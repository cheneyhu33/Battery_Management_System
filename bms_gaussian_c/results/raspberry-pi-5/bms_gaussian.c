/* bms_gaussian.c — 3D thermal model of one 20 mm cell with a central
 * Gaussian heat source.
 *
 * This is a standalone C port of the model that runs on the Arduino UNO Q
 * (ArduinoApps/bms/python/main.py). The physics is identical — same grid,
 * same time step, same source normalisation, same boundary condition, same
 * stop criterion — so that the numbers from different boards can be compared
 * directly. It depends on nothing but libm.
 *
 *     rho*Cp dT/dt = k Laplacian(T) + q'''(r)
 *
 *     q'''(r) = A exp(-r^2 / (2 sigma^2))
 *     A chosen numerically so that  sum(q''' * W * dx^3) = Q_total
 *
 *     -k dT/dn = h (T - T_amb)     on all six faces
 *
 *     the run stops when the hottest node on the surface reaches T_LIMIT
 *
 * Discretisation: explicit FTCS on a vertex-centred finite volume grid. Each
 * node owns a control volume of 1 cell in the interior, 1/2 on a face, 1/4 on
 * an edge and 1/8 at a corner (W below), and sheds heat through every face it
 * exposes (CONV_COEF below). On interior nodes this reduces exactly to the
 * usual seven-point stencil.
 *
 * The same source builds on any board with a C compiler, and the numbers it
 * prints are identical everywhere; only the timing differs. That is the point
 * of the exercise.
 *
 * Build:  cc -O2 -o bms_gaussian bms_gaussian.c -lm
 *
 * Usage:
 *   ./bms_gaussian                 plain run: solve, report, print timing
 *   ./bms_gaussian --export DIR    also write the temperature history and
 *                                  mid-plane slices to CSV, for heatmaps
 *
 * Keep the two apart when comparing boards: a timing run must be a plain run,
 * because writing CSVs every step changes the workload.
 */

#define _POSIX_C_SOURCE 200809L   /* mkdir, clock_gettime */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>

/* ---------------------------------------------------------------- domain */

enum { N = 40, P = N + 1, PP = P * P };

static const double Lx = 20.0e-3;   /* cell edge [m] */
static const double Ly = 20.0e-3;
static const double Lz = 20.0e-3;

static const double dx = 20.0e-3 / 40;   /* = Lx / N; uniform grid [m] */
static const double dy = 20.0e-3 / 40;
static const double dz = 20.0e-3 / 40;

/* -------------------------------------------------------------- material */

static const double k   = 1.04;     /* thermal conductivity [W/(m K)] */
static const double rho = 1555.0;   /* density [kg/m^3] */
static const double Cp  = 1437.0;   /* specific heat [J/(kg K)] */

/* ----------------------------------------------------------- heat source */

static const double Q_total = 2.0;       /* total generation [W] */
static const double sigma   = 2.0e-3;    /* Gaussian radius [m] */

/* ------------------------------------------------------- boundary / stop */

static const double h     = 10.0;    /* convection coefficient [W/(m^2 K)] */
static const double T_amb = 25.0;    /* ambient air [C] */
static const double T_LIMIT = 60.0;  /* shutdown threshold on the surface [C] */

/* ---------------------------------------------------------------- state */

static size_t IDX(int i, int j, int k) { return ((size_t)i * P + j) * P + k; }

/* Steps at which a mid-plane snapshot is written in --export mode. The final
 * step, whenever it falls, is always written too. */
static const long SNAP[] = {0, 100, 500, 1000, 2000, 3000, 4000};
#define NSNAP ((int)(sizeof SNAP / sizeof SNAP[0]))

static int is_snapshot(long step)
{
    for (int n = 0; n < NSNAP; n++)
        if (SNAP[n] == step)
            return 1;
    return 0;
}

/* Hottest node on the six faces: what a surface-mounted sensor reads. */
static double surface_max_rise(const double *T)
{
    double s = 0.0;
    for (int j = 0; j < P; j++)
        for (int m = 0; m < P; m++) {
            double a = T[IDX(0, j, m)], b = T[IDX(N, j, m)];
            if (a > s) s = a;
            if (b > s) s = b;
        }
    for (int i = 0; i < P; i++)
        for (int m = 0; m < P; m++) {
            double a = T[IDX(i, 0, m)], b = T[IDX(i, N, m)];
            if (a > s) s = a;
            if (b > s) s = b;
        }
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++) {
            double a = T[IDX(i, j, 0)], b = T[IDX(i, j, N)];
            if (a > s) s = a;
            if (b > s) s = b;
        }
    return s;
}

/* ------------------------------------------------------- frame stream
 *
 * The UNO Q web page gets its data from this stream; the Python side is a
 * relay and holds no physics at all. Little-endian, no padding:
 *
 *   char   magic[4]   "BMSF"
 *   int32  nodes      P
 *   int32  step
 *   int32  final      1 on the last frame of the run
 *   double t, t_max, t_center, t_surf, energy, t_amb, t_limit
 *   float  xy[P*P]    T[i][j][MI]   temperature rise, index i*P + j
 *   float  xz[P*P]    T[i][MI][k]
 *   float  yz[P*P]    T[MI][j][k]
 *
 * All temperatures here are rises above ambient; the receiver adds t_amb
 * back for display. stdout carries the binary stream, so in this mode every
 * human-readable line goes to stderr instead.
 */
static void emit_frame(FILE *f, long step, int final, double dt,
                       const double *T, const double *W, double cell_volume,
                       double t_surf)
{
    const size_t nn = (size_t)P * P * P;
    double tmax = T[0], energy = 0.0;
    for (size_t c = 0; c < nn; c++) {
        if (T[c] > tmax) tmax = T[c];
        energy += W[c] * T[c];
    }
    energy *= rho * Cp * cell_volume;

    const int nodes = P, st = (int)step;
    const double hdr[7] = {step * dt, tmax, T[IDX(N / 2, N / 2, N / 2)],
                           t_surf, energy, T_amb, T_LIMIT};
    float buf[P * P];

    fwrite("BMSF", 1, 4, f);
    fwrite(&nodes, sizeof nodes, 1, f);
    fwrite(&st, sizeof st, 1, f);
    fwrite(&final, sizeof final, 1, f);
    fwrite(hdr, sizeof hdr, 1, f);

    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++)
            buf[i * P + j] = (float)T[IDX(i, j, N / 2)];
    fwrite(buf, sizeof buf, 1, f);

    for (int i = 0; i < P; i++)
        for (int m = 0; m < P; m++)
            buf[i * P + m] = (float)T[IDX(i, N / 2, m)];
    fwrite(buf, sizeof buf, 1, f);

    for (int j = 0; j < P; j++)
        for (int m = 0; m < P; m++)
            buf[j * P + m] = (float)T[IDX(N / 2, j, m)];
    fwrite(buf, sizeof buf, 1, f);

    fflush(f);
}

/* One row of the time history. Energy is summed on demand because the extra
 * pass costs time and would otherwise perturb the timing comparison. */
static void write_history(FILE *f, long step, double dt, const double *T,
                          const double *W, const double cell_volume,
                          const double *av, const double *ex)
{
    double tmin = T[0], tmax = T[0], sur = 0.0, energy = 0.0;
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++)
            for (int m = 0; m < P; m++) {
                const size_t c = IDX(i, j, m);
                const double v = T[c];
                if (v < tmin) tmin = v;
                if (v > tmax) tmax = v;
                energy += W[c] * v;
                if (ex[i] > 0.0 || ex[j] > 0.0 || ex[m] > 0.0)
                    if (v > sur) sur = v;
            }
    (void)av;
    fprintf(f, "%ld,%.6f,%.10f,%.10f,%.10f,%.10f,%.10f\n", step, step * dt,
            T_amb + T[IDX(N / 2, N / 2, N / 2)], T_amb + sur, T_amb + tmin,
            T_amb + tmax, energy * rho * Cp * cell_volume);
}

/* The z = centre plane, one row per node, absolute temperature. */
static void write_slice(FILE *f, long step, const double *T)
{
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++)
            fprintf(f, "%ld,%d,%d,%.10f\n", step, i, j,
                    T_amb + T[IDX(i, j, N / 2)]);
}

int main(int argc, char **argv)
{
    FILE *out = stdout;
    const char *export_dir = NULL;
    int stream_mode = 0;
    long stream_every = 10;

    for (int a = 1; a < argc; a++) {
        if (strcmp(argv[a], "--export") == 0 && a + 1 < argc) {
            export_dir = argv[++a];
        } else if (strcmp(argv[a], "--stream") == 0) {
            stream_mode = 1;
            if (a + 1 < argc && argv[a + 1][0] != '-')
                stream_every = strtol(argv[++a], NULL, 10);
        } else if (strcmp(argv[a], "-h") == 0 || strcmp(argv[a], "--help") == 0) {
            fprintf(out, "usage: %s [--export DIR] [--stream [STEPS]]\n", argv[0]);
            return 0;
        } else {
            fprintf(stderr,
                    "unknown argument: %s\nusage: %s [--export DIR] [--stream [STEPS]]\n",
                    argv[a], argv[0]);
            return 2;
        }
    }
    if (stream_mode && stream_every < 1)
        stream_every = 1;

    /* The frame stream owns stdout, so keep the log on stderr. */
    if (stream_mode)
        out = stderr;

    FILE *hist = NULL, *slice = NULL;
    if (export_dir) {
        if (mkdir(export_dir, 0777) != 0 && errno != EEXIST) {
            fprintf(stderr, "cannot create %s: %s\n", export_dir, strerror(errno));
            return 1;
        }
        char path[1024];
        snprintf(path, sizeof path, "%s/temperature_history.csv", export_dir);
        hist = fopen(path, "w");
        if (!hist) {
            fprintf(stderr, "cannot write %s: %s\n", path, strerror(errno));
            return 1;
        }
        snprintf(path, sizeof path, "%s/midplane_slices.csv", export_dir);
        slice = fopen(path, "w");
        if (!slice) {
            fprintf(stderr, "cannot write %s: %s\n", path, strerror(errno));
            return 1;
        }
        fprintf(hist, "step,time_s,center_C,surface_max_C,min_C,max_C,energy_J\n");
        fprintf(slice, "step,i,j,T_C\n");
    }

    const double alpha = k / (rho * Cp);

    /* Explicit stability bound. The worst node is a corner: three conducting
     * neighbours but three convecting faces carrying only an eighth of a
     * cell of heat capacity. The source does not enter this bound. */
    const double stability_limit =
        1.0 / (2.0 * alpha * (1.0 / (dx * dx) + 1.0 / (dy * dy) + 1.0 / (dz * dz))
               + 2.0 * h / (rho * Cp) * (1.0 / dx + 1.0 / dy + 1.0 / dz));
    const double dt = 0.90 * stability_limit;

    /* Per-axis volume fraction (1 interior, 1/2 at the two ends) and the
     * matching exposed-face count (0 interior, 1 at the ends). */
    double av[P], ex[P], sa[P];
    for (int i = 0; i < P; i++) {
        av[i] = (i == 0 || i == N) ? 0.5 : 1.0;
        ex[i] = (i == 0 || i == N) ? 1.0 : 0.0;
        sa[i] = ex[i] / av[i];          /* 2 on a face node, 0 inside */
    }

    const double rx = alpha * dt / (dx * dx);
    const double ry = alpha * dt / (dy * dy);
    const double rz = alpha * dt / (dz * dz);

    /* Convective loss coefficient: -(h/(rho*Cp*dx)) * CONV_COEF * theta */
    const double conv_base = h * dt / (rho * Cp * dx);

    const size_t nn = (size_t)P * P * P;
    const double cell_volume = dx * dy * dz;

    double *T = malloc(nn * sizeof *T);
    double *Tn = malloc(nn * sizeof *Tn);
    double *heat = malloc(nn * sizeof *heat);
    double *W = malloc(nn * sizeof *W);
    if (!T || !Tn || !heat || !W) {
        fprintf(stderr, "out of memory\n");
        return 1;
    }

    /* Gaussian shape, then normalised on the discrete sum. The analytic
     * integral would be wrong here: the Gaussian is truncated at the faces
     * and the boundary control volumes are fractions of a cell. Dividing by
     * the discrete sum makes sum(q''' * W * dx^3) exactly Q_total. */
    double shape_sum = 0.0;
    for (int i = 0; i < P; i++) {
        const double xc = i * dx - 0.5 * Lx;
        for (int j = 0; j < P; j++) {
            const double yc = j * dy - 0.5 * Ly;
            for (int m = 0; m < P; m++) {
                const double zc = m * dz - 0.5 * Lz;
                const size_t c = IDX(i, j, m);
                const double sh =
                    exp(-(xc * xc + yc * yc + zc * zc) / (2.0 * sigma * sigma));
                W[c] = av[i] * av[j] * av[m];
                T[c] = sh;               /* park the shape here for now */
                shape_sum += sh * W[c];
            }
        }
    }
    const double source_norm = shape_sum * cell_volume;   /* effective volume */

    double q_peak = 0.0;
    for (size_t c = 0; c < nn; c++) {
        const double qv = Q_total / source_norm * T[c];
        if (qv > q_peak)
            q_peak = qv;
        heat[c] = qv / (rho * Cp) * dt;   /* temperature rise per step [K] */
        T[c] = 0.0;                        /* start at ambient */
    }

    /* ------------------------------------------------------------ banner */
    fprintf(out, "=======================================================================\n");
    fprintf(out, "bms - single-cell 3D heat diffusion, central Gaussian heat source\n");
    fprintf(out, "=======================================================================\n");
    fprintf(out, "Thermal properties\n");
    fprintf(out, "  k     = %.4f W/(m K)\n", k);
    fprintf(out, "  rho   = %.1f kg/m^3\n", rho);
    fprintf(out, "  Cp    = %.1f J/(kg K)\n", Cp);
    fprintf(out, "  alpha = %.4e m^2/s\n", alpha);
    fprintf(out, "dx    = %.4e m\n", dx);
    fprintf(out, "dt    = %.6e s  (0.90 x stability limit %.6e s)\n", dt, stability_limit);
    fprintf(out, "r     = %.6f   (pure-conduction limit 1/6 = %.6f)\n", rx, 1.0 / 6.0);
    fprintf(out, "-----------------------------------------------------------------------\n");
    fprintf(out, "Cell      : %.1f x %.1f x %.1f mm cube\n", Lx * 1e3, Ly * 1e3, Lz * 1e3);
    fprintf(out, "Grid      : %d^3 = %zu nodes  (dx = %.2f mm)\n", P, nn, dx * 1e3);
    fprintf(out, "Ambient   : %.1f C\n", T_amb);
    fprintf(out, "Cooling   : natural convection, h = %.1f W/(m^2 K)\n", h);
    fprintf(out, "Biot      : Bi = h*L/k = %.3f\n", h * Lx / k);
    fprintf(out, "Heat source\n");
    fprintf(out, "  total   : %.2f W  (%.0f W/kg)\n", Q_total, Q_total / (rho * Lx * Ly * Lz));
    fprintf(out, "  shape   : Gaussian about the centre, sigma = %.2f mm\n", sigma * 1e3);
    fprintf(out, "  peak    : %.4e W/m^3  (%.3f K/s at the core, not everywhere)\n",
           q_peak, q_peak / (rho * Cp));
    fprintf(out, "  spread  : effective source volume = %.2f%% of the cell\n",
           source_norm / (Lx * Ly * Lz) * 100.0);
    fprintf(out, "Shutdown\n");
    fprintf(out, "  trigger : surface >= %.1f C\n", T_LIMIT);
    fprintf(out, "  action  : stop advancing the field\n");
    fprintf(out, "Time constant rho*Cp*V/(h*A) = %.2f s\n",
           rho * Cp * Lx * Ly * Lz / (h * 6.0 * Lx * Lx));
    fprintf(out, "=======================================================================\n");
    fprintf(out, "Initial condition: cell at %.1f C, central heat source ON\n", T_amb);
    fprintf(out, "  source  : %.2f W, Gaussian at the centre, sigma = %.1f mm\n",
           Q_total, sigma * 1e3);
    fprintf(out, "  peak    : %.4e W/m^3 at the core\n", q_peak);
    fprintf(out, "  stop    : surface >= %.1f C\n", T_LIMIT);
    fprintf(out, "-----------------------------------------------------------------------\n");

    /* -------------------------------------------------------------- solve */
    long step = 0;
    double t_surf = 0.0, t_core = 0.0;

    if (hist) {
        write_history(hist, 0, dt, T, W, cell_volume, av, ex);
        write_slice(slice, 0, T);
    }

    if (stream_mode)
        emit_frame(stdout, 0, 0, dt, T, W, cell_volume, 0.0);

    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC, &t_start);

    for (;;) {
        for (int i = 0; i < P; i++) {
            const double rc_x = rx / av[i];
            for (int j = 0; j < P; j++) {
                const double rc_y = ry / av[j];
                const double cxy = conv_base * (sa[i] + sa[j]);
                for (int m = 0; m < P; m++) {
                    const size_t c = IDX(i, j, m);
                    const double t = T[c];
                    double v = t;

                    if (i > 0) v += rc_x * (T[c - PP] - t);
                    if (i < N) v += rc_x * (T[c + PP] - t);
                    if (j > 0) v += rc_y * (T[c - P] - t);
                    if (j < N) v += rc_y * (T[c + P] - t);
                    if (m > 0) v += (rz / av[m]) * (T[c - 1] - t);
                    if (m < N) v += (rz / av[m]) * (T[c + 1] - t);

                    v -= (cxy + conv_base * sa[m]) * t;
                    v += heat[c];

                    Tn[c] = v;
                }
            }
        }

        double *swap = T;
        T = Tn;
        Tn = swap;
        step++;

        t_surf = surface_max_rise(T);

        if (hist)
            write_history(hist, step, dt, T, W, cell_volume, av, ex);
        if (slice && is_snapshot(step))
            write_slice(slice, step, T);
        if (stream_mode && step % stream_every == 0)
            emit_frame(stdout, step, 0, dt, T, W, cell_volume, t_surf);

        if (step % 500 == 0) {
            t_core = T_amb + T[IDX(N / 2, N / 2, N / 2)];
            fprintf(out, "[step %6ld] t = %10.3f s | T_core = %9.4f C | T_surf = %9.4f C\n",
                   step, step * dt, t_core, T_amb + t_surf);
        }

        if (T_amb + t_surf >= T_LIMIT)
            break;
    }

    clock_gettime(CLOCK_MONOTONIC, &t_end);
    const double wall =
        (double)(t_end.tv_sec - t_start.tv_sec)
        + 1e-9 * (double)(t_end.tv_nsec - t_start.tv_nsec);

    /* The final field is always part of the export, even when the trip step
     * is not one of the scheduled snapshots. */
    if (slice && !is_snapshot(step))
        write_slice(slice, step, T);

    /* The last frame is flagged so the receiver can tell a finished run from
     * one that is still going. */
    if (stream_mode)
        emit_frame(stdout, step, 1, dt, T, W, cell_volume, t_surf);

    /* ------------------------------------------------------------ report */
    t_surf = T_amb + t_surf;
    t_core = T_amb + T[IDX(N / 2, N / 2, N / 2)];

    double stored = 0.0, peak = 0.0;
    for (int i = 0; i < P; i++)
        for (int j = 0; j < P; j++)
            for (int m = 0; m < P; m++) {
                const size_t c = IDX(i, j, m);
                stored += W[c] * T[c];
                if (T[c] > peak)
                    peak = T[c];
            }
    stored *= rho * Cp * cell_volume;

    fprintf(out, "-----------------------------------------------------------------------\n");
    fprintf(out, "  == LIMIT REACHED  t = %.3f s   (step %ld)\n", step * dt, step);
    fprintf(out, "     surface      : %.4f C   (limit %.1f C)\n", t_surf, T_LIMIT);
    fprintf(out, "     centre       : %.4f C   (%+.4f K above the surface)\n",
           t_core, t_core - t_surf);
    fprintf(out, "     -> simulation frozen at this frame\n");
    fprintf(out, "-----------------------------------------------------------------------\n");
    fprintf(out, "Peak cell temperature   : %.4f C\n", T_amb + peak);
    fprintf(out, "Energy stored above amb : %.4f J\n", stored);
    fprintf(out, "Energy generated        : %.4f J\n", Q_total * step * dt);
    fprintf(out, "Energy convected away   : %.4f J\n", Q_total * step * dt - stored);
    fprintf(out, "-----------------------------------------------------------------------\n");
    fprintf(out, "Timing\n");
    fprintf(out, "  wall    : %.6f s\n", wall);
    fprintf(out, "  per step: %.4f ms\n", wall / (double)step * 1e3);
    fprintf(out, "  workload: %s\n",
           (export_dir && stream_mode)
               ? "solve + CSV export + frame stream (not comparable with a plain run)"
           : export_dir
               ? "solve + CSV export (not comparable with a plain run)"
           : stream_mode
               ? "solve + frame stream (not comparable with a plain run)"
               : "plain solve, no file I/O inside the loop");

    if (hist)
        fclose(hist);
    if (slice)
        fclose(slice);

    free(T);
    free(Tn);
    free(heat);
    free(W);
    return 0;
}
