#!/usr/bin/env python3
"""Figures 6.3, 6.4 and 6.6 of the EE3180 report, rendered from the model's own output.

Nothing here is hand-edited: every curve comes from a computed field.

    fig6_3  Gaussian-pulse convergence test. The model profile is produced by
            this script, running the same vertex-centred finite-volume update
            as bms_gaussian.c (h = 0, no source); the exact solution is
            Eq. (6.3). The convergence panel uses the same three grids as
            Table 6.2, and reproduces its errors (5.24e-2 / 1.40e-2 / 3.52e-3).

    fig6_4  Core and surface temperature, read from
            results/uno-q/data/temperature_history.csv.

    fig6_6  The frame at the trip, read from results/uno-q/data/midplane_slices.csv.

    python3 report_figures.py [RESULTS_DIR] [OUT_DIR]

Needs numpy and matplotlib only.
"""
from pathlib import Path
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.ticker import FixedLocator, NullFormatter

# ---------------------------------------------------------------- style
BG = "#f7f8fc"          # figure face, the light grey-blue of the report figures
INK = "#172738"         # titles and axis text
MUTED = "#64748b"       # notes
GRID = "#e3e7ef"
CORE = "#b54812"        # core / hot
SURF = "#147a89"        # surface / cool
GATE = "#7d8791"        # the 60 C gate

# ---------------------------------------------------------------- model constants
K_COND = 1.04
RHO = 1555.0
CP = 1437.0
ALPHA = K_COND / (RHO * CP)
L_M = 20.0e-3
P_GRID = 41


# ======================================================================
# 1. Gaussian-pulse verification (Fig. 6.3)
# ======================================================================
def pulse_run(n, a0=35.0, sigma0=1.0e-3, sigma_end=2.5e-3):
    """Same discrete operator as bms_gaussian.c with h = 0 and no source.

    A node's temperature moves by (alpha*dt/dx^2)/av * sum(neighbour - self),
    where av = 1/2 on the outer planes; that halves the heat capacity of a
    boundary node, which is the zero-flux (mirror) condition. Returns the grid,
    the field at the final time and that time.
    """
    p = n + 1
    dx = L_M / n
    dt = 0.90 / (2.0 * ALPHA * 3.0 / (dx * dx))     # 0.9 x the corner bound
    t_end = (sigma_end ** 2 - sigma0 ** 2) / (2.0 * ALPHA)
    # t_end/dt is exactly 70 at N = 40 and 280 at N = 80; a rounding hair above
    # the integer would make ceil take one extra step and land past the time the
    # test is defined at, so nudge it down before rounding up.
    steps = int(np.ceil(t_end / dt - 1e-9))
    t_final = steps * dt
    r = ALPHA * dt / (dx * dx)

    x = np.arange(p) * dx - 0.5 * L_M
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    T = a0 * np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / (2.0 * sigma0 ** 2))

    def axis_term(T, r, axis):
        d = np.zeros_like(T)
        sl = [slice(None)] * 3
        # interior: r*(T[i-1] + T[i+1] - 2T[i])
        inner = [slice(None)] * 3
        inner[axis] = slice(1, -1)
        lo = [slice(None)] * 3; lo[axis] = slice(0, -2)
        hi = [slice(None)] * 3; hi[axis] = slice(2, None)
        d[tuple(inner)] += r * (T[tuple(lo)] + T[tuple(hi)] - 2.0 * T[tuple(inner)])
        # ends: (r/av)*(neighbour - self) with av = 1/2 -> 2r
        first = [slice(None)] * 3; first[axis] = 0
        second = [slice(None)] * 3; second[axis] = 1
        d[tuple(first)] += 2.0 * r * (T[tuple(second)] - T[tuple(first)])
        last = [slice(None)] * 3; last[axis] = -1
        prev = [slice(None)] * 3; prev[axis] = -2
        d[tuple(last)] += 2.0 * r * (T[tuple(prev)] - T[tuple(last)])
        return d

    for _ in range(steps):
        T = T + axis_term(T, r, 0) + axis_term(T, r, 1) + axis_term(T, r, 2)

    return x, T, t_final


def pulse_reference(t, a0=35.0, sigma0=1.0e-3):
    w2 = sigma0 ** 2 + 2.0 * ALPHA * t
    return a0 * (sigma0 ** 2 / w2) ** 1.5, w2


def figure_6_3(out):
    grids = [(20, 1.00), (40, 0.50), (80, 0.25)]     # (N, dx in mm)
    errors, peaks, exact = [], [], []
    profile = None
    for n, dx_mm in grids:
        x, T, t = pulse_run(n)
        amp, w2 = pulse_reference(t)
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        ref = amp * np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / (2.0 * w2))
        errors.append(float(np.abs(T - ref).max()))
        peaks.append(float(T.max()))
        exact.append(amp)
        if n == 40:
            mid = n // 2
            profile = ((x + 0.5 * L_M) * 1e3,
                       T[:, mid, mid].copy(),
                       ref[:, mid, mid].copy(),
                       t)

    t_final = profile[3]
    with plt.rc_context({"font.family": "serif",
                         "axes.edgecolor": "#333333",
                         "axes.labelcolor": "#111111",
                         "text.color": "#111111",
                         "xtick.color": "#111111",
                         "ytick.color": "#111111"}):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

        # ---- (a) the profile, model against exact -----------------------
        xs, model, ref, _ = profile
        x0 = np.linspace(0.0, 20.0, 400)
        t0 = 35.0 * np.exp(-(((x0 - 10.0) * 1e-3) ** 2) / (2.0 * (1.0e-3) ** 2))
        ax1.plot(x0, t0, color="#8a8a8a", lw=1.6, label="t = 0")
        ax1.plot(xs, ref, color="#9a9a9a", lw=4.0, alpha=0.75,
                 label="analytical (%.2f s)" % t_final)
        ax1.plot(xs, model, "o", ms=3.2, color="#1f77b4",
                 label="model, N = 40")
        ax1.set_xlabel("x [mm] (through centre)")
        ax1.set_ylabel("temperature rise [K]")
        ax1.set_xlim(0, 20)
        ax1.set_ylim(-1, 36)
        ax1.grid(color="#e6e6e6", lw=0.8)
        ax1.legend(loc="center left", frameon=True, fontsize=11)
        ax1.set_title("(a) pulse profile", fontsize=15)

        # ---- (b) convergence -------------------------------------------
        dx = np.array([g[1] for g in grids])
        err = np.array(errors)
        ax2.loglog(dx, err, "o-", color="#1f77b4", ms=9, lw=1.8,
                   label="max |error|")
        ax2.loglog(dx, err[-1] * (dx / dx[-1]) ** 2, "k--", lw=1.3,
                   label="slope 2")
        ax2.set_xlabel("dx [mm]")
        ax2.set_ylabel("max |error| [K]")
        # Explicit ticks: the default log locator prints 3e-4 and 4e-1 on top
        # of each other at this width.
        ax2.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1.0]))
        ax2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
            lambda v, _: ("%g" % v)))
        ax2.xaxis.set_minor_formatter(NullFormatter())
        ax2.set_xlim(0.20, 1.25)
        ax2.grid(color="#e6e6e6", lw=0.8, which="major")
        ax2.legend(loc="upper left", frameon=True, fontsize=12)
        ax2.set_title("(b) grid convergence", fontsize=15)

        fig.tight_layout()
        fig.savefig(out / "fig6_3_grid_convergence.png", dpi=210,
                    facecolor="white")
        plt.close(fig)

    for (n, dx_mm), e in zip(grids, errors):
        print("  N=%2d  dx=%.2f mm  max|err| = %.4e K" % (n, dx_mm, e))
    return errors


# ======================================================================
# 2. Core and surface temperature (Fig. 6.4)
# ======================================================================
def read_history(path):
    t, core, surf = [], [], []
    with open(path) as f:
        for row in csv.DictReader(f):
            t.append(float(row["time_s"]))
            core.append(float(row["center_C"]))
            surf.append(float(row["surface_max_C"]))
    return np.array(t), np.array(core), np.array(surf)


def figure_6_4(out, data_dir):
    t, core, surf = read_history(data_dir / "temperature_history.csv")
    trip = len(t) - 1
    t_trip, c_trip, s_trip = t[trip], core[trip], surf[trip]
    gap = c_trip - s_trip

    fig = plt.figure(figsize=(10.4, 5.8), facecolor=BG)
    ax = fig.add_axes((0.075, 0.20, 0.905, 0.62), facecolor="white")

    ax.plot(t, core, color=CORE, lw=2.6, solid_capstyle="round",
            label="Core (centre node)")
    ax.plot(t, surf, color=SURF, lw=2.6, solid_capstyle="round",
            label="Hottest surface node (the sensor)")
    ax.axhline(60.0, color=GATE, lw=1.3, ls=(0, (4, 3)))
    ax.axvline(t_trip, color=GATE, lw=1.0, ls=(0, (1, 2)))

    ax.set_xlim(0, t_trip * 1.02)
    ax.set_ylim(24, 105)
    ax.set_xlabel("Simulated time (s)", color=INK, fontsize=13)
    ax.set_ylabel("Temperature (C)", color=INK, fontsize=13)
    ax.tick_params(colors=INK, labelsize=11)
    ax.grid(color=GRID, lw=1.0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#c9d0dc")

    leg = ax.legend(loc="upper left", frameon=False, fontsize=12.5)
    for text in leg.get_texts():
        text.set_color(INK)

    # The trip annotation sits in the empty wedge between the two curves. No
    # leader line: a straight segment drawn across the plot from here to the
    # trip point reads as a stray trend line.
    ax.text(0.50, 0.62,
            "trip  %.2f s\nsurface  %.2f C\ncore  %.2f C"
            % (t_trip, s_trip, c_trip),
            transform=ax.transAxes, fontsize=13, color=INK, linespacing=1.5,
            va="center", ha="left",
            bbox=dict(boxstyle="round,pad=0.45", facecolor="white",
                      edgecolor="#dfe4ee", lw=1.0))

    fig.text(0.012, 0.055,
             "The gap is the sensor lag: %.1f K at the trip, and it is still widening."
             % gap, fontsize=10.5, color=MUTED)
    fig.text(0.012, 0.018,
             "Central Gaussian source (sigma = 2 mm, 2 W) in a 20 mm cell; "
             "mid-plane z = 10 mm. Gate: the surface reaches 60 C.",
             fontsize=10.5, color=MUTED)
    fig.text(0.5, 0.925, "Surface sensor versus core, up to the trip",
             ha="center", va="center", fontsize=20, fontweight="bold", color=INK)

    fig.savefig(out / "fig6_4_core_vs_surface.png", dpi=210, facecolor=BG)
    plt.close(fig)
    print("  trip %.2f s   core %.2f C   surface %.2f C   gap %.2f K"
          % (t_trip, c_trip, s_trip, gap))


# ======================================================================
# 3. The frame at the trip (Fig. 6.6)
# ======================================================================
def read_slice(path, step):
    P = P_GRID
    field = np.full((P, P), np.nan)
    with open(path) as f:
        for row in csv.DictReader(f):
            if int(row["step"]) != step:
                continue
            field[int(row["i"]), int(row["j"])] = float(row["T_C"])
    return field


def figure_6_6(out, data_dir):
    with open(data_dir / "midplane_slices.csv") as f:
        steps = sorted({int(r["step"]) for r in csv.DictReader(f)})
    step = steps[-1]
    field = read_slice(data_dir / "midplane_slices.csv", step)
    P = P_GRID
    g = np.linspace(0.0, 20.0, P)
    mid = P // 2
    centreline = field[:, mid]
    x_mm = g

    vmin, vmax = float(np.nanmin(field)), float(np.nanmax(field))
    t_trip = 387.86

    fig = plt.figure(figsize=(12.8, 5.78), facecolor=BG)
    # The colourbar gets its own row under the left panel: placed afterwards
    # with add_axes it collides with the footnotes at this figure height.
    gs = fig.add_gridspec(2, 2, height_ratios=(1.0, 0.055),
                          width_ratios=(1.0, 1.30),
                          left=0.055, right=0.975, top=0.80, bottom=0.205,
                          wspace=0.22, hspace=0.62)

    # ---- left: the mid-plane field ------------------------------------
    axL = fig.add_subplot(gs[0, 0])
    mesh = axL.pcolormesh(g, g, field.T, cmap="inferno", shading="gouraud",
                          vmin=vmin, vmax=vmax)
    axL.contour(g, g, field.T, levels=[60, 70, 80, 90, 100],
                colors="white", linewidths=0.9, alpha=0.75)
    axL.set_aspect("equal")
    axL.set_xlabel("x (mm)", color=INK, fontsize=12)
    axL.set_ylabel("y (mm)", color=INK, fontsize=12)
    axL.set_title("Mid-plane, %.2f to %.2f C" % (vmin, vmax),
                  fontsize=13, color=INK)
    axL.tick_params(colors=INK, labelsize=10)

    cax = fig.add_subplot(gs[1, 0])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.set_label("Temperature (C)", color=INK, fontsize=11, labelpad=6)
    cb.ax.tick_params(colors=INK, labelsize=10)
    cb.outline.set_edgecolor("#c9d0dc")

    # ---- right: the centreline ----------------------------------------
    axR = fig.add_subplot(gs[0, 1])
    axR.plot(x_mm, centreline, color=CORE, lw=2.6, solid_capstyle="round")
    axR.plot([10.0], [vmax], "o", ms=9, color=INK, zorder=5)
    axR.axhline(60.0, color=SURF, lw=1.4, ls=(0, (4, 3)))
    axR.set_xlim(0, 20)
    axR.set_ylim(56, 106)
    axR.set_xlabel("x (mm)", color=INK, fontsize=12)
    axR.set_ylabel("Temperature (C)", color=INK, fontsize=12)
    axR.set_title("Centreline through y = z = 10 mm", fontsize=13, color=INK)
    axR.tick_params(colors=INK, labelsize=10)
    axR.grid(color=GRID, lw=1.0)
    axR.set_axisbelow(True)
    for side in ("top", "right"):
        axR.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axR.spines[side].set_color("#c9d0dc")

    # Annotations placed where the curve is not: right of the peak but above
    # it, and on the gate line well right of the point where the curve crosses.
    axR.annotate("core  %.2f C" % vmax, xy=(10.0, vmax),
                 xytext=(11.4, 103.4), fontsize=11.5, color=INK,
                 va="bottom", ha="left")
    axR.text(13.0, 60.7, "60 C gate", fontsize=11.5, color=SURF,
             va="bottom", ha="left")

    fig.text(0.5, 0.925, "At the moment the sensor trips",
             ha="center", va="center", fontsize=20, fontweight="bold", color=INK)
    fig.text(0.012, 0.045,
             "The surface reads %.2f C while the core is %.2f C: a gap of %.1f K."
             % (np.nanmin(centreline), vmax, vmax - np.nanmin(centreline)),
             fontsize=10.5, color=MUTED)
    fig.text(0.012, 0.012,
             "Central Gaussian source (sigma = 2 mm, 2 W) in a 20 mm cell; "
             "mid-plane z = 10 mm. Gate: the surface reaches 60 C.",
             fontsize=10.5, color=MUTED)

    fig.savefig(out / "fig6_6_trip_frame.png", dpi=200, facecolor=BG)
    plt.close(fig)
    print("  step %d   mid-plane %.2f to %.2f C   surface %.4f C"
          % (step, vmin, vmax, centreline[0]))


def main(argv):
    root = Path(__file__).resolve().parent
    data = Path(argv[1]).resolve() if len(argv) > 1 else root / "results/uno-q/data"
    out = Path(argv[2]).resolve() if len(argv) > 2 else root / "report_figures"
    out.mkdir(parents=True, exist_ok=True)

    print("fig 6.3 (pulse convergence)")
    figure_6_3(out)
    print("fig 6.4 (core vs surface)")
    figure_6_4(out, data)
    print("fig 6.6 (trip frame)")
    figure_6_6(out, data)
    print("written to", out)


if __name__ == "__main__":
    main(sys.argv)
