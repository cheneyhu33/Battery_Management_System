# Cell thermal model with a central Gaussian heat source — portable C

One C source, one set of physics, buildable on any board with a C compiler.
Run it on each board, compare the numbers (which must be identical) and the
wall time (which is the measurement).

## Files

| File | What it is |
| --- | --- |
| `bms_gaussian.c` | the model — C99, links only against libm |
| `collect.sh` | runs the model on any board and packages the result |
| `bms_handoff.zip` | what you send to another board: `HANDOFF.md` + the two files above |
| `compare_boards.py` | compares result folders — physics and timing |
| `plot_heatmaps.py` | renders one board's CSVs into 2D figures |
| `plot_heatmap_3d.py` | renders the `--volume` field into 3D figures |
| `fig4_3_uniform.py` | renders Fig. 4.3 — the `--uniform` case, in the style of Fig. 4.1 |
| `fig4_2_field.py` | renders Fig. 4.2 — mid-plane, cut-away and the depth-slice strip |
| `fig_snap10_preview.py` | the mid-plane every 40 s from the dense-snapshot export |
| `cases/uniform-source/` | the `--uniform` run behind Fig. 4.3: its CSVs and log |
| `cases/gaussian-snap10/` | the Gaussian run exported every 10 s (40 frames, not 8) |
| `results/README.md` | the one-folder-per-board convention |
| `results/uno-q/` | **the UNO Q's results** — runs, `data/`, `plots/`, `environment.txt` |
| `main.py` | the UNO Q web relay — no physics, just the frame stream |

Every board's output lives in its own folder under `results/`, named after the
board. Ours is `results/uno-q/`; a classmate's arrives as
`results/<their-board>/`. Nothing is ever merged between folders.

## The model

```
rho*Cp dT/dt = k Laplacian(T) + q'''(r)
q'''(r) = A exp(-r^2 / (2 sigma^2)),   sum(q''' * W * dx^3) = Q_total
-k dT/dn = h (T - T_amb)   on all six faces
run stops when the hottest surface node reaches T_LIMIT
```

| Parameter | Value |
| --- | --- |
| Cell | 20 × 20 × 20 mm cube |
| Grid | 41 × 41 × 41 = 68 921 nodes, dx = 0.5 mm |
| k / rho / Cp | 1.04 W/(m K) / 1555 kg/m³ / 1437 J/(kg K) |
| Q_total | 2.0 W |
| sigma | 2.0 mm (peak q''' = 1.5873e7 W/m³) |
| h, T_amb | 10 W/(m² K), 25 °C |
| dt | 0.08018666 s (0.90 × stability limit) |
| Stop | surface ≥ 60 °C |

Vertex-centred finite volume: a node owns one cell in the interior, half a cell
on a face, a quarter on an edge and an eighth at a corner, and convects through
every face it exposes. On interior nodes this reduces exactly to the usual
seven-point FTCS stencil.

## Build and run

```bash
cc -std=c99 -O2 -o bms_gaussian bms_gaussian.c -lm

./bms_gaussian                 # timed run — this is the cross-board measurement
./bms_gaussian --export data   # same solve, plus the CSVs the 2D figures need
./bms_gaussian --volume vol    # same solve, plus the full 41^3 field for 3D
./bms_gaussian --stream 10     # same solve, plus binary frames on stdout

# Store a field frame every 10 s of simulated time instead of at the seven
# steps in SNAP. dt is untouched — this changes how often the field is
# stored, not how finely it is integrated.
./bms_gaussian --snap-dt 10 --export cases/gaussian-snap10/data \
                              --volume cases/gaussian-snap10/volume

# Same 2 W spread evenly instead of concentrated at the centre (Section 4.4).
# Physics, grid, time step and stop criterion are untouched: only the shape
# of q'''(r) changes, so it stays Q_total = 2 W either way.
./bms_gaussian --uniform
./bms_gaussian --uniform --export cases/uniform-source/data
```

For another architecture, add the cross-compiler and `-static`:

```bash
aarch64-linux-gnu-gcc -std=c99 -O2 -static -o bms_gaussian bms_gaussian.c -lm
```

## Why the modes are separate

A timing run must be a plain run. `--export` writes a row to a CSV on every
step; `--volume` writes 276 KB per snapshot; `--stream` sums the whole field and
writes 20 KB every ten steps. All three are real extra work. The program prints
which workload it measured so they can never be mistaken for each other:

```
Timing
  wall    : 21.960526 s
  per step: 4.5401 ms
  workload: plain solve, no file I/O inside the loop
```

Compare boards with the plain run only. Use `--export` and `--volume`
separately, when the board's output is needed for the figures.

`--snap-dt` does not change the solve, only how much of it is written: it sets
the frame cadence, so an export run with `--snap-dt 10` writes five times as
many frames as the stock `SNAP` table and takes correspondingly longer. Leave
it off for anything that is timed.

## Comparing boards

Send the other board `bms_handoff.zip` (or its three files). They extract it and run

```sh
sh collect.sh
```

which compiles, does three plain runs, one `--export` run, and writes
`bms_result_<hostname>/` with the logs, the CSVs and `environment.txt`
(CPU model, cpufreq governor and limits, memory, compiler version, load
average). They send that folder back. Nothing is installed on their board and
nothing outside the current directory is touched.

Rename it into its own slot before doing anything else, so the two boards can
never be confused later:

```bash
mv ~/Downloads/bms_result_beaglebone \
   results/beaglebone-ai64
```

Then, here:

```bash
python3 compare_boards.py results/uno-q results/beaglebone-ai64
```

The script enforces the rule that makes the comparison valid: it reads every
plain run of every board and exits non-zero if any physical quantity differs.
Only `wall` and `per step` are allowed to change. It prints

```
trip step      4837
trip time      387.8630 s
surface        60.0008 C
core           101.7759 C
core-surface   +41.7751 K
stored         632.5409 J
identical on all 2 boards

board                                          mean s    ms/step       vs first
UNO-Q                                         21.7402     4.4946          1.00x
Mac                                            0.3610     0.0774   60.2x faster
```

and, for boards that exported data, whether the slice CSVs came out
byte-identical.

If you would rather look by hand, the rule is the same: everything above the
`Timing` line must match digit for digit, and only `wall` / `per step` may
differ. If anything else differs, the boards are not running the same model,
and the timing comparison is meaningless.

For heatmaps, run `./bms_gaussian --export <dir>` on each board, copy the two
CSVs into `data/`, and run `plot_heatmaps.py`. The figures are then directly
comparable, having come from identical code and an identical renderer. A folder
that came back from another board is rendered in place, into its own `plots/`:

```bash
python3 plot_heatmaps.py results/beaglebone-ai64/data results/beaglebone-ai64/plots
```

Note that different architectures can round the last printed decimal
differently. Between the UNO Q and an Apple-silicon Mac the two CSVs differ by
at most `1e-10` °C — one unit in the last printed place, from summation order
in the energy column, not from a different solution. The slices are
byte-identical.

## Results on the Arduino UNO Q

| Quantity | Value |
| --- | --- |
| Trip time | 387.863 s of simulated time (6.46 min) |
| Trip step | 4837 |
| Surface at trip | 60.0008 °C |
| Core at trip | 101.7759 °C |
| Core − surface | +41.7751 K |
| Energy generated | 775.7258 J |
| Energy stored | 632.5409 J |
| Energy convected | 143.1849 J |
| Wall time, plain run | 21.740 s (4.495 ms/step), mean of three runs |
| Wall time, `--export` | 27.061 s (5.594 ms/step) |

The three plain runs were 21.634 / 21.957 / 21.629 s — a spread of 1.5 %, so
quote the mean and keep the spread in mind before reading anything into small
differences between boards.

For reference the same source compiled natively on an Apple-silicon Mac takes
0.361 s (0.077 ms/step, mean of three runs), 60.8× faster. That ratio is the
whole point of running it on more than one board.

The physics numbers do not depend on the board. The 41.8 K gap is a property
of the model, not of the hardware: the surface sensor reaches 60 °C while the
core is already at 101.8 °C.

## Figures

1. `01_temperature_evolution` — eight mid-plane snapshots on one shared scale.
2. `02_final_field_and_profile` — the frame the run froze on, and its
   centreline against the 60 °C gate.
3. `03_core_vs_surface` — the two temperatures through the whole run. The core
   jumps in the first 20 s (the Gaussian has almost no thermal mass) and then
   climbs steadily; the surface follows far behind.
4. `04_centreline_evolution` — the same centreline at every snapshot.

### 3D figures

`plot_heatmap_3d.py` reads a `--volume` dump instead of the slice CSV, so it
has the whole cube to work with:

```bash
./bms_gaussian --volume results/uno-q/volume      # once, to produce the field
python3 plot_heatmap_3d.py results/uno-q/volume results/uno-q/plots_3d
```

It needs numpy, matplotlib and scikit-image (marching cubes for the
isosurface). If the machine has none of them:

```bash
uv run --quiet --with numpy --with matplotlib --with scikit-image \
    python3 plot_heatmap_3d.py results/uno-q/volume results/uno-q/plots_3d
```

5. `05_3d_isosurface_final` — the 80 °C isosurface: everything inside it is
   hotter than the set point the surface sensor is protecting. In a 20 mm cube
   it is a sphere about 7 mm across, sitting well inside the skin that reads
   60 °C.
6. `06_3d_evolution` — the three mid-planes at six times on one shared colour
   scale, so the growth of the hot core is a colour change, not a rescaling.
7. `07_3d_midplanes` — the cut-away view: the three planes through the centre,
   where the core meets.

The 3D renderer draws surfaces rather than volumes: matplotlib's depth sorting
is unreliable for nested translucent shells, so the figures use one opaque
isosurface at a time and orthogonal planes, both of which it draws correctly.

## The UNO Q web page

The app framework on the UNO Q only runs Python on the MPU, so `main.py` is a
relay: it launches the compiled `bms_gaussian` as a child process, reads its
binary frame stream from stdout, and forwards each frame to the browser. It
contains no physics — no grid, no coefficients, no state — so the page shows
exactly what the C model computed.

To deploy it, the app directory looks like this:

```
ArduinoApps/bms/
  assets/                 unchanged (the page, the colour map, the timeline)
  python/main.py          the relay
  python/bms_gaussian.c   the model, for the record
  python/bin/bms_gaussian the model, cross-compiled for this board
```

The binary is built on the host and pushed, because the board has no Linux C
compiler:

```bash
zig cc -target aarch64-linux-musl -O2 -static -s -o bms_gaussian bms_gaussian.c -lm
adb push bms_gaussian /home/arduino/ArduinoApps/bms/python/bin/
adb shell chmod 755 /home/arduino/ArduinoApps/bms/python/bin/bms_gaussian
```

`main.py` starts the model on its first `loop()` call rather than at import:
the framework takes about twenty seconds to bring the web server up, and the
model takes about as long, so starting it earlier would mean it finished
before the page could be opened.

The web page runs the model at board speed — the whole 387.9 s of simulated
time in about 22 s of wall time — and then holds the final frame. The frozen
state is the same one the plain run reports.

There is therefore one implementation of the model, `bms_gaussian.c`, and
three ways to read it: plain run for timing, `--export` for heatmaps,
`--stream` for the live page.
