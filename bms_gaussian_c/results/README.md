# Results, one folder per board

Every board that runs the model gets **its own folder here, named after the
board**. Nothing is ever merged, moved between folders, or overwritten by
another board's run.

```
results/
  uno-q/                 Arduino UNO Q          <- ours
  raspberry-pi-5/        Raspberry Pi 5         <- arrives from the classmate
  beaglebone-ai64/       BeagleBone AI-64       <- arrives from the classmate
  <board>/               anything else, same layout
```

## What is inside a folder

Identical layout on every board, because `collect.sh` produces it:

| File | What it is |
| --- | --- |
| `run_1.txt` `run_2.txt` `run_3.txt` | three plain runs — the timing measurement |
| `run_export.txt` | the run that also wrote the CSVs (timing not comparable) |
| `data/temperature_history.csv` | per-step core / surface / min / max / energy |
| `data/midplane_slices.csv` | 8 snapshots × 3 mid-planes, full field |
| `environment.txt` | CPU, clock governor, memory, compiler, binary and source hashes |
| `bms_gaussian.c` `collect.sh` | the source and script that actually ran |
| `plots/` | the 2D figures — in `uno-q/` only, because they come out byte-identical on every board (see below) |
| `volume/` | the full P³ field at the snapshot steps, for the 3D figures |
| `plots_3d/` | the 3D figures drawn from `volume/` |

`volume/` exists only where someone ran `--volume`. On the UNO Q it was produced
by a host build, because the binary deployed on the board was built before that
mode existed and the board has no C compiler — see
`uno-q/volume/PROVENANCE.md`. The field itself is board-independent: the plain
run prints the same numbers everywhere and the slice CSVs are byte-identical.

`environment.txt` is what identifies a folder: it carries the board model,
architecture, CPU name and the hashes of the binary and the source. If a
returned folder has no `environment.txt`, its provenance is unknown — ask for
it again rather than trusting the numbers.

## One set of figures is enough

The boards compute the same temperatures and the figures are drawn by the same
script from those numbers, so the figures are not merely similar — they are the
same file. All three boards produce these hashes:

| Figure | md5 |
| --- | --- |
| `01_temperature_evolution.png` | `bc87dc6242937679b943e71c837dc0fc` |
| `02_final_field_and_profile.png` | `06e8112839ecf70d07ddf351299e2ca9` |
| `03_core_vs_surface.png` | `76afb12853f15c1fec6c0d585dea9d76` |
| `04_centreline_evolution.png` | `964954f7304f4f55d56c774f29e36f65` |

Storing three copies of the same four files adds nothing, so only
`uno-q/plots/` is kept. The `data/` CSVs stay in **every** folder: they are the
record that that board actually ran, and any board's figures can be rebuilt
from its own data with one command:

```bash
python3 plot_heatmaps.py results/raspberry-pi-5/data results/raspberry-pi-5/plots
```

The history CSVs are not byte-identical — the UNO Q's differs from the other
two in the last printed digit, up to `1e-10` (a summation-order effect, noted
in `bms_gaussian_c/README.md`) — but that is far below anything a raster at
210 dpi can show, so the images still come out identical. If a re-rendered set
ever fails to match the hashes above, the data behind it is different, and that
is worth investigating rather than ignoring.

## How our folder differs from a returned one

Ours was produced by running a **cross-compiled** binary rather than compiling
on the board, because the UNO Q has no Linux C compiler. `environment.txt`
says so explicitly:

```
method    : the board has no Linux C compiler, so this binary was
            cross-compiled on macOS with zig 0.17.0 and run directly.
            command: zig cc -target aarch64-linux-musl -O2 -static -s -lm
```

That is the only difference. The model, the grid, the arithmetic and the
output are the same, and the exported slices were byte-identical to the ones
produced by a native build on an Apple-silicon Mac.

## Two revisions of the source are in here

The `bms_gaussian.c` that was sent out is **not** the one at the top of this
folder any more.

| Revision | md5 | Where it is |
| --- | --- | --- |
| sent out / what ran | `b4ca229122dd130fcd2d8540fbb637c5` | `uno-q/bms_gaussian.c`, `raspberry-pi-5/bms_gaussian.c`, `beaglebone-ai64/bms_gaussian.c`, `bms_handoff.zip` |
| current | `9593c495e740e7f406a6c4a741661180` | `bms_gaussian_c/bms_gaussian.c` |

The only change is the addition of `--volume DIR`, which dumps the full P³
field for the 3D figures. No physics, no timing path, no output format was
touched: the plain run of the current source prints the same `Timing`-and-above
text byte for byte as the stored `run_1.txt`, and the exported CSVs are
identical apart from the known last-digit `energy_J` rounding.

So a folder that carries `b4ca2291...` is a valid result. Its numbers are the
same as a run of the current source. Keep the copy inside each folder as the
record of what actually produced it — don't overwrite it with the newer file.

`bms_handoff.zip` was deliberately left at the revision already sent to the
classmate, so everyone runs the identical source. It will be rebuilt with
`--volume` only if a new board joins after the 3D figures are wanted from it.

## Comparing

```bash
python3 compare_boards.py results/uno-q results/raspberry-pi-5 results/beaglebone-ai64
```

Add `LABEL=` to name a folder that has no `environment.txt`, e.g.
`"BeagleBone=/path/to/folder"`. The script exits non-zero if any physical
quantity differs, because then the boards were not running the same model and
the timing comparison means nothing.
