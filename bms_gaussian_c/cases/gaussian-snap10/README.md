# Dense-snapshot case — the same run, 40 field frames instead of 8

The Gaussian case of `results/uno-q/`, unchanged in every physical parameter,
exported with a finer **frame** cadence. The time step `dt` is not touched: it
stays at 0.90 × the stability limit (0.080187 s), because `dt` is a numerical
knob and changing it adds no data.

What changes is how often the field is stored. The default export writes a
frame only at the seven steps in `SNAP` (`0, 100, 500, 1000, 2000, 3000,
4000`), which is 8 frames over the 387.9 s run — the early part of the
transient, where most of the structure forms, gets one frame at 8 s and the
next at 40 s. `--snap-dt` stores every 10 s instead:

| | `results/uno-q/` | here |
| --- | --- | --- |
| field frames | 8 | **40** |
| step between frames | 100 → 1000 (uneven) | 125 (even) |
| cadence in seconds | 8, 40, 80, 160, … | 10.02 s throughout |
| `midplane_slices.csv` | 321 kB | 2.7 MB |
| `volume/` | 2.2 MB | 11 MB |

## How it was produced

```bash
cc -std=c99 -O2 -o bms_gaussian bms_gaussian.c -lm
./bms_gaussian --snap-dt 10 \
    --export cases/gaussian-snap10/data \
    --volume cases/gaussian-snap10/volume
```

`run.txt` is that run's full log. `--snap-dt` divides the requested cadence by
`dt` and rounds to whole steps, so 10 s becomes 125 steps (10.02 s); a cadence
shorter than one step stores every step. Without `--snap-dt` the program uses
the old `SNAP` table, and the default export is unchanged — `midplane_slices.csv`
for the stock run is byte-identical to `results/uno-q/data/`.

The trip frame is always written whatever the cadence, so this export ends with
`t = 387.863 s, core 101.776 °C, skin 60.001 °C` — the same numbers as
`results/uno-q/`.

Provenance: host build on the Mac. The model is board-independent; a plain run
prints the same temperatures as every board's stored `run_*.txt`. The history
CSV differs from the board's copy in one last digit on one row
(`energy_J`, 1 ULP), which is compiler rounding, not a model difference.
