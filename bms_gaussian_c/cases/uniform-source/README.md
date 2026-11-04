# Uniform-source case — Section 4.4, Fig. 4.3

The same cell, the same 2 W, the same grid, time step and stop criterion as
`results/uno-q/`; only the shape of the source changes. Instead of the central
Gaussian, `q'''` is spread evenly through the cell:

```
q'''(r) = Q_total / V = 2.5e5 W/m^3   (every node, every step)
```

This is a **parameter case, not a board**: it is not a run of a different
machine, so it lives here rather than under `results/` (which is one folder per
board).

## How it was produced

Host build on the Mac — the model is board-independent, and its plain run
prints the same numbers as every board's stored `run_*.txt`.

```bash
cc -std=c99 -O2 -o bms_gaussian bms_gaussian.c -lm
./bms_gaussian --uniform --export cases/uniform-source/data
```

`run.txt` is that run's full log. The source is `bms_gaussian.c` at the top of
this project with `--uniform` added (the default path is untouched).

## Result

| Quantity | Uniform | Gaussian (`results/uno-q`) |
| --- | --- | --- |
| cutoff | 402.778 s (step 5023) | 387.863 s (step 4837) |
| core at cutoff | 61.607 °C | 101.776 °C |
| surface at cutoff | 60.005 °C | 60.001 °C |
| core – surface | **1.602 K** | **41.775 K** |
| energy stored | 626.7 J | 632.5 J |
| energy generated | 805.6 J | 775.7 J |

These are the numbers in Table 4.2 of the report (402.8 s, 61.6 °C, 1.6 K),
recomputed from this run.

## Figure

```bash
python3 fig4_3_uniform.py        # -> report_figures/fig4_3_uniform.png
```

Fig. 4.3 is the same two panels as Fig. 4.1 — (a) core and hottest-surface
temperature up to the cutoff, (b) centreline profile at the snapshot steps —
drawn with the same colours, fonts and layout, so the two figures can be read
side by side.
