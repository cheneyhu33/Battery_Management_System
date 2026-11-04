# StarFive VisionFive 2 — RISC-V benchmark

| | |
| --- | --- |
| Board | StarFive VisionFive 2 (JH7110) |
| CPU | 4 × SiFive U74, rv64imafdc, in-order, 1.5 GHz |
| OS | Debian trixie, kernel 6.12.5-starfive |
| Compiler | gcc 12.2.0 (Debian 12.2.0-14), built on the board |
| Memory | 3.8 GiB |

## Runs

Eight plain runs, in two configurations:

| Runs | Governor | Wall (s) | ms/step |
| --- | --- | ---: | ---: |
| 1–3 | `ondemand` (default) | 46.98, 46.01, 42.62 | 9.713, 9.513, 8.812 |
| 4–8 | `performance`, pinned 1.500 GHz | 42.67, 42.83, 43.10, 43.94, 43.69 | 8.821, 8.854, 8.911, 9.085, 9.033 |

Median of the pinned runs: **43.10 s**, **8.91 ms/step**, RTF **9.0**.

The `ondemand` runs are dominated by the ramp from 750 MHz to 1.5 GHz over
the first tens of seconds; the third one, started warm, already matches the
pinned clock. The pinned runs are the steady-state measurement.

Thermals: 60.5 °C at the end of the pinned series. No throttling.

## Identity with the other boards

The export run (`run_export.txt`, 54.2 s) reproduced the cut-off exactly:

```
t = 387.863 s (step 4837)   surface 60.0008 °C   centre 101.7759 °C
stored 632.5409 J
```

`data/midplane_slices.csv` is **byte-identical** to `results/uno-q/data/`
(md5 `ff89056cf33c13ee2fc793499f61017a`), as it already was across the Pi 5,
the BeagleBone and the UNO Q.

`data/temperature_history.csv` differs on two of its 4839 rows, each by one
unit in the last printed digit:

| step | column | UNO Q | VisionFive 2 |
| ---: | --- | --- | --- |
| 2348 | `center_C` | 84.9731332506 | 84.9731332505 |
| 3681 | `energy_J` | 508.6357339214 | 508.6357339213 |

That is a difference of 1e-10 K and 1e-9 J, from summation order, not from a
different solution. The same one-ULP effect was already seen between the UNO Q
and a Mac.
