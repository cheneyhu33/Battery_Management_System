# Comparing the same model on two boards

This archive is a 3D thermal-diffusion model of a single battery cell, written
in C. It runs until the hottest surface node reaches 60 °C, then stops.

The point is simple: **run this exact same source on your board and on ours,
and check that the physics comes out identical and see how the wall time
compares.**

## What you need to do

```sh
sh collect.sh
```

That is the whole job. It compiles, does three plain (timed) runs, one run that
also exports the CSVs, and records your CPU environment. The result is a folder
called `bms_result_<your-hostname>/`. **Send back that whole folder.**

## What you need

- A C compiler. Most boards have one. If not:
  - Debian / Ubuntu / Raspberry Pi OS: `sudo apt install build-essential`
  - Fedora: `sudo dnf install gcc`
- A POSIX shell (`sh` is enough — nothing bash-specific)
- A few hundred KB of disk

Nothing is installed system-wide and nothing outside the current directory is
touched, except for the compiled `bms_gaussian` binary left in the directory
you ran the script from.

## Things not to do

| Don't | Why |
| --- | --- |
| **Don't** edit `bms_gaussian.c` | It is no longer the same model, and the comparison is worthless |
| **Don't** add `-ffast-math` or `-march=native` | They change the floating-point results, so the physics numbers stop matching |
| Don't hand-roll the compile command | `collect.sh` pins `-std=c99 -O2` so every board builds identically |
| Don't run heavy work alongside it | It inflates the timing; the script records the load average so this can be checked |
| Don't Ctrl-C partway | It takes roughly 30 s to 3 min; slower boards take longer |

## What comes back

```
bms_result_<hostname>/
  run_1.txt  run_2.txt  run_3.txt   three plain timed runs, full stdout
  run_export.txt                    the CSV-export run (its time is not compared)
  data/                             the CSVs the heatmaps are drawn from
  environment.txt                   CPU model, cpufreq governor and limits, memory,
                                    compiler version, binary size, load average
  bms_gaussian.c  collect.sh        the source you actually ran, for the record
```

## How we will use it

Every number above the `Timing` line in `run_*.txt` — trip step 4837, surface
60.0008 °C, core 101.7759 °C, energies 775.7258 / 632.5409 / 143.1849 J — must
match **digit for digit**. That is what "the same code" means. Only `wall` and
`per step` are allowed to differ; those two are the board's result.

`data/` goes straight into the same plotting script, so your heatmaps and ours
are produced by one renderer and can be put side by side.
