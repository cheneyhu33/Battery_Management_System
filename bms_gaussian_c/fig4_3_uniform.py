#!/usr/bin/env python3
"""Figure 4.3 of the EE3180 report: the uniform-source case, drawn in the same
style as Fig. 4.1 (Section 4.2).

    (a) core and hottest-surface temperature up to the 60 C cutoff
    (b) temperature along the centreline at the snapshot steps

The same two panels, the same colours and the same layout as Fig. 4.1, so the
two figures can be read side by side: there the core runs 41.8 K ahead of the
sensor, here only 1.6 K.

Data come from bms_gaussian.c run with --uniform --export:

    ./bms_gaussian --uniform --export cases/uniform-source/data

    python3 fig4_3_uniform.py [DATA_DIR] [OUT_PNG]

Needs numpy and matplotlib only.
"""
from pathlib import Path
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DATA = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    ROOT / "cases" / "uniform-source" / "data"
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else \
    ROOT / "report_figures" / "fig4_3_uniform.png"

P = 41                      # nodes per axis
L_MM = 20.0                 # cell edge [mm]
CORE = "#b2182b"            # the core colour of Fig. 4.1
SURF = "#2166ac"            # the surface colour of Fig. 4.1
GATE = "#5f6b76"            # the 60 C cutoff line
GRID = "#ebebeb"

# ---------------------------------------------------------------- data
t, core, surf = [], [], []
with open(DATA / "temperature_history.csv") as f:
    for row in csv.DictReader(f):
        t.append(float(row["time_s"]))
        core.append(float(row["center_C"]))
        surf.append(float(row["surface_max_C"]))
t, core, surf = np.array(t), np.array(core), np.array(surf)

slices = {}
with open(DATA / "midplane_slices.csv") as f:
    for row in csv.DictReader(f):
        n = int(row["step"])
        if n not in slices:
            slices[n] = np.full((P, P), np.nan)
        slices[n][int(row["i"]), int(row["j"])] = float(row["T_C"])

steps = sorted(slices)[1:]      # step 0 is the uniform 25 C start: Fig. 4.1(b)
                                # does not plot it either
trip = steps[-1]
trip_t, trip_core, trip_surf = t[-1], core[-1], surf[-1]

x = np.linspace(0.0, L_MM, P)

# ---------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.66, 4.65))

# ---- (a) core and surface up to the cutoff -------------------------------
ax1.plot(t, core, color=CORE, lw=2.6, solid_capstyle="round",
         label="core (centre node)")
ax1.plot(t, surf, color=SURF, lw=2.6, solid_capstyle="round",
         label="hottest surface node")
ax1.axhline(60.0, color=GATE, lw=1.3, ls=(0, (5, 4)))
ax1.text(4.0, 61.0, "60 °C cutoff", color="#3f4a55", fontsize=16, va="bottom")

ax1.set_xlim(0, trip_t * 1.02)
ax1.set_ylim(24, 66)
ax1.set_xlabel("time (s)")
ax1.set_ylabel("temperature (°C)")
ax1.grid(color=GRID, lw=1.0)
ax1.set_axisbelow(True)
ax1.legend(loc="lower right", frameon=False, handlelength=2.0,
           borderaxespad=1.0)
ax1.set_title("(a)", loc="left", fontsize=18)

# ---- (b) centreline at the snapshot steps --------------------------------
colors = plt.cm.viridis(np.linspace(.05, .88, len(steps)))
for n, c in zip(steps, colors):
    ax2.plot(x, slices[n][P // 2, :], color=c, lw=2.4,
             label="{:.0f} s".format(t[n]))

ax2.set_xlim(0, L_MM)
ax2.set_ylim(20, 100)       # headroom for the legend, as in Fig. 4.1(b)
ax2.set_xlabel("x along centreline (mm)")
ax2.set_ylabel("temperature (°C)")
ax2.set_xticks([0, 5, 10, 15, 20])
ax2.grid(color=GRID, lw=1.0)
ax2.set_axisbelow(True)
ax2.legend(loc="upper left", frameon=False, ncol=2, handlelength=2.0,
           columnspacing=1.6, borderaxespad=0.6)
ax2.set_title("(b)", loc="left", fontsize=18)

fig.subplots_adjust(left=0.085, right=0.985, top=0.86, bottom=0.175,
                    wspace=0.22)

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=150)
plt.close(fig)

print("trip  %.2f s   surface %.2f C   core %.2f C   gap %.2f K"
      % (trip_t, trip_surf, trip_core, trip_core - trip_surf))
print("written to", OUT)
