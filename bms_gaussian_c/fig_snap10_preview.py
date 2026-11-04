#!/usr/bin/env python3
"""Preview of what the denser export buys: the mid-plane every 40 s.

Reads cases/gaussian-snap10/data (40 frames at 10.02 s) and draws ten of them
on one colour scale, from 0 s to the trip. The old export had eight frames in
total, at 8, 40, 80, 160, 240, 320 and 388 s.

    python3 fig_snap10_preview.py [DATA_DIR] [OUT_PNG]

Needs numpy and matplotlib only.
"""
from pathlib import Path
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize

P = 41
L_MM = 20.0
EVERY_S = 40.0

ROOT = Path(__file__).resolve().parent
DATA = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    ROOT / "cases" / "gaussian-snap10" / "data"
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else \
    ROOT / "cases" / "gaussian-snap10" / "preview.png"

hist = np.genfromtxt(DATA / "temperature_history.csv", delimiter=",", names=True)
sl = np.genfromtxt(DATA / "midplane_slices.csv", delimiter=",", names=True)

steps = np.unique(sl["step"])
times = {int(s): float(hist["time_s"][np.where(hist["step"] == s)[0][0]])
         for s in steps}

# The 0 s frame is the whole cell at ambient: one flat colour, no information.
# Start at the first frame with structure and end on the cut-off frame.
trip = int(steps[-1])
picks, last_t = [], 0.0
for s in steps[:-1]:
    if times[int(s)] - last_t >= EVERY_S - 1e-9:
        picks.append(int(s))
        last_t = times[int(s)]
picks = picks[:9] + [trip]

frames = {s: sl[sl["step"] == s]["T_C"].reshape(P, P) for s in picks}
vmin = min(f.min() for f in frames.values())
vmax = hist["center_C"].max()

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8.5,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

fig, axes = plt.subplots(2, 5, figsize=(6.54, 2.95))
fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.30,
                    wspace=0.30, hspace=0.46)
norm = Normalize(vmin=float(vmin), vmax=float(vmax))
for n, (ax, s) in enumerate(zip(axes.flat, picks)):
    im = ax.imshow(frames[s], origin="lower", extent=(0, L_MM, 0, L_MM),
                   cmap=cm.inferno, norm=norm, interpolation="nearest")
    ax.set_xticks([0, 10, 20])
    ax.set_yticks([0, 10, 20])
    ax.tick_params(labelsize=6, pad=1)
    if n >= 5:
        ax.set_xlabel("x (mm)", fontsize=7, labelpad=1)
    else:
        ax.set_xticklabels([])
    if n % 5 == 0:
        ax.set_ylabel("y (mm)", fontsize=7, labelpad=1)
    else:
        ax.set_yticklabels([])
    ax.set_title("%.0f s" % times[s], fontsize=7.5, pad=2)

cax = fig.add_axes((0.35, 0.155, 0.33, 0.024))
cb = fig.colorbar(im, cax=cax, orientation="horizontal")
cb.set_label("temperature (°C)", fontsize=7, labelpad=2)
cb.ax.tick_params(labelsize=6, pad=1)

fig.text(0.075, 0.045,
         "Mid-plane z = 10 mm, one frame every 40 s, ending at the cut-off. "
         "40 frames at 10.02 s; the old table gave 8, unevenly spaced.",
         fontsize=6.8, color="#3a3a3a")

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=300)
plt.close(fig)
print("frames used:", [(s, round(times[s], 1)) for s in picks])
print("written to", OUT)
