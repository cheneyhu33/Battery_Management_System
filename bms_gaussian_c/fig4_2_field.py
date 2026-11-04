#!/usr/bin/env python3
"""Figure 4.2 of the EE3180 report: the temperature field at the cutoff, seen
from the outside in.

    (a) the mid-plane z = 10 mm with isotherms in C
    (b) the three mid-planes through the centre, cut away
    (c) planes parallel to the x-y face at z = 0, 2.5, 5, 7.5 and 10 mm on
        one shared colour scale: the 101.8 C core shows up as 60 C on the skin

The source is an isotropic Gaussian in a cubic cell, so the three orthogonal
mid-planes are identical (max |xy - xz| = 0 exactly). The planes in (c) are
the ones that actually differ, and they are the ones a surface sensor sees.

Run the model first:

    ./bms_gaussian --volume results/uno-q/volume

    python3 fig4_2_field.py [VOLUME_DIR] [OUT_PNG]

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
T_AMB = 25.0
VIEW = (22, -50)                # elev, azim, as in plot_heatmap_3d.py
DEPTHS = (0.0, 2.5, 5.0, 7.5, 10.0)

ROOT = Path(__file__).resolve().parent
VOL = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else \
    ROOT / "results" / "uno-q" / "volume"
OUT = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else \
    ROOT / "report_figures" / "fig4_2_field.png"

G = np.linspace(0.0, L_MM, P)   # node positions [mm]


# ---------------------------------------------------------------- data
with open(VOL / "volume_index.csv") as f:
    rows = list(csv.DictReader(f))
if not rows:
    raise SystemExit("no rows in volume_index.csv")
last = rows[-1]
step = int(last["step"])
trip_t = float(last["time_s"])
trip_core = float(last["center_C"])
trip_surf = float(last["surface_C"])

raw = np.fromfile(VOL / last["file"], dtype="<f4")
if raw.size != P * P * P:
    raise SystemExit("%s: %d values, expected %d"
                     % (last["file"], raw.size, P * P * P))
vol = raw.reshape(P, P, P).astype(np.float64) + T_AMB

c = P // 2
mid = vol[:, :, c]              # the z = centre plane, (x, y)

# ---------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.8,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

cmap = cm.inferno
norm = Normalize(vmin=float(vol.min()), vmax=float(vol.max()))

fig = plt.figure(figsize=(6.54, 5.5))
gs = fig.add_gridspec(2, 5, height_ratios=[1.5, 1.0],
                      left=0.075, right=0.965, top=0.955, bottom=0.215,
                      wspace=0.42, hspace=0.30)

# ---- (a) the mid-plane with isotherms ------------------------------------
ax = fig.add_subplot(gs[0, :2])
ax.imshow(mid, origin="lower", extent=(0, L_MM, 0, L_MM), cmap=cmap,
          norm=norm, interpolation="nearest")
levels = np.arange(60.0, 101.0, 10.0)
cs = ax.contour(G, G, mid, levels=levels, colors="white", linewidths=0.7)
ax.clabel(cs, fmt="%d", fontsize=6, inline=True)
ax.set_xlabel("x (mm)")
ax.set_ylabel("y (mm)")
ax.set_xticks([0, 5, 10, 15, 20])
ax.set_yticks([0, 5, 10, 15, 20])
ax.set_title("(a)  mid-plane  z = 10 mm", loc="left", fontsize=8.5, pad=4)

# ---- (b) the three mid-planes, cut away ----------------------------------
ax3 = fig.add_subplot(gs[0, 2:], projection="3d")
X, Y = np.meshgrid(G, G, indexing="ij")
face = cmap(norm(vol))
for surf in (ax3.plot_surface(X, Y, np.full_like(X, G[c]),
                              facecolors=face[:, :, c]),
             ax3.plot_surface(X, np.full_like(X, G[c]), Y,
                              facecolors=face[:, c, :]),
             ax3.plot_surface(np.full_like(X, G[c]), X, Y,
                              facecolors=face[c, :, :])):
    surf.set_edgecolor("none")
    surf.set_antialiased(False)
    surf.set_alpha(1.0)
for a, b in [((0, 0, 0), (20, 0, 0)), ((0, 20, 0), (20, 20, 0)),
             ((0, 0, 20), (20, 0, 20)), ((0, 20, 20), (20, 20, 20)),
             ((0, 0, 0), (0, 20, 0)), ((20, 0, 0), (20, 20, 0)),
             ((0, 0, 20), (0, 20, 20)), ((20, 0, 20), (20, 20, 20)),
             ((0, 0, 0), (0, 0, 20)), ((20, 0, 0), (20, 0, 20)),
             ((0, 20, 0), (0, 20, 20)), ((20, 20, 0), (20, 20, 20))]:
    ax3.plot(*zip(a, b), color="#c8ced6", lw=0.5)
ax3.set_xlim(0, L_MM)
ax3.set_ylim(0, L_MM)
ax3.set_zlim(0, L_MM)
ax3.set_box_aspect((1, 1, 1))
ax3.view_init(*VIEW)
ax3.set_xticks([0, 10, 20])
ax3.set_yticks([0, 10, 20])
ax3.set_zticks([0, 10, 20])
ax3.tick_params(labelsize=6, pad=0)
ax3.set_xlabel("x (mm)", fontsize=7, labelpad=2)
ax3.set_ylabel("y (mm)", fontsize=7, labelpad=2)
ax3.set_zlabel("z (mm)", fontsize=7, labelpad=3)
ax3.set_title("(b)  three mid-planes", loc="left", fontsize=8.5, pad=0)

# ---- (c) planes parallel to the face, outside in -------------------------
ims = []
for n, z in enumerate(DEPTHS):
    axc = fig.add_subplot(gs[1, n])
    k = int(round(z / (L_MM / (P - 1))))
    plane = vol[:, :, k]
    ims.append(axc.imshow(plane, origin="lower", extent=(0, L_MM, 0, L_MM),
                          cmap=cmap, norm=norm, interpolation="nearest"))
    axc.set_xticks([0, 10, 20])
    axc.set_yticks([0, 10, 20])
    axc.tick_params(labelsize=6, pad=1)
    axc.set_xlabel("x (mm)", fontsize=7, labelpad=1)
    if n == 0:
        axc.set_ylabel("y (mm)", fontsize=7, labelpad=1)
    else:
        axc.set_yticklabels([])
    axc.set_title("z = %.1f mm\n%.1f–%.1f °C"
                  % (z, plane.min(), plane.max()),
                  fontsize=6.8, pad=3)
    if n == 0:
        axc.text(-0.34, 1.45, "(c)  planes parallel to the face, outside in",
                 transform=axc.transAxes, fontsize=8.5, va="top")

# one colour bar for (a)-(c): every panel is on the same scale
cax = fig.add_axes((0.33, 0.105, 0.40, 0.019))
cb = fig.colorbar(ims[-1], cax=cax, orientation="horizontal")
cb.set_label("temperature (°C)", fontsize=7.5, labelpad=2)
cb.ax.tick_params(labelsize=6.5, pad=1)

fig.text(0.075, 0.022,
         "t = %.1f s (step %d): surface %.2f °C, core %.2f °C, gap %.1f K — "
         "the core is not where the sensor is."
         % (trip_t, step, trip_surf, trip_core, trip_core - trip_surf),
         fontsize=6.8, color="#3a3a3a")

OUT.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT, dpi=300)
plt.close(fig)

print("trip %.2f s  surface %.2f C  core %.2f C  gap %.2f K"
      % (trip_t, trip_surf, trip_core, trip_core - trip_surf))
for z in DEPTHS:
    k = int(round(z / (L_MM / (P - 1))))
    print("  z=%4.1f mm   %.2f .. %.2f C"
          % (z, vol[:, :, k].min(), vol[:, :, k].max()))
print("written to", OUT)
