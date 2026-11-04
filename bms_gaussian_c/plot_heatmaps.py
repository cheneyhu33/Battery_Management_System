#!/usr/bin/env python3
"""Render the heatmaps from the CSVs the C model exported.

This script knows nothing about the board it ran on: it reads one data/ folder
and writes one plots/ folder. The same CSVs come out of every board that runs
bms_gaussian.c, so the same figures can be produced for each board and compared
side by side.

Run: python3 plot_heatmaps.py                          (results/uno-q: data -> plots)
     python3 plot_heatmaps.py DIR DIR/plots            (e.g. a returned result folder)
(needs numpy and matplotlib)
"""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
DEFAULT_BOARD = ROOT / "results" / "uno-q"
DATA = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_BOARD / "data"
PLOTS = Path(sys.argv[2]).resolve() if len(sys.argv) > 2 else DEFAULT_BOARD / "plots"
PLOTS.mkdir(parents=True, exist_ok=True)

hist = np.genfromtxt(DATA / "temperature_history.csv", delimiter=",", names=True)
sl = np.genfromtxt(DATA / "midplane_slices.csv", delimiter=",", names=True)

P = 41                       # nodes per axis
L_MM = 20.0                  # cell edge [mm]
steps = np.unique(sl["step"])

core = hist["center_C"]
surf = hist["surface_max_C"]
trip_step = int(hist["step"][-1])
trip_t = float(hist["time_s"][-1])
trip_core = float(hist["center_C"][-1])
trip_surf = float(hist["surface_max_C"][-1])

assert hist["step"][0] == 0 and hist["step"][-1] == trip_step
assert np.isfinite(core).all() and np.isfinite(surf).all()
assert np.all(core >= surf - 1e-9), "core must not be cooler than the surface"
assert np.all(np.diff(core) > 0), "core temperature must rise monotonically"
assert trip_surf >= 60.0, "the run must stop at or above the limit"

def plane(n):
    return sl[sl["step"] == n]["T_C"].reshape(P, P)

# Every slice must be symmetric about both centre lines.
for n in steps:
    M = plane(n)
    assert np.max(np.abs(M - M.T)) < 1e-9
    assert np.max(np.abs(M - M[::-1, :])) < 1e-9

meta = {
    "source": "bms_gaussian.c --export, run on the board",
    "temperature_history_sha256": hashlib.sha256(
        (DATA / "temperature_history.csv").read_bytes()).hexdigest(),
    "midplane_slices_sha256": hashlib.sha256(
        (DATA / "midplane_slices.csv").read_bytes()).hexdigest(),
    "trip_step": trip_step,
    "trip_time_s": trip_t,
    "surface_at_trip_C": trip_surf,
    "core_at_trip_C": trip_core,
    "core_minus_surface_K": trip_core - trip_surf,
    "peak_C": float(hist["max_C"].max()),
    "min_C_at_trip": float(hist["min_C"][-1]),
    "energy_stored_J": float(hist["energy_J"][-1]),
    "energy_generated_J": 2.0 * trip_t,
    "snapshot_steps": [int(s) for s in steps],
    "cell_mm": L_MM,
    "nodes": P,
    "validation": "finite; core >= surface; core rises monotonically; "
                  "every slice symmetric under both axes; final surface >= 60 C",
}
try:
    (DATA / "validation.json").write_text(json.dumps(meta, indent=2) + "\n")
except OSError:
    # A returned result folder may be read-only; the figures are what matter.
    (PLOTS / "validation.json").write_text(json.dumps(meta, indent=2) + "\n")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
    "axes.labelsize": 10, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "#f6f8fc", "axes.facecolor": "white",
    "text.color": "#182738", "axes.labelcolor": "#182738",
    "xtick.color": "#415469", "ytick.color": "#415469",
    "savefig.facecolor": "#f6f8fc", "pdf.fonttype": 42,
})
FOOTER = ("Central Gaussian source (sigma = 2 mm, 2 W) in a 20 mm cell; "
          "mid-plane z = 10 mm. Gate: the surface reaches 60 C.")


def save(fig, name):
    fig.savefig(PLOTS / (name + ".png"), dpi=210)
    fig.savefig(PLOTS / (name + ".pdf"))
    plt.close(fig)


extent = (0, L_MM, 0, L_MM)

# ---- 1. evolution on a shared scale ---------------------------------------
panels = [s for s in steps if s in (0, 100, 500, 1000, 2000, 3000, 4000, trip_step)]
panels = panels[:8] if len(panels) > 8 else panels
fig, axes = plt.subplots(2, 4, figsize=(15, 7.6), layout="constrained")
fig.suptitle("Heat builds up in the core while the surface lags behind",
             fontsize=17, fontweight="bold", y=.98)
for ax, n in zip(axes.flat, panels):
    img = ax.imshow(plane(n), origin="lower", extent=extent,
                    vmin=25.0, vmax=trip_core, cmap="inferno",
                    interpolation="nearest")
    row = np.where(hist["step"] == n)[0][0]
    ax.set_title("step {}\n{:.1f} s   core {:.2f} C".format(
        int(n), float(hist["time_s"][row]), float(hist["center_C"][row])), fontsize=10)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
fig.colorbar(img, ax=list(axes.flat), label="Temperature (C), shared scale",
             shrink=.86, pad=.02)
fig.get_layout_engine().set(rect=(0, .07, 1, .80))
fig.text(.02, .02, FOOTER, fontsize=9, color="#526479")
save(fig, "01_temperature_evolution")

# ---- 2. final field, zoomed, plus its centreline ---------------------------
fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.6), layout="constrained")
fig.suptitle("At the moment the sensor trips", fontsize=17, fontweight="bold", y=.98)
last = plane(trip_step)
img = axes[0].imshow(last, origin="lower", extent=extent,
                     vmin=float(last.min()), vmax=float(last.max()),
                     cmap="inferno", interpolation="nearest")
levels = np.linspace(float(last.min()), float(last.max()), 8)[1:-1]
axes[0].contour(np.linspace(0, L_MM, P), np.linspace(0, L_MM, P), last,
                levels=levels, colors="white", alpha=.55, linewidths=.6)
axes[0].set(xlabel="x (mm)", ylabel="y (mm)",
            title="Mid-plane, {:.2f} to {:.2f} C".format(last.min(), last.max()))
fig.colorbar(img, ax=axes[0], label="Temperature (C)",
             orientation="horizontal", shrink=.86, pad=.20)

x = np.linspace(0, L_MM, P)
axes[1].plot(x, last[P // 2, :], color="#b64614", lw=2.6)
axes[1].axhline(60.0, color="#137b8a", ls="--", lw=1.4)
axes[1].text(1.0, 61.5, "60 C gate", color="#137b8a", fontsize=9)
axes[1].scatter([L_MM / 2], [trip_core], color="#182738", s=42, zorder=3)
axes[1].annotate("core {:.2f} C".format(trip_core),
                 xy=(L_MM / 2, trip_core), xytext=(.55, .86),
                 textcoords="axes fraction",
                 arrowprops={"arrowstyle": "-", "color": "#526479"}, fontsize=11)
axes[1].set(xlabel="x (mm)", ylabel="Temperature (C)",
            title="Centreline through y = z = 10 mm",
            xlim=(0, L_MM), ylim=(float(last.min()) - 2, trip_core + 3))
axes[1].grid(alpha=.18)
fig.get_layout_engine().set(rect=(0, .12, 1, .74))
fig.text(.02, .03,
         "The surface reads {:.2f} C while the core is {:.2f} C: a gap of "
         "{:.1f} K.\n".format(trip_surf, trip_core, trip_core - trip_surf) + FOOTER,
         fontsize=9, color="#526479")
save(fig, "02_final_field_and_profile")

# ---- 3. history -----------------------------------------------------------
fig, ax = plt.subplots(figsize=(10.4, 5.8), layout="constrained")
fig.suptitle("Surface sensor versus core, up to the trip",
             fontsize=17, fontweight="bold", y=.98)
t = hist["time_s"]
ax.plot(t, core, lw=2.6, color="#b64614", label="Core (centre node)")
ax.plot(t, surf, lw=2.6, color="#137b8a", label="Hottest surface node (the sensor)")
ax.axhline(60.0, color="#7a8699", ls="--", lw=1.2)
ax.axvline(trip_t, color="#7a8699", ls=":", lw=1.2)
ax.scatter([trip_t, trip_t], [trip_surf, trip_core], color="#182738", s=45, zorder=3)
ax.annotate("trip  {:.2f} s\nsurface {:.2f} C\ncore {:.2f} C".format(
    trip_t, trip_surf, trip_core), xy=(trip_t, trip_core), xytext=(.42, .30),
    textcoords="axes fraction",
    arrowprops={"arrowstyle": "-", "color": "#526479"}, fontsize=10.5)
ax.set(xlabel="Simulated time (s)", ylabel="Temperature (C)", xlim=(0, trip_t))
ax.grid(alpha=.18)
ax.legend(frameon=False, loc="upper left")
fig.get_layout_engine().set(rect=(0, .11, 1, .75))
fig.text(.02, .03,
         "The gap is the sensor lag: {:.1f} K at the trip, and it is still "
         "widening.\n".format(trip_core - trip_surf) + FOOTER,
         fontsize=9, color="#526479")
save(fig, "03_core_vs_surface")

# ---- 4. centreline evolution ---------------------------------------------
fig, ax = plt.subplots(figsize=(9.8, 6.0), layout="constrained")
fig.suptitle("The core heats up first and the profile stays steep",
             fontsize=17, fontweight="bold", y=.98)
colors = plt.cm.viridis(np.linspace(.05, .88, len(panels)))
for n, c in zip(panels, colors):
    row = np.where(hist["step"] == n)[0][0]
    ax.plot(x, plane(n)[P // 2, :], lw=2.1, color=c,
            label="{:.0f} s (step {})".format(float(hist["time_s"][row]), int(n)))
ax.set(xlabel="x (mm)", ylabel="Temperature (C)", xlim=(0, L_MM))
ax.grid(alpha=.18)
ax.legend(frameon=False, fontsize=9)
fig.get_layout_engine().set(rect=(0, .11, 1, .74))
fig.text(.02, .03, FOOTER, fontsize=9, color="#526479")
save(fig, "04_centreline_evolution")

print(json.dumps(meta, indent=2))
print("Plots generated:", *[f.name for f in sorted(PLOTS.glob("*.png"))], sep="\n")
