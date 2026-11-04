#!/usr/bin/env python3
"""3D heatmaps from the volume dump the model writes with --volume.

The mid-plane CSV is a single slice, so the older figures are 2D. This reads
the full P^3 field instead and draws it in three ways:

    05  the isosurface of "everything hotter than N C" (the hot volume)
    06  three mid-planes at six times (how the hot core grows)
    07  three mid-planes through the centre, the cut-away view

    python3 plot_heatmap_3d.py VOLUME_DIR [OUT_DIR]

VOLUME_DIR is what `./bms_gaussian --volume DIR` writes:
    volume_index.csv     step,time_s,file,center_C,surface_C
    vol_XXXXXX.f32       P^3 little-endian float32, rise above ambient,
                         index order ((i*P)+j)*P+k, i.e. numpy [i,j,k]

Needs numpy, matplotlib and scikit-image (for marching cubes).
"""
from pathlib import Path
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize, to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

P = 41
L_MM = 20.0                     # cell edge
T_AMB = 25.0
G = np.linspace(0.0, L_MM, P)   # node positions in mm
VIEW = (22, -50)                # elev, azim
EVOL_STEPS = (500, 1000, 2000, 3000, 4000)


def load_index(vol_dir):
    with open(vol_dir / "volume_index.csv") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        r["step"] = int(r["step"])
        r["time_s"] = float(r["time_s"])
        r["center_C"] = float(r["center_C"])
        r["surface_C"] = float(r["surface_C"])
    return rows


def load_volume(vol_dir, row):
    raw = np.fromfile(vol_dir / row["file"], dtype="<f4")
    if raw.size != P * P * P:
        raise SystemExit("%s: %d values, expected %d"
                         % (row["file"], raw.size, P * P * P))
    return raw.reshape(P, P, P).astype(np.float64) + T_AMB


def draw_box(ax):
    """The 20 mm cube the field lives in, so the slices have a frame."""
    e = L_MM
    for a, b in [((0, 0, 0), (e, 0, 0)), ((0, e, 0), (e, e, 0)),
                 ((0, 0, e), (e, 0, e)), ((0, e, e), (e, e, e)),
                 ((0, 0, 0), (0, e, 0)), ((e, 0, 0), (e, e, 0)),
                 ((0, 0, e), (0, e, e)), ((e, 0, e), (e, e, e)),
                 ((0, 0, 0), (0, 0, e)), ((e, 0, 0), (e, 0, e)),
                 ((0, e, 0), (0, e, e)), ((e, e, 0), (e, e, e))]:
        ax.plot(*zip(a, b), color="#9aa7b8", lw=0.7, alpha=0.9)


def plane(ax, facecolors, kind, offset):
    """One coloured mid-plane. shade=False keeps the colour literal: with
    shading on, matplotlib darkens faces by their normal and the same
    temperature reads as a different colour on each plane."""
    X, Y = np.meshgrid(G, G, indexing="ij")
    if kind == "yz":            # constant x, axes (y, z)
        surf = ax.plot_surface(np.full_like(X, offset), X, Y,
                               facecolors=facecolors)
    elif kind == "xz":          # constant y, axes (x, z)
        surf = ax.plot_surface(X, np.full_like(X, offset), Y,
                               facecolors=facecolors)
    else:                       # "xy", constant z, axes (x, y)
        surf = ax.plot_surface(X, Y, np.full_like(X, offset),
                               facecolors=facecolors)
    surf.set_edgecolor("none")
    surf.set_antialiased(False)
    surf.set_alpha(1.0)
    return surf


def draw_isosurface(ax, vol, level, color):
    """One opaque shell: everything hotter than `level`. A single surface has
    no depth-sorting problem, unlike nested translucent ones."""
    from skimage import measure
    if not (vol.min() < level < vol.max() - 1.0):
        return False
    verts, faces, _, _ = measure.marching_cubes(vol, level=level)
    # One RGBA row per face: matplotlib's shading needs a 2-D array here and
    # raises AxisError if handed a single colour tuple.
    facecolors = np.empty((len(faces), 4))
    facecolors[:] = to_rgba(color)
    # Poly3DCollection wants one entry per polygon, so index with faces:
    # verts alone is a per-vertex list and makes the shading code fail.
    mesh = Poly3DCollection(verts[faces] * (L_MM / (P - 1)), facecolors=facecolors,
                            edgecolor="none", alpha=1.0, shade=True)
    ax.add_collection3d(mesh)
    return True


def draw_midplanes(ax, vol, cmap, norm):
    """The three planes through the centre, which is what a cut-away shows."""
    c = P // 2
    face = cmap(norm(vol))
    plane(ax, face[:, c, :], "xz", G[c])
    plane(ax, face[c, :, :], "yz", G[c])
    plane(ax, face[:, :, c], "xy", G[c])


def style_3d(ax, title, labels=True):
    ax.set_xlim(0, L_MM)
    ax.set_ylim(0, L_MM)
    ax.set_zlim(0, L_MM)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(*VIEW)
    if labels:
        ax.set_xlabel("x (mm)", fontsize=8, labelpad=6)
        ax.set_ylabel("y (mm)", fontsize=8, labelpad=6)
        ax.set_zlabel("z (mm)", fontsize=8, labelpad=8)
    for setter in (ax.set_xticks, ax.set_yticks, ax.set_zticks):
        setter([0, 10, 20])
    ax.tick_params(labelsize=7, pad=1)
    ax.set_title(title, fontsize=10, pad=0)


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    vol_dir = Path(argv[1]).resolve()
    out_dir = (Path(argv[2]).resolve() if len(argv) > 2
               else vol_dir.parent / "plots_3d")
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_index(vol_dir)
    if not rows:
        raise SystemExit("no rows in volume_index.csv")
    last = rows[-1]

    cmap = cm.inferno
    # Two scales: the final-state figures use the range the field actually
    # occupies, which is what makes the structure visible. The evolution
    # panels share one global scale so their colours stay comparable.
    global_norm = Normalize(vmin=T_AMB, vmax=max(r["center_C"] for r in rows))

    # ---- 1. the frozen final state as a stack of slices ------------------
    vol = load_volume(vol_dir, last)
    norm = Normalize(vmin=float(vol.min()) - 1.0, vmax=float(vol.max()) + 0.5)
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    # A level a little above the skin: the shell is the part of the cell that
    # is dangerously hot, and it sits clearly inside the 20 mm cube.
    level = float(np.floor((last["surface_C"] + last["center_C"]) / 20.0) * 10.0)
    fig = plt.figure(figsize=(7.4, 7.0), facecolor="#f6f8fc")
    ax = fig.add_subplot(111, projection="3d")
    ok = draw_isosurface(ax, vol, level, cmap(global_norm(level)))
    draw_box(ax)
    style_3d(ax, "t = %.1f s   (step %d)\ncore %.1f C, skin %.1f C"
             % (last["time_s"], last["step"], last["center_C"], last["surface_C"]))
    # How big the hot body actually is, measured from the grid rather than
    # eyeballed: the bounding box of the nodes at or above `level`.
    above = np.argwhere(vol >= level)
    span = (above.max(0) - above.min(0)) * (L_MM / (P - 1)) if len(above) else None
    ax.text2D(0.02, 0.94,
              "%d C isosurface%s:\neverything inside it is hotter"
              % (level, "  (%.0f mm across)" % span.max() if span is not None else ""),
              transform=ax.transAxes, fontsize=9, color="#334155")
    fig.tight_layout()
    fig.savefig(out_dir / "05_3d_isosurface_final.png", dpi=210)
    fig.savefig(out_dir / "05_3d_isosurface_final.pdf")
    plt.close(fig)
    print("05_3d_isosurface_final: %d C%s" % (level, "" if ok else " (no surface)"))

    # ---- 2. how the hot core grows, one colour scale throughout ----------
    picks = [r for r in rows if r["step"] in EVOL_STEPS]
    if last not in picks:
        picks.append(last)
    fig = plt.figure(figsize=(11.6, 7.6), facecolor="#f6f8fc")
    for n_ax, row in enumerate(picks, start=1):
        ax = fig.add_subplot(2, 3, n_ax, projection="3d")
        draw_midplanes(ax, load_volume(vol_dir, row), cmap, global_norm)
        draw_box(ax)
        style_3d(ax, "t = %.0f s   core %.1f C"
                 % (row["time_s"], row["center_C"]), labels=False)
    fig.suptitle("Growth of the hot core   (mid-planes, same temperature scale everywhere)",
                 fontsize=11, y=0.97)
    sm = cm.ScalarMappable(norm=global_norm, cmap=cmap)
    cax = fig.add_axes((0.905, 0.18, 0.014, 0.62))
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("temperature (C)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.subplots_adjust(left=0.03, right=0.88, top=0.88, bottom=0.04,
                        wspace=0.02, hspace=0.16)
    fig.savefig(out_dir / "06_3d_evolution.png", dpi=210,
                bbox_inches="tight")
    fig.savefig(out_dir / "06_3d_evolution.pdf", bbox_inches="tight")
    plt.close(fig)
    print("06_3d_evolution: %d panels" % len(picks))

    # ---- 3. three orthogonal mid-planes through the centre ---------------
    fig = plt.figure(figsize=(7.4, 7.0), facecolor="#f6f8fc")
    ax = fig.add_subplot(111, projection="3d")
    draw_midplanes(ax, vol, cmap, norm)
    draw_box(ax)
    style_3d(ax, "Mid-planes through the centre at t = %.1f s" % last["time_s"])
    cb = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cb.set_label("temperature (C)", fontsize=9)
    cb.ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "07_3d_midplanes.png", dpi=210)
    fig.savefig(out_dir / "07_3d_midplanes.pdf")
    plt.close(fig)
    print("07_3d_midplanes")

    print("written to", out_dir)


if __name__ == "__main__":
    main(sys.argv)
