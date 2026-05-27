#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Debug visualization for DiMES mesh components.

Usage
-----
Point BASE_PATH at the stem of the saved mesh files, e.g.:

    BASE_PATH = "/home/user/mesh_repo/DiMES_3D"

The script looks for files named  <BASE_PATH>_<component>.npz  and:
  - prints a summary table  (component, nodes, triangles, centroid, area)
  - shows a colour-coded 3-D plot with component labels at their centroids

Requirements: numpy, matplotlib
"""

import glob
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# ── configuration ─────────────────────────────────────────────────────────────
BASE_PATH = "/home/cappellil/GITRmProject/mesh_repo/DiMES_3D"   # edit this
ALPHA     = 0.55          # surface transparency
MAX_TRIS  = 4000          # downsample for rendering speed (None = show all)
# ──────────────────────────────────────────────────────────────────────────────


def load_components(base_path):
    pattern = f"{base_path}_*.npz"
    files = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"No component files found matching: {pattern}")

    components = {}
    for f in files:
        stem = os.path.basename(f)[len(os.path.basename(base_path)) + 1:-4]
        data = np.load(f, allow_pickle=True)
        components[stem] = {
            "nodes":     data["nodes"],
            "triangles": data["triangles"].astype(np.int64),
        }
    return components


def _tri_areas(nodes, tris):
    if len(tris) == 0:
        return np.array([])
    v0 = nodes[tris[:, 0]]
    v1 = nodes[tris[:, 1]]
    v2 = nodes[tris[:, 2]]
    cross = np.cross(v1 - v0, v2 - v0)
    return 0.5 * np.linalg.norm(cross, axis=1)


def component_stats(nodes, tris):
    centroid = nodes.mean(axis=0) if len(nodes) else np.zeros(3)
    area = _tri_areas(nodes, tris).sum()
    return centroid, area


def print_summary(components):
    hdr = f"{'Component':<28} {'Nodes':>7} {'Tris':>7}   {'cx':>8} {'cy':>8} {'cz':>8}   {'Area':>10}"
    print("\n" + hdr)
    print("─" * len(hdr))
    for name, d in components.items():
        c, a = component_stats(d["nodes"], d["triangles"])
        print(f"{name:<28} {len(d['nodes']):>7} {len(d['triangles']):>7}   "
              f"{c[0]:8.3f} {c[1]:8.3f} {c[2]:8.3f}   {a:10.4f}")
    print()


def visualize(components, alpha=ALPHA, max_tris=MAX_TRIS):
    cmap   = cm.get_cmap("tab20", max(len(components), 1))
    fig    = plt.figure(figsize=(13, 8))
    ax     = fig.add_subplot(111, projection="3d")
    legend = []

    all_nodes = np.vstack([d["nodes"] for d in components.values()
                           if len(d["nodes"]) > 0])

    for i, (name, d) in enumerate(components.items()):
        nodes = d["nodes"]
        tris  = d["triangles"]
        color = cmap(i)

        if len(tris) == 0:
            continue

        # optional downsampling for rendering speed
        if max_tris is not None and len(tris) > max_tris:
            idx  = np.random.choice(len(tris), max_tris, replace=False)
            tris = tris[idx]

        polys = nodes[tris]
        col   = Poly3DCollection(polys, alpha=alpha,
                                 facecolor=color, edgecolor="none")
        ax.add_collection3d(col)
        legend.append(plt.Rectangle((0, 0), 1, 1, fc=color, label=name))

        centroid, _ = component_stats(nodes, d["triangles"])
        ax.text(centroid[0], centroid[1], centroid[2],
                name, fontsize=6, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", alpha=0.6, lw=0))

    # axis limits
    mn, mx = all_nodes.min(axis=0), all_nodes.max(axis=0)
    pad = 0.05 * (mx - mn).max()
    ax.set_xlim(mn[0] - pad, mx[0] + pad)
    ax.set_ylim(mn[1] - pad, mx[1] + pad)
    ax.set_zlim(mn[2] - pad, mx[2] + pad)

    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    ax.set_title(f"DiMES mesh components  ({len(components)} surfaces)")
    ax.legend(handles=legend, loc="upper right", fontsize=6,
              ncol=max(1, len(components) // 12))
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        BASE_PATH = sys.argv[1]

    components = load_components(BASE_PATH)
    print_summary(components)
    visualize(components)
