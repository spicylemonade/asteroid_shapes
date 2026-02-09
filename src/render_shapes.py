#!/usr/bin/env python3
"""
Render 3D figures of asteroid shape models from .obj files.

Reads converged pipeline results and generates:
  - Individual 3-view figures for each asteroid shape
  - A composite gallery figure showing all shapes in a grid

References:
  - Pipeline results: results/pipeline_results.json
  - Spin vectors: results/spin_vectors.csv
  - Shape files: results/shapes/*.obj
"""

import json
import math
import os
import sys

import numpy as np

# Use non-interactive backend before importing pyplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Import load_obj from mesh_utils in the same directory
sys.path.insert(0, os.path.dirname(__file__))
from mesh_utils import load_obj


# ============================================================
# Paths
# ============================================================
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
SHAPES_DIR = os.path.join(REPO_ROOT, "results", "shapes")
SPIN_CSV = os.path.join(REPO_ROOT, "results", "spin_vectors.csv")
PIPELINE_JSON = os.path.join(REPO_ROOT, "results", "pipeline_results.json")
FIGURES_DIR = os.path.join(REPO_ROOT, "figures")


# ============================================================
# Helper: set up a 3D axes for shape rendering
# ============================================================
def _setup_3d_ax(ax, vertices, elev, azim, title=None, fontsize=10):
    """Configure a 3D axes with black background and equal aspect."""
    ax.set_facecolor("black")
    ax.set_axis_off()
    ax.view_init(elev=elev, azim=azim)

    # Equal aspect ratio
    centroid = vertices.mean(axis=0)
    max_range = (vertices.max(axis=0) - vertices.min(axis=0)).max() / 2.0
    ax.set_xlim(centroid[0] - max_range, centroid[0] + max_range)
    ax.set_ylim(centroid[1] - max_range, centroid[1] + max_range)
    ax.set_zlim(centroid[2] - max_range, centroid[2] + max_range)

    if title:
        ax.set_title(title, fontsize=fontsize, color="white", pad=2)


def _add_mesh(ax, vertices, faces, color="0.7", edgecolor="0.4", alpha=0.9):
    """Add a triangular mesh to a 3D axes."""
    tri_verts = vertices[faces]
    mesh = Poly3DCollection(
        tri_verts, alpha=alpha, facecolor=color, edgecolor=edgecolor, linewidths=0.1
    )
    ax.add_collection3d(mesh)


# ============================================================
# Render individual asteroid shape (3 views)
# ============================================================
def render_individual(
    asteroid_id, vertices, faces, period_hours, lam, beta, chi2_reduced, out_path
):
    """Render a single asteroid shape from 3 viewing angles and save to *out_path*."""
    fig = plt.figure(figsize=(12, 4), facecolor="black")

    views = [
        (0, 0, "Equatorial 0\u00b0"),
        (0, 90, "Equatorial 90\u00b0"),
        (90, 0, "Polar"),
    ]

    for idx, (elev, azim, label) in enumerate(views):
        ax = fig.add_subplot(1, 3, idx + 1, projection="3d", facecolor="black")
        _add_mesh(ax, vertices, faces)
        _setup_3d_ax(ax, vertices, elev=elev, azim=azim, title=label, fontsize=10)

    suptitle = (
        f"Asteroid {asteroid_id}   |   "
        f"P = {period_hours:.4f} h   |   "
        f"\u03bb = {lam:.1f}\u00b0, \u03b2 = {beta:.1f}\u00b0   |   "
        f"\u03c7\u00b2_red = {chi2_reduced:.2f}"
    )
    fig.suptitle(suptitle, fontsize=10, color="white", y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=300, facecolor="black", bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Render gallery
# ============================================================
def render_gallery(entries, out_path):
    """Render a composite gallery figure.

    *entries* is a list of dicts with keys:
        asteroid_id, vertices, faces
    """
    n = len(entries)
    if n == 0:
        print("No shapes to render in gallery.")
        return

    # Determine grid dimensions (prefer ~5 columns)
    ncols = min(5, n)
    nrows = math.ceil(n / ncols)

    cell_w = 3.0
    cell_h = 3.0
    fig = plt.figure(
        figsize=(ncols * cell_w, nrows * cell_h + 0.6), facecolor="black"
    )

    for idx, entry in enumerate(entries):
        ax = fig.add_subplot(nrows, ncols, idx + 1, projection="3d", facecolor="black")
        _add_mesh(ax, entry["vertices"], entry["faces"])
        _setup_3d_ax(
            ax,
            entry["vertices"],
            elev=20,
            azim=30,
            title=str(entry["asteroid_id"]),
            fontsize=10,
        )

    fig.suptitle("Asteroid Shape Gallery", fontsize=12, color="white", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=300, facecolor="black", bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)

    # Load pipeline results
    with open(PIPELINE_JSON, "r") as f:
        pipeline = json.load(f)

    # Collect converged entries (status starts with "converged")
    converged = {
        k: v
        for k, v in pipeline.items()
        if v.get("status", "").startswith("converged")
    }
    print(f"Found {len(converged)} converged asteroids in pipeline results.")

    gallery_entries = []

    for ast_id, info in sorted(converged.items(), key=lambda x: int(x[0])):
        shape_file_rel = info.get("shape_file", "")
        shape_file = os.path.join(REPO_ROOT, shape_file_rel)

        if not os.path.isfile(shape_file):
            print(f"  [{ast_id}] Shape file not found, skipping: {shape_file}")
            continue

        try:
            vertices, faces = load_obj(shape_file)
        except Exception as exc:
            print(f"  [{ast_id}] Error loading OBJ: {exc}")
            continue

        if len(vertices) == 0 or len(faces) == 0:
            print(f"  [{ast_id}] Empty mesh, skipping.")
            continue

        period_hours = info.get("period_hours", 0.0)
        lam = info.get("spin_lam_deg", 0.0)
        beta = info.get("spin_beta_deg", 0.0)
        chi2_reduced = info.get("chi2_reduced", 0.0)

        out_path = os.path.join(FIGURES_DIR, f"{ast_id}_shape.png")
        print(f"  [{ast_id}] Rendering individual figure -> {out_path}")
        render_individual(
            ast_id, vertices, faces, period_hours, lam, beta, chi2_reduced, out_path
        )

        gallery_entries.append(
            {"asteroid_id": ast_id, "vertices": vertices, "faces": faces}
        )

    # Gallery
    if gallery_entries:
        gallery_path = os.path.join(FIGURES_DIR, "shape_gallery.png")
        print(f"\nRendering gallery ({len(gallery_entries)} shapes) -> {gallery_path}")
        render_gallery(gallery_entries, gallery_path)
    else:
        print("No shapes rendered; gallery skipped.")

    print("Done.")


if __name__ == "__main__":
    main()
