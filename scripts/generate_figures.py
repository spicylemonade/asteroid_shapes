#!/usr/bin/env python3
"""
Item 023: Generate 3D visualizations of asteroid shapes.

Produce publication-quality 3D rendered figures showing shapes from 3 viewing angles.
Uses matplotlib 3D projection with surface shading.
"""

import sys
import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from lci_engine.forward_model import load_mesh_obj, compute_facet_properties

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')
FIGURES_DIR = os.path.join(os.path.dirname(__file__), '..', 'figures')

os.makedirs(FIGURES_DIR, exist_ok=True)

def log(msg):
    print(msg, flush=True)


def render_mesh(vertices, faces, title, filepath, spin_axis=None,
                elev_azim_pairs=None, figsize=(18, 6)):
    """Render a 3D mesh from multiple viewing angles."""
    if elev_azim_pairs is None:
        elev_azim_pairs = [(20, 30), (20, 120), (80, 0)]  # equatorial x2, polar x1

    normals, areas, centroids = compute_facet_properties(vertices, faces)

    # Normalize to unit scale
    centroid = np.mean(vertices, axis=0)
    v_centered = vertices - centroid
    max_r = np.max(np.linalg.norm(v_centered, axis=1))
    if max_r > 0:
        v_centered = v_centered / max_r

    fig = plt.figure(figsize=figsize)

    for idx, (elev, azim) in enumerate(elev_azim_pairs):
        ax = fig.add_subplot(1, 3, idx + 1, projection='3d')

        # Light source direction
        light_dir = np.array([1, 0.5, 0.8])
        light_dir /= np.linalg.norm(light_dir)

        # Compute shading
        normals_c, _, _ = compute_facet_properties(v_centered, faces)
        shading = np.maximum(normals_c @ light_dir, 0.0)
        # Add ambient
        shading = 0.2 + 0.8 * shading

        # Build polygon collection
        polygons = []
        colors = []
        for fi, face in enumerate(faces):
            polygon = v_centered[face]
            polygons.append(polygon)
            s = shading[fi]
            colors.append((s * 0.7, s * 0.65, s * 0.6, 1.0))

        mesh = Poly3DCollection(polygons, linewidths=0.1, edgecolors=(0.3, 0.3, 0.3, 0.3))
        mesh.set_facecolor(colors)
        ax.add_collection3d(mesh)

        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_zlim(-1.2, 1.2)
        ax.view_init(elev=elev, azim=azim)

        view_labels = ['Equatorial View 1', 'Equatorial View 2', 'Polar View']
        ax.set_title(view_labels[idx], fontsize=11)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')

        # Draw spin axis if provided
        if spin_axis is not None:
            ax.quiver(0, 0, 0, spin_axis[0]*1.3, spin_axis[1]*1.3, spin_axis[2]*1.3,
                     color='red', arrow_length_ratio=0.1, linewidth=2, label='Spin axis')

        ax.set_aspect('equal')

    fig.suptitle(title, fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout()
    plt.savefig(filepath, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    log(f"  Saved: {filepath}")


def main():
    log("Generating 3D shape visualizations...")

    # 1. Ground truth shapes
    gt_shapes = {
        '433_Eros': {'title': '433 Eros - Ground Truth (NEAR Shoemaker)', 'spin': [0, 0, 1]},
        '216_Kleopatra': {'title': '216 Kleopatra - Ground Truth (Radar+AO)', 'spin': [0, 0, 1]},
        '25143_Itokawa': {'title': '25143 Itokawa - Ground Truth (Hayabusa)', 'spin': [0, 0, 1]},
    }

    for key, info in gt_shapes.items():
        obj_path = os.path.join(RESULTS_DIR, f"ground_truth_{key}.obj")
        if os.path.exists(obj_path):
            log(f"\nRendering ground truth: {key}")
            v, f = load_mesh_obj(obj_path)
            filepath = os.path.join(FIGURES_DIR, f"ground_truth_{key}_shape.png")
            render_mesh(v, f, info['title'], filepath, spin_axis=np.array(info['spin']))

    # 2. Blind validation results
    blind_shapes = {
        '433_Eros': 'Recovered Shape (Blind Test)',
        '216_Kleopatra': 'Recovered Shape (Blind Test)',
        '25143_Itokawa': 'Recovered Shape (Blind Test)',
    }

    # Load blind validation results for spin axis info
    blind_results_path = os.path.join(RESULTS_DIR, 'blind_validation_results.json')
    blind_data = {}
    if os.path.exists(blind_results_path):
        with open(blind_results_path) as f:
            blind_data = json.load(f)

    for key, subtitle in blind_shapes.items():
        obj_path = os.path.join(RESULTS_DIR, f"blind_test_{key}.obj")
        if os.path.exists(obj_path):
            log(f"\nRendering blind test: {key}")
            v, f = load_mesh_obj(obj_path)

            # Get recovered spin axis
            spin = np.array([0, 0, 1])
            if key in blind_data and 'recovered' in blind_data[key]:
                r = blind_data[key]['recovered']
                lam = np.radians(r.get('pole_lambda_deg', 0))
                bet = np.radians(r.get('pole_beta_deg', 45))
                spin = np.array([np.cos(bet)*np.cos(lam), np.cos(bet)*np.sin(lam), np.sin(bet)])

            title = f"{key.replace('_', ' ')} - {subtitle}"
            filepath = os.path.join(FIGURES_DIR, f"blind_test_{key}_shape.png")
            render_mesh(v, f, title, filepath, spin_axis=spin)

    # 3. Side-by-side comparison
    for key in ['433_Eros', '216_Kleopatra', '25143_Itokawa']:
        gt_path = os.path.join(RESULTS_DIR, f"ground_truth_{key}.obj")
        bt_path = os.path.join(RESULTS_DIR, f"blind_test_{key}.obj")
        if os.path.exists(gt_path) and os.path.exists(bt_path):
            log(f"\nRendering comparison: {key}")
            gt_v, gt_f = load_mesh_obj(gt_path)
            bt_v, bt_f = load_mesh_obj(bt_path)

            fig = plt.figure(figsize=(14, 6))

            for col, (verts, fcs, label) in enumerate([
                (gt_v, gt_f, 'Ground Truth'),
                (bt_v, bt_f, 'LCI Pipeline')
            ]):
                ax = fig.add_subplot(1, 2, col + 1, projection='3d')

                centroid = np.mean(verts, axis=0)
                v_c = verts - centroid
                max_r = np.max(np.linalg.norm(v_c, axis=1))
                if max_r > 0:
                    v_c /= max_r

                light = np.array([1, 0.5, 0.8])
                light /= np.linalg.norm(light)
                normals_c, _, _ = compute_facet_properties(v_c, fcs)
                shading = 0.2 + 0.8 * np.maximum(normals_c @ light, 0.0)

                polys = [v_c[face] for face in fcs]
                colors = [(s*0.7, s*0.65, s*0.6, 1.0) for s in shading]
                mesh = Poly3DCollection(polys, linewidths=0.1, edgecolors=(0.3,0.3,0.3,0.3))
                mesh.set_facecolor(colors)
                ax.add_collection3d(mesh)

                ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_zlim(-1.2, 1.2)
                ax.view_init(elev=20, azim=30)
                ax.set_title(label, fontsize=12)
                ax.set_aspect('equal')

            fig.suptitle(f"{key.replace('_', ' ')} - Shape Comparison", fontsize=14, fontweight='bold')
            plt.tight_layout()
            filepath = os.path.join(FIGURES_DIR, f"comparison_{key}.png")
            plt.savefig(filepath, dpi=200, bbox_inches='tight', facecolor='white')
            plt.close(fig)
            log(f"  Saved: {filepath}")

    # 4. Candidate shapes (if available)
    candidate_objs = []
    for f in os.listdir(RESULTS_DIR):
        if f.startswith('candidate_') and f.endswith('.obj'):
            candidate_objs.append(f)

    for obj_name in sorted(candidate_objs)[:10]:
        log(f"\nRendering candidate: {obj_name}")
        v, f = load_mesh_obj(os.path.join(RESULTS_DIR, obj_name))
        name = obj_name.replace('.obj', '').replace('candidate_', '').replace('_', ' ')
        filepath = os.path.join(FIGURES_DIR, f"{obj_name.replace('.obj', '_shape.png')}")
        render_mesh(v, f, f"Asteroid {name} - New Shape Model", filepath)

    log(f"\n{'='*60}")
    log(f"Figure generation complete. Files in {FIGURES_DIR}/")
    log(f"Total figures: {len(os.listdir(FIGURES_DIR))}")


if __name__ == '__main__':
    main()
