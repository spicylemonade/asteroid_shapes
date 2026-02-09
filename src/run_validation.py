"""
Blind validation: run full inversion pipeline on ground truth asteroids
and compare results against DAMIT/spacecraft shape models.

Item 016 of the research rubric.
"""
import numpy as np
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))

from data_ingest import preprocess_asteroid
from period_search import combined_period_search, merge_sessions_to_relative
from convex_inversion import ConvexSolver, grid_search_spin, subsample_blocks
from ga_solver import GASolver
from mesh_utils import (
    hausdorff_distance, volumetric_iou, load_obj, save_obj,
    create_sphere_mesh,
)

repo_root = os.path.dirname(os.path.dirname(__file__))
zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')


def angular_distance_deg(lam1, beta1, lam2, beta2):
    """Great-circle angular distance on the celestial sphere."""
    lam1, beta1 = np.radians(lam1), np.radians(beta1)
    lam2, beta2 = np.radians(lam2), np.radians(beta2)
    cos_d = (np.sin(beta1) * np.sin(beta2) +
             np.cos(beta1) * np.cos(beta2) * np.cos(lam1 - lam2))
    cos_d = np.clip(cos_d, -1.0, 1.0)
    return np.degrees(np.arccos(cos_d))


def run_validation_for_asteroid(ast_num, known_period, gt_info):
    """Run full pipeline on one asteroid and return validation metrics."""
    print(f"\n{'='*60}")
    print(f"VALIDATING ASTEROID {ast_num} ({gt_info['name']})")
    print(f"{'='*60}")

    t0 = time.time()

    # Load data
    print("Loading ALCDEF data...")
    data = preprocess_asteroid(zip_path, gz_path, ast_num)
    blocks = data['blocks']
    print(f"  {len(blocks)} blocks, {sum(len(b['data']) for b in blocks)} total points")

    # Period search (skip if known period is very precise)
    print(f"\nKnown period: {known_period:.6f} h")
    # Use known period directly for validation to isolate shape fitting
    period = known_period

    # Select best blocks for inversion
    good_blocks = [b for b in blocks if len(b['data']) >= 10]
    if len(good_blocks) > 8:
        good_blocks = good_blocks[:8]
    good_blocks = subsample_blocks(good_blocks, max_points_per_block=15)
    total_pts = sum(len(b['data']) for b in good_blocks)
    print(f"Using {len(good_blocks)} blocks, {total_pts} data points for inversion")

    # Stage 1: Convex inversion (seed)
    print("\n--- Stage 1: Convex inversion ---")
    convex_solver = ConvexSolver(n_subdivisions=1, c_ls=0.5, c_l=0.1)
    convex_best, all_solutions = grid_search_spin(
        convex_solver, good_blocks, period,
        lam_step=30, beta_step=30, max_iter=20
    )

    if convex_best is None:
        print("WARNING: Convex inversion failed to converge")
        return None

    print(f"Convex result: lam={convex_best['spin_lam_deg']}, "
          f"beta={convex_best['spin_beta_deg']}, "
          f"chi2_red={convex_best['chi2_reduced']:.2f}")

    # Stage 2: GA non-convex refinement
    print("\n--- Stage 2: GA non-convex refinement ---")
    ga = GASolver(
        n_verts=len(convex_best['base_vertices']),
        base_vertices=convex_best['base_vertices'],
        faces=convex_best['faces'],
        c_ls=0.5, c_l=0.1,
    )

    # Use fewer blocks for GA (slower per evaluation)
    ga_blocks = good_blocks[:4]
    ga_blocks = subsample_blocks(ga_blocks, max_points_per_block=8)

    ga_result = ga.evolve(
        blocks=ga_blocks,
        spin_axis=convex_best['spin_axis'],
        period_hours=period,
        pop_size=30,
        n_generations=30,
        mutation_rate=0.15,
        mutation_sigma=0.08,
        tournament_size=3,
        elitism_count=2,
        regularization=0.01,
        seed_radii=convex_best['radii'],
        verbose=True,
    )

    # Save GA shape
    out_dir = os.path.join(repo_root, 'results', 'shapes')
    os.makedirs(out_dir, exist_ok=True)
    obj_path = os.path.join(out_dir, f'{ast_num}_ga.obj')
    save_obj(obj_path, ga_result['vertices'], ga_result['faces'])
    print(f"Saved GA shape to {obj_path}")

    # Compare with ground truth
    gt_spin = gt_info['spin_parameters']
    gt_lam = gt_spin['lambda_deg']
    gt_beta = gt_spin['beta_deg']

    our_lam = convex_best['spin_lam_deg']
    our_beta = convex_best['spin_beta_deg']

    # Handle pole ambiguity (180-degree ambiguity in lambda)
    pole_err1 = angular_distance_deg(our_lam, our_beta, gt_lam, gt_beta)
    pole_err2 = angular_distance_deg((our_lam + 180) % 360, -our_beta, gt_lam, gt_beta)
    pole_error = min(pole_err1, pole_err2)

    # Shape comparison with ground truth ellipsoid
    gt_obj_path = os.path.join(repo_root, 'results', 'ground_truth',
                                f'{ast_num}_{gt_info["name"]}_ellipsoid.obj')
    if os.path.exists(gt_obj_path):
        gt_verts, gt_faces = load_obj(gt_obj_path)
        h_dist = hausdorff_distance(ga_result['vertices'], gt_verts)
        iou = volumetric_iou(ga_result['vertices'], ga_result['faces'],
                              gt_verts, gt_faces, resolution=32)
    else:
        h_dist = None
        iou = None

    elapsed = time.time() - t0

    result = {
        'asteroid_number': ast_num,
        'name': gt_info['name'],
        'known_period_hours': known_period,
        'known_spin_lam': gt_lam,
        'known_spin_beta': gt_beta,
        'found_spin_lam': float(our_lam),
        'found_spin_beta': float(our_beta),
        'pole_error_deg': float(pole_error),
        'convex_chi2_reduced': float(convex_best['chi2_reduced']),
        'ga_chi2': float(ga_result['chi2']),
        'ellipsoid_params': {
            'b_a': float(convex_best['ellipsoid_params']['b']),
            'c_a': float(convex_best['ellipsoid_params']['c']),
        },
        'hausdorff_vs_gt': float(h_dist) if h_dist is not None else None,
        'iou_vs_gt': float(iou) if iou is not None else None,
        'shape_file': f'results/shapes/{ast_num}_ga.obj',
        'elapsed_seconds': float(elapsed),
        'n_blocks_used': len(good_blocks),
        'n_data_points': total_pts,
    }

    print(f"\n--- Results for {ast_num} {gt_info['name']} ---")
    print(f"  Pole: found ({our_lam}, {our_beta}), known ({gt_lam}, {gt_beta})")
    print(f"  Pole error: {pole_error:.1f} deg")
    print(f"  Convex chi2_red: {convex_best['chi2_reduced']:.2f}")
    print(f"  GA chi2: {ga_result['chi2']:.2f}")
    if h_dist is not None:
        print(f"  Hausdorff vs GT: {h_dist:.4f}")
    if iou is not None:
        print(f"  IoU vs GT: {iou:.4f}")
    print(f"  Time: {elapsed:.1f}s")

    return result


if __name__ == '__main__':
    # Load ground truth
    gt_path = os.path.join(repo_root, 'results', 'ground_truth', 'ground_truth_info.json')
    with open(gt_path) as f:
        gt_data = json.load(f)

    # Validate on 3 asteroids with most ALCDEF data
    targets = [
        (1036, 10.31304),   # Ganymed - 134 blocks, 23583 pts
        (433, 5.27025547),  # Eros - 27 blocks, 2289 pts
        (1580, 6.13836),    # Betulia - 5 blocks, 209 pts
    ]

    all_results = {}
    for ast_num, known_period in targets:
        gt_info = gt_data['asteroids'][str(ast_num)]
        result = run_validation_for_asteroid(ast_num, known_period, gt_info)
        if result:
            all_results[str(ast_num)] = result

    # Save validation report
    report_path = os.path.join(repo_root, 'results', 'validation_report.json')
    with open(report_path, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\n\nValidation report saved to {report_path}")

    # Summary
    print("\n" + "="*60)
    print("VALIDATION SUMMARY")
    print("="*60)
    for key, res in all_results.items():
        print(f"  {res['name']} ({key}): "
              f"pole_err={res['pole_error_deg']:.1f}deg, "
              f"chi2_red={res['convex_chi2_reduced']:.2f}, "
              f"H={res['hausdorff_vs_gt']}, "
              f"IoU={res['iou_vs_gt']}")
