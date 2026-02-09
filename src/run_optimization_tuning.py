"""
Recursive optimization loop: tune pipeline parameters to improve validation.

Item 017 of the research rubric.

Systematically adjusts:
  - Spin-axis grid resolution (lam_step, beta_step)
  - GA mutation rate and sigma
  - Regularization strength
  - Number of data points used

Each tuning iteration is documented with parameter values and resulting metrics.
"""
import numpy as np
import json
import os
import sys
import time

sys.path.insert(0, os.path.dirname(__file__))

from data_ingest import preprocess_asteroid
from convex_inversion import ConvexSolver, grid_search_spin, subsample_blocks
from ga_solver import GASolver
from mesh_utils import hausdorff_distance, volumetric_iou, load_obj

repo_root = os.path.dirname(os.path.dirname(__file__))
zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')


def angular_distance_deg(lam1, beta1, lam2, beta2):
    lam1, beta1 = np.radians(lam1), np.radians(beta1)
    lam2, beta2 = np.radians(lam2), np.radians(beta2)
    cos_d = (np.sin(beta1) * np.sin(beta2) +
             np.cos(beta1) * np.cos(beta2) * np.cos(lam1 - lam2))
    cos_d = np.clip(cos_d, -1.0, 1.0)
    return np.degrees(np.arccos(cos_d))


def run_single_iteration(ast_num, known_period, gt_info, params):
    """Run pipeline with given parameters and return metrics."""
    data = preprocess_asteroid(zip_path, gz_path, ast_num)
    blocks = data['blocks']

    good_blocks = [b for b in blocks if len(b['data']) >= params.get('min_block_size', 10)]
    max_blocks = params.get('max_blocks', 8)
    if len(good_blocks) > max_blocks:
        good_blocks = good_blocks[:max_blocks]
    good_blocks = subsample_blocks(good_blocks,
                                    max_points_per_block=params.get('max_pts_per_block', 15))
    total_pts = sum(len(b['data']) for b in good_blocks)

    # Convex inversion
    convex_solver = ConvexSolver(n_subdivisions=1, c_ls=params['c_ls'], c_l=params['c_l'])
    convex_best, _ = grid_search_spin(
        convex_solver, good_blocks, known_period,
        lam_step=params['lam_step'], beta_step=params['beta_step'],
        max_iter=params.get('max_iter', 20)
    )

    if convex_best is None:
        return None

    # GA refinement
    ga = GASolver(
        n_verts=len(convex_best['base_vertices']),
        base_vertices=convex_best['base_vertices'],
        faces=convex_best['faces'],
        c_ls=params['c_ls'], c_l=params['c_l'],
    )

    ga_blocks = good_blocks[:params.get('ga_n_blocks', 4)]
    ga_blocks = subsample_blocks(ga_blocks,
                                  max_points_per_block=params.get('ga_max_pts', 8))

    ga_result = ga.evolve(
        blocks=ga_blocks,
        spin_axis=convex_best['spin_axis'],
        period_hours=known_period,
        pop_size=params.get('ga_pop_size', 30),
        n_generations=params.get('ga_n_gen', 30),
        mutation_rate=params.get('ga_mutation_rate', 0.15),
        mutation_sigma=params.get('ga_mutation_sigma', 0.08),
        regularization=params.get('ga_regularization', 0.01),
        seed_radii=convex_best['radii'],
        verbose=False,
    )

    # Evaluate
    gt_spin = gt_info['spin_parameters']
    our_lam = float(convex_best['spin_lam_deg'])
    our_beta = float(convex_best['spin_beta_deg'])
    pole_err1 = angular_distance_deg(our_lam, our_beta,
                                      gt_spin['lambda_deg'], gt_spin['beta_deg'])
    pole_err2 = angular_distance_deg((our_lam + 180) % 360, -our_beta,
                                      gt_spin['lambda_deg'], gt_spin['beta_deg'])
    pole_error = min(pole_err1, pole_err2)

    gt_obj_path = os.path.join(repo_root, 'results', 'ground_truth',
                                f'{ast_num}_{gt_info["name"]}_ellipsoid.obj')
    h_dist, iou = None, None
    if os.path.exists(gt_obj_path):
        gt_verts, gt_faces = load_obj(gt_obj_path)
        h_dist = float(hausdorff_distance(ga_result['vertices'], gt_verts))
        iou = float(volumetric_iou(ga_result['vertices'], ga_result['faces'],
                                     gt_verts, gt_faces, resolution=32))

    return {
        'spin_lam': our_lam,
        'spin_beta': our_beta,
        'pole_error_deg': float(pole_error),
        'convex_chi2_red': float(convex_best['chi2_reduced']),
        'ga_chi2': float(ga_result['chi2']),
        'hausdorff': h_dist,
        'iou': iou,
        'n_data_pts': total_pts,
    }


if __name__ == '__main__':
    # Load ground truth
    gt_path = os.path.join(repo_root, 'results', 'ground_truth', 'ground_truth_info.json')
    with open(gt_path) as f:
        gt_data = json.load(f)

    # Primary test target: Ganymed (most data, best initial result)
    test_target = (1036, 10.31304)
    gt_info = gt_data['asteroids']['1036']

    optimization_log = []

    # ============================================================
    # Iteration 1: Baseline parameters (from initial validation)
    # ============================================================
    print("="*60)
    print("ITERATION 1: Baseline parameters")
    print("="*60)
    params_1 = {
        'lam_step': 30, 'beta_step': 30,
        'c_ls': 0.5, 'c_l': 0.1,
        'max_blocks': 8, 'max_pts_per_block': 15,
        'min_block_size': 10,
        'ga_n_blocks': 4, 'ga_max_pts': 8,
        'ga_pop_size': 30, 'ga_n_gen': 30,
        'ga_mutation_rate': 0.15, 'ga_mutation_sigma': 0.08,
        'ga_regularization': 0.01,
        'max_iter': 20,
    }

    t0 = time.time()
    result_1 = run_single_iteration(*test_target, gt_info, params_1)
    elapsed_1 = time.time() - t0

    iter_1 = {
        'iteration': 1,
        'description': 'Baseline parameters',
        'parameters': params_1,
        'results': result_1,
        'elapsed_seconds': elapsed_1,
    }
    optimization_log.append(iter_1)
    print(f"  Pole error: {result_1['pole_error_deg']:.1f} deg")
    print(f"  Hausdorff: {result_1['hausdorff']}")
    print(f"  IoU: {result_1['iou']}")
    print(f"  Time: {elapsed_1:.1f}s")

    # ============================================================
    # Iteration 2: Finer spin grid + more data + higher GA pop
    # ============================================================
    print("\n" + "="*60)
    print("ITERATION 2: Finer spin grid (15 deg), more data, larger GA population")
    print("="*60)
    params_2 = {
        'lam_step': 15, 'beta_step': 15,
        'c_ls': 0.5, 'c_l': 0.1,
        'max_blocks': 12, 'max_pts_per_block': 20,
        'min_block_size': 8,
        'ga_n_blocks': 5, 'ga_max_pts': 10,
        'ga_pop_size': 40, 'ga_n_gen': 40,
        'ga_mutation_rate': 0.2, 'ga_mutation_sigma': 0.1,
        'ga_regularization': 0.005,
        'max_iter': 25,
    }

    t0 = time.time()
    result_2 = run_single_iteration(*test_target, gt_info, params_2)
    elapsed_2 = time.time() - t0

    iter_2 = {
        'iteration': 2,
        'description': 'Finer spin grid (15 deg steps), more data points, larger GA population, lower regularization',
        'parameters': params_2,
        'results': result_2,
        'elapsed_seconds': elapsed_2,
    }
    optimization_log.append(iter_2)
    print(f"  Pole error: {result_2['pole_error_deg']:.1f} deg")
    print(f"  Hausdorff: {result_2['hausdorff']}")
    print(f"  IoU: {result_2['iou']}")
    print(f"  Time: {elapsed_2:.1f}s")

    # ============================================================
    # Iteration 3: Different scattering weights
    # ============================================================
    print("\n" + "="*60)
    print("ITERATION 3: Adjusted scattering law weights (more Lommel-Seeliger)")
    print("="*60)
    params_3 = {
        'lam_step': 15, 'beta_step': 15,
        'c_ls': 0.7, 'c_l': 0.05,
        'max_blocks': 12, 'max_pts_per_block': 20,
        'min_block_size': 8,
        'ga_n_blocks': 5, 'ga_max_pts': 10,
        'ga_pop_size': 40, 'ga_n_gen': 40,
        'ga_mutation_rate': 0.2, 'ga_mutation_sigma': 0.1,
        'ga_regularization': 0.005,
        'max_iter': 25,
    }

    t0 = time.time()
    result_3 = run_single_iteration(*test_target, gt_info, params_3)
    elapsed_3 = time.time() - t0

    iter_3 = {
        'iteration': 3,
        'description': 'Adjusted scattering weights: c_ls=0.7, c_l=0.05 (more LS-dominant)',
        'parameters': params_3,
        'results': result_3,
        'elapsed_seconds': elapsed_3,
    }
    optimization_log.append(iter_3)
    print(f"  Pole error: {result_3['pole_error_deg']:.1f} deg")
    print(f"  Hausdorff: {result_3['hausdorff']}")
    print(f"  IoU: {result_3['iou']}")
    print(f"  Time: {elapsed_3:.1f}s")

    # Save optimization log
    log_path = os.path.join(repo_root, 'results', 'optimization_log.json')
    with open(log_path, 'w') as f:
        json.dump(optimization_log, f, indent=2)
    print(f"\nOptimization log saved to {log_path}")

    # Summary table
    print("\n" + "="*60)
    print("OPTIMIZATION SUMMARY")
    print("="*60)
    print(f"{'Iter':<6} {'Pole Err':<10} {'Chi2_red':<10} {'Hausdorff':<10} {'IoU':<8} {'Time(s)':<8}")
    print("-" * 52)
    for entry in optimization_log:
        r = entry['results']
        print(f"  {entry['iteration']:<4} {r['pole_error_deg']:<10.1f} "
              f"{r['convex_chi2_red']:<10.2f} "
              f"{r['hausdorff']:<10.4f} {r['iou']:<8.4f} "
              f"{entry['elapsed_seconds']:<8.1f}")
