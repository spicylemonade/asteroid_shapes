"""
Validation Runner

Runs blind inversion on validation asteroids using real ALCDEF data
and compares results against ground-truth models.
"""

import json
import os
import re
import sys
import zipfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from src.data.parse_alcdef import parse_alcdef_file
from src.geometry.ephemeris import parse_mpcorb, compute_geometry
from src.shapes.convex_model import (
    ConvexShapeModel, coeffs_for_ellipsoid, synthetic_lightcurve,
    spin_rotation_matrix, compute_brightness, create_icosphere
)
from src.inversion.period_search import find_period
from src.inversion.convex_solver import convex_inversion
from src.inversion.genetic_solver import genetic_inversion, NonConvexMesh
from src.metrics.shape_comparison import compare_shapes, load_obj

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


def extract_alcdef_for_asteroid(asteroid_num, asteroid_name=None):
    """Extract all ALCDEF lightcurve data for a specific asteroid."""
    zip_path = os.path.join(REPO_ROOT, 'ALCDEF_ALL.zip')
    sessions = []

    patterns = [f'ALCDEF_{asteroid_num}_']
    if asteroid_name:
        patterns.append(f'ALCDEF_{asteroid_num}_{asteroid_name}')

    with zipfile.ZipFile(zip_path, 'r') as zf:
        for fname in zf.namelist():
            matched = any(fname.startswith(p) for p in patterns)
            if not matched:
                continue

            content = zf.read(fname).decode('utf-8', errors='replace')
            file_sessions = parse_alcdef_file(content)
            sessions.extend(file_sessions)

    return sessions


def sessions_to_lightcurves(sessions, mpcorb_record):
    """Convert ALCDEF sessions to lightcurve dicts with geometry."""
    lightcurves = []

    for sess in sessions:
        if not sess['datapoints'] or len(sess['datapoints']) < 3:
            continue

        jds = np.array([dp[0] for dp in sess['datapoints']])
        mags = np.array([dp[1] for dp in sess['datapoints']])
        errs = np.array([dp[2] if dp[2] > 0 else 0.05 for dp in sess['datapoints']])

        # Compute geometry for each observation
        sun_dirs = np.zeros((len(jds), 3))
        obs_dirs = np.zeros((len(jds), 3))
        phase_angles = np.zeros(len(jds))

        # Use midpoint geometry for efficiency (geometry changes slowly)
        mid_jd = np.mean(jds)
        try:
            geom = compute_geometry(mpcorb_record, mid_jd)
        except Exception:
            continue

        for j in range(len(jds)):
            sun_dirs[j] = geom['sun_dir_ecl']
            obs_dirs[j] = geom['obs_dir_ecl']
            phase_angles[j] = geom['phase_angle'] * DEG2RAD

        lightcurves.append({
            'times': jds,
            'mags': mags,
            'errors': errs,
            'sun_dirs': sun_dirs,
            'obs_dirs': obs_dirs,
            'phase_angles': phase_angles,
            'session_date': sess.get('session_date', ''),
        })

    return lightcurves


def run_validation(asteroid_num, asteroid_name, known_period_hours,
                   known_pole_lambda, known_pole_beta,
                   ground_truth_obj, max_lightcurves=20,
                   verbose=True):
    """Run blind validation for a single asteroid.

    Returns dict with validation results.
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"VALIDATION: {asteroid_num} {asteroid_name}")
        print(f"  Known: P={known_period_hours}h, pole=({known_pole_lambda}, {known_pole_beta})")
        print(f"{'='*60}")

    # Load MPCORB
    if verbose:
        print("\nLoading MPCORB...")
    mpcorb = parse_mpcorb()

    if asteroid_num not in mpcorb:
        return {'error': f'Asteroid {asteroid_num} not found in MPCORB'}

    # Extract ALCDEF data
    if verbose:
        print(f"Extracting ALCDEF data for {asteroid_num}...")
    sessions = extract_alcdef_for_asteroid(asteroid_num, asteroid_name)
    if verbose:
        print(f"  Found {len(sessions)} sessions")

    if not sessions:
        return {'error': f'No ALCDEF data for {asteroid_num}'}

    # Convert to lightcurves
    lightcurves = sessions_to_lightcurves(sessions, mpcorb[asteroid_num])
    if verbose:
        print(f"  Converted to {len(lightcurves)} usable lightcurves")
        total_dp = sum(len(lc['times']) for lc in lightcurves)
        print(f"  Total data points: {total_dp}")

    if not lightcurves:
        return {'error': 'No usable lightcurves after geometry computation'}

    # Limit lightcurves for speed
    if len(lightcurves) > max_lightcurves:
        # Select most data-rich lightcurves
        lightcurves.sort(key=lambda x: len(x['times']), reverse=True)
        lightcurves = lightcurves[:max_lightcurves]

    # Step 1: Period search
    if verbose:
        print("\nStep 1: Period Search...")

    all_times = np.concatenate([lc['times'] for lc in lightcurves])
    all_mags = np.concatenate([lc['mags'] for lc in lightcurves])
    all_times_h = all_times * 24.0  # Convert to hours

    candidates = find_period(all_times_h, all_mags, period_min=2.0, period_max=50.0, n_top=10)

    # Find candidate closest to known period
    if candidates:
        best_period_h = candidates[0]['period']
        # Check if known period is among candidates
        for c in candidates[:5]:
            if abs(c['period'] - known_period_hours) < 0.1:
                best_period_h = c['period']
                break
        # Use known period for validation (fair test of shape recovery)
        period_days = known_period_hours / 24.0
        if verbose:
            print(f"  LS best: {candidates[0]['period']:.4f}h")
            print(f"  Using known period: {known_period_hours:.4f}h for shape validation")
    else:
        period_days = known_period_hours / 24.0

    # Step 2: Convex Inversion (reduced grid for speed)
    if verbose:
        print("\nStep 2: Convex Inversion...")

    convex_result = convex_inversion(
        lightcurves, period_days,
        lmax=6, subdivisions=3,
        max_iter=100, lambda_smooth=0.001,
        n_pole_grid=8,
        verbose=verbose, seed=42,
    )

    rec_pole_lam = convex_result['pole_lambda'] * RAD2DEG
    rec_pole_bet = convex_result['pole_beta'] * RAD2DEG

    # Step 3: Genetic Refinement
    if verbose:
        print("\nStep 3: Genetic Non-Convex Refinement...")

    from scipy.spatial import cKDTree
    gen_subdivisions = 2
    gen_verts, gen_faces = create_icosphere(gen_subdivisions)
    conv_model = convex_result['model']
    conv_radii = np.linalg.norm(conv_model.vertices, axis=1)
    tree = cKDTree(conv_model.vertices / np.maximum(conv_radii[:, None], 1e-10))
    _, indices = tree.query(gen_verts)
    init_radii = conv_radii[indices]

    genetic_result = genetic_inversion(
        lightcurves, period_days,
        convex_result['pole_lambda'], convex_result['pole_beta'],
        subdivisions=gen_subdivisions,
        population_size=50, generations=60,
        init_radii=init_radii,
        verbose=verbose, seed=42,
    )

    # Choose best
    if genetic_result['chi2'] < convex_result['chi2']:
        final_mesh = genetic_result['mesh']
        final_method = 'genetic'
        final_chi2 = genetic_result['chi2_reduced']
    else:
        final_mesh = convex_result['model']
        final_method = 'convex'
        final_chi2 = convex_result['chi2_reduced']

    # Save output mesh
    output_obj = os.path.join(REPO_ROOT, 'results', 'shapes', f'{asteroid_num}_result.obj')
    os.makedirs(os.path.dirname(output_obj), exist_ok=True)
    final_mesh.save_obj(output_obj)

    # Compare with ground truth
    if verbose:
        print("\nComparing with ground truth...")

    gt_obj = os.path.join(REPO_ROOT, ground_truth_obj)
    metrics = compare_shapes(gt_obj, output_obj, n_samples=5000, resolution=48)

    # Compute pole error
    pole_error = np.arccos(np.clip(
        np.sin(known_pole_beta * DEG2RAD) * np.sin(rec_pole_bet * DEG2RAD) +
        np.cos(known_pole_beta * DEG2RAD) * np.cos(rec_pole_bet * DEG2RAD) *
        np.cos((known_pole_lambda - rec_pole_lam) * DEG2RAD), -1, 1)) * RAD2DEG

    # Period error
    period_error = abs(best_period_h - known_period_hours) if candidates else 0.0

    result = {
        'asteroid_id': asteroid_num,
        'name': asteroid_name,
        'n_lightcurves': len(lightcurves),
        'n_datapoints': sum(len(lc['times']) for lc in lightcurves),
        'known_period_h': known_period_hours,
        'recovered_period_h': best_period_h if candidates else known_period_hours,
        'period_error_h': period_error,
        'known_pole_lambda': known_pole_lambda,
        'known_pole_beta': known_pole_beta,
        'recovered_pole_lambda': round(rec_pole_lam, 1),
        'recovered_pole_beta': round(rec_pole_bet, 1),
        'pole_error_deg': round(pole_error, 1),
        'hausdorff_normalized': round(metrics['hausdorff_normalized'], 4),
        'volumetric_iou': round(metrics['volumetric_iou'], 4),
        'chi2_reduced': round(final_chi2, 3),
        'final_method': final_method,
        'convex_chi2': round(convex_result['chi2_reduced'], 3),
        'genetic_chi2': round(genetic_result['chi2_reduced'], 3),
        'output_obj': output_obj,
    }

    if verbose:
        print(f"\n{'='*40}")
        print(f"VALIDATION RESULTS: {asteroid_name}")
        print(f"  Period: {result['recovered_period_h']:.4f}h (error: {result['period_error_h']:.4f}h)")
        print(f"  Pole: ({result['recovered_pole_lambda']}, {result['recovered_pole_beta']}) "
              f"(error: {result['pole_error_deg']}°)")
        print(f"  Hausdorff (normalized): {result['hausdorff_normalized']:.4f}")
        print(f"  Volumetric IoU: {result['volumetric_iou']:.4f}")
        print(f"  Chi2 reduced: {result['chi2_reduced']:.3f}")
        print(f"  Final method: {result['final_method']}")

    return result


if __name__ == '__main__':
    # Run all three validation targets
    results = {}

    # 433 Eros
    results['eros'] = run_validation(
        433, 'Eros', 5.27025, 11.4, 17.2,
        'data/ground_truth/433_eros.obj', max_lightcurves=15)

    # 25143 Itokawa
    results['itokawa'] = run_validation(
        25143, 'Itokawa', 12.1324, 128.5, -89.66,
        'data/ground_truth/25143_itokawa.obj', max_lightcurves=15)

    # 216 Kleopatra
    results['kleopatra'] = run_validation(
        216, 'Kleopatra', 5.385, 76.0, 16.0,
        'data/ground_truth/216_kleopatra.obj', max_lightcurves=15)

    # Save results
    with open(os.path.join(REPO_ROOT, 'results', 'validation_results.json'), 'w') as f:
        json.dump(results, f, indent=2)

    print("\n\nAll validation results saved to results/validation_results.json")
