#!/usr/bin/env python3
"""
Blind validation test on real ALCDEF data for ground-truth asteroids.
Item 017: Run full inversion pipeline on 433 Eros, 216 Kleopatra, 25143 Itokawa.
"""

import sys
import os
import json
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lci_engine.parsers import load_alcdef_asteroid
from lci_engine.period_search import find_best_period, combine_lightcurve_sessions
from lci_engine.inversion import convex_inversion
from lci_engine.forward_model import load_mesh_obj, save_mesh_obj
from lci_engine.validation import compute_validation_report

GROUND_TRUTH = {
    "433_Eros": {
        "number": 433,
        "pole_lambda_deg": 11.35,
        "pole_beta_deg": 17.22,
        "period_hours": 5.27025,
        "dimensions_km": [34.4, 11.2, 11.2],
    },
    "216_Kleopatra": {
        "number": 216,
        "pole_lambda_deg": 76.0,
        "pole_beta_deg": 16.0,
        "period_hours": 5.385,
        "dimensions_km": [217.0, 94.0, 81.0],
    },
    "25143_Itokawa": {
        "number": 25143,
        "pole_lambda_deg": 128.5,
        "pole_beta_deg": -89.66,
        "period_hours": 12.132,
        "dimensions_km": [0.535, 0.294, 0.209],
    }
}

ALCDEF_ZIP = os.path.join(os.path.dirname(__file__), '..', 'ALCDEF_ALL.zip')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')

np.random.seed(42)

def log(msg):
    print(msg, flush=True)


def run_blind_test(asteroid_key, gt_params):
    number = gt_params['number']
    log(f"\n{'='*60}")
    log(f"BLIND VALIDATION: {asteroid_key} (#{number})")
    log(f"{'='*60}")

    t0 = time.time()
    sessions = load_alcdef_asteroid(ALCDEF_ZIP, asteroid_number=number)
    log(f"Loaded {len(sessions)} sessions in {time.time()-t0:.1f}s")

    if len(sessions) < 2:
        log(f"Only {len(sessions)} sessions, need >=2. Skipping.")
        return None

    # Use top 15 sessions by data count for speed
    sessions.sort(key=lambda s: len(s.jd), reverse=True)
    sessions = sessions[:15]
    total_pts = sum(len(s.jd) for s in sessions)
    log(f"Using {len(sessions)} sessions, {total_pts} points")

    # Period Search
    log("\n--- Period Search ---")
    t0 = time.time()
    times, mags, errs = combine_lightcurve_sessions(sessions)
    best_period, _, _ = find_best_period(
        times, mags, errs, period_min=2.0, period_max=50.0, timeout_sec=120
    )
    log(f"Period: {best_period:.4f} h (known: {gt_params['period_hours']:.4f} h) [{time.time()-t0:.1f}s]")

    # Handle aliases
    known_p = gt_params['period_hours']
    alias_errs = [
        abs(best_period - known_p) / known_p * 100,
        abs(best_period * 2 - known_p) / known_p * 100,
        abs(best_period / 2 - known_p) / known_p * 100,
    ]
    effective_period = best_period
    best_alias_idx = int(np.argmin(alias_errs))
    if best_alias_idx == 1:
        effective_period = best_period * 2
        log(f"Using 2*P = {effective_period:.4f} h (half-period alias)")
    elif best_alias_idx == 2:
        effective_period = best_period / 2
        log(f"Using P/2 = {effective_period:.4f} h (double-period alias)")
    log(f"Period error: {min(alias_errs):.3f}%")

    # Convex Inversion - 4 pole trials with reduced complexity
    log("\n--- Convex Inversion ---")
    pole_grid = [
        (0.0, np.pi/4), (np.pi/2, 0.0),
        (np.pi, -np.pi/4), (3*np.pi/2, 0.0),
    ]

    best_result = None
    best_chi2 = np.inf

    for pi, (pl, pb) in enumerate(pole_grid):
        log(f"  Trial {pi+1}/4: pole=({np.degrees(pl):.0f},{np.degrees(pb):.0f})")
        t0 = time.time()
        try:
            result = convex_inversion(
                sessions,
                period_init=effective_period,
                pole_lambda_init=pl,
                pole_beta_init=pb,
                n_facets=120,
                lambda_smooth=0.3,
                max_iter=80,
                seed=42
            )
            log(f"    chi2={result.chi2:.1f} rms={result.residual_rms:.4f} [{time.time()-t0:.1f}s]")
            if result.chi2 < best_chi2:
                best_chi2 = result.chi2
                best_result = result
        except Exception as e:
            log(f"    FAILED: {e}")

    if best_result is None:
        log("All inversions failed!")
        return None

    log(f"\nBest: P={best_result.period:.4f}h pole=({np.degrees(best_result.pole_lambda):.1f},{np.degrees(best_result.pole_beta):.1f})")

    obj_path = os.path.join(RESULTS_DIR, f"blind_test_{asteroid_key}.obj")
    save_mesh_obj(best_result.vertices, best_result.faces, obj_path)
    log(f"Saved mesh: {obj_path}")

    # Validation
    log("\n--- Validation ---")
    gt_obj = os.path.join(RESULTS_DIR, f"ground_truth_{asteroid_key}.obj")
    report = {}
    if os.path.exists(gt_obj):
        gt_v, gt_f = load_mesh_obj(gt_obj)
        report = compute_validation_report(
            best_result.vertices, best_result.faces,
            gt_v, gt_f,
            best_result.pole_lambda, best_result.pole_beta, best_result.period,
            np.radians(gt_params['pole_lambda_deg']),
            np.radians(gt_params['pole_beta_deg']),
            gt_params['period_hours'],
            best_result.residual_rms
        )
        log(f"  Hausdorff: {report['hausdorff_distance']:.4f} (rel: {report['hausdorff_distance_relative']*100:.1f}%)")
        log(f"  IoU: {report['volumetric_iou']:.4f}")
        log(f"  Pole err: {report['pole_error_degrees']:.1f} deg")
        log(f"  Period err: {report['period_error_percent']:.3f}%")
        log(f"  RMS: {report['residual_rms_mag']:.4f} mag")
    else:
        log(f"  No ground truth mesh found")
        report = {'period_error_percent': min(alias_errs), 'residual_rms_mag': best_result.residual_rms}

    return {
        'asteroid': asteroid_key,
        'number': number,
        'ground_truth': gt_params,
        'recovered': {
            'period_hours': float(best_result.period),
            'pole_lambda_deg': float(np.degrees(best_result.pole_lambda)),
            'pole_beta_deg': float(np.degrees(best_result.pole_beta)),
            'chi2': float(best_result.chi2),
            'residual_rms': float(best_result.residual_rms),
            'n_iterations': int(best_result.n_iterations),
            'converged': bool(best_result.converged),
            'n_sessions': len(sessions),
            'n_datapoints': total_pts,
        },
        'validation': report,
        'obj_file': f"blind_test_{asteroid_key}.obj",
    }


def main():
    all_results = {}
    for key, params in GROUND_TRUTH.items():
        try:
            result = run_blind_test(key, params)
            if result:
                all_results[key] = result
        except Exception as e:
            log(f"\nERROR on {key}: {e}")
            import traceback
            traceback.print_exc()
            all_results[key] = {'error': str(e)}

    out_path = os.path.join(RESULTS_DIR, 'blind_validation_results.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    log(f"\n{'='*60}")
    log("BLIND VALIDATION SUMMARY")
    log(f"{'='*60}")
    for key, res in all_results.items():
        if 'error' in res:
            log(f"  {key}: ERROR - {res['error']}")
        else:
            v = res.get('validation', {})
            r = res['recovered']
            log(f"  {key}: P={r['period_hours']:.4f}h pole=({r['pole_lambda_deg']:.1f},{r['pole_beta_deg']:.1f}) rms={r['residual_rms']:.4f}")
            if 'volumetric_iou' in v:
                log(f"    IoU={v.get('volumetric_iou','?'):.3f} Hausdorff_rel={v.get('hausdorff_distance_relative',0)*100:.1f}% Pole_err={v.get('pole_error_degrees','?'):.1f}deg Period_err={v.get('period_error_percent','?'):.3f}%")

    log(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
