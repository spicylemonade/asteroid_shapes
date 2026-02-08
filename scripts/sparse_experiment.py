#!/usr/bin/env python3
"""
Item 021: Sparse-data inversion experiment.

Select 3 asteroids with dense ALCDEF coverage. Subsample their lightcurves
to simulate sparse survey conditions. Run sparse inversion and compare
recovered parameters against dense-data results.
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
from lci_engine.sparse_inversion import subsample_to_sparse, sparse_inversion

ALCDEF_ZIP = os.path.join(os.path.dirname(__file__), '..', 'ALCDEF_ALL.zip')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')

np.random.seed(42)

def log(msg):
    print(msg, flush=True)


# Use 3 asteroids with good dense coverage from our ground truth set
TARGETS = [
    {"key": "433_Eros", "number": 433, "known_period": 5.27025},
    {"key": "216_Kleopatra", "number": 216, "known_period": 5.385},
    {"key": "1943_Anteros", "number": 1943, "known_period": None},  # NEO candidate
]


def run_dense_vs_sparse(target):
    key = target['key']
    number = target['number']
    known_period = target['known_period']

    log(f"\n{'='*60}")
    log(f"SPARSE EXPERIMENT: {key} (#{number})")
    log(f"{'='*60}")

    # Load data
    sessions = load_alcdef_asteroid(ALCDEF_ZIP, asteroid_number=number)
    log(f"Loaded {len(sessions)} sessions")
    if len(sessions) < 3:
        log("Too few sessions, skipping")
        return None

    sessions.sort(key=lambda s: len(s.jd), reverse=True)
    sessions = sessions[:15]
    total_pts = sum(len(s.jd) for s in sessions)
    log(f"Using {len(sessions)} sessions, {total_pts} points")

    # Stage 1: Dense period search
    log("\n--- Dense Period Search ---")
    t0 = time.time()
    times, mags, errs = combine_lightcurve_sessions(sessions)
    dense_period, _, _ = find_best_period(
        times, mags, errs, period_min=2.0, period_max=50.0, timeout_sec=120
    )
    # Handle P/2 alias if known period available
    if known_period:
        alias_errs = [
            abs(dense_period - known_period) / known_period,
            abs(dense_period * 2 - known_period) / known_period,
            abs(dense_period / 2 - known_period) / known_period,
        ]
        best_alias = int(np.argmin(alias_errs))
        if best_alias == 1:
            dense_period = dense_period * 2
        elif best_alias == 2:
            dense_period = dense_period / 2
    log(f"Dense period: {dense_period:.4f} h [{time.time()-t0:.1f}s]")

    # Stage 2: Dense inversion (1 pole trial for speed)
    log("\n--- Dense Inversion ---")
    t0 = time.time()
    try:
        dense_result = convex_inversion(
            sessions,
            period_init=dense_period,
            pole_lambda_init=0.0,
            pole_beta_init=np.pi/4,
            n_facets=100,
            lambda_smooth=0.3,
            max_iter=60,
            seed=42
        )
        log(f"Dense: P={dense_result.period:.4f}h pole=({np.degrees(dense_result.pole_lambda):.1f},{np.degrees(dense_result.pole_beta):.1f}) rms={dense_result.residual_rms:.4f} [{time.time()-t0:.1f}s]")
    except Exception as e:
        log(f"Dense inversion failed: {e}")
        return None

    # Stage 3: Create sparse subsample
    log("\n--- Sparse Subsample ---")
    sparse_data = subsample_to_sparse(
        sessions, n_points_per_apparition=50, n_apparitions=3, seed=42
    )
    log(f"Sparse points: {len(sparse_data)} (from {total_pts} dense)")

    if len(sparse_data) < 20:
        log("Too few sparse points, skipping")
        return None

    # Stage 4: Sparse inversion
    log("\n--- Sparse Inversion ---")
    t0 = time.time()
    try:
        sparse_result = sparse_inversion(
            sparse_data,
            period_init=dense_period,
            pole_lambda_init=0.0,
            pole_beta_init=np.pi/4,
            H_init=10.0,
            n_facets=80,
            lambda_smooth=0.5,
            max_iter=60,
            seed=42
        )
        sparse_period = sparse_result['period']
        sparse_lambda = sparse_result['pole_lambda']
        sparse_beta = sparse_result['pole_beta']
        log(f"Sparse: P={sparse_period:.4f}h pole=({np.degrees(sparse_lambda):.1f},{np.degrees(sparse_beta):.1f}) chi2={sparse_result['chi2']:.1f} [{time.time()-t0:.1f}s]")
    except Exception as e:
        log(f"Sparse inversion failed: {e}")
        return None

    # Stage 5: Compare
    log("\n--- Comparison ---")
    # Period comparison
    if known_period:
        dense_period_err = min(
            abs(dense_result.period - known_period) / known_period * 100,
            abs(dense_result.period * 2 - known_period) / known_period * 100,
            abs(dense_result.period / 2 - known_period) / known_period * 100,
        )
        sparse_period_err = min(
            abs(sparse_period - known_period) / known_period * 100,
            abs(sparse_period * 2 - known_period) / known_period * 100,
            abs(sparse_period / 2 - known_period) / known_period * 100,
        )
    else:
        # Compare sparse vs dense
        dense_period_err = 0.0
        sparse_period_err = abs(sparse_period - dense_result.period) / dense_result.period * 100

    # Pole comparison (sparse vs dense)
    from lci_engine.validation import pole_direction_error
    pole_diff = pole_direction_error(
        sparse_lambda, sparse_beta,
        dense_result.pole_lambda, dense_result.pole_beta
    )

    log(f"  Dense period error: {dense_period_err:.3f}%")
    log(f"  Sparse period error: {sparse_period_err:.3f}%")
    log(f"  Pole difference (sparse vs dense): {pole_diff:.1f} deg")
    log(f"  Dense RMS: {dense_result.residual_rms:.4f} mag")

    return {
        'asteroid': key,
        'number': number,
        'n_dense_points': total_pts,
        'n_sparse_points': len(sparse_data),
        'n_apparitions_sparse': 3,
        'dense': {
            'period_hours': float(dense_result.period),
            'pole_lambda_deg': float(np.degrees(dense_result.pole_lambda)),
            'pole_beta_deg': float(np.degrees(dense_result.pole_beta)),
            'residual_rms': float(dense_result.residual_rms),
            'period_error_pct': float(dense_period_err),
        },
        'sparse': {
            'period_hours': float(sparse_period),
            'pole_lambda_deg': float(np.degrees(sparse_lambda)),
            'pole_beta_deg': float(np.degrees(sparse_beta)),
            'chi2': float(sparse_result['chi2']),
            'period_error_pct': float(sparse_period_err),
        },
        'degradation': {
            'period_error_increase_pct': float(sparse_period_err - dense_period_err),
            'pole_difference_deg': float(pole_diff),
        },
        'known_period': known_period,
    }


def main():
    all_results = {}
    for target in TARGETS:
        try:
            result = run_dense_vs_sparse(target)
            if result:
                all_results[target['key']] = result
        except Exception as e:
            log(f"\nERROR on {target['key']}: {e}")
            import traceback
            traceback.print_exc()

    out_path = os.path.join(RESULTS_DIR, 'sparse_experiment_results.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    log(f"\n{'='*60}")
    log("SPARSE EXPERIMENT SUMMARY")
    log(f"{'='*60}")
    for key, res in all_results.items():
        d = res['dense']
        s = res['sparse']
        deg = res['degradation']
        log(f"  {key}:")
        log(f"    Dense:  P={d['period_hours']:.4f}h err={d['period_error_pct']:.3f}%")
        log(f"    Sparse: P={s['period_hours']:.4f}h err={s['period_error_pct']:.3f}%")
        log(f"    Pole diff: {deg['pole_difference_deg']:.1f} deg")
        log(f"    Points: {res['n_dense_points']} dense -> {res['n_sparse_points']} sparse")

    log(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
