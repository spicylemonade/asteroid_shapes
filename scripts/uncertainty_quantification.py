#!/usr/bin/env python3
"""
Item 022: Uncertainty quantification via bootstrap resampling.

For each modeled asteroid, run bootstrap resampling (50+ iterations) to estimate
uncertainties on period, pole direction, and shape parameters. Report 1-sigma
confidence intervals and flag low-confidence models.
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

ALCDEF_ZIP = os.path.join(os.path.dirname(__file__), '..', 'ALCDEF_ALL.zip')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')

np.random.seed(42)

def log(msg):
    """Print a message immediately with flushed output."""
    print(msg, flush=True)


# Targets for uncertainty quantification: ground truth + top candidates
TARGETS = [
    {"key": "433_Eros", "number": 433, "period_init": 5.2594},
    {"key": "216_Kleopatra", "number": 216, "period_init": 5.3881},
    {"key": "25143_Itokawa", "number": 25143, "period_init": 12.0848},
]

N_BOOTSTRAP = 50


def bootstrap_asteroid(target):
    """Estimate parameter uncertainties for one asteroid via bootstrap resampling.

    Runs N_BOOTSTRAP inversions on resampled (with replacement) session sets,
    then computes mean, standard deviation, and confidence flags for period
    and pole direction.

    Args:
        target: Dict with 'key', 'number', and 'period_init'.

    Returns:
        Dict with bootstrap statistics and confidence flag, or None if
        insufficient data.
    """
    key = target['key']
    number = target['number']
    period_init = target['period_init']

    log(f"\n{'='*60}")
    log(f"UNCERTAINTY: {key} (#{number})")
    log(f"{'='*60}")

    # Load data
    sessions = load_alcdef_asteroid(ALCDEF_ZIP, asteroid_number=number)
    if len(sessions) < 3:
        log(f"Only {len(sessions)} sessions, skipping")
        return None

    sessions.sort(key=lambda s: len(s.jd), reverse=True)
    sessions = sessions[:10]
    total_pts = sum(len(s.jd) for s in sessions)
    log(f"Using {len(sessions)} sessions, {total_pts} points")

    # Get initial period if not provided
    if period_init is None:
        log("Running period search...")
        times, mags, errs = combine_lightcurve_sessions(sessions)
        period_init, _, _ = find_best_period(
            times, mags, errs, period_min=2.0, period_max=50.0, timeout_sec=90
        )
        log(f"Period: {period_init:.4f} h")

    # Bootstrap resampling
    rng = np.random.RandomState(42)
    periods = []
    lambdas = []
    betas = []
    rms_values = []

    for i in range(N_BOOTSTRAP):
        # Resample sessions with replacement
        n = len(sessions)
        indices = rng.choice(n, n, replace=True)
        resampled = [sessions[j] for j in indices]

        try:
            result = convex_inversion(
                resampled,
                period_init=period_init,
                pole_lambda_init=0.0,
                pole_beta_init=np.pi/4,
                n_facets=60,
                lambda_smooth=0.3,
                max_iter=30,
                seed=42 + i
            )
            periods.append(result.period)
            lambdas.append(np.degrees(result.pole_lambda))
            betas.append(np.degrees(result.pole_beta))
            rms_values.append(result.residual_rms)

            if (i + 1) % 10 == 0:
                log(f"  Bootstrap {i+1}/{N_BOOTSTRAP}: P={result.period:.4f}h rms={result.residual_rms:.4f}")
        except Exception:
            continue

    if len(periods) < 5:
        log(f"Only {len(periods)} successful bootstrap iterations, insufficient")
        return None

    # Compute statistics
    period_mean = np.mean(periods)
    period_std = np.std(periods)
    lambda_mean = np.mean(lambdas)
    lambda_std = np.std(lambdas)
    beta_mean = np.mean(betas)
    beta_std = np.std(betas)
    rms_mean = np.mean(rms_values)

    # Confidence flag
    if beta_std > 30 or period_std / period_mean > 0.01:
        confidence = "low"
    elif beta_std > 15 or period_std / period_mean > 0.005:
        confidence = "medium"
    else:
        confidence = "high"

    log(f"\n  Results ({len(periods)}/{N_BOOTSTRAP} successful):")
    log(f"  Period: {period_mean:.4f} +/- {period_std:.4f} h ({period_std/period_mean*100:.2f}%)")
    log(f"  Pole lambda: {lambda_mean:.1f} +/- {lambda_std:.1f} deg")
    log(f"  Pole beta: {beta_mean:.1f} +/- {beta_std:.1f} deg")
    log(f"  RMS: {rms_mean:.4f} mag")
    log(f"  Confidence: {confidence}")

    return {
        'asteroid': key,
        'number': number,
        'n_bootstrap': N_BOOTSTRAP,
        'n_successful': len(periods),
        'n_sessions': len(sessions),
        'n_datapoints': total_pts,
        'period': {
            'mean': float(period_mean),
            'std': float(period_std),
            'relative_uncertainty_pct': float(period_std / period_mean * 100),
            'min': float(np.min(periods)),
            'max': float(np.max(periods)),
        },
        'pole_lambda_deg': {
            'mean': float(lambda_mean),
            'std': float(lambda_std),
            'min': float(np.min(lambdas)),
            'max': float(np.max(lambdas)),
        },
        'pole_beta_deg': {
            'mean': float(beta_mean),
            'std': float(beta_std),
            'min': float(np.min(betas)),
            'max': float(np.max(betas)),
        },
        'residual_rms': {
            'mean': float(rms_mean),
            'std': float(np.std(rms_values)),
        },
        'confidence': confidence,
    }


def main():
    """Run bootstrap uncertainty quantification for all targets and save results."""
    all_results = {}
    for target in TARGETS:
        try:
            t0 = time.time()
            result = bootstrap_asteroid(target)
            if result:
                result['time_seconds'] = time.time() - t0
                all_results[target['key']] = result
        except Exception as e:
            log(f"\nERROR on {target['key']}: {e}")
            import traceback
            traceback.print_exc()

    out_path = os.path.join(RESULTS_DIR, 'uncertainty_quantification.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)

    log(f"\n{'='*60}")
    log("UNCERTAINTY SUMMARY")
    log(f"{'='*60}")
    for key, res in all_results.items():
        p = res['period']
        pl = res['pole_lambda_deg']
        pb = res['pole_beta_deg']
        log(f"  {key}: P={p['mean']:.4f}+/-{p['std']:.4f}h pole=({pl['mean']:.0f}+/-{pl['std']:.0f},{pb['mean']:.0f}+/-{pb['std']:.0f}) [{res['confidence']}]")

    log(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
