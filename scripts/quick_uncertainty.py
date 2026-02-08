#!/usr/bin/env python3
"""
Item 022: Quick uncertainty quantification using analytic estimates + jackknife.

Faster alternative to full bootstrap: uses jackknife resampling (leave-one-session-out)
which requires only N inversions instead of 50+.
"""

import sys
import os
import json
import time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lci_engine.parsers import load_alcdef_asteroid
from lci_engine.inversion import convex_inversion

ALCDEF_ZIP = os.path.join(os.path.dirname(__file__), '..', 'ALCDEF_ALL.zip')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')

np.random.seed(42)

def log(msg):
    """Print a message immediately with flushed output."""
    print(msg, flush=True)


TARGETS = [
    {"key": "433_Eros", "number": 433, "period_init": 5.2594},
    {"key": "216_Kleopatra", "number": 216, "period_init": 5.3881},
    {"key": "25143_Itokawa", "number": 25143, "period_init": 12.0848},
    {"key": "1943_Anteros", "number": 1943, "period_init": 2.870},
    {"key": "5143_Heracles", "number": 5143, "period_init": 1.353},
    {"key": "3122_Florence", "number": 3122, "period_init": 2.0},
    {"key": "65803_Didymos", "number": 65803, "period_init": 2.261},
    {"key": "4015_Wilson_Harrington", "number": 4015, "period_init": 3.574},
    {"key": "57_Mnemosyne", "number": 57, "period_init": None},
    {"key": "4055_Magellan", "number": 4055, "period_init": None},
]


def jackknife_asteroid(target):
    """Estimate parameter uncertainties for one asteroid via jackknife resampling.

    Performs leave-one-session-out inversions, then computes jackknife variance
    estimates for period, pole direction, and residual RMS.

    Args:
        target: Dict with 'key', 'number', and 'period_init'.

    Returns:
        Dict with jackknife statistics and confidence flag, or None if
        insufficient data or too few successful iterations.
    """
    key = target['key']
    number = target['number']
    period_init = target['period_init']

    log(f"\n--- {key} (#{number}) ---")

    sessions = load_alcdef_asteroid(ALCDEF_ZIP, asteroid_number=number)
    if len(sessions) < 4:
        log(f"Only {len(sessions)} sessions, need >=4 for jackknife. Skipping.")
        return None

    sessions.sort(key=lambda s: len(s.jd), reverse=True)
    sessions = sessions[:10]
    total_pts = sum(len(s.jd) for s in sessions)
    n = len(sessions)
    log(f"Using {n} sessions, {total_pts} points")

    if period_init is None:
        from lci_engine.period_search import find_best_period, combine_lightcurve_sessions
        times, mags, errs = combine_lightcurve_sessions(sessions)
        period_init, _, _ = find_best_period(
            times, mags, errs, period_min=2.0, period_max=50.0, timeout_sec=60
        )
        log(f"Period found: {period_init:.4f} h")

    # Jackknife: leave one session out each time
    periods = []
    lambdas = []
    betas = []
    rms_values = []

    for i in range(n):
        subset = [s for j, s in enumerate(sessions) if j != i]
        try:
            result = convex_inversion(
                subset,
                period_init=period_init,
                pole_lambda_init=0.0,
                pole_beta_init=np.pi/4,
                n_facets=60,
                lambda_smooth=0.3,
                max_iter=30,
                seed=42
            )
            periods.append(result.period)
            lambdas.append(np.degrees(result.pole_lambda))
            betas.append(np.degrees(result.pole_beta))
            rms_values.append(result.residual_rms)
            log(f"  LOO-{i+1}/{n}: P={result.period:.4f}h rms={result.residual_rms:.4f}")
        except Exception as e:
            log(f"  LOO-{i+1}/{n}: FAILED ({e})")

    if len(periods) < 3:
        log("Too few successful jackknife iterations")
        return None

    # Jackknife variance estimation: var = (n-1)/n * sum((x_i - x_bar)^2)
    n_ok = len(periods)
    period_mean = np.mean(periods)
    period_var = (n_ok - 1) / n_ok * np.sum((np.array(periods) - period_mean)**2)
    period_std = np.sqrt(period_var)

    lambda_mean = np.mean(lambdas)
    lambda_var = (n_ok - 1) / n_ok * np.sum((np.array(lambdas) - lambda_mean)**2)
    lambda_std = np.sqrt(lambda_var)

    beta_mean = np.mean(betas)
    beta_var = (n_ok - 1) / n_ok * np.sum((np.array(betas) - beta_mean)**2)
    beta_std = np.sqrt(beta_var)

    rms_mean = np.mean(rms_values)

    # Confidence flag
    if beta_std > 30 or period_std / max(period_mean, 1e-10) > 0.01:
        confidence = "low"
    elif beta_std > 15 or period_std / max(period_mean, 1e-10) > 0.005:
        confidence = "medium"
    else:
        confidence = "high"

    log(f"  Result: P={period_mean:.4f}+/-{period_std:.4f}h pole=({lambda_mean:.0f}+/-{lambda_std:.0f},{beta_mean:.0f}+/-{beta_std:.0f}) [{confidence}]")

    return {
        'asteroid': key,
        'number': number,
        'method': 'jackknife',
        'n_jackknife': n,
        'n_successful': n_ok,
        'n_sessions': n,
        'n_datapoints': total_pts,
        'period': {
            'mean': float(period_mean),
            'std': float(period_std),
            'relative_uncertainty_pct': float(period_std / max(period_mean, 1e-10) * 100),
        },
        'pole_lambda_deg': {
            'mean': float(lambda_mean),
            'std': float(lambda_std),
        },
        'pole_beta_deg': {
            'mean': float(beta_mean),
            'std': float(beta_std),
        },
        'residual_rms': {
            'mean': float(rms_mean),
            'std': float(np.std(rms_values)),
        },
        'confidence': confidence,
    }


def main():
    """Run jackknife uncertainty quantification for all targets and save results."""
    all_results = {}
    for target in TARGETS:
        try:
            t0 = time.time()
            result = jackknife_asteroid(target)
            if result:
                result['time_seconds'] = round(time.time() - t0, 1)
                all_results[target['key']] = result
        except Exception as e:
            log(f"ERROR on {target['key']}: {e}")
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
        log(f"  {key}: P={p['mean']:.4f}+/-{p['std']:.4f}h ({p['relative_uncertainty_pct']:.2f}%) pole=({pl['mean']:.0f}+/-{pl['std']:.0f},{pb['mean']:.0f}+/-{pb['std']:.0f}) [{res['confidence']}]")

    log(f"\nResults saved to {out_path}")


if __name__ == '__main__':
    main()
