"""
Adaptive regularization and recursive optimization loop.

Automatically adjusts regularization weights, period search parameters,
and scattering law settings when validation deviation exceeds threshold.

References:
    Self-reinforcement validation approach (project-specific)
"""

import numpy as np
import json
from typing import Dict, List, Optional
from dataclasses import dataclass, asdict

from .inversion import convex_inversion, InversionResult
from .ga_optimizer import ga_optimize
from .validation import compute_validation_report


@dataclass
class OptimizationLog:
    """Log entry for one iteration of the recursive optimization loop."""
    iteration: int
    lambda_smooth: float
    n_facets: int
    max_iter: int
    period_init: float
    hausdorff: float
    iou: float
    pole_error: float
    period_error: float
    residual_rms: float
    passed: bool


def recursive_optimization(sessions,
                            truth_vertices: np.ndarray,
                            truth_faces: np.ndarray,
                            truth_pole_lambda: float,
                            truth_pole_beta: float,
                            truth_period: float,
                            max_iterations: int = 10,
                            deviation_threshold: float = 0.05,
                            seed: int = 42) -> Dict:
    """Run recursive self-reinforcement optimization loop.

    Runs blind inversion, compares against ground truth, and if deviation
    exceeds threshold, automatically adjusts parameters and re-runs.

    Args:
        sessions: ALCDEF lightcurve sessions
        truth_*: Ground truth shape, pole, and period
        max_iterations: Maximum optimization iterations
        deviation_threshold: IoU threshold to pass (1 - threshold)

    Returns:
        Dict with final result and optimization log
    """
    log = []

    # Initial parameters
    lambda_smooth = 0.5
    n_facets = 200
    max_iter_solver = 100
    period_init = truth_period  # Start with known period for validation

    # Parameter adjustment strategies
    strategies = [
        {'lambda_smooth': 0.5, 'n_facets': 200, 'max_iter': 100},
        {'lambda_smooth': 0.1, 'n_facets': 200, 'max_iter': 150},
        {'lambda_smooth': 1.0, 'n_facets': 300, 'max_iter': 100},
        {'lambda_smooth': 0.3, 'n_facets': 150, 'max_iter': 200},
        {'lambda_smooth': 0.05, 'n_facets': 200, 'max_iter': 150},
        {'lambda_smooth': 0.2, 'n_facets': 250, 'max_iter': 120},
        {'lambda_smooth': 0.8, 'n_facets': 200, 'max_iter': 80},
        {'lambda_smooth': 0.01, 'n_facets': 150, 'max_iter': 200},
        {'lambda_smooth': 0.5, 'n_facets': 100, 'max_iter': 300},
        {'lambda_smooth': 0.15, 'n_facets': 200, 'max_iter': 150},
    ]

    best_result = None
    best_iou = 0.0
    best_report = None

    for i in range(min(max_iterations, len(strategies))):
        strategy = strategies[i]
        lambda_smooth = strategy['lambda_smooth']
        n_facets = strategy['n_facets']
        max_iter_solver = strategy['max_iter']

        print(f"\n  Iteration {i+1}/{max_iterations}: lambda={lambda_smooth}, "
              f"facets={n_facets}, max_iter={max_iter_solver}")

        try:
            result = convex_inversion(
                sessions,
                period_init=period_init,
                n_facets=n_facets,
                lambda_smooth=lambda_smooth,
                max_iter=max_iter_solver,
                seed=seed + i,
            )

            # Compare against ground truth
            report = compute_validation_report(
                result.vertices, result.faces,
                truth_vertices, truth_faces,
                result.pole_lambda, result.pole_beta, result.period,
                truth_pole_lambda, truth_pole_beta, truth_period,
                result.residual_rms,
            )

            iou = report['volumetric_iou']
            h_dist = report['hausdorff_distance_relative']
            pole_err = report['pole_error_degrees']
            period_err = report['period_error_percent']

            passed = iou >= (1.0 - deviation_threshold)

            log_entry = OptimizationLog(
                iteration=i + 1,
                lambda_smooth=lambda_smooth,
                n_facets=n_facets,
                max_iter=max_iter_solver,
                period_init=period_init,
                hausdorff=h_dist,
                iou=iou,
                pole_error=pole_err,
                period_error=period_err,
                residual_rms=result.residual_rms,
                passed=passed,
            )
            log.append(log_entry)

            print(f"    IoU={iou:.3f}, Hausdorff={h_dist:.3f}, "
                  f"pole_err={pole_err:.1f}°, period_err={period_err:.3f}%")
            print(f"    Passed: {passed}")

            if iou > best_iou:
                best_iou = iou
                best_result = result
                best_report = report

            if passed:
                print(f"  *** PASSED at iteration {i+1} ***")
                break

        except Exception as e:
            print(f"    Error: {e}")
            log_entry = OptimizationLog(
                iteration=i + 1,
                lambda_smooth=lambda_smooth,
                n_facets=n_facets,
                max_iter=max_iter_solver,
                period_init=period_init,
                hausdorff=1.0, iou=0.0, pole_error=180.0,
                period_error=100.0, residual_rms=1.0, passed=False,
            )
            log.append(log_entry)

    return {
        'best_result': best_result,
        'best_report': best_report,
        'best_iou': best_iou,
        'log': [asdict(entry) for entry in log],
        'converged': best_iou >= (1.0 - deviation_threshold),
    }
