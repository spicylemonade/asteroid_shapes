"""
Hybrid dense+sparse data fusion module (ADAM-inspired).

Combines dense lightcurve data (ALCDEF) with sparse photometric data
in a unified inversion framework with configurable weighting.

References:
    Viikinkoski et al. 2015 (Viikinkoski2015 in sources.bib)
"""

import numpy as np
from scipy.optimize import minimize
from typing import List, Dict, Optional

from .forward_model import (
    create_triaxial_ellipsoid, compute_facet_properties,
    ecliptic_to_body_frame, lommel_seeliger
)
from .inversion import (
    _compute_model_lightcurves, build_adjacency, prepare_observations,
    InversionResult
)
from .sparse_inversion import sparse_chi_squared


def hybrid_chi_squared(params: np.ndarray,
                        normals: np.ndarray,
                        base_vertices: np.ndarray,
                        faces: np.ndarray,
                        dense_observations: List[dict],
                        dense_mags: List[np.ndarray],
                        dense_errs: List[np.ndarray],
                        sparse_data: List[Dict],
                        epoch_jd: float,
                        w_dense: float = 1.0,
                        w_sparse: float = 1.0,
                        lambda_smooth: float = 0.1,
                        adjacency: Optional[np.ndarray] = None) -> float:
    """Unified chi-squared combining dense and sparse data.

    params: [pole_lambda, pole_beta, period, H_abs, log_areas...]
    """
    n_facets = len(normals)
    pole_lambda = params[0]
    pole_beta = params[1]
    period = params[2]
    H_abs = params[3]
    log_areas = params[4:]
    areas = np.exp(log_areas)

    # Dense component
    chi2_dense = 0.0
    if dense_observations:
        model_mags = _compute_model_lightcurves(
            areas, normals, base_vertices, faces, dense_observations,
            pole_lambda, pole_beta, period, epoch_jd
        )
        for obs_m, mod_m, err in zip(dense_mags, model_mags, dense_errs):
            if len(obs_m) < 2:
                continue
            offset = np.mean(obs_m) - np.mean(mod_m)
            weights = 1.0 / np.maximum(err, 0.001) ** 2
            chi2_dense += np.sum(weights * (obs_m - mod_m - offset) ** 2)

    # Sparse component
    chi2_sparse = 0.0
    if sparse_data:
        period_days = period / 24.0
        for dp in sparse_data:
            t = dp['time']
            phase = 2.0 * np.pi * (t - epoch_jd) / period_days
            sun_body, obs_body = ecliptic_to_body_frame(
                dp['sun_ecl'], dp['obs_ecl'], pole_lambda, pole_beta, phase
            )
            mu0 = normals @ sun_body
            mu = normals @ obs_body
            visible = (mu0 > 0) & (mu > 0)
            if np.any(visible):
                scatter = lommel_seeliger(mu0[visible], mu[visible])
                brightness = np.sum(scatter * areas[visible])
            else:
                brightness = 1e-10
            model_mag = H_abs - 2.5 * np.log10(max(brightness, 1e-30))
            err = max(dp['uncertainty'], 0.01)
            chi2_sparse += ((model_mag - dp['reduced_mag']) / err) ** 2

    # Total with weighting
    chi2_total = w_dense * chi2_dense + w_sparse * chi2_sparse

    # Regularization
    if adjacency is not None and lambda_smooth > 0:
        for i, j in adjacency:
            chi2_total += lambda_smooth * (log_areas[i] - log_areas[j]) ** 2

    return chi2_total


def hybrid_inversion(sessions,
                      sparse_data: List[Dict] = None,
                      period_init: float = None,
                      pole_lambda_init: float = 0.0,
                      pole_beta_init: float = 0.0,
                      H_init: float = 10.0,
                      w_dense: float = 1.0,
                      w_sparse: float = 0.5,
                      n_facets: int = 150,
                      lambda_smooth: float = 0.3,
                      max_iter: int = 150,
                      seed: int = 42) -> InversionResult:
    """Run hybrid dense+sparse inversion.

    Iterative: sparse data constrains pole, dense data refines shape.
    """
    rng = np.random.RandomState(seed)

    # Prepare dense observations
    dense_obs, dense_mags, dense_errs = prepare_observations(sessions)

    if period_init is None:
        from .period_search import find_best_period, combine_lightcurve_sessions
        times, mags, errs = combine_lightcurve_sessions(sessions)
        period_init, _, _ = find_best_period(times, mags, errs, timeout_sec=120)

    # Create mesh
    n_lat = max(6, int(np.sqrt(n_facets / 2)))
    n_lon = max(12, int(n_facets / n_lat))
    vertices, faces = create_triaxial_ellipsoid(1.0, 1.0, 1.0, n_lat, n_lon)
    normals, areas, _ = compute_facet_properties(vertices, faces)
    adjacency = build_adjacency(faces)
    n_actual = len(faces)

    epoch_jd = dense_obs[0]['times'][0] if dense_obs else sparse_data[0]['time']

    log_areas_init = np.log(areas) + 0.02 * rng.randn(n_actual)
    params_init = np.concatenate([
        [pole_lambda_init, pole_beta_init, period_init, H_init],
        log_areas_init
    ])

    bounds = (
        [(0, 2 * np.pi), (-np.pi / 2, np.pi / 2),
         (max(1.0, period_init - 1.0), period_init + 1.0),
         (H_init - 5.0, H_init + 5.0)]
        + [(-4, 4)] * n_actual
    )

    result = minimize(
        hybrid_chi_squared,
        params_init,
        args=(normals, vertices, faces, dense_obs, dense_mags, dense_errs,
              sparse_data or [], epoch_jd, w_dense, w_sparse, lambda_smooth, adjacency),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': max_iter, 'ftol': 1e-8, 'disp': False}
    )

    opt = result.x
    opt_areas = np.exp(opt[4:])

    # Compute residuals
    model_mags = _compute_model_lightcurves(
        opt_areas, normals, vertices, faces, dense_obs,
        opt[0], opt[1], opt[2], epoch_jd
    )
    residuals = []
    for obs_m, mod_m in zip(dense_mags, model_mags):
        offset = np.mean(obs_m) - np.mean(mod_m)
        residuals.extend(obs_m - mod_m - offset)
    residual_rms = np.sqrt(np.mean(np.array(residuals) ** 2)) if residuals else 0.0

    return InversionResult(
        pole_lambda=opt[0], pole_beta=opt[1], period=opt[2],
        vertices=vertices, faces=faces, facet_areas=opt_areas,
        chi2=result.fun, residual_rms=residual_rms,
        n_iterations=result.nit, converged=result.success,
    )
