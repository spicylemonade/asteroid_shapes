"""
Sparse data inversion module for survey-quality photometry.

Extends convex inversion to handle sparse photometric data from surveys
like Gaia DR3, Pan-STARRS, and ZTF, with absolute magnitude calibration.

References:
    Durech et al. 2009 (Durech2009 in sources.bib)
    Cellino et al. 2009 (Cellino2009 in sources.bib)
"""

import numpy as np
from scipy.optimize import minimize
from typing import Tuple, List, Dict, Optional

from .forward_model import (
    create_triaxial_ellipsoid, compute_facet_properties,
    ecliptic_to_body_frame, lommel_seeliger
)
from .inversion import build_adjacency


def hg_phase_function(phase_angle: float, G: float = 0.15) -> float:
    """IAU H-G phase function for absolute magnitude calibration.

    Phi(alpha) = (1-G)*Phi1(alpha) + G*Phi2(alpha)

    Reference: Bowell et al. 1989
    """
    alpha = np.radians(phase_angle)

    # Simplified Bowell et al. approximation
    A1 = np.exp(-3.33 * (np.tan(alpha / 2.0)) ** 0.63)
    A2 = np.exp(-1.87 * (np.tan(alpha / 2.0)) ** 1.22)

    return (1.0 - G) * A1 + G * A2


def reduce_magnitude(V_obs: float, r_au: float, delta_au: float,
                      phase_deg: float, G: float = 0.15) -> float:
    """Convert observed magnitude to reduced magnitude.

    V_red = V_obs - 5*log10(r*delta) + 2.5*log10(Phi(alpha))

    Args:
        V_obs: Observed apparent magnitude
        r_au: Heliocentric distance (AU)
        delta_au: Geocentric distance (AU)
        phase_deg: Phase angle (degrees)
        G: Slope parameter
    """
    dist_correction = 5.0 * np.log10(r_au * delta_au)
    phase_correction = -2.5 * np.log10(max(hg_phase_function(phase_deg, G), 1e-10))
    return V_obs - dist_correction + phase_correction


def sparse_chi_squared(params: np.ndarray,
                        normals: np.ndarray,
                        faces: np.ndarray,
                        sparse_data: List[Dict],
                        lambda_smooth: float = 0.1,
                        adjacency: Optional[np.ndarray] = None) -> float:
    """Chi-squared for sparse photometric data.

    params: [pole_lambda, pole_beta, period, H_abs, log_areas...]

    Each sparse data point has: time, reduced_magnitude, sun_ecl, obs_ecl, uncertainty.
    """
    n_facets = len(normals)
    pole_lambda = params[0]
    pole_beta = params[1]
    period = params[2]
    H_abs = params[3]
    log_areas = params[4:]
    areas = np.exp(log_areas)

    period_days = period / 24.0

    # Vectorize: precompute static rotation matrix
    from .forward_model import rotation_matrix_z, rotation_matrix_y
    R_beta = rotation_matrix_y(-(np.pi / 2 - pole_beta))
    R_lambda = rotation_matrix_z(-pole_lambda)
    R_static = R_beta @ R_lambda

    t0 = sparse_data[0]['time']
    chi2 = 0.0

    for dp in sparse_data:
        t = dp['time']
        phase = 2.0 * np.pi * (t - t0) / period_days
        cos_p, sin_p = np.cos(-phase), np.sin(-phase)

        sun_s = R_static @ dp['sun_ecl']
        obs_s = R_static @ dp['obs_ecl']

        sun_body = np.array([cos_p*sun_s[0]-sin_p*sun_s[1],
                             sin_p*sun_s[0]+cos_p*sun_s[1], sun_s[2]])
        obs_body = np.array([cos_p*obs_s[0]-sin_p*obs_s[1],
                             sin_p*obs_s[0]+cos_p*obs_s[1], obs_s[2]])

        mu0 = normals @ sun_body
        mu = normals @ obs_body
        visible = (mu0 > 0) & (mu > 0)

        if np.any(visible):
            denom = np.maximum(mu0[visible] + mu[visible], 1e-10)
            brightness = np.sum(mu0[visible] / denom * areas[visible])
        else:
            brightness = 1e-10

        model_mag = H_abs - 2.5 * np.log10(max(brightness, 1e-30))
        err = max(dp['uncertainty'], 0.01)
        chi2 += ((model_mag - dp['reduced_mag']) / err) ** 2

    # Smoothness regularization (vectorized)
    if adjacency is not None and len(adjacency) > 0 and lambda_smooth > 0:
        diffs = log_areas[adjacency[:, 0]] - log_areas[adjacency[:, 1]]
        chi2 += lambda_smooth * np.sum(diffs ** 2)

    return chi2


def sparse_inversion(sparse_data: List[Dict],
                      period_init: float,
                      pole_lambda_init: float = 0.0,
                      pole_beta_init: float = 0.0,
                      H_init: float = 10.0,
                      n_facets: int = 100,
                      lambda_smooth: float = 0.5,
                      max_iter: int = 100,
                      seed: int = 42) -> Dict:
    """Run sparse data inversion.

    Args:
        sparse_data: List of dicts with keys: time, reduced_mag, sun_ecl, obs_ecl, uncertainty
        period_init: Initial period guess (hours)
        pole_lambda_init, pole_beta_init: Initial pole (radians)
        H_init: Initial absolute magnitude
        n_facets: Number of mesh facets
        lambda_smooth: Smoothness regularization
        max_iter: Maximum iterations
        seed: Random seed

    Returns:
        Dict with results: period, pole, shape, etc.
    """
    rng = np.random.RandomState(seed)

    n_lat = max(5, int(np.sqrt(n_facets / 2)))
    n_lon = max(10, int(n_facets / n_lat))
    vertices, faces = create_triaxial_ellipsoid(1.0, 1.0, 1.0, n_lat, n_lon)
    normals, areas, _ = compute_facet_properties(vertices, faces)
    adjacency = build_adjacency(faces)
    n_actual = len(faces)

    log_areas_init = np.log(areas) + 0.01 * rng.randn(n_actual)
    params_init = np.concatenate([
        [pole_lambda_init, pole_beta_init, period_init, H_init],
        log_areas_init
    ])

    bounds = (
        [(0, 2 * np.pi), (-np.pi / 2, np.pi / 2),
         (max(1.0, period_init - 2.0), period_init + 2.0),
         (H_init - 3.0, H_init + 3.0)]
        + [(-3, 3)] * n_actual
    )

    result = minimize(
        sparse_chi_squared,
        params_init,
        args=(normals, faces, sparse_data, lambda_smooth, adjacency),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': max_iter, 'ftol': 1e-6, 'disp': False}
    )

    opt_lambda = result.x[0]
    opt_beta = result.x[1]
    opt_period = result.x[2]
    opt_H = result.x[3]
    opt_areas = np.exp(result.x[4:])

    return {
        'pole_lambda': opt_lambda,
        'pole_beta': opt_beta,
        'period': opt_period,
        'H_abs': opt_H,
        'facet_areas': opt_areas,
        'vertices': vertices,
        'faces': faces,
        'chi2': result.fun,
        'converged': result.success,
        'n_iterations': result.nit,
    }


def subsample_to_sparse(sessions, n_points_per_apparition: int = 50,
                          n_apparitions: int = 3, seed: int = 42) -> List[Dict]:
    """Subsample dense lightcurve sessions to simulate sparse survey data.

    Groups sessions by approximate apparition (time clusters), then randomly
    samples points from each apparition.
    """
    rng = np.random.RandomState(seed)

    # Collect all data points
    all_points = []
    for s in sessions:
        for j in range(len(s.jd)):
            all_points.append({
                'time': s.jd[j],
                'mag': s.mag[j],
                'uncertainty': s.mag_err[j],
                'phase_angle': s.phase_angle,
                'pabl': s.pabl,
                'pabb': s.pabb,
            })

    if not all_points:
        return []

    # Sort by time
    all_points.sort(key=lambda x: x['time'])
    times = np.array([p['time'] for p in all_points])

    # Cluster into apparitions (gaps > 30 days)
    apparitions = []
    current = [0]
    for i in range(1, len(times)):
        if times[i] - times[i - 1] > 30:
            apparitions.append(current)
            current = [i]
        else:
            current.append(i)
    apparitions.append(current)

    # Select up to n_apparitions
    if len(apparitions) > n_apparitions:
        selected = rng.choice(len(apparitions), n_apparitions, replace=False)
        apparitions = [apparitions[i] for i in sorted(selected)]

    # Subsample from each apparition
    sparse_data = []
    for app_indices in apparitions:
        n_sample = min(n_points_per_apparition, len(app_indices))
        sampled = rng.choice(app_indices, n_sample, replace=False)

        for idx in sampled:
            p = all_points[idx]
            # Create approximate geometry
            phase_rad = np.radians(abs(p['phase_angle']))
            pabl_rad = np.radians(p['pabl'])
            pabb_rad = np.radians(p['pabb'])

            pab = np.array([
                np.cos(pabb_rad) * np.cos(pabl_rad),
                np.cos(pabb_rad) * np.sin(pabl_rad),
                np.sin(pabb_rad)
            ])

            half_phase = phase_rad / 2.0
            sun_ecl = pab * np.cos(half_phase) + np.array([0, 0, 1]) * np.sin(half_phase)
            obs_ecl = pab * np.cos(half_phase) - np.array([0, 0, 1]) * np.sin(half_phase)
            sun_ecl /= np.linalg.norm(sun_ecl)
            obs_ecl /= np.linalg.norm(obs_ecl)

            sparse_data.append({
                'time': p['time'],
                'reduced_mag': p['mag'],  # Approximate: use as reduced mag
                'sun_ecl': sun_ecl,
                'obs_ecl': obs_ecl,
                'uncertainty': p['uncertainty'],
            })

    return sparse_data


def bootstrap_uncertainty(sparse_data: List[Dict],
                           period_init: float,
                           pole_lambda_init: float,
                           pole_beta_init: float,
                           H_init: float,
                           n_bootstrap: int = 50,
                           seed: int = 42) -> Dict:
    """Estimate parameter uncertainties via bootstrap resampling.

    Returns:
        Dict with mean, std for period, pole_lambda, pole_beta
    """
    rng = np.random.RandomState(seed)
    n = len(sparse_data)

    periods = []
    lambdas = []
    betas = []

    for i in range(n_bootstrap):
        # Resample with replacement
        indices = rng.choice(n, n, replace=True)
        resampled = [sparse_data[j] for j in indices]

        try:
            result = sparse_inversion(
                resampled, period_init,
                pole_lambda_init, pole_beta_init, H_init,
                n_facets=50, lambda_smooth=1.0, max_iter=30, seed=seed + i
            )
            periods.append(result['period'])
            lambdas.append(result['pole_lambda'])
            betas.append(result['pole_beta'])
        except Exception:
            continue

    if len(periods) < 3:
        return {
            'period_mean': period_init, 'period_std': 0,
            'pole_lambda_mean': pole_lambda_init, 'pole_lambda_std': 0,
            'pole_beta_mean': pole_beta_init, 'pole_beta_std': 0,
            'n_successful': len(periods),
        }

    return {
        'period_mean': np.mean(periods),
        'period_std': np.std(periods),
        'pole_lambda_mean': np.mean(lambdas),
        'pole_lambda_std': np.std(lambdas),
        'pole_beta_mean': np.mean(betas),
        'pole_beta_std': np.std(betas),
        'n_successful': len(periods),
    }
