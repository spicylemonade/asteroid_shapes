"""
Convex lightcurve inversion solver (Kaasalainen-Torppa method).

Implements gradient-based minimization of chi-squared between observed
and model lightcurves, optimizing spin axis, period, and convex shape.

No reliance on external inversion libraries for core logic.
numpy/scipy for linear algebra permitted.

References:
    Kaasalainen & Torppa 2001 (Kaasalainen2001a)
    Kaasalainen, Torppa & Muinonen 2001 (Kaasalainen2001b)
"""

import numpy as np
from scipy.optimize import minimize
from typing import Tuple, List, Optional
from dataclasses import dataclass

from .forward_model import (
    create_triaxial_ellipsoid, compute_facet_properties,
    ecliptic_to_body_frame, lommel_seeliger,
    compute_synthetic_brightness, brightness_to_magnitude,
    save_mesh_obj
)


@dataclass
class InversionResult:
    """Result of lightcurve inversion."""
    pole_lambda: float  # ecliptic longitude (radians)
    pole_beta: float  # ecliptic latitude (radians)
    period: float  # sidereal period (hours)
    vertices: np.ndarray  # shape mesh vertices
    faces: np.ndarray  # shape mesh faces
    facet_areas: np.ndarray  # optimized facet areas
    chi2: float  # final chi-squared
    residual_rms: float  # RMS of magnitude residuals
    n_iterations: int
    converged: bool


def _compute_model_lightcurves(facet_areas: np.ndarray,
                                normals: np.ndarray,
                                base_vertices: np.ndarray,
                                faces: np.ndarray,
                                observations: List[dict],
                                pole_lambda: float,
                                pole_beta: float,
                                period_hours: float,
                                epoch_jd: float) -> List[np.ndarray]:
    """Compute model lightcurves for all observation sessions.

    Vectorized over rotation phases for speed.
    """
    period_days = period_hours / 24.0
    model_mags = []

    # Precompute the static part of rotation matrix (pole alignment)
    from .forward_model import rotation_matrix_z, rotation_matrix_y
    R_beta = rotation_matrix_y(-(np.pi / 2 - pole_beta))
    R_lambda = rotation_matrix_z(-pole_lambda)
    R_static = R_beta @ R_lambda  # (3, 3)

    for obs in observations:
        times = obs['times']
        sun_ecl = obs['sun_ecl']
        obs_ecl = obs['obs_ecl']

        # Pre-rotate sun/obs by static matrix
        sun_static = R_static @ sun_ecl  # (3,)
        obs_static = R_static @ obs_ecl  # (3,)

        phases = 2.0 * np.pi * (times - epoch_jd) / period_days
        n_times = len(phases)

        # Vectorized spin rotation: compute cos/sin for all phases
        cos_p = np.cos(-phases)
        sin_p = np.sin(-phases)

        # For each time, R_spin @ sun_static where R_spin rotates around Z
        # R_spin = [[cos, -sin, 0], [sin, cos, 0], [0, 0, 1]]
        sun_body_x = cos_p * sun_static[0] - sin_p * sun_static[1]
        sun_body_y = sin_p * sun_static[0] + cos_p * sun_static[1]
        sun_body_z = np.full(n_times, sun_static[2])
        # sun_body shape: (n_times, 3)
        sun_body_all = np.column_stack([sun_body_x, sun_body_y, sun_body_z])

        obs_body_x = cos_p * obs_static[0] - sin_p * obs_static[1]
        obs_body_y = sin_p * obs_static[0] + cos_p * obs_static[1]
        obs_body_z = np.full(n_times, obs_static[2])
        obs_body_all = np.column_stack([obs_body_x, obs_body_y, obs_body_z])

        # Compute mu0 and mu for all times and all facets at once
        # normals: (M, 3), sun_body_all: (n_times, 3)
        # mu0: (n_times, M) = sun_body_all @ normals.T
        mu0_all = sun_body_all @ normals.T  # (n_times, M)
        mu_all = obs_body_all @ normals.T   # (n_times, M)

        # Lommel-Seeliger: mu0 / (mu0 + mu), only where both > 0
        visible = (mu0_all > 0) & (mu_all > 0)
        denom = mu0_all + mu_all
        denom = np.where(visible, np.maximum(denom, 1e-10), 1.0)
        scatter = np.where(visible, mu0_all / denom, 0.0)

        # brightness = sum over facets of scatter * area
        brightness = scatter @ facet_areas  # (n_times,)
        brightness = np.maximum(brightness, 1e-30)
        model_mags.append(-2.5 * np.log10(brightness))

    return model_mags


def _chi_squared(params: np.ndarray,
                  normals: np.ndarray,
                  base_vertices: np.ndarray,
                  faces: np.ndarray,
                  observations: List[dict],
                  observed_mags: List[np.ndarray],
                  observed_errs: List[np.ndarray],
                  epoch_jd: float,
                  lambda_smooth: float = 0.1,
                  adjacency: Optional[np.ndarray] = None) -> float:
    """Compute total chi-squared loss for inversion.

    params layout: [pole_lambda, pole_beta, period, log_areas...]
    """
    n_facets = len(normals)
    pole_lambda = params[0]
    pole_beta = params[1]
    period = params[2]
    log_areas = params[3:]

    # Ensure positive areas
    areas = np.exp(log_areas)

    # Compute model lightcurves
    model_mags = _compute_model_lightcurves(
        areas, normals, base_vertices, faces,
        observations, pole_lambda, pole_beta, period, epoch_jd
    )

    # Chi-squared: per-session with scaling factor
    chi2 = 0.0
    for i, (obs_m, mod_m, err) in enumerate(zip(observed_mags, model_mags, observed_errs)):
        if len(obs_m) < 2:
            continue
        # Optimal scaling: shift model to match observed mean
        offset = np.mean(obs_m) - np.mean(mod_m)
        mod_shifted = mod_m + offset
        weights = 1.0 / np.maximum(err, 0.001) ** 2
        chi2 += np.sum(weights * (obs_m - mod_shifted) ** 2)

    # Regularization: smoothness of facet areas (vectorized)
    if adjacency is not None and len(adjacency) > 0 and lambda_smooth > 0:
        diffs = log_areas[adjacency[:, 0]] - log_areas[adjacency[:, 1]]
        chi2 += lambda_smooth * np.sum(diffs ** 2)

    return chi2


def build_adjacency(faces: np.ndarray) -> np.ndarray:
    """Build adjacency list of face pairs sharing an edge."""
    from collections import defaultdict

    edge_to_face = defaultdict(list)
    for fi, face in enumerate(faces):
        for k in range(3):
            edge = tuple(sorted([face[k], face[(k + 1) % 3]]))
            edge_to_face[edge].append(fi)

    pairs = []
    for edge, face_list in edge_to_face.items():
        for i in range(len(face_list)):
            for j in range(i + 1, len(face_list)):
                pairs.append([face_list[i], face_list[j]])

    return np.array(pairs) if pairs else np.empty((0, 2), dtype=int)


def prepare_observations(sessions, default_sun_ecl=None, default_obs_ecl=None):
    """Convert ALCDEF sessions to observation dicts for inversion.

    Since ALCDEF provides phase angle but not full 3D geometry,
    we approximate sun/observer vectors from the Phase Angle Bisector.
    """
    observations = []
    observed_mags = []
    observed_errs = []

    for session in sessions:
        if len(session.jd) < 3:
            continue

        # Approximate geometry from PAB and phase angle
        phase_rad = np.radians(abs(session.phase_angle))
        pabl_rad = np.radians(session.pabl)
        pabb_rad = np.radians(session.pabb)

        # PAB direction (midpoint between sun and observer as seen from asteroid)
        pab = np.array([
            np.cos(pabb_rad) * np.cos(pabl_rad),
            np.cos(pabb_rad) * np.sin(pabl_rad),
            np.sin(pabb_rad)
        ])

        # Sun and observer are offset by phase_angle/2 from PAB
        # Simplified: assume both in same plane
        half_phase = phase_rad / 2.0
        sun_ecl = pab * np.cos(half_phase) + np.array([0, 0, 1]) * np.sin(half_phase)
        obs_ecl = pab * np.cos(half_phase) - np.array([0, 0, 1]) * np.sin(half_phase)
        sun_ecl /= np.linalg.norm(sun_ecl)
        obs_ecl /= np.linalg.norm(obs_ecl)

        observations.append({
            'times': session.jd,
            'sun_ecl': sun_ecl,
            'obs_ecl': obs_ecl,
        })

        # Normalize magnitudes per session
        mags = session.mag - np.mean(session.mag)
        observed_mags.append(mags)
        observed_errs.append(session.mag_err)

    return observations, observed_mags, observed_errs


def convex_inversion(sessions,
                      period_init: float = None,
                      pole_lambda_init: float = None,
                      pole_beta_init: float = None,
                      n_facets: int = 200,
                      lambda_smooth: float = 0.5,
                      max_iter: int = 200,
                      seed: int = 42) -> InversionResult:
    """Run convex lightcurve inversion.

    Args:
        sessions: List of LightcurveSession objects
        period_init: Initial period guess (hours). If None, run period search.
        pole_lambda_init: Initial pole longitude (radians)
        pole_beta_init: Initial pole latitude (radians)
        n_facets: Number of mesh facets
        lambda_smooth: Smoothness regularization weight
        max_iter: Maximum L-BFGS iterations
        seed: Random seed

    Returns:
        InversionResult with optimized parameters
    """
    rng = np.random.RandomState(seed)

    # Prepare observations
    observations, observed_mags, observed_errs = prepare_observations(sessions)

    if not observations:
        raise ValueError("No valid observations for inversion")

    # Period search if needed
    if period_init is None:
        from .period_search import find_best_period, combine_lightcurve_sessions
        times, mags, errs = combine_lightcurve_sessions(sessions)
        period_init, _, _ = find_best_period(times, mags, errs, timeout_sec=120)

    # Initial pole estimate (if not given, try a grid)
    if pole_lambda_init is None:
        pole_lambda_init = 0.0
    if pole_beta_init is None:
        pole_beta_init = np.pi / 4.0

    # Create initial mesh
    n_lat = max(8, int(np.sqrt(n_facets / 2)))
    n_lon = max(16, int(n_facets / n_lat))
    vertices, faces = create_triaxial_ellipsoid(1.0, 1.0, 1.0, n_lat, n_lon)
    normals, areas, centroids = compute_facet_properties(vertices, faces)
    n_actual_facets = len(faces)

    # Build adjacency for regularization
    adjacency = build_adjacency(faces)

    # Reference epoch
    epoch_jd = observations[0]['times'][0]

    # Initial parameters
    log_areas_init = np.log(areas) + 0.05 * rng.randn(n_actual_facets)
    params_init = np.concatenate([
        [pole_lambda_init, pole_beta_init, period_init],
        log_areas_init
    ])

    # Bounds: pole lambda [0, 2pi], beta [-pi/2, pi/2], period [P-1, P+1]
    bounds = (
        [(0, 2 * np.pi), (-np.pi / 2, np.pi / 2),
         (max(1.0, period_init - 1.0), period_init + 1.0)]
        + [(-5, 5)] * n_actual_facets
    )

    # Optimize
    result = minimize(
        _chi_squared,
        params_init,
        args=(normals, vertices, faces, observations, observed_mags,
              observed_errs, epoch_jd, lambda_smooth, adjacency),
        method='L-BFGS-B',
        bounds=bounds,
        options={'maxiter': max_iter, 'ftol': 1e-8, 'disp': False}
    )

    # Extract results
    opt_lambda = result.x[0]
    opt_beta = result.x[1]
    opt_period = result.x[2]
    opt_log_areas = result.x[3:]
    opt_areas = np.exp(opt_log_areas)

    # Scale vertices by optimized areas
    normals_unit, base_areas, _ = compute_facet_properties(vertices, faces)
    scale_factors = np.sqrt(opt_areas / np.maximum(base_areas, 1e-10))

    # Apply area scaling to vertices (approximate: scale along normal direction)
    new_vertices = vertices.copy()
    vertex_scales = np.ones(len(vertices))
    for fi, face in enumerate(faces):
        for vi in face:
            vertex_scales[vi] *= scale_factors[fi] ** (1.0 / 3.0)

    # Normalize
    vertex_scales /= np.mean(vertex_scales)
    norms = np.linalg.norm(new_vertices, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    directions = new_vertices / norms
    new_vertices = directions * (norms.ravel() * vertex_scales)[:, np.newaxis]

    # Compute residual RMS
    model_mags = _compute_model_lightcurves(
        opt_areas, normals, vertices, faces, observations,
        opt_lambda, opt_beta, opt_period, epoch_jd
    )

    residuals = []
    for obs_m, mod_m in zip(observed_mags, model_mags):
        offset = np.mean(obs_m) - np.mean(mod_m)
        residuals.extend(obs_m - mod_m - offset)

    residual_rms = np.sqrt(np.mean(np.array(residuals) ** 2))

    return InversionResult(
        pole_lambda=opt_lambda,
        pole_beta=opt_beta,
        period=opt_period,
        vertices=new_vertices,
        faces=faces,
        facet_areas=opt_areas,
        chi2=result.fun,
        residual_rms=residual_rms,
        n_iterations=result.nit,
        converged=result.success,
    )
