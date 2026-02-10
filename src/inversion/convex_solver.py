"""
Convex Inversion Solver (Kaasalainen-Torppa Method)

Implements Levenberg-Marquardt optimization of spherical harmonics
coefficients and pole orientation to fit observed lightcurves.
"""

import numpy as np
from ..shapes.convex_model import (
    ConvexShapeModel, synthetic_lightcurve, spin_rotation_matrix,
    compute_brightness, coeffs_for_ellipsoid
)

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


def compute_chi_squared(model, lightcurves, pole_lambda, pole_beta, period,
                        phi0, t0, lambda_smooth=0.001, lambda_convex=0.0):
    """Compute chi-squared between observed and synthetic lightcurves.

    Parameters
    ----------
    model : ConvexShapeModel
    lightcurves : list of dicts, each with:
        'times': array of JDs
        'mags': array of observed magnitudes
        'errors': array of uncertainties
        'sun_dirs': array (N,3) unit vectors toward Sun (ecliptic)
        'obs_dirs': array (N,3) unit vectors toward observer (ecliptic)
        'phase_angles': array (N,) phase angles in radians
    pole_lambda, pole_beta : float, pole direction (radians)
    period : float, rotation period (days)
    phi0 : float, initial phase (radians)
    t0 : float, reference epoch (JD)
    lambda_smooth : float, smoothness regularization weight
    lambda_convex : float, convexity penalty weight

    Returns
    -------
    chi2 : float, total chi-squared
    chi2_data : float, data chi-squared only
    n_data : int, total number of data points
    """
    chi2_data = 0.0
    n_data = 0

    for lc in lightcurves:
        times = lc['times']
        obs_mags = lc['mags']
        errors = lc['errors']
        sun_dirs = lc['sun_dirs']
        obs_dirs = lc['obs_dirs']
        phase_angles = lc['phase_angles']

        # Compute synthetic lightcurve
        syn_mags = synthetic_lightcurve(
            model, times, pole_lambda, pole_beta, period, phi0, t0,
            sun_dirs, obs_dirs, phase_angles
        )

        # Per-session offset (analytically solve)
        offset = np.mean(obs_mags - syn_mags)

        # Chi-squared for this lightcurve
        residuals = (obs_mags - syn_mags - offset) / errors
        chi2_data += np.sum(residuals**2)
        n_data += len(times)

    # Smoothness regularization (penalize high-order coefficients)
    chi2_smooth = 0.0
    if lambda_smooth > 0:
        idx = 0
        for l in range(model.lmax + 1):
            for m in range(-l, l + 1):
                if idx < len(model.coeffs):
                    chi2_smooth += l**2 * (l + 1)**2 * model.coeffs[idx]**2
                idx += 1
        chi2_smooth *= lambda_smooth

    chi2 = chi2_data + chi2_smooth
    return chi2, chi2_data, n_data


def numerical_gradient(model, lightcurves, pole_lambda, pole_beta, period,
                       phi0, t0, lambda_smooth=0.001, eps=1e-6):
    """Compute gradient of chi-squared w.r.t. shape coefficients using finite differences."""
    n_params = model.n_coeffs
    grad = np.zeros(n_params)
    base_coeffs = model.coeffs.copy()

    chi2_0, _, _ = compute_chi_squared(
        model, lightcurves, pole_lambda, pole_beta, period, phi0, t0,
        lambda_smooth)

    for i in range(n_params):
        model.coeffs[i] = base_coeffs[i] + eps
        model._update_mesh()
        chi2_plus, _, _ = compute_chi_squared(
            model, lightcurves, pole_lambda, pole_beta, period, phi0, t0,
            lambda_smooth)
        grad[i] = (chi2_plus - chi2_0) / eps
        model.coeffs[i] = base_coeffs[i]

    model.coeffs = base_coeffs
    model._update_mesh()
    return grad


def convex_inversion(lightcurves, period, pole_lambda_init=None, pole_beta_init=None,
                     phi0=0.0, t0=2451545.0, lmax=6, subdivisions=3,
                     max_iter=200, lambda_smooth=0.001,
                     pole_search_grid=True, n_pole_grid=12,
                     verbose=True, seed=42):
    """Run convex inversion to determine shape and pole.

    Parameters
    ----------
    lightcurves : list of lightcurve dicts
    period : float, rotation period in days
    pole_lambda_init, pole_beta_init : float or None, initial pole guess (radians)
    phi0 : float, initial rotation phase
    t0 : float, reference epoch
    lmax : int, max spherical harmonics degree
    subdivisions : int, icosphere subdivisions
    max_iter : int, max iterations
    lambda_smooth : float, smoothness weight
    pole_search_grid : bool, search over pole grid
    n_pole_grid : int, number of grid points per axis
    verbose : bool
    seed : int, random seed

    Returns
    -------
    result : dict with keys 'model', 'pole_lambda', 'pole_beta', 'period',
             'phi0', 'chi2', 'chi2_reduced', 'converged'
    """
    np.random.seed(seed)

    best_chi2 = np.inf
    best_result = None

    # Grid search over pole orientations
    if pole_search_grid and pole_lambda_init is None:
        pole_lambdas = np.linspace(0, 2 * np.pi, n_pole_grid, endpoint=False)
        pole_betas = np.linspace(-np.pi / 2, np.pi / 2, n_pole_grid // 2 + 1)
        pole_grid = [(l, b) for l in pole_lambdas for b in pole_betas]
    elif pole_lambda_init is not None:
        pole_grid = [(pole_lambda_init, pole_beta_init or 0.0)]
    else:
        pole_grid = [(0.0, 0.0)]

    if verbose:
        print(f"Convex inversion: period={period*24:.4f}h, "
              f"lmax={lmax}, {len(pole_grid)} pole trials")

    for trial_idx, (plam, pbet) in enumerate(pole_grid):
        model = ConvexShapeModel(lmax=lmax, subdivisions=subdivisions)

        # Initialize with slightly elongated shape
        init_coeffs = coeffs_for_ellipsoid(1.3, 1.0, 0.9, lmax)
        model.set_coefficients(init_coeffs)

        # Levenberg-Marquardt style optimization
        damping = 1.0
        current_chi2, _, n_data = compute_chi_squared(
            model, lightcurves, plam, pbet, period, phi0, t0, lambda_smooth)

        for iteration in range(max_iter):
            grad = numerical_gradient(
                model, lightcurves, plam, pbet, period, phi0, t0,
                lambda_smooth)

            # Gradient descent step with damping
            step = -grad / (np.linalg.norm(grad) + damping)
            step_size = 0.01 / (1 + iteration * 0.01)

            new_coeffs = model.coeffs + step_size * step
            model.set_coefficients(new_coeffs)

            new_chi2, chi2_data, _ = compute_chi_squared(
                model, lightcurves, plam, pbet, period, phi0, t0,
                lambda_smooth)

            if new_chi2 < current_chi2:
                current_chi2 = new_chi2
                damping *= 0.5
            else:
                # Revert
                model.set_coefficients(model.coeffs - step_size * step)
                damping *= 2.0

            if damping > 1e6:
                break

        # Try small pole refinement
        for d_lam in [-0.05, 0, 0.05]:
            for d_bet in [-0.05, 0, 0.05]:
                test_lam = plam + d_lam
                test_bet = np.clip(pbet + d_bet, -np.pi/2, np.pi/2)
                test_chi2, _, _ = compute_chi_squared(
                    model, lightcurves, test_lam, test_bet, period, phi0, t0,
                    lambda_smooth)
                if test_chi2 < current_chi2:
                    current_chi2 = test_chi2
                    plam = test_lam
                    pbet = test_bet

        if current_chi2 < best_chi2:
            best_chi2 = current_chi2
            chi2_reduced = current_chi2 / max(n_data - model.n_coeffs - 4, 1)
            best_result = {
                'model': model,
                'coeffs': model.coeffs.copy(),
                'pole_lambda': plam,
                'pole_beta': pbet,
                'period': period,
                'phi0': phi0,
                'chi2': current_chi2,
                'chi2_reduced': chi2_reduced,
                'n_data': n_data,
                'converged': True,
            }

        if verbose and (trial_idx + 1) % max(1, len(pole_grid) // 10) == 0:
            print(f"  Trial {trial_idx+1}/{len(pole_grid)}: "
                  f"chi2={current_chi2:.2f}, best={best_chi2:.2f}")

    if verbose and best_result:
        print(f"Best solution: pole=({best_result['pole_lambda']*RAD2DEG:.1f}, "
              f"{best_result['pole_beta']*RAD2DEG:.1f}), "
              f"chi2_red={best_result['chi2_reduced']:.3f}")

    return best_result
