"""
Sparse Data Inversion Module

Implements inversion using sparse photometric points (individual magnitude
measurements, not dense lightcurves) following the methodology of
Durech et al. (2010) \cite{durech2009}.
"""

import numpy as np
from ..shapes.convex_model import (
    ConvexShapeModel, spin_rotation_matrix, compute_brightness,
    coeffs_for_ellipsoid
)

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


def color_correction(mag, filter_name, target_filter='V'):
    """Apply approximate color-index correction between filters.

    Based on typical asteroid color indices (Bowell et al. 1989).
    """
    corrections = {
        'B': -0.85,  # B-V typical for S-type
        'V': 0.0,
        'R': 0.38,   # V-R typical
        'I': 0.70,   # V-I typical
        'C': 0.0,    # Clear, approximate as V
        'G': -0.2,   # Gaia G band
        'g': -0.35,  # SDSS g
        'r': 0.25,   # SDSS r
        'i': 0.45,   # SDSS i
        'z': 0.55,   # SDSS z
    }
    corr = corrections.get(filter_name, 0.0) - corrections.get(target_filter, 0.0)
    return mag - corr


def sparse_chi_squared(model, sparse_data, pole_lambda, pole_beta, period,
                       phi0, t0, lambda_smooth=0.001):
    """Compute chi-squared for sparse photometric data.

    sparse_data : list of dicts with:
        'jd': float, Julian Date
        'mag': float, observed magnitude (calibrated)
        'error': float, uncertainty
        'sun_dir': array (3,), unit vector toward Sun (ecliptic)
        'obs_dir': array (3,), unit vector toward observer (ecliptic)
        'phase_angle': float, phase angle (radians)
        'filter': str, filter name
    """
    chi2 = 0.0
    n_data = len(sparse_data)

    for dp in sparse_data:
        R = spin_rotation_matrix(pole_lambda, pole_beta, period, phi0,
                                 dp['jd'], t0)
        sun_body = R @ dp['sun_dir']
        obs_body = R @ dp['obs_dir']

        flux = compute_brightness(model, sun_body, obs_body, dp['phase_angle'])
        flux = max(flux, 1e-30)
        mag_model = -2.5 * np.log10(flux)

        # Correct observed magnitude for filter difference
        mag_obs = color_correction(dp['mag'], dp.get('filter', 'V'))

        residual = (mag_obs - mag_model) / dp['error']
        chi2 += residual**2

    # Smoothness regularization
    chi2_smooth = 0.0
    if lambda_smooth > 0:
        idx = 0
        for l in range(model.lmax + 1):
            for m in range(-l, l + 1):
                if idx < len(model.coeffs):
                    chi2_smooth += l**2 * (l + 1)**2 * model.coeffs[idx]**2
                idx += 1
        chi2_smooth *= lambda_smooth

    return chi2 + chi2_smooth, chi2, n_data


def sparse_inversion(sparse_data, period, pole_lambda_init=None, pole_beta_init=None,
                     phi0=0.0, t0=2451545.0, lmax=6, subdivisions=3,
                     max_iter=150, lambda_smooth=0.005,
                     n_pole_grid=10, verbose=True, seed=42):
    """Run sparse data inversion following Durech et al. (2010).

    Parameters
    ----------
    sparse_data : list of dicts with sparse photometric points
    period : float, rotation period in days
    """
    np.random.seed(seed)

    best_chi2 = np.inf
    best_result = None

    # Pole grid
    if pole_lambda_init is None:
        pole_lambdas = np.linspace(0, 2*np.pi, n_pole_grid, endpoint=False)
        pole_betas = np.linspace(-np.pi/2, np.pi/2, n_pole_grid // 2 + 1)
        pole_grid = [(l, b) for l in pole_lambdas for b in pole_betas]
    else:
        pole_grid = [(pole_lambda_init, pole_beta_init or 0.0)]

    if verbose:
        print(f"Sparse inversion: {len(sparse_data)} points, "
              f"period={period*24:.4f}h, {len(pole_grid)} pole trials")

    for trial_idx, (plam, pbet) in enumerate(pole_grid):
        model = ConvexShapeModel(lmax=lmax, subdivisions=subdivisions)
        init_coeffs = coeffs_for_ellipsoid(1.2, 1.0, 0.9, lmax)
        model.set_coefficients(init_coeffs)

        current_chi2, _, n_data = sparse_chi_squared(
            model, sparse_data, plam, pbet, period, phi0, t0, lambda_smooth)

        # Gradient descent
        eps = 1e-5
        damping = 1.0
        for iteration in range(max_iter):
            grad = np.zeros(model.n_coeffs)
            base_coeffs = model.coeffs.copy()

            for j in range(model.n_coeffs):
                model.coeffs[j] = base_coeffs[j] + eps
                model._update_mesh()
                chi2_plus, _, _ = sparse_chi_squared(
                    model, sparse_data, plam, pbet, period, phi0, t0, lambda_smooth)
                grad[j] = (chi2_plus - current_chi2) / eps
                model.coeffs[j] = base_coeffs[j]

            model.coeffs = base_coeffs
            model._update_mesh()

            step = -grad / (np.linalg.norm(grad) + damping)
            step_size = 0.005 / (1 + iteration * 0.01)

            new_coeffs = model.coeffs + step_size * step
            model.set_coefficients(new_coeffs)

            new_chi2, _, _ = sparse_chi_squared(
                model, sparse_data, plam, pbet, period, phi0, t0, lambda_smooth)

            if new_chi2 < current_chi2:
                current_chi2 = new_chi2
                damping *= 0.7
            else:
                model.set_coefficients(base_coeffs)
                damping *= 2.0

            if damping > 1e6:
                break

        if current_chi2 < best_chi2:
            best_chi2 = current_chi2
            best_result = {
                'model': model,
                'coeffs': model.coeffs.copy(),
                'pole_lambda': plam,
                'pole_beta': pbet,
                'period': period,
                'phi0': phi0,
                'chi2': current_chi2,
                'chi2_reduced': current_chi2 / max(n_data - model.n_coeffs - 4, 1),
                'n_data': n_data,
                'converged': True,
            }

        if verbose and (trial_idx + 1) % max(1, len(pole_grid) // 5) == 0:
            print(f"  Trial {trial_idx+1}/{len(pole_grid)}: chi2={current_chi2:.2f}")

    if verbose and best_result:
        print(f"Best sparse solution: pole=({best_result['pole_lambda']*RAD2DEG:.1f}, "
              f"{best_result['pole_beta']*RAD2DEG:.1f}), "
              f"chi2_red={best_result['chi2_reduced']:.3f}")

    return best_result
