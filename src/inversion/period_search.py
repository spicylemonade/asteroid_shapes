"""
Period Search Module

Implements Lomb-Scargle periodogram and Phase Dispersion Minimization (PDM)
for asteroid rotation period determination from lightcurve data.
"""

import numpy as np
from scipy.signal import lombscargle


def lomb_scargle_period(times, mags, errors=None,
                        period_min=2.0, period_max=100.0,
                        n_frequencies=50000):
    """Compute Lomb-Scargle periodogram to find rotation period.

    Parameters
    ----------
    times : array, observation times (JD or hours)
    mags : array, magnitudes
    errors : array or None, magnitude uncertainties
    period_min : float, minimum period to search (same units as times)
    period_max : float, maximum period to search
    n_frequencies : int, number of frequency steps

    Returns
    -------
    periods : array, candidate periods sorted by power (descending)
    powers : array, corresponding normalized LS powers
    all_periods : array, full period grid
    all_powers : array, full power spectrum
    """
    times = np.asarray(times, dtype=np.float64)
    mags = np.asarray(mags, dtype=np.float64)

    # Subtract mean
    mags_centered = mags - np.mean(mags)

    # For asteroids, the lightcurve has 2 maxima/minima per rotation
    # So we search at twice the frequency (half the period)
    freq_min = 1.0 / period_max
    freq_max = 1.0 / period_min

    # Angular frequencies
    angular_freqs = np.linspace(2 * np.pi * freq_min,
                                2 * np.pi * freq_max,
                                n_frequencies)

    # Compute periodogram
    power = lombscargle(times, mags_centered, angular_freqs, normalize=True)

    # Convert to periods
    periods = 2 * np.pi / angular_freqs

    # Sort by power (descending)
    sort_idx = np.argsort(power)[::-1]

    return periods[sort_idx], power[sort_idx], periods, power


def phase_dispersion_minimization(times, mags, errors=None,
                                  period_min=2.0, period_max=100.0,
                                  n_periods=50000, n_bins=10):
    """Phase Dispersion Minimization (PDM) period search.

    Parameters
    ----------
    times : array, observation times
    mags : array, magnitudes
    errors : array or None, uncertainties
    period_min, period_max : float, period range to search
    n_periods : int, number of trial periods
    n_bins : int, number of phase bins

    Returns
    -------
    periods : array, candidate periods sorted by PDM statistic (ascending)
    theta : array, PDM theta statistic values (lower = better)
    all_periods : array, full period grid
    all_theta : array, full theta spectrum
    """
    times = np.asarray(times, dtype=np.float64)
    mags = np.asarray(mags, dtype=np.float64)

    total_variance = np.var(mags)
    if total_variance < 1e-30:
        trial_periods = np.linspace(period_min, period_max, n_periods)
        return trial_periods, np.ones(n_periods), trial_periods, np.ones(n_periods)

    trial_periods = np.linspace(period_min, period_max, n_periods)
    theta = np.zeros(n_periods)

    for i, P in enumerate(trial_periods):
        # Phase fold
        phases = (times / P) % 1.0

        # Bin the data
        bin_edges = np.linspace(0, 1, n_bins + 1)
        bin_variances = []
        bin_counts = []

        for b in range(n_bins):
            mask = (phases >= bin_edges[b]) & (phases < bin_edges[b + 1])
            n_in_bin = np.sum(mask)
            if n_in_bin > 1:
                bin_variances.append(np.var(mags[mask]) * n_in_bin)
                bin_counts.append(n_in_bin)

        if sum(bin_counts) > n_bins:
            pooled_var = sum(bin_variances) / (sum(bin_counts) - len(bin_counts))
            theta[i] = pooled_var / total_variance
        else:
            theta[i] = 1.0

    # Sort by theta (ascending - lower is better)
    sort_idx = np.argsort(theta)

    return trial_periods[sort_idx], theta[sort_idx], trial_periods, theta


def find_period(times, mags, errors=None,
                period_min=2.0, period_max=100.0,
                n_top=10, use_double_period=True):
    """Find rotation period using both LS and PDM methods.

    Parameters
    ----------
    times : array, observation times
    mags : array, magnitudes
    errors : array or None
    period_min, period_max : float, search range
    n_top : int, number of top candidates to return
    use_double_period : bool, if True, also search at double the LS period
        (asteroid lightcurves are double-peaked)

    Returns
    -------
    candidates : list of dicts with keys 'period', 'ls_power', 'pdm_theta', 'method'
    """
    # Lomb-Scargle
    ls_periods, ls_powers, ls_all_p, ls_all_pow = lomb_scargle_period(
        times, mags, errors, period_min, period_max)

    # PDM
    pdm_periods, pdm_theta, pdm_all_p, pdm_all_theta = phase_dispersion_minimization(
        times, mags, errors, period_min, period_max)

    candidates = []

    # Top LS periods
    for i in range(min(n_top, len(ls_periods))):
        p = ls_periods[i]
        candidates.append({
            'period': p,
            'ls_power': ls_powers[i],
            'pdm_theta': 1.0,
            'method': 'LS',
        })
        # Also try double period (full rotation)
        if use_double_period:
            p2 = 2.0 * p
            if period_min <= p2 <= period_max:
                candidates.append({
                    'period': p2,
                    'ls_power': ls_powers[i] * 0.9,
                    'pdm_theta': 1.0,
                    'method': 'LS_double',
                })

    # Top PDM periods
    for i in range(min(n_top, len(pdm_periods))):
        p = pdm_periods[i]
        candidates.append({
            'period': p,
            'ls_power': 0.0,
            'pdm_theta': pdm_theta[i],
            'method': 'PDM',
        })

    # Cross-validate: for each candidate, get PDM theta
    for cand in candidates:
        if cand['method'] != 'PDM':
            # Find closest PDM period
            idx = np.argmin(np.abs(pdm_all_p - cand['period']))
            cand['pdm_theta'] = pdm_all_theta[idx]
        if cand['method'] != 'LS' and cand['method'] != 'LS_double':
            idx = np.argmin(np.abs(ls_all_p - cand['period']))
            cand['ls_power'] = ls_all_pow[idx]

    # Score: combine LS power (higher better) and PDM theta (lower better)
    for cand in candidates:
        cand['score'] = cand['ls_power'] * (1.0 / max(cand['pdm_theta'], 0.01))

    # Sort by score descending
    candidates.sort(key=lambda x: x['score'], reverse=True)

    # Remove near-duplicates
    unique = []
    for cand in candidates:
        is_dup = False
        for u in unique:
            if abs(cand['period'] - u['period']) < 0.005 * u['period']:
                is_dup = True
                break
        if not is_dup:
            unique.append(cand)

    return unique[:n_top]
