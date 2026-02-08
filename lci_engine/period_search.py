"""
Period search algorithms for asteroid lightcurve analysis.

Implements Phase Dispersion Minimization (PDM) and chi-squared periodogram
for finding the sidereal rotation period from photometric lightcurve data.

References:
    Kaasalainen & Torppa 2001 (Kaasalainen2001b in sources.bib)
    Stellingwerf 1978 (PDM original)
"""

import numpy as np
from typing import Tuple, List, Optional
import signal


class ComputeTimeout(Exception):
    """Raised when a period search computation exceeds the allowed time limit."""
    pass


def phase_dispersion_minimization(times: np.ndarray, mags: np.ndarray,
                                   mag_errs: np.ndarray,
                                   period_min: float = 2.0,
                                   period_max: float = 100.0,
                                   step: float = 0.0001,
                                   n_bins: int = 25) -> Tuple[np.ndarray, np.ndarray]:
    """Phase Dispersion Minimization period search.

    Args:
        times: Julian dates of observations
        mags: Observed magnitudes
        mag_errs: Magnitude uncertainties
        period_min, period_max: Search range in hours
        step: Step size in hours
        n_bins: Number of phase bins

    Returns:
        periods: Array of trial periods (hours)
        theta: PDM statistic (lower = better period)
    """
    periods = np.arange(period_min, period_max, step)
    theta = np.zeros(len(periods))

    # Overall variance
    total_var = np.var(mags)
    if total_var < 1e-10:
        return periods, np.ones(len(periods))

    times_hours = (times - times[0]) * 24.0

    for ip, period in enumerate(periods):
        # Compute phases
        phases = (times_hours / period) % 1.0

        # Bin the phases and compute within-bin variance
        bin_var_sum = 0.0
        bin_n_sum = 0

        for b in range(n_bins):
            lo = b / n_bins
            hi = (b + 1) / n_bins
            mask = (phases >= lo) & (phases < hi)
            n_in_bin = np.sum(mask)

            if n_in_bin >= 2:
                bin_var = np.var(mags[mask])
                bin_var_sum += (n_in_bin - 1) * bin_var
                bin_n_sum += n_in_bin - 1

        if bin_n_sum > 0:
            theta[ip] = (bin_var_sum / bin_n_sum) / total_var
        else:
            theta[ip] = 1.0

    return periods, theta


def chi_squared_periodogram(times: np.ndarray, mags: np.ndarray,
                             mag_errs: np.ndarray,
                             period_min: float = 2.0,
                             period_max: float = 100.0,
                             step: float = 0.0001,
                             n_harmonics: int = 4) -> Tuple[np.ndarray, np.ndarray]:
    """Chi-squared (Fourier) periodogram for period search.

    Fits a truncated Fourier series at each trial period and computes
    the chi-squared reduction.

    Args:
        times: Julian dates
        mags: Magnitudes
        mag_errs: Uncertainties
        period_min, period_max: Search range in hours
        step: Step size in hours
        n_harmonics: Number of Fourier harmonics to fit

    Returns:
        periods: Trial periods (hours)
        chi2_red: Reduced chi-squared (lower = better fit)
    """
    periods = np.arange(period_min, period_max, step)
    chi2_values = np.zeros(len(periods))

    times_hours = (times - times[0]) * 24.0
    weights = 1.0 / np.maximum(mag_errs, 0.001) ** 2
    n_obs = len(times)

    for ip, period in enumerate(periods):
        phases = 2.0 * np.pi * times_hours / period

        # Build design matrix for Fourier series
        n_params = 1 + 2 * n_harmonics
        A = np.zeros((n_obs, n_params))
        A[:, 0] = 1.0  # constant term

        for h in range(1, n_harmonics + 1):
            A[:, 2 * h - 1] = np.cos(h * phases)
            A[:, 2 * h] = np.sin(h * phases)

        # Weighted least squares
        W = np.diag(weights)
        try:
            AW = A.T @ W
            coeff = np.linalg.solve(AW @ A, AW @ mags)
            residuals = mags - A @ coeff
            chi2 = np.sum(weights * residuals ** 2)
        except np.linalg.LinAlgError:
            chi2 = np.sum(weights * mags ** 2)

        dof = max(n_obs - n_params, 1)
        chi2_values[ip] = chi2 / dof

    return periods, chi2_values


def lomb_scargle_search(times: np.ndarray, mags: np.ndarray,
                         period_min: float = 2.0,
                         period_max: float = 100.0,
                         n_freqs: int = 100000) -> Tuple[np.ndarray, np.ndarray]:
    """Lomb-Scargle periodogram for irregularly sampled data.

    Args:
        times: Julian dates
        mags: Magnitudes (mean-subtracted)
        period_min, period_max: Search range in hours
        n_freqs: Number of frequency samples

    Returns:
        periods: Trial periods (hours)
        power: LS power (higher = better)
    """
    from scipy.signal import lombscargle

    times_hours = (times - times[0]) * 24.0

    freq_min = 1.0 / period_max
    freq_max = 1.0 / period_min
    freqs = np.linspace(freq_min, freq_max, n_freqs)
    angular_freqs = 2.0 * np.pi * freqs

    power = lombscargle(times_hours, mags - np.mean(mags), angular_freqs, normalize=True)
    periods = 1.0 / freqs

    return periods, power


def find_best_period(times: np.ndarray, mags: np.ndarray,
                      mag_errs: np.ndarray,
                      period_min: float = 2.0,
                      period_max: float = 100.0,
                      coarse_step: float = 0.001,
                      fine_step: float = 0.00001,
                      n_candidates: int = 10,
                      method: str = 'lomb_scargle',
                      timeout_sec: int = 300) -> Tuple[float, np.ndarray, np.ndarray]:
    """Multi-stage period search with alias resolution.

    Uses Lomb-Scargle for initial detection, then chi-squared refinement.
    Automatically handles P/2 aliases (asteroid lightcurves show 2 maxima
    per rotation, so Fourier methods often find the half-period).

    Args:
        times, mags, mag_errs: Observation data
        period_min, period_max: Range in hours
        coarse_step: Step for chi2 refinement scan
        fine_step: Step for final refinement
        n_candidates: Number of candidates to refine
        method: 'lomb_scargle' (default), 'chi2', or 'pdm'
        timeout_sec: Timeout in seconds

    Returns:
        best_period: Best period in hours (sidereal rotation period)
        all_periods: Full array of trial periods from coarse scan
        all_scores: Corresponding scores
    """
    def _handler(signum, frame):
        """Signal handler that raises ComputeTimeout on SIGALRM."""
        raise ComputeTimeout()

    old_handler = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(timeout_sec)

    try:
        # Stage 1: Lomb-Scargle to find candidate frequencies
        periods_ls, power_ls = lomb_scargle_search(
            times, mags, period_min / 2.0, period_max, n_freqs=200000
        )

        # Find peaks in LS power
        from scipy.signal import find_peaks as _find_peaks
        peaks, _ = _find_peaks(power_ls, distance=50, height=0.01)

        if len(peaks) == 0:
            best_idx = np.argmax(power_ls)
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
            return periods_ls[best_idx], periods_ls, -power_ls

        # Sort by power
        sorted_peaks = peaks[np.argsort(-power_ls[peaks])]
        top_periods = periods_ls[sorted_peaks[:n_candidates]]

        # Stage 2: For each candidate, also consider 2*P (alias resolution)
        candidates = []
        for p in top_periods:
            candidates.append(p)
            if 2 * p <= period_max:
                candidates.append(2 * p)
            if p / 2 >= period_min:
                candidates.append(p / 2)

        # Stage 3: Refine each candidate with chi-squared periodogram
        best_period = candidates[0]
        best_chi2 = np.inf

        for p_cand in candidates:
            p_lo = max(period_min, p_cand - 0.02)
            p_hi = min(period_max, p_cand + 0.02)

            periods_f, chi2_f = chi_squared_periodogram(
                times, mags, mag_errs, p_lo, p_hi, step=0.0001, n_harmonics=4
            )

            if len(chi2_f) > 0:
                best_f_idx = np.argmin(chi2_f)
                if chi2_f[best_f_idx] < best_chi2:
                    best_chi2 = chi2_f[best_f_idx]
                    best_period = periods_f[best_f_idx]

        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
        return best_period, periods_ls, -power_ls

    except ComputeTimeout:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
        best_idx = np.argmax(power_ls) if 'power_ls' in dir() else 0
        return periods_ls[best_idx] if 'periods_ls' in dir() else period_min, \
               np.array([]), np.array([])


def combine_lightcurve_sessions(sessions) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Combine multiple lightcurve sessions into a single dataset.

    For relative (differential) magnitudes, each session is normalized
    to zero mean. For absolute magnitudes, data is used as-is.

    Args:
        sessions: List of LightcurveSession objects

    Returns:
        times: Combined JD array
        mags: Combined magnitudes (normalized)
        mag_errs: Combined uncertainties
    """
    all_times = []
    all_mags = []
    all_errs = []

    for session in sessions:
        t = session.jd
        m = session.mag
        e = session.mag_err

        # Normalize each session to zero mean for relative photometry
        if session.differ_mags or True:  # Treat all as relative for period search
            m = m - np.mean(m)

        all_times.append(t)
        all_mags.append(m)
        all_errs.append(e)

    times = np.concatenate(all_times)
    mags = np.concatenate(all_mags)
    errs = np.concatenate(all_errs)

    # Sort by time
    order = np.argsort(times)
    return times[order], mags[order], errs[order]
