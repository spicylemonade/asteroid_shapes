"""
Period search engine using Lomb-Scargle periodogram and Phase Dispersion Minimization (PDM).

References:
  - Lomb (1976), Scargle (1982) [lomb1976, scargle1982 in sources.bib]
  - Stellingwerf (1978) [stellingwerf1978 in sources.bib]
"""
import numpy as np
from scipy.signal import lombscargle


def lomb_scargle_search(times, mags, errors=None, period_min=0.5, period_max=100.0,
                        n_periods=50000, top_n=10):
    """Lomb-Scargle periodogram for asteroid rotation period search.

    Asteroids typically have double-peaked lightcurves, so we search for
    P_rot = 2 * P_LS (the Lomb-Scargle period is half the rotation period).

    Parameters
    ----------
    times : array-like
        Julian dates of observations
    mags : array-like
        Magnitudes
    errors : array-like, optional
        Magnitude uncertainties
    period_min : float
        Minimum search period in hours
    period_max : float
        Maximum search period in hours
    n_periods : int
        Number of trial periods
    top_n : int
        Number of top candidates to return

    Returns
    -------
    list of dict
        Top candidates with keys: period_hours, power, confidence
    """
    times = np.asarray(times, dtype=np.float64)
    mags = np.asarray(mags, dtype=np.float64)

    # Convert times to hours
    t_hours = (times - times[0]) * 24.0

    # Normalize magnitudes
    mags_norm = mags - np.mean(mags)

    # Angular frequency range (rad/hour)
    freq_min = 2 * np.pi / period_max
    freq_max = 2 * np.pi / period_min

    # Trial frequencies
    freqs = np.linspace(freq_min, freq_max, n_periods)

    # Compute Lomb-Scargle periodogram
    power = lombscargle(t_hours, mags_norm, freqs, normalize=True)

    # Convert to periods
    periods = 2 * np.pi / freqs

    # Find peaks
    candidates = _find_peaks(periods, power, top_n)

    # For asteroids, check both P and 2P (double-peaked lightcurves)
    all_candidates = []
    for c in candidates:
        all_candidates.append({
            'period_hours': c['period'],
            'power': c['power'],
            'confidence': c['power'] / max(power) if max(power) > 0 else 0,
            'method': 'lomb_scargle',
        })
        # Also add 2*P candidate (rotation period for double-peaked LC)
        all_candidates.append({
            'period_hours': 2 * c['period'],
            'power': c['power'] * 0.95,  # Slightly lower priority
            'confidence': c['power'] / max(power) * 0.95 if max(power) > 0 else 0,
            'method': 'lomb_scargle_2x',
        })

    all_candidates.sort(key=lambda x: -x['confidence'])
    return all_candidates[:top_n]


def phase_dispersion_minimization(times, mags, errors=None, period_min=0.5,
                                   period_max=100.0, n_periods=20000, n_bins=10,
                                   top_n=10):
    """Phase Dispersion Minimization (PDM) for period search.

    PDM folds the data at trial periods and measures the dispersion within
    phase bins. The minimum dispersion corresponds to the best period.

    Parameters
    ----------
    times : array-like
        Julian dates
    mags : array-like
        Magnitudes
    errors : array-like, optional
        Uncertainties (used for weighting)
    period_min, period_max : float
        Search range in hours
    n_periods : int
        Number of trial periods
    n_bins : int
        Number of phase bins
    top_n : int
        Number of top candidates

    Returns
    -------
    list of dict
        Top candidates with period_hours, theta, confidence
    """
    times = np.asarray(times, dtype=np.float64)
    mags = np.asarray(mags, dtype=np.float64)

    t_hours = (times - times[0]) * 24.0
    total_var = np.var(mags)
    if total_var == 0:
        return []

    # Trial periods
    periods = np.linspace(period_min, period_max, n_periods)
    thetas = np.zeros(n_periods)

    for i, P in enumerate(periods):
        phases = (t_hours / P) % 1.0
        bin_edges = np.linspace(0, 1, n_bins + 1)

        bin_var_sum = 0.0
        n_total = 0
        for b in range(n_bins):
            mask = (phases >= bin_edges[b]) & (phases < bin_edges[b + 1])
            if np.sum(mask) > 1:
                bin_var_sum += np.var(mags[mask]) * np.sum(mask)
                n_total += np.sum(mask)

        if n_total > n_bins:
            thetas[i] = (bin_var_sum / n_total) / total_var
        else:
            thetas[i] = 1.0

    # Find minima (best periods have lowest theta)
    candidates = _find_valleys(periods, thetas, top_n)

    return [{
        'period_hours': c['period'],
        'theta': c['value'],
        'confidence': 1.0 - c['value'],  # Higher is better
        'method': 'pdm',
    } for c in candidates]


def combined_period_search(times, mags, errors=None, period_min=0.5,
                            period_max=100.0, top_n=5):
    """Combined period search using both Lomb-Scargle and PDM.

    Returns candidates ranked by combined score.
    """
    ls_candidates = lomb_scargle_search(times, mags, errors, period_min, period_max,
                                         top_n=top_n * 2)
    pdm_candidates = phase_dispersion_minimization(times, mags, errors, period_min,
                                                     period_max, top_n=top_n * 2)

    # Combine: for each LS candidate, check if PDM agrees within 1%
    combined = []
    for ls in ls_candidates:
        best_pdm_match = None
        for pdm in pdm_candidates:
            if abs(ls['period_hours'] - pdm['period_hours']) / ls['period_hours'] < 0.01:
                if best_pdm_match is None or pdm['confidence'] > best_pdm_match['confidence']:
                    best_pdm_match = pdm
                    break

        combined_score = ls['confidence']
        if best_pdm_match:
            combined_score = (ls['confidence'] + best_pdm_match['confidence']) / 2.0
            combined_score *= 1.5  # Boost for agreement

        combined.append({
            'period_hours': ls['period_hours'],
            'ls_confidence': ls['confidence'],
            'pdm_confidence': best_pdm_match['confidence'] if best_pdm_match else 0,
            'combined_score': min(combined_score, 1.0),
            'ls_pdm_agree': best_pdm_match is not None,
        })

    # Add PDM-only candidates
    for pdm in pdm_candidates:
        is_new = True
        for c in combined:
            if abs(c['period_hours'] - pdm['period_hours']) / pdm['period_hours'] < 0.01:
                is_new = False
                break
        if is_new:
            combined.append({
                'period_hours': pdm['period_hours'],
                'ls_confidence': 0,
                'pdm_confidence': pdm['confidence'],
                'combined_score': pdm['confidence'] * 0.8,
                'ls_pdm_agree': False,
            })

    combined.sort(key=lambda x: -x['combined_score'])
    return combined[:top_n]


def _find_peaks(x, y, top_n=10):
    """Find local maxima in y, return top N."""
    peaks = []
    for i in range(1, len(y) - 1):
        if y[i] > y[i - 1] and y[i] > y[i + 1]:
            peaks.append({'period': x[i], 'power': y[i]})
    peaks.sort(key=lambda p: -p['power'])
    return peaks[:top_n]


def _find_valleys(x, y, top_n=10):
    """Find local minima in y, return top N."""
    valleys = []
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] < y[i + 1]:
            valleys.append({'period': x[i], 'value': y[i]})
    valleys.sort(key=lambda v: v['value'])
    return valleys[:top_n]


# ============================================================
# Validation
# ============================================================

def merge_sessions_to_relative(blocks):
    """Merge multi-session lightcurve blocks into relative magnitude format.

    Each session has different comparison stars. Normalize each session
    to zero mean, then combine all data. This removes inter-session offsets
    while preserving intra-session variations due to rotation.

    Parameters
    ----------
    blocks : list of dict
        Each with 'data' (Nx3 array: JD, mag, err)

    Returns
    -------
    jd, mag, err : arrays
        Combined, sorted, mean-subtracted data
    """
    all_jd = []
    all_mag = []
    all_err = []

    for block in blocks:
        d = block['data']
        if len(d) < 3:
            continue
        # Subtract session mean to normalize
        session_mean = np.mean(d[:, 1])
        all_jd.extend(d[:, 0])
        all_mag.extend(d[:, 1] - session_mean)
        all_err.extend(d[:, 2])

    if not all_jd:
        return np.array([]), np.array([]), np.array([])

    jd = np.array(all_jd)
    mag = np.array(all_mag)
    err = np.array(all_err)

    # Sort by time
    order = np.argsort(jd)
    return jd[order], mag[order], err[order]


if __name__ == '__main__':
    import os
    import sys
    import json
    sys.path.insert(0, os.path.dirname(__file__))
    from data_ingest import preprocess_asteroid

    repo_root = os.path.dirname(os.path.dirname(__file__))
    zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
    gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')

    # Known periods from LCDB:
    # 433 Eros: 5.270 hours
    # 1036 Ganymed: 10.314 hours
    known_periods = {
        433: 5.270,
        1036: 10.314,
    }

    results = {}

    for ast_num, known_P in known_periods.items():
        print(f"\n{'='*60}")
        print(f"Asteroid {ast_num} - Known period: {known_P} hours")
        print(f"{'='*60}")

        data = preprocess_asteroid(zip_path, gz_path, ast_num)

        # Merge sessions with relative magnitude normalization
        all_jd, all_mag, all_err = merge_sessions_to_relative(data['blocks'])

        print(f"Total data points: {len(all_jd)}")
        print(f"Time span: {(all_jd[-1]-all_jd[0]):.1f} days")

        # Run combined search - focus on asteroid rotation periods (2-50 hours)
        candidates = combined_period_search(all_jd, all_mag, all_err,
                                             period_min=2.0, period_max=50.0,
                                             top_n=20)

        # Filter out obvious 24h aliases
        filtered = [c for c in candidates
                    if not (23.5 < c['period_hours'] < 24.5 or
                            47.0 < c['period_hours'] < 49.0)]
        if not filtered:
            filtered = candidates

        print(f"\nTop 10 period candidates (24h aliases removed):")
        best_match = None
        for i, c in enumerate(filtered[:10]):
            match = abs(c['period_hours'] - known_P) / known_P < 0.002
            flag = " <-- MATCH!" if match else ""
            if match and best_match is None:
                best_match = c
            print(f"  {i+1}. P={c['period_hours']:.4f} h, score={c['combined_score']:.3f}, "
                  f"LS={c['ls_confidence']:.3f}, PDM={c['pdm_confidence']:.3f}, "
                  f"agree={c['ls_pdm_agree']}{flag}")

        if best_match:
            diff = abs(best_match['period_hours'] - known_P)
            print(f"\nBest match: {best_match['period_hours']:.4f} h "
                  f"(error: {diff:.4f} h, {diff/known_P*100:.4f}%)")
            results[ast_num] = {
                'known_period': known_P,
                'found_period': best_match['period_hours'],
                'error_hours': diff,
                'status': 'PASS' if diff < 0.001 else 'CLOSE',
            }
        else:
            # Check half-period
            for c in filtered[:10]:
                if abs(c['period_hours'] * 2 - known_P) / known_P < 0.002:
                    diff = abs(c['period_hours'] * 2 - known_P)
                    print(f"\nHalf-period match found: {c['period_hours']:.4f} h "
                          f"(2x = {c['period_hours']*2:.4f} h vs known {known_P} h)")
                    results[ast_num] = {
                        'known_period': known_P,
                        'found_period': c['period_hours'] * 2,
                        'error_hours': diff,
                        'status': 'HALF_PERIOD_MATCH',
                    }
                    break
            else:
                print(f"\nWARNING: Known period {known_P} h not found in top candidates")
                diffs = [abs(c['period_hours'] - known_P) for c in filtered[:10]]
                closest_idx = np.argmin(diffs)
                closest = filtered[closest_idx]
                print(f"Closest: {closest['period_hours']:.4f} h "
                      f"(diff: {diffs[closest_idx]:.4f} h)")
                results[ast_num] = {
                    'known_period': known_P,
                    'found_period': closest['period_hours'],
                    'error_hours': diffs[closest_idx],
                    'status': 'NOT_FOUND',
                }

    # Save results
    os.makedirs(os.path.join(repo_root, 'results'), exist_ok=True)
    with open(os.path.join(repo_root, 'results', 'period_search_validation.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to results/period_search_validation.json")
