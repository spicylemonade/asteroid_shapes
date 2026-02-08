#!/usr/bin/env python3
"""
Run asteroid lightcurve inversion pipeline on top 10 candidate asteroids.

Item 020 of the research rubric: Execute the validated inversion pipeline
on the top 10 candidates from the target selection list.

For each candidate:
  1. Load ALCDEF data
  2. Perform period search (timeout 120s, range 2-50h)
  3. Run convex inversion with 2 pole trials
  4. Save best .obj mesh
  5. Record period, pole, chi2, rms, converged status

Outputs:
  - results/candidate_NNNN_Name.obj  (one per asteroid)
  - results/candidate_inversion_results.json  (full results)
  - results/candidates_modeled.csv  (summary table)
"""

import sys
import os
import json
import csv
import time
import traceback
import signal
import numpy as np

# Add repo root to path
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, REPO_ROOT)

from lci_engine.parsers import load_alcdef_asteroid
from lci_engine.period_search import find_best_period, combine_lightcurve_sessions
from lci_engine.inversion import convex_inversion
from lci_engine.forward_model import save_mesh_obj

np.random.seed(42)

ALCDEF_ZIP = os.path.join(REPO_ROOT, 'ALCDEF_ALL.zip')
RESULTS_DIR = os.path.join(REPO_ROOT, 'results')

# Top 10 candidates from target_candidates_top50.csv
TOP10_CANDIDATES = [
    {"rank": 1,  "number": 1943,  "name": "Anteros",        "num_lcs": 117},
    {"rank": 2,  "number": 5143,  "name": "Heracles",       "num_lcs": 83},
    {"rank": 3,  "number": 3122,  "name": "Florence",       "num_lcs": 72},
    {"rank": 4,  "number": 65803, "name": "Didymos",        "num_lcs": 70},
    {"rank": 5,  "number": 4015,  "name": "Wilson-Harrington", "num_lcs": 64},
    {"rank": 6,  "number": 57,    "name": "Mnemosyne",      "num_lcs": 64},
    {"rank": 7,  "number": 4055,  "name": "Magellan",       "num_lcs": 62},
    {"rank": 8,  "number": 185,   "name": "Eunike",         "num_lcs": 56},
    {"rank": 9,  "number": 13553, "name": "Masaakikoyama",  "num_lcs": 55},
    {"rank": 10, "number": 887,   "name": "Alinda",         "num_lcs": 54},
]

# Pole trial grid: 2 trials as specified
POLE_TRIALS = [
    (0.0,       np.pi / 4),    # (0, pi/4)
    (np.pi,    -np.pi / 4),    # (pi, -pi/4)
]

# Inversion parameters
N_FACETS = 120
LAMBDA_SMOOTH = 0.3
MAX_ITER = 80
PERIOD_MIN = 2.0
PERIOD_MAX = 50.0
PERIOD_TIMEOUT = 120  # seconds
MAX_SESSIONS = 20  # cap sessions per asteroid for performance


def log(msg):
    """Print with flush for unbuffered output."""
    print(msg, flush=True)


class AsteroidTimeout(Exception):
    """Raised when processing a single asteroid exceeds the per-asteroid time limit."""
    pass


def _timeout_handler(signum, frame):
    """Signal handler that raises AsteroidTimeout on SIGALRM."""
    raise AsteroidTimeout("Asteroid processing timed out")


def run_single_candidate(candidate, per_asteroid_timeout=480):
    """Run the full inversion pipeline on a single candidate asteroid.

    Args:
        candidate: dict with keys 'number', 'name', 'num_lcs', 'rank'
        per_asteroid_timeout: max seconds per asteroid (default 480 = 8 min)

    Returns:
        dict with all results, or dict with 'error' key on failure
    """
    number = candidate['number']
    name = candidate['name']
    label = f"{number}_{name}"
    log(f"\n{'='*70}")
    log(f"  CANDIDATE {candidate['rank']}/10: {number} {name} ({candidate['num_lcs']} LCs)")
    log(f"{'='*70}")

    result_record = {
        'rank': candidate['rank'],
        'number': number,
        'name': name,
        'expected_lightcurves': candidate['num_lcs'],
    }

    # Set per-asteroid timeout
    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(per_asteroid_timeout)

    try:
        # Step 1: Load ALCDEF data
        t0 = time.time()
        log(f"  Loading ALCDEF data for #{number}...")
        sessions = load_alcdef_asteroid(ALCDEF_ZIP, asteroid_number=number)
        load_time = time.time() - t0
        log(f"  Loaded {len(sessions)} sessions in {load_time:.1f}s")

        if len(sessions) < 2:
            msg = f"Only {len(sessions)} sessions found, need >= 2. Skipping."
            log(f"  SKIP: {msg}")
            result_record['error'] = msg
            result_record['status'] = 'skipped'
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
            return result_record

        # Sort by data count and cap for performance
        sessions.sort(key=lambda s: len(s.jd), reverse=True)
        sessions_used = sessions[:MAX_SESSIONS]
        total_pts = sum(len(s.jd) for s in sessions_used)
        log(f"  Using {len(sessions_used)}/{len(sessions)} sessions, {total_pts} data points")

        result_record['sessions_loaded'] = len(sessions)
        result_record['sessions_used'] = len(sessions_used)
        result_record['total_datapoints'] = total_pts

        # Step 2: Period search
        log(f"  --- Period Search (range {PERIOD_MIN}-{PERIOD_MAX}h, timeout {PERIOD_TIMEOUT}s) ---")
        t0 = time.time()
        times, mags, errs = combine_lightcurve_sessions(sessions_used)
        best_period, _, _ = find_best_period(
            times, mags, errs,
            period_min=PERIOD_MIN,
            period_max=PERIOD_MAX,
            timeout_sec=PERIOD_TIMEOUT
        )
        period_time = time.time() - t0
        log(f"  Best period: {best_period:.4f} h (search took {period_time:.1f}s)")
        result_record['period_search_seconds'] = round(period_time, 1)
        result_record['period_found'] = round(float(best_period), 6)

        # Step 3: Convex inversion with 2 pole trials
        log(f"  --- Convex Inversion ({len(POLE_TRIALS)} pole trials, n_facets={N_FACETS}, "
            f"lambda={LAMBDA_SMOOTH}, max_iter={MAX_ITER}) ---")

        best_result = None
        best_chi2 = np.inf
        trial_results = []

        for ti, (pl, pb) in enumerate(POLE_TRIALS):
            pl_deg = np.degrees(pl)
            pb_deg = np.degrees(pb)
            log(f"    Trial {ti+1}/{len(POLE_TRIALS)}: pole=({pl_deg:.0f} deg, {pb_deg:.0f} deg)")
            t0 = time.time()
            try:
                inv_result = convex_inversion(
                    sessions_used,
                    period_init=best_period,
                    pole_lambda_init=pl,
                    pole_beta_init=pb,
                    n_facets=N_FACETS,
                    lambda_smooth=LAMBDA_SMOOTH,
                    max_iter=MAX_ITER,
                    seed=42,
                )
                inv_time = time.time() - t0
                log(f"      chi2={inv_result.chi2:.1f}  rms={inv_result.residual_rms:.4f}  "
                    f"converged={inv_result.converged}  iters={inv_result.n_iterations}  "
                    f"[{inv_time:.1f}s]")

                trial_info = {
                    'pole_lambda_init_deg': round(pl_deg, 1),
                    'pole_beta_init_deg': round(pb_deg, 1),
                    'chi2': round(float(inv_result.chi2), 2),
                    'rms': round(float(inv_result.residual_rms), 5),
                    'converged': bool(inv_result.converged),
                    'n_iterations': int(inv_result.n_iterations),
                    'period_refined': round(float(inv_result.period), 6),
                    'pole_lambda_deg': round(float(np.degrees(inv_result.pole_lambda)), 2),
                    'pole_beta_deg': round(float(np.degrees(inv_result.pole_beta)), 2),
                    'time_seconds': round(inv_time, 1),
                }
                trial_results.append(trial_info)

                if inv_result.chi2 < best_chi2:
                    best_chi2 = inv_result.chi2
                    best_result = inv_result

            except Exception as e:
                inv_time = time.time() - t0
                log(f"      FAILED: {e} [{inv_time:.1f}s]")
                trial_results.append({
                    'pole_lambda_init_deg': round(pl_deg, 1),
                    'pole_beta_init_deg': round(pb_deg, 1),
                    'error': str(e),
                    'time_seconds': round(inv_time, 1),
                })

        result_record['pole_trials'] = trial_results

        if best_result is None:
            msg = "All inversion trials failed"
            log(f"  ERROR: {msg}")
            result_record['error'] = msg
            result_record['status'] = 'failed'
            signal.alarm(0)
            signal.signal(signal.SIGALRM, old_handler)
            return result_record

        # Step 4: Save best .obj mesh
        # Build filename: candidate_NNNN_Name.obj
        safe_name = name.replace(' ', '_').replace('-', '_')
        obj_filename = f"candidate_{number}_{safe_name}.obj"
        obj_path = os.path.join(RESULTS_DIR, obj_filename)
        save_mesh_obj(best_result.vertices, best_result.faces, obj_path)
        log(f"  Saved mesh: {obj_filename}")
        log(f"    Vertices: {len(best_result.vertices)}, Faces: {len(best_result.faces)}")

        # Step 5: Record final results
        result_record['status'] = 'converged' if best_result.converged else 'completed'
        result_record['period_hours'] = round(float(best_result.period), 6)
        result_record['pole_lambda_deg'] = round(float(np.degrees(best_result.pole_lambda)), 2)
        result_record['pole_beta_deg'] = round(float(np.degrees(best_result.pole_beta)), 2)
        result_record['chi2'] = round(float(best_result.chi2), 2)
        result_record['residual_rms'] = round(float(best_result.residual_rms), 5)
        result_record['converged'] = bool(best_result.converged)
        result_record['n_iterations'] = int(best_result.n_iterations)
        result_record['n_vertices'] = int(len(best_result.vertices))
        result_record['n_faces'] = int(len(best_result.faces))
        result_record['obj_file'] = obj_filename

        log(f"\n  RESULT: P={result_record['period_hours']:.4f}h  "
            f"pole=({result_record['pole_lambda_deg']:.1f}, {result_record['pole_beta_deg']:.1f}) deg  "
            f"chi2={result_record['chi2']:.1f}  rms={result_record['residual_rms']:.4f}  "
            f"converged={result_record['converged']}")

    except AsteroidTimeout:
        log(f"  TIMEOUT: Processing exceeded {per_asteroid_timeout}s limit. Skipping.")
        result_record['error'] = f'Timeout after {per_asteroid_timeout}s'
        result_record['status'] = 'timeout'
    except Exception as e:
        log(f"  ERROR: {e}")
        traceback.print_exc(file=sys.stdout)
        result_record['error'] = str(e)
        result_record['status'] = 'error'
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

    return result_record


def save_results_json(all_results):
    """Save full results to JSON."""
    out_path = os.path.join(RESULTS_DIR, 'candidate_inversion_results.json')
    with open(out_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    log(f"\nFull results saved to: {out_path}")
    return out_path


def save_summary_csv(all_results):
    """Save summary CSV."""
    out_path = os.path.join(RESULTS_DIR, 'candidates_modeled.csv')
    fieldnames = [
        'rank', 'number', 'name', 'status', 'sessions_used', 'total_datapoints',
        'period_hours', 'pole_lambda_deg', 'pole_beta_deg',
        'chi2', 'residual_rms', 'converged', 'n_iterations', 'obj_file'
    ]

    with open(out_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        for rec in all_results:
            row = {k: rec.get(k, '') for k in fieldnames}
            writer.writerow(row)

    log(f"Summary CSV saved to: {out_path}")
    return out_path


def main():
    """Run the full inversion pipeline on all top-10 candidate asteroids and save results."""
    log("=" * 70)
    log("  ASTEROID LIGHTCURVE INVERSION PIPELINE - TOP 10 CANDIDATES")
    log("  Item 020 of research rubric")
    log("=" * 70)
    log(f"  ALCDEF archive: {ALCDEF_ZIP}")
    log(f"  Results dir:    {RESULTS_DIR}")
    log(f"  Random seed:    42")
    log(f"  Period range:   {PERIOD_MIN}-{PERIOD_MAX} h")
    log(f"  Pole trials:    {len(POLE_TRIALS)}")
    log(f"  n_facets:       {N_FACETS}")
    log(f"  lambda_smooth:  {LAMBDA_SMOOTH}")
    log(f"  max_iter:       {MAX_ITER}")
    log(f"  Max sessions:   {MAX_SESSIONS}")
    log("")

    os.makedirs(RESULTS_DIR, exist_ok=True)
    global_start = time.time()

    all_results = []
    success_count = 0
    converged_count = 0

    for candidate in TOP10_CANDIDATES:
        result = run_single_candidate(candidate, per_asteroid_timeout=480)
        all_results.append(result)

        status = result.get('status', 'unknown')
        if status in ('converged', 'completed'):
            success_count += 1
        if result.get('converged', False):
            converged_count += 1

        # Save intermediate results after each asteroid (in case of crash)
        save_results_json(all_results)

        elapsed = time.time() - global_start
        log(f"  [Total elapsed: {elapsed:.0f}s, {success_count}/{len(all_results)} successful]")

    # Final save
    json_path = save_results_json(all_results)
    csv_path = save_summary_csv(all_results)

    total_time = time.time() - global_start

    # Print summary
    log(f"\n{'='*70}")
    log("  PIPELINE SUMMARY")
    log(f"{'='*70}")
    log(f"  Total time:        {total_time:.0f}s ({total_time/60:.1f} min)")
    log(f"  Candidates run:    {len(all_results)}")
    log(f"  Successful models: {success_count}")
    log(f"  Converged models:  {converged_count}")
    log("")

    log(f"  {'Rank':<5} {'Number':<8} {'Name':<20} {'Status':<12} {'Period(h)':<10} "
        f"{'Pole(l,b)':<16} {'RMS':<8} {'Chi2':<10}")
    log(f"  {'-'*5} {'-'*8} {'-'*20} {'-'*12} {'-'*10} {'-'*16} {'-'*8} {'-'*10}")
    for rec in all_results:
        status = rec.get('status', 'unknown')
        period_str = f"{rec['period_hours']:.4f}" if 'period_hours' in rec else '-'
        pole_str = (f"({rec['pole_lambda_deg']:.0f},{rec['pole_beta_deg']:.0f})"
                    if 'pole_lambda_deg' in rec else '-')
        rms_str = f"{rec['residual_rms']:.4f}" if 'residual_rms' in rec else '-'
        chi2_str = f"{rec['chi2']:.1f}" if 'chi2' in rec else '-'
        log(f"  {rec['rank']:<5} {rec['number']:<8} {rec['name']:<20} {status:<12} "
            f"{period_str:<10} {pole_str:<16} {rms_str:<8} {chi2_str:<10}")

    log(f"\n  Output files:")
    log(f"    JSON: {json_path}")
    log(f"    CSV:  {csv_path}")
    for rec in all_results:
        if 'obj_file' in rec:
            log(f"    OBJ:  results/{rec['obj_file']}")

    log(f"\nDone.")


if __name__ == '__main__':
    main()
