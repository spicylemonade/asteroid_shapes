"""
Uncertainty quantification for asteroid shape solutions via jackknife bootstrap.

For the top 10 converged asteroids (by chi2_reduced), performs simplified
bootstrap analysis by re-running convex inversion with different subsets of
lightcurve blocks (jackknife: remove ~20% each time). Computes spin axis
uncertainty as the standard deviation of pole solutions across bootstrap samples.

Classifies each solution as:
  HIGH   (<5 deg uncertainty)
  MEDIUM (5-15 deg uncertainty)
  LOW    (>15 deg uncertainty)

References:
  - Kaasalainen et al. (2001) [Kaasalainen2001a in sources.bib]
  - Durech et al. (2009) [Durech2009 in sources.bib]
"""
import numpy as np
import json
import csv
import os
import sys
import time
import math
import zipfile
import gc

sys.path.insert(0, os.path.dirname(__file__))
from data_ingest import preprocess_asteroid
from convex_inversion import ConvexSolver, grid_search_spin, subsample_blocks

# ============================================================
# Configuration
# ============================================================

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ZIP_PATH = os.path.join(REPO_ROOT, 'ALCDEF_ALL.zip')
GZ_PATH = os.path.join(REPO_ROOT, 'MPCORB.DAT.gz')
RESULTS_DIR = os.path.join(REPO_ROOT, 'results')

PIPELINE_RESULTS_PATH = os.path.join(RESULTS_DIR, 'pipeline_results.json')
SPIN_VECTORS_PATH = os.path.join(RESULTS_DIR, 'spin_vectors.csv')
OUTPUT_PATH = os.path.join(RESULTS_DIR, 'uncertainty_report.csv')

N_TOP_ASTEROIDS = 10
N_BOOTSTRAP = 3
JACKKNIFE_DROP_FRACTION = 0.20

# Coarse grid search settings for bootstrap resamples (fast)
BOOTSTRAP_LAM_STEP = 90
BOOTSTRAP_BETA_STEP = 45
BOOTSTRAP_MAX_ITER = 10

np.random.seed(42)


# ============================================================
# Angular separation on the sphere
# ============================================================

def angular_separation_deg(lam1, beta1, lam2, beta2):
    """Compute angular separation between two sky positions (all in degrees).

    Uses the Vincenty formula for numerical stability.
    """
    lam1_r = math.radians(lam1)
    beta1_r = math.radians(beta1)
    lam2_r = math.radians(lam2)
    beta2_r = math.radians(beta2)

    dlam = lam2_r - lam1_r

    num = math.sqrt(
        (math.cos(beta2_r) * math.sin(dlam)) ** 2 +
        (math.cos(beta1_r) * math.sin(beta2_r) -
         math.sin(beta1_r) * math.cos(beta2_r) * math.cos(dlam)) ** 2
    )
    den = (math.sin(beta1_r) * math.sin(beta2_r) +
           math.cos(beta1_r) * math.cos(beta2_r) * math.cos(dlam))

    sep = math.atan2(num, den)
    return math.degrees(abs(sep))


def pole_scatter_deg(solutions):
    """Compute the standard deviation of pole directions across solutions.

    Converts pole directions to Cartesian unit vectors, computes the mean
    direction, and returns the RMS angular deviation from the mean.

    Parameters
    ----------
    solutions : list of dict
        Each dict must have 'spin_lam_deg' and 'spin_beta_deg'.

    Returns
    -------
    float
        Standard deviation of pole direction in degrees.
    """
    if len(solutions) < 2:
        return 0.0

    # Convert to Cartesian
    vecs = []
    for sol in solutions:
        lam = math.radians(sol['spin_lam_deg'])
        beta = math.radians(sol['spin_beta_deg'])
        x = math.cos(beta) * math.cos(lam)
        y = math.cos(beta) * math.sin(lam)
        z = math.sin(beta)
        vecs.append(np.array([x, y, z]))

    vecs = np.array(vecs)
    mean_vec = np.mean(vecs, axis=0)
    mean_norm = np.linalg.norm(mean_vec)

    if mean_norm < 1e-10:
        # Poles are scattered all over the sphere
        return 180.0

    mean_vec /= mean_norm

    # Compute angular deviations from mean
    deviations = []
    for v in vecs:
        cos_angle = np.clip(np.dot(v, mean_vec), -1.0, 1.0)
        angle_deg = math.degrees(math.acos(cos_angle))
        deviations.append(angle_deg)

    return float(np.std(deviations))


def classify_confidence(uncertainty_deg):
    """Classify solution confidence based on pole uncertainty."""
    if uncertainty_deg < 5.0:
        return "HIGH"
    elif uncertainty_deg <= 15.0:
        return "MEDIUM"
    else:
        return "LOW"


# ============================================================
# Jackknife block selection
# ============================================================

def jackknife_block_subsets(blocks, n_samples, drop_fraction=0.20):
    """Generate jackknife subsets of lightcurve blocks.

    Each subset drops approximately drop_fraction of the blocks.

    Parameters
    ----------
    blocks : list
        Full list of lightcurve blocks.
    n_samples : int
        Number of jackknife subsets to generate.
    drop_fraction : float
        Fraction of blocks to drop in each subset.

    Returns
    -------
    list of list
        Each inner list is a subset of blocks.
    """
    n_blocks = len(blocks)
    n_drop = max(1, int(round(n_blocks * drop_fraction)))

    subsets = []
    for i in range(n_samples):
        # Deterministic but different for each sample
        rng = np.random.RandomState(42 + i)
        drop_indices = set(rng.choice(n_blocks, size=min(n_drop, n_blocks - 1), replace=False))
        subset = [b for j, b in enumerate(blocks) if j not in drop_indices]
        if len(subset) >= 2:
            subsets.append(subset)

    return subsets


# ============================================================
# Main uncertainty analysis
# ============================================================

def load_pipeline_results():
    """Load and return pipeline results sorted by chi2_reduced."""
    with open(PIPELINE_RESULTS_PATH) as f:
        results = json.load(f)

    # Filter to converged asteroids and sort by chi2_reduced
    converged = []
    for key, val in results.items():
        if val.get('status') in ('converged', 'converged_convex_only'):
            converged.append(val)

    converged.sort(key=lambda x: x['chi2_reduced'])
    return converged


def load_spin_vectors():
    """Load spin vectors CSV into a dict keyed by asteroid_number."""
    spin_dict = {}
    with open(SPIN_VECTORS_PATH) as f:
        reader = csv.DictReader(f)
        for row in reader:
            ast_num = int(row['asteroid_number'])
            spin_dict[ast_num] = {
                'period_hours': float(row['period_hours']),
                'spin_lam_deg': float(row['spin_lam_deg']),
                'spin_beta_deg': float(row['spin_beta_deg']),
                'chi2_reduced': float(row['chi2_reduced']),
            }
    return spin_dict


def run_bootstrap_for_asteroid(ast_num, period_hours, original_lam, original_beta):
    """Run jackknife bootstrap analysis for a single asteroid.

    Parameters
    ----------
    ast_num : int
        Asteroid number.
    period_hours : float
        Best-fit period from pipeline.
    original_lam : float
        Original spin longitude (deg).
    original_beta : float
        Original spin latitude (deg).

    Returns
    -------
    dict with bootstrap results or None on failure.
    """
    print(f"  Loading ALCDEF data for asteroid {ast_num}...")
    try:
        data = preprocess_asteroid(ZIP_PATH, GZ_PATH, ast_num)
    except Exception as e:
        print(f"    ERROR loading data: {e}")
        return None

    all_blocks = data['blocks']
    print(f"    Found {len(all_blocks)} lightcurve blocks")

    if len(all_blocks) < 3:
        print(f"    Insufficient blocks for jackknife ({len(all_blocks)} < 3)")
        # Return a result with the original solution and high uncertainty marker
        return {
            'solutions': [{'spin_lam_deg': original_lam, 'spin_beta_deg': original_beta}],
            'pole_uncertainty_deg': 999.0,
            'n_bootstrap_samples': 0,
        }

    # Filter to blocks with enough data
    good_blocks = [b for b in all_blocks if len(b['data']) >= 5]
    if len(good_blocks) < 3:
        good_blocks = all_blocks[:max(3, len(all_blocks))]

    # Limit to manageable number of blocks for speed
    if len(good_blocks) > 10:
        good_blocks = good_blocks[:10]

    print(f"    Using {len(good_blocks)} blocks for jackknife analysis")

    # Generate jackknife subsets
    subsets = jackknife_block_subsets(good_blocks, N_BOOTSTRAP, JACKKNIFE_DROP_FRACTION)
    print(f"    Generated {len(subsets)} jackknife subsets")

    # Include the original solution
    solutions = [{'spin_lam_deg': original_lam, 'spin_beta_deg': original_beta}]

    # Create a single solver instance (reuse across bootstrap samples)
    solver = ConvexSolver(n_subdivisions=1, c_ls=0.5, c_l=0.1)

    for i, subset in enumerate(subsets):
        print(f"    Bootstrap sample {i+1}/{len(subsets)}: "
              f"{len(subset)} blocks, "
              f"{sum(len(b['data']) for b in subset)} raw points...")
        t0 = time.time()

        try:
            # Subsample for speed
            sub_blocks = subsample_blocks(subset, max_points_per_block=10)
            inv_pts = sum(len(b['data']) for b in sub_blocks)

            if inv_pts < 10:
                print(f"      Skipping (only {inv_pts} points after subsampling)")
                continue

            best, _ = grid_search_spin(
                solver, sub_blocks, period_hours,
                lam_step=BOOTSTRAP_LAM_STEP,
                beta_step=BOOTSTRAP_BETA_STEP,
                max_iter=BOOTSTRAP_MAX_ITER
            )

            elapsed = time.time() - t0

            if best is not None:
                solutions.append({
                    'spin_lam_deg': best['spin_lam_deg'],
                    'spin_beta_deg': best['spin_beta_deg'],
                })
                print(f"      Result: lam={best['spin_lam_deg']:.1f}, "
                      f"beta={best['spin_beta_deg']:.1f}, "
                      f"chi2_red={best['chi2_reduced']:.2f} "
                      f"({elapsed:.1f}s)")
            else:
                print(f"      Inversion failed ({elapsed:.1f}s)")

        except Exception as e:
            print(f"      ERROR: {e}")
            continue

    # Compute pole scatter
    uncertainty = pole_scatter_deg(solutions)
    n_samples = len(solutions) - 1  # Excluding the original

    print(f"    Pole uncertainty: {uncertainty:.2f} deg "
          f"({len(solutions)} solutions including original)")

    return {
        'solutions': solutions,
        'pole_uncertainty_deg': uncertainty,
        'n_bootstrap_samples': n_samples,
    }


def main():
    """Run uncertainty quantification on top converged asteroids."""
    print("=" * 70)
    print("UNCERTAINTY QUANTIFICATION FOR ASTEROID SHAPE SOLUTIONS")
    print("=" * 70)
    print(f"Pipeline results: {PIPELINE_RESULTS_PATH}")
    print(f"Spin vectors:     {SPIN_VECTORS_PATH}")
    print(f"Output:           {OUTPUT_PATH}")
    print(f"N bootstrap:      {N_BOOTSTRAP}")
    print(f"Jackknife drop:   {JACKKNIFE_DROP_FRACTION*100:.0f}%")
    print(f"Grid steps:       lam={BOOTSTRAP_LAM_STEP}, beta={BOOTSTRAP_BETA_STEP}")
    print()

    # Load pipeline results
    converged = load_pipeline_results()
    print(f"Found {len(converged)} converged asteroids")

    # Select top N by chi2_reduced
    top_asteroids = converged[:N_TOP_ASTEROIDS]
    print(f"Analyzing top {len(top_asteroids)} by chi2_reduced:")
    for ast in top_asteroids:
        print(f"  #{ast['asteroid_number']}: chi2_red={ast['chi2_reduced']:.4f}, "
              f"P={ast['period_hours']:.4f}h")

    # Load spin vectors for cross-reference
    spin_dict = load_spin_vectors()

    # Run bootstrap analysis
    print("\n" + "=" * 70)
    print("RUNNING JACKKNIFE BOOTSTRAP ANALYSIS")
    print("=" * 70)

    report_rows = []
    t_total = time.time()

    for i, ast_info in enumerate(top_asteroids):
        ast_num = ast_info['asteroid_number']
        period_hours = ast_info['period_hours']
        spin_lam = ast_info['spin_lam_deg']
        spin_beta = ast_info['spin_beta_deg']
        chi2_red = ast_info['chi2_reduced']

        print(f"\n[{i+1}/{len(top_asteroids)}] Asteroid {ast_num} "
              f"(chi2_red={chi2_red:.4f}, P={period_hours:.4f}h)")
        print("-" * 60)

        result = run_bootstrap_for_asteroid(
            ast_num, period_hours, spin_lam, spin_beta
        )

        if result is not None:
            uncertainty = result['pole_uncertainty_deg']
            confidence = classify_confidence(uncertainty)
            n_samples = result['n_bootstrap_samples']

            # Compute mean pole from all solutions for the report
            if len(result['solutions']) > 1:
                mean_lam = np.mean([s['spin_lam_deg'] for s in result['solutions']])
                mean_beta = np.mean([s['spin_beta_deg'] for s in result['solutions']])
            else:
                mean_lam = spin_lam
                mean_beta = spin_beta

            row = {
                'asteroid_number': ast_num,
                'period_hours': round(period_hours, 6),
                'spin_lam_deg': round(mean_lam, 2),
                'spin_beta_deg': round(mean_beta, 2),
                'pole_uncertainty_deg': round(uncertainty, 2),
                'confidence_level': confidence,
                'n_bootstrap_samples': n_samples,
                'chi2_reduced': round(chi2_red, 4),
            }
        else:
            row = {
                'asteroid_number': ast_num,
                'period_hours': round(period_hours, 6),
                'spin_lam_deg': round(spin_lam, 2),
                'spin_beta_deg': round(spin_beta, 2),
                'pole_uncertainty_deg': -1.0,
                'confidence_level': 'FAILED',
                'n_bootstrap_samples': 0,
                'chi2_reduced': round(chi2_red, 4),
            }

        report_rows.append(row)
        print(f"  => Confidence: {row['confidence_level']}, "
              f"Uncertainty: {row['pole_uncertainty_deg']} deg")

        gc.collect()

    # Save report
    print("\n" + "=" * 70)
    print("SAVING RESULTS")
    print("=" * 70)

    fieldnames = [
        'asteroid_number', 'period_hours', 'spin_lam_deg', 'spin_beta_deg',
        'pole_uncertainty_deg', 'confidence_level', 'n_bootstrap_samples',
        'chi2_reduced'
    ]

    os.makedirs(RESULTS_DIR, exist_ok=True)
    with open(OUTPUT_PATH, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)

    print(f"Report saved to {OUTPUT_PATH}")

    elapsed_total = time.time() - t_total
    print(f"\nTotal elapsed time: {elapsed_total:.1f}s")

    # Summary
    print("\n" + "=" * 70)
    print("UNCERTAINTY SUMMARY")
    print("=" * 70)
    print(f"{'Asteroid':>10s} {'Period(h)':>10s} {'Lam':>7s} {'Beta':>7s} "
          f"{'Uncert':>8s} {'Confidence':>12s} {'N_boot':>7s} {'Chi2_red':>10s}")
    print("-" * 80)
    for row in report_rows:
        print(f"{row['asteroid_number']:>10d} {row['period_hours']:>10.4f} "
              f"{row['spin_lam_deg']:>7.1f} {row['spin_beta_deg']:>7.1f} "
              f"{row['pole_uncertainty_deg']:>8.2f} {row['confidence_level']:>12s} "
              f"{row['n_bootstrap_samples']:>7d} {row['chi2_reduced']:>10.4f}")

    n_high = sum(1 for r in report_rows if r['confidence_level'] == 'HIGH')
    n_med = sum(1 for r in report_rows if r['confidence_level'] == 'MEDIUM')
    n_low = sum(1 for r in report_rows if r['confidence_level'] == 'LOW')
    n_fail = sum(1 for r in report_rows if r['confidence_level'] == 'FAILED')

    print(f"\nHIGH confidence:   {n_high}")
    print(f"MEDIUM confidence: {n_med}")
    print(f"LOW confidence:    {n_low}")
    if n_fail:
        print(f"FAILED:            {n_fail}")

    print("\nDone.")


if __name__ == '__main__':
    main()
