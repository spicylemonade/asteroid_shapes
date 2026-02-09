"""
Sparse data inversion module for asteroid lightcurve inversion.
Handles survey-like sparse photometry (Gaia/ZTF/Pan-STARRS style)
with <100 data points across multiple apparitions.

References:
  - Durech et al. (2010) [Durech2009, Durech2010 in sources.bib]
  - Cellino et al. (2015) [Cellino2015 in sources.bib]
  - Viikinkoski et al. (2015) [Viikinkoski2015 in sources.bib]
"""
import numpy as np
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from convex_inversion import ConvexSolver, grid_search_spin, subsample_blocks


def create_sparse_subset(blocks, max_points=80, min_apparitions=3, seed=42):
    """Create a sparse data subset from dense lightcurves.

    Simulates survey-like sparse sampling by taking 1-2 points
    per session spread across multiple apparitions.

    Parameters
    ----------
    blocks : list of preprocessed LC blocks
    max_points : int, max total points
    min_apparitions : int, minimum different observing epochs (>120 day gaps)
    seed : int

    Returns
    -------
    sparse_blocks : list of single-point blocks
    """
    rng = np.random.RandomState(seed)

    # Group blocks by apparition (sessions >120 days apart)
    apparitions = []
    current_app = []
    last_jd = None

    sorted_blocks = sorted(blocks, key=lambda b: b['data'][0, 0])

    for block in sorted_blocks:
        jd = block['data'][0, 0]
        if last_jd is not None and (jd - last_jd) > 120:
            if current_app:
                apparitions.append(current_app)
            current_app = []
        current_app.append(block)
        last_jd = jd

    if current_app:
        apparitions.append(current_app)

    if len(apparitions) < min_apparitions:
        # Relax: treat each block as its own apparition
        apparitions = [[b] for b in sorted_blocks]

    # Sample points per apparition
    points_per_app = max(1, max_points // len(apparitions))
    sparse_blocks = []

    for app in apparitions:
        # Sample a few points from each session in this apparition
        for block in app:
            n = len(block['data'])
            n_sample = min(points_per_app, n, 3)
            if n_sample < 1:
                continue
            indices = rng.choice(n, n_sample, replace=False)
            indices.sort()

            new_block = {
                'meta': block['meta'],
                'data': block['data'][indices],
                'geometry': [block['geometry'][i] for i in indices],
            }
            if 'reduced_mag' in block:
                new_block['reduced_mag'] = block['reduced_mag'][indices]
            sparse_blocks.append(new_block)

            if sum(len(b['data']) for b in sparse_blocks) >= max_points:
                break
        if sum(len(b['data']) for b in sparse_blocks) >= max_points:
            break

    total_pts = sum(len(b['data']) for b in sparse_blocks)
    n_app = len(apparitions)
    print(f"Sparse subset: {total_pts} points across {n_app} apparitions from {len(sparse_blocks)} blocks")

    return sparse_blocks


def sparse_inversion(blocks, period_hours, lam_step=60, beta_step=30,
                      n_subdivisions=1, c_ls=0.5, c_l=0.1):
    """Run sparse data inversion.

    Uses coarse grid search with convex model on sparse data points.

    Parameters
    ----------
    blocks : list of sparse LC blocks (absolute photometry)
    period_hours : float, rotation period
    lam_step, beta_step : spin axis grid resolution

    Returns
    -------
    best solution dict, all solutions
    """
    solver = ConvexSolver(n_subdivisions=n_subdivisions, c_ls=c_ls, c_l=c_l)

    # Run grid search
    best, all_solutions = grid_search_spin(solver, blocks, period_hours,
                                            lam_step=lam_step, beta_step=beta_step,
                                            max_iter=20)

    return best, all_solutions


def data_fusion_inversion(dense_blocks, sparse_blocks, period_hours,
                           dense_weight=1.0, sparse_weight=0.5,
                           lam_step=60, beta_step=30):
    """Hybrid dense + sparse data fusion inversion.

    Dense lightcurves contribute shape detail.
    Sparse absolute photometry constrains pole orientation.

    Parameters
    ----------
    dense_blocks : list of dense LC blocks
    sparse_blocks : list of sparse LC blocks
    period_hours : float
    dense_weight : float, weight for dense data
    sparse_weight : float, weight for sparse data

    Returns
    -------
    best solution dict
    """
    solver = ConvexSolver(n_subdivisions=1, c_ls=0.5, c_l=0.1)

    # Combine dense and sparse blocks with appropriate weighting
    # Scale error bars to adjust relative contribution
    combined = []

    for block in dense_blocks:
        new_block = {
            'meta': block['meta'],
            'data': block['data'].copy(),
            'geometry': block['geometry'],
        }
        # Increase weight of dense data by reducing errors
        new_block['data'][:, 2] /= dense_weight
        combined.append(new_block)

    for block in sparse_blocks:
        new_block = {
            'meta': block['meta'],
            'data': block['data'].copy(),
            'geometry': block['geometry'],
        }
        new_block['data'][:, 2] /= sparse_weight
        combined.append(new_block)

    combined = subsample_blocks(combined, max_points_per_block=15)

    best, all_solutions = grid_search_spin(solver, combined, period_hours,
                                            lam_step=lam_step, beta_step=beta_step,
                                            max_iter=20)

    return best, all_solutions


if __name__ == '__main__':
    import json
    from data_ingest import preprocess_asteroid
    from mesh_utils import save_obj

    repo_root = os.path.dirname(os.path.dirname(__file__))
    zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
    gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')

    # Test on 1036 Ganymed
    ast_num = 1036
    known_period = 10.314

    print(f"Loading data for asteroid {ast_num}...")
    data = preprocess_asteroid(zip_path, gz_path, ast_num)
    blocks = data['blocks']

    # Test sparse inversion
    print("\n=== Sparse Inversion ===")
    sparse = create_sparse_subset(blocks, max_points=60, min_apparitions=3)
    sparse_best, _ = sparse_inversion(sparse, known_period, lam_step=90, beta_step=45)

    if sparse_best:
        print(f"Sparse result: lam={sparse_best['spin_lam_deg']}, "
              f"beta={sparse_best['spin_beta_deg']}, "
              f"chi2_red={sparse_best['chi2_reduced']:.2f}")

    # Test dense inversion
    print("\n=== Dense Inversion ===")
    dense_sub = subsample_blocks([b for b in blocks if len(b['data']) >= 10][:5],
                                  max_points_per_block=15)
    solver = ConvexSolver(n_subdivisions=1, c_ls=0.5, c_l=0.1)
    dense_best, _ = grid_search_spin(solver, dense_sub, known_period,
                                      lam_step=90, beta_step=45, max_iter=20)

    if dense_best:
        print(f"Dense result: lam={dense_best['spin_lam_deg']}, "
              f"beta={dense_best['spin_beta_deg']}, "
              f"chi2_red={dense_best['chi2_reduced']:.2f}")

    # Test fusion
    print("\n=== Fusion Inversion ===")
    fusion_best, _ = data_fusion_inversion(dense_sub, sparse, known_period,
                                            lam_step=90, beta_step=45)

    if fusion_best:
        print(f"Fusion result: lam={fusion_best['spin_lam_deg']}, "
              f"beta={fusion_best['spin_beta_deg']}, "
              f"chi2_red={fusion_best['chi2_reduced']:.2f}")

    # Compare pole accuracies
    # Known Ganymed pole: approximately lam=18, beta=-47 (Durech et al.)
    known_lam, known_beta = 18, -47
    results = {}
    for name, sol in [('sparse', sparse_best), ('dense', dense_best), ('fusion', fusion_best)]:
        if sol:
            # Angular distance
            dlam = sol['spin_lam_deg'] - known_lam
            dbeta = sol['spin_beta_deg'] - known_beta
            ang_dist = np.sqrt(dlam**2 + dbeta**2)
            results[name] = {
                'spin_lam': float(sol['spin_lam_deg']),
                'spin_beta': float(sol['spin_beta_deg']),
                'chi2_reduced': float(sol['chi2_reduced']),
                'pole_error_deg': float(ang_dist),
            }
            print(f"\n{name}: pole error = {ang_dist:.1f} deg")

    # Save results
    with open(os.path.join(repo_root, 'results', 'sparse_fusion_comparison.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print("\nResults saved to results/sparse_fusion_comparison.json")
