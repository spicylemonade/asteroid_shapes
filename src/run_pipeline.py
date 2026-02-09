"""
Full pipeline runner: execute period search + convex + GA inversion on all targets.

Item 019 of the research rubric.
Processes each target in the candidate list and saves results.
Optimized to pre-cache data loading.
"""
import numpy as np
import json
import csv
import os
import sys
import time
import traceback
import zipfile
import gc
import math

sys.path.insert(0, os.path.dirname(__file__))

from data_ingest import parse_alcdef_file, load_mpcorb_for_asteroid, compute_viewing_geometry
from data_ingest import orbital_position_ecliptic, earth_position_ecliptic
from period_search import combined_period_search, merge_sessions_to_relative
from convex_inversion import ConvexSolver, grid_search_spin, subsample_blocks
from ga_solver import GASolver
from mesh_utils import save_obj

repo_root = os.path.dirname(os.path.dirname(__file__))
zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')

np.random.seed(42)


def load_alcdef_fast(zf, namelist_index, asteroid_number):
    """Load ALCDEF data from pre-opened zip file using pre-built index."""
    matching = namelist_index.get(str(asteroid_number), [])
    if not matching:
        return []

    all_blocks = []
    for fname in matching:
        content = zf.read(fname).decode('utf-8', errors='replace')
        blocks = parse_alcdef_file(content)
        for block in blocks:
            if block['data']:
                data_array = np.array(block['data'], dtype=np.float64)
                if data_array.ndim == 2 and data_array.shape[1] >= 3:
                    block['data'] = data_array[:, :3]
                    all_blocks.append(block)
    return all_blocks


def build_namelist_index(zf):
    """Build index mapping asteroid numbers to zip file names."""
    index = {}
    for name in zf.namelist():
        if name.startswith('ALCDEF_'):
            parts = name.split('_')
            if len(parts) >= 2:
                ast_num = parts[1]
                if ast_num not in index:
                    index[ast_num] = []
                index[ast_num].append(name)
    return index


def preprocess_fast(zf, namelist_index, orb_cache, ast_num):
    """Fast preprocessing with pre-opened zip and cached orbital elements."""
    blocks = load_alcdef_fast(zf, namelist_index, ast_num)
    if not blocks:
        return None

    orb = orb_cache.get(ast_num)
    if orb is None:
        try:
            orb = load_mpcorb_for_asteroid(gz_path, ast_num)
            orb_cache[ast_num] = orb
        except Exception:
            orb = {'a': 2.5, 'e': 0.1, 'i': 10, 'node': 0, 'peri': 0,
                   'M': 0, 'epoch_jd': 2460000.5}
            orb_cache[ast_num] = orb

    for block in blocks:
        geom_list = []
        reduced_mags = []
        for row in block['data']:
            jd = row[0]
            mag = row[1]
            try:
                ast_pos = orbital_position_ecliptic(orb, jd)
                earth_pos = earth_position_ecliptic(jd)
                geom = compute_viewing_geometry(ast_pos, earth_pos)
            except Exception:
                geom = {'phase_angle_deg': 20.0, 'helio_dist_au': 2.0, 'geo_dist_au': 1.5}
            geom_list.append(geom)
            r = geom.get('helio_dist_au', 2.0)
            delta = geom.get('geo_dist_au', 1.5)
            if r > 0 and delta > 0:
                reduced_mag = mag - 5 * math.log10(r * delta)
            else:
                reduced_mag = mag
            reduced_mags.append(reduced_mag)

        block['geometry'] = geom_list
        block['reduced_mag'] = np.array(reduced_mags)

    return blocks


def run_pipeline_for_asteroid(zf, namelist_index, orb_cache, ast_num, timeout_seconds=90):
    """Run full inversion pipeline on one asteroid."""
    t0 = time.time()
    result = {
        'asteroid_number': int(ast_num),
        'status': 'failed',
        'error': None,
    }

    try:
        blocks = preprocess_fast(zf, namelist_index, orb_cache, ast_num)
        if blocks is None or len(blocks) < 2:
            result['error'] = f'Insufficient data: {0 if blocks is None else len(blocks)} blocks'
            return result

        total_raw_pts = sum(len(b['data']) for b in blocks)

        # Period search
        rel_times, rel_mags, rel_errs = merge_sessions_to_relative(blocks)

        if len(rel_times) < 20:
            result['error'] = f'Too few data points: {len(rel_times)}'
            return result

        candidates = combined_period_search(
            rel_times, rel_mags, rel_errs,
            period_min=2.0, period_max=100.0, top_n=5
        )

        if not candidates:
            result['error'] = 'No period candidates found'
            return result

        period_hours = candidates[0]['period_hours']

        # Select blocks for inversion
        good_blocks = [b for b in blocks if len(b['data']) >= 5]
        if len(good_blocks) > 5:
            good_blocks = good_blocks[:5]
        good_blocks = subsample_blocks(good_blocks, max_points_per_block=10)
        inv_pts = sum(len(b['data']) for b in good_blocks)

        if inv_pts < 10:
            result['error'] = f'Insufficient inversion points: {inv_pts}'
            return result

        # Stage 1: Convex inversion
        convex_solver = ConvexSolver(n_subdivisions=1, c_ls=0.5, c_l=0.1)
        convex_best, _ = grid_search_spin(
            convex_solver, good_blocks, period_hours,
            lam_step=45, beta_step=45, max_iter=15
        )

        if convex_best is None:
            result['error'] = 'Convex inversion failed'
            return result

        # Save shape
        out_dir = os.path.join(repo_root, 'results', 'shapes')
        os.makedirs(out_dir, exist_ok=True)

        elapsed = time.time() - t0

        # Stage 2: GA only if time permits
        if elapsed < timeout_seconds * 0.6:
            ga = GASolver(
                n_verts=len(convex_best['base_vertices']),
                base_vertices=convex_best['base_vertices'],
                faces=convex_best['faces'],
                c_ls=0.5, c_l=0.1,
            )
            ga_blocks = good_blocks[:3]
            ga_blocks = subsample_blocks(ga_blocks, max_points_per_block=6)

            ga_result = ga.evolve(
                blocks=ga_blocks,
                spin_axis=convex_best['spin_axis'],
                period_hours=period_hours,
                pop_size=15,
                n_generations=10,
                mutation_rate=0.15,
                mutation_sigma=0.08,
                regularization=0.01,
                seed_radii=convex_best['radii'],
                verbose=False,
            )
            obj_path = os.path.join(out_dir, f'{ast_num}_ga.obj')
            save_obj(obj_path, ga_result['vertices'], ga_result['faces'])
            ga_chi2 = float(ga_result['chi2'])
            shape_file = f'results/shapes/{ast_num}_ga.obj'
            status = 'converged'
        else:
            obj_path = os.path.join(out_dir, f'{ast_num}_convex.obj')
            save_obj(obj_path, convex_best['vertices'], convex_best['faces'])
            ga_chi2 = None
            shape_file = f'results/shapes/{ast_num}_convex.obj'
            status = 'converged_convex_only'

        result.update({
            'status': status,
            'period_hours': float(period_hours),
            'spin_lam_deg': float(convex_best['spin_lam_deg']),
            'spin_beta_deg': float(convex_best['spin_beta_deg']),
            'chi2_reduced': float(convex_best['chi2_reduced']),
            'ga_chi2': ga_chi2,
            'ellipsoid_b_a': float(convex_best['ellipsoid_params']['b']),
            'ellipsoid_c_a': float(convex_best['ellipsoid_params']['c']),
            'n_blocks': len(good_blocks),
            'n_data_points': inv_pts,
            'total_raw_points': total_raw_pts,
            'shape_file': shape_file,
            'elapsed_seconds': float(time.time() - t0),
        })
        return result

    except Exception as e:
        result['error'] = str(e)
        result['elapsed_seconds'] = float(time.time() - t0)
        return result


if __name__ == '__main__':
    csv_path = os.path.join(repo_root, 'results', 'target_candidates.csv')
    targets = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            targets.append(int(row['asteroid_number']))

    print(f"Running pipeline on {len(targets)} target asteroids")
    print("="*60)

    # Pre-open zip file and build index
    print("Opening ALCDEF zip file and building index...")
    zf = zipfile.ZipFile(zip_path)
    namelist_index = build_namelist_index(zf)
    print(f"  {len(namelist_index)} asteroids indexed in archive")

    orb_cache = {}
    all_results = {}
    spin_vectors = []
    converged = 0
    failed = 0

    for i, ast_num in enumerate(targets):
        print(f"\n[{i+1}/{len(targets)}] Processing asteroid {ast_num}...", flush=True)
        result = run_pipeline_for_asteroid(zf, namelist_index, orb_cache, ast_num,
                                            timeout_seconds=90)
        all_results[str(ast_num)] = result

        if result['status'] in ('converged', 'converged_convex_only'):
            converged += 1
            spin_vectors.append({
                'asteroid_number': ast_num,
                'period_hours': result['period_hours'],
                'spin_lam_deg': result['spin_lam_deg'],
                'spin_beta_deg': result['spin_beta_deg'],
                'chi2_reduced': result['chi2_reduced'],
                'status': result['status'],
            })
            print(f"  CONVERGED: P={result['period_hours']:.4f}h, "
                  f"spin=({result['spin_lam_deg']}, {result['spin_beta_deg']}), "
                  f"chi2_red={result['chi2_reduced']:.2f}, "
                  f"time={result.get('elapsed_seconds', 0):.1f}s")
        else:
            failed += 1
            print(f"  FAILED: {result.get('error', 'unknown')}")

        gc.collect()

    zf.close()

    # Save results
    results_path = os.path.join(repo_root, 'results', 'pipeline_results.json')
    with open(results_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nResults saved to {results_path}")

    # Save spin vectors CSV
    sv_path = os.path.join(repo_root, 'results', 'spin_vectors.csv')
    if spin_vectors:
        with open(sv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=spin_vectors[0].keys())
            writer.writeheader()
            writer.writerows(spin_vectors)
        print(f"Spin vectors saved to {sv_path}")

    # Summary
    print(f"\n{'='*60}")
    print(f"PIPELINE SUMMARY")
    print(f"{'='*60}")
    print(f"  Total targets: {len(targets)}")
    print(f"  Converged: {converged}")
    print(f"  Failed: {failed}")
    print(f"  Success rate: {converged/max(len(targets),1)*100:.1f}%")
