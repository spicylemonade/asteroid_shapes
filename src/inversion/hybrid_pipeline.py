"""
Hybrid Multi-Method Inversion Pipeline

Combines convex inversion, genetic non-convex refinement, and sparse
data fusion into a unified pipeline.
"""

import json
import os
import numpy as np

from .period_search import find_period
from .convex_solver import convex_inversion
from .genetic_solver import genetic_inversion
from .sparse_solver import sparse_inversion
from ..shapes.convex_model import ConvexShapeModel

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


def load_config(config_path=None):
    """Load pipeline configuration from YAML or return defaults."""
    defaults = {
        'period_search': {
            'min_hours': 2.0,
            'max_hours': 100.0,
            'n_frequencies': 50000,
        },
        'convex_inversion': {
            'lmax': 6,
            'subdivisions': 3,
            'max_iter': 150,
            'lambda_smooth': 0.001,
            'n_pole_grid': 12,
        },
        'genetic_refinement': {
            'enabled': True,
            'subdivisions': 2,
            'population_size': 50,
            'generations': 80,
            'mutation_rate': 0.15,
            'lambda_smooth': 0.01,
        },
        'sparse_fusion': {
            'enabled': False,
            'weight': 0.1,
        },
        'output': {
            'save_obj': True,
            'save_spin': True,
        },
    }

    if config_path and os.path.exists(config_path):
        try:
            import yaml
            with open(config_path) as f:
                user_config = yaml.safe_load(f)
            if user_config:
                _deep_update(defaults, user_config)
        except ImportError:
            pass

    return defaults


def _deep_update(base, update):
    for k, v in update.items():
        if isinstance(v, dict) and k in base and isinstance(base[k], dict):
            _deep_update(base[k], v)
        else:
            base[k] = v


def run_hybrid_pipeline(lightcurves, sparse_data=None, known_period=None,
                        config=None, output_dir=None, asteroid_id=None,
                        verbose=True, seed=42):
    """Run the full hybrid inversion pipeline.

    Parameters
    ----------
    lightcurves : list of dicts
        Dense lightcurve data with 'times','mags','errors','sun_dirs','obs_dirs','phase_angles'
    sparse_data : list of dicts or None
        Sparse photometric points
    known_period : float or None
        Known rotation period in hours (if available)
    config : dict or None
        Pipeline configuration
    output_dir : str or None
        Directory for output files
    asteroid_id : str or None
        Asteroid identifier for filenames
    verbose : bool
    seed : int

    Returns
    -------
    result : dict with 'convex_result', 'genetic_result', 'final_model',
             'period_hours', 'pole_lambda_deg', 'pole_beta_deg', 'chi2_reduced'
    """
    if config is None:
        config = load_config()

    result = {
        'asteroid_id': asteroid_id,
        'converged': False,
    }

    # ===== Step 1: Period Search =====
    if verbose:
        print(f"\n{'='*60}")
        print(f"HYBRID PIPELINE: {asteroid_id or 'Unknown'}")
        print(f"{'='*60}")

    if known_period is not None:
        period_days = known_period / 24.0
        if verbose:
            print(f"\nUsing known period: {known_period:.4f} hours")
    else:
        if verbose:
            print(f"\nStep 1: Period Search")

        # Combine all lightcurve data for period search
        all_times = []
        all_mags = []
        for lc in lightcurves:
            all_times.extend(lc['times'])
            all_mags.extend(lc['mags'])

        all_times = np.array(all_times)
        all_mags = np.array(all_mags)

        # Convert to hours for period search
        pcfg = config['period_search']
        candidates = find_period(
            all_times * 24.0, all_mags,
            period_min=pcfg['min_hours'],
            period_max=pcfg['max_hours'],
            n_top=5,
        )

        if candidates:
            period_days = candidates[0]['period'] / 24.0
            if verbose:
                print(f"  Best period: {candidates[0]['period']:.4f} hours "
                      f"(score={candidates[0]['score']:.3f})")
        else:
            period_days = 6.0 / 24.0
            if verbose:
                print("  No period found, using default 6h")

    result['period_hours'] = period_days * 24.0
    result['period_days'] = period_days

    # ===== Step 2: Convex Inversion =====
    if verbose:
        print(f"\nStep 2: Convex Inversion")

    ccfg = config['convex_inversion']
    convex_result = convex_inversion(
        lightcurves, period_days,
        lmax=ccfg['lmax'],
        subdivisions=ccfg['subdivisions'],
        max_iter=ccfg['max_iter'],
        lambda_smooth=ccfg['lambda_smooth'],
        n_pole_grid=ccfg['n_pole_grid'],
        verbose=verbose,
        seed=seed,
    )

    result['convex_result'] = {
        'pole_lambda_deg': convex_result['pole_lambda'] * RAD2DEG,
        'pole_beta_deg': convex_result['pole_beta'] * RAD2DEG,
        'chi2': float(convex_result['chi2']),
        'chi2_reduced': float(convex_result['chi2_reduced']),
    }

    # ===== Step 3: Genetic Non-Convex Refinement =====
    gcfg = config['genetic_refinement']
    if gcfg['enabled']:
        if verbose:
            print(f"\nStep 3: Genetic Non-Convex Refinement")

        # Use convex solution as initial radii
        convex_model = convex_result['model']
        init_radii = np.linalg.norm(convex_model.vertices, axis=1)
        # Resample to genetic solver mesh resolution
        from ..shapes.convex_model import create_icosphere
        gen_verts, gen_faces = create_icosphere(gcfg['subdivisions'])
        # Interpolate radii from convex to genetic mesh
        from scipy.spatial import cKDTree
        tree = cKDTree(convex_model.vertices / np.linalg.norm(convex_model.vertices, axis=1, keepdims=True))
        _, indices = tree.query(gen_verts)
        init_radii_gen = np.linalg.norm(convex_model.vertices[indices], axis=1)

        genetic_result = genetic_inversion(
            lightcurves, period_days,
            convex_result['pole_lambda'], convex_result['pole_beta'],
            phi0=convex_result['phi0'],
            subdivisions=gcfg['subdivisions'],
            population_size=gcfg['population_size'],
            generations=gcfg['generations'],
            mutation_rate=gcfg['mutation_rate'],
            init_radii=init_radii_gen,
            lambda_smooth=gcfg['lambda_smooth'],
            verbose=verbose,
            seed=seed,
        )

        result['genetic_result'] = {
            'chi2': float(genetic_result['chi2']),
            'chi2_reduced': float(genetic_result['chi2_reduced']),
        }

        # Choose best result
        if genetic_result['chi2'] < convex_result['chi2']:
            final_mesh = genetic_result['mesh']
            final_chi2 = genetic_result['chi2_reduced']
            result['final_method'] = 'genetic'
            if verbose:
                print(f"  Genetic solver improved chi2: "
                      f"{convex_result['chi2']:.2f} -> {genetic_result['chi2']:.2f}")
        else:
            final_mesh = convex_result['model']
            final_chi2 = convex_result['chi2_reduced']
            result['final_method'] = 'convex'
    else:
        final_mesh = convex_result['model']
        final_chi2 = convex_result['chi2_reduced']
        result['final_method'] = 'convex'
        genetic_result = None

    # ===== Step 4: Sparse Data Fusion (optional) =====
    scfg = config['sparse_fusion']
    if scfg['enabled'] and sparse_data:
        if verbose:
            print(f"\nStep 4: Sparse Data Fusion")
        # Run sparse inversion starting from current solution
        sparse_result = sparse_inversion(
            sparse_data, period_days,
            pole_lambda_init=convex_result['pole_lambda'],
            pole_beta_init=convex_result['pole_beta'],
            verbose=verbose, seed=seed,
        )
        result['sparse_result'] = {
            'chi2': float(sparse_result['chi2']),
            'chi2_reduced': float(sparse_result['chi2_reduced']),
        }

    # ===== Finalize =====
    result['pole_lambda_deg'] = convex_result['pole_lambda'] * RAD2DEG
    result['pole_beta_deg'] = convex_result['pole_beta'] * RAD2DEG
    result['chi2_reduced'] = float(final_chi2)
    result['converged'] = final_chi2 < 10.0

    # Save outputs
    if output_dir and asteroid_id:
        os.makedirs(output_dir, exist_ok=True)

        if config['output']['save_obj']:
            obj_path = os.path.join(output_dir, f"{asteroid_id}.obj")
            final_mesh.save_obj(obj_path)

        if config['output']['save_spin']:
            spin_path = os.path.join(output_dir, f"{asteroid_id}_spin.json")
            spin_data = {
                'asteroid_id': asteroid_id,
                'pole_lambda_deg': result['pole_lambda_deg'],
                'pole_beta_deg': result['pole_beta_deg'],
                'period_hours': result['period_hours'],
                'chi2_reduced': result['chi2_reduced'],
                'method': result['final_method'],
            }
            with open(spin_path, 'w') as f:
                json.dump(spin_data, f, indent=2)

    if verbose:
        print(f"\nResult: period={result['period_hours']:.4f}h, "
              f"pole=({result['pole_lambda_deg']:.1f}, {result['pole_beta_deg']:.1f}), "
              f"chi2_red={result['chi2_reduced']:.3f}, "
              f"method={result.get('final_method', 'N/A')}")

    return result
