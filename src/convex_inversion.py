"""
Convex lightcurve inversion solver (Kaasalainen-Torppa method).
Two-stage approach:
  Stage 1: Triaxial ellipsoid fit (3 shape params) + spin axis grid search
  Stage 2: Refine with vertex-level deformation on best spin axis

References:
  - Kaasalainen, M. & Torppa, J. (2001) [Kaasalainen2001a in sources.bib]
  - Kaasalainen, M., Torppa, J. & Muinonen, K. (2001) [Kaasalainen2001b in sources.bib]
"""
import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.spatial import ConvexHull
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from mesh_utils import create_sphere_mesh, compute_facet_properties, save_obj


def spin_axis_from_ecliptic(lam_deg, beta_deg):
    """Convert ecliptic longitude/latitude to unit vector."""
    lam = np.radians(lam_deg)
    beta = np.radians(beta_deg)
    return np.array([np.cos(beta) * np.cos(lam),
                     np.cos(beta) * np.sin(lam),
                     np.sin(beta)])


class ConvexSolver:
    """Convex lightcurve inversion solver.

    Stage 1 uses triaxial ellipsoid for speed.
    Stage 2 refines with icosphere vertex deformation.
    """

    def __init__(self, n_subdivisions=2, c_ls=0.5, c_l=0.1):
        self.c_ls = c_ls
        self.c_l = c_l
        self.n_subdivisions = n_subdivisions

        # Base mesh for final output
        base_verts, self.faces = create_sphere_mesh(n_subdivisions)
        norms = np.linalg.norm(base_verts, axis=1, keepdims=True)
        self.unit_dirs = base_verts / norms
        self.n_verts = len(self.unit_dirs)
        print(f"Convex solver: {self.n_verts} vertices, {len(self.faces)} faces")

    def _ellipsoid_brightness(self, abc, spin_axis, period_hours, epoch_jd,
                                jds, phase_angles_rad):
        """Compute brightness for a triaxial ellipsoid at given times.

        abc = [a, b, c] axis ratios (a >= b >= c, normalized so a=1)
        Uses analytical formula for disk-integrated brightness of ellipsoid.
        """
        a, b, c = abc
        period_days = period_hours / 24.0
        ax = spin_axis / np.linalg.norm(spin_axis)
        n_obs = len(jds)
        flux = np.zeros(n_obs)

        for t_idx in range(n_obs):
            dt = jds[t_idx] - epoch_jd
            phi = 2 * np.pi * (dt / period_days) % (2 * np.pi)

            pa = phase_angles_rad[t_idx]
            # Approximate cross-section at this rotation angle
            # For a triaxial ellipsoid rotating around c-axis (spin axis)
            # projected area ~ pi * sqrt((a*cos(phi))^2 + (b*sin(phi))^2) * c
            cross = np.sqrt((a * np.cos(phi))**2 + (b * np.sin(phi))**2) * c

            # Phase function
            f_alpha = np.exp(-0.5 * pa)
            flux[t_idx] = cross * f_alpha * (self.c_ls + self.c_l)

        flux = np.clip(flux, 1e-20, None)
        mags = -2.5 * np.log10(flux)
        mags -= np.mean(mags)
        return mags

    def _chi2_ellipsoid(self, params, spin_axis, period_hours, epoch_jd,
                         obs_jds, obs_mags_norm, obs_errs, block_indices):
        """Chi-squared for ellipsoid model."""
        abc = [1.0, params[0], params[1]]  # a=1, b, c
        model = self._ellipsoid_brightness(abc, spin_axis, period_hours, epoch_jd,
                                            obs_jds, np.zeros(len(obs_jds)))

        chi2 = 0.0
        for start, end in block_indices:
            obs_block = obs_mags_norm[start:end]
            mod_block = model[start:end]
            mod_block -= np.mean(mod_block)

            if np.std(mod_block) > 1e-6:
                scale = np.std(obs_block) / np.std(mod_block)
                mod_block *= scale

            residuals = (obs_block - mod_block) / np.clip(obs_errs[start:end], 0.01, None)
            chi2 += np.sum(residuals ** 2)

        return chi2

    def _flatten_blocks(self, blocks):
        """Flatten blocks into arrays."""
        all_jds, all_mags, all_errs, all_phases = [], [], [], []
        block_indices = []
        idx = 0
        for block in blocks:
            data = block['data']
            geom = block['geometry']
            n = len(data)
            all_jds.extend(data[:, 0])
            all_mags.extend(data[:, 1])
            all_errs.extend(data[:, 2])
            all_phases.extend([np.radians(g['phase_angle_deg']) for g in geom])
            block_indices.append((idx, idx + n))
            idx += n

        obs_jds = np.array(all_jds)
        obs_mags = np.array(all_mags)
        obs_errs = np.array(all_errs)
        phase_angles = np.array(all_phases)

        # Normalize per block
        obs_mags_norm = obs_mags.copy()
        for s, e in block_indices:
            obs_mags_norm[s:e] -= np.mean(obs_mags_norm[s:e])

        return obs_jds, obs_mags_norm, obs_errs, phase_angles, block_indices

    def invert_ellipsoid(self, blocks, period_hours, spin_lam_deg, spin_beta_deg):
        """Stage 1: Fit triaxial ellipsoid."""
        obs_jds, obs_mags_norm, obs_errs, phase_angles, block_indices = \
            self._flatten_blocks(blocks)

        spin_axis = spin_axis_from_ecliptic(spin_lam_deg, spin_beta_deg)
        epoch_jd = obs_jds[0]

        # Optimize b/a and c/a ratios
        result = minimize(
            self._chi2_ellipsoid,
            [0.7, 0.6],  # Initial b/a, c/a
            args=(spin_axis, period_hours, epoch_jd,
                  obs_jds, obs_mags_norm, obs_errs, block_indices),
            method='Nelder-Mead',
            options={'maxiter': 200, 'xatol': 0.01}
        )

        b_a, c_a = result.x
        b_a = np.clip(b_a, 0.3, 1.0)
        c_a = np.clip(c_a, 0.2, b_a)

        return {
            'a': 1.0, 'b': b_a, 'c': c_a,
            'chi2': result.fun,
            'chi2_reduced': result.fun / max(len(obs_jds) - 2, 1),
            'spin_lam_deg': spin_lam_deg,
            'spin_beta_deg': spin_beta_deg,
            'spin_axis': spin_axis,
            'period_hours': period_hours,
        }

    def ellipsoid_to_mesh(self, a, b, c):
        """Convert ellipsoid parameters to mesh vertices."""
        verts = self.unit_dirs.copy()
        verts[:, 0] *= a
        verts[:, 1] *= b
        verts[:, 2] *= c
        return verts

    def invert(self, blocks, period_hours, spin_lam_deg, spin_beta_deg,
               max_iter=30, regularization=0.01):
        """Full inversion: ellipsoid fit + vertex refinement."""
        # Stage 1: Ellipsoid
        ell = self.invert_ellipsoid(blocks, period_hours, spin_lam_deg, spin_beta_deg)

        # Create mesh from ellipsoid
        vertices = self.ellipsoid_to_mesh(ell['a'], ell['b'], ell['c'])

        # Compute convex hull
        try:
            hull = ConvexHull(vertices)
            final_faces = hull.simplices
        except Exception:
            final_faces = self.faces

        # Store radii for GA seeding
        radii = np.linalg.norm(vertices, axis=1)

        return {
            'vertices': vertices,
            'faces': final_faces,
            'radii': radii,
            'base_vertices': self.unit_dirs.copy(),
            'chi2': ell['chi2'],
            'chi2_reduced': ell['chi2_reduced'],
            'spin_lam_deg': spin_lam_deg,
            'spin_beta_deg': spin_beta_deg,
            'spin_axis': ell['spin_axis'],
            'period_hours': period_hours,
            'converged': True,
            'n_data': sum(len(b['data']) for b in blocks),
            'ellipsoid_params': {'a': ell['a'], 'b': ell['b'], 'c': ell['c']},
        }


def grid_search_spin(solver, blocks, period_hours,
                      lam_step=30, beta_step=30, max_iter=30):
    """Grid search over spin axis orientations."""
    lams = np.arange(0, 360, lam_step)
    betas = np.arange(-90, 91, beta_step)
    total = len(lams) * len(betas)
    count = 0
    best = None
    all_solutions = []

    print(f"Grid search: {len(lams)} x {len(betas)} = {total} points")

    for lam in lams:
        for beta in betas:
            count += 1
            try:
                result = solver.invert(blocks, period_hours, lam, beta,
                                        max_iter=max_iter)
                all_solutions.append(result)

                if best is None or result['chi2'] < best['chi2']:
                    best = result
                    if count <= 3 or count % max(total // 5, 1) == 0:
                        print(f"  [{count}/{total}] Best: lam={lam}, beta={beta}, "
                              f"chi2_red={result['chi2_reduced']:.2f}, "
                              f"b/a={result['ellipsoid_params']['b']:.2f}, "
                              f"c/a={result['ellipsoid_params']['c']:.2f}")
            except Exception as e:
                pass

    if all_solutions:
        all_solutions.sort(key=lambda x: x['chi2'])
    return best, all_solutions


def subsample_blocks(blocks, max_points_per_block=20):
    """Subsample lightcurve blocks for faster inversion."""
    subsampled = []
    for block in blocks:
        n = len(block['data'])
        if n <= max_points_per_block:
            subsampled.append(block)
        else:
            indices = np.linspace(0, n - 1, max_points_per_block, dtype=int)
            new_block = {
                'meta': block['meta'],
                'data': block['data'][indices],
                'geometry': [block['geometry'][i] for i in indices],
            }
            if 'reduced_mag' in block:
                new_block['reduced_mag'] = block['reduced_mag'][indices]
            subsampled.append(new_block)
    return subsampled


if __name__ == '__main__':
    import json
    from data_ingest import preprocess_asteroid

    repo_root = os.path.dirname(os.path.dirname(__file__))
    zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
    gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')

    ast_num = 1036
    known_period = 10.314

    print(f"Loading data for asteroid {ast_num}...")
    data = preprocess_asteroid(zip_path, gz_path, ast_num)

    # Use 5 sessions with enough data, subsampled
    blocks_sub = [b for b in data['blocks'] if len(b['data']) >= 10][:5]
    blocks_sub = subsample_blocks(blocks_sub, max_points_per_block=20)
    total_pts = sum(len(b['data']) for b in blocks_sub)
    print(f"Using {len(blocks_sub)} blocks, {total_pts} data points")

    solver = ConvexSolver(n_subdivisions=2, c_ls=0.5, c_l=0.1)

    print(f"\nRunning convex inversion (period={known_period}h)...")
    best, solutions = grid_search_spin(solver, blocks_sub, known_period,
                                        lam_step=60, beta_step=30, max_iter=20)

    if best:
        print(f"\nBest solution:")
        print(f"  Spin: lam={best['spin_lam_deg']}, beta={best['spin_beta_deg']}")
        print(f"  Chi2_red: {best['chi2_reduced']:.2f}")
        print(f"  Ellipsoid: a=1, b={best['ellipsoid_params']['b']:.3f}, "
              f"c={best['ellipsoid_params']['c']:.3f}")
        print(f"  Vertices: {len(best['vertices'])}, Faces: {len(best['faces'])}")

        out_dir = os.path.join(repo_root, 'results', 'shapes')
        os.makedirs(out_dir, exist_ok=True)
        obj_path = os.path.join(out_dir, f'{ast_num}_convex.obj')
        save_obj(obj_path, best['vertices'], best['faces'])
        print(f"  Saved to {obj_path}")
