"""
Genetic algorithm non-convex shape solver for asteroid lightcurve inversion.

Encodes vertex radial displacements from a convex seed mesh as the genome,
allowing concavities (inward deformation) and protrusions (outward deformation).
Fitness is evaluated using the self-shadowing forward photometric model from
scattering.py, which is critical for correctly modelling non-convex shapes
where one part of the body may cast a shadow on another.

The GA operators -- tournament selection, two-point crossover, Gaussian
mutation, and elitism -- follow the SAGE approach described in:

    Bartczak, P. & Dudzinski, G. (2018), "SAGE modelling of global
    self-shadowing effects in asteroid shape reconstruction from
    lightcurves", MNRAS, 473, 5085-5098.  [Bartczak2018 in sources.bib]

References:
  - Bartczak & Dudzinski (2018) [Bartczak2018 in sources.bib]
  - Kaasalainen & Torppa (2001) [Kaasalainen2001a in sources.bib]
"""
import numpy as np
import os
import sys
import copy

sys.path.insert(0, os.path.dirname(__file__))

from mesh_utils import create_sphere_mesh, compute_facet_properties, save_obj
from scattering import (
    compute_brightness,
    compute_shadow_mask_vectorized,
    _rotate_vertices,
    build_bvh,
)
from convex_inversion import ConvexSolver, subsample_blocks


# ============================================================
# Helper: build deformed mesh from genome
# ============================================================

def _apply_genome(unit_dirs, radii_genome, faces):
    """Deform a unit-sphere mesh by applying per-vertex radii.

    Parameters
    ----------
    unit_dirs : (N, 3) unit direction vectors (from icosphere)
    radii_genome : (N,) radial distance per vertex
    faces : (M, 3) triangle connectivity

    Returns
    -------
    vertices : (N, 3)
    normals : (M, 3)
    areas : (M,)
    centers : (M, 3)
    """
    vertices = unit_dirs * radii_genome[:, np.newaxis]
    normals, areas, centers = compute_facet_properties(vertices, faces)
    return vertices, normals, areas, centers


# ============================================================
# Smoothness regularisation
# ============================================================

def _build_adjacency(n_verts, faces):
    """Build a vertex-adjacency list from a triangle mesh.

    Returns a list-of-sets where adjacency[i] contains the indices of all
    vertices sharing an edge with vertex i.
    """
    adjacency = [set() for _ in range(n_verts)]
    for f in faces:
        for a, b in [(f[0], f[1]), (f[1], f[2]), (f[2], f[0])]:
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def _smoothness_penalty(radii, adjacency):
    """Laplacian smoothness penalty: sum of squared differences between
    each vertex radius and the mean of its neighbours.

    A small penalty encourages locally smooth shapes and prevents
    high-frequency noise in the genome from producing unrealistic
    spiky surfaces.
    """
    penalty = 0.0
    for i, neighbours in enumerate(adjacency):
        if not neighbours:
            continue
        neighbour_mean = np.mean([radii[j] for j in neighbours])
        penalty += (radii[i] - neighbour_mean) ** 2
    return penalty


# ============================================================
# GASolver
# ============================================================

class GASolver:
    """Genetic-algorithm solver for non-convex asteroid shape inversion.

    The genome of each individual is an array of per-vertex radii that
    are applied to the unit directions of an icosphere mesh.  The convex
    seed from ``ConvexSolver`` provides the initial population centre;
    mutation then explores both inward (concavity) and outward
    perturbations.

    Parameters
    ----------
    n_verts : int
        Number of mesh vertices (must match ``base_vertices``).
    base_vertices : (N, 3) array
        Unit-direction vectors of the seed icosphere.
    faces : (M, 3) array of int
        Triangle connectivity (0-indexed).
    c_ls : float
        Lommel-Seeliger scattering weight.
    c_l : float
        Lambert scattering weight.
    """

    def __init__(self, n_verts, base_vertices, faces, c_ls=0.5, c_l=0.1):
        self.n_verts = n_verts
        # Ensure unit directions
        norms = np.linalg.norm(base_vertices, axis=1, keepdims=True)
        self.unit_dirs = base_vertices / np.clip(norms, 1e-10, None)
        self.faces = np.asarray(faces, dtype=np.int32)
        self.c_ls = c_ls
        self.c_l = c_l

        # Pre-build vertex adjacency for smoothness regularisation
        self.adjacency = _build_adjacency(self.n_verts, self.faces)

    # ----------------------------------------------------------
    # Fitness evaluation
    # ----------------------------------------------------------

    def fitness(self, individual, blocks, spin_axis, period_hours, epoch_jd,
                regularization=0.01):
        """Evaluate chi-squared fitness of a single genome.

        Parameters
        ----------
        individual : (N,) array
            Per-vertex radii (the genome).
        blocks : list of dict
            Preprocessed lightcurve blocks (subsampled for speed).
        spin_axis : (3,) unit vector
            Rotation pole in ecliptic coordinates.
        period_hours : float
            Sidereal rotation period.
        epoch_jd : float
            Julian Date of zero rotation phase.
        regularization : float
            Weight of the Laplacian smoothness penalty.

        Returns
        -------
        float
            Total chi-squared (lower is better).
        """
        radii = individual
        vertices, normals, areas, centers = _apply_genome(
            self.unit_dirs, radii, self.faces
        )

        period_days = period_hours / 24.0
        chi2 = 0.0

        for block in blocks:
            data = block['data']          # (n, 3): JD, mag, err
            geom = block['geometry']      # list of dicts per point
            n_pts = len(data)

            obs_mags = data[:, 1].copy()
            obs_errs = data[:, 2].copy()
            obs_mags -= np.mean(obs_mags)  # zero-mean relative magnitudes

            model_flux = np.zeros(n_pts)
            for k in range(n_pts):
                jd_k = data[k, 0]
                dt = jd_k - epoch_jd
                rotation_angle = (2.0 * np.pi * dt / period_days) % (2.0 * np.pi)

                rot_verts = _rotate_vertices(vertices, spin_axis, rotation_angle)
                rot_normals, rot_areas, rot_centers = compute_facet_properties(
                    rot_verts, self.faces
                )

                # Sun and observer directions from geometry metadata
                # We approximate sun/obs directions from the phase angle
                # stored in the geometry dict, following the same convention
                # as convex_inversion.py.
                phase_rad = np.radians(geom[k]['phase_angle_deg'])

                # Sun direction: +x axis (arbitrary reference, absorbed by
                # the zero-mean normalisation of relative magnitudes).
                sun_dir = np.array([1.0, 0.0, 0.0], dtype=np.float64)
                # Observer direction: rotated from sun direction by phase angle
                obs_dir = np.array([np.cos(phase_rad), np.sin(phase_rad), 0.0],
                                   dtype=np.float64)

                flux = compute_brightness(
                    rot_verts, self.faces, rot_normals, rot_areas, rot_centers,
                    sun_dir, obs_dir, bvh_root=None,
                    c_ls=self.c_ls, c_l=self.c_l,
                    use_shadow=True,
                )
                model_flux[k] = flux

            # Convert model fluxes to relative magnitudes
            model_flux = np.clip(model_flux, 1e-20, None)
            model_mags = -2.5 * np.log10(model_flux)
            model_mags -= np.mean(model_mags)

            # Scale model amplitude to match observed amplitude (free
            # albedo / distance factor absorbed per block).
            std_obs = np.std(obs_mags)
            std_mod = np.std(model_mags)
            if std_mod > 1e-8:
                model_mags *= std_obs / std_mod

            residuals = (obs_mags - model_mags) / np.clip(obs_errs, 0.01, None)
            chi2 += np.sum(residuals ** 2)

        # Add smoothness regularisation penalty
        chi2 += regularization * _smoothness_penalty(radii, self.adjacency)

        return chi2

    # ----------------------------------------------------------
    # Selection
    # ----------------------------------------------------------

    @staticmethod
    def tournament_selection(population, fitnesses, tournament_size=3):
        """Select one individual via tournament selection.

        Parameters
        ----------
        population : list of (N,) arrays
        fitnesses : (P,) array
        tournament_size : int

        Returns
        -------
        (N,) array  -- copy of the winner's genome.
        """
        indices = np.random.choice(len(population), size=tournament_size,
                                   replace=False)
        best_idx = indices[np.argmin(fitnesses[indices])]
        return population[best_idx].copy()

    # ----------------------------------------------------------
    # Crossover
    # ----------------------------------------------------------

    @staticmethod
    def crossover(parent_a, parent_b):
        """Two-point crossover.

        Two random cut points are chosen along the genome; the segment
        between them is swapped between the two parents to produce two
        offspring.

        Parameters
        ----------
        parent_a, parent_b : (N,) arrays

        Returns
        -------
        child_a, child_b : (N,) arrays
        """
        n = len(parent_a)
        pt1, pt2 = sorted(np.random.choice(n, size=2, replace=False))
        child_a = parent_a.copy()
        child_b = parent_b.copy()
        child_a[pt1:pt2] = parent_b[pt1:pt2]
        child_b[pt1:pt2] = parent_a[pt1:pt2]
        return child_a, child_b

    # ----------------------------------------------------------
    # Mutation
    # ----------------------------------------------------------

    @staticmethod
    def mutation(individual, mutation_rate=0.1, sigma=0.05, r_min=0.3):
        """Gaussian mutation of individual vertex radii.

        Each gene (vertex radius) has an independent probability
        ``mutation_rate`` of being perturbed by a zero-mean Gaussian
        with standard deviation ``sigma``.  A lower clamp ``r_min``
        prevents the vertex from collapsing through the origin.

        Parameters
        ----------
        individual : (N,) array
        mutation_rate : float in [0, 1]
        sigma : float
        r_min : float
            Minimum allowed radius.

        Returns
        -------
        (N,) array  -- mutated genome (in-place modified copy).
        """
        child = individual.copy()
        mask = np.random.rand(len(child)) < mutation_rate
        child[mask] += np.random.randn(np.sum(mask)) * sigma
        child = np.clip(child, r_min, None)
        return child

    # ----------------------------------------------------------
    # Main evolution loop
    # ----------------------------------------------------------

    def evolve(self, blocks, spin_axis, period_hours,
               pop_size=50, n_generations=100, mutation_rate=0.1,
               mutation_sigma=0.05, tournament_size=3,
               elitism_count=2, regularization=0.01,
               seed_radii=None, verbose=True):
        """Run the genetic algorithm.

        Parameters
        ----------
        blocks : list of dict
            Lightcurve observation blocks (should be subsampled via
            ``subsample_blocks`` for speed).
        spin_axis : (3,) array
            Unit rotation pole vector.
        period_hours : float
            Sidereal rotation period in hours.
        pop_size : int
            Population size.
        n_generations : int
            Number of generations.
        mutation_rate : float
            Per-gene probability of Gaussian mutation.
        mutation_sigma : float
            Standard deviation of the Gaussian mutation noise.
        tournament_size : int
            Number of contestants in tournament selection.
        elitism_count : int
            Number of elite individuals carried unchanged into the next
            generation.
        regularization : float
            Weight of smoothness regularisation in the fitness function.
        seed_radii : (N,) array or None
            Radii from the convex seed solution.  If *None*, all radii
            default to 1.0 (unit sphere).
        verbose : bool
            Print progress each generation.

        Returns
        -------
        dict with keys
            'vertices' : (N, 3)   best-fit deformed mesh vertices
            'faces'    : (M, 3)   triangle connectivity
            'chi2'     : float    best chi-squared value
            'log'      : list     best fitness per generation
        """
        spin_axis = np.asarray(spin_axis, dtype=np.float64)
        spin_axis /= np.linalg.norm(spin_axis)
        epoch_jd = blocks[0]['data'][0, 0]  # first observation as epoch

        # ---- initialise population ----
        if seed_radii is not None:
            base_radii = np.asarray(seed_radii, dtype=np.float64)
        else:
            base_radii = np.ones(self.n_verts, dtype=np.float64)

        population = []
        # First individual is the seed itself (elitism preserves it if good)
        population.append(base_radii.copy())
        for _ in range(pop_size - 1):
            noise = np.random.randn(self.n_verts) * mutation_sigma
            ind = np.clip(base_radii + noise, 0.3, None)
            population.append(ind)

        # ---- evaluate initial population ----
        fitnesses = np.array([
            self.fitness(ind, blocks, spin_axis, period_hours, epoch_jd,
                         regularization=regularization)
            for ind in population
        ])

        convergence_log = []
        best_ever_fitness = np.min(fitnesses)
        best_ever_idx = int(np.argmin(fitnesses))
        best_ever_genome = population[best_ever_idx].copy()

        if verbose:
            print(f"GA initialised: pop_size={pop_size}, "
                  f"n_verts={self.n_verts}, n_gen={n_generations}")
            print(f"  Initial best chi2 = {best_ever_fitness:.4f}")

        # ---- generational loop ----
        for gen in range(n_generations):
            # Sort by fitness (ascending = best first)
            order = np.argsort(fitnesses)
            sorted_pop = [population[i] for i in order]
            sorted_fit = fitnesses[order]

            new_population = []
            # Elitism: carry forward best individuals unchanged
            for e in range(min(elitism_count, pop_size)):
                new_population.append(sorted_pop[e].copy())

            # Fill the rest via selection, crossover, mutation
            while len(new_population) < pop_size:
                parent_a = self.tournament_selection(
                    population, fitnesses, tournament_size
                )
                parent_b = self.tournament_selection(
                    population, fitnesses, tournament_size
                )
                child_a, child_b = self.crossover(parent_a, parent_b)
                child_a = self.mutation(child_a, mutation_rate, mutation_sigma)
                child_b = self.mutation(child_b, mutation_rate, mutation_sigma)

                new_population.append(child_a)
                if len(new_population) < pop_size:
                    new_population.append(child_b)

            population = new_population

            # Evaluate new population
            fitnesses = np.array([
                self.fitness(ind, blocks, spin_axis, period_hours, epoch_jd,
                             regularization=regularization)
                for ind in population
            ])

            gen_best = np.min(fitnesses)
            convergence_log.append(float(gen_best))

            if gen_best < best_ever_fitness:
                best_ever_fitness = gen_best
                best_ever_idx = int(np.argmin(fitnesses))
                best_ever_genome = population[best_ever_idx].copy()

            if verbose:
                print(f"  Gen {gen + 1:4d}/{n_generations}  "
                      f"best={gen_best:.4f}  overall_best={best_ever_fitness:.4f}")

        # ---- build output mesh from best genome ----
        best_verts, _, _, _ = _apply_genome(
            self.unit_dirs, best_ever_genome, self.faces
        )

        return {
            'vertices': best_verts,
            'faces': self.faces,
            'radii': best_ever_genome,
            'chi2': float(best_ever_fitness),
            'log': convergence_log,
        }


# ============================================================
# Quick integration test
# ============================================================

if __name__ == '__main__':
    repo_root = os.path.dirname(os.path.dirname(__file__))
    sys.path.insert(0, os.path.join(repo_root, 'src'))

    from data_ingest import preprocess_asteroid

    zip_path = os.path.join(repo_root, 'ALCDEF_ALL.zip')
    gz_path = os.path.join(repo_root, 'MPCORB.DAT.gz')
    ast_num = 1036          # Ganymed
    known_period = 10.314   # hours (LCDB)

    # ---- Load and subsample data ----
    print(f"Loading ALCDEF data for asteroid {ast_num} (Ganymed)...")
    data = preprocess_asteroid(zip_path, gz_path, ast_num)

    blocks = [b for b in data['blocks'] if len(b['data']) >= 10][:5]
    blocks = subsample_blocks(blocks, max_points_per_block=10)
    total_pts = sum(len(b['data']) for b in blocks)
    print(f"Using {len(blocks)} blocks, {total_pts} data points "
          f"(subsampled, max 10/block)")

    # ---- Stage 1: Convex seed (42-vertex mesh, 1 subdivision) ----
    print("\n--- Stage 1: Convex inversion (seed) ---")
    convex_solver = ConvexSolver(n_subdivisions=1, c_ls=0.5, c_l=0.1)

    # Use a coarse spin grid to keep the test fast
    from convex_inversion import grid_search_spin
    convex_best, _ = grid_search_spin(
        convex_solver, blocks, known_period,
        lam_step=60, beta_step=30, max_iter=10,
    )

    if convex_best is None:
        print("Convex inversion did not converge -- using unit sphere seed.")
        base_verts, faces = create_sphere_mesh(n_subdivisions=1)
        seed_radii = np.ones(len(base_verts))
        norms = np.linalg.norm(base_verts, axis=1, keepdims=True)
        unit_dirs = base_verts / norms
        spin_axis = np.array([0.0, 0.0, 1.0])
        convex_chi2 = float('inf')
    else:
        print(f"Convex seed: spin=({convex_best['spin_lam_deg']}, "
              f"{convex_best['spin_beta_deg']}), "
              f"chi2_red={convex_best['chi2_reduced']:.2f}")
        seed_radii = convex_best['radii']
        unit_dirs = convex_best['base_vertices']
        faces = convex_best['faces']
        spin_axis = convex_best['spin_axis']
        convex_chi2 = convex_best['chi2']

    n_verts = len(unit_dirs)

    # ---- Stage 2: GA non-convex refinement ----
    print(f"\n--- Stage 2: GA non-convex refinement ({n_verts} vertices) ---")
    ga = GASolver(
        n_verts=n_verts,
        base_vertices=unit_dirs,
        faces=faces,
        c_ls=0.5,
        c_l=0.1,
    )

    result = ga.evolve(
        blocks=blocks,
        spin_axis=spin_axis,
        period_hours=known_period,
        pop_size=20,
        n_generations=20,
        mutation_rate=0.1,
        mutation_sigma=0.05,
        tournament_size=3,
        elitism_count=2,
        regularization=0.01,
        seed_radii=seed_radii,
        verbose=True,
    )

    # ---- Save output ----
    out_dir = os.path.join(repo_root, 'results', 'shapes')
    os.makedirs(out_dir, exist_ok=True)
    obj_path = os.path.join(out_dir, '1036_ga.obj')
    save_obj(obj_path, result['vertices'], result['faces'])
    print(f"\nSaved GA shape model to {obj_path}")

    # ---- Report chi2 improvement ----
    ga_chi2 = result['chi2']
    print(f"\nChi-squared comparison:")
    print(f"  Convex seed chi2 = {convex_chi2:.4f}")
    print(f"  GA best chi2     = {ga_chi2:.4f}")
    if convex_chi2 > 0 and np.isfinite(convex_chi2):
        improvement = (convex_chi2 - ga_chi2) / convex_chi2 * 100
        print(f"  Improvement       = {improvement:+.2f}%")
    print(f"\nConvergence log (best chi2 per generation):")
    for i, val in enumerate(result['log']):
        print(f"  gen {i + 1:3d}: {val:.4f}")
