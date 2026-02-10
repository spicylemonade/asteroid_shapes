"""
Genetic/Evolutionary Non-Convex Shape Solver (SAGE-Inspired)

Implements a genetic algorithm for asteroid shape inversion that allows
non-convex shapes (concavities, bifurcations). Based on the methodology
described in Bartczak & Dudzinski (2018) SAGE.
"""

import numpy as np
from ..shapes.convex_model import (
    create_icosphere, lommel_seeliger_lambert, spin_rotation_matrix
)

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


class NonConvexMesh:
    """Vertex-based mesh allowing non-convex shapes."""

    def __init__(self, subdivisions=2):
        self.base_verts, self.faces = create_icosphere(subdivisions)
        self.n_vertices = len(self.base_verts)
        self.n_facets = len(self.faces)
        # Directions for radial parameterization
        self.directions = self.base_verts.copy()
        # Radii (one per vertex)
        self.radii = np.ones(self.n_vertices)
        self._update()

    def set_radii(self, radii):
        self.radii = np.clip(radii, 0.1, 10.0)
        self._update()

    def _update(self):
        self.vertices = self.directions * self.radii[:, np.newaxis]
        v0 = self.vertices[self.faces[:, 0]]
        v1 = self.vertices[self.faces[:, 1]]
        v2 = self.vertices[self.faces[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)
        self.facet_areas = 0.5 * np.linalg.norm(cross, axis=1)
        norms = np.linalg.norm(cross, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-30)
        self.facet_normals = cross / norms
        self.facet_centroids = (v0 + v1 + v2) / 3.0

    def save_obj(self, filepath):
        with open(filepath, 'w') as f:
            f.write("# Non-convex asteroid shape model (genetic solver)\n")
            for v in self.vertices:
                f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            for face in self.faces:
                f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")


def compute_brightness_nonconvex(mesh, sun_dir, obs_dir, phase_angle=None):
    """Compute brightness for non-convex mesh (simplified, no ray-casting)."""
    mu0 = np.dot(mesh.facet_normals, sun_dir)
    mu = np.dot(mesh.facet_normals, obs_dir)
    visible = (mu > 0) & (mu0 > 0)

    if phase_angle is None:
        cos_alpha = np.clip(np.dot(sun_dir, obs_dir), -1, 1)
        phase_angle = np.arccos(cos_alpha)

    scatter = lommel_seeliger_lambert(mu, mu0, phase_angle)
    return np.sum(mesh.facet_areas[visible] * scatter[visible])


def synthetic_lc_nonconvex(mesh, times, pole_lambda, pole_beta, period,
                           phi0, t0, sun_dirs, obs_dirs, phase_angles):
    """Compute synthetic lightcurve for non-convex mesh."""
    n = len(times)
    flux = np.zeros(n)
    for i in range(n):
        R = spin_rotation_matrix(pole_lambda, pole_beta, period, phi0, times[i], t0)
        sun_body = R @ sun_dirs[i]
        obs_body = R @ obs_dirs[i]
        flux[i] = compute_brightness_nonconvex(mesh, sun_body, obs_body, phase_angles[i])
    flux = np.maximum(flux, 1e-30)
    return -2.5 * np.log10(flux)


def fitness_function(radii, mesh, lightcurves, pole_lambda, pole_beta,
                     period, phi0, t0, lambda_smooth=0.01, lambda_regularity=0.01):
    """Compute fitness (negative chi-squared) for genetic algorithm."""
    mesh.set_radii(radii)

    chi2 = 0.0
    n_data = 0
    for lc in lightcurves:
        syn = synthetic_lc_nonconvex(
            mesh, lc['times'], pole_lambda, pole_beta, period,
            phi0, t0, lc['sun_dirs'], lc['obs_dirs'], lc['phase_angles'])
        offset = np.mean(lc['mags'] - syn)
        residuals = (lc['mags'] - syn - offset) / lc['errors']
        chi2 += np.sum(residuals**2)
        n_data += len(lc['times'])

    # Smoothness: penalize difference between adjacent vertices
    smooth_penalty = 0.0
    for face in mesh.faces:
        for j in range(3):
            i1, i2 = face[j], face[(j+1) % 3]
            smooth_penalty += (radii[i1] - radii[i2])**2
    smooth_penalty *= lambda_smooth

    # Regularity: penalize variation in facet areas
    areas = mesh.facet_areas
    mean_area = np.mean(areas)
    regularity = lambda_regularity * np.sum((areas - mean_area)**2) / (mean_area**2 + 1e-10)

    total = chi2 + smooth_penalty + regularity
    return -total, chi2, n_data


def genetic_inversion(lightcurves, period, pole_lambda, pole_beta,
                      phi0=0.0, t0=2451545.0, subdivisions=2,
                      population_size=50, generations=100,
                      mutation_rate=0.15, crossover_rate=0.7,
                      elite_fraction=0.1, init_radii=None,
                      lambda_smooth=0.01, lambda_regularity=0.01,
                      verbose=True, seed=42):
    """Run genetic algorithm for non-convex shape inversion.

    Parameters
    ----------
    lightcurves : list of dicts with 'times','mags','errors','sun_dirs','obs_dirs','phase_angles'
    period : float, rotation period in days
    pole_lambda, pole_beta : float, pole direction (radians)
    init_radii : array or None, initial vertex radii (from convex solution)
    """
    rng = np.random.RandomState(seed)
    mesh = NonConvexMesh(subdivisions=subdivisions)
    n_genes = mesh.n_vertices

    # Initialize population
    population = np.zeros((population_size, n_genes))
    if init_radii is not None:
        # Seed from convex solution with perturbations
        for i in range(population_size):
            noise = rng.normal(0, 0.05, n_genes)
            population[i] = np.clip(init_radii + noise, 0.1, 10.0)
    else:
        for i in range(population_size):
            population[i] = 1.0 + rng.normal(0, 0.2, n_genes)
            population[i] = np.clip(population[i], 0.1, 10.0)

    n_elite = max(1, int(population_size * elite_fraction))
    best_fitness = -np.inf
    best_radii = None
    best_chi2 = np.inf

    if verbose:
        print(f"Genetic solver: {n_genes} vertices, pop={population_size}, "
              f"gen={generations}")

    for gen in range(generations):
        # Evaluate fitness
        fitnesses = np.zeros(population_size)
        chi2s = np.zeros(population_size)
        for i in range(population_size):
            fit, chi2, n_data = fitness_function(
                population[i], mesh, lightcurves, pole_lambda, pole_beta,
                period, phi0, t0, lambda_smooth, lambda_regularity)
            fitnesses[i] = fit
            chi2s[i] = chi2

        # Track best
        best_idx = np.argmax(fitnesses)
        if fitnesses[best_idx] > best_fitness:
            best_fitness = fitnesses[best_idx]
            best_radii = population[best_idx].copy()
            best_chi2 = chi2s[best_idx]

        if verbose and (gen + 1) % max(1, generations // 10) == 0:
            print(f"  Gen {gen+1}/{generations}: best_chi2={best_chi2:.2f}, "
                  f"mean_chi2={np.mean(chi2s):.2f}")

        # Selection (tournament)
        sorted_idx = np.argsort(fitnesses)[::-1]

        new_pop = np.zeros_like(population)
        # Elitism
        for i in range(n_elite):
            new_pop[i] = population[sorted_idx[i]]

        # Tournament selection + crossover + mutation
        for i in range(n_elite, population_size):
            # Tournament selection (size 3)
            t1, t2, t3 = rng.choice(population_size, 3, replace=False)
            parent1_idx = max(t1, t2, t3, key=lambda x: fitnesses[x])
            t1, t2, t3 = rng.choice(population_size, 3, replace=False)
            parent2_idx = max(t1, t2, t3, key=lambda x: fitnesses[x])

            parent1 = population[parent1_idx]
            parent2 = population[parent2_idx]

            # Crossover
            if rng.random() < crossover_rate:
                # Uniform crossover
                mask = rng.random(n_genes) < 0.5
                child = np.where(mask, parent1, parent2)
            else:
                child = parent1.copy()

            # Mutation
            mutate_mask = rng.random(n_genes) < mutation_rate
            child[mutate_mask] += rng.normal(0, 0.1, np.sum(mutate_mask))
            child = np.clip(child, 0.1, 10.0)

            new_pop[i] = child

        population = new_pop

    # Final result
    mesh.set_radii(best_radii)
    n_data = sum(len(lc['times']) for lc in lightcurves)

    return {
        'mesh': mesh,
        'radii': best_radii,
        'pole_lambda': pole_lambda,
        'pole_beta': pole_beta,
        'period': period,
        'chi2': best_chi2,
        'chi2_reduced': best_chi2 / max(n_data - n_genes, 1),
        'n_data': n_data,
        'converged': True,
    }
