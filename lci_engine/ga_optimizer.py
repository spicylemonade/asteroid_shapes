"""
Genetic/Evolutionary Algorithm for non-convex shape optimization (SAGE-inspired).

Implements evolutionary shape optimization where the genome encodes vertex
displacements from a convex hull, with mutation, crossover, and selection
operators to explore non-convex shapes.

References:
    Bartczak & Dudzinski 2018 (Bartczak2018 in sources.bib)
"""

import numpy as np
from typing import Tuple, List, Optional
from dataclasses import dataclass

from .forward_model import (
    compute_facet_properties, ecliptic_to_body_frame, lommel_seeliger,
    brightness_to_magnitude, create_triaxial_ellipsoid
)
from .inversion import prepare_observations, _compute_model_lightcurves


@dataclass
class GAResult:
    """Result from genetic algorithm shape optimization."""
    vertices: np.ndarray
    faces: np.ndarray
    pole_lambda: float
    pole_beta: float
    period: float
    fitness: float
    generation: int
    chi2: float
    residual_rms: float


class Individual:
    """A single shape solution in the GA population."""

    def __init__(self, displacements: np.ndarray, pole_lambda: float,
                 pole_beta: float, period: float):
        self.displacements = displacements.copy()
        self.pole_lambda = pole_lambda
        self.pole_beta = pole_beta
        self.period = period
        self.fitness = -np.inf

    def copy(self):
        ind = Individual(self.displacements, self.pole_lambda,
                         self.pole_beta, self.period)
        ind.fitness = self.fitness
        return ind


def _apply_displacements(base_vertices: np.ndarray,
                          displacements: np.ndarray) -> np.ndarray:
    """Apply radial displacements to base mesh vertices."""
    norms = np.linalg.norm(base_vertices, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    directions = base_vertices / norms
    new_radii = norms.ravel() * (1.0 + displacements)
    new_radii = np.maximum(new_radii, 0.1)
    return directions * new_radii[:, np.newaxis]


def _evaluate_fitness(individual: Individual,
                       base_vertices: np.ndarray,
                       faces: np.ndarray,
                       observations: List[dict],
                       observed_mags: List[np.ndarray],
                       observed_errs: List[np.ndarray],
                       epoch_jd: float,
                       lambda_smooth: float = 0.01) -> float:
    """Compute fitness (negative chi-squared) for an individual."""
    vertices = _apply_displacements(base_vertices, individual.displacements)
    normals, areas, _ = compute_facet_properties(vertices, faces)

    # Compute model lightcurves
    model_mags = _compute_model_lightcurves(
        areas, normals, base_vertices, faces, observations,
        individual.pole_lambda, individual.pole_beta,
        individual.period, epoch_jd
    )

    # Chi-squared
    chi2 = 0.0
    for obs_m, mod_m, err in zip(observed_mags, model_mags, observed_errs):
        if len(obs_m) < 2:
            continue
        offset = np.mean(obs_m) - np.mean(mod_m)
        weights = 1.0 / np.maximum(err, 0.001) ** 2
        chi2 += np.sum(weights * (obs_m - offset - mod_m) ** 2)

    # Smoothness regularization on displacements
    disp = individual.displacements
    smooth_penalty = lambda_smooth * np.sum(np.diff(disp) ** 2)

    return -(chi2 + smooth_penalty)


def _mutate(individual: Individual, rng: np.random.RandomState,
            mutation_rate: float = 0.1, mutation_scale: float = 0.05) -> Individual:
    """Apply mutation operators to an individual."""
    child = individual.copy()
    n = len(child.displacements)

    # Vertex perturbation: randomly perturb a fraction of vertices
    mask = rng.random(n) < mutation_rate
    child.displacements[mask] += mutation_scale * rng.randn(np.sum(mask))

    # Concavity introduction: occasionally push vertices inward
    if rng.random() < 0.05:
        idx = rng.randint(0, n)
        radius = max(1, n // 20)
        lo = max(0, idx - radius)
        hi = min(n, idx + radius)
        child.displacements[lo:hi] -= 0.1 * rng.random()

    # Local smoothing: average some vertices with neighbors
    if rng.random() < 0.1:
        smoothed = child.displacements.copy()
        smoothed[1:-1] = 0.5 * child.displacements[1:-1] + \
                         0.25 * child.displacements[:-2] + \
                         0.25 * child.displacements[2:]
        child.displacements = smoothed

    # Small pole perturbation
    child.pole_lambda += 0.01 * rng.randn()
    child.pole_beta += 0.01 * rng.randn()
    child.pole_beta = np.clip(child.pole_beta, -np.pi / 2, np.pi / 2)

    # Small period perturbation
    child.period += 0.001 * rng.randn()
    child.period = max(1.0, child.period)

    return child


def _crossover(parent_a: Individual, parent_b: Individual,
               rng: np.random.RandomState) -> Individual:
    """Uniform crossover between two parents."""
    n = len(parent_a.displacements)
    mask = rng.random(n) > 0.5
    child_disp = np.where(mask, parent_a.displacements, parent_b.displacements)

    # Average spin parameters
    lam = (parent_a.pole_lambda + parent_b.pole_lambda) / 2
    bet = (parent_a.pole_beta + parent_b.pole_beta) / 2
    per = (parent_a.period + parent_b.period) / 2

    child = Individual(child_disp, lam, bet, per)
    return child


def _tournament_select(population: List[Individual],
                        rng: np.random.RandomState, k: int = 3) -> Individual:
    """Tournament selection: pick best of k random individuals."""
    candidates = rng.choice(len(population), size=min(k, len(population)), replace=False)
    best = max(candidates, key=lambda i: population[i].fitness)
    return population[best].copy()


def ga_optimize(sessions,
                period_init: float,
                pole_lambda_init: float,
                pole_beta_init: float,
                convex_vertices: np.ndarray = None,
                convex_faces: np.ndarray = None,
                population_size: int = 50,
                n_generations: int = 100,
                mutation_rate: float = 0.15,
                mutation_scale: float = 0.05,
                lambda_smooth: float = 0.01,
                seed: int = 42,
                verbose: bool = True) -> GAResult:
    """Run genetic algorithm for non-convex shape optimization.

    Args:
        sessions: List of LightcurveSession objects
        period_init: Initial period (hours)
        pole_lambda_init, pole_beta_init: Initial pole (radians)
        convex_vertices, convex_faces: Initial convex shape (optional)
        population_size: GA population size
        n_generations: Number of generations
        mutation_rate: Fraction of vertices to mutate
        mutation_scale: Standard deviation of mutations
        lambda_smooth: Smoothness regularization weight
        seed: Random seed

    Returns:
        GAResult with optimized non-convex shape
    """
    rng = np.random.RandomState(seed)

    # Prepare observations
    observations, observed_mags, observed_errs = prepare_observations(sessions)
    if not observations:
        raise ValueError("No valid observations")

    epoch_jd = observations[0]['times'][0]

    # Base mesh
    if convex_vertices is None or convex_faces is None:
        n_lat, n_lon = 12, 24
        convex_vertices, convex_faces = create_triaxial_ellipsoid(1, 1, 1, n_lat, n_lon)

    n_vertices = len(convex_vertices)

    # Initialize population
    population = []
    for _ in range(population_size):
        disp = 0.05 * rng.randn(n_vertices)
        lam = pole_lambda_init + 0.1 * rng.randn()
        bet = pole_beta_init + 0.1 * rng.randn()
        bet = np.clip(bet, -np.pi / 2, np.pi / 2)
        per = period_init + 0.01 * rng.randn()
        population.append(Individual(disp, lam, bet, per))

    # Add initial convex solution as one member
    population[0] = Individual(
        np.zeros(n_vertices), pole_lambda_init, pole_beta_init, period_init
    )

    # Evaluate initial population
    for ind in population:
        ind.fitness = _evaluate_fitness(
            ind, convex_vertices, convex_faces, observations,
            observed_mags, observed_errs, epoch_jd, lambda_smooth
        )

    best_ever = max(population, key=lambda x: x.fitness).copy()

    # Evolution loop
    for gen in range(n_generations):
        new_population = []

        # Elitism: keep top 10%
        sorted_pop = sorted(population, key=lambda x: x.fitness, reverse=True)
        n_elite = max(2, population_size // 10)
        new_population.extend([ind.copy() for ind in sorted_pop[:n_elite]])

        # Generate rest through crossover + mutation
        while len(new_population) < population_size:
            parent_a = _tournament_select(population, rng)
            parent_b = _tournament_select(population, rng)
            child = _crossover(parent_a, parent_b, rng)
            child = _mutate(child, rng, mutation_rate, mutation_scale)
            child.fitness = _evaluate_fitness(
                child, convex_vertices, convex_faces, observations,
                observed_mags, observed_errs, epoch_jd, lambda_smooth
            )
            new_population.append(child)

        population = new_population

        # Track best
        gen_best = max(population, key=lambda x: x.fitness)
        if gen_best.fitness > best_ever.fitness:
            best_ever = gen_best.copy()

        if verbose and (gen + 1) % 10 == 0:
            print(f"  Gen {gen+1}/{n_generations}: best_fitness={best_ever.fitness:.2f}")

    # Build final mesh
    final_vertices = _apply_displacements(convex_vertices, best_ever.displacements)
    normals, areas, _ = compute_facet_properties(final_vertices, convex_faces)

    # Compute final residual RMS
    model_mags = _compute_model_lightcurves(
        areas, normals, convex_vertices, convex_faces, observations,
        best_ever.pole_lambda, best_ever.pole_beta, best_ever.period, epoch_jd
    )

    residuals = []
    for obs_m, mod_m in zip(observed_mags, model_mags):
        offset = np.mean(obs_m) - np.mean(mod_m)
        residuals.extend(obs_m - mod_m - offset)

    residual_rms = np.sqrt(np.mean(np.array(residuals) ** 2))

    return GAResult(
        vertices=final_vertices,
        faces=convex_faces,
        pole_lambda=best_ever.pole_lambda,
        pole_beta=best_ever.pole_beta,
        period=best_ever.period,
        fitness=best_ever.fitness,
        generation=n_generations,
        chi2=-best_ever.fitness,
        residual_rms=residual_rms,
    )
