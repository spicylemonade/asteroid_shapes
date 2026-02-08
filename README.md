# Asteroid Lightcurve Inversion Pipeline

A Python pipeline for deriving 3D shape models, spin axes, and rotation periods of asteroids from photometric lightcurve data. The pipeline ingests ALCDEF photometry observations and MPCORB orbital elements, searches for sidereal rotation periods, and solves the lightcurve inversion problem using convex optimization, genetic algorithms, sparse-data techniques, and hybrid dense+sparse fusion.

## Table of Contents

- [Project Description](#project-description)
- [Installation](#installation)
- [Usage Guide](#usage-guide)
- [Directory Structure](#directory-structure)
- [Dependencies](#dependencies)
- [References](#references)

---

## Project Description

Asteroid lightcurve inversion is the process of reconstructing a three-dimensional shape model and spin state from time-series brightness measurements (lightcurves). As an asteroid rotates, its apparent brightness changes due to varying cross-sectional area and illumination geometry. By fitting a forward model of reflected sunlight to observed lightcurves, the inverse problem yields the asteroid's convex (or non-convex) shape, spin axis orientation, and sidereal rotation period.

This pipeline implements:

- **Convex inversion** (Kaasalainen & Torppa 2001): L-BFGS-B optimization of chi-squared between observed and synthetic lightcurves over spin axis, period, and triangulated convex shape facet areas with smoothness regularization.
- **Genetic algorithm for non-convex shapes** (SAGE-inspired, Bartczak & Dudzinski 2018): Evolutionary optimization with vertex-displacement genomes, mutation operators (vertex perturbation, concavity introduction, local smoothing), uniform crossover, and tournament selection.
- **Sparse data inversion** (Durech et al. 2009, Cellino et al. 2009): Extension for survey-quality photometry with absolute magnitude calibration via the IAU H-G phase function and bootstrap uncertainty estimation.
- **Hybrid dense+sparse fusion** (ADAM-inspired, Viikinkoski et al. 2015): Unified chi-squared combining dense ALCDEF lightcurves with sparse survey photometry under configurable weighting.
- **Adaptive regularization loop**: Recursive self-reinforcement that automatically adjusts regularization weights, mesh resolution, and solver iterations when validation metrics fall below threshold.

The pipeline was validated on three ground-truth asteroids (433 Eros, 216 Kleopatra, 25143 Itokawa) and applied to produce new shape models for 10 candidate near-Earth and main-belt asteroids not previously modeled in the DAMIT database.

---

## Installation

### Prerequisites

- Python 3.8 or later
- pip (Python package manager)

### Install dependencies

```bash
pip install numpy scipy matplotlib
```

Optionally, for enhanced FITS file handling and coordinate transformations:

```bash
pip install astropy
```

### Clone the repository

```bash
git clone <repository-url>
cd repo
```

### Data files

The pipeline requires two data archives in the repository root:

| File | Description | Source |
|------|-------------|--------|
| `ALCDEF_ALL.zip` | ALCDEF photometric lightcurve archive (~135 MB) | [alcdef.org](https://alcdef.org) |
| `MPCORB.DAT.gz` | Minor Planet Center orbital elements (~91 MB) | [minorplanetcenter.net](https://minorplanetcenter.net/iau/MPCORB.html) |

These files are tracked via Git LFS and should be present after cloning.

### Verify installation

```python
from lci_engine import __version__
print(__version__)  # 0.1.0

from lci_engine.parsers import load_alcdef_asteroid
sessions = load_alcdef_asteroid('ALCDEF_ALL.zip', asteroid_number=433)
print(f"Loaded {len(sessions)} sessions for 433 Eros")
```

---

## Usage Guide

### Quick start: invert a single asteroid

```python
import numpy as np
from lci_engine.parsers import load_alcdef_asteroid
from lci_engine.period_search import find_best_period, combine_lightcurve_sessions
from lci_engine.inversion import convex_inversion
from lci_engine.forward_model import save_mesh_obj

# 1. Load ALCDEF lightcurve data
sessions = load_alcdef_asteroid('ALCDEF_ALL.zip', asteroid_number=1943)
print(f"Loaded {len(sessions)} lightcurve sessions")

# 2. Find the rotation period
times, mags, errs = combine_lightcurve_sessions(sessions)
best_period, periods, scores = find_best_period(
    times, mags, errs,
    period_min=2.0,     # hours
    period_max=50.0,    # hours
    timeout_sec=120,
)
print(f"Best period: {best_period:.4f} hours")

# 3. Run convex lightcurve inversion
result = convex_inversion(
    sessions,
    period_init=best_period,
    pole_lambda_init=0.0,         # initial pole longitude (radians)
    pole_beta_init=np.pi / 4,    # initial pole latitude (radians)
    n_facets=200,                 # mesh resolution
    lambda_smooth=0.3,            # regularization weight
    max_iter=100,                 # L-BFGS-B iterations
)

print(f"Period:  {result.period:.4f} h")
print(f"Pole:    lambda={np.degrees(result.pole_lambda):.1f} deg, "
      f"beta={np.degrees(result.pole_beta):.1f} deg")
print(f"RMS:     {result.residual_rms:.4f} mag")

# 4. Save the shape model as Wavefront OBJ
save_mesh_obj(result.vertices, result.faces, 'results/asteroid_1943.obj')
```

### Running the full candidate pipeline

To process all top-10 candidate asteroids at once:

```bash
python scripts/run_candidates.py
```

This script performs period search and convex inversion for each candidate, saves `.obj` meshes to `results/`, and produces a summary JSON and CSV.

### Running blind validation

To reproduce the ground-truth validation on Eros, Kleopatra, and Itokawa:

```bash
python scripts/blind_validation.py
```

### Sparse data experiment

To compare dense versus sparse inversion performance:

```bash
python scripts/sparse_experiment.py
```

### Uncertainty quantification

Bootstrap or jackknife uncertainty estimation:

```bash
python scripts/uncertainty_quantification.py   # bootstrap (50 iterations)
python scripts/quick_uncertainty.py             # jackknife (faster)
```

### Generating 3D shape figures

To render publication-quality 3D visualizations of all shape models:

```bash
python scripts/generate_figures.py
```

Figures are saved as PNG files in the `figures/` directory.

### Using the genetic algorithm for non-convex refinement

```python
from lci_engine.ga_optimizer import ga_optimize

ga_result = ga_optimize(
    sessions,
    period_init=best_period,
    pole_lambda_init=result.pole_lambda,
    pole_beta_init=result.pole_beta,
    convex_vertices=result.vertices,
    convex_faces=result.faces,
    population_size=50,
    n_generations=100,
)
print(f"GA fitness: {ga_result.fitness:.2f}, RMS: {ga_result.residual_rms:.4f}")
```

### Using hybrid dense+sparse fusion

```python
from lci_engine.hybrid_fusion import hybrid_inversion
from lci_engine.sparse_inversion import subsample_to_sparse

sparse_data = subsample_to_sparse(sessions, n_points_per_apparition=50)

hybrid_result = hybrid_inversion(
    sessions,
    sparse_data=sparse_data,
    period_init=best_period,
    w_dense=1.0,
    w_sparse=0.5,
)
```

---

## Directory Structure

```
repo/
|-- README.md                   Project documentation (this file)
|-- sources.bib                 BibTeX references for all consulted papers
|-- ALCDEF_ALL.zip              ALCDEF photometric lightcurve archive (Git LFS)
|-- MPCORB.DAT.gz              MPC orbital elements database (Git LFS)
|-- blind_test_eros.py          Quick standalone Eros inversion test
|
|-- lci_engine/                 Core Python package
|   |-- __init__.py             Package metadata and version
|   |-- parsers.py              ALCDEF and MPCORB data ingestion and parsing
|   |-- forward_model.py        Convex shape representation, Lommel-Seeliger
|   |                           scattering, synthetic lightcurve generation
|   |-- period_search.py        PDM, chi-squared periodogram, Lomb-Scargle
|   |                           period search with alias resolution
|   |-- inversion.py            Convex lightcurve inversion solver
|   |                           (Kaasalainen-Torppa L-BFGS-B method)
|   |-- ga_optimizer.py         Genetic algorithm for non-convex shape
|   |                           optimization (SAGE-inspired)
|   |-- sparse_inversion.py     Sparse survey-data inversion with H-G
|   |                           calibration and bootstrap uncertainty
|   |-- hybrid_fusion.py        Hybrid dense+sparse data fusion (ADAM-inspired)
|   |-- adaptive_loop.py        Adaptive regularization and recursive
|   |                           optimization loop
|   |-- validation.py           Hausdorff distance, volumetric IoU, pole/period
|   |                           error metrics, and validation reporting
|
|-- scripts/                    Experiment and analysis scripts
|   |-- run_candidates.py       Run inversion on top-10 candidate asteroids
|   |-- blind_validation.py     Blind validation on ground-truth asteroids
|   |-- sparse_experiment.py    Dense vs. sparse inversion comparison
|   |-- uncertainty_quantification.py  Bootstrap uncertainty estimation
|   |-- quick_uncertainty.py    Jackknife uncertainty estimation (faster)
|   |-- generate_figures.py     3D shape rendering and figure generation
|
|-- results/                    Output directory for shape models and reports
|   |-- *.obj                   Wavefront OBJ shape model files
|   |-- *.json                  Structured result and validation reports
|   |-- *.csv                   Summary tables
|
|-- figures/                    Output directory for rendered shape figures
|   |-- *.png                   3D shape visualizations
|
|-- .archivara/                 Project orchestration metadata
```

---

## Dependencies

### Required

| Package | Version | Purpose |
|---------|---------|---------|
| **numpy** | >= 1.20 | Array operations, linear algebra, mesh computation |
| **scipy** | >= 1.7 | L-BFGS-B optimization, Lomb-Scargle periodogram, spatial algorithms (cKDTree, Delaunay, ConvexHull) |
| **matplotlib** | >= 3.4 | 3D shape visualization and figure generation |

### Optional

| Package | Version | Purpose |
|---------|---------|---------|
| **astropy** | >= 5.0 | FITS file handling, coordinate frame transformations, time utilities |

### Python version

- Python >= 3.8

All core inversion logic is implemented from scratch using only NumPy and SciPy for numerical operations. No external lightcurve inversion libraries are required.

---

## References

Key papers underpinning this pipeline (see `sources.bib` for full BibTeX entries):

- Kaasalainen, M. & Torppa, J. (2001). Optimization methods for asteroid lightcurve inversion. *Icarus*, 153, 24-36.
- Kaasalainen, M., Torppa, J. & Muinonen, K. (2001). Optimization methods for asteroid lightcurve inversion II. *Icarus*, 153, 37-51.
- Bartczak, P. & Dudzinski, G. (2018). Shaping asteroids with genetic evolution (SAGE). *MNRAS*, 473, 5050-5065.
- Durech, J. et al. (2009). Asteroid models from combined sparse and dense photometric data. *A&A*, 493, 291-297.
- Cellino, A. et al. (2009). Asteroid photometric inversion from sparse data. *Icarus*, 204, 497-504.
- Viikinkoski, M. et al. (2015). ADAM: All-Data Asteroid Modeling. *A&A*, 576, A8.
- Hapke, B. (1993). *Theory of Reflectance and Emittance Spectroscopy*. Cambridge University Press.
- Bowell, E. et al. (1989). Application of photometric models to asteroids. In *Asteroids II*, University of Arizona Press.
