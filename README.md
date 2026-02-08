# Asteroid Lightcurve Inversion Pipeline

A hybrid lightcurve inversion pipeline for determining 3D shape models and spin states of asteroids from photometric observations. Implements convex inversion (Kaasalainen & Torppa 2001), genetic non-convex refinement (SAGE-inspired, Bartczak & Dudzinski 2018), and sparse data fusion (Durech et al. 2010).

## Features

- **ALCDEF Data Parsing**: Reads and summarizes all 24,643 lightcurve files from ALCDEF archives
- **MPCORB Orbital Elements**: Parses MPC orbital element catalogues for observation geometry computation
- **Convex Inversion**: Spherical harmonics shape parameterization with Levenberg-Marquardt optimization
- **Genetic Non-Convex Solver**: Vertex-based mesh evolution for contact binaries and bifurcated shapes
- **Sparse Data Handling**: Supports individual magnitude measurements from survey photometry
- **Period Search**: Lomb-Scargle periodogram and Phase Dispersion Minimization (PDM)
- **Shape Metrics**: Hausdorff distance and Volumetric IoU for model comparison
- **Target Selection**: Automated candidate identification cross-referencing MPCORB, ALCDEF, and DAMIT

## Installation

```bash
pip install -r requirements.txt
```

Dependencies: numpy, scipy, matplotlib, seaborn, astropy, pandas, PyYAML, Pillow, pytest

## Data Requirements

Place these files in the repository root:
- `ALCDEF_ALL.zip` — ALCDEF photometric database archive
- `MPCORB.DAT.gz` — Minor Planet Center orbital elements catalogue

## Usage

```bash
# Run full pipeline (parse data, select targets, inversion)
python run_pipeline.py --all

# Individual steps
python run_pipeline.py --parse-alcdef    # Parse ALCDEF data
python run_pipeline.py --select          # Select candidate targets
python run_pipeline.py --validate        # Run validation on ground-truth asteroids
python run_pipeline.py --batch           # Batch inversion on 50 candidates
```

### Validation

Run blind tests on asteroids with known shapes (433 Eros, 25143 Itokawa, 216 Kleopatra):

```bash
python -m src.inversion.run_validation
```

### Running Tests

```bash
python -m pytest tests/ -v
```

## Project Structure

```
├── run_pipeline.py              # Main entry point
├── requirements.txt             # Python dependencies
├── sources.bib                  # BibTeX references
├── ALCDEF_ALL.zip               # ALCDEF photometric data (LFS)
├── MPCORB.DAT.gz                # MPC orbital elements (LFS)
├── src/
│   ├── data/
│   │   └── parse_alcdef.py      # ALCDEF archive parser
│   ├── geometry/
│   │   └── ephemeris.py         # MPCORB parser & Keplerian geometry
│   ├── shapes/
│   │   └── convex_model.py      # Convex shape model & forward model
│   ├── inversion/
│   │   ├── period_search.py     # Period determination (LS + PDM)
│   │   ├── convex_solver.py     # Convex inversion solver
│   │   ├── genetic_solver.py    # Genetic non-convex solver
│   │   ├── sparse_solver.py     # Sparse photometry solver
│   │   ├── hybrid_pipeline.py   # Multi-method pipeline
│   │   └── run_validation.py    # Validation test runner
│   ├── metrics/
│   │   └── shape_comparison.py  # Hausdorff distance & IoU
│   └── targets/
│       └── selector.py          # Candidate selection engine
├── data/
│   └── ground_truth/            # Validation shape models (Eros, Itokawa, Kleopatra)
├── results/                     # Output CSV, JSON, reports
│   ├── shapes/                  # Output .obj shape models
│   ├── validation_report.md     # Validation results vs published methods
│   ├── batch_run_log.csv        # Batch inversion progress
│   └── final_candidates.csv     # Ranked new shape models
├── figures/                     # Plots and visualizations
│   └── shapes/                  # Per-asteroid multi-view renders
└── tests/                       # Unit tests (55 tests)
```

## Results Summary

### Validation (Ground-Truth Asteroids)

| Asteroid | Data Points | IoU | Hausdorff (norm.) | Chi-squared |
|----------|-------------|-----|-------------------|-------------|
| 216 Kleopatra | 66 | 0.574 | 0.175 | 1.838 |
| 25143 Itokawa | 33 | 0.364 | 0.574 | 11.927 |
| 433 Eros | 7 | 0.164 | 1.116 | 2.508 |

### Batch Run (50 Candidates)

- 26 out of 50 candidates converged (chi-squared < 5.0)
- All 26 converged models have confidence scores > 0.7
- Convergence rate (52%) consistent with published rates (Durech et al. 2016: 40-60%)

## References

See `sources.bib` for full bibliography. Key references:

- Kaasalainen, M. & Torppa, J. (2001), Icarus 153, 24-36
- Kaasalainen, M., Torppa, J. & Muinonen, K. (2001), Icarus 153, 37-51
- Durech, J. et al. (2010), A&A 513, A46
- Bartczak, P. & Dudzinski, G. (2018), MNRAS 473, 5085-5098
