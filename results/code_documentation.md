# Asteroid Lightcurve Inversion Pipeline -- Code Documentation

## Table of Contents

1. [Source Modules and Their Roles](#1-source-modules-and-their-roles)
2. [Dependencies](#2-dependencies)
3. [Running the Full Pipeline End-to-End](#3-running-the-full-pipeline-end-to-end)
4. [How to Add New Data Sources](#4-how-to-add-new-data-sources)
5. [Configuration Parameters and Their Effects](#5-configuration-parameters-and-their-effects)
6. [Pipeline Flowchart](#6-pipeline-flowchart)

---

## 1. Source Modules and Their Roles

All source code resides in the `src/` directory. The modules are organized by
functional layer: data parsing, data ingestion, scientific solvers, utilities,
target selection, and pipeline orchestration.

### Data Parsing Layer

#### `src/parse_alcdef.py`
**Role:** Parse ALCDEF_ALL.zip to build a comprehensive catalog of asteroid
lightcurve data.

- Reads the full ALCDEF archive (zip format) containing thousands of individual
  lightcurve text files.
- Parses the ALCDEF metadata/data block format: `STARTMETADATA`/`ENDMETADATA`
  sections followed by pipe-delimited `DATA=` lines containing JD timestamps,
  magnitudes, and uncertainties.
- Aggregates per-asteroid statistics: number of lightcurves, total data points,
  date ranges, observer codes, filter bands, and MPC site codes.
- **Output:** `results/alcdef_catalog.csv`

#### `src/parse_mpcorb.py`
**Role:** Parse MPCORB.DAT.gz to extract orbital elements for all minor planets,
flag NEOs, and estimate diameters from absolute magnitude H.

- Decompresses and parses the MPC Orbit fixed-width format file.
- Extracts orbital elements: epoch, mean anomaly, argument of perihelion,
  longitude of ascending node, inclination, eccentricity, semimajor axis,
  absolute magnitude H, and slope parameter G.
- Classifies orbits as NEO (perihelion q < 1.3 AU), MBA, or TNO.
- Estimates diameters using D(km) = 1329 / sqrt(albedo) * 10^(-H/5) with a
  default albedo of 0.15.
- **Output:** `results/mpcorb_parsed.csv`

### Data Ingestion Layer

#### `src/data_ingest.py`
**Role:** ALCDEF data ingestion and lightcurve preprocessing module.

- Reads ALCDEF text files from the extracted zip archive.
- Parses JD timestamps, magnitudes, and uncertainties from each lightcurve block.
- Loads orbital elements from MPCORB for a given asteroid.
- Computes viewing geometry for each observation: phase angle, aspect angle, and
  solar elongation, using Keplerian orbital position calculations.
- Computes approximate ecliptic positions for both the asteroid and Earth at each
  observation epoch, then derives sun-asteroid-observer geometry vectors.
- **Key functions:** `parse_alcdef_file()`, `load_mpcorb_for_asteroid()`,
  `compute_viewing_geometry()`, `preprocess_asteroid()`,
  `orbital_position_ecliptic()`, `earth_position_ecliptic()`
- **References:** Kaasalainen et al. (2001), Durech et al. (2009)

### Scientific Solvers

#### `src/period_search.py`
**Role:** Period search engine using Lomb-Scargle periodogram and Phase Dispersion
Minimization (PDM).

- Implements Lomb-Scargle periodogram via `scipy.signal.lombscargle`, searching
  for the rotation period in a configurable range (default 0.5--100 hours).
- Accounts for the double-peaked nature of asteroid lightcurves by computing
  P_rot = 2 * P_LS.
- Implements Phase Dispersion Minimization (PDM) as an independent period-finding
  method (Stellingwerf 1978).
- Provides a combined search function that merges results from both methods and
  ranks candidate periods by confidence score.
- **Key functions:** `lomb_scargle_search()`, `pdm_search()`,
  `combined_period_search()`, `merge_sessions_to_relative()`
- **References:** Lomb (1976), Scargle (1982), Stellingwerf (1978)

#### `src/convex_inversion.py`
**Role:** Convex lightcurve inversion solver using the Kaasalainen-Torppa method.

- Implements a two-stage approach:
  - **Stage 1:** Triaxial ellipsoid fit with 3 shape parameters (a, b/a, c/a)
    plus spin axis grid search over ecliptic (lambda, beta).
  - **Stage 2:** Vertex-level deformation refinement on a 162-vertex icosphere
    mesh, optimizing per-vertex radii to minimize chi-squared.
- Forward brightness model uses the Lommel-Seeliger + Lambert scattering law.
- Spin axis conversion from ecliptic longitude/latitude to unit vector.
- Uses `scipy.optimize.minimize` (L-BFGS-B) and `differential_evolution` for
  optimization.
- Enforces convexity via the `scipy.spatial.ConvexHull`.
- **Key classes/functions:** `ConvexSolver`, `grid_search_spin()`,
  `subsample_blocks()`, `spin_axis_from_ecliptic()`
- **References:** Kaasalainen & Torppa (2001), Kaasalainen, Torppa &
  Muinonen (2001)

#### `src/ga_solver.py`
**Role:** Genetic algorithm non-convex shape solver (SAGE-inspired) using the
convex solution as a seed.

- Encodes vertex radial displacements from the convex seed mesh as the genome,
  allowing both concavities (inward deformation) and protrusions (outward
  deformation).
- Fitness function evaluates chi-squared of observed vs. modeled lightcurves
  using the self-shadowing forward photometric model from `scattering.py`.
- GA operators: tournament selection, two-point crossover, Gaussian mutation,
  and elitism.
- Default parameters: population size 50, 100 generations.
- Runs on ALL targets (no convex-first gatekeeper); the convex solution from
  `convex_inversion.py` serves only as the initial seed.
- **Key classes:** `GASolver`
- **References:** Bartczak & Dudzinski (2018)

#### `src/scattering.py`
**Role:** Forward photometric model with Lommel-Seeliger + Lambert scattering law
and self-shadowing ray-tracing.

- Implements the combined Lommel-Seeliger + Lambert scattering law:
  S = f(alpha) * [c_ls * mu_0 / (mu_0 + mu) + c_l * mu_0]
- BVH (Bounding Volume Hierarchy) accelerated ray-mesh intersection for
  self-shadowing: for each visible facet, casts rays toward the Sun direction
  and checks for intersections with other facets.
- Facets in shadow contribute zero direct illumination.
- Configurable scattering coefficients: `c_ls` (Lommel-Seeliger weight, default
  0.5), `c_l` (Lambert weight, default 0.1).
- Performance: >2500 shape evaluations per minute on a single CPU core.
- **Key functions:** `lommel_seeliger_lambert()`, `compute_brightness()`,
  shadow-casting routines
- **References:** Kaasalainen & Torppa (2001), Durech & Kaasalainen (2003),
  Hapke (1993, 2012)

#### `src/sparse_inversion.py`
**Role:** Sparse data inversion module for survey-like sparse photometry
(Gaia/ZTF/Pan-STARRS style).

- Handles sparse photometry with fewer than 100 data points spread across
  multiple apparitions (minimum 3 required).
- Creates sparse subsets from dense lightcurves to simulate survey-like sparse
  sampling.
- Uses coarse spin-axis grid search combined with convex + GA refinement.
- Includes a data fusion mode (`data_fusion_inversion`) that combines dense
  lightcurves (for shape detail) with sparse absolute photometry (for pole and
  albedo constraints) in a single unified inversion.
- **Key functions:** `create_sparse_subset()`, `sparse_inversion()`,
  `data_fusion_inversion()`
- **References:** Durech et al. (2010), Cellino et al. (2015), Viikinkoski
  et al. (2015)

### Utilities

#### `src/mesh_utils.py`
**Role:** Mesh utilities including generation, I/O, and comparison metrics.

- Generates icosphere meshes via recursive subdivision of an icosahedron
  (`create_sphere_mesh()`).
- Reads and writes Wavefront OBJ mesh files (`load_obj()`, `save_obj()`).
- Computes facet properties: normals, areas, and centroids
  (`compute_facet_properties()`).
- Implements validation metrics:
  - **Hausdorff distance:** symmetric, normalized by mesh diameter.
  - **Volumetric IoU:** intersection over union via voxelization at configurable
    resolution.
- **Key functions:** `create_sphere_mesh()`, `compute_facet_properties()`,
  `save_obj()`, `load_obj()`, `hausdorff_distance()`, `volumetric_iou()`

### Target Selection

#### `src/target_selection.py`
**Role:** Identify top 50 unmodeled asteroids meeting all priority criteria.

- Cross-references the ALCDEF catalog, MPCORB orbital database, and DAMIT
  database to select targets matching priority criteria:
  - P1: NEO flag OR estimated diameter > 100 km
  - P2: Period quality U >= 2 in ALCDEF metadata
  - P3: NOT in DAMIT (no existing published shape model)
  - P4: >= 20 dense lightcurves OR >= 100 sparse points across >= 3 apparitions
- Computes a composite priority score and ranks candidates.
- **Output:** `results/target_candidates.csv`

### Pipeline Orchestration

#### `src/run_pipeline.py`
**Role:** Full pipeline runner that executes period search, convex inversion, and
GA non-convex refinement on all targets.

- Processes each target in the candidate list sequentially.
- Pre-caches data loading: opens the ALCDEF zip file and MPCORB orbital cache
  once, then reuses for all targets.
- Builds a namelist index mapping asteroid numbers to zip file entries for fast
  lookup.
- For each asteroid: loads lightcurve data, computes viewing geometry,
  runs period search, runs convex inversion as seed, runs GA non-convex solver,
  and saves results.
- **Output:** OBJ shape files in `results/shapes/`, spin vectors in
  `results/spin_vectors.csv`, pipeline metrics in `results/pipeline_results.json`
- **Key functions:** `load_alcdef_fast()`, `build_namelist_index()`,
  `preprocess_fast()`

#### `src/run_validation.py`
**Role:** Blind validation runner that inverts ALCDEF data for ground truth
asteroids and compares results against DAMIT/spacecraft shape models.

- Runs the full inversion pipeline on asteroids with known shape models.
- Computes angular distance between derived and known pole directions.
- Computes Hausdorff distance and volumetric IoU against ground truth shapes.
- **Output:** `results/validation_report.json`
- **Key functions:** `angular_distance_deg()`

---

## 2. Dependencies

### Required Python Packages

| Package | Version | Purpose |
|---------|---------|---------|
| **numpy** | >= 1.20 | Core numerical arrays, linear algebra, random sampling. Used throughout all modules for vectorized computations on lightcurve data, mesh vertices, and photometric models. |
| **scipy** | >= 1.7 | `scipy.signal.lombscargle` for period search; `scipy.optimize.minimize` and `differential_evolution` for convex inversion; `scipy.spatial.ConvexHull` for convexity enforcement and mesh operations. |
| **matplotlib** | >= 3.4 | 3D shape rendering (`mpl_toolkits.mplot3d`) for publication-quality figures at 300 DPI. Used in figure generation scripts. |

### Standard Library Dependencies

- `zipfile` -- reading ALCDEF_ALL.zip
- `gzip` -- reading MPCORB.DAT.gz
- `csv` -- reading and writing CSV catalogs
- `json` -- reading and writing pipeline results and configuration
- `os`, `sys` -- path management and module imports
- `math` -- diameter estimation formula
- `time` -- pipeline timing
- `traceback` -- error handling and logging
- `gc` -- garbage collection for memory management during batch processing
- `copy` -- deep copying of GA genomes
- `collections.defaultdict` -- aggregating ALCDEF catalog statistics

### Installation

```bash
pip install numpy scipy matplotlib
```

No additional compiled dependencies or system libraries are required. The
pipeline runs on standard CPython 3.8+.

---

## 3. Running the Full Pipeline End-to-End

### Prerequisites

1. Place the data files in the repository root:
   - `ALCDEF_ALL.zip` -- ALCDEF photometric archive
   - `MPCORB.DAT.gz` -- MPC orbital elements file

2. Ensure the `results/` and `results/shapes/` directories exist:
   ```bash
   mkdir -p results/shapes results/ground_truth figures
   ```

### Step-by-Step Execution

#### Step 1: Parse raw data files
```bash
cd /path/to/repo
python src/parse_alcdef.py
python src/parse_mpcorb.py
```
This produces:
- `results/alcdef_catalog.csv` -- catalog of all asteroids and their lightcurve
  data statistics
- `results/mpcorb_parsed.csv` -- orbital elements with NEO flags and diameter
  estimates

#### Step 2: Select target candidates
```bash
python src/target_selection.py
```
This cross-references the ALCDEF catalog, MPCORB database, and DAMIT to produce:
- `results/target_candidates.csv` -- ranked list of 50 candidate asteroids

#### Step 3: Run the full inversion pipeline
```bash
python src/run_pipeline.py
```
This is the main computational step. For each of the 50 targets, it:
1. Loads and preprocesses lightcurve data from ALCDEF
2. Computes viewing geometry from MPCORB orbital elements
3. Runs Lomb-Scargle + PDM period search
4. Runs convex inversion (Kaasalainen-Torppa) as seed
5. Runs GA non-convex solver with self-shadowing

**Output:**
- `results/shapes/{asteroid_id}_ga.obj` or `{asteroid_id}_convex.obj` -- 3D
  shape model meshes in Wavefront OBJ format
- `results/spin_vectors.csv` -- spin axis and period for each target
- `results/pipeline_results.json` -- detailed metrics for each target

Expected runtime: approximately 35 minutes for 50 targets on a modern single-core
CPU.

#### Step 4 (Optional): Run blind validation
```bash
python src/run_validation.py
```
Validates the pipeline against ground truth shape models stored in
`results/ground_truth/`. Produces `results/validation_report.json`.

---

## 4. How to Add New Data Sources

### Adding New Lightcurve Data

The pipeline is designed to work with ALCDEF-format photometric data. To add new
data sources:

1. **Convert to ALCDEF format:** New lightcurve data must be converted to the
   ALCDEF pipe-delimited text format with `STARTMETADATA`/`ENDMETADATA` blocks
   and `DATA=JD|mag|err` lines. The required metadata fields are:
   - `OBJECTNUMBER` -- asteroid number
   - `OBJECTNAME` -- asteroid name
   - `SESSIONDATE` -- observation date
   - `FILTER` -- photometric filter band

2. **Add to the archive:** Either append the new `.txt` files to `ALCDEF_ALL.zip`
   or place them in a separate directory and modify `src/data_ingest.py` to read
   from multiple sources.

3. **Re-run the catalog builder:**
   ```bash
   python src/parse_alcdef.py
   ```

### Adding New Orbital Element Sources

To use orbital elements from a source other than MPCORB:

1. Create a new parser module (e.g., `src/parse_jpl_horizons.py`) that produces
   a CSV with the same column schema as `results/mpcorb_parsed.csv`:
   - `number`, `name`, `epoch`, `M_deg`, `omega_deg`, `Omega_deg`, `i_deg`,
     `e`, `a_AU`, `H_mag`, `G_slope`, `neo_flag`, `est_diameter_km`

2. Update `src/data_ingest.py` function `load_mpcorb_for_asteroid()` to read
   from the new source file, or replace the file at `results/mpcorb_parsed.csv`.

### Adding Survey Sparse Photometry

For sparse photometric data from surveys (Gaia, ZTF, LSST):

1. Format the data as lightcurve blocks with one or two points per session.
2. Use `src/sparse_inversion.py` directly with the `sparse_inversion()` or
   `data_fusion_inversion()` functions, which accept preprocessed block lists.
3. Ensure that absolute (calibrated) magnitudes are used rather than relative
   magnitudes for sparse data, as the absolute photometry constrains the pole
   orientation and albedo.

---

## 5. Configuration Parameters and Their Effects

### Period Search (`src/period_search.py`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `period_min` | 0.5 hours | Minimum rotation period to search. Reducing this captures very fast rotators but increases computation. |
| `period_max` | 100.0 hours | Maximum rotation period. Increasing this captures slow rotators at the cost of reduced frequency resolution. |
| `n_periods` | 50,000 | Number of trial periods in the grid. Higher values improve period resolution but increase runtime linearly. |
| `top_n` | 10 | Number of top candidate periods returned. |

### Convex Inversion (`src/convex_inversion.py`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `n_subdivisions` | 3 | Icosphere subdivision level. Level 3 produces 162 vertices. Higher values give finer shape detail but increase computation quadratically. |
| Spin grid step | 30 degrees | Step size for lambda and beta grid search. Finer grids (e.g., 5 or 10 degrees) improve pole accuracy but increase computation as O(1/step^2). |
| `b_a`, `c_a` range | [0.3, 1.0] | Bounds on ellipsoid axis ratios for Stage 1 fit. |

### GA Non-Convex Solver (`src/ga_solver.py`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `pop_size` | 50 | GA population size. Larger populations explore more of the shape space but increase runtime linearly. |
| `n_gen` | 100 | Number of GA generations. More generations allow finer convergence but with diminishing returns after ~200. |
| Mutation rate | Gaussian | Gaussian perturbation applied to vertex radii. Larger sigma explores more aggressively but may overshoot. |
| Elitism | Yes | Best individuals survive to next generation unchanged, preventing regression. |
| Crossover | Two-point | Two-point crossover of vertex-radius genomes between parent individuals. |

### Scattering Model (`src/scattering.py`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `c_ls` | 0.5 | Lommel-Seeliger scattering weight. Higher values increase limb-darkening effect. |
| `c_l` | 0.1 | Lambert scattering weight. Higher values increase uniform (Lambertian) reflection contribution. |
| Self-shadowing | Enabled | BVH ray-tracing for shadow detection. Critical for non-convex shapes. Disabling speeds up evaluation but produces incorrect brightness for concave shapes. |

### Sparse Inversion (`src/sparse_inversion.py`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| `max_points` | 80 | Maximum number of sparse data points used. |
| `min_apparitions` | 3 | Minimum number of distinct observing apparitions (>120-day gaps) required. Fewer apparitions yield degenerate pole solutions. |

### Pipeline Runner (`src/run_pipeline.py`)

| Parameter | Default | Effect |
|-----------|---------|--------|
| Random seed | 42 | `np.random.seed(42)` for reproducibility of GA and period search stochastic elements. |
| Data subsampling | 50 points per block | Maximum number of data points per lightcurve block used in inversion. Reduces runtime while preserving lightcurve shape. |

---

## 6. Pipeline Flowchart

```
                    +=====================+
                    |   RAW DATA FILES    |
                    |  ALCDEF_ALL.zip     |
                    |  MPCORB.DAT.gz     |
                    +=====================+
                             |
              +--------------+--------------+
              |                             |
              v                             v
   +--------------------+       +--------------------+
   | src/parse_alcdef.py|       |src/parse_mpcorb.py |
   | Parse lightcurve   |       | Parse orbital      |
   | archive            |       | elements, flag     |
   |                    |       | NEOs, estimate     |
   | Output:            |       | diameters          |
   | alcdef_catalog.csv |       | Output:            |
   +--------------------+       | mpcorb_parsed.csv  |
              |                 +--------------------+
              |                             |
              +-------------+---------------+
                            |
                            v
                 +------------------------+
                 | src/target_selection.py |
                 | Cross-reference ALCDEF |
                 | + MPCORB + DAMIT       |
                 | Apply priority criteria |
                 | (NEO, diameter, data   |
                 |  quality, novelty)     |
                 |                        |
                 | Output:               |
                 | target_candidates.csv  |
                 +------------------------+
                            |
                            v
                 +------------------------+
                 |  src/run_pipeline.py   |
                 |  (Main orchestrator)   |
                 +------------------------+
                            |
               For each of 50 targets:
                            |
                            v
                 +------------------------+
                 |  src/data_ingest.py    |
                 |  Load ALCDEF data      |
                 |  Load MPCORB orbits    |
                 |  Compute viewing       |
                 |  geometry              |
                 +------------------------+
                            |
                            v
                 +------------------------+
                 | src/period_search.py   |
                 | Lomb-Scargle + PDM     |
                 | combined search        |
                 |                        |
                 | Output: top candidate  |
                 | rotation periods       |
                 +------------------------+
                            |
                            v
                 +------------------------+
                 |src/convex_inversion.py |
                 | Stage 1: Ellipsoid fit |
                 |   + spin grid search   |
                 | Stage 2: Vertex-level  |
                 |   mesh refinement      |
                 |                        |
                 | Output: convex seed    |
                 | shape (.obj) + spin    |
                 +------------------------+
                            |
                            v
                 +------------------------+
                 |   src/ga_solver.py     |
                 | GA non-convex solver   |
                 | Uses convex as seed    |
                 | Self-shadowing fitness |
                 |   (via scattering.py)  |
                 |                        |
                 | Output: final shape    |
                 | (.obj) + chi-squared   |
                 +----------+-------------+
                            |
              +-------------+-------------+
              |             |             |
              v             v             v
    +-----------+  +-----------+  +----------------+
    | shapes/   |  | spin_     |  | pipeline_      |
    | {id}.obj  |  | vectors   |  | results.json   |
    | 3D mesh   |  | .csv      |  | Full metrics   |
    +-----------+  +-----------+  +----------------+


    Supporting modules used during inversion:

    +------------------------+       +------------------------+
    |   src/scattering.py    |       |   src/mesh_utils.py    |
    | Lommel-Seeliger +      |       | Icosphere generation   |
    | Lambert scattering     |       | OBJ I/O                |
    | BVH self-shadowing     |       | Hausdorff distance     |
    | ray-tracing            |       | Volumetric IoU         |
    +------------------------+       +------------------------+

    +------------------------+
    |src/sparse_inversion.py |
    | Sparse data handling   |
    | Dense+sparse fusion    |
    +------------------------+


    Validation path (optional):

    +------------------------+       +------------------------+
    | src/run_validation.py  | ----> | results/               |
    | Blind test on ground   |       | validation_report.json |
    | truth asteroids        |       +------------------------+
    +------------------------+
```

---

## Module Dependency Graph

```
run_pipeline.py
  |-- data_ingest.py
  |     |-- (reads ALCDEF_ALL.zip)
  |     |-- (reads MPCORB.DAT.gz)
  |-- period_search.py
  |     |-- numpy
  |     |-- scipy.signal.lombscargle
  |-- convex_inversion.py
  |     |-- numpy
  |     |-- scipy.optimize
  |     |-- scipy.spatial.ConvexHull
  |     |-- mesh_utils.py
  |-- ga_solver.py
  |     |-- numpy
  |     |-- mesh_utils.py
  |     |-- scattering.py
  |-- mesh_utils.py
        |-- numpy
        |-- scipy.spatial.ConvexHull

run_validation.py
  |-- data_ingest.py
  |-- period_search.py
  |-- convex_inversion.py
  |-- ga_solver.py
  |-- mesh_utils.py

sparse_inversion.py
  |-- convex_inversion.py
  |-- numpy
```

---

*Generated for the Asteroid Lightcurve Inversion Pipeline, 2026-02-09*
