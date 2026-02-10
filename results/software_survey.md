# Survey of Asteroid Lightcurve Inversion Software

## 1. convexinv / conjgradinv (Kaasalainen & Durech)

| Property | Details |
|---|---|
| **Language** | Originally Fortran (Kaasalainen); converted to C by Durech. Repository is 87.4% C, 11.8% Fortran. |
| **License** | CC Attribution 4.0 International |
| **Algorithm** | Convex shape inversion using Levenberg-Marquardt optimization. `convexinv` parameterizes the shape with spherical harmonics and optimizes all parameters (shape, spin axis, period, scattering). `conjgradinv` uses conjugate-gradient optimization of facet areas directly and is intended for final-stage refinement. Mathematical uniqueness of the convex solution was proved by Lamberg & Kaasalainen (2001). |
| **Input format** | Text file with lightcurves. Each data row contains: Julian Date (epoch), relative intensity, Sun-to-asteroid Cartesian ecliptic vector (x, y, z in AU), observer-to-asteroid Cartesian ecliptic vector (x, y, z in AU). A configuration/parameter file specifies initial spin axis (ecliptic longitude/latitude), period, spherical harmonics order, scattering law parameters, and optimization controls. |
| **Output format** | Convex polyhedral shape model (vertices and triangular facets), sidereal rotation period, spin axis direction (ecliptic longitude and latitude), scattering law parameters. Shape can be exported as vertex/facet lists. |
| **Known limitations** | Restricted to convex shapes only; all concavities are filled in. Levenberg-Marquardt converges to local minima, so a grid search over period and pole direction is required. Needs lightcurves from multiple apparitions (different viewing/illumination geometries) for a unique solution. Does not handle non-principal-axis rotation (tumbling) in the base version. Cannot determine absolute size from lightcurves alone. |
| **Source code** | Publicly available. GitHub mirror: [mkretlow/DAMIT-convex](https://github.com/mkretlow/DAMIT-convex). Official download: [DAMIT software page](https://damit.cuni.cz/pages/software_download). Includes `convexinv`, `conjgradinv`, and `lcgenerator` (synthetic lightcurve generator). |
| **Also includes** | `lcgenerator` tool for computing synthetic lightcurves from a given shape model (the forward/direct problem). |

### Extensions

- Modified versions handle: (i) disk-resolved data and non-convex features, (ii) YORP-induced spin-rate changes, and (iii) excited (tumbling) rotation states.
- The algorithm underpins the Asteroids@home distributed computing project and MPO LCInvert.

---

## 2. SAGE (Shaping Asteroids with Genetic Evolution) -- Bartczak & Dudzinski

| Property | Details |
|---|---|
| **Language** | Not publicly disclosed. The companion AsteroidFarm project by co-author Dudzinski is in C++, suggesting SAGE is likely C++ or Fortran. |
| **License** | Not publicly released as open-source software. Available by contacting the authors at Adam Mickiewicz University, Poznan, Poland. |
| **Algorithm** | Genetic (evolutionary) algorithm for non-convex shape inversion. Starts from two spheres with random spin axis orientation and approximate period. Random mutations are applied to shape vertices, period, spin axis, separation, and size ratio, generating a population each generation. Synthetic lightcurves are computed for each candidate and compared to observations via chi-squared. The best-fitting model seeds the next generation. Iteration continues until chi-squared converges. A weighting scheme prevents trapping in local minima. Shape is represented by 242 vertices on fixed radial rays; the mesh is refined using Catmull-Clark subdivision to 3,842 vertices / 7,680 faces. |
| **Input format** | Disk-integrated photometric lightcurves (relative magnitudes vs. time with geometric metadata). A radar module was added in 2016 to also accept radar delay-Doppler images. |
| **Output format** | Non-convex polyhedral shape model (vertices and triangular facets), rotation period, spin axis orientation. For binary asteroids, also outputs component separation and size ratio. |
| **Known limitations** | Computationally expensive due to the stochastic global search (many forward models per generation over many generations). Requires high-quality lightcurves at a range of solar phase angles for reliable non-convex features. Non-convex features become degenerate at low phase angles. Source code is not publicly available. |
| **Source code** | **Not publicly available.** Must be obtained by contacting the authors. |
| **Key references** | Bartczak & Dudzinski (2018), MNRAS 473, 5050. Applied to (90) Antiope, (809) Lundia, (433) Eros, (9) Metis. |

---

## 3. KOALA (Knitted Occultation, Adaptive-optics, and Lightcurve Analysis) -- Carry et al.

| Property | Details |
|---|---|
| **Language** | Not publicly disclosed. Based on Kaasalainen's inversion framework, likely Fortran/C. |
| **License** | Not publicly released. Academic use; available by contacting the authors (Benoit Carry, Observatoire de la Cote d'Azur). |
| **Algorithm** | Multi-data inversion combining three observation types simultaneously: (1) disk-resolved adaptive optics images (boundary contour extraction), (2) optical lightcurves, and (3) stellar occultation silhouettes. Uses the same convex/non-convex optimization as the Kaasalainen framework but adds contour-fitting for AO images and occultation chords. A thermophysical module was added in 2012 (with Delbo and Durech) to also use thermal infrared photometry. |
| **Input format** | Optical lightcurves (same format as convexinv), adaptive optics images (extracted boundary contours), stellar occultation chord timings, and optionally thermal infrared photometry. |
| **Output format** | 3D polyhedral shape model with absolute size, spin axis direction, rotation period, geometric albedo, and thermal inertia. |
| **Known limitations** | Relies on boundary contour extraction from AO images, which can be noisy or ambiguous. Requires multiple high-quality data types (AO, occultations, lightcurves) which are available for only a limited number of asteroids. Largely superseded by ADAM, which handles raw pixel data directly without contour extraction. |
| **Source code** | **Not publicly available.** No public GitHub or other repository found. |
| **Validation** | Validated by ESA Rosetta flyby of (21) Lutetia: spin axis accurate to 2 degrees, diameter to 2%, volume to 10%. |
| **Key references** | Carry et al. (2010, 2012), Planetary and Space Science 66, 200-212. |

---

## 4. ADAM (All-Data Asteroid Modeling) -- Viikinkoski, Kaasalainen & Durech

| Property | Details |
|---|---|
| **Language** | Primarily C (95.8%), with MATLAB (1.2%) and Python (1.1%) utilities for visualization and plotting. Requires libraries: KissFFT, Iniparser, Wcstools, LAPACK/BLAS. |
| **License** | CC Attribution 4.0 International |
| **Algorithm** | General non-convex shape reconstruction using subdivision surfaces. Handles all disk-resolved data types via uniform 2D Fourier transform comparison (comparing FT samples of model projections with data images). Supports Catmull-Clark and butterfly subdivision of control meshes. Optimization uses gradient-based methods. Regularization includes mesh smoothness, facet area homogeneity, and inertia constraints. Can fit albedo variegation (surface albedo maps). Supports Hapke and other scattering laws. |
| **Input format** | Configuration (.ini) files specifying data sources. Accepts any combination of: lightcurves (same format as convexinv), adaptive optics images, HST/FGS data, stellar occultation chord timings, range-Doppler radar images, disk-resolved thermal images. Initial shape can be an ellipsoid or a scaled convex model from lightcurve inversion. |
| **Output format** | Non-convex polyhedral shape model (vertices and triangular facets), spin state, scattering parameters, albedo maps. Shape files compatible with standard 3D formats (OBJ-style vertex/face lists). Includes MATLAB and Python scripts for displaying shapes, plotting projections and occultation fits. |
| **Known limitations** | Requires disk-resolved data (AO, radar, etc.) for non-convex reconstruction; lightcurves alone are insufficient for detailed non-convex features. Computational cost increases exponentially with subdivision depth. Setup requires significant expertise (no GUI; command-line with config files). Authors state "we do not offer any user support: the files are presented as is." Linux/gcc build environment assumed. |
| **Source code** | Publicly available on GitHub: [matvii/ADAM](https://github.com/matvii/ADAM). Also registered at ASCL: [ascl.net/1502.004](https://ascl.net/1502.004). |
| **Key references** | Viikinkoski, Kaasalainen & Durech (2015), A&A 576, A8. |

---

## 5. MPO LCInvert -- Brian D. Warner (Bdw Publishing)

| Property | Details |
|---|---|
| **Language** | Windows application (compiled binary). Internal implementation based on the Kaasalainen & Durech C code. |
| **License** | Commercial (proprietary). Sold by Bdw Publishing. Educational licenses available. Individual purchase; pricing on the Bdw Publishing website. |
| **Algorithm** | Implementation of the Kaasalainen et al. (2001) convex lightcurve inversion method. Includes period spectrum search, pole direction search (chi-squared grid), and convex shape optimization. Essentially a GUI wrapper around the convexinv/conjgradinv algorithms. |
| **Input format** | ALCDEF (Asteroid Lightcurve Data Exchange Format) files, which use a FITS-like keyword=value structure with data rows of JD, magnitude, magnitude error, and airmass. Also accepts lightcurves generated by MPO Canopus. |
| **Output format** | 3D convex shape model (visualized in the application), spin axis coordinates, rotation period. Can export shape model files. |
| **Known limitations** | Convex shapes only (same fundamental limitation as convexinv). Windows-only (runs on Mac via virtual machine). Commercial/closed-source, so algorithms cannot be inspected or modified. Requires understanding of lightcurve data quality requirements. Aimed primarily at amateur astronomers and citizen scientists, though used in professional research (Minor Planet Bulletin). |
| **Source code** | **Not available.** Closed-source commercial software. |
| **URL** | [minplanobs.org](https://minplanobs.org/BdwPub/php/mpolcinvert.php) / [bdwpublishing.com](http://bdwpublishing.com) |
| **Platform** | Windows (Intel); runs on Mac via Parallels or similar VM. Electronic download only. |
| **Key references** | Warner (2006), "A Practical Guide to Lightcurve Photometry and Analysis" (Springer). |

---

## 6. Durech's Lightcurve Inversion Codes and Related Repositories

### 6a. Asteroids@home / PeriodSearch

| Property | Details |
|---|---|
| **Language** | C++ (89.5%), C (7.7%), CUDA (1.7%), with shell scripts. Multiple build targets: CPU (SSE3, AVX/FMA), NVIDIA CUDA, AMD OpenCL. |
| **License** | GPL-3.0 (GNU General Public License v3) |
| **Algorithm** | Distributed implementation of the convex lightcurve inversion (same Kaasalainen/Durech algorithm as convexinv). The computationally intensive period-scanning step is parallelized: the period range is divided into independent work units distributed to volunteer computers via the BOINC framework. Each work unit evaluates the chi-squared fit over a sub-range of trial periods. Results are validated by redundant computation. |
| **Input format** | Lightcurve data in the convexinv format (JD, intensity, Sun/observer vectors). Work units are packaged by the BOINC server. |
| **Output format** | Convex shape models with spin axis direction and rotation period. Results published in peer-reviewed journals and submitted to DAMIT. |
| **Known limitations** | Convex shapes only. Requires the BOINC infrastructure for distributed execution. Individual work units test narrow period ranges, so the full solution requires aggregation. The CUDA/OpenCL versions require compatible GPU hardware. |
| **Source code** | GitHub: [AsteroidsAtHome/PeriodSearch](https://github.com/AsteroidsAtHome/PeriodSearch). CUDA 12 port: [JStateson/CUDA12_PeriodSearch](https://github.com/JStateson/CUDA12_PeriodSearch). |
| **URL** | Project website: [asteroidsathome.net](https://asteroidsathome.net/boinc/) |

### 6b. DAMIT-convex (GitHub mirror)

| Property | Details |
|---|---|
| **Language** | C (87.4%), Fortran (11.8%) |
| **License** | CC Attribution 4.0 International |
| **Description** | GitHub mirror of the official DAMIT download, containing `convexinv`, `conjgradinv`, `lcgenerator`, Fortran originals, PDF manuals, and example data. Maintained by M. Kretlow. |
| **Source code** | [mkretlow/DAMIT-convex](https://github.com/mkretlow/DAMIT-convex) |

### 6c. Asteroids-MDSM (Multi-Data Shape Modeling tools)

| Property | Details |
|---|---|
| **Language** | Python (requires numpy, pyglet) |
| **License** | Not specified |
| **Description** | Collection of tools and scripts for multi-data 3D shape modeling work, including a simple 3D shape file viewer for asteroids/comets. Utility scripts rather than a full inversion package. |
| **Source code** | [mkretlow/Asteroids-MDSM](https://github.com/mkretlow/Asteroids-MDSM) |

---

## 7. Python Implementations of Asteroid Shape Inversion

### 7a. sbpy (Small-Body Planetary Astronomy)

| Property | Details |
|---|---|
| **Language** | Python (Astropy affiliated package) |
| **License** | BSD 3-Clause (open source) |
| **Algorithm** | The `sbpy.shape.Kaasalainen` module is planned to implement the Kaasalainen convex inversion method natively in Python. The `sbpy.shape.Lightcurve` class provides period fitting (Lomb-Scargle), Fourier coefficient fitting, and spin pole orientation estimation. Can load shape models in VRML and OBJ formats and compute illumination/viewing geometry for surface facets. |
| **Input format** | Standard Python/Astropy data structures. Can read asteroid shapes in VRML, OBJ formats. Lightcurve data as time-series arrays. |
| **Output format** | Astropy-compatible Python objects; shape models, period estimates, pole solutions. |
| **Known limitations** | The lightcurve inversion module (`sbpy.shape.Kaasalainen`) is still **under development / not yet fully implemented** as of the latest documentation. Core shape loading and geometry calculations are functional, but full inversion capability may not be available yet. |
| **Source code** | GitHub: [NASA-Planetary-Science/sbpy](https://github.com/NASA-Planetary-Science/sbpy). Docs: [sbpy.readthedocs.io](https://sbpy.readthedocs.io/). Install via pip or conda-forge. |
| **Funding** | NASA PDART Grants 80NSSC18K0987 and 80NSSC22K0143. |

### 7b. Bayesian Lightcurve Inversion (Muinonen et al. 2020)

| Property | Details |
|---|---|
| **Language** | Not specified in the paper; methodology described for implementation with standard MCMC tools (e.g., Python + emcee). |
| **License** | Academic (methods published; no standalone code repository identified). |
| **Algorithm** | Bayesian inference framework for asteroid lightcurve inversion. Uses random-walk MCMC and importance sampling with a novel proposal density based on virtual observations generating virtual least-squares solutions. Simultaneously estimates rotation period, pole orientation, convex shape, and photometric phase-curve parameters with full posterior probability distributions. |
| **Input format** | Photometric lightcurves (demonstrated on Gaia DR2 photometry). |
| **Output format** | Posterior probability distributions for all parameters; convex shape model, spin state, period. |
| **Known limitations** | Computationally much more expensive than maximum-likelihood approaches (many MCMC samples needed). No standalone public code release identified. Convex shapes only in the published formulation. |
| **Key references** | Muinonen et al. (2020), A&A 642, A138. |

### 7c. Deep Learning Shape Inversion (Tang et al. 2025)

| Property | Details |
|---|---|
| **Language** | Python (PyTorch-based deep learning framework). |
| **License** | Not yet determined; paper published February 2025, code availability not confirmed. |
| **Algorithm** | Deep neural network mapping from photometric lightcurve data directly to 3D point cloud representations of asteroid shapes. Uses a modified PoinTR architecture with: (1) 1D DGCNN for feature extraction, (2) farthest point sampling and k-nearest neighbors for downsampling, (3) geometry-aware transformer for capturing global context and local details. A separate module predicts concave regions on the convex hull using deviation between non-convex asteroid lightcurves and their convex hull lightcurves. |
| **Input format** | Photometric lightcurves (normalized). Trained on synthetic lightcurves generated from known 3D models. |
| **Output format** | 3D point cloud of the asteroid shape. Concavity detection map (IoU = 0.89). |
| **Known limitations** | Requires large training datasets of synthetic lightcurves from known shapes. Generalization to real-world sparse/noisy data is still being validated. Tested on Lowell Observatory data for asteroids (3337) Milo and (1289) Kuta. New approach; not yet widely adopted or independently validated. Code availability unclear. |
| **Key references** | Tang et al. (2025), A&A 696, A55. arXiv: [2502.16455](https://arxiv.org/abs/2502.16455). |

### 7d. asteroid_lightcurve (jalalirs)

| Property | Details |
|---|---|
| **Language** | Python (90.3%), CSS (9.7%) |
| **License** | MIT |
| **Algorithm** | **Not an inversion tool.** Forward problem only: renders synthetic lightcurves from existing 3D shape models (.obj format). Provides orbital simulation modes. Retrieves observed lightcurves from the ALCP Database for comparison. |
| **Source code** | [jalalirs/asteroid_lightcurve](https://github.com/jalalirs/asteroid_lightcurve) |
| **Note** | Useful for validation/visualization but does not perform shape reconstruction. NASA Space Apps Challenge 2021 submission. |

---

## 8. DAMIT Database Tools and Utilities

### 8a. DAMIT Web Interface and Database

| Property | Details |
|---|---|
| **URL** | [damit.cuni.cz](https://damit.cuni.cz/) (formerly astro.troja.mff.cuni.cz/projects/damit/) |
| **Language** | PHP web interface, MySQL backend |
| **License** | Models and data freely available. Software tools under GPL/CC-BY 4.0. |
| **Contents** | 10,755 asteroids with 16,091 shape models and 5 tumblers (as of latest count). For each asteroid: polyhedral shape model (vertices/facets), sidereal rotation period, spin axis direction (ecliptic lambda/beta), scattering law parameters (c, a, d, k), and the photometric lightcurve data used for inversion. |
| **Data formats** | Shape models as text files (vertex coordinates + face indices); lightcurves in the convexinv 8-column format (JD, intensity, Sun xyz, observer xyz); PNG rendered images; downloadable parameter tables. |
| **Maintainer** | Josef Durech, Astronomical Institute, Charles University, Prague. |

### 8b. Software Available from DAMIT

The DAMIT software download page provides:

1. **convexinv** -- Full-parameter convex inversion using spherical harmonics (C).
2. **conjgradinv** -- Shape-only optimization using facet areas directly (C).
3. **lcgenerator** -- Synthetic lightcurve generator (forward problem) from a given shape and spin state (C).
4. **Fortran originals** -- The original Kaasalainen Fortran source code.
5. **PDF manuals** -- Documentation for each program.
6. **Example lightcurves** -- Test data for verifying installation.

All available under GPL and/or CC-BY 4.0.

### 8c. Related Utility: 3D Asteroid Catalogue

| Property | Details |
|---|---|
| **URL** | [3d-asteroids.space](https://3d-asteroids.space/) |
| **Description** | Interactive web-based catalogue with 3D models, orbital/physical parameters, and current orbital positions. Most lightcurve-inversion models sourced from DAMIT. Allows visual inspection of shape models in the browser. |

### 8d. ALCDEF (Asteroid Lightcurve Data Exchange Format)

| Property | Details |
|---|---|
| **URL** | [alcdef.org](https://alcdef.org/) |
| **Description** | Standardized format for exchanging raw asteroid time-series photometry. FITS-like keyword=value metadata with pipe/semicolon/tab-delimited data rows (JD, magnitude, mag error, airmass). Used as the primary input format for MPO LCInvert and a common source of lightcurves for inversion. |
| **Related tool** | [dnl-blkv/alcdef2json](https://github.com/dnl-blkv/alcdef2json) -- Python converter from ALCDEF to JSON. |

---

## 9. Additional Tools

### 9a. AsteroidFarm (Dudzinski)

| Property | Details |
|---|---|
| **Language** | C++ (93.6%), Lua (4.0%), GLSL (1.9%) |
| **License** | Not specified |
| **Algorithm** | Collection of algorithms for asteroid shape modeling. Currently supports synthetic observation generation (lightcurves, adaptive optics, radar delay-Doppler). SAGE method integration listed as TODO/in development. |
| **Source code** | [perkun/AsteroidFarm](https://github.com/perkun/AsteroidFarm) (Grzegorz Dudzinski, co-author of SAGE). 90 commits, actively developed. |
| **Note** | Not yet a complete inversion tool; primarily a forward-modeling and utility framework with SAGE integration planned. |

### 9b. SHAPE (Ostro et al.)

| Property | Details |
|---|---|
| **Language** | C |
| **License** | Academic; available from authors. |
| **Algorithm** | Radar + lightcurve inversion for non-convex shapes. Uses delay-Doppler radar images, continuous-wave spectra, and lightcurves. Vertex-based or spherical harmonics shape representation with Levenberg-Marquardt optimization. The original and most widely used tool for radar-based asteroid modeling. |
| **Known limitations** | Requires radar data (Arecibo or Goldstone), which is available for only a few hundred near-Earth asteroids. Complex setup. Not publicly distributed via GitHub. |
| **Note** | Included for completeness as it is a foundational tool in the field, predating ADAM and KOALA. |

---

## Comparative Summary Table

| Software | Language | License | Convex | Non-convex | Multi-data | Source Available | GUI |
|---|---|---|---|---|---|---|---|
| **convexinv** | C/Fortran | CC-BY 4.0 | Yes | No | No | Yes (GitHub + DAMIT) | No |
| **SAGE** | C++ (likely) | Not released | No | Yes | Lightcurves + radar | No | No |
| **KOALA** | C/Fortran (likely) | Not released | Yes | Partial | AO + occultations + LC + thermal | No | No |
| **ADAM** | C + MATLAB/Python | CC-BY 4.0 | Yes | Yes | AO + radar + occultations + LC + thermal | Yes (GitHub) | No |
| **MPO LCInvert** | Windows binary | Commercial | Yes | No | No | No (closed source) | Yes |
| **PeriodSearch** | C++/CUDA/OpenCL | GPL-3.0 | Yes | No | No | Yes (GitHub) | Minimal |
| **sbpy** | Python | BSD 3-Clause | Planned | No | No | Yes (GitHub) | No |
| **Deep Learning (Tang)** | Python/PyTorch | Unknown | Yes | Partial | No | Unclear | No |
| **AsteroidFarm** | C++ | Unspecified | Planned | Planned (SAGE) | Planned | Yes (GitHub) | No |
| **SHAPE** | C | Academic | No | Yes | Radar + LC | By request | No |

---

## Key Findings

1. **convexinv remains the workhorse.** The Kaasalainen/Durech convex inversion code is the most widely used tool, with over 16,000 models in DAMIT. It is freely available, well-documented, and has been validated against spacecraft imagery.

2. **ADAM is the most capable open-source multi-data tool.** It handles the widest variety of data types and produces non-convex models. Its CC-BY 4.0 license and GitHub availability make it the best option for researchers needing multi-data non-convex inversion.

3. **SAGE is powerful but inaccessible.** The genetic algorithm approach produces detailed non-convex shapes from lightcurves alone, but the code is not publicly available.

4. **KOALA has been largely superseded by ADAM.** ADAM inherits KOALA's contour-fitting capabilities while adding direct pixel-level fitting and more data types.

5. **Python ecosystem is still maturing.** sbpy's Kaasalainen module is under development. There is no production-ready, pure-Python lightcurve inversion package today. The deep learning approach (Tang et al. 2025) is promising but new and unvalidated at scale.

6. **MPO LCInvert fills a niche for amateurs.** It provides a GUI-based workflow for observers who want to attempt inversion of their own data without command-line expertise, but it is closed-source and convex-only.

7. **Distributed computing has scaled the problem.** Asteroids@home has applied convex inversion to hundreds of thousands of asteroids using BOINC, and the GPU-accelerated (CUDA/OpenCL) versions dramatically reduce computation time.

---

## Sources

- [mkretlow/DAMIT-convex (GitHub)](https://github.com/mkretlow/DAMIT-convex)
- [DAMIT Software Download](https://damit.cuni.cz/pages/software_download)
- [matvii/ADAM (GitHub)](https://github.com/matvii/ADAM)
- [ASCL - ADAM](https://ascl.net/1502.004)
- [AsteroidsAtHome/PeriodSearch (GitHub)](https://github.com/AsteroidsAtHome/PeriodSearch)
- [JStateson/CUDA12_PeriodSearch (GitHub)](https://github.com/JStateson/CUDA12_PeriodSearch)
- [perkun/AsteroidFarm (GitHub)](https://github.com/perkun/AsteroidFarm)
- [NASA-Planetary-Science/sbpy (GitHub)](https://github.com/NASA-Planetary-Science/sbpy)
- [sbpy Documentation](https://sbpy.readthedocs.io/en/latest/about.html)
- [MPO LCInvert](https://minplanobs.org/BdwPub/php/mpolcinvert.php)
- [Asteroids@home](https://asteroidsathome.net/boinc/)
- [DAMIT Database](https://damit.cuni.cz/)
- [3D Asteroid Catalogue](https://3d-asteroids.space/)
- [ALCDEF Database](https://alcdef.org/)
- [mkretlow/Asteroids-MDSM (GitHub)](https://github.com/mkretlow/Asteroids-MDSM)
- [jalalirs/asteroid_lightcurve (GitHub)](https://github.com/jalalirs/asteroid_lightcurve)
- [Bartczak & Dudzinski (2018) - SAGE](https://arxiv.org/abs/1904.08940)
- [Carry et al. (2012) - KOALA](https://www.sciencedirect.com/science/article/abs/pii/S0032063311003916)
- [Viikinkoski et al. (2015) - ADAM](https://www.aanda.org/articles/aa/full_html/2015/04/aa25259-14/aa25259-14.html)
- [Muinonen et al. (2020) - Bayesian Inversion](https://www.aanda.org/articles/aa/full_html/2020/10/aa38036-20/aa38036-20.html)
- [Tang et al. (2025) - Deep Learning Inversion](https://arxiv.org/abs/2502.16455)
- [Durech et al. (2010) - DAMIT](https://www.aanda.org/articles/aa/full_html/2010/05/aa12693-09/aa12693-09.html)
- [Asteroids@home paper](https://arxiv.org/abs/1511.08640)
