# Comprehensive Literature Review: Asteroid Light Curve Inversion Methods

## Date: 2026-02-09
## Purpose: Phase 1 research for asteroid lightcurve inversion pipeline

---

## Table of Contents

1. [Convex Inversion Method (Kaasalainen and Torppa 2001)](#1-convex-inversion-method)
2. [SAGE Genetic Evolution for Non-Convex Shapes (Bartczak and Dudzinski 2018)](#2-sage-genetic-evolution)
3. [Sparse Photometric Inversion (Durech et al. 2009-2010)](#3-sparse-photometric-inversion)
4. [ADAM All-Data Asteroid Modeling (Viikinkoski et al. 2015)](#4-adam-all-data-asteroid-modeling)
5. [Self-Shadowing Ray-Tracing in Photometric Models](#5-self-shadowing-ray-tracing)
6. [Scattering Laws: Lommel-Seeliger, Lambert, Hapke](#6-scattering-laws)
7. [Cellino et al. Sparse Inversion Work](#7-cellino-sparse-inversion)
8. [Muinonen et al. Scattering Law and Phase Function Work](#8-muinonen-scattering-law)
9. [Software Tool Survey](#9-software-tool-survey)
10. [Gap Analysis and Pipeline Implications](#10-gap-analysis)

---

## 1. Convex Inversion Method

### Key Papers

- **Kaasalainen and Torppa (2001a)**: "Optimization Methods for Asteroid Lightcurve Inversion. I. Shape Determination," *Icarus*, 153, 24-36.
- **Kaasalainen, Torppa, and Muinonen (2001b)**: "Optimization Methods for Asteroid Lightcurve Inversion. II. The Complete Inverse Problem," *Icarus*, 153, 37-51.
- **Kaasalainen, Torppa, and Piironen (2001c)**: "Models of Twenty Asteroids from Photometric Data," *Icarus*, 153, 52-65.

### Algorithmic Overview

The convex inversion method is the foundational technique for deriving asteroid 3D shape models and spin states from disk-integrated photometry (lightcurves). The method proceeds as follows:

**Step 1: Shape Parameterization via Spherical Harmonics**

The asteroid surface is parameterized through Gaussian surface density on the unit sphere. The surface curvature function (Gaussian image) is expanded in a series of real spherical harmonics:

    G(theta, phi) = sum_{l=0}^{l_max} sum_{m=-l}^{l} c_{lm} Y_{lm}(theta, phi)

where G is the Gaussian surface density, Y_{lm} are real spherical harmonics of degree l and order m, and c_{lm} are the coefficients to be optimized. Typical values used are l_max = 6-8 for the spherical harmonics order.

The Minkowski stability of the solution requires G >= 0 everywhere, which enforces convexity. The actual shape (vertex coordinates) is then recovered from the Gaussian image by a Minkowski reconstruction.

Alternatively, the shape can be represented directly as facet areas of a triangulated convex polyhedron with approximately 500 or more facets.

**Step 2: Scattering Law (Forward Model)**

For each triangular facet k with area A_k and outward normal n_k, the observed brightness is:

    L_model = sum_{k: visible and illuminated} A_k * S(mu_k, mu0_k, alpha)

where:
- mu_k = cos(theta_emission) = max(0, n_k . e_obs)
- mu0_k = cos(theta_incidence) = max(0, n_k . e_sun)
- alpha = phase angle (Sun-asteroid-observer angle)
- S is the scattering function (typically combined Lommel-Seeliger + Lambert)

The combined scattering law is:

    S(mu, mu0, alpha) = f(alpha) * [(1 - c_L) * mu0/(mu0 + mu) + c_L * mu0]

where c_L is the Lambert weight (typically 0.1 for dark asteroids) and f(alpha) is the phase function.

**Step 3: Levenberg-Marquardt Optimization**

The objective function is the total chi-squared:

    chi^2_total = chi^2_data + lambda_reg * chi^2_reg

where:
- chi^2_data = sum_i [(L_obs_i - L_model_i) / sigma_i]^2  over all observations
- chi^2_reg is a regularization term enforcing smoothness of the Gaussian surface density
- lambda_reg is the regularization weight

The parameter vector p = {c_lm, lambda_pole, beta_pole, P_rot, phi_0, c_L, ...} is optimized using the Levenberg-Marquardt algorithm:

    (J^T J + lambda_LM * diag(J^T J)) * delta_p = J^T * r

where J is the Jacobian matrix of partial derivatives of residuals with respect to parameters, r is the residual vector, and lambda_LM is the damping parameter.

**Step 4: Pole Search via Grid**

The spin axis direction (lambda, beta) is searched over a coarse grid (typically 5-degree steps in ecliptic longitude lambda and latitude beta). For each trial pole, the Levenberg-Marquardt optimization is run on the remaining parameters. The best-fit pole is the one minimizing chi^2_total.

### Limitations

1. **Convexity constraint**: By definition, the method cannot resolve concavities. The resulting model is the convex hull of the true shape.
2. **Phase angle dependence**: Nonconvex features are only detectable at high solar phase angles (alpha > 60 degrees), rare for main-belt asteroids (Durech and Kaasalainen 2003).
3. **Data requirements**: Reliable models require lightcurves from at least 3-4 apparitions with diverse viewing geometries.
4. **Scattering law degeneracy**: Shape and scattering law parameters are partially degenerate.
5. **No size information**: Convex inversion of relative photometry cannot determine absolute size.

### Open-Source Implementation

- **convexinv** (C): Available from https://github.com/mkretlow/DAMIT-convex
- **conjgradinv** (C): Companion program using conjugate gradient optimization with direct facet area parameterization
- **MPO LCInvert** (Windows, commercial): GUI wrapper around the Kaasalainen/Durech algorithms

---

## 2. SAGE Genetic Evolution for Non-Convex Shapes

### Key Paper

- **Bartczak and Dudzinski (2018)**: "Shaping asteroid models using genetic evolution (SAGE)," *MNRAS*, 473, 5050-5065.
- **Bartczak et al. (2014)**: "A new non-convex model of the binary asteroid 90 Antiope," *MNRAS*, 443, 1802-1809.

### Algorithmic Overview

SAGE (Shaping Asteroid models using Genetic Evolution) is a genetic algorithm approach that produces non-convex asteroid shapes directly from photometric lightcurve data.

**Step 1: Mesh Construction via Vertex Rays**

Starting from the icosahedron (12 vertices), Catmull-Clark subdivision is applied iteratively:
- After 2 subdivisions: 242 vertices
- After 4 subdivisions: 3842 vertices, 7680 triangular facets

Each vertex position is parameterized as a radial distance along a ray from the center. The genome of each individual encodes these radial distances.

**Step 2: Genetic Algorithm Operations**

- **Population initialization**: N shapes (N >= 50) generated by random perturbation of a seed shape
- **Fitness evaluation**: Synthetic lightcurves computed via ray-tracing (with self-shadowing) and compared to observations
- **Selection**: Tournament or rank-based selection
- **Crossover**: Parent genomes (vertex radial distances) combined to produce offspring
- **Mutation**: Random perturbation of vertex positions, constrained for mesh quality
- **Iteration**: 500+ generations until convergence

**Step 3: Fitness Function**

    fitness = 1 / chi^2 = 1 / sum_i [(L_obs_i - L_model_i)^2 / sigma_i^2]

where L_model_i is computed via full ray-tracing on the non-convex mesh.

### Key Assumptions

- Uniform geometric albedo
- Homogeneous mass distribution
- Principal axis rotation (no tumbling)

### Validation

Validated on 433 Eros (NEAR Shoemaker ground truth), 9 Metis, and binary 90 Antiope. Eros non-convex model showed excellent agreement with spacecraft shape.

### Limitations

1. **Computational cost**: Ray-tracing per fitness evaluation, orders of magnitude more expensive than convex inversion
2. **Non-uniqueness**: Multiple non-convex shapes can produce identical lightcurves
3. **Local optima**: GA may converge to suboptimal solutions
4. **Data requirements**: Non-convex features need high phase-angle data
5. **Albedo uniformity**: Real albedo variegation can mimic shape features

### Software Availability

**SAGE is NOT publicly available.** Developed at Adam Mickiewicz University, Poznan, Poland. Contact authors directly. Funded by EU Horizon 2020 grant 687378.

---

## 3. Sparse Photometric Inversion

### Key Papers

- **Kaasalainen (2004)**: *A&A*, 422, L39-L42. First proof that sparse photometry can solve the inversion problem.
- **Durech et al. (2009)**: "Asteroid models from combined sparse and dense photometric data," *A&A*, 493, 291-297.
- **Durech et al. (2010)**: "DAMIT: a database of asteroid models," *A&A*, 513, A46.
- **Durech et al. (2016)**: "Asteroid models from the Lowell Photometric Database," *A&A*, 587, A48.
- **Durech et al. (2018)**: "Asteroid models from Lowell + WISE data," *A&A*, 617, A57.
- **Hanus et al. (2011)**: "Pole-latitude distribution study," *A&A*, 530, A134.
- **Hanus et al. (2013)**: "Combined dense and sparse photometry," *A&A*, 551, A67.

### Method Description

Sparse photometric inversion adapts the convex inversion framework to handle datasets with only a few measurements per night, as produced by all-sky surveys.

**Key differences from dense lightcurve inversion:**

1. **Absolute calibration is essential**: Sparse data points from different nights must be on the same photometric scale.
2. **Coarse grid search**: Period searched with step size ~10^{-5} to 10^{-6} hours; pole searched on 5-degree grid
3. **Data requirements**: Typically >= 100 calibrated measurements across >= 3 apparitions
4. **Combined dense + sparse**: Most powerful approach uses dense lightcurves for shape detail plus sparse points for pole and brightness constraints

### Survey Data Characteristics

| Survey | Typical accuracy | Points per asteroid | Time baseline |
|--------|-----------------|--------------------|----|
| Lowell/MPC archival | 0.1-0.2 mag | hundreds | 1998-2011 |
| Pan-STARRS | ~0.02 mag | tens to hundreds | 2010-present |
| Gaia | ~0.01 mag | ~70 per asteroid | 2014-present |
| ZTF | ~0.02 mag | tens to hundreds | 2018-present |
| ATLAS | ~0.02 mag | hundreds | 2015-present |
| LSST/Rubin | ~0.01 mag | thousands (projected) | 2025-2035 |

### Results Scale

- DAMIT database: >16,000 models as of 2024
- Durech et al. (2016): 328 new models from Lowell alone
- Asteroids@home BOINC project: distributed computing for mass inversion

### Limitations

1. Low shape resolution from sparse data
2. Severe period aliasing
3. Photometric systematics across surveys
4. Convex-only models

---

## 4. ADAM All-Data Asteroid Modeling

### Key Papers

- **Viikinkoski, Kaasalainen, and Durech (2015)**: *A&A*, 576, A8.
- **Viikinkoski et al. (2017)**: *A&A*, 607, A117.

### Method Description

ADAM simultaneously inverts: disk-integrated photometry, AO images, stellar occultations, radar delay-Doppler, thermal IR, and interferometric data.

**Key innovation**: Uniform handling of all disk-resolved data via 2D Fourier transform projections.

### Shape Representations

1. **Octantoids**: Spherical harmonics-based, smooth curved surfaces, parameterized as r(theta,phi) = sum a_{lm} Y_{lm}. Can represent non-star-shaped surfaces.
2. **Subdivision surfaces**: Control-point meshes with Loop/Catmull-Clark subdivision. Sharper local features.

**Reliability check**: If both representations converge to similar shapes, data are sufficient.

### Regularization Functions

- eta: Smoothness (adjacent facet normal differences)
- gamma_1: Angular constraint (sharp angle prevention)
- gamma_2: Convex regularization (optional concavity suppression)
- gamma_3: Area homogeneity (degenerate facet prevention)

### Total Objective

    chi^2_total = sum_d w_d * chi^2_d + sum_r lambda_r * R_r

Optimization via Levenberg-Marquardt.

### Limitations

1. Requires disk-resolved data for reliable non-convex results
2. Initialization-sensitive (starts from convex seed)
3. Few hundred asteroids have sufficient multi-data coverage
4. Computationally expensive

### Open-Source

GitHub: https://github.com/matvii/ADAM (Fortran + C)

---

## 5. Self-Shadowing Ray-Tracing in Asteroid Photometric Models

### Key Reference

- **Durech and Kaasalainen (2003)**: "Photometric signatures of highly nonconvex and binary asteroids," *A&A*, 404, 709-714.

### Ray-Tracing Algorithm

For each facet k on a non-convex mesh:

1. **Illumination check**: Cast ray from centroid toward Sun. Intersection with another facet means shadow.
2. **Visibility check**: Cast ray from centroid toward observer. Intersection means occlusion.
3. **Brightness**: Only facets that are both illuminated AND visible contribute:

       L_total = sum_{k: visible AND illuminated} A_k * S(mu_k, mu0_k, alpha)

### BVH Acceleration

- Binary tree of axis-aligned bounding boxes (AABBs) enclosing facet groups
- Ray-box intersection test before descending to children
- Average complexity: O(N_facets * log(N_facets)) per evaluation
- Construction via Surface Area Heuristic (SAH)

### Phase Angle Dependence

- alpha < 30 deg: Self-shadowing minimal
- alpha > 60 deg: Shadowing from concavities becomes pronounced
- Main-belt asteroids: max alpha ~20-30 deg (small effect)
- NEOs: alpha can exceed 90 deg (essential for accurate modeling)

### Key Finding (Durech and Kaasalainen 2003)

- Nonconvexity measure defined for known asteroid shapes
- Only binary/bifurcated shapes resolvable for main-belt asteroids
- Topologically star-like asteroids modeled as convex bodies
- Observations at alpha > 60 deg critical for non-convex recovery

---

## 6. Scattering Laws: Lommel-Seeliger, Lambert, Hapke

### 6.1 Lambert (Lambertian) Scattering

    r_Lambert = (A / pi) * mu_0

where A is albedo, mu_0 = cos(incidence angle). Radiance independent of emission angle. Unrealistic for most asteroids.

### 6.2 Lommel-Seeliger Scattering

Derived from radiative transfer (Chandrasekhar 1960):

    r_LS = (w_0 / (4*pi)) * mu_0 / (mu_0 + mu) * p(alpha)

Properties: valid for low-albedo surfaces, single scattering only, computationally simple, default in convex inversion codes.

### 6.3 Combined Lommel-Seeliger + Lambert

Standard for lightcurve inversion:

    S(mu, mu_0, alpha) = f(alpha) * [(1 - c_L) * mu_0/(mu_0 + mu) + c_L * mu_0]

Typical values: c_L ~ 0.1 (dark asteroids), c_L ~ 0.5 (bright asteroids).

### 6.4 Hapke Scattering Model (Hapke 1993, 2012)

Full bidirectional reflectance:

    r(i, e, alpha) = (w/(4*pi)) * (mu_0/(mu_0+mu)) * [p(alpha)*B_SH(alpha) + M(mu_0,mu)] * B_CB(alpha) * S(theta_bar,i,e,alpha)

Components:
- w: single scattering albedo
- p(alpha): Henyey-Greenstein phase function, p_HG = (1-g^2)/(1+2g*cos(alpha)+g^2)^{3/2}
- B_SH: Shadow-Hiding Opposition Effect, B_SH = 1 + B_S0/(1+tan(alpha/2)/h_S)
- M: Multiple scattering via H-functions, M = H(mu_0,w)*H(mu,w) - 1
- B_CB: Coherent Backscatter Opposition Effect
- S(theta_bar): Macroscopic roughness correction

Parameters (5-7 free): w, g, B_S0, h_S, B_C0, h_C, theta_bar

**Limitation**: Many correlated parameters, values may not correspond to physical properties, computationally expensive, degenerate with shape in inversion.

### Practical Recommendation

Combined LS + Lambert with linear phase coefficient is sufficient for lightcurve inversion. Hapke parameters are better constrained from disk-resolved data.

---

## 7. Cellino et al. Sparse Inversion Work

### Key Papers

- **Cellino et al. (2009)**: *A&A*, 506, 935-954. Genetic inversion of Hipparcos data.
- **Cellino et al. (2015)**: *PSS*, 118, 221-226. LS scattering on ellipsoids.
- **Santana-Ros et al. (2015)**: *MNRAS*, 450, 333-341. Testing Gaia inversion.

### Method

Genetic algorithm with triaxial ellipsoid shape model:
- Parameters: P, (lambda, beta), axial ratios (a/c, b/c), phi_0, phase slope
- Later extended to "Cellinoid" (8 ellipsoidal octants)
- Lommel-Seeliger scattering law on ellipsoidal surfaces

### Application to Gaia

Primary method for Gaia DPAC asteroid photometry inversion. Applied to 22,815 objects from Gaia DR3.

### Limitations

1. Low shape resolution (ellipsoid/Cellinoid only)
2. Failures for quasi-spherical or pole-on objects
3. Speed vs. accuracy trade-off

---

## 8. Muinonen et al. Scattering Law and Phase Function Work

### Key Papers

- **Muinonen et al. (2010)**: *Icarus*, 209, 542-555. H, G1, G2 phase function.
- **Muinonen et al. (2009)**: *M&PS*, 44, 1937-1946. Linear-exponential modeling.
- **Muinonen and Lumme (2015)**: *A&A*, 584, A23. LS ellipsoidal brightness.
- **Muinonen et al. (2020)**: *A&A*, 642, A138. Bayesian lightcurve inversion.

### H, G1, G2 Phase Function (IAU standard since 2012)

    V(alpha) = H - 2.5*log10[G1*Phi_1(alpha) + G2*Phi_2(alpha) + (1-G1-G2)*Phi_3(alpha)]

- H: absolute magnitude
- G1, G2: phase parameters (G1+G2 <= 1, both >= 0)
- Phi_1, Phi_2, Phi_3: cubic spline basis functions
- Simplified H, G12 two-parameter version for sparse data

### Linear-Exponential Phase Curve

    V(alpha) = V(0) + b*alpha - a_1*exp(-alpha/d_1)

### Bayesian Lightcurve Inversion (2020)

Full MCMC treatment: posterior distributions for all parameters, formal uncertainty quantification, Lommel-Seeliger scattering on triaxial ellipsoids.

### Taxonomic Implications

G12 parameters correlate with taxonomy: C-complex (low), S-complex (intermediate), E-type (high).

---

## 9. Software Tool Survey

### 9.1 convexinv / conjgradinv (DAMIT tools)

| Property | Value |
|----------|-------|
| Developer | Kaasalainen (Fortran original), Durech (C port) |
| Language | C |
| License | Freely available (academic use) |
| Source | https://github.com/mkretlow/DAMIT-convex |
| Sparse data | Yes |
| Non-convex | No |
| Self-shadowing | No |
| Optimization | Levenberg-Marquardt |

### 9.2 MPO LCInvert

| Property | Value |
|----------|-------|
| Developer | Brian D. Warner |
| Language | Windows application (proprietary) |
| License | Commercial |
| Sparse data | Yes (Catalina + dense) |
| Non-convex | No |
| Self-shadowing | No |
| Key feature | Guided Inversion Wizard GUI |

### 9.3 SAGE

| Property | Value |
|----------|-------|
| Developer | Bartczak, Dudzinski (AMU Poznan) |
| Language | Not publicly documented |
| License | NOT publicly available |
| Sparse data | Not documented |
| Non-convex | YES |
| Self-shadowing | YES |
| Mesh | 3842 vertices, 7680 facets |
| Optimization | Genetic algorithm |

### 9.4 KOALA

| Property | Value |
|----------|-------|
| Developer | Benoit Carry et al. |
| Language | Not publicly documented |
| License | NOT publicly available |
| Non-convex | Yes (multi-data) |
| Self-shadowing | Yes |
| Data types | Lightcurves + AO + occultations |
| Validation | Rosetta/Lutetia: spin to 2 deg, diameter to 2% |

### 9.5 ADAM

| Property | Value |
|----------|-------|
| Developer | Matti Viikinkoski |
| Language | Fortran + C |
| License | Open source |
| Source | https://github.com/matvii/ADAM |
| Sparse data | Yes |
| Non-convex | Yes |
| Self-shadowing | Yes |
| Data types | All (lightcurves, AO, occultations, radar, thermal, ALMA) |

### 9.6 LSST/Rubin Observatory

| Property | Value |
|----------|-------|
| Developer | Rubin Observatory team |
| Language | Python |
| License | Open source (GPL) |
| Shape inversion | NOT in core pipeline (as of 2026) |
| Products | Orbits, lightcurves, colors, rotation periods |
| Projected yield | ~5M MBAs, ~127K NEOs |
| Gap | Shape inversion is community value-added product |

### Comparison Table

| Feature | convexinv | LCInvert | SAGE | KOALA | ADAM | Rubin |
|---------|-----------|----------|------|-------|------|-------|
| Open source | Yes | No | No | No | Yes | Yes* |
| Non-convex | No | No | Yes | Yes | Yes | N/A |
| Self-shadow | No | No | Yes | Yes | Yes | N/A |
| Sparse data | Yes | Yes | ? | No | Yes | Yes* |
| Multi-data | No | No | Partial | Yes | Yes | No |
| Language | C | Windows | ? | ? | Fortran | Python |

---

## 10. Gap Analysis and Pipeline Implications

### Identified Gaps

1. **No open-source non-convex photometry-only inverter**: SAGE is not public. ADAM and KOALA need disk-resolved data.
2. **No unified dense+sparse non-convex pipeline**: Existing tools are either convex+sparse or non-convex+dense.
3. **No integrated BVH-accelerated ray-tracer for asteroid photometry**: Performance optimization not well-documented.
4. **LSST/Rubin shape model gap**: Unprecedented data but no inversion pipeline.
5. **Uncertainty quantification**: Only Muinonen et al. (2020) provide Bayesian UQ.

### Our Pipeline Should Implement

1. Convex inversion (Kaasalainen and Torppa 2001) as initial seed
2. GA non-convex solver (SAGE-inspired, Bartczak and Dudzinski 2018) using convex seed
3. Self-shadowing ray-tracer with BVH acceleration
4. Combined LS + Lambert scattering as default
5. Sparse + dense data fusion (Durech et al. 2009)
6. Bootstrap/jackknife uncertainty quantification
7. H, G1, G2 phase function (Muinonen et al. 2010) for absolute photometry

---

## BibTeX References

All 30 BibTeX entries are provided in `/sources.bib`. See that file for complete citation data.

---

*Compiled 2026-02-09, Phase 1: Problem Analysis and Literature Review.*
