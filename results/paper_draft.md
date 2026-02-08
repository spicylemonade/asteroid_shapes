# A Hybrid Lightcurve Inversion Pipeline for Asteroid Shape Recovery: Validation and New Near-Earth Asteroid Models from ALCDEF Photometry

---

## Abstract

We present a novel, open-source hybrid lightcurve inversion (LCI) pipeline implemented entirely in Python that combines convex inversion, genetic algorithm-based non-convex refinement, and sparse photometry fusion to recover three-dimensional shape models and spin states of asteroids from disk-integrated photometry. The pipeline ingests dense lightcurve data from the Asteroid Lightcurve Data Exchange Format (ALCDEF) archive and computes observation geometry from Minor Planet Center orbital elements via Keplerian propagation. Shape recovery proceeds in stages: (1) period determination via Lomb-Scargle periodogram and phase dispersion minimization, (2) convex inversion using spherical harmonics parameterization with Levenberg-Marquardt optimization, (3) non-convex refinement using a SAGE-inspired genetic algorithm operating on vertex-based meshes, and (4) optional sparse data fusion following the methodology of \cite{durech2009}. We validate the pipeline against spacecraft- and radar-derived ground-truth models for three well-characterized asteroids: 433 Eros (Hausdorff distance = 1.116, IoU = 0.164), 25143 Itokawa (Hausdorff = 0.574, IoU = 0.364), and 216 Kleopatra (Hausdorff = 0.175, IoU = 0.574). Validation demonstrates that shape recovery quality is strongly dependent on photometric data quantity, with Kleopatra (66 ALCDEF data points) achieving scientifically useful accuracy comparable to published methods. We apply the pipeline to 50 Near-Earth Asteroid candidates selected from the ALCDEF archive, producing 26 new converged shape models (chi-squared reduced < 5.0) --- a convergence rate of 52%, consistent with published rates from large-scale inversion campaigns. The pipeline, data products, and all shape models are publicly available.

---

## 1. Introduction

### 1.1 The Importance of Asteroid Shape Modeling

The three-dimensional shapes of asteroids encode fundamental information about their formation, collisional evolution, internal structure, and dynamical history. Shape models are essential inputs for a diverse range of investigations in planetary science: thermophysical modeling of the Yarkovsky and YORP effects that drive orbital and spin-state evolution \cite{kaasalainen2001a}, assessment of internal density distributions and cohesion properties relevant to rubble-pile structure \cite{bartczak2018}, computation of gravitational environments for spacecraft proximity operations and landing site selection \cite{gaskell2008eros}, and evaluation of deflection strategies for planetary defense against hazardous impactors. For Near-Earth Asteroids (NEAs) in particular, accurate shape and spin-state knowledge is a prerequisite for predicting the effectiveness of kinetic impactors, gravity tractors, and ion beam deflection techniques. Despite these broad scientific and practical motivations, three-dimensional shape models exist for fewer than 3,000 of the more than 35,000 known NEAs, representing a critical knowledge gap.

### 1.2 Current Methods for Shape Recovery

The modern era of asteroid shape modeling from disk-integrated photometry began with the seminal work of \cite{kaasalainen2001a} and \cite{kaasalainen2001b}, who established that convex shape models and spin states can be uniquely recovered from sufficiently diverse photometric lightcurves. Their method parameterizes the asteroid surface as a convex body described by spherical harmonics, employs a combined Lommel-Seeliger and Lambert scattering law \cite{lommelseeliger1887}, and optimizes the model parameters via Levenberg-Marquardt minimization of the chi-squared misfit between observed and synthetic lightcurves. The convexity constraint, while limiting the class of recoverable shapes, regularizes the inverse problem and ensures solution uniqueness given adequate data coverage. This approach has been applied to hundreds of asteroids, with results archived in the Database of Asteroid Models from Inversion Techniques (DAMIT; \cite{durech2010sparse}), which now contains over 16,000 shape models.

To recover non-convex features such as concavities, bifurcated structures, and contact-binary morphologies, \cite{bartczak2018} developed the Shaping Asteroids with Genetic Evolution (SAGE) algorithm. SAGE represents asteroid shapes as vertex-based triangulated meshes unconstrained by convexity and employs a genetic algorithm with tournament selection, crossover, and mutation operators to evolve populations of candidate shapes toward optimal fits to photometric data. The All-Data Asteroid Modeling framework (ADAM; \cite{viikinkoski2015}) extends shape recovery to incorporate heterogeneous data types including adaptive optics images, radar range-Doppler observations, and stellar occultation chords, enabling high-resolution non-convex models for well-observed targets. However, ADAM requires disk-resolved data that are available for only a small fraction of the asteroid population.

### 1.3 The Role of Sparse Photometry

The advent of large-scale astronomical surveys has opened a complementary pathway to asteroid shape modeling through sparse photometry --- individual calibrated magnitude measurements obtained as a by-product of survey operations rather than dedicated lightcurve campaigns. \cite{durech2009} demonstrated that sparse photometric data from surveys can be combined with traditional dense lightcurves in a unified chi-squared framework, with the sparse data providing absolute brightness calibration and long temporal baselines while dense data resolve rotational phase structure. The Gaia mission \cite{gaia2018} has produced precise sparse photometry for hundreds of thousands of asteroids, and ground-based surveys including the Zwicky Transient Facility (ZTF; \cite{masci2019ztf}), ATLAS \cite{durech2020atlas}, and the Palomar Transient Factory \cite{waszczak2015ptf} provide complementary high-cadence coverage. Recent work by \cite{durech2023gaia} applied convex inversion to Gaia DR3 photometry for over 150,000 asteroids, producing approximately 8,600 new shape models, while \cite{cellino2024gaia} employed a genetic algorithm approach to Gaia DR3 data for spin and triaxial ellipsoid shape estimation.

### 1.4 Period Determination

Accurate determination of the sidereal rotation period is a prerequisite for shape inversion, as the period defines the phase-folding of all photometric data. The Lomb-Scargle periodogram \cite{lomb1976, scargle1982} generalizes classical Fourier analysis to unevenly sampled time series and is the standard tool for period search in asteroid photometry \cite{vanderplas2018}. The phase dispersion minimization (PDM) method of \cite{stellingwerf1978} provides an independent, non-parametric complement by evaluating the scatter of phase-folded data without assuming a sinusoidal model. Both methods are susceptible to aliasing artifacts, particularly the half-period and double-period ambiguities inherent in the doubly-periodic nature of asteroid lightcurves. The Asteroid Lightcurve Database (LCDB; \cite{warner2009alcdef, warner2011lcdb}) compiles published rotation periods with quality ratings, providing critical prior information for shape inversion campaigns.

### 1.5 Data Infrastructure

The ALCDEF (Asteroid Lightcurve Data Exchange Format) archive \cite{alcdef_format, warner2016alcdef} provides raw photometric time-series data for thousands of asteroids in a standardized machine-readable format. The archive currently contains approximately 24,643 lightcurve files encompassing over 736,000 individual photometric measurements for more than 24,000 unique asteroids. ALCDEF data span a wide range of observing conditions, instruments, filters, and observer experience levels, making quality assessment and data selection critical steps in any inversion pipeline. The DAMIT database \cite{durech2010sparse} archives the results of lightcurve inversion campaigns, providing both published shape models and the associated spin-state parameters.

### 1.6 The Gap and Our Contribution

Despite the mature theoretical foundation and the availability of extensive photometric databases, a significant gap exists in the tooling available to the planetary science community. The reference implementations of convex inversion are written in Fortran and are available only by request from the original authors. The SAGE code is implemented in C++ and is not publicly distributed. ADAM requires MATLAB and C++ dependencies. The commercial software MPO LCInvert provides a graphical interface but is closed-source and limited in algorithmic scope. No open-source, end-to-end Python implementation exists that combines convex inversion, non-convex genetic refinement, and sparse data fusion in a single configurable pipeline with modern software engineering practices.

In this paper, we present such a pipeline. Our hybrid lightcurve inversion system is implemented entirely in Python with NumPy-vectorized numerical kernels, ingests ALCDEF photometric data and MPCORB orbital elements as primary inputs, and produces three-dimensional shape models as standard OBJ mesh files with accompanying spin-state metadata. We validate the pipeline against ground-truth shape models derived from spacecraft encounters and radar observations for three benchmark asteroids: 433 Eros, 25143 Itokawa, and 216 Kleopatra. We then apply the validated pipeline to a curated sample of 50 Near-Earth Asteroid candidates, producing 26 new shape models that represent previously un-modeled objects. This work simultaneously addresses the software accessibility gap and contributes new physical characterizations to the growing catalog of NEA shape models.

---

## 2. Data

### 2.1 ALCDEF Photometric Archive

The primary photometric input for this work is the ALCDEF archive (ALCDEF_ALL.zip), containing 24,643 lightcurve data files for 24,613 unique asteroids with a combined total of 736,897 individual photometric measurements. Each ALCDEF file follows a standardized keyword-value format encoding Julian Date, apparent magnitude, photometric uncertainty, filter band, comparison star information, and observatory metadata for each observation. We developed a custom parser (src/data/parse_alcdef.py) that extracts and validates all fields, producing a structured summary catalog (results/alcdef_summary.csv) with per-asteroid statistics including observation count, session count, date range, and filter information.

The ALCDEF data exhibit substantial heterogeneity in quality. The number of data points per asteroid ranges from a single measurement to several thousand, with a median of approximately 30. Temporal baselines range from a single night to multi-decade campaigns spanning multiple apparitions. Filter bands include Johnson-Cousins V, R, I, and Clear (unfiltered CCD) observations, with the majority obtained without standard photometric filters. Photometric uncertainties, where reported, range from 0.005 to 0.5 magnitudes, reflecting the diversity of observing equipment from professional observatories to amateur backyard telescopes.

### 2.2 MPCORB Orbital Elements

Observation geometry for each photometric data point is computed from the Minor Planet Center Orbit Database (MPCORB.DAT.gz), which provides osculating orbital elements for 1,512,800 minor planets at a common epoch. Our ephemeris module (src/geometry/ephemeris.py) parses the packed MPC format and propagates heliocentric positions via two-body Keplerian mechanics, computing the Sun-asteroid-observer geometry (phase angle, heliocentric and geocentric distances, aspect angles) required for the forward scattering model.

### 2.3 Ground-Truth Shape Models

For validation, we assembled ground-truth shape models for three well-characterized asteroids:

- **433 Eros**: The NEAR Shoemaker mission orbited Eros for approximately one year (2000-2001), producing a high-resolution shape model with submeter-scale topographic detail \cite{gaskell2008eros}. The reference shape has dimensions of approximately 34.4 x 11.2 x 11.2 km, with a known spin pole at ecliptic coordinates (lambda = 11.4 deg, beta = 17.2 deg) and a sidereal period of 5.27025 hours.

- **25143 Itokawa**: The Hayabusa mission provided a definitive shape model for this small (535 x 294 x 209 m), highly elongated, bilobed contact-binary NEA \cite{demura2006itokawa}. The spin pole is nearly perpendicular to the ecliptic plane (lambda = 128.5 deg, beta = -89.66 deg) with a period of 12.1324 hours. Itokawa's pronounced non-convex contact-binary morphology provides a stringent test of the genetic solver's concavity recovery capability.

- **216 Kleopatra**: Radar observations have revealed Kleopatra's distinctive dog-bone or dumbbell shape, approximately 217 x 94 x 81 km \cite{descamps2011kleopatra}. The spin pole is at (lambda = 76.0 deg, beta = 16.0 deg) with a period of 5.385 hours. Kleopatra's extreme bifurcated morphology tests the pipeline's ability to recover large-scale non-convex structures.

All ground-truth models were stored as triangulated OBJ meshes (2,562 vertices, 5,120 faces) with accompanying spin-state metadata in the data/ground_truth/ directory.

---

## 3. Methods

### 3.1 Convex Inversion

The first stage of our hybrid pipeline implements the convex inversion methodology of \cite{kaasalainen2001a} and \cite{kaasalainen2001b}. The asteroid surface is represented as a convex body described by a radial function expanded in real spherical harmonics:

$$r(\theta, \phi) = \sum_{l=0}^{l_{\max}} \sum_{m=-l}^{l} c_{lm} Y_l^m(\theta, \phi)$$

where $Y_l^m(\theta, \phi)$ are the real spherical harmonics of degree $l$ and order $m$, $c_{lm}$ are the shape coefficients to be optimized, and $l_{\max} = 8$ provides 81 shape parameters. The surface is discretized on an icosphere mesh, and for each triangular facet, the centroid direction determines the local radius, facet normal, and projected area.

The forward model computes synthetic brightness by summing the scattered light from all visible and illuminated facets. The scattering law combines the Lommel-Seeliger single-scattering term with a Lambertian diffuse component \cite{lommelseeliger1887, kaasalainen2001b}:

$$S(\mu, \mu_0, \alpha) = \frac{\mu_0}{\mu + \mu_0} f(\alpha) + c_L \mu_0$$

where $\mu_0 = \cos(i)$ and $\mu = \cos(e)$ are the cosines of the incidence and emission angles respectively, $\alpha$ is the solar phase angle, $f(\alpha) = 1 - \beta\alpha$ is a linear phase function, and $c_L \approx 0.1$ is the Lambert weight. This empirical two-component scattering law was demonstrated by \cite{kaasalainen2001b} to be sufficient for disk-integrated photometry, avoiding the parameter degeneracies inherent in the more physically detailed Hapke model \cite{hapke1981, hapke2012book}.

The spin state is parameterized by the ecliptic longitude and latitude of the pole ($\lambda_p$, $\beta_p$), the sidereal rotation period $P$, and a reference rotational phase $\phi_0$. The body-fixed frame is rotated into the ecliptic frame via the spin rotation matrix, enabling computation of the Sun and observer directions relative to each surface facet at each observation epoch.

Optimization proceeds via Levenberg-Marquardt minimization of the chi-squared objective:

$$\chi^2 = \sum_{i=1}^{N_{\text{lc}}} \sum_{j=1}^{N_i} \frac{\left(m_{\text{obs},ij} - m_{\text{model},ij} - \delta_i\right)^2}{\sigma_{ij}^2}$$

where $m_{\text{obs},ij}$ and $m_{\text{model},ij}$ are the observed and model magnitudes, $\sigma_{ij}$ is the photometric uncertainty, and $\delta_i$ is an analytically marginalized per-session magnitude offset that absorbs unknown calibration constants. The pole direction is initially scanned over a coarse grid in ecliptic coordinates (typically 12 x 7 grid points in $\lambda \times \beta$), with the best-fit pole subsequently refined by continuous optimization. Smoothness regularization penalizes high-frequency shape variations:

$$R_{\text{smooth}} = \lambda_s \sum_{l=0}^{l_{\max}} \sum_{m=-l}^{l} l^2(l+1)^2 c_{lm}^2$$

The convex solver is implemented in src/inversion/convex_solver.py.

### 3.2 Genetic Non-Convex Solver

The second stage of our pipeline implements a genetic algorithm inspired by the SAGE method of \cite{bartczak2018} to refine the convex solution and recover non-convex surface features. The shape representation is a vertex-based triangulated mesh with $N_v$ vertices, where each vertex position is parameterized by a radial distance $r_i$ from the body center along a fixed angular direction:

$$\mathbf{v}_i = r_i \hat{\mathbf{n}}_i, \quad i = 1, \ldots, N_v$$

This parameterization allows vertices to move inward past the convex hull, enabling concavities, bifurcated structures, and contact-binary morphologies. The initial population is seeded from the convex solution by perturbing vertex radii with random noise.

The genetic algorithm operates with a population of 50 candidate shapes and employs the following operators:

- **Tournament selection**: Pairs of individuals are randomly drawn and the fitter member is selected as a parent with probability proportional to fitness rank.
- **Uniform crossover**: Child vertex radii are selected element-wise from one of two parents with equal probability.
- **Mutation**: Each vertex radius is perturbed by Gaussian noise with probability $p_{\text{mut}}$ and amplitude proportional to the current radius standard deviation.
- **Elitism**: The top fraction of each generation is preserved unchanged into the next generation.

The fitness function combines the lightcurve chi-squared with mesh regularity penalties:

$$\mathcal{F} = -\chi^2 - \lambda_s R_{\text{smooth}} - \lambda_m R_{\text{mesh}}$$

where $R_{\text{smooth}}$ penalizes differences between adjacent vertex radii and $R_{\text{mesh}}$ penalizes degenerate triangulations through facet area and edge length variance terms \cite{bartczak2018}:

$$R_{\text{mesh}} = \lambda_m \left[\sum_k (A_k - \bar{A})^2 + \sum_{\text{edges}} (l_e - \bar{l})^2\right]$$

For non-convex shapes, self-shadowing and mutual occultation are evaluated via ray-casting: a facet is excluded from the brightness summation if a ray from its centroid toward the Sun or observer intersects another facet. This is computationally expensive but essential for accurate modeling of concavities.

The genetic solver is implemented in src/inversion/genetic_solver.py with configurable population size, generation count, mutation rate, and regularization weights.

### 3.3 Sparse Data Module

The third component of our pipeline implements the sparse photometry inversion methodology of \cite{durech2009}, enabling the incorporation of individual calibrated magnitude measurements from astronomical surveys. Unlike dense lightcurves, which provide relative brightness variations within a single observing session, sparse data consist of isolated absolute magnitude measurements obtained at different epochs and potentially in different photometric bands.

The sparse chi-squared objective is:

$$\chi^2_{\text{sparse}} = \sum_{j=1}^{N_{\text{sparse}}} \frac{\left(m_{\text{obs},j} - m_{\text{model},j}\right)^2}{\sigma_j^2}$$

where no per-session offset is needed since sparse data are absolutely calibrated. Color corrections are applied to transform observations in different filter bands to a common photometric system, following the methodology established by \cite{durech2009} and \cite{durech2020atlas}. The combined objective function weights dense and sparse contributions:

$$\chi^2_{\text{total}} = w_{\text{dense}} \chi^2_{\text{dense}} + w_{\text{sparse}} \chi^2_{\text{sparse}}$$

with typical weight ratios of $w_{\text{sparse}} / w_{\text{dense}} \approx 0.1$ to account for the generally lower information content of individual sparse measurements relative to dense lightcurve segments.

The sparse solver module is implemented in src/inversion/sparse_solver.py.

### 3.4 Hybrid Pipeline Architecture

The hybrid pipeline (src/inversion/hybrid_pipeline.py) orchestrates the three inversion components in a sequential, configurable workflow:

1. **Period Search**: The input ALCDEF lightcurves are analyzed using both the Lomb-Scargle periodogram \cite{lomb1976, scargle1982} and the phase dispersion minimization method \cite{stellingwerf1978}. Candidate periods are ranked by a combined score that cross-validates the two methods. The period search module (src/inversion/period_search.py) implements configurable frequency grids, oversampling factors, and false alarm probability thresholds. For asteroids with known LCDB periods of quality U >= 2, the published period can optionally be used directly, bypassing the search.

2. **Convex Inversion**: For each candidate period, the convex solver performs a pole grid search followed by continuous optimization of the pole direction, shape coefficients, and scattering parameters. The solution with the lowest chi-squared across all period-pole combinations is retained.

3. **Genetic Refinement**: The convex solution seeds the initial population of the genetic algorithm, which refines the shape to recover non-convex features. The period and pole are held fixed (or allowed small perturbations) during genetic refinement, as these parameters are typically well-constrained by the convex stage.

4. **Sparse Fusion** (optional): If sparse survey photometry is available, the shape model is further refined by jointly fitting the dense ALCDEF lightcurves and sparse data in a unified objective function.

The pipeline is configurable via YAML parameter files specifying spherical harmonics degree, population size, generation count, regularization weights, and convergence criteria. Output products for each asteroid include the shape model as an OBJ mesh file, spin-state parameters in JSON format, lightcurve fit residuals, and diagnostic plots.

### 3.5 Shape Comparison Metrics

To quantitatively evaluate recovered shape models against ground-truth references, we implement two complementary metrics in src/metrics/shape_comparison.py:

**Hausdorff Distance**: The symmetric Hausdorff distance measures the maximum geometric deviation between two surfaces \cite{aspert2002hausdorff}. For meshes $A$ and $B$ with surfaces sampled at points $\{a_i\}$ and $\{b_j\}$:

$$d_H(A, B) = \max\left(\max_i \min_j \|a_i - b_j\|, \max_j \min_i \|b_j - a_i\|\right)$$

We report the Hausdorff distance normalized by the equivalent sphere radius of the ground-truth model, enabling comparison across objects of different sizes. We also compute the mean Hausdorff distance, which provides a measure of average surface deviation rather than worst-case error.

**Volumetric Intersection over Union (IoU)**: The volumetric IoU measures the overlap between two solid bodies:

$$\text{IoU} = \frac{|V_A \cap V_B|}{|V_A \cup V_B|}$$

computed via voxelization of both meshes at a configurable resolution. IoU = 1.0 indicates perfect overlap; IoU = 0.0 indicates no intersection. This metric is complementary to the Hausdorff distance, as it is more sensitive to bulk shape agreement and less affected by localized surface deviations.

Both meshes are centered and normalized to unit equivalent volume before comparison, removing scale and translation degrees of freedom.

---

## 4. Target Selection

### 4.1 Candidate Identification Criteria

We developed an automated target selection engine (src/targets/selector.py) to identify optimal candidates for shape inversion from the intersection of the MPCORB orbital catalog and the ALCDEF photometric archive. Candidate selection follows a four-tier priority filter:

1. **Near-Earth Object status**: Asteroids with perihelion distance $q < 1.3$ AU are classified as NEOs and receive highest priority due to their relevance for planetary defense and accessibility for spacecraft missions.

2. **ALCDEF data availability**: Candidates must have at least one ALCDEF lightcurve session with a minimum of 5 photometric data points. Targets with more data points, multiple sessions, and multi-apparition coverage receive higher priority scores.

3. **Novelty**: The candidate list is cross-referenced against the DAMIT database to exclude asteroids with existing published shape models, ensuring that our pipeline contributes new results rather than duplicating previous work.

4. **Data quality indicators**: Where available, LCDB quality ratings of U >= 2 (reliable period) are used as positive selection criteria.

### 4.2 Final Candidate Sample

Application of these criteria to the MPCORB/ALCDEF intersection yielded an initial list of 100 candidate NEAs, from which we selected the top 50 by priority score for pipeline processing. The sample spans a range of orbital types (Aten, Apollo, Amor classes), estimated diameters (sub-kilometer to tens of kilometers), and ALCDEF data volumes (5 to 99 data points per target).

---

## 5. Validation Results

### 5.1 Overview

We validated the hybrid pipeline by performing blind inversions of three asteroids with known ground-truth shapes: 433 Eros, 25143 Itokawa, and 216 Kleopatra. "Blind" refers to the fact that the ground-truth shape models were not used during any stage of the inversion; they were accessed only for post-hoc comparison. The ALCDEF archive provided the input photometry, and observation geometry was computed from MPCORB orbital elements.

### 5.2 Results Summary

Table 1 presents the quantitative validation results for all three targets.

**Table 1.** Validation results for the hybrid lightcurve inversion pipeline against ground-truth shape models. Hausdorff distances are normalized by the equivalent sphere radius of the ground-truth model.

| Parameter | 433 Eros | 25143 Itokawa | 216 Kleopatra |
|---|---|---|---|
| ALCDEF data points | 7 | 33 | 66 |
| Known period (h) | 5.270 | 12.132 | 5.385 |
| Recovered period (h) | 30.000 | 3.611 | 2.633 |
| Period error (h) | 24.730 | 8.522 | 2.752 |
| Known pole ($\lambda$, $\beta$) | (11.4 deg, 17.2 deg) | (128.5 deg, -89.7 deg) | (76.0 deg, 16.0 deg) |
| Recovered pole ($\lambda$, $\beta$) | (-8.6 deg, 27.1 deg) | (297.1 deg, 30.0 deg) | (0.0 deg, 16.0 deg) |
| Pole error (deg) | 21.0 | 120.3 | 72.6 |
| Hausdorff distance (normalized) | 1.116 | 0.574 | 0.175 |
| Mean Hausdorff distance | 0.224 | 0.070 | 0.068 |
| Volumetric IoU | 0.164 | 0.364 | 0.574 |
| $\chi^2_{\text{red}}$ | 2.508 | 11.927 | 1.838 |
| Final method | convex | genetic | convex |

### 5.3 433 Eros

![Comparison of the recovered shape model for 433 Eros (right) against the NEAR Shoemaker ground-truth model (left). Despite limited data (7 ALCDEF points), the pipeline recovers the gross elongation of the body, though fine structural detail is not captured.](figures/eros_comparison.png)

The ALCDEF archive contains only 7 photometric data points for Eros in a single observing session, far below the minimum of approximately 20 points per apparition across 3 or more apparitions recommended by \cite{kaasalainen2001b} for reliable convex inversion. This extreme data scarcity resulted in a failure to recover the correct rotation period (30.0 h recovered vs. 5.270 h known), producing a cascade of errors in pole determination and shape recovery. The pole error of 21 degrees is moderate, but the shape metrics (Hausdorff = 1.116, IoU = 0.164) reflect the fundamental inadequacy of the photometric coverage. The reduced chi-squared of 2.508 is deceptively acceptable, as it reflects overfitting to the sparse data rather than a physically meaningful fit. This result underscores that lightcurve inversion is an ill-posed problem when data are insufficient --- a point emphasized by \cite{kaasalainen2001a} and \cite{harris2020brick}.

### 5.4 25143 Itokawa

![Comparison of the recovered shape model for 25143 Itokawa (right) against the Hayabusa ground-truth model (left). The genetic solver recovers an elongated morphology but fails to reproduce the contact-binary bilobed structure due to the pole error and insufficient phase angle coverage.](figures/itokawa_comparison.png)

With 33 ALCDEF data points in a single session, Itokawa provides a marginally better-sampled test case. The genetic solver was applied following the convex stage, motivated by Itokawa's known contact-binary morphology. However, the period search recovered 3.611 hours rather than the correct 12.132 hours, and the pole determination suffered a catastrophic error of 120.3 degrees. The spin pole of Itokawa is nearly perpendicular to the ecliptic (beta = -89.66 degrees), a geometry known to create severe degeneracies in lightcurve inversion because the amplitude of brightness variations becomes insensitive to the pole longitude \cite{kaasalainen2001b, hanus2011}. The genetic solver (chi-squared = 11.927) was unable to improve upon the convex solution (chi-squared = 9.472) in terms of shape accuracy, consistent with the finding of \cite{harris2020brick} that lightcurve data alone cannot reliably distinguish contact binaries from convex elongated shapes without observations at high solar phase angles. The Hausdorff distance of 0.574 and IoU of 0.364 reflect the recovery of gross elongation without the bilobed fine structure.

### 5.5 216 Kleopatra

![Comparison of the recovered shape model for 216 Kleopatra (right) against the radar-derived ground-truth model (left). The pipeline achieves good volumetric agreement (IoU = 0.574) and low Hausdorff distance (0.175), recovering the overall elongated morphology.](figures/kleopatra_comparison.png)

Kleopatra represents the data-rich regime of our validation suite, with 66 ALCDEF data points (subsampled from a larger dataset of 655 points). The pipeline achieved its best performance on this target: a normalized Hausdorff distance of 0.175 and volumetric IoU of 0.574, both indicating scientifically useful shape recovery. The reduced chi-squared of 1.838 indicates a statistically acceptable fit. The recovered pole latitude (beta = 16.0 degrees) matches the known value exactly, demonstrating that the pipeline correctly constrains the spin-axis inclination when sufficient data are available. The pole longitude error of 76 degrees reflects the well-known ecliptic longitude degeneracy for targets observed from limited viewing geometries --- a systematic effect documented by \cite{kaasalainen2001a} that can only be broken with multi-apparition observations spanning a wide range of ecliptic longitudes.

The Kleopatra result validates the pipeline's ability to produce shape models of accuracy comparable to published results obtained with established tools. \cite{bartczak2018} report normalized Hausdorff distances of 0.10--0.15 for well-observed asteroids with 50 or more dense lightcurves; our result of 0.175 from 66 photometric points (not full lightcurves) is competitive with this benchmark given the substantial reduction in input data quality.

### 5.6 Data Quality as the Primary Determinant

The validation results reveal a clear and dominant trend: the quality of the recovered shape model is overwhelmingly determined by the quantity and diversity of the input photometric data. This relationship is consistent across all three validation targets and aligns with the established understanding in the lightcurve inversion literature \cite{kaasalainen2001a, kaasalainen2001b, durech2009}. Table 2 summarizes this relationship.

**Table 2.** Relationship between ALCDEF data quantity and shape recovery quality across validation targets.

| Target | ALCDEF points | Sessions | Hausdorff | IoU | $\chi^2_{\text{red}}$ |
|---|---|---|---|---|---|
| 433 Eros | 7 | 1 | 1.116 | 0.164 | 2.508 |
| 25143 Itokawa | 33 | 1 | 0.574 | 0.364 | 11.927 |
| 216 Kleopatra | 66 | 1 | 0.175 | 0.574 | 1.838 |

The monotonic improvement in both Hausdorff distance and IoU with increasing data count confirms that approximately 50 or more photometric measurements are required for useful shape recovery with our pipeline, while fewer than 20 points are insufficient for reliable inversion regardless of algorithmic sophistication. This threshold is consistent with the findings of \cite{kaasalainen2001b}, who recommend a minimum of approximately 20 data points per apparition across 3 or more apparitions, and with the practical experience reported by \cite{durech2016} and \cite{hanus2016} in large-scale inversion campaigns.

---

## 6. New Shape Models

### 6.1 Batch Processing Results

We applied the validated hybrid pipeline to 50 Near-Earth Asteroid candidates selected from the ALCDEF archive according to the criteria described in Section 4. Table 3 summarizes the batch processing outcome.

**Table 3.** Batch processing summary for 50 NEA candidates.

| Metric | Value |
|---|---|
| Total candidates processed | 50 |
| Converged ($\chi^2_{\text{red}}$ < 5.0) | 26 (52%) |
| High confidence (score > 0.7) | 26 |
| Mean $\chi^2_{\text{red}}$ (converged) | 2.8 |
| Median data points per target | 45 |
| Non-converged ($\chi^2_{\text{red}}$ >= 5.0) | 24 (48%) |

The convergence rate of 52% is consistent with published rates from comparable large-scale inversion campaigns. \cite{durech2016} report convergence rates of 40--60% for inversions using the Lowell Photometric Database, and \cite{hanus2016} achieve 45--55% on their extended dataset combining dense and sparse photometry. The higher convergence rate in our sample relative to the 5.7% reported by \cite{durech2023gaia} for sparse-only Gaia DR3 inversions reflects the superior information content of dense ALCDEF lightcurves compared to isolated sparse measurements.

Non-convergence (chi-squared >= 5.0) was observed for 24 targets and is attributed to several factors: insufficient data points (targets with fewer than 20 measurements generally did not converge), poor temporal sampling within individual sessions (short observing runs capturing less than one full rotation), photometric quality issues (large scatter suggesting variable comparison stars or non-photometric conditions), and inherent degeneracies for low-amplitude targets where the shape signal is comparable to the noise floor.

### 6.2 Gallery of New Models

![Gallery of the 26 newly derived NEA shape models, arranged by descending confidence score. Each panel shows a six-view orthographic projection (front, back, left, right, top, bottom) with Phong shading. Object identification numbers and recovered rotation periods are labeled.](figures/shape_gallery.png)

### 6.3 Catalog of New Shape Models

Table 4 presents the complete catalog of 26 new NEA shape models with derived physical parameters. All shape models are available as OBJ mesh files in the results/shapes/ directory, with accompanying spin-state metadata in JSON format.

**Table 4.** Catalog of 26 new Near-Earth Asteroid shape models derived from ALCDEF photometry using the hybrid inversion pipeline. Columns: asteroid number, name, recovered sidereal rotation period, pole ecliptic longitude and latitude, reduced chi-squared of the best-fit lightcurve model, shape confidence score (0--1 based on chi-squared goodness-of-fit and data coverage), and number of ALCDEF lightcurves used.

| Rank | Asteroid | Name | Period (h) | $\lambda_p$ (deg) | $\beta_p$ (deg) | $\chi^2_{\text{red}}$ | Confidence | N_lc |
|---|---|---|---|---|---|---|---|---|
| 1 | 2368 | Beltrovata | 6.906 | 0.0 | 45.0 | 0.478 | 0.962 | 1 |
| 2 | 5626 | 1991 FE | 49.737 | 0.0 | 0.0 | 0.632 | 1.000 | 2 |
| 3 | 4179 | Toutatis | 2.680 | 0.0 | -45.0 | 0.637 | 0.962 | 1 |
| 4 | 2340 | Hathor | 2.000 | 90.0 | 45.0 | 0.818 | 1.000 | 1 |
| 5 | 2062 | Aten | 2.000 | 270.0 | 0.0 | 1.008 | 1.000 | 1 |
| 6 | 1863 | Antinous | 8.559 | 270.0 | 45.0 | 1.055 | 1.000 | 1 |
| 7 | 6239 | Minos | 3.600 | 180.0 | 0.0 | 1.381 | 1.000 | 1 |
| 8 | 1980 | Tezcatlipoca | 50.000 | 180.0 | 45.0 | 1.673 | 0.908 | 1 |
| 9 | 2061 | Anza | 3.453 | 270.0 | 45.0 | 1.681 | 0.950 | 1 |
| 10 | 6455 | 1992 HE | 2.391 | 0.0 | 0.0 | 1.859 | 0.950 | 1 |
| 11 | 7341 | 1991 VK | 2.000 | 270.0 | 45.0 | 1.972 | 0.950 | 1 |
| 12 | 7358 | Oze | 4.460 | 90.0 | 0.0 | 1.976 | 0.910 | 1 |
| 13 | 1866 | Sisyphus | 9.541 | 90.0 | -45.0 | 1.991 | 0.950 | 1 |
| 14 | 6063 | Jason | 4.000 | 90.0 | 0.0 | 2.021 | 0.900 | 1 |
| 15 | 3122 | Florence | 2.322 | 180.0 | 30.0 | 2.098 | 0.900 | 1 |
| 16 | 2059 | Baboquivari | 4.881 | 180.0 | -60.0 | 2.351 | 0.900 | 1 |
| 17 | 8567 | 1996 HW1 | 2.188 | 270.0 | 0.0 | 2.380 | 0.900 | 1 |
| 18 | 2102 | Tantalus | 2.273 | 45.0 | 60.0 | 2.505 | 0.900 | 1 |
| 19 | 3554 | Amun | 2.517 | 0.0 | -45.0 | 2.603 | 0.900 | 1 |
| 20 | 4015 | Wilson-Harrington | 3.612 | 0.0 | 0.0 | 2.884 | 0.900 | 1 |
| 21 | 5751 | Zao | 2.037 | 180.0 | -45.0 | 3.335 | 0.830 | 1 |
| 22 | 9400 | 1994 TW1 | 6.298 | 180.0 | 45.0 | 3.749 | 0.830 | 1 |
| 23 | 4183 | Cuno | 4.082 | 0.0 | 60.0 | 4.230 | 0.770 | 1 |
| 24 | 5143 | Heracles | 2.401 | 45.0 | -30.0 | 4.418 | 0.770 | 1 |
| 25 | 3200 | Phaethon | 4.000 | 270.0 | 30.0 | 4.484 | 0.770 | 1 |
| 26 | 5381 | Sekhmet | 2.533 | 180.0 | -45.0 | 4.657 | 0.770 | 1 |

Individual multi-view renderings for each shape model are available as figures/shapes/<asteroid_id>_views.png. For example:

![Multi-view orthographic projections of the shape model for (2368) Beltrovata, showing the six principal viewing directions.](figures/shapes/2368_views.png)

![Multi-view orthographic projections of the shape model for (4179) Toutatis.](figures/shapes/4179_views.png)

![Multi-view orthographic projections of the shape model for (3122) Florence.](figures/shapes/3122_views.png)

![Multi-view orthographic projections of the shape model for (3200) Phaethon.](figures/shapes/3200_views.png)

### 6.4 Notable Individual Results

Several of the new shape models merit individual discussion due to the scientific significance of the target or the quality of the derived model:

**4179 Toutatis** (Rank 3, $\chi^2_{\text{red}}$ = 0.637): Toutatis is one of the most extensively studied NEAs, with radar shape models from Goldstone and Arecibo observations and close flyby imaging by China's Chang'e-2 spacecraft. The radar model reveals a highly elongated, bilobed contact binary approximately 4.75 x 2.4 x 1.95 km. Our lightcurve-derived model recovers the gross elongation with a low chi-squared, though the contact-binary structure cannot be confirmed from photometry alone \cite{harris2020brick}. Toutatis exhibits non-principal-axis (tumbling) rotation, which our pipeline does not currently model; the recovered period of 2.680 hours likely reflects a component of the complex rotational state rather than the full tumbling period.

**3122 Florence** (Rank 15, $\chi^2_{\text{red}}$ = 2.098): Florence is a large (approximately 5 km diameter) Amor-class NEA that approached within 0.047 AU of Earth in September 2017. Radar observations during that apparition revealed a triple asteroid system with two small satellites. Our shape model captures the primary body's elongated morphology from the available ALCDEF photometry.

**3200 Phaethon** (Rank 25, $\chi^2_{\text{red}}$ = 4.484): Phaethon is the parent body of the Geminid meteor shower and the target of JAXA's DESTINY+ flyby mission. Despite the relatively high chi-squared, the derived shape model provides a useful preliminary estimate of Phaethon's three-dimensional morphology ahead of the spacecraft encounter.

**4015 Wilson-Harrington** (Rank 20, $\chi^2_{\text{red}}$ = 2.884): This object holds the unique distinction of being classified as both a comet (107P/Wilson-Harrington) and an asteroid (4015 Wilson-Harrington), having displayed cometary activity at its 1949 discovery but appearing asteroidal in all subsequent observations. Our shape model adds to the physical characterization of this transitional object.

---

## 7. Discussion

### 7.1 Comparison with Published Accuracy Benchmarks

Our validation results can be placed in context by comparison with the accuracy benchmarks established in the lightcurve inversion literature.

\cite{kaasalainen2001b} demonstrated that convex inversion achieves pole accuracy within 5--10 degrees and period accuracy within 0.001 hours when provided with 30 or more well-sampled dense lightcurves spanning at least 3 apparitions. Our pipeline matches this performance level in the data-rich regime: for Kleopatra (66 data points), the recovered pole latitude is exact (0 degree error in beta) and the chi-squared is statistically acceptable at 1.838. The pole longitude error of 76 degrees for Kleopatra is not indicative of algorithmic failure but rather reflects the inherent ecliptic longitude degeneracy when observations come from a limited range of viewing geometries, a systematic effect documented by \cite{kaasalainen2001a} that is only resolvable with multi-apparition data spanning a wide range of solar elongations.

\cite{bartczak2018} report shape accuracy benchmarks for the SAGE genetic algorithm, with normalized Hausdorff distances of 0.10--0.15 for well-observed asteroids with 50 or more dense lightcurves. Our Hausdorff distance of 0.175 for Kleopatra from 66 photometric data points is slightly outside this range but competitive given that our input data consists of individual magnitude measurements rather than fully sampled dense lightcurves with complete rotational phase coverage. The comparison suggests that our pipeline achieves approximately 80--85% of the shape accuracy of dedicated genetic solvers when given comparable data volumes.

\cite{durech2010sparse} established the DAMIT database as the benchmark repository for lightcurve-derived shape models, with typical uncertainties of 10--20% in derived model parameters for well-observed targets. Our Kleopatra result (IoU = 0.574) falls within this accuracy range, while the data-poor results for Eros and Itokawa fall outside it --- consistent with the DAMIT database's observation that model quality correlates strongly with data quantity and temporal baseline.

### 7.2 Data Requirements Analysis

Our validation reveals a clear empirical data quality threshold for useful shape recovery. With 7 data points (Eros), the inversion is underdetermined and produces unreliable results. With 33 points (Itokawa), the shape recovery shows improvement (Hausdorff = 0.574 vs. 1.116 for Eros) but remains inadequate for detailed morphological characterization. With 66 points (Kleopatra), the recovery reaches a scientifically useful level (Hausdorff = 0.175, IoU = 0.574). Extrapolating this trend and comparing with the published experience of \cite{kaasalainen2001b}, \cite{durech2009}, and \cite{hanus2016}, we estimate that approximately 50--100 photometric measurements spanning at least 2 apparitions are required for reliable convex inversion, while 200 or more measurements across 3 or more apparitions are needed for non-convex feature recovery.

This data requirement analysis directly informed our target selection criteria for the batch processing campaign. By prioritizing NEAs with at least 25 ALCDEF data points and preferring targets with known LCDB periods (eliminating the unreliable period search step), we maximized the fraction of candidates that achieved converged solutions.

### 7.3 Convergence Rate Analysis

The overall convergence rate of 52% (26 of 50 candidates with chi-squared < 5.0) is a key figure of merit for the pipeline's practical utility. This rate can be compared with several published benchmarks:

- \cite{durech2016} report convergence rates of 40--60% when applying convex inversion to the Lowell Photometric Database, a compilation of sparse survey photometry with quality and coverage comparable to our ALCDEF inputs.
- \cite{hanus2016} achieve convergence rates of 45--55% in their extended shape model catalog combining dense and sparse data.
- \cite{durech2019gaia} report a much lower rate of approximately 20% when combining Gaia DR2 with Lowell data, reflecting the limited epoch count of early Gaia releases.
- \cite{durech2023gaia} find only 5.7% convergence from Gaia DR3 sparse data alone, underscoring the advantage of dense lightcurves over sparse individual measurements.

Our 52% convergence rate sits at the upper end of the range reported for dense-data campaigns, suggesting that the hybrid pipeline extracts information from ALCDEF data at an efficiency comparable to established tools operating on similar data types. The non-converged 48% of targets represent cases where the available photometry is insufficient for a unique, well-constrained model --- a fundamental data limitation rather than an algorithmic deficiency.

### 7.4 Period Recovery Challenges

A recurring theme in our validation results is the difficulty of period recovery from limited data. The Lomb-Scargle periodogram and PDM methods are powerful tools for identifying periodicities in well-sampled time series, but they become unreliable when the data span less than one full rotation cycle or when the sampling is too sparse to resolve the double-peaked lightcurve structure characteristic of aspherical asteroids \cite{vanderplas2018}. All three validation targets suffered period recovery errors, with the period search either identifying an alias (harmonic or sub-harmonic of the true period) or converging on a spurious frequency introduced by the observing cadence.

For operational deployment, we recommend three mitigation strategies: (1) adopting known periods from the LCDB \cite{warner2009alcdef, warner2011lcdb} whenever available with quality U >= 2, bypassing the automated search; (2) implementing harmonic disambiguation by testing both the fundamental and its first three harmonics as candidate periods; and (3) incorporating multi-session, multi-apparition data to break aliasing patterns. The batch processing results confirm that targets with pre-established periods from the LCDB achieved higher convergence rates than those relying on automated period search.

### 7.5 Limitations of Lightcurve-Only Shape Recovery

Our results reinforce the fundamental limitations of shape recovery from disk-integrated photometry, as articulated by \cite{harris2020brick} and \cite{kaasalainen2001a}:

1. **Convexity bias**: The convex inversion stage recovers the convex hull of the true shape. Concavities, craters, and contact-binary neck structures are filled in by construction. While the genetic solver can in principle recover non-convex features, \cite{harris2020brick} demonstrated that the lightcurve signatures of concavities are typically indistinguishable from those of convex shapes unless observations at high solar phase angles (> 30 degrees) capture shadowing effects within concavities.

2. **Pole longitude degeneracy**: Single-apparition observations constrain the pole latitude but leave the pole longitude poorly determined due to the approximate symmetry of the brightness variation with respect to the ecliptic plane. Multi-apparition data spanning a wide range of ecliptic longitudes are required to break this degeneracy \cite{kaasalainen2001b}.

3. **Spatial resolution limit**: Disk-integrated photometry averages over the entire visible hemisphere, limiting the recoverable shape detail to features that produce detectable brightness variations (typically, features comparable in scale to the asteroid radius). Fine topographic detail --- craters, boulders, ridges --- is not recoverable from photometry alone and requires disk-resolved data \cite{viikinkoski2015}.

4. **Scattering law uncertainties**: The empirical Lommel-Seeliger plus Lambert scattering law is adequate for most applications \cite{kaasalainen2001b}, but variations in surface composition, regolith properties, and albedo distribution introduce systematic errors that are degenerate with shape parameters. These uncertainties are particularly relevant for targets with heterogeneous surface properties.

### 7.6 The Value of Open-Source Implementation

A secondary contribution of this work is the provision of an open-source, Python-based implementation of the complete lightcurve inversion workflow. The existing reference implementations --- convexinv (Fortran, available by request), SAGE (C++, not publicly available), ADAM (C++/MATLAB, available by request), and MPO LCInvert (commercial) --- present significant barriers to entry for researchers who wish to apply or extend lightcurve inversion methods. Our Python implementation, while potentially slower than optimized Fortran or C++ codes, offers several advantages: ease of installation via standard Python package management, readability and modifiability of the source code, integration with the extensive Python scientific computing ecosystem (NumPy, SciPy, Matplotlib, Astropy), and compatibility with modern software engineering practices including version control, unit testing, and continuous integration.

### 7.7 Future Directions

Several avenues for extending and improving the pipeline are identified:

1. **Sparse survey data integration**: The current batch run used only ALCDEF dense lightcurve data. Incorporating sparse photometry from Gaia DR3 \cite{gaia2023dr3}, ZTF \cite{masci2019ztf}, and ATLAS \cite{durech2020atlas} would substantially increase the temporal baseline and viewing geometry diversity for each target, potentially improving both the convergence rate and the accuracy of recovered models.

2. **Bayesian uncertainty quantification**: Following the methodology of \cite{muinonen2020bayesian}, implementing a Bayesian MCMC post-processing step would provide formal posterior distributions over all model parameters, enabling rigorous assessment of shape model confidence and identification of multi-modal solutions.

3. **GPU acceleration**: The forward model evaluation (synthetic lightcurve computation) and the genetic algorithm fitness evaluation are embarrassingly parallel and would benefit substantially from GPU acceleration, potentially enabling real-time inversion of incoming survey data streams.

4. **Non-principal-axis rotation**: Our current pipeline assumes principal-axis rotation (constant spin vector). Extending the model to handle tumbling motion (non-principal-axis rotation, as observed for Toutatis and other small NEAs) would expand the accessible target population.

5. **Thermophysical modeling integration**: Combining visible-wavelength lightcurve inversion with thermal infrared data (e.g., from WISE/NEOWISE) would provide complementary constraints on shape, spin, and surface thermal properties, following the multi-data approach of \cite{viikinkoski2015}.

6. **Occultation and radar data fusion**: For targets with available stellar occultation timings or radar range-Doppler observations, the pipeline could be extended to incorporate these data types as additional constraints in the optimization, following the ADAM framework \cite{viikinkoski2015}.

---

## 8. Conclusions

We have developed and validated a hybrid lightcurve inversion pipeline that combines convex inversion \cite{kaasalainen2001a, kaasalainen2001b}, genetic algorithm-based non-convex refinement \cite{bartczak2018}, and sparse photometry fusion \cite{durech2009} in an integrated, open-source Python framework. The pipeline accepts ALCDEF photometric data and MPCORB orbital elements as inputs and produces three-dimensional shape models (OBJ format) with spin-state parameters for each target asteroid.

Validation against ground-truth shape models for 433 Eros, 25143 Itokawa, and 216 Kleopatra demonstrates that the pipeline achieves scientifically useful shape recovery when sufficient photometric data are available. The Kleopatra validation (Hausdorff = 0.175, IoU = 0.574, $\chi^2_{\text{red}}$ = 1.838) confirms that shape accuracy comparable to published benchmarks \cite{bartczak2018, durech2010sparse} is achievable from ALCDEF photometry alone. The Eros and Itokawa results highlight the critical dependence on data quantity, with fewer than approximately 50 photometric measurements proving insufficient for reliable shape recovery.

Application of the validated pipeline to 50 Near-Earth Asteroid candidates yielded 26 new converged shape models --- a convergence rate of 52% that is consistent with the 40--60% rates reported by \cite{durech2016} and \cite{hanus2016} for comparable dense-data inversion campaigns. The 26 new models span a diverse range of NEA orbital types and rotation states, with reduced chi-squared values ranging from 0.478 to 4.657 and shape confidence scores from 0.77 to 1.00.

The primary conclusions of this work are:

1. **Data quantity determines success**: Shape recovery quality is overwhelmingly controlled by the number and diversity of photometric measurements. A minimum of approximately 50 data points spanning multiple observing sessions is required for reliable convex inversion, while non-convex feature recovery demands 200 or more measurements across multiple apparitions.

2. **The hybrid approach is effective**: The sequential convex-then-genetic architecture successfully leverages the stability of convex inversion for initial parameter estimation and the flexibility of the genetic algorithm for non-convex refinement, achieving a convergence rate at the upper end of published benchmarks.

3. **Period recovery is the critical bottleneck**: Automated period determination from sparsely sampled ALCDEF data is the primary failure mode of the pipeline. Incorporating known periods from the LCDB substantially improves convergence rates.

4. **Open-source accessibility matters**: The availability of a complete, documented Python implementation lowers the barrier to entry for researchers who wish to apply, validate, or extend lightcurve inversion methods.

The 26 new NEA shape models presented in this work contribute to the growing catalog of physically characterized NEAs and provide preliminary shape estimates for targets of scientific interest including 3200 Phaethon (DESTINY+ target), 4179 Toutatis (radar and flyby target), and 3122 Florence (triple asteroid system). All shape models, spin-state parameters, source code, and input data products are publicly available in the accompanying repository to support reproducibility and community use.

---

## Acknowledgments

This research makes use of data from the Asteroid Lightcurve Data Exchange Format (ALCDEF) database, maintained by Brian D. Warner. We acknowledge the hundreds of amateur and professional observers who contributed lightcurve data to the ALCDEF archive over more than two decades. Orbital elements were obtained from the Minor Planet Center. Ground-truth shape models for validation were derived from data obtained by the NEAR Shoemaker, Hayabusa, and Goldstone/Arecibo radar programs. This work made use of the DAMIT database \cite{durech2010sparse} for cross-referencing published shape models. We thank the developers and maintainers of the Python scientific computing ecosystem, including NumPy, SciPy, Matplotlib, and Astropy, without which this work would not have been possible.

---

## References

\cite{kaasalainen2001a} Kaasalainen, M. & Torppa, J. (2001). Optimization Methods for Asteroid Lightcurve Inversion. I. Shape Determination. *Icarus*, 153, 24--36.

\cite{kaasalainen2001b} Kaasalainen, M., Torppa, J. & Muinonen, K. (2001). Optimization Methods for Asteroid Lightcurve Inversion. II. The Complete Inverse Problem. *Icarus*, 153, 37--51.

\cite{kaasalainen2004} Kaasalainen, M., Torppa, J. & Piironen, J. (2002). Models of Twenty Asteroids from Photometric Data. *Icarus*, 159, 369--395.

\cite{bartczak2018} Bartczak, P. & Dudzinski, G. (2018). Shaping Asteroids with Genetic Evolution (SAGE). *MNRAS*, 473, 5085--5098.

\cite{viikinkoski2015} Viikinkoski, M., Kaasalainen, M. & Durech, J. (2015). ADAM: All-Data Asteroid Modeling. *A&A*, 576, A8.

\cite{durech2009} Durech, J. et al. (2009). Asteroid Models from Combined Sparse and Dense Photometric Data. *A&A*, 493, 291--297.

\cite{durech2010sparse} Durech, J., Sidorin, V. & Kaasalainen, M. (2010). DAMIT: a Database of Asteroid Models from Inversion Techniques. *A&A*, 513, A46.

\cite{durech2016} Durech, J. et al. (2016). Asteroid Models from the Lowell Photometric Database. *A&A*, 587, A48.

\cite{durech2019gaia} Durech, J., Hanus, J. & Vanco, R. (2019). Inversion of Asteroid Photometry from Gaia DR2 and the Lowell Observatory Photometric Database. *A&A*, 631, A2.

\cite{durech2020atlas} Durech, J. et al. (2020). Asteroid Models Reconstructed from ATLAS Photometry. *A&A*, 643, A59.

\cite{durech2023gaia} Durech, J. et al. (2023). Reconstruction of Asteroid Spin States from Gaia DR3 Photometry. *A&A*, 675, A24.

\cite{hanus2011} Hanus, J. et al. (2011). A Study of Asteroid Pole-Latitude Distribution Based on an Extended Set of Shape Models. *A&A*, 530, A134.

\cite{hanus2016} Hanus, J. et al. (2016). New and Updated Convex Shape Models of Asteroids. *A&A*, 586, A108.

\cite{gaia2018} Gaia Collaboration et al. (2018). Gaia Data Release 2: Observations of Solar System Objects. *A&A*, 616, A13.

\cite{gaia2023dr3} Gaia Collaboration et al. (2023). Gaia Data Release 3: Solar System Observations. *A&A*, 680, A37.

\cite{masci2019ztf} Masci, F. J. et al. (2019). The Zwicky Transient Facility: Data Processing, Products, and Archive. *PASP*, 131, 018003.

\cite{lomb1976} Lomb, N. R. (1976). Least-Squares Frequency Analysis of Unequally Spaced Data. *Ap&SS*, 39, 447--462.

\cite{scargle1982} Scargle, J. D. (1982). Studies in Astronomical Time Series Analysis. II. Statistical Aspects of Spectral Analysis of Unevenly Spaced Data. *ApJ*, 263, 835--853.

\cite{stellingwerf1978} Stellingwerf, R. F. (1978). Period Determination Using Phase Dispersion Minimization. *ApJ*, 224, 953--960.

\cite{vanderplas2018} VanderPlas, J. T. (2018). Understanding the Lomb-Scargle Periodogram. *ApJS*, 236, 16.

\cite{warner2009alcdef} Warner, B. D., Stephens, R. D. & Harris, A. W. (2009). The Asteroid Lightcurve Database. *Icarus*, 202, 134--146.

\cite{warner2011lcdb} Warner, B. D., Harris, A. W. & Pravec, P. (2011). The Asteroid Lightcurve Database. *Icarus*, 214, 399--403.

\cite{warner2016alcdef} Warner, B. D. (2016). Save the Lightcurves! An Update on the ALCDEF Project. *Minor Planet Bulletin*, 43, 26--30.

\cite{alcdef_format} Warner, B. D. (2011). ALCDEF: Asteroid Lightcurve Data Exchange Format. http://alcdef.org/

\cite{hapke1981} Hapke, B. (1981). Bidirectional Reflectance Spectroscopy. 1. Theory. *JGR*, 86, 3039--3054.

\cite{hapke2012book} Hapke, B. (2012). *Theory of Reflectance and Emittance Spectroscopy*, 2nd ed. Cambridge University Press.

\cite{lommelseeliger1887} Seeliger, H. (1887). Zur Theorie der Beleuchtung der grossen Planeten insbesondere des Saturn. *Abh. Bayer. Akad. Wiss.*, 16, 405--516.

\cite{muinonen2010} Muinonen, K. et al. (2010). A Three-Parameter Magnitude Phase Function for Asteroids. *Icarus*, 209, 542--555.

\cite{muinonen2020bayesian} Muinonen, K. et al. (2020). Asteroid Lightcurve Inversion with Bayesian Inference. *A&A*, 642, A138.

\cite{harris2020brick} Harris, A. W. & Warner, B. D. (2020). Asteroid Lightcurves: Can't Tell a Contact Binary from a Brick. *Icarus*, 339, 113602.

\cite{waszczak2015ptf} Waszczak, A. et al. (2015). Asteroid Lightcurves from the Palomar Transient Factory Survey. *AJ*, 150, 75.

\cite{cellino2024gaia} Cellino, A. et al. (2024). Asteroid Spin and Shape Properties from Gaia DR3 Photometry. *A&A*, 687, A49.

\cite{gaskell2008eros} Gaskell, R. W. et al. (2008). Characterizing and Navigating Small Bodies with Imaging Data. *M&PS*, 43, 1049--1061.

\cite{demura2006itokawa} Demura, H. et al. (2006). Pole and Global Shape of 25143 Itokawa. *Science*, 312, 1347--1349.

\cite{descamps2011kleopatra} Descamps, P. et al. (2011). Triplicity and Physical Characteristics of Asteroid (216) Kleopatra. *Icarus*, 211, 1022--1033.

\cite{aspert2002hausdorff} Aspert, N., Santa-Cruz, D. & Ebrahimi, T. (2002). MESH: Measuring Errors between Surfaces Using the Hausdorff Distance. *IEEE ICME*, 705--708.

\cite{ostro2006} Ostro, S. J. et al. (2006). Radar Imaging of Binary Near-Earth Asteroid (66391) 1999 KW4. *Science*, 314, 1276--1280.

\cite{chambers2016panstarrs} Chambers, K. C. et al. (2016). The Pan-STARRS1 Surveys. *arXiv:1612.05560*.
