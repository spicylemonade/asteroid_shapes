# Literature Review: Asteroid Light Curve Inversion Methods

## 1. Convex Inversion (Kaasalainen & Torppa, 2001)

### Method
The foundational convex inversion method (Kaasalainen2001a, Kaasalainen2001b) parameterizes asteroid shapes as convex polyhedra defined by spherical harmonics coefficients or triangulated surface facets. The shape is constrained to be convex, which guarantees mathematical uniqueness.

### Key Equations
- **Forward model**: For a convex shape with facets of area $A_i$ and normal $\hat{n}_i$, the total scattered flux is:
  $$F = \sum_{i \in \text{visible}} A_i \cdot S(\mu_0^i, \mu^i, \alpha)$$
  where $\mu_0 = \cos(i)$ (incidence), $\mu = \cos(e)$ (emission), $\alpha$ = phase angle
- **Scattering law** (Lommel-Seeliger + Lambert):
  $$S(\mu_0, \mu, \alpha) = f(\alpha) \left[ c_{\rm LS} \frac{\mu_0}{\mu_0 + \mu} + c_{\rm L} \mu_0 \right]$$
- **Optimization**: Levenberg-Marquardt minimization of:
  $$\chi^2 = \sum_j \frac{(F_{\rm obs}^j - F_{\rm model}^j)^2}{\sigma_j^2}$$
- **Spin search**: Grid search over ecliptic longitude $\lambda \in [0°, 360°]$ and latitude $\beta \in [-90°, 90°]$ with ~5° steps, refined with gradient descent.

### Algorithmic Steps
1. Parameterize shape as spherical harmonics (order ~6-8) or ~500+ facet polyhedron
2. For each trial spin axis ($\lambda$, $\beta$) and period $P$:
   a. Compute viewing/illumination geometry for each observation epoch
   b. Render model brightness using LS+Lambert scattering
   c. Optimize shape parameters to minimize $\chi^2$
3. Select global minimum across spin grid

### Limitations
- Cannot resolve concavities (produces convex hull approximation)
- Multiple local minima in period/spin space
- Requires dense multi-apparition data for reliable solutions

---

## 2. SAGE: Genetic Evolution for Non-Convex Shapes (Bartczak & Dudziński, 2018)

### Method
SAGE (Bartczak2018) uses a genetic algorithm to evolve non-convex asteroid shape models. The genome encodes vertex positions of a triangulated mesh, allowing arbitrary concavities.

### Key Principles
- **Genome**: Vertex radii/displacements from a seed mesh (convex or spherical)
- **Fitness function**: $\chi^2$ between observed and modeled lightcurves
- **Operators**: Tournament selection, two-point crossover, Gaussian mutation
- **Forward model**: Must include ray-tracing for self-shadowing on non-convex shapes

### Algorithmic Steps
1. Initialize population of ~100 individuals (random perturbations of seed shape)
2. For each individual: compute all modeled lightcurves via ray-tracing
3. Evaluate fitness ($\chi^2$)
4. Apply selection, crossover, mutation
5. Repeat for ~500-1000 generations until convergence
6. Output best-fit non-convex shape

### Limitations
- Computationally expensive (ray-tracing per fitness evaluation)
- May converge to local minima without sufficient population diversity
- Requires good period determination as input

---

## 3. Sparse Data Inversion (Ďurech et al., 2009, 2010)

### Method
Sparse photometric inversion (Durech2009, Durech2010) extends convex inversion to handle survey-like data with few measurements per night, spread over multiple apparitions.

### Key Adaptations
- Uses **absolute (calibrated) magnitudes** rather than relative lightcurves
- Coarse spin-axis grid search with convex shape fitting
- Combines sparse survey data with dense lightcurves when available

### Challenges
- Period ambiguity is severe with sparse data
- Shape models from sparse-only data are crude (angular, unrealistic)
- Dense data supplement dramatically improves shape quality

### Success Criteria
- Typically requires >100 sparse points across ≥3 apparitions
- High-amplitude asteroids converge more reliably
- Pole accuracy ~10-30° with sparse-only; ~5-10° with sparse+dense

---

## 4. ADAM: All-Data Asteroid Modeling (Viikinkoski et al., 2015)

### Method
ADAM (Viikinkoski2015) is a unified framework for fusing multiple data types: lightcurves, adaptive optics images, stellar occultation timings, radar delay-Doppler, and thermal IR.

### Key Innovation
- All disk-resolved data types handled via 2D Fourier transform of projected shape
- Non-convex and non-starlike shapes supported
- Weighting scheme balances contributions from different data modalities

### Implementation
- Core in MATLAB with C subroutines for performance
- Available on GitHub (github.com/matvii/ADAM)

---

## 5. Self-Shadowing Ray-Tracing

### Physical Motivation
For non-convex shapes, facets can cast shadows on other facets. The Lommel-Seeliger + Lambert model assumes all visible facets are illuminated, which is incorrect for concave regions (craters, bifurcated contact binaries).

### Implementation Approach
- For each visible facet, cast a ray toward the Sun direction
- Test for intersection with other mesh facets
- Facets in shadow contribute zero direct illumination
- Use BVH (Bounding Volume Hierarchy) for O(log N) intersection testing

### Impact
- At high phase angles, shadowing effects become pronounced (Durech2003)
- Essential for modeling contact binaries and deeply cratered bodies
- Brightness difference can exceed 1-5% for elongated/concave shapes

---

## 6. Scattering Laws

### Lommel-Seeliger (LS)
$$S_{\rm LS} = \frac{\omega}{8\pi} \frac{\mu_0}{\mu_0 + \mu} f(\alpha)$$
- Single-scattering model, valid for dark (low albedo) surfaces
- Applicable to C-type asteroids (Muinonen2015)

### Lambert
$$S_{\rm L} = A_L \mu_0$$
- Ideal diffuse surface, valid for high-albedo objects
- Rarely used alone for asteroids

### Combined LS + Lambert
$$S = f(\alpha) \left[ c_{\rm LS} \frac{\mu_0}{\mu_0 + \mu} + c_{\rm L} \mu_0 \right]$$
- Standard for lightcurve inversion (Kaasalainen2001a)
- Weight $c_{\rm LS}$ typically 0.5-0.9

### Hapke Model (Hapke1993, Hapke2012)
- Multi-parameter physically motivated model
- Parameters: single-scattering albedo $\omega$, asymmetry factor $g$, roughness $\bar{\theta}$, opposition effect (SHOE, CBOE)
- More accurate but computationally expensive and parameter-degenerate

---

## 7. Survey of Existing Tools

| Tool | Language | Non-Convex | Sparse | Self-Shadowing | License |
|------|----------|-----------|--------|---------------|---------|
| MPO LCInvert | C++ | No | No | No | Commercial |
| SAGE | Fortran | Yes | No | Yes (ray-tracing) | Not publicly released |
| KOALA | Fortran/C | Partial | No | No | Not public |
| ADAM | MATLAB/C | Yes | Partial | Via rendering | Public (GitHub) |
| convexinv | Fortran | No | Yes | No | Public (DAMIT) |
| LSST/Rubin | Python | No | Yes (planned) | No | Open source |

### Gap Analysis
Our pipeline must fill these gaps:
1. **No existing tool combines**: non-convex shapes + sparse data + self-shadowing in one integrated pipeline
2. **SAGE**: powerful but not publicly available and not designed for sparse data
3. **ADAM**: handles multiple data types but focused on resolved data fusion
4. **convexinv**: mature sparse handling but convex-only
5. **No tool** implements the "convex seed → GA non-convex refinement" workflow with self-shadowing
