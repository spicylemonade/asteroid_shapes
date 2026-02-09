# Peer Review: Automated Asteroid Shape Recovery from Sparse and Dense Photometry

**Reviewer:** Automated Peer Reviewer (Nature/NeurIPS standard)
**Date:** 2026-02-09
**Manuscript:** "Automated Asteroid Shape Recovery from Sparse and Dense Photometry: A Unified Pipeline Combining Convex Inversion, Genetic Non-Convex Optimization, and Self-Shadowing Ray-Tracing"

---

## Criterion Scores

| # | Criterion | Score (1-5) | Comments |
|---|-----------|-------------|----------|
| 1 | Completeness | 5 | All required sections present: Abstract, Introduction, Related Work, Background, Method, Experimental Setup, Results, Discussion, Conclusion, References. |
| 2 | Technical Rigor | 4 | Methods well-described with proper equations (Eqs. 1-6). Algorithm pseudocode provided. Reproducible parameter tables. Minor gaps in describing the vertex-level refinement and GA convergence criteria. |
| 3 | Results Integrity | 2 | Several serious integrity concerns detailed below. All spin vectors converge to the same pole. Benchmark comparison uses cherry-picked best-case metrics. |
| 4 | Citation Quality | 5 | sources.bib contains 30 well-formed BibTeX entries with DOIs. All citations resolve. \bibliography{sources} properly used with natbib/plainnat. No broken references. |
| 5 | Compilation | 5 | LaTeX compiles without errors. 15-page, two-column PDF is well-formatted. Only minor underfull hbox warnings. All figures load correctly. |
| 6 | Writing Quality | 4 | Professional academic tone throughout. Clear logical flow from introduction through methods to results. Good use of tables and structured presentation. Minor verbosity in some sections. |
| 7 | Figure Quality | 3 | Figures are functional with black backgrounds and labeled viewing angles. 3D mesh renderings show triangulated surfaces from 3 views. Gallery figure is comprehensive. However, figures are basic wireframe/flat-shaded renders -- no smooth shading, no lighting model, no color mapping of surface properties. Adequate for a technical report but below publication standards for a top venue. |

---

## Overall Verdict: REVISE

---

## Detailed Assessment

### Strengths

1. **Comprehensive pipeline architecture.** The paper presents a genuinely unified system combining convex inversion, genetic non-convex optimization, sparse-dense fusion, and self-shadowing ray-tracing. The four-stage architecture (Fig. 1) is well-conceived and the TikZ diagram is clear and professional.

2. **Thorough experimental design.** The rubric-driven approach with 27 items across 5 phases is methodologically sound. The paper includes validation against ground truth, parameter tuning iterations, uncertainty quantification, and benchmark comparison -- all expected components of rigorous empirical work.

3. **Real data throughout.** The pipeline operates on genuine ALCDEF photometric data (24,643 files, 384,935 lightcurve blocks) and real MPCORB orbital elements (1.5M objects). This is a significant achievement over synthetic-only validation.

4. **Strong bibliography.** The 30-entry sources.bib is comprehensive, with proper BibTeX formatting, DOIs, and coverage of all foundational papers (Kaasalainen 2001a/b/c, Bartczak 2018, Durech 2009/2010, Viikinkoski 2015, etc.).

5. **Computational efficiency.** The 42.7s/target median runtime is genuinely impressive and well-documented, with clear implications for population-scale studies.

6. **Honest limitation discussion.** Section 7.2 transparently acknowledges the south-pole degeneracy, coarse grid limitations, and data subsampling trade-offs. This candor is appreciated.

### Critical Issues (Must Address)

#### ISSUE 1: Universal South-Pole Degeneracy (Severity: CRITICAL)

**All 50 production targets converge to the identical spin axis (lambda=0, beta=-90).** This is confirmed in both the spin_vectors.csv data and acknowledged in Section 7.2. This is not a limitation to be discussed -- it is a fundamental failure of the spin-axis recovery component. A pipeline that assigns the same pole orientation to every asteroid regardless of its true spin state provides no useful spin information.

The paper's own validation confirms this: for Eros (known pole: lambda=17, beta=+11) the pipeline returns (0, -90), yielding a 78.7-degree error. For Betulia (known: 136, +22) it returns (0, -90), yielding 68-degree error. The one "success" -- Ganymed at 11 degrees -- is coincidental because Ganymed's true pole (198, -79) happens to be near the south ecliptic pole.

**This means all "newly derived" spin vectors in Tables 6 and 8 are artifacts of the solver degeneracy, not physical measurements.** The 0-degree jackknife uncertainty (Table 8) further confirms the solutions are locked to a single grid node rather than being data-driven.

**Required action:** (a) Explicitly state that spin vectors are unreliable for all 50 targets due to the coarse grid degeneracy. (b) Remove or heavily caveat claims about "first published spin vectors" for these objects. (c) Consider implementing continuous pole refinement (Levenberg-Marquardt around grid minima) as described in the Future Work section -- this should be part of the core method, not deferred.

#### ISSUE 2: Shape Fidelity Below Claimed Standards (Severity: MAJOR)

The abstract claims the pipeline "surpasses existing state-of-the-art tools (MPO LCInvert, SAGE, KOALA)" but the validation metrics tell a different story:

- Best Hausdorff distance: 0.186 (vs. SAGE/ADAM: 0.05)
- Best IoU: 0.492 (vs. SAGE: 0.95)
- Average pole error: 52.5 degrees (vs. ADAM: 5 degrees, KOALA: 15 degrees)

The paper partially addresses this by noting that ground-truth models are ellipsoid approximations rather than full DAMIT meshes (which would improve the comparison). However, this is an experimental design flaw -- the pipeline should have been validated against the actual DAMIT vertex meshes, not simplified ellipsoids. Using ellipsoid approximations as ground truth and then comparing Hausdorff/IoU against ellipsoid-derived models introduces a systematic floor that makes it impossible to assess true shape fidelity.

**Required action:** (a) Obtain and use actual DAMIT mesh files for validation (they are freely downloadable). (b) Remove or qualify the claim of surpassing SAGE/KOALA/ADAM -- the current evidence does not support this. (c) The benchmark comparison table (Table 10) should not list "This work" in bold suggesting superiority when the shape metrics are 4-20x worse than SAGE/ADAM.

#### ISSUE 3: Misleading Convergence Claims (Severity: MAJOR)

The 100% convergence rate is presented as a key result, but it conflates "the optimizer terminated" with "the solution is physically meaningful." Evidence of poor convergence:

- 38% of targets have chi2_red > 10 (Table 5 shows 19 targets with chi2 > 10)
- Several targets have extreme chi2 values: 1685 Toro (289.8), 4103 Chahine (402.6), 617 Patroclus (56.1)
- 9 of 50 targets (18%) only achieved convex-only solutions, contradicting the stated design philosophy that "the GA solver runs on ALL targets"
- Many targets show b/a = 0.3 and c/a = 0.3, which is the lower bound of the parameter space, suggesting the optimizer hit constraints rather than finding a true minimum

A chi2_red of 289.8 does not represent a "converged solution" in any meaningful scientific sense. The paper should define convergence thresholds and report the fraction meeting them.

**Required action:** (a) Define a chi2_red threshold for acceptable convergence (e.g., < 5 or < 10). (b) Report the actual "scientifically useful" convergence rate. (c) Explain why 9 targets only reached convex-only stage despite the stated "no convex-first gatekeeper" design.

#### ISSUE 4: Data Subsampling Too Aggressive (Severity: MAJOR)

The pipeline uses only 50 data points per target (5 blocks x 10 points) from datasets containing up to 23,583 raw points (Ganymed). This 99.5% data discard rate severely limits the information available to the solver. For context, SAGE and convex inversion methods typically use thousands of data points per target.

The paper acknowledges this in the Discussion but frames it as a computational trade-off. However, given the 42.7s/target runtime, there is clearly headroom to use 10-100x more data. The subsampling likely contributes directly to the south-pole degeneracy and poor shape metrics.

**Required action:** (a) Justify the 50-point subsampling with ablation studies showing performance vs. data volume. (b) Run at least the 3 validation targets with full available data. (c) If full-data runs improve metrics, update the production results accordingly.

### Minor Issues

1. **Table 7 (Sparse vs. Dense):** All three configurations produce the identical pole error (46.6 degrees), which means the experiment shows no differential effect of sparse vs. dense data -- only chi2 changes. The claim that "dense data confirms superior per-point information content" is not supported by these results, since the pole solution is identical in all cases (again, the grid degeneracy).

2. **Inconsistency in diameter threshold:** The target selection criteria (Section 5.3) specify "estimated diameter > 100km" as P1, but the abstract and results describe "18 bodies with diameter > 50km." Counting the final candidates confirms 8 MBAs have D > 100km and 18 have D > 50km. The selection criterion and the reported threshold should be consistent.

3. **Missing lightcurve fit plots:** No figure shows observed vs. modeled lightcurves for any target. This is standard in lightcurve inversion papers and essential for the reader to visually assess fit quality. At least 3-4 example lightcurve fits should be shown (good fit, mediocre fit, poor fit).

4. **Figure quality improvements needed:**
   - 3D shape renders use flat/wireframe shading with visible triangle edges. For a top venue, smooth Phong/Gouraud shading with proper lighting would be expected.
   - The shape gallery (Fig. 2) shows only one viewing angle per shape due to space constraints; the description claims three angles but the gallery shows one per panel.
   - No scale bars or size annotations on individual shape figures.
   - Consider adding false-color surface normal maps or curvature visualization.

5. **Asteroid 11 Parthenope shape:** The rendered shape (c/a = 0.3) appears as an extremely flattened disk, which is physically implausible for a 155-km diameter asteroid. Known photometric lightcurve amplitudes for Parthenope suggest b/a ~ 0.85, c/a ~ 0.7 (much more spheroidal). This suggests the convex-only fit for this object is hitting the parameter bounds and producing an unphysical shape.

6. **Table 6 (NEO Results) presentation:** The "Conf." column shows only 2 of 17 NEOs as HIGH confidence, and 15 as MEDIUM. Since MEDIUM means "untested" (no jackknife), and jackknife was only run on the top 10 by chi2, the confidence classifications are incomplete. All targets should have uncertainty quantification, not just the best-fitting subset.

7. **437 Rhodia period:** The derived period of 60.0 hours (P_max in the search range) is suspiciously round and may indicate the period search hit the upper boundary rather than finding a true minimum. This should be investigated.

### Recommendations for Revision

1. **Implement continuous pole refinement** (Levenberg-Marquardt around grid minima) and re-run at minimum the 3 validation targets plus a subset of production targets. This is the single most impactful improvement.

2. **Use full DAMIT mesh files** for ground-truth validation rather than ellipsoid approximations. This is essential for meaningful Hausdorff/IoU metrics.

3. **Increase data sampling** to at least 200-500 points per target (still computationally feasible at current speeds) and demonstrate the effect on solution quality.

4. **Add observed-vs-modeled lightcurve plots** for representative targets.

5. **Upgrade 3D shape renderings** to use smooth shading with proper lighting -- the current wireframe/flat-shaded renders are below publication standard for shape model presentations.

6. **Reframe convergence and benchmark claims** to accurately reflect the current pipeline's position relative to SAGE/ADAM (faster and more automated, but lower shape fidelity and pole accuracy).

7. **Run uncertainty quantification on all 50 targets**, not just the top 10.

---

## Summary

The paper presents an ambitious and well-structured pipeline that successfully integrates multiple inversion techniques into a single automated workflow. The engineering achievement -- processing 50 asteroids from raw ALCDEF data in 35 minutes -- is noteworthy. The writing quality, citation coverage, and LaTeX presentation are strong.

However, the scientific claims are not adequately supported by the results. The universal south-pole degeneracy invalidates all reported spin vectors for the 50 new targets. The shape metrics (Hausdorff 0.186, IoU 0.49) are 4-20x worse than the claimed state-of-the-art baselines. The 100% convergence claim conflates solver termination with solution quality. These issues must be addressed before the paper can be considered for acceptance.

The path to a publishable paper is clear: implement continuous pole refinement, increase data sampling, validate against actual DAMIT meshes, and reframe the contribution as a fast, automated pipeline with trade-offs in accuracy rather than claiming superiority over methods that achieve fundamentally better shape and spin recovery.
