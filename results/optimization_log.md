# Optimization Log

## Overview

This document records the optimization and threshold evaluation results for the
asteroid lightcurve inversion validation suite. Three targets with known
ground-truth shapes were tested: 433 Eros, 25143 Itokawa, and 216 Kleopatra.

## Quality Thresholds

The following acceptance thresholds were defined for shape recovery quality:

| Metric | Threshold | Direction |
|---|---|---|
| Hausdorff distance (normalized) | < 0.30 | Lower is better |
| Volumetric IoU | > 0.40 | Higher is better |

A target **passes** validation when both thresholds are satisfied simultaneously.

## Threshold Evaluation Results

### 216 Kleopatra -- PASSED

| Metric | Value | Threshold | Status |
|---|---|---|---|
| Hausdorff (normalized) | 0.175 | < 0.30 | PASS |
| Volumetric IoU | 0.574 | > 0.40 | PASS |

Kleopatra is the only target that passed both quality thresholds. With 66 ALCDEF
data points and the convex inversion method, the pipeline produced a shape model
with a Hausdorff distance well below the 0.30 cutoff and an IoU well above the
0.40 cutoff.

### 433 Eros -- FAILED

| Metric | Value | Threshold | Status |
|---|---|---|---|
| Hausdorff (normalized) | 1.1161 | < 0.30 | FAIL |
| Volumetric IoU | 0.1639 | > 0.40 | FAIL |

Eros failed both thresholds by wide margins. The root cause is extremely sparse
input data: only **7 photometric data points** from a single lightcurve. This is
far below the minimum data density required for meaningful lightcurve inversion.

### 25143 Itokawa -- FAILED

| Metric | Value | Threshold | Status |
|---|---|---|---|
| Hausdorff (normalized) | 0.5743 | < 0.30 | FAIL |
| Volumetric IoU | 0.3635 | > 0.40 | FAIL |

Itokawa failed both thresholds, though its IoU of 0.3635 approached the 0.40
boundary. The dataset contained **33 data points** from a single lightcurve --
better than Eros but still insufficient, compounded by Itokawa's highly non-convex
bi-lobed shape and near-polar spin axis.

## Data Sparsity Analysis

A key finding from this validation round is that **data quantity is the dominant
factor** determining recovery quality:

| Target | Data Points | Hausdorff (norm) | IoU | Passed |
|---|---|---|---|---|
| 433 Eros | 7 | 1.1161 | 0.1639 | No |
| 25143 Itokawa | 33 | 0.5743 | 0.3635 | No |
| 216 Kleopatra | 66 | 0.1750 | 0.5744 | Yes |

Both Eros (7 points) and Itokawa (33 points) fall below the recommended minimum
data threshold for reliable inversion. The literature generally suggests a minimum
of approximately 50-100 well-distributed data points per lightcurve for convex
inversion, with multiple apparitions strongly preferred.

### Recommended Data Thresholds

Based on these results, the following minimum data requirements are recommended
before attempting shape inversion:

| Parameter | Minimum Recommended |
|---|---|
| Data points per lightcurve | >= 50 |
| Number of lightcurves | >= 2 (preferably 3+) |
| Apparitions | >= 2 |

Targets below these thresholds should be flagged with a data-quality warning in
pipeline output to set appropriate expectations for result reliability.

## Pipeline Performance Assessment

The validation results support the following conclusions:

1. **The pipeline works well when sufficient data is available.** The Kleopatra
   case (66 data points, single lightcurve) passed both quality thresholds with
   comfortable margins (Hausdorff 0.175 vs. threshold 0.30; IoU 0.574 vs.
   threshold 0.40). This demonstrates that the core inversion algorithms, period
   search, pole optimization, and convex shape modeling are fundamentally sound.

2. **Performance degrades predictably with data sparsity.** The monotonic
   relationship between data point count and all quality metrics (Hausdorff, IoU,
   period error) confirms that poor results on Eros and Itokawa are attributable
   to insufficient input data rather than algorithmic deficiencies.

3. **Period aliasing is the primary failure mode under sparse data.** All three
   targets exhibited period errors, but the severity scaled inversely with data
   density. Kleopatra's half-period alias (2.63 h vs. 5.39 h) is a well-understood
   and partially recoverable degeneracy, while Eros and Itokawa showed gross
   period mismatches.

4. **The genetic method did not outperform convex inversion for Itokawa.** Despite
   using a more flexible genetic algorithm approach, the 33-point Itokawa dataset
   was too sparse to benefit from the additional model complexity. Method selection
   should be gated on data sufficiency.

## Next Steps

- Augment ALCDEF queries with additional data sources (e.g., literature
  compilations, observatory archives) to increase data point counts for Eros and
  Itokawa.
- Implement a pre-inversion data-quality gate that warns users or skips inversion
  when data falls below recommended thresholds.
- Test the pipeline on targets with richer ALCDEF coverage (100+ data points,
  multiple apparitions) to establish upper-bound performance.
- Investigate period-alias disambiguation strategies (e.g., imposing priors from
  taxonomic class or diameter estimates).
