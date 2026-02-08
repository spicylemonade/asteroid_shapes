# Validation Report: 433 Eros

## Summary

This report documents the validation test for asteroid **433 Eros**, a well-studied
near-Earth asteroid whose true shape, rotation period, and pole orientation are known
from the NEAR Shoemaker mission. The purpose of this test is to evaluate how
accurately the lightcurve inversion pipeline can recover these known physical
parameters from publicly available ALCDEF photometric data.

## Data Used

| Parameter | Value |
|---|---|
| Asteroid ID | 433 |
| Number of lightcurves | 1 |
| Number of data points | 7 |
| Inversion method | convex |

The input data consisted of a single lightcurve containing only **7 photometric
data points** sourced from the ALCDEF database. This is an exceptionally sparse
dataset, far below what is typically required for reliable lightcurve inversion.

## Period Search Results

| Parameter | Value |
|---|---|
| Known rotation period | 5.27025 h |
| Recovered rotation period | 30.0 h |
| Period error | 24.7298 h |

The recovered period of 30.0 h deviates drastically from the known value of
5.27025 h, yielding an absolute error of approximately 24.73 h. With only 7 data
points spanning a single lightcurve, the period search landscape is severely
under-constrained, and the algorithm converged on an incorrect local minimum.

## Shape Recovery Results

### Pole Direction

| Parameter | Known | Recovered | Error |
|---|---|---|---|
| Ecliptic longitude (lambda) | 11.4 deg | -8.6 deg | -- |
| Ecliptic latitude (beta) | 17.2 deg | 27.1 deg | -- |
| Angular separation | -- | -- | 21.0 deg |

The pole direction error of 21.0 deg is moderate in isolation, but this result
should be interpreted cautiously given the gross period mismatch.

### Shape Metrics

| Metric | Value |
|---|---|
| Hausdorff distance (normalized) | 1.1161 |
| Hausdorff distance (mean) | 0.2242 |
| Volumetric IoU | 0.1639 |
| Reduced chi-squared | 2.508 |

The normalized Hausdorff distance of 1.1161 and IoU of 0.1639 indicate very poor
shape recovery. The reconstructed convex model bears little resemblance to the
true elongated shape of Eros.

## Comparison Figure

![433 Eros shape comparison](figures/eros_comparison.png)

## Discussion of Limitations

The extremely poor performance on Eros is directly attributable to the severe
data sparsity. With only **7 photometric data points** from a single lightcurve,
the inverse problem is drastically under-determined:

- **Period search failure**: Seven points cannot adequately sample even a single
  full rotation cycle of a 5.27 h period object. The period search algorithm
  lacks the frequency resolution needed and settled on a spurious 30.0 h alias.

- **Shape degeneracy**: Convex inversion requires dense phase-angle and
  rotational-phase coverage to constrain facet reflectances. With 7 points, the
  shape solution is essentially unconstrained, yielding an IoU of only 0.164.

- **ALCDEF data context**: The ALCDEF database entry for Eros used in this test
  contained a fragmentary observation set. Eros has extensive photometric data in
  the literature, but the specific ALCDEF records available for this automated
  pipeline run were minimal.

For comparison, the Kleopatra validation (66 data points) achieved an IoU of
0.574, demonstrating that the pipeline performs substantially better when
sufficient data is available. The Eros result should be regarded as a lower-bound
stress test rather than a representative measure of pipeline capability.
