# Validation Report: 25143 Itokawa

## Summary

This report documents the validation test for asteroid **25143 Itokawa**, a small
rubble-pile near-Earth asteroid whose detailed shape model is known from the
Hayabusa mission. The purpose of this test is to evaluate how accurately the
lightcurve inversion pipeline can recover the known physical parameters from
publicly available ALCDEF photometric data.

## Data Used

| Parameter | Value |
|---|---|
| Asteroid ID | 25143 |
| Number of lightcurves | 1 |
| Number of data points | 33 |
| Inversion method | genetic |

The input data consisted of a single lightcurve containing **33 photometric data
points** sourced from the ALCDEF database. While this represents more data than
the Eros test case, it remains sparse for reliable shape inversion, particularly
for a highly irregular body like Itokawa.

## Period Search Results

| Parameter | Value |
|---|---|
| Known rotation period | 12.1324 h |
| Recovered rotation period | 3.6109 h |
| Period error | 8.5215 h |

The recovered period of 3.6109 h is significantly shorter than the known value of
12.1324 h, with an absolute error of approximately 8.52 h. The recovered period
is close to a 1/3 sub-harmonic of the true period, suggesting the algorithm
locked onto an aliased frequency. This is a common failure mode when rotational
phase coverage is incomplete.

## Shape Recovery Results

### Pole Direction

| Parameter | Known | Recovered | Error |
|---|---|---|---|
| Ecliptic longitude (lambda) | 128.5 deg | 297.1 deg | -- |
| Ecliptic latitude (beta) | -89.66 deg | 30.0 deg | -- |
| Angular separation | -- | -- | 120.3 deg |

The pole direction error of 120.3 deg is very large, indicating near-opposite
hemisphere placement of the recovered pole relative to the known solution. The
true pole of Itokawa is nearly aligned with the south ecliptic pole (beta =
-89.66 deg), a geometry that is particularly difficult to recover from sparse
single-apparition data.

### Shape Metrics

| Metric | Value |
|---|---|
| Hausdorff distance (normalized) | 0.5743 |
| Hausdorff distance (mean) | 0.0703 |
| Volumetric IoU | 0.3635 |
| Reduced chi-squared | 11.927 |

The IoU of 0.3635 indicates partial but insufficient shape recovery. The high
reduced chi-squared of 11.927 confirms a poor fit between the model and observed
lightcurves. The genetic algorithm method was employed here, but with only 33 data
points, it could not adequately explore the parameter space.

## Comparison Figure

![25143 Itokawa shape comparison](figures/itokawa_comparison.png)

## Discussion of Limitations

Itokawa presents a particularly challenging case for lightcurve inversion due to
several compounding factors:

- **Sparse ALCDEF data**: With only **33 data points** from a single lightcurve,
  the dataset falls well below the density needed for robust inversion. While
  better than the Eros case (7 points), it remains insufficient for a body with
  Itokawa's complexity.

- **Highly non-convex shape**: Itokawa's distinctive bi-lobed contact-binary
  morphology produces lightcurve features (shadowing, concavity effects) that
  cannot be captured by convex or low-resolution genetic algorithm models. The
  genetic method was selected automatically, but even this more flexible approach
  could not compensate for the data deficit.

- **Near-polar spin axis**: Itokawa's pole latitude of -89.66 deg means the spin
  axis is nearly perpendicular to the ecliptic plane. This geometry creates
  viewing conditions where single-apparition observations sample a very limited
  range of aspect angles, making pole determination especially difficult.

- **Period aliasing**: The recovered 3.6109 h period (roughly 1/3 of the true
  12.1324 h period) is a classic aliasing artifact arising from insufficient
  temporal coverage.

For comparison, the Kleopatra validation test achieved an IoU of 0.574 with 66
data points, demonstrating that the pipeline is capable of meaningful shape
recovery when more data is available.
