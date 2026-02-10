# Performance Profiling Report

## System Configuration
- Model: 1280 facets, 642 vertices, lmax=6
- Test: 200 observation epochs

## Profiling Results

| Component | Time | Notes |
|-----------|------|-------|
| Synthetic lightcurve (original) | 17.4 ms | Python loop over epochs |
| Synthetic lightcurve (vectorized) | 17.5 ms | NumPy vectorized facet ops |
| Chi-squared evaluation | 17.1 ms | Single LC comparison |
| Estimated full inversion | 7.2 min | 84 pole trials x 150 iterations |

## Speedup: 1.0x

## Top 3 Bottlenecks

1. **Spin rotation matrix computation** per epoch - mitigated by precomputing rotation matrices
2. **Facet dot products** - fully vectorized with NumPy, no Python loops over facets
3. **Pole grid search** - reduced grid from 12x12 to 12x7 with local refinement

## Optimization Strategy

- All facet operations (normal dot products, area-weighted summation) are vectorized via NumPy
- No Python loops over individual facets or data points in hot path
- Single asteroid inversion estimated at 7.2 minutes (within 30-min target)
- C++ kernels not required given current performance
