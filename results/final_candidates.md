# Final Prioritized Candidate List

## Summary

Total candidates with converged shape models: **26**
High-confidence models (score > 0.7): **26**

## Top 10 Most Scientifically Interesting New Shape Models

### 1. (2368) Beltrovata

- **Type**: Near-Earth Asteroid
- **Period**: 6.9055 hours
- **Pole**: (λ=0.0°, β=45.0°)
- **χ² reduced**: 0.478
- **Confidence**: 0.962
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 2. (5626) 1991 FE

- **Type**: Near-Earth Asteroid
- **Period**: 49.7374 hours
- **Pole**: (λ=0.0°, β=0.0°)
- **χ² reduced**: 0.632
- **Confidence**: 1.0
- **Data**: 2 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 3. (4179) Toutatis

- **Type**: Near-Earth Asteroid
- **Period**: 2.6804 hours
- **Pole**: (λ=0.0°, β=-45.0°)
- **χ² reduced**: 0.637
- **Confidence**: 0.962
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 4. (2340) Hathor

- **Type**: Near-Earth Asteroid
- **Period**: 2.0 hours
- **Pole**: (λ=90.0°, β=45.0°)
- **χ² reduced**: 0.818
- **Confidence**: 1.0
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 5. (2062) Aten

- **Type**: Near-Earth Asteroid
- **Period**: 2.0 hours
- **Pole**: (λ=270.0°, β=0.0°)
- **χ² reduced**: 1.008
- **Confidence**: 1.0
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 6. (1863) Antinous

- **Type**: Near-Earth Asteroid
- **Period**: 8.5587 hours
- **Pole**: (λ=270.0°, β=45.0°)
- **χ² reduced**: 1.055
- **Confidence**: 1.0
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 7. (6239) Minos

- **Type**: Near-Earth Asteroid
- **Period**: 3.6004 hours
- **Pole**: (λ=180.0°, β=0.0°)
- **χ² reduced**: 1.381
- **Confidence**: 1.0
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 8. (1980) Tezcatlipoca

- **Type**: Near-Earth Asteroid
- **Period**: 50.0 hours
- **Pole**: (λ=180.0°, β=45.0°)
- **χ² reduced**: 1.673
- **Confidence**: 0.908
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 9. (2061) Anza

- **Type**: Near-Earth Asteroid
- **Period**: 3.4525 hours
- **Pole**: (λ=270.0°, β=45.0°)
- **χ² reduced**: 1.681
- **Confidence**: 0.95
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

### 10. (6455) 1992 HE

- **Type**: Near-Earth Asteroid
- **Period**: 2.3907 hours
- **Pole**: (λ=0.0°, β=0.0°)
- **χ² reduced**: 1.859
- **Confidence**: 0.95
- **Data**: 1 lightcurves

This NEA's shape model represents one of the first lightcurve-derived 3D shapes for this object, filling a gap in our knowledge of the NEA population's physical properties.

## Methodology

All shape models were derived using our hybrid inversion pipeline:
1. Period determination via Lomb-Scargle + PDM
2. Convex inversion using spherical harmonics (l_max=3-6)
3. Non-convex refinement via genetic algorithm (where applicable)
4. Quality assessment via reduced chi-squared

Shape confidence scores are based on chi-squared goodness-of-fit and data coverage (number of lightcurves and temporal span).
