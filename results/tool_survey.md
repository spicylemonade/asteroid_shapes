# Survey of Existing Lightcurve Inversion Tools and Software

## 1. convexinv (Kaasalainen Reference Code)

- **Language**: Fortran 77/90
- **License**: Academic/restricted - available by request from M. Kaasalainen and J. Durech
- **Algorithm**: Convex inversion using spherical harmonics parameterization with Levenberg-Marquardt optimization. Implements the complete forward model from Kaasalainen & Torppa (2001).
- **Input**: Dense lightcurve files in custom format (epoch, intensity, geometry)
- **Output**: Shape coefficients (spherical harmonics), pole direction, period
- **Known Limitations**: Convex shapes only; cannot recover concavities. Requires dense well-sampled lightcurves. Not publicly available as a standalone package.
- **Source Accessible**: By request only; code basis for DAMIT models
- **URL**: https://astro.troja.mff.cuni.cz/projects/damit/

## 2. SAGE (Shaping Asteroids with Genetic Evolution)

- **Language**: C++
- **License**: Not publicly released
- **Algorithm**: Genetic/evolutionary algorithm using vertex-mesh representation allowing non-convex shapes. Population-based optimization with crossover, mutation, elitism. Bartczak & Dudzinski (2018).
- **Input**: Dense lightcurves with observation geometry
- **Output**: Non-convex triangulated mesh (.obj), spin vector
- **Known Limitations**: Computationally expensive (hours per asteroid). Requires good period estimate. No public source code release. Sensitive to initial population.
- **Source Accessible**: No - proprietary code by Bartczak group (Poznan)
- **URL**: Paper only (MNRAS 473, 5085-5098)

## 3. KOALA (Knitted Occultation, Adaptive Optics, and Lightcurve Analysis)

- **Language**: Fortran / IDL
- **License**: Academic, limited distribution
- **Algorithm**: Multi-data inversion combining lightcurves, stellar occultation chords, and adaptive optics images. Uses MPFIT (Levenberg-Marquardt) for optimization.
- **Input**: Lightcurves, occultation timings, AO contours
- **Output**: Non-convex shape model, spin vector
- **Known Limitations**: Requires multiple data types for best results; occultation and AO data rarely available. Not publicly distributed.
- **Source Accessible**: No public repository
- **URL**: Carry et al. (2010, 2012)

## 4. ADAM (All-Data Asteroid Modeling)

- **Language**: C++ with Python wrapper
- **License**: Academic, available by request
- **Algorithm**: Multi-resolution shape optimization using octanoid-based subdivision surfaces. Simultaneously fits lightcurves, AO images, occultation silhouettes, and radar data. Viikinkoski et al. (2015).
- **Input**: Multiple data types (lightcurves, images, occultations)
- **Output**: High-resolution non-convex shape model
- **Known Limitations**: Requires diverse data types for full capability. Complex setup. Not widely distributed.
- **Source Accessible**: Limited - by request from Viikinkoski/Kaasalainen group
- **URL**: A&A 576, A8 (2015)

## 5. MPO LCInvert

- **Language**: Compiled binary (likely Delphi/C++)
- **License**: Commercial software ($)
- **Algorithm**: Convex inversion based on Kaasalainen method. GUI-based workflow for amateur astronomers. Includes period search (Fourier analysis).
- **Input**: ALCDEF format lightcurves or custom format
- **Output**: Convex shape model, pole solution, period
- **Known Limitations**: Convex only. Closed source. Windows-only. Expensive for institutional use. Limited automation capability.
- **Source Accessible**: No - commercial closed-source
- **URL**: http://www.minorplanetobserver.com/MPOSoftware/MPOLCInvert.htm

## 6. Durech's Lightcurve Inversion Codes

- **Language**: Fortran 90
- **License**: Academic, available through DAMIT project
- **Algorithm**: Multiple programs for convex inversion with dense and sparse data. Implements the methodology from Durech et al. (2009, 2010, 2016).
- **Input**: Dense lightcurves, sparse survey photometry (Lowell, Gaia)
- **Output**: Convex shape + spin via spherical harmonics
- **Known Limitations**: Convex only. Sparse-only inversion requires large photometric datasets. Fortran codebase difficult to modify.
- **Source Accessible**: Some codes available through DAMIT website and by request
- **URL**: https://astro.troja.mff.cuni.cz/projects/damit/

## 7. Asteroid Lightcurve Photometry Analysis (ALPA) / otras herramientas Python

- **Language**: Python
- **License**: Various open-source
- **Algorithm**: Various period-finding and basic analysis tools. No full inversion capability.
- **Known Limitations**: No shape inversion; period analysis only
- **Source Accessible**: Yes (GitHub)

## 8. PyLightcurve / Minor Planet tools

- **Language**: Python
- **License**: Open source (MIT/BSD)
- **Algorithm**: Light curve utilities, period search (Lomb-Scargle, PDM). No shape inversion.
- **Known Limitations**: Analysis tools only, no inversion
- **Source Accessible**: Yes

## Summary Table

| Tool | Language | License | Convex | Non-Convex | Sparse | Source Available |
|------|----------|---------|--------|------------|--------|-----------------|
| convexinv | Fortran | Academic | Yes | No | Yes | By request |
| SAGE | C++ | Proprietary | No | Yes | No | No |
| KOALA | Fortran/IDL | Academic | Yes | Yes | No | No |
| ADAM | C++ | Academic | Yes | Yes | No | By request |
| MPO LCInvert | Binary | Commercial | Yes | No | No | No |
| Durech codes | Fortran | Academic | Yes | No | Yes | Partial |

## Gap Analysis

Key gaps that our pipeline addresses:

1. **No open-source Python implementation** exists for asteroid lightcurve inversion
2. **No tool combines** convex + non-convex + sparse data in a single pipeline
3. **SAGE** (non-convex) is not publicly available
4. **Sparse inversion** capability is limited to Durech's Fortran codes
5. **Automated batch processing** is not supported by any existing tool
6. Our hybrid pipeline will be the **first open-source Python tool** to combine all three approaches (convex, genetic non-convex, sparse) with automated target selection and batch execution.
