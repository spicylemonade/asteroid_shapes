# Ground Truth Shape Models

This directory contains reference shape models for validation asteroids.
Models are approximated from published spacecraft and radar-derived dimensions.

## 433 Eros

- **Source**: NEAR Shoemaker mission; Gaskell et al. (2008) \cite{gaskell2008eros}
- **Dimensions**: 34.4 x 11.2 x 11.2 km
- **Shape**: Elongated, saddle-like concavity
- **Spin**: lambda=11.4 deg, beta=17.2 deg, P=5.27025 h
- **Files**: `433_eros.obj`, `433_eros_spin.json`

## 25143 Itokawa

- **Source**: Hayabusa mission; Demura et al. (2006) \cite{demura2006itokawa}
- **Dimensions**: 535 x 294 x 209 m
- **Shape**: Contact binary (bilobed, "sea otter")
- **Spin**: lambda=128.5 deg, beta=-89.66 deg, P=12.1324 h
- **Files**: `25143_itokawa.obj`, `25143_itokawa_spin.json`

## 216 Kleopatra

- **Source**: Radar observations; Descamps et al. (2011) \cite{descamps2011kleopatra}
- **Dimensions**: 217 x 94 x 81 km
- **Shape**: Dog-bone / dumbbell, metallic M-type
- **Spin**: lambda=76 deg, beta=16 deg, P=5.385 h
- **Files**: `216_kleopatra.obj`, `216_kleopatra_spin.json`

## Notes

Shape models are generated as icosphere-based meshes (2562 vertices, 5120 faces)
with dimensions and morphology approximating the published radar/spacecraft models.
These serve as ground-truth references for validation of our lightcurve inversion pipeline.
