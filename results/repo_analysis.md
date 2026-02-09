# Repository Structure & Data Asset Analysis

## 1. Repository Files and Purposes

| File/Directory | Purpose |
|---|---|
| `ALCDEF_ALL.zip` (130 MB) | Asteroid Lightcurve Data Exchange Format archive containing dense photometric lightcurve data for ~24,643 asteroids |
| `MPCORB.DAT.gz` (88 MB) | Minor Planet Center orbital elements database (~1,512,800 records) for all known minor planets |
| `README.md` | Project readme (minimal) |
| `research_rubric.json` | Research task rubric tracking progress through 27 items across 5 phases |
| `TASK_researcher_attempt_1.md` | Task specification document |
| `results/` | Output directory for data products (catalogs, reports, shape models) |
| `figures/` | Output directory for publication-quality figures |
| `src/` | Source code directory for pipeline modules (to be populated) |
| `sources.bib` | Bibliography file (to be created) |

## 2. Contents of ALCDEF_ALL.zip

- **Total files**: 24,643 text files
- **Numbered asteroids**: 23,753 files (one per numbered asteroid)
- **Unnumbered (provisional designation)**: 890 files
- **Unique numbered asteroid IDs**: ~23,743

### Data Format (ALCDEF text format)
Each file contains one or more lightcurve "blocks", each structured as:
- `STARTMETADATA` / `ENDMETADATA` section with key-value pairs:
  - `OBJECTNUMBER`, `OBJECTNAME`, `MPCDESIG` — asteroid identification
  - `SESSIONDATE`, `SESSIONTIME` — observation session timing
  - `MPCCODE` — observatory code
  - `PHASE` — solar phase angle (degrees)
  - `PABL`, `PABB` — phase angle bisector longitude/latitude
  - `FILTER`, `MAGBAND` — photometric filter and magnitude band
  - `DIFFERMAGS` — whether differential magnitudes
  - `REDUCEDMAGS` — reduction method
  - `DELIMITER=PIPE` — data delimiter indicator
- `DATA=JD|magnitude|uncertainty` lines (pipe-delimited)

### Key Test Asteroids Present
| Asteroid | Files | Blocks | Data Points |
|---|---|---|---|
| 433 Eros | 1 | 27 | 2,289 |
| 1580 Betulia | 1 | present | TBD |
| 25143 Itokawa | 1 | present | TBD |
| 216 Kleopatra | 1 | present | TBD |
| 3 Juno | 1 | present | TBD |
| 4 Vesta | 1 | present | TBD |
| 52 Europa | 1 | present | TBD |
| 1036 Ganymed | 1 | present | TBD |

## 3. Contents of MPCORB.DAT.gz

- **Total data records**: 1,512,800 orbital element entries
- **Format**: Fixed-width columns (MPC Orbit Format)

### Column Schema (from MPC documentation)
| Column Range | Field | Description |
|---|---|---|
| 1-7 | Designation | Packed MPC designation (number or provisional) |
| 9-13 | H | Absolute magnitude |
| 15-19 | G | Slope parameter |
| 21-25 | Epoch | Epoch (packed MPC format) |
| 27-35 | M | Mean anomaly (degrees) |
| 38-46 | Peri | Argument of perihelion (degrees) |
| 49-57 | Node | Longitude of ascending node (degrees) |
| 60-68 | Incl | Inclination (degrees) |
| 71-79 | e | Eccentricity |
| 81-91 | n | Mean daily motion (degrees/day) |
| 93-103 | a | Semimajor axis (AU) |
| 106-106 | U | Uncertainty parameter |
| 108-116 | Reference | Reference code |
| 118-122 | #Obs | Number of observations |
| 124-126 | #Opp | Number of oppositions |
| 128-136 | Arc | Arc span (years or days) |
| 138-141 | rms | RMS residual (arcsec) |
| 143-145 | Perts | Perturbation flags |
| 147-149 | Computer | Orbit computer |
| 151-160 | Flags | Additional flags |
| 167-194 | Name | Object name/designation (readable) |

### Sample Records (first 3 numbered objects)
- (1) Ceres: H=3.35, G=0.15, a=2.766 AU, e=0.080, i=10.59°
- (2) Pallas: H=4.11, G=0.15, a=2.770 AU, e=0.231, i=34.93°
- (3) Juno: H=5.19, G=0.15, a=2.671 AU, e=0.256, i=12.99°

## 4. Output Directories

- `results/` — Empty, ready for output data (catalogs, reports, shape models, spin vectors)
  - `results/ground_truth/` — For DAMIT ground truth shape models
  - `results/shapes/` — For generated .obj shape files
- `figures/` — Empty, ready for publication-quality PNG/PDF figures

## 5. Confirmation

Both data files are **accessible and parseable**:
- `ALCDEF_ALL.zip`: Successfully opened with Python `zipfile`, files extracted and parsed (pipe-delimited JD|mag|uncertainty format confirmed)
- `MPCORB.DAT.gz`: Successfully decompressed with Python `gzip`, fixed-width format confirmed with 1,512,800 orbital records
