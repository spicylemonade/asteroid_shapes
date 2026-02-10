# Repository Inventory

## 1. ALCDEF_ALL.zip

- **File size**: 130 MB (compressed)
- **Total files**: 24,643 lightcurve data files
- **Naming convention**: `ALCDEF_<number>_<name>.txt` where `<number>` is the asteroid number (0 for unnumbered objects using provisional designation) and `<name>` is the asteroid name or provisional designation (e.g., `ALCDEF_433_Eros.txt`, `ALCDEF_0_1994_CJ1.txt`)
- **File format**: Plain text with ALCDEF (Asteroid Lightcurve Data Exchange Format) structure

### Sample Record Format (all ALCDEF fields observed)

Each file contains one or more lightcurve session blocks with the following structure:

**Metadata block** (between `STARTMETADATA` and `ENDMETADATA`):

| Field | Description | Example |
|-------|-------------|---------|
| `UPDATED` | Date of last update | `2018-12-07 00:00:00` |
| `REVISEDDATA` | Whether data has been revised | `FALSE` |
| `ALLOWSHARING` | Data sharing permission | `TRUE` |
| `SUBMITPDS` | Submitted to PDS | `TRUE` |
| `LCBLOCKID` | Unique lightcurve block identifier | `16365` |
| `OBJECTNUMBER` | Asteroid number (0 if unnumbered) | `0` |
| `OBJECTNAME` | Asteroid name or designation | `1994 CJ1` |
| `MPCDESIG` | MPC provisional designation | `1994 CJ1` |
| `SESSIONDATE` | Observation session date | `2014-06-27` |
| `SESSIONTIME` | Session start time (UTC) | `05:47:50` |
| `CONTACTNAME` | Observer contact name | `Brian D. Warner` |
| `CONTACTINFO` | Contact details | email + address |
| `OBSERVERS` | Observer name(s) | `B. D. Warner` |
| `FACILITY` | Observatory/facility name | full facility string |
| `MPCCODE` | MPC observatory code | `U82` |
| `TELESCOPE` | Telescope description | `0.50-m f/8.1 R-C` |
| `DETECTOR` | Detector/CCD type | `FLI-1001E CCD` |
| `OBSLONGITUDE` | Observatory longitude | `-116.38486` |
| `OBSLATITUDE` | Observatory latitude | `+34.27250` |
| `EXPOSURE` | Exposure time in seconds (-99 = unknown) | `-99` |
| `EXPJD` | When JD applies relative to exposure | `MID` |
| `OBJECTRA` | Right ascension of object | `14:17.9` |
| `OBJECTDEC` | Declination of object | `+15 16` |
| `EQUINOX` | Coordinate equinox | `J2000.0` |
| `PHASE` | Solar phase angle (degrees) | `+66.58` |
| `PABL` | Phase Angle Bisector longitude | `+241.0` |
| `PABB` | Phase Angle Bisector latitude | `+17.4` |
| `FILTER` | Photometric filter used | `C` (clear), `V`, `R`, etc. |
| `MAGBAND` | Magnitude band | `V` |
| `CICORRECTION` | Color index correction applied | `FALSE` |
| `CIBAND` | Color index band | (may be blank) |
| `DIFFERMAGS` | Whether differential magnitudes | `FALSE` |
| `STANDARD` | Calibration standard type | `INTERNAL` |
| `LTCAPP` | Light-time correction applied | `NONE` |
| `REDUCEDMAGS` | Reduced magnitude system | `NONE` |
| `COMMENT` | Free-text comments | (multiple allowed) |
| `DELIMITER` | Data column delimiter | `PIPE` |
| `COMPNAMEn` | Comparison star identifier | positional string |
| `COMPRAn` | Comparison star RA | `14:17:47.05` |
| `COMPDECn` | Comparison star Dec | `+15:19:40.1` |
| `COMPMAGn` | Comparison star magnitude | `+16.614` |
| `COMPMAGBANDn` | Comparison mag band | `V` |
| `COMPCIn` | Comparison color index | `+0.418` |
| `COMPCIBANDn` | Comparison CI band | `NONE` |

**Data block** (after `ENDMETADATA`):

Each data line has format: `DATA=<JD>|<magnitude>|<mag_error>`

- JD: Julian Date (mid-exposure)
- magnitude: Observed magnitude
- mag_error: Magnitude uncertainty

Example: `DATA=2456835.739711|18.582|0.091`

## 2. MPCORB.DAT.gz

- **File size**: 88 MB (compressed)
- **Total lines**: 1,512,845
- **Total orbit records**: 1,512,800 (after 43-line header + 2 separator lines)
- **Epoch**: K25BL (packed MPC format; corresponds to 2025)

### Column Definitions (MPC Orbit Format)

| Column Range | Field | Description |
|---|---|---|
| 1-7 | `Des'n` | Packed MPC designation (number or provisional) |
| 9-13 | `H` | Absolute magnitude |
| 15-19 | `G` | Slope parameter |
| 21-25 | `Epoch` | Epoch (packed MPC format) |
| 27-35 | `M` | Mean anomaly (degrees) |
| 37-46 | `Peri.` | Argument of perihelion (degrees) |
| 48-57 | `Node` | Longitude of ascending node (degrees) |
| 59-68 | `Incl.` | Inclination (degrees) |
| 70-79 | `e` | Eccentricity |
| 81-91 | `n` | Mean daily motion (degrees/day) |
| 93-103 | `a` | Semimajor axis (AU) |
| 106-106 | `U` | Uncertainty parameter |
| 108-116 | `Reference` | Reference identifier |
| 118-122 | `#Obs` | Number of observations |
| 124-126 | `#Opp` | Number of oppositions |
| 128-136 | `Arc` | Arc length (years or days) |
| 138-141 | `rms` | RMS residual (arcseconds) |
| 143-145 | `Perts` | Coarse perturbers indicator |
| 147-149 | `Perts` | Precise perturbers indicator |
| 151-160 | `Computer` | Orbit computer name |
| 162-165 | `Hex flags` | Hex digit flags |
| 167-194 | `Name/Desig` | Readable name or designation |
| 195-202 | `Last obs` | Date of last observation (YYYYMMDD) |

### Sample Records

```
00001    3.35  0.15 K25BL 231.53975   73.29974   80.24963   10.58789  0.0795763  0.21429712   2.7656157    (1) Ceres
00002    4.11  0.15 K25BL 211.52977  310.93340  172.88859   34.92833  0.2306430  0.21379713   2.7699258    (2) Pallas
00433   11.16  0.15 K25BL ...        ...        ...         ...       ...        ...          1.4583 ...   (433) Eros
```

## 3. Output Directories

| Directory | Contents |
|-----------|----------|
| `figures/` | Empty (0 files) - designated for plot output |
| `results/` | Empty (0 files) - designated for data/report output |

## 4. Source Code

**No existing source code found.** The repository contains only:

- `README.md` - Minimal project description ("# asteroid_shapes")
- `research_rubric.json` - Research plan/rubric with 27 items across 5 phases
- `TASK_researcher_attempt_1.md` - Task description document
- `ALCDEF_ALL.zip` - Lightcurve data archive (130 MB)
- `MPCORB.DAT.gz` - Orbital elements database (88 MB)

No `src/`, `tests/`, or Python/C++ source files exist. The entire inversion pipeline must be built from scratch.
