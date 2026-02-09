"""
Fetch and store ground truth asteroid shape models from the DAMIT database
(Database of Asteroid Models from Inversion Techniques).

This module provides:
  1. Hardcoded ground truth spin parameters for well-studied asteroids
     whose lightcurve data also exists in our ALCDEF dataset.
  2. Ellipsoid-approximation .obj model generation for shape validation.
  3. JSON export of ground truth parameters to results/ground_truth/.

Ground truth sources:
  - DAMIT: https://astro.troja.mff.cuni.cz/projects/damit/
    (Durech et al. 2010, A&A 513, A46)
  - NEAR Shoemaker mission (Miller et al. 2002) for 433 Eros
  - Radar + lightcurve inversions (Kaasalainen et al. 2001, 2002, 2004)
  - 3D Asteroid Catalogue: https://3d-asteroids.space/

References in sources.bib:
  - Durech2010, Kaasalainen2001a, Kaasalainen2001b
"""

import json
import os
import sys
import numpy as np
from datetime import datetime

# ---------------------------------------------------------------------------
# Repository paths
# ---------------------------------------------------------------------------
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GROUND_TRUTH_DIR = os.path.join(REPO_ROOT, "results", "ground_truth")
SRC_DIR = os.path.join(REPO_ROOT, "src")

sys.path.insert(0, SRC_DIR)
from mesh_utils import save_obj


# ============================================================
# Ground truth data from DAMIT and spacecraft missions
# ============================================================

# Each entry documents published spin parameters and shape information.
# Fields:
#   asteroid_number: IAU number
#   name: common name
#   damit_model_id: DAMIT internal model ID (if known)
#   spin_lambda_deg: ecliptic longitude of spin pole (J2000), degrees
#   spin_beta_deg:   ecliptic latitude of spin pole (J2000), degrees
#   period_hours:    sidereal rotation period, hours
#   period_uncertainty_hours: uncertainty on period
#   pole_uncertainty_deg: approximate pole direction uncertainty, degrees
#   axis_ratios:     [a/a, b/a, c/a] triaxial ellipsoid ratios (a >= b >= c)
#   dimensions_km:   [a, b, c] in km if known from spacecraft/radar
#   n_vertices:      number of vertices in DAMIT convex model (typical)
#   n_faces:         number of faces in DAMIT convex model (typical)
#   source:          primary reference for the values
#   notes:           additional context
#   in_alcdef:       whether this asteroid has lightcurve data in our ALCDEF_ALL.zip

DAMIT_GROUND_TRUTH = {
    433: {
        "asteroid_number": 433,
        "name": "Eros",
        "damit_model_id": None,  # Eros has spacecraft model; DAMIT may list LC inversion model
        "spin_lambda_deg": 17.22,
        "spin_beta_deg": 11.35,
        "period_hours": 5.27025547,
        "period_uncertainty_hours": 1e-8,
        "pole_uncertainty_deg": 0.01,
        "axis_ratios": [1.0, 0.33, 0.33],  # 34.4 x 11.2 x 11.3 km -> b/a~0.33
        "dimensions_km": [34.4, 11.2, 11.3],
        "equivalent_diameter_km": 16.84,
        "n_vertices": 7790,  # Gaskell high-res model
        "n_faces": 15572,
        "spectral_type": "S",
        "albedo": 0.25,
        "source": "Miller et al. 2002, Icarus 155, 3; NEAR Shoemaker",
        "notes": (
            "Ground truth from NEAR Shoemaker spacecraft. "
            "Pole: RA=11.3695 deg, Dec=17.2273 deg (J2000 equatorial). "
            "Ecliptic: lambda=17.22, beta=11.35. "
            "Lightcurve inversion by Durech et al. recovers lambda~17, beta~8. "
            "Shape model from Gaskell (PDS, 2008). "
            "DAMIT convex LC-inversion model may have ~500-1000 vertices."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 27,
        "alcdef_datapoints": 2289,
    },
    1036: {
        "asteroid_number": 1036,
        "name": "Ganymed",
        "damit_model_id": 1849,
        "spin_lambda_deg": 198.0,
        "spin_beta_deg": -79.0,
        "period_hours": 10.31304,
        "period_uncertainty_hours": 0.00005,
        "pole_uncertainty_deg": 5.0,
        "axis_ratios": [1.0, 0.976, 0.929],  # 42 x 41 x 39 km
        "dimensions_km": [42.0, 41.0, 39.0],
        "equivalent_diameter_km": 37.675,
        "n_vertices": 1022,  # typical DAMIT convex model
        "n_faces": 2040,
        "spectral_type": "S",
        "albedo": 0.238,
        "source": "Hanus et al. 2016; DAMIT Model 1849",
        "notes": (
            "DAMIT Model 1849. Multiple pole solutions in literature: "
            "(214,-73), (190,-78), (198,-79). Kaasalainen et al. (2001): "
            "lambda=208, beta=-76. Adopted value from Hanus et al. "
            "Retrograde rotator. Largest NEO (Amor group)."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 134,
        "alcdef_datapoints": 23583,
    },
    15: {
        "asteroid_number": 15,
        "name": "Eunomia",
        "damit_model_id": 1838,
        "spin_lambda_deg": 0.0,
        "spin_beta_deg": -68.0,
        "period_hours": 6.082752,
        "period_uncertainty_hours": 0.000002,
        "pole_uncertainty_deg": 5.0,
        "axis_ratios": [1.0, 0.85, 0.75],  # approximate from lightcurve amplitude
        "dimensions_km": [357.0, 255.0, 212.0],  # Vernazza et al. 2021 (SPHERE)
        "equivalent_diameter_km": 268.0,
        "n_vertices": 1022,
        "n_faces": 2040,
        "spectral_type": "S",
        "albedo": 0.248,
        "source": "Viikinkoski et al. 2017, A&A 607, A117; DAMIT Model 1838",
        "notes": (
            "DAMIT Model 1838. Pole solutions: Kaasalainen (2002): (3,-67); "
            "Viikinkoski (2017): (0,-68); Vernazza (2021): (356,-70). "
            "Retrograde rotator. Largest S-type MBA."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 7,
        "alcdef_datapoints": 2889,
    },
    52: {
        "asteroid_number": 52,
        "name": "Europa",
        "damit_model_id": 135,
        "spin_lambda_deg": 252.0,
        "spin_beta_deg": 38.0,
        "period_hours": 5.629959,
        "period_uncertainty_hours": 0.00005,
        "pole_uncertainty_deg": 5.0,
        "axis_ratios": [1.0, 0.87, 0.66],  # 379 x 330 x 249 km
        "dimensions_km": [379.0, 330.0, 249.0],
        "equivalent_diameter_km": 303.92,
        "n_vertices": 1022,
        "n_faces": 2040,
        "spectral_type": "C",
        "albedo": 0.058,
        "source": "Michalowski et al. 2004; Merline et al. 2013; DAMIT Model 135",
        "notes": (
            "DAMIT Model 135. Prograde rotator. "
            "Pole: Michalowski (2004): (252,+38); Merline (2013): (255,+35). "
            "Second largest C-type in main belt after Hygiea."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 1,
        "alcdef_datapoints": 636,
    },
    1580: {
        "asteroid_number": 1580,
        "name": "Betulia",
        "damit_model_id": 204,
        "spin_lambda_deg": 136.0,
        "spin_beta_deg": 22.0,
        "period_hours": 6.13836,
        "period_uncertainty_hours": 0.00005,
        "pole_uncertainty_deg": 5.0,
        "axis_ratios": [1.0, 0.83, 0.73],  # approx from radar model
        "dimensions_km": [6.59, 5.47, 4.81],  # Magri et al. 2007
        "equivalent_diameter_km": 5.39,
        "n_vertices": 1022,
        "n_faces": 2040,
        "spectral_type": "C",
        "albedo": 0.077,
        "source": "Kaasalainen et al. 2004; Magri et al. 2007, Icarus 186, 152",
        "notes": (
            "DAMIT Model 204. Confirmed by radar (Magri et al. 2007). "
            "Notable concavity in southern hemisphere. Triple-peaked lightcurve. "
            "NEO (Amor group)."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 5,
        "alcdef_datapoints": 209,
    },
    951: {
        "asteroid_number": 951,
        "name": "Gaspra",
        "damit_model_id": None,  # Spacecraft model; may have DAMIT LC model
        "spin_lambda_deg": 20.0,
        "spin_beta_deg": 22.0,
        "period_hours": 7.042,
        "period_uncertainty_hours": 0.001,
        "pole_uncertainty_deg": 5.0,
        "axis_ratios": [1.0, 0.75, 0.60],  # 18.2 x 10.5 x 8.9 km
        "dimensions_km": [18.2, 10.5, 8.9],
        "equivalent_diameter_km": 12.2,
        "n_vertices": 32768,  # Thomas et al. Galileo model
        "n_faces": 65536,
        "spectral_type": "S",
        "albedo": 0.22,
        "source": "Thomas et al. 1994, Icarus 107, 23; Galileo flyby",
        "notes": (
            "Ground truth from Galileo spacecraft (1991 flyby). "
            "Ecliptic pole: lambda~20, beta~+22 (±5 deg). "
            "Period = 7.042 h (0.2934197 d). "
            "Shape model from Galileo images (Thomas et al. 1994)."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 11,
        "alcdef_datapoints": 306,
    },
    243: {
        "asteroid_number": 243,
        "name": "Ida",
        "damit_model_id": None,  # Spacecraft model from Galileo
        "spin_lambda_deg": 73.0,
        "spin_beta_deg": -53.0,
        "period_hours": 4.633632,
        "period_uncertainty_hours": 0.000007,
        "pole_uncertainty_deg": 10.0,
        "axis_ratios": [1.0, 0.45, 0.40],  # 59.8 x 25.4 x 18.6 km
        "dimensions_km": [59.8, 25.4, 18.6],
        "equivalent_diameter_km": 31.4,
        "n_vertices": 32768,
        "n_faces": 65536,
        "spectral_type": "S",
        "albedo": 0.238,
        "source": "Belton et al. 1996, Icarus 120, 1; Galileo flyby",
        "notes": (
            "Ground truth from Galileo spacecraft (1993 flyby). "
            "Two pole solutions: ecliptic (73,-53) and (254,-55). "
            "Retrograde rotator. Has satellite Dactyl. "
            "Highly elongated shape."
        ),
        "in_alcdef": True,
        "alcdef_lightcurves": 5,
        "alcdef_datapoints": 356,
    },
}


def get_ground_truth(asteroid_number):
    """Retrieve ground truth parameters for an asteroid.

    Parameters
    ----------
    asteroid_number : int
        IAU asteroid number.

    Returns
    -------
    dict or None
        Ground truth dictionary if available, else None.
    """
    return DAMIT_GROUND_TRUTH.get(asteroid_number)


def list_ground_truth_asteroids():
    """Return list of asteroid numbers with ground truth data."""
    return sorted(DAMIT_GROUND_TRUTH.keys())


# ============================================================
# Ellipsoid approximation mesh generation
# ============================================================

def generate_ellipsoid_obj(a, b, c, n_lat=32, n_lon=64):
    """Generate a triangulated ellipsoid mesh.

    Parameters
    ----------
    a, b, c : float
        Semi-axis lengths along x, y, z respectively.
    n_lat : int
        Number of latitude divisions.
    n_lon : int
        Number of longitude divisions.

    Returns
    -------
    vertices : (N, 3) ndarray
    faces : (M, 3) ndarray of int (0-indexed)
    """
    vertices = []
    faces = []

    # Generate vertices
    for i in range(n_lat + 1):
        theta = np.pi * i / n_lat  # 0 to pi (pole to pole)
        for j in range(n_lon):
            phi = 2.0 * np.pi * j / n_lon  # 0 to 2*pi
            x = a * np.sin(theta) * np.cos(phi)
            y = b * np.sin(theta) * np.sin(phi)
            z = c * np.cos(theta)
            vertices.append([x, y, z])

    vertices = np.array(vertices, dtype=np.float64)

    # Generate faces (triangle strips)
    for i in range(n_lat):
        for j in range(n_lon):
            # Current row start index
            curr = i * n_lon + j
            next_j = i * n_lon + (j + 1) % n_lon
            # Next row
            below = (i + 1) * n_lon + j
            below_next = (i + 1) * n_lon + (j + 1) % n_lon

            # Two triangles per quad
            faces.append([curr, below, next_j])
            faces.append([next_j, below, below_next])

    faces = np.array(faces, dtype=np.int32)
    return vertices, faces


def generate_ground_truth_ellipsoid(asteroid_number, output_dir=None):
    """Generate an ellipsoid .obj approximation for a ground truth asteroid.

    The ellipsoid uses the published axis ratios from the ground truth data.
    The model is unit-volume normalized (consistent with DAMIT convention).

    Parameters
    ----------
    asteroid_number : int
        IAU asteroid number.
    output_dir : str, optional
        Directory to save the .obj file. Defaults to results/ground_truth/.

    Returns
    -------
    dict with keys: 'vertices', 'faces', 'obj_path', 'axis_ratios'
    """
    gt = get_ground_truth(asteroid_number)
    if gt is None:
        raise ValueError(f"No ground truth data for asteroid {asteroid_number}")

    if output_dir is None:
        output_dir = GROUND_TRUTH_DIR

    # Use axis ratios: a=1, b=b/a, c=c/a
    ratios = gt["axis_ratios"]
    a_ratio, b_ratio, c_ratio = ratios

    # Normalize to unit volume: V = (4/3)*pi*a*b*c = 1
    # => a*b*c = 3/(4*pi)
    vol_factor = (3.0 / (4.0 * np.pi)) ** (1.0 / 3.0)
    product = (a_ratio * b_ratio * c_ratio) ** (1.0 / 3.0)
    scale = vol_factor / product

    a = a_ratio * scale
    b = b_ratio * scale
    c = c_ratio * scale

    vertices, faces = generate_ellipsoid_obj(a, b, c, n_lat=32, n_lon=64)

    os.makedirs(output_dir, exist_ok=True)
    obj_path = os.path.join(output_dir, f"{asteroid_number}_{gt['name']}_ellipsoid.obj")
    save_obj(obj_path, vertices, faces)

    return {
        "vertices": vertices,
        "faces": faces,
        "obj_path": obj_path,
        "axis_ratios": ratios,
        "scale_factors": [a, b, c],
        "n_vertices": len(vertices),
        "n_faces": len(faces),
    }


# ============================================================
# JSON export
# ============================================================

def export_ground_truth_json(output_dir=None):
    """Export all ground truth parameters to a JSON file.

    Parameters
    ----------
    output_dir : str, optional
        Defaults to results/ground_truth/.

    Returns
    -------
    str : path to the saved JSON file
    """
    if output_dir is None:
        output_dir = GROUND_TRUTH_DIR

    os.makedirs(output_dir, exist_ok=True)

    # Build export data (exclude numpy arrays)
    export = {
        "metadata": {
            "description": (
                "Ground truth asteroid shape and spin parameters from DAMIT "
                "(Database of Asteroid Models from Inversion Techniques) and "
                "spacecraft missions. Used for blind validation of lightcurve "
                "inversion pipeline."
            ),
            "damit_url": "https://astro.troja.mff.cuni.cz/projects/damit/",
            "damit_reference": "Durech, J. et al. 2010, A&A 513, A46",
            "generated_at": datetime.utcnow().isoformat() + "Z",
            "num_asteroids": len(DAMIT_GROUND_TRUTH),
        },
        "asteroids": {},
    }

    for num, gt in sorted(DAMIT_GROUND_TRUTH.items()):
        export["asteroids"][str(num)] = {
            "asteroid_number": gt["asteroid_number"],
            "name": gt["name"],
            "damit_model_id": gt["damit_model_id"],
            "spin_parameters": {
                "lambda_deg": gt["spin_lambda_deg"],
                "beta_deg": gt["spin_beta_deg"],
                "period_hours": gt["period_hours"],
                "period_uncertainty_hours": gt["period_uncertainty_hours"],
                "pole_uncertainty_deg": gt["pole_uncertainty_deg"],
            },
            "shape": {
                "axis_ratios": gt["axis_ratios"],
                "dimensions_km": gt.get("dimensions_km"),
                "equivalent_diameter_km": gt.get("equivalent_diameter_km"),
                "n_vertices_damit": gt["n_vertices"],
                "n_faces_damit": gt["n_faces"],
            },
            "physical": {
                "spectral_type": gt.get("spectral_type"),
                "albedo": gt.get("albedo"),
            },
            "data_availability": {
                "in_alcdef": gt["in_alcdef"],
                "alcdef_lightcurves": gt.get("alcdef_lightcurves"),
                "alcdef_datapoints": gt.get("alcdef_datapoints"),
            },
            "source": gt["source"],
            "notes": gt["notes"],
        }

    json_path = os.path.join(output_dir, "ground_truth_info.json")
    with open(json_path, "w") as f:
        json.dump(export, f, indent=2)

    print(f"Ground truth JSON saved to {json_path}")
    return json_path


def export_summary_table():
    """Print a summary table of all ground truth asteroids."""
    header = (
        f"{'Number':>6}  {'Name':<12}  {'DAMIT ID':>8}  "
        f"{'Lambda':>7}  {'Beta':>6}  {'Period(h)':>10}  "
        f"{'a/a':>4} {'b/a':>5} {'c/a':>5}  "
        f"{'ALCDEF LCs':>10}  {'Source'}"
    )
    print(header)
    print("-" * len(header))

    for num in sorted(DAMIT_GROUND_TRUTH.keys()):
        gt = DAMIT_GROUND_TRUTH[num]
        mid = gt["damit_model_id"]
        mid_str = str(mid) if mid is not None else "N/A"
        ratios = gt["axis_ratios"]
        print(
            f"{gt['asteroid_number']:>6}  {gt['name']:<12}  {mid_str:>8}  "
            f"{gt['spin_lambda_deg']:>7.1f}  {gt['spin_beta_deg']:>6.1f}  "
            f"{gt['period_hours']:>10.6f}  "
            f"{ratios[0]:>4.2f} {ratios[1]:>5.2f} {ratios[2]:>5.2f}  "
            f"{gt.get('alcdef_lightcurves', 'N/A'):>10}  "
            f"{gt['source'][:50]}"
        )


# ============================================================
# Main execution
# ============================================================

if __name__ == "__main__":
    print("=" * 80)
    print("DAMIT Ground Truth Retrieval")
    print("=" * 80)
    print()

    # 1. Print summary table
    print("Ground truth asteroids with ALCDEF data:")
    print()
    export_summary_table()
    print()

    # 2. Export JSON
    print("Exporting ground truth parameters to JSON...")
    json_path = export_ground_truth_json()
    print()

    # 3. Generate ellipsoid approximation models
    print("Generating ellipsoid approximation .obj models...")
    print()

    # Generate for the 3 primary validation targets
    primary_targets = [433, 1036, 15]
    additional_targets = [52, 1580, 951, 243]

    for ast_num in primary_targets + additional_targets:
        gt = DAMIT_GROUND_TRUTH[ast_num]
        result = generate_ground_truth_ellipsoid(ast_num)
        print(
            f"  {ast_num:>5} {gt['name']:<12}: "
            f"{result['n_vertices']} vertices, {result['n_faces']} faces "
            f"-> {result['obj_path']}"
        )

    print()
    print("=" * 80)
    print("Summary of primary validation targets:")
    print("=" * 80)
    print()

    for ast_num in primary_targets:
        gt = DAMIT_GROUND_TRUTH[ast_num]
        print(f"  [{ast_num}] {gt['name']}:")
        print(f"    DAMIT Model ID:  {gt['damit_model_id'] or 'N/A (spacecraft model)'}")
        print(f"    Spin pole:       lambda = {gt['spin_lambda_deg']:.2f} deg, "
              f"beta = {gt['spin_beta_deg']:.2f} deg")
        print(f"    Period:          {gt['period_hours']:.8f} hours")
        print(f"    Axis ratios:     a:b:c = 1.00 : {gt['axis_ratios'][1]:.3f} : "
              f"{gt['axis_ratios'][2]:.3f}")
        if gt.get("dimensions_km"):
            dims = gt["dimensions_km"]
            print(f"    Dimensions:      {dims[0]} x {dims[1]} x {dims[2]} km")
        print(f"    Source:          {gt['source']}")
        print(f"    ALCDEF data:     {gt['alcdef_lightcurves']} lightcurves, "
              f"{gt['alcdef_datapoints']} points")
        print()

    print("Ground truth retrieval complete.")
    print(f"Files saved to: {GROUND_TRUTH_DIR}/")
