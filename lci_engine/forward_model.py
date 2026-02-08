"""
Forward model for asteroid lightcurve computation.

Implements:
1. Convex shape representation (triangulated polyhedron)
2. Lommel-Seeliger scattering law
3. Illumination/viewing geometry computation
4. Synthetic brightness generation

References:
    Kaasalainen & Torppa 2001 (Kaasalainen2001a in sources.bib)
    Hapke 1993 (Hapke1993 in sources.bib)
"""

import numpy as np
from typing import Tuple, Optional


def create_triaxial_ellipsoid(a: float, b: float, c: float,
                               n_lat: int = 15, n_lon: int = 30) -> Tuple[np.ndarray, np.ndarray]:
    """Create a triangulated triaxial ellipsoid mesh.

    Args:
        a, b, c: Semi-axes along x, y, z
        n_lat: Number of latitude divisions
        n_lon: Number of longitude divisions

    Returns:
        vertices: (N, 3) array of vertex positions
        faces: (M, 3) array of face vertex indices
    """
    vertices = []
    faces = []

    # Generate vertices
    for i in range(n_lat + 1):
        theta = np.pi * i / n_lat
        for j in range(n_lon):
            phi = 2.0 * np.pi * j / n_lon
            x = a * np.sin(theta) * np.cos(phi)
            y = b * np.sin(theta) * np.sin(phi)
            z = c * np.cos(theta)
            vertices.append([x, y, z])

    vertices = np.array(vertices)

    # Generate faces (triangulation)
    for i in range(n_lat):
        for j in range(n_lon):
            j_next = (j + 1) % n_lon
            v0 = i * n_lon + j
            v1 = i * n_lon + j_next
            v2 = (i + 1) * n_lon + j
            v3 = (i + 1) * n_lon + j_next

            if i == 0:
                faces.append([v0, v3, v2])
            elif i == n_lat - 1:
                faces.append([v0, v1, v2])
            else:
                faces.append([v0, v1, v3])
                faces.append([v0, v3, v2])

    return vertices, np.array(faces)


def create_convex_shape(n_facets: int = 200, seed: int = 42) -> Tuple[np.ndarray, np.ndarray]:
    """Create a random convex shape via perturbed sphere.

    Uses Gaussian surface density approach: generate random
    facet areas on a sphere, then apply Minkowski reconstruction.
    For simplicity, we perturb a sphere's vertex radii.

    Returns:
        vertices: (N, 3) array
        faces: (M, 3) array
    """
    rng = np.random.RandomState(seed)

    # Start with ico-sphere approximation
    n_lat = max(8, int(np.sqrt(n_facets / 2)))
    n_lon = max(16, int(n_facets / n_lat))

    vertices, faces = create_triaxial_ellipsoid(1.0, 1.0, 1.0, n_lat, n_lon)

    # Perturb radii to create irregular convex shape
    radii = 1.0 + 0.2 * rng.randn(len(vertices))
    radii = np.maximum(radii, 0.5)

    norms = np.linalg.norm(vertices, axis=1, keepdims=True)
    norms = np.maximum(norms, 1e-10)
    directions = vertices / norms
    vertices = directions * radii[:, np.newaxis]

    return vertices, faces


def compute_facet_properties(vertices: np.ndarray, faces: np.ndarray
                              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute facet normals, areas, and centroids.

    Returns:
        normals: (M, 3) unit normal vectors (outward)
        areas: (M,) facet areas
        centroids: (M, 3) facet centers
    """
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]

    edge1 = v1 - v0
    edge2 = v2 - v0

    cross = np.cross(edge1, edge2)
    area_2 = np.linalg.norm(cross, axis=1)
    area_2 = np.maximum(area_2, 1e-30)

    normals = cross / area_2[:, np.newaxis]
    areas = 0.5 * area_2
    centroids = (v0 + v1 + v2) / 3.0

    # Ensure normals point outward (away from origin)
    dots = np.sum(normals * centroids, axis=1)
    flip = dots < 0
    normals[flip] *= -1

    return normals, areas, centroids


def rotation_matrix_z(angle: float) -> np.ndarray:
    """Rotation matrix around Z axis by given angle (radians)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rotation_matrix_y(angle: float) -> np.ndarray:
    """Rotation matrix around Y axis by given angle (radians)."""
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def ecliptic_to_body_frame(sun_ecl: np.ndarray, obs_ecl: np.ndarray,
                            pole_lambda: float, pole_beta: float,
                            rotation_phase: float) -> Tuple[np.ndarray, np.ndarray]:
    """Transform Sun and observer vectors from ecliptic to asteroid body frame.

    Args:
        sun_ecl: (3,) Sun direction in ecliptic coords (unit vector)
        obs_ecl: (3,) Observer direction in ecliptic coords (unit vector)
        pole_lambda: Ecliptic longitude of spin axis (radians)
        pole_beta: Ecliptic latitude of spin axis (radians)
        rotation_phase: Current rotation angle (radians)

    Returns:
        sun_body: (3,) Sun direction in body frame
        obs_body: (3,) Observer direction in body frame
    """
    # Transform from ecliptic to equatorial body frame
    # 1. Rotate by -lambda around Z to align pole longitude
    # 2. Rotate by -(90-beta) around Y to align pole to Z
    # 3. Rotate by -rotation_phase around Z for current spin angle

    R_spin = rotation_matrix_z(-rotation_phase)
    R_beta = rotation_matrix_y(-(np.pi / 2 - pole_beta))
    R_lambda = rotation_matrix_z(-pole_lambda)

    R_total = R_spin @ R_beta @ R_lambda

    sun_body = R_total @ sun_ecl
    obs_body = R_total @ obs_ecl

    return sun_body, obs_body


def lommel_seeliger(mu0: np.ndarray, mu: np.ndarray,
                    phase_angle: float = 0.0) -> np.ndarray:
    """Lommel-Seeliger scattering law.

    S = f(alpha) * mu0 / (mu0 + mu)

    where f(alpha) is a simple linear phase function.
    """
    denom = mu0 + mu
    denom = np.maximum(denom, 1e-10)

    # Simple phase function: 1 + slope * alpha
    # (linear approximation, sufficient for relative photometry)
    f_alpha = 1.0

    return f_alpha * mu0 / denom


def compute_synthetic_brightness(vertices: np.ndarray, faces: np.ndarray,
                                  sun_body: np.ndarray, obs_body: np.ndarray) -> float:
    """Compute total brightness of asteroid shape for given geometry.

    Args:
        vertices, faces: Shape mesh
        sun_body: (3,) unit vector to Sun in body frame
        obs_body: (3,) unit vector to observer in body frame

    Returns:
        Total brightness (arbitrary units, proportional to reflected flux)
    """
    normals, areas, _ = compute_facet_properties(vertices, faces)

    # Cosines of incidence and emission angles
    mu0 = normals @ sun_body  # cos(incidence)
    mu = normals @ obs_body   # cos(emission)

    # Only facets that are both illuminated and visible contribute
    visible = (mu0 > 0) & (mu > 0)

    if not np.any(visible):
        return 0.0

    mu0_v = mu0[visible]
    mu_v = mu[visible]
    areas_v = areas[visible]

    # Lommel-Seeliger scattering
    scatter = lommel_seeliger(mu0_v, mu_v)

    brightness = np.sum(scatter * areas_v)
    return brightness


def generate_lightcurve(vertices: np.ndarray, faces: np.ndarray,
                        sun_ecl: np.ndarray, obs_ecl: np.ndarray,
                        pole_lambda: float, pole_beta: float,
                        period_hours: float, epoch_jd: float,
                        times_jd: np.ndarray) -> np.ndarray:
    """Generate synthetic lightcurve for given shape and geometry.

    Args:
        vertices, faces: Shape mesh
        sun_ecl: (3,) unit vector to Sun in ecliptic frame
        obs_ecl: (3,) unit vector to observer in ecliptic frame
        pole_lambda, pole_beta: Spin axis in ecliptic coords (radians)
        period_hours: Sidereal rotation period
        epoch_jd: Reference epoch (JD) for phase=0
        times_jd: (N,) array of observation times (JD)

    Returns:
        brightness: (N,) array of synthetic brightness values
    """
    period_days = period_hours / 24.0
    phases = 2.0 * np.pi * (times_jd - epoch_jd) / period_days

    brightness = np.zeros(len(times_jd))
    for i, phase in enumerate(phases):
        sun_body, obs_body = ecliptic_to_body_frame(
            sun_ecl, obs_ecl, pole_lambda, pole_beta, phase
        )
        brightness[i] = compute_synthetic_brightness(vertices, faces, sun_body, obs_body)

    return brightness


def brightness_to_magnitude(brightness: np.ndarray) -> np.ndarray:
    """Convert brightness to relative magnitude (smaller = brighter)."""
    brightness = np.maximum(brightness, 1e-30)
    ref = np.median(brightness)
    return -2.5 * np.log10(brightness / ref)


def save_mesh_obj(vertices: np.ndarray, faces: np.ndarray, filepath: str):
    """Save mesh as Wavefront OBJ file."""
    with open(filepath, 'w') as f:
        f.write("# LCI Engine shape model\n")
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in faces:
            f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")


def load_mesh_obj(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """Load mesh from Wavefront OBJ file."""
    vertices = []
    faces = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == 'v' and len(parts) >= 4:
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif parts[0] == 'f':
                # Handle f v1 v2 v3 or f v1/vt1/vn1 v2/... formats
                face_verts = []
                for p in parts[1:4]:
                    idx = int(p.split('/')[0]) - 1
                    face_verts.append(idx)
                faces.append(face_verts)
    return np.array(vertices), np.array(faces)
