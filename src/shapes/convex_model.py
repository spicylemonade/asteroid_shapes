"""
Convex Shape Representation and Synthetic Lightcurve Forward Model

Implements:
- Spherical harmonics shape parameterization (degree/order up to 8)
- Triangulated mesh generation from coefficients
- Facet normal and area computation
- Lommel-Seeliger + Lambert scattering law
- Synthetic brightness calculation
"""

import numpy as np
from scipy.special import sph_harm

DEG2RAD = np.pi / 180.0
RAD2DEG = 180.0 / np.pi


def create_icosphere(subdivisions=3):
    """Create an icosphere mesh with given subdivision level.

    Returns vertices (Nx3) and faces (Mx3, indices).
    """
    # Start with icosahedron
    t = (1.0 + np.sqrt(5.0)) / 2.0
    vertices = [
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ]
    vertices = np.array(vertices, dtype=np.float64)
    vertices = vertices / np.linalg.norm(vertices, axis=1, keepdims=True)

    faces = [
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ]
    faces = np.array(faces, dtype=np.int64)

    # Subdivide
    for _ in range(subdivisions):
        edge_midpoints = {}
        new_faces = []
        new_verts = list(vertices)

        def get_midpoint(i, j):
            key = (min(i, j), max(i, j))
            if key in edge_midpoints:
                return edge_midpoints[key]
            mid = (new_verts[i] + new_verts[j]) / 2.0
            mid = mid / np.linalg.norm(mid)
            idx = len(new_verts)
            new_verts.append(mid)
            edge_midpoints[key] = idx
            return idx

        for tri in faces:
            a, b, c = tri
            ab = get_midpoint(a, b)
            bc = get_midpoint(b, c)
            ca = get_midpoint(c, a)
            new_faces.extend([
                [a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]
            ])

        vertices = np.array(new_verts)
        faces = np.array(new_faces)

    return vertices, faces


def real_spherical_harmonics(l, m, theta, phi):
    """Compute real spherical harmonics Y_l^m(theta, phi).

    theta: colatitude [0, pi]
    phi: longitude [0, 2pi)
    """
    if m > 0:
        return np.sqrt(2) * np.real(sph_harm(m, l, phi, theta))
    elif m == 0:
        return np.real(sph_harm(0, l, phi, theta))
    else:
        return np.sqrt(2) * np.imag(sph_harm(-m, l, phi, theta))


def spherical_harmonics_radius(coeffs, theta, phi, lmax=8):
    """Compute radius from spherical harmonics coefficients.

    coeffs: flat array of coefficients, ordered as (l=0,m=0), (l=1,m=-1), (l=1,m=0), (l=1,m=1), ...
    """
    r = np.zeros_like(theta, dtype=np.float64)
    idx = 0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            if idx >= len(coeffs):
                break
            r += coeffs[idx] * real_spherical_harmonics(l, m, theta, phi)
            idx += 1
    return r


def coeffs_for_ellipsoid(a_axis, b_axis, c_axis, lmax=8):
    """Compute approximate spherical harmonics coefficients for a triaxial ellipsoid.

    Uses the fact that for an ellipsoid with semi-axes (a, b, c) aligned with
    body axes, the radius function is approximately:
    r(theta, phi) = 1/sqrt((sin(theta)cos(phi)/a)^2 + (sin(theta)sin(phi)/b)^2 + (cos(theta)/c)^2)
    """
    # We'll use a least-squares fit on the icosphere points
    verts, _ = create_icosphere(subdivisions=3)
    theta = np.arccos(np.clip(verts[:, 2], -1, 1))
    phi = np.arctan2(verts[:, 1], verts[:, 0]) % (2 * np.pi)

    # True radii for ellipsoid
    r_true = 1.0 / np.sqrt(
        (np.sin(theta) * np.cos(phi) / a_axis)**2 +
        (np.sin(theta) * np.sin(phi) / b_axis)**2 +
        (np.cos(theta) / c_axis)**2
    )

    # Build design matrix
    n_coeffs = (lmax + 1)**2
    A = np.zeros((len(theta), n_coeffs))
    idx = 0
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            A[:, idx] = real_spherical_harmonics(l, m, theta, phi)
            idx += 1

    # Least-squares fit
    coeffs, _, _, _ = np.linalg.lstsq(A, r_true, rcond=None)
    return coeffs


class ConvexShapeModel:
    """Convex asteroid shape model using spherical harmonics."""

    def __init__(self, lmax=8, subdivisions=3):
        self.lmax = lmax
        self.n_coeffs = (lmax + 1)**2
        self.subdivisions = subdivisions

        # Create base icosphere
        self._base_verts, self.faces = create_icosphere(subdivisions)
        self.n_vertices = len(self._base_verts)
        self.n_facets = len(self.faces)

        # Spherical coordinates of base vertices (unit sphere)
        self._theta = np.arccos(np.clip(self._base_verts[:, 2], -1, 1))
        self._phi = np.arctan2(self._base_verts[:, 1], self._base_verts[:, 0]) % (2 * np.pi)

        # Precompute spherical harmonics basis
        self._basis = np.zeros((self.n_vertices, self.n_coeffs))
        idx = 0
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                self._basis[:, idx] = real_spherical_harmonics(l, m, self._theta, self._phi)
                idx += 1

        # Default coefficients (sphere)
        self.coeffs = np.zeros(self.n_coeffs)
        self.coeffs[0] = np.sqrt(4 * np.pi)  # Y_0^0 normalization gives r=1

        self._update_mesh()

    def set_coefficients(self, coeffs):
        """Set spherical harmonics coefficients and update mesh."""
        self.coeffs = np.array(coeffs[:self.n_coeffs], dtype=np.float64)
        self._update_mesh()

    def _update_mesh(self):
        """Recompute mesh vertices, normals, and areas from current coefficients."""
        radii = self._basis @ self.coeffs
        # Ensure positive radii
        radii = np.maximum(radii, 0.01)
        self.vertices = self._base_verts * radii[:, np.newaxis]

        # Compute facet normals and areas
        v0 = self.vertices[self.faces[:, 0]]
        v1 = self.vertices[self.faces[:, 1]]
        v2 = self.vertices[self.faces[:, 2]]
        cross = np.cross(v1 - v0, v2 - v0)
        self.facet_areas = 0.5 * np.linalg.norm(cross, axis=1)
        # Avoid division by zero
        norms = np.linalg.norm(cross, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-30)
        self.facet_normals = cross / norms

        # Facet centroids
        self.facet_centroids = (v0 + v1 + v2) / 3.0

    def get_mesh(self):
        """Return vertices and faces for mesh export."""
        return self.vertices.copy(), self.faces.copy()

    def save_obj(self, filepath):
        """Save mesh as OBJ file."""
        with open(filepath, 'w') as f:
            f.write("# Asteroid convex shape model\n")
            for v in self.vertices:
                f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
            for face in self.faces:
                f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")


def lommel_seeliger_lambert(mu, mu0, alpha, c_lambert=0.1, phase_coeff=0.01):
    """Combined Lommel-Seeliger + Lambert scattering law.

    Parameters
    ----------
    mu : array, cos(emission angle)
    mu0 : array, cos(incidence angle)
    alpha : float, phase angle in radians
    c_lambert : float, Lambert scattering weight
    phase_coeff : float, linear phase coefficient

    Returns
    -------
    Scattered intensity per facet
    """
    denom = mu + mu0
    denom = np.maximum(denom, 1e-10)

    # Phase function
    f_phase = 1.0 - phase_coeff * np.abs(alpha)

    # Lommel-Seeliger term
    ls = mu0 / denom * f_phase

    # Lambert term
    lambert = c_lambert * mu0

    return ls + lambert


def compute_brightness(shape_model, sun_dir, obs_dir, phase_angle=None):
    """Compute total reflected brightness for given geometry.

    Parameters
    ----------
    shape_model : ConvexShapeModel
    sun_dir : array (3,), unit vector from asteroid toward Sun (body frame)
    obs_dir : array (3,), unit vector from asteroid toward observer (body frame)
    phase_angle : float, phase angle in radians (computed from directions if None)

    Returns
    -------
    Total brightness (proportional to flux)
    """
    normals = shape_model.facet_normals
    areas = shape_model.facet_areas

    # Cosines of incidence and emission angles
    mu0 = np.dot(normals, sun_dir)    # cos(incidence)
    mu = np.dot(normals, obs_dir)     # cos(emission)

    # Visibility mask: facet must face both Sun and observer
    visible = (mu > 0) & (mu0 > 0)

    if phase_angle is None:
        cos_alpha = np.clip(np.dot(sun_dir, obs_dir), -1, 1)
        phase_angle = np.arccos(cos_alpha)

    # Scattering
    scatter = lommel_seeliger_lambert(mu, mu0, phase_angle)

    # Total brightness
    brightness = np.sum(areas[visible] * scatter[visible])

    return brightness


def spin_rotation_matrix(pole_lambda, pole_beta, period, phi0, t, t0):
    """Compute rotation matrix from ecliptic to body-fixed frame.

    Parameters
    ----------
    pole_lambda : float, pole ecliptic longitude (radians)
    pole_beta : float, pole ecliptic latitude (radians)
    period : float, sidereal rotation period (days)
    phi0 : float, rotational phase at t0 (radians)
    t : float, current epoch (JD)
    t0 : float, reference epoch (JD)

    Returns
    -------
    R : 3x3 rotation matrix (ecliptic -> body-fixed)
    """
    # Spin phase
    phase = phi0 + 2.0 * np.pi * (t - t0) / period
    phase = phase % (2.0 * np.pi)

    # Spin axis in ecliptic coordinates
    sz = np.array([np.cos(pole_beta) * np.cos(pole_lambda),
                   np.cos(pole_beta) * np.sin(pole_lambda),
                   np.sin(pole_beta)])

    # Build body frame: z_body = spin axis
    # x_body: arbitrary direction perpendicular to spin axis
    if abs(sz[2]) < 0.9:
        sx = np.cross(np.array([0, 0, 1]), sz)
    else:
        sx = np.cross(np.array([1, 0, 0]), sz)
    sx = sx / np.linalg.norm(sx)
    sy = np.cross(sz, sx)

    # Pole orientation matrix (ecliptic -> body at phase=0)
    R_pole = np.array([sx, sy, sz])  # rows are body axes in ecliptic coords

    # Rotation about spin axis
    c, s = np.cos(phase), np.sin(phase)
    R_spin = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])

    # Combined: ecliptic -> body-fixed
    R = R_spin @ R_pole
    return R


def synthetic_lightcurve(shape_model, times, pole_lambda, pole_beta, period, phi0, t0,
                         sun_dirs_ecl, obs_dirs_ecl, phase_angles):
    """Compute synthetic lightcurve for a sequence of observation epochs.

    Parameters
    ----------
    shape_model : ConvexShapeModel
    times : array of JDs
    pole_lambda, pole_beta : pole direction in radians
    period : rotation period in days
    phi0 : initial phase in radians
    t0 : reference epoch JD
    sun_dirs_ecl : array (N, 3), unit vectors toward Sun in ecliptic
    obs_dirs_ecl : array (N, 3), unit vectors toward observer in ecliptic
    phase_angles : array (N,), phase angles in radians

    Returns
    -------
    magnitudes : array (N,), relative magnitudes (smaller = brighter)
    """
    n_obs = len(times)
    flux = np.zeros(n_obs)

    for i in range(n_obs):
        R = spin_rotation_matrix(pole_lambda, pole_beta, period, phi0,
                                 times[i], t0)
        # Transform directions to body frame
        sun_body = R @ sun_dirs_ecl[i]
        obs_body = R @ obs_dirs_ecl[i]

        flux[i] = compute_brightness(shape_model, sun_body, obs_body,
                                     phase_angles[i])

    # Convert to relative magnitudes
    flux = np.maximum(flux, 1e-30)
    mag = -2.5 * np.log10(flux)

    return mag
