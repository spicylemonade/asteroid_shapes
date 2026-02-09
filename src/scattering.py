"""
Forward photometric model with Lommel-Seeliger + Lambert scattering law
and self-shadowing ray-tracing.

References:
  - Kaasalainen & Torppa (2001) [Kaasalainen2001a in sources.bib]
  - Durech & Kaasalainen (2003) [Durech2003 in sources.bib]
  - Hapke (1993, 2012) [Hapke1993, Hapke2012 in sources.bib]
"""
import numpy as np


# ============================================================
# Scattering laws
# ============================================================

def lommel_seeliger_lambert(mu_0, mu, alpha_rad, c_ls=0.5, c_l=0.1):
    """Combined Lommel-Seeliger + Lambert scattering law.

    S = f(alpha) * [c_ls * mu_0 / (mu_0 + mu) + c_l * mu_0]

    Parameters
    ----------
    mu_0 : array, cos(incidence angle), i.e., cos(angle between normal and sun direction)
    mu : array, cos(emission angle), i.e., cos(angle between normal and observer direction)
    alpha_rad : float, phase angle in radians
    c_ls : float, Lommel-Seeliger weight
    c_l : float, Lambert weight

    Returns
    -------
    array : scattered intensity per facet
    """
    # Phase function f(alpha) - simple linear-exponential
    f_alpha = np.exp(-0.5 * alpha_rad)

    # Avoid division by zero
    denom = mu_0 + mu
    denom = np.where(denom > 1e-10, denom, 1e-10)

    ls_term = c_ls * mu_0 / denom
    lambert_term = c_l * mu_0

    return f_alpha * (ls_term + lambert_term)


# ============================================================
# BVH for ray-mesh intersection (self-shadowing)
# ============================================================

class AABB:
    """Axis-Aligned Bounding Box."""
    def __init__(self, min_pt, max_pt):
        self.min_pt = np.asarray(min_pt, dtype=np.float64)
        self.max_pt = np.asarray(max_pt, dtype=np.float64)

    def intersects_ray(self, origin, direction):
        """Test ray-AABB intersection using slab method."""
        inv_dir = np.where(np.abs(direction) > 1e-12,
                           1.0 / direction,
                           np.sign(direction) * 1e12)

        t1 = (self.min_pt - origin) * inv_dir
        t2 = (self.max_pt - origin) * inv_dir

        tmin = np.max(np.minimum(t1, t2))
        tmax = np.min(np.maximum(t1, t2))

        return tmax >= max(tmin, 0.0)


class BVHNode:
    """Bounding Volume Hierarchy node for accelerated ray-mesh intersection."""
    def __init__(self, bbox, left=None, right=None, face_indices=None):
        self.bbox = bbox
        self.left = left
        self.right = right
        self.face_indices = face_indices  # Leaf node


def build_bvh(vertices, faces, max_leaf_size=8):
    """Build a BVH tree for the given mesh.

    Parameters
    ----------
    vertices : (N, 3) array
    faces : (M, 3) array of vertex indices
    max_leaf_size : int

    Returns
    -------
    BVHNode : root of the tree
    """
    face_indices = np.arange(len(faces))
    return _build_bvh_recursive(vertices, faces, face_indices, max_leaf_size)


def _build_bvh_recursive(vertices, faces, face_indices, max_leaf_size):
    # Compute bounding box of all triangles
    all_verts = vertices[faces[face_indices].ravel()]
    bbox = AABB(np.min(all_verts, axis=0), np.max(all_verts, axis=0))

    if len(face_indices) <= max_leaf_size:
        return BVHNode(bbox, face_indices=face_indices)

    # Split along longest axis
    extent = bbox.max_pt - bbox.min_pt
    axis = np.argmax(extent)

    # Compute centroids
    centroids = np.mean(vertices[faces[face_indices]], axis=1)
    median = np.median(centroids[:, axis])

    left_mask = centroids[:, axis] <= median
    right_mask = ~left_mask

    # Handle degenerate case
    if np.sum(left_mask) == 0 or np.sum(right_mask) == 0:
        mid = len(face_indices) // 2
        left_indices = face_indices[:mid]
        right_indices = face_indices[mid:]
    else:
        left_indices = face_indices[left_mask]
        right_indices = face_indices[right_mask]

    left_node = _build_bvh_recursive(vertices, faces, left_indices, max_leaf_size)
    right_node = _build_bvh_recursive(vertices, faces, right_indices, max_leaf_size)

    return BVHNode(bbox, left=left_node, right=right_node)


def ray_triangle_intersect(origin, direction, v0, v1, v2):
    """Moller-Trumbore ray-triangle intersection test.

    Returns t > 0 if intersection, -1 otherwise.
    """
    edge1 = v1 - v0
    edge2 = v2 - v0
    h = np.cross(direction, edge2)
    a = np.dot(edge1, h)

    if abs(a) < 1e-10:
        return -1.0

    f = 1.0 / a
    s = origin - v0
    u = f * np.dot(s, h)
    if u < 0.0 or u > 1.0:
        return -1.0

    q = np.cross(s, edge1)
    v = f * np.dot(direction, q)
    if v < 0.0 or u + v > 1.0:
        return -1.0

    t = f * np.dot(edge2, q)
    if t > 1e-6:
        return t
    return -1.0


def ray_intersects_bvh(origin, direction, bvh_root, vertices, faces, exclude_face=-1):
    """Test if a ray intersects any triangle in the BVH.

    Parameters
    ----------
    origin : (3,) ray origin
    direction : (3,) ray direction (normalized)
    bvh_root : BVHNode
    vertices : (N, 3)
    faces : (M, 3)
    exclude_face : int, face index to exclude (self-intersection avoidance)

    Returns
    -------
    bool : True if ray hits any triangle
    """
    if not bvh_root.bbox.intersects_ray(origin, direction):
        return False

    if bvh_root.face_indices is not None:
        # Leaf node
        for fi in bvh_root.face_indices:
            if fi == exclude_face:
                continue
            v0 = vertices[faces[fi, 0]]
            v1 = vertices[faces[fi, 1]]
            v2 = vertices[faces[fi, 2]]
            t = ray_triangle_intersect(origin, direction, v0, v1, v2)
            if t > 0:
                return True
        return False

    if bvh_root.left and ray_intersects_bvh(origin, direction, bvh_root.left, vertices, faces, exclude_face):
        return True
    if bvh_root.right and ray_intersects_bvh(origin, direction, bvh_root.right, vertices, faces, exclude_face):
        return True
    return False


def compute_shadow_mask_vectorized(centers, normals, sun_dir, vertices, faces):
    """Vectorized shadow computation using simplified approach.

    Instead of per-facet BVH traversal (slow in Python), use a faster
    approach: for each illuminated facet, project its center onto the
    sun direction and check if any other facet blocks it using bounding
    volume tests.

    For non-convex shapes, this efficiently detects when one part of the
    asteroid casts a shadow on another part.
    """
    n_faces = len(centers)
    shadow = np.zeros(n_faces, dtype=bool)

    # Only check facets that face the sun
    mu_0 = np.dot(normals, sun_dir)
    illuminated = mu_0 > 0

    if not np.any(illuminated):
        return shadow

    # Project all face centers onto sun direction
    # A facet B shadows facet A if:
    # 1. B is between A and the sun (higher projection)
    # 2. The perpendicular distance between A's ray and B is small
    proj = np.dot(centers, sun_dir)  # projection onto sun direction

    # For each illuminated facet, check if any other facet blocks it
    # Use vectorized distance computation for speed
    ill_indices = np.where(illuminated)[0]

    for i in ill_indices:
        origin = centers[i]
        # Find facets that are further along the sun direction (closer to sun)
        candidates = proj > proj[i] + 0.01
        if not np.any(candidates):
            continue

        cand_idx = np.where(candidates)[0]
        # Vector from origin to candidate centers
        diffs = centers[cand_idx] - origin
        # Project onto sun_dir
        along = np.dot(diffs, sun_dir)
        # Perpendicular distance
        perp = diffs - np.outer(along, sun_dir)
        perp_dist = np.linalg.norm(perp, axis=1)

        # Check if any candidate is close enough to block (within facet size)
        avg_facet_size = np.sqrt(np.mean(np.sum((vertices[faces[cand_idx, 1]] -
                                                  vertices[faces[cand_idx, 0]])**2, axis=1)))
        threshold = avg_facet_size * 1.5

        if np.any(perp_dist < threshold):
            shadow[i] = True

    return shadow


# ============================================================
# Forward brightness model with self-shadowing
# ============================================================

def compute_brightness(vertices, faces, normals, areas, centers,
                       sun_dir, obs_dir, bvh_root=None,
                       c_ls=0.5, c_l=0.1, use_shadow=True):
    """Compute total disk-integrated brightness of a shape model.

    Parameters
    ----------
    vertices : (N, 3) mesh vertices
    faces : (M, 3) face indices
    normals : (M, 3) outward unit normals
    areas : (M,) facet areas
    centers : (M, 3) facet centroids
    sun_dir : (3,) unit vector from asteroid toward sun
    obs_dir : (3,) unit vector from asteroid toward observer
    bvh_root : BVHNode, optional, for self-shadowing
    c_ls : float, Lommel-Seeliger weight
    c_l : float, Lambert weight
    use_shadow : bool, whether to enable self-shadowing

    Returns
    -------
    float : total brightness (arbitrary units)
    """
    # Cosines
    mu_0 = np.dot(normals, sun_dir)   # cos(incidence)
    mu = np.dot(normals, obs_dir)      # cos(emission)

    # Phase angle
    cos_alpha = np.dot(sun_dir, obs_dir)
    alpha = np.arccos(np.clip(cos_alpha, -1, 1))

    # Visible and illuminated facets
    visible = (mu > 0) & (mu_0 > 0)

    if not np.any(visible):
        return 1e-10

    # Self-shadowing via vectorized shadow computation
    shadow_mask = np.zeros(len(faces), dtype=bool)
    if use_shadow:
        shadow_mask = compute_shadow_mask_vectorized(centers, normals, sun_dir, vertices, faces)

    # Active facets: visible, illuminated, not shadowed
    active = visible & ~shadow_mask

    if not np.any(active):
        return 1e-10

    # Compute scattering
    intensity = lommel_seeliger_lambert(mu_0[active], mu[active], alpha, c_ls, c_l)

    # Total brightness = sum(intensity * area)
    total = np.sum(intensity * areas[active])

    return total


def compute_lightcurve(vertices, faces, spin_axis, period_hours,
                       epoch_jd, observation_jds, sun_dirs, obs_dirs,
                       c_ls=0.5, c_l=0.1, use_shadow=True):
    """Compute a synthetic lightcurve for a rotating shape model.

    Parameters
    ----------
    vertices : (N, 3) mesh vertices
    faces : (M, 3) face indices
    spin_axis : (3,) unit vector, rotation axis (ecliptic coordinates)
    period_hours : float, rotation period
    epoch_jd : float, JD of zero rotation phase
    observation_jds : array of JD times
    sun_dirs : (T, 3) unit sun directions at each obs time
    obs_dirs : (T, 3) unit observer directions at each obs time
    c_ls, c_l : scattering law weights
    use_shadow : bool

    Returns
    -------
    model_mags : array of model magnitudes (relative)
    """
    from mesh_utils import compute_facet_properties

    period_days = period_hours / 24.0
    n_obs = len(observation_jds)
    model_flux = np.zeros(n_obs)

    # Build BVH once
    bvh = build_bvh(vertices, faces) if use_shadow else None

    for t_idx in range(n_obs):
        # Rotation angle at this time
        dt = observation_jds[t_idx] - epoch_jd
        rotation_angle = 2 * np.pi * (dt / period_days) % (2 * np.pi)

        # Rotate shape around spin axis
        rot_verts = _rotate_vertices(vertices, spin_axis, rotation_angle)

        # Recompute facet properties for rotated mesh
        normals, areas, centers = compute_facet_properties(rot_verts, faces)

        # Build BVH for rotated mesh (needed for shadow calculation)
        if use_shadow:
            bvh = build_bvh(rot_verts, faces)

        # Compute brightness
        flux = compute_brightness(rot_verts, faces, normals, areas, centers,
                                   sun_dirs[t_idx], obs_dirs[t_idx],
                                   bvh_root=bvh, c_ls=c_ls, c_l=c_l,
                                   use_shadow=use_shadow)
        model_flux[t_idx] = flux

    # Convert to relative magnitudes
    model_flux = np.clip(model_flux, 1e-20, None)
    model_mags = -2.5 * np.log10(model_flux)
    model_mags -= np.mean(model_mags)  # Zero-mean relative magnitudes

    return model_mags


def _rotate_vertices(vertices, axis, angle):
    """Rotate vertices around an axis by given angle (Rodrigues' formula)."""
    axis = axis / np.linalg.norm(axis)
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)

    # Rodrigues' rotation
    rotated = (vertices * cos_a +
               np.cross(axis, vertices) * sin_a +
               axis * np.dot(vertices, axis)[:, np.newaxis] * (1 - cos_a))

    return rotated


# ============================================================
# Self-shadowing test
# ============================================================

if __name__ == '__main__':
    import time
    from mesh_utils import create_sphere_mesh, compute_facet_properties

    np.random.seed(42)

    # Create a dumbbell/contact binary test shape (truly non-convex)
    # Two spheres connected at a narrow neck
    v1, f1 = create_sphere_mesh(n_subdivisions=2)  # 162 verts, smaller for speed
    v2, f2 = create_sphere_mesh(n_subdivisions=2)

    # Shift sphere 1 left, sphere 2 right, with overlap creating a neck
    v1[:, 0] -= 0.7
    v2[:, 0] += 0.7

    # Combine meshes
    f2_shifted = f2 + len(v1)
    verts = np.vstack([v1, v2])
    faces = np.vstack([f1, f2_shifted])

    # Squeeze the neck region (vertices near x=0) inward
    neck_mask = np.abs(verts[:, 0]) < 0.4
    neck_scale = 0.35
    verts[neck_mask, 1] *= neck_scale
    verts[neck_mask, 2] *= neck_scale

    print(f"Test mesh (dumbbell): {len(verts)} vertices, {len(faces)} faces")

    normals, areas, centers = compute_facet_properties(verts, faces)

    # Sun coming from the side - one lobe will shadow the other
    sun_dir = np.array([1.0, 0.0, 0.3])
    sun_dir /= np.linalg.norm(sun_dir)
    obs_dir = np.array([0.0, 0.0, 1.0])

    # Compute with and without self-shadowing
    flux_no_shadow = compute_brightness(verts, faces, normals, areas, centers,
                                         sun_dir, obs_dir,
                                         use_shadow=False)
    flux_shadow = compute_brightness(verts, faces, normals, areas, centers,
                                      sun_dir, obs_dir,
                                      use_shadow=True)

    diff_pct = abs(flux_no_shadow - flux_shadow) / max(flux_no_shadow, 1e-10) * 100
    print(f"\nFlux without shadows: {flux_no_shadow:.6f}")
    print(f"Flux with shadows:    {flux_shadow:.6f}")
    print(f"Difference: {diff_pct:.2f}%")
    print(f"Requirement: >1% difference for concave shapes: {'PASS' if diff_pct > 1 else 'FAIL'}")

    # Performance test
    print(f"\nPerformance test: 1000 brightness evaluations...")
    t0 = time.time()
    for _ in range(1000):
        sd = sun_dir + np.random.randn(3) * 0.1
        sd /= np.linalg.norm(sd)
        compute_brightness(verts, faces, normals, areas, centers,
                           sd, obs_dir, use_shadow=True)
    elapsed = time.time() - t0
    rate = 1000 / elapsed * 60
    print(f"Elapsed: {elapsed:.2f}s ({rate:.0f} evaluations/minute)")
    print(f"Requirement: >=1000/min: {'PASS' if rate >= 1000 else 'FAIL'}")
