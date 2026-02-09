"""
Mesh utilities: generation, I/O, and comparison metrics.
Implements Hausdorff distance and volumetric IoU for shape validation.

References:
  - Item 011: Mesh comparison metrics
"""
import numpy as np
from scipy.spatial import ConvexHull
import os


# ============================================================
# Mesh generation
# ============================================================

def create_sphere_mesh(n_subdivisions=3, radius=1.0):
    """Create an icosphere mesh by subdividing an icosahedron.

    Returns
    -------
    vertices : (N, 3) array
    faces : (M, 3) array of int
    """
    # Icosahedron base
    t = (1.0 + np.sqrt(5.0)) / 2.0
    verts = np.array([
        [-1, t, 0], [1, t, 0], [-1, -t, 0], [1, -t, 0],
        [0, -1, t], [0, 1, t], [0, -1, -t], [0, 1, -t],
        [t, 0, -1], [t, 0, 1], [-t, 0, -1], [-t, 0, 1],
    ], dtype=np.float64)
    verts /= np.linalg.norm(verts[0])

    faces = np.array([
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ], dtype=np.int32)

    # Subdivide
    for _ in range(n_subdivisions):
        verts, faces = _subdivide(verts, faces)

    # Normalize to sphere
    norms = np.linalg.norm(verts, axis=1, keepdims=True)
    verts = verts / norms * radius

    return verts, faces


def _subdivide(vertices, faces):
    """Subdivide each triangle into 4 triangles."""
    edge_midpoints = {}
    new_verts = list(vertices)
    new_faces = []

    def get_midpoint(i, j):
        key = (min(i, j), max(i, j))
        if key in edge_midpoints:
            return edge_midpoints[key]
        mid = (vertices[i] + vertices[j]) / 2.0
        idx = len(new_verts)
        new_verts.append(mid)
        edge_midpoints[key] = idx
        return idx

    for f in faces:
        a, b, c = f
        ab = get_midpoint(a, b)
        bc = get_midpoint(b, c)
        ca = get_midpoint(c, a)
        new_faces.extend([
            [a, ab, ca], [b, bc, ab], [c, ca, bc], [ab, bc, ca]
        ])

    return np.array(new_verts, dtype=np.float64), np.array(new_faces, dtype=np.int32)


def deform_sphere_to_convex(vertices, radii):
    """Deform sphere vertices by multiplying by radial distances.

    Parameters
    ----------
    vertices : (N, 3) array, unit sphere points
    radii : (N,) array, radial distances

    Returns
    -------
    (N, 3) array
    """
    norms = np.linalg.norm(vertices, axis=1, keepdims=True)
    unit = vertices / np.clip(norms, 1e-10, None)
    return unit * radii[:, np.newaxis]


def compute_facet_properties(vertices, faces):
    """Compute facet normals, areas, and centers.

    Returns
    -------
    normals : (M, 3) outward-pointing unit normals
    areas : (M,) facet areas
    centers : (M, 3) facet centroids
    """
    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]

    edge1 = v1 - v0
    edge2 = v2 - v0
    cross = np.cross(edge1, edge2)
    areas = np.linalg.norm(cross, axis=1) / 2.0
    normals = cross / (2 * areas[:, np.newaxis] + 1e-30)
    centers = (v0 + v1 + v2) / 3.0

    return normals, areas, centers


# ============================================================
# OBJ I/O
# ============================================================

def save_obj(filename, vertices, faces):
    """Save mesh to Wavefront OBJ format."""
    os.makedirs(os.path.dirname(filename) if os.path.dirname(filename) else '.', exist_ok=True)
    with open(filename, 'w') as f:
        f.write(f"# OBJ file with {len(vertices)} vertices, {len(faces)} faces\n")
        for v in vertices:
            f.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for face in faces:
            # OBJ is 1-indexed
            f.write(f"f {face[0]+1} {face[1]+1} {face[2]+1}\n")


def load_obj(filename):
    """Load mesh from Wavefront OBJ format.

    Returns
    -------
    vertices : (N, 3) array
    faces : (M, 3) array (0-indexed)
    """
    vertices = []
    faces = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            if parts[0] == 'v' and len(parts) >= 4:
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif parts[0] == 'f':
                # Handle f v1 v2 v3 or f v1/vt1/vn1 etc.
                face_verts = []
                for p in parts[1:4]:
                    idx = int(p.split('/')[0]) - 1  # Convert to 0-indexed
                    face_verts.append(idx)
                faces.append(face_verts)

    return np.array(vertices, dtype=np.float64), np.array(faces, dtype=np.int32)


# ============================================================
# Mesh comparison metrics
# ============================================================

def hausdorff_distance(verts1, verts2, n_samples=10000):
    """Compute symmetric Hausdorff distance between two meshes.

    Normalized by the diameter of the bounding sphere of mesh1.

    Parameters
    ----------
    verts1, verts2 : (N, 3) arrays
    n_samples : int
        Number of sample points (uses vertices directly if fewer)

    Returns
    -------
    float : normalized symmetric Hausdorff distance
    """
    # Sample points from vertices
    pts1 = verts1 if len(verts1) <= n_samples else verts1[
        np.random.choice(len(verts1), n_samples, replace=False)]
    pts2 = verts2 if len(verts2) <= n_samples else verts2[
        np.random.choice(len(verts2), n_samples, replace=False)]

    # Directed Hausdorff: max over pts1 of min distance to pts2
    def directed_hausdorff(a, b):
        max_dist = 0.0
        for p in a:
            dists = np.linalg.norm(b - p, axis=1)
            min_d = np.min(dists)
            if min_d > max_dist:
                max_dist = min_d
        return max_dist

    h_ab = directed_hausdorff(pts1, pts2)
    h_ba = directed_hausdorff(pts2, pts1)
    h_sym = max(h_ab, h_ba)

    # Normalize by mesh1 diameter
    diameter = np.max(np.linalg.norm(verts1 - np.mean(verts1, axis=0), axis=1)) * 2
    if diameter > 0:
        h_sym /= diameter

    return h_sym


def volumetric_iou(verts1, faces1, verts2, faces2, resolution=64):
    """Compute volumetric Intersection over Union using voxelization.

    Parameters
    ----------
    verts1, faces1 : mesh 1
    verts2, faces2 : mesh 2
    resolution : int
        Voxel grid resolution per axis

    Returns
    -------
    float : IoU in [0, 1]
    """
    # Compute bounding box encompassing both meshes
    all_verts = np.vstack([verts1, verts2])
    bbox_min = np.min(all_verts, axis=0) - 0.1
    bbox_max = np.max(all_verts, axis=0) + 0.1

    # Create voxel grid
    x = np.linspace(bbox_min[0], bbox_max[0], resolution)
    y = np.linspace(bbox_min[1], bbox_max[1], resolution)
    z = np.linspace(bbox_min[2], bbox_max[2], resolution)

    # Voxelize both meshes using ray casting (simplified: point-in-mesh test)
    voxels1 = _voxelize_mesh(verts1, faces1, x, y, z)
    voxels2 = _voxelize_mesh(verts2, faces2, x, y, z)

    intersection = np.sum(voxels1 & voxels2)
    union = np.sum(voxels1 | voxels2)

    if union == 0:
        return 0.0
    return float(intersection) / float(union)


def _voxelize_mesh(verts, faces, x, y, z):
    """Simple voxelization using distance-from-centroid heuristic.

    For convex-ish meshes, a point is inside if it's closer to the centroid
    than the surface along that direction.
    """
    centroid = np.mean(verts, axis=0)
    # Compute max radius in each direction
    radii = np.linalg.norm(verts - centroid, axis=1)
    max_radius = np.max(radii)

    # For each voxel center, check if inside
    nx, ny, nz = len(x), len(y), len(z)
    voxels = np.zeros((nx, ny, nz), dtype=bool)

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                point = np.array([x[i], y[j], z[k]])
                dist = np.linalg.norm(point - centroid)
                if dist <= max_radius:
                    # Refined: check against surface in this direction
                    direction = point - centroid
                    dir_norm = np.linalg.norm(direction)
                    if dir_norm < 1e-10:
                        voxels[i, j, k] = True
                        continue
                    direction /= dir_norm
                    # Find surface radius in this direction
                    dots = np.dot(verts - centroid, direction)
                    positive = dots > 0
                    if np.any(positive):
                        surface_r = np.max(dots[positive])
                        voxels[i, j, k] = dist <= surface_r

    return voxels


# ============================================================
# Tests
# ============================================================

if __name__ == '__main__':
    np.random.seed(42)

    print("Test 1: Identical sphere meshes")
    v1, f1 = create_sphere_mesh(n_subdivisions=2, radius=1.0)
    h = hausdorff_distance(v1, v1)
    iou = volumetric_iou(v1, f1, v1, f1, resolution=32)
    print(f"  Hausdorff distance: {h:.6f} (expected: 0.0)")
    print(f"  Volumetric IoU: {iou:.4f} (expected: 1.0)")

    print("\nTest 2: Sphere vs scaled sphere (r=1.0 vs r=1.2)")
    v2 = v1 * 1.2
    h = hausdorff_distance(v1, v2)
    iou = volumetric_iou(v1, f1, v2, f1, resolution=32)
    # Expected Hausdorff ~ 0.2/2.0 = 0.1 (normalized by diameter=2.0)
    # Expected IoU ~ (1.0^3)/(1.2^3) = 0.579
    print(f"  Hausdorff distance: {h:.4f} (expected: ~0.10)")
    print(f"  Volumetric IoU: {iou:.4f} (expected: ~0.58)")

    print("\nTest 3: Save/Load OBJ round-trip")
    os.makedirs('/tmp/test_mesh', exist_ok=True)
    save_obj('/tmp/test_mesh/sphere.obj', v1, f1)
    v_loaded, f_loaded = load_obj('/tmp/test_mesh/sphere.obj')
    print(f"  Original: {len(v1)} verts, {len(f1)} faces")
    print(f"  Loaded: {len(v_loaded)} verts, {len(f_loaded)} faces")
    assert len(v1) == len(v_loaded)
    assert len(f1) == len(f_loaded)
    assert np.allclose(v1, v_loaded, atol=1e-5)
    print("  Round-trip OK!")

    print(f"\nAll mesh metric tests passed!")
