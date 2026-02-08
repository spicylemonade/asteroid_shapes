"""
Shape Comparison Metrics

Implements Hausdorff distance and Volumetric IoU for comparing
3D asteroid shape models.
"""

import numpy as np


def load_obj(filepath):
    """Load an OBJ mesh file.

    Returns
    -------
    vertices : array (N, 3)
    faces : array (M, 3) of vertex indices (0-indexed)
    """
    vertices = []
    faces = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('v '):
                parts = line.split()
                vertices.append([float(parts[1]), float(parts[2]), float(parts[3])])
            elif line.startswith('f '):
                parts = line.split()
                # OBJ faces are 1-indexed, may have v/vt/vn format
                face = []
                for p in parts[1:4]:
                    idx = int(p.split('/')[0]) - 1
                    face.append(idx)
                faces.append(face)

    return np.array(vertices, dtype=np.float64), np.array(faces, dtype=np.int64)


def sample_points_on_mesh(vertices, faces, n_samples=10000, seed=42):
    """Uniformly sample points on mesh surface.

    Uses area-weighted random sampling on triangles.
    """
    rng = np.random.RandomState(seed)

    v0 = vertices[faces[:, 0]]
    v1 = vertices[faces[:, 1]]
    v2 = vertices[faces[:, 2]]

    # Triangle areas
    cross = np.cross(v1 - v0, v2 - v0)
    areas = 0.5 * np.linalg.norm(cross, axis=1)
    total_area = np.sum(areas)

    if total_area < 1e-30:
        return vertices[:n_samples] if len(vertices) >= n_samples else vertices

    # Sample triangles proportional to area
    probs = areas / total_area
    tri_indices = rng.choice(len(faces), size=n_samples, p=probs)

    # Random barycentric coordinates
    r1 = rng.random(n_samples)
    r2 = rng.random(n_samples)
    sqrt_r1 = np.sqrt(r1)
    u = 1 - sqrt_r1
    v = sqrt_r1 * (1 - r2)
    w = sqrt_r1 * r2

    points = (u[:, None] * v0[tri_indices] +
              v[:, None] * v1[tri_indices] +
              w[:, None] * v2[tri_indices])

    return points


def hausdorff_distance(mesh1_verts, mesh1_faces, mesh2_verts, mesh2_faces,
                       n_samples=10000, seed=42):
    """Compute Hausdorff distance between two meshes.

    Returns
    -------
    dict with:
        'forward': max distance from mesh1 to mesh2
        'backward': max distance from mesh2 to mesh1
        'symmetric': max of forward and backward
        'mean_forward': mean distance from mesh1 to mesh2
        'mean_backward': mean distance from mesh2 to mesh1
    """
    pts1 = sample_points_on_mesh(mesh1_verts, mesh1_faces, n_samples, seed)
    pts2 = sample_points_on_mesh(mesh2_verts, mesh2_faces, n_samples, seed + 1)

    # Forward: for each point in mesh1, find closest in mesh2
    forward_dists = _point_to_set_distances(pts1, pts2)
    backward_dists = _point_to_set_distances(pts2, pts1)

    return {
        'forward': np.max(forward_dists),
        'backward': np.max(backward_dists),
        'symmetric': max(np.max(forward_dists), np.max(backward_dists)),
        'mean_forward': np.mean(forward_dists),
        'mean_backward': np.mean(backward_dists),
    }


def _point_to_set_distances(points_a, points_b, batch_size=1000):
    """Compute minimum distance from each point in A to closest point in B."""
    n = len(points_a)
    min_dists = np.full(n, np.inf)

    for i in range(0, n, batch_size):
        batch = points_a[i:i + batch_size]
        # Compute distances to all points in B
        diff = batch[:, None, :] - points_b[None, :, :]  # (batch, M, 3)
        dists = np.sqrt(np.sum(diff**2, axis=2))  # (batch, M)
        min_dists[i:i + batch_size] = np.min(dists, axis=1)

    return min_dists


def voxelize_mesh(vertices, faces, resolution=64):
    """Voxelize a mesh at given resolution.

    Returns a 3D boolean array where True indicates interior voxels.
    """
    # Compute bounding box
    vmin = np.min(vertices, axis=0)
    vmax = np.max(vertices, axis=0)
    extent = vmax - vmin
    max_extent = np.max(extent) * 1.1  # 10% padding

    center = (vmin + vmax) / 2.0
    half = max_extent / 2.0
    origin = center - half

    voxel_size = max_extent / resolution

    # Sample points on a grid and check inside/outside using ray casting
    grid = np.zeros((resolution, resolution, resolution), dtype=bool)

    # Use surface point sampling + filling approach
    # First mark surface voxels
    pts = sample_points_on_mesh(vertices, faces, n_samples=resolution**2 * 2)
    for pt in pts:
        ix = int((pt[0] - origin[0]) / voxel_size)
        iy = int((pt[1] - origin[1]) / voxel_size)
        iz = int((pt[2] - origin[2]) / voxel_size)
        if 0 <= ix < resolution and 0 <= iy < resolution and 0 <= iz < resolution:
            grid[ix, iy, iz] = True

    # Simple ray-casting along Z-axis for inside/outside
    # For each (x,y) column, determine inside regions
    for ix in range(resolution):
        for iy in range(resolution):
            # Find z-crossings in this column
            col = grid[ix, iy, :]
            crossings = np.where(col)[0]
            if len(crossings) >= 2:
                # Fill between first and last crossing
                grid[ix, iy, crossings[0]:crossings[-1] + 1] = True

    return grid, origin, voxel_size


def volumetric_iou(mesh1_verts, mesh1_faces, mesh2_verts, mesh2_faces,
                   resolution=64):
    """Compute Volumetric Intersection over Union between two meshes.

    Parameters
    ----------
    mesh1_verts, mesh1_faces : first mesh
    mesh2_verts, mesh2_faces : second mesh
    resolution : voxel grid resolution

    Returns
    -------
    iou : float in [0, 1]
    """
    # Align bounding boxes - use common bounding box
    all_verts = np.vstack([mesh1_verts, mesh2_verts])
    vmin = np.min(all_verts, axis=0)
    vmax = np.max(all_verts, axis=0)
    max_extent = np.max(vmax - vmin) * 1.1

    center = (vmin + vmax) / 2.0

    # Normalize both meshes to same scale
    def normalize(verts):
        c = (np.min(verts, axis=0) + np.max(verts, axis=0)) / 2.0
        scale = np.max(np.max(verts, axis=0) - np.min(verts, axis=0))
        if scale < 1e-10:
            scale = 1.0
        return (verts - c) / scale

    norm1 = normalize(mesh1_verts)
    norm2 = normalize(mesh2_verts)

    grid1, _, _ = voxelize_mesh(norm1, mesh1_faces, resolution)
    grid2, _, _ = voxelize_mesh(norm2, mesh2_faces, resolution)

    intersection = np.sum(grid1 & grid2)
    union = np.sum(grid1 | grid2)

    if union == 0:
        return 0.0

    return float(intersection) / float(union)


def compare_shapes(mesh1_path, mesh2_path, n_samples=10000, resolution=64, seed=42):
    """Compare two shape models from OBJ files.

    Returns
    -------
    dict with Hausdorff distances and volumetric IoU
    """
    v1, f1 = load_obj(mesh1_path)
    v2, f2 = load_obj(mesh2_path)

    # Normalize both to unit equivalent radius for fair comparison
    def normalize_to_unit(verts):
        center = np.mean(verts, axis=0)
        verts_c = verts - center
        mean_r = np.mean(np.linalg.norm(verts_c, axis=1))
        if mean_r > 1e-10:
            verts_c = verts_c / mean_r
        return verts_c, mean_r

    v1_norm, r1 = normalize_to_unit(v1)
    v2_norm, r2 = normalize_to_unit(v2)

    hausdorff = hausdorff_distance(v1_norm, f1, v2_norm, f2, n_samples, seed)
    iou = volumetric_iou(v1_norm, f1, v2_norm, f2, resolution)

    # Also compute hausdorff as fraction of equivalent radius
    equiv_radius = (r1 + r2) / 2.0

    return {
        'hausdorff_symmetric': hausdorff['symmetric'],
        'hausdorff_forward': hausdorff['forward'],
        'hausdorff_backward': hausdorff['backward'],
        'hausdorff_mean': (hausdorff['mean_forward'] + hausdorff['mean_backward']) / 2.0,
        'hausdorff_normalized': hausdorff['symmetric'],  # already normalized
        'volumetric_iou': iou,
        'equiv_radius_1': r1,
        'equiv_radius_2': r2,
    }
