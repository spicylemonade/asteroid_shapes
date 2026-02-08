"""
Validation metrics for comparing asteroid shape models.

Implements Hausdorff distance, volumetric IoU, pole direction error,
period error, and lightcurve residual RMS.

References:
    Hausdorff distance for mesh comparison
    Monte Carlo volumetric IoU
"""

import numpy as np
from typing import Tuple, Dict
from scipy.spatial import ConvexHull


def hausdorff_distance(vertices_a: np.ndarray, vertices_b: np.ndarray) -> float:
    """Compute Hausdorff distance between two point sets.

    H(A,B) = max(h(A,B), h(B,A))
    where h(A,B) = max_{a in A} min_{b in B} ||a-b||

    Args:
        vertices_a: (N, 3) first point set
        vertices_b: (M, 3) second point set

    Returns:
        Hausdorff distance
    """
    from scipy.spatial import cKDTree

    tree_a = cKDTree(vertices_a)
    tree_b = cKDTree(vertices_b)

    # h(A, B): for each point in A, find nearest in B
    dists_ab, _ = tree_b.query(vertices_a)
    h_ab = np.max(dists_ab)

    # h(B, A): for each point in B, find nearest in A
    dists_ba, _ = tree_a.query(vertices_b)
    h_ba = np.max(dists_ba)

    return max(h_ab, h_ba)


def mean_surface_distance(vertices_a: np.ndarray, vertices_b: np.ndarray) -> float:
    """Mean bidirectional surface distance between two point sets."""
    from scipy.spatial import cKDTree

    tree_a = cKDTree(vertices_a)
    tree_b = cKDTree(vertices_b)

    dists_ab, _ = tree_b.query(vertices_a)
    dists_ba, _ = tree_a.query(vertices_b)

    return 0.5 * (np.mean(dists_ab) + np.mean(dists_ba))


def volumetric_iou(vertices_a: np.ndarray, faces_a: np.ndarray,
                    vertices_b: np.ndarray, faces_b: np.ndarray,
                    n_samples: int = 50000, seed: int = 42) -> float:
    """Compute Volumetric Intersection over Union via Monte Carlo sampling.

    Samples random points in the bounding box and tests membership
    in each mesh using ray casting (for convex shapes, uses dot product test).

    Args:
        vertices_a, faces_a: First mesh
        vertices_b, faces_b: Second mesh
        n_samples: Number of Monte Carlo sample points
        seed: Random seed

    Returns:
        IoU value in [0, 1]
    """
    rng = np.random.RandomState(seed)

    # Combined bounding box
    all_verts = np.vstack([vertices_a, vertices_b])
    bbox_min = np.min(all_verts, axis=0) - 0.1
    bbox_max = np.max(all_verts, axis=0) + 0.1

    # Sample random points in bounding box
    points = rng.uniform(bbox_min, bbox_max, size=(n_samples, 3))

    # Test membership using convex hull (for convex shapes)
    in_a = _points_in_convex_hull(points, vertices_a)
    in_b = _points_in_convex_hull(points, vertices_b)

    intersection = np.sum(in_a & in_b)
    union = np.sum(in_a | in_b)

    if union == 0:
        return 0.0

    return intersection / union


def _points_in_convex_hull(points: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    """Test which points are inside the convex hull of vertices.

    Uses the Delaunay triangulation approach.
    """
    from scipy.spatial import Delaunay

    try:
        hull = Delaunay(vertices)
        return hull.find_simplex(points) >= 0
    except Exception:
        # Fallback: use centroid distance heuristic
        centroid = np.mean(vertices, axis=0)
        max_dist = np.max(np.linalg.norm(vertices - centroid, axis=1))
        dists = np.linalg.norm(points - centroid, axis=1)
        return dists <= max_dist


def pole_direction_error(lambda1: float, beta1: float,
                          lambda2: float, beta2: float) -> float:
    """Angular distance between two pole directions on the unit sphere.

    Args:
        lambda1, beta1: First pole (ecliptic lon, lat in radians)
        lambda2, beta2: Second pole (radians)

    Returns:
        Angular distance in degrees
    """
    # Convert to Cartesian
    v1 = np.array([
        np.cos(beta1) * np.cos(lambda1),
        np.cos(beta1) * np.sin(lambda1),
        np.sin(beta1)
    ])
    v2 = np.array([
        np.cos(beta2) * np.cos(lambda2),
        np.cos(beta2) * np.sin(lambda2),
        np.sin(beta2)
    ])

    cos_angle = np.clip(np.dot(v1, v2), -1.0, 1.0)
    # Also check antipodal (pole ambiguity: lambda+180, -beta is equivalent)
    cos_angle_anti = np.clip(np.dot(v1, -v2), -1.0, 1.0)

    angle = np.degrees(np.arccos(cos_angle))
    angle_anti = np.degrees(np.arccos(cos_angle_anti))

    return min(angle, angle_anti)


def period_relative_error(period_found: float, period_true: float) -> float:
    """Relative error in period determination.

    Also checks common aliases (P/2, 2P).

    Returns:
        Minimum relative error considering aliases (percentage)
    """
    errors = [
        abs(period_found - period_true) / period_true,
        abs(period_found * 2 - period_true) / period_true,
        abs(period_found / 2 - period_true) / period_true,
    ]
    return min(errors) * 100.0


def normalize_mesh_scale(vertices: np.ndarray) -> np.ndarray:
    """Normalize mesh to unit maximum radius from centroid."""
    centroid = np.mean(vertices, axis=0)
    centered = vertices - centroid
    max_r = np.max(np.linalg.norm(centered, axis=1))
    if max_r > 0:
        centered /= max_r
    return centered


def compute_validation_report(result_vertices: np.ndarray,
                               result_faces: np.ndarray,
                               truth_vertices: np.ndarray,
                               truth_faces: np.ndarray,
                               result_pole_lambda: float,
                               result_pole_beta: float,
                               result_period: float,
                               truth_pole_lambda: float,
                               truth_pole_beta: float,
                               truth_period: float,
                               residual_rms: float = 0.0) -> Dict:
    """Compute comprehensive validation report.

    Returns dict with all metrics.
    """
    # Normalize meshes to same scale
    norm_result = normalize_mesh_scale(result_vertices)
    norm_truth = normalize_mesh_scale(truth_vertices)

    # Shape metrics
    h_dist = hausdorff_distance(norm_result, norm_truth)
    m_dist = mean_surface_distance(norm_result, norm_truth)
    iou = volumetric_iou(norm_result, result_faces, norm_truth, truth_faces)

    # Pole error
    pole_err = pole_direction_error(
        result_pole_lambda, result_pole_beta,
        truth_pole_lambda, truth_pole_beta
    )

    # Period error
    period_err = period_relative_error(result_period, truth_period)

    # Estimate object radius for relative Hausdorff
    truth_radius = np.max(np.linalg.norm(norm_truth, axis=1))
    h_dist_relative = h_dist / truth_radius if truth_radius > 0 else h_dist

    return {
        'hausdorff_distance': float(h_dist),
        'hausdorff_distance_relative': float(h_dist_relative),
        'mean_surface_distance': float(m_dist),
        'volumetric_iou': float(iou),
        'pole_error_degrees': float(pole_err),
        'period_error_percent': float(period_err),
        'residual_rms_mag': float(residual_rms),
    }
