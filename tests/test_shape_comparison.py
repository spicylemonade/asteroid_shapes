"""Tests for src.metrics.shape_comparison module."""

import numpy as np
import pytest

from src.metrics.shape_comparison import (
    hausdorff_distance,
    load_obj,
    sample_points_on_mesh,
    volumetric_iou,
)
from src.shapes.convex_model import create_icosphere


class TestSamplePoints:
    def test_correct_number_of_samples(self):
        verts, faces = create_icosphere(2)
        pts = sample_points_on_mesh(verts, faces, n_samples=500)
        assert pts.shape == (500, 3)

    def test_points_on_unit_sphere(self):
        """Sampled points on a unit icosphere should have norm ~ 1."""
        verts, faces = create_icosphere(3)
        pts = sample_points_on_mesh(verts, faces, n_samples=1000)
        norms = np.linalg.norm(pts, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=0.05)


class TestHausdorffDistance:
    def test_identical_meshes(self):
        """Hausdorff of identical meshes should have small mean distances."""
        verts, faces = create_icosphere(2)
        result = hausdorff_distance(verts, faces, verts, faces, n_samples=2000)
        # Mean distances should be very small; max (Hausdorff) has sampling noise
        assert result['mean_forward'] < 0.1
        assert result['mean_backward'] < 0.1

    def test_scaled_meshes(self):
        """Scaling one mesh should give bounded Hausdorff."""
        verts, faces = create_icosphere(2)
        verts_scaled = verts * 1.1
        result = hausdorff_distance(verts, faces, verts_scaled, faces, n_samples=1000)
        # The Hausdorff should be approximately 0.1 (the scaling offset)
        assert 0.05 < result['symmetric'] < 0.3

    def test_symmetry(self):
        """Hausdorff should be symmetric."""
        verts1, faces1 = create_icosphere(2)
        verts2 = verts1 * 1.2
        r1 = hausdorff_distance(verts1, faces1, verts2, faces1, n_samples=500)
        r2 = hausdorff_distance(verts2, faces1, verts1, faces1, n_samples=500)
        assert abs(r1['symmetric'] - r2['symmetric']) < 0.05


class TestVolumetricIoU:
    def test_identical_meshes(self):
        """IoU of identical meshes should be ~1.0."""
        verts, faces = create_icosphere(2)
        iou = volumetric_iou(verts, faces, verts, faces, resolution=32)
        assert iou > 0.8

    def test_different_shapes_lower_iou(self):
        """Very different shapes should have lower IoU than identical ones."""
        verts, faces = create_icosphere(2)
        # Create a highly elongated version
        verts_elongated = verts.copy()
        verts_elongated[:, 0] *= 3.0  # stretch along x by 3x
        iou_same = volumetric_iou(verts, faces, verts, faces, resolution=32)
        iou_diff = volumetric_iou(verts, faces, verts_elongated, faces, resolution=32)
        assert iou_diff < iou_same


class TestLoadObj:
    def test_load_obj(self, tmp_path):
        """Load a simple OBJ file."""
        obj_content = (
            "v 0.0 0.0 0.0\n"
            "v 1.0 0.0 0.0\n"
            "v 0.0 1.0 0.0\n"
            "v 0.0 0.0 1.0\n"
            "f 1 2 3\n"
            "f 1 2 4\n"
            "f 1 3 4\n"
            "f 2 3 4\n"
        )
        obj_path = tmp_path / "test.obj"
        obj_path.write_text(obj_content)
        verts, faces = load_obj(str(obj_path))
        assert verts.shape == (4, 3)
        assert faces.shape == (4, 3)
        assert faces.min() == 0  # 0-indexed
