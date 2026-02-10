"""Tests for src.inversion.genetic_solver module."""

import numpy as np
import pytest

from src.inversion.genetic_solver import (
    NonConvexMesh,
    compute_brightness_nonconvex,
)


class TestNonConvexMesh:
    def test_default_sphere(self):
        """Default mesh should be a unit sphere."""
        mesh = NonConvexMesh(subdivisions=1)
        norms = np.linalg.norm(mesh.vertices, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-10)

    def test_set_radii(self):
        """Setting radii should change vertex positions."""
        mesh = NonConvexMesh(subdivisions=1)
        new_radii = np.ones(mesh.n_vertices) * 2.0
        mesh.set_radii(new_radii)
        norms = np.linalg.norm(mesh.vertices, axis=1)
        np.testing.assert_allclose(norms, 2.0, atol=1e-10)

    def test_clamp_radii(self):
        """Radii should be clamped to [0.1, 10.0]."""
        mesh = NonConvexMesh(subdivisions=1)
        radii = np.array([0.001, 100.0] + [1.0] * (mesh.n_vertices - 2))
        mesh.set_radii(radii)
        assert mesh.radii[0] == 0.1
        assert mesh.radii[1] == 10.0

    def test_facet_normals_and_areas(self):
        """Facets should have positive areas and unit normals."""
        mesh = NonConvexMesh(subdivisions=2)
        assert np.all(mesh.facet_areas > 0)
        norms = np.linalg.norm(mesh.facet_normals, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_save_obj(self, tmp_path):
        """Should save valid OBJ file."""
        mesh = NonConvexMesh(subdivisions=1)
        path = tmp_path / "mesh.obj"
        mesh.save_obj(str(path))
        content = path.read_text()
        assert content.startswith("# Non-convex")
        v_lines = [l for l in content.splitlines() if l.startswith('v ')]
        f_lines = [l for l in content.splitlines() if l.startswith('f ')]
        assert len(v_lines) == mesh.n_vertices
        assert len(f_lines) == mesh.n_facets


class TestBrightnessNonConvex:
    def test_positive_brightness(self):
        """Illuminated sphere should have positive brightness."""
        mesh = NonConvexMesh(subdivisions=2)
        sun_dir = np.array([1.0, 0.0, 0.0])
        obs_dir = np.array([0.8, 0.6, 0.0])
        b = compute_brightness_nonconvex(mesh, sun_dir, obs_dir)
        assert b > 0

    def test_brightness_varies_with_shape(self):
        """Different shapes should give different brightness."""
        mesh1 = NonConvexMesh(subdivisions=2)
        mesh2 = NonConvexMesh(subdivisions=2)
        radii = np.ones(mesh2.n_vertices) * 1.5
        mesh2.set_radii(radii)

        sun_dir = np.array([1.0, 0.0, 0.0])
        obs_dir = np.array([0.8, 0.6, 0.0])
        b1 = compute_brightness_nonconvex(mesh1, sun_dir, obs_dir)
        b2 = compute_brightness_nonconvex(mesh2, sun_dir, obs_dir)
        assert b1 != b2
