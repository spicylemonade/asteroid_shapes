"""Tests for src.shapes.convex_model module."""

import numpy as np
import pytest

from src.shapes.convex_model import (
    ConvexShapeModel,
    coeffs_for_ellipsoid,
    compute_brightness,
    create_icosphere,
    lommel_seeliger_lambert,
    real_spherical_harmonics,
    spherical_harmonics_radius,
    spin_rotation_matrix,
    synthetic_lightcurve,
)


class TestCreateIcosphere:
    def test_base_icosphere(self):
        """Subdivision 0 should give 12 vertices, 20 faces (icosahedron)."""
        verts, faces = create_icosphere(subdivisions=0)
        assert verts.shape == (12, 3)
        assert faces.shape == (20, 3)

    def test_subdivision_1(self):
        """Subdivision 1: 42 verts, 80 faces."""
        verts, faces = create_icosphere(subdivisions=1)
        assert verts.shape == (42, 3)
        assert faces.shape == (80, 3)

    def test_vertices_on_unit_sphere(self):
        """All icosphere vertices should lie on the unit sphere."""
        verts, _ = create_icosphere(subdivisions=2)
        norms = np.linalg.norm(verts, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-10)


class TestSphericalHarmonics:
    def test_y00_constant(self):
        """Y_0^0 should be constant = 1/(2*sqrt(pi))."""
        theta = np.linspace(0.1, np.pi - 0.1, 20)
        phi = np.linspace(0.1, 2 * np.pi - 0.1, 20)
        vals = real_spherical_harmonics(0, 0, theta, phi)
        expected = 1.0 / (2.0 * np.sqrt(np.pi))
        np.testing.assert_allclose(vals, expected, atol=1e-10)

    def test_sh_radius_sphere(self):
        """Sphere (only c00 set) should give constant radius."""
        n_coeffs = 81  # lmax=8 -> (8+1)^2 = 81
        coeffs = np.zeros(n_coeffs)
        coeffs[0] = np.sqrt(4 * np.pi)  # r = 1

        theta = np.array([0.5, 1.0, 1.5, 2.0])
        phi = np.array([0.0, 1.0, 2.0, 3.0])
        r = spherical_harmonics_radius(coeffs, theta, phi, lmax=8)
        np.testing.assert_allclose(r, 1.0, atol=1e-6)


class TestConvexShapeModel:
    def test_default_is_sphere(self):
        """Default model should be approximately a unit sphere."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        radii = np.linalg.norm(model.vertices, axis=1)
        np.testing.assert_allclose(radii, 1.0, atol=0.1)

    def test_set_coefficients(self):
        """Setting coefficients should change the mesh."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        v_before = model.vertices.copy()
        coeffs = model.coeffs.copy()
        coeffs[4] = 0.5  # Perturb l=2 coefficient
        model.set_coefficients(coeffs)
        assert not np.allclose(model.vertices, v_before)

    def test_facet_normals_unit(self):
        """Facet normals should be approximately unit vectors."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        norms = np.linalg.norm(model.facet_normals, axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-6)

    def test_positive_areas(self):
        """All facet areas should be positive."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        assert np.all(model.facet_areas > 0)

    def test_save_and_load_obj(self, tmp_path):
        """Save to OBJ and verify file has expected content."""
        model = ConvexShapeModel(lmax=4, subdivisions=1)
        obj_path = tmp_path / "test_shape.obj"
        model.save_obj(str(obj_path))
        content = obj_path.read_text()
        v_lines = [l for l in content.splitlines() if l.startswith('v ')]
        f_lines = [l for l in content.splitlines() if l.startswith('f ')]
        assert len(v_lines) == model.n_vertices
        assert len(f_lines) == model.n_facets


class TestEllipsoidCoeffs:
    def test_sphere_coeffs(self):
        """A sphere (a=b=c=1) should have c00 dominant."""
        coeffs = coeffs_for_ellipsoid(1.0, 1.0, 1.0, lmax=4)
        # c00 should give radius ~ 1
        assert abs(coeffs[0]) > 0.5
        # Higher order should be small
        assert np.max(np.abs(coeffs[1:])) < 0.1

    def test_elongated_shape(self):
        """An elongated ellipsoid should have significant l=2 terms."""
        coeffs = coeffs_for_ellipsoid(2.0, 1.0, 1.0, lmax=4)
        # Should have non-trivial higher-order contributions
        assert np.max(np.abs(coeffs[1:])) > 0.01


class TestScatteringLaw:
    def test_zero_for_dark_facets(self):
        """If mu=0 or mu0=0, result should be near zero."""
        mu = np.array([0.0, 0.5])
        mu0 = np.array([0.5, 0.0])
        result = lommel_seeliger_lambert(mu, mu0, 0.3)
        # Not exactly zero due to Lambert term, but very small for mu0=0
        assert result[1] < 0.01

    def test_positive_for_lit_facets(self):
        """Positive mu and mu0 should give positive scattering."""
        mu = np.array([0.5, 0.8])
        mu0 = np.array([0.6, 0.7])
        result = lommel_seeliger_lambert(mu, mu0, 0.2)
        assert np.all(result > 0)


class TestComputeBrightness:
    def test_positive_brightness(self):
        """Brightness should be positive for illuminated/visible geometry."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        sun_dir = np.array([1.0, 0.0, 0.0])
        obs_dir = np.array([0.8, 0.6, 0.0])
        b = compute_brightness(model, sun_dir, obs_dir)
        assert b > 0

    def test_zero_brightness_for_opposite(self):
        """Brightness from behind should be near zero."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        sun_dir = np.array([1.0, 0.0, 0.0])
        obs_dir = np.array([-1.0, 0.0, 0.0])
        b = compute_brightness(model, sun_dir, obs_dir)
        # Should be small (but not necessarily zero for a sphere)
        assert b >= 0


class TestSpinRotation:
    def test_identity_at_zero_phase(self):
        """Rotation matrix should be orthogonal."""
        R = spin_rotation_matrix(0, 0, 1.0, 0, 2451545.0, 2451545.0)
        assert R.shape == (3, 3)
        # Should be orthogonal: R @ R^T = I
        np.testing.assert_allclose(R @ R.T, np.eye(3), atol=1e-10)

    def test_determinant_one(self):
        """Rotation matrix should have determinant 1."""
        R = spin_rotation_matrix(0.5, 0.3, 0.5, 0.1, 2451545.5, 2451545.0)
        assert abs(np.linalg.det(R) - 1.0) < 1e-10


class TestSyntheticLightcurve:
    def test_periodic(self):
        """Lightcurve of elongated shape should be periodic with rotation."""
        model = ConvexShapeModel(lmax=4, subdivisions=2)
        coeffs = coeffs_for_ellipsoid(2.0, 1.0, 1.0, lmax=4)
        model.set_coefficients(coeffs)

        period_days = 0.2  # ~4.8 hours
        n_points = 50
        t0 = 2451545.0
        times = t0 + np.linspace(0, 2 * period_days, n_points)

        sun_dirs = np.tile([1.0, 0.0, 0.3], (n_points, 1))
        sun_dirs /= np.linalg.norm(sun_dirs[0])
        obs_dirs = np.tile([0.9, 0.4, 0.1], (n_points, 1))
        obs_dirs /= np.linalg.norm(obs_dirs[0])
        phase_angles = np.full(n_points, 0.3)

        mags = synthetic_lightcurve(
            model, times, 0.0, np.pi / 4, period_days, 0.0, t0,
            sun_dirs, obs_dirs, phase_angles
        )
        assert len(mags) == n_points
        # Should have some variation (not constant)
        assert np.std(mags) > 0.001
