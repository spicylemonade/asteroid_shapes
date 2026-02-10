"""Tests for src.geometry.ephemeris module."""

import numpy as np
import pytest

from src.geometry.ephemeris import (
    calendar_to_jd,
    compute_aspect_angles,
    compute_geometry,
    earth_position,
    kepler_position,
    solve_kepler,
    unpack_mpc_epoch,
)


class TestCalendarToJD:
    def test_j2000_epoch(self):
        """J2000.0 = 2000-01-01.5 TT = JD 2451545.0."""
        jd = calendar_to_jd(2000, 1, 1.5)
        assert abs(jd - 2451545.0) < 0.01

    def test_known_date(self):
        """Test a well-known conversion: 2024-01-01 = JD 2460310.5."""
        jd = calendar_to_jd(2024, 1, 1.5)
        assert abs(jd - 2460310.5) < 1.0  # within 1 day tolerance


class TestUnpackMPCEpoch:
    def test_k25a5(self):
        """K25A5 = 2025-Jan-05."""
        jd = unpack_mpc_epoch("K25A5")
        expected = calendar_to_jd(2025, 10, 5)
        assert abs(jd - expected) < 1.0

    def test_invalid_returns_default(self):
        jd = unpack_mpc_epoch("???")
        assert jd == 2451545.0

    def test_j2000_packed(self):
        """K00A1 = 2000-Jan-01."""
        jd = unpack_mpc_epoch("K00A1")
        expected = calendar_to_jd(2000, 10, 1)
        assert abs(jd - expected) < 1.0


class TestSolveKepler:
    def test_circular_orbit(self):
        """For e=0, E should equal M."""
        M = 1.5
        E = solve_kepler(M, 0.0)
        assert abs(E - M) < 1e-10

    def test_moderate_eccentricity(self):
        """Verify M = E - e*sin(E) for e=0.3."""
        M = 2.0
        e = 0.3
        E = solve_kepler(M, e)
        assert abs(M - (E - e * np.sin(E))) < 1e-10

    def test_high_eccentricity(self):
        """Check convergence for e=0.9."""
        M = 0.5
        e = 0.9
        E = solve_kepler(M, e)
        assert abs(M - (E - e * np.sin(E))) < 1e-10


class TestKeplerPosition:
    def test_circular_orbit_distance(self):
        """Circular orbit should have constant distance = a."""
        a, e = 2.0, 0.0
        for M in [0, 1, 2, 3, 4, 5]:
            pos = kepler_position(a, e, 0, 0, 0, M)
            assert abs(np.linalg.norm(pos) - a) < 1e-8

    def test_perihelion_distance(self):
        """At M=0 (perihelion), distance should be a*(1-e)."""
        a, e = 2.5, 0.2
        pos = kepler_position(a, e, 0, 0, 0, 0.0)
        expected_r = a * (1 - e)
        assert abs(np.linalg.norm(pos) - expected_r) < 1e-6


class TestEarthPosition:
    def test_distance_approximately_1au(self):
        """Earth should be approximately 1 AU from the Sun."""
        jd = 2451545.0  # J2000
        pos = earth_position(jd)
        r = np.linalg.norm(pos)
        assert 0.98 < r < 1.02

    def test_different_epochs(self):
        """Earth position should change over time."""
        pos1 = earth_position(2451545.0)
        pos2 = earth_position(2451545.0 + 180)  # ~6 months later
        dist = np.linalg.norm(pos1 - pos2)
        assert dist > 1.0  # Should be roughly 2 AU apart


class TestComputeGeometry:
    def test_basic_geometry(self):
        """Test geometry computation with a simple asteroid record."""
        rec = {
            'a': 1.5, 'e': 0.2, 'incl': 0.1, 'node': 0.5,
            'peri': 1.0, 'M0': 0.5, 'n': 0.01, 'epoch_jd': 2451545.0,
        }
        result = compute_geometry(rec, 2451545.0)
        assert 'phase_angle' in result
        assert 'solar_elongation' in result
        assert 0 <= result['phase_angle'] <= 180
        assert 0 <= result['solar_elongation'] <= 180
        assert result['r_helio'] > 0
        assert result['delta_geo'] > 0

    def test_sun_dir_unit_vector(self):
        rec = {
            'a': 2.0, 'e': 0.1, 'incl': 0.05, 'node': 0.3,
            'peri': 0.8, 'M0': 1.0, 'n': 0.005, 'epoch_jd': 2451545.0,
        }
        result = compute_geometry(rec, 2451600.0)
        sun_norm = np.linalg.norm(result['sun_dir_ecl'])
        obs_norm = np.linalg.norm(result['obs_dir_ecl'])
        assert abs(sun_norm - 1.0) < 1e-10
        assert abs(obs_norm - 1.0) < 1e-10


class TestComputeAspectAngles:
    def test_pole_aligned_with_z(self):
        """If pole is at ecliptic north, sub-observer lat should match geometry."""
        sun_dir = np.array([1, 0, 0], dtype=float)
        obs_dir = np.array([0, 0, 1], dtype=float)
        result = compute_aspect_angles(sun_dir, obs_dir,
                                       pole_lambda=0.0, pole_beta=90.0)
        assert abs(result['sub_obs_lat'] - 90.0) < 1.0
