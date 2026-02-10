"""Tests for src.inversion.period_search module."""

import numpy as np
import pytest

from src.inversion.period_search import (
    find_period,
    lomb_scargle_period,
    phase_dispersion_minimization,
)


class TestLombScargle:
    def test_detects_known_period(self):
        """Should detect a strong sinusoidal signal."""
        np.random.seed(42)
        true_period = 5.0
        times = np.sort(np.random.uniform(0, 50, 200))
        signal = np.sin(2 * np.pi * times / true_period)
        mags = signal + np.random.normal(0, 0.1, len(times))

        periods, powers, _, _ = lomb_scargle_period(
            times, mags, period_min=2.0, period_max=20.0, n_frequencies=10000)

        # The best period should be close to true period
        assert abs(periods[0] - true_period) < 0.5

    def test_returns_correct_shapes(self):
        times = np.linspace(0, 10, 50)
        mags = np.sin(times) + 0.01 * np.random.randn(50)
        periods, powers, all_p, all_pow = lomb_scargle_period(
            times, mags, n_frequencies=1000)
        assert len(periods) == len(powers) == 1000
        assert len(all_p) == len(all_pow) == 1000


class TestPDM:
    def test_detects_known_period(self):
        """PDM should find minimum theta near true period."""
        np.random.seed(42)
        true_period = 8.0
        times = np.sort(np.random.uniform(0, 80, 150))
        signal = np.sin(2 * np.pi * times / true_period)
        mags = signal + np.random.normal(0, 0.1, len(times))

        periods, theta, _, _ = phase_dispersion_minimization(
            times, mags, period_min=2.0, period_max=20.0, n_periods=5000)

        # Best period should be near the true period
        assert abs(periods[0] - true_period) < 1.0

    def test_constant_signal(self):
        """Constant magnitudes should give theta ~ 1.0 everywhere."""
        times = np.linspace(0, 10, 50)
        mags = np.ones(50) * 15.0
        periods, theta, _, _ = phase_dispersion_minimization(
            times, mags, n_periods=100)
        # All theta should be 1.0 since there's no variance
        np.testing.assert_allclose(theta, 1.0, atol=0.01)


class TestFindPeriod:
    def test_returns_candidates(self):
        """find_period should return a list of candidate dicts."""
        np.random.seed(42)
        times = np.sort(np.random.uniform(0, 20, 100))
        mags = np.sin(2 * np.pi * times / 4.0) + 0.1 * np.random.randn(100)

        candidates = find_period(times, mags, period_min=2.0, period_max=20.0, n_top=5)
        assert len(candidates) > 0
        assert len(candidates) <= 5
        for c in candidates:
            assert 'period' in c
            assert 'score' in c
            assert 'method' in c
            assert c['period'] > 0

    def test_no_duplicates(self):
        """Candidates should not have near-duplicate periods."""
        np.random.seed(42)
        times = np.sort(np.random.uniform(0, 30, 100))
        mags = np.sin(2 * np.pi * times / 6.0) + 0.1 * np.random.randn(100)

        candidates = find_period(times, mags, n_top=10)
        periods = [c['period'] for c in candidates]
        for i in range(len(periods)):
            for j in range(i + 1, len(periods)):
                # No two periods within 0.5% of each other
                assert abs(periods[i] - periods[j]) >= 0.005 * min(periods[i], periods[j])
