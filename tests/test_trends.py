"""Spec for ``fit_trend`` and the ``trend_with_significance`` stub."""

from __future__ import annotations

import numpy as np
import pytest

from correlation_trends.trends import fit_trend, trend_with_significance


def _ar1(n: int, r: float, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    a = np.zeros(n)
    a[0] = rng.standard_normal()
    for i in range(1, n):
        a[i] = r * a[i - 1] + rng.standard_normal()
    return a


def test_fit_trend_recovers_planted_slope() -> None:
    t = np.arange(500.0)
    x = 0.5 * t + 3 + np.random.default_rng(0).standard_normal(500)
    slope, intercept = fit_trend(t, x)
    assert slope == pytest.approx(0.5, abs=0.01)
    assert intercept == pytest.approx(3.0, abs=1.0)


def test_trend_significance_recovers_slope() -> None:
    t = np.arange(1000.0)
    x = 0.02 * t + _ar1(1000, 0.8, seed=1)
    res = trend_with_significance(t, x, dt=1.0)
    assert res.slope == pytest.approx(0.02, abs=0.01)


def test_autocorrelation_inflates_the_standard_error() -> None:
    # red-noise residuals => honest SE larger, honest p larger, N_eff < N
    t = np.arange(2000.0)
    x = 0.005 * t + _ar1(2000, 0.9, seed=2)
    res = trend_with_significance(t, x, dt=1.0)
    assert res.se_eff > res.se
    assert res.p_eff > res.p_naive
    assert res.n_eff < len(x)
    # slope/SE ratios: naive over-confident (larger magnitude) than honest
    assert abs(res.t_naive) > abs(res.t_eff)


def test_slope_over_se_matches_ratio() -> None:
    t = np.arange(800.0)
    x = 0.01 * t + _ar1(800, 0.7, seed=3)
    res = trend_with_significance(t, x, dt=1.0)
    assert res.t_naive == pytest.approx(res.slope / res.se)
    assert res.t_eff == pytest.approx(res.slope / res.se_eff)
