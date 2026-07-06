"""Spec for the xarray seasonal-cycle helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
import pytest

from correlation_trends.seasonal import seasonal_climatology, remove_seasonal_cycle


def _monthly_series(
    years: int = 20, amp: float = 3.0, mean: float = 15.0, seed: int = 0
) -> xr.DataArray:
    t = pd.date_range("2000-01-01", periods=12 * years, freq="MS")
    rng = np.random.default_rng(seed)
    annual = amp * np.cos(2 * np.pi * (t.month.values - 1) / 12)
    data = mean + annual + 0.2 * rng.standard_normal(t.size)
    return xr.DataArray(data, coords={"TIME": t}, dims="TIME")


def test_climatology_has_twelve_months_and_right_amplitude() -> None:
    da = _monthly_series(amp=3.0)
    clim = seasonal_climatology(da)
    assert clim.sizes["month"] == 12
    assert clim.max().item() - clim.min().item() == pytest.approx(6.0, abs=0.3)


def test_removing_seasonal_cycle_preserves_mean_and_removes_cycle() -> None:
    da = _monthly_series(amp=5.0, mean=15.0)
    des = remove_seasonal_cycle(da)
    # the overall mean is kept (only the seasonal departures are removed)
    assert float(des.mean()) == pytest.approx(float(da.mean()), abs=1e-6)
    # the residual monthly climatology is flat (no annual cycle left)
    resid_clim = des.groupby("TIME.month").mean()
    assert resid_clim.max().item() - resid_clim.min().item() < 1e-6
    # variance is much reduced once the big annual cycle is gone
    assert float(des.var()) < 0.2 * float(da.var())
