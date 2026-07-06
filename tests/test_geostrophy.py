"""Spec for the geostrophic-transport helpers, on the real RAPID hydrography.

These tests read ``ts_gridded.nc`` (via amocatlas / the local ``data`` dir) and
require ``gsw``. They are skipped if the data cannot be loaded (e.g. offline with
no cached file).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("gsw")
from correlation_trends.data_io import load_amoc, load_ts_gridded
from correlation_trends.geostrophy import interior_geostrophic_transport


@pytest.fixture(scope="module")
def ts():
    try:
        return load_ts_gridded()
    except Exception as exc:  # pragma: no cover - environment dependent
        pytest.skip(f"ts_gridded unavailable: {exc}")


def test_interior_transport_is_southward_and_finite(ts) -> None:
    trans = interior_geostrophic_transport(ts.isel(TIME=slice(0, 2000)))
    v = trans.values
    assert np.isfinite(v).all()
    assert trans.attrs["units"] == "Sv"
    assert np.nanmean(v) < 0  # interior flow is southward


def test_interior_transport_tracks_trans_umo(ts) -> None:
    n = 2000
    trans = interior_geostrophic_transport(ts.isel(TIME=slice(0, n)))
    _, _, series = load_amoc()
    umo = series["TRANS_UMO"][:n]
    m = np.isfinite(trans.values) & np.isfinite(umo)
    r = np.corrcoef(trans.values[m], umo[m])[0, 1]
    assert r > 0.6  # upper mid-ocean estimate tracks TRANS_UMO well
