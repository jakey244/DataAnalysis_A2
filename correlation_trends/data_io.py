"""Load AMOC datasets via ``amocatlas`` (worked helpers).

- :func:`load_amoc` — RAPID 26N transports (MOC and components).
- :func:`load_ts_gridded` — RAPID 26N gridded boundary T/S (for geostrophy).
- :func:`load_47n` — NOAC 47N basin-wide AMOC transport (second dataset).
"""

from __future__ import annotations

import numpy as np
import xarray as xr
from numpy.typing import NDArray


def load_amoc() -> tuple[NDArray, float, dict[str, NDArray[np.float64]]]:
    """Load the standardised RAPID 26N transports via ``amocatlas``.

    Returns
    -------
    time : numpy.ndarray
        ``datetime64`` time axis.
    dt : float
        Median sample spacing, in days.
    series : dict of str to numpy.ndarray
        ``MOC``, ``TRANS_FC`` (Gulf Stream / Florida Current), ``TRANS_EKMAN``,
        and ``TRANS_UMO`` (upper-mid-ocean geostrophic), each gap-filled by
        linear interpolation over non-finite values.

    Notes
    -----
    Requires the ``amocatlas`` package and network access on first use.
    """
    from amocatlas import read

    ds = read.rapid()  # standardised names on the TIME coordinate
    t = ds["TIME"].values
    dt = float(np.median(np.diff(t)) / np.timedelta64(1, "D"))

    def v(name: str) -> NDArray[np.float64]:
        a = ds[name].values.astype(float)
        good = np.isfinite(a)
        idx = np.arange(a.size)
        a[~good] = np.interp(idx[~good], idx[good], a[good])
        return a

    return t, dt, {k: v(k) for k in ("MOC", "TRANS_FC", "TRANS_EKMAN", "TRANS_UMO")}


def load_ts_gridded(data_dir: str | None = "data") -> xr.Dataset:
    """Load the standardised RAPID 26N gridded T/S (``ts_gridded.nc``).

    Parameters
    ----------
    data_dir : str or None, default "data"
        Directory holding ``ts_gridded.nc``. If the file is present it is used
        directly; otherwise ``amocatlas`` downloads it (~485 MB).

    Returns
    -------
    xarray.Dataset
        Standardised dataset with ``TEMP_WEST``, ``PSAL_WEST``, ``TEMP_EAST``,
        ``PSAL_EAST`` (and the MAR boundaries) on ``PRESSURE`` x ``TIME``.
    """
    from amocatlas import read

    ds = read.rapid(
        file_list=["ts_gridded.nc"], transport_only=False, data_dir=data_dir
    )
    return ds[0] if isinstance(ds, list) else ds


def load_moc_sigma0_26n(
    data_dir: str | None = "data",
) -> tuple[NDArray, NDArray[np.float64]]:
    """Load RAPID 26N AMOC strength in density (sigma0) coordinates, ``MOC_SIGMA0``.

    Parameters
    ----------
    data_dir : str or None, default "data"
        Directory holding ``meridional_transports.nc`` (downloaded if absent).

    Returns
    -------
    time : numpy.ndarray
        ``datetime64`` time axis (~10-daily, 2004-2024).
    moc_sigma0 : numpy.ndarray
        AMOC strength (Sv) in sigma0 coordinates at 26N.
    """
    from amocatlas import read

    ds = read.rapid(
        file_list=["meridional_transports.nc"], transport_only=False, data_dir=data_dir
    )
    ds = ds[0] if isinstance(ds, list) else ds
    return ds["TIME"].values, ds["MOC_SIGMA0"].values.astype(float)


def load_47n() -> tuple[NDArray, NDArray[np.float64]]:
    """Load the NOAC 47N basin-wide AMOC volume transport, ``MOC_SIGMA0``.

    Returns
    -------
    time : numpy.ndarray
        ``datetime64`` time axis (monthly, 1993-2018).
    moc_sigma0 : numpy.ndarray
        AMOC volume transport (Sv) at 47N in sigma0 coordinates.

    Notes
    -----
    Requires network access on first use (``amocatlas.read.read_47n``). The raw
    column ``"Trans vol [Sv]"`` is the standardised ``MOC_SIGMA0`` (unfiltered).
    """
    from amocatlas import read

    out = read.read_47n()
    ds = out[0] if isinstance(out, (list, tuple)) else out
    return ds["TIME"].values, ds["Trans vol [Sv]"].values.astype(float)
