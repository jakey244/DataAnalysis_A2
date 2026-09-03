"""Interior geostrophic (thermal-wind) transport from boundary T/S profiles.

This ports the Lecture 3 method: convert practical T/S to TEOS-10, form the
dynamic height (geostrophic streamfunction) at the western and eastern
boundaries, difference them for the zonally-integrated geostrophic transport per
unit depth, and integrate over depth. Applied to the RAPID 26N gridded
hydrography (``ts_gridded.nc`` via :func:`correlation_trends.data_io.load_ts_gridded`).

Requires ``gsw`` (TEOS-10).

Notes
-----
The estimate here is the *baroclinic interior* transport between the two
boundaries. It captures the sign and gross variability of RAPID's published
``TRANS_UMO`` but is not identical: the official product also splits the interior
at the Mid-Atlantic Ridge and applies a mass-balance (external-transport)
adjustment, which removes variance this simple estimate retains.
"""

from __future__ import annotations

import os

import gsw
import numpy as np
import xarray as xr
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid

# (temperature var, salinity var, representative longitude) for each boundary
WEST = ("TEMP_WEST", "PSAL_WEST", -76.74)
EAST = ("TEMP_EAST", "PSAL_EAST", -16.23)


def to_teos10(
    temp: ArrayLike, salt: ArrayLike, pres: ArrayLike, lon: float, lat: float
) -> tuple[NDArray, NDArray]:
    """Convert practical temperature/salinity to TEOS-10 SA and CT.

    Parameters
    ----------
    temp : array_like
        In-situ temperature (degC).
    salt : array_like
        Practical salinity (PSS-78).
    pres : array_like
        Pressure (dbar), broadcastable to ``temp``.
    lon, lat : float
        Location, for the salinity composition correction.

    Returns
    -------
    SA : numpy.ndarray
        Absolute salinity (g/kg).
    CT : numpy.ndarray
        Conservative temperature (degC).
    """
    n_time = len(temp[0, :])
    n_pres = len(temp[:, 0])
    SA = np.empty([n_time, n_pres])
    CT = np.empty([n_time, n_pres])
    for i in range(n_time):
        SA[i] = gsw.SA_from_SP(
            np.asarray(salt[:, i], float), np.asarray(pres, float), lon, lat
        )
        CT[i] = gsw.CT_from_t(
            SA[i], np.asarray(temp[:, i], float), np.asarray(pres, float)
        )
    # SA = gsw.SA_from_SP(np.asarray(salt, float), np.asarray(pres, float), lon, lat)
    # CT = gsw.CT_from_t(SA, np.asarray(temp, float), np.asarray(pres, float))
    return SA, CT


def dynamic_height(
    SA: ArrayLike, CT: ArrayLike, pres: ArrayLike, p_ref: float = 4820.0
) -> NDArray:
    """Geostrophic streamfunction (dynamic height anomaly), referenced to ``p_ref``.

    Parameters
    ----------
    SA, CT : array_like
        Absolute salinity and conservative temperature, pressure along axis 0.
    pres : array_like
        Pressure (dbar), same shape as ``SA`` (or broadcastable), increasing
        along axis 0.
    p_ref : float, default 4820.0
        Reference pressure (assumed level of no motion).

    Returns
    -------
    numpy.ndarray
        Dynamic height anomaly (m2/s2), zero at ``p_ref``.
    """
    return gsw.geo_strf_dyn_height(
        np.asarray(SA, float),
        np.asarray(CT, float),
        np.asarray(pres, float),
        p_ref=p_ref,
        axis=0,
    )


def _fill_profiles(a: NDArray, p: NDArray) -> NDArray:
    """Linear-interpolate NaNs along pressure (axis 0), per time column."""
    a = np.array(a, float)
    for j in range(a.shape[1]):
        col = a[:, j]
        good = np.isfinite(col)
        if good.sum() >= 2:
            a[:, j] = np.interp(p, p[good], col[good])
    return a


def interior_geostrophic_transport(
    ds: xr.Dataset,
    lat: float = 26.5,
    p_ref: float = 4820.0,
    z_max: float = 1100.0,
    west: tuple = WEST,
    east: tuple = EAST,
) -> xr.DataArray:
    """Upper mid-ocean geostrophic transport from boundary T/S profiles.

    Follows the RAPID definition (McCarthy et al., 2015): the mid-ocean
    geostrophic transport per unit depth, ``(phi_east - phi_west) / f`` referenced
    to ``p_ref``, integrated **from the surface down to the depth of the AMOC
    maximum** (``z_max``, ~1100 m when northward AAIW is present). Integrating the
    full water column instead would fold in the deep NADW/AABW and is *not* the
    upper mid-ocean transport.

    Parameters
    ----------
    ds : xarray.Dataset
        Standardised ``ts_gridded`` dataset (``TEMP_WEST``, ``PSAL_WEST``,
        ``TEMP_EAST``, ``PSAL_EAST`` on ``PRESSURE`` x ``TIME``).
    lat : float, default 26.5
        Latitude for the Coriolis parameter.
    p_ref : float, default 4820.0
        Reference pressure (level of no motion), the deepest common level.
    z_max : float, default 1100.0
        Depth (m) of the upper integration limit -- the depth of the AMOC
        maximum. RAPID uses a time-varying value (~1100 m with AAIW, ~700 m
        without); a fixed value is a good approximation for teaching.
    west, east : tuple
        ``(temp_var, salt_var, longitude)`` for each boundary.

    Returns
    -------
    xarray.DataArray
        Transport (Sv) on the ``TIME`` coordinate. Negative is southward.

    Notes
    -----
    The zonally-integrated geostrophic transport per unit depth is
    ``(phi_east - phi_west) / f`` (the horizontal distance cancels between the
    thermal-wind velocity and its zonal integral). This single east-west estimate
    still differs from the published ``TRANS_UMO`` (correlation ~0.77): the
    official product also splits the interior at the Mid-Atlantic Ridge, applies a
    mass-balance adjustment, and uses the time-varying AMOC-maximum depth.
    """
    p = ds["PRESSURE"].values.astype(float)

    def boundary_phi(spec: tuple) -> NDArray:
        tname, sname, lon = spec
        T = _fill_profiles(ds[tname].values, p)
        S = _fill_profiles(ds[sname].values, p)
        p2 = p[:, None] * np.ones_like(T)
        SA, CT = to_teos10(T, S, p2, lon, lat)
        return dynamic_height(SA, CT, p2, p_ref=p_ref)

    phi_w = boundary_phi(west)
    phi_e = boundary_phi(east)
    f = gsw.f(lat)
    transport_per_depth = (phi_e - phi_w) / f  # m2/s
    depth = -gsw.z_from_p(p, lat)  # m, increasing downward
    upper = depth <= z_max  # surface -> AMOC-max depth
    trapezoid = getattr(np, "trapezoid", np.trapz)  # numpy>=2 renamed trapz
    transport = trapezoid(transport_per_depth[upper], depth[upper], axis=0) / 1e6

    return xr.DataArray(
        transport,
        coords={"TIME": ds["TIME"].values},
        dims="TIME",
        name="interior_geostrophic_transport",
        attrs={
            "units": "Sv",
            "long_name": "Upper mid-ocean geostrophic transport",
            "reference_pressure": p_ref,
            "integration_depth_m": z_max,
        },
    )


def geostrophy_jakob(
    ts,
    pref=4820.0,
    z_max=1100.0,
    lat=26.5,
    west_lon=-76.74,
    east_lon=16.23,
    i_want_to_use_saved_files=True,
) -> xr.DataArray:
    """Upper mid-ocean geostrophic transport from boundary T/S profiles.

    Follows the RAPID definition (McCarthy et al., 2015): the mid-ocean
    geostrophic transport per unit depth, ``(phi_east - phi_west) / f`` referenced
    to ``p_ref``, integrated **from the surface down to the depth of the AMOC
    maximum** (``z_max``, ~1100 m when northward AAIW is present). Integrating the
    full water column instead would fold in the deep NADW/AABW and is *not* the
    upper mid-ocean transport.

    Parameters
    ----------
    ds : xarray.Dataset
        Standardised ``ts_gridded`` dataset (``TEMP_WEST``, ``PSAL_WEST``,
        ``TEMP_EAST``, ``PSAL_EAST`` on ``PRESSURE`` x ``TIME``).
    pref : float, default 4820.0
        Reference pressure (level of no motion), the deepest common level.
    lat : float, default 26.5
        Latitude for the Coriolis parameter.
    z_max : float, default 1100.0
        Depth (m) of the upper integration limit -- the depth of the AMOC
        maximum. RAPID uses a time-varying value (~1100 m with AAIW, ~700 m
        without); a fixed value is a good approximation for teaching.
    west_lon, east_lon : float
        longitude coordinates of boundaries, default valid for RAPID (26°N)
    i_want_to_use_saved_files : bool, default True
        whether to use/write .npy files after the TEOS conversion,
        TEOS conversion of the whole timeseries and vertical column takes about 5 mins
        (on my local machine)

    Returns
    -------
    xarray.DataArray
        Transport (Sv) on the ``TIME`` coordinate. Negative is southward.
    """
    # define the vertical pressure levels
    p = ts.PRESSURE.values

    # define the coriolis parameter at given latitude
    f = gsw.f(lat)

    # converting to TEOS-10 units
    # check if the data files exist and read them, if they dont exist,
    # convert and write .npy files for them to make the calculation quicker (takes around 5 mins for me locally)
    if (
        os.path.isfile("../data/TS_WEST.npy")
        and os.path.isfile("../data/TS_EAST.npy")
        and i_want_to_use_saved_files
    ):
        TS_WEST = np.load("../data/TS_WEST.npy")
        TS_EAST = np.load("../data/TS_EAST.npy")
    else:
        # defining salinity and temperature arrays at the boundaries for further calculation
        PSAL_WEST = _fill_profiles(ts.PSAL_WEST.values, p)
        PSAL_EAST = _fill_profiles(ts.PSAL_EAST.values, p)
        T_WEST = _fill_profiles(ts.TEMP_WEST.values, p)
        T_EAST = _fill_profiles(ts.TEMP_EAST.values, p)
        # create data directory in the repo if it doesnt exist (not necessary with my gitignore settings but still cleaner)
        os.makedirs("../data", exist_ok=True)
        TS_WEST = to_teos10(temp=T_WEST, salt=PSAL_WEST, pres=p, lat=lat, lon=west_lon)
        TS_EAST = to_teos10(temp=T_EAST, salt=PSAL_EAST, pres=p, lat=lat, lon=east_lon)
        if i_want_to_use_saved_files:
            np.save("../data/TS_WEST.npy", TS_WEST)
            np.save("../data/TS_EAST.npy", TS_EAST)

    # calculate the dynamic heights from converted profiles
    dynamic_height_west_nan = dynamic_height(
        TS_WEST[0].T, TS_WEST[1].T, p, pref
    )  # have to transpose the arrays, otherwise the pressure is on the wrong axis
    dynamic_height_east_nan = dynamic_height(TS_EAST[0].T, TS_EAST[1].T, p, pref)
    # linearly interpolate the dynamic height profile time series
    dynamic_height_west = _fill_profiles(dynamic_height_west_nan, p)
    dynamic_height_east = _fill_profiles(dynamic_height_east_nan, p)

    psi_difference = dynamic_height_east - dynamic_height_west

    # zonally-integrated geostrophic transport per unit depth, m^2/s
    thermal_wind = psi_difference / f

    # integrate until reference pressure
    depth = -gsw.z_from_p(p, lat)  # m, increasing downward
    relevant_depths = depth <= z_max
    transport = trapezoid(
        y=thermal_wind[relevant_depths, :], x=depth[relevant_depths], axis=0
    )
    transport = transport / 1e6  # converted to Sverdrup unit

    return xr.DataArray(
        transport,
        coords={"TIME": ts["TIME"].values},
        dims="TIME",
        name="interior_geostrophic_transport",
        attrs={
            "units": "Sv",
            "long_name": "Upper mid-ocean geostrophic transport",
            "reference_pressure": pref,
            "integration_depth_m": z_max,
        },
    )
