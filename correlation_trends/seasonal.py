"""Remove the seasonal (annual) cycle with xarray ``groupby`` (worked helper).

The seasonal cycle is a large, deterministic signal. If two series share it, their
cross-correlation peaks near +/-12 months and their trends can be biased -- so it
is usually removed before correlation or trend analysis. The idiom is a
month-of-year ``groupby``: build the monthly climatology, then subtract it.
"""

from __future__ import annotations

import xarray as xr


def seasonal_climatology(da: xr.DataArray, group: str = "TIME.month") -> xr.DataArray:
    """Monthly climatology: the mean annual cycle.

    Parameters
    ----------
    da : xarray.DataArray
        Series with a datetime ``TIME`` coordinate.
    group : str, default "TIME.month"
        Grouping key. ``"TIME.month"`` gives a 12-value climatology; use
        ``"TIME.dayofyear"`` for a daily climatology.

    Returns
    -------
    xarray.DataArray
        The group-mean (e.g. 12 monthly means), indexed by the group label.
    """
    return da.groupby(group).mean()


def remove_seasonal_cycle(da: xr.DataArray, group: str = "TIME.month") -> xr.DataArray:
    """Return the series with its mean annual cycle removed, **mean preserved**.

    Subtract the monthly climatology (which removes the seasonal *departures* and
    the overall mean) and then add the overall mean back, so the deseasonalised
    series sits at the same level as the original -- only the seasonal wiggle is
    gone, not the mean.

    Parameters
    ----------
    da : xarray.DataArray
        Series with a datetime ``TIME`` coordinate.
    group : str, default "TIME.month"
        Grouping key passed to :func:`seasonal_climatology`.

    Returns
    -------
    xarray.DataArray
        ``da`` with the seasonal cycle removed, on the original ``TIME`` axis,
        retaining the original overall mean.

    Examples
    --------
    >>> clim = da.groupby("TIME.month").mean()
    >>> deseasonalised = da.groupby("TIME.month") - clim + da.mean()
    """
    # clim = seasonal_climatology(da, group)
    # then subtract it and add da.mean() back (keep the series' overall level)
    raise NotImplementedError("remove_seasonal_cycle")
