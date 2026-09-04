"""Remove the seasonal (annual) cycle with xarray ``groupby`` (worked helper).

The seasonal cycle is a large, deterministic signal. If two series share it, their
cross-correlation peaks near +/-12 months and their trends can be biased -- so it
is usually removed before correlation or trend analysis. The idiom is a
month-of-year ``groupby``: build the monthly climatology, then subtract it.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def seasonal_climatology(
    da: xr.DataArray, group: str = "TIME.month", mode: str = "mean"
) -> xr.DataArray:
    """Monthly climatology: the mean annual cycle.

    Parameters
    ----------
    da : xarray.DataArray
        Series with a datetime ``TIME`` coordinate.
    group : str, default "TIME.month"
        Grouping key. ``"TIME.month"`` gives a 12-value climatology; use
        ``"TIME.dayofyear"`` for a daily climatology.
    mode : str, default "mean"
        Calculate seasonal cycle based on monthly mean or median.

    Returns
    -------
    xarray.DataArray
        The group-mean (e.g. 12 monthly means), indexed by the group label.
    """
    if mode == "mean":
        out = da.groupby(group).mean()
    elif mode == "median":
        out = da.groupby(group).median()
    else:
        raise ValueError
    return out


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
    clim = seasonal_climatology(da, group)
    deseasonalized = da.groupby(group) - clim + da.mean()
    return deseasonalized


def plot_seasonal_cycle(ds):
    """Plots the seasonal cycle calculated by remove_seasonal_cycle
    
        Parameters
        ----------
        da : xarray.DataArray
            Series with a datetime ``TIME`` coordinate.
        
        Returns
        -------
        Nothing
            calls matplotlib.pyplot.plot()
    """
    clim = seasonal_climatology(ds)

    months = np.arange(12) + 1
    plt.figure(figsize=[8, 4])
    plt.plot(months, clim.values)
    plt.grid()
    plt.xlabel("Month")
    plt.xticks(np.arange(12) + 1, np.arange(12) + 1)
    try:
        plt.ylabel(f"{ds.name} transport (Sv)")
    except:
        plt.ylabel("Transport (Sv)")
