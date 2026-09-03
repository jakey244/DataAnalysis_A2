"""Time-domain characterisation helpers for Assignment 1.

``summary_stats``, ``seasonal_cycle`` and ``decorrelation_timescale`` characterise a
series in the time domain. Each docstring states exactly what to return.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import linregress


def summary_stats(values: np.ndarray) -> dict[str, float]:
    """Basic descriptive statistics of a (possibly gappy) series.

    Parameters
    ----------
    values : numpy.ndarray, shape (N,)
        The series, which may contain ``NaN``.

    Returns
    -------
    stats : dict[str, float]
        Dictionary with keys:

        ``n``
            Total number of samples (including gaps).
        ``n_missing``
            Number of ``NaN`` values.
        ``mean``, ``std``, ``median``, ``min``, ``max``
            Computed over the finite values (NaN-aware).
        ``range``
            ``max - min``.

    Notes
    -----
    Use the NaN-aware reductions (``numpy.nanmean`` etc.) so gaps do not poison the
    statistics. Decide and document whether ``std`` uses ``ddof=0`` or ``1``.
    """

    n = len(values)
    n_missing = int(np.isnan(values).sum().values)
    mean = np.nanmean(values)
    std = np.nanstd(values, ddof=0)  # ddof = 0
    median = np.nanmedian(values)
    min = np.nanmin(values)
    max = np.nanmax(values)
    range = max - min

    statistics = {
        "n": n,
        "n_missing": n_missing,
        "mean": mean,
        "std": std,
        "median": median,
        "min": min,
        "max": max,
        "range": range,
    }
    return statistics


def detrending(values: np.ndarray, time: np.ndarray) -> xr.DataArray:
    """
    Takes a timeseries and returns the detrended series.
    """
    d = np.asarray(values, dtype=float)
    predictor = np.arange(len(d))
    notnan = np.isfinite(d)
    d = d[notnan]
    predictor = predictor[notnan]
    
    linregress_stats = linregress(predictor, d)
    intercept = linregress_stats.intercept
    slope = linregress_stats.slope
    out = d - d.mean()
    for i in range(len(d)):
        out[i] = d[i] - intercept - slope * predictor[i]

    time = time[notnan]
    out_xarray = xr.DataArray(data=out, dims=["TIME"], coords=dict(TIME=time))
    return out_xarray


def seasonal_cycle(
    time: np.ndarray,
    values: np.ndarray,
    by: str = "month",
) -> tuple[pd.DataFrame, np.ndarray]:
    """Climatological seasonal cycle (mean and median per calendar period) and deseasonalized data.

    Parameters
    ----------
    time : numpy.ndarray, dtype ``datetime64``
        Time coordinate, same length as ``values``.
    values : numpy.ndarray, shape (N,)
        The series.
    by : {"month", "dayofyear"}, optional
        Grouping period. Default ``"month"`` (a 12-point climatology).

    Returns
    -------
    clim : pandas.DataFrame
        Indexed by the period (e.g. months 1–12), with columns ``"mean"`` and ``"median"``.
    deseasonalized : numpy.ndarray
        The input ``values`` with the seasonal mean subtracted for each period.

    Notes
    -----
    - The deseasonalized data is computed as:
      ``deseasonalized = values - clim["mean"].reindex(indexer).values``
    - The seasonal range is ``clim["mean"].max() - clim["mean"].min()``.
    """
    s = pd.Series(values, index=pd.to_datetime(time))

    # Group by the specified period
    if by == "month":
        grouped = s.groupby(s.index.month)
    elif by == "dayofyear":
        grouped = s.groupby(s.index.dayofyear)
    else:
        raise ValueError("'by' must be either 'month' or 'dayofyear'")

    # Aggregate by mean and median
    clim = grouped.agg(["mean", "median"])

    # Compute deseasonalized data
    if by == "month":
        indexer = s.index.month
    else:  # dayofyear
        indexer = s.index.dayofyear

    # Align the seasonal means with the original data
    seasonal_means = clim["mean"].reindex(indexer).values
    deseasonalized = values - seasonal_means

    return clim, deseasonalized, seasonal_means, indexer


def decorrelation_timescale(
    values: np.ndarray,
    dt: float,
) -> tuple[float, float]:
    """Integral decorrelation timescale and effective degrees of freedom.

    Parameters
    ----------
    values : numpy.ndarray, shape (N,)
        Evenly sampled series with **no gaps** (interpolate first if needed).
    dt : float
        Sampling interval (same time units you want the timescale in).

    Returns
    -------
    integral_scale : float
        Integral timescale ``tau`` (in units of ``dt``).
    ndof : float
        Effective degrees of freedom, ``N * dt / tau - 1``.

    Notes
    -----
    Integral-timescale method (after E. Frajka-Williams' ``calc_ndof.m``):

    1. Remove the mean.
    2. Form the **normalised autocovariance** ``R`` (autocorrelation, ``R(0) = 1``)
       at lags ``0, dt, 2*dt, ...`` (the non-negative lags).
    3. Integrate ``R`` from lag 0 **until its first zero crossing**, by the
       trapezoidal rule:
       ``tau = sum_i dt * (R[i] + R[i+1]) / 2`` while ``R[i] >= 0``.
       This one-sided integral (not doubled) is the integral timescale.
    4. ``ndof = N * dt / tau - 1``  (``N`` = number of samples).

    A white series gives ``tau`` of order ``dt`` and ``ndof`` close to ``N``; a
    strongly autocorrelated series gives a large ``tau`` and few independent samples.

    TODO (student): implement the autocovariance, the zero-crossing integral, and
    return ``(integral_scale, ndof)``.

    ### this function is not used here!
    """
    raise NotImplementedError("Implement the integral-timescale d.o.f. estimate.")
