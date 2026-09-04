"""Linear trends and their significance for autocorrelated series.

Worked helper: :func:`fit_trend`. Student stub: :func:`trend_with_significance`.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike
from scipy import stats

from .correlation import effective_dof


class TrendResult(NamedTuple):
    """Result of :func:`trend_with_significance`.

    Attributes
    ----------
    slope, intercept : float
        Least-squares fit coefficients (``x ~ slope * t + intercept``).
    se : float
        Naive OLS standard error of the slope (assumes independent residuals).
    se_eff : float
        Standard error inflated for autocorrelation, ``se * sqrt(N / N_eff)``.
    p_naive, p_eff : float
        Two-sided p-values for ``slope = 0`` using ``se`` (with ``N - 2`` d.o.f.)
        and ``se_eff`` (with ``N_eff - 2`` d.o.f.) respectively.
    n_eff : float
        Effective sample size ``N / (1 + 2 sum rho_k)`` from the residuals.
    t_naive, t_eff : float
        Slope in units of its standard error (``slope / se`` and ``slope / se_eff``)
        -- how many sigma the slope sits from zero. Significant at 95% needs
        ``|t| > ~1.96``.
    """

    slope: float
    intercept: float
    se: float
    se_eff: float
    p_naive: float
    p_eff: float
    n_eff: float
    t_naive: float
    t_eff: float


def fit_trend(t: ArrayLike, x: ArrayLike) -> tuple[float, float]:
    """Least-squares straight-line fit.

    Parameters
    ----------
    t : array_like
        Predictor (e.g. time).
    x : array_like
        Response series.

    Returns
    -------
    slope, intercept : float
        Coefficients of ``x ~ slope * t + intercept``.
    """
    slope, intercept = np.polyfit(np.asarray(t, float), np.asarray(x, float), 1)
    return float(slope), float(intercept)


def trend_with_significance(t: ArrayLike, x: ArrayLike, dt: float) -> TrendResult:
    """Slope, its standard error, ``slope/SE``, and autocorrelation-aware p-values.

    The naive standard error assumes independent residuals. For a red-noise
    series that is optimistic: neighbouring residuals are correlated, so the
    record carries fewer independent samples than ``N``. The honest version
    replaces ``N`` by the effective sample size ``N_eff`` derived from the
    two-sided integral timescale of the residuals, inflating the SE by
    ``sqrt(N / N_eff)``.

    Parameters
    ----------
    t : array_like
        Predictor (time), same units implied by ``dt``.
    x : array_like
        Response series.
    dt : float
        Sample spacing (present for interface symmetry; ``N_eff`` here is a
        sample count and does not depend on ``dt``).

    Returns
    -------
    TrendResult
        Slope, intercept, naive and effective SE, p-values, ``N_eff``, and the
        ``slope/SE`` ratios ``t_naive`` and ``t_eff``.

    Notes
    -----
    ``sigma`` in ``slope/SE`` is the standard error *of the slope itself* -- how
    much the fitted slope would wobble on resampling -- not the scatter of the
    data about the line. Significant at 95% means ``|slope| > ~1.96 * SE``.
    """
    try:
        x = x.values
    except:
        print("Not xarray format")
    notnan = np.isfinite(x)
    x = x[notnan]
    t = t[notnan]
    N = x.size
    slope, intercept = fit_trend(t, x)
    # print(slope,intercept)
    resid = x - (slope * t + intercept)
    # print(resid)
    N_eff = effective_dof(resid, dt)
    stderr = np.sqrt(np.sum(resid**2) / (N - 2) / np.sum((t - t.mean()) ** 2))
    effective_stderr = stderr * np.sqrt(N / N_eff)
    t_naive = slope / stderr
    t_eff = slope / effective_stderr
    p_naive = 2 * stats.t.sf(np.abs(t_naive), N_eff * 2)
    p_eff = 2 * stats.t.sf(np.abs(t_eff), N_eff)

    res = TrendResult(
        slope=slope,
        intercept=intercept,
        se=stderr,
        se_eff=effective_stderr,
        p_naive=p_naive,
        p_eff=p_eff,
        n_eff=N_eff,
        t_naive=t_naive,
        t_eff=t_eff,
    )
    return res


def confband(
    t: ArrayLike, x: ArrayLike, dt: float, alpha: float = 0.05
) -> tuple[np.ndarray, np.ndarray]:
    trend_result = trend_with_significance(t, x, dt)
    slope = trend_result.slope
    intercept = trend_result.intercept
    slope_se = trend_result.se_eff  # standard error of the slope
    EDOF = trend_result.n_eff
    n = t.size
    t_mean = t.mean()
    ss_t = ((t - t_mean) ** 2).sum()

    x_fit = intercept + slope * t

    # t critical value
    t_crit = stats.t.ppf(1 - alpha / 2, df=EDOF)

    # SE of the mean response using slope_se
    se_mean = slope_se * np.sqrt(ss_t) * np.sqrt(1 / n + (t - t_mean) ** 2 / ss_t)

    ci_lower = x_fit - t_crit * se_mean
    ci_upper = x_fit + t_crit * se_mean
    return ci_lower, ci_upper

def confband_not_timeseries(
    x: ArrayLike, y: ArrayLike, alpha: float = 0.05
) -> tuple[np.ndarray, np.ndarray]:
    trend_result = stats.linregress(x,y)
    slope = trend_result.slope
    intercept = trend_result.intercept
    EDOF = np.min(
        [effective_dof(x,dt = 0.5), effective_dof(y,dt = 0.5),]
        )
    print(EDOF)
    N = x.size
    slope_se = trend_result.stderr  * np.sqrt(N / EDOF)  # standard error of the slope
    x_mean = x.mean()
    ss_x = ((x - x_mean) ** 2).sum()

    y_fit = intercept + slope * x

    # t critical value
    t_crit = stats.t.ppf(1 - alpha / 2, df=EDOF)

    # SE of the mean response using slope_se
    se_mean = slope_se * np.sqrt(ss_x) * np.sqrt(1 / N + (x - x_mean) ** 2 / ss_x)

    ci_lower = y_fit - t_crit * se_mean
    ci_upper = y_fit + t_crit * se_mean
    return ci_lower, ci_upper

