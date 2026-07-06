"""Autocorrelation, decorrelation timescale, and cross-correlation.

Worked helpers: :func:`autocorr`, :func:`integral_timescale`, :func:`effective_dof`.
Student stub: :func:`cross_correlation`.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import signal


def autocorr(x: ArrayLike, biased: bool = True) -> NDArray[np.float64]:
    """One-sided sample autocorrelation ``R(tau)`` with ``R(0) = 1``.

    Parameters
    ----------
    x : array_like
        Input series. Non-finite values are dropped before estimation.
    biased : bool, default True
        If True, normalise each lag by ``N`` (the biased estimator: a triangular
        ``(N - tau)/N`` taper, low variance, positive-definite). If False,
        normalise by ``N - tau`` (unbiased in the mean, but the variance grows
        without bound as ``tau -> N``).

    Returns
    -------
    numpy.ndarray
        ``R(tau)`` for lags ``0, 1, ..., N-1``.

    Notes
    -----
    The biased and unbiased estimates agree at small lag (many pairs) and diverge
    at large lag; the biased tail is damped toward zero while the unbiased tail
    fans out. See the lecture figure comparing white noise, AR(1), and a sinusoid.
    """
    d = np.asarray(x, dtype=float)
    d = d[np.isfinite(d)]
    d = d - d.mean()
    n = d.size
    var = d.var()
    full = np.correlate(d, d, "full")[n - 1 :]  # lags 0 .. n-1
    if biased:
        return full / (n * var)
    return full / ((n - np.arange(n)) * var)


def integral_timescale(x: ArrayLike, dt: float, biased: bool = True) -> float:
    """Integral (decorrelation) timescale ``T*`` by trapezoid to the first zero.

    ``T* = int_0^inf R(tau) dtau``, estimated as the trapezoidal sum of the
    autocorrelation truncated at its first zero crossing. Truncating there keeps
    essentially all the real area while discarding the noisy large-lag tail; for
    the biased estimator this is close to the full integral (the tail is damped).

    Parameters
    ----------
    x : array_like
        Input series.
    dt : float
        Sample spacing, in the time units you want ``T*`` expressed in.
    biased : bool, default True
        Estimator passed to :func:`autocorr`. Biased is recommended.

    Returns
    -------
    float
        The one-sided integral timescale ``T*`` (same units as ``dt``).
    """
    R = autocorr(x, biased=biased)
    tau = 0.0
    i = 0
    while i + 1 < len(R) and R[i] >= 0:
        tau += dt * (R[i] + R[i + 1]) / 2.0
        i += 1
    return tau


def effective_dof(x: ArrayLike, dt: float, biased: bool = True) -> float:
    """Effective (equivalent) degrees of freedom ``EDOF`` for a series.

    Following Emery & Thomson, the raw degrees of freedom are
    ``DOF = record / T*`` and the *equivalent* degrees of freedom are
    ``EDOF = DOF / 2 = record / (2 T*)``. Use ``EDOF`` — not ``DOF`` — for error
    bars, confidence intervals, and significance tests: the factor of two comes
    from the two-sided variance sum ``1 + 2 sum rho_k``.

    Parameters
    ----------
    x : array_like
        Input series.
    dt : float
        Sample spacing.
    biased : bool, default True
        Estimator passed through to :func:`integral_timescale`.

    Returns
    -------
    float
        ``EDOF = record / (2 T*)``, where ``record = N * dt``.
    """
    # tstar = integral_timescale(x, dt, biased=biased)
    # record = N * dt  (N = number of finite samples)
    # EDOF = record / (2 * T*)   -- the factor of 2 is DOF -> EDOF
    raise NotImplementedError("effective_dof")


def cross_correlation(
    x: ArrayLike, y: ArrayLike
) -> tuple[NDArray[np.int_], NDArray[np.float64]]:
    """Normalised cross-correlation of two series, with the lag axis.

    Both series are standardised (subtract mean, divide by standard deviation)
    and correlated at every integer lag.

    Parameters
    ----------
    x, y : array_like
        Two series of equal length.

    Returns
    -------
    lags : numpy.ndarray
        Integer lags, from ``-(N-1)`` to ``N-1``.
    r : numpy.ndarray
        Cross-correlation at each lag; ``r`` peaks at the lag that best aligns
        the two series.

    Notes
    -----
    Sign convention: the peak sits at the lag by which ``y`` is shifted relative
    to ``x``. If ``y`` lags ``x`` by ``k`` samples (``x`` leads) the peak is at
    ``-k``; so a **negative** peak lag means ``x`` leads ``y``, a positive one
    means ``y`` leads ``x``.
    """
    # xa = (x - np.mean(x)) / np.std(x)
    # ya = (y - np.mean(y)) / np.std(y)
    # then signal.correlate(xa, ya, mode="full") / N and signal.correlation_lags(...)
    raise NotImplementedError("cross_correlation")
