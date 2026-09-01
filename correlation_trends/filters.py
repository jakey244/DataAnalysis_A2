"""Filter frequency-response helpers.

``nyquist_frequency`` and ``tukey_lowpass`` are worked helpers;
``butterworth_squared_response`` returns the filter's squared frequency response
(the closed form is given in its docstring).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def tukey_lowpass(
    values: np.ndarray,
    window: int,
    center: bool = True,
) -> np.ndarray:
    """Low-pass filter via a Tukey-windowed rolling mean (pandas).

    Mirrors the lab pattern
    ``series.rolling(window, center=True, win_type="tukey", min_periods=1).mean()``
    — a tapered-cosine running mean that looks like a boxcar but is "less jaggedy"
    and has cleaner spectral (frequency-response) properties.

    Parameters
    ----------
    values : numpy.ndarray, shape (N,)
        Evenly sampled series.
    window : int
        Window length in samples (e.g. 20 for a ~10-day filter on a 12-hour grid).
    center : bool, optional
        Centre the window (zero phase). Default ``True``.

    Returns
    -------
    filtered : numpy.ndarray, shape (N,)
        The Tukey-smoothed series.

    Notes
    -----
    Uses scipy's default Tukey shape parameter (``alpha = 0.5``). Compare its
    transfer function with a boxcar of the same length to see the reduced sidelobes.
    """
    s = pd.Series(np.asarray(values, dtype="float64"))
    out = s.rolling(window, center=center, win_type="tukey", min_periods=1).mean()
    return out.to_numpy()


def nyquist_frequency(dt_days: float) -> float:
    """Nyquist frequency for a regular grid.

    Parameters
    ----------
    dt_days : float
        Sample spacing in days.

    Returns
    -------
    f_nyquist : float
        Nyquist frequency in cycles per day, ``1 / (2 * dt_days)``.

    Examples
    --------
    >>> nyquist_frequency(0.5)
    1.0
    """
    return 1.0 / (2.0 * dt_days)


def butterworth_squared_response(
    freq: np.ndarray,
    f_cut: float,
    order: int = 5,
    zero_phase: bool = True,
) -> np.ndarray:
    """Squared magnitude response ``|H(f)|**2`` of a low-pass Butterworth filter.

    Parameters
    ----------
    freq : numpy.ndarray, shape (M,)
        Frequencies (cycles per day) at which to evaluate the response.
    f_cut : float
        Half-power cutoff frequency (cycles per day) — where ``|H|**2`` equals 0.5.
    order : int, optional
        Filter order ``n``. Default 5.
    zero_phase : bool, optional
        If ``True``, model a zero-phase (forward-backward, ``filtfilt``)
        application, which squares the magnitude response and so doubles the
        effective order. Default ``True``.

    Returns
    -------
    h2 : numpy.ndarray, shape (M,)
        Squared magnitude response, in ``[0, 1]``.

    Notes
    -----
    The single-pass Butterworth squared response is

    .. math::

        |H(f)|^2 = \\frac{1}{1 + (f / f_\\mathrm{cut})^{2n}}.

    A zero-phase application multiplies the magnitude by itself, i.e. squares the
    expression above (equivalently doubles the exponent ``2n -> 4n``). At
    ``f = f_cut`` the single-pass response is exactly 0.5 — useful as a test.

    TODO (student): implement and return ``h2``.
    """
    raise NotImplementedError("Implement |H(f)|**2 from the formula in the docstring.")
