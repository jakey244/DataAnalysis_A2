"""Answer key for the correlation-and-trends assignment (Lecture 4).

Worked helpers (``autocorr``, ``integral_timescale``, ``effective_dof``, and the
``data_io`` loader) are provided so students focus on the two estimators left as
stubs in the student version: ``cross_correlation`` and ``trend_with_significance``.
The accompanying ``pytest`` checks encode the behaviour they must satisfy.
"""

from .correlation import (
    autocorr,
    integral_timescale,
    effective_dof,
    cross_correlation,
)
from .trends import fit_trend, trend_with_significance, TrendResult
from .data_io import load_amoc, load_ts_gridded, load_47n, load_moc_sigma0_26n
from .geostrophy import (
    to_teos10,
    dynamic_height,
    interior_geostrophic_transport,
)
from .seasonal import seasonal_climatology, remove_seasonal_cycle

__all__ = [
    "autocorr",
    "integral_timescale",
    "effective_dof",
    "cross_correlation",
    "fit_trend",
    "trend_with_significance",
    "TrendResult",
    "load_amoc",
    "load_ts_gridded",
    "load_47n",
    "load_moc_sigma0_26n",
    "to_teos10",
    "dynamic_height",
    "interior_geostrophic_transport",
    "seasonal_climatology",
    "remove_seasonal_cycle",
]
