from __future__ import annotations

import numpy as np
import pytest

from amocatlas import read

from correlation_trends.analysis import detrending


def detrending_detrends_a_known_series() -> None:
    t = read.rapid().MOC["TIME"][:100]
    n = t.size
    x = np.random.rand(n) + np.arange(n) * 0.0005 + 10
    x_detrended = detrending(x,t)
    assert x_detrended == pytest.approx(x, abs=0.02)