""" TEST FOR THE DETRENDING FUNCTION"""
from __future__ import annotations

import numpy as np
import pytest

from amocatlas import read

from correlation_trends.analysis import detrending


def test_detrending_detrends_a_known_series() -> None:
    t = read.rapid().MOC["TIME"][:10000]
    n = t.size
    x = np.random.rand(n)*10 +5 
    x_res = x - x.mean()
    x_slope = x + np.arange(n)
    x_detrended = detrending(x_slope,t)
    assert np.allclose(x_detrended.values,x_res, rtol = 0.1, atol = 0.1)