==============================
Assignment 2 — starter package
==============================

*(This folder is generated from ``../solution2/`` — do not edit it by hand;
edit the solution and re-run ``generate_assignment.py``.)*

Implement the stubbed functions in ``correlation_trends/`` until the tests pass.

Layout
------

- ``correlation_trends/`` — the package. Worked helpers (the data loaders, the
  autocorrelation, the integral timescale, the seasonal climatology, and the
  geostrophy routines) are provided. You implement the functions that
  ``raise NotImplementedError``: ``correlation.effective_dof``,
  ``correlation.cross_correlation``, ``seasonal.remove_seasonal_cycle``, and
  ``trends.trend_with_significance``.
- ``tests/`` — ``pytest`` checks that pin each function's contract.

Run::

    pip install -r requirements.txt
    ruff format . && ruff check .    # format, then lint
    pytest -q                        # implement until the tests are green

Each stub's docstring states exactly what to return.
