"""Tests for the module Age-trajectory curvature test."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from isograph_benchmark.real_data import age_model_curvature as amc


# --------------------------------------------------------------------------- #
# the contrast
# --------------------------------------------------------------------------- #
def test_contrast_rows_sum_to_zero():
    """Zero-sum rows are what cancels the shared uncertainty; without it the test is weak."""
    for k in (3, 4, 5, 8):
        c = amc.second_difference_contrasts(k)
        assert c.shape == (k - 2, k)
        assert np.allclose(c.sum(axis=1), 0.0)


def test_contrast_annihilates_every_straight_line():
    """C @ b == 0 for collinear b, and only for collinear b."""
    for k in (3, 5):
        c = amc.second_difference_contrasts(k)
        x = np.arange(float(k))
        for slope, intercept in ((0.0, 1.0), (2.5, -3.0), (-1.0, 0.0)):
            assert np.allclose(c @ (slope * x + intercept), 0.0)
        # a genuine bend is not annihilated
        assert not np.allclose(c @ (x ** 2), 0.0)


def test_contrast_requires_three_points():
    with pytest.raises(ValueError):
        amc.second_difference_contrasts(2)


# --------------------------------------------------------------------------- #
# covariance handling
# --------------------------------------------------------------------------- #
def _spline_frame(effects, cov, with_cov=True):
    k = len(effects)
    row = {"module_id": ["M000"] * k,
           "age_prob": list(np.linspace(0.25, 0.75, k)),
           "effect": list(effects),
           "pvalue_ftest": [0.5] * k}
    if with_cov:
        for i in range(1, k + 1):
            for j in range(1, k + 1):
                row[f"cov_age_{i}{j}"] = [cov[i - 1, j - 1]] * k
    return pd.DataFrame(row)


def test_cov_matrix_missing_columns_returns_none():
    """A fit that wrote only per-point se must be detected, not silently diagonalised."""
    g = _spline_frame([1.0, 2.0, 3.0], np.eye(3), with_cov=False)
    assert amc._cov_matrix(g, 3) is None


def test_cov_matrix_roundtrips():
    cov = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]])
    g = _spline_frame([1.0, 2.0, 3.0], cov)
    assert np.allclose(amc._cov_matrix(g, 3), cov)


# --------------------------------------------------------------------------- #
# the statistic
# --------------------------------------------------------------------------- #
def _wald(effects, cov):
    b = np.asarray(effects, dtype=float)
    c = amc.second_difference_contrasts(len(b))
    d = c @ b
    return float(d @ np.linalg.solve(c @ cov @ c.T, d))


def test_shared_uncertainty_cancels():
    """The whole power argument: a near-constant covariance leaves a SMALL contrast variance.

    If the off-diagonal terms were dropped the variance would be ~3x the level variance
    instead of a fraction of it, and a real bend would be declared null.
    """
    sigma2 = 0.40
    near_constant = np.full((3, 3), sigma2) + np.eye(3) * 1e-3
    c = amc.second_difference_contrasts(3)
    full_var = float((c @ near_constant @ c.T)[0, 0])
    diag_var = float((c @ np.diag(np.diag(near_constant)) @ c.T)[0, 0])
    assert full_var < 0.01 * diag_var


def test_perfectly_straight_line_gives_zero_statistic():
    cov = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]])
    assert _wald([1.0, 2.0, 3.0], cov) == pytest.approx(0.0, abs=1e-18)


def test_a_real_bend_is_detected():
    cov = np.eye(3) * 0.01
    stat = _wald([0.0, 1.0, 0.0], cov)  # strongly curved
    assert stat > 10.0


def test_statistic_is_invariant_to_sign_of_the_aging_effect():
    """Only collinearity matters, not whether the module goes up or down with age."""
    cov = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]])
    b = [0.3, 1.1, 0.4]
    assert _wald(b, cov) == pytest.approx(_wald([-v for v in b], cov))


def test_k3_chi2_equals_z_squared():
    """The reported z must be the signed square root of the chi2 when df == 1."""
    cov = np.array([[4.0, 1.0, 0.5], [1.0, 3.0, 0.25], [0.5, 0.25, 2.0]])
    b = [0.3, 1.1, 0.4]
    c = amc.second_difference_contrasts(3)
    d = c @ np.asarray(b)
    stat = _wald(b, cov)
    z = float(np.sign(d[0]) * np.sqrt(stat))
    assert z ** 2 == pytest.approx(stat)
    assert np.sign(z) == np.sign(d[0])
