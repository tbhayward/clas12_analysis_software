#!/usr/bin/env python3
"""Small algebra and recovery checks for the no-flip model."""

from __future__ import annotations

import numpy as np
import pandas as pd

from fit_no_nucleon_flip_amplitudes import (
    amplitude_bilinears,
    common_observable_y_limits,
    fit_one_bin,
    fit_second_sector_at_fraction,
    predict_ratios,
)
from helicity_bilinear_transforms import SOURCE_PARAMETERS, transform_values


def test_forward_change_of_basis() -> None:
    amplitudes = np.asarray([1.2 + 0.3j, -0.4 + 0.8j, 0.7 - 0.2j])
    epsilon = 0.73
    direct = amplitude_bilinears(amplitudes, epsilon)
    through_ratios = transform_values(predict_ratios(amplitudes, epsilon))
    np.testing.assert_allclose(through_ratios, direct, rtol=1e-12, atol=1e-12)


def test_scale_and_phase_invariance() -> None:
    amplitudes = np.asarray([0.8 - 0.2j, 0.3 + 1.1j, -0.5 + 0.4j])
    epsilon = 0.61
    baseline = predict_ratios(amplitudes, epsilon)
    changed = predict_ratios(4.7 * np.exp(1.3j) * amplitudes, epsilon)
    np.testing.assert_allclose(changed, baseline, rtol=1e-12, atol=1e-12)


def test_exact_recovery() -> None:
    amplitudes = np.asarray([1.0 + 0.0j, 0.45 + 0.35j, -0.30 + 0.60j])
    epsilon = 0.68
    measured = predict_ratios(amplitudes, epsilon)
    covariance = np.diag(np.asarray([0.025, 0.025, 0.030, 0.030, 0.035, 0.025, 0.030]) ** 2)
    fit = fit_one_bin(measured, covariance, epsilon, starts=20, rng=np.random.default_rng(7))
    assert fit.success
    assert fit.chi2 < 1.0e-16
    np.testing.assert_allclose(fit.prediction, measured, rtol=1e-9, atol=1e-9)
    assert fit.dof == 3


def test_common_observable_y_limits() -> None:
    rows = []
    for x_index in range(2):
        for observable_index, observable in enumerate(SOURCE_PARAMETERS):
            rows.append({
                "x_index": x_index,
                "observable": observable,
                "measured": -0.2 + 0.1 * x_index + 0.01 * observable_index,
                "fit_uncertainty": 0.03,
                "predicted": 0.15 - 0.01 * observable_index,
            })
    detail = pd.DataFrame(rows)
    limits = common_observable_y_limits(detail)
    assert set(limits) == set(SOURCE_PARAMETERS)
    for observable, (lower, upper) in limits.items():
        points = detail.loc[detail["observable"] == observable]
        assert lower <= 0.0 <= upper
        assert lower <= np.min(points["measured"] - points["fit_uncertainty"])
        assert upper >= np.max(points["measured"] + points["fit_uncertainty"])
        assert lower <= np.min(points["predicted"])
        assert upper >= np.max(points["predicted"])


def test_fixed_second_sector_fraction_recovery() -> None:
    epsilon = 0.66
    sector_a = predict_ratios(np.asarray([1.0, 0.35 + 0.20j, -0.25 + 0.55j]), epsilon)
    sector_b = predict_ratios(np.asarray([0.7 - 0.1j, -0.5 + 0.8j, 0.4 + 0.2j]), epsilon)
    fraction = 0.32
    measured = (1.0 - fraction) * sector_a + fraction * sector_b
    covariance = np.diag(np.full(7, 0.025) ** 2)
    fit = fit_second_sector_at_fraction(
        measured,
        covariance,
        epsilon,
        sector_a,
        fraction,
        starts=20,
        rng=np.random.default_rng(19),
    )
    assert fit.success
    assert fit.chi2 < 1.0e-14
    np.testing.assert_allclose(fit.combined_prediction, measured, rtol=1e-8, atol=1e-8)


def test_zero_second_sector_is_anchor() -> None:
    epsilon = 0.71
    anchor = predict_ratios(np.asarray([1.0, 0.4 + 0.1j, -0.2 + 0.6j]), epsilon)
    measured = anchor + np.asarray([0.01, 0.0, -0.01, 0.02, 0.0, -0.02, 0.01])
    covariance = np.diag(np.full(7, 0.03) ** 2)
    fit = fit_second_sector_at_fraction(
        measured, covariance, epsilon, anchor, 0.0, starts=2, rng=np.random.default_rng(3)
    )
    np.testing.assert_allclose(fit.combined_prediction, anchor)
    np.testing.assert_allclose(fit.residual, measured - anchor)


def main() -> None:
    test_forward_change_of_basis()
    test_scale_and_phase_invariance()
    test_exact_recovery()
    test_common_observable_y_limits()
    test_fixed_second_sector_fraction_recovery()
    test_zero_second_sector_is_anchor()
    print("All no-flip amplitude tests passed.")


if __name__ == "__main__":
    main()
