#!/usr/bin/env python3
"""Small algebra and recovery checks for the no-flip model."""

from __future__ import annotations

import numpy as np

from fit_no_nucleon_flip_amplitudes import (
    amplitude_bilinears,
    fit_one_bin,
    predict_ratios,
)
from helicity_bilinear_transforms import transform_values


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


def main() -> None:
    test_forward_change_of_basis()
    test_scale_and_phase_invariance()
    test_exact_recovery()
    print("All no-flip amplitude tests passed.")


if __name__ == "__main__":
    main()
