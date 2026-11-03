#!/usr/bin/env python3
"""Test a parity-related three-amplitude, no-nucleon-flip model.

This is a deliberately small second-stage fit to the seven ratios already
extracted in each kinematic bin.  It does not repeat the event-level fit.

The retained independent amplitudes are

    L  = T_{0,+}^{0,+},   P = T_{+1,+}^{0,+},   M = T_{-1,+}^{0,+}.

All amplitudes with recoil-neutron helicity opposite to the initial-proton
helicity are set to zero.  Parity supplies the corresponding initial-proton
helicity - amplitudes.  A common complex phase and a common scale cancel in
the seven ratios, leaving four real fit parameters per bin.

After the original fit, the script can hold that best one-sector prediction
fixed and profile the denominator fraction carried by a second incoherent
recoil-helicity sector.  This is a conditional diagnostic of how much added
sector strength is needed to repair the one-sector discrepancies.  It is not
an unconstrained two-sector amplitude extraction.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.optimize import least_squares
from scipy.stats import chi2

from helicity_bilinear_transforms import (
    HELICITY_LABELS,
    HELICITY_PARAMETERS,
    SOURCE_PARAMETERS,
    TRANSFORM,
)


SOURCE_LABELS = {
    "u1": r"$F_{UU}^{\cos\phi}/F_{UU}$",
    "u2": r"$F_{UU}^{\cos2\phi}/F_{UU}$",
    "lu1": r"$F_{LU}^{\sin\phi}/F_{UU}$",
    "ul1": r"$F_{UL}^{\sin\phi}/F_{UU}$",
    "ul2": r"$F_{UL}^{\sin2\phi}/F_{UU}$",
    "ll0": r"$F_{LL}/F_{UU}$",
    "ll1": r"$F_{LL}^{\cos\phi}/F_{UU}$",
}

XB_LABELS = (
    r"$0.10<x_B<0.25$",
    r"$0.25<x_B<0.35$",
    r"$0.35<x_B<0.45$",
    r"$0.45<x_B<0.60$",
)
XB_STYLES = (
    ("#0072B2", "o"),
    ("#D55E00", "s"),
    ("#009E73", "^"),
    ("#CC79A7", "D"),
)
AMPLITUDE_NAMES = ("L", "P", "M")


@dataclass
class BinFit:
    amplitudes: np.ndarray
    prediction: np.ndarray
    residual: np.ndarray
    marginal_pull: np.ndarray
    chi2: float
    rank: int
    dof: int
    p_value: float
    reference: str
    success: bool
    message: str


@dataclass
class SecondSectorFit:
    """Best added-sector solution at one fixed denominator fraction."""

    fraction: float
    amplitudes: np.ndarray
    sector_prediction: np.ndarray
    combined_prediction: np.ndarray
    residual: np.ndarray
    marginal_pull: np.ndarray
    chi2: float
    rank: int
    dof: int
    p_value: float
    reference: str
    success: bool
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit the restricted three-complex-amplitude model to all seven RGC ratios."
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=Path("output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv"),
    )
    parser.add_argument(
        "--input-json",
        type=Path,
        default=Path("output/asymmetry_extraction/nominal/json/structure_function_ratios.json"),
        help="Production fit JSON containing the full seven-ratio statistical covariance.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/no_nucleon_flip_amplitude_test"),
    )
    parser.add_argument(
        "--default-epsilon",
        type=float,
        default=None,
        help="Fallback epsilon only when mean_epsilon is absent from the CSV.",
    )
    parser.add_argument(
        "--include-diagonal-systematics",
        action="store_true",
        help="Add marginal point-to-point systematic magnitudes in quadrature. This is approximate.",
    )
    parser.add_argument("--starts", type=int, default=40)
    parser.add_argument(
        "--skip-second-sector-profile",
        action="store_true",
        help="Run only the original one-sector study and preserve its historical behavior.",
    )
    parser.add_argument(
        "--profile-points",
        type=int,
        default=51,
        help="Number of equally spaced second-sector fractions between zero and the maximum.",
    )
    parser.add_argument(
        "--profile-max-fraction",
        type=float,
        default=1.0,
        help="Largest second-sector denominator fraction to test.",
    )
    parser.add_argument(
        "--profile-starts",
        type=int,
        default=8,
        help="Multistart count for the added-sector fit at each fixed fraction.",
    )
    parser.add_argument(
        "--profile-acceptance-alpha",
        type=float,
        default=0.05,
        help="Fixed-fraction goodness-of-fit rejection level used to quote the minimum fraction.",
    )
    parser.add_argument("--seed", type=int, default=24601)
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def amplitude_bilinears(amplitudes: np.ndarray, epsilon: float) -> np.ndarray:
    """Return the seven normalized bilinears in HELICITY_PARAMETERS order."""
    L, P, M = np.asarray(amplitudes, dtype=complex)
    denominator = 0.5 * (abs(P) ** 2 + abs(M) ** 2) + epsilon * abs(L) ** 2
    if not np.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("Amplitude denominator must be positive.")
    return np.asarray(
        [
            np.imag(L * np.conjugate(P)),
            -np.imag(L * np.conjugate(M)),
            np.imag(M * np.conjugate(P)),
            abs(P) ** 2 - abs(M) ** 2,
            np.real(L * np.conjugate(P)),
            -np.real(L * np.conjugate(M)),
            np.real(M * np.conjugate(P)),
        ],
        dtype=float,
    ) / denominator


def predict_ratios(amplitudes: np.ndarray, epsilon: float) -> np.ndarray:
    """Return (u1,u2,lu1,ul1,ul2,ll0,ll1) predicted by the model."""
    return np.linalg.solve(TRANSFORM, amplitude_bilinears(amplitudes, epsilon))


def amplitudes_from_chart(parameters: np.ndarray, reference: int) -> np.ndarray:
    """Use one real unit amplitude to remove the unobservable scale and phase."""
    parameters = np.asarray(parameters, dtype=float)
    amplitudes = np.empty(3, dtype=complex)
    amplitudes[reference] = 1.0 + 0.0j
    other = [index for index in range(3) if index != reference]
    amplitudes[other[0]] = parameters[0] + 1j * parameters[1]
    amplitudes[other[1]] = parameters[2] + 1j * parameters[3]
    return amplitudes


def canonicalize(amplitudes: np.ndarray) -> np.ndarray:
    """Normalize sum |T|^2 to one and make the largest amplitude real positive."""
    output = np.asarray(amplitudes, dtype=complex).copy()
    norm = np.sqrt(np.sum(np.abs(output) ** 2))
    if norm == 0.0:
        return output
    output /= norm
    anchor = int(np.argmax(np.abs(output)))
    output *= np.exp(-1j * np.angle(output[anchor]))
    return output


def covariance_whitener(covariance: np.ndarray) -> tuple[np.ndarray, int]:
    covariance = 0.5 * (np.asarray(covariance, dtype=float) + np.asarray(covariance, dtype=float).T)
    values, vectors = np.linalg.eigh(covariance)
    tolerance = max(float(np.max(values)), 1.0) * 1.0e-12
    keep = values > tolerance
    if not np.any(keep):
        raise ValueError("Covariance has no positive eigenvalues.")
    whitener = (vectors[:, keep] / np.sqrt(values[keep])).T
    return whitener, int(np.count_nonzero(keep))


def fit_one_bin(
    measured: np.ndarray,
    covariance: np.ndarray,
    epsilon: float,
    starts: int,
    rng: np.random.Generator,
) -> BinFit:
    whitener, rank = covariance_whitener(covariance)
    best = None
    best_reference = 0
    for reference in range(3):
        initial_points = [np.zeros(4)]
        initial_points.extend(rng.normal(0.0, 1.25, size=(max(starts - 1, 0), 4)))
        for initial in initial_points:
            result = least_squares(
                lambda pars: whitener @ (
                    measured - predict_ratios(amplitudes_from_chart(pars, reference), epsilon)
                ),
                initial,
                method="trf",
                max_nfev=5000,
                xtol=1.0e-12,
                ftol=1.0e-12,
                gtol=1.0e-12,
            )
            objective = float(np.dot(result.fun, result.fun))
            if best is None or objective < best[0]:
                best = (objective, result)
                best_reference = reference
    assert best is not None
    objective, result = best
    amplitudes = amplitudes_from_chart(result.x, best_reference)
    prediction = predict_ratios(amplitudes, epsilon)
    residual = measured - prediction
    marginal_sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    marginal_pull = np.divide(
        residual,
        marginal_sigma,
        out=np.full_like(residual, np.nan),
        where=marginal_sigma > 0.0,
    )
    dof = max(rank - 4, 0)
    return BinFit(
        amplitudes=canonicalize(amplitudes),
        prediction=prediction,
        residual=residual,
        marginal_pull=marginal_pull,
        chi2=objective,
        rank=rank,
        dof=dof,
        p_value=float(chi2.sf(objective, dof)) if dof > 0 else np.nan,
        reference=AMPLITUDE_NAMES[best_reference],
        success=bool(result.success),
        message=str(result.message),
    )


def fit_second_sector_at_fraction(
    measured: np.ndarray,
    covariance: np.ndarray,
    epsilon: float,
    anchor_prediction: np.ndarray,
    fraction: float,
    starts: int,
    rng: np.random.Generator,
) -> SecondSectorFit:
    """Fit the shape of a second incoherent sector at fixed strength fraction.

    ``fraction`` is D_B/(D_A+D_B), where D is the contribution to the
    unpolarized denominator.  The first-sector prediction is held fixed.
    For fraction > 0 the problem can be mapped exactly onto ``fit_one_bin``
    with a shifted pseudo-measurement and a rescaled covariance.
    """
    measured = np.asarray(measured, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    anchor_prediction = np.asarray(anchor_prediction, dtype=float)
    fraction = float(fraction)
    if fraction < 0.0 or fraction > 1.0:
        raise ValueError("Second-sector fraction must lie between zero and one.")
    whitener, rank = covariance_whitener(covariance)
    dof = max(rank - 4, 0)
    if fraction == 0.0:
        combined = anchor_prediction.copy()
        residual = measured - combined
        objective = float(np.dot(whitener @ residual, whitener @ residual))
        sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        pull = np.divide(residual, sigma, out=np.full_like(residual, np.nan), where=sigma > 0.0)
        return SecondSectorFit(
            fraction=0.0,
            amplitudes=np.full(3, np.nan + 1j * np.nan),
            sector_prediction=np.full_like(measured, np.nan),
            combined_prediction=combined,
            residual=residual,
            marginal_pull=pull,
            chi2=objective,
            rank=rank,
            dof=dof,
            p_value=float(chi2.sf(objective, dof)) if dof > 0 else np.nan,
            reference="none",
            success=True,
            message="No added sector at fraction zero.",
        )
    pseudo_measured = (measured - (1.0 - fraction) * anchor_prediction) / fraction
    pseudo_covariance = covariance / fraction**2
    sector_fit = fit_one_bin(pseudo_measured, pseudo_covariance, epsilon, starts, rng)
    combined = (1.0 - fraction) * anchor_prediction + fraction * sector_fit.prediction
    residual = measured - combined
    objective = float(np.dot(whitener @ residual, whitener @ residual))
    sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    pull = np.divide(residual, sigma, out=np.full_like(residual, np.nan), where=sigma > 0.0)
    return SecondSectorFit(
        fraction=fraction,
        amplitudes=sector_fit.amplitudes,
        sector_prediction=sector_fit.prediction,
        combined_prediction=combined,
        residual=residual,
        marginal_pull=pull,
        chi2=objective,
        rank=rank,
        dof=dof,
        p_value=float(chi2.sf(objective, dof)) if dof > 0 else np.nan,
        reference=sector_fit.reference,
        success=sector_fit.success,
        message=sector_fit.message,
    )


def source_covariances(path: Path) -> dict[int, np.ndarray]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    output: dict[int, np.ndarray] = {}
    for result in payload["results"]:
        fit = result["variants"]["nominal"]
        order = list(fit["parameter_order"])
        full = np.asarray(fit["covariance"], dtype=float)
        indices = [order.index(name) for name in SOURCE_PARAMETERS]
        output[int(result["bin_number"])] = full[np.ix_(indices, indices)]
    return output


def diagonal_stat_covariance(row: pd.Series) -> np.ndarray:
    errors = np.asarray([float(row[f"{name}_stat"]) for name in SOURCE_PARAMETERS])
    return np.diag(errors**2)


def point_to_point_variances(row: pd.Series) -> np.ndarray:
    variances = []
    for name in SOURCE_PARAMETERS:
        candidates = (
            f"{name}_point_to_point_systematic",
            f"{name}_point_to_point_systematic_diagonal",
        )
        value = next((float(row[key]) for key in candidates if key in row and pd.notna(row[key])), 0.0)
        variances.append(value**2)
    return np.asarray(variances)


def epsilon_for_row(row: pd.Series, fallback: float | None) -> float:
    for column in ("mean_epsilon", "epsilon"):
        if column in row and pd.notna(row[column]):
            return float(row[column])
    if fallback is None:
        raise KeyError("CSV has no mean_epsilon column. Supply --default-epsilon explicitly.")
    return float(fallback)


def fit_frame(
    frame: pd.DataFrame,
    covariances: dict[int, np.ndarray],
    default_epsilon: float | None,
    include_systematics: bool,
    starts: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    summary_rows: list[dict] = []
    detail_rows: list[dict] = []
    for _, row in frame.sort_values("bin_number").iterrows():
        bin_number = int(row["bin_number"])
        measured = np.asarray([float(row[name]) for name in SOURCE_PARAMETERS])
        covariance = covariances.get(bin_number, diagonal_stat_covariance(row)).copy()
        if include_systematics:
            covariance += np.diag(point_to_point_variances(row))
        epsilon = epsilon_for_row(row, default_epsilon)
        fit = fit_one_bin(measured, covariance, epsilon, starts, rng)
        amplitudes = fit.amplitudes
        common = {
            "bin_number": bin_number,
            "x_index": int(row["x_index"]),
            "t_index": int(row["t_index"]),
            "mean_minus_tprime_gev2": float(row["mean_minus_tprime_gev2"]),
            "mean_epsilon_used": epsilon,
        }
        summary = {
            **common,
            "chi2": fit.chi2,
            "covariance_rank": fit.rank,
            "model_parameter_count": 4,
            "dof": fit.dof,
            "chi2_per_dof": fit.chi2 / fit.dof if fit.dof else np.nan,
            "p_value": fit.p_value,
            "fit_reference_chart": fit.reference,
            "fit_success": fit.success,
            "fit_message": fit.message,
        }
        for name, amplitude in zip(AMPLITUDE_NAMES, amplitudes):
            summary[f"{name}_real_normalized"] = float(amplitude.real)
            summary[f"{name}_imag_normalized"] = float(amplitude.imag)
            summary[f"{name}_magnitude_normalized"] = float(abs(amplitude))
            summary[f"{name}_strength_fraction"] = float(abs(amplitude) ** 2)
            summary[f"{name}_phase_deg"] = float(np.degrees(np.angle(amplitude)))
        summary_rows.append(summary)
        sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        for index, name in enumerate(SOURCE_PARAMETERS):
            detail_rows.append(
                {
                    **common,
                    "observable": name,
                    "measured": measured[index],
                    "fit_uncertainty": sigma[index],
                    "predicted": fit.prediction[index],
                    "residual": fit.residual[index],
                    "marginal_pull": fit.marginal_pull[index],
                }
            )
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def profile_second_sector(
    frame: pd.DataFrame,
    covariances: dict[int, np.ndarray],
    anchor_detail: pd.DataFrame,
    default_epsilon: float | None,
    include_systematics: bool,
    fractions: np.ndarray,
    starts: int,
    acceptance_alpha: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Profile the strength of a second recoil-helicity sector in every bin."""
    rng = np.random.default_rng(seed)
    summary_rows: list[dict] = []
    scan_rows: list[dict] = []
    diagnostic_rows: list[dict] = []
    fractions = np.unique(np.asarray(fractions, dtype=float))
    for _, row in frame.sort_values("bin_number").iterrows():
        bin_number = int(row["bin_number"])
        measured = np.asarray([float(row[name]) for name in SOURCE_PARAMETERS])
        covariance = covariances.get(bin_number, diagonal_stat_covariance(row)).copy()
        if include_systematics:
            covariance += np.diag(point_to_point_variances(row))
        epsilon = epsilon_for_row(row, default_epsilon)
        anchor_points = anchor_detail.loc[anchor_detail["bin_number"] == bin_number].set_index("observable")
        anchor_prediction = np.asarray([float(anchor_points.loc[name, "predicted"]) for name in SOURCE_PARAMETERS])
        common = {
            "bin_number": bin_number,
            "x_index": int(row["x_index"]),
            "t_index": int(row["t_index"]),
            "mean_minus_tprime_gev2": float(row["mean_minus_tprime_gev2"]),
            "mean_epsilon_used": epsilon,
        }
        fits: list[SecondSectorFit] = []
        for fraction in fractions:
            fit = fit_second_sector_at_fraction(
                measured, covariance, epsilon, anchor_prediction, float(fraction), starts, rng
            )
            fits.append(fit)
            scan = {
                **common,
                "second_sector_fraction": fit.fraction,
                "chi2": fit.chi2,
                "covariance_rank": fit.rank,
                "fixed_fraction_parameter_count": 4 if fit.fraction > 0.0 else 0,
                "fixed_fraction_dof": fit.dof,
                "chi2_per_fixed_fraction_dof": fit.chi2 / fit.dof if fit.dof else np.nan,
                "fixed_fraction_p_value": fit.p_value,
                "fit_reference_chart": fit.reference,
                "fit_success": fit.success,
                "fit_message": fit.message,
            }
            for name, amplitude in zip(AMPLITUDE_NAMES, fit.amplitudes):
                scan[f"sector_B_{name}_real_normalized"] = float(np.real(amplitude))
                scan[f"sector_B_{name}_imag_normalized"] = float(np.imag(amplitude))
            scan_rows.append(scan)
        best = min(fits, key=lambda item: item.chi2)
        positive_fit = next((item for item in fits if item.fraction > 0.0), fits[0])
        threshold = float(chi2.ppf(1.0 - acceptance_alpha, positive_fit.dof)) if positive_fit.dof else np.nan
        acceptable = [item for item in fits if np.isfinite(threshold) and item.chi2 <= threshold]
        minimum = min(acceptable, key=lambda item: item.fraction) if acceptable else None
        selected = minimum if minimum is not None else best
        original = fits[0]
        summary_rows.append({
            **common,
            "original_one_sector_chi2": original.chi2,
            "fixed_fraction_dof": positive_fit.dof,
            "acceptance_alpha": acceptance_alpha,
            "acceptance_chi2_threshold": threshold,
            "minimum_acceptable_fraction": minimum.fraction if minimum is not None else np.nan,
            "minimum_acceptable_fraction_found": minimum is not None,
            "best_profile_fraction": best.fraction,
            "best_profile_chi2": best.chi2,
            "best_profile_delta_chi2": original.chi2 - best.chi2,
            "selected_fraction": selected.fraction,
            "selected_solution": "minimum_acceptable" if minimum is not None else "best_profile_no_acceptable_grid_point",
            "profile_grid_spacing": float(np.min(np.diff(fractions))) if len(fractions) > 1 else np.nan,
        })
        sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
        for index, name in enumerate(SOURCE_PARAMETERS):
            diagnostic_rows.append({
                **common,
                "observable": name,
                "measured": measured[index],
                "fit_uncertainty": sigma[index],
                "anchor_prediction": anchor_prediction[index],
                "selected_prediction": selected.combined_prediction[index],
                "anchor_residual": measured[index] - anchor_prediction[index],
                "selected_residual": selected.residual[index],
                "selected_marginal_pull": selected.marginal_pull[index],
                "selected_fraction": selected.fraction,
                "selected_solution": "minimum_acceptable" if minimum is not None else "best_profile_no_acceptable_grid_point",
            })
    return pd.DataFrame(summary_rows), pd.DataFrame(scan_rows), pd.DataFrame(diagnostic_rows)


def _style_axis(axis) -> None:
    axis.axhline(0.0, color="0.55", ls="--", lw=1.2, zorder=0)
    axis.set_xlim(0.0, 1.20)
    axis.minorticks_on()
    axis.tick_params(direction="in", top=True, right=True, which="both")


def _nice_tick_step(value: float) -> float:
    """Return a readable axis step near ``value``."""
    if not np.isfinite(value) or value <= 0.0:
        return 0.1
    power = 10.0 ** np.floor(np.log10(value))
    fraction = value / power
    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= candidate:
            return float(candidate * power)
    return float(10.0 * power)


def common_observable_y_limits(detail: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Choose one y-axis range per observable for every x_B canvas.

    The limits include every measured point with its displayed uncertainty,
    every fitted prediction and zero.  Rounding to a readable step keeps the
    canvases directly comparable without clipping an uncertainty bar.
    """
    limits: dict[str, tuple[float, float]] = {}
    for name in SOURCE_PARAMETERS:
        points = detail.loc[detail["observable"] == name]
        measured = points["measured"].to_numpy(float)
        uncertainty = points["fit_uncertainty"].to_numpy(float)
        predicted = points["predicted"].to_numpy(float)
        values = np.concatenate((
            measured - uncertainty,
            measured + uncertainty,
            predicted,
            np.asarray([0.0]),
        ))
        values = values[np.isfinite(values)]
        if values.size == 0:
            limits[name] = (-1.0, 1.0)
            continue
        low = float(np.min(values))
        high = float(np.max(values))
        span = max(high - low, 0.05)
        low -= 0.08 * span
        high += 0.08 * span
        step = _nice_tick_step((high - low) / 6.0)
        lower = float(np.floor(low / step) * step)
        upper = float(np.ceil(high / step) * step)
        if lower == upper:
            lower -= step
            upper += step
        limits[name] = (lower, upper)
    return limits


def plot_measured_vs_model(detail: pd.DataFrame, output_dir: Path) -> list[Path]:
    paths = []
    y_limits = common_observable_y_limits(detail)
    for x_index in sorted(detail["x_index"].unique()):
        subset = detail.loc[detail["x_index"] == x_index]
        color, marker = XB_STYLES[int(x_index)]
        fig, axes = plt.subplots(4, 2, figsize=(16, 20))
        axes = axes.ravel()
        for panel, name in enumerate(SOURCE_PARAMETERS):
            axis = axes[panel]
            points = subset.loc[subset["observable"] == name].sort_values("t_index")
            x = points["mean_minus_tprime_gev2"].to_numpy(float)
            y = points["measured"].to_numpy(float)
            error = points["fit_uncertainty"].to_numpy(float)
            model = points["predicted"].to_numpy(float)
            _style_axis(axis)
            axis.set_ylim(*y_limits[name])
            axis.errorbar(
                x, y, yerr=error, fmt=marker, mfc="white", mec=color, mew=2,
                ecolor=color, capsize=3, label="Measured ratio", zorder=2,
            )
            axis.plot(x, model, color="black", lw=2, marker=".", label="Three-amplitude fit")
            axis.set_ylabel(SOURCE_LABELS[name], fontsize=15, labelpad=12)
            axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$", fontsize=12)
            axis.text(.04, .90, f"({chr(97 + panel)})", transform=axis.transAxes,
                      fontsize=14, fontweight="bold")
        axes[-1].axis("off")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(.94, .05), frameon=False)
        title_index = int(x_index)
        title = XB_LABELS[title_index] if title_index < len(XB_LABELS) else f"x bin {title_index}"
        fig.suptitle(
            r"No-nucleon-flip three-amplitude test: " + title,
            fontsize=20,
        )
        fig.tight_layout(rect=(.03, .04, .98, .96))
        path = output_dir / f"xB_bin_{title_index + 1}_measured_vs_three_amplitude_model.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_goodness(summary: pd.DataFrame, output_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    for x_index, (color, marker) in enumerate(XB_STYLES):
        points = summary.loc[summary["x_index"] == x_index].sort_values("t_index")
        if points.empty:
            continue
        x = points["mean_minus_tprime_gev2"]
        axes[0].plot(x, points["chi2_per_dof"], marker=marker, color=color,
                     mfc="white", mew=2, label=XB_LABELS[x_index])
        axes[1].plot(x, points["p_value"], marker=marker, color=color,
                     mfc="white", mew=2)
    axes[0].axhline(1.0, color="0.4", ls="--")
    axes[0].set_ylabel(r"$\chi^2/\mathrm{dof}$")
    axes[1].axhline(0.05, color="0.4", ls="--")
    axes[1].set_ylabel("fit p-value")
    axes[1].set_ylim(-.02, 1.02)
    for axis in axes:
        axis.set_xlim(0.0, 1.20)
        axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$")
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Goodness of fit of the restricted three-amplitude model", fontsize=16)
    fig.tight_layout()
    path = output_dir / "three_amplitude_goodness_of_fit.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_pull_heatmap(detail: pd.DataFrame, output_dir: Path) -> Path:
    ordered = detail.sort_values(["x_index", "t_index", "observable"])
    bin_order = ordered[["bin_number", "x_index", "t_index"]].drop_duplicates()
    matrix = np.full((len(SOURCE_PARAMETERS), len(bin_order)), np.nan)
    for column, bin_number in enumerate(bin_order["bin_number"]):
        points = detail.loc[detail["bin_number"] == bin_number].set_index("observable")
        for row, name in enumerate(SOURCE_PARAMETERS):
            if name in points.index:
                matrix[row, column] = float(points.loc[name, "marginal_pull"])
    figure, axis = plt.subplots(figsize=(15, 5.5))
    image = axis.imshow(matrix, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    axis.set_yticks(range(len(SOURCE_PARAMETERS)), [SOURCE_LABELS[name] for name in SOURCE_PARAMETERS])
    tick_labels = [f"x{int(x)+1},t{int(t)+1}" for x, t in zip(bin_order["x_index"], bin_order["t_index"])]
    axis.set_xticks(range(len(tick_labels)), tick_labels, rotation=60, ha="right")
    axis.set_title("Marginal residuals: (measurement - model) / marginal uncertainty")
    figure.colorbar(image, ax=axis, label="marginal pull")
    figure.tight_layout()
    path = output_dir / "three_amplitude_marginal_pull_heatmap.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def plot_amplitudes(summary: pd.DataFrame, output_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=(17, 9), sharex=True)
    for x_index, (color, marker) in enumerate(XB_STYLES):
        points = summary.loc[summary["x_index"] == x_index].sort_values("t_index")
        if points.empty:
            continue
        x = points["mean_minus_tprime_gev2"]
        for column, name in enumerate(AMPLITUDE_NAMES):
            axes[0, column].plot(x, points[f"{name}_strength_fraction"], marker=marker,
                                 color=color, mfc="white", mew=1.8, label=XB_LABELS[x_index])
            axes[1, column].plot(x, points[f"{name}_phase_deg"], marker=marker,
                                 color=color, mfc="white", mew=1.8)
    for column, name in enumerate(AMPLITUDE_NAMES):
        axes[0, column].set_title(name + r" strength $|T|^2/\sum|T|^2$")
        axes[0, column].set_ylim(-.05, 1.05)
        axes[1, column].set_title(name + " phase in displayed convention")
        axes[1, column].set_ylim(-190, 190)
        axes[1, column].set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$")
    axes[0, 0].set_ylabel("strength fraction")
    axes[1, 0].set_ylabel("phase (degrees)")
    for axis in axes.ravel():
        axis.set_xlim(0.0, 1.20)
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle("Representative amplitude solution (common phase and scale fixed)", fontsize=16)
    fig.tight_layout()
    path = output_dir / "three_amplitude_representative_solution.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def second_sector_y_limits(detail: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Common per-observable limits for all second-sector repair canvases."""
    renamed = detail.rename(columns={"anchor_prediction": "predicted"}).copy()
    extra = detail.rename(columns={"selected_prediction": "predicted"}).copy()
    return common_observable_y_limits(pd.concat([renamed, extra], ignore_index=True))


def plot_second_sector_repairs(detail: pd.DataFrame, output_dir: Path) -> list[Path]:
    paths: list[Path] = []
    y_limits = second_sector_y_limits(detail)
    for x_index in sorted(detail["x_index"].unique()):
        subset = detail.loc[detail["x_index"] == x_index]
        color, marker = XB_STYLES[int(x_index)]
        fig, axes = plt.subplots(4, 2, figsize=(16, 20))
        axes = axes.ravel()
        for panel, name in enumerate(SOURCE_PARAMETERS):
            axis = axes[panel]
            points = subset.loc[subset["observable"] == name].sort_values("t_index")
            x = points["mean_minus_tprime_gev2"].to_numpy(float)
            _style_axis(axis)
            axis.set_ylim(*y_limits[name])
            axis.errorbar(
                x, points["measured"], yerr=points["fit_uncertainty"], fmt=marker,
                mfc="white", mec=color, mew=2, ecolor=color, capsize=3,
                label="Measured ratio", zorder=3,
            )
            axis.plot(x, points["anchor_prediction"], color="0.25", lw=2, ls="--",
                      marker=".", label="Fixed one-sector fit")
            axis.plot(x, points["selected_prediction"], color="#D55E00", lw=2,
                      marker=".", label="Profiled second-sector solution")
            axis.set_ylabel(SOURCE_LABELS[name], fontsize=15, labelpad=12)
            axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$", fontsize=12)
            axis.text(.04, .90, f"({chr(97 + panel)})", transform=axis.transAxes,
                      fontsize=14, fontweight="bold")
        axes[-1].axis("off")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(.95, .05), frameon=False)
        title_index = int(x_index)
        title = XB_LABELS[title_index] if title_index < len(XB_LABELS) else f"x bin {title_index}"
        fig.suptitle("Fixed first sector plus minimum profiled second sector: " + title, fontsize=20)
        fig.tight_layout(rect=(.03, .04, .98, .96))
        path = output_dir / f"xB_bin_{title_index + 1}_second_sector_profile.png"
        fig.savefig(path, dpi=200)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_second_sector_chi2_profiles(scan: pd.DataFrame, summary: pd.DataFrame, output_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(15, 11), sharex=True)
    axes = axes.ravel()
    for x_index, axis in enumerate(axes):
        subset = scan.loc[scan["x_index"] == x_index]
        for t_index in sorted(subset["t_index"].unique()):
            points = subset.loc[subset["t_index"] == t_index].sort_values("second_sector_fraction")
            mean_t = float(points["mean_minus_tprime_gev2"].iloc[0])
            axis.plot(points["second_sector_fraction"], points["chi2"],
                      lw=1.8, label=rf"$-t'={mean_t:.2f}$")
        threshold_points = summary.loc[summary["x_index"] == x_index, "acceptance_chi2_threshold"]
        if not threshold_points.empty:
            axis.axhline(float(threshold_points.iloc[0]), color="0.35", ls="--",
                         label=r"fixed-$r_B$ 5% threshold")
        axis.set_title(XB_LABELS[x_index])
        axis.set_xlabel(r"second-sector denominator fraction $r_B$")
        axis.set_ylabel(r"profile $\chi^2$")
        axis.set_xlim(0.0, float(scan["second_sector_fraction"].max()))
        axis.set_ylim(bottom=0.0)
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
        axis.legend(frameon=False, fontsize=8, ncol=2)
    fig.suptitle("Conditional second-sector strength profiles", fontsize=17)
    fig.tight_layout()
    path = output_dir / "second_sector_chi2_profiles.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def plot_second_sector_summary(summary: pd.DataFrame, output_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    for x_index, (color, marker) in enumerate(XB_STYLES):
        points = summary.loc[summary["x_index"] == x_index].sort_values("t_index")
        if points.empty:
            continue
        x = points["mean_minus_tprime_gev2"]
        axes[0].plot(x, points["minimum_acceptable_fraction"], marker=marker, color=color,
                     mfc="white", mew=2, label=XB_LABELS[x_index])
        axes[1].plot(x, points["best_profile_delta_chi2"], marker=marker, color=color,
                     mfc="white", mew=2)
    axes[0].set_ylabel(r"minimum acceptable $r_B$")
    axes[0].set_ylim(-.03, 1.03)
    axes[1].set_ylabel(r"$\chi^2_A-\min_{r_B}\chi^2_{A+B}$")
    axes[1].set_ylim(bottom=0.0)
    for axis in axes:
        axis.set_xlim(0.0, 1.20)
        axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$")
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("How much second-sector strength is required conditionally?", fontsize=16)
    fig.tight_layout()
    path = output_dir / "second_sector_minimum_strength_summary.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return path


def main() -> None:
    args = parse_args()
    if args.profile_points < 2:
        raise ValueError("--profile-points must be at least 2.")
    if not 0.0 < args.profile_max_fraction <= 1.0:
        raise ValueError("--profile-max-fraction must be greater than zero and no larger than one.")
    if not 0.0 < args.profile_acceptance_alpha < 1.0:
        raise ValueError("--profile-acceptance-alpha must lie strictly between zero and one.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.input_csv)
    covariances = source_covariances(args.input_json) if args.input_json.exists() else {}
    summary, detail = fit_frame(
        frame,
        covariances,
        args.default_epsilon,
        args.include_diagonal_systematics,
        args.starts,
        args.seed,
    )
    summary_path = args.output_dir / "no_flip_three_amplitude_fit_summary.csv"
    detail_path = args.output_dir / "no_flip_three_amplitude_observable_diagnostics.csv"
    summary.to_csv(summary_path, index=False)
    detail.to_csv(detail_path, index=False)
    figure_paths: list[Path] = []
    if not args.skip_plots:
        figure_paths.extend(plot_measured_vs_model(detail, args.output_dir))
        figure_paths.append(plot_goodness(summary, args.output_dir))
        figure_paths.append(plot_pull_heatmap(detail, args.output_dir))
        figure_paths.append(plot_amplitudes(summary, args.output_dir))
    profile_files: list[Path] = []
    profile_manifest = None
    if not args.skip_second_sector_profile:
        fractions = np.linspace(0.0, args.profile_max_fraction, args.profile_points)
        profile_summary, profile_scan, profile_detail = profile_second_sector(
            frame,
            covariances,
            detail,
            args.default_epsilon,
            args.include_diagonal_systematics,
            fractions,
            args.profile_starts,
            args.profile_acceptance_alpha,
            args.seed + 1000003,
        )
        profile_summary_path = args.output_dir / "second_sector_profile_summary.csv"
        profile_scan_path = args.output_dir / "second_sector_profile_scan.csv"
        profile_detail_path = args.output_dir / "second_sector_profile_observable_diagnostics.csv"
        profile_summary.to_csv(profile_summary_path, index=False)
        profile_scan.to_csv(profile_scan_path, index=False)
        profile_detail.to_csv(profile_detail_path, index=False)
        profile_files.extend([profile_summary_path, profile_scan_path, profile_detail_path])
        if not args.skip_plots:
            profile_files.extend(plot_second_sector_repairs(profile_detail, args.output_dir))
            profile_files.append(plot_second_sector_chi2_profiles(profile_scan, profile_summary, args.output_dir))
            profile_files.append(plot_second_sector_summary(profile_summary, args.output_dir))
        profile_manifest = {
            "definition": "r_B=D_B/(D_A+D_B), with D_s=0.5*(|P_s|^2+|M_s|^2)+epsilon*|L_s|^2",
            "prediction": "f_total=(1-r_B)*f_A+r_B*f_B",
            "conditioning": "The original best one-sector prediction f_A is fixed. Only the four observable shape parameters of sector B are refitted at each r_B.",
            "interpretation": "The minimum acceptable r_B is a grid-based conditional diagnostic. It is not a model-independent helicity-flip fraction or a confidence interval.",
            "recoil_sector_note": "For these recoil-unpolarized observables, a pure flip-only three-amplitude sector and a pure no-flip-only three-amplitude sector have the same ratio manifold. Sector B is therefore labeled generically as a second recoil-helicity sector.",
            "profile_points": int(args.profile_points),
            "profile_max_fraction": float(args.profile_max_fraction),
            "profile_starts": int(args.profile_starts),
            "acceptance_alpha": float(args.profile_acceptance_alpha),
            "files": list(map(str, profile_files)),
        }
    manifest = {
        "model": {
            "independent_amplitudes": {
                "L": "T_{0,+}^{0,+}",
                "P": "T_{+1,+}^{0,+}",
                "M": "T_{-1,+}^{0,+}",
            },
            "assumption": "All nucleon-helicity-flip amplitudes are zero; parity supplies the lambda=- amplitudes.",
            "real_parameters_per_bin": 4,
            "denominator": "0.5*(|P|^2+|M|^2)+epsilon*|L|^2; the common C cancels.",
            "interpretation": "A satisfactory fit means flip amplitudes are not required by these data. It does not prove that they vanish.",
        },
        "covariance": "full statistical covariance from JSON" if covariances else "diagonal statistical errors from CSV",
        "diagonal_systematics_included": bool(args.include_diagonal_systematics),
        "files": [str(summary_path), str(detail_path), *map(str, figure_paths), *map(str, profile_files)],
    }
    if profile_manifest is not None:
        manifest["second_sector_profile"] = profile_manifest
    (args.output_dir / "no_flip_three_amplitude_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    profile_note = " plus the conditional second-sector profile" if profile_manifest is not None else ""
    print(f"Wrote {len(summary)} bin fits{profile_note} to {args.output_dir}")


if __name__ == "__main__":
    main()
