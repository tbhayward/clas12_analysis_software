#!/usr/bin/env python3
"""Search for a sparse six-amplitude description of the seven RGC ratios.

The six independent complex amplitudes are represented as two recoil-helicity
sectors for initial-proton helicity +1/2,

    L+ = T_{0,+}^{0,+},   P+ = T_{+1,+}^{0,+},   M+ = T_{-1,+}^{0,+},
    L- = T_{0,+}^{0,-},   P- = T_{+1,+}^{0,-},   M- = T_{-1,+}^{0,-}.

Parity supplies the corresponding initial-proton-helicity -1/2 amplitudes.
Only ratios are fitted, so a common phase and scale are unobservable.  For
each candidate model one nonzero complex amplitude is therefore fixed to
1+0i.  A request for seven active real components then means the fixed real
anchor plus six fitted real components, with the other five of the original
twelve real/imaginary components set to zero in that phase convention.

The same component support is required in every kinematic bin.  Component
values are still fitted independently in each bin.  This is a sparse effective
model search, not a unique extraction of physical helicity amplitudes.
"""

from __future__ import annotations

import argparse
import itertools
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

from helicity_bilinear_transforms import SOURCE_PARAMETERS, TRANSFORM


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

AMPLITUDE_NAMES = ("Lp", "Pp", "Mp", "Lm", "Pm", "Mm")
AMPLITUDE_PDF_NAMES = {
    "Lp": "T_{0,+}^{0,+}",
    "Pp": "T_{+1,+}^{0,+}",
    "Mp": "T_{-1,+}^{0,+}",
    "Lm": "T_{0,+}^{0,-}",
    "Pm": "T_{+1,+}^{0,-}",
    "Mm": "T_{-1,+}^{0,-}",
}
COMPONENT_NAMES = tuple(
    f"{amplitude}_{part}"
    for amplitude in AMPLITUDE_NAMES
    for part in ("real", "imag")
)
COMPONENT_LABELS = {
    f"{name}_real": rf"$\mathrm{{Re}}\,{AMPLITUDE_PDF_NAMES[name]}$"
    for name in AMPLITUDE_NAMES
}
COMPONENT_LABELS.update({
    f"{name}_imag": rf"$\mathrm{{Im}}\,{AMPLITUDE_PDF_NAMES[name]}$"
    for name in AMPLITUDE_NAMES
})


@dataclass(frozen=True)
class SparseSupport:
    anchor: str
    active_components: tuple[str, ...]

    @property
    def fitted_components(self) -> tuple[str, ...]:
        return tuple(
            component for component in self.active_components
            if component != f"{self.anchor}_real"
        )

    @property
    def support_id(self) -> str:
        return self.anchor + "__" + "__".join(self.active_components)


@dataclass
class SparseFit:
    amplitudes: np.ndarray
    prediction: np.ndarray
    residual: np.ndarray
    marginal_pull: np.ndarray
    chi2: float
    covariance_rank: int
    model_rank: int
    dof: int
    p_value: float
    success: bool
    message: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Find sparse supports of the six parity-independent complex amplitudes "
            "that best reproduce all seven RGC structure-function ratios."
        )
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
        help="Production-fit JSON containing the full seven-ratio statistical covariance.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/sparse_helicity_amplitude_search"),
    )
    parser.add_argument(
        "--active-components",
        type=int,
        default=7,
        help=(
            "Number of nonzero real/imaginary components among the original twelve. "
            "One is the fixed real anchor, so the number fitted is one smaller."
        ),
    )
    parser.add_argument(
        "--anchor",
        choices=AMPLITUDE_NAMES,
        default="Lp",
        help="Amplitude fixed to 1+0i when --search-all-anchors is not used.",
    )
    parser.add_argument(
        "--search-all-anchors",
        action="store_true",
        help="Repeat the support search with each of the six amplitudes as the phase/scale anchor.",
    )
    parser.add_argument("--screen-starts", type=int, default=2)
    parser.add_argument("--refine-starts", type=int, default=20)
    parser.add_argument(
        "--refine-models",
        type=int,
        default=30,
        help="Number of best screened supports to refit with --refine-starts.",
    )
    parser.add_argument(
        "--max-screen-models",
        type=int,
        default=None,
        help="Optional debugging limit. Omit for the complete combinatorial search.",
    )
    parser.add_argument(
        "--include-diagonal-systematics",
        action="store_true",
        help="Add point-to-point systematic variances diagonally to the statistical covariance.",
    )
    parser.add_argument(
        "--default-epsilon",
        type=float,
        default=None,
        help="Fallback epsilon only when mean_epsilon is absent from the CSV.",
    )
    parser.add_argument("--seed", type=int, default=260827683)
    parser.add_argument("--skip-loo", action="store_true", help="Skip leave-one-observable-out checks.")
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def six_amplitude_bilinears(amplitudes: np.ndarray, epsilon: float) -> np.ndarray:
    """Return the seven normalized bilinears for both recoil-helicity sectors."""
    values = np.asarray(amplitudes, dtype=complex)
    if values.shape != (6,):
        raise ValueError("Six complex amplitudes are required.")
    Lp, Pp, Mp, Lm, Pm, Mm = values
    denominator = (
        0.5 * (abs(Pp) ** 2 + abs(Mp) ** 2 + abs(Pm) ** 2 + abs(Mm) ** 2)
        + epsilon * (abs(Lp) ** 2 + abs(Lm) ** 2)
    )
    if not np.isfinite(denominator) or denominator <= 0.0:
        raise ValueError("The unpolarized amplitude denominator must be positive.")
    lp = Lp * np.conjugate(Pp) + Lm * np.conjugate(Pm)
    lm = Lp * np.conjugate(Mp) + Lm * np.conjugate(Mm)
    mp = Mp * np.conjugate(Pp) + Mm * np.conjugate(Pm)
    return np.asarray(
        [
            np.imag(lp),
            -np.imag(lm),
            np.imag(mp),
            abs(Pp) ** 2 + abs(Pm) ** 2 - abs(Mp) ** 2 - abs(Mm) ** 2,
            np.real(lp),
            -np.real(lm),
            np.real(mp),
        ],
        dtype=float,
    ) / denominator


def predict_ratios(amplitudes: np.ndarray, epsilon: float) -> np.ndarray:
    """Return ratios in (u1,u2,lu1,ul1,ul2,ll0,ll1) order."""
    return np.linalg.solve(TRANSFORM, six_amplitude_bilinears(amplitudes, epsilon))


def amplitudes_from_parameters(parameters: np.ndarray, support: SparseSupport) -> np.ndarray:
    """Build six amplitudes in a fixed phase/scale convention."""
    parameters = np.asarray(parameters, dtype=float)
    if parameters.size != len(support.fitted_components):
        raise ValueError("Parameter count does not match the sparse support.")
    components = dict.fromkeys(COMPONENT_NAMES, 0.0)
    components[f"{support.anchor}_real"] = 1.0
    components[f"{support.anchor}_imag"] = 0.0
    for name, value in zip(support.fitted_components, parameters):
        components[name] = float(value)
    return np.asarray(
        [components[f"{name}_real"] + 1j * components[f"{name}_imag"] for name in AMPLITUDE_NAMES],
        dtype=complex,
    )


def displayed_amplitudes(amplitudes: np.ndarray) -> np.ndarray:
    """Normalize sum |T|^2 to one without changing the anchor phase convention."""
    values = np.asarray(amplitudes, dtype=complex).copy()
    norm = float(np.sqrt(np.sum(np.abs(values) ** 2)))
    return values / norm if norm > 0.0 else values


def covariance_whitener(covariance: np.ndarray) -> tuple[np.ndarray, int]:
    covariance = np.asarray(covariance, dtype=float)
    covariance = 0.5 * (covariance + covariance.T)
    values, vectors = np.linalg.eigh(covariance)
    tolerance = max(float(np.max(values)), 1.0) * 1.0e-12
    keep = values > tolerance
    if not np.any(keep):
        raise ValueError("Covariance has no positive eigenvalues.")
    return (vectors[:, keep] / np.sqrt(values[keep])).T, int(np.count_nonzero(keep))


def _initial_points(count: int, size: int, rng: np.random.Generator) -> list[np.ndarray]:
    points = [np.zeros(size)]
    if count > 1:
        points.append(rng.normal(0.0, 0.35, size=size))
    if count > 2:
        points.extend(rng.normal(0.0, 1.25, size=(count - 2, size)))
    return points[:count]


def fit_sparse_bin(
    measured: np.ndarray,
    covariance: np.ndarray,
    epsilon: float,
    support: SparseSupport,
    starts: int,
    rng: np.random.Generator,
    used_indices: np.ndarray | None = None,
) -> SparseFit:
    """Fit one support in one bin, optionally withholding some observables."""
    measured = np.asarray(measured, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if used_indices is None:
        used_indices = np.arange(measured.size)
    used_indices = np.asarray(used_indices, dtype=int)
    used_covariance = covariance[np.ix_(used_indices, used_indices)]
    whitener, covariance_rank = covariance_whitener(used_covariance)
    size = len(support.fitted_components)
    best = None
    for initial in _initial_points(max(starts, 1), size, rng):
        result = least_squares(
            lambda pars: whitener @ (
                measured[used_indices]
                - predict_ratios(amplitudes_from_parameters(pars, support), epsilon)[used_indices]
            ),
            initial,
            method="trf",
            max_nfev=4000,
            xtol=1.0e-11,
            ftol=1.0e-11,
            gtol=1.0e-11,
        )
        objective = float(np.dot(result.fun, result.fun))
        if best is None or objective < best[0]:
            best = (objective, result)
    assert best is not None
    objective, result = best
    amplitudes = amplitudes_from_parameters(result.x, support)
    prediction = predict_ratios(amplitudes, epsilon)
    residual = measured - prediction
    marginal_sigma = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    marginal_pull = np.divide(
        residual,
        marginal_sigma,
        out=np.full_like(residual, np.nan),
        where=marginal_sigma > 0.0,
    )
    singular_values = np.linalg.svd(result.jac, compute_uv=False)
    tolerance = max(result.jac.shape) * (singular_values[0] if singular_values.size else 0.0) * 1.0e-10
    model_rank = int(np.count_nonzero(singular_values > tolerance))
    dof = max(covariance_rank - model_rank, 0)
    return SparseFit(
        amplitudes=displayed_amplitudes(amplitudes),
        prediction=prediction,
        residual=residual,
        marginal_pull=marginal_pull,
        chi2=objective,
        covariance_rank=covariance_rank,
        model_rank=model_rank,
        dof=dof,
        p_value=float(chi2.sf(objective, dof)) if dof > 0 else np.nan,
        success=bool(result.success),
        message=str(result.message),
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


def point_to_point_errors(row: pd.Series) -> np.ndarray:
    errors = []
    for name in SOURCE_PARAMETERS:
        candidates = (
            f"{name}_point_to_point_systematic",
            f"{name}_point_to_point_systematic_diagonal",
        )
        value = next((float(row[key]) for key in candidates if key in row and pd.notna(row[key])), 0.0)
        errors.append(abs(value))
    return np.asarray(errors)


def epsilon_for_row(row: pd.Series, fallback: float | None) -> float:
    for column in ("mean_epsilon", "epsilon"):
        if column in row and pd.notna(row[column]):
            return float(row[column])
    if fallback is None:
        raise KeyError("CSV has no mean_epsilon column. Supply --default-epsilon explicitly.")
    return float(fallback)


def prepared_bins(
    frame: pd.DataFrame,
    covariances: dict[int, np.ndarray],
    default_epsilon: float | None,
    include_systematics: bool,
) -> list[dict]:
    bins = []
    for _, row in frame.sort_values("bin_number").iterrows():
        covariance = covariances.get(int(row["bin_number"]), diagonal_stat_covariance(row)).copy()
        systematic = point_to_point_errors(row)
        if include_systematics:
            covariance += np.diag(systematic**2)
        bins.append({
            "row": row,
            "measured": np.asarray([float(row[name]) for name in SOURCE_PARAMETERS]),
            "covariance": covariance,
            "stat_error": np.asarray([float(row[f"{name}_stat"]) for name in SOURCE_PARAMETERS]),
            "systematic_error": systematic,
            "epsilon": epsilon_for_row(row, default_epsilon),
        })
    return bins


def generate_supports(active_components: int, anchors: tuple[str, ...]) -> list[SparseSupport]:
    if not 2 <= active_components <= 11:
        raise ValueError("--active-components must be between 2 and 11 after fixing phase and scale.")
    supports = []
    for anchor in anchors:
        anchor_real = f"{anchor}_real"
        forbidden_imag = f"{anchor}_imag"
        available = [name for name in COMPONENT_NAMES if name not in (anchor_real, forbidden_imag)]
        for selected in itertools.combinations(available, active_components - 1):
            active = (anchor_real, *selected)
            supports.append(SparseSupport(anchor=anchor, active_components=tuple(active)))
    return supports


def evaluate_support(
    support: SparseSupport,
    bins: list[dict],
    starts: int,
    seed: int,
) -> tuple[dict, list[SparseFit]]:
    rng = np.random.default_rng(seed)
    fits = [
        fit_sparse_bin(
            item["measured"], item["covariance"], item["epsilon"], support, starts, rng
        )
        for item in bins
    ]
    total_chi2 = float(sum(fit.chi2 for fit in fits))
    total_model_rank = int(sum(fit.model_rank for fit in fits))
    total_covariance_rank = int(sum(fit.covariance_rank for fit in fits))
    total_dof = max(total_covariance_rank - total_model_rank, 0)
    number_measurements = total_covariance_rank
    aic = total_chi2 + 2.0 * total_model_rank
    bic = total_chi2 + total_model_rank * np.log(max(number_measurements, 1))
    row = {
        "support_id": support.support_id,
        "anchor": support.anchor,
        "anchor_pdf_name": AMPLITUDE_PDF_NAMES[support.anchor],
        "active_component_count": len(support.active_components),
        "fitted_parameter_count_per_bin": len(support.fitted_components),
        "active_components": ",".join(support.active_components),
        "zero_components": ",".join(name for name in COMPONENT_NAMES if name not in support.active_components),
        "total_chi2": total_chi2,
        "total_covariance_rank": total_covariance_rank,
        "total_model_rank": total_model_rank,
        "total_dof": total_dof,
        "chi2_per_dof": total_chi2 / total_dof if total_dof else np.nan,
        "global_p_value": float(chi2.sf(total_chi2, total_dof)) if total_dof else np.nan,
        "aic": aic,
        "bic": bic,
        "failed_bins": int(sum(not fit.success for fit in fits)),
        "rank_deficient_bins": int(sum(fit.model_rank < len(support.fitted_components) for fit in fits)),
    }
    return row, fits


def fit_details(support: SparseSupport, bins: list[dict], fits: list[SparseFit]) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    detail_rows = []
    for item, fit in zip(bins, fits):
        row = item["row"]
        common = {
            "bin_number": int(row["bin_number"]),
            "x_index": int(row["x_index"]),
            "t_index": int(row["t_index"]),
            "mean_minus_tprime_gev2": float(row["mean_minus_tprime_gev2"]),
            "mean_epsilon_used": item["epsilon"],
            "support_id": support.support_id,
        }
        summary = {
            **common,
            "chi2": fit.chi2,
            "covariance_rank": fit.covariance_rank,
            "model_rank": fit.model_rank,
            "dof": fit.dof,
            "chi2_per_dof": fit.chi2 / fit.dof if fit.dof else np.nan,
            "p_value": fit.p_value,
            "fit_success": fit.success,
            "fit_message": fit.message,
        }
        for name, amplitude in zip(AMPLITUDE_NAMES, fit.amplitudes):
            summary[f"{name}_real_normalized"] = float(amplitude.real)
            summary[f"{name}_imag_normalized"] = float(amplitude.imag)
            summary[f"{name}_magnitude_normalized"] = float(abs(amplitude))
            summary[f"{name}_phase_deg"] = float(np.degrees(np.angle(amplitude)))
        summary_rows.append(summary)
        total_error = np.sqrt(np.maximum(np.diag(item["covariance"]), 0.0))
        for index, observable in enumerate(SOURCE_PARAMETERS):
            detail_rows.append({
                **common,
                "observable": observable,
                "measured": item["measured"][index],
                "stat_error": item["stat_error"][index],
                "point_to_point_systematic": item["systematic_error"][index],
                "fit_uncertainty": total_error[index],
                "predicted": fit.prediction[index],
                "residual": fit.residual[index],
                "marginal_pull": fit.marginal_pull[index],
            })
    return pd.DataFrame(summary_rows), pd.DataFrame(detail_rows)


def leave_one_out(
    support: SparseSupport,
    bins: list[dict],
    starts: int,
    seed: int,
) -> pd.DataFrame:
    """Fit six ratios and predict the seventh for the selected support."""
    rng = np.random.default_rng(seed)
    rows = []
    all_indices = np.arange(len(SOURCE_PARAMETERS))
    for item in bins:
        source_row = item["row"]
        marginal_sigma = np.sqrt(np.maximum(np.diag(item["covariance"]), 0.0))
        for held_index, observable in enumerate(SOURCE_PARAMETERS):
            used = all_indices[all_indices != held_index]
            fit = fit_sparse_bin(
                item["measured"], item["covariance"], item["epsilon"], support, starts, rng, used
            )
            residual = item["measured"][held_index] - fit.prediction[held_index]
            rows.append({
                "bin_number": int(source_row["bin_number"]),
                "x_index": int(source_row["x_index"]),
                "t_index": int(source_row["t_index"]),
                "mean_minus_tprime_gev2": float(source_row["mean_minus_tprime_gev2"]),
                "held_out_observable": observable,
                "measured": item["measured"][held_index],
                "predicted_without_observable": fit.prediction[held_index],
                "held_out_residual": residual,
                "held_out_marginal_pull": residual / marginal_sigma[held_index]
                if marginal_sigma[held_index] > 0.0 else np.nan,
                "training_chi2": fit.chi2,
                "training_covariance_rank": fit.covariance_rank,
                "training_model_rank": fit.model_rank,
                "fit_success": fit.success,
            })
    return pd.DataFrame(rows)


def _style_axis(axis) -> None:
    axis.axhline(0.0, color="0.55", ls="--", lw=1.2, zorder=0)
    axis.set_xlim(0.0, 1.20)
    axis.minorticks_on()
    axis.tick_params(direction="in", top=True, right=True, which="both")


def _nice_tick_step(value: float) -> float:
    if not np.isfinite(value) or value <= 0.0:
        return 0.1
    power = 10.0 ** np.floor(np.log10(value))
    fraction = value / power
    for candidate in (1.0, 2.0, 2.5, 5.0, 10.0):
        if fraction <= candidate:
            return float(candidate * power)
    return float(10.0 * power)


def common_y_limits(detail: pd.DataFrame) -> dict[str, tuple[float, float]]:
    limits = {}
    for observable in SOURCE_PARAMETERS:
        points = detail.loc[detail["observable"] == observable]
        measured = points["measured"].to_numpy(float)
        uncertainty = points["fit_uncertainty"].to_numpy(float)
        predicted = points["predicted"].to_numpy(float)
        values = np.concatenate((measured - uncertainty, measured + uncertainty, predicted, [0.0]))
        values = values[np.isfinite(values)]
        low, high = float(np.min(values)), float(np.max(values))
        span = max(high - low, 0.05)
        step = _nice_tick_step(1.16 * span / 6.0)
        limits[observable] = (
            float(np.floor((low - 0.08 * span) / step) * step),
            float(np.ceil((high + 0.08 * span) / step) * step),
        )
    return limits


def plot_best_fit(detail: pd.DataFrame, output_dir: Path) -> list[Path]:
    paths = []
    y_limits = common_y_limits(detail)
    for x_index in sorted(detail["x_index"].unique()):
        subset = detail.loc[detail["x_index"] == x_index]
        color, marker = XB_STYLES[int(x_index)]
        figure, axes = plt.subplots(4, 2, figsize=(16, 20))
        axes = axes.ravel()
        for panel, observable in enumerate(SOURCE_PARAMETERS):
            axis = axes[panel]
            points = subset.loc[subset["observable"] == observable].sort_values("t_index")
            x = points["mean_minus_tprime_gev2"].to_numpy(float)
            y = points["measured"].to_numpy(float)
            stat = points["stat_error"].to_numpy(float)
            syst = points["point_to_point_systematic"].to_numpy(float)
            _style_axis(axis)
            axis.set_ylim(*y_limits[observable])
            width = 0.045
            for x_value, y_value, error in zip(x, y, syst):
                axis.add_patch(Rectangle(
                    (x_value - width / 2.0, y_value - error), width, 2.0 * error,
                    facecolor=color, edgecolor="none", alpha=0.18, zorder=1,
                ))
            axis.errorbar(
                x, y, yerr=stat, fmt=marker, mfc="white", mec=color, mew=2,
                ecolor=color, capsize=3, label="Measured ratio (stat.)", zorder=3,
            )
            axis.plot(x, points["predicted"], color="black", lw=2, marker=".",
                      label="Best sparse model", zorder=2)
            axis.set_ylabel(SOURCE_LABELS[observable], fontsize=15, labelpad=12)
            axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$", fontsize=12)
            axis.text(.04, .90, f"({chr(97 + panel)})", transform=axis.transAxes,
                      fontsize=14, fontweight="bold")
        axes[-1].axis("off")
        handles, labels = axes[0].get_legend_handles_labels()
        handles.insert(0, Rectangle((0, 0), 1, 1, facecolor=color, alpha=0.18))
        labels.insert(0, "Point-to-point syst.")
        figure.legend(handles, labels, loc="lower right", bbox_to_anchor=(.95, .05), frameon=False)
        title = XB_LABELS[int(x_index)] if int(x_index) < len(XB_LABELS) else f"x bin {x_index}"
        figure.suptitle("Best common sparse-amplitude support: " + title, fontsize=20)
        figure.tight_layout(rect=(.03, .04, .98, .96))
        path = output_dir / f"xB_bin_{int(x_index)+1}_best_sparse_model.png"
        figure.savefig(path, dpi=200)
        plt.close(figure)
        paths.append(path)
    return paths


def plot_goodness(summary: pd.DataFrame, output_dir: Path) -> Path:
    figure, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    for x_index, (color, marker) in enumerate(XB_STYLES):
        points = summary.loc[summary["x_index"] == x_index].sort_values("t_index")
        if points.empty:
            continue
        x = points["mean_minus_tprime_gev2"]
        axes[0].plot(x, points["chi2_per_dof"], marker=marker, color=color,
                     mfc="white", mew=2, label=XB_LABELS[x_index])
        axes[1].plot(x, points["p_value"], marker=marker, color=color, mfc="white", mew=2)
    axes[0].axhline(1.0, color="0.4", ls="--")
    axes[0].set_ylabel(r"$\chi^2/\mathrm{dof}$")
    axes[1].axhline(0.05, color="0.4", ls="--")
    axes[1].set_ylabel("fit p-value")
    axes[1].set_ylim(-0.02, 1.02)
    for axis in axes:
        axis.set_xlim(0.0, 1.20)
        axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^2)$")
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
    axes[0].legend(frameon=False, fontsize=9)
    figure.suptitle("Goodness of fit of the selected sparse support", fontsize=16)
    figure.tight_layout()
    path = output_dir / "best_sparse_model_goodness_of_fit.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def plot_model_ranking(ranking: pd.DataFrame, output_dir: Path) -> Path:
    points = ranking.loc[ranking["stage"] == "refined"].nsmallest(30, "total_chi2").copy()
    points = points.sort_values("total_chi2", ascending=True)
    figure, axis = plt.subplots(figsize=(11, 8))
    labels = [f"{row.anchor}: {row.active_components}" for row in points.itertuples()]
    y = np.arange(len(points))
    axis.barh(y, points["total_chi2"], color="#4C78A8")
    axis.set_yticks(y, labels, fontsize=7)
    axis.invert_yaxis()
    axis.set_xlabel(r"total $\chi^2$ across all bins")
    axis.set_title("Refined sparse-support ranking")
    figure.tight_layout()
    path = output_dir / "refined_sparse_support_ranking.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def plot_component_frequency(ranking: pd.DataFrame, output_dir: Path) -> Path:
    top = ranking.loc[ranking["stage"] == "refined"].nsmallest(20, "total_chi2")
    counts = {
        component: sum(component in value.split(",") for value in top["active_components"])
        for component in COMPONENT_NAMES
    }
    figure, axis = plt.subplots(figsize=(12, 5.5))
    names = list(COMPONENT_NAMES)
    axis.bar(range(len(names)), [counts[name] for name in names], color="#59A14F")
    axis.set_xticks(range(len(names)), [COMPONENT_LABELS[name] for name in names], rotation=45, ha="right")
    axis.set_ylabel("appearances among best 20 refined supports")
    axis.set_ylim(0, max(len(top), 1))
    axis.set_title("Support stability diagnostic")
    figure.tight_layout()
    path = output_dir / "top_model_component_frequency.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def plot_pull_heatmap(detail: pd.DataFrame, output_dir: Path) -> Path:
    bin_order = detail[["bin_number", "x_index", "t_index"]].drop_duplicates().sort_values(
        ["x_index", "t_index"]
    )
    matrix = np.full((len(SOURCE_PARAMETERS), len(bin_order)), np.nan)
    for column, bin_number in enumerate(bin_order["bin_number"]):
        points = detail.loc[detail["bin_number"] == bin_number].set_index("observable")
        for row, observable in enumerate(SOURCE_PARAMETERS):
            matrix[row, column] = float(points.loc[observable, "marginal_pull"])
    figure, axis = plt.subplots(figsize=(15, 5.5))
    image = axis.imshow(matrix, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    axis.set_yticks(range(len(SOURCE_PARAMETERS)), [SOURCE_LABELS[name] for name in SOURCE_PARAMETERS])
    labels = [f"x{int(x)+1},t{int(t)+1}" for x, t in zip(bin_order["x_index"], bin_order["t_index"])]
    axis.set_xticks(range(len(labels)), labels, rotation=60, ha="right")
    axis.set_title("Best sparse model marginal residuals")
    figure.colorbar(image, ax=axis, label="(measurement - model) / marginal uncertainty")
    figure.tight_layout()
    path = output_dir / "best_sparse_model_pull_heatmap.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def plot_loo_heatmap(loo: pd.DataFrame, output_dir: Path) -> Path:
    bin_order = loo[["bin_number", "x_index", "t_index"]].drop_duplicates().sort_values(
        ["x_index", "t_index"]
    )
    matrix = np.full((len(SOURCE_PARAMETERS), len(bin_order)), np.nan)
    for column, bin_number in enumerate(bin_order["bin_number"]):
        points = loo.loc[loo["bin_number"] == bin_number].set_index("held_out_observable")
        for row, observable in enumerate(SOURCE_PARAMETERS):
            matrix[row, column] = float(points.loc[observable, "held_out_marginal_pull"])
    figure, axis = plt.subplots(figsize=(15, 5.5))
    image = axis.imshow(matrix, cmap="RdBu_r", vmin=-4, vmax=4, aspect="auto")
    axis.set_yticks(range(len(SOURCE_PARAMETERS)), [SOURCE_LABELS[name] for name in SOURCE_PARAMETERS])
    labels = [f"x{int(x)+1},t{int(t)+1}" for x, t in zip(bin_order["x_index"], bin_order["t_index"])]
    axis.set_xticks(range(len(labels)), labels, rotation=60, ha="right")
    axis.set_title("Leave-one-observable-out prediction residuals")
    figure.colorbar(image, ax=axis, label="held-out marginal pull")
    figure.tight_layout()
    path = output_dir / "leave_one_observable_out_pull_heatmap.png"
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return path


def main() -> None:
    args = parse_args()
    if args.screen_starts < 1 or args.refine_starts < 1:
        raise ValueError("Multistart counts must be positive.")
    if args.refine_models < 1:
        raise ValueError("--refine-models must be positive.")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.input_csv)
    missing = [name for name in SOURCE_PARAMETERS if name not in frame.columns]
    if missing:
        raise KeyError(f"Input CSV is missing ratio columns: {missing}")
    covariances = source_covariances(args.input_json) if args.input_json.exists() else {}
    bins = prepared_bins(
        frame, covariances, args.default_epsilon, args.include_diagonal_systematics
    )
    anchors = AMPLITUDE_NAMES if args.search_all_anchors else (args.anchor,)
    supports = generate_supports(args.active_components, anchors)
    if args.max_screen_models is not None:
        supports = supports[:args.max_screen_models]
    print(f"Screening {len(supports)} sparse supports across {len(bins)} bins...")
    screening_rows = []
    support_lookup = {support.support_id: support for support in supports}
    for index, support in enumerate(supports):
        row, _ = evaluate_support(
            support, bins, args.screen_starts, args.seed + 1009 * index
        )
        row["stage"] = "screen"
        screening_rows.append(row)
        if (index + 1) % 50 == 0 or index + 1 == len(supports):
            print(f"  screened {index + 1}/{len(supports)}")
    screening = pd.DataFrame(screening_rows).sort_values("total_chi2")
    refine_ids = screening.head(min(args.refine_models, len(screening)))["support_id"].tolist()
    print(f"Refining the best {len(refine_ids)} supports...")
    refined_rows = []
    refined_fits: dict[str, list[SparseFit]] = {}
    for index, support_id in enumerate(refine_ids):
        support = support_lookup[support_id]
        row, fits = evaluate_support(
            support, bins, args.refine_starts, args.seed + 1000003 + 7919 * index
        )
        row["stage"] = "refined"
        refined_rows.append(row)
        refined_fits[support_id] = fits
    refined = pd.DataFrame(refined_rows).sort_values("total_chi2")
    ranking = pd.concat([screening, refined], ignore_index=True)
    ranking_path = args.output_dir / "sparse_support_ranking.csv"
    ranking.to_csv(ranking_path, index=False)
    best_row = refined.iloc[0]
    best_support = support_lookup[str(best_row["support_id"])]
    best_fits = refined_fits[best_support.support_id]
    summary, detail = fit_details(best_support, bins, best_fits)
    summary_path = args.output_dir / "best_sparse_model_bin_summary.csv"
    detail_path = args.output_dir / "best_sparse_model_observable_diagnostics.csv"
    summary.to_csv(summary_path, index=False)
    detail.to_csv(detail_path, index=False)
    files = [ranking_path, summary_path, detail_path]
    loo_path = None
    loo = None
    if not args.skip_loo:
        print("Running leave-one-observable-out checks for the selected support...")
        loo = leave_one_out(
            best_support, bins, args.refine_starts, args.seed + 2000003
        )
        loo_path = args.output_dir / "leave_one_observable_out_diagnostics.csv"
        loo.to_csv(loo_path, index=False)
        files.append(loo_path)
    if not args.skip_plots:
        files.extend(plot_best_fit(detail, args.output_dir))
        files.append(plot_goodness(summary, args.output_dir))
        files.append(plot_model_ranking(ranking, args.output_dir))
        files.append(plot_component_frequency(ranking, args.output_dir))
        files.append(plot_pull_heatmap(detail, args.output_dir))
        if loo is not None:
            files.append(plot_loo_heatmap(loo, args.output_dir))
    manifest = {
        "purpose": (
            "Sparse effective-amplitude model search. This does not uniquely extract physical "
            "helicity amplitudes or establish that omitted components vanish."
        ),
        "amplitudes": AMPLITUDE_PDF_NAMES,
        "component_convention": (
            "For each candidate support, the anchor amplitude is fixed to 1+0i to remove the "
            "common scale and phase that cancel in the seven ratios."
        ),
        "active_component_count": args.active_components,
        "fitted_parameter_count_per_bin": args.active_components - 1,
        "same_support_in_all_bins": True,
        "independent_values_in_each_bin": True,
        "searched_anchors": list(anchors),
        "screened_support_count": len(supports),
        "refined_support_count": len(refined),
        "best_support": {
            "support_id": best_support.support_id,
            "anchor": best_support.anchor,
            "anchor_pdf_name": AMPLITUDE_PDF_NAMES[best_support.anchor],
            "active_components": list(best_support.active_components),
            "zero_components": [
                name for name in COMPONENT_NAMES if name not in best_support.active_components
            ],
            "total_chi2": float(best_row["total_chi2"]),
            "total_dof": int(best_row["total_dof"]),
            "global_p_value": float(best_row["global_p_value"]),
        },
        "covariance": (
            "full statistical covariance from JSON"
            if covariances else "diagonal statistical covariance from CSV"
        ),
        "diagonal_systematics_included": bool(args.include_diagonal_systematics),
        "leave_one_observable_out_note": (
            "Held-out pulls use the held observable's marginal uncertainty. They are a practical "
            "predictive diagnostic, not a full uncertainty band on the nonlinear prediction."
        ),
        "files": list(map(str, files)),
    }
    manifest_path = args.output_dir / "sparse_helicity_amplitude_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("Best support:")
    print(f"  anchor: {best_support.anchor} = {AMPLITUDE_PDF_NAMES[best_support.anchor]} = 1+0i")
    print(f"  active: {', '.join(best_support.active_components)}")
    print(f"  total chi2/dof: {best_row['total_chi2']:.3f}/{int(best_row['total_dof'])}")
    print(f"Wrote sparse-amplitude study to {args.output_dir}")


if __name__ == "__main__":
    main()
