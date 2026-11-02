#!/usr/bin/env python3
"""Convert the production seven-parameter RGC fit to helicity bilinears.

This is deliberately a change of basis of the existing simultaneous fit, not
a second independent minimization.  Because the transformation is linear and
invertible, it has the same minimum likelihood and preserves all statistical
information when the production covariance matrices are supplied.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

from helicity_bilinear_transforms import (
    HELICITY_DESCRIPTIONS,
    HELICITY_LABELS,
    HELICITY_PARAMETERS,
    SOURCE_PARAMETERS,
    TRANSFORM,
    transform_covariance,
    transform_values,
)


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

COMPARABLE_PANELS = ((0, 1), (2,), (3,), (4, 5), (6,))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Transform the seven RGC structure-function ratios into seven normalized helicity bilinears."
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
        help="Detailed production JSON used to recover the exact parameter order and covariance.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/helicity_bilinear_extraction"),
    )
    parser.add_argument("--skip-plots", action="store_true")
    return parser.parse_args()


def source_covariances(path: Path) -> dict[int, np.ndarray]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    output: dict[int, np.ndarray] = {}
    for result in payload["results"]:
        fit = result["variants"]["nominal"]
        order = list(fit["parameter_order"])
        covariance = np.asarray(fit["covariance"], dtype=float)
        indices = [order.index(name) for name in SOURCE_PARAMETERS]
        output[int(result["bin_number"])] = covariance[np.ix_(indices, indices)]
    return output


def diagonal_covariance_from_frame(row: pd.Series) -> np.ndarray:
    errors = np.asarray([float(row[f"{name}_stat"]) for name in SOURCE_PARAMETERS])
    return np.diag(errors**2)


def propagated_diagonal_uncertainty(row: pd.Series, suffix: str) -> np.ndarray:
    errors = np.asarray(
        [float(row.get(f"{name}_{suffix}", 0.0)) for name in SOURCE_PARAMETERS],
        dtype=float,
    )
    return np.sqrt((TRANSFORM**2) @ (errors**2))


def _nice_step(span: float) -> float:
    if span <= 0.25:
        return 0.05
    if span <= 0.60:
        return 0.10
    if span <= 1.50:
        return 0.20
    return 0.50


def standardized_limits(frame: pd.DataFrame) -> dict[int, tuple[float, float]]:
    """Set common limits for like bilinears, using only this canvas's data."""
    limits: dict[int, tuple[float, float]] = {}
    for panel_group in COMPARABLE_PANELS:
        lows, highs = [0.0], [0.0]
        for panel in panel_group:
            name = HELICITY_PARAMETERS[panel]
            y = frame[name].to_numpy(dtype=float)
            stat = frame[f"{name}_stat"].to_numpy(dtype=float)
            systematic_name = f"{name}_point_to_point_systematic_diagonal"
            syst = frame[systematic_name].to_numpy(dtype=float) if systematic_name in frame else np.zeros_like(y)
            extent = np.maximum(stat, syst)
            lows.extend((y - extent).tolist())
            highs.extend((y + extent).tolist())
        low, high = min(lows), max(highs)
        padding = max(0.08 * (high - low), 0.02)
        step = _nice_step(high - low + 2.0 * padding)
        low = step * np.floor((low - padding) / step)
        high = step * np.ceil((high + padding) / step)
        for panel in panel_group:
            limits[panel] = (low, high)
    return limits


def draw_systematic_boxes(axis, x, y, syst, color, width=.045, label=None):
    for index, (x_value, y_value, error) in enumerate(zip(x, y, syst)):
        axis.add_patch(
            Rectangle(
                (x_value - width / 2.0, y_value - error), width, 2.0 * error,
                facecolor=color, edgecolor="none", alpha=.18,
                label=label if index == 0 else None, zorder=1,
            )
        )


def style_axis(axis, panel: int, limits: dict[int, tuple[float, float]]):
    name = HELICITY_PARAMETERS[panel]
    axis.axhline(0.0, color="0.55", ls="--", lw=1.4, zorder=0)
    axis.set_xlim(0.0, 1.20)
    axis.set_ylim(*limits[panel])
    axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^{2})$", fontsize=12)
    axis.set_ylabel(HELICITY_LABELS[name], fontsize=8.5, labelpad=8)
    axis.text(0.04, 0.91, f"({chr(97 + panel)})", transform=axis.transAxes,
              fontsize=14, fontweight="bold")
    axis.minorticks_on()
    axis.tick_params(direction="in", top=True, right=True, which="both")


def plot_x_group(frame: pd.DataFrame, x_index: int, output_dir: Path) -> Path:
    subset = frame.loc[frame["x_index"].astype(int) == x_index].sort_values("t_index")
    x = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
    limits = standardized_limits(subset)
    color, marker = XB_STYLES[x_index]
    fig, axes = plt.subplots(4, 2, figsize=(18, 22), sharex=False)
    axes = axes.ravel()
    for panel, (axis, name) in enumerate(zip(axes, HELICITY_PARAMETERS)):
        y = subset[name].to_numpy(dtype=float)
        stat = subset[f"{name}_stat"].to_numpy(dtype=float)
        style_axis(axis, panel, limits)
        systematic_column = f"{name}_point_to_point_systematic_diagonal"
        if systematic_column in subset:
            systematic = subset[systematic_column].to_numpy(dtype=float)
            draw_systematic_boxes(
                axis, x, y, systematic, color,
                label="P-to-p syst. (diagonal propagation)",
            )
        axis.errorbar(x, y, yerr=stat, fmt=marker, ms=7, mfc="white",
                      mec=color, mew=2, ecolor=color, capsize=4,
                      label="RGC transformed amplitude combination", zorder=2)
    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.94, 0.045), frameon=False)
    fig.suptitle(
        r"CLAS12 $ep\rightarrow e(n)\pi^+$ helicity-amplitude combinations: " + XB_LABELS[x_index],
        fontsize=20,
    )
    fig.text(
        0.5, 0.015,
        r"$T_{\mu,\lambda}^{0\sigma}$: photon helicity $\mu$, initial-proton helicity $\lambda$, "
        r"pion helicity 0, recoil-neutron helicity $\sigma$; sum is over $\sigma=\pm\frac{1}{2}$.",
        ha="center", fontsize=11,
    )
    fig.tight_layout(rect=(0.03, 0.05, 0.98, 0.96))
    path = output_dir / f"xB_bin_{x_index + 1}_seven_helicity_bilinears.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def plot_all_x_groups(frame: pd.DataFrame, output_dir: Path) -> Path:
    limits = standardized_limits(frame)
    fig, axes = plt.subplots(4, 2, figsize=(18, 22), sharex=False)
    axes = axes.ravel()
    for panel, (axis, name) in enumerate(zip(axes, HELICITY_PARAMETERS)):
        style_axis(axis, panel, limits)
        for x_index, (color, marker) in enumerate(XB_STYLES):
            subset = frame.loc[frame["x_index"].astype(int) == x_index].sort_values("t_index")
            x = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
            y = subset[name].to_numpy(dtype=float)
            stat = subset[f"{name}_stat"].to_numpy(dtype=float)
            systematic_column = f"{name}_point_to_point_systematic_diagonal"
            if systematic_column in subset:
                syst = subset[systematic_column].to_numpy(dtype=float)
                draw_systematic_boxes(axis, x, y, syst, color)
            axis.errorbar(
                x, y, yerr=stat, fmt=marker, ms=6.5, mfc="white", mec=color,
                mew=1.9, ecolor=color, capsize=3, label=XB_LABELS[x_index], zorder=2,
            )
    axes[-1].axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    axes[-1].legend(handles, labels, loc="center", frameon=False, fontsize=12,
                    title=r"$x_B$ bins", title_fontsize=13)
    axes[-1].text(
        .5, .20,
        "Statistical error bars; colored rectangles are\n"
        "point-to-point systematics (diagonal propagation)",
        ha="center", va="center", transform=axes[-1].transAxes, fontsize=11,
    )
    fig.suptitle(
        r"CLAS12 $ep\rightarrow e(n)\pi^+$ helicity-amplitude combinations: all $x_B$ bins",
        fontsize=20,
    )
    fig.text(
        0.5, 0.015,
        r"$T_{\mu,\lambda}^{0\sigma}$: photon helicity $\mu$, initial-proton helicity $\lambda$, "
        r"pion helicity 0, recoil-neutron helicity $\sigma$; sum is over $\sigma=\pm\frac{1}{2}$.",
        ha="center", fontsize=11,
    )
    fig.tight_layout(rect=(0.03, 0.05, 0.98, 0.96))
    path = output_dir / "all_xB_seven_helicity_bilinears.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def main() -> int:
    args = parse_args()
    output = args.output_dir.expanduser().resolve()
    tables = output / "tables"
    covariance_dir = output / "covariance"
    plots = output / "plots"
    for directory in (output, tables, covariance_dir, plots):
        directory.mkdir(parents=True, exist_ok=True)

    frame = pd.read_csv(args.input_csv.expanduser().resolve()).sort_values("bin_number").reset_index(drop=True)
    missing = [name for name in SOURCE_PARAMETERS if name not in frame.columns]
    if missing:
        raise RuntimeError(f"Missing source columns: {missing}")

    exact_covariances: dict[int, np.ndarray] = {}
    covariance_source = "diagonal approximation from marginal statistical errors"
    json_path = args.input_json.expanduser().resolve()
    if json_path.is_file():
        exact_covariances = source_covariances(json_path)
        covariance_source = "exact nominal Minuit covariance from production JSON"

    values = frame.loc[:, SOURCE_PARAMETERS].to_numpy(dtype=float)
    transformed = transform_values(values)
    result = frame.copy()
    for index, name in enumerate(HELICITY_PARAMETERS):
        result[name] = transformed[:, index]

    for row_index, row in result.iterrows():
        bin_number = int(row["bin_number"])
        source_covariance = exact_covariances.get(bin_number)
        if source_covariance is None:
            source_covariance = diagonal_covariance_from_frame(row)
        helicity_covariance = transform_covariance(source_covariance)
        helicity_correlation = np.divide(
            helicity_covariance,
            np.sqrt(np.outer(np.diag(helicity_covariance), np.diag(helicity_covariance))),
            out=np.zeros_like(helicity_covariance),
            where=np.outer(np.diag(helicity_covariance), np.diag(helicity_covariance)) > 0.0,
        )
        np.save(covariance_dir / f"bin_{bin_number:02d}_helicity_covariance.npy", helicity_covariance)
        np.save(covariance_dir / f"bin_{bin_number:02d}_helicity_correlation.npy", helicity_correlation)
        for index, name in enumerate(HELICITY_PARAMETERS):
            result.at[row_index, f"{name}_stat"] = np.sqrt(max(helicity_covariance[index, index], 0.0))

        # The production CSV stores systematic magnitudes, not their signed
        # covariance.  These columns therefore use explicitly labelled
        # diagonal propagation.  Exact systematic propagation should transform
        # each signed alternative-fit vector before taking its magnitude.
        for suffix in (
            "radiation_systematic",
            "channel_selection_systematic",
            "bin_migration_systematic",
            "point_to_point_systematic",
        ):
            propagated = propagated_diagonal_uncertainty(row, suffix)
            for index, name in enumerate(HELICITY_PARAMETERS):
                result.at[row_index, f"{name}_{suffix}_diagonal"] = propagated[index]

    csv_path = tables / "helicity_bilinear_combinations.csv"
    result.to_csv(csv_path, index=False)
    definitions_path = tables / "helicity_bilinear_definitions.json"
    definitions_path.write_text(
        json.dumps(
            {
                "source_parameter_order": SOURCE_PARAMETERS,
                "helicity_parameter_order": HELICITY_PARAMETERS,
                "transformation_matrix": TRANSFORM.tolist(),
                "covariance_rule": "V_H = M V_F M^T",
                "covariance_source": covariance_source,
                "definitions": {
                    name: {
                        "matplotlib_label": HELICITY_LABELS[name],
                        "description": HELICITY_DESCRIPTIONS[name],
                    }
                    for name in HELICITY_PARAMETERS
                },
                "important": (
                    "These are normalized bilinear combinations, not individual complex amplitudes. "
                    "Labels retain the PDF convention T_{mu,lambda}^{0 sigma}."
                ),
            },
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )

    plot_paths = []
    if not args.skip_plots:
        for x_index in range(4):
            plot_paths.append(str(plot_x_group(result, x_index, plots)))
        plot_paths.append(str(plot_all_x_groups(result, plots)))

    manifest = {
        "input_csv": str(args.input_csv.expanduser().resolve()),
        "input_json": str(json_path) if json_path.is_file() else None,
        "output_csv": str(csv_path),
        "definitions": str(definitions_path),
        "covariance_directory": str(covariance_dir),
        "plots": plot_paths,
        "method": "exact linear change of basis of the seven-parameter production fit",
    }
    (output / "helicity_bilinear_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote {csv_path}")
    print(f"Statistical covariance: {covariance_source}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
