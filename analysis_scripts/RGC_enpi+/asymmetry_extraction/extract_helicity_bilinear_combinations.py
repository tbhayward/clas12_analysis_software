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


def plot_x_group(frame: pd.DataFrame, x_index: int, output_dir: Path) -> Path:
    subset = frame.loc[frame["x_index"].astype(int) == x_index].sort_values("t_index")
    x = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
    fig, axes = plt.subplots(4, 2, figsize=(17, 22), sharex=True)
    axes = axes.ravel()
    for panel, (axis, name) in enumerate(zip(axes, HELICITY_PARAMETERS), start=1):
        y = subset[name].to_numpy(dtype=float)
        stat = subset[f"{name}_stat"].to_numpy(dtype=float)
        axis.axhline(0.0, color="0.55", ls="--", lw=1.4)
        systematic_column = f"{name}_point_to_point_systematic_diagonal"
        if systematic_column in subset:
            systematic = subset[systematic_column].to_numpy(dtype=float)
            width = 0.045
            axis.bar(x, 2.0 * systematic, bottom=-systematic, width=width,
                     color="#9ecae1", alpha=0.48, edgecolor="none",
                     label="P-to-p syst. (diagonal propagation)")
        axis.errorbar(x, y, yerr=stat, fmt="o", ms=7, mfc="white",
                      mec="#0072B2", mew=2, ecolor="#0072B2", capsize=4,
                      label="RGC transformed fit")
        axis.set_title(HELICITY_LABELS[name], fontsize=11, pad=9)
        axis.set_ylabel("normalized bilinear")
        axis.text(0.04, 0.91, f"({chr(96 + panel)})", transform=axis.transAxes,
                  fontsize=14, fontweight="bold")
        axis.minorticks_on()
        axis.tick_params(direction="in", top=True, right=True, which="both")
    axes[-1].axis("off")
    for axis in axes[-2:]:
        if axis.axison:
            axis.set_xlabel(r"$-t'\;(\mathrm{GeV}^{2})$")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.94, 0.045), frameon=False)
    fig.suptitle(
        r"CLAS12 $ep\rightarrow e(n)\pi^+$ normalized helicity bilinears: " + XB_LABELS[x_index],
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
