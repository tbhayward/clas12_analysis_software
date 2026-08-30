#!/usr/bin/env python3
"""
Compare photon-FD-sector dependence of the DVCS cross section with the
independent pi0-BDT tag-and-probe photon-efficiency correction.

The six DVCS CSVs are supplied explicitly at runtime.  The efficiency input is
the `nominal_fd_corrections.csv` table written by
derive_photon_efficiency_from_pi0_bdt.py.

The comparison deliberately uses sector-normalized shapes.  This removes an
irrelevant overall normalization and asks whether the independently observed
six-sector patterns are correlated.

Example
-------
python compare_dvcs_photon_sector_to_efficiency.py \
    --sec1 phot_sec1.csv --sec2 phot_sec2.csv --sec3 phot_sec3.csv \
    --sec4 phot_sec4.csv --sec5 phot_sec5.csv --sec6 phot_sec6.csv \
    --efficiency-table output/photon_efficiency_bdt/fa18_inb/full/nominal_fd_corrections.csv \
    --period "Fa18 Inb"
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_CROSS_SECTION_COLUMN = (
    "cross sections, ep->epg, exp, Fa18 Inb, unpol"
)
DEFAULT_GTHETA_COLUMN = "g_theta, Fa18 Inb"
DEFAULT_DETAILED_THETA_MIN = 5.5
DEFAULT_DETAILED_THETA_MAX = 30.0
DEFAULT_DETAILED_THETA_BINS = 7

# These match the current photon-efficiency correction bins.
MATCHED_THETA_BINS = (
    (5.5, 10.0),
    (10.0, 15.0),
    (15.0, 35.0),
)
HIGH_ENERGY_BINS = (
    (2.0, 3.0),
    (3.0, 9.5),
)

VOLUME_COLUMNS = (
    ("xBmin", "xBmax"),
    ("Q2min", "Q2max"),
    ("t_abs_min", "t_abs_max"),
    ("phimin", "phimax"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare six photon-sector DVCS cross sections with the "
            "sector-dependent tag-and-probe photon-efficiency correction."
        )
    )

    for sector in range(1, 7):
        parser.add_argument(
            f"--sec{sector}",
            type=Path,
            required=True,
            help=f"DVCS photon-sector-{sector} CSV.",
        )
    #endfor

    parser.add_argument(
        "--efficiency-table",
        type=Path,
        required=True,
        help=(
            "nominal_fd_corrections.csv written by the photon-efficiency "
            "script."
        ),
    )
    parser.add_argument(
        "--period",
        default="Fa18 Inb",
        help="Human-readable period label used in plot titles.",
    )
    parser.add_argument(
        "--cross-section-column",
        default=DEFAULT_CROSS_SECTION_COLUMN,
        help="DVCS CSV column containing the unpolarized cross section.",
    )
    parser.add_argument(
        "--gtheta-column",
        default=DEFAULT_GTHETA_COLUMN,
        help="DVCS CSV column containing photon polar angle in degrees.",
    )
    parser.add_argument(
        "--detailed-theta-min",
        type=float,
        default=DEFAULT_DETAILED_THETA_MIN,
    )
    parser.add_argument(
        "--detailed-theta-max",
        type=float,
        default=DEFAULT_DETAILED_THETA_MAX,
    )
    parser.add_argument(
        "--detailed-theta-bins",
        type=int,
        default=DEFAULT_DETAILED_THETA_BINS,
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output") / "photon_sector_dvcs_closure",
    )

    return parser.parse_args()


def tuple_value_unc(cell: object) -> Tuple[float, float]:
    """
    Parse cells such as "(value,stat,sys)" or "(value,stat)".

    Only the first uncertainty is used here.  The closure plots are intended
    as a statistical/shape diagnostic rather than a final systematic error
    propagation.
    """
    if cell is None:
        return np.nan, np.nan
    #endif

    if isinstance(cell, float) and np.isnan(cell):
        return np.nan, np.nan
    #endif

    if isinstance(cell, (int, float, np.integer, np.floating)):
        return float(cell), np.nan
    #endif

    text = str(cell).strip()
    if not text:
        return np.nan, np.nan
    #endif

    try:
        parsed = ast.literal_eval(text)
    except Exception:
        try:
            return float(text), np.nan
        except Exception:
            return np.nan, np.nan
        #endtry
    #endtry

    if isinstance(parsed, (tuple, list)):
        value = float(parsed[0]) if len(parsed) >= 1 else np.nan
        unc = float(parsed[1]) if len(parsed) >= 2 else np.nan
        return value, unc
    #endif

    try:
        return float(parsed), np.nan
    except Exception:
        return np.nan, np.nan
    #endtry


def bin_volume(df: pd.DataFrame) -> np.ndarray:
    volume = np.ones(len(df), dtype=float)

    for low_col, high_col in VOLUME_COLUMNS:
        if low_col not in df.columns or high_col not in df.columns:
            raise KeyError(
                f"Required DVCS bin-volume columns missing: "
                f"{low_col}, {high_col}"
            )
        #endif

        low = pd.to_numeric(df[low_col], errors="coerce").to_numpy()
        high = pd.to_numeric(df[high_col], errors="coerce").to_numpy()
        width = high - low

        # phi is stored in degrees in these CSVs.  Integrating a differential
        # cross section requires the angular width in radians.
        if low_col == "phimin":
            width = np.radians(width)
        #endif

        volume *= width
    #endfor

    return volume


def load_sector_csv(
    path: Path,
    sector: int,
    cross_section_column: str,
    gtheta_column: str,
) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    #endif

    df = pd.read_csv(path, low_memory=False)

    required = {
        cross_section_column,
        gtheta_column,
        *[name for pair in VOLUME_COLUMNS for name in pair],
    }
    missing = sorted(required.difference(df.columns))

    if missing:
        raise KeyError(
            f"{path} is missing required columns:\n  "
            + "\n  ".join(missing)
        )
    #endif

    values = np.full(len(df), np.nan)
    uncertainties = np.full(len(df), np.nan)

    for irow, cell in enumerate(df[cross_section_column].to_numpy()):
        values[irow], uncertainties[irow] = tuple_value_unc(cell)
    #endfor

    out = pd.DataFrame(
        {
            "sector": sector,
            "g_theta": pd.to_numeric(
                df[gtheta_column],
                errors="coerce",
            ).to_numpy(),
            "sigma": values,
            "sigma_unc": uncertainties,
            "volume": bin_volume(df),
        }
    )

    return out


def integrate_theta_range(
    df: pd.DataFrame,
    theta_low: float,
    theta_high: float,
) -> Tuple[float, float, float, int]:
    """
    Bin-volume-weighted integral over the underlying 4D DVCS bins.

    sigma_int = sum_i sigma_i * DeltaV_i.
    Statistical uncertainties are propagated in quadrature.
    The returned theta is volume*cross-section weighted when possible.
    """
    mask = (
        np.isfinite(df["g_theta"].to_numpy())
        & np.isfinite(df["sigma"].to_numpy())
        & np.isfinite(df["volume"].to_numpy())
        & (df["volume"].to_numpy() > 0.0)
        & (df["g_theta"].to_numpy() >= theta_low)
        & (df["g_theta"].to_numpy() < theta_high)
    )

    if not np.any(mask):
        return np.nan, np.nan, np.nan, 0
    #endif

    sigma = df.loc[mask, "sigma"].to_numpy(dtype=float)
    unc = df.loc[mask, "sigma_unc"].to_numpy(dtype=float)
    volume = df.loc[mask, "volume"].to_numpy(dtype=float)
    theta = df.loc[mask, "g_theta"].to_numpy(dtype=float)

    contribution = sigma * volume
    value = float(np.sum(contribution))

    finite_unc = np.isfinite(unc)
    uncertainty = float(
        np.sqrt(
            np.sum(
                (unc[finite_unc] * volume[finite_unc]) ** 2
            )
        )
    )

    positive = np.isfinite(contribution) & (contribution > 0.0)
    if np.any(positive):
        theta_mean = float(
            np.average(
                theta[positive],
                weights=contribution[positive],
            )
        )
    else:
        theta_mean = float(np.mean(theta))
    #endif

    return value, uncertainty, theta_mean, int(np.sum(mask))


def normalize_sector_values(
    values: np.ndarray,
    uncertainties: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Normalize six sector values to their finite arithmetic mean.

    The uncertainty shown is the numerator uncertainty divided by the same
    mean.  Correlations from using the measured sector mean as denominator are
    intentionally not modeled because this is a closure-shape diagnostic.
    """
    finite = np.isfinite(values) & (values > 0.0)

    normalized = np.full_like(values, np.nan, dtype=float)
    normalized_unc = np.full_like(values, np.nan, dtype=float)

    if not np.any(finite):
        return normalized, normalized_unc
    #endif

    mean = float(np.mean(values[finite]))
    normalized[finite] = values[finite] / mean

    good_unc = finite & np.isfinite(uncertainties)
    normalized_unc[good_unc] = uncertainties[good_unc] / mean

    return normalized, normalized_unc


def build_detailed_theta_projection(
    sectors: Mapping[int, pd.DataFrame],
    theta_edges: np.ndarray,
) -> pd.DataFrame:
    rows = []

    for ibin in range(len(theta_edges) - 1):
        low = float(theta_edges[ibin])
        high = float(theta_edges[ibin + 1])

        values = np.full(6, np.nan)
        uncertainties = np.full(6, np.nan)
        theta_means = np.full(6, np.nan)
        counts = np.zeros(6, dtype=int)

        for sector in range(1, 7):
            (
                values[sector - 1],
                uncertainties[sector - 1],
                theta_means[sector - 1],
                counts[sector - 1],
            ) = integrate_theta_range(
                sectors[sector],
                low,
                high,
            )
        #endfor

        norm, norm_unc = normalize_sector_values(
            values,
            uncertainties,
        )

        finite_theta = np.isfinite(theta_means)
        plot_theta = (
            float(np.mean(theta_means[finite_theta]))
            if np.any(finite_theta)
            else 0.5 * (low + high)
        )

        for sector in range(1, 7):
            rows.append(
                {
                    "theta_bin": ibin,
                    "theta_min": low,
                    "theta_max": high,
                    "theta_mean": plot_theta,
                    "sector": sector,
                    "sigma": values[sector - 1],
                    "sigma_unc": uncertainties[sector - 1],
                    "sector_ratio": norm[sector - 1],
                    "sector_ratio_unc": norm_unc[sector - 1],
                    "n_bins": counts[sector - 1],
                }
            )
        #endfor
    #endfor

    return pd.DataFrame(rows)


def build_matched_theta_cross_section(
    sectors: Mapping[int, pd.DataFrame],
) -> pd.DataFrame:
    rows = []

    for theta_low, theta_high in MATCHED_THETA_BINS:
        values = np.full(6, np.nan)
        uncertainties = np.full(6, np.nan)

        for sector in range(1, 7):
            values[sector - 1], uncertainties[sector - 1], _, _ = (
                integrate_theta_range(
                    sectors[sector],
                    theta_low,
                    theta_high,
                )
            )
        #endfor

        norm, norm_unc = normalize_sector_values(
            values,
            uncertainties,
        )

        for sector in range(1, 7):
            rows.append(
                {
                    "theta_min_deg": theta_low,
                    "theta_max_deg": theta_high,
                    "sector": sector,
                    "cross_section": values[sector - 1],
                    "cross_section_unc": uncertainties[sector - 1],
                    "cross_section_sector_shape": norm[sector - 1],
                    "cross_section_sector_shape_unc": norm_unc[
                        sector - 1
                    ],
                }
            )
        #endfor
    #endfor

    return pd.DataFrame(rows)


def load_efficiency_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    #endif

    df = pd.read_csv(path)

    required = {
        "theta_min_deg",
        "theta_max_deg",
        "sector",
        "Epred_min_GeV",
        "Epred_max_GeV",
        "data_over_aaogen",
        "data_over_aaogen_unc",
        "data_over_clasdis",
        "data_over_clasdis_unc",
    }
    missing = sorted(required.difference(df.columns))

    if missing:
        raise KeyError(
            f"{path} is missing required efficiency columns:\n  "
            + "\n  ".join(missing)
        )
    #endif

    for column in required.difference({"sector"}):
        df[column] = pd.to_numeric(df[column], errors="coerce")
    #endfor
    df["sector"] = pd.to_numeric(
        df["sector"],
        errors="coerce",
    ).astype("Int64")

    return df


def efficiency_sector_shape(
    efficiency: pd.DataFrame,
    theta_low: float,
    theta_high: float,
    energy_low: float,
    energy_high: float,
    mc_name: str,
) -> Tuple[np.ndarray, np.ndarray]:
    ratio_col = f"data_over_{mc_name}"
    unc_col = f"data_over_{mc_name}_unc"

    subset = efficiency[
        np.isclose(efficiency["theta_min_deg"], theta_low)
        & np.isclose(efficiency["theta_max_deg"], theta_high)
        & np.isclose(efficiency["Epred_min_GeV"], energy_low)
        & np.isclose(efficiency["Epred_max_GeV"], energy_high)
    ]

    values = np.full(6, np.nan)
    uncertainties = np.full(6, np.nan)

    for sector in range(1, 7):
        row = subset[subset["sector"] == sector]
        if len(row) != 1:
            continue
        #endif
        values[sector - 1] = float(row.iloc[0][ratio_col])
        uncertainties[sector - 1] = float(row.iloc[0][unc_col])
    #endfor

    return normalize_sector_values(values, uncertainties)


def pearson_finite(
    x: np.ndarray,
    y: np.ndarray,
) -> float:
    finite = np.isfinite(x) & np.isfinite(y)

    if np.sum(finite) < 3:
        return np.nan
    #endif

    x_use = x[finite]
    y_use = y[finite]

    if np.std(x_use) == 0.0 or np.std(y_use) == 0.0:
        return np.nan
    #endif

    return float(np.corrcoef(x_use, y_use)[0, 1])


def plot_detailed_cross_section(
    detailed: pd.DataFrame,
    period: str,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.0, 8.0),
        sharex=True,
    )

    for sector in range(1, 7):
        sub = detailed[detailed["sector"] == sector]

        axes[0].errorbar(
            sub["theta_mean"],
            sub["sigma"],
            yerr=sub["sigma_unc"],
            marker="o",
            linestyle="-",
            capsize=2,
            label=f"S{sector}",
        )
        axes[1].errorbar(
            sub["theta_mean"],
            sub["sector_ratio"],
            yerr=sub["sector_ratio_unc"],
            marker="o",
            linestyle="-",
            capsize=2,
            label=f"S{sector}",
        )
    #endfor

    axes[0].set_ylabel("Integrated cross section (arb. projection units)")
    axes[0].set_title(period)
    axes[0].grid(alpha=0.2)
    axes[0].legend(ncol=3, fontsize=8)

    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1].set_xlabel(r"Photon polar angle $\theta_\gamma$ (deg)")
    axes[1].set_ylabel("Sector / six-sector mean")
    axes[1].grid(alpha=0.2)

    fig.suptitle(
        "DVCS photon-FD-sector dependence",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_matched_comparison(
    cross_section: pd.DataFrame,
    efficiency: pd.DataFrame,
    mc_name: str,
    period: str,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(
        len(MATCHED_THETA_BINS),
        1,
        figsize=(9.0, 11.0),
        squeeze=False,
    )
    sectors = np.arange(1, 7)

    for irow, (theta_low, theta_high) in enumerate(
        MATCHED_THETA_BINS
    ):
        ax = axes[irow, 0]
        xs = cross_section[
            np.isclose(
                cross_section["theta_min_deg"],
                theta_low,
            )
            & np.isclose(
                cross_section["theta_max_deg"],
                theta_high,
            )
        ].sort_values("sector")

        xs_shape = xs[
            "cross_section_sector_shape"
        ].to_numpy(dtype=float)
        xs_unc = xs[
            "cross_section_sector_shape_unc"
        ].to_numpy(dtype=float)

        ax.errorbar(
            sectors,
            xs_shape,
            yerr=xs_unc,
            marker="o",
            linestyle="-",
            capsize=2,
            linewidth=1.7,
            label="DVCS cross section",
        )

        correlation_lines = []

        for energy_low, energy_high in HIGH_ENERGY_BINS:
            eff_shape, eff_unc = efficiency_sector_shape(
                efficiency,
                theta_low,
                theta_high,
                energy_low,
                energy_high,
                mc_name,
            )
            r = pearson_finite(xs_shape, eff_shape)
            correlation_lines.append(
                f"{energy_low:g}-{energy_high:g} GeV: r={r:.2f}"
            )

            ax.errorbar(
                sectors,
                eff_shape,
                yerr=eff_unc,
                marker="s",
                linestyle="--",
                capsize=2,
                label=(
                    f"Efficiency {energy_low:g}-"
                    f"{energy_high:g} GeV"
                ),
            )
        #endfor

        ax.axhline(1.0, linestyle=":", linewidth=1.0)
        ax.set_ylim(0.2, 1.8)
        ax.set_xticks(range(1, 7))
        ax.set_xlabel("Photon FD sector")
        ax.set_ylabel("Quantity / six-sector mean")
        ax.set_title(
            f"{theta_low:g} <= theta_gamma < "
            f"{theta_high:g} deg"
        )
        ax.grid(alpha=0.2)
        ax.text(
            0.02,
            0.04,
            "\n".join(correlation_lines),
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
        )

        if irow == 0:
            ax.legend(fontsize=8, ncol=3)
        #endif
    #endfor

    fig.suptitle(
        f"{period}: DVCS sector shape vs DATA/{mc_name.upper()} "
        "photon-efficiency sector shape",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def write_comparison_table(
    cross_section: pd.DataFrame,
    efficiency: pd.DataFrame,
    output_path: Path,
) -> None:
    rows = []

    for theta_low, theta_high in MATCHED_THETA_BINS:
        xs = cross_section[
            np.isclose(
                cross_section["theta_min_deg"],
                theta_low,
            )
            & np.isclose(
                cross_section["theta_max_deg"],
                theta_high,
            )
        ].sort_values("sector")

        for _, xs_row in xs.iterrows():
            sector = int(xs_row["sector"])

            row = {
                "theta_min_deg": theta_low,
                "theta_max_deg": theta_high,
                "sector": sector,
                "dvcs_cross_section": xs_row["cross_section"],
                "dvcs_cross_section_unc": xs_row[
                    "cross_section_unc"
                ],
                "dvcs_sector_shape": xs_row[
                    "cross_section_sector_shape"
                ],
                "dvcs_sector_shape_unc": xs_row[
                    "cross_section_sector_shape_unc"
                ],
            }

            for energy_low, energy_high in HIGH_ENERGY_BINS:
                for mc_name in ["aaogen", "clasdis"]:
                    shape, shape_unc = efficiency_sector_shape(
                        efficiency,
                        theta_low,
                        theta_high,
                        energy_low,
                        energy_high,
                        mc_name,
                    )
                    key = (
                        f"{mc_name}_{energy_low:g}_to_"
                        f"{energy_high:g}_GeV"
                    ).replace(".", "p")
                    row[key] = shape[sector - 1]
                    row[f"{key}_unc"] = shape_unc[sector - 1]
                #endfor
            #endfor

            rows.append(row)
        #endfor
    #endfor

    pd.DataFrame(rows).to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sector_paths = {
        sector: getattr(args, f"sec{sector}")
        for sector in range(1, 7)
    }

    sectors = {
        sector: load_sector_csv(
            sector_paths[sector],
            sector,
            args.cross_section_column,
            args.gtheta_column,
        )
        for sector in range(1, 7)
    }

    theta_edges = np.linspace(
        args.detailed_theta_min,
        args.detailed_theta_max,
        args.detailed_theta_bins + 1,
    )

    detailed = build_detailed_theta_projection(
        sectors,
        theta_edges,
    )
    matched = build_matched_theta_cross_section(sectors)
    efficiency = load_efficiency_table(args.efficiency_table)

    detailed.to_csv(
        args.output_dir / "dvcs_sector_theta_projection.csv",
        index=False,
    )
    matched.to_csv(
        args.output_dir / "dvcs_sector_matched_theta_bins.csv",
        index=False,
    )

    plot_detailed_cross_section(
        detailed,
        args.period,
        args.output_dir
        / "01_dvcs_photon_sector_theta_projection.png",
    )

    for mc_name in ["aaogen", "clasdis"]:
        plot_matched_comparison(
            matched,
            efficiency,
            mc_name,
            args.period,
            args.output_dir
            / f"02_dvcs_vs_efficiency_{mc_name}.png",
        )
    #endfor

    write_comparison_table(
        matched,
        efficiency,
        args.output_dir
        / "dvcs_efficiency_sector_shape_comparison.csv",
    )

    print("Wrote:")
    print(
        "  "
        + str(
            args.output_dir
            / "01_dvcs_photon_sector_theta_projection.png"
        )
    )
    print(
        "  "
        + str(
            args.output_dir
            / "02_dvcs_vs_efficiency_aaogen.png"
        )
    )
    print(
        "  "
        + str(
            args.output_dir
            / "02_dvcs_vs_efficiency_clasdis.png"
        )
    )
    print(
        "  "
        + str(
            args.output_dir
            / "dvcs_efficiency_sector_shape_comparison.csv"
        )
    )


if __name__ == "__main__":
    main()
