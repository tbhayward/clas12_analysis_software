#!/usr/bin/env python3
"""
CLAS12 DVCS high-luminosity reach study.

Reads the released RGA unpolarized DVCS cross sections and, optionally,
the published RGA beam-spin asymmetry (BSA) points.  The first-stage study
focuses on the condition |t|/Q^2 < 0.2 and projects statistical precision
for luminosity factors 1x, 2x, 5x, and 10x.

Expected layout:
  higher_level_dvcs_cross_section/
    import/
      clasdb_E214M1.txt
      rga_prl_bsa.txt
    high_luminosity/
      high_luminosity_dvcs_reach.py
"""

from pathlib import Path
import argparse
import math
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LUMI_FACTORS = (1, 2, 5, 10)
RATIO_LINES = (0.1, 0.2, 0.3)
T_EDGES = np.arange(0.10, 1.11, 0.10)
RELSTAT_THRESHOLDS = (0.10, 0.15, 0.20, 0.30)


def parse_args():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--xs",
        type=Path,
        default=here.parent / "import" / "clasdb_E214M1.txt",
    )
    parser.add_argument(
        "--bsa",
        type=Path,
        default=here.parent / "import" / "rga_prl_bsa.txt",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=here / "output",
    )
    parser.add_argument("--clean-ratio", type=float, default=0.20)
    parser.add_argument(
        "--include-bsa",
        action="store_true",
        help="Also make a BSA-vs-XS kinematic-coverage comparison.",
    )
    return parser.parse_args()


def load_xs(path):
    rows = []

    with path.open() as f:
        for line in f:
            fields = line.split()
            if len(fields) < 8 or re.fullmatch(r"\d+", fields[0]) is None:
                continue
            #endif

            rows.append(
                [
                    int(fields[0]),
                    float(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                    float(fields[4]),
                    float(fields[5]),
                    float(fields[6]),
                    float(fields[7]),
                ]
            )
        #endfor

    df = pd.DataFrame(
        rows,
        columns=[
            "bin",
            "xB",
            "Q2",
            "t",
            "phi",
            "sigma",
            "stat",
            "syst",
        ],
    )

    if len(df) == 0:
        raise RuntimeError(f"No cross-section rows parsed from {path}")
    #endif

    df["abs_t"] = np.abs(df["t"])
    df["t_over_Q2"] = df["abs_t"] / df["Q2"]
    df["rel_stat"] = df["stat"] / np.abs(df["sigma"])

    # Statistical-power proxy only. This is NOT the raw selected event count.
    df["N_eff"] = np.where(
        np.isfinite(df["rel_stat"]) & (df["rel_stat"] > 0),
        1.0 / df["rel_stat"]**2,
        np.nan,
    )

    return df


def load_bsa(path):
    df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=["phi", "Q2", "xB", "t", "Eb", "A", "sigA"],
    )
    df["abs_t"] = np.abs(df["t"])
    df["t_over_Q2"] = df["abs_t"] / df["Q2"]
    return df


def xs_cells(df):
    """
    The published XS 'bin' index identifies an (xB,Q2,t) cell; phi points
    within the same bin share that cell.
    """
    return (
        df.groupby("bin", as_index=False)
        .agg(
            xB=("xB", "median"),
            Q2=("Q2", "median"),
            abs_t=("abs_t", "median"),
            t_over_Q2=("t_over_Q2", "median"),
            n_phi=("phi", "size"),
            median_rel_stat=("rel_stat", "median"),
            sum_N_eff=("N_eff", "sum"),
        )
    )


def plot_q2_t_map(xs, clean_ratio, outdir):
    cells = xs_cells(xs)
    clean = cells["t_over_Q2"] < clean_ratio

    fig, ax = plt.subplots(figsize=(9, 7))

    ax.scatter(
        cells.loc[~clean, "abs_t"],
        cells.loc[~clean, "Q2"],
        s=24,
        alpha=0.40,
        label=rf"$|t|/Q^2 \geq {clean_ratio:.1f}$",
    )
    ax.scatter(
        cells.loc[clean, "abs_t"],
        cells.loc[clean, "Q2"],
        s=30,
        alpha=0.85,
        label=rf"$|t|/Q^2 < {clean_ratio:.1f}$",
    )

    tmax = max(1.05, 1.03 * cells["abs_t"].max())
    tgrid = np.linspace(0.0, tmax, 500)

    for ratio in RATIO_LINES:
        ax.plot(
            tgrid,
            tgrid / ratio,
            linewidth=1.5,
            label=rf"$|t|/Q^2={ratio:.1f}$",
        )
    #endfor

    ax.set_xlim(0.0, tmax)
    ax.set_ylim(0.8, max(6.0, 1.05 * cells["Q2"].max()))
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title("CLAS12 RGA DVCS cross-section kinematic reach")
    ax.grid(alpha=0.2)
    ax.legend(ncol=2, fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "01_q2_vs_t_ratio_boundaries.png", dpi=250)
    plt.close(fig)


def plot_clean_stat_map(xs, clean_ratio, outdir):
    clean = xs[xs["t_over_Q2"] < clean_ratio].copy()

    fig, ax = plt.subplots(figsize=(9, 7))
    sc = ax.scatter(
        clean["abs_t"],
        clean["Q2"],
        c=100.0 * clean["rel_stat"],
        s=22,
        alpha=0.75,
    )
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("Relative statistical uncertainty (%)")

    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title(rf"4D points satisfying $|t|/Q^2<{clean_ratio:.1f}$")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "02_clean_region_statistical_precision.png", dpi=250)
    plt.close(fig)


def make_t_summary(xs, clean_ratio):
    clean = xs[xs["t_over_Q2"] < clean_ratio].copy()
    clean["t_interval"] = pd.cut(
        clean["abs_t"],
        bins=T_EDGES,
        right=False,
        include_lowest=True,
    )

    rows = []

    for interval, g in clean.groupby("t_interval", observed=False):
        if len(g) == 0:
            continue
        #endif

        cells = xs_cells(g)

        row = {
            "t_low": interval.left,
            "t_high": interval.right,
            "n_4d_points": len(g),
            "n_xB_Q2_t_cells": len(cells),
            "median_phi_points_per_cell": cells["n_phi"].median(),
            "sum_N_eff_1x": g["N_eff"].sum(),
        }

        for factor in LUMI_FACTORS:
            rel = g["rel_stat"] / math.sqrt(factor)
            row[f"median_rel_stat_{factor}x"] = np.nanmedian(rel)
            row[f"p90_rel_stat_{factor}x"] = np.nanpercentile(rel, 90)
            row[f"sum_N_eff_{factor}x"] = factor * g["N_eff"].sum()
        #endfor

        rows.append(row)
    #endfor

    return pd.DataFrame(rows)


def plot_t_summary(summary, outdir):
    centers = 0.5 * (summary["t_low"] + summary["t_high"])

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.bar(
        centers,
        summary["n_4d_points"],
        width=0.085,
        alpha=0.8,
    )
    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Number of 4D cross-section points")
    ax.set_title(r"Measured points satisfying $|t|/Q^2<0.2$")
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "03_clean_points_vs_t.png", dpi=250)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 6))

    for factor in LUMI_FACTORS:
        ax.plot(
            centers,
            100.0 * summary[f"median_rel_stat_{factor}x"],
            marker="o",
            linewidth=1.8,
            label=f"{factor}x luminosity",
        )
    #endfor

    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel("Median relative statistical uncertainty (%)")
    ax.set_title(r"Projected precision for $|t|/Q^2<0.2$")
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "04_projected_precision_vs_t.png", dpi=250)
    plt.close(fig)


def make_usability_summary(xs, clean_ratio):
    """
    Exploratory threshold study.

    This deliberately does NOT declare a CFF fit "possible." It asks how
    many measured 4D points and (xB,Q2,t) cells reach various statistical
    precision thresholds as luminosity increases.
    """
    clean = xs[xs["t_over_Q2"] < clean_ratio].copy()
    rows = []

    for threshold in RELSTAT_THRESHOLDS:
        for factor in LUMI_FACTORS:
            temp = clean.copy()
            temp["projected_rel_stat"] = temp["rel_stat"] / math.sqrt(factor)
            temp["usable_point"] = temp["projected_rel_stat"] <= threshold

            per_cell = (
                temp.groupby("bin", as_index=False)
                .agg(
                    xB=("xB", "median"),
                    Q2=("Q2", "median"),
                    abs_t=("abs_t", "median"),
                    n_phi_total=("phi", "size"),
                    n_phi_usable=("usable_point", "sum"),
                )
            )

            for min_phi in (4, 6, 8):
                usable_cells = per_cell["n_phi_usable"] >= min_phi

                rows.append(
                    {
                        "rel_stat_threshold": threshold,
                        "luminosity_factor": factor,
                        "min_usable_phi_points": min_phi,
                        "n_clean_4d_points": len(temp),
                        "n_usable_4d_points": int(temp["usable_point"].sum()),
                        "n_clean_cells": len(per_cell),
                        "n_usable_cells": int(usable_cells.sum()),
                        "max_abs_t_usable": (
                            per_cell.loc[usable_cells, "abs_t"].max()
                            if usable_cells.any()
                            else np.nan
                        ),
                        "p90_abs_t_usable": (
                            np.percentile(
                                per_cell.loc[usable_cells, "abs_t"],
                                90,
                            )
                            if usable_cells.any()
                            else np.nan
                        ),
                    }
                )
            #endfor
        #endfor
    #endfor

    return pd.DataFrame(rows)


def plot_usability(summary, outdir):
    # Main workshop diagnostic: require six phi points and show several
    # point-level statistical precision requirements.
    s = summary[summary["min_usable_phi_points"] == 6]

    fig, ax = plt.subplots(figsize=(9, 6))

    for threshold in RELSTAT_THRESHOLDS:
        g = s[np.isclose(s["rel_stat_threshold"], threshold)]
        ax.plot(
            g["luminosity_factor"],
            g["n_usable_cells"],
            marker="o",
            linewidth=1.8,
            label=f"{100*threshold:.0f}% per phi point",
        )
    #endfor

    ax.set_xlabel("Luminosity relative to current RGA dataset")
    ax.set_ylabel(r"Usable $(x_B,Q^2,t)$ cells")
    ax.set_title(
        r"Statistically usable cells inside $|t|/Q^2<0.2$"
        "\n(require at least six phi points)"
    )
    ax.set_xticks(LUMI_FACTORS)
    ax.grid(alpha=0.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "05_usable_cells_vs_luminosity.png", dpi=250)
    plt.close(fig)


def plot_bsa_xs_coverage(xs, bsa, clean_ratio, outdir):
    """
    Optional diagnostic only. BSA points do not carry the same bin IDs as
    the cross-section table, so this is a coverage comparison rather than
    an attempted point-by-point matching.
    """
    xsc = xs_cells(xs)
    xs_clean = xsc["t_over_Q2"] < clean_ratio
    bsa_clean = bsa["t_over_Q2"] < clean_ratio

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.scatter(
        bsa.loc[~bsa_clean, "abs_t"],
        bsa.loc[~bsa_clean, "Q2"],
        s=10,
        alpha=0.15,
        label="BSA points outside preferred region",
    )
    ax.scatter(
        bsa.loc[bsa_clean, "abs_t"],
        bsa.loc[bsa_clean, "Q2"],
        s=12,
        alpha=0.35,
        label="BSA points inside preferred region",
    )
    ax.scatter(
        xsc.loc[xs_clean, "abs_t"],
        xsc.loc[xs_clean, "Q2"],
        s=35,
        marker="x",
        label="XS cells inside preferred region",
    )

    tmax = max(xsc["abs_t"].max(), bsa["abs_t"].max())
    tgrid = np.linspace(0.0, 1.03 * tmax, 500)
    ax.plot(
        tgrid,
        tgrid / clean_ratio,
        linewidth=1.8,
        label=rf"$|t|/Q^2={clean_ratio:.1f}$",
    )

    ax.set_xlabel(r"$|t|$ (GeV$^2$)")
    ax.set_ylabel(r"$Q^2$ (GeV$^2$)")
    ax.set_title("Published CLAS12 RGA BSA and cross-section coverage")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "06_bsa_xs_coverage.png", dpi=250)
    plt.close(fig)


def print_summary(xs, clean_ratio):
    cells = xs_cells(xs)
    clean = xs[xs["t_over_Q2"] < clean_ratio]
    clean_cells = xs_cells(clean)

    print()
    print("CLAS12 DVCS high-luminosity reach")
    print("---------------------------------")
    print(f"Cross-section 4D points: {len(xs)}")
    print(f"Unique (xB,Q2,t) cells: {len(cells)}")
    print(
        f"|t|/Q2 < {clean_ratio:.2f}: "
        f"{len(clean)} / {len(xs)} 4D points "
        f"({100*len(clean)/len(xs):.1f}%)"
    )
    print(
        f"|t|/Q2 < {clean_ratio:.2f}: "
        f"{len(clean_cells)} / {len(cells)} cells "
        f"({100*len(clean_cells)/len(cells):.1f}%)"
    )
    print(
        "Median relative statistical uncertainty in preferred region: "
        f"{100*np.nanmedian(clean['rel_stat']):.2f}%"
    )
    print(
        "90th-percentile relative statistical uncertainty: "
        f"{100*np.nanpercentile(clean['rel_stat'], 90):.2f}%"
    )
    print(
        "Largest measured |t| satisfying the ratio cut: "
        f"{clean_cells['abs_t'].max():.3f} GeV^2"
    )
    print()
    print(
        "N_eff=(sigma/stat)^2 is a statistical-power proxy, "
        "not the raw event count."
    )
    print()


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    xs = load_xs(args.xs)

    print_summary(xs, args.clean_ratio)

    # Preserve the full point-level table for the eventual PARTONS/CFF study.
    for factor in LUMI_FACTORS:
        xs[f"stat_{factor}x"] = xs["stat"] / math.sqrt(factor)
        xs[f"rel_stat_{factor}x"] = xs["rel_stat"] / math.sqrt(factor)
        xs[f"N_eff_{factor}x"] = factor * xs["N_eff"]
    #endfor

    xs["preferred_ratio_region"] = xs["t_over_Q2"] < args.clean_ratio
    xs.to_csv(args.output / "xs_pointwise_luminosity_projection.csv", index=False)

    cells = xs_cells(xs)
    cells.to_csv(args.output / "xs_cell_kinematics.csv", index=False)

    t_summary = make_t_summary(xs, args.clean_ratio)
    t_summary.to_csv(args.output / "xs_clean_t_summary.csv", index=False)

    usability = make_usability_summary(xs, args.clean_ratio)
    usability.to_csv(args.output / "xs_usability_vs_luminosity.csv", index=False)

    plot_q2_t_map(xs, args.clean_ratio, args.output)
    plot_clean_stat_map(xs, args.clean_ratio, args.output)
    plot_t_summary(t_summary, args.output)
    plot_usability(usability, args.output)

    if args.include_bsa:
        bsa = load_bsa(args.bsa)
        bsa["preferred_ratio_region"] = bsa["t_over_Q2"] < args.clean_ratio
        bsa.to_csv(args.output / "bsa_pointwise_ratio_study.csv", index=False)
        plot_bsa_xs_coverage(xs, bsa, args.clean_ratio, args.output)

        print(f"BSA 4D points read: {len(bsa)}")
        print(
            f"BSA points with |t|/Q2 < {args.clean_ratio:.2f}: "
            f"{(bsa['t_over_Q2'] < args.clean_ratio).sum()}"
        )
    #endif

    print(f"Outputs written to {args.output}")


if __name__ == "__main__":
    main()
