#!/usr/bin/env python3

import os
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def read_neupane_csv(path):
    df = pd.read_csv(path)

    out = pd.DataFrame()
    out["runnum"] = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    out["neupane_charge"] = pd.to_numeric(df.iloc[:, 1], errors="coerce")

    out = out.dropna(subset=["runnum", "neupane_charge"])
    out["runnum"] = out["runnum"].astype(int)

    return out
#enddef


def read_hayward_csv(path):
    df = pd.read_csv(path, header=None)

    out = pd.DataFrame()
    out["runnum"] = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    out["hayward_run_scaler"] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    out["hayward_hel_scaler"] = (
        pd.to_numeric(df.iloc[:, 2], errors="coerce")
        + pd.to_numeric(df.iloc[:, 3], errors="coerce")
    )

    out = out.dropna(subset=["runnum", "hayward_run_scaler", "hayward_hel_scaler"])
    out["runnum"] = out["runnum"].astype(int)

    return out
#enddef


def read_run_period_charge_file(path):
    df = pd.read_csv(path, header=None)

    out = pd.DataFrame()
    out["runnum"] = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    out["run_scaler"] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    out["hel_scaler"] = (
        pd.to_numeric(df.iloc[:, 2], errors="coerce")
        + pd.to_numeric(df.iloc[:, 3], errors="coerce")
    )

    out = out.dropna(subset=["runnum", "run_scaler", "hel_scaler"])
    out["runnum"] = out["runnum"].astype(int)

    return out
#enddef


def compute_percent_difference(df, numerator_column, denominator_column, output_column):
    local = df.copy()

    local = local[
        (local[numerator_column] != 0.0)
        & (local[denominator_column] != 0.0)
    ].copy()

    local[output_column] = 100.0 * (
        local[numerator_column] - local[denominator_column]
    ) / local[denominator_column]

    local = local[np.isfinite(local[output_column].to_numpy())].copy()

    return local
#enddef


def print_neupane_hayward_summary(local, hayward_label, print_large_differences):
    percent_diffs = local["percent_difference"].to_numpy()

    mean_percent_difference = np.mean(percent_diffs)
    median_percent_difference = np.median(percent_diffs)

    print("")
    print("Charge cross-check: Hayward {} vs Neupane RUN::Scaler".format(hayward_label))
    print("Common nonzero-charge runs used in comparison:", len(local))
    print("Mean percent difference   = {:.6f}%".format(mean_percent_difference))
    print("Median percent difference = {:.6f}%".format(median_percent_difference))
    print("Percent difference definition:")
    print("  100 * (Hayward - Neupane) / Neupane")
    print("")

    if print_large_differences:
        large = local[np.abs(local["percent_difference"]) > 1.0].copy()

        if len(large) == 0:
            print("No runs differ by more than 1 percent for Hayward {}.".format(hayward_label))
        else:
            print("Runs differing by more than 1 percent for Hayward {}:".format(hayward_label))
            print("runnum, neupane_charge, hayward_charge, percent_difference")
            for _, row in large.iterrows():
                print(
                    "{}, {:.6f}, {:.6f}, {:.6f}%".format(
                        int(row["runnum"]),
                        row["neupane_charge"],
                        row["hayward_charge"],
                        row["percent_difference"],
                    )
                )
            #endfor
        #endif
    #endif
#enddef


def print_run_period_summary(label, local, outlier_threshold_percent):
    percent_diffs = local["percent_difference"].to_numpy()

    mean_percent_difference = np.mean(percent_diffs)
    rms_percent_difference = math.sqrt(np.mean(percent_diffs * percent_diffs))

    print("")
    print("RUN::Scaler vs HEL::Scaler charge cross-check: {}".format(label))
    print("Nonzero-charge runs used in comparison:", len(local))
    print("Mean percent difference = {:.6f}%".format(mean_percent_difference))
    print("RMS percent difference  = {:.6f}%".format(rms_percent_difference))
    print("Percent difference definition:")
    print("  100 * (HEL::Scaler - RUN::Scaler) / RUN::Scaler")

    outliers = local[np.abs(local["percent_difference"]) > outlier_threshold_percent].copy()

    print("")
    print("Outlier threshold for {}: |percent difference| > {:.3f}%".format(label, outlier_threshold_percent))

    if len(outliers) == 0:
        print("No outlier runs found for {}.".format(label))
    else:
        print("Outlier runs for {}:".format(label))
        print("runnum, RUN::Scaler, HEL::Scaler, percent_difference")
        for _, row in outliers.iterrows():
            print(
                "{}, {:.6f}, {:.6f}, {:.6f}%".format(
                    int(row["runnum"]),
                    row["run_scaler"],
                    row["hel_scaler"],
                    row["percent_difference"],
                )
            )
        #endfor
    #endif
#enddef


def plot_neupane_hayward_panel(ax_top, ax_bottom, merged_all, comparison_local, hayward_column, title):
    ax_top.plot(
        merged_all["runnum"],
        merged_all["neupane_charge"],
        marker="o",
        linestyle="none",
        markersize=4,
        label="Neupane, RUN::Scaler",
    )

    ax_top.plot(
        merged_all["runnum"],
        merged_all[hayward_column],
        marker="s",
        linestyle="none",
        markersize=4,
        label=title,
    )

    ax_top.set_ylabel("Accumulated charge (nC)")
    ax_top.set_title(title)
    ax_top.legend()
    ax_top.grid(True, alpha=0.3)

    ax_bottom.axhline(0.0, linewidth=1)
    ax_bottom.plot(
        comparison_local["runnum"],
        comparison_local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=4,
    )

    ax_bottom.set_xlabel("Run number")
    ax_bottom.set_ylabel("% diff.")
    ax_bottom.grid(True, alpha=0.3)
#enddef


def plot_run_period_panel(ax_top, ax_mid, ax_bottom, comparison_local, title):
    ax_top.plot(
        comparison_local["runnum"],
        comparison_local["run_scaler"],
        marker="o",
        linestyle="none",
        markersize=4,
        label="RUN::Scaler",
    )

    ax_top.plot(
        comparison_local["runnum"],
        comparison_local["hel_scaler"],
        marker="s",
        linestyle="none",
        markersize=4,
        label="HEL::Scaler",
    )

    ax_top.set_ylabel("Charge (nC)")
    ax_top.set_title(title)
    ax_top.legend(fontsize=8)
    ax_top.grid(True, alpha=0.3)

    ax_mid.axhline(0.0, linewidth=1)
    ax_mid.axhline(10.0, linewidth=1, linestyle="--")
    ax_mid.axhline(-10.0, linewidth=1, linestyle="--")
    ax_mid.plot(
        comparison_local["runnum"],
        comparison_local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=3,
    )

    ax_mid.set_ylabel("% diff.")
    ax_mid.grid(True, alpha=0.3)

    ax_bottom.axhline(0.0, linewidth=1)
    ax_bottom.plot(
        comparison_local["runnum"],
        comparison_local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=3,
    )

    ax_bottom.set_ylim(-10.0, 10.0)
    ax_bottom.set_xlabel("Run number")
    ax_bottom.set_ylabel("% diff.")
    ax_bottom.grid(True, alpha=0.3)
#enddef


def make_neupane_hayward_comparison_plot(neupane_path, hayward_path, output_path):
    neupane = read_neupane_csv(neupane_path)
    hayward = read_hayward_csv(hayward_path)

    merged = pd.merge(neupane, hayward, on="runnum", how="inner")
    merged = merged.sort_values("runnum").reset_index(drop=True)

    if len(merged) == 0:
        raise RuntimeError("No common run numbers found between Neupane and Hayward CSV files.")
    #endif

    left_base = pd.DataFrame()
    left_base["runnum"] = merged["runnum"]
    left_base["neupane_charge"] = merged["neupane_charge"]
    left_base["hayward_charge"] = merged["hayward_run_scaler"]

    right_base = pd.DataFrame()
    right_base["runnum"] = merged["runnum"]
    right_base["neupane_charge"] = merged["neupane_charge"]
    right_base["hayward_charge"] = merged["hayward_hel_scaler"]

    left = compute_percent_difference(
        left_base,
        "hayward_charge",
        "neupane_charge",
        "percent_difference",
    )

    right = compute_percent_difference(
        right_base,
        "hayward_charge",
        "neupane_charge",
        "percent_difference",
    )

    print_neupane_hayward_summary(
        left,
        "RUN::Scaler",
        print_large_differences=True,
    )

    print_neupane_hayward_summary(
        right,
        "HEL::Scaler",
        print_large_differences=False,
    )

    fig = plt.figure(
        figsize=(18, 8),
        constrained_layout=True,
    )

    outer = fig.add_gridspec(
        1,
        2,
        wspace=0.12,
    )

    left_grid = outer[0].subgridspec(
        2,
        1,
        height_ratios=[3, 1],
        hspace=0.05,
    )

    right_grid = outer[1].subgridspec(
        2,
        1,
        height_ratios=[3, 1],
        hspace=0.05,
    )

    ax_left_top = fig.add_subplot(left_grid[0])
    ax_left_bottom = fig.add_subplot(left_grid[1], sharex=ax_left_top)

    ax_right_top = fig.add_subplot(right_grid[0])
    ax_right_bottom = fig.add_subplot(right_grid[1], sharex=ax_right_top)

    plot_neupane_hayward_panel(
        ax_left_top,
        ax_left_bottom,
        merged,
        left,
        "hayward_run_scaler",
        "Hayward, RUN::Scaler",
    )

    plot_neupane_hayward_panel(
        ax_right_top,
        ax_right_bottom,
        merged,
        right,
        "hayward_hel_scaler",
        "Hayward, HEL::Scaler",
    )

    plt.setp(ax_left_top.get_xticklabels(), visible=False)
    plt.setp(ax_right_top.get_xticklabels(), visible=False)

    fig.suptitle("Run-by-run accumulated charge comparison", fontsize=16)

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print("")
    print("Saved:", output_path)
#enddef


def make_run_period_scaler_comparison_plot(run_period_files, output_path):
    outlier_threshold_percent = 10.0

    fig = plt.figure(
        figsize=(18, 8),
        constrained_layout=True,
    )

    outer = fig.add_gridspec(
        1,
        3,
        wspace=0.18,
    )

    for index, item in enumerate(run_period_files):
        label = item["label"]
        path = item["path"]

        local = read_run_period_charge_file(path)
        local = local.sort_values("runnum").reset_index(drop=True)

        comparison_local = compute_percent_difference(
            local,
            "hel_scaler",
            "run_scaler",
            "percent_difference",
        )

        if len(comparison_local) == 0:
            raise RuntimeError("No nonzero RUN::Scaler and HEL::Scaler entries found for {}.".format(label))
        #endif

        print_run_period_summary(
            label,
            comparison_local,
            outlier_threshold_percent,
        )

        inner = outer[index].subgridspec(
            3,
            1,
            height_ratios=[3, 1, 1],
            hspace=0.05,
        )

        ax_top = fig.add_subplot(inner[0])
        ax_mid = fig.add_subplot(inner[1], sharex=ax_top)
        ax_bottom = fig.add_subplot(inner[2], sharex=ax_top)

        plot_run_period_panel(
            ax_top,
            ax_mid,
            ax_bottom,
            comparison_local,
            label,
        )

        plt.setp(ax_top.get_xticklabels(), visible=False)
        plt.setp(ax_mid.get_xticklabels(), visible=False)
    #endfor

    fig.suptitle("RUN::Scaler vs HEL::Scaler accumulated charge comparison", fontsize=16)

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print("")
    print("Saved:", output_path)
#enddef


def main():
    base_dir = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity"

    neupane_path = os.path.join(base_dir, "krishnas_charges.csv")
    hayward_path = os.path.join(base_dir, "global.csv")

    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    make_neupane_hayward_comparison_plot(
        neupane_path,
        hayward_path,
        os.path.join(output_dir, "charge_cross_check.png"),
    )

    run_period_files = [
        {
            "label": "Fa18 Inb",
            "path": os.path.join(base_dir, "rga_fa18_inb.txt"),
        },
        {
            "label": "Fa18 Out",
            "path": os.path.join(base_dir, "rga_fa18_out.txt"),
        },
        {
            "label": "Sp19 Inb",
            "path": os.path.join(base_dir, "rga_sp19_inb.txt"),
        },
    ]

    make_run_period_scaler_comparison_plot(
        run_period_files,
        os.path.join(output_dir, "run_period_scaler_cross_check.png"),
    )
#enddef


if __name__ == "__main__":
    main()
#endif