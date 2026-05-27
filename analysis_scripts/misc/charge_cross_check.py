#!/usr/bin/env python3

import os
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
    out["hayward_charge_col2"] = pd.to_numeric(df.iloc[:, 1], errors="coerce")
    out["hayward_charge_col3_col4"] = (
        pd.to_numeric(df.iloc[:, 2], errors="coerce")
        + pd.to_numeric(df.iloc[:, 3], errors="coerce")
    )

    out = out.dropna(subset=["runnum", "hayward_charge_col2", "hayward_charge_col3_col4"])
    out["runnum"] = out["runnum"].astype(int)

    return out
#enddef


def print_summary_and_large_differences(merged, charge_column, label):
    local = merged.copy()

    local = local[local["neupane_charge"] != 0].copy()

    if len(local) == 0:
        raise RuntimeError("No common runs with nonzero Neupane charge for {}.".format(label))
    #endif

    local["percent_difference"] = 100.0 * (
        local[charge_column] - local["neupane_charge"]
    ) / local["neupane_charge"]

    finite = np.isfinite(local["percent_difference"].to_numpy())
    percent_diffs = local.loc[finite, "percent_difference"].to_numpy()

    mean_percent_difference = np.mean(percent_diffs)
    median_percent_difference = np.median(percent_diffs)

    print("")
    print("Charge cross-check: Hayward ({}) vs Neupane RUN::Scaler".format(label))
    print("Common runs:", len(local))
    print("Mean percent difference   = {:.6f}%".format(mean_percent_difference))
    print("Median percent difference = {:.6f}%".format(median_percent_difference))
    print("Percent difference definition:")
    print("  100 * (Hayward - Neupane) / Neupane")
    print("")

    large = local[np.abs(local["percent_difference"]) > 1.0].copy()

    if len(large) == 0:
        print("No runs differ by more than 1 percent for {}.".format(label))
    else:
        print("Runs differing by more than 1 percent for {}:".format(label))
        print("runnum, neupane_charge, hayward_charge, percent_difference")
        for _, row in large.iterrows():
            print(
                "{}, {:.6f}, {:.6f}, {:.6f}%".format(
                    int(row["runnum"]),
                    row["neupane_charge"],
                    row[charge_column],
                    row["percent_difference"],
                )
            )
        #endfor
    #endif

    return local
#enddef


def plot_panel(ax_top, ax_bottom, local, charge_column, title):
    ax_top.plot(
        local["runnum"],
        local["neupane_charge"],
        marker="o",
        linestyle="none",
        markersize=4,
        label="Neupane, RUN::Scaler",
    )

    ax_top.plot(
        local["runnum"],
        local[charge_column],
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
        local["runnum"],
        local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=4,
    )

    ax_bottom.set_xlabel("Run number")
    ax_bottom.set_ylabel("% diff.")
    ax_bottom.grid(True, alpha=0.3)
#enddef


def main():
    neupane_path = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/krishnas_charges.csv"
    hayward_path = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity/global.csv"

    output_path = "output/charge_cross_check.png"
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    neupane = read_neupane_csv(neupane_path)
    hayward = read_hayward_csv(hayward_path)

    merged = pd.merge(neupane, hayward, on="runnum", how="inner")
    merged = merged.sort_values("runnum").reset_index(drop=True)

    if len(merged) == 0:
        raise RuntimeError("No common run numbers found between Neupane and Hayward CSV files.")
    #endif

    left = print_summary_and_large_differences(
        merged,
        "hayward_charge_col2",
        "RUN::Scaler",
    )

    right = print_summary_and_large_differences(
        merged,
        "hayward_charge_col3_col4",
        "HEL::Scaler",
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

    plot_panel(
        ax_left_top,
        ax_left_bottom,
        left,
        "hayward_charge_col2",
        "Hayward, RUN::Scaler",
    )

    plot_panel(
        ax_right_top,
        ax_right_bottom,
        right,
        "hayward_charge_col3_col4",
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


if __name__ == "__main__":
    main()
#endif