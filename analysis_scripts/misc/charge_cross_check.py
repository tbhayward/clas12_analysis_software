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
    out["hayward_charge"] = (
        pd.to_numeric(df.iloc[:, 2], errors="coerce")
        + pd.to_numeric(df.iloc[:, 3], errors="coerce")
    )

    out = out.dropna(subset=["runnum", "hayward_charge"])
    out["runnum"] = out["runnum"].astype(int)

    return out
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

    merged = merged[(merged["neupane_charge"] != 0) | (merged["hayward_charge"] != 0)].copy()

    if len(merged) == 0:
        raise RuntimeError("Common runs exist, but all common charge values are zero.")
    #endif

    merged["percent_difference"] = 100.0 * (
        merged["hayward_charge"] - merged["neupane_charge"]
    ) / merged["neupane_charge"]

    finite = np.isfinite(merged["percent_difference"].to_numpy())
    percent_diffs = merged.loc[finite, "percent_difference"].to_numpy()

    mean_percent_difference = np.mean(percent_diffs)
    median_percent_difference = np.median(percent_diffs)

    print("")
    print("Charge cross-check: Hayward vs Neupane")
    print("Common runs:", len(merged))
    print("Mean percent difference   = {:.6f}%".format(mean_percent_difference))
    print("Median percent difference = {:.6f}%".format(median_percent_difference))
    print("")
    print("Percent difference definition:")
    print("  100 * (Hayward - Neupane) / Neupane")
    print("")

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(13, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax_top.plot(
        merged["runnum"],
        merged["neupane_charge"],
        marker="o",
        linestyle="none",
        markersize=4,
        label="Neupane",
    )

    ax_top.plot(
        merged["runnum"],
        merged["hayward_charge"],
        marker="s",
        linestyle="none",
        markersize=4,
        label="Hayward",
    )

    ax_top.set_ylabel("Accumulated charge (nC)")
    ax_top.set_title("Run-by-run accumulated charge comparison")
    ax_top.legend()
    ax_top.grid(True, alpha=0.3)

    ax_bottom.axhline(0.0, linewidth=1)
    ax_bottom.plot(
        merged["runnum"],
        merged["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=4,
    )

    ax_bottom.set_xlabel("Run number")
    ax_bottom.set_ylabel("% diff.")
    ax_bottom.grid(True, alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print("Saved:", output_path)
#enddef


if __name__ == "__main__":
    main()
#endif