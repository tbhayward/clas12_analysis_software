#!/usr/bin/env python3

import os
import csv
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def fatal(message):
    raise RuntimeError(message)
#enddef


def read_neupane_csv(path):
    if not os.path.exists(path):
        fatal("Neupane CSV does not exist: {}".format(path))
    #endif

    df = pd.read_csv(path)

    if df.shape[1] < 2:
        fatal("Neupane CSV must have at least 2 columns: {}".format(path))
    #endif

    out = pd.DataFrame()
    out["runnum"] = pd.to_numeric(df.iloc[:, 0], errors="coerce")
    out["neupane_charge"] = pd.to_numeric(df.iloc[:, 1], errors="coerce")

    out = out.dropna(subset=["runnum", "neupane_charge"])
    out["runnum"] = out["runnum"].astype(int)

    return out
#enddef


def read_hayward_csv(path):
    if not os.path.exists(path):
        fatal("Hayward CSV does not exist: {}".format(path))
    #endif

    df = pd.read_csv(path, header=None)

    if df.shape[1] < 4:
        fatal("Hayward CSV must have at least 4 columns: {}".format(path))
    #endif

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


def read_flat_charge_file(path):
    if not os.path.exists(path):
        fatal("Charge file does not exist: {}".format(path))
    #endif

    df = pd.read_csv(path, header=None, comment="#")

    if df.shape[1] < 4:
        fatal("Flat charge file must have at least 4 columns: {}".format(path))
    #endif

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


def read_sectioned_charge_file(path, section_header):
    if not os.path.exists(path):
        fatal("Sectioned charge file does not exist: {}".format(path))
    #endif

    rows = []
    active = False
    found_header = False

    with open(path, "r") as f:
        for raw_line in f:
            line = raw_line.strip()

            if line == "":
                continue
            #endif

            if line.startswith("#"):
                if line == section_header:
                    active = True
                    found_header = True
                    continue
                #endif

                if active:
                    break
                #endif

                continue
            #endif

            if not active:
                continue
            #endif

            fields = next(csv.reader([line]))

            if len(fields) < 4:
                fatal(
                    "Malformed numeric row in section '{}'. Expected at least 4 columns, got {}: {}".format(
                        section_header,
                        len(fields),
                        line,
                    )
                )
            #endif

            try:
                runnum = int(float(fields[0]))
                run_scaler = float(fields[1])
                hel_scaler = float(fields[2]) + float(fields[3])
            except Exception as exc:
                fatal(
                    "Could not parse numeric row in section '{}': {}\nError: {}".format(
                        section_header,
                        line,
                        exc,
                    )
                )
            #endtry

            rows.append(
                {
                    "runnum": runnum,
                    "run_scaler": run_scaler,
                    "hel_scaler": hel_scaler,
                }
            )
        #endfor
    #endwith

    if not found_header:
        fatal("Did not find section header '{}' in {}".format(section_header, path))
    #endif

    if len(rows) == 0:
        fatal("Section '{}' in {} contained no numeric rows.".format(section_header, path))
    #endif

    out = pd.DataFrame(rows)
    out = out.dropna(subset=["runnum", "run_scaler", "hel_scaler"])
    out["runnum"] = out["runnum"].astype(int)

    return out
#enddef


def read_charge_source(source):
    if source["kind"] == "flat":
        return read_flat_charge_file(source["path"])
    #endif

    if source["kind"] == "sectioned":
        return read_sectioned_charge_file(source["path"], source["section_header"])
    #endif

    if source["kind"] == "combined_sectioned":
        frames = []

        for section_header in source["section_headers"]:
            frames.append(read_sectioned_charge_file(source["path"], section_header))
        #endfor

        out = pd.concat(frames, ignore_index=True)
        out = out.sort_values("runnum").reset_index(drop=True)

        duplicate_runnums = out[out.duplicated(subset=["runnum"], keep=False)]["runnum"].unique()

        if len(duplicate_runnums) > 0:
            duplicate_string = ", ".join(str(int(x)) for x in duplicate_runnums)
            fatal(
                "Duplicate run numbers found while combining sections for '{}': {}".format(
                    source["label"],
                    duplicate_string,
                )
            )
        #endif

        return out
    #endif

    fatal("Unknown source kind '{}' for label '{}'.".format(source["kind"], source["label"]))
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


def split_outliers(local, outlier_threshold_percent):
    non_outliers = local[np.abs(local["percent_difference"]) <= outlier_threshold_percent].copy()
    outliers = local[np.abs(local["percent_difference"]) > outlier_threshold_percent].copy()

    return non_outliers, outliers
#enddef


def compute_mean_and_std(values):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    if len(arr) == 0:
        return float("nan"), float("nan")
    #endif

    mean_value = float(np.mean(arr))
    std_value = float(math.sqrt(np.mean((arr - mean_value) * (arr - mean_value))))

    return mean_value, std_value
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


def print_run_period_summary(label, local, outlier_threshold_percent, print_outliers):
    non_outliers, outliers = split_outliers(local, outlier_threshold_percent)

    clean_percent_diffs = non_outliers["percent_difference"].to_numpy()
    clean_mean, clean_std = compute_mean_and_std(clean_percent_diffs)

    print("")
    print("RUN::Scaler vs HEL::Scaler charge cross-check: {}".format(label))
    print("Nonzero-charge runs found:", len(local))
    print("Outlier definition: |percent difference| > {:.3f}%".format(outlier_threshold_percent))
    print("Outlier runs excluded from mean/std:", len(outliers))
    print("Non-outlier runs used in mean/std:", len(non_outliers))
    print("Mean percent difference = {:.6f}%".format(clean_mean))
    print("Std. dev. percent difference = {:.6f}%".format(clean_std))
    print("Percent difference definition:")
    print("  100 * (HEL::Scaler - RUN::Scaler) / RUN::Scaler")

    if print_outliers:
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
    else:
        if len(outliers) > 0:
            print("")
            print(
                "Note: {} has {} runs beyond |percent difference| > {:.3f}%, but detailed outlier printing is disabled for this panel.".format(
                    label,
                    len(outliers),
                    outlier_threshold_percent,
                )
            )
        #endif
    #endif
#enddef


def print_final_summary(clean_summary_values):
    rga_values = []
    rgc_values = []
    overall_values = []

    for item in clean_summary_values:
        values = item["percent_differences"]

        overall_values.extend(values)

        if item["group"] == "RGA":
            rga_values.extend(values)
        elif item["group"] == "RGC":
            rgc_values.extend(values)
        else:
            fatal("Unknown summary group '{}'.".format(item["group"]))
        #endif
    #endfor

    rga_mean, rga_std = compute_mean_and_std(rga_values)
    rgc_mean, rgc_std = compute_mean_and_std(rgc_values)
    overall_mean, overall_std = compute_mean_and_std(overall_values)

    print("")
    print("============================================================")
    print("Final summary, ignoring outliers")
    print("Outlier definition used here: |percent difference| > 10%")
    print("Percent difference definition:")
    print("  100 * (HEL::Scaler - RUN::Scaler) / RUN::Scaler")
    print("")
    print("Mean RGA percent difference             = {:.6f}%".format(rga_mean))
    print("Std. dev. RGA percent difference        = {:.6f}%".format(rga_std))
    print("Mean RGC percent difference             = {:.6f}%".format(rgc_mean))
    print("Std. dev. RGC percent difference        = {:.6f}%".format(rgc_std))
    print("Mean overall percent difference         = {:.6f}%".format(overall_mean))
    print("Std. dev. overall percent difference    = {:.6f}%".format(overall_std))
    print("============================================================")
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


def plot_run_period_panel_with_outlier_axis(ax_top, ax_mid, ax_bottom, comparison_local, title):
    ax_top.plot(
        comparison_local["runnum"],
        comparison_local["run_scaler"],
        marker="o",
        linestyle="none",
        markersize=3,
        label="RUN::Scaler",
    )

    ax_top.plot(
        comparison_local["runnum"],
        comparison_local["hel_scaler"],
        marker="s",
        linestyle="none",
        markersize=3,
        label="HEL::Scaler",
    )

    ax_top.set_ylabel("Charge (nC)")
    ax_top.set_title(title)
    ax_top.legend(fontsize=7)
    ax_top.grid(True, alpha=0.3)

    ax_mid.axhline(0.0, linewidth=1)
    ax_mid.axhline(10.0, linewidth=1, linestyle="--")
    ax_mid.axhline(-10.0, linewidth=1, linestyle="--")
    ax_mid.plot(
        comparison_local["runnum"],
        comparison_local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=2,
    )

    ax_mid.set_ylabel("% diff.")
    ax_mid.grid(True, alpha=0.3)

    ax_bottom.axhline(0.0, linewidth=1)
    ax_bottom.plot(
        comparison_local["runnum"],
        comparison_local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=2,
    )

    ax_bottom.set_ylim(-10.0, 10.0)
    ax_bottom.set_xlabel("Run number")
    ax_bottom.set_ylabel("% diff.")
    ax_bottom.grid(True, alpha=0.3)
#enddef


def plot_run_period_panel_zoom_only(ax_top, ax_bottom, comparison_local, title):
    ax_top.plot(
        comparison_local["runnum"],
        comparison_local["run_scaler"],
        marker="o",
        linestyle="none",
        markersize=3,
        label="RUN::Scaler",
    )

    ax_top.plot(
        comparison_local["runnum"],
        comparison_local["hel_scaler"],
        marker="s",
        linestyle="none",
        markersize=3,
        label="HEL::Scaler",
    )

    ax_top.set_ylabel("Charge (nC)")
    ax_top.set_title(title)
    ax_top.legend(fontsize=7)
    ax_top.grid(True, alpha=0.3)

    ax_bottom.axhline(0.0, linewidth=1)
    ax_bottom.plot(
        comparison_local["runnum"],
        comparison_local["percent_difference"],
        marker="o",
        linestyle="none",
        markersize=2,
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
        fatal("No common run numbers found between Neupane and Hayward CSV files.")
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


def make_run_period_scaler_comparison_plot(run_period_sources, output_path):
    outlier_threshold_percent = 10.0
    clean_summary_values = []

    fig = plt.figure(
        figsize=(24, 14),
        constrained_layout=True,
    )

    outer = fig.add_gridspec(
        2,
        3,
        wspace=0.18,
        hspace=0.22,
    )

    for index, source in enumerate(run_period_sources):
        label = source["label"]
        group = source["group"]
        use_outlier_axis = source["use_outlier_axis"]

        local = read_charge_source(source)
        local = local.sort_values("runnum").reset_index(drop=True)

        comparison_local = compute_percent_difference(
            local,
            "hel_scaler",
            "run_scaler",
            "percent_difference",
        )

        if len(comparison_local) == 0:
            fatal("No nonzero RUN::Scaler and HEL::Scaler entries found for {}.".format(label))
        #endif

        print_run_period_summary(
            label,
            comparison_local,
            outlier_threshold_percent,
            print_outliers=use_outlier_axis,
        )

        non_outliers, outliers = split_outliers(
            comparison_local,
            outlier_threshold_percent,
        )

        clean_summary_values.append(
            {
                "label": label,
                "group": group,
                "percent_differences": non_outliers["percent_difference"].to_numpy().tolist(),
            }
        )

        row = index // 3
        col = index % 3

        if use_outlier_axis:
            inner = outer[row, col].subgridspec(
                3,
                1,
                height_ratios=[3, 1, 1],
                hspace=0.05,
            )

            ax_top = fig.add_subplot(inner[0])
            ax_mid = fig.add_subplot(inner[1], sharex=ax_top)
            ax_bottom = fig.add_subplot(inner[2], sharex=ax_top)

            plot_run_period_panel_with_outlier_axis(
                ax_top,
                ax_mid,
                ax_bottom,
                comparison_local,
                label,
            )

            plt.setp(ax_top.get_xticklabels(), visible=False)
            plt.setp(ax_mid.get_xticklabels(), visible=False)
        else:
            inner = outer[row, col].subgridspec(
                2,
                1,
                height_ratios=[3, 1],
                hspace=0.05,
            )

            ax_top = fig.add_subplot(inner[0])
            ax_bottom = fig.add_subplot(inner[1], sharex=ax_top)

            plot_run_period_panel_zoom_only(
                ax_top,
                ax_bottom,
                comparison_local,
                label,
            )

            plt.setp(ax_top.get_xticklabels(), visible=False)
        #endif
    #endfor

    fig.suptitle("RUN::Scaler vs HEL::Scaler accumulated charge comparison", fontsize=18)

    fig.savefig(output_path, dpi=300)
    plt.close(fig)

    print("")
    print("Saved:", output_path)

    print_final_summary(clean_summary_values)
#enddef


def main():
    rga_charge_dir = "/u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/imports/integrated_luminosity"
    rgc_run_info_path = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

    neupane_path = os.path.join(rga_charge_dir, "krishnas_charges.csv")
    hayward_path = os.path.join(rga_charge_dir, "global.csv")

    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    make_neupane_hayward_comparison_plot(
        neupane_path,
        hayward_path,
        os.path.join(output_dir, "charge_cross_check.png"),
    )

    run_period_sources = [
        {
            "label": "RGA Fa18 Inb",
            "group": "RGA",
            "kind": "flat",
            "path": os.path.join(rga_charge_dir, "rga_fa18_inb.txt"),
            "use_outlier_axis": True,
        },
        {
            "label": "RGA Fa18 Out",
            "group": "RGA",
            "kind": "flat",
            "path": os.path.join(rga_charge_dir, "rga_fa18_out.txt"),
            "use_outlier_axis": True,
        },
        {
            "label": "RGA Sp19 Inb",
            "group": "RGA",
            "kind": "flat",
            "path": os.path.join(rga_charge_dir, "rga_sp19_inb.txt"),
            "use_outlier_axis": True,
        },
        {
            "label": "RGC Su22 Inb",
            "group": "RGC",
            "kind": "combined_sectioned",
            "path": rgc_run_info_path,
            "section_headers": [
                "# RGC Su22 NH3",
                "# RGC Su22 Inb ND3 Runs",
            ],
            "use_outlier_axis": False,
        },
        {
            "label": "RGC Fa22 Inb",
            "group": "RGC",
            "kind": "combined_sectioned",
            "path": rgc_run_info_path,
            "section_headers": [
                "# RGC Fa22 NH3",
                "# RGC Fa22 Inb ND3 Runs",
            ],
            "use_outlier_axis": False,
        },
        {
            "label": "RGC Sp23",
            "group": "RGC",
            "kind": "combined_sectioned",
            "path": rgc_run_info_path,
            "section_headers": [
                "# RGC Sp23 NH3",
                "# RGC Sp23 ND3",
            ],
            "use_outlier_axis": False,
        },
    ]

    make_run_period_scaler_comparison_plot(
        run_period_sources,
        os.path.join(output_dir, "run_period_scaler_cross_check.png"),
    )
#enddef


if __name__ == "__main__":
    main()
#endif