#!/usr/bin/env python3

import os
import csv
import argparse
import numpy as np
import matplotlib.pyplot as plt


DEFAULT_INPUT_CSV = "non_qa_gated_totals.csv"
DEFAULT_QA_INPUT_CSV = "qa_gated_totals.csv"

DEFAULT_OUTPUT_PNG = "output/rga_rgb_rgc_non_gated_totals.png"
DEFAULT_OUTPUT_SORTED_CSV = "output/non_qa_gated_totals_sorted_by_period.csv"

RUN_OUTLIER_SIGMA_THRESHOLD = 5.0

# If the gap between adjacent run periods is larger than this,
# compress the empty space down to this size on the plot.
MAX_EMPTY_GAP_TO_KEEP = 350.0

# Minimum horizontal separation between x tick labels in compressed coordinates.
# Increase this if x labels still overlap.
MIN_XTICK_SEPARATION = 180.0


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compare RUN::Scaler and HEL::Scaler charge totals for "
            "non-QA-gated and QA-gated CSV files."
        )
    )

    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT_CSV,
        help=f"Input non-QA-gated sectioned CSV file. Default: {DEFAULT_INPUT_CSV}",
    )

    parser.add_argument(
        "--qa-input",
        default=DEFAULT_QA_INPUT_CSV,
        help=f"Input QA-gated sectioned CSV file. Default: {DEFAULT_QA_INPUT_CSV}",
    )

    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT_PNG,
        help=f"Output non-QA-gated plot path. Default: {DEFAULT_OUTPUT_PNG}",
    )

    parser.add_argument(
        "--sorted-output",
        default=DEFAULT_OUTPUT_SORTED_CSV,
        help=f"Output sorted non-QA-gated CSV path. Default: {DEFAULT_OUTPUT_SORTED_CSV}",
    )

    parser.add_argument(
        "--numerator",
        choices=["HEL", "RUN"],
        default="HEL",
        help=(
            "Top of the percent-difference numerator. "
            "'HEL' gives HEL - RUN. "
            "'RUN' gives RUN - HEL. "
            "Default: HEL"
        ),
    )

    parser.add_argument(
        "--denominator",
        choices=["RUN", "HEL"],
        default="RUN",
        help=(
            "Denominator of the percent difference. "
            "'RUN' divides by RUN::Scaler. "
            "'HEL' divides by HEL sum. "
            "Default: RUN"
        ),
    )

    parser.add_argument(
        "--outlier-sigma",
        type=float,
        default=RUN_OUTLIER_SIGMA_THRESHOLD,
        help=(
            "Robust median/MAD outlier threshold in sigma-like units. "
            f"Default: {RUN_OUTLIER_SIGMA_THRESHOLD}"
        ),
    )

    parser.add_argument(
        "--max-gap",
        type=float,
        default=MAX_EMPTY_GAP_TO_KEEP,
        help=(
            "Maximum empty run-number gap to keep in compressed x-axis coordinates. "
            f"Default: {MAX_EMPTY_GAP_TO_KEEP}"
        ),
    )

    parser.add_argument(
        "--min-xtick-separation",
        type=float,
        default=MIN_XTICK_SEPARATION,
        help=(
            "Minimum separation between compressed x tick labels. "
            f"Default: {MIN_XTICK_SEPARATION}"
        ),
    )

    return parser.parse_args()
#enddef


def derive_output_path(path, suffix):
    base, ext = os.path.splitext(path)

    if ext == "":
        ext = ".png"
    #endif

    return f"{base}{suffix}{ext}"
#enddef


def compute_percent_difference(run_scaler, hel_sum, numerator_choice, denominator_choice):
    """
    Compute percent difference with independent numerator and denominator choices.

    numerator_choice = "HEL":

        numerator = HEL - RUN

    numerator_choice = "RUN":

        numerator = RUN - HEL

    denominator_choice = "RUN":

        denominator = RUN

    denominator_choice = "HEL":

        denominator = HEL
    """

    if numerator_choice == "HEL":
        numerator = hel_sum - run_scaler
    elif numerator_choice == "RUN":
        numerator = run_scaler - hel_sum
    else:
        raise RuntimeError(f"Unknown numerator choice: {numerator_choice}")
    #endif

    if denominator_choice == "RUN":
        denominator = run_scaler
    elif denominator_choice == "HEL":
        denominator = hel_sum
    else:
        raise RuntimeError(f"Unknown denominator choice: {denominator_choice}")
    #endif

    if denominator == 0.0:
        return np.nan
    #endif

    return 100.0 * numerator / denominator
#enddef


def get_percent_difference_definition_string(numerator_choice, denominator_choice):
    if numerator_choice == "HEL":
        numerator_string = "HEL sum - RUN::Scaler"
    elif numerator_choice == "RUN":
        numerator_string = "RUN::Scaler - HEL sum"
    else:
        numerator_string = "unknown numerator"
    #endif

    if denominator_choice == "RUN":
        denominator_string = "RUN::Scaler"
    elif denominator_choice == "HEL":
        denominator_string = "HEL sum"
    else:
        denominator_string = "unknown denominator"
    #endif

    return f"100 * ({numerator_string}) / {denominator_string}"
#enddef


def get_ylabel(numerator_choice, denominator_choice):
    if numerator_choice == "HEL":
        numerator_latex = (
            r"\left("
            r"\mathrm{HEL}^{+} + \mathrm{HEL}^{-} + \mathrm{HEL}^{\mathrm{unassigned}}"
            r"\right)"
            r" - \mathrm{RUN}"
        )
    else:
        numerator_latex = (
            r"\mathrm{RUN} - "
            r"\left("
            r"\mathrm{HEL}^{+} + \mathrm{HEL}^{-} + \mathrm{HEL}^{\mathrm{unassigned}}"
            r"\right)"
        )
    #endif

    if denominator_choice == "RUN":
        denominator_latex = r"\mathrm{RUN}"
    else:
        denominator_latex = (
            r"\left("
            r"\mathrm{HEL}^{+} + \mathrm{HEL}^{-} + \mathrm{HEL}^{\mathrm{unassigned}}"
            r"\right)"
        )
    #endif

    return (
        r"$100 \times "
        r"\left["
        + numerator_latex
        + r"\right] / "
        + denominator_latex
        + r"$ [%]"
    )
#enddef


def read_period_charge_file(path, numerator_choice, denominator_choice):
    """
    Read a sectioned charge CSV.

    Expected format:

        # RGA Fa18 Inb
        run,RUN::Scaler,HEL::Scaler(+),HEL::Scaler(-),HEL::Scaler(unassigned),0,0

    The file may also contain QA-removed rows like:

        run,0.0,0.0,0.0,0.0

    where the percent difference becomes NaN because the denominator is zero.
    """

    if not os.path.exists(path):
        raise RuntimeError(f"Input CSV does not exist: {path}")
    #endif

    periods = []
    current_period = None

    with open(path, "r", newline="") as f:
        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue
            #endif

            if line.startswith("#"):
                period_name = line.lstrip("#").strip()

                current_period = {
                    "name": period_name,
                    "rows": [],
                }

                periods.append(current_period)
                continue
            #endif

            if current_period is None:
                raise RuntimeError(
                    f"Found data line before any period header in {path}: {line}"
                )
            #endif

            parts = [x.strip() for x in line.split(",")]

            if len(parts) < 5:
                raise RuntimeError(
                    f"Expected at least 5 comma-separated columns in {path}, "
                    f"got {len(parts)}: {line}"
                )
            #endif

            runnum = int(parts[0])
            run_scaler = float(parts[1])
            hel_pos = float(parts[2])
            hel_neg = float(parts[3])
            hel_unassigned = float(parts[4])
            hel_sum = hel_pos + hel_neg + hel_unassigned

            percent_difference = compute_percent_difference(
                run_scaler,
                hel_sum,
                numerator_choice,
                denominator_choice,
            )

            qadb_fully_removed = run_scaler == 0.0

            current_period["rows"].append(
                {
                    "period": current_period["name"],
                    "runnum": runnum,
                    "run_scaler": run_scaler,
                    "hel_pos": hel_pos,
                    "hel_neg": hel_neg,
                    "hel_unassigned": hel_unassigned,
                    "hel_sum": hel_sum,
                    "percent_difference": percent_difference,
                    "qadb_fully_removed": qadb_fully_removed,
                    "original_parts": parts,
                }
            )
        #endfor
    #endwith

    return periods
#enddef


def flatten_rows(periods):
    rows = []

    for period in periods:
        for row in period["rows"]:
            rows.append(row)
        #endfor
    #endfor

    return rows
#enddef


def make_run_lookup(periods):
    lookup = {}

    for period in periods:
        for row in period["rows"]:
            lookup[row["runnum"]] = row
        #endfor
    #endfor

    return lookup
#enddef


def find_qadb_removed_runs(nonqa_periods, qa_periods):
    """
    Return run numbers whose QA-gated RUN::Scaler column is exactly zero.

    The non-QA-gated rows are used for plotting the red open circles, because
    fully removed QA rows do not have a meaningful QA-gated percent difference.
    """

    nonqa_lookup = make_run_lookup(nonqa_periods)
    removed = []

    for qa_row in flatten_rows(qa_periods):
        if qa_row["run_scaler"] != 0.0:
            continue
        #endif

        runnum = qa_row["runnum"]

        if runnum not in nonqa_lookup:
            removed.append(
                {
                    "runnum": runnum,
                    "period": qa_row["period"],
                    "nonqa_row": None,
                    "qa_row": qa_row,
                }
            )
            continue
        #endif

        removed.append(
            {
                "runnum": runnum,
                "period": qa_row["period"],
                "nonqa_row": nonqa_lookup[runnum],
                "qa_row": qa_row,
            }
        )
    #endfor

    return removed
#enddef


def write_sorted_period_charge_file(periods, output_path):
    """
    Write a copy of the input CSV with the same period headers and same row
    contents, but with numeric rows sorted by run number inside each period.
    """

    output_dir = os.path.dirname(output_path)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    #endif

    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")

        for period in periods:
            f.write(f"# {period['name']}\n")

            sorted_rows = sorted(
                period["rows"],
                key=lambda row: row["runnum"],
            )

            for row in sorted_rows:
                writer.writerow(row["original_parts"])
            #endfor
        #endfor
    #endwith
#enddef


def constant_fit(values):
    """
    Constant fit for unweighted points.
    """

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    n = len(arr)

    if n == 0:
        return np.nan, np.nan, np.nan, 0
    #endif

    mean = float(np.mean(arr))

    if n > 1:
        rms = float(np.sqrt(np.mean((arr - mean) ** 2)))
        err = float(np.std(arr, ddof=1) / np.sqrt(n))
    else:
        rms = 0.0
        err = np.nan
    #endif

    return mean, err, rms, n
#enddef


def robust_sigma_from_mad(values):
    """
    Estimate a robust sigma using the median absolute deviation.

        sigma ~= 1.4826 * median(|x - median(x)|)

    For Gaussian-distributed data, this estimates the ordinary standard deviation.
    """

    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]

    if len(arr) == 0:
        return np.nan, np.nan
    #endif

    median_value = float(np.median(arr))
    mad_value = float(np.median(np.abs(arr - median_value)))

    robust_sigma = 1.4826 * mad_value

    return median_value, robust_sigma
#enddef


def find_run_outliers_iterative_robust(
    rows,
    runnums,
    percent_differences,
    sigma_threshold,
):
    """
    Identify individual run-level outliers inside one period using iterative
    robust median/MAD clipping.

    Procedure:
      1. Start with all finite points in the period.
      2. Compute the median percent difference.
      3. Compute robust_sigma = 1.4826 * MAD.
      4. Remove points farther than sigma_threshold * robust_sigma.
      5. Repeat until no additional points are removed.
    """

    runnums = np.asarray(runnums, dtype=float)
    values = np.asarray(percent_differences, dtype=float)

    finite_mask = np.isfinite(values)
    keep_mask = finite_mask.copy()

    outlier_info = []
    iteration = 0

    while True:
        iteration += 1

        current_values = values[keep_mask]

        if len(current_values) < 3:
            break
        #endif

        median_value, robust_sigma = robust_sigma_from_mad(current_values)

        if not np.isfinite(robust_sigma) or robust_sigma == 0.0:
            break
        #endif

        current_indices = np.where(keep_mask)[0]
        z_values_current = (values[current_indices] - median_value) / robust_sigma
        bad_local = np.abs(z_values_current) > sigma_threshold

        if not np.any(bad_local):
            break
        #endif

        bad_indices = current_indices[bad_local]

        for idx in bad_indices:
            keep_mask[idx] = False

            outlier_info.append(
                {
                    "iteration": iteration,
                    "period": rows[idx]["period"],
                    "runnum": int(runnums[idx]),
                    "run_scaler": float(rows[idx]["run_scaler"]),
                    "hel_sum": float(rows[idx]["hel_sum"]),
                    "percent_difference": float(values[idx]),
                    "median_used": median_value,
                    "robust_sigma_used": robust_sigma,
                    "z_value": float((values[idx] - median_value) / robust_sigma),
                }
            )
        #endfor
    #endwhile

    return keep_mask, outlier_info
#enddef


def process_periods(
    periods,
    dataset_label,
    numerator_choice,
    denominator_choice,
    outlier_sigma_threshold,
):
    """
    Run robust outlier identification and constant fits for one dataset.
    """

    period_fit_info = []
    all_outliers = []

    print()
    print("======================================================================")
    print(f"{dataset_label}")
    print("======================================================================")
    print("Run-level outlier search:")
    print(
        f"  Outlier definition inside each period: "
        f"|robust median/MAD z| > {outlier_sigma_threshold:.1f}"
    )
    print(
        "  z_i = (run_i percent difference - median of current kept runs) "
        "/ (1.4826 * MAD of current kept runs)"
    )
    print()
    print("Percent difference definition:")
    print(f"  {get_percent_difference_definition_string(numerator_choice, denominator_choice)}")
    print("  HEL sum = HEL+ + HEL- + HEL(unassigned)")
    print(f"  numerator choice   = {numerator_choice}")
    print(f"  denominator choice = {denominator_choice}")

    any_run_outliers = False

    for period in periods:
        period_name = period["name"]
        rows = sorted(period["rows"], key=lambda row: row["runnum"])

        if len(rows) == 0:
            continue
        #endif

        runnums = np.array([row["runnum"] for row in rows], dtype=float)
        percent_differences = np.array(
            [row["percent_difference"] for row in rows],
            dtype=float,
        )

        keep_mask, outlier_info = find_run_outliers_iterative_robust(
            rows,
            runnums,
            percent_differences,
            outlier_sigma_threshold,
        )

        all_outliers.extend(outlier_info)

        if len(outlier_info) > 0:
            any_run_outliers = True
            print()
            print(f"Outlier runs removed from period average for {period_name}:")
            print(
                f"{'iter':>4} "
                f"{'runnum':>8} "
                f"{'RUN::Scaler [nC]':>18} "
                f"{'HEL sum [nC]':>18} "
                f"{'percent_diff [%]':>18} "
                f"{'median [%]':>14} "
                f"{'robust sigma [%]':>18} "
                f"{'z':>12}"
            )

            for item in outlier_info:
                print(
                    f"{item['iteration']:>4d} "
                    f"{item['runnum']:>8d} "
                    f"{item['run_scaler']:>18.6f} "
                    f"{item['hel_sum']:>18.6f} "
                    f"{item['percent_difference']:>18.7f} "
                    f"{item['median_used']:>14.7f} "
                    f"{item['robust_sigma_used']:>18.7f} "
                    f"{item['z_value']:>12.3f}"
                )
            #endfor
        #endif

        fit_value, fit_err, fit_rms, n_fit = constant_fit(
            percent_differences[keep_mask]
        )

        period_fit_info.append(
            {
                "name": period_name,
                "dataset_label": dataset_label,
                "fit_value": fit_value,
                "fit_err": fit_err,
                "fit_rms": fit_rms,
                "n_fit": n_fit,
                "n_total_finite": int(np.sum(np.isfinite(percent_differences))),
                "n_total_rows": len(rows),
                "n_run_outliers": len(outlier_info),
                "rows": rows,
                "runnums": runnums,
                "percent_differences": percent_differences,
                "keep_mask": keep_mask,
                "outlier_info": outlier_info,
            }
        )
    #endfor

    if not any_run_outliers:
        print()
        print(
            f"No individual run-level outliers found using "
            f"|robust z| > {outlier_sigma_threshold:.1f}."
        )
    #endif

    clean_period_fit_values = [
        item["fit_value"] for item in period_fit_info
        if np.isfinite(item["fit_value"])
    ]

    if len(clean_period_fit_values) > 0:
        overall_period_average = float(np.mean(clean_period_fit_values))
        overall_period_rms = float(
            np.sqrt(
                np.mean(
                    (
                        np.asarray(clean_period_fit_values, dtype=float)
                        - overall_period_average
                    ) ** 2
                )
            )
        )
    else:
        overall_period_average = np.nan
        overall_period_rms = np.nan
    #endif

    print()
    print("Constant fits to percent difference after removing run-level outliers:")
    print(f"  percent difference = {get_percent_difference_definition_string(numerator_choice, denominator_choice)}")
    print("  sum(HEL::Scaler) = HEL+ + HEL- + HEL(unassigned)")
    print("-" * 119)
    print(
        f"{'Period':<20} "
        f"{'N used':>8} "
        f"{'N finite':>8} "
        f"{'N rows':>8} "
        f"{'N out':>7} "
        f"{'constant [%]':>16} "
        f"{'err [%]':>14} "
        f"{'RMS [%]':>14}"
    )
    print("-" * 119)

    for item in period_fit_info:
        print(
            f"{item['name']:<20} "
            f"{item['n_fit']:>8d} "
            f"{item['n_total_finite']:>8d} "
            f"{item['n_total_rows']:>8d} "
            f"{item['n_run_outliers']:>7d} "
            f"{item['fit_value']:>16.7f} "
            f"{item['fit_err']:>14.7f} "
            f"{item['fit_rms']:>14.7f}"
        )
    #endfor

    print("-" * 119)
    print()
    print("Overall average of period constants after run-level outlier removal:")
    print(f"  average = {overall_period_average:.7f}%")
    print(f"  RMS     = {overall_period_rms:.7f}%")
    print()

    return {
        "dataset_label": dataset_label,
        "period_fit_info": period_fit_info,
        "all_outliers": all_outliers,
        "overall_period_average": overall_period_average,
        "overall_period_rms": overall_period_rms,
    }
#enddef


def build_compressed_x_mapping(period_fit_info, max_empty_gap_to_keep):
    """
    Build a piecewise-linear compressed x-axis.

    The run numbers inside each period keep their true internal spacing.
    Large gaps between periods are compressed down to max_empty_gap_to_keep.
    """

    valid_periods = []

    for item in period_fit_info:
        percent_differences = item["percent_differences"]
        finite_mask = np.isfinite(percent_differences)

        if not np.any(finite_mask):
            continue
        #endif

        runnums = item["runnums"][finite_mask]

        valid_periods.append(
            {
                "name": item["name"],
                "real_min": float(np.min(runnums)),
                "real_max": float(np.max(runnums)),
            }
        )
    #endfor

    if len(valid_periods) == 0:
        return {
            "valid_periods": [],
            "gap_breaks": [],
        }
    #endif

    current_x_start = valid_periods[0]["real_min"]
    valid_periods[0]["plot_min"] = current_x_start
    valid_periods[0]["plot_max"] = current_x_start + (
        valid_periods[0]["real_max"] - valid_periods[0]["real_min"]
    )

    gap_breaks = []

    for i_period in range(1, len(valid_periods)):
        previous = valid_periods[i_period - 1]
        current = valid_periods[i_period]

        real_gap = current["real_min"] - previous["real_max"]
        kept_gap = min(real_gap, max_empty_gap_to_keep)

        current["plot_min"] = previous["plot_max"] + kept_gap
        current["plot_max"] = current["plot_min"] + (
            current["real_max"] - current["real_min"]
        )

        if real_gap > max_empty_gap_to_keep:
            skipped_low = previous["real_max"]
            skipped_high = current["real_min"]
            boundary_x = 0.5 * (previous["plot_max"] + current["plot_min"])

            gap_breaks.append(
                {
                    "left_period": previous["name"],
                    "right_period": current["name"],
                    "real_gap": real_gap,
                    "kept_gap": kept_gap,
                    "skipped_low": skipped_low,
                    "skipped_high": skipped_high,
                    "plot_x": boundary_x,
                }
            )
        #endif
    #endfor

    return {
        "valid_periods": valid_periods,
        "gap_breaks": gap_breaks,
    }
#enddef


def compress_runnums_for_period(runnums, period_name, mapping_info):
    """
    Convert real run numbers to compressed plot coordinates for one period.
    """

    runnums = np.asarray(runnums, dtype=float)

    for period_info in mapping_info["valid_periods"]:
        if period_info["name"] != period_name:
            continue
        #endif

        return period_info["plot_min"] + (runnums - period_info["real_min"])
    #endfor

    return runnums.copy()
#enddef


def make_compressed_xticks(mapping_info, min_tick_separation):
    """
    Construct x-axis tick positions and labels for the compressed axis.
    """

    candidate_ticks = []

    for period_info in mapping_info["valid_periods"]:
        real_min = period_info["real_min"]
        real_max = period_info["real_max"]
        plot_min = period_info["plot_min"]

        real_span = real_max - real_min

        if real_span <= 0.0:
            real_ticks = [real_min]
        elif real_span < 400.0:
            real_ticks = [
                real_min,
                real_max,
            ]
        else:
            real_ticks = [
                real_min,
                0.5 * (real_min + real_max),
                real_max,
            ]
        #endif

        for real_tick in real_ticks:
            plot_tick = plot_min + (real_tick - real_min)

            candidate_ticks.append(
                {
                    "plot_tick": plot_tick,
                    "real_tick": real_tick,
                }
            )
        #endfor
    #endfor

    candidate_ticks = sorted(candidate_ticks, key=lambda item: item["plot_tick"])

    kept_ticks = []

    for item in candidate_ticks:
        if len(kept_ticks) == 0:
            kept_ticks.append(item)
            continue
        #endif

        previous = kept_ticks[-1]

        if item["plot_tick"] - previous["plot_tick"] >= min_tick_separation:
            kept_ticks.append(item)
        else:
            continue
        #endif
    #endfor

    tick_positions = [item["plot_tick"] for item in kept_ticks]
    tick_labels = [str(int(round(item["real_tick"]))) for item in kept_ticks]

    return tick_positions, tick_labels
#enddef


def plot_processed_dataset(
    processed,
    output_png,
    title,
    numerator_choice,
    denominator_choice,
    max_empty_gap_to_keep,
    min_xtick_separation,
    qadb_removed_runs=None,
):
    """
    Plot one processed dataset.

    If qadb_removed_runs is provided, red open circles are drawn around the
    corresponding non-QA-gated points.
    """

    period_fit_info = processed["period_fit_info"]
    overall_period_average = processed["overall_period_average"]

    output_dir = os.path.dirname(output_png)

    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    #endif

    mapping_info = build_compressed_x_mapping(
        period_fit_info,
        max_empty_gap_to_keep,
    )

    if len(mapping_info["gap_breaks"]) > 0:
        print(f"Compressed large run-number gaps on plot {output_png}:")
        for gap in mapping_info["gap_breaks"]:
            print(
                f"  {gap['left_period']} -> {gap['right_period']}: "
                f"real gap {gap['skipped_low']:.0f} to {gap['skipped_high']:.0f} "
                f"compressed to {gap['kept_gap']:.0f} run-number units"
            )
        #endfor
        print()
    #endif

    removed_run_set = set()

    if qadb_removed_runs is not None:
        removed_run_set = {
            item["runnum"] for item in qadb_removed_runs
            if item["nonqa_row"] is not None
            and np.isfinite(item["nonqa_row"]["percent_difference"])
        }
    #endif

    fig, ax = plt.subplots(figsize=(18, 7))

    all_plot_runs = []
    removed_plot_x = []
    removed_plot_y = []

    for i_period, item in enumerate(period_fit_info):
        period_name = item["name"]
        runnums = item["runnums"]
        percent_differences = item["percent_differences"]
        keep_mask = item["keep_mask"]
        finite_mask = np.isfinite(percent_differences)
        outlier_mask = finite_mask & (~keep_mask)

        if not np.any(finite_mask):
            continue
        #endif

        plot_runnums = compress_runnums_for_period(
            runnums,
            period_name,
            mapping_info,
        )

        fit_value = item["fit_value"]

        normal_points = ax.plot(
            plot_runnums[keep_mask],
            percent_differences[keep_mask],
            marker="o",
            linestyle="none",
            markersize=4,
            label=period_name,
        )

        period_color = normal_points[0].get_color()

        if np.any(outlier_mask):
            ax.plot(
                plot_runnums[outlier_mask],
                percent_differences[outlier_mask],
                marker="x",
                linestyle="none",
                markersize=7,
                markeredgewidth=1.4,
                color=period_color,
            )
        #endif

        for i_point in range(len(runnums)):
            runnum = int(runnums[i_point])

            if runnum not in removed_run_set:
                continue
            #endif

            if not np.isfinite(percent_differences[i_point]):
                continue
            #endif

            removed_plot_x.append(plot_runnums[i_point])
            removed_plot_y.append(percent_differences[i_point])
        #endfor

        xmin = float(np.min(plot_runnums[finite_mask]))
        xmax = float(np.max(plot_runnums[finite_mask]))

        ax.hlines(
            fit_value,
            xmin,
            xmax,
            linewidth=1.2,
            linestyle="-",
            color=period_color,
        )

        text_x = 0.5 * (xmin + xmax)

        if i_period % 2 == 0:
            text_y = -18.8
            text_va = "bottom"
        else:
            text_y = 18.8
            text_va = "top"
        #endif

        label_text = f"{period_name}: C = {fit_value:.3f}%"

        if item["n_run_outliers"] > 0:
            label_text += f" ({item['n_run_outliers']} out)"
        #endif

        ax.text(
            text_x,
            text_y,
            label_text,
            ha="center",
            va=text_va,
            fontsize=8,
            rotation=0,
            clip_on=False,
            bbox={
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.75,
                "pad": 1.5,
            },
        )

        all_plot_runs.extend(plot_runnums[finite_mask].tolist())
    #endfor

    if len(removed_plot_x) > 0:
        ax.plot(
            removed_plot_x,
            removed_plot_y,
            marker="o",
            linestyle="none",
            markersize=11,
            markerfacecolor="none",
            markeredgecolor="red",
            markeredgewidth=1.5,
            label="QADB fully removed",
        )
    #endif

    for gap in mapping_info["gap_breaks"]:
        gap_x = gap["plot_x"]

        ax.axvline(
            gap_x,
            color="gray",
            linewidth=0.9,
            alpha=0.8,
            linestyle="--",
        )

        ax.text(
            gap_x,
            0.96,
            "//",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=16,
            fontweight="bold",
            color="gray",
        )

        ax.text(
            gap_x,
            0.04,
            "//",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=16,
            fontweight="bold",
            color="gray",
        )
    #endfor

    ax.axhline(
        0.0,
        color="black",
        linewidth=0.8,
        linestyle="--",
    )

    if np.isfinite(overall_period_average):
        ax.axhline(
            overall_period_average,
            color="black",
            linewidth=0.9,
            linestyle=":",
        )
    #endif

    tick_positions, tick_labels = make_compressed_xticks(
        mapping_info,
        min_xtick_separation,
    )

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(
        tick_labels,
        rotation=35,
        ha="right",
        fontsize=8,
    )

    ax.set_xlabel("Run Number")
    ax.set_ylabel(get_ylabel(numerator_choice, denominator_choice))
    ax.set_ylim(-20.0, 20.0)

    if len(all_plot_runs) > 0:
        xmin_all = min(all_plot_runs)
        xmax_all = max(all_plot_runs)
        dx = xmax_all - xmin_all

        if dx == 0.0:
            ax.set_xlim(xmin_all - 1.0, xmax_all + 1.0)
        else:
            ax.set_xlim(xmin_all - 0.02 * dx, xmax_all + 0.02 * dx)
        #endif
    #endif

    ax.set_title(title)
    ax.grid(True, alpha=0.25)

    ax.legend(
        fontsize=8,
        ncol=2,
        loc="upper right",
        frameon=True,
    )

    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.close(fig)

    print(f"Saved plot to: {output_png}")
#enddef


def outlier_lookup(processed):
    lookup = {}

    for item in processed["all_outliers"]:
        lookup[item["runnum"]] = item
    #endfor

    return lookup
#enddef


def row_lookup_from_period_fit_info(period_fit_info):
    lookup = {}

    for period in period_fit_info:
        for row in period["rows"]:
            lookup[row["runnum"]] = row
        #endfor
    #endfor

    return lookup
#enddef


def print_qadb_removed_runs(qadb_removed_runs):
    print()
    print("======================================================================")
    print("Runs fully removed by QADB")
    print("======================================================================")
    print("Definition: QA-gated RUN::Scaler column is exactly zero.")
    print()

    if len(qadb_removed_runs) == 0:
        print("No fully QADB-removed runs found.")
        return
    #endif

    print(
        f"{'period':<20} "
        f"{'runnum':>8} "
        f"{'nonQA RUN [nC]':>18} "
        f"{'nonQA HEL sum [nC]':>20} "
        f"{'nonQA percent diff [%]':>24}"
    )

    for item in qadb_removed_runs:
        nonqa_row = item["nonqa_row"]
        qa_row = item["qa_row"]

        if nonqa_row is None:
            print(
                f"{qa_row['period']:<20} "
                f"{qa_row['runnum']:>8d} "
                f"{'missing':>18} "
                f"{'missing':>20} "
                f"{'missing':>24}"
            )
            continue
        #endif

        print(
            f"{nonqa_row['period']:<20} "
            f"{nonqa_row['runnum']:>8d} "
            f"{nonqa_row['run_scaler']:>18.6f} "
            f"{nonqa_row['hel_sum']:>20.6f} "
            f"{nonqa_row['percent_difference']:>24.7f}"
        )
    #endfor
#enddef


def print_method_specific_outliers(nonqa_processed, qa_processed, qadb_removed_runs):
    """
    Print runs that are outliers in non-QA but not QA, or QA but not non-QA.

    Fully QADB-removed runs are reported separately, because the QA-gated
    percent difference is undefined when the QA-gated charge is zero.
    """

    nonqa_outliers = outlier_lookup(nonqa_processed)
    qa_outliers = outlier_lookup(qa_processed)

    nonqa_rows = row_lookup_from_period_fit_info(nonqa_processed["period_fit_info"])
    qa_rows = row_lookup_from_period_fit_info(qa_processed["period_fit_info"])

    qadb_removed_set = {item["runnum"] for item in qadb_removed_runs}

    nonqa_only = sorted(
        runnum for runnum in nonqa_outliers.keys()
        if runnum not in qa_outliers
        and runnum not in qadb_removed_set
    )

    qa_only = sorted(
        runnum for runnum in qa_outliers.keys()
        if runnum not in nonqa_outliers
    )

    nonqa_outlier_qadb_removed = sorted(
        runnum for runnum in nonqa_outliers.keys()
        if runnum in qadb_removed_set
    )

    print()
    print("======================================================================")
    print("Method-specific scaler-bank discrepancies")
    print("======================================================================")
    print("These are runs that are robust median/MAD outliers in one dataset but not the other.")
    print("Fully QADB-removed runs are listed separately because their QA-gated percent")
    print("difference is undefined.")
    print()

    if len(nonqa_only) == 0:
        print("No runs are non-QA-gated outliers while remaining non-outliers in QA-gated data.")
    else:
        print("Runs that are outliers in NON-QA-gated data but not in QA-gated data:")
        print(
            f"{'runnum':>8} "
            f"{'period':<20} "
            f"{'nonQA diff [%]':>18} "
            f"{'nonQA z':>12} "
            f"{'QA diff [%]':>18} "
            f"{'QA RUN [nC]':>16} "
            f"{'QA HEL sum [nC]':>18}"
        )

        for runnum in nonqa_only:
            nonqa_out = nonqa_outliers[runnum]
            qa_row = qa_rows.get(runnum)

            if qa_row is None or not np.isfinite(qa_row["percent_difference"]):
                qa_diff_string = "undefined"
                qa_run_string = "missing"
                qa_hel_string = "missing"
            else:
                qa_diff_string = f"{qa_row['percent_difference']:.7f}"
                qa_run_string = f"{qa_row['run_scaler']:.6f}"
                qa_hel_string = f"{qa_row['hel_sum']:.6f}"
            #endif

            print(
                f"{runnum:>8d} "
                f"{nonqa_out['period']:<20} "
                f"{nonqa_out['percent_difference']:>18.7f} "
                f"{nonqa_out['z_value']:>12.3f} "
                f"{qa_diff_string:>18} "
                f"{qa_run_string:>16} "
                f"{qa_hel_string:>18}"
            )
        #endfor
    #endif

    print()

    if len(qa_only) == 0:
        print("No runs are QA-gated outliers while remaining non-outliers in non-QA-gated data.")
    else:
        print("Runs that are outliers in QA-gated data but not in NON-QA-gated data:")
        print(
            f"{'runnum':>8} "
            f"{'period':<20} "
            f"{'QA diff [%]':>18} "
            f"{'QA z':>12} "
            f"{'nonQA diff [%]':>18} "
            f"{'nonQA RUN [nC]':>18} "
            f"{'nonQA HEL sum [nC]':>20}"
        )

        for runnum in qa_only:
            qa_out = qa_outliers[runnum]
            nonqa_row = nonqa_rows.get(runnum)

            if nonqa_row is None or not np.isfinite(nonqa_row["percent_difference"]):
                nonqa_diff_string = "undefined"
                nonqa_run_string = "missing"
                nonqa_hel_string = "missing"
            else:
                nonqa_diff_string = f"{nonqa_row['percent_difference']:.7f}"
                nonqa_run_string = f"{nonqa_row['run_scaler']:.6f}"
                nonqa_hel_string = f"{nonqa_row['hel_sum']:.6f}"
            #endif

            print(
                f"{runnum:>8d} "
                f"{qa_out['period']:<20} "
                f"{qa_out['percent_difference']:>18.7f} "
                f"{qa_out['z_value']:>12.3f} "
                f"{nonqa_diff_string:>18} "
                f"{nonqa_run_string:>18} "
                f"{nonqa_hel_string:>20}"
            )
        #endfor
    #endif

    print()

    if len(nonqa_outlier_qadb_removed) == 0:
        print("No non-QA-gated outlier runs were fully removed by QADB.")
    else:
        print("Runs that were non-QA-gated outliers and then fully removed by QADB:")
        print(
            f"{'runnum':>8} "
            f"{'period':<20} "
            f"{'nonQA diff [%]':>18} "
            f"{'nonQA z':>12} "
            f"{'nonQA RUN [nC]':>18} "
            f"{'nonQA HEL sum [nC]':>20}"
        )

        for runnum in nonqa_outlier_qadb_removed:
            nonqa_out = nonqa_outliers[runnum]
            nonqa_row = nonqa_rows.get(runnum)

            print(
                f"{runnum:>8d} "
                f"{nonqa_out['period']:<20} "
                f"{nonqa_out['percent_difference']:>18.7f} "
                f"{nonqa_out['z_value']:>12.3f} "
                f"{nonqa_row['run_scaler']:>18.6f} "
                f"{nonqa_row['hel_sum']:>20.6f}"
            )
        #endfor
    #endif
#enddef


def main():
    args = parse_args()

    input_csv = args.input
    qa_input_csv = args.qa_input
    output_png = args.output
    output_sorted_csv = args.sorted_output
    numerator_choice = args.numerator
    denominator_choice = args.denominator
    outlier_sigma_threshold = args.outlier_sigma
    max_empty_gap_to_keep = args.max_gap
    min_xtick_separation = args.min_xtick_separation

    qadb_overlay_png = derive_output_path(output_png, "_qadb_removed_circled")
    qa_output_png = derive_output_path(output_png, "_qa_gated")

    qa_sorted_output_csv = derive_output_path(output_sorted_csv, "_qa_gated")

    nonqa_periods = read_period_charge_file(
        input_csv,
        numerator_choice,
        denominator_choice,
    )

    qa_periods = read_period_charge_file(
        qa_input_csv,
        numerator_choice,
        denominator_choice,
    )

    write_sorted_period_charge_file(
        nonqa_periods,
        output_sorted_csv,
    )

    write_sorted_period_charge_file(
        qa_periods,
        qa_sorted_output_csv,
    )

    qadb_removed_runs = find_qadb_removed_runs(
        nonqa_periods,
        qa_periods,
    )

    print_qadb_removed_runs(qadb_removed_runs)

    nonqa_processed = process_periods(
        nonqa_periods,
        "NON-QA-gated charge totals",
        numerator_choice,
        denominator_choice,
        outlier_sigma_threshold,
    )

    qa_processed = process_periods(
        qa_periods,
        "QA-gated charge totals",
        numerator_choice,
        denominator_choice,
        outlier_sigma_threshold,
    )

    print_method_specific_outliers(
        nonqa_processed,
        qa_processed,
        qadb_removed_runs,
    )

    plot_processed_dataset(
        nonqa_processed,
        output_png,
        "Non-QA-gated RUN::Scaler vs HEL::Scaler charge comparison by run period",
        numerator_choice,
        denominator_choice,
        max_empty_gap_to_keep,
        min_xtick_separation,
        qadb_removed_runs=None,
    )

    plot_processed_dataset(
        nonqa_processed,
        qadb_overlay_png,
        "Non-QA-gated charge comparison with QADB-removed runs circled",
        numerator_choice,
        denominator_choice,
        max_empty_gap_to_keep,
        min_xtick_separation,
        qadb_removed_runs=qadb_removed_runs,
    )

    plot_processed_dataset(
        qa_processed,
        qa_output_png,
        "QA-gated RUN::Scaler vs HEL::Scaler charge comparison by run period",
        numerator_choice,
        denominator_choice,
        max_empty_gap_to_keep,
        min_xtick_separation,
        qadb_removed_runs=None,
    )

    print(f"Saved sorted non-QA CSV to: {output_sorted_csv}")
    print(f"Saved sorted QA CSV to: {qa_sorted_output_csv}")
    print(f"Saved non-QA plot to: {output_png}")
    print(f"Saved non-QA plot with QADB-removed runs circled to: {qadb_overlay_png}")
    print(f"Saved QA-gated plot to: {qa_output_png}")
#enddef


if __name__ == "__main__":
    main()
#endif