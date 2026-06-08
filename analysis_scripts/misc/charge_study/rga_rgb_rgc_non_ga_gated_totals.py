#!/usr/bin/env python3

import os
import csv
import numpy as np
import matplotlib.pyplot as plt


INPUT_CSV = "non_qa_gated_totals.csv"
OUTPUT_PNG = "output/rga_rgb_rgc_non_gated_totals.png"
OUTPUT_SORTED_CSV = "output/non_qa_gated_totals_sorted_by_period.csv"

RUN_OUTLIER_SIGMA_THRESHOLD = 5.0

# If the gap between adjacent run periods is larger than this,
# compress the empty space down to this size on the plot.
MAX_EMPTY_GAP_TO_KEEP = 350.0


def read_period_charge_file(path):
    """
    Read a sectioned charge CSV.

    Expected format:

        # RGA Fa18 Inb
        run,RUN::Scaler,HEL::Scaler(+),HEL::Scaler(-),HEL::Scaler(unassigned),0,0
        ...

    Convention
    ----------
    Plotted quantity:

        percent_difference = 100 * (HEL::Scaler - RUN::Scaler) / RUN::Scaler

    where:

        HEL::Scaler = HEL+ + HEL- + HEL(unassigned)
    """

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
                    f"Found data line before any period header: {line}"
                )
            #endif

            parts = [x.strip() for x in line.split(",")]

            if len(parts) < 5:
                raise RuntimeError(
                    f"Expected at least 5 comma-separated columns, got {len(parts)}: {line}"
                )
            #endif

            runnum = int(parts[0])
            run_scaler = float(parts[1])
            hel_pos = float(parts[2])
            hel_neg = float(parts[3])
            hel_unassigned = float(parts[4])
            hel_sum = hel_pos + hel_neg + hel_unassigned

            if run_scaler == 0.0:
                percent_difference = np.nan
            else:
                percent_difference = 100.0 * (hel_sum - run_scaler) / run_scaler
            #endif

            current_period["rows"].append(
                {
                    "runnum": runnum,
                    "run_scaler": run_scaler,
                    "hel_pos": hel_pos,
                    "hel_neg": hel_neg,
                    "hel_unassigned": hel_unassigned,
                    "hel_sum": hel_sum,
                    "percent_difference": percent_difference,
                    "original_parts": parts,
                }
            )
        #endfor
    #endwith

    return periods


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


def constant_fit(values):
    """
    Constant fit for unweighted points.

    Returns
    -------
    mean : float
        Best-fit constant for equal point weights.
    err : float
        Standard error on the mean using sample std / sqrt(N).
    rms : float
        Population RMS scatter about the mean.
    n : int
        Number of finite points.
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


def find_run_outliers_leave_one_out(runnums, percent_differences, sigma_threshold):
    """
    Identify individual run-level outliers inside one period.

    For each run i, compare its percent difference y_i against the mean and RMS
    of all other runs in the same period:

        z_i = (y_i - mean_others) / RMS_others

    If |z_i| > sigma_threshold, the run is flagged as an outlier.
    """

    runnums = np.asarray(runnums, dtype=float)
    values = np.asarray(percent_differences, dtype=float)

    finite_mask = np.isfinite(values)
    keep_mask = finite_mask.copy()
    z_values = np.full(len(values), np.nan, dtype=float)
    outlier_info = []

    finite_indices = np.where(finite_mask)[0]

    if len(finite_indices) < 3:
        return keep_mask, outlier_info, z_values
    #endif

    for idx in finite_indices:
        other_indices = finite_indices[finite_indices != idx]
        other_values = values[other_indices]

        if len(other_values) < 2:
            continue
        #endif

        mean_others = float(np.mean(other_values))
        rms_others = float(np.sqrt(np.mean((other_values - mean_others) ** 2)))

        if rms_others == 0.0:
            continue
        #endif

        z_value = (values[idx] - mean_others) / rms_others
        z_values[idx] = z_value

        if abs(z_value) > sigma_threshold:
            keep_mask[idx] = False

            outlier_info.append(
                {
                    "runnum": int(runnums[idx]),
                    "percent_difference": float(values[idx]),
                    "mean_others": mean_others,
                    "rms_others": rms_others,
                    "z_value": float(z_value),
                }
            )
        #endif
    #endfor

    return keep_mask, outlier_info, z_values


def build_compressed_x_mapping(period_fit_info, max_empty_gap_to_keep):
    """
    Build a piecewise-linear compressed x-axis.

    The run numbers inside each period keep their true internal spacing.
    Large gaps between periods are compressed down to max_empty_gap_to_keep.

    Returns
    -------
    mapping_info : dict
        Contains per-period x offsets and boundary information.
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


def make_compressed_xticks(mapping_info):
    """
    Construct x-axis tick positions and labels for the compressed axis.
    """

    tick_positions = []
    tick_labels = []

    for period_info in mapping_info["valid_periods"]:
        real_min = period_info["real_min"]
        real_max = period_info["real_max"]
        plot_min = period_info["plot_min"]

        real_span = real_max - real_min

        if real_span <= 0.0:
            real_ticks = [real_min]
        elif real_span < 250.0:
            real_ticks = [
                real_min,
                real_max,
            ]
        elif real_span < 800.0:
            real_ticks = [
                real_min,
                0.5 * (real_min + real_max),
                real_max,
            ]
        else:
            real_ticks = [
                real_min,
                real_min + 0.33 * real_span,
                real_min + 0.67 * real_span,
                real_max,
            ]
        #endif

        for real_tick in real_ticks:
            plot_tick = plot_min + (real_tick - real_min)

            tick_positions.append(plot_tick)
            tick_labels.append(str(int(round(real_tick))))
        #endfor
    #endfor

    return tick_positions, tick_labels


def main():
    periods = read_period_charge_file(INPUT_CSV)

    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)

    write_sorted_period_charge_file(
        periods,
        OUTPUT_SORTED_CSV,
    )

    period_fit_info = []

    print()
    print("Run-level outlier search:")
    print(
        f"  Outlier definition inside each period: "
        f"|leave-one-out z| > {RUN_OUTLIER_SIGMA_THRESHOLD:.1f}"
    )
    print(
        "  z_i = (run_i percent difference - mean of other runs in same period) "
        "/ RMS of other runs in same period"
    )

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

        keep_mask, outlier_info, z_values = find_run_outliers_leave_one_out(
            runnums,
            percent_differences,
            RUN_OUTLIER_SIGMA_THRESHOLD,
        )

        if len(outlier_info) > 0:
            any_run_outliers = True
            print()
            print(f"Outlier runs removed from period average for {period_name}:")
            print(
                f"{'runnum':>8} "
                f"{'percent_diff [%]':>18} "
                f"{'mean_others [%]':>18} "
                f"{'RMS_others [%]':>18} "
                f"{'z':>12}"
            )

            for item in outlier_info:
                print(
                    f"{item['runnum']:>8d} "
                    f"{item['percent_difference']:>18.7f} "
                    f"{item['mean_others']:>18.7f} "
                    f"{item['rms_others']:>18.7f} "
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
                "fit_value": fit_value,
                "fit_err": fit_err,
                "fit_rms": fit_rms,
                "n_fit": n_fit,
                "n_total_finite": int(np.sum(np.isfinite(percent_differences))),
                "n_run_outliers": len(outlier_info),
                "rows": rows,
                "runnums": runnums,
                "percent_differences": percent_differences,
                "keep_mask": keep_mask,
                "outlier_info": outlier_info,
                "z_values": z_values,
            }
        )
    #endfor

    if not any_run_outliers:
        print()
        print(
            f"No individual run-level outliers found using "
            f"|z| > {RUN_OUTLIER_SIGMA_THRESHOLD:.1f}."
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
    print("  percent difference = 100 * (sum(HEL::Scaler) - RUN::Scaler) / RUN::Scaler")
    print("  sum(HEL::Scaler) = HEL+ + HEL- + HEL(unassigned)")
    print("-" * 104)
    print(
        f"{'Period':<20} "
        f"{'N used':>8} "
        f"{'N total':>8} "
        f"{'N out':>7} "
        f"{'constant [%]':>16} "
        f"{'err [%]':>14} "
        f"{'RMS [%]':>14}"
    )
    print("-" * 104)

    for item in period_fit_info:
        print(
            f"{item['name']:<20} "
            f"{item['n_fit']:>8d} "
            f"{item['n_total_finite']:>8d} "
            f"{item['n_run_outliers']:>7d} "
            f"{item['fit_value']:>16.7f} "
            f"{item['fit_err']:>14.7f} "
            f"{item['fit_rms']:>14.7f}"
        )
    #endfor

    print("-" * 104)
    print()
    print("Overall average of period constants after run-level outlier removal:")
    print(f"  average = {overall_period_average:.7f}%")
    print(f"  RMS     = {overall_period_rms:.7f}%")
    print()

    mapping_info = build_compressed_x_mapping(
        period_fit_info,
        MAX_EMPTY_GAP_TO_KEEP,
    )

    if len(mapping_info["gap_breaks"]) > 0:
        print("Compressed large run-number gaps on the plot:")
        for gap in mapping_info["gap_breaks"]:
            print(
                f"  {gap['left_period']} -> {gap['right_period']}: "
                f"real gap {gap['skipped_low']:.0f} to {gap['skipped_high']:.0f} "
                f"compressed to {gap['kept_gap']:.0f} run-number units"
            )
        #endfor
        print()
    #endif

    fig, ax = plt.subplots(figsize=(18, 7))

    all_plot_runs = []
    all_percent_differences = []

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

        ax.plot(
            plot_runnums[keep_mask],
            percent_differences[keep_mask],
            marker="o",
            linestyle="none",
            markersize=4,
            label=period_name,
        )

        if np.any(outlier_mask):
            ax.plot(
                plot_runnums[outlier_mask],
                percent_differences[outlier_mask],
                marker="x",
                linestyle="none",
                markersize=7,
                markeredgewidth=1.4,
            )
        #endif

        xmin = float(np.min(plot_runnums[finite_mask]))
        xmax = float(np.max(plot_runnums[finite_mask]))

        ax.hlines(
            fit_value,
            xmin,
            xmax,
            linewidth=1.2,
            linestyle="-",
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
        all_percent_differences.extend(percent_differences[finite_mask].tolist())
    #endfor

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

    tick_positions, tick_labels = make_compressed_xticks(mapping_info)

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, rotation=0)

    ax.set_xlabel("Run Number")
    ax.set_ylabel(
        r"$100 \times "
        r"\left["
        r"\left("
        r"\mathrm{HEL}^{+} + \mathrm{HEL}^{-} + \mathrm{HEL}^{\mathrm{unassigned}}"
        r"\right)"
        r" - \mathrm{RUN}"
        r"\right] / \mathrm{RUN}$ [%]"
    )

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

    ax.set_title(
        "Non-QA-gated RUN::Scaler vs HEL::Scaler charge comparison by run period"
    )

    ax.grid(True, alpha=0.25)

    ax.legend(
        fontsize=8,
        ncol=2,
        loc="upper right",
        frameon=True,
    )

    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300)
    plt.close(fig)

    print(f"Saved plot to: {OUTPUT_PNG}")
    print(f"Saved sorted CSV to: {OUTPUT_SORTED_CSV}")


if __name__ == "__main__":
    main()
#endif