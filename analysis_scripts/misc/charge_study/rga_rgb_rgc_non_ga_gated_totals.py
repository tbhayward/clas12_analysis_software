#!/usr/bin/env python3

import os
import csv
import numpy as np
import matplotlib.pyplot as plt


INPUT_CSV = "non_qa_gated_totals.csv"
OUTPUT_PNG = "output/rga_rgb_rgc_non_gated_totals.png"
OUTPUT_SORTED_CSV = "output/non_qa_gated_totals_sorted_by_period.csv"

OUTLIER_SIGMA_THRESHOLD = 5.0


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


def find_period_outliers(period_fit_info, sigma_threshold):
    """
    Identify period-level outliers using a leave-one-out comparison.

    For each period i, compare its fitted constant C_i against the mean and
    RMS of the fitted constants from all other periods:

        z_i = (C_i - mean_others) / RMS_others

    If |z_i| > sigma_threshold, the period is flagged as an outlier.

    If RMS_others is zero, no sigma test is applied for that period.
    """

    outlier_names = set()

    valid = [
        item for item in period_fit_info
        if np.isfinite(item["fit_value"])
    ]

    if len(valid) < 3:
        return outlier_names
    #endif

    for item in valid:
        others = [
            other["fit_value"] for other in valid
            if other["name"] != item["name"]
        ]

        others = np.asarray(others, dtype=float)
        mean_others = float(np.mean(others))
        rms_others = float(np.sqrt(np.mean((others - mean_others) ** 2)))

        if rms_others == 0.0:
            continue
        #endif

        z_value = (item["fit_value"] - mean_others) / rms_others

        item["leave_one_out_mean"] = mean_others
        item["leave_one_out_rms"] = rms_others
        item["leave_one_out_z"] = z_value

        if abs(z_value) > sigma_threshold:
            outlier_names.add(item["name"])
        #endif
    #endfor

    return outlier_names


def main():
    periods = read_period_charge_file(INPUT_CSV)

    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)

    write_sorted_period_charge_file(
        periods,
        OUTPUT_SORTED_CSV,
    )

    period_fit_info = []

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

        finite_mask = np.isfinite(percent_differences)

        if not np.any(finite_mask):
            period_fit_info.append(
                {
                    "name": period_name,
                    "fit_value": np.nan,
                    "fit_err": np.nan,
                    "fit_rms": np.nan,
                    "n_fit": 0,
                    "rows": rows,
                    "runnums": runnums,
                    "percent_differences": percent_differences,
                    "finite_mask": finite_mask,
                }
            )
            continue
        #endif

        fit_value, fit_err, fit_rms, n_fit = constant_fit(percent_differences)

        period_fit_info.append(
            {
                "name": period_name,
                "fit_value": fit_value,
                "fit_err": fit_err,
                "fit_rms": fit_rms,
                "n_fit": n_fit,
                "rows": rows,
                "runnums": runnums,
                "percent_differences": percent_differences,
                "finite_mask": finite_mask,
            }
        )
    #endfor

    outlier_period_names = find_period_outliers(
        period_fit_info,
        OUTLIER_SIGMA_THRESHOLD,
    )

    clean_period_fit_values = [
        item["fit_value"] for item in period_fit_info
        if np.isfinite(item["fit_value"])
        and item["name"] not in outlier_period_names
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
    print("Constant fits to percent difference:")
    print("  percent difference = 100 * (sum(HEL::Scaler) - RUN::Scaler) / RUN::Scaler")
    print("  sum(HEL::Scaler) = HEL+ + HEL- + HEL(unassigned)")
    print("-" * 112)
    print(
        f"{'Period':<20} "
        f"{'N':>5} "
        f"{'constant [%]':>16} "
        f"{'err [%]':>14} "
        f"{'RMS [%]':>14} "
        f"{'period-z':>12} "
        f"{'used in avg?':>14}"
    )
    print("-" * 112)

    for item in period_fit_info:
        period_z = item.get("leave_one_out_z", np.nan)

        if item["name"] in outlier_period_names:
            used_string = "NO"
        elif np.isfinite(item["fit_value"]):
            used_string = "YES"
        else:
            used_string = "NO"
        #endif

        print(
            f"{item['name']:<20} "
            f"{item['n_fit']:>5d} "
            f"{item['fit_value']:>16.7f} "
            f"{item['fit_err']:>14.7f} "
            f"{item['fit_rms']:>14.7f} "
            f"{period_z:>12.3f} "
            f"{used_string:>14}"
        )
    #endfor

    print("-" * 112)

    if len(outlier_period_names) > 0:
        print()
        print(
            f"Period-level outliers removed from overall average "
            f"using |z| > {OUTLIER_SIGMA_THRESHOLD:.1f}:"
        )

        for item in period_fit_info:
            if item["name"] not in outlier_period_names:
                continue
            #endif

            print(
                f"  {item['name']}: "
                f"C = {item['fit_value']:.7f}%, "
                f"other-period mean = {item.get('leave_one_out_mean', np.nan):.7f}%, "
                f"other-period RMS = {item.get('leave_one_out_rms', np.nan):.7f}%, "
                f"z = {item.get('leave_one_out_z', np.nan):.3f}"
            )
        #endfor
    else:
        print()
        print(
            f"No period-level outliers found using "
            f"|z| > {OUTLIER_SIGMA_THRESHOLD:.1f}."
        )
    #endif

    print()
    print("Overall average of period constants, after period-outlier removal:")
    print(f"  average = {overall_period_average:.7f}%")
    print(f"  RMS     = {overall_period_rms:.7f}%")
    print()

    fig, ax = plt.subplots(figsize=(18, 7))

    all_runs = []
    all_percent_differences = []

    for i_period, item in enumerate(period_fit_info):
        period_name = item["name"]
        runnums = item["runnums"]
        percent_differences = item["percent_differences"]
        finite_mask = item["finite_mask"]

        if not np.any(finite_mask):
            continue
        #endif

        fit_value = item["fit_value"]

        ax.plot(
            runnums[finite_mask],
            percent_differences[finite_mask],
            marker="o",
            linestyle="none",
            markersize=4,
            label=period_name,
        )

        xmin = float(np.min(runnums))
        xmax = float(np.max(runnums))

        if period_name in outlier_period_names:
            fit_linestyle = "--"
            fit_linewidth = 1.0
        else:
            fit_linestyle = "-"
            fit_linewidth = 1.2
        #endif

        ax.hlines(
            fit_value,
            xmin,
            xmax,
            linewidth=fit_linewidth,
            linestyle=fit_linestyle,
        )

        text_x = 0.5 * (xmin + xmax)

        if i_period % 2 == 0:
            text_y = -18.8
        else:
            text_y = -16.3
        #endif

        if period_name in outlier_period_names:
            label_text = f"{period_name}: C = {fit_value:.3f}% OUT"
        else:
            label_text = f"{period_name}: C = {fit_value:.3f}%"
        #endif

        ax.text(
            text_x,
            text_y,
            label_text,
            ha="center",
            va="bottom",
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

        all_runs.extend(runnums[finite_mask].tolist())
        all_percent_differences.extend(percent_differences[finite_mask].tolist())

        item["sorted_rows"] = item["rows"]
        item["max_runnum"] = int(np.max(runnums))
    #endfor

    sorted_periods_for_boundaries = [
        item for item in period_fit_info
        if "max_runnum" in item
    ]

    for i_period, item in enumerate(sorted_periods_for_boundaries[:-1]):
        this_max = item["max_runnum"]

        next_rows = sorted_periods_for_boundaries[i_period + 1]["sorted_rows"]
        next_min = min(row["runnum"] for row in next_rows)

        boundary_x = 0.5 * (this_max + next_min)

        ax.axvline(
            boundary_x,
            color="gray",
            linewidth=0.7,
            alpha=0.6,
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

    if len(all_runs) > 0:
        xmin_all = min(all_runs)
        xmax_all = max(all_runs)
        dx = xmax_all - xmin_all

        if dx == 0.0:
            ax.set_xlim(xmin_all - 1.0, xmax_all + 1.0)
        else:
            ax.set_xlim(xmin_all - 0.02 * dx, xmax_all + 0.02 * dx)
        #endif
    #endif

    ax.set_title("Non-QA-gated RUN::Scaler vs HEL::Scaler charge comparison by run period")

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