#!/usr/bin/env python3

import os
import csv
import numpy as np
import matplotlib.pyplot as plt


INPUT_CSV = "non_qa_gated_totals.csv"
OUTPUT_PNG = "output/rga_rgb_rgc_non_gated_totals.png"
OUTPUT_SORTED_CSV = "output/non_qa_gated_totals_sorted_by_period.csv"


def read_period_charge_file(path):
    """
    Read a sectioned charge CSV.

    Expected format:

        # RGA Fa18 Inb
        run,RUN::Scaler,HEL::Scaler(+),HEL::Scaler(-),HEL::Scaler(unassigned),0,0
        ...

    Returns
    -------
    periods : list of dict
        Each entry has:
            name : str
            rows : list of dict

    Convention
    ----------
    The plotted quantity is the same convention as the older scaler
    cross-check script:

        percent_difference = 100 * (HEL::Scaler - RUN::Scaler) / RUN::Scaler

    where here:

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

            extra_columns = parts[5:]

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
                    "extra_columns": extra_columns,
                }
            )
        #endfor
    #endwith

    return periods


def write_sorted_period_charge_file(periods, output_path):
    """
    Write a copy of the input CSV with the same period headers and same row
    contents, but with numeric rows sorted by run number inside each period.

    This intentionally preserves the original columns from the input rows.
    It does not add the computed HEL sum or percent difference to the CSV.
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
        Standard error on the mean using the sample RMS / sqrt(N).
    rms : float
        RMS scatter about the mean.
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


def main():
    periods = read_period_charge_file(INPUT_CSV)

    os.makedirs(os.path.dirname(OUTPUT_PNG), exist_ok=True)

    write_sorted_period_charge_file(
        periods,
        OUTPUT_SORTED_CSV,
    )

    fig, ax = plt.subplots(figsize=(16, 7))

    print()
    print("Constant fits to percent difference:")
    print("  percent difference = 100 * (sum(HEL::Scaler) - RUN::Scaler) / RUN::Scaler")
    print("  sum(HEL::Scaler) = HEL+ + HEL- + HEL(unassigned)")
    print("-" * 86)
    print(
        f"{'Period':<20} "
        f"{'N':>5} "
        f"{'constant [%]':>16} "
        f"{'err [%]':>14} "
        f"{'RMS [%]':>14}"
    )
    print("-" * 86)

    all_runs = []
    all_percent_differences = []

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
            print(
                f"{period_name:<20} "
                f"{0:>5d} "
                f"{'nan':>16} "
                f"{'nan':>14} "
                f"{'nan':>14}"
            )
            continue
        #endif

        fit_value, fit_err, fit_rms, n_fit = constant_fit(percent_differences)

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

        ax.hlines(
            fit_value,
            xmin,
            xmax,
            linewidth=1.2,
            linestyle="-",
        )

        text_x = 0.5 * (xmin + xmax)
        text_y = fit_value + 1.0

        if text_y > 18.0:
            text_y = fit_value - 2.5
        #endif

        ax.text(
            text_x,
            text_y,
            f"{period_name}\nC = {fit_value:.3f}%",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )

        print(
            f"{period_name:<20} "
            f"{n_fit:>5d} "
            f"{fit_value:>16.7f} "
            f"{fit_err:>14.7f} "
            f"{fit_rms:>14.7f}"
        )

        all_runs.extend(runnums[finite_mask].tolist())
        all_percent_differences.extend(percent_differences[finite_mask].tolist())

        period["sorted_rows"] = rows
        period["max_runnum"] = int(np.max(runnums))
    #endfor

    print("-" * 86)
    print()

    sorted_periods = [
        period for period in periods
        if "max_runnum" in period
    ]

    for i_period, period in enumerate(sorted_periods[:-1]):
        this_max = period["max_runnum"]

        next_rows = sorted_periods[i_period + 1]["sorted_rows"]
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