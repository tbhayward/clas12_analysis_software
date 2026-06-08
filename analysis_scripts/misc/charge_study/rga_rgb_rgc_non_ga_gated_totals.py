#!/usr/bin/env python3

import os
import csv
import numpy as np
import matplotlib.pyplot as plt


INPUT_CSV = "non_qa_gated_totals.csv"
OUTPUT_PNG = "output/rga_rgb_rgc_non_gated_totals.png"


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

            if hel_sum == 0.0:
                ratio = np.nan
            else:
                ratio = run_scaler / hel_sum
            #endif

            current_period["rows"].append(
                {
                    "runnum": runnum,
                    "run_scaler": run_scaler,
                    "hel_pos": hel_pos,
                    "hel_neg": hel_neg,
                    "hel_unassigned": hel_unassigned,
                    "hel_sum": hel_sum,
                    "ratio": ratio,
                }
            )
        #endfor
    #endwith

    return periods


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

    fig, ax = plt.subplots(figsize=(16, 7))

    print()
    print("Constant fits to RUN::Scaler / sum(HEL::Scaler):")
    print("-" * 72)
    print(f"{'Period':<20} {'N':>5} {'constant':>14} {'err':>14} {'RMS':>14}")
    print("-" * 72)

    all_runs = []
    all_ratios = []

    for period in periods:
        period_name = period["name"]
        rows = sorted(period["rows"], key=lambda row: row["runnum"])

        if len(rows) == 0:
            continue
        #endif

        runnums = np.array([row["runnum"] for row in rows], dtype=float)
        ratios = np.array([row["ratio"] for row in rows], dtype=float)

        finite_mask = np.isfinite(ratios)

        if not np.any(finite_mask):
            print(f"{period_name:<20} {'0':>5} {'nan':>14} {'nan':>14} {'nan':>14}")
            continue
        #endif

        fit_value, fit_err, fit_rms, n_fit = constant_fit(ratios)

        ax.plot(
            runnums[finite_mask],
            ratios[finite_mask],
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
        text_y = fit_value + 0.012

        if text_y > 1.17:
            text_y = fit_value - 0.025
        #endif

        ax.text(
            text_x,
            text_y,
            f"{period_name}\nC = {fit_value:.5f}",
            ha="center",
            va="bottom",
            fontsize=8,
            rotation=90,
        )

        print(
            f"{period_name:<20} "
            f"{n_fit:>5d} "
            f"{fit_value:>14.7f} "
            f"{fit_err:>14.7f} "
            f"{fit_rms:>14.7f}"
        )

        all_runs.extend(runnums[finite_mask].tolist())
        all_ratios.extend(ratios[finite_mask].tolist())

        period["sorted_rows"] = rows
        period["max_runnum"] = int(np.max(runnums))
    #endfor

    print("-" * 72)
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
        1.0,
        color="black",
        linewidth=0.8,
        linestyle="--",
    )

    ax.set_xlabel("Run Number")
    ax.set_ylabel(
        r"RUN::Scaler charge / "
        r"$\left(\mathrm{HEL}^{+} + \mathrm{HEL}^{-} + \mathrm{HEL}^{\mathrm{unassigned}}\right)$ charge"
    )

    ax.set_ylim(0.8, 1.2)

    if len(all_runs) > 0:
        xmin_all = min(all_runs)
        xmax_all = max(all_runs)
        dx = xmax_all - xmin_all
        ax.set_xlim(xmin_all - 0.02 * dx, xmax_all + 0.02 * dx)
    #endif

    ax.set_title("Non-QA-gated charge comparison by run period")

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


if __name__ == "__main__":
    main()
#endif