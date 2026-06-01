#!/usr/bin/env python3

import argparse
import os
import re
import numpy as np
import matplotlib.pyplot as plt


def parse_fit_file(path):
    """
    Parse a Mathematica-style text file containing entries like:

    enpiLowxBGEchi2FitsALUsinphi = {{x, y, dy}, {x, y, dy}, ...};

    Returns
    -------
    data : dict
        data[name] = numpy array with shape (N, 3)
        columns are x, y, dy.
    """

    with open(path, "r") as f:
        text = f.read()

    pattern = re.compile(r"([A-Za-z0-9_]+)\s*=\s*\{\{(.*?)\}\}\s*;", re.DOTALL)

    data = {}

    for match in pattern.finditer(text):
        name = match.group(1)
        body = match.group(2)

        rows = re.findall(r"\{([^{}]+)\}", "{{" + body + "}}")
        parsed_rows = []

        for row in rows:
            parts = [p.strip() for p in row.split(",")]

            if len(parts) != 3:
                raise ValueError(f"Malformed row in {path} for {name}: {row}")
            #endif

            parsed_rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
        #endfor

        data[name] = np.array(parsed_rows, dtype=float)
    #endfor

    if len(data) == 0:
        raise ValueError(f"No fit arrays found in file: {path}")
    #endif

    return data


def get_series(data, xb_tag, modulation):
    """
    Retrieve the array corresponding to an xB tag and modulation.

    Example key:
      enpiLowxBGEchi2FitsALUsinphi
    """

    key = f"enpi{xb_tag}xBGEchi2Fits{modulation}"

    if key not in data:
        available = "\n".join(sorted(data.keys()))
        raise KeyError(
            f"Missing required key:\n"
            f"  {key}\n\n"
            f"Available keys are:\n{available}"
        )
    #endif

    return data[key]


def build_24_point_arrays(data, modulation):
    """
    Build a 24-point array for a given modulation by concatenating the four
    xB bins in the order:
      Low, MidLow, MidHigh, High

    Each xB bin contributes six -t' points.
    """

    xb_tags = ["Low", "MidLow", "MidHigh", "High"]
    arrays = []

    for xb_tag in xb_tags:
        arr = get_series(data, xb_tag, modulation)

        if arr.shape != (6, 3):
            raise ValueError(
                f"Expected 6 rows and 3 columns for {xb_tag}, {modulation}, "
                f"but got shape {arr.shape}"
            )
        #endif

        arrays.append(arr)
    #endfor

    return np.vstack(arrays)


def compute_two_file_pulls(y_fixed, dy_fixed, y_free, dy_free):
    """
    Compute pulls for the two-file comparison mode using independent uncertainties:

      pull = (free - fixed) / sqrt(sigma_free^2 + sigma_fixed^2)
    """

    denom = np.sqrt(dy_fixed * dy_fixed + dy_free * dy_free)
    pulls = np.full_like(y_fixed, np.nan, dtype=float)

    good = denom > 0.0

    pulls[good] = (y_free[good] - y_fixed[good]) / denom[good]

    return pulls


def compute_three_file_pulls(arrays):
    """
    Compute pulls for the three-file run-period comparison mode.

    The point-by-point reference is the simple arithmetic mean of the three
    fitted values:

      mean_i = (y_i,Su22 + y_i,Fa22 + y_i,Sp23) / 3

    The pull for each run period is:

      pull_i,r = (y_i,r - mean_i) / sigma_i,r

    Parameters
    ----------
    arrays : list of numpy arrays
        Each entry has shape (24, 3), with columns x, y, dy.

    Returns
    -------
    pulls_by_period : list of numpy arrays
        One pull array per run period.
    mean_values : numpy array
        The point-by-point mean fitted value.
    """

    y_values = np.vstack([arr[:, 1] for arr in arrays])
    dy_values = np.vstack([arr[:, 2] for arr in arrays])

    mean_values = np.mean(y_values, axis=0)

    pulls_by_period = []

    for i_period in range(len(arrays)):
        pulls = np.full(24, np.nan, dtype=float)
        good = dy_values[i_period, :] > 0.0
        pulls[good] = (y_values[i_period, good] - mean_values[good]) / dy_values[i_period, good]
        pulls_by_period.append(pulls)
    #endfor

    return pulls_by_period, mean_values


def summarize_pulls(pulls):
    """
    Return summary statistics for one pull array.
    """

    finite = pulls[np.isfinite(pulls)]

    if finite.size == 0:
        return {
            "n": 0,
            "mean": np.nan,
            "median": np.nan,
            "mean_abs": np.nan,
            "median_abs": np.nan,
            "rms": np.nan,
            "sigma": np.nan,
        }
    #endif

    return {
        "n": finite.size,
        "mean": np.mean(finite),
        "median": np.median(finite),
        "mean_abs": np.mean(np.abs(finite)),
        "median_abs": np.median(np.abs(finite)),
        "rms": np.sqrt(np.mean(finite * finite)),
        "sigma": np.std(finite, ddof=1) if finite.size > 1 else 0.0,
    }


def print_two_file_pull_summary(modulation, pulls):
    """
    Print mean, median, mean absolute, median absolute, RMS, and sigma of pulls
    for the AUU fixed vs AUU free comparison.
    """

    stats = summarize_pulls(pulls)

    if stats["n"] == 0:
        print(f"{modulation}: no finite pulls")
        return
    #endif

    print("")
    print(f"{modulation}")
    print("-" * len(modulation))
    print(f"  comparison             = AUU free - AUU fixed")
    print(f"  number of points       = {stats['n']:d}")
    print(f"  mean pull              = {stats['mean']:.6f}")
    print(f"  median pull            = {stats['median']:.6f}")
    print(f"  mean |pull|            = {stats['mean_abs']:.6f}")
    print(f"  median |pull|          = {stats['median_abs']:.6f}")
    print(f"  RMS pull               = {stats['rms']:.6f}")
    print(f"  sigma pull             = {stats['sigma']:.6f}")


def print_three_file_pull_summary(modulation, labels, pulls_by_period):
    """
    Print pull summaries for each run period with respect to the point-by-point
    mean of the three run periods.
    """

    print("")
    print(f"{modulation}")
    print("-" * len(modulation))
    print("  comparison             = each run period - point-by-point mean of three")

    for label, pulls in zip(labels, pulls_by_period):
        stats = summarize_pulls(pulls)

        if stats["n"] == 0:
            print(f"  {label}: no finite pulls")
            continue
        #endif

        print(f"  {label}")
        print(f"    number of points     = {stats['n']:d}")
        print(f"    mean pull            = {stats['mean']:.6f}")
        print(f"    median pull          = {stats['median']:.6f}")
        print(f"    mean |pull|          = {stats['mean_abs']:.6f}")
        print(f"    median |pull|        = {stats['median_abs']:.6f}")
        print(f"    RMS pull             = {stats['rms']:.6f}")
        print(f"    sigma pull           = {stats['sigma']:.6f}")
    #endfor


def modulation_label(modulation):
    """
    Convert internal modulation names to plot labels.
    """

    labels = {
        "ALUsinphi": r"$F_{LU}^{\sin\phi}/F_{UU}$",
        "AULsinphi": r"$F_{UL}^{\sin\phi}/F_{UU}$",
        "AULsin2phi": r"$F_{UL}^{\sin 2\phi}/F_{UU}$",
        "ALL": r"$F_{LL}/F_{UU}$",
        "ALLcosphi": r"$F_{LL}^{\cos\phi}/F_{UU}$",
    }

    return labels.get(modulation, modulation)


def add_common_formatting(ax_top, ax_pull):
    """
    Add common axis formatting, xB separators, and xB labels.
    """

    ax_top.grid(True, alpha=0.3)
    ax_top.legend(loc="best", fontsize=11)

    ax_pull.axhline(0.0, linestyle="--", linewidth=1)
    ax_pull.axhline(1.0, linestyle=":", linewidth=1)
    ax_pull.axhline(-1.0, linestyle=":", linewidth=1)
    ax_pull.axhline(2.0, linestyle=":", linewidth=1)
    ax_pull.axhline(-2.0, linestyle=":", linewidth=1)

    ax_pull.set_ylabel("pull", fontsize=13)
    ax_pull.set_xlabel("point index", fontsize=14)
    ax_pull.set_xlim(0.5, 24.5)
    ax_pull.set_xticks(np.arange(1, 25, 1))
    ax_pull.grid(True, alpha=0.3)

    for xline in [6.5, 12.5, 18.5]:
        ax_top.axvline(xline, linestyle="--", linewidth=1, alpha=0.4)
        ax_pull.axvline(xline, linestyle="--", linewidth=1, alpha=0.4)
    #endfor

    xb_labels = [
        (3.5, r"$0.10<x_B<0.25$"),
        (9.5, r"$0.25<x_B<0.35$"),
        (15.5, r"$0.35<x_B<0.45$"),
        (21.5, r"$0.45<x_B<0.60$"),
    ]

    ymin, ymax = ax_top.get_ylim()
    y_text = ymax - 0.06 * (ymax - ymin)

    for xpos, label in xb_labels:
        ax_top.text(
            xpos,
            y_text,
            label,
            ha="center",
            va="top",
            fontsize=10,
        )
    #endfor


def make_two_file_comparison_plot(
    modulation,
    fixed_arr,
    free_arr,
    output_dir,
    output_format,
):
    """
    Make a two-panel comparison plot for one modulation in AUU fixed/free mode.
    """

    point_index = np.arange(1, 25)

    y_fixed = fixed_arr[:, 1]
    dy_fixed = fixed_arr[:, 2]

    y_free = free_arr[:, 1]
    dy_free = free_arr[:, 2]

    pulls = compute_two_file_pulls(y_fixed, dy_fixed, y_free, dy_free)

    fig, (ax_top, ax_pull) = plt.subplots(
        2,
        1,
        figsize=(12.0, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.25], "hspace": 0.05},
    )

    ax_top.errorbar(
        point_index - 0.08,
        y_fixed,
        yerr=dy_fixed,
        fmt="o",
        markersize=4,
        capsize=2,
        linewidth=1,
        label="AUU fixed",
    )

    ax_top.errorbar(
        point_index + 0.08,
        y_free,
        yerr=dy_free,
        fmt="s",
        markersize=4,
        capsize=2,
        linewidth=1,
        label="AUU free",
    )

    ax_top.set_ylabel(modulation_label(modulation), fontsize=14)
    ax_top.set_title(
        f"Exclusive $e p \\rightarrow e n \\pi^+$ comparison: {modulation_label(modulation)}",
        fontsize=15,
    )

    ax_pull.plot(point_index, pulls, "o", markersize=4)

    add_common_formatting(ax_top, ax_pull)

    safe_name = modulation.replace("/", "_").replace(" ", "_")
    output_path = os.path.join(output_dir, f"compare_AUU_fixed_free_{safe_name}.{output_format}")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print_two_file_pull_summary(modulation, pulls)
    print(f"  saved plot             = {output_path}")


def make_three_file_comparison_plot(
    modulation,
    arrays,
    labels,
    output_dir,
    output_format,
):
    """
    Make a two-panel comparison plot for one modulation in three-run-period mode.

    Top panel:
      Su22, Fa22, Sp23 fitted values.

    Bottom panel:
      Pulls with respect to the point-by-point mean of the three fitted values.
    """

    point_index = np.arange(1, 25)

    pulls_by_period, mean_values = compute_three_file_pulls(arrays)

    fig, (ax_top, ax_pull) = plt.subplots(
        2,
        1,
        figsize=(12.0, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.25], "hspace": 0.05},
    )

    x_offsets = [-0.16, 0.00, 0.16]
    marker_styles = ["o", "s", "^"]

    for arr, label, x_offset, marker_style in zip(arrays, labels, x_offsets, marker_styles):
        y = arr[:, 1]
        dy = arr[:, 2]

        ax_top.errorbar(
            point_index + x_offset,
            y,
            yerr=dy,
            fmt=marker_style,
            markersize=4,
            capsize=2,
            linewidth=1,
            label=label,
        )
    #endfor

    for pulls, label, marker_style in zip(pulls_by_period, labels, marker_styles):
        ax_pull.plot(
            point_index,
            pulls,
            marker_style,
            markersize=4,
            linestyle="None",
            label=label,
        )
    #endfor

    ax_top.plot(
        point_index,
        mean_values,
        linestyle="--",
        linewidth=1,
        label="mean",
    )

    ax_top.set_ylabel(modulation_label(modulation), fontsize=14)
    ax_top.set_title(
        f"Exclusive $e p \\rightarrow e n \\pi^+$ run-period comparison: {modulation_label(modulation)}",
        fontsize=15,
    )

    add_common_formatting(ax_top, ax_pull)

    safe_name = modulation.replace("/", "_").replace(" ", "_")
    output_path = os.path.join(output_dir, f"compare_run_periods_{safe_name}.{output_format}")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print_three_file_pull_summary(modulation, labels, pulls_by_period)
    print(f"  saved plot             = {output_path}")


def run_two_file_mode(args):
    """
    Run the original AUU fixed/free comparison mode.
    """

    fixed_data = parse_fit_file(args.fit_files[0])
    free_data = parse_fit_file(args.fit_files[1])

    modulations = [
        "ALUsinphi",
        "AULsinphi",
        "AULsin2phi",
        "ALL",
        "ALLcosphi",
    ]

    for modulation in modulations:
        fixed_arr = build_24_point_arrays(fixed_data, modulation)
        free_arr = build_24_point_arrays(free_data, modulation)

        make_two_file_comparison_plot(
            modulation=modulation,
            fixed_arr=fixed_arr,
            free_arr=free_arr,
            output_dir=args.output_dir,
            output_format=args.format,
        )
    #endfor


def run_three_file_mode(args):
    """
    Run the three-run-period comparison mode.
    """

    labels = ["Su22", "Fa22", "Sp23"]
    data_sets = [parse_fit_file(path) for path in args.fit_files]

    modulations = [
        "ALUsinphi",
        "AULsinphi",
        "AULsin2phi",
        "ALL",
        "ALLcosphi",
    ]

    for modulation in modulations:
        arrays = []

        for data in data_sets:
            arrays.append(build_24_point_arrays(data, modulation))
        #endfor

        make_three_file_comparison_plot(
            modulation=modulation,
            arrays=arrays,
            labels=labels,
            output_dir=args.output_dir,
            output_format=args.format,
        )
    #endfor


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare fit-parameter text files. "
            "With two files, the script compares AUU fixed vs AUU free. "
            "With three files, the script compares Su22, Fa22, and Sp23."
        )
    )

    parser.add_argument(
        "fit_files",
        nargs="+",
        help=(
            "Input fit text files. "
            "Pass two files for AUU fixed/free mode, or three files for Su22/Fa22/Sp23 mode."
        ),
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        default="output/UU_comparison",
        help="Directory where comparison plots will be written.",
    )

    parser.add_argument(
        "--format",
        default="png",
        choices=["png", "pdf"],
        help="Output plot format.",
    )

    args = parser.parse_args()

    if len(args.fit_files) not in [2, 3]:
        raise ValueError(
            "This script requires either exactly two fit files "
            "for AUU fixed/free mode or exactly three fit files "
            "for Su22/Fa22/Sp23 mode."
        )
    #endif

    os.makedirs(args.output_dir, exist_ok=True)

    if len(args.fit_files) == 2:
        run_two_file_mode(args)
    elif len(args.fit_files) == 3:
        run_three_file_mode(args)
    #endif


if __name__ == "__main__":
    main()
#endif