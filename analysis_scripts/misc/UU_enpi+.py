#!/usr/bin/env python3

import argparse
import os
import re
import math
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


def compute_pulls(y_fixed, dy_fixed, y_free, dy_free):
    """
    Compute pulls using independent uncertainties:

      pull = (free - fixed) / sqrt(sigma_free^2 + sigma_fixed^2)
    """

    denom = np.sqrt(dy_fixed * dy_fixed + dy_free * dy_free)
    pulls = np.full_like(y_fixed, np.nan, dtype=float)

    good = denom > 0.0

    pulls[good] = (y_free[good] - y_fixed[good]) / denom[good]

    return pulls


def print_pull_summary(modulation, pulls):
    """
    Print mean, median, mean absolute, median absolute, and RMS pull.
    """

    finite = pulls[np.isfinite(pulls)]

    if finite.size == 0:
        print(f"{modulation}: no finite pulls")
        return
    #endif

    mean_pull = np.mean(finite)
    median_pull = np.median(finite)
    mean_abs_pull = np.mean(np.abs(finite))
    median_abs_pull = np.median(np.abs(finite))
    rms_pull = np.sqrt(np.mean(finite * finite))
    sigma_pull = np.std(finite, ddof=1) if finite.size > 1 else 0.0

    print("")
    print(f"{modulation}")
    print("-" * len(modulation))
    print(f"  number of points       = {finite.size:d}")
    print(f"  mean pull              = {mean_pull:.6f}")
    print(f"  median pull            = {median_pull:.6f}")
    print(f"  mean |pull|            = {mean_abs_pull:.6f}")
    print(f"  median |pull|          = {median_abs_pull:.6f}")
    print(f"  RMS pull               = {rms_pull:.6f}")
    print(f"  sigma pull             = {sigma_pull:.6f}")


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


def make_comparison_plot(
    modulation,
    fixed_arr,
    free_arr,
    output_dir,
    output_format,
):
    """
    Make a two-panel comparison plot for one modulation.
    """

    point_index = np.arange(1, 25)

    y_fixed = fixed_arr[:, 1]
    dy_fixed = fixed_arr[:, 2]

    y_free = free_arr[:, 1]
    dy_free = free_arr[:, 2]

    pulls = compute_pulls(y_fixed, dy_fixed, y_free, dy_free)

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
    ax_top.set_title(f"Exclusive $e p \\rightarrow e n \\pi^+$ comparison: {modulation_label(modulation)}", fontsize=15)
    ax_top.grid(True, alpha=0.3)
    ax_top.legend(loc="best", fontsize=11)

    ax_pull.axhline(0.0, linestyle="--", linewidth=1)
    ax_pull.axhline(1.0, linestyle=":", linewidth=1)
    ax_pull.axhline(-1.0, linestyle=":", linewidth=1)
    ax_pull.axhline(2.0, linestyle=":", linewidth=1)
    ax_pull.axhline(-2.0, linestyle=":", linewidth=1)

    ax_pull.plot(point_index, pulls, "o", markersize=4)

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

    safe_name = modulation.replace("/", "_").replace(" ", "_")
    output_path = os.path.join(output_dir, f"compare_AUU_fixed_free_{safe_name}.{output_format}")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print_pull_summary(modulation, pulls)
    print(f"  saved plot             = {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compare fit-parameter text files with AUU fixed vs AUU free. "
            "Each input file should contain Mathematica-style arrays."
        )
    )

    parser.add_argument(
        "auu_fixed_file",
        help="Text file for the AUU-fixed fit results.",
    )

    parser.add_argument(
        "auu_free_file",
        help="Text file for the AUU-free fit results.",
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

    os.makedirs(args.output_dir, exist_ok=True)

    fixed_data = parse_fit_file(args.auu_fixed_file)
    free_data = parse_fit_file(args.auu_free_file)

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

        make_comparison_plot(
            modulation=modulation,
            fixed_arr=fixed_arr,
            free_arr=free_arr,
            output_dir=args.output_dir,
            output_format=args.format,
        )
    #endfor


if __name__ == "__main__":
    main()
#endif