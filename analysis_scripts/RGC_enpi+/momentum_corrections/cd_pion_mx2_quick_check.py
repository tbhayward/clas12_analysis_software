#!/usr/bin/env python3

"""
Quick CD-pion Mx2 diagnostic.

This script reads the RGC combined inbending NH3 e pi+ ROOT file and:

1. Selects only events whose pion is reconstructed in the Central Detector:
       detector == 2

2. Reports the percentage of all events satisfying:
       0.70 < Mx2 < 1.00
   that have a CD pion (detector == 2).

3. Reports how the CD subset of that Mx2-selected sample is distributed
   across each requested x interval:
       0.10 < x < 0.25
       0.25 < x < 0.35
       0.35 < x < 0.45
       0.45 < x < 0.60

4. Produces:
       cd_pion_mx2_overall.png
       cd_pion_mx2_by_x.png

The ROOT file can be supplied as a positional command-line argument. If it
is omitted, the configured FILE_PATH is used. The ROOT tree name can be
overridden with --tree-name; PhysicsEvents remains the default.

The ROOT file is streamed in chunks, so the full event sample is never held
in memory at once.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


# ============================================================================
# Configuration
# ============================================================================

FILE_PATH = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
    "paper_versions/rgc_combined_inb_NH3_epi+_mom_corrections.root"
)

TREE_NAME = "PhysicsEvents"

MX2_BRANCH = "Mx2"
X_BRANCH = "x"
PION_DETECTOR_BRANCH = "detector"

# CLAS12 detector code used by the pi+ momentum-correction analysis:
#   1 = Forward Detector
#   2 = Central Detector
CD_DETECTOR_VALUE = 2

MX2_CUT_MIN = 0.70
MX2_CUT_MAX = 1.00

X_BINS = [
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
]

# Full range displayed so the neutron-region peak and surrounding background
# can be judged rather than plotting only the accepted window.
MX2_PLOT_MIN = -0.50
MX2_PLOT_MAX = 2.50
N_MX2_BINS = 150

USE_LOG_Y = True
STEP_SIZE = "250 MB"

OUTPUT_DIRECTORY = Path("cd_pion_quick_check")
OVERALL_OUTPUT = OUTPUT_DIRECTORY / "cd_pion_mx2_overall.png"
X_BINNED_OUTPUT = OUTPUT_DIRECTORY / "cd_pion_mx2_by_x.png"


# ============================================================================
# Main processing
# ============================================================================

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Check the CD-pion fraction inside the configured Mx2 window."
        )
    )
    parser.add_argument(
        "input_root",
        nargs="?",
        type=Path,
        default=FILE_PATH,
        help=(
            "Input ROOT file. Defaults to "
            f"{FILE_PATH}"
        ),
    )
    parser.add_argument(
        "--tree-name",
        default=TREE_NAME,
        help=f"ROOT tree name (default: {TREE_NAME}).",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    input_root = arguments.input_root
    tree_name = arguments.tree_name

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)

    histogram_edges = np.linspace(
        MX2_PLOT_MIN,
        MX2_PLOT_MAX,
        N_MX2_BINS + 1,
    )

    overall_histogram = np.zeros(N_MX2_BINS, dtype=np.int64)

    x_histograms = np.zeros(
        (len(X_BINS), N_MX2_BINS),
        dtype=np.int64,
    )

    # Overall CD counts.
    n_all = 0
    n_all_mx2_selected = 0
    n_cd = 0
    n_cd_mx2_selected = 0

    # Counts in each x bin.
    n_all_mx2_selected_in_x_bin = np.zeros(
        len(X_BINS),
        dtype=np.int64,
    )
    n_cd_in_x_bin = np.zeros(len(X_BINS), dtype=np.int64)
    n_cd_mx2_selected_in_x_bin = np.zeros(
        len(X_BINS),
        dtype=np.int64,
    )

    with uproot.open(input_root) as root_file:
        if tree_name not in root_file:
            available = "\n".join(f"  {key}" for key in root_file.keys())
            raise KeyError(
                f"Tree '{tree_name}' was not found in:\n"
                f"{input_root}\n\n"
                f"Available ROOT objects:\n{available}"
            )

        tree = root_file[tree_name]

        required_branches = [
            MX2_BRANCH,
            X_BRANCH,
            PION_DETECTOR_BRANCH,
        ]

        available_branches = set(tree.keys())
        missing = [
            branch
            for branch in required_branches
            if branch not in available_branches
        ]

        if missing:
            branch_listing = "\n".join(
                f"  {branch}" for branch in sorted(available_branches)
            )
            raise KeyError(
                "Missing required branch or branches:\n"
                f"  {', '.join(missing)}\n\n"
                f"Available branches:\n{branch_listing}"
            )

        print(f"File: {input_root}")
        print(f"Tree: {tree_name}")
        print(f"Entries in tree: {tree.num_entries:,}")
        print(
            "CD-pion requirement: "
            f"{PION_DETECTOR_BRANCH} == {CD_DETECTOR_VALUE}"
        )
        print()

        for arrays in tree.iterate(
            expressions=required_branches,
            step_size=STEP_SIZE,
            library="np",
        ):
            mx2 = np.asarray(arrays[MX2_BRANCH])
            x = np.asarray(arrays[X_BRANCH])
            detector = np.asarray(arrays[PION_DETECTOR_BRANCH])

            finite_mx2 = np.isfinite(mx2)
            finite_x = np.isfinite(x)

            is_cd = detector == CD_DETECTOR_VALUE

            all_base = finite_mx2
            all_mx2_selected = (
                all_base
                & (mx2 > MX2_CUT_MIN)
                & (mx2 < MX2_CUT_MAX)
            )

            cd_base = is_cd & finite_mx2
            cd_mx2_selected = all_mx2_selected & is_cd

            n_all += np.count_nonzero(all_base)
            n_all_mx2_selected += np.count_nonzero(all_mx2_selected)
            n_cd += np.count_nonzero(cd_base)
            n_cd_mx2_selected += np.count_nonzero(cd_mx2_selected)

            chunk_overall_histogram, _ = np.histogram(
                mx2[cd_base],
                bins=histogram_edges,
            )
            overall_histogram += chunk_overall_histogram

            for index, (x_low, x_high) in enumerate(X_BINS):
                in_x_bin = (
                    finite_x
                    & (x > x_low)
                    & (x < x_high)
                )

                all_selected_in_x_bin = (
                    all_mx2_selected
                    & in_x_bin
                )

                cd_in_x_bin = (
                    cd_base
                    & in_x_bin
                )

                cd_selected_in_x_bin = (
                    cd_mx2_selected
                    & in_x_bin
                )

                n_all_mx2_selected_in_x_bin[index] += np.count_nonzero(
                    all_selected_in_x_bin
                )
                n_cd_in_x_bin[index] += np.count_nonzero(cd_in_x_bin)
                n_cd_mx2_selected_in_x_bin[index] += np.count_nonzero(
                    cd_selected_in_x_bin
                )

                chunk_x_histogram, _ = np.histogram(
                    mx2[cd_in_x_bin],
                    bins=histogram_edges,
                )
                x_histograms[index] += chunk_x_histogram

    print_summary(
        n_all=n_all,
        n_all_mx2_selected=n_all_mx2_selected,
        n_cd=n_cd,
        n_cd_mx2_selected=n_cd_mx2_selected,
        n_all_mx2_selected_in_x_bin=n_all_mx2_selected_in_x_bin,
        n_cd_in_x_bin=n_cd_in_x_bin,
        n_cd_mx2_selected_in_x_bin=n_cd_mx2_selected_in_x_bin,
    )

    make_overall_plot(
        histogram_edges=histogram_edges,
        histogram_counts=overall_histogram,
        n_all_mx2_selected=n_all_mx2_selected,
        n_cd=n_cd,
        n_cd_mx2_selected=n_cd_mx2_selected,
    )

    make_x_binned_plot(
        histogram_edges=histogram_edges,
        histogram_counts=x_histograms,
        n_all_mx2_selected_in_x_bin=n_all_mx2_selected_in_x_bin,
        n_cd_in_x_bin=n_cd_in_x_bin,
        n_cd_mx2_selected_in_x_bin=n_cd_mx2_selected_in_x_bin,
    )


# ============================================================================
# Printed results
# ============================================================================

def print_summary(
    n_all: int,
    n_all_mx2_selected: int,
    n_cd: int,
    n_cd_mx2_selected: int,
    n_all_mx2_selected_in_x_bin: np.ndarray,
    n_cd_in_x_bin: np.ndarray,
    n_cd_mx2_selected_in_x_bin: np.ndarray,
) -> None:
    print("Mx2-selected sample and CD fraction")
    print("-----------------------------------")
    print(f"All events with finite Mx2:                 {n_all:,}")
    print(
        f"All events with {MX2_CUT_MIN:.2f} < Mx2 "
        f"< {MX2_CUT_MAX:.2f}:       {n_all_mx2_selected:,}"
    )
    print(f"CD-pion events with finite Mx2:             {n_cd:,}")
    print(
        f"CD-pion events with {MX2_CUT_MIN:.2f} < Mx2 "
        f"< {MX2_CUT_MAX:.2f}: {n_cd_mx2_selected:,}"
    )

    if n_all_mx2_selected > 0:
        cd_fraction_of_window = (
            100.0 * n_cd_mx2_selected / n_all_mx2_selected
        )
        print(
            "Percentage of all Mx2-selected events that are CD: "
            f"{cd_fraction_of_window:.6f}%"
        )
    else:
        print(
            "Percentage of all Mx2-selected events that are CD: undefined"
        )

    print()
    print("CD fraction of the Mx2-selected sample in each x bin")
    print("----------------------------------------------------")

    for (
        x_low,
        x_high,
    ), n_all_selected_x, n_cd_total_x, n_cd_selected_x in zip(
        X_BINS,
        n_all_mx2_selected_in_x_bin,
        n_cd_in_x_bin,
        n_cd_mx2_selected_in_x_bin,
    ):
        cd_fraction_in_selected_x_bin = (
            100.0 * n_cd_selected_x / n_all_selected_x
            if n_all_selected_x > 0
            else np.nan
        )

        cd_mx2_pass_fraction = (
            100.0 * n_cd_selected_x / n_cd_total_x
            if n_cd_total_x > 0
            else np.nan
        )

        print(
            f"{x_low:.2f} < x < {x_high:.2f}: "
            f"{n_cd_selected_x:,} CD / {n_all_selected_x:,} all "
            f"Mx2-selected events; "
            f"{cd_fraction_in_selected_x_bin:.6f}% are CD; "
            f"{cd_mx2_pass_fraction:.6f}% of CD events in this x bin "
            f"pass the Mx2 window"
        )

    n_all_selected_in_listed_bins = int(
        np.sum(n_all_mx2_selected_in_x_bin)
    )
    n_cd_selected_in_listed_bins = int(
        np.sum(n_cd_mx2_selected_in_x_bin)
    )

    listed_bins_cd_fraction = (
        100.0
        * n_cd_selected_in_listed_bins
        / n_all_selected_in_listed_bins
        if n_all_selected_in_listed_bins > 0
        else np.nan
    )

    print()
    print(
        "All Mx2-selected events in the four listed x bins: "
        f"{n_all_selected_in_listed_bins:,}"
    )
    print(
        "CD Mx2-selected events in the four listed x bins: "
        f"{n_cd_selected_in_listed_bins:,}"
    )
    print(
        "Percentage of Mx2-selected events in those bins that are CD: "
        f"{listed_bins_cd_fraction:.6f}%"
    )


# ============================================================================
# Plotting
# ============================================================================

def configure_mx2_axis(axis: plt.Axes) -> None:
    axis.axvspan(
        MX2_CUT_MIN,
        MX2_CUT_MAX,
        alpha=0.20,
        label=(
            rf"${MX2_CUT_MIN:.2f}<M_x^2"
            rf"<{MX2_CUT_MAX:.2f}$"
        ),
    )

    axis.axvline(
        MX2_CUT_MIN,
        linestyle="--",
        linewidth=1.5,
    )
    axis.axvline(
        MX2_CUT_MAX,
        linestyle="--",
        linewidth=1.5,
    )

    # Neutron mass squared, useful as a visual reference for e pi+ missing mass.
    neutron_mass_gev = 0.939565
    neutron_mass_squared = neutron_mass_gev**2

    axis.axvline(
        neutron_mass_squared,
        linestyle=":",
        linewidth=1.5,
        label=rf"$m_n^2={neutron_mass_squared:.3f}$ GeV$^2$",
    )

    axis.set_xlim(MX2_PLOT_MIN, MX2_PLOT_MAX)
    axis.set_xlabel(r"$M_x^2(e\pi^+)$ (GeV$^2$)")
    axis.set_ylabel("Counts")
    axis.grid(axis="y", alpha=0.25)

    if USE_LOG_Y:
        axis.set_yscale("log")
        axis.set_ylim(bottom=0.8)


def make_overall_plot(
    histogram_edges: np.ndarray,
    histogram_counts: np.ndarray,
    n_all_mx2_selected: int,
    n_cd: int,
    n_cd_mx2_selected: int,
) -> None:
    centers = 0.5 * (
        histogram_edges[:-1] + histogram_edges[1:]
    )
    widths = np.diff(histogram_edges)

    figure, axis = plt.subplots(figsize=(10, 7))

    axis.bar(
        centers,
        histogram_counts,
        width=widths,
        align="center",
        linewidth=0,
    )

    configure_mx2_axis(axis)

    percentage = (
        100.0 * n_cd_mx2_selected / n_all_mx2_selected
        if n_all_mx2_selected > 0
        else np.nan
    )

    axis.text(
        0.97,
        0.95,
        (
            f"All events inside Mx2 window: {n_all_mx2_selected:,}\n"
            f"CD events inside Mx2 window: {n_cd_mx2_selected:,}\n"
            f"CD fraction of Mx2 window: {percentage:.3f}%"
        ),
        transform=axis.transAxes,
        horizontalalignment="right",
        verticalalignment="top",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "black",
            "alpha": 0.88,
        },
    )

    axis.set_title(
        r"CD-pion $M_x^2(e\pi^+)$ distribution"
    )
    axis.legend(loc="lower right")

    figure.tight_layout()
    figure.savefig(
        OVERALL_OUTPUT,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)

    print()
    print(f"Saved overall CD plot to: {OVERALL_OUTPUT}")


def make_x_binned_plot(
    histogram_edges: np.ndarray,
    histogram_counts: np.ndarray,
    n_all_mx2_selected_in_x_bin: np.ndarray,
    n_cd_in_x_bin: np.ndarray,
    n_cd_mx2_selected_in_x_bin: np.ndarray,
) -> None:
    centers = 0.5 * (
        histogram_edges[:-1] + histogram_edges[1:]
    )
    widths = np.diff(histogram_edges)

    figure, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14, 10),
        sharex=True,
    )

    for (
        axis,
        (x_low, x_high),
        counts,
        n_all_selected_x,
        n_cd_total_x,
        n_cd_selected_x,
    ) in zip(
        axes.ravel(),
        X_BINS,
        histogram_counts,
        n_all_mx2_selected_in_x_bin,
        n_cd_in_x_bin,
        n_cd_mx2_selected_in_x_bin,
    ):
        axis.bar(
            centers,
            counts,
            width=widths,
            align="center",
            linewidth=0,
        )

        configure_mx2_axis(axis)

        cd_fraction_in_selected_x_bin = (
            100.0 * n_cd_selected_x / n_all_selected_x
            if n_all_selected_x > 0
            else np.nan
        )

        cd_mx2_pass_fraction = (
            100.0 * n_cd_selected_x / n_cd_total_x
            if n_cd_total_x > 0
            else np.nan
        )

        axis.text(
            0.97,
            0.95,
            (
                f"All events in Mx2 window: {n_all_selected_x:,}\n"
                f"CD events in Mx2 window: {n_cd_selected_x:,}\n"
                f"CD fraction: {cd_fraction_in_selected_x_bin:.3f}%\n"
                f"CD Mx2 pass fraction: {cd_mx2_pass_fraction:.3f}%"
            ),
            transform=axis.transAxes,
            horizontalalignment="right",
            verticalalignment="top",
            bbox={
                "boxstyle": "round",
                "facecolor": "white",
                "edgecolor": "black",
                "alpha": 0.88,
            },
        )

        axis.set_title(
            rf"${x_low:.2f}<x<{x_high:.2f}$"
        )
        axis.legend(loc="lower right")

    figure.suptitle(
        r"CD-only $M_x^2(e\pi^+)$ distributions in each $x$ interval",
        fontsize=16,
    )

    figure.tight_layout(
        rect=(0.0, 0.0, 1.0, 0.96),
    )

    figure.savefig(
        X_BINNED_OUTPUT,
        dpi=200,
        bbox_inches="tight",
    )
    plt.close(figure)

    print(f"Saved CD x-binned canvas to: {X_BINNED_OUTPUT}")


if __name__ == "__main__":
    main()
