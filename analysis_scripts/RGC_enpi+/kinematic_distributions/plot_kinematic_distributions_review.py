#!/usr/bin/env python3
"""Publication-style RGC e n pi+ kinematic-distribution plots for reviewer response.

Reads nominal NH3 paper-version ROOT trees for Su22/Fa22/Sp23 and applies the
period- and (xB,-t')-dependent nominal +/-2 sigma Mx2 windows from channel
selection. No asymmetry-extraction quantities are fitted or modified.

In addition to the existing kinematic plots, this script produces pion-momentum
distributions in the final 24 (xB,-t') analysis bins and exports period-resolved
0.25-GeV momentum histograms.  The latter are designed to be combined directly
with the outputs of the beta_pid_2d_fit particle-misidentification study.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot

from matplotlib.colors import LogNorm


# =============================================================================
# Analysis configuration
# =============================================================================

PERIODS = (
    "su22",
    "fa22",
    "sp23",
)


XB_BINS = (
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
)


TP_BINS = (
    (0.05, 0.25),
    (0.25, 0.45),
    (0.45, 0.65),
    (0.65, 0.85),
    (0.85, 1.05),
    (1.05, 1.25),
)


PAPER_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions"
)


DEFAULT_INPUTS = {
    period:
        PAPER_DIR
        / f"rgc_{period}_inb_NH3_epi+_mom_corrections.root"

    for period in PERIODS
}


DEFAULT_CUTS = Path(
    "../channel_selection/output/"
    "channel_selection_mx2_fit_stability/"
    "final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts.json"
)


DEFAULT_OUT = Path(
    "output/kinematic_distributions_review"
)


TREE = "PhysicsEvents"

CHUNK = "250 MB"


CMAP = "turbo"


# =============================================================================
# Pion-momentum binning
#
# Deliberately identical to the portable beta_pid_2d_fit output:
#
#     0.50--5.00 GeV in 0.25-GeV intervals.
#
# This lets the two CSVs be joined without interpolation.
# =============================================================================

PION_P_MIN = 0.50
PION_P_MAX = 5.00
PION_P_WIDTH = 0.25


PION_P_EDGES = np.arange(
    PION_P_MIN,
    PION_P_MAX + PION_P_WIDTH,
    PION_P_WIDTH,
)


# =============================================================================
# ROOT branch aliases
# =============================================================================

ALIASES = {
    "xB": (
        "xB",
        "x",
        "xb",
        "x_b",
    ),

    "Q2": (
        "Q2",
        "q2",
    ),

    "tprime": (
        "tprime",
        "t_prime",
        "tp",
        "tPrime",
    ),

    "Mx2": (
        "Mx2",
        "mx2",
        "Mx2_epi",
        "Mx2_epip",
        "missing_mass_squared",
        "missing_mass2",
    ),

    "phi": (
        "phi",
        "phi1",
        "phi_h",
        "trento_phi",
    ),

    "W": (
        "W",
        "w",
    ),

    "y": (
        "y",
        "inelasticity",
    ),

    # Pion momentum.
    "p_p": (
        "p_p",
        "p_pi",
        "p_pip",
        "pion_p",
    ),
}


# =============================================================================
# ROOT helpers
# =============================================================================

def resolve(
    tree,
    logical,
    required=True,
):

    names = set(
        tree.keys()
    )


    for candidate in ALIASES[
        logical
    ]:

        if candidate in names:
            return candidate


    if required:

        raise KeyError(
            f"Could not resolve {logical}; "
            f"tried {ALIASES[logical]}"
        )


    return None


def phi_degrees(
    values,
    branch,
):

    values = np.asarray(
        values,
        dtype=float,
    )


    finite_values = values[
        np.isfinite(
            values
        )
    ]


    if not finite_values.size:
        return values


    if (
        np.nanmax(
            np.abs(
                finite_values
            )
        )
        <= 2 * np.pi + 0.25
    ):

        values = np.degrees(
            values
        )


    return np.mod(
        values,
        360.0,
    )


# =============================================================================
# Analysis-bin helpers
# =============================================================================

def bin_index(
    xb,
    minus_tprime,
):

    output = np.full(
        xb.shape,
        -1,
        dtype=np.int16,
    )


    for ix, (
        xb_low,
        xb_high,
    ) in enumerate(
        XB_BINS
    ):

        for it, (
            tp_low,
            tp_high,
        ) in enumerate(
            TP_BINS
        ):

            analysis_bin = (
                ix * len(
                    TP_BINS
                )
                + it
                + 1
            )


            mask = (
                (xb >= xb_low)
                & (xb < xb_high)
                & (minus_tprime >= tp_low)
                & (minus_tprime < tp_high)
            )


            output[
                mask
            ] = analysis_bin


    return output


def load_nominal_windows(
    path,
):

    payload = json.loads(
        path.read_text()
    )


    if (
        float(
            payload.get(
                "nominal_sigma_multiple",
                np.nan,
            )
        )
        != 2.0
    ):

        raise RuntimeError(
            f"Expected nominal 2-sigma windows in {path}"
        )


    windows = {}


    for period in PERIODS:

        windows[
            period
        ] = {
            int(
                row[
                    "bin_number"
                ]
            ):
                tuple(
                    map(
                        float,
                        row[
                            "nominal"
                        ],
                    )
                )

            for row in payload[
                "periods"
            ][
                period
            ]
        }


        if (
            len(
                windows[
                    period
                ]
            )
            != 24
        ):

            raise RuntimeError(
                f"Expected 24 nominal windows for {period}"
            )


    return windows


# =============================================================================
# Event collection
# =============================================================================

def collect(
    period,
    path,
    windows,
):

    print(
        f"[load] {period}: {path}",
        flush=True,
    )


    with uproot.open(
        path
    ) as root_file:

        tree = root_file[
            TREE
        ]


        branches = {
            key:
                resolve(
                    tree,
                    key,
                    required=(
                        key not in {
                            "W",
                            "y",
                        }
                    ),
                )

            for key in ALIASES
        }


        print(
            f"[load] {period}: pion momentum branch = "
            f"{branches['p_p']}",
            flush=True,
        )


        expressions = [
            value
            for value in branches.values()
            if value is not None
        ]


        accumulator = {
            key:
                []

            for key in (
                "xB",
                "Q2",
                "minus_tprime",
                "Mx2",
                "phi",
                "W",
                "y",
                "p_p",
                "analysis_bin",
            )
        }


        seen = 0
        selected = 0


        for chunk_index, arrays in enumerate(
            tree.iterate(
                expressions=expressions,
                step_size=CHUNK,
                library="np",
            ),
            start=1,
        ):

            xb = np.asarray(
                arrays[
                    branches[
                        "xB"
                    ]
                ],
                dtype=float,
            )


            q2 = np.asarray(
                arrays[
                    branches[
                        "Q2"
                    ]
                ],
                dtype=float,
            )


            minus_tprime = -np.asarray(
                arrays[
                    branches[
                        "tprime"
                    ]
                ],
                dtype=float,
            )


            mx2 = np.asarray(
                arrays[
                    branches[
                        "Mx2"
                    ]
                ],
                dtype=float,
            )


            phi = phi_degrees(
                arrays[
                    branches[
                        "phi"
                    ]
                ],
                branches[
                    "phi"
                ],
            )


            pion_momentum = np.asarray(
                arrays[
                    branches[
                        "p_p"
                    ]
                ],
                dtype=float,
            )


            if branches[
                "W"
            ]:

                w = np.asarray(
                    arrays[
                        branches[
                            "W"
                        ]
                    ],
                    dtype=float,
                )


            else:

                w = np.full(
                    xb.shape,
                    np.nan,
                )


            if branches[
                "y"
            ]:

                y = np.asarray(
                    arrays[
                        branches[
                            "y"
                        ]
                    ],
                    dtype=float,
                )


            else:

                y = np.full(
                    xb.shape,
                    np.nan,
                )


            analysis_bin = bin_index(
                xb,
                minus_tprime,
            )


            keep = (
                (analysis_bin > 0)

                & np.isfinite(
                    xb
                )

                & np.isfinite(
                    q2
                )

                & np.isfinite(
                    minus_tprime
                )

                & np.isfinite(
                    mx2
                )

                & np.isfinite(
                    phi
                )

                & np.isfinite(
                    pion_momentum
                )
            )


            # ================================================================
            # Final period- and analysis-bin-dependent nominal Mx2 windows.
            # ================================================================

            for analysis_bin_number in range(
                1,
                25,
            ):

                low, high = windows[
                    analysis_bin_number
                ]


                bin_mask = (
                    analysis_bin
                    == analysis_bin_number
                )


                keep[
                    bin_mask
                ] &= (
                    (
                        mx2[
                            bin_mask
                        ]
                        >= low
                    )

                    & (
                        mx2[
                            bin_mask
                        ]
                        <= high
                    )
                )


            # ================================================================
            # DIS cuts.
            # ================================================================

            keep &= (
                q2 > 1.0
            )


            if branches[
                "W"
            ]:

                keep &= (
                    w > 2.0
                )


            if branches[
                "y"
            ]:

                keep &= (
                    y < 0.80
                )


            seen += xb.size


            selected += int(
                np.count_nonzero(
                    keep
                )
            )


            values = {
                "xB":
                    xb,

                "Q2":
                    q2,

                "minus_tprime":
                    minus_tprime,

                "Mx2":
                    mx2,

                "phi":
                    phi,

                "W":
                    w,

                "y":
                    y,

                "p_p":
                    pion_momentum,

                "analysis_bin":
                    analysis_bin,
            }


            for key, values_array in values.items():

                accumulator[
                    key
                ].append(
                    values_array[
                        keep
                    ]
                )


            print(
                f"[load] {period}: "
                f"chunk {chunk_index}; "
                f"seen={seen:,}; "
                f"selected={selected:,}",

                flush=True,
            )


    result = {
        key:
            (
                np.concatenate(
                    value
                )

                if value

                else np.empty(
                    0
                )
            )

        for key, value in accumulator.items()
    }


    result[
        "period"
    ] = np.full(
        result[
            "xB"
        ].size,
        period,
    )


    print(
        f"[load] {period}: DONE; "
        f"selected={result['xB'].size:,}",

        flush=True,
    )


    return result


# =============================================================================
# Generic plotting helpers
# =============================================================================

def finite(
    values,
):

    return values[
        np.isfinite(
            values
        )
    ]


def hist2d(
    ax,
    x,
    y,
    bins,
    ranges,
    xlabel,
    ylabel,
):

    good = (
        np.isfinite(
            x
        )

        & np.isfinite(
            y
        )
    )


    histogram, x_edges, y_edges = np.histogram2d(
        x[
            good
        ],
        y[
            good
        ],
        bins=bins,
        range=ranges,
    )


    positive = histogram[
        histogram > 0
    ]


    norm = LogNorm(
        vmin=(
            max(
                1.0,
                float(
                    positive.min()
                ),
            )

            if positive.size

            else 1.0
        ),

        vmax=max(
            2.0,
            float(
                histogram.max()
            ),
        ),
    )


    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        histogram.T,
        shading="auto",
        norm=norm,
        cmap=CMAP,
    )


    ax.set_xlabel(
        xlabel
    )


    ax.set_ylabel(
        ylabel
    )


    return mesh


def save(
    fig,
    path,
):

    fig.savefig(
        path,
        dpi=220,
        bbox_inches="tight",
    )


    plt.close(
        fig
    )


    print(
        f"[plot] wrote {path}",
        flush=True,
    )


# =============================================================================
# Existing analysis-bin shape plot
# =============================================================================

def plot_variable_by_tprime_with_xb_overlays(
    data,
    key,
    edges,
    xlabel,
    path,
):

    """One row of -t' panels, with four xB bins overlaid."""

    fig, axes = plt.subplots(
        1,
        6,
        figsize=(
            19.0,
            4.1,
        ),
        sharex=True,
        sharey=True,
    )


    centers = 0.5 * (
        edges[:-1]
        + edges[1:]
    )


    for ax, (
        tp_low,
        tp_high,
    ) in zip(
        axes,
        TP_BINS,
    ):

        for xb_low, xb_high in XB_BINS:

            mask = (
                (
                    data[
                        "xB"
                    ]
                    >= xb_low
                )

                & (
                    data[
                        "xB"
                    ]
                    < xb_high
                )

                & (
                    data[
                        "minus_tprime"
                    ]
                    >= tp_low
                )

                & (
                    data[
                        "minus_tprime"
                    ]
                    < tp_high
                )

                & np.isfinite(
                    data[
                        key
                    ]
                )
            )


            values = data[
                key
            ][
                mask
            ]


            counts, _ = np.histogram(
                values,
                bins=edges,
            )


            total = counts.sum()


            density = (
                counts / total

                if total

                else counts.astype(
                    float
                )
            )


            ax.step(
                centers,
                density,
                where="mid",
                linewidth=1.35,
                label=(
                    fr"${xb_low:.2f}\leq x_B<{xb_high:.2f}$ "
                    fr"($N={total:,}$)"
                ),
            )


        ax.set_title(
            fr"${tp_low:.2f}\leq -t'<{tp_high:.2f}$"
        )


        ax.set_xlabel(
            xlabel
        )


        ax.tick_params(
            direction="in",
            top=True,
            right=True,
        )


        ax.legend(
            fontsize=6.5,
            frameon=False,
            loc="best",
        )


    axes[
        0
    ].set_ylabel(
        "Fraction of events / bin"
    )


    fig.suptitle(
        rf"{xlabel} distributions across the analysis "
        rf"$x_B$ and $-t'$ bins",
        fontsize=15,
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.92,
        )
    )


    save(
        fig,
        path,
    )


# =============================================================================
# NEW: pion-momentum distributions
# =============================================================================

def plot_pion_momentum_by_analysis_bin(
    data,
    path,
):

    """Six -t' panels, each containing the four xB pion-momentum distributions."""

    fig, axes = plt.subplots(
        1,
        6,
        figsize=(
            20.0,
            4.3,
        ),
        sharex=True,
        sharey=True,
    )


    # Slightly finer plotting bins than the exported contamination histogram.
    plot_edges = np.linspace(
        PION_P_MIN,
        PION_P_MAX,
        91,
    )


    centers = 0.5 * (
        plot_edges[:-1]
        + plot_edges[1:]
    )


    for tp_index, (
        ax,
        (
            tp_low,
            tp_high,
        ),
    ) in enumerate(
        zip(
            axes,
            TP_BINS,
        )
    ):

        for xb_index, (
            xb_low,
            xb_high,
        ) in enumerate(
            XB_BINS
        ):

            analysis_bin = (
                xb_index
                * len(
                    TP_BINS
                )
                + tp_index
                + 1
            )


            mask = (
                (
                    data[
                        "analysis_bin"
                    ]
                    == analysis_bin
                )

                & np.isfinite(
                    data[
                        "p_p"
                    ]
                )

                & (
                    data[
                        "p_p"
                    ]
                    >= PION_P_MIN
                )

                & (
                    data[
                        "p_p"
                    ]
                    < PION_P_MAX
                )
            )


            momentum = data[
                "p_p"
            ][
                mask
            ]


            counts, _ = np.histogram(
                momentum,
                bins=plot_edges,
            )


            total = int(
                counts.sum()
            )


            if total > 0:

                fraction = (
                    counts
                    / total
                )


            else:

                fraction = counts.astype(
                    float
                )


            ax.step(
                centers,
                fraction,
                where="mid",
                linewidth=1.35,
                label=(
                    fr"${xb_low:.2f}\leq x_B<{xb_high:.2f}$"
                    fr" ($N={total:,}$)"
                ),
            )


        ax.set_title(
            fr"${tp_low:.2f}\leq -t'<{tp_high:.2f}$"
        )


        ax.set_xlabel(
            r"$p_{\pi^+}$ (GeV)"
        )


        ax.set_xlim(
            PION_P_MIN,
            PION_P_MAX,
        )


        ax.tick_params(
            direction="in",
            top=True,
            right=True,
        )


        ax.legend(
            fontsize=6.5,
            frameon=False,
            loc="best",
        )


    axes[
        0
    ].set_ylabel(
        "Fraction of events / momentum bin"
    )


    fig.suptitle(
        r"$\pi^+$ momentum distributions in the 24 analysis bins",
        fontsize=15,
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.92,
        )
    )


    save(
        fig,
        path,
    )


# =============================================================================
# NEW: portable period-resolved pion-momentum histogram
# =============================================================================

def write_pion_momentum_histograms(
    data,
    path,
):

    """Write momentum populations for direct combination with PID study.

    Each row corresponds to

        period x analysis bin x pion-momentum interval.

    The momentum intervals exactly match beta_pid_2d_fit's portable
    contamination lookup: 0.50--5.00 GeV in 0.25-GeV bins.
    """

    rows = []


    for period in PERIODS:

        period_mask = (
            data[
                "period"
            ]
            == period
        )


        for xb_index, (
            xb_low,
            xb_high,
        ) in enumerate(
            XB_BINS
        ):

            for tp_index, (
                tp_low,
                tp_high,
            ) in enumerate(
                TP_BINS
            ):

                analysis_bin = (
                    xb_index
                    * len(
                        TP_BINS
                    )
                    + tp_index
                    + 1
                )


                analysis_mask = (
                    period_mask

                    & (
                        data[
                            "analysis_bin"
                        ]
                        == analysis_bin
                    )

                    & np.isfinite(
                        data[
                            "p_p"
                        ]
                    )

                    & (
                        data[
                            "p_p"
                        ]
                        >= PION_P_MIN
                    )

                    & (
                        data[
                            "p_p"
                        ]
                        < PION_P_MAX
                    )
                )


                momentum = data[
                    "p_p"
                ][
                    analysis_mask
                ]


                counts, _ = np.histogram(
                    momentum,
                    bins=PION_P_EDGES,
                )


                total = int(
                    counts.sum()
                )


                for momentum_index, count in enumerate(
                    counts
                ):

                    p_low = float(
                        PION_P_EDGES[
                            momentum_index
                        ]
                    )


                    p_high = float(
                        PION_P_EDGES[
                            momentum_index + 1
                        ]
                    )


                    fraction = (
                        float(
                            count
                        )
                        / total

                        if total > 0

                        else 0.0
                    )


                    rows.append(
                        {
                            "period":
                                period,

                            "analysis_bin":
                                analysis_bin,

                            "xB_bin":
                                xb_index + 1,

                            "tprime_bin":
                                tp_index + 1,

                            "xB_min":
                                xb_low,

                            "xB_max":
                                xb_high,

                            "minus_tprime_min_GeV2":
                                tp_low,

                            "minus_tprime_max_GeV2":
                                tp_high,

                            "p_min_GeV":
                                p_low,

                            "p_max_GeV":
                                p_high,

                            "p_center_GeV":
                                0.5
                                * (
                                    p_low
                                    + p_high
                                ),

                            "event_count":
                                int(
                                    count
                                ),

                            "analysis_bin_period_events_0p5_to_5":
                                total,

                            "fraction_of_analysis_bin_period":
                                fraction,
                        }
                    )


    fieldnames = [
        "period",
        "analysis_bin",
        "xB_bin",
        "tprime_bin",
        "xB_min",
        "xB_max",
        "minus_tprime_min_GeV2",
        "minus_tprime_max_GeV2",
        "p_min_GeV",
        "p_max_GeV",
        "p_center_GeV",
        "event_count",
        "analysis_bin_period_events_0p5_to_5",
        "fraction_of_analysis_bin_period",
    ]


    with path.open(
        "w",
        newline="",
    ) as output_file:

        writer = csv.DictWriter(
            output_file,
            fieldnames=fieldnames,
        )


        writer.writeheader()


        writer.writerows(
            rows
        )


    print(
        f"[csv] wrote {path}",
        flush=True,
    )


    # -------------------------------------------------------------------------
    # Useful terminal summary.
    # -------------------------------------------------------------------------

    print(
        "\n[pion momentum] period-resolved populations for contamination folding:",
        flush=True,
    )


    for analysis_bin in range(
        1,
        25,
    ):

        mask = (
            data[
                "analysis_bin"
            ]
            == analysis_bin
        )


        momentum = data[
            "p_p"
        ][
            mask
        ]


        momentum = momentum[
            np.isfinite(
                momentum
            )
        ]


        n_total = len(
            momentum
        )


        n_pid_range = int(
            np.count_nonzero(
                (
                    momentum
                    >= PION_P_MIN
                )

                & (
                    momentum
                    < PION_P_MAX
                )
            )
        )


        n_above_15 = int(
            np.count_nonzero(
                (
                    momentum
                    >= 1.50
                )

                & (
                    momentum
                    < PION_P_MAX
                )
            )
        )


        fraction_above_15 = (
            n_above_15
            / n_pid_range

            if n_pid_range

            else np.nan
        )


        print(
            f"  bin {analysis_bin:2d}: "
            f"N(all p)={n_total:7d}, "
            f"N(0.5-5)={n_pid_range:7d}, "
            f"N(1.5-5)={n_above_15:7d}, "
            f"fraction p>=1.5="
            f"{100.0 * fraction_above_15:6.2f}%",

            flush=True,
        )


# =============================================================================
# Main
# =============================================================================

def main():

    parser = argparse.ArgumentParser()


    parser.add_argument(
        "--cuts",
        type=Path,
        default=DEFAULT_CUTS,
    )


    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUT,
    )


    for period in PERIODS:

        parser.add_argument(
            f"--{period}",
            type=Path,
            default=DEFAULT_INPUTS[
                period
            ],
        )


    args = parser.parse_args()


    args.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )


    windows = load_nominal_windows(
        args.cuts
    )


    datasets = [
        collect(
            period,
            getattr(
                args,
                period,
            ),
            windows[
                period
            ],
        )

        for period in PERIODS
    ]


    data = {
        key:
            np.concatenate(
                [
                    dataset[
                        key
                    ]
                    for dataset in datasets
                ]
            )

        for key in datasets[
            0
        ]
    }


    n = data[
        "xB"
    ].size


    print(
        f"[summary] combined selected sample: {n:,} events",
        flush=True,
    )


    q2max = max(
        8.0,
        np.nanpercentile(
            data[
                "Q2"
            ],
            99.7,
        ),
    )


    wmax = (
        max(
            4.5,
            np.nanpercentile(
                data[
                    "W"
                ],
                99.7,
            ),
        )

        if np.isfinite(
            data[
                "W"
            ]
        ).any()

        else 4.5
    )


    # =========================================================================
    # 01: one-dimensional kinematics
    # =========================================================================

    one_d = [
        (
            "xB",
            np.linspace(
                0.10,
                0.60,
                51,
            ),
            r"$x_B$",
        ),

        (
            "Q2",
            np.linspace(
                1,
                q2max,
                55,
            ),
            r"$Q^2$ (GeV$^2$)",
        ),

        (
            "minus_tprime",
            np.linspace(
                0.05,
                1.25,
                49,
            ),
            r"$-t'$ (GeV$^2$)",
        ),

        (
            "phi",
            np.linspace(
                0,
                360,
                49,
            ),
            r"$\phi$ (deg)",
        ),
    ]


    if np.isfinite(
        data[
            "W"
        ]
    ).any():

        one_d.append(
            (
                "W",
                np.linspace(
                    2,
                    wmax,
                    50,
                ),
                r"$W$ (GeV)",
            )
        )


    if np.isfinite(
        data[
            "y"
        ]
    ).any():

        one_d.append(
            (
                "y",
                np.linspace(
                    0,
                    0.8,
                    49,
                ),
                r"$y$",
            )
        )


    fig, axes = plt.subplots(
        2,
        3,
        figsize=(
            11.5,
            7.0,
        ),
    )


    for ax, (
        key,
        bins,
        label,
    ) in zip(
        axes.flat,
        one_d,
    ):

        ax.hist(
            finite(
                data[
                    key
                ]
            ),
            bins=bins,
            histtype="step",
            linewidth=1.5,
        )


        ax.set_xlabel(
            label
        )


        ax.set_ylabel(
            "Events"
        )


        ax.tick_params(
            direction="in",
            top=True,
            right=True,
        )


    for ax in axes.flat[
        len(
            one_d
        ):
    ]:

        ax.axis(
            "off"
        )


    fig.suptitle(
        r"RGC $e n \pi^+$ kinematics after nominal exclusivity selection"
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.96,
        )
    )


    save(
        fig,
        args.output_dir
        / "01_kinematic_1d.png",
    )


    # =========================================================================
    # 02: core correlations
    # =========================================================================

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(
            14,
            4.1,
        ),
    )


    specs = [
        (
            data[
                "xB"
            ],
            data[
                "Q2"
            ],
            (
                50,
                55,
            ),
            (
                (
                    0.10,
                    0.60,
                ),
                (
                    1,
                    q2max,
                ),
            ),
            r"$x_B$",
            r"$Q^2$ (GeV$^2$)",
        ),

        (
            data[
                "xB"
            ],
            data[
                "minus_tprime"
            ],
            (
                50,
                48,
            ),
            (
                (
                    0.10,
                    0.60,
                ),
                (
                    0.05,
                    1.25,
                ),
            ),
            r"$x_B$",
            r"$-t'$ (GeV$^2$)",
        ),

        (
            data[
                "Q2"
            ],
            data[
                "minus_tprime"
            ],
            (
                55,
                48,
            ),
            (
                (
                    1,
                    q2max,
                ),
                (
                    0.05,
                    1.25,
                ),
            ),
            r"$Q^2$ (GeV$^2$)",
            r"$-t'$ (GeV$^2$)",
        ),
    ]


    for ax, spec in zip(
        axes,
        specs,
    ):

        mesh = hist2d(
            ax,
            *spec,
        )


        fig.colorbar(
            mesh,
            ax=ax,
            label="Events",
        )


    fig.suptitle(
        "Non-redundant kinematic correlations"
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.94,
        )
    )


    save(
        fig,
        args.output_dir
        / "02_core_correlations.png",
    )


    # =========================================================================
    # 03: phi correlations
    # =========================================================================

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(
            9,
            10,
        ),
        sharex=True,
    )


    phi_edges = np.linspace(
        0,
        360,
        49,
    )


    y_specs = [
        (
            data[
                "Q2"
            ],
            np.linspace(
                1,
                q2max,
                55,
            ),
            r"$Q^2$ (GeV$^2$)",
        ),

        (
            data[
                "xB"
            ],
            np.linspace(
                0.10,
                0.60,
                51,
            ),
            r"$x_B$",
        ),

        (
            data[
                "minus_tprime"
            ],
            np.linspace(
                0.05,
                1.25,
                49,
            ),
            r"$-t'$ (GeV$^2$)",
        ),
    ]


    for ax, (
        yy,
        y_edges,
        ylabel,
    ) in zip(
        axes,
        y_specs,
    ):

        good = (
            np.isfinite(
                data[
                    "phi"
                ]
            )

            & np.isfinite(
                yy
            )
        )


        histogram, x_edges, y_edges_hist = np.histogram2d(
            data[
                "phi"
            ][
                good
            ],
            yy[
                good
            ],
            bins=(
                phi_edges,
                y_edges,
            ),
        )


        positive = histogram[
            histogram > 0
        ]


        mesh = ax.pcolormesh(
            x_edges,
            y_edges_hist,
            histogram.T,
            shading="auto",
            norm=LogNorm(
                vmin=(
                    max(
                        1.0,
                        float(
                            positive.min()
                        ),
                    )

                    if positive.size

                    else 1.0
                ),

                vmax=max(
                    2.0,
                    float(
                        histogram.max()
                    ),
                ),
            ),
            cmap=CMAP,
        )


        ax.set_ylabel(
            ylabel
        )


        fig.colorbar(
            mesh,
            ax=ax,
            label="Events",
        )


    axes[
        -1
    ].set_xlabel(
        r"$\phi$ (deg)"
    )


    axes[
        -1
    ].set_xlim(
        0,
        360,
    )


    fig.suptitle(
        r"Common-event-sample correlations with $\phi$"
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.96,
        )
    )


    save(
        fig,
        args.output_dir
        / "03_phi_correlations.png",
    )


    # =========================================================================
    # 04: phi projections
    # =========================================================================

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(
            10.5,
            7.5,
        ),
        sharex=True,
        sharey=True,
    )


    for ax, (
        xb_low,
        xb_high,
    ) in zip(
        axes.flat,
        XB_BINS,
    ):

        mask = (
            (
                data[
                    "xB"
                ]
                >= xb_low
            )

            & (
                data[
                    "xB"
                ]
                < xb_high
            )
        )


        counts, edges = np.histogram(
            data[
                "phi"
            ][
                mask
            ],
            bins=phi_edges,
        )


        centers = 0.5 * (
            edges[:-1]
            + edges[1:]
        )


        density = (
            counts
            / counts.sum()

            if counts.sum()

            else counts.astype(
                float
            )
        )


        ax.step(
            centers,
            density,
            where="mid",
            linewidth=1.5,
        )


        ax.set_title(
            fr"${xb_low:.2f} \leq x_B < {xb_high:.2f}$  "
            fr"($N={counts.sum():,}$)"
        )


        ax.set_xlim(
            0,
            360,
        )


        ax.tick_params(
            direction="in",
            top=True,
            right=True,
        )


    for ax in axes[
        -1,
        :
    ]:

        ax.set_xlabel(
            r"$\phi$ (deg)"
        )


    for ax in axes[
        :,
        0
    ]:

        ax.set_ylabel(
            "Fraction of events / bin"
        )


    fig.suptitle(
        r"$\phi$ acceptance shape in the four analysis $x_B$ intervals"
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.95,
        )
    )


    save(
        fig,
        args.output_dir
        / "04_phi_projections.png",
    )


    # =========================================================================
    # 05: Q2 distributions
    # =========================================================================

    plot_variable_by_tprime_with_xb_overlays(
        data,
        "Q2",
        np.linspace(
            1,
            q2max,
            45,
        ),
        r"$Q^2$ (GeV$^2$)",
        args.output_dir
        / "05_Q2_by_analysis_bin.png",
    )


    # =========================================================================
    # 06: W distributions
    # =========================================================================

    if np.isfinite(
        data[
            "W"
        ]
    ).any():

        plot_variable_by_tprime_with_xb_overlays(
            data,
            "W",
            np.linspace(
                2,
                wmax,
                45,
            ),
            r"$W$ (GeV)",
            args.output_dir
            / "06_W_by_analysis_bin.png",
        )


    # =========================================================================
    # 07: Q2-W by xB
    # =========================================================================

    if np.isfinite(
        data[
            "W"
        ]
    ).any():

        fig, axes = plt.subplots(
            2,
            2,
            figsize=(
                10.5,
                8.0,
            ),
            sharex=True,
            sharey=True,
        )


        for ax, (
            xb_low,
            xb_high,
        ) in zip(
            axes.flat,
            XB_BINS,
        ):

            mask = (
                (
                    data[
                        "xB"
                    ]
                    >= xb_low
                )

                & (
                    data[
                        "xB"
                    ]
                    < xb_high
                )
            )


            mesh = hist2d(
                ax,
                data[
                    "Q2"
                ][
                    mask
                ],
                data[
                    "W"
                ][
                    mask
                ],
                (
                    50,
                    45,
                ),
                (
                    (
                        1,
                        q2max,
                    ),
                    (
                        2,
                        wmax,
                    ),
                ),
                r"$Q^2$ (GeV$^2$)",
                r"$W$ (GeV)",
            )


            ax.set_title(
                fr"${xb_low:.2f}\leq x_B<{xb_high:.2f}$  "
                fr"($N={np.count_nonzero(mask):,}$)"
            )


            fig.colorbar(
                mesh,
                ax=ax,
                label="Events",
            )


        fig.suptitle(
            r"$Q^2$--$W$ phase space in the four analysis $x_B$ intervals"
        )


        fig.tight_layout(
            rect=(
                0,
                0,
                1,
                0.95,
            )
        )


        save(
            fig,
            args.output_dir
            / "07_Q2_W_by_xB.png",
        )


    # =========================================================================
    # 08: NEW pion-momentum distributions
    #
    # Four panels = four xB bins.
    # Six curves/panel = six -t' bins.
    # =========================================================================

    plot_pion_momentum_by_analysis_bin(
        data,
        args.output_dir
        / "08_pion_momentum_by_analysis_bin.png",
    )


    # =========================================================================
    # NEW portable momentum histogram for kaon-contamination folding.
    # =========================================================================

    write_pion_momentum_histograms(
        data,
        args.output_dir
        / "pion_momentum_histograms_by_analysis_bin.csv",
    )


    # =========================================================================
    # Summary JSON
    # =========================================================================

    summary = {
        "cuts_json":
            str(
                args.cuts.resolve()
            ),

        "nominal_sigma_multiple":
            2.0,

        "total_selected_events":
            int(
                n
            ),

        "selected_by_period": {
            period:
                int(
                    np.count_nonzero(
                        data[
                            "period"
                        ]
                        == period
                    )
                )

            for period in PERIODS
        },

        "inputs": {
            period:
                str(
                    getattr(
                        args,
                        period,
                    ).resolve()
                )

            for period in PERIODS
        },

        "explicit_dis_cuts": {
            "Q2_min_gev2":
                1.0,

            "W_min_gev_if_branch_present":
                2.0,

            "y_max_if_branch_present":
                0.80,
        },

        "analysis_binning": {
            "xB":
                XB_BINS,

            "minus_tprime_gev2":
                TP_BINS,
        },

        "pion_momentum": {
            "branch":
                "p_p",

            "histogram_min_GeV":
                PION_P_MIN,

            "histogram_max_GeV":
                PION_P_MAX,

            "histogram_width_GeV":
                PION_P_WIDTH,

            "portable_histogram_csv":
                "pion_momentum_histograms_by_analysis_bin.csv",

            "purpose":
                (
                    "Period-resolved momentum populations for direct "
                    "combination with beta_pid_2d_fit kaon-contamination "
                    "lookup/replica outputs."
                ),
        },

        "plot_colormap_2d":
            CMAP,
    }


    summary_path = (
        args.output_dir
        / "kinematic_summary.json"
    )


    summary_path.write_text(
        json.dumps(
            summary,
            indent=2,
        )
    )


    print(
        f"[summary] wrote {summary_path}",
        flush=True,
    )


    print(
        "\n"
        "============================================================\n"
        " Kinematic-distribution study complete\n"
        "============================================================\n"
        "New pion-momentum products:\n"
        "  08_pion_momentum_by_analysis_bin.png\n"
        "  pion_momentum_histograms_by_analysis_bin.csv\n\n"
        "The CSV uses the same 0.25-GeV bins from 0.5 to 5.0 GeV as the\n"
        "beta_pid_2d_fit contamination lookup and retains Su22/Fa22/Sp23\n"
        "separately for direct period-dependent contamination folding.",

        flush=True,
    )


if __name__ == "__main__":

    main()