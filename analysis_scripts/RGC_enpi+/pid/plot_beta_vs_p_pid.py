#!/usr/bin/env python3

"""
RGC positive-hadron beta PID study.

For each RGC period this script produces:

1. A two-panel beta-vs-p diagnostic:
   - all positive FD hadron tracks before the chi2pid cut;
   - the same tracks after |chi2pid| < 3.5 and 0.5 < p < 5 GeV.

   The expected beta(p) curves for pi+, K+, and p are overlaid.

2. Momentum-sliced beta distributions:
   - p = 1.50--5.00 GeV;
   - 0.25-GeV momentum bins;
   - beta = 0.8--1.1;
   - no integrated panel;
   - four columns;
   - PNG output only.

Per the analysis prescription for this study, contamination below 1.50 GeV
is assumed negligible and is not quantified.

The beta spectra above 1.50 GeV are fit with pi/K/p components. The resulting
non-pion contribution inside an illustrative fitted pion beta band is used
only as a detector/PID overlap diagnostic. It is NOT a reinterpretation of
REC::Particle.chi2pid.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import uproot

from scipy.optimize import least_squares
from scipy.special import ndtr


# =============================================================================
# Configuration
# =============================================================================

INPUTS = {
    "Su22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/"
        "data/pass2/calibration/"
        "rgc_su22_inb_NH3_epi+X_calibration.root"
    ),

    "Fa22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/"
        "data/pass2/calibration/"
        "rgc_fa22_inb_NH3_epi+X_calibration.root"
    ),

    "Sp23": Path(
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/"
        "data/pass2/calibration/"
        "rgc_sp23_inb_NH3_epi+X_calibration.root"
    ),
}


OUTDIR = Path(
    "output/beta_pid_study"
)


POSITIVE_HADRON_PIDS = (
    211,
    321,
    2212,
)


MASS = {
    "pi": 0.13957039,
    "K": 0.493677,
    "p": 0.9382720813,
}


# =============================================================================
# 2D beta-vs-p configuration
# =============================================================================

PLOT_2D_P_RANGE = (
    0.5,
    5.0,
)

PLOT_2D_BETA_RANGE = (
    0.20,
    1.20,
)

PLOT_2D_P_BINS = 180
PLOT_2D_BETA_BINS = 200


# =============================================================================
# Momentum-slice configuration
#
# 1.50 -> 5.00 GeV in 0.25-GeV intervals gives 14 panels.
# =============================================================================

P_SLICE_EDGES = np.arange(
    1.50,
    5.00 + 0.25,
    0.25,
)


P_RANGES = [
    (
        float(P_SLICE_EDGES[index]),
        float(P_SLICE_EDGES[index + 1]),
    )
    for index in range(
        len(P_SLICE_EDGES) - 1
    )
]


SLICE_BETA_RANGE = (
    0.80,
    1.10,
)

N_SLICE_BETA_BINS = 240

MIN_FIT_ENTRIES = 200

MEAN_WINDOW = 0.025

SIGMA_MIN = 0.0015
SIGMA_MAX = 0.050


# =============================================================================
# ROOT helpers
# =============================================================================

def find_tree(
    root_file,
):

    for _, obj in root_file.items(
        recursive=True
    ):

        if isinstance(
            obj,
            uproot.behaviors.TTree.TTree,
        ):

            return obj

    raise RuntimeError(
        "No TTree found in ROOT file."
    )


def find_branch(
    tree,
    *aliases,
):

    names = {
        str(key).split(";")[0]
        for key in tree.keys()
    }

    for alias in aliases:

        if alias in names:
            return alias

    raise KeyError(
        "Could not find any of these branches:\n"
        f"{aliases}\n\n"
        "Available branches:\n"
        f"{sorted(names)}"
    )


# =============================================================================
# Physics
# =============================================================================

def beta_expected(
    momentum,
    mass,
):

    momentum = np.asarray(
        momentum,
        dtype=float,
    )

    return (
        momentum
        / np.sqrt(
            momentum * momentum
            + mass * mass
        )
    )


# =============================================================================
# Gaussian model
# =============================================================================

def gaussian_counts(
    x,
    area,
    mu,
    sigma,
    bin_width,
):

    normalization = (
        area
        * bin_width
        / (
            np.sqrt(
                2.0 * np.pi
            )
            * sigma
        )
    )

    return (
        normalization
        * np.exp(
            -0.5
            * (
                (
                    x - mu
                )
                / sigma
            ) ** 2
        )
    )


def mixture_model(
    parameters,
    x,
    bin_width,
):

    result = np.zeros_like(
        x,
        dtype=float,
    )

    for index in range(
        3
    ):

        result += gaussian_counts(
            x,

            parameters[
                3 * index
            ],

            parameters[
                3 * index + 1
            ],

            parameters[
                3 * index + 2
            ],

            bin_width,
        )

    return result


def gaussian_integral(
    area,
    mu,
    sigma,
    lower,
    upper,
):

    return (
        area
        * (
            ndtr(
                (
                    upper - mu
                )
                / sigma
            )
            -
            ndtr(
                (
                    lower - mu
                )
                / sigma
            )
        )
    )


# =============================================================================
# Load one period
# =============================================================================

def load_period(
    filename,
):

    with uproot.open(
        filename
    ) as root_file:

        tree = find_tree(
            root_file
        )


        pid_branch = find_branch(
            tree,
            "particle_pid",
            "pid",
        )


        p_branch = find_branch(
            tree,
            "p",
            "particle_p",
        )


        beta_branch = find_branch(
            tree,
            "particle_beta",
            "beta",
        )


        chi2pid_branch = find_branch(
            tree,
            "particle_chi2pid",
            "chi2pid",
        )


        status_branch = find_branch(
            tree,
            "particle_status",
            "status",
        )


        arrays = tree.arrays(
            [
                pid_branch,
                p_branch,
                beta_branch,
                chi2pid_branch,
                status_branch,
            ],

            library="np",
        )


    pid = np.asarray(
        arrays[
            pid_branch
        ]
    )


    momentum = np.asarray(
        arrays[
            p_branch
        ],
        dtype=float,
    )


    beta = np.asarray(
        arrays[
            beta_branch
        ],
        dtype=float,
    )


    chi2pid = np.asarray(
        arrays[
            chi2pid_branch
        ],
        dtype=float,
    )


    status = np.abs(
        np.asarray(
            arrays[
                status_branch
            ]
        )
    )


    positive_fd = (
        np.isin(
            pid,
            POSITIVE_HADRON_PIDS,
        )

        & (
            status >= 2000
        )

        & (
            status < 4000
        )

        & np.isfinite(
            momentum
        )

        & np.isfinite(
            beta
        )

        & np.isfinite(
            chi2pid
        )
    )


    return {
        "pid":
            pid[
                positive_fd
            ],

        "p":
            momentum[
                positive_fd
            ],

        "beta":
            beta[
                positive_fd
            ],

        "chi2pid":
            chi2pid[
                positive_fd
            ],
    }


# =============================================================================
# 2D beta-vs-p diagnostic
# =============================================================================

def plot_beta_vs_p(
    period,
    data,
):

    momentum = data[
        "p"
    ]

    beta = data[
        "beta"
    ]

    chi2pid = data[
        "chi2pid"
    ]


    before = (
        (
            momentum
            >= PLOT_2D_P_RANGE[0]
        )

        & (
            momentum
            < PLOT_2D_P_RANGE[1]
        )

        & (
            beta
            >= PLOT_2D_BETA_RANGE[0]
        )

        & (
            beta
            <= PLOT_2D_BETA_RANGE[1]
        )
    )


    after = (
        before

        & (
            np.abs(
                chi2pid
            )
            < 3.5
        )
    )


    fig, axes = plt.subplots(
        1,
        2,

        figsize=(
            14,
            5.5,
        ),
    )


    selections = [
        (
            before,
            "Before PID cuts",
        ),

        (
            after,
            r"After $|\chi^2_{\rm PID}|<3.5$, "
            r"$0.5<p<5.0$ GeV",
        ),
    ]


    p_curve = np.linspace(
        PLOT_2D_P_RANGE[0],
        PLOT_2D_P_RANGE[1],
        600,
    )


    for ax, (
        selection,
        title,
    ) in zip(
        axes,
        selections,
    ):


        hist = ax.hist2d(
            momentum[
                selection
            ],

            beta[
                selection
            ],

            bins=[
                PLOT_2D_P_BINS,
                PLOT_2D_BETA_BINS,
            ],

            range=[
                PLOT_2D_P_RANGE,
                PLOT_2D_BETA_RANGE,
            ],

            norm=LogNorm(
                vmin=1
            ),
        )


        fig.colorbar(
            hist[3],
            ax=ax,
            label="Counts",
        )


        # ---------------------------------------------------------------------
        # Theoretical beta(p) curves.
        #
        # These are especially useful at low momentum, where the proton and
        # kaon bands curve rapidly and therefore generate broad structures when
        # projected over a finite momentum interval.
        # ---------------------------------------------------------------------

        ax.plot(
            p_curve,

            beta_expected(
                p_curve,
                MASS[
                    "pi"
                ],
            ),

            linewidth=1.4,

            label=r"$\pi^+$",
        )


        ax.plot(
            p_curve,

            beta_expected(
                p_curve,
                MASS[
                    "K"
                ],
            ),

            linewidth=1.4,

            label=r"$K^+$",
        )


        ax.plot(
            p_curve,

            beta_expected(
                p_curve,
                MASS[
                    "p"
                ],
            ),

            linewidth=1.4,

            label=r"$p$",
        )


        ax.set_xlabel(
            r"$p$ (GeV)"
        )


        ax.set_ylabel(
            r"$\beta$"
        )


        ax.set_title(
            title
        )


        ax.legend(
            loc="lower right",
            fontsize=9,
        )


    fig.suptitle(
        f"{period} FD positive hadrons",
        fontsize=15,
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.94,
        )
    )


    output = (
        OUTDIR
        / f"beta_vs_p_{period.lower()}.png"
    )


    fig.savefig(
        output,

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


    print(
        f"[{period}] wrote {output}",
        flush=True,
    )


# =============================================================================
# Fit one 0.25-GeV momentum interval
# =============================================================================

def fit_beta_slice(
    beta_values,
    p_min,
    p_max,
):

    counts, edges = np.histogram(
        beta_values,

        bins=N_SLICE_BETA_BINS,

        range=SLICE_BETA_RANGE,
    )


    centers = 0.5 * (
        edges[:-1]
        + edges[1:]
    )


    bin_width = (
        edges[1]
        - edges[0]
    )


    n_entries = int(
        counts.sum()
    )


    if n_entries < MIN_FIT_ENTRIES:

        return (
            counts,
            edges,
            None,
        )


    p_center = 0.5 * (
        p_min
        + p_max
    )


    species = (
        "pi",
        "K",
        "p",
    )


    expected_mu = np.array(
        [
            beta_expected(
                p_center,
                MASS[
                    name
                ],
            )
            for name in species
        ]
    )


    # =========================================================================
    # Initial component populations
    # =========================================================================

    initial_areas = []


    for mu in expected_mu:

        near_peak = (
            np.abs(
                centers - mu
            )
            < 0.012
        )


        estimate = float(
            counts[
                near_peak
            ].sum()
        )


        initial_areas.append(
            max(
                estimate,
                0.02 * n_entries,
            )
        )


    initial_sigmas = np.array(
        [
            0.007,
            0.009,
            0.011,
        ]
    )


    initial_parameters = np.column_stack(
        (
            np.asarray(
                initial_areas
            ),

            expected_mu,

            initial_sigmas,
        )
    ).ravel()


    # =========================================================================
    # Bounds
    # =========================================================================

    lower_bounds = []
    upper_bounds = []


    for mu in expected_mu:

        lower_bounds.extend(
            [
                0.0,

                max(
                    SLICE_BETA_RANGE[0],
                    mu - MEAN_WINDOW,
                ),

                SIGMA_MIN,
            ]
        )


        upper_bounds.extend(
            [
                10.0
                * n_entries,

                min(
                    SLICE_BETA_RANGE[1],
                    mu + MEAN_WINDOW,
                ),

                SIGMA_MAX,
            ]
        )


    lower_bounds = np.asarray(
        lower_bounds
    )


    upper_bounds = np.asarray(
        upper_bounds
    )


    errors = np.sqrt(
        np.maximum(
            counts,
            1.0,
        )
    )


    def residual(
        parameters,
    ):

        prediction = mixture_model(
            parameters,
            centers,
            bin_width,
        )

        return (
            prediction
            - counts
        ) / errors


    try:

        fit = least_squares(
            residual,

            initial_parameters,

            bounds=(
                lower_bounds,
                upper_bounds,
            ),

            x_scale="jac",

            max_nfev=5000,
        )


    except Exception as error:

        warnings.warn(
            f"Fit failed for "
            f"{p_min:.2f} < p < {p_max:.2f} GeV: "
            f"{error}"
        )

        return (
            counts,
            edges,
            None,
        )


    parameters = fit.x


    prediction = mixture_model(
        parameters,
        centers,
        bin_width,
    )


    chi2 = float(
        np.sum(
            (
                (
                    counts
                    - prediction
                )
                / errors
            ) ** 2
        )
    )


    ndf = max(
        len(
            counts
        )
        - len(
            parameters
        ),
        1,
    )


    # =========================================================================
    # Components
    # =========================================================================

    components = {}


    for index, name in enumerate(
        species
    ):

        components[
            name
        ] = {
            "area":
                parameters[
                    3 * index
                ],

            "mu":
                parameters[
                    3 * index + 1
                ],

            "sigma":
                parameters[
                    3 * index + 2
                ],
        }


    pi_component = components[
        "pi"
    ]

    k_component = components[
        "K"
    ]

    p_component = components[
        "p"
    ]


    # =========================================================================
    # Peak separation
    #
    # This is what was previously labeled "S".
    #
    # A value of 3 means that the fitted peak centers differ by three times the
    # quadrature combination of their fitted widths.
    # =========================================================================

    separation_pi_k = (
        abs(
            pi_component[
                "mu"
            ]
            - k_component[
                "mu"
            ]
        )

        / np.sqrt(
            pi_component[
                "sigma"
            ] ** 2

            + k_component[
                "sigma"
            ] ** 2
        )
    )


    separation_k_p = (
        abs(
            k_component[
                "mu"
            ]
            - p_component[
                "mu"
            ]
        )

        / np.sqrt(
            k_component[
                "sigma"
            ] ** 2

            + p_component[
                "sigma"
            ] ** 2
        )
    )


    # =========================================================================
    # Fitted contamination in the pion beta region
    #
    # Use +/-3.5 fitted pion sigma only as an empirical beta-band overlap
    # diagnostic.
    # =========================================================================

    pion_lower = (
        pi_component[
            "mu"
        ]

        - 3.5
        * pi_component[
            "sigma"
        ]
    )


    pion_upper = (
        pi_component[
            "mu"
        ]

        + 3.5
        * pi_component[
            "sigma"
        ]
    )


    inside = {}


    for name in species:

        component = components[
            name
        ]


        inside[
            name
        ] = gaussian_integral(
            component[
                "area"
            ],

            component[
                "mu"
            ],

            component[
                "sigma"
            ],

            pion_lower,
            pion_upper,
        )


    total_inside = sum(
        inside.values()
    )


    if total_inside > 0:

        fitted_nonpion_fraction = (
            inside[
                "K"
            ]

            + inside[
                "p"
            ]
        ) / total_inside


    else:

        fitted_nonpion_fraction = np.nan


    fit_reliable = bool(
        fit.success

        and separation_pi_k
        >= 1.0
    )


    return (
        counts,
        edges,

        {
            "parameters":
                parameters,

            "components":
                components,

            "chi2_ndf":
                chi2 / ndf,

            "separation_piK":
                separation_pi_k,

            "separation_Kp":
                separation_k_p,

            "pion_window":
                (
                    pion_lower,
                    pion_upper,
                ),

            "inside":
                inside,

            "fitted_nonpion_fraction":
                fitted_nonpion_fraction,

            "fit_reliable":
                fit_reliable,
        },
    )


# =============================================================================
# Momentum-sliced beta distributions
# =============================================================================

def plot_momentum_slices(
    period,
    data,
):

    momentum = data[
        "p"
    ]

    beta = data[
        "beta"
    ]


    # 14 momentum bins on a 4x4 canvas.
    n_columns = 4

    n_rows = int(
        np.ceil(
            len(
                P_RANGES
            )
            / n_columns
        )
    )


    fig, axes = plt.subplots(
        n_rows,
        n_columns,

        figsize=(
            16,
            3.7 * n_rows,
        ),
    )


    axes = np.asarray(
        axes
    ).ravel()


    fit_rows = []
    overlap_rows = []


    for panel, (
        p_min,
        p_max,
    ) in enumerate(
        P_RANGES
    ):


        ax = axes[
            panel
        ]


        selection = (
            (
                momentum
                >= p_min
            )

            & (
                momentum
                < p_max
            )

            & (
                beta
                >= SLICE_BETA_RANGE[0]
            )

            & (
                beta
                <= SLICE_BETA_RANGE[1]
            )
        )


        beta_bin = beta[
            selection
        ]


        (
            counts,
            edges,
            result,
        ) = fit_beta_slice(
            beta_bin,
            p_min,
            p_max,
        )


        centers = 0.5 * (
            edges[:-1]
            + edges[1:]
        )


        bin_width = (
            edges[1]
            - edges[0]
        )


        # ---------------------------------------------------------------------
        # Data
        # ---------------------------------------------------------------------

        ax.step(
            centers,
            counts,

            where="mid",

            linewidth=1.0,

            label="Data",
        )


        # ---------------------------------------------------------------------
        # Expected beta positions at center of momentum interval.
        # ---------------------------------------------------------------------

        p_center = 0.5 * (
            p_min
            + p_max
        )


        expected_positions = {}


        for species in (
            "pi",
            "K",
            "p",
        ):

            expected = beta_expected(
                p_center,

                MASS[
                    species
                ],
            )


            expected_positions[
                species
            ] = expected


            ax.axvline(
                expected,

                linestyle=":",

                linewidth=1.0,
            )


        # ---------------------------------------------------------------------
        # Fit
        # ---------------------------------------------------------------------

        if result is not None:

            parameters = result[
                "parameters"
            ]


            total_prediction = mixture_model(
                parameters,
                centers,
                bin_width,
            )


            ax.plot(
                centers,
                total_prediction,

                linewidth=1.5,

                label="Total fit",
            )


            for index, species in enumerate(
                (
                    "pi",
                    "K",
                    "p",
                )
            ):

                component_prediction = gaussian_counts(
                    centers,

                    parameters[
                        3 * index
                    ],

                    parameters[
                        3 * index + 1
                    ],

                    parameters[
                        3 * index + 2
                    ],

                    bin_width,
                )


                ax.plot(
                    centers,
                    component_prediction,

                    linestyle="--",

                    linewidth=1.0,
                )


            nonpion_percent = (
                100.0
                * result[
                    "fitted_nonpion_fraction"
                ]
            )


            if result[
                "fit_reliable"
            ]:

                quality = (
                    "resolved"
                )

            else:

                quality = (
                    "fit caution"
                )


            # No more cryptic "S" labels.
            annotation = (
                f"N = {len(beta_bin):,}\n"

                f"pi/K sep. = "
                f"{result['separation_piK']:.2f}\n"

                f"K/p sep. = "
                f"{result['separation_Kp']:.2f}\n"

                f"non-pion overlap = "
                f"{nonpion_percent:.2f}%\n"

                f"{quality}"
            )


            ax.text(
                0.03,
                0.96,

                annotation,

                transform=ax.transAxes,

                va="top",
                ha="left",

                fontsize=7.5,
            )


            # =================================================================
            # Fit summary
            # =================================================================

            row = {
                "period":
                    period,

                "p_min_GeV":
                    p_min,

                "p_max_GeV":
                    p_max,

                "N":
                    len(
                        beta_bin
                    ),

                "chi2_ndf":
                    result[
                        "chi2_ndf"
                    ],

                "pi_K_separation":
                    result[
                        "separation_piK"
                    ],

                "K_p_separation":
                    result[
                        "separation_Kp"
                    ],

                "fit_reliable":
                    result[
                        "fit_reliable"
                    ],
            }


            for index, species in enumerate(
                (
                    "pi",
                    "K",
                    "p",
                )
            ):

                row[
                    f"{species}_expected_mu"
                ] = expected_positions[
                    species
                ]


                row[
                    f"{species}_area"
                ] = parameters[
                    3 * index
                ]


                row[
                    f"{species}_mu"
                ] = parameters[
                    3 * index + 1
                ]


                row[
                    f"{species}_sigma"
                ] = parameters[
                    3 * index + 2
                ]


            fit_rows.append(
                row
            )


            # =================================================================
            # Overlap summary
            # =================================================================

            overlap_rows.append(
                {
                    "period":
                        period,

                    "p_min_GeV":
                        p_min,

                    "p_max_GeV":
                        p_max,

                    "pion_fit_window_min":
                        result[
                            "pion_window"
                        ][0],

                    "pion_fit_window_max":
                        result[
                            "pion_window"
                        ][1],

                    "fit_pi_in_window":
                        result[
                            "inside"
                        ][
                            "pi"
                        ],

                    "fit_K_in_window":
                        result[
                            "inside"
                        ][
                            "K"
                        ],

                    "fit_p_in_window":
                        result[
                            "inside"
                        ][
                            "p"
                        ],

                    "fit_nonpi_fraction_in_pion_window":
                        result[
                            "fitted_nonpion_fraction"
                        ],

                    "fit_reliable":
                        result[
                            "fit_reliable"
                        ],
                }
            )


        else:

            ax.text(
                0.03,
                0.96,

                (
                    f"N = {len(beta_bin):,}\n"
                    "fit unavailable"
                ),

                transform=ax.transAxes,

                va="top",

                fontsize=8,
            )


        ax.set_title(
            f"{p_min:.2f} < p < {p_max:.2f} GeV",

            fontsize=10,
        )


        ax.set_xlim(
            *SLICE_BETA_RANGE
        )


        ax.set_xlabel(
            r"$\beta$"
        )


        ax.set_ylabel(
            "Counts"
        )


    # =========================================================================
    # Empty final two panels
    # =========================================================================

    for panel in range(
        len(
            P_RANGES
        ),
        len(
            axes
        ),
    ):

        axes[
            panel
        ].axis(
            "off"
        )


    fig.suptitle(
        f"{period} FD positive hadrons: "
        r"$\beta$ distributions in 0.25-GeV momentum bins",

        fontsize=15,
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.965,
        )
    )


    # PNG ONLY.
    output = (
        OUTDIR
        / f"beta_momentum_slices_{period.lower()}.png"
    )


    fig.savefig(
        output,

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


    print(
        f"[{period}] wrote {output}",
        flush=True,
    )


    # =========================================================================
    # Numerical summaries
    # =========================================================================

    fit_output = (
        OUTDIR
        / f"beta_fit_summary_{period.lower()}.csv"
    )


    pd.DataFrame(
        fit_rows
    ).to_csv(
        fit_output,

        index=False,
    )


    overlap_output = (
        OUTDIR
        / f"beta_overlap_summary_{period.lower()}.csv"
    )


    pd.DataFrame(
        overlap_rows
    ).to_csv(
        overlap_output,

        index=False,
    )


    print(
        f"[{period}] wrote {fit_output}",
        flush=True,
    )


    print(
        f"[{period}] wrote {overlap_output}",
        flush=True,
    )


# =============================================================================
# Main
# =============================================================================

def main():

    OUTDIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    for (
        period,
        filename,
    ) in INPUTS.items():


        if not filename.is_file():

            raise FileNotFoundError(
                f"Missing calibration ROOT file:\n"
                f"{filename}"
            )


        print(
            "\n"
            "============================================================\n"
            f" {period}\n"
            "============================================================",

            flush=True,
        )


        data = load_period(
            filename
        )


        print(
            f"[{period}] positive FD hadron tracks: "
            f"{len(data['p']):,}",

            flush=True,
        )


        # 2D beta-vs-p before/after PID.
        plot_beta_vs_p(
            period,
            data,
        )


        # Quantitative study begins at p = 1.50 GeV.
        plot_momentum_slices(
            period,
            data,
        )


if __name__ == "__main__":
    main()