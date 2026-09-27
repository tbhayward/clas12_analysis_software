#!/usr/bin/env python3

"""
RGC positive-hadron beta-band study in momentum slices.

Purpose
-------
Use the existing RGC calibration ROOT trees to examine the measured beta
distribution for all positive FD hadron tracks in bins of momentum.

The goal is to directly resolve the pi+, K+, and proton beta bands and quantify
their separation/overlap.

Important
---------
This is a PID diagnostic only.

The fitted +/-3.5 sigma pion interval below is an empirical beta-band overlap
diagnostic. It is NOT a replacement for, or reinterpretation of, the
production REC::Particle chi2pid cut.

Inputs
------
/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/
    rgc_su22_inb_NH3_epi+X_calibration.root
    rgc_fa22_inb_NH3_epi+X_calibration.root
    rgc_sp23_inb_NH3_epi+X_calibration.root

Outputs
-------
output/beta_pid_momentum_fits/

For each period:
    beta_momentum_slices_<period>.png
    beta_momentum_slices_<period>.pdf
    beta_fit_summary_<period>.csv
    beta_overlap_summary_<period>.csv
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

OUTDIR = Path("output/beta_pid_momentum_fits")


# Positive hadron PIDs stored in the existing calibration trees.
POSITIVE_HADRON_PIDS = (
    211,    # pi+
    321,    # K+
    2212,   # proton
)


# Particle masses in GeV.
MASS = {
    "pi": 0.13957039,
    "K": 0.493677,
    "p": 0.9382720813,
}


# Momentum bins.
#
# These intentionally resemble the layout of the old beta-vs-p diagnostic.
P_RANGES = [
    (0.50, 1.00),
    (1.00, 1.50),
    (1.50, 2.00),
    (2.00, 2.50),
    (2.50, 3.00),
    (3.00, 3.50),
    (3.50, 4.00),
    (4.00, 5.00),
]


BETA_RANGE = (0.20, 1.05)

N_BETA_BINS = 340

MIN_FIT_ENTRIES = 500


# Fit bounds.
#
# The expected beta values provide the initial positions for the pi/K/p peaks.
# The fitted means can move somewhat to accommodate the actual detector
# response.
MEAN_WINDOW = 0.035

SIGMA_MIN = 0.002
SIGMA_MAX = 0.080


# =============================================================================
# ROOT helpers
# =============================================================================

def find_tree(root_file):
    """
    Return the first TTree found in a ROOT file.
    """

    for _, obj in root_file.items(recursive=True):

        if isinstance(obj, uproot.behaviors.TTree.TTree):
            return obj

    raise RuntimeError(
        "No TTree was found in the ROOT file."
    )


def find_branch(tree, *aliases):
    """
    Find the first available branch from a list of possible names.
    """

    names = {
        str(key).split(";")[0]
        for key in tree.keys()
    }

    for alias in aliases:

        if alias in names:
            return alias

    raise KeyError(
        "Could not find any of the requested branches:\n"
        f"    {aliases}\n"
        "Available branches include:\n"
        f"    {sorted(names)}"
    )


# =============================================================================
# Physics helpers
# =============================================================================

def beta_expected(momentum, mass):
    """
    Relativistic beta for a particle with momentum p and mass m.

        beta = p / sqrt(p^2 + m^2)
    """

    momentum = np.asarray(
        momentum,
        dtype=float,
    )

    return momentum / np.sqrt(
        momentum * momentum + mass * mass
    )


# =============================================================================
# Fit model
# =============================================================================

def gaussian_counts(
    x,
    area,
    mu,
    sigma,
    bin_width,
):
    """
    Gaussian expressed in histogram counts.

    'area' corresponds approximately to the integrated number of particles
    belonging to that Gaussian component.
    """

    normalization = (
        area
        * bin_width
        / (
            np.sqrt(2.0 * np.pi)
            * sigma
        )
    )

    exponent = (
        -0.5
        * (
            (x - mu)
            / sigma
        ) ** 2
    )

    return normalization * np.exp(
        exponent
    )


def mixture_model(
    parameters,
    x,
    bin_width,
):
    """
    Three-Gaussian pi/K/p model.

    Parameter ordering:

        [Npi, mupi, sigmapi,
         NK,  muK,  sigmaK,
         Np,  mup,  sigmap]
    """

    result = np.zeros_like(
        x,
        dtype=float,
    )

    for index in range(3):

        area = parameters[
            3 * index
        ]

        mu = parameters[
            3 * index + 1
        ]

        sigma = parameters[
            3 * index + 2
        ]

        result += gaussian_counts(
            x,
            area,
            mu,
            sigma,
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
    """
    Integral of one fitted Gaussian between lower and upper beta.
    """

    upper_cdf = ndtr(
        (upper - mu)
        / sigma
    )

    lower_cdf = ndtr(
        (lower - mu)
        / sigma
    )

    return area * (
        upper_cdf
        - lower_cdf
    )


# =============================================================================
# Beta-distribution fitting
# =============================================================================

def fit_beta_distribution(
    beta_values,
    p_min,
    p_max,
):
    """
    Fit one momentum-bin beta distribution with pi/K/p Gaussians.

    At high momentum the pi/K peaks may become poorly separated. Those cases
    are retained in the output but explicitly flagged using the fitted
    separation metric.
    """

    counts, edges = np.histogram(
        beta_values,
        bins=N_BETA_BINS,
        range=BETA_RANGE,
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


    # -------------------------------------------------------------------------
    # Expected beta positions
    # -------------------------------------------------------------------------

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
                MASS[name],
            )
            for name in species
        ]
    )


    # -------------------------------------------------------------------------
    # Initial population estimates
    # -------------------------------------------------------------------------

    initial_areas = []

    for mu in expected_mu:

        near_peak = (
            np.abs(
                centers - mu
            )
            < 0.025
        )

        estimate = float(
            counts[
                near_peak
            ].sum()
        )

        initial_areas.append(
            max(
                estimate,
                0.03 * n_entries,
            )
        )


    initial_sigma = np.array(
        [
            0.012,
            0.014,
            0.016,
        ]
    )


    initial_parameters = np.column_stack(
        (
            np.asarray(
                initial_areas
            ),
            expected_mu,
            initial_sigma,
        )
    ).ravel()


    # -------------------------------------------------------------------------
    # Fit bounds
    # -------------------------------------------------------------------------

    lower_bounds = []
    upper_bounds = []

    for mu in expected_mu:

        lower_bounds.extend(
            [
                0.0,

                max(
                    BETA_RANGE[0],
                    mu - MEAN_WINDOW,
                ),

                SIGMA_MIN,
            ]
        )

        upper_bounds.extend(
            [
                10.0 * n_entries,

                min(
                    BETA_RANGE[1],
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


    # -------------------------------------------------------------------------
    # Fit
    # -------------------------------------------------------------------------

    errors = np.sqrt(
        np.maximum(
            counts,
            1.0,
        )
    )


    def residual(parameters):

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
            max_nfev=4000,
        )

    except Exception as error:

        warnings.warn(
            "Fit failed for "
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
        len(counts)
        - len(parameters),
        1,
    )


    # -------------------------------------------------------------------------
    # Store fitted components
    # -------------------------------------------------------------------------

    components = {}

    for index, name in enumerate(
        species
    ):

        components[name] = {
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


    # -------------------------------------------------------------------------
    # Separation metrics
    # -------------------------------------------------------------------------

    separation_pi_k = (
        abs(
            pi_component["mu"]
            - k_component["mu"]
        )
        / np.sqrt(
            pi_component["sigma"] ** 2
            + k_component["sigma"] ** 2
        )
    )


    separation_k_p = (
        abs(
            k_component["mu"]
            - p_component["mu"]
        )
        / np.sqrt(
            k_component["sigma"] ** 2
            + p_component["sigma"] ** 2
        )
    )


    # -------------------------------------------------------------------------
    # Empirical fitted overlap with a +/-3.5 sigma pion band
    #
    # Again: this is NOT being claimed to reproduce the production chi2pid
    # definition. It is simply a useful way of quantifying the fitted beta
    # overlap.
    # -------------------------------------------------------------------------

    pion_lower = (
        pi_component["mu"]
        - 3.5
        * pi_component["sigma"]
    )

    pion_upper = (
        pi_component["mu"]
        + 3.5
        * pi_component["sigma"]
    )


    inside = {}

    for name in species:

        component = components[
            name
        ]

        inside[name] = gaussian_integral(
            component["area"],
            component["mu"],
            component["sigma"],
            pion_lower,
            pion_upper,
        )


    total_inside = sum(
        inside.values()
    )


    if total_inside > 0:

        nonpion_fraction = (
            inside["K"]
            + inside["p"]
        ) / total_inside

    else:

        nonpion_fraction = np.nan


    # A very simple warning criterion.
    #
    # We should inspect the actual distributions before deciding what numerical
    # separation threshold is appropriate for a publication-level statement.
    reliable = bool(
        fit.success
        and separation_pi_k >= 1.0
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

            "nonpion_fraction":
                nonpion_fraction,

            "reliable":
                reliable,
        },
    )


# =============================================================================
# Input
# =============================================================================

def load_positive_fd_tracks(
    filename,
):
    """
    Load all positive FD hadron tracks represented in the calibration tree.

    No pion PID cut and no chi2pid cut are applied here.

    The existing calibration trees contain REC-assigned pi+, K+, and proton
    tracks, which are the positive hadron assignments retained by the original
    calibration skim.
    """

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


        momentum_branch = find_branch(
            tree,
            "p",
            "particle_p",
        )


        beta_branch = find_branch(
            tree,
            "particle_beta",
            "beta",
        )


        status_branch = find_branch(
            tree,
            "particle_status",
            "status",
        )


        arrays = tree.arrays(
            [
                pid_branch,
                momentum_branch,
                beta_branch,
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
            momentum_branch
        ],
        dtype=float,
    )


    beta = np.asarray(
        arrays[
            beta_branch
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


    selection = (
        np.isin(
            pid,
            POSITIVE_HADRON_PIDS,
        )

        & (
            status
            >= 2000
        )

        & (
            status
            < 4000
        )

        & np.isfinite(
            momentum
        )

        & np.isfinite(
            beta
        )

        & (
            momentum
            >= 0.5
        )

        & (
            momentum
            < 5.0
        )

        & (
            beta
            > BETA_RANGE[0]
        )

        & (
            beta
            < BETA_RANGE[1]
        )
    )


    return (
        momentum[
            selection
        ],

        beta[
            selection
        ],

        pid[
            selection
        ],
    )


# =============================================================================
# Plot one period
# =============================================================================

def analyze_period(
    period,
    filename,
):
    """
    Produce the complete 3x3 beta-distribution figure and CSV summaries.
    """

    print(
        f"[{period}] reading {filename}",
        flush=True,
    )


    (
        momentum,
        beta,
        assigned_pid,
    ) = load_positive_fd_tracks(
        filename
    )


    print(
        f"[{period}] selected positive FD tracks: "
        f"{len(momentum):,}",
        flush=True,
    )


    fig, axes = plt.subplots(
        3,
        3,
        figsize=(
            14.5,
            11.0,
        ),
        sharex=False,
        sharey=False,
    )


    axes = axes.ravel()


    fit_rows = []
    overlap_rows = []


    # =========================================================================
    # Integrated panel
    #
    # Do NOT fit the integrated distribution with three Gaussians because the
    # expected beta peak positions vary significantly over 0.5 < p < 5 GeV.
    # =========================================================================

    counts, edges = np.histogram(
        beta,
        bins=N_BETA_BINS,
        range=BETA_RANGE,
    )


    centers = 0.5 * (
        edges[:-1]
        + edges[1:]
    )


    axes[0].step(
        centers,
        counts,
        where="mid",
        linewidth=1.0,
    )


    axes[0].set_title(
        "Integrated\n"
        f"N = {len(beta):,}"
    )


    axes[0].set_xlabel(
        r"$\beta$"
    )


    axes[0].set_ylabel(
        "Counts"
    )


    axes[0].set_xlim(
        *BETA_RANGE
    )


    # =========================================================================
    # Momentum slices
    # =========================================================================

    for panel, (
        p_min,
        p_max,
    ) in enumerate(
        P_RANGES,
        start=1,
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
        )


        beta_bin = beta[
            selection
        ]


        (
            counts,
            edges,
            result,
        ) = fit_beta_distribution(
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
        # Expected beta positions
        # ---------------------------------------------------------------------

        p_center = 0.5 * (
            p_min
            + p_max
        )


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


            total_model = mixture_model(
                parameters,
                centers,
                bin_width,
            )


            ax.plot(
                centers,
                total_model,
                linewidth=1.5,
                label="3-Gaussian fit",
            )


            for index, species in enumerate(
                (
                    "pi",
                    "K",
                    "p",
                )
            ):

                area = parameters[
                    3 * index
                ]

                mu = parameters[
                    3 * index + 1
                ]

                sigma = parameters[
                    3 * index + 2
                ]


                component = gaussian_counts(
                    centers,
                    area,
                    mu,
                    sigma,
                    bin_width,
                )


                ax.plot(
                    centers,
                    component,
                    linestyle="--",
                    linewidth=1.0,
                )


            if result[
                "reliable"
            ]:

                quality_text = (
                    "resolved"
                )

            else:

                quality_text = (
                    "fit caution"
                )


            ax.text(
                0.03,
                0.96,

                (
                    f"N = {len(beta_bin):,}\n"

                    rf"$S_{{\pi K}}"
                    rf"={result['separation_piK']:.2f}$"
                    "\n"

                    rf"$S_{{Kp}}"
                    rf"={result['separation_Kp']:.2f}$"
                    "\n"

                    f"{quality_text}"
                ),

                transform=ax.transAxes,

                va="top",
                ha="left",

                fontsize=8,
            )


            # -----------------------------------------------------------------
            # Fit-summary table
            # -----------------------------------------------------------------

            row = {
                "period":
                    period,

                "p_min":
                    p_min,

                "p_max":
                    p_max,

                "N":
                    len(
                        beta_bin
                    ),

                "chi2_ndf":
                    result[
                        "chi2_ndf"
                    ],

                "separation_piK":
                    result[
                        "separation_piK"
                    ],

                "separation_Kp":
                    result[
                        "separation_Kp"
                    ],

                "fit_reliable":
                    result[
                        "reliable"
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


            # -----------------------------------------------------------------
            # Fitted-overlap table
            # -----------------------------------------------------------------

            overlap_rows.append(
                {
                    "period":
                        period,

                    "p_min":
                        p_min,

                    "p_max":
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
                            "nonpion_fraction"
                        ],

                    "fit_reliable":
                        result[
                            "reliable"
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


        # ---------------------------------------------------------------------
        # Panel formatting
        # ---------------------------------------------------------------------

        ax.set_title(
            f"{p_min:.2f} < p < {p_max:.2f} GeV"
        )


        ax.set_xlabel(
            r"$\beta$"
        )


        ax.set_ylabel(
            "Counts"
        )


        ax.set_xlim(
            *BETA_RANGE
        )


    # =========================================================================
    # Save figure
    # =========================================================================

    fig.suptitle(
        f"{period} FD: all positive hadron tracks",
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


    png_output = (
        OUTDIR
        / f"beta_momentum_slices_{period.lower()}.png"
    )


    pdf_output = (
        OUTDIR
        / f"beta_momentum_slices_{period.lower()}.pdf"
    )


    fig.savefig(
        png_output,
        dpi=200,
        bbox_inches="tight",
    )


    fig.savefig(
        pdf_output,
        bbox_inches="tight",
    )


    plt.close(
        fig
    )


    print(
        f"[{period}] wrote {png_output}",
        flush=True,
    )


    print(
        f"[{period}] wrote {pdf_output}",
        flush=True,
    )


    # =========================================================================
    # Save numerical summaries
    # =========================================================================

    fit_table = pd.DataFrame(
        fit_rows
    )


    fit_table.to_csv(
        OUTDIR
        / f"beta_fit_summary_{period.lower()}.csv",

        index=False,
    )


    overlap_table = pd.DataFrame(
        overlap_rows
    )


    overlap_table.to_csv(
        OUTDIR
        / f"beta_overlap_summary_{period.lower()}.csv",

        index=False,
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


        analyze_period(
            period,
            filename,
        )


if __name__ == "__main__":
    main()