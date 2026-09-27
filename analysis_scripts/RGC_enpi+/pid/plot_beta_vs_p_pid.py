#!/usr/bin/env python3

"""
RGC positive-hadron PID study using a native (p, beta) mixture model.

Purpose
-------
Estimate residual kaon/proton contamination in the actual pion PID sample
as a function of momentum, for later folding through the pion momentum
distributions in the 24 final analysis bins.

The key improvement over the previous projected-beta study is that every
track is evaluated at its own momentum:

    beta_h(p) = p / sqrt(p^2 + m_h^2)

rather than approximating all tracks in a finite momentum interval by a
single Gaussian centered at beta_h(p_center).

Two samples are fit independently:

    1. all_positive
       All positive FD tracks assigned by REC as pi+, K+, or proton.
       This is primarily a validation of the three-species mixture model.

    2. pion_selected
       REC PID == 211 and |chi2pid| < 3.5.
       This is the sample that matters for the contamination estimate.

The stored REC::Particle chi2pid value is used only for the hypothesis to
which it actually corresponds. In particular, the stored chi2pid of an
assigned kaon is NOT interpreted as its pion-hypothesis chi2pid.

Detector response
-----------------
For each species, clean REC-assigned control tracks are used to measure

    Delta beta_h = beta_measured - beta_expected_h(p)

versus momentum.

In narrow momentum bins the robust median gives the offset and the robust
MAD gives the Gaussian-core resolution. Smooth functions in 1/p are then
fit to the measured offset and resolution.

These response functions are fixed in the subsequent mixture fits. This
prevents unresolved high-momentum pi/K populations from being accommodated
by arbitrarily broad Gaussian components.

Mixture fit
-----------
Within each 0.25-GeV momentum interval, the conditional likelihood is

    P(beta_j | p_j) =
          f_pi G_pi(beta_j | p_j)
        + f_K  G_K (beta_j | p_j)
        + f_p  G_p (beta_j | p_j),

with

    f_pi + f_K + f_p = 1.

The three fractions are the only free parameters in each mixture fit.

Bootstrap replicas provide statistical uncertainties on the fitted species
fractions.

Important interpretation
------------------------
For the pion_selected sample, the fitted f_K is an empirical estimate of
the residual kaon-like component in the tracks that REC classified as pi+
and that passed the actual |chi2pid| < 3.5 requirement.

This is therefore the quantity intended for later folding with the momentum
distribution in each of the 24 final analysis bins.

No contamination estimate is assigned below p = 1.50 GeV in this study.

Outputs
-------
For each period:

    beta_vs_p_<period>.png
    beta_response_<period>.png
    mixture_slices_all_positive_<period>.png
    mixture_slices_pion_selected_<period>.png

    response_points_<period>.csv
    mixture_all_positive_<period>.csv
    mixture_pion_selected_<period>.csv
    kaon_contamination_<period>.csv

The last file is the compact f_K(p) result intended for the eventual
24-bin folding study.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import uproot

from scipy.optimize import minimize


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
    "output/beta_pid_2d_fit"
)


# =============================================================================
# Species definitions
# =============================================================================

SPECIES = (
    "pi",
    "K",
    "p",
)


PID_MAP = {
    "pi": 211,
    "K": 321,
    "p": 2212,
}


MASS = {
    "pi": 0.13957039,
    "K": 0.493677,
    "p": 0.9382720813,
}


# =============================================================================
# Analysis momentum range
# =============================================================================

P_MIN = 1.50
P_MAX = 5.00

P_SLICE_WIDTH = 0.25


P_SLICE_EDGES = np.arange(
    P_MIN,
    P_MAX + P_SLICE_WIDTH,
    P_SLICE_WIDTH,
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


# =============================================================================
# Response calibration
# =============================================================================

# Use somewhat finer bins for measuring detector response.
RESPONSE_P_BIN_WIDTH = 0.10


RESPONSE_P_EDGES = np.arange(
    0.50,
    P_MAX + RESPONSE_P_BIN_WIDTH,
    RESPONSE_P_BIN_WIDTH,
)


# Tight assigned-hypothesis requirement used only to measure the detector
# response around the clearly identified REC species bands.
#
# This is NOT the production pion selection.
RESPONSE_CHI2_MAX = 2.0


# Minimum entries required to retain a response point.
MIN_RESPONSE_ENTRIES = 150


# Robust trimming around the median, in robust sigma units.
RESPONSE_TRIM_NSIGMA = 4.0


# Polynomial order in x = 1/p for the smooth response functions.
RESPONSE_OFFSET_POLY_ORDER = 2
RESPONSE_SIGMA_POLY_ORDER = 2


# Conservative lower/upper limits on response width.
SIGMA_FLOOR = 0.0010
SIGMA_CEILING = 0.0500


# =============================================================================
# Mixture fits
# =============================================================================

MIN_MIXTURE_ENTRIES = 300


# Limit the number of tracks used in a single likelihood fit.
#
# The selection is deterministic through RNG_SEED. This keeps bootstrap
# runtime manageable without changing the sample definition.
MAX_FIT_EVENTS_PER_BIN = 150000


N_BOOTSTRAP = 50

RNG_SEED = 20260927


# =============================================================================
# Plotting
# =============================================================================

BETA_2D_RANGE = (
    0.20,
    1.20,
)


BETA_SLICE_RANGE = (
    0.80,
    1.10,
)


N_BETA_SLICE_BINS = 180


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
# Physics helpers
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


def gaussian_pdf(
    x,
    mu,
    sigma,
):

    sigma = np.clip(
        sigma,
        SIGMA_FLOOR,
        SIGMA_CEILING,
    )


    z = (
        x - mu
    ) / sigma


    return (
        np.exp(
            -0.5 * z * z
        )
        / (
            np.sqrt(
                2.0 * np.pi
            )
            * sigma
        )
    )


def robust_location_scale(
    values,
):

    values = np.asarray(
        values,
        dtype=float,
    )


    values = values[
        np.isfinite(
            values
        )
    ]


    if len(
        values
    ) == 0:

        return (
            np.nan,
            np.nan,
            0,
        )


    median = np.median(
        values
    )


    mad = np.median(
        np.abs(
            values - median
        )
    )


    sigma = (
        1.4826
        * mad
    )


    if (
        not np.isfinite(
            sigma
        )
        or sigma <= 0
    ):

        sigma = np.std(
            values
        )


    if (
        not np.isfinite(
            sigma
        )
        or sigma <= 0
    ):

        return (
            median,
            np.nan,
            len(
                values
            ),
        )


    keep = (
        np.abs(
            values - median
        )
        <= RESPONSE_TRIM_NSIGMA * sigma
    )


    trimmed = values[
        keep
    ]


    if len(
        trimmed
    ) < 10:

        trimmed = values


    location = np.median(
        trimmed
    )


    mad_trimmed = np.median(
        np.abs(
            trimmed - location
        )
    )


    scale = (
        1.4826
        * mad_trimmed
    )


    if (
        not np.isfinite(
            scale
        )
        or scale <= 0
    ):

        scale = np.std(
            trimmed
        )


    return (
        float(
            location
        ),

        float(
            scale
        ),

        int(
            len(
                trimmed
            )
        ),
    )


# =============================================================================
# Input
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
        ],
        dtype=int,
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
            tuple(
                PID_MAP.values()
            ),
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
# Response calibration
# =============================================================================

def measure_response_points(
    period,
    data,
):

    rows = []


    for species in SPECIES:

        pid_value = PID_MAP[
            species
        ]


        mass = MASS[
            species
        ]


        species_selection = (
            (
                data[
                    "pid"
                ]
                == pid_value
            )

            & (
                np.abs(
                    data[
                        "chi2pid"
                    ]
                )
                < RESPONSE_CHI2_MAX
            )
        )


        p_species = data[
            "p"
        ][
            species_selection
        ]


        beta_species = data[
            "beta"
        ][
            species_selection
        ]


        residual_species = (
            beta_species

            - beta_expected(
                p_species,
                mass,
            )
        )


        for index in range(
            len(
                RESPONSE_P_EDGES
            ) - 1
        ):

            p_low = RESPONSE_P_EDGES[
                index
            ]


            p_high = RESPONSE_P_EDGES[
                index + 1
            ]


            selection = (
                (
                    p_species
                    >= p_low
                )

                & (
                    p_species
                    < p_high
                )
            )


            values = residual_species[
                selection
            ]


            if len(
                values
            ) < MIN_RESPONSE_ENTRIES:

                continue


            (
                offset,
                sigma,
                n_used,
            ) = robust_location_scale(
                values
            )


            if (
                not np.isfinite(
                    offset
                )

                or not np.isfinite(
                    sigma
                )

                or sigma <= 0
            ):

                continue


            rows.append(
                {
                    "period":
                        period,

                    "species":
                        species,

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

                    "N_raw":
                        int(
                            len(
                                values
                            )
                        ),

                    "N_used":
                        n_used,

                    "delta_beta":
                        offset,

                    "sigma_beta":
                        sigma,
                }
            )


    return pd.DataFrame(
        rows
    )


# =============================================================================
# Smooth response functions
# =============================================================================

def fit_response_model(
    response_table,
):

    models = {}


    for species in SPECIES:

        table = response_table[
            response_table[
                "species"
            ]
            == species
        ].copy()


        if len(
            table
        ) < 4:

            raise RuntimeError(
                f"Not enough response points for {species}."
            )


        p = table[
            "p_center_GeV"
        ].to_numpy(
            dtype=float
        )


        x = 1.0 / p


        offset = table[
            "delta_beta"
        ].to_numpy(
            dtype=float
        )


        sigma = table[
            "sigma_beta"
        ].to_numpy(
            dtype=float
        )


        n = table[
            "N_used"
        ].to_numpy(
            dtype=float
        )


        weights = np.sqrt(
            np.maximum(
                n,
                1.0,
            )
        )


        offset_coefficients = np.polyfit(
            x,
            offset,

            deg=RESPONSE_OFFSET_POLY_ORDER,

            w=weights,
        )


        # Fit log(sigma) so the smooth model is intrinsically positive.
        sigma_coefficients = np.polyfit(
            x,
            np.log(
                np.clip(
                    sigma,
                    SIGMA_FLOOR,
                    SIGMA_CEILING,
                )
            ),

            deg=RESPONSE_SIGMA_POLY_ORDER,

            w=weights,
        )


        models[
            species
        ] = {
            "offset_coefficients":
                offset_coefficients,

            "sigma_coefficients":
                sigma_coefficients,
        }


    return models


def response_offset(
    momentum,
    species,
    models,
):

    momentum = np.asarray(
        momentum,
        dtype=float,
    )


    x = 1.0 / np.clip(
        momentum,
        0.2,
        None,
    )


    return np.polyval(
        models[
            species
        ][
            "offset_coefficients"
        ],

        x,
    )


def response_sigma(
    momentum,
    species,
    models,
):

    momentum = np.asarray(
        momentum,
        dtype=float,
    )


    x = 1.0 / np.clip(
        momentum,
        0.2,
        None,
    )


    log_sigma = np.polyval(
        models[
            species
        ][
            "sigma_coefficients"
        ],

        x,
    )


    sigma = np.exp(
        log_sigma
    )


    return np.clip(
        sigma,
        SIGMA_FLOOR,
        SIGMA_CEILING,
    )


def response_mean(
    momentum,
    species,
    models,
):

    return (
        beta_expected(
            momentum,
            MASS[
                species
            ],
        )

        + response_offset(
            momentum,
            species,
            models,
        )
    )


# =============================================================================
# Plot detector response
# =============================================================================

def plot_response(
    period,
    response_table,
    models,
):

    fig, axes = plt.subplots(
        2,
        1,

        figsize=(
            8,
            8,
        ),

        sharex=True,
    )


    p_curve = np.linspace(
        0.5,
        5.0,
        500,
    )


    for species in SPECIES:

        table = response_table[
            response_table[
                "species"
            ]
            == species
        ]


        axes[
            0
        ].plot(
            table[
                "p_center_GeV"
            ],

            table[
                "delta_beta"
            ],

            marker="o",

            linestyle="none",

            label=species,
        )


        axes[
            0
        ].plot(
            p_curve,

            response_offset(
                p_curve,
                species,
                models,
            ),

            linewidth=1.5,
        )


        axes[
            1
        ].plot(
            table[
                "p_center_GeV"
            ],

            table[
                "sigma_beta"
            ],

            marker="o",

            linestyle="none",

            label=species,
        )


        axes[
            1
        ].plot(
            p_curve,

            response_sigma(
                p_curve,
                species,
                models,
            ),

            linewidth=1.5,
        )


    axes[
        0
    ].axhline(
        0.0,

        linewidth=0.8,
    )


    axes[
        0
    ].set_ylabel(
        r"$\Delta\beta$"
    )


    axes[
        0
    ].set_title(
        f"{period}: measured TOF response offsets"
    )


    axes[
        0
    ].legend()


    axes[
        1
    ].set_xlabel(
        r"$p$ (GeV)"
    )


    axes[
        1
    ].set_ylabel(
        r"$\sigma_\beta$"
    )


    axes[
        1
    ].set_title(
        "Measured TOF response widths"
    )


    axes[
        1
    ].legend()


    fig.tight_layout()


    output = (
        OUTDIR
        / f"beta_response_{period.lower()}.png"
    )


    fig.savefig(
        output,

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# Two-dimensional beta-vs-p diagnostic
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

    pid = data[
        "pid"
    ]

    chi2pid = data[
        "chi2pid"
    ]


    base = (
        (
            momentum >= 0.5
        )

        & (
            momentum < 5.0
        )

        & (
            beta >= BETA_2D_RANGE[0]
        )

        & (
            beta <= BETA_2D_RANGE[1]
        )
    )


    pion_selected = (
        base

        & (
            pid == 211
        )

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
            base,
            "All positive FD tracks",
        ),

        (
            pion_selected,
            r"REC $\pi^+$, $|\chi^2_{\rm PID}|<3.5$",
        ),
    ]


    p_curve = np.linspace(
        0.5,
        5.0,
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
                180,
                200,
            ],

            range=[
                (
                    0.5,
                    5.0,
                ),

                BETA_2D_RANGE,
            ],

            norm=LogNorm(
                vmin=1
            ),
        )


        fig.colorbar(
            hist[
                3
            ],

            ax=ax,

            label="Counts",
        )


        for species, label in (
            (
                "pi",
                r"$\pi^+$",
            ),

            (
                "K",
                r"$K^+$",
            ),

            (
                "p",
                r"$p$",
            ),
        ):

            ax.plot(
                p_curve,

                beta_expected(
                    p_curve,
                    MASS[
                        species
                    ],
                ),

                linewidth=1.3,

                label=label,
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
        f"{period} positive-hadron PID",
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


# =============================================================================
# Mixture-fraction parameterization
# =============================================================================

def softmax_fractions(
    parameters,
):

    # Fix the pion log-weight to zero.
    logits = np.array(
        [
            0.0,
            parameters[
                0
            ],
            parameters[
                1
            ],
        ],
        dtype=float,
    )


    logits -= np.max(
        logits
    )


    weights = np.exp(
        logits
    )


    fractions = (
        weights
        / np.sum(
            weights
        )
    )


    return fractions


# =============================================================================
# Native (p, beta) mixture likelihood
# =============================================================================

def fit_mixture(
    momentum,
    beta,
    models,
    initial_fractions=None,
):

    momentum = np.asarray(
        momentum,
        dtype=float,
    )


    beta = np.asarray(
        beta,
        dtype=float,
    )


    if len(
        momentum
    ) < MIN_MIXTURE_ENTRIES:

        return None


    pdf_components = []


    for species in SPECIES:

        mu = response_mean(
            momentum,
            species,
            models,
        )


        sigma = response_sigma(
            momentum,
            species,
            models,
        )


        pdf = gaussian_pdf(
            beta,
            mu,
            sigma,
        )


        pdf_components.append(
            pdf
        )


    pdf_components = np.asarray(
        pdf_components
    )


    if initial_fractions is None:

        initial_fractions = np.array(
            [
                0.90,
                0.07,
                0.03,
            ]
        )


    initial_fractions = np.asarray(
        initial_fractions,
        dtype=float,
    )


    initial_fractions = np.clip(
        initial_fractions,
        1.0e-6,
        None,
    )


    initial_fractions /= np.sum(
        initial_fractions
    )


    initial_parameters = np.array(
        [
            np.log(
                initial_fractions[
                    1
                ]
                / initial_fractions[
                    0
                ]
            ),

            np.log(
                initial_fractions[
                    2
                ]
                / initial_fractions[
                    0
                ]
            ),
        ]
    )


    def nll(
        parameters,
    ):

        fractions = softmax_fractions(
            parameters
        )


        total_pdf = np.sum(
            fractions[
                :,
                None,
            ]

            * pdf_components,

            axis=0,
        )


        total_pdf = np.clip(
            total_pdf,
            1.0e-300,
            None,
        )


        return -np.sum(
            np.log(
                total_pdf
            )
        )


    fit = minimize(
        nll,

        initial_parameters,

        method="BFGS",

        options={
            "maxiter":
                1000,

            "gtol":
                1.0e-7,
        },
    )


    fractions = softmax_fractions(
        fit.x
    )


    # =========================================================================
    # Per-event posterior species probabilities
    # =========================================================================

    weighted_pdf = (
        fractions[
            :,
            None,
        ]

        * pdf_components
    )


    denominator = np.sum(
        weighted_pdf,
        axis=0,
    )


    denominator = np.clip(
        denominator,
        1.0e-300,
        None,
    )


    posterior = (
        weighted_pdf
        / denominator[
            None,
            :,
        ]
    )


    return {
        "success":
            bool(
                fit.success
            ),

        "message":
            str(
                fit.message
            ),

        "nll":
            float(
                fit.fun
            ),

        "fractions":
            fractions,

        "posterior":
            posterior,

        "pdf_components":
            pdf_components,
    }


# =============================================================================
# Bootstrap
# =============================================================================

def bootstrap_mixture(
    momentum,
    beta,
    models,
    nominal_fractions,
    rng,
):

    n = len(
        momentum
    )


    replicas = []


    for _ in range(
        N_BOOTSTRAP
    ):

        indices = rng.integers(
            0,
            n,

            size=n,
        )


        result = fit_mixture(
            momentum[
                indices
            ],

            beta[
                indices
            ],

            models,

            initial_fractions=nominal_fractions,
        )


        if result is None:

            continue


        if not np.all(
            np.isfinite(
                result[
                    "fractions"
                ]
            )
        ):

            continue


        replicas.append(
            result[
                "fractions"
            ]
        )


    if len(
        replicas
    ) < 5:

        return {
            "n_bootstrap":
                len(
                    replicas
                ),

            "std":
                np.full(
                    3,
                    np.nan,
                ),

            "p16":
                np.full(
                    3,
                    np.nan,
                ),

            "p84":
                np.full(
                    3,
                    np.nan,
                ),
        }


    replicas = np.asarray(
        replicas
    )


    return {
        "n_bootstrap":
            len(
                replicas
            ),

        "std":
            np.std(
                replicas,

                axis=0,

                ddof=1,
            ),

        "p16":
            np.percentile(
                replicas,
                16,

                axis=0,
            ),

        "p84":
            np.percentile(
                replicas,
                84,

                axis=0,
            ),
    }


# =============================================================================
# Fit all momentum intervals for one sample
# =============================================================================

def fit_sample_by_momentum(
    period,
    sample_name,
    data,
    sample_selection,
    models,
    rng,
):

    rows = []

    plot_payload = []


    for bin_index, (
        p_min,
        p_max,
    ) in enumerate(
        P_RANGES,
        start=1,
    ):

        selection = (
            sample_selection

            & (
                data[
                    "p"
                ]
                >= p_min
            )

            & (
                data[
                    "p"
                ]
                < p_max
            )

            & (
                data[
                    "beta"
                ]
                >= BETA_SLICE_RANGE[0]
            )

            & (
                data[
                    "beta"
                ]
                <= BETA_SLICE_RANGE[1]
            )
        )


        momentum = data[
            "p"
        ][
            selection
        ]


        beta = data[
            "beta"
        ][
            selection
        ]


        n_original = len(
            momentum
        )


        # ---------------------------------------------------------------------
        # Deterministic downsampling for likelihood speed.
        # ---------------------------------------------------------------------

        if (
            MAX_FIT_EVENTS_PER_BIN is not None

            and len(
                momentum
            ) > MAX_FIT_EVENTS_PER_BIN
        ):

            indices = rng.choice(
                len(
                    momentum
                ),

                size=MAX_FIT_EVENTS_PER_BIN,

                replace=False,
            )


            momentum_fit = momentum[
                indices
            ]


            beta_fit = beta[
                indices
            ]


        else:

            momentum_fit = momentum
            beta_fit = beta


        if len(
            momentum_fit
        ) < MIN_MIXTURE_ENTRIES:

            print(
                f"[{period}] {sample_name} "
                f"{p_min:.2f}-{p_max:.2f} GeV: "
                f"only {len(momentum_fit):,} events; skipping.",

                flush=True,
            )

            continue


        result = fit_mixture(
            momentum_fit,
            beta_fit,
            models,
        )


        if result is None:

            continue


        fractions = result[
            "fractions"
        ]


        bootstrap = bootstrap_mixture(
            momentum_fit,
            beta_fit,
            models,
            fractions,
            rng,
        )


        row = {
            "period":
                period,

            "sample":
                sample_name,

            "momentum_bin":
                bin_index,

            "p_min_GeV":
                p_min,

            "p_max_GeV":
                p_max,

            "p_mean_GeV":
                float(
                    np.mean(
                        momentum
                    )
                ),

            "N_total":
                n_original,

            "N_fit":
                len(
                    momentum_fit
                ),

            "fit_success":
                result[
                    "success"
                ],

            "fit_message":
                result[
                    "message"
                ],

            "nll":
                result[
                    "nll"
                ],

            "bootstrap_replicas":
                bootstrap[
                    "n_bootstrap"
                ],
        }


        for species_index, species in enumerate(
            SPECIES
        ):

            row[
                f"f_{species}"
            ] = fractions[
                species_index
            ]


            row[
                f"f_{species}_bootstrap_std"
            ] = bootstrap[
                "std"
            ][
                species_index
            ]


            row[
                f"f_{species}_p16"
            ] = bootstrap[
                "p16"
            ][
                species_index
            ]


            row[
                f"f_{species}_p84"
            ] = bootstrap[
                "p84"
            ][
                species_index
            ]


        rows.append(
            row
        )


        plot_payload.append(
            {
                "p_min":
                    p_min,

                "p_max":
                    p_max,

                "momentum":
                    momentum_fit,

                "beta":
                    beta_fit,

                "result":
                    result,

                "bootstrap":
                    bootstrap,
            }
        )


        print(
            f"[{period}] {sample_name} "
            f"{p_min:.2f}-{p_max:.2f} GeV | "
            f"N={n_original:,} | "
            f"pi={fractions[0]:.4f} | "
            f"K={fractions[1]:.4f} | "
            f"p={fractions[2]:.4f} | "
            f"sigma_K={bootstrap['std'][1]:.4f}",

            flush=True,
        )


    return (
        pd.DataFrame(
            rows
        ),

        plot_payload,
    )


# =============================================================================
# Plot native-fit momentum slices
# =============================================================================

def plot_mixture_slices(
    period,
    sample_name,
    payload,
    models,
):

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


    payload_by_range = {
        (
            item[
                "p_min"
            ],

            item[
                "p_max"
            ],
        ):
            item

        for item in payload
    }


    for panel, (
        p_min,
        p_max,
    ) in enumerate(
        P_RANGES
    ):

        ax = axes[
            panel
        ]


        key = (
            p_min,
            p_max,
        )


        if key not in payload_by_range:

            ax.axis(
                "off"
            )

            continue


        item = payload_by_range[
            key
        ]


        momentum = item[
            "momentum"
        ]


        beta = item[
            "beta"
        ]


        result = item[
            "result"
        ]


        bootstrap = item[
            "bootstrap"
        ]


        fractions = result[
            "fractions"
        ]


        # ---------------------------------------------------------------------
        # Data histogram
        # ---------------------------------------------------------------------

        counts, edges = np.histogram(
            beta,

            bins=N_BETA_SLICE_BINS,

            range=BETA_SLICE_RANGE,
        )


        centers = 0.5 * (
            edges[:-1]
            + edges[1:]
        )


        bin_width = (
            edges[
                1
            ]
            - edges[
                0
            ]
        )


        ax.step(
            centers,
            counts,

            where="mid",

            linewidth=1.0,

            label="Data",
        )


        # ---------------------------------------------------------------------
        # Projection of the native event-level model onto beta.
        #
        # For every event momentum, evaluate the fitted response density on
        # the beta grid, then average over the actual p distribution.
        # ---------------------------------------------------------------------

        total_curve = np.zeros_like(
            centers
        )


        for species_index, species in enumerate(
            SPECIES
        ):

            component_curve = np.zeros_like(
                centers
            )


            # Chunking avoids constructing an unnecessarily huge 2D array.
            chunk_size = 10000


            for start in range(
                0,
                len(
                    momentum
                ),
                chunk_size,
            ):

                stop = min(
                    start + chunk_size,
                    len(
                        momentum
                    ),
                )


                p_chunk = momentum[
                    start:stop
                ]


                mu = response_mean(
                    p_chunk,
                    species,
                    models,
                )


                sigma = response_sigma(
                    p_chunk,
                    species,
                    models,
                )


                z = (
                    centers[
                        None,
                        :
                    ]

                    - mu[
                        :,
                        None,
                    ]
                ) / sigma[
                    :,
                    None,
                ]


                pdf = (
                    np.exp(
                        -0.5
                        * z
                        * z
                    )

                    / (
                        np.sqrt(
                            2.0
                            * np.pi
                        )

                        * sigma[
                            :,
                            None,
                        ]
                    )
                )


                component_curve += np.sum(
                    pdf,

                    axis=0,
                )


            component_curve *= (
                fractions[
                    species_index
                ]

                * bin_width
            )


            total_curve += component_curve


            ax.plot(
                centers,
                component_curve,

                linestyle="--",

                linewidth=1.0,
            )


        ax.plot(
            centers,
            total_curve,

            linewidth=1.5,

            label="Mixture fit",
        )


        # ---------------------------------------------------------------------
        # Annotation
        # ---------------------------------------------------------------------

        k_error = bootstrap[
            "std"
        ][
            1
        ]


        p_error = bootstrap[
            "std"
        ][
            2
        ]


        annotation = (
            f"N = {len(momentum):,}\n"

            f"pi = "
            f"{100.0 * fractions[0]:.2f}%\n"

            f"K = "
            f"{100.0 * fractions[1]:.2f}%"
        )


        if np.isfinite(
            k_error
        ):

            annotation += (
                f" +/- "
                f"{100.0 * k_error:.2f}%"
            )


        annotation += (
            f"\np = "
            f"{100.0 * fractions[2]:.2f}%"
        )


        if np.isfinite(
            p_error
        ):

            annotation += (
                f" +/- "
                f"{100.0 * p_error:.2f}%"
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


        ax.set_title(
            f"{p_min:.2f} < p < {p_max:.2f} GeV",

            fontsize=10,
        )


        ax.set_xlim(
            *BETA_SLICE_RANGE
        )


        ax.set_xlabel(
            r"$\beta$"
        )


        ax.set_ylabel(
            "Counts"
        )


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


    if sample_name == "pion_selected":

        sample_title = (
            r"REC $\pi^+$ with "
            r"$|\chi^2_{\rm PID}|<3.5$"
        )


    else:

        sample_title = (
            "all positive FD tracks"
        )


    fig.suptitle(
        f"{period}: native $(p,\\beta)$ mixture fit, "
        f"{sample_title}",

        fontsize=14,
    )


    fig.tight_layout(
        rect=(
            0,
            0,
            1,
            0.965,
        )
    )


    output = (
        OUTDIR
        / (
            f"mixture_slices_"
            f"{sample_name}_"
            f"{period.lower()}.png"
        )
    )


    fig.savefig(
        output,

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# Compact kaon-contamination output
# =============================================================================

def write_kaon_contamination_table(
    period,
    pion_table,
):

    if len(
        pion_table
    ) == 0:

        return


    columns = [
        "period",
        "momentum_bin",
        "p_min_GeV",
        "p_max_GeV",
        "p_mean_GeV",
        "N_total",
        "f_K",
        "f_K_bootstrap_std",
        "f_K_p16",
        "f_K_p84",
        "f_p",
        "f_p_bootstrap_std",
        "fit_success",
    ]


    output_table = pion_table[
        columns
    ].copy()


    output = (
        OUTDIR
        / f"kaon_contamination_{period.lower()}.csv"
    )


    output_table.to_csv(
        output,

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


    master_rng = np.random.default_rng(
        RNG_SEED
    )


    for period_index, (
        period,
        filename,
    ) in enumerate(
        INPUTS.items()
    ):

        print(
            "\n"
            "============================================================\n"
            f" {period}\n"
            "============================================================",

            flush=True,
        )


        if not filename.is_file():

            raise FileNotFoundError(
                f"Missing calibration ROOT file:\n"
                f"{filename}"
            )


        data = load_period(
            filename
        )


        print(
            f"[{period}] positive FD tracks: "
            f"{len(data['p']):,}",

            flush=True,
        )


        # =====================================================================
        # Basic 2D diagnostic
        # =====================================================================

        plot_beta_vs_p(
            period,
            data,
        )


        # =====================================================================
        # Determine detector response
        # =====================================================================

        response_table = measure_response_points(
            period,
            data,
        )


        response_output = (
            OUTDIR
            / f"response_points_{period.lower()}.csv"
        )


        response_table.to_csv(
            response_output,

            index=False,
        )


        models = fit_response_model(
            response_table
        )


        plot_response(
            period,
            response_table,
            models,
        )


        # =====================================================================
        # Sample 1:
        #
        # all positive FD tracks.
        #
        # This is primarily the validation fit.
        # =====================================================================

        all_positive_selection = (
            (
                data[
                    "p"
                ]
                >= P_MIN
            )

            & (
                data[
                    "p"
                ]
                < P_MAX
            )
        )


        rng_all = np.random.default_rng(
            master_rng.integers(
                0,
                2**32 - 1,
            )
        )


        (
            all_table,
            all_payload,
        ) = fit_sample_by_momentum(
            period,
            "all_positive",
            data,
            all_positive_selection,
            models,
            rng_all,
        )


        all_output = (
            OUTDIR
            / f"mixture_all_positive_{period.lower()}.csv"
        )


        all_table.to_csv(
            all_output,

            index=False,
        )


        plot_mixture_slices(
            period,
            "all_positive",
            all_payload,
            models,
        )


        # =====================================================================
        # Sample 2:
        #
        # Actual pion-selected sample.
        #
        # THIS is the sample from which the residual kaon-like fraction will
        # ultimately be propagated through the final 24-bin momentum spectra.
        # =====================================================================

        pion_selected = (
            (
                data[
                    "pid"
                ]
                == 211
            )

            & (
                np.abs(
                    data[
                        "chi2pid"
                    ]
                )
                < 3.5
            )

            & (
                data[
                    "p"
                ]
                >= P_MIN
            )

            & (
                data[
                    "p"
                ]
                < P_MAX
            )
        )


        rng_pion = np.random.default_rng(
            master_rng.integers(
                0,
                2**32 - 1,
            )
        )


        (
            pion_table,
            pion_payload,
        ) = fit_sample_by_momentum(
            period,
            "pion_selected",
            data,
            pion_selected,
            models,
            rng_pion,
        )


        pion_output = (
            OUTDIR
            / f"mixture_pion_selected_{period.lower()}.csv"
        )


        pion_table.to_csv(
            pion_output,

            index=False,
        )


        plot_mixture_slices(
            period,
            "pion_selected",
            pion_payload,
            models,
        )


        write_kaon_contamination_table(
            period,
            pion_table,
        )


        print(
            f"\n[{period}] completed.",
            flush=True,
        )


    print(
        "\n"
        "============================================================\n"
        " PID 2D-fit study complete\n"
        "============================================================\n"
        f"Output directory: {OUTDIR}",

        flush=True,
    )


if __name__ == "__main__":
    main()