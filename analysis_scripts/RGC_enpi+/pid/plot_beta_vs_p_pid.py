#!/usr/bin/env python3

"""
RGC positive-hadron PID / particle-misidentification study.

This study has two distinct purposes.

(1) Reviewer-requested PID validation
-------------------------------------
Reproduce the beta-versus-momentum comparison requested by the reviewers:

    left:
        all positive FD hadrons before the chi2pid requirement

    right:
        all positive FD hadrons after
            |chi2pid| < 3.5
            0.5 < p < 5.0 GeV

The REC particle assignments pi+, K+, and proton are all retained in both
panels.  This figure belongs in the early PID section of the analysis note.

(2) Residual particle-misidentification study
---------------------------------------------
Estimate whether the sample actually USED as pion candidates,

    REC PID == 211
    |chi2pid| < 3.5
    0.5 < p < 5.0 GeV,

contains residual beta-band components statistically consistent with kaons
or protons.

The contamination fit is performed only for

    1.50 <= p < 5.00 GeV.

Per the analysis prescription, contamination below 1.50 GeV is taken to
be negligible.

IMPORTANT CONCEPTUAL POINT
--------------------------
The same native (p,beta) three-component mixture model is fit twice:

    A. all_positive:
       all positive FD tracks assigned by REC as pi+, K+, or proton.

       This is a validation fit.  The REC assignment is NOT used as the
       species identity in the likelihood.  The likelihood asks how much
       of the measured (p,beta) distribution follows each empirically
       calibrated species response.

    B. pion_selected:
       only tracks satisfying
           REC PID == 211
           |chi2pid| < 3.5.

       AFTER making that selection, the SAME three-component (p,beta)
       mixture model is fit again.

       Thus the fitted kaon fraction is not obtained by counting tracks
       assigned as kaons.  It asks whether any of the tracks that REC
       accepted as pion candidates have measured (p,beta) values that
       statistically follow the calibrated kaon response.

The REC selection and this purity test are not completely independent:
chi2pid itself uses detector PID information including TOF.  The result
should therefore be described as a data-driven residual-band decomposition
of the post-PID pion sample, not as an independent truth tag.

Detector response
-----------------
For each period and species,

    Delta beta_h = beta_measured - beta_expected_h(p),

where

    beta_expected_h(p) = p / sqrt(p^2 + m_h^2).

Clean REC-assigned control tracks are used to determine robust response
points in narrow momentum intervals.

The response is PERIOD SPECIFIC.  This is intentional: Su22 and
Fa22/Sp23 may have different reconstruction software and TOF calibrations.

Rather than imposing a global polynomial in 1/p, the response points are
smoothed directly as functions of momentum using weighted smoothing
splines.  The mixture fit therefore follows the measured detector response
more locally.

Mixture model
-------------
Within each 0.25-GeV momentum interval,

    P(beta_j | p_j) =
          f_pi G_pi(beta_j | p_j)
        + f_K  G_K (beta_j | p_j)
        + f_p  G_p (beta_j | p_j),

with

    f_pi = 1 - f_K - f_p,
    f_K >= 0,
    f_p >= 0,
    f_K + f_p <= 1.

The fit is performed directly with bounded physical fractions.  A zero
kaon or proton fraction is therefore a valid boundary solution and does
not require a softmax parameter to approach minus infinity.

Bootstrap
---------
Bootstrap replicas vary BOTH:

    * the response calibration;
    * the mixture-fit sample.

For every replica:
    1. bootstrap the clean response-control tracks;
    2. rebuild the period-specific response splines;
    3. bootstrap the mixture sample;
    4. refit f_pi, f_K, f_p.

The resulting f_K replica distribution is intended to be propagated later
through the actual momentum distributions of the final 24 analysis bins.

No change to the asymmetry extraction is made by this script.
"""

from pathlib import Path
import warnings

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import uproot

from scipy.interpolate import UnivariateSpline
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
# Production PID momentum region
# =============================================================================

PID_P_MIN = 0.50
PID_P_MAX = 5.00

CHI2PID_MAX = 3.5


# =============================================================================
# Contamination-fit region
# =============================================================================

CONTAMINATION_P_MIN = 1.50
CONTAMINATION_P_MAX = 5.00

P_SLICE_WIDTH = 0.25


P_SLICE_EDGES = np.arange(
    CONTAMINATION_P_MIN,
    CONTAMINATION_P_MAX + P_SLICE_WIDTH,
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

RESPONSE_P_MIN = 0.50
RESPONSE_P_MAX = 5.00
RESPONSE_P_BIN_WIDTH = 0.10


RESPONSE_P_EDGES = np.arange(
    RESPONSE_P_MIN,
    RESPONSE_P_MAX + RESPONSE_P_BIN_WIDTH,
    RESPONSE_P_BIN_WIDTH,
)


# Tight assigned-hypothesis control sample.
#
# This is only used to determine the measured detector response.
RESPONSE_CHI2_MAX = 2.0


MIN_RESPONSE_ENTRIES = 150

RESPONSE_TRIM_NSIGMA = 4.0


SIGMA_FLOOR = 0.0010
SIGMA_CEILING = 0.0500


# Smoothing strength.
#
# The spline target is approximately:
#
#     sum_i [ (y_i - spline_i) / error_i ]^2 ~ N
#
# so a value around 1 gives a modest smoothing appropriate to measured
# response points rather than exact interpolation of every statistical
# fluctuation.
SPLINE_SMOOTHING_SCALE_OFFSET = 1.0
SPLINE_SMOOTHING_SCALE_SIGMA = 1.0


# =============================================================================
# Mixture fits
# =============================================================================

MIN_MIXTURE_ENTRIES = 300


MAX_FIT_EVENTS_PER_BIN = 150000


# Full response+mixture bootstrap.
#
# 50 is a reasonable first production value.  Increase to 100 or 200 for
# the final quoted result after the study is stable.
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

def find_tree(root_file):

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


    if len(values) == 0:

        return (
            np.nan,
            np.nan,
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


    sigma = 1.4826 * mad


    if (
        not np.isfinite(sigma)
        or sigma <= 0
    ):

        sigma = np.std(
            values
        )


    if (
        not np.isfinite(sigma)
        or sigma <= 0
    ):

        return (
            median,
            np.nan,
            np.nan,
            np.nan,
            len(values),
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


    if len(trimmed) < 10:
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
        not np.isfinite(scale)
        or scale <= 0
    ):

        scale = np.std(
            trimmed
        )


    n = len(
        trimmed
    )


    # Approximate uncertainty of a median for a Gaussian-like core.
    location_error = (
        1.2533
        * scale
        / np.sqrt(
            max(n, 1)
        )
    )


    # Approximate statistical uncertainty on a Gaussian width.
    scale_error = (
        scale
        / np.sqrt(
            max(
                2 * (n - 1),
                1,
            )
        )
    )


    return (
        float(location),
        float(scale),
        float(location_error),
        float(scale_error),
        int(n),
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
            ],
            dtype=int,
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
# Response-control samples
# =============================================================================

def get_response_control_sample(
    data,
    species,
):

    selection = (
        (
            data[
                "pid"
            ]
            == PID_MAP[
                species
            ]
        )

        & (
            np.abs(
                data[
                    "chi2pid"
                ]
            )
            < RESPONSE_CHI2_MAX
        )

        & (
            data[
                "p"
            ]
            >= RESPONSE_P_MIN
        )

        & (
            data[
                "p"
            ]
            < RESPONSE_P_MAX
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


    return (
        momentum,
        beta,
    )


# =============================================================================
# Response-point measurement
# =============================================================================

def measure_species_response_points(
    period,
    species,
    momentum,
    beta,
):

    residual = (
        beta

        - beta_expected(
            momentum,
            MASS[
                species
            ],
        )
    )


    rows = []


    for index in range(
        len(RESPONSE_P_EDGES) - 1
    ):

        p_low = RESPONSE_P_EDGES[
            index
        ]


        p_high = RESPONSE_P_EDGES[
            index + 1
        ]


        selection = (
            (
                momentum >= p_low
            )

            & (
                momentum < p_high
            )
        )


        values = residual[
            selection
        ]


        if len(values) < MIN_RESPONSE_ENTRIES:
            continue


        (
            offset,
            sigma,
            offset_error,
            sigma_error,
            n_used,
        ) = robust_location_scale(
            values
        )


        if (
            not np.isfinite(offset)
            or not np.isfinite(sigma)
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
                    float(
                        np.mean(
                            momentum[
                                selection
                            ]
                        )
                    ),

                "N_raw":
                    int(
                        len(values)
                    ),

                "N_used":
                    n_used,

                "delta_beta":
                    offset,

                "delta_beta_error":
                    offset_error,

                "sigma_beta":
                    sigma,

                "sigma_beta_error":
                    sigma_error,
            }
        )


    return pd.DataFrame(
        rows
    )


def measure_response_points(
    period,
    data,
):

    tables = []


    for species in SPECIES:

        (
            momentum,
            beta,
        ) = get_response_control_sample(
            data,
            species,
        )


        table = measure_species_response_points(
            period,
            species,
            momentum,
            beta,
        )


        tables.append(
            table
        )


    return pd.concat(
        tables,
        ignore_index=True,
    )


# =============================================================================
# Smoothed response model
# =============================================================================

class SpeciesResponseModel:

    def __init__(
        self,
        momentum_points,
        offsets,
        offset_errors,
        sigmas,
        sigma_errors,
    ):

        order = np.argsort(
            momentum_points
        )


        p = np.asarray(
            momentum_points,
            dtype=float,
        )[
            order
        ]


        offsets = np.asarray(
            offsets,
            dtype=float,
        )[
            order
        ]


        offset_errors = np.asarray(
            offset_errors,
            dtype=float,
        )[
            order
        ]


        sigmas = np.asarray(
            sigmas,
            dtype=float,
        )[
            order
        ]


        sigma_errors = np.asarray(
            sigma_errors,
            dtype=float,
        )[
            order
        ]


        if len(p) < 4:

            raise RuntimeError(
                "Need at least four response points "
                "for spline response model."
            )


        offset_error_floor = max(
            np.nanmedian(
                offset_errors[
                    np.isfinite(
                        offset_errors
                    )
                    & (
                        offset_errors > 0
                    )
                ]
            ),

            1.0e-5,
        )


        sigma_error_floor = max(
            np.nanmedian(
                sigma_errors[
                    np.isfinite(
                        sigma_errors
                    )
                    & (
                        sigma_errors > 0
                    )
                ]
            ),

            1.0e-5,
        )


        offset_errors = np.where(
            np.isfinite(
                offset_errors
            )
            & (
                offset_errors > 0
            ),

            offset_errors,

            offset_error_floor,
        )


        sigma_errors = np.where(
            np.isfinite(
                sigma_errors
            )
            & (
                sigma_errors > 0
            ),

            sigma_errors,

            sigma_error_floor,
        )


        sigmas = np.clip(
            sigmas,
            SIGMA_FLOOR,
            SIGMA_CEILING,
        )


        # Fit log(sigma) so the smoothed width remains positive.
        log_sigma = np.log(
            sigmas
        )


        log_sigma_errors = (
            sigma_errors
            / sigmas
        )


        log_sigma_errors = np.clip(
            log_sigma_errors,
            1.0e-5,
            None,
        )


        self.p_min = float(
            np.min(p)
        )


        self.p_max = float(
            np.max(p)
        )


        self.offset_spline = UnivariateSpline(
            p,
            offsets,

            w=1.0 / offset_errors,

            s=(
                SPLINE_SMOOTHING_SCALE_OFFSET
                * len(p)
            ),

            k=min(
                3,
                len(p) - 1,
            ),

            ext=3,
        )


        self.log_sigma_spline = UnivariateSpline(
            p,
            log_sigma,

            w=1.0 / log_sigma_errors,

            s=(
                SPLINE_SMOOTHING_SCALE_SIGMA
                * len(p)
            ),

            k=min(
                3,
                len(p) - 1,
            ),

            ext=3,
        )


    def offset(
        self,
        momentum,
    ):

        momentum = np.asarray(
            momentum,
            dtype=float,
        )


        momentum_eval = np.clip(
            momentum,
            self.p_min,
            self.p_max,
        )


        return self.offset_spline(
            momentum_eval
        )


    def sigma(
        self,
        momentum,
    ):

        momentum = np.asarray(
            momentum,
            dtype=float,
        )


        momentum_eval = np.clip(
            momentum,
            self.p_min,
            self.p_max,
        )


        sigma = np.exp(
            self.log_sigma_spline(
                momentum_eval
            )
        )


        return np.clip(
            sigma,
            SIGMA_FLOOR,
            SIGMA_CEILING,
        )


def build_response_models_from_table(
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


        models[
            species
        ] = SpeciesResponseModel(
            table[
                "p_center_GeV"
            ].to_numpy(
                dtype=float
            ),

            table[
                "delta_beta"
            ].to_numpy(
                dtype=float
            ),

            table[
                "delta_beta_error"
            ].to_numpy(
                dtype=float
            ),

            table[
                "sigma_beta"
            ].to_numpy(
                dtype=float
            ),

            table[
                "sigma_beta_error"
            ].to_numpy(
                dtype=float
            ),
        )


    return models


def response_offset(
    momentum,
    species,
    models,
):

    return models[
        species
    ].offset(
        momentum
    )


def response_sigma(
    momentum,
    species,
    models,
):

    return models[
        species
    ].sigma(
        momentum
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
# Response bootstrap
# =============================================================================

def bootstrap_response_models(
    period,
    data,
    rng,
):

    tables = []


    for species in SPECIES:

        (
            momentum,
            beta,
        ) = get_response_control_sample(
            data,
            species,
        )


        if len(momentum) == 0:

            raise RuntimeError(
                f"No response-control tracks for {species}."
            )


        indices = rng.integers(
            0,
            len(momentum),

            size=len(momentum),
        )


        table = measure_species_response_points(
            period,
            species,
            momentum[
                indices
            ],
            beta[
                indices
            ],
        )


        tables.append(
            table
        )


    response_table = pd.concat(
        tables,
        ignore_index=True,
    )


    return build_response_models_from_table(
        response_table
    )


# =============================================================================
# Response plot
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
        RESPONSE_P_MIN,
        RESPONSE_P_MAX,
        700,
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
        ].errorbar(
            table[
                "p_center_GeV"
            ],

            table[
                "delta_beta"
            ],

            yerr=table[
                "delta_beta_error"
            ],

            marker="o",

            markersize=3,

            linestyle="none",

            capsize=1.5,

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
        ].errorbar(
            table[
                "p_center_GeV"
            ],

            table[
                "sigma_beta"
            ],

            yerr=table[
                "sigma_beta_error"
            ],

            marker="o",

            markersize=3,

            linestyle="none",

            capsize=1.5,

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
        f"{period}: period-specific TOF response offset"
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
        "Period-specific TOF response width"
    )


    axes[
        1
    ].legend()


    fig.tight_layout()


    fig.savefig(
        OUTDIR
        / f"beta_response_{period.lower()}.png",

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# 2D beta-versus-p plot helper
# =============================================================================

def draw_beta_panel(
    fig,
    ax,
    momentum,
    beta,
    selection,
    title,
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
                PID_P_MIN,
                PID_P_MAX,
            ),

            BETA_2D_RANGE,
        ],

        norm=LogNorm(
            vmin=1
        ),

        cmap="turbo",
    )


    fig.colorbar(
        hist[
            3
        ],

        ax=ax,

        label="Counts",
    )


    p_curve = np.linspace(
        PID_P_MIN,
        PID_P_MAX,
        700,
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

            linewidth=1.2,

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


# =============================================================================
# Reviewer-requested beta-vs-p figure
# =============================================================================

def plot_beta_vs_p_reviewer(
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


    base = (
        (
            momentum >= PID_P_MIN
        )

        & (
            momentum < PID_P_MAX
        )

        & (
            beta >= BETA_2D_RANGE[0]
        )

        & (
            beta <= BETA_2D_RANGE[1]
        )
    )


    after_pid = (
        base

        & (
            np.abs(
                chi2pid
            )
            < CHI2PID_MAX
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


    draw_beta_panel(
        fig,
        axes[
            0
        ],
        momentum,
        beta,
        base,
        "Before PID cut",
    )


    draw_beta_panel(
        fig,
        axes[
            1
        ],
        momentum,
        beta,
        after_pid,
        (
            r"After $|\chi^2_{\rm PID}|<3.5$, "
            r"$0.5<p<5.0$ GeV"
        ),
    )


    fig.suptitle(
        f"{period}: positive FD hadron PID",
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


    fig.savefig(
        OUTDIR
        / f"beta_vs_p_reviewer_{period.lower()}.png",

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# Pion-selected contamination-study beta-vs-p figure
# =============================================================================

def plot_beta_vs_p_pion_selected(
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
            momentum >= PID_P_MIN
        )

        & (
            momentum < PID_P_MAX
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
            < CHI2PID_MAX
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


    draw_beta_panel(
        fig,
        axes[
            0
        ],
        momentum,
        beta,
        base,
        "All positive FD tracks",
    )


    draw_beta_panel(
        fig,
        axes[
            1
        ],
        momentum,
        beta,
        pion_selected,
        (
            r"REC $\pi^+$ with "
            r"$|\chi^2_{\rm PID}|<3.5$"
        ),
    )


    fig.suptitle(
        f"{period}: pion-selected contamination sample",
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


    fig.savefig(
        OUTDIR
        / f"beta_vs_p_pion_selected_{period.lower()}.png",

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# Physical mixture likelihood
# =============================================================================

def fractions_from_parameters(
    parameters,
):

    f_k = float(
        parameters[
            0
        ]
    )


    f_p = float(
        parameters[
            1
        ]
    )


    f_pi = (
        1.0
        - f_k
        - f_p
    )


    return np.array(
        [
            f_pi,
            f_k,
            f_p,
        ],
        dtype=float,
    )


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


    if len(momentum) < MIN_MIXTURE_ENTRIES:
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


        pdf_components.append(
            gaussian_pdf(
                beta,
                mu,
                sigma,
            )
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
            ],
            dtype=float,
        )


    initial_fractions = np.asarray(
        initial_fractions,
        dtype=float,
    )


    initial_fractions = np.clip(
        initial_fractions,
        0.0,
        1.0,
    )


    initial_fractions /= np.sum(
        initial_fractions
    )


    x0 = np.array(
        [
            initial_fractions[
                1
            ],

            initial_fractions[
                2
            ],
        ],
        dtype=float,
    )


    # -------------------------------------------------------------------------
    # Physical constraints:
    #
    #     f_K >= 0
    #     f_p >= 0
    #     f_K + f_p <= 1
    # -------------------------------------------------------------------------

    bounds = [
        (
            0.0,
            1.0,
        ),

        (
            0.0,
            1.0,
        ),
    ]


    constraints = [
        {
            "type":
                "ineq",

            "fun":
                lambda x:
                    1.0
                    - x[
                        0
                    ]
                    - x[
                        1
                    ],
        }
    ]


    def nll(
        parameters,
    ):

        fractions = fractions_from_parameters(
            parameters
        )


        if np.any(
            fractions < 0
        ):

            return 1.0e100


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


    # Try several starts.  This is cheap because there are only two free
    # fractions and makes the boundary solution robust.
    starting_points = [
        x0,

        np.array(
            [
                0.0,
                0.0,
            ]
        ),

        np.array(
            [
                0.01,
                0.0,
            ]
        ),

        np.array(
            [
                0.05,
                0.01,
            ]
        ),

        np.array(
            [
                0.10,
                0.05,
            ]
        ),
    ]


    fits = []


    for start in starting_points:

        try:

            fit = minimize(
                nll,

                start,

                method="SLSQP",

                bounds=bounds,

                constraints=constraints,

                options={
                    "maxiter":
                        1000,

                    "ftol":
                        1.0e-10,
                },
            )


            if (
                np.isfinite(
                    fit.fun
                )

                and np.all(
                    np.isfinite(
                        fit.x
                    )
                )
            ):

                fits.append(
                    fit
                )


        except Exception:

            continue


    if len(fits) == 0:

        return None


    fit = min(
        fits,
        key=lambda item:
            item.fun,
    )


    fractions = fractions_from_parameters(
        fit.x
    )


    # Numerical cleanup at physical boundaries.
    fractions[
        np.abs(
            fractions
        )
        < 1.0e-12
    ] = 0.0


    fractions = np.clip(
        fractions,
        0.0,
        1.0,
    )


    fractions /= np.sum(
        fractions
    )


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
# Full response + mixture bootstrap
# =============================================================================

def bootstrap_mixture_full(
    period,
    data,
    momentum,
    beta,
    nominal_fractions,
    rng,
):

    n = len(
        momentum
    )


    replicas = []


    for replica_index in range(
        N_BOOTSTRAP
    ):

        # ---------------------------------------------------------------------
        # Rebuild detector response from a bootstrapped control sample.
        # ---------------------------------------------------------------------

        try:

            replica_models = bootstrap_response_models(
                period,
                data,
                rng,
            )


        except Exception as exc:

            warnings.warn(
                f"Response bootstrap replica "
                f"{replica_index} failed: {exc}"
            )

            continue


        # ---------------------------------------------------------------------
        # Bootstrap the mixture sample independently.
        # ---------------------------------------------------------------------

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

            replica_models,

            initial_fractions=nominal_fractions,
        )


        if result is None:
            continue


        fractions = result[
            "fractions"
        ]


        if not np.all(
            np.isfinite(
                fractions
            )
        ):

            continue


        replicas.append(
            fractions
        )


    if len(replicas) == 0:

        return {
            "replicas":
                np.empty(
                    (
                        0,
                        3,
                    )
                ),

            "n_bootstrap":
                0,

            "mean":
                np.full(
                    3,
                    np.nan,
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

            "p50":
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
        replicas,
        dtype=float,
    )


    if len(replicas) > 1:

        std = np.std(
            replicas,

            axis=0,

            ddof=1,
        )


    else:

        std = np.full(
            3,
            np.nan,
        )


    return {
        "replicas":
            replicas,

        "n_bootstrap":
            len(
                replicas
            ),

        "mean":
            np.mean(
                replicas,

                axis=0,
            ),

        "std":
            std,

        "p16":
            np.percentile(
                replicas,
                16,

                axis=0,
            ),

        "p50":
            np.percentile(
                replicas,
                50,

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
# Fit one sample in all momentum intervals
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

    replica_rows = []


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


        if (
            MAX_FIT_EVENTS_PER_BIN is not None

            and len(momentum)
            > MAX_FIT_EVENTS_PER_BIN
        ):

            indices = rng.choice(
                len(momentum),

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


        if len(momentum_fit) < MIN_MIXTURE_ENTRIES:

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

            print(
                f"[{period}] {sample_name} "
                f"{p_min:.2f}-{p_max:.2f} GeV: "
                "fit failed.",

                flush=True,
            )

            continue


        fractions = result[
            "fractions"
        ]


        bootstrap = bootstrap_mixture_full(
            period,
            data,
            momentum_fit,
            beta_fit,
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
                f"f_{species}_bootstrap_mean"
            ] = bootstrap[
                "mean"
            ][
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
                f"f_{species}_p50"
            ] = bootstrap[
                "p50"
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


        for replica_index, replica in enumerate(
            bootstrap[
                "replicas"
            ]
        ):

            replica_rows.append(
                {
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

                    "replica":
                        replica_index,

                    "f_pi":
                        replica[
                            0
                        ],

                    "f_K":
                        replica[
                            1
                        ],

                    "f_p":
                        replica[
                            2
                        ],
                }
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
            f"pi={fractions[0]:.5f} | "
            f"K={fractions[1]:.5f} | "
            f"p={fractions[2]:.5f} | "
            f"K bootstrap="
            f"{bootstrap['p50'][1]:.5f}"
            f" +{bootstrap['p84'][1] - bootstrap['p50'][1]:.5f}"
            f" -{bootstrap['p50'][1] - bootstrap['p16'][1]:.5f}",

            flush=True,
        )


    return (
        pd.DataFrame(
            rows
        ),

        pd.DataFrame(
            replica_rows
        ),

        plot_payload,
    )


# =============================================================================
# Mixture projection plots
# =============================================================================

def project_component_to_beta(
    centers,
    bin_width,
    momentum,
    species,
    fraction,
    models,
):

    curve = np.zeros_like(
        centers,
        dtype=float,
    )


    chunk_size = 10000


    for start in range(
        0,
        len(momentum),
        chunk_size,
    ):

        stop = min(
            start + chunk_size,
            len(momentum),
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


        curve += np.sum(
            pdf,
            axis=0,
        )


    curve *= (
        fraction
        * bin_width
    )


    return curve


def plot_mixture_slices(
    period,
    sample_name,
    payload,
    models,
):

    n_columns = 4


    n_rows = int(
        np.ceil(
            len(P_RANGES)
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


        fractions = item[
            "result"
        ][
            "fractions"
        ]


        bootstrap = item[
            "bootstrap"
        ]


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


        total_curve = np.zeros_like(
            centers
        )


        for species_index, species in enumerate(
            SPECIES
        ):

            component_curve = project_component_to_beta(
                centers,
                bin_width,
                momentum,
                species,
                fractions[
                    species_index
                ],
                models,
            )


            total_curve += component_curve


            ax.plot(
                centers,
                component_curve,

                linestyle="--",

                linewidth=1.0,

                label=species,
            )


        ax.plot(
            centers,
            total_curve,

            linewidth=1.5,

            label="Total",
        )


        k_median = bootstrap[
            "p50"
        ][
            1
        ]


        k_low = (
            k_median
            - bootstrap[
                "p16"
            ][
                1
            ]
        )


        k_high = (
            bootstrap[
                "p84"
            ][
                1
            ]
            - k_median
        )


        annotation = (
            f"N = {len(momentum):,}\n"

            f"pi = "
            f"{100.0 * fractions[0]:.3f}%\n"

            f"K = "
            f"{100.0 * fractions[1]:.3f}%\n"

            f"p = "
            f"{100.0 * fractions[2]:.3f}%"
        )


        if np.isfinite(
            k_median
        ):

            annotation += (
                "\n"
                f"K boot. = "
                f"{100.0 * k_median:.3f}"
                f" +{100.0 * k_high:.3f}"
                f" -{100.0 * k_low:.3f}%"
            )


        ax.text(
            0.03,
            0.96,

            annotation,

            transform=ax.transAxes,

            va="top",
            ha="left",

            fontsize=7.0,
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
        len(P_RANGES),
        len(axes),
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


    fig.savefig(
        OUTDIR
        / (
            f"mixture_slices_"
            f"{sample_name}_"
            f"{period.lower()}.png"
        ),

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# Kaon-contamination summary plot
# =============================================================================

def plot_kaon_contamination(
    period,
    pion_table,
):

    if len(pion_table) == 0:
        return


    p = pion_table[
        "p_mean_GeV"
    ].to_numpy(
        dtype=float
    )


    nominal = pion_table[
        "f_K"
    ].to_numpy(
        dtype=float
    )


    median = pion_table[
        "f_K_p50"
    ].to_numpy(
        dtype=float
    )


    p16 = pion_table[
        "f_K_p16"
    ].to_numpy(
        dtype=float
    )


    p84 = pion_table[
        "f_K_p84"
    ].to_numpy(
        dtype=float
    )


    fig, ax = plt.subplots(
        figsize=(
            8,
            5,
        )
    )


    ax.errorbar(
        p,
        100.0 * median,

        yerr=np.vstack(
            [
                100.0
                * (
                    median - p16
                ),

                100.0
                * (
                    p84 - median
                ),
            ]
        ),

        marker="o",

        linestyle="-",

        capsize=3,

        label="Full bootstrap",
    )


    ax.scatter(
        p,
        100.0 * nominal,

        marker="x",

        label="Nominal response fit",
    )


    ax.axhline(
        0.0,

        linewidth=0.8,
    )


    ax.set_xlim(
        CONTAMINATION_P_MIN,
        CONTAMINATION_P_MAX,
    )


    ax.set_xlabel(
        r"$p$ (GeV)"
    )


    ax.set_ylabel(
        r"Kaon-like fraction (\%)"
    )


    ax.set_title(
        f"{period}: residual kaon-like component "
        "of pion-selected sample"
    )


    ax.legend()


    fig.tight_layout()


    fig.savefig(
        OUTDIR
        / f"kaon_contamination_{period.lower()}.png",

        dpi=200,

        bbox_inches="tight",
    )


    plt.close(
        fig
    )


# =============================================================================
# Compact contamination table
# =============================================================================

def write_kaon_contamination_table(
    period,
    pion_table,
):

    if len(pion_table) == 0:
        return


    columns = [
        "period",
        "momentum_bin",
        "p_min_GeV",
        "p_max_GeV",
        "p_mean_GeV",
        "N_total",
        "f_K",
        "f_K_bootstrap_mean",
        "f_K_bootstrap_std",
        "f_K_p16",
        "f_K_p50",
        "f_K_p84",
        "f_p",
        "f_p_bootstrap_mean",
        "f_p_bootstrap_std",
        "f_p_p16",
        "f_p_p50",
        "f_p_p84",
        "fit_success",
    ]


    pion_table[
        columns
    ].to_csv(
        OUTDIR
        / f"kaon_contamination_{period.lower()}.csv",

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


    for period, filename in INPUTS.items():

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
        # Reviewer-requested PID figure.
        # =====================================================================

        plot_beta_vs_p_reviewer(
            period,
            data,
        )


        # =====================================================================
        # Separate pion-selected contamination-study figure.
        # =====================================================================

        plot_beta_vs_p_pion_selected(
            period,
            data,
        )


        # =====================================================================
        # Period-specific detector response.
        # =====================================================================

        response_table = measure_response_points(
            period,
            data,
        )


        response_table.to_csv(
            OUTDIR
            / f"response_points_{period.lower()}.csv",

            index=False,
        )


        models = build_response_models_from_table(
            response_table
        )


        plot_response(
            period,
            response_table,
            models,
        )


        # =====================================================================
        # Validation fit:
        # all positive tracks.
        # =====================================================================

        all_positive_selection = (
            (
                data[
                    "p"
                ]
                >= CONTAMINATION_P_MIN
            )

            & (
                data[
                    "p"
                ]
                < CONTAMINATION_P_MAX
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
            all_replicas,
            all_payload,
        ) = fit_sample_by_momentum(
            period,
            "all_positive",
            data,
            all_positive_selection,
            models,
            rng_all,
        )


        all_table.to_csv(
            OUTDIR
            / f"mixture_all_positive_{period.lower()}.csv",

            index=False,
        )


        all_replicas.to_csv(
            OUTDIR
            / f"mixture_all_positive_replicas_{period.lower()}.csv",

            index=False,
        )


        plot_mixture_slices(
            period,
            "all_positive",
            all_payload,
            models,
        )


        # =====================================================================
        # Physics sample for residual misidentification:
        #
        # First make the actual pion-candidate selection.
        #
        # THEN fit the surviving tracks with exactly the same pi/K/p
        # (p,beta) mixture model.
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
                < CHI2PID_MAX
            )

            & (
                data[
                    "p"
                ]
                >= CONTAMINATION_P_MIN
            )

            & (
                data[
                    "p"
                ]
                < CONTAMINATION_P_MAX
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
            pion_replicas,
            pion_payload,
        ) = fit_sample_by_momentum(
            period,
            "pion_selected",
            data,
            pion_selected,
            models,
            rng_pion,
        )


        pion_table.to_csv(
            OUTDIR
            / f"mixture_pion_selected_{period.lower()}.csv",

            index=False,
        )


        pion_replicas.to_csv(
            OUTDIR
            / f"mixture_pion_selected_replicas_{period.lower()}.csv",

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


        plot_kaon_contamination(
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
        " PID / particle-misidentification study complete\n"
        "============================================================\n"
        f"Output directory: {OUTDIR}",

        flush=True,
    )


if __name__ == "__main__":

    main()