#!/usr/bin/env python3

"""
RGC positive-hadron PID / particle-misidentification study.

Final purpose
-------------
1. Produce the reviewer-requested beta-versus-p PID validation plots.

2. Use a period-specific native (p,beta) mixture fit to estimate the
   residual kaon-like component of the actual pion-selected sample:

       REC PID == 211
       |chi2pid| < 3.5
       no explicit pion-momentum cut

3. Produce a compact, portable f_K(period,p) lookup and a coherent set of
   bootstrap replicas that can be loaded by the final kinematic-distribution
   code and folded through the pion momentum distribution in each of the
   24 final (x_B,-t') analysis bins.

Physics prescription
--------------------
For p < 1.50 GeV:

    f_K(p) = 0

by the adopted low-momentum prescription, where the pion/kaon TOF separation is strong.

For

    1.50 <= p < 10.50 GeV,

the pion-selected sample is fit in 0.25-GeV momentum intervals with

    P(beta_j | p_j)
      = f_pi G_pi(beta_j | p_j)
      + f_K  G_K (beta_j | p_j)
      + f_p  G_p(beta_j | p_j),

where every event is evaluated at its own momentum.

The previous 5-GeV upper boundary has been removed.  The fit is extended to
10.5 GeV, covering the physically populated RGC momentum range.  This 10.5-GeV
value is a numerical study boundary near the beam energy, not an analysis cut;
tracks above it are reported separately as pathological/overflow entries and
are not used to define the PID response.

Detector response
-----------------
For each period and species,

    Delta beta_h = beta_measured - beta_expected_h(p),

with

    beta_expected_h(p) = p / sqrt(p^2 + m_h^2).

The detector response is determined separately for Su22, Fa22, and Sp23.
This is intentional because the periods need not share identical
reconstruction/calibration conditions.

The response points are described by moderately smoothed splines.  The
smoothing is deliberately stronger than in the previous version so that
individual low-statistics high-p response points cannot generate sharp
spline excursions near 5 GeV.

Bootstrap
---------
The final contamination bootstrap is coherent across momentum intervals.

For bootstrap replica b:

    1. one fluctuated period-specific response model is generated;
    2. every pion-selected momentum interval is independently resampled;
    3. f_K is refit in every interval using that same response replica.

Thus replica b represents one complete possible f_K(p) curve.

This is important for the later 24-bin propagation.

The bootstrap replicas are parallelized over 8 worker processes.

Outputs
-------
For each period:

    beta_vs_p_reviewer_<period>.png
    beta_vs_p_pion_selected_<period>.png
    beta_response_<period>.png

    mixture_slices_all_positive_<period>.png
    mixture_slices_pion_selected_<period>.png
    kaon_contamination_<period>.png

    response_points_<period>.csv
    mixture_all_positive_<period>.csv
    mixture_pion_selected_<period>.csv

    kaon_contamination_lookup_<period>.csv
    kaon_contamination_replicas_<period>.csv
    kaon_contamination_replica_matrix_<period>.csv

Combined portable outputs:

    kaon_contamination_lookup_all_periods.csv
    kaon_contamination_replicas_all_periods.csv

The two combined files are the intended interface to the later
kinematic-distribution / 24-bin contamination calculation.

No change to the production asymmetry extraction is made here.
"""

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
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
    "output/beta_pid_2d_fit_full_momentum"
)


N_WORKERS = 8

N_BOOTSTRAP = 200

RNG_SEED = 20260927


# =============================================================================
# Species
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
# Production PID region
# =============================================================================

PID_P_MIN = 0.00
PID_P_MAX = 10.50

CHI2PID_MAX = 3.5


# =============================================================================
# Contamination region
# =============================================================================

CONTAMINATION_P_MIN = 1.50
CONTAMINATION_P_MAX = 10.50

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


# Explicit zero-contamination bins below 1.5 GeV.
LOW_P_EDGES = np.arange(
    PID_P_MIN,
    CONTAMINATION_P_MIN + P_SLICE_WIDTH,
    P_SLICE_WIDTH,
)


LOW_P_RANGES = [
    (
        float(LOW_P_EDGES[index]),
        float(LOW_P_EDGES[index + 1]),
    )
    for index in range(
        len(LOW_P_EDGES) - 1
    )
]


# =============================================================================
# Response calibration
# =============================================================================

RESPONSE_P_MIN = 0.20
RESPONSE_P_MAX = 10.50
RESPONSE_P_BIN_WIDTH = 0.10


RESPONSE_P_EDGES = np.arange(
    RESPONSE_P_MIN,
    RESPONSE_P_MAX + RESPONSE_P_BIN_WIDTH,
    RESPONSE_P_BIN_WIDTH,
)


RESPONSE_CHI2_MAX = 2.0

MIN_RESPONSE_ENTRIES = 150

RESPONSE_TRIM_NSIGMA = 4.0


SIGMA_FLOOR = 0.0010
SIGMA_CEILING = 0.0500


# Stronger smoothing than the previous iteration.
#
# The previous value was 1.0.  Four times the nominal smoothing suppresses
# high-p edge wiggles while retaining the broad measured momentum dependence.
SPLINE_SMOOTHING_SCALE_OFFSET = 4.0
SPLINE_SMOOTHING_SCALE_SIGMA = 4.0


# =============================================================================
# Mixture fitting
# =============================================================================

MIN_MIXTURE_ENTRIES = 300

MAX_FIT_EVENTS_PER_BIN = 150000


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
# Worker globals
#
# On the JLab Linux farm we use fork.  The parent loads the arrays once and
# the worker processes inherit them copy-on-write instead of serializing the
# full calibration sample for every bootstrap job.
# =============================================================================

_WORKER_PERIOD = None
_WORKER_DATA = None
_WORKER_PION_BIN_DATA = None
_WORKER_RESPONSE_TABLE = None


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


    sigma = (
        1.4826
        * mad
    )


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


    location_error = (
        1.2533
        * scale
        / np.sqrt(
            max(n, 1)
        )
    )


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
# Response samples
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


    return (
        data[
            "p"
        ][
            selection
        ],

        data[
            "beta"
        ][
            selection
        ],
    )


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


        tables.append(
            measure_species_response_points(
                period,
                species,
                momentum,
                beta,
            )
        )


    return pd.concat(
        tables,
        ignore_index=True,
    )


# =============================================================================
# Response model
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
                "Need at least four response points."
            )


        valid_offset_errors = offset_errors[
            np.isfinite(
                offset_errors
            )
            & (
                offset_errors > 0
            )
        ]


        valid_sigma_errors = sigma_errors[
            np.isfinite(
                sigma_errors
            )
            & (
                sigma_errors > 0
            )
        ]


        offset_error_floor = max(
            np.nanmedian(
                valid_offset_errors
            )
            if len(valid_offset_errors)
            else 1.0e-5,

            1.0e-5,
        )


        sigma_error_floor = max(
            np.nanmedian(
                valid_sigma_errors
            )
            if len(valid_sigma_errors)
            else 1.0e-5,

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
# Efficient response-replica construction
#
# We already measured robust response points and their statistical errors.
# For the final bootstrap we fluctuate those measured points rather than
# repeatedly resampling millions of control tracks.
#
# This retains the measured response uncertainty while making 200 coherent
# replicas computationally practical.
# =============================================================================

def fluctuate_response_table(
    response_table,
    rng,
):

    replica = response_table.copy()


    replica[
        "delta_beta"
    ] = rng.normal(
        response_table[
            "delta_beta"
        ].to_numpy(
            dtype=float
        ),

        response_table[
            "delta_beta_error"
        ].to_numpy(
            dtype=float
        ),
    )


    sigma_nominal = response_table[
        "sigma_beta"
    ].to_numpy(
        dtype=float
    )


    sigma_error = response_table[
        "sigma_beta_error"
    ].to_numpy(
        dtype=float
    )


    sigma_replica = rng.normal(
        sigma_nominal,
        sigma_error,
    )


    replica[
        "sigma_beta"
    ] = np.clip(
        sigma_replica,
        SIGMA_FLOOR,
        SIGMA_CEILING,
    )


    return replica


# =============================================================================
# Response plotting
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
# 2D beta-p plotting
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
            rf"$0<p<{PID_P_MAX:.1f}$ GeV"
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
# Mixture likelihood
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
                0.95,
                0.04,
                0.01,
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
                0.002,
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
    }


# =============================================================================
# Prepare momentum-bin samples
# =============================================================================

def prepare_pion_bin_data(
    data,
):

    base = (
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


    result = []


    for bin_index, (
        p_min,
        p_max,
    ) in enumerate(
        P_RANGES,
        start=1,
    ):

        selection = (
            base

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


        result.append(
            {
                "momentum_bin":
                    bin_index,

                "p_min":
                    p_min,

                "p_max":
                    p_max,

                "p":
                    np.asarray(
                        data[
                            "p"
                        ][
                            selection
                        ],
                        dtype=float,
                    ),

                "beta":
                    np.asarray(
                        data[
                            "beta"
                        ][
                            selection
                        ],
                        dtype=float,
                    ),
            }
        )


    return result


# =============================================================================
# Nominal sample fits
# =============================================================================

def fit_nominal_sample(
    period,
    sample_name,
    data,
    sample_selection,
    models,
    rng,
):

    rows = []

    payload = []


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

            and n_original
            > MAX_FIT_EVENTS_PER_BIN
        ):

            indices = rng.choice(
                n_original,

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
            continue


        fit = fit_mixture(
            momentum_fit,
            beta_fit,
            models,
        )


        if fit is None:
            continue


        fractions = fit[
            "fractions"
        ]


        rows.append(
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
                    fit[
                        "success"
                    ],

                "fit_message":
                    fit[
                        "message"
                    ],

                "nll":
                    fit[
                        "nll"
                    ],

                "f_pi":
                    fractions[
                        0
                    ],

                "f_K":
                    fractions[
                        1
                    ],

                "f_p":
                    fractions[
                        2
                    ],
            }
        )


        payload.append(
            {
                "p_min":
                    p_min,

                "p_max":
                    p_max,

                "momentum":
                    momentum_fit,

                "beta":
                    beta_fit,

                "fractions":
                    fractions,
            }
        )


        print(
            f"[{period}] {sample_name} "
            f"{p_min:.2f}-{p_max:.2f} GeV | "
            f"N={n_original:,} | "
            f"pi={fractions[0]:.5f} | "
            f"K={fractions[1]:.5f} | "
            f"p={fractions[2]:.5f}",

            flush=True,
        )


    return (
        pd.DataFrame(
            rows
        ),

        payload,
    )


# =============================================================================
# Parallel coherent bootstrap
# =============================================================================

def initialize_worker_globals(
    period,
    data,
    pion_bin_data,
    response_table,
):

    global _WORKER_PERIOD
    global _WORKER_DATA
    global _WORKER_PION_BIN_DATA
    global _WORKER_RESPONSE_TABLE


    _WORKER_PERIOD = period
    _WORKER_DATA = data
    _WORKER_PION_BIN_DATA = pion_bin_data
    _WORKER_RESPONSE_TABLE = response_table


def bootstrap_replica_worker(
    replica_index,
    seed,
):

    rng = np.random.default_rng(
        seed
    )


    # -------------------------------------------------------------------------
    # One detector-response realization for the ENTIRE f_K(p) curve.
    # -------------------------------------------------------------------------

    response_replica = fluctuate_response_table(
        _WORKER_RESPONSE_TABLE,
        rng,
    )


    models = build_response_models_from_table(
        response_replica
    )


    rows = []


    for item in _WORKER_PION_BIN_DATA:

        momentum = item[
            "p"
        ]


        beta = item[
            "beta"
        ]


        n = len(
            momentum
        )


        if n < MIN_MIXTURE_ENTRIES:

            rows.append(
                {
                    "period":
                        _WORKER_PERIOD,

                    "replica":
                        replica_index,

                    "momentum_bin":
                        item[
                            "momentum_bin"
                        ],

                    "p_min_GeV":
                        item[
                            "p_min"
                        ],

                    "p_max_GeV":
                        item[
                            "p_max"
                        ],

                    "f_pi":
                        np.nan,

                    "f_K":
                        np.nan,

                    "f_p":
                        np.nan,
                }
            )

            continue


        # ---------------------------------------------------------------------
        # Bootstrap the pion-selected tracks in this momentum interval.
        # ---------------------------------------------------------------------

        indices = rng.integers(
            0,
            n,

            size=n,
        )


        momentum_replica = momentum[
            indices
        ]


        beta_replica = beta[
            indices
        ]


        if (
            MAX_FIT_EVENTS_PER_BIN is not None

            and len(
                momentum_replica
            ) > MAX_FIT_EVENTS_PER_BIN
        ):

            sub_indices = rng.choice(
                len(
                    momentum_replica
                ),

                size=MAX_FIT_EVENTS_PER_BIN,

                replace=False,
            )


            momentum_replica = momentum_replica[
                sub_indices
            ]


            beta_replica = beta_replica[
                sub_indices
            ]


        fit = fit_mixture(
            momentum_replica,
            beta_replica,
            models,
        )


        if fit is None:

            fractions = np.full(
                3,
                np.nan,
            )


        else:

            fractions = fit[
                "fractions"
            ]


        rows.append(
            {
                "period":
                    _WORKER_PERIOD,

                "replica":
                    replica_index,

                "momentum_bin":
                    item[
                        "momentum_bin"
                    ],

                "p_min_GeV":
                    item[
                        "p_min"
                    ],

                "p_max_GeV":
                    item[
                        "p_max"
                    ],

                "f_pi":
                    fractions[
                        0
                    ],

                "f_K":
                    fractions[
                        1
                    ],

                "f_p":
                    fractions[
                        2
                    ],
            }
        )


    return rows


def run_parallel_bootstrap(
    period,
    data,
    pion_bin_data,
    response_table,
):

    print(
        f"\n[{period}] Starting {N_BOOTSTRAP} coherent bootstrap replicas "
        f"with {N_WORKERS} workers...",

        flush=True,
    )


    seed_sequence = np.random.SeedSequence(
        RNG_SEED
        + {
            "Su22":
                1000,

            "Fa22":
                2000,

            "Sp23":
                3000,
        }[
            period
        ]
    )


    child_sequences = seed_sequence.spawn(
        N_BOOTSTRAP
    )


    seeds = [
        int(
            sequence.generate_state(
                1
            )[
                0
            ]
        )
        for sequence in child_sequences
    ]


    all_rows = []


    # JLab ifarm is Linux; fork allows the large numpy arrays to be inherited
    # copy-on-write rather than repeatedly serialized.
    context = mp.get_context(
        "fork"
    )


    initialize_worker_globals(
        period,
        data,
        pion_bin_data,
        response_table,
    )


    with ProcessPoolExecutor(
        max_workers=N_WORKERS,
        mp_context=context,
    ) as executor:

        futures = {
            executor.submit(
                bootstrap_replica_worker,
                replica_index,
                seeds[
                    replica_index
                ],
            ):
                replica_index

            for replica_index in range(
                N_BOOTSTRAP
            )
        }


        completed = 0


        for future in as_completed(
            futures
        ):

            replica_index = futures[
                future
            ]


            try:

                rows = future.result()


            except Exception as exc:

                warnings.warn(
                    f"{period} bootstrap replica "
                    f"{replica_index} failed: {exc}"
                )

                continue


            all_rows.extend(
                rows
            )


            completed += 1


            if (
                completed == 1

                or completed % 10 == 0

                or completed == N_BOOTSTRAP
            ):

                print(
                    f"[{period}] bootstrap: "
                    f"{completed}/{N_BOOTSTRAP} replicas complete",

                    flush=True,
                )


    replicas = pd.DataFrame(
        all_rows
    )


    if len(replicas):

        replicas = replicas.sort_values(
            [
                "replica",
                "momentum_bin",
            ]
        ).reset_index(
            drop=True
        )


    return replicas


# =============================================================================
# Bootstrap summaries
# =============================================================================

def summarize_bootstrap(
    nominal_table,
    replicas,
):

    rows = []


    for _, nominal in nominal_table.iterrows():

        momentum_bin = int(
            nominal[
                "momentum_bin"
            ]
        )


        subset = replicas[
            replicas[
                "momentum_bin"
            ]
            == momentum_bin
        ]


        row = nominal.to_dict()


        for species in SPECIES:

            values = subset[
                f"f_{species}"
            ].to_numpy(
                dtype=float
            )


            values = values[
                np.isfinite(
                    values
                )
            ]


            if len(values) == 0:

                row[
                    f"f_{species}_bootstrap_mean"
                ] = np.nan


                row[
                    f"f_{species}_bootstrap_std"
                ] = np.nan


                row[
                    f"f_{species}_p16"
                ] = np.nan


                row[
                    f"f_{species}_p50"
                ] = np.nan


                row[
                    f"f_{species}_p84"
                ] = np.nan


            else:

                row[
                    f"f_{species}_bootstrap_mean"
                ] = float(
                    np.mean(
                        values
                    )
                )


                row[
                    f"f_{species}_bootstrap_std"
                ] = float(
                    np.std(
                        values,
                        ddof=1,
                    )
                ) if len(values) > 1 else np.nan


                (
                    row[
                        f"f_{species}_p16"
                    ],

                    row[
                        f"f_{species}_p50"
                    ],

                    row[
                        f"f_{species}_p84"
                    ],

                ) = np.percentile(
                    values,
                    [
                        16,
                        50,
                        84,
                    ],
                )


        row[
            "bootstrap_replicas"
        ] = int(
            subset[
                "replica"
            ].nunique()
        )


        rows.append(
            row
        )


    return pd.DataFrame(
        rows
    )


# =============================================================================
# Mixture plots
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
    summary_table=None,
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
            "fractions"
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


        annotation = (
            f"N = {len(momentum):,}\n"

            f"pi = "
            f"{100.0 * fractions[0]:.3f}%\n"

            f"K = "
            f"{100.0 * fractions[1]:.3f}%\n"

            f"p = "
            f"{100.0 * fractions[2]:.3f}%"
        )


        if (
            summary_table is not None
            and len(summary_table)
        ):

            match = summary_table[
                np.isclose(
                    summary_table[
                        "p_min_GeV"
                    ],
                    p_min,
                )

                & np.isclose(
                    summary_table[
                        "p_max_GeV"
                    ],
                    p_max,
                )
            ]


            if len(match):

                match = match.iloc[
                    0
                ]


                median = match[
                    "f_K_p50"
                ]


                p16 = match[
                    "f_K_p16"
                ]


                p84 = match[
                    "f_K_p84"
                ]


                annotation += (
                    "\n"
                    f"K boot. = "
                    f"{100.0 * median:.3f}"
                    f" +{100.0 * (p84 - median):.3f}"
                    f" -{100.0 * (median - p16):.3f}%"
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
# Portable contamination lookup
# =============================================================================

def build_lookup_table(
    period,
    pion_summary,
):

    rows = []


    # -------------------------------------------------------------------------
    # Explicit adopted zero-contamination region.
    # -------------------------------------------------------------------------

    for p_min, p_max in LOW_P_RANGES:

        rows.append(
            {
                "period":
                    period,

                "p_min_GeV":
                    p_min,

                "p_max_GeV":
                    p_max,

                "p_center_GeV":
                    0.5
                    * (
                        p_min
                        + p_max
                    ),

                "N_calibration":
                    np.nan,

                "f_K_nominal":
                    0.0,

                "f_K_mean":
                    0.0,

                "f_K_std":
                    0.0,

                "f_K_p16":
                    0.0,

                "f_K_p50":
                    0.0,

                "f_K_p84":
                    0.0,

                "prescription":
                    "assumed_zero_below_1p5_GeV",
            }
        )


    # -------------------------------------------------------------------------
    # Fitted region.
    # -------------------------------------------------------------------------

    for _, row in pion_summary.iterrows():

        rows.append(
            {
                "period":
                    period,

                "p_min_GeV":
                    row[
                        "p_min_GeV"
                    ],

                "p_max_GeV":
                    row[
                        "p_max_GeV"
                    ],

                "p_center_GeV":
                    row[
                        "p_mean_GeV"
                    ],

                "N_calibration":
                    row[
                        "N_total"
                    ],

                "f_K_nominal":
                    row[
                        "f_K"
                    ],

                "f_K_mean":
                    row[
                        "f_K_bootstrap_mean"
                    ],

                "f_K_std":
                    row[
                        "f_K_bootstrap_std"
                    ],

                "f_K_p16":
                    row[
                        "f_K_p16"
                    ],

                "f_K_p50":
                    row[
                        "f_K_p50"
                    ],

                "f_K_p84":
                    row[
                        "f_K_p84"
                    ],

                "prescription":
                    "native_p_beta_mixture_fit",
            }
        )


    lookup = pd.DataFrame(
        rows
    )


    return lookup.sort_values(
        "p_min_GeV"
    ).reset_index(
        drop=True
    )


def build_portable_replica_table(
    period,
    replicas,
):

    rows = []


    replica_ids = sorted(
        replicas[
            "replica"
        ].unique()
    )


    for replica_id in replica_ids:

        # ---------------------------------------------------------------------
        # Explicit zero region.
        # ---------------------------------------------------------------------

        for p_min, p_max in LOW_P_RANGES:

            rows.append(
                {
                    "period":
                        period,

                    "replica":
                        int(
                            replica_id
                        ),

                    "p_min_GeV":
                        p_min,

                    "p_max_GeV":
                        p_max,

                    "f_K":
                        0.0,

                    "prescription":
                        "assumed_zero_below_1p5_GeV",
                }
            )


        # ---------------------------------------------------------------------
        # Fitted region.
        # ---------------------------------------------------------------------

        subset = replicas[
            replicas[
                "replica"
            ]
            == replica_id
        ]


        for _, row in subset.iterrows():

            rows.append(
                {
                    "period":
                        period,

                    "replica":
                        int(
                            replica_id
                        ),

                    "p_min_GeV":
                        row[
                            "p_min_GeV"
                        ],

                    "p_max_GeV":
                        row[
                            "p_max_GeV"
                        ],

                    "f_K":
                        row[
                            "f_K"
                        ],

                    "prescription":
                        "native_p_beta_mixture_fit",
                }
            )


    result = pd.DataFrame(
        rows
    )


    return result.sort_values(
        [
            "replica",
            "p_min_GeV",
        ]
    ).reset_index(
        drop=True
    )


def write_replica_matrix(
    period,
    portable_replicas,
):

    matrix = portable_replicas.pivot(
        index="replica",
        columns=[
            "p_min_GeV",
            "p_max_GeV",
        ],
        values="f_K",
    )


    matrix.columns = [
        (
            f"fK_"
            f"{p_min:.2f}_"
            f"{p_max:.2f}_GeV"
        )
        for p_min, p_max in matrix.columns
    ]


    matrix = matrix.reset_index()


    matrix.insert(
        0,
        "period",
        period,
    )


    matrix.to_csv(
        OUTDIR
        / (
            f"kaon_contamination_replica_matrix_"
            f"{period.lower()}.csv"
        ),

        index=False,
    )


# =============================================================================
# Contamination plot
# =============================================================================

def plot_kaon_contamination(
    period,
    lookup,
):

    fitted = lookup[
        lookup[
            "prescription"
        ]
        == "native_p_beta_mixture_fit"
    ]


    fig, ax = plt.subplots(
        figsize=(
            8,
            5,
        )
    )


    p = fitted[
        "p_center_GeV"
    ].to_numpy(
        dtype=float
    )


    median = fitted[
        "f_K_p50"
    ].to_numpy(
        dtype=float
    )


    p16 = fitted[
        "f_K_p16"
    ].to_numpy(
        dtype=float
    )


    p84 = fitted[
        "f_K_p84"
    ].to_numpy(
        dtype=float
    )


    nominal = fitted[
        "f_K_nominal"
    ].to_numpy(
        dtype=float
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

        label="Bootstrap median and 68% interval",
    )


    ax.scatter(
        p,
        100.0 * nominal,

        marker="x",

        label="Nominal fit",
    )


    ax.plot(
        [
            PID_P_MIN,
            CONTAMINATION_P_MIN,
        ],

        [
            0.0,
            0.0,
        ],

        linewidth=2.0,

        label=r"Adopted $f_K=0$ below 1.5 GeV",
    )


    ax.axhline(
        0.0,

        linewidth=0.8,
    )


    ax.set_xlim(
        PID_P_MIN,
        PID_P_MAX,
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


    ax.legend(
        fontsize=8,
    )


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


    combined_lookup_tables = []

    combined_replica_tables = []


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


        # Full-momentum coverage diagnostic.  The 10.5-GeV endpoint is a
        # numerical PID-study boundary, not a physics-selection cut.
        pvals = data["p"]
        n_all = len(pvals)
        print(
            f"[{period}] momentum coverage: "
            f"p>5: {np.count_nonzero(pvals > 5.0):,} "
            f"({100.0*np.count_nonzero(pvals > 5.0)/max(n_all,1):.2f}%), "
            f"p>8: {np.count_nonzero(pvals > 8.0):,} "
            f"({100.0*np.count_nonzero(pvals > 8.0)/max(n_all,1):.2f}%), "
            f"p>=10.5 overflow: {np.count_nonzero(pvals >= PID_P_MAX):,}",
            flush=True,
        )


        # =====================================================================
        # Reviewer figure.
        # =====================================================================

        plot_beta_vs_p_reviewer(
            period,
            data,
        )


        # =====================================================================
        # Pion-selected study figure.
        # =====================================================================

        plot_beta_vs_p_pion_selected(
            period,
            data,
        )


        # =====================================================================
        # Nominal response.
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
        # All-positive validation.
        #
        # Nominal only: no need to spend 200 bootstrap replicas on a control
        # sample that is not propagated into the final contamination estimate.
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
            all_payload,
        ) = fit_nominal_sample(
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


        plot_mixture_slices(
            period,
            "all_positive",
            all_payload,
            models,
            summary_table=None,
        )


        # =====================================================================
        # Nominal pion-selected fits.
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
            pion_nominal,
            pion_payload,
        ) = fit_nominal_sample(
            period,
            "pion_selected",
            data,
            pion_selected,
            models,
            rng_pion,
        )


        # =====================================================================
        # Coherent 8-worker bootstrap.
        # =====================================================================

        pion_bin_data = prepare_pion_bin_data(
            data
        )


        replicas = run_parallel_bootstrap(
            period,
            data,
            pion_bin_data,
            response_table,
        )


        pion_summary = summarize_bootstrap(
            pion_nominal,
            replicas,
        )


        pion_summary.to_csv(
            OUTDIR
            / f"mixture_pion_selected_{period.lower()}.csv",

            index=False,
        )


        # =====================================================================
        # Updated mixture plots with bootstrap intervals.
        # =====================================================================

        plot_mixture_slices(
            period,
            "pion_selected",
            pion_payload,
            models,
            summary_table=pion_summary,
        )


        # =====================================================================
        # Portable downstream interface.
        # =====================================================================

        lookup = build_lookup_table(
            period,
            pion_summary,
        )


        portable_replicas = build_portable_replica_table(
            period,
            replicas,
        )


        lookup.to_csv(
            OUTDIR
            / f"kaon_contamination_lookup_{period.lower()}.csv",

            index=False,
        )


        portable_replicas.to_csv(
            OUTDIR
            / f"kaon_contamination_replicas_{period.lower()}.csv",

            index=False,
        )


        write_replica_matrix(
            period,
            portable_replicas,
        )


        plot_kaon_contamination(
            period,
            lookup,
        )


        combined_lookup_tables.append(
            lookup
        )


        combined_replica_tables.append(
            portable_replicas
        )


        print(
            f"\n[{period}] completed.",

            flush=True,
        )


    # =========================================================================
    # Combined files intended for the kinematic-distribution code.
    # =========================================================================

    combined_lookup = pd.concat(
        combined_lookup_tables,
        ignore_index=True,
    )


    combined_replicas = pd.concat(
        combined_replica_tables,
        ignore_index=True,
    )


    combined_lookup.to_csv(
        OUTDIR
        / "kaon_contamination_lookup_all_periods.csv",

        index=False,
    )


    combined_replicas.to_csv(
        OUTDIR
        / "kaon_contamination_replicas_all_periods.csv",

        index=False,
    )


    print(
        "\n"
        "============================================================\n"
        " PID / particle-misidentification study complete\n"
        "============================================================\n"
        f"Output directory:\n"
        f"  {OUTDIR}\n\n"
        "Primary downstream files:\n"
        "  kaon_contamination_lookup_all_periods.csv\n"
        "  kaon_contamination_replicas_all_periods.csv\n\n"
        "These are ready to be folded through the final 24-bin pion "
        "momentum distributions.",

        flush=True,
    )


if __name__ == "__main__":

    main()