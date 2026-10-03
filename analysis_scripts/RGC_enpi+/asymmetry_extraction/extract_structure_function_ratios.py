#!/usr/bin/env python3
"""
extract_structure_function_ratios.py

Standalone event-level asymmetry extraction for the RGC exclusive
e p -> e' n pi+ analysis.

Run this program from

    RGC_enpi+/asymmetry_extraction/

with sibling analysis directories

    ../channel_selection/
    ../dilution_factor/
    ../momentum_corrections/

Nominal statistical model
-------------------------
The three RGC run periods are fitted simultaneously in every (xB, -tprime)
kinematic bin.  The physics parameters are common to Su22, Fa22, and Sp23,
while the likelihood uses the period-dependent beam polarization, the signed
run-by-run target polarization, the run/helicity accumulated charges, and the
period/bin-dependent dilution factors.

The nominal estimator is an unbinned conditional likelihood.  For an event
with accepted kinematics x in period p, the probability of its observed
(run, beam-helicity) state is

    P(r,h | x,p) =
        Q(r,h) g(x; r,h,p) /
        sum_{r' in p, h'=+/-1} Q(r',h') g(x; r',h',p).

The sum in the denominator is restricted to the same run period.  This permits
the three periods to constrain one common set of structure-function ratios
without assuming that the detector acceptance is identical between periods.
Within one period, helicity- and target-state-independent acceptance is assumed.

The fitted structure-function ratios are

    u1   = F_UU^{cos(phi)}     / F_UU
    u2   = F_UU^{cos(2phi)}    / F_UU
    lu1  = F_LU^{sin(phi)}     / F_UU
    ul1  = F_UL^{sin(phi)}     / F_UU
    ul2  = F_UL^{sin(2phi)}    / F_UU
    ll0  = F_LL                / F_UU
    ll1  = F_LL^{cos(phi)}     / F_UU

The unpolarized cosine modulations u1 and u2 float in the nominal fit.

Target-axis treatments
----------------------
Three fits are performed in every kinematic bin.

  nominal
    P_parallel = P_t in the laboratory (beam-axis) frame.  No cos(theta_gamma)
    projection is applied.  These are the reported observables.

  photon_axis_projection
    P_L = P_t cos(theta_gamma), P_T = 0.  This is retained only as an
    interpretation study of the conversion to the virtual-photon axis.

  external_data_informed
    P_L = P_t cos(theta_gamma), P_T = P_t sin(theta_gamma).  Fixed
    transverse amplitudes are supplied from external HERMES measurements.
    After phi_S = 0, the two distinct HERMES sin(phi_h +/- phi_S) UT
    amplitudes become the same observed sin(phi_h) harmonic but retain
    different depolarization factors.  For this controlled leakage study
    they are assigned the same signed magnitude:

      F_UT^{sin(phi-phi_S)} / F_UU   = -0.117
      F_UT^{sin(phi+phi_S)} / F_UU   = -0.117
      F_UT^{sin(2phi-phi_S)} / F_UU  = +0.109
      F_LT^{cos(phi_S)} / F_UU       = -0.150
      F_LT^{cos(phi-phi_S)} / F_UU   = -0.150

    The sin(phi-phi_S) and sin(2phi-phi_S) UT inputs are based on the
    inverse-total-variance weighted HERMES exclusive-pi+ averages over
    0.07 < x_B < 0.35, omitting the lowest-x_B bin.  The additional
    sin(phi+phi_S) input is assigned the same signed value as the measured
    sin(phi-phi_S) input for this study, as requested.  Statistical and
    systematic uncertainties are combined in quadrature for the weighted
    averages. No factor of two is applied.

    The two LT inputs are rounded estimates from the open highest-z HERMES
    SIDIS pi+ points in 0.84 < z < 1.20. Those points are closest to the
    exclusive z -> 1 limit and are not included in the ordinary x and P_hT
    projections. They are documented as digitized/rounded estimates rather
    than exact tabulated values. No factor of two is applied.

After phi_S = 0, the two sin(phi_h +/- phi_S) terms are added coherently
with their proper event-by-event depolarization factors: the minus term has
unit coefficient in this normalization and the plus term carries r_b = B/A.
There is no fixed transverse sin(3phi) or cos(2phi) input. The external-data-informed fit is a
controlled leakage study, not a claim that the SIDIS amplitudes equal the
exclusive amplitudes at CLAS12 kinematics.

The photon-axis-projection and external-data-informed fits are retained as
interpretation studies.  Their shifts relative to the laboratory-frame nominal
result are written to the diagnostic tables and plots, but no target-axis
systematic is assigned and these shifts do not enter the point-to-point
systematic uncertainty.


Dilution-factor uncertainty
---------------------------
The production nominal dilution factor is Method 1, the direct five-target
calculation.  Its bootstrap statistical uncertainty is propagated in the
likelihood as a period/bin-dependent nuisance parameter.  The 4% thermal-
contraction uncertainty is a correlated scale uncertainty and is not folded
into the point-to-point systematic bands here.

For the radiation variation, the extraction uses the dilution factors
recalculated from the matched internal+external-ISR samples and their matched
exclusivity cuts.  Thus the radiation comparison consistently propagates the
change in event kinematics, event selection, and dilution factor through the
full asymmetry extraction.

Momentum-correction diagnostic
------------------------------
The nominal production files include the established electron and pi+ momentum
corrections.  A matched extraction is also performed with the corresponding
uncorrected ROOT files while retaining the same nominal channel-selection cuts,
run information, and nominal dilution factors.  This comparison is retained as
a correction-validation diagnostic only.  Turning the established correction
off is not an equally plausible analysis alternative, so

    abs(a_without_corrections - a_with_corrections)

is NOT assigned as a point-to-point systematic uncertainty and is NOT included
in the final systematic quadrature.  A residual momentum-correction uncertainty
would require a separate estimate of the uncertainty on the correction itself.
The signed comparison vector, covariance diagnostic, and Barlow readout are
still retained.

Channel-selection uncertainty
-----------------------------
The ordinary data are also refitted with the matched tight, nominal, and loose
Mx2 windows (mu +/- 1 sigma, 2 sigma, and 3 sigma, respectively). Each
extraction uses the corresponding dilution factor determined with the same
window. The loose cache
is built once from ROOT and the nominal and tight caches are derived from that
superset. For each observable, the recommended pointwise uncertainty is the RMS
of the tight-minus-nominal and loose-minus-nominal shifts. The complete signed
variation vectors and their outer-product covariance are retained for coherent
propagation across bins.

The radiation sample uses its separately recalculated Method-1 dilution factor
and bootstrap statistical uncertainty.  The correlated 4% dilution-factor
scale uncertainty is intentionally not included in the radiation difference;
it is imposed separately as a scale uncertainty on the final observables.

Bin-migration uncertainty
-------------------------
The established 24x24 MC migration matrix from the previous analysis-note
version is embedded below.  Matrix row i gives the fractional composition of
reconstructed bin i in terms of generated bins j.  Because the published table
is rounded to two decimal places, each row is renormalized to unit sum before
application; otherwise a constant observable would acquire an artificial 1--2%
shift from rounding alone.  For each of the five published polarized ratios,

    a_migrated[i] = sum_j M_normalized[i,j] * a_nominal[j]

    delta_migration[i] = 1.5 * abs(a_migrated[i] - a_nominal[i]).

The factor 1.5 retains the established conservative prescription for GEMC
underestimating the relevant detector resolution.  The two fitted unpolarized
modulations are nuisance quantities required by the simultaneous likelihood;
they are not publication observables and are excluded from migration and from
publication-level systematic averages.

The final published point-to-point systematic is therefore

    sqrt(delta_radiation^2
         + delta_channel_selection^2 + delta_migration^2).

Barlow consistency criterion
----------------------------
Every explicitly evaluated systematic variation is tested with the Barlow
criterion.  For a variation value a_v and nominal value a_0 with statistical
uncertainties sigma_v and sigma_0, this program defines

    sigma_B = sqrt(abs(sigma_v^2 - sigma_0^2)),
    B       = abs(a_v - a_0) / sigma_B.

A variation passes the criterion when B >= 1.  A zero denominator with a
nonzero shift is treated as B = infinity and therefore passes; a zero
denominator with a zero shift is reported as undefined.  The complete Barlow
values and pass/fail/undefined states are written for target-axis,
ISR, momentum-correction, and channel-selection variations.  The Barlow decision is diagnostic and
does not silently delete a requested systematic; both the raw variation and
the criterion result are retained.

Beam- and target-polarization uncertainties are intentionally ignored by this
program.  Their statistical and systematic components will be combined later
into correlated polarization scale systematics.

Parallelism
-----------
ROOT input is read once in a serial preprocessing pass and written to a compact
cache.  Independent kinematic-bin fits are then distributed across at most
eight worker processes.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
import time
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
import re
from pathlib import Path
import sys
from typing import Any, Iterable, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot

try:
    from iminuit import Minuit
except ImportError as exc:
    raise RuntimeError(
        "This program requires iminuit. Install it with: python -m pip install iminuit"
    ) from exc
# endtry


# =============================================================================
# Fixed analysis definitions
# =============================================================================

DIS_W_MIN_GEV = 2.0

PERIODS: tuple[str, ...] = ("su22", "fa22", "sp23")
PERIOD_LABELS: dict[str, str] = {
    "su22": "Su22",
    "fa22": "Fa22",
    "sp23": "Sp23",
}

PERIOD_COLORS: dict[str, str] = {
    "su22": "tab:orange",
    "fa22": "tab:blue",
    "sp23": "tab:green",
}
COMBINED_COLOR = "black"

VARIANT_LABELS: dict[str, str] = {
    "nominal": r"Nominal lab frame: $P_{\parallel}=P_t$",
    "photon_axis_projection": r"Photon-axis projection: $P_L=P_t\cos\theta_\gamma$",
    "external_data_informed": r"External-data-informed transverse leakage",
}
VARIANT_COLORS: dict[str, str] = {
    "nominal": "tab:blue",
    "photon_axis_projection": "tab:purple",
    "external_data_informed": "tab:orange",
}

# External transverse amplitudes used in the leakage fit.
#
# Exclusive UT reference:
#   A. Airapetian et al. (HERMES),
#   "Single-spin azimuthal asymmetry in exclusive electroproduction of pi+
#   mesons on transversely polarized protons,"
#   Phys. Lett. B 682 (2010) 345-350,
#   DOI: 10.1016/j.physletb.2009.11.039,
#   arXiv:0907.2596.
#   Official HERMES ASCII data file: ivana.AUTexclpi-1.dat.
#
# SIDIS LT reference:
#   A. Airapetian et al. (HERMES),
#   "Azimuthal single- and double-spin asymmetries in semi-inclusive
#   deep-inelastic lepton scattering by transversely polarized protons,"
#   JHEP 12 (2020) 010,
#   DOI: 10.1007/JHEP12(2020)010,
#   arXiv:2007.07755.
#
# The UT values below are direct inverse-total-variance weighted averages of
# exact tabulated exclusive-pi+ measurements. The LT values are rounded
# estimates from the open highest-z pi+ points (0.84 < z < 1.20), which are not
# included in the ordinary x and P_hT projections. No factor of two is applied
# to any of the four inputs.
EXTERNAL_TRANSVERSE_INPUTS: dict[str, dict[str, Any]] = {
    "ut_sin_phi": {
        "value": -0.117,
        "measured_central_value": -0.11729214289950972,
        "scale_factor": 1.0,
        "observable": "A_UT,l^{sin(phi-phi_S)}",
        "channel": "exclusive pi+",
        "kinematic_selection": {
            "averaging_policy": (
                "Inverse-total-variance weighted average over the HERMES "
                "x_B bins 0.07-0.10, 0.10-0.15, and 0.15-0.35; the "
                "0.03-0.07 bin is omitted. Statistical and systematic "
                "uncertainties are combined in quadrature."
            ),
            "xB_range_used": [0.07, 0.35],
            "weighted_mean_xB": 0.13848440446760085,
            "weighted_mean_uncertainty": 0.09362561740330531,
            "input_rows": [
                {
                    "xB_bin": [0.07, 0.10],
                    "mean_xB": 0.09,
                    "asymmetry": -0.071,
                    "stat_uncertainty": 0.180,
                    "syst_uncertainty": 0.022,
                },
                {
                    "xB_bin": [0.10, 0.15],
                    "mean_xB": 0.13,
                    "asymmetry": -0.093,
                    "stat_uncertainty": 0.130,
                    "syst_uncertainty": 0.029,
                },
                {
                    "xB_bin": [0.15, 0.35],
                    "mean_xB": 0.21,
                    "asymmetry": -0.219,
                    "stat_uncertainty": 0.180,
                    "syst_uncertainty": 0.065,
                },
            ],
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": (
                "A. Airapetian et al., Phys. Lett. B 682 (2010) 345-350"
            ),
            "title": (
                "Single-spin azimuthal asymmetry in exclusive "
                "electroproduction of pi+ mesons on transversely "
                "polarized protons"
            ),
            "doi": "10.1016/j.physletb.2009.11.039",
            "arxiv": "0907.2596",
            "data_source": "Official HERMES ASCII file ivana.AUTexclpi-1.dat",
            "data_note": (
                "Direct weighted average of exact tabulated exclusive-pi+ "
                "measurements; no doubling factor."
            ),
        },
    },
    "ut_sin_phi_plus": {
        "value": -0.117,
        "measured_central_value": -0.117,
        "scale_factor": 1.0,
        "observable": "A_UT,l^{sin(phi+phi_S)}",
        "channel": "exclusive pi+",
        "kinematic_selection": {
            "assignment_policy": (
                "Assigned the same signed magnitude as the HERMES-weighted "
                "A_UT,l^{sin(phi-phi_S)} input for the controlled transverse-"
                "leakage study.  The two terms remain distinct in the cross "
                "section because sin(phi+phi_S) carries the B/A = epsilon "
                "depolarization factor after phi_S = 0."
            ),
            "xB_range_used": [0.07, 0.35],
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": (
                "A. Airapetian et al., Phys. Lett. B 682 (2010) 345-350"
            ),
            "title": (
                "Single-spin azimuthal asymmetry in exclusive "
                "electroproduction of pi+ mesons on transversely "
                "polarized protons"
            ),
            "doi": "10.1016/j.physletb.2009.11.039",
            "arxiv": "0907.2596",
            "data_note": (
                "For this analysis study the signed value is set equal to "
                "the sin(phi-phi_S) input, -0.117; it is not presented as a "
                "new independent weighted average."
            ),
        },
    },
    "ut_sin_2phi": {
        "value": 0.109,
        "measured_central_value": 0.10930777687298002,
        "scale_factor": 1.0,
        "observable": "A_UT,l^{sin(2phi-phi_S)}",
        "channel": "exclusive pi+",
        "kinematic_selection": {
            "averaging_policy": (
                "Inverse-total-variance weighted average over the HERMES "
                "x_B bins 0.07-0.10, 0.10-0.15, and 0.15-0.35; the "
                "0.03-0.07 bin is omitted. Statistical and systematic "
                "uncertainties are combined in quadrature."
            ),
            "xB_range_used": [0.07, 0.35],
            "weighted_mean_xB": 0.13896698591627263,
            "weighted_mean_uncertainty": 0.09952058355881978,
            "input_rows": [
                {
                    "xB_bin": [0.07, 0.10],
                    "mean_xB": 0.09,
                    "asymmetry": -0.196,
                    "stat_uncertainty": 0.181,
                    "syst_uncertainty": 0.106,
                },
                {
                    "xB_bin": [0.10, 0.15],
                    "mean_xB": 0.13,
                    "asymmetry": 0.178,
                    "stat_uncertainty": 0.132,
                    "syst_uncertainty": 0.024,
                },
                {
                    "xB_bin": [0.15, 0.35],
                    "mean_xB": 0.21,
                    "asymmetry": 0.247,
                    "stat_uncertainty": 0.192,
                    "syst_uncertainty": 0.085,
                },
            ],
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": (
                "A. Airapetian et al., Phys. Lett. B 682 (2010) 345-350"
            ),
            "title": (
                "Single-spin azimuthal asymmetry in exclusive "
                "electroproduction of pi+ mesons on transversely "
                "polarized protons"
            ),
            "doi": "10.1016/j.physletb.2009.11.039",
            "arxiv": "0907.2596",
            "data_source": "Official HERMES ASCII file ivana.AUTexclpi-1.dat",
            "data_note": (
                "Direct weighted average of exact tabulated exclusive-pi+ "
                "measurements; no doubling factor."
            ),
        },
    },
    "lt_cos_0phi": {
        "value": -0.150,
        "measured_central_value": -0.150,
        "scale_factor": 1.0,
        "observable": "2<cos(phi_S)>/sqrt(2 epsilon (1-epsilon))",
        "channel": "SIDIS pi+",
        "kinematic_selection": {
            "z_bin": [0.84, 1.20],
            "x_acceptance": [0.023, 0.600],
            "selection_note": (
                "Open highest-z point in the one-dimensional pi+ z "
                "projection. This point is not included in the ordinary "
                "x and P_hT projections and is used as the closest "
                "available SIDIS proxy for z approximately 1."
            ),
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": "A. Airapetian et al., JHEP 12 (2020) 010",
            "title": (
                "Azimuthal single- and double-spin asymmetries in "
                "semi-inclusive deep-inelastic lepton scattering by "
                "transversely polarized protons"
            ),
            "doi": "10.1007/JHEP12(2020)010",
            "arxiv": "2007.07755",
            "figure": 31,
            "data_note": (
                "Digitized/rounded estimate from the open highest-z pi+ "
                "point; not an exact tabulated value; no doubling factor."
            ),
        },
    },
    "lt_cos_phi": {
        "value": -0.150,
        "measured_central_value": -0.150,
        "scale_factor": 1.0,
        "observable": "2<cos(phi-phi_S)>/sqrt(1-epsilon^2)",
        "channel": "SIDIS pi+",
        "kinematic_selection": {
            "z_bin": [0.84, 1.20],
            "x_acceptance": [0.023, 0.600],
            "selection_note": (
                "Open highest-z point in the one-dimensional pi+ z "
                "projection. This point is not included in the ordinary "
                "x and P_hT projections and is used as the closest "
                "available SIDIS proxy for z approximately 1."
            ),
        },
        "reference": {
            "collaboration": "HERMES",
            "citation": "A. Airapetian et al., JHEP 12 (2020) 010",
            "title": (
                "Azimuthal single- and double-spin asymmetries in "
                "semi-inclusive deep-inelastic lepton scattering by "
                "transversely polarized protons"
            ),
            "doi": "10.1007/JHEP12(2020)010",
            "arxiv": "2007.07755",
            "figure": 21,
            "data_note": (
                "Digitized/rounded estimate from the open highest-z pi+ "
                "point; not an exact tabulated value; no doubling factor."
            ),
        },
    },
}

PERIOD_INDEX: dict[str, int] = {
    period: index for index, period in enumerate(PERIODS)
}

BEAM_POLARIZATION: dict[str, float] = {
    "su22": 0.8384,
    "fa22": 0.8372,
    "sp23": 0.8040,
}

PROTON_MASS_GEV = 0.9382720813

BEAM_ENERGY_GEV: dict[str, float] = {
    "su22": 10.5473,
    "fa22": 10.5563,
    "sp23": 10.5593,
}

XB_BINS: tuple[tuple[float, float], ...] = (
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
)

MINUS_TPRIME_BINS_GEV2: tuple[tuple[float, float], ...] = (
    (0.05, 0.25),
    (0.25, 0.45),
    (0.45, 0.65),
    (0.65, 0.85),
    (0.85, 1.05),
    (1.05, 1.25),
)

NUMBER_OF_BINS = len(XB_BINS) * len(MINUS_TPRIME_BINS_GEV2)
MAXIMUM_WORKERS = 8
DEFAULT_TREE_NAME = "PhysicsEvents"
DEFAULT_CHUNK_SIZE = "250 MB"
DEFAULT_OUTPUT_DIR = Path("output/asymmetry_extraction")
DEFAULT_CACHE_PATH = DEFAULT_OUTPUT_DIR / "cache/selected_events.npz"

DEFAULT_RUN_INFO_CSV = Path("clas12_run_info.csv")
DEFAULT_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts.json"
)
DEFAULT_POLYNOMIAL_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "nominal/final_polynomial_only_cuts/tables/"
    "final_polynomial_only_mx2_cuts.json"
)

DEFAULT_DILUTION_DIR = Path(
    "../dilution_factor/output/dilution_factor_determination"
)

PAPER_VERSIONS_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
    "paper_versions"
)

DEFAULT_INPUTS: dict[str, Path] = {
    period: PAPER_VERSIONS_DIR
    / f"rgc_{period}_inb_NH3_epi+_mom_corrections.root"
    for period in PERIODS
}

DEFAULT_ISR_INPUTS: dict[str, Path] = {
    period: PAPER_VERSIONS_DIR
    / f"rgc_{period}_inb_NH3_epi+_ISR_externalISR_mom_corrections.root"
    for period in PERIODS
}

UNCORRECTED_INPUTS_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+"
)

DEFAULT_UNCORRECTED_INPUTS: dict[str, Path] = {
    period: UNCORRECTED_INPUTS_DIR
    / f"rgc_{period}_inb_NH3_epi+_2.root"
    for period in PERIODS
}

DEFAULT_CHANNEL_SELECTION_MANIFEST = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "channel_selection_manifest.json"
)
DEFAULT_ISR_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "isr/final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts.json"
)
FIT_VARIANTS: tuple[str, ...] = (
    "nominal",
    "photon_axis_projection",
    "external_data_informed",
)


PHYSICS_PARAMETERS: tuple[str, ...] = (
    "u1",
    "u2",
    "lu1",
    "ul1",
    "ul2",
    "ll0",
    "ll1",
)

PARAMETER_INITIAL_VALUES: dict[str, float] = {
    "u1": 0.0,
    "u2": 0.0,
    "lu1": 0.05,
    "ul1": 0.0,
    "ul2": 0.0,
    "ll0": 0.0,
    "ll1": 0.0,
}

PARAMETER_LIMITS: dict[str, tuple[float, float]] = {
    name: (-1.5, 1.5) for name in PHYSICS_PARAMETERS
}


PARAMETER_LABELS: dict[str, str] = {
    "u1": r"$F_{UU}^{\cos\phi}/F_{UU}$",
    "u2": r"$F_{UU}^{\cos2\phi}/F_{UU}$",
    "lu1": r"$F_{LU}^{\sin\phi}/F_{UU}$",
    "ul1": r"$A_{UL,\mathrm{lab}}^{\sin\phi}$",
    "ul2": r"$A_{UL,\mathrm{lab}}^{\sin2\phi}$",
    "ll0": r"$A_{LL,\mathrm{lab}}$",
    "ll1": r"$A_{LL,\mathrm{lab}}^{\cos\phi}$",
}

PARAMETER_Y_LIMITS: dict[str, tuple[float, float] | None] = {
    "u1": (-1.0, 1.0),
    "u2": (-1.0, 1.0),
    "lu1": (-0.5, 0.5),
    "ul1": (-0.5, 0.5),
    "ul2": (-0.5, 0.5),
    "ll0": (-1.0, 1.0),
    "ll1": (-1.0, 1.0),
}

# Standardized systematic-comparison axes.  Single-spin ratios share one
# scale, while unpolarized and double-spin ratios share the broader scale.
SINGLE_SPIN_PARAMETERS = frozenset(("lu1", "ul1", "ul2"))
UNPOLARIZED_DOUBLE_SPIN_PARAMETERS = frozenset(
    ("u1", "u2", "ll0", "ll1")
)
SYSTEMATIC_COMPARISON_SINGLE_SPIN_Y_LIMITS = (-0.5, 0.5)
SYSTEMATIC_COMPARISON_UNPOLARIZED_DOUBLE_Y_LIMITS = (-1.0, 1.0)
SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS = (1.0e-2, 1.0e1)
SYSTEMATIC_COMPARISON_FIGSIZE = (11, 10)
SYSTEMATIC_COMPARISON_DPI = 200
PUBLISHED_SYSTEMATIC_PARAMETERS: tuple[str, ...] = (
    "lu1", "ul1", "ul2", "ll0", "ll1",
)

# Established MC bin-migration composition matrix from the previous analysis
# note.  Row i is the generated-bin composition of reconstructed bin i.
# Values were published rounded to two decimal places, so rows are normalized
# before use.
MIGRATION_RESOLUTION_SCALE = 1.5
MIGRATION_MATRIX = np.asarray([
    [0.93,0.03,0,0,0,0,0.03,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0.11,0.80,0.04,0,0,0,0,0.04,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0.01,0.11,0.83,0.02,0,0,0,0,0.01,0.03,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0.10,0.82,0.03,0,0,0,0.01,0.03,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0.01,0.09,0.83,0.04,0,0,0,0,0.01,0.02,0,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0.01,0.11,0.84,0,0,0,0,0,0.01,0.02,0,0,0,0,0,0,0,0,0,0,0],
    [0.01,0,0,0,0,0,0.93,0.02,0,0,0,0,0.05,0,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0.16,0.70,0.02,0,0,0,0.04,0.07,0,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0.01,0.14,0.71,0.02,0,0,0,0.04,0.07,0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0.01,0.17,0.68,0.03,0,0,0.01,0.03,0.06,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0.01,0.14,0.70,0.05,0,0,0.01,0.03,0.04,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0,0.01,0.14,0.76,0,0,0.01,0.01,0.03,0.03,0,0,0,0,0,0],
    [0,0,0,0,0,0,0.01,0,0,0,0,0,0.94,0.03,0,0,0,0,0.02,0,0,0,0,0],
    [0,0,0,0,0,0,0,0.01,0,0,0,0,0.17,0.74,0.03,0,0,0,0.02,0.02,0,0,0,0],
    [0,0,0,0,0,0,0,0,0.01,0,0,0,0.01,0.16,0.72,0.03,0,0,0.01,0.02,0.03,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.19,0.70,0.04,0,0,0.01,0.02,0.02,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.17,0.69,0.04,0,0.01,0.01,0.01,0.03,0],
    [0,0,0,0,0,0,0,0,0,0,0,0.01,0,0,0,0.02,0.19,0.71,0,0,0.01,0.01,0.02,0.03],
    [0,0,0,0,0,0,0,0,0,0,0,0,0.03,0,0,0,0,0,0.89,0.06,0.01,0,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0.01,0,0,0,0.23,0.69,0.04,0.01,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0.01,0,0,0.02,0.24,0.67,0.05,0,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0,0,0,0.03,0.24,0.66,0.05,0],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.01,0,0,0.01,0.05,0.18,0.70,0.04],
    [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.02,0,0,0.01,0.06,0.25,0.65],
], dtype=np.float64)

# Physics-motivated 3x4 canvas layout:
#   top row:    UU (unpolarized) structure-function ratios
#   middle row: all fitted single-spin ratios (LU and UL)
#   bottom row: all fitted double-spin ratios (LL)
#
# The unused cells in the top and bottom rows are intentionally left blank.
PHYSICS_PANEL_ROWS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("UU", ("u1", "u2")),
    ("Single-spin", ("lu1", "ul1", "ul2")),
    ("Double-spin", ("ll0", "ll1")),
)

AGGREGATED_PANEL_ORDER: tuple[str, ...] = tuple(
    parameter
    for _, row_parameters in PHYSICS_PANEL_ROWS
    for parameter in row_parameters
)

PROBABILITY_FLOOR = 1.0e-300
CROSS_SECTION_FLOOR = 1.0e-10
INVALID_NLL = 1.0e100

BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "runnum": ("runnum", "run", "run_number", "RunNumber"),
    "helicity": ("helicity", "hel", "beam_helicity"),
    "xB": ("xB", "x", "xb", "x_b"),
    "tprime": ("tprime", "t_prime", "tp", "tPrime"),
    "t": ("t", "T", "mandelstam_t"),
    "W": ("W", "w", "invariant_mass_W"),
    "Mx2": (
        "Mx2",
        "mx2",
        "Mx2_epi",
        "Mx2_epip",
        "missing_mass_squared",
        "missing_mass2",
    ),
    "phi": ("phi", "phi1", "phi_h", "trento_phi"),
    "Q2": ("Q2", "q2"),
    "y": ("y", "inelasticity"),
    "e_p": ("e_p", "p_e", "electron_p"),
    "e_theta": ("e_theta", "theta_e", "electron_theta"),
    "DepA": ("DepA", "depA"),
    "DepB": ("DepB", "depB"),
    "DepC": ("DepC", "depC"),
    "DepV": ("DepV", "depV"),
    "DepW": ("DepW", "depW"),
}

PERIOD_DIAGNOSTIC_BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "e_phi": ("e_phi", "phi_e", "electron_phi"),
    "p_phi": ("p_phi", "pip_phi", "pi_phi", "hadron_phi"),
    "p_p": ("p_p", "pip_p", "pi_p", "hadron_p"),
    "p_theta": ("p_theta", "pip_theta", "pi_theta", "hadron_theta"),
}


# =============================================================================
# Data containers
# =============================================================================

@dataclass(frozen=True)
class RunRecord:
    period: str
    run: int
    charge_plus: float
    charge_minus: float
    target_polarization: float
    target_polarization_uncertainty: float


@dataclass(frozen=True)
class DilutionRecord:
    period: str
    bin_number: int
    x_index: int
    t_index: int
    value: float
    stat_uncertainty: float


@dataclass(frozen=True)
class CutRecord:
    period: str
    bin_number: int
    x_index: int
    t_index: int
    low_gev2: float
    high_gev2: float


# =============================================================================
# General utilities
# =============================================================================

def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    # endif
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    # endif
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    # endif
    if isinstance(value, np.integer):
        return int(value)
    # endif
    if isinstance(value, np.floating):
        value = float(value)
    # endif
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    # endif
    if isinstance(value, Path):
        return str(value)
    # endif
    return value


def write_json(path: Path, payload: Any) -> None:
    ensure_directory(path.parent)
    path.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


BARLOW_PASS_THRESHOLD = 1.0


def calculate_barlow_arrays(
    nominal_values: pd.Series | np.ndarray,
    variation_values: pd.Series | np.ndarray,
    nominal_stat: pd.Series | np.ndarray,
    variation_stat: pd.Series | np.ndarray,
) -> dict[str, np.ndarray]:
    """Return binwise Barlow diagnostics for one systematic variation.

    The denominator follows the standard correlated-sample prescription

        sqrt(abs(sigma_variation**2 - sigma_nominal**2)).

    Status codes are ``pass``, ``fail``, or ``undefined``.  Non-finite inputs
    are undefined.  A finite nonzero shift with exactly zero denominator has
    infinite Barlow significance and passes.
    """
    nominal_values_array = np.asarray(nominal_values, dtype=np.float64)
    variation_values_array = np.asarray(variation_values, dtype=np.float64)
    nominal_stat_array = np.asarray(nominal_stat, dtype=np.float64)
    variation_stat_array = np.asarray(variation_stat, dtype=np.float64)

    difference = variation_values_array - nominal_values_array
    variance_difference = variation_stat_array**2 - nominal_stat_array**2
    denominator = np.sqrt(np.abs(variance_difference))
    barlow = np.full(difference.shape, np.nan, dtype=np.float64)
    status = np.full(difference.shape, "undefined", dtype=object)

    finite_inputs = (
        np.isfinite(difference)
        & np.isfinite(nominal_stat_array)
        & np.isfinite(variation_stat_array)
        & (nominal_stat_array >= 0.0)
        & (variation_stat_array >= 0.0)
    )
    positive_denominator = finite_inputs & (denominator > 0.0)
    barlow[positive_denominator] = (
        np.abs(difference[positive_denominator])
        / denominator[positive_denominator]
    )

    zero_denominator_nonzero_shift = (
        finite_inputs & (denominator == 0.0) & (np.abs(difference) > 0.0)
    )
    barlow[zero_denominator_nonzero_shift] = np.inf

    defined = np.isfinite(barlow) | np.isinf(barlow)
    passed = defined & (barlow >= BARLOW_PASS_THRESHOLD)
    failed = defined & ~passed
    status[passed] = "pass"
    status[failed] = "fail"

    return {
        "difference": difference,
        "variance_difference": variance_difference,
        "denominator": denominator,
        "barlow": barlow,
        "passed": passed,
        "status": status,
    }


def summarize_barlow_status(
    frame: pd.DataFrame,
    *,
    effect: str,
    variation: str,
    parameter: str,
    barlow_column: str,
    status_column: str,
) -> dict[str, Any]:
    """Summarize one parameter/variation Barlow test."""
    status = frame[status_column].astype(str)
    values = pd.to_numeric(frame[barlow_column], errors="coerce").to_numpy()
    finite_values = values[np.isfinite(values)]
    return {
        "effect": effect,
        "variation": variation,
        "parameter": parameter,
        "threshold": BARLOW_PASS_THRESHOLD,
        "number_of_bins": int(len(frame)),
        "number_pass": int((status == "pass").sum()),
        "number_fail": int((status == "fail").sum()),
        "number_undefined": int((status == "undefined").sum()),
        "number_defined": int((status != "undefined").sum()),
        "pass_fraction_defined": (
            float((status == "pass").sum() / (status != "undefined").sum())
            if int((status != "undefined").sum()) > 0 else None
        ),
        "fail_fraction_defined": (
            float((status == "fail").sum() / (status != "undefined").sum())
            if int((status != "undefined").sum()) > 0 else None
        ),
        "all_defined_bins_pass": bool(
            (status == "pass").all()
        ),
        "minimum_finite_barlow": (
            float(np.min(finite_values)) if finite_values.size else None
        ),
        "maximum_finite_barlow": (
            float(np.max(finite_values)) if finite_values.size else None
        ),
    }


def write_barlow_summary_products(
    records: list[dict[str, Any]],
    output_dir: Path,
) -> dict[str, Any]:
    """Write parameter-level and effect-level Barlow summary products."""
    ensure_directory(output_dir)
    frame = pd.DataFrame(records)
    csv_path = output_dir / "barlow_criteria_summary.csv"
    effect_csv_path = output_dir / "barlow_criteria_effect_summary.csv"
    json_path = output_dir / "barlow_criteria_summary.json"
    text_path = output_dir / "barlow_criteria_summary.txt"
    frame.to_csv(csv_path, index=False)

    effect_records: list[dict[str, Any]] = []
    if not frame.empty:
        for (effect, variation), group in frame.groupby(
            ["effect", "variation"], sort=False
        ):
            number_pass = int(group["number_pass"].sum())
            number_fail = int(group["number_fail"].sum())
            number_undefined = int(group["number_undefined"].sum())
            finite_minima = pd.to_numeric(
                group["minimum_finite_barlow"], errors="coerce"
            ).dropna()
            finite_maxima = pd.to_numeric(
                group["maximum_finite_barlow"], errors="coerce"
            ).dropna()
            effect_records.append({
                "effect": str(effect),
                "variation": str(variation),
                "threshold": BARLOW_PASS_THRESHOLD,
                "number_of_parameters": int(len(group)),
                "number_of_bin_tests": int(group["number_of_bins"].sum()),
                "number_pass": number_pass,
                "number_fail": number_fail,
                "number_undefined": number_undefined,
                "number_defined": number_pass + number_fail,
                "pass_fraction_defined": (
                    float(number_pass / (number_pass + number_fail))
                    if (number_pass + number_fail) > 0 else None
                ),
                "fail_fraction_defined": (
                    float(number_fail / (number_pass + number_fail))
                    if (number_pass + number_fail) > 0 else None
                ),
                "all_defined_bins_pass": bool(number_fail == 0),
                "minimum_finite_barlow": (
                    float(finite_minima.min()) if not finite_minima.empty else None
                ),
                "maximum_finite_barlow": (
                    float(finite_maxima.max()) if not finite_maxima.empty else None
                ),
            })
        # endfor
    # endif
    effect_frame = pd.DataFrame(effect_records)
    effect_frame.to_csv(effect_csv_path, index=False)

    write_json(json_path, {
        "schema_version": 3,
        "definition": (
            "B = abs(variation - nominal) / "
            "sqrt(abs(sigma_variation^2 - sigma_nominal^2)); pass if B >= 1"
        ),
        "effect_summary": effect_records,
        "parameter_summary": records,
    })

    lines = [
        "Barlow criterion summary",
        "=========================",
        (
            "Definition: B = |variation - nominal| / "
            "sqrt(|sigma_variation^2 - sigma_nominal^2|); pass if B >= 1."
        ),
        "",
        "Effect-level summary",
        "--------------------",
    ]
    for record in effect_records:
        defined = int(record["number_defined"])
        pass_fraction = record["pass_fraction_defined"]
        fail_fraction = record["fail_fraction_defined"]
        pass_percent = 100.0 * pass_fraction if pass_fraction is not None else float("nan")
        fail_percent = 100.0 * fail_fraction if fail_fraction is not None else float("nan")
        lines.append(
            f"{record['effect']:22s} {record['variation']:32s}: "
            f"{record['number_pass']:3d}/{defined:3d} defined bins pass "
            f"({pass_percent:5.1f}%); "
            f"{record['number_fail']:3d}/{defined:3d} fail "
            f"({fail_percent:5.1f}%); "
            f"undefined={record['number_undefined']:3d}"
        )
    # endfor
    lines.extend(["", "Parameter-level summary", "-----------------------"])
    for record in records:
        defined = int(record["number_defined"])
        pass_fraction = record["pass_fraction_defined"]
        fail_fraction = record["fail_fraction_defined"]
        pass_percent = 100.0 * pass_fraction if pass_fraction is not None else float("nan")
        fail_percent = 100.0 * fail_fraction if fail_fraction is not None else float("nan")
        lines.append(
            f"{record['effect']:22s} {record['variation']:32s} "
            f"{record['parameter']:4s}: "
            f"{record['number_pass']:2d}/{defined:2d} pass ({pass_percent:5.1f}%); "
            f"{record['number_fail']:2d}/{defined:2d} fail ({fail_percent:5.1f}%); "
            f"undefined={record['number_undefined']:2d}"
        )
    # endfor
    text_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "csv": str(csv_path),
        "effect_csv": str(effect_csv_path),
        "json": str(json_path),
        "text": str(text_path),
        "effect_summary": effect_records,
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
        # endfor
    # endwith
    return digest.hexdigest()


def combined_bin_number(x_index: int, t_index: int) -> int:
    return x_index * len(MINUS_TPRIME_BINS_GEV2) + t_index + 1


def bin_indices(
    x_values: np.ndarray,
    minus_tprime_values: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_index = np.full(x_values.shape, -1, dtype=np.int16)
    t_index = np.full(minus_tprime_values.shape, -1, dtype=np.int16)

    for index, (low, high) in enumerate(XB_BINS):
        mask = (x_values >= low) & (x_values < high)
        x_index[mask] = index
    # endfor

    for index, (low, high) in enumerate(MINUS_TPRIME_BINS_GEV2):
        mask = (minus_tprime_values >= low) & (minus_tprime_values < high)
        t_index[mask] = index
    # endfor

    valid = (x_index >= 0) & (t_index >= 0)
    number = np.full(x_values.shape, -1, dtype=np.int16)
    number[valid] = (
        x_index[valid] * len(MINUS_TPRIME_BINS_GEV2)
        + t_index[valid]
        + 1
    )
    return x_index, t_index, number


def resolve_branch(
    tree: uproot.behaviors.TTree.TTree,
    aliases: Iterable[str],
) -> str:
    available = {str(key).split(";")[0] for key in tree.keys()}
    for alias in aliases:
        if alias in available:
            return alias
        # endif
    # endfor
    raise RuntimeError(
        f"None of the aliases {tuple(aliases)} was found. "
        f"Available branches include {sorted(available)[:120]}"
    )


def detect_angle_unit(values: np.ndarray, branch_name: str) -> str:
    finite = np.abs(values[np.isfinite(values)])
    if finite.size == 0:
        raise RuntimeError(f"Cannot determine units for empty branch {branch_name}.")
    # endif
    percentile = float(np.percentile(finite, 99.5))
    return "degrees" if percentile > 2.0 * math.pi + 0.25 else "radians"


def angle_to_radians(values: np.ndarray, branch_name: str) -> np.ndarray:
    unit = detect_angle_unit(values, branch_name)
    if unit == "degrees":
        return np.deg2rad(values)
    # endif
    return values


# =============================================================================
# Run information
# =============================================================================

def parse_run_info_csv(path: Path) -> dict[int, RunRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing run-information CSV: {path}")
    # endif

    period_lookup = {
        "Su22": "su22",
        "Fa22": "fa22",
        "Sp23": "sp23",
    }

    records: dict[int, RunRecord] = {}
    current_period: str | None = None
    current_target: str | None = None

    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8").splitlines(),
        start=1,
    ):
        line = raw_line.strip()
        if not line:
            continue
        # endif

        if line.startswith("#"):
            header = line[1:].strip().split()
            current_period = None
            current_target = None
            if len(header) >= 3 and header[0] == "RGC":
                current_period = period_lookup.get(header[1])
                current_target = header[2]
            # endif
            continue
        # endif

        if current_period not in PERIODS or current_target != "NH3":
            continue
        # endif

        fields = [field.strip() for field in line.split(",")]
        if len(fields) < 6:
            raise RuntimeError(
                f"Malformed NH3 run-info row at line {line_number}: {raw_line}"
            )
        # endif

        try:
            run = int(fields[0])
            charge_plus = float(fields[2])
            charge_minus = float(fields[3])
            target_polarization = float(fields[4])
            target_uncertainty = float(fields[-1])
        except ValueError as exc:
            raise RuntimeError(
                f"Invalid numeric field at line {line_number}: {raw_line}"
            ) from exc
        # endtry

        if run in records:
            raise RuntimeError(f"Duplicate NH3 run {run} in {path}.")
        # endif
        if not (
            math.isfinite(charge_plus)
            and math.isfinite(charge_minus)
            and charge_plus >= 0.0
            and charge_minus >= 0.0
        ):
            raise RuntimeError(
                f"Run {run} has invalid helicity charges: "
                f"Q+={charge_plus}, Q-={charge_minus}. Charges must be "
                "finite and nonnegative."
            )
        # endif

        disabled_by_qa = charge_plus == 0.0 and charge_minus == 0.0
        only_one_helicity_disabled = (charge_plus == 0.0) != (charge_minus == 0.0)
        if only_one_helicity_disabled:
            raise RuntimeError(
                f"Run {run} has only one zero helicity charge: "
                f"Q+={charge_plus}, Q-={charge_minus}. The current likelihood "
                "expects either two positive helicity charges or a fully "
                "QA-disabled run with both charges set to zero."
            )
        # endif

        # Fully QA-disabled runs are intentionally retained in the CSV with
        # Q+ = Q- = 0.  They are valid bookkeeping rows but must never enter
        # the likelihood.  A nonzero finite target polarization is required
        # only for active runs.
        if not disabled_by_qa and (
            not math.isfinite(target_polarization)
            or target_polarization == 0.0
        ):
            raise RuntimeError(
                f"Active run {run} has invalid target polarization "
                f"{target_polarization}."
            )
        # endif

        records[run] = RunRecord(
            period=current_period,
            run=run,
            charge_plus=charge_plus,
            charge_minus=charge_minus,
            target_polarization=target_polarization,
            target_polarization_uncertainty=target_uncertainty,
        )
    # endfor

    for period in PERIODS:
        if not any(
            record.period == period
            and record.charge_plus > 0.0
            and record.charge_minus > 0.0
            for record in records.values()
        ):
            raise RuntimeError(
                f"No active NH3 run records with positive Q+ and Q- were "
                f"found for {period}."
            )
        # endif
    # endfor

    return records


def run_state_arrays(
    run_records: Mapping[int, RunRecord],
) -> dict[str, dict[str, np.ndarray]]:
    result: dict[str, dict[str, np.ndarray]] = {}
    for period in PERIODS:
        period_records = sorted(
            (
                record
                for record in run_records.values()
                if (
                    record.period == period
                    and record.charge_plus > 0.0
                    and record.charge_minus > 0.0
                )
            ),
            key=lambda record: record.run,
        )
        result[period] = {
            "run": np.asarray([record.run for record in period_records], dtype=np.int32),
            "pt": np.asarray(
                [record.target_polarization for record in period_records],
                dtype=np.float64,
            ),
            "q_plus": np.asarray(
                [record.charge_plus for record in period_records],
                dtype=np.float64,
            ),
            "q_minus": np.asarray(
                [record.charge_minus for record in period_records],
                dtype=np.float64,
            ),
        }
    # endfor
    return result


# =============================================================================
# Channel-selection cuts
# =============================================================================

def read_json_tolerating_trailing_escaped_whitespace(path: Path) -> Any:
    """Read JSON, tolerating only the known trailing literal ``\\n`` writer artifact."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing JSON file: {path}")
    # endif
    raw_json = path.read_text(encoding="utf-8")
    try:
        return json.loads(raw_json)
    except json.JSONDecodeError as exc:
        cleaned_json = re.sub(r"(?:\\[nrt])+\s*$", "", raw_json)
        if cleaned_json == raw_json:
            raise
        # endif
        try:
            payload = json.loads(cleaned_json)
        except json.JSONDecodeError:
            raise exc
        # endtry
        print(
            f"[json] WARNING: ignored literal escaped trailing whitespace in {path}",
            flush=True,
        )
        return payload
    # endtry


def write_hybrid_fit_method_cut_jsons(
    carbon_path: Path,
    polynomial_path: Path,
    output_dir: Path,
) -> tuple[Path, Path]:
    """Build nominal +/-2sigma hybrid cuts for the fit-method decomposition.

    The two outputs isolate (1) the polynomial centroid at the carbon width and
    (2) the carbon centroid at the polynomial width.  Centers and half-widths
    are inferred directly from each period/bin nominal interval, so this stays
    compatible with the existing channel-selection JSON schema.
    """
    carbon = read_json_tolerating_trailing_escaped_whitespace(carbon_path)
    polynomial = read_json_tolerating_trailing_escaped_whitespace(polynomial_path)
    cperiods = carbon.get("periods")
    pperiods = polynomial.get("periods")
    if not isinstance(cperiods, dict) or not isinstance(pperiods, dict):
        raise RuntimeError("Fit-method cut JSONs must contain a 'periods' mapping.")
    # endif

    import copy
    mu_poly_sigma_c = copy.deepcopy(carbon)
    mu_c_sigma_poly = copy.deepcopy(carbon)
    for period in PERIODS:
        crows = cperiods.get(period)
        prows = pperiods.get(period)
        if not isinstance(crows, list) or not isinstance(prows, list):
            raise RuntimeError(f"Missing period rows for {period} in fit-method JSONs.")
        # endif
        pmap = {int(row["bin_number"]): row for row in prows}
        out_mu = mu_poly_sigma_c["periods"][period]
        out_sigma = mu_c_sigma_poly["periods"][period]
        for index, crow in enumerate(crows):
            bin_number = int(crow["bin_number"])
            prow = pmap.get(bin_number)
            if prow is None:
                raise RuntimeError(f"Polynomial cut JSON missing {period} bin {bin_number}.")
            # endif
            clo, chi = map(float, crow["nominal"])
            plo, phi = map(float, prow["nominal"])
            mu_c = 0.5 * (clo + chi)
            mu_p = 0.5 * (plo + phi)
            half_c = 0.5 * (chi - clo)
            half_p = 0.5 * (phi - plo)
            out_mu[index]["nominal"] = [mu_p - half_c, mu_p + half_c]
            out_sigma[index]["nominal"] = [mu_c - half_p, mu_c + half_p]
        # endfor
    # endfor

    ensure_directory(output_dir)
    mu_path = output_dir / "hybrid_mu_polynomial_sigma_carbon_cuts.json"
    sigma_path = output_dir / "hybrid_mu_carbon_sigma_polynomial_cuts.json"
    write_json(mu_path, mu_poly_sigma_c)
    write_json(sigma_path, mu_c_sigma_poly)
    return mu_path, sigma_path

def load_channel_cuts(
    path: Path,
    cut_label: str = "nominal",
) -> dict[tuple[str, int], CutRecord]:
    payload = read_json_tolerating_trailing_escaped_whitespace(path)

    period_payload = payload.get("periods")
    if not isinstance(period_payload, dict):
        raise RuntimeError(
            f"Cut JSON {path} does not contain a 'periods' mapping."
        )
    # endif

    cuts: dict[tuple[str, int], CutRecord] = {}
    for period in PERIODS:
        rows = period_payload.get(period)
        if not isinstance(rows, list):
            raise RuntimeError(f"Cut JSON has no row list for {period}.")
        # endif

        for row in rows:
            x_index = int(row["x_index"])
            t_index = int(row["t_index"])
            bin_number = int(
                row.get("bin_number", combined_bin_number(x_index, t_index))
            )
            interval = row.get(cut_label)
            if not isinstance(interval, list) or len(interval) != 2:
                raise RuntimeError(
                    f"Missing {cut_label} cut interval for {period}, bin {bin_number}."
                )
            # endif
            low, high = float(interval[0]), float(interval[1])
            if not (
                math.isfinite(low)
                and math.isfinite(high)
                and high > low
            ):
                raise RuntimeError(
                    f"Invalid {cut_label} cut for {period}, bin {bin_number}: "
                    f"[{low}, {high}]."
                )
            # endif
            cuts[(period, bin_number)] = CutRecord(
                period=period,
                bin_number=bin_number,
                x_index=x_index,
                t_index=t_index,
                low_gev2=low,
                high_gev2=high,
            )
        # endfor
    # endfor

    expected = len(PERIODS) * NUMBER_OF_BINS
    if len(cuts) != expected:
        raise RuntimeError(
            f"Expected {expected} period/bin cuts, found {len(cuts)}."
        )
    # endif
    return cuts


# =============================================================================
# Dilution factors
# =============================================================================

def find_isr_dilution_json(directory: Path) -> Path:
    """Locate the compact dilution-factor JSON for the radiation sample."""
    preferred = [
        directory / "isr/dilution_factors_production.json",
        directory / "isr/tables/dilution_factors_production.json",
    ]
    for path in preferred:
        if path.is_file():
            return path
        # endif
    # endfor
    candidates = sorted(
        (directory / "isr").rglob("dilution_factors_production*.json")
        if (directory / "isr").is_dir() else [],
        key=lambda path: (path.stat().st_mtime, path.name),
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            "Could not locate the radiation-sample dilution-factor JSON under "
            f"{directory / 'isr'}. Run determine_dilution_factor.py with its "
            "ISR diagnostic enabled first, or pass --isr-dilution-json."
        )
    # endif
    return candidates[0]


def find_default_dilution_json(directory: Path) -> Path:
    if not directory.is_dir():
        raise FileNotFoundError(
            f"Missing dilution-factor output directory: {directory}"
        )
    # endif

    preferred = [
        directory / "production/tables/dilution_factors_production.json",
        directory / "production/dilution_factors_production.json",
        directory / "analysis/dilution_factors_production.json",
    ]
    for path in preferred:
        if path.is_file():
            return path
        # endif
    # endfor

    candidates = sorted(
        directory.rglob("dilution_factors_production*.json"),
        key=lambda path: (path.stat().st_mtime, path.name),
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            "Could not locate dilution_factors_production*.json under "
            f"{directory}."
        )
    # endif
    return candidates[0]


def _extract_recommended_record(cut_payload: Mapping[str, Any]) -> tuple[float, float]:
    record = cut_payload.get("recommended")
    if not isinstance(record, Mapping):
        raise RuntimeError("Dilution cut payload has no 'recommended' record.")
    # endif

    value_keys = ("value", "central", "recommended_value")
    uncertainty_keys = (
        "stat_uncertainty",
        "statistical_uncertainty",
        "recommended_stat_uncertainty",
    )

    value = None
    for key in value_keys:
        if key in record:
            value = float(record[key])
            break
        # endif
    # endfor

    uncertainty = None
    for key in uncertainty_keys:
        if key in record:
            uncertainty = float(record[key])
            break
        # endif
    # endfor

    if value is None or uncertainty is None:
        raise RuntimeError(
            "Recommended dilution record lacks value/stat_uncertainty fields: "
            f"{dict(record)}"
        )
    # endif
    return value, uncertainty


def load_dilution_factors(
    path: Path,
    cut_label: str = "nominal",
) -> dict[tuple[str, int], DilutionRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing dilution-factor JSON: {path}")
    # endif

    payload = json.loads(path.read_text(encoding="utf-8"))
    periods = payload.get("periods")
    if not isinstance(periods, dict):
        raise RuntimeError(
            f"Dilution JSON {path} does not contain a 'periods' mapping."
        )
    # endif

    records: dict[tuple[str, int], DilutionRecord] = {}
    for period in PERIODS:
        period_payload = periods.get(period)
        if not isinstance(period_payload, Mapping):
            raise RuntimeError(f"Dilution JSON has no payload for {period}.")
        # endif
        bins = period_payload.get("bins")
        if not isinstance(bins, list):
            raise RuntimeError(f"Dilution JSON has no bin list for {period}.")
        # endif

        for row in bins:
            x_index = int(row["x_index"])
            t_index = int(row["t_index"])
            bin_number = int(
                row.get("bin_number", combined_bin_number(x_index, t_index))
            )
            cuts = row.get("cuts")
            if not isinstance(cuts, Mapping) or cut_label not in cuts:
                raise RuntimeError(
                    f"Dilution JSON lacks {cut_label} cut for {period}, "
                    f"bin {bin_number}."
                )
            # endif
            value, uncertainty = _extract_recommended_record(cuts[cut_label])
            if not (
                math.isfinite(value)
                and math.isfinite(uncertainty)
                and value > 0.0
                and uncertainty >= 0.0
            ):
                raise RuntimeError(
                    f"Invalid dilution factor for {period}, bin {bin_number}: "
                    f"{value} +/- {uncertainty}."
                )
            # endif
            records[(period, bin_number)] = DilutionRecord(
                period=period,
                bin_number=bin_number,
                x_index=x_index,
                t_index=t_index,
                value=value,
                stat_uncertainty=uncertainty,
            )
        # endfor
    # endfor

    expected = len(PERIODS) * NUMBER_OF_BINS
    if len(records) != expected:
        raise RuntimeError(
            f"Expected {expected} dilution records, found {len(records)}."
        )
    # endif
    return records


# =============================================================================
# Input ROOT handling and event cache
# =============================================================================

def parse_input_override(text: str) -> tuple[str, Path]:
    try:
        period_text, path_text = text.split("=", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Input overrides must have the form PERIOD=/path/file.root"
        ) from exc
    # endtry
    period = period_text.strip().lower()
    if period not in PERIODS:
        raise argparse.ArgumentTypeError(
            f"Unknown period {period!r}; expected one of {PERIODS}."
        )
    # endif
    return period, Path(path_text).expanduser()


def resolve_tree(
    root_file: uproot.reading.ReadOnlyDirectory,
    requested_name: str,
    path: Path,
):
    """Resolve a ROOT tree without relying on Uproot membership testing."""
    lookup_names: list[str] = [requested_name]

    for key in root_file.keys():
        key_text = str(key)
        base_name = key_text.split(";")[0]
        if base_name == requested_name:
            lookup_names.extend((key_text, base_name))
        # endif
    # endfor

    attempted: set[str] = set()
    lookup_errors: list[str] = []

    for name in lookup_names:
        if name in attempted:
            continue
        # endif
        attempted.add(name)

        try:
            obj = root_file[name]
        except Exception as exc:
            lookup_errors.append(f"{name!r}: {exc}")
            continue
        # endtry

        if hasattr(obj, "arrays") and hasattr(obj, "keys"):
            return obj
        # endif
    # endfor

    tree_like: list[tuple[str, Any]] = []
    for key in root_file.keys():
        key_text = str(key)
        try:
            obj = root_file[key_text]
        except Exception:
            continue
        # endtry

        if hasattr(obj, "arrays") and hasattr(obj, "keys"):
            tree_like.append((key_text, obj))
        # endif
    # endfor

    if len(tree_like) == 1:
        return tree_like[0][1]
    # endif

    available = [str(key) for key in root_file.keys()]
    details = "; ".join(lookup_errors[:5])
    raise RuntimeError(
        f"Tree {requested_name!r} not found in {path}. "
        f"Available keys: {available}. "
        f"Direct-lookup diagnostics: {details}"
    )


def resolve_tree_branches(
    path: Path,
    tree_name: str,
) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing ROOT input: {path}")
    # endif

    with uproot.open(path) as root_file:
        tree = resolve_tree(root_file, tree_name, path)
        return {
            logical_name: resolve_branch(tree, aliases)
            for logical_name, aliases in BRANCH_ALIASES.items()
        }
    # endwith


def compute_theta_gamma(
    electron_momentum_gev: np.ndarray,
    electron_theta_rad: np.ndarray,
    beam_energy_gev: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute the virtual-photon angle relative to the incident beam.

    With the incident electron along +z,

        q_T = p_e' sin(theta_e),
        q_z = E_beam - p_e' cos(theta_e),

        sin(theta_gamma) = q_T / |q|,
        cos(theta_gamma) = q_z / |q|.

    This is algebraically equivalent to the analysis-note expression for
    theta_gamma and guarantees sin^2 + cos^2 = 1 numerically.
    """
    q_transverse = electron_momentum_gev * np.sin(electron_theta_rad)
    q_longitudinal = (
        beam_energy_gev
        - electron_momentum_gev * np.cos(electron_theta_rad)
    )
    q_magnitude = np.hypot(q_transverse, q_longitudinal)
    valid = np.isfinite(q_magnitude) & (q_magnitude > 0.0)

    sin_theta = np.full(q_magnitude.shape, np.nan, dtype=np.float64)
    cos_theta = np.full(q_magnitude.shape, np.nan, dtype=np.float64)
    sin_theta[valid] = q_transverse[valid] / q_magnitude[valid]
    cos_theta[valid] = q_longitudinal[valid] / q_magnitude[valid]

    sin_theta = np.clip(sin_theta, -1.0, 1.0)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return sin_theta, cos_theta


def build_event_cache(
    input_paths: Mapping[str, Path],
    tree_name: str,
    chunk_size: str,
    run_records: Mapping[int, RunRecord],
    cuts: Mapping[tuple[str, int], CutRecord],
    cache_path: Path,
) -> dict[str, Any]:
    ensure_directory(cache_path.parent)
    cache_build_start = time.perf_counter()
    print(
        f"[cache] BUILD START -> {cache_path.resolve()}",
        flush=True,
    )
    print(
        f"[cache] tree={tree_name}; chunk_size={chunk_size}; "
        f"periods={', '.join(PERIODS)}",
        flush=True,
    )

    collected: dict[str, list[np.ndarray]] = {
        "period_index": [],
        "runnum": [],
        "helicity": [],
        "bin_number": [],
        "xB": [],
        "minus_tprime": [],
        "minus_t": [],
        "W": [],
        "Mx2": [],
        "phi": [],
        "Q2": [],
        "epsilon": [],
        "DepA": [],
        "DepB": [],
        "DepC": [],
        "DepV": [],
        "DepW": [],
        "sin_theta_gamma": [],
        "cos_theta_gamma": [],
        "rB": [],
        "rC": [],
        "rV": [],
        "rW": [],
    }

    period_statistics: dict[str, dict[str, int]] = {}

    for period in PERIODS:
        period_start = time.perf_counter()
        path = input_paths[period].expanduser().resolve()
        print(f"[cache] {period}: opening {path}", flush=True)
        if not path.is_file():
            raise FileNotFoundError(f"Missing ROOT input: {path}")
        # endif

        total_seen = 0
        total_selected = 0
        total_zero_helicity = 0
        total_missing_run = 0
        total_zero_charge_run_events = 0
        missing_run_event_counts: dict[int, int] = {}
        zero_charge_run_event_counts: dict[int, int] = {}

        # Iterate the resolved TTree object directly.  The external-radiation
        # files can carry a nonstandard key such as ``PhysicsEvents;1;1``;
        # reconstructing a ``path:PhysicsEvents`` source string loses that
        # resolved key and makes Uproot fail even though the tree was found.
        with uproot.open(path) as root_file:
            tree = resolve_tree(root_file, tree_name, path)
            branches = {
                logical_name: resolve_branch(tree, aliases)
                for logical_name, aliases in BRANCH_ALIASES.items()
            }
            expressions = list(dict.fromkeys(branches.values()))

            chunk_number = 0
            for arrays in tree.iterate(
                expressions=expressions,
                step_size=chunk_size,
                library="np",
            ):
                chunk_number += 1
                n_chunk = len(arrays[branches["runnum"]])
                total_seen += n_chunk
                if chunk_number == 1 or chunk_number % 10 == 0:
                    print(
                        f"[cache] {period}: chunk {chunk_number}; "
                        f"seen={total_seen:,}; selected={total_selected:,}; "
                        f"elapsed={time.perf_counter() - period_start:.1f} s",
                        flush=True,
                    )
                # endif

                runnum = np.asarray(
                    arrays[branches["runnum"]],
                    dtype=np.int64,
                )
                helicity = np.asarray(
                    arrays[branches["helicity"]],
                    dtype=np.int8,
                )
                xB = np.asarray(arrays[branches["xB"]], dtype=np.float64)
                minus_tprime = -np.asarray(
                    arrays[branches["tprime"]],
                    dtype=np.float64,
                )
                minus_t = -np.asarray(
                    arrays[branches["t"]],
                    dtype=np.float64,
                )
                w = np.asarray(
                    arrays[branches["W"]],
                    dtype=np.float64,
                )
                mx2 = np.asarray(arrays[branches["Mx2"]], dtype=np.float64)
                phi_raw = np.asarray(arrays[branches["phi"]], dtype=np.float64)
                phi = angle_to_radians(phi_raw, branches["phi"])
                phi = np.mod(phi, 2.0 * math.pi)
                q2 = np.asarray(arrays[branches["Q2"]], dtype=np.float64)
                y = np.asarray(arrays[branches["y"]], dtype=np.float64)
                gamma2 = 4.0 * PROTON_MASS_GEV**2 * xB**2 / q2
                epsilon = (1.0 - y - 0.25 * gamma2 * y**2) / (
                    1.0 - y + 0.5 * y**2 + 0.25 * gamma2 * y**2
                )
                e_p = np.asarray(arrays[branches["e_p"]], dtype=np.float64)
                e_theta_raw = np.asarray(
                    arrays[branches["e_theta"]],
                    dtype=np.float64,
                )
                e_theta = angle_to_radians(e_theta_raw, branches["e_theta"])

                dep_a = np.asarray(arrays[branches["DepA"]], dtype=np.float64)
                dep_b = np.asarray(arrays[branches["DepB"]], dtype=np.float64)
                dep_c = np.asarray(arrays[branches["DepC"]], dtype=np.float64)
                dep_v = np.asarray(arrays[branches["DepV"]], dtype=np.float64)
                dep_w = np.asarray(arrays[branches["DepW"]], dtype=np.float64)

                sin_theta_gamma, cos_theta_gamma = compute_theta_gamma(
                    e_p,
                    e_theta,
                    BEAM_ENERGY_GEV[period],
                )

                _, _, bin_number = bin_indices(xB, minus_tprime)

                known_run = np.fromiter(
                    (
                        int(run) in run_records
                        and run_records[int(run)].period == period
                        for run in runnum
                    ),
                    count=runnum.size,
                    dtype=bool,
                )
                active_run = np.fromiter(
                    (
                        int(run) in run_records
                        and run_records[int(run)].period == period
                        and run_records[int(run)].charge_plus > 0.0
                        and run_records[int(run)].charge_minus > 0.0
                        for run in runnum
                    ),
                    count=runnum.size,
                    dtype=bool,
                )
                zero_charge_run = known_run & ~active_run

                missing_mask = ~known_run
                total_missing_run += int(np.count_nonzero(missing_mask))
                total_zero_charge_run_events += int(
                    np.count_nonzero(zero_charge_run)
                )
                total_zero_helicity += int(np.count_nonzero(helicity == 0))

                if np.any(missing_mask):
                    missing_runs, missing_counts = np.unique(
                        runnum[missing_mask],
                        return_counts=True,
                    )
                    for missing_run, missing_count in zip(
                        missing_runs,
                        missing_counts,
                    ):
                        key = int(missing_run)
                        missing_run_event_counts[key] = (
                            missing_run_event_counts.get(key, 0)
                            + int(missing_count)
                        )
                    # endfor
                # endif

                if np.any(zero_charge_run):
                    disabled_runs, disabled_counts = np.unique(
                        runnum[zero_charge_run],
                        return_counts=True,
                    )
                    for disabled_run, disabled_count in zip(
                        disabled_runs,
                        disabled_counts,
                    ):
                        key = int(disabled_run)
                        zero_charge_run_event_counts[key] = (
                            zero_charge_run_event_counts.get(key, 0)
                            + int(disabled_count)
                        )
                    # endfor
                # endif

                finite = (
                    np.isfinite(xB)
                    & np.isfinite(minus_tprime)
                    & np.isfinite(minus_t)
                    & np.isfinite(w)
                    & np.isfinite(mx2)
                    & np.isfinite(phi)
                    & np.isfinite(q2)
                    & np.isfinite(epsilon)
                    & np.isfinite(sin_theta_gamma)
                    & np.isfinite(cos_theta_gamma)
                    & np.isfinite(dep_a)
                    & np.isfinite(dep_b)
                    & np.isfinite(dep_c)
                    & np.isfinite(dep_v)
                    & np.isfinite(dep_w)
                )
                # Re-impose the DIS W cut on the kinematics stored in the
                # final analysis trees.  These values include the established
                # momentum corrections, so this explicitly guarantees W > 2
                # GeV after correction rather than relying on the upstream
                # pre-correction event selection.
                base = (
                    finite
                    & active_run
                    & np.isin(helicity, (-1, 1))
                    & (bin_number >= 1)
                    & (w > DIS_W_MIN_GEV)
                    & (dep_a > 0.0)
                )

                selected = np.zeros(base.shape, dtype=bool)
                for current_bin in range(1, NUMBER_OF_BINS + 1):
                    cut = cuts[(period, current_bin)]
                    selected |= (
                        base
                        & (bin_number == current_bin)
                        & (mx2 >= cut.low_gev2)
                        & (mx2 < cut.high_gev2)
                    )
                # endfor

                if not np.any(selected):
                    continue
                # endif

                total_selected += int(np.count_nonzero(selected))
                collected["period_index"].append(
                    np.full(
                        int(np.count_nonzero(selected)),
                        PERIOD_INDEX[period],
                        dtype=np.int8,
                    )
                )
                collected["runnum"].append(runnum[selected].astype(np.int32))
                collected["helicity"].append(helicity[selected].astype(np.int8))
                collected["bin_number"].append(
                    bin_number[selected].astype(np.int16)
                )
                collected["xB"].append(xB[selected])
                collected["minus_tprime"].append(minus_tprime[selected])
                collected["minus_t"].append(minus_t[selected])
                collected["W"].append(w[selected])
                collected["Mx2"].append(mx2[selected])
                collected["phi"].append(phi[selected])
                collected["Q2"].append(q2[selected])
                collected["epsilon"].append(epsilon[selected])
                collected["DepA"].append(dep_a[selected])
                collected["DepB"].append(dep_b[selected])
                collected["DepC"].append(dep_c[selected])
                collected["DepV"].append(dep_v[selected])
                collected["DepW"].append(dep_w[selected])
                collected["sin_theta_gamma"].append(sin_theta_gamma[selected])
                collected["cos_theta_gamma"].append(cos_theta_gamma[selected])
                collected["rB"].append(dep_b[selected] / dep_a[selected])
                collected["rC"].append(dep_c[selected] / dep_a[selected])
                collected["rV"].append(dep_v[selected] / dep_a[selected])
                collected["rW"].append(dep_w[selected] / dep_a[selected])
            # endfor

            if missing_run_event_counts:
                details = ", ".join(
                    f"{run}: {count:,} events"
                    for run, count in sorted(missing_run_event_counts.items())
                )
                raise RuntimeError(
                    f"{period} ROOT input contains events from runs absent from "
                    f"the NH3 section of the run-information CSV: {details}."
                )
            # endif

            if zero_charge_run_event_counts:
                details = ", ".join(
                    f"{run}: {count:,} events"
                    for run, count in sorted(
                        zero_charge_run_event_counts.items()
                    )
                )
                raise RuntimeError(
                    f"{period} ROOT input contains events from QA-disabled runs "
                    f"whose CSV charges are Q+=Q-=0: {details}. Zero-charge rows "
                    "are allowed in clas12_run_info.csv, but no event from those "
                    "runs may enter the ROOT input used for asymmetry extraction."
                )
            # endif

        print(
            f"[cache] {period}: DONE; seen={total_seen:,}; "
            f"selected={total_selected:,}; "
            f"elapsed={time.perf_counter() - period_start:.1f} s",
            flush=True,
        )
        period_statistics[period] = {
            "events_seen": total_seen,
            "events_selected": total_selected,
            "helicity_zero_events_ignored": total_zero_helicity,
            "events_with_missing_run_info": total_missing_run,
            "events_from_zero_charge_runs": total_zero_charge_run_events,
            "active_runs_in_likelihood": sum(
                1
                for record in run_records.values()
                if (
                    record.period == period
                    and record.charge_plus > 0.0
                    and record.charge_minus > 0.0
                )
            ),
        }
    # endfor

    if not collected["runnum"]:
        raise RuntimeError("No events survived the nominal selection.")
    # endif

    cache = {
        name: np.concatenate(chunks)
        for name, chunks in collected.items()
    }
    print(
        f"[cache] compressing {cache["runnum"].size:,} selected events -> "
        f"{cache_path.resolve()}",
        flush=True,
    )
    np.savez_compressed(cache_path, **cache)
    print(
        f"[cache] BUILD DONE; selected={cache['runnum'].size:,}; "
        f"elapsed={time.perf_counter() - cache_build_start:.1f} s",
        flush=True,
    )

    return {
        "cache_path": str(cache_path.resolve()),
        "number_of_selected_events": int(cache["runnum"].size),
        "period_statistics": period_statistics,
    }



def derive_event_cache(
    *,
    source_cache_path: Path,
    cuts: Mapping[tuple[str, int], CutRecord],
    cache_path: Path,
) -> dict[str, Any]:
    """Filter a previously selected superset cache without rereading ROOT."""
    derive_start = time.perf_counter()
    print(
        f"[cache] DERIVE START: {source_cache_path.resolve()} -> "
        f"{cache_path.resolve()}",
        flush=True,
    )
    source = load_event_cache(source_cache_path)
    print(
        f"[cache] source cache loaded: {source['runnum'].size:,} events",
        flush=True,
    )
    selected = np.zeros(source["runnum"].shape, dtype=bool)
    period_index = np.asarray(source["period_index"], dtype=np.int8)
    bin_number = np.asarray(source["bin_number"], dtype=np.int16)
    mx2 = np.asarray(source["Mx2"], dtype=np.float64)

    for period in PERIODS:
        period_mask = period_index == PERIOD_INDEX[period]
        for current_bin in range(1, NUMBER_OF_BINS + 1):
            cut = cuts[(period, current_bin)]
            selected |= (
                period_mask
                & (bin_number == current_bin)
                & (mx2 >= cut.low_gev2)
                & (mx2 < cut.high_gev2)
            )
        # endfor
    # endfor

    if not np.any(selected):
        raise RuntimeError(
            f"No events survived cache-derived selection for {cache_path}."
        )
    # endif

    ensure_directory(cache_path.parent)
    filtered = {name: values[selected] for name, values in source.items()}
    np.savez_compressed(cache_path, **filtered)
    print(
        f"[cache] DERIVE DONE; selected={np.count_nonzero(selected):,}/"
        f"{selected.size:,}; elapsed={time.perf_counter() - derive_start:.1f} s",
        flush=True,
    )
    return {
        "cache_path": str(cache_path.resolve()),
        "number_of_selected_events": int(np.count_nonzero(selected)),
        "number_of_source_events": int(selected.size),
        "derived_from_cache": str(source_cache_path.resolve()),
    }


def load_event_cache(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing event cache: {path}")
    # endif
    with np.load(path, allow_pickle=False) as payload:
        result = {
            key: np.asarray(payload[key])
            for key in payload.files
        }
    # endwith

    sizes = {array.shape[0] for array in result.values()}
    if len(sizes) != 1:
        raise RuntimeError(
            f"Cache arrays have inconsistent lengths: {sorted(sizes)}"
        )
    # endif
    return result


# =============================================================================
# Likelihood
# =============================================================================

def external_data_informed_transverse_terms(
    phi: np.ndarray,
    r_b: np.ndarray,
    r_c: np.ndarray,
    r_v: np.ndarray,
    r_w: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return the fixed external-data-informed transverse leakage terms.

    The restricted harmonic mapping is

      F_UT^{sin(phi-phi_S)} / F_UU
      F_UT^{sin(phi+phi_S)} / F_UU
      F_UT^{sin(2phi-phi_S)} / F_UU
      F_LT^{cos(phi_S)} / F_UU
      F_LT^{cos(phi-phi_S)} / F_UU

    These are the complete fixed transverse inputs used in the controlled
    external-data-informed leakage study.  No UT or LT amplitude is fitted.

    With phi_S = 0, sin(phi-phi_S) and sin(phi+phi_S) both reduce
    to sin(phi), but they are not averaged: their cross-section
    contributions add.  In this convention the minus term carries no
    extra depolarization factor, while the plus term carries r_b = B/A.
    The A_UT,l^{sin(2phi-phi_S)} term carries r_v.
    """
    ut_sin_phi = float(EXTERNAL_TRANSVERSE_INPUTS["ut_sin_phi"]["value"])
    ut_sin_phi_plus = float(
        EXTERNAL_TRANSVERSE_INPUTS["ut_sin_phi_plus"]["value"]
    )
    ut_sin_2phi = float(
        EXTERNAL_TRANSVERSE_INPUTS["ut_sin_2phi"]["value"]
    )
    lt_cos_0phi = float(EXTERNAL_TRANSVERSE_INPUTS["lt_cos_0phi"]["value"])
    lt_cos_phi = float(EXTERNAL_TRANSVERSE_INPUTS["lt_cos_phi"]["value"])

    ut = (
        (ut_sin_phi + r_b * ut_sin_phi_plus) * np.sin(phi)
        + r_v * ut_sin_2phi * np.sin(2.0 * phi)
    )
    lt = (
        r_w * lt_cos_0phi
        + r_c * lt_cos_phi * np.cos(phi)
    )
    return ut, lt


def evaluate_cross_section_factor(
    *,
    variant: str,
    phi: np.ndarray,
    r_b: np.ndarray,
    r_c: np.ndarray,
    r_v: np.ndarray,
    r_w: np.ndarray,
    sin_theta_gamma: np.ndarray,
    cos_theta_gamma: np.ndarray,
    helicity: np.ndarray | float,
    beam_polarization: float,
    target_polarization: np.ndarray | float,
    dilution: float,
    transverse_scales: Mapping[str, float] | None,
    u1: float,
    u2: float,
    lu1: float,
    ul1: float,
    ul2: float,
    ll0: float,
    ll1: float,
) -> np.ndarray:
    h = np.asarray(helicity, dtype=np.float64)
    pt = np.asarray(target_polarization, dtype=np.float64)

    if variant == "photon_axis_projection":
        # Interpretation study: project the laboratory target polarization
        # onto the virtual-photon direction and neglect the transverse piece.
        p_longitudinal = pt * cos_theta_gamma
        p_transverse = np.zeros_like(pt, dtype=np.float64)
    elif variant == "external_data_informed":
        # Interpretation study in the virtual-photon basis, including the
        # geometrically induced transverse component with fixed external inputs.
        p_longitudinal = pt * cos_theta_gamma
        p_transverse = pt * sin_theta_gamma
    else:
        # Production nominal: report the observable for the polarization state
        # actually prepared in the laboratory, longitudinal to the beam line.
        # No cos(theta_gamma) projection is applied to the fitted target
        # polarization.  Photon-axis decompositions are diagnostic studies only.
        p_longitudinal = pt
        p_transverse = np.zeros_like(pt, dtype=np.float64)
    # endif

    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    cos_2phi = np.cos(2.0 * phi)
    sin_2phi = np.sin(2.0 * phi)

    unpolarized = 1.0 + r_v * u1 * cos_phi + r_b * u2 * cos_2phi
    beam_spin = h * beam_polarization * r_w * lu1 * sin_phi
    target_longitudinal = (
        dilution
        * p_longitudinal
        * (r_v * ul1 * sin_phi + r_b * ul2 * sin_2phi)
    )
    double_longitudinal = (
        h
        * beam_polarization
        * dilution
        * p_longitudinal
        * (r_c * ll0 + r_w * ll1 * cos_phi)
    )

    if variant == "external_data_informed":
        ut_fixed, lt_fixed = external_data_informed_transverse_terms(
            phi,
            r_b,
            r_c,
            r_v,
            r_w,
        )
    else:
        ut_fixed = np.zeros(phi.shape, dtype=np.float64)
        lt_fixed = np.zeros(phi.shape, dtype=np.float64)
    # endif

    # UT and LT amplitudes are not fitted observables in this analysis.
    # They enter only as fixed external inputs in the external-data-informed
    # target-axis leakage study.
    target_transverse = dilution * p_transverse * ut_fixed
    double_transverse = (
        h * beam_polarization * dilution * p_transverse * lt_fixed
    )

    return (
        unpolarized
        + beam_spin
        + target_longitudinal
        + double_longitudinal
        + target_transverse
        + double_transverse
    )


def make_bin_nll(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    bin_number: int,
    variant: str,
    active_periods: tuple[str, ...] = PERIODS,
    transverse_scales: Mapping[str, float] | None = None,
    target_sign_filter: int | None = None,
    run_ranges_filter: tuple[tuple[int, int], ...] | None = None,
    beam_polarization_scales: Mapping[str, float] | None = None,
    double_spin_products_by_run: Mapping[int, float] | None = None,
):
    beam_polarization_scales = dict(beam_polarization_scales or {})
    double_spin_products_by_run = dict(double_spin_products_by_run or {})
    mask = events["bin_number"] == bin_number
    if len(active_periods) != len(PERIODS):
        allowed_indices = np.asarray(
            [PERIOD_INDEX[period] for period in active_periods],
            dtype=np.int8,
        )
        mask &= np.isin(events["period_index"], allowed_indices)
    # endif

    if run_ranges_filter is not None:
        run_mask = np.zeros(mask.shape, dtype=bool)
        for run_low, run_high in run_ranges_filter:
            run_mask |= (events["runnum"] >= run_low) & (events["runnum"] <= run_high)
        # endfor
        mask &= run_mask
    # endif

    # Optional fixed-target-orientation diagnostic.  Determine the target sign
    # from the run-state polarization associated with each event; this keeps the
    # ordinary event selection untouched while restricting the conditional
    # likelihood to one physical target-polarization orientation.
    if target_sign_filter is not None:
        requested_sign = int(np.sign(target_sign_filter))
        if requested_sign == 0:
            raise ValueError("target_sign_filter must be +1, -1, or None.")
        # endif
        event_target_sign = np.zeros(mask.shape, dtype=np.int8)
        candidate_indices = np.flatnonzero(mask)
        for period in active_periods:
            period_candidates = candidate_indices[
                events["period_index"][candidate_indices] == PERIOD_INDEX[period]
            ]
            if period_candidates.size == 0:
                continue
            # endif
            state = run_states[period]
            run_to_sign = {
                int(run): int(np.sign(pt))
                for run, pt in zip(state["run"], state["pt"])
            }
            event_target_sign[period_candidates] = np.fromiter(
                (run_to_sign[int(run)] for run in events["runnum"][period_candidates]),
                count=period_candidates.size,
                dtype=np.int8,
            )
        # endfor
        mask &= event_target_sign == requested_sign
    # endif

    if not np.any(mask):
        raise RuntimeError(
            f"Bin {bin_number} has no selected events for {active_periods}."
        )
    # endif

    period_index = events["period_index"][mask].astype(np.int8, copy=False)
    runnum = events["runnum"][mask].astype(np.int32, copy=False)
    helicity = events["helicity"][mask].astype(np.int8, copy=False)
    phi = events["phi"][mask].astype(np.float64, copy=False)
    sin_theta = events["sin_theta_gamma"][mask].astype(np.float64, copy=False)
    cos_theta = events["cos_theta_gamma"][mask].astype(np.float64, copy=False)
    r_b = events["rB"][mask].astype(np.float64, copy=False)
    r_c = events["rC"][mask].astype(np.float64, copy=False)
    r_v = events["rV"][mask].astype(np.float64, copy=False)
    r_w = events["rW"][mask].astype(np.float64, copy=False)

    period_event_indices = {
        period: np.flatnonzero(period_index == PERIOD_INDEX[period])
        for period in active_periods
    }

    run_lookup = {
        period: {
            int(run): index
            for index, run in enumerate(run_states[period]["run"])
        }
        for period in active_periods
    }

    observed_run_indices: dict[str, np.ndarray] = {}
    for period in active_periods:
        indices = period_event_indices[period]
        observed_run_indices[period] = np.fromiter(
            (run_lookup[period][int(run)] for run in runnum[indices]),
            count=indices.size,
            dtype=np.int32,
        )
    # endfor

    dilution_central = {
        period: dilution_records[(period, bin_number)].value
        for period in active_periods
    }
    dilution_sigma = {
        period: dilution_records[(period, bin_number)].stat_uncertainty
        for period in active_periods
    }

    # Precompute all parameter-independent event arrays and run-state charge
    # moments once.  The previous implementation recomputed several of these
    # quantities at every Minuit function call.
    period_data: dict[str, dict[str, Any]] = {}
    for period in active_periods:
        indices = period_event_indices[period]
        event_phi = phi[indices]
        state = run_states[period]
        state_pt = state["pt"]
        state_q_plus = state["q_plus"]
        state_q_minus = state["q_minus"]
        observed_state_index = observed_run_indices[period]
        state_use = np.ones(state_pt.shape, dtype=bool)
        if run_ranges_filter is not None:
            state_run = state["run"]
            state_use = np.zeros(state_pt.shape, dtype=bool)
            for run_low, run_high in run_ranges_filter:
                state_use |= (state_run >= run_low) & (state_run <= run_high)
            # endfor
        # endif
        if target_sign_filter is not None:
            state_use &= np.sign(state_pt) == int(np.sign(target_sign_filter))
        # endif
        event_h = helicity[indices].astype(np.float64, copy=False)

        period_data[period] = {
            "phi": event_phi,
            "sin_phi": np.sin(event_phi),
            "cos_phi": np.cos(event_phi),
            "sin_2phi": np.sin(2.0 * event_phi),
            "cos_2phi": np.cos(2.0 * event_phi),
            "sin_3phi": np.sin(3.0 * event_phi),
            "r_b": r_b[indices],
            "r_c": r_c[indices],
            "r_v": r_v[indices],
            "r_w": r_w[indices],
            "sin_theta": sin_theta[indices],
            "cos_theta": cos_theta[indices],
            "h": event_h,
            "observed_pt": state_pt[observed_state_index],
            "observed_pbpt": np.asarray([
                double_spin_products_by_run.get(
                    int(run),
                    BEAM_POLARIZATION[period] * state_pt[state_index],
                )
                for run, state_index in zip(
                    runnum[indices], observed_state_index
                )
            ], dtype=np.float64),
            "double_charge_sum": float(np.sum([
                (state_q_plus[j] - state_q_minus[j])
                * double_spin_products_by_run.get(
                    int(state["run"][j]),
                    BEAM_POLARIZATION[period] * state_pt[j],
                )
                for j in np.flatnonzero(state_use)
            ])),
            "observed_charge": np.where(
                event_h > 0.0,
                state_q_plus[observed_state_index],
                state_q_minus[observed_state_index],
            ),
            "charge_sum": float(np.sum((state_q_plus + state_q_minus)[state_use])),
            "helicity_charge_sum": float(
                np.sum((state_q_plus - state_q_minus)[state_use])
            ),
            "target_charge_sum": float(
                np.sum((state_pt * (state_q_plus + state_q_minus))[state_use])
            ),
            "helicity_target_charge_sum": float(
                np.sum((state_pt * (state_q_plus - state_q_minus))[state_use])
            ),
        }
    # endfor

    def nll(
        u1: float,
        u2: float,
        lu1: float,
        ul1: float,
        ul2: float,
        ll0: float,
        ll1: float,
        f_su22: float,
        f_fa22: float,
        f_sp23: float,
    ) -> float:
        parameter_values = (u1, u2, lu1, ul1, ul2, ll0, ll1)
        if not all(math.isfinite(value) for value in parameter_values):
            return INVALID_NLL
        # endif

        dilution_by_period = {
            "su22": f_su22,
            "fa22": f_fa22,
            "sp23": f_sp23,
        }
        total = 0.0

        for period in active_periods:
            dilution = dilution_by_period[period]
            sigma = dilution_sigma[period]
            central = dilution_central[period]
            if not math.isfinite(dilution) or dilution <= 0.0:
                return INVALID_NLL
            # endif
            if sigma > 0.0:
                total += 0.5 * ((dilution - central) / sigma) ** 2
            elif abs(dilution - central) > 1.0e-14:
                return INVALID_NLL
            # endif

            data = period_data[period]
            event_phi = data["phi"]
            event_r_b = data["r_b"]
            event_r_c = data["r_c"]
            event_r_v = data["r_v"]
            event_r_w = data["r_w"]
            event_sin = data["sin_theta"]
            event_cos = data["cos_theta"]
            event_h = data["h"]

            numerator_factor = evaluate_cross_section_factor(
                variant=variant,
                phi=event_phi,
                r_b=event_r_b,
                r_c=event_r_c,
                r_v=event_r_v,
                r_w=event_r_w,
                sin_theta_gamma=event_sin,
                cos_theta_gamma=event_cos,
                helicity=event_h,
                beam_polarization=(BEAM_POLARIZATION[period] * beam_polarization_scales.get(period, 1.0)),
                target_polarization=data["observed_pt"],
                dilution=dilution,
                transverse_scales=transverse_scales,
                u1=u1,
                u2=u2,
                lu1=lu1,
                ul1=ul1,
                ul2=ul2,
                ll0=ll0,
                ll1=ll1,
            )
            if double_spin_products_by_run:
                # Replace only the longitudinal double-spin Pb*Pt factor.
                # LU continues to use the nominal beam polarization and UL
                # continues to use the nominal target polarization.
                nominal_pbpt = (
                    BEAM_POLARIZATION[period]
                    * beam_polarization_scales.get(period, 1.0)
                    * data["observed_pt"]
                )
                if variant == "photon_axis_projection" or variant == "external_data_informed":
                    replacement_longitudinal_geometry = event_cos
                else:
                    replacement_longitudinal_geometry = np.ones_like(event_phi)
                # endif
                replacement = (
                    event_h * dilution * replacement_longitudinal_geometry
                    * (data["observed_pbpt"] - nominal_pbpt)
                    * (event_r_c * ll0 + event_r_w * ll1 * data["cos_phi"])
                )
                numerator_factor = numerator_factor + replacement
            # endif

            if (
                np.any(~np.isfinite(numerator_factor))
                or np.any(numerator_factor <= CROSS_SECTION_FLOOR)
            ):
                return INVALID_NLL
            # endif

            unpolarized = (
                1.0
                + event_r_v * u1 * data["cos_phi"]
                + event_r_b * u2 * data["cos_2phi"]
            )
            beam_coefficient = (
                (BEAM_POLARIZATION[period] * beam_polarization_scales.get(period, 1.0))
                * event_r_w
                * lu1
                * data["sin_phi"]
            )

            if variant == "photon_axis_projection":
                longitudinal_geometry = event_cos
                transverse_geometry = np.zeros_like(event_phi)
            elif variant == "external_data_informed":
                longitudinal_geometry = event_cos
                transverse_geometry = event_sin
            else:
                # Production nominal: the target polarization is defined along
                # the incident beam in the laboratory frame.  Keep the same
                # event-by-event SIDIS depolarization factors, but do not apply
                # an additional cos(theta_gamma) target-axis projection.
                longitudinal_geometry = np.ones_like(event_phi)
                transverse_geometry = np.zeros_like(event_phi)
            # endif

            target_longitudinal_coefficient = (
                dilution
                * longitudinal_geometry
                * (
                    event_r_v * ul1 * data["sin_phi"]
                    + event_r_b * ul2 * data["sin_2phi"]
                )
            )
            double_longitudinal_coefficient = (
                (BEAM_POLARIZATION[period] * beam_polarization_scales.get(period, 1.0))
                * dilution
                * longitudinal_geometry
                * (
                    event_r_c * ll0
                    + event_r_w * ll1 * data["cos_phi"]
                )
            )

            # UT/LT terms are fixed external leakage inputs only; there
            # are no fitted transverse-target amplitudes.
            target_transverse_coefficient = np.zeros_like(event_phi)
            double_transverse_coefficient = np.zeros_like(event_phi)

            if variant == "external_data_informed":
                ut_sin_phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["ut_sin_phi"]["value"]
                )
                ut_sin_phi_plus = float(
                    EXTERNAL_TRANSVERSE_INPUTS["ut_sin_phi_plus"]["value"]
                )
                ut_sin_2phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["ut_sin_2phi"]["value"]
                )
                lt_cos_0phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["lt_cos_0phi"]["value"]
                )
                lt_cos_phi = float(
                    EXTERNAL_TRANSVERSE_INPUTS["lt_cos_phi"]["value"]
                )
                ut_fixed = (
                    (ut_sin_phi + event_r_b * ut_sin_phi_plus)
                    * data["sin_phi"]
                    + event_r_v * ut_sin_2phi * data["sin_2phi"]
                )
                lt_fixed = (
                    event_r_w * lt_cos_0phi
                    + event_r_c * lt_cos_phi * data["cos_phi"]
                )
                target_transverse_coefficient = (
                    target_transverse_coefficient
                    + dilution * transverse_geometry * ut_fixed
                )
                double_transverse_coefficient = (
                    double_transverse_coefficient
                    + (BEAM_POLARIZATION[period] * beam_polarization_scales.get(period, 1.0))
                    * dilution
                    * transverse_geometry
                    * lt_fixed
                )
            # endif

            target_coefficient = (
                target_longitudinal_coefficient
                + target_transverse_coefficient
            )
            double_coefficient = (
                double_longitudinal_coefficient
                + double_transverse_coefficient
            )

            denominator = (
                data["charge_sum"] * unpolarized
                + data["helicity_charge_sum"] * beam_coefficient
                + data["target_charge_sum"] * target_coefficient
                + (
                    data["double_charge_sum"]
                    * dilution
                    * longitudinal_geometry
                    * (event_r_c * ll0 + event_r_w * ll1 * data["cos_phi"])
                    if double_spin_products_by_run
                    else data["helicity_target_charge_sum"] * double_coefficient
                )
            )
            if (
                np.any(~np.isfinite(denominator))
                or np.any(denominator <= CROSS_SECTION_FLOOR)
            ):
                return INVALID_NLL
            # endif

            numerator = data["observed_charge"] * numerator_factor
            probability = numerator / denominator
            if (
                np.any(~np.isfinite(probability))
                or np.any(probability <= 0.0)
                or np.any(probability > 1.0 + 1.0e-10)
            ):
                return INVALID_NLL
            # endif
            total -= float(
                np.sum(np.log(np.maximum(probability, PROBABILITY_FLOOR)))
            )
        # endfor

        return total

    metadata = {
        "active_periods": list(active_periods),
        "target_sign_filter": target_sign_filter,
        "number_of_events": int(np.count_nonzero(mask)),
        "events_by_period": {
            period: int(period_event_indices.get(period, np.empty(0)).size)
            for period in PERIODS
        },
        "dilution_central": {
            period: (
                dilution_records[(period, bin_number)].value
                if period in active_periods else None
            )
            for period in PERIODS
        },
        "dilution_stat_uncertainty": {
            period: (
                dilution_records[(period, bin_number)].stat_uncertainty
                if period in active_periods else None
            )
            for period in PERIODS
        },
        # Publication-facing target-axis kinematic quantity.  Tabulate
        # <sin(theta_gamma)> directly rather than converting the sample to a
        # mean theta_gamma.
        "mean_sin_theta_gamma": float(np.mean(sin_theta)),
        "rms_sin_theta_gamma": float(np.std(sin_theta, ddof=1))
        if sin_theta.size > 1 else 0.0,
        "mean_cos_theta_gamma": float(np.mean(cos_theta)),
        "rms_cos_theta_gamma": float(np.std(cos_theta, ddof=1))
        if cos_theta.size > 1 else 0.0,
        "min_xB": float(np.min(events["xB"][mask])),
        "median_xB": float(np.median(events["xB"][mask])),
        "mean_xB": float(np.mean(events["xB"][mask])),
        "max_xB": float(np.max(events["xB"][mask])),
        "min_Q2_gev2": float(np.min(events["Q2"][mask])),
        "median_Q2_gev2": float(np.median(events["Q2"][mask])),
        "mean_Q2_gev2": float(np.mean(events["Q2"][mask])),
        "max_Q2_gev2": float(np.max(events["Q2"][mask])),
        "min_W_gev": float(np.min(events["W"][mask])),
        "median_W_gev": float(np.median(events["W"][mask])),
        "mean_W_gev": float(np.mean(events["W"][mask])),
        "max_W_gev": float(np.max(events["W"][mask])),
        "min_minus_t_gev2": float(np.min(events["minus_t"][mask])),
        "median_minus_t_gev2": float(np.median(events["minus_t"][mask])),
        "mean_minus_t_gev2": float(np.mean(events["minus_t"][mask])),
        "max_minus_t_gev2": float(np.max(events["minus_t"][mask])),
        "min_minus_tprime_gev2": float(np.min(events["minus_tprime"][mask])),
        "median_minus_tprime_gev2": float(np.median(events["minus_tprime"][mask])),
        "mean_minus_tprime_gev2": float(np.mean(events["minus_tprime"][mask])),
        "max_minus_tprime_gev2": float(np.max(events["minus_tprime"][mask])),
        "min_epsilon": float(np.min(events["epsilon"][mask])),
        "median_epsilon": float(np.median(events["epsilon"][mask])),
        "mean_epsilon": float(np.mean(events["epsilon"][mask])),
        "max_epsilon": float(np.max(events["epsilon"][mask])),
        "min_DepA": float(np.min(events["DepA"][mask])),
        "median_DepA": float(np.median(events["DepA"][mask])),
        "mean_DepA": float(np.mean(events["DepA"][mask])),
        "max_DepA": float(np.max(events["DepA"][mask])),
        "min_DepB": float(np.min(events["DepB"][mask])),
        "median_DepB": float(np.median(events["DepB"][mask])),
        "mean_DepB": float(np.mean(events["DepB"][mask])),
        "max_DepB": float(np.max(events["DepB"][mask])),
        "min_DepC": float(np.min(events["DepC"][mask])),
        "median_DepC": float(np.median(events["DepC"][mask])),
        "mean_DepC": float(np.mean(events["DepC"][mask])),
        "max_DepC": float(np.max(events["DepC"][mask])),
        "min_DepV": float(np.min(events["DepV"][mask])),
        "median_DepV": float(np.median(events["DepV"][mask])),
        "mean_DepV": float(np.mean(events["DepV"][mask])),
        "max_DepV": float(np.max(events["DepV"][mask])),
        "min_DepW": float(np.min(events["DepW"][mask])),
        "median_DepW": float(np.median(events["DepW"][mask])),
        "mean_DepW": float(np.mean(events["DepW"][mask])),
        "max_DepW": float(np.max(events["DepW"][mask])),
    }
    return nll, metadata


def minuit_result_payload(minuit: Minuit) -> dict[str, Any]:
    parameter_names = list(minuit.parameters)
    values = {
        name: float(minuit.values[name])
        for name in parameter_names
    }
    errors = {
        name: float(minuit.errors[name])
        for name in parameter_names
    }

    covariance = None
    correlation = None
    if minuit.covariance is not None:
        covariance = [
            [
                float(minuit.covariance[name_i, name_j])
                for name_j in parameter_names
            ]
            for name_i in parameter_names
        ]
        diagonal = np.sqrt(
            np.maximum(np.diag(np.asarray(covariance)), 0.0)
        )
        matrix = np.asarray(covariance, dtype=np.float64)
        denominator = np.outer(diagonal, diagonal)
        correlation_array = np.divide(
            matrix,
            denominator,
            out=np.zeros_like(matrix),
            where=denominator > 0.0,
        )
        correlation = correlation_array.tolist()
    # endif

    fmin = minuit.fmin
    return {
        "parameter_order": parameter_names,
        "values": values,
        "errors": errors,
        "covariance": covariance,
        "correlation": correlation,
        "minimum_nll": float(minuit.fval),
        "valid": bool(fmin.is_valid),
        "accurate_covariance": bool(fmin.has_accurate_covar),
        "positive_definite_covariance": bool(fmin.has_posdef_covar),
        "parameters_at_limit": bool(fmin.has_parameters_at_limit),
        "edm": float(fmin.edm),
        "edm_goal": float(fmin.edm_goal),
        "nfcn": int(fmin.nfcn),
    }


def fit_one_variant(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    bin_number: int,
    variant: str,
    active_periods: tuple[str, ...] = PERIODS,
    transverse_scales: Mapping[str, float] | None = None,
    initial_values: Mapping[str, float] | None = None,
    fixed_physics_parameters: Mapping[str, float] | None = None,
    target_sign_filter: int | None = None,
    run_ranges_filter: tuple[tuple[int, int], ...] | None = None,
    beam_polarization_scales: Mapping[str, float] | None = None,
    double_spin_products_by_run: Mapping[int, float] | None = None,
    fast_profile: bool = False,
) -> dict[str, Any]:
    nll, metadata = make_bin_nll(
        events,
        run_states,
        dilution_records,
        bin_number,
        variant,
        active_periods=active_periods,
        transverse_scales=transverse_scales,
        target_sign_filter=target_sign_filter,
        run_ranges_filter=run_ranges_filter,
        beam_polarization_scales=beam_polarization_scales,
        double_spin_products_by_run=double_spin_products_by_run,
    )

    initial = dict(PARAMETER_INITIAL_VALUES)
    if initial_values is not None:
        for name in PHYSICS_PARAMETERS:
            if name in initial_values:
                initial[name] = float(initial_values[name])
            # endif
        # endfor
    # endif
    effective_fixed_physics_parameters = dict(
        fixed_physics_parameters or {}
    )

    for name, value in effective_fixed_physics_parameters.items():
        if name not in PHYSICS_PARAMETERS:
            raise RuntimeError(
                f"Unknown fixed physics parameter {name!r}."
            )
        # endif
        initial[name] = float(value)
    # endfor

    initial.update(
        {
            f"f_{period}": dilution_records[(period, bin_number)].value
            for period in PERIODS
        }
    )

    def configured_minuit(start_values: Mapping[str, float]) -> Minuit:
        candidate = Minuit(nll, **dict(start_values))
        candidate.errordef = Minuit.LIKELIHOOD
        candidate.print_level = 0
        candidate.strategy = 0 if fast_profile else 1

        for parameter_name in PHYSICS_PARAMETERS:
            candidate.limits[parameter_name] = PARAMETER_LIMITS[
                parameter_name
            ]
            if (
                parameter_name in effective_fixed_physics_parameters
            ):
                candidate.values[parameter_name] = float(
                    effective_fixed_physics_parameters[parameter_name]
                )
                candidate.fixed[parameter_name] = True
            # endif
        # endfor
        for period in PERIODS:
            nuisance_name = f"f_{period}"
            central = dilution_records[(period, bin_number)].value
            sigma = dilution_records[(period, bin_number)].stat_uncertainty
            width = max(8.0 * sigma, 0.20 * central, 0.02)
            candidate.limits[nuisance_name] = (
                max(1.0e-6, central - width),
                central + width,
            )
            if period not in active_periods or sigma == 0.0:
                candidate.fixed[nuisance_name] = True
            # endif
        # endfor
        return candidate

    start_candidates: list[dict[str, float]] = [dict(initial)]
    zero_start = dict(initial)
    zero_start.update({name: 0.0 for name in PHYSICS_PARAMETERS})
    start_candidates.append(zero_start)

    # A deterministic midpoint start is useful for difficult low-statistics
    # bins without introducing run-to-run randomness.
    midpoint_start = dict(initial)
    midpoint_start.update(
        {
            name: 0.5 * float(initial[name])
            for name in PHYSICS_PARAMETERS
        }
    )
    start_candidates.append(midpoint_start)

    sign_reversed = dict(initial)
    for name in PHYSICS_PARAMETERS:
        if (
            name not in effective_fixed_physics_parameters
        ):
            sign_reversed[name] = -float(initial[name])
        # endif
    # endfor
    start_candidates.append(sign_reversed)

    bounded_offset = dict(initial)
    offset_pattern = {
        "u1": 0.20,
        "u2": -0.20,
        "lu1": 0.05,
        "ul1": -0.05,
        "ul2": 0.05,
        "ll0": 0.20,
        "ll1": -0.20,
    }
    for name, offset in offset_pattern.items():
        if (
            name not in effective_fixed_physics_parameters
        ):
            low, high = PARAMETER_LIMITS[name]
            bounded_offset[name] = float(
                np.clip(
                    float(initial[name]) + offset,
                    low + 1.0e-6,
                    high - 1.0e-6,
                )
            )
        # endif
    # endfor
    start_candidates.append(bounded_offset)

    # Profile scans start from the preceding constrained minimum.  Repeating the
    # full five-start production recovery and HESSE at every scan point is both
    # unnecessary and extremely expensive; only the profiled NLL is required.
    if fast_profile:
        start_candidates = start_candidates[:1]
    # endif

    attempted: list[Minuit] = []
    for start_values in start_candidates:
        candidate = configured_minuit(start_values)
        candidate.migrad(ncall=12000 if fast_profile else 50000)
        if not candidate.fmin.is_valid:
            if fast_profile:
                # One bounded recovery attempt.  Do not allow a pathological
                # profile point to consume hundreds of thousands of calls.
                candidate.strategy = 1
                candidate.simplex(ncall=8000)
                candidate.migrad(ncall=20000)
            else:
                candidate.simplex(ncall=50000)
                candidate.strategy = 2
                candidate.migrad(ncall=120000)
            # endif
        # endif
        if not fast_profile:
            candidate.hesse()
        # endif
        attempted.append(candidate)
    # endfor

    def candidate_quality(candidate: Minuit) -> tuple[int, float, float]:
        fmin = candidate.fmin
        trustworthy = (
            fmin.is_valid
            and fmin.has_accurate_covar
            and fmin.has_posdef_covar
            and not fmin.has_parameters_at_limit
            and math.isfinite(float(candidate.fval))
        )
        valid = fmin.is_valid and math.isfinite(float(candidate.fval))
        rank = 0 if trustworthy else (1 if valid else 2)
        edm = (
            float(fmin.edm)
            if math.isfinite(float(fmin.edm))
            else math.inf
        )
        return rank, float(candidate.fval), edm

    minuit = min(attempted, key=candidate_quality)

    result = minuit_result_payload(minuit)
    # Evaluation-only diagnostics: these do not alter Minuit, the likelihood,
    # starting values, limits, or the selected minimum.
    result["diagnostics"] = {
        "initial_values": {name: float(value) for name, value in initial.items()},
        "parameter_limits": {
            name: [float(PARAMETER_LIMITS[name][0]), float(PARAMETER_LIMITS[name][1])]
            for name in PHYSICS_PARAMETERS
        },
        "nll_probes_at_selected_minimum": _finite_nll_probes(
            nll, {name: float(minuit.values[name]) for name in minuit.parameters}
        ),
        "fmin": {
            "is_valid": bool(minuit.fmin.is_valid),
            "has_accurate_covar": bool(minuit.fmin.has_accurate_covar),
            "has_posdef_covar": bool(minuit.fmin.has_posdef_covar),
            "has_parameters_at_limit": bool(minuit.fmin.has_parameters_at_limit),
            "has_reached_call_limit": bool(minuit.fmin.has_reached_call_limit),
            "hesse_failed": bool(minuit.fmin.hesse_failed),
            "nfcn": int(minuit.fmin.nfcn),
            "edm": float(minuit.fmin.edm),
            "edm_goal": float(minuit.fmin.edm_goal),
        },
    }
    result.update(
        {
            "bin_number": bin_number,
            "variant": variant,
            "active_periods": list(active_periods),
            "target_sign_filter": target_sign_filter,
            "transverse_scales": (
                dict(transverse_scales)
                if transverse_scales is not None else None
            ),
            "fixed_physics_parameters": (
                dict(effective_fixed_physics_parameters)
                if effective_fixed_physics_parameters else None
            ),
            "metadata": metadata,
        }
    )
    return result


_WORKER_EVENTS: dict[str, np.ndarray] | None = None
_WORKER_RUN_STATES: dict[str, dict[str, np.ndarray]] | None = None
_WORKER_DILUTION_RECORDS: dict[tuple[str, int], DilutionRecord] | None = None


def initialize_fit_worker(
    cache_path_text: str,
    run_state_payload: dict[str, dict[str, list[float] | list[int]]],
    dilution_payload: dict[str, dict[str, dict[str, float | int]]],
) -> None:
    """Load immutable fit inputs once per worker process."""
    global _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS

    _WORKER_EVENTS = load_event_cache(Path(cache_path_text))
    _WORKER_RUN_STATES = {
        period: {
            key: np.asarray(values)
            for key, values in state.items()
        }
        for period, state in run_state_payload.items()
    }
    _WORKER_DILUTION_RECORDS = {
        (period, int(bin_text)): DilutionRecord(
            period=period,
            bin_number=int(bin_text),
            x_index=int(record["x_index"]),
            t_index=int(record["t_index"]),
            value=float(record["value"]),
            stat_uncertainty=float(record["stat_uncertainty"]),
        )
        for period, period_payload in dilution_payload.items()
        for bin_text, record in period_payload.items()
    }


def fit_bin_worker(
    bin_number: int,
    include_target_axis_study: bool = True,
    include_period_diagnostics: bool = True,
) -> dict[str, Any]:
    if (
        _WORKER_EVENTS is None
        or _WORKER_RUN_STATES is None
        or _WORKER_DILUTION_RECORDS is None
    ):
        raise RuntimeError("Fit worker was not initialized.")
    # endif

    events = _WORKER_EVENTS
    run_states = _WORKER_RUN_STATES
    dilution_records = _WORKER_DILUTION_RECORDS
    worker_start = time.perf_counter()

    def report(stage: str, detail: str = "") -> None:
        elapsed = time.perf_counter() - worker_start
        suffix = f" | {detail}" if detail else ""
        print(
            f"[worker bin {bin_number:02d}] {stage} "
            f"(elapsed {elapsed:8.1f} s){suffix}",
            flush=True,
        )
    # enddef

    report("START")

    report("nominal simultaneous fit: START")
    nominal = fit_one_variant(
        events,
        run_states,
        dilution_records,
        bin_number,
        "nominal",
    )
    report(
        "nominal simultaneous fit: DONE",
        f"valid={nominal['valid']}; EDM={nominal['edm']:.3e}",
    )
    variants = {"nominal": nominal}
    photon_axis_projection_shift: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }
    external_data_shift: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }
    full_three_fit_spread: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }
    systematics: dict[str, float | None] = {
        parameter: None for parameter in PHYSICS_PARAMETERS
    }

    if include_target_axis_study:
        report("target-axis photon_axis_projection: START")
        photon_axis_projection = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "photon_axis_projection",
            initial_values=nominal["values"],
        )
        report(
            "target-axis photon_axis_projection: DONE",
            f"valid={photon_axis_projection['valid']}; EDM={photon_axis_projection['edm']:.3e}",
        )
        report("target-axis external_data_informed: START")
        external_data_informed = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "external_data_informed",
            initial_values=nominal["values"],
        )
        report(
            "target-axis external_data_informed: DONE",
            f"valid={external_data_informed['valid']}; "
            f"EDM={external_data_informed['edm']:.3e}",
        )
        variants.update({
            "photon_axis_projection": photon_axis_projection,
            "external_data_informed": external_data_informed,
        })
        for parameter in PHYSICS_PARAMETERS:
            nominal_value = nominal["values"][parameter]
            photon_axis_projection_value = photon_axis_projection["values"][parameter]
            external_value = external_data_informed["values"][parameter]
            external_data_shift[parameter] = abs(
                external_value - nominal_value
            )

            photon_axis_projection_shift[parameter] = abs(
                photon_axis_projection_value - nominal_value
            )
            systematics[parameter] = max(
                photon_axis_projection_shift[parameter],
                external_data_shift[parameter],
            )
            full_three_fit_spread[parameter] = (
                max(nominal_value, photon_axis_projection_value, external_value)
                - min(nominal_value, photon_axis_projection_value, external_value)
            )
        # endfor
    # endif

    # Period-only and constrained fits are diagnostics only.  They are
    # intentionally evaluated for the production nominal extraction, where
    # they feed the final period-stability products, but are skipped for
    # systematic-variation samples.  This does not alter the simultaneous
    # extraction used for any quoted result or systematic variation.
    period_fits: dict[str, dict[str, Any]] = {}
    period_constraint_fits: dict[str, dict[str, dict[str, Any]]] = {
        "fix_u1": {}, "fix_u2": {}, "fix_u1_u2": {}
    }
    period_consistency: dict[str, dict[str, float]] = {}
    if include_period_diagnostics:
        # Independent nominal fits for period-consistency plots.  These do not
        # enter the quoted combined result or its target-axis systematic.
        period_fits = {}
        for period in PERIODS:
            report(f"period-only {period}: START")
            period_fits[period] = fit_one_variant(
                events,
                run_states,
                dilution_records,
                bin_number,
                "nominal",
                active_periods=(period,),
                initial_values=nominal["values"],
            )
            report(
                f"period-only {period}: DONE",
                f"valid={period_fits[period]['valid']}; "
                f"EDM={period_fits[period]['edm']:.3e}",
            )
        # endfor

        period_constraint_fits: dict[
            str,
            dict[str, dict[str, Any]],
        ] = {
            "fix_u1": {},
            "fix_u2": {},
            "fix_u1_u2": {},
        }
        for period in PERIODS:
            report(f"period constraint fix_u1 {period}: START")
            period_constraint_fits["fix_u1"][period] = fit_one_variant(
                events,
                run_states,
                dilution_records,
                bin_number,
                "nominal",
                active_periods=(period,),
                initial_values=period_fits[period]["values"],
                fixed_physics_parameters={
                    "u1": nominal["values"]["u1"],
                },
            )
            report(
                f"period constraint fix_u1 {period}: DONE",
                f"valid={period_constraint_fits['fix_u1'][period]['valid']}; "
                f"EDM={period_constraint_fits['fix_u1'][period]['edm']:.3e}",
            )
            report(f"period constraint fix_u2 {period}: START")
            period_constraint_fits["fix_u2"][period] = fit_one_variant(
                events,
                run_states,
                dilution_records,
                bin_number,
                "nominal",
                active_periods=(period,),
                initial_values=period_fits[period]["values"],
                fixed_physics_parameters={
                    "u2": nominal["values"]["u2"],
                },
            )
            report(
                f"period constraint fix_u2 {period}: DONE",
                f"valid={period_constraint_fits['fix_u2'][period]['valid']}; "
                f"EDM={period_constraint_fits['fix_u2'][period]['edm']:.3e}",
            )
            report(f"period constraint fix_u1_u2 {period}: START")
            period_constraint_fits["fix_u1_u2"][period] = fit_one_variant(
                events,
                run_states,
                dilution_records,
                bin_number,
                "nominal",
                active_periods=(period,),
                initial_values=period_fits[period]["values"],
                fixed_physics_parameters={
                    "u1": nominal["values"]["u1"],
                    "u2": nominal["values"]["u2"],
                },
            )
            report(
                f"period constraint fix_u1_u2 {period}: DONE",
                f"valid={period_constraint_fits['fix_u1_u2'][period]['valid']}; "
                f"EDM={period_constraint_fits['fix_u1_u2'][period]['edm']:.3e}",
            )
        # endfor

        # Quantify each period's tension with the simultaneous solution.  For each
        # period, compare its nominal NLL at the combined best-fit point with the
        # independently minimized period-only NLL.  The difference is nonnegative
        # up to minimizer precision and is a compact period-consistency diagnostic.
        period_consistency: dict[str, dict[str, float]] = {}
        for period in PERIODS:
            report(f"period consistency {period}: START")
            period_nll, _ = make_bin_nll(
                events,
                run_states,
                dilution_records,
                bin_number,
                "nominal",
                active_periods=(period,),
            )
            combined_nll_for_period = float(
                period_nll(
                    **{
                        name: nominal["values"][name]
                        for name in (
                            *PHYSICS_PARAMETERS,
                            "f_su22",
                            "f_fa22",
                            "f_sp23",
                        )
                    }
                )
            )
            period_minimum_nll = float(period_fits[period]["minimum_nll"])
            period_consistency[period] = {
                "nll_at_combined_solution": combined_nll_for_period,
                "period_only_minimum_nll": period_minimum_nll,
                "delta_nll": max(
                    0.0,
                    combined_nll_for_period - period_minimum_nll,
                ),
            }
            report(f"period consistency {period}: DONE")
        # endfor
    # endif

    report("DONE; returning result to parent")
    return {
        "bin_number": bin_number,
        "variants": variants,
        "period_fits": period_fits,
        "period_constraint_fits": period_constraint_fits,
        "period_consistency": period_consistency,
        "target_axis_study_envelope": systematics,
        "photon_axis_projection_shift": photon_axis_projection_shift,
        "external_data_shift": external_data_shift,
        "full_three_fit_spread": full_three_fit_spread,
        "external_transverse_inputs": (
            EXTERNAL_TRANSVERSE_INPUTS if include_target_axis_study else None
        ),
        "target_axis_study_performed": include_target_axis_study,
        "period_diagnostics_performed": include_period_diagnostics,
    }


# =============================================================================
# Outputs
# =============================================================================



def _nll_vector(nll: Any, values: Mapping[str, float]) -> float:
    """Evaluate the existing NLL without changing any fit configuration."""
    order = (*PHYSICS_PARAMETERS, "f_su22", "f_fa22", "f_sp23")
    return float(nll(*[float(values[name]) for name in order]))


def _finite_nll_probes(
    nll: Any,
    values: Mapping[str, float],
    *,
    fractional_step: float = 0.02,
    absolute_step: float = 0.01,
) -> dict[str, Any]:
    """Probe local NLL sensitivity using evaluations only (no minimization)."""
    center = _nll_vector(nll, values)
    probes: dict[str, Any] = {}
    for name in PHYSICS_PARAMETERS:
        low, high = PARAMETER_LIMITS[name]
        step = max(absolute_step, fractional_step * (high - low))
        minus_values = dict(values)
        plus_values = dict(values)
        minus_values[name] = max(low + 1.0e-8, float(values[name]) - step)
        plus_values[name] = min(high - 1.0e-8, float(values[name]) + step)
        minus_nll = _nll_vector(nll, minus_values)
        plus_nll = _nll_vector(nll, plus_values)
        probes[name] = {
            "step": float(step),
            "minus_value": float(minus_values[name]),
            "plus_value": float(plus_values[name]),
            "delta_nll_minus": float(minus_nll - center),
            "delta_nll_plus": float(plus_nll - center),
        }
    # endfor
    return {"center_nll": center, "parameters": probes}


def _period_state_counts(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    bin_number: int,
    period: str,
) -> dict[str, Any]:
    """Count observed helicity/target-polarization states for diagnostics."""
    mask = (
        (events["bin_number"] == bin_number)
        & (events["period_index"] == PERIOD_INDEX[period])
    )
    runs = events["runnum"][mask].astype(np.int32, copy=False)
    helicity = events["helicity"][mask].astype(np.int8, copy=False)
    lookup = {
        int(run): float(pt)
        for run, pt in zip(
            run_states[period]["run"], run_states[period]["pt"]
        )
    }
    pt = np.fromiter((lookup[int(run)] for run in runs), dtype=np.float64,
                     count=runs.size)
    target_sign = np.sign(pt).astype(np.int8, copy=False)
    counts = {}
    for h in (-1, 1):
        for t in (-1, 0, 1):
            counts[f"h{h:+d}_t{t:+d}"] = int(
                np.count_nonzero((helicity == h) & (target_sign == t))
            )
        # endfor
    # endfor
    return {
        "events": int(runs.size),
        "unique_runs": int(np.unique(runs).size),
        "helicity_plus": int(np.count_nonzero(helicity > 0)),
        "helicity_minus": int(np.count_nonzero(helicity < 0)),
        "target_plus": int(np.count_nonzero(target_sign > 0)),
        "target_minus": int(np.count_nonzero(target_sign < 0)),
        "state_counts": counts,
    }


def period_preflight_worker(task: dict[str, Any]) -> dict[str, Any]:
    """Diagnose one period/bin likelihood without running Minuit."""
    if (_WORKER_EVENTS is None or _WORKER_RUN_STATES is None
            or _WORKER_DILUTION_RECORDS is None):
        raise RuntimeError("Fit worker was not initialized.")
    bin_number = int(task["bin_number"])
    period = str(task["period"])
    nominal = task["nominal"]
    # This worker diagnoses a single run period at a time.  The period-only
    # fit workers pass ``active_periods=(period,)`` explicitly; mirror that
    # here rather than referring to an undefined outer-scope variable.
    active_periods = (period,)
    # make_bin_nll constructs the likelihood only; parameter fixing is a
    # Minuit configuration handled by fit_one_variant.  The preflight does not
    # run Minuit, so enforce the zero-UU diagnostic model by setting the probe
    # values below rather than passing an unsupported keyword here.
    nll, metadata = make_bin_nll(
        _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
        bin_number, "nominal", active_periods=active_periods,
    )
    values = dict(PARAMETER_INITIAL_VALUES)
    values.update({name: float(nominal["values"][name])
                   for name in PHYSICS_PARAMETERS})
    # Match the actual period-stability fit model exactly.
    values["u1"] = 0.0
    values["u2"] = 0.0
    values.update({
        f"f_{p}": float(_WORKER_DILUTION_RECORDS[(p, bin_number)].value)
        for p in PERIODS
    })
    probes = _finite_nll_probes(nll, values)
    response_by_parameter = {
        name: max(abs(item["delta_nll_minus"]), abs(item["delta_nll_plus"]))
        for name, item in probes["parameters"].items()
    }
    weakest = sorted(response_by_parameter, key=response_by_parameter.get)[:2]
    scans: dict[str, Any] = {}
    for name in weakest:
        low, high = PARAMETER_LIMITS[name]
        grid = np.linspace(low, high, 31)
        nll_values = []
        for value in grid:
            trial = dict(values)
            trial[name] = float(value)
            nll_values.append(_nll_vector(nll, trial))
        # endfor
        finite = np.asarray(nll_values, dtype=np.float64)
        finite_mask = np.isfinite(finite)
        minimum = float(np.min(finite[finite_mask])) if np.any(finite_mask) else math.inf
        scans[name] = {
            "grid": grid.tolist(),
            "delta_nll": [
                float(v - minimum) if math.isfinite(float(v)) and math.isfinite(minimum) else None
                for v in nll_values
            ],
        }
    # endfor
    payload = {
        "bin_number": bin_number,
        "period": period,
        "counts": _period_state_counts(
            _WORKER_EVENTS, _WORKER_RUN_STATES, bin_number, period
        ),
        "start_values": values,
        "nll_probes_at_start": probes,
        "weakest_parameters": weakest,
        "one_dimensional_scans": scans,
        "metadata": metadata,
    }
    flatness = min(response_by_parameter.values())
    print(
        f"[preflight] bin {bin_number:02d} | {period} | "
        f"N={payload['counts']['events']:,} | weakest local response={flatness:.3e}",
        flush=True,
    )
    return payload


def run_period_preflight_stage(
    *, tasks: list[dict[str, Any]], workers: int, cache_path: Path,
    run_state_payload: dict[str, Any], dilution_payload: dict[str, Any],
    output_dir: Path,
) -> list[dict[str, Any]]:
    """Run evaluation-only diagnostics before any period-only minimizations."""
    print("=" * 78, flush=True)
    print(f"[stage period preflight] START | tasks={len(tasks)} | max_workers={workers}", flush=True)
    print("=" * 78, flush=True)
    completed = []
    mp_context = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers, mp_context=mp_context,
        initializer=initialize_fit_worker,
        initargs=(str(cache_path), run_state_payload, dilution_payload),
    ) as executor:
        futures = {executor.submit(period_preflight_worker, task): task for task in tasks}
        for index, future in enumerate(as_completed(futures), start=1):
            completed.append(future.result())
            write_json(
                output_dir / "json" / "checkpoint_03a_period_preflight_partial.json",
                {"schema_version": 1, "completed": sorted(
                    completed, key=lambda x: (x["bin_number"], x["period"])
                )},
            )
            print(f"[stage period preflight] progress {index}/{len(tasks)}", flush=True)
        # endfor
    # endwith
    write_json(
        output_dir / "json" / "checkpoint_03a_period_preflight.json",
        {"schema_version": 1, "completed": sorted(
            completed, key=lambda x: (x["bin_number"], x["period"])
        )},
    )
    return completed


def _task_key(task: Mapping[str, Any]) -> str:
    # Treat an omitted optional task field and an explicit JSON null as the
    # same key. Worker results store absent period/constraint values as None,
    # while the submitted task dictionaries may omit those keys entirely.
    # Without this normalization, completed checkpoint tasks are loaded but
    # fail to match the corresponding pending tasks on resume.
    fields = ("kind", "bin_number", "period", "constraint")
    return "|".join(
        "" if task.get(key) is None else str(task.get(key))
        for key in fields
    )

def fit_stage_worker(task: dict[str, Any]) -> dict[str, Any]:
    """Execute exactly one existing fit task for staged production.

    This function changes only scheduling/granularity.  Every fit delegates to
    fit_one_variant with the same arguments used by the original bundled
    fit_bin_worker implementation.
    """
    if (
        _WORKER_EVENTS is None
        or _WORKER_RUN_STATES is None
        or _WORKER_DILUTION_RECORDS is None
    ):
        raise RuntimeError("Fit worker was not initialized.")
    # endif

    events = _WORKER_EVENTS
    run_states = _WORKER_RUN_STATES
    dilution_records = _WORKER_DILUTION_RECORDS
    bin_number = int(task["bin_number"])
    kind = str(task["kind"])
    period = task.get("period")
    constraint = task.get("constraint")
    nominal = task.get("nominal")
    period_fit = task.get("period_fit")
    t0 = time.perf_counter()

    label = f"bin {bin_number:02d} | {kind}"
    if period is not None:
        label += f" | {period}"
    # endif
    if constraint is not None:
        label += f" | {constraint}"
    # endif
    print(f"[stage fit START] {label}", flush=True)

    if kind == "nominal":
        fit = fit_one_variant(
            events, run_states, dilution_records, bin_number, "nominal"
        )
    elif kind == "nominal_zero_uu":
        fit = fit_one_variant(
            events, run_states, dilution_records, bin_number, "nominal",
            fixed_physics_parameters={"u1": 0.0, "u2": 0.0},
        )
    elif kind in {"photon_axis_projection", "external_data_informed"}:
        fit = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            kind,
            initial_values=nominal["values"],
        )
    elif kind in {"period_only", "period_only_zero_uu"}:
        fit = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=nominal["values"],
            fixed_physics_parameters=(
                {"u1": 0.0, "u2": 0.0}
                if kind == "period_only_zero_uu" else None
            ),
        )
    elif kind == "period_constraint":
        if constraint == "fix_u1":
            fixed = {"u1": nominal["values"]["u1"]}
        elif constraint == "fix_u2":
            fixed = {"u2": nominal["values"]["u2"]}
        elif constraint == "fix_u1_u2":
            fixed = {
                "u1": nominal["values"]["u1"],
                "u2": nominal["values"]["u2"],
            }
        else:
            raise ValueError(f"Unknown period constraint: {constraint}")
        # endif
        fit = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=period_fit["values"],
            fixed_physics_parameters=fixed,
        )
    else:
        raise ValueError(f"Unknown staged fit kind: {kind}")
    # endif

    elapsed = time.perf_counter() - t0
    print(
        f"[stage fit DONE ] {label} | {elapsed:.1f} s | "
        f"valid={fit['valid']} | EDM={fit['edm']:.3e}",
        flush=True,
    )
    return {
        "bin_number": bin_number,
        "kind": kind,
        "period": period,
        "constraint": constraint,
        "fit": fit,
        "elapsed_seconds": elapsed,
    }


def run_fit_stage(
    *,
    stage_name: str,
    tasks: list[dict[str, Any]],
    workers: int,
    cache_path: Path,
    run_state_payload: dict[str, Any],
    dilution_payload: dict[str, Any],
    output_dir: Path | None = None,
    checkpoint_tag: str | None = None,
) -> list[dict[str, Any]]:
    """Run one fit stage with incremental task-level checkpointing.

    Completed tasks are written immediately.  If the same stage is restarted,
    those exact completed fit payloads are reused and only missing tasks are
    submitted.  This changes scheduling only, never the fit itself.
    """
    partial_path = None
    completed_by_key: dict[str, dict[str, Any]] = {}
    if output_dir is not None and checkpoint_tag is not None:
        partial_path = output_dir / "json" / f"checkpoint_{checkpoint_tag}_partial.json"
        if partial_path.is_file():
            try:
                payload = json.loads(partial_path.read_text())
                for item in payload.get("completed", []):
                    key = _task_key(item)
                    completed_by_key[key] = item
                # endfor
                print(f"[stage {stage_name}] RESUME: loaded {len(completed_by_key)} completed tasks", flush=True)
            except Exception as exc:
                print(f"[stage {stage_name}] WARNING: could not read partial checkpoint: {exc}", flush=True)
            # endtry
        # endif
    # endif
    pending_tasks = [task for task in tasks if _task_key(task) not in completed_by_key]

    print("=" * 78, flush=True)
    print(
        f"[stage {stage_name}] START | tasks={len(tasks)} | pending={len(pending_tasks)} | "
        f"max_workers={workers}",
        flush=True,
    )
    print("=" * 78, flush=True)
    t0 = time.perf_counter()
    completed: list[dict[str, Any]] = list(completed_by_key.values())
    if not pending_tasks:
        print(f"[stage {stage_name}] all tasks already checkpointed.", flush=True)
        return sorted(completed, key=lambda x: (_task_key(x)))
    # endif
    mp_context = mp.get_context("spawn")
    executor = ProcessPoolExecutor(
        max_workers=workers,
        mp_context=mp_context,
        initializer=initialize_fit_worker,
        initargs=(str(cache_path), run_state_payload, dilution_payload),
    )
    try:
        futures = {
            executor.submit(fit_stage_worker, task): task
            for task in pending_tasks
        }
        total = len(tasks)
        already = len(completed_by_key)
        for index, future in enumerate(as_completed(futures), start=already + 1):
            task = futures[future]
            try:
                item = future.result()
            except Exception:
                print(
                    f"[stage {stage_name}] ERROR in bin "
                    f"{int(task['bin_number']):02d}, kind={task['kind']}, "
                    f"period={task.get('period')}, "
                    f"constraint={task.get('constraint')}",
                    flush=True,
                )
                raise
            # endtry
            completed.append(item)
            if partial_path is not None:
                write_json(
                    partial_path,
                    {
                        "schema_version": 1,
                        "stage": stage_name,
                        "updated_utc": datetime.now(timezone.utc).isoformat(),
                        "completed": sorted(completed, key=_task_key),
                    },
                )
            # endif
            print(
                f"[stage {stage_name}] progress {index}/{total} completed",
                flush=True,
            )
        # endfor
        print(
            f"[stage {stage_name}] all futures returned; shutting down workers...",
            flush=True,
        )
    finally:
        executor.shutdown(wait=True, cancel_futures=False)
    # endtry
    elapsed = time.perf_counter() - t0
    print(
        f"[stage {stage_name}] DONE | {elapsed:.1f} s",
        flush=True,
    )
    return completed


def initialize_result_shell(
    bin_number: int,
    nominal: dict[str, Any],
    include_target_axis_study: bool,
    include_period_diagnostics: bool,
) -> dict[str, Any]:
    empty = {parameter: None for parameter in PHYSICS_PARAMETERS}
    return {
        "bin_number": bin_number,
        "variants": {"nominal": nominal},
        "period_fits": {},
        "period_constraint_fits": {
            "fix_u1": {}, "fix_u2": {}, "fix_u1_u2": {}
        },
        "period_consistency": {},
        "target_axis_study_envelope": dict(empty),
        "photon_axis_projection_shift": dict(empty),
        "external_data_shift": dict(empty),
        "full_three_fit_spread": dict(empty),
        "external_transverse_inputs": (
            EXTERNAL_TRANSVERSE_INPUTS if include_target_axis_study else None
        ),
        "target_axis_study_performed": include_target_axis_study,
        "period_diagnostics_performed": include_period_diagnostics,
    }


def finish_target_axis_study_envelopes(result: dict[str, Any]) -> None:
    nominal = result["variants"]["nominal"]
    photon_axis_projection = result["variants"]["photon_axis_projection"]
    external_data_informed = result["variants"]["external_data_informed"]
    for parameter in PHYSICS_PARAMETERS:
        nominal_value = nominal["values"][parameter]
        photon_axis_projection_value = photon_axis_projection["values"][parameter]
        external_value = external_data_informed["values"][parameter]
        external_shift = abs(external_value - nominal_value)
        result["external_data_shift"][parameter] = external_shift
        projection_shift = abs(photon_axis_projection_value - nominal_value)
        result["photon_axis_projection_shift"][parameter] = projection_shift
        result["target_axis_study_envelope"][parameter] = max(
            projection_shift, external_shift
        )
        result["full_three_fit_spread"][parameter] = (
            max(nominal_value, photon_axis_projection_value, external_value)
            - min(nominal_value, photon_axis_projection_value, external_value)
        )
    # endfor


def finish_period_consistency(
    result: dict[str, Any],
    events: dict[str, np.ndarray],
    run_states: dict[str, dict[str, np.ndarray]],
    dilution_records: dict[tuple[str, int], DilutionRecord],
) -> None:
    """Reproduce the original non-minimizing period-consistency calculation."""
    bin_number = int(result["bin_number"])
    nominal = result["variants"]["nominal"]
    for period in PERIODS:
        period_nll, _ = make_bin_nll(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
        )
        combined_nll_for_period = float(
            period_nll(
                **{
                    name: nominal["values"][name]
                    for name in (
                        *PHYSICS_PARAMETERS,
                        "f_su22",
                        "f_fa22",
                        "f_sp23",
                    )
                }
            )
        )
        period_minimum_nll = float(
            result["period_fits"][period]["minimum_nll"]
        )
        result["period_consistency"][period] = {
            "nll_at_combined_solution": combined_nll_for_period,
            "period_only_minimum_nll": period_minimum_nll,
            "delta_nll": max(
                0.0, combined_nll_for_period - period_minimum_nll
            ),
        }
    # endfor


def write_stage_checkpoint(
    output_dir: Path,
    stage_name: str,
    results: list[dict[str, Any]],
) -> Path:
    path = output_dir / "json" / f"checkpoint_{stage_name}.json"
    write_json(
        path,
        {
            "schema_version": 1,
            "stage": stage_name,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "results": sorted(results, key=lambda item: item["bin_number"]),
        },
    )
    print(f"[checkpoint] wrote {path}", flush=True)
    return path

def flatten_fit_results(
    results: list[dict[str, Any]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for result in sorted(results, key=lambda item: item["bin_number"]):
        nominal = result["variants"]["nominal"]
        metadata = nominal["metadata"]
        row: dict[str, Any] = {
            "bin_number": result["bin_number"],
            "x_index": (result["bin_number"] - 1)
            // len(MINUS_TPRIME_BINS_GEV2),
            "t_index": (result["bin_number"] - 1)
            % len(MINUS_TPRIME_BINS_GEV2),
            "min_xB": metadata["min_xB"],
            "median_xB": metadata["median_xB"],
            "mean_xB": metadata["mean_xB"],
            "max_xB": metadata["max_xB"],
            "min_Q2_gev2": metadata["min_Q2_gev2"],
            "median_Q2_gev2": metadata["median_Q2_gev2"],
            "mean_Q2_gev2": metadata["mean_Q2_gev2"],
            "max_Q2_gev2": metadata["max_Q2_gev2"],
            "min_W_gev": metadata["min_W_gev"],
            "median_W_gev": metadata["median_W_gev"],
            "mean_W_gev": metadata["mean_W_gev"],
            "max_W_gev": metadata["max_W_gev"],
            "min_minus_t_gev2": metadata["min_minus_t_gev2"],
            "median_minus_t_gev2": metadata["median_minus_t_gev2"],
            "mean_minus_t_gev2": metadata["mean_minus_t_gev2"],
            "max_minus_t_gev2": metadata["max_minus_t_gev2"],
            "min_minus_tprime_gev2": metadata["min_minus_tprime_gev2"],
            "median_minus_tprime_gev2": metadata["median_minus_tprime_gev2"],
            "mean_minus_tprime_gev2": metadata[
                "mean_minus_tprime_gev2"
            ],
            "max_minus_tprime_gev2": metadata["max_minus_tprime_gev2"],
            "min_epsilon": metadata["min_epsilon"],
            "median_epsilon": metadata["median_epsilon"],
            "mean_epsilon": metadata["mean_epsilon"],
            "max_epsilon": metadata["max_epsilon"],
            "min_DepA": metadata["min_DepA"],
            "median_DepA": metadata["median_DepA"],
            "mean_DepA": metadata["mean_DepA"],
            "max_DepA": metadata["max_DepA"],
            "min_DepB": metadata["min_DepB"],
            "median_DepB": metadata["median_DepB"],
            "mean_DepB": metadata["mean_DepB"],
            "max_DepB": metadata["max_DepB"],
            "min_DepC": metadata["min_DepC"],
            "median_DepC": metadata["median_DepC"],
            "mean_DepC": metadata["mean_DepC"],
            "max_DepC": metadata["max_DepC"],
            "min_DepV": metadata["min_DepV"],
            "median_DepV": metadata["median_DepV"],
            "mean_DepV": metadata["mean_DepV"],
            "max_DepV": metadata["max_DepV"],
            "min_DepW": metadata["min_DepW"],
            "median_DepW": metadata["median_DepW"],
            "mean_DepW": metadata["mean_DepW"],
            "max_DepW": metadata["max_DepW"],
            "number_of_events": metadata["number_of_events"],
            "mean_sin_theta_gamma": metadata["mean_sin_theta_gamma"],
            "rms_sin_theta_gamma": metadata["rms_sin_theta_gamma"],
            "mean_cos_theta_gamma": metadata["mean_cos_theta_gamma"],
            "rms_cos_theta_gamma": metadata["rms_cos_theta_gamma"],
            "nominal_fit_valid": nominal["valid"],
            "nominal_minimum_nll": nominal["minimum_nll"],
            "nominal_edm": nominal["edm"],
        }

        for period in PERIODS:
            row[f"events_{period}"] = metadata["events_by_period"][period]
            row[f"dilution_{period}"] = metadata[
                "dilution_central"
            ][period]
            row[f"dilution_stat_{period}"] = metadata[
                "dilution_stat_uncertainty"
            ][period]
            row[f"fitted_dilution_{period}"] = nominal[
                "values"
            ][f"f_{period}"]
            row[f"fitted_dilution_error_{period}"] = nominal[
                "errors"
            ][f"f_{period}"]

            period_fit = result["period_fits"].get(period)
            if period_fit is None:
                row[f"period_fit_valid_{period}"] = np.nan
                row[f"period_fit_accurate_covariance_{period}"] = np.nan
                row[f"period_fit_positive_definite_covariance_{period}"] = np.nan
                row[f"period_fit_parameters_at_limit_{period}"] = np.nan
                row[f"period_fit_nll_{period}"] = np.nan
                row[f"period_fit_edm_{period}"] = np.nan
                row[f"period_delta_nll_{period}"] = np.nan
                for constraint_name in ("fix_u1", "fix_u2", "fix_u1_u2"):
                    prefix = f"{constraint_name}_period_fit"
                    row[f"{prefix}_valid_{period}"] = np.nan
                    row[f"{prefix}_accurate_covariance_{period}"] = np.nan
                    row[f"{prefix}_positive_definite_covariance_{period}"] = np.nan
                    row[f"{prefix}_parameters_at_limit_{period}"] = np.nan
                    row[f"{prefix}_edm_{period}"] = np.nan
                # endfor
                continue
            # endif
            row[f"period_fit_valid_{period}"] = period_fit["valid"]
            row[f"period_fit_accurate_covariance_{period}"] = period_fit[
                "accurate_covariance"
            ]
            row[f"period_fit_positive_definite_covariance_{period}"] = (
                period_fit["positive_definite_covariance"]
            )
            row[f"period_fit_parameters_at_limit_{period}"] = period_fit[
                "parameters_at_limit"
            ]
            row[f"period_fit_nll_{period}"] = period_fit["minimum_nll"]
            row[f"period_fit_edm_{period}"] = period_fit["edm"]
            row[f"period_delta_nll_{period}"] = result[
                "period_consistency"
            ][period]["delta_nll"]

            for constraint_name in ("fix_u1", "fix_u2", "fix_u1_u2"):
                constrained_fit = result["period_constraint_fits"][
                    constraint_name
                ][period]
                prefix = f"{constraint_name}_period_fit"
                row[f"{prefix}_valid_{period}"] = constrained_fit["valid"]
                row[f"{prefix}_accurate_covariance_{period}"] = (
                    constrained_fit["accurate_covariance"]
                )
                row[
                    f"{prefix}_positive_definite_covariance_{period}"
                ] = constrained_fit["positive_definite_covariance"]
                row[f"{prefix}_parameters_at_limit_{period}"] = (
                    constrained_fit["parameters_at_limit"]
                )
                row[f"{prefix}_edm_{period}"] = constrained_fit["edm"]
            # endfor
        # endfor

        ll0_value = float(nominal["values"]["ll0"])
        ll0_error = float(nominal["errors"]["ll0"])
        row["ll0_central_outside_unit_interval"] = abs(ll0_value) > 1.0
        row["ll0_stat_interval_crosses_unit_boundary"] = (
            ll0_value - ll0_error < -1.0
            or ll0_value + ll0_error > 1.0
        )
        row["ll0_distance_to_nearest_unit_boundary"] = (
            1.0 - abs(ll0_value)
        )

        for parameter in PHYSICS_PARAMETERS:
            row[parameter] = nominal["values"][parameter]
            row[f"{parameter}_stat"] = nominal["errors"][parameter]
            row[f"{parameter}_projection_sys"] = result[
                "photon_axis_projection_shift"
            ][parameter]
            row[f"{parameter}_external_data_sys"] = result[
                "external_data_shift"
            ][parameter]
            row[f"{parameter}_three_fit_full_spread"] = result[
                "full_three_fit_spread"
            ][parameter]
            row[f"{parameter}_target_axis_study_envelope"] = result[
                "target_axis_study_envelope"
            ][parameter]
            for variant in FIT_VARIANTS:
                variant_fit = result["variants"].get(variant)
                row[f"{parameter}_{variant}"] = (
                    variant_fit["values"][parameter]
                    if variant_fit is not None
                    else np.nan
                )
                row[f"{parameter}_stat_{variant}"] = (
                    variant_fit["errors"][parameter]
                    if variant_fit is not None
                    else np.nan
                )
            # endfor
            for period in PERIODS:
                period_fit = result["period_fits"].get(period)
                if period_fit is None:
                    row[f"{parameter}_{period}"] = np.nan
                    row[f"{parameter}_stat_{period}"] = np.nan
                    for constraint_name in ("fix_u1", "fix_u2", "fix_u1_u2"):
                        row[f"{parameter}_{constraint_name}_{period}"] = np.nan
                        row[f"{parameter}_{constraint_name}_stat_{period}"] = np.nan
                    # endfor
                    continue
                # endif
                row[f"{parameter}_{period}"] = period_fit[
                    "values"
                ][parameter]
                row[f"{parameter}_stat_{period}"] = period_fit[
                    "errors"
                ][parameter]

                for constraint_name in (
                    "fix_u1",
                    "fix_u2",
                    "fix_u1_u2",
                ):
                    constrained_fit = result["period_constraint_fits"][
                        constraint_name
                    ][period]
                    row[
                        f"{parameter}_{constraint_name}_{period}"
                    ] = constrained_fit["values"][parameter]
                    row[
                        f"{parameter}_{constraint_name}_stat_{period}"
                    ] = constrained_fit["errors"][parameter]
                # endfor
            # endfor
        # endfor
        rows.append(row)
    # endfor
    return pd.DataFrame(rows)


def apply_parameter_y_limits(ax: plt.Axes, parameter: str) -> None:
    limits = PARAMETER_Y_LIMITS[parameter]
    if limits is not None:
        ax.set_ylim(*limits)
    # endif


def systematic_comparison_y_limits(parameter: str) -> tuple[float, float]:
    """Return the common ratio/difference scale for a parameter family."""
    if parameter in SINGLE_SPIN_PARAMETERS:
        return SYSTEMATIC_COMPARISON_SINGLE_SPIN_Y_LIMITS
    # endif
    if parameter in UNPOLARIZED_DOUBLE_SPIN_PARAMETERS:
        return SYSTEMATIC_COMPARISON_UNPOLARIZED_DOUBLE_Y_LIMITS
    # endif
    raise KeyError(f"No systematic-comparison axis family for {parameter!r}.")


def systematic_size_y_limits(parameter: str) -> tuple[float, float]:
    """Return the common positive scale for assigned systematic magnitudes."""
    if parameter in SINGLE_SPIN_PARAMETERS:
        return (0.0, 0.2)
    # endif
    if parameter in UNPOLARIZED_DOUBLE_SPIN_PARAMETERS:
        return (0.0, 0.2)
    # endif
    raise KeyError(f"No systematic-size axis family for {parameter!r}.")


def configure_systematic_comparison_axes(
    axes: np.ndarray,
    parameter: str,
) -> None:
    """Apply identical axes to every three-panel systematic canvas."""
    axes[0].set_ylim(*systematic_comparison_y_limits(parameter))
    axes[1].set_ylim(*systematic_size_y_limits(parameter))
    axes[2].set_yscale("log")
    axes[2].set_ylim(*SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS)
    for ax in axes:
        ax.grid(alpha=0.25, which="both")
    # endfor


def values_for_fixed_log_axis(values: Any) -> tuple[np.ndarray, int, int]:
    """Clip positive values to the fixed 1e-2--1e1 logarithmic axis."""
    raw = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype=float)
    plotted = np.full(raw.shape, np.nan, dtype=float)
    positive = np.isfinite(raw) & (raw > 0.0)
    plotted[positive] = np.clip(
        raw[positive],
        SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS[0],
        SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS[1],
    )
    number_below = int(np.count_nonzero(
        positive & (raw < SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS[0])
    ))
    number_above = int(np.count_nonzero(
        positive & (raw > SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS[1])
    ))
    return plotted, number_below, number_above


def annotate_log_clipping(ax: plt.Axes, below: int, above: int) -> None:
    """Annotate values clipped by the standardized logarithmic range."""
    notes: list[str] = []
    if below:
        notes.append(
            f"{below} below {SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS[0]:g}"
        )
    # endif
    if above:
        notes.append(
            f"{above} above {SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS[1]:g}"
        )
    # endif
    if notes:
        ax.text(
            0.99, 0.03, "; ".join(notes), transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9,
            bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "0.7"},
        )
    # endif


def status_masks(status_values: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return pass, fail, and undefined masks from Barlow status strings."""
    status = pd.Series(status_values, dtype="object").fillna("undefined")
    passed = status.eq("pass").to_numpy(dtype=bool)
    failed = status.eq("fail").to_numpy(dtype=bool)
    undefined = ~(passed | failed)
    return passed, failed, undefined


def plot_systematic_status_markers(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    status_values: Any,
    *,
    label_prefix: str = "",
) -> None:
    """Plot filled pass markers, open fail markers, and x for undefined."""
    passed, failed, undefined = status_masks(status_values)
    prefix = f"{label_prefix}: " if label_prefix else ""
    if np.any(passed):
        ax.plot(
            x[passed], y[passed], "o", color="C0",
            label=f"{prefix}passes Barlow",
        )
    # endif
    if np.any(failed):
        line = ax.plot(
            x[failed], y[failed], "o", color="C0", markerfacecolor="none",
            label=f"{prefix}fails Barlow",
        )[0]
        line.set_markeredgewidth(1.5)
    # endif
    if np.any(undefined):
        ax.plot(
            x[undefined], y[undefined], "x", color="C0",
            label=f"{prefix}Barlow undefined",
        )
    # endif


def draw_systematic_comparison(
    axes: np.ndarray,
    *,
    x: np.ndarray,
    parameter: str,
    top_series: list[dict[str, Any]],
    systematic: Any,
    status: Any,
    systematic_definition: str,
    show_legends: bool = True,
    show_xlabel: bool = True,
) -> None:
    """Draw the common value, assigned-systematic, and normalized-size panels."""
    for series in top_series:
        axes[0].errorbar(
            x,
            np.asarray(series["values"], dtype=float),
            yerr=np.asarray(series["errors"], dtype=float),
            fmt=series.get("fmt", "o"),
            capsize=2,
            label=series["label"],
        )
    # endfor
    axes[0].set_ylabel(PARAMETER_LABELS[parameter])
    if show_legends:
        axes[0].legend()
    # endif

    systematic_array = np.asarray(systematic, dtype=float)
    plot_systematic_status_markers(
        axes[1], x, systematic_array, status,
    )
    axes[1].set_ylabel("Assigned systematic")
    axes[1].text(
        0.02, 0.95, systematic_definition,
        transform=axes[1].transAxes, ha="left", va="top", fontsize=9,
        bbox={"facecolor": "white", "alpha": 0.80, "edgecolor": "0.7"},
    )
    if show_legends:
        axes[1].legend()
    # endif

    nominal_errors = np.asarray(top_series[0]["errors"], dtype=float)
    ratio = np.divide(
        systematic_array,
        nominal_errors,
        out=np.full(systematic_array.shape, np.nan, dtype=float),
        where=np.isfinite(nominal_errors) & (nominal_errors > 0.0),
    )
    ratio_plot, below, above = values_for_fixed_log_axis(ratio)
    plot_systematic_status_markers(axes[2], x, ratio_plot, status)
    axes[2].axhline(
        1.0, linewidth=1.0, linestyle="--",
        label=r"$\delta_{\rm syst}/\sigma_{\rm stat}=1$",
    )
    annotate_log_clipping(axes[2], below, above)
    axes[2].set_ylabel(r"$\delta_{\rm syst}/\sigma_{\rm stat}$")
    if show_xlabel:
        axes[2].set_xlabel("Combined kinematic-bin number")
    # endif
    if show_legends:
        axes[2].legend()
    # endif
    configure_systematic_comparison_axes(axes, parameter)


def write_systematic_summary_canvas(
    *,
    output_dir: Path,
    filename_stem: str,
    draw_parameter: Any,
) -> dict[str, str]:
    """Write a 3xN overview canvas excluding unpolarized modulations."""
    ensure_directory(output_dir)
    ncols = len(PUBLISHED_SYSTEMATIC_PARAMETERS)
    fig, axes = plt.subplots(
        3, ncols, figsize=(4.8 * ncols, 10.5), sharex="col", squeeze=False,
    )
    for column, parameter in enumerate(PUBLISHED_SYSTEMATIC_PARAMETERS):
        draw_parameter(
            parameter,
            axes[:, column],
            show_legends=(column == 0),
            show_xlabel=True,
        )
        axes[0, column].set_title(PARAMETER_LABELS[parameter])
        if column != 0:
            for row in range(3):
                axes[row, column].set_ylabel("")
            # endfor
        # endif
    # endfor
    fig.tight_layout()
    png_path = output_dir / f"{filename_stem}.png"
    pdf_path = output_dir / f"{filename_stem}.pdf"
    fig.savefig(png_path, dpi=SYSTEMATIC_COMPARISON_DPI)
    fig.savefig(pdf_path)
    plt.close(fig)
    return {"png": str(png_path), "pdf": str(pdf_path)}


def period_fit_quality_mask(
    frame: pd.DataFrame,
    period: str,
) -> np.ndarray:
    """Return True only for trustworthy period-only diagnostic fits."""
    return (
        frame[f"period_fit_valid_{period}"].astype(bool).to_numpy()
        & frame[
            f"period_fit_accurate_covariance_{period}"
        ].astype(bool).to_numpy()
        & frame[
            f"period_fit_positive_definite_covariance_{period}"
        ].astype(bool).to_numpy()
        & ~frame[
            f"period_fit_parameters_at_limit_{period}"
        ].astype(bool).to_numpy()
    )


def plot_parameter_summaries(
    frame: pd.DataFrame,
    output_dir: Path,
    include_target_axis_uncertainty: bool,
) -> list[str]:
    """Write the original one-parameter, all-24-bin summary plots."""
    ensure_directory(output_dir)
    paths: list[str] = []
    bins = frame["bin_number"].to_numpy()

    for parameter in PHYSICS_PARAMETERS:
        fig, ax = plt.subplots(figsize=(14, 5.5))
        point_to_point_column = f"{parameter}_point_to_point_systematic"
        if point_to_point_column in frame.columns:
            point_to_point = pd.to_numeric(
                frame[point_to_point_column], errors="coerce"
            ).to_numpy(dtype=float)
            ax.bar(
                bins,
                point_to_point,
                width=0.70,
                bottom=0.0,
                alpha=0.28,
                label="Point-to-point systematic",
                zorder=0,
            )
        # endif
        ax.errorbar(
            bins,
            frame[parameter],
            yerr=frame[f"{parameter}_stat"],
            marker="o",
            linestyle="none",
            capsize=2,
            label="Statistical uncertainty",
        )
        if include_target_axis_uncertainty:
            target_axis_uncertainty = pd.to_numeric(
                frame[f"{parameter}_target_axis_study_envelope"],
                errors="coerce",
            ).to_numpy(dtype=float)
            if np.isfinite(target_axis_uncertainty).any():
                ax.errorbar(
                    bins,
                    frame[parameter],
                    yerr=target_axis_uncertainty,
                    marker="none",
                    linestyle="none",
                    capsize=4,
                    linewidth=1.0,
                    label="Target-axis systematic",
                )
            # endif
        # endif
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("Combined kinematic-bin number")
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_xticks(bins)
        apply_parameter_y_limits(ax, parameter)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        path = output_dir / f"{parameter}_summary.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def make_grouped_physics_canvas(
    *,
    figsize: tuple[float, float] = (18.0, 12.0),
) -> tuple[plt.Figure, dict[str, plt.Axes], list[plt.Axes]]:
    """
    Create the common 3x4 physics canvas used by the aggregate plots.

    Row 1 contains the two UU ratios, row 2 contains the three fitted
    single-spin ratios, and row 3 contains the two fitted double-spin ratios.
    The remaining five cells are left blank.
    """
    fig, axes = plt.subplots(
        3,
        4,
        figsize=figsize,
        sharex=True,
        squeeze=False,
    )

    panel_positions: dict[str, tuple[int, int]] = {
        "u1": (0, 0),
        "u2": (0, 1),
        "lu1": (1, 0),
        "ul1": (1, 1),
        "ul2": (1, 2),
        "ll0": (2, 0),
        "ll1": (2, 1),
    }

    axes_by_parameter: dict[str, plt.Axes] = {}
    axes_in_order: list[plt.Axes] = []
    for parameter in AGGREGATED_PANEL_ORDER:
        row_index, column_index = panel_positions[parameter]
        ax = axes[row_index, column_index]
        axes_by_parameter[parameter] = ax
        axes_in_order.append(ax)
    # endfor

    used_positions = set(panel_positions.values())
    for row_index in range(3):
        for column_index in range(4):
            if (row_index, column_index) not in used_positions:
                axes[row_index, column_index].axis("off")
            # endif
        # endfor
    # endfor

    return fig, axes_by_parameter, axes_in_order

def plot_aggregated_by_x(
    frame: pd.DataFrame,
    output_dir: Path,
    include_target_axis_uncertainty: bool,
) -> list[str]:
    """
    Write one polarization-grouped 3x4 canvas per xB bin.

    Top row: the two UU ratios.
    Middle row: all four single-spin ratios (LU, UL, UL, UT).
    Bottom row: all three double-spin ratios (LL, LL, LT).
    Unused cells are left blank.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        fig, axes_by_parameter, axes_in_order = make_grouped_physics_canvas()

        for parameter in AGGREGATED_PANEL_ORDER:
            ax = axes_by_parameter[parameter]

            x_values = subset["mean_minus_tprime_gev2"]
            point_to_point_column = f"{parameter}_point_to_point_systematic"
            if point_to_point_column in subset.columns:
                point_to_point = pd.to_numeric(
                    subset[point_to_point_column], errors="coerce"
                ).to_numpy(dtype=float)
                x_array = x_values.to_numpy(dtype=float)
                if x_array.size > 1:
                    bar_width = 0.55 * float(np.nanmin(np.diff(x_array)))
                else:
                    bar_width = 0.05
                # endif
                ax.bar(
                    x_array,
                    point_to_point,
                    width=bar_width,
                    bottom=0.0,
                    alpha=0.28,
                    label="Point-to-point systematic",
                    zorder=0,
                )
            # endif
            ax.errorbar(
                x_values,
                subset[parameter],
                yerr=subset[f"{parameter}_stat"],
                marker="o",
                linestyle="none",
                capsize=2,
                label="Statistical uncertainty",
            )
            if include_target_axis_uncertainty:
                target_axis_uncertainty = pd.to_numeric(
                    subset[f"{parameter}_target_axis_study_envelope"],
                    errors="coerce",
                ).to_numpy(dtype=float)
                if np.isfinite(target_axis_uncertainty).any():
                    ax.errorbar(
                        x_values,
                        subset[parameter],
                        yerr=target_axis_uncertainty,
                        marker="none",
                        linestyle="none",
                        capsize=4,
                        linewidth=1.0,
                        label="Target-axis systematic",
                    )
                # endif
            # endif
            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if parameter in PHYSICS_PANEL_ROWS[-1][1]:
                ax.set_xlabel(r"$-t^\\prime$ (GeV$^2$)")
            # endif
        # endfor

        handles, labels = axes_in_order[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower right",
            bbox_to_anchor=(0.96, 0.06),
        )
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_combined.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths

def plot_aggregated_by_period(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """
    Write polarization-grouped canvases comparing the simultaneous result with
    independent Su22, Fa22, and Sp23 diagnostic fits.

    Invalid or covariance-defective period fits are not drawn as ordinary
    error bars.  Their central values are marked with an x so they cannot be
    mistaken for precise measurements with zero or misleading Hesse errors.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        x_values = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
        fig, axes_by_parameter, axes_in_order = make_grouped_physics_canvas()

        for parameter in AGGREGATED_PANEL_ORDER:
            ax = axes_by_parameter[parameter]

            ax.errorbar(
                x_values,
                subset[parameter],
                yerr=subset[f"{parameter}_stat"],
                marker="o",
                linestyle="-",
                capsize=2,
                color=COMBINED_COLOR,
                label="Simultaneous",
            )

            for period in PERIODS:
                quality = period_fit_quality_mask(subset, period)
                values = subset[f"{parameter}_{period}"].to_numpy(dtype=float)
                errors = subset[
                    f"{parameter}_stat_{period}"
                ].to_numpy(dtype=float)

                if np.any(quality):
                    ax.errorbar(
                        x_values[quality],
                        values[quality],
                        yerr=errors[quality],
                        marker="o",
                        linestyle="none",
                        capsize=2,
                        color=PERIOD_COLORS[period],
                        label=PERIOD_LABELS[period],
                    )
                # endif

                invalid = ~quality
                if np.any(invalid):
                    ax.plot(
                        x_values[invalid],
                        values[invalid],
                        marker="x",
                        linestyle="none",
                        markersize=8,
                        markeredgewidth=1.8,
                        color=PERIOD_COLORS[period],
                        label=f"{PERIOD_LABELS[period]} invalid",
                    )
                # endif
            # endfor

            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if parameter in PHYSICS_PANEL_ROWS[-1][1]:
                ax.set_xlabel(r"$-t^\\prime$ (GeV$^2$)")
            # endif
        # endfor

        # Deduplicate legend entries produced in every panel.
        handles: list[Any] = []
        labels: list[str] = []
        for ax in axes_in_order:
            panel_handles, panel_labels = ax.get_legend_handles_labels()
            for handle, label in zip(panel_handles, panel_labels):
                if label not in labels:
                    handles.append(handle)
                    labels.append(label)
                # endif
            # endfor
        # endfor
        fig.legend(
            handles,
            labels,
            loc="lower right",
            bbox_to_anchor=(0.96, 0.06),
        )
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: period diagnostics",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_by_period.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def plot_target_axis_variants(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """Compare nominal, photon-axis-projection, and external-data-informed T fits."""
    ensure_directory(output_dir)
    paths: list[str] = []
    variants = tuple(VARIANT_LABELS)

    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        subset = frame.loc[frame["x_index"] == x_index].sort_values("t_index")
        x_values = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
        fig, axes_by_parameter, axes_in_order = make_grouped_physics_canvas()

        for parameter in AGGREGATED_PANEL_ORDER:
            ax = axes_by_parameter[parameter]

            for variant in variants:
                ax.errorbar(
                    x_values,
                    subset[f"{parameter}_{variant}"],
                    yerr=(
                        subset[f"{parameter}_stat"]
                        if variant == "nominal"
                        else None
                    ),
                    marker="o",
                    linestyle="-",
                    linewidth=1.0,
                    capsize=2,
                    color=VARIANT_COLORS[variant],
                    label=VARIANT_LABELS[variant],
                )
            # endfor

            # Draw the target-axis study envelope for diagnostic comparison as narrow, non-overlapping
            # gray rectangles centered on the nominal points.
            x_array = np.asarray(x_values, dtype=np.float64)
            y_array = subset[parameter].to_numpy(dtype=np.float64)
            sys_array = subset[
                f"{parameter}_target_axis_study_envelope"
            ].to_numpy(dtype=np.float64)
            if x_array.size > 1:
                separations = np.diff(np.sort(x_array))
                half_width = 0.10 * float(np.min(separations))
            else:
                half_width = 0.015
            # endif
            for point_index, (x_point, y_point, y_sys) in enumerate(
                zip(x_array, y_array, sys_array)
            ):
                ax.fill_between(
                    [x_point - half_width, x_point + half_width],
                    [y_point - y_sys, y_point - y_sys],
                    [y_point + y_sys, y_point + y_sys],
                    color="0.65",
                    alpha=0.45,
                    linewidth=0.0,
                    label=(
                        "Target-axis study envelope (diagnostic only)"
                        if point_index == 0 else None
                    ),
                    zorder=1,
                )
            # endfor

            ax.axhline(0.0, linewidth=0.8)
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            apply_parameter_y_limits(ax, parameter)
            ax.grid(alpha=0.25)
            if parameter in PHYSICS_PANEL_ROWS[-1][1]:
                ax.set_xlabel(r"$-t^\\prime$ (GeV$^2$)")
            # endif
        # endfor

        handles, labels = axes_in_order[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower right",
            bbox_to_anchor=(0.97, 0.055),
        )
        fig.suptitle(
            rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: target-axis variants",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
        path = output_dir / f"xB_bin_{x_index + 1}_target_axis_variants.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(str(path))
    # endfor
    return paths


def _period_stability_references(frame: pd.DataFrame, parameter: str) -> dict[str, np.ndarray]:
    """Return simultaneous zero-UU and independent-period weighted references."""
    n = len(frame)
    combined = frame[parameter].to_numpy(dtype=float)
    combined_error = frame[f"{parameter}_stat"].to_numpy(dtype=float)
    weighted_mean = np.full(n, np.nan, dtype=float)
    weighted_error = np.full(n, np.nan, dtype=float)
    for i in range(n):
        vals, weights = [], []
        for period in PERIODS:
            value = float(frame.iloc[i][f"{parameter}_{period}"])
            sigma = float(frame.iloc[i][f"{parameter}_stat_{period}"])
            if (period_fit_quality_mask(frame, period)[i]
                    and np.isfinite(value) and np.isfinite(sigma) and sigma > 0.0):
                vals.append(value); weights.append(1.0 / sigma**2)
            # endif
        # endfor
        if len(vals) >= 2:
            w = np.asarray(weights); v = np.asarray(vals)
            weighted_mean[i] = float(np.sum(w * v) / np.sum(w))
            weighted_error[i] = math.sqrt(1.0 / float(np.sum(w)))
        # endif
    # endfor
    return {"combined": combined, "combined_error": combined_error,
            "weighted_mean": weighted_mean, "weighted_error": weighted_error}


def plot_period_stability_published(frame: pd.DataFrame, output_dir: Path) -> list[str]:
    """Plot zero-UU period fits against the simultaneous zero-UU reference."""
    ensure_directory(output_dir)
    paths: list[str] = []
    bin_numbers = frame["bin_number"].to_numpy(dtype=int)
    offsets = {"su22": -0.18, "fa22": 0.0, "sp23": 0.18}
    markers = {"su22": "o", "fa22": "s", "sp23": "^"}
    period_colors = {"su22": "tab:orange", "fa22": "tab:blue", "sp23": "tab:green"}
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        refs = _period_stability_references(frame, parameter)
        combined, combined_error = refs["combined"], refs["combined_error"]
        fig = plt.figure(figsize=(15, 7.5))
        grid = fig.add_gridspec(2, 1, height_ratios=(2.2, 1.0), hspace=0.06)
        ax = fig.add_subplot(grid[0]); pull_ax = fig.add_subplot(grid[1], sharex=ax)
        cv = np.isfinite(combined) & np.isfinite(combined_error) & (combined_error > 0)
        ax.errorbar(bin_numbers[cv], combined[cv], yerr=combined_error[cv], color="black",
                    marker="o", ms=4.5, lw=1.0, capsize=2, label="Simultaneous fit", zorder=4)
        for period in PERIODS:
            values = frame[f"{parameter}_{period}"].to_numpy(dtype=float)
            errors = frame[f"{parameter}_stat_{period}"].to_numpy(dtype=float)
            valid = period_fit_quality_mask(frame, period) & np.isfinite(values) & np.isfinite(errors) & (errors > 0)
            x = bin_numbers.astype(float) + offsets[period]
            ax.errorbar(x[valid], values[valid], yerr=errors[valid], marker=markers[period], linestyle="none",
                        capsize=2, color=period_colors[period], label=PERIOD_LABELS[period])
            denom2 = errors**2 - combined_error**2
            pv = valid & cv & (denom2 > 0)
            pulls = np.full(len(frame), np.nan)
            pulls[pv] = (values[pv] - combined[pv]) / np.sqrt(denom2[pv])
            pull_ax.plot(x[pv], pulls[pv], color=period_colors[period], marker=markers[period], linestyle="none")
        # endfor
        ax.axhline(0.0, lw=0.8); ax.set_ylabel(PARAMETER_LABELS[parameter]); apply_parameter_y_limits(ax, parameter)
        ax.grid(alpha=0.25); ax.legend(ncol=4, loc="best"); ax.tick_params(labelbottom=False)
        for y, ls, lw in ((0,"-",1.0),(1,"--",0.6),(-1,"--",0.6),(2,":",0.6),(-2,":",0.6),(2.5,"--",0.8),(-2.5,"--",0.8)):
            pull_ax.axhline(y, color="black" if y == 0 else None, linestyle=ls, linewidth=lw)
        # endfor
        pull_ax.set_ylabel(r"Pull wrt. simultaneous"); pull_ax.set_xlabel("Combined kinematic-bin number")
        pull_ax.set_xticks(bin_numbers); pull_ax.set_xlim(0.4, NUMBER_OF_BINS + 0.6); pull_ax.grid(alpha=0.25)
        fig.suptitle(f"Run-period stability: {PARAMETER_LABELS[parameter]} ($u_1=u_2=0$)", y=0.995)
        fig.tight_layout(rect=(0,0,1,0.975))
        png_path = output_dir / f"period_stability_{parameter}_bins_01_24.png"
        fig.savefig(png_path, dpi=200); plt.close(fig); paths.append(str(png_path))
    # endfor
    return paths


def write_period_stability_excursion_table(frame: pd.DataFrame, output_path: Path, threshold: float = 2.5) -> pd.DataFrame:
    """Write period residuals exceeding threshold relative to simultaneous zero-UU fit."""
    rows: list[dict[str, Any]] = []
    estimator_rows: list[dict[str, Any]] = []
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        refs = _period_stability_references(frame, parameter)
        for i, source_row in frame.iterrows():
            combined = float(refs["combined"][i]); combined_err = float(refs["combined_error"][i])
            wmean = float(refs["weighted_mean"][i]); werr = float(refs["weighted_error"][i])
            estimator_rows.append({"bin_number": int(source_row["bin_number"]), "parameter": parameter,
                "simultaneous": combined, "simultaneous_stat": combined_err, "period_weighted_mean": wmean,
                "period_weighted_mean_stat": werr, "difference": wmean-combined})
            for period in PERIODS:
                value = float(source_row[f"{parameter}_{period}"]); sigma = float(source_row[f"{parameter}_stat_{period}"])
                if not (period_fit_quality_mask(frame, period)[i] and np.isfinite(value) and np.isfinite(sigma)
                        and sigma > 0 and np.isfinite(combined) and np.isfinite(combined_err)):
                    continue
                # endif
                denom2 = sigma**2 - combined_err**2
                if denom2 <= 0: continue
                pull = (value-combined)/math.sqrt(denom2)
                if abs(pull) <= threshold: continue
                rows.append({"bin_number": int(source_row["bin_number"]), "parameter": parameter,
                    "parameter_label": PARAMETER_LABELS[parameter], "period": period, "period_label": PERIOD_LABELS[period],
                    "value": value, "stat_uncertainty": sigma, "simultaneous_value": combined,
                    "simultaneous_stat_uncertainty": combined_err, "pull": pull, "abs_pull": abs(pull),
                    "period_weighted_mean": wmean, "period_weighted_mean_stat": werr,
                    "events_period": int(source_row[f"events_{period}"])})
            # endfor
        # endfor
    # endfor
    excursions = pd.DataFrame(rows)
    if not excursions.empty:
        excursions = excursions.sort_values(["abs_pull","parameter","bin_number","period"], ascending=[False,True,True,True]).reset_index(drop=True)
    # endif
    ensure_directory(output_path.parent); excursions.to_csv(output_path, index=False)
    pd.DataFrame(estimator_rows).to_csv(output_path.parent / "period_stability_reference_estimator_comparison.csv", index=False)
    return excursions




def write_period_stability_observable_likelihood_ratio_tests(
    results: list[dict[str, Any]],
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path,
) -> Path:
    """Exact profiled LRT for period dependence of each published amplitude.

    For a tested amplitude, H0 constrains that amplitude to be common to all
    three periods while the other four published polarized amplitudes are
    independently profiled in each period.  H1 is the fully period-separated
    model already fitted by the period-stability extraction.  Thus H1 has two
    additional parameters per bin and the asymptotic reference distribution is
    chi-square with 2 dof per bin (48 dof for 24 valid bins).

    The zero-UU period-stability model is used throughout: u1=u2=0.  Dilution
    factors are profiled with exactly the same Gaussian constraints as in the
    nominal period fits.
    """
    ensure_directory(output_dir)
    try:
        from scipy.stats import chi2 as chi2_distribution
    except Exception:
        chi2_distribution = None
    # endtry

    results_by_bin = {int(item["bin_number"]): item for item in results}
    rows: list[dict[str, Any]] = []

    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        total_stat = 0.0
        total_df = 0
        for bin_number in range(1, NUMBER_OF_BINS + 1):
            result = results_by_bin[bin_number]
            period_fits = result["period_fits"]

            period_nlls = {}
            for period in PERIODS:
                period_nlls[period], _ = make_bin_nll(
                    events, run_states, dilution_records, bin_number,
                    "nominal", active_periods=(period,),
                )
            # endfor

            # One common tested amplitude; all other published amplitudes are
            # period-specific.  The three active dilution factors are also
            # profiled.  u1 and u2 are fixed to zero by construction below.
            names = [parameter]
            for period in PERIODS:
                for other in PUBLISHED_SYSTEMATIC_PARAMETERS:
                    if other != parameter:
                        names.append(f"{other}_{period}")
                    # endif
                # endfor
            # endfor
            names.extend(f"f_{period}" for period in PERIODS)

            starts = {}
            common_weights = []
            common_values = []
            for period in PERIODS:
                fit = period_fits[period]
                value = float(fit["values"][parameter])
                error = float(fit["errors"][parameter])
                if np.isfinite(error) and error > 0.0:
                    common_values.append(value)
                    common_weights.append(1.0 / error**2)
                # endif
                for other in PUBLISHED_SYSTEMATIC_PARAMETERS:
                    if other != parameter:
                        starts[f"{other}_{period}"] = float(fit["values"][other])
                    # endif
                # endfor
                starts[f"f_{period}"] = float(
                    dilution_records[(period, bin_number)].value
                )
            # endfor
            if common_weights:
                starts[parameter] = float(
                    np.average(common_values, weights=common_weights)
                )
            else:
                starts[parameter] = float(
                    result["variants"]["nominal"]["values"][parameter]
                )
            # endif

            central_dilutions = {
                period: float(dilution_records[(period, bin_number)].value)
                for period in PERIODS
            }

            def joint_nll(*args: float) -> float:
                pars = dict(zip(names, args))
                total = 0.0
                for period in PERIODS:
                    physics = {"u1": 0.0, "u2": 0.0}
                    for published in PUBLISHED_SYSTEMATIC_PARAMETERS:
                        physics[published] = (
                            float(pars[parameter])
                            if published == parameter
                            else float(pars[f"{published}_{period}"])
                        )
                    # endfor
                    dilution_values = {
                        f"f_{p}": (
                            float(pars[f"f_{p}"])
                            if p == period else central_dilutions[p]
                        )
                        for p in PERIODS
                    }
                    total += float(period_nlls[period](
                        **physics, **dilution_values
                    ))
                # endfor
                return total

            start_vector = [float(starts[name]) for name in names]

            def configured_candidate(scale: float = 1.0) -> Minuit:
                candidate_start = list(start_vector)
                if scale != 1.0:
                    for i, name in enumerate(names):
                        if name.startswith("f_"):
                            continue
                        # endif
                        candidate_start[i] *= scale
                    # endfor
                # endif
                candidate = Minuit(joint_nll, *candidate_start, name=names)
                candidate.errordef = Minuit.LIKELIHOOD
                candidate.print_level = 0
                candidate.strategy = 1
                for name in names:
                    base = name.split("_", 1)[0]
                    if name in PUBLISHED_SYSTEMATIC_PARAMETERS:
                        candidate.limits[name] = PARAMETER_LIMITS[name]
                    elif base in PUBLISHED_SYSTEMATIC_PARAMETERS:
                        candidate.limits[name] = PARAMETER_LIMITS[base]
                    elif name.startswith("f_"):
                        period = name[2:]
                        record = dilution_records[(period, bin_number)]
                        central = float(record.value)
                        sigma = float(record.stat_uncertainty)
                        width = max(8.0 * sigma, 0.20 * central, 0.02)
                        candidate.limits[name] = (
                            max(1.0e-6, central - width), central + width
                        )
                        if sigma == 0.0:
                            candidate.fixed[name] = True
                        # endif
                    # endif
                # endfor
                candidate.migrad(ncall=120000)
                if not candidate.fmin.is_valid:
                    candidate.simplex(ncall=50000)
                    candidate.strategy = 2
                    candidate.migrad(ncall=180000)
                # endif
                candidate.hesse()
                return candidate

            candidates = [configured_candidate(1.0), configured_candidate(0.5)]
            valid_candidates = [
                candidate for candidate in candidates
                if candidate.fmin.is_valid and math.isfinite(float(candidate.fval))
            ]
            null_fit = min(
                valid_candidates if valid_candidates else candidates,
                key=lambda candidate: float(candidate.fval)
                if math.isfinite(float(candidate.fval)) else math.inf,
            )

            null_nll = float(null_fit.fval)
            alt_nll = float(sum(
                float(period_fits[period]["minimum_nll"])
                for period in PERIODS
            ))
            valid = (
                null_fit.fmin.is_valid
                and math.isfinite(null_nll)
                and math.isfinite(alt_nll)
            )
            stat = max(0.0, 2.0 * (null_nll - alt_nll)) if valid else math.nan
            df = 2 if valid else 0
            pvalue = (
                float(chi2_distribution.sf(stat, df))
                if valid and chi2_distribution is not None else math.nan
            )
            rows.append({
                "scope": "bin", "parameter": parameter,
                "bin_number": bin_number, "minus2logLambda": stat,
                "df": df, "pvalue": pvalue,
                "null_minimum_nll": null_nll,
                "alternative_minimum_nll": alt_nll,
                "null_fit_valid": bool(null_fit.fmin.is_valid),
                "null_fit_accurate_covariance": bool(null_fit.fmin.has_accurate_covar),
                "null_fit_positive_definite_covariance": bool(null_fit.fmin.has_posdef_covar),
                "null_fit_parameters_at_limit": bool(null_fit.fmin.has_parameters_at_limit),
                "null_fit_edm": float(null_fit.fmin.edm),
                "null_fit_nfcn": int(null_fit.fmin.nfcn),
            })
            if valid:
                total_stat += stat
                total_df += df
            # endif
            print(
                f"[observable LRT] {parameter} bin {bin_number:02d}: "
                f"-2dlnL={stat:.3f}, df={df}, p={pvalue:.4g}, "
                f"valid={bool(null_fit.fmin.is_valid)}",
                flush=True,
            )
        # endfor

        total_pvalue = (
            float(chi2_distribution.sf(total_stat, total_df))
            if total_df > 0 and chi2_distribution is not None else math.nan
        )
        rows.append({
            "scope": "all_24_bins", "parameter": parameter,
            "bin_number": np.nan, "minus2logLambda": total_stat,
            "df": total_df, "pvalue": total_pvalue,
            "null_minimum_nll": np.nan, "alternative_minimum_nll": np.nan,
            "null_fit_valid": True,
            "null_fit_accurate_covariance": np.nan,
            "null_fit_positive_definite_covariance": np.nan,
            "null_fit_parameters_at_limit": np.nan,
            "null_fit_edm": np.nan, "null_fit_nfcn": np.nan,
        })
        print(
            f"[observable LRT] {parameter} TOTAL: "
            f"-2dlnL={total_stat:.3f}, df={total_df}, p={total_pvalue:.6g}",
            flush=True,
        )
    # endfor

    path = output_dir / "period_stability_observable_likelihood_ratio_tests.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def write_period_stability_global_diagnostics(
    frame: pd.DataFrame,
    output_dir: Path,
    cuts: Mapping[tuple[str, int], CutRecord],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, str]:
    """Quantify global and structured run-period consistency for zero-UU fits."""
    ensure_directory(output_dir)
    try:
        from scipy.stats import chi2 as chi2_distribution, pearsonr, spearmanr
    except Exception:
        chi2_distribution = None
        pearsonr = spearmanr = None
    # endtry

    # Exact all-five-amplitude likelihood-ratio test.  For each bin the null
    # has five common polarized amplitudes; the alternative has five amplitudes
    # independently in each of three periods, hence 10 additional parameters.
    lrt_rows = []
    total_stat = 0.0
    total_df = 0
    for _, row in frame.sort_values("bin_number").iterrows():
        deltas = np.asarray(
            [float(row.get(f"period_delta_nll_{period}", np.nan)) for period in PERIODS],
            dtype=float,
        )
        valid = np.all(np.isfinite(deltas))
        stat = 2.0 * float(np.sum(deltas)) if valid else math.nan
        df = 2 * len(PUBLISHED_SYSTEMATIC_PARAMETERS) if valid else 0
        pvalue = (
            float(chi2_distribution.sf(stat, df))
            if valid and chi2_distribution is not None else math.nan
        )
        lrt_rows.append({
            "scope": "bin", "bin_number": int(row["bin_number"]),
            "minus2logLambda": stat, "df": df, "pvalue": pvalue,
        })
        if valid:
            total_stat += stat
            total_df += df
        # endif
    # endfor
    global_p = (
        float(chi2_distribution.sf(total_stat, total_df))
        if total_df > 0 and chi2_distribution is not None else math.nan
    )
    lrt_rows.append({
        "scope": "all_24_bins", "bin_number": np.nan,
        "minus2logLambda": total_stat, "df": total_df, "pvalue": global_p,
    })
    lrt_path = output_dir / "period_stability_global_likelihood_ratio.csv"
    pd.DataFrame(lrt_rows).to_csv(lrt_path, index=False)

    # Observable-by-observable heterogeneity test.  This is a Wald/Cochran-Q
    # test based on the three independent period estimators, not an LRT; it is
    # intentionally labelled as such so the statistical interpretation is clear.
    observable_rows = []
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        q_total = 0.0
        df_total = 0
        for _, row in frame.sort_values("bin_number").iterrows():
            vals, weights = [], []
            for period in PERIODS:
                value = float(row[f"{parameter}_{period}"])
                sigma = float(row[f"{parameter}_stat_{period}"])
                if np.isfinite(value) and np.isfinite(sigma) and sigma > 0.0:
                    vals.append(value)
                    weights.append(1.0 / sigma**2)
                # endif
            # endfor
            if len(vals) < 2:
                continue
            # endif
            v = np.asarray(vals); w = np.asarray(weights)
            mean = float(np.sum(w * v) / np.sum(w))
            q = float(np.sum(w * (v - mean)**2))
            df_q = len(vals) - 1
            observable_rows.append({
                "scope": "bin", "parameter": parameter,
                "bin_number": int(row["bin_number"]), "Q": q, "df": df_q,
                "pvalue": float(chi2_distribution.sf(q, df_q)) if chi2_distribution is not None else math.nan,
            })
            q_total += q
            df_total += df_q
        # endfor
        observable_rows.append({
            "scope": "all_24_bins", "parameter": parameter,
            "bin_number": np.nan, "Q": q_total, "df": df_total,
            "pvalue": float(chi2_distribution.sf(q_total, df_total)) if df_total and chi2_distribution is not None else math.nan,
        })
    # endfor
    observable_path = output_dir / "period_stability_observable_wald_tests.csv"
    pd.DataFrame(observable_rows).to_csv(observable_path, index=False)

    # Build one row per bin/observable/period with the simultaneous-reference
    # pull and experimental context needed to trace coherent patterns.
    context_rows = []
    pt_by_period = {}
    for period in PERIODS:
        state = run_states[period]
        qtot = np.asarray(state["q_plus"]) + np.asarray(state["q_minus"])
        qsum = float(np.sum(qtot))
        pt_by_period[period] = (
            float(np.sum(np.abs(np.asarray(state["pt"])) * qtot) / qsum)
            if qsum > 0.0 else math.nan
        )
    # endfor
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        refs = _period_stability_references(frame, parameter)
        for i, row in frame.iterrows():
            combined = float(refs["combined"][i])
            combined_err = float(refs["combined_error"][i])
            bin_number = int(row["bin_number"])
            for period in PERIODS:
                value = float(row[f"{parameter}_{period}"])
                sigma_stat = float(row[f"{parameter}_stat_{period}"])
                denom2 = sigma_stat**2 - combined_err**2
                pull = ((value - combined) / math.sqrt(denom2)) if denom2 > 0.0 else math.nan
                cut = cuts[(period, bin_number)]
                mu = 0.5 * (cut.low_gev2 + cut.high_gev2)
                mx2_sigma = 0.25 * (cut.high_gev2 - cut.low_gev2)  # nominal is mu +/- 2 sigma
                context_rows.append({
                    "bin_number": bin_number, "parameter": parameter, "period": period,
                    "period_value": value, "period_stat": sigma_stat,
                    "simultaneous_value": combined, "simultaneous_stat": combined_err,
                    "signed_difference": value - combined, "pull": pull,
                    "events": int(row[f"events_{period}"]),
                    "dilution": float(row[f"dilution_{period}"]),
                    "dilution_stat": float(row[f"dilution_stat_{period}"]),
                    "charge_weighted_mean_abs_target_polarization": pt_by_period[period],
                    "mx2_mu_gev2": mu, "mx2_sigma_gev2": mx2_sigma,
                    "mx2_low_gev2": cut.low_gev2, "mx2_high_gev2": cut.high_gev2,
                })
            # endfor
        # endfor
    # endfor
    context = pd.DataFrame(context_rows)
    context_path = output_dir / "period_stability_full_context.csv"
    context.to_csv(context_path, index=False)

    pull_rows = []
    for (parameter, period), group in context.groupby(["parameter", "period"]):
        x = group["pull"].to_numpy(dtype=float)
        x = x[np.isfinite(x)]
        pull_rows.append({
            "parameter": parameter, "period": period, "n": len(x),
            "mean_pull": float(np.mean(x)) if len(x) else math.nan,
            "rms_pull": float(np.sqrt(np.mean(x**2))) if len(x) else math.nan,
            "sample_std_pull": float(np.std(x, ddof=1)) if len(x) > 1 else math.nan,
            "n_positive": int(np.count_nonzero(x > 0.0)),
            "n_negative": int(np.count_nonzero(x < 0.0)),
            "n_abs_gt_2p5": int(np.count_nonzero(np.abs(x) > 2.5)),
        })
    # endfor
    pull_path = output_dir / "period_stability_signed_pull_summary.csv"
    pd.DataFrame(pull_rows).to_csv(pull_path, index=False)

    corr_rows = []
    covariates = ("mx2_sigma_gev2", "dilution", "dilution_stat", "events")
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        for period_key in (*PERIODS, "all_periods"):
            subset = context.loc[context["parameter"] == parameter]
            if period_key != "all_periods":
                subset = subset.loc[subset["period"] == period_key]
            # endif
            for covariate in covariates:
                x = subset[covariate].to_numpy(dtype=float)
                y = subset["pull"].to_numpy(dtype=float)
                good = np.isfinite(x) & np.isfinite(y)
                if np.count_nonzero(good) >= 3 and np.std(x[good]) > 0.0 and np.std(y[good]) > 0.0:
                    pr = pearsonr(x[good], y[good]) if pearsonr is not None else (math.nan, math.nan)
                    sr = spearmanr(x[good], y[good]) if spearmanr is not None else (math.nan, math.nan)
                    pearson_r, pearson_p = float(pr[0]), float(pr[1])
                    spearman_r, spearman_p = float(sr[0]), float(sr[1])
                else:
                    pearson_r = pearson_p = spearman_r = spearman_p = math.nan
                # endif
                corr_rows.append({
                    "parameter": parameter, "period": period_key, "covariate": covariate,
                    "n": int(np.count_nonzero(good)), "pearson_r": pearson_r,
                    "pearson_pvalue": pearson_p, "spearman_rho": spearman_r,
                    "spearman_pvalue": spearman_p,
                })
            # endfor
        # endfor
    # endfor
    corr_path = output_dir / "period_stability_pull_context_correlations.csv"
    pd.DataFrame(corr_rows).to_csv(corr_path, index=False)

    # Coherent multiplicative period-scale test.  Fit period = s * simultaneous
    # through the origin using the correlated residual variance sigma_p^2-sigma_c^2.
    scale_groups = {
        "UL": ("ul1", "ul2"),
        "LL": ("ll0", "ll1"),
        "all_target_polarized": ("ul1", "ul2", "ll0", "ll1"),
    }
    scale_rows = []
    for period in PERIODS:
        for group_name, parameters in scale_groups.items():
            xs, ys, ws = [], [], []
            for parameter in parameters:
                refs = _period_stability_references(frame, parameter)
                for i, row in frame.iterrows():
                    c = float(refs["combined"][i]); ce = float(refs["combined_error"][i])
                    y = float(row[f"{parameter}_{period}"]); ye = float(row[f"{parameter}_stat_{period}"])
                    variance = ye**2 - ce**2
                    if np.isfinite(c) and np.isfinite(y) and variance > 0.0:
                        xs.append(c); ys.append(y); ws.append(1.0 / variance)
                    # endif
                # endfor
            # endfor
            x = np.asarray(xs); y = np.asarray(ys); w = np.asarray(ws)
            denom = float(np.sum(w * x**2)) if len(x) else 0.0
            scale = float(np.sum(w * x * y) / denom) if denom > 0.0 else math.nan
            scale_err = math.sqrt(1.0 / denom) if denom > 0.0 else math.nan
            chi2 = float(np.sum(w * (y - scale * x)**2)) if denom > 0.0 else math.nan
            df_scale = max(0, len(x) - 1)
            scale_rows.append({
                "period": period, "group": group_name, "n": len(x), "scale": scale,
                "scale_stat": scale_err, "scale_minus_one_sigma": (scale - 1.0) / scale_err if scale_err > 0.0 else math.nan,
                "chi2": chi2, "df": df_scale,
                "chi2_pvalue": float(chi2_distribution.sf(chi2, df_scale)) if df_scale and chi2_distribution is not None else math.nan,
            })
        # endfor
    # endfor
    scale_path = output_dir / "period_stability_period_scale_tests.csv"
    pd.DataFrame(scale_rows).to_csv(scale_path, index=False)

    return {
        "global_lrt": str(lrt_path), "observable_wald": str(observable_path),
        "context": str(context_path), "pull_summary": str(pull_path),
        "correlations": str(corr_path), "scale_tests": str(scale_path),
    }



def _copy_run_states_with_charge_distortion(
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    distorted_period: str,
    epsilon: float,
    mode: str,
) -> dict[str, dict[str, np.ndarray]]:
    """Return run states with a controlled charge distortion.

    ``target_epoch`` applies Q_h -> Q_h (1 + epsilon*sign(Pt)) to both
    helicities in a target-polarization epoch.  This is the physically motivated
    Faraday-cup scale test: a charge calibration offset tied to target epoch/sign.

    ``double_spin`` applies Q_h -> Q_h (1 + epsilon*h*sign(Pt)).  This is the
    normalization mode maximally degenerate with a constant double-spin term and
    is retained as a useful response/control test.
    """
    if mode not in {"target_epoch", "double_spin"}:
        raise ValueError(f"Unknown charge-distortion mode: {mode}")
    # endif

    result: dict[str, dict[str, np.ndarray]] = {}
    for period in PERIODS:
        result[period] = {
            name: np.array(values, copy=True)
            for name, values in run_states[period].items()
        }
        if period != distorted_period:
            continue
        # endif

        sign_pt = np.sign(result[period]["pt"])
        if mode == "target_epoch":
            factor = 1.0 + epsilon * sign_pt
            result[period]["q_plus"] *= factor
            result[period]["q_minus"] *= factor
        else:
            result[period]["q_plus"] *= 1.0 + epsilon * sign_pt
            result[period]["q_minus"] *= 1.0 - epsilon * sign_pt
        # endif

        if np.any(result[period]["q_plus"] <= 0.0) or np.any(result[period]["q_minus"] <= 0.0):
            raise RuntimeError("Charge-distortion diagnostic produced nonpositive charge.")
        # endif
    # endfor
    return result


def write_charge_normalization_diagnostics(
    *,
    frame: pd.DataFrame,
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path,
    run_full_scan: bool = False,
) -> dict[str, str]:
    """Audit charge bookkeeping and optionally stress-test two normalization modes.

    The audit records the actual Q+ and Q- values supplied to the likelihood by
    run, period, and contiguous target-sign epoch.  The stress test then scans
    both a target-epoch Faraday-cup scale distortion, (1 + epsilon*s), and a
    double-spin distortion, (1 + epsilon*h*s), at +/-0.25, 0.5, 1, 2, and 5%.
    All 24 bins are refit with u1=u2=0 so the response of every published
    polarized amplitude can be compared directly.
    """
    ensure_directory(output_dir)

    # ---- run-level and period-level charge audit -------------------------
    run_rows: list[dict[str, Any]] = []
    period_rows: list[dict[str, Any]] = []
    epoch_rows: list[dict[str, Any]] = []
    for period in PERIODS:
        state = run_states[period]
        runs = np.asarray(state["run"], dtype=int)
        pt = np.asarray(state["pt"], dtype=float)
        qp = np.asarray(state["q_plus"], dtype=float)
        qm = np.asarray(state["q_minus"], dtype=float)
        qsum = qp + qm
        qdiff = qp - qm
        sign_pt = np.sign(pt)
        for r, ptt, qpp, qmm in zip(runs, pt, qp, qm):
            total = qpp + qmm
            run_rows.append({
                "period": period, "run": int(r), "target_polarization": float(ptt),
                "target_sign": int(np.sign(ptt)), "charge_plus": float(qpp),
                "charge_minus": float(qmm), "charge_total": float(total),
                "helicity_charge_asymmetry": float((qpp-qmm)/total) if total > 0 else math.nan,
            })
        # endfor
        total_q = float(np.sum(qsum))
        hel_q = float(np.sum(qdiff))
        ds_q = float(np.sum(sign_pt * qdiff))
        target_q = float(np.sum(sign_pt * qsum))
        period_rows.append({
            "period": period, "number_of_runs": len(runs), "charge_total": total_q,
            "helicity_charge_asymmetry": hel_q/total_q if total_q > 0 else math.nan,
            "target_sign_charge_asymmetry": target_q/total_q if total_q > 0 else math.nan,
            "double_spin_charge_asymmetry": ds_q/total_q if total_q > 0 else math.nan,
            "charge_weighted_mean_target_polarization": float(np.sum(pt*qsum)/total_q) if total_q > 0 else math.nan,
            "charge_weighted_mean_abs_target_polarization": float(np.sum(np.abs(pt)*qsum)/total_q) if total_q > 0 else math.nan,
        })

        # Contiguous target-sign blocks are a definition-free polarization epoch.
        if len(runs):
            start = 0
            for j in range(1, len(runs)+1):
                boundary = j == len(runs) or sign_pt[j] != sign_pt[j-1]
                if not boundary:
                    continue
                # endif
                sl = slice(start, j)
                tq = float(np.sum(qsum[sl]))
                epoch_rows.append({
                    "period": period, "epoch": len([x for x in epoch_rows if x["period"] == period]) + 1,
                    "run_low": int(runs[start]), "run_high": int(runs[j-1]),
                    "target_sign": int(sign_pt[start]), "number_of_runs": int(j-start),
                    "charge_total": tq,
                    "helicity_charge_asymmetry": float(np.sum(qdiff[sl])/tq) if tq > 0 else math.nan,
                    "double_spin_charge_asymmetry": float(np.sum(sign_pt[sl]*qdiff[sl])/tq) if tq > 0 else math.nan,
                    "charge_weighted_mean_target_polarization": float(np.sum(pt[sl]*qsum[sl])/tq) if tq > 0 else math.nan,
                })
                start = j
            # endfor
        # endif
    # endfor
    run_path = output_dir / "charge_audit_by_run.csv"
    period_path = output_dir / "charge_audit_by_period.csv"
    epoch_path = output_dir / "charge_audit_target_sign_epochs.csv"
    pd.DataFrame(run_rows).to_csv(run_path, index=False)
    pd.DataFrame(period_rows).to_csv(period_path, index=False)
    pd.DataFrame(epoch_rows).to_csv(epoch_path, index=False)

    if not run_full_scan:
        return {
            "run_audit": str(run_path), "period_audit": str(period_path),
            "epoch_audit": str(epoch_path), "response": "",
            "summary": "", "scan": "",
        }
    # endif

    # ---- controlled target-epoch and h*s charge distortions -------------
    eps_scan = (0.0025, 0.0050, 0.0100, 0.0200, 0.0500)
    eps0 = 0.0050
    modes = ("target_epoch", "double_spin")
    response_rows: list[dict[str, Any]] = []
    scan_rows: list[dict[str, Any]] = []
    fixed = {"u1": 0.0, "u2": 0.0}
    frame_by_bin = frame.set_index("bin_number")
    for mode in modes:
        for period in PERIODS:
            distorted_states = {
                eps: _copy_run_states_with_charge_distortion(run_states, period, eps, mode)
                for mag in eps_scan for eps in (-mag, +mag)
            }
            for bin_number in range(1, NUMBER_OF_BINS + 1):
                nominal_row = frame_by_bin.loc[bin_number]
                initial = {
                    name: float(nominal_row[f"{name}_{period}"])
                    for name in PHYSICS_PARAMETERS
                    if f"{name}_{period}" in nominal_row.index
                }
                fits: dict[float, dict[str, Any]] = {}
                for eps in sorted(distorted_states):
                    fit = fit_one_variant(
                        events, distorted_states[eps], dilution_records, bin_number, "nominal",
                        active_periods=(period,), initial_values=initial,
                        fixed_physics_parameters=fixed,
                    )
                    fits[eps] = fit
                    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
                        scan_rows.append({
                            "distortion_mode": mode, "period": period, "bin_number": bin_number,
                            "parameter": parameter, "epsilon": eps,
                            "percent_charge_distortion": 100.0*eps,
                            "value": float(fit["values"][parameter]),
                            "stat": float(fit["errors"][parameter]), "fit_valid": bool(fit["valid"]),
                        })
                    # endfor
                # endfor
                for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
                    vm = float(fits[-eps0]["values"][parameter])
                    vp = float(fits[+eps0]["values"][parameter])
                    derivative = (vp-vm)/(2.0*eps0)
                    slopes = []
                    for mag in eps_scan:
                        slopes.append((float(fits[+mag]["values"][parameter]) - float(fits[-mag]["values"][parameter]))/(2.0*mag))
                    # endfor
                    period_value = float(nominal_row[f"{parameter}_{period}"])
                    combined = float(nominal_row[parameter])
                    observed_shift = period_value-combined
                    required = observed_shift/derivative if abs(derivative) > 1.0e-12 else math.nan
                    response_rows.append({
                        "distortion_mode": mode, "period": period, "bin_number": bin_number,
                        "parameter": parameter, "epsilon_minus": -eps0, "value_minus": vm,
                        "epsilon_plus": eps0, "value_plus": vp, "dA_depsilon": derivative,
                        "slope_at_0p25pct": slopes[0], "slope_at_0p5pct": slopes[1],
                        "slope_at_1pct": slopes[2], "slope_at_2pct": slopes[3],
                        "slope_at_5pct": slopes[4],
                        "max_fractional_slope_change": (max(abs(x-derivative) for x in slopes)/abs(derivative)) if abs(derivative)>1e-12 else math.nan,
                        "nominal_period_value": period_value, "simultaneous_value": combined,
                        "observed_period_minus_simultaneous": observed_shift,
                        "epsilon_required_to_match_observed_shift": required,
                        "percent_charge_distortion_required": 100.0*required if np.isfinite(required) else math.nan,
                        "fit_minus_valid": bool(fits[-eps0]["valid"]),
                        "fit_plus_valid": bool(fits[+eps0]["valid"]),
                    })
                # endfor
            # endfor
        # endfor
    # endfor

    response = pd.DataFrame(response_rows)
    response_path = output_dir / "charge_distortion_response.csv"
    response.to_csv(response_path, index=False)
    scan_path = output_dir / "charge_distortion_scan.csv"
    pd.DataFrame(scan_rows).to_csv(scan_path, index=False)

    summary_rows: list[dict[str, Any]] = []
    for (mode, period, parameter), group in response.groupby(["distortion_mode", "period", "parameter"]):
        req = group["epsilon_required_to_match_observed_shift"].to_numpy(float)
        der = group["dA_depsilon"].to_numpy(float)
        req = req[np.isfinite(req)]
        der = der[np.isfinite(der)]
        summary_rows.append({
            "distortion_mode": mode, "period": period, "parameter": parameter, "n": len(group),
            "median_dA_depsilon": float(np.median(der)) if len(der) else math.nan,
            "mean_dA_depsilon": float(np.mean(der)) if len(der) else math.nan,
            "median_required_epsilon": float(np.median(req)) if len(req) else math.nan,
            "median_required_percent": float(100.0*np.median(req)) if len(req) else math.nan,
        })
    # endfor
    summary_path = output_dir / "charge_distortion_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    return {
        "run_audit": str(run_path), "period_audit": str(period_path),
        "epoch_audit": str(epoch_path), "response": str(response_path),
        "summary": str(summary_path), "scan": str(scan_path),
    }



def _copy_run_states_with_target_sign_total_scale(
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    period: str,
    target_sign: int,
    scale: float,
) -> dict[str, dict[str, np.ndarray]]:
    """Scale one of the six period x target-sign integrated charge totals.

    Both helicity charges for runs with the selected target sign are multiplied
    by the same factor.  This is the literal test of a Faraday-cup normalization
    error in one period/target-polarization sample.
    """
    if scale <= 0.0:
        raise ValueError("Target-sign charge scale must remain positive.")
    # endif
    result = {
        p: {name: np.array(values, copy=True) for name, values in state.items()}
        for p, state in run_states.items()
    }
    state = result[period]
    mask = np.sign(state["pt"]).astype(int) == int(np.sign(target_sign))
    state["q_plus"][mask] *= scale
    state["q_minus"][mask] *= scale
    return result


def write_six_charge_total_corrections(
    *, frame: pd.DataFrame, events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path, epsilon: float = 0.01,
) -> dict[str, str]:
    """Infer the correction to each of the six period x target-sign charges.

    One charge total is perturbed at a time by +/-epsilon.  For every perturbation
    both the affected period fit and the three-period simultaneous fit are rerun,
    so the motion of the reference mean is included rather than held fixed.
    Only ll0 is floated; all other physics amplitudes are fixed to the relevant
    nominal solution.  The reported correction is the local linear correction
    to the *currently used* charge total that would remove the weighted period-
    minus-simultaneous A_LL displacement if that one charge total were the sole
    source of the discrepancy.
    """
    ensure_directory(output_dir)
    frame_by_bin = frame.set_index("bin_number")
    rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []

    for period in PERIODS:
        state = run_states[period]
        for target_sign in (1, -1):
            sign_mask = np.sign(state["pt"]).astype(int) == target_sign
            nominal_charge = float(np.sum((state["q_plus"] + state["q_minus"])[sign_mask]))
            distorted = {
                e: _copy_run_states_with_target_sign_total_scale(
                    run_states, period, target_sign, 1.0 + e
                )
                for e in (-epsilon, +epsilon)
            }
            group_rows: list[dict[str, Any]] = []
            for bin_number in range(1, NUMBER_OF_BINS + 1):
                nominal = frame_by_bin.loc[bin_number]
                period_initial = {
                    name: float(nominal[f"{name}_{period}"])
                    for name in PHYSICS_PARAMETERS
                }
                combined_initial = {
                    name: float(nominal[name]) for name in PHYSICS_PARAMETERS
                }
                period_fixed = {
                    name: value for name, value in period_initial.items() if name != "ll0"
                }
                combined_fixed = {
                    name: value for name, value in combined_initial.items() if name != "ll0"
                }
                dvals: dict[float, float] = {}
                pvals: dict[float, float] = {}
                cvals: dict[float, float] = {}
                for e in (-epsilon, +epsilon):
                    pf = fit_one_variant(
                        events, distorted[e], dilution_records, bin_number, "nominal",
                        active_periods=(period,), initial_values=period_initial,
                        fixed_physics_parameters=period_fixed,
                    )
                    cf = fit_one_variant(
                        events, distorted[e], dilution_records, bin_number, "nominal",
                        initial_values=combined_initial,
                        fixed_physics_parameters=combined_fixed,
                    )
                    pvals[e] = float(pf["values"]["ll0"])
                    cvals[e] = float(cf["values"]["ll0"])
                    dvals[e] = pvals[e] - cvals[e]
                # endfor
                slope = (dvals[+epsilon] - dvals[-epsilon]) / (2.0 * epsilon)
                d0 = float(nominal[f"ll0_{period}"] - nominal["ll0"])
                required = -d0 / slope if abs(slope) > 1.0e-12 else math.nan
                pstat = float(nominal[f"ll0_stat_{period}"])
                cstat = float(nominal["ll0_stat"])
                nested_var = max(pstat * pstat - cstat * cstat, 0.0)
                row = {
                    "period": period, "target_sign": target_sign,
                    "bin_number": bin_number, "nominal_charge_total": nominal_charge,
                    "epsilon_probe": epsilon,
                    "period_ll0_minus": pvals[-epsilon], "combined_ll0_minus": cvals[-epsilon],
                    "period_ll0_plus": pvals[+epsilon], "combined_ll0_plus": cvals[+epsilon],
                    "nominal_period_minus_combined": d0,
                    "d_period_minus_combined_depsilon": slope,
                    "required_fractional_charge_correction": required,
                    "required_percent_charge_correction": 100.0 * required if np.isfinite(required) else math.nan,
                    "nested_variance": nested_var,
                }
                rows.append(row)
                group_rows.append(row)
            # endfor

            finite = [r for r in group_rows if np.isfinite(r["d_period_minus_combined_depsilon"])]
            num = den = 0.0
            for r in finite:
                var = r["nested_variance"]
                if var <= 0.0:
                    continue
                # endif
                w = 1.0 / var
                slope = r["d_period_minus_combined_depsilon"]
                d0 = r["nominal_period_minus_combined"]
                num += w * slope * d0
                den += w * slope * slope
            # endfor
            required_global = -num / den if den > 0.0 else math.nan
            summaries.append({
                "period": period, "target_sign": target_sign,
                "nominal_charge_total": nominal_charge,
                "required_fractional_charge_correction": required_global,
                "required_percent_charge_correction": 100.0 * required_global if np.isfinite(required_global) else math.nan,
                "corrected_charge_total": nominal_charge * (1.0 + required_global) if np.isfinite(required_global) else math.nan,
                "weighted_bins": sum(1 for r in finite if r["nested_variance"] > 0.0),
                "interpretation": "correction to currently used charge if this one target-sign total alone explains the A_LL period displacement",
            })
            print(
                f"[six-charge solve] {period} Pt sign {target_sign:+d}: "
                f"required correction={100.0*required_global:+.3f}%",
                flush=True,
            )
        # endfor
    # endfor

    bins_path = output_dir / "six_charge_total_corrections_by_bin.csv"
    summary_path = output_dir / "six_charge_total_corrections_summary.csv"
    pd.DataFrame(rows).to_csv(bins_path, index=False)
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    return {"bins": str(bins_path), "summary": str(summary_path)}



def _copy_run_states_with_period_target_ratio_scale(
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    period_scales: Mapping[str, float],
) -> dict[str, dict[str, np.ndarray]]:
    """Apply antisymmetric target-sign charge scales within each period.

    For period p with parameter delta_p,
        Q(Pt>0) -> Q(Pt>0) * (1 + delta_p)
        Q(Pt<0) -> Q(Pt<0) * (1 - delta_p).
    Thus only the physically relevant Pt+/Pt- charge ratio is changed; the
    nearly invisible common charge normalization of a period is not floated.
    """
    result = {
        p: {name: np.array(values, copy=True) for name, values in state.items()}
        for p, state in run_states.items()
    }
    for period, delta in period_scales.items():
        if abs(delta) >= 0.95:
            raise ValueError("Target-sign ratio distortion must keep both charge scales positive.")
        # endif
        state = result[period]
        signs = np.sign(state["pt"]).astype(int)
        plus = signs > 0
        minus = signs < 0
        state["q_plus"][plus] *= 1.0 + delta
        state["q_minus"][plus] *= 1.0 + delta
        state["q_plus"][minus] *= 1.0 - delta
        state["q_minus"][minus] *= 1.0 - delta
    # endfor
    return result


def write_simultaneous_target_charge_ratio_fit(
    *, frame: pd.DataFrame, events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path, epsilon: float = 0.01,
) -> dict[str, str]:
    """Simultaneously solve the three period Pt+/Pt- charge-ratio corrections.

    The fit has one physically relevant charge parameter per run period rather
    than six absolute charge scales.  A parameter delta_p applies +delta_p to
    the Pt>0 total and -delta_p to the Pt<0 total, so

        R_p' / R_p = (1 + delta_p) / (1 - delta_p),

    where R_p = Q(Pt>0)/Q(Pt<0).  All three deltas float simultaneously and
    every response includes the movement of the three-period simultaneous A_LL
    reference.  A +/-1% finite-difference response gives the simultaneous WLS
    solution.  The resulting point is then evaluated with full likelihood
    refits (not extrapolated) so the CSV reports whether the linear solution is
    actually valid.  If a solution leaves the few-percent regime, it should be
    regarded as evidence that this charge-ratio mechanism cannot plausibly
    account for the observed period displacement, not as a literal calibration.
    """
    ensure_directory(output_dir)
    frame_by_bin = frame.set_index("bin_number")
    labels = [f"{p}_Ptplus_over_Ptminus" for p in PERIODS]

    obs = []
    for bin_number in range(1, NUMBER_OF_BINS + 1):
        nominal = frame_by_bin.loc[bin_number]
        for period in PERIODS:
            d0 = float(nominal[f"ll0_{period}"] - nominal["ll0"])
            pstat = float(nominal[f"ll0_stat_{period}"])
            cstat = float(nominal["ll0_stat"])
            var = max(pstat * pstat - cstat * cstat, 1.0e-12)
            obs.append((bin_number, period, d0, var))
        # endfor
    # endfor

    y = np.asarray([x[2] for x in obs], dtype=float)
    var = np.asarray([x[3] for x in obs], dtype=float)
    sqrt_w = 1.0 / np.sqrt(var)
    response = np.zeros((len(obs), len(PERIODS)), dtype=float)
    probe_rows: list[dict[str, Any]] = []

    for j, changed_period in enumerate(PERIODS):
        distorted = {
            e: _copy_run_states_with_period_target_ratio_scale(
                run_states, {changed_period: e}
            )
            for e in (-epsilon, +epsilon)
        }
        combined_by_eps = {-epsilon: {}, +epsilon: {}}
        affected_by_eps = {-epsilon: {}, +epsilon: {}}

        for bin_number in range(1, NUMBER_OF_BINS + 1):
            nominal = frame_by_bin.loc[bin_number]
            combined_initial = {name: float(nominal[name]) for name in PHYSICS_PARAMETERS}
            combined_fixed = {name: value for name, value in combined_initial.items() if name != "ll0"}
            period_initial = {name: float(nominal[f"{name}_{changed_period}"]) for name in PHYSICS_PARAMETERS}
            period_fixed = {name: value for name, value in period_initial.items() if name != "ll0"}
            for e in (-epsilon, +epsilon):
                cf = fit_one_variant(
                    events, distorted[e], dilution_records, bin_number, "nominal",
                    initial_values=combined_initial, fixed_physics_parameters=combined_fixed,
                )
                pf = fit_one_variant(
                    events, distorted[e], dilution_records, bin_number, "nominal",
                    active_periods=(changed_period,), initial_values=period_initial,
                    fixed_physics_parameters=period_fixed,
                )
                combined_by_eps[e][bin_number] = float(cf["values"]["ll0"])
                affected_by_eps[e][bin_number] = float(pf["values"]["ll0"])
            # endfor
        # endfor

        for i, (bin_number, residual_period, d0, nested_var) in enumerate(obs):
            if residual_period == changed_period:
                dm = affected_by_eps[-epsilon][bin_number] - combined_by_eps[-epsilon][bin_number]
                dp = affected_by_eps[+epsilon][bin_number] - combined_by_eps[+epsilon][bin_number]
            else:
                nominal_period = float(frame_by_bin.loc[bin_number, f"ll0_{residual_period}"])
                dm = nominal_period - combined_by_eps[-epsilon][bin_number]
                dp = nominal_period - combined_by_eps[+epsilon][bin_number]
            # endif
            slope = (dp - dm) / (2.0 * epsilon)
            response[i, j] = slope
            probe_rows.append({
                "charge_ratio_parameter": labels[j], "changed_period": changed_period,
                "bin_number": bin_number, "residual_period": residual_period,
                "nominal_period_minus_combined": d0, "residual_minus": dm,
                "residual_plus": dp, "dresidual_ddelta": slope,
                "nested_variance": nested_var,
            })
        # endfor
        print(f"[target-charge-ratio] response {j+1}/3: {labels[j]}", flush=True)
    # endfor

    Aw = response * sqrt_w[:, None]
    yw = y * sqrt_w
    ata = Aw.T @ Aw
    aty = Aw.T @ yw
    covariance = np.linalg.pinv(ata, rcond=1.0e-12)
    delta = -covariance @ aty
    chi2_before = float(np.sum(yw * yw))
    linear_resid = y + response @ delta
    chi2_linear = float(np.sum((linear_resid * sqrt_w) ** 2))

    # Full-likelihood validation at the simultaneous three-parameter solution.
    # This is deliberately one exact evaluation rather than a costly nonlinear
    # optimizer over hundreds of repeated fits.  For a plausible few-percent
    # solution the +/-1% response is expected to be locally linear; the exact
    # validation quantifies any failure of that approximation.
    exact_valid = bool(np.all(np.abs(delta) < 0.25))
    exact_residuals = np.full_like(y, np.nan)
    chi2_exact = math.nan
    if exact_valid:
        distorted_solution = _copy_run_states_with_period_target_ratio_scale(
            run_states, {p: float(delta[j]) for j, p in enumerate(PERIODS)}
        )
        exact_period: dict[tuple[int, str], float] = {}
        exact_combined: dict[int, float] = {}
        for bin_number in range(1, NUMBER_OF_BINS + 1):
            nominal = frame_by_bin.loc[bin_number]
            combined_initial = {name: float(nominal[name]) for name in PHYSICS_PARAMETERS}
            combined_fixed = {name: value for name, value in combined_initial.items() if name != "ll0"}
            cf = fit_one_variant(
                events, distorted_solution, dilution_records, bin_number, "nominal",
                initial_values=combined_initial, fixed_physics_parameters=combined_fixed,
            )
            exact_combined[bin_number] = float(cf["values"]["ll0"])
            for period in PERIODS:
                period_initial = {name: float(nominal[f"{name}_{period}"]) for name in PHYSICS_PARAMETERS}
                period_fixed = {name: value for name, value in period_initial.items() if name != "ll0"}
                pf = fit_one_variant(
                    events, distorted_solution, dilution_records, bin_number, "nominal",
                    active_periods=(period,), initial_values=period_initial,
                    fixed_physics_parameters=period_fixed,
                )
                exact_period[(bin_number, period)] = float(pf["values"]["ll0"])
            # endfor
        # endfor
        for i, (bin_number, period, _, _) in enumerate(obs):
            exact_residuals[i] = exact_period[(bin_number, period)] - exact_combined[bin_number]
        # endfor
        chi2_exact = float(np.sum((exact_residuals * sqrt_w) ** 2))
    # endif

    solution_rows = []
    for j, period in enumerate(PERIODS):
        state = run_states[period]
        signs = np.sign(state["pt"]).astype(int)
        qplus = float(np.sum((state["q_plus"] + state["q_minus"])[signs > 0]))
        qminus = float(np.sum((state["q_plus"] + state["q_minus"])[signs < 0]))
        ratio0 = qplus / qminus if qminus > 0.0 else math.nan
        d = float(delta[j])
        ratio_factor = (1.0 + d) / (1.0 - d) if abs(d) < 1.0 else math.nan
        sigma = math.sqrt(max(float(covariance[j, j]), 0.0))
        solution_rows.append({
            "period": period, "charge_ratio_parameter": labels[j],
            "nominal_Q_Ptplus": qplus, "nominal_Q_Ptminus": qminus,
            "nominal_Qplus_over_Qminus": ratio0,
            "best_fit_delta": d, "best_fit_delta_percent": 100.0 * d,
            "Ptplus_charge_scale": 1.0 + d, "Ptminus_charge_scale": 1.0 - d,
            "Qplus_over_Qminus_scale_factor": ratio_factor,
            "Qplus_over_Qminus_percent_change": 100.0 * (ratio_factor - 1.0) if np.isfinite(ratio_factor) else math.nan,
            "linearized_delta_uncertainty": sigma,
            "linearized_delta_uncertainty_percent": 100.0 * sigma,
        })
    # endfor

    prior_rows = []
    for prior_sigma in (0.01, 0.02, 0.03, 0.05):
        hessian = ata + np.eye(len(PERIODS)) / (prior_sigma * prior_sigma)
        dprior = -np.linalg.solve(hessian, aty)
        data_resid = y + response @ dprior
        data_chi2 = float(np.sum((data_resid * sqrt_w) ** 2))
        prior_chi2 = float(np.sum((dprior / prior_sigma) ** 2))
        for j, period in enumerate(PERIODS):
            ratio_factor = (1.0 + dprior[j]) / (1.0 - dprior[j])
            prior_rows.append({
                "prior_delta_percent": 100.0 * prior_sigma, "period": period,
                "best_fit_delta_percent": 100.0 * float(dprior[j]),
                "Qplus_over_Qminus_percent_change": 100.0 * (ratio_factor - 1.0),
                "data_chi2": data_chi2, "prior_chi2": prior_chi2,
                "penalized_chi2": data_chi2 + prior_chi2,
                "nominal_data_chi2": chi2_before,
            })
        # endfor
    # endfor

    validation_rows = []
    for i, (bin_number, period, d0, nested_var) in enumerate(obs):
        validation_rows.append({
            "bin_number": bin_number, "period": period,
            "nominal_residual": d0, "linearized_corrected_residual": linear_resid[i],
            "exact_corrected_residual": exact_residuals[i], "nested_variance": nested_var,
        })
    # endfor

    summary = pd.DataFrame([{
        "n_residuals": len(y), "n_charge_ratio_parameters": len(PERIODS),
        "chi2_before": chi2_before, "chi2_linear_solution": chi2_linear,
        "chi2_exact_validation": chi2_exact,
        "exact_validation_performed": exact_valid,
        "max_abs_delta_percent": 100.0 * float(np.max(np.abs(delta))),
        "condition_number": float(np.linalg.cond(Aw)),
        "note": "three simultaneous antisymmetric Pt+/Pt- charge-ratio corrections; exact likelihood validation performed when all |delta|<25%",
    }])

    paths = {
        "solution": output_dir / "simultaneous_target_charge_ratio_solution.csv",
        "summary": output_dir / "simultaneous_target_charge_ratio_summary.csv",
        "response": output_dir / "simultaneous_target_charge_ratio_response.csv",
        "priors": output_dir / "simultaneous_target_charge_ratio_prior_profiles.csv",
        "validation": output_dir / "simultaneous_target_charge_ratio_validation.csv",
    }
    pd.DataFrame(solution_rows).to_csv(paths["solution"], index=False)
    summary.to_csv(paths["summary"], index=False)
    pd.DataFrame(probe_rows).to_csv(paths["response"], index=False)
    pd.DataFrame(prior_rows).to_csv(paths["priors"], index=False)
    pd.DataFrame(validation_rows).to_csv(paths["validation"], index=False)

    print(
        f"[target-charge-ratio] chi2 {chi2_before:.3f} -> {chi2_linear:.3f} (linear)",
        flush=True,
    )
    if exact_valid:
        print(f"[target-charge-ratio] exact validation chi2={chi2_exact:.3f}", flush=True)
    else:
        print("[target-charge-ratio] exact validation skipped: linear solution exceeds 25% in at least one period", flush=True)
    # endif
    for row in solution_rows:
        print(
            f"[target-charge-ratio] {row['period']}: delta={row['best_fit_delta_percent']:+.3f}% "
            f"=> Q(Pt+)/Q(Pt-) change={row['Qplus_over_Qminus_percent_change']:+.3f}%",
            flush=True,
        )
    # endfor
    return {key: str(value) for key, value in paths.items()}


def write_beam_polarization_response_diagnostic(
    *, frame: pd.DataFrame, events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path, epsilon: float = 0.01,
) -> dict[str, str]:
    """Test whether a period beam-polarization scale can explain A_LL.

    The affected period beam polarization is shifted by +/-1%.  The affected
    period and simultaneous fits are both rerun, so the simultaneous reference
    is allowed to move.  LU, LL and LLcosphi are floated together because all
    three contain P_b; UU and UL amplitudes are fixed.  The summary compares the
    beam-polarization correction inferred from A_LL with that inferred from LU.
    """
    ensure_directory(output_dir)
    frame_by_bin = frame.set_index("bin_number")
    parameters = ("lu1", "ll0", "ll1")
    rows: list[dict[str, Any]] = []

    for period in PERIODS:
        for bin_number in range(1, NUMBER_OF_BINS + 1):
            nominal = frame_by_bin.loc[bin_number]
            period_initial = {name: float(nominal[f"{name}_{period}"]) for name in PHYSICS_PARAMETERS}
            combined_initial = {name: float(nominal[name]) for name in PHYSICS_PARAMETERS}
            period_fixed = {
                name: value for name, value in period_initial.items() if name not in parameters
            }
            combined_fixed = {
                name: value for name, value in combined_initial.items() if name not in parameters
            }
            fits: dict[float, tuple[dict[str, Any], dict[str, Any]]] = {}
            for e in (-epsilon, +epsilon):
                scales = {period: 1.0 + e}
                pf = fit_one_variant(
                    events, run_states, dilution_records, bin_number, "nominal",
                    active_periods=(period,), initial_values=period_initial,
                    fixed_physics_parameters=period_fixed,
                    beam_polarization_scales=scales,
                )
                cf = fit_one_variant(
                    events, run_states, dilution_records, bin_number, "nominal",
                    initial_values=combined_initial,
                    fixed_physics_parameters=combined_fixed,
                    beam_polarization_scales=scales,
                )
                fits[e] = (pf, cf)
            # endfor
            for parameter in parameters:
                dm = float(fits[-epsilon][0]["values"][parameter] - fits[-epsilon][1]["values"][parameter])
                dp = float(fits[+epsilon][0]["values"][parameter] - fits[+epsilon][1]["values"][parameter])
                slope = (dp - dm) / (2.0 * epsilon)
                d0 = float(nominal[f"{parameter}_{period}"] - nominal[parameter])
                required = -d0 / slope if abs(slope) > 1.0e-12 else math.nan
                pstat = float(nominal[f"{parameter}_stat_{period}"])
                cstat = float(nominal[f"{parameter}_stat"])
                rows.append({
                    "period": period, "bin_number": bin_number, "parameter": parameter,
                    "nominal_beam_polarization": BEAM_POLARIZATION[period],
                    "epsilon_probe": epsilon, "nominal_period_minus_combined": d0,
                    "d_period_minus_combined_depsilon": slope,
                    "required_fractional_beam_polarization_correction": required,
                    "required_percent_beam_polarization_correction": 100.0 * required if np.isfinite(required) else math.nan,
                    "nested_variance": max(pstat*pstat - cstat*cstat, 0.0),
                })
            # endfor
        # endfor
        print(f"[beam-polarization response] {period} complete", flush=True)
    # endfor

    result = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for (period, parameter), group in result.groupby(["period", "parameter"]):
        num = den = 0.0
        for _, r in group.iterrows():
            var = float(r["nested_variance"])
            slope = float(r["d_period_minus_combined_depsilon"])
            d0 = float(r["nominal_period_minus_combined"])
            if var <= 0.0 or not np.isfinite(slope):
                continue
            # endif
            w = 1.0 / var
            num += w * slope * d0
            den += w * slope * slope
        # endfor
        req = -num / den if den > 0.0 else math.nan
        summary_rows.append({
            "period": period, "parameter": parameter,
            "required_fractional_beam_polarization_correction": req,
            "required_percent_beam_polarization_correction": 100.0 * req if np.isfinite(req) else math.nan,
            "corrected_beam_polarization": BEAM_POLARIZATION[period] * (1.0 + req) if np.isfinite(req) else math.nan,
        })
    # endfor
    summary = pd.DataFrame(summary_rows)
    # Put the LU and LL answers next to one another for the key consistency test.
    pivot = summary.pivot(index="period", columns="parameter", values="required_percent_beam_polarization_correction").reset_index()
    pivot = pivot.rename(columns={
        "lu1": "required_percent_from_LU_sinphi",
        "ll0": "required_percent_from_ALL",
        "ll1": "required_percent_from_ALL_cosphi",
    })
    summary_path = output_dir / "beam_polarization_required_corrections.csv"
    bins_path = output_dir / "beam_polarization_response_by_bin.csv"
    result.to_csv(bins_path, index=False)
    pivot.to_csv(summary_path, index=False)
    return {"bins": str(bins_path), "summary": str(summary_path)}


def write_flagged_epoch_refits(
    *, frame: pd.DataFrame, events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path, threshold: float = 2.5,
) -> str:
    """Fast target-polarization-epoch diagnostic for the constant A_LL term.

    The previous implementation floated all five polarized amplitudes in every
    flagged-bin/epoch fit.  Small target epochs can make those fits extremely
    slow and poorly constrained.  Here the diagnostic is deliberately surgical:
    every contiguous target-sign epoch is fit in every kinematic bin with ll0
    free while u1, u2 and the other four polarized amplitudes are fixed to the
    corresponding full-period solution.  Beam helicity still flips within an
    epoch, so ll0 remains directly identifiable.  This answers whether the
    period-wide ll0 displacement is present in both target signs or is localized
    to a particular target epoch without asking low-statistics subsets to
    redetermine unrelated harmonics.
    """
    ensure_directory(output_dir)
    frame_by_bin = frame.set_index("bin_number")

    # Build the same contiguous target-sign epochs used by the charge audit.
    epoch_defs: dict[str, list[tuple[int, int, int, int]]] = {}
    luminosity_rows: list[dict[str, Any]] = []
    for period in PERIODS:
        state = run_states[period]
        runs = np.asarray(state["run"], dtype=int)
        pt = np.asarray(state["pt"], dtype=float)
        qp = np.asarray(state["q_plus"], dtype=float)
        qm = np.asarray(state["q_minus"], dtype=float)
        sign_pt = np.sign(pt).astype(int)
        epochs: list[tuple[int, int, int, int]] = []
        if len(runs):
            start = 0
            epoch_index = 0
            for j in range(1, len(runs) + 1):
                boundary = j == len(runs) or sign_pt[j] != sign_pt[j - 1]
                if not boundary:
                    continue
                # endif

                epoch_index += 1
                sl = slice(start, j)
                target_sign = int(sign_pt[start])
                q_plus = float(np.sum(qp[sl]))
                q_minus = float(np.sum(qm[sl]))
                q_total = q_plus + q_minus
                # For a fixed target sign, h*s=+ corresponds to h=s.
                q_hs_plus = q_plus if target_sign > 0 else q_minus
                q_hs_minus = q_minus if target_sign > 0 else q_plus
                hs_asym = (
                    (q_hs_plus - q_hs_minus) / q_total
                    if q_total > 0.0 else math.nan
                )
                luminosity_rows.append({
                    "period": period, "epoch": epoch_index,
                    "run_low": int(runs[start]), "run_high": int(runs[j - 1]),
                    "target_sign": target_sign, "number_of_runs": int(j - start),
                    "charge_h_plus": q_plus, "charge_h_minus": q_minus,
                    "charge_hs_plus": q_hs_plus, "charge_hs_minus": q_hs_minus,
                    "charge_total": q_total,
                    "helicity_charge_asymmetry": (
                        (q_plus - q_minus) / q_total if q_total > 0.0 else math.nan
                    ),
                    "hs_charge_asymmetry": hs_asym,
                    "charge_weighted_mean_target_polarization": (
                        float(np.sum(pt[sl] * (qp[sl] + qm[sl])) / q_total)
                        if q_total > 0.0 else math.nan
                    ),
                })
                epochs.append((epoch_index, int(runs[start]), int(runs[j - 1]), target_sign))
                start = j
            # endfor
        # endif
        epoch_defs[period] = epochs
    # endfor

    luminosity_path = output_dir / "target_epoch_four_state_luminosity.csv"
    pd.DataFrame(luminosity_rows).to_csv(luminosity_path, index=False)

    # Fit only ll0 in each epoch.  Everything else is anchored to the existing
    # full-period solution for that same bin, making this both fast and directly
    # targeted at the observed anomaly.
    rows: list[dict[str, Any]] = []
    for period in PERIODS:
        for epoch_index, run_low, run_high, target_sign in epoch_defs[period]:
            for bin_number in range(1, NUMBER_OF_BINS + 1):
                nominal_row = frame_by_bin.loc[bin_number]
                mask = (
                    (events["bin_number"] == bin_number)
                    & (events["period_index"] == PERIOD_INDEX[period])
                    & (events["runnum"] >= run_low)
                    & (events["runnum"] <= run_high)
                )
                n_events = int(np.count_nonzero(mask))
                if n_events < 50:
                    continue
                # endif

                fixed = {"u1": 0.0, "u2": 0.0}
                for parameter in ("lu1", "ul1", "ul2", "ll1"):
                    fixed[parameter] = float(nominal_row[f"{parameter}_{period}"])
                # endfor
                initial = {
                    name: float(nominal_row[f"{name}_{period}"])
                    for name in PUBLISHED_SYSTEMATIC_PARAMETERS
                }
                try:
                    fit = fit_one_variant(
                        events, run_states, dilution_records, bin_number, "nominal",
                        active_periods=(period,), initial_values=initial,
                        fixed_physics_parameters=fixed,
                        run_ranges_filter=((run_low, run_high),),
                    )
                    ll0 = float(fit["values"]["ll0"])
                    ll0_stat = float(fit["errors"]["ll0"])
                    valid = bool(fit["valid"])
                    edm = float(fit["edm"])
                    error = ""
                except Exception as exc:
                    ll0 = ll0_stat = edm = math.nan
                    valid = False
                    error = str(exc)
                # endtry

                period_ll0 = float(nominal_row[f"ll0_{period}"])
                period_stat = float(nominal_row[f"ll0_stat_{period}"])
                combined_ll0 = float(nominal_row["ll0"])
                combined_stat = float(nominal_row["ll0_stat"])
                rows.append({
                    "period": period, "epoch": epoch_index,
                    "run_low": run_low, "run_high": run_high,
                    "target_sign": target_sign, "bin_number": bin_number,
                    "events": n_events, "valid": valid, "edm": edm, "error": error,
                    "ll0": ll0, "ll0_stat": ll0_stat,
                    "period_ll0": period_ll0, "period_ll0_stat": period_stat,
                    "combined_ll0": combined_ll0, "combined_ll0_stat": combined_stat,
                    "epoch_minus_period": ll0 - period_ll0 if np.isfinite(ll0) else math.nan,
                    "epoch_minus_combined": ll0 - combined_ll0 if np.isfinite(ll0) else math.nan,
                })
            # endfor
            print(
                f"[epoch ll0] {period} epoch {epoch_index}: runs {run_low}-{run_high}, "
                f"target sign {target_sign:+d} complete",
                flush=True,
            )
        # endfor
    # endfor

    path = output_dir / "target_epoch_ll0_refits.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return str(path)

def plot_period_stability(
    frame: pd.DataFrame,
    output_dir: Path,
) -> list[str]:
    """
    Compare fully free period fits with diagnostic fits fixing u1 only,
    u2 only, or both u1 and u2 to the simultaneous values.
    """
    ensure_directory(output_dir)
    paths: list[str] = []

    constraint_specs = (
        ("free", "Fully free period fit", "o", None),
        ("fix_u1", r"$u_1$ fixed", "s", "tab:purple"),
        ("fix_u2", r"$u_2$ fixed", "^", "tab:brown"),
        ("fix_u1_u2", r"$u_1,u_2$ fixed", "D", COMBINED_COLOR),
    )

    for period in PERIODS:
        for x_index, (x_low, x_high) in enumerate(XB_BINS):
            subset = frame.loc[
                frame["x_index"] == x_index
            ].sort_values("t_index")
            x_values = subset["mean_minus_tprime_gev2"].to_numpy(dtype=float)
            fig, axes_by_parameter, axes_in_order = make_grouped_physics_canvas()

            for parameter in AGGREGATED_PANEL_ORDER:
                ax = axes_by_parameter[parameter]

                for constraint_name, label, marker, color in constraint_specs:
                    if constraint_name == "free":
                        values = subset[f"{parameter}_{period}"]
                        errors = subset[f"{parameter}_stat_{period}"]
                        draw_color = PERIOD_COLORS[period]
                    else:
                        values = subset[
                            f"{parameter}_{constraint_name}_{period}"
                        ]
                        errors = subset[
                            f"{parameter}_{constraint_name}_stat_{period}"
                        ]
                        draw_color = color
                    # endif

                    ax.errorbar(
                        x_values,
                        values,
                        yerr=errors,
                        marker=marker,
                        linestyle="none",
                        capsize=2,
                        color=draw_color,
                        label=label,
                    )
                # endfor

                ax.axhline(0.0, linewidth=0.8)
                ax.set_ylabel(PARAMETER_LABELS[parameter])
                apply_parameter_y_limits(ax, parameter)
                ax.grid(alpha=0.25)
                if parameter in PHYSICS_PANEL_ROWS[-1][1]:
                    ax.set_xlabel(r"$-t^\\prime$ (GeV$^2$)")
                # endif
            # endfor

            handles, labels = axes_in_order[0].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="lower right",
                bbox_to_anchor=(0.97, 0.055),
            )
            fig.suptitle(
                rf"{PERIOD_LABELS[period]}, "
                rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: fit stability",
                y=0.995,
            )
            fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.98))
            path = (
                output_dir
                / f"{period}_xB_bin_{x_index + 1}_stability.png"
            )
            fig.savefig(path, dpi=180)
            plt.close(fig)
            paths.append(str(path))
        # endfor
    # endfor
    return paths


def plot_period_consistency_heatmap(
    frame: pd.DataFrame,
    output_dir: Path,
) -> str:
    """Plot Delta NLL for each period and combined kinematic bin."""
    ensure_directory(output_dir)
    matrix = np.asarray(
        [
            frame[f"period_delta_nll_{period}"].to_numpy(dtype=float)
            for period in PERIODS
        ],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(15, 3.8))
    image = ax.imshow(matrix, aspect="auto", origin="upper")
    ax.set_xticks(np.arange(NUMBER_OF_BINS))
    ax.set_xticklabels(frame["bin_number"].astype(int).tolist())
    ax.set_yticks(np.arange(len(PERIODS)))
    ax.set_yticklabels([PERIOD_LABELS[period] for period in PERIODS])
    ax.set_xlabel("Combined kinematic-bin number")
    ax.set_ylabel("Run period")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label(
        r"$\Delta\mathrm{NLL}_p="
        r"\mathrm{NLL}_p(\widehat{\theta}_{\mathrm{combined}})"
        r"-\mathrm{NLL}_{p,\min}$"
    )
    fig.tight_layout()
    path = output_dir / "period_consistency_delta_nll.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)

def write_latex_table(
    frame: pd.DataFrame,
    path: Path,
    *,
    include_target_axis_uncertainty: bool,
    sample_variant: str,
) -> None:
    """Write the publication-facing table for the five polarized observables.

    The two UU harmonics remain fitted nuisance quantities and are deliberately
    omitted.  Once the final point-to-point columns have been attached, the
    nominal table quotes statistical and total point-to-point uncertainties.
    Earlier/intermediate variant tables quote the uncertainty information that
    is actually available at that stage.
    """
    ensure_directory(path.parent)
    parameters = PUBLISHED_SYSTEMATIC_PARAMETERS
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\small",
        r"\begin{tabular}{rrrrrr}",
        r"\hline",
        (
            r"Bin & $F_{LU}^{\sin\phi}/F_{UU}$ & "
            r"$F_{UL}^{\sin\phi}/F_{UU}$ & "
            r"$F_{UL}^{\sin2\phi}/F_{UU}$ & "
            r"$F_{LL}/F_{UU}$ & "
            r"$F_{LL}^{\cos\phi}/F_{UU}$ \\"
        ),
        r"\hline",
    ]

    has_final_ptp = all(
        f"{parameter}_point_to_point_systematic" in frame.columns
        for parameter in parameters
    )
    for row in frame.itertuples(index=False):
        entries = [str(int(row.bin_number))]
        for parameter in parameters:
            value = float(getattr(row, parameter))
            stat = float(getattr(row, f"{parameter}_stat"))
            if has_final_ptp:
                ptp = float(getattr(row, f"{parameter}_point_to_point_systematic"))
                entries.append(rf"${value:.5f}\pm{stat:.5f}\pm{ptp:.5f}$")
            elif include_target_axis_uncertainty:
                axis = float(getattr(row, f"{parameter}_target_axis_study_envelope"))
                entries.append(rf"${value:.5f}\pm{stat:.5f}\pm{axis:.5f}$")
            else:
                entries.append(rf"${value:.5f}\pm{stat:.5f}$")
            # endif
        # endfor
        lines.append(" & ".join(entries) + r" \\")
    # endfor

    if has_final_ptp:
        caption = (
            r"Nominal simultaneous unbinned-likelihood results for the five "
            r"published polarized structure-function ratios. The first "
            r"uncertainty is statistical and includes the Gaussian-constrained "
            r"dilution-factor statistical uncertainty. The second is the total "
            r"point-to-point systematic uncertainty, combining radiation, "
            r"target-axis treatment, channel selection, and bin migration in "
            r"quadrature. Correlated polarization and dilution-model scale "
            r"uncertainties are not included."
        )
    elif include_target_axis_uncertainty:
        caption = (
            r"Intermediate nominal simultaneous unbinned-likelihood results "
            r"for the five published polarized structure-function ratios. "
            r"The first uncertainty is statistical and the second is the "
            r"target-axis treatment uncertainty."
        )
    else:
        caption = (
            rf"{sample_variant.upper()} simultaneous unbinned-likelihood "
            r"diagnostic results for the five polarized structure-function "
            r"ratios. The quoted uncertainty is statistical."
        )
    # endif

    lines.extend([
        r"\hline",
        r"\end{tabular}",
        rf"\caption{{{caption}}}",
        r"\end{table}",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")



def write_external_transverse_input_tables(
    csv_path: Path,
    json_path: Path,
) -> None:
    """Write the fixed external leakage inputs and full provenance."""
    rows: list[dict[str, Any]] = []
    for term, payload in EXTERNAL_TRANSVERSE_INPUTS.items():
        reference = payload["reference"]
        kinematics = payload["kinematic_selection"]
        rows.append(
            {
                "term": term,
                "input_value": payload["value"],
                "measured_central_value": payload[
                    "measured_central_value"
                ],
                "scale_factor": payload["scale_factor"],
                "observable": payload["observable"],
                "channel": payload["channel"],
                "kinematic_selection_json": json.dumps(
                    kinematics,
                    sort_keys=True,
                ),
                "citation": reference["citation"],
                "title": reference["title"],
                "doi": reference["doi"],
                "arxiv": reference["arxiv"],
                "figure": reference.get("figure"),
                "data_note": reference["data_note"],
            }
        )
    # endfor
    ensure_directory(csv_path.parent)
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    write_json(
        json_path,
        {
            "policy": (
                "The UT inputs are direct inverse-total-variance weighted "
                "averages of exclusive-pi+ measurements. The LT inputs are "
                "rounded estimates from the open highest-z SIDIS pi+ points. "
                "No factor of two is applied."
            ),
            "inputs": EXTERNAL_TRANSVERSE_INPUTS,
        },
    )



# =============================================================================
# Command line
# =============================================================================

def resolve_isr_cut_json(explicit_path: Path | None) -> Path:
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Explicit ISR cut JSON does not exist: {path}")
        return path
    # endif
    manifest = DEFAULT_CHANNEL_SELECTION_MANIFEST.expanduser().resolve()
    if manifest.is_file():
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        retained = payload.get("retained_analysis_products", {})
        isr_root = retained.get("isr")
        if isr_root:
            candidate = Path(isr_root) / "final_carbon_assisted_cuts/tables/final_carbon_assisted_mx2_cuts.json"
            if candidate.is_file():
                return candidate.resolve()
            # endif
        # endif
    # endif
    path = DEFAULT_ISR_CUT_JSON.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Could not resolve ISR channel-selection JSON: {path}")
    return path


def run_analysis_variant(
    *,
    sample_variant: str,
    input_paths: dict[str, Path],
    run_info_path: Path,
    cut_json_path: Path,
    dilution_json_path: Path,
    output_dir: Path,
    cache_path: Path,
    tree_name: str,
    chunk_size: str,
    workers: int,
    reuse_cache: bool,
    skip_plots: bool,
    include_target_axis_study: bool,
    include_period_diagnostics: bool = False,
    cut_label: str = "nominal",
    source_cache_path: Path | None = None,
    zero_uu_baseline: bool = False,
) -> dict[str, Any]:
    tables_dir = output_dir / "tables"
    json_dir = output_dir / "json"
    covariance_dir = output_dir / "covariance"
    plots_dir = output_dir / "plots"
    all_bins_plots_dir = plots_dir / "all_bins"
    aggregated_plots_dir = plots_dir / "aggregated"
    period_plots_dir = aggregated_plots_dir / "by_period"
    target_axis_variants_dir = plots_dir / "target_axis_variants"
    period_stability_dir = plots_dir / "period_stability"
    latex_dir = output_dir / "latex"
    for directory in (
        output_dir, tables_dir, json_dir, covariance_dir, plots_dir,
        all_bins_plots_dir, aggregated_plots_dir, period_plots_dir,
        period_stability_dir, latex_dir,
    ):
        ensure_directory(directory)
    # endfor
    if include_target_axis_study:
        ensure_directory(target_axis_variants_dir)
    # endif

    print("-" * 78)
    print(f"Starting {sample_variant} structure-function-ratio analysis")
    print("-" * 78)
    print(f"Channel cuts:         {cut_json_path}")
    print(f"Dilution factors:     {dilution_json_path}")
    print(f"Output directory:     {output_dir}")
    print(f"Selected-event cache: {cache_path}")
    print(f"Target-axis study:    {include_target_axis_study}")
    print(f"Period diagnostics:   {include_period_diagnostics}")
    print(f"Exclusivity window:   {cut_label}")

    run_records = parse_run_info_csv(run_info_path)
    run_states = run_state_arrays(run_records)
    cuts = load_channel_cuts(cut_json_path, cut_label=cut_label)
    dilution_records = load_dilution_factors(
        dilution_json_path, cut_label=cut_label
    )
    if reuse_cache:
        print(f"[cache] REUSE requested: {cache_path}", flush=True)
        events = load_event_cache(cache_path)
        print(
            f"[cache] REUSE DONE: loaded {events['runnum'].size:,} events",
            flush=True,
        )
        cache_summary = {
            "cache_path": str(cache_path),
            "number_of_selected_events": int(events["runnum"].size),
            "reused": True,
            "derived_from_cache": None,
        }
    elif source_cache_path is not None:
        cache_summary = derive_event_cache(
            source_cache_path=source_cache_path,
            cuts=cuts,
            cache_path=cache_path,
        )
        cache_summary["reused"] = False
        events = load_event_cache(cache_path)
    else:
        cache_summary = build_event_cache(
            input_paths=input_paths,
            tree_name=tree_name,
            chunk_size=chunk_size,
            run_records=run_records,
            cuts=cuts,
            cache_path=cache_path,
        )
        cache_summary["reused"] = False
        cache_summary["derived_from_cache"] = None
        events = load_event_cache(cache_path)
    # endif

    run_state_payload = {period: {key: values.tolist() for key, values in state.items()} for period, state in run_states.items()}
    dilution_payload = {period: {str(bin_number): {"x_index": record.x_index, "t_index": record.t_index, "value": record.value, "stat_uncertainty": record.stat_uncertainty} for (record_period, bin_number), record in dilution_records.items() if record_period == period} for period in PERIODS}
    results: list[dict[str, Any]] = []

    # Stage 1: the quoted simultaneous extraction.  These are exactly the
    # original nominal fit_one_variant calls, now returned and checkpointed
    # before any diagnostic minimizations begin.
    nominal_task_kind = "nominal_zero_uu" if zero_uu_baseline else "nominal"
    nominal_tasks = [
        {"kind": nominal_task_kind, "bin_number": bin_number}
        for bin_number in range(1, NUMBER_OF_BINS + 1)
    ]
    nominal_stage = run_fit_stage(
        stage_name=f"{sample_variant}: simultaneous",
        tasks=nominal_tasks,
        workers=workers,
        cache_path=cache_path,
        run_state_payload=run_state_payload,
        dilution_payload=dilution_payload,
        output_dir=output_dir, checkpoint_tag="01_simultaneous_tasks",
    )
    nominal_by_bin = {
        int(item["bin_number"]): item["fit"] for item in nominal_stage
    }
    for bin_number in range(1, NUMBER_OF_BINS + 1):
        results.append(
            initialize_result_shell(
                bin_number,
                nominal_by_bin[bin_number],
                include_target_axis_study,
                include_period_diagnostics,
            )
        )
    # endfor
    results_by_bin = {item["bin_number"]: item for item in results}
    write_stage_checkpoint(output_dir, "01_simultaneous", results)

    invalid_nominal = [
        bin_number
        for bin_number, fit in sorted(nominal_by_bin.items())
        if not bool(fit["valid"])
    ]
    if invalid_nominal:
        print(
            f"[{sample_variant}] WARNING: simultaneous fit invalid in bins "
            f"{invalid_nominal}. No automatic retry or fit-setting change "
            "has been applied.",
            flush=True,
        )
    # endif

    # Stage 2: target-axis alternatives.  The calls and starting values are
    # identical to the original bundled worker.
    if include_target_axis_study:
        target_tasks: list[dict[str, Any]] = []
        for bin_number in range(1, NUMBER_OF_BINS + 1):
            nominal = nominal_by_bin[bin_number]
            for kind in ("photon_axis_projection", "external_data_informed"):
                target_tasks.append({
                    "kind": kind,
                    "bin_number": bin_number,
                    "nominal": nominal,
                })
            # endfor
        # endfor
        target_stage = run_fit_stage(
            stage_name=f"{sample_variant}: target-axis",
            tasks=target_tasks,
            workers=workers,
            cache_path=cache_path,
            run_state_payload=run_state_payload,
            dilution_payload=dilution_payload,
            output_dir=output_dir, checkpoint_tag="02_target_axis_tasks",
        )
        for item in target_stage:
            result = results_by_bin[int(item["bin_number"])]
            result["variants"][item["kind"]] = item["fit"]
        # endfor
        for result in results:
            finish_target_axis_study_envelopes(result)
        # endfor
        write_stage_checkpoint(output_dir, "02_target_axis", results)
    # endif

    # Stages 3 and 4 are nominal-only diagnostics.  They remain complete final
    # outputs, but no longer block saving the production simultaneous fits.
    if include_period_diagnostics:
        period_tasks: list[dict[str, Any]] = []
        for bin_number in range(1, NUMBER_OF_BINS + 1):
            nominal = nominal_by_bin[bin_number]
            for period in PERIODS:
                period_tasks.append({
                    "kind": (
                        "period_only_zero_uu" if zero_uu_baseline else "period_only"
                    ),
                    "bin_number": bin_number,
                    "period": period,
                    "nominal": nominal,
                })
            # endfor
        # endfor
        # Evaluation-only preflight: inspect state populations and local NLL
        # sensitivity before launching the expensive period-only minimizations.
        # This is diagnostic only and does not alter any fit.
        preflight_tasks = [
            {"bin_number": task["bin_number"], "period": task["period"],
             "nominal": task["nominal"]}
            for task in period_tasks
        ]
        run_period_preflight_stage(
            tasks=preflight_tasks, workers=workers, cache_path=cache_path,
            run_state_payload=run_state_payload,
            dilution_payload=dilution_payload, output_dir=output_dir,
        )

        period_stage = run_fit_stage(
            stage_name=f"{sample_variant}: period-only",
            tasks=period_tasks,
            workers=workers,
            cache_path=cache_path,
            run_state_payload=run_state_payload,
            dilution_payload=dilution_payload,
            output_dir=output_dir, checkpoint_tag="03_period_only_tasks",
        )
        for item in period_stage:
            results_by_bin[int(item["bin_number"])]["period_fits"][
                item["period"]
            ] = item["fit"]
        # endfor
        write_stage_checkpoint(output_dir, "03_period_only", results)

        constraint_tasks: list[dict[str, Any]] = []
        if zero_uu_baseline:
            # In the dedicated zero-UU period-stability cross-check, u1 and u2
            # are already fixed to zero in every period-only fit.  The older
            # constraint variants would therefore be redundant minimizations.
            # Populate their compatibility slots with the same fit so the
            # generic table writer remains backward compatible.
            for bin_number in range(1, NUMBER_OF_BINS + 1):
                result = results_by_bin[bin_number]
                for period in PERIODS:
                    period_fit = result["period_fits"][period]
                    for constraint in ("fix_u1", "fix_u2", "fix_u1_u2"):
                        result["period_constraint_fits"][constraint][period] = period_fit
                    # endfor
                # endfor
            # endfor
        else:
            for bin_number in range(1, NUMBER_OF_BINS + 1):
                result = results_by_bin[bin_number]
                nominal = nominal_by_bin[bin_number]
                for period in PERIODS:
                    period_fit = result["period_fits"][period]
                    for constraint in ("fix_u1", "fix_u2", "fix_u1_u2"):
                        constraint_tasks.append({
                            "kind": "period_constraint",
                            "bin_number": bin_number,
                            "period": period,
                            "constraint": constraint,
                            "nominal": nominal,
                            "period_fit": period_fit,
                        })
                    # endfor
                # endfor
            # endfor
            constraint_stage = run_fit_stage(
                stage_name=f"{sample_variant}: period-constraints",
                tasks=constraint_tasks,
                workers=workers,
                cache_path=cache_path,
                run_state_payload=run_state_payload,
                dilution_payload=dilution_payload,
                output_dir=output_dir, checkpoint_tag="04_period_constraints_tasks",
            )
            for item in constraint_stage:
                results_by_bin[int(item["bin_number"])][
                    "period_constraint_fits"
                ][item["constraint"]][item["period"]] = item["fit"]
            # endfor
        # endif

        print(
            f"[{sample_variant}] computing period-consistency NLL diagnostics...",
            flush=True,
        )
        for result in results:
            finish_period_consistency(
                result, events, run_states, dilution_records
            )
        # endfor
        write_stage_checkpoint(output_dir, "04_period_constraints", results)
    # endif

    print(
        f"[{sample_variant}] all requested fit stages complete; "
        "assembling final tables and plots.",
        flush=True,
    )
    results.sort(key=lambda item: item["bin_number"])
    frame = flatten_fit_results(results)
    csv_path = tables_dir / "structure_function_ratios.csv"
    frame.to_csv(csv_path, index=False)
    detailed_json_path = json_dir / "structure_function_ratios.json"
    write_json(detailed_json_path, {
        "schema_version": 2,
        "analysis": "RGC exclusive enpi+ structure-function-ratio extraction",
        "sample_variant": sample_variant,
        "diagnostic_only": sample_variant != "nominal",
        "exclusivity_window": cut_label,
        "target_axis_study_performed": include_target_axis_study,
        "period_diagnostics_performed": include_period_diagnostics,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "beam_polarization": BEAM_POLARIZATION,
        "beam_energy_gev": BEAM_ENERGY_GEV,
        "cache": cache_summary,
        "inputs": {
            "root_files": {period: str(input_paths[period].expanduser().resolve()) for period in PERIODS},
            "run_info_csv": str(run_info_path), "run_info_csv_sha256": sha256_file(run_info_path),
            "channel_cut_json": str(cut_json_path), "channel_cut_json_sha256": sha256_file(cut_json_path),
            "dilution_json": str(dilution_json_path), "dilution_json_sha256": sha256_file(dilution_json_path),
        },
        "results": results,
    })
    for result in results:
        nominal = result["variants"]["nominal"]
        if nominal["covariance"] is not None:
            b = result["bin_number"]
            np.save(covariance_dir / f"bin_{b:02d}_nominal_covariance.npy", np.asarray(nominal["covariance"], dtype=np.float64))
            np.save(covariance_dir / f"bin_{b:02d}_nominal_correlation.npy", np.asarray(nominal["correlation"], dtype=np.float64))
        # endif
    # endfor
    latex_path = latex_dir / "structure_function_ratios.tex"
    write_latex_table(
        frame,
        latex_path,
        include_target_axis_uncertainty=include_target_axis_study,
        sample_variant=sample_variant,
    )
    plot_paths = {"all_bins": [], "aggregated": [], "aggregated_by_period": [], "target_axis_variants": [], "period_stability": []}
    if not skip_plots:
        print(f"[{sample_variant}] writing plots...", flush=True)
        plot_paths["all_bins"] = plot_parameter_summaries(
            frame,
            all_bins_plots_dir,
            include_target_axis_uncertainty=include_target_axis_study,
        )
        plot_paths["aggregated"] = plot_aggregated_by_x(
            frame,
            aggregated_plots_dir,
            include_target_axis_uncertainty=include_target_axis_study,
        )
        if include_period_diagnostics:
            print(f"[{sample_variant}] writing period-consistency diagnostics...", flush=True)
            plot_paths["aggregated_by_period"] = plot_aggregated_by_period(frame, period_plots_dir)
            plot_paths["aggregated_by_period"].append(
                plot_period_consistency_heatmap(frame, period_plots_dir)
            )
            plot_paths["period_stability"] = plot_period_stability(
                frame, period_stability_dir
            )
        # endif
        if include_target_axis_study:
            plot_paths["target_axis_variants"] = plot_target_axis_variants(frame, target_axis_variants_dir)
        # endif
    # endif
    manifest_path = output_dir / "analysis_variant_manifest.json"
    write_json(manifest_path, {"schema_version": 3, "sample_variant": sample_variant, "diagnostic_only": sample_variant != "nominal", "exclusivity_window": cut_label, "target_axis_study_performed": include_target_axis_study, "period_diagnostics_performed": include_period_diagnostics, "products": {"csv": str(csv_path), "detailed_json": str(detailed_json_path), "latex": str(latex_path), "covariance_directory": str(covariance_dir), "plots": plot_paths, "cache": str(cache_path)}})
    invalid_bins = frame.loc[~frame["nominal_fit_valid"].astype(bool), "bin_number"].astype(int).tolist()
    return {"sample_variant": sample_variant, "frame": frame, "results": results, "events": int(events["runnum"].size), "csv": str(csv_path), "json": str(detailed_json_path), "latex": str(latex_path), "manifest": str(manifest_path), "invalid_bins": invalid_bins}


def write_nominal_isr_comparison_products(
    nominal: pd.DataFrame,
    isr: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Write nominal-versus-ISR products using the common layout."""
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    covariance_dir = output_dir / "covariance"
    summary_dir = output_dir / "summary"
    for directory in (tables_dir, plots_dir, covariance_dir, summary_dir):
        ensure_directory(directory)
    # endfor

    keys = ["bin_number", "x_index", "t_index"]
    keep = keys + [
        "min_xB", "median_xB", "mean_xB", "max_xB", "min_Q2_gev2", "median_Q2_gev2", "mean_Q2_gev2", "max_Q2_gev2", "min_W_gev", "median_W_gev", "mean_W_gev", "max_W_gev", "min_minus_t_gev2", "median_minus_t_gev2", "mean_minus_t_gev2", "max_minus_t_gev2", "min_minus_tprime_gev2", "median_minus_tprime_gev2", "mean_minus_tprime_gev2", "max_minus_tprime_gev2", "min_epsilon", "median_epsilon", "mean_epsilon", "max_epsilon", "min_DepA", "median_DepA", "mean_DepA", "max_DepA", "min_DepB", "median_DepB", "mean_DepB", "max_DepB", "min_DepC", "median_DepC", "mean_DepC", "max_DepC", "min_DepV", "median_DepV", "mean_DepV", "max_DepV", "min_DepW", "median_DepW", "mean_DepW", "max_DepW",
        "number_of_events",
    ] + [
        item for parameter in PHYSICS_PARAMETERS
        for item in (parameter, f"{parameter}_stat")
    ]
    merged = nominal[keep].merge(
        isr[keep], on=keys, suffixes=("_nominal", "_isr"),
        validate="one_to_one",
    )

    barlow_records: list[dict[str, Any]] = []
    covariance_products: dict[str, str] = {}
    for parameter in PHYSICS_PARAMETERS:
        shift_column = f"{parameter}_isr_minus_nominal"
        systematic_column = f"{parameter}_absolute_isr_systematic"
        merged[shift_column] = (
            merged[f"{parameter}_isr"] - merged[f"{parameter}_nominal"]
        )
        merged[systematic_column] = merged[shift_column].abs()
        merged[f"{parameter}_isr_systematic_to_nominal_stat"] = np.divide(
            merged[systematic_column], merged[f"{parameter}_stat_nominal"],
        )
        barlow = calculate_barlow_arrays(
            merged[f"{parameter}_nominal"], merged[f"{parameter}_isr"],
            merged[f"{parameter}_stat_nominal"], merged[f"{parameter}_stat_isr"],
        )
        prefix = f"{parameter}_isr_barlow"
        merged[f"{prefix}_statistical_scale"] = barlow["denominator"]
        merged[f"{prefix}_value"] = barlow["barlow"]
        merged[f"{prefix}_passes"] = barlow["passed"]
        merged[f"{prefix}_status"] = barlow["status"]
        barlow_records.append(summarize_barlow_status(
            merged, effect="ISR", variation="isr_minus_nominal",
            parameter=parameter, barlow_column=f"{prefix}_value",
            status_column=f"{prefix}_status",
        ))
        shift = merged[shift_column].to_numpy(dtype=np.float64)
        covariance_path = covariance_dir / f"{parameter}_isr_covariance.npy"
        np.save(covariance_path, np.outer(shift, shift))
        covariance_products[parameter] = str(covariance_path)
    # endfor

    x = merged["bin_number"].to_numpy(dtype=float)
    definition = r"$\delta_{\rm ISR}=|a_{\rm ISR}-a_{\rm nominal}|$"

    def draw(parameter: str, axes: np.ndarray, *, show_legends: bool, show_xlabel: bool) -> None:
        draw_systematic_comparison(
            axes, x=x, parameter=parameter,
            top_series=[
                {"values": merged[f"{parameter}_nominal"], "errors": merged[f"{parameter}_stat_nominal"], "fmt": "o", "label": "Nominal"},
                {"values": merged[f"{parameter}_isr"], "errors": merged[f"{parameter}_stat_isr"], "fmt": "s", "label": "ISR"},
            ],
            systematic=merged[f"{parameter}_absolute_isr_systematic"],
            status=merged[f"{parameter}_isr_barlow_status"],
            systematic_definition=definition,
            show_legends=show_legends, show_xlabel=show_xlabel,
        )
    # enddef

    plot_paths: list[str] = []
    for parameter in PHYSICS_PARAMETERS:
        fig, axes = plt.subplots(3, 1, figsize=SYSTEMATIC_COMPARISON_FIGSIZE, sharex=True)
        draw(parameter, axes, show_legends=True, show_xlabel=True)
        fig.tight_layout()
        path = plots_dir / f"nominal_vs_isr_{parameter}.png"
        fig.savefig(path, dpi=SYSTEMATIC_COMPARISON_DPI)
        plt.close(fig)
        plot_paths.append(str(path))
    # endfor
    summary = write_systematic_summary_canvas(
        output_dir=summary_dir, filename_stem="isr_systematic_summary",
        draw_parameter=draw,
    )

    csv_path = tables_dir / "nominal_vs_isr_structure_function_ratios.csv"
    json_path = tables_dir / "nominal_vs_isr_structure_function_ratios.json"
    merged.to_csv(csv_path, index=False)
    write_json(json_path, {
        "schema_version": 4,
        "systematic_definition": "delta_ISR = abs(ISR - nominal)",
        "middle_panel": "Assigned systematic; filled marker passes Barlow, open marker fails Barlow.",
        "bottom_panel": "Assigned systematic divided by nominal statistical uncertainty.",
        "plot_axis_convention": "Top-panel axes are common by parameter family; every assigned-systematic panel uses 0 to 0.2; normalized-size axis is logarithmic from 1e-2 to 1e1.",
        "barlow_summary": barlow_records,
        "rows": merged.to_dict(orient="records"),
    })
    return {"csv": str(csv_path), "json": str(json_path), "plots_directory": str(plots_dir), "plots": plot_paths, "summary": summary, "covariance_directory": str(covariance_dir), "covariance": covariance_products, "barlow_summary": barlow_records}

def write_momentum_correction_comparison_products(
    corrected: pd.DataFrame,
    uncorrected: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Write corrected-versus-uncorrected products using the common layout."""
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    covariance_dir = output_dir / "covariance"
    summary_dir = output_dir / "summary"
    for directory in (tables_dir, plots_dir, covariance_dir, summary_dir):
        ensure_directory(directory)
    # endfor
    keys = ["bin_number", "x_index", "t_index"]
    keep = keys + ["min_xB", "median_xB", "mean_xB", "max_xB", "min_Q2_gev2", "median_Q2_gev2", "mean_Q2_gev2", "max_Q2_gev2", "min_W_gev", "median_W_gev", "mean_W_gev", "max_W_gev", "min_minus_t_gev2", "median_minus_t_gev2", "mean_minus_t_gev2", "max_minus_t_gev2", "min_minus_tprime_gev2", "median_minus_tprime_gev2", "mean_minus_tprime_gev2", "max_minus_tprime_gev2", "min_epsilon", "median_epsilon", "mean_epsilon", "max_epsilon", "min_DepA", "median_DepA", "mean_DepA", "max_DepA", "min_DepB", "median_DepB", "mean_DepB", "max_DepB", "min_DepC", "median_DepC", "mean_DepC", "max_DepC", "min_DepV", "median_DepV", "mean_DepV", "max_DepV", "min_DepW", "median_DepW", "mean_DepW", "max_DepW", "number_of_events"] + [item for parameter in PHYSICS_PARAMETERS for item in (parameter, f"{parameter}_stat")]
    merged = corrected[keep].merge(uncorrected[keep], on=keys, suffixes=("_corrected", "_uncorrected"), validate="one_to_one")
    barlow_records: list[dict[str, Any]] = []
    covariance_products: dict[str, str] = {}
    for parameter in PHYSICS_PARAMETERS:
        shift = merged[f"{parameter}_uncorrected"] - merged[f"{parameter}_corrected"]
        merged[f"{parameter}_uncorrected_minus_corrected"] = shift
        systematic_column = f"{parameter}_absolute_momentum_correction_systematic"
        merged[systematic_column] = shift.abs()
        merged[f"{parameter}_momentum_systematic_to_nominal_stat"] = np.divide(merged[systematic_column], merged[f"{parameter}_stat_corrected"])
        barlow = calculate_barlow_arrays(merged[f"{parameter}_corrected"], merged[f"{parameter}_uncorrected"], merged[f"{parameter}_stat_corrected"], merged[f"{parameter}_stat_uncorrected"])
        prefix = f"{parameter}_momentum_correction_barlow"
        merged[f"{prefix}_statistical_scale"] = barlow["denominator"]
        merged[f"{prefix}_value"] = barlow["barlow"]
        merged[f"{prefix}_passes"] = barlow["passed"]
        merged[f"{prefix}_status"] = barlow["status"]
        barlow_records.append(summarize_barlow_status(merged, effect="momentum_corrections", variation="uncorrected_minus_corrected", parameter=parameter, barlow_column=f"{prefix}_value", status_column=f"{prefix}_status"))
        covariance_path = covariance_dir / f"{parameter}_momentum_correction_covariance.npy"
        vector = shift.to_numpy(dtype=np.float64)
        np.save(covariance_path, np.outer(vector, vector))
        covariance_products[parameter] = str(covariance_path)
    # endfor
    x = merged["bin_number"].to_numpy(dtype=float)
    definition = r"$\delta_{\rm mom}=|a_{\rm uncorrected}-a_{\rm corrected}|$"
    def draw(parameter: str, axes: np.ndarray, *, show_legends: bool, show_xlabel: bool) -> None:
        draw_systematic_comparison(axes, x=x, parameter=parameter, top_series=[
            {"values": merged[f"{parameter}_corrected"], "errors": merged[f"{parameter}_stat_corrected"], "fmt": "o", "label": "Momentum corrected (nominal)"},
            {"values": merged[f"{parameter}_uncorrected"], "errors": merged[f"{parameter}_stat_uncorrected"], "fmt": "s", "label": "No momentum corrections"},
        ], systematic=merged[f"{parameter}_absolute_momentum_correction_systematic"], status=merged[f"{parameter}_momentum_correction_barlow_status"], systematic_definition=definition, show_legends=show_legends, show_xlabel=show_xlabel)
    # enddef
    plot_paths: list[str] = []
    for parameter in PHYSICS_PARAMETERS:
        fig, axes = plt.subplots(3, 1, figsize=SYSTEMATIC_COMPARISON_FIGSIZE, sharex=True)
        draw(parameter, axes, show_legends=True, show_xlabel=True)
        fig.tight_layout()
        path = plots_dir / f"momentum_corrections_{parameter}.png"
        fig.savefig(path, dpi=SYSTEMATIC_COMPARISON_DPI)
        plt.close(fig)
        plot_paths.append(str(path))
    # endfor
    summary = write_systematic_summary_canvas(output_dir=summary_dir, filename_stem="momentum_corrections_systematic_summary", draw_parameter=draw)
    csv_path = tables_dir / "momentum_correction_structure_function_ratios.csv"
    json_path = tables_dir / "momentum_correction_structure_function_ratios.json"
    merged.to_csv(csv_path, index=False)
    write_json(json_path, {"schema_version": 3, "systematic_definition": "delta_mom = abs(uncorrected - corrected)", "middle_panel": "Assigned systematic; filled marker passes Barlow, open marker fails Barlow.", "bottom_panel": "Assigned systematic divided by corrected nominal statistical uncertainty.", "plot_axis_convention": "Top-panel axes are common by parameter family; every assigned-systematic panel uses 0 to 0.2; normalized-size axis is logarithmic from 1e-2 to 1e1.", "barlow_summary": barlow_records, "rows": merged.to_dict(orient="records")})
    return {"csv": str(csv_path), "json": str(json_path), "plots_directory": str(plots_dir), "plots": plot_paths, "summary": summary, "covariance_directory": str(covariance_dir), "covariance": covariance_products, "barlow_summary": barlow_records}

def write_fit_method_diagnostic_products(
    nominal: pd.DataFrame,
    centroid_only: pd.DataFrame,
    width_only: pd.DataFrame,
    polynomial: pd.DataFrame,
    output_dir: Path,
    carbon_cut_json: Path | None = None,
    polynomial_cut_json: Path | None = None,
) -> dict[str, Any]:
    """Write the four-way fit-method decomposition diagnostic."""
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    for directory in (output_dir, tables_dir, plots_dir):
        ensure_directory(directory)
    # endfor

    keys = ["bin_number", "x_index", "t_index"]
    keep = keys + ["number_of_events"] + [
        item for parameter in PHYSICS_PARAMETERS
        for item in (parameter, f"{parameter}_stat")
    ]
    frames = {
        "carbon": nominal[keep],
        "mu_poly_sigma_carbon": centroid_only[keep],
        "mu_carbon_sigma_poly": width_only[keep],
        "polynomial": polynomial[keep],
    }
    merged = None
    for label, frame in frames.items():
        renamed = frame.rename(columns={
            column: f"{column}_{label}" for column in frame.columns if column not in keys
        })
        merged = renamed if merged is None else merged.merge(renamed, on=keys, validate="one_to_one")
    # endfor

    for parameter in PHYSICS_PARAMETERS:
        base = merged[f"{parameter}_carbon"]
        stat = merged[f"{parameter}_stat_carbon"]
        for label in ("mu_poly_sigma_carbon", "mu_carbon_sigma_poly", "polynomial"):
            shift = merged[f"{parameter}_{label}"] - base
            merged[f"{parameter}_{label}_minus_carbon"] = shift
            merged[f"{parameter}_{label}_absolute_shift"] = shift.abs()
            merged[f"{parameter}_{label}_shift_over_carbon_stat"] = np.divide(shift.abs(), stat)
        # endfor
    # endfor

    csv_path = tables_dir / "fit_method_mu_sigma_decomposition.csv"
    merged.to_csv(csv_path, index=False)
    json_path = tables_dir / "fit_method_mu_sigma_decomposition.json"
    write_json(json_path, {
        "schema_version": 2,
        "diagnostic_only": True,
        "systematic_assignment": "none",
        "selections": {
            "carbon": "A(mu_C,sigma_C)",
            "mu_poly_sigma_carbon": "A(mu_poly,sigma_C): centroid-only variation",
            "mu_carbon_sigma_poly": "A(mu_C,sigma_poly): width-only variation",
            "polynomial": "A(mu_poly,sigma_poly): full polynomial-only variation",
        },
        "dilution_factors": "identical nominal dilution factors in all four extractions",
        "rows": merged.to_dict(orient="records"),
    })

    bins = merged["bin_number"].to_numpy(dtype=float)
    plot_paths: list[str] = []
    labels = {
        "carbon": r"$A(\mu_C,\sigma_C)$",
        "mu_poly_sigma_carbon": r"$A(\mu_{\rm poly},\sigma_C)$",
        "mu_carbon_sigma_poly": r"$A(\mu_C,\sigma_{\rm poly})$",
        "polynomial": r"$A(\mu_{\rm poly},\sigma_{\rm poly})$",
    }
    markers = {"carbon": "o", "mu_poly_sigma_carbon": "^", "mu_carbon_sigma_poly": "v", "polynomial": "s"}
    offsets = {"carbon": -0.24, "mu_poly_sigma_carbon": -0.08, "mu_carbon_sigma_poly": 0.08, "polynomial": 0.24}
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        fig, ax = plt.subplots(figsize=(14, 5.5))
        for label in ("carbon", "mu_poly_sigma_carbon", "mu_carbon_sigma_poly", "polynomial"):
            ax.errorbar(
                bins + offsets[label], merged[f"{parameter}_{label}"],
                yerr=merged[f"{parameter}_stat_{label}"], marker=markers[label],
                linestyle="none", capsize=2, label=labels[label],
            )
        # endfor
        ax.axhline(0.0, linewidth=0.8)
        ax.set_xlabel("Combined kinematic-bin number")
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_xticks(bins)
        apply_parameter_y_limits(ax, parameter)
        ax.grid(alpha=0.25)
        ax.legend(ncol=2)
        ax.set_title("Fit-method decomposition (statistical uncertainties only)")
        fig.tight_layout()
        path = plots_dir / f"fit_method_decomposition_{parameter}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        plot_paths.append(str(path))
    # endfor

    # Relate the displacement of the fitted neutron-peak centroid to the
    # quality of the polynomial-only missing-mass fit.  This is diagnostic:
    # it is not used to assign or rescale any systematic uncertainty.
    centroid_quality_csv = None
    centroid_quality_plot = None
    if carbon_cut_json is not None and polynomial_cut_json is not None:
        carbon_cuts = read_json_tolerating_trailing_escaped_whitespace(carbon_cut_json)
        polynomial_cuts = read_json_tolerating_trailing_escaped_whitespace(polynomial_cut_json)

        def _mu_by_bin(payload: dict[str, Any]) -> dict[int, float]:
            rows = payload.get("flat_rows")
            if rows:
                return {int(row["bin_number"]): float(row.get("shared_mean_gev2", row.get("mu_gev2"))) for row in rows}
            # endif
            first_period = next(iter(payload["periods"].values()))
            return {int(row["bin_number"]): float(row["mu_gev2"]) for row in first_period}

        carbon_mu = _mu_by_bin(carbon_cuts)
        polynomial_mu = _mu_by_bin(polynomial_cuts)

        # channel_selection_mx2_fits.py writes the full polynomial fit table
        # two directory levels above final_polynomial_only_cuts/tables/.
        fit_table_path = polynomial_cut_json.parents[2] / "tables" / "mx2_peak_fit_results_v27.csv"
        if fit_table_path.is_file():
            fit_table = pd.read_csv(fit_table_path)
            selected = fit_table[(fit_table["stage"] == "after") & fit_table["is_recommended"].astype(bool)].copy()
            quality = (
                selected.groupby("bin_number", as_index=False)
                .agg(polynomial_chi2_ndf=("joint_chi2_ndf", "first"))
                .sort_values("bin_number")
            )
            quality["mu_carbon_gev2"] = quality["bin_number"].map(carbon_mu)
            quality["mu_polynomial_gev2"] = quality["bin_number"].map(polynomial_mu)
            quality["delta_mu_gev2"] = quality["mu_polynomial_gev2"] - quality["mu_carbon_gev2"]
            quality["absolute_delta_mu_gev2"] = quality["delta_mu_gev2"].abs()

            centroid_quality_csv = tables_dir / "polynomial_centroid_shift_vs_fit_quality.csv"
            quality.to_csv(centroid_quality_csv, index=False)

            fig, ax = plt.subplots(figsize=(9.0, 6.5))
            ax.scatter(quality["polynomial_chi2_ndf"], quality["absolute_delta_mu_gev2"], s=42, zorder=3)
            for row in quality.itertuples(index=False):
                ax.annotate(
                    f"bin {int(row.bin_number)}",
                    (row.polynomial_chi2_ndf, row.absolute_delta_mu_gev2),
                    xytext=(5, 4), textcoords="offset points", fontsize=8,
                )
            # endfor
            ax.set_xlabel(r"Polynomial-only fit $\chi^2/\mathrm{ndf}$")
            ax.set_ylabel(r"$|\mu_{\rm poly}-\mu_C|$ (GeV$^2$)")
            ax.set_title("Polynomial-only centroid displacement vs fit quality")
            ax.grid(alpha=0.25)
            fig.tight_layout()
            centroid_quality_plot = plots_dir / "polynomial_centroid_shift_vs_fit_quality.png"
            fig.savefig(centroid_quality_plot, dpi=180)
            plt.close(fig)
            plot_paths.append(str(centroid_quality_plot))
        else:
            print(
                "[fit-method-diagnostic] WARNING: polynomial fit-quality table not found at "
                f"{fit_table_path}; skipping centroid-shift-vs-chi2 plot.",
                flush=True,
            )
        # endif
    # endif

    return {
        "csv": str(csv_path), "json": str(json_path), "plots": plot_paths,
        "centroid_quality_csv": str(centroid_quality_csv) if centroid_quality_csv else None,
        "centroid_quality_plot": str(centroid_quality_plot) if centroid_quality_plot else None,
    }


def write_channel_selection_comparison_products(
    *,
    tight: pd.DataFrame,
    nominal: pd.DataFrame,
    loose: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Write channel-selection products using the common layout."""
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    covariance_dir = output_dir / "covariance"
    summary_dir = output_dir / "summary"
    for directory in (tables_dir, plots_dir, covariance_dir, summary_dir):
        ensure_directory(directory)
    # endfor
    keys = ["bin_number", "x_index", "t_index"]
    keep = keys + ["min_xB", "median_xB", "mean_xB", "max_xB", "min_Q2_gev2", "median_Q2_gev2", "mean_Q2_gev2", "max_Q2_gev2", "min_W_gev", "median_W_gev", "mean_W_gev", "max_W_gev", "min_minus_t_gev2", "median_minus_t_gev2", "mean_minus_t_gev2", "max_minus_t_gev2", "min_minus_tprime_gev2", "median_minus_tprime_gev2", "mean_minus_tprime_gev2", "max_minus_tprime_gev2", "min_epsilon", "median_epsilon", "mean_epsilon", "max_epsilon", "min_DepA", "median_DepA", "mean_DepA", "max_DepA", "min_DepB", "median_DepB", "mean_DepB", "max_DepB", "min_DepC", "median_DepC", "mean_DepC", "max_DepC", "min_DepV", "median_DepV", "mean_DepV", "max_DepV", "min_DepW", "median_DepW", "mean_DepW", "max_DepW", "number_of_events"] + [item for parameter in PHYSICS_PARAMETERS for item in (parameter, f"{parameter}_stat")]
    merged = nominal[keep].merge(tight[keep], on=keys, suffixes=("_nominal", "_tight"), validate="one_to_one").merge(loose[keep], on=keys, validate="one_to_one")
    merged = merged.rename(columns={column: f"{column}_loose" for column in keep if column not in keys})
    covariance_products: dict[str, str] = {}
    barlow_records: list[dict[str, Any]] = []
    for parameter in PHYSICS_PARAMETERS:
        delta_tight = merged[f"{parameter}_tight"] - merged[f"{parameter}_nominal"]
        delta_loose = merged[f"{parameter}_loose"] - merged[f"{parameter}_nominal"]
        merged[f"{parameter}_tight_minus_nominal"] = delta_tight
        merged[f"{parameter}_loose_minus_nominal"] = delta_loose
        syst = np.sqrt(0.5 * (delta_tight**2 + delta_loose**2))
        merged[f"{parameter}_channel_selection_rms_systematic"] = syst
        merged[f"{parameter}_channel_selection_systematic_to_nominal_stat"] = np.divide(syst, merged[f"{parameter}_stat_nominal"])
        for variation in ("tight", "loose"):
            barlow = calculate_barlow_arrays(merged[f"{parameter}_nominal"], merged[f"{parameter}_{variation}"], merged[f"{parameter}_stat_nominal"], merged[f"{parameter}_stat_{variation}"])
            prefix = f"{parameter}_{variation}_barlow"
            merged[f"{prefix}_statistical_scale"] = barlow["denominator"]
            merged[f"{prefix}_value"] = barlow["barlow"]
            merged[f"{prefix}_passes"] = barlow["passed"]
            merged[f"{prefix}_status"] = barlow["status"]
            barlow_records.append(summarize_barlow_status(merged, effect="channel_selection", variation=f"{variation}_minus_nominal", parameter=parameter, barlow_column=f"{prefix}_value", status_column=f"{prefix}_status"))
        # endfor
        tight_status = merged[f"{parameter}_tight_barlow_status"].astype(str)
        loose_status = merged[f"{parameter}_loose_barlow_status"].astype(str)
        merged[f"{parameter}_channel_selection_barlow_status"] = np.where((tight_status == "pass") | (loose_status == "pass"), "pass", np.where((tight_status == "fail") & (loose_status == "fail"), "fail", "undefined"))
        dt = delta_tight.to_numpy(dtype=np.float64)
        dl = delta_loose.to_numpy(dtype=np.float64)
        covariance_path = covariance_dir / f"{parameter}_channel_selection_covariance.npy"
        np.save(covariance_path, 0.5 * (np.outer(dt, dt) + np.outer(dl, dl)))
        covariance_products[parameter] = str(covariance_path)
    # endfor
    x = merged["bin_number"].to_numpy(dtype=float)
    definition = r"$\delta_{\rm ch}=\sqrt{(\Delta_{\rm tight}^2+\Delta_{\rm loose}^2)/2}$"
    def draw(parameter: str, axes: np.ndarray, *, show_legends: bool, show_xlabel: bool) -> None:
        draw_systematic_comparison(axes, x=x, parameter=parameter, top_series=[
            {"values": merged[f"{parameter}_nominal"], "errors": merged[f"{parameter}_stat_nominal"], "fmt": "o", "label": r"Nominal ($\mu\pm2\sigma$)"},
            {"values": merged[f"{parameter}_tight"], "errors": merged[f"{parameter}_stat_tight"], "fmt": "^", "label": r"Tight ($\mu\pm1\sigma$)"},
            {"values": merged[f"{parameter}_loose"], "errors": merged[f"{parameter}_stat_loose"], "fmt": "s", "label": r"Loose ($\mu\pm3\sigma$)"},
        ], systematic=merged[f"{parameter}_channel_selection_rms_systematic"], status=merged[f"{parameter}_channel_selection_barlow_status"], systematic_definition=definition, show_legends=show_legends, show_xlabel=show_xlabel)
    # enddef
    plot_paths: list[str] = []
    for parameter in PHYSICS_PARAMETERS:
        fig, axes = plt.subplots(3, 1, figsize=SYSTEMATIC_COMPARISON_FIGSIZE, sharex=True)
        draw(parameter, axes, show_legends=True, show_xlabel=True)
        fig.tight_layout()
        path = plots_dir / f"channel_selection_{parameter}.png"
        fig.savefig(path, dpi=SYSTEMATIC_COMPARISON_DPI)
        plt.close(fig)
        plot_paths.append(str(path))
    # endfor
    summary = write_systematic_summary_canvas(output_dir=summary_dir, filename_stem="channel_selection_systematic_summary", draw_parameter=draw)
    csv_path = tables_dir / "channel_selection_structure_function_ratios.csv"
    json_path = tables_dir / "channel_selection_structure_function_ratios.json"
    merged.to_csv(csv_path, index=False)
    write_json(json_path, {"schema_version": 4, "systematic_definition": "delta_ch = sqrt(((tight-nominal)^2 + (loose-nominal)^2)/2)", "combined_barlow_marker_rule": "Filled if either tight or loose variation passes Barlow; open only if both defined variations fail; x if unresolved/undefined.", "middle_panel": "Assigned RMS systematic with Barlow marker state.", "bottom_panel": "Assigned RMS systematic divided by nominal statistical uncertainty.", "plot_axis_convention": "Top-panel axes are common by parameter family; every assigned-systematic panel uses 0 to 0.2; normalized-size axis is logarithmic from 1e-2 to 1e1.", "barlow_summary": barlow_records, "rows": merged.to_dict(orient="records")})
    return {"csv": str(csv_path), "json": str(json_path), "plots_directory": str(plots_dir), "plots": plot_paths, "summary": summary, "covariance_directory": str(covariance_dir), "covariance": covariance_products, "barlow_summary": barlow_records}

def calculate_migration_systematics(
    nominal: pd.DataFrame,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
    """Calculate current-result migration shifts from the established matrix.

    The published matrix rows are rounded generated-bin compositions of each
    reconstructed bin.  Renormalizing each row is essential: a constant input
    observable must remain constant after migration despite the 0.98--1.01 row
    sums introduced by two-decimal-place publication rounding.

    Returns
    -------
    migrated
        Migration-smeared values for the five published polarized ratios.
    assigned
        1.5 * abs(migrated - nominal), the established conservative systematic.
    matrix_normalized
        Row-normalized matrix actually used in the calculation.
    """
    if len(nominal) != NUMBER_OF_BINS:
        raise ValueError(
            f"Migration calculation requires {NUMBER_OF_BINS} nominal bins; "
            f"received {len(nominal)}."
        )
    # endif
    ordered = nominal.sort_values("bin_number")
    expected_bins = np.arange(1, NUMBER_OF_BINS + 1, dtype=int)
    actual_bins = ordered["bin_number"].to_numpy(dtype=int)
    if not np.array_equal(actual_bins, expected_bins):
        raise ValueError(
            "Migration matrix assumes nominal bins are numbered consecutively "
            f"1..{NUMBER_OF_BINS}; received {actual_bins.tolist()}."
        )
    # endif

    row_sums = MIGRATION_MATRIX.sum(axis=1)
    if np.any(~np.isfinite(row_sums)) or np.any(row_sums <= 0.0):
        raise ValueError("Migration matrix contains a non-finite or empty row.")
    # endif
    matrix_normalized = MIGRATION_MATRIX / row_sums[:, None]

    # Sanity check: row normalization must preserve a constant observable.
    constant_check = matrix_normalized @ np.ones(NUMBER_OF_BINS, dtype=float)
    if not np.allclose(constant_check, 1.0, rtol=0.0, atol=1.0e-12):
        raise RuntimeError("Normalized migration matrix does not preserve constants.")
    # endif

    migrated: dict[str, np.ndarray] = {}
    assigned: dict[str, np.ndarray] = {}
    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        values = pd.to_numeric(ordered[parameter], errors="coerce").to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(
                f"Cannot calculate migration systematic for {parameter}: "
                "nominal values contain non-finite entries."
            )
        # endif
        smeared = matrix_normalized @ values
        migrated[parameter] = smeared
        assigned[parameter] = MIGRATION_RESOLUTION_SCALE * np.abs(smeared - values)
    # endfor
    return migrated, assigned, matrix_normalized


def write_migration_systematic_products(
    nominal: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Write migration matrix, per-bin shifts, covariance, and diagnostics."""
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    covariance_dir = output_dir / "covariance"
    for directory in (tables_dir, plots_dir, covariance_dir):
        ensure_directory(directory)
    # endfor

    ordered = nominal.sort_values("bin_number").reset_index(drop=True)
    migrated, assigned, matrix_normalized = calculate_migration_systematics(ordered)
    row_sums = MIGRATION_MATRIX.sum(axis=1)

    matrix_raw_path = tables_dir / "migration_matrix_published_rounded.csv"
    matrix_norm_path = tables_dir / "migration_matrix_row_normalized.csv"
    pd.DataFrame(MIGRATION_MATRIX).to_csv(matrix_raw_path, index=False)
    pd.DataFrame(matrix_normalized).to_csv(matrix_norm_path, index=False)

    columns: dict[str, Any] = {
        "bin_number": ordered["bin_number"].to_numpy(dtype=int),
        "x_index": ordered["x_index"].to_numpy(dtype=int),
        "t_index": ordered["t_index"].to_numpy(dtype=int),
        "migration_matrix_published_row_sum": row_sums,
    }
    covariance_products: dict[str, str] = {}
    plot_paths: list[str] = []
    x = ordered["bin_number"].to_numpy(dtype=float)

    for parameter in PUBLISHED_SYSTEMATIC_PARAMETERS:
        nominal_values = ordered[parameter].to_numpy(dtype=float)
        smeared = migrated[parameter]
        signed_base_shift = smeared - nominal_values
        signed_assigned_shift = MIGRATION_RESOLUTION_SCALE * signed_base_shift
        systematic = assigned[parameter]
        stat = ordered[f"{parameter}_stat"].to_numpy(dtype=float)

        columns[f"{parameter}_nominal"] = nominal_values
        columns[f"{parameter}_migrated"] = smeared
        columns[f"{parameter}_migration_signed_shift_unscaled"] = signed_base_shift
        columns[f"{parameter}_migration_signed_shift_scaled"] = signed_assigned_shift
        columns[f"{parameter}_migration_systematic"] = systematic
        columns[f"{parameter}_migration_systematic_to_stat"] = np.divide(
            systematic,
            stat,
            out=np.full_like(systematic, np.nan),
            where=stat > 0.0,
        )

        covariance_path = covariance_dir / f"{parameter}_migration_covariance.npy"
        np.save(covariance_path, np.outer(signed_assigned_shift, signed_assigned_shift))
        covariance_products[parameter] = str(covariance_path)

        fig, axes = plt.subplots(3, 1, figsize=SYSTEMATIC_COMPARISON_FIGSIZE, sharex=True)
        axes[0].errorbar(
            x, nominal_values, yerr=stat, marker="o", linestyle="none",
            capsize=2, label="Nominal extraction",
        )
        axes[0].plot(x, smeared, "s", label="Migration-smeared nominal values")
        axes[0].set_ylabel(PARAMETER_LABELS[parameter])
        apply_parameter_y_limits(axes[0], parameter)
        axes[0].legend()
        axes[0].grid(alpha=0.25)

        axes[1].bar(x, systematic, width=0.70, alpha=0.55)
        axes[1].set_ylabel("Assigned systematic")
        axes[1].text(
            0.02, 0.95,
            r"$\delta_{\rm migr}=1.5\,|M_{\rm row\,norm}a-a|$",
            transform=axes[1].transAxes, va="top",
        )
        axes[1].set_ylim(*systematic_size_y_limits(parameter))
        axes[1].grid(alpha=0.25)

        ratio = np.divide(
            systematic, stat, out=np.full_like(systematic, np.nan), where=stat > 0.0,
        )
        axes[2].plot(x, ratio, "o")
        axes[2].axhline(1.0, linewidth=0.8)
        axes[2].set_yscale("log")
        axes[2].set_ylim(*SYSTEMATIC_TO_STAT_RATIO_Y_LIMITS)
        axes[2].set_ylabel(r"$\delta_{\rm migr}/\sigma_{\rm stat}$")
        axes[2].set_xlabel("Combined kinematic-bin number")
        axes[2].grid(alpha=0.25)
        fig.tight_layout()
        path = plots_dir / f"migration_{parameter}.png"
        fig.savefig(path, dpi=SYSTEMATIC_COMPARISON_DPI)
        plt.close(fig)
        plot_paths.append(str(path))
    # endfor

    frame = pd.DataFrame(columns)
    csv_path = tables_dir / "migration_systematics.csv"
    json_path = tables_dir / "migration_systematics.json"
    frame.to_csv(csv_path, index=False)
    write_json(json_path, {
        "schema_version": 1,
        "published_parameters_only": list(PUBLISHED_SYSTEMATIC_PARAMETERS),
        "matrix_definition": (
            "Published row i gives generated-bin composition j of reconstructed bin i. "
            "Rows are renormalized before application because the published entries are "
            "rounded to two decimal places."
        ),
        "resolution_scale_factor": MIGRATION_RESOLUTION_SCALE,
        "systematic_definition": (
            "delta_migration[i] = 1.5 * abs(sum_j M_row_normalized[i,j] * "
            "a_nominal[j] - a_nominal[i])"
        ),
        "published_matrix_row_sums_before_normalization": row_sums.tolist(),
        "raw_matrix_csv": str(matrix_raw_path),
        "normalized_matrix_csv": str(matrix_norm_path),
        "covariance": covariance_products,
        "rows": frame.to_dict(orient="records"),
    })
    return {
        "csv": str(csv_path),
        "json": str(json_path),
        "raw_matrix_csv": str(matrix_raw_path),
        "normalized_matrix_csv": str(matrix_norm_path),
        "plots": plot_paths,
        "covariance": covariance_products,
    }


def attach_point_to_point_systematics(
    nominal: pd.DataFrame,
    isr: pd.DataFrame | None,
    momentum_uncorrected: pd.DataFrame | None,
    channel_tight: pd.DataFrame | None,
    channel_loose: pd.DataFrame | None,
) -> pd.DataFrame:
    """Attach source-by-source and final point-to-point uncertainties.

    For the five published polarized observables, the assigned total is

        delta_ptp = sqrt(delta_rad^2 + delta_axis^2
                         + delta_ch^2 + delta_migration^2).

    The target-axis study envelope and momentum-correction on/off difference are retained in diagnostic
    columns but are deliberately excluded from the assigned total.  The two UU
    modulations are fitted nuisance quantities and receive no migration term;
    they are excluded from publication-level systematic summaries.
    """
    output = nominal.copy().sort_values("bin_number").reset_index(drop=True)
    _, migration_assigned, _ = calculate_migration_systematics(output)

    for parameter in PHYSICS_PARAMETERS:
        nominal_values = pd.to_numeric(output[parameter], errors="coerce").to_numpy(dtype=float)

        if isr is None:
            radiation = np.zeros_like(nominal_values)
        else:
            isr_ordered = isr.sort_values("bin_number")
            radiation = np.abs(
                pd.to_numeric(isr_ordered[parameter], errors="coerce").to_numpy(dtype=float)
                - nominal_values
            )
        # endif

        target_axis = pd.to_numeric(
            output[f"{parameter}_target_axis_study_envelope"], errors="coerce"
        ).to_numpy(dtype=float)

        if momentum_uncorrected is None:
            momentum_diagnostic = np.zeros_like(nominal_values)
        else:
            momentum_ordered = momentum_uncorrected.sort_values("bin_number")
            momentum_diagnostic = np.abs(
                pd.to_numeric(momentum_ordered[parameter], errors="coerce").to_numpy(dtype=float)
                - nominal_values
            )
        # endif

        if channel_tight is None or channel_loose is None:
            channel = np.zeros_like(nominal_values)
        else:
            tight_ordered = channel_tight.sort_values("bin_number")
            loose_ordered = channel_loose.sort_values("bin_number")
            tight_shift = (
                pd.to_numeric(tight_ordered[parameter], errors="coerce").to_numpy(dtype=float)
                - nominal_values
            )
            loose_shift = (
                pd.to_numeric(loose_ordered[parameter], errors="coerce").to_numpy(dtype=float)
                - nominal_values
            )
            channel = np.sqrt(0.5 * (tight_shift**2 + loose_shift**2))
        # endif

        migration = (
            migration_assigned[parameter]
            if parameter in PUBLISHED_SYSTEMATIC_PARAMETERS
            else np.zeros_like(nominal_values)
        )
        total = np.sqrt(radiation**2 + channel**2 + migration**2)

        output[f"{parameter}_radiation_systematic"] = radiation
        output[f"{parameter}_target_axis_study_envelope"] = target_axis
        output[f"{parameter}_momentum_correction_difference_diagnostic"] = momentum_diagnostic
        # Keep the legacy column name for downstream compatibility, but make
        # its diagnostic-only status explicit in the new manifest and JSON.
        output[f"{parameter}_momentum_correction_systematic"] = momentum_diagnostic
        output[f"{parameter}_channel_selection_systematic"] = channel
        output[f"{parameter}_bin_migration_systematic"] = migration
        output[f"{parameter}_point_to_point_systematic"] = total
        stat = pd.to_numeric(output[f"{parameter}_stat"], errors="coerce").to_numpy(dtype=float)
        output[f"{parameter}_point_to_point_systematic_to_stat"] = np.divide(
            total,
            stat,
            out=np.full_like(total, np.nan),
            where=stat > 0.0,
        )
    # endfor
    return output


def write_target_axis_study_products(
    nominal_frame: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Write target-axis tables, plots, covariance, and Barlow diagnostics.

    The output table is assembled from a dictionary of complete columns and
    materialized once.  This avoids pandas DataFrame fragmentation from
    repeatedly inserting columns inside the parameter/variation loops.
    """
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    covariance_dir = output_dir / "covariance"
    summary_dir = output_dir / "summary"
    for directory in (tables_dir, plots_dir, covariance_dir, summary_dir):
        ensure_directory(directory)
    # endfor

    keys = [
        "bin_number",
        "x_index",
        "t_index",
        "min_xB",
        "mean_xB",
        "max_xB",
        "min_Q2_gev2",
        "mean_Q2_gev2",
        "max_Q2_gev2",
        "min_W_gev",
        "mean_W_gev",
        "max_W_gev",
        "min_minus_t_gev2",
        "mean_minus_t_gev2",
        "max_minus_t_gev2",
        "min_minus_tprime_gev2",
        "mean_minus_tprime_gev2",
        "max_minus_tprime_gev2",
    ]

    column_data: dict[str, Any] = {
        key: nominal_frame[key].to_numpy(copy=True)
        for key in keys
    }
    records: list[dict[str, Any]] = []
    covariance_products: dict[str, dict[str, str]] = {}

    for parameter in PHYSICS_PARAMETERS:
        nominal_values = nominal_frame[parameter].to_numpy(dtype=np.float64, copy=True)
        nominal_stat = nominal_frame[f"{parameter}_stat"].to_numpy(
            dtype=np.float64, copy=True
        )

        column_data[f"{parameter}_nominal"] = nominal_values
        column_data[f"{parameter}_stat_nominal"] = nominal_stat
        covariance_products[parameter] = {}

        variation_results: dict[str, dict[str, np.ndarray]] = {}
        for variation in ("photon_axis_projection", "external_data_informed"):
            variation_values = nominal_frame[f"{parameter}_{variation}"].to_numpy(
                dtype=np.float64, copy=True
            )
            variation_stat = nominal_frame[
                f"{parameter}_stat_{variation}"
            ].to_numpy(dtype=np.float64, copy=True)

            column_data[f"{parameter}_{variation}"] = variation_values
            column_data[f"{parameter}_stat_{variation}"] = variation_stat

            barlow = calculate_barlow_arrays(
                nominal_values,
                variation_values,
                nominal_stat,
                variation_stat,
            )
            prefix = f"{parameter}_{variation}_barlow"
            column_data[f"{prefix}_difference"] = np.asarray(
                barlow["difference"]
            )
            column_data[f"{prefix}_statistical_scale"] = np.asarray(
                barlow["denominator"]
            )
            column_data[f"{prefix}_value"] = np.asarray(barlow["barlow"])
            column_data[f"{prefix}_passes"] = np.asarray(barlow["passed"])
            column_data[f"{prefix}_status"] = np.asarray(
                barlow["status"], dtype=object
            )

            variation_results[variation] = {
                "difference": np.asarray(barlow["difference"]),
                "status": np.asarray(barlow["status"], dtype=object),
            }

            summary_frame = pd.DataFrame(
                {
                    f"{prefix}_value": column_data[f"{prefix}_value"],
                    f"{prefix}_status": column_data[f"{prefix}_status"],
                }
            )
            records.append(
                summarize_barlow_status(
                    summary_frame,
                    effect="target_axis",
                    variation=f"{variation}_minus_nominal",
                    parameter=parameter,
                    barlow_column=f"{prefix}_value",
                    status_column=f"{prefix}_status",
                )
            )

            vector = variation_results[variation]["difference"].astype(
                np.float64, copy=False
            )
            covariance_path = (
                covariance_dir / f"{parameter}_{variation}_covariance.npy"
            )
            np.save(covariance_path, np.outer(vector, vector))
            covariance_products[parameter][variation] = str(covariance_path)
        # endfor

        no_difference = variation_results["photon_axis_projection"]["difference"]
        ext_difference = variation_results["external_data_informed"][
            "difference"
        ]
        no_abs = np.abs(no_difference)
        ext_abs = np.abs(ext_difference)
        choose_no = no_abs >= ext_abs
        target_axis_study_envelope = np.maximum(no_abs, ext_abs)  # diagnostic envelope only

        systematic_to_stat = np.divide(
            target_axis_study_envelope,
            nominal_stat,
            out=np.full_like(target_axis_study_envelope, np.nan, dtype=np.float64),
            where=np.isfinite(nominal_stat) & (nominal_stat > 0.0),
        )

        no_status = variation_results["photon_axis_projection"]["status"]
        ext_status = variation_results["external_data_informed"]["status"]

        column_data[f"{parameter}_target_axis_study_envelope"] = (
            target_axis_study_envelope
        )
        column_data[
            f"{parameter}_target_axis_study_envelope_to_nominal_stat"
        ] = systematic_to_stat
        column_data[f"{parameter}_target_axis_selected_variation"] = np.where(
            choose_no, "photon_axis_projection", "external_data_informed"
        )
        column_data[f"{parameter}_target_axis_study_status"] = np.where(
            choose_no, no_status, ext_status
        )
    # endfor

    frame = pd.DataFrame(column_data)

    x = frame["bin_number"].to_numpy(dtype=float)
    definition = (
        r"$\delta_{\rm axis}=\max(|a_{\rm no\ proj}-a_{\rm nom}|,"
        r"|a_{\rm ext}-a_{\rm nom}|)$"
    )

    def draw(
        parameter: str,
        axes: np.ndarray,
        *,
        show_legends: bool,
        show_xlabel: bool,
    ) -> None:
        draw_systematic_comparison(
            axes,
            x=x,
            parameter=parameter,
            top_series=[
                {
                    "values": frame[f"{parameter}_nominal"],
                    "errors": frame[f"{parameter}_stat_nominal"],
                    "fmt": "o",
                    "label": "Nominal projection",
                },
                {
                    "values": frame[f"{parameter}_photon_axis_projection"],
                    "errors": frame[f"{parameter}_stat_photon_axis_projection"],
                    "fmt": "^",
                    "label": "Photon-axis projection",
                },
                {
                    "values": frame[f"{parameter}_external_data_informed"],
                    "errors": frame[
                        f"{parameter}_stat_external_data_informed"
                    ],
                    "fmt": "s",
                    "label": "External-data-informed",
                },
            ],
            systematic=frame[f"{parameter}_target_axis_study_envelope"],
            status=frame[f"{parameter}_target_axis_study_status"],
            systematic_definition=definition,
            show_legends=show_legends,
            show_xlabel=show_xlabel,
        )
        axes[1].set_ylabel("Diagnostic envelope")
        axes[2].set_ylabel("Envelope / nominal stat.")
    # enddef

    plot_paths: list[str] = []
    for parameter in PHYSICS_PARAMETERS:
        fig, axes = plt.subplots(
            3, 1, figsize=SYSTEMATIC_COMPARISON_FIGSIZE, sharex=True
        )
        draw(parameter, axes, show_legends=True, show_xlabel=True)
        fig.tight_layout()
        path = plots_dir / f"target_axis_{parameter}.png"
        fig.savefig(path, dpi=SYSTEMATIC_COMPARISON_DPI)
        plt.close(fig)
        plot_paths.append(str(path))
    # endfor

    summary = write_systematic_summary_canvas(
        output_dir=summary_dir,
        filename_stem="target_axis_study_envelope_summary",
        draw_parameter=draw,
    )
    csv_path = tables_dir / "target_axis_study_criteria.csv"
    json_path = tables_dir / "target_axis_study_criteria.json"
    frame.to_csv(csv_path, index=False)
    write_json(
        json_path,
        {
            "schema_version": 3,
            "study_envelope_definition": (
                "diagnostic envelope = max(abs(photon_axis_projection-nominal), "
                "abs(external_data_informed-nominal)); not assigned as a systematic"
            ),
            "selected_barlow_marker_rule": (
                "The marker state is taken from the variation that supplies "
                "the envelope in that bin."
            ),
            "middle_panel": (
                "Target-axis diagnostic envelope with Barlow marker state; not assigned as a systematic."
            ),
            "bottom_panel": (
                "Target-axis diagnostic envelope divided by nominal statistical "
                "uncertainty."
            ),
            "plot_axis_convention": (
                "Top and middle axes are common by parameter family; "
                "normalized-size axis is logarithmic from 1e-2 to 1e1."
            ),
            "summary": records,
            "rows": frame.to_dict(orient="records"),
        },
    )
    return {
        "csv": str(csv_path),
        "json": str(json_path),
        "plots_directory": str(plots_dir),
        "plots": plot_paths,
        "summary_canvas": summary,
        "covariance_directory": str(covariance_dir),
        "covariance": covariance_products,
        "barlow_summary": records,
    }


# =============================================================================
# RGA Diehl et al. exclusive-pi+ cross-check
# =============================================================================
# Supplemental material to S. Diehl et al., Phys. Lett. B 839 (2023) 137761.
# Q2-xB polygon vertices are connected in the order listed in supplemental
# Table 10; -t edges are supplemental Table 11.  Published points below are
# supplemental Tables 1--9.
RGA_Q2_XB_POLYGONS = {
    1: ((0.095,1.5),(0.21,3.2),(0.21,1.5)),
    2: ((0.21,1.5),(0.21,2.2),(0.30,2.7),(0.30,1.5)),
    3: ((0.21,2.2),(0.21,3.2),(0.30,4.56),(0.30,2.7)),
    4: ((0.30,1.5),(0.30,2.7),(0.37,3.2),(0.37,1.765),(0.332,1.5)),
    5: ((0.30,2.7),(0.30,4.57),(0.37,5.62),(0.37,3.2)),
    6: ((0.37,1.765),(0.37,3.2),(0.45,3.85),(0.45,2.48)),
    7: ((0.37,3.2),(0.37,5.62),(0.45,6.8),(0.45,3.85)),
    8: ((0.45,2.48),(0.45,3.85),(0.67,6.14),(0.63,5.15),(0.57,4.05),(0.50,3.05)),
    9: ((0.45,3.85),(0.45,6.8),(0.677,10.185),(0.7896,11.351),(0.75,9.52),(0.708,7.42),(0.67,6.14)),
}
RGA_MINUS_T_EDGES = {
    1:(0.01,0.07,0.12,0.21,0.35,0.60,0.90),
    2:(0.02,0.15,0.25,0.40,0.60,0.85),
    3:(0.02,0.15,0.25,0.40,0.66,0.90),
    4:(0.05,0.21,0.31,0.45,0.65,0.90),
    5:(0.05,0.23,0.34,0.55,0.90),
    6:(0.10,0.31,0.425,0.59,0.90),
    7:(0.10,0.34,0.46,0.63,0.90),
    8:(0.20,0.52,0.75,0.95,1.20),
    9:(0.20,0.62,0.85,1.00,1.20),
}
RGA_PUBLISHED = {
1:[(1.811,.175,.052,.710,.0292,.0117,.0095),(1.841,.183,.093,.730,.0766,.0109,.0072),(1.851,.182,.160,.721,.1111,.0137,.0092),(1.858,.180,.273,.709,.1073,.0129,.0102),(1.861,.178,.464,.700,.0695,.0123,.0059),(1.864,.177,.740,.694,.0560,.0136,.0070)],
2:[(1.856,.252,.111,.866,.0556,.0078,.0064),(1.884,.263,.195,.874,.1202,.0094,.0112),(1.881,.264,.317,.875,.1450,.0108,.0105),(1.875,.263,.492,.875,.1096,.0114,.0102),(1.873,.263,.716,.875,.1112,.0124,.0123)],
3:[(2.784,.248,.111,.661,.0670,.0126,.0048),(2.882,.260,.195,.668,.1220,.0157,.0112),(2.900,.260,.318,.663,.1070,.0166,.0099),(2.909,.259,.517,.657,.0960,.0158,.0077),(2.913,.259,.774,.655,.0580,.0209,.0119)],
4:[(2.028,.327,.171,.906,.0887,.0102,.0066),(2.090,.334,.257,.904,.1092,.0099,.0092),(2.091,.335,.375,.905,.1443,.0123,.0130),(2.084,.335,.542,.906,.1540,.0118,.0132),(2.097,.335,.765,.904,.1253,.0126,.0110)],
5:[(3.440,.327,.185,.700,.0561,.0155,.0093),(3.515,.335,.281,.701,.1217,.0155,.0119),(3.530,.336,.434,.699,.1886,.0182,.0138),(3.549,.336,.706,.695,.1015,.0178,.0064)],
6:[(2.518,.395,.260,.899,.1016,.0126,.0104),(2.629,.405,.365,.895,.1337,.0119,.0119),(2.646,.407,.502,.894,.1781,.0140,.0123),(2.644,.407,.729,.895,.1580,.0120,.0118)],
7:[(4.108,.398,.282,.707,.1210,.0174,.0133),(4.224,.408,.398,.704,.1670,.0187,.0127),(4.237,.409,.540,.704,.1282,.0220,.0125),(4.257,.410,.754,.702,.0991,.0202,.0096)],
8:[(3.332,.477,.432,.875,.1036,.0145,.0085),(3.594,.498,.632,.865,.1329,.0134,.0114),(3.734,.508,.845,.860,.1694,.0168,.0128),(3.781,.511,1.068,.858,.1569,.0175,.0087)],
9:[(5.065,.486,.499,.695,.1155,.0150,.0084),(5.444,.517,.734,.685,.1272,.0160,.0116),(5.673,.537,.924,.682,.1117,.0220,.0100),(5.842,.549,1.098,.677,.1386,.0190,.0126)],
}

def _points_in_polygon(x, y, vertices):
    from matplotlib.path import Path as MplPath
    pts = np.column_stack((x, y))

    # IMPORTANT: when Path(..., closed=True) is given an un-repeated vertex
    # list, Matplotlib treats the final supplied vertex as the CLOSEPOLY
    # placeholder and ignores its coordinates.  That made the three-vertex
    # Diehl panel 1 degenerate and also dropped the final physical vertex from
    # every other Diehl polygon.  Explicitly repeat the first vertex so the
    # CLOSEPOLY placeholder is the duplicate rather than a physical vertex.
    verts = np.asarray(vertices, dtype=float)
    closed_verts = np.vstack((verts, verts[0]))
    return MplPath(closed_verts, closed=True).contains_points(
        pts, radius=1e-12
    )


def _assign_rga_bins(events):
    x = np.asarray(events["xB"], dtype=float)
    q2 = np.asarray(events["Q2"], dtype=float)
    mt = np.asarray(events["minus_t"], dtype=float)
    w = np.asarray(events["W"], dtype=float)
    period_index = np.asarray(events["period_index"], dtype=np.int8)
    y = np.empty_like(x)
    for period in PERIODS:
        m = period_index == PERIOD_INDEX[period]
        y[m] = q2[m] / (2.0 * PROTON_MASS_GEV * BEAM_ENERGY_GEV[period] * x[m])
    # endfor
    base = np.isfinite(y) & (q2 > 1.5) & (w > 2.0) & (y < 0.75)
    panel = np.full(x.shape, -1, dtype=np.int16)
    subbin = np.full(x.shape, -1, dtype=np.int16)
    global_bin = np.full(x.shape, -1, dtype=np.int16)
    offset = 0
    mapping = {}
    for ip in range(1,10):
        pmask = base & _points_in_polygon(x, q2, RGA_Q2_XB_POLYGONS[ip])
        edges = RGA_MINUS_T_EDGES[ip]
        for it,(lo,hi) in enumerate(zip(edges[:-1], edges[1:]), start=1):
            g = offset + it
            m = pmask & (mt >= lo) & (mt < hi)
            panel[m], subbin[m], global_bin[m] = ip, it, g
            mapping[g] = (ip,it,lo,hi)
        # endfor
        offset += len(edges)-1
    # endfor
    return panel, subbin, global_bin, y, mapping


def _make_rga_uu_lu_nll(events, run_states, mask):
    """Three-parameter RGA validation fit in the Diehl bin itself.

    This comparison-only likelihood deliberately contains no target-spin terms.
    It floats only

        u1  = F_UU^{cos(phi)}  / F_UU
        u2  = F_UU^{cos(2phi)}/ F_UU
        lu1 = F_LU^{sin(phi)}  / F_UU

    using the RGC events that fall directly in the requested Diehl
    (Q2, xB, -t) bin.  No amplitudes are imported from the nominal 24-bin
    production extraction, and the production likelihood is not modified.
    """
    idx = np.flatnonzero(mask)
    if idx.size == 0:
        raise RuntimeError("Empty RGA comparison bin")
    # endif

    period_idx = events["period_index"][idx].astype(np.int8, copy=False)
    runnum = events["runnum"][idx].astype(np.int32, copy=False)
    helicity = events["helicity"][idx].astype(np.float64, copy=False)
    phi = events["phi"][idx].astype(np.float64, copy=False)
    r_b = events["rB"][idx].astype(np.float64, copy=False)
    r_v = events["rV"][idx].astype(np.float64, copy=False)
    r_w = events["rW"][idx].astype(np.float64, copy=False)

    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    cos_2phi = np.cos(2.0 * phi)

    run_lookup = {
        period: {
            int(run): i
            for i, run in enumerate(run_states[period]["run"])
        }
        for period in PERIODS
    }

    period_data = {}
    for period in PERIODS:
        loc = np.flatnonzero(period_idx == PERIOD_INDEX[period])
        if loc.size == 0:
            period_data[period] = {"loc": loc}
            continue
        # endif

        state = run_states[period]
        observed_state_index = np.fromiter(
            (run_lookup[period][int(run)] for run in runnum[loc]),
            count=loc.size,
            dtype=np.int32,
        )
        event_h = helicity[loc]

        period_data[period] = {
            "loc": loc,
            "h": event_h,
            "observed_charge": np.where(
                event_h > 0.0,
                state["q_plus"][observed_state_index],
                state["q_minus"][observed_state_index],
            ),
            "charge_sum": float(
                np.sum(state["q_plus"] + state["q_minus"])
            ),
            "helicity_charge_sum": float(
                np.sum(state["q_plus"] - state["q_minus"])
            ),
        }
    # endfor

    def nll(u1: float, u2: float, lu1: float) -> float:
        if not all(math.isfinite(v) for v in (u1, u2, lu1)):
            return INVALID_NLL
        # endif

        unpolarized = (
            1.0
            + r_v * u1 * cos_phi
            + r_b * u2 * cos_2phi
        )
        beam = np.empty_like(phi)
        for period in PERIODS:
            loc = period_data[period]["loc"]
            if loc.size == 0:
                continue
            # endif
            beam[loc] = (
                BEAM_POLARIZATION[period]
                * r_w[loc]
                * lu1
                * sin_phi[loc]
            )
        # endfor

        total = 0.0
        for period in PERIODS:
            data = period_data[period]
            loc = data["loc"]
            if loc.size == 0:
                continue
            # endif

            uu = unpolarized[loc]
            bb = beam[loc]
            hh = data["h"]

            factor = uu + hh * bb
            denominator = (
                data["charge_sum"] * uu
                + data["helicity_charge_sum"] * bb
            )
            if (
                np.any(~np.isfinite(factor))
                or np.any(factor <= CROSS_SECTION_FLOOR)
                or np.any(~np.isfinite(denominator))
                or np.any(denominator <= CROSS_SECTION_FLOOR)
            ):
                return INVALID_NLL
            # endif

            probability = data["observed_charge"] * factor / denominator
            if (
                np.any(~np.isfinite(probability))
                or np.any(probability <= 0.0)
                or np.any(probability > 1.0 + 1.0e-10)
            ):
                return INVALID_NLL
            # endif

            total -= float(
                np.sum(np.log(np.maximum(probability, PROBABILITY_FLOOR)))
            )
        # endfor
        return total

    return nll, idx


def _fit_rga_uu_lu(events, run_states, mask):
    """Fit (u1, u2, lu1) and return Minuit plus selected-event indices."""
    from iminuit import Minuit

    nll, idx = _make_rga_uu_lu_nll(events, run_states, mask)
    m = Minuit(
        nll,
        u1=PARAMETER_INITIAL_VALUES["u1"],
        u2=PARAMETER_INITIAL_VALUES["u2"],
        lu1=PARAMETER_INITIAL_VALUES["lu1"],
    )
    m.errordef = Minuit.LIKELIHOOD
    m.print_level = 0
    m.strategy = 1
    for name in ("u1", "u2", "lu1"):
        m.limits[name] = PARAMETER_LIMITS[name]
    # endfor

    # The fit is only three dimensional, so use a deterministic retry from zero
    # if the ordinary start does not give a trustworthy covariance.
    m.migrad()
    m.hesse()
    good = (
        bool(m.valid)
        and m.covariance is not None
        and all(
            np.isfinite(float(m.errors[name])) and float(m.errors[name]) > 0.0
            for name in ("u1", "u2", "lu1")
        )
    )
    if not good:
        # Comparison-only systematic samples can be less well conditioned.
        # Try several deterministic starts, including SIMPLEX recovery, before
        # declaring the Diehl-bin fit invalid.
        starts = (
            (0.0, 0.0, 0.0),
            (float(m.values["u1"]), float(m.values["u2"]), 0.0),
            (0.0, 0.0, float(m.values["lu1"])),
        )
        best = m
        best_fval = float(m.fval) if np.isfinite(float(m.fval)) else np.inf
        for u1_start, u2_start, lu1_start in starts:
            retry = Minuit(
                nll, u1=u1_start, u2=u2_start, lu1=lu1_start
            )
            retry.errordef = Minuit.LIKELIHOOD
            retry.print_level = 0
            retry.strategy = 2
            for name in ("u1", "u2", "lu1"):
                retry.limits[name] = PARAMETER_LIMITS[name]
            # endfor
            retry.migrad(ncall=20000)
            if not retry.valid:
                retry.simplex(ncall=20000)
                retry.migrad(ncall=20000)
            # endif
            retry.hesse()
            retry_good = (
                bool(retry.valid)
                and retry.covariance is not None
                and all(
                    np.isfinite(float(retry.errors[name]))
                    and float(retry.errors[name]) > 0.0
                    for name in ("u1", "u2", "lu1")
                )
            )
            retry_fval = (
                float(retry.fval)
                if np.isfinite(float(retry.fval))
                else np.inf
            )
            if retry_good and retry_fval < best_fval:
                best = retry
                best_fval = retry_fval
            # endif
        # endfor
        m = best
    # endif
    return m, idx



def _rga_fit_cache(events, run_states):
    """Fit all populated Diehl bins for one already-selected RGC cache."""
    panel, subbin, gbin, y, mapping = _assign_rga_bins(events)
    rows = []
    for g, (ip, it, lo, hi) in mapping.items():
        mask = gbin == g
        n = int(np.count_nonzero(mask))
        if n < 20:
            continue
        # endif
        m, idx = _fit_rga_uu_lu(events, run_states, mask)
        valid = (
            bool(m.valid)
            and m.covariance is not None
            and all(
                np.isfinite(float(m.errors[name]))
                and float(m.errors[name]) > 0.0
                for name in ("u1", "u2", "lu1")
            )
        )
        rows.append({
            "rga_panel": ip,
            "t_bin": it,
            "global_bin": g,
            "minus_t_low": lo,
            "minus_t_high": hi,
            "n": n,
            "mean_xB": float(np.mean(events["xB"][idx])),
            "mean_Q2": float(np.mean(events["Q2"][idx])),
            "mean_minus_t": float(np.mean(events["minus_t"][idx])),
            "mean_minus_tprime": float(np.mean(events["minus_tprime"][idx])),
            "mean_y": float(np.mean(y[idx])),
            "mean_epsilon": float(np.nanmean(events["epsilon"][idx])),
            "u1": float(m.values["u1"]),
            "u1_stat": float(m.errors["u1"]),
            "u2": float(m.values["u2"]),
            "u2_stat": float(m.errors["u2"]),
            "lu1": float(m.values["lu1"]),
            "lu1_stat": float(m.errors["lu1"]),
            "fit_valid": valid,
            "edm": float(m.fmin.edm),
        })
    # endfor
    return pd.DataFrame(rows), (panel, subbin, gbin, y, mapping)


def _rga_migration_remap(events, gbin, nominal_production_frame):
    """Event-weight the established LU migration uncertainty into Diehl bins.

    Every selected event already carries its original 24-bin production
    ``bin_number``.  The LU migration uncertainty assigned to that production
    bin is therefore averaged over the events entering each Diehl bin.
    """
    _, assigned, _ = calculate_migration_systematics(nominal_production_frame)
    lu_by_bin = np.asarray(assigned["lu1"], dtype=float)
    result = {}
    for g in np.unique(gbin[gbin > 0]):
        mask = gbin == g
        old_bins = np.asarray(events["bin_number"][mask], dtype=int)
        good = (old_bins >= 1) & (old_bins <= NUMBER_OF_BINS)
        if not np.any(good):
            result[int(g)] = np.nan
            continue
        # endif
        result[int(g)] = float(np.mean(lu_by_bin[old_bins[good] - 1]))
    # endfor
    return result


def _rga_rgc_compatibility(frame):
    """Direct RGC-RGA compatibility using stat+point-to-point errors.

    The Moller beam-polarization scale systematic is common-mode between the
    two measurements and therefore cancels from their relative comparison.  It
    is intentionally not included as either a nuisance parameter or a
    point-to-point uncertainty here.
    """
    use = (
        frame["fit_valid"].astype(bool).to_numpy()
        & np.isfinite(frame["lu1_rgc"].to_numpy(float))
        & np.isfinite(frame["lu1_rga"].to_numpy(float))
        & np.isfinite(frame["uncorr_sigma"].to_numpy(float))
        & (frame["uncorr_sigma"].to_numpy(float) > 0.0)
    )
    residual = (
        frame.loc[use, "lu1_rgc"].to_numpy(float)
        - frame.loc[use, "lu1_rga"].to_numpy(float)
    )
    sigma = frame.loc[use, "uncorr_sigma"].to_numpy(float)
    pulls = residual / sigma
    chi2 = float(np.sum(pulls**2))
    ndf = int(pulls.size)
    try:
        from scipy.stats import chi2 as chi2_distribution
        p_value = float(chi2_distribution.sf(chi2, ndf))
    except Exception:
        p_value = np.nan
    # endtry
    return {
        "chi2": chi2,
        "ndf": ndf,
        "chi2_per_ndf": float(chi2 / ndf) if ndf else np.nan,
        "p_value": p_value,
        "n_points": int(pulls.size),
        "pulls": pulls,
        "use_mask": use,
    }

def _rga_raw_cutflow(args, run_records, output_csv):
    """Raw-ROOT cut flow for every one of the 41 published Diehl bins.

    Diagnostic only: no production or RGA-cross-check selection is changed.
    The table is deliberately evaluated before the selected-event cache so we
    can identify exactly where apparently missing RGC comparison bins vanish.
    """
    cuts = load_channel_cuts(
        args.cut_json.expanduser().resolve(), cut_label="nominal"
    )
    inputs = {p: Path(v) for p, v in DEFAULT_INPUTS.items()}
    for period, path in args.input:
        inputs[period] = path
    # endfor

    stages = (
        "inside_Q2_xB_polygon",
        "Q2_gt_1p5",
        "W_gt_2",
        "y_lt_0p75",
        "inside_Diehl_t_bin",
        "inside_production_xB_tprime",
        "nominal_Mx2",
    )
    keys = [
        (period, ip, it, stage)
        for period in PERIODS
        for ip in range(1, 10)
        for it in range(1, len(RGA_MINUS_T_EDGES[ip]))
        for stage in stages
    ]
    counts = {key: 0 for key in keys}

    for period in PERIODS:
        path = inputs[period].expanduser().resolve()
        print(f"[RGA cutflow] scanning raw nominal tree: {period} {path}", flush=True)
        with uproot.open(path) as root_file:
            tree = resolve_tree(root_file, args.tree, path)
            wanted = {
                key: resolve_branch(tree, BRANCH_ALIASES[key])
                for key in ("xB", "Q2", "W", "t", "tprime", "Mx2", "y")
            }
            expressions = list(dict.fromkeys(wanted.values()))
            for arrays in tree.iterate(
                expressions=expressions,
                step_size=args.chunk_size,
                library="np",
            ):
                x = np.asarray(arrays[wanted["xB"]], float)
                q2 = np.asarray(arrays[wanted["Q2"]], float)
                w = np.asarray(arrays[wanted["W"]], float)
                mt = -np.asarray(arrays[wanted["t"]], float)
                mtp = -np.asarray(arrays[wanted["tprime"]], float)
                mx2 = np.asarray(arrays[wanted["Mx2"]], float)
                y = np.asarray(arrays[wanted["y"]], float)
                _, _, old_bin = bin_indices(x, mtp)

                for ip in range(1, 10):
                    polygon = _points_in_polygon(x, q2, RGA_Q2_XB_POLYGONS[ip])
                    m_q2 = polygon & (q2 > 1.5)
                    m_w = m_q2 & (w > 2.0)
                    m_y = m_w & (y < 0.75)
                    edges = RGA_MINUS_T_EDGES[ip]

                    for it in range(1, len(edges)):
                        counts[(period, ip, it, "inside_Q2_xB_polygon")] += int(np.count_nonzero(polygon))
                        counts[(period, ip, it, "Q2_gt_1p5")] += int(np.count_nonzero(m_q2))
                        counts[(period, ip, it, "W_gt_2")] += int(np.count_nonzero(m_w))
                        counts[(period, ip, it, "y_lt_0p75")] += int(np.count_nonzero(m_y))

                        m_t = m_y & (mt >= edges[it - 1]) & (mt < edges[it])
                        counts[(period, ip, it, "inside_Diehl_t_bin")] += int(np.count_nonzero(m_t))

                        m_prod = m_t & (old_bin >= 1)
                        counts[(period, ip, it, "inside_production_xB_tprime")] += int(np.count_nonzero(m_prod))

                        m_mx = np.zeros_like(m_prod)
                        for old in range(1, NUMBER_OF_BINS + 1):
                            channel_cut = cuts[(period, old)]
                            m_mx |= (
                                m_prod
                                & (old_bin == old)
                                & (mx2 >= channel_cut.low_gev2)
                                & (mx2 < channel_cut.high_gev2)
                            )
                        # endfor
                        counts[(period, ip, it, "nominal_Mx2")] += int(np.count_nonzero(m_mx))
                    # endfor
                # endfor
            # endfor
        # endwith
    # endfor

    rows = []
    global_bin = 0
    for ip in range(1, 10):
        edges = RGA_MINUS_T_EDGES[ip]
        for it in range(1, len(edges)):
            global_bin += 1
            row = {
                "global_bin": global_bin,
                "rga_panel": ip,
                "t_bin": it,
                "minus_t_low": float(edges[it - 1]),
                "minus_t_high": float(edges[it]),
            }
            for stage in stages:
                row[stage] = int(sum(
                    counts[(period, ip, it, stage)] for period in PERIODS
                ))
            # endfor
            # Useful survival fractions that immediately expose the killing cut.
            before_prod = row["inside_Diehl_t_bin"]
            after_prod = row["inside_production_xB_tprime"]
            after_mx = row["nominal_Mx2"]
            row["production_phase_space_survival"] = (
                after_prod / before_prod if before_prod else np.nan
            )
            row["Mx2_survival_after_production_phase_space"] = (
                after_mx / after_prod if after_prod else np.nan
            )
            rows.append(row)
        # endfor
    # endfor

    frame = pd.DataFrame(rows)
    frame.to_csv(output_csv, index=False)
    print("[RGA cutflow] 41-bin raw-event cut flow:", flush=True)
    print(frame.to_string(index=False), flush=True)
    return frame

def _rga_variant_cache_paths(args):
    root = args.output_dir.expanduser().resolve()
    return {
        "nominal": (
            args.cache.expanduser().resolve()
            if args.cache
            else root / "nominal/cache/selected_events.npz"
        ),
        "radiation": root / "isr/cache/selected_events.npz",
        "tight": root / "channel_selection/tight/cache/selected_events.npz",
        "loose": root / "channel_selection/loose/cache/selected_events.npz",
    }


def run_rga_cross_check(args):
    """Final comparison study in the Diehl et al. RGA binning.

    Production extraction is untouched.  This mode:
      * fits (u1,u2,lu1) in each populated Diehl bin;
      * repeats the fit for radiation and tight/loose channel variations;
      * remaps the established LU migration uncertainty event-by-event from
        the original 24-bin (xB,-t') scheme;
      * combines RGC point-to-point systematics in quadrature;
      * compares against published RGA statistical and systematic errors;
      * treats the 2.8% Moller scale systematic as common-mode and therefore
        excludes it from the relative RGA-RGC compatibility test;
      * writes pull distributions and a raw-ROOT 41-bin cut flow.
    """
    out = args.output_dir.expanduser().resolve() / "rga_cross_check"
    tables = out / "tables"
    plots = out / "plots"
    for directory in (out, tables, plots):
        ensure_directory(directory)
    # endfor

    run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    run_states = run_state_arrays(run_records)
    paths = _rga_variant_cache_paths(args)

    if not paths["nominal"].is_file():
        raise FileNotFoundError(
            "RGA cross-check requires the nominal selected-event cache: "
            f"{paths['nominal']}"
        )
    # endif

    print(f"[RGA cross-check] nominal cache: {paths['nominal']}", flush=True)
    nominal_events = load_event_cache(paths["nominal"])
    nominal_fit, assignment = _rga_fit_cache(nominal_events, run_states)
    panel, subbin, gbin, y, mapping = assignment
    print(
        f"[RGA cross-check] nominal Diehl-binned events: "
        f"{np.count_nonzero(gbin > 0):,}",
        flush=True,
    )

    # Production LU migration uncertainty, remapped into the Diehl bins.
    nominal_production_csv = (
        args.output_dir.expanduser().resolve()
        / "nominal/tables/structure_function_ratios.csv"
    )
    if not nominal_production_csv.is_file():
        raise FileNotFoundError(
            "Need the completed nominal production table to remap the "
            f"migration systematic: {nominal_production_csv}"
        )
    # endif
    nominal_production = pd.read_csv(nominal_production_csv)
    migration = _rga_migration_remap(
        nominal_events, gbin, nominal_production
    )

    variant_frames = {}
    for name in ("radiation", "tight", "loose"):
        path = paths[name]
        if not path.is_file():
            raise FileNotFoundError(
                f"RGA final systematic study requires the {name} cache: {path}. "
                "Run the ordinary production/systematic suite first."
            )
        # endif
        print(f"[RGA cross-check] fitting {name} cache: {path}", flush=True)
        events = load_event_cache(path)
        variant_frames[name], _ = _rga_fit_cache(events, run_states)
    # endfor

    def variant_records(frame):
        return {
            (int(row.rga_panel), int(row.t_bin)): row
            for row in frame.itertuples()
        }

    radiation_records = variant_records(variant_frames["radiation"])
    tight_records = variant_records(variant_frames["tight"])
    loose_records = variant_records(variant_frames["loose"])

    def valid_lu(records, key):
        row = records.get(key)
        if row is None or not bool(row.fit_valid):
            return np.nan
        # endif
        value = float(row.lu1)
        return value if np.isfinite(value) else np.nan

    def diagnostic(records, key, field, default=np.nan):
        row = records.get(key)
        return default if row is None else getattr(row, field)

    rows = []
    for row in nominal_fit.itertuples():
        ip, it = int(row.rga_panel), int(row.t_bin)
        pub = RGA_PUBLISHED[ip][it - 1]
        rga_q2, rga_x, rga_t, rga_eps, rga_val, rga_stat, rga_sys = pub
        key = (ip, it)

        radiation_value = valid_lu(radiation_records, key)
        tight_value = valid_lu(tight_records, key)
        loose_value = valid_lu(loose_records, key)

        rad = (
            abs(radiation_value - row.lu1)
            if bool(row.fit_valid) and np.isfinite(radiation_value)
            else np.nan
        )
        if (
            bool(row.fit_valid)
            and np.isfinite(tight_value)
            and np.isfinite(loose_value)
        ):
            dt = tight_value - row.lu1
            dl = loose_value - row.lu1
            channel = math.sqrt(0.5 * (dt * dt + dl * dl))
        else:
            channel = np.nan
        # endif

        migration_sys = float(migration.get(int(row.global_bin), np.nan))
        rgc_ptp = (
            math.sqrt(rad**2 + channel**2 + migration_sys**2)
            if np.all(np.isfinite([rad, channel, migration_sys]))
            else np.nan
        )
        uncorr_sigma = (
            math.sqrt(
                row.lu1_stat**2 + rgc_ptp**2 + rga_stat**2 + rga_sys**2
            )
            if bool(row.fit_valid) and np.isfinite(rgc_ptp)
            else np.nan
        )
        raw_delta = row.lu1 - rga_val
        raw_pull = (
            raw_delta / uncorr_sigma
            if np.isfinite(uncorr_sigma) and uncorr_sigma > 0.0
            else np.nan
        )

        rows.append({
            "rga_panel": ip,
            "t_bin": it,
            "global_bin": int(row.global_bin),
            "n_rgc": int(row.n),
            "mean_xB_rgc": row.mean_xB,
            "mean_Q2_rgc": row.mean_Q2,
            "mean_minus_t_rgc": row.mean_minus_t,
            "mean_minus_tprime_rgc": row.mean_minus_tprime,
            "mean_y_rgc": row.mean_y,
            "mean_epsilon_rgc": row.mean_epsilon,
            "u1_rgc": row.u1,
            "stat_u1_rgc": row.u1_stat,
            "u2_rgc": row.u2,
            "stat_u2_rgc": row.u2_stat,
            "lu1_rgc": row.lu1,
            "stat_rgc": row.lu1_stat,
            "radiation_rgc": rad,
            "channel_selection_rgc": channel,
            "migration_remap_rgc": migration_sys,
            "ptp_rgc": rgc_ptp,
            "fit_valid": bool(row.fit_valid),
            "edm": row.edm,
            "radiation_n": int(diagnostic(radiation_records, key, "n", 0)),
            "radiation_fit_valid": bool(diagnostic(
                radiation_records, key, "fit_valid", False
            )),
            "radiation_edm": float(diagnostic(radiation_records, key, "edm")),
            "radiation_lu1": float(diagnostic(radiation_records, key, "lu1")),
            "radiation_lu1_stat": float(diagnostic(
                radiation_records, key, "lu1_stat"
            )),
            "tight_n": int(diagnostic(tight_records, key, "n", 0)),
            "tight_fit_valid": bool(diagnostic(
                tight_records, key, "fit_valid", False
            )),
            "tight_edm": float(diagnostic(tight_records, key, "edm")),
            "tight_lu1": float(diagnostic(tight_records, key, "lu1")),
            "tight_lu1_stat": float(diagnostic(
                tight_records, key, "lu1_stat"
            )),
            "loose_n": int(diagnostic(loose_records, key, "n", 0)),
            "loose_fit_valid": bool(diagnostic(
                loose_records, key, "fit_valid", False
            )),
            "loose_edm": float(diagnostic(loose_records, key, "edm")),
            "loose_lu1": float(diagnostic(loose_records, key, "lu1")),
            "loose_lu1_stat": float(diagnostic(
                loose_records, key, "lu1_stat"
            )),
            "radiation_fallback_used": False,
            "channel_fallback_used": False,
            "xB_rga": rga_x,
            "Q2_rga": rga_q2,
            "minus_t_rga": rga_t,
            "epsilon_rga": rga_eps,
            "lu1_rga": rga_val,
            "stat_rga": rga_stat,
            "ptp_rga": rga_sys,
            "uncorr_sigma": uncorr_sigma,
            "delta_rgc_minus_rga": raw_delta,
            "raw_pull_with_ptp": raw_pull,
        })
    # endfor

    frame = pd.DataFrame(rows)

    # Last-resort comparison-only prescription: after the robust retries above,
    # copy any still-missing systematic component from the closest available
    # adjacent -t bin in the SAME Diehl Q2-xB panel.  Central values and
    # statistical errors are never replaced.  Every fallback is recorded.
    frame["radiation_rgc_fallback_source_t_bin"] = np.nan
    frame["channel_selection_rgc_fallback_source_t_bin"] = np.nan

    def fill_adjacent_systematic(column, flag_column, source_column):
        missing = frame.index[~np.isfinite(frame[column].to_numpy(float))]
        for idx in missing:
            panel = int(frame.at[idx, "rga_panel"])
            tbin = int(frame.at[idx, "t_bin"])
            candidates = frame[
                (frame["rga_panel"] == panel)
                & np.isfinite(frame[column])
            ].copy()
            if candidates.empty:
                continue
            # endif
            candidates["_distance"] = np.abs(
                candidates["t_bin"].astype(int) - tbin
            )
            source_idx = candidates.sort_values(
                ["_distance", "t_bin"]
            ).index[0]
            frame.at[idx, column] = float(frame.at[source_idx, column])
            frame.at[idx, flag_column] = True
            frame.at[idx, source_column] = int(
                frame.at[source_idx, "t_bin"]
            )
        # endfor

    fill_adjacent_systematic(
        "radiation_rgc",
        "radiation_fallback_used",
        "radiation_rgc_fallback_source_t_bin",
    )
    fill_adjacent_systematic(
        "channel_selection_rgc",
        "channel_fallback_used",
        "channel_selection_rgc_fallback_source_t_bin",
    )

    # Recompute the total point-to-point uncertainty after any fallback.
    frame["ptp_rgc"] = np.sqrt(
        frame["radiation_rgc"]**2
        + frame["channel_selection_rgc"]**2
        + frame["migration_remap_rgc"]**2
    )
    frame["uncorr_sigma"] = np.sqrt(
        frame["stat_rgc"]**2
        + frame["ptp_rgc"]**2
        + frame["stat_rga"]**2
        + frame["ptp_rga"]**2
    )
    frame["raw_pull_with_ptp"] = (
        frame["delta_rgc_minus_rga"] / frame["uncorr_sigma"]
    )

    compatibility = _rga_rgc_compatibility(frame)
    frame["pull"] = np.nan
    frame.loc[compatibility["use_mask"], "pull"] = compatibility["pulls"]

    csv = tables / "rga_rgc_final_comparison.csv"
    frame.to_csv(csv, index=False)

    # Compact statistics/systematics comparison table.
    summary_rows = []
    for label, stat_col, ptp_col in (
        ("RGC", "stat_rgc", "ptp_rgc"),
        ("RGA", "stat_rga", "ptp_rga"),
    ):
        for kind, column in (("stat", stat_col), ("ptp", ptp_col)):
            values = pd.to_numeric(frame[column], errors="coerce").to_numpy(float)
            values = values[np.isfinite(values)]
            summary_rows.append({
                "dataset": label,
                "uncertainty": kind,
                "n": int(values.size),
                "mean": float(np.mean(values)) if values.size else np.nan,
                "median": float(np.median(values)) if values.size else np.nan,
                "rms": float(np.sqrt(np.mean(values**2))) if values.size else np.nan,
                "max": float(np.max(values)) if values.size else np.nan,
            })
        # endfor
    # endfor
    uncertainty_summary = pd.DataFrame(summary_rows)
    uncertainty_summary.to_csv(
        tables / "rga_rgc_uncertainty_summary.csv", index=False
    )

    # Raw ROOT cut flow for the panel-1 investigation.
    cutflow = _rga_raw_cutflow(
        args, run_records, tables / "rga_41bin_cutflow.csv"
    )

    valid_pulls = frame["pull"].to_numpy(float)
    valid_pulls = valid_pulls[np.isfinite(valid_pulls)]
    pull_stats = {
        "n": int(valid_pulls.size),
        "mean": float(np.mean(valid_pulls)) if valid_pulls.size else np.nan,
        "sample_sigma": (
            float(np.std(valid_pulls, ddof=1))
            if valid_pulls.size > 1 else np.nan
        ),
        "rms": (
            float(np.sqrt(np.mean(valid_pulls**2)))
            if valid_pulls.size else np.nan
        ),
        "median": (
            float(np.median(valid_pulls)) if valid_pulls.size else np.nan
        ),
    }

    fallback_rows = frame[
        frame["radiation_fallback_used"]
        | frame["channel_fallback_used"]
    ]
    if len(fallback_rows):
        print(
            "\n[RGA cross-check] adjacent-bin systematic fallback usage:",
            flush=True,
        )
        print(
            fallback_rows[[
                "rga_panel", "t_bin",
                "radiation_fallback_used",
                "radiation_rgc_fallback_source_t_bin",
                "channel_fallback_used",
                "channel_selection_rgc_fallback_source_t_bin",
            ]].to_string(index=False),
            flush=True,
        )
    # endif

    print("\n[RGA cross-check] uncertainty summary", flush=True)
    print(uncertainty_summary.to_string(index=False), flush=True)
    print(
        "[RGA cross-check] direct compatibility (common Moller scale cancels): "
        f"chi2/ndf={compatibility['chi2']:.3f}/{compatibility['ndf']}="
        f"{compatibility['chi2_per_ndf']:.3f}, "
        f"p={compatibility['p_value']:.4g}",
        flush=True,
    )
    print(
        "[RGA cross-check] pull distribution: "
        f"N={pull_stats['n']}, mean={pull_stats['mean']:+.3f}, "
        f"sigma={pull_stats['sample_sigma']:.3f}, "
        f"RMS={pull_stats['rms']:.3f}",
        flush=True,
    )

    if not args.skip_plots:
        valid = frame[
            frame["fit_valid"]
            & np.isfinite(frame["ptp_rgc"])
        ].copy()

        # Final 3x3 overlay: inner bars are statistical, outer bars are
        # stat+point-to-point.  The common 2.8% Moller scale cancels from this
        # comparison.  The legend identifies datasets only; the nested error
        # bars make the statistical versus total point-to-point treatment clear.
        from matplotlib.lines import Line2D

        fig, axes = plt.subplots(3, 3, figsize=(12, 10), sharey=True)
        for ip, ax in enumerate(axes.flat, start=1):
            data = valid[valid.rga_panel == ip]
            pub = np.asarray(RGA_PUBLISHED[ip], float)
            rga_outer = np.sqrt(pub[:,5]**2 + pub[:,6]**2)
            ax.errorbar(
                pub[:,2], pub[:,4], yerr=rga_outer, fmt="o",
                capsize=2, color="tab:blue"
            )
            ax.errorbar(
                pub[:,2], pub[:,4], yerr=pub[:,5], fmt="none",
                capsize=4, color="tab:blue"
            )
            if len(data):
                rgc_outer = np.sqrt(data.stat_rgc**2 + data.ptp_rgc**2)
                ax.errorbar(
                    data.mean_minus_t_rgc, data.lu1_rgc,
                    yerr=rgc_outer, fmt="s", capsize=2,
                    color="tab:orange"
                )
                ax.errorbar(
                    data.mean_minus_t_rgc, data.lu1_rgc,
                    yerr=data.stat_rgc, fmt="none", capsize=4,
                    color="tab:orange"
                )
            # endif
            ax.axhline(0.0, lw=0.8)
            ax.set_ylim(-0.1, 0.3)
            ax.set_title(f"RGA $Q^2$-$x_B$ bin {ip}")
            ax.set_xlabel(r"$-t$ (GeV$^2$)")
            if ip in (1,4,7):
                ax.set_ylabel(
                    r"$F_{LU}^{\sin\phi}/F_{UU}"
                    r"=\sigma_{LT'}/\sigma_0$"
                )
            # endif
        # endfor

        # Hard-coded two-entry dataset legend so it is independent of whether
        # a particular panel contains RGC points.
        legend_ax = axes.flat[0]
        legend_ax.legend(
            handles=[
                Line2D(
                    [0], [0], marker="o", linestyle="none",
                    color="tab:blue", label="RGA"
                ),
                Line2D(
                    [0], [0], marker="s", linestyle="none",
                    color="tab:orange", label="RGC"
                ),
            ],
            fontsize=8,
        )
        fig.tight_layout()
        fig.savefig(plots / "rga_rgc_final_overlay.png", dpi=200)
        plt.close(fig)

        # Point-by-point statistical and point-to-point uncertainty comparison.
        fig, ax = plt.subplots(figsize=(11, 5.5))
        xplot = np.arange(len(valid))
        ax.plot(xplot, valid.stat_rgc, "o", label="RGC stat")
        ax.plot(xplot, valid.ptp_rgc, "s", label="RGC point-to-point")
        ax.plot(xplot, valid.stat_rga, "^", label="RGA stat")
        ax.plot(xplot, valid.ptp_rga, "v", label="RGA point-to-point")
        ax.set_xlabel("Matched RGA/RGC comparison point")
        ax.set_ylabel(r"Absolute uncertainty on $F_{LU}^{\sin\phi}/F_{UU}$")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots / "rga_rgc_uncertainty_comparison.png", dpi=200)
        plt.close(fig)

        # Source breakdown of the RGC point-to-point systematic.
        fig, ax = plt.subplots(figsize=(11, 5.5))
        ax.plot(xplot, valid.radiation_rgc, "o", label="Radiation")
        ax.plot(xplot, valid.channel_selection_rgc, "s", label="Channel selection")
        ax.plot(xplot, valid.migration_remap_rgc, "^", label="Migration remap")
        ax.plot(xplot, valid.ptp_rgc, "D", label="Total RGC point-to-point")
        ax.set_xlabel("Matched RGA/RGC comparison point")
        ax.set_ylabel("Absolute systematic uncertainty")
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots / "rgc_ptp_systematic_breakdown.png", dpi=200)
        plt.close(fig)

        # Pull histogram for the direct RGC-RGA comparison.
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        bins = np.linspace(-4.0, 4.0, 17)
        ax.hist(valid_pulls, bins=bins, histtype="stepfilled", alpha=0.65)
        ax.axvline(0.0, lw=1.0)
        ax.set_xlabel("RGC - RGA pull")
        ax.set_ylabel("Comparison points")
        ax.text(
            0.03, 0.97,
            f"N = {pull_stats['n']}\n"
            f"mean = {pull_stats['mean']:+.2f}\n"
            f"sigma = {pull_stats['sample_sigma']:.2f}\n"
            f"RMS = {pull_stats['rms']:.2f}\n"
            f"chi2/ndf = {compatibility['chi2_per_ndf']:.2f}",
            transform=ax.transAxes, va="top"
        )
        fig.tight_layout()
        fig.savefig(plots / "rga_rgc_pull_histogram.png", dpi=200)
        plt.close(fig)

        # Pull versus matched comparison point.
        fig, ax = plt.subplots(figsize=(11, 5.5))
        ax.axhline(0.0, lw=0.8)
        ax.axhline(+1.0, lw=0.7, ls="--")
        ax.axhline(-1.0, lw=0.7, ls="--")
        ax.plot(np.arange(len(valid_pulls)), valid_pulls, "o")
        ax.set_xlabel("Matched RGA/RGC comparison point")
        ax.set_ylabel("RGC - RGA pull")
        fig.tight_layout()
        fig.savefig(plots / "rga_rgc_pulls.png", dpi=200)
        plt.close(fig)
    # endif

    manifest = {
        "selection": {
            "Q2_min_gev2": 1.5,
            "W_min_gev": 2.0,
            "y_max": 0.75,
        },
        "fit_parameters": [
            "u1 = F_UU^{cos(phi)}/F_UU",
            "u2 = F_UU^{cos(2phi)}/F_UU",
            "lu1 = F_LU^{sin(phi)}/F_UU",
        ],
        "rgc_point_to_point": (
            "sqrt(radiation^2 + channel_selection^2 + "
            "migration_remap^2)"
        ),
        "radiation": "|LU_ISR - LU_nominal| in the Diehl binning",
        "channel_selection": (
            "sqrt(((LU_tight-LU_nominal)^2 + "
            "(LU_loose-LU_nominal)^2)/2) in the Diehl binning"
        ),
        "migration": (
            "event-weighted remap of the established production-bin LU "
            "migration systematic using each event's original (xB,-t') bin"
        ),
        "failed_variation_fit_fallback": (
            "After deterministic Minuit retries, any still-missing radiation "
            "or channel-selection systematic is copied from the closest "
            "available adjacent -t bin within the same Diehl Q2-xB panel. "
            "Fallback use and source t-bin are recorded explicitly. Nominal "
            "central values and statistical uncertainties are never replaced."
        ),
        "rga_point_to_point": (
            "published Diehl et al. systematic uncertainty for each point"
        ),
        "beam_polarization_scale": {
            "fractional_uncertainty": 0.028,
            "treatment_in_comparison": (
                "common-mode between RGA and RGC; cancels from the relative "
                "comparison and is not included in chi2 or pulls"
            ),
        },
        "compatibility": {
            key: value for key, value in compatibility.items()
            if key not in ("pulls", "use_mask")
        },
        "pull_distribution": pull_stats,
        "products": {
            "comparison_csv": str(csv),
            "uncertainty_summary_csv": str(
                tables / "rga_rgc_uncertainty_summary.csv"
            ),
            "raw_cutflow_csv": str(tables / "rga_41bin_cutflow.csv"),
            "plots_directory": str(plots),
        },
        "important_note": (
            "This entire study is isolated behind --rga-cross-check. "
            "The production extraction is unchanged."
        ),
    }
    write_json(out / "rga_cross_check_manifest.json", manifest)
    print(f"[RGA cross-check] wrote {csv}", flush=True)
    return 0



def fit_xb_integrated_tprime_zero_uu(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    x_index: int,
) -> dict[str, Any]:
    """Nominal unbinned likelihood integrated over all six -t' bins in one xB bin.

    Each original (xB,-t') bin keeps its own period-dependent dilution-factor
    nuisance parameters.  Only the five polarized physics amplitudes are shared
    across -t', and u1=u2=0 are fixed.  This is a diagnostic of statistical
    stability, not a replacement for the production 24-bin result.
    """
    t_count = len(MINUS_TPRIME_BINS_GEV2)
    bin_numbers = [
        combined_bin_number(x_index, t_index)
        for t_index in range(t_count)
    ]
    sub_nlls = {
        bin_number: make_bin_nll(
            events, run_states, dilution_records, bin_number, "nominal"
        )[0]
        for bin_number in bin_numbers
    }

    physics_names = list(PHYSICS_PARAMETERS)
    dilution_names = [
        f"f_{period}_bin{bin_number:02d}"
        for bin_number in bin_numbers
        for period in PERIODS
    ]
    names = physics_names + dilution_names

    start = [float(PARAMETER_INITIAL_VALUES[name]) for name in physics_names]
    for bin_number in bin_numbers:
        for period in PERIODS:
            start.append(float(dilution_records[(period, bin_number)].value))
        # endfor
    # endfor

    def combined_nll(*pars: float) -> float:
        values = dict(zip(names, pars))
        total = 0.0
        for bin_number in bin_numbers:
            total += sub_nlls[bin_number](
                values["u1"], values["u2"], values["lu1"],
                values["ul1"], values["ul2"], values["ll0"], values["ll1"],
                values[f"f_su22_bin{bin_number:02d}"],
                values[f"f_fa22_bin{bin_number:02d}"],
                values[f"f_sp23_bin{bin_number:02d}"],
            )
        # endfor
        return float(total)

    def configured_minuit(start_values: list[float]) -> Minuit:
        candidate = Minuit(combined_nll, *start_values, name=names)
        candidate.errordef = Minuit.LIKELIHOOD
        candidate.print_level = 0
        candidate.strategy = 1
        for name in physics_names:
            candidate.limits[name] = PARAMETER_LIMITS[name]
        # endfor
        candidate.values["u1"] = 0.0
        candidate.values["u2"] = 0.0
        candidate.fixed["u1"] = True
        candidate.fixed["u2"] = True
        for bin_number in bin_numbers:
            for period in PERIODS:
                name = f"f_{period}_bin{bin_number:02d}"
                record = dilution_records[(period, bin_number)]
                width = max(
                    8.0 * record.stat_uncertainty,
                    0.20 * record.value,
                    0.02,
                )
                candidate.limits[name] = (
                    max(1.0e-6, record.value - width),
                    record.value + width,
                )
                if record.stat_uncertainty == 0.0:
                    candidate.fixed[name] = True
                # endif
            # endfor
        # endfor
        return candidate

    # Two deterministic starts are sufficient for this diagnostic and avoid
    # multiplying the already-larger integrated likelihood cost unnecessarily.
    starts = [list(start)]
    zero_physics = list(start)
    for i, name in enumerate(physics_names):
        if name not in ("u1", "u2"):
            zero_physics[i] = 0.0
        # endif
    # endfor
    starts.append(zero_physics)

    attempted = []
    for start_values in starts:
        candidate = configured_minuit(start_values)
        candidate.migrad(ncall=80000)
        if not candidate.fmin.is_valid:
            candidate.simplex(ncall=40000)
            candidate.strategy = 2
            candidate.migrad(ncall=120000)
        # endif
        candidate.hesse()
        attempted.append(candidate)
    # endfor

    def quality(candidate: Minuit) -> tuple[int, float, float]:
        valid = candidate.fmin.is_valid and math.isfinite(float(candidate.fval))
        return (
            0 if valid else 1,
            float(candidate.fval) if math.isfinite(float(candidate.fval)) else math.inf,
            float(candidate.fmin.edm) if math.isfinite(float(candidate.fmin.edm)) else math.inf,
        )

    best = min(attempted, key=quality)
    values = {name: float(best.values[name]) for name in PHYSICS_PARAMETERS}
    errors = {name: float(best.errors[name]) for name in PHYSICS_PARAMETERS}
    return {
        "x_index": x_index,
        "xB_low": float(XB_BINS[x_index][0]),
        "xB_high": float(XB_BINS[x_index][1]),
        "xB_center": 0.5 * float(XB_BINS[x_index][0] + XB_BINS[x_index][1]),
        "bin_numbers": bin_numbers,
        "values": values,
        "errors": errors,
        "valid": bool(best.fmin.is_valid),
        "edm": float(best.fmin.edm),
        "fval": float(best.fval),
        "nfcn": int(best.fmin.nfcn),
    }


def run_xb_integrated_tprime_zero_uu_study(
    cache_path: Path,
    run_info_path: Path,
    dilution_json_path: Path,
    output_dir: Path,
) -> Path:
    """Run and write the four xB-only, t'-integrated nominal baseline fits."""
    print("=" * 78, flush=True)
    print("[xB-integrated] START: nominal unbinned likelihood, u1=u2=0", flush=True)
    print("=" * 78, flush=True)
    events = load_event_cache(cache_path)
    run_records = parse_run_info_csv(run_info_path)
    run_states = run_state_arrays(run_records)
    dilution_records = load_dilution_factors(dilution_json_path, cut_label="nominal")

    rows = []
    details = []
    for x_index in range(len(XB_BINS)):
        t0 = time.perf_counter()
        fit = fit_xb_integrated_tprime_zero_uu(
            events, run_states, dilution_records, x_index
        )
        elapsed = time.perf_counter() - t0
        row = {
            "x_index": x_index,
            "xB_low": fit["xB_low"],
            "xB_high": fit["xB_high"],
            "xB_center": fit["xB_center"],
            "fit_valid": fit["valid"],
            "edm": fit["edm"],
            "fval": fit["fval"],
            "nfcn": fit["nfcn"],
            "elapsed_seconds": elapsed,
        }
        for name in PUBLISHED_SYSTEMATIC_PARAMETERS:
            row[name] = fit["values"][name]
            row[f"{name}_error"] = fit["errors"][name]
        # endfor
        rows.append(row)
        details.append(fit)
        print(
            f"[xB-integrated] x bin {x_index + 1}/{len(XB_BINS)} DONE | "
            f"{elapsed:.1f} s | valid={fit['valid']} | EDM={fit['edm']:.3e}",
            flush=True,
        )
    # endfor

    tables = output_dir / "tables"
    ensure_directory(tables)
    csv_path = tables / "structure_function_ratios_xB_integrated_tprime.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    write_json(
        tables / "structure_function_ratios_xB_integrated_tprime_details.json",
        details,
    )
    print(f"[xB-integrated] wrote {csv_path}", flush=True)
    return csv_path




# =============================================================================
# Targeted run-period stability diagnostics
# =============================================================================

PERIOD_STABILITY_DIAGNOSTIC_BINS: tuple[int, ...] = (4, 6, 7, 14, 17, 18, 19, 21)
PERIOD_STABILITY_PROFILE_PARAMETERS: dict[int, str] = {
    7: "ll0",
    14: "ll1",
    17: "ll1",
    19: "ll1",
    21: "ll1",
}
PERIOD_STABILITY_2D_PROFILE_BINS: tuple[int, ...] = (7, 19)
COMBINED_PERIOD_KEY = "combined"


def parse_bin_list(value: str) -> tuple[int, ...]:
    bins = tuple(sorted({int(item.strip()) for item in value.split(",") if item.strip()}))
    if not bins or any(bin_number < 1 or bin_number > NUMBER_OF_BINS for bin_number in bins):
        raise argparse.ArgumentTypeError(f"Bins must be comma-separated integers in 1--{NUMBER_OF_BINS}.")
    # endif
    return bins


def _diagnostic_period_fit_worker(task: dict[str, Any]) -> dict[str, Any]:
    """Fit one period/bin and optionally profile one physics parameter."""
    if (_WORKER_EVENTS is None or _WORKER_RUN_STATES is None
            or _WORKER_DILUTION_RECORDS is None):
        raise RuntimeError("Diagnostic worker was not initialized.")
    # endif
    bin_number = int(task["bin_number"])
    period = str(task["period"])
    active_periods = PERIODS if period == COMBINED_PERIOD_KEY else (period,)
    parameter = task.get("profile_parameter")
    points = int(task.get("profile_points", 21))
    fit = fit_one_variant(
        _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
        bin_number, "nominal", active_periods=active_periods,
        fixed_physics_parameters={"u1": 0.0, "u2": 0.0},
    )
    result: dict[str, Any] = {
        "bin_number": bin_number, "period": period, "fit": fit,
        "profile_parameter": parameter, "profile": None,
    }
    if parameter is not None and fit["valid"]:
        center = float(fit["values"][parameter])
        sigma = float(fit["errors"][parameter])
        low_limit, high_limit = PARAMETER_LIMITS[parameter]
        half_width = max(4.0 * sigma, 0.20)
        low = max(low_limit, center - half_width)
        high = min(high_limit, center + half_width)
        grid = np.linspace(low, high, points)
        profile_nll: list[float] = []
        start = dict(fit["values"])
        for value in grid:
            profiled = fit_one_variant(
                _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
                bin_number, "nominal", active_periods=active_periods,
                initial_values=start,
                fixed_physics_parameters={"u1": 0.0, "u2": 0.0, parameter: float(value)},
            )
            profile_nll.append(float(profiled["minimum_nll"]))
            if profiled["valid"]:
                start = dict(profiled["values"])
            # endif
        # endfor
        minimum = min(profile_nll)
        result["profile"] = {
            "grid": grid.tolist(),
            "two_delta_nll": [2.0 * (value - minimum) for value in profile_nll],
        }
    # endif
    return result


def _charge_normalized_state_yields(
    events: Mapping[str, np.ndarray], run_states: Mapping[str, Mapping[str, np.ndarray]],
    bin_number: int, period: str, phi_edges: np.ndarray,
) -> dict[str, Any]:
    """Return raw charge-normalized phi yields for the four target/helicity states."""
    period_mask = (
        (events["bin_number"] == bin_number)
        & (events["period_index"] == PERIOD_INDEX[period])
    )
    runnum = events["runnum"][period_mask].astype(np.int32, copy=False)
    helicity = events["helicity"][period_mask].astype(np.int8, copy=False)
    phi = events["phi"][period_mask].astype(np.float64, copy=False)
    state = run_states[period]
    lookup = {int(run): i for i, run in enumerate(state["run"])}
    event_state = np.asarray([lookup[int(run)] for run in runnum], dtype=np.int32)
    event_target_sign = np.sign(state["pt"][event_state]).astype(np.int8)
    out: dict[str, Any] = {}
    for target_sign in (-1, 1):
        run_mask = np.sign(state["pt"]) == target_sign
        for h in (-1, 1):
            charge = float(np.sum((state["q_plus"] if h > 0 else state["q_minus"])[run_mask]))
            selected = (event_target_sign == target_sign) & (helicity == h)
            counts, _ = np.histogram(phi[selected], bins=phi_edges)
            if charge > 0.0:
                yields = counts.astype(float) / charge
                errors = np.sqrt(counts.astype(float)) / charge
            else:
                yields = np.full(counts.shape, np.nan)
                errors = np.full(counts.shape, np.nan)
            # endif
            out[f"t{target_sign:+d}_h{h:+d}"] = {
                "charge": charge, "counts": counts.tolist(),
                "yield": yields.tolist(), "error": errors.tolist(),
            }
        # endfor
    # endfor
    return out



def run_flagged_exclusivity_window_diagnostic(
    args: argparse.Namespace,
    root: Path,
    bins: tuple[int, ...],
    profile_map: Mapping[int, tuple[str, ...]],
    workers: int,
) -> Path:
    """Refit flagged period-stability bins with 1, 2, and 3 sigma Mx2 windows."""
    out_dir = root / "period_stability" / "diagnostics" / "period_excursions"
    ensure_directory(out_dir)
    run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    run_states = run_state_arrays(run_records)
    run_state_payload = {
        period: {key: value.tolist() for key, value in state.items()}
        for period, state in run_states.items()
    }

    dilution_path = (
        args.dilution_json.expanduser().resolve() if args.dilution_json
        else find_default_dilution_json(args.dilution_dir.expanduser().resolve()).resolve()
    )
    loose_cache = root / "channel_selection" / "loose" / "cache" / "selected_events.npz"
    if not loose_cache.is_file():
        print("[period window diagnostic] loose 3-sigma cache missing; building it from ROOT...", flush=True)
        inputs = {period: Path(value) for period, value in DEFAULT_INPUTS.items()}
        for period, path in args.input:
            inputs[period] = path
        # endfor
        loose_cuts = load_channel_cuts(args.cut_json.expanduser().resolve(), cut_label="loose")
        build_event_cache(
            input_paths=inputs, tree_name=args.tree, chunk_size=args.chunk_size,
            run_records=run_records, cuts=loose_cuts, cache_path=loose_cache,
        )
    # endif

    cache_by_window = {
        "loose": loose_cache,
        "nominal": root / "period_stability" / "cache" / "selected_events.npz",
        "tight": out_dir / "cache" / "selected_events_tight.npz",
    }
    if not cache_by_window["tight"].is_file():
        tight_cuts = load_channel_cuts(args.cut_json.expanduser().resolve(), cut_label="tight")
        derive_event_cache(
            source_cache_path=loose_cache, cuts=tight_cuts,
            cache_path=cache_by_window["tight"],
        )
    # endif

    all_rows = []
    for window in ("tight", "nominal", "loose"):
        dilution_records = load_dilution_factors(dilution_path, cut_label=window)
        dilution_payload = {
            period: {
                str(bin_number): {
                    "x_index": record.x_index, "t_index": record.t_index,
                    "value": record.value, "stat_uncertainty": record.stat_uncertainty,
                }
                for (record_period, bin_number), record in dilution_records.items()
                if record_period == period
            }
            for period in PERIODS
        }
        tasks = [
            {"bin_number": bin_number, "period": period, "profile_parameter": None}
            for bin_number in bins for period in (*PERIODS, COMBINED_PERIOD_KEY)
        ]
        results = []
        mp_context = mp.get_context("spawn")
        with ProcessPoolExecutor(
            max_workers=min(workers, len(tasks)), mp_context=mp_context,
            initializer=initialize_fit_worker,
            initargs=(str(cache_by_window[window]), run_state_payload, dilution_payload),
        ) as executor:
            futures = [executor.submit(_diagnostic_period_fit_worker, task) for task in tasks]
            for index, future in enumerate(as_completed(futures), start=1):
                results.append(future.result())
                print(f"[period window diagnostic:{window}] {index}/{len(tasks)}", flush=True)
            # endfor
        # endwith
        by_key = {(item["bin_number"], item["period"]): item["fit"] for item in results}
        for bin_number in bins:
            combined_fit = by_key[(bin_number, COMBINED_PERIOD_KEY)]
            for parameter in profile_map.get(bin_number, tuple()):
                combined = float(combined_fit["values"][parameter])
                combined_err = float(combined_fit["errors"][parameter])
                for period in PERIODS:
                    fit = by_key[(bin_number, period)]
                    value = float(fit["values"][parameter])
                    error = float(fit["errors"][parameter])
                    denom2 = error**2 - combined_err**2
                    pull = (value - combined) / math.sqrt(denom2) if denom2 > 0.0 else math.nan
                    all_rows.append({
                        "window": window, "bin_number": bin_number, "parameter": parameter,
                        "period": period, "value": value, "stat": error,
                        "simultaneous_value": combined, "simultaneous_stat": combined_err,
                        "pull_wrt_simultaneous": pull,
                        "events": int(fit["metadata"]["number_of_events"]),
                        "fit_valid": fit["valid"], "accurate_covariance": fit["accurate_covariance"],
                        "positive_definite_covariance": fit["positive_definite_covariance"],
                        "parameters_at_limit": fit["parameters_at_limit"], "edm": fit["edm"],
                    })
                # endfor
            # endfor
        # endfor
    # endfor

    frame = pd.DataFrame(all_rows)
    path = out_dir / "flagged_bins_exclusivity_window_period_refits.csv"
    frame.to_csv(path, index=False)

    plot_dir = out_dir / "plots" / "exclusivity_window"
    ensure_directory(plot_dir)
    window_offsets = {"tight": -0.18, "nominal": 0.0, "loose": 0.18}
    for (bin_number, parameter), group in frame.groupby(["bin_number", "parameter"]):
        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        xbase = {period: i for i, period in enumerate(PERIODS)}
        for window in ("tight", "nominal", "loose"):
            subset = group.loc[group["window"] == window]
            xs = [xbase[p] + window_offsets[window] for p in subset["period"]]
            ax.errorbar(xs, subset["value"], yerr=subset["stat"], marker="o", linestyle="none", capsize=3, label=window.capitalize())
        # endfor
        ax.set_xticks(range(len(PERIODS)), [PERIOD_LABELS[p] for p in PERIODS])
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_title(f"Bin {int(bin_number)}: exclusivity-window dependence of period fits")
        ax.grid(alpha=0.25); ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(plot_dir / f"bin_{int(bin_number):02d}_{parameter}_window_period_fits.png", dpi=180)
        plt.close(fig)
    # endfor
    return path

def run_period_stability_diagnostics(args: argparse.Namespace, root: Path, workers: int) -> int:
    """Run focused diagnostics for the largest period-stability excursions."""
    nominal_dir = root / "period_stability"
    cache_path = (
        args.cache.expanduser().resolve() if args.cache
        else nominal_dir / "cache/selected_events.npz"
    )
    if not cache_path.is_file():
        raise FileNotFoundError(
            f"Period-stability cache not found: {cache_path}. Run --period-stability-only first."
        )
    # endif
    run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    run_states = run_state_arrays(run_records)
    dilution_path = (
        args.dilution_json.expanduser().resolve() if args.dilution_json
        else find_default_dilution_json(args.dilution_dir.expanduser().resolve()).resolve()
    )
    dilution_records = load_dilution_factors(dilution_path, cut_label="nominal")
    cuts = load_channel_cuts(args.cut_json.expanduser().resolve(), cut_label="nominal")
    events = load_event_cache(cache_path)
    bins = args.diagnostic_bins
    flagged_map = getattr(args, "_period_stability_flagged_parameters", None)
    profile_map = (flagged_map if flagged_map is not None else
                   {b: ((PERIOD_STABILITY_PROFILE_PARAMETERS[b],) if b in PERIOD_STABILITY_PROFILE_PARAMETERS else tuple()) for b in bins})
    diagnostic_dir = nominal_dir / "diagnostics" / "period_excursions"
    plot_dir = diagnostic_dir / "plots"
    ensure_directory(diagnostic_dir)
    ensure_directory(plot_dir)

    run_state_payload = {period: {key: value.tolist() for key, value in state.items()} for period, state in run_states.items()}
    dilution_payload = {
        period: {
            str(bin_number): {
                "x_index": record.x_index,
                "t_index": record.t_index,
                "value": record.value,
                "stat_uncertainty": record.stat_uncertainty,
            }
            for (record_period, bin_number), record in dilution_records.items()
            if record_period == period
        }
        for period in PERIODS
    }
    tasks = []
    for bin_number in bins:
        parameters = profile_map.get(bin_number, tuple()) or (None,)
        for profile_parameter in parameters:
            for period in (*PERIODS, COMBINED_PERIOD_KEY):
                tasks.append({"bin_number": bin_number, "period": period,
                    "profile_parameter": profile_parameter, "profile_points": args.diagnostic_profile_points})
            # endfor
        # endfor
    # endfor
    results = []
    mp_context = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=min(workers, len(tasks)), mp_context=mp_context,
        initializer=initialize_fit_worker,
        initargs=(str(cache_path), run_state_payload, dilution_payload),
    ) as executor:
        futures = [executor.submit(_diagnostic_period_fit_worker, task) for task in tasks]
        for index, future in enumerate(as_completed(futures), start=1):
            results.append(future.result())
            print(f"[period diagnostics] completed {index}/{len(tasks)} period/bin tasks", flush=True)
        # endfor
    # endwith

    fit_by_key: dict[tuple[int, str], dict[str, Any]] = {}
    profile_by_key = {(item["bin_number"], item["period"], item.get("profile_parameter")): item for item in results}
    for item in results:
        fit_by_key.setdefault((item["bin_number"], item["period"]), item)
    # endfor
    rows = []
    phi_edges = np.linspace(0.0, 2.0 * math.pi, 13)
    phi_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
    raw_yields: dict[str, Any] = {}
    for bin_number in bins:
        # Raw four-state charge-normalized yields: three periods plus the
        # combined data sample.  The combined panel sums counts and charges
        # before forming each charge-normalized yield.
        fig, axes_grid = plt.subplots(2, 2, figsize=(11.5, 8.2), sharex=True)
        axes = axes_grid.ravel()
        period_payloads = {
            period: _charge_normalized_state_yields(
                events, run_states, bin_number, period, phi_edges
            )
            for period in PERIODS
        }
        combined_payload: dict[str, Any] = {}
        for state_key in period_payloads[PERIODS[0]]:
            charge = sum(period_payloads[p][state_key]["charge"] for p in PERIODS)
            counts = np.sum(
                [np.asarray(period_payloads[p][state_key]["counts"], dtype=float) for p in PERIODS],
                axis=0,
            )
            combined_payload[state_key] = {
                "charge": float(charge),
                "counts": counts.astype(int).tolist(),
                "yield": (counts / charge).tolist() if charge > 0.0 else [math.nan] * len(counts),
                "error": (np.sqrt(counts) / charge).tolist() if charge > 0.0 else [math.nan] * len(counts),
            }
        # endfor
        all_payloads = {**period_payloads, COMBINED_PERIOD_KEY: combined_payload}
        for ax, period in zip(axes, (*PERIODS, COMBINED_PERIOD_KEY)):
            payload = all_payloads[period]
            raw_yields[f"bin_{bin_number:02d}_{period}"] = payload
            state_specs = (
                (-1, -1, "o", r"$P_t<0,h=-1$"), (-1, 1, "s", r"$P_t<0,h=+1$"),
                (1, -1, "^", r"$P_t>0,h=-1$"), (1, 1, "D", r"$P_t>0,h=+1$"),
            )
            for target_sign, h, marker, label in state_specs:
                item = payload[f"t{target_sign:+d}_h{h:+d}"]
                ax.errorbar(phi_centers, item["yield"], yerr=item["error"], fmt=marker, ms=3.5, capsize=2, label=label)
            # endfor
            ax.set_title("Combined" if period == COMBINED_PERIOD_KEY else PERIOD_LABELS[period])
            ax.set_xlabel(r"$\phi$ (rad)")
            ax.grid(alpha=0.25)
        # endfor
        axes[0].set_ylabel("charge-normalized yield")
        axes[2].set_ylabel("charge-normalized yield")
        axes[-1].legend(fontsize=8, frameon=False)
        fig.suptitle(f"Bin {bin_number}: raw target/helicity-state yields")
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(plot_dir / f"bin_{bin_number:02d}_raw_state_yields.png", dpi=180)
        plt.close(fig)

        # Correlation matrices for the three independent period fits and
        # the simultaneous combined fit.  Invalid covariance estimates are
        # never rendered as an identity matrix.
        fig, axes_grid = plt.subplots(2, 2, figsize=(11.5, 9.0))
        axes = axes_grid.ravel()
        image = None
        for ax, period in zip(axes, (*PERIODS, COMBINED_PERIOD_KEY)):
            fit = fit_by_key[(bin_number, period)]["fit"]
            covariance_ok = (
                bool(fit.get("valid", False))
                and bool(fit.get("accurate_covariance", False))
                and bool(fit.get("positive_definite_covariance", False))
                and not bool(fit.get("parameters_at_limit", False))
                and fit.get("correlation") is not None
            )
            if not covariance_ok:
                ax.text(
                    0.5, 0.54, "Covariance unavailable",
                    ha="center", va="center", fontsize=12, fontweight="bold",
                    transform=ax.transAxes,
                )
                ax.text(
                    0.5, 0.44, "fit failed covariance-quality requirement",
                    ha="center", va="center", fontsize=9, transform=ax.transAxes,
                )
                ax.set_axis_off()
                continue
            # endif
            order = fit["parameter_order"]
            physics_indices = [order.index(name) for name in PHYSICS_PARAMETERS]
            corr = np.asarray(fit["correlation"], dtype=float)
            matrix = corr[np.ix_(physics_indices, physics_indices)]
            image = ax.imshow(matrix, vmin=-1.0, vmax=1.0, cmap="coolwarm")
            ax.set_xticks(
                range(len(PHYSICS_PARAMETERS)), PHYSICS_PARAMETERS,
                rotation=45, ha="right",
            )
            ax.set_yticks(range(len(PHYSICS_PARAMETERS)), PHYSICS_PARAMETERS)
            ax.set_title(
                "Combined" if period == COMBINED_PERIOD_KEY else PERIOD_LABELS[period]
            )
        # endfor
        if image is not None:
            # Reserve a dedicated colorbar axis outside the 2x2 grid so it
            # cannot overlap the Sp23/combined matrices.
            fig.subplots_adjust(left=0.07, right=0.88, bottom=0.10, top=0.90, hspace=0.30, wspace=0.28)
            cax = fig.add_axes([0.91, 0.18, 0.018, 0.64])
            fig.colorbar(image, cax=cax, label="correlation")
        # endif
        fig.suptitle(f"Bin {bin_number}: MLE parameter correlations")
        fig.savefig(plot_dir / f"bin_{bin_number:02d}_correlations.png", dpi=180)
        plt.close(fig)

        for period in (*PERIODS, COMBINED_PERIOD_KEY):
            fit = fit_by_key[(bin_number, period)]["fit"]
            if period in PERIODS:
                state = run_states[period]
                qtot = state["q_plus"] + state["q_minus"]
                qsum = float(np.sum(qtot))
                mean_abs_pt = float(np.sum(np.abs(state["pt"]) * qtot) / qsum) if qsum > 0 else math.nan
                cut = cuts[(period, bin_number)]
                beam_pol = BEAM_POLARIZATION[period]
                dilution = dilution_records[(period, bin_number)].value
                dilution_stat = dilution_records[(period, bin_number)].stat_uncertainty
                mx2_low, mx2_high = cut.low_gev2, cut.high_gev2
            else:
                mean_abs_pt = beam_pol = dilution = dilution_stat = mx2_low = mx2_high = math.nan
            # endif
            row = {
                "bin_number": bin_number, "period": period,
                "events": fit["metadata"]["number_of_events"],
                "beam_polarization": beam_pol,
                "charge_weighted_mean_abs_target_polarization": mean_abs_pt,
                "dilution": dilution, "dilution_stat": dilution_stat,
                "mx2_low_gev2": mx2_low, "mx2_high_gev2": mx2_high,
                "fit_valid": fit["valid"], "accurate_covariance": fit["accurate_covariance"],
                "positive_definite_covariance": fit["positive_definite_covariance"],
                "parameters_at_limit": fit["parameters_at_limit"], "edm": fit["edm"],
                "minimum_nll": fit.get("minimum_nll", math.nan),
                "u1_fixed_zero": True, "u2_fixed_zero": True,
            }
            for parameter in PHYSICS_PARAMETERS:
                row[parameter] = fit["values"][parameter]
                row[f"{parameter}_stat"] = fit["errors"][parameter]
                if fit["correlation"] is not None:
                    order = fit["parameter_order"]
                    row[f"corr_ll1_{parameter}"] = fit["correlation"][order.index("ll1")][order.index(parameter)]
                    row[f"corr_ll0_{parameter}"] = fit["correlation"][order.index("ll0")][order.index(parameter)]
                # endif
            # endfor
            rows.append(row)
        # endfor

        for profile_parameter in profile_map.get(bin_number, tuple()):
            fig, ax = plt.subplots(figsize=(7.2, 5.0))
            profiles = []
            for period in (*PERIODS, COMBINED_PERIOD_KEY):
                profile = profile_by_key[(bin_number, period, profile_parameter)]["profile"]
                if profile is None:
                    continue
                # endif
                grid = np.asarray(profile["grid"], dtype=float)
                delta = np.asarray(profile["two_delta_nll"], dtype=float)
                if period in PERIODS:
                    profiles.append((period, grid, delta))
                # endif
                color = "black" if period == COMBINED_PERIOD_KEY else PERIOD_COLORS[period]
                label = "Combined" if period == COMBINED_PERIOD_KEY else PERIOD_LABELS[period]
                ax.plot(grid, delta, marker="o", ms=3, color=color, label=label)
            # endfor
            ax.axhline(1.0, color="black", ls="--", lw=1.0, alpha=0.6)
            ax.axhline(4.0, color="black", ls=":", lw=1.0, alpha=0.6)
            ax.set_xlabel(PARAMETER_LABELS[profile_parameter])
            ax.set_ylabel(r"$2\Delta\mathrm{NLL}$")
            ax.set_ylim(bottom=0.0)
            ax.grid(alpha=0.25)
            ax.legend(frameon=False)
            ax.set_title(f"Bin {bin_number}: period profile likelihoods")
            fig.tight_layout()
            fig.savefig(plot_dir / f"bin_{bin_number:02d}_{profile_parameter}_profile_likelihood.png", dpi=180)
            plt.close(fig)

            # Common-amplitude profile test: interpolate each period profile onto
            # the overlap grid and compare the summed common minimum with the sum
            # of independent profile minima.  Other amplitudes are independently
            # profiled within each period at every fixed value.
            if len(profiles) == 3:
                common_low = max(grid[0] for _, grid, _ in profiles)
                common_high = min(grid[-1] for _, grid, _ in profiles)
                if common_high > common_low:
                    common_grid = np.linspace(common_low, common_high, 401)
                    summed = np.zeros_like(common_grid)
                    for _, grid, delta in profiles:
                        summed += np.interp(common_grid, grid, delta)
                    # endfor
                    common_stat = float(np.min(summed))
                    try:
                        from scipy.stats import chi2 as chi2_distribution
                        common_p = float(chi2_distribution.sf(common_stat, 2))
                    except Exception:
                        common_p = math.nan
                    # endtry
                    for row in rows:
                        if row["bin_number"] == bin_number:
                            row["profile_common_parameter"] = profile_parameter
                            row["profile_common_delta_minus2logL"] = common_stat
                            row["profile_common_pvalue_df2"] = common_p
                        # endif
                    # endfor
                # endif
            # endif
        # endif
    # endfor

    summary = pd.DataFrame(rows).sort_values(["bin_number", "period"])
    summary_path = diagnostic_dir / "period_excursion_diagnostics.csv"
    summary.to_csv(summary_path, index=False)
    write_json(diagnostic_dir / "period_excursion_profiles.json", {
        "bins": list(bins), "profile_parameters": {str(k): list(v) for k, v in profile_map.items()},
        "fit_results": results, "raw_state_yields": raw_yields,
    })
    window_path = run_flagged_exclusivity_window_diagnostic(
        args, root, tuple(bins), profile_map, workers
    )
    print("[period-stability-diagnostics] complete", flush=True)
    print(f"  Window refits: {window_path}", flush=True)
    print(f"  Summary: {summary_path}", flush=True)
    print(f"  Plots:   {plot_dir}", flush=True)
    return 0



# -----------------------------------------------------------------------------
# CLAS6 EG1b exclusive-pi+ cross-check
# -----------------------------------------------------------------------------

CLAS6_W_MIN = 2.0
CLAS6_W_MAX = 2.6
CLAS6_Q2_MIN = 1.0
CLAS6_Q2_MAX = 5.0
PION_MASS_GEV = 0.13957039
NEUTRON_MASS_GEV = 0.93956542
CLAS6_MATCH_VARIABLES = ("Q2", "W", "minus_tprime", "pion_p_lab")


def _clas6_t_and_tprime(w, q2, cos_theta):
    w = np.asarray(w, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    cth = np.asarray(cos_theta, dtype=float)
    q0 = (w*w - PROTON_MASS_GEV**2 - q2) / (2.0*w)
    qmag = np.sqrt(np.maximum(q0*q0 + q2, 0.0))
    epi = (w*w + PION_MASS_GEV**2 - NEUTRON_MASS_GEV**2) / (2.0*w)
    ppi = np.sqrt(np.maximum(epi*epi - PION_MASS_GEV**2, 0.0))
    t = -q2 + PION_MASS_GEV**2 - 2.0*(q0*epi - qmag*ppi*cth)
    tmin = -q2 + PION_MASS_GEV**2 - 2.0*(q0*epi - qmag*ppi)
    return t, tmin - t


def _exclusive_cos_theta_from_t(w, q2, t):
    """Recover cos(theta_pi*) for gamma*p -> pi+n from W,Q2,t."""
    w = np.asarray(w, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    t = np.asarray(t, dtype=float)
    q0 = (w*w - PROTON_MASS_GEV**2 - q2) / (2.0*w)
    qmag = np.sqrt(np.maximum(q0*q0 + q2, 0.0))
    epi = (w*w + PION_MASS_GEV**2 - NEUTRON_MASS_GEV**2) / (2.0*w)
    ppi = np.sqrt(np.maximum(epi*epi - PION_MASS_GEV**2, 0.0))
    denom = 2.0*qmag*ppi
    out = np.full(np.broadcast(w, q2, t).shape, np.nan, dtype=float)
    good = np.isfinite(denom) & (denom > 0.0)
    out[good] = (
        t[good] + q2[good] - PION_MASS_GEV**2 + 2.0*q0[good]*epi[good]
    ) / denom[good]
    out[good] = np.clip(out[good], -1.0, 1.0)
    return out


def _exclusive_pion_lab_momentum(beam_energy, w, q2, cos_theta):
    """Two-body pi+ lab momentum for a stationary proton target.

    The boost is along q.  This is exact for gamma*p -> pi+n at the supplied
    W,Q2,cos(theta*) and is used only as a phase-space/acceptance diagnostic.
    """
    e0 = np.asarray(beam_energy, dtype=float)
    w = np.asarray(w, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    cth = np.asarray(cos_theta, dtype=float)
    nu = (w*w + q2 - PROTON_MASS_GEV**2) / (2.0*PROTON_MASS_GEV)
    qmag_lab = np.sqrt(np.maximum(nu*nu + q2, 0.0))
    total_e = PROTON_MASS_GEV + nu
    beta = np.divide(qmag_lab, total_e, out=np.zeros_like(qmag_lab), where=total_e > 0.0)
    gamma = np.divide(total_e, w, out=np.full_like(total_e, np.nan), where=w > 0.0)
    epi_star = (w*w + PION_MASS_GEV**2 - NEUTRON_MASS_GEV**2) / (2.0*w)
    ppi_star = np.sqrt(np.maximum(epi_star*epi_star - PION_MASS_GEV**2, 0.0))
    epi_lab = gamma*(epi_star + beta*ppi_star*cth)
    ppi_lab = np.sqrt(np.maximum(epi_lab*epi_lab - PION_MASS_GEV**2, 0.0))
    # e0 is deliberately accepted/checked here: the same W,Q2 point can only
    # occur at a beam energy for which the electron kinematics are physical.
    physical = np.isfinite(e0) & (e0 > nu)
    return np.where(physical, ppi_lab, np.nan)


def _clas6_analysis_bin(xb, minus_tprime):
    xb_edges = np.asarray([0.10, 0.25, 0.35, 0.45, 0.60], dtype=float)
    tp_edges = np.asarray([0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25], dtype=float)
    ix = np.searchsorted(xb_edges, xb, side="right") - 1
    it = np.searchsorted(tp_edges, minus_tprime, side="right") - 1
    good = (ix >= 0) & (ix < 4) & (it >= 0) & (it < 6)
    out = np.full(np.asarray(xb).shape, -1, dtype=int)
    out[good] = ix[good]*6 + it[good] + 1
    return out


def _load_clas6_exclpip(path):
    names = [
        "E", "Code", "W_bin", "Q2_bin", "theta_bin", "phi_bin",
        "W", "Q2", "cos_theta", "phi", "epsilon",
        "sin_phi", "sin_2phi", "cos_phi", "cos_2phi",
        "A_LL", "A_LL_stat", "A_UL", "A_UL_stat",
    ]
    frame = pd.read_csv(
        path, sep=r"\s+", comment="#", skiprows=5, names=names,
        engine="python", on_bad_lines="skip",
    )
    for column in names:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    # endfor
    frame = frame.dropna(subset=["E", "Code", "W", "Q2", "cos_theta", "phi"])
    frame = frame[
        (frame["Code"] == 1)
        & (frame["W"] > CLAS6_W_MIN) & (frame["W"] < CLAS6_W_MAX)
        & (frame["Q2"] > CLAS6_Q2_MIN) & (frame["Q2"] < CLAS6_Q2_MAX)
    ].copy()
    frame["xB"] = frame["Q2"] / (frame["W"]**2 - PROTON_MASS_GEV**2 + frame["Q2"])
    t, tp = _clas6_t_and_tprime(frame["W"], frame["Q2"], frame["cos_theta"])
    frame["t"] = t
    frame["minus_tprime"] = tp
    frame["pion_p_lab"] = _exclusive_pion_lab_momentum(
        frame["E"].to_numpy(float), frame["W"].to_numpy(float),
        frame["Q2"].to_numpy(float), frame["cos_theta"].to_numpy(float),
    )
    frame["analysis_bin"] = _clas6_analysis_bin(frame["xB"].to_numpy(), tp)
    frame = frame[frame["analysis_bin"] > 0].copy()
    return frame


def _weighted_linear_fit(design, values, errors):
    x = np.asarray(design, dtype=float)
    y = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    good = np.all(np.isfinite(x), axis=1) & np.isfinite(y) & np.isfinite(e) & (e > 0.0)
    x, y, e = x[good], y[good], e[good]
    if len(y) < x.shape[1] + 1:
        return None
    # endif
    weight = 1.0/(e*e)
    normal = x.T @ (weight[:, None]*x)
    rhs = x.T @ (weight*y)
    try:
        cov = np.linalg.inv(normal)
        beta = cov @ rhs
    except np.linalg.LinAlgError:
        return None
    # endtry
    residual = y - x@beta
    chi2 = float(np.sum((residual/e)**2))
    return beta, np.sqrt(np.diag(cov)), cov, chi2, int(len(y)-len(beta)), int(len(y))


def _fit_clas6_bin(data, reflect_phi=True):
    """Fit published CLAS6 beam-axis asymmetries using each point's epsilon.

    The adopted comparison convention is phi_RGC = 2*pi - phi_CLAS6.  Thus the
    CLAS6 sine harmonics are reflected while constant/cosine harmonics are not.
    Each CLAS6 point retains its own published epsilon, so its beam-energy
    dependence enters the depolarization factors point by point.
    """
    eps = data["epsilon"].to_numpy(float)
    phi = data["phi"].to_numpy(float)
    if reflect_phi:
        phi = np.mod(2.0*math.pi - phi, 2.0*math.pi)
    # endif
    r_b = eps
    r_c = np.sqrt(np.maximum(1.0 - eps*eps, 0.0))
    r_v = np.sqrt(np.maximum(2.0*eps*(1.0 + eps), 0.0))
    r_w = np.sqrt(np.maximum(2.0*eps*(1.0 - eps), 0.0))
    ul_design = np.column_stack((r_v*np.sin(phi), r_b*np.sin(2.0*phi)))
    ll_design = np.column_stack((r_c, r_w*np.cos(phi)))
    ul = _weighted_linear_fit(ul_design, data["A_UL"], data["A_UL_stat"])
    ll = _weighted_linear_fit(ll_design, data["A_LL"], data["A_LL_stat"])
    if ul is None or ll is None:
        return None
    # endif
    return {
        "ul1": float(ul[0][0]), "ul1_stat": float(ul[1][0]),
        "ul2": float(ul[0][1]), "ul2_stat": float(ul[1][1]),
        "ll0": float(ll[0][0]), "ll0_stat": float(ll[1][0]),
        "ll1": float(ll[0][1]), "ll1_stat": float(ll[1][1]),
        "ul_chi2": ul[3], "ul_ndf": ul[4],
        "ll_chi2": ll[3], "ll_ndf": ll[4],
        "n_ul_points": ul[5], "n_ll_points": ll[5],
    }


def _add_rgc_crosscheck_kinematics(events):
    """Add exclusive reconstructed pion momentum to an RGC event dictionary."""
    out = {key: np.asarray(value) for key, value in events.items()}
    period_index = np.asarray(out["period_index"], dtype=int)
    beam = np.asarray([BEAM_ENERGY_GEV[PERIODS[i]] for i in period_index], dtype=float)
    cos_theta = _exclusive_cos_theta_from_t(
        out["W"], out["Q2"], -np.asarray(out["minus_t"], dtype=float),
    )
    out["cos_theta_star"] = cos_theta
    out["pion_p_lab"] = _exclusive_pion_lab_momentum(
        beam, out["W"], out["Q2"], cos_theta,
    )
    return out


def _common_envelope_for_bin(rgc_events, clas6_bin, bin_number):
    """Intersection of the populated RGC/CLAS6 ranges in four key variables."""
    rgc_mask = np.asarray(rgc_events["bin_number"]) == int(bin_number)
    envelope = {}
    for variable in CLAS6_MATCH_VARIABLES:
        r = np.asarray(rgc_events[variable], dtype=float)[rgc_mask]
        c = pd.to_numeric(clas6_bin[variable], errors="coerce").to_numpy(float)
        r = r[np.isfinite(r)]
        c = c[np.isfinite(c)]
        if len(r) == 0 or len(c) == 0:
            return None
        # endif
        low = max(float(np.min(r)), float(np.min(c)))
        high = min(float(np.max(r)), float(np.max(c)))
        if not np.isfinite(low) or not np.isfinite(high) or high <= low:
            return None
        # endif
        envelope[variable] = (low, high)
    # endfor
    return envelope


def _apply_envelope_to_rgc(events, bin_number, envelope):
    mask = np.asarray(events["bin_number"]) == int(bin_number)
    for variable, (low, high) in envelope.items():
        values = np.asarray(events[variable], dtype=float)
        mask &= np.isfinite(values) & (values >= low) & (values <= high)
    # endfor
    return mask


def _apply_envelope_to_clas6(frame, envelope):
    mask = np.ones(len(frame), dtype=bool)
    for variable, (low, high) in envelope.items():
        values = pd.to_numeric(frame[variable], errors="coerce").to_numpy(float)
        mask &= np.isfinite(values) & (values >= low) & (values <= high)
    # endfor
    return frame.loc[mask].copy()


def _kinematic_summary(prefix, rgc_events, rgc_mask, clas6_frame):
    result = {}
    variables = ("xB", "Q2", "W", "minus_tprime", "epsilon", "pion_p_lab")
    for variable in variables:
        r = np.asarray(rgc_events[variable], dtype=float)[rgc_mask]
        c = pd.to_numeric(clas6_frame[variable], errors="coerce").to_numpy(float)
        for dataset, values in (("rgc", r), ("clas6", c)):
            values = values[np.isfinite(values)]
            if len(values) == 0:
                continue
            # endif
            result[f"{prefix}_{variable}_{dataset}_mean"] = float(np.mean(values))
            result[f"{prefix}_{variable}_{dataset}_rms"] = float(np.std(values))
            result[f"{prefix}_{variable}_{dataset}_q10"] = float(np.quantile(values, 0.10))
            result[f"{prefix}_{variable}_{dataset}_median"] = float(np.median(values))
            result[f"{prefix}_{variable}_{dataset}_q90"] = float(np.quantile(values, 0.90))
        # endfor
    # endfor
    return result


_CLAS6_MATCH_ENVELOPES = None


def initialize_clas6_cross_check_worker(
    cache_path_text: str,
    run_state_payload: dict[str, dict[str, list[float] | list[int]]],
    dilution_payload: dict[str, dict[str, dict[str, float | int]]],
    match_envelopes: dict[int, dict[str, tuple[float, float]]] | None = None,
) -> None:
    """Load cache once per worker; optionally retain only common-envelope events."""
    initialize_fit_worker(cache_path_text, run_state_payload, dilution_payload)
    global _WORKER_EVENTS, _CLAS6_MATCH_ENVELOPES
    common = (
        (_WORKER_EVENTS["W"] > CLAS6_W_MIN)
        & (_WORKER_EVENTS["W"] < CLAS6_W_MAX)
        & (_WORKER_EVENTS["Q2"] > CLAS6_Q2_MIN)
        & (_WORKER_EVENTS["Q2"] < CLAS6_Q2_MAX)
    )
    _WORKER_EVENTS = {key: np.asarray(value)[common] for key, value in _WORKER_EVENTS.items()}
    _WORKER_EVENTS = _add_rgc_crosscheck_kinematics(_WORKER_EVENTS)
    _CLAS6_MATCH_ENVELOPES = match_envelopes
    if match_envelopes:
        keep = np.zeros(len(_WORKER_EVENTS["bin_number"]), dtype=bool)
        for bin_number, envelope in match_envelopes.items():
            keep |= _apply_envelope_to_rgc(_WORKER_EVENTS, int(bin_number), envelope)
        # endfor
        _WORKER_EVENTS = {key: np.asarray(value)[keep] for key, value in _WORKER_EVENTS.items()}
    # endif


def _clas6_rgc_fit_worker(bin_number: int) -> tuple[int, dict[str, Any]]:
    if _WORKER_EVENTS is None or _WORKER_RUN_STATES is None or _WORKER_DILUTION_RECORDS is None:
        raise RuntimeError("CLAS6 cross-check worker was not initialized.")
    # endif
    result = fit_one_variant(
        _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
        int(bin_number), "nominal",
    )
    return int(bin_number), result


def _run_clas6_rgc_parallel_fit(
    cache_path, run_state_payload, dilution_payload, eligible_bins, workers,
    match_envelopes=None, label="broad",
):
    fits = {}
    n_workers = max(1, min(int(workers), 8, len(eligible_bins)))
    print(
        f"[CLAS6 cross-check] {label}: fitting {len(eligible_bins)} bins with "
        f"{n_workers} workers", flush=True,
    )
    mp_context = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=n_workers, mp_context=mp_context,
        initializer=initialize_clas6_cross_check_worker,
        initargs=(str(cache_path), run_state_payload, dilution_payload, match_envelopes),
    ) as executor:
        futures = {executor.submit(_clas6_rgc_fit_worker, b): b for b in eligible_bins}
        for completed, future in enumerate(as_completed(futures), start=1):
            bin_number, fit = future.result()
            fits[bin_number] = fit
            print(
                f"[CLAS6 cross-check] {label}: completed bin {bin_number:02d} "
                f"({completed}/{len(eligible_bins)})", flush=True,
            )
        # endfor
    # endwith
    return fits


def _comparison_summary(frame, fit_tag):
    rows = []
    for parameter in ("ul1", "ul2", "ll0", "ll1"):
        column = f"{fit_tag}_{parameter}_pull"
        pulls = pd.to_numeric(frame.get(column, pd.Series(dtype=float)), errors="coerce").to_numpy(float)
        pulls = pulls[np.isfinite(pulls)]
        rows.append({
            "comparison": fit_tag, "parameter": parameter, "n": int(len(pulls)),
            "chi2": float(np.sum(pulls*pulls)), "ndf": int(len(pulls)),
            "chi2_per_ndf": float(np.mean(pulls*pulls)) if len(pulls) else np.nan,
            "pull_mean": float(np.mean(pulls)) if len(pulls) else np.nan,
            "pull_rms": float(np.sqrt(np.mean(pulls*pulls))) if len(pulls) else np.nan,
        })
    # endfor
    return rows


def run_clas6_cross_check(args, workers):
    """CLAS6 EG1b comparison with broad and common-envelope phase-space fits."""
    out = args.output_dir.expanduser().resolve() / "clas6_cross_check"
    tables, plots = out / "tables", out / "plots"
    for directory in (out, tables, plots):
        ensure_directory(directory)
    # endfor

    clas6_path = args.clas6_data.expanduser().resolve()
    if not clas6_path.is_file():
        raise FileNotFoundError(f"Missing CLAS6 data file: {clas6_path}")
    # endif
    clas6 = _load_clas6_exclpip(clas6_path)
    clas6.to_csv(tables / "clas6_points_common_W_Q2.csv", index=False)

    phi_c6 = clas6["phi"].to_numpy(float)
    trig_checks = {
        "sin_phi": np.sin(phi_c6), "sin_2phi": np.sin(2.0*phi_c6),
        "cos_phi": np.cos(phi_c6), "cos_2phi": np.cos(2.0*phi_c6),
    }
    print("\n[CLAS6 cross-check] exclpip.txt trigonometric-column check", flush=True)
    for column, calculated in trig_checks.items():
        tabulated = clas6[column].to_numpy(float)
        good = np.isfinite(tabulated) & np.isfinite(calculated)
        max_abs = float(np.max(np.abs(tabulated[good] - calculated[good]))) if np.any(good) else np.nan
        print(f"  {column:8s}: max |table - calculated| = {max_abs:.3e}", flush=True)
    # endfor

    cache_path = _rga_variant_cache_paths(args)["nominal"]
    if not cache_path.is_file():
        raise FileNotFoundError(f"Missing nominal selected-event cache: {cache_path}")
    # endif
    events = load_event_cache(cache_path)
    common = (
        (events["W"] > CLAS6_W_MIN) & (events["W"] < CLAS6_W_MAX)
        & (events["Q2"] > CLAS6_Q2_MIN) & (events["Q2"] < CLAS6_Q2_MAX)
    )
    cut_events = {key: np.asarray(value)[common] for key, value in events.items()}
    cut_events = _add_rgc_crosscheck_kinematics(cut_events)

    run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    run_states = run_state_arrays(run_records)
    dilution_path = (
        args.dilution_json.expanduser().resolve() if args.dilution_json
        else find_default_dilution_json(args.dilution_dir.expanduser().resolve()).resolve()
    )
    dilution_records = load_dilution_factors(dilution_path, "nominal")

    eligible_bins = []
    envelopes = {}
    matched_clas6 = {}
    for bin_number in range(1, NUMBER_OF_BINS + 1):
        rgc_bin_mask = np.asarray(cut_events["bin_number"]) == bin_number
        n_rgc = int(np.count_nonzero(rgc_bin_mask))
        c6 = clas6[clas6["analysis_bin"] == bin_number]
        if n_rgc == 0 or c6.empty:
            continue
        # endif
        envelope = _common_envelope_for_bin(cut_events, c6, bin_number)
        c6_match = _apply_envelope_to_clas6(c6, envelope) if envelope is not None else c6.iloc[0:0].copy()
        rgc_match = _apply_envelope_to_rgc(cut_events, bin_number, envelope) if envelope is not None else np.zeros(len(cut_events["bin_number"]), dtype=bool)
        n_rgc_match = int(np.count_nonzero(rgc_match))
        if envelope is not None and n_rgc_match > 0 and len(c6_match) >= 5:
            envelopes[bin_number] = envelope
            matched_clas6[bin_number] = c6_match
        # endif
        eligible_bins.append(bin_number)
        print(
            f"[CLAS6 cross-check] bin {bin_number:02d}: broad RGC={n_rgc:,}, "
            f"CLAS6={len(c6):,}; matched RGC={n_rgc_match:,}, CLAS6={len(c6_match):,}",
            flush=True,
        )
    # endfor

    run_state_payload = {
        period: {key: np.asarray(values).tolist() for key, values in state.items()}
        for period, state in run_states.items()
    }
    dilution_payload = {
        period: {
            str(bin_number): {
                "x_index": record.x_index, "t_index": record.t_index,
                "value": record.value, "stat_uncertainty": record.stat_uncertainty,
            }
            for (record_period, bin_number), record in dilution_records.items()
            if record_period == period
        }
        for period in PERIODS
    }

    broad_fits = _run_clas6_rgc_parallel_fit(
        cache_path, run_state_payload, dilution_payload, eligible_bins, workers,
        match_envelopes=None, label="broad overlap",
    )
    matched_bins = sorted(envelopes)
    matched_fits = _run_clas6_rgc_parallel_fit(
        cache_path, run_state_payload, dilution_payload, matched_bins, workers,
        match_envelopes=envelopes, label="common envelope",
    ) if matched_bins else {}

    envelope_rows = []
    rows = []
    for bin_number in eligible_bins:
        c6_broad = clas6[clas6["analysis_bin"] == bin_number].copy()
        broad_mask = np.asarray(cut_events["bin_number"]) == bin_number
        broad_c6_fit = _fit_clas6_bin(c6_broad, reflect_phi=True)
        if broad_c6_fit is None or bin_number not in broad_fits:
            continue
        # endif
        row = {
            "bin_number": bin_number,
            "n_rgc_broad": int(np.count_nonzero(broad_mask)),
            "n_clas6_broad": int(len(c6_broad)),
        }
        row.update(_kinematic_summary("broad", cut_events, broad_mask, c6_broad))

        if bin_number in envelopes:
            envelope = envelopes[bin_number]
            rgc_match_mask = _apply_envelope_to_rgc(cut_events, bin_number, envelope)
            c6_match = matched_clas6[bin_number]
            row["n_rgc_matched"] = int(np.count_nonzero(rgc_match_mask))
            row["n_clas6_matched"] = int(len(c6_match))
            row.update(_kinematic_summary("matched", cut_events, rgc_match_mask, c6_match))
            envelope_row = {"bin_number": bin_number}
            for variable, (low, high) in envelope.items():
                envelope_row[f"{variable}_min"] = low
                envelope_row[f"{variable}_max"] = high
            # endfor
            envelope_rows.append(envelope_row)
        else:
            row["n_rgc_matched"] = 0
            row["n_clas6_matched"] = 0
        # endif

        for fit_tag, rgc_fit, c6_data in (
            ("broad", broad_fits.get(bin_number), c6_broad),
            ("matched", matched_fits.get(bin_number), matched_clas6.get(bin_number)),
        ):
            if rgc_fit is None or c6_data is None or len(c6_data) < 5:
                continue
            # endif
            c6_fit = _fit_clas6_bin(c6_data, reflect_phi=True)
            if c6_fit is None:
                continue
            # endif
            row[f"{fit_tag}_rgc_fit_valid"] = bool(rgc_fit["valid"])
            row[f"{fit_tag}_rgc_edm"] = float(rgc_fit["edm"])
            for parameter in ("ul1", "ul2", "ll0", "ll1"):
                rgc_value = float(rgc_fit["values"][parameter])
                rgc_stat = float(rgc_fit["errors"][parameter])
                c6_value = float(c6_fit[parameter])
                c6_stat = float(c6_fit[f"{parameter}_stat"])
                sigma = math.hypot(rgc_stat, c6_stat)
                row[f"{fit_tag}_{parameter}_rgc"] = rgc_value
                row[f"{fit_tag}_{parameter}_rgc_stat"] = rgc_stat
                row[f"{fit_tag}_{parameter}_clas6"] = c6_value
                row[f"{fit_tag}_{parameter}_clas6_stat"] = c6_stat
                row[f"{fit_tag}_{parameter}_delta"] = rgc_value - c6_value
                row[f"{fit_tag}_{parameter}_pull"] = (rgc_value - c6_value)/sigma if sigma > 0.0 else np.nan
            # endfor
        # endfor
        rows.append(row)
    # endfor

    frame = pd.DataFrame(rows)
    frame.to_csv(tables / "clas6_rgc_phase_space_matched_comparison.csv", index=False)
    pd.DataFrame(envelope_rows).to_csv(tables / "clas6_rgc_common_envelopes.csv", index=False)

    summary_rows = _comparison_summary(frame, "broad") + _comparison_summary(frame, "matched")
    summary_frame = pd.DataFrame(summary_rows)
    summary_frame.to_csv(tables / "clas6_rgc_statistical_summary.csv", index=False)
    print("\n[CLAS6 cross-check] statistical-only summary (adopted phi reflection)", flush=True)
    print(summary_frame.to_string(index=False), flush=True)

    # Compact per-bin phase-space table for immediate inspection.
    phase_columns = [
        "bin_number", "n_rgc_broad", "n_clas6_broad", "n_rgc_matched", "n_clas6_matched",
        "broad_Q2_rgc_mean", "broad_Q2_clas6_mean", "matched_Q2_rgc_mean", "matched_Q2_clas6_mean",
        "broad_W_rgc_mean", "broad_W_clas6_mean", "matched_W_rgc_mean", "matched_W_clas6_mean",
        "broad_minus_tprime_rgc_mean", "broad_minus_tprime_clas6_mean",
        "matched_minus_tprime_rgc_mean", "matched_minus_tprime_clas6_mean",
        "broad_epsilon_rgc_mean", "broad_epsilon_clas6_mean",
        "matched_epsilon_rgc_mean", "matched_epsilon_clas6_mean",
        "broad_pion_p_lab_rgc_mean", "broad_pion_p_lab_clas6_mean",
        "matched_pion_p_lab_rgc_mean", "matched_pion_p_lab_clas6_mean",
    ]
    existing_phase_columns = [column for column in phase_columns if column in frame.columns]
    frame[existing_phase_columns].to_csv(tables / "clas6_rgc_kinematic_diagnostics.csv", index=False)

    if not args.skip_plots and not frame.empty:
        labels = {
            "ul1": r"$A_{UL,\mathrm{lab}}^{\sin\phi}$",
            "ul2": r"$A_{UL,\mathrm{lab}}^{\sin2\phi}$",
            "ll0": r"$A_{LL,\mathrm{lab}}$",
            "ll1": r"$A_{LL,\mathrm{lab}}^{\cos\phi}$",
        }
        for parameter, ylabel in labels.items():
            valid = frame[np.isfinite(pd.to_numeric(frame.get(f"broad_{parameter}_pull"), errors="coerce"))].copy()
            if valid.empty:
                continue
            # endif
            x = np.arange(len(valid))
            fig, ax = plt.subplots(figsize=(10.5, 5.5))
            # Keep the publication-level comparison deliberately simple.  The
            # common-envelope study is retained in the diagnostic tables, but the
            # primary figures show only the broad-overlap CLAS6 and RGC results.
            ax.errorbar(
                x - 0.08, valid[f"broad_{parameter}_clas6"],
                yerr=valid[f"broad_{parameter}_clas6_stat"], fmt="o", capsize=3,
                color="red", label=r"CLAS6 ($2\pi-\phi$)",
            )
            ax.errorbar(
                x + 0.08, valid[f"broad_{parameter}_rgc"],
                yerr=valid[f"broad_{parameter}_rgc_stat"], fmt="s", capsize=3,
                color="blue", label="CLAS12 RGC",
            )
            ax.axhline(0.0, lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(valid["bin_number"].astype(int))
            ax.set_xlabel("RGC analysis bin")
            ax.set_ylabel(ylabel)
            # Use a common fixed asymmetry range for every CLAS6--RGC amplitude
            # comparison so visual differences are directly comparable panel to panel.
            ax.set_ylim(-1.0, 1.0)
            ax.legend()
            fig.tight_layout()
            fig.savefig(plots / f"clas6_rgc_{parameter}.png", dpi=200)
            plt.close(fig)
        # endfor

        # Mean-kinematics diagnostics: broad and matched means for both experiments.
        kin_labels = {
            "Q2": r"$Q^2$ (GeV$^2$)", "W": r"$W$ (GeV)",
            "minus_tprime": r"$-t'$ (GeV$^2$)", "epsilon": r"$\epsilon$",
            "pion_p_lab": r"$p_{\pi}$ (GeV)",
        }
        for variable, ylabel in kin_labels.items():
            if f"broad_{variable}_rgc_mean" not in frame.columns:
                continue
            # endif
            x = np.arange(len(frame))
            fig, ax = plt.subplots(figsize=(10.5, 5.5))
            ax.plot(x, frame[f"broad_{variable}_clas6_mean"], "o", label="CLAS6 broad")
            ax.plot(x, frame[f"broad_{variable}_rgc_mean"], "s", label="RGC broad")
            if f"matched_{variable}_clas6_mean" in frame.columns:
                ax.plot(x, frame[f"matched_{variable}_clas6_mean"], "^", label="CLAS6 matched")
                ax.plot(x, frame[f"matched_{variable}_rgc_mean"], "D", label="RGC matched")
            # endif
            ax.set_xticks(x)
            ax.set_xticklabels(frame["bin_number"].astype(int))
            ax.set_xlabel("RGC analysis bin")
            ax.set_ylabel(ylabel)
            ax.legend(ncol=2)
            fig.tight_layout()
            fig.savefig(plots / f"clas6_rgc_kinematics_{variable}.png", dpi=200)
            plt.close(fig)
        # endfor
    # endif

    write_json(out / "clas6_cross_check_manifest.json", {
        "selection": {
            "broad": {"W_GeV": [CLAS6_W_MIN, CLAS6_W_MAX], "Q2_GeV2": [CLAS6_Q2_MIN, CLAS6_Q2_MAX], "channel_code": 1},
            "matched": (
                "Within each RGC (xB,-tprime) bin, both samples are additionally restricted "
                "to the intersection of their populated Q2, W, -tprime, and reconstructed "
                "exclusive pion lab-momentum ranges. No arbitrary detector-momentum cut is imposed."
            ),
        },
        "phi_convention": (
            "Adopt phi_RGC = 2*pi - phi_CLAS6 for this cross-check. This reverses the sine "
            "harmonics and leaves constant/cosine harmonics unchanged. The earlier diagnostic "
            "showed that this transformation resolves the sine-only sign disagreement."
        ),
        "depolarization_factors": (
            "CLAS6 uses its tabulated epsilon point by point: B/A=epsilon, C/A=sqrt(1-epsilon^2), "
            "V/A=sqrt(2 epsilon (1+epsilon)), W/A=sqrt(2 epsilon (1-epsilon)). RGC continues to "
            "use its own event-by-event DepB/DepA, DepC/DepA, DepV/DepA, and DepW/DepA values. "
            "Thus the different beam energies are not treated with a common depolarization factor."
        ),
        "pion_momentum": (
            "For both experiments p_pi(lab) is reconstructed from exclusive two-body gamma*p -> pi+n "
            "kinematics. CLAS6 uses tabulated E,W,Q2,cos(theta*); RGC obtains cos(theta*) from W,Q2,t "
            "and uses the period beam energy. It is a phase-space diagnostic/matching variable, not a "
            "new production selection."
        ),
        "uncertainties": "Statistical only; no RGC or CLAS6 systematic uncertainties are reproduced.",
        "binning": "Common RGC 24-bin (xB,-tprime) scheme.",
        "clas6_data": str(clas6_path),
        "rgc_cache": str(cache_path),
    })
    print(f"[CLAS6 cross-check] wrote upgraded study to {out}", flush=True)
    return 0



# -----------------------------------------------------------------------------
# Diagnostic: unbinned MLE versus 12-bin-in-phi chi2 extraction
# -----------------------------------------------------------------------------

MLE_CHI2_PARAMETERS: tuple[str, ...] = ("lu1", "ul1", "ul2", "ll0", "ll1")
MLE_CHI2_PHI_BINS = 12


def _binned_phi_chi2_fit(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    bin_number: int,
) -> dict[str, Any]:
    """Fit the five polarized amplitudes to 12 charge-normalized phi bins.

    The unpolarized cos(phi) and cos(2phi) amplitudes are fixed to zero.  For
    each run period and phi bin a free spin-independent normalization is
    profiled analytically, so acceptance and the unpolarized phi shape are not
    fitted as physics amplitudes.  The period/bin dilution factors are Gaussian
    constrained nuisances, exactly paralleling their statistical treatment in
    the nominal MLE.
    """
    phi_edges = np.linspace(0.0, 2.0 * math.pi, MLE_CHI2_PHI_BINS + 1)
    rows: list[dict[str, float]] = []

    for period in PERIODS:
        mask = (
            (events["bin_number"] == bin_number)
            & (events["period_index"] == PERIOD_INDEX[period])
        )
        phi = events["phi"][mask].astype(float, copy=False)
        runnum = events["runnum"][mask].astype(np.int32, copy=False)
        helicity = events["helicity"][mask].astype(np.int8, copy=False)
        state = run_states[period]
        lookup = {int(run): i for i, run in enumerate(state["run"])}
        event_state = np.asarray([lookup[int(run)] for run in runnum], dtype=np.int32)
        target_sign = np.sign(state["pt"][event_state]).astype(np.int8)

        dep_arrays = {
            "r_b": events["rB"][mask].astype(float, copy=False),
            "r_c": events["rC"][mask].astype(float, copy=False),
            "r_v": events["rV"][mask].astype(float, copy=False),
            "r_w": events["rW"][mask].astype(float, copy=False),
        }
        phi_index = np.searchsorted(phi_edges, phi, side="right") - 1
        phi_index = np.clip(phi_index, 0, MLE_CHI2_PHI_BINS - 1)

        for j in range(MLE_CHI2_PHI_BINS):
            in_phi = phi_index == j
            if not np.any(in_phi):
                continue
            # endif
            phi_mean = float(np.angle(np.mean(np.exp(1j * phi[in_phi]))) % (2.0 * math.pi))
            dep_mean = {name: float(np.mean(values[in_phi])) for name, values in dep_arrays.items()}

            for sgn in (-1, 1):
                run_mask = np.sign(state["pt"]) == sgn
                for h in (-1, 1):
                    charges = state["q_plus"] if h > 0 else state["q_minus"]
                    charge = float(np.sum(charges[run_mask]))
                    selected = in_phi & (target_sign == sgn) & (helicity == h)
                    count = int(np.count_nonzero(selected))
                    if charge <= 0.0 or count <= 0:
                        continue
                    # endif
                    # Charge-weighted polarization is the appropriate effective
                    # polarization for this charge-normalized spin-state yield.
                    qsel = charges[run_mask]
                    ptsel = state["pt"][run_mask]
                    pt_eff = float(np.sum(qsel * ptsel) / np.sum(qsel))
                    rows.append({
                        "period": period,
                        "phi_bin": float(j),
                        "phi": phi_mean,
                        "h": float(h),
                        "pt": pt_eff,
                        "yield": float(count / charge),
                        "error": float(math.sqrt(count) / charge),
                        **dep_mean,
                    })
                # endfor
            # endfor
        # endfor
    # endfor

    if not rows:
        raise RuntimeError(f"No binned chi2 data for analysis bin {bin_number}.")
    # endif

    frame = pd.DataFrame(rows)
    nuisance_names = [f"f_{period}" for period in PERIODS]
    names = list(MLE_CHI2_PARAMETERS) + nuisance_names
    starts = [0.0] * len(MLE_CHI2_PARAMETERS) + [
        dilution_records[(period, bin_number)].value for period in PERIODS
    ]

    def objective(*values: float) -> float:
        pars = dict(zip(names, values))
        total = 0.0
        # Profile one arbitrary spin-independent normalization independently in
        # every period/phi bin.  This makes the comparison insensitive to the
        # unpolarized phi distribution while retaining the polarized dependence.
        for (period, phi_bin), group in frame.groupby(["period", "phi_bin"], sort=False):
            f = float(pars[f"f_{period}"])
            phi_g = group["phi"].to_numpy(float)
            h_g = group["h"].to_numpy(float)
            pt_g = group["pt"].to_numpy(float)
            factor = evaluate_cross_section_factor(
                variant="nominal",
                phi=phi_g,
                r_b=group["r_b"].to_numpy(float),
                r_c=group["r_c"].to_numpy(float),
                r_v=group["r_v"].to_numpy(float),
                r_w=group["r_w"].to_numpy(float),
                sin_theta_gamma=np.zeros(len(group), dtype=float),
                cos_theta_gamma=np.ones(len(group), dtype=float),
                helicity=h_g,
                beam_polarization=BEAM_POLARIZATION[period],
                target_polarization=pt_g,
                dilution=f,
                transverse_scales=None,
                u1=0.0, u2=0.0,
                lu1=float(pars["lu1"]), ul1=float(pars["ul1"]),
                ul2=float(pars["ul2"]), ll0=float(pars["ll0"]),
                ll1=float(pars["ll1"]),
            )
            if np.any(~np.isfinite(factor)) or np.any(factor <= 0.0):
                return 1.0e100
            # endif
            y = group["yield"].to_numpy(float)
            err = group["error"].to_numpy(float)
            weight = 1.0 / (err * err)
            denominator = float(np.sum(weight * factor * factor))
            if denominator <= 0.0:
                return 1.0e100
            # endif
            norm = float(np.sum(weight * factor * y) / denominator)
            residual = (y - norm * factor) / err
            total += float(np.sum(residual * residual))
        # endfor

        # Statistical dilution-factor uncertainty: each period/bin dilution is
        # a Gaussian-constrained nuisance, just as in the unbinned likelihood.
        for period in PERIODS:
            record = dilution_records[(period, bin_number)]
            f = float(pars[f"f_{period}"])
            if f <= 0.0:
                return 1.0e100
            # endif
            if record.stat_uncertainty > 0.0:
                total += ((f - record.value) / record.stat_uncertainty) ** 2
            # endif
        # endfor
        return total

    m = Minuit(objective, *starts, name=names)
    m.errordef = Minuit.LEAST_SQUARES
    m.print_level = 0
    for name in MLE_CHI2_PARAMETERS:
        m.limits[name] = PARAMETER_LIMITS[name]
    # endfor
    for period in PERIODS:
        name = f"f_{period}"
        record = dilution_records[(period, bin_number)]
        width = max(8.0 * record.stat_uncertainty, 0.20 * record.value, 0.02)
        m.limits[name] = (max(1.0e-6, record.value - width), record.value + width)
        if record.stat_uncertainty == 0.0:
            m.fixed[name] = True
        # endif
    # endfor
    m.migrad()
    m.hesse()

    # There are four spin-state yields per period/phi bin and one profiled
    # normalization for each such group. Gaussian nuisance constraints add one
    # datum and one nuisance parameter apiece, so they cancel in the nominal ndf.
    groups = int(frame.groupby(["period", "phi_bin"]).ngroups)
    n_data = int(len(frame))
    ndf = max(0, n_data - groups - len(MLE_CHI2_PARAMETERS))
    return {
        "bin_number": int(bin_number),
        "values": {name: float(m.values[name]) for name in MLE_CHI2_PARAMETERS},
        "errors": {name: float(m.errors[name]) for name in MLE_CHI2_PARAMETERS},
        "dilution_values": {period: float(m.values[f"f_{period}"]) for period in PERIODS},
        "chi2": float(m.fval),
        "ndf": int(ndf),
        "chi2_per_ndf": float(m.fval / ndf) if ndf > 0 else math.nan,
        "valid": bool(m.fmin.is_valid),
        "number_of_state_yields": n_data,
        "number_of_profiled_normalizations": groups,
    }


def run_mle_vs_binned_chi2_diagnostic(
    *, cache_path: Path, run_info_path: Path, dilution_json_path: Path,
    output_dir: Path,
) -> None:
    """Compare polarized-only unbinned MLE and 12-bin phi chi2 fits."""
    print("[MLE-vs-chi2] START: u1=u2=0, 12 phi bins", flush=True)
    out = output_dir / "diagnostics" / "mle_vs_binned_chi2"
    plots = out / "plots"
    tables = out / "tables"
    ensure_directory(plots)
    ensure_directory(tables)
    events = load_event_cache(cache_path)
    run_states = run_state_arrays(parse_run_info_csv(run_info_path))
    dilution_records = load_dilution_factors(dilution_json_path, cut_label="nominal")

    rows: list[dict[str, Any]] = []
    for bin_number in range(1, NUMBER_OF_BINS + 1):
        mle = fit_one_variant(
            events, run_states, dilution_records, bin_number, "nominal",
            fixed_physics_parameters={"u1": 0.0, "u2": 0.0},
        )
        chi2 = _binned_phi_chi2_fit(
            events, run_states, dilution_records, bin_number,
        )
        row: dict[str, Any] = {
            "bin_number": bin_number,
            "mle_valid": mle["valid"], "chi2_valid": chi2["valid"],
            "chi2": chi2["chi2"], "ndf": chi2["ndf"],
            "chi2_per_ndf": chi2["chi2_per_ndf"],
        }
        for parameter in MLE_CHI2_PARAMETERS:
            row[f"mle_{parameter}"] = mle["values"][parameter]
            row[f"mle_{parameter}_stat"] = mle["errors"][parameter]
            row[f"chi2_{parameter}"] = chi2["values"][parameter]
            row[f"chi2_{parameter}_stat"] = chi2["errors"][parameter]
        # endfor
        rows.append(row)
        print(
            f"[MLE-vs-chi2] bin {bin_number:02d}/24 | "
            f"chi2/ndf={chi2['chi2_per_ndf']:.3f} | "
            f"MLE valid={mle['valid']} | chi2 valid={chi2['valid']}",
            flush=True,
        )
    # endfor

    result = pd.DataFrame(rows)
    result.to_csv(tables / "mle_vs_12bin_phi_chi2.csv", index=False)
    x = result["bin_number"].to_numpy(int)
    for parameter in MLE_CHI2_PARAMETERS:
        fig, ax = plt.subplots(figsize=(10.5, 5.5))
        ax.errorbar(
            x - 0.08, result[f"mle_{parameter}"],
            yerr=result[f"mle_{parameter}_stat"], fmt="o", capsize=3,
            color="red", label="Unbinned MLE",
        )
        ax.errorbar(
            x + 0.08, result[f"chi2_{parameter}"],
            yerr=result[f"chi2_{parameter}_stat"], fmt="s", capsize=3,
            color="blue", label=r"12-bin $\chi^2$ fit",
        )
        ax.axhline(0.0, lw=0.8, color="black")
        ax.set_xlabel("Analysis bin")
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_xticks(x)
        limits = PARAMETER_Y_LIMITS.get(parameter)
        if limits is not None:
            ax.set_ylim(*limits)
        # endif
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(plots / f"mle_vs_chi2_{parameter}.png", dpi=200)
        plt.close(fig)
    # endfor
    print(f"[MLE-vs-chi2] wrote {tables / 'mle_vs_12bin_phi_chi2.csv'}", flush=True)
    print(f"[MLE-vs-chi2] wrote five comparison plots to {plots}", flush=True)



def run_double_spin_target_split_diagnostic(
    *, cache_path: Path, run_info_path: Path, dilution_json_path: Path,
    output_dir: Path, skip_plots: bool = False,
) -> int:
    """Extract LL independently for the two target-polarization signs.

    A single target orientation cannot determine UU and UL separately.  For
    this diagnostic the nominal simultaneous-fit UU/UL amplitudes (u1, u2,
    ul1, ul2) are therefore fixed, while lu1, ll0 and ll1 are refitted using
    only one target-polarization sign at a time.  The comparison of interest
    is ll0 and ll1 between the P_t>0 and P_t<0 samples.
    """
    print("[double-spin target split] START", flush=True)
    out = output_dir / "diagnostics" / "double_spin_target_split"
    plots = out / "plots"
    tables = out / "tables"
    ensure_directory(plots)
    ensure_directory(tables)

    if not cache_path.is_file():
        raise FileNotFoundError(
            "The target-split diagnostic uses the nominal selected-event cache, "
            f"but it was not found: {cache_path}. Run the normal extraction once first."
        )
    # endif

    nominal_table_path = output_dir / "nominal" / "tables" / "structure_function_ratios.csv"
    if not nominal_table_path.is_file():
        raise FileNotFoundError(
            "The target-split diagnostic fixes the UU/UL denominator terms to "
            "the completed nominal simultaneous fit, but that table was not found: "
            f"{nominal_table_path}. Run the normal extraction once first."
        )
    # endif
    nominal_frame = pd.read_csv(nominal_table_path).set_index("bin_number")
    required_columns = ("u1", "u2", "ul1", "ul2")
    missing_columns = [name for name in required_columns if name not in nominal_frame.columns]
    if missing_columns:
        raise RuntimeError(
            "Nominal table is missing columns required by the target-split diagnostic: "
            + ", ".join(missing_columns)
        )
    # endif

    events = load_event_cache(cache_path)
    run_states = run_state_arrays(parse_run_info_csv(run_info_path))
    dilution_records = load_dilution_factors(dilution_json_path, cut_label="nominal")

    rows: list[dict[str, Any]] = []
    for bin_number in range(1, NUMBER_OF_BINS + 1):
        nominal_row = nominal_frame.loc[bin_number]
        fixed_denominator = {
            name: float(nominal_row[name])
            for name in required_columns
        }
        initial_values = {
            name: float(nominal_row[name])
            for name in PHYSICS_PARAMETERS
            if name in nominal_frame.columns
        }

        fits: dict[int, dict[str, Any]] = {}
        for target_sign in (1, -1):
            fit = fit_one_variant(
                events, run_states, dilution_records, bin_number, "nominal",
                initial_values=initial_values,
                fixed_physics_parameters=fixed_denominator,
                target_sign_filter=target_sign,
            )
            fits[target_sign] = fit
        # endfor

        row: dict[str, Any] = {
            "bin_number": bin_number,
            **{f"fixed_{name}": value for name, value in fixed_denominator.items()},
        }
        for target_sign, tag in ((1, "target_plus"), (-1, "target_minus")):
            fit = fits[target_sign]
            row[f"{tag}_valid"] = bool(fit["valid"])
            row[f"{tag}_events"] = int(fit["metadata"]["number_of_events"])
            # Only LU and LL are independently refitted in a single target state.
            for parameter in ("lu1", "ll0", "ll1"):
                row[f"{tag}_{parameter}"] = float(fit["values"][parameter])
                row[f"{tag}_{parameter}_stat"] = float(fit["errors"][parameter])
            # endfor
        # endfor
        for parameter in ("ll0", "ll1"):
            vp = row[f"target_plus_{parameter}"]
            vm = row[f"target_minus_{parameter}"]
            ep = row[f"target_plus_{parameter}_stat"]
            em = row[f"target_minus_{parameter}_stat"]
            sigma = math.sqrt(ep * ep + em * em)
            row[f"delta_{parameter}"] = vp - vm
            row[f"pull_{parameter}"] = (vp - vm) / sigma if sigma > 0.0 else math.nan
        # endfor
        rows.append(row)
        print(
            f"[double-spin target split] bin {bin_number:02d}/24 | "
            f"A_LL: +={row['target_plus_ll0']:+.4f}, -={row['target_minus_ll0']:+.4f} | "
            f"A_LL^cosphi: +={row['target_plus_ll1']:+.4f}, -={row['target_minus_ll1']:+.4f}",
            flush=True,
        )
    # endfor

    frame = pd.DataFrame(rows)
    csv_path = tables / "double_spin_target_split.csv"
    frame.to_csv(csv_path, index=False)

    x = frame["bin_number"].to_numpy(float)
    for parameter in (() if skip_plots else ("ll0", "ll1")):
        fig, ax = plt.subplots(figsize=(12, 6.5))
        ax.errorbar(
            x - 0.08, frame[f"target_plus_{parameter}"],
            yerr=frame[f"target_plus_{parameter}_stat"],
            fmt="o", capsize=3, label=r"$P_t>0$",
        )
        ax.errorbar(
            x + 0.08, frame[f"target_minus_{parameter}"],
            yerr=frame[f"target_minus_{parameter}_stat"],
            fmt="s", capsize=3, label=r"$P_t<0$",
        )
        ax.axhline(0.0, linewidth=1.0)
        ax.set_xlabel("Analysis bin")
        ax.set_ylabel(PARAMETER_LABELS[parameter])
        ax.set_xticks(np.arange(1, NUMBER_OF_BINS + 1))
        limits = PARAMETER_Y_LIMITS.get(parameter)
        if limits is not None:
            ax.set_ylim(*limits)
        # endif
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(plots / f"double_spin_target_split_{parameter}.png", dpi=200)
        plt.close(fig)
    # endfor

    for parameter in ("ll0", "ll1"):
        pulls = frame[f"pull_{parameter}"].to_numpy(float)
        pulls = pulls[np.isfinite(pulls)]
        if pulls.size:
            print(
                f"[double-spin target split] {parameter}: "
                f"mean pull={np.mean(pulls):+.3f}, RMS={np.sqrt(np.mean(pulls**2)):.3f}, "
                f"max |pull|={np.max(np.abs(pulls)):.3f}",
                flush=True,
            )
        # endif
    # endfor
    print(
        "[double-spin target split] fixed nominal amplitudes: u1, u2, ul1, ul2; "
        "refitted amplitudes: lu1, ll0, ll1",
        flush=True,
    )
    print(f"[double-spin target split] table: {csv_path}", flush=True)
    print(f"[double-spin target split] plots: {plots}", flush=True)
    return 0



SOLENOID_SPLIT_CONFIG = {
    "negative": {
        "label": r"Solenoid $-1$",
        "run_ranges": ((16043, 16772), (16843, 17183), (17477, 17811)),
        "active_periods": ("su22", "fa22", "sp23"),
    },
    "positive": {
        "label": r"Solenoid $+1$",
        "run_ranges": ((17185, 17408),),
        "active_periods": ("fa22",),
    },
}


def fit_xb_integrated_solenoid_tprime(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    t_index: int,
    solenoid_tag: str,
) -> dict[str, Any]:
    """Fit all seven amplitudes in one -t' bin, integrated over xB."""
    config = SOLENOID_SPLIT_CONFIG[solenoid_tag]
    active_periods = tuple(config["active_periods"])
    run_ranges = tuple(config["run_ranges"])
    bin_numbers = [
        combined_bin_number(x_index, t_index)
        for x_index in range(len(XB_BINS))
    ]
    sub_nlls = {
        bin_number: make_bin_nll(
            events, run_states, dilution_records, bin_number, "nominal",
            active_periods=active_periods,
            run_ranges_filter=run_ranges,
        )[0]
        for bin_number in bin_numbers
    }

    physics_names = list(PHYSICS_PARAMETERS)
    dilution_names = [
        f"f_{period}_bin{bin_number:02d}"
        for bin_number in bin_numbers
        for period in PERIODS
    ]
    names = physics_names + dilution_names
    start = [float(PARAMETER_INITIAL_VALUES[name]) for name in physics_names]
    for bin_number in bin_numbers:
        for period in PERIODS:
            start.append(float(dilution_records[(period, bin_number)].value))
        # endfor
    # endfor

    def combined_nll(*pars: float) -> float:
        values = dict(zip(names, pars))
        total = 0.0
        for bin_number in bin_numbers:
            dilution_args = [
                values[f"f_{period}_bin{bin_number:02d}"]
                for period in PERIODS
            ]
            total += sub_nlls[bin_number](
                values["u1"], values["u2"], values["lu1"],
                values["ul1"], values["ul2"], values["ll0"], values["ll1"],
                *dilution_args,
            )
        # endfor
        return float(total)

    def configured_minuit(start_values: list[float]) -> Minuit:
        candidate = Minuit(combined_nll, *start_values, name=names)
        candidate.errordef = Minuit.LIKELIHOOD
        candidate.print_level = 0
        candidate.strategy = 1
        for name in physics_names:
            candidate.limits[name] = PARAMETER_LIMITS[name]
        # endfor
        for bin_number in bin_numbers:
            for period in PERIODS:
                name = f"f_{period}_bin{bin_number:02d}"
                record = dilution_records[(period, bin_number)]
                width = max(8.0 * record.stat_uncertainty, 0.20 * record.value, 0.02)
                candidate.limits[name] = (
                    max(1.0e-6, record.value - width), record.value + width
                )
                if period not in active_periods or record.stat_uncertainty == 0.0:
                    candidate.fixed[name] = True
                # endif
            # endfor
        # endfor
        return candidate

    starts = [list(start)]
    zero_physics = list(start)
    for i, name in enumerate(physics_names):
        zero_physics[i] = 0.0
    # endfor
    starts.append(zero_physics)
    attempted = []
    for start_values in starts:
        candidate = configured_minuit(start_values)
        candidate.migrad(ncall=80000)
        if not candidate.fmin.is_valid:
            candidate.simplex(ncall=40000)
            candidate.strategy = 2
            candidate.migrad(ncall=120000)
        # endif
        candidate.hesse()
        attempted.append(candidate)
    # endfor

    def quality(candidate: Minuit) -> tuple[int, float, float]:
        valid = candidate.fmin.is_valid and math.isfinite(float(candidate.fval))
        return (
            0 if valid else 1,
            float(candidate.fval) if math.isfinite(float(candidate.fval)) else math.inf,
            float(candidate.fmin.edm) if math.isfinite(float(candidate.fmin.edm)) else math.inf,
        )

    best = min(attempted, key=quality)
    values = {name: float(best.values[name]) for name in PHYSICS_PARAMETERS}
    errors = {name: float(best.errors[name]) for name in PHYSICS_PARAMETERS}
    event_mask = np.isin(events["bin_number"], np.asarray(bin_numbers, dtype=np.int16))
    run_mask = np.zeros(event_mask.shape, dtype=bool)
    for run_low, run_high in run_ranges:
        run_mask |= (events["runnum"] >= run_low) & (events["runnum"] <= run_high)
    # endfor
    event_mask &= run_mask
    return {
        "t_index": t_index,
        "solenoid": solenoid_tag,
        "bin_numbers": bin_numbers,
        "number_of_events": int(np.count_nonzero(event_mask)),
        "values": values,
        "errors": errors,
        "valid": bool(best.fmin.is_valid),
        "edm": float(best.fmin.edm),
        "fval": float(best.fval),
    }


def _solenoid_split_worker(task: tuple[str, int, str, str, str]) -> dict[str, Any]:
    cache_text, run_info_text, dilution_text, solenoid_tag, t_index_text = task
    events = load_event_cache(Path(cache_text))
    run_states = run_state_arrays(parse_run_info_csv(Path(run_info_text)))
    dilution_records = load_dilution_factors(Path(dilution_text), cut_label="nominal")
    return fit_xb_integrated_solenoid_tprime(
        events, run_states, dilution_records, int(t_index_text), solenoid_tag
    )


def run_solenoid_split_diagnostic(
    *, cache_path: Path, run_info_path: Path, dilution_json_path: Path,
    output_dir: Path, workers: int, skip_plots: bool = False,
) -> int:
    """Compare all seven amplitudes for the two solenoid polarities in six -t' bins."""
    print("[solenoid split] START: xB-integrated, six -t' bins", flush=True)
    out = output_dir / "diagnostics" / "solenoid_split"
    plots = out / "plots"
    tables = out / "tables"
    ensure_directory(plots)
    ensure_directory(tables)
    if not cache_path.is_file():
        raise FileNotFoundError(
            f"Solenoid-split diagnostic requires the nominal cache: {cache_path}"
        )
    # endif

    tasks = [
        (str(cache_path), str(run_info_path), str(dilution_json_path), tag, str(t_index))
        for t_index in range(len(MINUS_TPRIME_BINS_GEV2))
        for tag in ("negative", "positive")
    ]
    results: list[dict[str, Any]] = []
    max_workers = max(1, min(int(workers), len(tasks)))
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_solenoid_split_worker, task): task for task in tasks}
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(
                f"[solenoid split] -t' bin {result['t_index'] + 1}/6 | "
                f"{result['solenoid']} | N={result['number_of_events']:,} | "
                f"valid={result['valid']}", flush=True,
            )
        # endfor
    # endwith

    by_key = {(item["t_index"], item["solenoid"]): item for item in results}
    rows: list[dict[str, Any]] = []
    for t_index, (t_low, t_high) in enumerate(MINUS_TPRIME_BINS_GEV2):
        neg = by_key[(t_index, "negative")]
        pos = by_key[(t_index, "positive")]
        row: dict[str, Any] = {
            "t_index": t_index,
            "t_bin": t_index + 1,
            "minus_tprime_low": float(t_low),
            "minus_tprime_high": float(t_high),
            "solenoid_negative_events": neg["number_of_events"],
            "solenoid_positive_events": pos["number_of_events"],
            "solenoid_negative_valid": neg["valid"],
            "solenoid_positive_valid": pos["valid"],
        }
        for parameter in PHYSICS_PARAMETERS:
            vn, en = neg["values"][parameter], neg["errors"][parameter]
            vp, ep = pos["values"][parameter], pos["errors"][parameter]
            sigma = math.sqrt(en * en + ep * ep)
            row[f"solenoid_negative_{parameter}"] = vn
            row[f"solenoid_negative_{parameter}_stat"] = en
            row[f"solenoid_positive_{parameter}"] = vp
            row[f"solenoid_positive_{parameter}_stat"] = ep
            row[f"delta_{parameter}"] = vp - vn
            row[f"pull_{parameter}"] = (vp - vn) / sigma if sigma > 0.0 else math.nan
        # endfor
        rows.append(row)
    # endfor

    frame = pd.DataFrame(rows)
    csv_path = tables / "solenoid_split.csv"
    frame.to_csv(csv_path, index=False)
    x = frame["t_bin"].to_numpy(float)
    if not skip_plots:
        for parameter in PHYSICS_PARAMETERS:
            fig, ax = plt.subplots(figsize=(10.5, 5.8))
            ax.errorbar(
                x - 0.06, frame[f"solenoid_negative_{parameter}"],
                yerr=frame[f"solenoid_negative_{parameter}_stat"],
                fmt="o", capsize=3, label=r"Solenoid $-1$",
            )
            ax.errorbar(
                x + 0.06, frame[f"solenoid_positive_{parameter}"],
                yerr=frame[f"solenoid_positive_{parameter}_stat"],
                fmt="s", capsize=3, label=r"Solenoid $+1$",
            )
            ax.axhline(0.0, linewidth=1.0)
            ax.set_xlabel(r"$-t'$ bin (integrated over $x_B$)")
            ax.set_ylabel(PARAMETER_LABELS[parameter])
            ax.set_xticks(np.arange(1, len(MINUS_TPRIME_BINS_GEV2) + 1))
            limits = PARAMETER_Y_LIMITS.get(parameter)
            if limits is not None:
                ax.set_ylim(*limits)
            # endif
            ax.legend(frameon=False)
            fig.tight_layout()
            fig.savefig(plots / f"solenoid_split_{parameter}.png", dpi=200)
            plt.close(fig)
        # endfor
    # endif

    for parameter in PHYSICS_PARAMETERS:
        pulls = frame[f"pull_{parameter}"].to_numpy(float)
        pulls = pulls[np.isfinite(pulls)]
        if pulls.size:
            print(
                f"[solenoid split] {parameter}: mean pull={np.mean(pulls):+.3f}, "
                f"RMS={np.sqrt(np.mean(pulls**2)):.3f}, "
                f"chi2/{pulls.size}={np.sum(pulls**2):.2f}/{pulls.size}",
                flush=True,
            )
        # endif
    # endfor
    print(f"[solenoid split] table: {csv_path}", flush=True)
    print(f"[solenoid split] plots: {plots}", flush=True)
    return 0


# =============================================================================
# Appendix MLE fit-quality diagnostics
# =============================================================================

def _appendix_profile_worker(task: dict[str, Any]) -> list[dict[str, Any]]:
    """Profile all seven parameters for one bin using one shared nominal MLE.

    One process owns one physics bin.  This avoids repeating the expensive nominal
    seven-parameter fit once per profile parameter while retaining warm starts
    along each one-dimensional profile scan.
    """
    if (_WORKER_EVENTS is None or _WORKER_RUN_STATES is None
            or _WORKER_DILUTION_RECORDS is None):
        raise RuntimeError("Appendix diagnostic worker was not initialized.")
    # endif
    bin_number = int(task["bin_number"])
    points = int(task["profile_points"])
    fit = fit_one_variant(
        _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
        bin_number, "nominal", active_periods=PERIODS,
    )
    results: list[dict[str, Any]] = []
    if not fit["valid"]:
        for parameter in PHYSICS_PARAMETERS:
            results.append({"bin_number": bin_number, "parameter": parameter,
                            "fit": fit, "profile": None})
        # endfor
        return results
    # endif

    for parameter in PHYSICS_PARAMETERS:
        center = float(fit["values"][parameter])
        sigma = max(float(fit["errors"][parameter]), 1.0e-3)
        low_limit, high_limit = PARAMETER_LIMITS[parameter]
        half_width = max(8.0 * sigma, 0.35)
        low = max(low_limit, center - half_width)
        high = min(high_limit, center + half_width)
        grid = np.linspace(low, high, points)
        nll_values = np.full(points, math.nan, dtype=float)
        valid_values = np.zeros(points, dtype=bool)

        # Walk outward from the nominal MLE independently on the two sides.
        # This gives every constrained fit a nearby warm start instead of
        # dragging a solution all the way from one edge of the profile to the other.
        center_index = int(np.argmin(np.abs(grid - center)))
        walk_orders = [
            list(range(center_index, -1, -1)),
            list(range(center_index + 1, points)),
        ]
        for walk_order in walk_orders:
            start_values = dict(fit["values"])
            for grid_index in walk_order:
                value = float(grid[grid_index])
                profiled = fit_one_variant(
                    _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
                    bin_number, "nominal", active_periods=PERIODS,
                    initial_values=start_values,
                    fixed_physics_parameters={parameter: value},
                    fast_profile=True,
                )
                nll_values[grid_index] = float(profiled["minimum_nll"])
                valid_values[grid_index] = bool(profiled["valid"])
                if profiled["valid"]:
                    start_values = dict(profiled["values"])
                # endif
            # endfor
        # endfor
        finite = np.isfinite(nll_values) & (np.asarray(nll_values) < INVALID_NLL / 10.0)
        minimum = min([float(fit["minimum_nll"])]
                      + [float(nll_values[i]) for i in range(points) if finite[i]])
        two_delta = [2.0 * (float(value) - minimum) if finite[i] else math.nan
                     for i, value in enumerate(nll_values)]
        results.append({
            "bin_number": bin_number, "parameter": parameter, "fit": fit,
            "profile": {"grid": grid.tolist(), "two_delta_nll": two_delta,
                        "valid": valid_values.tolist()},
        })
    # endfor
    return results


def _profile_crossing(grid: np.ndarray, delta: np.ndarray, center: float,
                      level: float, side: str) -> float | None:
    finite = np.isfinite(grid) & np.isfinite(delta)
    x = grid[finite]; y = delta[finite]
    if side == "low":
        mask = x <= center
        x = x[mask][::-1]; y = y[mask][::-1]
    else:
        mask = x >= center
        x = x[mask]; y = y[mask]
    # endif
    if x.size < 2:
        return None
    # endif
    for i in range(x.size - 1):
        y0, y1 = y[i], y[i + 1]
        if (y0 - level) * (y1 - level) <= 0.0 and y1 != y0:
            fraction = (level - y0) / (y1 - y0)
            return float(x[i] + fraction * (x[i + 1] - x[i]))
        # endif
    # endfor
    return None


def _appendix_conditional_prediction(
    events: Mapping[str, np.ndarray], run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord], bin_number: int,
    fit: Mapping[str, Any], phi_edges: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Observed and conditional-MLE expected populations in four spin categories."""
    values = fit["values"]
    categories = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    nbins = len(phi_edges) - 1
    observed = {cat: np.zeros(nbins, dtype=float) for cat in categories}
    expected = {cat: np.zeros(nbins, dtype=float) for cat in categories}
    variance = {cat: np.zeros(nbins, dtype=float) for cat in categories}

    for period in PERIODS:
        mask = ((events["bin_number"] == bin_number)
                & (events["period_index"] == PERIOD_INDEX[period]))
        if not np.any(mask):
            continue
        # endif
        phi = events["phi"][mask].astype(float, copy=False)
        r_b = events["rB"][mask].astype(float, copy=False)
        r_c = events["rC"][mask].astype(float, copy=False)
        r_v = events["rV"][mask].astype(float, copy=False)
        r_w = events["rW"][mask].astype(float, copy=False)
        sin_theta = events["sin_theta_gamma"][mask].astype(float, copy=False)
        cos_theta = events["cos_theta_gamma"][mask].astype(float, copy=False)
        runnum = events["runnum"][mask].astype(np.int32, copy=False)
        helicity = events["helicity"][mask].astype(np.int8, copy=False)
        state = run_states[period]
        lookup = {int(run): i for i, run in enumerate(state["run"])}
        observed_state = np.fromiter((lookup[int(run)] for run in runnum),
                                     count=runnum.size, dtype=np.int32)
        observed_sign = np.sign(state["pt"][observed_state]).astype(np.int8)
        dilution = float(values[f"f_{period}"])
        beam_pol = BEAM_POLARIZATION[period]
        bin_index = np.searchsorted(phi_edges, phi, side="right") - 1
        in_range = (bin_index >= 0) & (bin_index < nbins)

        weights_by_cat: dict[tuple[int, int], np.ndarray] = {}
        denominator = np.zeros(phi.shape, dtype=float)
        for target_sign, h in categories:
            state_mask = np.sign(state["pt"]) == target_sign
            charges = state["q_plus"] if h > 0 else state["q_minus"]
            cat_weight = np.zeros(phi.shape, dtype=float)
            for j in np.flatnonzero(state_mask):
                charge = float(charges[j])
                if charge <= 0.0:
                    continue
                # endif
                factor = evaluate_cross_section_factor(
                    variant="nominal", phi=phi, r_b=r_b, r_c=r_c, r_v=r_v, r_w=r_w,
                    sin_theta_gamma=sin_theta, cos_theta_gamma=cos_theta,
                    helicity=np.full(phi.shape, h, dtype=float),
                    beam_polarization=beam_pol,
                    target_polarization=np.full(phi.shape, float(state["pt"][j])),
                    dilution=dilution, transverse_scales=None,
                    **{name: float(values[name]) for name in PHYSICS_PARAMETERS},
                )
                cat_weight += charge * factor
            # endfor
            weights_by_cat[(target_sign, h)] = cat_weight
            denominator += cat_weight
        # endfor
        if np.any(denominator <= CROSS_SECTION_FLOOR):
            raise RuntimeError(f"Non-positive appendix diagnostic denominator in bin {bin_number}, {period}.")
        # endif

        for cat in categories:
            target_sign, h = cat
            obs_mask = in_range & (observed_sign == target_sign) & (helicity == h)
            observed[cat] += np.bincount(bin_index[obs_mask], minlength=nbins)[:nbins]
            probability = weights_by_cat[cat] / denominator
            for k in range(nbins):
                k_mask = in_range & (bin_index == k)
                p = probability[k_mask]
                expected[cat][k] += float(np.sum(p))
                variance[cat][k] += float(np.sum(p * (1.0 - p)))
            # endfor
        # endfor
    # endfor

    rows: list[dict[str, Any]] = []
    residual_values: list[float] = []
    for cat in categories:
        target_sign, h = cat
        for k in range(nbins):
            sigma = math.sqrt(max(variance[cat][k], 0.0))
            residual = ((observed[cat][k] - expected[cat][k]) / sigma
                        if sigma > 0.0 else math.nan)
            if math.isfinite(residual):
                residual_values.append(residual)
            # endif
            rows.append({
                "bin_number": bin_number, "target_sign": target_sign, "helicity": h,
                "phi_bin": k, "phi_low": float(phi_edges[k]), "phi_high": float(phi_edges[k + 1]),
                "phi_center": float(0.5 * (phi_edges[k] + phi_edges[k + 1])),
                "observed": float(observed[cat][k]), "expected": float(expected[cat][k]),
                "conditional_variance": float(variance[cat][k]), "residual": residual,
            })
        # endfor
    # endfor
    summary = {
        "bin_number": bin_number,
        "number_of_events": int(fit["metadata"]["number_of_events"]),
        "residual_mean": float(np.mean(residual_values)) if residual_values else math.nan,
        "residual_rms": float(np.sqrt(np.mean(np.square(residual_values)))) if residual_values else math.nan,
        "maximum_absolute_residual": float(np.max(np.abs(residual_values))) if residual_values else math.nan,
    }
    return rows, summary


def run_appendix_mle_diagnostics(args: argparse.Namespace, root: Path, workers: int,
                                 dilution_json_path: Path) -> int:
    """Produce appendix-ready diagnostics for the production seven-parameter joint MLE."""
    out_dir = root / "appendix_mle_diagnostics"
    conditional_dir = out_dir / "conditional_fit_quality"
    profile_dir = out_dir / "profile_likelihoods"
    table_dir = out_dir / "tables"
    for directory in (out_dir, conditional_dir, profile_dir, table_dir):
        ensure_directory(directory)
    # endfor
    cache_path = (args.cache.expanduser().resolve() if args.cache
                  else root / "nominal/cache/selected_events.npz")
    if not cache_path.is_file():
        raise FileNotFoundError(
            f"Nominal selected-event cache not found: {cache_path}. Run the nominal extraction first or pass --cache."
        )
    # endif
    events = load_event_cache(cache_path)
    run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    run_states = run_state_arrays(run_records)
    dilution_records = load_dilution_factors(dilution_json_path, cut_label="nominal")
    run_state_payload = {p: {k: v.tolist() for k, v in s.items()} for p, s in run_states.items()}
    dilution_payload = {
        p: {str(b): {"x_index": r.x_index, "t_index": r.t_index, "value": r.value,
                     "stat_uncertainty": r.stat_uncertainty}
            for (rp, b), r in dilution_records.items() if rp == p}
        for p in PERIODS
    }
    # Parallelize by physics bin, not by (bin, parameter).  Each worker performs
    # the nominal joint fit once and reuses it for all seven profile scans.
    # With 24 bins this still provides ample coarse-grained work for eight
    # processes while eliminating 6 redundant nominal fits per bin.
    tasks = [{"bin_number": b, "profile_points": args.appendix_profile_points}
             for b in range(1, NUMBER_OF_BINS + 1)]
    results: list[dict[str, Any]] = []
    mp_context = mp.get_context("spawn")
    n_appendix_workers = max(1, min(int(workers), 8, len(tasks)))
    print(
        f"[appendix MLE] profiling 24 bins x {len(PHYSICS_PARAMETERS)} parameters "
        f"with {n_appendix_workers} workers",
        flush=True,
    )
    with ProcessPoolExecutor(
        max_workers=n_appendix_workers, mp_context=mp_context,
        initializer=initialize_fit_worker,
        initargs=(str(cache_path), run_state_payload, dilution_payload),
    ) as executor:
        futures = [executor.submit(_appendix_profile_worker, task) for task in tasks]
        completed_profiles = 0
        for index, future in enumerate(as_completed(futures), start=1):
            bin_results = future.result()
            results.extend(bin_results)
            completed_profiles += len(bin_results)
            print(
                f"[appendix MLE] bins {index}/{len(tasks)}; "
                f"profiles {completed_profiles}/{len(tasks) * len(PHYSICS_PARAMETERS)}",
                flush=True,
            )
        # endfor
    # endwith
    by_key = {(r["bin_number"], r["parameter"]): r for r in results}

    profile_rows: list[dict[str, Any]] = []
    crossing_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    for b in range(1, NUMBER_OF_BINS + 1):
        fit = by_key[(b, PHYSICS_PARAMETERS[0])]["fit"]
        fit_rows.append({
            "bin_number": b, "number_of_events": fit["metadata"]["number_of_events"],
            "valid": fit["valid"], "accurate_covariance": fit["accurate_covariance"],
            "positive_definite_covariance": fit["positive_definite_covariance"],
            "parameters_at_limit": fit["parameters_at_limit"], "edm": fit["edm"],
            "minimum_nll": fit["minimum_nll"],
            **{name: fit["values"][name] for name in PHYSICS_PARAMETERS},
            **{f"{name}_stat": fit["errors"][name] for name in PHYSICS_PARAMETERS},
        })
        for par in PHYSICS_PARAMETERS:
            item = by_key[(b, par)]; profile = item["profile"]
            if profile is None:
                continue
            # endif
            grid = np.asarray(profile["grid"], dtype=float)
            delta = np.asarray(profile["two_delta_nll"], dtype=float)
            center = float(item["fit"]["values"][par])
            for x, y, valid in zip(grid, delta, profile["valid"]):
                profile_rows.append({"bin_number": b, "parameter": par,
                                     "parameter_value": x, "two_delta_nll": y,
                                     "profile_fit_valid": valid})
            # endfor
            crossing_rows.append({
                "bin_number": b, "parameter": par, "mle": center,
                "hesse_stat": float(item["fit"]["errors"][par]),
                "low_68": _profile_crossing(grid, delta, center, 1.0, "low"),
                "high_68": _profile_crossing(grid, delta, center, 1.0, "high"),
                "low_95": _profile_crossing(grid, delta, center, 3.84, "low"),
                "high_95": _profile_crossing(grid, delta, center, 3.84, "high"),
            })
        # endfor
    # endfor
    pd.DataFrame(fit_rows).to_csv(table_dir / "joint_fit_status.csv", index=False)
    pd.DataFrame(profile_rows).to_csv(table_dir / "profile_likelihood_points.csv", index=False)
    pd.DataFrame(crossing_rows).to_csv(table_dir / "profile_likelihood_crossings.csv", index=False)

    for par in PHYSICS_PARAMETERS:
        for x_index in range(len(XB_BINS)):
            fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.8), sharey=True)
            for t_index, ax in enumerate(axes.flat):
                b = combined_bin_number(x_index, t_index)
                item = by_key[(b, par)]
                profile = item["profile"]
                if profile is None:
                    ax.text(0.5, 0.5, "invalid fit", ha="center", va="center", transform=ax.transAxes)
                    continue
                # endif
                grid = np.asarray(profile["grid"], dtype=float)
                delta = np.asarray(profile["two_delta_nll"], dtype=float)
                center = float(item["fit"]["values"][par])
                error = float(item["fit"]["errors"][par])
                ax.plot(grid, delta, marker="o", ms=2.8, lw=1.2)
                ax.axhline(1.0, ls="--", lw=0.9)
                ax.axhline(3.84, ls=":", lw=0.9)
                ax.axvline(center, lw=1.0)
                ax.axvline(center - error, ls="--", lw=0.8)
                ax.axvline(center + error, ls="--", lw=0.8)
                if par in ("ll0", "ll1", "u1", "u2"):
                    ax.set_xlim(-1.0, 1.0)
                else:
                    ax.set_xlim(-0.5, 0.5)
                # endif
                ax.set_ylim(bottom=0.0, top=max(5.0, min(12.0, np.nanmax(delta) * 1.05)))
                ax.set_title(f"Bin {b}: $-t^\\prime$ bin {t_index + 1}")
                ax.grid(alpha=0.20)
                if t_index >= 3:
                    ax.set_xlabel(PARAMETER_LABELS[par])
                # endif
                if t_index % 3 == 0:
                    ax.set_ylabel(r"$-2\Delta\ln\mathcal{L}$")
                # endif
            # endfor
            xlow, xhigh = XB_BINS[x_index]
            fig.suptitle(f"Joint-fit profile likelihood: {PARAMETER_LABELS[par]},  {xlow:.3f} < $x_B$ < {xhigh:.3f}")
            fig.tight_layout(rect=(0, 0, 1, 0.96))
            fig.savefig(profile_dir / f"profile_{par}_xbin_{x_index + 1}.png", dpi=200)
            plt.close(fig)
        # endfor
    # endfor

    phi_edges = np.linspace(0.0, 2.0 * math.pi, args.appendix_phi_bins + 1)
    conditional_rows: list[dict[str, Any]] = []
    residual_summaries: list[dict[str, Any]] = []
    state_titles = {(-1, -1): r"$P_t<0,\ h=-1$", (-1, 1): r"$P_t<0,\ h=+1$",
                    (1, -1): r"$P_t>0,\ h=-1$", (1, 1): r"$P_t>0,\ h=+1$"}
    for b in range(1, NUMBER_OF_BINS + 1):
        fit = by_key[(b, PHYSICS_PARAMETERS[0])]["fit"]
        rows, summary = _appendix_conditional_prediction(
            events, run_states, dilution_records, b, fit, phi_edges)
        conditional_rows.extend(rows); residual_summaries.append(summary)
        frame = pd.DataFrame(rows)
        fig = plt.figure(figsize=(12.0, 8.8))
        outer = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.30)
        for index, cat in enumerate([(-1, -1), (-1, 1), (1, -1), (1, 1)]):
            inner = outer[index // 2, index % 2].subgridspec(2, 1, height_ratios=(3.2, 1.0), hspace=0.04)
            ax = fig.add_subplot(inner[0]); rax = fig.add_subplot(inner[1], sharex=ax)
            sub = frame[(frame.target_sign == cat[0]) & (frame.helicity == cat[1])]
            x = sub.phi_center.to_numpy(); obs = sub.observed.to_numpy(); exp = sub.expected.to_numpy()
            ax.errorbar(x, obs, yerr=np.sqrt(np.maximum(obs, 1.0)), fmt="o", ms=4, capsize=2, label="Observed")
            ax.step(x, exp, where="mid", lw=1.5, label="Conditional MLE")
            ax.set_title(state_titles[cat]); ax.set_ylabel("Events"); ax.grid(alpha=0.20); ax.legend(frameon=False, fontsize=8)
            rax.axhline(0.0, lw=0.8); rax.plot(x, sub.residual.to_numpy(), "o", ms=3.5)
            rax.set_ylim(-3.0, 3.0); rax.set_ylabel("Pull"); rax.set_xlabel(r"$\phi$ (rad)"); rax.grid(alpha=0.20)
            plt.setp(ax.get_xticklabels(), visible=False)
        # endfor
        fig.suptitle(f"Bin {b}: conditional-likelihood data/model diagnostic")
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(conditional_dir / f"conditional_fit_bin_{b:02d}.png", dpi=200)
        plt.close(fig)
    # endfor
    pd.DataFrame(conditional_rows).to_csv(table_dir / "conditional_phi_data_model.csv", index=False)
    pd.DataFrame(residual_summaries).to_csv(table_dir / "conditional_residual_summary.csv", index=False)
    write_json(out_dir / "appendix_mle_diagnostics_manifest.json", {
        "mode": "appendix_mle_diagnostics", "fit_model": "production nominal seven-parameter joint MLE",
        "parameters": list(PHYSICS_PARAMETERS), "periods": list(PERIODS),
        "profile_levels_two_delta_nll": {"68.3_percent": 1.0, "95_percent": 3.84},
        "profile_points": int(args.appendix_profile_points), "phi_bins": int(args.appendix_phi_bins),
        "cache": str(cache_path), "dilution_json": str(dilution_json_path),
        "conditional_plot_note": "Phi binning is diagnostic only. Expected populations are sums of event-by-event conditional spin-state probabilities at the observed event kinematics.",
    })
    print("[appendix MLE] complete", flush=True)
    print(f"  Output: {out_dir}", flush=True)
    print(f"  Conditional plots: {NUMBER_OF_BINS}", flush=True)
    print(f"  Profile canvases: {len(PHYSICS_PARAMETERS) * len(XB_BINS)}", flush=True)
    return 0


# =============================================================================
# Extended run-period A_LL diagnostics
# =============================================================================

def load_dis_pbpt_spreadsheets(base_dir: Path) -> dict[int, tuple[float, float]]:
    """Load independent DIS run-by-run PbPt products supplied by the RGC group."""
    files = {
        "su22": base_dir / "PbPt_NH3_summer22_F1F221_2025-06-24.xlsx",
        "fa22": base_dir / "PbPt_NH3_fall22_F1F221_2025-06-24.xlsx",
        "sp23": base_dir / "PbPt_NH3_spring23_F1F221_2025-06-24.xlsx",
    }
    result: dict[int, tuple[float, float]] = {}
    for period, path in files.items():
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing DIS PbPt spreadsheet for {period}: {path}"
            )
        # endif
        frame = pd.read_excel(path)
        required = {"run", "PbPt", "dPbPt"}
        if not required.issubset(frame.columns):
            raise RuntimeError(
                f"DIS PbPt spreadsheet {path} lacks columns {sorted(required)}."
            )
        # endif
        for row in frame.itertuples(index=False):
            run = int(getattr(row, "run"))
            pbpt = float(getattr(row, "PbPt"))
            dpbpt = float(getattr(row, "dPbPt"))
            if math.isfinite(pbpt):
                result[run] = (pbpt, dpbpt)
            # endif
        # endfor
    # endfor
    return result


def build_period_stability_diagnostic_cache(
    input_paths: Mapping[str, Path],
    tree_name: str,
    chunk_size: str,
    run_records: Mapping[int, RunRecord],
    cache_path: Path,
) -> dict[str, np.ndarray]:
    """Build a production-phase-space cache without any Mx2 requirement."""
    if cache_path.is_file():
        cached = load_event_cache(cache_path)
        needed = {"e_phi", "p_phi", "p_p", "p_theta", "y", "Mx2"}
        if needed.issubset(cached):
            print(f"[extended diagnostics] reusing {cache_path}", flush=True)
            return cached
        # endif
    # endif

    collected: dict[str, list[np.ndarray]] = {
        key: [] for key in (
            "period_index", "runnum", "helicity", "bin_number", "xB",
            "minus_tprime", "minus_t", "W", "Mx2", "phi", "Q2", "y",
            "epsilon", "DepA", "DepB", "DepC", "DepV", "DepW",
            "sin_theta_gamma", "cos_theta_gamma", "rB", "rC", "rV", "rW",
            "e_phi", "p_phi", "p_p", "p_theta",
        )
    }
    for period in PERIODS:
        path = input_paths[period].expanduser().resolve()
        print(f"[extended diagnostics] reading {period}: {path}", flush=True)
        with uproot.open(path) as root_file:
            tree = resolve_tree(root_file, tree_name, path)
            aliases = dict(BRANCH_ALIASES)
            aliases.update(PERIOD_DIAGNOSTIC_BRANCH_ALIASES)
            branches = {
                logical: resolve_branch(tree, choices)
                for logical, choices in aliases.items()
            }
            expressions = list(dict.fromkeys(branches.values()))
            for arrays in tree.iterate(
                expressions=expressions, step_size=chunk_size, library="np"
            ):
                runnum = np.asarray(arrays[branches["runnum"]], dtype=np.int64)
                helicity = np.asarray(arrays[branches["helicity"]], dtype=np.int8)
                x_b = np.asarray(arrays[branches["xB"]], dtype=float)
                mtp = -np.asarray(arrays[branches["tprime"]], dtype=float)
                mt = -np.asarray(arrays[branches["t"]], dtype=float)
                w = np.asarray(arrays[branches["W"]], dtype=float)
                mx2 = np.asarray(arrays[branches["Mx2"]], dtype=float)
                phi = np.mod(
                    angle_to_radians(
                        np.asarray(arrays[branches["phi"]], dtype=float),
                        branches["phi"],
                    ),
                    2.0 * math.pi,
                )
                q2 = np.asarray(arrays[branches["Q2"]], dtype=float)
                y = np.asarray(arrays[branches["y"]], dtype=float)
                e_p = np.asarray(arrays[branches["e_p"]], dtype=float)
                e_theta = angle_to_radians(
                    np.asarray(arrays[branches["e_theta"]], dtype=float),
                    branches["e_theta"],
                )
                dep_a = np.asarray(arrays[branches["DepA"]], dtype=float)
                dep_b = np.asarray(arrays[branches["DepB"]], dtype=float)
                dep_c = np.asarray(arrays[branches["DepC"]], dtype=float)
                dep_v = np.asarray(arrays[branches["DepV"]], dtype=float)
                dep_w = np.asarray(arrays[branches["DepW"]], dtype=float)
                gamma2 = 4.0 * PROTON_MASS_GEV**2 * x_b**2 / q2
                epsilon = (1.0 - y - 0.25 * gamma2 * y**2) / (
                    1.0 - y + 0.5 * y**2 + 0.25 * gamma2 * y**2
                )
                sin_tg, cos_tg = compute_theta_gamma(
                    e_p, e_theta, BEAM_ENERGY_GEV[period]
                )
                _, _, bins = bin_indices(x_b, mtp)
                known_active = np.fromiter(
                    (
                        int(run) in run_records
                        and run_records[int(run)].period == period
                        and run_records[int(run)].charge_plus > 0.0
                        and run_records[int(run)].charge_minus > 0.0
                        for run in runnum
                    ),
                    count=runnum.size,
                    dtype=bool,
                )
                finite = (
                    np.isfinite(x_b) & np.isfinite(mtp) & np.isfinite(mt)
                    & np.isfinite(w) & np.isfinite(mx2) & np.isfinite(phi)
                    & np.isfinite(q2) & np.isfinite(y) & np.isfinite(epsilon)
                    & np.isfinite(dep_a) & np.isfinite(dep_b) & np.isfinite(dep_c)
                    & np.isfinite(dep_v) & np.isfinite(dep_w)
                )
                selected = (
                    finite & known_active & np.isin(helicity, (-1, 1))
                    & (bins >= 1) & (w > DIS_W_MIN_GEV) & (dep_a > 0.0)
                )
                if not np.any(selected):
                    continue
                # endif
                def put(name: str, values: np.ndarray) -> None:
                    collected[name].append(np.asarray(values)[selected])
                # enddef
                put("period_index", np.full(runnum.shape, PERIOD_INDEX[period], dtype=np.int8))
                put("runnum", runnum.astype(np.int32))
                put("helicity", helicity)
                put("bin_number", bins.astype(np.int16))
                put("xB", x_b); put("minus_tprime", mtp); put("minus_t", mt)
                put("W", w); put("Mx2", mx2); put("phi", phi); put("Q2", q2); put("y", y)
                put("epsilon", epsilon); put("DepA", dep_a); put("DepB", dep_b)
                put("DepC", dep_c); put("DepV", dep_v); put("DepW", dep_w)
                put("sin_theta_gamma", sin_tg); put("cos_theta_gamma", cos_tg)
                put("rB", dep_b / dep_a); put("rC", dep_c / dep_a)
                put("rV", dep_v / dep_a); put("rW", dep_w / dep_a)
                for name in ("e_phi", "p_phi", "p_theta"):
                    raw = np.asarray(arrays[branches[name]], dtype=float)
                    put(name, angle_to_radians(raw, branches[name]))
                # endfor
                put("p_p", np.asarray(arrays[branches["p_p"]], dtype=float))
            # endfor
        # endwith
    # endfor
    cache = {key: np.concatenate(value) for key, value in collected.items()}
    ensure_directory(cache_path.parent)
    np.savez_compressed(cache_path, **cache)
    print(
        f"[extended diagnostics] cache contains {cache['runnum'].size:,} events",
        flush=True,
    )
    return cache


def select_nominal_mx2_from_diagnostic_cache(
    events: Mapping[str, np.ndarray],
    cuts: Mapping[tuple[str, int], CutRecord],
) -> dict[str, np.ndarray]:
    mask = np.zeros(events["runnum"].shape, dtype=bool)
    for period in PERIODS:
        pmask = events["period_index"] == PERIOD_INDEX[period]
        for bin_number in range(1, NUMBER_OF_BINS + 1):
            cut = cuts[(period, bin_number)]
            mask |= (
                pmask & (events["bin_number"] == bin_number)
                & (events["Mx2"] >= cut.low_gev2)
                & (events["Mx2"] < cut.high_gev2)
            )
        # endfor
    # endfor
    return {name: np.asarray(values)[mask] for name, values in events.items()}


def _fixed_other_physics(frame: pd.DataFrame, bin_number: int) -> dict[str, float]:
    row = frame.loc[frame["bin_number"].astype(int) == int(bin_number)].iloc[0]
    return {
        parameter: float(row[parameter])
        for parameter in PHYSICS_PARAMETERS
        if parameter != "ll0"
    }


def _unit_dilutions() -> dict[tuple[str, int], DilutionRecord]:
    result = {}
    for period in PERIODS:
        for b in range(1, NUMBER_OF_BINS + 1):
            xi = (b - 1) // len(MINUS_TPRIME_BINS_GEV2)
            ti = (b - 1) % len(MINUS_TPRIME_BINS_GEV2)
            result[(period, b)] = DilutionRecord(period, b, xi, ti, 1.0, 0.0)
        # endfor
    # endfor
    return result


def write_dis_pbpt_all_diagnostic(
    frame: pd.DataFrame,
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    pbpt: Mapping[int, tuple[float, float]],
    output_dir: Path,
) -> Path:
    ensure_directory(output_dir)
    products = {run: value[0] for run, value in pbpt.items()}
    rows = []
    for b in range(1, NUMBER_OF_BINS + 1):
        fixed = _fixed_other_physics(frame, b)
        initial = {"ll0": float(frame.loc[frame["bin_number"] == b, "ll0"].iloc[0])}
        for period in PERIODS:
            period_runs = set(np.asarray(run_states[period]["run"], dtype=int))
            missing = sorted(run for run in period_runs if run not in products)
            if missing:
                print(
                    f"[DIS PbPt] {period}: {len(missing)} run-info runs lack DIS PbPt; "
                    "only events in matched runs are retained",
                    flush=True,
                )
            # endif
            matched_ranges = tuple((run, run) for run in sorted(period_runs & products.keys()))
            try:
                fit = fit_one_variant(
                    events, run_states, dilution_records, b, "nominal",
                    active_periods=(period,), initial_values=initial,
                    fixed_physics_parameters=fixed,
                    run_ranges_filter=matched_ranges,
                    double_spin_products_by_run=products,
                )
                rows.append({
                    "bin_number": b, "period": period,
                    "ll0_dis_pbpt": fit["values"]["ll0"],
                    "ll0_dis_pbpt_stat": fit["errors"]["ll0"],
                    "valid": fit["valid"], "edm": fit["edm"],
                    "number_of_events": fit["metadata"]["number_of_events"],
                })
            except Exception as exc:
                rows.append({"bin_number": b, "period": period, "error": str(exc)})
            # endtry
        # endfor
    # endfor
    out = output_dir / "all_period_stability_dis_pbpt.csv"
    result = pd.DataFrame(rows)
    result.to_csv(out, index=False)
    if not result.empty:
        fig, ax = plt.subplots(figsize=(11.0, 5.2))
        for i, period in enumerate(PERIODS):
            sub = result[result["period"] == period]
            ax.errorbar(sub["bin_number"] + 0.12 * (i - 1), sub["ll0_dis_pbpt"], yerr=sub["ll0_dis_pbpt_stat"], fmt="o", label=period)
        # endfor
        ax.axhline(0.0, linewidth=1.0, linestyle="--")
        ax.set_xlabel("Kinematic bin")
        ax.set_ylabel(r"$A_{LL}$ using DIS run-by-run $P_bP_t$")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "all_period_stability_dis_pbpt.png", dpi=180)
        plt.close(fig)
    # endif
    return out


def write_inclusive_epi_control_diagnostic(
    frame: pd.DataFrame,
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    output_dir: Path,
) -> Path:
    """Fit inclusive e pi+ X A_raw/(Pb Pt) with no dilution-factor correction."""
    ensure_directory(output_dir)
    unit_f = _unit_dilutions()
    rows = []
    zero_fixed = {p: 0.0 for p in PHYSICS_PARAMETERS if p != "ll0"}
    for b in range(1, NUMBER_OF_BINS + 1):
        for period in PERIODS:
            try:
                fit = fit_one_variant(
                    events, run_states, unit_f, b, "nominal",
                    active_periods=(period,), fixed_physics_parameters=zero_fixed,
                )
                rows.append({
                    "bin_number": b, "period": period,
                    "all_over_pbpt": fit["values"]["ll0"],
                    "all_over_pbpt_stat": fit["errors"]["ll0"],
                    "valid": fit["valid"], "edm": fit["edm"],
                    "number_of_events": fit["metadata"]["number_of_events"],
                })
            except Exception as exc:
                rows.append({"bin_number": b, "period": period, "error": str(exc)})
            # endtry
        # endfor
    # endfor
    out = output_dir / "inclusive_epiX_all_over_pbpt.csv"
    result = pd.DataFrame(rows)
    result.to_csv(out, index=False)
    if not result.empty:
        fig, ax = plt.subplots(figsize=(11.0, 5.2))
        for i, period in enumerate(PERIODS):
            sub = result[result["period"] == period]
            ax.errorbar(sub["bin_number"] + 0.12 * (i - 1), sub["all_over_pbpt"], yerr=sub["all_over_pbpt_stat"], fmt="o", label=period)
        # endfor
        ax.axhline(0.0, linewidth=1.0, linestyle="--")
        ax.set_xlabel("Kinematic bin")
        ax.set_ylabel(r"Inclusive $e\pi^+X$: $A_{LL}^{\rm raw}/(P_bP_t)$ (no dilution correction)")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "inclusive_epiX_all_over_pbpt.png", dpi=180)
        plt.close(fig)
    # endif
    return out


def fd_sector_from_phi_rad_array(phi: np.ndarray) -> np.ndarray:
    deg = np.mod(np.degrees(np.asarray(phi, dtype=float)), 360.0)
    sector = np.full(deg.shape, -1, dtype=np.int8)
    sector[(deg >= 330.0) | (deg < 30.0)] = 1
    sector[(deg >= 30.0) & (deg < 90.0)] = 2
    sector[(deg >= 90.0) & (deg < 150.0)] = 3
    sector[(deg >= 150.0) & (deg < 210.0)] = 4
    sector[(deg >= 210.0) & (deg < 270.0)] = 5
    sector[(deg >= 270.0) & (deg < 330.0)] = 6
    return sector


def _inverse_variance_combine(values: list[tuple[float, float]]) -> tuple[float, float, int]:
    good = [(v, e) for v, e in values if math.isfinite(v) and math.isfinite(e) and e > 0.0]
    if not good:
        return math.nan, math.nan, 0
    # endif
    w = np.asarray([1.0 / e**2 for _, e in good], dtype=float)
    v = np.asarray([x for x, _ in good], dtype=float)
    return float(np.sum(w * v) / np.sum(w)), float(math.sqrt(1.0 / np.sum(w))), len(good)


def write_sector_all_diagnostics(
    frame: pd.DataFrame,
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path,
) -> Path:
    ensure_directory(output_dir)
    rows = []
    for branch, label in (("e_phi", "electron"), ("p_phi", "pion")):
        sectors = fd_sector_from_phi_rad_array(events[branch])
        for sector in range(1, 7):
            subset_mask = sectors == sector
            subset = {name: np.asarray(values)[subset_mask] for name, values in events.items()}
            for period in PERIODS:
                for ti in range(len(MINUS_TPRIME_BINS_GEV2)):
                    estimates = []
                    n_events = 0
                    for xi in range(len(XB_BINS)):
                        b = xi * len(MINUS_TPRIME_BINS_GEV2) + ti + 1
                        fixed = _fixed_other_physics(frame, b)
                        try:
                            fit = fit_one_variant(
                                subset, run_states, dilution_records, b, "nominal",
                                active_periods=(period,), fixed_physics_parameters=fixed,
                            )
                            if fit["valid"]:
                                estimates.append((fit["values"]["ll0"], fit["errors"]["ll0"]))
                                n_events += int(fit["metadata"]["number_of_events"])
                            # endif
                        except Exception:
                            pass
                        # endtry
                    # endfor
                    value, error, n_xb = _inverse_variance_combine(estimates)
                    rows.append({
                        "particle": label, "sector": sector, "period": period,
                        "t_index": ti, "minus_tprime_low_gev2": MINUS_TPRIME_BINS_GEV2[ti][0],
                        "minus_tprime_high_gev2": MINUS_TPRIME_BINS_GEV2[ti][1],
                        "ll0": value, "ll0_stat": error, "xB_bins_combined": n_xb,
                        "number_of_events": n_events,
                    })
                # endfor
            # endfor
        # endfor
    # endfor
    out = output_dir / "all_sector_dependence_xB_integrated.csv"
    result = pd.DataFrame(rows)
    result.to_csv(out, index=False)
    for particle in ("electron", "pion"):
        for ti in range(len(MINUS_TPRIME_BINS_GEV2)):
            fig, ax = plt.subplots(figsize=(8.0, 5.0))
            for period in PERIODS:
                sub = result[(result["particle"] == particle) & (result["period"] == period) & (result["t_index"] == ti)]
                ax.errorbar(sub["sector"], sub["ll0"], yerr=sub["ll0_stat"], fmt="o-", label=period)
            # endfor
            ax.axhline(0.0, linewidth=1.0, linestyle="--")
            ax.set_xlabel(f"{particle.capitalize()} FD sector")
            ax.set_ylabel(r"$x_B$-integrated $A_{LL}$")
            ax.set_xticks(range(1, 7))
            ax.legend(frameon=False)
            fig.tight_layout()
            fig.savefig(output_dir / f"all_{particle}_sector_tbin_{ti + 1}.png", dpi=180)
            plt.close(fig)
        # endfor
    # endfor
    return out


def write_run_number_all_diagnostic(
    frame: pd.DataFrame,
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    dilution_records: Mapping[tuple[str, int], DilutionRecord],
    output_dir: Path,
    runs_per_block: int = 8,
) -> Path:
    ensure_directory(output_dir)
    rows = []
    for period in PERIODS:
        runs = np.asarray(run_states[period]["run"], dtype=int)
        pts = np.asarray(run_states[period]["pt"], dtype=float)
        # Never allow a block to cross a target-polarization sign reversal.
        start = 0
        blocks = []
        while start < runs.size:
            sign = int(np.sign(pts[start]))
            stop = start
            while stop < runs.size and int(np.sign(pts[stop])) == sign and stop - start < runs_per_block:
                stop += 1
            # endwhile
            blocks.append((runs[start:stop], sign))
            start = stop
        # endwhile
        for ib, (block_runs, sign) in enumerate(blocks, start=1):
            estimates = []
            n_events = 0
            ranges = tuple((int(run), int(run)) for run in block_runs)
            for b in range(1, NUMBER_OF_BINS + 1):
                fixed = _fixed_other_physics(frame, b)
                try:
                    fit = fit_one_variant(
                        events, run_states, dilution_records, b, "nominal",
                        active_periods=(period,), fixed_physics_parameters=fixed,
                        run_ranges_filter=ranges,
                    )
                    if fit["valid"]:
                        estimates.append((fit["values"]["ll0"], fit["errors"]["ll0"]))
                        n_events += int(fit["metadata"]["number_of_events"])
                    # endif
                except Exception:
                    pass
                # endtry
            # endfor
            value, error, n_bins = _inverse_variance_combine(estimates)
            rows.append({
                "period": period, "block": ib, "run_min": int(block_runs.min()),
                "run_max": int(block_runs.max()), "mean_run": float(np.mean(block_runs)),
                "target_sign": sign, "ll0": value, "ll0_stat": error,
                "bins_combined": n_bins, "number_of_events": n_events,
            })
        # endfor
    # endfor
    out = output_dir / "all_vs_run_number_blocks.csv"
    result = pd.DataFrame(rows)
    result.to_csv(out, index=False)
    if not result.empty:
        fig, ax = plt.subplots(figsize=(11.0, 5.5))
        for period in PERIODS:
            sub = result[result["period"] == period]
            ax.errorbar(sub["mean_run"], sub["ll0"], yerr=sub["ll0_stat"], fmt="o", label=period)
        # endfor
        ax.axhline(0.0, linewidth=1.0, linestyle="--")
        ax.set_xlabel("Run number")
        ax.set_ylabel(r"Integrated $A_{LL}$")
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(output_dir / "all_vs_run_number_blocks.png", dpi=180)
        plt.close(fig)
    # endif
    return out


def write_charge_normalized_kinematic_distributions(
    events: Mapping[str, np.ndarray],
    run_states: Mapping[str, Mapping[str, np.ndarray]],
    output_dir: Path,
) -> Path:
    ensure_directory(output_dir)
    variables = {
        "Q2": ("Q2", 30), "W": ("W", 30), "y": ("y", 30),
        "phi": ("phi", 30), "p_pi": ("p_p", 30), "theta_pi": ("p_theta", 30),
        "xB": ("xB", 30), "minus_tprime": ("minus_tprime", 30),
    }
    rows = []
    for b in range(1, NUMBER_OF_BINS + 1):
        bmask = events["bin_number"] == b
        if not np.any(bmask):
            continue
        # endif
        for label, (branch, nbins) in variables.items():
            all_values = np.asarray(events[branch][bmask], dtype=float)
            finite = all_values[np.isfinite(all_values)]
            if finite.size < 2:
                continue
            # endif
            low, high = np.quantile(finite, [0.005, 0.995])
            if label == "phi":
                low, high = 0.0, 2.0 * math.pi
            # endif
            edges = np.linspace(low, high, nbins + 1)
            fig, ax = plt.subplots(figsize=(7.5, 5.0))
            for period in PERIODS:
                state = run_states[period]
                for helicity in (-1, 1):
                    mask = bmask & (events["period_index"] == PERIOD_INDEX[period]) & (events["helicity"] == helicity)
                    vals = np.asarray(events[branch][mask], dtype=float)
                    charge = float(np.sum(state["q_plus"] if helicity > 0 else state["q_minus"]))
                    hist, _ = np.histogram(vals[np.isfinite(vals)], bins=edges)
                    rate = hist / charge if charge > 0.0 else np.full(hist.shape, np.nan)
                    centers = 0.5 * (edges[:-1] + edges[1:])
                    ax.step(centers, rate, where="mid", label=f"{period} h={helicity:+d}")
                    rows.append({
                        "bin_number": b, "variable": label, "period": period,
                        "helicity": helicity, "charge": charge,
                        "mean": float(np.mean(vals)) if vals.size else math.nan,
                        "std": float(np.std(vals, ddof=1)) if vals.size > 1 else math.nan,
                        "events": int(vals.size),
                    })
                # endfor
            # endfor
            ax.set_xlabel(label)
            ax.set_ylabel("Counts / accumulated charge")
            ax.legend(frameon=False, fontsize=8, ncol=2)
            fig.tight_layout()
            fig.savefig(output_dir / f"bin_{b:02d}_{label}_charge_normalized.png", dpi=160)
            plt.close(fig)
        # endfor
    # endfor
    out = output_dir / "kinematic_distribution_summary.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    return out

def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Fit nominal and systematic-variation RGC exclusive-pi+ "
            "structure-function ratios."
        )
    )
    parser.add_argument("--tree", default=DEFAULT_TREE_NAME)
    parser.add_argument("--input", action="append", default=[], type=parse_input_override, metavar="PERIOD=FILE")
    parser.add_argument("--isr-input", action="append", default=[], type=parse_input_override, metavar="PERIOD=FILE")
    parser.add_argument(
        "--uncorrected-input", action="append", default=[],
        type=parse_input_override, metavar="PERIOD=FILE",
        help=(
            "Override a no-momentum-corrections ROOT input. May be repeated "
            "for PERIOD=FILE."
        ),
    )
    parser.add_argument("--run-info-csv", type=Path, default=DEFAULT_RUN_INFO_CSV)
    parser.add_argument("--cut-json", type=Path, default=DEFAULT_CUT_JSON)
    parser.add_argument(
        "--polynomial-cut-json", type=Path, default=DEFAULT_POLYNOMIAL_CUT_JSON,
        help="Polynomial-only cut JSON produced by channel_selection_mx2_fits.py.",
    )
    parser.add_argument(
        "--fit-method-diagnostic", action="store_true",
        help=(
            "Run only the statistical comparison of the nominal carbon-assisted "
            "A(mu_C,sigma_C) extraction with centroid-only A(mu_poly,sigma_C), "
            "width-only A(mu_C,sigma_poly), and full polynomial "
            "A(mu_poly,sigma_poly) extractions. Writes five 24-bin polarized "
            "comparison canvases plus CSV/JSON tables and assigns no systematic."
        ),
    )
    parser.add_argument("--isr-cut-json", type=Path, default=None)
    parser.add_argument("--dilution-json", type=Path, default=None)
    parser.add_argument(
        "--isr-dilution-json",
        type=Path,
        default=None,
        help=(
            "Compact dilution-factor JSON recalculated from the matched "
            "internal+external-ISR samples. By default this is resolved from "
            "the ISR output of determine_dilution_factor.py."
        ),
    )
    parser.add_argument("--dilution-dir", type=Path, default=DEFAULT_DILUTION_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache", type=Path, default=None, help="Legacy nominal-cache override.")
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--disable-isr", action="store_true")
    parser.add_argument(
        "--disable-momentum-corrections", action="store_true",
        help="Skip the momentum-corrected versus uncorrected study.",
    )
    parser.add_argument(
        "--disable-channel-selection", action="store_true",
        help="Skip the tight/nominal/loose exclusivity-window study.",
    )
    parser.add_argument("--chunk-size", default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--workers", type=int, default=MAXIMUM_WORKERS)
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--appendix-mle-diagnostics", action="store_true",
        help=("Run only appendix-ready fit-quality diagnostics for the production "
              "seven-parameter joint MLE: 24 conditional phi data/model plots and "
              "profile-likelihood canvases for all seven amplitudes."))
    parser.add_argument("--appendix-profile-points", type=int, default=51,
        help="Fixed-parameter points per appendix profile likelihood (default: 17).")
    parser.add_argument("--appendix-phi-bins", type=int, default=18,
        help="Diagnostic phi bins for conditional data/model plots (default: 12).")
    parser.add_argument(
        "--period-stability-only", action="store_true",
        help=(
            "Run only the zero-UU (u1=u2=0) simultaneous MLE plus the three "
            "independent zero-UU period-only MLE extractions in each of the 24 "
            "bins, then make the five polarized run-period comparison/pull plots "
            "and write a table of every |pull|>2.5 excursion. No systematic or "
            "target-axis studies are run."
        ),
    )
    parser.add_argument(
        "--period-full-charge-scan", action="store_true",
        help=(
            "With --period-stability-only, additionally run the expensive full "
            "target-epoch and h*s charge-distortion scan. Off by default."
        ),
    )
    parser.add_argument(
        "--period-epoch-refits", action="store_true",
        help=(
            "With --period-stability-only, additionally run all target-polarization "
            "epoch A_LL refits. Off by default after the dedicated validation study."
        ),
    )
    parser.add_argument(
        "--period-flagged-diagnostics", action="store_true",
        help=(
            "With --period-stability-only, additionally rerun the expensive flagged-bin "
            "profiles/exclusivity-window diagnostics. Off by default once validated."
        ),
    )
    parser.add_argument(
        "--period-stability-plot-only", action="store_true",
        help=(
            "Regenerate the five published run-period stability PNGs directly "
            "from an existing period-stability structure_function_ratios.csv; "
            "do not read ROOT files or rerun any MLE fits."
        ),
    )
    parser.add_argument(
        "--period-stability-csv", type=Path, default=None,
        help=(
            "CSV to use with --period-stability-plot-only. By default uses "
            "OUTPUT_DIR/period_stability/tables/structure_function_ratios.csv."
        ),
    )
    parser.add_argument(
        "--period-stability-diagnostics", action="store_true",
        help=(
            "Run targeted nominal period-only diagnostics for the largest "
            "run-period excursions using the existing period-stability cache. "
            "Writes 2x2 raw phi-state yields and correlation matrices (Su22, Fa22, "
            "Sp23, combined), period and combined profile likelihoods, and a "
            "diagnostic summary table; no systematics."
        ),
    )
    parser.add_argument(
        "--diagnostic-bins", type=parse_bin_list,
        default=PERIOD_STABILITY_DIAGNOSTIC_BINS,
        help="Comma-separated bins for --period-stability-diagnostics "
             "(default: 4,6,7,14,17,18,19,21).",
    )
    parser.add_argument(
        "--diagnostic-profile-points", type=int, default=21,
        help="Number of fixed-parameter points per period profile (default: 21).",
    )
    parser.add_argument(
        "--baseline-zero-uu-only", action="store_true",
        help="Run only the baseline nominal likelihood with u1=u2=0 fixed; skip all systematic studies.",
    )
    parser.add_argument(
        "--clas6-cross-check", action="store_true",
        help=(
            "Run one statistical-only CLAS6 EG1b exclusive-pi+ cross-check. "
            "Applies 2<W<2.6 GeV and 1<Q2<5 GeV^2 to both datasets, "
            "re-fits RGC events in the common 24-bin scheme, and fits the "
            "tabulated CLAS6 A_UL/A_LL phi dependences without systematics."
        ),
    )
    parser.add_argument(
        "--clas6-data", type=Path,
        default=Path(__file__).resolve().parent / "exclpip.txt",
        help="Path to the CLAS6 EG1b exclpip.txt table (default: beside this script).",
    )
    parser.add_argument(
        "--double-spin-target-split", action="store_true",
        help=(
            "Run only the nominal fixed-target-orientation double-spin stability "
            "test. Fits P_t>0 and P_t<0 samples independently in all 24 bins "
            "with the full seven-term nominal likelihood and compares A_LL and "
            "A_LL^cos(phi). Uses the existing nominal selected-event cache; no "
            "systematic or alternative-method studies are run."
        ),
    )
    parser.add_argument(
        "--solenoid-split", action="store_true",
        help=(
            "Run only the solenoid-polarity stability diagnostic. Integrates "
            "over xB, fits all seven amplitudes independently for solenoid -1 "
            "and +1 in the six -tprime bins, and uses the existing nominal cache. "
            "No systematic or alternative-method studies are run."
        ),
    )
    parser.add_argument(
        "--rga-cross-check", action="store_true",
        help=(
            "Run only the Diehl et al. RGA exclusive-pi+ cross-check. "
            "Uses the published Q2-xB polygons and -t bins, imposes "
            "Q2>1.5 GeV^2, W>2 GeV, y<0.75, and fits only "
            "(F_UU^cos(phi), F_UU^cos(2phi), F_LU^sin(phi))/F_UU "
            "directly in each Diehl bin. Production extraction is unchanged."
        ),
    )
    return parser


def main() -> int:
    args = build_argument_parser().parse_args()
    workers = max(
        1,
        min(int(args.workers), MAXIMUM_WORKERS, 8, os.cpu_count() or 1, NUMBER_OF_BINS),
    )
    if args.clas6_cross_check:
        return run_clas6_cross_check(args, workers)
    # endif
    if args.rga_cross_check:
        return run_rga_cross_check(args)
    # endif
    root = args.output_dir.expanduser().resolve()
    if args.period_stability_diagnostics:
        return run_period_stability_diagnostics(args, root, workers)
    # endif
    if args.period_stability_plot_only:
        csv_path = (
            args.period_stability_csv.expanduser().resolve()
            if args.period_stability_csv is not None
            else root / "period_stability/tables/structure_function_ratios.csv"
        )
        if not csv_path.is_file():
            raise FileNotFoundError(
                f"Period-stability CSV not found: {csv_path}"
            )
        # endif
        frame = pd.read_csv(csv_path)
        stability_dir = root / "period_stability/plots/period_stability_published"
        paths = plot_period_stability_published(frame, stability_dir)
        print("[period-stability-plot-only] complete", flush=True)
        print(f"  Input:  {csv_path}", flush=True)
        print(f"  Plots:  {stability_dir}", flush=True)
        print(f"  Wrote:  {len(paths)} polarized stability PNGs", flush=True)
        return 0
    # endif
    if args.baseline_zero_uu_only or args.period_stability_only:
        args.disable_isr = True
        args.disable_momentum_corrections = True
        args.disable_channel_selection = True
    # endif
    nominal_dir = root / (
        "period_stability" if args.period_stability_only
        else ("nominal_zero_uu" if args.baseline_zero_uu_only else "nominal")
    )
    isr_dir = root / "isr"
    momentum_dir = root / "momentum_corrections"
    channel_dir = root / "channel_selection"
    diagnostics_dir = root / "diagnostics"
    for directory in (
        root, nominal_dir, isr_dir, momentum_dir, channel_dir, diagnostics_dir
    ):
        ensure_directory(directory)
    # endfor

    nominal_inputs = {period: Path(value) for period, value in DEFAULT_INPUTS.items()}
    for period, path in args.input:
        nominal_inputs[period] = path
    # endfor
    nominal_dilution = (
        args.dilution_json.expanduser().resolve()
        if args.dilution_json
        else find_default_dilution_json(
            args.dilution_dir.expanduser().resolve()
        ).resolve()
    )
    if args.appendix_mle_diagnostics:
        return run_appendix_mle_diagnostics(args, root, workers, nominal_dilution)
    # endif
    if args.double_spin_target_split:
        return run_double_spin_target_split_diagnostic(
            cache_path=(
                args.cache.expanduser().resolve()
                if args.cache
                else root / "nominal/cache/selected_events.npz"
            ),
            run_info_path=args.run_info_csv.expanduser().resolve(),
            dilution_json_path=nominal_dilution,
            output_dir=root,
            skip_plots=args.skip_plots,
        )
    # endif
    if args.solenoid_split:
        return run_solenoid_split_diagnostic(
            cache_path=(
                args.cache.expanduser().resolve()
                if args.cache
                else root / "nominal/cache/selected_events.npz"
            ),
            run_info_path=args.run_info_csv.expanduser().resolve(),
            dilution_json_path=nominal_dilution,
            output_dir=root,
            workers=workers,
            skip_plots=args.skip_plots,
        )
    # endif

    isr_dilution = None
    if not args.disable_isr and not args.fit_method_diagnostic:
        isr_dilution = (
            args.isr_dilution_json.expanduser().resolve()
            if args.isr_dilution_json
            else find_isr_dilution_json(
                args.dilution_dir.expanduser().resolve()
            ).resolve()
        )
    # endif

    nominal_cache = (
        args.cache.expanduser().resolve()
        if args.cache
        else nominal_dir / "cache/selected_events.npz"
    )

    channel_loose_cache = channel_dir / "loose/cache/selected_events.npz"

    print("=" * 78, flush=True)
    print("RGC enpi+ structure-function-ratio extraction: startup", flush=True)
    print("=" * 78, flush=True)
    print(f"Output root:          {root}", flush=True)
    print(f"Workers:              {workers} (maximum {MAXIMUM_WORKERS})", flush=True)
    print(f"Reuse cache:          {args.reuse_cache}", flush=True)
    print(f"Skip plots:           {args.skip_plots}", flush=True)
    print(f"ISR study enabled:    {not args.disable_isr and not args.fit_method_diagnostic}", flush=True)
    print(
        f"Momentum study:       {not args.disable_momentum_corrections and not args.fit_method_diagnostic}",
        flush=True,
    )
    print(
        f"Channel-selection:    {not args.disable_channel_selection and not args.fit_method_diagnostic}",
        flush=True,
    )
    print(f"Nominal cache:        {nominal_cache}", flush=True)
    if not args.disable_channel_selection and not args.fit_method_diagnostic:
        print(f"Loose superset cache: {channel_loose_cache}", flush=True)
        if args.reuse_cache:
            print(
                "[startup] --reuse-cache: existing loose cache will be loaded; "
                "no ROOT cache rebuild will be performed.",
                flush=True,
            )
        else:
            print(
                "[startup] Clean run: building the loose (3 sigma) superset "
                "cache FIRST. The nominal (2 sigma) cache will then be derived "
                "from it without rereading the nominal ROOT trees.",
                flush=True,
            )
        # endif
    else:
        print(
            "[startup] Channel-selection study disabled: nominal cache will be "
            "built directly from the nominal ROOT inputs.",
            flush=True,
        )
    # endif
    print("=" * 78, flush=True)

    nominal_source_cache: Path | None = None
    channel_loose_cache_ready = False
    if not args.disable_channel_selection and not args.fit_method_diagnostic:
        if args.reuse_cache:
            print(
                f"[startup/cache] checking loose superset cache: "
                f"{channel_loose_cache}",
                flush=True,
            )
            if not channel_loose_cache.is_file():
                raise FileNotFoundError(
                    "--reuse-cache was requested, but the loose channel-selection "
                    f"cache is missing: {channel_loose_cache}"
                )
            # endif
            print(
                "[startup/cache] loose superset cache exists; it will be reused.",
                flush=True,
            )
            channel_loose_cache_ready = True
        else:
            print(
                "[startup/cache] preparing loose (3 sigma) superset cache "
                "before nominal fitting...",
                flush=True,
            )
            run_records_for_cache = parse_run_info_csv(
                args.run_info_csv.expanduser().resolve()
            )
            loose_cuts_for_cache = load_channel_cuts(
                args.cut_json.expanduser().resolve(), cut_label="loose"
            )
            build_event_cache(
                input_paths=nominal_inputs,
                tree_name=args.tree,
                chunk_size=args.chunk_size,
                run_records=run_records_for_cache,
                cuts=loose_cuts_for_cache,
                cache_path=channel_loose_cache,
            )
            channel_loose_cache_ready = True
            print(
                "[startup/cache] loose superset cache is ready; proceeding to "
                "the nominal analysis.",
                flush=True,
            )
        # endif
        nominal_source_cache = channel_loose_cache
    # endif

    nominal_result = run_analysis_variant(
        sample_variant="nominal",
        input_paths=nominal_inputs,
        run_info_path=args.run_info_csv.expanduser().resolve(),
        cut_json_path=args.cut_json.expanduser().resolve(),
        dilution_json_path=nominal_dilution,
        output_dir=nominal_dir,
        cache_path=nominal_cache,
        tree_name=args.tree,
        chunk_size=args.chunk_size,
        workers=workers,
        reuse_cache=args.reuse_cache,
        skip_plots=args.skip_plots,
        include_target_axis_study=(
            not args.baseline_zero_uu_only and not args.period_stability_only
            and not args.fit_method_diagnostic
        ),
        include_period_diagnostics=(
            not args.baseline_zero_uu_only and not args.fit_method_diagnostic
        ),
        cut_label="nominal",
        source_cache_path=(None if args.reuse_cache else nominal_source_cache),
        zero_uu_baseline=(args.baseline_zero_uu_only or args.period_stability_only),
    )

    if args.period_stability_only:
        stability_dir = nominal_dir / "plots/period_stability_published"
        paths = []
        if not args.skip_plots:
            paths = plot_period_stability_published(
                nominal_result["frame"], stability_dir
            )
        # endif
        excursion_path = nominal_dir / "tables/period_stability_excursions_gt_2p5sigma.csv"
        excursions = write_period_stability_excursion_table(
            nominal_result["frame"], excursion_path, threshold=2.5
        )
        print("[period-stability-only] complete (u1=u2=0)", flush=True)
        print(f"  Results:    {nominal_result['csv']}", flush=True)
        print(f"  Plots:      {stability_dir}", flush=True)
        print(f"  Excursions: {excursion_path} ({len(excursions)} rows)", flush=True)
        if paths:
            print(f"  Wrote:   {len(paths)} polarized stability PNGs", flush=True)
        # endif

        # Global/structured consistency diagnostics are always written, even if
        # no single residual exceeds the 2.5-sigma inspection threshold.
        stability_run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
        stability_run_states = run_state_arrays(stability_run_records)
        stability_cuts = load_channel_cuts(args.cut_json.expanduser().resolve(), cut_label="nominal")
        global_products = write_period_stability_global_diagnostics(
            nominal_result["frame"], nominal_dir / "diagnostics" / "global_consistency",
            stability_cuts, stability_run_states,
        )
        print(f"  Global LRT: {global_products['global_lrt']}", flush=True)
        print(f"  Observable Wald tests: {global_products['observable_wald']}", flush=True)

        # Exact observable-by-observable profiled likelihood-ratio tests.
        # For each observable, only that amplitude is constrained to be common
        # across periods under H0; the other four polarized amplitudes remain
        # period-specific and are profiled in both hypotheses.
        observable_lrt_path = write_period_stability_observable_likelihood_ratio_tests(
            nominal_result["results"], load_event_cache(nominal_cache),
            stability_run_states,
            load_dilution_factors(nominal_dilution, cut_label="nominal"),
            nominal_dir / "diagnostics" / "global_consistency",
        )
        print(f"  Observable exact LRTs: {observable_lrt_path}", flush=True)
        print(f"  Period scale tests: {global_products['scale_tests']}", flush=True)

        # Audit the actual run/helicity charge bookkeeping and measure the
        # response to a controlled h*sign(Pt)-correlated charge distortion.
        # This directly tests whether a small relative Faraday-cup/exposure
        # error could preferentially generate the constant A_LL period pattern.
        stability_events = load_event_cache(nominal_cache)
        stability_dilutions = load_dilution_factors(nominal_dilution, cut_label="nominal")
        charge_products = write_charge_normalization_diagnostics(
            frame=nominal_result["frame"], events=stability_events,
            run_states=stability_run_states, dilution_records=stability_dilutions,
            output_dir=nominal_dir / "diagnostics" / "charge_normalization",
            run_full_scan=args.period_full_charge_scan,
        )
        print(f"  Charge audit: {charge_products['period_audit']}", flush=True)
        if args.period_full_charge_scan:
            print(f"  Full charge response: {charge_products['response']}", flush=True)
            print(f"  Full charge summary: {charge_products['summary']}", flush=True)
        # endif

        target_charge_ratio = write_simultaneous_target_charge_ratio_fit(
            frame=nominal_result["frame"], events=stability_events,
            run_states=stability_run_states, dilution_records=stability_dilutions,
            output_dir=nominal_dir / "diagnostics" / "charge_normalization",
        )
        print(f"  Target-charge-ratio fit: {target_charge_ratio['solution']}", flush=True)
        print(f"  Target-charge-ratio priors: {target_charge_ratio['priors']}", flush=True)

        beam_products = write_beam_polarization_response_diagnostic(
            frame=nominal_result["frame"], events=stability_events,
            run_states=stability_run_states, dilution_records=stability_dilutions,
            output_dir=nominal_dir / "diagnostics" / "beam_polarization",
        )
        print(f"  Beam-polarization response: {beam_products['summary']}", flush=True)

        # Extended final run-period diagnostics: inclusive e pi+ X control,
        # independent DIS PbPt, detector-sector dependence, chronological
        # stability, and charge-normalized kinematic population comparisons.
        extended_dir = nominal_dir / "diagnostics" / "extended_period_checks"
        diagnostic_cache_path = nominal_dir / "cache" / "period_stability_no_mx2_events.npz"
        diagnostic_events = build_period_stability_diagnostic_cache(
            nominal_inputs, args.tree, args.chunk_size, stability_run_records,
            diagnostic_cache_path,
        )
        inclusive_path = write_inclusive_epi_control_diagnostic(
            nominal_result["frame"], diagnostic_events, stability_run_states,
            extended_dir / "inclusive_epiX",
        )
        print(f"  Inclusive e pi+ X A_raw/(PbPt), no dilution correction: {inclusive_path}", flush=True)

        dis_pbpt = load_dis_pbpt_spreadsheets(Path(__file__).resolve().parent)
        dis_pbpt_path = write_dis_pbpt_all_diagnostic(
            nominal_result["frame"], stability_events, stability_run_states,
            stability_dilutions, dis_pbpt, extended_dir / "dis_pbpt",
        )
        print(f"  DIS PbPt A_LL stability: {dis_pbpt_path}", flush=True)

        nominal_diagnostic_events = select_nominal_mx2_from_diagnostic_cache(
            diagnostic_events, stability_cuts
        )
        sector_path = write_sector_all_diagnostics(
            nominal_result["frame"], nominal_diagnostic_events, stability_run_states,
            stability_dilutions, extended_dir / "sectors",
        )
        print(f"  Sector A_LL dependence: {sector_path}", flush=True)

        run_path = write_run_number_all_diagnostic(
            nominal_result["frame"], stability_events, stability_run_states,
            stability_dilutions, extended_dir / "run_number",
        )
        print(f"  Run-number A_LL dependence: {run_path}", flush=True)

        distribution_path = write_charge_normalized_kinematic_distributions(
            nominal_diagnostic_events, stability_run_states,
            extended_dir / "kinematic_distributions",
        )
        print(f"  Charge-normalized kinematics: {distribution_path}", flush=True)

        if args.period_epoch_refits:
            epoch_refits = write_flagged_epoch_refits(
                frame=nominal_result["frame"], events=stability_events,
                run_states=stability_run_states, dilution_records=stability_dilutions,
                output_dir=nominal_dir / "diagnostics" / "run_epochs", threshold=2.5,
            )
            print(f"  Target-epoch refits: {epoch_refits}", flush=True)
        # endif

        if args.period_flagged_diagnostics and not excursions.empty and not args.skip_plots:
            flagged_bins = tuple(sorted(set(excursions["bin_number"].astype(int))))
            flagged_map = {int(b): tuple(sorted(set(g["parameter"].astype(str)))) for b, g in excursions.groupby("bin_number")}
            args.diagnostic_bins = flagged_bins
            args._period_stability_flagged_parameters = flagged_map
            print(f"[period-stability-only] automatically diagnosing flagged bins {flagged_bins}", flush=True)
            run_period_stability_diagnostics(args, root, workers)
        # endif
        return 0
    # endif

    if args.baseline_zero_uu_only:
        run_xb_integrated_tprime_zero_uu_study(
            cache_path=nominal_cache,
            run_info_path=args.run_info_csv.expanduser().resolve(),
            dilution_json_path=nominal_dilution,
            output_dir=nominal_dir,
        )
        print(f"[baseline-zero-uu] complete: {nominal_dir}", flush=True)
        return 0
    # endif

    if args.fit_method_diagnostic:
        polynomial_cut_json = args.polynomial_cut_json.expanduser().resolve()
        if not polynomial_cut_json.is_file():
            raise FileNotFoundError(
                "Fit-method diagnostic requires the polynomial-only cut JSON: "
                f"{polynomial_cut_json}. Run channel_selection_mx2_fits.py first."
            )
        # endif
        diagnostic_dir = root / "fit_method_diagnostic"
        generated_cuts_dir = diagnostic_dir / "generated_cuts"
        centroid_cut_json, width_cut_json = write_hybrid_fit_method_cut_jsons(
            carbon_path=args.cut_json.expanduser().resolve(),
            polynomial_path=polynomial_cut_json,
            output_dir=generated_cuts_dir,
        )

        variants = {
            "centroid_only_extraction": ("fit_method_centroid_only", centroid_cut_json),
            "width_only_extraction": ("fit_method_width_only", width_cut_json),
            "polynomial_extraction": ("fit_method_polynomial", polynomial_cut_json),
        }
        diagnostic_results: dict[str, Any] = {}
        for dirname, (sample_variant, cut_path) in variants.items():
            variant_dir = diagnostic_dir / dirname
            diagnostic_results[dirname] = run_analysis_variant(
                sample_variant=sample_variant,
                input_paths=nominal_inputs,
                run_info_path=args.run_info_csv.expanduser().resolve(),
                cut_json_path=cut_path,
                dilution_json_path=nominal_dilution,
                output_dir=variant_dir,
                cache_path=variant_dir / "cache/selected_events.npz",
                tree_name=args.tree,
                chunk_size=args.chunk_size,
                workers=workers,
                reuse_cache=args.reuse_cache,
                skip_plots=True,
                include_target_axis_study=False,
                include_period_diagnostics=False,
                cut_label="nominal",
            )
        # endfor
        products = write_fit_method_diagnostic_products(
            nominal=nominal_result["frame"],
            centroid_only=diagnostic_results["centroid_only_extraction"]["frame"],
            width_only=diagnostic_results["width_only_extraction"]["frame"],
            polynomial=diagnostic_results["polynomial_extraction"]["frame"],
            output_dir=diagnostic_dir,
            carbon_cut_json=args.cut_json.expanduser().resolve(),
            polynomial_cut_json=polynomial_cut_json,
        )
        print("[fit-method-diagnostic] complete", flush=True)
        print(f"  CSV:   {products['csv']}", flush=True)
        print(f"  Plots: {diagnostic_dir / 'plots'}", flush=True)
        return 0
    # endif

    target_axis_study = write_target_axis_study_products(
        nominal_result["frame"], diagnostics_dir / "target_axis"
    )

    # Collaboration-facing estimator cross-check: compare the same five
    # polarized amplitudes with u1=u2 fixed to zero in an unbinned MLE and a
    # conventional 12-bin-in-phi chi2 fit. This diagnostic is nominal-only and
    # does not participate in any systematic or alternative-method studies.
    run_mle_vs_binned_chi2_diagnostic(
        cache_path=nominal_cache,
        run_info_path=args.run_info_csv.expanduser().resolve(),
        dilution_json_path=nominal_dilution,
        output_dir=root,
    )

    isr_result = None
    isr_comparison = None
    if not args.disable_isr:
        isr_inputs = {period: Path(value) for period, value in DEFAULT_ISR_INPUTS.items()}
        for period, path in args.isr_input:
            isr_inputs[period] = path
        # endfor
        isr_cut = resolve_isr_cut_json(args.isr_cut_json)
        isr_result = run_analysis_variant(
            sample_variant="isr",
            input_paths=isr_inputs,
            run_info_path=args.run_info_csv.expanduser().resolve(),
            cut_json_path=isr_cut,
            dilution_json_path=isr_dilution,
            output_dir=isr_dir,
            cache_path=isr_dir / "cache/selected_events.npz",
            tree_name=args.tree,
            chunk_size=args.chunk_size,
            workers=workers,
            reuse_cache=args.reuse_cache,
            skip_plots=args.skip_plots,
            include_target_axis_study=False,
            cut_label="nominal",
        )
        isr_comparison = write_nominal_isr_comparison_products(
            nominal_result["frame"], isr_result["frame"], diagnostics_dir / "isr"
        )
    # endif

    momentum_result = None
    momentum_comparison = None
    if not args.disable_momentum_corrections:
        uncorrected_inputs = {
            period: Path(value)
            for period, value in DEFAULT_UNCORRECTED_INPUTS.items()
        }
        for period, path in args.uncorrected_input:
            uncorrected_inputs[period] = path
        # endfor
        uncorrected_dir = momentum_dir / "uncorrected"
        momentum_result = run_analysis_variant(
            sample_variant="momentum_uncorrected",
            input_paths=uncorrected_inputs,
            run_info_path=args.run_info_csv.expanduser().resolve(),
            cut_json_path=args.cut_json.expanduser().resolve(),
            dilution_json_path=nominal_dilution,
            output_dir=uncorrected_dir,
            cache_path=uncorrected_dir / "cache/selected_events.npz",
            tree_name=args.tree,
            chunk_size=args.chunk_size,
            workers=workers,
            reuse_cache=args.reuse_cache,
            skip_plots=args.skip_plots,
            include_target_axis_study=False,
            cut_label="nominal",
        )
        momentum_comparison = write_momentum_correction_comparison_products(
            corrected=nominal_result["frame"],
            uncorrected=momentum_result["frame"],
            output_dir=diagnostics_dir / "momentum_corrections",
        )
    # endif

    channel_results: dict[str, Any] = {}
    channel_comparison = None
    if not args.disable_channel_selection:
        loose_dir = channel_dir / "loose"
        nominal_window_dir = channel_dir / "nominal"
        tight_dir = channel_dir / "tight"

        # The loose cache was built once before the production nominal fit. It
        # is a strict superset, so all narrower selections are derived without
        # any additional ROOT pass.
        loose_result = run_analysis_variant(
            sample_variant="channel_selection_loose",
            input_paths=nominal_inputs,
            run_info_path=args.run_info_csv.expanduser().resolve(),
            cut_json_path=args.cut_json.expanduser().resolve(),
            dilution_json_path=nominal_dilution,
            output_dir=loose_dir,
            cache_path=channel_loose_cache,
            tree_name=args.tree,
            chunk_size=args.chunk_size,
            workers=workers,
            reuse_cache=channel_loose_cache_ready,
            skip_plots=args.skip_plots,
            include_target_axis_study=False,
            cut_label="loose",
        )
        # The production nominal result already uses the same nominal ROOT
        # inputs, nominal dilution factors, and nominal (2 sigma) exclusivity
        # selection.  Reuse that fitted result rather than repeating the same
        # simultaneous likelihood fit a second time.  This is an exact reuse,
        # not a change to the extraction.
        print(
            "[channel_selection] reusing production nominal (2 sigma) fit "
            "for the channel-selection comparison.",
            flush=True,
        )
        nominal_window_result = nominal_result

        tight_result = run_analysis_variant(
            sample_variant="channel_selection_tight",
            input_paths=nominal_inputs,
            run_info_path=args.run_info_csv.expanduser().resolve(),
            cut_json_path=args.cut_json.expanduser().resolve(),
            dilution_json_path=nominal_dilution,
            output_dir=tight_dir,
            cache_path=tight_dir / "cache/selected_events.npz",
            tree_name=args.tree,
            chunk_size=args.chunk_size,
            workers=workers,
            reuse_cache=args.reuse_cache,
            skip_plots=args.skip_plots,
            include_target_axis_study=False,
            cut_label="tight",
            source_cache_path=channel_loose_cache,
        )
        channel_results = {
            "tight": tight_result,
            "nominal": nominal_window_result,
            "loose": loose_result,
        }
        channel_comparison = write_channel_selection_comparison_products(
            tight=tight_result["frame"],
            nominal=nominal_window_result["frame"],
            loose=loose_result["frame"],
            output_dir=diagnostics_dir / "channel_selection",
        )
    # endif

    migration_products = write_migration_systematic_products(
        nominal_result["frame"], diagnostics_dir / "bin_migration"
    )

    nominal_with_systematics = attach_point_to_point_systematics(
        nominal=nominal_result["frame"],
        isr=(isr_result["frame"] if isr_result is not None else None),
        momentum_uncorrected=(
            momentum_result["frame"] if momentum_result is not None else None
        ),
        channel_tight=(
            channel_results["tight"]["frame"] if channel_results else None
        ),
        channel_loose=(
            channel_results["loose"]["frame"] if channel_results else None
        ),
    )
    nominal_result["frame"] = nominal_with_systematics
    point_to_point_csv = diagnostics_dir / "point_to_point_systematics.csv"
    nominal_with_systematics.to_csv(point_to_point_csv, index=False)
    # Also update the production CSV so downstream consumers receive the
    # source-by-source and total point-to-point systematic columns directly.
    nominal_with_systematics.to_csv(Path(nominal_result["csv"]), index=False)
    # Rewrite the publication-facing nominal LaTeX table now that the final
    # point-to-point systematic has been assembled.
    write_latex_table(
        nominal_with_systematics,
        Path(nominal_result["latex"]),
        include_target_axis_uncertainty=False,
        sample_variant="nominal",
    )

    if not args.skip_plots:
        plot_parameter_summaries(
            nominal_with_systematics,
            nominal_dir / "plots/all_bins",
            include_target_axis_uncertainty=False,
        )
        plot_aggregated_by_x(
            nominal_with_systematics,
            nominal_dir / "plots/aggregated",
            include_target_axis_uncertainty=False,
        )
    # endif

    all_barlow_records: list[dict[str, Any]] = list(
        target_axis_study["barlow_summary"]
    )
    if isr_comparison is not None:
        all_barlow_records.extend(isr_comparison["barlow_summary"])
    # endif
    if momentum_comparison is not None:
        all_barlow_records.extend(momentum_comparison["barlow_summary"])
    # endif
    if channel_comparison is not None:
        all_barlow_records.extend(channel_comparison["barlow_summary"])
    # endif
    # Publication-level Barlow summaries exclude the two fitted UU nuisance
    # modulations.  Their detailed diagnostic files remain available in each
    # source-specific directory.
    published_barlow_records = [
        record for record in all_barlow_records
        if record.get("parameter") in PUBLISHED_SYSTEMATIC_PARAMETERS
    ]
    barlow_summary_products = write_barlow_summary_products(
        published_barlow_records, diagnostics_dir / "barlow"
    )

    manifest_path = root / "asymmetry_extraction_manifest.json"
    write_json(manifest_path, {
        "schema_version": 7,
        "production_policy": (
            "Nominal results are production. ISR/external-ISR results are a "
            "separate radiation diagnostic. Momentum corrections are evaluated "
            "by comparing the nominal corrected extraction with matched "
            "uncorrected ROOT files. Channel selection is evaluated "
            "with matched tight, nominal, and loose dilution factors. The "
            "recommended per-bin channel-selection uncertainty is the RMS of "
            "the tight-minus-nominal and loose-minus-nominal shifts, while the "
            "complete variation vectors are retained as coherent alternatives. "
            "The momentum-correction on/off comparison is diagnostic only and is "
            "not assigned as a systematic. For the five published polarized ratios, "
            "radiation, channel-selection RMS, and 1.5-scaled bin-migration uncertainties "
            "are combined in quadrature and drawn as bars from y=0. Target-axis "
            "variants are retained as interpretation studies and are not assigned "
            "as systematic uncertainties. "
            "The two fitted UU modulations are excluded from publication-level "
            "systematic summaries."
        ),
        "nominal": {
            key: value for key, value in nominal_result.items()
            if key not in ("frame", "results")
        },
        "isr": (
            {
                key: value for key, value in isr_result.items()
                if key not in ("frame", "results")
            }
            if isr_result else None
        ),
        "target_axis_study": target_axis_study,
        "nominal_isr_comparison": isr_comparison,
        "momentum_corrections": (
            {
                key: value for key, value in momentum_result.items()
                if key not in ("frame", "results")
            }
            if momentum_result else None
        ),
        "momentum_correction_comparison": momentum_comparison,
        "momentum_correction_policy": (
            "Diagnostic only; abs(uncorrected-corrected) is not included in the "
            "assigned point-to-point systematic."
        ),
        "bin_migration": migration_products,
        "channel_selection": {
            label: {
                key: value for key, value in result.items()
                if key not in ("frame", "results")
            }
            for label, result in channel_results.items()
        } if channel_results else None,
        "channel_selection_comparison": channel_comparison,
        "point_to_point_systematics_csv": str(point_to_point_csv),
        "point_to_point_systematic_definition": (
            "sqrt(radiation^2 + channel_selection_RMS^2 + bin_migration^2), "
            "for published polarized observables; target-axis variants are "
            "interpretation studies only"
        ),
        "barlow_summary": barlow_summary_products,
    })

    print("\nStructure-function-ratio study complete.")
    print(f"  Nominal:           {nominal_dir}")
    if isr_result:
        print(f"  ISR:               {isr_dir}")
        print(f"  ISR diagnostics:   {diagnostics_dir / 'isr'}")
    # endif
    if momentum_result:
        print(f"  Momentum study:    {momentum_dir}")
        print(
            f"  Momentum diagnostics: "
            f"{diagnostics_dir / 'momentum_corrections'}"
        )
    # endif
    if channel_results:
        print(f"  Channel selection: {channel_dir}")
        print(
            f"  Channel diagnostics: {diagnostics_dir / 'channel_selection'}"
        )
    # endif
    print(f"  Diagnostics:       {diagnostics_dir}")
    print(f"  Barlow readout:    {barlow_summary_products['text']}")
    print(f"  Manifest:          {manifest_path}")

    total_pass = sum(record["number_pass"] for record in published_barlow_records)
    total_fail = sum(record["number_fail"] for record in published_barlow_records)
    total_undefined = sum(
        record["number_undefined"] for record in published_barlow_records
    )
    total_defined = total_pass + total_fail
    total_pass_percent = (
        100.0 * total_pass / total_defined if total_defined > 0 else float("nan")
    )
    total_fail_percent = (
        100.0 * total_fail / total_defined if total_defined > 0 else float("nan")
    )
    print(
        "  Barlow bins:       "
        f"{total_pass}/{total_defined} defined bins pass "
        f"({total_pass_percent:.1f}%); "
        f"{total_fail}/{total_defined} fail "
        f"({total_fail_percent:.1f}%); "
        f"undefined={total_undefined}"
    )
    for record in barlow_summary_products["effect_summary"]:
        defined = int(record["number_defined"])
        pass_fraction = record["pass_fraction_defined"]
        fail_fraction = record["fail_fraction_defined"]
        pass_percent = 100.0 * pass_fraction if pass_fraction is not None else float("nan")
        fail_percent = 100.0 * fail_fraction if fail_fraction is not None else float("nan")
        print(
            f"    {record['effect']} / {record['variation']}: "
            f"{record['number_pass']}/{defined} defined bins pass "
            f"({pass_percent:.1f}%); "
            f"{record['number_fail']}/{defined} fail "
            f"({fail_percent:.1f}%); "
            f"undefined={record['number_undefined']}"
        )
    # endfor

    invalid = nominal_result["invalid_bins"]
    if isr_result is not None:
        invalid += isr_result["invalid_bins"]
    # endif
    if momentum_result is not None:
        invalid += momentum_result["invalid_bins"]
    # endif
    for result in channel_results.values():
        invalid += result["invalid_bins"]
    # endfor
    return 2 if invalid else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
