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

PERIODS: tuple[str, ...] = ("su22", "fa22", "sp23")
PERIOD_LABELS: dict[str, str] = {
    "su22": "Su22",
    "fa22": "Fa22",
    "sp23": "Sp23",
}

PERIOD_COLORS: dict[str, str] = {
    "su22": "tab:orange",
    "fa22": "tab:green",
    "sp23": "tab:red",
}
COMBINED_COLOR = "tab:blue"

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

def load_channel_cuts(
    path: Path,
    cut_label: str = "nominal",
) -> dict[tuple[str, int], CutRecord]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing channel-selection cut JSON: {path}")
    # endif

    payload = json.loads(path.read_text(encoding="utf-8"))
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
                base = (
                    finite
                    & active_run
                    & np.isin(helicity, (-1, 1))
                    & (bin_number >= 1)
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
):
    mask = events["bin_number"] == bin_number
    if len(active_periods) != len(PERIODS):
        allowed_indices = np.asarray(
            [PERIOD_INDEX[period] for period in active_periods],
            dtype=np.int8,
        )
        mask &= np.isin(events["period_index"], allowed_indices)
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
            "observed_charge": np.where(
                event_h > 0.0,
                state_q_plus[observed_state_index],
                state_q_minus[observed_state_index],
            ),
            "charge_sum": float(np.sum(state_q_plus + state_q_minus)),
            "helicity_charge_sum": float(
                np.sum(state_q_plus - state_q_minus)
            ),
            "target_charge_sum": float(
                np.sum(state_pt * (state_q_plus + state_q_minus))
            ),
            "helicity_target_charge_sum": float(
                np.sum(state_pt * (state_q_plus - state_q_minus))
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
                beam_polarization=BEAM_POLARIZATION[period],
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
                BEAM_POLARIZATION[period]
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
                BEAM_POLARIZATION[period]
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
                    + BEAM_POLARIZATION[period]
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
                + data["helicity_target_charge_sum"] * double_coefficient
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
) -> dict[str, Any]:
    nll, metadata = make_bin_nll(
        events,
        run_states,
        dilution_records,
        bin_number,
        variant,
        active_periods=active_periods,
        transverse_scales=transverse_scales,
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
        candidate.strategy = 1

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

    attempted: list[Minuit] = []
    for start_values in start_candidates:
        candidate = configured_minuit(start_values)
        candidate.migrad(ncall=50000)
        if not candidate.fmin.is_valid:
            candidate.simplex(ncall=50000)
            candidate.strategy = 2
            candidate.migrad(ncall=120000)
        # endif
        candidate.hesse()
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
    nll, metadata = make_bin_nll(
        _WORKER_EVENTS, _WORKER_RUN_STATES, _WORKER_DILUTION_RECORDS,
        bin_number, "nominal", active_periods=(period,),
    )
    values = dict(PARAMETER_INITIAL_VALUES)
    values.update({name: float(nominal["values"][name])
                   for name in PHYSICS_PARAMETERS})
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
    return "|".join(str(task.get(k, "")) for k in
                    ("kind", "bin_number", "period", "constraint"))

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
    elif kind in {"photon_axis_projection", "external_data_informed"}:
        fit = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            kind,
            initial_values=nominal["values"],
        )
    elif kind == "period_only":
        fit = fit_one_variant(
            events,
            run_states,
            dilution_records,
            bin_number,
            "nominal",
            active_periods=(period,),
            initial_values=nominal["values"],
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
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
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
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
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
                ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
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
                    ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
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
    nominal_tasks = [
        {"kind": "nominal", "bin_number": bin_number}
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
                    "kind": "period_only",
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
    # Tiny positive radius includes boundary points without materially changing bins.
    return MplPath(np.asarray(vertices, dtype=float), closed=True).contains_points(pts, radius=1e-12)


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
        retry = Minuit(nll, u1=0.0, u2=0.0, lu1=0.0)
        retry.errordef = Minuit.LIKELIHOOD
        retry.print_level = 0
        retry.strategy = 2
        for name in ("u1", "u2", "lu1"):
            retry.limits[name] = PARAMETER_LIMITS[name]
        # endfor
        retry.migrad()
        retry.hesse()
        if (
            bool(retry.valid)
            and retry.covariance is not None
            and all(
                np.isfinite(float(retry.errors[name]))
                and float(retry.errors[name]) > 0.0
                for name in ("u1", "u2", "lu1")
            )
        ):
            m = retry
        # endif
    # endif
    return m, idx


def run_rga_cross_check(args):
    """Comparison-only extraction in the exact Diehl et al. binning."""
    out = args.output_dir.expanduser().resolve() / "rga_cross_check"
    tables = out / "tables"
    plots = out / "plots"
    cache_dir = out / "cache"
    for directory in (out, tables, plots, cache_dir):
        ensure_directory(directory)
    # endfor

    # Reuse/build the ordinary nominal selected-event sample.  The RGC
    # production exclusivity selection is preserved; only the comparison DIS
    # cuts and Diehl binning are changed.
    nominal_cache = (
        args.cache.expanduser().resolve()
        if args.cache
        else args.output_dir.expanduser().resolve()
        / "nominal/cache/selected_events.npz"
    )
    run_records = parse_run_info_csv(args.run_info_csv.expanduser().resolve())
    run_states = run_state_arrays(run_records)

    if nominal_cache.is_file():
        print(
            f"[RGA cross-check] loading nominal selected-event cache: "
            f"{nominal_cache}",
            flush=True,
        )
        events = load_event_cache(nominal_cache)
    else:
        print(
            "[RGA cross-check] nominal cache missing; building it directly.",
            flush=True,
        )
        inputs = {p: Path(v) for p, v in DEFAULT_INPUTS.items()}
        for period, path in args.input:
            inputs[period] = path
        # endfor
        cuts = load_channel_cuts(
            args.cut_json.expanduser().resolve(),
            cut_label="nominal",
        )
        build_event_cache(
            inputs,
            args.tree,
            args.chunk_size,
            run_records,
            cuts,
            nominal_cache,
        )
        events = load_event_cache(nominal_cache)
    # endif

    panel, subbin, gbin, y, mapping = _assign_rga_bins(events)
    print(
        "[RGA cross-check] fit model: "
        "(u1, u2, lu1) = "
        "(F_UU^cos(phi), F_UU^cos(2phi), F_LU^sin(phi))/F_UU; "
        "no production-bin amplitudes are imported.",
        flush=True,
    )
    print(
        f"[RGA cross-check] events after Diehl cuts and binning: "
        f"{np.count_nonzero(gbin > 0):,}",
        flush=True,
    )

    rows = []
    for g, (ip, it, lo, hi) in mapping.items():
        mask = gbin == g
        n = int(np.count_nonzero(mask))
        if n < 20:
            print(
                f"[RGA cross-check] panel {ip} t-bin {it}: "
                f"only {n} events; skipping",
                flush=True,
            )
            continue
        # endif

        m, idx = _fit_rga_uu_lu(events, run_states, mask)
        pub = RGA_PUBLISHED[ip][it - 1]

        fit_valid = (
            bool(m.valid)
            and m.covariance is not None
            and all(
                np.isfinite(float(m.errors[name]))
                and float(m.errors[name]) > 0.0
                for name in ("u1", "u2", "lu1")
            )
        )
        rgc_u1 = float(m.values["u1"])
        rgc_u2 = float(m.values["u2"])
        rgc_lu = float(m.values["lu1"])
        stat_u1 = float(m.errors["u1"])
        stat_u2 = float(m.errors["u2"])
        rgc_stat = float(m.errors["lu1"])

        correlation = np.full((3, 3), np.nan)
        if m.covariance is not None:
            names = ("u1", "u2", "lu1")
            cov = np.asarray(
                [
                    [float(m.covariance[a, b]) for b in names]
                    for a in names
                ],
                dtype=float,
            )
            sig = np.sqrt(np.maximum(np.diag(cov), 0.0))
            correlation = np.divide(
                cov,
                np.outer(sig, sig),
                out=np.full_like(cov, np.nan),
                where=np.outer(sig, sig) > 0.0,
            )
        # endif

        rga_q2, rga_x, rga_t, rga_eps, rga_val, rga_stat, rga_sys = pub
        if fit_valid:
            combined = math.sqrt(
                rgc_stat**2 + rga_stat**2 + rga_sys**2
            )
            delta = rgc_lu - rga_val
            pull = delta / combined if combined > 0.0 else np.nan
        else:
            combined = np.nan
            delta = np.nan
            pull = np.nan
        # endif

        epsilon_values = (
            events["epsilon"][idx]
            if "epsilon" in events
            else np.full(idx.size, np.nan)
        )
        rows.append(
            {
                "rga_panel": ip,
                "t_bin": it,
                "minus_t_low": lo,
                "minus_t_high": hi,
                "n_rgc": n,
                "mean_xB_rgc": float(np.mean(events["xB"][idx])),
                "mean_Q2_rgc": float(np.mean(events["Q2"][idx])),
                "mean_minus_t_rgc": float(
                    np.mean(events["minus_t"][idx])
                ),
                "mean_y_rgc": float(np.mean(y[idx])),
                "mean_epsilon_rgc": float(
                    np.nanmean(epsilon_values)
                ),
                "u1_rgc": rgc_u1,
                "stat_u1_rgc": stat_u1,
                "u2_rgc": rgc_u2,
                "stat_u2_rgc": stat_u2,
                "lu1_rgc": rgc_lu,
                "stat_rgc": rgc_stat,
                "fit_valid": fit_valid,
                "edm": float(m.fmin.edm),
                "corr_u1_u2": float(correlation[0, 1]),
                "corr_u1_lu1": float(correlation[0, 2]),
                "corr_u2_lu1": float(correlation[1, 2]),
                "xB_rga": rga_x,
                "Q2_rga": rga_q2,
                "minus_t_rga": rga_t,
                "epsilon_rga": rga_eps,
                "lu1_rga": rga_val,
                "stat_rga": rga_stat,
                "sys_rga": rga_sys,
                "delta_rgc_minus_rga": delta,
                "combined_uncertainty": combined,
                "pull": pull,
            }
        )

        status = "OK" if fit_valid else "INVALID"
        print(
            f"[RGA cross-check] panel {ip} t{it}: N={n:,}, "
            f"u1={rgc_u1:+.4f}+/-{stat_u1:.4f}, "
            f"u2={rgc_u2:+.4f}+/-{stat_u2:.4f}, "
            f"LU={rgc_lu:+.4f}+/-{rgc_stat:.4f}; "
            f"RGA={rga_val:+.4f} [{status}]",
            flush=True,
        )
    # endfor

    frame = pd.DataFrame(rows)
    csv = tables / "rga_cross_check.csv"
    frame.to_csv(csv, index=False)

    valid_frame = frame[frame["fit_valid"]].copy()

    if not args.skip_plots:
        fig, axes = plt.subplots(
            3, 3, figsize=(12, 10), sharey=True
        )
        for ip, ax in enumerate(axes.flat, start=1):
            data = valid_frame[valid_frame.rga_panel == ip]
            pub = np.asarray(RGA_PUBLISHED[ip], float)
            ax.errorbar(
                pub[:, 2],
                pub[:, 4],
                yerr=np.hypot(pub[:, 5], pub[:, 6]),
                fmt="o",
                label="RGA (Diehl et al.)",
            )
            if len(data):
                ax.errorbar(
                    data.mean_minus_t_rgc,
                    data.lu1_rgc,
                    yerr=data.stat_rgc,
                    fmt="s",
                    label="RGC cross-check",
                )
            # endif
            ax.axhline(0, lw=0.8)
            ax.set_title(f"RGA $Q^2$-$x_B$ bin {ip}")
            ax.set_xlabel(r"$-t$ (GeV$^2$)")
            if ip in (1, 4, 7):
                ax.set_ylabel(
                    r"$F_{LU}^{\sin\phi}/F_{UU}"
                    r"=\sigma_{LT'}/\sigma_0$"
                )
            # endif
            if ip == 1:
                ax.legend(fontsize=8)
            # endif
        # endfor
        fig.tight_layout()
        fig.savefig(
            plots / "rga_rgc_lu_overlay.png",
            dpi=200,
        )
        plt.close(fig)

        fig, axes = plt.subplots(
            3, 3, figsize=(12, 10), sharey=True
        )
        for ip, ax in enumerate(axes.flat, start=1):
            data = valid_frame[valid_frame.rga_panel == ip]
            ax.axhline(0, lw=0.8)
            if len(data):
                ax.plot(
                    data.mean_minus_t_rgc,
                    data.pull,
                    "o",
                )
            # endif
            ax.set_title(f"RGA $Q^2$-$x_B$ bin {ip}")
            ax.set_xlabel(r"$-t$ (GeV$^2$)")
            if ip in (1, 4, 7):
                ax.set_ylabel("RGC - RGA pull")
            # endif
        # endfor
        fig.tight_layout()
        fig.savefig(
            plots / "rga_rgc_lu_pulls.png",
            dpi=200,
        )
        plt.close(fig)
    # endif

    write_json(
        out / "rga_cross_check_manifest.json",
        {
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
            "observable": (
                "F_LU^{sin(phi)}/F_UU = sigma_LT'/sigma_0"
            ),
            "number_of_comparison_points": int(len(frame)),
            "number_of_valid_comparison_points": int(
                np.count_nonzero(frame["fit_valid"])
            ),
            "important_note": (
                "Comparison-only three-parameter UU+LU fit performed "
                "directly in each Diehl (Q2,xB,-t) bin. No amplitudes "
                "or dilution factors are imported from the nominal "
                "24-bin production fit. Target-spin UL/LL terms are "
                "omitted in this first cross-check implementation. "
                "The production extraction is unchanged."
            ),
        },
    )
    print(
        f"[RGA cross-check] wrote {csv}",
        flush=True,
    )
    return 0


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
    if args.rga_cross_check:
        return run_rga_cross_check(args)
    # endif
    workers = max(
        1,
        min(int(args.workers), MAXIMUM_WORKERS, os.cpu_count() or 1, NUMBER_OF_BINS),
    )
    root = args.output_dir.expanduser().resolve()
    nominal_dir = root / "nominal"
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
    isr_dilution = None
    if not args.disable_isr:
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
    print(f"ISR study enabled:    {not args.disable_isr}", flush=True)
    print(
        f"Momentum study:       {not args.disable_momentum_corrections}",
        flush=True,
    )
    print(
        f"Channel-selection:    {not args.disable_channel_selection}",
        flush=True,
    )
    print(f"Nominal cache:        {nominal_cache}", flush=True)
    if not args.disable_channel_selection:
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
    if not args.disable_channel_selection:
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
        include_target_axis_study=True,
        include_period_diagnostics=True,
        cut_label="nominal",
        source_cache_path=(None if args.reuse_cache else nominal_source_cache),
    )

    target_axis_study = write_target_axis_study_products(
        nominal_result["frame"], diagnostics_dir / "target_axis"
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
