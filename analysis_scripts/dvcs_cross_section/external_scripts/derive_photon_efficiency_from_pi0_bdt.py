#!/usr/bin/env python3
"""
First-stage pi0-BDT tag-and-probe photon-efficiency extraction.

Purpose
-------
Use the already-trained single-photon pi0 BDT as an aggressive TAG-purity
selector, then measure the probability to reconstruct the kinematically
predicted partner photon from pi0 -> gamma gamma.

The nominal first-stage result is only a function of predicted probe energy,
with FT and FD treated separately.  Sector/theta dependence is deliberately
left for a later iteration.

This version also performs the validation studies needed before interpreting
the efficiency:
    1) denominator closure: compare the missing-vector prediction against the
       reconstructed truth-qualified sister photon in MC;
    2) angular-association plateau: scan the predicted/reconstructed opening
       angle and determine where correct-pair recovery saturates;
    3) extract the numerator from fits to M(gamma gamma), in bins of the SAME
       predicted probe energy used for the denominator;
    4) repeat those fits for 8, 10, and 12 degree association windows;
    5) build a fine probe-energy response diagnostic for AAOgen and CLASDIS.

The raw efficiency in earlier versions was already binned in predicted probe
energy for BOTH numerator and denominator.  The new fit-based numerator keeps
that definition; it replaces raw numerator counts with fitted pi0 yields.

The current ROOT schema does not store generated sister-photon four-vectors.
Therefore the denominator-closure study uses the reconstructed sister photon
after truth qualification.  This tests the missing-vector bin assignment plus
detector resolution; it is not a generator-level resolution study.

Core philosophy
---------------
DENOMINATOR (e'p'gamma X tree):
    * reconstructed electron + proton + one reconstructed photon (the tag)
    * tag passes the same loose single-photon BDT preselection
    * tag BDT score >= an aggressive threshold (default 0.80)
    * tag energy >= 0.4 GeV (the reconstructed-photon threshold)
    * modest low-Mx2(ep) requirement to suppress obviously heavier hadronic
      systems (default Mx2(ep) < 0.25 GeV^2)
    * NO tight Mx2(epgamma_tag) requirement in the nominal sample: the
      aggressiveness is intentionally supplied by the BDT score, not by
      kinematic sculpting
    * the missing four-vector after e'p'gamma is interpreted as the predicted
      probe photon
    * predicted probe E >= 0.4 GeV
    * predicted probe direction must point to FT or FD acceptance

NUMERATOR (e'p'gammagamma X tree):
    * build each reconstructed pair as two DIRECTED tag -> probe trials
    * apply the SAME tag-side selection as the denominator
    * require the second reconstructed photon to be geometrically compatible
      with the predicted missing-photon direction
    * DO NOT require a BDT score on the probe: that would fold classifier
      efficiency into photon reconstruction efficiency

Then:
    epsilon_gamma(E_pred) = N_reconstructed_probe(E_pred) /
                            N_predicted_probe(E_pred)

This is done independently for:
    * real data
    * AAOgen (nominal exclusive-pi0 MC reference)
    * CLASDIS (independent hadronic-environment cross-check)

and the first correction diagnostic is:
    C(E_pred) = epsilon_data / epsilon_AAOgen

Important
---------
The e'pgamma X and e'pgammagamma X inputs for a given source MUST represent
the same processed exposure/sample before the resulting efficiency is treated
as physical.  The code does not luminosity-normalize one topology to the other.

This first version is intentionally aggressive and diagnostic-heavy.  It is
not yet the final photon-efficiency correction implementation.
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

import joblib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import uproot


# =============================================================================
# Constants / configuration
# =============================================================================

START_TIME = time.perf_counter()

PROTON_MASS = 0.9382720813

FEATURES = [
    "p2_p",
    "p2_theta",
    "open_angle_ep2",
    "open_angle_p1p2",
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    "Emiss2",
    "theta_gamma_gamma",
    "pTmiss",
]

REGIONS = {
    "FT": 0,
    "FD": 1,
}

BEAM_ENERGIES = {
    "sp18_inb": 10.594,
    "sp18_out": 10.594,
    "fa18_inb": 10.604,
    "fa18_out": 10.604,
    "sp19_inb": 10.200,
}

BASE_DIRS = {
    period: Path(f"/work/clas12/thayward/pi0_BDT/training_sample/{period}")
    for period in BEAM_ENERGIES
}

DEFAULT_ENERGY_BINS = np.asarray(
    [0.40, 1.00, 2.00, 3.00, 9.50],
    dtype=float,
)

# BDT preselection -- kept consistent with training.
OPEN_ANGLE_MIN = 7.0
LOOSE_MX2_MIN = -0.5
LOOSE_MX2_MAX = 0.5

# Aggressive nominal tag-and-probe selection.
#
# IMPORTANT: the aggressiveness comes from the BDT purity requirement, not
# from imposing a high tag-energy threshold.  Requiring E_tag >= 2 GeV would
# bias the directed tag-and-probe phase space: in particular it would remove
# high-energy probes whose partner photon is the lower-energy daughter.  Keep
# the tag threshold at the ordinary reconstructed-photon threshold so the
# probe can populate the full accessible energy spectrum.
TAG_SCORE_MIN = 0.80
TAG_ENERGY_MIN = 0.40
PROBE_ENERGY_MIN = 0.40
PROBE_ENERGY_MAX = 9.50

# Low Mx2(ep) is kept only as a modest hadronic-environment cleanup.
# The dominant purity requirement is the BDT score.
MX2_EP_MIN = -0.10
MX2_EP_MAX = 0.25

# Retained only for plotting / optional future studies.  The nominal v3
# selection does NOT cut on Mx2(epgamma_tag), because that strongly sculpted
# the predicted-probe energy and removed the FT probe population.
MISSING_MASS2_MIN = -0.08
MISSING_MASS2_MAX = 0.08

# Predicted probe acceptance.  Initial broad detector definitions only.
FT_THETA_MIN = 2.0
FT_THETA_MAX = 5.5
FD_THETA_MIN = 5.5
FD_THETA_MAX = 35.0

# Numerator: reconstructed partner must point near the predicted probe.
PROBE_MATCH_ANGLE_MAX = 3.0

# Pair mass window is a DIAGNOSTIC only in v1 and is not required in the
# numerator, because doing so would add an extra reconstruction-resolution
# requirement that is absent from the denominator.
PI0_MASS_WINDOW = (0.110, 0.160)

# Angular-association validation scan.  The nominal 3-degree requirement is
# NOT changed by this study; the scan is used to determine whether a looser
# matching window reaches a stable plateau.
ASSOCIATION_ANGLE_SCAN = np.asarray(
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 15.0],
    dtype=float,
)

ASSOCIATION_ENERGY_BINS = (
    (0.40, 1.00),
    (1.00, 2.00),
    (2.00, 9.50),
)

# Association windows used for the fitted-yield stability test.
# 10 degrees is the nominal working point; 8 and 12 degrees bracket it as
# matching/association variations.
MGG_FIT_ASSOCIATION_ANGLES = (8.0, 10.0, 12.0)
MGG_FIT_NOMINAL_ANGLE = 10.0

# Diagnostic scan of the minimum BDT tag score.  This deliberately extends
# below the nominal 0.80 working point and very close to 1.0 so we can test
# whether denominator contamination is responsible for an anomalously low
# apparent photon efficiency in data.
BDT_EFFICIENCY_SCAN_THRESHOLDS = np.asarray(
    [0.50, 0.60, 0.70, 0.80, 0.85, 0.90, 0.925, 0.95, 0.975],
    dtype=float,
)

# Coarse probe-polar-angle dependence for the correction.
# FT remains one angular region.  FD is split into three broad bins,
# motivated by the strong score-dependent theta migration seen in AAOgen.
FD_PROBE_THETA_SPLITS = (10.0, 15.0)

# Fit range deliberately wider than the nominal pi0 window so the continuum
# is constrained by data rather than by a hard mass-window count.
MGG_FIT_RANGE = (0.070, 0.240)
MGG_FIT_BINS = 85

# Fine response binning used only as a migration diagnostic.
RESPONSE_ENERGY_EDGES = np.asarray(
    [0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 2.00, 2.50, 3.50, 5.00, 7.00, 9.50],
    dtype=float,
)


# =============================================================================
# Small containers
# =============================================================================

@dataclass
class DirectedSample:
    predicted_energy: np.ndarray
    predicted_theta: np.ndarray
    tag_score: np.ndarray
    tag_energy: np.ndarray
    mx2_ep: np.ndarray
    missing_mass2: np.ndarray
    region: np.ndarray
    partner_match_angle: Optional[np.ndarray] = None
    partner_energy: Optional[np.ndarray] = None
    pair_mass: Optional[np.ndarray] = None
    tag_parent_pid: Optional[np.ndarray] = None


@dataclass
class EfficiencyCurve:
    edges: np.ndarray
    centers: np.ndarray
    denominator: np.ndarray
    numerator: np.ndarray
    efficiency: np.ndarray
    uncertainty: np.ndarray


@dataclass
class MassFitResult:
    success: bool
    yield_pi0: float
    yield_uncertainty: float
    mean: float
    sigma: float
    sigma_core: float
    sigma_tail: float
    core_fraction: float
    chi2_ndf: float
    x_centers: np.ndarray
    counts: np.ndarray
    model_counts: np.ndarray
    signal_counts: np.ndarray
    background_counts: np.ndarray


@dataclass
class FitEfficiencyCurve:
    edges: np.ndarray
    centers: np.ndarray
    denominator: np.ndarray
    numerator_yield: np.ndarray
    numerator_uncertainty: np.ndarray
    efficiency: np.ndarray
    uncertainty: np.ndarray
    fit_mean: np.ndarray
    fit_sigma: np.ndarray
    fit_chi2_ndf: np.ndarray


# =============================================================================
# CLI / paths
# =============================================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="First-stage pi0-BDT tag-and-probe photon efficiency."
    )

    parser.add_argument(
        "--period",
        choices=sorted(BEAM_ENERGIES),
        default="fa18_inb",
    )
    parser.add_argument(
        "--sample-tag",
        default="temp",
        help="Input suffix tag. 'temp' -> *_temp.root; 'full' -> *.root.",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=None,
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help=(
            "Directory containing FT/pi0_bdt_model.joblib and "
            "FD/pi0_bdt_model.joblib. Default is "
            "output/pi0_bdt_truth/<period>/<sample-tag>."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
    )

    parser.add_argument("--tag-score-min", type=float, default=TAG_SCORE_MIN)
    parser.add_argument("--tag-energy-min", type=float, default=TAG_ENERGY_MIN)
    parser.add_argument("--probe-energy-min", type=float, default=PROBE_ENERGY_MIN)
    parser.add_argument("--probe-energy-max", type=float, default=PROBE_ENERGY_MAX)

    parser.add_argument("--mx2-ep-min", type=float, default=MX2_EP_MIN)
    parser.add_argument("--mx2-ep-max", type=float, default=MX2_EP_MAX)
    parser.add_argument(
        "--missing-mass2-min",
        type=float,
        default=MISSING_MASS2_MIN,
    )
    parser.add_argument(
        "--missing-mass2-max",
        type=float,
        default=MISSING_MASS2_MAX,
    )
    parser.add_argument(
        "--probe-match-angle-max",
        type=float,
        default=PROBE_MATCH_ANGLE_MAX,
    )

    parser.add_argument("--ft-theta-min", type=float, default=FT_THETA_MIN)
    parser.add_argument("--ft-theta-max", type=float, default=FT_THETA_MAX)
    parser.add_argument("--fd-theta-min", type=float, default=FD_THETA_MIN)
    parser.add_argument("--fd-theta-max", type=float, default=FD_THETA_MAX)

    parser.add_argument(
        "--energy-bins",
        default=",".join(f"{x:g}" for x in DEFAULT_ENERGY_BINS),
        help="Comma-separated predicted-probe-energy bin edges in GeV.",
    )

    return parser.parse_args()


def sample_suffix(sample_tag: str) -> str:
    tag = sample_tag.strip()

    if tag.lower() in {"", "full", "none"}:
        return ""
    #endif

    return tag if tag.startswith("_") else "_" + tag


def default_model_dir(args: argparse.Namespace) -> Path:
    if args.model_dir is not None:
        return args.model_dir
    #endif

    return (
        Path("output")
        / "pi0_bdt_truth"
        / args.period
        / args.sample_tag
    )


def default_output_dir(args: argparse.Namespace) -> Path:
    if args.output_dir is not None:
        return args.output_dir
    #endif

    return (
        Path("output")
        / "photon_efficiency_bdt"
        / args.period
        / args.sample_tag
    )


def build_paths(args: argparse.Namespace) -> Dict[str, Path]:
    base = args.base_dir if args.base_dir is not None else BASE_DIRS[args.period]
    suffix = sample_suffix(args.sample_tag)
    prefix = f"rga_{args.period}"

    return {
        "data_epg": base / f"data_{prefix}_epgammaX{suffix}.root",
        "data_epgg": base / f"data_{prefix}_epgammagammaX{suffix}.root",
        "aaogen_epg": base / f"aaogen_{prefix}_epgammaX{suffix}.root",
        "aaogen_epgg": base / f"aaogen_{prefix}_epgammagammaX{suffix}.root",
        "clasdis_epg": base / f"clasdis_{prefix}_epgammaX{suffix}.root",
        "clasdis_epgg": base / f"clasdis_{prefix}_epgammagammaX{suffix}.root",
    }


# =============================================================================
# Utilities
# =============================================================================

def progress(message: str) -> None:
    elapsed = time.perf_counter() - START_TIME
    print(f"[+{elapsed:8.1f} s] {message}", flush=True)


def first_tree_name(path: Path) -> str:
    with uproot.open(path) as root_file:
        for key, obj in root_file.items():
            if hasattr(obj, "arrays") and hasattr(obj, "num_entries"):
                return key.split(";")[0]
            #endif
        #endfor
    #endwith

    raise RuntimeError(f"No TTree-like object found in {path}")


def available_branches(path: Path) -> set[str]:
    tree_name = first_tree_name(path)

    with uproot.open(path) as root_file:
        return {
            name.split(";")[0]
            for name in root_file[tree_name].keys()
        }
    #endwith


def read_arrays(
    path: Path,
    required: Sequence[str],
    optional: Sequence[str] = (),
) -> Dict[str, np.ndarray]:
    tree_name = first_tree_name(path)
    available = available_branches(path)

    missing = [branch for branch in required if branch not in available]

    if missing:
        raise RuntimeError(
            f"{path} is missing required branches:\n  "
            + "\n  ".join(missing)
        )
    #endif

    branches = list(required) + [
        branch
        for branch in optional
        if branch in available and branch not in required
    ]

    progress(f"LOAD: {path}")
    start = time.perf_counter()

    with uproot.open(path) as root_file:
        raw = root_file[tree_name].arrays(branches, library="np")
    #endwith

    arrays = {branch: np.asarray(raw[branch]) for branch in branches}
    n = len(next(iter(arrays.values()))) if arrays else 0

    progress(
        f"LOAD COMPLETE: {n:,} rows, {len(branches)} branches "
        f"({time.perf_counter() - start:.1f} s)"
    )
    return arrays


def angle_to_radians(values: np.ndarray) -> np.ndarray:
    """
    Convert stored reconstruction angles to radians.

    The epgamma / epgammagamma ROOT trees used by this analysis store the
    reconstructed spherical angles in radians:
      theta ~ O(0.1-1.0)
      phi   ~ 0..2*pi

    Earlier versions used a magnitude heuristic that treated phi values above
    3.5 as degrees and therefore applied np.radians() to an angle that was
    already in radians.  That compressed phi by a factor of ~57.3 and broke the
    missing-probe direction while leaving the missing energy unchanged.

    For these trees, return the stored angle directly.
    """
    return np.asarray(values, dtype=float)


def angle_to_degrees(values: np.ndarray) -> np.ndarray:
    """
    Convert the angular BDT features to degrees.

    This helper is used only for theta-like BDT inputs, not phi.  Their stored
    values are in radians and occupy the small-angle CLAS12 range, so values
    below 3.5 are converted to degrees.  Do not use this helper for azimuth.
    """
    arr = np.asarray(values, dtype=float)
    finite = arr[np.isfinite(arr)]

    if len(finite) == 0:
        return arr
    #endif

    if np.nanpercentile(np.abs(finite), 99.0) <= 3.5:
        return np.degrees(arr)
    #endif

    return arr



def spherical_to_cartesian_explicit(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    theta_units: str,
    phi_units: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert spherical coordinates using an EXPLICIT unit hypothesis.

    This exists only for the geometry-debug study.  It intentionally avoids
    the automatic branch-level angle-unit heuristic so that we can determine
    unambiguously which convention reproduces the reconstructed partner
    photon in e'p'gammagamma events.
    """
    p = np.asarray(momentum, dtype=float)
    theta_arr = np.asarray(theta, dtype=float)
    phi_arr = np.asarray(phi, dtype=float)

    if theta_units == "deg":
        theta_rad = np.radians(theta_arr)
    elif theta_units == "rad":
        theta_rad = theta_arr
    else:
        raise ValueError(f"Unsupported theta_units={theta_units}")
    #endif

    if phi_units == "deg":
        phi_rad = np.radians(phi_arr)
    elif phi_units == "rad":
        phi_rad = phi_arr
    else:
        raise ValueError(f"Unsupported phi_units={phi_units}")
    #endif

    px = p * np.sin(theta_rad) * np.cos(phi_rad)
    py = p * np.sin(theta_rad) * np.sin(phi_rad)
    pz = p * np.cos(theta_rad)
    return px, py, pz


def vector_phi_deg(
    px: np.ndarray,
    py: np.ndarray,
) -> np.ndarray:
    return np.degrees(np.arctan2(py, px))


def wrapped_delta_phi_deg(
    phi_a_deg: np.ndarray,
    phi_b_deg: np.ndarray,
) -> np.ndarray:
    return (
        np.asarray(phi_a_deg, dtype=float)
        - np.asarray(phi_b_deg, dtype=float)
        + 180.0
    ) % 360.0 - 180.0


def spherical_to_cartesian(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta_rad = angle_to_radians(theta)
    phi_rad = angle_to_radians(phi)
    p = np.asarray(momentum, dtype=float)

    px = p * np.sin(theta_rad) * np.cos(phi_rad)
    py = p * np.sin(theta_rad) * np.sin(phi_rad)
    pz = p * np.cos(theta_rad)
    return px, py, pz


def vector_theta_deg(
    px: np.ndarray,
    py: np.ndarray,
    pz: np.ndarray,
) -> np.ndarray:
    p = np.sqrt(px * px + py * py + pz * pz)
    cos_theta = np.divide(
        pz,
        p,
        out=np.full_like(p, np.nan, dtype=float),
        where=p > 0.0,
    )
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))


def opening_angle_deg(
    ax: np.ndarray,
    ay: np.ndarray,
    az: np.ndarray,
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
) -> np.ndarray:
    adotb = ax * bx + ay * by + az * bz
    amag = np.sqrt(ax * ax + ay * ay + az * az)
    bmag = np.sqrt(bx * bx + by * by + bz * bz)
    denom = amag * bmag

    cos_angle = np.divide(
        adotb,
        denom,
        out=np.full_like(adotb, np.nan, dtype=float),
        where=denom > 0.0,
    )
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    return np.degrees(np.arccos(cos_angle))




def assign_fd_sector(
    phi_deg: np.ndarray,
) -> np.ndarray:
    """
    Assign CLAS12 FD sector from predicted probe azimuth.

    Sector convention:
      S1: 330-360 or 0-30 deg
      S2: 30-90 deg
      S3: 90-150 deg
      S4: 150-210 deg
      S5: 210-270 deg
      S6: 270-330 deg

    Returns 1-6, or -1 for non-finite values.
    """
    phi = np.asarray(phi_deg, dtype=float)
    wrapped = np.mod(phi, 360.0)
    sector = np.full(len(wrapped), -1, dtype=np.int8)

    finite = np.isfinite(wrapped)
    sector[
        finite
        & ((wrapped >= 330.0) | (wrapped < 30.0))
    ] = 1
    sector[
        finite
        & (wrapped >= 30.0)
        & (wrapped < 90.0)
    ] = 2
    sector[
        finite
        & (wrapped >= 90.0)
        & (wrapped < 150.0)
    ] = 3
    sector[
        finite
        & (wrapped >= 150.0)
        & (wrapped < 210.0)
    ] = 4
    sector[
        finite
        & (wrapped >= 210.0)
        & (wrapped < 270.0)
    ] = 5
    sector[
        finite
        & (wrapped >= 270.0)
        & (wrapped < 330.0)
    ] = 6

    return sector


def probe_theta_bins(
    args: argparse.Namespace,
) -> list[Tuple[str, int, float, float]]:
    """
    Return the coarse probe-theta bins used by the correction study.

    FT is kept as one bin.  FD is divided into 5.5-10, 10-15, and 15-35 deg
    for the default acceptance.  The outer limits continue to respect the
    command-line FT/FD acceptance bounds.
    """
    split1, split2 = FD_PROBE_THETA_SPLITS

    if not (
        args.fd_theta_min < split1
        and split1 < split2
        and split2 < args.fd_theta_max
    ):
        raise ValueError(
            "FD theta acceptance must contain the fixed coarse splits "
            f"{split1:g} and {split2:g} deg."
        )
    #endif

    return [
        (
            f"FT {args.ft_theta_min:g}-{args.ft_theta_max:g} deg",
            REGIONS["FT"],
            float(args.ft_theta_min),
            float(args.ft_theta_max),
        ),
        (
            f"FD {args.fd_theta_min:g}-{split1:g} deg",
            REGIONS["FD"],
            float(args.fd_theta_min),
            float(split1),
        ),
        (
            f"FD {split1:g}-{split2:g} deg",
            REGIONS["FD"],
            float(split1),
            float(split2),
        ),
        (
            f"FD {split2:g}-{args.fd_theta_max:g} deg",
            REGIONS["FD"],
            float(split2),
            float(args.fd_theta_max),
        ),
    ]


def assign_region(
    theta_deg: np.ndarray,
    args: argparse.Namespace,
) -> np.ndarray:
    theta = np.asarray(theta_deg, dtype=float)
    region = np.full(len(theta), -1, dtype=np.int8)

    ft = (
        np.isfinite(theta)
        & (theta >= args.ft_theta_min)
        & (theta < args.ft_theta_max)
    )
    fd = (
        np.isfinite(theta)
        & (theta >= args.fd_theta_min)
        & (theta < args.fd_theta_max)
    )

    region[ft] = REGIONS["FT"]
    region[fd] = REGIONS["FD"]
    return region


def feature_matrix_epg(arrays: Mapping[str, np.ndarray]) -> np.ndarray:
    columns = []

    for feature in FEATURES:
        values = np.asarray(arrays[feature], dtype=float)

        if feature in {"p2_theta", "theta_gamma_gamma"}:
            values = angle_to_degrees(values)
        #endif

        columns.append(values)
    #endfor

    return np.column_stack(columns)


def feature_matrix_epgg(
    arrays: Mapping[str, np.ndarray],
    tag_index: int,
) -> np.ndarray:
    if tag_index == 1:
        mapping = {
            "p2_p": "p2_p",
            "p2_theta": "p2_theta",
            "open_angle_ep2": "open_angle_ep2",
            "open_angle_p1p2": "open_angle_p1p2",
            "Mx2": "gamma1_epgamma_Mx2",
            "Mx2_1": "gamma1_ep_Mx2",
            "Mx2_2": "gamma1_egamma_Mx2",
            "Emiss2": "gamma1_Emiss2",
            "theta_gamma_gamma": "gamma1_theta_gamma_gamma",
            "pTmiss": "gamma1_pTmiss",
        }
    elif tag_index == 2:
        mapping = {
            "p2_p": "p3_p",
            "p2_theta": "p3_theta",
            "open_angle_ep2": "open_angle_ep3",
            "open_angle_p1p2": "open_angle_p1p3",
            "Mx2": "gamma2_epgamma_Mx2",
            "Mx2_1": "gamma2_ep_Mx2",
            "Mx2_2": "gamma2_egamma_Mx2",
            "Emiss2": "gamma2_Emiss2",
            "theta_gamma_gamma": "gamma2_theta_gamma_gamma",
            "pTmiss": "gamma2_pTmiss",
        }
    else:
        raise ValueError("tag_index must be 1 or 2")
    #endif

    columns = []

    for feature in FEATURES:
        values = np.asarray(arrays[mapping[feature]], dtype=float)

        if feature in {"p2_theta", "theta_gamma_gamma"}:
            values = angle_to_degrees(values)
        #endif

        columns.append(values)
    #endfor

    return np.column_stack(columns)


def score_model(model, X: np.ndarray) -> np.ndarray:
    return np.asarray(model.predict_proba(X)[:, 1], dtype=float)


def load_models(model_dir: Path) -> Dict[str, object]:
    models = {}

    for region in REGIONS:
        path = model_dir / region / "pi0_bdt_model.joblib"

        if not path.exists():
            raise FileNotFoundError(f"Missing BDT model bundle: {path}")
        #endif

        bundle = joblib.load(path)
        model = bundle["model"] if isinstance(bundle, dict) else bundle
        models[region] = model
        progress(f"MODEL LOADED: {path}")
    #endfor

    return models


# =============================================================================
# Required ROOT branches
# =============================================================================

EPG_REQUIRED = sorted(
    set(
        FEATURES
        + [
            "detector2",
            "e_p",
            "e_theta",
            "e_phi",
            "p1_p",
            "p1_theta",
            "p1_phi",
            "p2_phi",
        ]
    )
)

EPG_OPTIONAL_TRUTH = [
    "matching_gamma_pid",
    "gamma_mcindex",
    "gamma_parent_pid",
    "gamma_grandparent_pid",
]

EPGG_REQUIRED = sorted(
    set(
        [
            "detector_gamma1",
            "detector_gamma2",
            "e_p",
            "e_theta",
            "e_phi",
            "p1_p",
            "p1_theta",
            "p1_phi",
            "p2_p",
            "p2_theta",
            "p2_phi",
            "p3_p",
            "p3_theta",
            "p3_phi",
            "open_angle_ep2",
            "open_angle_ep3",
            "open_angle_p1p2",
            "open_angle_p1p3",
            "gamma1_epgamma_Mx2",
            "gamma1_ep_Mx2",
            "gamma1_egamma_Mx2",
            "gamma1_Emiss2",
            "gamma1_theta_gamma_gamma",
            "gamma1_pTmiss",
            "gamma2_epgamma_Mx2",
            "gamma2_ep_Mx2",
            "gamma2_egamma_Mx2",
            "gamma2_Emiss2",
            "gamma2_theta_gamma_gamma",
            "gamma2_pTmiss",
            "Mh_gammagamma",
        ]
    )
)

EPGG_OPTIONAL_TRUTH = [
    "matching_gamma1_pid",
    "gamma1_mcindex",
    "gamma1_parent_pid",
    "gamma1_grandparent_pid",
    "matching_gamma2_pid",
    "gamma2_mcindex",
    "gamma2_parent_pid",
    "gamma2_grandparent_pid",
]


# =============================================================================
# Missing-probe reconstruction
# =============================================================================

def predicted_probe_from_tag(
    arrays: Mapping[str, np.ndarray],
    beam_energy: float,
    tag_p_branch: str,
    tag_theta_branch: str,
    tag_phi_branch: str,
) -> Dict[str, np.ndarray]:
    epx, epy, epz = spherical_to_cartesian(
        arrays["e_p"],
        arrays["e_theta"],
        arrays["e_phi"],
    )
    ppx, ppy, ppz = spherical_to_cartesian(
        arrays["p1_p"],
        arrays["p1_theta"],
        arrays["p1_phi"],
    )
    gpx, gpy, gpz = spherical_to_cartesian(
        arrays[tag_p_branch],
        arrays[tag_theta_branch],
        arrays[tag_phi_branch],
    )

    e_p = np.asarray(arrays["e_p"], dtype=float)
    proton_p = np.asarray(arrays["p1_p"], dtype=float)
    gamma_p = np.asarray(arrays[tag_p_branch], dtype=float)

    electron_energy = e_p
    proton_energy = np.sqrt(proton_p * proton_p + PROTON_MASS * PROTON_MASS)
    gamma_energy = gamma_p

    missing_energy = (
        beam_energy
        + PROTON_MASS
        - electron_energy
        - proton_energy
        - gamma_energy
    )

    missing_px = -epx - ppx - gpx
    missing_py = -epy - ppy - gpy
    missing_pz = beam_energy - epz - ppz - gpz
    missing_p = np.sqrt(
        missing_px * missing_px
        + missing_py * missing_py
        + missing_pz * missing_pz
    )
    missing_mass2 = missing_energy * missing_energy - missing_p * missing_p
    missing_theta = vector_theta_deg(
        missing_px,
        missing_py,
        missing_pz,
    )

    return {
        "energy": missing_energy,
        "px": missing_px,
        "py": missing_py,
        "pz": missing_pz,
        "p": missing_p,
        "mass2": missing_mass2,
        "theta": missing_theta,
        "phi": vector_phi_deg(
            missing_px,
            missing_py,
        ),
    }



def predicted_probe_from_tag_explicit_units(
    arrays: Mapping[str, np.ndarray],
    beam_energy: float,
    tag_p_branch: str,
    tag_theta_branch: str,
    tag_phi_branch: str,
    theta_units: str,
    phi_units: str,
) -> Dict[str, np.ndarray]:
    """
    Missing-photon four-vector under a specified angular-unit hypothesis.
    """
    epx, epy, epz = spherical_to_cartesian_explicit(
        arrays["e_p"],
        arrays["e_theta"],
        arrays["e_phi"],
        theta_units,
        phi_units,
    )
    ppx, ppy, ppz = spherical_to_cartesian_explicit(
        arrays["p1_p"],
        arrays["p1_theta"],
        arrays["p1_phi"],
        theta_units,
        phi_units,
    )
    gpx, gpy, gpz = spherical_to_cartesian_explicit(
        arrays[tag_p_branch],
        arrays[tag_theta_branch],
        arrays[tag_phi_branch],
        theta_units,
        phi_units,
    )

    e_p = np.asarray(arrays["e_p"], dtype=float)
    proton_p = np.asarray(arrays["p1_p"], dtype=float)
    gamma_p = np.asarray(arrays[tag_p_branch], dtype=float)

    electron_energy = e_p
    proton_energy = np.sqrt(
        proton_p * proton_p + PROTON_MASS * PROTON_MASS
    )
    gamma_energy = gamma_p

    missing_energy = (
        beam_energy
        + PROTON_MASS
        - electron_energy
        - proton_energy
        - gamma_energy
    )

    missing_px = -epx - ppx - gpx
    missing_py = -epy - ppy - gpy
    missing_pz = beam_energy - epz - ppz - gpz
    missing_p = np.sqrt(
        missing_px * missing_px
        + missing_py * missing_py
        + missing_pz * missing_pz
    )

    return {
        "energy": missing_energy,
        "px": missing_px,
        "py": missing_py,
        "pz": missing_pz,
        "p": missing_p,
        "theta": vector_theta_deg(
            missing_px,
            missing_py,
            missing_pz,
        ),
        "phi": vector_phi_deg(
            missing_px,
            missing_py,
        ),
    }


ANGLE_HYPOTHESES = {
    "theta=rad, phi=rad": ("rad", "rad"),
    "theta=deg, phi=deg": ("deg", "deg"),
    "theta=deg, phi=rad": ("deg", "rad"),
    "theta=rad, phi=deg": ("rad", "deg"),
}


def geometry_hypothesis_diagnostics(
    arrays: Mapping[str, np.ndarray],
    beam_energy: float,
    tag_index: int,
    base_mask: np.ndarray,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    Compare several explicit angular-unit hypotheses before imposing any
    predicted FT/FD acceptance or angular matching.

    The correct convention should give a narrow peak near zero in the full 3D
    opening angle and in Delta-theta / Delta-phi for true exclusive pi0 pairs.
    """
    if tag_index == 1:
        tag_p = "p2_p"
        tag_theta = "p2_theta"
        tag_phi = "p2_phi"
        probe_p = "p3_p"
        probe_theta = "p3_theta"
        probe_phi = "p3_phi"
    elif tag_index == 2:
        tag_p = "p3_p"
        tag_theta = "p3_theta"
        tag_phi = "p3_phi"
        probe_p = "p2_p"
        probe_theta = "p2_theta"
        probe_phi = "p2_phi"
    else:
        raise ValueError("tag_index must be 1 or 2")
    #endif

    out: Dict[str, Dict[str, np.ndarray]] = {}

    for label, (theta_units, phi_units) in ANGLE_HYPOTHESES.items():
        pred = predicted_probe_from_tag_explicit_units(
            arrays,
            beam_energy,
            tag_p,
            tag_theta,
            tag_phi,
            theta_units,
            phi_units,
        )
        rpx, rpy, rpz = spherical_to_cartesian_explicit(
            arrays[probe_p],
            arrays[probe_theta],
            arrays[probe_phi],
            theta_units,
            phi_units,
        )
        reco_theta = vector_theta_deg(rpx, rpy, rpz)
        reco_phi = vector_phi_deg(rpx, rpy)

        opening = opening_angle_deg(
            pred["px"],
            pred["py"],
            pred["pz"],
            rpx,
            rpy,
            rpz,
        )
        delta_theta = pred["theta"] - reco_theta
        delta_phi = wrapped_delta_phi_deg(pred["phi"], reco_phi)

        finite = (
            np.asarray(base_mask, dtype=bool)
            & np.isfinite(opening)
            & np.isfinite(delta_theta)
            & np.isfinite(delta_phi)
            & np.isfinite(pred["energy"])
        )

        out[label] = {
            "mask": finite,
            "opening": opening,
            "delta_theta": delta_theta,
            "delta_phi": delta_phi,
            "pred_theta": pred["theta"],
            "pred_phi": pred["phi"],
            "pred_energy": pred["energy"],
            "reco_theta": reco_theta,
            "reco_phi": reco_phi,
        }
    #endfor

    return out


# =============================================================================
# Directed denominator / numerator construction
# =============================================================================

def score_epg_tags(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, object],
) -> np.ndarray:
    X = feature_matrix_epg(arrays)
    detector = np.asarray(arrays["detector2"], dtype=np.int64)

    valid = np.all(np.isfinite(X), axis=1)
    valid &= X[:, FEATURES.index("open_angle_ep2")] > OPEN_ANGLE_MIN
    valid &= X[:, FEATURES.index("Mx2")] > LOOSE_MX2_MIN
    valid &= X[:, FEATURES.index("Mx2")] < LOOSE_MX2_MAX

    scores = np.full(len(X), np.nan, dtype=float)

    for region_name, detector_value in REGIONS.items():
        mask = valid & (detector == detector_value)

        if np.any(mask):
            scores[mask] = score_model(models[region_name], X[mask])
        #endif
    #endfor

    return scores


def denominator_from_epg(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, object],
    beam_energy: float,
    args: argparse.Namespace,
) -> Tuple[DirectedSample, Dict[str, np.ndarray]]:
    X = feature_matrix_epg(arrays)
    score = score_epg_tags(arrays, models)
    pred = predicted_probe_from_tag(
        arrays,
        beam_energy,
        "p2_p",
        "p2_theta",
        "p2_phi",
    )

    mx2_ep = X[:, FEATURES.index("Mx2_1")]
    branch_missing_mass2 = X[:, FEATURES.index("Mx2")]
    tag_energy = np.asarray(arrays["p2_p"], dtype=float)
    region = assign_region(pred["theta"], args)

    base = np.isfinite(score)
    base &= score >= args.tag_score_min
    base &= tag_energy >= args.tag_energy_min
    base &= np.isfinite(mx2_ep)
    base &= mx2_ep >= args.mx2_ep_min
    base &= mx2_ep < args.mx2_ep_max
    # Mx2(epgamma_tag) is deliberately NOT cut here.  Keep it as a
    # diagnostic variable only; the aggressive purity selection is the BDT.
    base &= np.isfinite(branch_missing_mass2)
    base &= np.isfinite(pred["energy"])
    base &= pred["energy"] >= args.probe_energy_min
    base &= pred["energy"] < args.probe_energy_max
    base &= region >= 0

    tag_parent = None

    if "gamma_parent_pid" in arrays:
        tag_parent = np.asarray(arrays["gamma_parent_pid"], dtype=np.int64)[base]
    #endif

    sample = DirectedSample(
        predicted_energy=pred["energy"][base],
        predicted_theta=pred["theta"][base],
        tag_score=score[base],
        tag_energy=tag_energy[base],
        mx2_ep=mx2_ep[base],
        missing_mass2=branch_missing_mass2[base],
        region=region[base],
        tag_parent_pid=tag_parent,
    )

    diagnostics = {
        "all_score": score,
        "all_tag_energy": tag_energy,
        "all_mx2_ep": mx2_ep,
        "all_missing_mass2": branch_missing_mass2,
        "all_pred_energy": pred["energy"],
        "all_pred_theta": pred["theta"],
        "all_pred_phi": pred["phi"],
        "selected_mask": base,
        # Compare the branch Mx2(epgamma) against a direct four-vector check.
        "computed_missing_mass2": pred["mass2"],
    }

    # Keep truth labels on the unfiltered denominator population so AAOgen
    # can be restricted to genuinely truth-matched reconstructed photons.
    if "matching_gamma_pid" in arrays:
        diagnostics["all_matching_gamma_pid"] = np.asarray(
            arrays["matching_gamma_pid"],
            dtype=np.int64,
        )
    #endif
    if "gamma_mcindex" in arrays:
        diagnostics["all_gamma_mcindex"] = np.asarray(
            arrays["gamma_mcindex"],
            dtype=np.int64,
        )
    #endif
    if "gamma_parent_pid" in arrays:
        diagnostics["all_gamma_parent_pid"] = np.asarray(
            arrays["gamma_parent_pid"],
            dtype=np.int64,
        )
    #endif

    return sample, diagnostics


def score_epgg_tag_direction(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, object],
    tag_index: int,
) -> np.ndarray:
    X = feature_matrix_epgg(arrays, tag_index)

    detector_branch = (
        "detector_gamma1"
        if tag_index == 1
        else "detector_gamma2"
    )
    detector = np.asarray(arrays[detector_branch], dtype=np.int64)

    valid = np.all(np.isfinite(X), axis=1)
    valid &= X[:, FEATURES.index("open_angle_ep2")] > OPEN_ANGLE_MIN
    valid &= X[:, FEATURES.index("Mx2")] > LOOSE_MX2_MIN
    valid &= X[:, FEATURES.index("Mx2")] < LOOSE_MX2_MAX

    scores = np.full(len(X), np.nan, dtype=float)

    for region_name, detector_value in REGIONS.items():
        mask = valid & (detector == detector_value)

        if np.any(mask):
            scores[mask] = score_model(models[region_name], X[mask])
        #endif
    #endfor

    return scores


def numerator_direction_from_epgg(
    arrays: Mapping[str, np.ndarray],
    models: Mapping[str, object],
    beam_energy: float,
    args: argparse.Namespace,
    tag_index: int,
) -> Tuple[DirectedSample, Dict[str, np.ndarray]]:
    if tag_index == 1:
        tag_p = "p2_p"
        tag_theta = "p2_theta"
        tag_phi = "p2_phi"
        probe_p = "p3_p"
        probe_theta = "p3_theta"
        probe_phi = "p3_phi"
        probe_detector = "detector_gamma2"
        tag_parent_branch = "gamma1_parent_pid"
        tag_pid_branch = "matching_gamma1_pid"
        tag_mcindex_branch = "gamma1_mcindex"
        probe_pid_branch = "matching_gamma2_pid"
        probe_mcindex_branch = "gamma2_mcindex"
        probe_parent_branch = "gamma2_parent_pid"
    elif tag_index == 2:
        tag_p = "p3_p"
        tag_theta = "p3_theta"
        tag_phi = "p3_phi"
        probe_p = "p2_p"
        probe_theta = "p2_theta"
        probe_phi = "p2_phi"
        probe_detector = "detector_gamma1"
        tag_parent_branch = "gamma2_parent_pid"
        tag_pid_branch = "matching_gamma2_pid"
        tag_mcindex_branch = "gamma2_mcindex"
        probe_pid_branch = "matching_gamma1_pid"
        probe_mcindex_branch = "gamma1_mcindex"
        probe_parent_branch = "gamma1_parent_pid"
    else:
        raise ValueError("tag_index must be 1 or 2")
    #endif

    X = feature_matrix_epgg(arrays, tag_index)
    scores = score_epgg_tag_direction(
        arrays,
        models,
        tag_index,
    )
    pred = predicted_probe_from_tag(
        arrays,
        beam_energy,
        tag_p,
        tag_theta,
        tag_phi,
    )

    actual_px, actual_py, actual_pz = spherical_to_cartesian(
        arrays[probe_p],
        arrays[probe_theta],
        arrays[probe_phi],
    )
    match_angle = opening_angle_deg(
        pred["px"],
        pred["py"],
        pred["pz"],
        actual_px,
        actual_py,
        actual_pz,
    )

    mx2_ep = X[:, FEATURES.index("Mx2_1")]
    branch_missing_mass2 = X[:, FEATURES.index("Mx2")]
    tag_energy = np.asarray(arrays[tag_p], dtype=float)
    partner_energy = np.asarray(arrays[probe_p], dtype=float)
    actual_detector = np.asarray(arrays[probe_detector], dtype=np.int64)
    predicted_region = assign_region(pred["theta"], args)

    # Build the numerator cut flow explicitly.  The nominal purity selection
    # is intentionally dominated by the BDT score; Mx2(epgamma_tag) is only
    # retained as a diagnostic.
    masks: Dict[str, np.ndarray] = {}

    masks["all"] = np.ones(len(scores), dtype=bool)
    masks["finite_bdt"] = np.isfinite(scores)
    masks["score"] = masks["finite_bdt"] & (scores >= args.tag_score_min)
    masks["tag_energy"] = masks["score"] & (tag_energy >= args.tag_energy_min)
    masks["mx2_ep"] = (
        masks["tag_energy"]
        & np.isfinite(mx2_ep)
        & (mx2_ep >= args.mx2_ep_min)
        & (mx2_ep < args.mx2_ep_max)
    )
    masks["probe_energy"] = (
        masks["mx2_ep"]
        & np.isfinite(pred["energy"])
        & (pred["energy"] >= args.probe_energy_min)
        & (pred["energy"] < args.probe_energy_max)
    )
    # This is the last stage BEFORE any predicted angular-acceptance
    # requirement.  Geometry-debug plots use it deliberately so that a bad
    # predicted direction cannot preselect the diagnostic sample.
    masks["pre_acceptance"] = masks["probe_energy"]

    masks["probe_acceptance"] = (
        masks["pre_acceptance"]
        & (predicted_region >= 0)
    )
    masks["finite_match"] = (
        masks["probe_acceptance"]
        & np.isfinite(match_angle)
    )
    masks["same_detector"] = (
        masks["finite_match"]
        & (actual_detector == predicted_region)
    )
    masks["partner_energy"] = (
        masks["same_detector"]
        & (partner_energy >= args.probe_energy_min)
    )
    masks["angular_match"] = (
        masks["partner_energy"]
        & (match_angle < args.probe_match_angle_max)
    )

    base = masks["angular_match"]

    tag_parent = None

    if tag_parent_branch in arrays:
        tag_parent = np.asarray(
            arrays[tag_parent_branch],
            dtype=np.int64,
        )[base]
    #endif

    sample = DirectedSample(
        predicted_energy=pred["energy"][base],
        predicted_theta=pred["theta"][base],
        tag_score=scores[base],
        tag_energy=tag_energy[base],
        mx2_ep=mx2_ep[base],
        missing_mass2=branch_missing_mass2[base],
        region=predicted_region[base],
        partner_match_angle=match_angle[base],
        partner_energy=partner_energy[base],
        pair_mass=np.asarray(arrays["Mh_gammagamma"], dtype=float)[base],
        tag_parent_pid=tag_parent,
    )

    reco_probe_theta = angle_to_degrees(
        np.asarray(arrays[probe_theta], dtype=float)
    )
    reco_probe_phi = angle_to_degrees(
        np.asarray(arrays[probe_phi], dtype=float)
    )

    diagnostics = {
        "tag_index": np.asarray([tag_index], dtype=np.int64),
        "scores": scores,
        "tag_energy": tag_energy,
        "mx2_ep": mx2_ep,
        "missing_mass2": branch_missing_mass2,
        "pred_energy": pred["energy"],
        "pred_theta": pred["theta"],
        "pred_phi": pred["phi"],
        "pred_region": predicted_region,
        "actual_detector": actual_detector,
        "partner_energy": partner_energy,
        "partner_theta": reco_probe_theta,
        "partner_phi": reco_probe_phi,
        "match_angle": match_angle,
        "pair_mass": np.asarray(arrays["Mh_gammagamma"], dtype=float),
    }

    truth_branches = {
        "tag_truth_pid": tag_pid_branch,
        "tag_truth_mcindex": tag_mcindex_branch,
        "tag_truth_parent_pid": tag_parent_branch,
        "probe_truth_pid": probe_pid_branch,
        "probe_truth_mcindex": probe_mcindex_branch,
        "probe_truth_parent_pid": probe_parent_branch,
    }

    for output_name, branch_name in truth_branches.items():
        if branch_name in arrays:
            diagnostics[output_name] = np.asarray(
                arrays[branch_name],
                dtype=np.int64,
            )
        #endif
    #endfor

    for name, mask in masks.items():
        diagnostics[f"mask_{name}"] = mask
    #endfor

    return sample, diagnostics


def print_numerator_cutflow(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
) -> None:
    tag_index = int(diagnostics["tag_index"][0])
    ordered = [
        ("all", "all directed pairs"),
        ("finite_bdt", "finite BDT / training preselection"),
        ("score", "tag BDT score cut"),
        ("tag_energy", "tag energy cut"),
        ("mx2_ep", "Mx2(ep) cleanup"),
        ("probe_energy", "physical probe-energy range"),
        ("pre_acceptance", "pre-acceptance geometry sample"),
        ("probe_acceptance", "predicted FT/FD acceptance"),
        ("finite_match", "finite predicted/reco angle"),
        ("same_detector", "reco probe in predicted detector"),
        ("partner_energy", "reco probe above threshold"),
        ("angular_match", "angular match"),
    ]

    progress(
        f"{source_name.upper()} EPGG tag-direction {tag_index} cut flow:"
    )

    previous = None

    for key, label in ordered:
        count = int(np.sum(diagnostics[f"mask_{key}"]))

        if previous is None or previous == 0:
            frac = 1.0 if previous is None else 0.0
        else:
            frac = count / previous
        #endif

        print(
            f"    {label:38s}: {count:10,d}"
            + ("" if previous is None else f"  ({100.0*frac:6.2f}% of previous)")
        )
        previous = count
    #endfor


def concatenate_numerator_diagnostics(
    first: Mapping[str, np.ndarray],
    second: Mapping[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    keys = [
        "scores",
        "tag_energy",
        "mx2_ep",
        "missing_mass2",
        "pred_energy",
        "pred_theta",
        "pred_phi",
        "pred_region",
        "actual_detector",
        "partner_energy",
        "partner_theta",
        "partner_phi",
        "match_angle",
        "pair_mass",
    ]

    optional_keys = [
        "tag_truth_pid",
        "tag_truth_mcindex",
        "tag_truth_parent_pid",
        "probe_truth_pid",
        "probe_truth_mcindex",
        "probe_truth_parent_pid",
    ]

    for key in optional_keys:
        if key in first and key in second:
            keys.append(key)
        #endif
    #endfor
    out = {
        key: np.concatenate([first[key], second[key]])
        for key in keys
    }

    mask_keys = [
        key
        for key in first
        if key.startswith("mask_") and key in second
    ]

    for key in mask_keys:
        out[key] = np.concatenate([first[key], second[key]])
    #endfor

    return out


def concatenate_directed(
    first: DirectedSample,
    second: DirectedSample,
) -> DirectedSample:
    def cat_optional(
        a: Optional[np.ndarray],
        b: Optional[np.ndarray],
    ) -> Optional[np.ndarray]:
        if a is None or b is None:
            return None
        #endif
        return np.concatenate([a, b])

    return DirectedSample(
        predicted_energy=np.concatenate(
            [first.predicted_energy, second.predicted_energy]
        ),
        predicted_theta=np.concatenate(
            [first.predicted_theta, second.predicted_theta]
        ),
        tag_score=np.concatenate([first.tag_score, second.tag_score]),
        tag_energy=np.concatenate([first.tag_energy, second.tag_energy]),
        mx2_ep=np.concatenate([first.mx2_ep, second.mx2_ep]),
        missing_mass2=np.concatenate(
            [first.missing_mass2, second.missing_mass2]
        ),
        region=np.concatenate([first.region, second.region]),
        partner_match_angle=cat_optional(
            first.partner_match_angle,
            second.partner_match_angle,
        ),
        partner_energy=cat_optional(
            first.partner_energy,
            second.partner_energy,
        ),
        pair_mass=cat_optional(first.pair_mass, second.pair_mass),
        tag_parent_pid=cat_optional(
            first.tag_parent_pid,
            second.tag_parent_pid,
        ),
    )



# =============================================================================
# Denominator-closure and association-plateau validation
# =============================================================================

def truth_pair_mask_for_source(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
) -> Optional[np.ndarray]:
    """
    Candidate-level truth qualification for the reconstructed sister photon.

    AAOgen:
        The generator is exclusive ep pi0, so two distinct truth-matched
        photons are the two pi0 daughters.

    CLASDIS:
        Require both reconstructed photons to be truth-matched photons whose
        immediate parent PID is 111, plus distinct MC::Particle indices.
        The current ROOT schema does not store the parent row/index, so this
        cannot prove that both photons came from the exact same pi0 object.
        It is therefore a truth-enriched CLASDIS diagnostic, not an exact
        same-parent truth label.

    Data:
        No truth mask is available.
    """
    required = [
        "tag_truth_pid",
        "tag_truth_mcindex",
        "probe_truth_pid",
        "probe_truth_mcindex",
    ]

    if any(key not in diagnostics for key in required):
        return None
    #endif

    tag_pid = np.asarray(diagnostics["tag_truth_pid"], dtype=np.int64)
    probe_pid = np.asarray(diagnostics["probe_truth_pid"], dtype=np.int64)
    tag_mcindex = np.asarray(
        diagnostics["tag_truth_mcindex"],
        dtype=np.int64,
    )
    probe_mcindex = np.asarray(
        diagnostics["probe_truth_mcindex"],
        dtype=np.int64,
    )

    mask = (
        (tag_pid == 22)
        & (probe_pid == 22)
        & (tag_mcindex >= 0)
        & (probe_mcindex >= 0)
        & (tag_mcindex != probe_mcindex)
    )

    if source_name == "clasdis":
        if (
            "tag_truth_parent_pid" not in diagnostics
            or "probe_truth_parent_pid" not in diagnostics
        ):
            return None
        #endif

        tag_parent = np.asarray(
            diagnostics["tag_truth_parent_pid"],
            dtype=np.int64,
        )
        probe_parent = np.asarray(
            diagnostics["probe_truth_parent_pid"],
            dtype=np.int64,
        )
        mask &= (tag_parent == 111) & (probe_parent == 111)
    #endif

    return mask


def reconstructed_partner_region(
    diagnostics: Mapping[str, np.ndarray],
    args: argparse.Namespace,
) -> np.ndarray:
    theta = np.asarray(diagnostics["partner_theta"], dtype=float)
    return assign_region(theta, args)


def denominator_closure_base_mask(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    args: argparse.Namespace,
) -> Optional[np.ndarray]:
    truth_mask = truth_pair_mask_for_source(
        source_name,
        diagnostics,
    )

    if truth_mask is None:
        return None
    #endif

    # Use the tag-side denominator selection through physical predicted probe
    # energy, but DO NOT require predicted angular acceptance.  Otherwise the
    # FT/FD migration test would preselect away failures of the prediction.
    base = (
        np.asarray(diagnostics["mask_probe_energy"], dtype=bool)
        & truth_mask
        & np.isfinite(diagnostics["partner_energy"])
        & np.isfinite(diagnostics["partner_theta"])
        & (
            np.asarray(diagnostics["partner_energy"], dtype=float)
            >= args.probe_energy_min
        )
    )

    return base


def plot_denominator_closure_migrations(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Compare missing-vector probe assignment to the reconstructed truth-qualified
    sister photon.  This is the strongest closure test available with the
    current ROOT schema because generated sister four-vectors are not stored.
    """
    base = denominator_closure_base_mask(
        source_name.lower(),
        diagnostics,
        args,
    )

    if base is None or np.sum(base) == 0:
        progress(
            f"SKIP: {source_name} denominator closure -- "
            "truth-qualified sister information unavailable"
        )
        return
    #endif

    pred_e = np.asarray(diagnostics["pred_energy"], dtype=float)[base]
    reco_e = np.asarray(diagnostics["partner_energy"], dtype=float)[base]
    pred_theta = np.asarray(diagnostics["pred_theta"], dtype=float)[base]
    reco_theta = np.asarray(diagnostics["partner_theta"], dtype=float)[base]

    pred_region = assign_region(pred_theta, args)
    reco_region = assign_region(reco_theta, args)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 10.0))

    h = axes[0, 0].hist2d(
        reco_e,
        pred_e,
        bins=[
            np.linspace(args.probe_energy_min, args.probe_energy_max, 70),
            np.linspace(args.probe_energy_min, args.probe_energy_max, 70),
        ],
        cmin=1,
    )
    axes[0, 0].plot(
        [args.probe_energy_min, args.probe_energy_max],
        [args.probe_energy_min, args.probe_energy_max],
        linestyle="--",
        linewidth=1.0,
    )
    axes[0, 0].set_xlabel(r"Reconstructed sister $E_\gamma$ (GeV)")
    axes[0, 0].set_ylabel(r"Predicted missing-vector $E_\gamma$ (GeV)")
    axes[0, 0].set_title("Probe-energy closure")
    fig.colorbar(h[3], ax=axes[0, 0], label="Candidates")

    h = axes[0, 1].hist2d(
        reco_theta,
        pred_theta,
        bins=[
            np.linspace(2.0, 35.0, 70),
            np.linspace(2.0, 35.0, 70),
        ],
        cmin=1,
    )
    axes[0, 1].plot(
        [2.0, 35.0],
        [2.0, 35.0],
        linestyle="--",
        linewidth=1.0,
    )
    axes[0, 1].set_xlabel(
        r"Reconstructed sister $\theta_\gamma$ (deg)"
    )
    axes[0, 1].set_ylabel(
        r"Predicted missing-vector $\theta_\gamma$ (deg)"
    )
    axes[0, 1].set_title("Probe-angle closure")
    fig.colorbar(h[3], ax=axes[0, 1], label="Candidates")

    # Region migration.  Include outside acceptance explicitly so a bad
    # prediction cannot disappear from the matrix.
    region_labels = ["FT", "FD", "outside"]
    pred_region3 = np.where(pred_region < 0, 2, pred_region)
    reco_region3 = np.where(reco_region < 0, 2, reco_region)
    region_matrix = np.zeros((3, 3), dtype=float)

    for true_bin in range(3):
        denom = np.sum(reco_region3 == true_bin)

        for pred_bin in range(3):
            if denom > 0:
                region_matrix[pred_bin, true_bin] = (
                    np.sum(
                        (pred_region3 == pred_bin)
                        & (reco_region3 == true_bin)
                    )
                    / denom
                )
            #endif
        #endfor
    #endfor

    im = axes[1, 0].imshow(
        region_matrix,
        origin="lower",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
    )
    axes[1, 0].set_xticks(range(3), region_labels)
    axes[1, 0].set_yticks(range(3), region_labels)
    axes[1, 0].set_xlabel("Reconstructed sister region")
    axes[1, 0].set_ylabel("Predicted region")
    axes[1, 0].set_title("FT/FD migration; columns normalized")

    for iy in range(3):
        for ix in range(3):
            axes[1, 0].text(
                ix,
                iy,
                f"{100.0*region_matrix[iy, ix]:.1f}%",
                ha="center",
                va="center",
            )
        #endfor
    #endfor

    fig.colorbar(im, ax=axes[1, 0], label="Fraction")

    # Energy-bin migration in the exact nominal broad bins.
    energy_edges = np.asarray(
        [args.probe_energy_min, 1.0, 2.0, args.probe_energy_max],
        dtype=float,
    )
    pred_bin = np.digitize(pred_e, energy_edges) - 1
    reco_bin = np.digitize(reco_e, energy_edges) - 1
    nbin = len(energy_edges) - 1
    energy_matrix = np.zeros((nbin, nbin), dtype=float)

    for true_bin in range(nbin):
        denom = np.sum(reco_bin == true_bin)

        for predicted_bin in range(nbin):
            if denom > 0:
                energy_matrix[predicted_bin, true_bin] = (
                    np.sum(
                        (pred_bin == predicted_bin)
                        & (reco_bin == true_bin)
                    )
                    / denom
                )
            #endif
        #endfor
    #endfor

    labels = [
        f"{energy_edges[i]:g}-{energy_edges[i+1]:g}"
        for i in range(nbin)
    ]
    im = axes[1, 1].imshow(
        energy_matrix,
        origin="lower",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
    )
    axes[1, 1].set_xticks(range(nbin), labels)
    axes[1, 1].set_yticks(range(nbin), labels)
    axes[1, 1].set_xlabel("Reconstructed sister E bin (GeV)")
    axes[1, 1].set_ylabel("Predicted E bin (GeV)")
    axes[1, 1].set_title("Energy-bin migration; columns normalized")

    for iy in range(nbin):
        for ix in range(nbin):
            axes[1, 1].text(
                ix,
                iy,
                f"{100.0*energy_matrix[iy, ix]:.1f}%",
                ha="center",
                va="center",
            )
        #endfor
    #endfor

    fig.colorbar(im, ax=axes[1, 1], label="Fraction")

    fig.suptitle(
        f"{source_name}: denominator probe-assignment closure",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)


def plot_denominator_closure_resolution(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Resolution summaries versus the reconstructed truth-qualified sister.
    """
    base = denominator_closure_base_mask(
        source_name.lower(),
        diagnostics,
        args,
    )

    if base is None or np.sum(base) == 0:
        return
    #endif

    pred_e = np.asarray(diagnostics["pred_energy"], dtype=float)[base]
    reco_e = np.asarray(diagnostics["partner_energy"], dtype=float)[base]
    pred_theta = np.asarray(diagnostics["pred_theta"], dtype=float)[base]
    reco_theta = np.asarray(diagnostics["partner_theta"], dtype=float)[base]

    frac_de = np.divide(
        pred_e - reco_e,
        reco_e,
        out=np.full_like(pred_e, np.nan),
        where=reco_e > 0.0,
    )
    dtheta = pred_theta - reco_theta

    energy_edges = np.asarray(
        [0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 2.00, 2.50, 3.50, 5.00, 7.00, 9.50],
        dtype=float,
    )
    theta_edges = np.asarray(
        [2.0, 3.0, 4.0, 5.5, 7.0, 9.0, 12.0, 16.0, 21.0, 27.0, 35.0],
        dtype=float,
    )

    def summarize(
        x: np.ndarray,
        y: np.ndarray,
        edges: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        centers = 0.5 * (edges[:-1] + edges[1:])
        median = np.full(len(centers), np.nan)
        half68 = np.full(len(centers), np.nan)

        for ibin in range(len(centers)):
            mask = (
                np.isfinite(x)
                & np.isfinite(y)
                & (x >= edges[ibin])
                & (x < edges[ibin + 1])
            )
            values = y[mask]

            if len(values) < 20:
                continue
            #endif

            q16, q50, q84 = np.quantile(
                values,
                [0.16, 0.50, 0.84],
            )
            median[ibin] = q50
            half68[ibin] = 0.5 * (q84 - q16)
        #endfor

        return centers, median, half68

    e_centers, e_med, e_width = summarize(
        reco_e,
        frac_de,
        energy_edges,
    )
    t_centers, t_med, t_width = summarize(
        reco_theta,
        dtheta,
        theta_edges,
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))

    axes[0].errorbar(
        e_centers,
        e_med,
        yerr=e_width,
        marker="o",
        linestyle="-",
        capsize=2,
    )
    axes[0].axhline(0.0, linestyle="--", linewidth=1.0)
    axes[0].set_xlabel(r"Reconstructed sister $E_\gamma$ (GeV)")
    axes[0].set_ylabel(
        r"$(E_{\rm pred}-E_{\rm reco})/E_{\rm reco}$"
    )
    axes[0].set_title("Median and central 68% half-width")

    axes[1].errorbar(
        t_centers,
        t_med,
        yerr=t_width,
        marker="o",
        linestyle="-",
        capsize=2,
    )
    axes[1].axhline(0.0, linestyle="--", linewidth=1.0)
    axes[1].set_xlabel(
        r"Reconstructed sister $\theta_\gamma$ (deg)"
    )
    axes[1].set_ylabel(
        r"$\theta_{\rm pred}-\theta_{\rm reco}$ (deg)"
    )
    axes[1].set_title("Median and central 68% half-width")

    for ax in axes:
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: denominator prediction resolution",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    save_figure(fig, output_path)


def association_plateau_curves(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    region_value: int,
    energy_range: Tuple[float, float],
) -> Optional[Dict[str, np.ndarray]]:
    truth_mask = truth_pair_mask_for_source(
        source_name,
        diagnostics,
    )

    if truth_mask is None:
        return None
    #endif

    pred_e = np.asarray(diagnostics["pred_energy"], dtype=float)
    pred_region = np.asarray(diagnostics["pred_region"], dtype=np.int64)
    match_angle = np.asarray(diagnostics["match_angle"], dtype=float)

    # Reconstructed partner already exists by construction of epgg.  Start
    # from the tag-side denominator-like selection and require the partner to
    # be above threshold.  Do not impose the nominal angular cut.
    base = (
        np.asarray(diagnostics["mask_partner_energy"], dtype=bool)
        & (pred_region == region_value)
        & (pred_e >= energy_range[0])
        & (pred_e < energy_range[1])
        & np.isfinite(match_angle)
    )

    correct = base & truth_mask
    wrong = base & ~truth_mask

    n_correct_total = int(np.sum(correct))

    recovery = np.full(len(ASSOCIATION_ANGLE_SCAN), np.nan)
    wrong_fraction = np.full(len(ASSOCIATION_ANGLE_SCAN), np.nan)

    for iangle, angle_max in enumerate(ASSOCIATION_ANGLE_SCAN):
        correct_pass = np.sum(
            correct & (match_angle < angle_max)
        )
        wrong_pass = np.sum(
            wrong & (match_angle < angle_max)
        )
        all_pass = correct_pass + wrong_pass

        if n_correct_total > 0:
            recovery[iangle] = correct_pass / n_correct_total
        #endif

        if all_pass > 0:
            wrong_fraction[iangle] = wrong_pass / all_pass
        #endif
    #endfor

    return {
        "angle": ASSOCIATION_ANGLE_SCAN.copy(),
        "recovery": recovery,
        "wrong_fraction": wrong_fraction,
        "n_correct_total": np.asarray([n_correct_total]),
    }


def plot_association_plateau(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Correct-pair recovery and wrong-candidate fraction versus angular window.

    For AAOgen this is a clean candidate-level truth study because the
    generator contains the exclusive pi0 -> gamma gamma final state.

    For CLASDIS the current ancestry schema only proves both photons have
    parent PID 111, not that they share the same parent row.  Treat the
    resulting wrong-fraction estimate as diagnostic.
    """
    if truth_pair_mask_for_source(
        source_name.lower(),
        diagnostics,
    ) is None:
        progress(
            f"SKIP: {source_name} association plateau -- truth unavailable"
        )
        return
    #endif

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.5))
    region_names = {0: "FT", 1: "FD"}

    for icol, region_value in enumerate((0, 1)):
        for energy_range in ASSOCIATION_ENERGY_BINS:
            result = association_plateau_curves(
                source_name.lower(),
                diagnostics,
                args,
                region_value,
                energy_range,
            )

            if result is None:
                continue
            #endif

            label = (
                f"{energy_range[0]:g}-{energy_range[1]:g} GeV"
            )

            axes[0, icol].plot(
                result["angle"],
                result["recovery"],
                marker="o",
                label=label,
            )
            axes[1, icol].plot(
                result["angle"],
                result["wrong_fraction"],
                marker="o",
                label=label,
            )
        #endfor

        axes[0, icol].axvline(
            args.probe_match_angle_max,
            linestyle="--",
            linewidth=1.0,
            label=f"nominal {args.probe_match_angle_max:g} deg",
        )
        axes[1, icol].axvline(
            args.probe_match_angle_max,
            linestyle="--",
            linewidth=1.0,
        )

        axes[0, icol].set_title(
            f"{region_names[region_value]}: correct-pair recovery"
        )
        axes[0, icol].set_xlabel(
            r"Association window $\Delta\Omega_{\max}$ (deg)"
        )
        axes[0, icol].set_ylabel(
            "Fraction of truth-qualified pairs recovered"
        )
        axes[0, icol].set_ylim(0.0, 1.05)

        axes[1, icol].set_title(
            f"{region_names[region_value]}: wrong-candidate fraction"
        )
        axes[1, icol].set_xlabel(
            r"Association window $\Delta\Omega_{\max}$ (deg)"
        )
        axes[1, icol].set_ylabel(
            "Wrong candidates / accepted candidates"
        )
        axes[1, icol].set_ylim(0.0, 1.0)
    #endfor

    for ax in axes.flat:
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    #endfor

    fig.suptitle(
        f"{source_name}: angular-association plateau study",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_mgg_angle_scan(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Show how the reconstructed pi0 peak and combinatorial continuum evolve as
    the angular association window is loosened.  This is deliberately NOT the
    final numerator extraction; the final numerator should come from fits to
    these Mgg spectra.
    """
    pred_region = np.asarray(
        diagnostics["pred_region"],
        dtype=np.int64,
    )
    match_angle = np.asarray(
        diagnostics["match_angle"],
        dtype=float,
    )
    pair_mass = np.asarray(
        diagnostics["pair_mass"],
        dtype=float,
    )
    base = (
        np.asarray(diagnostics["mask_partner_energy"], dtype=bool)
        & np.isfinite(match_angle)
        & np.isfinite(pair_mass)
    )

    scan_angles = [3.0, 5.0, 8.0, 12.0]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    region_names = {0: "FT", 1: "FD"}

    for iax, region_value in enumerate((0, 1)):
        region_base = base & (pred_region == region_value)

        for angle_max in scan_angles:
            mask = region_base & (match_angle < angle_max)
            masses = pair_mass[mask]

            axes[iax].hist(
                masses,
                bins=np.linspace(0.0, 0.40, 121),
                histtype="step",
                linewidth=1.2,
                label=f"< {angle_max:g} deg (N={len(masses):,})",
            )
        #endfor

        axes[iax].axvspan(
            PI0_MASS_WINDOW[0],
            PI0_MASS_WINDOW[1],
            alpha=0.08,
        )
        axes[iax].set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
        axes[iax].set_ylabel("Candidate pairs")
        axes[iax].set_yscale("log")
        axes[iax].set_title(
            f"{region_names[region_value]}: pair-mass evolution"
        )
        axes[iax].grid(alpha=0.2)
        axes[iax].legend(fontsize=8)
    #endfor

    fig.suptitle(
        f"{source_name}: pi0 mass peak versus association window",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    save_figure(fig, output_path)




# =============================================================================
# M(gamma gamma) fit numerator and probe-energy response
# =============================================================================

def mgg_model(
    x: np.ndarray,
    amplitude_core: float,
    amplitude_tail: float,
    mean: float,
    sigma_core: float,
    delta_sigma: float,
    b0: float,
    b1: float,
    b2: float,
) -> np.ndarray:
    """
    Common-mean double Gaussian pi0 signal plus quadratic continuum.

    sigma_tail is parameterized as sigma_core + delta_sigma so the broad
    component is guaranteed to remain broader than the core.
    """
    sigma_tail = sigma_core + delta_sigma
    core = amplitude_core * np.exp(
        -0.5 * ((x - mean) / sigma_core) ** 2
    )
    tail = amplitude_tail * np.exp(
        -0.5 * ((x - mean) / sigma_tail) ** 2
    )
    dx = x - 0.135
    background = b0 + b1 * dx + b2 * dx * dx
    return core + tail + background


def fit_mgg_spectrum(
    masses: np.ndarray,
) -> MassFitResult:
    masses = np.asarray(masses, dtype=float)
    masses = masses[
        np.isfinite(masses)
        & (masses >= MGG_FIT_RANGE[0])
        & (masses < MGG_FIT_RANGE[1])
    ]

    hist_edges = np.linspace(
        MGG_FIT_RANGE[0],
        MGG_FIT_RANGE[1],
        MGG_FIT_BINS + 1,
    )
    counts, _ = np.histogram(masses, bins=hist_edges)
    counts = counts.astype(float)
    centers = 0.5 * (hist_edges[:-1] + hist_edges[1:])
    bin_width = hist_edges[1] - hist_edges[0]

    def failed_result() -> MassFitResult:
        return MassFitResult(
            success=False,
            yield_pi0=np.nan,
            yield_uncertainty=np.nan,
            mean=np.nan,
            sigma=np.nan,
            sigma_core=np.nan,
            sigma_tail=np.nan,
            core_fraction=np.nan,
            chi2_ndf=np.nan,
            x_centers=centers,
            counts=counts,
            model_counts=np.full_like(centers, np.nan),
            signal_counts=np.full_like(centers, np.nan),
            background_counts=np.full_like(centers, np.nan),
        )

    if len(masses) < 100 or np.sum(counts) <= 0:
        return failed_result()
    #endif

    sideband = (
        (centers < 0.110)
        | (centers > 0.165)
    )
    background_guess = (
        float(np.median(counts[sideband]))
        if np.any(sideband)
        else float(np.median(counts))
    )
    peak_region = (
        (centers >= 0.120)
        & (centers <= 0.150)
    )
    peak_guess = max(
        float(np.max(counts[peak_region])) - background_guess,
        1.0,
    )

    # Start with a dominant narrow core and a smaller broad component.
    p0 = [
        0.75 * peak_guess,       # core amplitude
        0.25 * peak_guess,       # tail amplitude
        0.135,                   # common mean
        0.007,                   # core width
        0.008,                   # tail-core width difference
        max(background_guess, 0.0),
        0.0,
        0.0,
    ]

    amplitude_upper = max(10.0 * np.max(counts), 10.0)

    lower = [
        0.0,        # core amplitude
        0.0,        # tail amplitude
        0.125,      # mean
        0.0025,     # core width
        0.0005,     # tail-core separation
        0.0,        # background constant
        -1.0e7,
        -1.0e8,
    ]
    upper = [
        amplitude_upper,
        amplitude_upper,
        0.145,
        0.018,      # core width
        0.035,      # additional tail width
        amplitude_upper,
        1.0e7,
        1.0e8,
    ]

    sigma_y = np.sqrt(np.maximum(counts, 1.0))

    try:
        popt, pcov = curve_fit(
            mgg_model,
            centers,
            counts,
            p0=p0,
            sigma=sigma_y,
            absolute_sigma=True,
            bounds=(lower, upper),
            maxfev=100000,
        )

        (
            amplitude_core,
            amplitude_tail,
            mean,
            sigma_core,
            delta_sigma,
            b0,
            b1,
            b2,
        ) = popt
        sigma_tail = sigma_core + delta_sigma

        model_counts = mgg_model(centers, *popt)
        core_counts = amplitude_core * np.exp(
            -0.5 * ((centers - mean) / sigma_core) ** 2
        )
        tail_counts = amplitude_tail * np.exp(
            -0.5 * ((centers - mean) / sigma_tail) ** 2
        )
        signal_counts = core_counts + tail_counts

        dx = centers - 0.135
        background_counts = b0 + b1 * dx + b2 * dx * dx

        norm = np.sqrt(2.0 * np.pi) / bin_width
        core_yield = amplitude_core * sigma_core * norm
        tail_yield = amplitude_tail * sigma_tail * norm
        yield_pi0 = core_yield + tail_yield

        if yield_pi0 > 0.0:
            core_fraction = core_yield / yield_pi0
            sigma_effective = np.sqrt(
                (
                    core_yield * sigma_core * sigma_core
                    + tail_yield * sigma_tail * sigma_tail
                )
                / yield_pi0
            )
        else:
            core_fraction = np.nan
            sigma_effective = np.nan
        #endif

        # Full covariance propagation for
        # Y = norm * [A_c*s_c + A_t*(s_c + delta)].
        gradient = np.zeros(len(popt), dtype=float)
        gradient[0] = norm * sigma_core
        gradient[1] = norm * sigma_tail
        gradient[2] = 0.0
        gradient[3] = norm * (
            amplitude_core + amplitude_tail
        )
        gradient[4] = norm * amplitude_tail
        variance = float(gradient @ pcov @ gradient)
        yield_uncertainty = np.sqrt(max(variance, 0.0))

        residual = (counts - model_counts) / sigma_y
        ndf = max(len(counts) - len(popt), 1)
        chi2_ndf = float(
            np.sum(residual * residual) / ndf
        )

        success = (
            np.isfinite(yield_pi0)
            and yield_pi0 >= 0.0
            and np.isfinite(yield_uncertainty)
            and 0.125 <= mean <= 0.145
            and 0.0025 <= sigma_core <= 0.018
            and sigma_tail > sigma_core
            and sigma_tail <= 0.053
        )

        return MassFitResult(
            success=bool(success),
            yield_pi0=float(yield_pi0),
            yield_uncertainty=float(yield_uncertainty),
            mean=float(mean),
            sigma=float(sigma_effective),
            sigma_core=float(sigma_core),
            sigma_tail=float(sigma_tail),
            core_fraction=float(core_fraction),
            chi2_ndf=chi2_ndf,
            x_centers=centers,
            counts=counts,
            model_counts=model_counts,
            signal_counts=signal_counts,
            background_counts=background_counts,
        )
    except Exception as exc:
        progress(f"WARNING: Mgg fit failed: {exc}")
        return failed_result()
    #endtry



def mgg_fit_candidate_mask(
    diagnostics: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    region_value: int,
    e_low: float,
    e_high: float,
    angle_max: float,
) -> np.ndarray:
    """
    Candidate mask for the fit-based numerator.

    Important:
      * region and energy bin are defined by the PREDICTED probe;
      * no probe BDT requirement;
      * do not require actual_detector == predicted_region;
      * reconstructed partner only needs to exist above threshold and lie
        within the deliberately loose association cone.
    """
    pred_energy = np.asarray(
        diagnostics["pred_energy"],
        dtype=float,
    )
    pred_region = np.asarray(
        diagnostics["pred_region"],
        dtype=np.int64,
    )
    partner_energy = np.asarray(
        diagnostics["partner_energy"],
        dtype=float,
    )
    match_angle = np.asarray(
        diagnostics["match_angle"],
        dtype=float,
    )

    return (
        np.asarray(diagnostics["mask_finite_match"], dtype=bool)
        & np.isfinite(partner_energy)
        & (partner_energy >= args.probe_energy_min)
        & (pred_region == region_value)
        & (pred_energy >= e_low)
        & (pred_energy < e_high)
        & np.isfinite(match_angle)
        & (match_angle < angle_max)
    )


def fit_efficiency_curve(
    denominator: DirectedSample,
    diagnostics: Mapping[str, np.ndarray],
    region_value: int,
    edges: np.ndarray,
    angle_max: float,
    args: argparse.Namespace,
) -> Tuple[FitEfficiencyCurve, list[MassFitResult]]:
    den_values = denominator.predicted_energy[
        denominator.region == region_value
    ]
    den, _ = np.histogram(den_values, bins=edges)
    den = den.astype(float)

    nbin = len(edges) - 1
    numerator_yield = np.full(nbin, np.nan)
    numerator_uncertainty = np.full(nbin, np.nan)
    fit_mean = np.full(nbin, np.nan)
    fit_sigma = np.full(nbin, np.nan)
    fit_chi2_ndf = np.full(nbin, np.nan)
    fits: list[MassFitResult] = []

    pair_mass = np.asarray(
        diagnostics["pair_mass"],
        dtype=float,
    )

    for ibin in range(nbin):
        mask = mgg_fit_candidate_mask(
            diagnostics,
            args,
            region_value,
            edges[ibin],
            edges[ibin + 1],
            angle_max,
        )
        fit = fit_mgg_spectrum(pair_mass[mask])
        fits.append(fit)

        if fit.success:
            numerator_yield[ibin] = fit.yield_pi0
            numerator_uncertainty[ibin] = fit.yield_uncertainty
            fit_mean[ibin] = fit.mean
            fit_sigma[ibin] = fit.sigma
            fit_chi2_ndf[ibin] = fit.chi2_ndf
        #endif
    #endfor

    efficiency = np.divide(
        numerator_yield,
        den,
        out=np.full(nbin, np.nan),
        where=den > 0.0,
    )
    uncertainty = np.divide(
        numerator_uncertainty,
        den,
        out=np.full(nbin, np.nan),
        where=den > 0.0,
    )

    return (
        FitEfficiencyCurve(
            edges=edges,
            centers=0.5 * (edges[:-1] + edges[1:]),
            denominator=den,
            numerator_yield=numerator_yield,
            numerator_uncertainty=numerator_uncertainty,
            efficiency=efficiency,
            uncertainty=uncertainty,
            fit_mean=fit_mean,
            fit_sigma=fit_sigma,
            fit_chi2_ndf=fit_chi2_ndf,
        ),
        fits,
    )



def denominator_counts_for_score_threshold(
    diagnostics: Mapping[str, np.ndarray],
    region_value: int,
    edges: np.ndarray,
    score_threshold: float,
    args: argparse.Namespace,
    require_truth_matched_photon: bool = False,
    theta_low: Optional[float] = None,
    theta_high: Optional[float] = None,
    fd_sector: Optional[int] = None,
) -> np.ndarray:
    """
    Rebuild the epgamma denominator directly from the pre-cut diagnostic
    arrays for an arbitrary BDT threshold.

    This does NOT use the nominal denominator DirectedSample because that
    object has already been filtered at args.tag_score_min and therefore
    cannot support a scan below the nominal score cut.
    """
    score = np.asarray(diagnostics["all_score"], dtype=float)
    tag_energy = np.asarray(
        diagnostics["all_tag_energy"],
        dtype=float,
    )
    mx2_ep = np.asarray(
        diagnostics["all_mx2_ep"],
        dtype=float,
    )
    missing_mass2 = np.asarray(
        diagnostics["all_missing_mass2"],
        dtype=float,
    )
    pred_energy = np.asarray(
        diagnostics["all_pred_energy"],
        dtype=float,
    )
    pred_theta = np.asarray(
        diagnostics["all_pred_theta"],
        dtype=float,
    )
    region = assign_region(pred_theta, args)

    mask = (
        np.isfinite(score)
        & (score >= score_threshold)
        & np.isfinite(tag_energy)
        & (tag_energy >= args.tag_energy_min)
        & np.isfinite(mx2_ep)
        & (mx2_ep >= args.mx2_ep_min)
        & (mx2_ep < args.mx2_ep_max)
        & np.isfinite(missing_mass2)
        & np.isfinite(pred_energy)
        & (pred_energy >= args.probe_energy_min)
        & (pred_energy < args.probe_energy_max)
        & (region == region_value)
    )

    if theta_low is not None:
        mask &= pred_theta >= float(theta_low)
    #endif
    if theta_high is not None:
        mask &= pred_theta < float(theta_high)
    #endif

    if fd_sector is not None:
        pred_phi = np.asarray(
            diagnostics["all_pred_phi"],
            dtype=float,
        )
        mask &= assign_fd_sector(pred_phi) == int(fd_sector)
    #endif

    if require_truth_matched_photon:
        if (
            "all_matching_gamma_pid" not in diagnostics
            or "all_gamma_mcindex" not in diagnostics
        ):
            raise RuntimeError(
                "Truth-matched denominator requested, but epgamma truth "
                "branches are unavailable."
            )
        #endif
        mask &= (
            np.asarray(
                diagnostics["all_matching_gamma_pid"],
                dtype=np.int64,
            ) == 22
        )
        mask &= (
            np.asarray(
                diagnostics["all_gamma_mcindex"],
                dtype=np.int64,
            ) >= 0
        )
    #endif

    counts, _ = np.histogram(
        pred_energy[mask],
        bins=edges,
    )
    return counts.astype(float)


def fit_numerator_for_score_threshold(
    diagnostics: Mapping[str, np.ndarray],
    region_value: int,
    edges: np.ndarray,
    score_threshold: float,
    angle_max: float,
    args: argparse.Namespace,
    require_truth_matched_tag: bool = False,
    theta_low: Optional[float] = None,
    theta_high: Optional[float] = None,
    fd_sector: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, list[MassFitResult]]:
    """
    Fit the pi0 numerator in each E_pred bin at an arbitrary minimum tag score.

    Selection is rebuilt from the unfiltered epgammagamma diagnostics so the
    scan is independent of the nominal args.tag_score_min.
    """
    scores = np.asarray(diagnostics["scores"], dtype=float)
    tag_energy = np.asarray(
        diagnostics["tag_energy"],
        dtype=float,
    )
    mx2_ep = np.asarray(
        diagnostics["mx2_ep"],
        dtype=float,
    )
    pred_energy = np.asarray(
        diagnostics["pred_energy"],
        dtype=float,
    )
    pred_region = np.asarray(
        diagnostics["pred_region"],
        dtype=np.int64,
    )
    partner_energy = np.asarray(
        diagnostics["partner_energy"],
        dtype=float,
    )
    match_angle = np.asarray(
        diagnostics["match_angle"],
        dtype=float,
    )
    pair_mass = np.asarray(
        diagnostics["pair_mass"],
        dtype=float,
    )

    common = (
        np.isfinite(scores)
        & (scores >= score_threshold)
        & np.isfinite(tag_energy)
        & (tag_energy >= args.tag_energy_min)
        & np.isfinite(mx2_ep)
        & (mx2_ep >= args.mx2_ep_min)
        & (mx2_ep < args.mx2_ep_max)
        & np.isfinite(pred_energy)
        & (pred_energy >= args.probe_energy_min)
        & (pred_energy < args.probe_energy_max)
        & (pred_region == region_value)
        & np.isfinite(partner_energy)
        & (partner_energy >= args.probe_energy_min)
        & np.isfinite(match_angle)
        & (match_angle < angle_max)
        & np.isfinite(pair_mass)
    )

    if theta_low is not None:
        common &= np.asarray(
            diagnostics["pred_theta"],
            dtype=float,
        ) >= float(theta_low)
    #endif
    if theta_high is not None:
        common &= np.asarray(
            diagnostics["pred_theta"],
            dtype=float,
        ) < float(theta_high)
    #endif

    if fd_sector is not None:
        pred_phi = np.asarray(
            diagnostics["pred_phi"],
            dtype=float,
        )
        common &= assign_fd_sector(pred_phi) == int(fd_sector)
    #endif

    if require_truth_matched_tag:
        if (
            "tag_truth_pid" not in diagnostics
            or "tag_truth_mcindex" not in diagnostics
        ):
            raise RuntimeError(
                "Truth-matched numerator requested, but epgammagamma truth "
                "branches are unavailable."
            )
        #endif
        common &= (
            np.asarray(
                diagnostics["tag_truth_pid"],
                dtype=np.int64,
            ) == 22
        )
        common &= (
            np.asarray(
                diagnostics["tag_truth_mcindex"],
                dtype=np.int64,
            ) >= 0
        )
    #endif

    nbin = len(edges) - 1
    yields = np.full(nbin, np.nan)
    uncertainties = np.full(nbin, np.nan)
    fits: list[MassFitResult] = []

    for ibin in range(nbin):
        mask = (
            common
            & (pred_energy >= edges[ibin])
            & (pred_energy < edges[ibin + 1])
        )
        fit = fit_mgg_spectrum(pair_mass[mask])
        fits.append(fit)

        if fit.success:
            yields[ibin] = fit.yield_pi0
            uncertainties[ibin] = fit.yield_uncertainty
        #endif
    #endfor

    return yields, uncertainties, fits


def build_bdt_efficiency_threshold_scan(
    denominator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    numerator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    """
    Scan apparent fitted efficiency versus minimum BDT tag score.

    Returned arrays have shape (N_threshold, N_energy_bin).
    """
    results: Dict[
        str,
        Dict[str, Dict[str, np.ndarray]],
    ] = {}

    n_threshold = len(BDT_EFFICIENCY_SCAN_THRESHOLDS)
    n_bin = len(edges) - 1

    for source_name in ["data", "aaogen", "clasdis"]:
        results[source_name] = {}

        for region_name, region_value in REGIONS.items():
            den = np.zeros((n_threshold, n_bin), dtype=float)
            num = np.full((n_threshold, n_bin), np.nan)
            num_unc = np.full((n_threshold, n_bin), np.nan)
            eff = np.full((n_threshold, n_bin), np.nan)
            eff_unc = np.full((n_threshold, n_bin), np.nan)

            for ithr, threshold in enumerate(
                BDT_EFFICIENCY_SCAN_THRESHOLDS
            ):
                den[ithr] = denominator_counts_for_score_threshold(
                    denominator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    float(threshold),
                    args,
                )

                (
                    num[ithr],
                    num_unc[ithr],
                    _,
                ) = fit_numerator_for_score_threshold(
                    numerator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    float(threshold),
                    MGG_FIT_NOMINAL_ANGLE,
                    args,
                )

                eff[ithr] = np.divide(
                    num[ithr],
                    den[ithr],
                    out=np.full(n_bin, np.nan),
                    where=den[ithr] > 0.0,
                )
                eff_unc[ithr] = np.divide(
                    num_unc[ithr],
                    den[ithr],
                    out=np.full(n_bin, np.nan),
                    where=den[ithr] > 0.0,
                )
            #endfor

            results[source_name][region_name] = {
                "denominator": den,
                "numerator": num,
                "numerator_uncertainty": num_unc,
                "efficiency": eff,
                "efficiency_uncertainty": eff_unc,
            }
        #endfor
    #endfor

    return results



def build_aaogen_truth_tag_threshold_scan(
    denominator_diagnostics: Mapping[str, np.ndarray],
    numerator_diagnostics: Mapping[str, np.ndarray],
    edges: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, Dict[str, np.ndarray]]:
    """
    AAOgen-only threshold scan requiring the reconstructed tag photon to be
    truth matched to a generated photon in BOTH denominator and numerator.

    Because AAOgen is exclusive ep pi0, a truth-matched generated photon is
    a pi0 daughter by construction.  This isolates BDT/phase-space selection
    effects from unmatched/nonphoton contamination in the AAOgen tag sample.
    """
    results: Dict[str, Dict[str, np.ndarray]] = {}
    n_threshold = len(BDT_EFFICIENCY_SCAN_THRESHOLDS)
    n_bin = len(edges) - 1

    for region_name, region_value in REGIONS.items():
        den = np.zeros((n_threshold, n_bin), dtype=float)
        num = np.full((n_threshold, n_bin), np.nan)
        num_unc = np.full((n_threshold, n_bin), np.nan)
        eff = np.full((n_threshold, n_bin), np.nan)
        eff_unc = np.full((n_threshold, n_bin), np.nan)

        for ithr, threshold in enumerate(
            BDT_EFFICIENCY_SCAN_THRESHOLDS
        ):
            den[ithr] = denominator_counts_for_score_threshold(
                denominator_diagnostics,
                region_value,
                edges,
                float(threshold),
                args,
                require_truth_matched_photon=True,
            )

            (
                num[ithr],
                num_unc[ithr],
                _,
            ) = fit_numerator_for_score_threshold(
                numerator_diagnostics,
                region_value,
                edges,
                float(threshold),
                MGG_FIT_NOMINAL_ANGLE,
                args,
                require_truth_matched_tag=True,
            )

            eff[ithr] = np.divide(
                num[ithr],
                den[ithr],
                out=np.full(n_bin, np.nan),
                where=den[ithr] > 0.0,
            )
            eff_unc[ithr] = np.divide(
                num_unc[ithr],
                den[ithr],
                out=np.full(n_bin, np.nan),
                where=den[ithr] > 0.0,
            )
        #endfor

        results[region_name] = {
            "denominator": den,
            "numerator": num,
            "numerator_uncertainty": num_unc,
            "efficiency": eff,
            "efficiency_uncertainty": eff_unc,
        }
    #endfor

    return results




def build_theta_binned_bdt_scan(
    denominator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    numerator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, Dict[str, Dict[str, np.ndarray]]]:
    """
    Repeat the apparent-efficiency threshold scan inside coarse probe-theta
    bins.  This controls the dominant phase-space migration induced by the
    BDT score and directly prepares the eventual E_pred x theta_pred
    correction structure.
    """
    results: Dict[str, Dict[str, Dict[str, np.ndarray]]] = {}
    theta_bins = probe_theta_bins(args)
    n_threshold = len(BDT_EFFICIENCY_SCAN_THRESHOLDS)
    n_energy = len(edges) - 1

    for source_name in ["data", "aaogen", "clasdis"]:
        results[source_name] = {}

        for label, region_value, theta_low, theta_high in theta_bins:
            den = np.zeros((n_threshold, n_energy), dtype=float)
            num = np.full((n_threshold, n_energy), np.nan)
            num_unc = np.full((n_threshold, n_energy), np.nan)
            eff = np.full((n_threshold, n_energy), np.nan)
            eff_unc = np.full((n_threshold, n_energy), np.nan)

            for ithr, threshold in enumerate(
                BDT_EFFICIENCY_SCAN_THRESHOLDS
            ):
                den[ithr] = denominator_counts_for_score_threshold(
                    denominator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    float(threshold),
                    args,
                    theta_low=theta_low,
                    theta_high=theta_high,
                )

                (
                    num[ithr],
                    num_unc[ithr],
                    _,
                ) = fit_numerator_for_score_threshold(
                    numerator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    float(threshold),
                    MGG_FIT_NOMINAL_ANGLE,
                    args,
                    theta_low=theta_low,
                    theta_high=theta_high,
                )

                eff[ithr] = np.divide(
                    num[ithr],
                    den[ithr],
                    out=np.full(n_energy, np.nan),
                    where=den[ithr] > 0.0,
                )
                eff_unc[ithr] = np.divide(
                    num_unc[ithr],
                    den[ithr],
                    out=np.full(n_energy, np.nan),
                    where=den[ithr] > 0.0,
                )
            #endfor

            results[source_name][label] = {
                "denominator": den,
                "numerator": num,
                "numerator_uncertainty": num_unc,
                "efficiency": eff,
                "efficiency_uncertainty": eff_unc,
                "theta_low": np.asarray([theta_low], dtype=float),
                "theta_high": np.asarray([theta_high], dtype=float),
            }
        #endfor
    #endfor

    return results


def plot_theta_binned_nominal_data_mc_ratio(
    results: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Nominal DATA/MC correction versus E_pred in one FT and three FD theta bins.
    """
    theta_bins = probe_theta_bins(args)

    # Use the scan point nearest the nominal requested BDT threshold.
    ithr = int(
        np.argmin(
            np.abs(
                BDT_EFFICIENCY_SCAN_THRESHOLDS
                - float(args.tag_score_min)
            )
        )
    )
    used_threshold = float(
        BDT_EFFICIENCY_SCAN_THRESHOLDS[ithr]
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(11.5, 8.5),
        squeeze=False,
    )

    centers = 0.5 * (edges[:-1] + edges[1:])

    for iax, (label, _, _, _) in enumerate(theta_bins):
        ax = axes.flat[iax]
        data = results["data"][label]

        for mc_name in ["aaogen", "clasdis"]:
            mc = results[mc_name][label]

            data_eff = data["efficiency"][ithr]
            data_unc = data["efficiency_uncertainty"][ithr]
            mc_eff = mc["efficiency"][ithr]
            mc_unc = mc["efficiency_uncertainty"][ithr]

            valid = (
                np.isfinite(data_eff)
                & np.isfinite(data_unc)
                & np.isfinite(mc_eff)
                & np.isfinite(mc_unc)
                & (data_eff > 0.0)
                & (mc_eff > 0.0)
            )

            ratio = np.full_like(data_eff, np.nan)
            ratio_unc = np.full_like(data_eff, np.nan)
            ratio[valid] = data_eff[valid] / mc_eff[valid]
            ratio_unc[valid] = ratio[valid] * np.sqrt(
                (data_unc[valid] / data_eff[valid]) ** 2
                + (mc_unc[valid] / mc_eff[valid]) ** 2
            )

            ax.errorbar(
                centers,
                ratio,
                yerr=ratio_unc,
                marker="o",
                linestyle="-",
                capsize=2,
                label=f"DATA / {mc_name.upper()}",
            )
        #endfor

        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.set_title(label)
        ax.set_xlabel(r"Predicted probe energy $E_{\rm pred}$ (GeV)")
        ax.set_ylabel("Photon-efficiency ratio DATA / MC")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8)
    #endfor

    fig.suptitle(
        "Probe-theta-binned photon-efficiency correction "
        f"(BDT >= {used_threshold:.2f}, 10 deg association)",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_fd_theta_binned_bdt_ratio_scan(
    results: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    DATA/MC ratio versus BDT threshold after controlling the probe theta.

    Rows are the three coarse FD theta bins; columns are E_pred bins.  If the
    broad-bin threshold dependence was largely caused by theta migration,
    these curves should be appreciably flatter than the old FD-integrated
    threshold scan.
    """
    fd_bins = [
        item
        for item in probe_theta_bins(args)
        if item[1] == REGIONS["FD"]
    ]
    n_energy = len(edges) - 1

    fig, axes = plt.subplots(
        len(fd_bins),
        n_energy,
        figsize=(5.0 * n_energy, 4.1 * len(fd_bins)),
        squeeze=False,
        sharex=True,
    )

    for irow, (label, _, _, _) in enumerate(fd_bins):
        data = results["data"][label]

        for ibin in range(n_energy):
            ax = axes[irow, ibin]

            for mc_name in ["aaogen", "clasdis"]:
                mc = results[mc_name][label]

                data_eff = data["efficiency"][:, ibin]
                data_unc = data["efficiency_uncertainty"][:, ibin]
                mc_eff = mc["efficiency"][:, ibin]
                mc_unc = mc["efficiency_uncertainty"][:, ibin]

                valid = (
                    np.isfinite(data_eff)
                    & np.isfinite(data_unc)
                    & np.isfinite(mc_eff)
                    & np.isfinite(mc_unc)
                    & (data_eff > 0.0)
                    & (mc_eff > 0.0)
                )

                ratio = np.full_like(data_eff, np.nan)
                ratio_unc = np.full_like(data_eff, np.nan)
                ratio[valid] = data_eff[valid] / mc_eff[valid]
                ratio_unc[valid] = ratio[valid] * np.sqrt(
                    (data_unc[valid] / data_eff[valid]) ** 2
                    + (mc_unc[valid] / mc_eff[valid]) ** 2
                )

                ax.errorbar(
                    BDT_EFFICIENCY_SCAN_THRESHOLDS,
                    ratio,
                    yerr=ratio_unc,
                    marker="o",
                    linestyle="-",
                    capsize=2,
                    label=f"DATA / {mc_name.upper()}",
                )
            #endfor

            ax.axhline(1.0, linestyle="--", linewidth=1.0)
            ax.axvline(
                float(args.tag_score_min),
                linestyle=":",
                linewidth=1.0,
            )
            ax.set_title(
                f"{label}; "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("Minimum tag BDT score")
            ax.set_ylabel("Efficiency ratio DATA / MC")
            ax.set_xlim(
                BDT_EFFICIENCY_SCAN_THRESHOLDS[0] - 0.01,
                0.99,
            )
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "FD DATA/MC efficiency ratio versus BDT score "
        "inside coarse probe-theta bins",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)



def build_fd_theta_sector_nominal_scan(
    denominator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    numerator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, Dict[str, Dict[int, Dict[str, np.ndarray]]]]:
    """
    Build the nominal FD correction in E_pred x theta_pred x sector.

    Only the nominal BDT threshold is used here.  The purpose is to establish
    whether sector-to-sector response differences are large enough to retain
    in the final photon-efficiency correction.
    """
    results: Dict[
        str,
        Dict[str, Dict[int, Dict[str, np.ndarray]]],
    ] = {}
    fd_bins = [
        item
        for item in probe_theta_bins(args)
        if item[1] == REGIONS["FD"]
    ]
    n_energy = len(edges) - 1

    for source_name in ["data", "aaogen", "clasdis"]:
        results[source_name] = {}

        for label, region_value, theta_low, theta_high in fd_bins:
            results[source_name][label] = {}

            for sector in range(1, 7):
                den = denominator_counts_for_score_threshold(
                    denominator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    float(args.tag_score_min),
                    args,
                    theta_low=theta_low,
                    theta_high=theta_high,
                    fd_sector=sector,
                )
                num, num_unc, _ = fit_numerator_for_score_threshold(
                    numerator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    float(args.tag_score_min),
                    MGG_FIT_NOMINAL_ANGLE,
                    args,
                    theta_low=theta_low,
                    theta_high=theta_high,
                    fd_sector=sector,
                )

                eff = np.divide(
                    num,
                    den,
                    out=np.full(n_energy, np.nan),
                    where=den > 0.0,
                )
                eff_unc = np.divide(
                    num_unc,
                    den,
                    out=np.full(n_energy, np.nan),
                    where=den > 0.0,
                )

                results[source_name][label][sector] = {
                    "denominator": den,
                    "numerator": num,
                    "numerator_uncertainty": num_unc,
                    "efficiency": eff,
                    "efficiency_uncertainty": eff_unc,
                }
            #endfor
        #endfor
    #endfor

    return results


def plot_fd_sector_dependence_nominal(
    results: Mapping[
        str,
        Mapping[str, Mapping[int, Mapping[str, np.ndarray]]],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Sector dependence of the nominal DATA/MC correction.

    Rows are coarse FD theta bins, columns are E_pred bins, and the horizontal
    axis is CLAS12 sector.  This compactly shows whether the final correction
    needs explicit sector dependence.
    """
    fd_bins = [
        item
        for item in probe_theta_bins(args)
        if item[1] == REGIONS["FD"]
    ]
    n_energy = len(edges) - 1
    sectors = np.arange(1, 7, dtype=float)

    fig, axes = plt.subplots(
        len(fd_bins),
        n_energy,
        figsize=(5.0 * n_energy, 4.1 * len(fd_bins)),
        squeeze=False,
        sharex=True,
    )

    for irow, (label, _, _, _) in enumerate(fd_bins):
        for ibin in range(n_energy):
            ax = axes[irow, ibin]

            for mc_name in ["aaogen", "clasdis"]:
                ratio = np.full(6, np.nan)
                ratio_unc = np.full(6, np.nan)

                for isec, sector in enumerate(range(1, 7)):
                    data = results["data"][label][sector]
                    mc = results[mc_name][label][sector]

                    data_eff = data["efficiency"][ibin]
                    data_unc = data["efficiency_uncertainty"][ibin]
                    mc_eff = mc["efficiency"][ibin]
                    mc_unc = mc["efficiency_uncertainty"][ibin]

                    if (
                        np.isfinite(data_eff)
                        and np.isfinite(data_unc)
                        and np.isfinite(mc_eff)
                        and np.isfinite(mc_unc)
                        and data_eff > 0.0
                        and mc_eff > 0.0
                    ):
                        ratio[isec] = data_eff / mc_eff
                        ratio_unc[isec] = ratio[isec] * np.sqrt(
                            (data_unc / data_eff) ** 2
                            + (mc_unc / mc_eff) ** 2
                        )
                    #endif
                #endfor

                ax.errorbar(
                    sectors,
                    ratio,
                    yerr=ratio_unc,
                    marker="o",
                    linestyle="-",
                    capsize=2,
                    label=f"DATA / {mc_name.upper()}",
                )
            #endfor

            ax.axhline(1.0, linestyle="--", linewidth=1.0)
            ax.set_title(
                f"{label}; "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("FD sector")
            ax.set_ylabel("Photon-efficiency ratio DATA / MC")
            ax.set_xticks(range(1, 7))
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "Nominal FD photon-efficiency correction versus sector "
        f"(BDT >= {args.tag_score_min:.2f}, 10 deg association)",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)


def plot_fd_sector_averaged_correction(
    results: Mapping[
        str,
        Mapping[str, Mapping[int, Mapping[str, np.ndarray]]],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Show all six sector-specific correction curves versus E_pred, one panel
    per coarse FD theta bin.  DATA/AAOgen is used for the main sector-shape
    view; DATA/CLASDIS remains available in the companion sector canvas and
    printed numerical summary.
    """
    fd_bins = [
        item
        for item in probe_theta_bins(args)
        if item[1] == REGIONS["FD"]
    ]
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, axes = plt.subplots(
        1,
        len(fd_bins),
        figsize=(5.2 * len(fd_bins), 4.8),
        squeeze=False,
    )

    for iax, (label, _, _, _) in enumerate(fd_bins):
        ax = axes[0, iax]

        for sector in range(1, 7):
            data = results["data"][label][sector]
            mc = results["aaogen"][label][sector]

            data_eff = data["efficiency"]
            data_unc = data["efficiency_uncertainty"]
            mc_eff = mc["efficiency"]
            mc_unc = mc["efficiency_uncertainty"]

            valid = (
                np.isfinite(data_eff)
                & np.isfinite(data_unc)
                & np.isfinite(mc_eff)
                & np.isfinite(mc_unc)
                & (data_eff > 0.0)
                & (mc_eff > 0.0)
            )
            ratio = np.full_like(data_eff, np.nan)
            ratio_unc = np.full_like(data_eff, np.nan)
            ratio[valid] = data_eff[valid] / mc_eff[valid]
            ratio_unc[valid] = ratio[valid] * np.sqrt(
                (data_unc[valid] / data_eff[valid]) ** 2
                + (mc_unc[valid] / mc_eff[valid]) ** 2
            )

            ax.errorbar(
                centers,
                ratio,
                yerr=ratio_unc,
                marker="o",
                linestyle="-",
                capsize=2,
                label=f"S{sector}",
            )
        #endfor

        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.set_title(label)
        ax.set_xlabel(r"Predicted probe energy $E_{\rm pred}$ (GeV)")
        ax.set_ylabel("DATA / AAOgen efficiency ratio")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, ncol=2)
    #endfor

    fig.suptitle(
        "FD sector-specific photon-efficiency corrections",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    save_figure(fig, output_path)



def plot_predicted_vs_reconstructed_sector_migration(
    numerator_diagnostics_by_source: Mapping[
        str,
        Mapping[str, np.ndarray],
    ],
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Compare the sector assigned from the missing-vector predicted probe phi
    with the sector of the actually reconstructed partner photon.

    Each predicted-sector row is normalized to unity.  A strong diagonal is
    required before large sector-dependent correction factors are interpreted
    as detector effects rather than sector-migration artifacts.
    """
    sources = ["data", "aaogen", "clasdis"]
    fig, axes = plt.subplots(
        1,
        len(sources),
        figsize=(5.2 * len(sources), 4.8),
        squeeze=False,
    )

    im = None

    for iax, source_name in enumerate(sources):
        ax = axes[0, iax]
        d = numerator_diagnostics_by_source[source_name]

        score = np.asarray(d["scores"], dtype=float)
        pred_theta = np.asarray(d["pred_theta"], dtype=float)
        pred_phi = np.asarray(d["pred_phi"], dtype=float)
        partner_phi = np.asarray(d["partner_phi"], dtype=float)
        partner_energy = np.asarray(d["partner_energy"], dtype=float)
        match_angle = np.asarray(d["match_angle"], dtype=float)
        mx2_ep = np.asarray(d["mx2_ep"], dtype=float)
        tag_energy = np.asarray(d["tag_energy"], dtype=float)

        pred_sector = assign_fd_sector(pred_phi)
        reco_sector = assign_fd_sector(partner_phi)

        mask = (
            np.isfinite(score)
            & (score >= args.tag_score_min)
            & np.isfinite(tag_energy)
            & (tag_energy >= args.tag_energy_min)
            & np.isfinite(mx2_ep)
            & (mx2_ep >= args.mx2_ep_min)
            & (mx2_ep < args.mx2_ep_max)
            & np.isfinite(pred_theta)
            & (pred_theta >= args.fd_theta_min)
            & (pred_theta < args.fd_theta_max)
            & np.isfinite(partner_energy)
            & (partner_energy >= args.probe_energy_min)
            & np.isfinite(match_angle)
            & (match_angle < MGG_FIT_NOMINAL_ANGLE)
            & (pred_sector >= 1)
            & (reco_sector >= 1)
        )

        matrix = np.zeros((6, 6), dtype=float)

        for pred_sec in range(1, 7):
            for reco_sec in range(1, 7):
                matrix[pred_sec - 1, reco_sec - 1] = np.sum(
                    mask
                    & (pred_sector == pred_sec)
                    & (reco_sector == reco_sec)
                )
            #endfor
        #endfor

        row_sum = np.sum(matrix, axis=1, keepdims=True)
        matrix_norm = np.divide(
            matrix,
            row_sum,
            out=np.zeros_like(matrix),
            where=row_sum > 0.0,
        )

        im = ax.imshow(
            matrix_norm,
            origin="lower",
            vmin=0.0,
            vmax=1.0,
            aspect="equal",
        )

        for irow in range(6):
            for icol in range(6):
                ax.text(
                    icol,
                    irow,
                    f"{100.0 * matrix_norm[irow, icol]:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8,
                )
            #endfor
        #endfor

        diagonal = np.trace(matrix)
        total = np.sum(matrix)
        agreement = diagonal / total if total > 0.0 else np.nan

        ax.set_title(
            f"{source_name.upper()}\n"
            f"same-sector fraction = {agreement:.3f}"
        )
        ax.set_xlabel("Reconstructed partner sector")
        ax.set_ylabel("Predicted probe sector")
        ax.set_xticks(np.arange(6))
        ax.set_xticklabels(np.arange(1, 7))
        ax.set_yticks(np.arange(6))
        ax.set_yticklabels(np.arange(1, 7))
    #endfor

    if im is not None:
        fig.colorbar(
            im,
            ax=axes.ravel().tolist(),
            fraction=0.025,
            pad=0.02,
            label="Row-normalized fraction",
        )
    #endif

    fig.suptitle(
        "FD predicted-versus-reconstructed probe sector migration "
        "(nominal BDT and 10 deg association)",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 0.97, 0.94])
    save_figure(fig, output_path)


def plot_fd_sector_statistics(
    results: Mapping[
        str,
        Mapping[str, Mapping[int, Mapping[str, np.ndarray]]],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Show the nominal tag-denominator population in every sector cell.

    This is primarily a fit-stability diagnostic: sector corrections with very
    small populations should not be interpreted as true detector response.
    """
    fd_bins = [
        item
        for item in probe_theta_bins(args)
        if item[1] == REGIONS["FD"]
    ]
    n_energy = len(edges) - 1
    sectors = np.arange(1, 7)

    fig, axes = plt.subplots(
        len(fd_bins),
        n_energy,
        figsize=(4.4 * n_energy, 4.0 * len(fd_bins)),
        squeeze=False,
        sharex=True,
    )

    for irow, (label, _, _, _) in enumerate(fd_bins):
        for ibin in range(n_energy):
            ax = axes[irow, ibin]

            for source_name in ["data", "aaogen", "clasdis"]:
                counts = np.asarray(
                    [
                        results[source_name][label][sector][
                            "denominator"
                        ][ibin]
                        for sector in range(1, 7)
                    ],
                    dtype=float,
                )
                ax.plot(
                    sectors,
                    counts,
                    marker="o",
                    linestyle="-",
                    label=source_name.upper(),
                )
            #endfor

            ax.set_yscale("log")
            ax.set_title(
                f"{label}; "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("FD sector")
            ax.set_ylabel("Tag denominator count")
            ax.set_xticks(range(1, 7))
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "Nominal FD sector statistics",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)


def plot_fd_sector_raw_efficiencies(
    results: Mapping[
        str,
        Mapping[str, Mapping[int, Mapping[str, np.ndarray]]],
    ],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Plot fitted apparent efficiency before DATA/MC division.

    This identifies whether an anomalous correction originates predominantly
    in DATA, AAOgen, CLASDIS, or a low-statistics Mgg fit.
    """
    fd_bins = [
        item
        for item in probe_theta_bins(args)
        if item[1] == REGIONS["FD"]
    ]
    n_energy = len(edges) - 1
    sectors = np.arange(1, 7)

    fig, axes = plt.subplots(
        len(fd_bins),
        n_energy,
        figsize=(4.4 * n_energy, 4.0 * len(fd_bins)),
        squeeze=False,
        sharex=True,
    )

    for irow, (label, _, _, _) in enumerate(fd_bins):
        for ibin in range(n_energy):
            ax = axes[irow, ibin]

            for source_name in ["data", "aaogen", "clasdis"]:
                eff = np.asarray(
                    [
                        results[source_name][label][sector][
                            "efficiency"
                        ][ibin]
                        for sector in range(1, 7)
                    ],
                    dtype=float,
                )
                unc = np.asarray(
                    [
                        results[source_name][label][sector][
                            "efficiency_uncertainty"
                        ][ibin]
                        for sector in range(1, 7)
                    ],
                    dtype=float,
                )

                ax.errorbar(
                    sectors,
                    eff,
                    yerr=unc,
                    marker="o",
                    linestyle="-",
                    capsize=2,
                    label=source_name.upper(),
                )
            #endfor

            ax.set_title(
                f"{label}; "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("FD sector")
            ax.set_ylabel(r"Fitted $N_{\pi^0}/N_{\rm tag}$")
            ax.set_xticks(range(1, 7))
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "Nominal apparent photon efficiency versus FD sector",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)


def plot_aaogen_truth_tag_score_phase_space(
    diagnostics: Mapping[str, np.ndarray],
    edges: np.ndarray,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Diagnose why truth-matched AAOgen efficiency changes with BDT threshold.

    Focus on FD, where the final correction is currently best constrained.
    Within each broad E_pred bin, overlay unit-normalized distributions of
    predicted probe theta and tag energy for DISJOINT tag-score ranges:
      0.50-0.60, 0.60-0.70, 0.70-0.80, and >=0.80.

    A corresponding numerical table of means for E_pred, theta_pred, E_tag,
    Mx2(ep), and Mx2(epgamma_tag) is printed to stdout.
    """
    score = np.asarray(diagnostics["all_score"], dtype=float)
    tag_energy = np.asarray(
        diagnostics["all_tag_energy"],
        dtype=float,
    )
    mx2_ep = np.asarray(
        diagnostics["all_mx2_ep"],
        dtype=float,
    )
    mx2_epgamma = np.asarray(
        diagnostics["all_missing_mass2"],
        dtype=float,
    )
    pred_energy = np.asarray(
        diagnostics["all_pred_energy"],
        dtype=float,
    )
    pred_theta = np.asarray(
        diagnostics["all_pred_theta"],
        dtype=float,
    )
    truth_pid = np.asarray(
        diagnostics["all_matching_gamma_pid"],
        dtype=np.int64,
    )
    truth_mcindex = np.asarray(
        diagnostics["all_gamma_mcindex"],
        dtype=np.int64,
    )
    region = assign_region(pred_theta, args)

    base = (
        np.isfinite(score)
        & np.isfinite(tag_energy)
        & (tag_energy >= args.tag_energy_min)
        & np.isfinite(mx2_ep)
        & (mx2_ep >= args.mx2_ep_min)
        & (mx2_ep < args.mx2_ep_max)
        & np.isfinite(mx2_epgamma)
        & np.isfinite(pred_energy)
        & (pred_energy >= args.probe_energy_min)
        & (pred_energy < args.probe_energy_max)
        & np.isfinite(pred_theta)
        & (region == REGIONS["FD"])
        & (truth_pid == 22)
        & (truth_mcindex >= 0)
    )

    score_ranges = [
        (0.50, 0.60, "0.50-0.60"),
        (0.60, 0.70, "0.60-0.70"),
        (0.70, 0.80, "0.70-0.80"),
        (0.80, np.inf, ">=0.80"),
    ]

    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(5.1 * nbin, 8.0),
        squeeze=False,
    )

    progress(
        "AAOgen FD truth-matched tag phase-space means by disjoint BDT range:"
    )

    for ibin in range(nbin):
        e_mask = (
            base
            & (pred_energy >= edges[ibin])
            & (pred_energy < edges[ibin + 1])
        )

        # Use common axes within each E_pred bin so score-range shape changes
        # are visually meaningful.
        theta_values = pred_theta[e_mask]
        tag_values = tag_energy[e_mask]

        if np.any(np.isfinite(theta_values)):
            theta_lo, theta_hi = np.nanpercentile(
                theta_values,
                [0.5, 99.5],
            )
        else:
            theta_lo, theta_hi = 5.5, 35.0
        #endif

        if np.any(np.isfinite(tag_values)):
            tag_lo, tag_hi = np.nanpercentile(
                tag_values,
                [0.5, 99.5],
            )
        else:
            tag_lo, tag_hi = args.tag_energy_min, 9.5
        #endif

        theta_bins = np.linspace(theta_lo, theta_hi, 45)
        tag_bins = np.linspace(tag_lo, tag_hi, 45)

        print(
            f"  FD Epred {edges[ibin]:g}-{edges[ibin+1]:g} GeV:"
        )

        for low, high, label in score_ranges:
            score_mask = (
                e_mask
                & (score >= low)
                & (score < high)
            )

            n = int(np.sum(score_mask))
            if n == 0:
                print(f"    score {label:9s}: N=0")
                continue
            #endif

            axes[0, ibin].hist(
                pred_theta[score_mask],
                bins=theta_bins,
                histtype="step",
                density=True,
                linewidth=1.5,
                label=label,
            )
            axes[1, ibin].hist(
                tag_energy[score_mask],
                bins=tag_bins,
                histtype="step",
                density=True,
                linewidth=1.5,
                label=label,
            )

            def mean_or_nan(values: np.ndarray) -> float:
                finite = np.isfinite(values)
                return (
                    float(np.mean(values[finite]))
                    if np.any(finite)
                    else np.nan
                )
            #enddef

            print(
                f"    score {label:9s}: N={n:9,d}, "
                f"<Epred>={mean_or_nan(pred_energy[score_mask]):.3f} GeV, "
                f"<theta_pred>={mean_or_nan(pred_theta[score_mask]):.3f} deg, "
                f"<Etag>={mean_or_nan(tag_energy[score_mask]):.3f} GeV, "
                f"<Mx2(ep)>={mean_or_nan(mx2_ep[score_mask]):.4f} GeV^2, "
                f"<Mx2(epgamma)>={mean_or_nan(mx2_epgamma[score_mask]):.4f} GeV^2"
            )
        #endfor

        axes[0, ibin].set_title(
            f"FD, {edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
        )
        axes[0, ibin].set_xlabel(
            r"Predicted probe $\theta$ (deg)"
        )
        axes[0, ibin].set_ylabel("Unit-normalized density")
        axes[0, ibin].grid(alpha=0.2)

        axes[1, ibin].set_xlabel(r"Tag photon energy $E_{\rm tag}$ (GeV)")
        axes[1, ibin].set_ylabel("Unit-normalized density")
        axes[1, ibin].grid(alpha=0.2)

        if ibin == 0:
            axes[0, ibin].legend(
                title="BDT score",
                fontsize=8,
                title_fontsize=8,
            )
        #endif
    #endfor

    fig.suptitle(
        "AAOgen truth-matched FD tags: phase-space evolution with BDT score",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_aaogen_truth_tag_threshold_scan(
    standard_scan: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    truth_scan: Mapping[str, Mapping[str, np.ndarray]],
    edges: np.ndarray,
    output_path: Path,
) -> None:
    """
    Compare ordinary AAOgen and truth-matched-tag AAOgen efficiency versus
    BDT threshold.  If both retain similar threshold dependence, the effect
    is selection/phase-space correlation rather than AAOgen tag impurity.
    """
    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(5.0 * nbin, 8.5),
        squeeze=False,
        sharex=True,
    )

    for irow, region_name in enumerate(["FT", "FD"]):
        standard = standard_scan["aaogen"][region_name]
        truth = truth_scan[region_name]

        for ibin in range(nbin):
            ax = axes[irow, ibin]
            ax.errorbar(
                BDT_EFFICIENCY_SCAN_THRESHOLDS,
                standard["efficiency"][:, ibin],
                yerr=standard["efficiency_uncertainty"][:, ibin],
                marker="o",
                linestyle="-",
                capsize=2,
                label="AAOgen all selected tags",
            )
            ax.errorbar(
                BDT_EFFICIENCY_SCAN_THRESHOLDS,
                truth["efficiency"][:, ibin],
                yerr=truth["efficiency_uncertainty"][:, ibin],
                marker="s",
                linestyle="--",
                capsize=2,
                label="AAOgen truth-matched tags",
            )
            ax.axvline(0.80, linestyle=":", linewidth=1.0)
            ax.set_title(
                f"{region_name}, "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("Minimum tag BDT score")
            ax.set_ylabel(r"Fitted $N_{\pi^0}/N_{\rm tag}$")
            ax.grid(alpha=0.2)
            ax.set_xlim(
                BDT_EFFICIENCY_SCAN_THRESHOLDS[0] - 0.01,
                0.99,
            )

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "AAOgen BDT-threshold dependence: all tags versus truth-matched tags",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_data_truth_aaogen_ratio_scan(
    standard_scan: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    truth_scan: Mapping[str, Mapping[str, np.ndarray]],
    edges: np.ndarray,
    output_path: Path,
) -> None:
    """
    DATA / truth-matched-AAOgen apparent-efficiency ratio versus BDT threshold.

    Persistence of threshold dependence here demonstrates that the effect
    cannot be explained by unmatched/nonphoton AAOgen denominator tags.
    """
    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(5.0 * nbin, 8.5),
        squeeze=False,
        sharex=True,
    )

    for irow, region_name in enumerate(["FT", "FD"]):
        data = standard_scan["data"][region_name]
        truth = truth_scan[region_name]

        for ibin in range(nbin):
            ax = axes[irow, ibin]
            data_eff = data["efficiency"][:, ibin]
            data_unc = data["efficiency_uncertainty"][:, ibin]
            mc_eff = truth["efficiency"][:, ibin]
            mc_unc = truth["efficiency_uncertainty"][:, ibin]

            valid = (
                np.isfinite(data_eff)
                & np.isfinite(data_unc)
                & np.isfinite(mc_eff)
                & np.isfinite(mc_unc)
                & (data_eff > 0.0)
                & (mc_eff > 0.0)
            )

            ratio = np.full_like(data_eff, np.nan, dtype=float)
            ratio_unc = np.full_like(data_eff, np.nan, dtype=float)
            ratio[valid] = data_eff[valid] / mc_eff[valid]
            ratio_unc[valid] = ratio[valid] * np.sqrt(
                (data_unc[valid] / data_eff[valid]) ** 2
                + (mc_unc[valid] / mc_eff[valid]) ** 2
            )

            ax.errorbar(
                BDT_EFFICIENCY_SCAN_THRESHOLDS,
                ratio,
                yerr=ratio_unc,
                marker="o",
                linestyle="-",
                capsize=2,
                label="DATA / truth-matched AAOgen",
            )
            ax.axhline(1.0, linestyle="--", linewidth=1.0)
            ax.axvline(0.80, linestyle=":", linewidth=1.0)
            ax.set_title(
                f"{region_name}, "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("Minimum tag BDT score")
            ax.set_ylabel("Apparent efficiency ratio DATA / MC")
            ax.grid(alpha=0.2)
            ax.set_xlim(
                BDT_EFFICIENCY_SCAN_THRESHOLDS[0] - 0.01,
                0.99,
            )

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "DATA / truth-matched AAOgen efficiency versus minimum tag BDT score",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_bdt_efficiency_threshold_scan(
    results: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    edges: np.ndarray,
    output_path: Path,
) -> None:
    """
    Primary diagnostic: apparent efficiency versus minimum BDT tag score.

    Rows are FT/FD; columns are E_pred bins.  Data, AAOgen, and CLASDIS are
    shown together so a purity-driven rise in data can be compared directly
    with the MC behavior.
    """
    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(5.0 * nbin, 8.5),
        squeeze=False,
        sharex=True,
    )

    for irow, region_name in enumerate(["FT", "FD"]):
        for ibin in range(nbin):
            ax = axes[irow, ibin]

            for source_name in ["data", "aaogen", "clasdis"]:
                block = results[source_name][region_name]
                ax.errorbar(
                    BDT_EFFICIENCY_SCAN_THRESHOLDS,
                    block["efficiency"][:, ibin],
                    yerr=block["efficiency_uncertainty"][:, ibin],
                    marker="o",
                    linestyle="-",
                    capsize=2,
                    label=source_name.upper(),
                )
            #endfor

            ax.axvline(
                0.80,
                linestyle="--",
                linewidth=1.0,
            )
            ax.set_title(
                f"{region_name}, "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("Minimum tag BDT score")
            ax.set_ylabel(
                r"Fitted $N_{\pi^0}/N_{\rm tag}$"
            )
            ax.set_xlim(
                BDT_EFFICIENCY_SCAN_THRESHOLDS[0] - 0.01,
                0.99,
            )
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "Apparent photon efficiency versus minimum tag BDT score "
        "(10 deg association)",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)



def plot_bdt_threshold_data_mc_ratio(
    results: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    edges: np.ndarray,
    output_path: Path,
) -> None:
    """
    Plot the data/MC apparent-efficiency ratio versus minimum BDT tag score.

    This is the key working-point diagnostic.  Absolute apparent efficiencies
    may change with tag-score threshold because the selected pi0 population
    changes, especially for CLASDIS where genuine non-pi0 photons are present.
    What matters for the correction is whether the DATA/MC ratio is stable.

    Ratios are shown separately for AAOgen and CLASDIS in each FT/FD and
    E_pred bin.  Statistical uncertainties from the fitted numerator yields
    are propagated in quadrature.
    """
    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(5.0 * nbin, 8.5),
        squeeze=False,
        sharex=True,
    )

    for irow, region_name in enumerate(["FT", "FD"]):
        data_block = results["data"][region_name]

        for ibin in range(nbin):
            ax = axes[irow, ibin]

            data_eff = data_block["efficiency"][:, ibin]
            data_unc = data_block["efficiency_uncertainty"][:, ibin]

            for mc_name in ["aaogen", "clasdis"]:
                mc_block = results[mc_name][region_name]
                mc_eff = mc_block["efficiency"][:, ibin]
                mc_unc = mc_block["efficiency_uncertainty"][:, ibin]

                valid = (
                    np.isfinite(data_eff)
                    & np.isfinite(data_unc)
                    & np.isfinite(mc_eff)
                    & np.isfinite(mc_unc)
                    & (data_eff > 0.0)
                    & (mc_eff > 0.0)
                )

                ratio = np.full_like(data_eff, np.nan, dtype=float)
                ratio_unc = np.full_like(data_eff, np.nan, dtype=float)

                ratio[valid] = (
                    data_eff[valid] / mc_eff[valid]
                )

                rel_data = np.zeros_like(data_eff, dtype=float)
                rel_mc = np.zeros_like(mc_eff, dtype=float)
                rel_data[valid] = (
                    data_unc[valid] / data_eff[valid]
                )
                rel_mc[valid] = (
                    mc_unc[valid] / mc_eff[valid]
                )
                ratio_unc[valid] = ratio[valid] * np.sqrt(
                    rel_data[valid] * rel_data[valid]
                    + rel_mc[valid] * rel_mc[valid]
                )

                ax.errorbar(
                    BDT_EFFICIENCY_SCAN_THRESHOLDS,
                    ratio,
                    yerr=ratio_unc,
                    marker="o",
                    linestyle="-",
                    capsize=2,
                    label=f"DATA / {mc_name.upper()}",
                )
            #endfor

            ax.axhline(
                1.0,
                linestyle="--",
                linewidth=1.0,
            )
            ax.axvline(
                0.80,
                linestyle=":",
                linewidth=1.0,
            )
            ax.set_title(
                f"{region_name}, "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("Minimum tag BDT score")
            ax.set_ylabel("Apparent efficiency ratio DATA / MC")
            ax.set_xlim(
                BDT_EFFICIENCY_SCAN_THRESHOLDS[0] - 0.01,
                0.99,
            )
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "DATA/MC photon-efficiency ratio versus minimum tag BDT score "
        "(10 deg association)",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_data_bdt_threshold_retention(
    results: Mapping[
        str,
        Mapping[str, Mapping[str, np.ndarray]],
    ],
    edges: np.ndarray,
    output_path: Path,
) -> None:
    """
    Data-only retention diagnostic.

    Denominator and fitted pi0 numerator are independently normalized to their
    values at the loosest scanned threshold.  If denominator contamination is
    important, tightening the score should suppress the denominator faster
    than the fitted pi0 numerator.
    """
    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(5.0 * nbin, 8.5),
        squeeze=False,
        sharex=True,
        sharey=True,
    )

    for irow, region_name in enumerate(["FT", "FD"]):
        block = results["data"][region_name]

        for ibin in range(nbin):
            ax = axes[irow, ibin]
            den = block["denominator"][:, ibin]
            num = block["numerator"][:, ibin]

            den0 = den[0]
            num0 = num[0]

            den_ret = np.divide(
                den,
                den0,
                out=np.full_like(den, np.nan),
                where=den0 > 0.0,
            )
            num_ret = np.divide(
                num,
                num0,
                out=np.full_like(num, np.nan),
                where=np.isfinite(num0) & (num0 > 0.0),
            )

            ax.plot(
                BDT_EFFICIENCY_SCAN_THRESHOLDS,
                den_ret,
                marker="o",
                linestyle="-",
                label="Tag denominator",
            )
            ax.plot(
                BDT_EFFICIENCY_SCAN_THRESHOLDS,
                num_ret,
                marker="s",
                linestyle="--",
                label=r"Fitted $\pi^0$ numerator",
            )
            ax.axvline(
                0.80,
                linestyle=":",
                linewidth=1.0,
            )
            ax.set_title(
                f"{region_name}, "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel("Minimum tag BDT score")
            ax.set_ylabel("Retention relative to score >= 0.50")
            ax.set_xlim(
                BDT_EFFICIENCY_SCAN_THRESHOLDS[0] - 0.01,
                0.99,
            )
            ax.set_ylim(0.0, 1.05)
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        "Data: denominator versus fitted-pi0 numerator retention "
        "under tighter BDT selection",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_mgg_fit_examples(
    source_name: str,
    fits_by_region: Mapping[str, list[MassFitResult]],
    edges: np.ndarray,
    output_path: Path,
    angle_max: float,
) -> None:
    nbin = len(edges) - 1
    fig, axes = plt.subplots(
        2,
        nbin,
        figsize=(4.8 * nbin, 8.2),
        squeeze=False,
    )

    for irow, region_name in enumerate(["FT", "FD"]):
        fits = fits_by_region[region_name]

        for ibin in range(nbin):
            ax = axes[irow, ibin]
            fit = fits[ibin]

            ax.step(
                fit.x_centers,
                fit.counts,
                where="mid",
                linewidth=1.2,
                label="Candidates",
            )

            if fit.success:
                ax.plot(
                    fit.x_centers,
                    fit.model_counts,
                    linewidth=1.5,
                    label="Total fit",
                )
                ax.plot(
                    fit.x_centers,
                    fit.background_counts,
                    linestyle="--",
                    linewidth=1.2,
                    label="Background",
                )
                ax.plot(
                    fit.x_centers,
                    fit.signal_counts,
                    linestyle=":",
                    linewidth=1.2,
                    label=r"$\pi^0$ double Gaussian",
                )
                fit_text = (
                    f"Npi0={fit.yield_pi0:,.0f}\\n"
                    f"mu={1000.0*fit.mean:.1f} MeV\\n"
                    f"sigma_core={1000.0*fit.sigma_core:.1f} MeV\\n"
                    f"sigma_tail={1000.0*fit.sigma_tail:.1f} MeV\\n"
                    f"f_core={fit.core_fraction:.2f}\\n"
                    f"chi2/ndf={fit.chi2_ndf:.2f}"
                )
            else:
                fit_text = "fit failed / insufficient statistics"
            #endif

            ax.text(
                0.03,
                0.96,
                fit_text,
                transform=ax.transAxes,
                va="top",
                fontsize=8,
            )
            ax.set_title(
                f"{region_name}, "
                f"{edges[ibin]:g} <= Epred < {edges[ibin+1]:g} GeV"
            )
            ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
            ax.set_ylabel("Candidate pairs / bin")
            ax.grid(alpha=0.2)

            if irow == 0 and ibin == 0:
                ax.legend(fontsize=8)
            #endif
        #endfor
    #endfor

    fig.suptitle(
        f"{source_name}: M(gamma gamma) fits, "
        f"association < {angle_max:g} deg",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_fit_efficiency_angle_stability(
    source_name: str,
    curves_by_angle: Mapping[float, Mapping[str, FitEfficiencyCurve]],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9))

    for iax, region_name in enumerate(["FT", "FD"]):
        ax = axes[iax]

        for angle_max in MGG_FIT_ASSOCIATION_ANGLES:
            curve = curves_by_angle[angle_max][region_name]
            ax.errorbar(
                curve.centers,
                curve.efficiency,
                yerr=curve.uncertainty,
                marker="o",
                linestyle="-",
                capsize=2,
                label=f"{angle_max:g} deg",
            )
        #endfor

        ax.set_xlabel(r"Predicted probe energy $E_{\rm pred}$ (GeV)")
        ax.set_ylabel(
            r"Fitted $N_{\pi^0}$ / tag denominator"
        )
        ax.set_title(region_name)
        ax.grid(alpha=0.2)
        ax.legend()
    #endfor

    fig.suptitle(
        f"{source_name}: fitted-efficiency association-angle stability",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.94])
    save_figure(fig, output_path)


def plot_fit_efficiency_data_mc_ratio(
    fit_curves: Mapping[
        str,
        Mapping[str, FitEfficiencyCurve],
    ],
    output_path: Path,
) -> None:
    """
    Nominal 10-degree fitted efficiencies and data/MC ratios.

    With partial data or CLASDIS epgg files the absolute ratios are not yet
    physically interpretable; the plot is nevertheless useful once exposures
    are complete/matched.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))

    for icol, region_name in enumerate(["FT", "FD"]):
        ax_eff = axes[0, icol]
        ax_ratio = axes[1, icol]

        for source_name in ["data", "aaogen", "clasdis"]:
            curve = fit_curves[source_name][region_name]
            ax_eff.errorbar(
                curve.centers,
                curve.efficiency,
                yerr=curve.uncertainty,
                marker="o",
                linestyle="-",
                capsize=2,
                label=source_name.upper(),
            )
        #endfor

        data_curve = fit_curves["data"][region_name]

        for mc_name in ["aaogen", "clasdis"]:
            mc_curve = fit_curves[mc_name][region_name]
            ratio = np.divide(
                data_curve.efficiency,
                mc_curve.efficiency,
                out=np.full_like(data_curve.efficiency, np.nan),
                where=(
                    np.isfinite(mc_curve.efficiency)
                    & (mc_curve.efficiency > 0.0)
                ),
            )

            rel_data = np.divide(
                data_curve.uncertainty,
                data_curve.efficiency,
                out=np.zeros_like(data_curve.efficiency),
                where=(
                    np.isfinite(data_curve.efficiency)
                    & (data_curve.efficiency > 0.0)
                ),
            )
            rel_mc = np.divide(
                mc_curve.uncertainty,
                mc_curve.efficiency,
                out=np.zeros_like(mc_curve.efficiency),
                where=(
                    np.isfinite(mc_curve.efficiency)
                    & (mc_curve.efficiency > 0.0)
                ),
            )
            ratio_unc = ratio * np.sqrt(
                rel_data * rel_data
                + rel_mc * rel_mc
            )

            ax_ratio.errorbar(
                data_curve.centers,
                ratio,
                yerr=ratio_unc,
                marker="o",
                linestyle="-",
                capsize=2,
                label=f"data / {mc_name}",
            )
        #endfor

        ax_eff.set_title(f"{region_name}: fitted efficiency")
        ax_eff.set_xlabel(r"Predicted probe energy $E_{\rm pred}$ (GeV)")
        ax_eff.set_ylabel(r"Fitted $N_{\pi^0}/N_{\rm tag}$")
        ax_eff.grid(alpha=0.2)
        ax_eff.legend()

        ax_ratio.axhline(1.0, linestyle="--", linewidth=1.0)
        ax_ratio.set_title(f"{region_name}: data / MC")
        ax_ratio.set_xlabel(r"Predicted probe energy $E_{\rm pred}$ (GeV)")
        ax_ratio.set_ylabel("Efficiency ratio")
        ax_ratio.grid(alpha=0.2)
        ax_ratio.legend()
    #endfor

    fig.suptitle(
        "Nominal 10 deg M(gamma gamma)-fit efficiency",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_probe_energy_response_fine(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Fine response of E_pred against the reconstructed truth-qualified sister.

    NOTE: generated sister-photon kinematics are not present in the current
    ROOT schema.  The x-axis is therefore reconstructed sister energy after
    truth qualification, not E_true.  This is intentionally labeled as such
    rather than pretending that a generator-level response is available.
    """
    base = denominator_closure_base_mask(
        source_name.lower(),
        diagnostics,
        args,
    )

    if base is None or np.sum(base) == 0:
        return
    #endif

    pred_e = np.asarray(
        diagnostics["pred_energy"],
        dtype=float,
    )[base]
    sister_e = np.asarray(
        diagnostics["partner_energy"],
        dtype=float,
    )[base]

    edges = RESPONSE_ENERGY_EDGES
    matrix = np.zeros(
        (len(edges) - 1, len(edges) - 1),
        dtype=float,
    )

    sister_bin = np.digitize(sister_e, edges) - 1
    pred_bin = np.digitize(pred_e, edges) - 1

    for true_like_bin in range(len(edges) - 1):
        denom = np.sum(sister_bin == true_like_bin)

        for predicted_bin in range(len(edges) - 1):
            if denom > 0:
                matrix[predicted_bin, true_like_bin] = (
                    np.sum(
                        (pred_bin == predicted_bin)
                        & (sister_bin == true_like_bin)
                    )
                    / denom
                )
            #endif
        #endfor
    #endfor

    labels = [
        f"{edges[i]:g}-{edges[i+1]:g}"
        for i in range(len(edges) - 1)
    ]

    fig, ax = plt.subplots(figsize=(11.5, 9.0))
    im = ax.imshow(
        matrix,
        origin="lower",
        vmin=0.0,
        vmax=1.0,
        aspect="auto",
    )
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlabel(
        "Reconstructed truth-qualified sister E bin (GeV)"
    )
    ax.set_ylabel("Predicted probe E bin (GeV)")
    ax.set_title(
        f"{source_name}: fine probe-energy response; columns normalized"
    )
    fig.colorbar(im, ax=ax, label="Fraction")

    fig.tight_layout()
    save_figure(fig, output_path)



# =============================================================================
# Efficiency calculation
# =============================================================================

def efficiency_curve(
    denominator: DirectedSample,
    numerator: DirectedSample,
    region_value: int,
    edges: np.ndarray,
) -> EfficiencyCurve:
    den_values = denominator.predicted_energy[
        denominator.region == region_value
    ]
    num_values = numerator.predicted_energy[
        numerator.region == region_value
    ]

    den, _ = np.histogram(den_values, bins=edges)
    num, _ = np.histogram(num_values, bins=edges)

    den = den.astype(float)
    num = num.astype(float)

    efficiency = np.divide(
        num,
        den,
        out=np.full_like(den, np.nan),
        where=den > 0,
    )

    # First-stage statistical estimate.  When the topology bookkeeping is
    # perfectly subset-like, this is the familiar binomial uncertainty.
    variance = np.divide(
        efficiency * (1.0 - efficiency),
        den,
        out=np.full_like(den, np.nan),
        where=(den > 0) & (efficiency >= 0.0) & (efficiency <= 1.0),
    )
    uncertainty = np.sqrt(np.clip(variance, 0.0, None))

    centers = 0.5 * (edges[:-1] + edges[1:])

    return EfficiencyCurve(
        edges=edges,
        centers=centers,
        denominator=den,
        numerator=num,
        efficiency=efficiency,
        uncertainty=uncertainty,
    )


def ratio_curve(
    data: EfficiencyCurve,
    mc: EfficiencyCurve,
) -> Tuple[np.ndarray, np.ndarray]:
    ratio = np.divide(
        data.efficiency,
        mc.efficiency,
        out=np.full_like(data.efficiency, np.nan),
        where=np.isfinite(mc.efficiency) & (mc.efficiency > 0.0),
    )

    rel_data = np.divide(
        data.uncertainty,
        data.efficiency,
        out=np.zeros_like(data.efficiency),
        where=np.isfinite(data.efficiency) & (data.efficiency > 0.0),
    )
    rel_mc = np.divide(
        mc.uncertainty,
        mc.efficiency,
        out=np.zeros_like(mc.efficiency),
        where=np.isfinite(mc.efficiency) & (mc.efficiency > 0.0),
    )

    uncertainty = ratio * np.sqrt(rel_data * rel_data + rel_mc * rel_mc)
    return ratio, uncertainty


# =============================================================================
# Plotting
# =============================================================================

def save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    progress(f"PLOT: {path}")


def plot_pure_sample_inputs(
    source_name: str,
    arrays: Mapping[str, np.ndarray],
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    selected = diagnostics["selected_mask"]

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8))

    specs = [
        (
            diagnostics["all_score"],
            np.linspace(0.0, 1.0, 61),
            r"Tag BDT $\pi^0$ score",
            args.tag_score_min,
        ),
        (
            diagnostics["all_tag_energy"],
            np.linspace(0.0, 9.5, 76),
            r"$E_{\rm tag}$ (GeV)",
            args.tag_energy_min,
        ),
        (
            diagnostics["all_mx2_ep"],
            np.linspace(-0.5, 1.5, 100),
            r"$M_X^2(ep)$ (GeV$^2$)",
            args.mx2_ep_max,
        ),
        (
            diagnostics["all_missing_mass2"],
            np.linspace(-0.5, 0.5, 100),
            r"$M_X^2(ep\gamma_{\rm tag})$ (GeV$^2$)",
            None,
        ),
        (
            diagnostics["all_pred_energy"],
            np.linspace(0.0, 9.5, 76),
            r"Predicted $E_{\rm probe}$ (GeV)",
            args.probe_energy_min,
        ),
        (
            diagnostics["all_pred_theta"],
            np.linspace(0.0, 40.0, 80),
            r"Predicted $\theta_{\rm probe}$ (deg)",
            args.ft_theta_max,
        ),
    ]

    for ax, (values, bins, xlabel, marker) in zip(axes.flat, specs):
        finite = np.isfinite(values)

        if np.any(finite):
            ax.hist(
                values[finite],
                bins=bins,
                density=True,
                histtype="step",
                linewidth=1.4,
                label=f"Before aggressive cuts (N={np.sum(finite):,})",
            )
        #endif

        selected_values = values[selected & finite]

        if len(selected_values) > 0:
            ax.hist(
                selected_values,
                bins=bins,
                density=True,
                histtype="step",
                linewidth=1.7,
                label=f"Selected tags (N={len(selected_values):,})",
            )
        #endif

        if marker is not None:
            ax.axvline(marker, linestyle="--", linewidth=1.0)
        #endif

        ax.set_xlabel(xlabel)
        ax.set_ylabel("Normalized density")
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7)
    #endfor

    fig.suptitle(
        f"{source_name}: BDT-driven tag-and-probe denominator selection",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)



def print_geometry_hypothesis_summary(
    source_name: str,
    tag_index: int,
    results: Mapping[str, Mapping[str, np.ndarray]],
) -> None:
    progress(
        f"{source_name.upper()} EPGG tag-direction {tag_index} "
        "geometry hypotheses:"
    )

    for label, values in results.items():
        mask = values["mask"]
        opening = values["opening"][mask]
        dtheta = values["delta_theta"][mask]
        dphi = values["delta_phi"][mask]

        if len(opening) == 0:
            print(f"    {label:24s}: no finite trials")
            continue
        #endif

        print(
            f"    {label:24s}: "
            f"N={len(opening):9,d}, "
            f"median 3D={np.nanmedian(opening):8.3f} deg, "
            f"median |dtheta|={np.nanmedian(np.abs(dtheta)):8.3f} deg, "
            f"median |dphi|={np.nanmedian(np.abs(dphi)):8.3f} deg"
        )
    #endfor


def plot_geometry_hypotheses(
    source_name: str,
    tag_index: int,
    results: Mapping[str, Mapping[str, np.ndarray]],
    output_path: Path,
) -> None:
    """
    Unit-convention debugging plot.  No predicted FT/FD acceptance cut is
    applied here.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))

    for label, values in results.items():
        mask = values["mask"]
        opening = values["opening"][mask]
        dtheta = values["delta_theta"][mask]
        dphi = values["delta_phi"][mask]
        pred_theta = values["pred_theta"][mask]

        axes[0, 0].hist(
            opening,
            bins=np.linspace(0.0, 60.0, 181),
            histtype="step",
            linewidth=1.35,
            label=label,
        )
        axes[0, 1].hist(
            dtheta,
            bins=np.linspace(-60.0, 60.0, 241),
            histtype="step",
            linewidth=1.35,
            label=label,
        )
        axes[1, 0].hist(
            dphi,
            bins=np.linspace(-180.0, 180.0, 241),
            histtype="step",
            linewidth=1.35,
            label=label,
        )
        axes[1, 1].hist(
            pred_theta,
            bins=np.linspace(0.0, 60.0, 181),
            histtype="step",
            linewidth=1.35,
            label=label,
        )
    #endfor

    axes[0, 0].set_xlabel(
        r"$\angle(\gamma_{\rm pred},\gamma_{\rm reco})$ (deg)"
    )
    axes[0, 0].set_ylabel("Directed trials")
    axes[0, 0].set_yscale("log")

    axes[0, 1].set_xlabel(
        r"$\theta_{\rm pred}-\theta_{\rm reco}$ (deg)"
    )
    axes[0, 1].set_ylabel("Directed trials")
    axes[0, 1].set_yscale("log")

    axes[1, 0].set_xlabel(
        r"$\phi_{\rm pred}-\phi_{\rm reco}$ (deg)"
    )
    axes[1, 0].set_ylabel("Directed trials")
    axes[1, 0].set_yscale("log")

    axes[1, 1].set_xlabel(
        r"Predicted $\theta_{\rm probe}$ (deg)"
    )
    axes[1, 1].set_ylabel("Directed trials")
    axes[1, 1].set_yscale("log")

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=2,
        fontsize=8,
    )

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: geometry-unit hypotheses, tag direction {tag_index}",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.89])
    save_figure(fig, output_path)


def plot_raw_angle_branches(
    source_name: str,
    arrays: Mapping[str, np.ndarray],
    output_path: Path,
) -> None:
    """
    Show the raw stored angular branch ranges without any conversion.  This
    makes radians-vs-degrees immediately visible.
    """
    branches = [
        ("e_theta", r"$e'$ theta"),
        ("e_phi", r"$e'$ phi"),
        ("p1_theta", r"$p'$ theta"),
        ("p1_phi", r"$p'$ phi"),
        ("p2_theta", r"$\gamma_1$ theta"),
        ("p2_phi", r"$\gamma_1$ phi"),
        ("p3_theta", r"$\gamma_2$ theta"),
        ("p3_phi", r"$\gamma_2$ phi"),
    ]

    fig, axes = plt.subplots(2, 4, figsize=(16.0, 7.5))

    for ax, (branch, label) in zip(axes.flat, branches):
        values = np.asarray(arrays[branch], dtype=float)
        values = values[np.isfinite(values)]

        if len(values) > 0:
            lo, hi = np.nanpercentile(values, [0.2, 99.8])

            if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
                lo = np.nanmin(values)
                hi = np.nanmax(values)
            #endif

            ax.hist(
                values,
                bins=np.linspace(lo, hi, 120),
                histtype="step",
                linewidth=1.4,
            )
            ax.set_title(
                f"{branch}\n"
                f"p1={np.nanpercentile(values, 1):.3g}, "
                f"p50={np.nanmedian(values):.3g}, "
                f"p99={np.nanpercentile(values, 99):.3g}",
                fontsize=9,
            )
        #endif

        ax.set_xlabel(f"raw {label}")
        ax.set_ylabel("Rows")
        ax.set_yscale("log")
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: raw stored angular branches",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)



def binomial_efficiency(
    passed: np.ndarray,
    total: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    passed = np.asarray(passed, dtype=float)
    total = np.asarray(total, dtype=float)

    efficiency = np.full_like(total, np.nan, dtype=float)
    uncertainty = np.full_like(total, np.nan, dtype=float)

    valid = total > 0.0
    efficiency[valid] = passed[valid] / total[valid]
    uncertainty[valid] = np.sqrt(
        efficiency[valid]
        * (1.0 - efficiency[valid])
        / total[valid]
    )

    return efficiency, uncertainty


def binned_matching_efficiency(
    diagnostics: Mapping[str, np.ndarray],
    variable: np.ndarray,
    bin_edges: np.ndarray,
    region: int,
    angle_cut: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Conditional efficiency of the final angular association cut.

    Denominator:
      - predicted probe lies in requested FT/FD region,
      - reconstructed partner is in the same detector region,
      - reconstructed partner is above the nominal energy threshold.

    Numerator:
      denominator plus DeltaOmega(pred,reco) < angle_cut.

    This intentionally isolates ONLY the angular-matching requirement.
    """
    base = (
        diagnostics["mask_partner_energy"]
        & (diagnostics["pred_region"] == region)
        & np.isfinite(variable)
        & np.isfinite(diagnostics["match_angle"])
    )

    passed = base & (diagnostics["match_angle"] < angle_cut)

    total_counts, _ = np.histogram(
        variable[base],
        bins=bin_edges,
    )
    pass_counts, _ = np.histogram(
        variable[passed],
        bins=bin_edges,
    )

    efficiency, uncertainty = binomial_efficiency(
        pass_counts,
        total_counts,
    )

    return efficiency, uncertainty, total_counts


def binned_containment_angles(
    diagnostics: Mapping[str, np.ndarray],
    variable: np.ndarray,
    bin_edges: np.ndarray,
    region: int,
    quantiles: Tuple[float, ...] = (0.68, 0.90, 0.95, 0.99),
) -> Dict[float, np.ndarray]:
    """
    Angular containment of DeltaOmega(pred,reco) before applying the final
    angular-match requirement.
    """
    out = {
        q: np.full(len(bin_edges) - 1, np.nan, dtype=float)
        for q in quantiles
    }

    base = (
        diagnostics["mask_partner_energy"]
        & (diagnostics["pred_region"] == region)
        & np.isfinite(variable)
        & np.isfinite(diagnostics["match_angle"])
    )

    for ibin in range(len(bin_edges) - 1):
        in_bin = (
            base
            & (variable >= bin_edges[ibin])
            & (variable < bin_edges[ibin + 1])
        )
        angles = diagnostics["match_angle"][in_bin]

        if len(angles) < 20:
            continue
        #endif

        for q in quantiles:
            out[q][ibin] = np.quantile(angles, q)
        #endfor
    #endfor

    return out


def combine_direction_diagnostics(
    first: Mapping[str, np.ndarray],
    second: Mapping[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    """
    Concatenate only the event-wise arrays needed by the matching study.
    """
    keys = [
        "pred_energy",
        "pred_theta",
        "pred_region",
        "tag_energy",
        "partner_energy",
        "match_angle",
        "mask_partner_energy",
    ]

    return {
        key: np.concatenate(
            [
                np.asarray(first[key]),
                np.asarray(second[key]),
            ]
        )
        for key in keys
    }


def plot_matching_efficiency_study(
    source_name: str,
    numdiag1: Mapping[str, np.ndarray],
    numdiag2: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Matching efficiency versus predicted probe energy and theta.

    Direction 1 and direction 2 are kept separate, and their combination is
    shown as well.  FT and FD occupy separate columns.
    """
    energy_edges = np.asarray(
        [0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.5, 5.0, 7.0, 9.5],
        dtype=float,
    )
    theta_edges_by_region = {
        0: np.asarray([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]),
        1: np.asarray([5.5, 7.0, 9.0, 12.0, 16.0, 21.0, 27.0, 35.0]),
    }
    region_names = {0: "FT", 1: "FD"}

    combined = combine_direction_diagnostics(
        numdiag1,
        numdiag2,
    )
    samples = [
        ("tag direction 1", numdiag1),
        ("tag direction 2", numdiag2),
        ("combined", combined),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))

    for icol, region in enumerate((0, 1)):
        # Energy dependence.
        centers = 0.5 * (energy_edges[:-1] + energy_edges[1:])

        for label, diag in samples:
            eff, err, counts = binned_matching_efficiency(
                diag,
                np.asarray(diag["pred_energy"], dtype=float),
                energy_edges,
                region,
                args.probe_match_angle_max,
            )
            valid = np.isfinite(eff) & (counts > 0)

            axes[0, icol].errorbar(
                centers[valid],
                eff[valid],
                yerr=err[valid],
                marker="o",
                linestyle="-",
                linewidth=1.1,
                markersize=4,
                label=label,
            )
        #endfor

        axes[0, icol].set_title(
            f"{region_names[region]}: matching efficiency vs predicted energy"
        )
        axes[0, icol].set_xlabel(
            r"Predicted $E_{\rm probe}$ (GeV)"
        )
        axes[0, icol].set_ylabel(
            rf"$P(\Delta\Omega<{args.probe_match_angle_max:g}^\circ)$"
        )
        axes[0, icol].set_ylim(0.0, 1.05)

        # Theta dependence.
        theta_edges = theta_edges_by_region[region]
        theta_centers = 0.5 * (
            theta_edges[:-1] + theta_edges[1:]
        )

        for label, diag in samples:
            eff, err, counts = binned_matching_efficiency(
                diag,
                np.asarray(diag["pred_theta"], dtype=float),
                theta_edges,
                region,
                args.probe_match_angle_max,
            )
            valid = np.isfinite(eff) & (counts > 0)

            axes[1, icol].errorbar(
                theta_centers[valid],
                eff[valid],
                yerr=err[valid],
                marker="o",
                linestyle="-",
                linewidth=1.1,
                markersize=4,
                label=label,
            )
        #endfor

        axes[1, icol].set_title(
            f"{region_names[region]}: matching efficiency vs predicted theta"
        )
        axes[1, icol].set_xlabel(
            r"Predicted $\theta_{\rm probe}$ (deg)"
        )
        axes[1, icol].set_ylabel(
            rf"$P(\Delta\Omega<{args.probe_match_angle_max:g}^\circ)$"
        )
        axes[1, icol].set_ylim(0.0, 1.05)
    #endfor

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    #endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=3,
        fontsize=9,
    )
    fig.suptitle(
        f"{source_name}: conditional angular-matching efficiency",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.89])
    save_figure(fig, output_path)


def plot_matching_containment_study(
    source_name: str,
    numdiag1: Mapping[str, np.ndarray],
    numdiag2: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    68%, 90%, 95%, and 99% containment angle versus predicted E and theta.
    Uses both tag directions combined.
    """
    energy_edges = np.asarray(
        [0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.5, 5.0, 7.0, 9.5],
        dtype=float,
    )
    theta_edges_by_region = {
        0: np.asarray([2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]),
        1: np.asarray([5.5, 7.0, 9.0, 12.0, 16.0, 21.0, 27.0, 35.0]),
    }
    region_names = {0: "FT", 1: "FD"}
    quantiles = (0.68, 0.90, 0.95, 0.99)

    combined = combine_direction_diagnostics(
        numdiag1,
        numdiag2,
    )

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))

    for icol, region in enumerate((0, 1)):
        energy_centers = 0.5 * (
            energy_edges[:-1] + energy_edges[1:]
        )
        energy_q = binned_containment_angles(
            combined,
            combined["pred_energy"],
            energy_edges,
            region,
            quantiles,
        )

        for q in quantiles:
            values = energy_q[q]
            valid = np.isfinite(values)
            axes[0, icol].plot(
                energy_centers[valid],
                values[valid],
                marker="o",
                linewidth=1.2,
                markersize=4,
                label=f"{int(round(100*q))}% containment",
            )
        #endfor

        axes[0, icol].axhline(
            args.probe_match_angle_max,
            linestyle="--",
            linewidth=1.0,
            label=f"nominal {args.probe_match_angle_max:g} deg cut",
        )
        axes[0, icol].set_title(
            f"{region_names[region]}: angular containment vs energy"
        )
        axes[0, icol].set_xlabel(
            r"Predicted $E_{\rm probe}$ (GeV)"
        )
        axes[0, icol].set_ylabel(
            r"$\Delta\Omega$ containment angle (deg)"
        )

        theta_edges = theta_edges_by_region[region]
        theta_centers = 0.5 * (
            theta_edges[:-1] + theta_edges[1:]
        )
        theta_q = binned_containment_angles(
            combined,
            combined["pred_theta"],
            theta_edges,
            region,
            quantiles,
        )

        for q in quantiles:
            values = theta_q[q]
            valid = np.isfinite(values)
            axes[1, icol].plot(
                theta_centers[valid],
                values[valid],
                marker="o",
                linewidth=1.2,
                markersize=4,
                label=f"{int(round(100*q))}% containment",
            )
        #endfor

        axes[1, icol].axhline(
            args.probe_match_angle_max,
            linestyle="--",
            linewidth=1.0,
            label=f"nominal {args.probe_match_angle_max:g} deg cut",
        )
        axes[1, icol].set_title(
            f"{region_names[region]}: angular containment vs theta"
        )
        axes[1, icol].set_xlabel(
            r"Predicted $\theta_{\rm probe}$ (deg)"
        )
        axes[1, icol].set_ylabel(
            r"$\Delta\Omega$ containment angle (deg)"
        )
    #endfor

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    #endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=5,
        fontsize=8,
    )
    fig.suptitle(
        f"{source_name}: reconstructed-probe angular containment",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.89])
    save_figure(fig, output_path)


def plot_matching_angle_slices(
    source_name: str,
    numdiag1: Mapping[str, np.ndarray],
    numdiag2: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Direct DeltaOmega distributions in broad predicted-energy and theta bins.
    Uses both tag directions combined and shows FT/FD separately.
    """
    combined = combine_direction_diagnostics(
        numdiag1,
        numdiag2,
    )

    energy_slices = [
        (0.4, 1.0),
        (1.0, 2.0),
        (2.0, 9.5),
    ]
    theta_slices = {
        0: [(2.0, 3.2), (3.2, 4.4), (4.4, 5.5)],
        1: [(5.5, 10.0), (10.0, 20.0), (20.0, 35.0)],
    }
    region_names = {0: "FT", 1: "FD"}

    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))

    for icol, region in enumerate((0, 1)):
        base = (
            combined["mask_partner_energy"]
            & (combined["pred_region"] == region)
            & np.isfinite(combined["match_angle"])
        )

        for lo, hi in energy_slices:
            mask = (
                base
                & (combined["pred_energy"] >= lo)
                & (combined["pred_energy"] < hi)
            )
            angles = combined["match_angle"][mask]

            axes[0, icol].hist(
                angles,
                bins=np.linspace(0.0, 12.0, 121),
                histtype="step",
                linewidth=1.4,
                density=True,
                label=f"{lo:g}-{hi:g} GeV (N={len(angles):,})",
            )
        #endfor

        axes[0, icol].axvline(
            args.probe_match_angle_max,
            linestyle="--",
            linewidth=1.0,
            label=f"{args.probe_match_angle_max:g} deg cut",
        )
        axes[0, icol].set_title(
            f"{region_names[region]}: DeltaOmega by predicted energy"
        )
        axes[0, icol].set_xlabel(
            r"$\Delta\Omega_{\rm pred,reco}$ (deg)"
        )
        axes[0, icol].set_ylabel("Normalized density")
        axes[0, icol].set_yscale("log")
        axes[0, icol].legend(fontsize=7)

        for lo, hi in theta_slices[region]:
            mask = (
                base
                & (combined["pred_theta"] >= lo)
                & (combined["pred_theta"] < hi)
            )
            angles = combined["match_angle"][mask]

            axes[1, icol].hist(
                angles,
                bins=np.linspace(0.0, 12.0, 121),
                histtype="step",
                linewidth=1.4,
                density=True,
                label=f"{lo:g}-{hi:g} deg (N={len(angles):,})",
            )
        #endfor

        axes[1, icol].axvline(
            args.probe_match_angle_max,
            linestyle="--",
            linewidth=1.0,
            label=f"{args.probe_match_angle_max:g} deg cut",
        )
        axes[1, icol].set_title(
            f"{region_names[region]}: DeltaOmega by predicted theta"
        )
        axes[1, icol].set_xlabel(
            r"$\Delta\Omega_{\rm pred,reco}$ (deg)"
        )
        axes[1, icol].set_ylabel("Normalized density")
        axes[1, icol].set_yscale("log")
        axes[1, icol].legend(fontsize=7)
    #endfor

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: angular-resolution slices before final match cut",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_tag_direction_energy_ordering(
    source_name: str,
    epgg: Mapping[str, np.ndarray],
    numdiag1: Mapping[str, np.ndarray],
    numdiag2: Mapping[str, np.ndarray],
    output_path: Path,
) -> None:
    """
    Diagnose whether p2/p3 ordering explains the direction-dependent matching
    efficiency.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))

    p2 = np.asarray(epgg["p2_p"], dtype=float)
    p3 = np.asarray(epgg["p3_p"], dtype=float)
    finite = np.isfinite(p2) & np.isfinite(p3)

    axes[0, 0].hist2d(
        p2[finite],
        p3[finite],
        bins=[
            np.linspace(0.4, 8.0, 80),
            np.linspace(0.4, 8.0, 80),
        ],
    )
    axes[0, 0].plot(
        [0.4, 8.0],
        [0.4, 8.0],
        linestyle="--",
        linewidth=1.0,
    )
    axes[0, 0].set_xlabel(r"$E_{\gamma_1}=p2_p$ (GeV)")
    axes[0, 0].set_ylabel(r"$E_{\gamma_2}=p3_p$ (GeV)")
    axes[0, 0].set_title("Raw reconstructed photon-energy ordering")

    frac_p2_higher = (
        np.sum(finite & (p2 >= p3)) / np.sum(finite)
        if np.sum(finite) > 0
        else np.nan
    )
    axes[0, 0].text(
        0.04,
        0.96,
        f"P(p2 >= p3) = {frac_p2_higher:.3f}",
        transform=axes[0, 0].transAxes,
        va="top",
    )

    for label, diag in [
        ("tag direction 1: probe=p3", numdiag1),
        ("tag direction 2: probe=p2", numdiag2),
    ]:
        base = diag["mask_partner_energy"]
        energies = diag["pred_energy"][base]
        axes[0, 1].hist(
            energies,
            bins=np.linspace(0.4, 9.5, 80),
            histtype="step",
            linewidth=1.4,
            density=True,
            label=label,
        )

        angles = diag["match_angle"][base]
        axes[1, 0].hist(
            angles,
            bins=np.linspace(0.0, 12.0, 121),
            histtype="step",
            linewidth=1.4,
            density=True,
            label=label,
        )

        tag_energy = diag["tag_energy"][base]
        partner_energy = diag["partner_energy"][base]
        denom = tag_energy + partner_energy
        sharing = np.full_like(denom, np.nan, dtype=float)
        valid = np.isfinite(denom) & (denom > 0.0)
        sharing[valid] = partner_energy[valid] / denom[valid]

        axes[1, 1].hist(
            sharing[np.isfinite(sharing)],
            bins=np.linspace(0.0, 1.0, 80),
            histtype="step",
            linewidth=1.4,
            density=True,
            label=label,
        )
    #endfor

    axes[0, 1].set_xlabel(
        r"Predicted probe energy (GeV)"
    )
    axes[0, 1].set_ylabel("Normalized density")
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_title("Probe-energy distributions by tag direction")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].set_xlabel(
        r"$\Delta\Omega_{\rm pred,reco}$ (deg)"
    )
    axes[1, 0].set_ylabel("Normalized density")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Angular resolution by tag direction")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].set_xlabel(
        r"$E_{\rm probe}/(E_{\rm tag}+E_{\rm probe})$"
    )
    axes[1, 1].set_ylabel("Normalized density")
    axes[1, 1].set_title("Photon energy sharing by tag direction")
    axes[1, 1].legend(fontsize=8)

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: tag-direction ordering diagnostics",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def print_matching_efficiency_summary(
    source_name: str,
    numdiag1: Mapping[str, np.ndarray],
    numdiag2: Mapping[str, np.ndarray],
    args: argparse.Namespace,
) -> None:
    progress(
        f"{source_name.upper()} conditional angular-matching summary:"
    )

    for region, region_name in [(0, "FT"), (1, "FD")]:
        for label, diag in [
            ("direction 1", numdiag1),
            ("direction 2", numdiag2),
        ]:
            base = (
                diag["mask_partner_energy"]
                & (diag["pred_region"] == region)
                & np.isfinite(diag["match_angle"])
            )
            n_total = int(np.sum(base))
            n_pass = int(
                np.sum(
                    base
                    & (
                        diag["match_angle"]
                        < args.probe_match_angle_max
                    )
                )
            )
            efficiency = (
                n_pass / n_total
                if n_total > 0
                else np.nan
            )

            print(
                f"    {region_name:2s} {label:11s}: "
                f"{n_pass:9,d}/{n_total:9,d} = "
                f"{efficiency:7.4f}"
            )
        #endfor
    #endfor



def plot_pre_match_probe_diagnostics(
    source_name: str,
    diagnostics: Mapping[str, np.ndarray],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    """
    Plot predicted-vs-reconstructed probe behavior BEFORE predicted angular
    acceptance and before the final angular-match requirement.
    """
    base = diagnostics["mask_pre_acceptance"]
    accepted = diagnostics["mask_probe_acceptance"]

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9.0))

    values_all = diagnostics["match_angle"][base]
    values_all = values_all[np.isfinite(values_all)]
    axes[0, 0].hist(
        values_all,
        bins=np.linspace(0.0, 60.0, 181),
        histtype="step",
        linewidth=1.5,
        label=f"all pre-acceptance (N={len(values_all):,})",
    )

    values_acc = diagnostics["match_angle"][accepted]
    values_acc = values_acc[np.isfinite(values_acc)]
    axes[0, 0].hist(
        values_acc,
        bins=np.linspace(0.0, 60.0, 181),
        histtype="step",
        linewidth=1.5,
        label=f"predicted FT/FD accepted (N={len(values_acc):,})",
    )

    axes[0, 0].axvline(
        args.probe_match_angle_max,
        linestyle="--",
        linewidth=1.0,
        label=f"nominal match < {args.probe_match_angle_max:g} deg",
    )
    axes[0, 0].set_xlabel(
        r"$\angle(\gamma_{\rm pred},\gamma_{\rm reco})$ (deg)"
    )
    axes[0, 0].set_ylabel("Directed trials")
    axes[0, 0].set_yscale("log")
    axes[0, 0].legend(fontsize=8)

    finite = (
        base
        & np.isfinite(diagnostics["pred_energy"])
        & np.isfinite(diagnostics["partner_energy"])
    )
    x = diagnostics["pred_energy"][finite]
    y = diagnostics["partner_energy"][finite]

    axes[0, 1].hist2d(
        x,
        y,
        bins=[np.linspace(0.4, 9.5, 70), np.linspace(0.4, 9.5, 70)],
    )
    axes[0, 1].plot(
        [0.4, 9.5],
        [0.4, 9.5],
        linestyle="--",
        linewidth=1.0,
    )
    axes[0, 1].set_xlabel(r"Predicted $E_{\rm probe}$ (GeV)")
    axes[0, 1].set_ylabel(r"Reconstructed $E_{\rm probe}$ (GeV)")

    finite_res = finite & (diagnostics["pred_energy"] > 0.0)
    residual = (
        diagnostics["partner_energy"][finite_res]
        - diagnostics["pred_energy"][finite_res]
    ) / diagnostics["pred_energy"][finite_res]

    axes[1, 0].hist(
        residual[np.isfinite(residual)],
        bins=np.linspace(-1.0, 1.0, 120),
        histtype="step",
        linewidth=1.5,
    )
    axes[1, 0].set_xlabel(
        r"$(E_{\rm reco}-E_{\rm pred})/E_{\rm pred}$"
    )
    axes[1, 0].set_ylabel("Directed trials")

    pair_mask = base & np.isfinite(diagnostics["pair_mass"])
    axes[1, 1].hist(
        diagnostics["pair_mass"][pair_mask],
        bins=np.linspace(0.0, 0.8, 120),
        histtype="step",
        linewidth=1.5,
    )
    axes[1, 1].axvspan(
        PI0_MASS_WINDOW[0],
        PI0_MASS_WINDOW[1],
        alpha=0.12,
        label=r"$\pi^0$ mass window (diagnostic only)",
    )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    axes[1, 1].set_ylabel("Directed trials")
    axes[1, 1].legend(fontsize=8)

    for ax in axes.flat:
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: probe prediction validation before final matching",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)



def plot_partner_matching(
    source_name: str,
    numerator: DirectedSample,
    output_path: Path,
) -> None:
    if numerator.partner_match_angle is None:
        return
    #endif

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6))

    axes[0].hist(
        numerator.partner_match_angle,
        bins=np.linspace(0.0, 3.0, 61),
        histtype="step",
        linewidth=1.5,
    )
    axes[0].set_xlabel(
        r"$\angle(\gamma_{\rm pred},\gamma_{\rm reco})$ (deg)"
    )

    if numerator.partner_energy is not None:
        residual = (
            numerator.partner_energy - numerator.predicted_energy
        ) / numerator.predicted_energy
        finite = np.isfinite(residual)

        axes[1].hist(
            residual[finite],
            bins=np.linspace(-0.6, 0.6, 100),
            histtype="step",
            linewidth=1.5,
        )
        axes[1].set_xlabel(
            r"$(E_{\rm reco}-E_{\rm pred})/E_{\rm pred}$"
        )
    #endif

    if numerator.pair_mass is not None:
        axes[2].hist(
            numerator.pair_mass,
            bins=np.linspace(0.0, 0.8, 120),
            histtype="step",
            linewidth=1.5,
        )
        axes[2].axvspan(
            PI0_MASS_WINDOW[0],
            PI0_MASS_WINDOW[1],
            alpha=0.12,
            label=r"$\pi^0$ mass window (diagnostic only)",
        )
        axes[2].set_yscale("log")
        axes[2].set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
        axes[2].legend(fontsize=8)
    #endif

    for ax in axes:
        ax.set_ylabel("Directed trials")
        ax.grid(alpha=0.2)
    #endfor

    fig.suptitle(
        f"{source_name}: reconstructed-probe matching diagnostics",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_efficiencies(
    curves: Mapping[str, Mapping[str, EfficiencyCurve]],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), sharex="col")

    for column, region_name in enumerate(["FT", "FD"]):
        ax_eff = axes[0, column]
        ax_ratio = axes[1, column]

        for source_name in ["data", "aaogen", "clasdis"]:
            curve = curves[source_name][region_name]

            ax_eff.errorbar(
                curve.centers,
                curve.efficiency,
                yerr=curve.uncertainty,
                fmt="o-",
                markersize=4,
                linewidth=1.3,
                label=source_name.upper(),
            )
        #endfor

        data_curve = curves["data"][region_name]
        aaogen_curve = curves["aaogen"][region_name]
        clasdis_curve = curves["clasdis"][region_name]

        ratio_aao, ratio_aao_err = ratio_curve(
            data_curve,
            aaogen_curve,
        )
        ratio_dis, ratio_dis_err = ratio_curve(
            data_curve,
            clasdis_curve,
        )

        ax_ratio.errorbar(
            data_curve.centers,
            ratio_aao,
            yerr=ratio_aao_err,
            fmt="o-",
            markersize=4,
            linewidth=1.3,
            label="Data / AAOgen",
        )
        ax_ratio.errorbar(
            data_curve.centers,
            ratio_dis,
            yerr=ratio_dis_err,
            fmt="s--",
            markersize=4,
            linewidth=1.2,
            label="Data / CLASDIS",
        )

        ax_eff.set_title(region_name)
        ax_eff.set_ylabel(r"$\epsilon_\gamma$")
        ax_eff.set_ylim(0.0, 1.10)
        ax_eff.grid(alpha=0.2)
        ax_eff.legend(fontsize=8)

        ax_ratio.axhline(1.0, linestyle="--", linewidth=1.0)
        ax_ratio.set_xlabel(r"Predicted $E_{\rm probe}$ (GeV)")
        ax_ratio.set_ylabel(r"$\epsilon_{\rm data}/\epsilon_{\rm MC}$")
        ax_ratio.set_ylim(0.5, 1.5)
        ax_ratio.grid(alpha=0.2)
        ax_ratio.legend(fontsize=8)
    #endfor

    fig.suptitle(
        "First-stage aggressive pi0-BDT tag-and-probe photon efficiency",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)


def plot_counts(
    curves: Mapping[str, Mapping[str, EfficiencyCurve]],
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.0), sharex=True)

    for row, region_name in enumerate(["FT", "FD"]):
        for source_name in ["data", "aaogen", "clasdis"]:
            curve = curves[source_name][region_name]

            axes[row, 0].step(
                curve.centers,
                curve.denominator,
                where="mid",
                label=source_name.upper(),
            )
            axes[row, 1].step(
                curve.centers,
                curve.numerator,
                where="mid",
                label=source_name.upper(),
            )
        #endfor

        axes[row, 0].set_title(f"{region_name}: predicted probes")
        axes[row, 1].set_title(f"{region_name}: reconstructed probes")

        for ax in axes[row]:
            ax.set_yscale("log")
            ax.set_ylabel("Directed trials")
            ax.grid(alpha=0.2)
            ax.legend(fontsize=8)
        #endfor
    #endfor

    for ax in axes[-1]:
        ax.set_xlabel(r"Predicted $E_{\rm probe}$ (GeV)")
    #endfor

    fig.suptitle(
        "Tag-and-probe raw numerator / denominator statistics",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.96])
    save_figure(fig, output_path)


def plot_clasdis_tag_truth_purity(
    denominator: DirectedSample,
    output_path: Path,
    edges: np.ndarray,
) -> None:
    if denominator.tag_parent_pid is None:
        return
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)

    for ax, (region_name, region_value) in zip(axes, REGIONS.items()):
        mask = denominator.region == region_value
        energy = denominator.predicted_energy[mask]
        parent = denominator.tag_parent_pid[mask]

        total, _ = np.histogram(energy, bins=edges)
        true_pi0, _ = np.histogram(
            energy[parent == 111],
            bins=edges,
        )

        purity = np.divide(
            true_pi0,
            total,
            out=np.full(len(total), np.nan, dtype=float),
            where=total > 0,
        )

        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.plot(centers, purity, "o-")
        ax.set_title(region_name)
        ax.set_xlabel(r"Predicted $E_{\rm probe}$ (GeV)")
        ax.set_ylim(0.0, 1.05)
        ax.grid(alpha=0.2)
    #endfor

    axes[0].set_ylabel(
        r"Fraction of selected tags with truth parent $\pi^0$"
    )

    fig.suptitle(
        "CLASDIS truth check of aggressive tag-sample purity",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


def plot_threshold_scan(
    arrays_by_source: Mapping[str, Mapping[str, np.ndarray]],
    models: Mapping[str, object],
    beam_energy: float,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    """
    Integrated diagnostic only: how aggressive tag-score selection changes
    retained data statistics and CLASDIS truth purity.
    """
    thresholds = np.asarray([0.60, 0.70, 0.75, 0.80, 0.85, 0.90])
    original = args.tag_score_min

    data_counts = {"FT": [], "FD": []}
    clasdis_purity = {"FT": [], "FD": []}

    for threshold in thresholds:
        args.tag_score_min = float(threshold)

        data_den, _ = denominator_from_epg(
            arrays_by_source["data"],
            models,
            beam_energy,
            args,
        )
        dis_den, _ = denominator_from_epg(
            arrays_by_source["clasdis"],
            models,
            beam_energy,
            args,
        )

        for region_name, region_value in REGIONS.items():
            data_counts[region_name].append(
                np.sum(data_den.region == region_value)
            )

            mask = dis_den.region == region_value

            if (
                dis_den.tag_parent_pid is not None
                and np.sum(mask) > 0
            ):
                purity = np.mean(
                    dis_den.tag_parent_pid[mask] == 111
                )
            else:
                purity = np.nan
            #endif

            clasdis_purity[region_name].append(purity)
        #endfor
    #endfor

    args.tag_score_min = original

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))

    for region_name in ["FT", "FD"]:
        counts = np.asarray(data_counts[region_name], dtype=float)

        if counts[0] > 0:
            retained = counts / counts[0]
        else:
            retained = np.full_like(counts, np.nan)
        #endif

        axes[0].plot(
            thresholds,
            retained,
            "o-",
            label=region_name,
        )
        axes[1].plot(
            thresholds,
            clasdis_purity[region_name],
            "o-",
            label=region_name,
        )
    #endfor

    axes[0].set_xlabel("Minimum tag BDT score")
    axes[0].set_ylabel("Data denominator retained / value at 0.60")
    axes[1].set_xlabel("Minimum tag BDT score")
    axes[1].set_ylabel(
        r"CLASDIS selected-tag truth-$\pi^0$ fraction"
    )

    for ax in axes:
        ax.grid(alpha=0.2)
        ax.legend()
    #endfor

    fig.suptitle(
        "Aggressiveness scan for the tag definition",
        y=0.99,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])
    save_figure(fig, output_path)


# =============================================================================
# Source workflow
# =============================================================================

def process_source(
    source_name: str,
    epg_path: Path,
    epgg_path: Path,
    models: Mapping[str, object],
    beam_energy: float,
    args: argparse.Namespace,
    output_dir: Path,
) -> Tuple[
    DirectedSample,
    DirectedSample,
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
]:
    truth_epg = EPG_OPTIONAL_TRUTH if source_name != "data" else []
    truth_epgg = EPGG_OPTIONAL_TRUTH if source_name != "data" else []

    epg = read_arrays(
        epg_path,
        EPG_REQUIRED,
        truth_epg,
    )
    epgg = read_arrays(
        epgg_path,
        EPGG_REQUIRED,
        truth_epgg,
    )

    denominator, diagnostics = denominator_from_epg(
        epg,
        models,
        beam_energy,
        args,
    )

    num1, numdiag1 = numerator_direction_from_epgg(
        epgg,
        models,
        beam_energy,
        args,
        tag_index=1,
    )
    num2, numdiag2 = numerator_direction_from_epgg(
        epgg,
        models,
        beam_energy,
        args,
        tag_index=2,
    )

    print_numerator_cutflow(source_name, numdiag1)
    print_numerator_cutflow(source_name, numdiag2)

    numerator = concatenate_directed(num1, num2)
    numerator_diagnostics = concatenate_numerator_diagnostics(
        numdiag1,
        numdiag2,
    )

    # The earlier geometry, closure, association-plateau, and raw-input
    # canvases served their development purpose and are no longer produced.
    # Keep the numerical cut flow above, but reserve output plots for the
    # diagnostics that directly inform the final tag-and-probe correction.

    progress(
        f"{source_name.upper()}: denominator directed probes="
        f"{len(denominator.predicted_energy):,}; "
        f"numerator reconstructed probes="
        f"{len(numerator.predicted_energy):,}"
    )

    return (
        denominator,
        numerator,
        diagnostics,
        epg,
        numerator_diagnostics,
    )


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    edges = np.asarray(
        [float(token) for token in args.energy_bins.split(",")],
        dtype=float,
    )

    if len(edges) < 2 or np.any(np.diff(edges) <= 0.0):
        raise ValueError("--energy-bins must be strictly increasing.")
    #endif

    paths = build_paths(args)
    model_dir = default_model_dir(args)
    output_dir = default_output_dir(args)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 92)
    print("FIRST-STAGE PI0-BDT TAG-AND-PROBE PHOTON EFFICIENCY")
    print(
        f"period={args.period}, sample_tag={args.sample_tag}, "
        f"Ebeam={BEAM_ENERGIES[args.period]:.3f} GeV"
    )
    print(
        "Aggressive tag definition: "
        f"score>={args.tag_score_min:.2f}, "
        f"Etag>={args.tag_energy_min:.2f} GeV, "
        f"{args.mx2_ep_min:.3f}<=Mx2(ep)<{args.mx2_ep_max:.3f} GeV^2; "
        "NO nominal Mx2(epgamma_tag) cut"
    )
    print(
        "Probe definition: "
        f"{args.probe_energy_min:.2f}<=Epred<{args.probe_energy_max:.2f} GeV, "
        f"partner angular match<{args.probe_match_angle_max:.1f} deg"
    )
    print(
        "Fit numerator: common-mean double Gaussian pi0 signal + "
        "quadratic background; nominal association=10 deg, "
        "variations=8/12 deg"
    )
    print("=" * 92)

    for key, path in paths.items():
        print(f"INPUT {key}: {path}")

        if not path.exists():
            raise FileNotFoundError(path)
        #endif
    #endfor

    models = load_models(model_dir)
    beam_energy = BEAM_ENERGIES[args.period]

    denominators: Dict[str, DirectedSample] = {}
    numerators: Dict[str, DirectedSample] = {}
    epg_arrays: Dict[str, Mapping[str, np.ndarray]] = {}
    numerator_diagnostics_by_source: Dict[
        str,
        Dict[str, np.ndarray],
    ] = {}
    denominator_diagnostics_by_source: Dict[
        str,
        Dict[str, np.ndarray],
    ] = {}

    for source_name in ["data", "aaogen", "clasdis"]:
        (
            denominator,
            numerator,
            denominator_diagnostics,
            epg,
            numerator_diagnostics,
        ) = process_source(
            source_name,
            paths[f"{source_name}_epg"],
            paths[f"{source_name}_epgg"],
            models,
            beam_energy,
            args,
            output_dir,
        )
        denominators[source_name] = denominator
        numerators[source_name] = numerator
        epg_arrays[source_name] = epg
        denominator_diagnostics_by_source[source_name] = (
            denominator_diagnostics
        )
        numerator_diagnostics_by_source[source_name] = (
            numerator_diagnostics
        )
    #endfor

    curves: Dict[str, Dict[str, EfficiencyCurve]] = {}

    for source_name in ["data", "aaogen", "clasdis"]:
        curves[source_name] = {}

        for region_name, region_value in REGIONS.items():
            curves[source_name][region_name] = efficiency_curve(
                denominators[source_name],
                numerators[source_name],
                region_value,
                edges,
            )

            curve = curves[source_name][region_name]
            finite = np.isfinite(curve.efficiency)

            if np.any(
                finite
                & (
                    (curve.efficiency < 0.0)
                    | (curve.efficiency > 1.0)
                )
            ):
                progress(
                    f"WARNING: {source_name} {region_name} has raw "
                    "numerator/denominator bins outside [0,1]. "
                    "This indicates pair combinatorics or inconsistent "
                    "topology exposure and must be investigated."
                )
            #endif
        #endfor
    #endfor

    plot_clasdis_tag_truth_purity(
        denominators["clasdis"],
        output_dir / "05_clasdis_tag_truth_purity.png",
        edges,
    )

    # ------------------------------------------------------------------
    # M(gamma gamma)-fit numerator.
    # Both numerator and denominator remain binned in E_pred.  The only
    # conceptual change here is replacing raw reconstructed-pair counts by
    # the fitted pi0 signal yield.
    # ------------------------------------------------------------------
    fitted_curves_by_source_angle: Dict[
        str,
        Dict[float, Dict[str, FitEfficiencyCurve]],
    ] = {}

    nominal_fit_curves: Dict[
        str,
        Dict[str, FitEfficiencyCurve],
    ] = {}

    for source_name in ["data", "aaogen", "clasdis"]:
        fitted_curves_by_source_angle[source_name] = {}

        for angle_max in MGG_FIT_ASSOCIATION_ANGLES:
            fitted_curves_by_source_angle[source_name][angle_max] = {}
            fits_for_nominal: Dict[str, list[MassFitResult]] = {}

            for region_name, region_value in REGIONS.items():
                curve, fits = fit_efficiency_curve(
                    denominators[source_name],
                    numerator_diagnostics_by_source[source_name],
                    region_value,
                    edges,
                    angle_max,
                    args,
                )
                fitted_curves_by_source_angle[
                    source_name
                ][angle_max][region_name] = curve

                if np.isclose(
                    angle_max,
                    MGG_FIT_NOMINAL_ANGLE,
                ):
                    fits_for_nominal[region_name] = fits
                #endif
            #endfor

            if np.isclose(
                angle_max,
                MGG_FIT_NOMINAL_ANGLE,
            ):
                nominal_fit_curves[source_name] = (
                    fitted_curves_by_source_angle[
                        source_name
                    ][angle_max]
                )
                plot_mgg_fit_examples(
                    source_name.upper(),
                    fits_for_nominal,
                    edges,
                    output_dir
                    / f"11_{source_name}_mgg_fit_examples_10deg.png",
                    angle_max,
                )
            #endif
        #endfor

        if source_name == "data":
            plot_fit_efficiency_angle_stability(
                source_name.upper(),
                fitted_curves_by_source_angle[source_name],
                output_dir
                / "12_data_fitted_efficiency_angle_stability.png",
            )
        #endif
    #endfor

    plot_fit_efficiency_data_mc_ratio(
        nominal_fit_curves,
        output_dir / "14_mgg_fit_efficiency_and_data_mc_ratio.png",
    )

    # ------------------------------------------------------------------
    # BDT-threshold diagnostic.
    # The DATA/MC ratio scan is retained as the primary working-point test.
    # Earlier raw-retention and apparent-efficiency canvases are no longer
    # produced because they are superseded by the ratio diagnostic.
    # ------------------------------------------------------------------
    progress(
        "Building apparent-efficiency scan versus minimum BDT tag score..."
    )
    bdt_efficiency_scan = build_bdt_efficiency_threshold_scan(
        denominator_diagnostics_by_source,
        numerator_diagnostics_by_source,
        edges,
        args,
    )
    plot_bdt_threshold_data_mc_ratio(
        bdt_efficiency_scan,
        edges,
        output_dir / "17_bdt_threshold_data_mc_ratio.png",
    )

    progress(
        "Building AAOgen truth-matched-tag BDT-threshold diagnostic..."
    )
    aaogen_truth_tag_scan = build_aaogen_truth_tag_threshold_scan(
        denominator_diagnostics_by_source["aaogen"],
        numerator_diagnostics_by_source["aaogen"],
        edges,
        args,
    )
    plot_aaogen_truth_tag_threshold_scan(
        bdt_efficiency_scan,
        aaogen_truth_tag_scan,
        edges,
        output_dir / "18_aaogen_truth_tag_bdt_threshold_scan.png",
    )
    plot_aaogen_truth_tag_score_phase_space(
        denominator_diagnostics_by_source["aaogen"],
        edges,
        args,
        output_dir / "20_aaogen_truth_tag_score_phase_space.png",
    )

    progress(
        "Building coarse Epred x theta_pred photon-efficiency diagnostics..."
    )
    theta_binned_scan = build_theta_binned_bdt_scan(
        denominator_diagnostics_by_source,
        numerator_diagnostics_by_source,
        edges,
        args,
    )
    plot_theta_binned_nominal_data_mc_ratio(
        theta_binned_scan,
        edges,
        args,
        output_dir / "21_theta_binned_nominal_data_mc_ratio.png",
    )
    plot_fd_theta_binned_bdt_ratio_scan(
        theta_binned_scan,
        edges,
        args,
        output_dir / "22_fd_theta_binned_bdt_ratio_scan.png",
    )

    progress(
        "Building nominal FD Epred x theta_pred x sector corrections..."
    )
    fd_sector_scan = build_fd_theta_sector_nominal_scan(
        denominator_diagnostics_by_source,
        numerator_diagnostics_by_source,
        edges,
        args,
    )
    plot_fd_sector_dependence_nominal(
        fd_sector_scan,
        edges,
        args,
        output_dir / "23_fd_sector_dependence_nominal.png",
    )
    plot_fd_sector_averaged_correction(
        fd_sector_scan,
        edges,
        args,
        output_dir / "24_fd_sector_specific_corrections.png",
    )
    plot_predicted_vs_reconstructed_sector_migration(
        numerator_diagnostics_by_source,
        args,
        output_dir / "25_fd_sector_migration.png",
    )
    plot_fd_sector_statistics(
        fd_sector_scan,
        edges,
        args,
        output_dir / "26_fd_sector_statistics.png",
    )
    plot_fd_sector_raw_efficiencies(
        fd_sector_scan,
        edges,
        args,
        output_dir / "27_fd_sector_raw_efficiencies.png",
    )

    print("\n" + "=" * 92)
    print("SUMMARY")
    print("=" * 92)

    for region_name in ["FT", "FD"]:
        print(f"\n{region_name}")

        for source_name in ["data", "aaogen", "clasdis"]:
            curve = curves[source_name][region_name]
            den_total = int(np.sum(curve.denominator))
            num_total = int(np.sum(curve.numerator))
            raw = num_total / den_total if den_total > 0 else np.nan
            print(
                f"  {source_name:8s}: denominator={den_total:10,d}  "
                f"numerator={num_total:10,d}  raw integrated={raw:.4f}"
            )
        #endfor
    #endfor

    print("\nNominal 10-degree Mgg-fit efficiencies:")

    for region_name in ["FT", "FD"]:
        print(f"  {region_name}")

        for source_name in ["data", "aaogen", "clasdis"]:
            curve = nominal_fit_curves[source_name][region_name]
            values = []

            for ibin in range(len(curve.centers)):
                values.append(
                    f"{curve.edges[ibin]:g}-{curve.edges[ibin+1]:g}:"
                    f"{curve.efficiency[ibin]:.4f}"
                )
            #endfor

            print(
                f"    {source_name:8s}: "
                + ", ".join(values)
            )
        #endfor
    #endfor

    print("\nBDT-threshold apparent-efficiency scan:")
    for region_name in ["FT", "FD"]:
        print(f"  {region_name}")
        for source_name in ["data", "aaogen", "clasdis"]:
            block = bdt_efficiency_scan[source_name][region_name]
            values = []
            for ithr, threshold in enumerate(
                BDT_EFFICIENCY_SCAN_THRESHOLDS
            ):
                finite = np.isfinite(block["efficiency"][ithr])
                if np.any(finite):
                    den = np.sum(
                        block["denominator"][ithr][finite]
                    )
                    num = np.nansum(
                        block["numerator"][ithr][finite]
                    )
                    integrated = (
                        num / den if den > 0.0 else np.nan
                    )
                else:
                    integrated = np.nan
                #endif
                values.append(
                    f"{threshold:.3f}:{integrated:.4f}"
                )
            #endfor
            print(
                f"    {source_name:8s}: "
                + ", ".join(values)
            )
        #endfor
    #endfor

    print("\nBDT-threshold DATA/MC ratio scan:")
    for region_name in ["FT", "FD"]:
        print(f"  {region_name}")
        for ibin in range(len(edges) - 1):
            label = (
                f"{edges[ibin]:g}-{edges[ibin+1]:g} GeV"
            )
            pieces = []
            for mc_name in ["aaogen", "clasdis"]:
                data_eff = (
                    bdt_efficiency_scan["data"][region_name]
                    ["efficiency"][:, ibin]
                )
                mc_eff = (
                    bdt_efficiency_scan[mc_name][region_name]
                    ["efficiency"][:, ibin]
                )
                ratio = np.divide(
                    data_eff,
                    mc_eff,
                    out=np.full_like(data_eff, np.nan),
                    where=(
                        np.isfinite(data_eff)
                        & np.isfinite(mc_eff)
                        & (mc_eff > 0.0)
                    ),
                )
                values = ", ".join(
                    f"{thr:.3f}:{val:.3f}"
                    for thr, val in zip(
                        BDT_EFFICIENCY_SCAN_THRESHOLDS,
                        ratio,
                    )
                    if np.isfinite(val)
                )
                pieces.append(
                    f"DATA/{mc_name.upper()} [{values}]"
                )
            #endfor
            print(f"    {label}: " + " ; ".join(pieces))
        #endfor
    #endfor

    print("\nAAOgen truth-matched-tag threshold scan:")
    for region_name in ["FT", "FD"]:
        print(f"  {region_name}")
        for ibin in range(len(edges) - 1):
            standard = (
                bdt_efficiency_scan["aaogen"][region_name]
                ["efficiency"][:, ibin]
            )
            truth = (
                aaogen_truth_tag_scan[region_name]
                ["efficiency"][:, ibin]
            )
            pieces = []
            for ithr, threshold in enumerate(
                BDT_EFFICIENCY_SCAN_THRESHOLDS
            ):
                if np.isfinite(standard[ithr]) and np.isfinite(truth[ithr]):
                    pieces.append(
                        f"{threshold:.3f}:"
                        f"{standard[ithr]:.3f}->{truth[ithr]:.3f}"
                    )
                #endif
            #endfor
            print(
                f"    {edges[ibin]:g}-{edges[ibin+1]:g} GeV: "
                + ", ".join(pieces)
            )
        #endfor
    #endfor

    print("\nNominal coarse Epred x theta_pred DATA/MC correction:")
    nominal_index = int(
        np.argmin(
            np.abs(
                BDT_EFFICIENCY_SCAN_THRESHOLDS
                - float(args.tag_score_min)
            )
        )
    )
    for label, _, _, _ in probe_theta_bins(args):
        print(f"  {label}")
        data_block = theta_binned_scan["data"][label]
        for mc_name in ["aaogen", "clasdis"]:
            mc_block = theta_binned_scan[mc_name][label]
            data_eff = data_block["efficiency"][nominal_index]
            mc_eff = mc_block["efficiency"][nominal_index]
            ratio = np.divide(
                data_eff,
                mc_eff,
                out=np.full_like(data_eff, np.nan),
                where=(
                    np.isfinite(data_eff)
                    & np.isfinite(mc_eff)
                    & (mc_eff > 0.0)
                ),
            )
            values = ", ".join(
                f"{edges[i]:g}-{edges[i+1]:g}:{ratio[i]:.3f}"
                for i in range(len(ratio))
                if np.isfinite(ratio[i])
            )
            print(f"    DATA/{mc_name.upper()}: {values}")
        #endfor
    #endfor

    print("\nNominal FD sector-specific DATA/MC correction:")
    for label, region_value, theta_low, theta_high in probe_theta_bins(args):
        if region_value != REGIONS["FD"]:
            continue
        #endif
        print(f"  {label}")
        for sector in range(1, 7):
            pieces = []
            for mc_name in ["aaogen", "clasdis"]:
                data_eff = (
                    fd_sector_scan["data"][label][sector]["efficiency"]
                )
                mc_eff = (
                    fd_sector_scan[mc_name][label][sector]["efficiency"]
                )
                ratio = np.divide(
                    data_eff,
                    mc_eff,
                    out=np.full_like(data_eff, np.nan),
                    where=(
                        np.isfinite(data_eff)
                        & np.isfinite(mc_eff)
                        & (mc_eff > 0.0)
                    ),
                )
                values = ", ".join(
                    f"{edges[i]:g}-{edges[i+1]:g}:{ratio[i]:.3f}"
                    for i in range(len(ratio))
                    if np.isfinite(ratio[i])
                )
                pieces.append(
                    f"DATA/{mc_name.upper()} [{values}]"
                )
            #endfor
            print(f"    S{sector}: " + " ; ".join(pieces))
        #endfor
    #endfor

    print("\nOutputs retained after diagnostic pruning:")
    print(f"  {output_dir / '05_clasdis_tag_truth_purity.png'}")
    print("  11_<source>_mgg_fit_examples_10deg.png")
    print("  12_data_fitted_efficiency_angle_stability.png")
    print("  14_mgg_fit_efficiency_and_data_mc_ratio.png")
    print("  17_bdt_threshold_data_mc_ratio.png")
    print("  18_aaogen_truth_tag_bdt_threshold_scan.png")
    print("  20_aaogen_truth_tag_score_phase_space.png")
    print("  21_theta_binned_nominal_data_mc_ratio.png")
    print("  22_fd_theta_binned_bdt_ratio_scan.png")
    print("  23_fd_sector_dependence_nominal.png")
    print("  24_fd_sector_specific_corrections.png")
    print("  25_fd_sector_migration.png")
    print("  26_fd_sector_statistics.png")
    print("  27_fd_sector_raw_efficiencies.png")
    print("\nThis is a first-stage diagnostic extraction, not yet a final correction.")


if __name__ == "__main__":
    main()
