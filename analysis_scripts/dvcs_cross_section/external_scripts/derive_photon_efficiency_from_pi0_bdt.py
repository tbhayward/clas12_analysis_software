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
    [0.40, 1.00, 2.00, 9.50],
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
        "selected_mask": base,
        # Compare the branch Mx2(epgamma) against a direct four-vector check.
        "computed_missing_mass2": pred["mass2"],
    }
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
    elif tag_index == 2:
        tag_p = "p3_p"
        tag_theta = "p3_theta"
        tag_phi = "p3_phi"
        probe_p = "p2_p"
        probe_theta = "p2_theta"
        probe_phi = "p2_phi"
        probe_detector = "detector_gamma1"
        tag_parent_branch = "gamma2_parent_pid"
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

    diagnostics = {
        "tag_index": np.asarray([tag_index], dtype=np.int64),
        "scores": scores,
        "tag_energy": tag_energy,
        "mx2_ep": mx2_ep,
        "missing_mass2": branch_missing_mass2,
        "pred_energy": pred["energy"],
        "pred_theta": pred["theta"],
        "pred_region": predicted_region,
        "actual_detector": actual_detector,
        "partner_energy": partner_energy,
        "match_angle": match_angle,
        "pair_mass": np.asarray(arrays["Mh_gammagamma"], dtype=float),
    }

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
        "pred_region",
        "actual_detector",
        "partner_energy",
        "match_angle",
        "pair_mass",
    ]
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
) -> Tuple[DirectedSample, DirectedSample, Dict[str, np.ndarray], Dict[str, np.ndarray]]:
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

    progress(
        f"{source_name.upper()}: denominator directed probes="
        f"{len(denominator.predicted_energy):,}; "
        f"numerator reconstructed probes="
        f"{len(numerator.predicted_energy):,}"
    )

    plot_pure_sample_inputs(
        source_name.upper(),
        epg,
        diagnostics,
        output_dir / f"01_{source_name}_tag_selection.png",
        args,
    )
    plot_pre_match_probe_diagnostics(
        source_name.upper(),
        numerator_diagnostics,
        output_dir / f"02_{source_name}_probe_prediction_before_match.png",
        args,
    )
    plot_partner_matching(
        source_name.upper(),
        numerator,
        output_dir / f"02b_{source_name}_matched_probe_diagnostics.png",
    )

    return denominator, numerator, diagnostics, epg


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

    for source_name in ["data", "aaogen", "clasdis"]:
        denominator, numerator, _, epg = process_source(
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

    plot_counts(
        curves,
        output_dir / "03_raw_probe_counts.png",
    )
    plot_efficiencies(
        curves,
        output_dir / "04_efficiency_and_data_mc_ratio.png",
    )
    plot_clasdis_tag_truth_purity(
        denominators["clasdis"],
        output_dir / "05_clasdis_tag_truth_purity.png",
        edges,
    )
    plot_threshold_scan(
        epg_arrays,
        models,
        beam_energy,
        args,
        output_dir / "06_tag_score_threshold_scan.png",
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

    print("\nOutputs:")
    print(f"  {output_dir / '03_raw_probe_counts.png'}")
    print(f"  {output_dir / '04_efficiency_and_data_mc_ratio.png'}")
    print(f"  {output_dir / '05_clasdis_tag_truth_purity.png'}")
    print(f"  {output_dir / '06_tag_score_threshold_scan.png'}")
    print("\nThis is a first-stage diagnostic extraction, not yet a final correction.")


if __name__ == "__main__":
    main()
