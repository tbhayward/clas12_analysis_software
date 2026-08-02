#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v24_fitter_preflight_and_failure_audit.py

High-energy photon efficiency extraction using an observed PASS count and a
FAIL-only BH/DVCS + exclusive-pi0 template fit.

For each selected epgamma opportunity:
  PASS = a validated native epgammagamma partner is reconstructed above 2 GeV;
  FAIL = no accepted partner is found.

The PASS yield is measured directly from data. Only the FAIL sample is fitted:

  D_FAIL = N_pi0_FAIL T_pi0_FAIL + N_BH/DVCS_FAIL T_DVCS.

The data efficiency is then

  epsilon_data = N_PASS / (N_PASS + N_pi0_FAIL),

while epsilon_MC is the direct AAOGEN PASS fraction. The fit uses one total
FAIL count and normalized shape likelihoods, avoiding the repeated extended
normalization that invalidated the earlier simultaneous fit. Independent
per-variable fits are also produced; their spread is reported as a model
diagnostic/systematic. Underflow and overflow bins are included in every fit.

The tag-dependent z, t1, and electron-tag opening-angle cuts remain disabled.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import importlib.util
import inspect
import traceback
import math
import os
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize, minimize_scalar
from scipy.special import xlogy

try:
    import uproot
except ImportError as exc:
    raise SystemExit("This script requires uproot: python3 -m pip install uproot") from exc
#endif



TREE_NAME = "PhysicsEvents"
DEFAULT_OUTPUT_DIR = "output/photon_efficiency_study"
DEFAULT_STEP_SIZE = "200 MB"
MAX_WORKERS = 8

PROTON_MASS_GEV = 0.9382720813
ELECTRON_MASS_GEV = 0.00051099895


@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy_GeV: float
    epg_data: str
    epgg_data: str
    dvcs_mc: str
    pi0_epg_mc: str
    pi0_epgg_mc: str


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        "fa18_inb", "Fa18 Inb", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "fa18_out", "Fa18 Out", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "sp18_inb", "Sp18 Inb", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp18_out", "Sp18 Out", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp19_inb", "Sp19 Inb", 10.200,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root",
    ),
)


ALIASES: Mapping[str, Tuple[str, ...]] = {
    "runnum": ("runnum", "run", "run_number"),
    "eventnum": ("evnum", "eventnum", "event", "event_number"),
    "e_p": ("e_p", "p_e", "electron_p"),
    "e_theta": ("e_theta", "theta_e", "electron_theta"),
    "e_phi": ("e_phi", "phi_e", "electron_phi"),
    "p_p": ("p1_p", "p_p", "proton_p"),
    "p_theta": ("p1_theta", "p_theta", "proton_theta"),
    "p_phi": ("p1_phi", "p_phi", "proton_phi"),
    "tag_E": ("p2_p", "g1_p", "g1_E", "gamma1_p", "photon1_p"),
    "tag_theta": ("p2_theta", "g1_theta", "gamma1_theta", "photon1_theta"),
    "tag_phi": ("p2_phi", "g1_phi", "gamma1_phi", "photon1_phi"),
    "tag_detector": ("detector2", "g1_detector", "gamma1_detector"),
    "proton_detector": ("detector1", "p_detector", "proton_detector"),
    "fiducial_status": ("fiducial_status",),
    "Q2": ("Q2", "q2"),
    "W": ("W", "w"),
    "y": ("y", "inelasticity"),
    "z": ("z",),
    "t1": ("t1", "t", "minus_t"),
    "open_angle_ep2": ("open_angle_ep2", "open_angle_eg", "open_angle_ephoton"),
    "Delta_phi": ("Delta_phi", "delta_phi", "dphi"),
    "theta_cm": ("theta", "theta2", "theta_2"),
    "theta_gamma_gamma": ("theta_gamma_gamma", "theta_pi0_pi0"),
    "pTmiss": ("pTmiss", "ptmiss", "pT_miss"),
    "Emiss2": ("Emiss2", "E_miss2", "Emiss_sq"),
    "Mx2": ("Mx2", "Mx2_epg", "Mx2_epgamma", "Mx2_eppi0", "Mx2_epi0"),
    "Mx2_1": ("Mx2_1", "Mx2_ep", "Mx2_x1", "Mx2_proton", "Mx2_p"),
    "Mx2_2": ("Mx2_2", "Mx2_egamma", "Mx2_gamma", "Mx2_pi0", "Mx2_x2"),
    # Native epgammagamma reconstruction.
    "pi0_p": ("p2_p", "pi0_p"),
    "pi0_theta": ("p2_theta", "pi0_theta"),
    "pi0_phi": ("p2_phi", "pi0_phi"),
    "pi0_mass": ("Mh_gammagamma",),
    "open1": ("open_angle_egamma1",),
    "open2": ("open_angle_egamma2",),
    "trento1": ("gamma_phi1",),
    "trento2": ("gamma_phi2",),
    "native_gamma_detector1": ("detector_gamma1",),
    "native_gamma_detector2": ("detector_gamma2",),
}


@dataclass(frozen=True)
class FitVariable:
    key: str
    label: str
    bins: int
    low: float
    high: float


FIT_VARIABLES: Tuple[FitVariable, ...] = (
    FitVariable("Delta_phi", r"$\Delta\phi$ (rad)", 100, 2.84159, 3.44159),
    FitVariable("theta_cm", r"$\theta_{p\gamma}^{\rm CM}$ (rad)", 100, 2.0, math.pi),
    FitVariable("theta_gamma_gamma", r"$\theta_{\gamma\gamma}$ (deg)", 120, 0.0, 50.0),
    FitVariable("pTmiss", r"$p_T^{\rm miss}$ (GeV)", 125, 0.0, 0.5),
    FitVariable("Emiss2", r"$E_{\rm miss}$ (GeV)", 120, 1.5, 4.5),
    FitVariable("Mx2", r"$M_x^2$ (GeV$^2$)", 100, -0.03, 0.03),
    FitVariable("Mx2_2", r"$M_{x2}^2$ (GeV$^2$)", 125, -1.0, 9.0),
)

FRACTION_DRIVERS: Tuple[str, ...] = (
    "theta_gamma_gamma",
)

TOPOLOGY_INFO: Mapping[str, Tuple[str, int, int]] = {
    "CD_FT": ("(CD, FT)", 2, 0),
    "CD_FD": ("(CD, FD)", 2, 1),
    "FD_FD": ("(FD, FD)", 1, 1),
}


@dataclass
class OpportunityRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    tag_E: np.ndarray
    tag_theta: np.ndarray
    tag_phi: np.ndarray
    tag_detector: np.ndarray
    proton_detector: np.ndarray
    probe_E: np.ndarray
    probe_p: np.ndarray
    probe_theta: np.ndarray
    probe_phi: np.ndarray
    probe_m2: np.ndarray
    probe_E_minus_p: np.ndarray
    probe_detector: np.ndarray
    probe_sector: np.ndarray
    topology_code: np.ndarray
    Delta_phi: np.ndarray
    theta_cm: np.ndarray
    theta_gamma_gamma: np.ndarray
    pTmiss: np.ndarray
    Emiss2: np.ndarray
    Mx2: np.ndarray
    Mx2_2: np.ndarray

    def size(self) -> int:
        return int(self.tag_E.size)


@dataclass
class PhotonCandidateRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    photon_E: np.ndarray
    photon_theta: np.ndarray
    photon_phi: np.ndarray
    photon_detector: np.ndarray

    def size(self) -> int:
        return int(self.photon_E.size)


@dataclass
class IdentityCatalog:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray

    def size(self) -> int:
        return int(self.e_p.size)


@dataclass
class EpggRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    g1_E: np.ndarray
    g1_theta: np.ndarray
    g1_phi: np.ndarray
    g2_E: np.ndarray
    g2_theta: np.ndarray
    g2_phi: np.ndarray
    solution_valid: np.ndarray
    solution_closure: np.ndarray
    detector1: np.ndarray
    detector2: np.ndarray

    def size(self) -> int:
        return int(self.g1_E.shape[0])


@dataclass
class PartnerMatchResult:
    matched: np.ndarray
    catalog_confirmed: np.ndarray
    partner_index: np.ndarray
    probe_angle_deg: np.ndarray
    probe_relative_E: np.ndarray
    pair_mass_GeV: np.ndarray
    score: np.ndarray
    partner_E: np.ndarray
    tag_angle_deg: np.ndarray
    tag_relative_E: np.ndarray
    solution_closure: np.ndarray
    summary: Dict[str, object]


@dataclass
class EfficiencyRow:
    period: str
    period_label: str
    detector: str
    sector: int
    observed_data: float
    expected_data: float
    expected_data_err: float
    found_data: float
    raw_efficiency_data: float
    efficiency_data: float
    efficiency_data_err: float
    expected_mc: float
    found_mc: float
    raw_efficiency_mc: float
    efficiency_mc: float
    efficiency_mc_err: float
    raw_scale_factor: float
    scale_factor: float
    scale_factor_err: float
    population_valid: bool
    population_failure_reasons: str
    fit_success: bool
    fit_quality_warning: bool
    fit_warning_reasons: str
    fit_fraction_pi0: float
    fit_fraction_pi0_err: float
    fit_deviance: float
    fit_ndf: int
    fit_reduced_deviance: float


_EXCLUSIVITY_MODULE = None


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Derive data/MC scale factors for reconstructing photons above 2 GeV."
    )
    parser.add_argument("--period", action="append", choices=[p.key for p in PERIODS])
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--exclusivity-fitter",
        default=None,
        help="Optional explicit path to plot_exclusivity_data_dvcs_pi0_mc.py.",
    )

    parser.add_argument("--tag-E-min", type=float, default=0.40)
    parser.add_argument("--probe-E-min", type=float, default=2.00)
    parser.add_argument("--probe-E-max", type=float, default=9.50)
    parser.add_argument("--probe-m2-abs-max", type=float, default=0.10)
    parser.add_argument("--probe-E-minus-p-abs-max", type=float, default=0.10)

    parser.add_argument("--ft-theta-min", type=float, default=2.5)
    parser.add_argument("--ft-theta-max", type=float, default=5.0)
    parser.add_argument("--fd-theta-min", type=float, default=5.0)
    parser.add_argument("--fd-theta-max", type=float, default=35.0)

    parser.add_argument("--Q2-min", type=float, default=1.0)
    parser.add_argument("--W-min", type=float, default=2.0)
    parser.add_argument("--y-max", type=float, default=0.8)
    parser.add_argument("--z-min", type=float, default=0.65)
    parser.add_argument("--minus-t-max", type=float, default=1.0)
    parser.add_argument("--open-angle-min-deg", type=float, default=5.0)
    parser.add_argument("--require-fiducial-status-111", action="store_true")

    parser.add_argument("--tag-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument("--tag-match-relative-E-max", type=float, default=0.35)
    parser.add_argument("--probe-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument("--probe-match-relative-E-max", type=float, default=0.35)
    parser.add_argument("--pi0-mass-min", type=float, default=0.10)
    parser.add_argument("--pi0-mass-max", type=float, default=0.17)
    parser.add_argument("--mc-signature-decimals", type=int, default=10)
    parser.add_argument("--native-closure-max", type=float, default=0.10)
    parser.add_argument("--native-minimum-photon-E", type=float, default=0.20)
    parser.add_argument("--found-probe-E-min", type=float, default=2.00)
    parser.add_argument(
        "--require-probe-residual-match",
        action="store_true",
        help="Also require the truth partner to pass probe residual windows.",
    )


    parser.add_argument("--fit-min-counts", type=int, default=100)
    parser.add_argument(
        "--fit-max-reduced-deviance",
        type=float,
        default=5.0,
        help=(
            "Issue a prominent warning if the combined detector-category "
            "template fit has deviance/ndf above this value. The fit is still "
            "retained and the efficiency is still calculated. Use a "
            "non-positive value to disable this warning threshold."
        ),
    )
    parser.add_argument(
        "--fit-fraction-boundary-margin",
        type=float,
        default=1.0e-4,
        help=(
            "Abort when the fitted pi0 fraction lies this close to 0 or 1. "
            "The default rejects boundary-dominated fits such as 0.999991."
        ),
    )
    parser.add_argument(
        "--yield-consistency-tolerance",
        type=float,
        default=1.0e-6,
        help=(
            "Absolute numerical tolerance allowed when checking that found "
            "counts do not exceed expected counts."
        ),
    )
    parser.add_argument("--exclusivity-max-shift-bins", type=float, default=5.0)
    parser.add_argument("--exclusivity-max-smear-bins", type=float, default=8.0)
    parser.add_argument("--exclusivity-shift-prior-bins", type=float, default=2.0)
    parser.add_argument("--exclusivity-smear-prior-bins", type=float, default=3.0)
    parser.add_argument("--disable-exclusivity-nuisance-penalties", action="store_true")
    parser.add_argument("--exclusivity-dvcs-core-containment", type=float, default=0.90)
    parser.add_argument("--exclusivity-dvcs-fraction-containment", type=float, default=0.95)
    parser.add_argument("--exclusivity-pi0-core-containment", type=float, default=0.90)
    parser.add_argument("--exclusivity-pi0-fraction-containment", type=float, default=0.95)
    parser.add_argument("--exclusivity-outside-overshoot-penalty", type=float, default=0.25)
    parser.add_argument("--exclusivity-emiss2-mean-order-penalty", type=float, default=25.0)
    return parser.parse_args()


def selected_periods(args: argparse.Namespace) -> List[PeriodConfig]:
    if not args.period:
        return list(PERIODS)
    # endif
    requested = set(args.period)
    return [period for period in PERIODS if period.key in requested]


def require_tree(path: str) -> Tuple[int, List[str]]:
    if not Path(path).is_file():
        raise FileNotFoundError(f"Missing input ROOT file: {path}")
    # endif
    with uproot.open(path) as root_file:
        if TREE_NAME not in root_file:
            raise KeyError(f"Tree '{TREE_NAME}' is absent from {path}")
        # endif
        tree = root_file[TREE_NAME]
        return int(tree.num_entries), [str(key) for key in tree.keys()]
    # endwith


def resolve_branches(
    path: str,
    logical_names: Sequence[str],
    optional: Sequence[str] = (),
) -> Dict[str, Optional[str]]:
    _, keys = require_tree(path)
    available = set(keys)
    optional_set = set(optional)
    resolved: Dict[str, Optional[str]] = {}
    missing: List[str] = []

    for logical in logical_names:
        aliases = ALIASES.get(logical, (logical,))
        branch = next((candidate for candidate in aliases if candidate in available), None)
        resolved[logical] = branch
        if branch is None and logical not in optional_set:
            missing.append(f"{logical}: tried {', '.join(aliases)}")
        # endif
    # endfor

    if missing:
        raise KeyError(
            f"Missing required branches in {path}:\n  " + "\n  ".join(missing)
        )
    # endif
    return resolved


def finite_array(
    arrays: Mapping[str, np.ndarray],
    branch: Optional[str],
    default: float = math.nan,
    dtype=np.float64,
) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    if branch is None:
        return np.full(n, default, dtype=dtype)
    # endif
    return np.asarray(arrays[branch], dtype=dtype)


def spherical_to_cartesian(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sin_theta = np.sin(theta)
    return (
        momentum * sin_theta * np.cos(phi),
        momentum * sin_theta * np.sin(phi),
        momentum * np.cos(theta),
    )


def massive_four_vector(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    mass: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical_to_cartesian(momentum, theta, phi)
    energy = np.sqrt(np.maximum(momentum * momentum + mass * mass, 0.0))
    return energy, px, py, pz


def photon_four_vector(
    energy: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical_to_cartesian(energy, theta, phi)
    return energy, px, py, pz


def reconstruct_probe(
    beam_energy: float,
    e_p: np.ndarray,
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    p_p: np.ndarray,
    p_theta: np.ndarray,
    p_phi: np.ndarray,
    tag_E: np.ndarray,
    tag_theta: np.ndarray,
    tag_phi: np.ndarray,
) -> Dict[str, np.ndarray]:
    e_E, e_px, e_py, e_pz = massive_four_vector(
        e_p, e_theta, e_phi, ELECTRON_MASS_GEV
    )
    p_E, p_px, p_py, p_pz = massive_four_vector(
        p_p, p_theta, p_phi, PROTON_MASS_GEV
    )
    g_E, g_px, g_py, g_pz = photon_four_vector(tag_E, tag_theta, tag_phi)

    probe_E = beam_energy + PROTON_MASS_GEV - e_E - p_E - g_E
    probe_px = -e_px - p_px - g_px
    probe_py = -e_py - p_py - g_py
    probe_pz = beam_energy - e_pz - p_pz - g_pz
    probe_p = np.sqrt(np.maximum(probe_px**2 + probe_py**2 + probe_pz**2, 0.0))
    probe_pt = np.sqrt(np.maximum(probe_px**2 + probe_py**2, 0.0))
    probe_theta = np.arctan2(probe_pt, probe_pz)
    probe_phi = np.mod(np.arctan2(probe_py, probe_px), 2.0 * math.pi)
    return {
        "E": probe_E,
        "p": probe_p,
        "theta": probe_theta,
        "phi": probe_phi,
        "m2": probe_E**2 - probe_p**2,
        "E_minus_p": probe_E - probe_p,
    }


def fd_sector_from_phi(phi_rad: np.ndarray) -> np.ndarray:
    phi_deg = np.mod(np.degrees(phi_rad), 360.0)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)
    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def classify_probe(
    theta_deg: np.ndarray,
    phi_rad: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    detector = np.full(theta_deg.shape, -1, dtype=np.int16)
    detector[
        (theta_deg >= args.ft_theta_min)
        & (theta_deg < args.ft_theta_max)
    ] = 0
    detector[
        (theta_deg >= args.fd_theta_min)
        & (theta_deg < args.fd_theta_max)
    ] = 1
    sector = fd_sector_from_phi(phi_rad)
    sector[detector != 1] = 0
    return detector, sector


def topology_code(
    proton_detector: np.ndarray,
    tag_detector: np.ndarray,
) -> np.ndarray:
    code = np.full(proton_detector.shape, -1, dtype=np.int8)
    code[(proton_detector == 2) & (tag_detector == 0)] = 0  # CD_FT
    code[(proton_detector == 2) & (tag_detector == 1)] = 1  # CD_FD
    code[(proton_detector == 1) & (tag_detector == 1)] = 2  # FD_FD
    return code


def topology_key(code: int) -> str:
    return {0: "CD_FT", 1: "CD_FD", 2: "FD_FD"}[int(code)]


def opening_angle_deg(
    theta1: np.ndarray,
    phi1: np.ndarray,
    theta2: np.ndarray,
    phi2: np.ndarray,
) -> np.ndarray:
    cosine = (
        np.sin(theta1) * np.sin(theta2) * np.cos(phi1 - phi2)
        + np.cos(theta1) * np.cos(theta2)
    )
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def concat(parts: Sequence[np.ndarray], dtype=np.float64) -> np.ndarray:
    if not parts:
        return np.asarray([], dtype=dtype)
    # endif
    return np.concatenate(parts).astype(dtype, copy=False)


def subset_opportunities(records: OpportunityRecords, mask: np.ndarray) -> OpportunityRecords:
    return OpportunityRecords(
        **{
            name: np.asarray(getattr(records, name))[mask]
            for name in records.__dataclass_fields__
        }
    )


def exact_duplicate_keep_mask(
    arrays: Sequence[np.ndarray],
    decimals: int = 10,
) -> np.ndarray:
    if not arrays or len(arrays[0]) == 0:
        return np.ones(0, dtype=bool)
    # endif
    rounded = [
        np.round(np.asarray(values, dtype=np.float64), decimals)
        for values in arrays
    ]
    matrix = np.column_stack(rounded)
    _, first = np.unique(matrix, axis=0, return_index=True)
    keep = np.zeros(matrix.shape[0], dtype=bool)
    keep[np.sort(first)] = True
    return keep



CUT_STAGE_ORDER: Tuple[str, ...] = (
    "tree_entries",
    "finite_kinematics",
    "Q2",
    "W",
    "y",
    "z_tag",
    "minus_t1_tag",
    "open_angle_ep2_tag",
    "fiducial_status",
    "tag_energy",
    "probe_finite",
    "probe_energy_min",
    "probe_energy_max",
    "probe_m2",
    "probe_E_minus_p",
    "probe_acceptance",
    "supported_topology",
    "selected_before_deduplication",
    "selected_after_deduplication",
)


def cut_stage_definitions(
    resolved: Mapping[str, Optional[str]],
    args: argparse.Namespace,
) -> Dict[str, Dict[str, object]]:
    """Return human-readable definitions for every sequential opportunity cut."""
    return {
        "tree_entries": {
            "label": "Input ROOT-tree entries",
            "condition": "all entries",
            "branch": None,
            "applied": True,
            "scope": "input",
        },
        "finite_kinematics": {
            "label": "Finite positive reconstructed e, p, and tag photon",
            "condition": (
                "finite e/p/gamma momentum and angles; "
                "p_e > 0, p_p > 0, E_tag > 0"
            ),
            "branch": "e_p,e_theta,e_phi,p1_p,p1_theta,p1_phi,p2_p,p2_theta,p2_phi",
            "applied": True,
            "scope": "common event/tag",
        },
        "Q2": {
            "label": "DIS virtuality",
            "condition": f"Q2 > {args.Q2_min:g} GeV^2",
            "branch": resolved.get("Q2"),
            "applied": resolved.get("Q2") is not None,
            "scope": "common event",
        },
        "W": {
            "label": "DIS invariant mass",
            "condition": f"W > {args.W_min:g} GeV",
            "branch": resolved.get("W"),
            "applied": resolved.get("W") is not None,
            "scope": "common event",
        },
        "y": {
            "label": "DIS inelasticity",
            "condition": f"y < {args.y_max:g}",
            "branch": resolved.get("y"),
            "applied": resolved.get("y") is not None,
            "scope": "common event",
        },
        "z_tag": {
            "label": "Tag-photon z requirement",
            "condition": (
                f"NOT APPLIED: z(tag photon) > {args.z_min:g} "
                "was removed from the nominal tag-and-probe selection"
            ),
            "branch": resolved.get("z"),
            "applied": False,
            "scope": "REMOVED — tag-dependent DVCS cut",
        },
        "minus_t1_tag": {
            "label": "Tag-photon -t1 requirement",
            "condition": (
                f"NOT APPLIED: -t1(tag photon) < {args.minus_t_max:g} GeV^2 "
                "was removed from the nominal tag-and-probe selection"
            ),
            "branch": resolved.get("t1"),
            "applied": False,
            "scope": "REMOVED — tag-dependent DVCS cut",
        },
        "open_angle_ep2_tag": {
            "label": "Electron–tag-photon opening angle",
            "condition": (
                f"NOT APPLIED: open_angle(e, tag photon) > "
                f"{args.open_angle_min_deg:g} deg was removed from the "
                "nominal tag-and-probe selection"
            ),
            "branch": resolved.get("open_angle_ep2"),
            "applied": False,
            "scope": "REMOVED — tag-dependent DVCS cut",
        },
        "fiducial_status": {
            "label": "Combined fiducial status",
            "condition": "fiducial_status == 111",
            "branch": resolved.get("fiducial_status"),
            "applied": (
                bool(args.require_fiducial_status_111)
                and resolved.get("fiducial_status") is not None
            ),
            "scope": "optional reconstructed-particle fiducial",
        },
        "tag_energy": {
            "label": "Reconstructed tag-photon energy",
            "condition": f"E_tag >= {args.tag_E_min:g} GeV",
            "branch": resolved.get("tag_E"),
            "applied": True,
            "scope": "tag",
        },
        "probe_finite": {
            "label": "Finite physical missing-photon four-vector",
            "condition": (
                "finite E_X, |p_X|, theta_X, phi_X, m_X^2, E_X-|p_X|; "
                "E_X > 0 and |p_X| > 0"
            ),
            "branch": "derived from beam + target - e - p - tag",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_energy_min": {
            "label": "Predicted probe minimum energy",
            "condition": f"E_X >= {args.probe_E_min:g} GeV",
            "branch": "derived E_X",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_energy_max": {
            "label": "Predicted probe maximum energy",
            "condition": f"E_X < {args.probe_E_max:g} GeV",
            "branch": "derived E_X",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_m2": {
            "label": "Predicted probe mass-squared consistency",
            "condition": f"|m_X^2| < {args.probe_m2_abs_max:g} GeV^2",
            "branch": "derived m_X^2",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_E_minus_p": {
            "label": "Predicted probe massless-energy consistency",
            "condition": (
                f"|E_X - |p_X|| < "
                f"{args.probe_E_minus_p_abs_max:g} GeV"
            ),
            "branch": "derived E_X-|p_X|",
            "applied": True,
            "scope": "predicted probe",
        },
        "probe_acceptance": {
            "label": "Predicted probe points into FT or FD",
            "condition": (
                f"FT: {args.ft_theta_min:g} <= theta_X < "
                f"{args.ft_theta_max:g} deg; FD: "
                f"{args.fd_theta_min:g} <= theta_X < "
                f"{args.fd_theta_max:g} deg"
            ),
            "branch": "derived theta_X and phi_X",
            "applied": True,
            "scope": "predicted probe detector",
        },
        "supported_topology": {
            "label": "Supported reconstructed proton/tag topology",
            "condition": "proton and tag detector codes map to a supported topology",
            "branch": "detector1, detector2",
            "applied": True,
            "scope": "bookkeeping only; proton topologies later combined",
        },
        "selected_before_deduplication": {
            "label": "Selected directed opportunities before MC deduplication",
            "condition": "all preceding requirements",
            "branch": None,
            "applied": True,
            "scope": "selected sample",
        },
        "selected_after_deduplication": {
            "label": "Selected directed opportunities after MC deduplication",
            "condition": "exact duplicate removal for MC; unchanged for data",
            "branch": None,
            "applied": True,
            "scope": "final fit denominator population",
        },
    }


def stage_mask_with_optional_branch(
    previous: np.ndarray,
    arrays: Mapping[str, np.ndarray],
    branch: Optional[str],
    predicate,
) -> np.ndarray:
    """Apply one branch cut sequentially; absent optional branches leave the mask unchanged."""
    if branch is None:
        return previous.copy()
    # endif
    values = finite_array(arrays, branch)
    return previous & np.isfinite(values) & predicate(values)


def finalize_cutflow_report(cutflow: Dict[str, object]) -> Dict[str, object]:
    """Add incremental and total survival rates to the raw cumulative counters."""
    definitions = cutflow["stage_definitions"]
    counts = cutflow["stage_counts"]
    rows: List[Dict[str, object]] = []
    total = int(counts.get("tree_entries", 0))
    previous_count = total

    for index, stage in enumerate(CUT_STAGE_ORDER):
        count = int(counts.get(stage, 0))
        definition = definitions[stage]
        applied = bool(definition["applied"])
        reference = total if index == 0 else previous_count
        incremental = (
            float(count / reference)
            if reference > 0
            else math.nan
        )
        total_survival = (
            float(count / total)
            if total > 0
            else math.nan
        )
        rows.append(
            {
                "order": index,
                "stage": stage,
                "label": definition["label"],
                "condition": definition["condition"],
                "branch": definition["branch"],
                "scope": definition["scope"],
                "applied": applied,
                "count": count,
                "rejected_at_stage": (
                    0 if index == 0 else max(previous_count - count, 0)
                ),
                "incremental_survival": incremental,
                "incremental_rejection": (
                    1.0 - incremental
                    if math.isfinite(incremental)
                    else math.nan
                ),
                "total_survival": total_survival,
            }
        )
        previous_count = count
    # endfor

    cutflow["stages"] = rows
    cutflow["largest_incremental_loss"] = max(
        rows[1:],
        key=lambda row: (
            row["incremental_rejection"]
            if math.isfinite(row["incremental_rejection"])
            else -1.0
        ),
    ) if len(rows) > 1 else None
    return cutflow


def common_global_mask(
    arrays: Mapping[str, np.ndarray],
    resolved: Mapping[str, Optional[str]],
    args: argparse.Namespace,
) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    mask = np.ones(n, dtype=bool)

    def lower(key: str, threshold: float) -> None:
        nonlocal mask
        branch = resolved.get(key)
        if branch is not None:
            values = finite_array(arrays, branch)
            mask &= np.isfinite(values) & (values > threshold)
        # endif

    def upper(key: str, threshold: float) -> None:
        nonlocal mask
        branch = resolved.get(key)
        if branch is not None:
            values = finite_array(arrays, branch)
            mask &= np.isfinite(values) & (values < threshold)
        # endif

    lower("Q2", args.Q2_min)
    lower("W", args.W_min)
    upper("y", args.y_max)
    lower("z", args.z_min)

    # Tag-dependent z, t1, and electron-tag opening-angle cuts are
    # intentionally omitted from the nominal opportunity definition.

    if args.require_fiducial_status_111 and resolved.get("fiducial_status") is not None:
        status = finite_array(arrays, resolved["fiducial_status"])
        mask &= np.isfinite(status) & (status == 111)
    # endif
    return mask


def empty_opportunity_store() -> Dict[str, List[np.ndarray]]:
    return {name: [] for name in OpportunityRecords.__dataclass_fields__}


def append_opportunity_store(
    store: Dict[str, List[np.ndarray]],
    values: Mapping[str, np.ndarray],
    mask: np.ndarray,
) -> None:
    for name in store:
        store[name].append(np.asarray(values[name])[mask])
    # endfor


def finalize_opportunity_store(store: Mapping[str, Sequence[np.ndarray]]) -> OpportunityRecords:
    integer_fields = {
        "runnum", "eventnum", "tag_detector", "proton_detector",
        "probe_detector", "probe_sector", "topology_code",
    }
    payload = {}
    for name, parts in store.items():
        dtype = np.int64 if name in {"runnum", "eventnum"} else (
            np.int16 if name in integer_fields else np.float64
        )
        payload[name] = concat(parts, dtype=dtype)
    # endfor
    return OpportunityRecords(**payload)


def empty_candidate_store() -> Dict[str, List[np.ndarray]]:
    return {name: [] for name in PhotonCandidateRecords.__dataclass_fields__}


def append_candidate_store(
    store: Dict[str, List[np.ndarray]],
    values: Mapping[str, np.ndarray],
    mask: np.ndarray,
) -> None:
    for name in store:
        store[name].append(np.asarray(values[name])[mask])
    # endfor


def finalize_candidate_store(
    store: Mapping[str, Sequence[np.ndarray]],
) -> PhotonCandidateRecords:
    integer_fields = {"runnum", "eventnum", "photon_detector"}
    payload = {}
    for name, parts in store.items():
        dtype = np.int64 if name in {"runnum", "eventnum"} else (
            np.int16 if name in integer_fields else np.float64
        )
        payload[name] = concat(parts, dtype=dtype)
    # endfor
    return PhotonCandidateRecords(**payload)


def subset_candidates(
    records: PhotonCandidateRecords,
    mask: np.ndarray,
) -> PhotonCandidateRecords:
    return PhotonCandidateRecords(
        **{
            name: np.asarray(getattr(records, name))[mask]
            for name in records.__dataclass_fields__
        }
    )


def read_opportunities(
    path: str,
    beam_energy: float,
    role: str,
    args: argparse.Namespace,
    deduplicate: bool,
) -> Tuple[OpportunityRecords, PhotonCandidateRecords, Dict[str, object]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "tag_E", "tag_theta", "tag_phi",
        "tag_detector", "proton_detector",
        "fiducial_status", "Q2", "W", "y", "z", "t1", "open_angle_ep2",
        "Delta_phi", "theta_cm", "theta_gamma_gamma", "pTmiss",
        "Emiss2", "Mx2", "Mx2_1", "Mx2_2",
    )
    optional = (
        "runnum", "eventnum", "tag_detector", "proton_detector",
        "fiducial_status", "Q2", "W", "y", "z", "t1", "open_angle_ep2",
        "theta_cm", "Mx2_1",
    )
    resolved = resolve_branches(path, logical, optional=optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    store = empty_opportunity_store()
    candidate_store = empty_candidate_store()
    stage_definitions = cut_stage_definitions(resolved, args)
    cutflow: Dict[str, object] = {
        "role": role,
        "path": path,
        "stage_definitions": stage_definitions,
        "stage_counts": {stage: 0 for stage in CUT_STAGE_ORDER},
        "detector_counts": {
            stage: {"FT": 0, "FD": 0}
            for stage in (
                "probe_finite",
                "probe_energy_min",
                "probe_energy_max",
                "probe_m2",
                "probe_E_minus_p",
                "probe_acceptance",
                "supported_topology",
                "selected_before_deduplication",
                "selected_after_deduplication",
            )
        },
        "partner_candidates_finite_and_E_above_threshold": 0,
    }
    stage_counts = cutflow["stage_counts"]

    seen = 0
    log(f"Reading {role} from {path}")
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        # endif
        if args.max_events is not None and seen + n > args.max_events:
            n = args.max_events - seen
            arrays = {key: value[:n] for key, value in arrays.items()}
        # endif
        seen += n
        stage_counts["tree_entries"] += int(n)

        def arr(key: str, default: float = math.nan, dtype=np.float64) -> np.ndarray:
            return finite_array(arrays, resolved.get(key), default=default, dtype=dtype)

        e_p, e_theta, e_phi = arr("e_p"), arr("e_theta"), arr("e_phi")
        p_p, p_theta, p_phi = arr("p_p"), arr("p_theta"), arr("p_phi")
        tag_E = arr("tag_E")
        tag_theta = arr("tag_theta")
        tag_phi = np.mod(arr("tag_phi"), 2.0 * math.pi)
        tag_detector = arr("tag_detector", default=-1, dtype=np.int16)
        proton_detector = arr("proton_detector", default=-1, dtype=np.int16)

        finite = (
            np.isfinite(e_p) & np.isfinite(e_theta) & np.isfinite(e_phi)
            & np.isfinite(p_p) & np.isfinite(p_theta) & np.isfinite(p_phi)
            & np.isfinite(tag_E) & np.isfinite(tag_theta) & np.isfinite(tag_phi)
            & (e_p > 0.0) & (p_p > 0.0) & (tag_E > 0.0)
        )
        stage_counts["finite_kinematics"] += int(np.count_nonzero(finite))

        mask_Q2 = stage_mask_with_optional_branch(
            finite, arrays, resolved.get("Q2"),
            lambda values: values > args.Q2_min,
        )
        stage_counts["Q2"] += int(np.count_nonzero(mask_Q2))

        mask_W = stage_mask_with_optional_branch(
            mask_Q2, arrays, resolved.get("W"),
            lambda values: values > args.W_min,
        )
        stage_counts["W"] += int(np.count_nonzero(mask_W))

        mask_y = stage_mask_with_optional_branch(
            mask_W, arrays, resolved.get("y"),
            lambda values: values < args.y_max,
        )
        stage_counts["y"] += int(np.count_nonzero(mask_y))

        # These three DVCS-style cuts are tag dependent. The reconstructed
        # tag may be the low-energy pi0 daughter, whereas the efficiency probe
        # is the predicted high-energy photon. They are therefore retained as
        # explicit audit stages but intentionally do not alter the mask.
        mask_z = mask_y.copy()
        stage_counts["z_tag"] += int(np.count_nonzero(mask_z))

        mask_t1 = mask_z.copy()
        stage_counts["minus_t1_tag"] += int(np.count_nonzero(mask_t1))

        mask_open = mask_t1.copy()
        stage_counts["open_angle_ep2_tag"] += int(np.count_nonzero(mask_open))

        if (
            args.require_fiducial_status_111
            and resolved.get("fiducial_status") is not None
        ):
            status = finite_array(arrays, resolved["fiducial_status"])
            mask_fiducial = (
                mask_open
                & np.isfinite(status)
                & (status == 111)
            )
        else:
            mask_fiducial = mask_open.copy()
        # endif
        stage_counts["fiducial_status"] += int(
            np.count_nonzero(mask_fiducial)
        )

        tag_mask = mask_fiducial & (tag_E >= args.tag_E_min)
        stage_counts["tag_energy"] += int(np.count_nonzero(tag_mask))

        # The partner pool intentionally uses every finite reconstructed photon
        # above threshold, independent of the tag-specific opportunity cuts.
        candidate_mask = finite & (tag_E >= args.tag_E_min)
        candidate_values = {
            "runnum": arr("runnum", default=-1, dtype=np.int64),
            "eventnum": arr("eventnum", default=-1, dtype=np.int64),
            "e_p": e_p, "e_theta": e_theta, "e_phi": e_phi,
            "p_p": p_p, "p_theta": p_theta, "p_phi": p_phi,
            "photon_E": tag_E,
            "photon_theta": tag_theta,
            "photon_phi": tag_phi,
            "photon_detector": tag_detector,
        }
        append_candidate_store(candidate_store, candidate_values, candidate_mask)
        cutflow["partner_candidates_finite_and_E_above_threshold"] += int(
            np.count_nonzero(candidate_mask)
        )

        probe = reconstruct_probe(
            beam_energy,
            e_p, e_theta, e_phi,
            p_p, p_theta, p_phi,
            tag_E, tag_theta, tag_phi,
        )
        probe_theta_deg = np.degrees(probe["theta"])
        probe_detector, probe_sector = classify_probe(
            probe_theta_deg, probe["phi"], args
        )

        def accumulate_detector_stage(stage: str, mask: np.ndarray) -> None:
            cutflow["detector_counts"][stage]["FT"] += int(
                np.count_nonzero(mask & (probe_detector == 0))
            )
            cutflow["detector_counts"][stage]["FD"] += int(
                np.count_nonzero(mask & (probe_detector == 1))
            )
        # enddef

        probe_finite = (
            tag_mask
            & np.isfinite(probe["E"]) & np.isfinite(probe["p"])
            & np.isfinite(probe["theta"]) & np.isfinite(probe["phi"])
            & np.isfinite(probe["m2"]) & np.isfinite(probe["E_minus_p"])
            & (probe["E"] > 0.0) & (probe["p"] > 0.0)
        )
        stage_counts["probe_finite"] += int(np.count_nonzero(probe_finite))
        accumulate_detector_stage("probe_finite", probe_finite)

        probe_E_min = probe_finite & (probe["E"] >= args.probe_E_min)
        stage_counts["probe_energy_min"] += int(
            np.count_nonzero(probe_E_min)
        )
        accumulate_detector_stage("probe_energy_min", probe_E_min)

        probe_E_window = probe_E_min & (probe["E"] < args.probe_E_max)
        stage_counts["probe_energy_max"] += int(
            np.count_nonzero(probe_E_window)
        )
        accumulate_detector_stage("probe_energy_max", probe_E_window)

        probe_m2 = (
            probe_E_window
            & (np.abs(probe["m2"]) < args.probe_m2_abs_max)
        )
        stage_counts["probe_m2"] += int(np.count_nonzero(probe_m2))
        accumulate_detector_stage("probe_m2", probe_m2)

        photon_like = (
            probe_m2
            & (
                np.abs(probe["E_minus_p"])
                < args.probe_E_minus_p_abs_max
            )
        )
        stage_counts["probe_E_minus_p"] += int(
            np.count_nonzero(photon_like)
        )
        accumulate_detector_stage("probe_E_minus_p", photon_like)

        accepted = photon_like & (probe_detector >= 0)
        stage_counts["probe_acceptance"] += int(np.count_nonzero(accepted))
        accumulate_detector_stage("probe_acceptance", accepted)

        topo = topology_code(proton_detector, tag_detector)
        selected = accepted & (topo >= 0)
        stage_counts["supported_topology"] += int(np.count_nonzero(selected))
        stage_counts["selected_before_deduplication"] += int(
            np.count_nonzero(selected)
        )
        accumulate_detector_stage("supported_topology", selected)
        accumulate_detector_stage("selected_before_deduplication", selected)

        values = {
            "runnum": arr("runnum", default=-1, dtype=np.int64),
            "eventnum": arr("eventnum", default=-1, dtype=np.int64),
            "e_p": e_p, "e_theta": e_theta, "e_phi": e_phi,
            "p_p": p_p, "p_theta": p_theta, "p_phi": p_phi,
            "tag_E": tag_E, "tag_theta": tag_theta, "tag_phi": tag_phi,
            "tag_detector": tag_detector,
            "proton_detector": proton_detector,
            "probe_E": probe["E"], "probe_p": probe["p"],
            "probe_theta": probe["theta"], "probe_phi": probe["phi"],
            "probe_m2": probe["m2"],
            "probe_E_minus_p": probe["E_minus_p"],
            "probe_detector": probe_detector,
            "probe_sector": probe_sector,
            "topology_code": topo,
            "Delta_phi": arr("Delta_phi"),
            "theta_cm": arr("theta_cm"),
            "theta_gamma_gamma": arr("theta_gamma_gamma"),
            "pTmiss": arr("pTmiss"),
            "Emiss2": arr("Emiss2"),
            "Mx2": arr("Mx2"),
            "Mx2_2": arr("Mx2_2"),
        }
        append_opportunity_store(store, values, selected)
    # endfor

    records = finalize_opportunity_store(store)
    candidates = finalize_candidate_store(candidate_store)
    if deduplicate and candidates.size() > 0:
        candidate_keep = exact_duplicate_keep_mask(
            [
                candidates.e_p, candidates.e_theta, candidates.e_phi,
                candidates.p_p, candidates.p_theta, candidates.p_phi,
                candidates.photon_E, candidates.photon_theta,
                candidates.photon_phi,
            ],
            decimals=args.mc_signature_decimals,
        )
        candidates = subset_candidates(candidates, candidate_keep)
    # endif
    if deduplicate and records.size() > 0:
        keep = exact_duplicate_keep_mask(
            [
                records.e_p, records.e_theta, records.e_phi,
                records.p_p, records.p_theta, records.p_phi,
                records.tag_E, records.tag_theta, records.tag_phi,
            ],
            decimals=args.mc_signature_decimals,
        )
        records = subset_opportunities(records, keep)
    # endif
    stage_counts["selected_after_deduplication"] = records.size()
    cutflow["detector_counts"]["selected_after_deduplication"]["FT"] = int(
        np.count_nonzero(records.probe_detector == 0)
    )
    cutflow["detector_counts"]["selected_after_deduplication"]["FD"] = int(
        np.count_nonzero(records.probe_detector == 1)
    )
    cutflow["branch_mapping"] = resolved
    cutflow["partner_candidate_records"] = candidates.size()
    cutflow = finalize_cutflow_report(cutflow)
    return records, candidates, cutflow


def unit_vector(theta: float, phi: float) -> np.ndarray:
    return np.asarray(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ],
        dtype=float,
    )

def q_trento_basis(
    beam_energy: float,
    electron_p: float,
    electron_theta: float,
    electron_phi: float,
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    electron_direction = unit_vector(electron_theta, electron_phi)
    electron_vector = electron_p * electron_direction
    q_vector = np.asarray([0.0, 0.0, beam_energy]) - electron_vector
    q_magnitude = float(np.linalg.norm(q_vector))
    if q_magnitude <= 1.0e-12:
        return None
    # endif
    q_hat = q_vector / q_magnitude

    electron_transverse = (
        electron_direction
        - float(np.dot(electron_direction, q_hat)) * q_hat
    )
    transverse_magnitude = float(np.linalg.norm(electron_transverse))
    if transverse_magnitude <= 1.0e-12:
        return None
    # endif

    x_hat = electron_transverse / transverse_magnitude
    y_hat = np.cross(q_hat, x_hat)
    y_hat /= float(np.linalg.norm(y_hat))
    return q_hat, x_hat, y_hat, electron_direction

def solve_q_polar_angles(
    electron_direction: np.ndarray,
    q_hat: np.ndarray,
    x_hat: np.ndarray,
    trento_phi: float,
    opening_angle: float,
) -> List[float]:
    """
    Solve for the photon polar angle about q.

    With n_gamma = sin(theta_q)(cos(phi)x + sin(phi)y) + cos(theta_q)q,
    the stored electron-photon opening angle gives

        A sin(theta_q) + B cos(theta_q) = cos(alpha),

    where A = (e_hat dot x_hat) cos(phi) and B = e_hat dot q_hat.
    """
    A = float(np.dot(electron_direction, x_hat)) * math.cos(trento_phi)
    B = float(np.dot(electron_direction, q_hat))
    C = math.cos(opening_angle)
    radius = math.hypot(A, B)
    if radius <= 1.0e-12:
        return []
    # endif

    ratio = C / radius
    if ratio < -1.0 - 1.0e-9 or ratio > 1.0 + 1.0e-9:
        return []
    # endif
    ratio = float(np.clip(ratio, -1.0, 1.0))

    phase = math.atan2(A, B)
    offset = math.acos(ratio)
    candidates: List[float] = []
    for raw in (phase + offset, phase - offset):
        theta = raw % (2.0 * math.pi)
        if theta > math.pi:
            theta = 2.0 * math.pi - theta
        # endif
        if (
            0.0 <= theta <= math.pi
            and all(abs(theta - existing) > 1.0e-9 for existing in candidates)
        ):
            candidates.append(theta)
        # endif
    # endfor
    return candidates

def lab_direction_from_trento(
    theta_q: float,
    trento_phi: float,
    q_hat: np.ndarray,
    x_hat: np.ndarray,
    y_hat: np.ndarray,
) -> np.ndarray:
    direction = (
        math.sin(theta_q)
        * (
            math.cos(trento_phi) * x_hat
            + math.sin(trento_phi) * y_hat
        )
        + math.cos(theta_q) * q_hat
    )
    direction /= float(np.linalg.norm(direction))
    return direction

def direction_to_angles(direction: np.ndarray) -> Tuple[float, float]:
    theta = math.acos(float(np.clip(direction[2], -1.0, 1.0)))
    phi = math.atan2(float(direction[1]), float(direction[0]))
    return theta, float(np.mod(phi, 2.0 * math.pi))

def reconstruct_daughter_solutions(
    beam_energy: float,
    electron_p: float,
    electron_theta: float,
    electron_phi: float,
    pi0_p: float,
    pi0_theta: float,
    pi0_phi: float,
    pi0_mass: float,
    opening1: float,
    opening2: float,
    trento1: float,
    trento2: float,
    args: argparse.Namespace,
) -> List[Tuple[float, float, float, float, float, float, float]]:
    basis = q_trento_basis(
        beam_energy,
        electron_p,
        electron_theta,
        electron_phi,
    )
    if basis is None:
        return []
    # endif
    q_hat, x_hat, y_hat, electron_direction = basis

    theta1_candidates = solve_q_polar_angles(
        electron_direction,
        q_hat,
        x_hat,
        trento1,
        opening1,
    )
    theta2_candidates = solve_q_polar_angles(
        electron_direction,
        q_hat,
        x_hat,
        trento2,
        opening2,
    )
    if not theta1_candidates or not theta2_candidates:
        return []
    # endif

    pi0_direction = unit_vector(pi0_theta, pi0_phi)
    pi0_vector = pi0_p * pi0_direction
    pi0_energy = math.sqrt(max(pi0_p**2 + pi0_mass**2, 0.0))
    target = np.asarray([pi0_energy, *pi0_vector], dtype=float)
    normalization = max(pi0_energy, 1.0e-12)

    solutions = []
    for theta1_q in theta1_candidates:
        direction1 = lab_direction_from_trento(
            theta1_q, trento1, q_hat, x_hat, y_hat
        )
        theta1_lab, phi1_lab = direction_to_angles(direction1)

        for theta2_q in theta2_candidates:
            direction2 = lab_direction_from_trento(
                theta2_q, trento2, q_hat, x_hat, y_hat
            )
            theta2_lab, phi2_lab = direction_to_angles(direction2)

            response = np.asarray(
                [
                    [1.0, 1.0],
                    [direction1[0], direction2[0]],
                    [direction1[1], direction2[1]],
                    [direction1[2], direction2[2]],
                ],
                dtype=float,
            )
            energies, _, rank, _ = np.linalg.lstsq(
                response,
                target,
                rcond=None,
            )
            if (
                rank < 2
                or energies[0] <= args.native_minimum_photon_E
                or energies[1] <= args.native_minimum_photon_E
            ):
                continue
            # endif

            closure = float(
                np.linalg.norm(response @ energies - target)
                / normalization
            )
            if closure > args.native_closure_max:
                continue
            # endif

            solutions.append(
                (
                    closure,
                    float(energies[0]),
                    theta1_lab,
                    phi1_lab,
                    float(energies[1]),
                    theta2_lab,
                    phi2_lab,
                )
            )
        # endfor
    # endfor

    solutions.sort(key=lambda item: item[0])
    return solutions[:4]

def read_epgg(
    path: str,
    beam_energy: float,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[EpggRecords, Dict[str, int]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "pi0_p", "pi0_theta", "pi0_phi", "pi0_mass",
        "open1", "open2", "trento1", "trento2",
        "native_gamma_detector1", "native_gamma_detector2",
    )
    optional_identity = (
        ("runnum", "eventnum")
        if "AAOGEN" in mode.upper()
        else ()
    )
    resolved = resolve_branches(path, logical, optional=optional_identity)
    expressions = sorted(
        {branch for branch in resolved.values() if branch is not None}
    )

    scalar_store: Dict[str, List[np.ndarray]] = {
        name: [] for name in (
            "runnum", "eventnum",
            "e_p", "e_theta", "e_phi",
            "p_p", "p_theta", "p_phi",
            "detector1", "detector2",
        )
    }
    solution_store: Dict[str, List[np.ndarray]] = {
        name: [] for name in (
            "g1_E", "g1_theta", "g1_phi",
            "g2_E", "g2_theta", "g2_phi",
            "solution_valid", "solution_closure",
        )
    }
    counts = {
        "tree_entries": 0,
        "events_with_at_least_one_solution": 0,
        "valid_solution_count": 0,
    }

    seen = 0
    log(f"Reading {mode} epgammagamma records from {path}")
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        # endif
        if args.max_events is not None and seen + n > args.max_events:
            n = args.max_events - seen
            arrays = {key: values[:n] for key, values in arrays.items()}
        # endif
        seen += n
        counts["tree_entries"] += int(n)

        def arr(name: str, default=math.nan, dtype=np.float64):
            return finite_array(arrays, resolved.get(name), default, dtype)

        e_p = arr("e_p")
        e_theta = arr("e_theta")
        e_phi = np.mod(arr("e_phi"), 2.0 * math.pi)
        p_p = arr("p_p")
        p_theta = arr("p_theta")
        p_phi = np.mod(arr("p_phi"), 2.0 * math.pi)
        pi0_p = arr("pi0_p")
        pi0_theta = arr("pi0_theta")
        pi0_phi = np.mod(arr("pi0_phi"), 2.0 * math.pi)
        pi0_mass = arr("pi0_mass")
        open1 = np.deg2rad(arr("open1"))
        open2 = np.deg2rad(arr("open2"))
        trento1 = np.mod(arr("trento1"), 2.0 * math.pi)
        trento2 = np.mod(arr("trento2"), 2.0 * math.pi)

        max_solutions = 4
        g1_E = np.full((n, max_solutions), np.nan)
        g1_theta = np.full((n, max_solutions), np.nan)
        g1_phi = np.full((n, max_solutions), np.nan)
        g2_E = np.full((n, max_solutions), np.nan)
        g2_theta = np.full((n, max_solutions), np.nan)
        g2_phi = np.full((n, max_solutions), np.nan)
        valid = np.zeros((n, max_solutions), dtype=bool)
        closure = np.full((n, max_solutions), np.nan)

        for index in range(n):
            values = (
                e_p[index], e_theta[index], e_phi[index],
                pi0_p[index], pi0_theta[index], pi0_phi[index],
                pi0_mass[index], open1[index], open2[index],
                trento1[index], trento2[index],
            )
            if not all(math.isfinite(float(value)) for value in values):
                continue
            # endif

            solutions = reconstruct_daughter_solutions(
                beam_energy,
                float(e_p[index]),
                float(e_theta[index]),
                float(e_phi[index]),
                float(pi0_p[index]),
                float(pi0_theta[index]),
                float(pi0_phi[index]),
                float(pi0_mass[index]),
                float(open1[index]),
                float(open2[index]),
                float(trento1[index]),
                float(trento2[index]),
                args,
            )
            for solution_index, solution in enumerate(solutions):
                (
                    closure_value,
                    energy1, theta1, phi1,
                    energy2, theta2, phi2,
                ) = solution
                closure[index, solution_index] = closure_value
                g1_E[index, solution_index] = energy1
                g1_theta[index, solution_index] = theta1
                g1_phi[index, solution_index] = phi1
                g2_E[index, solution_index] = energy2
                g2_theta[index, solution_index] = theta2
                g2_phi[index, solution_index] = phi2
                valid[index, solution_index] = True
            # endfor
        # endfor

        counts["events_with_at_least_one_solution"] += int(
            np.count_nonzero(np.any(valid, axis=1))
        )
        counts["valid_solution_count"] += int(np.count_nonzero(valid))

        scalar_values = {
            "runnum": arr("runnum", -1, np.int64),
            "eventnum": arr("eventnum", -1, np.int64),
            "e_p": e_p,
            "e_theta": e_theta,
            "e_phi": e_phi,
            "p_p": p_p,
            "p_theta": p_theta,
            "p_phi": p_phi,
            "detector1": arr("native_gamma_detector1", -1, np.int16),
            "detector2": arr("native_gamma_detector2", -1, np.int16),
        }
        for name in scalar_store:
            scalar_store[name].append(np.asarray(scalar_values[name]))
        # endfor

        solution_values = {
            "g1_E": g1_E,
            "g1_theta": g1_theta,
            "g1_phi": g1_phi,
            "g2_E": g2_E,
            "g2_theta": g2_theta,
            "g2_phi": g2_phi,
            "solution_valid": valid,
            "solution_closure": closure,
        }
        for name in solution_store:
            solution_store[name].append(solution_values[name])
        # endfor
    # endfor

    payload = {}
    for name, parts in scalar_store.items():
        if parts:
            payload[name] = np.concatenate(parts)
        else:
            dtype = (
                np.int64
                if name in {"runnum", "eventnum"}
                else np.int16
                if name in {"detector1", "detector2"}
                else np.float64
            )
            payload[name] = np.asarray([], dtype=dtype)
        # endif
    # endfor
    for name, parts in solution_store.items():
        if parts:
            payload[name] = np.concatenate(parts, axis=0)
        else:
            dtype = bool if name == "solution_valid" else np.float64
            payload[name] = np.empty((0, 4), dtype=dtype)
        # endif
    # endfor

    return EpggRecords(**payload), counts

def read_identity_catalog(
    path: str,
    role: str,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[IdentityCatalog, Dict[str, object]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
    )
    optional = ("runnum", "eventnum") if mode == "mc" else ()
    resolved = resolve_branches(path, logical, optional=optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    store = {name: [] for name in IdentityCatalog.__dataclass_fields__}
    seen = 0
    log(f"Reading {role} identity catalog from {path}")

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        # endif
        if args.max_events is not None and seen + n > args.max_events:
            n = args.max_events - seen
            arrays = {key: value[:n] for key, value in arrays.items()}
        # endif
        seen += n

        for name in store:
            branch = resolved.get(name)
            default = -1 if name in {"runnum", "eventnum"} else math.nan
            dtype = np.int64 if name in {"runnum", "eventnum"} else np.float64
            store[name].append(finite_array(arrays, branch, default=default, dtype=dtype))
        # endfor
    # endfor

    catalog = IdentityCatalog(
        runnum=concat(store["runnum"], np.int64),
        eventnum=concat(store["eventnum"], np.int64),
        e_p=concat(store["e_p"]),
        e_theta=concat(store["e_theta"]),
        e_phi=concat(store["e_phi"]),
        p_p=concat(store["p_p"]),
        p_theta=concat(store["p_theta"]),
        p_phi=concat(store["p_phi"]),
    )

    keys = identity_keys_data(catalog) if mode == "data" else identity_keys_mc(
        catalog, args.mc_signature_decimals
    )
    unique_count = int(np.unique(structured_keys(keys)).size) if catalog.size() else 0
    diagnostics = {
        "role": role,
        "mode": mode,
        "path": path,
        "tree_entries_read": int(seen),
        "catalog_records": catalog.size(),
        "unique_identity_keys": unique_count,
    }
    return catalog, diagnostics


def identity_keys_data(records) -> np.ndarray:
    return np.column_stack(
        (
            np.asarray(records.runnum, dtype=np.int64),
            np.asarray(records.eventnum, dtype=np.int64),
        )
    )


def identity_keys_mc(records, decimals: int) -> np.ndarray:
    values = (
        records.e_p, records.e_theta, records.e_phi,
        records.p_p, records.p_theta, records.p_phi,
    )
    return np.column_stack(
        [np.round(np.asarray(value, dtype=np.float64), decimals) for value in values]
    )


def identity_keys_epgg_data(records: EpggRecords) -> np.ndarray:
    return np.column_stack((records.runnum.astype(np.int64), records.eventnum.astype(np.int64)))


def identity_keys_epgg_mc(records: EpggRecords, decimals: int) -> np.ndarray:
    return np.column_stack([
        np.round(records.e_p, decimals),
        np.round(records.e_theta, decimals),
        np.round(records.e_phi, decimals),
        np.round(records.p_p, decimals),
        np.round(records.p_theta, decimals),
        np.round(records.p_phi, decimals),
    ])


def structured_keys(matrix: np.ndarray) -> np.ndarray:
    matrix = np.ascontiguousarray(matrix)
    dtype = np.dtype([(f"f{index}", matrix.dtype) for index in range(matrix.shape[1])])
    return matrix.view(dtype).reshape(-1)


def group_indices(keys: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    structured = structured_keys(keys)
    _, inverse, counts = np.unique(structured, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind="stable")
    offsets = np.empty(counts.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    return inverse, order, offsets


def photon_pair_mass(
    E1: float,
    theta1: float,
    phi1: float,
    E2: float,
    theta2: float,
    phi2: float,
) -> float:
    cosine = (
        math.sin(theta1) * math.sin(theta2) * math.cos(phi1 - phi2)
        + math.cos(theta1) * math.cos(theta2)
    )
    mass2 = 2.0 * E1 * E2 * (1.0 - max(-1.0, min(1.0, cosine)))
    return math.sqrt(max(mass2, 0.0))


def match_truth_partners(
    opportunities: OpportunityRecords,
    epgg: EpggRecords,
    mode: str,
    args: argparse.Namespace,
) -> PartnerMatchResult:
    n = opportunities.size()
    matched = np.zeros(n, dtype=bool)
    catalog_confirmed = np.zeros(n, dtype=bool)
    partner_index = np.full(n, -1, dtype=np.int64)
    probe_angle = np.full(n, np.nan)
    probe_relative_E = np.full(n, np.nan)
    pair_mass = np.full(n, np.nan)
    score = np.full(n, np.nan)
    partner_E = np.full(n, np.nan)
    tag_angle = np.full(n, np.nan)
    tag_relative_E = np.full(n, np.nan)
    closure = np.full(n, np.nan)

    opportunity_keys = identity_keys_data(opportunities) if mode == "data" else identity_keys_mc(opportunities, args.mc_signature_decimals)
    native_keys = identity_keys_epgg_data(epgg) if mode == "data" else identity_keys_epgg_mc(epgg, args.mc_signature_decimals)

    native_structured = structured_keys(native_keys)
    unique, inverse, counts = np.unique(native_structured, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind="stable")
    offsets = np.empty(counts.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])
    group_by_key = {unique[i].tobytes(): i for i in range(unique.size)}

    identity_matches = valid_opportunities = tag_matches = partner_threshold = probe_passes = 0
    records_tested = solutions_tested = 0

    for oi, key in enumerate(structured_keys(opportunity_keys)):
        group = group_by_key.get(key.tobytes())
        if group is None:
            continue
        # endif
        catalog_confirmed[oi] = True
        identity_matches += 1
        best = None
        any_valid = False
        members = order[offsets[group]:offsets[group + 1]]
        for ni in members:
            records_tested += 1
            for si in range(epgg.solution_valid.shape[1]):
                if not epgg.solution_valid[ni, si]:
                    continue
                # endif
                any_valid = True
                solutions_tested += 1
                assignments = (
                    (epgg.g1_E[ni,si], epgg.g1_theta[ni,si], epgg.g1_phi[ni,si],
                     epgg.g2_E[ni,si], epgg.g2_theta[ni,si], epgg.g2_phi[ni,si]),
                    (epgg.g2_E[ni,si], epgg.g2_theta[ni,si], epgg.g2_phi[ni,si],
                     epgg.g1_E[ni,si], epgg.g1_theta[ni,si], epgg.g1_phi[ni,si]),
                )
                for tE,tth,tph,pE,pth,pph in assignments:
                    ta = float(opening_angle_deg(
                        np.asarray([opportunities.tag_theta[oi]]), np.asarray([opportunities.tag_phi[oi]]),
                        np.asarray([tth]), np.asarray([tph]))[0])
                    tr = abs(opportunities.tag_E[oi]-tE)/max(tE,1e-12)
                    ts = (ta/max(args.tag_match_angle_max_deg,1e-12))**2 + (tr/max(args.tag_match_relative_E_max,1e-12))**2
                    if best is None or ts < best[0]:
                        pa = float(opening_angle_deg(
                            np.asarray([opportunities.probe_theta[oi]]), np.asarray([opportunities.probe_phi[oi]]),
                            np.asarray([pth]), np.asarray([pph]))[0])
                        pr = abs(opportunities.probe_E[oi]-pE)/max(pE,1e-12)
                        pm = photon_pair_mass(float(opportunities.tag_E[oi]),float(opportunities.tag_theta[oi]),float(opportunities.tag_phi[oi]),float(pE),float(pth),float(pph))
                        best=(ts,ta,tr,pa,pr,pm,float(pE),int(ni),float(epgg.solution_closure[ni,si]))
                    # endif
                # endfor
            # endfor
        # endfor
        if any_valid:
            valid_opportunities += 1
        # endif
        if best is None:
            continue
        # endif
        ts,ta,tr,pa,pr,pm,pE,ni,cl = best
        partner_index[oi]=ni; probe_angle[oi]=pa; probe_relative_E[oi]=pr
        pair_mass[oi]=pm; score[oi]=ts; partner_E[oi]=pE
        tag_angle[oi]=ta; tag_relative_E[oi]=tr; closure[oi]=cl
        if not (ta < args.tag_match_angle_max_deg and tr < args.tag_match_relative_E_max):
            continue
        # endif
        tag_matches += 1
        if pE < args.found_probe_E_min:
            continue
        # endif
        partner_threshold += 1
        probe_ok = pa < args.probe_match_angle_max_deg and pr < args.probe_match_relative_E_max
        probe_passes += int(probe_ok)
        if args.require_probe_residual_match and not probe_ok:
            continue
        # endif
        matched[oi]=True
    # endfor

    summary={
        "mode":mode,"opportunities":n,"identity_matches":identity_matches,
        "valid_solution_opportunities":valid_opportunities,"tag_matched":tag_matches,
        "partner_above_threshold":partner_threshold,"probe_residual_passes":probe_passes,
        "matched_opportunities":int(np.count_nonzero(matched)),
        "epgammagamma_records_tested":records_tested,"daughter_solutions_tested":solutions_tested,
        "fraction_identity_matched":identity_matches/n if n else math.nan,
        "fraction_matched":float(np.count_nonzero(matched))/n if n else math.nan,
        "conditional_valid_solution_given_identity":valid_opportunities/identity_matches if identity_matches else math.nan,
        "conditional_tag_match_given_valid_solution":tag_matches/valid_opportunities if valid_opportunities else math.nan,
        "require_probe_residual_match":bool(args.require_probe_residual_match),
    }
    return PartnerMatchResult(
        matched,catalog_confirmed,partner_index,probe_angle,probe_relative_E,pair_mass,score,
        partner_E,tag_angle,tag_relative_E,closure,summary
    )



def category_mask(
    records: OpportunityRecords,
    detector: str,
    sector: int,
) -> np.ndarray:
    if detector == "FT":
        return records.probe_detector == 0
    # endif
    return (records.probe_detector == 1) & (records.probe_sector == sector)


def variable_values(records: OpportunityRecords, key: str) -> np.ndarray:
    return np.asarray(getattr(records, key), dtype=np.float64)


def histograms_for_mask(
    records: OpportunityRecords,
    mask: np.ndarray,
) -> Dict[str, np.ndarray]:
    histograms: Dict[str, np.ndarray] = {}
    for variable in FIT_VARIABLES:
        values = variable_values(records, variable.key)[mask]
        values = values[np.isfinite(values)]
        histograms[variable.key], _ = np.histogram(
            values,
            bins=variable.bins,
            range=(variable.low, variable.high),
        )
    # endfor
    return histograms







def plot_expected_probe_diagnostics(
    path: Path,
    period_label: str,
    samples: Mapping[str, OpportunityRecords],
) -> None:
    variables = (
        ("tag_E", np.linspace(0.4, 9.5, 92), r"$E_{\rm tag}$ (GeV)"),
        ("probe_E", np.linspace(2.0, 9.5, 76), r"$E_{\rm probe}^{\rm pred}$ (GeV)"),
        ("probe_theta", np.linspace(math.radians(2.5), math.radians(35.0), 80), r"$\theta_{\rm probe}^{\rm pred}$ (rad)"),
        ("probe_phi", np.linspace(0.0, 2.0 * math.pi, 73), r"$\phi_{\rm probe}^{\rm pred}$ (rad)"),
        ("probe_m2", np.linspace(-0.1, 0.1, 101), r"$m_X^2$ (GeV$^2$)"),
        ("probe_E_minus_p", np.linspace(-0.1, 0.1, 101), r"$E_X-|\vec p_X|$ (GeV)"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for axis, (name, bins, label) in zip(axes.flat, variables):
        for sample, records in samples.items():
            values = np.asarray(getattr(records, name), dtype=float)
            values = values[np.isfinite(values)]
            if values.size == 0:
                continue
            # endif
            axis.hist(
                values,
                bins=bins,
                weights=np.ones(values.size) / values.size,
                histtype="step",
                linewidth=1.3,
                label=f"{sample} ({values.size:,})",
            )
        # endfor
        axis.set_xlabel(label)
        axis.set_ylabel("Fraction / bin")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    # endfor
    fig.suptitle(f"{period_label}: selected expected-probe opportunities")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)

def plot_matching_residuals(
    path: Path,
    period_label: str,
    data_match: PartnerMatchResult,
    mc_match: PartnerMatchResult,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    panels = (
        ("tag_angle_deg", "Tag angular residual (deg)", 0.0, 5.0),
        ("tag_relative_E", "Tag relative energy residual", 0.0, 0.5),
        ("probe_angle_deg", "Predicted-X to partner angle (deg)", 0.0, 20.0),
        ("probe_relative_E", "Predicted-X to partner relative energy", 0.0, 1.5),
        ("partner_E", "Reconstructed partner energy (GeV)", 0.0, 10.0),
        ("solution_closure", "Native daughter closure", 0.0, 0.10),
    )
    for axis,(field,label,low,high) in zip(axes.flat,panels):
        for sample,result in (("Data",data_match),("AAOGEN MC",mc_match)):
            values=np.asarray(getattr(result,field),dtype=float)
            values=values[np.isfinite(values)]
            if values.size:
                hist,edges=np.histogram(values,bins=np.linspace(low,high,101))
                if hist.sum()>0: hist=hist/hist.sum()
                centers=0.5*(edges[:-1]+edges[1:])
                axis.step(centers,hist,where="mid",label=sample)
            # endif
        # endfor
        axis.set_xlabel(label); axis.set_ylabel("Normalized entries")
        axis.set_xlim(low,high); axis.grid(alpha=0.25); axis.legend(fontsize=8)
    # endfor
    fig.suptitle(f"{period_label}: native truth-partner matching")
    fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(path,dpi=180); plt.close(fig)

def plot_cutflow(path: Path, title: str, cutflow: Mapping[str, object]) -> None:
    rows = cutflow.get("stages", [])
    if not rows:
        return
    # endif
    labels = [str(row["stage"]) for row in rows]
    counts = np.asarray([float(row["count"]) for row in rows])
    total_survival = np.asarray(
        [float(row["total_survival"]) for row in rows]
    )
    incremental = np.asarray(
        [float(row["incremental_survival"]) for row in rows]
    )

    fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True)
    x = np.arange(len(rows))

    axes[0].step(x, np.maximum(counts, 0.5), where="mid", linewidth=1.5)
    axes[0].scatter(x, np.maximum(counts, 0.5), s=25)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Cumulative entries")
    axes[0].grid(alpha=0.25)

    axes[1].plot(x, 100.0 * total_survival, marker="o", label="Total survival")
    axes[1].plot(x, 100.0 * incremental, marker="s", label="Incremental survival")
    axes[1].set_ylabel("Survival (%)")
    axes[1].set_ylim(0.0, 105.0)
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=45, ha="right")

    suspect_stages = {"z_tag", "minus_t1_tag", "open_angle_ep2_tag"}
    for axis in axes:
        for index, stage in enumerate(labels):
            if stage in suspect_stages:
                axis.axvspan(index - 0.45, index + 0.45, alpha=0.10)
            # endif
        # endfor
    # endfor

    largest = cutflow.get("largest_incremental_loss")
    subtitle = ""
    if isinstance(largest, Mapping):
        subtitle = (
            f"Largest incremental loss: {largest['stage']} — "
            f"{100.0 * float(largest['incremental_rejection']):.2f}% rejected"
        )
    # endif
    fig.suptitle(f"{title}\n{subtitle}")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)

def write_cutflow_csv(path: Path, cutflow: Mapping[str, object]) -> None:
    """Write the finalized sequential cutflow table to CSV."""
    rows = list(cutflow.get("stages", []))
    if not rows:
        return
    # endif

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "order",
        "stage",
        "label",
        "condition",
        "branch",
        "scope",
        "applied",
        "count",
        "rejected_at_stage",
        "incremental_survival",
        "incremental_rejection",
        "total_survival",
    ]

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: row.get(field)
                    for field in fieldnames
                }
            )

def write_cutflow_text(path: Path, cutflow: Mapping[str, object]) -> None:
    lines = [
        f"Role: {cutflow.get('role')}",
        f"Path: {cutflow.get('path')}",
        "",
        (
            "Columns: cumulative count | incremental survival | total survival "
            "| rejected at stage | cut description"
        ),
        "",
    ]
    for row in cutflow.get("stages", []):
        branch = row.get("branch")
        branch_text = f"; branch={branch}" if branch else ""
        applied_text = "applied" if row.get("applied") else "NOT APPLIED (branch/option absent)"
        lines.append(
            f"{int(row['order']):2d} {row['stage']:<34s} "
            f"{int(row['count']):>12,d} | "
            f"{100.0 * float(row['incremental_survival']):>8.3f}% | "
            f"{100.0 * float(row['total_survival']):>8.3f}% | "
            f"rejected={int(row['rejected_at_stage']):>10,d} | "
            f"{row['condition']} [{applied_text}; scope={row['scope']}"
            f"{branch_text}]"
        )
    # endfor

    lines.extend(
        [
            "",
            "Predicted-probe detector populations:",
            json.dumps(cutflow.get("detector_counts", {}), indent=2),
            "",
            "NOTE: z_tag, minus_t1_tag, and open_angle_ep2_tag are intentionally "
            "NOT APPLIED. They remain in the stage table with 100% incremental "
            "survival so the removed selection is explicit in the provenance.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def emit_cutflow_diagnostics(
    period_dir: Path,
    period_label: str,
    named_cutflows: Mapping[str, Mapping[str, object]],
) -> None:
    audit_dir = period_dir / "cutflow_diagnostics"
    audit_dir.mkdir(parents=True, exist_ok=True)

    for sample_name, cutflow in named_cutflows.items():
        safe = sample_name.replace(" ", "_").replace("/", "_")
        write_cutflow_csv(audit_dir / f"{safe}_cutflow.csv", cutflow)
        write_cutflow_text(audit_dir / f"{safe}_cutflow.txt", cutflow)
        plot_cutflow(
            audit_dir / f"{safe}_cutflow.png",
            f"{period_label}: {sample_name} opportunity cutflow",
            cutflow,
        )

        largest = cutflow.get("largest_incremental_loss")
        if isinstance(largest, Mapping):
            log(
                f"{period_label} {sample_name} cutflow: "
                f"{int(largest['count']):,} survive {largest['stage']}; "
                f"incremental survival="
                f"{100.0 * float(largest['incremental_survival']):.3f}%, "
                f"total survival="
                f"{100.0 * float(largest['total_survival']):.3f}%"
            )
        # endif

        detector_counts = cutflow.get("detector_counts", {})
        for row in cutflow.get("stages", []):
            detector_suffix = ""
            if row["stage"] in detector_counts:
                stage_detector_counts = detector_counts[row["stage"]]
                detector_suffix = (
                    f"; FT={int(stage_detector_counts.get('FT', 0)):,}, "
                    f"FD={int(stage_detector_counts.get('FD', 0)):,}"
                )
            # endif
            log(
                f"{period_label} {sample_name}: "
                f"{row['stage']}: {int(row['count']):,} "
                f"({100.0 * float(row['incremental_survival']):.3f}% of prior, "
                f"{100.0 * float(row['total_survival']):.3f}% of input)"
                f"{detector_suffix} — {row['condition']}"
            )


@dataclass
class FailOnlyEfficiencyRow:
    period: str
    period_label: str
    detector: str
    sector: int
    data_pass_count: int
    data_fail_count: int
    pi0_mc_pass_count: int
    pi0_mc_fail_count: int
    dvcs_mc_count: int
    fit_fraction_pi0_fail: float
    fit_fraction_pi0_fail_stat_err: float
    fit_fraction_pi0_fail_model_spread: float
    fitted_pi0_fail: float
    fitted_pi0_fail_stat_err: float
    fitted_pi0_fail_model_spread: float
    efficiency_data: float
    efficiency_data_stat_err: float
    efficiency_data_model_err: float
    efficiency_mc: float
    efficiency_mc_err: float
    scale_factor: float
    scale_factor_stat_err: float
    scale_factor_model_err: float
    fit_success: bool
    fit_message: str
    fit_deviance: float
    fit_ndf: int
    fit_reduced_deviance: float
    fit_quality_warning: bool
    per_variable_fractions: str

@dataclass
class FractionFit:
    success: bool
    message: str
    fraction: float
    stat_error: float
    nll: float
    deviance: float
    ndf: int
    reduced_deviance: float
    per_variable: Dict[str, Dict[str, float]]
    payload: Dict[str, Dict[str, object]]
    warnings: List[str]

FAIL_FIT_DRIVERS: Tuple[str,...] = FRACTION_DRIVERS

def parse_args() -> argparse.Namespace:
    p=argparse.ArgumentParser(description='FAIL-only photon efficiency template extraction')
    p.add_argument('--period',action='append',choices=[x.key for x in PERIODS])
    p.add_argument('--workers',type=int,default=5)
    p.add_argument('--step-size',default=DEFAULT_STEP_SIZE); p.add_argument('--max-events',type=int,default=None)
    p.add_argument('--output-dir',default='output/photon_efficiency_study/fail_only_template_fit')
    p.add_argument('--tag-E-min',type=float,default=.40); p.add_argument('--probe-E-min',type=float,default=2.0); p.add_argument('--probe-E-max',type=float,default=9.5)
    p.add_argument('--probe-m2-abs-max',type=float,default=.10); p.add_argument('--probe-E-minus-p-abs-max',type=float,default=.10)
    p.add_argument('--ft-theta-min',type=float,default=2.5); p.add_argument('--ft-theta-max',type=float,default=5.0); p.add_argument('--fd-theta-min',type=float,default=5.0); p.add_argument('--fd-theta-max',type=float,default=35.0)
    p.add_argument('--Q2-min',type=float,default=1.0); p.add_argument('--W-min',type=float,default=2.0); p.add_argument('--y-max',type=float,default=.8)
    p.add_argument('--z-min',type=float,default=.65); p.add_argument('--minus-t-max',type=float,default=1.0); p.add_argument('--open-angle-min-deg',type=float,default=5.0); p.add_argument('--require-fiducial-status-111',action='store_true')
    p.add_argument('--tag-match-angle-max-deg',type=float,default=3.0); p.add_argument('--tag-match-relative-E-max',type=float,default=.35); p.add_argument('--probe-match-angle-max-deg',type=float,default=3.0); p.add_argument('--probe-match-relative-E-max',type=float,default=.35)
    p.add_argument('--mc-signature-decimals',type=int,default=10); p.add_argument('--native-closure-max',type=float,default=.10); p.add_argument('--native-minimum-photon-E',type=float,default=.20); p.add_argument('--found-probe-E-min',type=float,default=2.0); p.add_argument('--require-probe-residual-match',action='store_true')
    p.add_argument(
        '--fit-drivers',
        nargs='+',
        default=list(FAIL_FIT_DRIVERS),
        choices=[v.key for v in FIT_VARIABLES],
        help=(
            'Variables used in the nominal FAIL-only fit. The default is '
            'theta_gamma_gamma only. All seven observables are still plotted '
            'as diagnostics.'
        ),
    )
    p.add_argument('--template-pseudocount',type=float,default=.25); p.add_argument('--fit-min-data-fail',type=int,default=100); p.add_argument('--fit-min-template-counts',type=int,default=100)
    p.add_argument('--fit-max-reduced-deviance',type=float,default=5.0)
    p.add_argument('--fit-max-variable-spread',type=float,default=.15)
    p.add_argument(
        '--exclusivity-fitter-script',
        default=str(
            Path(__file__).resolve().parent.parent
            / 'plot_exclusivity_data_dvcs_pi0_mc.py'
        ),
        help=(
            'Path to the validated exclusivity-template script. Its exact '
            'fit_shared_two_templates implementation is used for the FAIL fit.'
        ),
    )
    p.add_argument('--max-shift-bins', type=float, default=12.0)
    p.add_argument('--max-smear-bins', type=float, default=20.0)
    p.add_argument('--shift-prior-bins', type=float, default=4.0)
    p.add_argument('--smear-prior-bins', type=float, default=8.0)
    p.add_argument('--disable-nuisance-penalties', action='store_true')
    p.add_argument('--fit-core-containment', type=float, default=0.90)
    p.add_argument('--fit-fraction-containment', type=float, default=0.95)
    p.add_argument('--pi0-support-core-containment', type=float, default=0.90)
    p.add_argument('--pi0-support-fraction-containment', type=float, default=0.95)
    p.add_argument('--outside-overshoot-penalty-weight', type=float, default=0.25)
    p.add_argument('--emiss2-mean-order-penalty-weight', type=float, default=25.0)
    return p.parse_args()

def category_stem(detector,sector): return 'ft' if detector=='FT' else f'fd_sector_{sector}'
def category_title(detector,sector): return 'FT' if detector=='FT' else f'FD sector {sector}'
def finite_values(r,key,mask):
    x=np.asarray(getattr(r,key),float)[mask]; return x[np.isfinite(x)]

def histogram_with_flow(r,mask,var):
    x=finite_values(r,var.key,mask); core,_=np.histogram(x,bins=var.bins,range=(var.low,var.high))
    return np.concatenate(([np.count_nonzero(x<var.low)],core,[np.count_nonzero(x>=var.high)])).astype(float)

def shape(counts,pseudo):
    a=np.asarray(counts,float)+max(float(pseudo),0.0); return a/a.sum()

def multinomial_nll(f,data,p0,p1):
    prob=np.maximum(f*p0+(1-f)*p1,1e-300); return float(-np.sum(xlogy(data,prob)))

def fit_one_fraction(data,pi0,dvcs):
    obj=lambda z: multinomial_nll(float(z),data,pi0,dvcs)
    res=minimize_scalar(obj,bounds=(1e-6,1-1e-6),method='bounded',options={'xatol':1e-10})
    f=float(res.x); h=1e-4; h=min(h,.45*f,.45*(1-f)); h=max(h,1e-7)
    curv=(obj(f+h)-2*obj(f)+obj(f-h))/(h*h); err=math.sqrt(1/curv) if curv>0 and math.isfinite(curv) else math.nan
    mu=data.sum()*(f*pi0+(1-f)*dvcs); dev=poisson_deviance(data,mu); ndf=max(data.size-1,1)
    return f,err,float(res.fun),dev,ndf,bool(res.success),str(res.message)

_EXCLUSIVITY_FITTER_CACHE = {}


def resolve_exclusivity_fitter_path(script_path: str) -> Path:
    """Resolve the validated fitter path with explicit, deterministic fallbacks."""
    requested = Path(script_path).expanduser()
    candidates = [
        requested,
        Path(__file__).resolve().parent.parent
        / "plot_exclusivity_data_dvcs_pi0_mc.py",
        Path(__file__).resolve().parent
        / "plot_exclusivity_data_dvcs_pi0_mc.py",
    ]

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
        # endif
    # endfor

    attempted = "\n  ".join(str(candidate.resolve()) for candidate in candidates)
    raise FileNotFoundError(
        "Validated exclusivity fitter was not found. Tried:\n  " + attempted
    )


def load_exclusivity_fitter(script_path: str):
    """
    Load and validate the exact template-morphing implementation.

    The module is registered in sys.modules before exec_module(), which is
    required by dataclasses and by imported code that resolves its own module.
    Each worker process caches the successfully imported module.
    """
    path = resolve_exclusivity_fitter_path(script_path)
    cache_key = str(path)
    if cache_key in _EXCLUSIVITY_FITTER_CACHE:
        return _EXCLUSIVITY_FITTER_CACHE[cache_key]
    # endif

    module_name = "validated_exclusivity_fitter"
    specification = importlib.util.spec_from_file_location(module_name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"Could not create an import specification for {path}")
    # endif

    module = importlib.util.module_from_spec(specification)
    sys.modules[module_name] = module
    try:
        specification.loader.exec_module(module)
    except Exception:
        sys.modules.pop(module_name, None)
        raise
    # endtry

    required_symbols = (
        "fit_shared_two_templates",
        "VARIABLES",
        "TopologyConfig",
    )
    missing = [name for name in required_symbols if not hasattr(module, name)]
    if missing:
        raise RuntimeError(
            f"{path} is missing required fitter symbols: " + ", ".join(missing)
        )
    # endif

    expected_parameters = {
        "data_histograms",
        "dvcs_histograms",
        "pi0_histograms",
        "topology",
        "max_shift_bins",
        "max_smear_bins",
        "min_counts",
        "fraction_variable_branches",
        "shift_prior_bins",
        "smear_prior_bins",
        "use_nuisance_penalties",
        "core_containment",
        "fraction_containment",
        "pi0_support_core_containment",
        "pi0_support_fraction_containment",
        "pi0_core_calibration",
        "outside_overshoot_penalty_weight",
        "emiss2_mean_order_penalty_weight",
    }
    signature = inspect.signature(module.fit_shared_two_templates)
    available_parameters = set(signature.parameters)
    missing_parameters = sorted(expected_parameters - available_parameters)
    if missing_parameters:
        raise RuntimeError(
            "Validated fit_shared_two_templates() interface is incompatible. "
            "Missing parameter(s): " + ", ".join(missing_parameters)
            + f". Resolved fitter: {path}"
        )
    # endif

    variables = list(module.VARIABLES)
    variable_branches = {variable.branch for variable in variables}
    missing_drivers = sorted(
        set(FAIL_FIT_DRIVERS) - variable_branches
    )
    if missing_drivers:
        raise RuntimeError(
            "Validated fitter does not define nominal FAIL driver(s): "
            + ", ".join(missing_drivers)
        )
    # endif

    # Lightweight interface exercise: instantiate the topology object now,
    # before any ROOT files are read.
    module.TopologyConfig(
        key="FAIL_ONLY_PREFLIGHT",
        label="FAIL-only preflight",
        detector1=-1,
        detector2=-1,
    )

    module._resolved_fitter_path = str(path)
    module._validated_fit_signature = str(signature)
    _EXCLUSIVITY_FITTER_CACHE[cache_key] = module
    return module


def preflight_exclusivity_fitter(args: argparse.Namespace) -> Dict[str, object]:
    """Validate and report the external fitter before expensive processing."""
    module = load_exclusivity_fitter(args.exclusivity_fitter_script)
    payload = {
        "resolved_path": module._resolved_fitter_path,
        "fit_signature": module._validated_fit_signature,
        "variables": [variable.branch for variable in module.VARIABLES],
        "nominal_fraction_drivers": list(args.fit_drivers),
        "module_name": module.__name__,
    }
    log(
        "Validated exclusivity fitter interface: "
        f"{payload['resolved_path']} | signature "
        f"{payload['fit_signature']}"
    )
    return payload


def exclusivity_variable_values(
    records: OpportunityRecords,
    branch: str,
    mask: np.ndarray,
) -> np.ndarray:
    """
    Map this script's stored quantities onto the validated fitter's conventions.

    In particular, theta_gamma_gamma is stored here in degrees but is fitted in
    radians by plot_exclusivity_data_dvcs_pi0_mc.py.
    """
    source_branch = "theta_cm" if branch == "theta" else branch
    values = finite_values(records, source_branch, mask)

    if branch == "theta_gamma_gamma":
        values = np.radians(values)
    # endif
    return values


def build_exclusivity_histograms(
    fitter_module,
    records: OpportunityRecords,
    mask: np.ndarray,
) -> Dict[str, np.ndarray]:
    histograms: Dict[str, np.ndarray] = {}
    for variable in fitter_module.VARIABLES:
        values = exclusivity_variable_values(
            records,
            variable.branch,
            mask,
        )
        counts, _ = np.histogram(
            values,
            bins=variable.bins,
            range=(variable.xmin, variable.xmax),
        )
        histograms[variable.branch] = counts.astype(float)
    # endfor
    return histograms


def fit_fail_shapes(data, dvcs, pi0, data_fail, dvcs_mask, pi0_fail, args):
    """
    Fit the FAIL sample using the exact shape-morphing model from the validated
    exclusivity determination.

    The PASS count remains fixed and external to this fit. Only the FAIL pi0
    fraction is extracted. Both the BH/DVCS and AAOGEN FAIL templates receive
    the same variable-dependent shift/smearing treatment, nuisance priors,
    support-region handling, and physical penalties used by the exclusivity
    script.
    """
    n_fail = int(np.count_nonzero(data_fail))
    if n_fail < args.fit_min_data_fail:
        raise RuntimeError(
            f"insufficient FAIL data: {n_fail} < {args.fit_min_data_fail}"
        )
    # endif

    fitter = load_exclusivity_fitter(args.exclusivity_fitter_script)

    data_histograms = build_exclusivity_histograms(
        fitter, data, data_fail
    )
    dvcs_histograms = build_exclusivity_histograms(
        fitter, dvcs, dvcs_mask
    )
    pi0_histograms = build_exclusivity_histograms(
        fitter, pi0, pi0_fail
    )

    topology = fitter.TopologyConfig(
        key="FAIL_ONLY",
        label="FAIL-only combined proton topologies",
        detector1=-1,
        detector2=-1,
    )

    summary = fitter.fit_shared_two_templates(
        data_histograms=data_histograms,
        dvcs_histograms=dvcs_histograms,
        pi0_histograms=pi0_histograms,
        topology=topology,
        max_shift_bins=float(args.max_shift_bins),
        max_smear_bins=float(args.max_smear_bins),
        min_counts=int(args.fit_min_template_counts),
        fraction_variable_branches=tuple(args.fit_drivers),
        shift_prior_bins=float(args.shift_prior_bins),
        smear_prior_bins=float(args.smear_prior_bins),
        use_nuisance_penalties=not args.disable_nuisance_penalties,
        core_containment=float(args.fit_core_containment),
        fraction_containment=float(args.fit_fraction_containment),
        pi0_support_core_containment=float(
            args.pi0_support_core_containment
        ),
        pi0_support_fraction_containment=float(
            args.pi0_support_fraction_containment
        ),
        pi0_core_calibration=None,
        outside_overshoot_penalty_weight=float(
            args.outside_overshoot_penalty_weight
        ),
        emiss2_mean_order_penalty_weight=float(
            args.emiss2_mean_order_penalty_weight
        ),
    )

    if not summary.success:
        raise RuntimeError(
            "Validated exclusivity-template fitter failed: "
            + str(summary.message)
        )
    # endif

    fraction = float(summary.f_pi0)
    stat_error = float(summary.f_pi0_err)
    variable_results = summary.variable_results or {}

    payload: Dict[str, Dict[str, object]] = {}
    per_variable: Dict[str, Dict[str, float]] = {}

    # Retain every successfully profiled variable for diagnostics, while the
    # requested fit drivers alone determine the shared fraction.
    for branch, result in variable_results.items():
        if (
            result.fit_data_counts is None
            or result.model_counts is None
            or result.dvcs_component_counts is None
            or result.pi0_component_counts is None
        ):
            continue
        # endif

        fitted_variable = next(
            variable
            for variable in fitter.VARIABLES
            if variable.branch == branch
        )
        local_key = "theta_cm" if branch == "theta" else branch
        local_variable = next(
            (
                variable
                for variable in FIT_VARIABLES
                if variable.key == local_key
            ),
            FitVariable(
                local_key,
                fitted_variable.label,
                fitted_variable.bins,
                fitted_variable.xmin,
                fitted_variable.xmax,
            ),
        )

        data_counts = np.asarray(result.fit_data_counts, dtype=float)
        model = np.asarray(result.model_counts, dtype=float)
        dvcs_component = np.asarray(
            result.dvcs_component_counts,
            dtype=float,
        )
        pi0_component = np.asarray(
            result.pi0_component_counts,
            dtype=float,
        )
        residual = data_counts - model
        pull = residual / np.sqrt(np.maximum(model, 1.0))
        deviance = float(result.deviance)
        ndf = max(int(result.ndf), 1)

        payload[local_key] = {
            "variable": local_variable,
            "data": data_counts,
            "model": model,
            "pi0_component": pi0_component,
            "dvcs_component": dvcs_component,
            "residual": residual,
            "pull": pull,
            "combined_deviance": deviance,
            "combined_ndf": ndf,
            "combined_reduced_deviance": deviance / ndf,
            "pi0_raw": pi0_histograms[branch],
            "dvcs_raw": dvcs_histograms[branch],
            "pi0_shape": np.asarray(
                result.transformed_pi0_shape,
                dtype=float,
            ),
            "dvcs_shape": np.asarray(
                result.transformed_dvcs_shape,
                dtype=float,
            ),
        }

        per_variable[local_key] = {
            "fraction": fraction,
            "stat_error": stat_error,
            "nll": math.nan,
            "deviance_at_individual_optimum": deviance,
            "ndf": ndf,
            "reduced_deviance_at_individual_optimum": deviance / ndf,
            "success": bool(result.success),
            "message": str(result.message),
            "deviance_at_combined_fraction": deviance,
            "reduced_deviance_at_combined_fraction": deviance / ndf,
            "fraction_minus_combined": 0.0,
            "shift": float(result.shift),
            "shift_error": float(result.shift_err),
            "smear": float(result.sigma_add),
            "smear_error": float(result.sigma_add_err),
            "pi0_shift": float(result.pi0_shift),
            "pi0_shift_error": float(result.pi0_shift_err),
            "pi0_smear": float(result.pi0_sigma_add),
            "pi0_smear_error": float(result.pi0_sigma_add_err),
            "morph_label": str(result.morph_label),
            "fit_region_data_counts": float(result.fit_region_data_counts),
            "fit_region_model_counts": float(result.fit_region_model_counts),
            "full_range_model_to_data": float(
                result.full_range_model_to_data
            ),
        }
    # endfor

    reduced_deviance = (
        float(summary.deviance) / int(summary.ndf)
        if int(summary.ndf) > 0
        else math.nan
    )
    warnings: List[str] = []
    if (
        args.fit_max_reduced_deviance > 0.0
        and math.isfinite(reduced_deviance)
        and reduced_deviance > args.fit_max_reduced_deviance
    ):
        warnings.append(
            f"D/ndf={reduced_deviance:.3f} exceeds "
            f"{args.fit_max_reduced_deviance:g}"
        )
    # endif

    return FractionFit(
        True,
        str(summary.message),
        fraction,
        stat_error,
        math.nan,
        float(summary.deviance),
        int(summary.ndf),
        reduced_deviance,
        per_variable,
        payload,
        warnings,
    )


def poisson_deviance(obs,mu):
    obs=np.asarray(obs,float); mu=np.maximum(np.asarray(mu,float),1e-12); term=np.where(obs>0,mu-obs+obs*np.log(obs/mu),mu); return float(2*np.sum(term))

def binomial_efficiency(passed,total):
    if total<=0 or passed<0 or passed>total:return math.nan,math.nan
    e=passed/total; return e,math.sqrt(max(e*(1-e)/total,0))

def efficiency_from_pass_fail(npass,nfail,nfail_err):
    den=npass+nfail
    if den<=0:return math.nan,math.nan
    e=npass/den; dp=nfail/(den*den); df=-npass/(den*den); var=dp*dp*max(npass,1.0)+df*df*nfail_err*nfail_err; return e,math.sqrt(max(var,0))

def ratio_errors(a,ae,b,be):
    if not all(map(math.isfinite,[a,ae,b,be])) or b<=0:return math.nan,math.nan
    r=a/b; return r,abs(r)*math.sqrt((ae/a)**2+(be/b)**2) if a>0 else math.nan

def plot_shape_diagnostics(out,period,det,sec,data,dvcs,pi0,dp,df,pp,pf,dm):
    out.mkdir(parents=True,exist_ok=True)
    for var in FIT_VARIABLES:
        fig,ax=plt.subplots(2,2,figsize=(15,10)); panels=[(ax[0,0],[('Data PASS',data,dp),('Data FAIL',data,df)],'Data PASS vs FAIL'),(ax[0,1],[('Matched data',data,dp),('AAOGEN PASS',pi0,pp)],'PASS closure'),(ax[1,0],[('Data FAIL',data,df),('AAOGEN FAIL',pi0,pf),('DVCSGEN',dvcs,dm)],'FAIL templates'),(ax[1,1],[('Data PASS',data,dp),('AAOGEN PASS',pi0,pp)],'PASS validation')]
        for a,curves,title in panels:
            for lab,r,m in curves:
                vals=finite_values(r,var.key,m); c,e=np.histogram(vals,bins=var.bins,range=(var.low,var.high)); dens=c/c.sum() if c.sum() else c.astype(float); ctr=(e[:-1]+e[1:])/2; a.step(ctr,dens,where='mid',label=f'{lab} ({np.count_nonzero(m):,})')
            a.set_title(title); a.set_xlabel(var.label); a.set_ylabel('Fraction / bin'); a.grid(alpha=.25); a.legend(fontsize=8)
        fig.suptitle(f'{period}, {category_title(det,sec)}: {var.label}'); fig.tight_layout(rect=(0,0,1,.95)); fig.savefig(out/f'{var.key}_pass_fail_shapes.png',dpi=180); plt.close(fig)

def flow_bin_labels(variable: FitVariable, size: Optional[int] = None) -> List[str]:
    count = variable.bins if size is None else int(size)
    if count == variable.bins + 2:
        labels = ["UF"]
        labels.extend(str(index + 1) for index in range(variable.bins))
        labels.append("OF")
        return labels
    # endif
    return [str(index + 1) for index in range(count)]


def plot_one_fail_fit_variable(
    path: Path,
    period: str,
    detector: str,
    sector: int,
    key: str,
    fit: FractionFit,
) -> None:
    """
    Plot one FAIL-fit observable with stacked components and a residual panel.
    """
    item = fit.payload[key]
    variable = item["variable"]
    data_counts = np.asarray(item["data"], dtype=float)
    pi0_component = np.asarray(item["pi0_component"], dtype=float)
    dvcs_component = np.asarray(item["dvcs_component"], dtype=float)
    model = np.asarray(item["model"], dtype=float)
    residual = np.asarray(item["residual"], dtype=float)

    x = np.arange(data_counts.size, dtype=float)
    labels = flow_bin_labels(variable, data_counts.size)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(11, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.2]},
    )

    axes[0].bar(
        x,
        dvcs_component,
        width=0.92,
        label="Fitted BH/DVCS FAIL",
        alpha=0.65,
    )
    axes[0].bar(
        x,
        pi0_component,
        width=0.92,
        bottom=dvcs_component,
        label=r"Fitted $\pi^0$ FAIL",
        alpha=0.65,
    )
    axes[0].step(
        x,
        model,
        where="mid",
        linewidth=1.5,
        label="Total fit",
    )
    axes[0].errorbar(
        x,
        data_counts,
        yerr=np.sqrt(np.maximum(data_counts, 1.0)),
        fmt=".",
        markersize=5,
        capsize=1.5,
        label="Data FAIL",
    )
    axes[0].set_ylabel("Events / bin")
    axes[0].set_title(
        f"{period}, {category_title(detector, sector)}: {variable.label}\n"
        f"combined f_pi0={fit.fraction:.5f}; "
        f"individual f_pi0={fit.per_variable[key]['fraction']:.5f}; "
        f"D/ndf={item['combined_reduced_deviance']:.3f}"
    )
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend(fontsize=9)

    axes[1].axhline(0.0, linewidth=1.0)
    axes[1].bar(x, residual, width=0.75)
    axes[1].set_ylabel("Data - fit")
    axes[1].set_xlabel(
        f"{variable.label} bin (UF = underflow, OF = overflow)"
    )
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=90 if len(labels) > 24 else 0)

    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fail_fit(path, period, detector, sector, fit):
    """
    Retain a compact all-driver canvas for quick inspection.

    Dedicated stacked/residual figures are written separately for each driver.
    """
    keys = list(fit.payload)
    fig, axes = plt.subplots(
        2,
        len(keys),
        figsize=(6 * len(keys), 8),
        squeeze=False,
    )

    for column, key in enumerate(keys):
        item = fit.payload[key]
        variable = item["variable"]
        x = np.arange(item["data"].size)

        axes[0, column].bar(
            x,
            item["dvcs_component"],
            width=0.92,
            label="BH/DVCS",
            alpha=0.65,
        )
        axes[0, column].bar(
            x,
            item["pi0_component"],
            width=0.92,
            bottom=item["dvcs_component"],
            label=r"$\pi^0$ FAIL",
            alpha=0.65,
        )
        axes[0, column].step(
            x,
            item["model"],
            where="mid",
            linewidth=1.4,
            label="Total",
        )
        axes[0, column].errorbar(
            x,
            item["data"],
            yerr=np.sqrt(np.maximum(item["data"], 1.0)),
            fmt=".",
            label="Data",
        )
        axes[0, column].set_title(
            f"{variable.label}\n"
            f"D/ndf={item['combined_reduced_deviance']:.2f}"
        )
        axes[0, column].set_ylabel("Events / bin")
        axes[0, column].legend(fontsize=8)
        axes[0, column].grid(axis="y", alpha=0.25)

        axes[1, column].axhline(0.0, linewidth=1.0)
        axes[1, column].bar(
            x,
            item["residual"],
            width=0.75,
        )
        axes[1, column].set_ylabel("Data - fit")
        axes[1, column].set_xlabel("bin: UF, regular bins, OF")
        axes[1, column].grid(axis="y", alpha=0.25)
    # endfor

    fig.suptitle(
        f"{period}, {category_title(detector, sector)} FAIL-only fit: "
        f"f_pi0={fit.fraction:.4f} +/- {fit.stat_error:.4f}, "
        f"D/ndf={fit.reduced_deviance:.3f}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fit_diagnostic_summary(
    path: Path,
    period: str,
    detector: str,
    sector: int,
    fit: FractionFit,
) -> None:
    """
    Summarize per-variable deviance and independently preferred pi0 fractions.
    """
    keys = list(fit.per_variable)
    labels = [
        next(variable.label for variable in FIT_VARIABLES if variable.key == key)
        for key in keys
    ]
    x = np.arange(len(keys), dtype=float)

    combined_deviances = np.asarray(
        [
            fit.per_variable[key]["deviance_at_combined_fraction"]
            for key in keys
        ],
        dtype=float,
    )
    combined_reduced = np.asarray(
        [
            fit.per_variable[key]["reduced_deviance_at_combined_fraction"]
            for key in keys
        ],
        dtype=float,
    )
    fractions = np.asarray(
        [fit.per_variable[key]["fraction"] for key in keys],
        dtype=float,
    )
    fraction_errors = np.asarray(
        [fit.per_variable[key]["stat_error"] for key in keys],
        dtype=float,
    )

    fig, axes = plt.subplots(2, 1, figsize=(11, 9))

    axes[0].bar(x, combined_deviances)
    for index, (deviance, reduced) in enumerate(
        zip(combined_deviances, combined_reduced)
    ):
        axes[0].text(
            index,
            deviance,
            f"D/ndf={reduced:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    # endfor
    axes[0].set_xticks(x, labels, rotation=25, ha="right")
    axes[0].set_ylabel("Deviance contribution")
    axes[0].set_title(
        "Contribution of each fit observable at the combined fraction"
    )
    axes[0].grid(axis="y", alpha=0.25)

    axes[1].errorbar(
        x,
        fractions,
        yerr=fraction_errors,
        fmt="o",
        capsize=3,
        label="Independent one-variable fit",
    )
    axes[1].axhline(
        fit.fraction,
        linewidth=1.5,
        label=f"Combined fit = {fit.fraction:.5f}",
    )
    if math.isfinite(fit.stat_error):
        axes[1].axhspan(
            fit.fraction - fit.stat_error,
            fit.fraction + fit.stat_error,
            alpha=0.15,
            label="Combined statistical uncertainty",
        )
    # endif
    axes[1].set_xticks(x, labels, rotation=25, ha="right")
    axes[1].set_ylabel(r"Fitted $f_{\pi^0,\mathrm{FAIL}}$")
    axes[1].set_ylim(0.0, 1.0)
    axes[1].grid(alpha=0.25)
    axes[1].legend(fontsize=9)

    fig.suptitle(
        f"{period}, {category_title(detector, sector)} fit-diagnostic summary"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fit_diagnostics_payload(
    period: str,
    detector: str,
    sector: int,
    fit: FractionFit,
    n_pass: int,
    n_fail: int,
    n_pi0_mc_pass: int,
    n_pi0_mc_fail: int,
    n_dvcs_mc: int,
) -> Dict[str, object]:
    """
    Return a JSON-serializable record of every fit input and prediction.
    """
    variables: Dict[str, object] = {}
    for key, item in fit.payload.items():
        variable = item["variable"]
        variables[key] = {
            "label": variable.label,
            "range": [variable.low, variable.high],
            "regular_bins": variable.bins,
            "bin_convention": "underflow, regular bins, overflow",
            "data_fail": np.asarray(item["data"], dtype=float).tolist(),
            "pi0_template_raw": np.asarray(
                item["pi0_raw"], dtype=float
            ).tolist(),
            "dvcs_template_raw": np.asarray(
                item["dvcs_raw"], dtype=float
            ).tolist(),
            "pi0_template_probability": np.asarray(
                item["pi0_shape"], dtype=float
            ).tolist(),
            "dvcs_template_probability": np.asarray(
                item["dvcs_shape"], dtype=float
            ).tolist(),
            "fit_total": np.asarray(item["model"], dtype=float).tolist(),
            "fit_pi0_fail": np.asarray(
                item["pi0_component"], dtype=float
            ).tolist(),
            "fit_bh_dvcs_fail": np.asarray(
                item["dvcs_component"], dtype=float
            ).tolist(),
            "residual_data_minus_fit": np.asarray(
                item["residual"], dtype=float
            ).tolist(),
            "pull": np.asarray(item["pull"], dtype=float).tolist(),
            "deviance_at_combined_fraction": item["combined_deviance"],
            "ndf_at_combined_fraction": item["combined_ndf"],
            "reduced_deviance_at_combined_fraction": item[
                "combined_reduced_deviance"
            ],
            "individual_fit": fit.per_variable[key],
        }
    # endfor

    return {
        "period": period,
        "detector": detector,
        "sector": sector,
        "category": category_title(detector, sector),
        "counts": {
            "data_pass": n_pass,
            "data_fail": n_fail,
            "pi0_mc_pass": n_pi0_mc_pass,
            "pi0_mc_fail": n_pi0_mc_fail,
            "dvcs_mc": n_dvcs_mc,
        },
        "combined_fit": {
            "success": fit.success,
            "message": fit.message,
            "fraction_pi0_fail": fit.fraction,
            "fraction_pi0_fail_stat_error": fit.stat_error,
            "nll": fit.nll,
            "deviance": fit.deviance,
            "ndf": fit.ndf,
            "reduced_deviance": fit.reduced_deviance,
            "warnings": fit.warnings,
        },
        "per_variable": fit.per_variable,
        "variables": variables,
    }



def write_components(pathstem, fit):
    pathstem.mkdir(parents=True, exist_ok=True)
    written: Dict[str, str] = {}

    for key, item in fit.payload.items():
        path = pathstem / f"{key}_components.csv"
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "bin_index",
                    "bin_type",
                    "data_fail",
                    "fit_total",
                    "fit_pi0_fail",
                    "fit_bh_dvcs",
                    "residual_data_minus_fit",
                    "pull",
                    "poisson_deviance_contribution",
                ]
            )
            for index in range(item["data"].size):
                observed = float(item["data"][index])
                expected = max(float(item["model"][index]), 1.0e-12)
                deviance_contribution = 2.0 * (
                    expected
                    - observed
                    + (
                        observed * math.log(observed / expected)
                        if observed > 0.0
                        else 0.0
                    )
                )
                if item["data"].size == item["variable"].bins + 2:
                    if index == 0:
                        bin_type = "underflow"
                    elif index == item["data"].size - 1:
                        bin_type = "overflow"
                    else:
                        bin_type = "regular"
                    # endif
                else:
                    bin_type = "regular"
                # endif

                writer.writerow(
                    [
                        index,
                        bin_type,
                        observed,
                        expected,
                        item["pi0_component"][index],
                        item["dvcs_component"][index],
                        item["residual"][index],
                        item["pull"][index],
                        deviance_contribution,
                    ]
                )
            # endfor
        # endwith
        written[key] = str(path)
    # endfor

    return written




def ratio_with_poisson_uncertainty(
    numerator_counts: np.ndarray,
    denominator_counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Ratio of unit-normalized shapes with uncertainties from the raw counts.

    The central value is (n_i/N_n)/(d_i/N_d). The displayed uncertainty uses
    the raw per-bin Poisson terms rather than incorrectly treating normalized
    fractions as event counts.
    """
    numerator_counts = np.asarray(numerator_counts, dtype=float)
    denominator_counts = np.asarray(denominator_counts, dtype=float)
    numerator_total = float(np.sum(numerator_counts))
    denominator_total = float(np.sum(denominator_counts))

    ratio = np.full_like(numerator_counts, np.nan, dtype=float)
    error = np.full_like(numerator_counts, np.nan, dtype=float)
    if numerator_total <= 0.0 or denominator_total <= 0.0:
        return ratio, error
    # endif

    numerator_fraction = numerator_counts / numerator_total
    denominator_fraction = denominator_counts / denominator_total
    valid = denominator_fraction > 0.0
    ratio[valid] = numerator_fraction[valid] / denominator_fraction[valid]

    positive = valid & (numerator_counts > 0.0) & (denominator_counts > 0.0)
    error[positive] = ratio[positive] * np.sqrt(
        1.0 / numerator_counts[positive]
        + 1.0 / denominator_counts[positive]
    )
    zero_numerator = valid & (numerator_counts <= 0.0)
    error[zero_numerator] = 0.0
    return ratio, error


def plot_theta_control_diagnostics(
    path: Path,
    period: str,
    detector: str,
    sector: int,
    data: OpportunityRecords,
    pi0: OpportunityRecords,
    data_pass: np.ndarray,
    pi0_pass: np.ndarray,
    pi0_fail: np.ndarray,
) -> None:
    """
    Validate the theta_gamma_gamma control shapes before interpreting the FAIL fit.

    Panels:
      1. matched PASS data versus AAOGEN PASS;
      2. PASS-data / AAOGEN-PASS ratio;
      3. AAOGEN PASS versus AAOGEN FAIL;
      4. AAOGEN-PASS / AAOGEN-FAIL ratio.

    Coarser binning is used here to expose broad shape differences without
    allowing very fine statistical fluctuations to dominate the comparison.
    """
    variable = next(
        item for item in FIT_VARIABLES
        if item.key == "theta_gamma_gamma"
    )
    bins = 40
    low = variable.low
    high = variable.high
    edges = np.linspace(low, high, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])

    data_pass_counts, _ = np.histogram(
        finite_values(data, "theta_gamma_gamma", data_pass),
        bins=edges,
    )
    pi0_pass_counts, _ = np.histogram(
        finite_values(pi0, "theta_gamma_gamma", pi0_pass),
        bins=edges,
    )
    pi0_fail_counts, _ = np.histogram(
        finite_values(pi0, "theta_gamma_gamma", pi0_fail),
        bins=edges,
    )

    def normalized(counts: np.ndarray) -> np.ndarray:
        total = float(np.sum(counts))
        return counts / total if total > 0.0 else np.zeros_like(counts, dtype=float)
    # enddef

    data_pass_norm = normalized(data_pass_counts.astype(float))
    pi0_pass_norm = normalized(pi0_pass_counts.astype(float))
    pi0_fail_norm = normalized(pi0_fail_counts.astype(float))

    pass_ratio, pass_ratio_err = ratio_with_poisson_uncertainty(
        data_pass_counts,
        pi0_pass_counts,
    )
    mc_ratio, mc_ratio_err = ratio_with_poisson_uncertainty(
        pi0_pass_counts,
        pi0_fail_counts,
    )

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0, 0].step(
        centers,
        data_pass_norm,
        where="mid",
        label=f"Matched data PASS ({int(np.count_nonzero(data_pass)):,})",
    )
    axes[0, 0].step(
        centers,
        pi0_pass_norm,
        where="mid",
        label=f"AAOGEN PASS ({int(np.count_nonzero(pi0_pass)):,})",
    )
    axes[0, 0].set_title("PASS control-shape closure")
    axes[0, 0].set_ylabel("Fraction / bin")
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].grid(alpha=0.25)

    axes[0, 1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[0, 1].errorbar(
        centers,
        pass_ratio,
        yerr=pass_ratio_err,
        fmt=".",
        capsize=2,
    )
    axes[0, 1].set_title("Matched data PASS / AAOGEN PASS")
    axes[0, 1].set_ylabel("Shape ratio")
    axes[0, 1].grid(alpha=0.25)

    axes[1, 0].step(
        centers,
        pi0_pass_norm,
        where="mid",
        label=f"AAOGEN PASS ({int(np.count_nonzero(pi0_pass)):,})",
    )
    axes[1, 0].step(
        centers,
        pi0_fail_norm,
        where="mid",
        label=f"AAOGEN FAIL ({int(np.count_nonzero(pi0_fail)):,})",
    )
    axes[1, 0].set_title("AAOGEN reconstruction-bias comparison")
    axes[1, 0].set_xlabel(variable.label)
    axes[1, 0].set_ylabel("Fraction / bin")
    axes[1, 0].legend(fontsize=9)
    axes[1, 0].grid(alpha=0.25)

    axes[1, 1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1, 1].errorbar(
        centers,
        mc_ratio,
        yerr=mc_ratio_err,
        fmt=".",
        capsize=2,
    )
    axes[1, 1].set_title("AAOGEN PASS / AAOGEN FAIL")
    axes[1, 1].set_xlabel(variable.label)
    axes[1, 1].set_ylabel("Shape ratio")
    axes[1, 1].grid(alpha=0.25)

    for axis in axes.flat:
        axis.set_xlim(low, high)
    # endfor

    fig.suptitle(
        f"{period}, {category_title(detector, sector)}: "
        r"$\theta_{\gamma\gamma}$ control diagnostics"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def process_period(period,args_dict):
    args = argparse.Namespace(**args_dict)
    # Import and validate once per worker, before any ROOT I/O.
    fitter_preflight = preflight_exclusivity_fitter(args)
    pdir = Path(args.output_dir) / period.key
    for directory in [
        pdir,
        pdir / "shape_diagnostics",
        pdir / "fail_fits",
        pdir / "fit_component_tables",
        pdir / "fit_diagnostic_summaries",
        pdir / "fit_diagnostic_json",
        pdir / "theta_control_diagnostics",
    ]:
        directory.mkdir(parents=True, exist_ok=True)
    # endfor
    data,_,dc=read_opportunities(period.epg_data,period.beam_energy_GeV,'epgamma data',args,False); dvcs,_,bc=read_opportunities(period.dvcs_mc,period.beam_energy_GeV,'DVCSGEN epgamma MC',args,True); pi0,_,pc=read_opportunities(period.pi0_epg_mc,period.beam_energy_GeV,'AAOGEN-as-epgamma MC',args,True)
    cuts={'data':dc,'dvcs_mc':bc,'pi0_mc':pc}; json.dump(cuts,open(pdir/'opportunity_cutflows.json','w'),indent=2); emit_cutflow_diagnostics(pdir,period.label,cuts)
    egd,egd_diag=read_epgg(period.epgg_data,period.beam_energy_GeV,'epgammagamma data',args); egm,egm_diag=read_epgg(period.pi0_epgg_mc,period.beam_energy_GeV,'AAOGEN epgammagamma MC',args)
    dm=match_truth_partners(data,egd,'data',args); pm=match_truth_partners(pi0,egm,'mc',args)
    plot_expected_probe_diagnostics(pdir/'expected_probe_diagnostics.png',period.label,{'Data':data,'DVCSGEN':dvcs,'AAOGEN':pi0}); plot_matching_residuals(pdir/'matching_residuals.png',period.label,dm,pm)
    rows = []
    meta = {}
    failures: Dict[str, object] = {}
    for det,sec in [('FT',0)]+[('FD',s) for s in range(1,7)]:
        stem=category_stem(det,sec); dcat=category_mask(data,det,sec); bcat=category_mask(dvcs,det,sec); pcat=category_mask(pi0,det,sec); dpass=dcat&dm.matched; dfail=dcat&~dm.matched; ppass=pcat&pm.matched; pfail=pcat&~pm.matched
        plot_shape_diagnostics(
            pdir / 'shape_diagnostics' / stem,
            period.label,
            det,
            sec,
            data,
            dvcs,
            pi0,
            dpass,
            dfail,
            ppass,
            pfail,
            bcat,
        )
        theta_control_plot = (
            pdir
            / "theta_control_diagnostics"
            / f"{stem}_theta_gamma_gamma_control.png"
        )
        plot_theta_control_diagnostics(
            theta_control_plot,
            period.label,
            det,
            sec,
            data,
            pi0,
            dpass,
            ppass,
            pfail,
        )
        try:
            fit = fit_fail_shapes(data, dvcs, pi0, dfail, bcat, pfail, args)
        except Exception as error:
            failure_traceback = traceback.format_exc()
            failure_payload = {
                "period": period.key,
                "period_label": period.label,
                "detector": det,
                "sector": sec,
                "category": category_title(det, sec),
                "exception_type": type(error).__name__,
                "exception": str(error),
                "traceback": failure_traceback,
                "resolved_fitter_path": fitter_preflight["resolved_path"],
                "fit_signature": fitter_preflight["fit_signature"],
                "counts": {
                    "data_pass": int(np.count_nonzero(dpass)),
                    "data_fail": int(np.count_nonzero(dfail)),
                    "pi0_mc_pass": int(np.count_nonzero(ppass)),
                    "pi0_mc_fail": int(np.count_nonzero(pfail)),
                    "dvcs_mc": int(np.count_nonzero(bcat)),
                },
                "theta_control_plot": str(theta_control_plot),
            }
            failures[stem] = failure_payload
            failure_path = (
                pdir / "fit_diagnostic_json" / f"{stem}_fit_failure.json"
            )
            with open(failure_path, "w", encoding="utf-8") as handle:
                json.dump(failure_payload, handle, indent=2)
            # endwith
            log(
                f"WARNING {period.label} {category_title(det, sec)} "
                f"fit failed: {type(error).__name__}: {error}. "
                f"Full traceback: {failure_path}"
            )
            continue
        # endtry
        plot_fail_fit(
            pdir / "fail_fits" / f"{stem}_fail_fit.png",
            period.label,
            det,
            sec,
            fit,
        )

        variable_plot_paths: Dict[str, str] = {}
        for key in fit.payload:
            variable_plot_path = (
                pdir
                / "fail_fits"
                / f"{stem}_{key}_stacked_residual.png"
            )
            plot_one_fail_fit_variable(
                variable_plot_path,
                period.label,
                det,
                sec,
                key,
                fit,
            )
            variable_plot_paths[key] = str(variable_plot_path)
        # endfor

        summary_plot_path = (
            pdir
            / "fit_diagnostic_summaries"
            / f"{stem}_deviance_and_fraction_summary.png"
        )
        plot_fit_diagnostic_summary(
            summary_plot_path,
            period.label,
            det,
            sec,
            fit,
        )

        component_table_paths = write_components(
            pdir / "fit_component_tables" / stem,
            fit,
        )
        npass=int(np.count_nonzero(dpass)); nfail=int(np.count_nonzero(dfail)); npi0fail=fit.fraction*nfail; npi0fail_stat=fit.stat_error*nfail; vals=[v['fraction'] for v in fit.per_variable.values()]; spread=max(abs(v-fit.fraction) for v in vals) if vals else math.nan; npi0fail_model=spread*nfail
        ed,edstat=efficiency_from_pass_fail(npass,npi0fail,npi0fail_stat); edmodel=(npass/((npass+npi0fail)**2))*npi0fail_model if (npass+npi0fail)>0 else math.nan; em,emerr=binomial_efficiency(int(np.count_nonzero(ppass)),int(np.count_nonzero(pcat))); sf,sfstat=ratio_errors(ed,edstat,em,emerr); sfmodel=abs(sf)*(edmodel/ed) if ed>0 and math.isfinite(edmodel) else math.nan
        row=FailOnlyEfficiencyRow(period.key,period.label,det,sec,npass,nfail,int(np.count_nonzero(ppass)),int(np.count_nonzero(pfail)),int(np.count_nonzero(bcat)),fit.fraction,fit.stat_error,spread,npi0fail,npi0fail_stat,npi0fail_model,ed,edstat,edmodel,em,emerr,sf,sfstat,sfmodel,fit.success,fit.message,fit.deviance,fit.ndf,fit.reduced_deviance,bool(fit.warnings),json.dumps({k:v['fraction'] for k,v in fit.per_variable.items()},sort_keys=True))
        rows.append(asdict(row))

        diagnostic_payload = fit_diagnostics_payload(
            period.label,
            det,
            sec,
            fit,
            npass,
            nfail,
            int(np.count_nonzero(ppass)),
            int(np.count_nonzero(pfail)),
            int(np.count_nonzero(bcat)),
        )
        diagnostic_payload["outputs"] = {
            "combined_fit_plot": str(
                pdir / "fail_fits" / f"{stem}_fail_fit.png"
            ),
            "per_variable_fit_plots": variable_plot_paths,
            "summary_plot": str(summary_plot_path),
            "component_tables": component_table_paths,
            "theta_control_plot": str(theta_control_plot),
        }

        diagnostic_json_path = (
            pdir / "fit_diagnostic_json" / f"{stem}_fit_diagnostics.json"
        )
        with open(
            diagnostic_json_path,
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(diagnostic_payload, handle, indent=2)
        # endwith

        meta[stem] = {
            "row": asdict(row),
            "per_variable": fit.per_variable,
            "warnings": fit.warnings,
            "diagnostic_json": str(diagnostic_json_path),
            "outputs": diagnostic_payload["outputs"],
        }
        log(
            f"{period.label} {category_title(det, sec)}: "
            f"PASS={npass:,}, FAIL={nfail:,}, "
            f"f_pi0_FAIL={fit.fraction:.4f}, "
            f"N_pi0_FAIL={npi0fail:.1f}, "
            f"eps_data={ed:.4f}, eps_MC={em:.4f}, "
            f"S={sf:.4f}, D/ndf={fit.reduced_deviance:.2f}"
        )
        for key, variable_result in fit.per_variable.items():
            log(
                f"  {period.label} {category_title(det, sec)} {key}: "
                f"individual f_pi0_FAIL="
                f"{variable_result['fraction']:.5f} +/- "
                f"{variable_result['stat_error']:.5f}; "
                f"deviance contribution at combined fraction="
                f"{variable_result['deviance_at_combined_fraction']:.2f}; "
                f"D/ndf="
                f"{variable_result['reduced_deviance_at_combined_fraction']:.2f}"
            )
        # endfor
    write_rows_csv(pdir / "fail_only_results.csv", rows)
    period_payload = {
        "rows": rows,
        "categories": meta,
        "fit_failures": failures,
        "fit_summary": {
            "successful_categories": len(rows),
            "failed_categories": len(failures),
            "expected_categories": 7,
        },
        "fitter_preflight": fitter_preflight,
        "matching": {"data": dm.summary, "mc": pm.summary},
        "epgg": {"data": egd_diag, "mc": egm_diag},
    }
    with open(pdir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(period_payload, handle, indent=2)
    # endwith

    if not rows:
        raise RuntimeError(
            f"{period.label}: all seven detector-category fits failed. "
            f"See {pdir / 'metadata.json'} and fit_diagnostic_json/*_fit_failure.json"
        )
    # endif
    return period.key, rows, period_payload

def write_rows_csv(path,rows):
    if not rows:return
    with open(path,'w',newline='') as h: w=csv.DictWriter(h,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def preflight(periods):
    out=[]
    for p in periods:
        d=asdict(p)
        for role in ['epg_data','epgg_data','dvcs_mc','pi0_epg_mc','pi0_epgg_mc']:
            n,b=require_tree(getattr(p,role)); d[role+'_entries']=n; log(f'Preflight {p.label} {role}: {n:,} entries')
        out.append(d)
    return out

def plot_all_period_results(path, rows):
    periods=[p.label for p in PERIODS]; cats=['FT']+[f'FD S{s}' for s in range(1,7)]; x=np.arange(7,dtype=float); offsets=np.linspace(-.24,.24,len(periods))
    fig,axes=plt.subplots(3,1,figsize=(14,14),sharex=True)
    for off,label in zip(offsets,periods):
        selected=[r for r in rows if r['period_label']==label]; selected.sort(key=lambda r:(0 if r['detector']=='FT' else 1,r['sector']))
        if len(selected)!=7: continue
        axes[0].errorbar(x+off,[r['efficiency_data'] for r in selected],yerr=[r['efficiency_data_stat_err'] for r in selected],fmt='o',capsize=2,label=label)
        axes[1].errorbar(x+off,[r['efficiency_mc'] for r in selected],yerr=[r['efficiency_mc_err'] for r in selected],fmt='o',capsize=2,label=label)
        axes[2].errorbar(x+off,[r['scale_factor'] for r in selected],yerr=[r['scale_factor_stat_err'] for r in selected],fmt='o',capsize=2,label=label)
    axes[0].set_ylabel(r'$\epsilon_{\mathrm{data}}$'); axes[1].set_ylabel(r'$\epsilon_{\mathrm{MC}}$'); axes[2].set_ylabel(r'$S_\gamma$'); axes[2].axhline(1,ls='--')
    axes[0].set_ylim(0,1.05); axes[1].set_ylim(0,1.05); vals=[r['scale_factor'] for r in rows if math.isfinite(r['scale_factor'])]; axes[2].set_ylim(0,max(1.25,1.2*max(vals)) if vals else 1.25)
    for a in axes: a.grid(alpha=.25); a.legend(ncol=3,fontsize=8)
    axes[2].set_xticks(x,cats); fig.suptitle('FAIL-only photon efficiencies and data/MC scale factors'); fig.tight_layout(rect=(0,0,1,.96)); fig.savefig(path,dpi=180); plt.close(fig)

def main():
    args = parse_args()
    periods = selected_periods(args)
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Resolve/import/validate before the much more expensive ROOT preflight.
    fitter_preflight = preflight_exclusivity_fitter(args)
    manifest = {
        "args": vars(args),
        "fitter_preflight": fitter_preflight,
        "periods": preflight(periods),
    }
    with open(out / "input_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    workers = max(1, min(args.workers, MAX_WORKERS, len(periods)))
    rows: List[Dict[str, object]] = []
    metadata: Dict[str, object] = {}
    period_failures: Dict[str, object] = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(process_period, period, vars(args)): period
            for period in periods
        }
        for future in concurrent.futures.as_completed(futures):
            period = futures[future]
            try:
                key, period_rows, period_metadata = future.result()
            except Exception as error:
                failure_payload = {
                    "period": period.key,
                    "period_label": period.label,
                    "exception_type": type(error).__name__,
                    "exception": str(error),
                    "traceback": traceback.format_exc(),
                }
                period_failures[period.key] = failure_payload
                log(
                    f"WARNING period {period.label} failed: "
                    f"{type(error).__name__}: {error}"
                )
                continue
            # endtry
            rows.extend(period_rows)
            metadata[key] = period_metadata
            log(f"Completed {key}")
        # endfor
    # endwith

    order = {period.key: index for index, period in enumerate(PERIODS)}
    rows.sort(
        key=lambda row: (
            order[row["period"]],
            0 if row["detector"] == "FT" else 1,
            row["sector"],
        )
    )
    write_rows_csv(out / "photon_efficiency_fail_only_results.csv", rows)

    final_payload = {
        "rows": rows,
        "metadata": metadata,
        "period_failures": period_failures,
        "fitter_preflight": fitter_preflight,
        "summary": {
            "successful_category_rows": len(rows),
            "expected_category_rows": 7 * len(periods),
            "failed_periods": len(period_failures),
        },
    }
    with open(
        out / "photon_efficiency_fail_only_results.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(final_payload, handle, indent=2)
    # endwith

    if not rows:
        raise RuntimeError(
            "All detector-category fits failed; refusing to write an empty "
            "summary plot. Inspect photon_efficiency_fail_only_results.json "
            "and per-period fit failure JSON files."
        )
    # endif

    plot_all_period_results(
        out / "all_periods_integrated_efficiencies.png",
        rows,
    )
    log(f"Wrote outputs to {out}")
    return 0

if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception as e: print(f'FATAL ERROR: {e}',file=sys.stderr,flush=True); raise
