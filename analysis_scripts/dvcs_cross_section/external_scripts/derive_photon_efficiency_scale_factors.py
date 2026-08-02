#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v16_cutflow_csv_fix.py

Production extraction of the data/MC efficiency scale factor for reconstructing
the high-energy photon in exclusive pi0 tag-and-probe events.

This revision removes the three inappropriate tag-dependent DVCS cuts
(z(tag), -t1(tag), and electron-tag opening angle) from the nominal
opportunity selection. It retains the predicted-probe energy, massless
four-vector consistency, and FT/FD acceptance requirements. It also writes
a fully sequential opportunity cut audit before any template fit is attempted.

Poor reduced deviance is recorded and printed as a warning, but does not abort the extraction. Population inconsistencies such as N_found > N_expected also no longer abort diagnostic production: the official efficiency is set to NaN while raw ratios, component tables, pull plots, and yield-diagnosis plots are retained. Mathematical failures such as non-convergence, invalid parameters, boundary solutions, missing component fits, or unphysical yields remain fatal. Each cut is reported with its mathematical
definition, branch name, cumulative count, incremental survival fraction,
total survival fraction, rejection count, and FT/FD population where
applicable. The audit is written even if a later fit aborts.

The selected epgammaX denominator is decomposed in data with one validated
BH/DVCS + AAOGEN template fit per photon detector category. CD-FD and FD-FD
proton topologies are combined before fitting. Native epgammagamma records provide the found
probe. Their daughter photons are reconstructed with the validated Trento-basis
transformation, stored electron-photon opening angles, and four-vector closure
to the reconstructed pi0. The measured epgamma tag identifies one daughter;
the opposite daughter is the found high-energy probe.

Outputs are epsilon_data, epsilon_MC, and S_gamma = epsilon_data/epsilon_MC for
FT and FD sectors 1--6 in every RGA period.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import importlib.util
import json
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
    "Delta_phi",
    "theta_gamma_gamma",
    "pTmiss",
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


def locate_exclusivity_fitter(explicit_path: Optional[str]) -> Path:
    if explicit_path is not None:
        candidate = Path(explicit_path).expanduser().resolve()
        if not candidate.is_file():
            raise FileNotFoundError(f"Explicit exclusivity fitter not found: {candidate}")
        # endif
        return candidate

    script_dir = Path(__file__).resolve().parent
    cwd = Path.cwd().resolve()
    search_dirs = (
        script_dir,
        script_dir.parent,
        script_dir.parent / "external_scripts",
        cwd,
        cwd.parent,
        cwd / "external_scripts",
        cwd.parent / "external_scripts",
    )
    candidates: List[Path] = []
    for directory in search_dirs:
        candidates.append(directory / "plot_exclusivity_data_dvcs_pi0_mc.py")
    # endfor
    for directory in search_dirs:
        candidates.extend(sorted(directory.glob("plot_exclusivity_data_dvcs_pi0_mc(*).py"), reverse=True))
    # endfor

    source = next((candidate for candidate in candidates if candidate.is_file()), None)
    if source is None:
        raise FileNotFoundError(
            "Could not locate plot_exclusivity_data_dvcs_pi0_mc.py beside this "
            "script or one directory above it. Use --exclusivity-fitter."
        )
    # endif
    return source


def load_exclusivity_fitter(args: argparse.Namespace):
    global _EXCLUSIVITY_MODULE
    if _EXCLUSIVITY_MODULE is not None:
        return _EXCLUSIVITY_MODULE
    # endif

    source = locate_exclusivity_fitter(args.exclusivity_fitter)
    spec = importlib.util.spec_from_file_location("_photon_eff_exclusivity", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Unable to load exclusivity fitter from {source}")
    # endif
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    original_aliases = {
        variable.branch: tuple(getattr(variable, "aliases", ()))
        for variable in module.VARIABLES
    }
    module.VARIABLES = tuple(
        module.VariableConfig(
            "theta" if variable.key == "theta_cm" else variable.key,
            variable.label,
            variable.bins,
            variable.low,
            variable.high,
            aliases=original_aliases.get(
                "theta" if variable.key == "theta_cm" else variable.key,
                (),
            ),
        )
        for variable in FIT_VARIABLES
    )
    _EXCLUSIVITY_MODULE = module
    log(f"Loaded validated exclusivity fitter from {source}")
    return module


def mapped_histograms(histograms: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    return {
        "Delta_phi": np.asarray(histograms["Delta_phi"], dtype=np.float64),
        "theta": np.asarray(histograms["theta_cm"], dtype=np.float64),
        "theta_gamma_gamma": np.asarray(histograms["theta_gamma_gamma"], dtype=np.float64),
        "pTmiss": np.asarray(histograms["pTmiss"], dtype=np.float64),
        "Emiss2": np.asarray(histograms["Emiss2"], dtype=np.float64),
        "Mx2": np.asarray(histograms["Mx2"], dtype=np.float64),
        "Mx2_2": np.asarray(histograms["Mx2_2"], dtype=np.float64),
    }


def run_template_fit(
    data_histograms: Mapping[str, np.ndarray],
    dvcs_histograms: Mapping[str, np.ndarray],
    pi0_histograms: Mapping[str, np.ndarray],
    topology: str,
    args: argparse.Namespace,
):
    module = load_exclusivity_fitter(args)
    label, detector1, detector2 = TOPOLOGY_INFO[topology]
    topology_config = module.TopologyConfig(
        topology,
        label,
        detector1,
        detector2,
    )
    return module.fit_shared_two_templates(
        mapped_histograms(data_histograms),
        mapped_histograms(dvcs_histograms),
        mapped_histograms(pi0_histograms),
        topology_config,
        max_shift_bins=args.exclusivity_max_shift_bins,
        max_smear_bins=args.exclusivity_max_smear_bins,
        min_counts=args.fit_min_counts,
        fraction_variable_branches=tuple(
            "theta" if key == "theta_cm" else key
            for key in FRACTION_DRIVERS
        ),
        shift_prior_bins=args.exclusivity_shift_prior_bins,
        smear_prior_bins=args.exclusivity_smear_prior_bins,
        use_nuisance_penalties=not args.disable_exclusivity_nuisance_penalties,
        core_containment=args.exclusivity_dvcs_core_containment,
        fraction_containment=args.exclusivity_dvcs_fraction_containment,
        pi0_support_core_containment=args.exclusivity_pi0_core_containment,
        pi0_support_fraction_containment=args.exclusivity_pi0_fraction_containment,
        pi0_core_calibration=None,
        outside_overshoot_penalty_weight=args.exclusivity_outside_overshoot_penalty,
        emiss2_mean_order_penalty_weight=args.exclusivity_emiss2_mean_order_penalty,
    )


def plot_template_fit(
    path: Path,
    title: str,
    summary,
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    variable_results = summary.variable_results or {}
    for axis, variable in zip(axes.flat, FIT_VARIABLES):
        key = "theta" if variable.key == "theta_cm" else variable.key
        result = variable_results.get(key)
        if result is None or not result.success or result.fit_data_counts is None:
            axis.text(0.5, 0.5, "No valid fit", ha="center", va="center", transform=axis.transAxes)
            axis.set_axis_off()
            continue
        # endif
        edges = np.linspace(variable.low, variable.high, variable.bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        data = np.asarray(result.fit_data_counts)
        model = np.asarray(result.model_counts)
        dvcs = np.asarray(result.dvcs_component_counts)
        pi0 = np.asarray(result.pi0_component_counts)
        axis.errorbar(centers, data, yerr=np.sqrt(np.maximum(data, 1.0)), fmt="o", ms=2.5, label="Data")
        axis.step(centers, model, where="mid", linewidth=1.5, label="Total fit")
        axis.step(centers, dvcs, where="mid", linewidth=1.1, linestyle="--", label="BH/DVCS")
        axis.step(centers, pi0, where="mid", linewidth=1.1, linestyle=":", label=r"$\pi^0$")
        axis.set_xlabel(variable.label)
        axis.set_ylabel("Counts / bin")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    # endfor
    axes.flat[-1].axis("off")
    fig.suptitle(
        f"{title}\n"
        rf"$f_{{\pi^0}}={summary.f_pi0:.4f}\pm{summary.f_pi0_err:.4f}$, "
        rf"$D/\mathrm{{ndf}}={summary.deviance:.1f}/{summary.ndf}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def write_template_fit_components(
    output_dir: Path,
    category_stem: str,
    summary,
) -> Dict[str, str]:
    """
    Write one CSV per fit variable containing data, total model, BH/DVCS,
    pi0 component, residual, pull, and Poisson deviance contribution.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    written: Dict[str, str] = {}
    variable_results = getattr(summary, "variable_results", None) or {}

    for variable in FIT_VARIABLES:
        key = "theta" if variable.key == "theta_cm" else variable.key
        result = variable_results.get(key)
        if (
            result is None
            or getattr(result, "fit_data_counts", None) is None
            or getattr(result, "model_counts", None) is None
        ):
            continue
        # endif

        edges = np.linspace(variable.low, variable.high, variable.bins + 1)
        data = np.asarray(result.fit_data_counts, dtype=float)
        model = np.asarray(result.model_counts, dtype=float)
        dvcs = np.asarray(result.dvcs_component_counts, dtype=float)
        pi0 = np.asarray(result.pi0_component_counts, dtype=float)
        variance = np.maximum(model, 1.0)
        residual = data - model
        pull = residual / np.sqrt(variance)

        deviance_contribution = np.zeros_like(data, dtype=float)
        positive = (data > 0.0) & (model > 0.0)
        deviance_contribution[positive] = 2.0 * (
            model[positive] - data[positive]
            + data[positive] * np.log(data[positive] / model[positive])
        )
        zero_data = (data <= 0.0) & (model > 0.0)
        deviance_contribution[zero_data] = 2.0 * model[zero_data]

        path = output_dir / f"{category_stem}_{key}_components.csv"
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "bin_low",
                    "bin_high",
                    "bin_center",
                    "data",
                    "total_model",
                    "bh_dvcs_component",
                    "pi0_component",
                    "residual",
                    "pull",
                    "poisson_deviance_contribution",
                ]
            )
            for index in range(data.size):
                writer.writerow(
                    [
                        edges[index],
                        edges[index + 1],
                        0.5 * (edges[index] + edges[index + 1]),
                        data[index],
                        model[index],
                        dvcs[index],
                        pi0[index],
                        residual[index],
                        pull[index],
                        deviance_contribution[index],
                    ]
                )
            # endfor
        # endwith
        written[key] = str(path)
    # endfor
    return written


def plot_template_fit_pulls(
    path: Path,
    title: str,
    summary,
) -> None:
    """Plot residual pulls for every fitted variable."""
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    variable_results = getattr(summary, "variable_results", None) or {}

    for axis, variable in zip(axes.flat, FIT_VARIABLES):
        key = "theta" if variable.key == "theta_cm" else variable.key
        result = variable_results.get(key)
        if (
            result is None
            or getattr(result, "fit_data_counts", None) is None
            or getattr(result, "model_counts", None) is None
        ):
            axis.text(
                0.5,
                0.5,
                "No valid fit",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
            axis.set_axis_off()
            continue
        # endif

        edges = np.linspace(variable.low, variable.high, variable.bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        data = np.asarray(result.fit_data_counts, dtype=float)
        model = np.asarray(result.model_counts, dtype=float)
        pull = (data - model) / np.sqrt(np.maximum(model, 1.0))

        axis.axhline(0.0, linewidth=1.0)
        axis.axhline(3.0, linestyle="--", linewidth=0.8)
        axis.axhline(-3.0, linestyle="--", linewidth=0.8)
        axis.plot(centers, pull, marker="o", linestyle="none", markersize=2.5)
        axis.set_xlabel(variable.label)
        axis.set_ylabel(r"$(N_{\rm data}-N_{\rm fit})/\sqrt{N_{\rm fit}}$")
        axis.grid(alpha=0.25)
    # endfor

    axes.flat[-1].axis("off")
    reduced = (
        float(summary.deviance) / int(summary.ndf)
        if int(summary.ndf) > 0
        else math.inf
    )
    fig.suptitle(
        f"{title}\n"
        rf"$f_{{\pi^0}}={summary.f_pi0:.4f}\pm{summary.f_pi0_err:.4f}$, "
        rf"$D/\mathrm{{ndf}}={reduced:.3f}$"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_category_yield_diagnostics(
    path: Path,
    period_label: str,
    detector: str,
    sector: int,
    observed_data: float,
    expected_data: float,
    expected_data_err: float,
    found_data: float,
    expected_mc: float,
    found_mc: float,
    fit_details: Mapping[str, object],
    population_failures: Sequence[str],
) -> None:
    """Visualize the numerator/denominator bookkeeping and fit decomposition."""
    category_label = detector + (
        f" sector {sector}" if detector == "FD" else ""
    )
    fraction = float(fit_details.get("fraction_pi0", math.nan))
    fraction_err = float(fit_details.get("fraction_pi0_err", math.nan))
    reduced = float(fit_details.get("reduced_deviance", math.nan))
    expected_dvcs = observed_data - expected_data

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    labels_data = [
        "Observed\ndata",
        "Fitted\nBH/DVCS",
        "Fitted\npi0",
        "Found\npartner",
    ]
    values_data = [
        observed_data,
        expected_dvcs,
        expected_data,
        found_data,
    ]
    errors_data = [
        math.sqrt(max(observed_data, 1.0)),
        expected_data_err,
        expected_data_err,
        math.sqrt(max(found_data, 1.0)),
    ]
    x_data = np.arange(len(labels_data))
    axes[0, 0].bar(x_data, values_data, yerr=errors_data, capsize=3)
    axes[0, 0].set_xticks(x_data, labels_data)
    axes[0, 0].set_ylabel("Events")
    axes[0, 0].set_title("Data yield decomposition")
    axes[0, 0].grid(axis="y", alpha=0.25)

    labels_mc = ["Expected\nAAOGEN", "Found\npartner"]
    values_mc = [expected_mc, found_mc]
    x_mc = np.arange(len(labels_mc))
    axes[0, 1].bar(x_mc, values_mc)
    axes[0, 1].set_xticks(x_mc, labels_mc)
    axes[0, 1].set_ylabel("Events")
    axes[0, 1].set_title("AAOGEN reconstruction counts")
    axes[0, 1].grid(axis="y", alpha=0.25)

    raw_data_eff = (
        found_data / expected_data if expected_data > 0.0 else math.nan
    )
    raw_mc_eff = found_mc / expected_mc if expected_mc > 0.0 else math.nan
    ratios = [fraction, raw_data_eff, raw_mc_eff]
    ratio_labels = [
        r"Fit $f_{\pi^0}$",
        r"Raw $\epsilon_{\rm data}$",
        r"$\epsilon_{\rm MC}$",
    ]
    x_ratio = np.arange(len(ratios))
    axes[1, 0].bar(x_ratio, ratios)
    axes[1, 0].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1, 0].set_xticks(x_ratio, ratio_labels)
    axes[1, 0].set_ylabel("Fraction")
    finite_ratios = [value for value in ratios if math.isfinite(value)]
    ymax = max([1.1] + [1.15 * value for value in finite_ratios])
    axes[1, 0].set_ylim(0.0, ymax)
    axes[1, 0].grid(axis="y", alpha=0.25)

    axes[1, 1].axis("off")
    warning_lines = list(fit_details.get("warning_reasons", []) or [])
    status_lines = [
        f"Observed data entries: {observed_data:,.0f}",
        f"Fitted pi0 yield: {expected_data:,.1f} +/- {expected_data_err:,.1f}",
        f"Found data partners: {found_data:,.0f}",
        f"Expected AAOGEN opportunities: {expected_mc:,.0f}",
        f"Found AAOGEN partners: {found_mc:,.0f}",
        "",
        f"f_pi0 = {fraction:.6f} +/- {fraction_err:.6f}",
        f"Deviance / ndf = {reduced:.4f}",
        f"Raw data efficiency = {raw_data_eff:.6f}",
        f"MC efficiency = {raw_mc_eff:.6f}",
    ]
    if population_failures:
        status_lines.extend(
            ["", "POPULATION INVALID:"] + [f"- {item}" for item in population_failures]
        )
    else:
        status_lines.extend(["", "Population consistency: PASS"])
    # endif
    if warning_lines:
        status_lines.extend(
            ["", "FIT WARNINGS:"] + [f"- {item}" for item in warning_lines]
        )
    # endif
    axes[1, 1].text(
        0.02,
        0.98,
        "\n".join(status_lines),
        ha="left",
        va="top",
        transform=axes[1, 1].transAxes,
        family="monospace",
        fontsize=9,
    )

    fig.suptitle(f"{period_label}: {category_label} yield diagnosis")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def fit_category(
    data: OpportunityRecords,
    dvcs: OpportunityRecords,
    pi0: OpportunityRecords,
    detector: str,
    sector: int,
    fit_dir: Path,
    period_label: str,
    args: argparse.Namespace,
) -> Tuple[float, float, float, float, int, bool, Dict[str, object]]:
    """
    Fit one photon-detector category after combining all supported proton
    topologies.

    The photon-efficiency correction depends on the predicted probe detector,
    not on whether the proton was reconstructed in the CD or FD.  Therefore:

      * FT contains every selected opportunity whose predicted probe is FT;
      * FD sector N contains every selected opportunity whose predicted probe
        is in FD sector N.

    CD-FD and FD-FD populations are merged before histogramming and before the
    single two-template fit.  No topology can silently disappear from the data
    denominator while remaining in the matched numerator.
    """
    data_mask = category_mask(data, detector, sector)
    dvcs_mask = category_mask(dvcs, detector, sector)
    pi0_mask = category_mask(pi0, detector, sector)

    n_data = int(np.count_nonzero(data_mask))
    n_dvcs = int(np.count_nonzero(dvcs_mask))
    n_pi0 = int(np.count_nonzero(pi0_mask))
    category_label = detector + (f"_sector_{sector}" if detector == "FD" else "")

    details: Dict[str, object] = {
        "fit_scope": "combined_proton_topologies",
        "photon_detector": detector,
        "photon_sector": int(sector),
        "data": n_data,
        "dvcs": n_dvcs,
        "pi0": n_pi0,
        "included_topology_codes": sorted(
            {
                int(code)
                for code in np.concatenate(
                    (
                        data.topology_code[data_mask],
                        dvcs.topology_code[dvcs_mask],
                        pi0.topology_code[pi0_mask],
                    )
                )
                if int(code) >= 0
            }
        ),
    }

    if min(n_data, n_dvcs, n_pi0) < args.fit_min_counts:
        raise RuntimeError(
            f"{period_label} {category_label}: combined detector-category fit "
            f"has insufficient support: data={n_data:,}, DVCSGEN={n_dvcs:,}, "
            f"AAOGEN={n_pi0:,}; each must be >= {args.fit_min_counts:,}. "
            "No efficiency or scale factor will be computed."
        )
    # endif

    # The validated fitter requires a TopologyConfig.  Its topology object is
    # used to configure/label the shared fit; the histograms supplied here are
    # already the combined detector-category population.  Use the matching
    # photon topology as the representative configuration.
    representative_topology = "CD_FT" if detector == "FT" else "CD_FD"
    summary = run_template_fit(
        histograms_for_mask(data, data_mask),
        histograms_for_mask(dvcs, dvcs_mask),
        histograms_for_mask(pi0, pi0_mask),
        representative_topology,
        args,
    )

    plot_path = fit_dir / (
        "ft_combined_proton_topologies.png"
        if detector == "FT"
        else f"fd_sector_{sector}_combined_proton_topologies.png"
    )
    category_stem = (
        "ft_combined_proton_topologies"
        if detector == "FT"
        else f"fd_sector_{sector}_combined_proton_topologies"
    )
    category_title = (
        f"{period_label}: predicted probe {detector}"
        + (f" sector {sector}" if detector == "FD" else "")
        + ", combined proton topologies"
    )
    plot_template_fit(plot_path, category_title, summary)

    pull_plot_path = fit_dir / f"{category_stem}_pulls.png"
    plot_template_fit_pulls(pull_plot_path, category_title, summary)

    component_dir = fit_dir / "fit_component_tables"
    component_tables = write_template_fit_components(
        component_dir,
        category_stem,
        summary,
    )

    fraction = float(getattr(summary, "f_pi0", math.nan))
    fraction_err = float(getattr(summary, "f_pi0_err", math.nan))
    deviance = float(getattr(summary, "deviance", math.nan))
    ndf = int(getattr(summary, "ndf", 0))
    reduced_deviance = deviance / ndf if ndf > 0 else math.inf
    summary_success = bool(getattr(summary, "success", False))
    message = str(getattr(summary, "message", ""))

    details.update(
        {
            "success": summary_success,
            "message": message,
            "fraction_pi0": fraction,
            "fraction_pi0_err": fraction_err,
            "deviance": deviance,
            "ndf": ndf,
            "reduced_deviance": reduced_deviance,
            "representative_fitter_topology": representative_topology,
            "fit_plot": str(plot_path),
            "pull_plot": str(pull_plot_path),
            "component_tables": component_tables,
        }
    )

    failure_reasons: List[str] = []
    if not summary_success:
        failure_reasons.append(f"fitter reported failure: {message}")
    # endif
    if not math.isfinite(fraction):
        failure_reasons.append("fitted pi0 fraction is non-finite")
    # endif
    if not math.isfinite(fraction_err) or fraction_err < 0.0:
        failure_reasons.append("fitted pi0-fraction uncertainty is invalid")
    # endif
    if not math.isfinite(deviance) or ndf <= 0:
        failure_reasons.append(
            f"invalid fit deviance/ndf: deviance={deviance}, ndf={ndf}"
        )
    # endif

    margin = max(float(args.fit_fraction_boundary_margin), 0.0)
    if math.isfinite(fraction) and (
        fraction <= margin or fraction >= 1.0 - margin
    ):
        failure_reasons.append(
            f"fitted pi0 fraction {fraction:.8g} is within {margin:g} "
            "of a physical boundary"
        )
    # endif

    warning_reasons: List[str] = []
    max_reduced_deviance = float(args.fit_max_reduced_deviance)
    if (
        max_reduced_deviance > 0.0
        and math.isfinite(reduced_deviance)
        and reduced_deviance > max_reduced_deviance
    ):
        warning_reasons.append(
            f"deviance/ndf={reduced_deviance:.4g} exceeds "
            f"{max_reduced_deviance:g}"
        )
    # endif

    variable_results = getattr(summary, "variable_results", None) or {}
    invalid_variables = []
    for variable in FIT_VARIABLES:
        key = "theta" if variable.key == "theta_cm" else variable.key
        result = variable_results.get(key)
        if result is None or not bool(getattr(result, "success", False)):
            invalid_variables.append(key)
        # endif
    # endfor
    if invalid_variables:
        failure_reasons.append(
            "invalid component fits for: " + ", ".join(invalid_variables)
        )
    # endif

    details["fit_quality_warning"] = bool(warning_reasons)
    details["warning_reasons"] = warning_reasons
    details["passes_reduced_deviance_threshold"] = not bool(warning_reasons)
    details["fit_max_reduced_deviance_threshold"] = max_reduced_deviance

    if warning_reasons:
        log(
            f"WARNING: {period_label} {category_label}: combined template fit "
            f"retained despite poor goodness of fit: "
            + "; ".join(warning_reasons)
            + f". Inspect diagnostic plot: {plot_path}"
        )
    # endif

    if failure_reasons:
        details["failure_reasons"] = failure_reasons
        raise RuntimeError(
            f"{period_label} {category_label}: combined template fit rejected. "
            + "; ".join(failure_reasons)
            + f". Diagnostic plot: {plot_path}"
        )
    # endif

    expected = fraction * n_data
    # Treat the fitted fraction uncertainty and the finite data count as
    # independent contributions, consistent with the previous implementation.
    expected_err = math.sqrt(
        (n_data * fraction_err) ** 2 + fraction**2 * n_data
    )

    if not math.isfinite(expected) or expected <= 0.0:
        raise RuntimeError(
            f"{period_label} {category_label}: fitted expected data yield is "
            f"invalid ({expected})."
        )
    # endif
    if not math.isfinite(expected_err) or expected_err < 0.0:
        raise RuntimeError(
            f"{period_label} {category_label}: fitted expected-yield "
            f"uncertainty is invalid ({expected_err})."
        )
    # endif

    details["expected_data"] = expected
    details["expected_data_err"] = expected_err
    return (
        expected,
        expected_err,
        fraction,
        deviance,
        ndf,
        True,
        details,
    )



def efficiency_error_ratio(
    numerator: float,
    numerator_err: float,
    denominator: float,
    denominator_err: float,
) -> Tuple[float, float]:
    if denominator <= 0.0:
        return math.nan, math.nan
    # endif
    efficiency = numerator / denominator
    error = math.sqrt(
        (numerator_err / denominator)**2
        + (numerator * denominator_err / denominator**2)**2
    )
    return efficiency, error


def binomial_efficiency(found: float, expected: float) -> Tuple[float, float]:
    if expected <= 0.0:
        return math.nan, math.nan
    # endif
    efficiency = found / expected
    variance = max(efficiency * (1.0 - efficiency) / expected, 0.0)
    return efficiency, math.sqrt(variance)


def scale_factor(
    data_efficiency: float,
    data_error: float,
    mc_efficiency: float,
    mc_error: float,
) -> Tuple[float, float]:
    if (
        not math.isfinite(data_efficiency)
        or not math.isfinite(mc_efficiency)
        or mc_efficiency <= 0.0
    ):
        return math.nan, math.nan
    # endif
    value = data_efficiency / mc_efficiency
    relative_variance = 0.0
    if data_efficiency > 0.0 and math.isfinite(data_error):
        relative_variance += (data_error / data_efficiency)**2
    # endif
    if math.isfinite(mc_error):
        relative_variance += (mc_error / mc_efficiency)**2
    # endif
    return value, abs(value) * math.sqrt(relative_variance)


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



def plot_period_efficiencies(
    path: Path,
    period_label: str,
    rows: Sequence[EfficiencyRow],
) -> None:
    labels = ["FT"] + [f"FD S{sector}" for sector in range(1, 7)]
    x = np.arange(len(labels), dtype=float)
    data_eff = np.asarray([row.efficiency_data for row in rows])
    data_err = np.asarray([row.efficiency_data_err for row in rows])
    mc_eff = np.asarray([row.efficiency_mc for row in rows])
    mc_err = np.asarray([row.efficiency_mc_err for row in rows])
    scale = np.asarray([row.scale_factor for row in rows])
    scale_err = np.asarray([row.scale_factor_err for row in rows])

    fig, axes = plt.subplots(2, 1, figsize=(11, 9), sharex=True)
    axes[0].errorbar(x - 0.08, data_eff, yerr=data_err, fmt="o", label="Data (valid)")
    axes[0].errorbar(x + 0.08, mc_eff, yerr=mc_err, fmt="s", label="AAOGEN MC")

    invalid_mask = np.asarray([not row.population_valid for row in rows], dtype=bool)
    raw_invalid = np.asarray([row.raw_efficiency_data for row in rows], dtype=float)
    axes[0].plot(
        x[invalid_mask] - 0.08,
        raw_invalid[invalid_mask],
        linestyle="none",
        marker="x",
        markersize=8,
        label="Data raw ratio (invalid)",
    )
    axes[0].set_ylabel("Photon efficiency")
    finite_values = np.concatenate(
        (
            data_eff[np.isfinite(data_eff)],
            mc_eff[np.isfinite(mc_eff)],
            raw_invalid[np.isfinite(raw_invalid)],
        )
    )
    axes[0].set_ylim(
        0.0,
        max(1.05, 1.15 * float(np.max(finite_values)))
        if finite_values.size
        else 1.05,
    )
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].errorbar(x, scale, yerr=scale_err, fmt="o")
    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1].set_ylabel(r"$\epsilon_{\rm data}/\epsilon_{\rm MC}$")
    axes[1].set_xticks(x, labels)
    axes[1].grid(alpha=0.25)

    fig.suptitle(f"{period_label}: integrated photon-efficiency results")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def validate_efficiency_population(
    period_label: str,
    detector: str,
    sector: int,
    expected_data: float,
    found_data: float,
    expected_mc: float,
    found_mc: float,
    fit_success: bool,
    args: argparse.Namespace,
) -> List[str]:
    """
    Return population-consistency failures without aborting.

    An invalid population is retained for diagnosis, but its official
    efficiency and scale factor are set to NaN. Raw ratios are still recorded.
    """
    category_label = detector + (
        f" sector {sector}" if detector == "FD" else ""
    )
    tolerance = max(float(args.yield_consistency_tolerance), 0.0)

    failures: List[str] = []
    if not fit_success:
        failures.append("combined data template fit was not successful")
    # endif
    if not math.isfinite(expected_data) or expected_data <= 0.0:
        failures.append(f"expected data yield is invalid ({expected_data})")
    # endif
    if not math.isfinite(found_data) or found_data < 0.0:
        failures.append(f"found data yield is invalid ({found_data})")
    # endif
    if not math.isfinite(expected_mc) or expected_mc <= 0.0:
        failures.append(f"expected MC count is invalid ({expected_mc})")
    # endif
    if not math.isfinite(found_mc) or found_mc < 0.0:
        failures.append(f"found MC count is invalid ({found_mc})")
    # endif
    if (
        math.isfinite(expected_data)
        and math.isfinite(found_data)
        and found_data > expected_data + tolerance
    ):
        failures.append(
            f"found data count {found_data:.6g} exceeds fitted expected "
            f"pi0 yield {expected_data:.6g}"
        )
    # endif
    if (
        math.isfinite(expected_mc)
        and math.isfinite(found_mc)
        and found_mc > expected_mc + tolerance
    ):
        failures.append(
            f"found MC count {found_mc:.6g} exceeds expected MC count "
            f"{expected_mc:.6g}"
        )
    # endif

    if failures:
        log(
            f"WARNING: {period_label} {category_label}: population consistency "
            f"failed; official efficiency and scale factor will be NaN, but all "
            f"raw counts, ratios, and diagnostic plots will be retained: "
            + "; ".join(failures)
        )
    # endif
    return failures




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
        # endfor
    # endwith



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
        # endfor
    # endfor


def process_period(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    fit_dir = period_dir / "template_fits"
    category_diagnostic_dir = period_dir / "category_diagnostics"
    period_dir.mkdir(parents=True, exist_ok=True)
    fit_dir.mkdir(parents=True, exist_ok=True)
    category_diagnostic_dir.mkdir(parents=True, exist_ok=True)

    data, _, data_cutflow = read_opportunities(
        period.epg_data, period.beam_energy_GeV, "epgamma data", args, deduplicate=False
    )
    dvcs, _, dvcs_cutflow = read_opportunities(
        period.dvcs_mc, period.beam_energy_GeV, "DVCSGEN epgamma MC", args, deduplicate=True
    )
    pi0, _, pi0_cutflow = read_opportunities(
        period.pi0_epg_mc, period.beam_energy_GeV, "AAOGEN-as-epgamma MC", args, deduplicate=True
    )

    # Write the complete sequential audit before any native-tree reconstruction
    # or template fit. Therefore a later intentional fit failure cannot erase
    # the information needed to diagnose the low statistics.
    initial_cutflows = {
        "data": data_cutflow,
        "dvcs_mc": dvcs_cutflow,
        "pi0_mc": pi0_cutflow,
    }
    with open(
        period_dir / "opportunity_cutflows.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(initial_cutflows, handle, indent=2)
    # endwith
    emit_cutflow_diagnostics(period_dir, period.label, initial_cutflows)

    epgg_data, epgg_data_diag = read_epgg(
        period.epgg_data, period.beam_energy_GeV, "epgammagamma data", args
    )
    epgg_mc, epgg_mc_diag = read_epgg(
        period.pi0_epgg_mc, period.beam_energy_GeV, "AAOGEN epgammagamma MC", args
    )

    data_match = match_truth_partners(data, epgg_data, "data", args)
    mc_match = match_truth_partners(pi0, epgg_mc, "mc", args)

    for sample_label, result in (("data", data_match), ("AAOGEN MC", mc_match)):
        summary = result.summary
        log(
            f"{period.label} {sample_label}: identity={summary['identity_matches']:,}/"
            f"{summary['opportunities']:,}, valid|identity="
            f"{summary['conditional_valid_solution_given_identity']:.3%}, "
            f"tag|valid={summary['conditional_tag_match_given_valid_solution']:.3%}, "
            f"partner E>{args.found_probe_E_min:g} GeV="
            f"{summary['partner_above_threshold']:,}, matched="
            f"{summary['matched_opportunities']:,}"
        )
    # endfor

    plot_expected_probe_diagnostics(
        period_dir / "expected_probe_diagnostics.png",
        period.label,
        {"Data": data, "DVCSGEN": dvcs, "AAOGEN": pi0},
    )
    plot_matching_residuals(
        period_dir / "matching_residuals.png",
        period.label,
        data_match,
        mc_match,
    )

    rows: List[EfficiencyRow] = []
    fit_metadata: Dict[str, object] = {}
    categories = [("FT", 0)] + [("FD", sector) for sector in range(1, 7)]
    for detector, sector in categories:
        (
            expected_data,
            expected_data_err,
            fraction_pi0,
            fit_deviance,
            fit_ndf,
            fit_success,
            fit_details,
        ) = fit_category(
            data, dvcs, pi0,
            detector, sector,
            fit_dir, period.label, args,
        )

        data_category = category_mask(data, detector, sector)
        mc_category = category_mask(pi0, detector, sector)
        found_data = float(np.count_nonzero(data_category & data_match.matched))
        expected_mc = float(np.count_nonzero(mc_category))
        found_mc = float(np.count_nonzero(mc_category & mc_match.matched))

        population_failures = validate_efficiency_population(
            period.label,
            detector,
            sector,
            expected_data,
            found_data,
            expected_mc,
            found_mc,
            fit_success,
            args,
        )
        population_valid = not population_failures

        raw_efficiency_data = (
            found_data / expected_data
            if math.isfinite(expected_data) and expected_data > 0.0
            else math.nan
        )
        raw_efficiency_mc = (
            found_mc / expected_mc
            if math.isfinite(expected_mc) and expected_mc > 0.0
            else math.nan
        )
        raw_correction = (
            raw_efficiency_data / raw_efficiency_mc
            if (
                math.isfinite(raw_efficiency_data)
                and math.isfinite(raw_efficiency_mc)
                and raw_efficiency_mc > 0.0
            )
            else math.nan
        )

        if population_valid:
            efficiency_data, efficiency_data_err = efficiency_error_ratio(
                found_data,
                math.sqrt(max(found_data, 1.0)),
                expected_data,
                expected_data_err,
            )
            efficiency_mc, efficiency_mc_err = binomial_efficiency(
                found_mc,
                expected_mc,
            )
            correction, correction_err = scale_factor(
                efficiency_data,
                efficiency_data_err,
                efficiency_mc,
                efficiency_mc_err,
            )
        else:
            efficiency_data = math.nan
            efficiency_data_err = math.nan
            efficiency_mc = raw_efficiency_mc
            efficiency_mc_err = (
                binomial_efficiency(found_mc, expected_mc)[1]
                if math.isfinite(raw_efficiency_mc)
                and 0.0 <= raw_efficiency_mc <= 1.0
                else math.nan
            )
            correction = math.nan
            correction_err = math.nan
        # endif

        observed_data = float(np.count_nonzero(data_category))
        fit_details["population_valid"] = population_valid
        fit_details["population_failure_reasons"] = population_failures
        fit_details["observed_data"] = observed_data
        fit_details["found_data"] = found_data
        fit_details["expected_mc"] = expected_mc
        fit_details["found_mc"] = found_mc
        fit_details["raw_efficiency_data"] = raw_efficiency_data
        fit_details["raw_efficiency_mc"] = raw_efficiency_mc
        fit_details["raw_scale_factor"] = raw_correction

        category_stem = (
            "ft" if detector == "FT" else f"fd_sector_{sector}"
        )
        category_plot_path = (
            category_diagnostic_dir / f"{category_stem}_yield_diagnosis.png"
        )
        plot_category_yield_diagnostics(
            category_plot_path,
            period.label,
            detector,
            sector,
            observed_data,
            expected_data,
            expected_data_err,
            found_data,
            expected_mc,
            found_mc,
            fit_details,
            population_failures,
        )
        fit_details["yield_diagnostic_plot"] = str(category_plot_path)

        row = EfficiencyRow(
            period=period.key,
            period_label=period.label,
            detector=detector,
            sector=sector,
            observed_data=observed_data,
            expected_data=expected_data,
            expected_data_err=expected_data_err,
            found_data=found_data,
            raw_efficiency_data=raw_efficiency_data,
            efficiency_data=efficiency_data,
            efficiency_data_err=efficiency_data_err,
            expected_mc=expected_mc,
            found_mc=found_mc,
            raw_efficiency_mc=raw_efficiency_mc,
            efficiency_mc=efficiency_mc,
            efficiency_mc_err=efficiency_mc_err,
            raw_scale_factor=raw_correction,
            scale_factor=correction,
            scale_factor_err=correction_err,
            population_valid=population_valid,
            population_failure_reasons="; ".join(population_failures),
            fit_success=fit_success,
            fit_quality_warning=bool(
                fit_details.get("fit_quality_warning", False)
            ),
            fit_warning_reasons="; ".join(
                fit_details.get("warning_reasons", []) or []
            ),
            fit_fraction_pi0=fraction_pi0,
            fit_fraction_pi0_err=(
                expected_data_err / max(float(np.count_nonzero(data_category)), 1.0)
                if math.isfinite(expected_data_err)
                else math.nan
            ),
            fit_deviance=fit_deviance,
            fit_ndf=fit_ndf,
            fit_reduced_deviance=(
                fit_deviance / fit_ndf if fit_ndf > 0 else math.inf
            ),
        )
        rows.append(row)
        fit_metadata[f"{detector}_{sector}"] = fit_details

        status = "VALID" if population_valid else "INVALID-DIAGNOSTIC"
        log(
            f"{period.label} {detector}"
            + (f" sector {sector}" if detector == "FD" else "")
            + f" [{status}]: Nobs_data={observed_data:.0f}, "
            + f"Nexp_data={expected_data:.1f}, Nfound_data={found_data:.0f}, "
            + f"raw_eps_data={raw_efficiency_data:.4f}, "
            + f"Nexp_MC={expected_mc:.0f}, Nfound_MC={found_mc:.0f}, "
            + f"eps_MC={raw_efficiency_mc:.4f}, "
            + f"raw_S_gamma={raw_correction:.4f}"
        )
    # endfor

    write_rows_csv(
        period_dir / "category_fit_and_efficiency_diagnostics.csv",
        [asdict(row) for row in rows],
    )

    plot_period_efficiencies(
        period_dir / "integrated_efficiencies.png",
        period.label,
        rows,
    )

    metadata = {
        "period": asdict(period),
        "opportunity_cutflows": {
            "data": data_cutflow,
            "dvcs_mc": dvcs_cutflow,
            "pi0_mc": pi0_cutflow,
        },
        "epgammagamma_reconstruction": {
            "data": epgg_data_diag,
            "mc": epgg_mc_diag,
        },
        "matching": {
            "data": data_match.summary,
            "mc": mc_match.summary,
        },
        "template_fits": fit_metadata,
    }
    with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    # endwith
    with open(period_dir / "opportunity_cutflows.json", "w", encoding="utf-8") as handle:
        json.dump(metadata["opportunity_cutflows"], handle, indent=2)
    # endwith
    with open(period_dir / "matching_summary.json", "w", encoding="utf-8") as handle:
        json.dump(metadata["matching"], handle, indent=2)
    # endwith

    return period.key, [asdict(row) for row in rows], metadata


def write_rows_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    # endif
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def plot_all_periods(
    efficiency_path: Path,
    scale_path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    periods = [period.label for period in PERIODS]
    categories = [("FT", 0)] + [("FD", sector) for sector in range(1, 7)]
    labels = ["FT"] + [f"FD S{sector}" for sector in range(1, 7)]
    x = np.arange(len(categories), dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(13, 10), sharex=True)
    offsets = np.linspace(-0.24, 0.24, len(periods))
    for offset, period_label in zip(offsets, periods):
        selected = [
            next(
                (
                    row for row in rows
                    if row["period_label"] == period_label
                    and row["detector"] == detector
                    and int(row["sector"]) == sector
                ),
                None,
            )
            for detector, sector in categories
        ]
        data = np.asarray([
            row["efficiency_data"] if row is not None else math.nan
            for row in selected
        ])
        data_err = np.asarray([
            row["efficiency_data_err"] if row is not None else math.nan
            for row in selected
        ])
        mc = np.asarray([
            row["efficiency_mc"] if row is not None else math.nan
            for row in selected
        ])
        mc_err = np.asarray([
            row["efficiency_mc_err"] if row is not None else math.nan
            for row in selected
        ])
        axes[0].errorbar(x + offset, data, yerr=data_err, fmt="o", ms=4, label=period_label)
        axes[1].errorbar(x + offset, mc, yerr=mc_err, fmt="s", ms=4, label=period_label)
    # endfor
    axes[0].set_ylabel(r"$\epsilon_{\rm data}$")
    axes[1].set_ylabel(r"$\epsilon_{\rm MC}$")
    for axis in axes:
        axis.set_ylim(0.0, 1.05)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8, ncol=3)
    # endfor
    axes[1].set_xticks(x, labels)
    fig.suptitle("Integrated photon efficiencies by period")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(efficiency_path, dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(13, 6))
    for offset, period_label in zip(offsets, periods):
        selected = [
            next(
                (
                    row for row in rows
                    if row["period_label"] == period_label
                    and row["detector"] == detector
                    and int(row["sector"]) == sector
                ),
                None,
            )
            for detector, sector in categories
        ]
        values = np.asarray([
            row["scale_factor"] if row is not None else math.nan
            for row in selected
        ])
        errors = np.asarray([
            row["scale_factor_err"] if row is not None else math.nan
            for row in selected
        ])
        axis.errorbar(x + offset, values, yerr=errors, fmt="o", ms=4, label=period_label)
    # endfor
    axis.axhline(1.0, linestyle="--", linewidth=1.0)
    axis.set_xticks(x, labels)
    axis.set_ylabel(r"$S_\gamma=\epsilon_{\rm data}/\epsilon_{\rm MC}$")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8, ncol=3)
    axis.set_title("Integrated photon-efficiency scale factors")
    fig.tight_layout()
    fig.savefig(scale_path, dpi=180)
    plt.close(fig)


def preflight(periods: Sequence[PeriodConfig], args: argparse.Namespace) -> Dict[str, object]:
    manifest: Dict[str, object] = {
        "script": Path(__file__).name,
        "created_unix_time": time.time(),
        "arguments": vars(args),
        "periods": [],
    }
    for period in periods:
        period_payload = asdict(period)
        for role in (
            "epg_data", "epgg_data", "dvcs_mc", "pi0_epg_mc", "pi0_epgg_mc"
        ):
            path = getattr(period, role)
            entries, keys = require_tree(path)
            period_payload[f"{role}_entries"] = entries
            period_payload[f"{role}_branches"] = keys
            log(f"Preflight {period.label} {role}: {entries:,} entries")
        # endfor
        manifest["periods"].append(period_payload)
    # endfor

    # Fail early if the fitter cannot be loaded.
    load_exclusivity_fitter(args)
    return manifest


def validate_runtime_helpers() -> None:
    """Fail immediately if a required output helper was accidentally omitted."""
    required_helpers = {
        "write_cutflow_csv": write_cutflow_csv,
        "write_cutflow_text": write_cutflow_text,
        "plot_cutflow": plot_cutflow,
        "emit_cutflow_diagnostics": emit_cutflow_diagnostics,
        "write_rows_csv": write_rows_csv,
    }
    missing = [
        name
        for name, function in required_helpers.items()
        if not callable(function)
    ]
    if missing:
        raise RuntimeError(
            "Required runtime helper(s) are unavailable: "
            + ", ".join(missing)
        )
    # endif


def main() -> int:
    validate_runtime_helpers()
    args = parse_args()
    periods = selected_periods(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = preflight(periods, args)
    with open(output_dir / "input_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    workers = max(1, min(int(args.workers), MAX_WORKERS, len(periods)))
    log(
        f"Selected {len(periods)} period(s); using {workers} worker(s). "
        f"Directed opportunity: E_tag >= {args.tag_E_min:g} GeV and "
        f"E_probe >= {args.probe_E_min:g} GeV."
    )

    rows: List[Dict[str, object]] = []
    metadata: Dict[str, object] = {}
    args_dict = vars(args)

    if workers == 1:
        for period in periods:
            key, period_rows, payload = process_period(period, args_dict)
            rows.extend(period_rows)
            metadata[key] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(process_period, period, args_dict): period.key
                for period in periods
            }
            for future in concurrent.futures.as_completed(futures):
                key, period_rows, payload = future.result()
                rows.extend(period_rows)
                metadata[key] = payload
                log(f"Completed {key}.")
            # endfor
        # endwith
    # endif

    period_order = {period.key: index for index, period in enumerate(PERIODS)}
    rows.sort(
        key=lambda row: (
            period_order[row["period"]],
            0 if row["detector"] == "FT" else 1,
            int(row["sector"]),
        )
    )

    write_rows_csv(output_dir / "photon_efficiency_scale_factors.csv", rows)
    with open(
        output_dir / "photon_efficiency_scale_factors.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "rows": rows,
                "period_metadata": metadata,
            },
            handle,
            indent=2,
        )
    # endwith

    plot_all_periods(
        output_dir / "all_periods_integrated_efficiencies.png",
        output_dir / "all_periods_integrated_scale_factors.png",
        rows,
    )

    log(f"Wrote photon-efficiency outputs to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
