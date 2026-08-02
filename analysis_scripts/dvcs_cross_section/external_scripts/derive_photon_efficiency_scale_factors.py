#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors.py

RGA high-energy photon-efficiency extraction using the newly reprocessed
multi-candidate epgamma trees.

Each epgamma entry contains one real reconstructed photon candidate. Entries
sharing an event identity therefore provide the actual reconstructed photon
list for that event:

  * data identity: exact (runnum, evnum);
  * AAOGEN identity: reconstructed electron/proton kinematic signature.

For every directed epgamma tag opportunity, the script calculates

    p_X = k + p_target - k' - p' - p_gamma_tag

and requires the predicted X to be a photon with E_X >= 2 GeV in FT or FD
acceptance. No upper bound is imposed on E_tag.

The found-probe decision no longer attempts to reconstruct daughter photons
from the pi0 tree. Instead it searches the other real epgamma entries in the
same event/signature group for a photon matching X. The epgammagamma tree is
used only as an event/signature-level confirmation that the event entered the
reconstructed-pi0 sample. The matched tag-partner pair is also required to lie
inside the configured pi0 invariant-mass window.

Data denominator:
  BH/DVCS + AAOGEN two-template fit to selected epgammaX opportunities.

Data numerator:
  selected opportunities with a matching partner photon in the same data
  event and with that event present in the epgammagamma data tree.

MC denominator:
  all selected AAOGEN epgammaX opportunities.

MC numerator:
  selected AAOGEN opportunities with a matching partner photon in the same
  reconstructed e/p signature group and with that signature present in the
  AAOGEN epgammagamma tree.

The final correction is S_gamma = epsilon_data / epsilon_MC for FT and FD
sectors 1--6.
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
class PartnerMatchResult:
    matched: np.ndarray
    catalog_confirmed: np.ndarray
    partner_index: np.ndarray
    probe_angle_deg: np.ndarray
    probe_relative_E: np.ndarray
    pair_mass_GeV: np.ndarray
    score: np.ndarray
    summary: Dict[str, object]


@dataclass
class EfficiencyRow:
    period: str
    period_label: str
    detector: str
    sector: int
    expected_data: float
    expected_data_err: float
    found_data: float
    efficiency_data: float
    efficiency_data_err: float
    expected_mc: float
    found_mc: float
    efficiency_mc: float
    efficiency_mc_err: float
    scale_factor: float
    scale_factor_err: float
    fit_success: bool
    fit_fraction_pi0: float
    fit_fraction_pi0_err: float
    fit_deviance: float
    fit_ndf: int


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

    parser.add_argument("--tag-match-angle-max-deg", type=float, default=1.0)
    parser.add_argument("--tag-match-relative-E-max", type=float, default=0.15)
    parser.add_argument("--probe-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument("--probe-match-relative-E-max", type=float, default=0.35)
    parser.add_argument("--pi0-mass-min", type=float, default=0.10)
    parser.add_argument("--pi0-mass-max", type=float, default=0.17)
    parser.add_argument("--mc-signature-decimals", type=int, default=10)


    parser.add_argument("--fit-min-counts", type=int, default=100)
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

    if resolved.get("t1") is not None:
        t1 = finite_array(arrays, resolved["t1"])
        mask &= np.isfinite(t1) & ((-t1) < args.minus_t_max)
    # endif

    if resolved.get("open_angle_ep2") is not None:
        opening = finite_array(arrays, resolved["open_angle_ep2"])
        mask &= np.isfinite(opening) & (opening > args.open_angle_min_deg)
    # endif

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
    cutflow = {
        "role": role,
        "path": path,
        "tree_entries": 0,
        "finite_kinematics": 0,
        "global_cuts": 0,
        "tag_energy": 0,
        "probe_finite": 0,
        "probe_energy": 0,
        "probe_photon_hypothesis": 0,
        "probe_acceptance": 0,
        "supported_topology": 0,
        "selected_before_deduplication": 0,
        "selected_after_deduplication": 0,
    }

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
        cutflow["tree_entries"] += int(n)

        def arr(key: str, default: float = math.nan, dtype=np.float64) -> np.ndarray:
            return finite_array(arrays, resolved.get(key), default=default, dtype=dtype)

        e_p, e_theta, e_phi = arr("e_p"), arr("e_theta"), arr("e_phi")
        p_p, p_theta, p_phi = arr("p_p"), arr("p_theta"), arr("p_phi")
        tag_E, tag_theta, tag_phi = arr("tag_E"), arr("tag_theta"), np.mod(arr("tag_phi"), 2.0 * math.pi)
        tag_detector = arr("tag_detector", default=-1, dtype=np.int16)
        proton_detector = arr("proton_detector", default=-1, dtype=np.int16)

        finite = (
            np.isfinite(e_p) & np.isfinite(e_theta) & np.isfinite(e_phi)
            & np.isfinite(p_p) & np.isfinite(p_theta) & np.isfinite(p_phi)
            & np.isfinite(tag_E) & np.isfinite(tag_theta) & np.isfinite(tag_phi)
            & (e_p > 0.0) & (p_p > 0.0) & (tag_E > 0.0)
        )
        cutflow["finite_kinematics"] += int(np.count_nonzero(finite))

        global_mask = finite & common_global_mask(arrays, resolved, args)
        cutflow["global_cuts"] += int(np.count_nonzero(global_mask))

        tag_mask = global_mask & (tag_E >= args.tag_E_min)
        cutflow["tag_energy"] += int(np.count_nonzero(tag_mask))

        # The partner search must use every reconstructed photon candidate in
        # the event, not only candidates whose own missing-X hypothesis passes
        # the expected-probe selection.
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
        append_candidate_store(candidate_store, candidate_values, tag_mask)

        probe = reconstruct_probe(
            beam_energy,
            e_p, e_theta, e_phi,
            p_p, p_theta, p_phi,
            tag_E, tag_theta, tag_phi,
        )
        probe_finite = (
            tag_mask
            & np.isfinite(probe["E"]) & np.isfinite(probe["p"])
            & np.isfinite(probe["theta"]) & np.isfinite(probe["phi"])
            & np.isfinite(probe["m2"]) & np.isfinite(probe["E_minus_p"])
            & (probe["E"] > 0.0) & (probe["p"] > 0.0)
        )
        cutflow["probe_finite"] += int(np.count_nonzero(probe_finite))

        probe_energy = (
            probe_finite
            & (probe["E"] >= args.probe_E_min)
            & (probe["E"] < args.probe_E_max)
        )
        cutflow["probe_energy"] += int(np.count_nonzero(probe_energy))

        photon_like = (
            probe_energy
            & (np.abs(probe["m2"]) < args.probe_m2_abs_max)
            & (np.abs(probe["E_minus_p"]) < args.probe_E_minus_p_abs_max)
        )
        cutflow["probe_photon_hypothesis"] += int(np.count_nonzero(photon_like))

        probe_theta_deg = np.degrees(probe["theta"])
        probe_detector, probe_sector = classify_probe(
            probe_theta_deg, probe["phi"], args
        )
        accepted = photon_like & (probe_detector >= 0)
        cutflow["probe_acceptance"] += int(np.count_nonzero(accepted))

        topo = topology_code(proton_detector, tag_detector)
        selected = accepted & (topo >= 0)
        cutflow["supported_topology"] += int(np.count_nonzero(selected))
        cutflow["selected_before_deduplication"] += int(np.count_nonzero(selected))

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
    cutflow["selected_after_deduplication"] = records.size()
    cutflow["branch_mapping"] = resolved
    cutflow["partner_candidate_records"] = candidates.size()
    return records, candidates, cutflow


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


def match_partner_photons(
    opportunities: OpportunityRecords,
    candidates: PhotonCandidateRecords,
    catalog: IdentityCatalog,
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

    opportunity_keys = identity_keys_data(opportunities) if mode == "data" else identity_keys_mc(
        opportunities, args.mc_signature_decimals
    )
    candidate_keys = identity_keys_data(candidates) if mode == "data" else identity_keys_mc(
        candidates, args.mc_signature_decimals
    )
    catalog_keys = identity_keys_data(catalog) if mode == "data" else identity_keys_mc(
        catalog, args.mc_signature_decimals
    )

    candidate_inverse, candidate_order, candidate_offsets = group_indices(candidate_keys)
    candidate_structured = structured_keys(candidate_keys)
    candidate_unique, candidate_first = np.unique(candidate_structured, return_index=True)
    candidate_group_by_key = {
        candidate_unique[index].tobytes(): int(candidate_inverse[candidate_first[index]])
        for index in range(candidate_unique.size)
    }
    catalog_key_set = {key.tobytes() for key in np.unique(structured_keys(catalog_keys))}

    shared_candidate_groups = 0
    catalog_groups = 0
    candidate_tests = 0
    exact_tag_candidates_skipped = 0
    mass_window_rejections = 0

    for opportunity_index, key in enumerate(structured_keys(opportunity_keys)):
        key_bytes = key.tobytes()
        if key_bytes not in catalog_key_set:
            continue
        # endif
        catalog_confirmed[opportunity_index] = True
        catalog_groups += 1
        group_index = candidate_group_by_key.get(key_bytes)
        if group_index is None:
            continue
        # endif
        shared_candidate_groups += 1
        members = candidate_order[
            candidate_offsets[group_index]:candidate_offsets[group_index + 1]
        ]
        best = None
        for candidate_index in members:
            candidate_tests += 1
            tag_angle = float(opening_angle_deg(
                np.asarray([opportunities.tag_theta[opportunity_index]]),
                np.asarray([opportunities.tag_phi[opportunity_index]]),
                np.asarray([candidates.photon_theta[candidate_index]]),
                np.asarray([candidates.photon_phi[candidate_index]]),
            )[0])
            tag_relE = abs(
                opportunities.tag_E[opportunity_index] - candidates.photon_E[candidate_index]
            ) / max(candidates.photon_E[candidate_index], 1.0e-12)
            if tag_angle < 1.0e-6 and tag_relE < 1.0e-10:
                exact_tag_candidates_skipped += 1
                continue
            # endif

            angle = float(opening_angle_deg(
                np.asarray([opportunities.probe_theta[opportunity_index]]),
                np.asarray([opportunities.probe_phi[opportunity_index]]),
                np.asarray([candidates.photon_theta[candidate_index]]),
                np.asarray([candidates.photon_phi[candidate_index]]),
            )[0])
            relative_E = abs(
                opportunities.probe_E[opportunity_index] - candidates.photon_E[candidate_index]
            ) / max(candidates.photon_E[candidate_index], 1.0e-12)
            if (
                candidates.photon_E[candidate_index] < args.probe_E_min
                or angle > args.probe_match_angle_max_deg
                or relative_E > args.probe_match_relative_E_max
            ):
                continue
            # endif

            mass = photon_pair_mass(
                float(opportunities.tag_E[opportunity_index]),
                float(opportunities.tag_theta[opportunity_index]),
                float(opportunities.tag_phi[opportunity_index]),
                float(candidates.photon_E[candidate_index]),
                float(candidates.photon_theta[candidate_index]),
                float(candidates.photon_phi[candidate_index]),
            )
            if not (args.pi0_mass_min <= mass <= args.pi0_mass_max):
                mass_window_rejections += 1
                continue
            # endif

            candidate_score = (
                (angle / max(args.probe_match_angle_max_deg, 1.0e-12)) ** 2
                + (relative_E / max(args.probe_match_relative_E_max, 1.0e-12)) ** 2
            )
            if best is None or candidate_score < best[0]:
                best = (candidate_score, int(candidate_index), angle, relative_E, mass)
            # endif
        # endfor

        if best is not None:
            matched[opportunity_index] = True
            score[opportunity_index] = best[0]
            partner_index[opportunity_index] = best[1]
            probe_angle[opportunity_index] = best[2]
            probe_relative_E[opportunity_index] = best[3]
            pair_mass[opportunity_index] = best[4]
        # endif
    # endfor

    return PartnerMatchResult(
        matched=matched,
        catalog_confirmed=catalog_confirmed,
        partner_index=partner_index,
        probe_angle_deg=probe_angle,
        probe_relative_E=probe_relative_E,
        pair_mass_GeV=pair_mass,
        score=score,
        summary={
            "mode": mode,
            "opportunities": n,
            "candidate_records": candidates.size(),
            "epgammagamma_catalog_records": catalog.size(),
            "catalog_confirmed_opportunities": int(np.count_nonzero(catalog_confirmed)),
            "catalog_confirmed_fraction": float(np.mean(catalog_confirmed)) if n else math.nan,
            "opportunities_with_candidate_group": shared_candidate_groups,
            "candidate_tests": candidate_tests,
            "exact_tag_candidates_skipped": exact_tag_candidates_skipped,
            "pi0_mass_window_rejections": mass_window_rejections,
            "matched_opportunities": int(np.count_nonzero(matched)),
            "matched_fraction": float(np.mean(matched)) if n else math.nan,
            "matched_fraction_given_catalog": (
                float(np.count_nonzero(matched) / np.count_nonzero(catalog_confirmed))
                if np.count_nonzero(catalog_confirmed) else math.nan
            ),
        },
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
    category_data = category_mask(data, detector, sector)
    category_dvcs = category_mask(dvcs, detector, sector)
    category_pi0 = category_mask(pi0, detector, sector)

    expected = 0.0
    expected_variance = 0.0
    deviance = 0.0
    ndf = 0
    fit_success = True
    weighted_fraction = 0.0
    weighted_total = 0.0
    details: Dict[str, object] = {}

    for code in (0, 1, 2):
        topology = topology_key(code)
        data_mask = category_data & (data.topology_code == code)
        dvcs_mask = category_dvcs & (dvcs.topology_code == code)
        pi0_mask = category_pi0 & (pi0.topology_code == code)

        n_data = int(np.count_nonzero(data_mask))
        n_dvcs = int(np.count_nonzero(dvcs_mask))
        n_pi0 = int(np.count_nonzero(pi0_mask))
        if min(n_data, n_dvcs, n_pi0) < args.fit_min_counts:
            details[topology] = {
                "success": False,
                "reason": "insufficient support",
                "data": n_data,
                "dvcs": n_dvcs,
                "pi0": n_pi0,
            }
            if n_data > 0:
                fit_success = False
            # endif
            continue
        # endif

        summary = run_template_fit(
            histograms_for_mask(data, data_mask),
            histograms_for_mask(dvcs, dvcs_mask),
            histograms_for_mask(pi0, pi0_mask),
            topology,
            args,
        )
        success = bool(summary.success) and math.isfinite(summary.f_pi0)
        details[topology] = {
            "success": success,
            "message": str(summary.message),
            "data": n_data,
            "dvcs": n_dvcs,
            "pi0": n_pi0,
            "fraction_pi0": float(summary.f_pi0),
            "fraction_pi0_err": float(summary.f_pi0_err),
            "deviance": float(summary.deviance),
            "ndf": int(summary.ndf),
        }
        if not success:
            fit_success = False
            continue
        # endif

        fraction = float(summary.f_pi0)
        fraction_err = (
            float(summary.f_pi0_err)
            if math.isfinite(summary.f_pi0_err)
            else math.sqrt(max(fraction * (1.0 - fraction) / max(n_data, 1), 0.0))
        )
        yield_pi0 = fraction * n_data
        yield_err = math.sqrt((n_data * fraction_err)**2 + fraction**2 * n_data)
        expected += yield_pi0
        expected_variance += yield_err**2
        weighted_fraction += fraction * n_data
        weighted_total += n_data
        deviance += float(summary.deviance)
        ndf += int(summary.ndf)

        plot_template_fit(
            fit_dir / f"{detector.lower()}_sector_{sector}_{topology}.png",
            f"{period_label}: predicted probe {detector}"
            + (f" sector {sector}" if detector == "FD" else "")
            + f", {TOPOLOGY_INFO[topology][0]}",
            summary,
        )
    # endfor

    combined_fraction = weighted_fraction / weighted_total if weighted_total > 0 else math.nan
    return (
        expected,
        math.sqrt(expected_variance),
        combined_fraction,
        deviance,
        ndf,
        fit_success,
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
    variables = (
        ("probe_angle_deg", np.linspace(0.0, 6.0, 81), "Probe angular residual (deg)"),
        ("probe_relative_E", np.linspace(0.0, 0.70, 81), "Probe relative energy residual"),
        ("pair_mass_GeV", np.linspace(0.08, 0.20, 81), r"Matched $M_{\gamma\gamma}$ (GeV)"),
        ("score", np.linspace(0.0, 2.0, 81), "Partner-match score"),
    )
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for axis, (name, bins, label) in zip(axes.flat, variables):
        for sample, result in (("Data", data_match), ("AAOGEN", mc_match)):
            values = np.asarray(getattr(result, name))
            values = values[result.matched & np.isfinite(values)]
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
        axis.legend()
    # endfor
    fig.suptitle(f"{period_label}: accepted partner-photon matches")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


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
    axes[0].errorbar(x - 0.08, data_eff, yerr=data_err, fmt="o", label="Data")
    axes[0].errorbar(x + 0.08, mc_eff, yerr=mc_err, fmt="s", label="AAOGEN MC")
    axes[0].set_ylabel("Photon efficiency")
    axes[0].set_ylim(0.0, 1.05)
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


def process_period(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    fit_dir = period_dir / "template_fits"
    period_dir.mkdir(parents=True, exist_ok=True)
    fit_dir.mkdir(parents=True, exist_ok=True)

    data, data_candidates, data_cutflow = read_opportunities(
        period.epg_data, period.beam_energy_GeV, "epgamma data", args, deduplicate=False
    )
    dvcs, _, dvcs_cutflow = read_opportunities(
        period.dvcs_mc, period.beam_energy_GeV, "DVCSGEN epgamma MC", args, deduplicate=True
    )
    pi0, pi0_candidates, pi0_cutflow = read_opportunities(
        period.pi0_epg_mc, period.beam_energy_GeV, "AAOGEN-as-epgamma MC", args, deduplicate=True
    )

    epgg_data, epgg_data_diag = read_identity_catalog(
        period.epgg_data, "epgammagamma data", "data", args
    )
    epgg_mc, epgg_mc_diag = read_identity_catalog(
        period.pi0_epgg_mc, "AAOGEN epgammagamma MC", "mc", args
    )

    data_match = match_partner_photons(
        data, data_candidates, epgg_data, "data", args
    )
    mc_match = match_partner_photons(
        pi0, pi0_candidates, epgg_mc, "mc", args
    )

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

        row = EfficiencyRow(
            period=period.key,
            period_label=period.label,
            detector=detector,
            sector=sector,
            expected_data=expected_data,
            expected_data_err=expected_data_err,
            found_data=found_data,
            efficiency_data=efficiency_data,
            efficiency_data_err=efficiency_data_err,
            expected_mc=expected_mc,
            found_mc=found_mc,
            efficiency_mc=efficiency_mc,
            efficiency_mc_err=efficiency_mc_err,
            scale_factor=correction,
            scale_factor_err=correction_err,
            fit_success=fit_success,
            fit_fraction_pi0=fraction_pi0,
            fit_fraction_pi0_err=(
                expected_data_err / max(float(np.count_nonzero(data_category)), 1.0)
                if math.isfinite(expected_data_err)
                else math.nan
            ),
            fit_deviance=fit_deviance,
            fit_ndf=fit_ndf,
        )
        rows.append(row)
        fit_metadata[f"{detector}_{sector}"] = fit_details

        log(
            f"{period.label} {detector}"
            + (f" sector {sector}" if detector == "FD" else "")
            + f": Nexp_data={expected_data:.1f}, Nfound_data={found_data:.0f}, "
            + f"eps_data={efficiency_data:.4f}, Nexp_MC={expected_mc:.0f}, "
            + f"Nfound_MC={found_mc:.0f}, eps_MC={efficiency_mc:.4f}, "
            + f"S_gamma={correction:.4f}"
        )
    # endfor

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
        "epgammagamma_identity_catalogs": {
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


def main() -> int:
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
