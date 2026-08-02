#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v18_binomial_helper_fix.py

Pass/fail simultaneous extraction of the RGA high-energy photon reconstruction
efficiency.

This is a clean replacement for the earlier inclusive BH/DVCS + pi0 denominator
fit. The selected epgamma opportunity sample is split first:

    PASS: a native epgammagamma record is found, its two daughter photons are
          reconstructed, the epgamma tag is identified, and the opposite
          daughter has E_gamma >= 2 GeV;

    FAIL: the same selected epgamma opportunity has no accepted reconstructed
          partner.

For each run period and predicted photon detector category (FT and FD sectors
1--6), the script performs a simultaneous extended Poisson fit to the PASS and
FAIL histograms. The model is

    PASS = N_pi0 * epsilon_data * T_pi0_pass

    FAIL = N_pi0 * (1 - epsilon_data) * T_pi0_fail
         + N_DVCS * T_DVCS_fail.

The AAOGEN PASS and FAIL templates are obtained with the same validated native
epgammagamma truth-partner matcher used for data. The DVCSGEN template enters
only the FAIL sample: native epgammagamma matching defines a reconstructed pi0
partner, so a genuine BH/DVCS event has no corresponding signal PASS component.

The shared fit parameters are N_pi0, N_DVCS, and epsilon_data. Therefore the
fitted total pi0 yield can never be smaller than its fitted PASS component.
The MC efficiency is measured directly as N_AAOGEN_pass / N_AAOGEN_total, and
the correction is S_gamma = epsilon_data / epsilon_MC.

The nominal fit drivers are Delta_phi, theta_gamma_gamma, and pTmiss. All seven
stored exclusivity variables are nevertheless plotted in four-panel shape
diagnostics:

    1. data PASS versus data FAIL;
    2. matched data versus AAOGEN PASS;
    3. data FAIL versus DVCSGEN and AAOGEN FAIL;
    4. data PASS versus AAOGEN PASS, with DVCSGEN shown only as a reference.

The tag-dependent cuts z(tag), -t1(tag), and electron-tag opening angle remain
disabled. The retained opportunity requirements are Q2, W, y, tag E >= 0.4 GeV,
predicted probe E >= 2 GeV, predicted-probe massless consistency, and predicted
FT/FD acceptance.

Output directory:
    output/photon_efficiency_study/pass_fail_simultaneous_fit/

The script writes complete cutflows, shape diagnostics, simultaneous-fit
overlays and pulls, bin-by-bin component tables, closure summaries, per-period
CSV/JSON products, and all-period efficiency/scale-factor canvases.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
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
from scipy.optimize import minimize
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





@dataclass
class PassFailEfficiencyRow:
    period: str
    period_label: str
    detector: str
    sector: int

    data_pass_count: int
    data_fail_count: int
    data_total_count: int

    pi0_mc_pass_count: int
    pi0_mc_fail_count: int
    pi0_mc_total_count: int
    dvcs_mc_count: int

    fitted_pi0_total: float
    fitted_pi0_total_err: float
    fitted_dvcs_total: float
    fitted_dvcs_total_err: float
    fitted_pi0_pass: float
    fitted_pi0_fail: float

    efficiency_data: float
    efficiency_data_err: float
    efficiency_mc: float
    efficiency_mc_err: float
    scale_factor: float
    scale_factor_err: float

    fit_success: bool
    fit_status: str
    fit_message: str
    fit_nll: float
    fit_deviance: float
    fit_ndf: int
    fit_reduced_deviance: float
    fit_quality_warning: bool
    fit_warning_reasons: str
    covariance_valid: bool

    closure_pass_minus_fit: float
    closure_pass_pull: float
    closure_total_data_minus_model: float


@dataclass
class SimultaneousFitResult:
    success: bool
    status: str
    message: str
    n_pi0: float
    n_dvcs: float
    epsilon_data: float
    covariance: np.ndarray
    errors: np.ndarray
    nll: float
    deviance: float
    ndf: int
    reduced_deviance: float
    covariance_valid: bool
    warning_reasons: List[str]
    driver_payload: Dict[str, Dict[str, object]]


PASS_FAIL_DRIVER_KEYS: Tuple[str, ...] = FRACTION_DRIVERS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract data/MC photon-efficiency scale factors with a simultaneous "
            "PASS/FAIL pi0 tag-and-probe fit."
        )
    )
    parser.add_argument("--period", action="append", choices=[p.key for p in PERIODS])
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument(
        "--output-dir",
        default=(
            "output/photon_efficiency_study/"
            "pass_fail_simultaneous_fit"
        ),
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
    # Retained only because the cutflow definitions document the disabled cuts.
    parser.add_argument("--z-min", type=float, default=0.65)
    parser.add_argument("--minus-t-max", type=float, default=1.0)
    parser.add_argument("--open-angle-min-deg", type=float, default=5.0)
    parser.add_argument("--require-fiducial-status-111", action="store_true")

    parser.add_argument("--tag-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument("--tag-match-relative-E-max", type=float, default=0.35)
    parser.add_argument("--probe-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument("--probe-match-relative-E-max", type=float, default=0.35)
    parser.add_argument("--mc-signature-decimals", type=int, default=10)
    parser.add_argument("--native-closure-max", type=float, default=0.10)
    parser.add_argument("--native-minimum-photon-E", type=float, default=0.20)
    parser.add_argument("--found-probe-E-min", type=float, default=2.00)
    parser.add_argument(
        "--require-probe-residual-match",
        action="store_true",
        help=(
            "Require the reconstructed partner to satisfy the predicted-probe "
            "angular and relative-energy residual windows. By default the "
            "validated identity and tag match define PASS."
        ),
    )

    parser.add_argument(
        "--fit-drivers",
        nargs="+",
        default=list(PASS_FAIL_DRIVER_KEYS),
        choices=[variable.key for variable in FIT_VARIABLES],
        help=(
            "Variables used in the simultaneous PASS/FAIL likelihood. All seven "
            "variables are still plotted."
        ),
    )
    parser.add_argument("--fit-min-data-pass", type=int, default=50)
    parser.add_argument("--fit-min-data-fail", type=int, default=100)
    parser.add_argument("--fit-min-template-counts", type=int, default=100)
    parser.add_argument(
        "--fit-max-reduced-deviance",
        type=float,
        default=5.0,
        help=(
            "Warning threshold only. Poor goodness of fit does not stop output "
            "production."
        ),
    )
    parser.add_argument(
        "--template-pseudocount",
        type=float,
        default=0.25,
        help=(
            "Small symmetric pseudocount added to each template bin before "
            "normalization. It prevents zero-probability bins without changing "
            "the template normalization appreciably."
        ),
    )
    parser.add_argument(
        "--fit-maxiter",
        type=int,
        default=2000,
    )
    parser.add_argument(
        "--fail-on-fit-failure",
        action="store_true",
        help=(
            "Abort when a category fit is mathematically invalid. By default "
            "the category is marked invalid and the remaining diagnostics run."
        ),
    )
    return parser.parse_args()


def category_stem(detector: str, sector: int) -> str:
    return "ft" if detector == "FT" else f"fd_sector_{sector}"


def category_title(detector: str, sector: int) -> str:
    return "FT" if detector == "FT" else f"FD sector {sector}"


def finite_values(records: OpportunityRecords, key: str, mask: np.ndarray) -> np.ndarray:
    values = np.asarray(getattr(records, key), dtype=float)[mask]
    return values[np.isfinite(values)]


def normalized_histogram(
    records: OpportunityRecords,
    mask: np.ndarray,
    variable: FitVariable,
    pseudocount: float,
) -> Tuple[np.ndarray, np.ndarray]:
    values = finite_values(records, variable.key, mask)
    counts, edges = np.histogram(
        values,
        bins=variable.bins,
        range=(variable.low, variable.high),
    )
    smoothed = counts.astype(float) + max(float(pseudocount), 0.0)
    total = float(np.sum(smoothed))
    shape = (
        smoothed / total
        if total > 0.0
        else np.full(variable.bins, 1.0 / variable.bins)
    )
    return shape, counts.astype(float)


def density_histogram(
    records: OpportunityRecords,
    mask: np.ndarray,
    variable: FitVariable,
) -> Tuple[np.ndarray, np.ndarray]:
    values = finite_values(records, variable.key, mask)
    counts, edges = np.histogram(
        values,
        bins=variable.bins,
        range=(variable.low, variable.high),
    )
    widths = np.diff(edges)
    total = float(np.sum(counts))
    density = (
        counts / (total * widths)
        if total > 0.0
        else np.zeros_like(counts, dtype=float)
    )
    return density.astype(float), edges


def poisson_nll(observed: np.ndarray, expected: np.ndarray) -> float:
    expected = np.maximum(np.asarray(expected, dtype=float), 1.0e-12)
    observed = np.asarray(observed, dtype=float)
    return float(np.sum(expected - xlogy(observed, expected)))


def poisson_deviance(observed: np.ndarray, expected: np.ndarray) -> float:
    expected = np.maximum(np.asarray(expected, dtype=float), 1.0e-12)
    observed = np.asarray(observed, dtype=float)
    terms = expected - observed + xlogy(observed, observed / expected)
    terms = np.where(observed > 0.0, terms, expected)
    return float(2.0 * np.sum(terms))


def numerical_hessian(
    function,
    point: np.ndarray,
    lower: np.ndarray,
    upper: np.ndarray,
) -> np.ndarray:
    point = np.asarray(point, dtype=float)
    n = point.size
    steps = np.maximum(np.abs(point) * 2.0e-4, np.asarray([1.0, 1.0, 2.0e-5]))
    hessian = np.zeros((n, n), dtype=float)
    f0 = float(function(point))

    for i in range(n):
        hi = min(
            steps[i],
            0.45 * max(point[i] - lower[i], steps[i]),
            0.45 * max(upper[i] - point[i], steps[i]),
        )
        hi = max(hi, 1.0e-7)
        plus = point.copy()
        minus = point.copy()
        plus[i] += hi
        minus[i] -= hi
        hessian[i, i] = (
            function(plus) - 2.0 * f0 + function(minus)
        ) / (hi * hi)

        for j in range(i + 1, n):
            hj = min(
                steps[j],
                0.45 * max(point[j] - lower[j], steps[j]),
                0.45 * max(upper[j] - point[j], steps[j]),
            )
            hj = max(hj, 1.0e-7)
            pp = point.copy()
            pm = point.copy()
            mp = point.copy()
            mm = point.copy()
            pp[i] += hi
            pp[j] += hj
            pm[i] += hi
            pm[j] -= hj
            mp[i] -= hi
            mp[j] += hj
            mm[i] -= hi
            mm[j] -= hj
            mixed = (
                function(pp) - function(pm) - function(mp) + function(mm)
            ) / (4.0 * hi * hj)
            hessian[i, j] = mixed
            hessian[j, i] = mixed
        # endfor
    # endfor
    return hessian


def build_driver_payload(
    data: OpportunityRecords,
    dvcs: OpportunityRecords,
    pi0: OpportunityRecords,
    data_pass: np.ndarray,
    data_fail: np.ndarray,
    pi0_pass: np.ndarray,
    pi0_fail: np.ndarray,
    dvcs_mask: np.ndarray,
    args: argparse.Namespace,
) -> Dict[str, Dict[str, object]]:
    variable_by_key = {variable.key: variable for variable in FIT_VARIABLES}
    payload: Dict[str, Dict[str, object]] = {}

    for key in args.fit_drivers:
        variable = variable_by_key[key]
        data_pass_counts, edges = np.histogram(
            finite_values(data, key, data_pass),
            bins=variable.bins,
            range=(variable.low, variable.high),
        )
        data_fail_counts, _ = np.histogram(
            finite_values(data, key, data_fail),
            bins=variable.bins,
            range=(variable.low, variable.high),
        )
        pi0_pass_shape, pi0_pass_raw = normalized_histogram(
            pi0, pi0_pass, variable, args.template_pseudocount
        )
        pi0_fail_shape, pi0_fail_raw = normalized_histogram(
            pi0, pi0_fail, variable, args.template_pseudocount
        )
        dvcs_fail_shape, dvcs_raw = normalized_histogram(
            dvcs, dvcs_mask, variable, args.template_pseudocount
        )

        payload[key] = {
            "variable": variable,
            "edges": edges,
            "data_pass": data_pass_counts.astype(float),
            "data_fail": data_fail_counts.astype(float),
            "pi0_pass_shape": pi0_pass_shape,
            "pi0_fail_shape": pi0_fail_shape,
            "dvcs_fail_shape": dvcs_fail_shape,
            "pi0_pass_raw": pi0_pass_raw,
            "pi0_fail_raw": pi0_fail_raw,
            "dvcs_raw": dvcs_raw,
        }
    # endfor
    return payload


def simultaneous_pass_fail_fit(
    driver_payload: Dict[str, Dict[str, object]],
    data_pass_count: int,
    data_fail_count: int,
    args: argparse.Namespace,
) -> SimultaneousFitResult:
    if data_pass_count < args.fit_min_data_pass:
        raise RuntimeError(
            f"PASS data support {data_pass_count} is below "
            f"--fit-min-data-pass={args.fit_min_data_pass}"
        )
    # endif
    if data_fail_count < args.fit_min_data_fail:
        raise RuntimeError(
            f"FAIL data support {data_fail_count} is below "
            f"--fit-min-data-fail={args.fit_min_data_fail}"
        )
    # endif

    for key, item in driver_payload.items():
        for template_name in ("pi0_pass_raw", "pi0_fail_raw", "dvcs_raw"):
            support = int(np.sum(item[template_name]))
            if support < args.fit_min_template_counts:
                raise RuntimeError(
                    f"{key} {template_name} support {support} is below "
                    f"--fit-min-template-counts={args.fit_min_template_counts}"
                )
            # endif
        # endfor
    # endfor

    total_data = float(data_pass_count + data_fail_count)
    lower = np.asarray(
        [max(float(data_pass_count), 1.0e-6), 0.0, 1.0e-6],
        dtype=float,
    )
    upper = np.asarray(
        [
            max(10.0 * total_data, lower[0] + 1.0),
            max(10.0 * total_data, 1.0),
            1.0 - 1.0e-6,
        ],
        dtype=float,
    )

    observed_pass_fraction = data_pass_count / max(total_data, 1.0)
    initial_epsilon = float(np.clip(observed_pass_fraction, 0.05, 0.95))
    initial_pi0 = max(
        float(data_pass_count) / max(initial_epsilon, 1.0e-6),
        float(data_pass_count) + 1.0,
    )
    initial_dvcs = max(total_data - initial_pi0, 0.1 * total_data, 1.0)
    x0 = np.asarray(
        [
            min(initial_pi0, upper[0] * 0.8),
            min(initial_dvcs, upper[1] * 0.8),
            initial_epsilon,
        ],
        dtype=float,
    )
    x0 = np.minimum(np.maximum(x0, lower + 1.0e-7), upper - 1.0e-7)

    def objective(parameters: np.ndarray) -> float:
        n_pi0, n_dvcs, efficiency = parameters
        if (
            n_pi0 < lower[0]
            or n_dvcs < 0.0
            or not (0.0 < efficiency < 1.0)
        ):
            return 1.0e100
        # endif

        value = 0.0
        for item in driver_payload.values():
            expected_pass = (
                n_pi0 * efficiency * item["pi0_pass_shape"]
            )
            expected_fail = (
                n_pi0 * (1.0 - efficiency) * item["pi0_fail_shape"]
                + n_dvcs * item["dvcs_fail_shape"]
            )
            value += poisson_nll(item["data_pass"], expected_pass)
            value += poisson_nll(item["data_fail"], expected_fail)
        # endfor
        return float(value)

    result = minimize(
        objective,
        x0,
        method="L-BFGS-B",
        bounds=list(zip(lower, upper)),
        options={
            "maxiter": int(args.fit_maxiter),
            "ftol": 1.0e-12,
            "gtol": 1.0e-8,
            "maxls": 50,
        },
    )

    point = np.asarray(result.x, dtype=float)
    n_pi0, n_dvcs, efficiency = point

    covariance_valid = False
    covariance = np.full((3, 3), np.nan)
    errors = np.full(3, np.nan)
    try:
        hessian = numerical_hessian(objective, point, lower, upper)
        eigenvalues = np.linalg.eigvalsh(hessian)
        if np.all(np.isfinite(hessian)) and np.all(eigenvalues > 0.0):
            covariance = np.linalg.inv(hessian)
            diagonal = np.diag(covariance)
            if np.all(diagonal >= 0.0) and np.all(np.isfinite(diagonal)):
                errors = np.sqrt(diagonal)
                covariance_valid = True
            # endif
        # endif
    except Exception:
        covariance_valid = False
    # endtry

    deviance = 0.0
    number_of_bins = 0
    fitted_payload: Dict[str, Dict[str, object]] = {}
    for key, item in driver_payload.items():
        expected_pass = n_pi0 * efficiency * item["pi0_pass_shape"]
        pi0_fail_component = (
            n_pi0 * (1.0 - efficiency) * item["pi0_fail_shape"]
        )
        dvcs_fail_component = n_dvcs * item["dvcs_fail_shape"]
        expected_fail = pi0_fail_component + dvcs_fail_component

        deviance += poisson_deviance(item["data_pass"], expected_pass)
        deviance += poisson_deviance(item["data_fail"], expected_fail)
        number_of_bins += int(item["data_pass"].size + item["data_fail"].size)

        fitted_payload[key] = {
            **item,
            "expected_pass": expected_pass,
            "expected_fail": expected_fail,
            "pi0_fail_component": pi0_fail_component,
            "dvcs_fail_component": dvcs_fail_component,
        }
    # endfor

    ndf = max(number_of_bins - 3, 1)
    reduced_deviance = deviance / ndf
    warning_reasons: List[str] = []

    if (
        args.fit_max_reduced_deviance > 0.0
        and reduced_deviance > args.fit_max_reduced_deviance
    ):
        warning_reasons.append(
            f"deviance/ndf={reduced_deviance:.4g} exceeds "
            f"{args.fit_max_reduced_deviance:g}"
        )
    # endif
    if not covariance_valid:
        warning_reasons.append("numerical covariance is not positive definite")
    # endif
    if efficiency < 1.0e-4 or efficiency > 1.0 - 1.0e-4:
        warning_reasons.append(
            f"epsilon_data={efficiency:.6g} is near a physical boundary"
        )
    # endif

    success = bool(
        result.success
        and np.all(np.isfinite(point))
        and n_pi0 >= data_pass_count - 1.0e-6
        and n_dvcs >= 0.0
        and 0.0 < efficiency < 1.0
    )
    status = "success" if success else "failed"
    message = str(result.message)

    return SimultaneousFitResult(
        success=success,
        status=status,
        message=message,
        n_pi0=float(n_pi0),
        n_dvcs=float(n_dvcs),
        epsilon_data=float(efficiency),
        covariance=covariance,
        errors=errors,
        nll=float(result.fun),
        deviance=float(deviance),
        ndf=int(ndf),
        reduced_deviance=float(reduced_deviance),
        covariance_valid=covariance_valid,
        warning_reasons=warning_reasons,
        driver_payload=fitted_payload,
    )


def plot_shape_diagnostics(
    output_dir: Path,
    period_label: str,
    detector: str,
    sector: int,
    data: OpportunityRecords,
    dvcs: OpportunityRecords,
    pi0: OpportunityRecords,
    data_pass: np.ndarray,
    data_fail: np.ndarray,
    pi0_pass: np.ndarray,
    pi0_fail: np.ndarray,
    dvcs_mask: np.ndarray,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    title_category = category_title(detector, sector)

    for variable in FIT_VARIABLES:
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))

        panels = (
            (
                axes[0, 0],
                (
                    ("Data PASS", data, data_pass),
                    ("Data FAIL", data, data_fail),
                ),
                "Data PASS versus FAIL",
            ),
            (
                axes[0, 1],
                (
                    ("Matched data", data, data_pass),
                    ("AAOGEN PASS", pi0, pi0_pass),
                ),
                "Matched data versus AAOGEN PASS",
            ),
            (
                axes[1, 0],
                (
                    ("Data FAIL", data, data_fail),
                    ("DVCSGEN", dvcs, dvcs_mask),
                    ("AAOGEN FAIL", pi0, pi0_fail),
                ),
                "FAIL population templates",
            ),
            (
                axes[1, 1],
                (
                    ("Data PASS", data, data_pass),
                    ("AAOGEN PASS", pi0, pi0_pass),
                    ("DVCSGEN reference", dvcs, dvcs_mask),
                ),
                "PASS control population",
            ),
        )

        for axis, curves, panel_title in panels:
            for label, records, mask in curves:
                density, edges = density_histogram(records, mask, variable)
                centers = 0.5 * (edges[:-1] + edges[1:])
                axis.step(
                    centers,
                    density,
                    where="mid",
                    linewidth=1.3,
                    label=f"{label} ({int(np.count_nonzero(mask)):,})",
                )
            # endfor
            axis.set_title(panel_title)
            axis.set_xlabel(variable.label)
            axis.set_ylabel("Unit-normalized density")
            axis.grid(alpha=0.25)
            axis.legend(fontsize=8)
        # endfor

        fig.suptitle(
            f"{period_label}, {title_category}: {variable.label}\n"
            "All observables are the original epgamma tag-based quantities"
        )
        fig.tight_layout(rect=(0, 0, 1, 0.94))
        fig.savefig(
            output_dir / f"{variable.key}_pass_fail_shapes.png",
            dpi=180,
        )
        plt.close(fig)
    # endfor


def plot_simultaneous_fit(
    path: Path,
    period_label: str,
    detector: str,
    sector: int,
    fit: SimultaneousFitResult,
) -> None:
    driver_keys = list(fit.driver_payload)
    n_columns = len(driver_keys)
    fig, axes = plt.subplots(
        4,
        n_columns,
        figsize=(6.0 * n_columns, 14),
        squeeze=False,
    )

    for column, key in enumerate(driver_keys):
        item = fit.driver_payload[key]
        variable = item["variable"]
        edges = np.asarray(item["edges"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])

        data_pass = np.asarray(item["data_pass"], dtype=float)
        expected_pass = np.asarray(item["expected_pass"], dtype=float)
        data_fail = np.asarray(item["data_fail"], dtype=float)
        expected_fail = np.asarray(item["expected_fail"], dtype=float)
        pi0_fail = np.asarray(item["pi0_fail_component"], dtype=float)
        dvcs_fail = np.asarray(item["dvcs_fail_component"], dtype=float)

        axes[0, column].errorbar(
            centers,
            data_pass,
            yerr=np.sqrt(np.maximum(data_pass, 1.0)),
            fmt=".",
            label="Data PASS",
        )
        axes[0, column].step(
            centers,
            expected_pass,
            where="mid",
            linewidth=1.5,
            label="Fitted pi0 PASS",
        )
        axes[0, column].set_title(f"PASS: {variable.label}")
        axes[0, column].set_ylabel("Events / bin")
        axes[0, column].grid(alpha=0.25)
        axes[0, column].legend(fontsize=8)

        pass_pull = (
            data_pass - expected_pass
        ) / np.sqrt(np.maximum(expected_pass, 1.0))
        axes[1, column].axhline(0.0, linewidth=1.0)
        axes[1, column].axhline(3.0, linestyle="--", linewidth=0.8)
        axes[1, column].axhline(-3.0, linestyle="--", linewidth=0.8)
        axes[1, column].plot(centers, pass_pull, ".", markersize=3)
        axes[1, column].set_ylabel("PASS pull")
        axes[1, column].grid(alpha=0.25)

        axes[2, column].errorbar(
            centers,
            data_fail,
            yerr=np.sqrt(np.maximum(data_fail, 1.0)),
            fmt=".",
            label="Data FAIL",
        )
        axes[2, column].step(
            centers,
            expected_fail,
            where="mid",
            linewidth=1.5,
            label="Total FAIL fit",
        )
        axes[2, column].step(
            centers,
            pi0_fail,
            where="mid",
            linewidth=1.0,
            label="pi0 FAIL",
        )
        axes[2, column].step(
            centers,
            dvcs_fail,
            where="mid",
            linewidth=1.0,
            label="BH/DVCS FAIL",
        )
        axes[2, column].set_ylabel("Events / bin")
        axes[2, column].grid(alpha=0.25)
        axes[2, column].legend(fontsize=8)

        fail_pull = (
            data_fail - expected_fail
        ) / np.sqrt(np.maximum(expected_fail, 1.0))
        axes[3, column].axhline(0.0, linewidth=1.0)
        axes[3, column].axhline(3.0, linestyle="--", linewidth=0.8)
        axes[3, column].axhline(-3.0, linestyle="--", linewidth=0.8)
        axes[3, column].plot(centers, fail_pull, ".", markersize=3)
        axes[3, column].set_ylabel("FAIL pull")
        axes[3, column].set_xlabel(variable.label)
        axes[3, column].grid(alpha=0.25)
    # endfor

    category = category_title(detector, sector)
    error_epsilon = (
        fit.errors[2]
        if fit.errors.size >= 3 and math.isfinite(fit.errors[2])
        else math.nan
    )
    fig.suptitle(
        f"{period_label}, {category}: simultaneous PASS/FAIL fit\n"
        rf"$N_{{\pi^0}}={fit.n_pi0:.1f}$, "
        rf"$N_{{\rm DVCS}}={fit.n_dvcs:.1f}$, "
        rf"$\epsilon_{{\rm data}}={fit.epsilon_data:.5f}"
        rf"\pm{error_epsilon:.5f}$, "
        rf"$D/\mathrm{{ndf}}={fit.reduced_deviance:.3f}$"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_fit_component_tables(
    output_dir: Path,
    stem: str,
    fit: SimultaneousFitResult,
) -> Dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: Dict[str, str] = {}

    for key, item in fit.driver_payload.items():
        edges = np.asarray(item["edges"], dtype=float)
        data_pass = np.asarray(item["data_pass"], dtype=float)
        expected_pass = np.asarray(item["expected_pass"], dtype=float)
        data_fail = np.asarray(item["data_fail"], dtype=float)
        expected_fail = np.asarray(item["expected_fail"], dtype=float)
        pi0_fail = np.asarray(item["pi0_fail_component"], dtype=float)
        dvcs_fail = np.asarray(item["dvcs_fail_component"], dtype=float)

        path = output_dir / f"{stem}_{key}_components.csv"
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "bin_low",
                    "bin_high",
                    "bin_center",
                    "data_pass",
                    "fit_pass_pi0",
                    "pass_pull",
                    "data_fail",
                    "fit_fail_total",
                    "fit_fail_pi0",
                    "fit_fail_bh_dvcs",
                    "fail_pull",
                    "pass_deviance_contribution",
                    "fail_deviance_contribution",
                ]
            )
            for index in range(data_pass.size):
                pass_mu = max(expected_pass[index], 1.0e-12)
                fail_mu = max(expected_fail[index], 1.0e-12)
                pass_obs = data_pass[index]
                fail_obs = data_fail[index]
                pass_dev = 2.0 * (
                    pass_mu - pass_obs
                    + (
                        pass_obs * math.log(pass_obs / pass_mu)
                        if pass_obs > 0.0
                        else 0.0
                    )
                )
                fail_dev = 2.0 * (
                    fail_mu - fail_obs
                    + (
                        fail_obs * math.log(fail_obs / fail_mu)
                        if fail_obs > 0.0
                        else 0.0
                    )
                )
                writer.writerow(
                    [
                        edges[index],
                        edges[index + 1],
                        0.5 * (edges[index] + edges[index + 1]),
                        pass_obs,
                        pass_mu,
                        (pass_obs - pass_mu) / math.sqrt(max(pass_mu, 1.0)),
                        fail_obs,
                        fail_mu,
                        pi0_fail[index],
                        dvcs_fail[index],
                        (fail_obs - fail_mu) / math.sqrt(max(fail_mu, 1.0)),
                        pass_dev,
                        fail_dev,
                    ]
                )
            # endfor
        # endwith
        paths[key] = str(path)
    # endfor
    return paths


def plot_closure_summary(
    path: Path,
    period_label: str,
    detector: str,
    sector: int,
    row: PassFailEfficiencyRow,
) -> None:
    category = category_title(detector, sector)
    labels = [
        "Data PASS",
        "Fitted pi0 PASS",
        "Data FAIL",
        "Fitted pi0 FAIL",
        "Fitted BH/DVCS FAIL",
    ]
    values = [
        row.data_pass_count,
        row.fitted_pi0_pass,
        row.data_fail_count,
        row.fitted_pi0_fail,
        row.fitted_dvcs_total,
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    x = np.arange(len(labels))
    axes[0].bar(x, values)
    axes[0].set_xticks(x, labels, rotation=25, ha="right")
    axes[0].set_ylabel("Events")
    axes[0].set_title("Fitted event decomposition")
    axes[0].grid(axis="y", alpha=0.25)

    ratio_labels = [
        r"$\epsilon_{\rm data}$",
        r"$\epsilon_{\rm MC}$",
        r"$S_\gamma$",
    ]
    ratio_values = [
        row.efficiency_data,
        row.efficiency_mc,
        row.scale_factor,
    ]
    ratio_errors = [
        row.efficiency_data_err,
        row.efficiency_mc_err,
        row.scale_factor_err,
    ]
    xr = np.arange(len(ratio_labels))
    axes[1].errorbar(
        xr,
        ratio_values,
        yerr=ratio_errors,
        fmt="o",
        capsize=3,
    )
    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1].set_xticks(xr, ratio_labels)
    axes[1].set_ylabel("Efficiency or ratio")
    finite = [value for value in ratio_values if math.isfinite(value)]
    axes[1].set_ylim(0.0, max([1.1] + [1.2 * value for value in finite]))
    axes[1].grid(alpha=0.25)

    fig.suptitle(
        f"{period_label}, {category}: PASS/FAIL closure\n"
        f"PASS data - fitted pi0 PASS = {row.closure_pass_minus_fit:.2f}, "
        f"PASS pull = {row.closure_pass_pull:.3f}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def binomial_efficiency(
    passed: float,
    total: float,
) -> Tuple[float, float]:
    """
    Return a binomial efficiency and its Gaussian statistical uncertainty.

    Parameters
    ----------
    passed
        Number of successful trials.
    total
        Total number of trials.

    Returns
    -------
    efficiency, uncertainty
        The efficiency passed / total and
        sqrt(efficiency * (1 - efficiency) / total).

    Notes
    -----
    Invalid populations return (nan, nan). A tiny floating-point tolerance is
    allowed at the physical boundaries, but genuinely unphysical inputs such
    as passed < 0 or passed > total are rejected.
    """
    if (
        not math.isfinite(passed)
        or not math.isfinite(total)
        or total <= 0.0
        or passed < 0.0
        or passed > total + 1.0e-9
    ):
        return math.nan, math.nan
    # endif

    efficiency = passed / total
    efficiency = min(max(efficiency, 0.0), 1.0)
    variance = efficiency * (1.0 - efficiency) / total

    return efficiency, math.sqrt(max(variance, 0.0))



def uncertainty_ratio(
    numerator: float,
    numerator_error: float,
    denominator: float,
    denominator_error: float,
) -> Tuple[float, float]:
    if (
        not math.isfinite(numerator)
        or not math.isfinite(denominator)
        or denominator <= 0.0
    ):
        return math.nan, math.nan
    # endif
    ratio = numerator / denominator
    if numerator <= 0.0:
        return ratio, math.nan
    # endif
    relative_variance = (
        (numerator_error / numerator) ** 2
        + (denominator_error / denominator) ** 2
    )
    return ratio, abs(ratio) * math.sqrt(max(relative_variance, 0.0))




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

def process_period(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    cutflow_dir = period_dir / "cutflow_diagnostics"
    shape_root = period_dir / "shape_diagnostics"
    fit_root = period_dir / "simultaneous_fits"
    closure_root = period_dir / "closure_diagnostics"
    table_root = period_dir / "fit_component_tables"

    for directory in (
        period_dir,
        cutflow_dir,
        shape_root,
        fit_root,
        closure_root,
        table_root,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    # endfor

    data, _, data_cutflow = read_opportunities(
        period.epg_data,
        period.beam_energy_GeV,
        "epgamma data",
        args,
        deduplicate=False,
    )
    dvcs, _, dvcs_cutflow = read_opportunities(
        period.dvcs_mc,
        period.beam_energy_GeV,
        "DVCSGEN epgamma MC",
        args,
        deduplicate=True,
    )
    pi0, _, pi0_cutflow = read_opportunities(
        period.pi0_epg_mc,
        period.beam_energy_GeV,
        "AAOGEN-as-epgamma MC",
        args,
        deduplicate=True,
    )

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
        period.epgg_data,
        period.beam_energy_GeV,
        "epgammagamma data",
        args,
    )
    epgg_mc, epgg_mc_diag = read_epgg(
        period.pi0_epgg_mc,
        period.beam_energy_GeV,
        "AAOGEN epgammagamma MC",
        args,
    )

    data_match = match_truth_partners(data, epgg_data, "data", args)
    pi0_match = match_truth_partners(pi0, epgg_mc, "mc", args)

    for sample_label, result in (
        ("data", data_match),
        ("AAOGEN MC", pi0_match),
    ):
        summary = result.summary
        log(
            f"{period.label} {sample_label}: identity="
            f"{summary['identity_matches']:,}/{summary['opportunities']:,}, "
            f"valid|identity="
            f"{summary['conditional_valid_solution_given_identity']:.3%}, "
            f"tag|valid="
            f"{summary['conditional_tag_match_given_valid_solution']:.3%}, "
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
        pi0_match,
    )

    rows: List[PassFailEfficiencyRow] = []
    category_metadata: Dict[str, object] = {}
    categories = [("FT", 0)] + [("FD", sector) for sector in range(1, 7)]

    for detector, sector in categories:
        stem = category_stem(detector, sector)
        title = category_title(detector, sector)

        data_category = category_mask(data, detector, sector)
        dvcs_category = category_mask(dvcs, detector, sector)
        pi0_category = category_mask(pi0, detector, sector)

        data_pass = data_category & data_match.matched
        data_fail = data_category & ~data_match.matched
        pi0_pass = pi0_category & pi0_match.matched
        pi0_fail = pi0_category & ~pi0_match.matched

        data_pass_count = int(np.count_nonzero(data_pass))
        data_fail_count = int(np.count_nonzero(data_fail))
        pi0_pass_count = int(np.count_nonzero(pi0_pass))
        pi0_fail_count = int(np.count_nonzero(pi0_fail))
        dvcs_count = int(np.count_nonzero(dvcs_category))

        shape_dir = shape_root / stem
        plot_shape_diagnostics(
            shape_dir,
            period.label,
            detector,
            sector,
            data,
            dvcs,
            pi0,
            data_pass,
            data_fail,
            pi0_pass,
            pi0_fail,
            dvcs_category,
        )

        driver_payload = build_driver_payload(
            data,
            dvcs,
            pi0,
            data_pass,
            data_fail,
            pi0_pass,
            pi0_fail,
            dvcs_category,
            args,
        )

        fit_failure: Optional[str] = None
        try:
            fit = simultaneous_pass_fail_fit(
                driver_payload,
                data_pass_count,
                data_fail_count,
                args,
            )
        except Exception as exc:
            fit_failure = str(exc)
            fit = SimultaneousFitResult(
                success=False,
                status="failed",
                message=fit_failure,
                n_pi0=math.nan,
                n_dvcs=math.nan,
                epsilon_data=math.nan,
                covariance=np.full((3, 3), np.nan),
                errors=np.full(3, np.nan),
                nll=math.nan,
                deviance=math.nan,
                ndf=0,
                reduced_deviance=math.nan,
                covariance_valid=False,
                warning_reasons=[fit_failure],
                driver_payload={},
            )
            if args.fail_on_fit_failure:
                raise
            # endif
            log(
                f"WARNING: {period.label} {title}: simultaneous fit failed; "
                f"diagnostics retained and official values set to NaN: "
                f"{fit_failure}"
            )
        # endtry

        if fit.driver_payload:
            fit_plot = fit_root / f"{stem}_simultaneous_fit.png"
            plot_simultaneous_fit(
                fit_plot,
                period.label,
                detector,
                sector,
                fit,
            )
            component_tables = write_fit_component_tables(
                table_root,
                stem,
                fit,
            )
        else:
            fit_plot = None
            component_tables = {}
        # endif

        efficiency_mc, efficiency_mc_err = binomial_efficiency(
            float(pi0_pass_count),
            float(pi0_pass_count + pi0_fail_count),
        )

        efficiency_data = (
            fit.epsilon_data if fit.success else math.nan
        )
        efficiency_data_err = (
            float(fit.errors[2])
            if (
                fit.success
                and fit.errors.size >= 3
                and math.isfinite(fit.errors[2])
            )
            else math.nan
        )
        scale, scale_err = uncertainty_ratio(
            efficiency_data,
            efficiency_data_err,
            efficiency_mc,
            efficiency_mc_err,
        )

        fitted_pi0_pass = (
            fit.n_pi0 * fit.epsilon_data
            if fit.success
            else math.nan
        )
        fitted_pi0_fail = (
            fit.n_pi0 * (1.0 - fit.epsilon_data)
            if fit.success
            else math.nan
        )
        fitted_pi0_error = (
            float(fit.errors[0])
            if fit.errors.size >= 1 and math.isfinite(fit.errors[0])
            else math.nan
        )
        fitted_dvcs_error = (
            float(fit.errors[1])
            if fit.errors.size >= 2 and math.isfinite(fit.errors[1])
            else math.nan
        )

        pass_difference = (
            data_pass_count - fitted_pi0_pass
            if math.isfinite(fitted_pi0_pass)
            else math.nan
        )
        pass_pull = (
            pass_difference / math.sqrt(max(fitted_pi0_pass, 1.0))
            if math.isfinite(pass_difference)
            and math.isfinite(fitted_pi0_pass)
            and fitted_pi0_pass > 0.0
            else math.nan
        )
        total_model = (
            fit.n_pi0 + fit.n_dvcs if fit.success else math.nan
        )
        total_difference = (
            data_pass_count + data_fail_count - total_model
            if math.isfinite(total_model)
            else math.nan
        )

        row = PassFailEfficiencyRow(
            period=period.key,
            period_label=period.label,
            detector=detector,
            sector=sector,
            data_pass_count=data_pass_count,
            data_fail_count=data_fail_count,
            data_total_count=data_pass_count + data_fail_count,
            pi0_mc_pass_count=pi0_pass_count,
            pi0_mc_fail_count=pi0_fail_count,
            pi0_mc_total_count=pi0_pass_count + pi0_fail_count,
            dvcs_mc_count=dvcs_count,
            fitted_pi0_total=fit.n_pi0,
            fitted_pi0_total_err=fitted_pi0_error,
            fitted_dvcs_total=fit.n_dvcs,
            fitted_dvcs_total_err=fitted_dvcs_error,
            fitted_pi0_pass=fitted_pi0_pass,
            fitted_pi0_fail=fitted_pi0_fail,
            efficiency_data=efficiency_data,
            efficiency_data_err=efficiency_data_err,
            efficiency_mc=efficiency_mc,
            efficiency_mc_err=efficiency_mc_err,
            scale_factor=scale,
            scale_factor_err=scale_err,
            fit_success=fit.success,
            fit_status=fit.status,
            fit_message=fit.message,
            fit_nll=fit.nll,
            fit_deviance=fit.deviance,
            fit_ndf=fit.ndf,
            fit_reduced_deviance=fit.reduced_deviance,
            fit_quality_warning=bool(fit.warning_reasons),
            fit_warning_reasons="; ".join(fit.warning_reasons),
            covariance_valid=fit.covariance_valid,
            closure_pass_minus_fit=pass_difference,
            closure_pass_pull=pass_pull,
            closure_total_data_minus_model=total_difference,
        )
        rows.append(row)

        closure_plot = closure_root / f"{stem}_closure.png"
        plot_closure_summary(
            closure_plot,
            period.label,
            detector,
            sector,
            row,
        )

        metadata_entry = {
            "category": title,
            "counts": {
                "data_pass": data_pass_count,
                "data_fail": data_fail_count,
                "pi0_mc_pass": pi0_pass_count,
                "pi0_mc_fail": pi0_fail_count,
                "dvcs_mc": dvcs_count,
            },
            "fit": {
                "success": fit.success,
                "status": fit.status,
                "message": fit.message,
                "n_pi0": fit.n_pi0,
                "n_dvcs": fit.n_dvcs,
                "epsilon_data": fit.epsilon_data,
                "errors": fit.errors.tolist(),
                "covariance": fit.covariance.tolist(),
                "covariance_valid": fit.covariance_valid,
                "nll": fit.nll,
                "deviance": fit.deviance,
                "ndf": fit.ndf,
                "reduced_deviance": fit.reduced_deviance,
                "warning_reasons": fit.warning_reasons,
                "drivers": list(args.fit_drivers),
            },
            "efficiency": {
                "data": efficiency_data,
                "data_error": efficiency_data_err,
                "mc": efficiency_mc,
                "mc_error": efficiency_mc_err,
                "scale_factor": scale,
                "scale_factor_error": scale_err,
            },
            "closure": {
                "fitted_pi0_pass": fitted_pi0_pass,
                "data_pass_minus_fit": pass_difference,
                "pass_pull": pass_pull,
                "total_data_minus_model": total_difference,
            },
            "outputs": {
                "shape_directory": str(shape_dir),
                "fit_plot": str(fit_plot) if fit_plot else None,
                "component_tables": component_tables,
                "closure_plot": str(closure_plot),
            },
        }
        category_metadata[stem] = metadata_entry

        warning_text = (
            " | WARN: " + "; ".join(fit.warning_reasons)
            if fit.warning_reasons
            else ""
        )
        log(
            f"{period.label} {title}: data PASS/FAIL="
            f"{data_pass_count:,}/{data_fail_count:,}, "
            f"AAOGEN PASS/FAIL={pi0_pass_count:,}/{pi0_fail_count:,}, "
            f"Npi0={fit.n_pi0:.1f}, NDVCS={fit.n_dvcs:.1f}, "
            f"eps_data={efficiency_data:.5f}, eps_MC={efficiency_mc:.5f}, "
            f"S_gamma={scale:.5f}, D/ndf={fit.reduced_deviance:.3f}"
            f"{warning_text}"
        )
    # endfor

    row_dicts = [asdict(row) for row in rows]
    write_rows_csv(
        period_dir / "pass_fail_efficiency_results.csv",
        row_dicts,
    )

    metadata = {
        "period": asdict(period),
        "arguments": vars(args),
        "opportunity_cutflows": initial_cutflows,
        "epgammagamma_reconstruction": {
            "data": epgg_data_diag,
            "mc": epgg_mc_diag,
        },
        "matching": {
            "data": data_match.summary,
            "mc": pi0_match.summary,
        },
        "categories": category_metadata,
        "rows": row_dicts,
    }
    with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    # endwith
    with open(period_dir / "matching_summary.json", "w", encoding="utf-8") as handle:
        json.dump(metadata["matching"], handle, indent=2)
    # endwith

    plot_period_pass_fail_results(
        period_dir / "integrated_efficiencies_and_scale_factors.png",
        period.label,
        rows,
    )

    return period.key, row_dicts, metadata


def plot_period_pass_fail_results(
    path: Path,
    title: str,
    rows: Sequence[PassFailEfficiencyRow],
) -> None:
    labels = [
        "FT" if row.detector == "FT" else f"FD S{row.sector}"
        for row in rows
    ]
    x = np.arange(len(rows), dtype=float)

    data_eff = np.asarray([row.efficiency_data for row in rows], dtype=float)
    data_err = np.asarray([row.efficiency_data_err for row in rows], dtype=float)
    mc_eff = np.asarray([row.efficiency_mc for row in rows], dtype=float)
    mc_err = np.asarray([row.efficiency_mc_err for row in rows], dtype=float)
    scale = np.asarray([row.scale_factor for row in rows], dtype=float)
    scale_err = np.asarray([row.scale_factor_err for row in rows], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    axes[0].errorbar(
        x - 0.07,
        data_eff,
        yerr=data_err,
        fmt="o",
        capsize=3,
        label="Data",
    )
    axes[0].errorbar(
        x + 0.07,
        mc_eff,
        yerr=mc_err,
        fmt="s",
        capsize=3,
        label="AAOGEN MC",
    )
    axes[0].set_ylabel("Photon efficiency")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].errorbar(
        x,
        scale,
        yerr=scale_err,
        fmt="o",
        capsize=3,
    )
    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1].set_ylabel(r"$S_\gamma=\epsilon_{\rm data}/\epsilon_{\rm MC}$")
    finite = scale[np.isfinite(scale)]
    axes[1].set_ylim(
        0.0,
        max(1.25, 1.2 * float(np.max(finite)))
        if finite.size
        else 1.25,
    )
    axes[1].grid(alpha=0.25)
    axes[1].set_xticks(x, labels)

    fig.suptitle(title)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_all_period_results(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    periods = [period.label for period in PERIODS]
    categories = ["FT"] + [f"FD S{sector}" for sector in range(1, 7)]
    x = np.arange(len(categories), dtype=float)
    offsets = np.linspace(-0.24, 0.24, len(periods))

    fig, axes = plt.subplots(3, 1, figsize=(14, 14), sharex=True)

    for offset, period_label in zip(offsets, periods):
        selected = [
            row for row in rows
            if row["period_label"] == period_label
        ]
        selected.sort(
            key=lambda row: (
                0 if row["detector"] == "FT" else 1,
                int(row["sector"]),
            )
        )
        if len(selected) != 7:
            continue
        # endif

        axes[0].errorbar(
            x + offset,
            [row["efficiency_data"] for row in selected],
            yerr=[row["efficiency_data_err"] for row in selected],
            fmt="o",
            capsize=2,
            label=period_label,
        )
        axes[1].errorbar(
            x + offset,
            [row["efficiency_mc"] for row in selected],
            yerr=[row["efficiency_mc_err"] for row in selected],
            fmt="o",
            capsize=2,
            label=period_label,
        )
        axes[2].errorbar(
            x + offset,
            [row["scale_factor"] for row in selected],
            yerr=[row["scale_factor_err"] for row in selected],
            fmt="o",
            capsize=2,
            label=period_label,
        )
    # endfor

    axes[0].set_ylabel(r"$\epsilon_{\rm data}$")
    axes[1].set_ylabel(r"$\epsilon_{\rm MC}$")
    axes[2].set_ylabel(r"$S_\gamma$")
    axes[2].axhline(1.0, linestyle="--", linewidth=1.0)

    axes[0].set_ylim(0.0, 1.05)
    axes[1].set_ylim(0.0, 1.05)
    finite_scale = np.asarray(
        [
            row["scale_factor"]
            for row in rows
            if math.isfinite(float(row["scale_factor"]))
        ],
        dtype=float,
    )
    axes[2].set_ylim(
        0.0,
        max(1.25, 1.2 * float(np.max(finite_scale)))
        if finite_scale.size
        else 1.25,
    )

    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend(ncol=3, fontsize=8)
    # endfor
    axes[2].set_xticks(x, categories)

    fig.suptitle("Integrated photon efficiencies and data/MC scale factors")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_rows_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    # endif
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def preflight(periods: Sequence[PeriodConfig]) -> List[Dict[str, object]]:
    manifest: List[Dict[str, object]] = []
    for period in periods:
        entry = asdict(period)
        for role in (
            "epg_data",
            "epgg_data",
            "dvcs_mc",
            "pi0_epg_mc",
            "pi0_epgg_mc",
        ):
            path = getattr(period, role)
            entries, branches = require_tree(path)
            entry[f"{role}_entries"] = entries
            entry[f"{role}_branches"] = branches
            log(f"Preflight {period.label} {role}: {entries:,} entries")
        # endfor
        manifest.append(entry)
    # endfor
    return manifest


def validate_runtime_helpers() -> None:
    required = {
        "read_opportunities": read_opportunities,
        "read_epgg": read_epgg,
        "match_truth_partners": match_truth_partners,
        "emit_cutflow_diagnostics": emit_cutflow_diagnostics,
        "plot_expected_probe_diagnostics": plot_expected_probe_diagnostics,
        "plot_matching_residuals": plot_matching_residuals,
        "binomial_efficiency": binomial_efficiency,
        "uncertainty_ratio": uncertainty_ratio,
    }
    missing = [name for name, value in required.items() if not callable(value)]
    if missing:
        raise RuntimeError(
            "Required runtime helpers are unavailable: " + ", ".join(missing)
        )
    # endif


def main() -> int:
    validate_runtime_helpers()
    args = parse_args()
    periods = selected_periods(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "script": Path(__file__).name,
        "created_unix_time": time.time(),
        "arguments": vars(args),
        "periods": preflight(periods),
    }
    with open(output_dir / "input_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    workers = max(1, min(int(args.workers), MAX_WORKERS, len(periods)))
    log(
        f"PASS/FAIL SIMULTANEOUS FIT: {len(periods)} period(s), "
        f"{workers} worker(s). Directed opportunity: E_tag >= "
        f"{args.tag_E_min:g} GeV, predicted E_probe >= "
        f"{args.probe_E_min:g} GeV. Fit drivers: "
        + ", ".join(args.fit_drivers)
    )

    all_rows: List[Dict[str, object]] = []
    metadata: Dict[str, object] = {}

    if workers == 1:
        for period in periods:
            key, rows, payload = process_period(period, vars(args))
            metadata[key] = payload
            all_rows.extend(rows)
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers
        ) as executor:
            futures = {
                executor.submit(process_period, period, vars(args)): period.key
                for period in periods
            }
            for future in concurrent.futures.as_completed(futures):
                key, rows, payload = future.result()
                metadata[key] = payload
                all_rows.extend(rows)
                log(f"Completed {key}.")
            # endfor
        # endwith
    # endif

    period_order = {period.key: index for index, period in enumerate(PERIODS)}
    all_rows.sort(
        key=lambda row: (
            period_order.get(row["period"], 999),
            0 if row["detector"] == "FT" else 1,
            int(row["sector"]),
        )
    )

    write_rows_csv(
        output_dir / "photon_efficiency_pass_fail_results.csv",
        all_rows,
    )
    with open(
        output_dir / "photon_efficiency_pass_fail_results.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "arguments": vars(args),
                "rows": all_rows,
                "period_metadata": metadata,
            },
            handle,
            indent=2,
        )
    # endwith

    plot_all_period_results(
        output_dir / "all_periods_integrated_efficiencies.png",
        all_rows,
    )
    log(f"Wrote PASS/FAIL photon-efficiency outputs to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
