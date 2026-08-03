#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v38_template_fit_directory_fix.py

Stepwise photon-efficiency study.

The accepted epgamma population is split into events represented by exactly
one reconstructed-photon entry and events represented by exactly two
reconstructed-photon entries. Events represented by three or more photon
entries are rejected in full.

For each run period and accepted multiplicity, the after-cuts data
distributions are fitted simultaneously as a mixture of:

    DVCSGEN epgamma MC
    AAOGEN-as-epgamma pi0 background MC

using the validated shared two-template morphing fitter from
plot_exclusivity_data_dvcs_pi0_mc.py.

The five one-dimensional projections fitted simultaneously are:

    Delta_phi
    theta_gamma_gamma
    pTmiss
    Emiss2
    Mx2_2

A single pi0 fraction is shared by all five projections. Each projection has
independent DVCSGEN and AAOGEN shape-nuisance parameters.

Morphing model:

  * Delta_phi, Emiss2, and Mx2_2:
      additive translation plus symmetric Gaussian broadening;
  * theta_gamma_gamma and pTmiss:
      translation and Gaussian broadening in log(x + epsilon), preserving the
      physical lower boundary;
  * DVCSGEN and AAOGEN receive separate shift and broadening parameters;
  * Gaussian nuisance penalties regularize large shifts and broadenings;
  * the Emiss2 and Mx2_2 morphed DVCSGEN means are required to remain below
    the corresponding AAOGEN means.

The mixture is fitted to the raw histogram counts with a binned Poisson
deviance. Each template is normalized internally, so the fit determines the
composition of the selected data population rather than an absolute MC
normalization.

The five projections are marginals of the same event sample. Their combined
likelihood is therefore a composite likelihood rather than a statistically
independent five-dimensional likelihood. The fitted fraction is meaningful as
a common shape-composition estimator, while its curvature uncertainty does not
include correlations among the five projections.

Data products are written under:

  output/photon_efficiency_study/data/
      shape_comparison/
      template_fits/
          <period>_one_photon_template_fit.png
          <period>_two_photon_template_fit.png
          template_fit_results.csv
          template_fit_results.json
          all_periods_fitted_pi0_fraction.png

The completed AAOGEN MC photon-efficiency machinery remains under:

  output/photon_efficiency_study/mc/

and can be rerun with --run-mc-efficiency.
"""


from __future__ import annotations

import argparse
import importlib.util
import inspect
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
    FitVariable(
        "Delta_phi",
        r"$\Delta\phi$ (rad)",
        120,
        0.5 * math.pi,
        1.5 * math.pi,
    ),
    FitVariable(
        "theta_cm",
        r"$\theta_{p\gamma}^{\rm CM}$ (rad)",
        120,
        1.0,
        math.pi,
    ),
    FitVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$ (deg)",
        120,
        0.0,
        40.0,
    ),
    FitVariable(
        "pTmiss",
        r"$p_T^{\rm miss}$ (GeV)",
        125,
        0.0,
        2.0,
    ),
    FitVariable(
        "Emiss2",
        r"$E_{\rm miss}$ (GeV)",
        120,
        0.0,
        9.0,
    ),
    FitVariable(
        "Mx2",
        r"$M_x^2$ (GeV$^2$)",
        100,
        -0.03,
        0.03,
    ),
    FitVariable(
        "Mx2_2",
        r"$M_{x2}^2$ (GeV$^2$)",
        160,
        0.0,
        16.0,
    ),
)


DATA_SHAPE_VARIABLES: Tuple[FitVariable, ...] = (
    FitVariable(
        "Delta_phi",
        r"$\Delta\phi$ (rad)",
        120,
        0.5 * math.pi,
        1.5 * math.pi,
    ),
    FitVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$ (deg)",
        120,
        0.0,
        40.0,
    ),
    FitVariable(
        "pTmiss",
        r"$p_T^{\rm miss}$ (GeV)",
        125,
        0.0,
        2.0,
    ),
    FitVariable(
        "Emiss2",
        r"$E_{\rm miss}$ (GeV)",
        120,
        0.0,
        9.0,
    ),
    FitVariable(
        "Mx2_2",
        r"$M_{x2}^2$ (GeV$^2$)",
        160,
        0.0,
        16.0,
    ),
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




def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)




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
    "photon_multiplicity_lt_3",
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
        "photon_multiplicity_lt_3": {
            "label": "Event photon-entry multiplicity",
            "condition": (
                "retain only events with one or two reconstructed-photon "
                "entries; reject multiplicity >= 3"
            ),
            "branch": (
                "runnum,evnum for data; rounded reconstructed e/p "
                "kinematics for MC"
            ),
            "applied": True,
            "scope": "hard event-level cut",
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

    entry_stop = (
        min(require_tree(path)[0], int(args.max_events))
        if args.max_events is not None
        else require_tree(path)[0]
    )
    mc_identity_resolved = {
        logical_name: resolved[logical_name]
        for logical_name in identity_logical_names("mc")
    }
    (
        _all_mc_signatures,
        one_mc_signatures,
        two_mc_signatures,
        rejected_mc_signatures,
        multiplicity_cut,
    ) = scan_event_multiplicity(
        path,
        "mc",
        mc_identity_resolved,
        args,
        entry_stop,
    )

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
        "hard_multiplicity_cut": {
            "condition": (
                "retain event multiplicity < 3; reject every event with "
                "multiplicity >= 3"
            ),
            **multiplicity_cut,
        },
    }
    stage_counts = cutflow["stage_counts"]

    seen = 0
    log(f"Reading {role} from {path}")
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
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

        identity_arrays = {
            "e_p": e_p,
            "e_theta": e_theta,
            "e_phi": e_phi,
            "p_p": p_p,
            "p_theta": p_theta,
            "p_phi": p_phi,
        }
        (
            one_multiplicity_keep,
            two_multiplicity_keep,
            _valid_identity,
        ) = multiplicity_entry_masks(
            identity_arrays,
            "mc",
            one_mc_signatures,
            two_mc_signatures,
            args,
        )
        multiplicity_keep = (
            one_multiplicity_keep | two_multiplicity_keep
        )
        stage_counts["photon_multiplicity_lt_3"] += int(
            np.count_nonzero(multiplicity_keep)
        )

        finite = multiplicity_keep & (
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


@dataclass
class MCEfficiencyRow:
    period: str
    period_label: str
    topology: str
    opportunities: int
    reconstructed: int
    missed: int
    efficiency: float
    efficiency_error: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Produce before-cuts and after-cuts normalized epgamma shape "
            "comparisons, with the completed AAOGEN MC efficiency available "
            "as an optional stage."
        )
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[period.key for period in PERIODS],
        help="Restrict processing to one or more run periods.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=8,
        help=f"Number of period workers; hard maximum is {MAX_WORKERS}.",
    )
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional debugging cap per input ROOT tree.",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Common output root containing top-level data/ and mc/ folders.",
    )
    parser.add_argument(
        "--shape-log-y",
        action="store_true",
        help="Also write logarithmic-y normalized shape canvases.",
    )

    parser.add_argument(
        "--data-Emiss-min",
        type=float,
        default=1.0,
        help="After-cuts requirement: Emiss2 greater than this value (GeV).",
    )
    parser.add_argument(
        "--data-Mx2-2-min",
        type=float,
        default=2.0,
        help="After-cuts requirement: Mx2_2 greater than this value (GeV^2).",
    )
    parser.add_argument(
        "--data-pTmiss-min",
        type=float,
        default=0.5,
        help="After-cuts requirement: pTmiss greater than this value (GeV).",
    )

    parser.add_argument(
        "--run-mc-efficiency",
        action="store_true",
        help=(
            "Also rerun the previously completed AAOGEN-only MC photon "
            "efficiency stage under the top-level mc/ directory."
        ),
    )

    # Options retained by the optional AAOGEN MC efficiency stage.
    parser.add_argument("--tag-E-min", type=float, default=0.40)
    parser.add_argument("--probe-E-min", type=float, default=2.00)
    parser.add_argument("--probe-E-max", type=float, default=9.50)
    parser.add_argument("--probe-m2-abs-max", type=float, default=0.10)
    parser.add_argument(
        "--probe-E-minus-p-abs-max",
        type=float,
        default=0.10,
    )
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
    parser.add_argument(
        "--require-fiducial-status-111",
        action="store_true",
    )
    parser.add_argument("--tag-match-angle-max-deg", type=float, default=3.0)
    parser.add_argument(
        "--tag-match-relative-E-max",
        type=float,
        default=0.35,
    )
    parser.add_argument(
        "--probe-match-angle-max-deg",
        type=float,
        default=3.0,
    )
    parser.add_argument(
        "--probe-match-relative-E-max",
        type=float,
        default=0.35,
    )
    parser.add_argument("--mc-signature-decimals", type=int, default=10)
    parser.add_argument("--native-closure-max", type=float, default=0.10)
    parser.add_argument(
        "--native-minimum-photon-E",
        type=float,
        default=0.20,
    )
    parser.add_argument("--found-probe-E-min", type=float, default=2.00)
    parser.add_argument(
        "--require-probe-residual-match",
        action="store_true",
    )
    parser.add_argument("--energy-bins", type=int, default=15)
    parser.add_argument("--theta-bins", type=int, default=13)
    parser.add_argument("--phi-bins", type=int, default=36)
    parser.add_argument(
        "--minimum-bin-opportunities",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--skip-data-template-fits",
        action="store_true",
        help="Skip the shared morphed DVCSGEN+AAOGEN fits.",
    )
    parser.add_argument(
        "--allow-data-fit-failures",
        action="store_true",
        help=(
            "Retain successful categories and serialize failed categories "
            "instead of aborting the run."
        ),
    )
    parser.add_argument(
        "--exclusivity-fitter-script",
        default=str(
            Path(__file__).resolve().parent.parent
            / "plot_exclusivity_data_dvcs_pi0_mc.py"
        ),
        help=(
            "Path to the validated exclusivity template-morphing script. "
            "The parent dvcs_cross_section directory is used by default."
        ),
    )
    parser.add_argument(
        "--fit-max-shift-bins",
        type=float,
        default=12.0,
        help="Maximum absolute additive shift in histogram bins.",
    )
    parser.add_argument(
        "--fit-max-smear-bins",
        type=float,
        default=20.0,
        help="Maximum additional Gaussian broadening in histogram bins.",
    )
    parser.add_argument(
        "--fit-min-counts",
        type=int,
        default=100,
        help="Minimum in-range data count required in every active projection.",
    )
    parser.add_argument(
        "--fit-shift-prior-bins",
        type=float,
        default=4.0,
        help="Gaussian additive-shift prior width in histogram bins.",
    )
    parser.add_argument(
        "--fit-smear-prior-bins",
        type=float,
        default=8.0,
        help="Half-Gaussian broadening-prior width in histogram bins.",
    )
    parser.add_argument(
        "--disable-fit-nuisance-penalties",
        action="store_true",
        help="Disable shift and broadening nuisance penalties.",
    )
    parser.add_argument(
        "--fit-outside-overshoot-penalty",
        type=float,
        default=0.25,
        help="One-sided model-overshoot penalty outside active support.",
    )
    parser.add_argument(
        "--fit-emiss2-mean-order-penalty",
        type=float,
        default=25.0,
        help=(
            "Compatibility argument passed to the validated fitter. The "
            "current fitter applies hard mean ordering for Emiss2 and Mx2_2."
        ),
    )
    parser.add_argument(
        "--fit-dpi",
        type=int,
        default=180,
        help="Resolution of template-fit diagnostic plots.",
    )
    return parser.parse_args()


def binomial_efficiency(
    reconstructed: int,
    opportunities: int,
) -> Tuple[float, float]:
    if opportunities <= 0:
        return math.nan, math.nan
    # endif
    if reconstructed < 0 or reconstructed > opportunities:
        raise ValueError(
            "Invalid binomial population: "
            f"reconstructed={reconstructed}, opportunities={opportunities}"
        )
    # endif

    efficiency = reconstructed / opportunities
    uncertainty = math.sqrt(
        max(efficiency * (1.0 - efficiency) / opportunities, 0.0)
    )
    return efficiency, uncertainty


def topology_mask(
    records: OpportunityRecords,
    topology: str,
) -> np.ndarray:
    if topology == "FT":
        return records.probe_detector == 0
    # endif
    if topology == "FD":
        return records.probe_detector == 1
    # endif
    raise ValueError(f"Unsupported topology: {topology}")


def efficiency_row(
    period: PeriodConfig,
    topology: str,
    mask: np.ndarray,
    matched: np.ndarray,
) -> MCEfficiencyRow:
    opportunities = int(np.count_nonzero(mask))
    reconstructed = int(np.count_nonzero(mask & matched))
    efficiency, uncertainty = binomial_efficiency(
        reconstructed,
        opportunities,
    )
    return MCEfficiencyRow(
        period=period.key,
        period_label=period.label,
        topology=topology,
        opportunities=opportunities,
        reconstructed=reconstructed,
        missed=opportunities - reconstructed,
        efficiency=efficiency,
        efficiency_error=uncertainty,
    )


def binned_efficiency(
    values: np.ndarray,
    selected: np.ndarray,
    matched: np.ndarray,
    edges: np.ndarray,
) -> List[Dict[str, object]]:
    values = np.asarray(values, dtype=float)
    selected = np.asarray(selected, dtype=bool)
    matched = np.asarray(matched, dtype=bool)

    finite = selected & np.isfinite(values)
    denominator, _ = np.histogram(values[finite], bins=edges)
    numerator, _ = np.histogram(values[finite & matched], bins=edges)

    rows: List[Dict[str, object]] = []
    for index in range(edges.size - 1):
        efficiency, uncertainty = binomial_efficiency(
            int(numerator[index]),
            int(denominator[index]),
        )
        rows.append(
            {
                "bin_index": index,
                "bin_low": float(edges[index]),
                "bin_high": float(edges[index + 1]),
                "bin_center": float(
                    0.5 * (edges[index] + edges[index + 1])
                ),
                "opportunities": int(denominator[index]),
                "reconstructed": int(numerator[index]),
                "missed": int(denominator[index] - numerator[index]),
                "efficiency": efficiency,
                "efficiency_error": uncertainty,
            }
        )
    # endfor
    return rows


def write_dict_rows(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    # endif

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def plot_matching_residuals_mc_only(
    path: Path,
    period_label: str,
    match: PartnerMatchResult,
) -> None:
    panels = (
        ("tag_angle_deg", "Tag angular residual (deg)", 0.0, 5.0),
        (
            "tag_relative_E",
            "Tag relative energy residual",
            0.0,
            0.5,
        ),
        (
            "probe_angle_deg",
            "Predicted probe to partner angle (deg)",
            0.0,
            20.0,
        ),
        (
            "probe_relative_E",
            "Predicted probe to partner relative energy",
            0.0,
            1.5,
        ),
        (
            "partner_E",
            "Reconstructed partner energy (GeV)",
            0.0,
            10.0,
        ),
        (
            "solution_closure",
            "Native daughter-solution closure",
            0.0,
            0.10,
        ),
    )

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for axis, (field, label, low, high) in zip(axes.flat, panels):
        values = np.asarray(getattr(match, field), dtype=float)
        values = values[np.isfinite(values)]
        counts, edges = np.histogram(
            values,
            bins=np.linspace(low, high, 101),
        )
        if counts.sum() > 0:
            counts = counts / counts.sum()
        # endif
        centers = 0.5 * (edges[:-1] + edges[1:])
        axis.step(centers, counts, where="mid")
        axis.set_xlabel(label)
        axis.set_ylabel("Normalized entries")
        axis.set_xlim(low, high)
        axis.grid(alpha=0.25)
    # endfor

    fig.suptitle(f"{period_label}: AAOGEN truth-partner matching")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_efficiency_dependence(
    path: Path,
    title: str,
    variable_label: str,
    topology_rows: Mapping[str, Sequence[Mapping[str, object]]],
    minimum_opportunities: int,
) -> None:
    fig, axis = plt.subplots(figsize=(11, 7))

    for topology, rows in topology_rows.items():
        centers = np.asarray(
            [row["bin_center"] for row in rows],
            dtype=float,
        )
        efficiencies = np.asarray(
            [row["efficiency"] for row in rows],
            dtype=float,
        )
        uncertainties = np.asarray(
            [row["efficiency_error"] for row in rows],
            dtype=float,
        )
        support = np.asarray(
            [row["opportunities"] for row in rows],
            dtype=int,
        )
        valid = (
            np.isfinite(efficiencies)
            & np.isfinite(uncertainties)
            & (support >= minimum_opportunities)
        )
        axis.errorbar(
            centers[valid],
            efficiencies[valid],
            yerr=uncertainties[valid],
            fmt="o-",
            capsize=2,
            linewidth=1.0,
            markersize=4,
            label=topology,
        )
    # endfor

    axis.set_xlabel(variable_label)
    axis.set_ylabel(r"$\epsilon_{\mathrm{MC}}$")
    axis.set_ylim(0.0, 1.05)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fd_sector_diagnostics(
    path: Path,
    period_label: str,
    records: OpportunityRecords,
    matched: np.ndarray,
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for sector in range(1, 7):
        mask = (
            (records.probe_detector == 1)
            & (records.probe_sector == sector)
        )
        opportunities = int(np.count_nonzero(mask))
        reconstructed = int(np.count_nonzero(mask & matched))
        efficiency, uncertainty = binomial_efficiency(
            reconstructed,
            opportunities,
        )
        rows.append(
            {
                "sector": sector,
                "opportunities": opportunities,
                "reconstructed": reconstructed,
                "missed": opportunities - reconstructed,
                "efficiency": efficiency,
                "efficiency_error": uncertainty,
            }
        )
    # endfor

    x = np.arange(1, 7)
    efficiency = np.asarray(
        [row["efficiency"] for row in rows],
        dtype=float,
    )
    uncertainty = np.asarray(
        [row["efficiency_error"] for row in rows],
        dtype=float,
    )

    fig, axis = plt.subplots(figsize=(10, 6))
    axis.errorbar(
        x,
        efficiency,
        yerr=uncertainty,
        fmt="o",
        capsize=3,
    )
    axis.set_xticks(x)
    axis.set_xlabel("Predicted FD photon sector")
    axis.set_ylabel(r"$\epsilon_{\mathrm{MC}}$")
    axis.set_ylim(0.0, 1.05)
    axis.grid(alpha=0.25)
    axis.set_title(
        f"{period_label}: FD-sector diagnostic photon efficiency"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return rows



RAW_SHAPE_SAMPLE_INFO: Tuple[Tuple[str, str, str], ...] = (
    ("data", "Data (BH/DVCS + pi0)", "epg_data"),
    ("dvcsgen", "DVCSGEN", "dvcs_mc"),
    ("aaogen", "AAOGEN as epgamma", "pi0_epg_mc"),
)


def resolve_branches_from_keys(
    path: str,
    available_keys: Sequence[str],
    logical_names: Sequence[str],
) -> Dict[str, str]:
    """
    Resolve logical branches from an already opened tree branch list.

    This avoids reopening every ROOT file inside resolve_branches().
    """
    available = set(available_keys)
    resolved: Dict[str, str] = {}
    missing: List[str] = []

    for logical_name in logical_names:
        aliases = ALIASES.get(logical_name, (logical_name,))
        branch = next(
            (candidate for candidate in aliases if candidate in available),
            None,
        )
        if branch is None:
            missing.append(
                f"{logical_name}: tried {', '.join(aliases)}"
            )
        else:
            resolved[logical_name] = branch
        # endif
    # endfor

    if missing:
        raise KeyError(
            f"Missing required branches in {path}:\n  "
            + "\n  ".join(missing)
        )
    # endif
    return resolved


def empty_stage_histograms() -> Dict[str, Dict[str, object]]:
    histograms: Dict[str, Dict[str, object]] = {}
    for variable in DATA_SHAPE_VARIABLES:
        edges = np.linspace(
            variable.low,
            variable.high,
            variable.bins + 1,
        )
        histograms[variable.key] = {
            "edges": edges,
            "counts": np.zeros(variable.bins, dtype=np.int64),
            "selected_entries": 0,
            "in_range_entries": 0,
            "underflow_entries": 0,
            "overflow_entries": 0,
        }
    # endfor
    return histograms


def update_stage_histograms(
    stage_histograms: Dict[str, Dict[str, object]],
    arrays_by_logical_name: Mapping[str, np.ndarray],
    mask: np.ndarray,
) -> None:
    for variable in DATA_SHAPE_VARIABLES:
        selected = np.asarray(
            arrays_by_logical_name[variable.key],
            dtype=float,
        )[mask]
        histogram = stage_histograms[variable.key]
        edges = np.asarray(histogram["edges"], dtype=float)
        counts, _ = np.histogram(selected, bins=edges)

        histogram["counts"] += counts.astype(np.int64)
        histogram["selected_entries"] += int(selected.size)
        histogram["in_range_entries"] += int(np.sum(counts))
        histogram["underflow_entries"] += int(
            np.count_nonzero(selected < variable.low)
        )
        histogram["overflow_entries"] += int(
            np.count_nonzero(selected >= variable.high)
        )
    # endfor


def finalize_stage_histograms(
    stage_histograms: Dict[str, Dict[str, object]],
) -> Dict[str, Dict[str, object]]:
    finalized: Dict[str, Dict[str, object]] = {}

    for variable in DATA_SHAPE_VARIABLES:
        histogram = stage_histograms[variable.key]
        counts = np.asarray(histogram["counts"], dtype=np.int64)
        in_range = int(histogram["in_range_entries"])
        selected = int(histogram["selected_entries"])

        finalized[variable.key] = {
            "edges": np.asarray(histogram["edges"], dtype=float),
            "counts": counts,
            "unit_area": (
                counts.astype(float) / in_range
                if in_range > 0
                else np.zeros_like(counts, dtype=float)
            ),
            "selected_entries": selected,
            "in_range_entries": in_range,
            "underflow_entries": int(
                histogram["underflow_entries"]
            ),
            "overflow_entries": int(
                histogram["overflow_entries"]
            ),
            "in_range_fraction_of_selected": (
                in_range / selected if selected > 0 else math.nan
            ),
            "underflow_fraction_of_selected": (
                int(histogram["underflow_entries"]) / selected
                if selected > 0
                else math.nan
            ),
            "overflow_fraction_of_selected": (
                int(histogram["overflow_entries"]) / selected
                if selected > 0
                else math.nan
            ),
        }
    # endfor

    return finalized


def uint64_mix(values: np.ndarray) -> np.ndarray:
    """
    SplitMix64 finalizer used to build a deterministic 128-bit MC signature.
    """
    values = np.asarray(values, dtype=np.uint64)
    values = values ^ (values >> np.uint64(30))
    values = values * np.uint64(0xBF58476D1CE4E5B9)
    values = values ^ (values >> np.uint64(27))
    values = values * np.uint64(0x94D049BB133111EB)
    values = values ^ (values >> np.uint64(31))
    return values


def data_identity_signature(
    arrays_by_logical_name: Mapping[str, np.ndarray],
) -> Tuple[np.ndarray, np.ndarray]:
    runnum = np.asarray(
        arrays_by_logical_name["runnum"],
        dtype=float,
    )
    eventnum = np.asarray(
        arrays_by_logical_name["eventnum"],
        dtype=float,
    )
    valid = np.isfinite(runnum) & np.isfinite(eventnum)

    signature = np.empty(
        int(np.count_nonzero(valid)),
        dtype=[("runnum", "<i8"), ("eventnum", "<i8")],
    )
    signature["runnum"] = np.rint(runnum[valid]).astype(np.int64)
    signature["eventnum"] = np.rint(eventnum[valid]).astype(np.int64)
    return signature, valid


def mc_identity_signature(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    decimals: int,
) -> Tuple[np.ndarray, np.ndarray]:
    logical_names = (
        "e_p",
        "e_theta",
        "e_phi",
        "p_p",
        "p_theta",
        "p_phi",
    )
    columns = [
        np.asarray(arrays_by_logical_name[name], dtype=float)
        for name in logical_names
    ]

    valid = np.ones(columns[0].size, dtype=bool)
    for column in columns:
        valid &= np.isfinite(column)
    # endfor

    scale = float(10 ** int(decimals))
    quantized = [
        np.rint(column[valid] * scale).astype(np.int64)
        for column in columns
    ]

    count = int(np.count_nonzero(valid))
    h1 = np.full(count, np.uint64(0x243F6A8885A308D3))
    h2 = np.full(count, np.uint64(0x13198A2E03707344))

    for index, column in enumerate(quantized):
        unsigned = column.view(np.uint64)
        mixed_a = uint64_mix(
            unsigned
            + np.uint64(
                (0x9E3779B97F4A7C15 * (index + 1))
                & 0xFFFFFFFFFFFFFFFF
            )
        )
        mixed_b = uint64_mix(
            unsigned
            ^ np.uint64(
                (0xD1B54A32D192ED03 * (index + 1))
                & 0xFFFFFFFFFFFFFFFF
            )
        )
        h1 = uint64_mix(h1 ^ mixed_a)
        h2 = uint64_mix(h2 + mixed_b)
    # endfor

    signature = np.empty(
        count,
        dtype=[("hash1", "<u8"), ("hash2", "<u8")],
    )
    signature["hash1"] = h1
    signature["hash2"] = h2
    return signature, valid


def photon_multiplicity_summary(
    signatures: Sequence[np.ndarray],
    period: PeriodConfig,
    sample_key: str,
    sample_label: str,
    entries_read: int,
    valid_identity_entries: int,
    args: argparse.Namespace,
) -> Dict[str, object]:
    if signatures:
        all_signatures = np.concatenate(signatures)
    else:
        if sample_key == "data":
            all_signatures = np.empty(
                0,
                dtype=[("runnum", "<i8"), ("eventnum", "<i8")],
            )
        else:
            all_signatures = np.empty(
                0,
                dtype=[("hash1", "<u8"), ("hash2", "<u8")],
            )
        # endif
    # endif

    if all_signatures.size > 0:
        _, multiplicities = np.unique(
            all_signatures,
            return_counts=True,
        )
    else:
        multiplicities = np.empty(0, dtype=np.int64)
    # endif

    exact_counts = {
        multiplicity: int(
            np.count_nonzero(multiplicities == multiplicity)
        )
        for multiplicity in range(1, 7)
    }
    above_six = int(np.count_nonzero(multiplicities > 6))
    unique_events = int(multiplicities.size)
    accounted_entries = int(np.sum(multiplicities))

    return {
        "period": period.key,
        "period_label": period.label,
        "sample": sample_key,
        "sample_label": sample_label,
        "identity_method": (
            "exact runnum + evnum"
            if sample_key == "data"
            else (
                "128-bit hash of rounded reconstructed electron/proton "
                f"kinematics ({args.mc_signature_decimals} decimals)"
            )
        ),
        "tree_entries_read": int(entries_read),
        "valid_identity_entries": int(valid_identity_entries),
        "unique_events": unique_events,
        "events_with_1_photon": exact_counts[1],
        "events_with_2_photons": exact_counts[2],
        "events_with_3_photons": exact_counts[3],
        "events_with_4_photons": exact_counts[4],
        "events_with_5_photons": exact_counts[5],
        "events_with_6_photons": exact_counts[6],
        "events_with_more_than_6_photons": above_six,
        "maximum_photons_in_event": (
            int(np.max(multiplicities))
            if multiplicities.size > 0
            else 0
        ),
        "photon_entries_accounted_for": accounted_entries,
    }



def identity_logical_names(sample_key: str) -> Tuple[str, ...]:
    if sample_key == "data":
        return ("runnum", "eventnum")
    # endif
    return (
        "e_p",
        "e_theta",
        "e_phi",
        "p_p",
        "p_theta",
        "p_phi",
    )


def build_identity_signature(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    sample_key: str,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    if sample_key == "data":
        return data_identity_signature(arrays_by_logical_name)
    # endif
    return mc_identity_signature(
        arrays_by_logical_name,
        args.mc_signature_decimals,
    )


def scan_event_multiplicity(
    path: str,
    sample_key: str,
    resolved: Mapping[str, str],
    args: argparse.Namespace,
    entry_stop: int,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Dict[str, int],
]:
    """
    Scan identity branches and classify events by photon-entry multiplicity.

    Returns:
      all unique event signatures,
      exact-one-photon event signatures,
      exact-two-photon event signatures,
      rejected event signatures with multiplicity >= 3,
      entry/event accounting.
    """
    logical_names = identity_logical_names(sample_key)
    expressions = sorted(
        {resolved[logical_name] for logical_name in logical_names}
    )
    signature_chunks: List[np.ndarray] = []
    valid_identity_entries = 0
    entries_scanned = 0

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        chunk_size = len(next(iter(arrays.values()))) if arrays else 0
        entries_scanned += chunk_size

        logical_arrays = {
            logical_name: np.asarray(
                arrays[resolved[logical_name]],
                dtype=float,
            )
            for logical_name in logical_names
        }
        signatures, valid_identity = build_identity_signature(
            logical_arrays,
            sample_key,
            args,
        )
        valid_identity_entries += int(np.count_nonzero(valid_identity))
        signature_chunks.append(signatures)
    # endfor

    if signature_chunks:
        all_entry_signatures = np.concatenate(signature_chunks)
    else:
        if sample_key == "data":
            all_entry_signatures = np.empty(
                0,
                dtype=[("runnum", "<i8"), ("eventnum", "<i8")],
            )
        else:
            all_entry_signatures = np.empty(
                0,
                dtype=[("hash1", "<u8"), ("hash2", "<u8")],
            )
        # endif
    # endif

    if all_entry_signatures.size > 0:
        unique_signatures, multiplicities = np.unique(
            all_entry_signatures,
            return_counts=True,
        )
        one_photon_signatures = unique_signatures[multiplicities == 1]
        two_photon_signatures = unique_signatures[multiplicities == 2]
        rejected_signatures = unique_signatures[multiplicities >= 3]
        rejected_entries = int(
            np.sum(multiplicities[multiplicities >= 3])
        )
    else:
        unique_signatures = all_entry_signatures.copy()
        one_photon_signatures = all_entry_signatures.copy()
        two_photon_signatures = all_entry_signatures.copy()
        rejected_signatures = all_entry_signatures.copy()
        multiplicities = np.empty(0, dtype=np.int64)
        rejected_entries = 0
    # endif

    accounting = {
        "entries_scanned": int(entries_scanned),
        "valid_identity_entries": int(valid_identity_entries),
        "unique_events": int(unique_signatures.size),
        "one_photon_events": int(one_photon_signatures.size),
        "two_photon_events": int(two_photon_signatures.size),
        "rejected_events_multiplicity_ge_3": int(
            rejected_signatures.size
        ),
        "rejected_entries_multiplicity_ge_3": rejected_entries,
        "retained_events_multiplicity_lt_3": int(
            one_photon_signatures.size + two_photon_signatures.size
        ),
        "retained_entries_multiplicity_lt_3": int(
            one_photon_signatures.size
            + 2 * two_photon_signatures.size
        ),
        "maximum_multiplicity": (
            int(np.max(multiplicities))
            if multiplicities.size > 0
            else 0
        ),
    }
    return (
        unique_signatures,
        one_photon_signatures,
        two_photon_signatures,
        rejected_signatures,
        accounting,
    )


def multiplicity_entry_masks(
    arrays_by_logical_name: Mapping[str, np.ndarray],
    sample_key: str,
    one_photon_signatures: np.ndarray,
    two_photon_signatures: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return full-length masks for exact-one and exact-two photon events.
    """
    signatures, valid_identity = build_identity_signature(
        arrays_by_logical_name,
        sample_key,
        args,
    )

    one_valid = np.zeros(signatures.size, dtype=bool)
    two_valid = np.zeros(signatures.size, dtype=bool)

    if signatures.size > 0:
        if one_photon_signatures.size > 0:
            one_valid = np.isin(signatures, one_photon_signatures)
        # endif
        if two_photon_signatures.size > 0:
            two_valid = np.isin(signatures, two_photon_signatures)
        # endif
    # endif

    one_full = np.zeros(valid_identity.size, dtype=bool)
    two_full = np.zeros(valid_identity.size, dtype=bool)
    one_full[valid_identity] = one_valid
    two_full[valid_identity] = two_valid

    return one_full, two_full, valid_identity



def process_data_shape_sample(
    period: PeriodConfig,
    sample_key: str,
    sample_label: str,
    period_attribute: str,
    args_dict: Mapping[str, object],
) -> Tuple[str, str, Dict[str, object]]:
    """
    Determine event multiplicities, then fill exact-one and exact-two shapes.
    """
    args = argparse.Namespace(**args_dict)
    path = getattr(period, period_attribute)

    total_entries, available_keys = require_tree(path)

    shape_logical_names = [
        *(variable.key for variable in DATA_SHAPE_VARIABLES),
        "open_angle_ep2",
    ]
    identity_names = list(identity_logical_names(sample_key))
    logical_names = list(
        dict.fromkeys(shape_logical_names + identity_names)
    )
    resolved = resolve_branches_from_keys(
        path,
        available_keys,
        logical_names,
    )
    expressions = sorted(set(resolved.values()))

    entry_stop = (
        min(total_entries, int(args.max_events))
        if args.max_events is not None
        else total_entries
    )

    (
        all_signatures,
        one_photon_signatures,
        two_photon_signatures,
        rejected_signatures,
        multiplicity_cut,
    ) = scan_event_multiplicity(
        path,
        sample_key,
        resolved,
        args,
        entry_stop,
    )

    multiplicity = photon_multiplicity_summary(
        [all_signatures],
        period,
        sample_key,
        sample_label,
        multiplicity_cut["entries_scanned"],
        multiplicity_cut["valid_identity_entries"],
        args,
    )

    one_before_histograms = empty_stage_histograms()
    one_after_histograms = empty_stage_histograms()
    two_before_histograms = empty_stage_histograms()
    two_after_histograms = empty_stage_histograms()

    entries_read = 0
    valid_identity_entries = 0

    one_finite = 0
    one_open_angle = 0
    one_after_emiss = 0
    one_after_mx2 = 0
    one_after_ptmiss = 0

    two_finite = 0
    two_open_angle = 0
    two_after_emiss = 0
    two_after_mx2 = 0
    two_after_ptmiss = 0

    log(
        f"{period.label} {sample_label}: exact-one events="
        f"{multiplicity_cut['one_photon_events']:,}, exact-two events="
        f"{multiplicity_cut['two_photon_events']:,}, rejected >=3="
        f"{multiplicity_cut['rejected_events_multiplicity_ge_3']:,}"
    )

    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        entry_stop=entry_stop,
        library="np",
    ):
        chunk_size = len(next(iter(arrays.values()))) if arrays else 0
        entries_read += chunk_size

        logical_arrays = {
            logical_name: np.asarray(
                arrays[resolved[logical_name]],
                dtype=float,
            )
            for logical_name in logical_names
        }

        (
            one_multiplicity_mask,
            two_multiplicity_mask,
            valid_identity,
        ) = multiplicity_entry_masks(
            logical_arrays,
            sample_key,
            one_photon_signatures,
            two_photon_signatures,
            args,
        )
        valid_identity_entries += int(np.count_nonzero(valid_identity))

        common_finite = np.isfinite(logical_arrays["open_angle_ep2"])
        for variable in DATA_SHAPE_VARIABLES:
            common_finite &= np.isfinite(logical_arrays[variable.key])
        # endfor

        one_before_mask = one_multiplicity_mask & common_finite
        one_finite += int(np.count_nonzero(one_before_mask))
        one_before_mask &= logical_arrays["open_angle_ep2"] > 5.0
        one_open_angle += int(np.count_nonzero(one_before_mask))

        one_after_mask = one_before_mask.copy()
        one_after_mask &= logical_arrays["Emiss2"] > args.data_Emiss_min
        one_after_emiss += int(np.count_nonzero(one_after_mask))
        one_after_mask &= logical_arrays["Mx2_2"] > args.data_Mx2_2_min
        one_after_mx2 += int(np.count_nonzero(one_after_mask))
        one_after_mask &= logical_arrays["pTmiss"] > args.data_pTmiss_min
        one_after_ptmiss += int(np.count_nonzero(one_after_mask))

        two_before_mask = two_multiplicity_mask & common_finite
        two_finite += int(np.count_nonzero(two_before_mask))
        two_before_mask &= logical_arrays["open_angle_ep2"] > 5.0
        two_open_angle += int(np.count_nonzero(two_before_mask))

        two_after_mask = two_before_mask.copy()
        two_after_mask &= logical_arrays["Emiss2"] > args.data_Emiss_min
        two_after_emiss += int(np.count_nonzero(two_after_mask))
        two_after_mask &= logical_arrays["Mx2_2"] > args.data_Mx2_2_min
        two_after_mx2 += int(np.count_nonzero(two_after_mask))
        two_after_mask &= logical_arrays["pTmiss"] > args.data_pTmiss_min
        two_after_ptmiss += int(np.count_nonzero(two_after_mask))

        update_stage_histograms(
            one_before_histograms,
            logical_arrays,
            one_before_mask,
        )
        update_stage_histograms(
            one_after_histograms,
            logical_arrays,
            one_after_mask,
        )
        update_stage_histograms(
            two_before_histograms,
            logical_arrays,
            two_before_mask,
        )
        update_stage_histograms(
            two_after_histograms,
            logical_arrays,
            two_after_mask,
        )
    # endfor

    if entries_read != multiplicity_cut["entries_scanned"]:
        raise RuntimeError(
            f"{period.label} {sample_label}: multiplicity prepass read "
            f"{multiplicity_cut['entries_scanned']:,} entries but the shape "
            f"pass read {entries_read:,}."
        )
    # endif

    payload = {
        "period": period.key,
        "period_label": period.label,
        "sample": sample_key,
        "sample_label": sample_label,
        "path": path,
        "tree_entries": total_entries,
        "entries_read": entries_read,
        "resolved_branches": resolved,
        "hard_multiplicity_cut": {
            "condition": (
                "retain exact multiplicity 1 or 2; reject multiplicity >= 3"
            ),
            **multiplicity_cut,
        },
        "cutflow": {
            "tree_entries_read": entries_read,
            "valid_identity_entries": valid_identity_entries,
            "one_photon": {
                "finite_required_branches": one_finite,
                "open_angle_ep2_gt_5_deg": one_open_angle,
                "after_Emiss2_gt_min": one_after_emiss,
                "after_Mx2_2_gt_min": one_after_mx2,
                "after_pTmiss_gt_min": one_after_ptmiss,
            },
            "two_photon": {
                "finite_required_branches": two_finite,
                "open_angle_ep2_gt_5_deg": two_open_angle,
                "after_Emiss2_gt_min": two_after_emiss,
                "after_Mx2_2_gt_min": two_after_mx2,
                "after_pTmiss_gt_min": two_after_ptmiss,
            },
        },
        "one_photon": {
            "before_histograms": finalize_stage_histograms(
                one_before_histograms
            ),
            "after_histograms": finalize_stage_histograms(
                one_after_histograms
            ),
        },
        "two_photon": {
            "before_histograms": finalize_stage_histograms(
                two_before_histograms
            ),
            "after_histograms": finalize_stage_histograms(
                two_after_histograms
            ),
        },
        "photon_multiplicity": multiplicity,
    }
    return period.key, sample_key, payload


def plot_before_after_shape_canvas(
    path: Path,
    period_label: str,
    multiplicity_label: str,
    sample_payloads: Mapping[str, Mapping[str, object]],
    multiplicity_key: str,
    args: argparse.Namespace,
) -> None:
    """
    Write one 2x5 canvas for one accepted event-multiplicity category.
    """
    fig, axes = plt.subplots(2, 5, figsize=(24, 10))

    for row_index, stage_name in enumerate(
        ("before_histograms", "after_histograms")
    ):
        stage_label = "Before cuts" if row_index == 0 else "After cuts"

        for column_index, variable in enumerate(DATA_SHAPE_VARIABLES):
            axis = axes[row_index, column_index]

            for sample_key, sample_label, _ in RAW_SHAPE_SAMPLE_INFO:
                histogram = sample_payloads[sample_key][
                    multiplicity_key
                ][stage_name][variable.key]
                edges = np.asarray(histogram["edges"], dtype=float)
                centers = 0.5 * (edges[:-1] + edges[1:])
                normalized = np.asarray(
                    histogram["unit_area"],
                    dtype=float,
                )

                axis.step(
                    centers,
                    normalized,
                    where="mid",
                    linewidth=1.35,
                    label=(
                        f"{sample_label} "
                        f"({histogram['in_range_entries']:,})"
                    ),
                )
            # endfor

            axis.set_xlim(variable.low, variable.high)
            axis.set_xlabel(variable.label)
            axis.set_ylabel("Fraction / bin")
            axis.grid(alpha=0.25)
            axis.legend(fontsize=7)
            axis.set_title(f"{stage_label}: {variable.label}")

            if args.shape_log_y:
                axis.set_yscale("log")
                positive: List[float] = []
                for sample_key, _, _ in RAW_SHAPE_SAMPLE_INFO:
                    values = np.asarray(
                        sample_payloads[sample_key][
                            multiplicity_key
                        ][stage_name][variable.key]["unit_area"],
                        dtype=float,
                    )
                    positive.extend(values[values > 0.0].tolist())
                # endfor
                if positive:
                    axis.set_ylim(
                        bottom=max(min(positive) * 0.5, 1.0e-8)
                    )
                # endif
            # endif
        # endfor
    # endfor

    fig.suptitle(
        f"{period_label}: {multiplicity_label} epgamma shape comparison\n"
        r"Both rows: open_angle_ep2 $>5^\circ$. "
        f"Bottom row additionally: Emiss2 > {args.data_Emiss_min:g} GeV, "
        f"Mx2_2 > {args.data_Mx2_2_min:g} GeV$^2$, "
        f"pTmiss > {args.data_pTmiss_min:g} GeV. "
        "Each sample is independently unit normalized."
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def serializable_histograms(
    histograms: Mapping[str, Mapping[str, object]],
) -> Dict[str, Dict[str, object]]:
    output: Dict[str, Dict[str, object]] = {}
    for variable in DATA_SHAPE_VARIABLES:
        histogram = histograms[variable.key]
        output[variable.key] = {
            key: (
                value.tolist()
                if isinstance(value, np.ndarray)
                else value
            )
            for key, value in histogram.items()
        }
    # endfor
    return output



def multiplicity_categories() -> Tuple[Tuple[str, str], ...]:
    return (
        ("1", "events_with_1_photon"),
        ("2", "events_with_2_photons"),
        ("3", "events_with_3_photons"),
        ("4", "events_with_4_photons"),
        ("5", "events_with_5_photons"),
        ("6", "events_with_6_photons"),
        (">6", "events_with_more_than_6_photons"),
    )


def plot_period_photon_multiplicity(
    path: Path,
    period_label: str,
    rows: Sequence[Mapping[str, object]],
) -> None:
    categories = multiplicity_categories()
    labels = [label for label, _ in categories]
    x = np.arange(len(categories), dtype=float)

    sample_rows = {
        str(row["sample"]): row
        for row in rows
    }
    offsets = {
        "data": -0.24,
        "dvcsgen": 0.0,
        "aaogen": 0.24,
    }
    width = 0.22

    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    for sample_key, sample_label, _ in RAW_SHAPE_SAMPLE_INFO:
        if sample_key not in sample_rows:
            raise RuntimeError(
                f"{period_label}: missing multiplicity row for "
                f"{sample_key}."
            )
        # endif

        row = sample_rows[sample_key]
        totals = np.asarray(
            [int(row[field_name]) for _, field_name in categories],
            dtype=float,
        )
        unique_events = int(row["unique_events"])
        percentages = (
            100.0 * totals / unique_events
            if unique_events > 0
            else np.zeros_like(totals)
        )

        axes[0].bar(
            x + offsets[sample_key],
            totals,
            width=width,
            label=sample_label,
        )
        axes[1].bar(
            x + offsets[sample_key],
            percentages,
            width=width,
            label=sample_label,
        )
    # endfor

    axes[0].set_ylabel("Unique events")
    axes[0].set_title("Total event count")
    axes[0].grid(axis="y", alpha=0.25)
    axes[0].legend()

    axes[1].set_ylabel("Percentage of unique events (%)")
    axes[1].set_xlabel("Number of reconstructed-photon entries")
    axes[1].set_title("Percentage of unique events")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylim(bottom=0.0)
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend()

    fig.suptitle(
        f"{period_label}: epgamma photon-entry multiplicity"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_all_periods_photon_multiplicity(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    categories = multiplicity_categories()
    labels = [label for label, _ in categories]
    x = np.arange(len(categories), dtype=float)

    offsets = {
        "data": -0.24,
        "dvcsgen": 0.0,
        "aaogen": 0.24,
    }
    width = 0.22

    selected_periods = [
        period
        for period in PERIODS
        if any(str(row["period"]) == period.key for row in rows)
    ]

    fig, axes = plt.subplots(
        len(selected_periods),
        2,
        figsize=(18, 4.2 * len(selected_periods)),
        squeeze=False,
    )

    for period_index, period in enumerate(selected_periods):
        period_rows = {
            str(row["sample"]): row
            for row in rows
            if str(row["period"]) == period.key
        }

        total_axis = axes[period_index, 0]
        percent_axis = axes[period_index, 1]

        for sample_key, sample_label, _ in RAW_SHAPE_SAMPLE_INFO:
            if sample_key not in period_rows:
                raise RuntimeError(
                    f"{period.label}: missing multiplicity row for "
                    f"{sample_key}."
                )
            # endif

            row = period_rows[sample_key]
            totals = np.asarray(
                [int(row[field_name]) for _, field_name in categories],
                dtype=float,
            )
            unique_events = int(row["unique_events"])
            percentages = (
                100.0 * totals / unique_events
                if unique_events > 0
                else np.zeros_like(totals)
            )

            total_axis.bar(
                x + offsets[sample_key],
                totals,
                width=width,
                label=sample_label,
            )
            percent_axis.bar(
                x + offsets[sample_key],
                percentages,
                width=width,
                label=sample_label,
            )
        # endfor

        total_axis.set_title(f"{period.label}: total events")
        total_axis.set_ylabel("Unique events")
        total_axis.set_xticks(x, labels)
        total_axis.grid(axis="y", alpha=0.25)

        percent_axis.set_title(f"{period.label}: percentage")
        percent_axis.set_ylabel("Percentage of unique events (%)")
        percent_axis.set_xticks(x, labels)
        percent_axis.set_ylim(bottom=0.0)
        percent_axis.grid(axis="y", alpha=0.25)

        if period_index == 0:
            total_axis.legend(fontsize=8)
            percent_axis.legend(fontsize=8)
        # endif
    # endfor

    axes[-1, 0].set_xlabel("Number of reconstructed-photon entries")
    axes[-1, 1].set_xlabel("Number of reconstructed-photon entries")

    fig.suptitle(
        "epgamma photon-entry multiplicity by run period"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def write_photon_multiplicity_outputs(
    data_root: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    ordered_fields = [
        "period",
        "period_label",
        "sample",
        "sample_label",
        "identity_method",
        "tree_entries_read",
        "valid_identity_entries",
        "unique_events",
        "events_with_1_photon",
        "events_with_2_photons",
        "events_with_3_photons",
        "events_with_4_photons",
        "events_with_5_photons",
        "events_with_6_photons",
        "events_with_more_than_6_photons",
        "maximum_photons_in_event",
        "photon_entries_accounted_for",
    ]

    with open(
        data_root / "photon_multiplicity_summary.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered_fields)
        writer.writeheader()
        writer.writerows(rows)
    # endwith

    with open(
        data_root / "photon_multiplicity_summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(list(rows), handle, indent=2)
    # endwith

    plot_dir = data_root / "photon_multiplicity"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for period in PERIODS:
        period_rows = [
            row
            for row in rows
            if str(row["period"]) == period.key
        ]
        if not period_rows:
            continue
        # endif

        plot_period_photon_multiplicity(
            plot_dir / f"{period.key}_photon_multiplicity.png",
            period.label,
            period_rows,
        )
    # endfor

    plot_all_periods_photon_multiplicity(
        plot_dir / "all_periods_photon_multiplicity.png",
        rows,
    )
    # endwith



_EXCLUSIVITY_FITTER_MODULE = None


def resolve_exclusivity_fitter_path(
    configured_path: str,
) -> Path:
    requested = Path(configured_path).expanduser()
    candidates = (
        requested,
        Path(__file__).resolve().parent.parent
        / "plot_exclusivity_data_dvcs_pi0_mc.py",
        Path(__file__).resolve().parent
        / "plot_exclusivity_data_dvcs_pi0_mc.py",
    )

    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.exists():
            return resolved
        # endif
    # endfor

    attempted = "\n  ".join(
        str(candidate.resolve()) for candidate in candidates
    )
    raise FileNotFoundError(
        "Validated exclusivity fitter was not found. Tried:\n  "
        + attempted
    )


def load_exclusivity_fitter(
    configured_path: str,
):
    global _EXCLUSIVITY_FITTER_MODULE

    resolved_path = resolve_exclusivity_fitter_path(configured_path)
    if (
        _EXCLUSIVITY_FITTER_MODULE is not None
        and Path(_EXCLUSIVITY_FITTER_MODULE.__file__).resolve()
        == resolved_path
    ):
        return _EXCLUSIVITY_FITTER_MODULE
    # endif

    module_name = (
        "_photon_efficiency_validated_exclusivity_fitter"
    )
    spec = importlib.util.spec_from_file_location(
        module_name,
        resolved_path,
    )
    if spec is None or spec.loader is None:
        raise ImportError(
            f"Could not create import specification for {resolved_path}"
        )
    # endif

    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)

    required_symbols = (
        "VariableConfig",
        "TopologyConfig",
        "FitResult",
        "SharedFitSummary",
        "fit_shared_two_templates",
    )
    missing = [
        name for name in required_symbols
        if not hasattr(module, name)
    ]
    if missing:
        raise AttributeError(
            f"{resolved_path} is missing required fitter symbols: "
            + ", ".join(missing)
        )
    # endif

    signature = inspect.signature(module.fit_shared_two_templates)
    required_parameters = (
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
    )
    absent_parameters = [
        name for name in required_parameters
        if name not in signature.parameters
    ]
    if absent_parameters:
        raise TypeError(
            "Validated fit_shared_two_templates() has an incompatible "
            "interface. Missing parameters: "
            + ", ".join(absent_parameters)
        )
    # endif

    module.VARIABLES = tuple(
        module.VariableConfig(
            branch=variable.key,
            label=variable.label,
            bins=variable.bins,
            xmin=variable.low,
            xmax=variable.high,
            aliases=(),
        )
        for variable in DATA_SHAPE_VARIABLES
    )

    _EXCLUSIVITY_FITTER_MODULE = module
    log(
        "Loaded validated exclusivity template fitter from "
        f"{resolved_path}"
    )
    log(
        "Validated fitter signature: "
        f"{signature}"
    )
    return module


def fit_histogram_inputs(
    sample_payloads: Mapping[str, Mapping[str, object]],
    multiplicity_key: str,
) -> Tuple[
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
    Dict[str, np.ndarray],
]:
    def extract(
        sample_key: str,
    ) -> Dict[str, np.ndarray]:
        return {
            variable.key: np.asarray(
                sample_payloads[sample_key][multiplicity_key][
                    "after_histograms"
                ][variable.key]["counts"],
                dtype=np.float64,
            )
            for variable in DATA_SHAPE_VARIABLES
        }

    return extract("data"), extract("dvcsgen"), extract("aaogen")


def fit_result_to_dict(
    result,
) -> Dict[str, object]:
    scalar_fields = (
        "success",
        "message",
        "f_pi0",
        "f_pi0_err",
        "shift",
        "shift_err",
        "sigma_add",
        "sigma_add_err",
        "sigma_right",
        "sigma_right_err",
        "pi0_shift",
        "pi0_shift_err",
        "pi0_sigma_add",
        "pi0_sigma_add_err",
        "pi0_sigma_right",
        "pi0_sigma_right_err",
        "deviance",
        "ndf",
        "data_total",
        "morph_label",
        "excluded_data_counts",
        "excluded_model_counts",
        "excluded_excess_counts",
        "excluded_excess_fraction",
        "fit_region_data_counts",
        "fit_region_model_counts",
        "fit_region_model_to_data",
        "fit_region_closure_difference",
        "full_range_model_counts",
        "full_range_model_to_data",
    )
    return {
        field: getattr(result, field)
        for field in scalar_fields
        if hasattr(result, field)
    }


def shared_fit_to_dict(
    summary,
) -> Dict[str, object]:
    return {
        "success": bool(summary.success),
        "message": str(summary.message),
        "f_pi0": float(summary.f_pi0),
        "f_pi0_err": float(summary.f_pi0_err),
        "deviance": float(summary.deviance),
        "ndf": int(summary.ndf),
        "deviance_per_ndf": (
            float(summary.deviance) / int(summary.ndf)
            if int(summary.ndf) > 0
            else math.nan
        ),
        "fraction_variables": list(summary.fraction_variables),
        "variables": {
            variable_key: fit_result_to_dict(result)
            for variable_key, result in (
                summary.variable_results or {}
            ).items()
        },
    }


def run_one_shared_template_fit(
    period: PeriodConfig,
    multiplicity_key: str,
    multiplicity_label: str,
    sample_payloads: Mapping[str, Mapping[str, object]],
    args: argparse.Namespace,
):
    fitter = load_exclusivity_fitter(
        args.exclusivity_fitter_script
    )
    data_histograms, dvcs_histograms, pi0_histograms = (
        fit_histogram_inputs(
            sample_payloads,
            multiplicity_key,
        )
    )

    insufficient = []
    for variable in DATA_SHAPE_VARIABLES:
        data_count = int(np.sum(data_histograms[variable.key]))
        dvcs_count = int(np.sum(dvcs_histograms[variable.key]))
        pi0_count = int(np.sum(pi0_histograms[variable.key]))
        if (
            data_count < args.fit_min_counts
            or dvcs_count <= 0
            or pi0_count <= 0
        ):
            insufficient.append(
                f"{variable.key}: data={data_count}, "
                f"DVCSGEN={dvcs_count}, AAOGEN={pi0_count}"
            )
        # endif
    # endfor

    if insufficient:
        raise RuntimeError(
            f"{period.label} {multiplicity_label}: insufficient support "
            "for shared template fit:\n  "
            + "\n  ".join(insufficient)
        )
    # endif

    topology = fitter.TopologyConfig(
        key=f"{period.key}_{multiplicity_key}",
        label=f"{period.label}, {multiplicity_label}",
        detector1=-1,
        detector2=-1,
    )

    summary = fitter.fit_shared_two_templates(
        data_histograms=data_histograms,
        dvcs_histograms=dvcs_histograms,
        pi0_histograms=pi0_histograms,
        topology=topology,
        max_shift_bins=args.fit_max_shift_bins,
        max_smear_bins=args.fit_max_smear_bins,
        min_counts=args.fit_min_counts,
        fraction_variable_branches=[
            variable.key for variable in DATA_SHAPE_VARIABLES
        ],
        shift_prior_bins=args.fit_shift_prior_bins,
        smear_prior_bins=args.fit_smear_prior_bins,
        use_nuisance_penalties=(
            not args.disable_fit_nuisance_penalties
        ),
        core_containment=0.90,
        fraction_containment=0.95,
        pi0_support_core_containment=0.90,
        pi0_support_fraction_containment=0.95,
        pi0_core_calibration=None,
        outside_overshoot_penalty_weight=(
            args.fit_outside_overshoot_penalty
        ),
        emiss2_mean_order_penalty_weight=(
            args.fit_emiss2_mean_order_penalty
        ),
    )

    if not summary.success:
        raise RuntimeError(
            f"{period.label} {multiplicity_label}: shared template fit "
            f"failed: {summary.message}"
        )
    # endif

    return summary, data_histograms


def plot_shared_template_fit(
    path: Path,
    period_label: str,
    multiplicity_label: str,
    summary,
    data_histograms: Mapping[str, np.ndarray],
    args: argparse.Namespace,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        2,
        5,
        figsize=(24, 10),
        sharex="col",
        gridspec_kw={"height_ratios": [3.0, 1.2]},
    )

    for column, variable in enumerate(DATA_SHAPE_VARIABLES):
        fit_result = summary.variable_results[variable.key]
        if not fit_result.success:
            raise RuntimeError(
                f"{period_label} {multiplicity_label} "
                f"{variable.key}: variable fit failed: "
                f"{fit_result.message}"
            )
        # endif

        top = axes[0, column]
        bottom = axes[1, column]

        data = np.asarray(
            data_histograms[variable.key],
            dtype=float,
        )
        edges = np.linspace(
            variable.low,
            variable.high,
            variable.bins + 1,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        model = np.asarray(
            fit_result.model_counts,
            dtype=float,
        )
        dvcs_component = np.asarray(
            fit_result.dvcs_component_counts,
            dtype=float,
        )
        pi0_component = np.asarray(
            fit_result.pi0_component_counts,
            dtype=float,
        )

        top.errorbar(
            centers,
            data,
            yerr=np.sqrt(np.maximum(data, 1.0)),
            fmt="o",
            markersize=2.8,
            linewidth=0.8,
            capsize=1.5,
            label="Data",
        )
        top.step(
            centers,
            model,
            where="mid",
            linewidth=1.8,
            label="Total fit",
        )
        top.step(
            centers,
            dvcs_component,
            where="mid",
            linewidth=1.3,
            linestyle="--",
            label="Morphed DVCSGEN",
        )
        top.step(
            centers,
            pi0_component,
            where="mid",
            linewidth=1.3,
            linestyle=":",
            label=r"Morphed AAOGEN $\pi^0$",
        )
        top.set_ylabel("Entries / bin")
        top.set_title(variable.label)
        top.grid(alpha=0.25)
        top.legend(fontsize=7)

        ratio = np.full_like(data, np.nan, dtype=float)
        ratio_error = np.full_like(data, np.nan, dtype=float)
        valid = model > 0.0
        ratio[valid] = data[valid] / model[valid]
        ratio_error[valid] = (
            np.sqrt(np.maximum(data[valid], 1.0))
            / model[valid]
        )

        bottom.errorbar(
            centers[valid],
            ratio[valid],
            yerr=ratio_error[valid],
            fmt="o",
            markersize=2.5,
            linewidth=0.7,
            capsize=1.2,
        )
        bottom.axhline(1.0, linewidth=1.0)
        bottom.set_ylim(0.0, 2.0)
        bottom.set_xlabel(variable.label)
        bottom.set_ylabel("Data / fit")
        bottom.grid(alpha=0.25)

        if args.shape_log_y:
            top.set_yscale("log")
        # endif
    # endfor

    fig.suptitle(
        f"{period_label}: {multiplicity_label} shared morphed-template fit\n"
        rf"$f_{{\pi^0}}={summary.f_pi0:.4f}"
        rf"\pm{summary.f_pi0_err:.4f}$, "
        f"deviance/ndf="
        f"{summary.deviance:.1f}/{summary.ndf}"
        + (
            f"={summary.deviance / summary.ndf:.2f}"
            if summary.ndf > 0
            else ""
        )
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(path, dpi=args.fit_dpi)
    plt.close(fig)


def flatten_shared_fit_row(
    period: PeriodConfig,
    multiplicity_key: str,
    multiplicity_label: str,
    summary,
) -> Dict[str, object]:
    row: Dict[str, object] = {
        "period": period.key,
        "period_label": period.label,
        "multiplicity": multiplicity_key,
        "multiplicity_label": multiplicity_label,
        "fit_success": bool(summary.success),
        "fit_message": str(summary.message),
        "f_pi0": float(summary.f_pi0),
        "f_pi0_error": float(summary.f_pi0_err),
        "deviance": float(summary.deviance),
        "ndf": int(summary.ndf),
        "deviance_per_ndf": (
            float(summary.deviance) / int(summary.ndf)
            if int(summary.ndf) > 0
            else math.nan
        ),
    }

    for variable in DATA_SHAPE_VARIABLES:
        result = summary.variable_results[variable.key]
        prefix = variable.key
        for field in (
            "shift",
            "shift_err",
            "sigma_add",
            "sigma_add_err",
            "sigma_right",
            "sigma_right_err",
            "pi0_shift",
            "pi0_shift_err",
            "pi0_sigma_add",
            "pi0_sigma_add_err",
            "pi0_sigma_right",
            "pi0_sigma_right_err",
            "deviance",
            "ndf",
            "morph_label",
            "full_range_model_to_data",
        ):
            row[f"{prefix}_{field}"] = getattr(result, field)
        # endfor
    # endfor

    return row


def plot_all_period_fitted_pi0_fractions(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    period_labels = [period.label for period in PERIODS]
    x = np.arange(len(period_labels), dtype=float)
    offsets = {
        "one_photon": -0.08,
        "two_photon": 0.08,
    }
    display_labels = {
        "one_photon": "One-photon events",
        "two_photon": "Two-photon events",
    }

    fig, axis = plt.subplots(figsize=(12, 7))

    for multiplicity_key in ("one_photon", "two_photon"):
        selected = []
        for period in PERIODS:
            match = next(
                (
                    row for row in rows
                    if row["period"] == period.key
                    and row["multiplicity"] == multiplicity_key
                    and bool(row["fit_success"])
                ),
                None,
            )
            selected.append(match)
        # endfor

        values = np.asarray(
            [
                row["f_pi0"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )
        errors = np.asarray(
            [
                row["f_pi0_error"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )

        axis.errorbar(
            x + offsets[multiplicity_key],
            values,
            yerr=errors,
            fmt="o",
            capsize=3,
            label=display_labels[multiplicity_key],
        )
    # endfor

    axis.set_xticks(x, period_labels, rotation=20, ha="right")
    axis.set_ylabel(r"Fitted $\pi^0$ fraction")
    axis.set_ylim(0.0, 1.0)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(
        "Shared five-projection morphed-template fits"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_template_fit_outputs(
    fit_dir: Path,
    fit_rows: Sequence[Mapping[str, object]],
    fit_payloads: Mapping[str, object],
) -> None:
    fit_dir.mkdir(parents=True, exist_ok=True)

    if fit_rows:
        fieldnames = list(fit_rows[0].keys())
        with open(
            fit_dir / "template_fit_results.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
            )
            writer.writeheader()
            writer.writerows(fit_rows)
        # endwith
    # endif

    with open(
        fit_dir / "template_fit_results.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "fit_model": {
                    "shared_fraction": True,
                    "fraction_variables": [
                        variable.key
                        for variable in DATA_SHAPE_VARIABLES
                    ],
                    "templates": [
                        "DVCSGEN",
                        "AAOGEN-as-epgamma pi0",
                    ],
                    "data_population": "after_cuts",
                    "multiplicity_categories": [
                        "one_photon",
                        "two_photon",
                    ],
                    "projection_likelihood_note": (
                        "The five one-dimensional projections contain the "
                        "same events and are combined as a composite "
                        "likelihood; projection correlations are not included "
                        "in the curvature uncertainty."
                    ),
                },
                "rows": list(fit_rows),
                "fits": fit_payloads,
            },
            handle,
            indent=2,
        )
    # endwith

    if fit_rows:
        plot_all_period_fitted_pi0_fractions(
            fit_dir / "all_periods_fitted_pi0_fraction.png",
            fit_rows,
        )
    # endif



def run_data_shape_stage(
    periods: Sequence[PeriodConfig],
    args: argparse.Namespace,
    workers: int,
) -> Dict[str, object]:
    """
    Run 15 independent file tasks with no nested multiprocessing.
    """
    data_root = Path(args.output_dir) / "data"
    shape_dir = data_root / "shape_comparison"
    shape_dir.mkdir(parents=True, exist_ok=True)

    tasks = [
        (period, sample_key, sample_label, period_attribute)
        for period in periods
        for sample_key, sample_label, period_attribute
        in RAW_SHAPE_SAMPLE_INFO
    ]
    file_workers = max(
        1,
        min(int(workers), MAX_WORKERS, len(tasks)),
    )

    log(
        f"DATA SHAPE + MULTIPLICITY STAGE: {len(tasks)} files, "
        f"{file_workers} worker(s). Each ROOT file is read once."
    )

    results: Dict[str, Dict[str, Dict[str, object]]] = {
        period.key: {} for period in periods
    }

    if file_workers == 1:
        for period, sample_key, sample_label, period_attribute in tasks:
            key, returned_sample, payload = process_data_shape_sample(
                period,
                sample_key,
                sample_label,
                period_attribute,
                vars(args),
            )
            results[key][returned_sample] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=file_workers
        ) as executor:
            future_map = {
                executor.submit(
                    process_data_shape_sample,
                    period,
                    sample_key,
                    sample_label,
                    period_attribute,
                    vars(args),
                ): (period.key, sample_key)
                for (
                    period,
                    sample_key,
                    sample_label,
                    period_attribute,
                ) in tasks
            }
            for future in concurrent.futures.as_completed(future_map):
                key, returned_sample, payload = future.result()
                results[key][returned_sample] = payload
                log(
                    f"Completed {payload['period_label']} "
                    f"{payload['sample_label']}."
                )
            # endfor
        # endwith
    # endif

    period_order = {
        period.key: index for index, period in enumerate(PERIODS)
    }
    sample_order = {
        sample_key: index
        for index, (sample_key, _, _) in enumerate(
            RAW_SHAPE_SAMPLE_INFO
        )
    }

    multiplicity_rows: List[Dict[str, object]] = []
    audit: Dict[str, object] = {}
    fit_rows: List[Dict[str, object]] = []
    fit_payloads: Dict[str, object] = {}
    fit_dir = data_root / "template_fits"
    if not args.skip_data_template_fits:
        fit_dir.mkdir(parents=True, exist_ok=True)
    # endif

    for period in periods:
        sample_payloads = results[period.key]
        missing_samples = [
            sample_key
            for sample_key, _, _ in RAW_SHAPE_SAMPLE_INFO
            if sample_key not in sample_payloads
        ]
        if missing_samples:
            raise RuntimeError(
                f"{period.label}: missing processed samples: "
                + ", ".join(missing_samples)
            )
        # endif

        one_photon_plot_path = (
            shape_dir
            / f"{period.key}_one_photon_before_after_cuts.png"
        )
        two_photon_plot_path = (
            shape_dir
            / f"{period.key}_two_photon_before_after_cuts.png"
        )

        plot_before_after_shape_canvas(
            one_photon_plot_path,
            period.label,
            "exactly one reconstructed-photon entry",
            sample_payloads,
            "one_photon",
            args,
        )
        plot_before_after_shape_canvas(
            two_photon_plot_path,
            period.label,
            "exactly two reconstructed-photon entries",
            sample_payloads,
            "two_photon",
            args,
        )

        period_fit_payload: Dict[str, object] = {}
        if not args.skip_data_template_fits:
            fit_categories = (
                (
                    "one_photon",
                    "exactly one reconstructed-photon entry",
                ),
                (
                    "two_photon",
                    "exactly two reconstructed-photon entries",
                ),
            )

            for multiplicity_key, multiplicity_label in fit_categories:
                try:
                    summary, data_histograms = (
                        run_one_shared_template_fit(
                            period,
                            multiplicity_key,
                            multiplicity_label,
                            sample_payloads,
                            args,
                        )
                    )

                    fit_plot_path = (
                        fit_dir
                        / f"{period.key}_{multiplicity_key}"
                        "_template_fit.png"
                    )
                    plot_shared_template_fit(
                        fit_plot_path,
                        period.label,
                        multiplicity_label,
                        summary,
                        data_histograms,
                        args,
                    )

                    fit_row = flatten_shared_fit_row(
                        period,
                        multiplicity_key,
                        multiplicity_label,
                        summary,
                    )
                    fit_rows.append(fit_row)
                    period_fit_payload[multiplicity_key] = {
                        "plot": str(fit_plot_path),
                        **shared_fit_to_dict(summary),
                    }

                    log(
                        f"{period.label} {multiplicity_label}: "
                        f"f_pi0={summary.f_pi0:.6f} +/- "
                        f"{summary.f_pi0_err:.6f}, "
                        f"deviance/ndf={summary.deviance:.2f}/"
                        f"{summary.ndf}"
                    )
                except Exception as exc:
                    failure = {
                        "success": False,
                        "message": str(exc),
                    }
                    period_fit_payload[multiplicity_key] = failure
                    if not args.allow_data_fit_failures:
                        raise
                    # endif
                    log(
                        f"WARNING: {period.label} {multiplicity_label} "
                        f"template fit failed: {exc}"
                    )
                # endtry
            # endfor
        # endif

        fit_payloads[period.key] = period_fit_payload

        period_audit: Dict[str, object] = {
            "plots": {
                "one_photon": str(one_photon_plot_path),
                "two_photon": str(two_photon_plot_path),
            },
            "template_fits": period_fit_payload,
            "samples": {},
        }
        for sample_key, _, _ in RAW_SHAPE_SAMPLE_INFO:
            payload = sample_payloads[sample_key]
            multiplicity_rows.append(
                payload["photon_multiplicity"]
            )
            period_audit["samples"][sample_key] = {
                "sample_label": payload["sample_label"],
                "path": payload["path"],
                "tree_entries": payload["tree_entries"],
                "entries_read": payload["entries_read"],
                "resolved_branches": payload["resolved_branches"],
                "hard_multiplicity_cut": payload[
                    "hard_multiplicity_cut"
                ],
                "cutflow": payload["cutflow"],
                "photon_multiplicity": payload[
                    "photon_multiplicity"
                ],
                "one_photon": {
                    "before_histograms": serializable_histograms(
                        payload["one_photon"]["before_histograms"]
                    ),
                    "after_histograms": serializable_histograms(
                        payload["one_photon"]["after_histograms"]
                    ),
                },
                "two_photon": {
                    "before_histograms": serializable_histograms(
                        payload["two_photon"]["before_histograms"]
                    ),
                    "after_histograms": serializable_histograms(
                        payload["two_photon"]["after_histograms"]
                    ),
                },
            }
        # endfor
        audit[period.key] = period_audit
    # endfor

    multiplicity_rows.sort(
        key=lambda row: (
            period_order.get(str(row["period"]), 999),
            sample_order.get(str(row["sample"]), 999),
        )
    )
    write_photon_multiplicity_outputs(
        data_root,
        multiplicity_rows,
    )

    if not args.skip_data_template_fits:
        write_template_fit_outputs(
            fit_dir,
            fit_rows,
            fit_payloads,
        )
    # endif

    with open(
        data_root / "shape_comparison_audit.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(audit, handle, indent=2)
    # endwith

    for row in multiplicity_rows:
        log(
            f"{row['period_label']} {row['sample_label']}: "
            f"events by photon multiplicity "
            f"1={row['events_with_1_photon']:,}, "
            f"2={row['events_with_2_photons']:,}, "
            f"3={row['events_with_3_photons']:,}, "
            f"4={row['events_with_4_photons']:,}, "
            f"5={row['events_with_5_photons']:,}, "
            f"6={row['events_with_6_photons']:,}, "
            f">6={row['events_with_more_than_6_photons']:,}"
        )
    # endfor

    return {
        "periods": audit,
        "photon_multiplicity_rows": multiplicity_rows,
        "template_fit_rows": fit_rows,
        "template_fit_payloads": fit_payloads,
    }



def process_period_mc_only(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = (
        Path(args.output_dir)
        / "mc"
        / period.key
    )
    period_dir.mkdir(parents=True, exist_ok=True)

    opportunities, _, cutflow = read_opportunities(
        period.pi0_epg_mc,
        period.beam_energy_GeV,
        "AAOGEN epgamma opportunities",
        args,
        deduplicate=True,
    )
    emit_cutflow_diagnostics(
        period_dir,
        period.label,
        {"aaogen_epgamma": cutflow},
    )

    epgg, epgg_diagnostics = read_epgg(
        period.pi0_epgg_mc,
        period.beam_energy_GeV,
        "AAOGEN",
        args,
    )
    match = match_truth_partners(
        opportunities,
        epgg,
        "mc",
        args,
    )

    summary = match.summary
    log(
        f"{period.label} AAOGEN: opportunities="
        f"{summary['opportunities']:,}, identity="
        f"{summary['identity_matches']:,}, valid daughter solution="
        f"{summary['valid_solution_opportunities']:,}, tag matched="
        f"{summary['tag_matched']:,}, partner above "
        f"{args.found_probe_E_min:g} GeV="
        f"{summary['partner_above_threshold']:,}, reconstructed="
        f"{summary['matched_opportunities']:,}"
    )

    plot_expected_probe_diagnostics(
        period_dir / "selected_opportunity_kinematics.png",
        period.label,
        {"AAOGEN": opportunities},
    )
    plot_matching_residuals_mc_only(
        period_dir / "matching_residuals.png",
        period.label,
        match,
    )

    official_rows: List[MCEfficiencyRow] = []
    for topology in ("FT", "FD"):
        mask = topology_mask(opportunities, topology)
        row = efficiency_row(
            period,
            topology,
            mask,
            match.matched,
        )
        official_rows.append(row)
        log(
            f"{period.label} {topology}: "
            f"N_reconstructed={row.reconstructed:,}, "
            f"N_opportunities={row.opportunities:,}, "
            f"epsilon_MC={row.efficiency:.6f} +/- "
            f"{row.efficiency_error:.6f}"
        )
    # endfor

    official_dicts = [asdict(row) for row in official_rows]
    write_dict_rows(
        period_dir / "integrated_mc_efficiency.csv",
        official_dicts,
    )

    topology_masks = {
        topology: topology_mask(opportunities, topology)
        for topology in ("FT", "FD")
    }

    energy_edges = np.linspace(
        args.probe_E_min,
        args.probe_E_max,
        args.energy_bins + 1,
    )
    theta_edges = np.linspace(
        args.ft_theta_min,
        args.fd_theta_max,
        args.theta_bins + 1,
    )
    phi_edges = np.linspace(
        0.0,
        2.0 * math.pi,
        args.phi_bins + 1,
    )

    dependence_specs = {
        "probe_energy": (
            opportunities.probe_E,
            energy_edges,
            "Predicted probe energy (GeV)",
        ),
        "probe_theta": (
            np.degrees(opportunities.probe_theta),
            theta_edges,
            "Predicted probe theta (deg)",
        ),
        "probe_phi": (
            np.mod(opportunities.probe_phi, 2.0 * math.pi),
            phi_edges,
            "Predicted probe phi (rad)",
        ),
    }

    binned_outputs: Dict[str, Dict[str, List[Dict[str, object]]]] = {}
    for name, (values, edges, label) in dependence_specs.items():
        topology_rows: Dict[str, List[Dict[str, object]]] = {}
        for topology, mask in topology_masks.items():
            rows = binned_efficiency(
                values,
                mask,
                match.matched,
                edges,
            )
            topology_rows[topology] = rows
            write_dict_rows(
                period_dir / f"mc_efficiency_vs_{name}_{topology.lower()}.csv",
                rows,
            )
        # endfor

        binned_outputs[name] = topology_rows
        plot_efficiency_dependence(
            period_dir / f"mc_efficiency_vs_{name}.png",
            f"{period.label}: AAOGEN photon efficiency",
            label,
            topology_rows,
            args.minimum_bin_opportunities,
        )
    # endfor

    sector_rows = plot_fd_sector_diagnostics(
        period_dir / "fd_sector_diagnostic_efficiency.png",
        period.label,
        opportunities,
        match.matched,
    )
    write_dict_rows(
        period_dir / "fd_sector_diagnostic_efficiency.csv",
        sector_rows,
    )

    metadata = {
        "period": {
            "key": period.key,
            "label": period.label,
            "beam_energy_GeV": period.beam_energy_GeV,
            "pi0_epg_mc": period.pi0_epg_mc,
            "pi0_epgg_mc": period.pi0_epgg_mc,
        },
        "arguments": vars(args),
        "opportunity_cutflow": cutflow,
        "epgammagamma_reconstruction": epgg_diagnostics,
        "matching": summary,
        "integrated_efficiencies": official_dicts,
        "fd_sector_diagnostics": sector_rows,
        "binned_efficiencies": binned_outputs,
    }
    with open(
        period_dir / "mc_efficiency_metadata.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(metadata, handle, indent=2)
    # endwith

    return period.key, official_dicts, metadata


def preflight_enabled_stages(
    periods: Sequence[PeriodConfig],
    args: argparse.Namespace,
) -> List[Dict[str, object]]:
    manifest: List[Dict[str, object]] = []

    for period in periods:
        period_entry: Dict[str, object] = {
            "key": period.key,
            "label": period.label,
            "beam_energy_GeV": period.beam_energy_GeV,
        }

        roles = ["epg_data", "dvcs_mc", "pi0_epg_mc"]
        if args.run_mc_efficiency:
            roles.append("pi0_epgg_mc")
        # endif

        for role in roles:
            path = getattr(period, role)
            entries, branches = require_tree(path)
            period_entry[role] = path
            period_entry[f"{role}_entries"] = entries
            period_entry[f"{role}_branches"] = branches
            log(
                f"Preflight {period.label} {role}: {entries:,} entries"
            )
        # endfor

        manifest.append(period_entry)
    # endfor

    return manifest



def plot_all_period_integrated(
    path: Path,
    rows: Sequence[Mapping[str, object]],
) -> None:
    period_labels = [period.label for period in PERIODS]
    x = np.arange(len(period_labels), dtype=float)
    offsets = {"FT": -0.08, "FD": 0.08}

    fig, axis = plt.subplots(figsize=(12, 7))
    for topology in ("FT", "FD"):
        selected = []
        for label in period_labels:
            match = next(
                (
                    row for row in rows
                    if row["period_label"] == label
                    and row["topology"] == topology
                ),
                None,
            )
            selected.append(match)
        # endfor

        efficiencies = np.asarray(
            [
                row["efficiency"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )
        errors = np.asarray(
            [
                row["efficiency_error"] if row is not None else math.nan
                for row in selected
            ],
            dtype=float,
        )
        axis.errorbar(
            x + offsets[topology],
            efficiencies,
            yerr=errors,
            fmt="o",
            capsize=3,
            label=topology,
        )
    # endfor

    axis.set_xticks(x, period_labels, rotation=20, ha="right")
    axis.set_ylabel(r"$\epsilon_{\mathrm{MC}}$")
    axis.set_ylim(0.0, 1.05)
    axis.grid(alpha=0.25)
    axis.legend()
    axis.set_title(
        "AAOGEN high-energy photon reconstruction efficiency"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    periods = selected_periods(args)
    output_root = Path(args.output_dir)
    data_root = output_root / "data"
    mc_root = output_root / "mc"

    output_root.mkdir(parents=True, exist_ok=True)
    data_root.mkdir(parents=True, exist_ok=True)
    mc_root.mkdir(parents=True, exist_ok=True)

    workers = max(
        1,
        min(int(args.workers), MAX_WORKERS),
    )

    fitter_preflight: Dict[str, object] = {
        "enabled": not args.skip_data_template_fits,
    }
    if not args.skip_data_template_fits:
        fitter = load_exclusivity_fitter(
            args.exclusivity_fitter_script
        )
        fitter_preflight.update(
            {
                "resolved_path": str(
                    Path(fitter.__file__).resolve()
                ),
                "fit_signature": str(
                    inspect.signature(
                        fitter.fit_shared_two_templates
                    )
                ),
                "variables": [
                    {
                        "branch": variable.branch,
                        "label": variable.label,
                        "bins": variable.bins,
                        "xmin": variable.xmin,
                        "xmax": variable.xmax,
                    }
                    for variable in fitter.VARIABLES
                ],
            }
        )
    # endif

    manifest = {
        "script": Path(__file__).name,
        "mode": (
            "Fast data 2x5 shape comparison and photon multiplicity "
            "diagnostics with preserved optional AAOGEN MC efficiency"
        ),
        "created_unix_time": time.time(),
        "arguments": vars(args),
        "fitter_preflight": fitter_preflight,
        "periods": preflight_enabled_stages(periods, args),
    }
    with open(
        output_root / "input_manifest.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    data_metadata = run_data_shape_stage(
        periods,
        args,
        workers,
    )

    mc_rows: List[Dict[str, object]] = []
    mc_metadata: Dict[str, object] = {}

    if args.run_mc_efficiency:
        mc_workers = max(
            1,
            min(workers, len(periods)),
        )
        log(
            f"OPTIONAL AAOGEN MC EFFICIENCY STAGE: "
            f"{len(periods)} period(s), {mc_workers} worker(s)."
        )

        if mc_workers == 1:
            for period in periods:
                key, rows, payload = process_period_mc_only(
                    period,
                    vars(args),
                )
                mc_metadata[key] = payload
                mc_rows.extend(rows)
            # endfor
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=mc_workers
            ) as executor:
                future_map = {
                    executor.submit(
                        process_period_mc_only,
                        period,
                        vars(args),
                    ): period.key
                    for period in periods
                }
                for future in concurrent.futures.as_completed(future_map):
                    key, rows, payload = future.result()
                    mc_metadata[key] = payload
                    mc_rows.extend(rows)
                # endfor
            # endwith
        # endif

        period_order = {
            period.key: index
            for index, period in enumerate(PERIODS)
        }
        topology_order = {"FT": 0, "FD": 1}
        mc_rows.sort(
            key=lambda row: (
                period_order.get(str(row["period"]), 999),
                topology_order.get(str(row["topology"]), 999),
            )
        )

        write_dict_rows(
            mc_root / "photon_efficiency_mc.csv",
            mc_rows,
        )
        with open(
            mc_root / "photon_efficiency_mc.json",
            "w",
            encoding="utf-8",
        ) as handle:
            json.dump(
                {
                    "mode": "AAOGEN MC only",
                    "arguments": vars(args),
                    "rows": mc_rows,
                    "period_metadata": mc_metadata,
                },
                handle,
                indent=2,
            )
        # endwith

        plot_all_period_integrated(
            mc_root / "all_periods_integrated_mc_efficiency.png",
            mc_rows,
        )
    # endif

    with open(
        output_root / "stepwise_study_summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {
                "arguments": vars(args),
                "data_shape_and_multiplicity_stage": data_metadata,
                "mc_efficiency_rows": mc_rows,
                "mc_efficiency_metadata": mc_metadata,
            },
            handle,
            indent=2,
        )
    # endwith

    log(f"Wrote photon-efficiency study outputs to {output_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
