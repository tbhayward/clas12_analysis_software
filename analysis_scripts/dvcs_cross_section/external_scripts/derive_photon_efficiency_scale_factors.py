#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v8.py

Exploratory RGA photon-efficiency tag-and-probe study for the DVCS analysis.

The script measures epsilon_b = N_pass,b/(N_pass,b+N_fail,b) separately in
data and AAOGEN MC, where b is a predicted probe-photon bin. Native eppi0
events provide two directed passing probes. One-photon epgamma candidates
provide the failed-probe population after a simultaneous DVCS/BH plus pi0
template decomposition.

The failed-probe decomposition follows the stable exclusivity-study strategy:
Delta_phi, theta_gamma_gamma and pTmiss are fit simultaneously with one shared
pi0 fraction. Detector-category nuisance morphologies are calibrated on the
aggregate FT or FD-sector sample and then fixed in the differential E_gamma and
theta_gamma fits. Delta_phi uses additive shift/smearing; the positive-definite
theta_gamma_gamma and pTmiss projections use log-space shift/smearing.

The production extraction is deliberately integrated over energy and polar angle,
retaining only FT and FD sectors 1--6. A detailed passing-sample audit records
every event-level and directed-probe rejection stage for data and AAOGEN MC.
Native photon directions are recovered with a cone-constrained analytic solver
that does not use gamma_phi1/2, since those are Trento rather than lab angles. Processing is parallelized by run period with a hard maximum of seven
workers.

The script derives S_gamma,b = epsilon_data,b/epsilon_MC,b. It does not yet
propagate that scale factor into DVCS acceptance or pi0 migration.
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

try:
    import uproot
except ImportError as exc:
    raise SystemExit("This script requires uproot: python3 -m pip install uproot") from exc
#endif

try:
    from scipy.optimize import minimize
    from scipy.ndimage import gaussian_filter1d
except ImportError as exc:
    raise SystemExit("This script requires scipy: python3 -m pip install scipy") from exc
#endif


TREE_NAME = "PhysicsEvents"
DEFAULT_OUTPUT_DIR = "output/photon_efficiency_study"
DEFAULT_STEP_SIZE = "200 MB"
PROTON_MASS_GEV = 0.9382720813
ELECTRON_MASS_GEV = 0.00051099895
MAX_WORKERS = 7


@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy_GeV: float
    epg_data: str
    dvcs_mc: str
    pi0_as_epg_mc: str
    epgg_data: str
    epgg_mc: str


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        "fa18_inb", "Fa18 Inb", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_50nA_10604MeV_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "fa18_out", "Fa18 Out", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_50nA_10604MeV_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "sp19_inb", "Sp19 Inb", 10.200,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_50nA_10200MeV_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root",
    ),
    PeriodConfig(
        "sp18_inb", "Sp18 Inb", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_50nA_10594MeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_50nA_10594MeV_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp18_out", "Sp18 Out", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_45nA_10594MeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_45nA_10594MeV_epgamma.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root",
    ),
)


# Branch aliases are intentionally broad.  The script prints the chosen mapping
# and exits with a precise error when a required physical quantity is absent.
ALIASES: Mapping[str, Tuple[str, ...]] = {
    "runnum": ("runnum", "run", "run_number"),
    "eventnum": ("eventnum", "event", "event_number"),
    "e_p": ("e_p", "p_e", "electron_p"),
    "e_theta": ("e_theta", "theta_e", "electron_theta"),
    "e_phi": ("e_phi", "phi_e", "electron_phi"),
    "p_p": ("p1_p", "p_p", "proton_p"),
    "p_theta": ("p1_theta", "p_theta", "proton_theta"),
    "p_phi": ("p1_phi", "p_phi", "proton_phi"),
    "g1_E": ("p2_p", "g1_p", "g1_E", "gamma1_p", "photon1_p"),
    "g1_theta": ("p2_theta", "g1_theta", "gamma1_theta", "photon1_theta"),
    "g1_phi": ("p2_phi", "g1_phi", "gamma1_phi", "photon1_phi"),
    "g1_detector": ("detector2", "g1_detector", "gamma1_detector"),
    "g2_E": ("p3_p", "g2_p", "g2_E", "gamma2_p", "photon2_p"),
    "g2_theta": ("p3_theta", "g2_theta", "gamma2_theta", "photon2_theta"),
    "g2_phi": ("p3_phi", "g2_phi", "gamma2_phi", "photon2_phi"),
    "g2_detector": ("detector3", "g2_detector", "gamma2_detector"),
    "pi0_p": ("p2_p", "pi0_p"),
    "pi0_theta": ("p2_theta", "pi0_theta"),
    "pi0_phi": ("p2_phi", "pi0_phi"),
    "gamma1_phi_native": ("gamma_phi1",),
    "gamma2_phi_native": ("gamma_phi2",),
    "open_angle_egamma1": ("open_angle_egamma1",),
    "open_angle_egamma2": ("open_angle_egamma2",),
    "gamma1_detector_native": ("detector_gamma1",),
    "gamma2_detector_native": ("detector_gamma2",),
    "Mh_gammagamma": ("Mh_gammagamma",),
    "proton_detector": ("detector1", "p_detector", "proton_detector"),
    "Mx2_1": ("Mx2_1", "Mx2_ep", "Mx2_x1", "Mx2_proton", "Mx2_p"),
    "Mx2": ("Mx2", "Mx2_epg", "Mx2_epgamma", "Mx2_eppi0", "Mx2_epi0"),
    "pTmiss": ("pTmiss", "ptmiss", "pT_miss"),
    "Delta_phi": ("Delta_phi", "delta_phi", "dphi"),
    "theta_gamma_gamma": ("theta_gamma_gamma", "theta_pi0_pi0"),
    "fiducial_status": ("fiducial_status",),
}


@dataclass
class TrialArrays:
    E: np.ndarray
    theta_deg: np.ndarray
    phi_rad: np.ndarray
    detector: np.ndarray
    sector: np.ndarray
    mx2_ep: np.ndarray
    delta_phi: np.ndarray
    theta_gamma_gamma: np.ndarray
    pTmiss: np.ndarray
    weight: np.ndarray

    @staticmethod
    def empty() -> "TrialArrays":
        return TrialArrays(
            E=np.empty(0, dtype=np.float32),
            theta_deg=np.empty(0, dtype=np.float32),
            phi_rad=np.empty(0, dtype=np.float32),
            detector=np.empty(0, dtype=np.int8),
            sector=np.empty(0, dtype=np.int8),
            mx2_ep=np.empty(0, dtype=np.float32),
            delta_phi=np.empty(0, dtype=np.float32),
            theta_gamma_gamma=np.empty(0, dtype=np.float32),
            pTmiss=np.empty(0, dtype=np.float32),
            weight=np.empty(0, dtype=np.float32),
        )

    def size(self) -> int:
        return int(self.E.size)


PASS_EVENT_STAGES = (
    "all_events", "basic_quality", "mx2", "pTmiss", "Delta_phi",
    "native_inputs_finite", "longitudinal_energy_solution",
    "transverse_cone_solution", "detector_compatible_solution", "closure_pass",
)

PASS_PROBE_STAGES = (
    "two_trials_from_closure_events", "prediction_finite", "energy_range",
    "photon_like_E_minus_p", "photon_like_m2", "predicted_detector",
    "angle_match", "energy_match", "final_pass",
)

DIAGNOSTIC_HIST_EDGES = {
    "Mx2": np.linspace(-0.30, 0.30, 121),
    "pTmiss": np.linspace(0.0, 1.0, 101),
    "abs_Delta_phi_minus_pi": np.linspace(0.0, 1.2, 121),
    "closure": np.linspace(0.0, 0.20, 101),
    "ambiguity": np.linspace(0.0, 0.20, 101),
    "transverse_mismatch": np.linspace(0.0, 0.50, 101),
    "energy_fraction_g1": np.linspace(0.0, 1.0, 101),
    "n_solutions": np.arange(-0.5, 8.5, 1.0),
    "g1_E": np.linspace(0.0, 10.0, 101),
    "g2_E": np.linspace(0.0, 10.0, 101),
    "g1_theta_deg": np.linspace(0.0, 40.0, 101),
    "g2_theta_deg": np.linspace(0.0, 40.0, 101),
    "pred_E": np.linspace(0.0, 10.0, 101),
    "pred_theta_deg": np.linspace(0.0, 40.0, 101),
    "pred_m2": np.linspace(-0.5, 0.5, 121),
    "pred_E_minus_p": np.linspace(-0.8, 0.8, 121),
    "match_angle_deg": np.linspace(0.0, 15.0, 121),
    "relative_E_residual": np.linspace(0.0, 2.0, 121),
}


def empty_pass_diagnostics(sample: str, path: str) -> Dict[str, object]:
    return {
        "sample": sample,
        "path": path,
        "event_cutflow": {stage: 0 for stage in PASS_EVENT_STAGES},
        "probe_cutflow": {stage: 0 for stage in PASS_PROBE_STAGES},
        "histograms": {
            key: {"edges": edges.tolist(), "counts": np.zeros(len(edges)-1, dtype=np.int64).tolist()}
            for key, edges in DIAGNOSTIC_HIST_EDGES.items()
        },
        "detector_pairs_after_closure": {},
        "final_trials_by_category": {"FT": 0, **{f"FD sector {i}": 0 for i in range(1, 7)}},
    }


def diag_add_count(diag: Dict[str, object], group: str, stage: str, count: int) -> None:
    diag[group][stage] = int(diag[group].get(stage, 0)) + int(count)


def diag_fill(diag: Dict[str, object], key: str, values: np.ndarray, mask: Optional[np.ndarray] = None) -> None:
    values = np.asarray(values, dtype=float)
    if mask is not None:
        values = values[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return
    edges = np.asarray(diag["histograms"][key]["edges"], dtype=float)
    counts, _ = np.histogram(values, bins=edges)
    old = np.asarray(diag["histograms"][key]["counts"], dtype=np.int64)
    diag["histograms"][key]["counts"] = (old + counts).tolist()


@dataclass(frozen=True)
class BinDefinition:
    bin_id: int
    detector: str
    sector: int
    E_low: float
    E_high: float
    theta_low_deg: float
    theta_high_deg: float


@dataclass
class FitSummary:
    success: bool
    n_pi0: float
    n_pi0_err: float
    n_dvcs: float
    n_dvcs_err: float
    deviance: float
    ndf: int
    message: str
    data_counts: np.ndarray
    model_counts: np.ndarray
    pi0_counts: np.ndarray
    dvcs_counts: np.ndarray
    edges: np.ndarray


@dataclass(frozen=True)
class TemplateMorphology:
    pi0_shift: float = 0.0
    pi0_sigma: float = 0.0
    dvcs_shift: float = 0.0
    dvcs_sigma: float = 0.0
    success: bool = True
    deviance: float = math.nan
    ndf: int = 0


@dataclass
class EfficiencyRow:
    period: str
    period_label: str
    bin_id: int
    detector: str
    sector: int
    E_low: float
    E_high: float
    theta_low_deg: float
    theta_high_deg: float
    pass_data: float
    pass_data_err: float
    fail_data: float
    fail_data_err: float
    efficiency_data: float
    efficiency_data_err: float
    pass_mc: float
    pass_mc_err: float
    fail_mc: float
    fail_mc_err: float
    efficiency_mc: float
    efficiency_mc_err: float
    scale_factor: float
    scale_factor_err: float
    template_fit_success: bool
    template_fit_deviance: float
    template_fit_ndf: int


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Derive RGA photon-efficiency data/MC scale factors with pi0 tag-and-probe.",
    )
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--period", action="append", choices=[p.key for p in PERIODS])
    parser.add_argument("--workers", type=int, default=7,
                        help="Maximum period-level worker processes; capped at 7. "
                             "At most one worker is used per selected period.")
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument("--inspect-branches", action="store_true")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Debug-only maximum events read from each tree.")

    parser.add_argument("--ft-theta-min", type=float, default=2.5)
    parser.add_argument("--ft-theta-max", type=float, default=5.0)
    parser.add_argument("--fd-theta-min", type=float, default=5.0)
    parser.add_argument("--fd-theta-max", type=float, default=35.0)
    parser.add_argument("--probe-E-min", type=float, default=0.35)
    parser.add_argument("--probe-E-max", type=float, default=9.5)

    parser.add_argument("--ft-energy-bins", type=int, default=6)
    parser.add_argument("--fd-energy-bins", type=int, default=4)
    parser.add_argument("--fd-theta-bins", type=int, default=3)
    parser.add_argument("--min-probes-per-bin", type=int, default=250)

    parser.add_argument("--mx2-ep-min", type=float, default=-0.12)
    parser.add_argument("--mx2-ep-max", type=float, default=0.20)
    parser.add_argument("--mx2-ep-bins", type=int, default=80)
    parser.add_argument("--pass-mx2-abs-max", type=float, default=0.05,
                        help="Broad |Mx2| cut for direct epgamma-gamma pass trials, if Mx2 exists.")
    parser.add_argument("--pass-pTmiss-max", type=float, default=0.30,
                        help="Broad pTmiss cut for direct epgamma-gamma pass trials, if pTmiss exists.")
    parser.add_argument("--pass-delta-phi-window", type=float, default=0.30,
                        help="Broad native eppi0 proton-pi0 coplanarity requirement |Delta_phi-pi|. Set to a negative value to disable.")
    parser.add_argument("--disable-template-morphing", action="store_true",
                        help="Disable shared detector/sector nuisance morphing in the multi-projection fit.")
    parser.add_argument("--fit-min-counts", type=int, default=100,
                        help="Minimum data and template entries required for a differential multi-projection fit.")
    parser.add_argument("--probe-match-angle-max-deg", type=float, default=2.0,
                        help="Maximum angle between predicted and observed probe in epgamma-gamma events.")
    parser.add_argument("--probe-match-relative-E-max", type=float, default=0.35,
                        help="Maximum |Eobs-Epred|/Epred for a passing probe.")
    parser.add_argument("--pass-photon-closure-max", type=float, default=0.03,
                        help="Maximum normalized four-vector closure residual for reconstructed native eppi0 photons.")
    parser.add_argument("--cone-energy-denominator-min", type=float, default=1.0e-5,
                        help="Minimum |cos(alpha1)-cos(alpha2)| for the analytic photon-energy solution.")
    parser.add_argument("--cone-min-photon-energy", type=float, default=0.05,
                        help="Minimum physical photon energy (GeV) in the native eppi0 cone solution.")
    parser.add_argument("--cone-mass-residual-weight", type=float, default=0.05,
                        help="Relative weight of the diphoton-mass residual when choosing between mirror solutions.")
    parser.add_argument("--diagnostic-max-examples", type=int, default=0,
                        help="Reserved for future event-level diagnostic dumps; zero writes no event records.")
    parser.add_argument("--skip-efficiency-extraction", action="store_true",
                        help="Run the passing-sample audit and stop before epgamma template fits.")
    parser.add_argument("--require-fiducial-status", type=int, default=None,
                        help="Optional exact fiducial_status requirement.")

    parser.add_argument("--epg-data", default=None,
                        help="Override epgamma data path; valid only with one --period.")
    parser.add_argument("--dvcs-mc", default=None)
    parser.add_argument("--pi0-as-epg-mc", default=None)
    parser.add_argument("--epgg-data", default=None)
    parser.add_argument("--epgg-mc", default=None)
    return parser.parse_args()


def require_file(path: str) -> None:
    if not Path(path).is_file():
        raise FileNotFoundError(f"Input ROOT file does not exist: {path}")
    # endif


def tree_keys(path: str) -> List[str]:
    require_file(path)
    with uproot.open(path) as root_file:
        if TREE_NAME not in root_file:
            raise KeyError(f"Tree '{TREE_NAME}' not found in {path}")
        # endif
        return [str(key) for key in root_file[TREE_NAME].keys()]
    # endwith


def resolve_branches(path: str, logical_names: Sequence[str], optional: Sequence[str] = ()) -> Dict[str, Optional[str]]:
    available = set(tree_keys(path))
    optional_set = set(optional)
    resolved: Dict[str, Optional[str]] = {}
    missing: List[str] = []
    for logical in logical_names:
        aliases = ALIASES.get(logical, (logical,))
        branch = next((name for name in aliases if name in available), None)
        resolved[logical] = branch
        if branch is None and logical not in optional_set:
            missing.append(f"{logical}: tried {aliases}")
        # endif
    # endfor
    if missing:
        raise KeyError(
            f"Required branches are missing from {path}:\n  " + "\n  ".join(missing)
            + "\nAvailable branches include:\n  " + ", ".join(sorted(available)[:200])
        )
    # endif
    return resolved


def spherical_to_cartesian(p: np.ndarray, theta: np.ndarray, phi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    st = np.sin(theta)
    return p * st * np.cos(phi), p * st * np.sin(phi), p * np.cos(theta)


def four_vector_from_p_theta_phi(
    p: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    mass: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical_to_cartesian(p, theta, phi)
    energy = np.sqrt(np.maximum(p * p + mass * mass, 0.0))
    return energy, px, py, pz


def photon_four_vector(E: np.ndarray, theta: np.ndarray, phi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical_to_cartesian(E, theta, phi)
    return E, px, py, pz


def predicted_probe(
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
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ee, epx, epy, epz = four_vector_from_p_theta_phi(e_p, e_theta, e_phi, ELECTRON_MASS_GEV)
    pe, ppx, ppy, ppz = four_vector_from_p_theta_phi(p_p, p_theta, p_phi, PROTON_MASS_GEV)
    ge, gpx, gpy, gpz = photon_four_vector(tag_E, tag_theta, tag_phi)

    pred_E = beam_energy + PROTON_MASS_GEV - ee - pe - ge
    pred_px = -epx - ppx - gpx
    pred_py = -epy - ppy - gpy
    pred_pz = beam_energy - epz - ppz - gpz
    pred_p = np.sqrt(np.maximum(pred_px * pred_px + pred_py * pred_py + pred_pz * pred_pz, 0.0))
    pred_theta = np.arctan2(np.sqrt(pred_px * pred_px + pred_py * pred_py), pred_pz)
    pred_phi = np.arctan2(pred_py, pred_px)
    pred_m2 = pred_E * pred_E - pred_p * pred_p
    return pred_E, pred_theta, pred_phi, pred_m2, pred_p


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


def classify_predicted_detector(theta_deg: np.ndarray, phi_rad: np.ndarray, args: argparse.Namespace) -> Tuple[np.ndarray, np.ndarray]:
    detector = np.full(theta_deg.shape, -1, dtype=np.int16)
    detector[(theta_deg >= args.ft_theta_min) & (theta_deg < args.ft_theta_max)] = 0
    detector[(theta_deg >= args.fd_theta_min) & (theta_deg < args.fd_theta_max)] = 1
    sector = fd_sector_from_phi(phi_rad)
    sector[detector != 1] = 0
    return detector, sector


def concatenate_trials(chunks: Sequence[TrialArrays]) -> TrialArrays:
    if not chunks:
        return TrialArrays.empty()
    # endif
    return TrialArrays(
        E=np.concatenate([chunk.E for chunk in chunks]).astype(np.float32, copy=False),
        theta_deg=np.concatenate([chunk.theta_deg for chunk in chunks]).astype(np.float32, copy=False),
        phi_rad=np.concatenate([chunk.phi_rad for chunk in chunks]).astype(np.float32, copy=False),
        detector=np.concatenate([chunk.detector for chunk in chunks]).astype(np.int8, copy=False),
        sector=np.concatenate([chunk.sector for chunk in chunks]).astype(np.int8, copy=False),
        mx2_ep=np.concatenate([chunk.mx2_ep for chunk in chunks]).astype(np.float32, copy=False),
        delta_phi=np.concatenate([chunk.delta_phi for chunk in chunks]).astype(np.float32, copy=False),
        theta_gamma_gamma=np.concatenate([chunk.theta_gamma_gamma for chunk in chunks]).astype(np.float32, copy=False),
        pTmiss=np.concatenate([chunk.pTmiss for chunk in chunks]).astype(np.float32, copy=False),
        weight=np.concatenate([chunk.weight for chunk in chunks]).astype(np.float32, copy=False),
    )


def finite_array(arrays: Mapping[str, np.ndarray], branch: Optional[str], default: float = math.nan) -> np.ndarray:
    if branch is None:
        first = next(iter(arrays.values()))
        return np.full(len(first), default, dtype=float)
    # endif
    return np.asarray(arrays[branch], dtype=float)


def basic_quality_mask(arrays: Mapping[str, np.ndarray], resolved: Mapping[str, Optional[str]], args: argparse.Namespace) -> np.ndarray:
    e_p = finite_array(arrays, resolved["e_p"])
    p_p = finite_array(arrays, resolved["p_p"])
    mask = np.isfinite(e_p) & np.isfinite(p_p) & (e_p > 0.0) & (p_p > 0.0)
    if args.require_fiducial_status is not None and resolved.get("fiducial_status") is not None:
        status = np.asarray(arrays[resolved["fiducial_status"]], dtype=int)
        mask &= status == args.require_fiducial_status
    # endif
    return mask


def pass_sample_mask(arrays: Mapping[str, np.ndarray], resolved: Mapping[str, Optional[str]], args: argparse.Namespace) -> np.ndarray:
    mask = basic_quality_mask(arrays, resolved, args)
    if resolved.get("Mx2") is not None:
        mx2 = finite_array(arrays, resolved["Mx2"])
        mask &= np.isfinite(mx2) & (np.abs(mx2) < args.pass_mx2_abs_max)
    # endif
    if resolved.get("pTmiss") is not None:
        pt = finite_array(arrays, resolved["pTmiss"])
        mask &= np.isfinite(pt) & (pt < args.pass_pTmiss_max)
    # endif
    if args.pass_delta_phi_window is not None and args.pass_delta_phi_window >= 0.0 and resolved.get("Delta_phi") is not None:
        dphi = finite_array(arrays, resolved["Delta_phi"])
        mask &= np.isfinite(dphi) & (np.abs(dphi - math.pi) < args.pass_delta_phi_window)
    # endif
    return mask



def _unit_vector(theta: float, phi: float) -> np.ndarray:
    st = math.sin(theta)
    return np.asarray([st * math.cos(phi), st * math.sin(phi), math.cos(theta)], dtype=float)


def _orthonormal_basis_about(axis: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors transverse to a unit axis."""
    axis = np.asarray(axis, dtype=float)
    if abs(axis[2]) < 0.9:
        reference = np.asarray([0.0, 0.0, 1.0], dtype=float)
    else:
        reference = np.asarray([1.0, 0.0, 0.0], dtype=float)
    u = np.cross(reference, axis)
    norm_u = float(np.linalg.norm(u))
    if norm_u <= 1.0e-14:
        reference = np.asarray([0.0, 1.0, 0.0], dtype=float)
        u = np.cross(reference, axis)
        norm_u = float(np.linalg.norm(u))
    u /= norm_u
    v = np.cross(axis, u)
    v /= float(np.linalg.norm(v))
    return u, v


def _detector_compatible(theta_rad: float, detector: int, args: argparse.Namespace) -> bool:
    theta_deg = math.degrees(theta_rad)
    if detector == 0:
        return args.ft_theta_min <= theta_deg < args.ft_theta_max
    if detector == 1:
        return args.fd_theta_min <= theta_deg < args.fd_theta_max
    return False


def reconstruct_native_eppi0_photons(
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    pi0_p: np.ndarray,
    pi0_theta: np.ndarray,
    pi0_phi: np.ndarray,
    opening1: np.ndarray,
    opening2: np.ndarray,
    pi0_mass: np.ndarray,
    detector1: np.ndarray,
    detector2: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, ...]:
    """Recover the two photons using cone constraints around the electron.

    For photon i, the measured electron-photon opening angle alpha_i places its
    direction on a cone around the measured electron direction.  The measured
    pi0 four-vector then fixes the photon energies through the total-energy and
    electron-axis momentum equations.  The transverse pi0 momentum determines
    the two azimuths around the electron axis up to a mirror ambiguity.

    When finite detector resolution makes the transverse triangle slightly
    inconsistent, the nearest triangle boundary is used and the resulting
    normalized four-vector mismatch is retained as the closure diagnostic.
    """
    n_events = len(e_theta)
    g1_E = np.full(n_events, np.nan, dtype=float)
    g1_theta = np.full(n_events, np.nan, dtype=float)
    g1_phi = np.full(n_events, np.nan, dtype=float)
    g2_E = np.full(n_events, np.nan, dtype=float)
    g2_theta = np.full(n_events, np.nan, dtype=float)
    g2_phi = np.full(n_events, np.nan, dtype=float)
    closure = np.full(n_events, np.nan, dtype=float)
    ambiguity = np.full(n_events, np.nan, dtype=float)
    transverse_mismatch = np.full(n_events, np.nan, dtype=float)
    energy_fraction_g1 = np.full(n_events, np.nan, dtype=float)
    n_solutions = np.zeros(n_events, dtype=np.int16)
    input_finite = np.zeros(n_events, dtype=bool)
    longitudinal_solution = np.zeros(n_events, dtype=bool)
    transverse_solution = np.zeros(n_events, dtype=bool)
    detector_solution = np.zeros(n_events, dtype=bool)
    closure_pass = np.zeros(n_events, dtype=bool)

    for i in range(n_events):
        values = (
            e_theta[i], e_phi[i], pi0_p[i], pi0_theta[i], pi0_phi[i],
            opening1[i], opening2[i], pi0_mass[i], detector1[i], detector2[i],
        )
        if not all(math.isfinite(float(v)) for v in values):
            continue
        if pi0_p[i] <= 0.0 or pi0_mass[i] <= 0.0:
            continue
        if not (0.0 < opening1[i] < math.pi and 0.0 < opening2[i] < math.pi):
            continue
        input_finite[i] = True

        ehat = _unit_vector(float(e_theta[i]), float(e_phi[i]))
        npi = _unit_vector(float(pi0_theta[i]), float(pi0_phi[i]))
        pvec = float(pi0_p[i]) * npi
        epi = math.sqrt(float(pi0_p[i]) ** 2 + float(pi0_mass[i]) ** 2)
        p_parallel = float(np.dot(pvec, ehat))

        c1 = math.cos(float(opening1[i]))
        c2 = math.cos(float(opening2[i]))
        denom = c1 - c2
        if abs(denom) <= args.cone_energy_denominator_min:
            continue
        E1 = (p_parallel - epi * c2) / denom
        E2 = epi - E1
        if E1 <= args.cone_min_photon_energy or E2 <= args.cone_min_photon_energy:
            continue
        longitudinal_solution[i] = True
        energy_fraction_g1[i] = E1 / epi

        u, v = _orthonormal_basis_about(ehat)
        p_perp_vec = pvec - p_parallel * ehat
        px = float(np.dot(p_perp_vec, u))
        py = float(np.dot(p_perp_vec, v))
        R = math.hypot(px, py)
        psi = math.atan2(py, px) if R > 1.0e-14 else 0.0
        A = E1 * math.sin(float(opening1[i]))
        B = E2 * math.sin(float(opening2[i]))
        if A <= 1.0e-14 or B <= 1.0e-14:
            continue

        lower = abs(A - B)
        upper = A + B
        target_R = min(max(R, lower), upper)
        transverse_mismatch[i] = abs(R - target_R) / max(epi, 1.0e-12)

        if target_R <= 1.0e-14:
            # Degenerate transverse configuration: choose an arbitrary axis;
            # detector matching and closure will decide whether it is usable.
            delta1 = 0.0
            delta2 = math.pi if A >= B else 0.0
        else:
            cos_delta1 = (target_R * target_R + A * A - B * B) / (2.0 * target_R * A)
            cos_delta2 = (target_R * target_R + B * B - A * A) / (2.0 * target_R * B)
            delta1 = math.acos(float(np.clip(cos_delta1, -1.0, 1.0)))
            delta2 = math.acos(float(np.clip(cos_delta2, -1.0, 1.0)))
        transverse_solution[i] = True

        candidates: List[Tuple[float, float, float, float, float, float, float]] = []
        for sign in (+1.0, -1.0):
            beta1 = psi + sign * delta1
            beta2 = psi - sign * delta2
            n1 = c1 * ehat + math.sin(float(opening1[i])) * (math.cos(beta1) * u + math.sin(beta1) * v)
            n2 = c2 * ehat + math.sin(float(opening2[i])) * (math.cos(beta2) * u + math.sin(beta2) * v)
            n1 /= float(np.linalg.norm(n1))
            n2 /= float(np.linalg.norm(n2))
            th1 = math.acos(float(np.clip(n1[2], -1.0, 1.0)))
            ph1 = math.atan2(float(n1[1]), float(n1[0]))
            th2 = math.acos(float(np.clip(n2[2], -1.0, 1.0)))
            ph2 = math.atan2(float(n2[1]), float(n2[0]))
            if not (_detector_compatible(th1, int(detector1[i]), args) and
                    _detector_compatible(th2, int(detector2[i]), args)):
                continue
            four_residual = np.asarray([
                E1 + E2 - epi,
                *(E1 * n1 + E2 * n2 - pvec),
            ], dtype=float)
            residual = float(np.linalg.norm(four_residual) / max(epi, 1.0e-12))
            mass2 = 2.0 * E1 * E2 * (1.0 - float(np.dot(n1, n2)))
            mass_residual = abs(mass2 - float(pi0_mass[i]) ** 2) / max(float(pi0_mass[i]) ** 2, 1.0e-12)
            score = math.hypot(residual, args.cone_mass_residual_weight * mass_residual)
            candidates.append((score, residual, mass_residual, th1, ph1, th2, ph2))

        n_solutions[i] = len(candidates)
        if not candidates:
            continue
        detector_solution[i] = True
        candidates.sort(key=lambda item: item[0])
        best = candidates[0]
        closure[i] = best[1]
        ambiguity[i] = candidates[1][0] - best[0] if len(candidates) > 1 else math.inf
        g1_E[i], g2_E[i] = E1, E2
        g1_theta[i], g1_phi[i] = best[3], best[4]
        g2_theta[i], g2_phi[i] = best[5], best[6]
        closure_pass[i] = best[1] <= args.pass_photon_closure_max

    return (
        g1_E, g1_theta, g1_phi, g2_E, g2_theta, g2_phi,
        closure, ambiguity, transverse_mismatch, energy_fraction_g1,
        n_solutions, input_finite, longitudinal_solution, transverse_solution,
        detector_solution, closure_pass,
    )

def read_pass_trials(path: str, beam_energy: float, args: argparse.Namespace,
                     sample_name: str) -> Tuple[TrialArrays, Dict[str, object]]:
    """Read native eppi0 trees, construct directed probes, and audit every loss."""
    logical = [
        "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi",
        "pi0_p", "pi0_theta", "pi0_phi",
        "open_angle_egamma1", "open_angle_egamma2", "gamma1_detector_native",
        "gamma2_detector_native", "Mh_gammagamma", "Mx2_1", "Mx2", "pTmiss",
        "Delta_phi", "fiducial_status",
    ]
    optional = ["Mx2_1", "Mx2", "pTmiss", "Delta_phi", "fiducial_status"]
    resolved = resolve_branches(path, logical, optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    chunks: List[TrialArrays] = []
    seen = 0
    diag = empty_pass_diagnostics(sample_name, path)

    log(f"Reading native eppi0 PASS trials from {path}")
    for arrays in uproot.iterate(f"{path}:{TREE_NAME}", expressions=expressions,
                                 step_size=args.step_size, library="np"):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        if args.max_events is not None and seen + n > args.max_events:
            keep = args.max_events - seen
            arrays = {key: value[:keep] for key, value in arrays.items()}
            n = keep
        seen += n
        diag_add_count(diag, "event_cutflow", "all_events", n)

        quality = basic_quality_mask(arrays, resolved, args)
        mx2 = finite_array(arrays, resolved.get("Mx2"))
        ptmiss_all = finite_array(arrays, resolved.get("pTmiss"))
        delta_phi_all = finite_array(arrays, resolved.get("Delta_phi"))
        mx2_stage = quality.copy()
        if resolved.get("Mx2") is not None:
            mx2_stage &= np.isfinite(mx2) & (np.abs(mx2) < args.pass_mx2_abs_max)
        pt_stage = mx2_stage.copy()
        if resolved.get("pTmiss") is not None:
            pt_stage &= np.isfinite(ptmiss_all) & (ptmiss_all < args.pass_pTmiss_max)
        dphi_stage = pt_stage.copy()
        if args.pass_delta_phi_window is not None and args.pass_delta_phi_window >= 0.0 and resolved.get("Delta_phi") is not None:
            dphi_stage &= np.isfinite(delta_phi_all) & (np.abs(delta_phi_all - math.pi) < args.pass_delta_phi_window)
        base = dphi_stage

        diag_add_count(diag, "event_cutflow", "basic_quality", np.count_nonzero(quality))
        diag_add_count(diag, "event_cutflow", "mx2", np.count_nonzero(mx2_stage))
        diag_add_count(diag, "event_cutflow", "pTmiss", np.count_nonzero(pt_stage))
        diag_add_count(diag, "event_cutflow", "Delta_phi", np.count_nonzero(dphi_stage))
        diag_fill(diag, "Mx2", mx2, quality)
        diag_fill(diag, "pTmiss", ptmiss_all, mx2_stage)
        diag_fill(diag, "abs_Delta_phi_minus_pi", np.abs(delta_phi_all-math.pi), pt_stage)

        e_p = finite_array(arrays, resolved["e_p"]); e_theta = finite_array(arrays, resolved["e_theta"]); e_phi = finite_array(arrays, resolved["e_phi"])
        p_p = finite_array(arrays, resolved["p_p"]); p_theta = finite_array(arrays, resolved["p_theta"]); p_phi = finite_array(arrays, resolved["p_phi"])
        pi0_p = finite_array(arrays, resolved["pi0_p"]); pi0_theta = finite_array(arrays, resolved["pi0_theta"]); pi0_phi = finite_array(arrays, resolved["pi0_phi"])
        opening1 = finite_array(arrays, resolved["open_angle_egamma1"]); opening2 = finite_array(arrays, resolved["open_angle_egamma2"])
        detector1_obs = finite_array(arrays, resolved["gamma1_detector_native"], default=-1.0).astype(int)
        detector2_obs = finite_array(arrays, resolved["gamma2_detector_native"], default=-1.0).astype(int)
        pi0_mass = finite_array(arrays, resolved["Mh_gammagamma"])
        mx2_ep = finite_array(arrays, resolved.get("Mx2_1"))
        theta_gg_all = np.full(n, np.nan, dtype=float)

        rec = reconstruct_native_eppi0_photons(
            e_theta, e_phi, pi0_p, pi0_theta, pi0_phi, opening1, opening2,
            pi0_mass, detector1_obs, detector2_obs, args)
        (g1_E, g1_theta, g1_phi, g2_E, g2_theta, g2_phi, closure, ambiguity,
         transverse_mismatch, energy_fraction_g1, n_solutions, input_finite,
         longitudinal_solution, transverse_solution, detector_solution,
         closure_pass) = rec

        for stage, mask in (("native_inputs_finite", input_finite),
                            ("longitudinal_energy_solution", longitudinal_solution),
                            ("transverse_cone_solution", transverse_solution),
                            ("detector_compatible_solution", detector_solution),
                            ("closure_pass", closure_pass)):
            diag_add_count(diag, "event_cutflow", stage, np.count_nonzero(base & mask))
        diag_fill(diag, "closure", closure, base & detector_solution)
        diag_fill(diag, "ambiguity", ambiguity, base & detector_solution & np.isfinite(ambiguity))
        diag_fill(diag, "transverse_mismatch", transverse_mismatch, base & longitudinal_solution)
        diag_fill(diag, "energy_fraction_g1", energy_fraction_g1, base & longitudinal_solution)
        diag_fill(diag, "n_solutions", n_solutions, base & transverse_solution)
        diag_fill(diag, "g1_E", g1_E, base & closure_pass)
        diag_fill(diag, "g2_E", g2_E, base & closure_pass)
        diag_fill(diag, "g1_theta_deg", np.degrees(g1_theta), base & closure_pass)
        diag_fill(diag, "g2_theta_deg", np.degrees(g2_theta), base & closure_pass)

        pair_keys = np.char.add(np.char.add(detector1_obs.astype(str), "-"), detector2_obs.astype(str))
        for key in np.unique(pair_keys[base & closure_pass]):
            diag["detector_pairs_after_closure"][str(key)] = int(diag["detector_pairs_after_closure"].get(str(key), 0)) + int(np.count_nonzero((pair_keys == key) & base & closure_pass))

        pred_sets = (
            (*predicted_probe(beam_energy, e_p, e_theta, e_phi, p_p, p_theta, p_phi, g1_E, g1_theta, g1_phi), g2_E, g2_theta, g2_phi),
            (*predicted_probe(beam_energy, e_p, e_theta, e_phi, p_p, p_theta, p_phi, g2_E, g2_theta, g2_phi), g1_E, g1_theta, g1_phi),
        )
        for pred_E, pred_th, pred_ph, pred_m2, pred_p, obs_E, obs_th, obs_ph in pred_sets:
            seed = base & closure_pass
            diag_add_count(diag, "probe_cutflow", "two_trials_from_closure_events", np.count_nonzero(seed))
            theta_deg = np.degrees(pred_th)
            detector, sector = classify_predicted_detector(theta_deg, pred_ph, args)
            cos_opening = np.sin(pred_th)*np.sin(obs_th)*np.cos(pred_ph-obs_ph)+np.cos(pred_th)*np.cos(obs_th)
            match_angle_deg = np.degrees(np.arccos(np.clip(cos_opening, -1.0, 1.0)))
            relative_E_residual = np.abs(obs_E-pred_E)/np.maximum(pred_E, 1.0e-9)
            finite = seed & np.isfinite(pred_E) & np.isfinite(pred_p) & np.isfinite(theta_deg) & np.isfinite(pred_ph) & np.isfinite(obs_E) & np.isfinite(obs_th) & np.isfinite(obs_ph) & np.isfinite(match_angle_deg) & np.isfinite(relative_E_residual)
            energy = finite & (pred_E >= args.probe_E_min) & (pred_E < args.probe_E_max)
            ep = energy & (pred_p > 0.0) & (np.abs(pred_E-pred_p) < 0.30)
            m2 = ep & (np.abs(pred_m2) < 0.20)
            accepted = m2 & (detector >= 0)
            angle = accepted & (match_angle_deg < args.probe_match_angle_max_deg)
            eres = angle & (relative_E_residual < args.probe_match_relative_E_max)
            for stage, mask in (("prediction_finite", finite), ("energy_range", energy), ("photon_like_E_minus_p", ep),
                                ("photon_like_m2", m2), ("predicted_detector", accepted),
                                ("angle_match", angle), ("energy_match", eres), ("final_pass", eres)):
                diag_add_count(diag, "probe_cutflow", stage, np.count_nonzero(mask))
            diag_fill(diag, "pred_E", pred_E, finite); diag_fill(diag, "pred_theta_deg", theta_deg, finite)
            diag_fill(diag, "pred_m2", pred_m2, ep); diag_fill(diag, "pred_E_minus_p", pred_E-pred_p, energy)
            diag_fill(diag, "match_angle_deg", match_angle_deg, accepted)
            diag_fill(diag, "relative_E_residual", relative_E_residual, angle)
            for cat, cmask in [("FT", eres & (detector == 0))] + [(f"FD sector {j}", eres & (detector == 1) & (sector == j)) for j in range(1,7)]:
                diag["final_trials_by_category"][cat] += int(np.count_nonzero(cmask))
            chunks.append(TrialArrays(
                E=pred_E[eres].astype(np.float32, copy=False), theta_deg=theta_deg[eres].astype(np.float32, copy=False),
                phi_rad=pred_ph[eres].astype(np.float32, copy=False), detector=detector[eres].astype(np.int8, copy=False),
                sector=sector[eres].astype(np.int8, copy=False), mx2_ep=mx2_ep[eres].astype(np.float32, copy=False),
                delta_phi=delta_phi_all[eres].astype(np.float32, copy=False), theta_gamma_gamma=theta_gg_all[eres].astype(np.float32, copy=False),
                pTmiss=ptmiss_all[eres].astype(np.float32, copy=False), weight=np.ones(np.count_nonzero(eres), dtype=np.float32)))
    result = concatenate_trials(chunks)
    log(f"Native eppi0 audit for {Path(path).name}: events={seen:,}, closure events={diag['event_cutflow']['closure_pass']:,}, final directed trials={result.size():,}")
    return result, diag

def read_fail_trials(path: str, beam_energy: float, args: argparse.Namespace) -> TrialArrays:
    logical = [
        "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi",
        "g1_E", "g1_theta", "g1_phi", "Mx2_1", "Delta_phi",
        "theta_gamma_gamma", "pTmiss", "fiducial_status",
    ]
    optional = ["fiducial_status"]
    resolved = resolve_branches(path, logical, optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    chunks: List[TrialArrays] = []
    seen = 0

    log(f"Reading FAIL candidates from {path}")
    for arrays in uproot.iterate(f"{path}:{TREE_NAME}", expressions=expressions,
                                 step_size=args.step_size, library="np"):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        # endif
        if args.max_events is not None and seen + n > args.max_events:
            keep = args.max_events - seen
            arrays = {key: value[:keep] for key, value in arrays.items()}
            n = keep
        # endif
        seen += n

        base = basic_quality_mask(arrays, resolved, args)
        e_p = finite_array(arrays, resolved["e_p"])
        e_theta = finite_array(arrays, resolved["e_theta"])
        e_phi = finite_array(arrays, resolved["e_phi"])
        p_p = finite_array(arrays, resolved["p_p"])
        p_theta = finite_array(arrays, resolved["p_theta"])
        p_phi = finite_array(arrays, resolved["p_phi"])
        tag_E = finite_array(arrays, resolved["g1_E"])
        tag_theta = finite_array(arrays, resolved["g1_theta"])
        tag_phi = finite_array(arrays, resolved["g1_phi"])
        mx2_ep = finite_array(arrays, resolved["Mx2_1"])

        pred_E, pred_th, pred_ph, pred_m2, pred_p = predicted_probe(
            beam_energy, e_p, e_theta, e_phi, p_p, p_theta, p_phi,
            tag_E, tag_theta, tag_phi,
        )
        theta_deg = np.degrees(pred_th)
        detector, sector = classify_predicted_detector(theta_deg, pred_ph, args)
        good = (
            base
            & np.isfinite(pred_E) & np.isfinite(pred_p)
            & np.isfinite(theta_deg) & np.isfinite(pred_ph)
            & np.isfinite(mx2_ep)
            & (pred_E >= args.probe_E_min) & (pred_E < args.probe_E_max)
            & (pred_p > 0.0)
            & (np.abs(pred_E - pred_p) < 0.30)
            & (np.abs(pred_m2) < 0.20)
            & (detector >= 0)
        )
        chunks.append(TrialArrays(
            E=pred_E[good].astype(np.float32, copy=False),
            theta_deg=theta_deg[good].astype(np.float32, copy=False),
            phi_rad=pred_ph[good].astype(np.float32, copy=False),
            detector=detector[good].astype(np.int8, copy=False),
            sector=sector[good].astype(np.int8, copy=False),
            mx2_ep=mx2_ep[good].astype(np.float32, copy=False),
            delta_phi=finite_array(arrays, resolved.get("Delta_phi"))[good].astype(np.float32, copy=False),
            theta_gamma_gamma=finite_array(arrays, resolved.get("theta_gamma_gamma"))[good].astype(np.float32, copy=False),
            pTmiss=finite_array(arrays, resolved.get("pTmiss"))[good].astype(np.float32, copy=False),
            weight=np.ones(np.count_nonzero(good), dtype=np.float32),
        ))
    # endfor
    return concatenate_trials(chunks)


def strict_quantile_edges(values: np.ndarray, n_bins: int, low: float, high: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values >= low) & (values < high)]
    if values.size == 0 or n_bins <= 1:
        return np.asarray([low, high], dtype=float)
    # endif
    quantiles = np.linspace(0.0, 1.0, n_bins + 1)
    edges = np.quantile(values, quantiles)
    edges[0] = low
    edges[-1] = high
    edges = np.unique(edges)
    if edges.size < 2:
        return np.asarray([low, high], dtype=float)
    # endif
    return edges


def build_bins(pass_mc: TrialArrays, fail_mc: TrialArrays, args: argparse.Namespace) -> List[BinDefinition]:
    """Build the deliberately coarse first-pass detector binning.

    The study is integrated over the full configured photon-energy and
    polar-angle acceptance.  The only retained categories are FT and the six
    FD sectors.  The MC arguments are kept in the signature so later versions
    can restore adaptive differential binning without changing callers.
    """
    del pass_mc, fail_mc
    bins: List[BinDefinition] = [
        BinDefinition(
            0, "FT", 0,
            float(args.probe_E_min), float(args.probe_E_max),
            float(args.ft_theta_min), float(args.ft_theta_max),
        )
    ]
    for sector in range(1, 7):
        bins.append(
            BinDefinition(
                sector, "FD", sector,
                float(args.probe_E_min), float(args.probe_E_max),
                float(args.fd_theta_min), float(args.fd_theta_max),
            )
        )
    return bins


def bin_mask(trials: TrialArrays, definition: BinDefinition) -> np.ndarray:
    if definition.detector == "FT":
        detector_mask = trials.detector == 0
    else:
        detector_mask = (trials.detector == 1) & (trials.sector == definition.sector)
    # endif
    return (
        detector_mask
        & (trials.E >= definition.E_low) & (trials.E < definition.E_high)
        & (trials.theta_deg >= definition.theta_low_deg)
        & (trials.theta_deg < definition.theta_high_deg)
    )


def assign_bin_ids(trials: TrialArrays, definitions: Sequence[BinDefinition]) -> np.ndarray:
    """Assign every trial to at most one adaptive probe bin in one vectorized pass.

    The bin definitions are disjoint by construction.  This replaces repeated
    full-array scans for every bin with one scan per detector/sector category.
    """
    assigned = np.full(trials.size(), -1, dtype=np.int16)
    if trials.size() == 0:
        return assigned
    # endif

    categories: Dict[Tuple[str, int], List[BinDefinition]] = {}
    for definition in definitions:
        categories.setdefault((definition.detector, definition.sector), []).append(definition)
    # endfor

    for (detector_name, sector), category_bins in categories.items():
        if detector_name == "FT":
            category_mask = trials.detector == 0
        else:
            category_mask = (trials.detector == 1) & (trials.sector == sector)
        # endif
        indices = np.flatnonzero(category_mask)
        if indices.size == 0:
            continue
        # endif
        energy = trials.E[indices]
        theta = trials.theta_deg[indices]
        for definition in category_bins:
            local = (
                (energy >= definition.E_low) & (energy < definition.E_high)
                & (theta >= definition.theta_low_deg)
                & (theta < definition.theta_high_deg)
            )
            assigned[indices[local]] = definition.bin_id
        # endfor
    # endfor
    return assigned


def weighted_counts_by_bin(trials: TrialArrays, bin_ids: np.ndarray, n_bins: int) -> np.ndarray:
    valid = (bin_ids >= 0) & (bin_ids < n_bins)
    return np.bincount(
        bin_ids[valid].astype(np.int64, copy=False),
        weights=trials.weight[valid],
        minlength=n_bins,
    ).astype(float, copy=False)


def bulk_mx2_histograms(
    trials: TrialArrays,
    bin_ids: np.ndarray,
    n_probe_bins: int,
    mx2_edges: np.ndarray,
) -> np.ndarray:
    """Fill all probe-bin Mx2 histograms in one np.bincount operation."""
    n_mx2_bins = len(mx2_edges) - 1
    mx2_index = np.searchsorted(mx2_edges, trials.mx2_ep, side="right") - 1
    valid = (
        (bin_ids >= 0) & (bin_ids < n_probe_bins)
        & np.isfinite(trials.mx2_ep)
        & (mx2_index >= 0) & (mx2_index < n_mx2_bins)
    )
    flat_index = (
        bin_ids[valid].astype(np.int64, copy=False) * n_mx2_bins
        + mx2_index[valid].astype(np.int64, copy=False)
    )
    counts = np.bincount(
        flat_index,
        weights=trials.weight[valid],
        minlength=n_probe_bins * n_mx2_bins,
    )
    return counts.reshape(n_probe_bins, n_mx2_bins).astype(float, copy=False)


def poisson_deviance(data: np.ndarray, model: np.ndarray) -> float:
    data = np.asarray(data, dtype=float)
    model = np.clip(np.asarray(model, dtype=float), 1e-12, None)
    terms = model - data
    positive = data > 0.0
    terms[positive] += data[positive] * np.log(data[positive] / model[positive])
    return float(2.0 * np.sum(terms))


def morph_template_histogram(counts: np.ndarray, edges: np.ndarray, shift: float, sigma: float) -> np.ndarray:
    """Shift and Gaussian-smear a histogram while preserving its normalization."""
    counts = np.asarray(counts, dtype=float)
    total = float(counts.sum())
    if total <= 0.0:
        return np.zeros_like(counts, dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    shifted = np.interp(centers - shift, centers, counts, left=0.0, right=0.0)
    mean_width = float(np.mean(np.diff(edges)))
    if sigma > 0.0 and mean_width > 0.0:
        shifted = gaussian_filter1d(shifted, sigma / mean_width, mode="constant", cval=0.0)
    shifted = np.clip(shifted, 0.0, None)
    norm = float(shifted.sum())
    return shifted * (total / norm) if norm > 0.0 else np.zeros_like(counts, dtype=float)


def fit_shared_template_morphology(
    data_counts: np.ndarray,
    dvcs_counts: np.ndarray,
    pi0_counts: np.ndarray,
    edges: np.ndarray,
) -> TemplateMorphology:
    """Fit one shared shift/smearing model for a period and detector category."""
    data = np.asarray(data_counts, dtype=float)
    dvcs = np.asarray(dvcs_counts, dtype=float)
    pi0 = np.asarray(pi0_counts, dtype=float)
    if data.sum() <= 0.0 or dvcs.sum() <= 0.0 or pi0.sum() <= 0.0:
        return TemplateMorphology(success=False)
    total = float(data.sum())

    def objective(x: np.ndarray) -> float:
        n_pi0, n_dvcs = np.exp(x[:2])
        pi0_shape = morph_template_histogram(pi0, edges, x[2], x[3])
        dvcs_shape = morph_template_histogram(dvcs, edges, x[4], x[5])
        pi0_shape /= max(float(pi0_shape.sum()), 1e-12)
        dvcs_shape /= max(float(dvcs_shape.sum()), 1e-12)
        model = n_pi0 * pi0_shape + n_dvcs * dvcs_shape
        # Weak regularization prevents pathological over-morphing in sparse categories.
        penalty = 0.5 * ((x[2] / 0.02) ** 2 + (x[4] / 0.02) ** 2 +
                         (x[3] / 0.025) ** 2 + (x[5] / 0.025) ** 2)
        return 0.5 * poisson_deviance(data, model) + penalty

    x0 = np.asarray([math.log(max(0.25 * total, 1.0)), math.log(max(0.75 * total, 1.0)),
                     0.0, 0.004, 0.0, 0.004], dtype=float)
    bounds = [
        (math.log(1e-6), math.log(max(10.0 * total, 1.0))),
        (math.log(1e-6), math.log(max(10.0 * total, 1.0))),
        (-0.04, 0.04), (0.0, 0.05), (-0.04, 0.04), (0.0, 0.05),
    ]
    result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)
    x = result.x
    pi0_shape = morph_template_histogram(pi0, edges, x[2], x[3]); pi0_shape /= max(pi0_shape.sum(), 1e-12)
    dvcs_shape = morph_template_histogram(dvcs, edges, x[4], x[5]); dvcs_shape /= max(dvcs_shape.sum(), 1e-12)
    model = math.exp(x[0]) * pi0_shape + math.exp(x[1]) * dvcs_shape
    return TemplateMorphology(float(x[2]), float(x[3]), float(x[4]), float(x[5]),
                              bool(result.success), poisson_deviance(data, model),
                              max(int(np.count_nonzero(data + model)) - 6, 0))


def fit_two_template_histograms(
    data_counts: np.ndarray,
    dvcs_counts: np.ndarray,
    pi0_counts: np.ndarray,
    edges: np.ndarray,
    morphology: Optional[TemplateMorphology] = None,
) -> FitSummary:
    data_counts = np.asarray(data_counts, dtype=float)
    dvcs_counts = np.asarray(dvcs_counts, dtype=float)
    pi0_counts = np.asarray(pi0_counts, dtype=float)

    if data_counts.sum() <= 0 or dvcs_counts.sum() <= 0 or pi0_counts.sum() <= 0:
        zeros = np.zeros_like(data_counts, dtype=float)
        return FitSummary(False, math.nan, math.nan, math.nan, math.nan, math.nan, 0,
                          "Empty data or template histogram", data_counts.astype(float), zeros,
                          zeros, zeros, edges)
    # endif

    morphology = morphology or TemplateMorphology()
    dvcs_morphed = morph_template_histogram(dvcs_counts, edges, morphology.dvcs_shift, morphology.dvcs_sigma)
    pi0_morphed = morph_template_histogram(pi0_counts, edges, morphology.pi0_shift, morphology.pi0_sigma)
    dvcs_shape = dvcs_morphed / max(float(dvcs_morphed.sum()), 1e-12)
    pi0_shape = pi0_morphed / max(float(pi0_morphed.sum()), 1e-12)
    total = float(data_counts.sum())

    def objective(log_yields: np.ndarray) -> float:
        yields = np.exp(log_yields)
        model = yields[0] * pi0_shape + yields[1] * dvcs_shape
        return 0.5 * poisson_deviance(data_counts, model)

    initial = np.log([max(0.25 * total, 1.0), max(0.75 * total, 1.0)])
    result = minimize(objective, initial, method="L-BFGS-B",
                      bounds=[(math.log(1e-6), math.log(max(10.0 * total, 1.0)))] * 2)
    yields = np.exp(result.x)
    n_pi0, n_dvcs = float(yields[0]), float(yields[1])
    model = n_pi0 * pi0_shape + n_dvcs * dvcs_shape

    # Numerical Hessian in log-yield coordinates; transform to yield errors.
    def hessian_2d(x: np.ndarray, step: float = 1e-4) -> np.ndarray:
        h = np.zeros((2, 2), dtype=float)
        f0 = objective(x)
        for i in range(2):
            ei = np.zeros(2); ei[i] = step
            h[i, i] = (objective(x + ei) - 2.0 * f0 + objective(x - ei)) / (step * step)
            for j in range(i + 1, 2):
                ej = np.zeros(2); ej[j] = step
                h[i, j] = h[j, i] = (
                    objective(x + ei + ej) - objective(x + ei - ej)
                    - objective(x - ei + ej) + objective(x - ei - ej)
                ) / (4.0 * step * step)
            # endfor
        # endfor
        return h

    try:
        covariance_log = np.linalg.inv(hessian_2d(result.x))
        errors = yields * np.sqrt(np.maximum(np.diag(covariance_log), 0.0))
        n_pi0_err, n_dvcs_err = float(errors[0]), float(errors[1])
    except np.linalg.LinAlgError:
        n_pi0_err = math.nan
        n_dvcs_err = math.nan
    # endtry

    return FitSummary(
        bool(result.success), n_pi0, n_pi0_err, n_dvcs, n_dvcs_err,
        poisson_deviance(data_counts, model), max(int(np.count_nonzero(data_counts + model)) - 2, 0),
        str(result.message), data_counts.astype(float), model,
        n_pi0 * pi0_shape, n_dvcs * dvcs_shape, edges,
    )


def efficiency_and_error(n_pass: float, pass_err: float, n_fail: float, fail_err: float) -> Tuple[float, float]:
    total = n_pass + n_fail
    if not (math.isfinite(total) and total > 0.0):
        return math.nan, math.nan
    # endif
    efficiency = n_pass / total
    dp = n_fail / (total * total)
    df = -n_pass / (total * total)
    variance = dp * dp * pass_err * pass_err + df * df * fail_err * fail_err
    return efficiency, math.sqrt(max(variance, 0.0))


def scale_factor_and_error(ed: float, sed: float, em: float, sem: float) -> Tuple[float, float]:
    if not (math.isfinite(ed) and math.isfinite(em) and em > 0.0):
        return math.nan, math.nan
    # endif
    sf = ed / em
    variance = (sed / em) ** 2 + (ed * sem / (em * em)) ** 2
    return sf, math.sqrt(max(variance, 0.0))


def draw_fit_panel(ax: plt.Axes, title: str, fit: FitSummary) -> None:
    centers = 0.5 * (fit.edges[:-1] + fit.edges[1:])
    ax.errorbar(centers, fit.data_counts, yerr=np.sqrt(fit.data_counts), fmt="o", markersize=2.2, label="Data epgamma")
    ax.step(fit.edges[:-1], fit.model_counts, where="post", label="Total fit")
    ax.step(fit.edges[:-1], fit.pi0_counts, where="post", label="pi0 as epgamma")
    ax.step(fit.edges[:-1], fit.dvcs_counts, where="post", label="DVCS/BH")
    ax.set_title(title, fontsize=8)
    ratio = fit.deviance / fit.ndf if fit.ndf > 0 else math.nan
    ax.text(0.98, 0.96, f"Npi0={fit.n_pi0:.0f}\nD/ndf={ratio:.2f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=7)


def plot_aggregated_fits(path: Path, period_label: str, definitions: Sequence[BinDefinition],
                         fits: Sequence[FitSummary], detector: str, sector: int) -> None:
    selected = [(d, f) for d, f in zip(definitions, fits)
                if d.detector == detector and d.sector == sector]
    if not selected:
        return
    ncols = 3
    nrows = int(math.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 3.8 * nrows), squeeze=False, sharex=True)
    for ax, (definition, fit) in zip(axes.flat, selected):
        if detector == "FT":
            subtitle = f"E {definition.E_low:.2f}-{definition.E_high:.2f} GeV"
        else:
            subtitle = (f"E {definition.E_low:.2f}-{definition.E_high:.2f} GeV, "
                        f"theta {definition.theta_low_deg:.1f}-{definition.theta_high_deg:.1f} deg")
        draw_fit_panel(ax, subtitle, fit)
    for ax in axes.flat[len(selected):]:
        ax.set_visible(False)
    for ax in axes[-1, :]:
        if ax.get_visible(): ax.set_xlabel(r"$M_X^2(ep)$ (GeV$^2$)")
    for row in axes:
        if row[0].get_visible(): row[0].set_ylabel("Counts")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    category = "FT" if detector == "FT" else f"FD sector {sector}"
    fig.suptitle(f"{period_label}: {category} fail-sample template fits", fontsize=15)
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.965))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_fit(path: Path, title: str, fit: FitSummary) -> None:
    centers = 0.5 * (fit.edges[:-1] + fit.edges[1:])
    widths = np.diff(fit.edges)
    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    ax.errorbar(centers, fit.data_counts, yerr=np.sqrt(fit.data_counts), fmt="o",
                markersize=3.5, label="Data epgamma")
    ax.step(fit.edges[:-1], fit.model_counts, where="post", label="Total fit")
    ax.step(fit.edges[:-1], fit.pi0_counts, where="post", label="pi0 as epgamma")
    ax.step(fit.edges[:-1], fit.dvcs_counts, where="post", label="DVCS/BH")
    ax.set_xlabel(r"$M_X^2(ep)$ (GeV$^2$)")
    ax.set_ylabel(f"Counts / {np.mean(widths):.4g} GeV$^2$")
    ax.set_title(title)
    ax.legend(frameon=False)
    text = (
        f"Npi0 = {fit.n_pi0:.1f} +/- {fit.n_pi0_err:.1f}\n"
        f"NDVCS = {fit.n_dvcs:.1f} +/- {fit.n_dvcs_err:.1f}\n"
        f"deviance/ndf = {fit.deviance:.1f}/{fit.ndf}"
    )
    ax.text(0.98, 0.96, text, transform=ax.transAxes, ha="right", va="top")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_scale_factors(path: Path, period_label: str, rows: Sequence[EfficiencyRow], detector: str, sector: int) -> None:
    selected = [row for row in rows if row.detector == detector and row.sector == sector]
    if not selected:
        return
    # endif
    x = np.arange(len(selected), dtype=float)
    labels = []
    for row in selected:
        if detector == "FT":
            labels.append(f"{row.E_low:.2f}-{row.E_high:.2f}")
        else:
            labels.append(
                f"E {row.E_low:.2f}-{row.E_high:.2f}\n"
                f"th {row.theta_low_deg:.1f}-{row.theta_high_deg:.1f}"
            )
        # endif
    # endfor

    fig, ax = plt.subplots(figsize=(max(9.0, 0.65 * len(selected)), 5.5))
    ax.errorbar(x, [r.scale_factor for r in selected],
                yerr=[r.scale_factor_err for r in selected], fmt="o")
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel(r"$S_\gamma=\epsilon_{data}/\epsilon_{MC}$")
    ax.set_xlabel("Predicted probe bin")
    category = "FT" if detector == "FT" else f"FD sector {sector}"
    ax.set_title(f"{period_label}: photon-efficiency scale factor, {category}")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_efficiencies(path: Path, period_label: str, rows: Sequence[EfficiencyRow], detector: str, sector: int) -> None:
    selected = [row for row in rows if row.detector == detector and row.sector == sector]
    if not selected:
        return
    # endif
    x = np.arange(len(selected), dtype=float)
    fig, ax = plt.subplots(figsize=(max(9.0, 0.65 * len(selected)), 5.5))
    ax.errorbar(x - 0.08, [r.efficiency_data for r in selected],
                yerr=[r.efficiency_data_err for r in selected], fmt="o", label="Data")
    ax.errorbar(x + 0.08, [r.efficiency_mc for r in selected],
                yerr=[r.efficiency_mc_err for r in selected], fmt="s", label="AAOGEN MC")
    ax.set_ylim(0.0, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels([str(r.bin_id) for r in selected])
    ax.set_xlabel("Probe-bin ID")
    ax.set_ylabel("Tag-and-probe photon efficiency")
    category = "FT" if detector == "FT" else f"FD sector {sector}"
    ax.set_title(f"{period_label}: {category}")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)



@dataclass(frozen=True)
class FitVariableSpec:
    key: str
    attr: str
    label: str
    bins: int
    low: float
    high: float
    log_morph: bool = False


FIT_VARIABLES: Tuple[FitVariableSpec, ...] = (
    FitVariableSpec("Delta_phi", "delta_phi", r"$\Delta\phi$ (rad)", 100, 2.84159, 3.44159, False),
    FitVariableSpec("theta_gamma_gamma", "theta_gamma_gamma", r"$\theta_{\gamma\gamma}$ (rad)", 120, 0.0, 3.0, True),
    FitVariableSpec("pTmiss", "pTmiss", r"$p_T^{\mathrm{miss}}$ (GeV)", 125, 0.0, 0.5, True),
)


@dataclass
class ProjectionResult:
    spec: FitVariableSpec
    edges: np.ndarray
    data_counts: np.ndarray
    model_counts: np.ndarray
    pi0_counts: np.ndarray
    dvcs_counts: np.ndarray
    deviance: float
    ndf: int


@dataclass
class MultiFitSummary:
    success: bool
    n_pi0: float
    n_pi0_err: float
    n_dvcs: float
    n_dvcs_err: float
    fraction_pi0: float
    fraction_pi0_err: float
    deviance: float
    ndf: int
    message: str
    projections: Dict[str, ProjectionResult]


@dataclass
class MultiMorphology:
    params: Dict[str, Tuple[float, float, float, float]]
    success: bool = True
    deviance: float = math.nan
    ndf: int = 0


def fit_edges() -> Dict[str, np.ndarray]:
    return {spec.key: np.linspace(spec.low, spec.high, spec.bins + 1) for spec in FIT_VARIABLES}


def common_fit_mask(trials: TrialArrays) -> np.ndarray:
    mask = np.ones(trials.size(), dtype=bool)
    for spec in FIT_VARIABLES:
        values = np.asarray(getattr(trials, spec.attr), dtype=float)
        mask &= np.isfinite(values) & (values >= spec.low) & (values < spec.high)
    return mask


def bulk_variable_histograms(trials: TrialArrays, bin_ids: np.ndarray, n_probe_bins: int,
                             spec: FitVariableSpec, edges: np.ndarray,
                             require_common_support: bool = True) -> np.ndarray:
    values = np.asarray(getattr(trials, spec.attr), dtype=float)
    value_index = np.searchsorted(edges, values, side="right") - 1
    valid = (bin_ids >= 0) & (bin_ids < n_probe_bins) & np.isfinite(values)
    if require_common_support:
        valid &= common_fit_mask(trials)
    valid &= (value_index >= 0) & (value_index < spec.bins)
    flat = bin_ids[valid].astype(np.int64, copy=False) * spec.bins + value_index[valid].astype(np.int64, copy=False)
    counts = np.bincount(flat, weights=trials.weight[valid], minlength=n_probe_bins * spec.bins)
    return counts.reshape(n_probe_bins, spec.bins).astype(float, copy=False)


def morph_fit_template(counts: np.ndarray, edges: np.ndarray, shift: float, sigma: float,
                       log_morph: bool) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    if counts.sum() <= 0.0:
        return np.zeros_like(counts)
    if not log_morph:
        return morph_template_histogram(counts, edges, shift, sigma)
    centers = 0.5 * (edges[:-1] + edges[1:])
    epsilon = max(0.25 * float(np.mean(np.diff(edges))), 1.0e-5)
    z = np.log(np.maximum(centers - edges[0] + epsilon, epsilon))
    z_uniform = np.linspace(float(z.min()), float(z.max()), len(z))
    shape_z = np.interp(z_uniform, z, counts, left=0.0, right=0.0)
    dz = float(np.mean(np.diff(z_uniform))) if len(z_uniform) > 1 else 1.0
    shifted = np.interp(z_uniform - shift, z_uniform, shape_z, left=0.0, right=0.0)
    if sigma > 0.0 and dz > 0.0:
        shifted = gaussian_filter1d(shifted, sigma / dz, mode="constant", cval=0.0)
    result = np.interp(z, z_uniform, shifted, left=0.0, right=0.0)
    result = np.clip(result, 0.0, None)
    return result * (counts.sum() / result.sum()) if result.sum() > 0.0 else np.zeros_like(counts)


def fit_category_multi_morphology(data_hists: Mapping[str, np.ndarray],
                                  dvcs_hists: Mapping[str, np.ndarray],
                                  pi0_hists: Mapping[str, np.ndarray],
                                  edges_by_var: Mapping[str, np.ndarray],
                                  disabled: bool = False) -> MultiMorphology:
    if disabled:
        return MultiMorphology({spec.key: (0.0, 0.0, 0.0, 0.0) for spec in FIT_VARIABLES})
    active = [spec for spec in FIT_VARIABLES if data_hists[spec.key].sum() > 0 and
              dvcs_hists[spec.key].sum() > 0 and pi0_hists[spec.key].sum() > 0]
    if len(active) < 2:
        return MultiMorphology({spec.key: (0.0, 0.0, 0.0, 0.0) for spec in FIT_VARIABLES}, False)
    # x = logit(f_pi0), then pi0 shift/sigma and DVCS shift/sigma for every variable.
    x0 = [math.log(0.25 / 0.75)]
    bounds = [(-8.0, 8.0)]
    for spec in active:
        if spec.log_morph:
            x0 += [0.0, 0.03, 0.0, 0.03]
            bounds += [(-0.35, 0.35), (0.0, 0.40), (-0.35, 0.35), (0.0, 0.40)]
        else:
            width = float(np.mean(np.diff(edges_by_var[spec.key])))
            x0 += [0.0, 1.0 * width, 0.0, 1.0 * width]
            bounds += [(-12*width, 12*width), (0.0, 12*width), (-12*width, 12*width), (0.0, 12*width)]
    def objective(x: np.ndarray) -> float:
        f = 1.0 / (1.0 + math.exp(-float(x[0])))
        total = 0.0
        pos = 1
        for spec in active:
            ps, psg, ds, dsg = map(float, x[pos:pos+4]); pos += 4
            edges = edges_by_var[spec.key]
            pshape = morph_fit_template(pi0_hists[spec.key], edges, ps, psg, spec.log_morph)
            dshape = morph_fit_template(dvcs_hists[spec.key], edges, ds, dsg, spec.log_morph)
            pshape /= max(pshape.sum(), 1e-12); dshape /= max(dshape.sum(), 1e-12)
            n = float(data_hists[spec.key].sum())
            total += 0.5 * poisson_deviance(data_hists[spec.key], n * (f * pshape + (1-f) * dshape))
            # Same weak nuisance regularization philosophy as the exclusivity study.
            if spec.log_morph:
                total += 0.5*((ps/0.15)**2 + (ds/0.15)**2 + (psg/0.20)**2 + (dsg/0.20)**2)
            else:
                w = float(np.mean(np.diff(edges)))
                total += 0.5*((ps/(4*w))**2 + (ds/(4*w))**2 + (psg/(5*w))**2 + (dsg/(5*w))**2)
        return total
    result = minimize(objective, np.asarray(x0), method="L-BFGS-B", bounds=bounds, options={"maxiter": 800})
    params = {spec.key: (0.0, 0.0, 0.0, 0.0) for spec in FIT_VARIABLES}
    pos = 1
    for spec in active:
        params[spec.key] = tuple(map(float, result.x[pos:pos+4])); pos += 4
    ndf = sum(int(np.count_nonzero(data_hists[s.key])) for s in active) - len(result.x)
    return MultiMorphology(params, bool(result.success), 2.0 * float(result.fun), max(ndf, 0))


def fit_multi_projection_bin(data_hists: Mapping[str, np.ndarray],
                             dvcs_hists: Mapping[str, np.ndarray],
                             pi0_hists: Mapping[str, np.ndarray],
                             edges_by_var: Mapping[str, np.ndarray],
                             morphology: MultiMorphology, total_candidates: float,
                             min_counts: int) -> MultiFitSummary:
    active = [s for s in FIT_VARIABLES if data_hists[s.key].sum() >= min_counts and
              dvcs_hists[s.key].sum() >= min_counts and pi0_hists[s.key].sum() >= min_counts]
    if len(active) < 2 or total_candidates <= 0.0:
        return MultiFitSummary(False, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan,
                               math.nan, 0, "Insufficient common-support statistics", {})
    shapes = {}
    for spec in active:
        ps, psg, ds, dsg = morphology.params.get(spec.key, (0.0,0.0,0.0,0.0))
        pshape = morph_fit_template(pi0_hists[spec.key], edges_by_var[spec.key], ps, psg, spec.log_morph)
        dshape = morph_fit_template(dvcs_hists[spec.key], edges_by_var[spec.key], ds, dsg, spec.log_morph)
        pshape /= max(pshape.sum(), 1e-12); dshape /= max(dshape.sum(), 1e-12)
        shapes[spec.key] = (pshape, dshape)
    def objective_scalar(f: float) -> float:
        f = min(max(float(f), 1e-6), 1.0-1e-6)
        total = 0.0
        for spec in active:
            pshape, dshape = shapes[spec.key]
            n = float(data_hists[spec.key].sum())
            total += 0.5 * poisson_deviance(data_hists[spec.key], n*(f*pshape + (1-f)*dshape))
        return total
    result = minimize(lambda x: objective_scalar(float(x[0])), np.asarray([0.25]),
                      method="L-BFGS-B", bounds=[(1e-6, 1.0-1e-6)])
    f = float(result.x[0])
    h = min(1e-3, 0.25*min(f,1-f))
    curvature = (objective_scalar(f+h)-2*objective_scalar(f)+objective_scalar(f-h))/(h*h) if h>0 else math.nan
    ferr = math.sqrt(1.0/curvature) if math.isfinite(curvature) and curvature>0 else math.nan
    projections = {}
    dev = 0.0; npoints = 0
    for spec in active:
        pshape, dshape = shapes[spec.key]
        n = float(data_hists[spec.key].sum())
        pc = n*f*pshape; dc = n*(1-f)*dshape; model = pc+dc
        one_dev = poisson_deviance(data_hists[spec.key], model)
        projections[spec.key] = ProjectionResult(spec, edges_by_var[spec.key], data_hists[spec.key], model, pc, dc, one_dev, max(int(np.count_nonzero(data_hists[spec.key]))-1,0))
        dev += one_dev; npoints += int(np.count_nonzero(data_hists[spec.key]))
    npi0 = f*total_candidates; ndvcs=(1-f)*total_candidates
    return MultiFitSummary(bool(result.success), npi0, ferr*total_candidates, ndvcs, ferr*total_candidates,
                           f, ferr, dev, max(npoints-1,0), str(result.message), projections)


def draw_projection_panel(ax, title: str, fit: MultiFitSummary, variable_key: str) -> None:
    projection = fit.projections.get(variable_key)
    if projection is None:
        ax.text(0.5,0.5,"insufficient statistics",ha="center",va="center",transform=ax.transAxes); ax.set_title(title,fontsize=8); return
    centers=0.5*(projection.edges[:-1]+projection.edges[1:])
    ax.errorbar(centers, projection.data_counts, yerr=np.sqrt(projection.data_counts), fmt="o", markersize=2.2, label="Data epgamma")
    ax.step(projection.edges[:-1], projection.model_counts, where="post", label="Total fit")
    ax.step(projection.edges[:-1], projection.pi0_counts, where="post", label="pi0 as epgamma")
    ax.step(projection.edges[:-1], projection.dvcs_counts, where="post", label="DVCS/BH")
    ratio=fit.deviance/fit.ndf if fit.ndf>0 else math.nan
    ax.text(0.98,0.96,f"fpi0={fit.fraction_pi0:.3f}\nD/ndf={ratio:.2f}",transform=ax.transAxes,ha="right",va="top",fontsize=7)
    ax.set_title(title,fontsize=8)


def plot_aggregated_multi_fits(path: Path, period_label: str, definitions: Sequence[BinDefinition],
                               fits: Sequence[MultiFitSummary], detector: str, sector: int,
                               variable_key: str) -> None:
    selected=[(d,f) for d,f in zip(definitions,fits) if d.detector==detector and d.sector==sector]
    if not selected: return
    ncols=3; nrows=int(math.ceil(len(selected)/ncols))
    fig,axes=plt.subplots(nrows,ncols,figsize=(15,3.8*nrows),squeeze=False,sharex=True)
    for ax,(d,f) in zip(axes.flat,selected):
        title=f"E {d.E_low:.2f}-{d.E_high:.2f} GeV" if detector=="FT" else f"E {d.E_low:.2f}-{d.E_high:.2f} GeV, theta {d.theta_low_deg:.1f}-{d.theta_high_deg:.1f} deg"
        draw_projection_panel(ax,title,f,variable_key)
    for ax in axes.flat[len(selected):]: ax.set_visible(False)
    spec=next(s for s in FIT_VARIABLES if s.key==variable_key)
    for ax in axes[-1,:]:
        if ax.get_visible(): ax.set_xlabel(spec.label)
    for row in axes:
        if row[0].get_visible(): row[0].set_ylabel("Counts")
    handles,labels=axes.flat[0].get_legend_handles_labels()
    category="FT" if detector=="FT" else f"FD sector {sector}"
    fig.suptitle(f"{period_label}: {category} epgamma multi-projection fits — {variable_key}",fontsize=15)
    if handles: fig.legend(handles,labels,loc="upper center",ncol=4,frameon=False,bbox_to_anchor=(0.5,0.965))
    fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(path,dpi=180); plt.close(fig)

def plot_category_fit_summary(path: Path, period_label: str,
                              definition: BinDefinition,
                              fit: MultiFitSummary) -> None:
    """Write one compact 1x3 fit canvas for a detector category."""
    fig, axes = plt.subplots(1, len(FIT_VARIABLES), figsize=(18, 5.2), squeeze=False)
    category = "FT" if definition.detector == "FT" else f"FD sector {definition.sector}"
    for ax, spec in zip(axes[0], FIT_VARIABLES):
        draw_projection_panel(ax, spec.label, fit, spec.key)
        ax.set_xlabel(spec.label)
        ax.set_ylabel("Counts")
    handles, labels = axes[0][0].get_legend_handles_labels()
    ratio = fit.deviance / fit.ndf if fit.ndf > 0 else math.nan
    fig.suptitle(
        f"{period_label}: {category} integrated epgamma template fit  "
        f"(fpi0={fit.fraction_pi0:.4f}, deviance/ndf={ratio:.2f})",
        fontsize=15,
    )
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=4,
                   frameon=False, bbox_to_anchor=(0.5, 0.94))
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_integrated_efficiency_summary(path: Path, period_label: str,
                                       rows: Sequence[EfficiencyRow]) -> None:
    """Compare all seven integrated detector categories in one canvas."""
    ordered = sorted(rows, key=lambda r: (0 if r.detector == "FT" else 1, r.sector))
    if not ordered:
        return
    labels = ["FT" if r.detector == "FT" else f"FD S{r.sector}" for r in ordered]
    x = np.arange(len(ordered), dtype=float)
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    axes[0].errorbar(x - 0.08, [r.efficiency_data for r in ordered],
                     yerr=[r.efficiency_data_err for r in ordered],
                     fmt="o", label="Data")
    axes[0].errorbar(x + 0.08, [r.efficiency_mc for r in ordered],
                     yerr=[r.efficiency_mc_err for r in ordered],
                     fmt="s", label="AAOGEN MC")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].set_ylabel("Tag-and-probe photon efficiency")
    axes[0].legend(frameon=False)
    axes[0].set_title(f"{period_label}: integrated photon-efficiency study")

    axes[1].errorbar(x, [r.scale_factor for r in ordered],
                     yerr=[r.scale_factor_err for r in ordered], fmt="o")
    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1].set_ylabel(r"$S_\gamma=\epsilon_{data}/\epsilon_{MC}$")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels)
    axes[1].set_xlabel("Photon detector category")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _fraction_sequence(cutflow: Mapping[str, int], stages: Sequence[str]) -> np.ndarray:
    first = max(int(cutflow.get(stages[0], 0)), 1)
    return np.asarray([int(cutflow.get(stage, 0))/first for stage in stages], dtype=float)


def plot_pass_cutflow(path: Path, period_label: str, data_diag: Mapping[str, object], mc_diag: Mapping[str, object]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 10))
    for ax, group, stages, title in ((axes[0], "event_cutflow", PASS_EVENT_STAGES, "Native eppi0 event cut flow"),
                                     (axes[1], "probe_cutflow", PASS_PROBE_STAGES, "Directed tag-and-probe cut flow")):
        x=np.arange(len(stages)); width=0.38
        for offset, diag, label in ((-width/2, data_diag, "Data"),(width/2, mc_diag, "AAOGEN MC")):
            vals=_fraction_sequence(diag[group], stages)
            ax.bar(x+offset, vals, width, label=label)
            for xx,v,stage in zip(x+offset, vals, stages):
                ax.text(xx, max(v,1e-6)*1.08, f"{diag[group][stage]:,}", rotation=90, ha="center", va="bottom", fontsize=7)
        ax.set_yscale("log"); ax.set_ylim(1e-6, 2.5); ax.set_ylabel("Fraction of initial sample")
        ax.set_xticks(x); ax.set_xticklabels([v.replace("_"," ") for v in stages], rotation=35, ha="right")
        ax.set_title(title); ax.legend(frameon=False)
    fig.suptitle(f"{period_label}: passing-sample audit")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def plot_pass_diagnostic_histograms(path: Path, period_label: str, data_diag: Mapping[str, object], mc_diag: Mapping[str, object], keys: Sequence[str], title: str) -> None:
    ncols=3; nrows=math.ceil(len(keys)/ncols)
    fig, axes=plt.subplots(nrows,ncols,figsize=(15,4.2*nrows)); axes=np.atleast_1d(axes).ravel()
    for ax,key in zip(axes,keys):
        for diag,label in ((data_diag,"Data"),(mc_diag,"AAOGEN MC")):
            h=diag["histograms"][key]; edges=np.asarray(h["edges"]); counts=np.asarray(h["counts"],dtype=float)
            if counts.sum()>0: counts/=counts.sum()
            ax.stairs(counts,edges,label=label)
        ax.set_xlabel(key.replace("_"," ")); ax.set_ylabel("Unit-normalized counts"); ax.legend(frameon=False)
    for ax in axes[len(keys):]: ax.axis("off")
    fig.suptitle(f"{period_label}: {title}"); fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)


def write_pass_audit(period_dir: Path, period_label: str, data_diag: Mapping[str, object], mc_diag: Mapping[str, object]) -> None:
    audit_dir=period_dir/"pass_sample_audit"; audit_dir.mkdir(parents=True,exist_ok=True)
    with open(audit_dir/"pass_sample_cutflow.json","w",encoding="utf-8") as h: json.dump({"data":data_diag,"mc":mc_diag},h,indent=2)
    plot_pass_cutflow(audit_dir/"cutflow.png",period_label,data_diag,mc_diag)
    plot_pass_diagnostic_histograms(audit_dir/"event_selection_variables.png",period_label,data_diag,mc_diag,
        ["Mx2","pTmiss","abs_Delta_phi_minus_pi"],"broad ep-pi0 selection")
    plot_pass_diagnostic_histograms(audit_dir/"native_photon_reconstruction.png",period_label,data_diag,mc_diag,
        ["closure","transverse_mismatch","energy_fraction_g1","ambiguity","n_solutions","g1_E","g2_E","g1_theta_deg","g2_theta_deg"],"native photon recovery")
    plot_pass_diagnostic_histograms(audit_dir/"directed_probe_matching.png",period_label,data_diag,mc_diag,
        ["pred_E","pred_theta_deg","pred_m2","pred_E_minus_p","match_angle_deg","relative_E_residual"],"predicted-probe matching")


def process_period(period: PeriodConfig, args_dict: Mapping[str, object]) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    period_start = time.perf_counter()
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    aggregate_fit_dir = period_dir / "epgamma_template_fits"
    plot_dir = period_dir / "plots"
    aggregate_fit_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    pass_data, pass_data_diag = read_pass_trials(period.epgg_data, period.beam_energy_GeV, args, "data")
    pass_mc, pass_mc_diag = read_pass_trials(period.epgg_mc, period.beam_energy_GeV, args, "mc")
    write_pass_audit(period_dir, period.label, pass_data_diag, pass_mc_diag)
    if args.skip_efficiency_extraction:
        metadata = {"period": asdict(period), "pass_sample_audit": {"data": pass_data_diag, "mc": pass_mc_diag},
                    "runtime_seconds": time.perf_counter() - period_start}
        with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2)
        return period.key, [], metadata

    fail_data = read_fail_trials(period.epg_data, period.beam_energy_GeV, args)
    fail_dvcs_mc = read_fail_trials(period.dvcs_mc, period.beam_energy_GeV, args)
    fail_pi0_mc = read_fail_trials(period.pi0_as_epg_mc, period.beam_energy_GeV, args)

    log(
        f"{period.label}: pass data={pass_data.size():,}, pass MC={pass_mc.size():,}, "
        f"fail data candidates={fail_data.size():,}, DVCS template={fail_dvcs_mc.size():,}, "
        f"pi0 template={fail_pi0_mc.size():,}"
    )

    bins = build_bins(pass_mc, fail_pi0_mc, args)
    rows: List[EfficiencyRow] = []
    n_bins = len(bins)
    edges_by_var = fit_edges()

    log(f"{period.label}: assigning trials to {n_bins} integrated detector categories")
    pass_data_ids = assign_bin_ids(pass_data, bins)
    pass_mc_ids = assign_bin_ids(pass_mc, bins)
    fail_data_ids = assign_bin_ids(fail_data, bins)
    fail_dvcs_ids = assign_bin_ids(fail_dvcs_mc, bins)
    fail_pi0_ids = assign_bin_ids(fail_pi0_mc, bins)

    pass_data_counts = weighted_counts_by_bin(pass_data, pass_data_ids, n_bins)
    pass_mc_counts = weighted_counts_by_bin(pass_mc, pass_mc_ids, n_bins)
    fail_data_total_counts = weighted_counts_by_bin(fail_data, fail_data_ids, n_bins)
    fail_pi0_counts = weighted_counts_by_bin(fail_pi0_mc, fail_pi0_ids, n_bins)

    data_hists = {spec.key: bulk_variable_histograms(fail_data, fail_data_ids, n_bins, spec, edges_by_var[spec.key]) for spec in FIT_VARIABLES}
    dvcs_hists = {spec.key: bulk_variable_histograms(fail_dvcs_mc, fail_dvcs_ids, n_bins, spec, edges_by_var[spec.key]) for spec in FIT_VARIABLES}
    pi0_hists = {spec.key: bulk_variable_histograms(fail_pi0_mc, fail_pi0_ids, n_bins, spec, edges_by_var[spec.key]) for spec in FIT_VARIABLES}

    morphology_by_category: Dict[Tuple[str, int], MultiMorphology] = {}
    for detector_name, sector in [("FT", 0)] + [("FD", value) for value in range(1, 7)]:
        indices = [d.bin_id for d in bins if d.detector == detector_name and d.sector == sector]
        if not indices:
            continue
        morphology_by_category[(detector_name, sector)] = fit_category_multi_morphology(
            {spec.key: data_hists[spec.key][indices].sum(axis=0) for spec in FIT_VARIABLES},
            {spec.key: dvcs_hists[spec.key][indices].sum(axis=0) for spec in FIT_VARIABLES},
            {spec.key: pi0_hists[spec.key][indices].sum(axis=0) for spec in FIT_VARIABLES},
            edges_by_var, args.disable_template_morphing,
        )

    fits: List[MultiFitSummary] = []
    for definition in bins:
        index = definition.bin_id
        n_pass_data = float(pass_data_counts[index])
        n_pass_mc = float(pass_mc_counts[index])
        n_fail_mc = float(fail_pi0_counts[index])
        fit = fit_multi_projection_bin(
            {spec.key: data_hists[spec.key][index] for spec in FIT_VARIABLES},
            {spec.key: dvcs_hists[spec.key][index] for spec in FIT_VARIABLES},
            {spec.key: pi0_hists[spec.key][index] for spec in FIT_VARIABLES},
            edges_by_var, morphology_by_category.get((definition.detector, definition.sector), MultiMorphology({})),
            float(fail_data_total_counts[index]), args.fit_min_counts,
        )
        fits.append(fit)
        n_fail_data = fit.n_pi0
        n_fail_data_err = fit.n_pi0_err
        pass_data_err = math.sqrt(max(n_pass_data, 0.0))
        pass_mc_err = math.sqrt(max(n_pass_mc, 0.0))
        fail_mc_err = math.sqrt(max(n_fail_mc, 0.0))
        eff_data, eff_data_err = efficiency_and_error(n_pass_data, pass_data_err, n_fail_data, n_fail_data_err)
        eff_mc, eff_mc_err = efficiency_and_error(n_pass_mc, pass_mc_err, n_fail_mc, fail_mc_err)
        sf, sf_err = scale_factor_and_error(eff_data, eff_data_err, eff_mc, eff_mc_err)
        rows.append(EfficiencyRow(
            period.key, period.label, definition.bin_id, definition.detector, definition.sector,
            definition.E_low, definition.E_high, definition.theta_low_deg, definition.theta_high_deg,
            n_pass_data, pass_data_err, n_fail_data, n_fail_data_err, eff_data, eff_data_err,
            n_pass_mc, pass_mc_err, n_fail_mc, fail_mc_err, eff_mc, eff_mc_err, sf, sf_err,
            fit.success, fit.deviance, fit.ndf,
        ))

    for definition, fit in zip(bins, fits):
        category_name = "FT" if definition.detector == "FT" else f"FD_sector_{definition.sector}"
        plot_category_fit_summary(
            aggregate_fit_dir / f"{category_name}.png",
            period.label, definition, fit,
        )

    plot_integrated_efficiency_summary(
        plot_dir / "integrated_efficiencies_and_scale_factors.png",
        period.label, rows,
    )

    metadata = {
        "period": asdict(period),
        "pass_sample_audit": {"data": pass_data_diag, "mc": pass_mc_diag},
        "counts": {
            "pass_data_directed_trials": pass_data.size(),
            "pass_mc_directed_trials": pass_mc.size(),
            "fail_data_candidates": fail_data.size(),
            "fail_dvcs_mc_candidates": fail_dvcs_mc.size(),
            "fail_pi0_mc_candidates": fail_pi0_mc.size(),
            "n_bins": len(rows),
        },
        "multi_projection_template_morphology": {
            f"{detector}_s{sector}": {
                "success": value.success, "deviance": value.deviance, "ndf": value.ndf,
                "parameters": {key: list(params) for key, params in value.params.items()},
            }
            for (detector, sector), value in morphology_by_category.items()
        },
        "runtime_seconds": time.perf_counter() - period_start,
    }
    with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    # endwith
    return period.key, [asdict(row) for row in rows], metadata


def apply_path_overrides(periods: List[PeriodConfig], args: argparse.Namespace) -> List[PeriodConfig]:
    overrides = [args.epg_data, args.dvcs_mc, args.pi0_as_epg_mc, args.epgg_data, args.epgg_mc]
    if any(value is not None for value in overrides):
        if len(periods) != 1:
            raise ValueError("Input-path overrides require exactly one --period.")
        # endif
        p = periods[0]
        periods[0] = PeriodConfig(
            p.key, p.label, p.beam_energy_GeV,
            args.epg_data or p.epg_data,
            args.dvcs_mc or p.dvcs_mc,
            args.pi0_as_epg_mc or p.pi0_as_epg_mc,
            args.epgg_data or p.epgg_data,
            args.epgg_mc or p.epgg_mc,
        )
    # endif
    return periods


def inspect_selected_branches(periods: Sequence[PeriodConfig]) -> None:
    for period in periods:
        print(f"\n[{period.label}]")
        for label, path in (
            ("epgamma data", period.epg_data),
            ("DVCS MC", period.dvcs_mc),
            ("pi0-as-epgamma MC", period.pi0_as_epg_mc),
            ("epgamma-gamma data", period.epgg_data),
            ("epgamma-gamma MC", period.epgg_mc),
        ):
            print(f"\n  {label}: {path}")
            try:
                print("    " + ", ".join(tree_keys(path)))
            except Exception as exc:
                print(f"    ERROR: {exc}")
            # endtry
        # endfor
    # endfor


def write_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    # endif
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def main() -> int:
    args = parse_args()
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    # endif
    args.workers = min(args.workers, MAX_WORKERS)

    selected_keys = set(args.period or [p.key for p in PERIODS])
    periods = [p for p in PERIODS if p.key in selected_keys]
    periods = apply_path_overrides(periods, args)

    if args.inspect_branches:
        inspect_selected_branches(periods)
        return 0
    # endif

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for period in periods:
        for path in (period.epg_data, period.dvcs_mc, period.pi0_as_epg_mc,
                     period.epgg_data, period.epgg_mc):
            require_file(path)
        # endfor
    # endfor

    effective_workers = min(args.workers, len(periods))
    log(
        f"Selected {len(periods)} period(s); using {effective_workers} period worker(s) "
        f"(requested={args.workers}, hard maximum={MAX_WORKERS})"
    )

    args_dict = vars(args).copy()
    all_rows: List[Dict[str, object]] = []
    metadata: Dict[str, object] = {}

    if effective_workers == 1:
        for period in periods:
            key, rows, period_metadata = process_period(period, args_dict)
            all_rows.extend(rows)
            metadata[key] = period_metadata
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=effective_workers) as executor:
            futures = {
                executor.submit(process_period, period, args_dict): period.key
                for period in periods
            }
            for future in concurrent.futures.as_completed(futures):
                key, rows, period_metadata = future.result()
                all_rows.extend(rows)
                metadata[key] = period_metadata
            # endfor
        # endwith
    # endif

    all_rows.sort(key=lambda row: (str(row["period"]), int(row["bin_id"])))
    write_rows(output_dir / "photon_efficiency_scale_factors.csv", all_rows)

    json_payload = {
        "schema_version": 2,
        "description": "Exploratory RGA pi0 tag-and-probe photon-efficiency data/MC scale factors.",
        "formula": "epsilon = N_pass / (N_pass + N_fail); S_gamma = epsilon_data / epsilon_mc",
        "arguments": args_dict,
        "period_metadata": metadata,
        "bins": all_rows,
    }
    with open(output_dir / "photon_efficiency_scale_factors.json", "w", encoding="utf-8") as handle:
        json.dump(json_payload, handle, indent=2, allow_nan=True)
    # endwith

    provenance = {
        "script": Path(__file__).name,
        "created_unix_time": time.time(),
        "tree_name": TREE_NAME,
        "periods": [asdict(period) for period in periods],
        "notes": [
            "Each epgamma-gamma event contributes two directed passing probes.",
            "Each fitted pi0-as-epgamma event contributes one failed probe.",
            "The extraction is integrated over photon energy and polar angle within FT and each FD sector.",
            "A complete data/MC passing-sample cut-flow audit is written for every period.",
            "No equal-efficiency approximation is used.",
            "Scale factors are not yet propagated to DVCS acceptance or pi0 migration.",
            "Period processing is parallelized with a hard maximum of seven workers.",
            "Trials are assigned to adaptive bins once and all three discriminator histograms are filled in bulk.",
        ],
    }
    with open(output_dir / "provenance.json", "w", encoding="utf-8") as handle:
        json.dump(provenance, handle, indent=2)
    # endwith

    log(f"Wrote {output_dir / 'photon_efficiency_scale_factors.csv'}")
    log(f"Wrote {output_dir / 'photon_efficiency_scale_factors.json'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise
    # endtry
