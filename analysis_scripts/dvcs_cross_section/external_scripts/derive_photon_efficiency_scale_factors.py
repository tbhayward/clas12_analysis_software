#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v29.py

RGA photon-reconstruction data/MC study based on exclusive pi0 control samples.

The ordinary reconstructed epgamma data sample is decomposed into genuine
DVCS/BH and pi0->gamma gamma events reconstructed with only one photon by
fitting the stored epgamma exclusivity distributions with DVCSGEN and
AAOGEN-as-epgamma templates. No missing-photon four-vector is constructed
before that fit.

The epgammagamma data and AAOGEN samples provide the corresponding population
where both pi0 decay photons are reconstructed. The photons are energy ordered:
the leading photon must satisfy the production DVCS threshold
(E_gamma1 >= 2 GeV), while the second reconstructed photon must satisfy the
loose threshold (E_gamma2 >= 0.4 GeV).

For each leading-photon detector category,

    epsilon_data = N_pi0(epgammagamma, data)
                   / [N_pi0(epgammagamma, data)
                      + N_pi0(epgamma, data from template fit)]

    epsilon_mc = N_pi0(epgammagamma, AAOGEN)
                 / [N_pi0(epgammagamma, AAOGEN)
                    + N_pi0(epgamma, AAOGEN)]

and S_gamma = epsilon_data / epsilon_mc.

This is directly an event-migration efficiency for recovering the additional
pi0 decay photon conditional on a reconstructed leading photon above the DVCS
threshold. Its data/MC ratio constrains photon-reconstruction mismodelling; it
is not by itself a direct measurement of the absolute probability to lose the
sole photon in a true DVCS event.
"""

from __future__ import annotations

import argparse
import importlib.util
import concurrent.futures
import csv
import hashlib
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
    from scipy.optimize import minimize, minimize_scalar
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
    "theta_cm": ("theta", "theta2", "theta_2"),
    "Emiss2": ("Emiss2", "E_miss2", "Emiss_sq"),
    "Mx2_2": ("Mx2_2", "Mx2_egamma", "Mx2_gamma", "Mx2_pi0", "Mx2_x2"),
    "fiducial_status": ("fiducial_status",),
    "t1": ("t1", "t", "minus_t"),
    "open_angle_ep2": ("open_angle_ep2", "open_angle_eg", "open_angle_ephoton"),
    "detector1": ("detector1", "proton_detector", "p_detector"),
    "Q2": ("Q2", "q2"),
    "W": ("W", "w"),
    "y": ("y", "inelasticity"),
    "z": ("z",),
    # Optional AAOGEN truth-level probe information.  These aliases are deliberately
    # broad; truth closure is skipped with an explicit reason when no complete set exists.
    "gen_g1_E": ("gen_p2_p", "mc_p2_p", "generated_gamma1_p", "gen_gamma1_p"),
    "gen_g1_theta": ("gen_p2_theta", "mc_p2_theta", "generated_gamma1_theta", "gen_gamma1_theta"),
    "gen_g1_phi": ("gen_p2_phi", "mc_p2_phi", "generated_gamma1_phi", "gen_gamma1_phi"),
    "gen_g2_E": ("gen_p3_p", "mc_p3_p", "generated_gamma2_p", "gen_gamma2_p"),
    "gen_g2_theta": ("gen_p3_theta", "mc_p3_theta", "generated_gamma2_theta", "gen_gamma2_theta"),
    "gen_g2_phi": ("gen_p3_phi", "mc_p3_phi", "generated_gamma2_phi", "gen_gamma2_phi"),
}


@dataclass
class TrialArrays:
    E: np.ndarray
    tag_E: np.ndarray
    theta_deg: np.ndarray
    phi_rad: np.ndarray
    detector: np.ndarray
    sector: np.ndarray
    mx2_ep: np.ndarray
    delta_phi: np.ndarray
    theta_gamma_gamma: np.ndarray
    pTmiss: np.ndarray
    theta_cm: np.ndarray
    Emiss2: np.ndarray
    Mx2: np.ndarray
    Mx2_2: np.ndarray
    proton_detector: np.ndarray
    weight: np.ndarray

    @staticmethod
    def empty() -> "TrialArrays":
        return TrialArrays(
            E=np.empty(0, dtype=np.float32),
            tag_E=np.empty(0, dtype=np.float32),
            theta_deg=np.empty(0, dtype=np.float32),
            phi_rad=np.empty(0, dtype=np.float32),
            detector=np.empty(0, dtype=np.int8),
            sector=np.empty(0, dtype=np.int8),
            mx2_ep=np.empty(0, dtype=np.float32),
            delta_phi=np.empty(0, dtype=np.float32),
            theta_gamma_gamma=np.empty(0, dtype=np.float32),
            pTmiss=np.empty(0, dtype=np.float32),
            theta_cm=np.empty(0, dtype=np.float32),
            Emiss2=np.empty(0, dtype=np.float32),
            Mx2=np.empty(0, dtype=np.float32),
            Mx2_2=np.empty(0, dtype=np.float32),
            proton_detector=np.empty(0, dtype=np.int8),
            weight=np.empty(0, dtype=np.float32),
        )

    def size(self) -> int:
        return int(self.E.size)


PASS_EVENT_STAGES = (
    "all_events", "basic_quality", "mx2", "pTmiss", "Delta_phi", "production_event_cuts",
    "native_inputs_finite", "longitudinal_energy_solution",
    "transverse_cone_solution", "detector_compatible_solution", "closure_pass",
)

PASS_PROBE_STAGES = (
    "two_trials_from_closure_events", "prediction_finite", "tag_energy_range", "probe_energy_range",
    "photon_like_E_minus_p", "photon_like_m2", "predicted_detector",
    "global_cuts", "angle_match", "energy_match", "final_pass",
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
            key: {"edges": edges.tolist(), "counts": np.zeros(len(edges)-1, dtype=float).tolist()}
            for key, edges in DIAGNOSTIC_HIST_EDGES.items()
        },
        "detector_pairs_after_closure": {},
        "final_trials_by_category": {"FT": 0, **{f"FD sector {i}": 0 for i in range(1, 7)}},
        "matching_scans": {},
        "mirror_diagnostics": {
            "events_with_one_solution": 0,
            "events_with_two_solutions": 0,
            "events_with_same_category_pair": 0,
            "events_with_different_category_pair": 0,
            "photon1_category_migration": [[0.0 for _ in range(7)] for _ in range(7)],
            "photon2_category_migration": [[0.0 for _ in range(7)] for _ in range(7)],
            "category_labels": ["FT", "FD S1", "FD S2", "FD S3", "FD S4", "FD S5", "FD S6"],
        },
    }


def diag_add_count(diag: Dict[str, object], group: str, stage: str, count: float) -> None:
    diag[group][stage] = float(diag[group].get(stage, 0.0)) + float(count)


def category_code(detector: np.ndarray, sector: np.ndarray) -> np.ndarray:
    """Map FT to 0 and FD sectors 1--6 to their sector number."""
    code = np.full(np.asarray(detector).shape, -1, dtype=np.int16)
    code[np.asarray(detector) == 0] = 0
    fd = np.asarray(detector) == 1
    code[fd] = np.asarray(sector, dtype=np.int16)[fd]
    return code


def weighted_count(mask: np.ndarray, weights: Optional[np.ndarray] = None) -> float:
    mask = np.asarray(mask, dtype=bool)
    if weights is None:
        return float(np.count_nonzero(mask))
    return float(np.sum(np.asarray(weights, dtype=float)[mask]))


def add_matching_scan_count(diag: Dict[str, object], key: str, detector: np.ndarray,
                            sector: np.ndarray, mask: np.ndarray,
                            weights: Optional[np.ndarray] = None) -> None:
    entry = diag["matching_scans"].setdefault(
        key, {"all": 0.0, "FT": 0.0, **{f"FD sector {i}": 0.0 for i in range(1, 7)}}
    )
    entry["all"] += weighted_count(mask, weights)
    entry["FT"] += weighted_count(mask & (detector == 0), weights)
    for i in range(1, 7):
        entry[f"FD sector {i}"] += weighted_count(mask & (detector == 1) & (sector == i), weights)


def diag_fill(diag: Dict[str, object], key: str, values: np.ndarray,
              mask: Optional[np.ndarray] = None,
              weights: Optional[np.ndarray] = None) -> None:
    values = np.asarray(values, dtype=float)
    selected = np.ones(values.shape, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    selected &= np.isfinite(values)
    values = values[selected]
    hist_weights = None if weights is None else np.asarray(weights, dtype=float)[selected]
    if values.size == 0:
        return
    edges = np.asarray(diag["histograms"][key]["edges"], dtype=float)
    counts, _ = np.histogram(values, bins=edges, weights=hist_weights)
    old = np.asarray(diag["histograms"][key]["counts"], dtype=float)
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
    parser.add_argument("--preflight-only", action="store_true",
                        help="Validate all inputs, audit observed-photon energy support, write preflight outputs, and exit before event processing.")
    parser.add_argument("--preflight-stable-age-min", type=float, default=10.0,
                        help="Minimum file age in minutes required by preflight to consider a ROOT file complete/stable; set negative to disable.")
    parser.add_argument("--preflight-min-entries", type=int, default=1000,
                        help="Minimum PhysicsEvents entries required for each one-photon loose input.")
    parser.add_argument("--min-fit-support-fraction", type=float, default=0.0,
                        help=("Optional minimum common fit-support fraction. The default 0 records coverage without rejecting categories; "
                              "use a positive value only for explicit validation studies."))
    parser.add_argument("--max-events", type=int, default=None,
                        help="Debug-only maximum events read from each tree.")

    parser.add_argument("--ft-theta-min", type=float, default=2.5)
    parser.add_argument("--ft-theta-max", type=float, default=5.0)
    parser.add_argument("--fd-theta-min", type=float, default=5.0)
    parser.add_argument("--fd-theta-max", type=float, default=35.0)
    parser.add_argument("--tag-E-min", type=float, default=0.40,
                        help="Minimum reconstructed tag-photon energy (GeV).")
    parser.add_argument("--tag-E-max", type=float, default=9.5,
                        help="Maximum reconstructed tag-photon energy (GeV).")
    parser.add_argument("--probe-E-min", type=float, default=2.0,
                        help="Minimum predicted probe-photon energy (GeV), matching the DVCS photon threshold.")
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
    parser.add_argument("--fit-support-containment", type=float, default=0.995,
                        help="MC-defined central containment used while calibrating shared template morphologies.")
    parser.add_argument("--constraint-variable-weight", type=float, default=0.35,
                        help="Legacy option retained for command-line compatibility; the exact exclusivity fitter now controls auxiliary profiling.")
    parser.add_argument("--exclusivity-max-shift-bins", type=float, default=12.0)
    parser.add_argument("--exclusivity-max-smear-bins", type=float, default=20.0)
    parser.add_argument("--exclusivity-shift-prior-bins", type=float, default=4.0)
    parser.add_argument("--exclusivity-smear-prior-bins", type=float, default=8.0)
    parser.add_argument("--exclusivity-dvcs-core-containment", type=float, default=0.90)
    parser.add_argument("--exclusivity-dvcs-fraction-containment", type=float, default=0.95)
    parser.add_argument("--exclusivity-pi0-core-containment", type=float, default=0.90)
    parser.add_argument("--exclusivity-pi0-fraction-containment", type=float, default=0.95)
    parser.add_argument("--exclusivity-outside-overshoot-penalty", type=float, default=0.25)
    parser.add_argument("--exclusivity-emiss2-mean-order-penalty", type=float, default=25.0)
    parser.add_argument("--disable-exclusivity-nuisance-penalties", action="store_true")
    parser.add_argument("--disable-fit-variant-systematic", action="store_true",
                        help="Disable leave-one-driver-out and expanded-fit model-systematic variants.")
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
    parser.add_argument("--enable-production-global-cuts", action="store_true",
                        help=("Opt in to the additional (-t1), z, and predicted electron-photon opening-angle cuts. "
                              "They are disabled by default because they were responsible for the order-of-magnitude "
                              "loss of one-photon candidates relative to the original study."))
    parser.add_argument("--disable-production-global-cuts", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--global-t-abs-max", type=float, default=1.0,
                        help="Production global requirement (-t1) < this value when t1 exists.")
    parser.add_argument("--global-z-min", type=float, default=0.65,
                        help="Common production requirement z > this value.")
    parser.add_argument("--global-open-angle-min-deg", type=float, default=5.0,
                        help="Production electron-photon opening-angle requirement.")
    parser.add_argument("--enable-global-dis-cuts", action="store_true",
                        help="Additionally require Q2>1 GeV2, W>2 GeV and y<0.8 when all branches exist.")
    parser.add_argument("--input-hash-bytes", type=int, default=8*1024*1024,
                        help="Bytes sampled from the start and end of each ROOT file for identity checks; 0 hashes the full file.")
    parser.add_argument("--allow-identical-period-inputs", action="store_true",
                        help="Allow nominally different period inputs with identical fingerprints.")
    parser.add_argument("--match-angle-scan-deg", default="1,2,3,4,5",
                        help="Comma-separated angular matching thresholds used for stability scans.")
    parser.add_argument("--match-energy-scan", default="0.20,0.35,0.50,1.00",
                        help="Comma-separated relative-energy matching thresholds used for stability scans.")
    parser.add_argument("--predicted-m2-scan", default="0.05,0.10,0.20,0.40,-1",
                        help="Comma-separated |m2_pred| thresholds; negative means no m2 cut.")
    parser.add_argument(
        "--enable-matching-scans", action="store_true",
        help=("Enable the expensive multidimensional matching-cut scan diagnostics. "
              "Disabled by default because they do not affect the nominal efficiency result."),
    )
    parser.add_argument("--mirror-policy", choices=("best", "half_weight", "unique_sector"), default="best",
                        help="Treatment of the two cone mirror solutions in passing probes.")
    parser.add_argument("--skip-truth-closure", action="store_true",
                        help="Skip optional AAOGEN truth-level closure even if truth branches exist.")

    parser.add_argument("--epg-data", default=None,
                        help="Override epgamma data path; valid only with one --period.")
    parser.add_argument("--dvcs-mc", default=None)
    parser.add_argument("--pi0-as-epg-mc", default=None)
    parser.add_argument("--epgg-data", default=None)
    parser.add_argument("--epgg-mc", default=None)
    return parser.parse_args()



def parse_float_list(text: str) -> List[float]:
    values = []
    for token in str(text).split(","):
        token = token.strip()
        if token:
            values.append(float(token))
    return values


def file_identity(path: str, sample_bytes: int = 8 * 1024 * 1024) -> Dict[str, object]:
    """Return a stable, inexpensive identity record for period-input validation."""
    resolved = Path(path).resolve()
    stat = resolved.stat()
    digest = hashlib.sha256()
    with resolved.open("rb") as handle:
        if sample_bytes == 0 or stat.st_size <= 2 * sample_bytes:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
            method = "sha256_full"
        else:
            digest.update(handle.read(sample_bytes))
            handle.seek(max(stat.st_size - sample_bytes, 0))
            digest.update(handle.read(sample_bytes))
            digest.update(str(stat.st_size).encode("ascii"))
            method = f"sha256_first_last_{sample_bytes}_bytes_plus_size"
    with uproot.open(f"{path}:{TREE_NAME}") as tree:
        entries = int(tree.num_entries)
        uuid = str(getattr(tree.file, "uuid", ""))
    return {
        "path": path, "resolved_path": str(resolved), "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns), "tree_entries": entries,
        "root_uuid": uuid, "fingerprint_method": method, "fingerprint": digest.hexdigest(),
    }


def validate_period_input_identities(periods: Sequence[PeriodConfig], args: argparse.Namespace) -> Dict[str, object]:
    records: Dict[str, object] = {}
    roles = ("epg_data", "dvcs_mc", "pi0_as_epg_mc", "epgg_data", "epgg_mc")
    duplicates = []
    by_role: Dict[str, Dict[Tuple[int, int, str], List[str]]] = {role: {} for role in roles}
    for period in periods:
        records[period.key] = {}
        for role in roles:
            ident = file_identity(getattr(period, role), args.input_hash_bytes)
            records[period.key][role] = ident
            key = (int(ident["size_bytes"]), int(ident["tree_entries"]), str(ident["fingerprint"]))
            by_role[role].setdefault(key, []).append(period.key)
    for role, groups in by_role.items():
        for keys in groups.values():
            if len(keys) > 1:
                duplicates.append({"role": role, "periods": keys})
    if duplicates and not args.allow_identical_period_inputs:
        detail = "; ".join(f"{d['role']}: {','.join(d['periods'])}" for d in duplicates)
        raise RuntimeError("Identical ROOT inputs detected across nominally different periods: " + detail +
                           ". Use --allow-identical-period-inputs only when intentional.")
    return {"files": records, "duplicates": duplicates}


def electron_photon_opening_deg(e_theta: np.ndarray, e_phi: np.ndarray,
                                g_theta: np.ndarray, g_phi: np.ndarray) -> np.ndarray:
    cosine = (np.sin(e_theta) * np.sin(g_theta) * np.cos(e_phi - g_phi)
              + np.cos(e_theta) * np.cos(g_theta))
    return np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))


def production_event_mask(arrays: Mapping[str, np.ndarray], resolved: Mapping[str, Optional[str]],
                          args: argparse.Namespace) -> np.ndarray:
    """Common event-level preselection shared by epgamma and epgammagamma.

    Only quantities with the same physical meaning in both reconstructed
    topologies are applied here.  Topology-specific exclusivity quantities
    (Mx2, pTmiss, Delta_phi, photon-pair matching, and cone closure) are not
    allowed to enter the nominal efficiency numerator or denominator.
    """
    n = len(next(iter(arrays.values())))
    mask = np.ones(n, dtype=bool)
    if args.disable_production_global_cuts:
        return mask
    if resolved.get("t1") is not None:
        t1 = finite_array(arrays, resolved.get("t1"))
        mask &= np.isfinite(t1) & ((-t1) < args.global_t_abs_max)
    if args.enable_global_dis_cuts:
        required = (resolved.get("Q2"), resolved.get("W"), resolved.get("y"))
        if all(branch is not None for branch in required):
            q2 = finite_array(arrays, resolved.get("Q2"))
            w = finite_array(arrays, resolved.get("W"))
            y = finite_array(arrays, resolved.get("y"))
            mask &= (
                np.isfinite(q2) & np.isfinite(w) & np.isfinite(y)
                & (q2 > 1.0) & (w > 2.0) & (y < 0.8)
            )
        else:
            raise RuntimeError("--enable-global-dis-cuts requested but Q2/W/y are not all present")
    return mask


def inspect_truth_closure_availability(path: str, args: argparse.Namespace) -> Dict[str, object]:
    logical = ["gen_g1_E", "gen_g1_theta", "gen_g1_phi", "gen_g2_E", "gen_g2_theta", "gen_g2_phi"]
    if args.skip_truth_closure:
        return {"available": False, "status": "disabled_by_argument", "resolved_branches": {}}
    try:
        resolved = resolve_branches(path, logical, optional=logical)
    except Exception as exc:
        return {"available": False, "status": f"branch_inspection_failed: {exc}", "resolved_branches": {}}
    available = all(resolved.get(name) is not None for name in logical)
    return {
        "available": available,
        "status": "available_for_future_direct_truth_matching" if available else "truth_four_vectors_not_present_in_this_tree",
        "resolved_branches": resolved,
        "note": "The current reconstructed AAOGEN eppi0 tree is audited for truth-branch availability. A direct generated-denominator closure requires these branches or a matched generated tree.",
    }

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
        tag_E=np.concatenate([chunk.tag_E for chunk in chunks]).astype(np.float32, copy=False),
        theta_deg=np.concatenate([chunk.theta_deg for chunk in chunks]).astype(np.float32, copy=False),
        phi_rad=np.concatenate([chunk.phi_rad for chunk in chunks]).astype(np.float32, copy=False),
        detector=np.concatenate([chunk.detector for chunk in chunks]).astype(np.int8, copy=False),
        sector=np.concatenate([chunk.sector for chunk in chunks]).astype(np.int8, copy=False),
        mx2_ep=np.concatenate([chunk.mx2_ep for chunk in chunks]).astype(np.float32, copy=False),
        delta_phi=np.concatenate([chunk.delta_phi for chunk in chunks]).astype(np.float32, copy=False),
        theta_gamma_gamma=np.concatenate([chunk.theta_gamma_gamma for chunk in chunks]).astype(np.float32, copy=False),
        pTmiss=np.concatenate([chunk.pTmiss for chunk in chunks]).astype(np.float32, copy=False),
        theta_cm=np.concatenate([chunk.theta_cm for chunk in chunks]).astype(np.float32, copy=False),
        Emiss2=np.concatenate([chunk.Emiss2 for chunk in chunks]).astype(np.float32, copy=False),
        Mx2=np.concatenate([chunk.Mx2 for chunk in chunks]).astype(np.float32, copy=False),
        Mx2_2=np.concatenate([chunk.Mx2_2 for chunk in chunks]).astype(np.float32, copy=False),
        proton_detector=np.concatenate([chunk.proton_detector for chunk in chunks]).astype(np.int8, copy=False),
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
    active_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, ...]:
    """Recover both detector-compatible cone-mirror solutions when available.

    The first solution index is the nominal best-score solution.  The second is
    the alternate mirror.  Storing both solutions allows the downstream
    ``best``, ``half_weight`` and ``unique_sector`` prescriptions to be applied
    without rerunning the cone reconstruction.
    """
    n_events = len(e_theta)
    g1_E = np.full(n_events, np.nan, dtype=float)
    g2_E = np.full(n_events, np.nan, dtype=float)
    solution_g1_theta = np.full((n_events, 2), np.nan, dtype=float)
    solution_g1_phi = np.full((n_events, 2), np.nan, dtype=float)
    solution_g2_theta = np.full((n_events, 2), np.nan, dtype=float)
    solution_g2_phi = np.full((n_events, 2), np.nan, dtype=float)
    solution_closure = np.full((n_events, 2), np.nan, dtype=float)
    solution_score = np.full((n_events, 2), np.nan, dtype=float)
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

    if active_mask is None:
        active_indices = range(n_events)
    else:
        active_indices = np.flatnonzero(np.asarray(active_mask, dtype=bool))

    for i in active_indices:
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
        g1_E[i], g2_E[i] = E1, E2

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
            four_residual = np.asarray([E1 + E2 - epi, *(E1 * n1 + E2 * n2 - pvec)], dtype=float)
            residual = float(np.linalg.norm(four_residual) / max(epi, 1.0e-12))
            mass2 = 2.0 * E1 * E2 * (1.0 - float(np.dot(n1, n2)))
            mass_residual = abs(mass2 - float(pi0_mass[i]) ** 2) / max(float(pi0_mass[i]) ** 2, 1.0e-12)
            score = math.hypot(residual, args.cone_mass_residual_weight * mass_residual)
            candidates.append((score, residual, mass_residual, th1, ph1, th2, ph2))

        candidates.sort(key=lambda item: item[0])
        n_solutions[i] = min(len(candidates), 2)
        if not candidates:
            continue
        detector_solution[i] = True
        for j, candidate in enumerate(candidates[:2]):
            solution_score[i, j] = candidate[0]
            solution_closure[i, j] = candidate[1]
            solution_g1_theta[i, j], solution_g1_phi[i, j] = candidate[3], candidate[4]
            solution_g2_theta[i, j], solution_g2_phi[i, j] = candidate[5], candidate[6]
        closure[i] = solution_closure[i, 0]
        ambiguity[i] = (solution_score[i, 1] - solution_score[i, 0]
                        if n_solutions[i] > 1 else math.inf)
        closure_pass[i] = solution_closure[i, 0] <= args.pass_photon_closure_max

    return (
        g1_E, g2_E, solution_g1_theta, solution_g1_phi,
        solution_g2_theta, solution_g2_phi, solution_closure, solution_score,
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
        "Delta_phi", "fiducial_status", "t1", "detector1", "Q2", "W", "y", "z",
    ]
    optional = ["Mx2_1", "Mx2", "pTmiss", "Delta_phi", "fiducial_status", "t1", "detector1", "Q2", "W", "y"]
    resolved = resolve_branches(path, logical, optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    chunks: List[TrialArrays] = []
    diag = empty_pass_diagnostics(sample_name, path)
    seen = 0
    angle_scan = parse_float_list(args.match_angle_scan_deg) if args.enable_matching_scans else []
    energy_scan = parse_float_list(args.match_energy_scan) if args.enable_matching_scans else []
    m2_scan = parse_float_list(args.predicted_m2_scan) if args.enable_matching_scans else []

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
        base = quality & production_event_mask(arrays, resolved, args)

        # Apply only the same leading-photon opening-angle requirement used for
        # the one-photon sample.  The stored first photon is energy ordered.
        opening1_deg = finite_array(arrays, resolved["open_angle_egamma1"])
        if not args.disable_production_global_cuts:
            base &= np.isfinite(opening1_deg) & (opening1_deg > args.global_open_angle_min_deg)

        diag_add_count(diag, "event_cutflow", "basic_quality", np.count_nonzero(quality))
        diag_add_count(diag, "event_cutflow", "mx2", np.count_nonzero(quality))
        diag_add_count(diag, "event_cutflow", "pTmiss", np.count_nonzero(quality))
        diag_add_count(diag, "event_cutflow", "Delta_phi", np.count_nonzero(quality))
        diag_add_count(diag, "event_cutflow", "production_event_cuts", np.count_nonzero(base))
        diag_fill(diag, "Mx2", mx2, quality)
        diag_fill(diag, "pTmiss", ptmiss_all, quality)
        diag_fill(diag, "abs_Delta_phi_minus_pi", np.abs(delta_phi_all-math.pi), quality)

        e_p = finite_array(arrays, resolved["e_p"]); e_theta = finite_array(arrays, resolved["e_theta"]); e_phi = finite_array(arrays, resolved["e_phi"])
        p_p = finite_array(arrays, resolved["p_p"]); p_theta = finite_array(arrays, resolved["p_theta"]); p_phi = finite_array(arrays, resolved["p_phi"])
        pi0_p = finite_array(arrays, resolved["pi0_p"]); pi0_theta = finite_array(arrays, resolved["pi0_theta"]); pi0_phi = finite_array(arrays, resolved["pi0_phi"])
        opening1 = np.deg2rad(finite_array(arrays, resolved["open_angle_egamma1"]))
        opening2 = np.deg2rad(finite_array(arrays, resolved["open_angle_egamma2"]))
        detector1_obs = finite_array(arrays, resolved["gamma1_detector_native"], default=-1.0).astype(int)
        detector2_obs = finite_array(arrays, resolved["gamma2_detector_native"], default=-1.0).astype(int)
        pi0_mass = finite_array(arrays, resolved["Mh_gammagamma"])
        mx2_ep = finite_array(arrays, resolved.get("Mx2_1"))
        theta_gg_all = np.full(n, np.nan, dtype=float)

        rec = reconstruct_native_eppi0_photons(
            e_theta, e_phi, pi0_p, pi0_theta, pi0_phi, opening1, opening2,
            pi0_mass, detector1_obs, detector2_obs, args, active_mask=base)
        (g1_E, g2_E, sol_g1_theta, sol_g1_phi, sol_g2_theta, sol_g2_phi,
         sol_closure, sol_score, closure, ambiguity, transverse_mismatch,
         energy_fraction_g1, n_solutions, input_finite, longitudinal_solution,
         transverse_solution, detector_solution, closure_pass) = rec

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
        diag_fill(diag, "g1_theta_deg", np.degrees(sol_g1_theta[:, 0]), base & closure_pass)
        diag_fill(diag, "g2_theta_deg", np.degrees(sol_g2_theta[:, 0]), base & closure_pass)

        pair_keys = np.char.add(np.char.add(detector1_obs.astype(str), "-"), detector2_obs.astype(str))
        for key in np.unique(pair_keys[base & closure_pass]):
            diag["detector_pairs_after_closure"][str(key)] = int(diag["detector_pairs_after_closure"].get(str(key), 0)) + int(np.count_nonzero((pair_keys == key) & base & closure_pass))

        # Mirror-category diagnostics.  Category code 0 is FT and 1--6 are FD sectors.
        det11, sec11 = classify_predicted_detector(np.degrees(sol_g1_theta[:, 0]), sol_g1_phi[:, 0], args)
        det21, sec21 = classify_predicted_detector(np.degrees(sol_g2_theta[:, 0]), sol_g2_phi[:, 0], args)
        det12, sec12 = classify_predicted_detector(np.degrees(sol_g1_theta[:, 1]), sol_g1_phi[:, 1], args)
        det22, sec22 = classify_predicted_detector(np.degrees(sol_g2_theta[:, 1]), sol_g2_phi[:, 1], args)
        cat11, cat21 = category_code(det11, sec11), category_code(det21, sec21)
        cat12, cat22 = category_code(det12, sec12), category_code(det22, sec22)
        one = base & closure_pass & (n_solutions == 1)
        two = base & closure_pass & (n_solutions >= 2)
        same_pair = two & (cat11 == cat12) & (cat21 == cat22)
        different_pair = two & ~same_pair
        md = diag["mirror_diagnostics"]
        md["events_with_one_solution"] += int(np.count_nonzero(one))
        md["events_with_two_solutions"] += int(np.count_nonzero(two))
        md["events_with_same_category_pair"] += int(np.count_nonzero(same_pair))
        md["events_with_different_category_pair"] += int(np.count_nonzero(different_pair))
        for a, b, matrix_key in ((cat11, cat12, "photon1_category_migration"),
                                 (cat21, cat22, "photon2_category_migration")):
            matrix = np.asarray(md[matrix_key], dtype=float)
            valid = two & (a >= 0) & (a <= 6) & (b >= 0) & (b <= 6)
            np.add.at(matrix, (a[valid], b[valid]), 1.0)
            md[matrix_key] = matrix.tolist()

        # Select one or both mirror geometries according to the requested policy.
        if args.mirror_policy == "best":
            variants = [(0, base & detector_solution & (n_solutions >= 1), np.ones(n, dtype=float))]
        elif args.mirror_policy == "half_weight":
            w0 = np.where(n_solutions >= 2, 0.5, 1.0)
            w1 = np.where(n_solutions >= 2, 0.5, 0.0)
            variants = [
                (0, base & detector_solution & (n_solutions >= 1), w0),
                (1, base & detector_solution & (n_solutions >= 2), w1),
            ]
        elif args.mirror_policy == "unique_sector":
            unique = base & detector_solution & ((n_solutions == 1) | same_pair)
            variants = [(0, unique, np.ones(n, dtype=float))]
        else:
            raise RuntimeError(f"Unsupported mirror policy: {args.mirror_policy}")

        for solution_index, solution_seed, solution_weight in variants:
            g1_theta = sol_g1_theta[:, solution_index]; g1_phi = sol_g1_phi[:, solution_index]
            g2_theta = sol_g2_theta[:, solution_index]; g2_phi = sol_g2_phi[:, solution_index]
            # The stored photons are energy ordered.  The nominal numerator is
            # the directly reconstructed two-photon category.  No inferred
            # missing-photon closure, E-p, m2, angular-match, or energy-match
            # requirement is imposed, because those conditions have no
            # counterpart in the one-photon denominator.
            seed = solution_seed
            lead_theta_deg = np.degrees(g1_theta)
            lead_detector, lead_sector = classify_predicted_detector(
                lead_theta_deg, g1_phi, args
            )
            proton_detector = finite_array(
                arrays, resolved.get("detector1"), default=-1.0
            ).astype(np.int16)

            finite = (
                seed
                & np.isfinite(g1_E) & np.isfinite(g2_E)
                & np.isfinite(g1_theta) & np.isfinite(g1_phi)
                & np.isfinite(g2_theta) & np.isfinite(g2_phi)
            )
            leading_energy = (
                finite
                & (g1_E >= args.probe_E_min)
                & (g1_E < args.probe_E_max)
            )
            second_energy = (
                leading_energy
                & (g2_E >= args.tag_E_min)
                & (g2_E < args.tag_E_max)
            )
            supported_topology = (
                ((proton_detector == 2) & (lead_detector == 0))
                | ((proton_detector == 2) & (lead_detector == 1))
                | ((proton_detector == 1) & (lead_detector == 1))
            )
            detector_ok = second_energy & (lead_detector >= 0) & supported_topology
            final_pass = detector_ok

            for stage, mask in (
                ("prediction_finite", finite),
                ("tag_energy_range", leading_energy),
                ("probe_energy_range", second_energy),
                ("photon_like_E_minus_p", second_energy),
                ("photon_like_m2", second_energy),
                ("predicted_detector", detector_ok),
                ("global_cuts", detector_ok),
                ("angle_match", detector_ok),
                ("energy_match", final_pass),
                ("final_pass", final_pass),
            ):
                diag_add_count(
                    diag, "probe_cutflow", stage,
                    weighted_count(mask, solution_weight),
                )

            for cat, cmask in [("FT", final_pass & (lead_detector == 0))] + [
                (
                    f"FD sector {j}",
                    final_pass & (lead_detector == 1) & (lead_sector == j),
                )
                for j in range(1, 7)
            ]:
                diag["final_trials_by_category"][cat] += weighted_count(
                    cmask, solution_weight
                )

            chunks.append(
                TrialArrays(
                    # Binning/category variables describe the actual leading
                    # photon shared by the epgamma and epgammagamma samples.
                    E=g1_E[final_pass].astype(np.float32, copy=False),
                    tag_E=g1_E[final_pass].astype(np.float32, copy=False),
                    theta_deg=lead_theta_deg[final_pass].astype(np.float32, copy=False),
                    phi_rad=g1_phi[final_pass].astype(np.float32, copy=False),
                    detector=lead_detector[final_pass].astype(np.int8, copy=False),
                    sector=lead_sector[final_pass].astype(np.int8, copy=False),
                    mx2_ep=mx2_ep[final_pass].astype(np.float32, copy=False),
                    delta_phi=delta_phi_all[final_pass].astype(np.float32, copy=False),
                    theta_gamma_gamma=theta_gg_all[final_pass].astype(np.float32, copy=False),
                    pTmiss=ptmiss_all[final_pass].astype(np.float32, copy=False),
                    theta_cm=np.full(np.count_nonzero(final_pass), np.nan, dtype=np.float32),
                    Emiss2=np.full(np.count_nonzero(final_pass), np.nan, dtype=np.float32),
                    Mx2=np.full(np.count_nonzero(final_pass), np.nan, dtype=np.float32),
                    Mx2_2=np.full(np.count_nonzero(final_pass), np.nan, dtype=np.float32),
                    proton_detector=proton_detector[final_pass].astype(np.int8, copy=False),
                    weight=solution_weight[final_pass].astype(np.float32, copy=False),
                )
            )
    result = concatenate_trials(chunks)
    log(f"Native eppi0 audit for {Path(path).name}: events={seen:,}, closure events={diag['event_cutflow']['closure_pass']:,.0f}, final directed trial weight={float(np.sum(result.weight)):,.1f}")
    return result, diag

def _new_one_photon_cutflow(label: str, path: str) -> Dict[str, object]:
    stages = [
        "all_tree_entries", "basic_quality", "production_event_cuts",
        "leading_photon_energy", "leading_photon_detector", "selected",
    ]
    return {
        "label": label,
        "path": path,
        "stages": {stage: 0.0 for stage in stages},
        "selected_by_category": {
            "FT": 0.0,
            **{f"FD sector {i}": 0.0 for i in range(1, 7)},
        },
    }


def _add_cutflow(diag: Dict[str, object], stage: str, mask_or_count) -> None:
    value = float(mask_or_count) if np.isscalar(mask_or_count) else float(np.count_nonzero(mask_or_count))
    diag["stages"][stage] = float(diag["stages"].get(stage, 0.0)) + value


def read_one_photon_candidates(
    path: str,
    beam_energy: float,
    args: argparse.Namespace,
    sample_label: str = "one-photon sample",
) -> Tuple[TrialArrays, Dict[str, object]]:
    """
    Read the ordinary reconstructed epgamma population directly.

    No missing-photon four-vector is constructed and no event is assumed to be
    pi0 before the template decomposition.  The measured leading photon is the
    photon that would enter the production DVCS candidate.  The same selection
    is applied to data, DVCSGEN, and AAOGEN-as-epgamma samples.
    """
    del beam_energy  # retained in the call signature for backward compatibility

    logical = [
        "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi",
        "g1_E", "g1_theta", "g1_phi", "g1_detector",
        "Mx2_1", "Delta_phi", "theta_gamma_gamma", "pTmiss",
        "theta_cm", "Emiss2", "Mx2", "Mx2_2",
        "fiducial_status", "t1", "open_angle_ep2", "proton_detector", "Q2", "W", "y", "z",
    ]
    optional = [
        "g1_detector", "theta_cm", "Emiss2", "Mx2", "Mx2_2",
        "fiducial_status", "t1", "open_angle_ep2", "proton_detector", "Q2", "W", "y",
    ]
    resolved = resolve_branches(path, logical, optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    chunks: List[TrialArrays] = []
    seen = 0
    diag = _new_one_photon_cutflow(sample_label, path)

    log(f"Reading {sample_label} from {path}")
    for arrays in uproot.iterate(
        f"{path}:{TREE_NAME}",
        expressions=expressions,
        step_size=args.step_size,
        library="np",
    ):
        n = len(next(iter(arrays.values())))
        if args.max_events is not None and seen >= args.max_events:
            break
        if args.max_events is not None and seen + n > args.max_events:
            keep = args.max_events - seen
            arrays = {key: value[:keep] for key, value in arrays.items()}
            n = keep
        seen += n
        _add_cutflow(diag, "all_tree_entries", n)

        quality = basic_quality_mask(arrays, resolved, args)
        production = quality & production_event_mask(arrays, resolved, args)
        if resolved.get("open_angle_ep2") is not None and not args.disable_production_global_cuts:
            opening = finite_array(arrays, resolved.get("open_angle_ep2"))
            production &= np.isfinite(opening) & (opening > args.global_open_angle_min_deg)
        _add_cutflow(diag, "basic_quality", quality)
        _add_cutflow(diag, "production_event_cuts", production)

        g_E = finite_array(arrays, resolved["g1_E"])
        g_theta = finite_array(arrays, resolved["g1_theta"])
        g_phi = finite_array(arrays, resolved["g1_phi"])
        theta_deg = np.degrees(g_theta)

        # The epgamma template fit represents the actual DVCS candidate sample.
        # Therefore the measured photon, not an inferred partner, must satisfy
        # the production photon threshold.
        energy_ok = (
            production
            & np.isfinite(g_E)
            & (g_E >= args.probe_E_min)
            & (g_E < args.probe_E_max)
            & np.isfinite(theta_deg)
            & np.isfinite(g_phi)
        )
        _add_cutflow(diag, "leading_photon_energy", energy_ok)

        detector, sector = classify_predicted_detector(theta_deg, g_phi, args)
        proton_detector = finite_array(
            arrays, resolved.get("proton_detector"), default=-1.0
        ).astype(np.int16)
        supported_topology = (
            ((proton_detector == 2) & (detector == 0))
            | ((proton_detector == 2) & (detector == 1))
            | ((proton_detector == 1) & (detector == 1))
        )
        detector_ok = energy_ok & (detector >= 0) & supported_topology
        _add_cutflow(diag, "leading_photon_detector", energy_ok & (detector >= 0))
        _add_cutflow(diag, "selected", detector_ok)

        for name, cmask in [("FT", detector_ok & (detector == 0))] + [
            (f"FD sector {j}", detector_ok & (detector == 1) & (sector == j))
            for j in range(1, 7)
        ]:
            diag["selected_by_category"][name] += float(np.count_nonzero(cmask))

        good = detector_ok
        chunks.append(
            TrialArrays(
                E=g_E[good].astype(np.float32, copy=False),
                tag_E=g_E[good].astype(np.float32, copy=False),
                theta_deg=theta_deg[good].astype(np.float32, copy=False),
                phi_rad=g_phi[good].astype(np.float32, copy=False),
                detector=detector[good].astype(np.int8, copy=False),
                sector=sector[good].astype(np.int8, copy=False),
                mx2_ep=finite_array(arrays, resolved["Mx2_1"])[good].astype(np.float32, copy=False),
                delta_phi=finite_array(arrays, resolved.get("Delta_phi"))[good].astype(np.float32, copy=False),
                theta_gamma_gamma=finite_array(arrays, resolved.get("theta_gamma_gamma"))[good].astype(np.float32, copy=False),
                pTmiss=finite_array(arrays, resolved.get("pTmiss"))[good].astype(np.float32, copy=False),
                theta_cm=finite_array(arrays, resolved.get("theta_cm"))[good].astype(np.float32, copy=False),
                Emiss2=finite_array(arrays, resolved.get("Emiss2"))[good].astype(np.float32, copy=False),
                Mx2=finite_array(arrays, resolved.get("Mx2"))[good].astype(np.float32, copy=False),
                Mx2_2=finite_array(arrays, resolved.get("Mx2_2"))[good].astype(np.float32, copy=False),
                proton_detector=proton_detector[good].astype(np.int8, copy=False),
                weight=np.ones(np.count_nonzero(good), dtype=np.float32),
            )
        )

    return concatenate_trials(chunks), diag


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



def weighted_sumw2_by_bin(trials: TrialArrays, bin_ids: np.ndarray, n_bins: int) -> np.ndarray:
    valid = (bin_ids >= 0) & (bin_ids < n_bins)
    return np.bincount(
        bin_ids[valid].astype(np.int64, copy=False),
        weights=np.square(trials.weight[valid]),
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
    fraction_driver: bool = True
    fit_weight: float = 1.0


FIT_VARIABLES: Tuple[FitVariableSpec, ...] = (
    # Use the validated exclusivity-template supports.  Broad distributions are
    # produced separately as shape diagnostics and are not used to drive the fit.
    FitVariableSpec("Delta_phi", "delta_phi", r"$\Delta\phi$ (rad)", 100, 2.84159, 3.44159, False, True, 1.0),
    FitVariableSpec("theta_gamma_gamma", "theta_gamma_gamma", r"$\theta_{\gamma\gamma}$ (rad)", 120, 0.0, 3.0, False, True, 1.0),
    FitVariableSpec("pTmiss", "pTmiss", r"$p_T^{\mathrm{miss}}$ (GeV)", 125, 0.0, 0.5, False, True, 1.0),
    FitVariableSpec("theta_cm", "theta_cm", r"$\theta_{p\gamma}^{\mathrm{CM}}$ (rad)", 100, 2.0, math.pi, False, False, 0.0),
    FitVariableSpec("Emiss2", "Emiss2", r"$E_{\mathrm{miss}}$ (GeV)", 100, -1.0, 2.0, False, False, 0.0),
    FitVariableSpec("Mx2", "Mx2", r"$M_x^2$ (GeV$^2$)", 100, -0.03, 0.03, False, False, 0.0),
    FitVariableSpec("Mx2_2", "Mx2_2", r"$M_{x2}^2$ (GeV$^2$)", 125, -1.0, 4.0, False, False, 0.0),
)
FRACTION_DRIVER_KEYS = tuple(spec.key for spec in FIT_VARIABLES if spec.fraction_driver)


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
    """Common finite population for all fraction-driving projections.

    Histogram-range cuts must not silently redefine the physical FAIL sample.
    Finite values outside the nominal plotting support are folded into the edge
    bins by ``bulk_variable_histograms``.
    """
    mask = np.ones(trials.size(), dtype=bool)
    for spec in FIT_VARIABLES:
        if spec.fraction_driver:
            mask &= np.isfinite(np.asarray(getattr(trials, spec.attr), dtype=float))
    return mask


def bulk_variable_histograms(trials: TrialArrays, bin_ids: np.ndarray, n_probe_bins: int,
                             spec: FitVariableSpec, edges: np.ndarray,
                             require_common_support: bool = True) -> np.ndarray:
    values = np.asarray(getattr(trials, spec.attr), dtype=float)
    value_index = np.searchsorted(edges, values, side="right") - 1
    valid = (bin_ids >= 0) & (bin_ids < n_probe_bins) & np.isfinite(values)
    valid &= (values >= spec.low) & (values < spec.high)
    if require_common_support:
        valid &= common_fit_mask(trials)
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


def mc_containment_mask(counts: np.ndarray, containment: float) -> np.ndarray:
    counts = np.asarray(counts, dtype=float)
    mask = np.ones(counts.size, dtype=bool)
    total = float(np.sum(counts))
    if total <= 0.0 or containment >= 0.999999:
        return mask
    cdf = np.cumsum(np.clip(counts, 0.0, None)) / total
    tail = 0.5 * (1.0 - containment)
    lo = int(np.searchsorted(cdf, tail, side="left"))
    hi = int(np.searchsorted(cdf, 1.0 - tail, side="left"))
    mask[:] = False
    mask[max(0, lo):min(counts.size, hi + 1)] = True
    if np.count_nonzero(mask) < 5:
        mask[:] = True
    return mask


def fit_category_multi_morphology(data_hists: Mapping[str, np.ndarray],
                                  dvcs_hists: Mapping[str, np.ndarray],
                                  pi0_hists: Mapping[str, np.ndarray],
                                  edges_by_var: Mapping[str, np.ndarray],
                                  disabled: bool = False, support_containment: float = 0.995,
                                  constraint_weight: float = 0.35) -> MultiMorphology:
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
            model = n * (f * pshape + (1-f) * dshape)
            support = mc_containment_mask(dvcs_hists[spec.key] + pi0_hists[spec.key], support_containment)
            weight = 1.0 if spec.fraction_driver else constraint_weight
            total += weight * 0.5 * poisson_deviance(data_hists[spec.key][support], model[support])
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
            weight = spec.fit_weight
            total += weight * 0.5 * poisson_deviance(data_hists[spec.key], n*(f*pshape + (1-f)*dshape))
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
    """Write fit projections with standardized signed-pull residual panels."""
    fig, axes = plt.subplots(2, len(FIT_VARIABLES), figsize=(4.7 * len(FIT_VARIABLES), 8.0), squeeze=False,
                             gridspec_kw={"height_ratios": [3.0, 1.2]}, sharex="col")
    category = "FT" if definition.detector == "FT" else f"FD sector {definition.sector}"
    for col, spec in enumerate(FIT_VARIABLES):
        ax = axes[0, col]
        draw_projection_panel(ax, spec.label, fit, spec.key)
        ax.set_ylabel("Counts")
        projection = fit.projections.get(spec.key)
        rax = axes[1, col]
        if projection is not None:
            centers = 0.5 * (projection.edges[:-1] + projection.edges[1:])
            sigma = np.sqrt(np.maximum(projection.model_counts, 1.0))
            pull = (projection.data_counts - projection.model_counts) / sigma
            rax.axhline(0.0, linewidth=1.0)
            rax.axhline(2.0, linestyle="--", linewidth=0.8)
            rax.axhline(-2.0, linestyle="--", linewidth=0.8)
            rax.plot(centers, pull, "o", markersize=2.2)
            rax.set_ylim(-6.0, 6.0)
            rax.text(0.98, 0.92, f"D/ndf={projection.deviance/max(projection.ndf,1):.2f}",
                     transform=rax.transAxes, ha="right", va="top", fontsize=7)
        else:
            rax.text(0.5, 0.5, "insufficient statistics", ha="center", va="center", transform=rax.transAxes)
        rax.set_xlabel(spec.label)
        rax.set_ylabel("Pull")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    ratio = fit.deviance / fit.ndf if fit.ndf > 0 else math.nan
    fig.suptitle(f"{period_label}: {category} integrated epgamma template fit  "
                 f"(fpi0={fit.fraction_pi0:.4f}, deviance/ndf={ratio:.2f})", fontsize=15)
    if handles:
        fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.94))
    fig.tight_layout(rect=(0, 0, 1, 0.89))
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


def plot_all_period_scale_factor_summary(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    period_order = [p.key for p in PERIODS]
    labels = ["FT"] + [f"FD S{i}" for i in range(1, 7)]
    x = np.arange(7, dtype=float)
    offsets = np.linspace(-0.18, 0.18, len(period_order))
    fig, ax = plt.subplots(figsize=(12, 6.8))
    for offset, key in zip(offsets, period_order):
        selected = [r for r in rows if str(r["period"]) == key]
        if not selected:
            continue
        selected.sort(key=lambda r: (0 if str(r["detector"]) == "FT" else 1, int(r["sector"])))
        label = next(p.label for p in PERIODS if p.key == key)
        ax.errorbar(x + offset, [float(r["scale_factor"]) for r in selected],
                    yerr=[float(r["scale_factor_err"]) for r in selected], fmt="o", label=label)
    ax.axhline(1.0, linestyle="--", linewidth=1.0)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_xlabel("Photon detector category")
    ax.set_ylabel(r"$S_\gamma=\epsilon_{data}/\epsilon_{MC}$")
    ax.set_title("Integrated photon-efficiency scale factors by run period")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def parse_matching_scan_key(key: str) -> Tuple[float, float, float]:
    try:
        angle = float(key.split("angle_", 1)[1].split("_energy_", 1)[0])
        energy = float(key.split("_energy_", 1)[1].split("_m2_", 1)[0])
        m2 = float(key.rsplit("_m2_", 1)[1])
        return angle, energy, m2
    except Exception as exc:
        raise ValueError(f"Malformed matching-scan key: {key}") from exc


def matching_scan_category_value(entry: object, category: str) -> float:
    if isinstance(entry, Mapping):
        return float(entry.get(category, 0.0))
    # Backward compatibility with v10/v11 JSON, where each key held only total counts.
    return float(entry) if category == "all" else 0.0


def matching_scan_summary(diag: Mapping[str, object], nominal_energy: float = 0.35,
                          nominal_m2: float = 0.20) -> Dict[str, object]:
    categories = ["all", "FT"] + [f"FD sector {i}" for i in range(1, 7)]
    output: Dict[str, object] = {
        "nominal_relative_energy_cut": nominal_energy,
        "nominal_predicted_m2_cut_GeV2": nominal_m2,
        "angle_scan": {},
    }
    for key, entry in diag.get("matching_scans", {}).items():
        angle, energy, m2 = parse_matching_scan_key(str(key))
        if not (math.isclose(energy, nominal_energy, rel_tol=0.0, abs_tol=1.0e-9)
                and math.isclose(m2, nominal_m2, rel_tol=0.0, abs_tol=1.0e-9)):
            continue
        output["angle_scan"][f"{angle:g}"] = {
            category: matching_scan_category_value(entry, category) for category in categories
        }
    return output


def plot_matching_scan_summary(path: Path, period_label: str,
                               data_diag: Mapping[str, object],
                               mc_diag: Mapping[str, object]) -> None:
    """Plot data/AAOGEN pass-count ratios for all detector categories."""
    categories = ["all", "FT"] + [f"FD sector {i}" for i in range(1, 7)]
    labels = {"all": "All", "FT": "FT", **{f"FD sector {i}": f"FD S{i}" for i in range(1, 7)}}
    dscan = data_diag.get("matching_scans", {})
    mscan = mc_diag.get("matching_scans", {})
    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    drew = False
    for category in categories:
        points = []
        for key, dentry in dscan.items():
            angle, energy, m2 = parse_matching_scan_key(str(key))
            if not (math.isclose(energy, 0.35, abs_tol=1.0e-9)
                    and math.isclose(m2, 0.20, abs_tol=1.0e-9)):
                continue
            dcount = matching_scan_category_value(dentry, category)
            mcount = matching_scan_category_value(mscan.get(key, {}), category)
            if mcount > 0.0:
                points.append((angle, dcount / mcount))
        if points:
            points.sort()
            ax.plot([p[0] for p in points], [p[1] for p in points], "o-", label=labels[category])
            drew = True
    if not drew:
        plt.close(fig)
        return
    ax.set_xlabel("Predicted-observed angular match threshold (deg)")
    ax.set_ylabel("Passing-probe weighted-count ratio, data/AAOGEN")
    ax.set_title(f"{period_label}: matching-cut stability (energy residual < 0.35, |m²| < 0.20 GeV²)")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def plot_mirror_diagnostics(path: Path, period_label: str,
                            data_diag: Mapping[str, object],
                            mc_diag: Mapping[str, object]) -> None:
    """Show mirror multiplicity and category migration for data and AAOGEN."""
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    labels = ["One solution", "Two solutions", "Same categories", "Different categories"]
    x = np.arange(len(labels), dtype=float)
    width = 0.38
    for offset, diag, sample in ((-width/2, data_diag, "Data"), (width/2, mc_diag, "AAOGEN MC")):
        md = diag.get("mirror_diagnostics", {})
        values = [float(md.get("events_with_one_solution", 0.0)),
                  float(md.get("events_with_two_solutions", 0.0)),
                  float(md.get("events_with_same_category_pair", 0.0)),
                  float(md.get("events_with_different_category_pair", 0.0))]
        axes[0, 0].bar(x + offset, values, width, label=sample)
    axes[0, 0].set_yscale("log")
    axes[0, 0].set_xticks(x); axes[0, 0].set_xticklabels(labels, rotation=25, ha="right")
    axes[0, 0].set_ylabel("Events")
    axes[0, 0].set_title("Mirror-solution event counts")
    axes[0, 0].legend(frameon=False)

    category_labels = ["FT"] + [f"S{i}" for i in range(1, 7)]
    panels = [
        (axes[0, 1], data_diag, "photon1_category_migration", "Data photon 1"),
        (axes[0, 2], data_diag, "photon2_category_migration", "Data photon 2"),
        (axes[1, 1], mc_diag, "photon1_category_migration", "AAOGEN photon 1"),
        (axes[1, 2], mc_diag, "photon2_category_migration", "AAOGEN photon 2"),
    ]
    for ax, diag, key, title in panels:
        matrix = np.asarray(diag.get("mirror_diagnostics", {}).get(key, np.zeros((7, 7))), dtype=float)
        row_sum = matrix.sum(axis=1, keepdims=True)
        normalized = np.divide(matrix, row_sum, out=np.zeros_like(matrix), where=row_sum > 0)
        image = ax.imshow(normalized, origin="lower", vmin=0.0, vmax=1.0, aspect="auto")
        ax.set_xticks(range(7)); ax.set_xticklabels(category_labels)
        ax.set_yticks(range(7)); ax.set_yticklabels(category_labels)
        ax.set_xlabel("Alternate mirror category"); ax.set_ylabel("Best mirror category")
        ax.set_title(title)
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Row-normalized fraction")
    axes[1, 0].axis("off")
    for row, diag, sample in ((0, data_diag, "Data"), (1, mc_diag, "AAOGEN MC")):
        md = diag.get("mirror_diagnostics", {})
        two = float(md.get("events_with_two_solutions", 0.0))
        different = float(md.get("events_with_different_category_pair", 0.0))
        text = (f"{sample}\nTwo-solution events: {two:,.0f}\n"
                f"Different category pair: {different:,.0f}\n"
                f"Fraction: {different/two:.4f}" if two > 0 else f"{sample}\nNo two-solution events")
        target = axes[1, 0] if row == 0 else axes[1, 0]
        target.text(0.05, 0.75 - 0.42*row, text, transform=target.transAxes, va="top", fontsize=11)
    fig.suptitle(f"{period_label}: cone-mirror category ambiguity")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def efficiency_from_counts_simple(npass: float, nfail: float) -> float:
    denom = npass + nfail
    return npass / denom if denom > 0.0 else math.nan


def plot_matching_scale_factor_scan(path: Path, period_label: str,
                                    data_diag: Mapping[str, object],
                                    mc_diag: Mapping[str, object],
                                    rows: Sequence[EfficiencyRow]) -> None:
    """Plot S_gamma versus angular matching threshold using nominal fail yields."""
    row_map = {("FT" if r.detector == "FT" else f"FD sector {r.sector}"): r for r in rows}
    categories = ["FT"] + [f"FD sector {i}" for i in range(1, 7)]
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.5), squeeze=False, sharex=True)
    for ax, category in zip(axes.flat, categories):
        row = row_map.get(category)
        if row is None:
            ax.axis("off"); continue
        points = []
        for key, dentry in data_diag.get("matching_scans", {}).items():
            angle, energy, m2 = parse_matching_scan_key(str(key))
            if not (math.isclose(energy, 0.35, abs_tol=1.0e-9)
                    and math.isclose(m2, 0.20, abs_tol=1.0e-9)):
                continue
            pd = matching_scan_category_value(dentry, category)
            pm = matching_scan_category_value(mc_diag.get("matching_scans", {}).get(key, {}), category)
            ed = efficiency_from_counts_simple(pd, row.fail_data)
            em = efficiency_from_counts_simple(pm, row.fail_mc)
            sf = ed / em if em > 0.0 else math.nan
            points.append((angle, sf))
        points.sort()
        ax.plot([p[0] for p in points], [p[1] for p in points], "o-")
        ax.axhline(1.0, linestyle="--", linewidth=0.8)
        ax.axvline(2.0, linestyle=":", linewidth=0.8)
        ax.set_title("FT" if category == "FT" else category.replace("sector ", "S"))
        ax.set_ylabel(r"$S_\gamma$")
        ax.set_xlabel("Angle threshold (deg)")
    axes[1, 3].axis("off")
    fig.suptitle(f"{period_label}: scale-factor stability versus angular match\n"
                 "Nominal fitted fail yields held fixed; energy residual < 0.35, |m²| < 0.20 GeV²")
    fig.tight_layout(); fig.savefig(path, dpi=180); plt.close(fig)


def plot_matching_stability_maps(path: Path, period_label: str,
                                 data_diag: Mapping[str, object],
                                 mc_diag: Mapping[str, object],
                                 rows: Sequence[EfficiencyRow]) -> None:
    """Plot S_gamma maps versus angle and predicted-mass threshold at nominal energy cut."""
    row_map = {("FT" if r.detector == "FT" else f"FD sector {r.sector}"): r for r in rows}
    categories = ["FT"] + [f"FD sector {i}" for i in range(1, 7)]
    all_keys = list(data_diag.get("matching_scans", {}).keys())
    angles = sorted({parse_matching_scan_key(str(k))[0] for k in all_keys
                     if math.isclose(parse_matching_scan_key(str(k))[1], 0.35, abs_tol=1.0e-9)})
    m2_values = sorted({parse_matching_scan_key(str(k))[2] for k in all_keys
                        if math.isclose(parse_matching_scan_key(str(k))[1], 0.35, abs_tol=1.0e-9)})
    fig, axes = plt.subplots(2, 4, figsize=(18, 9), squeeze=False)
    images = []
    for ax, category in zip(axes.flat, categories):
        row = row_map.get(category)
        matrix = np.full((len(m2_values), len(angles)), np.nan, dtype=float)
        if row is not None:
            for iy, m2 in enumerate(m2_values):
                for ix, angle in enumerate(angles):
                    key = f"angle_{angle:g}_energy_0.35_m2_{m2:g}"
                    pd = matching_scan_category_value(data_diag.get("matching_scans", {}).get(key, {}), category)
                    pm = matching_scan_category_value(mc_diag.get("matching_scans", {}).get(key, {}), category)
                    ed = efficiency_from_counts_simple(pd, row.fail_data)
                    em = efficiency_from_counts_simple(pm, row.fail_mc)
                    matrix[iy, ix] = ed / em if em > 0.0 else math.nan
        if matrix.size == 0 or len(angles) == 0 or len(m2_values) == 0:
            ax.text(0.5, 0.5, "Matching scan disabled", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue
        image = ax.imshow(matrix, origin="lower", aspect="auto", vmin=0.5, vmax=1.2)
        images.append(image)
        ax.set_xticks(range(len(angles))); ax.set_xticklabels([f"{v:g}" for v in angles])
        ax.set_yticks(range(len(m2_values))); ax.set_yticklabels(["none" if v < 0 else f"{v:g}" for v in m2_values])
        ax.set_xlabel("Angle threshold (deg)"); ax.set_ylabel(r"$|m^2_{pred}|$ threshold (GeV$^2$)")
        ax.set_title("FT" if category == "FT" else category.replace("sector ", "S"))
    axes[1, 3].axis("off")
    if images:
        fig.colorbar(images[0], ax=[a for a in axes.flat if a.get_visible()], fraction=0.025, pad=0.02,
                     label=r"$S_\gamma$ (nominal fail yields fixed)")
    fig.suptitle(f"{period_label}: matching-cut scale-factor stability maps (energy residual < 0.35)")
    fig.subplots_adjust(top=0.90, right=0.92, wspace=0.28, hspace=0.30)
    fig.savefig(path, dpi=180); plt.close(fig)

def _fraction_sequence(cutflow: Mapping[str, float], stages: Sequence[str]) -> np.ndarray:
    first = max(float(cutflow.get(stages[0], 0.0)), 1.0)
    return np.asarray([float(cutflow.get(stage, 0.0))/first for stage in stages], dtype=float)


def plot_pass_cutflow(path: Path, period_label: str, data_diag: Mapping[str, object], mc_diag: Mapping[str, object]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(13, 10))
    for ax, group, stages, title in ((axes[0], "event_cutflow", PASS_EVENT_STAGES, "Native eppi0 event cut flow"),
                                     (axes[1], "probe_cutflow", PASS_PROBE_STAGES, "Directed tag-and-probe cut flow")):
        x=np.arange(len(stages)); width=0.38
        for offset, diag, label in ((-width/2, data_diag, "Data"),(width/2, mc_diag, "AAOGEN MC")):
            vals=_fraction_sequence(diag[group], stages)
            ax.bar(x+offset, vals, width, label=label)
            for xx,v,stage in zip(x+offset, vals, stages):
                ax.text(xx, max(v,1e-6)*1.08, f"{float(diag[group][stage]):,.1f}", rotation=90, ha="center", va="bottom", fontsize=7)
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
    plot_matching_scan_summary(audit_dir/"matching_cut_scan.png", period_label, data_diag, mc_diag)
    plot_mirror_diagnostics(audit_dir/"mirror_ambiguity.png", period_label, data_diag, mc_diag)


def fit_variant_model_systematic(data_hists, dvcs_hists, pi0_hists, edges_by_var,
                                 morphology, total_candidates, min_counts):
    """Return nominal fit and a conservative yield spread from credible projection variants."""
    nominal = fit_multi_projection_bin(data_hists, dvcs_hists, pi0_hists, edges_by_var,
                                       morphology, total_candidates, min_counts)
    variants = {"nominal": nominal}
    for omitted in FRACTION_DRIVER_KEYS:
        saved = []
        for spec in FIT_VARIABLES:
            if spec.key == omitted:
                saved.append((spec, data_hists[spec.key].copy(), dvcs_hists[spec.key].copy(), pi0_hists[spec.key].copy()))
                data_hists[spec.key] = np.zeros_like(data_hists[spec.key])
                dvcs_hists[spec.key] = np.zeros_like(dvcs_hists[spec.key])
                pi0_hists[spec.key] = np.zeros_like(pi0_hists[spec.key])
        variants[f"omit_{omitted}"] = fit_multi_projection_bin(data_hists, dvcs_hists, pi0_hists,
                                                                 edges_by_var, morphology,
                                                                 total_candidates, min_counts)
        for spec, da, dv, pi in saved:
            data_hists[spec.key], dvcs_hists[spec.key], pi0_hists[spec.key] = da, dv, pi
    good = [v.n_pi0 for v in variants.values() if v.success and math.isfinite(v.n_pi0)]
    spread = max((abs(v - nominal.n_pi0) for v in good), default=0.0) if nominal.success else math.nan
    return nominal, spread, variants


def plot_combined_fd_fit(path: Path, period_label: str, fit: MultiFitSummary) -> None:
    definition = BinDefinition(-1, "FD", 0, 0.0, 0.0, 0.0, 0.0)
    plot_category_fit_summary(path, period_label + " (all six FD sectors combined)", definition, fit)




_EXCLUSIVITY_FITTER_MODULE = None


def load_exclusivity_fitter_module():
    """Load the established exclusivity two-template fitter from nearby paths.

    During standalone development this script may live under ``external_scripts``
    while the validated exclusivity program remains one directory higher.  Search
    a small, deterministic set of nearby locations rather than requiring the two
    files to be colocated.
    """
    global _EXCLUSIVITY_FITTER_MODULE
    if _EXCLUSIVITY_FITTER_MODULE is not None:
        return _EXCLUSIVITY_FITTER_MODULE

    script_dir = Path(__file__).resolve().parent
    cwd = Path.cwd().resolve()
    search_dirs = [
        script_dir,
        script_dir.parent,
        script_dir.parent / "external_scripts",
        cwd,
        cwd.parent,
        cwd / "external_scripts",
        cwd.parent / "external_scripts",
    ]

    # Preserve order while removing duplicate resolved directories.
    unique_dirs = []
    seen_dirs = set()
    for directory in search_dirs:
        resolved = directory.resolve()
        if resolved not in seen_dirs:
            seen_dirs.add(resolved)
            unique_dirs.append(resolved)

    candidates = []
    canonical_name = "plot_exclusivity_data_dvcs_pi0_mc.py"
    versioned_pattern = "plot_exclusivity_data_dvcs_pi0_mc(*).py"
    for directory in unique_dirs:
        # Prefer the canonical repository filename in every search directory.
        candidates.append(directory / canonical_name)
    for directory in unique_dirs:
        # Fall back to versioned development copies, newest filename first.
        candidates.extend(sorted(directory.glob(versioned_pattern), reverse=True))

    source = next((candidate for candidate in candidates if candidate.is_file()), None)
    if source is None:
        searched = "\n  - ".join(str(directory) for directory in unique_dirs)
        raise FileNotFoundError(
            "Could not locate the validated exclusivity fitter "
            f"{canonical_name}. Searched:\n  - {searched}"
        )

    spec = importlib.util.spec_from_file_location("_photon_eff_exclusivity_fitter", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load exclusivity fitter from {source}")
    module = importlib.util.module_from_spec(spec)
    import sys
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)

    # The validated fitter normally uses the production-DVCS histogram
    # definitions stored in its module-level VARIABLES tuple.  This efficiency
    # study intentionally uses broader supports, so the fitter's VariableConfig
    # objects must be replaced with definitions that exactly match the
    # histograms supplied below.  Merely changing FIT_VARIABLES in this script
    # is insufficient: the fitter builds its morphing matrices from its own
    # VariableConfig.bins/xmin/xmax values.  A mismatch (for example a 120-bin
    # Delta_phi histogram with the fitter's native 100-bin Delta_phi config)
    # produces a 100x100 transport matrix multiplied by a 120-element template.
    fitter_variables = []
    aliases_by_branch = {
        variable.branch: tuple(getattr(variable, "aliases", ()))
        for variable in getattr(module, "VARIABLES", ())
    }
    for fit_spec in FIT_VARIABLES:
        branch = "theta" if fit_spec.key == "theta_cm" else fit_spec.key
        fitter_variables.append(
            module.VariableConfig(
                branch,
                fit_spec.label,
                int(fit_spec.bins),
                float(fit_spec.low),
                float(fit_spec.high),
                aliases=aliases_by_branch.get(branch, ()),
            )
        )
    module.VARIABLES = tuple(fitter_variables)

    # Validate the adapter immediately so future changes fail with a readable
    # message instead of a deep scipy/numpy matrix-dimension traceback.
    expected = {
        ("theta" if item.key == "theta_cm" else item.key): int(item.bins)
        for item in FIT_VARIABLES
    }
    configured = {item.branch: int(item.bins) for item in module.VARIABLES}
    if configured != expected:
        raise RuntimeError(
            "Internal exclusivity-fitter histogram configuration mismatch: "
            f"expected {expected}, configured {configured}"
        )

    _EXCLUSIVITY_FITTER_MODULE = module
    return module


def _exact_histogram_mapping(histograms: Mapping[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Translate this script's branch aliases to the exclusivity fitter keys."""
    return {
        "Delta_phi": np.asarray(histograms["Delta_phi"], dtype=np.float64),
        "theta": np.asarray(histograms["theta_cm"], dtype=np.float64),
        "theta_gamma_gamma": np.asarray(histograms["theta_gamma_gamma"], dtype=np.float64),
        "pTmiss": np.asarray(histograms["pTmiss"], dtype=np.float64),
        "Emiss2": np.asarray(histograms["Emiss2"], dtype=np.float64),
        "Mx2": np.asarray(histograms["Mx2"], dtype=np.float64),
        "Mx2_2": np.asarray(histograms["Mx2_2"], dtype=np.float64),
    }


def _current_key_from_exact(key: str) -> str:
    return "theta_cm" if key == "theta" else key


def convert_exact_summary(summary, data_hists: Mapping[str, np.ndarray]) -> MultiFitSummary:
    """Convert the established fitter result to this script's output schema."""
    projections: Dict[str, ProjectionResult] = {}
    if summary.variable_results:
        for exact_key, result in summary.variable_results.items():
            current_key = _current_key_from_exact(exact_key)
            spec = next((item for item in FIT_VARIABLES if item.key == current_key), None)
            if spec is None or not result.success or result.model_counts is None:
                continue
            edges = np.linspace(spec.low, spec.high, spec.bins + 1)
            projections[current_key] = ProjectionResult(
                spec=spec,
                edges=edges,
                data_counts=np.asarray(result.fit_data_counts, dtype=np.float64),
                model_counts=np.asarray(result.model_counts, dtype=np.float64),
                pi0_counts=np.asarray(result.pi0_component_counts, dtype=np.float64),
                dvcs_counts=np.asarray(result.dvcs_component_counts, dtype=np.float64),
                deviance=float(result.deviance),
                ndf=int(result.ndf),
            )
    total_candidates = float(np.sum(np.asarray(data_hists["Delta_phi"], dtype=np.float64)))
    fraction = float(summary.f_pi0)
    fraction_err = float(summary.f_pi0_err) if math.isfinite(summary.f_pi0_err) else math.nan
    return MultiFitSummary(
        success=bool(summary.success),
        n_pi0=fraction * total_candidates,
        n_pi0_err=(fraction_err * total_candidates if math.isfinite(fraction_err) else math.sqrt(max(fraction * total_candidates, 0.0))),
        n_dvcs=(1.0 - fraction) * total_candidates,
        n_dvcs_err=(fraction_err * total_candidates if math.isfinite(fraction_err) else math.sqrt(max((1.0-fraction) * total_candidates, 0.0))),
        fraction_pi0=fraction,
        fraction_pi0_err=fraction_err,
        deviance=float(summary.deviance),
        ndf=int(summary.ndf),
        message=str(summary.message),
        projections=projections,
    )


def run_exact_exclusivity_fit(
    data_hists: Mapping[str, np.ndarray],
    dvcs_hists: Mapping[str, np.ndarray],
    pi0_hists: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    fraction_drivers: Optional[Sequence[str]] = None,
    topology_key: str = "CD_FT",
) -> Tuple[MultiFitSummary, object]:
    """Run the exact support/profile fitter used by the exclusivity analysis."""
    module = load_exclusivity_fitter_module()
    topology_map = {
        "CD_FT": (2, 0),
        "CD_FD": (2, 1),
        "FD_FD": (1, 1),
    }
    if topology_key not in topology_map:
        raise ValueError(f"Unsupported production topology: {topology_key}")
    detector1, detector2 = topology_map[topology_key]
    topology = module.TopologyConfig(
        topology_key,
        TOPOLOGY_LABELS.get(topology_key, topology_key),
        detector1,
        detector2,
    )
    active_drivers = tuple(FRACTION_DRIVER_KEYS if fraction_drivers is None else fraction_drivers)
    exact_drivers = tuple("theta" if key == "theta_cm" else key for key in active_drivers)

    mapped_data = _exact_histogram_mapping(data_hists)
    mapped_dvcs = _exact_histogram_mapping(dvcs_hists)
    mapped_pi0 = _exact_histogram_mapping(pi0_hists)
    configured_bins = {item.branch: int(item.bins) for item in module.VARIABLES}
    for sample_name, mapping in (
        ("data", mapped_data),
        ("DVCS MC", mapped_dvcs),
        ("pi0 MC", mapped_pi0),
    ):
        for branch, expected_bins in configured_bins.items():
            actual_bins = int(np.asarray(mapping[branch]).size)
            if actual_bins != expected_bins:
                raise RuntimeError(
                    "Histogram/fitter bin mismatch before exclusivity fit: "
                    f"{sample_name} {branch} has {actual_bins} bins, but the "
                    f"configured fitter expects {expected_bins}."
                )

    summary = module.fit_shared_two_templates(
        mapped_data,
        mapped_dvcs,
        mapped_pi0,
        topology,
        max_shift_bins=args.exclusivity_max_shift_bins,
        max_smear_bins=args.exclusivity_max_smear_bins,
        min_counts=args.fit_min_counts,
        fraction_variable_branches=exact_drivers,
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
    return convert_exact_summary(summary, data_hists), summary


def fit_frozen_exclusivity_shapes(
    data_hists: Mapping[str, np.ndarray],
    reference_exact_summary,
    drivers: Sequence[str],
) -> MultiFitSummary:
    """Fit only f_pi0 using morphology frozen by the high-statistics reference fit.

    This is the simultaneous-sector implementation: the combined-FD fit learns
    the nuisance morphology once; each sector then profiles only its mixture
    fraction. It is much faster and prevents random sector nuisance solutions.
    """
    module = load_exclusivity_fitter_module()
    results = reference_exact_summary.variable_results or {}
    usable = []
    for key in drivers:
        exact_key = "theta" if key == "theta_cm" else key
        result = results.get(exact_key)
        if result is None or not result.success or result.transformed_dvcs_shape is None or result.transformed_pi0_shape is None:
            continue
        current_key = _current_key_from_exact(exact_key)
        data = np.asarray(data_hists[current_key], dtype=np.float64)
        usable.append((current_key, result, data))
    if not usable:
        return MultiFitSummary(False, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, 0, "No usable frozen projections", {})

    def objective(fraction: float) -> float:
        value = 0.0
        for _, result, data in usable:
            dv = np.asarray(result.transformed_dvcs_shape, dtype=np.float64)
            pi = np.asarray(result.transformed_pi0_shape, dtype=np.float64)
            shape = (1.0 - fraction) * dv + fraction * pi
            total = float(np.sum(data))
            if total <= 0.0 or float(np.sum(shape)) <= 0.0:
                continue
            model = total * shape / float(np.sum(shape))
            value += module.poisson_deviance(data, model)
        return value

    opt = minimize_scalar(objective, bounds=(0.0, 1.0), method="bounded", options={"xatol": 2e-5, "maxiter": 150})
    fraction = float(opt.x) if opt.success else math.nan
    if not math.isfinite(fraction):
        return MultiFitSummary(False, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, 0, "Frozen-shape fraction fit failed", {})
    step = 1e-3
    ferr = math.nan
    if step < fraction < 1.0-step:
        curvature = (objective(fraction-step) - 2*objective(fraction) + objective(fraction+step)) / step**2
        if curvature > 0.0:
            ferr = math.sqrt(2.0/curvature)
    projections: Dict[str, ProjectionResult] = {}
    total_deviance = 0.0
    total_ndf = 0
    for spec in FIT_VARIABLES:
        exact_key = "theta" if spec.key == "theta_cm" else spec.key
        result = results.get(exact_key)
        if result is None or not result.success or result.transformed_dvcs_shape is None:
            continue
        data = np.asarray(data_hists[spec.key], dtype=np.float64)
        dv = np.asarray(result.transformed_dvcs_shape, dtype=np.float64)
        pi = np.asarray(result.transformed_pi0_shape, dtype=np.float64)
        total = float(np.sum(data))
        shape = (1.0-fraction)*dv + fraction*pi
        norm = float(np.sum(shape))
        if total <= 0.0 or norm <= 0.0:
            continue
        scale = total / norm
        model = scale*shape
        dvcomp = scale*(1.0-fraction)*dv
        picomp = scale*fraction*pi
        dev = module.poisson_deviance(data, model)
        ndf = max(0, int(np.count_nonzero(data + model)) - 1)
        if spec.key in drivers:
            total_deviance += dev
            total_ndf += ndf
        projections[spec.key] = ProjectionResult(spec, np.linspace(spec.low, spec.high, spec.bins+1), data, model, picomp, dvcomp, dev, ndf)
    total_candidates = float(np.sum(np.asarray(data_hists["Delta_phi"], dtype=np.float64)))
    return MultiFitSummary(True, fraction*total_candidates, (ferr*total_candidates if math.isfinite(ferr) else math.sqrt(max(fraction*total_candidates,0.0))), (1-fraction)*total_candidates, (ferr*total_candidates if math.isfinite(ferr) else math.sqrt(max((1-fraction)*total_candidates,0.0))), fraction, ferr, total_deviance, total_ndf, "Exact exclusivity morphology frozen from combined category", projections)


def frozen_variant_systematic(data_hists, reference_exact_summary):
    nominal = fit_frozen_exclusivity_shapes(data_hists, reference_exact_summary, FRACTION_DRIVER_KEYS)
    variants = {"nominal": nominal}
    for omitted in FRACTION_DRIVER_KEYS:
        active = tuple(key for key in FRACTION_DRIVER_KEYS if key != omitted)
        variants[f"omit_{omitted}"] = fit_frozen_exclusivity_shapes(data_hists, reference_exact_summary, active)
    # The fit projections cover only their configured histogram supports.  The
    # leave-one-variable-out variation is therefore retained in *fraction*
    # units and converted to a yield only after multiplying by the full number
    # of selected failing candidates in the detector category.
    good = [fit.fraction_pi0 for fit in variants.values()
            if fit.success and math.isfinite(fit.fraction_pi0)]
    spread = max((abs(value - nominal.fraction_pi0) for value in good), default=0.0)         if nominal.success else math.nan
    return nominal, spread, variants

def summarize_selected_one_photon_support(label: str, trials: TrialArrays, args: argparse.Namespace) -> Dict[str, float]:
    """Summarize selected one-photon candidates without imposing a low-tag quota.

    Tags in 0.4 <= E_tag < E_probe,min are useful diagnostics, but their
    fraction is sample dependent. In particular, the DVCS-generator sample is
    only a background-template shape and is not expected to contain a large
    low-energy-tag population. The only fatal condition is an empty selected
    sample.
    """
    total = float(np.sum(trials.weight))
    low_mask = (
        np.isfinite(trials.tag_E)
        & (trials.tag_E >= args.tag_E_min)
        & (trials.tag_E < args.probe_E_min)
    )
    low = float(np.sum(trials.weight[low_mask]))
    fraction = low / total if total > 0.0 else math.nan
    log(
        f"{label}: selected one-photon total={total:,.1f}, "
        f"0.4-to-probe-threshold tags={low:,.1f}, fraction={fraction:.4f}"
    )
    if total <= 0.0:
        raise RuntimeError(
            f"{label}: no selected one-photon candidates remain after the common "
            "tag/probe selection"
        )
    return {
        "total_weight": total,
        "low_tag_weight": low,
        "low_tag_fraction": fraction,
        "low_tag_fraction_is_diagnostic_only": True,
    }


def category_fit_support_fraction(trials: TrialArrays, bin_ids: np.ndarray, bin_id: int) -> float:
    selected = bin_ids == bin_id
    denominator = float(np.sum(trials.weight[selected]))
    if denominator <= 0.0:
        return 0.0
    supported = selected & common_fit_mask(trials)
    return float(np.sum(trials.weight[supported])) / denominator



def plot_one_photon_cutflows(path: Path, period_label: str,
                             diagnostics: Sequence[Mapping[str, object]]) -> None:
    stages = ["all_tree_entries", "basic_quality", "production_event_cuts",
              "leading_photon_energy", "leading_photon_detector", "selected"]
    fig, axes = plt.subplots(len(diagnostics), 1, figsize=(13, 4.2 * len(diagnostics)), squeeze=False)
    for ax, diag in zip(axes[:, 0], diagnostics):
        values = np.asarray([float(diag["stages"].get(stage, 0.0)) for stage in stages])
        x = np.arange(len(stages))
        ax.plot(x, values, marker="o")
        ax.set_yscale("log")
        ax.set_xticks(x); ax.set_xticklabels(stages, rotation=35, ha="right")
        ax.set_ylabel("Candidates")
        ax.set_title(str(diag.get("label", "one-photon sample")))
        ax.grid(alpha=0.25)
        for xx, value in zip(x, values):
            if value > 0:
                ax.annotate(f"{value:,.0f}", (xx, value), xytext=(0, 5), textcoords="offset points",
                            ha="center", fontsize=8)
    fig.suptitle(f"{period_label}: one-photon candidate cutflow")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _shape_hist(trials: TrialArrays, attr: str, edges: np.ndarray,
                selection: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:
    values = np.asarray(getattr(trials, attr), dtype=float)
    weights = np.asarray(trials.weight, dtype=float)
    finite = np.isfinite(values)
    if selection is not None:
        finite &= selection
    total = float(np.sum(weights[finite]))
    inside = finite & (values >= edges[0]) & (values < edges[-1])
    counts = np.histogram(values[inside], bins=edges, weights=weights[inside])[0].astype(float)
    integral = float(np.sum(counts))
    return (counts / integral if integral > 0 else counts), (integral / total if total > 0 else math.nan)


def plot_shape_comparison_canvas(path: Path, period_label: str, category_label: str,
                                 data: TrialArrays, dvcs: TrialArrays, pi0: TrialArrays,
                                 category_masks: Tuple[np.ndarray, np.ndarray, np.ndarray],
                                 broad: bool = False) -> None:
    if broad:
        specs = [
            ("delta_phi", r"$\Delta\phi$ (rad)", 120, 0.0, 2*math.pi),
            ("theta_gamma_gamma", r"$\theta_{\gamma\gamma}$ (rad)", 120, 0.0, math.pi),
            ("pTmiss", r"$p_T^{\mathrm{miss}}$ (GeV)", 120, 0.0, 3.0),
            ("theta_cm", r"$\theta_{p\gamma}^{\mathrm{CM}}$ (rad)", 100, 0.0, math.pi),
            ("Emiss2", r"$E_{\mathrm{miss}}$ (GeV)", 120, -5.0, 5.0),
            ("Mx2", r"$M_x^2$ (GeV$^2$)", 120, -2.0, 2.0),
            ("Mx2_2", r"$M_{x2}^2$ (GeV$^2$)", 120, -2.0, 6.0),
        ]
    else:
        specs = [(sp.attr, sp.label, sp.bins, sp.low, sp.high) for sp in FIT_VARIABLES]
    fig, axes = plt.subplots(2, 4, figsize=(19, 8.5), squeeze=False)
    axes = axes.ravel()
    samples = [("Data one-photon", data, category_masks[0]),
               ("DVCS/BH template", dvcs, category_masks[1]),
               ("pi0 template", pi0, category_masks[2])]
    for ax, (attr, label, bins, low, high) in zip(axes, specs):
        edges = np.linspace(low, high, bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        for sample_label, trials, mask in samples:
            shape, containment = _shape_hist(trials, attr, edges, mask)
            ax.step(centers, shape, where="mid", linewidth=1.3,
                    label=f"{sample_label} ({100*containment:.1f}% in range)" if math.isfinite(containment) else sample_label)
        ax.set_xlabel(label); ax.set_ylabel("Unit-normalized entries")
        ax.grid(alpha=0.2); ax.legend(fontsize=7)
    for ax in axes[len(specs):]:
        ax.axis("off")
    kind = "broad diagnostic ranges" if broad else "validated fit ranges"
    fig.suptitle(f"{period_label}: {category_label} one-photon shape comparison — {kind}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_all_period_efficiency_summary(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    period_order = [p.key for p in PERIODS]
    labels = ["FT"] + [f"FD S{i}" for i in range(1, 7)]
    x = np.arange(7, dtype=float)
    offsets = np.linspace(-0.18, 0.18, len(period_order))
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    for offset, key in zip(offsets, period_order):
        selected = [r for r in rows if str(r["period"]) == key]
        if not selected:
            continue
        selected.sort(key=lambda r: (0 if str(r["detector"]) == "FT" else 1, int(r["sector"])))
        label = next(p.label for p in PERIODS if p.key == key)
        axes[0].errorbar(x + offset, [float(r["efficiency_data"]) for r in selected],
                         yerr=[float(r["efficiency_data_err"]) for r in selected], fmt="o", label=label)
        axes[1].errorbar(x + offset, [float(r["efficiency_mc"]) for r in selected],
                         yerr=[float(r["efficiency_mc_err"]) for r in selected], fmt="o", label=label)
    axes[0].set_ylabel(r"$\epsilon_{data}$"); axes[1].set_ylabel(r"$\epsilon_{MC}$")
    for ax in axes:
        ax.set_ylim(0.0, 1.05); ax.grid(alpha=0.25); ax.legend(frameon=False, ncol=3)
    axes[1].set_xticks(x); axes[1].set_xticklabels(labels); axes[1].set_xlabel("Photon detector category")
    fig.suptitle("Integrated photon efficiencies by run period")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=180)
    plt.close(fig)



TOPOLOGY_LABELS: Mapping[str, str] = {
    "CD_FT": "(CD, FT)",
    "CD_FD": "(CD, FD)",
    "FD_FD": "(FD, FD)",
}


def trial_topology_mask(trials: TrialArrays, topology_key: str) -> np.ndarray:
    """Return the production-analysis proton/photon topology mask."""
    if topology_key == "CD_FT":
        return (trials.proton_detector == 2) & (trials.detector == 0)
    if topology_key == "CD_FD":
        return (trials.proton_detector == 2) & (trials.detector == 1)
    if topology_key == "FD_FD":
        return (trials.proton_detector == 1) & (trials.detector == 1)
    raise ValueError(f"Unsupported topology key: {topology_key}")


def topology_keys_for_definition(definition: BinDefinition) -> Tuple[str, ...]:
    return ("CD_FT",) if definition.detector == "FT" else ("CD_FD", "FD_FD")


def histograms_for_mask(
    trials: TrialArrays,
    mask: np.ndarray,
    edges_by_var: Mapping[str, np.ndarray],
) -> Dict[str, np.ndarray]:
    """Fill all template projections from one explicitly common event mask."""
    result: Dict[str, np.ndarray] = {}
    selected_weights = trials.weight[mask]
    for spec in FIT_VARIABLES:
        values = np.asarray(getattr(trials, spec.attr), dtype=float)[mask]
        finite = np.isfinite(values) & np.isfinite(selected_weights)
        result[spec.key] = np.histogram(
            values[finite],
            bins=edges_by_var[spec.key],
            weights=selected_weights[finite],
        )[0].astype(float)
    return result


def weighted_count_for_mask(trials: TrialArrays, mask: np.ndarray) -> Tuple[float, float]:
    weights = np.asarray(trials.weight[mask], dtype=float)
    return float(np.sum(weights)), float(np.sum(weights * weights))


def exact_fit_variants(
    data_hists: Mapping[str, np.ndarray],
    dvcs_hists: Mapping[str, np.ndarray],
    pi0_hists: Mapping[str, np.ndarray],
    args: argparse.Namespace,
    topology_key: str,
) -> Tuple[MultiFitSummary, float, Dict[str, MultiFitSummary], object]:
    """Fit Delta_phi, theta_gamma_gamma and pTmiss and assess driver stability."""
    nominal, exact = run_exact_exclusivity_fit(
        data_hists, dvcs_hists, pi0_hists, args, FRACTION_DRIVER_KEYS,
        topology_key=topology_key,
    )
    variants: Dict[str, MultiFitSummary] = {"nominal": nominal}
    for omitted in FRACTION_DRIVER_KEYS:
        active = tuple(key for key in FRACTION_DRIVER_KEYS if key != omitted)
        variant, _ = run_exact_exclusivity_fit(
            data_hists, dvcs_hists, pi0_hists, args, active,
            topology_key=topology_key,
        )
        variants[f"omit_{omitted}"] = variant
    good = [
        fit.fraction_pi0 for fit in variants.values()
        if fit.success and math.isfinite(fit.fraction_pi0)
    ]
    spread = (
        max((abs(value - nominal.fraction_pi0) for value in good), default=0.0)
        if nominal.success else math.nan
    )
    if args.disable_fit_variant_systematic:
        spread = 0.0
    return nominal, spread, variants, exact


def process_period(period: PeriodConfig, args_dict: Mapping[str, object]) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    period_start = time.perf_counter()
    args = argparse.Namespace(**args_dict)
    args.period_key = period.key
    period_dir = Path(args.output_dir) / period.key
    aggregate_fit_dir = period_dir / "epgamma_template_fits"
    plot_dir = period_dir / "plots"
    aggregate_fit_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    pass_data, pass_data_diag = read_pass_trials(period.epgg_data, period.beam_energy_GeV, args, "data")
    pass_mc, pass_mc_diag = read_pass_trials(period.epgg_mc, period.beam_energy_GeV, args, "mc")
    write_pass_audit(period_dir, period.label, pass_data_diag, pass_mc_diag)
    if args.skip_efficiency_extraction:
        metadata = {
            "period": asdict(period),
            "pass_sample_audit": {"data": pass_data_diag, "mc": pass_mc_diag},
            "mirror_policy": args.mirror_policy,
            "matching_scan_summary": {
                "data": matching_scan_summary(pass_data_diag),
                "mc": matching_scan_summary(pass_mc_diag),
            },
            "runtime_seconds": time.perf_counter() - period_start,
        }
        with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2)
        return period.key, [], metadata

    data_one_photon_candidates, data_one_photon_diag = read_one_photon_candidates(
        period.epg_data, period.beam_energy_GeV, args, "data one-photon candidates"
    )
    dvcs_background_template, dvcs_template_diag = read_one_photon_candidates(
        period.dvcs_mc, period.beam_energy_GeV, args, "DVCS/BH background template"
    )
    pi0_missing_probe_mc, pi0_missing_diag = read_one_photon_candidates(
        period.pi0_as_epg_mc, period.beam_energy_GeV, args, "pi0 missing-probe MC/template"
    )

    plot_one_photon_cutflows(
        period_dir / "one_photon_candidate_cutflow.png", period.label,
        [data_one_photon_diag, dvcs_template_diag, pi0_missing_diag],
    )

    one_photon_support_audit = {
        "data_one_photon_candidates": summarize_selected_one_photon_support(
            f"{period.label} data one-photon candidates", data_one_photon_candidates, args
        ),
        "dvcs_background_template": summarize_selected_one_photon_support(
            f"{period.label} DVCS background template", dvcs_background_template, args
        ),
        "pi0_missing_probe_mc": summarize_selected_one_photon_support(
            f"{period.label} pi0 missing-probe MC", pi0_missing_probe_mc, args
        ),
    }

    log(
        f"{period.label}: pass data weight={float(np.sum(pass_data.weight)):,.1f} "
        f"({pass_data.size():,} stored), pass MC weight={float(np.sum(pass_mc.weight)):,.1f} "
        f"({pass_mc.size():,} stored), "
        f"data one-photon candidates={data_one_photon_candidates.size():,}, "
        f"DVCS background-template events={dvcs_background_template.size():,}, "
        f"pi0 missing-probe/template events={pi0_missing_probe_mc.size():,}"
    )

    energy_support = {
        "pass_data": summarize_tag_regions(pass_data),
        "pass_mc": summarize_tag_regions(pass_mc),
        "data_one_photon_candidates": summarize_tag_regions(data_one_photon_candidates),
        "dvcs_background_template": summarize_tag_regions(dvcs_background_template),
        "pi0_missing_probe_mc": summarize_tag_regions(pi0_missing_probe_mc),
    }
    plot_selected_energy_support(
        period_dir / "selected_energy_support.png", period.label,
        pass_data, pass_mc, data_one_photon_candidates, dvcs_background_template, pi0_missing_probe_mc, args,
    )


    bins = build_bins(pass_mc, pi0_missing_probe_mc, args)
    rows: List[EfficiencyRow] = []
    n_bins = len(bins)
    edges_by_var = fit_edges()

    log(f"{period.label}: assigning trials to {n_bins} integrated detector categories")
    pass_data_ids = assign_bin_ids(pass_data, bins)
    pass_mc_ids = assign_bin_ids(pass_mc, bins)
    data_one_photon_ids = assign_bin_ids(data_one_photon_candidates, bins)
    dvcs_template_ids = assign_bin_ids(dvcs_background_template, bins)
    pi0_missing_probe_ids = assign_bin_ids(pi0_missing_probe_mc, bins)

    shape_dir = period_dir / "shape_comparisons"
    shape_dir.mkdir(parents=True, exist_ok=True)

    # Reproduce the production-analysis topology separation.  In particular,
    # FT means (CD,FT); FD is the sum of distinct (CD,FD) and (FD,FD) fits.
    for topology_key in ("CD_FT", "CD_FD", "FD_FD"):
        data_top = trial_topology_mask(data_one_photon_candidates, topology_key)
        dvcs_top = trial_topology_mask(dvcs_background_template, topology_key)
        pi0_top = trial_topology_mask(pi0_missing_probe_mc, topology_key)
        if not (np.any(data_top) and np.any(dvcs_top) and np.any(pi0_top)):
            continue
        label = topology_key.lower()
        plot_shape_comparison_canvas(
            shape_dir / f"{label}_fit_ranges.png",
            period.label,
            TOPOLOGY_LABELS[topology_key],
            data_one_photon_candidates,
            dvcs_background_template,
            pi0_missing_probe_mc,
            (data_top, dvcs_top, pi0_top),
            broad=False,
        )
        plot_shape_comparison_canvas(
            shape_dir / f"{label}_broad_ranges.png",
            period.label,
            TOPOLOGY_LABELS[topology_key],
            data_one_photon_candidates,
            dvcs_background_template,
            pi0_missing_probe_mc,
            (data_top, dvcs_top, pi0_top),
            broad=True,
        )

    pass_data_counts = weighted_counts_by_bin(pass_data, pass_data_ids, n_bins)
    pass_mc_counts = weighted_counts_by_bin(pass_mc, pass_mc_ids, n_bins)
    pass_data_sumw2 = weighted_sumw2_by_bin(pass_data, pass_data_ids, n_bins)
    pass_mc_sumw2 = weighted_sumw2_by_bin(pass_mc, pass_mc_ids, n_bins)

    # Learn the nuisance morphology separately for each production topology.
    topology_references: Dict[str, object] = {}
    topology_reference_fits: Dict[str, MultiFitSummary] = {}
    topology_reference_variants: Dict[str, object] = {}
    for topology_key in ("CD_FT", "CD_FD", "FD_FD"):
        if topology_key == "CD_FT":
            category_mask_data = data_one_photon_candidates.detector == 0
            category_mask_dvcs = dvcs_background_template.detector == 0
            category_mask_pi0 = pi0_missing_probe_mc.detector == 0
        else:
            category_mask_data = data_one_photon_candidates.detector == 1
            category_mask_dvcs = dvcs_background_template.detector == 1
            category_mask_pi0 = pi0_missing_probe_mc.detector == 1

        data_mask = category_mask_data & trial_topology_mask(data_one_photon_candidates, topology_key)
        dvcs_mask = category_mask_dvcs & trial_topology_mask(dvcs_background_template, topology_key)
        pi0_mask = category_mask_pi0 & trial_topology_mask(pi0_missing_probe_mc, topology_key)

        data_h = histograms_for_mask(data_one_photon_candidates, data_mask, edges_by_var)
        dvcs_h = histograms_for_mask(dvcs_background_template, dvcs_mask, edges_by_var)
        pi0_h = histograms_for_mask(pi0_missing_probe_mc, pi0_mask, edges_by_var)

        if min(np.sum(data_h[key]) for key in FRACTION_DRIVER_KEYS) < args.fit_min_counts:
            continue
        if min(np.sum(dvcs_h[key]) for key in FRACTION_DRIVER_KEYS) < args.fit_min_counts:
            continue
        if min(np.sum(pi0_h[key]) for key in FRACTION_DRIVER_KEYS) < args.fit_min_counts:
            continue

        reference_fit, _, reference_variants, reference_exact = exact_fit_variants(
            data_h, dvcs_h, pi0_h, args, topology_key
        )
        topology_reference_fits[topology_key] = reference_fit
        topology_reference_variants[topology_key] = {
            name: {
                "success": value.success,
                "fraction_pi0": value.fraction_pi0,
                "fraction_pi0_err": value.fraction_pi0_err,
                "deviance": value.deviance,
                "ndf": value.ndf,
            }
            for name, value in reference_variants.items()
        }
        if reference_fit.success:
            topology_references[topology_key] = reference_exact
            definition = BinDefinition(
                -1,
                "FT" if topology_key == "CD_FT" else "FD",
                0,
                args.probe_E_min,
                args.probe_E_max,
                args.ft_theta_min if topology_key == "CD_FT" else args.fd_theta_min,
                args.ft_theta_max if topology_key == "CD_FT" else args.fd_theta_max,
            )
            plot_category_fit_summary(
                aggregate_fit_dir / f"{topology_key}_integrated.png",
                f"{period.label} {TOPOLOGY_LABELS[topology_key]}",
                definition,
                reference_fit,
            )

    fit_model_systematics: Dict[str, float] = {}
    fit_variants_metadata: Dict[str, object] = {}
    fit_support_coverage: Dict[str, Dict[str, float]] = {}
    category_fit_summaries: Dict[str, MultiFitSummary] = {}

    for definition in bins:
        index = definition.bin_id
        category_key = "FT" if definition.detector == "FT" else f"FD_sector_{definition.sector}"
        topology_details: Dict[str, object] = {}

        n_probe_missing_data = 0.0
        n_probe_missing_data_variance = 0.0
        n_probe_missing_mc = 0.0
        all_success = True
        weighted_fraction_numerator = 0.0
        weighted_fraction_denominator = 0.0
        total_deviance = 0.0
        total_ndf = 0
        maximum_model_spread = 0.0

        for topology_key in topology_keys_for_definition(definition):
            data_mask = (
                (data_one_photon_ids == index)
                & trial_topology_mask(data_one_photon_candidates, topology_key)
            )
            dvcs_mask = (
                (dvcs_template_ids == index)
                & trial_topology_mask(dvcs_background_template, topology_key)
            )
            pi0_mask = (
                (pi0_missing_probe_ids == index)
                & trial_topology_mask(pi0_missing_probe_mc, topology_key)
            )

            n_candidates, n_candidates_sumw2 = weighted_count_for_mask(
                data_one_photon_candidates, data_mask
            )
            n_pi0_mc_top, _ = weighted_count_for_mask(pi0_missing_probe_mc, pi0_mask)
            n_probe_missing_mc += n_pi0_mc_top

            data_h = histograms_for_mask(data_one_photon_candidates, data_mask, edges_by_var)
            dvcs_h = histograms_for_mask(dvcs_background_template, dvcs_mask, edges_by_var)
            pi0_h = histograms_for_mask(pi0_missing_probe_mc, pi0_mask, edges_by_var)

            support = {
                "data_Delta_phi": float(np.sum(data_h["Delta_phi"])),
                "data_theta_gamma_gamma": float(np.sum(data_h["theta_gamma_gamma"])),
                "data_pTmiss": float(np.sum(data_h["pTmiss"])),
                "dvcs_Delta_phi": float(np.sum(dvcs_h["Delta_phi"])),
                "pi0_Delta_phi": float(np.sum(pi0_h["Delta_phi"])),
            }

            reference_exact = topology_references.get(topology_key)
            if reference_exact is None:
                fit = MultiFitSummary(
                    False, math.nan, math.nan, math.nan, math.nan,
                    math.nan, math.nan, math.nan, 0,
                    f"No valid integrated reference fit for {topology_key}", {},
                )
                spread = math.nan
                variants = {"nominal": fit}
            elif topology_key == "CD_FT":
                # FT has one integrated output category, so use the full exact
                # topology fit rather than refitting frozen shapes.
                fit, spread, variants, _ = exact_fit_variants(
                    data_h, dvcs_h, pi0_h, args, topology_key
                )
            else:
                fit, spread, variants = frozen_variant_systematic(
                    data_h, reference_exact
                )

            if not fit.success or not math.isfinite(fit.fraction_pi0):
                all_success = False
                assigned = math.nan
                assigned_err = math.nan
            else:
                fraction = float(fit.fraction_pi0)
                fraction_err = (
                    float(fit.fraction_pi0_err)
                    if math.isfinite(fit.fraction_pi0_err) else 0.0
                )
                spread_value = float(spread) if math.isfinite(spread) else 0.0
                assigned = fraction * n_candidates
                fit_term = n_candidates * fraction_err
                model_term = n_candidates * spread_value
                count_term = fraction * math.sqrt(max(n_candidates_sumw2, 0.0))
                assigned_err = math.sqrt(
                    fit_term * fit_term
                    + model_term * model_term
                    + count_term * count_term
                )
                n_probe_missing_data += assigned
                n_probe_missing_data_variance += assigned_err * assigned_err
                weighted_fraction_numerator += fraction * n_candidates
                weighted_fraction_denominator += n_candidates
                total_deviance += fit.deviance
                total_ndf += fit.ndf
                maximum_model_spread = max(maximum_model_spread, spread_value)

            topology_details[topology_key] = {
                "label": TOPOLOGY_LABELS[topology_key],
                "selected_data_one_photon_candidates": n_candidates,
                "selected_pi0_mc_one_photon_candidates": n_pi0_mc_top,
                "support": support,
                "fit_success": fit.success,
                "fraction_pi0": fit.fraction_pi0,
                "fraction_pi0_err": fit.fraction_pi0_err,
                "assigned_pi0_data_yield": assigned,
                "assigned_pi0_data_yield_err": assigned_err,
                "model_spread_fraction": spread,
                "deviance": fit.deviance,
                "ndf": fit.ndf,
                "variants": {
                    name: {
                        "success": value.success,
                        "fraction_pi0": value.fraction_pi0,
                        "fraction_pi0_err": value.fraction_pi0_err,
                        "deviance": value.deviance,
                        "ndf": value.ndf,
                    }
                    for name, value in variants.items()
                },
            }

            if fit.success:
                plot_category_fit_summary(
                    aggregate_fit_dir / f"{category_key}_{topology_key}.png",
                    f"{period.label} {TOPOLOGY_LABELS[topology_key]}",
                    definition,
                    fit,
                )

        if not all_success:
            n_probe_missing_data = math.nan
            n_probe_missing_data_err = math.nan
        else:
            n_probe_missing_data_err = math.sqrt(n_probe_missing_data_variance)

        n_pass_data = float(pass_data_counts[index])
        n_pass_mc = float(pass_mc_counts[index])
        pass_data_err = math.sqrt(max(float(pass_data_sumw2[index]), 0.0))
        pass_mc_err = math.sqrt(max(float(pass_mc_sumw2[index]), 0.0))
        probe_missing_mc_err = math.sqrt(max(n_probe_missing_mc, 0.0))

        eff_data, eff_data_err = efficiency_and_error(
            n_pass_data, pass_data_err, n_probe_missing_data, n_probe_missing_data_err
        )
        eff_mc, eff_mc_err = efficiency_and_error(
            n_pass_mc, pass_mc_err, n_probe_missing_mc, probe_missing_mc_err
        )
        sf, sf_err = scale_factor_and_error(
            eff_data, eff_data_err, eff_mc, eff_mc_err
        )

        combined_fraction = (
            weighted_fraction_numerator / weighted_fraction_denominator
            if weighted_fraction_denominator > 0.0 else math.nan
        )
        combined_fit = MultiFitSummary(
            all_success,
            n_probe_missing_data,
            n_probe_missing_data_err,
            weighted_fraction_denominator - n_probe_missing_data
            if all_success else math.nan,
            math.nan,
            combined_fraction,
            math.nan,
            total_deviance,
            total_ndf,
            "Topology-separated production-style fit",
            {},
        )
        category_fit_summaries[category_key] = combined_fit
        fit_model_systematics[category_key] = maximum_model_spread
        fit_variants_metadata[category_key] = topology_details
        fit_support_coverage[category_key] = {
            key: value
            for topology in topology_details.values()
            for key, value in {
                f"{topology['label']} {name}": val
                for name, val in topology["support"].items()
            }.items()
        }

        rows.append(
            EfficiencyRow(
                period.key,
                period.label,
                definition.bin_id,
                definition.detector,
                definition.sector,
                definition.E_low,
                definition.E_high,
                definition.theta_low_deg,
                definition.theta_high_deg,
                n_pass_data,
                pass_data_err,
                n_probe_missing_data,
                n_probe_missing_data_err,
                eff_data,
                eff_data_err,
                n_pass_mc,
                pass_mc_err,
                n_probe_missing_mc,
                probe_missing_mc_err,
                eff_mc,
                eff_mc_err,
                sf,
                sf_err,
                all_success,
                total_deviance,
                total_ndf,
            )
        )

        log(
            f"{period.label} {category_key}: "
            f"topology-separated f_pi0={combined_fraction:.6g}, "
            f"assigned pi0 one-photon yield={n_probe_missing_data:,.1f}, "
            f"epgammagamma data={n_pass_data:,.1f}, "
            f"epgammagamma MC={n_pass_mc:,.1f}, "
            f"pi0-as-epgamma MC={n_probe_missing_mc:,.1f}"
        )

    plot_integrated_efficiency_summary(
        plot_dir / "integrated_efficiencies_and_scale_factors.png",
        period.label, rows,
    )
    audit_dir = period_dir / "pass_sample_audit"
    plot_matching_scale_factor_scan(
        audit_dir / "matching_cut_scale_factor_scan.png", period.label,
        pass_data_diag, pass_mc_diag, rows,
    )
    plot_matching_stability_maps(
        audit_dir / "matching_cut_stability_maps.png", period.label,
        pass_data_diag, pass_mc_diag, rows,
    )

    metadata = {
        "period": asdict(period),
        "pass_sample_audit": {"data": pass_data_diag, "mc": pass_mc_diag},
        "counts": {
            "epgammagamma_data_events_stored": pass_data.size(),
            "epgammagamma_mc_events_stored": pass_mc.size(),
            "epgammagamma_data_event_weight": float(np.sum(pass_data.weight)),
            "epgammagamma_mc_event_weight": float(np.sum(pass_mc.weight)),
            "data_one_photon_candidates": data_one_photon_candidates.size(),
            "dvcs_background_template_candidates": dvcs_background_template.size(),
            "pi0_missing_probe_mc_candidates": pi0_missing_probe_mc.size(),
            "n_bins": len(rows),
        },
        "input_file_identity": args.input_identity_records.get(period.key, {}),
        "selected_one_photon_support_audit": one_photon_support_audit,
        "one_photon_cutflow": {
            "data": data_one_photon_diag,
            "dvcs_background_template": dvcs_template_diag,
            "pi0_missing_probe_mc": pi0_missing_diag,
        },
        "fit_support_coverage": fit_support_coverage,
        "energy_support": energy_support,
        "mirror_policy": args.mirror_policy,
        "matching_scan_summary": {
            "data": matching_scan_summary(pass_data_diag),
            "mc": matching_scan_summary(pass_mc_diag),
        },
        "aaogen_truth_closure": inspect_truth_closure_availability(period.epgg_mc, args),
        "fit_architecture": {
            "fraction_drivers": list(FRACTION_DRIVER_KEYS),
            "diagnostic_only_variables": [s.key for s in FIT_VARIABLES if not s.fraction_driver],
            "constraint_variables": [],
            "support_fit_source": "plot_exclusivity_data_dvcs_pi0_mc.py exact fitter",
            "FD_morphology": "separate exact reference fits for (CD,FD) and (FD,FD), frozen independently for sector fractions",
            "FT_morphology": "exact exclusivity fit for the production (CD,FT) topology",
            "yield_normalization": "fitted pi0 fraction multiplied by the full selected data epgamma candidate count only after finite common-support validation",
            "mc_missing_probe_source": "production eppi0_bkg_aaogen_norad_*_epgamma.root sample from load_trees.cpp; no template fit is used for the MC one-photon pi0 count",
            "dvcs_mc_role": "shape-only DVCS/BH background template for the data one-photon decomposition; never enters the MC efficiency denominator",
            "tag_energy_policy": "0.4--probe-threshold tag fractions are recorded as diagnostics only; no minimum fraction is imposed",
            "fit_support_policy": "Delta_phi, theta_gamma_gamma and pTmiss determine f_pi0 on one common population; all other projections are diagnostic only",
            "common_selection_policy": "epgamma and epgammagamma share basic quality, -t, leading-photon opening angle, optional DIS cuts, leading-photon E>=2 GeV, and the same production topology definitions; topology-specific exclusivity and photon-matching cuts are diagnostic only",
            "production_global_cuts_enabled": bool(not args.disable_production_global_cuts),
        },
        "fit_model_systematics_on_fraction_pi0": fit_model_systematics,
        "fit_variants": fit_variants_metadata,
        "topology_reference_fits": {
            key: {
                "success": fit.success,
                "fraction_pi0": fit.fraction_pi0,
                "fraction_pi0_err": fit.fraction_pi0_err,
                "deviance": fit.deviance,
                "ndf": fit.ndf,
                "variants": topology_reference_variants.get(key, {}),
            }
            for key, fit in topology_reference_fits.items()
        },
        "fit_diagnostics": {
            category_key: {
                "mixture_fit_converged": fit.success,
                "fit_message": fit.message,
                "combined_fraction_pi0": fit.fraction_pi0,
                "combined_deviance": fit.deviance,
                "combined_ndf": fit.ndf,
                "topology_components": fit_variants_metadata.get(category_key, {}),
                "assigned_pi0_missing_probe_yield": float(rows[i].fail_data),
                "assigned_pi0_missing_probe_yield_err": float(rows[i].fail_data_err),
            }
            for i, (category_key, fit) in enumerate(category_fit_summaries.items())
        },
        "selected_energy_support": energy_support,
        "exclusivity_fitter_reference": {
            "module": str(Path(load_exclusivity_fitter_module().__file__).resolve()),
            "topology_reference_success": {
                key: bool(fit.success) for key, fit in topology_reference_fits.items()
            },
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



def summarize_tag_regions(trials: TrialArrays, split_GeV: float = 2.0) -> Dict[str, float]:
    """Summarize selected directed trials by observed-tag energy."""
    finite = np.isfinite(trials.tag_E) & np.isfinite(trials.weight)
    low = finite & (trials.tag_E < split_GeV)
    high = finite & (trials.tag_E >= split_GeV)
    return {
        "stored_trials": int(trials.size()),
        "weighted_total": float(np.sum(trials.weight[finite])),
        "weighted_low_tag": float(np.sum(trials.weight[low])),
        "weighted_high_tag": float(np.sum(trials.weight[high])),
        "low_tag_fraction": float(np.sum(trials.weight[low]) / np.sum(trials.weight[finite])) if np.sum(trials.weight[finite]) > 0 else math.nan,
    }


def plot_selected_energy_support(
    path: Path,
    period_label: str,
    pass_data: TrialArrays,
    pass_mc: TrialArrays,
    data_one_photon: TrialArrays,
    dvcs_template: TrialArrays,
    pi0_missing_probe: TrialArrays,
    args: argparse.Namespace,
) -> None:
    """Compare the actual leading-photon energy support of all five samples."""
    samples = [
        ("epgammagamma data", pass_data),
        ("epgammagamma AAOGEN", pass_mc),
        ("epgamma data", data_one_photon),
        ("epgamma DVCSGEN", dvcs_template),
        ("epgamma AAOGEN", pi0_missing_probe),
    ]
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    edges = np.linspace(args.probe_E_min, args.probe_E_max, 76)
    for label, trial in samples:
        if trial.size() == 0:
            continue
        ax.hist(
            trial.E,
            bins=edges,
            weights=trial.weight,
            histtype="step",
            density=True,
            linewidth=1.4,
            label=f"{label} (N={np.sum(trial.weight):,.0f})",
        )
    ax.axvline(args.probe_E_min, linestyle="--", linewidth=1.0)
    ax.set_xlabel("Leading reconstructed photon energy (GeV)")
    ax.set_ylabel("Unit-normalized entries")
    ax.set_title(f"{period_label}: common leading-photon support")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def audit_one_photon_file(path: str, role: str, period: PeriodConfig,
                          args: argparse.Namespace) -> Dict[str, object]:
    """Audit a production epgamma ROOT file without running the physics extraction."""
    record: Dict[str, object] = {
        "period": period.key,
        "period_label": period.label,
        "role": role,
        "path": path,
        "exists": Path(path).is_file(),
        "ok": False,
    }
    if not record["exists"]:
        record["status"] = "missing"
        return record
    stat = Path(path).stat()
    age_minutes = (time.time() - stat.st_mtime) / 60.0
    record.update({"size_bytes": int(stat.st_size), "mtime_unix": stat.st_mtime,
                   "age_minutes": age_minutes})
    try:
        with uproot.open(path) as root_file:
            if TREE_NAME not in root_file:
                record["status"] = f"missing_tree_{TREE_NAME}"
                return record
            tree = root_file[TREE_NAME]
            record["entries"] = int(tree.num_entries)
            keys = {str(k).split(";")[0] for k in tree.keys()}
            branch = next((name for name in ALIASES["g1_E"] if name in keys), None)
            record["observed_photon_energy_branch"] = branch
            if branch is None:
                record["status"] = "missing_observed_photon_energy_branch"
                return record
        counts = {"finite": 0, "below_tag_min": 0, "tag_to_2": 0, "above_2": 0}
        hist_edges = np.linspace(0.0, max(args.tag_E_max, 10.0), 101)
        hist = np.zeros(hist_edges.size - 1, dtype=np.int64)
        for arrays in uproot.iterate(f"{path}:{TREE_NAME}", expressions=[branch],
                                     step_size=args.step_size, library="np"):
            values = np.asarray(arrays[branch], dtype=float)
            values = values[np.isfinite(values)]
            counts["finite"] += int(values.size)
            counts["below_tag_min"] += int(np.count_nonzero(values < args.tag_E_min))
            counts["tag_to_2"] += int(np.count_nonzero((values >= args.tag_E_min) & (values < 2.0)))
            counts["above_2"] += int(np.count_nonzero(values >= 2.0))
            hist += np.histogram(values, bins=hist_edges)[0]
        record.update(counts)
        record["fraction_tag_to_2"] = counts["tag_to_2"] / counts["finite"] if counts["finite"] else math.nan
        record["hist_edges_GeV"] = hist_edges.tolist()
        record["hist_counts"] = hist.tolist()
        checks = {
            "minimum_entries": int(record["entries"]) >= args.preflight_min_entries,
            "stable_age": args.preflight_stable_age_min < 0 or age_minutes >= args.preflight_stable_age_min,
        }
        record["checks"] = checks
        record["ok"] = all(checks.values())
        record["status"] = "ok" if record["ok"] else "failed_checks"
    except Exception as exc:
        record["status"] = "error"
        record["error"] = repr(exc)
    return record


def write_preflight_outputs(output_dir: Path, periods: Sequence[PeriodConfig],
                            args: argparse.Namespace) -> Dict[str, object]:
    """Write a complete production-input readiness and energy-support audit."""
    preflight_dir = output_dir / "preflight"
    preflight_dir.mkdir(parents=True, exist_ok=True)
    records: List[Dict[str, object]] = []
    for period in periods:
        roles = (("data", period.epg_data), ("dvcs_mc", period.dvcs_mc),
                 ("pi0_as_epg_mc", period.pi0_as_epg_mc))
        for role, path in roles:
            log(f"Preflight: {period.label} {role}: {path}")
            records.append(audit_one_photon_file(path, role, period, args))
    payload = {
        "created_unix_time": time.time(),
        "tag_energy_min_GeV": args.tag_E_min,
        "probe_energy_min_GeV": args.probe_E_min,
        "all_ok": bool(records) and all(bool(r.get("ok")) for r in records),
        "records": records,
    }
    with open(preflight_dir / "production_epgamma_preflight.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, allow_nan=True)
    flat_rows = []
    for r in records:
        flat_rows.append({k: v for k, v in r.items()
                          if k not in {"hist_edges_GeV", "hist_counts", "checks"}} |
                         {f"check_{k}": v for k, v in r.get("checks", {}).items()})
    write_rows(preflight_dir / "production_epgamma_preflight.csv", flat_rows)

    for period in periods:
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        for role in ("data", "dvcs_mc", "pi0_as_epg_mc"):
            rec = next((r for r in records if r["period"] == period.key and r["role"] == role), None)
            if not rec or "hist_counts" not in rec:
                continue
            edges = np.asarray(rec["hist_edges_GeV"], dtype=float)
            counts = np.asarray(rec["hist_counts"], dtype=float)
            total = counts.sum()
            if total > 0:
                counts /= total
            ax.stairs(counts, edges, label=role)
        ax.axvline(args.tag_E_min, linestyle="--", linewidth=1.0)
        ax.axvline(2.0, linestyle="--", linewidth=1.0)
        ax.set_xlim(0.0, min(args.tag_E_max, 6.0))
        ax.set_xlabel("Observed photon energy (GeV)")
        ax.set_ylabel("Fraction per bin")
        ax.set_title(f"{period.label}: production epgamma input support")
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(preflight_dir / f"{period.key}_observed_photon_energy_support.png", dpi=180)
        plt.close(fig)
    return payload

def main() -> int:
    args = parse_args()
    if args.probe_E_min < 2.0:
        log("WARNING: --probe-E-min is below the production DVCS leading-photon threshold.")
    log(
        f"Direct sample definition: leading reconstructed photon "
        f"{args.probe_E_min:g} <= E_gamma1 < {args.probe_E_max:g} GeV; "
        f"second photon in epgammagamma "
        f"{args.tag_E_min:g} <= E_gamma2 < {args.tag_E_max:g} GeV"
    )
    if args.workers < 1:
        raise ValueError("--workers must be at least 1")
    # endif
    args.workers = min(args.workers, MAX_WORKERS)

    selected_keys = set(args.period or [p.key for p in PERIODS])
    periods = [p for p in PERIODS if p.key in selected_keys]
    periods = apply_path_overrides(periods, args)

    log("Input ROOT-tree roles:")
    for period in periods:
        log(f"  {period.label} epgamma data mixture: {period.epg_data}")
        log(f"  {period.label} DVCS/BH epgamma template: {period.dvcs_mc}")
        log(f"  {period.label} AAOGEN pi0-as-epgamma template/MC one-photon count: {period.pi0_as_epg_mc}")
        log(f"  {period.label} epgammagamma data two-photon count: {period.epgg_data}")
        log(f"  {period.label} AAOGEN epgammagamma MC two-photon count: {period.epgg_mc}")

    if args.inspect_branches:
        inspect_selected_branches(periods)
        return 0
    # endif

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    preflight_payload = write_preflight_outputs(output_dir, periods, args)
    if args.preflight_only:
        log(f"Preflight complete: all_ok={preflight_payload['all_ok']}")
        return 0 if preflight_payload["all_ok"] else 2

    if not preflight_payload["all_ok"]:
        raise RuntimeError(
            "Production epgamma input preflight failed. Inspect output/photon_efficiency_study/preflight; "
            "use --preflight-only while files are still being produced."
        )

    for period in periods:
        for path in (period.epg_data, period.dvcs_mc, period.pi0_as_epg_mc,
                     period.epgg_data, period.epgg_mc):
            require_file(path)
        # endfor
    # endfor

    identity_report = validate_period_input_identities(periods, args)
    args.input_identity_records = identity_report["files"]
    with open(output_dir / "input_file_identity.json", "w", encoding="utf-8") as handle:
        json.dump(identity_report, handle, indent=2)

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
    plot_all_period_scale_factor_summary(output_dir / "all_periods_integrated_scale_factors.png", all_rows)
    plot_all_period_efficiency_summary(output_dir / "all_periods_integrated_efficiencies.png", all_rows)

    json_payload = {
        "schema_version": 3,
        "description": "RGA exclusive-pi0 event-migration photon-reconstruction data/MC scale factors.",
        "formula": "epsilon = N_pi0_epgammagamma / (N_pi0_epgammagamma + N_pi0_epgamma); S_gamma = epsilon_data / epsilon_mc",
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
            "Each accepted epgammagamma event contributes one two-photon event. The energy-ordered leading photon satisfies probe-E-min and the second reconstructed photon satisfies tag-E-min.",
            "The nominal thresholds are E_gamma1>2 GeV and E_gamma2>0.40 GeV, matching the production DVCS photon requirement and the loose second-photon control-sample threshold.",
            "The one-photon pi0 yield in data is obtained only from the DVCSGEN plus AAOGEN template decomposition of the ordinary epgamma distributions.",
            "The epgamma data, DVCSGEN, and AAOGEN-as-epgamma samples use the same measured leading-photon requirement; no predicted missing photon is constructed.",
            "The extraction is integrated in categories of the actual leading reconstructed photon: FT and each FD sector.",
            "A complete data/MC passing-sample cut-flow audit is written for every period.",
            "Additional production-style (-t1), z, and predicted electron-photon opening-angle cuts are opt-in via --enable-production-global-cuts and are disabled by default. Sp18 Out sector-quality exclusions are intentionally not applied in this diagnostic study.",
            "Input ROOT identity records and duplicate-period checks are written before processing.",
            "Matching-cut scans are optional (--enable-matching-scans) because they are expensive and do not affect the nominal result; scale-factor stability, mirror-category migration, and fit residual/deviance diagnostics remain enabled.",
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
