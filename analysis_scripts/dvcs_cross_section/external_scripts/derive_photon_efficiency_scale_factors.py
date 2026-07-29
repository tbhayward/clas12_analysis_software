#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v2.py

Exploratory RGA photon-efficiency tag-and-probe study for the DVCS analysis.

Physics definition
------------------
The script measures, separately in data and AAOGEN Monte Carlo,

    epsilon_b = N_pass,b / (N_pass,b + N_fail,b),

where b denotes the predicted probe-photon detector/kinematic bin.

* An exclusive epgamma-gamma event contributes two directed PASS trials:
  gamma_1 is used as the tag and gamma_2 as the probe, and vice versa.
* An epgamma event attributed to pi0 production contributes one FAIL trial:
  the reconstructed photon is the tag and the missing photon is predicted from

      p_probe^pred = p_beam + p_target - p_e - p_p - p_tag.

No equal-efficiency assumption for the two pi0 photons is made. Every directed
trial is assigned using the kinematics of its own predicted probe.

The residual photon-efficiency correction is the data/MC scale factor

    S_gamma,b = epsilon_data,b / epsilon_MC,b.

The first implementation uses:

* FT: adaptive bins in predicted E_gamma, integrated over the configured FT
  angular range.
* FD: separate sectors 1--6, with adaptive rectangular bins in predicted
  E_gamma and theta_gamma.
* epgamma data FAIL yields: a non-negative binned Poisson template fit of the
  Mx2_1 = M_X^2(ep) distribution using DVCSGEN and pi0-as-DVCS AAOGEN.
* epgamma MC FAIL yields: direct pi0-as-DVCS AAOGEN counts.
* epgamma-gamma PASS yields: direct exclusive-pi0 data and AAOGEN counts after
  configurable broad exclusivity cuts.  This is deliberately simple in v1;
  later iterations can replace it with a simultaneous pass/fail nuisance fit.

Important limitations of v1
---------------------------
1. The script measures a tag-conditioned efficiency.  AAOGEN truth closure is
   not possible unless generated/truth photon branches are available; v1 writes
   all ingredients needed for reconstructed tag-and-probe closure.
2. Direct epgamma-gamma data are treated as pi0 signal after the supplied broad
   cuts.  The output includes Mx2_1 diagnostics so residual backgrounds can be
   evaluated before production use.
3. The FT/FD prediction uses configurable theta windows.  It intentionally does
   not attempt a detailed FT face or ECAL local-coordinate projection yet.
4. This script derives scale factors only.  It does not yet modify DVCSGEN
   acceptance or the AAOGEN one-photon/two-photon migration ratio.

Dependencies
------------
Python 3, numpy, scipy, matplotlib, uproot.

Typical usage
-------------

  python3 derive_photon_efficiency_scale_factors_v2.py \
      --period fa18_inb --workers 1

Run all periods with up to seven workers:

  python3 derive_photon_efficiency_scale_factors_v2.py --workers 7

Use --inspect-branches first if the local trees use unexpected branch names.
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
    "proton_detector": ("detector1", "p_detector", "proton_detector"),
    "Mx2_1": ("Mx2_1", "Mx2_ep", "Mx2_x1", "Mx2_proton", "Mx2_p"),
    "Mx2": ("Mx2", "Mx2_epg", "Mx2_epgamma", "Mx2_eppi0", "Mx2_epi0"),
    "pTmiss": ("pTmiss", "ptmiss", "pT_miss"),
    "Delta_phi": ("Delta_phi", "delta_phi", "dphi"),
    "theta_pi0_pi0": ("theta_pi0_pi0", "theta_gamma_gamma"),
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
            weight=np.empty(0, dtype=np.float32),
        )

    def size(self) -> int:
        return int(self.E.size)


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
    fail_fit_success: bool
    fail_fit_deviance: float
    fail_fit_ndf: int


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
    parser.add_argument("--ft-theta-max", type=float, default=4.5)
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
    parser.add_argument("--pass-delta-phi-window", type=float, default=0.35,
                        help="Require |Delta_phi-pi| below this value, if branch exists.")
    parser.add_argument("--probe-match-angle-max-deg", type=float, default=2.0,
                        help="Maximum angle between predicted and observed probe in epgamma-gamma events.")
    parser.add_argument("--probe-match-relative-E-max", type=float, default=0.35,
                        help="Maximum |Eobs-Epred|/Epred for a passing probe.")
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
    if resolved.get("Delta_phi") is not None:
        dphi = finite_array(arrays, resolved["Delta_phi"])
        mask &= np.isfinite(dphi) & (np.abs(dphi - math.pi) < args.pass_delta_phi_window)
    # endif
    return mask


def read_pass_trials(path: str, beam_energy: float, args: argparse.Namespace) -> TrialArrays:
    logical = [
        "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi",
        "g1_E", "g1_theta", "g1_phi", "g2_E", "g2_theta", "g2_phi",
        "Mx2_1", "Mx2", "pTmiss", "Delta_phi", "fiducial_status",
    ]
    optional = ["Mx2_1", "Mx2", "pTmiss", "Delta_phi", "fiducial_status"]
    resolved = resolve_branches(path, logical, optional)
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    chunks: List[TrialArrays] = []
    seen = 0

    log(f"Reading PASS trials from {path}")
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

        base = pass_sample_mask(arrays, resolved, args)
        e_p = finite_array(arrays, resolved["e_p"])
        e_theta = finite_array(arrays, resolved["e_theta"])
        e_phi = finite_array(arrays, resolved["e_phi"])
        p_p = finite_array(arrays, resolved["p_p"])
        p_theta = finite_array(arrays, resolved["p_theta"])
        p_phi = finite_array(arrays, resolved["p_phi"])
        g1_E = finite_array(arrays, resolved["g1_E"])
        g1_theta = finite_array(arrays, resolved["g1_theta"])
        g1_phi = finite_array(arrays, resolved["g1_phi"])
        g2_E = finite_array(arrays, resolved["g2_E"])
        g2_theta = finite_array(arrays, resolved["g2_theta"])
        g2_phi = finite_array(arrays, resolved["g2_phi"])
        mx2_ep = finite_array(arrays, resolved.get("Mx2_1"))

        # Directed trial 1: gamma1 tag, gamma2 probe.  The bin coordinates use
        # the kinematically predicted probe, not the observed probe.
        pred1_E, pred1_th, pred1_ph, pred1_m2, pred1_p = predicted_probe(
            beam_energy, e_p, e_theta, e_phi, p_p, p_theta, p_phi,
            g1_E, g1_theta, g1_phi,
        )
        # Directed trial 2: gamma2 tag, gamma1 probe.
        pred2_E, pred2_th, pred2_ph, pred2_m2, pred2_p = predicted_probe(
            beam_energy, e_p, e_theta, e_phi, p_p, p_theta, p_phi,
            g2_E, g2_theta, g2_phi,
        )

        for pred_E, pred_th, pred_ph, pred_m2, pred_p, obs_E, obs_th, obs_ph in (
            (pred1_E, pred1_th, pred1_ph, pred1_m2, pred1_p, g2_E, g2_theta, g2_phi),
            (pred2_E, pred2_th, pred2_ph, pred2_m2, pred2_p, g1_E, g1_theta, g1_phi),
        ):
            theta_deg = np.degrees(pred_th)
            detector, sector = classify_predicted_detector(theta_deg, pred_ph, args)
            cos_opening = (
                np.sin(pred_th) * np.sin(obs_th) * np.cos(pred_ph - obs_ph)
                + np.cos(pred_th) * np.cos(obs_th)
            )
            match_angle_deg = np.degrees(np.arccos(np.clip(cos_opening, -1.0, 1.0)))
            relative_E_residual = np.abs(obs_E - pred_E) / np.maximum(pred_E, 1e-9)
            good = (
                base
                & np.isfinite(pred_E) & np.isfinite(pred_p)
                & np.isfinite(theta_deg) & np.isfinite(pred_ph)
                & np.isfinite(obs_E) & np.isfinite(obs_th) & np.isfinite(obs_ph)
                & np.isfinite(match_angle_deg) & np.isfinite(relative_E_residual)
                & (pred_E >= args.probe_E_min) & (pred_E < args.probe_E_max)
                & (pred_p > 0.0)
                & (np.abs(pred_E - pred_p) < 0.30)
                & (np.abs(pred_m2) < 0.20)
                & (match_angle_deg < args.probe_match_angle_max_deg)
                & (relative_E_residual < args.probe_match_relative_E_max)
                & (detector >= 0)
            )
            chunks.append(TrialArrays(
                E=pred_E[good].astype(np.float32, copy=False),
                theta_deg=theta_deg[good].astype(np.float32, copy=False),
                phi_rad=pred_ph[good].astype(np.float32, copy=False),
                detector=detector[good].astype(np.int8, copy=False),
                sector=sector[good].astype(np.int8, copy=False),
                mx2_ep=mx2_ep[good].astype(np.float32, copy=False),
                weight=np.ones(np.count_nonzero(good), dtype=np.float32),
            ))
        # endfor
    # endfor
    return concatenate_trials(chunks)


def read_fail_trials(path: str, beam_energy: float, args: argparse.Namespace) -> TrialArrays:
    logical = [
        "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi",
        "g1_E", "g1_theta", "g1_phi", "Mx2_1", "fiducial_status",
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
    all_E = np.concatenate([pass_mc.E, fail_mc.E])
    all_theta = np.concatenate([pass_mc.theta_deg, fail_mc.theta_deg])
    all_detector = np.concatenate([pass_mc.detector, fail_mc.detector])
    all_sector = np.concatenate([pass_mc.sector, fail_mc.sector])

    bins: List[BinDefinition] = []
    next_id = 0

    ft_mask = all_detector == 0
    ft_edges = strict_quantile_edges(all_E[ft_mask], args.ft_energy_bins,
                                     args.probe_E_min, args.probe_E_max)
    for elo, ehi in zip(ft_edges[:-1], ft_edges[1:]):
        bins.append(BinDefinition(next_id, "FT", 0, float(elo), float(ehi),
                                  args.ft_theta_min, args.ft_theta_max))
        next_id += 1
    # endfor

    for sector in range(1, 7):
        smask = (all_detector == 1) & (all_sector == sector)
        e_edges = strict_quantile_edges(all_E[smask], args.fd_energy_bins,
                                        args.probe_E_min, args.probe_E_max)
        for elo, ehi in zip(e_edges[:-1], e_edges[1:]):
            emask = smask & (all_E >= elo) & (all_E < ehi)
            t_edges = strict_quantile_edges(all_theta[emask], args.fd_theta_bins,
                                            args.fd_theta_min, args.fd_theta_max)
            for tlo, thi in zip(t_edges[:-1], t_edges[1:]):
                count = int(np.count_nonzero(emask & (all_theta >= tlo) & (all_theta < thi)))
                if count < args.min_probes_per_bin and len(t_edges) > 2:
                    continue
                # endif
                bins.append(BinDefinition(next_id, "FD", sector, float(elo), float(ehi),
                                          float(tlo), float(thi)))
                next_id += 1
            # endfor
        # endfor
    # endfor
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


def fit_two_template_histograms(
    data_counts: np.ndarray,
    dvcs_counts: np.ndarray,
    pi0_counts: np.ndarray,
    edges: np.ndarray,
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

    dvcs_shape = dvcs_counts.astype(float) / dvcs_counts.sum()
    pi0_shape = pi0_counts.astype(float) / pi0_counts.sum()
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


def process_period(period: PeriodConfig, args_dict: Mapping[str, object]) -> Tuple[str, List[Dict[str, object]], Dict[str, object]]:
    period_start = time.perf_counter()
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    fit_dir = period_dir / "fail_fits"
    plot_dir = period_dir / "plots"
    fit_dir.mkdir(parents=True, exist_ok=True)
    plot_dir.mkdir(parents=True, exist_ok=True)

    pass_data = read_pass_trials(period.epgg_data, period.beam_energy_GeV, args)
    pass_mc = read_pass_trials(period.epgg_mc, period.beam_energy_GeV, args)
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
    mx2_edges = np.linspace(args.mx2_ep_min, args.mx2_ep_max, args.mx2_ep_bins + 1)

    log(f"{period.label}: assigning trials to {n_bins} adaptive bins")
    pass_data_ids = assign_bin_ids(pass_data, bins)
    pass_mc_ids = assign_bin_ids(pass_mc, bins)
    fail_data_ids = assign_bin_ids(fail_data, bins)
    fail_dvcs_ids = assign_bin_ids(fail_dvcs_mc, bins)
    fail_pi0_ids = assign_bin_ids(fail_pi0_mc, bins)

    pass_data_counts = weighted_counts_by_bin(pass_data, pass_data_ids, n_bins)
    pass_mc_counts = weighted_counts_by_bin(pass_mc, pass_mc_ids, n_bins)
    fail_pi0_counts = weighted_counts_by_bin(fail_pi0_mc, fail_pi0_ids, n_bins)

    data_mx2_hists = bulk_mx2_histograms(fail_data, fail_data_ids, n_bins, mx2_edges)
    dvcs_mx2_hists = bulk_mx2_histograms(fail_dvcs_mc, fail_dvcs_ids, n_bins, mx2_edges)
    pi0_mx2_hists = bulk_mx2_histograms(fail_pi0_mc, fail_pi0_ids, n_bins, mx2_edges)

    for definition in bins:
        index = definition.bin_id
        n_pass_data = float(pass_data_counts[index])
        n_pass_mc = float(pass_mc_counts[index])
        n_fail_mc = float(fail_pi0_counts[index])

        fit = fit_two_template_histograms(
            data_mx2_hists[index],
            dvcs_mx2_hists[index],
            pi0_mx2_hists[index],
            mx2_edges,
        )
        n_fail_data = fit.n_pi0
        n_fail_data_err = fit.n_pi0_err

        pass_data_err = math.sqrt(max(n_pass_data, 0.0))
        pass_mc_err = math.sqrt(max(n_pass_mc, 0.0))
        fail_mc_err = math.sqrt(max(n_fail_mc, 0.0))

        eff_data, eff_data_err = efficiency_and_error(
            n_pass_data, pass_data_err, n_fail_data, n_fail_data_err
        )
        eff_mc, eff_mc_err = efficiency_and_error(
            n_pass_mc, pass_mc_err, n_fail_mc, fail_mc_err
        )
        sf, sf_err = scale_factor_and_error(eff_data, eff_data_err, eff_mc, eff_mc_err)

        row = EfficiencyRow(
            period.key, period.label, definition.bin_id, definition.detector,
            definition.sector, definition.E_low, definition.E_high,
            definition.theta_low_deg, definition.theta_high_deg,
            n_pass_data, pass_data_err, n_fail_data, n_fail_data_err,
            eff_data, eff_data_err, n_pass_mc, pass_mc_err, n_fail_mc, fail_mc_err,
            eff_mc, eff_mc_err, sf, sf_err, fit.success, fit.deviance, fit.ndf,
        )
        rows.append(row)

        fit_name = f"bin_{definition.bin_id:03d}_{definition.detector}_s{definition.sector}.png"
        title = (
            f"{period.label}, {definition.detector} sector {definition.sector}, "
            f"E {definition.E_low:.2f}-{definition.E_high:.2f} GeV, "
            f"theta {definition.theta_low_deg:.1f}-{definition.theta_high_deg:.1f} deg"
        )
        plot_fit(fit_dir / fit_name, title, fit)
    # endfor

    plot_scale_factors(plot_dir / "scale_factors_FT.png", period.label, rows, "FT", 0)
    plot_efficiencies(plot_dir / "efficiencies_FT.png", period.label, rows, "FT", 0)
    for sector in range(1, 7):
        plot_scale_factors(plot_dir / f"scale_factors_FD_sector_{sector}.png",
                           period.label, rows, "FD", sector)
        plot_efficiencies(plot_dir / f"efficiencies_FD_sector_{sector}.png",
                          period.label, rows, "FD", sector)
    # endfor

    metadata = {
        "period": asdict(period),
        "counts": {
            "pass_data_directed_trials": pass_data.size(),
            "pass_mc_directed_trials": pass_mc.size(),
            "fail_data_candidates": fail_data.size(),
            "fail_dvcs_mc_candidates": fail_dvcs_mc.size(),
            "fail_pi0_mc_candidates": fail_pi0_mc.size(),
            "n_bins": len(rows),
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
            "FT is binned only in predicted probe energy.",
            "FD is binned by predicted sector, energy and polar angle.",
            "No equal-efficiency approximation is used.",
            "Scale factors are not yet propagated to DVCS acceptance or pi0 migration.",
            "Period processing is parallelized with a hard maximum of seven workers.",
            "Trials are assigned to adaptive bins once and Mx2 histograms are filled in bulk.",
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
