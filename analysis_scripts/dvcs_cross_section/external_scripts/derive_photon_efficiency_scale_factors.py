#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v9_epgammagamma_truth_partner_audit.py

Compact truth-partner audit for the RGA photon-efficiency study.

Purpose
-------
This release does not fit templates and does not calculate efficiencies or
scale factors.  It validates the central tag-and-probe geometry before the
production extraction is restored.

For each directed epgamma opportunity:

    1. use the measured epgamma photon as the tag;
    2. calculate the missing four-vector X;
    3. find the corresponding native epgammagamma record;
    4. reconstruct both daughter-photon laboratory four-vectors from the
       quantities stored in the native pi0 tree;
    5. identify which reconstructed daughter is the epgamma tag;
    6. define the opposite daughter as the true reconstructed partner;
    7. compare X directly with that partner.

The native tree stores Trento azimuths gamma_phi1 and gamma_phi2, not laboratory
azimuths.  This script converts them correctly by constructing the Trento basis
around the virtual-photon direction.  The polar angle about q is solved from
the stored electron-photon opening angle.  Photon energies are then obtained
from four-vector closure to the stored reconstructed pi0.

No DVCS kinematic, exclusivity, fiducial, detector-acceptance, or pi0-mass cuts
are imposed.  The only opportunity requirements are finite reconstructed
electron/proton/tag kinematics, E_tag >= 0.4 GeV, and predicted E_X >= 2 GeV.

Identity matching
-----------------
Data:
    exact (runnum, evnum)

AAOGEN MC:
    rounded reconstructed electron/proton signature
    (e_p, e_theta, e_phi, p_p, p_theta, p_phi)

Outputs
-------
output/photon_efficiency_study/truth_partner_audit/

    truth_partner_audit_summary.csv
    truth_partner_audit_summary.json
    input_manifest.json
    <period>/
        data_truth_partner_residuals.png
        aaogen_truth_partner_residuals.png
        data_truth_partner_sample.csv
        aaogen_truth_partner_sample.csv
        metadata.json

Only compact exact counters, histogram bin contents, and a small deterministic
reservoir sample are written.  No full per-opportunity table is produced.
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
from typing import Dict, List, Mapping, Optional, Sequence, Tuple

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", f"/tmp/matplotlib-{os.getuid()}")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import uproot


TREE_NAME = "PhysicsEvents"
MAX_WORKERS = 8
DEFAULT_STEP_SIZE = "200 MB"
DEFAULT_OUTPUT_DIR = "output/photon_efficiency_study/truth_partner_audit"

PROTON_MASS_GEV = 0.9382720813
ELECTRON_MASS_GEV = 0.00051099895


@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy_GeV: float
    epg_data: str
    epgg_data: str
    pi0_epg_mc: str
    pi0_epgg_mc: str


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        "fa18_inb", "Fa18 Inb", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "fa18_out", "Fa18 Out", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_50nA_10604MeV.root",
    ),
    PeriodConfig(
        "sp18_inb", "Sp18 Inb", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_50nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp18_out", "Sp18 Out", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_45nA_10594MeV.root",
    ),
    PeriodConfig(
        "sp19_inb", "Sp19 Inb", 10.200,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
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
    "g_E": ("p2_p", "g1_p", "g1_E", "gamma1_p", "photon1_p"),
    "g_theta": ("p2_theta", "g1_theta", "gamma1_theta", "photon1_theta"),
    "g_phi": ("p2_phi", "g1_phi", "gamma1_phi", "photon1_phi"),
    "pi0_p": ("p2_p", "pi0_p"),
    "pi0_theta": ("p2_theta", "pi0_theta"),
    "pi0_phi": ("p2_phi", "pi0_phi"),
    "pi0_mass": ("Mh_gammagamma",),
    "open1": ("open_angle_egamma1",),
    "open2": ("open_angle_egamma2",),
    "trento1": ("gamma_phi1",),
    "trento2": ("gamma_phi2",),
    "det1": ("detector_gamma1",),
    "det2": ("detector_gamma2",),
}


@dataclass
class EpgRecords:
    runnum: np.ndarray
    eventnum: np.ndarray
    e_p: np.ndarray
    e_theta: np.ndarray
    e_phi: np.ndarray
    p_p: np.ndarray
    p_theta: np.ndarray
    p_phi: np.ndarray
    g_E: np.ndarray
    g_theta: np.ndarray
    g_phi: np.ndarray
    probe_E: np.ndarray
    probe_theta: np.ndarray
    probe_phi: np.ndarray
    opportunity: np.ndarray

    def size(self) -> int:
        return int(self.g_E.size)


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


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Audit epgammagamma-derived truth partners for epgamma tags."
    )
    parser.add_argument("--period", action="append", choices=[p.key for p in PERIODS])
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--tag-E-min", type=float, default=0.40)
    parser.add_argument("--probe-E-min", type=float, default=2.00)
    parser.add_argument("--mc-signature-decimals", type=int, default=10)
    parser.add_argument("--tag-angle-max-deg", type=float, default=3.0)
    parser.add_argument("--tag-relative-E-max", type=float, default=0.35)
    parser.add_argument("--closure-max", type=float, default=0.10)
    parser.add_argument("--minimum-photon-E", type=float, default=0.20)
    parser.add_argument("--diagnostic-sample-size", type=int, default=5000)
    parser.add_argument("--diagnostic-seed", type=int, default=20260801)
    return parser.parse_args()


def selected_periods(args: argparse.Namespace) -> List[PeriodConfig]:
    if not args.period:
        return list(PERIODS)
    # endif
    requested = set(args.period)
    return [period for period in PERIODS if period.key in requested]


def require_tree(path: str) -> Tuple[int, List[str]]:
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    # endif
    with uproot.open(path) as root_file:
        if TREE_NAME not in root_file:
            raise KeyError(f"{TREE_NAME} absent from {path}")
        # endif
        tree = root_file[TREE_NAME]
        return int(tree.num_entries), [str(key) for key in tree.keys()]
    # endwith


def resolve(
    path: str,
    logical: Sequence[str],
    optional: Sequence[str] = (),
) -> Dict[str, Optional[str]]:
    _, keys = require_tree(path)
    available = set(keys)
    optional_set = set(optional)
    result: Dict[str, Optional[str]] = {}
    missing: List[str] = []
    for name in logical:
        branch = next(
            (
                candidate
                for candidate in ALIASES.get(name, (name,))
                if candidate in available
            ),
            None,
        )
        result[name] = branch
        if branch is None and name not in optional_set:
            missing.append(f"{name}: tried {ALIASES.get(name, (name,))}")
        # endif
    # endfor
    if missing:
        raise KeyError(
            f"Missing branches in {path}:\n  " + "\n  ".join(missing)
        )
    # endif
    return result


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


def unit_vector(theta: float, phi: float) -> np.ndarray:
    return np.asarray(
        [
            math.sin(theta) * math.cos(phi),
            math.sin(theta) * math.sin(phi),
            math.cos(theta),
        ],
        dtype=float,
    )


def spherical(
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


def massive(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
    mass: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical(momentum, theta, phi)
    return np.sqrt(np.maximum(momentum**2 + mass**2, 0.0)), px, py, pz


def reconstruct_probe(
    beam_energy: float,
    e_p: np.ndarray,
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    p_p: np.ndarray,
    p_theta: np.ndarray,
    p_phi: np.ndarray,
    g_E: np.ndarray,
    g_theta: np.ndarray,
    g_phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    e_E, e_px, e_py, e_pz = massive(
        e_p, e_theta, e_phi, ELECTRON_MASS_GEV
    )
    p_E, p_px, p_py, p_pz = massive(
        p_p, p_theta, p_phi, PROTON_MASS_GEV
    )
    g_px, g_py, g_pz = spherical(g_E, g_theta, g_phi)

    probe_E = beam_energy + PROTON_MASS_GEV - e_E - p_E - g_E
    probe_px = -e_px - p_px - g_px
    probe_py = -e_py - p_py - g_py
    probe_pz = beam_energy - e_pz - p_pz - g_pz
    probe_pt = np.sqrt(np.maximum(probe_px**2 + probe_py**2, 0.0))
    probe_theta = np.arctan2(probe_pt, probe_pz)
    probe_phi = np.mod(np.arctan2(probe_py, probe_px), 2.0 * math.pi)
    return probe_E, probe_theta, probe_phi


def read_epg(
    path: str,
    beam_energy: float,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[EpgRecords, Dict[str, int]]:
    logical = (
        "runnum", "eventnum",
        "e_p", "e_theta", "e_phi",
        "p_p", "p_theta", "p_phi",
        "g_E", "g_theta", "g_phi",
    )
    resolved = resolve(path, logical, optional=("runnum", "eventnum"))
    expressions = sorted(
        {branch for branch in resolved.values() if branch is not None}
    )
    store: Dict[str, List[np.ndarray]] = {
        name: [] for name in EpgRecords.__dataclass_fields__
    }
    counts = {
        "tree_entries": 0,
        "finite_epgamma_entries": 0,
        "photon_candidates_E_above_threshold": 0,
        "directed_opportunities_probe_E_above_threshold": 0,
    }

    seen = 0
    log(f"Reading {mode} epgamma records from {path}")
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
        g_E = arr("g_E")
        g_theta = arr("g_theta")
        g_phi = np.mod(arr("g_phi"), 2.0 * math.pi)

        finite = (
            np.isfinite(e_p)
            & np.isfinite(e_theta)
            & np.isfinite(e_phi)
            & np.isfinite(p_p)
            & np.isfinite(p_theta)
            & np.isfinite(p_phi)
            & np.isfinite(g_E)
            & np.isfinite(g_theta)
            & np.isfinite(g_phi)
            & (e_p > 0.0)
            & (p_p > 0.0)
            & (g_E > 0.0)
        )
        counts["finite_epgamma_entries"] += int(np.count_nonzero(finite))

        candidate = finite & (g_E >= args.tag_E_min)
        counts["photon_candidates_E_above_threshold"] += int(
            np.count_nonzero(candidate)
        )

        probe_E, probe_theta, probe_phi = reconstruct_probe(
            beam_energy,
            e_p, e_theta, e_phi,
            p_p, p_theta, p_phi,
            g_E, g_theta, g_phi,
        )
        opportunity = (
            candidate
            & np.isfinite(probe_E)
            & np.isfinite(probe_theta)
            & np.isfinite(probe_phi)
            & (probe_E >= args.probe_E_min)
        )
        counts["directed_opportunities_probe_E_above_threshold"] += int(
            np.count_nonzero(opportunity)
        )

        values = {
            "runnum": arr("runnum", -1, np.int64),
            "eventnum": arr("eventnum", -1, np.int64),
            "e_p": e_p,
            "e_theta": e_theta,
            "e_phi": e_phi,
            "p_p": p_p,
            "p_theta": p_theta,
            "p_phi": p_phi,
            "g_E": g_E,
            "g_theta": g_theta,
            "g_phi": g_phi,
            "probe_E": probe_E,
            "probe_theta": probe_theta,
            "probe_phi": probe_phi,
            "opportunity": opportunity,
        }
        for name in store:
            store[name].append(np.asarray(values[name])[candidate])
        # endfor
    # endfor

    payload = {}
    for name, parts in store.items():
        if parts:
            payload[name] = np.concatenate(parts)
        else:
            dtype = (
                bool
                if name == "opportunity"
                else np.int64
                if name in {"runnum", "eventnum"}
                else np.float64
            )
            payload[name] = np.asarray([], dtype=dtype)
        # endif
    # endfor
    return EpgRecords(**payload), counts


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
                or energies[0] <= args.minimum_photon_E
                or energies[1] <= args.minimum_photon_E
            ):
                continue
            # endif

            closure = float(
                np.linalg.norm(response @ energies - target)
                / normalization
            )
            if closure > args.closure_max:
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
        "open1", "open2", "trento1", "trento2", "det1", "det2",
    )
    resolved = resolve(path, logical, optional=("runnum", "eventnum"))
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
            "detector1": arr("det1", -1, np.int16),
            "detector2": arr("det2", -1, np.int16),
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


def keys_epg(records: EpgRecords, mode: str, decimals: int) -> np.ndarray:
    if mode == "data":
        return np.column_stack(
            (records.runnum.astype(np.int64), records.eventnum.astype(np.int64))
        )
    # endif
    return np.column_stack(
        [
            np.round(records.e_p, decimals),
            np.round(records.e_theta, decimals),
            np.round(records.e_phi, decimals),
            np.round(records.p_p, decimals),
            np.round(records.p_theta, decimals),
            np.round(records.p_phi, decimals),
        ]
    )


def keys_epgg(records: EpggRecords, mode: str, decimals: int) -> np.ndarray:
    if mode == "data":
        return np.column_stack(
            (records.runnum.astype(np.int64), records.eventnum.astype(np.int64))
        )
    # endif
    return np.column_stack(
        [
            np.round(records.e_p, decimals),
            np.round(records.e_theta, decimals),
            np.round(records.e_phi, decimals),
            np.round(records.p_p, decimals),
            np.round(records.p_theta, decimals),
            np.round(records.p_phi, decimals),
        ]
    )


def structured(matrix: np.ndarray) -> np.ndarray:
    matrix = np.ascontiguousarray(matrix)
    dtype = np.dtype([(f"f{index}", matrix.dtype) for index in range(matrix.shape[1])])
    return matrix.view(dtype).reshape(-1)


def angular_distance_deg(
    theta1: float,
    phi1: float,
    theta2: float,
    phi2: float,
) -> float:
    cosine = (
        math.sin(theta1) * math.sin(theta2) * math.cos(phi1 - phi2)
        + math.cos(theta1) * math.cos(theta2)
    )
    return math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))


def update_reservoir(
    reservoir: List[Dict[str, object]],
    row: Dict[str, object],
    seen: int,
    capacity: int,
    rng: np.random.Generator,
) -> None:
    if capacity <= 0:
        return
    # endif
    if len(reservoir) < capacity:
        reservoir.append(row)
        return
    # endif
    replacement_index = int(rng.integers(0, seen))
    if replacement_index < capacity:
        reservoir[replacement_index] = row
    # endif


def audit_truth_partners(
    epg: EpgRecords,
    epgg: EpggRecords,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[List[Dict[str, object]], Dict[str, object], Dict[str, object]]:
    epg_keys = structured(keys_epg(epg, mode, args.mc_signature_decimals))
    epgg_keys = structured(keys_epgg(epgg, mode, args.mc_signature_decimals))

    epgg_unique, epgg_inverse, epgg_counts = np.unique(
        epgg_keys,
        return_inverse=True,
        return_counts=True,
    )
    epgg_order = np.argsort(epgg_inverse, kind="stable")
    epgg_offsets = np.empty(epgg_counts.size + 1, dtype=np.int64)
    epgg_offsets[0] = 0
    np.cumsum(epgg_counts, out=epgg_offsets[1:])
    epgg_group_by_key = {
        epgg_unique[index].tobytes(): index
        for index in range(epgg_unique.size)
    }

    bins = {
        "tag_angle_deg": np.linspace(0.0, 20.0, 101),
        "tag_relative_E": np.linspace(0.0, 1.0, 101),
        "probe_angle_deg": np.linspace(0.0, 30.0, 121),
        "probe_relative_E": np.linspace(0.0, 3.0, 121),
        "solution_closure": np.linspace(0.0, args.closure_max, 101),
    }
    counts = {
        key: np.zeros(len(edges) - 1, dtype=np.int64)
        for key, edges in bins.items()
    }

    summary = {
        "mode": mode,
        "directed_opportunities": int(np.count_nonzero(epg.opportunity)),
        "opportunities_with_identity_match": 0,
        "opportunities_with_valid_daughter_solution": 0,
        "opportunities_with_tag_identified": 0,
        "opportunities_passing_tag_match_cuts": 0,
        "candidate_epgammagamma_records_tested": 0,
        "candidate_daughter_solutions_tested": 0,
    }

    sample_capacity = max(0, int(args.diagnostic_sample_size))
    rng = np.random.default_rng(
        int(args.diagnostic_seed) + (0 if mode == "data" else 1_000_003)
    )
    sample: List[Dict[str, object]] = []
    accepted_seen = 0

    opportunity_indices = np.flatnonzero(epg.opportunity)
    for opportunity_number, epg_index in enumerate(opportunity_indices, start=1):
        key = epg_keys[epg_index].tobytes()
        group = epgg_group_by_key.get(key)
        if group is None:
            continue
        # endif
        summary["opportunities_with_identity_match"] += 1
        members = epgg_order[epgg_offsets[group]:epgg_offsets[group + 1]]

        best = None
        any_valid = False
        for epgg_index in members:
            summary["candidate_epgammagamma_records_tested"] += 1
            for solution_index in range(epgg.solution_valid.shape[1]):
                if not epgg.solution_valid[epgg_index, solution_index]:
                    continue
                # endif
                any_valid = True
                summary["candidate_daughter_solutions_tested"] += 1

                daughter_values = (
                    (
                        epgg.g1_E[epgg_index, solution_index],
                        epgg.g1_theta[epgg_index, solution_index],
                        epgg.g1_phi[epgg_index, solution_index],
                        epgg.g2_E[epgg_index, solution_index],
                        epgg.g2_theta[epgg_index, solution_index],
                        epgg.g2_phi[epgg_index, solution_index],
                        1,
                    ),
                    (
                        epgg.g2_E[epgg_index, solution_index],
                        epgg.g2_theta[epgg_index, solution_index],
                        epgg.g2_phi[epgg_index, solution_index],
                        epgg.g1_E[epgg_index, solution_index],
                        epgg.g1_theta[epgg_index, solution_index],
                        epgg.g1_phi[epgg_index, solution_index],
                        2,
                    ),
                )

                for (
                    tag_E, tag_theta, tag_phi,
                    partner_E, partner_theta, partner_phi,
                    tag_choice,
                ) in daughter_values:
                    tag_angle = angular_distance_deg(
                        float(epg.g_theta[epg_index]),
                        float(epg.g_phi[epg_index]),
                        float(tag_theta),
                        float(tag_phi),
                    )
                    tag_relative_E = abs(
                        float(epg.g_E[epg_index]) - float(tag_E)
                    ) / max(float(tag_E), 1.0e-12)
                    tag_score = (
                        (tag_angle / max(args.tag_angle_max_deg, 1.0e-12))**2
                        + (
                            tag_relative_E
                            / max(args.tag_relative_E_max, 1.0e-12)
                        )**2
                    )
                    if best is None or tag_score < best[0]:
                        probe_angle = angular_distance_deg(
                            float(epg.probe_theta[epg_index]),
                            float(epg.probe_phi[epg_index]),
                            float(partner_theta),
                            float(partner_phi),
                        )
                        probe_relative_E = abs(
                            float(epg.probe_E[epg_index]) - float(partner_E)
                        ) / max(float(partner_E), 1.0e-12)
                        best = (
                            tag_score,
                            tag_angle,
                            tag_relative_E,
                            probe_angle,
                            probe_relative_E,
                            float(partner_E),
                            float(partner_theta),
                            float(partner_phi),
                            int(epgg_index),
                            int(solution_index),
                            int(tag_choice),
                            float(epgg.solution_closure[epgg_index, solution_index]),
                        )
                    # endif
                # endfor
            # endfor
        # endfor

        if any_valid:
            summary["opportunities_with_valid_daughter_solution"] += 1
        # endif
        if best is None:
            continue
        # endif
        summary["opportunities_with_tag_identified"] += 1

        (
            _,
            tag_angle,
            tag_relative_E,
            probe_angle,
            probe_relative_E,
            partner_E,
            partner_theta,
            partner_phi,
            matched_epgg_index,
            solution_index,
            tag_choice,
            closure,
        ) = best

        passes_tag = (
            tag_angle < args.tag_angle_max_deg
            and tag_relative_E < args.tag_relative_E_max
        )
        if not passes_tag:
            continue
        # endif
        summary["opportunities_passing_tag_match_cuts"] += 1
        accepted_seen += 1

        residual_values = {
            "tag_angle_deg": tag_angle,
            "tag_relative_E": tag_relative_E,
            "probe_angle_deg": probe_angle,
            "probe_relative_E": probe_relative_E,
            "solution_closure": closure,
        }
        for name, value in residual_values.items():
            edges = bins[name]
            bin_index = int(np.searchsorted(edges, value, side="right") - 1)
            if 0 <= bin_index < counts[name].size:
                counts[name][bin_index] += 1
            # endif
        # endfor

        row = {
            "epgamma_opportunity_index": int(epg_index),
            "matched_epgammagamma_index": matched_epgg_index,
            "solution_index": solution_index,
            "tag_choice": tag_choice,
            "tag_E": float(epg.g_E[epg_index]),
            "truth_partner_E": partner_E,
            "predicted_probe_E": float(epg.probe_E[epg_index]),
            "tag_angle_deg": tag_angle,
            "tag_relative_E": tag_relative_E,
            "probe_angle_deg": probe_angle,
            "probe_relative_E": probe_relative_E,
            "solution_closure": closure,
        }
        update_reservoir(
            sample,
            row,
            accepted_seen,
            sample_capacity,
            rng,
        )
    # endfor

    directed = summary["directed_opportunities"]
    for key in (
        "opportunities_with_identity_match",
        "opportunities_with_valid_daughter_solution",
        "opportunities_with_tag_identified",
        "opportunities_passing_tag_match_cuts",
    ):
        summary[f"fraction_{key}"] = (
            summary[key] / directed if directed > 0 else math.nan
        )
    # endfor

    histogram_payload = {
        name: {
            "edges": bins[name].tolist(),
            "counts": counts[name].tolist(),
        }
        for name in bins
    }
    histogram_payload["accepted_entries"] = int(
        summary["opportunities_passing_tag_match_cuts"]
    )
    summary["diagnostic_rows_written"] = len(sample)
    return sample, summary, histogram_payload


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    # endif
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def plot_histograms(
    path: Path,
    title: str,
    payload: Mapping[str, object],
) -> None:
    panels = (
        ("tag_angle_deg", "Tag angular residual (deg)"),
        ("tag_relative_E", "Tag relative energy residual"),
        ("probe_angle_deg", "Predicted-X to truth-partner angle (deg)"),
        ("probe_relative_E", "Predicted-X to truth-partner relative energy"),
        ("solution_closure", "Daughter reconstruction closure"),
    )
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for axis, (name, label) in zip(axes.flat, panels):
        edges = np.asarray(payload[name]["edges"], dtype=float)
        values = np.asarray(payload[name]["counts"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        axis.step(centers, values, where="mid", linewidth=1.3)
        axis.set_xlabel(label)
        axis.set_ylabel("Entries")
        axis.grid(alpha=0.25)
    # endfor
    axes.flat[-1].axis("off")
    fig.suptitle(
        f"{title}\nAccepted tag matches: "
        f"{int(payload.get('accepted_entries', 0)):,}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def process_period(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    period_dir.mkdir(parents=True, exist_ok=True)

    data_epg, data_epg_counts = read_epg(
        period.epg_data,
        period.beam_energy_GeV,
        "data",
        args,
    )
    mc_epg, mc_epg_counts = read_epg(
        period.pi0_epg_mc,
        period.beam_energy_GeV,
        "AAOGEN",
        args,
    )
    data_epgg, data_epgg_counts = read_epgg(
        period.epgg_data,
        period.beam_energy_GeV,
        "data",
        args,
    )
    mc_epgg, mc_epgg_counts = read_epgg(
        period.pi0_epgg_mc,
        period.beam_energy_GeV,
        "AAOGEN",
        args,
    )

    data_sample, data_summary, data_histograms = audit_truth_partners(
        data_epg,
        data_epgg,
        "data",
        args,
    )
    mc_sample, mc_summary, mc_histograms = audit_truth_partners(
        mc_epg,
        mc_epgg,
        "mc",
        args,
    )

    write_csv(period_dir / "data_truth_partner_sample.csv", data_sample)
    write_csv(period_dir / "aaogen_truth_partner_sample.csv", mc_sample)
    plot_histograms(
        period_dir / "data_truth_partner_residuals.png",
        f"{period.label}: data epgammagamma truth partner",
        data_histograms,
    )
    plot_histograms(
        period_dir / "aaogen_truth_partner_residuals.png",
        f"{period.label}: AAOGEN epgammagamma truth partner",
        mc_histograms,
    )

    payload = {
        "period": asdict(period),
        "data_epgamma_read_counts": data_epg_counts,
        "mc_epgamma_read_counts": mc_epg_counts,
        "data_epgammagamma_reconstruction_counts": data_epgg_counts,
        "mc_epgammagamma_reconstruction_counts": mc_epgg_counts,
        "data_truth_partner_summary": data_summary,
        "mc_truth_partner_summary": mc_summary,
        "data_histograms": data_histograms,
        "mc_histograms": mc_histograms,
    }
    with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    # endwith

    log(
        f"{period.label} data: opportunities="
        f"{data_summary['directed_opportunities']:,}, identity match="
        f"{data_summary['fraction_opportunities_with_identity_match']:.3%}, "
        f"valid daughter solution="
        f"{data_summary['fraction_opportunities_with_valid_daughter_solution']:.3%}, "
        f"tag matched="
        f"{data_summary['fraction_opportunities_passing_tag_match_cuts']:.3%}"
    )
    log(
        f"{period.label} AAOGEN: opportunities="
        f"{mc_summary['directed_opportunities']:,}, identity match="
        f"{mc_summary['fraction_opportunities_with_identity_match']:.3%}, "
        f"valid daughter solution="
        f"{mc_summary['fraction_opportunities_with_valid_daughter_solution']:.3%}, "
        f"tag matched="
        f"{mc_summary['fraction_opportunities_passing_tag_match_cuts']:.3%}"
    )
    return period.key, payload


def main() -> int:
    args = parse_args()
    periods = selected_periods(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "script": Path(__file__).name,
        "created_unix_time": time.time(),
        "arguments": vars(args),
        "periods": [],
    }
    for period in periods:
        item = asdict(period)
        for role in ("epg_data", "epgg_data", "pi0_epg_mc", "pi0_epgg_mc"):
            entries, branches = require_tree(getattr(period, role))
            item[f"{role}_entries"] = entries
            item[f"{role}_branches"] = branches
        # endfor
        manifest["periods"].append(item)
    # endfor

    with open(output_dir / "input_manifest.json", "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    workers = max(1, min(int(args.workers), MAX_WORKERS, len(periods)))
    log(
        f"EPGAMMAGAMMA TRUTH-PARTNER AUDIT: {len(periods)} period(s), "
        f"{workers} worker(s). Only E_tag >= {args.tag_E_min:g} GeV and "
        f"predicted E_X >= {args.probe_E_min:g} GeV are required."
    )

    metadata: Dict[str, object] = {}
    if workers == 1:
        for period in periods:
            key, payload = process_period(period, vars(args))
            metadata[key] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers
        ) as executor:
            futures = [
                executor.submit(process_period, period, vars(args))
                for period in periods
            ]
            for future in concurrent.futures.as_completed(futures):
                key, payload = future.result()
                metadata[key] = payload
                log(f"Completed {key}.")
            # endfor
        # endwith
    # endif

    rows = []
    for period in periods:
        payload = metadata[period.key]
        for sample, key in (
            ("Data", "data_truth_partner_summary"),
            ("AAOGEN MC", "mc_truth_partner_summary"),
        ):
            rows.append(
                {
                    "period": period.key,
                    "period_label": period.label,
                    "sample": sample,
                    **payload[key],
                }
            )
        # endfor
    # endfor

    write_csv(output_dir / "truth_partner_audit_summary.csv", rows)
    with open(
        output_dir / "truth_partner_audit_summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(
            {"rows": rows, "period_metadata": metadata},
            handle,
            indent=2,
        )
    # endwith

    log(f"Wrote truth-partner audit to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
