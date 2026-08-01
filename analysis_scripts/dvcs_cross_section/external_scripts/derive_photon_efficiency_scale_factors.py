#!/usr/bin/env python3
"""
photon_efficiency_step1_expected_probe_audit_v1.py

Step 1 of the restarted RGA photon-efficiency study.

Purpose
-------
Construct and validate the ep -> ep gamma X expected-probe population from the
three loose-threshold epgamma samples:

  1. real data;
  2. DVCSGEN/BH Monte Carlo;
  3. AAOGEN exclusive-pi0 Monte Carlo processed as epgamma.

For every selected epgamma entry, the measured photon is treated as the tag and
the missing four-vector

    X = k + p_target - k' - p' - gamma_tag

is reconstructed directly from four-momentum conservation.

This first step deliberately performs NO template fit, NO epgamma/epgammagamma
matching, NO efficiency calculation, and NO pi0-fraction extraction.  Its only
job is to verify that all three input samples produce sensible and comparably
defined expected-probe populations before the next stage introduces exclusive-
pi0 isolation.

Nominal directed opportunity requirements
-----------------------------------------
  E_tag >= 0.40 GeV
  E_probe >= 2.00 GeV

There is no E_tag < 2 GeV requirement.  Events in which both photons would have
energies above 2 GeV remain eligible.

Outputs
-------
output/photon_efficiency_study/step1_expected_probe_audit/

  input_manifest.json
  all_periods_cutflow.csv
  all_periods_expected_probe_counts.csv
  <period>/
      branch_mapping.json
      cutflow.json
      expected_probe_records.npz
      expected_probe_kinematics.png
      photon_hypothesis_diagnostics.png
      topology_breakdown.png

The NPZ files contain only the selected expected-probe records needed for
validation and later development.  They are not production corrections.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import uproot
except ImportError as exc:
    raise SystemExit(
        "This script requires numpy, matplotlib, and uproot.\n"
        "Install with: python3 -m pip install numpy matplotlib uproot"
    ) from exc
#endif


TREE_NAME = "PhysicsEvents"
OUTPUT_DIR = Path("output/photon_efficiency_study/step1_expected_probe_audit")
MAX_WORKERS = 8

PROTON_MASS_GEV = 0.9382720813
ELECTRON_MASS_GEV = 0.00051099895


@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy_GeV: float
    data: str
    dvcs_mc: str
    pi0_mc: str


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        "fa18_inb",
        "Fa18 Inb",
        10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
    ),
    PeriodConfig(
        "fa18_out",
        "Fa18 Out",
        10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/bkg_rga_fa18_out_epgamma_0.40GeV.root",
    ),
    PeriodConfig(
        "sp19_inb",
        "Sp19 Inb",
        10.200,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
    ),
    PeriodConfig(
        "sp18_inb",
        "Sp18 Inb",
        10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
    ),
    PeriodConfig(
        "sp18_out",
        "Sp18 Out",
        10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/"
        "dvcsgen_files_greater_than_0.35GeV/bkg_rga_sp18_out_epgamma_0.40GeV.root",
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
    "g_detector": ("detector2", "g1_detector", "gamma1_detector"),
    "p_detector": ("detector1", "p_detector", "proton_detector"),
    "fiducial_status": ("fiducial_status",),
    "Q2": ("Q2", "q2"),
    "W": ("W", "w"),
    "y": ("y", "inelasticity"),
    "z": ("z",),
    "t1": ("t1", "t", "minus_t"),
    "Delta_phi": ("Delta_phi", "delta_phi", "dphi"),
    "theta_gamma_gamma": ("theta_gamma_gamma", "theta_pi0_pi0"),
    "pTmiss": ("pTmiss", "ptmiss", "pT_miss"),
    "Emiss2": ("Emiss2", "E_miss2", "Emiss_sq"),
    "Mx2": ("Mx2", "Mx2_epg", "Mx2_epgamma"),
    "Mx2_2": ("Mx2_2", "Mx2_egamma", "Mx2_gamma"),
}


REQUIRED_BRANCHES: Tuple[str, ...] = (
    "e_p", "e_theta", "e_phi",
    "p_p", "p_theta", "p_phi",
    "g_E", "g_theta", "g_phi",
)

OPTIONAL_BRANCHES: Tuple[str, ...] = (
    "runnum", "eventnum", "g_detector", "p_detector", "fiducial_status",
    "Q2", "W", "y", "z", "t1",
    "Delta_phi", "theta_gamma_gamma", "pTmiss", "Emiss2", "Mx2", "Mx2_2",
)


@dataclass
class SampleResult:
    sample: str
    total_entries: int
    finite_four_vectors: int
    common_global_cuts: int
    tag_energy: int
    predicted_probe_finite: int
    predicted_probe_energy: int
    predicted_probe_acceptance: int
    selected_expected_probes: int
    selected_ft: int
    selected_fd: int
    selected_outside_acceptance: int


def log(message: str) -> None:
    timestamp = time.strftime("%H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Step 1: construct and validate epgammaX expected-probe populations "
            "for data, DVCSGEN, and AAOGEN-as-epgamma."
        )
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[period.key for period in PERIODS],
        help="Run only selected period(s). May be repeated. Default: all periods.",
    )
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--step-size", default="200 MB")
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--output-dir", default=str(OUTPUT_DIR))

    parser.add_argument("--tag-E-min", type=float, default=0.40)
    parser.add_argument("--probe-E-min", type=float, default=2.00)
    parser.add_argument("--probe-E-max", type=float, default=9.50)

    parser.add_argument("--ft-theta-min", type=float, default=2.5)
    parser.add_argument("--ft-theta-max", type=float, default=5.0)
    parser.add_argument("--fd-theta-min", type=float, default=5.0)
    parser.add_argument("--fd-theta-max", type=float, default=35.0)

    parser.add_argument("--Q2-min", type=float, default=1.0)
    parser.add_argument("--W-min", type=float, default=2.0)
    parser.add_argument("--y-max", type=float, default=0.8)
    parser.add_argument("--z-min", type=float, default=0.65)
    parser.add_argument("--minus-t-max", type=float, default=1.0)

    parser.add_argument(
        "--disable-global-cuts",
        action="store_true",
        help="Disable Q2, W, y, z, and -t cuts for a diagnostic comparison.",
    )
    parser.add_argument(
        "--require-fiducial-status-111",
        action="store_true",
        help="Optionally require fiducial_status == 111 when the branch exists.",
    )
    return parser.parse_args()


def selected_periods(args: argparse.Namespace) -> List[PeriodConfig]:
    if not args.period:
        return list(PERIODS)
    # endif
    wanted = set(args.period)
    return [period for period in PERIODS if period.key in wanted]


def require_tree(path: str) -> Tuple[int, List[str]]:
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"Input ROOT file does not exist: {path}")
    # endif

    with uproot.open(path) as root_file:
        if TREE_NAME not in root_file:
            raise KeyError(f"Tree '{TREE_NAME}' does not exist in {path}")
        # endif
        tree = root_file[TREE_NAME]
        return int(tree.num_entries), [str(key) for key in tree.keys()]
    # endwith


def resolve_branches(path: str) -> Dict[str, Optional[str]]:
    _, keys = require_tree(path)
    available = set(keys)
    resolved: Dict[str, Optional[str]] = {}
    missing: List[str] = []

    for logical in REQUIRED_BRANCHES + OPTIONAL_BRANCHES:
        branch = next(
            (candidate for candidate in ALIASES[logical] if candidate in available),
            None,
        )
        resolved[logical] = branch
        if logical in REQUIRED_BRANCHES and branch is None:
            missing.append(
                f"{logical}: tried {', '.join(ALIASES[logical])}"
            )
        # endif
    # endfor

    if missing:
        raise KeyError(
            f"Missing required branches in {path}:\n  "
            + "\n  ".join(missing)
        )
    # endif

    return resolved


def finite_array(
    arrays: Mapping[str, np.ndarray],
    branch: Optional[str],
    default: float = math.nan,
) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    if branch is None:
        return np.full(n, default, dtype=np.float64)
    # endif
    return np.asarray(arrays[branch], dtype=np.float64)


def spherical_to_cartesian(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sin_theta = np.sin(theta)
    px = momentum * sin_theta * np.cos(phi)
    py = momentum * sin_theta * np.sin(phi)
    pz = momentum * np.cos(theta)
    return px, py, pz


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


def reconstruct_missing_probe(
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
    g_E, g_px, g_py, g_pz = photon_four_vector(
        tag_E, tag_theta, tag_phi
    )

    probe_E = beam_energy + PROTON_MASS_GEV - e_E - p_E - g_E
    probe_px = -e_px - p_px - g_px
    probe_py = -e_py - p_py - g_py
    probe_pz = beam_energy - e_pz - p_pz - g_pz

    probe_p = np.sqrt(
        np.maximum(
            probe_px * probe_px + probe_py * probe_py + probe_pz * probe_pz,
            0.0,
        )
    )
    probe_pt = np.sqrt(
        np.maximum(probe_px * probe_px + probe_py * probe_py, 0.0)
    )
    probe_theta = np.arctan2(probe_pt, probe_pz)
    probe_phi = np.arctan2(probe_py, probe_px)
    probe_m2 = probe_E * probe_E - probe_p * probe_p
    probe_E_minus_p = probe_E - probe_p

    return {
        "probe_E": probe_E,
        "probe_px": probe_px,
        "probe_py": probe_py,
        "probe_pz": probe_pz,
        "probe_p": probe_p,
        "probe_pt": probe_pt,
        "probe_theta": probe_theta,
        "probe_phi": probe_phi,
        "probe_m2": probe_m2,
        "probe_E_minus_p": probe_E_minus_p,
    }


def classify_probe(
    theta_deg: np.ndarray,
    phi_rad: np.ndarray,
    args: argparse.Namespace,
) -> Tuple[np.ndarray, np.ndarray]:
    detector = np.full(theta_deg.shape, -1, dtype=np.int8)
    sector = np.zeros(theta_deg.shape, dtype=np.int8)

    ft = (
        np.isfinite(theta_deg)
        & (theta_deg >= args.ft_theta_min)
        & (theta_deg < args.ft_theta_max)
    )
    fd = (
        np.isfinite(theta_deg)
        & (theta_deg >= args.fd_theta_min)
        & (theta_deg < args.fd_theta_max)
    )

    detector[ft] = 0
    detector[fd] = 1

    phi_deg = np.degrees(phi_rad)
    wrapped = (phi_deg + 30.0) % 360.0
    sector[fd] = (np.floor(wrapped[fd] / 60.0).astype(np.int8) + 1)

    return detector, sector


def common_global_mask(
    arrays: Mapping[str, np.ndarray],
    resolved: Mapping[str, Optional[str]],
    args: argparse.Namespace,
) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    mask = np.ones(n, dtype=bool)

    if args.disable_global_cuts:
        return mask
    # endif

    def apply_lower(logical: str, threshold: float) -> None:
        nonlocal mask
        branch = resolved[logical]
        if branch is not None:
            values = finite_array(arrays, branch)
            mask &= np.isfinite(values) & (values > threshold)
        # endif

    def apply_upper(logical: str, threshold: float) -> None:
        nonlocal mask
        branch = resolved[logical]
        if branch is not None:
            values = finite_array(arrays, branch)
            mask &= np.isfinite(values) & (values < threshold)
        # endif

    apply_lower("Q2", args.Q2_min)
    apply_lower("W", args.W_min)
    apply_upper("y", args.y_max)
    apply_lower("z", args.z_min)

    if resolved["t1"] is not None:
        t1 = finite_array(arrays, resolved["t1"])
        mask &= np.isfinite(t1) & ((-t1) < args.minus_t_max)
    # endif

    if args.require_fiducial_status_111 and resolved["fiducial_status"] is not None:
        fiducial = finite_array(arrays, resolved["fiducial_status"])
        mask &= np.isfinite(fiducial) & (fiducial == 111)
    # endif

    return mask


def append_selected(
    storage: Dict[str, List[np.ndarray]],
    name: str,
    values: np.ndarray,
    mask: np.ndarray,
) -> None:
    storage.setdefault(name, []).append(
        np.asarray(values[mask], dtype=np.float32)
    )


def concatenate_storage(
    storage: Mapping[str, List[np.ndarray]],
) -> Dict[str, np.ndarray]:
    result: Dict[str, np.ndarray] = {}
    for name, chunks in storage.items():
        if chunks:
            result[name] = np.concatenate(chunks)
        else:
            result[name] = np.empty(0, dtype=np.float32)
        # endif
    # endfor
    return result


def process_sample(
    path: str,
    sample: str,
    beam_energy: float,
    args: argparse.Namespace,
) -> Tuple[SampleResult, Dict[str, np.ndarray], Dict[str, Optional[str]]]:
    entries, _ = require_tree(path)
    resolved = resolve_branches(path)

    expressions = sorted(
        {branch for branch in resolved.values() if branch is not None}
    )

    counts = {
        "total_entries": 0,
        "finite_four_vectors": 0,
        "common_global_cuts": 0,
        "tag_energy": 0,
        "predicted_probe_finite": 0,
        "predicted_probe_energy": 0,
        "predicted_probe_acceptance": 0,
        "selected_expected_probes": 0,
        "selected_ft": 0,
        "selected_fd": 0,
        "selected_outside_acceptance": 0,
    }
    storage: Dict[str, List[np.ndarray]] = {}

    log(f"Reading {sample}: {path}")

    seen = 0
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
            keep = args.max_events - seen
            arrays = {key: value[:keep] for key, value in arrays.items()}
            n = keep
        # endif
        seen += n
        counts["total_entries"] += n

        e_p = finite_array(arrays, resolved["e_p"])
        e_theta = finite_array(arrays, resolved["e_theta"])
        e_phi = finite_array(arrays, resolved["e_phi"])
        p_p = finite_array(arrays, resolved["p_p"])
        p_theta = finite_array(arrays, resolved["p_theta"])
        p_phi = finite_array(arrays, resolved["p_phi"])
        tag_E = finite_array(arrays, resolved["g_E"])
        tag_theta = finite_array(arrays, resolved["g_theta"])
        tag_phi = finite_array(arrays, resolved["g_phi"])

        finite_four_vectors = (
            np.isfinite(e_p)
            & np.isfinite(e_theta)
            & np.isfinite(e_phi)
            & np.isfinite(p_p)
            & np.isfinite(p_theta)
            & np.isfinite(p_phi)
            & np.isfinite(tag_E)
            & np.isfinite(tag_theta)
            & np.isfinite(tag_phi)
            & (e_p > 0.0)
            & (p_p > 0.0)
            & (tag_E > 0.0)
        )
        counts["finite_four_vectors"] += int(
            np.count_nonzero(finite_four_vectors)
        )

        global_mask = finite_four_vectors & common_global_mask(
            arrays, resolved, args
        )
        counts["common_global_cuts"] += int(np.count_nonzero(global_mask))

        tag_mask = global_mask & (tag_E >= args.tag_E_min)
        counts["tag_energy"] += int(np.count_nonzero(tag_mask))

        probe = reconstruct_missing_probe(
            beam_energy,
            e_p,
            e_theta,
            e_phi,
            p_p,
            p_theta,
            p_phi,
            tag_E,
            tag_theta,
            tag_phi,
        )

        probe_finite = (
            tag_mask
            & np.isfinite(probe["probe_E"])
            & np.isfinite(probe["probe_p"])
            & np.isfinite(probe["probe_theta"])
            & np.isfinite(probe["probe_phi"])
            & np.isfinite(probe["probe_m2"])
            & (probe["probe_E"] > 0.0)
            & (probe["probe_p"] > 0.0)
        )
        counts["predicted_probe_finite"] += int(
            np.count_nonzero(probe_finite)
        )

        probe_energy = (
            probe_finite
            & (probe["probe_E"] >= args.probe_E_min)
            & (probe["probe_E"] < args.probe_E_max)
        )
        counts["predicted_probe_energy"] += int(
            np.count_nonzero(probe_energy)
        )

        probe_theta_deg = np.degrees(probe["probe_theta"])
        detector, sector = classify_probe(
            probe_theta_deg, probe["probe_phi"], args
        )
        accepted = probe_energy & (detector >= 0)
        outside = probe_energy & (detector < 0)

        counts["predicted_probe_acceptance"] += int(
            np.count_nonzero(accepted)
        )
        counts["selected_expected_probes"] += int(
            np.count_nonzero(accepted)
        )
        counts["selected_ft"] += int(
            np.count_nonzero(accepted & (detector == 0))
        )
        counts["selected_fd"] += int(
            np.count_nonzero(accepted & (detector == 1))
        )
        counts["selected_outside_acceptance"] += int(
            np.count_nonzero(outside)
        )

        append_selected(storage, "tag_E", tag_E, accepted)
        append_selected(storage, "tag_theta_deg", np.degrees(tag_theta), accepted)
        append_selected(storage, "tag_phi_rad", tag_phi, accepted)
        append_selected(storage, "probe_E", probe["probe_E"], accepted)
        append_selected(storage, "probe_p", probe["probe_p"], accepted)
        append_selected(storage, "probe_theta_deg", probe_theta_deg, accepted)
        append_selected(storage, "probe_phi_rad", probe["probe_phi"], accepted)
        append_selected(storage, "probe_pt", probe["probe_pt"], accepted)
        append_selected(storage, "probe_m2", probe["probe_m2"], accepted)
        append_selected(
            storage, "probe_E_minus_p", probe["probe_E_minus_p"], accepted
        )
        append_selected(storage, "detector", detector, accepted)
        append_selected(storage, "sector", sector, accepted)

        for logical in (
            "Delta_phi", "theta_gamma_gamma", "pTmiss",
            "Emiss2", "Mx2", "Mx2_2", "Q2", "W", "y", "z", "t1",
        ):
            values = finite_array(arrays, resolved[logical])
            append_selected(storage, logical, values, accepted)
        # endfor
    # endfor

    records = concatenate_storage(storage)
    result = SampleResult(sample=sample, **counts)

    if result.total_entries != min(entries, args.max_events or entries):
        log(
            f"WARNING: {sample} expected {entries:,} entries but processed "
            f"{result.total_entries:,}."
        )
    # endif

    return result, records, resolved


def normalized_hist(
    axis: plt.Axes,
    values: np.ndarray,
    bins: np.ndarray,
    label: str,
) -> None:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return
    # endif
    weights = np.ones(finite.size, dtype=float) / float(finite.size)
    axis.hist(
        finite,
        bins=bins,
        weights=weights,
        histtype="step",
        linewidth=1.4,
        label=f"{label} ({finite.size:,})",
    )


def plot_expected_probe_kinematics(
    path: Path,
    period_label: str,
    records_by_sample: Mapping[str, Mapping[str, np.ndarray]],
) -> None:
    variables = (
        ("tag_E", np.linspace(0.4, 9.5, 92), r"Tag photon energy (GeV)"),
        ("probe_E", np.linspace(2.0, 9.5, 76), r"Predicted probe energy (GeV)"),
        ("tag_theta_deg", np.linspace(0.0, 40.0, 81), r"Tag photon $\theta$ (deg)"),
        ("probe_theta_deg", np.linspace(0.0, 40.0, 81), r"Predicted probe $\theta$ (deg)"),
        ("tag_phi_rad", np.linspace(-math.pi, math.pi, 73), r"Tag photon $\phi$ (rad)"),
        ("probe_phi_rad", np.linspace(-math.pi, math.pi, 73), r"Predicted probe $\phi$ (rad)"),
    )

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for axis, (name, bins, xlabel) in zip(axes.flat, variables):
        for sample, records in records_by_sample.items():
            normalized_hist(axis, records[name], bins, sample)
        # endfor
        axis.set_xlabel(xlabel)
        axis.set_ylabel("Fraction / bin")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    # endfor

    fig.suptitle(
        f"{period_label}: Step 1 expected-probe kinematics\n"
        r"$E_{\rm tag}\geq0.4$ GeV, "
        r"$E_{\rm probe}^{\rm pred}\geq2$ GeV, FT/FD angular acceptance"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_photon_hypothesis(
    path: Path,
    period_label: str,
    records_by_sample: Mapping[str, Mapping[str, np.ndarray]],
) -> None:
    variables = (
        ("probe_E_minus_p", np.linspace(-1.0, 1.0, 121), r"$E_X-|\vec p_X|$ (GeV)"),
        ("probe_m2", np.linspace(-2.0, 2.0, 121), r"$m_X^2$ (GeV$^2$)"),
        ("probe_pt", np.linspace(0.0, 1.5, 101), r"Predicted probe $p_T$ (GeV)"),
        ("Delta_phi", np.linspace(2.75, 3.50, 101), r"$\Delta\phi$ (rad)"),
        ("theta_gamma_gamma", np.linspace(0.0, 3.1, 101), r"$\theta_{\gamma\gamma}$ (rad)"),
        ("pTmiss", np.linspace(0.0, 0.6, 101), r"$p_T^{\rm miss}$ (GeV)"),
        ("Emiss2", np.linspace(-1.0, 2.0, 101), r"$E_{\rm miss}$ diagnostic"),
        ("Mx2_2", np.linspace(-1.0, 4.0, 101), r"$M_{x2}^2$ (GeV$^2$)"),
    )

    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    for axis, (name, bins, xlabel) in zip(axes.flat, variables):
        for sample, records in records_by_sample.items():
            normalized_hist(axis, records[name], bins, sample)
        # endfor
        axis.set_xlabel(xlabel)
        axis.set_ylabel("Fraction / bin")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    # endfor

    fig.suptitle(
        f"{period_label}: Step 1 photon-hypothesis and exclusivity diagnostics\n"
        "No cuts on the displayed variables"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_topology_breakdown(
    path: Path,
    period_label: str,
    records_by_sample: Mapping[str, Mapping[str, np.ndarray]],
) -> None:
    labels = ["FT", "FD S1", "FD S2", "FD S3", "FD S4", "FD S5", "FD S6"]
    x = np.arange(len(labels), dtype=float)
    width = 0.24

    fig, axis = plt.subplots(figsize=(12, 6))
    for offset_index, (sample, records) in enumerate(records_by_sample.items()):
        detector = records["detector"].astype(int)
        sector = records["sector"].astype(int)
        counts = [int(np.count_nonzero(detector == 0))]
        counts.extend(
            int(np.count_nonzero((detector == 1) & (sector == value)))
            for value in range(1, 7)
        )
        # endfor
        total = max(sum(counts), 1)
        fractions = np.asarray(counts, dtype=float) / float(total)
        offset = (offset_index - 1) * width
        axis.bar(x + offset, fractions, width=width, label=sample)
    # endfor

    axis.set_xticks(x, labels)
    axis.set_ylabel("Fraction of selected expected probes")
    axis.set_title(f"{period_label}: predicted-probe detector categories")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def process_period(
    period: PeriodConfig,
    args_dict: Mapping[str, object],
) -> Tuple[str, Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    period_dir.mkdir(parents=True, exist_ok=True)

    sample_paths = {
        "Data": period.data,
        "DVCSGEN": period.dvcs_mc,
        "AAOGEN pi0": period.pi0_mc,
    }

    results: Dict[str, SampleResult] = {}
    records_by_sample: Dict[str, Dict[str, np.ndarray]] = {}
    mapping: Dict[str, Dict[str, Optional[str]]] = {}

    for sample, path in sample_paths.items():
        result, records, resolved = process_sample(
            path, sample, period.beam_energy_GeV, args
        )
        results[sample] = result
        records_by_sample[sample] = records
        mapping[sample] = resolved
    # endfor

    with open(period_dir / "branch_mapping.json", "w", encoding="utf-8") as handle:
        json.dump(mapping, handle, indent=2)
    # endwith

    with open(period_dir / "cutflow.json", "w", encoding="utf-8") as handle:
        json.dump(
            {sample: asdict(result) for sample, result in results.items()},
            handle,
            indent=2,
        )
    # endwith

    # Save one compact NPZ per period with sample-prefixed arrays.
    npz_payload: Dict[str, np.ndarray] = {}
    for sample, records in records_by_sample.items():
        token = sample.lower().replace(" ", "_")
        for name, values in records.items():
            npz_payload[f"{token}__{name}"] = values
        # endfor
    # endfor
    np.savez_compressed(
        period_dir / "expected_probe_records.npz",
        **npz_payload,
    )

    plot_expected_probe_kinematics(
        period_dir / "expected_probe_kinematics.png",
        period.label,
        records_by_sample,
    )
    plot_photon_hypothesis(
        period_dir / "photon_hypothesis_diagnostics.png",
        period.label,
        records_by_sample,
    )
    plot_topology_breakdown(
        period_dir / "topology_breakdown.png",
        period.label,
        records_by_sample,
    )

    payload = {
        "period": asdict(period),
        "results": {
            sample: asdict(result) for sample, result in results.items()
        },
    }
    return period.key, payload


def write_summary_csvs(
    output_dir: Path,
    payload_by_period: Mapping[str, Mapping[str, object]],
) -> None:
    cutflow_fields = [
        "period", "sample", "total_entries", "finite_four_vectors",
        "common_global_cuts", "tag_energy", "predicted_probe_finite",
        "predicted_probe_energy", "predicted_probe_acceptance",
        "selected_expected_probes",
    ]
    count_fields = [
        "period", "sample", "selected_expected_probes",
        "selected_ft", "selected_fd", "selected_outside_acceptance",
    ]

    with open(
        output_dir / "all_periods_cutflow.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=cutflow_fields)
        writer.writeheader()
        for period_key, payload in payload_by_period.items():
            for sample, result in payload["results"].items():
                row = {"period": period_key, "sample": sample}
                row.update({field: result[field] for field in cutflow_fields[2:]})
                writer.writerow(row)
            # endfor
        # endfor
    # endwith

    with open(
        output_dir / "all_periods_expected_probe_counts.csv",
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=count_fields)
        writer.writeheader()
        for period_key, payload in payload_by_period.items():
            for sample, result in payload["results"].items():
                row = {"period": period_key, "sample": sample}
                row.update({field: result[field] for field in count_fields[2:]})
                writer.writerow(row)
            # endfor
        # endfor
    # endwith


def main() -> int:
    args = parse_args()
    periods = selected_periods(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    workers = max(1, min(int(args.workers), MAX_WORKERS, len(periods)))
    log(
        f"Step 1 only: constructing epgammaX expected-probe populations for "
        f"{len(periods)} period(s) with {workers} worker(s)."
    )
    log(
        f"Directed opportunity: E_tag >= {args.tag_E_min:g} GeV; "
        f"{args.probe_E_min:g} <= E_probe < {args.probe_E_max:g} GeV."
    )
    log("No template fit, event matching, or efficiency extraction will run.")

    manifest = {
        "script": Path(__file__).name,
        "created_unix_time": time.time(),
        "step": 1,
        "purpose": "Construct and validate epgammaX expected-probe populations.",
        "arguments": vars(args),
        "periods": [asdict(period) for period in periods],
    }
    with open(
        output_dir / "input_manifest.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(manifest, handle, indent=2)
    # endwith

    payload_by_period: Dict[str, Dict[str, object]] = {}
    args_dict = vars(args)

    if workers == 1:
        for period in periods:
            key, payload = process_period(period, args_dict)
            payload_by_period[key] = payload
        # endfor
    else:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(process_period, period, args_dict): period.key
                for period in periods
            }
            for future in as_completed(futures):
                key, payload = future.result()
                payload_by_period[key] = payload
                log(f"Completed {key}.")
            # endfor
        # endwith
    # endif

    ordered_payload = {
        period.key: payload_by_period[period.key] for period in periods
    }
    write_summary_csvs(output_dir, ordered_payload)

    with open(
        output_dir / "summary.json",
        "w",
        encoding="utf-8",
    ) as handle:
        json.dump(ordered_payload, handle, indent=2)
    # endwith

    log(f"Wrote Step-1 outputs to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
