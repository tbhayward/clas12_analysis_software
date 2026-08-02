#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v7_compact_matching_audit.py

Temporary no-cuts matching audit for the newly reprocessed multi-photon epgamma
trees.

This release deliberately does NOT calculate photon efficiencies, template
fractions, or scale factors.  It is designed to answer one question first:

    When an epgamma entry predicts a missing photon with E_X >= 2 GeV, is there
    another reconstructed photon entry in the same underlying event whose
    measured four-vector is compatible with X?

Opportunity definition
----------------------
Only the following requirements are applied:

    * finite electron, proton, and tag-photon kinematics;
    * E_tag >= 0.40 GeV;
    * finite missing four-vector;
    * E_X >= 2.00 GeV.

No Q2, W, y, z, t, fiducial, opening-angle, missing-mass, E-|p|,
detector-acceptance, topology, pi0-mass, or epgammagamma-catalog requirement is
applied.

Partner pool
------------
Every finite reconstructed photon entry with E_gamma >= 0.40 GeV is retained.
The current tag entry is excluded.  Among the remaining photons in the same
underlying event, the script records the nearest candidate in combined angular
and relative-energy residual.

Data identity:

    (runnum, evnum)

AAOGEN identity:

    rounded reconstructed electron/proton signature
    (e_p, e_theta, e_phi, p_p, p_theta, p_phi)

The epgammagamma trees are read only as identity catalogs for diagnostics.  They
do not gate the matching.

Outputs
-------
output/photon_efficiency_study/matching_audit/

    matching_audit_summary.csv
    matching_audit_summary.json
    input_manifest.json
    <period>/
        data_matching_residuals.png
        aaogen_matching_residuals.png
        data_nearest_match_sample.csv
        aaogen_nearest_match_sample.csv
        metadata.json

This is intentionally an audit release.  It never writes one row per
opportunity.  Exact counters and residual histograms are accumulated in memory,
while only a small deterministic reservoir sample is written for manual
inspection.  Once event-level photon matching is validated, the physics cuts
and BH/DVCS + pi0 template extraction can be restored in a later production
release.
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
DEFAULT_OUTPUT_DIR = "output/photon_efficiency_study/matching_audit"
DEFAULT_STEP_SIZE = "200 MB"
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
class IdentityCatalog:
    keys: np.ndarray
    entries: int
    unique_keys: int


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="No-cuts event-level partner-photon matching audit."
    )
    parser.add_argument("--period", action="append", choices=[p.key for p in PERIODS])
    parser.add_argument("--workers", type=int, default=5)
    parser.add_argument("--step-size", default=DEFAULT_STEP_SIZE)
    parser.add_argument("--max-events", type=int, default=None)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--tag-E-min", type=float, default=0.40)
    parser.add_argument("--probe-E-min", type=float, default=2.00)
    parser.add_argument("--mc-signature-decimals", type=int, default=10)
    parser.add_argument("--score-angle-scale-deg", type=float, default=3.0)
    parser.add_argument("--score-relative-E-scale", type=float, default=0.35)
    parser.add_argument(
        "--diagnostic-sample-size",
        type=int,
        default=5000,
        help=(
            "Maximum number of nearest-match rows written per sample and "
            "period. Exact summary counters and histograms still use every "
            "opportunity."
        ),
    )
    parser.add_argument(
        "--diagnostic-seed",
        type=int,
        default=20260801,
        help="Base seed for deterministic reservoir sampling.",
    )
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


def resolve(path: str, logical: Sequence[str], optional: Sequence[str] = ()) -> Dict[str, Optional[str]]:
    _, keys = require_tree(path)
    available = set(keys)
    optional_set = set(optional)
    result: Dict[str, Optional[str]] = {}
    missing = []
    for name in logical:
        branch = next((candidate for candidate in ALIASES.get(name, (name,)) if candidate in available), None)
        result[name] = branch
        if branch is None and name not in optional_set:
            missing.append(name)
        # endif
    # endfor
    if missing:
        raise KeyError(f"Missing branches in {path}: {missing}")
    # endif
    return result


def finite_array(arrays: Mapping[str, np.ndarray], branch: Optional[str], default=math.nan, dtype=np.float64) -> np.ndarray:
    n = len(next(iter(arrays.values())))
    if branch is None:
        return np.full(n, default, dtype=dtype)
    # endif
    return np.asarray(arrays[branch], dtype=dtype)


def spherical(momentum: np.ndarray, theta: np.ndarray, phi: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    st = np.sin(theta)
    return momentum * st * np.cos(phi), momentum * st * np.sin(phi), momentum * np.cos(theta)


def massive(momentum: np.ndarray, theta: np.ndarray, phi: np.ndarray, mass: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    px, py, pz = spherical(momentum, theta, phi)
    return np.sqrt(momentum**2 + mass**2), px, py, pz


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
    e_E, e_px, e_py, e_pz = massive(e_p, e_theta, e_phi, ELECTRON_MASS_GEV)
    p_E, p_px, p_py, p_pz = massive(p_p, p_theta, p_phi, PROTON_MASS_GEV)
    g_px, g_py, g_pz = spherical(g_E, g_theta, g_phi)

    E = beam_energy + PROTON_MASS_GEV - e_E - p_E - g_E
    px = -e_px - p_px - g_px
    py = -e_py - p_py - g_py
    pz = beam_energy - e_pz - p_pz - g_pz
    pt = np.sqrt(np.maximum(px**2 + py**2, 0.0))
    theta = np.arctan2(pt, pz)
    phi = np.mod(np.arctan2(py, px), 2.0 * math.pi)
    return E, theta, phi


def read_epg(path: str, beam_energy: float, mode: str, args: argparse.Namespace) -> Tuple[EpgRecords, Dict[str, int]]:
    logical = ("runnum", "eventnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi", "g_E", "g_theta", "g_phi")
    resolved = resolve(path, logical, optional=("runnum", "eventnum"))
    expressions = sorted({branch for branch in resolved.values() if branch is not None})

    store: Dict[str, List[np.ndarray]] = {name: [] for name in EpgRecords.__dataclass_fields__}
    counts = {
        "tree_entries": 0,
        "finite_epgamma_entries": 0,
        "photon_candidates_E_above_0p4": 0,
        "directed_opportunities_probe_E_above_2": 0,
    }

    seen = 0
    log(f"Reading {mode} epgamma records from {path}")
    for arrays in uproot.iterate(f"{path}:{TREE_NAME}", expressions=expressions, step_size=args.step_size, library="np"):
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

        def arr(name: str, default=math.nan, dtype=np.float64) -> np.ndarray:
            return finite_array(arrays, resolved.get(name), default, dtype)

        e_p, e_theta, e_phi = arr("e_p"), arr("e_theta"), np.mod(arr("e_phi"), 2.0 * math.pi)
        p_p, p_theta, p_phi = arr("p_p"), arr("p_theta"), np.mod(arr("p_phi"), 2.0 * math.pi)
        g_E, g_theta, g_phi = arr("g_E"), arr("g_theta"), np.mod(arr("g_phi"), 2.0 * math.pi)

        finite = (
            np.isfinite(e_p) & np.isfinite(e_theta) & np.isfinite(e_phi)
            & np.isfinite(p_p) & np.isfinite(p_theta) & np.isfinite(p_phi)
            & np.isfinite(g_E) & np.isfinite(g_theta) & np.isfinite(g_phi)
            & (e_p > 0.0) & (p_p > 0.0) & (g_E > 0.0)
        )
        counts["finite_epgamma_entries"] += int(np.count_nonzero(finite))
        candidate = finite & (g_E >= args.tag_E_min)
        counts["photon_candidates_E_above_0p4"] += int(np.count_nonzero(candidate))

        probe_E, probe_theta, probe_phi = reconstruct_probe(
            beam_energy, e_p, e_theta, e_phi, p_p, p_theta, p_phi, g_E, g_theta, g_phi
        )
        opportunity = (
            candidate
            & np.isfinite(probe_E)
            & np.isfinite(probe_theta)
            & np.isfinite(probe_phi)
            & (probe_E >= args.probe_E_min)
        )
        counts["directed_opportunities_probe_E_above_2"] += int(np.count_nonzero(opportunity))

        values = {
            "runnum": arr("runnum", -1, np.int64),
            "eventnum": arr("eventnum", -1, np.int64),
            "e_p": e_p, "e_theta": e_theta, "e_phi": e_phi,
            "p_p": p_p, "p_theta": p_theta, "p_phi": p_phi,
            "g_E": g_E, "g_theta": g_theta, "g_phi": g_phi,
            "probe_E": probe_E, "probe_theta": probe_theta, "probe_phi": probe_phi,
            "opportunity": opportunity,
        }
        for name in store:
            store[name].append(np.asarray(values[name])[candidate])
        # endfor
    # endfor

    payload = {}
    for name, parts in store.items():
        if not parts:
            dtype = bool if name == "opportunity" else (np.int64 if name in {"runnum", "eventnum"} else np.float64)
            payload[name] = np.asarray([], dtype=dtype)
        else:
            payload[name] = np.concatenate(parts)
        # endif
    # endfor

    records = EpgRecords(**payload)

    # Remove exact duplicate MC rows only. Distinct photons from one event remain.
    if mode == "mc" and records.size() > 0:
        matrix = np.column_stack([
            np.round(records.e_p, args.mc_signature_decimals),
            np.round(records.e_theta, args.mc_signature_decimals),
            np.round(records.e_phi, args.mc_signature_decimals),
            np.round(records.p_p, args.mc_signature_decimals),
            np.round(records.p_theta, args.mc_signature_decimals),
            np.round(records.p_phi, args.mc_signature_decimals),
            np.round(records.g_E, args.mc_signature_decimals),
            np.round(records.g_theta, args.mc_signature_decimals),
            np.round(records.g_phi, args.mc_signature_decimals),
        ])
        _, first = np.unique(matrix, axis=0, return_index=True)
        keep = np.zeros(records.size(), dtype=bool)
        keep[np.sort(first)] = True
        records = EpgRecords(**{name: np.asarray(getattr(records, name))[keep] for name in records.__dataclass_fields__})
        counts["after_exact_duplicate_removal"] = records.size()
    # endif

    return records, counts


def keys_for_records(records: EpgRecords, mode: str, decimals: int) -> np.ndarray:
    if mode == "data":
        return np.column_stack((records.runnum.astype(np.int64), records.eventnum.astype(np.int64)))
    # endif
    return np.column_stack([
        np.round(records.e_p, decimals),
        np.round(records.e_theta, decimals),
        np.round(records.e_phi, decimals),
        np.round(records.p_p, decimals),
        np.round(records.p_theta, decimals),
        np.round(records.p_phi, decimals),
    ])


def structured(matrix: np.ndarray) -> np.ndarray:
    matrix = np.ascontiguousarray(matrix)
    dtype = np.dtype([(f"f{i}", matrix.dtype) for i in range(matrix.shape[1])])
    return matrix.view(dtype).reshape(-1)


def read_identity_catalog(path: str, mode: str, args: argparse.Namespace) -> IdentityCatalog:
    logical = ("runnum", "eventnum", "e_p", "e_theta", "e_phi", "p_p", "p_theta", "p_phi")
    resolved = resolve(path, logical, optional=("runnum", "eventnum"))
    expressions = sorted({branch for branch in resolved.values() if branch is not None})
    parts: List[np.ndarray] = []

    for arrays in uproot.iterate(f"{path}:{TREE_NAME}", expressions=expressions, step_size=args.step_size, library="np"):
        def arr(name: str, default=math.nan, dtype=np.float64) -> np.ndarray:
            return finite_array(arrays, resolved.get(name), default, dtype)
        if mode == "data":
            matrix = np.column_stack((arr("runnum", -1, np.int64), arr("eventnum", -1, np.int64)))
        else:
            matrix = np.column_stack([
                np.round(arr("e_p"), args.mc_signature_decimals),
                np.round(arr("e_theta"), args.mc_signature_decimals),
                np.round(arr("e_phi"), args.mc_signature_decimals),
                np.round(arr("p_p"), args.mc_signature_decimals),
                np.round(arr("p_theta"), args.mc_signature_decimals),
                np.round(arr("p_phi"), args.mc_signature_decimals),
            ])
        # endif
        parts.append(matrix)
    # endfor

    matrix = np.concatenate(parts, axis=0) if parts else np.empty((0, 2 if mode == "data" else 6))
    keys = np.unique(structured(matrix))
    return IdentityCatalog(keys=keys, entries=int(matrix.shape[0]), unique_keys=int(keys.size))


def angular_distance_deg(theta1: float, phi1: float, theta2: float, phi2: float) -> float:
    cosine = (
        math.sin(theta1) * math.sin(theta2) * math.cos(phi1 - phi2)
        + math.cos(theta1) * math.cos(theta2)
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def pair_mass(E1: float, theta1: float, phi1: float, E2: float, theta2: float, phi2: float) -> float:
    opening = math.radians(angular_distance_deg(theta1, phi1, theta2, phi2))
    return math.sqrt(max(2.0 * E1 * E2 * (1.0 - math.cos(opening)), 0.0))


def residual_histogram_specs() -> Dict[str, np.ndarray]:
    return {
        "nearest_angle_deg": np.linspace(0.0, 60.0, 121),
        "nearest_relative_E": np.linspace(0.0, 3.0, 121),
        "nearest_pair_mass_GeV": np.linspace(0.0, 1.0, 121),
        "nearest_score": np.linspace(0.0, 100.0, 121),
    }


def update_reservoir(
    reservoir: List[Dict[str, object]],
    row: Dict[str, object],
    seen_rows: int,
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

    replacement_index = int(rng.integers(0, seen_rows))
    if replacement_index < capacity:
        reservoir[replacement_index] = row
    # endif


def audit_matches(
    records: EpgRecords,
    catalog: IdentityCatalog,
    mode: str,
    args: argparse.Namespace,
) -> Tuple[List[Dict[str, object]], Dict[str, object], Dict[str, object]]:
    """
    Audit every directed opportunity without retaining one Python dictionary per
    opportunity.

    Exact counters and histogram bin contents include the full sample.  A small
    deterministic reservoir sample is retained only for manual CSV inspection.
    """
    keys = keys_for_records(records, mode, args.mc_signature_decimals)
    key_structured = structured(keys)
    unique, inverse, counts = np.unique(
        key_structured,
        return_inverse=True,
        return_counts=True,
    )
    order = np.argsort(inverse, kind="stable")
    offsets = np.empty(counts.size + 1, dtype=np.int64)
    offsets[0] = 0
    np.cumsum(counts, out=offsets[1:])

    catalog_set = set(item.tobytes() for item in catalog.keys)
    opportunity_indices = np.flatnonzero(records.opportunity)

    sample_capacity = max(0, int(args.diagnostic_sample_size))
    seed_offset = 0 if mode == "data" else 1_000_003
    rng = np.random.default_rng(int(args.diagnostic_seed) + seed_offset)
    sampled_rows: List[Dict[str, object]] = []

    histogram_edges = residual_histogram_specs()
    histogram_counts = {
        key: np.zeros(len(edges) - 1, dtype=np.int64)
        for key, edges in histogram_edges.items()
    }

    opportunities_in_catalog = 0
    opportunities_with_non_tag = 0
    total_group_size = 0
    total_non_tag_candidates = 0
    nearest_rows_seen = 0

    # Reservoirs used only to estimate medians without retaining all residuals.
    median_reservoir_capacity = max(20_000, sample_capacity)
    angle_reservoir: List[float] = []
    relE_reservoir: List[float] = []
    residual_seen = 0

    angle_cuts = (1.0, 3.0, 5.0, 10.0, 20.0)
    relE_cuts = (0.15, 0.35, 0.50, 1.00)
    threshold_counts = {
        (angle_cut, relE_cut): 0
        for angle_cut in angle_cuts
        for relE_cut in relE_cuts
    }

    for index in opportunity_indices:
        group = int(inverse[index])
        members = order[offsets[group]:offsets[group + 1]]
        in_catalog = key_structured[index].tobytes() in catalog_set
        opportunities_in_catalog += int(in_catalog)
        total_group_size += int(members.size)

        best = None
        non_tag_candidates = 0
        for candidate in members:
            same_tag = (
                abs(records.g_E[candidate] - records.g_E[index]) < 1.0e-10
                and abs(records.g_theta[candidate] - records.g_theta[index]) < 1.0e-10
                and abs(records.g_phi[candidate] - records.g_phi[index]) < 1.0e-10
            )
            if same_tag:
                continue
            # endif

            non_tag_candidates += 1
            angle = angular_distance_deg(
                float(records.probe_theta[index]),
                float(records.probe_phi[index]),
                float(records.g_theta[candidate]),
                float(records.g_phi[candidate]),
            )
            relative_energy = abs(
                float(records.probe_E[index] - records.g_E[candidate])
            ) / max(float(records.g_E[candidate]), 1.0e-12)
            score = (
                (
                    angle
                    / max(args.score_angle_scale_deg, 1.0e-12)
                ) ** 2
                + (
                    relative_energy
                    / max(args.score_relative_E_scale, 1.0e-12)
                ) ** 2
            )
            invariant_mass = pair_mass(
                float(records.g_E[index]),
                float(records.g_theta[index]),
                float(records.g_phi[index]),
                float(records.g_E[candidate]),
                float(records.g_theta[candidate]),
                float(records.g_phi[candidate]),
            )
            if best is None or score < best[0]:
                best = (
                    score,
                    int(candidate),
                    angle,
                    relative_energy,
                    invariant_mass,
                )
            # endif
        # endfor

        total_non_tag_candidates += non_tag_candidates
        if best is not None:
            opportunities_with_non_tag += 1
            nearest_rows_seen += 1
            residual_seen += 1

            nearest_values = {
                "nearest_angle_deg": float(best[2]),
                "nearest_relative_E": float(best[3]),
                "nearest_pair_mass_GeV": float(best[4]),
                "nearest_score": float(best[0]),
            }

            for key, value in nearest_values.items():
                edges = histogram_edges[key]
                bin_index = int(np.searchsorted(edges, value, side="right") - 1)
                if 0 <= bin_index < histogram_counts[key].size:
                    histogram_counts[key][bin_index] += 1
                # endif
            # endfor

            for angle_cut in angle_cuts:
                for relE_cut in relE_cuts:
                    if best[2] < angle_cut and best[3] < relE_cut:
                        threshold_counts[(angle_cut, relE_cut)] += 1
                    # endif
                # endfor
            # endfor

            if len(angle_reservoir) < median_reservoir_capacity:
                angle_reservoir.append(float(best[2]))
                relE_reservoir.append(float(best[3]))
            else:
                replacement_index = int(rng.integers(0, residual_seen))
                if replacement_index < median_reservoir_capacity:
                    angle_reservoir[replacement_index] = float(best[2])
                    relE_reservoir[replacement_index] = float(best[3])
                # endif
            # endif
        # endif

        row = {
            "opportunity_index": int(index),
            "group_size": int(members.size),
            "non_tag_candidates": int(non_tag_candidates),
            "in_epgammagamma_catalog": bool(in_catalog),
            "tag_E": float(records.g_E[index]),
            "predicted_probe_E": float(records.probe_E[index]),
            "has_non_tag_candidate": best is not None,
            "nearest_candidate_index": int(best[1]) if best is not None else -1,
            "nearest_angle_deg": float(best[2]) if best is not None else math.nan,
            "nearest_relative_E": float(best[3]) if best is not None else math.nan,
            "nearest_pair_mass_GeV": float(best[4]) if best is not None else math.nan,
            "nearest_score": float(best[0]) if best is not None else math.nan,
        }
        update_reservoir(
            sampled_rows,
            row,
            seen_rows=int(index) + 1,
            capacity=sample_capacity,
            rng=rng,
        )
    # endfor

    directed_opportunities = int(opportunity_indices.size)
    summary: Dict[str, object] = {
        "mode": mode,
        "photon_candidate_records": records.size(),
        "directed_opportunities": directed_opportunities,
        "opportunities_in_epgammagamma_catalog": opportunities_in_catalog,
        "opportunities_with_non_tag_candidate": opportunities_with_non_tag,
        "fraction_with_non_tag_candidate": (
            opportunities_with_non_tag / directed_opportunities
            if directed_opportunities > 0
            else math.nan
        ),
        "group_size_mean": (
            total_group_size / directed_opportunities
            if directed_opportunities > 0
            else math.nan
        ),
        "non_tag_candidates_per_opportunity_mean": (
            total_non_tag_candidates / directed_opportunities
            if directed_opportunities > 0
            else math.nan
        ),
        "nearest_angle_deg_median": (
            float(np.median(angle_reservoir))
            if angle_reservoir
            else math.nan
        ),
        "nearest_relative_E_median": (
            float(np.median(relE_reservoir))
            if relE_reservoir
            else math.nan
        ),
        "median_estimate_reservoir_size": len(angle_reservoir),
        "diagnostic_rows_written": len(sampled_rows),
        "diagnostic_rows_total_population": directed_opportunities,
    }

    for angle_cut in angle_cuts:
        for relE_cut in relE_cuts:
            key = (
                f"fraction_angle_lt_{angle_cut:g}"
                f"_and_relE_lt_{relE_cut:g}"
            )
            summary[key] = (
                threshold_counts[(angle_cut, relE_cut)]
                / directed_opportunities
                if directed_opportunities > 0
                else math.nan
            )
        # endfor
    # endfor

    histogram_payload: Dict[str, object] = {
        key: {
            "edges": edges.tolist(),
            "counts": histogram_counts[key].tolist(),
        }
        for key, edges in histogram_edges.items()
    }
    histogram_payload["nearest_match_entries"] = nearest_rows_seen

    return sampled_rows, summary, histogram_payload



def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        return
    # endif
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    # endwith


def plot_residuals(
    path: Path,
    title: str,
    histogram_payload: Mapping[str, object],
) -> None:
    panels = (
        ("nearest_angle_deg", "Nearest angular residual (deg)"),
        ("nearest_relative_E", "Nearest relative energy residual"),
        ("nearest_pair_mass_GeV", r"Tag-partner mass (GeV)"),
        ("nearest_score", "Nearest matching score"),
    )

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    for axis, (key, label) in zip(axes.flat, panels):
        payload = histogram_payload[key]
        edges = np.asarray(payload["edges"], dtype=float)
        counts = np.asarray(payload["counts"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        axis.step(centers, counts, where="mid", linewidth=1.3)
        axis.set_xlabel(label)
        axis.set_ylabel("Entries")
        axis.grid(alpha=0.25)
    # endfor

    fig.suptitle(
        f"{title}\n"
        f"Nearest-match entries: "
        f"{int(histogram_payload.get('nearest_match_entries', 0)):,}"
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=180)
    plt.close(fig)



def process_period(period: PeriodConfig, args_dict: Mapping[str, object]) -> Tuple[str, Dict[str, object]]:
    args = argparse.Namespace(**args_dict)
    period_dir = Path(args.output_dir) / period.key
    period_dir.mkdir(parents=True, exist_ok=True)

    data, data_counts = read_epg(period.epg_data, period.beam_energy_GeV, "data", args)
    mc, mc_counts = read_epg(period.pi0_epg_mc, period.beam_energy_GeV, "mc", args)
    data_catalog = read_identity_catalog(period.epgg_data, "data", args)
    mc_catalog = read_identity_catalog(period.pi0_epgg_mc, "mc", args)

    data_rows, data_summary, data_histograms = audit_matches(
        data, data_catalog, "data", args
    )
    mc_rows, mc_summary, mc_histograms = audit_matches(
        mc, mc_catalog, "mc", args
    )

    write_csv(period_dir / "data_nearest_match_sample.csv", data_rows)
    write_csv(period_dir / "aaogen_nearest_match_sample.csv", mc_rows)
    plot_residuals(
        period_dir / "data_matching_residuals.png",
        f"{period.label}: data no-cuts matching",
        data_histograms,
    )
    plot_residuals(
        period_dir / "aaogen_matching_residuals.png",
        f"{period.label}: AAOGEN no-cuts matching",
        mc_histograms,
    )

    payload = {
        "period": asdict(period),
        "data_read_counts": data_counts,
        "mc_read_counts": mc_counts,
        "data_catalog": asdict(data_catalog) | {"keys": None},
        "mc_catalog": asdict(mc_catalog) | {"keys": None},
        "data_matching": data_summary,
        "mc_matching": mc_summary,
        "data_residual_histograms": data_histograms,
        "mc_residual_histograms": mc_histograms,
    }
    with open(period_dir / "metadata.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
    # endwith

    log(
        f"{period.label} data: opportunities={data_summary['directed_opportunities']:,}, "
        f"with another photon={data_summary['fraction_with_non_tag_candidate']:.3%}, "
        f"median dtheta={data_summary['nearest_angle_deg_median']:.3f} deg, "
        f"median dE/E={data_summary['nearest_relative_E_median']:.3f}"
    )
    log(
        f"{period.label} AAOGEN: opportunities={mc_summary['directed_opportunities']:,}, "
        f"with another photon={mc_summary['fraction_with_non_tag_candidate']:.3%}, "
        f"median dtheta={mc_summary['nearest_angle_deg_median']:.3f} deg, "
        f"median dE/E={mc_summary['nearest_relative_E_median']:.3f}"
    )
    return period.key, payload


def main() -> int:
    args = parse_args()
    periods = selected_periods(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "script": Path(__file__).name,
        "arguments": vars(args),
        "created_unix_time": time.time(),
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
        f"NO-CUTS MATCHING AUDIT: {len(periods)} period(s), {workers} worker(s). "
        f"Only E_tag >= {args.tag_E_min:g} GeV and predicted E_X >= "
        f"{args.probe_E_min:g} GeV are required. "
        f"At most {args.diagnostic_sample_size:,} diagnostic rows per sample "
        f"and period will be written."
    )

    metadata: Dict[str, object] = {}
    if workers == 1:
        for period in periods:
            key, payload = process_period(period, vars(args))
            metadata[key] = payload
        # endfor
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [executor.submit(process_period, period, vars(args)) for period in periods]
            for future in concurrent.futures.as_completed(futures):
                key, payload = future.result()
                metadata[key] = payload
                log(f"Completed {key}.")
            # endfor
        # endwith
    # endif

    summary_rows = []
    for period in periods:
        payload = metadata[period.key]
        for sample in ("data", "mc"):
            summary = payload[f"{sample}_matching"]
            row = {
                "period": period.key,
                "period_label": period.label,
                "sample": "Data" if sample == "data" else "AAOGEN MC",
                **summary,
            }
            summary_rows.append(row)
        # endfor
    # endfor

    write_csv(output_dir / "matching_audit_summary.csv", summary_rows)
    with open(output_dir / "matching_audit_summary.json", "w", encoding="utf-8") as handle:
        json.dump({"rows": summary_rows, "period_metadata": metadata}, handle, indent=2)
    # endwith

    log(f"Wrote no-cuts matching audit to {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr, flush=True)
        raise
    # endtry
