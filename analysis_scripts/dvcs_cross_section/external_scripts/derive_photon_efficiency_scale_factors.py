#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v039.py

Development script for the relative data/MC photon-reconstruction efficiency
measurement in CLAS12 RGA.

CURRENT USER-FACING WORKFLOW
----------------------------
1. Denominator composition:
   real epgamma data are decomposed into exclusive-pi0 missing-photon tags and
   genuine BH/DVCS epgamma using aaogen and dvcsgen templates.

2. Photon-efficiency extraction:
   reconstructed companion photons are associated with those tags in data and
   MC and the relative efficiency scale factor is formed,

       SF_gamma = epsilon_data / epsilon_MC.

The old "Stage I" resolution study is no longer a user-visible analysis stage.
Its only surviving role is an INTERNAL aaogen MC association kernel. Because
runnum/evnum are not useful identifiers in MC, aaogen epgamma and reconstructed
eppi0 skims are matched through reconstructed electron/proton Cartesian
momenta. That association supplies the MC reconstructed-probe numerator. No
routine Stage-I plots, resolution tables, summaries, or NPZ files are written.

Data companion-photon association uses exact (runnum, evnum) identity.

The current wagon-derived data support the reference extraction only for
0.4 <= E_probe^pred < 2 GeV. The result therefore remains a PRELIMINARY
overlap/reference extraction until the loose nSidis-derived data reproduce it
below 2 GeV and permit extension above 2 GeV.

The nominal denominator model remains aaogen pi0 + dvcsgen BH/DVCS. Mixed data
remain diagnostic-only. In FT, a coarse 0.40--1.20 / 1.20--2.00 GeV diagnostic
also floats the deterministic mixed-event wrong-photon template as a third
component with the two-component morph nuisance values held fixed. This does
not feed the nominal efficiency extraction.

Routine outputs:
    output/photon_efficiency/summary/
    output/photon_efficiency/stage2/
    output/photon_efficiency/stage3/

Default --plot-mode compact keeps only high-value plots. --plot-mode full
restores development/debug canvases.

Typical full-statistics run:
    python derive_photon_efficiency_scale_factors_v039.py --max-entries 0

One period:
    python derive_photon_efficiency_scale_factors_v039.py \
        --max-entries 0 --period fa18_out
"""

from __future__ import annotations

import argparse
import csv
import gc
import json
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy.spatial import cKDTree
from scipy.optimize import minimize, minimize_scalar
from scipy.ndimage import gaussian_filter1d, shift as ndimage_shift


# -----------------------------------------------------------------------------
# Physics constants
# -----------------------------------------------------------------------------

M_E = 0.00051099895000
M_P = 0.93827208816
M_PI0 = 0.1349768
TWO_PI = 2.0 * math.pi


# Probe-energy binning used for Stage-Ib resolution studies.
# Confirmed photon detector-code convention for reconstructed eppi0 trees.
# Plot aesthetics aligned with the main exclusivity-determination script.
SAMPLE_COLORS = {
    "data": "black",
    "dvcs_mc": "tab:blue",
    "pi0_mc": "tab:red",
    "fit": "tab:green",
}
SAMPLE_LABELS = {
    "data": r"$e'p'\gamma$ data",
    "dvcs_mc": "BH/DVCS MC",
    "pi0_mc": r"$e'p'\pi^0$ MC contribution",
}

PHOTON_DETECTOR_FT = 0
PHOTON_DETECTOR_FD = 1

PROBE_ENERGY_EDGES = np.asarray(
    [0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 2.00, 3.00, 4.50, 6.00, 10.00],
    dtype=float,
)


# -----------------------------------------------------------------------------
# Period configuration
# -----------------------------------------------------------------------------

@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy: float
    epgamma_data: str
    eppi0_data: str
    epgamma_dvcs_mc: str
    epgamma_pi0_mc: str
    eppi0_pi0_mc: str


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
        "sp19_inb", "Sp19 Inb", 10.1998,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_50nA_10200MeV.root",
    ),
)


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def wrap_phi(phi: np.ndarray | float) -> np.ndarray | float:
    return (phi + math.pi) % TWO_PI - math.pi


def angle_between(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom <= 0.0:
        return float("nan")
    #endif
    c = float(np.dot(a, b) / denom)
    return math.acos(max(-1.0, min(1.0, c)))


def direction(theta: float, phi: float) -> np.ndarray:
    st = math.sin(theta)
    return np.asarray([st * math.cos(phi), st * math.sin(phi), math.cos(theta)], dtype=float)


def cartesian_from_spherical(p: np.ndarray, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    st = np.sin(theta)
    return np.column_stack((
        p * st * np.cos(phi),
        p * st * np.sin(phi),
        p * np.cos(theta),
    ))


def infer_angle_unit(arrays: Dict[str, np.ndarray], requested: str) -> str:
    if requested in ("rad", "deg"):
        return requested
    #endif

    probes: List[np.ndarray] = []
    for name in ("e_theta", "p1_theta", "p2_theta", "e_phi", "p1_phi", "p2_phi"):
        if name in arrays:
            v = np.asarray(arrays[name], dtype=float)
            v = v[np.isfinite(v)]
            if v.size:
                probes.append(np.abs(v[: min(v.size, 10000)]))
            #endif
        #endif
    #endfor

    if not probes:
        raise RuntimeError("Could not infer angle units: no angular branches were available.")
    #endif

    vmax = float(np.nanpercentile(np.concatenate(probes), 99.0))
    return "deg" if vmax > 7.0 else "rad"


def to_radians(values: np.ndarray, unit: str) -> np.ndarray:
    return np.deg2rad(values) if unit == "deg" else np.asarray(values, dtype=float)


def finite_mask(*arrays: np.ndarray) -> np.ndarray:
    out = np.ones(len(arrays[0]), dtype=bool)
    for arr in arrays:
        out &= np.isfinite(arr)
    #endfor
    return out


def find_tree(root_file: uproot.ReadOnlyDirectory, requested: Optional[str]) -> Tuple[str, object]:
    if requested:
        if requested not in root_file:
            raise KeyError(f"Requested tree '{requested}' not found. Keys: {root_file.keys()[:20]}")
        #endif
        return requested, root_file[requested]
    #endif

    for key, classname in root_file.classnames().items():
        if "TTree" in classname or "RNTuple" in classname:
            clean = key.split(";")[0]
            return clean, root_file[clean]
        #endif
    #endfor

    raise RuntimeError(f"No TTree/RNTuple found. Keys: {root_file.keys()[:20]}")


def read_branches(
    path: str,
    required: Sequence[str],
    optional: Sequence[str],
    tree_name: Optional[str],
    max_entries: int,
) -> Tuple[Dict[str, np.ndarray], str, int]:
    with uproot.open(path) as root_file:
        found_name, tree = find_tree(root_file, tree_name)
        available = set(tree.keys())

        missing = sorted(set(required) - available)
        if missing:
            raise KeyError(
                f"{path}\nTree '{found_name}' is missing required branches:\n  "
                + "\n  ".join(missing)
            )
        #endif

        chosen = list(required) + [b for b in optional if b in available and b not in required]
        entry_stop = None if max_entries == 0 else min(max_entries, int(tree.num_entries))
        arrays = tree.arrays(chosen, entry_stop=entry_stop, library="np")
        arrays = {k: np.asarray(v) for k, v in arrays.items()}
        return arrays, found_name, int(tree.num_entries)


# -----------------------------------------------------------------------------
# Schema adapters
# -----------------------------------------------------------------------------

EPG_REQUIRED = (
    "e_p", "e_theta", "e_phi",
    "p1_p", "p1_theta", "p1_phi",
    "p2_p", "p2_theta", "p2_phi",
)

EPG_OPTIONAL_PI0_MC = (
    "detector2",
    "pTmiss",
    "theta_gamma_gamma",
)

EPG_OPTIONAL_DVCS_MC = (
    "pTmiss",
    "theta_gamma_gamma",
)

EPG_OPTIONAL_DATA = (
    "runnum",
    "evnum",
    "detector2",
    "pTmiss",
    "theta_gamma_gamma",
)

EPPIO_REQUIRED = (
    "e_p", "e_theta", "e_phi",
    "p1_p", "p1_theta", "p1_phi",
    "p2_p", "p2_theta", "p2_phi",
    "Mh_gammagamma",
)

EPPIO_OPTIONAL_PI0_MC = (
    "detector_gamma1",
    "detector_gamma2",
    "Emiss2",
    "pTmiss",
    "theta_pi0_pi0",
)

EPPIO_OPTIONAL_DATA = (
    "runnum",
    "evnum",
    "detector_gamma1",
    "detector_gamma2",
    "Emiss2",
    "pTmiss",
    "theta_pi0_pi0",
)


@dataclass
class EPGammaSample:
    electron_p3: np.ndarray
    proton_p3: np.ndarray
    tag_p3: np.ndarray
    tag_energy: np.ndarray
    raw: Dict[str, np.ndarray]
    angle_unit: str


@dataclass
class EPPi0Sample:
    electron_p3: np.ndarray
    proton_p3: np.ndarray
    pi0_p3: np.ndarray
    pi0_p: np.ndarray
    pi0_mass: np.ndarray
    pi0_theta: np.ndarray
    raw: Dict[str, np.ndarray]
    angle_unit: str


def _optional_raw_only(
    arrays: Dict[str, np.ndarray],
    required: Sequence[str],
) -> Dict[str, np.ndarray]:
    required_set = set(required)
    return {
        key: np.asarray(value)
        for key, value in arrays.items()
        if key not in required_set
    }


def extract_epgamma(arrays: Dict[str, np.ndarray], angle_mode: str) -> EPGammaSample:
    unit = infer_angle_unit(arrays, angle_mode)

    e_theta = to_radians(arrays["e_theta"], unit)
    e_phi = to_radians(arrays["e_phi"], unit)
    p_theta = to_radians(arrays["p1_theta"], unit)
    p_phi = to_radians(arrays["p1_phi"], unit)
    g_theta = to_radians(arrays["p2_theta"], unit)
    g_phi = to_radians(arrays["p2_phi"], unit)

    e_p = np.asarray(arrays["e_p"], dtype=float)
    p_p = np.asarray(arrays["p1_p"], dtype=float)
    g_e = np.asarray(arrays["p2_p"], dtype=float)

    return EPGammaSample(
        electron_p3=cartesian_from_spherical(e_p, e_theta, e_phi),
        proton_p3=cartesian_from_spherical(p_p, p_theta, p_phi),
        tag_p3=cartesian_from_spherical(g_e, g_theta, g_phi),
        tag_energy=g_e,
        raw=_optional_raw_only(arrays, EPG_REQUIRED),
        angle_unit=unit,
    )


def extract_eppi0(arrays: Dict[str, np.ndarray], angle_mode: str) -> EPPi0Sample:
    unit = infer_angle_unit(arrays, angle_mode)

    e_theta = to_radians(arrays["e_theta"], unit)
    e_phi = to_radians(arrays["e_phi"], unit)
    p_theta = to_radians(arrays["p1_theta"], unit)
    p_phi = to_radians(arrays["p1_phi"], unit)
    pi_theta = to_radians(arrays["p2_theta"], unit)
    pi_phi = to_radians(arrays["p2_phi"], unit)

    e_p = np.asarray(arrays["e_p"], dtype=float)
    p_p = np.asarray(arrays["p1_p"], dtype=float)
    pi_p = np.asarray(arrays["p2_p"], dtype=float)
    pi_mass = np.asarray(arrays["Mh_gammagamma"], dtype=float)

    return EPPi0Sample(
        electron_p3=cartesian_from_spherical(e_p, e_theta, e_phi),
        proton_p3=cartesian_from_spherical(p_p, p_theta, p_phi),
        pi0_p3=cartesian_from_spherical(pi_p, pi_theta, pi_phi),
        pi0_p=pi_p,
        pi0_mass=pi_mass,
        pi0_theta=pi_theta,
        raw=_optional_raw_only(arrays, EPPIO_REQUIRED),
        angle_unit=unit,
    )


# -----------------------------------------------------------------------------
# MC event matching through reconstructed e/p kinematics
# -----------------------------------------------------------------------------

@dataclass
class MatchResult:
    epg_index: np.ndarray
    eppi0_index: np.ndarray
    nearest_distance: np.ndarray
    max_component_delta: np.ndarray


def match_parent_kinematics(
    epg: EPGammaSample,
    eppi0: EPPi0Sample,
    tag_min: float,
    tag_max: float,
    component_tolerance: float,
    nearest_distance_max: float,
    kdtree_workers: int = 1,
    query_chunk_size: int = 500_000,
) -> MatchResult:
    """
    Same MC parent match as before, executed in bounded epgamma query chunks.
    """
    good_pi0 = (
        np.all(np.isfinite(eppi0.electron_p3), axis=1)
        & np.all(np.isfinite(eppi0.proton_p3), axis=1)
    )
    good_epg = (
        np.all(np.isfinite(epg.electron_p3), axis=1)
        & np.all(np.isfinite(epg.proton_p3), axis=1)
        & np.isfinite(epg.tag_energy)
        & (epg.tag_energy >= tag_min)
        & (epg.tag_energy < tag_max)
    )

    pi0_indices = np.flatnonzero(good_pi0)
    epg_indices = np.flatnonzero(good_epg)

    if pi0_indices.size == 0 or epg_indices.size == 0:
        return MatchResult(
            epg_index=np.asarray([], dtype=np.int64),
            eppi0_index=np.asarray([], dtype=np.int64),
            nearest_distance=np.asarray([], dtype=float),
            max_component_delta=np.asarray([], dtype=float),
        )
    #endif

    pi0_search = np.empty((pi0_indices.size, 6), dtype=np.float32)
    pi0_search[:, :3] = (
        eppi0.electron_p3[pi0_indices] / component_tolerance
    ).astype(np.float32, copy=False)
    pi0_search[:, 3:] = (
        eppi0.proton_p3[pi0_indices] / component_tolerance
    ).astype(np.float32, copy=False)

    tree = cKDTree(pi0_search, compact_nodes=True, balanced_tree=True)
    del pi0_search

    chunk_size = max(10_000, int(query_chunk_size))
    epg_parts: List[np.ndarray] = []
    pi0_parts: List[np.ndarray] = []
    distance_parts: List[np.ndarray] = []
    delta_parts: List[np.ndarray] = []

    for start in range(0, epg_indices.size, chunk_size):
        stop = min(start + chunk_size, epg_indices.size)
        idx = epg_indices[start:stop]

        query = np.empty((len(idx), 6), dtype=np.float32)
        query[:, :3] = (
            epg.electron_p3[idx] / component_tolerance
        ).astype(np.float32, copy=False)
        query[:, 3:] = (
            epg.proton_p3[idx] / component_tolerance
        ).astype(np.float32, copy=False)

        distance, local_index = tree.query(
            query,
            k=1,
            workers=max(1, int(kdtree_workers)),
        )
        del query

        candidate_pi0 = pi0_indices[local_index]
        de = np.abs(
            epg.electron_p3[idx] - eppi0.electron_p3[candidate_pi0]
        )
        dp = np.abs(
            epg.proton_p3[idx] - eppi0.proton_p3[candidate_pi0]
        )
        max_delta = np.maximum(np.max(de, axis=1), np.max(dp, axis=1))

        accepted = (
            np.isfinite(distance)
            & (distance <= nearest_distance_max)
            & (max_delta <= component_tolerance)
        )
        if np.any(accepted):
            epg_parts.append(idx[accepted])
            pi0_parts.append(candidate_pi0[accepted])
            distance_parts.append(distance[accepted])
            delta_parts.append(max_delta[accepted])
        #endif
    #endfor

    if not epg_parts:
        return MatchResult(
            epg_index=np.asarray([], dtype=np.int64),
            eppi0_index=np.asarray([], dtype=np.int64),
            nearest_distance=np.asarray([], dtype=float),
            max_component_delta=np.asarray([], dtype=float),
        )
    #endif

    return MatchResult(
        epg_index=np.concatenate(epg_parts),
        eppi0_index=np.concatenate(pi0_parts),
        nearest_distance=np.concatenate(distance_parts),
        max_component_delta=np.concatenate(delta_parts),
    )



def packed_event_keys(runnum: np.ndarray, evnum: np.ndarray) -> np.ndarray:
    """Pack nonnegative (runnum, evnum) data identifiers into uint64 keys."""
    run = np.asarray(runnum, dtype=np.int64)
    ev = np.asarray(evnum, dtype=np.int64)
    if np.any(run < 0) or np.any(ev < 0):
        raise ValueError("Data runnum/evnum must be nonnegative for exact matching.")
    #endif
    if np.any(run > np.iinfo(np.uint32).max) or np.any(ev > np.iinfo(np.uint32).max):
        raise ValueError("Data runnum/evnum exceeds the uint32 packing range.")
    #endif
    return (run.astype(np.uint64) << np.uint64(32)) | ev.astype(np.uint64)



@dataclass
class DataAssociationResult:
    """
    Detailed exact-event Stage-III data association.

    The current eppi0 tree stores the reconstructed pi0 four-vector but not the
    individual daughter-photon four-vectors. Therefore photon-level identity is
    tested kinematically: removing the epgamma tag from P_pi0 must leave a
    positive-energy, approximately massless companion, and that companion is
    optionally required to agree loosely with the probe predicted from epgamma
    exclusivity.
    """
    best_match: MatchResult
    stage_lookup: Dict[str, np.ndarray]
    counters: Dict[str, int]
    best_diagnostics: Dict[str, np.ndarray]


def _lookup_from_indices(n: int, indices: np.ndarray) -> np.ndarray:
    out = np.zeros(int(n), dtype=bool)
    if len(indices):
        out[np.asarray(indices, dtype=np.int64)] = True
    #endif
    return out


def match_data_event_candidates(
    period: PeriodConfig,
    epg: EPGammaSample,
    eppi0: EPPi0Sample,
    tag_remainder_m2_max: float,
    reco_probe_energy_min: float,
    probe_angle_max_deg: float,
    probe_frac_energy_max: float,
) -> DataAssociationResult:
    """
    Exact-event Stage-III DATA association.

    For data, (runnum, evnum) is the complete parent-event identity. No
    electron/proton kinematic matching is performed. MC continues to use the
    independent e/p kinematic matcher developed in Stage I.

    For every epgamma tag candidate, all same-event reconstructed eppi0
    candidates are considered. Candidate association is determined entirely
    from the photon side:

      same_event
          At least one eppi0 candidate has identical runnum and evnum.

      positive_remainder
          P_pi0 - k_tag has positive energy and non-zero momentum.

      tag_mass_shell
          The remainder is photon-like:
              |(P_pi0-k_tag)^2| <= tag_remainder_m2_max.

      probe_energy
          The reconstructed companion energy is above the photon threshold.

      probe_pred_consistent
          The reconstructed companion agrees loosely with the independent
          epgamma-predicted probe in opening angle and fractional energy.

    The best candidate is selected by the photon remainder mass shell, then by
    predicted/reconstructed companion agreement. No data e/p momentum
    difference enters the score or any hard selection.
    """
    for name in ("runnum", "evnum"):
        if name not in epg.raw or name not in eppi0.raw:
            raise KeyError(
                f"Stage III exact data matching requires '{name}' in both "
                "epgamma and eppi0 data trees."
            )
        #endif
    #endfor

    n_epg = len(epg.tag_energy)
    n_pi0 = len(eppi0.pi0_mass)
    empty_match = MatchResult(
        epg_index=np.asarray([], dtype=np.int64),
        eppi0_index=np.asarray([], dtype=np.int64),
        nearest_distance=np.asarray([], dtype=float),
        max_component_delta=np.asarray([], dtype=float),
    )
    stage_names = (
        "same_event",
        "positive_remainder",
        "tag_mass_shell",
        "probe_energy",
        "probe_pred_consistent",
    )
    empty_lookup = {
        name: np.zeros(n_epg, dtype=bool)
        for name in stage_names
    }

    def empty_result(extra_counters: Optional[Dict[str, int]] = None) -> DataAssociationResult:
        counters = {
            "epgamma_entries": int(n_epg),
            "eppi0_entries": int(n_pi0),
            **{name: 0 for name in stage_names},
            "same_event_candidate_pairs": 0,
            "positive_remainder_pairs": 0,
            "tag_mass_shell_pairs": 0,
            "probe_energy_pairs": 0,
            "probe_pred_consistent_pairs": 0,
            "accepted_best_matches": 0,
        }
        if extra_counters:
            counters.update(extra_counters)
        #endif
        return DataAssociationResult(
            best_match=empty_match,
            stage_lookup={k: v.copy() for k, v in empty_lookup.items()},
            counters=counters,
            best_diagnostics={},
        )
    #enddef

    if n_epg == 0 or n_pi0 == 0:
        return empty_result()
    #endif

    key_g = packed_event_keys(epg.raw["runnum"], epg.raw["evnum"])
    key_p = packed_event_keys(eppi0.raw["runnum"], eppi0.raw["evnum"])

    # Only the eppi0 keys need sorting. Every epgamma candidate then gets the
    # contiguous same-event eppi0 slice by two vectorized searchsorted calls.
    order_p = np.argsort(key_p, kind="mergesort")
    sorted_p = key_p[order_p]
    left = np.searchsorted(sorted_p, key_g, side="left")
    right = np.searchsorted(sorted_p, key_g, side="right")
    counts = (right - left).astype(np.int64)

    same_lookup = counts > 0
    total_pairs = int(np.sum(counts))
    if total_pairs == 0:
        result = empty_result({
            "same_event": int(np.count_nonzero(same_lookup)),
        })
        result.stage_lookup["same_event"] = same_lookup
        return result
    #endif

    # Expand only same-event candidate pairs.
    epg_rep = np.repeat(np.arange(n_epg, dtype=np.int64), counts)
    left_rep = np.repeat(left.astype(np.int64), counts)
    prefix_before = np.repeat((np.cumsum(counts) - counts), counts)
    offsets = np.arange(total_pairs, dtype=np.int64) - prefix_before
    pi_rep = order_p[left_rep + offsets]

    # Reconstructed companion: P_pi0 - k_tag.
    tag3 = epg.tag_p3[epg_rep]
    pi3 = eppi0.pi0_p3[pi_rep]
    tag_E = np.linalg.norm(tag3, axis=1)
    pi_p = np.linalg.norm(pi3, axis=1)
    pi_m = np.asarray(eppi0.pi0_mass[pi_rep], dtype=float)
    pi_E = np.sqrt(np.maximum(0.0, pi_p * pi_p + pi_m * pi_m))

    reco_E = pi_E - tag_E
    reco3 = pi3 - tag3
    reco_p = np.linalg.norm(reco3, axis=1)
    reco_m2 = reco_E * reco_E - reco_p * reco_p

    positive = (
        np.isfinite(reco_E)
        & np.isfinite(reco_p)
        & (reco_E > 0.0)
        & (reco_p > 0.0)
    )
    mass_shell = (
        positive
        & np.isfinite(reco_m2)
        & (np.abs(reco_m2) <= tag_remainder_m2_max)
    )
    above_threshold = mass_shell & (reco_E >= reco_probe_energy_min)

    # Independently predicted companion from the epgamma kinematics.
    e3 = epg.electron_p3[epg_rep]
    p3 = epg.proton_p3[epg_rep]
    e_pmag = np.linalg.norm(e3, axis=1)
    p_pmag = np.linalg.norm(p3, axis=1)
    e_E = np.sqrt(e_pmag * e_pmag + M_E * M_E)
    p_E = np.sqrt(p_pmag * p_pmag + M_P * M_P)

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    initial_E = period.beam_energy + M_P
    initial_p3 = np.asarray([0.0, 0.0, beam_p], dtype=float)

    pred_E = initial_E - e_E - p_E - tag_E
    pred3 = initial_p3[None, :] - e3 - p3 - tag3
    pred_p = np.linalg.norm(pred3, axis=1)

    dot = np.sum(pred3 * reco3, axis=1)
    denom = pred_p * reco_p
    with np.errstate(invalid="ignore", divide="ignore"):
        cos_open = np.clip(dot / denom, -1.0, 1.0)
        delta_theta_deg = np.degrees(np.arccos(cos_open))
        frac_delta_E = (reco_E - pred_E) / pred_E
    #endwith

    pred_consistent = (
        above_threshold
        & np.isfinite(pred_E)
        & (pred_E > 0.0)
        & np.isfinite(delta_theta_deg)
        & (delta_theta_deg <= probe_angle_max_deg)
        & np.isfinite(frac_delta_E)
        & (np.abs(frac_delta_E) <= probe_frac_energy_max)
    )

    pair_stage_masks = {
        "positive_remainder": positive,
        "tag_mass_shell": mass_shell,
        "probe_energy": above_threshold,
        "probe_pred_consistent": pred_consistent,
    }
    stage_lookup: Dict[str, np.ndarray] = {
        "same_event": same_lookup.copy(),
    }
    for name, pair_mask in pair_stage_masks.items():
        stage_lookup[name] = _lookup_from_indices(
            n_epg,
            np.unique(epg_rep[pair_mask]),
        )
    #endfor

    # Select one reconstructed-pi0 candidate per epgamma tag. Only photon-side
    # quantities enter the score.
    eligible_idx = np.flatnonzero(above_threshold)
    if eligible_idx.size == 0:
        best_match = empty_match
        best_diag: Dict[str, np.ndarray] = {}
    else:
        angle_term = np.nan_to_num(
            delta_theta_deg[eligible_idx] / max(probe_angle_max_deg, 1.0e-9),
            nan=1.0e6,
            posinf=1.0e6,
            neginf=1.0e6,
        )
        energy_term = np.nan_to_num(
            np.abs(frac_delta_E[eligible_idx])
            / max(probe_frac_energy_max, 1.0e-9),
            nan=1.0e6,
            posinf=1.0e6,
            neginf=1.0e6,
        )
        score = (
            np.abs(reco_m2[eligible_idx])
            / max(tag_remainder_m2_max, 1.0e-12)
            + 0.05 * angle_term
            + 0.05 * energy_term
        )

        gvalid = epg_rep[eligible_idx]
        order_best = np.lexsort((score, gvalid))
        sorted_eligible = eligible_idx[order_best]
        sorted_g = epg_rep[sorted_eligible]
        first = np.ones(len(sorted_eligible), dtype=bool)
        if len(sorted_eligible) > 1:
            first[1:] = sorted_g[1:] != sorted_g[:-1]
        #endif
        chosen = sorted_eligible[first]

        best_match = MatchResult(
            epg_index=np.asarray(epg_rep[chosen], dtype=np.int64),
            eppi0_index=np.asarray(pi_rep[chosen], dtype=np.int64),
            nearest_distance=np.zeros(len(chosen), dtype=float),
            max_component_delta=np.zeros(len(chosen), dtype=float),
        )
        best_diag = {
            "epg_index": np.asarray(epg_rep[chosen], dtype=np.int64),
            "eppi0_index": np.asarray(pi_rep[chosen], dtype=np.int64),
            "reco_probe_energy": np.asarray(reco_E[chosen], dtype=float),
            "pred_probe_energy": np.asarray(pred_E[chosen], dtype=float),
            "reco_probe_mass2": np.asarray(reco_m2[chosen], dtype=float),
            "probe_delta_theta_deg": np.asarray(
                delta_theta_deg[chosen], dtype=float
            ),
            "probe_delta_E_over_E": np.asarray(
                frac_delta_E[chosen], dtype=float
            ),
            "passes_pred_consistency": np.asarray(
                pred_consistent[chosen], dtype=bool
            ),
        }
    #endif

    counters = {
        "epgamma_entries": int(n_epg),
        "eppi0_entries": int(n_pi0),
        "same_event": int(np.count_nonzero(stage_lookup["same_event"])),
        "positive_remainder": int(
            np.count_nonzero(stage_lookup["positive_remainder"])
        ),
        "tag_mass_shell": int(
            np.count_nonzero(stage_lookup["tag_mass_shell"])
        ),
        "probe_energy": int(
            np.count_nonzero(stage_lookup["probe_energy"])
        ),
        "probe_pred_consistent": int(
            np.count_nonzero(stage_lookup["probe_pred_consistent"])
        ),
        "same_event_candidate_pairs": total_pairs,
        "positive_remainder_pairs": int(np.count_nonzero(positive)),
        "tag_mass_shell_pairs": int(np.count_nonzero(mass_shell)),
        "probe_energy_pairs": int(np.count_nonzero(above_threshold)),
        "probe_pred_consistent_pairs": int(
            np.count_nonzero(pred_consistent)
        ),
        "accepted_best_matches": int(len(best_match.epg_index)),
    }

    return DataAssociationResult(
        best_match=best_match,
        stage_lookup=stage_lookup,
        counters=counters,
        best_diagnostics=best_diag,
    )


def stage1_success_lookup(
    n_epgamma: int,
    arrays: Dict[str, np.ndarray],
    clean_mask: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (any_matched, clean_reconstructed_probe) booleans by epgamma row."""
    matched = np.zeros(int(n_epgamma), dtype=bool)
    clean = np.zeros(int(n_epgamma), dtype=bool)
    if arrays and "epg_index" in arrays:
        idx = np.asarray(arrays["epg_index"], dtype=np.int64)
        matched[idx] = True
        if len(clean_mask) != len(idx):
            raise ValueError("Stage-I clean mask length does not match epg_index.")
        #endif
        clean[idx[np.asarray(clean_mask, dtype=bool)]] = True
    #endif
    return matched, clean


# -----------------------------------------------------------------------------
# Stage-I / Stage-Ib pair construction
# -----------------------------------------------------------------------------

def build_stage1_arrays(
    period: PeriodConfig,
    epg: EPGammaSample,
    eppi0: EPPi0Sample,
    matches: MatchResult,
) -> Tuple[Dict[str, np.ndarray], Dict[str, int]]:
    """
    Vectorized construction of all Stage-I quantities.

    For each epgamma tag matched to an eppi0 candidate through the reconstructed
    electron/proton parent kinematics:

        P_probe^reco = P_pi0^reco - k_tag^reco
        P_probe^pred = P_beam + P_target - P_e - P_p - k_tag^reco

    Mh_gammagamma supplies the reconstructed pi0 mass event by event.

    This routine is vectorized because it is applied to O(10^5-10^6) pairs per
    period in full-statistics production.
    """
    i = np.asarray(matches.epg_index, dtype=np.int64)
    j = np.asarray(matches.eppi0_index, dtype=np.int64)
    n = len(i)

    if n == 0:
        return {}, {
            "matched_parent_pairs": 0,
            "invalid_pi0_mass": 0,
            "nonphysical_reco_probe_energy": 0,
            "nonphysical_predicted_probe": 0,
            "accepted_stage1_pairs": 0,
        }
    #endif

    e3 = np.asarray(epg.electron_p3[i], dtype=float)
    p3 = np.asarray(epg.proton_p3[i], dtype=float)
    tag3 = np.asarray(epg.tag_p3[i], dtype=float)
    pi3 = np.asarray(eppi0.pi0_p3[j], dtype=float)

    e_pmag = np.linalg.norm(e3, axis=1)
    p_pmag = np.linalg.norm(p3, axis=1)
    tag_E = np.linalg.norm(tag3, axis=1)
    pi_pmag = np.linalg.norm(pi3, axis=1)
    pi_mass = np.asarray(eppi0.pi0_mass[j], dtype=float)

    e_E = np.sqrt(e_pmag * e_pmag + M_E * M_E)
    p_E = np.sqrt(p_pmag * p_pmag + M_P * M_P)
    pi_E = np.sqrt(np.maximum(0.0, pi_pmag * pi_pmag + pi_mass * pi_mass))

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    initial_E = period.beam_energy + M_P
    initial_p3 = np.asarray([0.0, 0.0, beam_p], dtype=float)

    # Reconstructed companion photon from the reconstructed pi0 and tag.
    reco_E = pi_E - tag_E
    reco3 = pi3 - tag3
    reco_p = np.linalg.norm(reco3, axis=1)
    reco_m2 = reco_E * reco_E - reco_p * reco_p
    reco_E_minus_p = reco_E - reco_p

    # Missing/probe photon predicted independently of the eppi0 reconstruction.
    pred_E = initial_E - e_E - p_E - tag_E
    pred3 = initial_p3[None, :] - e3 - p3 - tag3
    pred_p = np.linalg.norm(pred3, axis=1)
    pred_m2 = pred_E * pred_E - pred_p * pred_p
    pred_E_minus_p = pred_E - pred_p

    # Complete reconstructed e p pi0 exclusivity closure.  This is diagnostic
    # only and is deliberately NOT used to define the clean tag association.
    miss_E = initial_E - e_E - p_E - pi_E
    miss3 = initial_p3[None, :] - e3 - p3 - pi3
    miss_p = np.linalg.norm(miss3, axis=1)
    miss_pT = np.hypot(miss3[:, 0], miss3[:, 1])
    miss_m2 = miss_E * miss_E - miss_p * miss_p

    # Angular quantities.
    with np.errstate(invalid="ignore", divide="ignore"):
        pred_unit = pred3 / pred_p[:, None]
        reco_unit = reco3 / reco_p[:, None]
        cos_open = np.sum(pred_unit * reco_unit, axis=1)
        cos_open = np.clip(cos_open, -1.0, 1.0)
        opening_deg = np.degrees(np.arccos(cos_open))

        pred_theta_deg = np.degrees(
            np.arccos(np.clip(pred_unit[:, 2], -1.0, 1.0))
        )
        reco_theta_deg = np.degrees(
            np.arccos(np.clip(reco_unit[:, 2], -1.0, 1.0))
        )
        pred_phi_deg = np.degrees(np.arctan2(pred3[:, 1], pred3[:, 0]))
        reco_phi_deg = np.degrees(np.arctan2(reco3[:, 1], reco3[:, 0]))

        delta_theta_deg = reco_theta_deg - pred_theta_deg
        delta_phi_deg = (
            (reco_phi_deg - pred_phi_deg + 180.0) % 360.0
        ) - 180.0

        delta_E = reco_E - pred_E
        frac_delta_E = delta_E / pred_E
    #endwith

    # Exact CLAS12 FD sector convention supplied for this analysis:
    #   S1: [330,360) U [0,30)
    #   S2: [ 30, 90)
    #   S3: [ 90,150)
    #   S4: [150,210)
    #   S5: [210,270)
    #   S6: [270,330)
    #
    # Convert predicted lab phi to [0,360) and assign the wrapped ranges
    # explicitly rather than relying on an approximate sector-center formula.
    pred_phi_360 = pred_phi_deg % 360.0
    pred_sector = np.full(len(pred_phi_360), -1, dtype=np.int16)

    pred_sector[
        (pred_phi_360 >= 330.0) | (pred_phi_360 < 30.0)
    ] = 1
    pred_sector[
        (pred_phi_360 >= 30.0) & (pred_phi_360 < 90.0)
    ] = 2
    pred_sector[
        (pred_phi_360 >= 90.0) & (pred_phi_360 < 150.0)
    ] = 3
    pred_sector[
        (pred_phi_360 >= 150.0) & (pred_phi_360 < 210.0)
    ] = 4
    pred_sector[
        (pred_phi_360 >= 210.0) & (pred_phi_360 < 270.0)
    ] = 5
    pred_sector[
        (pred_phi_360 >= 270.0) & (pred_phi_360 < 330.0)
    ] = 6

    detector_tag = np.full(n, -999, dtype=np.int16)
    if "detector2" in epg.raw:
        detector_tag = np.asarray(epg.raw["detector2"][i], dtype=np.int16)
    #endif

    arrays = {
        "epg_index": i,
        "eppi0_index": j,
        "parent_distance": np.asarray(matches.nearest_distance, dtype=float),
        "parent_max_component_delta": np.asarray(matches.max_component_delta, dtype=float),
        "tag_energy": tag_E,
        "pi0_mass": pi_mass,
        "pi0_energy": pi_E,
        "reco_probe_energy": reco_E,
        "reco_probe_p": reco_p,
        "reco_probe_mass2": reco_m2,
        "reco_probe_E_minus_p": reco_E_minus_p,
        "pred_probe_energy": pred_E,
        "pred_probe_p": pred_p,
        "pred_probe_mass2": pred_m2,
        "pred_probe_E_minus_p": pred_E_minus_p,
        "probe_delta_E": delta_E,
        "probe_delta_E_over_E": frac_delta_E,
        "probe_delta_theta_deg": delta_theta_deg,
        "probe_delta_phi_deg": delta_phi_deg,
        "probe_opening_residual_deg": opening_deg,
        "pred_probe_theta_deg": pred_theta_deg,
        "pred_probe_phi_deg": pred_phi_deg,
        "pred_probe_sector": pred_sector,
        "reco_probe_theta_deg": reco_theta_deg,
        "reco_probe_phi_deg": reco_phi_deg,
        "exclusivity_missing_energy": miss_E,
        "exclusivity_missing_p": miss_p,
        "exclusivity_missing_pT": miss_pT,
        "exclusivity_missing_mass2": miss_m2,
        "detector_tag_epgamma": detector_tag,
    }

    finite_common = (
        np.isfinite(pi_mass)
        & np.isfinite(reco_E)
        & np.isfinite(reco_p)
        & np.isfinite(pred_E)
        & np.isfinite(pred_p)
    )
    arrays = {name: value[finite_common] for name, value in arrays.items()}

    counters = {
        "matched_parent_pairs": int(n),
        "invalid_pi0_mass": int(np.count_nonzero(~np.isfinite(pi_mass) | (pi_mass <= 0.0))),
        "nonphysical_reco_probe_energy": int(np.count_nonzero(~np.isfinite(reco_E) | (reco_E <= 0.0))),
        "nonphysical_predicted_probe": int(np.count_nonzero(~np.isfinite(pred_E) | (pred_p <= 0.0))),
        "accepted_stage1_pairs": int(np.count_nonzero(finite_common)),
    }
    return arrays, counters


def build_clean_association_mask(
    arrays: Dict[str, np.ndarray],
    mgg_min: float,
    mgg_max: float,
    remainder_mass2_max: float,
    reco_probe_energy_min: float,
) -> np.ndarray:
    """
    Define a clean reconstructed-pi0 tag association.

    Crucially, this uses only reconstructed-side quantities.  It does NOT cut on
    the predicted probe, the predicted/reconstructed angular residual, or the
    ep-pi0 exclusivity residual.  Therefore the subsequent probe-resolution
    study is not artificially narrowed by its own selection.

    Defaults:
      0.10 < M_gg < 0.17 GeV
      |(P_pi0 - k_tag)^2| < 1e-3 GeV^2
      E_probe^reco >= 0.40 GeV
    """
    return (
        np.isfinite(arrays["pi0_mass"])
        & (arrays["pi0_mass"] >= mgg_min)
        & (arrays["pi0_mass"] <= mgg_max)
        & np.isfinite(arrays["reco_probe_mass2"])
        & (np.abs(arrays["reco_probe_mass2"]) <= remainder_mass2_max)
        & np.isfinite(arrays["reco_probe_energy"])
        & (arrays["reco_probe_energy"] >= reco_probe_energy_min)
        & (arrays["reco_probe_p"] > 0.0)
    )


def quantile_or_nan(values: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    #endif
    return float(np.quantile(values, q))


def resolution_rows(
    period: PeriodConfig,
    arrays: Dict[str, np.ndarray],
    clean: np.ndarray,
    ft_theta_max: float,
    min_count: int,
) -> List[Dict[str, object]]:
    """
    Measure probe prediction resolution versus predicted energy and detector
    region.  These are diagnostics/candidate matching radii, not final cuts.

    Region assignment is based on the PREDICTED probe direction so that the same
    definition can later be applied when the probe is not reconstructed:
      FT-like: theta_pred <= ft_theta_max
      FD S1..S6: theta_pred > ft_theta_max, sector from exact wrapped phi ranges

    The exact FT/FD geometric/fiducial acceptance will be refined later.
    """
    E = arrays["pred_probe_energy"]
    opening = arrays["probe_opening_residual_deg"]
    dE = arrays["probe_delta_E"]
    frac = arrays["probe_delta_E_over_E"]
    theta = arrays["pred_probe_theta_deg"]
    sector = arrays["pred_probe_sector"]

    region_masks: List[Tuple[str, np.ndarray]] = [
        ("all", np.ones(len(E), dtype=bool)),
        ("FT_like", theta <= ft_theta_max),
    ]
    for s in range(1, 7):
        region_masks.append(
            (f"FD_S{s}", (theta > ft_theta_max) & (sector == s))
        )
    #endfor

    rows: List[Dict[str, object]] = []
    for region_name, region_mask in region_masks:
        for ibin in range(len(PROBE_ENERGY_EDGES) - 1):
            lo = float(PROBE_ENERGY_EDGES[ibin])
            hi = float(PROBE_ENERGY_EDGES[ibin + 1])
            mask = (
                clean
                & region_mask
                & np.isfinite(E)
                & (E >= lo)
                & (E < hi)
                & np.isfinite(opening)
            )
            n = int(np.count_nonzero(mask))

            row: Dict[str, object] = {
                "period": period.key,
                "label": period.label,
                "region": region_name,
                "energy_low_GeV": lo,
                "energy_high_GeV": hi,
                "energy_center_GeV": 0.5 * (lo + hi),
                "count": n,
            }

            if n >= min_count:
                a = opening[mask]
                de = dE[mask]
                fr = frac[mask]
                row.update({
                    "angle_q50_deg": quantile_or_nan(a, 0.50),
                    "angle_q68_deg": quantile_or_nan(a, 0.68),
                    "angle_q90_deg": quantile_or_nan(a, 0.90),
                    "angle_q95_deg": quantile_or_nan(a, 0.95),
                    "angle_q99_deg": quantile_or_nan(a, 0.99),
                    "deltaE_median_GeV": quantile_or_nan(de, 0.50),
                    "deltaE_q16_GeV": quantile_or_nan(de, 0.16),
                    "deltaE_q84_GeV": quantile_or_nan(de, 0.84),
                    "frac_deltaE_median": quantile_or_nan(fr, 0.50),
                    "frac_deltaE_q16": quantile_or_nan(fr, 0.16),
                    "frac_deltaE_q84": quantile_or_nan(fr, 0.84),
                })
            else:
                row.update({
                    "angle_q50_deg": float("nan"),
                    "angle_q68_deg": float("nan"),
                    "angle_q90_deg": float("nan"),
                    "angle_q95_deg": float("nan"),
                    "angle_q99_deg": float("nan"),
                    "deltaE_median_GeV": float("nan"),
                    "deltaE_q16_GeV": float("nan"),
                    "deltaE_q84_GeV": float("nan"),
                    "frac_deltaE_median": float("nan"),
                    "frac_deltaE_q16": float("nan"),
                    "frac_deltaE_q84": float("nan"),
                })
            #endif

            rows.append(row)
        #endfor
    #endfor

    return rows


def write_rows_csv(rows: List[Dict[str, object]], path: Path) -> None:
    """
    Write a list of dictionaries whose schemas may legitimately differ.

    Some diagnostic columns (for example closure-test fields) exist only for
    selected rows.  Therefore deriving fieldnames from rows[0] is unsafe and
    can fail late, after all expensive physics processing has completed.

    Column order is deterministic: first occurrence of each key while scanning
    rows in their current order.
    """
    if not rows:
        return
    #endif

    fieldnames: List[str] = []
    seen: set[str] = set()

    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
            #endif
        #endfor
    #endfor

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            extrasaction="ignore",
            restval="",
        )
        writer.writeheader()
        writer.writerows(rows)
    #endwith

# -----------------------------------------------------------------------------
# Stage-II: real epgamma denominator composition and discriminator stability
# -----------------------------------------------------------------------------

STAGE2_DISCRIMINATORS: Tuple[str, ...] = (
    "mx2_ep_1d",
    "probe_m2_1d",
    "mx2_ep_x_probe_m2",
    "mx2_ep_x_pTmiss",
    "mx2_ep_x_theta_gg",
)

STAGE2_NOMINAL_DISCRIMINATOR = "mx2_ep_x_probe_m2"
STAGE2_PRODUCTION_MODEL = "shared_morphed_2d"
STAGE2_DRIVER_DISCRIMINATORS: Tuple[str, ...] = (
    "mx2_ep_x_probe_m2",
    "mx2_ep_x_pTmiss",
)
STAGE2_VALIDATION_DISCRIMINATORS: Tuple[str, ...] = (
    "mx2_ep_x_theta_gg",
    "mx2_ep_1d",
    "probe_m2_1d",
)


def stage2_energy_edges(max_energy: float) -> np.ndarray:
    """
    Return Stage-II probe-energy edges, truncated at the supported data endpoint.
    """
    base = PROBE_ENERGY_EDGES[PROBE_ENERGY_EDGES < max_energy]
    edges = np.concatenate((base, np.asarray([max_energy], dtype=float)))
    edges = np.unique(edges)
    return edges[edges >= PROBE_ENERGY_EDGES[0]]


def parse_mix_shifts(value: str) -> List[int]:
    shifts: List[int] = []
    for token in str(value).split(","):
        token = token.strip()
        if not token:
            continue
        #endif
        shift = int(token)
        if shift == 0:
            raise ValueError("--mix-shifts may not contain zero.")
        #endif
        shifts.append(shift)
    #endfor
    if not shifts:
        raise ValueError("--mix-shifts must contain at least one nonzero integer.")
    #endif
    return shifts


def predicted_region_arrays(
    pred_theta_deg: np.ndarray,
    pred_phi_deg: np.ndarray,
    ft_theta_max: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (sector, region_code) from the predicted photon direction.

    region_code:
      0 = outside/unclassified
      1 = FT, theta <= ft_theta_max
      2..7 = FD sectors 1..6, theta > ft_theta_max
    """
    phi360 = np.asarray(pred_phi_deg, dtype=float) % 360.0
    theta = np.asarray(pred_theta_deg, dtype=float)

    sector = np.full(len(theta), -1, dtype=np.int16)
    sector[(phi360 >= 330.0) | (phi360 < 30.0)] = 1
    sector[(phi360 >= 30.0) & (phi360 < 90.0)] = 2
    sector[(phi360 >= 90.0) & (phi360 < 150.0)] = 3
    sector[(phi360 >= 150.0) & (phi360 < 210.0)] = 4
    sector[(phi360 >= 210.0) & (phi360 < 270.0)] = 5
    sector[(phi360 >= 270.0) & (phi360 < 330.0)] = 6

    region = np.zeros(len(theta), dtype=np.int16)
    finite = np.isfinite(theta) & np.isfinite(phi360)
    region[finite & (theta <= ft_theta_max)] = 1
    fd = finite & (theta > ft_theta_max) & (sector >= 1)
    region[fd] = sector[fd] + 1

    return sector, region


def build_epgamma_denominator_features(
    period: PeriodConfig,
    epg: EPGammaSample,
    tag_min: float,
    tag_max: float,
    ft_theta_max: float,
    mixed_tag_shift: int = 0,
) -> Dict[str, np.ndarray]:
    """
    Same Stage-II kinematics with lower peak memory and a float32 persistent
    histogram feature store.
    """
    e3 = epg.electron_p3
    p3 = epg.proton_p3
    tag_source = epg.tag_p3
    n = len(epg.tag_energy)

    if mixed_tag_shift and n:
        shift = int(mixed_tag_shift) % n
        if shift == 0:
            shift = max(1, n // 2)
        #endif
        tag_index = (np.arange(n, dtype=np.int64) - shift) % n
        tag_E = epg.tag_energy[tag_index]
    else:
        tag_index = None
        tag_E = epg.tag_energy
    #endif

    e_E = np.sqrt(np.einsum("ij,ij->i", e3, e3) + M_E * M_E)
    p_E = np.sqrt(np.einsum("ij,ij->i", p3, p3) + M_P * M_P)

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    initial_E = period.beam_energy + M_P
    epmiss_E = initial_E - e_E - p_E
    del e_E, p_E

    px = -e3[:, 0] - p3[:, 0]
    py = -e3[:, 1] - p3[:, 1]
    pz = beam_p - e3[:, 2] - p3[:, 2]

    epmiss_m2 = epmiss_E * epmiss_E - (px * px + py * py + pz * pz)
    pred_E = epmiss_E - tag_E

    if tag_index is None:
        px -= tag_source[:, 0]
        py -= tag_source[:, 1]
        pz -= tag_source[:, 2]
    else:
        px -= tag_source[tag_index, 0]
        py -= tag_source[tag_index, 1]
        pz -= tag_source[tag_index, 2]
    #endif

    pred_p2 = px * px + py * py + pz * pz
    pred_p = np.sqrt(pred_p2)
    pred_m2 = pred_E * pred_E - pred_p2

    with np.errstate(invalid="ignore", divide="ignore"):
        pred_theta_deg = np.degrees(
            np.arccos(
                np.clip(
                    pz / np.where(pred_p > 0.0, pred_p, np.nan),
                    -1.0,
                    1.0,
                )
            )
        )
        pred_phi_deg = np.degrees(np.arctan2(py, px))
    #endwith

    sector, _region = predicted_region_arrays(
        pred_theta_deg,
        pred_phi_deg,
        ft_theta_max=ft_theta_max,
    )

    valid = (
        np.isfinite(tag_E)
        & (tag_E >= tag_min)
        & (tag_E < tag_max)
        & np.isfinite(epmiss_m2)
        & np.isfinite(pred_E)
        & np.isfinite(pred_m2)
        & np.isfinite(pred_theta_deg)
        & np.isfinite(pred_phi_deg)
    )

    out = {
        "ep_missing_mass2": np.asarray(epmiss_m2, dtype=np.float32),
        "pred_probe_energy": np.asarray(pred_E, dtype=np.float32),
        "pred_probe_mass2": np.asarray(pred_m2, dtype=np.float32),
        "pred_probe_theta_deg": np.asarray(pred_theta_deg, dtype=np.float32),
        "pred_probe_sector": np.asarray(sector, dtype=np.int16),
        "valid_tag": np.asarray(valid, dtype=bool),
    }

    if not mixed_tag_shift:
        for raw_name, out_name in (
            ("pTmiss", "stored_pTmiss"),
            ("theta_gamma_gamma", "stored_theta_gamma_gamma"),
        ):
            if raw_name in epg.raw:
                out[out_name] = np.asarray(
                    epg.raw[raw_name], dtype=np.float32
                )
            #endif
        #endfor
    else:
        # For the mixed-event wrong-photon template the skim-level pTmiss
        # belongs to the original photon and is therefore invalid. Recompute
        # the epgamma missing transverse momentum after substituting the
        # deterministically shifted photon. This makes the mixed template
        # usable in the same pTmiss production driver without borrowing the
        # original event's correlated value.
        out["stored_pTmiss"] = np.asarray(
            np.sqrt(px * px + py * py),
            dtype=np.float32,
        )
    #endif

    if "runnum" in epg.raw:
        out["runnum"] = np.asarray(epg.raw["runnum"])
    #endif
    if "evnum" in epg.raw:
        out["evnum"] = np.asarray(epg.raw["evnum"])
    #endif

    return out


def stage2_region_mask(
    features: Dict[str, np.ndarray],
    region_name: str,
    ft_theta_max: float,
) -> np.ndarray:
    theta = features["pred_probe_theta_deg"]
    sector = features["pred_probe_sector"]

    if region_name == "all":
        return np.ones(len(theta), dtype=bool)
    #endif
    if region_name == "FT":
        return theta <= ft_theta_max
    #endif
    if region_name == "FD_all":
        return (theta > ft_theta_max) & (sector >= 1) & (sector <= 6)
    #endif
    if region_name.startswith("FD_S"):
        s = int(region_name[-1])
        return (theta > ft_theta_max) & (sector == s)
    #endif
    raise ValueError(f"Unknown Stage-II region '{region_name}'")


def stage2_fit_mask(
    feat: Dict[str, np.ndarray],
    region_name: str,
    ft_theta_max: float,
    elo: float,
    ehi: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
) -> np.ndarray:
    return (
        feat["valid_tag"]
        & stage2_region_mask(feat, region_name, ft_theta_max)
        & (feat["pred_probe_energy"] >= elo)
        & (feat["pred_probe_energy"] < ehi)
        & (feat["ep_missing_mass2"] >= mm2_min)
        & (feat["ep_missing_mass2"] < mm2_max)
        & (np.abs(feat["pred_probe_mass2"]) < probe_m2_max)
    )


def discriminator_available(
    discriminator: str,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
) -> bool:
    if discriminator == "mx2_ep_x_pTmiss":
        return all("stored_pTmiss" in f for f in (data_f, pi0_f, dvcs_f))
    #endif
    if discriminator == "mx2_ep_x_theta_gg":
        return all("stored_theta_gamma_gamma" in f for f in (data_f, pi0_f, dvcs_f))
    #endif
    return True


def histogram_for_discriminator(
    discriminator: str,
    feat: Dict[str, np.ndarray],
    mask: np.ndarray,
    *,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    bins_1d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
) -> np.ndarray:
    """
    Build one histogram for the requested discriminator.

    The same baseline energy/region/MM2/probe-M2 support mask is applied before
    constructing every discriminator so the extracted fractions are compared on
    as nearly identical an event sample as possible.
    """
    if discriminator == "mx2_ep_1d":
        h, _ = np.histogram(
            feat["ep_missing_mass2"][mask],
            bins=np.linspace(mm2_min, mm2_max, bins_1d + 1),
        )
        return h.astype(float)
    #endif

    if discriminator == "probe_m2_1d":
        h, _ = np.histogram(
            feat["pred_probe_mass2"][mask],
            bins=np.linspace(-probe_m2_max, probe_m2_max, bins_1d + 1),
        )
        return h.astype(float)
    #endif

    if discriminator == "mx2_ep_x_probe_m2":
        h, _, _ = np.histogram2d(
            feat["ep_missing_mass2"][mask],
            feat["pred_probe_mass2"][mask],
            bins=(
                np.linspace(mm2_min, mm2_max, mm2_bins_2d + 1),
                np.linspace(-probe_m2_max, probe_m2_max, probe_m2_bins_2d + 1),
            ),
        )
        return h.astype(float)
    #endif

    if discriminator == "mx2_ep_x_pTmiss":
        local = mask & np.isfinite(feat["stored_pTmiss"])
        h, _, _ = np.histogram2d(
            feat["ep_missing_mass2"][local],
            feat["stored_pTmiss"][local],
            bins=(
                np.linspace(mm2_min, mm2_max, mm2_bins_2d + 1),
                np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            ),
        )
        return h.astype(float)
    #endif

    if discriminator == "mx2_ep_x_theta_gg":
        local = mask & np.isfinite(feat["stored_theta_gamma_gamma"])
        h, _, _ = np.histogram2d(
            feat["ep_missing_mass2"][local],
            feat["stored_theta_gamma_gamma"][local],
            bins=(
                np.linspace(mm2_min, mm2_max, mm2_bins_2d + 1),
                np.linspace(0.0, theta_max, theta_bins + 1),
            ),
        )
        return h.astype(float)
    #endif

    raise ValueError(f"Unknown Stage-II discriminator '{discriminator}'")


@dataclass
class TwoComponentFitResult:
    success: bool
    message: str
    data_count: int
    aaogen_count: int
    dvcsgen_count: int
    pi0_yield: float = float("nan")
    dvcs_yield: float = float("nan")
    pi0_fraction: float = float("nan")
    dvcs_fraction: float = float("nan")
    nll: float = float("nan")
    poisson_deviance: float = float("nan")
    ndof: int = 0


@dataclass
class ThreeComponentDiagnosticResult:
    success: bool
    pi0_fraction: float = float("nan")
    dvcs_fraction: float = float("nan")
    mixed_fraction: float = float("nan")
    pi0_yield: float = float("nan")
    mixed_yield: float = float("nan")


def normalized_template(hist: np.ndarray) -> np.ndarray:
    x = np.asarray(hist, dtype=float).ravel()
    eps = 1.0e-12
    return (x + eps) / np.sum(x + eps)


def poisson_deviance_quality(
    data: np.ndarray,
    mu: np.ndarray,
    n_parameters: int,
) -> Tuple[float, int]:
    data = np.asarray(data, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1.0e-12, None)

    positive = data > 0.0
    terms = np.zeros_like(data, dtype=float)
    terms[positive] = data[positive] * np.log(data[positive] / mu[positive])
    deviance = 2.0 * float(np.sum(mu - data + terms))

    populated = int(np.count_nonzero((data > 0.0) | (mu > 1.0e-10)))
    ndof = max(1, populated - n_parameters)
    return deviance, ndof


def fit_two_component_poisson(
    data_hist: np.ndarray,
    pi0_hist: np.ndarray,
    dvcs_hist: np.ndarray,
) -> TwoComponentFitResult:
    """
    Extended binned Poisson fit:

      mu_i = N_pi0 T_pi0_i + N_DVCS T_DVCS_i.

    aaogen and dvcsgen provide local shapes only.
    """
    data = np.asarray(data_hist, dtype=float).ravel()
    p_raw = np.asarray(pi0_hist, dtype=float).ravel()
    d_raw = np.asarray(dvcs_hist, dtype=float).ravel()

    nd = int(round(float(np.sum(data))))
    np0 = int(round(float(np.sum(p_raw))))
    ndv = int(round(float(np.sum(d_raw))))

    if nd <= 0 or np0 <= 0 or ndv <= 0:
        return TwoComponentFitResult(
            False, "one or more empty histograms", nd, np0, ndv
        )
    #endif

    p = normalized_template(p_raw)
    d = normalized_template(d_raw)
    T = np.column_stack((p, d))

    def objective(yields: np.ndarray) -> float:
        mu = np.clip(T @ yields, 1.0e-12, None)
        return float(np.sum(mu - data * np.log(mu)))
    #enddef

    result = minimize(
        objective,
        np.asarray([0.70 * nd, 0.30 * nd], dtype=float),
        method="L-BFGS-B",
        bounds=((0.0, None), (0.0, None)),
        options={"maxiter": 400, "ftol": 1.0e-11},
    )

    y = np.clip(np.asarray(result.x, dtype=float), 0.0, None)
    total = float(np.sum(y))
    if total <= 0.0:
        return TwoComponentFitResult(
            False, "fit returned zero total yield", nd, np0, ndv
        )
    #endif

    mu = T @ y
    deviance, ndof = poisson_deviance_quality(data, mu, n_parameters=2)

    return TwoComponentFitResult(
        success=bool(result.success),
        message=str(result.message),
        data_count=nd,
        aaogen_count=np0,
        dvcsgen_count=ndv,
        pi0_yield=float(y[0]),
        dvcs_yield=float(y[1]),
        pi0_fraction=float(y[0] / total),
        dvcs_fraction=float(y[1] / total),
        nll=float(result.fun),
        poisson_deviance=deviance,
        ndof=ndof,
    )


def fit_three_component_diagnostic(
    data_hist: np.ndarray,
    pi0_hist: np.ndarray,
    dvcs_hist: np.ndarray,
    mixed_hist: np.ndarray,
) -> ThreeComponentDiagnosticResult:
    """
    Diagnostic only: allow one mixed-event component to float and record the
    induced pi0-fraction shift. Never used for the nominal denominator.
    """
    data = np.asarray(data_hist, dtype=float).ravel()
    if (
        np.sum(data) <= 0
        or np.sum(pi0_hist) <= 0
        or np.sum(dvcs_hist) <= 0
        or np.sum(mixed_hist) <= 0
    ):
        return ThreeComponentDiagnosticResult(False)
    #endif

    T = np.column_stack((
        normalized_template(pi0_hist),
        normalized_template(dvcs_hist),
        normalized_template(mixed_hist),
    ))
    nd = float(np.sum(data))

    def objective(yields: np.ndarray) -> float:
        mu = np.clip(T @ yields, 1.0e-12, None)
        return float(np.sum(mu - data * np.log(mu)))
    #enddef

    result = minimize(
        objective,
        np.asarray([0.70 * nd, 0.29 * nd, 0.01 * nd], dtype=float),
        method="L-BFGS-B",
        bounds=((0.0, None), (0.0, None), (0.0, None)),
        options={"maxiter": 400, "ftol": 1.0e-11},
    )
    y = np.clip(np.asarray(result.x, dtype=float), 0.0, None)
    total = float(np.sum(y))
    if total <= 0.0:
        return ThreeComponentDiagnosticResult(False)
    #endif

    return ThreeComponentDiagnosticResult(
        success=bool(result.success),
        pi0_fraction=float(y[0] / total),
        dvcs_fraction=float(y[1] / total),
        mixed_fraction=float(y[2] / total),
        pi0_yield=float(y[0]),
        mixed_yield=float(y[2]),
    )


def run_stage2_discriminator_fits(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    bins_1d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
    min_data_count: int,
    min_template_count: int,
) -> List[Dict[str, object]]:
    """
    Fit all available discriminators independently in each
    period / region / energy bin.

    The comparison among extracted pi0 fractions is the central Stage-II
    robustness test.
    """
    edges = stage2_energy_edges(max_probe_energy)
    regions = ["all", "FT", "FD_all"] + [f"FD_S{s}" for s in range(1, 7)]
    rows: List[Dict[str, object]] = []

    for region_name in regions:
        for ibin in range(len(edges) - 1):
            elo = float(edges[ibin])
            ehi = float(edges[ibin + 1])

            masks = {
                "data": stage2_fit_mask(
                    data_f, region_name, ft_theta_max, elo, ehi,
                    mm2_min, mm2_max, probe_m2_max
                ),
                "pi0": stage2_fit_mask(
                    pi0_f, region_name, ft_theta_max, elo, ehi,
                    mm2_min, mm2_max, probe_m2_max
                ),
                "dvcs": stage2_fit_mask(
                    dvcs_f, region_name, ft_theta_max, elo, ehi,
                    mm2_min, mm2_max, probe_m2_max
                ),
            }

            raw_counts = {
                key: int(np.count_nonzero(mask))
                for key, mask in masks.items()
            }

            for disc in STAGE2_DISCRIMINATORS:
                row: Dict[str, object] = {
                    "period": period.key,
                    "label": period.label,
                    "beam_energy_GeV": period.beam_energy,
                    "region": region_name,
                    "energy_low_GeV": elo,
                    "energy_high_GeV": ehi,
                    "energy_center_GeV": 0.5 * (elo + ehi),
                    "discriminator": disc,
                    "is_nominal_discriminator": int(
                        disc == STAGE2_NOMINAL_DISCRIMINATOR
                    ),
                    "data_count_pre_discriminator": raw_counts["data"],
                    "aaogen_count_pre_discriminator": raw_counts["pi0"],
                    "dvcsgen_count_pre_discriminator": raw_counts["dvcs"],
                }

                if not discriminator_available(disc, data_f, pi0_f, dvcs_f):
                    row.update({
                        "fit_success": 0,
                        "fit_message": "required diagnostic branch unavailable",
                        "data_count": 0,
                        "aaogen_count": 0,
                        "dvcsgen_count": 0,
                        "pi0_yield": float("nan"),
                        "dvcs_yield": float("nan"),
                        "pi0_fraction": float("nan"),
                        "dvcs_fraction": float("nan"),
                        "nll": float("nan"),
                        "poisson_deviance": float("nan"),
                        "ndof": 0,
                        "deviance_per_ndof": float("nan"),
                    })
                    rows.append(row)
                    continue
                #endif

                hd = histogram_for_discriminator(
                    disc, data_f, masks["data"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    probe_m2_max=probe_m2_max,
                    mm2_bins_2d=mm2_bins_2d,
                    probe_m2_bins_2d=probe_m2_bins_2d,
                    bins_1d=bins_1d,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )
                hp = histogram_for_discriminator(
                    disc, pi0_f, masks["pi0"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    probe_m2_max=probe_m2_max,
                    mm2_bins_2d=mm2_bins_2d,
                    probe_m2_bins_2d=probe_m2_bins_2d,
                    bins_1d=bins_1d,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )
                hv = histogram_for_discriminator(
                    disc, dvcs_f, masks["dvcs"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    probe_m2_max=probe_m2_max,
                    mm2_bins_2d=mm2_bins_2d,
                    probe_m2_bins_2d=probe_m2_bins_2d,
                    bins_1d=bins_1d,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )

                counts = {
                    "data": int(round(float(np.sum(hd)))),
                    "pi0": int(round(float(np.sum(hp)))),
                    "dvcs": int(round(float(np.sum(hv)))),
                }
                row.update({
                    "data_count": counts["data"],
                    "aaogen_count": counts["pi0"],
                    "dvcsgen_count": counts["dvcs"],
                })

                enough = (
                    counts["data"] >= min_data_count
                    and counts["pi0"] >= min_template_count
                    and counts["dvcs"] >= min_template_count
                )
                if enough:
                    fit = fit_two_component_poisson(hd, hp, hv)
                    row.update({
                        "fit_success": int(fit.success),
                        "fit_message": fit.message,
                        "pi0_yield": fit.pi0_yield,
                        "dvcs_yield": fit.dvcs_yield,
                        "pi0_fraction": fit.pi0_fraction,
                        "dvcs_fraction": fit.dvcs_fraction,
                        "nll": fit.nll,
                        "poisson_deviance": fit.poisson_deviance,
                        "ndof": fit.ndof,
                        "deviance_per_ndof": (
                            fit.poisson_deviance / fit.ndof
                            if fit.ndof > 0 else float("nan")
                        ),
                    })
                else:
                    row.update({
                        "fit_success": 0,
                        "fit_message": "insufficient data/template statistics",
                        "pi0_yield": float("nan"),
                        "dvcs_yield": float("nan"),
                        "pi0_fraction": float("nan"),
                        "dvcs_fraction": float("nan"),
                        "nll": float("nan"),
                        "poisson_deviance": float("nan"),
                        "ndof": 0,
                        "deviance_per_ndof": float("nan"),
                    })
                #endif

                rows.append(row)
            #endfor
        #endfor
    #endfor

    return rows



@dataclass
class SharedMorphedFitResult:
    success: bool
    message: str
    pi0_fraction: float = float("nan")
    objective: float = float("nan")
    poisson_deviance: float = float("nan")
    ndof: int = 0
    nuisance: Optional[Dict[str, float]] = None


def morph_template_second_axis(hist: np.ndarray, shift_bins: float, sigma_bins: float) -> np.ndarray:
    """Shift and additionally smear a 2D template along its second discriminator axis."""
    h = np.asarray(hist, dtype=float)
    if h.ndim != 2:
        raise ValueError("Template morphing expects a 2D histogram.")
    out = h.copy()
    if sigma_bins > 1.0e-8:
        out = gaussian_filter1d(out, sigma=float(sigma_bins), axis=1, mode="nearest")
    #endif
    if abs(shift_bins) > 1.0e-8:
        out = ndimage_shift(out, shift=(0.0, float(shift_bins)), order=1, mode="nearest", prefilter=False)
    #endif
    out = np.clip(out, 0.0, None)
    total = float(np.sum(out))
    if total <= 0.0:
        return normalized_template(h).reshape(h.shape)
    #endif
    return out / total


def fit_control_shift_smear(data_values: np.ndarray, mc_values: np.ndarray, lo: float, hi: float, bins: int) -> Dict[str, float]:
    """Calibrate a modest MC shift/additional-smearing from a clean control projection."""
    d = np.asarray(data_values, dtype=float)
    m = np.asarray(mc_values, dtype=float)
    d = d[np.isfinite(d) & (d >= lo) & (d < hi)]
    m = m[np.isfinite(m) & (m >= lo) & (m < hi)]
    if len(d) < 500 or len(m) < 500:
        return {"success": 0, "shift_bins": 0.0, "sigma_bins": 0.0, "deviance_per_ndof": float("nan")}
    #endif
    hd, _ = np.histogram(d, bins=np.linspace(lo, hi, bins + 1))
    hm, _ = np.histogram(m, bins=np.linspace(lo, hi, bins + 1))
    hd = hd.astype(float)
    hm2 = hm.astype(float)[None, :]
    nd = float(np.sum(hd))

    def objective(x: np.ndarray) -> float:
        t = morph_template_second_axis(hm2, x[0], x[1]).ravel()
        mu = np.clip(nd * t, 1.0e-12, None)
        return float(np.sum(mu - hd * np.log(mu)))
    #enddef

    res = minimize(objective, np.asarray([0.0, 0.5]), method="L-BFGS-B",
                   bounds=((-4.0, 4.0), (0.0, 4.0)), options={"maxiter": 250, "ftol": 1e-10})
    t = morph_template_second_axis(hm2, res.x[0], res.x[1]).ravel()
    dev, ndof = poisson_deviance_quality(hd, nd * t, 2)
    return {"success": int(res.success), "shift_bins": float(res.x[0]), "sigma_bins": float(res.x[1]),
            "deviance_per_ndof": float(dev / ndof), "data_count": int(np.sum(hd)), "mc_count": int(np.sum(hm))}


def robust_control_summary(values: np.ndarray) -> Dict[str, float]:
    """Robust location/width/tail summary for one control distribution."""
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {
            "count": 0,
            "q16": float("nan"),
            "q50": float("nan"),
            "q84": float("nan"),
            "sigma_core": float("nan"),
            "tail_2sigma_self": float("nan"),
            "tail_3sigma_self": float("nan"),
        }
    #endif
    q16, q50, q84 = np.quantile(x, [0.16, 0.50, 0.84])
    sigma = 0.5 * float(q84 - q16)
    if sigma > 0.0:
        tail2 = float(np.mean(np.abs(x - q50) > 2.0 * sigma))
        tail3 = float(np.mean(np.abs(x - q50) > 3.0 * sigma))
    else:
        tail2 = tail3 = float("nan")
    #endif
    return {
        "count": int(x.size),
        "q16": float(q16),
        "q50": float(q50),
        "q84": float(q84),
        "sigma_core": sigma,
        "tail_2sigma_self": tail2,
        "tail_3sigma_self": tail3,
    }


def parse_float_edges(value: str, name: str) -> np.ndarray:
    vals = [float(tok.strip()) for tok in str(value).split(",") if tok.strip()]
    if len(vals) < 2:
        raise ValueError(f"{name} requires at least two comma-separated edges.")
    #endif
    arr = np.asarray(vals, dtype=float)
    if not np.all(np.isfinite(arr)) or np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} edges must be finite and strictly increasing.")
    #endif
    return arr


def eppi0_control_base_arrays(sample: EPPi0Sample) -> Dict[str, np.ndarray]:
    pi0_energy = np.sqrt(np.maximum(0.0, sample.pi0_p * sample.pi0_p + sample.pi0_mass * sample.pi0_mass))
    out = {
        "mgg": np.asarray(sample.pi0_mass, dtype=float),
        "pi0_energy": np.asarray(pi0_energy, dtype=float),
        "pi0_theta_deg": np.degrees(np.asarray(sample.pi0_theta, dtype=float)),
    }
    for raw_name in ("pTmiss", "Emiss2", "theta_pi0_pi0"):
        if raw_name in sample.raw:
            out[raw_name] = np.asarray(sample.raw[raw_name], dtype=float)
        #endif
    #endfor
    for raw_name in ("detector_gamma1", "detector_gamma2"):
        if raw_name in sample.raw:
            out[raw_name] = np.asarray(sample.raw[raw_name], dtype=np.int16)
        #endif
    #endfor
    return out


def eppi0_control_selection_masks(
    arrays: Dict[str, np.ndarray],
    mgg_min: float,
    mgg_max: float,
    emiss_abs_max: float,
) -> Dict[str, np.ndarray]:
    finite_mass = np.isfinite(arrays["mgg"])
    baseline = finite_mass
    mass = baseline & (arrays["mgg"] >= mgg_min) & (arrays["mgg"] <= mgg_max)
    final = mass.copy()
    if "Emiss2" in arrays:
        final &= np.isfinite(arrays["Emiss2"]) & (np.abs(arrays["Emiss2"]) <= emiss_abs_max)
    #endif
    return {
        "baseline": baseline,
        "mass_window": mass,
        "mass_plus_emiss": final,
    }


def control_ptmiss_row(
    period: PeriodConfig,
    selection: str,
    group: str,
    data_values: np.ndarray,
    mc_values: np.ndarray,
) -> Dict[str, object]:
    d = robust_control_summary(data_values)
    m = robust_control_summary(mc_values)
    sigma_add = float("nan")
    if np.isfinite(d["sigma_core"]) and np.isfinite(m["sigma_core"]):
        sigma_add = math.sqrt(max(0.0, d["sigma_core"]**2 - m["sigma_core"]**2))
    #endif

    # Tail fractions against the MC core are especially useful: they separate
    # a genuine resolution-width difference from an excess data tail.
    mc_tail2_data = mc_tail3_data = mc_tail2_mc = mc_tail3_mc = float("nan")
    if np.isfinite(m["q50"]) and np.isfinite(m["sigma_core"]) and m["sigma_core"] > 0.0:
        dv = np.asarray(data_values, dtype=float); dv = dv[np.isfinite(dv)]
        mv = np.asarray(mc_values, dtype=float); mv = mv[np.isfinite(mv)]
        mc_tail2_data = float(np.mean(np.abs(dv - m["q50"]) > 2.0*m["sigma_core"])) if dv.size else float("nan")
        mc_tail3_data = float(np.mean(np.abs(dv - m["q50"]) > 3.0*m["sigma_core"])) if dv.size else float("nan")
        mc_tail2_mc = float(np.mean(np.abs(mv - m["q50"]) > 2.0*m["sigma_core"])) if mv.size else float("nan")
        mc_tail3_mc = float(np.mean(np.abs(mv - m["q50"]) > 3.0*m["sigma_core"])) if mv.size else float("nan")
    #endif

    return {
        "period": period.key,
        "label": period.label,
        "selection": selection,
        "group": group,
        "data_count": d["count"],
        "mc_count": m["count"],
        "data_q16_GeV": d["q16"],
        "data_q50_GeV": d["q50"],
        "data_q84_GeV": d["q84"],
        "mc_q16_GeV": m["q16"],
        "mc_q50_GeV": m["q50"],
        "mc_q84_GeV": m["q84"],
        "median_shift_data_minus_mc_GeV": d["q50"] - m["q50"] if np.isfinite(d["q50"]) and np.isfinite(m["q50"]) else float("nan"),
        "data_sigma_core_GeV": d["sigma_core"],
        "mc_sigma_core_GeV": m["sigma_core"],
        "sigma_add_core_GeV": sigma_add,
        "data_tail_2sigma_self": d["tail_2sigma_self"],
        "mc_tail_2sigma_self": m["tail_2sigma_self"],
        "data_tail_3sigma_self": d["tail_3sigma_self"],
        "mc_tail_3sigma_self": m["tail_3sigma_self"],
        "data_tail_2sigma_about_mc_core": mc_tail2_data,
        "mc_tail_2sigma_about_mc_core": mc_tail2_mc,
        "data_tail_3sigma_about_mc_core": mc_tail3_data,
        "mc_tail_3sigma_about_mc_core": mc_tail3_mc,
    }


def build_pi0_control_validation(
    period: PeriodConfig,
    eppi0_data: EPPi0Sample,
    eppi0_mc: EPPi0Sample,
    mgg_min: float,
    mgg_max: float,
    emiss_abs_max: float,
    ft_theta_max: float,
    pi0_energy_edges: np.ndarray,
    ptmiss_max: float,
    ptmiss_bins: int,
) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    """
    Reconstructed-eppi0 data/aaogen control validation.

    IMPORTANT: v031 does NOT feed this calibration into the denominator fit.
    The previous full-range Gaussian-smearing calibration saturated at its hard
    bound and failed goodness-of-fit, so it is retained only as a diagnostic.
    """
    da = eppi0_control_base_arrays(eppi0_data)
    ma = eppi0_control_base_arrays(eppi0_mc)
    ds = eppi0_control_selection_masks(da, mgg_min, mgg_max, emiss_abs_max)
    ms = eppi0_control_selection_masks(ma, mgg_min, mgg_max, emiss_abs_max)

    summary: Dict[str, object] = {
        "period": period.key,
        "label": period.label,
        "control_prior_applied_to_denominator": 0,
        "mgg_min_GeV": mgg_min,
        "mgg_max_GeV": mgg_max,
        "emiss_abs_max_GeV": emiss_abs_max,
        "pi0_energy_edges_GeV": [float(x) for x in pi0_energy_edges],
        "individual_daughter_photon_energy_available": 0,
        "energy_binning_note": "Control energy bins use reconstructed pi0 energy because daughter photon energies are not stored in the eppi0 tree.",
        "upstream_eppi0_mass_window_note": "The eppi0 skim already requires 0.11 < Mgg < 0.16 GeV; the control Mgg plot is therefore a tight diagnostic of the retained peak, not an independent cleaning cut.",
    }

    rows: List[Dict[str, object]] = []
    if "pTmiss" not in da or "pTmiss" not in ma:
        summary["pTmiss_available"] = 0
        return summary, rows
    #endif
    summary["pTmiss_available"] = 1

    # Progressive selection rows expose whether a discrepancy is background/selection driven.
    for selection in ("baseline", "mass_window", "mass_plus_emiss"):
        rows.append(control_ptmiss_row(
            period, selection, "all",
            da["pTmiss"][ds[selection]],
            ma["pTmiss"][ms[selection]],
        ))
    #endfor

    final_d = ds["mass_plus_emiss"]
    final_m = ms["mass_plus_emiss"]

    # A topology proxy based on reconstructed pi0 direction. This is NOT called FT/FD:
    # individual daughter directions are not available as simple lab theta branches.
    for group, dm, mm in (
        ("pi0_forward_theta_le_ft_boundary", final_d & (da["pi0_theta_deg"] <= ft_theta_max), final_m & (ma["pi0_theta_deg"] <= ft_theta_max)),
        ("pi0_wide_theta_gt_ft_boundary", final_d & (da["pi0_theta_deg"] > ft_theta_max), final_m & (ma["pi0_theta_deg"] > ft_theta_max)),
    ):
        rows.append(control_ptmiss_row(period, "mass_plus_emiss", group, da["pTmiss"][dm], ma["pTmiss"][mm]))
    #endfor

    for ib in range(len(pi0_energy_edges)-1):
        lo, hi = float(pi0_energy_edges[ib]), float(pi0_energy_edges[ib+1])
        dm = final_d & (da["pi0_energy"] >= lo) & (da["pi0_energy"] < hi)
        mm = final_m & (ma["pi0_energy"] >= lo) & (ma["pi0_energy"] < hi)
        rows.append(control_ptmiss_row(
            period, "mass_plus_emiss", f"pi0_energy_{lo:g}_{hi:g}_GeV",
            da["pTmiss"][dm], ma["pTmiss"][mm],
        ))
    #endfor

    # Photon detector topology. Confirmed convention:
    #   detector_gamma == 0 -> Forward Tagger (FT)
    #   detector_gamma == 1 -> Forward Detector calorimeter (FD)
    #
    # Pair codes are sorted so (FT,FD) and (FD,FT) are one topology.
    if (
        all(name in da for name in ("detector_gamma1", "detector_gamma2"))
        and all(name in ma for name in ("detector_gamma1", "detector_gamma2"))
    ):
        dp1 = np.minimum(da["detector_gamma1"], da["detector_gamma2"])
        dp2 = np.maximum(da["detector_gamma1"], da["detector_gamma2"])
        mp1 = np.minimum(ma["detector_gamma1"], ma["detector_gamma2"])
        mp2 = np.maximum(ma["detector_gamma1"], ma["detector_gamma2"])

        topology_specs = (
            (
                "FT_FT",
                final_d & (dp1 == PHOTON_DETECTOR_FT) & (dp2 == PHOTON_DETECTOR_FT),
                final_m & (mp1 == PHOTON_DETECTOR_FT) & (mp2 == PHOTON_DETECTOR_FT),
            ),
            (
                "FT_FD",
                final_d & (dp1 == PHOTON_DETECTOR_FT) & (dp2 == PHOTON_DETECTOR_FD),
                final_m & (mp1 == PHOTON_DETECTOR_FT) & (mp2 == PHOTON_DETECTOR_FD),
            ),
            (
                "FD_FD",
                final_d & (dp1 == PHOTON_DETECTOR_FD) & (dp2 == PHOTON_DETECTOR_FD),
                final_m & (mp1 == PHOTON_DETECTOR_FD) & (mp2 == PHOTON_DETECTOR_FD),
            ),
            (
                "FT_containing",
                final_d & ((dp1 == PHOTON_DETECTOR_FT) | (dp2 == PHOTON_DETECTOR_FT)),
                final_m & ((mp1 == PHOTON_DETECTOR_FT) | (mp2 == PHOTON_DETECTOR_FT)),
            ),
        )

        summary["detector_code_convention"] = {
            "0": "FT",
            "1": "FD",
        }
        summary["control_topology_groups"] = [
            "FT_FT", "FT_FD", "FD_FD", "FT_containing"
        ]

        for group, dm, mm in topology_specs:
            rows.append(
                control_ptmiss_row(
                    period,
                    "mass_plus_emiss",
                    group,
                    da["pTmiss"][dm],
                    ma["pTmiss"][mm],
                )
            )
        #endfor

        # Also retain exact observed pair-code rows as a transparent audit trail.
        pairs = sorted(
            set(zip(dp1[final_d].tolist(), dp2[final_d].tolist()))
            | set(zip(mp1[final_m].tolist(), mp2[final_m].tolist()))
        )
        summary["observed_detector_code_pairs"] = [
            [int(a), int(b)] for a, b in pairs
        ]
        for a, b in pairs:
            dm = final_d & (dp1 == a) & (dp2 == b)
            mm = final_m & (mp1 == a) & (mp2 == b)
            rows.append(
                control_ptmiss_row(
                    period,
                    "mass_plus_emiss",
                    f"detector_pair_{int(a)}_{int(b)}",
                    da["pTmiss"][dm],
                    ma["pTmiss"][mm],
                )
            )
        #endfor
    #endif

    # Retain the old global full-range morph fit strictly as a diagnostic.
    full_fit = fit_control_shift_smear(
        da["pTmiss"][final_d],
        ma["pTmiss"][final_m],
        0.0, ptmiss_max, ptmiss_bins,
    )
    full_fit["used_as_denominator_prior"] = 0
    full_fit["saturated_sigma_bound"] = int(
        np.isfinite(full_fit.get("sigma_bins", float("nan")))
        and abs(float(full_fit.get("sigma_bins", 0.0)) - 4.0) < 1.0e-6
    )
    summary["legacy_full_range_gaussian_morph_diagnostic"] = full_fit

    final_all = next((r for r in rows if r["selection"] == "mass_plus_emiss" and r["group"] == "all"), None)
    if final_all is not None:
        summary["robust_core_global"] = {
            k: v for k, v in final_all.items()
            if k not in ("period", "label", "selection", "group")
        }
    #endif
    return summary, rows


def normalized_histogram(values: np.ndarray, edges: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    h, _ = np.histogram(vals, bins=edges)
    h = h.astype(float)
    if np.sum(h) > 0.0:
        h /= np.sum(h)
    #endif
    return h, 0.5*(edges[:-1] + edges[1:])


def make_pi0_control_canvases(
    period: PeriodConfig,
    eppi0_data: EPPi0Sample,
    eppi0_mc: EPPi0Sample,
    rows: List[Dict[str, object]],
    outdir: Path,
    mgg_min: float,
    mgg_max: float,
    emiss_abs_max: float,
    ft_theta_max: float,
    pi0_energy_edges: np.ndarray,
    ptmiss_max: float,
) -> None:
    da = eppi0_control_base_arrays(eppi0_data)
    ma = eppi0_control_base_arrays(eppi0_mc)
    ds = eppi0_control_selection_masks(da, mgg_min, mgg_max, emiss_abs_max)
    ms = eppi0_control_selection_masks(ma, mgg_min, mgg_max, emiss_abs_max)
    fd = ds["mass_plus_emiss"]; fm = ms["mass_plus_emiss"]

    fig, ax = plt.subplots(2, 3, figsize=(16.0, 9.5))

    # Mgg before control cuts.
    for arr, mask, label in ((da, ds["baseline"], "Data"), (ma, ms["baseline"], "aaogen")):
        h, c = normalized_histogram(arr["mgg"][mask], np.linspace(0.108, 0.162, 271))
        ax[0,0].step(c, h, where="mid", label=label)
    #endfor
    ax[0,0].axvline(mgg_min, linestyle="--"); ax[0,0].axvline(mgg_max, linestyle="--")
    ax[0,0].set_xlim(0.108, 0.162); ax[0,0].set_xlabel(r"$M_{\gamma\gamma}$ (GeV)"); ax[0,0].set_ylabel("Unit-normalized entries"); ax[0,0].set_title(r"Retained $M_{\gamma\gamma}$ peak (upstream skim: 0.11–0.16 GeV)"); ax[0,0].legend(); ax[0,0].grid(alpha=.18)

    if "pTmiss" in da and "pTmiss" in ma:
        edges = np.linspace(0.0, ptmiss_max, 241)
        for arr, mask, label in ((da, fd, "Data"), (ma, fm, "aaogen")):
            h, c = normalized_histogram(arr["pTmiss"][mask], edges)
            ax[0,1].step(c, h, where="mid", label=label)
        #endfor
        ax[0,1].set_xlabel(r"$p_{T,\mathrm{miss}}$ (GeV)"); ax[0,1].set_ylabel("Unit-normalized entries"); ax[0,1].set_title("pTmiss full control range"); ax[0,1].legend(); ax[0,1].grid(alpha=.18)

        final = next((r for r in rows if r["selection"]=="mass_plus_emiss" and r["group"]=="all"), None)
        if final and np.isfinite(final.get("mc_q50_GeV", np.nan)) and np.isfinite(final.get("mc_sigma_core_GeV", np.nan)):
            center=float(final["mc_q50_GeV"]); sig=max(float(final["mc_sigma_core_GeV"]),1e-4)
            lo=max(0.0, center-4.0*sig); hi=min(ptmiss_max, center+4.0*sig)
        else:
            lo,hi=0.0,min(ptmiss_max,0.4)
        #endif
        edges_core=np.linspace(lo,hi,281)
        hd,c=normalized_histogram(da["pTmiss"][fd],edges_core); hm,_=normalized_histogram(ma["pTmiss"][fm],edges_core)
        ax[0,2].step(c,hd,where="mid",label="Data"); ax[0,2].step(c,hm,where="mid",label="aaogen")
        ax[0,2].set_xlabel(r"$p_{T,\mathrm{miss}}$ (GeV)"); ax[0,2].set_ylabel("Unit-normalized entries"); ax[0,2].set_title("pTmiss robust-core zoom"); ax[0,2].legend(); ax[0,2].grid(alpha=.18)

        ratio=np.divide(hd,hm,out=np.full_like(hd,np.nan),where=hm>0)
        ax[1,2].plot(c,ratio,marker=".",linestyle="none",ms=3); ax[1,2].axhline(1.0,linestyle="--",linewidth=.8)
        ax[1,2].set_xlabel(r"$p_{T,\mathrm{miss}}$ (GeV)"); ax[1,2].set_ylabel("Data / aaogen (normalized)"); ax[1,2].set_title("pTmiss core shape ratio"); ax[1,2].grid(alpha=.18)
    else:
        ax[0,1].text(.5,.5,"pTmiss unavailable",ha="center",va="center"); ax[0,1].set_axis_off()
        ax[0,2].set_axis_off(); ax[1,2].set_axis_off()
    #endif

    if "Emiss2" in da and "Emiss2" in ma:
        edges=np.linspace(-1.5,1.5,241)
        for arr, mask, label in ((da, ds["mass_window"], "Data"), (ma, ms["mass_window"], "aaogen")):
            h,c=normalized_histogram(arr["Emiss2"][mask],edges); ax[1,0].step(c,h,where="mid",label=label)
        #endfor
        ax[1,0].axvline(-emiss_abs_max,linestyle="--"); ax[1,0].axvline(emiss_abs_max,linestyle="--")
        ax[1,0].set_xlabel(r"$E_{\mathrm{miss}}$ (GeV)"); ax[1,0].set_ylabel("Unit-normalized entries"); ax[1,0].set_title("Loose exclusivity control"); ax[1,0].legend(); ax[1,0].grid(alpha=.18)
    else:
        ax[1,0].text(.5,.5,"Emiss2 unavailable",ha="center",va="center"); ax[1,0].set_axis_off()
    #endif

    edges=np.linspace(max(0.0,pi0_energy_edges[0]-0.5),pi0_energy_edges[-1],241)
    for arr, mask, label in ((da,fd,"Data"),(ma,fm,"aaogen")):
        h,c=normalized_histogram(arr["pi0_energy"][mask],edges); ax[1,1].step(c,h,where="mid",label=label)
    #endfor
    for edge in pi0_energy_edges: ax[1,1].axvline(edge,linewidth=.5,alpha=.35)
    #endfor
    ax[1,1].set_xlabel(r"$E_{\pi^0}^{\mathrm{reco}}$ (GeV)"); ax[1,1].set_ylabel("Unit-normalized entries"); ax[1,1].set_title("Control energy coverage"); ax[1,1].legend(); ax[1,1].grid(alpha=.18)

    fig.suptitle(f"{period.label}: reconstructed eppi0 data/aaogen control validation",fontsize=14)
    safe_finalize_figure(fig,outdir/"canvas_pi0_control_validation.png",rect=(0,0,1,.95)); plt.close(fig)

    # Dependence canvas from robust summary rows.
    fig, ax = plt.subplots(1,2,figsize=(14.0,5.4))
    energy_rows=[]
    for r in rows:
        m=re.fullmatch(r"pi0_energy_([0-9.]+)_([0-9.]+)_GeV",str(r["group"]))
        if m and r["selection"]=="mass_plus_emiss":
            energy_rows.append((0.5*(float(m.group(1))+float(m.group(2))),r))
        #endif
    #endfor
    energy_rows.sort(key=lambda z:z[0])
    if energy_rows:
        x=np.asarray([z[0] for z in energy_rows])
        ax[0].plot(x,[z[1]["data_sigma_core_GeV"] for z in energy_rows],marker="o",label="Data")
        ax[0].plot(x,[z[1]["mc_sigma_core_GeV"] for z in energy_rows],marker="o",label="aaogen")
        ax[1].plot(x,[z[1]["sigma_add_core_GeV"] for z in energy_rows],marker="o",label=r"$\sigma_{\mathrm{add}}$")
    #endif
    ax[0].set_xlabel(r"$E_{\pi^0}^{\mathrm{reco}}$ (GeV)"); ax[0].set_ylabel(r"robust $p_{T,\mathrm{miss}}$ core width (GeV)"); ax[0].set_title("Core width vs reconstructed pi0 energy"); ax[0].legend(); ax[0].grid(alpha=.18)
    ax[1].set_xlabel(r"$E_{\pi^0}^{\mathrm{reco}}$ (GeV)"); ax[1].set_ylabel(r"inferred $\sigma_{\mathrm{add}}$ (GeV)"); ax[1].set_title("Extra core width suggested by data"); ax[1].grid(alpha=.18)
    fig.suptitle(f"{period.label}: pTmiss control-core dependence",fontsize=14)
    safe_finalize_figure(fig,outdir/"canvas_pi0_control_dependence.png",rect=(0,0,1,.94)); plt.close(fig)


    # Detector-topology control canvas using the confirmed photon detector codes.
    topology_names = ("FD_FD", "FT_FD", "FT_FT", "FT_containing")
    topology_rows = {
        str(r["group"]): r
        for r in rows
        if r["selection"] == "mass_plus_emiss"
        and str(r["group"]) in topology_names
    }

    if topology_rows and "pTmiss" in da and "pTmiss" in ma:
        if all(name in da for name in ("detector_gamma1", "detector_gamma2")) and all(
            name in ma for name in ("detector_gamma1", "detector_gamma2")
        ):
            dp1 = np.minimum(da["detector_gamma1"], da["detector_gamma2"])
            dp2 = np.maximum(da["detector_gamma1"], da["detector_gamma2"])
            mp1 = np.minimum(ma["detector_gamma1"], ma["detector_gamma2"])
            mp2 = np.maximum(ma["detector_gamma1"], ma["detector_gamma2"])

            topology_masks = {
                "FD_FD": (
                    fd & (dp1 == PHOTON_DETECTOR_FD) & (dp2 == PHOTON_DETECTOR_FD),
                    fm & (mp1 == PHOTON_DETECTOR_FD) & (mp2 == PHOTON_DETECTOR_FD),
                ),
                "FT_FD": (
                    fd & (dp1 == PHOTON_DETECTOR_FT) & (dp2 == PHOTON_DETECTOR_FD),
                    fm & (mp1 == PHOTON_DETECTOR_FT) & (mp2 == PHOTON_DETECTOR_FD),
                ),
                "FT_FT": (
                    fd & (dp1 == PHOTON_DETECTOR_FT) & (dp2 == PHOTON_DETECTOR_FT),
                    fm & (mp1 == PHOTON_DETECTOR_FT) & (mp2 == PHOTON_DETECTOR_FT),
                ),
            }

            fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.4))
            for ax_top, topology in zip(axes, ("FD_FD", "FT_FD", "FT_FT")):
                dm, mm = topology_masks[topology]
                row = topology_rows.get(topology, {})
                center = float(row.get("mc_q50_GeV", 0.15))
                sig = max(float(row.get("mc_sigma_core_GeV", 0.05)), 1.0e-3)
                lo = max(0.0, center - 5.0 * sig)
                hi = min(ptmiss_max, center + 7.0 * sig)
                if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
                    lo, hi = 0.0, min(ptmiss_max, 0.6)
                #endif
                edges_top = np.linspace(lo, hi, 241)
                hd, c = normalized_histogram(da["pTmiss"][dm], edges_top)
                hm, _ = normalized_histogram(ma["pTmiss"][mm], edges_top)
                ax_top.step(c, hd, where="mid", label="Data")
                ax_top.step(c, hm, where="mid", label="aaogen")
                ax_top.set_xlabel(r"$p_{T,\mathrm{miss}}$ (GeV)")
                ax_top.set_ylabel("Unit-normalized entries")
                ax_top.set_title(
                    f"{topology.replace('_', '–')}\n"
                    f"N(data)={int(np.count_nonzero(dm)):,}, "
                    f"N(MC)={int(np.count_nonzero(mm)):,}"
                )
                ax_top.legend(fontsize=8)
                ax_top.grid(alpha=0.18)
            #endfor

            fig.suptitle(
                rf"{period.label}: reconstructed $ep\pi^0$ control by photon detector topology "
                "(0 = FT, 1 = FD)",
                fontsize=14,
            )
            safe_finalize_figure(
                fig,
                outdir / "canvas_pi0_control_detector_topology.png",
                rect=(0, 0, 1, 0.93),
            )
            plt.close(fig)
        #endif
    #endif


def fit_shared_morphed_composition(
    histograms: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    pi0_control: Dict[str, object],
    nuisance_shift_prior: float,
    nuisance_sigma_prior: float,
    max_shift_bins: float,
    max_sigma_bins: float,
) -> SharedMorphedFitResult:
    """
    Shared two-driver composition fit with profile-initialized template morphing.

    The earlier simultaneous L-BFGS-B fit frequently remained exactly at the
    zero-smearing boundary even when explicit profile scans strongly preferred
    non-zero smearing.  v021 therefore uses:

      1. zero-centered weak nuisance priors;
      2. one-coordinate-at-a-time grid profiling with f_pi0 re-profiled at
         every candidate point;
      3. repeated sweeps until stable;
      4. optional continuous L-BFGS-B refinement starting from the non-zero
         profiled solution.

    The reconstructed-eppi0 control remains diagnostic-only.
    """
    names = [n for n in STAGE2_DRIVER_DISCRIMINATORS if n in histograms]
    if not names:
        return SharedMorphedFitResult(False, "no available driver discriminators")
    #endif

    data_hists = {
        name: np.asarray(histograms[name][0], dtype=float)
        for name in names
    }

    # Nuisance state is stored explicitly by (driver, component, nuisance).
    nuisance_state: Dict[str, float] = {}
    for name in names:
        for component in ("pi0", "dvcs"):
            nuisance_state[f"{name}_{component}_shift_bins"] = 0.0
            nuisance_state[f"{name}_{component}_sigma_bins"] = 0.0
        #endfor
    #endfor

    shift_grid = np.arange(
        -float(max_shift_bins),
        float(max_shift_bins) + 0.5,
        1.0,
        dtype=float,
    )
    sigma_grid = np.arange(
        0.0,
        float(max_sigma_bins) + 0.25,
        0.5,
        dtype=float,
    )

    def morphed_templates(
        state: Dict[str, float],
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        tp: Dict[str, np.ndarray] = {}
        td: Dict[str, np.ndarray] = {}
        for name in names:
            _, hp, hv = histograms[name]
            tp[name] = morph_template_second_axis(
                hp,
                state[f"{name}_pi0_shift_bins"],
                state[f"{name}_pi0_sigma_bins"],
            )
            td[name] = morph_template_second_axis(
                hv,
                state[f"{name}_dvcs_shift_bins"],
                state[f"{name}_dvcs_sigma_bins"],
            )
        #endfor
        return tp, td
    #enddef

    def nuisance_penalty(state: Dict[str, float]) -> float:
        penalty = 0.0
        for key, value in state.items():
            width = (
                nuisance_sigma_prior
                if key.endswith("_sigma_bins")
                else nuisance_shift_prior
            )
            penalty += 0.5 * (float(value) / float(width)) ** 2
        #endfor
        return float(penalty)
    #enddef

    def profiled_objective(
        state: Dict[str, float],
    ) -> Tuple[float, float]:
        tp, td = morphed_templates(state)
        f_prof, nll = _profile_fraction_for_templates(
            data_hists,
            tp,
            td,
        )
        return f_prof, float(nll + nuisance_penalty(state))
    #enddef

    f_best, objective_best = profiled_objective(nuisance_state)

    # Alternating coordinate-profile sweeps.
    max_sweeps = 4
    for sweep in range(max_sweeps):
        changed = False

        for name in names:
            for component in ("pi0", "dvcs"):
                for nuisance in ("shift", "sigma"):
                    key = f"{name}_{component}_{nuisance}_bins"
                    grid = shift_grid if nuisance == "shift" else sigma_grid

                    current_value = float(nuisance_state[key])
                    candidate_values = np.unique(
                        np.concatenate((grid, np.asarray([current_value])))
                    )

                    best_local_value = current_value
                    best_local_f = f_best
                    best_local_objective = objective_best

                    for candidate in candidate_values:
                        trial = dict(nuisance_state)
                        trial[key] = float(candidate)
                        f_trial, objective_trial = profiled_objective(trial)

                        if objective_trial < best_local_objective - 1.0e-8:
                            best_local_value = float(candidate)
                            best_local_f = float(f_trial)
                            best_local_objective = float(objective_trial)
                        #endif
                    #endfor

                    if abs(best_local_value - current_value) > 1.0e-10:
                        nuisance_state[key] = best_local_value
                        f_best = best_local_f
                        objective_best = best_local_objective
                        changed = True
                    #endif
                #endfor
            #endfor
        #endfor

        if not changed:
            break
        #endif
    #endfor

    # Record the profile-grid initialization before continuous refinement.
    profile_initial = dict(nuisance_state)
    profile_initial_f = float(f_best)
    profile_initial_objective = float(objective_best)

    # Continuous local refinement starts from the profile solution rather than
    # from the problematic all-zero boundary.
    x0 = [math.log(profile_initial_f / (1.0 - profile_initial_f))]
    bounds = [(-8.0, 8.0)]
    ordered_keys: List[str] = []

    for name in names:
        for component in ("pi0", "dvcs"):
            for nuisance in ("shift", "sigma"):
                key = f"{name}_{component}_{nuisance}_bins"
                ordered_keys.append(key)
                x0.append(profile_initial[key])
                if nuisance == "shift":
                    bounds.append((-max_shift_bins, max_shift_bins))
                else:
                    bounds.append((0.0, max_sigma_bins))
                #endif
            #endfor
        #endfor
    #endfor

    def continuous_objective(x: np.ndarray) -> float:
        f = 1.0 / (1.0 + math.exp(-float(x[0])))
        state = {
            key: float(value)
            for key, value in zip(ordered_keys, x[1:])
        }
        tp, td = morphed_templates(state)

        total_nll = 0.0
        for name in names:
            data = data_hists[name].ravel()
            nd = float(np.sum(data))
            mu = np.clip(
                nd * (
                    f * tp[name].ravel()
                    + (1.0 - f) * td[name].ravel()
                ),
                1.0e-12,
                None,
            )
            total_nll += float(
                np.sum(mu - data * np.log(mu))
            ) / len(names)
        #endfor

        return float(total_nll + nuisance_penalty(state))
    #enddef

    res = minimize(
        continuous_objective,
        np.asarray(x0, dtype=float),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 500, "ftol": 1.0e-10, "maxls": 40},
    )

    # Use the continuous result only if it actually improves the profiled
    # coordinate solution. This prevents a local optimizer regression.
    use_continuous = (
        bool(res.success)
        and np.isfinite(res.fun)
        and float(res.fun) <= profile_initial_objective + 1.0e-8
    )

    if use_continuous:
        f_final = 1.0 / (1.0 + math.exp(-float(res.x[0])))
        final_state = {
            key: float(value)
            for key, value in zip(ordered_keys, res.x[1:])
        }
        objective_final = float(res.fun)
        message = (
            "profile-grid initialization + continuous refinement: "
            + str(res.message)
        )
    else:
        f_final = profile_initial_f
        final_state = profile_initial
        objective_final = profile_initial_objective
        message = (
            "profile-grid/alternating solution retained; continuous refinement "
            "did not improve the objective"
        )
    #endif

    nuisance: Dict[str, float] = {
        "control_prior_applied": 0.0,
        "profile_initialized": 1.0,
        "profile_initial_pi0_fraction": profile_initial_f,
        "profile_initial_objective": profile_initial_objective,
        "continuous_refinement_used": float(use_continuous),
    }
    nuisance.update(final_state)

    # Preserve the grid-initial values for direct debugging.
    for key, value in profile_initial.items():
        nuisance[f"{key}_profile_initial"] = float(value)
    #endfor

    total_dev = 0.0
    total_ndof = 0

    for name in names:
        hd, hp, hv = histograms[name]
        tp = morph_template_second_axis(
            hp,
            final_state[f"{name}_pi0_shift_bins"],
            final_state[f"{name}_pi0_sigma_bins"],
        ).ravel()
        td = morph_template_second_axis(
            hv,
            final_state[f"{name}_dvcs_shift_bins"],
            final_state[f"{name}_dvcs_sigma_bins"],
        ).ravel()
        data = np.asarray(hd, dtype=float).ravel()
        nd = float(np.sum(data))
        mu = nd * (f_final * tp + (1.0 - f_final) * td)

        dev, ndof = poisson_deviance_quality(
            data,
            mu,
            1 + 4 * len(names),
        )
        nuisance[f"{name}_deviance_per_ndof"] = float(dev / ndof)
        total_dev += dev
        total_ndof += ndof
    #endfor

    return SharedMorphedFitResult(
        True,
        message,
        float(f_final),
        objective_final,
        float(total_dev),
        int(total_ndof),
        nuisance,
    )




def build_stage2_template_mixture_closure(
    period: PeriodConfig,
    region: str,
    energy_low: float,
    energy_high: float,
    histograms: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
    nominal_fit: SharedMorphedFitResult,
    pi0_control: Dict[str, object],
    truth_fractions: Sequence[float],
    nuisance_shift_prior: float,
    nuisance_sigma_prior: float,
    max_shift_bins: float,
    max_sigma_bins: float,
) -> List[Dict[str, object]]:
    """
    Deterministic Stage-II template-mixture closure using the *real-data best-fit
    morphed shapes* as the pseudo-data generator.

    This is deliberately harder than mixing the raw templates with themselves.
    For every production driver we first apply the nuisance shifts/smearings
    found in the nominal real-data fit, inject a controlled pi0 fraction, and
    then rerun the complete shared-morphed fit from its ordinary zero-nuisance
    initialization.  Any failure to recover the injected fraction therefore
    exposes optimizer/boundary/identifiability problems inside the model family.
    """
    if (
        region not in ("FT", "FD_all")
        or not nominal_fit.success
        or not nominal_fit.nuisance
        or not histograms
    ):
        return []
    #endif

    names = [n for n in STAGE2_DRIVER_DISCRIMINATORS if n in histograms]
    if not names:
        return []
    #endif

    injected_templates: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for name in names:
        _, hp, hv = histograms[name]
        tp = morph_template_second_axis(
            hp,
            float(nominal_fit.nuisance.get(f"{name}_pi0_shift_bins", 0.0)),
            float(nominal_fit.nuisance.get(f"{name}_pi0_sigma_bins", 0.0)),
        )
        td = morph_template_second_axis(
            hv,
            float(nominal_fit.nuisance.get(f"{name}_dvcs_shift_bins", 0.0)),
            float(nominal_fit.nuisance.get(f"{name}_dvcs_sigma_bins", 0.0)),
        )
        injected_templates[name] = (tp, td)
    #endfor

    rows: List[Dict[str, object]] = []
    for truth in truth_fractions:
        truth = float(truth)
        if not (0.0 < truth < 1.0):
            continue
        #endif

        pseudo: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for name in names:
            hd0, hp0, hv0 = histograms[name]
            nd = max(1.0, float(np.sum(hd0)))
            tp, td = injected_templates[name]
            pseudo_data = nd * (truth * tp + (1.0 - truth) * td)
            pseudo[name] = (pseudo_data, hp0, hv0)
        #endfor

        refit = fit_shared_morphed_composition(
            pseudo,
            pi0_control,
            nuisance_shift_prior,
            nuisance_sigma_prior,
            max_shift_bins,
            max_sigma_bins,
        )

        row: Dict[str, object] = {
            "period": period.key,
            "label": period.label,
            "beam_energy_GeV": period.beam_energy,
            "region": region,
            "energy_low_GeV": float(energy_low),
            "energy_high_GeV": float(energy_high),
            "energy_center_GeV": 0.5 * (float(energy_low) + float(energy_high)),
            "closure_model": "real-data-best-fit-morphed-Asimov",
            "n_drivers": len(names),
            "injected_pi0_fraction": truth,
            "fit_success": int(refit.success),
            "recovered_pi0_fraction": float(refit.pi0_fraction),
            "closure_bias": (
                float(refit.pi0_fraction - truth)
                if np.isfinite(refit.pi0_fraction)
                else float("nan")
            ),
            "abs_closure_bias": (
                abs(float(refit.pi0_fraction - truth))
                if np.isfinite(refit.pi0_fraction)
                else float("nan")
            ),
            "refit_deviance_per_ndof": (
                float(refit.poisson_deviance / refit.ndof)
                if refit.ndof else float("nan")
            ),
            "nominal_real_data_pi0_fraction": float(nominal_fit.pi0_fraction),
            "nominal_real_data_deviance_per_ndof": (
                float(nominal_fit.poisson_deviance / nominal_fit.ndof)
                if nominal_fit.ndof else float("nan")
            ),
        }

        # Record the actual injected morph state and the refitted state.  This
        # lets us distinguish a composition bias from a nuisance degeneracy.
        for name in names:
            for component in ("pi0", "dvcs"):
                for nuisance in ("shift", "sigma"):
                    key = f"{name}_{component}_{nuisance}_bins"
                    row[f"injected_{key}"] = float(
                        nominal_fit.nuisance.get(key, 0.0)
                    )
                    row[f"recovered_{key}"] = float(
                        (refit.nuisance or {}).get(key, float("nan"))
                    )
                #endfor
            #endfor
        #endfor

        rows.append(row)
    #endfor

    return rows


def make_stage2_template_mixture_closure_canvas(
    closure_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """Compact cross-period closure summary for FT and FD_all."""
    if not closure_rows:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(15.8, 10.0))

    for icol, region in enumerate(("FT", "FD_all")):
        rr_region = [
            r for r in closure_rows
            if str(r.get("region", "")) == region
            and int(r.get("fit_success", 0)) == 1
        ]

        # Top: recovered versus injected fraction. Each period is a separate
        # line/marker family; all energy bins are intentionally overlaid.
        ax = axes[0, icol]
        for period in PERIODS:
            rr = [r for r in rr_region if str(r.get("period", "")) == period.key]
            if not rr:
                continue
            #endif
            ax.scatter(
                [float(r["injected_pi0_fraction"]) for r in rr],
                [float(r["recovered_pi0_fraction"]) for r in rr],
                s=24,
                alpha=0.70,
                label=period.label,
            )
        #endfor
        ax.plot([0.0, 1.0], [0.0, 1.0], color="black", linestyle="--", linewidth=1.0)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel(r"injected $f_{\pi^0}$")
        ax.set_ylabel(r"recovered $f_{\pi^0}$")
        ax.set_title("FT" if region == "FT" else "FD all")
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8, frameon=False)

        # Bottom: worst absolute closure bias across injected fractions in each
        # energy bin. This makes energy-localized identifiability failures easy
        # to see without another forest of curves.
        ax = axes[1, icol]
        for period in PERIODS:
            rr = [r for r in rr_region if str(r.get("period", "")) == period.key]
            if not rr:
                continue
            #endif
            grouped: Dict[float, List[float]] = {}
            for row in rr:
                e = float(row["energy_center_GeV"])
                grouped.setdefault(e, []).append(abs(float(row["closure_bias"])))
            #endfor
            xs = sorted(grouped)
            ys = [max(grouped[x]) for x in xs]
            ax.plot(xs, ys, marker="o", linewidth=1.1, label=period.label)
        #endfor
        ax.axhline(0.01, color="black", linestyle=":", linewidth=1.0, label="1% bias")
        ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        ax.set_ylabel(r"max injected-fraction $|\Delta f_{\pi^0}|$")
        ax.set_ylim(0.0, 0.10)
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8, frameon=False)
    #endfor

    fig.suptitle(
        "Stage-II template-mixture closure with real-data best-fit morphing",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        Path(outdir) / "canvas_template_mixture_closure.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)


def run_stage2_ft_coarse_shared_morphed_fits(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    pi0_control: Dict[str, object],
    ft_theta_max: float,
    energy_edges: Sequence[float],
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
    min_data_count: int,
    min_template_count: int,
    nuisance_shift_prior: float,
    nuisance_sigma_prior: float,
    max_shift_bins: float,
    max_sigma_bins: float,
    closure_truth_fractions: Sequence[float],
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """
    Diagnostic coarse-energy FT version of the nominal shared-morphed Stage-II
    fit.

    Physics model, discriminators, nuisance priors, and template morphing are
    identical to the nominal fit. Only the FT predicted-probe energy binning is
    coarsened. These results are *not* used by Stage III.
    """
    edges = np.asarray(energy_edges, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or np.any(~np.isfinite(edges)):
        raise ValueError("FT coarse energy edges must be a finite 1D array.")
    #endif
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("FT coarse energy edges must be strictly increasing.")
    #endif

    rows: List[Dict[str, object]] = []
    closure_rows: List[Dict[str, object]] = []

    for ib in range(len(edges) - 1):
        elo = float(edges[ib])
        ehi = float(edges[ib + 1])

        masks = {
            key: stage2_fit_mask(
                feat,
                "FT",
                ft_theta_max,
                elo,
                ehi,
                mm2_min,
                mm2_max,
                probe_m2_max,
            )
            for key, feat in (
                ("data", data_f),
                ("pi0", pi0_f),
                ("dvcs", dvcs_f),
            )
        }

        row: Dict[str, object] = {
            "period": period.key,
            "label": period.label,
            "beam_energy_GeV": period.beam_energy,
            "region": "FT",
            "energy_low_GeV": elo,
            "energy_high_GeV": ehi,
            "energy_center_GeV": 0.5 * (elo + ehi),
            "fit_model": "shared_morphed_2d_FT_coarse_energy",
            "data_candidate_count": int(np.count_nonzero(masks["data"])),
            "aaogen_candidate_count": int(np.count_nonzero(masks["pi0"])),
            "dvcsgen_candidate_count": int(np.count_nonzero(masks["dvcs"])),
        }

        hist: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for disc in STAGE2_DRIVER_DISCRIMINATORS:
            if not discriminator_available(disc, data_f, pi0_f, dvcs_f):
                continue
            #endif
            hd = histogram_for_discriminator(
                disc, data_f, masks["data"],
                mm2_min=mm2_min, mm2_max=mm2_max,
                probe_m2_max=probe_m2_max,
                mm2_bins_2d=mm2_bins_2d,
                probe_m2_bins_2d=probe_m2_bins_2d,
                bins_1d=90,
                ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins,
                theta_max=theta_max, theta_bins=theta_bins,
            )
            hp = histogram_for_discriminator(
                disc, pi0_f, masks["pi0"],
                mm2_min=mm2_min, mm2_max=mm2_max,
                probe_m2_max=probe_m2_max,
                mm2_bins_2d=mm2_bins_2d,
                probe_m2_bins_2d=probe_m2_bins_2d,
                bins_1d=90,
                ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins,
                theta_max=theta_max, theta_bins=theta_bins,
            )
            hv = histogram_for_discriminator(
                disc, dvcs_f, masks["dvcs"],
                mm2_min=mm2_min, mm2_max=mm2_max,
                probe_m2_max=probe_m2_max,
                mm2_bins_2d=mm2_bins_2d,
                probe_m2_bins_2d=probe_m2_bins_2d,
                bins_1d=90,
                ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins,
                theta_max=theta_max, theta_bins=theta_bins,
            )

            if (
                float(np.sum(hd)) >= min_data_count
                and float(np.sum(hp)) >= min_template_count
                and float(np.sum(hv)) >= min_template_count
            ):
                hist[disc] = (hd, hp, hv)
            #endif
        #endfor

        fit = fit_shared_morphed_composition(
            hist,
            pi0_control,
            nuisance_shift_prior,
            nuisance_sigma_prior,
            max_shift_bins,
            max_sigma_bins,
        )

        row.update({
            "fit_success": int(fit.success),
            "fit_message": fit.message,
            "n_drivers": len(hist),
            "pi0_fraction": float(fit.pi0_fraction),
            "dvcs_fraction": (
                1.0 - float(fit.pi0_fraction)
                if np.isfinite(fit.pi0_fraction)
                else float("nan")
            ),
            "objective": float(fit.objective),
            "poisson_deviance": float(fit.poisson_deviance),
            "ndof": int(fit.ndof),
            "deviance_per_ndof": (
                float(fit.poisson_deviance / fit.ndof)
                if fit.ndof else float("nan")
            ),
        })
        if fit.nuisance:
            row.update(fit.nuisance)
        #endif
        rows.append(row)

        if hist and fit.success:
            closure_rows.extend(
                build_stage2_template_mixture_closure(
                    period,
                    "FT",
                    elo,
                    ehi,
                    hist,
                    fit,
                    pi0_control,
                    closure_truth_fractions,
                    nuisance_shift_prior,
                    nuisance_sigma_prior,
                    max_shift_bins,
                    max_sigma_bins,
                )
            )
        #endif
    #endfor

    return rows, closure_rows


def make_ft_coarse_composition_canvas(
    period: PeriodConfig,
    fine_rows: List[Dict[str, object]],
    coarse_rows: List[Dict[str, object]],
    coarse_closure_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Compare the nominal fine-bin FT composition with the diagnostic coarse-bin
    FT composition and show closure quality of the coarse fits.
    """
    if not coarse_rows:
        return
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(15.2, 5.7))

    fine = sorted(
        [
            r for r in fine_rows
            if str(r.get("region", "")) == "FT"
            and int(r.get("fit_success", 0)) == 1
        ],
        key=lambda r: float(r["energy_center_GeV"]),
    )
    coarse = sorted(
        [r for r in coarse_rows if int(r.get("fit_success", 0)) == 1],
        key=lambda r: float(r["energy_low_GeV"]),
    )

    if fine:
        axes[0].plot(
            [float(r["energy_center_GeV"]) for r in fine],
            [float(r["pi0_fraction"]) for r in fine],
            marker="o",
            linewidth=1.1,
            label="nominal fine bins",
        )
    #endif

    for i, row in enumerate(coarse):
        elo = float(row["energy_low_GeV"])
        ehi = float(row["energy_high_GeV"])
        f = float(row["pi0_fraction"])
        axes[0].hlines(
            f,
            elo,
            ehi,
            linewidth=3.0,
            label="coarse FT fit" if i == 0 else None,
        )
        axes[0].plot(
            [0.5 * (elo + ehi)],
            [f],
            marker="s",
            markersize=6,
        )
    #endfor

    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0].set_ylabel(r"fitted $f_{\pi^0}$")
    axes[0].set_title("FT composition: fine versus coarse energy bins")
    axes[0].grid(alpha=0.20)
    axes[0].legend(frameon=False, fontsize=8)

    # Worst closure bias per coarse energy bin.
    grouped: Dict[Tuple[float, float], List[float]] = {}
    for row in coarse_closure_rows:
        if int(row.get("fit_success", 0)) != 1:
            continue
        #endif
        key = (
            float(row["energy_low_GeV"]),
            float(row["energy_high_GeV"]),
        )
        grouped.setdefault(key, []).append(
            abs(float(row["closure_bias"]))
        )
    #endfor

    for (elo, ehi), values in sorted(grouped.items()):
        axes[1].bar(
            0.5 * (elo + ehi),
            max(values),
            width=0.85 * (ehi - elo),
            alpha=0.65,
        )
    #endfor
    axes[1].axhline(
        0.01,
        color="black",
        linestyle=":",
        linewidth=1.0,
        label="1% bias",
    )
    axes[1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1].set_ylabel(r"max injected-fraction $|\Delta f_{\pi^0}|$")
    axes[1].set_title("Coarse-bin FT template-mixture closure")
    axes[1].grid(alpha=0.20)
    axes[1].legend(frameon=False, fontsize=8)

    fig.suptitle(
        f"{period.label}: diagnostic coarse-energy FT composition",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        Path(outdir) / "canvas_ft_coarse_energy_composition.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)



@dataclass
class SharedThreeComponentDiagnosticResult:
    success: bool
    message: str = ""
    pi0_fraction: float = float("nan")
    dvcs_fraction: float = float("nan")
    mixed_fraction: float = float("nan")
    objective: float = float("nan")
    poisson_deviance: float = float("nan")
    ndof: int = 0
    driver_pi0_fraction: Optional[Dict[str, float]] = None
    driver_mixed_fraction: Optional[Dict[str, float]] = None


def _softmax_three_from_two(x: Sequence[float]) -> Tuple[float, float, float]:
    """Stable 3-component softmax with BH/DVCS as the zero-logit reference."""
    a = float(np.clip(x[0], -40.0, 40.0))
    b = float(np.clip(x[1], -40.0, 40.0))
    ea = math.exp(a)
    eb = math.exp(b)
    den = 1.0 + ea + eb
    return ea / den, 1.0 / den, eb / den


def _logits_from_three_fractions(
    f_pi0: float,
    f_dvcs: float,
    f_mixed: float,
) -> np.ndarray:
    eps = 1.0e-8
    p = max(eps, float(f_pi0))
    d = max(eps, float(f_dvcs))
    m = max(eps, float(f_mixed))
    norm = p + d + m
    p /= norm
    d /= norm
    m /= norm
    return np.asarray([math.log(p / d), math.log(m / d)], dtype=float)


def _fixed_morph_three_component_templates(
    histograms: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    two_component_fit: SharedMorphedFitResult,
) -> Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Return data plus fixed-morphed pi0/DVCS and raw mixed templates."""
    out: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    if not two_component_fit.nuisance:
        return out
    #endif
    for name in STAGE2_DRIVER_DISCRIMINATORS:
        if name not in histograms:
            continue
        #endif
        hd, hp, hv, hm = histograms[name]
        if min(float(np.sum(hp)), float(np.sum(hv)), float(np.sum(hm))) <= 0.0:
            continue
        #endif
        tp = morph_template_second_axis(
            hp,
            float(two_component_fit.nuisance.get(f"{name}_pi0_shift_bins", 0.0)),
            float(two_component_fit.nuisance.get(f"{name}_pi0_sigma_bins", 0.0)),
        )
        td = morph_template_second_axis(
            hv,
            float(two_component_fit.nuisance.get(f"{name}_dvcs_shift_bins", 0.0)),
            float(two_component_fit.nuisance.get(f"{name}_dvcs_sigma_bins", 0.0)),
        )
        tm = normalized_template(hm).reshape(hm.shape)
        out[name] = (np.asarray(hd, dtype=float), tp, td, tm)
    #endfor
    return out


def fit_shared_three_component_fixed_morph(
    histograms: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    two_component_fit: SharedMorphedFitResult,
) -> SharedThreeComponentDiagnosticResult:
    """
    Shared coarse-FT diagnostic fit with BH/DVCS + pi0 + mixed-event photons.

    The pi0 and BH/DVCS morph nuisance values are fixed to the corresponding
    two-component coarse-FT fit. Only the common composition fractions float.
    This tests a missing uncorrelated/wrong-photon population without giving the
    templates additional shape freedom that could worsen the known degeneracy.
    """
    if not two_component_fit.success or not two_component_fit.nuisance:
        return SharedThreeComponentDiagnosticResult(False, "two-component coarse fit unavailable")
    #endif

    shaped = _fixed_morph_three_component_templates(histograms, two_component_fit)
    names = [n for n in STAGE2_DRIVER_DISCRIMINATORS if n in shaped]
    if not names:
        return SharedThreeComponentDiagnosticResult(False, "no complete production drivers")
    #endif

    seed_mixed = 0.02
    seed_pi0 = (1.0 - seed_mixed) * float(two_component_fit.pi0_fraction)
    seed_dvcs = max(1.0e-6, 1.0 - seed_mixed - seed_pi0)
    x0 = _logits_from_three_fractions(seed_pi0, seed_dvcs, seed_mixed)

    def objective(x: np.ndarray) -> float:
        fp, fd, fm = _softmax_three_from_two(x)
        total = 0.0
        for name in names:
            hd, tp, td, tm = shaped[name]
            data = hd.ravel()
            nd = float(np.sum(data))
            mu = np.clip(nd * (fp * tp.ravel() + fd * td.ravel() + fm * tm.ravel()), 1.0e-12, None)
            total += float(np.sum(mu - data * np.log(mu))) / len(names)
        #endfor
        return float(total)
    #enddef

    res = minimize(
        objective,
        x0,
        method="L-BFGS-B",
        bounds=((-12.0, 12.0), (-12.0, 12.0)),
        options={"maxiter": 400, "ftol": 1.0e-11, "maxls": 40},
    )
    if not res.success or not np.isfinite(res.fun):
        return SharedThreeComponentDiagnosticResult(False, str(res.message))
    #endif

    fp, fd, fm = _softmax_three_from_two(res.x)
    total_dev = 0.0
    total_ndof = 0
    driver_pi0: Dict[str, float] = {}
    driver_mixed: Dict[str, float] = {}

    for name in names:
        hd, tp, td, tm = shaped[name]
        data = hd.ravel()
        nd = float(np.sum(data))
        mu = nd * (fp * tp.ravel() + fd * td.ravel() + fm * tm.ravel())
        dev, ndof = poisson_deviance_quality(data, mu, n_parameters=2)
        total_dev += dev
        total_ndof += ndof

        def local_objective(x: np.ndarray) -> float:
            lp, ld, lm = _softmax_three_from_two(x)
            lmu = np.clip(nd * (lp * tp.ravel() + ld * td.ravel() + lm * tm.ravel()), 1.0e-12, None)
            return float(np.sum(lmu - data * np.log(lmu)))
        #enddef

        lres = minimize(
            local_objective,
            res.x,
            method="L-BFGS-B",
            bounds=((-12.0, 12.0), (-12.0, 12.0)),
            options={"maxiter": 250, "ftol": 1.0e-11},
        )
        if lres.success and np.all(np.isfinite(lres.x)):
            lp, _, lm = _softmax_three_from_two(lres.x)
            driver_pi0[name] = float(lp)
            driver_mixed[name] = float(lm)
        #endif
    #endfor

    return SharedThreeComponentDiagnosticResult(
        True,
        str(res.message),
        float(fp),
        float(fd),
        float(fm),
        float(res.fun),
        float(total_dev),
        int(total_ndof),
        driver_pi0_fraction=driver_pi0,
        driver_mixed_fraction=driver_mixed,
    )


def build_ft_coarse_three_component_closure(
    period: PeriodConfig,
    energy_low: float,
    energy_high: float,
    histograms: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]],
    two_component_fit: SharedMorphedFitResult,
    truth_points: Sequence[Tuple[float, float]],
) -> List[Dict[str, object]]:
    """Asimov identifiability test for the fixed-morph 3-component FT fit."""
    shaped = _fixed_morph_three_component_templates(histograms, two_component_fit)
    names = [n for n in STAGE2_DRIVER_DISCRIMINATORS if n in shaped]
    if not names:
        return []
    #endif

    rows: List[Dict[str, object]] = []
    for truth_pi0, truth_mixed in truth_points:
        truth_pi0 = float(truth_pi0)
        truth_mixed = float(truth_mixed)
        truth_dvcs = 1.0 - truth_pi0 - truth_mixed
        if min(truth_pi0, truth_mixed, truth_dvcs) < 0.0:
            continue
        #endif

        pseudo: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        for name in names:
            hd, tp, td, tm = shaped[name]
            _, hp0, hv0, hm0 = histograms[name]
            nd = max(1.0, float(np.sum(hd)))
            pseudo_data = nd * (truth_pi0 * tp + truth_dvcs * td + truth_mixed * tm)
            pseudo[name] = (pseudo_data, hp0, hv0, hm0)
        #endfor

        refit = fit_shared_three_component_fixed_morph(pseudo, two_component_fit)
        rows.append({
            "period": period.key,
            "label": period.label,
            "region": "FT",
            "energy_low_GeV": float(energy_low),
            "energy_high_GeV": float(energy_high),
            "energy_center_GeV": 0.5 * (float(energy_low) + float(energy_high)),
            "closure_model": "coarse_FT_three_component_fixed_morph_Asimov",
            "injected_pi0_fraction": truth_pi0,
            "injected_dvcs_fraction": truth_dvcs,
            "injected_mixed_fraction": truth_mixed,
            "fit_success": int(refit.success),
            "recovered_pi0_fraction": float(refit.pi0_fraction),
            "recovered_dvcs_fraction": float(refit.dvcs_fraction),
            "recovered_mixed_fraction": float(refit.mixed_fraction),
            "pi0_closure_bias": float(refit.pi0_fraction - truth_pi0) if np.isfinite(refit.pi0_fraction) else float("nan"),
            "mixed_closure_bias": float(refit.mixed_fraction - truth_mixed) if np.isfinite(refit.mixed_fraction) else float("nan"),
            "max_abs_component_bias": max(abs(float(refit.pi0_fraction - truth_pi0)), abs(float(refit.mixed_fraction - truth_mixed))) if np.isfinite(refit.pi0_fraction) and np.isfinite(refit.mixed_fraction) else float("nan"),
        })
    #endfor
    return rows


def run_ft_coarse_three_component_diagnostic(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    mixed_f: Dict[str, np.ndarray],
    coarse_rows: List[Dict[str, object]],
    shift: int,
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
    min_data_count: int,
    min_template_count: int,
) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    """Run the 3-component diagnostic only in the coarse FT energy bins."""
    rows: List[Dict[str, object]] = []
    closure_rows: List[Dict[str, object]] = []
    truth_points = ((0.20, 0.00), (0.50, 0.00), (0.80, 0.00), (0.20, 0.10), (0.50, 0.10), (0.20, 0.25), (0.50, 0.25))

    for coarse in [r for r in coarse_rows if int(r.get("fit_success", 0)) == 1]:
        elo = float(coarse["energy_low_GeV"])
        ehi = float(coarse["energy_high_GeV"])
        masks = {
            key: stage2_fit_mask(feat, "FT", ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)
            for key, feat in (("data", data_f), ("pi0", pi0_f), ("dvcs", dvcs_f), ("mixed", mixed_f))
        }

        hist4: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
        for disc in STAGE2_DRIVER_DISCRIMINATORS:
            if disc not in ("mx2_ep_x_probe_m2", "mx2_ep_x_pTmiss"):
                continue
            #endif
            if not discriminator_available(disc, data_f, pi0_f, dvcs_f):
                continue
            #endif
            # The mixed feature store explicitly recomputes pTmiss for the
            # substituted photon, so both production drivers are valid here.
            hd = histogram_for_discriminator(disc, data_f, masks["data"], mm2_min=mm2_min, mm2_max=mm2_max, probe_m2_max=probe_m2_max, mm2_bins_2d=mm2_bins_2d, probe_m2_bins_2d=probe_m2_bins_2d, bins_1d=90, ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins, theta_max=theta_max, theta_bins=theta_bins)
            hp = histogram_for_discriminator(disc, pi0_f, masks["pi0"], mm2_min=mm2_min, mm2_max=mm2_max, probe_m2_max=probe_m2_max, mm2_bins_2d=mm2_bins_2d, probe_m2_bins_2d=probe_m2_bins_2d, bins_1d=90, ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins, theta_max=theta_max, theta_bins=theta_bins)
            hv = histogram_for_discriminator(disc, dvcs_f, masks["dvcs"], mm2_min=mm2_min, mm2_max=mm2_max, probe_m2_max=probe_m2_max, mm2_bins_2d=mm2_bins_2d, probe_m2_bins_2d=probe_m2_bins_2d, bins_1d=90, ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins, theta_max=theta_max, theta_bins=theta_bins)
            hm = histogram_for_discriminator(disc, mixed_f, masks["mixed"], mm2_min=mm2_min, mm2_max=mm2_max, probe_m2_max=probe_m2_max, mm2_bins_2d=mm2_bins_2d, probe_m2_bins_2d=probe_m2_bins_2d, bins_1d=90, ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins, theta_max=theta_max, theta_bins=theta_bins)
            if (
                float(np.sum(hd)) >= min_data_count
                and float(np.sum(hp)) >= min_template_count
                and float(np.sum(hv)) >= min_template_count
                and float(np.sum(hm)) >= min_template_count
            ):
                hist4[disc] = (hd, hp, hv, hm)
            #endif
        #endfor
        if not hist4:
            continue
        #endif

        nuisance = {
            key: float(value)
            for key, value in coarse.items()
            if (key.endswith("_shift_bins") or key.endswith("_sigma_bins")) and np.isfinite(float(value))
        }
        two_fit = SharedMorphedFitResult(
            True,
            str(coarse.get("fit_message", "")),
            float(coarse["pi0_fraction"]),
            float(coarse.get("objective", np.nan)),
            float(coarse.get("poisson_deviance", np.nan)),
            int(coarse.get("ndof", 0)),
            nuisance,
        )
        fit3 = fit_shared_three_component_fixed_morph(hist4, two_fit)
        if not fit3.success:
            continue
        #endif

        d2 = float(coarse.get("deviance_per_ndof", np.nan))
        d3 = float(fit3.poisson_deviance / fit3.ndof) if fit3.ndof else float("nan")
        dpi0 = fit3.driver_pi0_fraction or {}
        dmix = fit3.driver_mixed_fraction or {}
        ppi0 = float(dpi0.get("mx2_ep_x_probe_m2", np.nan))
        tpi0 = float(dpi0.get("mx2_ep_x_pTmiss", np.nan))
        pmix = float(dmix.get("mx2_ep_x_probe_m2", np.nan))
        tmix = float(dmix.get("mx2_ep_x_pTmiss", np.nan))

        rows.append({
            "period": period.key,
            "label": period.label,
            "region": "FT",
            "energy_low_GeV": elo,
            "energy_high_GeV": ehi,
            "energy_center_GeV": 0.5 * (elo + ehi),
            "mix_shift": int(shift),
            "fit_model": "coarse_FT_three_component_fixed_two_component_morph",
            "two_component_pi0_fraction": float(coarse["pi0_fraction"]),
            "three_component_pi0_fraction": float(fit3.pi0_fraction),
            "three_component_dvcs_fraction": float(fit3.dvcs_fraction),
            "three_component_mixed_fraction": float(fit3.mixed_fraction),
            "delta_pi0_fraction_vs_two_component": float(fit3.pi0_fraction - float(coarse["pi0_fraction"])),
            "two_component_deviance_per_ndof": d2,
            "three_component_deviance_per_ndof": d3,
            "deviance_per_ndof_improvement": d2 - d3 if np.isfinite(d2) and np.isfinite(d3) else float("nan"),
            "probe_m2_driver_pi0_fraction": ppi0,
            "ptmiss_driver_pi0_fraction": tpi0,
            "driver_abs_delta_pi0_fraction": abs(ppi0 - tpi0) if np.isfinite(ppi0) and np.isfinite(tpi0) else float("nan"),
            "probe_m2_driver_mixed_fraction": pmix,
            "ptmiss_driver_mixed_fraction": tmix,
            "driver_abs_delta_mixed_fraction": abs(pmix - tmix) if np.isfinite(pmix) and np.isfinite(tmix) else float("nan"),
        })
        for crow in build_ft_coarse_three_component_closure(period, elo, ehi, hist4, two_fit, truth_points):
            crow["mix_shift"] = int(shift)
            closure_rows.append(crow)
        #endfor
    #endfor
    return rows, closure_rows


def summarize_ft_coarse_three_component_shifts(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    """Summarize median and deterministic-shift envelope of 3-component fits."""
    grouped: Dict[Tuple[str, float, float], List[Dict[str, object]]] = {}
    for row in rows:
        key = (str(row["period"]), float(row["energy_low_GeV"]), float(row["energy_high_GeV"]))
        grouped.setdefault(key, []).append(row)
    #endfor
    out: List[Dict[str, object]] = []
    fields = ("three_component_pi0_fraction", "three_component_mixed_fraction", "three_component_deviance_per_ndof", "deviance_per_ndof_improvement", "driver_abs_delta_pi0_fraction", "driver_abs_delta_mixed_fraction")
    for key, rr in grouped.items():
        first = rr[0]
        o: Dict[str, object] = {
            "period": key[0], "label": first["label"], "region": "FT",
            "energy_low_GeV": key[1], "energy_high_GeV": key[2],
            "energy_center_GeV": 0.5 * (key[1] + key[2]),
            "n_mix_shifts": len(rr),
            "two_component_pi0_fraction": float(first["two_component_pi0_fraction"]),
            "two_component_deviance_per_ndof": float(first["two_component_deviance_per_ndof"]),
        }
        for field_name in fields:
            x = np.asarray([float(r.get(field_name, np.nan)) for r in rr], dtype=float)
            x = x[np.isfinite(x)]
            o[f"median_{field_name}"] = float(np.median(x)) if x.size else float("nan")
            o[f"min_{field_name}"] = float(np.min(x)) if x.size else float("nan")
            o[f"max_{field_name}"] = float(np.max(x)) if x.size else float("nan")
        #endfor
        out.append(o)
    #endfor
    return out


def make_ft_coarse_three_component_canvas(
    period: PeriodConfig,
    summary_rows: List[Dict[str, object]],
    closure_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    if not summary_rows:
        return
    #endif
    rr = sorted(summary_rows, key=lambda r: float(r["energy_low_GeV"]))
    x = np.asarray([float(r["energy_center_GeV"]) for r in rr], dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(15.6, 10.0))

    axes[0,0].plot(x, [float(r["two_component_pi0_fraction"]) for r in rr], marker="o", linewidth=1.2, label=r"2-comp $f_{\pi^0}$")
    axes[0,0].plot(x, [float(r["median_three_component_pi0_fraction"]) for r in rr], marker="s", linewidth=1.2, label=r"3-comp median $f_{\pi^0}$")
    axes[0,0].fill_between(x, [float(r["min_three_component_pi0_fraction"]) for r in rr], [float(r["max_three_component_pi0_fraction"]) for r in rr], alpha=0.18, label="mix-shift envelope")
    axes[0,0].set_ylim(0,1); axes[0,0].set_ylabel(r"$f_{\pi^0}$"); axes[0,0].set_title("Composition shift from third component"); axes[0,0].grid(alpha=.2); axes[0,0].legend(frameon=False, fontsize=8)

    axes[0,1].plot(x, [float(r["median_three_component_mixed_fraction"]) for r in rr], marker="o", linewidth=1.2)
    axes[0,1].fill_between(x, [float(r["min_three_component_mixed_fraction"]) for r in rr], [float(r["max_three_component_mixed_fraction"]) for r in rr], alpha=.18)
    axes[0,1].set_ylim(0,1); axes[0,1].set_ylabel("mixed/wrong-photon fraction"); axes[0,1].set_title("Preferred uncorrelated-photon component"); axes[0,1].grid(alpha=.2)

    axes[1,0].plot(x, [float(r["two_component_deviance_per_ndof"]) for r in rr], marker="o", linewidth=1.2, label="2-component")
    axes[1,0].plot(x, [float(r["median_three_component_deviance_per_ndof"]) for r in rr], marker="s", linewidth=1.2, label="3-component")
    axes[1,0].set_ylabel("Poisson deviance / ndof"); axes[1,0].set_title("Fit-quality comparison"); axes[1,0].grid(alpha=.2); axes[1,0].legend(frameon=False, fontsize=8)

    grouped: Dict[Tuple[float,float], List[float]] = {}
    for row in closure_rows:
        if int(row.get("fit_success",0)) != 1: continue
        key = (float(row["energy_low_GeV"]), float(row["energy_high_GeV"]))
        v = float(row.get("max_abs_component_bias", np.nan))
        if np.isfinite(v): grouped.setdefault(key, []).append(v)
    closure_max = [max(grouped.get((float(r["energy_low_GeV"]), float(r["energy_high_GeV"])), [float("nan")])) for r in rr]
    axes[1,1].plot(x, [float(r["median_driver_abs_delta_pi0_fraction"]) for r in rr], marker="o", linewidth=1.2, label=r"driver $|\Delta f_{\pi^0}|$")
    axes[1,1].plot(x, closure_max, marker="s", linewidth=1.2, label="max 3-comp closure bias")
    axes[1,1].axhline(.01, color="black", linestyle=":", linewidth=1, label="1%")
    axes[1,1].set_ylabel("absolute fraction difference"); axes[1,1].set_title("Identifiability / driver consistency"); axes[1,1].grid(alpha=.2); axes[1,1].legend(frameon=False, fontsize=8)
    for ax in axes.ravel(): ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    fig.suptitle(f"{period.label}: coarse-FT 3-component wrong-photon diagnostic", fontsize=14)
    safe_finalize_figure(fig, Path(outdir)/"canvas_ft_coarse_three_component_diagnostic.png", rect=(0,0,1,.94)); plt.close(fig)

def run_stage2_shared_morphed_fits(period: PeriodConfig, data_f: Dict[str,np.ndarray], pi0_f: Dict[str,np.ndarray],
                                   dvcs_f: Dict[str,np.ndarray], pi0_control: Dict[str,Dict[str,float]],
                                   ft_theta_max: float, max_probe_energy: float, mm2_min: float, mm2_max: float,
                                   probe_m2_max: float, mm2_bins_2d: int, probe_m2_bins_2d: int,
                                   ptmiss_max: float, ptmiss_bins: int, theta_max: float, theta_bins: int,
                                   min_data_count: int, min_template_count: int, nuisance_shift_prior: float,
                                   nuisance_sigma_prior: float, max_shift_bins: float, max_sigma_bins: float,
                                   closure_truth_fractions: Sequence[float]) -> Tuple[List[Dict[str,object]], List[Dict[str,object]]]:
    edges=stage2_energy_edges(max_probe_energy); regions=["all","FT","FD_all"]+[f"FD_S{s}" for s in range(1,7)]
    rows=[]
    closure_rows: List[Dict[str,object]] = []
    for region in regions:
        for ib in range(len(edges)-1):
            elo,ehi=float(edges[ib]),float(edges[ib+1])
            masks={k:stage2_fit_mask(f,region,ft_theta_max,elo,ehi,mm2_min,mm2_max,probe_m2_max)
                   for k,f in (("data",data_f),("pi0",pi0_f),("dvcs",dvcs_f))}
            row={"period":period.key,"label":period.label,"beam_energy_GeV":period.beam_energy,"region":region,
                 "energy_low_GeV":elo,"energy_high_GeV":ehi,"energy_center_GeV":0.5*(elo+ehi),
                 "fit_model":"shared_morphed_2d",
                 "data_candidate_count":int(np.count_nonzero(masks["data"])),
                 "aaogen_candidate_count":int(np.count_nonzero(masks["pi0"])),
                 "dvcsgen_candidate_count":int(np.count_nonzero(masks["dvcs"]))}
            hist={}
            for disc in STAGE2_DRIVER_DISCRIMINATORS:
                if not discriminator_available(disc,data_f,pi0_f,dvcs_f): continue
                hd=histogram_for_discriminator(disc,data_f,masks["data"],mm2_min=mm2_min,mm2_max=mm2_max,probe_m2_max=probe_m2_max,
                    mm2_bins_2d=mm2_bins_2d,probe_m2_bins_2d=probe_m2_bins_2d,bins_1d=90,ptmiss_max=ptmiss_max,ptmiss_bins=ptmiss_bins,theta_max=theta_max,theta_bins=theta_bins)
                hp=histogram_for_discriminator(disc,pi0_f,masks["pi0"],mm2_min=mm2_min,mm2_max=mm2_max,probe_m2_max=probe_m2_max,
                    mm2_bins_2d=mm2_bins_2d,probe_m2_bins_2d=probe_m2_bins_2d,bins_1d=90,ptmiss_max=ptmiss_max,ptmiss_bins=ptmiss_bins,theta_max=theta_max,theta_bins=theta_bins)
                hv=histogram_for_discriminator(disc,dvcs_f,masks["dvcs"],mm2_min=mm2_min,mm2_max=mm2_max,probe_m2_max=probe_m2_max,
                    mm2_bins_2d=mm2_bins_2d,probe_m2_bins_2d=probe_m2_bins_2d,bins_1d=90,ptmiss_max=ptmiss_max,ptmiss_bins=ptmiss_bins,theta_max=theta_max,theta_bins=theta_bins)
                if np.sum(hd)>=min_data_count and np.sum(hp)>=min_template_count and np.sum(hv)>=min_template_count: hist[disc]=(hd,hp,hv)
            #endfor
            fit=fit_shared_morphed_composition(hist,pi0_control,nuisance_shift_prior,nuisance_sigma_prior,max_shift_bins,max_sigma_bins)
            row.update({"fit_success":int(fit.success),"fit_message":fit.message,"n_drivers":len(hist),"pi0_fraction":fit.pi0_fraction,
                        "dvcs_fraction":1.0-fit.pi0_fraction if np.isfinite(fit.pi0_fraction) else float("nan"),"objective":fit.objective,
                        "poisson_deviance":fit.poisson_deviance,"ndof":fit.ndof,"deviance_per_ndof":fit.poisson_deviance/fit.ndof if fit.ndof else float("nan")})
            if fit.nuisance: row.update(fit.nuisance)
            # Meaningful deterministic closure: inject the *real-data best-fit
            # morphed shapes* at controlled pi0 fractions and rerun the complete
            # Stage-II shared fit from its normal initialization.  Restrict this
            # to FT and FD_all to directly test the two production topologies
            # without multiplying the sector-level runtime.
            if region in ("FT", "FD_all") and hist and fit.success:
                closure_rows.extend(
                    build_stage2_template_mixture_closure(
                        period,
                        region,
                        elo,
                        ehi,
                        hist,
                        fit,
                        pi0_control,
                        closure_truth_fractions,
                        nuisance_shift_prior,
                        nuisance_sigma_prior,
                        max_shift_bins,
                        max_sigma_bins,
                    )
                )
            #endif
            rows.append(row)
        #endfor
    #endfor
    return rows, closure_rows




def _profile_fraction_for_templates(
    data_hists: Dict[str, np.ndarray],
    pi0_templates: Dict[str, np.ndarray],
    dvcs_templates: Dict[str, np.ndarray],
) -> Tuple[float, float]:
    """
    Profile only the shared pi0 fraction for fixed morphed templates.

    A bounded one-dimensional minimization is substantially faster and more
    robust here than running L-BFGS-B on a logit parameter thousands of times
    during nuisance scans.
    """
    names = [n for n in STAGE2_DRIVER_DISCRIMINATORS if n in data_hists]
    if not names:
        return float("nan"), float("nan")
    #endif

    prepared = []
    for name in names:
        data = np.asarray(data_hists[name], dtype=float).ravel()
        tp = np.asarray(pi0_templates[name], dtype=float).ravel()
        td = np.asarray(dvcs_templates[name], dtype=float).ravel()
        nd = float(np.sum(data))
        prepared.append((data, tp, td, nd))
    #endfor

    def objective_f(f: float) -> float:
        ff = float(np.clip(f, 3.0e-4, 1.0 - 3.0e-4))
        total = 0.0
        for data, tp, td, nd in prepared:
            mu = np.clip(
                nd * (ff * tp + (1.0 - ff) * td),
                1.0e-12,
                None,
            )
            total += float(np.sum(mu - data * np.log(mu))) / len(prepared)
        #endfor
        return total
    #enddef

    result = minimize_scalar(
        objective_f,
        bounds=(3.0e-4, 1.0 - 3.0e-4),
        method="bounded",
        options={"xatol": 1.0e-7, "maxiter": 120},
    )
    return float(result.x), float(result.fun)



def run_stage2_morph_profile_scans(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
    min_data_count: int,
    min_template_count: int,
    shift_grid: np.ndarray,
    sigma_grid: np.ndarray,
) -> List[Dict[str, object]]:
    """
    Explicit one-nuisance-at-a-time profile scans in the all-region sample.

    At every scan point:
      * one morph nuisance is fixed to the requested value;
      * all other morph nuisances are fixed at zero;
      * the shared pi0 fraction is profiled.

    This directly diagnoses whether zero is genuinely preferred by the
    likelihood or whether the multidimensional L-BFGS-B fit is failing to move.
    """
    edges = stage2_energy_edges(max_probe_energy)
    rows: List[Dict[str, object]] = []

    for ib in range(len(edges) - 1):
        elo = float(edges[ib])
        ehi = float(edges[ib + 1])

        masks = {
            key: stage2_fit_mask(
                feat,
                "all",
                ft_theta_max,
                elo,
                ehi,
                mm2_min,
                mm2_max,
                probe_m2_max,
            )
            for key, feat in (
                ("data", data_f),
                ("pi0", pi0_f),
                ("dvcs", dvcs_f),
            )
        }

        hist: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for disc in STAGE2_DRIVER_DISCRIMINATORS:
            if not discriminator_available(disc, data_f, pi0_f, dvcs_f):
                continue
            #endif
            hd = histogram_for_discriminator(
                disc,
                data_f,
                masks["data"],
                mm2_min=mm2_min,
                mm2_max=mm2_max,
                probe_m2_max=probe_m2_max,
                mm2_bins_2d=mm2_bins_2d,
                probe_m2_bins_2d=probe_m2_bins_2d,
                bins_1d=90,
                ptmiss_max=ptmiss_max,
                ptmiss_bins=ptmiss_bins,
                theta_max=theta_max,
                theta_bins=theta_bins,
            )
            hp = histogram_for_discriminator(
                disc,
                pi0_f,
                masks["pi0"],
                mm2_min=mm2_min,
                mm2_max=mm2_max,
                probe_m2_max=probe_m2_max,
                mm2_bins_2d=mm2_bins_2d,
                probe_m2_bins_2d=probe_m2_bins_2d,
                bins_1d=90,
                ptmiss_max=ptmiss_max,
                ptmiss_bins=ptmiss_bins,
                theta_max=theta_max,
                theta_bins=theta_bins,
            )
            hv = histogram_for_discriminator(
                disc,
                dvcs_f,
                masks["dvcs"],
                mm2_min=mm2_min,
                mm2_max=mm2_max,
                probe_m2_max=probe_m2_max,
                mm2_bins_2d=mm2_bins_2d,
                probe_m2_bins_2d=probe_m2_bins_2d,
                bins_1d=90,
                ptmiss_max=ptmiss_max,
                ptmiss_bins=ptmiss_bins,
                theta_max=theta_max,
                theta_bins=theta_bins,
            )
            if (
                np.sum(hd) >= min_data_count
                and np.sum(hp) >= min_template_count
                and np.sum(hv) >= min_template_count
            ):
                hist[disc] = (hd, hp, hv)
            #endif
        #endfor

        if not hist:
            continue
        #endif

        data_hists = {name: h[0] for name, h in hist.items()}
        base_pi0 = {
            name: normalized_template(h[1]).reshape(h[1].shape)
            for name, h in hist.items()
        }
        base_dvcs = {
            name: normalized_template(h[2]).reshape(h[2].shape)
            for name, h in hist.items()
        }

        # Reference objective at every nuisance set to zero.
        f0, obj0 = _profile_fraction_for_templates(
            data_hists,
            base_pi0,
            base_dvcs,
        )

        for disc in hist.keys():
            for component in ("pi0", "dvcs"):
                for nuisance, grid in (
                    ("shift", shift_grid),
                    ("sigma", sigma_grid),
                ):
                    for value in np.asarray(grid, dtype=float):
                        tp = {name: arr.copy() for name, arr in base_pi0.items()}
                        td = {name: arr.copy() for name, arr in base_dvcs.items()}

                        if component == "pi0":
                            tp[disc] = morph_template_second_axis(
                                hist[disc][1],
                                float(value) if nuisance == "shift" else 0.0,
                                float(value) if nuisance == "sigma" else 0.0,
                            )
                        else:
                            td[disc] = morph_template_second_axis(
                                hist[disc][2],
                                float(value) if nuisance == "shift" else 0.0,
                                float(value) if nuisance == "sigma" else 0.0,
                            )
                        #endif

                        f_prof, objective = _profile_fraction_for_templates(
                            data_hists,
                            tp,
                            td,
                        )
                        rows.append({
                            "period": period.key,
                            "label": period.label,
                            "region": "all",
                            "energy_low_GeV": elo,
                            "energy_high_GeV": ehi,
                            "energy_center_GeV": 0.5 * (elo + ehi),
                            "driver": disc,
                            "component": component,
                            "nuisance": nuisance,
                            "fixed_value_bins": float(value),
                            "profiled_pi0_fraction": f_prof,
                            "objective": objective,
                            "delta_objective_from_all_zero": (
                                objective - obj0
                                if np.isfinite(objective) and np.isfinite(obj0)
                                else float("nan")
                            ),
                            "all_zero_profiled_pi0_fraction": f0,
                            "all_zero_objective": obj0,
                        })
                    #endfor
                #endfor
            #endfor
        #endfor
    #endfor

    return rows


def make_morph_profile_canvas(
    period: PeriodConfig,
    profile_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Visualize the most important nuisance profile:
    aaogen pTmiss additional smearing.

    Left: Δ objective vs fixed smearing for every Stage-II energy bin.
    Right: profile-preferred smearing and associated f_pi0 vs energy.
    """
    target = [
        r for r in profile_rows
        if r["driver"] == "mx2_ep_x_pTmiss"
        and r["component"] == "pi0"
        and r["nuisance"] == "sigma"
        and np.isfinite(r["delta_objective_from_all_zero"])
    ]
    if not target:
        return
    #endif

    groups: Dict[Tuple[float, float], List[Dict[str, object]]] = {}
    for row in target:
        key = (
            float(row["energy_low_GeV"]),
            float(row["energy_high_GeV"]),
        )
        groups.setdefault(key, []).append(row)
    #endfor

    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6))
    best_rows = []

    for (elo, ehi), rr in sorted(groups.items()):
        rr = sorted(rr, key=lambda r: float(r["fixed_value_bins"]))
        x = np.asarray([r["fixed_value_bins"] for r in rr], dtype=float)
        y = np.asarray([r["delta_objective_from_all_zero"] for r in rr], dtype=float)
        # Display relative to the best scanned point so the visual minimum is clear.
        y_rel = y - np.nanmin(y)
        axes[0].plot(
            x,
            y_rel,
            marker="o",
            ms=3,
            label=f"{elo:g}–{ehi:g} GeV",
        )
        best = rr[int(np.nanargmin(y))]
        best_rows.append(best)
    #endfor

    axes[0].axvline(0.0, linestyle="--", linewidth=0.9)
    axes[0].set_xlabel("Fixed aaogen $p_{T,\\mathrm{miss}}$ smearing (bins)")
    axes[0].set_ylabel(r"$\Delta$ objective from scan minimum")
    axes[0].set_title("Explicit nuisance profile")
    axes[0].legend(fontsize=7, ncol=2)
    axes[0].grid(alpha=0.18)

    best_rows = sorted(best_rows, key=lambda r: float(r["energy_center_GeV"]))
    x_energy = np.asarray(
        [r["energy_center_GeV"] for r in best_rows],
        dtype=float,
    )
    axes[1].plot(
        x_energy,
        [r["fixed_value_bins"] for r in best_rows],
        marker="o",
        label="Preferred smearing",
    )
    axes[1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1].set_ylabel("Profile-preferred smearing (bins)")
    axes[1].set_title("Best scanned aaogen $p_{T,\\mathrm{miss}}$ smearing")
    axes[1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: shared-fit morphing diagnostic",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_morph_nuisance_profiles.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)


def make_shared_fit_canvas(period: PeriodConfig, shared_rows: List[Dict[str,object]], raw_rows: List[Dict[str,object]],
                           pi0_control: Dict[str,object], outdir: Path) -> None:
    """Compact shared-fit summary; reconstructed-pi0 control is diagnostic only."""
    rr=[r for r in shared_rows if r["region"]=="all" and int(r.get("fit_success",0))]
    if not rr: return
    rr=sorted(rr,key=lambda r:float(r["energy_center_GeV"])); x=np.asarray([r["energy_center_GeV"] for r in rr])
    fig,ax=plt.subplots(2,2,figsize=(13.5,9.5))
    ax[0,0].plot(x,[r["pi0_fraction"] for r in rr],marker="o",label="shared morphed fit")
    for disc,label in (("mx2_ep_x_probe_m2", r"raw $M_X^2\otimes M_{\mathrm{probe}}^2$"),("mx2_ep_x_pTmiss", r"raw $M_X^2\otimes p_{T,\mathrm{miss}}$")):
        q=sorted([r for r in raw_rows if r["region"]=="all" and r["discriminator"]==disc and int(r["fit_success"])],key=lambda r:float(r["energy_center_GeV"]))
        if q: ax[0,0].plot([r["energy_center_GeV"] for r in q],[r["pi0_fraction"] for r in q],marker=".",linestyle="--",label=label)
    #endfor
    ax[0,0].set_ylabel(r"$f_{\pi^0}$"); ax[0,0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"); ax[0,0].set_ylim(0,1); ax[0,0].legend(fontsize=8); ax[0,0].grid(alpha=.18)
    ax[0,1].plot(x,[r["deviance_per_ndof"] for r in rr],marker="o")
    ax[0,1].set_ylabel("combined Poisson deviance / ndof"); ax[0,1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"); ax[0,1].grid(alpha=.18)
    for key,label in (
        ("mx2_ep_x_probe_m2_pi0_shift_bins", "probe-M2 pi0 shift"),
        ("mx2_ep_x_probe_m2_pi0_sigma_bins", "probe-M2 pi0 smear"),
        ("mx2_ep_x_probe_m2_dvcs_shift_bins", "probe-M2 DVCS shift"),
        ("mx2_ep_x_probe_m2_dvcs_sigma_bins", "probe-M2 DVCS smear"),
        ("mx2_ep_x_pTmiss_pi0_shift_bins", "pTmiss pi0 shift"),
        ("mx2_ep_x_pTmiss_pi0_sigma_bins", "pTmiss pi0 smear"),
        ("mx2_ep_x_pTmiss_dvcs_shift_bins", "pTmiss DVCS shift"),
        ("mx2_ep_x_pTmiss_dvcs_sigma_bins", "pTmiss DVCS smear"),
    ):
        values=np.asarray([r.get(key,np.nan) for r in rr],dtype=float)
        if np.any(np.isfinite(values)): ax[1,0].plot(x,values,marker=".",label=label)
        #endif
    #endfor
    ax[1,0].axhline(0,linewidth=.8); ax[1,0].set_ylabel("template morph nuisance (bins)"); ax[1,0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"); ax[1,0].legend(fontsize=7); ax[1,0].grid(alpha=.18)
    for truth in (20,50,80):
        values=np.asarray([r.get(f"closure_bias_{truth}",np.nan) for r in rr],dtype=float)
        if np.any(np.isfinite(values)): ax[1,1].plot(x,values,marker=".",label=f"truth {truth/100:.1f}")
        #endif
    #endfor
    ax[1,1].axhline(0,linewidth=.8); ax[1,1].set_ylabel(r"closure $f_{\mathrm{fit}}-f_{\mathrm{true}}$"); ax[1,1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"); handles_11, labels_11 = ax[1,1].get_legend_handles_labels(); ax[1,1].legend(fontsize=8) if handles_11 else None; ax[1,1].grid(alpha=.18)
    core=pi0_control.get("robust_core_global",{}) if isinstance(pi0_control,dict) else {}
    sadd=core.get("sigma_add_core_GeV",float("nan")) if isinstance(core,dict) else float("nan")
    fig.suptitle(
        f"{period.label}: shared morphed Stage-II composition (each FT/FD region fitted independently); "
        f"eppi0 control sigma_add(core)={sadd:.4f} GeV; control prior NOT applied",
        fontsize=12.2,
    )
    safe_finalize_figure(fig,outdir/"canvas_shared_morphed_composition.png",rect=(0,0,1,.95)); plt.close(fig)



def make_ft_fd_composition_canvas(
    period: PeriodConfig,
    shared_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Make the fitted composition distinction between FT and FD unmistakable.

    Left panel:
      FT pi0 and BH/DVCS fractions.

    Right panel:
      FD-all BH/DVCS fraction plus the six individual FD-sector BH/DVCS
      fractions. Since BH/DVCS and pi0 sum to unity, this directly displays the
      expected larger BH/DVCS contribution in FT without conflating it with a
      template-shape problem.
    """
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.7))

    ft = sorted(
        [
            r for r in shared_rows
            if r["region"] == "FT"
            and int(r["fit_success"]) == 1
            and np.isfinite(r["pi0_fraction"])
        ],
        key=lambda r: float(r["energy_center_GeV"]),
    )

    if ft:
        x = [r["energy_center_GeV"] for r in ft]
        fpi0 = [r["pi0_fraction"] for r in ft]
        fbh = [r["dvcs_fraction"] for r in ft]
        axes[0].plot(x, fbh, marker="o", label="BH/DVCS")
        axes[0].plot(x, fpi0, marker="o", label=r"$\pi^0$")
    #endif
    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0].set_ylabel("Fitted fraction")
    axes[0].set_title("FT: independent fitted composition")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.18)

    fd_all = sorted(
        [
            r for r in shared_rows
            if r["region"] == "FD_all"
            and int(r["fit_success"]) == 1
            and np.isfinite(r["dvcs_fraction"])
        ],
        key=lambda r: float(r["energy_center_GeV"]),
    )
    if fd_all:
        axes[1].plot(
            [r["energy_center_GeV"] for r in fd_all],
            [r["dvcs_fraction"] for r in fd_all],
            marker="o",
            linewidth=2.0,
            label="FD all",
        )
    #endif

    for sector in range(1, 7):
        region = f"FD_S{sector}"
        rr = sorted(
            [
                r for r in shared_rows
                if r["region"] == region
                and int(r["fit_success"]) == 1
                and np.isfinite(r["dvcs_fraction"])
            ],
            key=lambda r: float(r["energy_center_GeV"]),
        )
        if rr:
            axes[1].plot(
                [r["energy_center_GeV"] for r in rr],
                [r["dvcs_fraction"] for r in rr],
                marker=".",
                linewidth=1.0,
                label=f"S{sector}",
            )
        #endif
    #endfor

    axes[1].set_ylim(0.0, 1.02)
    axes[1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1].set_ylabel("Fitted BH/DVCS fraction")
    axes[1].set_title("FD: BH/DVCS fraction by sector")
    axes[1].legend(fontsize=8, ncol=2)
    axes[1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: Stage-II fitted composition by predicted detector region",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_ft_fd_composition.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)


def shared_driver_projection(
    region: str,
    row: Dict[str, object],
    discriminator: str,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
    display_bin_factor: int = 1,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """
    Build a one-dimensional second-axis projection using the actual shared-fit
    region and fraction.

    Returns data, final morphed pi0/BH-DVCS components, the total final model,
    and the corresponding *starting/unmorphed* pi0/BH-DVCS components.  The raw
    templates make it visually obvious what the morphing changed.
    """
    display_bin_factor = max(1, int(display_bin_factor))
    display_ptmiss_bins = int(ptmiss_bins) * display_bin_factor
    display_probe_m2_bins = int(probe_m2_bins_2d) * display_bin_factor

    elo = float(row["energy_low_GeV"])
    ehi = float(row["energy_high_GeV"])

    masks = {
        key: stage2_fit_mask(
            feat,
            region,
            ft_theta_max,
            elo,
            ehi,
            mm2_min,
            mm2_max,
            probe_m2_max,
        )
        for key, feat in (
            ("data", data_f),
            ("pi0", pi0_f),
            ("dvcs", dvcs_f),
        )
    }

    hd = histogram_for_discriminator(
        discriminator,
        data_f,
        masks["data"],
        mm2_min=mm2_min,
        mm2_max=mm2_max,
        probe_m2_max=probe_m2_max,
        mm2_bins_2d=mm2_bins_2d,
        probe_m2_bins_2d=display_probe_m2_bins,
        bins_1d=90,
        ptmiss_max=ptmiss_max,
        ptmiss_bins=display_ptmiss_bins,
        theta_max=theta_max,
        theta_bins=theta_bins,
    )
    hp = histogram_for_discriminator(
        discriminator,
        pi0_f,
        masks["pi0"],
        mm2_min=mm2_min,
        mm2_max=mm2_max,
        probe_m2_max=probe_m2_max,
        mm2_bins_2d=mm2_bins_2d,
        probe_m2_bins_2d=display_probe_m2_bins,
        bins_1d=90,
        ptmiss_max=ptmiss_max,
        ptmiss_bins=display_ptmiss_bins,
        theta_max=theta_max,
        theta_bins=theta_bins,
    )
    hv = histogram_for_discriminator(
        discriminator,
        dvcs_f,
        masks["dvcs"],
        mm2_min=mm2_min,
        mm2_max=mm2_max,
        probe_m2_max=probe_m2_max,
        mm2_bins_2d=mm2_bins_2d,
        probe_m2_bins_2d=display_probe_m2_bins,
        bins_1d=90,
        ptmiss_max=ptmiss_max,
        ptmiss_bins=display_ptmiss_bins,
        theta_max=theta_max,
        theta_bins=theta_bins,
    )

    raw_pi0_shape = normalized_template(hp).reshape(hp.shape)
    raw_dvcs_shape = normalized_template(hv).reshape(hv.shape)

    tp = morph_template_second_axis(
        hp,
        display_bin_factor
        * float(row.get(f"{discriminator}_pi0_shift_bins", 0.0)),
        display_bin_factor
        * float(row.get(f"{discriminator}_pi0_sigma_bins", 0.0)),
    )
    td = morph_template_second_axis(
        hv,
        display_bin_factor
        * float(row.get(f"{discriminator}_dvcs_shift_bins", 0.0)),
        display_bin_factor
        * float(row.get(f"{discriminator}_dvcs_sigma_bins", 0.0)),
    )

    data_proj = np.sum(hd, axis=0).astype(float)
    pi0_shape = np.sum(tp, axis=0).astype(float)
    dvcs_shape = np.sum(td, axis=0).astype(float)
    raw_pi0_proj = np.sum(raw_pi0_shape, axis=0).astype(float)
    raw_dvcs_proj = np.sum(raw_dvcs_shape, axis=0).astype(float)

    pi0_shape /= max(np.sum(pi0_shape), 1.0e-30)
    dvcs_shape /= max(np.sum(dvcs_shape), 1.0e-30)
    raw_pi0_proj /= max(np.sum(raw_pi0_proj), 1.0e-30)
    raw_dvcs_proj /= max(np.sum(raw_dvcs_proj), 1.0e-30)

    nd = float(np.sum(data_proj))
    fpi0 = float(row["pi0_fraction"])
    pi0_component = nd * fpi0 * pi0_shape
    dvcs_component = nd * (1.0 - fpi0) * dvcs_shape
    raw_pi0_component = nd * fpi0 * raw_pi0_proj
    raw_dvcs_component = nd * (1.0 - fpi0) * raw_dvcs_proj
    model = pi0_component + dvcs_component

    if discriminator == "mx2_ep_x_pTmiss":
        edges = np.linspace(0.0, ptmiss_max, display_ptmiss_bins + 1)
    else:
        edges = np.linspace(
            -probe_m2_max,
            probe_m2_max,
            display_probe_m2_bins + 1,
        )
    #endif
    centers = 0.5 * (edges[:-1] + edges[1:])

    return (
        centers,
        data_proj,
        pi0_component,
        dvcs_component,
        model,
        raw_pi0_component,
        raw_dvcs_component,
    )


def make_ft_fd_fit_overlay_canvas(
    period: PeriodConfig,
    shared_rows: List[Dict[str, object]],
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    outdir: Path,
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins_2d: int,
    probe_m2_bins_2d: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
) -> None:
    """
    Explicitly show what is being fit in FT versus FD.

    Rows:
      top    = FT
      bottom = all FD sectors combined

    Columns:
      representative predicted-probe energy bins.

    The pTmiss projection is used because it currently provides the clearest
    shape leverage and is where morphing is most strongly preferred.
    """
    targets = ((0.40, 0.50), (0.80, 1.00), (1.50, 2.00))
    fig, axes = plt.subplots(2, 3, figsize=(18.5, 10.2))

    for irow, region in enumerate(("FT", "FD_all")):
        for icol, (elo, ehi) in enumerate(targets):
            ax = axes[irow, icol]
            candidates = [
                r for r in shared_rows
                if r["region"] == region
                and int(r["fit_success"]) == 1
                and abs(float(r["energy_low_GeV"]) - elo) < 1.0e-9
                and abs(float(r["energy_high_GeV"]) - ehi) < 1.0e-9
            ]
            if not candidates:
                ax.text(
                    0.5,
                    0.5,
                    "No successful fit",
                    ha="center",
                    va="center",
                )
                ax.set_axis_off()
                continue
            #endif

            row = candidates[0]
            (
                centers,
                data,
                pi0_c,
                dvcs_c,
                model,
                raw_pi0_c,
                raw_dvcs_c,
            ) = shared_driver_projection(
                region,
                row,
                "mx2_ep_x_pTmiss",
                data_f,
                pi0_f,
                dvcs_f,
                ft_theta_max,
                mm2_min,
                mm2_max,
                probe_m2_max,
                mm2_bins_2d,
                probe_m2_bins_2d,
                ptmiss_max,
                ptmiss_bins,
                theta_max,
                theta_bins,
                display_bin_factor=2,
            )

            if not np.allclose(
                model,
                dvcs_c + pi0_c,
                rtol=1.0e-12,
                atol=1.0e-10,
            ):
                raise RuntimeError(
                    f"{period.label} {region} {elo:.2f}-{ehi:.2f} GeV: "
                    "display-model component closure failed."
                )
            #endif

            # Match the visual language of the main exclusivity fitter:
            # black points = data; blue = BH/DVCS; red = pi0; green = total fit.
            # Thin dashed curves show the starting/unmorphed templates; solid
            # curves show the final morphed components used in the shared fit.
            ax.errorbar(
                centers,
                data,
                yerr=np.sqrt(np.maximum(data, 1.0)),
                fmt="o",
                ms=2.8,
                linewidth=0.8,
                capsize=0.0,
                color=SAMPLE_COLORS["data"],
                label="Data",
                zorder=6,
            )
            ax.step(
                centers,
                raw_dvcs_c,
                where="mid",
                color=SAMPLE_COLORS["dvcs_mc"],
                linewidth=0.75,
                linestyle="--",
                alpha=0.75,
                label="BH/DVCS pre-morph shape (reference)",
                zorder=1,
            )
            ax.step(
                centers,
                raw_pi0_c,
                where="mid",
                color=SAMPLE_COLORS["pi0_mc"],
                linewidth=0.75,
                linestyle="--",
                alpha=0.75,
                label=r"$\pi^0$ pre-morph shape (reference)",
                zorder=1,
            )
            ax.step(
                centers,
                dvcs_c,
                where="mid",
                color=SAMPLE_COLORS["dvcs_mc"],
                linewidth=1.6,
                label="BH/DVCS fitted component",
                zorder=3,
            )
            ax.step(
                centers,
                pi0_c,
                where="mid",
                color=SAMPLE_COLORS["pi0_mc"],
                linewidth=1.6,
                label=r"$\pi^0$ fitted component",
                zorder=3,
            )
            ax.step(
                centers,
                model,
                where="mid",
                color=SAMPLE_COLORS["fit"],
                linewidth=2.0,
                label=r"total fit = fitted BH/DVCS + fitted $\pi^0$",
                zorder=4,
            )
            ax.set_xlabel(r"$p_{T,\mathrm{miss}}$ (GeV)")
            ax.set_ylabel("Entries / bin")
            if region == "FT":
                ax.set_xlim(0.0, 0.30)
            else:
                ax.set_xlim(0.0, ptmiss_max)
            #endif

            fbh = float(row["dvcs_fraction"])
            fpi0 = float(row["pi0_fraction"])
            smear_pi0 = float(
                row.get("mx2_ep_x_pTmiss_pi0_sigma_bins", 0.0)
            )
            smear_dvcs = float(
                row.get("mx2_ep_x_pTmiss_dvcs_sigma_bins", 0.0)
            )
            ax.set_title(
                f"{region}, {elo:.2f}–{ehi:.2f} GeV\n"
                f"BH/DVCS={fbh:.2f}, $\\pi^0$={fpi0:.2f}; "
                f"smear=({smear_dvcs:.1f},{smear_pi0:.1f}) bins"
            )
            ax.text(
                0.985,
                0.025,
                "solid blue + solid red = green\n"
                "dashed curves are pre-morph references only",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=6.7,
                alpha=0.75,
            )
            ax.grid(alpha=0.20)
        #endfor
    #endfor

    axes[0, 0].legend(fontsize=7.5, frameon=False, ncol=2)
    fig.suptitle(
        f"{period.label}: explicit FT versus FD Stage-II $p_{{T,\\mathrm{{miss}}}}$ fits",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_ft_fd_fit_overlays.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)


def closure_rows_for_shared_fit(period: PeriodConfig, shared_rows: List[Dict[str,object]]) -> List[Dict[str,object]]:
    """Lightweight bookkeeping closure target grid; full event bootstrap is deferred until final extraction."""
    rows=[]
    for r in shared_rows:
        if r["region"] != "all" or not int(r["fit_success"]): continue
        for truth in (0.2,0.5,0.8):
            rows.append({"period":period.key,"region":r["region"],"energy_low_GeV":r["energy_low_GeV"],"energy_high_GeV":r["energy_high_GeV"],
                         "closure_truth_fraction":truth,"note":"reserved target for pseudo-data closure in next statistical-validation step"})
        #endfor
    #endfor
    return rows

def discriminator_spread_rows(
    fit_rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    """
    Collapse the discriminator fit table into one stability row per
    period/region/energy bin.

    Spread quantities exclude unavailable/failed discriminators.
    """
    groups: Dict[Tuple[str, str, float, float], List[Dict[str, object]]] = {}
    for row in fit_rows:
        key = (
            str(row["period"]),
            str(row["region"]),
            float(row["energy_low_GeV"]),
            float(row["energy_high_GeV"]),
        )
        groups.setdefault(key, []).append(row)
    #endfor

    out: List[Dict[str, object]] = []
    for key, rows in groups.items():
        good = [
            r for r in rows
            if int(r["fit_success"]) == 1
            and np.isfinite(r["pi0_fraction"])
        ]
        nominal = next(
            (
                r for r in good
                if str(r["discriminator"]) == STAGE2_NOMINAL_DISCRIMINATOR
            ),
            None,
        )

        values = np.asarray(
            [float(r["pi0_fraction"]) for r in good],
            dtype=float,
        )
        if values.size:
            fmin = float(np.min(values))
            fmax = float(np.max(values))
            spread = fmax - fmin
            half_range = 0.5 * spread
            std = float(np.std(values))
        else:
            fmin = fmax = spread = half_range = std = float("nan")
        #endif

        row0 = rows[0]
        out.append({
            "period": key[0],
            "label": row0["label"],
            "region": key[1],
            "energy_low_GeV": key[2],
            "energy_high_GeV": key[3],
            "energy_center_GeV": float(row0["energy_center_GeV"]),
            "n_successful_discriminators": len(good),
            "nominal_pi0_fraction": (
                float(nominal["pi0_fraction"]) if nominal is not None else float("nan")
            ),
            "min_pi0_fraction": fmin,
            "max_pi0_fraction": fmax,
            "full_spread_pi0_fraction": spread,
            "half_range_pi0_fraction": half_range,
            "std_pi0_fraction": std,
            "nominal_deviance_per_ndof": (
                float(nominal["deviance_per_ndof"])
                if nominal is not None else float("nan")
            ),
        })
    #endfor

    return out


def run_mixed_shift_diagnostic(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    mixed_f: Dict[str, np.ndarray],
    fit_rows: List[Dict[str, object]],
    shift: int,
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    mm2_bins: int,
    probe_m2_bins: int,
    min_template_count: int,
) -> List[Dict[str, object]]:
    """
    Repeat only the nominal 2D discriminator fit with one extra mixed-data
    component. This is a stress test, not a production denominator result.
    """
    nominal_rows = [
        r for r in fit_rows
        if str(r["discriminator"]) == STAGE2_NOMINAL_DISCRIMINATOR
        and int(r["fit_success"]) == 1
    ]

    mm2_edges = np.linspace(mm2_min, mm2_max, mm2_bins + 1)
    pm2_edges = np.linspace(-probe_m2_max, probe_m2_max, probe_m2_bins + 1)

    rows: List[Dict[str, object]] = []
    for nominal in nominal_rows:
        region = str(nominal["region"])
        elo = float(nominal["energy_low_GeV"])
        ehi = float(nominal["energy_high_GeV"])

        md = stage2_fit_mask(data_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)
        mp = stage2_fit_mask(pi0_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)
        mv = stage2_fit_mask(dvcs_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)
        mm = stage2_fit_mask(mixed_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)

        if int(np.count_nonzero(mm)) < min_template_count:
            continue
        #endif

        hd, _, _ = np.histogram2d(
            data_f["ep_missing_mass2"][md],
            data_f["pred_probe_mass2"][md],
            bins=(mm2_edges, pm2_edges),
        )
        hp, _, _ = np.histogram2d(
            pi0_f["ep_missing_mass2"][mp],
            pi0_f["pred_probe_mass2"][mp],
            bins=(mm2_edges, pm2_edges),
        )
        hv, _, _ = np.histogram2d(
            dvcs_f["ep_missing_mass2"][mv],
            dvcs_f["pred_probe_mass2"][mv],
            bins=(mm2_edges, pm2_edges),
        )
        hm, _, _ = np.histogram2d(
            mixed_f["ep_missing_mass2"][mm],
            mixed_f["pred_probe_mass2"][mm],
            bins=(mm2_edges, pm2_edges),
        )

        fit = fit_three_component_diagnostic(hd, hp, hv, hm)
        if not fit.success:
            continue
        #endif

        nominal_pi0 = float(nominal["pi0_fraction"])
        rows.append({
            "period": period.key,
            "label": period.label,
            "region": region,
            "energy_low_GeV": elo,
            "energy_high_GeV": ehi,
            "energy_center_GeV": float(nominal["energy_center_GeV"]),
            "mix_shift": int(shift),
            "nominal_pi0_fraction": nominal_pi0,
            "diagnostic_pi0_fraction": fit.pi0_fraction,
            "diagnostic_dvcs_fraction": fit.dvcs_fraction,
            "diagnostic_mixed_fraction": fit.mixed_fraction,
            "delta_pi0_fraction": fit.pi0_fraction - nominal_pi0,
        })
    #endfor

    return rows


def representative_fit_projection(
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    row: Dict[str, object],
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    bins: int = 90,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Return M_X^2(ep) data and fitted component projections for one nominal
    mx2_ep_x_probe_m2 fit row.
    """
    region = str(row["region"])
    elo = float(row["energy_low_GeV"])
    ehi = float(row["energy_high_GeV"])

    md = stage2_fit_mask(data_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)
    mp = stage2_fit_mask(pi0_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)
    mv = stage2_fit_mask(dvcs_f, region, ft_theta_max, elo, ehi, mm2_min, mm2_max, probe_m2_max)

    edges = np.linspace(mm2_min, mm2_max, bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hd, _ = np.histogram(data_f["ep_missing_mass2"][md], bins=edges)
    hp, _ = np.histogram(pi0_f["ep_missing_mass2"][mp], bins=edges)
    hv, _ = np.histogram(dvcs_f["ep_missing_mass2"][mv], bins=edges)

    hp = hp.astype(float)
    hv = hv.astype(float)
    hp /= max(float(np.sum(hp)), 1.0)
    hv /= max(float(np.sum(hv)), 1.0)

    pi0_component = float(row["pi0_yield"]) * hp
    dvcs_component = float(row["dvcs_yield"]) * hv
    model = pi0_component + dvcs_component
    return centers, hd.astype(float), pi0_component, dvcs_component, model


def unit_hist(ax, values: np.ndarray, bins: np.ndarray, label: str) -> None:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    weights = np.ones(len(vals), dtype=float) / max(len(vals), 1)
    ax.hist(
        vals,
        bins=bins,
        weights=weights,
        histtype="step",
        linewidth=1.15,
        label=label,
    )


def sanitize_figure_text(fig, strip_math: bool = False) -> None:
    """
    Sanitize every Matplotlib Text artist before rendering.

    Development-time plotting bugs have occasionally introduced doubled LaTeX
    backslashes such as ``\\otimes``.  Collapse those to a single backslash.

    If strip_math=True, remove dollar delimiters as a last-resort plain-text
    fallback.  This changes presentation only; physics arrays and fit results
    are untouched.
    """
    from matplotlib.text import Text

    for artist in fig.findobj(match=Text):
        value = artist.get_text()
        if not isinstance(value, str) or not value:
            continue
        #endif

        cleaned = value
        while "\\\\" in cleaned:
            cleaned = cleaned.replace("\\\\", "\\")
        #endwhile

        if strip_math:
            cleaned = cleaned.replace("$", "")
        #endif

        if cleaned != value:
            artist.set_text(cleaned)
        #endif
    #endfor


def safe_finalize_figure(
    fig,
    output_path: Path,
    rect: Tuple[float, float, float, float] = (0.0, 0.0, 1.0, 0.96),
) -> None:
    """
    Best-effort figure renderer.

    A plot-label or layout problem must never terminate a multi-period physics
    job after the numerical results have already been computed.

    Attempts:
      1. sanitize MathText, tight_layout, save;
      2. sanitize again, manual subplot spacing, save;
      3. strip MathText dollar delimiters and save as plain text;
      4. log and skip this one canvas.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sanitize_figure_text(fig, strip_math=False)

    try:
        fig.tight_layout(rect=rect)
    except Exception as exc:
        log(
            f"WARNING: tight_layout failed for {output_path.name}: {exc}. "
            "Using manual subplot spacing."
        )
        fig.subplots_adjust(
            left=0.08,
            right=0.97,
            bottom=0.10,
            top=max(0.50, rect[3] - 0.02),
            wspace=0.28,
            hspace=0.30,
        )
    #endtry

    try:
        fig.savefig(output_path, dpi=180)
        return
    except Exception as exc:
        log(
            f"WARNING: first render failed for {output_path.name}: {exc}. "
            "Retrying after text sanitization."
        )
    #endtry

    sanitize_figure_text(fig, strip_math=False)
    fig.subplots_adjust(
        left=0.08,
        right=0.97,
        bottom=0.10,
        top=max(0.50, rect[3] - 0.02),
        wspace=0.28,
        hspace=0.30,
    )

    try:
        fig.savefig(output_path, dpi=180)
        return
    except Exception as exc:
        log(
            f"WARNING: sanitized render failed for {output_path.name}: {exc}. "
            "Retrying with plain-text labels."
        )
    #endtry

    sanitize_figure_text(fig, strip_math=True)

    try:
        fig.savefig(output_path, dpi=180)
    except Exception as exc:
        log(
            f"WARNING: skipping unsavable canvas {output_path.name}: {exc}"
        )
    #endtry



def make_stage2_canvases(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    fit_rows: List[Dict[str, object]],
    spread_rows: List[Dict[str, object]],
    mixed_rows: List[Dict[str, object]],
    outdir: Path,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
) -> None:
    """
    Compact Stage-II diagnostic suite emphasizing discriminator stability.
    """
    common = {}
    for key, feat in (("data", data_f), ("pi0", pi0_f), ("dvcs", dvcs_f)):
        common[key] = (
            feat["valid_tag"]
            & (feat["pred_probe_energy"] >= 0.40)
            & (feat["pred_probe_energy"] < max_probe_energy)
            & (feat["ep_missing_mass2"] >= mm2_min)
            & (feat["ep_missing_mass2"] < mm2_max)
            & (np.abs(feat["pred_probe_mass2"]) < probe_m2_max)
        )
    #endfor

    # --------------------------------------------------------------
    # Canvas 1: discriminator shapes.
    # --------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 10.2))

    for key, feat, label in (
        ("data", data_f, "Data"),
        ("pi0", pi0_f, "aaogen $\\pi^0$"),
        ("dvcs", dvcs_f, "dvcsgen BH/DVCS"),
    ):
        unit_hist(
            axes[0, 0],
            feat["ep_missing_mass2"][common[key]],
            np.linspace(mm2_min, mm2_max, 220),
            label,
        )
        unit_hist(
            axes[0, 1],
            feat["pred_probe_mass2"][common[key]],
            np.linspace(-0.06, 0.06, 241),
            label,
        )
    #endfor

    axes[0, 0].axvline(M_PI0 * M_PI0, linestyle="--", linewidth=1.0)
    axes[0, 0].set_xlabel(r"$M_X^2(ep)$ (GeV$^2$)")
    axes[0, 0].set_ylabel("Unit-normalized entries")
    axes[0, 0].set_title("$ep$ missing mass squared")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.18)

    axes[0, 1].set_xlim(-0.06, 0.06)
    axes[0, 1].set_xlabel(r"$(P_{\mathrm{probe}}^{\mathrm{pred}})^2$ (GeV$^2$)")
    axes[0, 1].set_ylabel("Unit-normalized entries")
    axes[0, 1].set_title("Predicted probe mass squared")
    axes[0, 1].grid(alpha=0.18)

    if all("stored_pTmiss" in f for f in (data_f, pi0_f, dvcs_f)):
        vals_all = np.concatenate([
            f["stored_pTmiss"][common[k]]
            for k, f in (("data", data_f), ("pi0", pi0_f), ("dvcs", dvcs_f))
        ])
        vals_all = vals_all[np.isfinite(vals_all)]
        hi = max(0.2, min(1.5, float(np.percentile(vals_all, 99.5)))) if vals_all.size else 1.0
        for key, feat, label in (
            ("data", data_f, "Data"),
            ("pi0", pi0_f, "aaogen $\\pi^0$"),
            ("dvcs", dvcs_f, "dvcsgen BH/DVCS"),
        ):
            unit_hist(
                axes[1, 0],
                feat["stored_pTmiss"][common[key]],
                np.linspace(0.0, hi, 220),
                label,
            )
        #endfor
        axes[1, 0].set_xlabel(r"stored $p_{T,\mathrm{miss}}$ (GeV)")
        axes[1, 0].set_ylabel("Unit-normalized entries")
        axes[1, 0].set_title(r"Skim-level $p_{T,\mathrm{miss}}$ diagnostic")
    else:
        axes[1, 0].text(0.5, 0.5, "pTmiss unavailable", ha="center", va="center")
        axes[1, 0].set_axis_off()
    #endif
    axes[1, 0].grid(alpha=0.18)

    if all("stored_theta_gamma_gamma" in f for f in (data_f, pi0_f, dvcs_f)):
        vals_all = np.concatenate([
            f["stored_theta_gamma_gamma"][common[k]]
            for k, f in (("data", data_f), ("pi0", pi0_f), ("dvcs", dvcs_f))
        ])
        vals_all = vals_all[np.isfinite(vals_all)]
        theta_display_max = 4.0
        for key, feat, label in (
            ("data", data_f, "Data"),
            ("pi0", pi0_f, "aaogen $\\pi^0$"),
            ("dvcs", dvcs_f, "dvcsgen BH/DVCS"),
        ):
            unit_hist(
                axes[1, 1],
                feat["stored_theta_gamma_gamma"][common[key]],
                np.linspace(0.0, theta_display_max, 221),
                label,
            )
        #endfor
        axes[1, 1].set_xlim(0.0, 4.0)
        axes[1, 1].set_xlabel(r"stored $\theta_{\gamma\gamma}$ (deg)")
        axes[1, 1].set_ylabel("Unit-normalized entries")
        axes[1, 1].set_title(r"$\theta_{\gamma\gamma}$ diagnostic")
    else:
        axes[1, 1].text(
            0.5, 0.5, "theta_gamma_gamma unavailable",
            ha="center", va="center"
        )
        axes[1, 1].set_axis_off()
    #endif
    axes[1, 1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: Stage II discriminator shapes "
        f"(beam {period.beam_energy:.4f} GeV)",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_discriminator_shapes.png",
        rect=(0, 0, 1, 0.96),
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # Canvas 2: all-region discriminator comparison.
    # --------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 10.2))

    good_all = [
        r for r in fit_rows
        if r["region"] == "all"
        and int(r["fit_success"]) == 1
        and np.isfinite(r["pi0_fraction"])
    ]

    for disc in STAGE2_DISCRIMINATORS:
        rr = [r for r in good_all if r["discriminator"] == disc]
        if not rr:
            continue
        #endif
        axes[0, 0].plot(
            [r["energy_center_GeV"] for r in rr],
            [r["pi0_fraction"] for r in rr],
            marker="o", ms=3, label=disc,
        )
        axes[0, 1].plot(
            [r["energy_center_GeV"] for r in rr],
            [r["deviance_per_ndof"] for r in rr],
            marker="o", ms=3, label=disc,
        )
    #endfor

    axes[0, 0].set_ylim(0.0, 1.05)
    axes[0, 0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0, 0].set_ylabel("Fitted $\\pi^0$ fraction")
    axes[0, 0].set_title("Discriminator stability")
    axes[0, 0].legend(fontsize=7)
    axes[0, 0].grid(alpha=0.18)

    axes[0, 1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0, 1].set_ylabel("Poisson deviance / ndof")
    axes[0, 1].set_title("Absolute template goodness")
    axes[0, 1].legend(fontsize=7)
    axes[0, 1].grid(alpha=0.18)

    rrspread = [
        r for r in spread_rows
        if r["region"] == "all"
        and np.isfinite(r["full_spread_pi0_fraction"])
    ]
    if rrspread:
        x = [r["energy_center_GeV"] for r in rrspread]
        axes[1, 0].plot(
            x,
            [r["full_spread_pi0_fraction"] for r in rrspread],
            marker="o",
            label="full discriminator spread",
        )
        axes[1, 0].plot(
            x,
            [r["half_range_pi0_fraction"] for r in rrspread],
            marker="o",
            label="half range",
        )
        axes[1, 0].legend(fontsize=8)
    #endif
    axes[1, 0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1, 0].set_ylabel(r"spread in $f_{\pi^0}$")
    axes[1, 0].set_title("Model dependence of extracted composition")
    axes[1, 0].grid(alpha=0.18)

    # Mixed-event sensitivity for nominal discriminator only.
    shifts = sorted(set(int(r["mix_shift"]) for r in mixed_rows))
    for shift in shifts:
        rr = [
            r for r in mixed_rows
            if r["region"] == "all"
            and int(r["mix_shift"]) == shift
            and np.isfinite(r["delta_pi0_fraction"])
        ]
        if not rr:
            continue
        #endif
        axes[1, 1].plot(
            [r["energy_center_GeV"] for r in rr],
            [r["delta_pi0_fraction"] for r in rr],
            marker="o", ms=3, label=f"shift {shift}",
        )
    #endfor
    axes[1, 1].axhline(0.0, linestyle="--", linewidth=1.0)
    axes[1, 1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1, 1].set_ylabel(r"$\Delta f_{\pi^0}$")
    axes[1, 1].set_title("Diagnostic mixed-event sensitivity")
    axes[1, 1].legend(fontsize=7)
    axes[1, 1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: Stage II discriminator robustness",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_discriminator_stability.png",
        rect=(0, 0, 1, 0.96),
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # Canvas 3: requested FT / FD summary using the nominal discriminator,
    # with discriminator spread indicated as a band.
    # --------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5.6))

    def draw_region(ax, region: str, label: str) -> None:
        rr_nom = [
            r for r in fit_rows
            if r["region"] == region
            and r["discriminator"] == STAGE2_NOMINAL_DISCRIMINATOR
            and int(r["fit_success"]) == 1
            and np.isfinite(r["pi0_fraction"])
        ]
        rr_sp = {
            (float(r["energy_low_GeV"]), float(r["energy_high_GeV"])): r
            for r in spread_rows if r["region"] == region
        }

        if not rr_nom:
            return
        #endif

        x = np.asarray([r["energy_center_GeV"] for r in rr_nom])
        y = np.asarray([r["pi0_fraction"] for r in rr_nom])
        low = []
        high = []
        for r in rr_nom:
            key = (float(r["energy_low_GeV"]), float(r["energy_high_GeV"]))
            sr = rr_sp.get(key)
            low.append(
                float(sr["min_pi0_fraction"])
                if sr is not None and np.isfinite(sr["min_pi0_fraction"])
                else float(r["pi0_fraction"])
            )
            high.append(
                float(sr["max_pi0_fraction"])
                if sr is not None and np.isfinite(sr["max_pi0_fraction"])
                else float(r["pi0_fraction"])
            )
        #endfor

        ax.plot(x, y, marker="o", ms=3, label=label)
        ax.fill_between(x, low, high, alpha=0.18)
    #enddef

    draw_region(axes[0], "FT", "FT")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0].set_ylabel("Nominal fitted $\\pi^0$ fraction")
    axes[0].set_title("FT")
    axes[0].grid(alpha=0.18)

    for s in range(1, 7):
        draw_region(axes[1], f"FD_S{s}", f"S{s}")
    #endfor
    axes[1].set_ylim(0.0, 1.05)
    axes[1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1].set_ylabel("Nominal fitted $\\pi^0$ fraction")
    axes[1].set_title("FD sectors")
    axes[1].legend(ncol=2, fontsize=8)
    axes[1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: nominal denominator purity "
        "(bands = discriminator spread)",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_pi0_purity_by_region.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)

    # --------------------------------------------------------------
    # Canvas 4: representative nominal 2D fits.
    # --------------------------------------------------------------
    targets = ((0.40, 0.50), (0.80, 1.00), (1.50, 2.00))
    fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.2))

    for ax, (elo, ehi) in zip(axes, targets):
        candidates = [
            r for r in fit_rows
            if r["region"] == "all"
            and r["discriminator"] == STAGE2_NOMINAL_DISCRIMINATOR
            and int(r["fit_success"]) == 1
            and abs(float(r["energy_low_GeV"]) - elo) < 1.0e-9
            and abs(float(r["energy_high_GeV"]) - ehi) < 1.0e-9
        ]
        if not candidates:
            ax.text(0.5, 0.5, "No successful fit", ha="center", va="center")
            ax.set_axis_off()
            continue
        #endif

        row = candidates[0]
        centers, hd, hp, hv, model = representative_fit_projection(
            data_f, pi0_f, dvcs_f, row,
            ft_theta_max=ft_theta_max,
            mm2_min=mm2_min,
            mm2_max=mm2_max,
            probe_m2_max=probe_m2_max,
        )
        ax.errorbar(
            centers,
            hd,
            yerr=np.sqrt(np.maximum(hd, 1.0)),
            fmt="o",
            ms=2.5,
            label="Data",
        )
        ax.step(centers, model, where="mid", linewidth=1.4, label="Total fit")
        ax.step(centers, hp, where="mid", linewidth=1.0, label="$\\pi^0$")
        ax.step(centers, hv, where="mid", linewidth=1.0, label="BH/DVCS")
        ax.set_xlabel(r"$M_X^2(ep)$ (GeV$^2$)")
        ax.set_ylabel("Entries / bin")
        ax.set_title(
            f"{elo:.2f}–{ehi:.2f} GeV\n"
            f"$f_{{\\pi^0}}$={float(row['pi0_fraction']):.3f}, "
            f"D/ndof={float(row['deviance_per_ndof']):.2f}"
        )
        ax.grid(alpha=0.18)
    #endfor
    axes[0].legend(fontsize=8)

    fig.suptitle(
        f"{period.label}: representative nominal "
        f"{STAGE2_NOMINAL_DISCRIMINATOR} fits",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_representative_nominal_fits.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)



# -----------------------------------------------------------------------------
# Stage III: preliminary reference efficiency and data/MC scale factor
# -----------------------------------------------------------------------------

def provisional_ratio_error(num: float, den: float) -> float:
    """
    Counting-only uncertainty for num/den used strictly as a development
    diagnostic.  The fitted data denominator and shared-event correlations are
    not included; final uncertainties require event-level bootstrap replicas.
    """
    if not np.isfinite(num) or not np.isfinite(den) or den <= 0.0 or num < 0.0:
        return float("nan")
    #endif
    if num == 0.0:
        return 1.0 / den
    #endif
    return float(math.sqrt(num) / den)


def provisional_binomial_error(num: int, den: int) -> float:
    if den <= 0 or num < 0:
        return float("nan")
    #endif
    eff = float(num) / float(den)
    if not (0.0 <= eff <= 1.0):
        return float("nan")
    #endif
    return float(math.sqrt(max(eff * (1.0 - eff), 0.0) / den))


def build_stage3_reference_rows(
    period: PeriodConfig,
    shared_rows: List[Dict[str, object]],
    data_f: Dict[str, np.ndarray],
    aaogen_f: Dict[str, np.ndarray],
    data_stage_lookup: Dict[str, np.ndarray],
    data_mass_shell_success: np.ndarray,
    data_final_success: np.ndarray,
    mc_any_match: np.ndarray,
    mc_clean_success: np.ndarray,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
) -> List[Dict[str, object]]:
    """
    Construct the Stage-III reference extraction with a transparent data
    association cut-flow.

    The data numerator is reported in two forms:

      mass-shell-only:
          same-event eppi0 candidate whose tag removal leaves
          a massless, above-threshold companion;

      final:
          mass-shell-only plus the loose predicted/reconstructed companion
          consistency gate.

    The production-like reference efficiency uses ``data_final_success``.
    Keeping the mass-shell-only alternative makes the association sensitivity
    explicit and is especially useful for diagnosing epsilon_data > 1.
    """
    edges = stage2_energy_edges(max_probe_energy)
    regions = ["all", "FT", "FD_all"] + [f"FD_S{s}" for s in range(1, 7)]
    shared_map = {
        (
            str(r["region"]),
            round(float(r["energy_low_GeV"]), 9),
            round(float(r["energy_high_GeV"]), 9),
        ): r
        for r in shared_rows
    }

    required_stages = (
        "same_event",
        "positive_remainder",
        "tag_mass_shell",
        "probe_energy",
        "probe_pred_consistent",
    )
    for stage in required_stages:
        if stage not in data_stage_lookup:
            raise KeyError(f"Missing Stage-III data association stage '{stage}'.")
        #endif
    #endfor

    rows: List[Dict[str, object]] = []
    for region in regions:
        for ib in range(len(edges) - 1):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])
            key = (region, round(elo, 9), round(ehi, 9))
            fit_row = shared_map.get(key)

            md = stage2_fit_mask(
                data_f, region, ft_theta_max, elo, ehi,
                mm2_min, mm2_max, probe_m2_max,
            )
            mm = stage2_fit_mask(
                aaogen_f, region, ft_theta_max, elo, ehi,
                mm2_min, mm2_max, probe_m2_max,
            )

            n_data_tags = int(np.count_nonzero(md))
            stage_counts = {
                stage: int(np.count_nonzero(md & data_stage_lookup[stage]))
                for stage in required_stages
            }
            n_data_mass_shell = int(
                np.count_nonzero(md & data_mass_shell_success)
            )
            n_data_success = int(
                np.count_nonzero(md & data_final_success)
            )

            n_mc_tags = int(np.count_nonzero(mm))
            n_mc_any = int(np.count_nonzero(mm & mc_any_match))
            n_mc_success = int(np.count_nonzero(mm & mc_clean_success))

            fpi0 = float("nan")
            fit_success = 0
            fit_dev = float("nan")
            stage2_data_count = -1
            stage2_aaogen_count = -1
            if fit_row is not None:
                fit_success = int(fit_row.get("fit_success", 0))
                fpi0 = float(fit_row.get("pi0_fraction", float("nan")))
                fit_dev = float(fit_row.get("deviance_per_ndof", float("nan")))
                stage2_data_count = int(
                    fit_row.get("data_candidate_count", -1)
                )
                stage2_aaogen_count = int(
                    fit_row.get("aaogen_candidate_count", -1)
                )
            #endif

            data_den = (
                fpi0 * n_data_tags
                if fit_success and np.isfinite(fpi0) and fpi0 > 0.0
                else float("nan")
            )

            data_eff_mass_shell = (
                float(n_data_mass_shell) / data_den
                if np.isfinite(data_den) and data_den > 0.0
                else float("nan")
            )
            data_eff = (
                float(n_data_success) / data_den
                if np.isfinite(data_den) and data_den > 0.0
                else float("nan")
            )
            mc_eff = (
                float(n_mc_success) / float(n_mc_tags)
                if n_mc_tags > 0 else float("nan")
            )
            sf = (
                data_eff / mc_eff
                if np.isfinite(data_eff)
                and np.isfinite(mc_eff)
                and mc_eff > 0.0
                else float("nan")
            )
            sf_mass_shell = (
                data_eff_mass_shell / mc_eff
                if np.isfinite(data_eff_mass_shell)
                and np.isfinite(mc_eff)
                and mc_eff > 0.0
                else float("nan")
            )

            data_err = provisional_ratio_error(
                float(n_data_success), data_den
            )
            mc_err = provisional_binomial_error(n_mc_success, n_mc_tags)
            if (
                np.isfinite(sf)
                and sf != 0.0
                and np.isfinite(data_err)
                and np.isfinite(mc_err)
                and np.isfinite(data_eff)
                and data_eff != 0.0
                and np.isfinite(mc_eff)
                and mc_eff != 0.0
            ):
                sf_err = abs(sf) * math.sqrt(
                    (data_err / data_eff) ** 2
                    + (mc_err / mc_eff) ** 2
                )
            else:
                sf_err = float("nan")
            #endif

            status = "ok"
            if fit_row is None or fit_success != 1:
                status = "no_successful_stage2_fit"
            elif (
                (stage2_data_count >= 0 and stage2_data_count != n_data_tags)
                or (
                    stage2_aaogen_count >= 0
                    and stage2_aaogen_count != n_mc_tags
                )
            ):
                status = "stage2_stage3_candidate_support_mismatch"
            elif n_mc_tags <= 0:
                status = "no_aaogen_denominator"
            elif stage_counts["same_event"] <= 0:
                status = "no_same_event_overlap"
            # Exact runnum+evnum defines data parent identity in v026.
            # Parent kinematic differences are diagnostic only.
            elif stage_counts["tag_mass_shell"] <= 0:
                status = "no_kinematic_tag_identity"
            elif n_data_mass_shell <= 0:
                status = "no_above_threshold_mass_shell_probe"
            elif n_data_success <= 0:
                status = "mass_shell_probe_but_no_pred_consistent_probe"
            elif np.isfinite(data_eff) and data_eff > 1.05:
                status = (
                    "data_efficiency_above_unity_check_association_or_denominator"
                )
            elif np.isfinite(mc_eff) and mc_eff > 1.0 + 1.0e-12:
                status = "mc_efficiency_above_unity_internal_error"
            #endif

            row = {
                "period": period.key,
                "label": period.label,
                "beam_energy_GeV": period.beam_energy,
                "data_mode": "wagon_reference",
                "region": region,
                "energy_low_GeV": elo,
                "energy_high_GeV": ehi,
                "energy_center_GeV": 0.5 * (elo + ehi),
                "stage2_fit_success": fit_success,
                "stage2_deviance_per_ndof": fit_dev,
                "fitted_pi0_fraction": fpi0,
                "stage2_data_candidate_count": stage2_data_count,
                "stage2_aaogen_candidate_count": stage2_aaogen_count,
                "stage2_minus_stage3_data_candidate_count": (
                    stage2_data_count - n_data_tags
                    if stage2_data_count >= 0 else -999999
                ),
                "stage2_minus_stage3_aaogen_candidate_count": (
                    stage2_aaogen_count - n_mc_tags
                    if stage2_aaogen_count >= 0 else -999999
                ),
                "data_tag_candidates": n_data_tags,
                "data_fitted_pi0_denominator": data_den,
                "data_mass_shell_reconstructed_probe": n_data_mass_shell,
                "data_clean_reconstructed_probe": n_data_success,
                "data_efficiency_mass_shell_only": data_eff_mass_shell,
                "data_efficiency": data_eff,
                "data_efficiency_counting_error_provisional": data_err,
                "mc_tag_candidates": n_mc_tags,
                "mc_any_eppi0_parent_match": n_mc_any,
                "mc_clean_reconstructed_probe": n_mc_success,
                "mc_efficiency": mc_eff,
                "mc_efficiency_binomial_error_provisional": mc_err,
                "data_over_mc_scale_factor_mass_shell_only": sf_mass_shell,
                "data_over_mc_scale_factor": sf,
                "scale_factor_counting_error_provisional": sf_err,
                "data_mass_shell_success_fraction_of_tags_raw": (
                    n_data_mass_shell / n_data_tags
                    if n_data_tags else float("nan")
                ),
                "data_clean_success_fraction_of_tags_raw": (
                    n_data_success / n_data_tags
                    if n_data_tags else float("nan")
                ),
                "mc_any_match_fraction_of_tags": (
                    n_mc_any / n_mc_tags if n_mc_tags else float("nan")
                ),
                "uncertainty_model": (
                    "provisional counting-only; fitted-denominator and "
                    "shared-event correlations omitted; final result requires "
                    "event bootstrap"
                ),
                "status": status,
            }

            # Full sequential cut-flow, both counts and raw-tag fractions.
            previous = n_data_tags
            for stage in required_stages:
                count = stage_counts[stage]
                row[f"data_cutflow_{stage}"] = count
                row[f"data_cutflow_{stage}_fraction_of_raw_tags"] = (
                    count / n_data_tags if n_data_tags else float("nan")
                )
                row[f"data_cutflow_{stage}_conditional_fraction"] = (
                    count / previous if previous else float("nan")
                )
                previous = count
            #endfor

            rows.append(row)
        #endfor
    #endfor
    return rows


def _stage3_rows_for_region(
    rows: List[Dict[str, object]],
    region: str,
) -> List[Dict[str, object]]:
    return sorted(
        [
            r for r in rows
            if r["region"] == region
            and np.isfinite(r.get("energy_center_GeV", np.nan))
        ],
        key=lambda r: float(r["energy_center_GeV"]),
    )



def _qfinite(values: np.ndarray, q: float) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.quantile(arr, q)) if arr.size else float("nan")



def make_stage3_canvases(
    period: PeriodConfig,
    rows: List[Dict[str, object]],
    outdir: Path,
    compact: bool = True,
) -> None:
    """
    Compact Stage-III production-development output.

    Only three canvases are retained:
      1. data and MC efficiencies;
      2. data/MC scale factors;
      3. compact numerator/denominator bookkeeping.

    The detailed parent-delta, multiple-electron, cut-flow, conditional-
    retention, and residual canvases were useful during development but are
    intentionally removed from the normal production output.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Needed by the scale-factor canvas in BOTH compact and full modes.
    ft = _stage3_rows_for_region(rows, "FT")

    if not compact:
        # ------------------------------------------------------------------
        # 1. Efficiencies: FT and FD sectors.
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.8))
        if ft:
            x = np.asarray([r["energy_center_GeV"] for r in ft], dtype=float)
            axes[0].errorbar(
                x,
                [r["data_efficiency"] for r in ft],
                yerr=[
                    r["data_efficiency_counting_error_provisional"]
                    for r in ft
                ],
                marker="o",
                linestyle="-",
                label="Data",
            )
            axes[0].errorbar(
                x,
                [r["mc_efficiency"] for r in ft],
                yerr=[
                    r["mc_efficiency_binomial_error_provisional"]
                    for r in ft
                ],
                marker="o",
                linestyle="--",
                label="aaogen MC",
            )
            # Retain mass-shell-only as a small diagnostic overlay, not a separate
            # canvas.
            axes[0].plot(
                x,
                [r["data_efficiency_mass_shell_only"] for r in ft],
                marker=".",
                linestyle=":",
                linewidth=0.9,
                label="Data, mass-shell only",
            )
        #endif

        finite_eff = [
            float(r[key])
            for r in rows
            for key in (
                "data_efficiency",
                "data_efficiency_mass_shell_only",
                "mc_efficiency",
            )
            if np.isfinite(r.get(key, np.nan))
        ]
        eff_ymax = (
            max(1.15, min(2.0, 1.10 * max(finite_eff)))
            if finite_eff else 1.15
        )
        axes[0].set_ylim(0.0, eff_ymax)
        axes[0].axhline(1.0, linestyle="--", linewidth=0.8)
        axes[0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        axes[0].set_ylabel("Photon reconstruction efficiency")
        axes[0].set_title("FT")
        axes[0].legend(fontsize=8)
        axes[0].grid(alpha=0.18)

        for s in range(1, 7):
            rr = _stage3_rows_for_region(rows, f"FD_S{s}")
            if not rr:
                continue
            #endif
            x = [r["energy_center_GeV"] for r in rr]
            axes[1].plot(
                x,
                [r["data_efficiency"] for r in rr],
                marker="o",
                linewidth=1.0,
                label=f"Data S{s}",
            )
            axes[1].plot(
                x,
                [r["mc_efficiency"] for r in rr],
                marker=".",
                linestyle="--",
                linewidth=0.8,
                label=f"MC S{s}",
            )
        #endfor
        axes[1].axhline(1.0, linestyle="--", linewidth=0.8)
        axes[1].set_ylim(0.0, eff_ymax)
        axes[1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        axes[1].set_ylabel("Photon reconstruction efficiency")
        axes[1].set_title("FD sectors")
        axes[1].legend(fontsize=7, ncol=3)
        axes[1].grid(alpha=0.18)

        fig.suptitle(
            f"{period.label}: Stage-III reference photon efficiencies",
            fontsize=13,
        )
        safe_finalize_figure(
            fig,
            outdir / "canvas_reference_efficiencies.png",
            rect=(0, 0, 1, 0.93),
        )
        plt.close(fig)

    #endif

    # ------------------------------------------------------------------
    # 2. Data/MC scale factors.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.8))
    if ft:
        x = [r["energy_center_GeV"] for r in ft]
        axes[0].errorbar(
            x,
            [r["data_over_mc_scale_factor"] for r in ft],
            yerr=[
                r["scale_factor_counting_error_provisional"]
                for r in ft
            ],
            marker="o",
            linestyle="-",
            label="Final",
        )
        axes[0].plot(
            x,
            [r["data_over_mc_scale_factor_mass_shell_only"] for r in ft],
            marker=".",
            linestyle=":",
            linewidth=0.9,
            label="Mass-shell only",
        )
    #endif
    axes[0].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0].set_ylabel(r"$\epsilon_{\mathrm{data}}/\epsilon_{\mathrm{MC}}$")
    axes[0].set_title("FT")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.18)

    for s in range(1, 7):
        rr = _stage3_rows_for_region(rows, f"FD_S{s}")
        if rr:
            axes[1].plot(
                [r["energy_center_GeV"] for r in rr],
                [r["data_over_mc_scale_factor"] for r in rr],
                marker="o",
                linewidth=1.0,
                label=f"S{s}",
            )
        #endif
    #endfor
    axes[1].axhline(1.0, linestyle="--", linewidth=1.0)
    axes[1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1].set_ylabel(r"$\epsilon_{\mathrm{data}}/\epsilon_{\mathrm{MC}}$")
    axes[1].set_title("FD sectors")
    axes[1].legend(fontsize=8, ncol=2)
    axes[1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: preliminary photon-efficiency data/MC scale factor",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_reference_scale_factors.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)

    if not compact:
        # ------------------------------------------------------------------
        # 3. Minimal bookkeeping: fitted pi0 denominator and final numerator.
        # ------------------------------------------------------------------
        fig, axes = plt.subplots(1, 2, figsize=(15.0, 5.8))
        for ax, region, title in (
            (axes[0], "FT", "FT"),
            (axes[1], "FD_all", "FD all"),
        ):
            rr = _stage3_rows_for_region(rows, region)
            if rr:
                x = [r["energy_center_GeV"] for r in rr]
                ax.plot(
                    x,
                    [r["data_fitted_pi0_denominator"] for r in rr],
                    marker="o",
                    label=r"Fitted data $\pi^0$ denominator",
                )
                ax.plot(
                    x,
                    [r["data_clean_reconstructed_probe"] for r in rr],
                    marker="o",
                    label="Final reconstructed-probe numerator",
                )
                ax.plot(
                    x,
                    [r["mc_clean_reconstructed_probe"] for r in rr],
                    marker=".",
                    linestyle="--",
                    label="aaogen reconstructed probe",
                )
            #endif
            ax.set_yscale("log")
            ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
            ax.set_ylabel("Directed tag-probe candidates")
            ax.set_title(title)
            ax.legend(fontsize=8)
            ax.grid(alpha=0.18)
        #endfor

        fig.suptitle(
            f"{period.label}: Stage-III numerator/denominator bookkeeping",
            fontsize=13,
        )
        safe_finalize_figure(
            fig,
            outdir / "canvas_reference_bookkeeping.png",
            rect=(0, 0, 1, 0.93),
        )
        plt.close(fig)


    #endif
def summarize_stage3(
    period: PeriodConfig,
    rows: List[Dict[str, object]],
    data_match_counters: Dict[str, int],
    max_probe_energy: float,
    use_pred_consistency: bool,
    association_settings: Dict[str, float],
) -> Dict[str, object]:
    good = [
        r for r in rows
        if r["region"] in ("FT", "FD_all")
        and r.get("status") == "ok"
        and np.isfinite(r.get("data_over_mc_scale_factor", np.nan))
    ]
    status_counts: Dict[str, int] = {}
    for row in rows:
        status = str(row.get("status", "unknown"))
        status_counts[status] = status_counts.get(status, 0) + 1
    #endfor

    return {
        "period": period.key,
        "label": period.label,
        "beam_energy_GeV": period.beam_energy,
        "data_mode": "wagon_reference",
        "probe_energy_min_GeV": 0.40,
        "probe_energy_max_GeV": max_probe_energy,
        "purpose": (
            "Stage-III association-validation reference extraction for overlap "
            "testing; future loose-nSidis production data will replace/extend "
            "the data side"
        ),
        "parent_consistency_diagnostic": "Active 0.002-GeV-style Cartesian gate retained for this run; independent tolerance scan quantifies its appropriateness.",
        "association_schema_limit": (
            "Current eppi0 trees do not store daughter-photon full three-vectors. "
            "Tag identity is therefore tested by exact event identity, e/p "
            "compatibility, and the massless remainder P_pi0-k_tag."
        ),
        "data_parent_identity": "exact runnum+evnum; no data e/p kinematic matching",
        "use_predicted_probe_consistency_in_numerator": bool(
            use_pred_consistency
        ),
        "association_settings": association_settings,
        "uncertainty_status": (
            "provisional counting-only diagnostics; final statistical covariance "
            "requires event-level bootstrap"
        ),
        "data_exact_event_matching": data_match_counters,
        "status_counts": status_counts,
        "n_reference_rows": len(rows),
        "n_good_ft_fdall_rows": len(good),
        "median_good_scale_factor": (
            float(
                np.median(
                    [r["data_over_mc_scale_factor"] for r in good]
                )
            )
            if good else float("nan")
        ),
    }


def summarize_stage2(
    period: PeriodConfig,
    fit_rows: List[Dict[str, object]],
    spread_rows: List[Dict[str, object]],
    mixed_rows: List[Dict[str, object]],
    data_total: int,
    pi0_total: int,
    dvcs_total: int,
    max_probe_energy: float,
) -> Dict[str, object]:
    nominal_good = [
        r for r in fit_rows
        if r["region"] == "all"
        and r["discriminator"] == STAGE2_NOMINAL_DISCRIMINATOR
        and int(r["fit_success"]) == 1
        and np.isfinite(r["pi0_yield"])
    ]

    pi0_yield = (
        float(np.sum([r["pi0_yield"] for r in nominal_good]))
        if nominal_good else float("nan")
    )
    dvcs_yield = (
        float(np.sum([r["dvcs_yield"] for r in nominal_good]))
        if nominal_good else float("nan")
    )
    total = pi0_yield + dvcs_yield if nominal_good else float("nan")

    all_spreads = [
        float(r["full_spread_pi0_fraction"])
        for r in spread_rows
        if r["region"] == "all"
        and np.isfinite(r["full_spread_pi0_fraction"])
    ]

    mix_delta = [
        abs(float(r["delta_pi0_fraction"]))
        for r in mixed_rows
        if r["region"] == "all"
        and np.isfinite(r["delta_pi0_fraction"])
    ]

    discriminator_summary: Dict[str, Dict[str, float]] = {}
    for disc in STAGE2_DISCRIMINATORS:
        rr = [
            r for r in fit_rows
            if r["region"] == "all"
            and r["discriminator"] == disc
            and int(r["fit_success"]) == 1
        ]
        if rr:
            discriminator_summary[disc] = {
                "successful_bins": int(len(rr)),
                "median_pi0_fraction": float(
                    np.median([r["pi0_fraction"] for r in rr])
                ),
                "median_deviance_per_ndof": float(
                    np.median([r["deviance_per_ndof"] for r in rr])
                ),
            }
        #endif
    #endfor

    return {
        "period": period.key,
        "label": period.label,
        "beam_energy_GeV": period.beam_energy,
        "stage2_probe_energy_max_GeV": max_probe_energy,
        "nominal_discriminator": STAGE2_NOMINAL_DISCRIMINATOR,
        "epgamma_data_entries_total": int(data_total),
        "aaogen_epgamma_entries_total": int(pi0_total),
        "dvcsgen_epgamma_entries_total": int(dvcs_total),
        "successful_nominal_all_region_energy_fits": int(len(nominal_good)),
        "summed_nominal_pi0_yield_in_fit_windows": pi0_yield,
        "summed_nominal_dvcs_yield_in_fit_windows": dvcs_yield,
        "summed_nominal_pi0_fraction": (
            pi0_yield / total
            if nominal_good and total > 0.0 else float("nan")
        ),
        "summed_nominal_dvcs_fraction": (
            dvcs_yield / total
            if nominal_good and total > 0.0 else float("nan")
        ),
        "median_nominal_deviance_per_ndof": (
            float(np.median([r["deviance_per_ndof"] for r in nominal_good]))
            if nominal_good else float("nan")
        ),
        "median_all_region_discriminator_spread": (
            float(np.median(all_spreads)) if all_spreads else float("nan")
        ),
        "max_all_region_discriminator_spread": (
            max(all_spreads) if all_spreads else float("nan")
        ),
        "median_abs_pi0_fraction_shift_from_mixed_component": (
            float(np.median(mix_delta)) if mix_delta else float("nan")
        ),
        "max_abs_pi0_fraction_shift_from_mixed_component": (
            max(mix_delta) if mix_delta else float("nan")
        ),
        "discriminator_summary": discriminator_summary,
    }


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def save_hist(
    values: np.ndarray,
    path: Path,
    xlabel: str,
    title: str,
    bins: int | np.ndarray = 120,
    xlim: Optional[Tuple[float, float]] = None,
    logy: bool = False,
) -> None:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    ax.hist(values, bins=bins, histtype="step", linewidth=1.4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    ax.set_title(title)
    if xlim is not None:
        ax.set_xlim(*xlim)
    #endif
    if logy:
        ax.set_yscale("log")
    #endif
    ax.grid(alpha=0.20)
    safe_finalize_figure(fig, Path(path), rect=(0, 0, 1, 1))
    plt.close(fig)


def save_hist2d(
    x: np.ndarray,
    y: np.ndarray,
    path: Path,
    xlabel: str,
    ylabel: str,
    title: str,
    bins: Tuple[int, int] = (120, 120),
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    logz: bool = False,
) -> None:
    from matplotlib.colors import LogNorm

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    norm = LogNorm() if logz and x.size else None
    h = ax.hist2d(x, y, bins=bins, norm=norm)
    fig.colorbar(h[3], ax=ax, label="Entries")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if xlim is not None:
        ax.set_xlim(*xlim)
    #endif
    if ylim is not None:
        ax.set_ylim(*ylim)
    #endif
    ax.grid(alpha=0.12)
    safe_finalize_figure(fig, Path(path), rect=(0, 0, 1, 1))
    plt.close(fig)


def finite_percentile(values: np.ndarray, q: float, default: float) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return default
    #endif
    return float(np.percentile(values, q))


def symmetric_zoom(values: np.ndarray, q: float, floor: float, ceiling: Optional[float] = None) -> float:
    values = np.asarray(values, dtype=float)
    values = np.abs(values[np.isfinite(values)])
    if values.size == 0:
        lim = floor
    else:
        lim = max(floor, float(np.percentile(values, q)))
    #endif
    if ceiling is not None:
        lim = min(lim, ceiling)
    #endif
    return lim


def make_full_and_zoom_hist(
    values: np.ndarray,
    outdir: Path,
    stem: str,
    xlabel: str,
    title: str,
    zoom_xlim: Tuple[float, float],
    bins: int = 160,
    logy_full: bool = True,
) -> None:
    save_hist(
        values,
        outdir / f"{stem}_full.png",
        xlabel,
        title + " (full range)",
        bins=bins,
        xlim=None,
        logy=logy_full,
    )
    save_hist(
        values,
        outdir / f"{stem}_zoom.png",
        xlabel,
        title + " (core zoom)",
        bins=bins,
        xlim=zoom_xlim,
        logy=False,
    )



def _finite(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return values[np.isfinite(values)]


def _hist_panel(
    ax,
    values: np.ndarray,
    xlabel: str,
    bins: int = 260,
    xlim: Optional[Tuple[float, float]] = None,
    logy: bool = False,
) -> None:
    values = _finite(values)
    ax.hist(values, bins=bins, histtype="step", linewidth=1.3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Entries")
    if xlim is not None:
        ax.set_xlim(*xlim)
    #endif
    if logy:
        ax.set_yscale("log")
    #endif
    ax.grid(alpha=0.18)


def _hist2d_panel(
    ax,
    x: np.ndarray,
    y: np.ndarray,
    xlabel: str,
    ylabel: str,
    bins: Tuple[int, int] = (160, 160),
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    logz: bool = True,
):
    from matplotlib.colors import LogNorm
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    norm = LogNorm() if logz and np.any(mask) else None
    h = ax.hist2d(x[mask], y[mask], bins=bins, norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    #endif
    if ylim is not None:
        ax.set_ylim(*ylim)
    #endif
    ax.grid(alpha=0.10)
    return h[3]


def make_compact_canvases(
    period: PeriodConfig,
    arrays: Dict[str, np.ndarray],
    clean: np.ndarray,
    rows: List[Dict[str, object]],
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_mass2_max: float,
) -> None:
    title = period.label
    outdir = Path(".")  # overwritten by caller through matplotlib path helper


def make_plots(
    period: PeriodConfig,
    arrays: Dict[str, np.ndarray],
    clean: np.ndarray,
    rows: List[Dict[str, object]],
    outdir: Path,
    assoc_mgg_min: float,
    assoc_mgg_max: float,
    assoc_mass2_max: float,
) -> None:
    if len(arrays["pred_probe_energy"]) == 0:
        log(f"{period.label}: no accepted pairs; skipping plots.")
        return
    #endif

    for old_png in outdir.glob("*.png"):
        old_png.unlink()
    #endfor

    title = period.label
    all_mask = np.ones(len(clean), dtype=bool)

    # ------------------------------------------------------------------
    # Canvas 1: define and validate the reconstructed-side association.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))

    vals = arrays["pi0_mass"]
    axes[0, 0].hist(vals[np.isfinite(vals)], bins=350, histtype="step", linewidth=1.2)
    axes[0, 0].axvline(assoc_mgg_min, linestyle="--")
    axes[0, 0].axvline(assoc_mgg_max, linestyle="--")
    axes[0, 0].set_xlim(0.08, 0.20)
    axes[0, 0].set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
    axes[0, 0].set_ylabel("Entries")
    axes[0, 0].set_title("Reconstructed $\\pi^0$ mass window")
    axes[0, 0].grid(alpha=0.18)

    m2 = arrays["reco_probe_mass2"]
    finite = np.isfinite(m2)
    axes[0, 1].hist(m2[finite], bins=500, histtype="step", linewidth=1.2)
    axes[0, 1].axvline(-assoc_mass2_max, linestyle="--")
    axes[0, 1].axvline(+assoc_mass2_max, linestyle="--")
    axes[0, 1].set_xlim(-0.02, 0.02)
    axes[0, 1].set_yscale("log")
    axes[0, 1].set_xlabel(r"$(P_{\pi^0}^{reco}-k_{tag}^{reco})^2$ (GeV$^2$)")
    axes[0, 1].set_ylabel("Entries")
    axes[0, 1].set_title("Companion-photon mass shell")
    axes[0, 1].grid(alpha=0.18)

    bins_e = np.linspace(0.0, 7.0, 141)
    axes[1, 0].hist(
        arrays["reco_probe_energy"][all_mask],
        bins=bins_e, histtype="step", linewidth=1.1, label="All parent matches"
    )
    axes[1, 0].hist(
        arrays["reco_probe_energy"][clean],
        bins=bins_e, histtype="step", linewidth=1.3, label="Clean association"
    )
    axes[1, 0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{reco}}$ (GeV)")
    axes[1, 0].set_ylabel("Entries")
    axes[1, 0].set_title("Companion-photon energy")
    axes[1, 0].legend()
    axes[1, 0].grid(alpha=0.18)

    bins_a = np.linspace(0.0, 180.0, 361)
    axes[1, 1].hist(
        arrays["probe_opening_residual_deg"][all_mask],
        bins=bins_a, histtype="step", linewidth=1.0, label="All parent matches"
    )
    axes[1, 1].hist(
        arrays["probe_opening_residual_deg"][clean],
        bins=bins_a, histtype="step", linewidth=1.3, label="Clean association"
    )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlim(0.0, 180.0)
    axes[1, 1].set_xlabel(r"$\Delta\alpha_{\rm probe}$ (deg)")
    axes[1, 1].set_ylabel("Entries")
    axes[1, 1].set_title("Prediction residual: diagnostic only")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.18)

    nall = len(clean)
    nclean = int(np.count_nonzero(clean))
    fig.suptitle(
        f"{title}: clean reconstructed $\\pi^0$ tag association "
        f"({nclean:,}/{nall:,} = {100.0*nclean/max(nall,1):.1f}%)",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_clean_tag_association.png",
        rect=(0, 0, 1, 0.96),
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # Canvas 2: resolution in the clean sample.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))
    from matplotlib.colors import LogNorm

    x = arrays["pred_probe_energy"][clean]
    y = arrays["reco_probe_energy"][clean]
    finite = np.isfinite(x) & np.isfinite(y)
    h = axes[0, 0].hist2d(
        x[finite], y[finite], bins=(180, 180),
        range=((0.0, 7.0), (0.0, 7.0)), norm=LogNorm()
    )
    fig.colorbar(h[3], ax=axes[0, 0], label="Entries")
    axes[0, 0].plot([0, 7], [0, 7], linestyle="--", linewidth=1.0)
    axes[0, 0].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[0, 0].set_ylabel(r"$E_{\mathrm{probe}}^{\mathrm{reco}}$ (GeV)")
    axes[0, 0].set_title("Predicted vs reconstructed energy")

    core_a = arrays["probe_opening_residual_deg"][clean]
    core_a = core_a[np.isfinite(core_a)]
    axes[0, 1].hist(core_a, bins=600, range=(0.0, 20.0), histtype="step", linewidth=1.25)
    axes[0, 1].set_xlabel(r"$\Delta\alpha_{\rm probe}$ (deg)")
    axes[0, 1].set_ylabel("Entries")
    axes[0, 1].set_title("Angular residual (fine core binning)")
    axes[0, 1].grid(alpha=0.18)

    de = arrays["probe_delta_E"][clean]
    de = de[np.isfinite(de)]
    axes[1, 0].hist(de, bins=500, range=(-1.0, 1.0), histtype="step", linewidth=1.25)
    axes[1, 0].set_xlabel(r"$E_{\mathrm{reco}}-E_{\mathrm{pred}}$ (GeV)")
    axes[1, 0].set_ylabel("Entries")
    axes[1, 0].set_title("Energy residual (fine core binning)")
    axes[1, 0].grid(alpha=0.18)

    E = arrays["pred_probe_energy"][clean]
    A = arrays["probe_opening_residual_deg"][clean]
    finite = np.isfinite(E) & np.isfinite(A)
    h = axes[1, 1].hist2d(
        E[finite], A[finite], bins=(160, 160),
        range=((0.4, 7.0), (0.0, 30.0)), norm=LogNorm()
    )
    fig.colorbar(h[3], ax=axes[1, 1], label="Entries")
    axes[1, 1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1, 1].set_ylabel(r"$\Delta\alpha_{\rm probe}$ (deg)")
    axes[1, 1].set_title("Angular resolution vs predicted energy")

    fig.suptitle(f"{title}: clean-sample probe prediction resolution", fontsize=15)
    safe_finalize_figure(
        fig,
        outdir / "canvas_clean_probe_resolution.png",
        rect=(0, 0, 1, 0.96),
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # Canvas 3: energy-dependent quantiles, all + detector-region split.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))
    region_order = ["all", "FT_like"] + [f"FD_S{s}" for s in range(1, 7)]

    for region in region_order:
        rr = [r for r in rows if r["region"] == region and np.isfinite(r["angle_q68_deg"])]
        if not rr:
            continue
        #endif
        xx = np.asarray([r["energy_center_GeV"] for r in rr])
        axes[0, 0].plot(xx, [r["angle_q68_deg"] for r in rr], marker="o", ms=3, label=region)
        axes[0, 1].plot(xx, [r["angle_q95_deg"] for r in rr], marker="o", ms=3, label=region)
        axes[1, 0].plot(xx, [r["angle_q99_deg"] for r in rr], marker="o", ms=3, label=region)
    #endfor

    for ax, q in (
        (axes[0, 0], "68%"),
        (axes[0, 1], "95%"),
        (axes[1, 0], "99%"),
    ):
        ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        ax.set_ylabel(r"$\Delta\alpha$ quantile (deg)")
        ax.set_title(f"{q} containment")
        ax.grid(alpha=0.18)
    #endfor
    axes[0, 1].legend(fontsize=8, ncol=2)

    # Counts by energy for all clean probes.
    rr = [r for r in rows if r["region"] == "all"]
    axes[1, 1].step(
        [r["energy_center_GeV"] for r in rr],
        [r["count"] for r in rr],
        where="mid",
    )
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    axes[1, 1].set_ylabel("Clean entries")
    axes[1, 1].set_title("Statistics by predicted energy")
    axes[1, 1].grid(alpha=0.18)

    fig.suptitle(
        f"{title}: candidate angular matching containment by predicted region",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_matching_resolution_quantiles.png",
        rect=(0, 0, 1, 0.96),
    )
    plt.close(fig)

    # ------------------------------------------------------------------
    # Canvas 4: check that clean association did not impose exclusivity.
    # ------------------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))
    variables = (
        ("exclusivity_missing_energy", r"$E_{\rm miss}(ep\pi^0)$ (GeV)", (-1.5, 1.5)),
        ("exclusivity_missing_p", r"$|\vec p_{\rm miss}(ep\pi^0)|$ (GeV)", (0.0, 1.5)),
        ("exclusivity_missing_pT", r"$p_{T,\rm miss}(ep\pi^0)$ (GeV)", (0.0, 0.75)),
        ("pred_probe_mass2", r"$(p_{\rm probe}^{pred})^2$ (GeV$^2$)", (-1.0, 1.0)),
    )
    for ax, (name, xlabel, lim) in zip(axes.flat, variables):
        vals = arrays[name][clean]
        vals = vals[np.isfinite(vals)]
        ax.hist(vals, bins=350, range=lim, histtype="step", linewidth=1.2)
        ax.set_xlim(*lim)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Entries")
        ax.grid(alpha=0.18)
    #endfor
    fig.suptitle(
        f"{title}: exclusivity diagnostics after reconstructed-side association",
        fontsize=15,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_exclusivity_after_clean_association.png",
        rect=(0, 0, 1, 0.96),
    )
    plt.close(fig)
# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

def percentile_or_nan(x: np.ndarray, q: float) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.percentile(x, q)) if x.size else float("nan")


def rms_or_nan(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(np.sqrt(np.mean(x * x))) if x.size else float("nan")


def summarize_period(
    period: PeriodConfig,
    epg_entries_total: int,
    eppi0_entries_total: int,
    epg_loaded: int,
    eppi0_loaded: int,
    epg_angle_unit: str,
    eppi0_angle_unit: str,
    matches: MatchResult,
    arrays: Dict[str, np.ndarray],
    counters: Dict[str, int],
    clean: Optional[np.ndarray] = None,
) -> Dict[str, object]:
    n = len(arrays["pred_probe_energy"])
    summary: Dict[str, object] = {
        "period": period.key,
        "label": period.label,
        "beam_energy_GeV": period.beam_energy,
        "epgamma_entries_total": epg_entries_total,
        "eppi0_entries_total": eppi0_entries_total,
        "epgamma_entries_loaded": epg_loaded,
        "eppi0_entries_loaded": eppi0_loaded,
        "epgamma_angle_unit": epg_angle_unit,
        "eppi0_angle_unit": eppi0_angle_unit,
        **counters,
        "stage1_pair_fraction_of_parent_matches": (
            n / len(matches.epg_index) if len(matches.epg_index) else float("nan")
        ),
        "median_parent_max_component_delta_GeV": percentile_or_nan(
            arrays["parent_max_component_delta"], 50.0
        ),
        "p99_parent_max_component_delta_GeV": percentile_or_nan(
            arrays["parent_max_component_delta"], 99.0
        ),
        "median_Mh_gammagamma_GeV": percentile_or_nan(
            arrays["pi0_mass"], 50.0
        ),
        "median_reco_probe_mass2_GeV2": percentile_or_nan(
            arrays["reco_probe_mass2"], 50.0
        ),
        "p95_abs_reco_probe_mass2_GeV2": percentile_or_nan(
            np.abs(arrays["reco_probe_mass2"]), 95.0
        ),
        "median_reco_probe_E_minus_p_GeV": percentile_or_nan(
            arrays["reco_probe_E_minus_p"], 50.0
        ),
        "median_pred_probe_mass2_GeV2": percentile_or_nan(
            arrays["pred_probe_mass2"], 50.0
        ),
        "p95_abs_pred_probe_mass2_GeV2": percentile_or_nan(
            np.abs(arrays["pred_probe_mass2"]), 95.0
        ),
        "median_probe_direction_residual_deg": percentile_or_nan(
            arrays["probe_opening_residual_deg"], 50.0
        ),
        "p68_probe_direction_residual_deg": percentile_or_nan(
            arrays["probe_opening_residual_deg"], 68.0
        ),
        "p95_probe_direction_residual_deg": percentile_or_nan(
            arrays["probe_opening_residual_deg"], 95.0
        ),
        "median_probe_deltaE_GeV": percentile_or_nan(
            arrays["probe_delta_E"], 50.0
        ),
        "rms_probe_deltaE_GeV": rms_or_nan(
            arrays["probe_delta_E"]
        ),
        "median_probe_fractional_deltaE": percentile_or_nan(
            arrays["probe_delta_E_over_E"], 50.0
        ),
        "median_eppi0_missing_energy_GeV": percentile_or_nan(
            arrays["exclusivity_missing_energy"], 50.0
        ),
        "median_eppi0_missing_p_GeV": percentile_or_nan(
            arrays["exclusivity_missing_p"], 50.0
        ),
        "median_eppi0_missing_pT_GeV": percentile_or_nan(
            arrays["exclusivity_missing_pT"], 50.0
        ),
        "median_eppi0_missing_mass2_GeV2": percentile_or_nan(
            arrays["exclusivity_missing_mass2"], 50.0
        ),
    }

    if clean is not None and len(clean) == len(arrays["pred_probe_energy"]):
        nclean = int(np.count_nonzero(clean))
        summary["clean_association_pairs"] = nclean
        summary["clean_association_fraction"] = (
            nclean / len(clean) if len(clean) else float("nan")
        )
        if nclean:
            summary["clean_angle_q50_deg"] = quantile_or_nan(
                arrays["probe_opening_residual_deg"][clean], 0.50
            )
            summary["clean_angle_q68_deg"] = quantile_or_nan(
                arrays["probe_opening_residual_deg"][clean], 0.68
            )
            summary["clean_angle_q95_deg"] = quantile_or_nan(
                arrays["probe_opening_residual_deg"][clean], 0.95
            )
            summary["clean_angle_q99_deg"] = quantile_or_nan(
                arrays["probe_opening_residual_deg"][clean], 0.99
            )
        #endif
    #endif

    return summary


def write_summary_csv(summaries: List[Dict[str, object]], path: Path) -> None:
    if not summaries:
        return
    #endif

    keys: List[str] = []
    seen = set()
    for row in summaries:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
            #endif
        #endfor
    #endfor

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(summaries)



def read_csv_rows(path: Path) -> List[Dict[str, object]]:
    path = Path(path)
    if not path.exists():
        return []
    #endif
    with path.open(newline="") as f:
        return [dict(row) for row in csv.DictReader(f)]
    #endwith


def run_internal_self_tests(output_root: Path) -> None:
    """
    Fast startup tests for the two classes of late failures encountered during
    development: heterogeneous CSV row schemas and malformed MathText.
    """
    test_dir = Path(output_root) / ".photon_efficiency_selftest"
    test_dir.mkdir(parents=True, exist_ok=True)
    csv_path = test_dir / "rows.csv"
    png_path = test_dir / "mathtext.png"

    try:
        rows = [
            {"period": "selftest", "region": "FD_S1", "pi0_fraction": 0.70},
            {
                "period": "selftest",
                "region": "all",
                "pi0_fraction": 0.80,
                "closure_truth_20": 0.20,
                "closure_fit_20": 0.21,
                "closure_bias_20": 0.01,
            },
        ]
        write_rows_csv(rows, csv_path)
        loaded = read_csv_rows(csv_path)

        if len(loaded) != 2:
            raise RuntimeError("CSV self-test returned wrong row count.")
        #endif
        if "closure_bias_20" not in loaded[0]:
            raise RuntimeError("CSV self-test did not preserve union schema.")
        #endif
        if loaded[1].get("closure_bias_20", "") == "":
            raise RuntimeError("CSV self-test lost closure value.")
        #endif

        fig, ax = plt.subplots(figsize=(4.0, 3.0))
        ax.plot(
            [0.0, 1.0],
            [0.0, 1.0],
            label=r"raw $M_X^2\\otimes M_{probe}^2$",
        )
        ax.set_ylabel(r"$f_{\\pi^0}$")
        ax.legend()
        safe_finalize_figure(fig, png_path, rect=(0, 0, 1, 1))
        plt.close(fig)

        if not png_path.exists() or png_path.stat().st_size <= 0:
            raise RuntimeError("MathText-safe renderer self-test failed.")
        #endif

        log(
            "Internal self-tests passed: heterogeneous CSV writing and "
            "MathText-safe rendering."
        )
    finally:
        for candidate in (csv_path, png_path):
            try:
                if candidate.exists():
                    candidate.unlink()
                #endif
            except OSError:
                pass
            #endtry
        #endfor
        try:
            if test_dir.exists():
                test_dir.rmdir()
            #endif
        except OSError:
            pass
        #endtry
    #endtry


def aggregate_existing_outputs(
    selected: Sequence[PeriodConfig],
    stage2_outroot: Path,
    stage3_outroot: Path,
    include_stage3: bool,
    summary_outroot: Path,
) -> None:
    """Rebuild aggregate denominator/efficiency outputs from per-period files."""
    order = {p.key: i for i, p in enumerate(selected)}

    stage2_summaries: List[Dict[str, object]] = []
    disc_rows: List[Dict[str, object]] = []
    spread_rows: List[Dict[str, object]] = []
    mixed_rows: List[Dict[str, object]] = []
    shared_rows: List[Dict[str, object]] = []
    closure_rows: List[Dict[str, object]] = []
    ft_coarse_rows: List[Dict[str, object]] = []
    ft_coarse_closure_rows: List[Dict[str, object]] = []
    ft3_rows: List[Dict[str, object]] = []
    ft3_closure_rows: List[Dict[str, object]] = []
    ft3_summary_rows: List[Dict[str, object]] = []
    profile_rows: List[Dict[str, object]] = []
    control_rows: List[Dict[str, object]] = []
    stage3_rows: List[Dict[str, object]] = []
    stage3_summaries: List[Dict[str, object]] = []
    missing: List[str] = []

    for period in selected:
        p2 = Path(stage2_outroot) / period.key
        s2 = p2 / "stage2_summary.json"
        if s2.exists():
            with s2.open() as f:
                stage2_summaries.append(json.load(f))
            #endwith
        else:
            missing.append(str(s2))
        #endif

        disc_rows.extend(read_csv_rows(p2 / "denominator_discriminator_fits.csv"))
        spread_rows.extend(read_csv_rows(p2 / "denominator_discriminator_spread.csv"))
        mixed_rows.extend(read_csv_rows(p2 / "mixed_component_diagnostics.csv"))
        shared_rows.extend(read_csv_rows(p2 / "denominator_shared_morphed_fits.csv"))
        closure_rows.extend(read_csv_rows(p2 / "template_mixture_closure.csv"))
        ft_coarse_rows.extend(read_csv_rows(p2 / "ft_coarse_shared_morphed_fits.csv"))
        ft_coarse_closure_rows.extend(
            read_csv_rows(p2 / "ft_coarse_template_mixture_closure.csv")
        )
        ft3_rows.extend(
            read_csv_rows(p2 / "ft_coarse_three_component_by_mix_shift.csv")
        )
        ft3_closure_rows.extend(
            read_csv_rows(p2 / "ft_coarse_three_component_closure.csv")
        )
        ft3_summary_rows.extend(
            read_csv_rows(p2 / "ft_coarse_three_component_summary.csv")
        )
        profile_rows.extend(read_csv_rows(p2 / "morph_nuisance_profiles.csv"))
        control_rows.extend(read_csv_rows(p2 / "pi0_control_validation.csv"))

        if include_stage3:
            p3 = Path(stage3_outroot) / period.key
            s3 = p3 / "stage3_summary.json"
            if s3.exists():
                with s3.open() as f:
                    stage3_summaries.append(json.load(f))
                #endwith
            else:
                missing.append(str(s3))
            #endif
            stage3_rows.extend(
                read_csv_rows(p3 / "reference_efficiency_scale_factors.csv")
            )
        #endif
    #endfor

    if missing:
        raise FileNotFoundError(
            "--aggregate-only is missing required per-period files:\\n  "
            + "\\n  ".join(missing)
        )
    #endif

    def sort_rows(rows: List[Dict[str, object]]) -> None:
        rows.sort(
            key=lambda r: (
                order.get(str(r.get("period", "")), 999),
                str(r.get("region", "")),
                float(r.get("energy_low_GeV", "nan")),
                str(r.get("discriminator", "")),
                str(r.get("mix_shift", "")),
            )
        )
    #enddef

    stage2_summaries.sort(
        key=lambda r: order.get(str(r.get("period", "")), 999)
    )
    write_summary_csv(stage2_summaries, Path(stage2_outroot) / "stage2_summary.csv")

    for rows, filename in (
        (disc_rows, "denominator_discriminator_fits.csv"),
        (spread_rows, "denominator_discriminator_spread.csv"),
        (mixed_rows, "mixed_component_diagnostics.csv"),
        (shared_rows, "denominator_shared_morphed_fits.csv"),
        (closure_rows, "template_mixture_closure.csv"),
        (ft_coarse_rows, "ft_coarse_shared_morphed_fits.csv"),
        (ft_coarse_closure_rows, "ft_coarse_template_mixture_closure.csv"),
        (ft3_rows, "ft_coarse_three_component_by_mix_shift.csv"),
        (ft3_closure_rows, "ft_coarse_three_component_closure.csv"),
        (ft3_summary_rows, "ft_coarse_three_component_summary.csv"),
        (profile_rows, "morph_nuisance_profiles.csv"),
        (control_rows, "pi0_control_validation.csv"),
    ):
        if rows:
            sort_rows(rows)
            write_rows_csv(rows, Path(stage2_outroot) / filename)
        #endif
    #endfor

    if shared_rows and disc_rows:
        make_cross_period_stage2_composition_canvas(
            shared_rows, disc_rows, Path(stage2_outroot)
        )
    #endif

    if include_stage3 and stage3_rows:
        sort_rows(stage3_rows)
        write_rows_csv(
            stage3_rows,
            Path(stage3_outroot) / "reference_efficiency_scale_factors.csv",
        )
        make_cross_period_stage3_scale_factor_canvas(
            stage3_rows, Path(stage3_outroot)
        )
    #endif

    if include_stage3 and stage3_summaries:
        stage3_summaries.sort(
            key=lambda r: order.get(str(r.get("period", "")), 999)
        )
        write_summary_csv(
            stage3_summaries, Path(stage3_outroot) / "stage3_summary.csv"
        )
    #endif

    make_summary_dashboard(
        shared_rows,
        disc_rows,
        closure_rows,
        stage3_rows,
        ft3_summary_rows,
        Path(summary_outroot),
    )

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Photon-efficiency development: denominator composition plus preliminary "
            "data/MC photon-efficiency scale-factor extraction."
        )
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[p.key for p in PERIODS],
        help="Run only this period. May be supplied more than once. Default: all periods.",
    )
    parser.add_argument(
        "--output-dir",
        default="output/photon_efficiency",
        help="Photon-efficiency output root. Default: output/photon_efficiency.",
    )
    parser.add_argument(
        "--tree",
        default=None,
        help="TTree name. Default: auto-detect first TTree/RNTuple.",
    )
    parser.add_argument(
        "--max-entries",
        type=int,
        default=500_000,
        help="Maximum entries loaded from each ROOT tree. 0 means all. Default: 500000.",
    )
    parser.add_argument(
        "--angles",
        choices=("auto", "rad", "deg"),
        default="auto",
        help="Input angle units. Default: auto.",
    )
    parser.add_argument(
        "--tag-min",
        type=float,
        default=2.0,
        help="Minimum reconstructed epgamma tag energy in GeV. Default: 2.0.",
    )
    parser.add_argument(
        "--tag-max",
        type=float,
        default=9.5,
        help="Exclusive upper reconstructed epgamma tag energy in GeV. Default: 9.5.",
    )
    parser.add_argument(
        "--parent-component-tol",
        type=float,
        default=0.002,
        help=(
            "Maximum allowed absolute mismatch in each electron/proton Cartesian "
            "momentum component (GeV). Default: 0.002."
        ),
    )
    parser.add_argument(
        "--parent-distance-max",
        type=float,
        default=2.0,
        help=(
            "Maximum six-dimensional KD-tree distance after scaling each momentum "
            "component by --parent-component-tol. Default: 2.0."
        ),
    )
    parser.add_argument(
        "--max-tag-reco-angle",
        type=float,
        default=1.0,
        help=(
            "Deprecated in v004 and ignored. Stage I no longer applies an "
            "angular association cut; the complete residual distribution is "
            "retained for diagnosis."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help=(
            "Maximum number of independent period worker processes. Hard-capped "
            "at 8. Full-statistics jobs are memory/I/O intensive, so the default "
            "is 2 rather than five simultaneous periods."
        ),
    )
    parser.add_argument(
        "--kdtree-workers",
        type=int,
        default=1,
        help=(
            "Threads used inside each scipy cKDTree query. Default: 1. Keeping "
            "this at 1 avoids CPU oversubscription when periods run in parallel."
        ),
    )
    parser.add_argument(
        "--kdtree-query-chunk",
        type=int,
        default=500000,
        help=(
            "Maximum epgamma candidates per KD-tree query chunk. Default: 500000."
        ),
    )
    parser.add_argument(
        "--plot-mode",
        choices=("compact", "full"),
        default="compact",
        help=(
            "Plot-output mode. 'compact' writes the physics-summary canvases "
            "and a reduced diagnostic set; 'full' additionally restores the "
            "legacy development canvases. Default: compact."
        ),
    )
    parser.add_argument(
        "--save-npz",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--assoc-mgg-min",
        type=float,
        default=0.11,
        help="Lower clean-association M_gammagamma edge (GeV). Default: 0.11 (matches the upstream eppi0 skim).",
    )
    parser.add_argument(
        "--assoc-mgg-max",
        type=float,
        default=0.16,
        help="Upper clean-association M_gammagamma edge (GeV). Default: 0.16 (matches the upstream eppi0 skim).",
    )
    parser.add_argument(
        "--assoc-remainder-mass2-max",
        type=float,
        default=1.0e-3,
        help=(
            "Require |(P_pi0-k_tag)^2| below this value (GeV^2) for a clean "
            "tag association. Default: 1e-3."
        ),
    )
    parser.add_argument(
        "--assoc-probe-energy-min",
        type=float,
        default=0.40,
        help=(
            "Minimum reconstructed companion-photon energy (GeV) for the clean "
            "association sample. Default: 0.40."
        ),
    )
    parser.add_argument(
        "--ft-theta-max",
        type=float,
        default=5.5,
        help=(
            "Predicted-probe polar-angle boundary (deg). FT is theta_pred <= this "
            "value; FD is above it. Default: 5.5."
        ),
    )
    parser.add_argument(
        "--min-resolution-count",
        type=int,
        default=100,
        help=(
            "Minimum clean entries required to report quantiles for one "
            "energy/region bin. Default: 100."
        ),
    )

    # Stage-II denominator-composition controls.
    parser.add_argument(
        "--stage1-only",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--stage2-output-dir",
        default="output/photon_efficiency/stage2",
        help="Stage-II output directory. Default: output/photon_efficiency/stage2.",
    )
    parser.add_argument(
        "--den-fit-mm2-min",
        type=float,
        default=-0.08,
        help="Lower M_X^2(ep) fit edge (GeV^2). Default: -0.08.",
    )
    parser.add_argument(
        "--den-fit-mm2-max",
        type=float,
        default=0.10,
        help="Upper M_X^2(ep) fit edge (GeV^2). Default: 0.11 (matches the upstream eppi0 skim).",
    )
    parser.add_argument(
        "--den-fit-probe-m2-max",
        type=float,
        default=0.60,
        help=(
            "Fit |M_X^2(epgamma_tag)| up to this value (GeV^2). "
            "Default: 0.60."
        ),
    )
    parser.add_argument(
        "--den-fit-mm2-bins",
        type=int,
        default=48,
        help="Number of M_X^2(ep) bins in each Stage-II template fit. Default: 48.",
    )
    parser.add_argument(
        "--den-fit-probe-m2-bins",
        type=int,
        default=48,
        help=(
            "Number of predicted-probe M^2 bins in each Stage-II template fit. "
            "Default: 48."
        ),
    )
    parser.add_argument(
        "--den-min-data-count",
        type=int,
        default=250,
        help=(
            "Minimum data entries inside a region/energy fit window before a "
            "Stage-II component fit is attempted. Default: 250."
        ),
    )
    parser.add_argument(
        "--den-min-template-count",
        type=int,
        default=100,
        help=(
            "Minimum aaogen and dvcsgen template entries in a region/energy fit "
            "window. Default: 100."
        ),
    )
    parser.add_argument(
        "--den-probe-energy-max",
        type=float,
        default=2.0,
        help=(
            "Maximum predicted probe energy (GeV) used for Stage-II denominator "
            "fits. The real epgamma skim has an observed missing-energy endpoint "
            "near 2 GeV, so 2.0 is the default supported upper boundary."
        ),
    )
    parser.add_argument(
        "--mix-shifts",
        default="7919,15401,32749",
        help=(
            "Comma-separated deterministic cyclic shifts used only for the "
            "mixed-data wrong-tag diagnostic. Mixed events are NOT part of the "
            "nominal denominator fit. Default: 7919,15401,32749."
        ),
    )
    parser.add_argument(
        "--disc-1d-bins",
        type=int,
        default=90,
        help=(
            "Histogram bins for 1D Stage-II discriminator fits. Default: 90."
        ),
    )
    parser.add_argument(
        "--disc-ptmiss-max",
        type=float,
        default=1.0,
        help=(
            "Upper stored pTmiss edge (GeV) used by the optional "
            "M_X^2(ep) x pTmiss discriminator. Default: 1.0."
        ),
    )
    parser.add_argument(
        "--disc-ptmiss-bins",
        type=int,
        default=40,
        help=(
            "pTmiss bins for the optional 2D discriminator fit. Default: 40."
        ),
    )
    parser.add_argument(
        "--disc-theta-max",
        type=float,
        default=40.0,
        help=(
            "Upper theta_gamma_gamma edge (deg) used by the optional "
            "M_X^2(ep) x theta_gamma_gamma discriminator. Default: 40."
        ),
    )
    parser.add_argument(
        "--disc-theta-bins",
        type=int,
        default=40,
        help=(
            "theta_gamma_gamma bins for the optional 2D discriminator fit. "
            "Default: 40."
        ),
    )
    parser.add_argument(
        "--control-mgg-min",
        type=float,
        default=0.11,
        help="Lower reconstructed-pi0 control M_gammagamma edge (GeV). Default: 0.11 (matches the upstream eppi0 skim).",
    )
    parser.add_argument(
        "--control-mgg-max",
        type=float,
        default=0.16,
        help="Upper reconstructed-pi0 control M_gammagamma edge (GeV). Default: 0.16 (matches the upstream eppi0 skim).",
    )
    parser.add_argument(
        "--control-emiss-abs-max",
        type=float,
        default=0.50,
        help=(
            "Loose |Emiss2| upper bound (GeV) for the final reconstructed-pi0 "
            "control sample. pTmiss itself is never cut when calibrating pTmiss. "
            "Default: 0.50."
        ),
    )
    parser.add_argument(
        "--control-pi0-energy-edges",
        default="2.4,3.0,4.0,5.0,7.0,10.0",
        help=(
            "Comma-separated reconstructed pi0-energy edges (GeV) for control "
            "validation. Individual daughter-photon energies are not stored in "
            "the eppi0 tree, so v018 does not pretend these are photon-energy bins."
        ),
    )
    parser.add_argument("--morph-shift-prior-bins", type=float, default=1.0, help="Gaussian prior width for template shifts in histogram bins. Default: 1.0.")
    parser.add_argument("--morph-sigma-prior-bins", type=float, default=1.0, help="Gaussian prior width for additional template smearing in histogram bins. Default: 1.0.")
    parser.add_argument("--morph-max-shift-bins", type=float, default=4.0, help="Hard bound on template shifts in histogram bins. Default: 4.0.")
    parser.add_argument("--morph-max-sigma-bins", type=float, default=6.0, help="Hard bound on additional template smearing in histogram bins. Default: 6.0.")
    parser.add_argument(
        "--ft-coarse-energy-edges",
        default="0.40,1.20,2.00",
        help=(
            "Comma-separated predicted-probe energy edges (GeV) for the "
            "diagnostic coarse-bin FT composition study. This does not alter "
            "the nominal Stage-II or Stage-III binning. "
            "Default: 0.40,1.20,2.00."
        ),
    )
    parser.add_argument(
        "--stage2-closure-truth-fractions",
        default="0.2,0.5,0.8",
        help=(
            "Comma-separated injected pi0 fractions for the FT/FD_all Stage-II "
            "morphed-template closure study. Default: 0.2,0.5,0.8."
        ),
    )
    parser.add_argument(
        "--morph-profile-shift-grid",
        default="-4,-3,-2,-1,0,1,2,3,4",
        help=(
            "Comma-separated fixed shift values (histogram bins) used for "
            "all-region nuisance profile scans. Default: -4,-3,-2,-1,0,1,2,3,4."
        ),
    )
    parser.add_argument(
        "--morph-profile-sigma-grid",
        default="0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6",
        help=(
            "Comma-separated fixed additional-smearing values (histogram bins) "
            "used for all-region nuisance profile scans. "
            "Default: 0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6."
        ),
    )
    parser.add_argument(
        "--skip-stage3",
        action="store_true",
        help=(
            "Skip the preliminary Stage-III efficiency/scale-factor reference "
            "extraction. Stage III runs by default whenever Stage II runs."
        ),
    )
    parser.add_argument(
        "--stage3-tag-remainder-m2-max",
        type=float,
        default=1.0e-3,
        help=(
            "Maximum |(P_pi0-k_tag)^2| (GeV^2) for the current-data "
            "kinematic tag-in-pi0 identity test. Default: 1e-3."
        ),
    )
    parser.add_argument(
        "--stage3-probe-angle-max-deg",
        type=float,
        default=10.0,
        help=(
            "Loose maximum opening angle (deg) between the reconstructed "
            "companion P_pi0-k_tag and the epgamma-predicted probe. This is an "
            "association diagnostic/gate, not an exclusivity cut. Default: 10."
        ),
    )
    parser.add_argument(
        "--stage3-probe-frac-energy-max",
        type=float,
        default=0.50,
        help=(
            "Loose maximum |E_reco-E_pred|/E_pred for the final current-data "
            "probe-consistency association. Default: 0.50."
        ),
    )
    parser.add_argument(
        "--stage3-use-pred-consistency",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Use the loose predicted-vs-reconstructed companion consistency "
            "gate in the Stage-III data numerator. Disable with "
            "--no-stage3-use-pred-consistency for a mass-shell-only diagnostic. "
            "Default: enabled."
        ),
    )
    parser.add_argument(
        "--stage3-output-dir",
        default="output/photon_efficiency/stage3",
        help=(
            "Stage-III output directory. Default: "
            "output/photon_efficiency/stage3."
        ),
    )
    parser.add_argument(
        "--aggregate-only",
        action="store_true",
        help=(
            "Do not reread ROOT files. Rebuild top-level aggregate CSV/JSON "
            "outputs from completed per-period files."
        ),
    )
    return parser.parse_args()



def process_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
    output_dir: str,
    stage2_output_dir: str,
) -> Dict[str, object]:
    """
    Process one run period completely.

    This top-level function is compatible with ProcessPoolExecutor and avoids
    relying on fork-specific behavior.
    """
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    tree_name = args_dict["tree"]
    max_entries = int(args_dict["max_entries"])
    angle_mode = str(args_dict["angles"])
    tag_min = float(args_dict["tag_min"])
    tag_max = float(args_dict["tag_max"])
    parent_component_tol = float(args_dict["parent_component_tol"])
    parent_distance_max = float(args_dict["parent_distance_max"])
    kdtree_workers = int(args_dict["kdtree_workers"])

    outroot = Path(output_dir)

    t0 = time.perf_counter()

    log(f"{period.label}: reading aaogen epgamma sample.")
    epg_arrays, epg_tree, epg_total = read_branches(
        period.epgamma_pi0_mc,
        EPG_REQUIRED,
        EPG_OPTIONAL_PI0_MC,
        tree_name,
        max_entries,
    )
    epg = extract_epgamma(epg_arrays, angle_mode)
    del epg_arrays

    log(
        f"{period.label}: epgamma tree '{epg_tree}', "
        f"loaded {len(epg.tag_energy):,}/{epg_total:,} entries, "
        f"angles interpreted as {epg.angle_unit}."
    )

    log(f"{period.label}: reading reconstructed aaogen eppi0 sample.")
    pi_arrays, pi_tree, pi_total = read_branches(
        period.eppi0_pi0_mc,
        EPPIO_REQUIRED,
        EPPIO_OPTIONAL_PI0_MC,
        tree_name,
        max_entries,
    )
    eppi0 = extract_eppi0(pi_arrays, angle_mode)
    del pi_arrays

    log(
        f"{period.label}: eppi0 tree '{pi_tree}', "
        f"loaded {len(eppi0.pi0_p):,}/{pi_total:,} entries, "
        f"angles interpreted as {eppi0.angle_unit}."
    )

    log(f"{period.label}: reading reconstructed eppi0 data control sample.")
    pi_data_arrays, pi_data_tree, pi_data_total = read_branches(
        period.eppi0_data, EPPIO_REQUIRED, EPPIO_OPTIONAL_DATA, tree_name, max_entries
    )
    eppi0_data = extract_eppi0(pi_data_arrays, angle_mode)
    del pi_data_arrays
    gc.collect()
    log(f"{period.label}: eppi0 data tree '{pi_data_tree}', loaded {len(eppi0_data.pi0_p):,}/{pi_data_total:,} entries.")

    log(f"{period.label}: matching e/p parent kinematics between MC skims.")
    matches = match_parent_kinematics(
        epg,
        eppi0,
        tag_min=tag_min,
        tag_max=tag_max,
        component_tolerance=parent_component_tol,
        nearest_distance_max=parent_distance_max,
        kdtree_workers=kdtree_workers,
        query_chunk_size=int(args_dict["kdtree_query_chunk"]),
    )
    log(f"{period.label}: accepted {len(matches.epg_index):,} parent matches.")

    log(f"{period.label}: vectorized construction of reconstructed and predicted probes.")
    pair_np, counters = build_stage1_arrays(
        period,
        epg,
        eppi0,
        matches,
    )

    clean = build_clean_association_mask(
        pair_np,
        mgg_min=float(args_dict["assoc_mgg_min"]),
        mgg_max=float(args_dict["assoc_mgg_max"]),
        remainder_mass2_max=float(args_dict["assoc_remainder_mass2_max"]),
        reco_probe_energy_min=float(args_dict["assoc_probe_energy_min"]),
    )

    mc_association_summary = {
        "period": period.key,
        "label": period.label,
        "aaogen_epgamma_loaded": int(len(epg.tag_energy)),
        "aaogen_eppi0_loaded": int(len(eppi0.pi0_p)),
        "accepted_parent_matches": int(len(matches.epg_index)),
        "accepted_probe_pairs": int(counters["accepted_stage1_pairs"]),
        "clean_reconstructed_probes": int(np.count_nonzero(clean)),
        "invalid_pi0_mass": int(counters["invalid_pi0_mass"]),
        "nonphysical_reco_probe_energy": int(
            counters["nonphysical_reco_probe_energy"]
        ),
    }
    log(
        f"{period.label}: internal MC tag/probe association: "
        f"{len(matches.epg_index):,} parent matches; "
        f"{np.count_nonzero(clean):,} clean reconstructed companion probes."
    )

    stage2_summary: Optional[Dict[str, object]] = None
    stage2_rows: List[Dict[str, object]] = []
    stage2_spread_rows: List[Dict[str, object]] = []
    mixed_diag_rows: List[Dict[str, object]] = []
    profile_rows: List[Dict[str, object]] = []
    shared_rows: List[Dict[str, object]] = []
    closure_rows: List[Dict[str, object]] = []
    ft_coarse_rows: List[Dict[str, object]] = []
    ft_coarse_closure_rows: List[Dict[str, object]] = []
    ft_coarse_three_component_rows: List[Dict[str, object]] = []
    ft_coarse_three_component_closure_rows: List[Dict[str, object]] = []
    ft_coarse_three_component_summary_rows: List[Dict[str, object]] = []
    control_rows: List[Dict[str, object]] = []
    stage3_rows: List[Dict[str, object]] = []
    stage3_summary: Optional[Dict[str, object]] = None

    if True:
        stage2_dir = Path(stage2_output_dir) / period.key
        stage2_dir.mkdir(parents=True, exist_ok=True)
        # Clear stale canvases exactly once BEFORE generating any Stage-II plots.
        # Older versions deleted *.png inside make_stage2_canvases(), which
        # erased the new control-validation and shared-fit canvases.
        for old_png in stage2_dir.glob("*.png"):
            old_png.unlink()
        #endfor

        log(f"{period.label}: Stage II reading real epgamma data.")
        data_arrays, data_tree, data_total = read_branches(
            period.epgamma_data,
            EPG_REQUIRED,
            EPG_OPTIONAL_DATA,
            tree_name,
            max_entries,
        )
        data_epg = extract_epgamma(data_arrays, angle_mode)
        del data_arrays
        log(
            f"{period.label}: data epgamma tree '{data_tree}', "
            f"loaded {len(data_epg.tag_energy):,}/{data_total:,} entries."
        )

        log(f"{period.label}: Stage II reading dvcsgen epgamma MC.")
        dv_arrays, dv_tree, dv_total = read_branches(
            period.epgamma_dvcs_mc,
            EPG_REQUIRED,
            EPG_OPTIONAL_DVCS_MC,
            tree_name,
            max_entries,
        )
        dv_epg = extract_epgamma(dv_arrays, angle_mode)
        del dv_arrays
        log(
            f"{period.label}: dvcsgen tree '{dv_tree}', "
            f"loaded {len(dv_epg.tag_energy):,}/{dv_total:,} entries."
        )

        ft_theta_max = float(args_dict["ft_theta_max"])
        max_probe_energy = float(args_dict["den_probe_energy_max"])

        # Reuse aaogen epgamma already resident from Stage I.
        pi0_f = build_epgamma_denominator_features(
            period, epg, tag_min, tag_max, ft_theta_max
        )
        data_f = build_epgamma_denominator_features(
            period, data_epg, tag_min, tag_max, ft_theta_max
        )
        dvcs_f = build_epgamma_denominator_features(
            period, dv_epg, tag_min, tag_max, ft_theta_max
        )
        del dv_epg
        gc.collect()

        available = [
            disc for disc in STAGE2_DISCRIMINATORS
            if discriminator_available(disc, data_f, pi0_f, dvcs_f)
        ]
        log(
            f"{period.label}: Stage II discriminator study for "
            f"0.4 <= E_probe^pred < {max_probe_energy:g} GeV: "
            + ", ".join(available)
        )

        stage2_rows = run_stage2_discriminator_fits(
            period,
            data_f=data_f,
            pi0_f=pi0_f,
            dvcs_f=dvcs_f,
            ft_theta_max=ft_theta_max,
            max_probe_energy=max_probe_energy,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
            mm2_bins_2d=int(args_dict["den_fit_mm2_bins"]),
            probe_m2_bins_2d=int(args_dict["den_fit_probe_m2_bins"]),
            bins_1d=int(args_dict["disc_1d_bins"]),
            ptmiss_max=float(args_dict["disc_ptmiss_max"]),
            ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
            theta_max=float(args_dict["disc_theta_max"]),
            theta_bins=int(args_dict["disc_theta_bins"]),
            min_data_count=int(args_dict["den_min_data_count"]),
            min_template_count=int(args_dict["den_min_template_count"]),
        )
        write_rows_csv(
            stage2_rows,
            stage2_dir / "denominator_discriminator_fits.csv",
        )

        control_energy_edges = parse_float_edges(
            str(args_dict["control_pi0_energy_edges"]),
            "--control-pi0-energy-edges",
        )
        pi0_control, control_rows = build_pi0_control_validation(
            period,
            eppi0_data,
            eppi0,
            mgg_min=float(args_dict["control_mgg_min"]),
            mgg_max=float(args_dict["control_mgg_max"]),
            emiss_abs_max=float(args_dict["control_emiss_abs_max"]),
            ft_theta_max=ft_theta_max,
            pi0_energy_edges=control_energy_edges,
            ptmiss_max=float(args_dict["disc_ptmiss_max"]),
            ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
        )
        write_rows_csv(control_rows, stage2_dir / "pi0_control_validation.csv")
        with (stage2_dir / "pi0_control_calibration.json").open("w") as f:
            json.dump(pi0_control, f, indent=2, allow_nan=True)
        #endwith
        if not compact_plot_enabled(args_dict):
            make_pi0_control_canvases(
                period,
                eppi0_data,
                eppi0,
                control_rows,
                stage2_dir,
                mgg_min=float(args_dict["control_mgg_min"]),
                mgg_max=float(args_dict["control_mgg_max"]),
                emiss_abs_max=float(args_dict["control_emiss_abs_max"]),
                ft_theta_max=ft_theta_max,
                pi0_energy_edges=control_energy_edges,
                ptmiss_max=float(args_dict["disc_ptmiss_max"]),
            )
        #endif
        log(
            f"{period.label}: reconstructed-eppi0 control validation written; "
            "control calibration is diagnostic-only and is not used as a morph prior."
        )
        del eppi0
        gc.collect()
        # Keep the object through the transfer-factor calculation below; it is
        # still modest compared with the epgamma MC samples.
        shared_rows, closure_rows = run_stage2_shared_morphed_fits(
            period, data_f, pi0_f, dvcs_f, pi0_control, ft_theta_max, max_probe_energy,
            float(args_dict["den_fit_mm2_min"]), float(args_dict["den_fit_mm2_max"]), float(args_dict["den_fit_probe_m2_max"]),
            int(args_dict["den_fit_mm2_bins"]), int(args_dict["den_fit_probe_m2_bins"]), float(args_dict["disc_ptmiss_max"]),
            int(args_dict["disc_ptmiss_bins"]), float(args_dict["disc_theta_max"]), int(args_dict["disc_theta_bins"]),
            int(args_dict["den_min_data_count"]), int(args_dict["den_min_template_count"]),
            float(args_dict["morph_shift_prior_bins"]), float(args_dict["morph_sigma_prior_bins"]),
            float(args_dict["morph_max_shift_bins"]), float(args_dict["morph_max_sigma_bins"]),
            parse_float_edges(
                str(args_dict["stage2_closure_truth_fractions"]),
                "--stage2-closure-truth-fractions",
            ),
        )
        write_rows_csv(shared_rows, stage2_dir / "denominator_shared_morphed_fits.csv")
        write_rows_csv(closure_rows, stage2_dir / "template_mixture_closure.csv")
        if not compact_plot_enabled(args_dict):
            make_stage2_template_mixture_closure_canvas(
                closure_rows,
                stage2_dir,
            )
        #endif

        ft_coarse_edges = parse_float_edges(
            str(args_dict["ft_coarse_energy_edges"]),
            "--ft-coarse-energy-edges",
        )
        ft_coarse_rows, ft_coarse_closure_rows = (
            run_stage2_ft_coarse_shared_morphed_fits(
                period,
                data_f,
                pi0_f,
                dvcs_f,
                pi0_control,
                ft_theta_max,
                ft_coarse_edges,
                float(args_dict["den_fit_mm2_min"]),
                float(args_dict["den_fit_mm2_max"]),
                float(args_dict["den_fit_probe_m2_max"]),
                int(args_dict["den_fit_mm2_bins"]),
                int(args_dict["den_fit_probe_m2_bins"]),
                float(args_dict["disc_ptmiss_max"]),
                int(args_dict["disc_ptmiss_bins"]),
                float(args_dict["disc_theta_max"]),
                int(args_dict["disc_theta_bins"]),
                int(args_dict["den_min_data_count"]),
                int(args_dict["den_min_template_count"]),
                float(args_dict["morph_shift_prior_bins"]),
                float(args_dict["morph_sigma_prior_bins"]),
                float(args_dict["morph_max_shift_bins"]),
                float(args_dict["morph_max_sigma_bins"]),
                parse_float_edges(
                    str(args_dict["stage2_closure_truth_fractions"]),
                    "--stage2-closure-truth-fractions",
                ),
            )
        )
        write_rows_csv(
            ft_coarse_rows,
            stage2_dir / "ft_coarse_shared_morphed_fits.csv",
        )
        write_rows_csv(
            ft_coarse_closure_rows,
            stage2_dir / "ft_coarse_template_mixture_closure.csv",
        )
        if not compact_plot_enabled(args_dict):
            make_ft_coarse_composition_canvas(
                period,
                shared_rows,
                ft_coarse_rows,
                ft_coarse_closure_rows,
                stage2_dir,
            )
        #endif


        log(
            f"{period.label}: wrote per-period shared-morphed numerical results "
            f"({len(shared_rows)} rows) before canvas rendering."
        )
        if not compact_plot_enabled(args_dict):
            make_shared_fit_canvas(period, shared_rows, stage2_rows, pi0_control, stage2_dir)
        #endif
        make_ft_fd_composition_canvas(
            period,
            shared_rows,
            stage2_dir,
        )
        make_ft_fd_fit_overlay_canvas(
            period,
            shared_rows,
            data_f,
            pi0_f,
            dvcs_f,
            stage2_dir,
            ft_theta_max=ft_theta_max,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
            mm2_bins_2d=int(args_dict["den_fit_mm2_bins"]),
            probe_m2_bins_2d=int(args_dict["den_fit_probe_m2_bins"]),
            ptmiss_max=float(args_dict["disc_ptmiss_max"]),
            ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
            theta_max=float(args_dict["disc_theta_max"]),
            theta_bins=int(args_dict["disc_theta_bins"]),
        )

        shift_profile_grid = np.asarray(
            parse_float_edges(
                str(args_dict["morph_profile_shift_grid"]),
                "--morph-profile-shift-grid",
            ),
            dtype=float,
        )
        sigma_profile_grid = np.asarray(
            parse_float_edges(
                str(args_dict["morph_profile_sigma_grid"]),
                "--morph-profile-sigma-grid",
            ),
            dtype=float,
        )
        profile_rows = run_stage2_morph_profile_scans(
            period,
            data_f=data_f,
            pi0_f=pi0_f,
            dvcs_f=dvcs_f,
            ft_theta_max=ft_theta_max,
            max_probe_energy=max_probe_energy,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
            mm2_bins_2d=int(args_dict["den_fit_mm2_bins"]),
            probe_m2_bins_2d=int(args_dict["den_fit_probe_m2_bins"]),
            ptmiss_max=float(args_dict["disc_ptmiss_max"]),
            ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
            theta_max=float(args_dict["disc_theta_max"]),
            theta_bins=int(args_dict["disc_theta_bins"]),
            min_data_count=int(args_dict["den_min_data_count"]),
            min_template_count=int(args_dict["den_min_template_count"]),
            shift_grid=shift_profile_grid,
            sigma_grid=sigma_profile_grid,
        )
        write_rows_csv(
            profile_rows,
            stage2_dir / "morph_nuisance_profiles.csv",
        )
        make_morph_profile_canvas(
            period,
            profile_rows,
            stage2_dir,
        )
        log(
            f"{period.label}: wrote explicit all-region morph nuisance profiles "
            f"({len(profile_rows)} scan rows)."
        )

        stage2_spread_rows = discriminator_spread_rows(stage2_rows)
        write_rows_csv(
            stage2_spread_rows,
            stage2_dir / "denominator_discriminator_spread.csv",
        )

        # Mixed-event stress tests remain nominal-discriminator-only and are
        # processed one shift at a time to bound peak memory.
        for shift in parse_mix_shifts(str(args_dict["mix_shifts"])):
            log(
                f"{period.label}: mixed-data diagnostic shift {shift} "
                "(nominal discriminator stress test only)."
            )
            mixed_f = build_epgamma_denominator_features(
                period,
                data_epg,
                tag_min=tag_min,
                tag_max=tag_max,
                ft_theta_max=ft_theta_max,
                mixed_tag_shift=shift,
            )
            mixed_diag_rows.extend(
                run_mixed_shift_diagnostic(
                    period,
                    data_f=data_f,
                    pi0_f=pi0_f,
                    dvcs_f=dvcs_f,
                    mixed_f=mixed_f,
                    fit_rows=stage2_rows,
                    shift=shift,
                    ft_theta_max=ft_theta_max,
                    mm2_min=float(args_dict["den_fit_mm2_min"]),
                    mm2_max=float(args_dict["den_fit_mm2_max"]),
                    probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
                    mm2_bins=int(args_dict["den_fit_mm2_bins"]),
                    probe_m2_bins=int(args_dict["den_fit_probe_m2_bins"]),
                    min_template_count=int(args_dict["den_min_template_count"]),
                )
            )
            three_rows_shift, three_closure_shift = (
                run_ft_coarse_three_component_diagnostic(
                    period,
                    data_f=data_f,
                    pi0_f=pi0_f,
                    dvcs_f=dvcs_f,
                    mixed_f=mixed_f,
                    coarse_rows=ft_coarse_rows,
                    shift=shift,
                    ft_theta_max=ft_theta_max,
                    mm2_min=float(args_dict["den_fit_mm2_min"]),
                    mm2_max=float(args_dict["den_fit_mm2_max"]),
                    probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
                    mm2_bins_2d=int(args_dict["den_fit_mm2_bins"]),
                    probe_m2_bins_2d=int(args_dict["den_fit_probe_m2_bins"]),
                    ptmiss_max=float(args_dict["disc_ptmiss_max"]),
                    ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
                    theta_max=float(args_dict["disc_theta_max"]),
                    theta_bins=int(args_dict["disc_theta_bins"]),
                    min_data_count=int(args_dict["den_min_data_count"]),
                    min_template_count=int(args_dict["den_min_template_count"]),
                )
            )
            ft_coarse_three_component_rows.extend(three_rows_shift)
            ft_coarse_three_component_closure_rows.extend(three_closure_shift)
            del mixed_f
        #endfor

        write_rows_csv(
            mixed_diag_rows,
            stage2_dir / "mixed_component_diagnostics.csv",
        )

        ft_coarse_three_component_summary_rows = summarize_ft_coarse_three_component_shifts(
            ft_coarse_three_component_rows
        )
        write_rows_csv(ft_coarse_three_component_rows, stage2_dir / "ft_coarse_three_component_by_mix_shift.csv")
        write_rows_csv(ft_coarse_three_component_closure_rows, stage2_dir / "ft_coarse_three_component_closure.csv")
        write_rows_csv(ft_coarse_three_component_summary_rows, stage2_dir / "ft_coarse_three_component_summary.csv")
        if not compact_plot_enabled(args_dict):
            make_ft_coarse_three_component_canvas(
                period,
                ft_coarse_three_component_summary_rows,
                ft_coarse_three_component_closure_rows,
                stage2_dir,
            )
        #endif

        if not compact_plot_enabled(args_dict):
            make_stage2_canvases(
                period,
                data_f=data_f,
                pi0_f=pi0_f,
                dvcs_f=dvcs_f,
                fit_rows=stage2_rows,
                spread_rows=stage2_spread_rows,
                mixed_rows=mixed_diag_rows,
                outdir=stage2_dir,
                ft_theta_max=ft_theta_max,
                max_probe_energy=max_probe_energy,
                mm2_min=float(args_dict["den_fit_mm2_min"]),
                mm2_max=float(args_dict["den_fit_mm2_max"]),
                probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
            )
        #endif

        stage2_summary = summarize_stage2(
            period,
            fit_rows=stage2_rows,
            spread_rows=stage2_spread_rows,
            mixed_rows=mixed_diag_rows,
            data_total=data_total,
            pi0_total=epg_total,
            dvcs_total=dv_total,
            max_probe_energy=max_probe_energy,
        )
        with (stage2_dir / "stage2_summary.json").open("w") as f:
            json.dump(stage2_summary, f, indent=2, allow_nan=True)
        #endwith

        log(
            f"{period.label}: Stage II completed "
            f"{stage2_summary['successful_nominal_all_region_energy_fits']} "
            "successful nominal all-region fits."
        )

        if not bool(args_dict.get("skip_stage3", False)):
            stage3_dir = Path(args_dict["stage3_output_dir"]) / period.key
            stage3_dir.mkdir(parents=True, exist_ok=True)
            for old_png in stage3_dir.glob("*.png"):
                old_png.unlink()
            #endfor

            log(
                f"{period.label}: Stage III exact runnum/evnum matching of "
                "epgamma data tags to reconstructed eppi0 data candidates."
            )
            data_assoc = match_data_event_candidates(
                period,
                data_epg,
                eppi0_data,
                tag_remainder_m2_max=float(
                    args_dict["stage3_tag_remainder_m2_max"]
                ),
                reco_probe_energy_min=float(args_dict["assoc_probe_energy_min"]),
                probe_angle_max_deg=float(
                    args_dict["stage3_probe_angle_max_deg"]
                ),
                probe_frac_energy_max=float(
                    args_dict["stage3_probe_frac_energy_max"]
                ),
            )

            # Best mass-shell/threshold candidate lookup. The final numerator
            # optionally adds the loose predicted/reconstructed consistency gate.
            data_mass_shell_lookup = data_assoc.stage_lookup["probe_energy"]
            if bool(args_dict["stage3_use_pred_consistency"]):
                data_final_lookup = data_assoc.stage_lookup[
                    "probe_pred_consistent"
                ]
            else:
                data_final_lookup = data_mass_shell_lookup.copy()
            #endif

            # MC uses the internal aaogen e/p parent matcher.
            mc_any_lookup, mc_clean_lookup = stage1_success_lookup(
                len(epg.tag_energy),
                pair_np,
                clean,
            )

            stage3_rows = build_stage3_reference_rows(
                period,
                shared_rows=shared_rows,
                data_f=data_f,
                aaogen_f=pi0_f,
                data_stage_lookup=data_assoc.stage_lookup,
                data_mass_shell_success=data_mass_shell_lookup,
                data_final_success=data_final_lookup,
                mc_any_match=mc_any_lookup,
                mc_clean_success=mc_clean_lookup,
                ft_theta_max=ft_theta_max,
                max_probe_energy=max_probe_energy,
                mm2_min=float(args_dict["den_fit_mm2_min"]),
                mm2_max=float(args_dict["den_fit_mm2_max"]),
                probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
            )
            write_rows_csv(
                stage3_rows,
                stage3_dir / "reference_efficiency_scale_factors.csv",
            )
            make_stage3_canvases(
                period,
                stage3_rows,
                stage3_dir,
                compact=compact_plot_enabled(args_dict),
            )
            stage3_summary = summarize_stage3(
                period,
                stage3_rows,
                data_match_counters={
                    **data_assoc.counters,
                    "mass_shell_above_threshold_associations": int(
                        np.count_nonzero(data_mass_shell_lookup)
                    ),
                    "final_numerator_associations": int(
                        np.count_nonzero(data_final_lookup)
                    ),
                },
                max_probe_energy=max_probe_energy,
                use_pred_consistency=bool(
                    args_dict["stage3_use_pred_consistency"]
                ),
                association_settings={
                    "tag_remainder_m2_max_GeV2": float(
                        args_dict["stage3_tag_remainder_m2_max"]
                    ),
                    "reco_probe_energy_min_GeV": float(
                        args_dict["assoc_probe_energy_min"]
                    ),
                    "probe_angle_max_deg": float(
                        args_dict["stage3_probe_angle_max_deg"]
                    ),
                    "probe_frac_energy_max": float(
                        args_dict["stage3_probe_frac_energy_max"]
                    ),
                },
            )
            with (stage3_dir / "stage3_summary.json").open("w") as f:
                json.dump(stage3_summary, f, indent=2, allow_nan=True)
            #endwith
            log(
                f"{period.label}: Stage III wrote {len(stage3_rows)} reference "
                "efficiency rows. These are validation results, not the final "
                "production scale factors."
            )
        #endif
    #endif

    summary["wall_time_s"] = float(time.perf_counter() - t0)
    with (period_dir / "stage1_summary.json").open("w") as f:
        json.dump(summary, f, indent=2, allow_nan=True)
    #endwith

    log(f"{period.label}: total worker wall time = {summary['wall_time_s']:.1f} s.")

    return {
        "mc_association_summary": mc_association_summary,
        "stage2_summary": stage2_summary,
        "stage2_rows": stage2_rows,
        "stage2_spread_rows": stage2_spread_rows,
        "mixed_diag_rows": mixed_diag_rows,
        "profile_rows": profile_rows,
        "shared_rows": shared_rows,
        "closure_rows": closure_rows,
        "ft_coarse_rows": ft_coarse_rows,
        "ft_coarse_closure_rows": ft_coarse_closure_rows,
        "ft_coarse_three_component_rows": ft_coarse_three_component_rows,
        "ft_coarse_three_component_closure_rows": ft_coarse_three_component_closure_rows,
        "ft_coarse_three_component_summary_rows": ft_coarse_three_component_summary_rows,
        "control_rows": control_rows,
        "stage3_rows": stage3_rows,
        "stage3_summary": stage3_summary,
    }



def make_cross_period_stage3_scale_factor_canvas(
    stage3_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """All-period Stage-III scale-factor summary for FT and FD_all."""
    if not stage3_rows:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.8))

    region_specs = (
        ("FT", "FT", 1.8),
        ("FD_all", "FD all", 1.35),
    )

    for ax, region, title, ymax in (
        (axes[0], *region_specs[0]),
        (axes[1], *region_specs[1]),
    ):
        for period in PERIODS:
            rr = sorted(
                [
                    r for r in stage3_rows
                    if str(r.get("period", "")) == period.key
                    and str(r.get("region", "")) == region
                    and np.isfinite(
                        float(r.get("data_over_mc_scale_factor", np.nan))
                    )
                ],
                key=lambda r: float(r["energy_center_GeV"]),
            )
            if not rr:
                continue
            #endif

            x = np.asarray(
                [float(r["energy_center_GeV"]) for r in rr],
                dtype=float,
            )
            y = np.asarray(
                [float(r["data_over_mc_scale_factor"]) for r in rr],
                dtype=float,
            )
            yerr = np.asarray(
                [
                    float(
                        r.get(
                            "scale_factor_counting_error_provisional",
                            np.nan,
                        )
                    )
                    for r in rr
                ],
                dtype=float,
            )
            ax.errorbar(
                x,
                y,
                yerr=yerr,
                marker="o",
                ms=4.0,
                linewidth=1.1,
                capsize=0.0,
                label=period.label,
            )
        #endfor

        ax.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
        ax.set_ylim(0.0, ymax)
        ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        ax.set_ylabel(
            r"$SF_\gamma=\epsilon_{\mathrm{data}}/\epsilon_{\mathrm{MC}}$"
        )
        ax.set_title(title)
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8, frameon=False)
    #endfor

    fig.suptitle(
        "RGA photon-efficiency reference scale factors by run period",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        Path(outdir) / "canvas_cross_period_stage3_scale_factors.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)


def make_cross_period_stage2_composition_canvas(
    shared_rows: List[Dict[str, object]],
    discriminator_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Cross-period Stage-II composition summary.

    Top: shared morphed f_pi0 versus predicted probe energy.
    Bottom: absolute disagreement between the two production-driver fits.
    """
    if not shared_rows:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)

    driver_map: Dict[Tuple[str, str, float, float, str], Dict[str, object]] = {}
    for row in discriminator_rows:
        disc = str(row.get("discriminator", ""))
        if (
            int(row.get("fit_success", 0)) != 1
            or disc not in STAGE2_DRIVER_DISCRIMINATORS
        ):
            continue
        #endif
        key = (
            str(row.get("period", "")),
            str(row.get("region", "")),
            round(float(row.get("energy_low_GeV", np.nan)), 9),
            round(float(row.get("energy_high_GeV", np.nan)), 9),
            disc,
        )
        driver_map[key] = row
    #endfor

    fig, axes = plt.subplots(2, 2, figsize=(15.8, 10.0), sharex="col")

    for icol, region in enumerate(("FT", "FD_all")):
        for period in PERIODS:
            rr = sorted(
                [
                    r for r in shared_rows
                    if str(r.get("period", "")) == period.key
                    and str(r.get("region", "")) == region
                    and int(r.get("fit_success", 0)) == 1
                    and np.isfinite(float(r.get("pi0_fraction", np.nan)))
                ],
                key=lambda r: float(r["energy_center_GeV"]),
            )
            if not rr:
                continue
            #endif

            x = np.asarray(
                [float(r["energy_center_GeV"]) for r in rr],
                dtype=float,
            )
            axes[0, icol].plot(
                x,
                [float(r["pi0_fraction"]) for r in rr],
                marker="o",
                ms=4.0,
                linewidth=1.1,
                label=period.label,
            )

            dx = []
            dy = []
            for row in rr:
                elo = round(float(row["energy_low_GeV"]), 9)
                ehi = round(float(row["energy_high_GeV"]), 9)
                rp = driver_map.get(
                    (
                        period.key,
                        region,
                        elo,
                        ehi,
                        "mx2_ep_x_probe_m2",
                    )
                )
                rt = driver_map.get(
                    (
                        period.key,
                        region,
                        elo,
                        ehi,
                        "mx2_ep_x_pTmiss",
                    )
                )
                if rp is None or rt is None:
                    continue
                #endif
                fp = float(rp.get("pi0_fraction", np.nan))
                ft = float(rt.get("pi0_fraction", np.nan))
                if np.isfinite(fp) and np.isfinite(ft):
                    dx.append(float(row["energy_center_GeV"]))
                    dy.append(abs(fp - ft))
                #endif
            #endfor

            if dx:
                axes[1, icol].plot(
                    dx,
                    dy,
                    marker="o",
                    ms=4.0,
                    linewidth=1.1,
                    label=period.label,
                )
            #endif
        #endfor

        axes[0, icol].set_ylim(0.0, 1.02)
        axes[0, icol].set_ylabel(r"shared fitted $f_{\pi^0}$")
        axes[0, icol].set_title("FT" if region == "FT" else "FD all")
        axes[0, icol].grid(alpha=0.20)
        axes[0, icol].legend(fontsize=8, frameon=False)

        axes[1, icol].set_ylim(0.0, 0.70 if region == "FT" else 0.25)
        axes[1, icol].set_xlabel(
            r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
        )
        axes[1, icol].set_ylabel(
            r"$|\Delta f_{\pi^0}|$ between production drivers"
        )
        axes[1, icol].grid(alpha=0.20)
    #endfor

    fig.suptitle(
        "RGA Stage-II denominator composition and driver consistency",
        fontsize=14,
    )
    safe_finalize_figure(
        fig,
        Path(outdir) / "canvas_cross_period_stage2_composition.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)




def make_cross_period_closure_summary_canvas(
    closure_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """Worst absolute template-mixture closure bias by period and energy."""
    if not closure_rows:
        return
    #endif

    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.8))

    for ax, region, title, ymax in (
        (axes[0], "FT", "FT", 0.30),
        (axes[1], "FD_all", "FD all", 0.20),
    ):
        for period in PERIODS:
            grouped: Dict[Tuple[float, float], List[float]] = {}
            for row in closure_rows:
                if (
                    str(row.get("period", "")) != period.key
                    or str(row.get("region", "")) != region
                    or int(row.get("fit_success", 0)) != 1
                ):
                    continue
                #endif
                bias = float(row.get("closure_bias", np.nan))
                if not np.isfinite(bias):
                    continue
                #endif
                key = (
                    float(row["energy_low_GeV"]),
                    float(row["energy_high_GeV"]),
                )
                grouped.setdefault(key, []).append(abs(bias))
            #endfor

            if grouped:
                keys = sorted(grouped)
                ax.plot(
                    [0.5 * (a + b) for a, b in keys],
                    [max(grouped[k]) for k in keys],
                    marker="o",
                    linewidth=1.1,
                    label=period.label,
                )
            #endif
        #endfor

        ax.axhline(
            0.01, color="black", linestyle=":", linewidth=1.0, label="1% bias"
        )
        ax.axhline(
            0.05, color="black", linestyle="--", linewidth=0.9,
            alpha=0.7, label="5% bias"
        )
        ax.set_ylim(0.0, ymax)
        ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        ax.set_ylabel(r"max closure $|\Delta f_{\pi^0}|$")
        ax.set_title(title)
        ax.grid(alpha=0.20)
        ax.legend(fontsize=8, frameon=False)
    #endfor

    fig.suptitle("RGA denominator-fit template-mixture closure", fontsize=14)
    safe_finalize_figure(
        fig,
        Path(outdir) / "canvas_summary_template_closure.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)


def make_summary_dashboard(
    stage2_shared_rows: List[Dict[str, object]],
    stage2_discriminator_rows: List[Dict[str, object]],
    closure_rows: List[Dict[str, object]],
    stage3_rows: List[Dict[str, object]],
    ft_three_component_summary_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Create a deliberately small physics dashboard for routine inspection.

    The summary directory contains only high-value cross-period plots:
      1. Stage-III data/MC photon-efficiency scale factors;
      2. Stage-II pi0 composition + driver disagreement;
      3. coarse-FT three-component diagnostic across periods;
      4. template-mixture closure summary.

    No physics calculation is performed here; only already-computed rows are
    visualized.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Reuse the existing cross-period summary makers, but render into summary/.
    if stage3_rows:
        make_cross_period_stage3_scale_factor_canvas(stage3_rows, outdir)
    #endif

    if stage2_shared_rows and stage2_discriminator_rows:
        make_cross_period_stage2_composition_canvas(
            stage2_shared_rows,
            stage2_discriminator_rows,
            outdir,
        )
    #endif

    if closure_rows:
        make_cross_period_closure_summary_canvas(closure_rows, outdir)
    #endif

    if ft_three_component_summary_rows:
        order = {p.key: i for i, p in enumerate(PERIODS)}
        rr = sorted(
            ft_three_component_summary_rows,
            key=lambda r: (
                float(r.get("energy_low_GeV", np.nan)),
                order.get(str(r.get("period", "")), 999),
            ),
        )

        fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.7))

        for period in PERIODS:
            pr = sorted(
                [
                    r for r in rr
                    if str(r.get("period", "")) == period.key
                ],
                key=lambda r: float(r["energy_center_GeV"]),
            )
            if not pr:
                continue
            #endif

            x = np.asarray(
                [float(r["energy_center_GeV"]) for r in pr],
                dtype=float,
            )

            axes[0].plot(
                x,
                [float(r["two_component_pi0_fraction"]) for r in pr],
                marker="o",
                linestyle="--",
                linewidth=1.0,
                alpha=0.70,
            )
            axes[0].plot(
                x,
                [
                    float(r["median_three_component_pi0_fraction"])
                    for r in pr
                ],
                marker="s",
                linewidth=1.2,
                label=period.label,
            )

            axes[1].plot(
                x,
                [
                    float(r["median_three_component_mixed_fraction"])
                    for r in pr
                ],
                marker="o",
                linewidth=1.2,
                label=period.label,
            )

            axes[2].plot(
                x,
                [
                    float(r["median_driver_abs_delta_pi0_fraction"])
                    for r in pr
                ],
                marker="o",
                linewidth=1.2,
                label=period.label,
            )
        #endfor

        axes[0].set_title(r"FT $f_{\pi^0}$: dashed 2-comp, solid 3-comp")
        axes[0].set_ylabel(r"$f_{\pi^0}$")
        axes[0].set_ylim(0.0, 1.0)

        axes[1].set_title("FT mixed/wrong-photon fraction")
        axes[1].set_ylabel(r"$f_{\mathrm{mixed}}$")
        axes[1].set_ylim(0.0, 0.35)

        axes[2].set_title(r"FT driver disagreement after 3-comp fit")
        axes[2].set_ylabel(r"$|\Delta f_{\pi^0}|$")
        axes[2].set_ylim(bottom=0.0)

        for ax in axes:
            ax.set_xlabel(
                r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
            )
            ax.grid(alpha=0.20)
        #endfor
        axes[1].legend(fontsize=8, frameon=False)

        fig.suptitle(
            "RGA coarse-FT three-component diagnostic",
            fontsize=14,
        )
        safe_finalize_figure(
            fig,
            Path(outdir) / "canvas_summary_ft_three_component.png",
            rect=(0, 0, 1, 0.93),
        )
        plt.close(fig)
    #endif


def compact_plot_enabled(args_dict: Dict[str, object]) -> bool:
    return str(args_dict.get("plot_mode", "compact")) == "compact"


def main() -> int:
    args = parse_args()

    if args.max_entries < 0:
        raise ValueError("--max-entries must be >= 0.")
    #endif
    if (
        not args.stage1_only
        and not args.skip_stage3
        and args.max_entries != 0
    ):
        raise ValueError(
            "Stage III requires --max-entries 0 because numerator associations "
            "span independently ordered epgamma/eppi0 files. Use --skip-stage3 "
            "for truncated Stage-I/II development runs."
        )
    #endif
    if not (0.0 < args.tag_min < args.tag_max):
        raise ValueError("Require 0 < --tag-min < --tag-max.")
    #endif
    if args.parent_component_tol <= 0.0:
        raise ValueError("--parent-component-tol must be > 0.")
    #endif
    if args.stage3_tag_remainder_m2_max <= 0.0:
        raise ValueError("--stage3-tag-remainder-m2-max must be > 0.")
    #endif
    if args.stage3_probe_angle_max_deg <= 0.0:
        raise ValueError("--stage3-probe-angle-max-deg must be > 0.")
    #endif
    if args.stage3_probe_frac_energy_max <= 0.0:
        raise ValueError("--stage3-probe-frac-energy-max must be > 0.")
    #endif
    if args.workers < 1:
        raise ValueError("--workers must be >= 1.")
    #endif
    if args.kdtree_workers < 1:
        raise ValueError("--kdtree-workers must be >= 1.")
    #endif
    if args.kdtree_query_chunk < 10000:
        raise ValueError("--kdtree-query-chunk must be >= 10000.")
    #endif
    closure_truth_fractions = parse_float_edges(
        args.stage2_closure_truth_fractions,
        "--stage2-closure-truth-fractions",
    )
    ft_coarse_energy_edges = np.asarray(
        parse_float_edges(
            args.ft_coarse_energy_edges,
            "--ft-coarse-energy-edges",
        ),
        dtype=float,
    )
    if (
        ft_coarse_energy_edges.size < 2
        or np.any(~np.isfinite(ft_coarse_energy_edges))
        or np.any(np.diff(ft_coarse_energy_edges) <= 0.0)
    ):
        raise ValueError(
            "--ft-coarse-energy-edges must contain at least two finite, "
            "strictly increasing values."
        )
    #endif
    if (
        closure_truth_fractions.size == 0
        or np.any(closure_truth_fractions <= 0.0)
        or np.any(closure_truth_fractions >= 1.0)
    ):
        raise ValueError(
            "--stage2-closure-truth-fractions values must lie strictly between 0 and 1."
        )
    #endif
    if not (args.den_fit_mm2_min < args.den_fit_mm2_max):
        raise ValueError("Require --den-fit-mm2-min < --den-fit-mm2-max.")
    #endif
    if args.den_fit_probe_m2_max <= 0.0:
        raise ValueError("--den-fit-probe-m2-max must be > 0.")
    #endif
    if args.den_fit_mm2_bins < 8 or args.den_fit_probe_m2_bins < 8:
        raise ValueError("Stage-II histogram bin counts must each be >= 8.")
    #endif
    if not (0.40 < args.den_probe_energy_max <= 10.0):
        raise ValueError("--den-probe-energy-max must be > 0.40 and <= 10 GeV.")
    #endif
    if args.disc_1d_bins < 8:
        raise ValueError("--disc-1d-bins must be >= 8.")
    #endif
    if args.disc_ptmiss_max <= 0.0 or args.disc_ptmiss_bins < 8:
        raise ValueError(
            "--disc-ptmiss-max must be > 0 and --disc-ptmiss-bins >= 8."
        )
    #endif
    if args.disc_theta_max <= 0.0 or args.disc_theta_bins < 8:
        raise ValueError(
            "--disc-theta-max must be > 0 and --disc-theta-bins >= 8."
        )
    #endif
    parse_mix_shifts(args.mix_shifts)
    if not (args.control_mgg_min < args.control_mgg_max):
        raise ValueError("Require --control-mgg-min < --control-mgg-max.")
    #endif
    if args.control_emiss_abs_max <= 0.0:
        raise ValueError("--control-emiss-abs-max must be > 0.")
    #endif
    parse_float_edges(args.control_pi0_energy_edges, "--control-pi0-energy-edges")

    requested_workers = min(int(args.workers), 8)

    selected = [
        p for p in PERIODS
        if args.period is None or p.key in set(args.period)
    ]

    outroot = Path(args.output_dir)
    outroot.mkdir(parents=True, exist_ok=True)

    stage2_outroot = Path(args.stage2_output_dir)
    stage2_outroot.mkdir(parents=True, exist_ok=True)
    stage3_outroot = Path(args.stage3_output_dir)
    if not args.skip_stage3:
        stage3_outroot.mkdir(parents=True, exist_ok=True)
    #endif

    run_internal_self_tests(outroot)

    if args.aggregate_only:
        aggregate_existing_outputs(
            selected=selected,
            stage2_outroot=stage2_outroot,
            stage3_outroot=stage3_outroot,
            include_stage3=(not args.skip_stage3),
            summary_outroot=outroot / "summary",
        )
        return 0
    #endif

    log(
        "Internal MC tag definition: "
        f"{args.tag_min:g} <= E_tag < {args.tag_max:g} GeV."
    )
    if args.skip_stage3:
        log(
            "Efficiency extraction disabled by --skip-stage3; "
            "denominator-composition study only."
        )
    else:
        log(
            "Efficiency extraction enabled: forming a PRELIMINARY wagon-reference "
            "data/MC photon-efficiency ratio below 2 GeV."
        )
    #endif

    preflight(
        selected,
        include_stage2=True,
        include_stage3=(not args.skip_stage3),
    )

    n_processes = min(requested_workers, len(selected))
    log(
        f"Period-level parallelism: {n_processes} process(es); "
        f"cKDTree workers/process = {args.kdtree_workers}; "
        f"KD query chunk = {args.kdtree_query_chunk:,}."
    )

    provenance = {
        "script": Path(__file__).name,
        "stage": 3 if not args.skip_stage3 else 2,
        "purpose": (
            "internal aaogen MC tag/probe association + real-data denominator composition + "
            "preliminary wagon-reference efficiency scale factor"
        ),
        "arguments": vars(args),
        "effective_period_workers": n_processes,
        "worker_cap": 8,
        "periods": [asdict(p) for p in selected],
        "schema_assumptions": {
            "epgamma": "p1=proton, p2=photon/tag",
            "eppi0": (
                "p1=proton, p2=reconstructed pi0 momentum; Mh_gammagamma is "
                "the event-by-event reconstructed pi0 mass; reconstructed probe "
                "is formed directly as P_pi0_reco - k_tag_reco"
            ),
        },
        "performance_notes": {
            "parallelization": "independent run periods use separate processes",
            "mc_association_output": "internal only; no routine pair-level files are written",
            "detector_region_definition": (
                "FT-like: theta_pred <= 5.5 deg by default; FD sectors use exact "
                "wrapped phi intervals [330,30), [30,90), [90,150), [150,210), "
                "[210,270), [270,330)"
            ),
            "kdtree": (
                "search coordinates use float32; final physical matching cut "
                "uses original float64 momentum components"
            ),
            "nested_parallelism": (
                "cKDTree defaults to one thread/process to avoid oversubscription"
            ),
            "full_statistics_worker_default": (
                "2 period processes to reduce RAM/page-cache and filesystem thrashing"
            ),
            "kdtree_query_chunk": int(args.kdtree_query_chunk),
            "plot_mode": str(args.plot_mode),
            "memory_model": (
                "sample-specific ROOT branch reads; slim sample objects; float32 "
                "persistent Stage-II histogram features; prompt release of large "
                "temporary/sample objects"
            ),
            "stage3_reference_model": (
                "Data numerator uses exact runnum+evnum epgamma-to-eppi0 association; "
                "MC numerator uses the internal aaogen e/p kinematic matcher. Data denominator "
                "is fitted Stage-II pi0 fraction times tag count; aaogen denominator is "
                "the corresponding exclusive-pi0 tag count. Stage-III uncertainties are "
                "provisional counting-only diagnostics pending event-level bootstrap."
            ),
            "stage2_denominator_model": (
                "real epgamma data are fit locally with floated aaogen-pi0 and "
                "dvcsgen BH/DVCS templates under multiple discriminator choices. "
                "No generator relative normalization is imposed. Shared morphing "
                "uses zero-centered weak priors in v031; reconstructed-eppi0 control "
                "results are diagnostic and are not injected as nuisance centers. "
                "Mixed-data wrong-tag templates remain diagnostic-only stress tests. A diagnostic coarse-energy FT fit uses 0.40--1.20 and 1.20--2.00 GeV bins. The existing deterministic mixed-event wrong-photon construction is floated as a third component in this coarse FT diagnostic only, with the two-component pi0/BH-DVCS morph nuisance values held fixed."
            ),
        },
    }

    args_dict = vars(args).copy()
    mc_association_summaries: List[Dict[str, object]] = []
    stage2_summaries: List[Dict[str, object]] = []
    all_stage2_rows: List[Dict[str, object]] = []
    all_stage2_spread_rows: List[Dict[str, object]] = []
    all_mixed_diag_rows: List[Dict[str, object]] = []
    all_shared_rows: List[Dict[str, object]] = []
    all_closure_rows: List[Dict[str, object]] = []
    all_ft_coarse_rows: List[Dict[str, object]] = []
    all_ft_coarse_closure_rows: List[Dict[str, object]] = []
    all_ft_coarse_three_component_rows: List[Dict[str, object]] = []
    all_ft_coarse_three_component_closure_rows: List[Dict[str, object]] = []
    all_ft_coarse_three_component_summary_rows: List[Dict[str, object]] = []
    all_profile_rows: List[Dict[str, object]] = []
    all_control_rows: List[Dict[str, object]] = []
    all_stage3_rows: List[Dict[str, object]] = []
    stage3_summaries: List[Dict[str, object]] = []

    if n_processes == 1:
        for period in selected:
            result = process_period(
                period,
                args_dict,
                str(outroot),
                str(stage2_outroot),
            )
            mc_association_summaries.append(
                result["mc_association_summary"]
            )
            if result["stage2_summary"] is not None:
                stage2_summaries.append(result["stage2_summary"])
                all_stage2_rows.extend(result["stage2_rows"])
                all_stage2_spread_rows.extend(result["stage2_spread_rows"])
                all_mixed_diag_rows.extend(result["mixed_diag_rows"])
                all_profile_rows.extend(result.get("profile_rows", []))
                all_shared_rows.extend(result.get("shared_rows", []))
                all_closure_rows.extend(result.get("closure_rows", []))
                all_ft_coarse_rows.extend(result.get("ft_coarse_rows", []))
                all_ft_coarse_closure_rows.extend(result.get("ft_coarse_closure_rows", []))
                all_ft_coarse_three_component_rows.extend(result.get("ft_coarse_three_component_rows", []))
                all_ft_coarse_three_component_closure_rows.extend(result.get("ft_coarse_three_component_closure_rows", []))
                all_ft_coarse_three_component_summary_rows.extend(result.get("ft_coarse_three_component_summary_rows", []))
                all_control_rows.extend(result.get("control_rows", []))
                all_stage3_rows.extend(result.get("stage3_rows", []))
                if result.get("stage3_summary") is not None:
                    stage3_summaries.append(result["stage3_summary"])
                #endif
            #endif
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            future_to_period = {
                executor.submit(
                    process_period,
                    period,
                    args_dict,
                    str(outroot),
                    str(stage2_outroot),
                ): period
                for period in selected
            }

            for future in as_completed(future_to_period):
                period = future_to_period[future]
                try:
                    result = future.result()
                    summaries.append(result["summary"])
                    all_resolution_rows.extend(result["resolution_rows"])
                    if result["stage2_summary"] is not None:
                        stage2_summaries.append(result["stage2_summary"])
                        all_stage2_rows.extend(result["stage2_rows"])
                        all_stage2_spread_rows.extend(result.get("stage2_spread_rows", []))
                        all_mixed_diag_rows.extend(result.get("mixed_diag_rows", []))
                        all_profile_rows.extend(result.get("profile_rows", []))
                        all_shared_rows.extend(result.get("shared_rows", []))
                        all_closure_rows.extend(result.get("closure_rows", []))
                        all_ft_coarse_rows.extend(result.get("ft_coarse_rows", []))
                        all_ft_coarse_closure_rows.extend(result.get("ft_coarse_closure_rows", []))
                        all_ft_coarse_three_component_rows.extend(result.get("ft_coarse_three_component_rows", []))
                        all_ft_coarse_three_component_closure_rows.extend(result.get("ft_coarse_three_component_closure_rows", []))
                        all_ft_coarse_three_component_summary_rows.extend(result.get("ft_coarse_three_component_summary_rows", []))
                        all_control_rows.extend(result.get("control_rows", []))
                        all_stage3_rows.extend(result.get("stage3_rows", []))
                        if result.get("stage3_summary") is not None:
                            stage3_summaries.append(result["stage3_summary"])
                        #endif
                    #endif
                except Exception as exc:
                    for other in future_to_period:
                        other.cancel()
                    #endfor
                    raise RuntimeError(
                        f"Parallel photon-efficiency processing failed for {period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    order = {p.key: i for i, p in enumerate(selected)}
    # Internal MC association produces no standalone aggregate files.

    if True:
        stage2_summaries.sort(
            key=lambda row: order.get(str(row["period"]), 999)
        )
        write_summary_csv(
            stage2_summaries,
            stage2_outroot / "stage2_summary.csv",
        )

        if all_stage2_rows:
            region_order2 = {
                "all": 0,
                "FT": 1,
                "FD_all": 2,
                **{f"FD_S{s}": s + 2 for s in range(1, 7)},
            }
            all_stage2_rows.sort(
                key=lambda r: (
                    order.get(str(r["period"]), 999),
                    region_order2.get(str(r["region"]), 999),
                    float(r["energy_low_GeV"]),
                )
            )
            write_rows_csv(
                all_stage2_rows,
                stage2_outroot / "denominator_discriminator_fits.csv",
            )
        #endif

        if all_stage2_spread_rows:
            all_stage2_spread_rows.sort(
                key=lambda r: (
                    order.get(str(r["period"]), 999),
                    str(r["region"]),
                    float(r["energy_low_GeV"]),
                )
            )
            write_rows_csv(
                all_stage2_spread_rows,
                stage2_outroot / "denominator_discriminator_spread.csv",
            )
        #endif

        if all_mixed_diag_rows:
            all_mixed_diag_rows.sort(
                key=lambda r: (
                    order.get(str(r["period"]), 999),
                    str(r["region"]),
                    float(r["energy_low_GeV"]),
                    int(r["mix_shift"]),
                )
            )
            write_rows_csv(
                all_mixed_diag_rows,
                stage2_outroot / "mixed_component_diagnostics.csv",
            )
        #endif

        if all_shared_rows:
            all_shared_rows.sort(key=lambda r: (order.get(str(r["period"]),999), str(r["region"]), float(r["energy_low_GeV"])))
            write_rows_csv(all_shared_rows, stage2_outroot / "denominator_shared_morphed_fits.csv")
        #endif

        if all_closure_rows:
            all_closure_rows.sort(
                key=lambda r: (
                    order.get(str(r.get("period", "")), 999),
                    str(r.get("region", "")),
                    float(r.get("energy_low_GeV", "nan")),
                    float(r.get("injected_pi0_fraction", "nan")),
                )
            )
            write_rows_csv(
                all_closure_rows,
                stage2_outroot / "template_mixture_closure.csv",
            )
            make_stage2_template_mixture_closure_canvas(
                all_closure_rows,
                stage2_outroot,
            )
        #endif

        if all_ft_coarse_rows:
            all_ft_coarse_rows.sort(
                key=lambda r: (
                    order.get(str(r.get("period", "")), 999),
                    float(r.get("energy_low_GeV", "nan")),
                )
            )
            write_rows_csv(
                all_ft_coarse_rows,
                stage2_outroot / "ft_coarse_shared_morphed_fits.csv",
            )
        #endif

        if all_ft_coarse_closure_rows:
            all_ft_coarse_closure_rows.sort(
                key=lambda r: (
                    order.get(str(r.get("period", "")), 999),
                    float(r.get("energy_low_GeV", "nan")),
                    float(r.get("injected_pi0_fraction", "nan")),
                )
            )
            write_rows_csv(
                all_ft_coarse_closure_rows,
                stage2_outroot / "ft_coarse_template_mixture_closure.csv",
            )
        #endif


        if all_ft_coarse_three_component_rows:
            write_rows_csv(all_ft_coarse_three_component_rows, stage2_outroot / "ft_coarse_three_component_by_mix_shift.csv")
        #endif
        if all_ft_coarse_three_component_closure_rows:
            write_rows_csv(all_ft_coarse_three_component_closure_rows, stage2_outroot / "ft_coarse_three_component_closure.csv")
        #endif
        if all_ft_coarse_three_component_summary_rows:
            write_rows_csv(all_ft_coarse_three_component_summary_rows, stage2_outroot / "ft_coarse_three_component_summary.csv")
        #endif

        if all_profile_rows:
            all_profile_rows.sort(
                key=lambda r: (
                    order.get(str(r.get("period", "")), 999),
                    float(r.get("energy_low_GeV", "nan")),
                    str(r.get("driver", "")),
                    str(r.get("component", "")),
                    str(r.get("nuisance", "")),
                    float(r.get("fixed_value_bins", "nan")),
                )
            )
            write_rows_csv(
                all_profile_rows,
                stage2_outroot / "morph_nuisance_profiles.csv",
            )
        #endif

        if all_control_rows:
            all_control_rows.sort(key=lambda r: (order.get(str(r.get("period","")),999), str(r.get("selection","")), str(r.get("group",""))))
            write_rows_csv(all_control_rows, stage2_outroot / "pi0_control_validation.csv")
        #endif

        if all_shared_rows and all_stage2_rows:
            make_cross_period_stage2_composition_canvas(
                all_shared_rows,
                all_stage2_rows,
                stage2_outroot,
            )
        #endif

        with (stage2_outroot / "stage2_summary.json").open("w") as f:
            json.dump(
                {
                    "method": (
                        "Local two-component extended-Poisson fits compare several "
                        "independent discriminator choices. aaogen models exclusive-"
                        "pi0 missing-photon tags and dvcsgen models genuine BH/DVCS "
                        "epgamma; both normalizations float independently. The "
                        "shared morphed fit uses zero-centered weak nuisance priors; "
                        "profile-grid initialized alternating morph fits plus explicit all-region "
                        "one-nuisance scans test whether "
                        "nonzero morphing is actually preferred. The reconstructed-eppi0 "
                        "pTmiss control study uses the confirmed 0=FT, 1=FD photon "
                        "detector convention and remains diagnostic-only until its "
                        "data/MC core and tails are validated. Deterministic "
                        "mixed-data tag samples remain diagnostic-only stress tests."
                    ),
                    "period_summaries": stage2_summaries,
                },
                f,
                indent=2,
                allow_nan=True,
            )
        #endwith

    if not args.skip_stage3:
        order = {p.key: i for i, p in enumerate(selected)}
        if all_stage3_rows:
            region_order3 = {
                "all": 0,
                "FT": 1,
                "FD_all": 2,
                **{f"FD_S{s}": s + 2 for s in range(1, 7)},
            }
            all_stage3_rows.sort(
                key=lambda r: (
                    order.get(str(r.get("period", "")), 999),
                    region_order3.get(str(r.get("region", "")), 999),
                    float(r.get("energy_low_GeV", "nan")),
                )
            )
            write_rows_csv(
                all_stage3_rows,
                stage3_outroot / "reference_efficiency_scale_factors.csv",
            )
            make_cross_period_stage3_scale_factor_canvas(
                all_stage3_rows,
                stage3_outroot,
            )
        #endif
        if stage3_summaries:
            stage3_summaries.sort(
                key=lambda r: order.get(str(r.get("period", "")), 999)
            )
            write_summary_csv(
                stage3_summaries,
                stage3_outroot / "stage3_summary.csv",
            )
            summary_outroot = outroot / "summary"
        make_summary_dashboard(
            all_shared_rows,
            all_stage2_rows,
            all_closure_rows,
            all_stage3_rows,
            all_ft_coarse_three_component_summary_rows,
            summary_outroot,
        )

        with (stage3_outroot / "stage3_summary.json").open("w") as f:
                json.dump(
                    {
                        "status": "preliminary wagon-reference extraction",
                        "final_production_result": False,
                        "period_summaries": stage3_summaries,
                        "future_validation": (
                            "The future loose-nSidis data extraction must reproduce "
                            "this reference below 2 GeV before extending above 2 GeV."
                        ),
                    },
                    f,
                    indent=2,
                    allow_nan=True,
                )
            #endwith
        #endif

    log(f"Done. Stage-I outputs are in {outroot}.")
    if True:
        log(f"Stage-II outputs are in {stage2_outroot}.")
    #endif
    if not args.skip_stage3:
        log(f"Stage-III reference outputs are in {stage3_outroot}.")
    #endif
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
