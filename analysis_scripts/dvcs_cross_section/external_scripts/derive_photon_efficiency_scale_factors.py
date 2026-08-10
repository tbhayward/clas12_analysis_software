#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v004.py

Stage-I development script for a relative data/MC photon-reconstruction
efficiency measurement in CLAS12 RGA.

THIS VERSION DOES NOT EXTRACT A PHOTON EFFICIENCY OR SCALE FACTOR.

Its purpose is to validate the central tag-and-probe kinematics using the
exclusive-pi0 aaogen samples:

    aaogen reconstructed as epgamma:
        .../bkg_rga_<period>_epgamma_0.40GeV.root

    aaogen reconstructed as eppi0:
        .../rec_aaogen_norad_<period>_...root

Because runnum/evnum are not useful identifiers in the MC, the two files are
matched through the reconstructed electron and proton Cartesian momenta.

For every matched epgamma -> eppi0 pair:

  1. The epgamma photon is treated as the directed tag.
  2. The missing/probe photon is predicted from four-momentum conservation,

         p_probe^pred = p_beam + p_target - p_e - p_p - p_tag .

  3. The reconstructed pi0 four-vector is formed directly from the eppi0 tree:
       * p2_p, p2_theta, p2_phi provide the reconstructed pi0 momentum;
       * Mh_gammagamma provides the event-by-event reconstructed pi0 mass;
       * E_pi0 = sqrt(|p_pi0|^2 + Mh_gammagamma^2).

  4. The reconstructed companion photon is then obtained without any Trento
     inversion:

         p_probe^reco = p_pi0^reco - p_tag^reco .

     If the epgamma tag is genuinely one daughter of the reconstructed pi0,
     this remainder should itself be photon-like:
         E_probe^reco > 0,
         (p_probe^reco)^2 ~ 0,
         E_probe^reco - |p_probe^reco| ~ 0.

  5. Predicted-vs-reconstructed probe residuals and full exclusive eppi0
     closure diagnostics are written and plotted.

The stored gamma_phi1/2 and open_angle_egamma1/2 branches are retained only as
optional future diagnostics; they are NOT used to construct the Stage-I probe.

Important schema assumptions in v001
------------------------------------
epgamma tree:
    p1_* = proton
    p2_* = reconstructed photon/tag

eppi0 tree:
    p1_* = proton
    p2_* = reconstructed pi0 momentum
    Mh_gammagamma = event-by-event reconstructed gamma-gamma mass
    gamma_phi1/2 and open_angle_egamma1/2 are auxiliary only in v004.

These assumptions are checked with diagnostic closures, not silently treated
as established facts.  If the actual producer uses another mapping, change
the adapters in `extract_epgamma()` / `extract_eppi0()` rather than changing
the physics calculations downstream.

Outputs
-------
All products go below:

    output/photon_efficiency/stage1/

including:
    stage1_summary.csv
    stage1_summary.json
    <period>/stage1_matched_pairs.npz
    <period>/*.png

Dependencies
------------
    numpy
    scipy
    uproot
    matplotlib

Typical usage
-------------
Quick orientation run over the first 500k entries of each relevant file:

    python derive_photon_efficiency_scale_factors_v004.py

Run one period:

    python derive_photon_efficiency_scale_factors_v004.py --period fa18_inb

Run all available entries:

    python derive_photon_efficiency_scale_factors_v004.py --max-entries 0

Use up to eight worker processes (default):

    python derive_photon_efficiency_scale_factors_v004.py --max-entries 0 --workers 8

If the ROOT angles are known explicitly:

    python derive_photon_efficiency_scale_factors_v004.py --angles rad

The Stage-I defaults intentionally require a reconstructed tag energy
E_tag >= 2 GeV, while no efficiency denominator is formed yet.
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
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from scipy.spatial import cKDTree


# -----------------------------------------------------------------------------
# Physics constants
# -----------------------------------------------------------------------------

M_E = 0.00051099895000
M_P = 0.93827208816
M_PI0 = 0.1349768
TWO_PI = 2.0 * math.pi


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
        "sp19_inb", "Sp19 Inb", 10.200,
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

EPG_OPTIONAL = (
    "runnum", "evnum", "detector1", "detector2",
    "Mx2", "Mx2_1", "Mx2_2", "Emiss2", "pTmiss",
)

EPPIO_REQUIRED = (
    "e_p", "e_theta", "e_phi",
    "p1_p", "p1_theta", "p1_phi",
    "p2_p", "p2_theta", "p2_phi",
    "Mh_gammagamma",
)

EPPIO_OPTIONAL = (
    "runnum", "evnum", "detector_gamma1", "detector_gamma2",
    "gamma_phi1", "gamma_phi2",
    "open_angle_egamma1", "open_angle_egamma2",
    "Mx2", "Mx2_1", "Mx2_2",
    "Emiss2", "pTmiss", "fiducial_status",
)


@dataclass
class EPGammaSample:
    electron_p3: np.ndarray
    proton_p3: np.ndarray
    tag_p3: np.ndarray
    electron_p: np.ndarray
    proton_p: np.ndarray
    tag_energy: np.ndarray
    electron_theta: np.ndarray
    electron_phi: np.ndarray
    proton_theta: np.ndarray
    proton_phi: np.ndarray
    tag_theta: np.ndarray
    tag_phi: np.ndarray
    raw: Dict[str, np.ndarray]
    angle_unit: str


@dataclass
class EPPi0Sample:
    electron_p3: np.ndarray
    proton_p3: np.ndarray
    pi0_p3: np.ndarray
    electron_p: np.ndarray
    proton_p: np.ndarray
    pi0_p: np.ndarray
    pi0_mass: np.ndarray
    electron_theta: np.ndarray
    electron_phi: np.ndarray
    proton_theta: np.ndarray
    proton_phi: np.ndarray
    pi0_theta: np.ndarray
    pi0_phi: np.ndarray
    raw: Dict[str, np.ndarray]
    angle_unit: str


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
        electron_p=e_p,
        proton_p=p_p,
        tag_energy=g_e,
        electron_theta=e_theta,
        electron_phi=e_phi,
        proton_theta=p_theta,
        proton_phi=p_phi,
        tag_theta=g_theta,
        tag_phi=g_phi,
        raw=arrays,
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
        electron_p=e_p,
        proton_p=p_p,
        pi0_p=pi_p,
        pi0_mass=pi_mass,
        electron_theta=e_theta,
        electron_phi=e_phi,
        proton_theta=p_theta,
        proton_phi=p_phi,
        pi0_theta=pi_theta,
        pi0_phi=pi_phi,
        raw=arrays,
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
) -> MatchResult:
    """
    Match reconstructed parent (electron, proton) states.

    Features are the six Cartesian momentum components:
        e_px,e_py,e_pz,p_px,p_py,p_pz.

    Dividing by component_tolerance makes the KD-tree metric dimensionless.
    A match must pass BOTH:
        Euclidean scaled distance <= nearest_distance_max
        max absolute Cartesian component difference <= component_tolerance

    The second condition prevents a moderately bad match in one coordinate
    from being hidden by the six-dimensional Euclidean metric.
    """
    x_pi0 = np.column_stack((eppi0.electron_p3, eppi0.proton_p3))
    x_epg = np.column_stack((epg.electron_p3, epg.proton_p3))

    good_pi0 = np.all(np.isfinite(x_pi0), axis=1)
    good_epg = (
        np.all(np.isfinite(x_epg), axis=1)
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

    # Use float32 search coordinates to reduce KD-tree memory and bandwidth.
    # The final physical component-difference cut is still evaluated with the
    # original float64 momenta below.
    pi0_search = np.asarray(
        x_pi0[pi0_indices] / component_tolerance, dtype=np.float32
    )
    epg_search = np.asarray(
        x_epg[epg_indices] / component_tolerance, dtype=np.float32
    )

    tree = cKDTree(pi0_search, compact_nodes=True, balanced_tree=True)
    distance, local_index = tree.query(
        epg_search,
        k=1,
        workers=max(1, int(kdtree_workers)),
    )

    candidate_pi0 = pi0_indices[local_index]
    deltas = np.abs(x_epg[epg_indices] - x_pi0[candidate_pi0])
    max_delta = np.max(deltas, axis=1)

    accepted = (
        np.isfinite(distance)
        & (distance <= nearest_distance_max)
        & (max_delta <= component_tolerance)
    )

    return MatchResult(
        epg_index=epg_indices[accepted],
        eppi0_index=candidate_pi0[accepted],
        nearest_distance=distance[accepted],
        max_component_delta=max_delta[accepted],
    )


# -----------------------------------------------------------------------------
# Stage-I pair construction
# -----------------------------------------------------------------------------

@dataclass
class PairArrays:
    epg_index: List[int]
    eppi0_index: List[int]
    parent_distance: List[float]
    parent_max_component_delta: List[float]

    tag_energy: List[float]

    pi0_mass: List[float]
    pi0_energy: List[float]

    reco_probe_energy: List[float]
    reco_probe_p: List[float]
    reco_probe_mass2: List[float]
    reco_probe_E_minus_p: List[float]

    pred_probe_energy: List[float]
    pred_probe_p: List[float]
    pred_probe_mass2: List[float]
    pred_probe_E_minus_p: List[float]

    probe_delta_E: List[float]
    probe_delta_E_over_E: List[float]
    probe_delta_theta_deg: List[float]
    probe_delta_phi_deg: List[float]
    probe_opening_residual_deg: List[float]

    exclusivity_missing_energy: List[float]
    exclusivity_missing_p: List[float]
    exclusivity_missing_pT: List[float]
    exclusivity_missing_mass2: List[float]

    detector_tag_epgamma: List[int]

    @classmethod
    def empty(cls) -> "PairArrays":
        return cls(*([[] for _ in cls.__annotations__]))


def event_four_vector(momentum3: np.ndarray, mass: float) -> np.ndarray:
    p2 = float(np.dot(momentum3, momentum3))
    E = math.sqrt(max(0.0, p2 + mass * mass))
    return np.asarray([E, momentum3[0], momentum3[1], momentum3[2]], dtype=float)


def photon_four_vector(momentum3: np.ndarray) -> np.ndarray:
    E = float(np.linalg.norm(momentum3))
    return np.asarray([E, momentum3[0], momentum3[1], momentum3[2]], dtype=float)


def detector_value(raw: Dict[str, np.ndarray], branch: str, index: int) -> int:
    if branch not in raw:
        return -999
    #endif
    try:
        return int(raw[branch][index])
    except Exception:
        return -999


def build_stage1_pairs(
    period: PeriodConfig,
    epg: EPGammaSample,
    eppi0: EPPi0Sample,
    matches: MatchResult,
) -> Tuple[PairArrays, Dict[str, int]]:
    """
    Build direct reconstructed-probe and missing-probe quantities.

    The reconstructed companion photon is

        probe_reco = pi0_reco - tag_reco,

    where the reconstructed pi0 energy uses the measured Mh_gammagamma.

    No angular matching cut is imposed in Stage I.  We want the complete
    distribution, including failures, so that association quality is visible
    rather than hidden by an upstream cut.
    """
    out = PairArrays.empty()
    counters = {
        "matched_parent_pairs": int(len(matches.epg_index)),
        "invalid_pi0_mass": 0,
        "nonphysical_reco_probe_energy": 0,
        "nonphysical_predicted_probe": 0,
        "accepted_stage1_pairs": 0,
    }

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    beam4 = np.asarray([period.beam_energy, 0.0, 0.0, beam_p], dtype=float)
    target4 = np.asarray([M_P, 0.0, 0.0, 0.0], dtype=float)

    for k in range(len(matches.epg_index)):
        i = int(matches.epg_index[k])
        j = int(matches.eppi0_index[k])

        pi0_mass = float(eppi0.pi0_mass[j])
        if not math.isfinite(pi0_mass) or pi0_mass <= 0.0:
            counters["invalid_pi0_mass"] += 1
            continue
        #endif

        electron4 = event_four_vector(epg.electron_p3[i], M_E)
        proton4 = event_four_vector(epg.proton_p3[i], M_P)
        tag4 = photon_four_vector(epg.tag_p3[i])

        pi0_p3 = np.asarray(eppi0.pi0_p3[j], dtype=float)
        pi0_p = float(np.linalg.norm(pi0_p3))
        pi0_energy = math.sqrt(max(0.0, pi0_p * pi0_p + pi0_mass * pi0_mass))
        pi04 = np.asarray(
            [pi0_energy, pi0_p3[0], pi0_p3[1], pi0_p3[2]],
            dtype=float,
        )

        # Reconstructed companion photon candidate from reconstructed pi0.
        reco_probe4 = pi04 - tag4
        Ereco = float(reco_probe4[0])
        preco3 = reco_probe4[1:4]
        preco = float(np.linalg.norm(preco3))
        m2reco = Ereco * Ereco - preco * preco

        if not np.all(np.isfinite(reco_probe4)):
            counters["nonphysical_reco_probe_energy"] += 1
            continue
        #endif
        if Ereco <= 0.0:
            counters["nonphysical_reco_probe_energy"] += 1
        #endif

        # Probe predicted only from beam/target/e/p/tag.
        pred_probe4 = beam4 + target4 - electron4 - proton4 - tag4
        Epred = float(pred_probe4[0])
        ppred3 = pred_probe4[1:4]
        ppred = float(np.linalg.norm(ppred3))
        m2pred = Epred * Epred - ppred * ppred

        if not np.all(np.isfinite(pred_probe4)) or ppred <= 0.0:
            counters["nonphysical_predicted_probe"] += 1
            continue
        #endif

        # Full ep -> ep pi0 exclusivity closure.
        miss4 = beam4 + target4 - electron4 - proton4 - pi04
        Emiss = float(miss4[0])
        pmiss3 = miss4[1:4]
        pmiss = float(np.linalg.norm(pmiss3))
        pTmiss = float(math.hypot(pmiss3[0], pmiss3[1]))
        Mmiss2 = Emiss * Emiss - pmiss * pmiss

        # Direction residuals are meaningful only when the reconstructed
        # remainder has nonzero momentum.
        if preco > 0.0:
            npred = ppred3 / ppred
            nreco = preco3 / preco

            theta_pred = math.acos(max(-1.0, min(1.0, float(npred[2]))))
            phi_pred = math.atan2(float(npred[1]), float(npred[0]))
            theta_reco = math.acos(max(-1.0, min(1.0, float(nreco[2]))))
            phi_reco = math.atan2(float(nreco[1]), float(nreco[0]))

            dtheta = math.degrees(theta_reco - theta_pred)
            dphi = math.degrees(float(wrap_phi(phi_reco - phi_pred)))
            opening = math.degrees(angle_between(npred, nreco))
        else:
            dtheta = float("nan")
            dphi = float("nan")
            opening = float("nan")
        #endif

        out.epg_index.append(i)
        out.eppi0_index.append(j)
        out.parent_distance.append(float(matches.nearest_distance[k]))
        out.parent_max_component_delta.append(float(matches.max_component_delta[k]))

        out.tag_energy.append(float(epg.tag_energy[i]))

        out.pi0_mass.append(pi0_mass)
        out.pi0_energy.append(pi0_energy)

        out.reco_probe_energy.append(Ereco)
        out.reco_probe_p.append(preco)
        out.reco_probe_mass2.append(m2reco)
        out.reco_probe_E_minus_p.append(Ereco - preco)

        out.pred_probe_energy.append(Epred)
        out.pred_probe_p.append(ppred)
        out.pred_probe_mass2.append(m2pred)
        out.pred_probe_E_minus_p.append(Epred - ppred)

        out.probe_delta_E.append(Ereco - Epred)
        out.probe_delta_E_over_E.append(
            (Ereco - Epred) / Epred if Epred != 0.0 else float("nan")
        )
        out.probe_delta_theta_deg.append(dtheta)
        out.probe_delta_phi_deg.append(dphi)
        out.probe_opening_residual_deg.append(opening)

        out.exclusivity_missing_energy.append(Emiss)
        out.exclusivity_missing_p.append(pmiss)
        out.exclusivity_missing_pT.append(pTmiss)
        out.exclusivity_missing_mass2.append(Mmiss2)

        out.detector_tag_epgamma.append(
            detector_value(epg.raw, "detector2", i)
        )

        counters["accepted_stage1_pairs"] += 1
    #endfor

    return out, counters


def pair_arrays_to_numpy(pairs: PairArrays) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    integer_names = {"epg_index", "eppi0_index", "detector_tag_epgamma"}
    for name in pairs.__annotations__:
        dtype = np.int64 if name in integer_names else float
        out[name] = np.asarray(getattr(pairs, name), dtype=dtype)
    #endfor
    return out


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
    fig.tight_layout()
    fig.savefig(path, dpi=180)
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
    fig.tight_layout()
    fig.savefig(path, dpi=180)
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


def make_plots(period: PeriodConfig, arrays: Dict[str, np.ndarray], outdir: Path) -> None:
    if len(arrays["pred_probe_energy"]) == 0:
        log(f"{period.label}: no accepted pairs; skipping plots.")
        return
    #endif

    title = period.label

    save_hist(
        arrays["parent_max_component_delta"],
        outdir / "parent_match_max_cartesian_delta.png",
        r"max $|\Delta p_{e/p,x/y/z}|$ (GeV)",
        f"{title}: e/p kinematic match quality",
        bins=120,
        logy=True,
    )

    # Reconstructed pi0 mass directly from Mh_gammagamma.
    save_hist(
        arrays["pi0_mass"],
        outdir / "pi0_Mh_gammagamma.png",
        r"$M_{\gamma\gamma}$ (GeV)",
        f"{title}: reconstructed $M_{{\\gamma\\gamma}}$ used for $P_{{\\pi^0}}$",
        bins=160,
        xlim=(0.0, 0.30),
    )

    # Photon-likeness of P_pi0 - k_tag.
    reco_m2_zoom = symmetric_zoom(arrays["reco_probe_mass2"], 99.0, 0.02, 1.0)
    make_full_and_zoom_hist(
        arrays["reco_probe_mass2"],
        outdir,
        "reco_probe_mass2",
        r"$(P_{\pi^0}^{reco}-k_{tag}^{reco})^2$ (GeV$^2$)",
        f"{title}: reconstructed companion photon mass shell",
        (-reco_m2_zoom, reco_m2_zoom),
    )

    reco_ep_zoom = symmetric_zoom(arrays["reco_probe_E_minus_p"], 99.0, 0.02, 1.0)
    make_full_and_zoom_hist(
        arrays["reco_probe_E_minus_p"],
        outdir,
        "reco_probe_E_minus_p",
        r"$E_{\rm probe}^{reco}-|\vec p_{\rm probe}^{reco}|$ (GeV)",
        f"{title}: reconstructed companion photon mass-shell residual",
        (-reco_ep_zoom, reco_ep_zoom),
    )

    # Predicted missing photon mass shell.
    pred_m2_zoom = symmetric_zoom(arrays["pred_probe_mass2"], 99.0, 0.02, 1.0)
    make_full_and_zoom_hist(
        arrays["pred_probe_mass2"],
        outdir,
        "pred_probe_mass2",
        r"$(p_{\rm probe}^{pred})^2$ (GeV$^2$)",
        f"{title}: predicted missing-photon mass shell",
        (-pred_m2_zoom, pred_m2_zoom),
    )

    pred_ep_zoom = symmetric_zoom(arrays["pred_probe_E_minus_p"], 99.0, 0.02, 1.0)
    make_full_and_zoom_hist(
        arrays["pred_probe_E_minus_p"],
        outdir,
        "pred_probe_E_minus_p",
        r"$E_{\rm probe}^{pred}-|\vec p_{\rm probe}^{pred}|$ (GeV)",
        f"{title}: predicted photon mass-shell residual",
        (-pred_ep_zoom, pred_ep_zoom),
    )

    # Probe energy comparison. Full physically broad range plus a log-z copy.
    emax = max(
        10.0,
        finite_percentile(
            np.concatenate((arrays["pred_probe_energy"], arrays["reco_probe_energy"])),
            99.9,
            10.0,
        ),
    )
    emax = min(emax, 15.0)

    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["reco_probe_energy"],
        outdir / "probe_energy_pred_vs_reco_full_linear.png",
        r"$E_{\rm probe}^{pred}$ (GeV)",
        r"$E_{\rm probe}^{reco}$ (GeV)",
        f"{title}: predicted vs reconstructed probe energy",
        bins=(140, 140),
        xlim=(-1.0, emax),
        ylim=(-1.0, emax),
        logz=False,
    )
    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["reco_probe_energy"],
        outdir / "probe_energy_pred_vs_reco_full_logz.png",
        r"$E_{\rm probe}^{pred}$ (GeV)",
        r"$E_{\rm probe}^{reco}$ (GeV)",
        f"{title}: predicted vs reconstructed probe energy",
        bins=(140, 140),
        xlim=(-1.0, emax),
        ylim=(-1.0, emax),
        logz=True,
    )

    de_zoom = symmetric_zoom(arrays["probe_delta_E"], 99.0, 0.10, 3.0)
    make_full_and_zoom_hist(
        arrays["probe_delta_E"],
        outdir,
        "probe_delta_E",
        r"$E_{\rm probe}^{reco}-E_{\rm probe}^{pred}$ (GeV)",
        f"{title}: probe energy residual",
        (-de_zoom, de_zoom),
    )

    rel_zoom = symmetric_zoom(arrays["probe_delta_E_over_E"], 99.0, 0.10, 3.0)
    make_full_and_zoom_hist(
        arrays["probe_delta_E_over_E"],
        outdir,
        "probe_fractional_delta_E",
        r"$(E_{\rm reco}-E_{\rm pred})/E_{\rm pred}$",
        f"{title}: fractional probe energy residual",
        (-rel_zoom, rel_zoom),
    )

    # Full angular range is always shown; zoom is additional.
    save_hist(
        arrays["probe_opening_residual_deg"],
        outdir / "probe_direction_opening_residual_full.png",
        r"$\Delta\alpha(\gamma_{\rm probe}^{pred},\gamma_{\rm probe}^{reco})$ (deg)",
        f"{title}: probe direction residual (full range)",
        bins=180,
        xlim=(0.0, 180.0),
        logy=True,
    )
    opening_zoom = max(
        1.0,
        min(20.0, finite_percentile(arrays["probe_opening_residual_deg"], 95.0, 5.0)),
    )
    save_hist(
        arrays["probe_opening_residual_deg"],
        outdir / "probe_direction_opening_residual_zoom.png",
        r"$\Delta\alpha(\gamma_{\rm probe}^{pred},\gamma_{\rm probe}^{reco})$ (deg)",
        f"{title}: probe direction residual (core zoom)",
        bins=160,
        xlim=(0.0, opening_zoom),
        logy=False,
    )

    dtheta_zoom = symmetric_zoom(arrays["probe_delta_theta_deg"], 99.0, 1.0, 30.0)
    make_full_and_zoom_hist(
        arrays["probe_delta_theta_deg"],
        outdir,
        "probe_delta_theta",
        r"$\theta_{\rm reco}-\theta_{\rm pred}$ (deg)",
        f"{title}: probe polar-angle residual",
        (-dtheta_zoom, dtheta_zoom),
    )

    dphi_zoom = symmetric_zoom(arrays["probe_delta_phi_deg"], 99.0, 1.0, 60.0)
    make_full_and_zoom_hist(
        arrays["probe_delta_phi_deg"],
        outdir,
        "probe_delta_phi",
        r"$\phi_{\rm reco}-\phi_{\rm pred}$ (deg)",
        f"{title}: probe azimuth residual",
        (-dphi_zoom, dphi_zoom),
    )

    # Full exclusive eppi0 closure.
    emiss_zoom = symmetric_zoom(arrays["exclusivity_missing_energy"], 99.0, 0.05, 2.0)
    make_full_and_zoom_hist(
        arrays["exclusivity_missing_energy"],
        outdir,
        "eppi0_exclusivity_missing_energy",
        r"$E_{\rm miss}(ep\pi^0)$ (GeV)",
        f"{title}: full reconstructed $ep\\pi^0$ energy closure",
        (-emiss_zoom, emiss_zoom),
    )

    save_hist(
        arrays["exclusivity_missing_p"],
        outdir / "eppi0_exclusivity_missing_p_full.png",
        r"$|\vec p_{\rm miss}(ep\pi^0)|$ (GeV)",
        f"{title}: full reconstructed $ep\\pi^0$ momentum closure",
        bins=160,
        logy=True,
    )
    pmiss_zoom = max(
        0.05,
        min(2.0, finite_percentile(arrays["exclusivity_missing_p"], 99.0, 0.5)),
    )
    save_hist(
        arrays["exclusivity_missing_p"],
        outdir / "eppi0_exclusivity_missing_p_zoom.png",
        r"$|\vec p_{\rm miss}(ep\pi^0)|$ (GeV)",
        f"{title}: full reconstructed $ep\\pi^0$ momentum closure (core zoom)",
        bins=160,
        xlim=(0.0, pmiss_zoom),
    )

    save_hist(
        arrays["exclusivity_missing_pT"],
        outdir / "eppi0_exclusivity_missing_pT.png",
        r"$p_{T,\rm miss}(ep\pi^0)$ (GeV)",
        f"{title}: reconstructed $ep\\pi^0$ transverse closure",
        bins=160,
        logy=True,
    )

    mm2_zoom = symmetric_zoom(arrays["exclusivity_missing_mass2"], 99.0, 0.02, 2.0)
    make_full_and_zoom_hist(
        arrays["exclusivity_missing_mass2"],
        outdir,
        "eppi0_exclusivity_missing_mass2",
        r"$M_X^2(ep\pi^0)$ (GeV$^2$)",
        f"{title}: reconstructed $ep\\pi^0$ missing mass squared",
        (-mm2_zoom, mm2_zoom),
    )

    # Residuals versus predicted probe energy; log-z makes sparse structures visible.
    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["probe_opening_residual_deg"],
        outdir / "probe_direction_residual_vs_energy_logz.png",
        r"$E_{\rm probe}^{pred}$ (GeV)",
        r"$\Delta\alpha_{\rm probe}$ (deg)",
        f"{title}: probe direction residual vs predicted energy",
        bins=(120, 120),
        xlim=(-1.0, emax),
        ylim=(0.0, 180.0),
        logz=True,
    )

    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["probe_delta_E_over_E"],
        outdir / "probe_fractional_energy_residual_vs_energy_logz.png",
        r"$E_{\rm probe}^{pred}$ (GeV)",
        r"$(E_{\rm reco}-E_{\rm pred})/E_{\rm pred}$",
        f"{title}: probe energy closure vs predicted energy",
        bins=(120, 120),
        xlim=(-1.0, emax),
        ylim=(-3.0, 3.0),
        logz=True,
    )

    save_hist2d(
        arrays["reco_probe_mass2"],
        arrays["probe_opening_residual_deg"],
        outdir / "tag_association_mass2_vs_direction_residual_logz.png",
        r"$(P_{\pi^0}^{reco}-k_{tag}^{reco})^2$ (GeV$^2$)",
        r"$\Delta\alpha_{\rm probe}$ (deg)",
        f"{title}: tag-association quality",
        bins=(120, 120),
        xlim=(-1.0, 1.0),
        ylim=(0.0, 180.0),
        logz=True,
    )


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


# -----------------------------------------------------------------------------
# Preflight and driver
# -----------------------------------------------------------------------------

def preflight(periods: Sequence[PeriodConfig]) -> None:
    log("Stage-I preflight: checking aaogen epgamma and aaogen eppi0 files.")
    missing: List[str] = []

    for p in periods:
        for label, path in (
            ("epgamma pi0 MC", p.epgamma_pi0_mc),
            ("eppi0 pi0 MC", p.eppi0_pi0_mc),
        ):
            if not Path(path).exists():
                missing.append(f"{p.label}: {label}: {path}")
            else:
                log(f"Preflight OK: {p.label}: {label}")
            #endif
        #endfor
    #endfor

    if missing:
        raise FileNotFoundError(
            "The following Stage-I input files do not exist:\n  " + "\n  ".join(missing)
        )
    #endif


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stage-I aaogen tag-and-probe kinematic validation."
    )
    parser.add_argument(
        "--period",
        action="append",
        choices=[p.key for p in PERIODS],
        help="Run only this period. May be supplied more than once. Default: all periods.",
    )
    parser.add_argument(
        "--output-dir",
        default="output/photon_efficiency/stage1",
        help="Stage-I output directory.",
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
        default=8,
        help=(
            "Maximum number of independent period worker processes. Hard-capped "
            "at 8. With all RGA periods there are only five independent jobs, so "
            "at most five processes run concurrently. Default: 8."
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
    return parser.parse_args()



def process_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
    output_dir: str,
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
    period_dir = outroot / period.key
    period_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.perf_counter()

    log(f"{period.label}: reading aaogen epgamma sample.")
    epg_arrays, epg_tree, epg_total = read_branches(
        period.epgamma_pi0_mc,
        EPG_REQUIRED,
        EPG_OPTIONAL,
        tree_name,
        max_entries,
    )
    epg = extract_epgamma(epg_arrays, angle_mode)

    log(
        f"{period.label}: epgamma tree '{epg_tree}', "
        f"loaded {len(epg.tag_energy):,}/{epg_total:,} entries, "
        f"angles interpreted as {epg.angle_unit}."
    )

    log(f"{period.label}: reading reconstructed aaogen eppi0 sample.")
    pi_arrays, pi_tree, pi_total = read_branches(
        period.eppi0_pi0_mc,
        EPPIO_REQUIRED,
        EPPIO_OPTIONAL,
        tree_name,
        max_entries,
    )
    eppi0 = extract_eppi0(pi_arrays, angle_mode)

    log(
        f"{period.label}: eppi0 tree '{pi_tree}', "
        f"loaded {len(eppi0.pi0_p):,}/{pi_total:,} entries, "
        f"angles interpreted as {eppi0.angle_unit}."
    )

    log(f"{period.label}: matching e/p parent kinematics between MC skims.")
    matches = match_parent_kinematics(
        epg,
        eppi0,
        tag_min=tag_min,
        tag_max=tag_max,
        component_tolerance=parent_component_tol,
        nearest_distance_max=parent_distance_max,
        kdtree_workers=kdtree_workers,
    )
    log(f"{period.label}: accepted {len(matches.epg_index):,} parent matches.")

    log(f"{period.label}: forming reconstructed pi0-minus-tag remainder and validating predicted probe.")
    pairs, counters = build_stage1_pairs(
        period,
        epg,
        eppi0,
        matches,
    )
    pair_np = pair_arrays_to_numpy(pairs)

    np.savez_compressed(period_dir / "stage1_matched_pairs.npz", **pair_np)
    make_plots(period, pair_np, period_dir)

    summary = summarize_period(
        period,
        epg_entries_total=epg_total,
        eppi0_entries_total=pi_total,
        epg_loaded=len(epg.tag_energy),
        eppi0_loaded=len(eppi0.pi0_p),
        epg_angle_unit=epg.angle_unit,
        eppi0_angle_unit=eppi0.angle_unit,
        matches=matches,
        arrays=pair_np,
        counters=counters,
    )
    summary["wall_time_s"] = float(time.perf_counter() - t0)

    with (period_dir / "stage1_summary.json").open("w") as f:
        json.dump(summary, f, indent=2, allow_nan=True)
    #endif

    log(
        f"{period.label}: Stage-I retained pairs = "
        f"{counters['accepted_stage1_pairs']:,}; "
        f"invalid Mh_gammagamma = {counters['invalid_pi0_mass']:,}; "
        f"nonpositive reconstructed remainder energy = "
        f"{counters['nonphysical_reco_probe_energy']:,}; "
        f"wall time = {summary['wall_time_s']:.1f} s."
    )

    return summary

def main() -> int:
    args = parse_args()

    if args.max_entries < 0:
        raise ValueError("--max-entries must be >= 0.")
    #endif
    if not (0.0 < args.tag_min < args.tag_max):
        raise ValueError("Require 0 < --tag-min < --tag-max.")
    #endif
    if args.parent_component_tol <= 0.0:
        raise ValueError("--parent-component-tol must be > 0.")
    #endif
    if args.workers < 1:
        raise ValueError("--workers must be >= 1.")
    #endif
    if args.kdtree_workers < 1:
        raise ValueError("--kdtree-workers must be >= 1.")
    #endif

    requested_workers = min(int(args.workers), 8)

    selected = [
        p for p in PERIODS
        if args.period is None or p.key in set(args.period)
    ]

    outroot = Path(args.output_dir)
    outroot.mkdir(parents=True, exist_ok=True)

    log(
        "Directed-tag Stage-I definition: "
        f"{args.tag_min:g} <= E_tag < {args.tag_max:g} GeV."
    )
    log(
        "No photon-efficiency denominator/numerator is formed in v004; "
        "this run validates exact e/p event matching, direct tag association, "
        "the reconstructed companion-photon remainder, and missing-photon prediction."
    )

    preflight(selected)

    n_processes = min(requested_workers, len(selected))
    log(
        f"Period-level parallelism: {n_processes} process(es); "
        f"cKDTree workers/process = {args.kdtree_workers}."
    )

    provenance = {
        "script": Path(__file__).name,
        "stage": 1,
        "purpose": "aaogen tag-and-probe kinematic validation",
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
            "kdtree": (
                "search coordinates use float32; final physical matching cut "
                "uses original float64 momentum components"
            ),
            "nested_parallelism": (
                "cKDTree defaults to one thread/process to avoid oversubscription"
            ),
        },
    }

    args_dict = vars(args).copy()
    summaries: List[Dict[str, object]] = []

    if n_processes == 1:
        for period in selected:
            summaries.append(process_period(period, args_dict, str(outroot)))
        #endfor
    else:
        with ProcessPoolExecutor(max_workers=n_processes) as executor:
            future_to_period = {
                executor.submit(
                    process_period,
                    period,
                    args_dict,
                    str(outroot),
                ): period
                for period in selected
            }

            for future in as_completed(future_to_period):
                period = future_to_period[future]
                try:
                    summaries.append(future.result())
                except Exception as exc:
                    for other in future_to_period:
                        other.cancel()
                    #endfor
                    raise RuntimeError(
                        f"Parallel Stage-I processing failed for {period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    order = {p.key: i for i, p in enumerate(selected)}
    summaries.sort(key=lambda row: order[str(row["period"])])

    write_summary_csv(summaries, outroot / "stage1_summary.csv")

    with (outroot / "stage1_summary.json").open("w") as f:
        json.dump(
            {
                "provenance": provenance,
                "period_summaries": summaries,
            },
            f,
            indent=2,
            allow_nan=True,
        )
    #endif

    log(f"Done. Stage-I outputs are in {outroot}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
