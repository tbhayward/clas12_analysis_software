#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v007.py

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

    python derive_photon_efficiency_scale_factors_v007.py

Run one period:

    python derive_photon_efficiency_scale_factors_v007.py --period fa18_inb

Run all available entries:

    python derive_photon_efficiency_scale_factors_v007.py --max-entries 0

Use up to eight worker processes (default):

    python derive_photon_efficiency_scale_factors_v007.py --max-entries 0 --workers 8

If the ROOT angles are known explicitly:

    python derive_photon_efficiency_scale_factors_v007.py --angles rad

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


# Probe-energy binning used for Stage-Ib resolution studies.
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
    if not rows:
        return
    #endif
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    #endwith
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
    axes[1, 0].set_xlabel(r"$E_{\rm probe}^{reco}$ (GeV)")
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
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / "canvas_clean_tag_association.png", dpi=180)
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
    axes[0, 0].set_xlabel(r"$E_{\rm probe}^{pred}$ (GeV)")
    axes[0, 0].set_ylabel(r"$E_{\rm probe}^{reco}$ (GeV)")
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
    axes[1, 0].set_xlabel(r"$E_{\rm reco}-E_{\rm pred}$ (GeV)")
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
    axes[1, 1].set_xlabel(r"$E_{\rm probe}^{pred}$ (GeV)")
    axes[1, 1].set_ylabel(r"$\Delta\alpha_{\rm probe}$ (deg)")
    axes[1, 1].set_title("Angular resolution vs predicted energy")

    fig.suptitle(f"{title}: clean-sample probe prediction resolution", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / "canvas_clean_probe_resolution.png", dpi=180)
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
        ax.set_xlabel(r"$E_{\rm probe}^{pred}$ (GeV)")
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
    axes[1, 1].set_xlabel(r"$E_{\rm probe}^{pred}$ (GeV)")
    axes[1, 1].set_ylabel("Clean entries")
    axes[1, 1].set_title("Statistics by predicted energy")
    axes[1, 1].grid(alpha=0.18)

    fig.suptitle(
        f"{title}: candidate angular matching containment by predicted region",
        fontsize=15,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / "canvas_matching_resolution_quantiles.png", dpi=180)
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
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(outdir / "canvas_exclusivity_after_clean_association.png", dpi=180)
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
    parser.add_argument(
        "--save-npz",
        action="store_true",
        help=(
            "Optionally save the large pair-level stage1_matched_pairs.npz file. "
            "Disabled by default."
        ),
    )
    parser.add_argument(
        "--assoc-mgg-min",
        type=float,
        default=0.10,
        help="Lower clean-association M_gammagamma edge (GeV). Default: 0.10.",
    )
    parser.add_argument(
        "--assoc-mgg-max",
        type=float,
        default=0.17,
        help="Upper clean-association M_gammagamma edge (GeV). Default: 0.17.",
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
            "Predicted-probe polar-angle boundary (deg) used for the Stage-Ib "
            "FT-like versus FD-like resolution diagnostic. FT-like is defined "
            "as theta_pred <= this value. Default: 5.5."
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

    res_rows = resolution_rows(
        period,
        pair_np,
        clean,
        ft_theta_max=float(args_dict["ft_theta_max"]),
        min_count=int(args_dict["min_resolution_count"]),
    )
    write_rows_csv(res_rows, period_dir / "probe_resolution_by_energy_region.csv")

    if bool(args_dict.get("save_npz", False)):
        np.savez_compressed(
            period_dir / "stage1_matched_pairs.npz",
            clean_association=clean,
            **pair_np,
        )
    else:
        old_npz = period_dir / "stage1_matched_pairs.npz"
        if old_npz.exists():
            old_npz.unlink()
        #endif
    #endif

    make_plots(
        period,
        pair_np,
        clean,
        res_rows,
        period_dir,
        assoc_mgg_min=float(args_dict["assoc_mgg_min"]),
        assoc_mgg_max=float(args_dict["assoc_mgg_max"]),
        assoc_mass2_max=float(args_dict["assoc_remainder_mass2_max"]),
    )

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
        clean=clean,
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

    return {"summary": summary, "resolution_rows": res_rows}

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
        "No data/MC photon-efficiency ratio is formed in v006. "
        "This run defines a reconstructed-side clean pi0 tag association and "
        "measures probe-prediction resolution versus energy and predicted region."
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
        "purpose": "aaogen clean tag association and probe-prediction resolution",
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
            "pair_level_npz": "disabled by default; enable with --save-npz",
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
        },
    }

    args_dict = vars(args).copy()
    summaries: List[Dict[str, object]] = []
    all_resolution_rows: List[Dict[str, object]] = []

    if n_processes == 1:
        for period in selected:
            result = process_period(period, args_dict, str(outroot))
            summaries.append(result["summary"])
            all_resolution_rows.extend(result["resolution_rows"])
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
                    result = future.result()
                    summaries.append(result["summary"])
                    all_resolution_rows.extend(result["resolution_rows"])
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
    if all_resolution_rows:
        region_order = {"all": 0, "FT_like": 1, **{f"FD_S{s}": s + 1 for s in range(1, 7)}}
        period_order = {p.key: i for i, p in enumerate(selected)}
        all_resolution_rows.sort(
            key=lambda r: (
                period_order.get(str(r["period"]), 999),
                region_order.get(str(r["region"]), 999),
                float(r["energy_low_GeV"]),
            )
        )
        write_rows_csv(
            all_resolution_rows,
            outroot / "probe_resolution_by_energy_region.csv",
        )
    #endif

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
