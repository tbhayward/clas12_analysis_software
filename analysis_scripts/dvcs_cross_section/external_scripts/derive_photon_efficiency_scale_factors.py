#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v017.py

Stage-I + Stage-II development script for a relative data/MC photon-reconstruction
efficiency measurement in CLAS12 RGA.

THIS VERSION DOES NOT YET FORM THE FINAL DATA/MC EFFICIENCY SCALE FACTOR.

Stage I validates the exclusive-pi0 tag-and-probe kinematics and derives
probe-prediction resolution. Stage II studies the real-data denominator using local two-component
aaogen-pi0 + dvcsgen BH/DVCS fits under several alternative discriminator
choices. The spread in the extracted pi0 fraction is used to test robustness
against template-shape mismodeling. Mixed-data wrong-tag samples remain
diagnostic stress tests, not nominal fit components.

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
    <period>/stage1_matched_pairs.npz   (only with --save-npz)
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

    python derive_photon_efficiency_scale_factors_v017.py

Run one period:

    python derive_photon_efficiency_scale_factors_v017.py --period fa18_inb

Run all available entries:

    python derive_photon_efficiency_scale_factors_v017.py --max-entries 0

Use up to eight worker processes (default):

    python derive_photon_efficiency_scale_factors_v017.py --max-entries 0 --workers 8

If the ROOT angles are known explicitly:

    python derive_photon_efficiency_scale_factors_v017.py --angles rad

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
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter1d, shift as ndimage_shift


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

EPG_OPTIONAL = (
    "runnum", "evnum", "detector1", "detector2",
    "Mx2", "Mx2_1", "Mx2_2", "Emiss2", "pTmiss",
    "theta_gamma_gamma",
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
    Vectorized Stage-II kinematics for an epgamma sample.

    The physical variables are recomputed from e, p, and the tag photon using
    the period-specific beam energy.

    Stored pTmiss / theta_gamma_gamma are preserved only for the unshifted
    physical samples. They are used as alternative discriminator diagnostics,
    not as nominal cuts.
    """
    e3 = np.asarray(epg.electron_p3, dtype=float)
    p3 = np.asarray(epg.proton_p3, dtype=float)
    tag3 = np.asarray(epg.tag_p3, dtype=float)

    if mixed_tag_shift:
        n = len(tag3)
        if n:
            shift = int(mixed_tag_shift) % n
            if shift == 0:
                shift = max(1, n // 2)
            #endif
            tag3 = np.roll(tag3, shift=shift, axis=0)
        #endif
    #endif

    e_pmag = np.linalg.norm(e3, axis=1)
    p_pmag = np.linalg.norm(p3, axis=1)
    tag_E = np.linalg.norm(tag3, axis=1)

    e_E = np.sqrt(e_pmag * e_pmag + M_E * M_E)
    p_E = np.sqrt(p_pmag * p_pmag + M_P * M_P)

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    initial_E = period.beam_energy + M_P
    initial_p3 = np.asarray([0.0, 0.0, beam_p], dtype=float)

    epmiss_E = initial_E - e_E - p_E
    epmiss3 = initial_p3[None, :] - e3 - p3
    epmiss_p = np.linalg.norm(epmiss3, axis=1)
    epmiss_m2 = epmiss_E * epmiss_E - epmiss_p * epmiss_p

    pred_E = epmiss_E - tag_E
    pred3 = epmiss3 - tag3
    pred_p = np.linalg.norm(pred3, axis=1)
    pred_m2 = pred_E * pred_E - pred_p * pred_p
    pred_E_minus_p = pred_E - pred_p
    pred_pT = np.hypot(pred3[:, 0], pred3[:, 1])

    with np.errstate(invalid="ignore", divide="ignore"):
        pred_unit = pred3 / pred_p[:, None]
        pred_theta_deg = np.degrees(
            np.arccos(np.clip(pred_unit[:, 2], -1.0, 1.0))
        )
        pred_phi_deg = np.degrees(np.arctan2(pred3[:, 1], pred3[:, 0]))
    #endwith

    sector, region = predicted_region_arrays(
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
        "tag_energy": tag_E,
        "ep_missing_mass2": epmiss_m2,
        "ep_missing_energy": epmiss_E,
        "pred_probe_energy": pred_E,
        "pred_probe_p": pred_p,
        "pred_probe_pT": pred_pT,
        "pred_probe_mass2": pred_m2,
        "pred_probe_E_minus_p": pred_E_minus_p,
        "pred_probe_theta_deg": pred_theta_deg,
        "pred_probe_phi_deg": pred_phi_deg,
        "pred_probe_sector": sector,
        "pred_region_code": region,
        "valid_tag": valid,
    }

    if not mixed_tag_shift:
        for raw_name, out_name in (
            ("Emiss2", "stored_Emiss2"),
            ("pTmiss", "stored_pTmiss"),
            ("theta_gamma_gamma", "stored_theta_gamma_gamma"),
        ):
            if raw_name in epg.raw:
                out[out_name] = np.asarray(epg.raw[raw_name], dtype=float)
            #endif
        #endfor
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
    regions = ["all", "FT"] + [f"FD_S{s}" for s in range(1, 7)]
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


def build_pi0_control_calibration(eppi0_data: EPPi0Sample, eppi0_mc: EPPi0Sample,
                                  ptmiss_max: float, ptmiss_bins: int) -> Dict[str, Dict[str, float]]:
    """Use reconstructed ep-pi0 data/aaogen to anchor the pi0 pTmiss morph nuisance."""
    out: Dict[str, Dict[str, float]] = {}
    if "pTmiss" in eppi0_data.raw and "pTmiss" in eppi0_mc.raw:
        out["pTmiss"] = fit_control_shift_smear(
            np.asarray(eppi0_data.raw["pTmiss"], dtype=float),
            np.asarray(eppi0_mc.raw["pTmiss"], dtype=float),
            0.0, ptmiss_max, ptmiss_bins,
        )
    else:
        out["pTmiss"] = {"success": 0, "shift_bins": 0.0, "sigma_bins": 0.0, "deviance_per_ndof": float("nan")}
    #endif
    return out


def fit_shared_morphed_composition(histograms: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]],
                                   pi0_control: Dict[str, Dict[str, float]],
                                   nuisance_shift_prior: float,
                                   nuisance_sigma_prior: float,
                                   max_shift_bins: float,
                                   max_sigma_bins: float) -> SharedMorphedFitResult:
    """
    Composite profile-likelihood fit with one shared pi0 fraction. Each 2D driver
    gets constrained shift/smear nuisances along its second axis for aaogen and
    dvcsgen. pTmiss aaogen is anchored to the clean reconstructed-pi0 control.
    """
    names = [n for n in STAGE2_DRIVER_DISCRIMINATORS if n in histograms]
    if not names:
        return SharedMorphedFitResult(False, "no available driver discriminators")
    #endif
    # x = [logit(f), then pi_shift,pi_sigma,dv_shift,dv_sigma for each driver]
    x0 = [math.log(0.7 / 0.3)]
    bounds = [(-8.0, 8.0)]
    centers: List[Tuple[float, float, float, float]] = []
    for name in names:
        ps, pg = 0.0, 0.0
        if name == "mx2_ep_x_pTmiss" and int(pi0_control.get("pTmiss", {}).get("success", 0)):
            ps = float(pi0_control["pTmiss"]["shift_bins"])
            pg = float(pi0_control["pTmiss"]["sigma_bins"])
        #endif
        centers.append((ps, pg, 0.0, 0.0))
        x0.extend([ps, pg, 0.0, 0.0])
        bounds.extend([(-max_shift_bins, max_shift_bins), (0.0, max_sigma_bins),
                       (-max_shift_bins, max_shift_bins), (0.0, max_sigma_bins)])
    #endfor

    def unpack_fraction(z: float) -> float:
        return 1.0 / (1.0 + math.exp(-float(z)))
    #enddef

    def objective(x: np.ndarray) -> float:
        f = unpack_fraction(x[0])
        total_nll = 0.0
        penalty = 0.0
        k = 1
        for idx, name in enumerate(names):
            hd, hp, hv = histograms[name]
            ps, pg, ds, dg = map(float, x[k:k+4]); k += 4
            tp = morph_template_second_axis(hp, ps, pg).ravel()
            td = morph_template_second_axis(hv, ds, dg).ravel()
            data = np.asarray(hd, dtype=float).ravel()
            nd = float(np.sum(data))
            mu = np.clip(nd * (f * tp + (1.0 - f) * td), 1.0e-12, None)
            # Average the projection likelihoods because they reuse the same events.
            total_nll += float(np.sum(mu - data * np.log(mu))) / len(names)
            c = centers[idx]
            penalty += 0.5 * ((ps-c[0])/nuisance_shift_prior)**2
            penalty += 0.5 * ((pg-c[1])/nuisance_sigma_prior)**2
            penalty += 0.5 * ((ds-c[2])/nuisance_shift_prior)**2
            penalty += 0.5 * ((dg-c[3])/nuisance_sigma_prior)**2
        #endfor
        return total_nll + penalty
    #enddef

    res = minimize(objective, np.asarray(x0, dtype=float), method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": 700, "ftol": 1.0e-10, "maxls": 40})
    f = unpack_fraction(res.x[0])
    nuisance: Dict[str, float] = {}
    total_dev, total_ndof = 0.0, 0
    k = 1
    for name in names:
        hd, hp, hv = histograms[name]
        ps, pg, ds, dg = map(float, res.x[k:k+4]); k += 4
        nuisance[f"{name}_pi0_shift_bins"] = ps
        nuisance[f"{name}_pi0_sigma_bins"] = pg
        nuisance[f"{name}_dvcs_shift_bins"] = ds
        nuisance[f"{name}_dvcs_sigma_bins"] = dg
        tp = morph_template_second_axis(hp, ps, pg).ravel()
        td = morph_template_second_axis(hv, ds, dg).ravel()
        data = np.asarray(hd, dtype=float).ravel(); nd = float(np.sum(data))
        dev, ndof = poisson_deviance_quality(data, nd*(f*tp+(1-f)*td), 1+4*len(names))
        nuisance[f"{name}_deviance_per_ndof"] = float(dev/ndof)
        total_dev += dev; total_ndof += ndof
    #endfor
    return SharedMorphedFitResult(bool(res.success), str(res.message), f, float(res.fun),
                                  float(total_dev), int(total_ndof), nuisance)


def run_stage2_shared_morphed_fits(period: PeriodConfig, data_f: Dict[str,np.ndarray], pi0_f: Dict[str,np.ndarray],
                                   dvcs_f: Dict[str,np.ndarray], pi0_control: Dict[str,Dict[str,float]],
                                   ft_theta_max: float, max_probe_energy: float, mm2_min: float, mm2_max: float,
                                   probe_m2_max: float, mm2_bins_2d: int, probe_m2_bins_2d: int,
                                   ptmiss_max: float, ptmiss_bins: int, theta_max: float, theta_bins: int,
                                   min_data_count: int, min_template_count: int, nuisance_shift_prior: float,
                                   nuisance_sigma_prior: float, max_shift_bins: float, max_sigma_bins: float) -> List[Dict[str,object]]:
    edges=stage2_energy_edges(max_probe_energy); regions=["all","FT"]+[f"FD_S{s}" for s in range(1,7)]
    rows=[]
    for region in regions:
        for ib in range(len(edges)-1):
            elo,ehi=float(edges[ib]),float(edges[ib+1])
            masks={k:stage2_fit_mask(f,region,ft_theta_max,elo,ehi,mm2_min,mm2_max,probe_m2_max)
                   for k,f in (("data",data_f),("pi0",pi0_f),("dvcs",dvcs_f))}
            row={"period":period.key,"label":period.label,"beam_energy_GeV":period.beam_energy,"region":region,
                 "energy_low_GeV":elo,"energy_high_GeV":ehi,"energy_center_GeV":0.5*(elo+ehi),
                 "fit_model":"shared_morphed_2d"}
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
            # Deterministic template-mixture closure in the high-statistics all-region sample.
            # This catches optimizer/boundary pathologies without multiplying runtime in every sector.
            if region == "all" and hist:
                for truth in (0.2, 0.5, 0.8):
                    pseudo = {}
                    for disc, (hd0, hp0, hv0) in hist.items():
                        nd0 = max(1.0, float(np.sum(hd0)))
                        tp0 = normalized_template(hp0).reshape(hp0.shape)
                        td0 = normalized_template(hv0).reshape(hv0.shape)
                        pseudo[disc] = (nd0 * (truth * tp0 + (1.0-truth) * td0), hp0, hv0)
                    #endfor
                    cf = fit_shared_morphed_composition(pseudo, pi0_control, nuisance_shift_prior, nuisance_sigma_prior, max_shift_bins, max_sigma_bins)
                    tag = str(int(round(100*truth)))
                    row[f"closure_truth_{tag}"] = truth
                    row[f"closure_fit_{tag}"] = cf.pi0_fraction
                    row[f"closure_bias_{tag}"] = cf.pi0_fraction-truth if np.isfinite(cf.pi0_fraction) else float("nan")
                #endfor
            #endif
            rows.append(row)
        #endfor
    #endfor
    return rows



def make_shared_fit_canvas(period: PeriodConfig, shared_rows: List[Dict[str,object]], raw_rows: List[Dict[str,object]],
                           pi0_control: Dict[str,Dict[str,float]], outdir: Path) -> None:
    """Compact production-fit summary: shared fraction, raw-driver comparison, fit quality, nuisances."""
    rr=[r for r in shared_rows if r["region"]=="all" and int(r.get("fit_success",0))]
    if not rr: return
    rr=sorted(rr,key=lambda r:float(r["energy_center_GeV"])); x=np.asarray([r["energy_center_GeV"] for r in rr])
    fig,ax=plt.subplots(2,2,figsize=(13.5,9.5))
    ax[0,0].plot(x,[r["pi0_fraction"] for r in rr],marker="o",label="shared morphed fit")
    for disc,label in (("mx2_ep_x_probe_m2", r"raw $M_X^2\otimes M_{\mathrm{probe}}^2$"),("mx2_ep_x_pTmiss", r"raw $M_X^2\otimes p_{T,\mathrm{miss}}$")):
        q=sorted([r for r in raw_rows if r["region"]=="all" and r["discriminator"]==disc and int(r["fit_success"])],key=lambda r:float(r["energy_center_GeV"]))
        if q: ax[0,0].plot([r["energy_center_GeV"] for r in q],[r["pi0_fraction"] for r in q],marker=".",linestyle="--",label=label)
    #endfor
    ax[0,0].set_ylabel(r"$f_{\pi^0}$"); ax[0,0].set_xlabel(r"$E_{probe}^{pred}$ (GeV)"); ax[0,0].set_ylim(0,1); ax[0,0].legend(fontsize=8); ax[0,0].grid(alpha=.18)
    ax[0,1].plot(x,[r["deviance_per_ndof"] for r in rr],marker="o")
    ax[0,1].set_ylabel("combined Poisson deviance / ndof"); ax[0,1].set_xlabel(r"$E_{probe}^{pred}$ (GeV)"); ax[0,1].grid(alpha=.18)
    for key, label in (
        ("mx2_ep_x_probe_m2_pi0_shift_bins", "probe-M2 pi0 shift"),
        ("mx2_ep_x_probe_m2_dvcs_shift_bins", "probe-M2 DVCS shift"),
        ("mx2_ep_x_pTmiss_pi0_shift_bins", "pTmiss pi0 shift"),
        ("mx2_ep_x_pTmiss_dvcs_shift_bins", "pTmiss DVCS shift"),
    ):
        values = np.asarray([r.get(key, np.nan) for r in rr], dtype=float)
        if np.any(np.isfinite(values)):
            ax[1, 0].plot(x, values, marker=".", label=label)
        #endif
    #endfor
    ax[1,0].axhline(0,linewidth=.8); ax[1,0].set_ylabel("template shift (bins)"); ax[1,0].set_xlabel(r"$E_{probe}^{pred}$ (GeV)"); ax[1,0].legend(fontsize=7); ax[1,0].grid(alpha=.18)
    for truth in (20, 50, 80):
        key = f"closure_bias_{truth}"
        values = np.asarray([r.get(key, np.nan) for r in rr], dtype=float)
        if np.any(np.isfinite(values)):
            ax[1, 1].plot(
                x,
                values,
                marker=".",
                label=f"truth {truth/100:.1f}",
            )
        #endif
    #endfor
    ax[1,1].axhline(0,linewidth=.8); ax[1,1].set_ylabel(r"closure $f_{fit}-f_{true}$"); ax[1,1].set_xlabel(r"$E_{probe}^{pred}$ (GeV)"); ax[1,1].legend(fontsize=8); ax[1,1].grid(alpha=.18)
    ctl=pi0_control.get("pTmiss",{}); fig.suptitle(f"{period.label}: shared morphed Stage-II composition; pi0 pTmiss control shift={ctl.get('shift_bins',float('nan')):.2f} bins, smear={ctl.get('sigma_bins',float('nan')):.2f} bins",fontsize=13)
    safe_finalize_figure(
        fig,
        outdir / "canvas_shared_morphed_composition.png",
        rect=(0, 0, 1, 0.95),
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
    for old in outdir.glob("*.png"):
        old.unlink()
    #endfor

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
            np.linspace(-probe_m2_max, probe_m2_max, 220),
            label,
        )
    #endfor

    axes[0, 0].axvline(M_PI0 * M_PI0, linestyle="--", linewidth=1.0)
    axes[0, 0].set_xlabel(r"$M_X^2(ep)$ (GeV$^2$)")
    axes[0, 0].set_ylabel("Unit-normalized entries")
    axes[0, 0].set_title("$ep$ missing mass squared")
    axes[0, 0].legend(fontsize=8)
    axes[0, 0].grid(alpha=0.18)

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
        hi = max(5.0, min(180.0, float(np.percentile(vals_all, 99.5)))) if vals_all.size else 40.0
        for key, feat, label in (
            ("data", data_f, "Data"),
            ("pi0", pi0_f, "aaogen $\\pi^0$"),
            ("dvcs", dvcs_f, "dvcsgen BH/DVCS"),
        ):
            unit_hist(
                axes[1, 1],
                feat["stored_theta_gamma_gamma"][common[key]],
                np.linspace(0.0, hi, 220),
                label,
            )
        #endfor
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
    stage1_outroot: Path,
    stage2_outroot: Path,
    include_stage2: bool,
) -> None:
    """Rebuild aggregate files from existing per-period outputs only."""
    order = {p.key: i for i, p in enumerate(selected)}

    stage1_summaries: List[Dict[str, object]] = []
    resolution_rows: List[Dict[str, object]] = []
    stage2_summaries: List[Dict[str, object]] = []
    disc_rows: List[Dict[str, object]] = []
    spread_rows: List[Dict[str, object]] = []
    mixed_rows: List[Dict[str, object]] = []
    shared_rows: List[Dict[str, object]] = []
    missing: List[str] = []

    for period in selected:
        p1 = Path(stage1_outroot) / period.key
        s1 = p1 / "stage1_summary.json"
        if s1.exists():
            with s1.open() as f:
                stage1_summaries.append(json.load(f))
            #endwith
        else:
            missing.append(str(s1))
        #endif

        resolution_rows.extend(
            read_csv_rows(p1 / "probe_resolution_by_energy_region.csv")
        )

        if include_stage2:
            p2 = Path(stage2_outroot) / period.key
            s2 = p2 / "stage2_summary.json"
            if s2.exists():
                with s2.open() as f:
                    stage2_summaries.append(json.load(f))
                #endwith
            else:
                missing.append(str(s2))
            #endif

            disc_rows.extend(
                read_csv_rows(p2 / "denominator_discriminator_fits.csv")
            )
            spread_rows.extend(
                read_csv_rows(p2 / "denominator_discriminator_spread.csv")
            )
            mixed_rows.extend(
                read_csv_rows(p2 / "mixed_component_diagnostics.csv")
            )
            shared_rows.extend(
                read_csv_rows(p2 / "denominator_shared_morphed_fits.csv")
            )
        #endif
    #endfor

    if missing:
        raise FileNotFoundError(
            "--aggregate-only is missing required per-period files:\n  "
            + "\n  ".join(missing)
        )
    #endif

    stage1_summaries.sort(
        key=lambda r: order.get(str(r.get("period", "")), 999)
    )
    write_summary_csv(
        stage1_summaries,
        Path(stage1_outroot) / "stage1_summary.csv",
    )

    if resolution_rows:
        resolution_rows.sort(
            key=lambda r: (
                order.get(str(r.get("period", "")), 999),
                str(r.get("region", "")),
                float(r.get("energy_low_GeV", "nan")),
            )
        )
        write_rows_csv(
            resolution_rows,
            Path(stage1_outroot) / "probe_resolution_by_energy_region.csv",
        )
    #endif

    with (Path(stage1_outroot) / "stage1_summary.json").open("w") as f:
        json.dump(
            {
                "aggregation_mode": "aggregate-only",
                "period_summaries": stage1_summaries,
            },
            f,
            indent=2,
            allow_nan=True,
        )
    #endwith

    if include_stage2:
        stage2_summaries.sort(
            key=lambda r: order.get(str(r.get("period", "")), 999)
        )
        write_summary_csv(
            stage2_summaries,
            Path(stage2_outroot) / "stage2_summary.csv",
        )

        def _sort(rows: List[Dict[str, object]]) -> None:
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

        for rows, filename in (
            (disc_rows, "denominator_discriminator_fits.csv"),
            (spread_rows, "denominator_discriminator_spread.csv"),
            (mixed_rows, "mixed_component_diagnostics.csv"),
            (shared_rows, "denominator_shared_morphed_fits.csv"),
        ):
            if rows:
                _sort(rows)
                write_rows_csv(rows, Path(stage2_outroot) / filename)
            #endif
        #endfor

        with (Path(stage2_outroot) / "stage2_summary.json").open("w") as f:
            json.dump(
                {
                    "aggregation_mode": "aggregate-only",
                    "period_summaries": stage2_summaries,
                },
                f,
                indent=2,
                allow_nan=True,
            )
        #endwith
    #endif

    log("Aggregate-only recovery completed successfully.")


# -----------------------------------------------------------------------------
# Preflight and driver
# -----------------------------------------------------------------------------

def preflight(
    periods: Sequence[PeriodConfig],
    include_stage2: bool,
) -> None:
    log("Preflight: checking required photon-efficiency input files.")
    missing: List[str] = []

    for p in periods:
        checks = [
            ("epgamma pi0 MC", p.epgamma_pi0_mc),
            ("eppi0 pi0 MC", p.eppi0_pi0_mc),
            ("eppi0 data control", p.eppi0_data),
        ]
        if include_stage2:
            checks.extend([
                ("epgamma data", p.epgamma_data),
                ("epgamma BH/DVCS MC", p.epgamma_dvcs_mc),
            ])
        #endif

        for label, path in checks:
            if not Path(path).exists():
                missing.append(f"{p.label}: {label}: {path}")
            else:
                log(f"Preflight OK: {p.label}: {label}")
            #endif
        #endfor
    #endfor

    if missing:
        raise FileNotFoundError(
            "The following required input files do not exist:\n  " + "\n  ".join(missing)
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

    # Stage-II denominator-composition controls.
    parser.add_argument(
        "--stage1-only",
        action="store_true",
        help=(
            "Run only the existing Stage-I/Stage-Ib validation. By default v008 "
            "also runs the Stage-II epgamma data-composition study."
        ),
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
        help="Upper M_X^2(ep) fit edge (GeV^2). Default: 0.10.",
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
    parser.add_argument("--morph-shift-prior-bins", type=float, default=1.0, help="Gaussian prior width for template shifts in histogram bins. Default: 1.0.")
    parser.add_argument("--morph-sigma-prior-bins", type=float, default=1.0, help="Gaussian prior width for additional template smearing in histogram bins. Default: 1.0.")
    parser.add_argument("--morph-max-shift-bins", type=float, default=4.0, help="Hard bound on template shifts in histogram bins. Default: 4.0.")
    parser.add_argument("--morph-max-sigma-bins", type=float, default=4.0, help="Hard bound on additional template smearing in histogram bins. Default: 4.0.")
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

    log(f"{period.label}: reading reconstructed eppi0 data control sample.")
    pi_data_arrays, pi_data_tree, pi_data_total = read_branches(
        period.eppi0_data, EPPIO_REQUIRED, EPPIO_OPTIONAL, tree_name, max_entries
    )
    eppi0_data = extract_eppi0(pi_data_arrays, angle_mode)
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
    log(
        f"{period.label}: Stage-I retained pairs = "
        f"{counters['accepted_stage1_pairs']:,}; "
        f"invalid Mh_gammagamma = {counters['invalid_pi0_mass']:,}; "
        f"nonpositive reconstructed remainder energy = "
        f"{counters['nonphysical_reco_probe_energy']:,}."
    )

    stage2_summary: Optional[Dict[str, object]] = None
    stage2_rows: List[Dict[str, object]] = []
    stage2_spread_rows: List[Dict[str, object]] = []
    mixed_diag_rows: List[Dict[str, object]] = []
    shared_rows: List[Dict[str, object]] = []

    if not bool(args_dict.get("stage1_only", False)):
        stage2_dir = Path(stage2_output_dir) / period.key
        stage2_dir.mkdir(parents=True, exist_ok=True)

        log(f"{period.label}: Stage II reading real epgamma data.")
        data_arrays, data_tree, data_total = read_branches(
            period.epgamma_data,
            EPG_REQUIRED,
            EPG_OPTIONAL,
            tree_name,
            max_entries,
        )
        data_epg = extract_epgamma(data_arrays, angle_mode)
        log(
            f"{period.label}: data epgamma tree '{data_tree}', "
            f"loaded {len(data_epg.tag_energy):,}/{data_total:,} entries."
        )

        log(f"{period.label}: Stage II reading dvcsgen epgamma MC.")
        dv_arrays, dv_tree, dv_total = read_branches(
            period.epgamma_dvcs_mc,
            EPG_REQUIRED,
            EPG_OPTIONAL,
            tree_name,
            max_entries,
        )
        dv_epg = extract_epgamma(dv_arrays, angle_mode)
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

        pi0_control = build_pi0_control_calibration(
            eppi0_data, eppi0, float(args_dict["disc_ptmiss_max"]), int(args_dict["disc_ptmiss_bins"])
        )
        with (stage2_dir / "pi0_control_calibration.json").open("w") as f:
            json.dump(pi0_control, f, indent=2, allow_nan=True)
        #endwith
        shared_rows = run_stage2_shared_morphed_fits(
            period, data_f, pi0_f, dvcs_f, pi0_control, ft_theta_max, max_probe_energy,
            float(args_dict["den_fit_mm2_min"]), float(args_dict["den_fit_mm2_max"]), float(args_dict["den_fit_probe_m2_max"]),
            int(args_dict["den_fit_mm2_bins"]), int(args_dict["den_fit_probe_m2_bins"]), float(args_dict["disc_ptmiss_max"]),
            int(args_dict["disc_ptmiss_bins"]), float(args_dict["disc_theta_max"]), int(args_dict["disc_theta_bins"]),
            int(args_dict["den_min_data_count"]), int(args_dict["den_min_template_count"]),
            float(args_dict["morph_shift_prior_bins"]), float(args_dict["morph_sigma_prior_bins"]),
            float(args_dict["morph_max_shift_bins"]), float(args_dict["morph_max_sigma_bins"]),
        )
        write_rows_csv(shared_rows, stage2_dir / "denominator_shared_morphed_fits.csv")
        log(
            f"{period.label}: wrote per-period shared-morphed numerical results "
            f"({len(shared_rows)} rows) before canvas rendering."
        )
        make_shared_fit_canvas(period, shared_rows, stage2_rows, pi0_control, stage2_dir)

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
            del mixed_f
        #endfor

        write_rows_csv(
            mixed_diag_rows,
            stage2_dir / "mixed_component_diagnostics.csv",
        )

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
    #endif

    summary["wall_time_s"] = float(time.perf_counter() - t0)
    with (period_dir / "stage1_summary.json").open("w") as f:
        json.dump(summary, f, indent=2, allow_nan=True)
    #endwith

    log(f"{period.label}: total worker wall time = {summary['wall_time_s']:.1f} s.")

    return {
        "summary": summary,
        "resolution_rows": res_rows,
        "stage2_summary": stage2_summary,
        "stage2_rows": stage2_rows,
        "stage2_spread_rows": stage2_spread_rows,
        "mixed_diag_rows": mixed_diag_rows,
        "shared_rows": shared_rows,
    }

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

    requested_workers = min(int(args.workers), 8)

    selected = [
        p for p in PERIODS
        if args.period is None or p.key in set(args.period)
    ]

    outroot = Path(args.output_dir)
    outroot.mkdir(parents=True, exist_ok=True)

    stage2_outroot = Path(args.stage2_output_dir)
    if not args.stage1_only:
        stage2_outroot.mkdir(parents=True, exist_ok=True)
    #endif

    run_internal_self_tests(outroot)

    if args.aggregate_only:
        aggregate_existing_outputs(
            selected=selected,
            stage1_outroot=outroot,
            stage2_outroot=stage2_outroot,
            include_stage2=not args.stage1_only,
        )
        return 0
    #endif

    log(
        "Directed-tag Stage-I definition: "
        f"{args.tag_min:g} <= E_tag < {args.tag_max:g} GeV."
    )
    log(
        "No data/MC photon-efficiency ratio is formed in v006. "
        "This run defines a reconstructed-side clean pi0 tag association and "
        "measures probe-prediction resolution versus energy and predicted region."
    )

    preflight(selected, include_stage2=not args.stage1_only)

    n_processes = min(requested_workers, len(selected))
    log(
        f"Period-level parallelism: {n_processes} process(es); "
        f"cKDTree workers/process = {args.kdtree_workers}."
    )

    provenance = {
        "script": Path(__file__).name,
        "stage": 1,
        "purpose": "Stage I aaogen resolution + Stage II real-data denominator composition",
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
            "stage2_denominator_model": (
                "real epgamma data are fit locally with floated aaogen-pi0 and "
                "dvcsgen BH/DVCS templates under multiple discriminator choices. "
                "No generator relative normalization is imposed. The nominal "
                "reference is mx2_ep_x_probe_m2; discriminator spread tests shape "
                "model dependence. Mixed-data wrong-tag templates remain "
                "diagnostic-only stress tests."
            ),
        },
    }

    args_dict = vars(args).copy()
    summaries: List[Dict[str, object]] = []
    all_resolution_rows: List[Dict[str, object]] = []
    stage2_summaries: List[Dict[str, object]] = []
    all_stage2_rows: List[Dict[str, object]] = []
    all_stage2_spread_rows: List[Dict[str, object]] = []
    all_mixed_diag_rows: List[Dict[str, object]] = []
    all_shared_rows: List[Dict[str, object]] = []

    if n_processes == 1:
        for period in selected:
            result = process_period(
                period,
                args_dict,
                str(outroot),
                str(stage2_outroot),
            )
            summaries.append(result["summary"])
            all_resolution_rows.extend(result["resolution_rows"])
            if result["stage2_summary"] is not None:
                stage2_summaries.append(result["stage2_summary"])
                all_stage2_rows.extend(result["stage2_rows"])
                all_stage2_spread_rows.extend(result["stage2_spread_rows"])
                all_mixed_diag_rows.extend(result["mixed_diag_rows"])
                all_shared_rows.extend(result.get("shared_rows", []))
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
                        all_shared_rows.extend(result.get("shared_rows", []))
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
    #endwith

    if not args.stage1_only:
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
                **{f"FD_S{s}": s + 1 for s in range(1, 7)},
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

        with (stage2_outroot / "stage2_summary.json").open("w") as f:
            json.dump(
                {
                    "method": (
                        "Local two-component extended-Poisson fits compare several "
                        "independent discriminator choices. aaogen models exclusive-"
                        "pi0 missing-photon tags and dvcsgen models genuine BH/DVCS "
                        "epgamma; both normalizations float independently. The "
                        "mx2_ep_x_probe_m2 result remains the nominal reference while "
                        "the spread among available discriminators quantifies model "
                        "dependence. Deterministic mixed-data tag samples are "
                        "diagnostic-only stress tests."
                    ),
                    "period_summaries": stage2_summaries,
                },
                f,
                indent=2,
                allow_nan=True,
            )
        #endwith

    log(f"Done. Stage-I outputs are in {outroot}.")
    if not args.stage1_only:
        log(f"Stage-II outputs are in {stage2_outroot}.")
    #endif
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
