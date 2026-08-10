#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors_v002.py

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

  3. The two reconstructed photons in the eppi0 event are recovered from the
     stored eppi0 quantities.  The supplied tree schema does not expose
     gamma1_p/gamma1_theta/gamma2_p/gamma2_theta explicitly.  Importantly,
     gamma_phi1/2 are TRento azimuths, not laboratory azimuths, while
     open_angle_egamma1/2 are lab-frame electron-photon opening angles.
     Therefore this script:
       * constructs the event-by-event Trento basis from the beam and measured
         scattered-electron momentum, with z_T parallel to q;
       * solves for each photon's polar angle about q using its Trento phi and
         measured lab-frame electron-photon opening angle;
       * transforms the resulting photon direction back into the lab frame;
       * uses the reconstructed pi0 four-vector (p2_p,p2_theta,p2_phi) and the
         pi0 -> gamma gamma mass-shell relation to recover photon energies;
       * chooses the physical pair of angular solutions with the best pi0
         four-vector closure.

  4. The reconstructed eppi0 photon closest in angle to the epgamma tag is
     identified as the same tag.  The other photon is the reconstructed probe.
  5. Predicted-vs-reconstructed probe residuals are written and plotted.

Important schema assumptions in v001
------------------------------------
epgamma tree:
    p1_* = proton
    p2_* = reconstructed photon/tag

eppi0 tree:
    p1_* = proton
    p2_* = reconstructed pi0
    gamma_phi1/2 = photon Trento azimuths about q
    open_angle_egamma1/2 = lab-frame electron-photon opening angles.

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

    python derive_photon_efficiency_scale_factors_v002.py

Run one period:

    python derive_photon_efficiency_scale_factors_v002.py --period fa18_inb

Run all available entries:

    python derive_photon_efficiency_scale_factors_v002.py --max-entries 0

If the ROOT angles are known explicitly:

    python derive_photon_efficiency_scale_factors_v002.py --angles rad

The Stage-I defaults intentionally require a reconstructed tag energy
E_tag >= 2 GeV, while no efficiency denominator is formed yet.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
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
    "gamma_phi1", "gamma_phi2",
    "open_angle_egamma1", "open_angle_egamma2",
)

EPPIO_OPTIONAL = (
    "runnum", "evnum", "detector_gamma1", "detector_gamma2",
    "Mh_gammagamma", "Mx2", "Mx2_1", "Mx2_2",
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
    electron_theta: np.ndarray
    electron_phi: np.ndarray
    proton_theta: np.ndarray
    proton_phi: np.ndarray
    pi0_theta: np.ndarray
    pi0_phi: np.ndarray
    gamma_phi1: np.ndarray
    gamma_phi2: np.ndarray
    open_angle_egamma1: np.ndarray
    open_angle_egamma2: np.ndarray
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

    return EPPi0Sample(
        electron_p3=cartesian_from_spherical(e_p, e_theta, e_phi),
        proton_p3=cartesian_from_spherical(p_p, p_theta, p_phi),
        pi0_p3=cartesian_from_spherical(pi_p, pi_theta, pi_phi),
        electron_p=e_p,
        proton_p=p_p,
        pi0_p=pi_p,
        electron_theta=e_theta,
        electron_phi=e_phi,
        proton_theta=p_theta,
        proton_phi=p_phi,
        pi0_theta=pi_theta,
        pi0_phi=pi_phi,
        gamma_phi1=to_radians(arrays["gamma_phi1"], unit),
        gamma_phi2=to_radians(arrays["gamma_phi2"], unit),
        open_angle_egamma1=to_radians(arrays["open_angle_egamma1"], unit),
        open_angle_egamma2=to_radians(arrays["open_angle_egamma2"], unit),
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

    tree = cKDTree(x_pi0[pi0_indices] / component_tolerance)
    distance, local_index = tree.query(
        x_epg[epg_indices] / component_tolerance,
        k=1,
        workers=-1,
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
# Recovery of the two photons stored only indirectly in the eppi0 tree
# -----------------------------------------------------------------------------

def build_trento_basis(
    beam_p3: np.ndarray,
    electron_p3: np.ndarray,
) -> Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Construct an event-by-event Trento-like orthonormal basis in the lab.

    z_T is parallel to q = k - k'.

    y_T is normal to the lepton plane, using k x k'.  x_T = y_T x z_T.
    This is the standard right-handed choice appropriate for a Trento azimuth
    measured about q relative to the lepton plane.

    The exact sign convention used by the producer is validated downstream by
    reconstructed-pi0 four-vector closure; a wrong y_T sign would manifest as
    a large and systematic closure failure.
    """
    q = np.asarray(beam_p3, dtype=float) - np.asarray(electron_p3, dtype=float)
    qmag = float(np.linalg.norm(q))
    if qmag <= 1.0e-14:
        return None
    #endif
    zhat = q / qmag

    normal = np.cross(np.asarray(beam_p3, dtype=float), np.asarray(electron_p3, dtype=float))
    nmag = float(np.linalg.norm(normal))
    if nmag <= 1.0e-14:
        return None
    #endif
    yhat = normal / nmag

    xhat = np.cross(yhat, zhat)
    xmag = float(np.linalg.norm(xhat))
    if xmag <= 1.0e-14:
        return None
    #endif
    xhat /= xmag

    # Re-orthogonalize y to protect against floating-point drift.
    yhat = np.cross(zhat, xhat)
    yhat /= max(float(np.linalg.norm(yhat)), 1.0e-30)

    return xhat, yhat, zhat


def trento_polar_solutions_from_electron_opening(
    electron_hat: np.ndarray,
    trento_phi: float,
    opening_angle_lab: float,
    xhat: np.ndarray,
    yhat: np.ndarray,
    zhat: np.ndarray,
) -> List[float]:
    """
    Solve for the photon polar angle theta_{gamma q} about q.

    In the Trento basis,

      n_gamma =
          sin(theta) cos(phi_T) x_T
        + sin(theta) sin(phi_T) y_T
        + cos(theta)            z_T.

    The measured lab electron-photon opening angle alpha supplies

      e_hat dot n_gamma = cos(alpha).

    Thus

      A cos(theta) + B sin(theta) = C,

    with
      A = e_hat dot z_T
      B = cos(phi_T) e_hat dot x_T + sin(phi_T) e_hat dot y_T
      C = cos(alpha).

    There may be two physical theta solutions; both are retained and the
    reconstructed-pi0 four-vector closure resolves the ambiguity.
    """
    if not math.isfinite(trento_phi) or not math.isfinite(opening_angle_lab):
        return []
    #endif

    A = float(np.dot(electron_hat, zhat))
    B = (
        math.cos(trento_phi) * float(np.dot(electron_hat, xhat))
        + math.sin(trento_phi) * float(np.dot(electron_hat, yhat))
    )
    C = math.cos(opening_angle_lab)
    R = math.hypot(A, B)

    if R <= 1.0e-14:
        return []
    #endif

    ratio = C / R
    if ratio < -1.0 - 1.0e-10 or ratio > 1.0 + 1.0e-10:
        return []
    #endif

    ratio = max(-1.0, min(1.0, ratio))
    delta = math.atan2(B, A)
    beta = math.acos(ratio)

    out: List[float] = []
    for raw in (delta + beta, delta - beta):
        # theta is a polar angle about q and must lie in [0, pi].
        # Add/subtract 2pi until we can test the physical representative.
        candidates = (raw, raw + TWO_PI, raw - TWO_PI)
        for candidate in candidates:
            if 0.0 <= candidate <= math.pi:
                if not any(abs(candidate - old) < 1.0e-10 for old in out):
                    out.append(candidate)
                #endif
                break
            #endif
        #endfor
    #endfor

    return out


def direction_from_trento(
    theta_q: float,
    phi_trento: float,
    xhat: np.ndarray,
    yhat: np.ndarray,
    zhat: np.ndarray,
) -> np.ndarray:
    st = math.sin(theta_q)
    n = (
        st * math.cos(phi_trento) * xhat
        + st * math.sin(phi_trento) * yhat
        + math.cos(theta_q) * zhat
    )
    norm = float(np.linalg.norm(n))
    if norm <= 1.0e-14:
        return np.asarray([float("nan")] * 3, dtype=float)
    #endif
    return n / norm


@dataclass
class RecoveredPi0Photons:
    success: bool
    n1: Optional[np.ndarray] = None
    n2: Optional[np.ndarray] = None
    theta_q1: float = float("nan")
    theta_q2: float = float("nan")
    phi_trento1: float = float("nan")
    phi_trento2: float = float("nan")
    lab_theta1: float = float("nan")
    lab_theta2: float = float("nan")
    lab_phi1: float = float("nan")
    lab_phi2: float = float("nan")
    energy1: float = float("nan")
    energy2: float = float("nan")
    mass_gg: float = float("nan")
    energy_closure: float = float("nan")
    momentum_closure: float = float("nan")
    opening1_closure: float = float("nan")
    opening2_closure: float = float("nan")
    score: float = float("inf")


def recover_pi0_photons(
    beam_p3: np.ndarray,
    electron_p3: np.ndarray,
    pi0_p3: np.ndarray,
    gamma_phi1_trento: float,
    gamma_phi2_trento: float,
    opening1_lab: float,
    opening2_lab: float,
) -> RecoveredPi0Photons:
    """
    Recover the two reconstructed photons using Trento azimuths and lab-frame
    electron-photon opening angles.

    Photon directions are first solved in the event-by-event Trento basis and
    transformed back to lab coordinates.

    For a candidate photon direction n_i and reconstructed pi0 P=(E,p),
    enforcing (P-k_i)^2=0 gives

        E_gamma_i = m_pi0^2 / [2 (E_pi0 - p_pi0 dot n_i)].

    All physical angular-solution pairs are tested.  The pair with the best
    pi0 four-vector closure is retained.
    """
    beam_p3 = np.asarray(beam_p3, dtype=float)
    electron_p3 = np.asarray(electron_p3, dtype=float)
    pi0_p3 = np.asarray(pi0_p3, dtype=float)

    e_mag = float(np.linalg.norm(electron_p3))
    p_pi = float(np.linalg.norm(pi0_p3))
    if e_mag <= 1.0e-14 or not math.isfinite(p_pi):
        return RecoveredPi0Photons(False)
    #endif
    electron_hat = electron_p3 / e_mag

    basis = build_trento_basis(beam_p3, electron_p3)
    if basis is None:
        return RecoveredPi0Photons(False)
    #endif
    xhat, yhat, zhat = basis

    E_pi = math.sqrt(max(0.0, p_pi * p_pi + M_PI0 * M_PI0))

    sol1 = trento_polar_solutions_from_electron_opening(
        electron_hat,
        gamma_phi1_trento,
        opening1_lab,
        xhat, yhat, zhat,
    )
    sol2 = trento_polar_solutions_from_electron_opening(
        electron_hat,
        gamma_phi2_trento,
        opening2_lab,
        xhat, yhat, zhat,
    )

    best = RecoveredPi0Photons(False)

    for theta_q1 in sol1:
        n1 = direction_from_trento(theta_q1, gamma_phi1_trento, xhat, yhat, zhat)
        if not np.all(np.isfinite(n1)):
            continue
        #endif
        denom1 = 2.0 * (E_pi - float(np.dot(pi0_p3, n1)))
        if denom1 <= 1.0e-12:
            continue
        #endif
        E1 = M_PI0 * M_PI0 / denom1

        for theta_q2 in sol2:
            n2 = direction_from_trento(theta_q2, gamma_phi2_trento, xhat, yhat, zhat)
            if not np.all(np.isfinite(n2)):
                continue
            #endif
            denom2 = 2.0 * (E_pi - float(np.dot(pi0_p3, n2)))
            if denom2 <= 1.0e-12:
                continue
            #endif
            E2 = M_PI0 * M_PI0 / denom2

            if (
                not math.isfinite(E1)
                or not math.isfinite(E2)
                or E1 <= 0.0
                or E2 <= 0.0
            ):
                continue
            #endif

            p_sum = E1 * n1 + E2 * n2
            energy_closure = (E1 + E2) - E_pi
            momentum_closure = float(np.linalg.norm(p_sum - pi0_p3))

            cos12 = max(-1.0, min(1.0, float(np.dot(n1, n2))))
            m2 = max(0.0, 2.0 * E1 * E2 * (1.0 - cos12))
            mass_gg = math.sqrt(m2)

            opening1_reco = angle_between(electron_hat, n1)
            opening2_reco = angle_between(electron_hat, n2)
            opening1_closure = opening1_reco - opening1_lab
            opening2_closure = opening2_reco - opening2_lab

            scale = max(E_pi, 0.25)
            score = (
                momentum_closure / scale
                + abs(energy_closure) / scale
                + 0.25 * abs(mass_gg - M_PI0) / M_PI0
                + 0.05 * abs(opening1_closure)
                + 0.05 * abs(opening2_closure)
            )

            if score < best.score:
                lab_theta1 = math.acos(max(-1.0, min(1.0, float(n1[2]))))
                lab_theta2 = math.acos(max(-1.0, min(1.0, float(n2[2]))))
                lab_phi1 = math.atan2(float(n1[1]), float(n1[0]))
                lab_phi2 = math.atan2(float(n2[1]), float(n2[0]))

                best = RecoveredPi0Photons(
                    True,
                    n1=n1,
                    n2=n2,
                    theta_q1=theta_q1,
                    theta_q2=theta_q2,
                    phi_trento1=gamma_phi1_trento,
                    phi_trento2=gamma_phi2_trento,
                    lab_theta1=lab_theta1,
                    lab_theta2=lab_theta2,
                    lab_phi1=lab_phi1,
                    lab_phi2=lab_phi2,
                    energy1=E1,
                    energy2=E2,
                    mass_gg=mass_gg,
                    energy_closure=energy_closure,
                    momentum_closure=momentum_closure,
                    opening1_closure=opening1_closure,
                    opening2_closure=opening2_closure,
                    score=score,
                )
            #endif
        #endfor
    #endfor

    return best


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
    reco_probe_energy: List[float]
    pred_probe_energy: List[float]
    pred_probe_p: List[float]
    pred_probe_mass2: List[float]
    pred_probe_E_minus_p: List[float]
    probe_delta_E: List[float]
    probe_delta_E_over_E: List[float]
    probe_delta_theta_deg: List[float]
    probe_delta_phi_deg: List[float]
    probe_opening_residual_deg: List[float]
    tag_reco_match_angle_deg: List[float]
    pi0_energy_closure: List[float]
    pi0_momentum_closure: List[float]
    pi0_mass_recovered: List[float]
    pi0_recovery_score: List[float]
    reco_gamma_gamma_opening_deg: List[float]
    detector_tag: List[int]
    detector_probe: List[int]

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
    max_tag_reco_angle_deg: float,
) -> Tuple[PairArrays, Dict[str, int]]:
    out = PairArrays.empty()
    counters = {
        "matched_parent_pairs": int(len(matches.epg_index)),
        "pi0_photon_recovery_failed": 0,
        "tag_reco_match_failed": 0,
        "nonphysical_predicted_probe": 0,
        "accepted_stage1_pairs": 0,
    }

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    beam4 = np.asarray([period.beam_energy, 0.0, 0.0, beam_p], dtype=float)
    target4 = np.asarray([M_P, 0.0, 0.0, 0.0], dtype=float)
    max_tag_angle = math.radians(max_tag_reco_angle_deg)

    for k in range(len(matches.epg_index)):
        i = int(matches.epg_index[k])
        j = int(matches.eppi0_index[k])

        recovered = recover_pi0_photons(
            beam_p3=beam4[1:4],
            electron_p3=eppi0.electron_p3[j],
            pi0_p3=eppi0.pi0_p3[j],
            gamma_phi1_trento=float(eppi0.gamma_phi1[j]),
            gamma_phi2_trento=float(eppi0.gamma_phi2[j]),
            opening1_lab=float(eppi0.open_angle_egamma1[j]),
            opening2_lab=float(eppi0.open_angle_egamma2[j]),
        )

        if not recovered.success:
            counters["pi0_photon_recovery_failed"] += 1
            continue
        #endif

        n1 = recovered.n1
        n2 = recovered.n2
        assert n1 is not None and n2 is not None
        ntag = epg.tag_p3[i] / max(float(np.linalg.norm(epg.tag_p3[i])), 1.0e-30)

        a1 = angle_between(ntag, n1)
        a2 = angle_between(ntag, n2)

        if not math.isfinite(a1) or not math.isfinite(a2):
            counters["tag_reco_match_failed"] += 1
            continue
        #endif

        if a1 <= a2:
            tag_id = 1
            tag_match = a1
            reco_probe_n = n2
            reco_probe_E = recovered.energy2
            reco_probe_theta = recovered.lab_theta2
            reco_probe_phi = recovered.lab_phi2
        else:
            tag_id = 2
            tag_match = a2
            reco_probe_n = n1
            reco_probe_E = recovered.energy1
            reco_probe_theta = recovered.lab_theta1
            reco_probe_phi = recovered.lab_phi1
        #endif

        if tag_match > max_tag_angle:
            counters["tag_reco_match_failed"] += 1
            continue
        #endif

        electron4 = event_four_vector(epg.electron_p3[i], M_E)
        proton4 = event_four_vector(epg.proton_p3[i], M_P)
        tag4 = photon_four_vector(epg.tag_p3[i])
        probe4 = beam4 + target4 - electron4 - proton4 - tag4

        Epred = float(probe4[0])
        ppred3 = probe4[1:4]
        ppred = float(np.linalg.norm(ppred3))
        m2pred = Epred * Epred - ppred * ppred

        if Epred <= 0.0 or ppred <= 0.0 or not np.all(np.isfinite(probe4)):
            counters["nonphysical_predicted_probe"] += 1
            continue
        #endif

        npred = ppred3 / ppred
        theta_pred = math.acos(max(-1.0, min(1.0, float(npred[2]))))
        phi_pred = math.atan2(float(npred[1]), float(npred[0]))

        dtheta = reco_probe_theta - theta_pred
        dphi = float(wrap_phi(reco_probe_phi - phi_pred))
        open_res = angle_between(npred, reco_probe_n)
        gg_open = angle_between(n1, n2)

        out.epg_index.append(i)
        out.eppi0_index.append(j)
        out.parent_distance.append(float(matches.nearest_distance[k]))
        out.parent_max_component_delta.append(float(matches.max_component_delta[k]))
        out.tag_energy.append(float(epg.tag_energy[i]))
        out.reco_probe_energy.append(float(reco_probe_E))
        out.pred_probe_energy.append(Epred)
        out.pred_probe_p.append(ppred)
        out.pred_probe_mass2.append(m2pred)
        out.pred_probe_E_minus_p.append(Epred - ppred)
        out.probe_delta_E.append(float(reco_probe_E - Epred))
        out.probe_delta_E_over_E.append(float((reco_probe_E - Epred) / Epred))
        out.probe_delta_theta_deg.append(math.degrees(dtheta))
        out.probe_delta_phi_deg.append(math.degrees(dphi))
        out.probe_opening_residual_deg.append(math.degrees(open_res))
        out.tag_reco_match_angle_deg.append(math.degrees(tag_match))
        out.pi0_energy_closure.append(float(recovered.energy_closure))
        out.pi0_momentum_closure.append(float(recovered.momentum_closure))
        out.pi0_mass_recovered.append(float(recovered.mass_gg))
        out.pi0_recovery_score.append(float(recovered.score))
        out.reco_gamma_gamma_opening_deg.append(math.degrees(gg_open))

        if tag_id == 1:
            out.detector_tag.append(detector_value(eppi0.raw, "detector_gamma1", j))
            out.detector_probe.append(detector_value(eppi0.raw, "detector_gamma2", j))
        else:
            out.detector_tag.append(detector_value(eppi0.raw, "detector_gamma2", j))
            out.detector_probe.append(detector_value(eppi0.raw, "detector_gamma1", j))
        #endif

        counters["accepted_stage1_pairs"] += 1
    #endfor

    return out, counters


def pair_arrays_to_numpy(pairs: PairArrays) -> Dict[str, np.ndarray]:
    out: Dict[str, np.ndarray] = {}
    integer_names = {"epg_index", "eppi0_index", "detector_tag", "detector_probe"}
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
    bins: Tuple[int, int] = (100, 100),
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
) -> None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    h = ax.hist2d(x, y, bins=bins)
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


def robust_symmetric_limit(values: np.ndarray, percentile: float = 99.0, floor: float = 1.0e-6) -> float:
    values = np.asarray(values, dtype=float)
    values = np.abs(values[np.isfinite(values)])
    if values.size == 0:
        return floor
    #endif
    return max(floor, float(np.percentile(values, percentile)))


def make_plots(period: PeriodConfig, arrays: Dict[str, np.ndarray], outdir: Path) -> None:
    if len(arrays["pred_probe_energy"]) == 0:
        log(f"{period.label}: no accepted pairs; skipping plots.")
        return
    #endif

    title = period.label

    save_hist(
        arrays["parent_max_component_delta"],
        outdir / "parent_match_max_cartesian_delta.png",
        r"max $|\Delta p_{e/p,x/y/z}|$ [GeV]",
        f"{title}: e/p kinematic match quality",
        bins=120,
        logy=True,
    )

    save_hist(
        arrays["tag_reco_match_angle_deg"],
        outdir / "tag_reco_angular_match.png",
        r"$\Delta\alpha(\gamma_{\rm tag}^{ep\gamma},\gamma_{\rm tag}^{ep\pi^0})$ [deg]",
        f"{title}: tag correspondence after e/p event matching",
        bins=120,
        logy=True,
    )

    emax = max(
        2.5,
        float(np.nanpercentile(
            np.concatenate((arrays["pred_probe_energy"], arrays["reco_probe_energy"])),
            99.5,
        )),
    )
    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["reco_probe_energy"],
        outdir / "probe_energy_pred_vs_reco.png",
        r"$E_{\rm probe}^{pred}$ [GeV]",
        r"$E_{\rm probe}^{reco}$ [GeV]",
        f"{title}: predicted vs reconstructed probe energy",
        bins=(120, 120),
        xlim=(0.0, emax),
        ylim=(0.0, emax),
    )

    de_lim = robust_symmetric_limit(arrays["probe_delta_E"], 99.0, 0.05)
    save_hist(
        arrays["probe_delta_E"],
        outdir / "probe_energy_residual.png",
        r"$E_{\rm probe}^{reco}-E_{\rm probe}^{pred}$ [GeV]",
        f"{title}: probe energy residual",
        bins=160,
        xlim=(-de_lim, de_lim),
    )

    rel_lim = min(2.0, robust_symmetric_limit(arrays["probe_delta_E_over_E"], 99.0, 0.05))
    save_hist(
        arrays["probe_delta_E_over_E"],
        outdir / "probe_fractional_energy_residual.png",
        r"$(E_{\rm reco}-E_{\rm pred})/E_{\rm pred}$",
        f"{title}: fractional probe energy residual",
        bins=160,
        xlim=(-rel_lim, rel_lim),
    )

    th_lim = min(5.0, robust_symmetric_limit(arrays["probe_delta_theta_deg"], 99.5, 0.1))
    ph_lim = min(8.0, robust_symmetric_limit(arrays["probe_delta_phi_deg"], 99.5, 0.1))

    save_hist(
        arrays["probe_delta_theta_deg"],
        outdir / "probe_theta_residual.png",
        r"$\theta_{\rm reco}-\theta_{\rm pred}$ [deg]",
        f"{title}: probe polar-angle residual",
        bins=160,
        xlim=(-th_lim, th_lim),
    )
    save_hist(
        arrays["probe_delta_phi_deg"],
        outdir / "probe_phi_residual.png",
        r"$\phi_{\rm reco}-\phi_{\rm pred}$ [deg]",
        f"{title}: probe azimuth residual",
        bins=160,
        xlim=(-ph_lim, ph_lim),
    )

    open_hi = min(
        10.0,
        max(0.25, float(np.nanpercentile(arrays["probe_opening_residual_deg"], 99.5))),
    )
    save_hist(
        arrays["probe_opening_residual_deg"],
        outdir / "probe_direction_opening_residual.png",
        r"$\Delta\alpha(\gamma_{\rm probe}^{pred},\gamma_{\rm probe}^{reco})$ [deg]",
        f"{title}: probe direction matching resolution",
        bins=160,
        xlim=(0.0, open_hi),
        logy=True,
    )

    m2lim = robust_symmetric_limit(arrays["pred_probe_mass2"], 99.0, 0.02)
    save_hist(
        arrays["pred_probe_mass2"],
        outdir / "predicted_probe_mass2.png",
        r"$(p_{\rm probe}^{pred})^2$ [GeV$^2$]",
        f"{title}: predicted missing-photon mass shell",
        bins=160,
        xlim=(-m2lim, m2lim),
    )

    eplim = robust_symmetric_limit(arrays["pred_probe_E_minus_p"], 99.0, 0.02)
    save_hist(
        arrays["pred_probe_E_minus_p"],
        outdir / "predicted_probe_E_minus_p.png",
        r"$E_{\rm probe}^{pred}-|\vec p_{\rm probe}^{pred}|$ [GeV]",
        f"{title}: predicted photon mass-shell residual",
        bins=160,
        xlim=(-eplim, eplim),
    )

    save_hist(
        arrays["pi0_mass_recovered"],
        outdir / "recovered_gamma_gamma_mass.png",
        r"$M_{\gamma\gamma}^{recovered}$ [GeV]",
        f"{title}: eppi0 photon-recovery closure",
        bins=140,
        xlim=(0.05, 0.22),
    )

    save_hist(
        arrays["pi0_momentum_closure"],
        outdir / "recovered_pi0_momentum_closure.png",
        r"$|\vec p_{\gamma1}+\vec p_{\gamma2}-\vec p_{\pi^0}|$ [GeV]",
        f"{title}: reconstructed-pi0 momentum closure",
        bins=140,
        logy=True,
    )

    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["probe_opening_residual_deg"],
        outdir / "probe_direction_residual_vs_energy.png",
        r"$E_{\rm probe}^{pred}$ [GeV]",
        r"$\Delta\alpha_{\rm probe}$ [deg]",
        f"{title}: probe direction resolution vs energy",
        bins=(100, 100),
        xlim=(0.0, emax),
        ylim=(0.0, open_hi),
    )

    save_hist2d(
        arrays["pred_probe_energy"],
        arrays["probe_delta_E_over_E"],
        outdir / "probe_fractional_energy_residual_vs_energy.png",
        r"$E_{\rm probe}^{pred}$ [GeV]",
        r"$(E_{\rm reco}-E_{\rm pred})/E_{\rm pred}$",
        f"{title}: probe energy closure vs energy",
        bins=(100, 100),
        xlim=(0.0, emax),
        ylim=(-rel_lim, rel_lim),
    )

    save_hist(
        arrays["reco_gamma_gamma_opening_deg"],
        outdir / "reco_gamma_gamma_opening_angle.png",
        r"$\Delta\alpha_{\gamma\gamma}$ [deg]",
        f"{title}: reconstructed pi0 photon opening angle",
        bins=140,
        logy=True,
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
        "median_tag_reco_match_deg": percentile_or_nan(
            arrays["tag_reco_match_angle_deg"], 50.0
        ),
        "p99_tag_reco_match_deg": percentile_or_nan(
            arrays["tag_reco_match_angle_deg"], 99.0
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
        "median_probe_deltaE_GeV": percentile_or_nan(arrays["probe_delta_E"], 50.0),
        "rms_probe_deltaE_GeV": rms_or_nan(arrays["probe_delta_E"]),
        "median_probe_fractional_deltaE": percentile_or_nan(
            arrays["probe_delta_E_over_E"], 50.0
        ),
        "median_abs_pred_probe_mass2_GeV2": percentile_or_nan(
            np.abs(arrays["pred_probe_mass2"]), 50.0
        ),
        "p95_abs_pred_probe_mass2_GeV2": percentile_or_nan(
            np.abs(arrays["pred_probe_mass2"]), 95.0
        ),
        "median_pi0_momentum_closure_GeV": percentile_or_nan(
            arrays["pi0_momentum_closure"], 50.0
        ),
        "median_recovered_Mgg_GeV": percentile_or_nan(
            arrays["pi0_mass_recovered"], 50.0
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
            "momentum component, GeV. Default: 0.002."
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
            "Maximum angle in degrees between epgamma tag and the corresponding "
            "recovered eppi0 photon. Default: 1.0 deg."
        ),
    )
    return parser.parse_args()


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
        "No photon-efficiency denominator/numerator is formed in v002; "
        "this run validates event matching and missing-photon prediction only."
    )

    preflight(selected)

    summaries: List[Dict[str, object]] = []
    provenance = {
        "script": Path(__file__).name,
        "stage": 1,
        "purpose": "aaogen tag-and-probe kinematic validation",
        "arguments": vars(args),
        "periods": [asdict(p) for p in selected],
        "schema_assumptions": {
            "epgamma": "p1=proton, p2=photon/tag",
            "eppi0": (
                "p1=proton, p2=pi0; gamma_phi1/2 are Trento azimuths and "
                "open_angle_egamma1/2 are lab e-gamma opening angles; individual "
                "photon lab directions/energies are recovered with the Trento "
                "basis plus the reconstructed pi0 four-vector"
            ),
        },
    }

    for period in selected:
        log(f"{period.label}: reading aaogen epgamma sample.")
        epg_arrays, epg_tree, epg_total = read_branches(
            period.epgamma_pi0_mc,
            EPG_REQUIRED,
            EPG_OPTIONAL,
            args.tree,
            args.max_entries,
        )
        epg = extract_epgamma(epg_arrays, args.angles)

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
            args.tree,
            args.max_entries,
        )
        eppi0 = extract_eppi0(pi_arrays, args.angles)

        log(
            f"{period.label}: eppi0 tree '{pi_tree}', "
            f"loaded {len(eppi0.pi0_p):,}/{pi_total:,} entries, "
            f"angles interpreted as {eppi0.angle_unit}."
        )

        log(f"{period.label}: matching e/p parent kinematics between MC skims.")
        matches = match_parent_kinematics(
            epg,
            eppi0,
            tag_min=args.tag_min,
            tag_max=args.tag_max,
            component_tolerance=args.parent_component_tol,
            nearest_distance_max=args.parent_distance_max,
        )
        log(f"{period.label}: accepted {len(matches.epg_index):,} parent matches.")

        log(f"{period.label}: recovering eppi0 photons and validating predicted probe.")
        pairs, counters = build_stage1_pairs(
            period,
            epg,
            eppi0,
            matches,
            max_tag_reco_angle_deg=args.max_tag_reco_angle,
        )
        pair_np = pair_arrays_to_numpy(pairs)

        period_dir = outroot / period.key
        period_dir.mkdir(parents=True, exist_ok=True)

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
        summaries.append(summary)

        with (period_dir / "stage1_summary.json").open("w") as f:
            json.dump(summary, f, indent=2, allow_nan=True)
        #endif

        log(
            f"{period.label}: Stage-I accepted pairs = "
            f"{counters['accepted_stage1_pairs']:,}; "
            f"photon-recovery failures = {counters['pi0_photon_recovery_failed']:,}; "
            f"tag-correspondence failures = {counters['tag_reco_match_failed']:,}."
        )
    #endfor

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

    log(f"Done. Stage-I outputs are in {outroot}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
