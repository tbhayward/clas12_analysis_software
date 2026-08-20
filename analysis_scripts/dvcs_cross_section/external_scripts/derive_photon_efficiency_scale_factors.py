#!/usr/bin/env python3
"""
---------------
For the nSidis denominator, the nominal pi0 fraction is now the simultaneous
shared fit to Delta_phi, pTmiss, and xF2 in each detector/energy bin.
The three one-variable fits are retained as model alternatives. The reported
composition systematic on epsilon_data and SF_gamma is the envelope obtained
by recalculating the full result with each alternative fraction. The RMS and
maximum |Delta f_pi0| are also written as diagnostics, but the envelope is the
assigned systematic because it is conservative and directly propagated through
the nonlinear 1/f_pi0 efficiency dependence.
derive_photon_efficiency_scale_factors_v067_wrapup_nominal.py

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

The production denominator/template selection now also applies Valerii
Klimenko's photon-efficiency exclusivity cuts: -0.231 < MM^2(epX) < 0.309
GeV^2, MM^2(e gamma X) > 1.4 GeV^2, |Delta phi(p,gamma)| < 5.7 deg, and
Angle(gamma,X) diagnostic.

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
component with the two-component morph nuisance values held fixed. This does
not feed the nominal efficiency extraction.

Routine outputs:
    output/photon_efficiency/summary/
    output/photon_efficiency/stage2/
    output/photon_efficiency/stage3/

Default --diagnostics selection runs only the focused discriminator-selection
study. --diagnostics full restores the expensive closure/profile/coarse-FT
suite. --plot-mode compact remains the default plotting mode.

Typical full-statistics run:
    python derive_photon_efficiency_scale_factors_v067_wrapup_nominal.py --max-entries 0

One period:
    python derive_photon_efficiency_scale_factors_v067_wrapup_nominal.py \
        --max-entries 0 --period fa18_out


LOOSE nSIDIS TRANSITION STUDY
-----------------------------
Fa18 Inb and Sp19 Inb now have loose nSidis-derived epgammaX and eppi0X data trees whose only
upstream event requirement is the nSidis electron selection.  These trees do
not impose the old DVCS/DVPi0P wagon missing-energy requirement.

Use:
    python derive_photon_efficiency_scale_factors_v067_wrapup_nominal.py \
        --max-entries 0 --period fa18_inb --nsidis-study

This isolated mode:
  * keeps and reads the original wagon epgamma/eppi0 trees;
  * checks exact runnum+evnum event recovery in the nSidis trees;
  * compares old and new epgamma kinematics below 2 GeV;
  * conditions old-wagon overlap checks on the nSidis parent requirement
    e_p > 2 GeV;
  * derives the FT and FD theta acceptance directly from real reconstructed
    photons (detector 0 = FT, detector 1 = FD), excluding the FT beam hole and
    out-of-acceptance predicted probes;
  * applies theta_ep > 5 deg identically to data and MC templates;
  * reconstructs the sharp M_X^2(ep pi0) exclusive core, scans M_X^2 + missing-momentum cuts, and quantifies rejection of the loose SIDIS-rich
    eppi0X population; E_miss is diagnostic-only;
  * scans the missing-photon mass-shell support without cutting pTmiss;
  * optionally runs the finalized M_X^2(ep) x pTmiss composition model below
    and above 2 GeV, including support-cut stability and actual fit projections;
  * does NOT form a new efficiency scale factor.

After inspecting those outputs, an optional diagnostic composition pilot is:
    --nsidis-pilot-fit

which uses the finalized M_X^2(ep) x pTmiss denominator model, compares the
old-wagon and nSidis results below 2 GeV, scans |M_probe^2| support, and extends
the diagnostic fit above 2 GeV. It remains diagnostic until the nSidis
reconstructed-pi0 exclusivity prescription is accepted.

"""

from __future__ import annotations

import argparse
import ast
import inspect
import textwrap
import csv
import gc
import json
import math
import os
import re
import sys
import time
import warnings
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
from scipy.optimize import brentq, minimize, minimize_scalar
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

THETA_EP_MIN_DEG = 5.0  # legacy diagnostic name
THETA_EGAMMA_MIN_DEG = 5.0

# Valerii Klimenko photon-efficiency exclusivity cuts (June 30, 2026):
#   -0.231 < MM^2(epX) < 0.309 GeV^2
#   MM^2(e gamma X) > 1.4 GeV^2
#   -5.7 deg < Delta phi(p,gamma) < 5.7 deg
#   Angle(gamma,X) is retained as a diagnostic only
# Here Delta phi(p,gamma) is the signed coplanarity residual relative to
# 180 deg, and X is the photon predicted from four-momentum conservation
# in e p -> e p gamma X.
VALERII_EGAMMA_MM2_MIN_GEV2 = 1.4
VALERII_DPHI_PGAMMA_MAX_DEG = 5.7
VALERII_ANGLE_EX_MAX_DEG = 9.2

PROBE_ENERGY_EDGES = np.asarray(
    # The original three sub-0.8 GeV bins produced a poorly constrained FT
    # pi0 fraction and correspondingly enormous profile-likelihood tails.
    # Use one 0.40-0.80 GeV bin.  The former 4.5-6 and 6-10 GeV bins are also
    # kept merged because of sparse FT statistics. stage2_energy_edges()
    # truncates the upper edge to the requested endpoint (9.5 GeV for nSidis).
    [0.40, 0.80, 1.00, 1.25, 1.50, 2.00, 3.00, 4.50, 10.00],
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
    nsidis_epgamma_data: Optional[str] = None
    nsidis_eppi0_data: Optional[str] = None
    clasdis_epgamma_mc: Optional[str] = None
    clasdis_eppi0_mc: Optional[str] = None


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        "fa18_inb", "Fa18 Inb", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/aaogen_rga_fa18_inb_eppi0_0.40GeV.root",
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_fa18_inb_epgamma.root"
        ),
        nsidis_eppi0_data=(
            "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
            "rga_fa18_inb_eppi0.root"
        ),
        clasdis_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "clasdis/rec_clasdis_rga_fa18_inb_epgammaX.root"
        ),
        clasdis_eppi0_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "clasdis/rec_clasdis_rga_fa18_inb_eppi0X.root"
        ),
    ),
    PeriodConfig(
        "fa18_out", "Fa18 Out", 10.604,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_fa18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/aaogen_rga_fa18_out_eppi0_0.40GeV.root",
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_fa18_out_epgamma.root"
        ),
        nsidis_eppi0_data=(
            "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
            "rga_fa18_out_eppi0.root"
        ),
        clasdis_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "clasdis/rec_clasdis_rga_fa18_out_epgammaX.root"
        ),
        clasdis_eppi0_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "clasdis/rec_clasdis_rga_fa18_out_eppi0X.root"
        ),
    ),
    PeriodConfig(
        "sp18_inb", "Sp18 Inb", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/aaogen_rga_sp18_inb_eppi0_0.40GeV.root",
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_sp18_inb_epgamma.root"
        ),
        nsidis_eppi0_data=(
            "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
            "rga_sp18_inb_eppi0.root"
        ),
    ),
    PeriodConfig(
        "sp18_out", "Sp18 Out", 10.594,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp18_out_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/aaogen_rga_sp18_out_eppi0_0.40GeV.root",
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_sp18_out_epgamma.root"
        ),
        nsidis_eppi0_data=(
            "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
            "rga_sp18_out_eppi0.root"
        ),
    ),
    PeriodConfig(
        "sp19_inb", "Sp19 Inb", 10.1998,
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/bkg_rga_sp19_inb_epgamma_0.40GeV.root",
        "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/dvcsgen_files_greater_than_0.40GeV/aaogen_rga_sp19_inb_eppi0_0.40GeV.root",
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_sp19_inb_epgamma.root"
        ),
        nsidis_eppi0_data=(
            "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/"
            "rga_sp19_inb_eppi0.root"
        ),
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
        p * st * np.cos(phi).astype(np.float32, copy=False),
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


def read_branches_with_progress(
    path: str,
    required: Sequence[str],
    optional: Sequence[str],
    tree_name: Optional[str],
    max_entries: int,
    *,
    progress_label: str,
    chunk_entries: int = 500_000,
) -> Tuple[Dict[str, np.ndarray], str, int]:
    """Read scalar ROOT branches in chunks with visible progress."""
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

        chosen = list(required) + [
            b for b in optional
            if b in available and b not in required
        ]

        total_entries = int(tree.num_entries)
        target_entries = (
            total_entries
            if max_entries == 0
            else min(int(max_entries), total_entries)
        )

        log(
            f"{progress_label}: reading {len(chosen)} branches, "
            f"{target_entries:,}/{total_entries:,} entries "
            f"in {int(chunk_entries):,}-entry chunks."
        )

        if target_entries == 0:
            return {}, found_name, total_entries
        #endif

        output: Dict[str, np.ndarray] = {}
        start = 0
        t0 = time.perf_counter()

        while start < target_entries:
            stop = min(
                start + int(chunk_entries),
                target_entries,
            )
            chunk = tree.arrays(
                chosen,
                entry_start=start,
                entry_stop=stop,
                library="np",
            )

            if not output:
                for name in chosen:
                    values = np.asarray(chunk[name])
                    output[name] = np.empty(
                        (target_entries,) + tuple(values.shape[1:]),
                        dtype=values.dtype,
                    )
                #endfor
            #endif

            for name in chosen:
                output[name][start:stop] = np.asarray(
                    chunk[name]
                )
            #endfor

            start = stop
            elapsed = time.perf_counter() - t0
            rate = (
                float(start) / elapsed
                if elapsed > 0.0
                else float("nan")
            )
            log(
                f"{progress_label}: "
                f"{start:,}/{target_entries:,} "
                f"({100.0 * start / target_entries:.1f}%), "
                f"{rate / 1.0e6:.2f} M entries/s."
            )
        #endwhile

        return output, found_name, total_entries




# -----------------------------------------------------------------------------
# Schema adapters
# -----------------------------------------------------------------------------

EPG_REQUIRED = (
    "e_p", "e_theta", "e_phi",
    "p1_p", "p1_theta", "p1_phi",
    "p2_p", "p2_theta", "p2_phi",
    "Mx2", "Mx2_1", "Mx2_2",
)

EPG_OPTIONAL_PI0_MC = (
    "detector2",
    "pTmiss",
    "pT",
    "theta_gamma_gamma",
    "Delta_phi",
    "Emiss2",
    "xF2",
)

EPG_OPTIONAL_DVCS_MC = (
    "pTmiss",
    "pT",
    "theta_gamma_gamma",
    "Delta_phi",
    "Emiss2",
    "xF2",
)

EPG_OPTIONAL_CLASDIS_MC = (
    "pTmiss",
    "pT",
    "theta_gamma_gamma",
    "Delta_phi",
    "Emiss2",
    "xF2",
)

EPG_OPTIONAL_DATA = (
    "runnum",
    "evnum",
    "detector2",
    "pTmiss",
    "pT",
    "theta_gamma_gamma",
    "Delta_phi",
    "Emiss2",
    "xF2",
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


# -----------------------------------------------------------------------------
# Stage-1 "grand diagnostic" inventory
# -----------------------------------------------------------------------------
#
# This is deliberately separate from the normal production branch list.  The
# full inventory is read only for --stage1-only, so the normal efficiency run
# does not pay the I/O/memory cost of dozens of extra branches.
GRAND_STAGE1_OPTIONAL_BRANCHES = (
    # lab/opening-angle kinematics
    "open_angle_ep",
    "open_angle_ep1",
    "open_angle_ep2",
    "open_angle_p1p2",
    # vertices
    "vz_e",
    "vz_p1",
    "vz_p2",
    # inclusive/exclusive kinematics
    "W",
    "x",
    "t",
    "t1",
    "t2",
    "tmin",
    "y",
    "Mx2",
    "Mx2_1",
    "Mx2_2",
    # SIDIS
    "z",
    "xF",
    "pT",
    "eta",
    "eta_gN",
    "xi",
    # SIDIS dihadron
    "z1",
    "z2",
    "xF1",
    "xF2",
    "xi1",
    "xi2",
    "Mh",
    "pT1",
    "pT2",
    "pTpT",
    "eta1",
    "eta2",
    "Delta_eta",
    "eta1_gN",
    "eta2_gN",
    # angles
    "phi1",
    "phi2",
    "Delta_phi",
    "phih",
    "phi",
    "phiR",
    "theta",
    # depolarization factors: support both common tree spellings
    "DepA",
    "DepB",
    "DepC",
    "DepV",
    "DepW",
    "Depolarization_A",
    "Depolarization_B",
    "Depolarization_C",
    "Depolarization_V",
    "Depolarization_W",
    # exclusivity
    "Emiss2",
    "theta_gamma_gamma",
    "pTmiss",
)

# Each tuple is:
#   (display_key, branch aliases, plot label, transform)
#
# "lab_angle_deg" converts the raw lab-angle branch according to the detected
# input angle convention.  Other analysis angles are kept in the ROOT-tree
# convention because they are derived observables rather than detector angles.
GRAND_STAGE1_GROUPS = {
    "lab_kinematics": (
        ("e_p", ("e_p",), r"$p_e$ (GeV)", "identity"),
        ("e_theta", ("e_theta",), r"$\theta_e$ (deg)", "lab_angle_deg"),
        ("e_phi", ("e_phi",), r"$\phi_e$ (deg)", "lab_angle_deg"),
        ("p1_p", ("p1_p",), r"$p_1$ (GeV)", "identity"),
        ("p1_theta", ("p1_theta",), r"$\theta_1$ (deg)", "lab_angle_deg"),
        ("p1_phi", ("p1_phi",), r"$\phi_1^{lab}$ (deg)", "lab_angle_deg"),
        ("p2_p", ("p2_p",), r"$p_2$ (GeV)", "identity"),
        ("p2_theta", ("p2_theta",), r"$\theta_2$ (deg)", "lab_angle_deg"),
        ("p2_phi", ("p2_phi",), r"$\phi_2^{lab}$ (deg)", "lab_angle_deg"),
        ("open_angle_ep", ("open_angle_ep",), "open_angle_ep", "identity"),
        ("open_angle_ep1", ("open_angle_ep1",), "open_angle_ep1", "identity"),
        ("open_angle_ep2", ("open_angle_ep2",), "open_angle_ep2", "identity"),
        ("open_angle_p1p2", ("open_angle_p1p2",), "open_angle_p1p2", "identity"),
    ),
    "inclusive_exclusive_kinematics": (
        ("W", ("W",), r"$W$ (GeV)", "identity"),
        ("x", ("x",), r"$x_B$", "identity"),
        ("t", ("t",), r"$t$ (GeV$^2$)", "identity"),
        ("t1", ("t1",), r"$t_1$ (GeV$^2$)", "identity"),
        ("t2", ("t2",), r"$t_2$ (GeV$^2$)", "identity"),
        ("tmin", ("tmin",), r"$t_{min}$ (GeV$^2$)", "identity"),
        ("y", ("y",), r"$y$", "identity"),
        ("Mx2", ("Mx2",), r"$M_X^2$ (GeV$^2$)", "identity"),
        ("Mx2_1", ("Mx2_1",), r"$M_{X,1}^2$ (GeV$^2$)", "identity"),
        ("Mx2_2", ("Mx2_2",), r"$M_{X,2}^2$ (GeV$^2$)", "identity"),
        ("selection_Mx2_ep", (), r"selection $M_X^2(ep)$ [recomputed] (GeV$^2$)", "derived_Mx2_ep"),
        (
            "selection_Mx2_epgamma",
            (),
            r"selection $M_X^2(ep\gamma)$ [recomputed] (GeV$^2$)",
            "derived_Mx2_epgamma",
        ),
    ),
    "sidis": (
        ("z", ("z",), r"$z$", "identity"),
        ("xF", ("xF",), r"$x_F$", "identity"),
        ("pT", ("pT",), r"$p_T$ (GeV)", "identity"),
        ("eta", ("eta",), r"$\eta$", "identity"),
        ("eta_gN", ("eta_gN",), r"$\eta_{gN}$", "identity"),
        ("xi", ("xi",), r"$\xi$", "identity"),
    ),
    "sidis_dihadron": (
        ("z1", ("z1",), r"$z_1$", "identity"),
        ("z2", ("z2",), r"$z_2$", "identity"),
        ("xF1", ("xF1",), r"$x_{F,1}$", "identity"),
        ("xF2", ("xF2",), r"$x_{F,2}$", "identity"),
        ("xi1", ("xi1",), r"$\xi_1$", "identity"),
        ("xi2", ("xi2",), r"$\xi_2$", "identity"),
        ("Mh", ("Mh",), r"$M_h$ (GeV)", "identity"),
        ("pT1", ("pT1",), r"$p_{T,1}$ (GeV)", "identity"),
        ("pT2", ("pT2",), r"$p_{T,2}$ (GeV)", "identity"),
        ("pTpT", ("pTpT",), r"$p_{T}p_{T}$", "identity"),
        ("eta1", ("eta1",), r"$\eta_1$", "identity"),
        ("eta2", ("eta2",), r"$\eta_2$", "identity"),
        ("Delta_eta", ("Delta_eta",), r"$\Delta\eta$", "identity"),
        ("eta1_gN", ("eta1_gN",), r"$\eta_{1,gN}$", "identity"),
        ("eta2_gN", ("eta2_gN",), r"$\eta_{2,gN}$", "identity"),
    ),
    "angles": (
        ("phi1", ("phi1",), r"$\phi_1$", "identity"),
        ("phi2", ("phi2",), r"$\phi_2$", "identity"),
        ("Delta_phi", ("Delta_phi",), r"$\Delta\phi$", "identity"),
        ("phih", ("phih", "phi"), r"$\phi_h$", "identity"),
        ("phiR", ("phiR",), r"$\phi_R$", "identity"),
        ("theta", ("theta",), r"$\theta$", "identity"),
    ),
    "depolarization": (
        ("DepA", ("DepA", "Depolarization_A"), r"$A$", "identity"),
        ("DepB", ("DepB", "Depolarization_B"), r"$B$", "identity"),
        ("DepC", ("DepC", "Depolarization_C"), r"$C$", "identity"),
        ("DepV", ("DepV", "Depolarization_V"), r"$V$", "identity"),
        ("DepW", ("DepW", "Depolarization_W"), r"$W_{dep}$", "identity"),
    ),
    "exclusivity": (
        ("Emiss2", ("Emiss2",), r"$E_{miss}$", "identity"),
        (
            "theta_gamma_gamma",
            ("theta_gamma_gamma",),
            r"$\theta_{\gamma\gamma}$",
            "identity",
        ),
        ("pTmiss", ("pTmiss",), r"$p_{T,miss}$ (GeV)", "identity"),
        ("selection_Mx2_ep", (), r"selection $M_X^2(ep)$ [recomputed] (GeV$^2$)", "derived_Mx2_ep"),
        (
            "selection_Mx2_epgamma",
            (),
            r"selection $M_X^2(ep\gamma)$ [recomputed] (GeV$^2$)",
            "derived_Mx2_epgamma",
        ),
    ),
}


def grand_stage1_optional_branches() -> Tuple[str, ...]:
    """Deduplicated optional branch list for the all-variable Stage-1 scan."""
    names = list(EPG_OPTIONAL_DATA)
    names.extend(GRAND_STAGE1_OPTIONAL_BRANCHES)
    return tuple(dict.fromkeys(names))


def grand_stage1_value(
    arrays: Dict[str, np.ndarray],
    features: Dict[str, np.ndarray],
    aliases: Sequence[str],
    transform: str,
    angle_unit: str,
) -> Optional[np.ndarray]:
    """Return one diagnostic variable, or None when the branch is unavailable."""
    if transform == "derived_Mx2_ep":
        return np.asarray(features["ep_missing_mass2"], dtype=float)
    #endif
    if transform == "derived_Mx2_epgamma":
        return np.asarray(features["pred_probe_mass2"], dtype=float)
    #endif

    branch = next(
        (name for name in aliases if name in arrays),
        None,
    )
    if branch is None:
        return None
    #endif

    values = np.asarray(arrays[branch], dtype=float)
    if transform == "lab_angle_deg":
        if angle_unit == "rad":
            values = np.degrees(values)
        #endif
    #endif
    return values


def grand_stage1_range(
    arrays_by_sample: Dict[str, Dict[str, np.ndarray]],
    features_by_sample: Dict[str, Dict[str, np.ndarray]],
    masks_by_sample: Dict[str, np.ndarray],
    aliases: Sequence[str],
    transform: str,
    angle_units: Dict[str, str],
) -> Optional[Tuple[float, float]]:
    """Common robust plotting range from all three samples."""
    pooled = []
    for sample_name in ("data", "pi0", "dvcs"):
        values = grand_stage1_value(
            arrays_by_sample[sample_name],
            features_by_sample[sample_name],
            aliases,
            transform,
            angle_units[sample_name],
        )
        if values is None:
            return None
        #endif
        mask = masks_by_sample[sample_name]
        local = values[mask & np.isfinite(values)]
        if local.size:
            # Bound the contribution of giant samples to percentile estimation.
            stride = max(1, local.size // 250_000)
            pooled.append(local[::stride])
        #endif
    #endfor

    if not pooled:
        return None
    #endif

    values = np.concatenate(pooled)
    if values.size == 0:
        return None
    #endif

    lo, hi = np.nanpercentile(values, [0.5, 99.5])
    if not np.isfinite(lo) or not np.isfinite(hi):
        return None
    #endif
    if hi <= lo:
        pad = max(1.0e-6, abs(lo) * 0.05 + 1.0e-6)
        lo -= pad
        hi += pad
    else:
        pad = 0.04 * (hi - lo)
        lo -= pad
        hi += pad
    #endif
    return float(lo), float(hi)


def grand_stage1_masks(
    feat: Dict[str, np.ndarray],
    region: str,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Minimal and nominal selections for the grand Stage-1 scan."""
    minimal = (
        np.asarray(feat["valid_tag"], dtype=bool)
        & stage2_region_mask(
            feat,
            region,
            ft_theta_max,
        )
        & np.isfinite(feat["electron_p"])
        & (feat["electron_p"] > 2.0)
        & np.isfinite(feat["theta_egamma_deg"])
        & (
            feat["theta_egamma_deg"]
            > THETA_EGAMMA_MIN_DEG
        )
        & np.isfinite(feat["pred_probe_energy"])
        & (feat["pred_probe_energy"] >= 0.40)
        & (feat["pred_probe_energy"] < max_probe_energy)
    )
    nominal = stage2_fit_mask(
        feat,
        region,
        ft_theta_max,
        0.40,
        max_probe_energy,
        mm2_min,
        mm2_max,
        probe_m2_max,
    )
    nominal &= (
        np.isfinite(feat["electron_p"])
        & (feat["electron_p"] > 2.0)
    )
    return minimal, nominal


def make_grand_stage1_diagnostics(
    period: PeriodConfig,
    arrays_by_sample: Dict[str, Dict[str, np.ndarray]],
    samples_by_name: Dict[str, EPGammaSample],
    features_by_sample: Dict[str, Dict[str, np.ndarray]],
    outdir: Path,
    *,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
) -> None:
    """
    Broad visual search for composition discriminators.

    For each variable group and detector region this writes:
      * unit-normalized 1D data/pi0/DVCS overlays under the NOMINAL
        exclusivity support;
      * E_gamma,tag correlations after the NOMINAL support, with one row per
        variable and columns = data, aaogen pi0, dvcsgen.

    All axes are linear.
    """
    outdir.mkdir(parents=True, exist_ok=True)

    sample_specs = (
        ("data", "data", "black"),
        ("pi0", r"$\pi^0$ MC", "tab:red"),
        ("dvcs", "BH/DVCS MC", "tab:blue"),
    )
    angle_units = {
        name: sample.angle_unit
        for name, sample in samples_by_name.items()
    }

    for region in ("FT", "FD_all"):
        minimal_masks = {}
        nominal_masks = {}
        for sample_name in ("data", "pi0", "dvcs"):
            minimal, nominal = grand_stage1_masks(
                features_by_sample[sample_name],
                region,
                ft_theta_max,
                max_probe_energy,
                mm2_min,
                mm2_max,
                probe_m2_max,
            )
            minimal_masks[sample_name] = minimal
            nominal_masks[sample_name] = nominal
        #endfor

        for group_name, group_specs in GRAND_STAGE1_GROUPS.items():
            # Keep only variables actually available in all three samples.
            available_specs = []
            ranges = {}
            for display_key, aliases, label, transform in group_specs:
                plot_range = grand_stage1_range(
                    arrays_by_sample,
                    features_by_sample,
                    minimal_masks,
                    aliases,
                    transform,
                    angle_units,
                )
                if plot_range is None:
                    continue
                #endif
                available_specs.append(
                    (display_key, aliases, label, transform)
                )
                ranges[display_key] = plot_range
            #endfor

            if not available_specs:
                continue
            #endif

            # ----------------------------------------------------------
            # 1D overlays after the nominal exclusivity support only.
            # The old "minimal" duplicates were intentionally removed.
            # ----------------------------------------------------------
            masks = nominal_masks
            nvars = len(available_specs)
            ncols = min(4, nvars)
            nrows = int(math.ceil(nvars / ncols))
            fig, axes = plt.subplots(
                nrows,
                ncols,
                figsize=(4.4 * ncols, 3.3 * nrows + 0.7),
                squeeze=False,
            )
            flat = axes.ravel()

            for ivar, (
                display_key,
                aliases,
                label,
                transform,
            ) in enumerate(available_specs):
                ax = flat[ivar]
                lo, hi = ranges[display_key]

                for sample_name, sample_label, color in sample_specs:
                    values = grand_stage1_value(
                        arrays_by_sample[sample_name],
                        features_by_sample[sample_name],
                        aliases,
                        transform,
                        angle_units[sample_name],
                    )
                    mask = (
                        masks[sample_name]
                        & np.isfinite(values)
                        & (values >= lo)
                        & (values < hi)
                    )
                    h, edges = np.histogram(
                        values[mask],
                        bins=80,
                        range=(lo, hi),
                    )
                    h = h.astype(float)
                    if np.sum(h) > 0.0:
                        h /= np.sum(h)
                    #endif
                    centers = 0.5 * (
                        edges[:-1] + edges[1:]
                    )
                    ax.step(
                        centers,
                        h,
                        where="mid",
                        linewidth=1.15,
                        color=color,
                        label=sample_label,
                    )
                #endfor

                ax.set_xlim(lo, hi)
                ax.set_ylim(bottom=0.0)
                ax.set_xlabel(label)
                ax.set_ylabel("unit-normalized")
                ax.grid(alpha=0.15)
                ax.set_title(display_key)
            #endfor

            for iax in range(
                len(available_specs),
                len(flat),
            ):
                flat[iax].set_axis_off()
            #endfor

            handles, labels = flat[0].get_legend_handles_labels()
            fig.legend(
                handles,
                labels,
                loc="upper center",
                ncol=3,
                frameon=False,
                bbox_to_anchor=(0.5, 0.935),
                fontsize=8.5,
            )
            fig.suptitle(
                f"{period.label}: {region} grand Stage-1 "
                f"{group_name.replace('_', ' ')} — nominal exclusivity",
                fontsize=12.5,
                y=0.99,
            )
            safe_finalize_figure(
                fig,
                outdir
                / (
                    f"grand_{group_name}_nominal_"
                    f"{period.key}_{region.lower()}.png"
                ),
                rect=(0.0, 0.0, 1.0, 0.88),
            )
            plt.close(fig)

            # ----------------------------------------------------------
            # 2D E_tag correlations after nominal exclusivity support.
            # One row per variable, columns = data/pi0/DVCS.
            # ----------------------------------------------------------
            nvars = len(available_specs)
            fig, axes = plt.subplots(
                nvars,
                3,
                figsize=(13.4, max(4.0, 2.55 * nvars + 0.8)),
                squeeze=False,
            )

            for irow, (
                display_key,
                aliases,
                label,
                transform,
            ) in enumerate(available_specs):
                lo, hi = ranges[display_key]

                # Shared z maximum across data/pi0/DVCS for this variable.
                histograms = {}
                zmax = 0.0
                for sample_name, _, _ in sample_specs:
                    values = grand_stage1_value(
                        arrays_by_sample[sample_name],
                        features_by_sample[sample_name],
                        aliases,
                        transform,
                        angle_units[sample_name],
                    )
                    feat = features_by_sample[sample_name]
                    etag = np.asarray(
                        feat["tag_energy"],
                        dtype=float,
                    )
                    mask = (
                        nominal_masks[sample_name]
                        & np.isfinite(etag)
                        & np.isfinite(values)
                        & (etag >= 0.4)
                        & (etag < 9.5)
                        & (values >= lo)
                        & (values < hi)
                    )
                    h2, xedges, yedges = np.histogram2d(
                        etag[mask],
                        values[mask],
                        bins=(60, 60),
                        range=((0.4, 9.5), (lo, hi)),
                    )
                    h2 = h2.astype(float)
                    if np.sum(h2) > 0.0:
                        h2 /= np.sum(h2)
                    #endif
                    histograms[sample_name] = (
                        h2,
                        xedges,
                        yedges,
                    )
                    zmax = max(
                        zmax,
                        float(np.max(h2))
                        if h2.size
                        else 0.0,
                    )
                #endfor

                for icol, (
                    sample_name,
                    sample_label,
                    _,
                ) in enumerate(sample_specs):
                    ax = axes[irow, icol]
                    h2, xedges, yedges = histograms[
                        sample_name
                    ]
                    mesh = ax.pcolormesh(
                        xedges,
                        yedges,
                        h2.T,
                        shading="auto",
                        vmin=0.0,
                        vmax=zmax if zmax > 0.0 else 1.0,
                    )
                    ax.set_xlim(0.4, 9.5)
                    ax.set_ylim(lo, hi)
                    ax.set_xlabel(
                        r"$E_{\gamma,\mathrm{tag}}$ (GeV)"
                    )
                    ax.set_ylabel(label)
                    if irow == 0:
                        ax.set_title(sample_label)
                    #endif
                #endfor
            #endfor

            fig.suptitle(
                f"{period.label}: {region} "
                f"$E_{{\\gamma,tag}}$ correlations — "
                f"{group_name.replace('_', ' ')}\n"
                "nominal exclusivity support",
                fontsize=12.5,
                y=0.995,
            )
            safe_finalize_figure(
                fig,
                outdir
                / (
                    f"grand_{group_name}_etag_correlations_"
                    f"{period.key}_{region.lower()}.png"
                ),
                rect=(0.0, 0.0, 1.0, 0.965),
            )
            plt.close(fig)
        #endfor
    #endfor




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


def subset_event_mapping(
    mapping: Dict[str, object],
    selector: np.ndarray,
    n_entries: int,
    *,
    preserve_metadata: bool,
) -> Dict[str, object]:
    """
    Safely subset a mixed event/metadata mapping.

    Event-length arrays are sliced by ``selector``.
    Scalar and other non-event metadata are either preserved verbatim or
    omitted according to ``preserve_metadata``.

    This is required for denominator feature stores because they intentionally
    contain scalar fixed-acceptance metadata alongside per-event arrays.
    """
    sel = np.asarray(selector)
    out: Dict[str, object] = {}

    for key, value in mapping.items():
        arr = np.asarray(value)

        if arr.ndim >= 1 and arr.shape[0] == int(n_entries):
            out[key] = np.asarray(arr[sel])
        elif preserve_metadata:
            # Keep true metadata exactly.  Scalars are converted back to their
            # Python/NumPy scalar form where possible so existing float(...)
            # consumers behave identically.
            if arr.ndim == 0:
                out[key] = arr.item()
            else:
                out[key] = value
            #endif
        #endif
    #endfor

    return out


def subset_feature_store(
    features: Dict[str, object],
    selector: np.ndarray,
    n_entries: int,
) -> Dict[str, object]:
    """
    Subset per-event feature arrays while preserving scalar detector metadata.
    """
    return subset_event_mapping(
        features,
        selector,
        n_entries,
        preserve_metadata=True,
    )


def subset_epgamma_sample(
    sample: EPGammaSample,
    selector: np.ndarray,
) -> EPGammaSample:
    """Return an EPGammaSample containing only selected entries."""
    sel = np.asarray(selector)
    return EPGammaSample(
        electron_p3=np.asarray(sample.electron_p3[sel]),
        proton_p3=np.asarray(sample.proton_p3[sel]),
        tag_p3=np.asarray(sample.tag_p3[sel]),
        tag_energy=np.asarray(sample.tag_energy[sel]),
        raw=subset_event_mapping(
            sample.raw,
            sel,
            len(sample.tag_energy),
            preserve_metadata=False,
        ),
        angle_unit=sample.angle_unit,
    )


def subset_eppi0_sample(
    sample: EPPi0Sample,
    selector: np.ndarray,
) -> EPPi0Sample:
    """Return an EPPi0Sample containing only selected entries."""
    sel = np.asarray(selector)
    return EPPi0Sample(
        electron_p3=np.asarray(sample.electron_p3[sel]),
        proton_p3=np.asarray(sample.proton_p3[sel]),
        pi0_p3=np.asarray(sample.pi0_p3[sel]),
        pi0_p=np.asarray(sample.pi0_p[sel]),
        pi0_mass=np.asarray(sample.pi0_mass[sel]),
        pi0_theta=np.asarray(sample.pi0_theta[sel]),
        raw=subset_event_mapping(
            sample.raw,
            sel,
            len(sample.pi0_mass),
            preserve_metadata=False,
        ),
        angle_unit=sample.angle_unit,
    )




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

    # Bulk reconstructed kinematics do not need float64 precision.
    # Converting here (rather than after feature construction) cuts the
    # memory traffic of the very large nSidis samples approximately in half.
    e_p = np.asarray(arrays["e_p"], dtype=np.float32)
    p_p = np.asarray(arrays["p1_p"], dtype=np.float32)
    g_e = np.asarray(arrays["p2_p"], dtype=np.float32)

    # Keep all non-vector branches, including the ROOT-tree missing-mass
    # observables.  Mx2/Mx2_1/Mx2_2 are required inputs, but they must remain
    # available in sample.raw because the denominator diagnostics use the
    # stored ROOT quantities directly.
    vector_only_branches = {
        "e_p", "e_theta", "e_phi",
        "p1_p", "p1_theta", "p1_phi",
        "p2_p", "p2_theta", "p2_phi",
    }
    raw = {
        key: np.asarray(value)
        for key, value in arrays.items()
        if key not in vector_only_branches
    }

    return EPGammaSample(
        electron_p3=cartesian_from_spherical(e_p, e_theta, e_phi),
        proton_p3=cartesian_from_spherical(p_p, p_theta, p_phi),
        tag_p3=cartesian_from_spherical(g_e, g_theta, g_phi),
        tag_energy=g_e,
        raw=raw,
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

    e_p = np.asarray(arrays["e_p"], dtype=np.float32)
    p_p = np.asarray(arrays["p1_p"], dtype=np.float32)
    pi_p = np.asarray(arrays["p2_p"], dtype=np.float32)
    pi_mass = np.asarray(arrays["Mh_gammagamma"], dtype=np.float32)

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
    empty_best_diagnostics = {
        "epg_index": np.asarray([], dtype=np.int64),
        "eppi0_index": np.asarray([], dtype=np.int64),
        "reco_probe_energy": np.asarray([], dtype=float),
        "pred_probe_energy": np.asarray([], dtype=float),
        "pred_probe_theta_deg": np.asarray([], dtype=float),
        "pi0_mass": np.asarray([], dtype=float),
        "reco_probe_mass2": np.asarray([], dtype=float),
        "probe_delta_theta_deg": np.asarray([], dtype=float),
        "probe_delta_E_over_E": np.asarray([], dtype=float),
        "passes_pred_consistency": np.asarray([], dtype=bool),
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
            best_diagnostics={
                key: value.copy()
                for key, value in empty_best_diagnostics.items()
            },
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
    with np.errstate(invalid="ignore", divide="ignore"):
        pred_theta_deg = np.degrees(
            np.arccos(
                np.clip(
                    pred3[:, 2] / pred_p,
                    -1.0,
                    1.0,
                )
            )
        )
    #endwith

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
        best_diag = {
            key: value.copy()
            for key, value in empty_best_diagnostics.items()
        }
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
            "pred_probe_theta_deg": np.asarray(
                pred_theta_deg[chosen], dtype=float
            ),
            "pi0_mass": np.asarray(pi_m[chosen], dtype=float),
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
    "theta_gg_1d",
    "mx2_ep_x_pTmiss",
    "mx2_ep_x_theta_gg",
)

# Production choice established by the Fa18 Inb driver study:
# M_X^2(ep) x pTmiss is the only nominal composition driver.
STAGE2_NOMINAL_DISCRIMINATOR = "mx2_ep_x_pTmiss"
STAGE2_PRODUCTION_MODEL = "morphed_mx2_ep_x_pTmiss"
STAGE2_DRIVER_DISCRIMINATORS: Tuple[str, ...] = (
    "mx2_ep_x_pTmiss",
)
STAGE2_VALIDATION_DISCRIMINATORS: Tuple[str, ...] = (
    "theta_gg_1d",
    "mx2_ep_x_theta_gg",
    "mx2_ep_1d",
)



def masked_finite_mean(values: np.ndarray, mask: np.ndarray) -> float:
    """Finite arithmetic mean of values selected by mask."""
    vv = np.asarray(values, dtype=float)
    mm = np.asarray(mask, dtype=bool)
    if vv.shape != mm.shape:
        raise ValueError(
            f"masked_finite_mean shape mismatch: values={vv.shape}, mask={mm.shape}"
        )
    #endif
    good = mm & np.isfinite(vv)
    if not np.any(good):
        return float("nan")
    #endif
    return float(np.mean(vv[good]))


def row_energy_coordinate(row: Dict[str, object]) -> float:
    """
    Plot binned energy points at the actual average predicted-probe energy.

    Legacy rows without an average-energy field fall back to the bin midpoint.
    """
    for key in (
        "mean_probe_energy_GeV",
        "data_mean_probe_energy_GeV",
        "plot_energy_GeV",
    ):
        try:
            value = float(row.get(key, float("nan")))
        except (TypeError, ValueError):
            continue
        #endtry
        if np.isfinite(value):
            return value
        #endif
    #endfor
    return 0.5 * (
        float(row["energy_low_GeV"])
        + float(row["energy_high_GeV"])
    )



def stage2_energy_edges(max_energy: float) -> np.ndarray:
    """
    Return Stage-II probe-energy edges, truncated at the supported data endpoint.
    """
    base = PROBE_ENERGY_EDGES[PROBE_ENERGY_EDGES < max_energy]
    edges = np.concatenate((base, np.asarray([max_energy], dtype=float)))
    edges = np.unique(edges)
    return edges[edges >= PROBE_ENERGY_EDGES[0]]


def nsidis_energy_edges_for_region(
    region: str,
    max_energy: float,
) -> np.ndarray:
    """
    Region-specific nSidis production binning.

    FD keeps the standard energy binning.
    FT merges the three lowest standard bins into 0.40-1.25 GeV because the
    low-energy FT numerator/denominator statistics are sparse and the user
    explicitly prefers one combined low-energy FT correction.
    """
    edges = stage2_energy_edges(max_energy)
    if str(region) != "FT":
        return edges
    #endif

    merged = [float(edges[0])]
    for edge in edges[1:]:
        edge = float(edge)

        # Merge the three low-energy bins into 0.40-1.25 GeV.
        if edge < 1.25 - 1.0e-12:
            continue
        #endif

        # Also merge 1.25-1.50 and 1.50-2.00 into 1.25-2.00 GeV.
        if 1.25 + 1.0e-12 < edge < 2.00 - 1.0e-12:
            continue
        #endif

        merged.append(edge)
    #endfor

    if merged[-1] < float(max_energy) - 1.0e-12:
        merged.append(float(max_energy))
    #endif
    return np.asarray(sorted(set(merged)), dtype=float)


@dataclass(frozen=True)
class PhotonAngularAcceptance:
    """
    Robust reconstructed-photon polar-angle support for one run period.

    The operational FT/FD boundaries are central-containment quantiles, not
    literal extrema.  Raw extrema are retained for diagnostics/provenance.
    """
    ft_theta_min_deg: float
    ft_theta_max_deg: float
    fd_theta_min_deg: float
    fd_theta_max_deg: float
    ft_count: int
    fd_count: int
    source: str
    containment: float
    quantile_low: float
    quantile_high: float
    ft_raw_min_deg: float
    ft_raw_max_deg: float
    fd_raw_min_deg: float
    fd_raw_max_deg: float
    ft_q001_deg: float
    ft_q999_deg: float
    fd_q001_deg: float
    fd_q999_deg: float

    def as_dict(self) -> Dict[str, object]:
        return {
            "ft_theta_min_deg": float(self.ft_theta_min_deg),
            "ft_theta_max_deg": float(self.ft_theta_max_deg),
            "fd_theta_min_deg": float(self.fd_theta_min_deg),
            "fd_theta_max_deg": float(self.fd_theta_max_deg),
            "ft_count": int(self.ft_count),
            "fd_count": int(self.fd_count),
            "source": str(self.source),
            "containment": float(self.containment),
            "quantile_low": float(self.quantile_low),
            "quantile_high": float(self.quantile_high),
            "ft_raw_min_deg": float(self.ft_raw_min_deg),
            "ft_raw_max_deg": float(self.ft_raw_max_deg),
            "fd_raw_min_deg": float(self.fd_raw_min_deg),
            "fd_raw_max_deg": float(self.fd_raw_max_deg),
            "ft_q001_deg": float(self.ft_q001_deg),
            "ft_q999_deg": float(self.ft_q999_deg),
            "fd_q001_deg": float(self.fd_q001_deg),
            "fd_q999_deg": float(self.fd_q999_deg),
        }


def photon_theta_deg_from_epgamma(sample: EPGammaSample) -> np.ndarray:
    """Polar angle of the actually reconstructed photon in an epgamma tree."""
    p3 = np.asarray(sample.tag_p3, dtype=float)
    p = np.sqrt(np.einsum("ij,ij->i", p3, p3))
    with np.errstate(invalid="ignore", divide="ignore"):
        c = p3[:, 2] / np.where(p > 0.0, p, np.nan)
        theta = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    #endwith
    return np.asarray(theta, dtype=float)


def infer_photon_angular_acceptance(
    sample: EPGammaSample,
    *,
    source: str,
    parent_mask: Optional[np.ndarray] = None,
    containment: float = 0.999,
) -> PhotonAngularAcceptance:
    """
    Fixed production photon angular acceptance.

      FT: 2.4 <= theta_gamma <= 5.0 deg
      FD: 6.0 <= theta_gamma <= 35.0 deg

    The arguments are retained for downstream interface compatibility only.
    No reconstructed-event extrema or percentiles determine these boundaries.
    """
    del sample, parent_mask, containment
    return PhotonAngularAcceptance(
        ft_theta_min_deg=2.4,
        ft_theta_max_deg=5.0,
        fd_theta_min_deg=6.0,
        fd_theta_max_deg=35.0,
        ft_count=0,
        fd_count=0,
        source=f"{source}; fixed production angular acceptance",
        containment=1.0,
        quantile_low=0.0,
        quantile_high=1.0,
        ft_raw_min_deg=2.4,
        ft_raw_max_deg=5.0,
        fd_raw_min_deg=6.0,
        fd_raw_max_deg=35.0,
        ft_q001_deg=2.4,
        ft_q999_deg=5.0,
        fd_q001_deg=6.0,
        fd_q999_deg=35.0,
    )



def attach_photon_angular_acceptance(
    features: Dict[str, np.ndarray],
    acceptance: PhotonAngularAcceptance,
) -> None:
    """
    Attach scalar angular-support metadata to a denominator feature store.

    Keeping the metadata inside the feature dict lets every existing
    stage2_fit_mask() caller automatically use the same period-specific
    acceptance without a large and error-prone signature refactor.
    """
    features["_ft_theta_min_deg"] = float(acceptance.ft_theta_min_deg)
    features["_ft_theta_max_deg"] = float(acceptance.ft_theta_max_deg)
    features["_fd_theta_min_deg"] = float(acceptance.fd_theta_min_deg)
    features["_fd_theta_max_deg"] = float(acceptance.fd_theta_max_deg)


def feature_angular_acceptance(
    features: Dict[str, np.ndarray],
    fallback_ft_theta_max: float,
) -> Tuple[float, float, float, float]:
    """Return attached period-specific support, or the legacy fallback."""
    if all(
        key in features
        for key in (
            "_ft_theta_min_deg",
            "_ft_theta_max_deg",
            "_fd_theta_min_deg",
            "_fd_theta_max_deg",
        )
    ):
        return (
            float(features["_ft_theta_min_deg"]),
            float(features["_ft_theta_max_deg"]),
            float(features["_fd_theta_min_deg"]),
            float(features["_fd_theta_max_deg"]),
        )
    #endif

    # Legacy fallback for internal/self-test paths that do not attach data-derived
    # acceptance. Normal production/nSidis processing attaches exact boundaries.
    return (0.0, float(fallback_ft_theta_max), float(fallback_ft_theta_max), 180.0)


def electron_proton_opening_angle_deg(
    electron_p3: np.ndarray,
    proton_p3: np.ndarray,
) -> np.ndarray:
    """Lab opening angle between two supplied three-momentum vectors."""
    e3 = np.asarray(electron_p3, dtype=float)
    p3 = np.asarray(proton_p3, dtype=float)
    en = np.sqrt(np.einsum("ij,ij->i", e3, e3))
    pn = np.sqrt(np.einsum("ij,ij->i", p3, p3))
    dot = np.einsum("ij,ij->i", e3, p3)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = dot / np.where((en > 0.0) & (pn > 0.0), en * pn, np.nan)
        theta = np.degrees(np.arccos(np.clip(c, -1.0, 1.0)))
    #endwith
    return np.asarray(theta, dtype=float)


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
) -> Dict[str, np.ndarray]:
    """
    Build denominator kinematics for real epgamma events. Missing-mass
    observables are read directly from the ROOT-tree Mx2/Mx2_1/Mx2_2 branches;
    no mixed-event/synthetic-photon construction is performed.
    """
    e3 = epg.electron_p3
    p3 = epg.proton_p3
    tag_source = epg.tag_p3
    n = len(epg.tag_energy)
    tag_E = epg.tag_energy

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

    # Use the missing-mass quantities already stored in the ROOT tree.
    # The processing convention is fixed: p1 = proton and p2 = photon.
    # Therefore ThreeParticles stores:
    #   Mx2_1 = MM^2(epX)
    #   Mx2_2 = MM^2(e gamma X)
    #   Mx2   = MM^2(ep gamma X)
    epmiss_m2 = np.asarray(epg.raw["Mx2_1"], dtype=float)
    egamma_miss_m2 = np.asarray(epg.raw["Mx2_2"], dtype=float)
    epgamma_miss_m2 = np.asarray(epg.raw["Mx2"], dtype=float)

    # Use the stored Delta_phi branch for the proton-photon coplanarity.
    # The tree value is centered near pi, so convert it to the signed residual
    # about 180 degrees in [-pi, pi).
    if "Delta_phi" in epg.raw:
        raw_delta_phi_tree = np.asarray(epg.raw["Delta_phi"], dtype=float)
        dphi_pgamma = (raw_delta_phi_tree - math.pi + math.pi) % TWO_PI - math.pi
        dphi_pgamma_deg = np.degrees(dphi_pgamma)
    else:
        proton_phi = np.arctan2(p3[:, 1], p3[:, 0])
        tag_phi = np.arctan2(tag_source[:, 1], tag_source[:, 0])
        raw_pair_delta_phi = proton_phi - tag_phi
        dphi_pgamma = (raw_pair_delta_phi - math.pi + math.pi) % TWO_PI - math.pi
        dphi_pgamma_deg = np.degrees(dphi_pgamma)
    #endif

    px -= tag_source[:, 0]
    py -= tag_source[:, 1]
    pz -= tag_source[:, 2]

    pred_p2 = px * px + py * py + pz * pz
    pred_p = np.sqrt(pred_p2)
    pred_m2 = pred_E * pred_E - pred_p2

    pred_p3 = np.column_stack((px, py, pz))

    # Angle(gamma,X): X is the inclusive hadronic system recoiling against
    # the scattered electron,
    #   P_X = P_beam + P_target - P_e.
    # Therefore X contains the recoil proton and the complete produced
    # hadronic/photon final state.  Compare its three-momentum direction to
    # the reconstructed tag-photon direction.
    ex_p3 = np.column_stack((
        -e3[:, 0],
        -e3[:, 1],
        beam_p - e3[:, 2],
    ))
    theta_gammaX_deg = electron_proton_opening_angle_deg(tag_source, ex_p3)

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
        "egamma_missing_mass2": np.asarray(egamma_miss_m2, dtype=np.float32),
        "epgamma_missing_mass2": np.asarray(epgamma_miss_m2, dtype=np.float32),
        "delta_phi_pgamma_deg": np.asarray(dphi_pgamma_deg, dtype=np.float32),
        "theta_gammaX_deg": np.asarray(theta_gammaX_deg, dtype=np.float32),
        "tag_energy": np.asarray(tag_E, dtype=np.float32),
        "pred_probe_energy": np.asarray(pred_E, dtype=np.float32),
        "pred_probe_mass2": np.asarray(pred_m2, dtype=np.float32),
        "pred_probe_theta_deg": np.asarray(pred_theta_deg, dtype=np.float32),
        "pred_probe_sector": np.asarray(sector, dtype=np.int16),
        "electron_p": np.asarray(
            np.sqrt(np.einsum("ij,ij->i", e3, e3)),
            dtype=np.float32,
        ),
        "theta_ep_deg": np.asarray(
            electron_proton_opening_angle_deg(e3, p3),
            dtype=np.float32,
        ),
        "theta_egamma_deg": np.asarray(
            electron_proton_opening_angle_deg(e3, tag_source),
            dtype=np.float32,
        ),
        "valid_tag": np.asarray(valid, dtype=bool),
    }

    # pTmiss(epgamma) is fully reconstructable from the measured four-vectors.
    # The loose nSidis trees are not required to store the old skim-level branch,
    # so always compute it and use it as a fallback.
    recomputed_ptmiss = np.asarray(
        np.sqrt(px * px + py * py),
        dtype=np.float32,
    )
    out["recomputed_pTmiss"] = recomputed_ptmiss

    if True:
        if "pTmiss" in epg.raw:
            out["stored_pTmiss"] = np.asarray(
                epg.raw["pTmiss"], dtype=np.float32
            )
            out["pTmiss_recomputed_fallback"] = np.zeros(n, dtype=bool)
        else:
            out["stored_pTmiss"] = recomputed_ptmiss
            out["pTmiss_recomputed_fallback"] = np.ones(n, dtype=bool)
        #endif

        if "theta_gamma_gamma" in epg.raw:
            out["stored_theta_gamma_gamma"] = np.asarray(
                epg.raw["theta_gamma_gamma"], dtype=np.float32
            )
        #endif
        if "Delta_phi" in epg.raw:
            raw_delta_phi = np.asarray(
                epg.raw["Delta_phi"], dtype=float
            )
            wrapped_delta_phi = np.abs(
                (raw_delta_phi + math.pi) % (2.0 * math.pi)
                - math.pi
            )
            out["stored_delta_phi_abs_rad"] = np.asarray(
                wrapped_delta_phi,
                dtype=np.float32,
            )
            out["stored_delta_phi_rad"] = np.asarray(
                np.mod(raw_delta_phi, 2.0 * math.pi),
                dtype=np.float32,
            )
        #endif
        if "Emiss2" in epg.raw:
            out["stored_Emiss2"] = np.asarray(
                epg.raw["Emiss2"],
                dtype=np.float32,
            )
        #endif
        if "xF2" in epg.raw:
            out["stored_xF2"] = np.asarray(
                epg.raw["xF2"],
                dtype=np.float32,
            )
        #endif
        if "pT" in epg.raw:
            out["stored_sidis_pT"] = np.asarray(
                epg.raw["pT"],
                dtype=np.float32,
            )
        #endif
    else:
        # The original stored pTmiss is invalid after mixed-event substitution.
        out["stored_pTmiss"] = recomputed_ptmiss
        out["pTmiss_recomputed_fallback"] = np.ones(n, dtype=bool)
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
    """
    Predicted-probe detector-region mask.

    Normal processing attaches fixed FT/FD polar-angle acceptance to every
    feature store.  This excludes the FT beam hole and the angular gap/outside
    acceptance instead of treating every theta <= 5.5 deg as reconstructable.
    """
    theta = features["pred_probe_theta_deg"]
    sector = features["pred_probe_sector"]
    ft_min, ft_max, fd_min, fd_max = feature_angular_acceptance(
        features, ft_theta_max
    )

    finite = np.isfinite(theta)
    ft = finite & (theta >= ft_min) & (theta <= ft_max)
    fd = (
        finite
        & (theta >= fd_min)
        & (theta <= fd_max)
        & (sector >= 1)
        & (sector <= 6)
    )

    if region_name == "all":
        return ft | fd
    #endif
    if region_name == "FT":
        return ft
    #endif
    if region_name == "FD_all":
        return fd
    #endif
    if region_name.startswith("FD_S"):
        s = int(region_name[-1])
        return fd & (sector == s)
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
    """
    Common denominator/template selection.

    The nominal exclusivity requirement is now the M_X^2(epX) window only.
    The lower edge is supplied by ``mm2_min`` and the upper edge is the
    period-specific CLASDIS FD valley between the missing-pi0 and missing-eta
    peaks.  That same upper edge is imposed in FT and FD.

    The previously tested Valerii-labelled MM^2(e gamma X), Delta phi, and
    Angle(gamma,X) requirements are not imposed because the labels/cut
    definitions in the presentation are not sufficiently self-consistent for
    production use.  Those quantities remain diagnostic variables.
    """
    theta_egamma = feat.get("theta_egamma_deg")
    if theta_egamma is None:
        theta_egamma_mask = np.ones(
            len(feat["pred_probe_energy"]), dtype=bool
        )
    else:
        theta_egamma_mask = (
            np.isfinite(theta_egamma)
            & (theta_egamma > THETA_EGAMMA_MIN_DEG)
        )
    #endif

    return (
        feat["valid_tag"]
        & theta_egamma_mask
        & stage2_region_mask(feat, region_name, ft_theta_max)
        & (feat["pred_probe_energy"] >= elo)
        & (feat["pred_probe_energy"] < ehi)
        & np.isfinite(feat["ep_missing_mass2"])
        & (feat["ep_missing_mass2"] >= mm2_min)
        & (feat["ep_missing_mass2"] < mm2_max)
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
    if discriminator in ("theta_gg_1d", "mx2_ep_x_theta_gg"):
        return all(
            "stored_theta_gamma_gamma" in f
            for f in (data_f, pi0_f, dvcs_f)
        )
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

    if discriminator == "theta_gg_1d":
        local = mask & np.isfinite(feat["stored_theta_gamma_gamma"])
        h, _ = np.histogram(
            feat["stored_theta_gamma_gamma"][local],
            bins=np.linspace(0.0, theta_max, theta_bins + 1),
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



def poisson_deviance_value(
    data: np.ndarray,
    mu: np.ndarray,
) -> float:
    """Raw binned Poisson deviance without attaching an ndof interpretation."""
    data = np.asarray(data, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1.0e-12, None)
    positive = data > 0.0
    terms = np.zeros_like(data, dtype=float)
    terms[positive] = data[positive] * np.log(data[positive] / mu[positive])
    return 2.0 * float(np.sum(mu - data + terms))


def binned_model_quality(
    data: np.ndarray,
    mu: np.ndarray,
    *,
    active_mu_threshold: float = 1.0,
    pearson_mu_threshold: float = 5.0,
) -> Dict[str, float]:
    """
    Transparent descriptive goodness metrics for a histogram/model pair.

    The old code reported Poisson deviance / ndof after counting essentially
    every 2D cell with a tiny nonzero model expectation as a degree of freedom.
    With finely binned sparse 2D histograms this can make a visibly poor
    projection appear to have D/ndof ~ 1.  The deviance itself is valid, but
    that reduced-number interpretation is not robust.

    Here:
      active bins = data > 0 OR model expectation >= active_mu_threshold;
      D/active-bin is a descriptive sparse-histogram metric;
      Pearson chi2/bin is reported only where model expectation >= 5.

    These are diagnostics, not p-values.
    """
    d = np.asarray(data, dtype=float)
    m = np.clip(np.asarray(mu, dtype=float), 1.0e-12, None)
    if d.shape != m.shape:
        raise ValueError("data/model shapes differ in binned_model_quality")
    #endif

    dev = poisson_deviance_value(d, m)
    active = (d > 0.0) | (m >= float(active_mu_threshold))
    n_active = int(np.count_nonzero(active))
    dev_active = poisson_deviance_value(d[active], m[active]) if n_active else float("nan")

    pearson_mask = m >= float(pearson_mu_threshold)
    n_pearson = int(np.count_nonzero(pearson_mask))
    if n_pearson:
        chi2 = float(np.sum(
            (d[pearson_mask] - m[pearson_mask]) ** 2 / m[pearson_mask]
        ))
    else:
        chi2 = float("nan")
    #endif

    return {
        "poisson_deviance": float(dev),
        "active_bins": n_active,
        "poisson_deviance_active": float(dev_active),
        "deviance_per_active_bin": (
            float(dev_active / n_active) if n_active else float("nan")
        ),
        "pearson_chi2_mu_ge5": float(chi2),
        "pearson_bins_mu_ge5": n_pearson,
        "pearson_chi2_per_bin_mu_ge5": (
            float(chi2 / n_pearson) if n_pearson else float("nan")
        ),
    }


def fit_quality_for_morphed_2d_model(
    data_hist: np.ndarray,
    pi0_hist: np.ndarray,
    dvcs_hist: np.ndarray,
    fit: SharedMorphedFitResult,
    discriminator: str,
) -> Dict[str, float]:
    """
    Compute 2D and both 1D-projection quality metrics for the ACTUAL fitted
    morphed model.  These metrics correspond directly to the plotted curves.
    """
    if not fit.success or not fit.nuisance:
        return {}
    #endif

    hp = np.asarray(pi0_hist, dtype=float)
    hv = np.asarray(dvcs_hist, dtype=float)
    hd = np.asarray(data_hist, dtype=float)

    tp = morph_template_second_axis(
        hp,
        float(fit.nuisance.get(f"{discriminator}_pi0_shift_bins", 0.0)),
        float(fit.nuisance.get(f"{discriminator}_pi0_sigma_bins", 0.0)),
    )
    td = morph_template_second_axis(
        hv,
        float(fit.nuisance.get(f"{discriminator}_dvcs_shift_bins", 0.0)),
        float(fit.nuisance.get(f"{discriminator}_dvcs_sigma_bins", 0.0)),
    )

    ndata = float(np.sum(hd))
    fpi0 = float(fit.pi0_fraction)
    mu2d = ndata * (fpi0 * tp + (1.0 - fpi0) * td)

    q2 = binned_model_quality(hd, mu2d)
    qx = binned_model_quality(np.sum(hd, axis=1), np.sum(mu2d, axis=1))
    qy = binned_model_quality(np.sum(hd, axis=0), np.sum(mu2d, axis=0))

    out: Dict[str, float] = {}
    for prefix, q in (
        ("quality_2d", q2),
        ("quality_mx2_projection", qx),
        ("quality_ptmiss_projection", qy),
    ):
        for key, value in q.items():
            out[f"{prefix}_{key}"] = value
        #endfor
    #endfor
    return out


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
    driver_names: Optional[Sequence[str]] = None,
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
    requested_drivers = (
        tuple(STAGE2_DRIVER_DISCRIMINATORS)
        if driver_names is None
        else tuple(str(x) for x in driver_names)
    )
    names = [n for n in requested_drivers if n in histograms]
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

    # The coordinate-profile sweeps revisit the same discrete morph states
    # many times. Cache each component template by its exact nuisance pair.
    morph_cache: Dict[
        Tuple[str, str, float, float],
        np.ndarray,
    ] = {}

    def cached_morph(
        name: str,
        component: str,
        hist: np.ndarray,
        shift_bins: float,
        sigma_bins: float,
    ) -> np.ndarray:
        key = (
            str(name),
            str(component),
            round(float(shift_bins), 12),
            round(float(sigma_bins), 12),
        )
        cached = morph_cache.get(key)
        if cached is not None:
            return cached
        #endif

        value = morph_template_second_axis(
            hist,
            float(shift_bins),
            float(sigma_bins),
        )
        morph_cache[key] = value
        return value
    #enddef

    def morphed_templates(
        state: Dict[str, float],
    ) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
        tp: Dict[str, np.ndarray] = {}
        td: Dict[str, np.ndarray] = {}
        for name in names:
            _, hp, hv = histograms[name]

            pi0_shift = state[f"{name}_pi0_shift_bins"]
            pi0_sigma = state[f"{name}_pi0_sigma_bins"]
            dvcs_shift = state[f"{name}_dvcs_shift_bins"]
            dvcs_sigma = state[f"{name}_dvcs_sigma_bins"]

            tp[name] = cached_morph(
                name,
                "pi0",
                hp,
                pi0_shift,
                pi0_sigma,
            )
            td[name] = cached_morph(
                name,
                "dvcs",
                hv,
                dvcs_shift,
                dvcs_sigma,
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

    if not np.isfinite(f_final):
        return SharedMorphedFitResult(
            False,
            "non-finite fitted pi0 fraction after profiling/refinement",
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


def run_stage2_shared_morphed_fits(period: PeriodConfig, data_f: Dict[str,np.ndarray], pi0_f: Dict[str,np.ndarray],
                                   dvcs_f: Dict[str,np.ndarray], pi0_control: Dict[str,Dict[str,float]],
                                   ft_theta_max: float, max_probe_energy: float, mm2_min: float, mm2_max: float,
                                   probe_m2_max: float, mm2_bins_2d: int, probe_m2_bins_2d: int,
                                   ptmiss_max: float, ptmiss_bins: int, theta_max: float, theta_bins: int,
                                   min_data_count: int, min_template_count: int, nuisance_shift_prior: float,
                                   nuisance_sigma_prior: float, max_shift_bins: float, max_sigma_bins: float,
                                   closure_truth_fractions: Sequence[float],
                                   run_closure: bool = True) -> Tuple[List[Dict[str,object]], List[Dict[str,object]]]:
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
                 "fit_model":STAGE2_PRODUCTION_MODEL,
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
            if run_closure and region in ("FT", "FD_all") and hist and fit.success:
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
    Profile the shared pi0 fraction for fixed templates.

    The extended-Poisson objective is convex in the single mixture fraction f.
    Instead of repeatedly evaluating it through minimize_scalar, solve the
    monotonic score equation exactly with Brent's method.  This is materially
    faster because this function is called thousands of times during morph
    coordinate profiling.

    The previous bounded minimize_scalar implementation remains as a numerical
    fallback for pathological inputs.  The objective, allowed fraction range,
    and returned NLL are unchanged.
    """
    # Use exactly the histogram keys supplied by the calling fit.  Earlier
    # versions intersected with STAGE2_DRIVER_DISCRIMINATORS, which silently
    # discarded diagnostic-only/custom driver names and produced NaN fractions
    # despite fit_success=True.
    names = [
        str(name)
        for name in data_hists
        if name in pi0_templates and name in dvcs_templates
    ]
    if not names:
        return float("nan"), float("nan")
    #endif

    f_lo = 3.0e-4
    f_hi = 1.0 - 3.0e-4

    prepared = []
    for name in names:
        data = np.asarray(data_hists[name], dtype=float).ravel()
        tp = np.asarray(pi0_templates[name], dtype=float).ravel()
        td = np.asarray(dvcs_templates[name], dtype=float).ravel()
        nd = float(np.sum(data))
        delta = tp - td
        prepared.append((data, tp, td, delta, nd))
    #endfor

    def objective_f(f: float) -> float:
        ff = float(np.clip(f, f_lo, f_hi))
        total = 0.0
        for data, tp, td, delta, nd in prepared:
            shape = td + ff * delta
            mu = np.clip(nd * shape, 1.0e-12, None)
            total += float(np.sum(mu - data * np.log(mu))) / len(prepared)
        #endfor
        return float(total)
    #enddef

    def score_f(f: float) -> float:
        """
        d(NLL)/df.  For mu_i = N [td_i + f (tp_i-td_i)],

            dNLL/df = sum_i N delta_i [1 - data_i/mu_i].

        The derivative is monotonic nondecreasing because the second derivative
        is positive semidefinite, so a bracketed root is the unique interior
        minimum.
        """
        ff = float(np.clip(f, f_lo, f_hi))
        total = 0.0
        for data, tp, td, delta, nd in prepared:
            shape = td + ff * delta
            mu = np.clip(nd * shape, 1.0e-12, None)
            total += float(
                np.sum(nd * delta * (1.0 - data / mu))
            ) / len(prepared)
        #endfor
        return float(total)
    #enddef

    try:
        score_lo = score_f(f_lo)
        score_hi = score_f(f_hi)

        if not np.isfinite(score_lo) or not np.isfinite(score_hi):
            raise FloatingPointError("non-finite profile score")
        #endif

        # Convex objective: if the derivative already has one sign across the
        # interval, the minimum is at the corresponding boundary.
        if score_lo >= 0.0:
            f_best = f_lo
        elif score_hi <= 0.0:
            f_best = f_hi
        else:
            f_best = float(
                brentq(
                    score_f,
                    f_lo,
                    f_hi,
                    xtol=1.0e-10,
                    rtol=1.0e-12,
                    maxiter=80,
                )
            )
        #endif

        return float(f_best), objective_f(f_best)

    except Exception:
        # Conservative fallback: preserve the prior implementation if the
        # score equation is numerically pathological for any unexpected input.
        result = minimize_scalar(
            objective_f,
            bounds=(f_lo, f_hi),
            method="bounded",
            options={"xatol": 1.0e-7, "maxiter": 120},
        )
        return float(result.x), float(result.fun)
    #endtry



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




SELECTION_CANDIDATE_MODELS: Tuple[Tuple[str, Tuple[str, ...]], ...] = (
    ("mx2_x_ptmiss", ("mx2_ep_x_pTmiss",)),
    ("mx2_x_theta_gg", ("mx2_ep_x_theta_gg",)),
    (
        "shared_ptmiss_theta",
        ("mx2_ep_x_pTmiss", "mx2_ep_x_theta_gg"),
    ),
)


def run_stage2_candidate_model_study(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    support_probe_m2_max: float,
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
) -> List[Dict[str, object]]:
    """
    Region-specific denominator-model selection study.

    Candidate models are evaluated only for FT and FD_all.  They are
    DIAGNOSTIC and do not feed Stage III in v045.

    Crucially, M_probe^2 is not one of the candidate production drivers here.
    It remains available as a high-resolution diagnostic.
    """
    rows: List[Dict[str, object]] = []
    edges = stage2_energy_edges(max_probe_energy)

    for region in ("FT", "FD_all"):
        for ib in range(len(edges) - 1):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])

            masks = {
                key: stage2_fit_mask(
                    feat,
                    region,
                    ft_theta_max,
                    elo,
                    ehi,
                    mm2_min,
                    mm2_max,
                    support_probe_m2_max,
                )
                for key, feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                )
            }

            for model_name, drivers in SELECTION_CANDIDATE_MODELS:
                hist: Dict[
                    str,
                    Tuple[np.ndarray, np.ndarray, np.ndarray],
                ] = {}

                for disc in drivers:
                    if not discriminator_available(
                        disc, data_f, pi0_f, dvcs_f
                    ):
                        continue
                    #endif

                    hd = histogram_for_discriminator(
                        disc,
                        data_f,
                        masks["data"],
                        mm2_min=mm2_min,
                        mm2_max=mm2_max,
                        probe_m2_max=support_probe_m2_max,
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
                        probe_m2_max=support_probe_m2_max,
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
                        probe_m2_max=support_probe_m2_max,
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

                row: Dict[str, object] = {
                    "period": period.key,
                    "label": period.label,
                    "region": region,
                    "energy_low_GeV": elo,
                    "energy_high_GeV": ehi,
                    "energy_center_GeV": 0.5 * (elo + ehi),
                    "candidate_model": model_name,
                    "drivers": "+".join(drivers),
                    "candidate_is_production": 0,
                }

                if len(hist) != len(drivers):
                    row.update({
                        "fit_success": 0,
                        "fit_message": "required driver unavailable or insufficient statistics",
                        "pi0_fraction": float("nan"),
                        "deviance_per_ndof": float("nan"),
                        "poisson_deviance": float("nan"),
                        "ndof": 0,
                    })
                    rows.append(row)
                    continue
                #endif

                fit = fit_shared_morphed_composition(
                    hist,
                    {},
                    nuisance_shift_prior,
                    nuisance_sigma_prior,
                    max_shift_bins,
                    max_sigma_bins,
                    driver_names=drivers,
                )
                row.update({
                    "fit_success": int(fit.success),
                    "fit_message": fit.message,
                    "pi0_fraction": float(fit.pi0_fraction),
                    "dvcs_fraction": (
                        1.0 - float(fit.pi0_fraction)
                        if np.isfinite(fit.pi0_fraction)
                        else float("nan")
                    ),
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
            #endfor
        #endfor
    #endfor

    return rows


def _component_projection_from_independent_fit(
    hd: np.ndarray,
    hp: np.ndarray,
    hv: np.ndarray,
    fit: TwoComponentFitResult,
    axis: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Project data and fitted raw-template components onto one 2D axis."""
    data = np.sum(hd, axis=axis)
    pshape = np.sum(hp, axis=axis)
    dshape = np.sum(hv, axis=axis)

    pshape = (
        pshape / np.sum(pshape)
        if np.sum(pshape) > 0 else pshape
    )
    dshape = (
        dshape / np.sum(dshape)
        if np.sum(dshape) > 0 else dshape
    )

    pi0_c = float(fit.pi0_yield) * pshape
    dvcs_c = float(fit.dvcs_yield) * dshape
    return data, pi0_c, dvcs_c, pi0_c + dvcs_c


def make_selection_driver_fit_canvases(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    outdir: Path,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    support_probe_m2_max: float,
    mm2_bins: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
    min_data_count: int,
    min_template_count: int,
) -> None:
    """
    Show the shapes that actually distinguish the two components.

    For every nominal energy bin and for FT/FD_all separately:
      column 1: M_X^2(ep) projection from the M_X^2 x pTmiss fit
      column 2: pTmiss projection from that same fit
      column 3: M_X^2(ep) projection from the M_X^2 x theta_gg fit
      column 4: theta_gg projection from that same fit

    This directly addresses the fact that previous canvases hid M_X^2 by
    projecting only onto the second axis.
    """
    edges = stage2_energy_edges(max_probe_energy)

    for region in ("FT", "FD_all"):
        nrows = len(edges) - 1
        fig, axes = plt.subplots(
            nrows,
            4,
            figsize=(19.0, max(9.0, 2.65 * nrows)),
            squeeze=False,
        )


        for ib in range(nrows):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])

            masks = {
                key: stage2_fit_mask(
                    feat,
                    region,
                    ft_theta_max,
                    elo,
                    ehi,
                    mm2_min,
                    mm2_max,
                    support_probe_m2_max,
                )
                for key, feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                )
            }

            for pair_index, disc in enumerate(
                ("mx2_ep_x_pTmiss", "mx2_ep_x_theta_gg")
            ):
                hd = histogram_for_discriminator(
                    disc, data_f, masks["data"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    probe_m2_max=support_probe_m2_max,
                    mm2_bins_2d=mm2_bins,
                    probe_m2_bins_2d=120,
                    bins_1d=90,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )
                hp = histogram_for_discriminator(
                    disc, pi0_f, masks["pi0"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    probe_m2_max=support_probe_m2_max,
                    mm2_bins_2d=mm2_bins,
                    probe_m2_bins_2d=120,
                    bins_1d=90,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )
                hv = histogram_for_discriminator(
                    disc, dvcs_f, masks["dvcs"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    probe_m2_max=support_probe_m2_max,
                    mm2_bins_2d=mm2_bins,
                    probe_m2_bins_2d=120,
                    bins_1d=90,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )

                fit = fit_two_component_poisson(hd, hp, hv)
                ax_mx = axes[ib, 2 * pair_index]
                ax_second = axes[ib, 2 * pair_index + 1]

                if not fit.success:
                    for ax in (ax_mx, ax_second):
                        ax.text(
                            0.5, 0.5, "fit unavailable",
                            transform=ax.transAxes,
                            ha="center", va="center",
                        )
                    #endfor
                    continue
                #endif

                # Mx2 is axis 0, so sum over axis 1.
                d_mx, p_mx, v_mx, m_mx = (
                    _component_projection_from_independent_fit(
                        hd, hp, hv, fit, axis=1
                    )
                )
                mx_edges = np.linspace(mm2_min, mm2_max, mm2_bins + 1)
                mx_centers = 0.5 * (mx_edges[:-1] + mx_edges[1:])

                # Second coordinate is axis 1, so sum over axis 0.
                d_2, p_2, v_2, m_2 = (
                    _component_projection_from_independent_fit(
                        hd, hp, hv, fit, axis=0
                    )
                )

                if disc == "mx2_ep_x_pTmiss":
                    second_edges = np.linspace(
                        0.0, ptmiss_max, ptmiss_bins + 1
                    )
                    second_label = r"$p_{T,\mathrm{miss}}$ (GeV)"
                    second_xlim = (
                        (0.0, 0.30)
                        if region == "FT"
                        else (0.0, ptmiss_max)
                    )
                    model_label = r"$M_X^2\otimes p_{T,\mathrm{miss}}$"
                else:
                    second_edges = np.linspace(
                        0.0, theta_max, theta_bins + 1
                    )
                    second_label = r"$\theta_{\gamma\gamma}$ (deg)"
                    second_xlim = (0.0, theta_max)
                    model_label = r"$M_X^2\otimes\theta_{\gamma\gamma}$"
                #endif
                second_centers = 0.5 * (
                    second_edges[:-1] + second_edges[1:]
                )

                for ax, xx, dd, pp, vv, mm, xlabel in (
                    (
                        ax_mx, mx_centers, d_mx, p_mx, v_mx, m_mx,
                        r"$M_X^2(ep)$ (GeV$^2$)",
                    ),
                    (
                        ax_second, second_centers, d_2, p_2, v_2, m_2,
                        second_label,
                    ),
                ):
                    ax.errorbar(
                        xx, dd,
                        yerr=np.sqrt(np.maximum(dd, 1.0)),
                        fmt="o", ms=2.0, linewidth=0.6,
                        color=SAMPLE_COLORS["data"],
                        label="data",
                    )
                    ax.step(
                        xx, vv, where="mid",
                        color=SAMPLE_COLORS["dvcs_mc"],
                        linewidth=1.1,
                        label="BH/DVCS",
                    )
                    ax.step(
                        xx, pp, where="mid",
                        color=SAMPLE_COLORS["pi0_mc"],
                        linewidth=1.1,
                        label=r"$\pi^0$",
                    )
                    ax.step(
                        xx, mm, where="mid",
                        color=SAMPLE_COLORS["fit"],
                        linewidth=1.4,
                        label="total fit",
                    )
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel("entries / bin")
                    ax.grid(alpha=0.18)
                #endfor

                ax_second.set_xlim(*second_xlim)
                ax_mx.set_title(
                    f"{elo:.2f}-{ehi:.2f} GeV; {model_label}\n"
                    + rf"$f_{{\pi^0}}={fit.pi0_fraction:.3f}$, "
                    + rf"$D/\mathrm{{ndof}}="
                    + (
                        f"{fit.poisson_deviance/fit.ndof:.2f}"
                        if fit.ndof else "nan"
                    ),
                    fontsize=8.0,
                )
                ax_second.set_title("same fit: second-axis projection", fontsize=8.0)
            #endfor

        #endfor

        axes[0, 0].legend(fontsize=6.5, frameon=False, ncol=2)
        fig.suptitle(
            f"{period.label}: discriminator-selection ACTUAL fit projections — "
            f"{region}\n"
            r"$M_X^2(ep)$ is shown explicitly in columns 1 and 3",
            fontsize=13,
        )
        safe_finalize_figure(
            fig,
            outdir / f"canvas_selection_fit_projections_{region.lower()}.png",
            rect=(0, 0, 1, 0.965),
        )
        plt.close(fig)


    #endfor


def make_candidate_model_summary_canvas(
    period: PeriodConfig,
    rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """Compact comparison of the three candidate production models."""
    if not rows:
        return
    #endif

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 9.0))
    labels = {
        "mx2_x_ptmiss": r"$M_X^2\otimes p_{T,\mathrm{miss}}$",
        "mx2_x_theta_gg": r"$M_X^2\otimes\theta_{\gamma\gamma}$",
        "shared_ptmiss_theta": (
            r"shared: $(M_X^2\otimes p_T)+(M_X^2\otimes\theta_{\gamma\gamma})$"
        ),
    }

    for irow, region in enumerate(("FT", "FD_all")):
        for model, label in labels.items():
            rr = sorted(
                [
                    r for r in rows
                    if r["region"] == region
                    and r["candidate_model"] == model
                    and int(r["fit_success"]) == 1
                ],
                key=lambda r: float(r["energy_center_GeV"]),
            )
            if not rr:
                continue
            #endif

            x = [float(r["energy_center_GeV"]) for r in rr]
            axes[irow, 0].plot(
                x,
                [float(r["pi0_fraction"]) for r in rr],
                marker="o",
                linewidth=1.0,
                label=label,
            )
            axes[irow, 1].plot(
                x,
                [float(r["deviance_per_ndof"]) for r in rr],
                marker="o",
                linewidth=1.0,
                label=label,
            )
        #endfor

        axes[irow, 0].set_ylim(0.0, 1.02)
        axes[irow, 0].set_ylabel(r"$f_{\pi^0}$")
        axes[irow, 0].set_title(f"{region}: extracted composition")
        axes[irow, 1].set_ylabel("Poisson deviance / ndof")
        axes[irow, 1].set_title(f"{region}: absolute template goodness")
        for ax in axes[irow]:
            ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
            ax.grid(alpha=0.18)
        #endfor
    #endfor

    axes[0, 0].legend(fontsize=7.2, frameon=False)
    axes[0, 1].legend(fontsize=7.2, frameon=False)
    fig.suptitle(
        f"{period.label}: denominator-driver candidate selection "
        "(DIAGNOSTIC; does not feed Stage III)",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_candidate_model_summary.png",
        rect=(0, 0, 1, 0.95),
    )
    plt.close(fig)


def make_actual_nominal_stage2_driver_canvases(
    period: PeriodConfig,
    shared_rows: List[Dict[str, object]],
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    outdir: Path,
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
) -> None:
    """
    High-value production diagnostic: show BOTH ACTUAL nominal Stage-II drivers
    in EVERY nominal energy bin.

    The shared fit uses the same f_pi0 simultaneously in:
      1. M_X^2(ep) x M_probe^2
      2. M_X^2(ep) x pTmiss

    This canvas projects each 2D driver onto its second coordinate.  The curves
    are reconstructed from the same shared fit row and morph nuisance values
    that feed Stage III.  These are therefore not example/independent fits.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    edges = stage2_energy_edges(max_probe_energy)

    row_map = {
        (
            str(r.get("region", "")),
            round(float(r.get("energy_low_GeV", np.nan)), 10),
            round(float(r.get("energy_high_GeV", np.nan)), 10),
        ): r
        for r in shared_rows
        if int(r.get("fit_success", 0)) == 1
    }

    for region in ("FT", "FD_all"):
        nrows = len(edges) - 1
        fig, axes = plt.subplots(
            nrows,
            2,
            figsize=(13.5, max(8.0, 2.8 * nrows)),
            squeeze=False,
        )

        for ib in range(nrows):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])
            row = row_map.get(
                (region, round(elo, 10), round(ehi, 10))
            )

            for j, disc in enumerate(STAGE2_DRIVER_DISCRIMINATORS):
                ax = axes[ib, j]
                if row is None:
                    ax.text(
                        0.5, 0.5, "no successful shared fit",
                        transform=ax.transAxes, ha="center", va="center"
                    )
                    ax.set_axis_off()
                    continue
                #endif

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
                    disc,
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
                    display_bin_factor=1,
                )

                if not np.allclose(
                    model, pi0_c + dvcs_c,
                    rtol=1.0e-12, atol=1.0e-10
                ):
                    raise RuntimeError(
                        f"{period.label} {region} {elo:.2f}-{ehi:.2f}: "
                        "nominal driver display components do not sum to model."
                    )
                #endif

                ax.errorbar(
                    centers,
                    data,
                    yerr=np.sqrt(np.maximum(data, 1.0)),
                    fmt="o",
                    ms=2.3,
                    linewidth=0.7,
                    capsize=0.0,
                    color=SAMPLE_COLORS["data"],
                    label="data",
                    zorder=6,
                )
                ax.step(
                    centers, raw_dvcs_c, where="mid",
                    color=SAMPLE_COLORS["dvcs_mc"],
                    linewidth=0.65, linestyle="--", alpha=0.55,
                    label="BH/DVCS pre-morph",
                )
                ax.step(
                    centers, raw_pi0_c, where="mid",
                    color=SAMPLE_COLORS["pi0_mc"],
                    linewidth=0.65, linestyle="--", alpha=0.55,
                    label=r"$\pi^0$ pre-morph",
                )
                ax.step(
                    centers, dvcs_c, where="mid",
                    color=SAMPLE_COLORS["dvcs_mc"],
                    linewidth=1.25,
                    label="BH/DVCS fitted",
                )
                ax.step(
                    centers, pi0_c, where="mid",
                    color=SAMPLE_COLORS["pi0_mc"],
                    linewidth=1.25,
                    label=r"$\pi^0$ fitted",
                )
                ax.step(
                    centers, model, where="mid",
                    color=SAMPLE_COLORS["fit"],
                    linewidth=1.6,
                    label="total shared fit",
                )

                if disc == "mx2_ep_x_probe_m2":
                    ax.set_xlabel(
                        r"$(P_{\mathrm{probe}}^{\mathrm{pred}})^2$ "
                        r"(GeV$^2$)"
                    )
                    driver_label = (
                        r"$M_X^2(ep)\otimes M_{\mathrm{probe}}^2$"
                    )
                else:
                    ax.set_xlabel(r"$p_{T,\mathrm{miss}}$ (GeV)")
                    if region == "FT":
                        ax.set_xlim(0.0, 0.30)
                    else:
                        ax.set_xlim(0.0, ptmiss_max)
                    #endif
                    driver_label = (
                        r"$M_X^2(ep)\otimes p_{T,\mathrm{miss}}$"
                    )
                #endif

                ax.set_ylabel("entries / bin")
                ax.set_title(
                    f"{elo:.2f}-{ehi:.2f} GeV; {driver_label}\n"
                    + rf"shared nominal $f_{{\pi^0}}={float(row['pi0_fraction']):.3f}$, "
                    + rf"$D/\mathrm{{ndof}}={float(row['deviance_per_ndof']):.2f}$",
                    fontsize=8.6,
                )
                ax.grid(alpha=0.18)
            #endfor
        #endfor

        axes[0, 0].legend(fontsize=6.7, frameon=False, ncol=2)
        fig.suptitle(
            f"{period.label}: ACTUAL nominal Stage-II shared-fit projections — "
            f"{region}\n"
            r"Both columns are used simultaneously to determine the same "
            r"$f_{\pi^0}$ that feeds the efficiency",
            fontsize=13,
        )
        safe_finalize_figure(
            fig,
            outdir / (
                "canvas_nominal_stage2_drivers_"
                + region.lower()
                + ".png"
            ),
            rect=(0, 0, 1, 0.965),
        )
        plt.close(fig)
    #endfor


def make_theta_gg_alternative_canvas(
    period: PeriodConfig,
    fit_rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Compare theta_gg alternatives against the two nominal-driver independent
    fits.  This is diagnostic only.

    theta_gg_1d is especially important because it tests theta_gg WITHOUT
    relying on M_X^2(ep).  mx2_ep_x_theta_gg tests whether theta_gg helps when
    combined with M_X^2(ep).
    """
    choices = (
        (
            "mx2_ep_x_probe_m2",
            r"$M_X^2\otimes M_{\mathrm{probe}}^2$",
        ),
        (
            "mx2_ep_x_pTmiss",
            r"$M_X^2\otimes p_{T,\mathrm{miss}}$",
        ),
        (
            "mx2_ep_x_theta_gg",
            r"$M_X^2\otimes\theta_{\gamma\gamma}$",
        ),
        (
            "theta_gg_1d",
            r"$\theta_{\gamma\gamma}$ only",
        ),
    )

    fig, axes = plt.subplots(2, 2, figsize=(15.0, 9.2))

    for irow, region in enumerate(("FT", "FD_all")):
        for disc, label in choices:
            rr = sorted(
                [
                    r for r in fit_rows
                    if str(r.get("region", "")) == region
                    and str(r.get("discriminator", "")) == disc
                    and int(r.get("fit_success", 0)) == 1
                    and np.isfinite(float(r.get("pi0_fraction", np.nan)))
                ],
                key=lambda r: float(r["energy_center_GeV"]),
            )
            if not rr:
                continue
            #endif

            x = [float(r["energy_center_GeV"]) for r in rr]
            axes[irow, 0].plot(
                x,
                [float(r["pi0_fraction"]) for r in rr],
                marker="o",
                linewidth=1.0,
                label=label,
            )
            axes[irow, 1].plot(
                x,
                [float(r["deviance_per_ndof"]) for r in rr],
                marker="o",
                linewidth=1.0,
                label=label,
            )
        #endfor

        axes[irow, 0].set_ylim(0.0, 1.02)
        axes[irow, 0].set_ylabel(r"independent-fit $f_{\pi^0}$")
        axes[irow, 0].set_title(f"{region}: extracted composition")
        axes[irow, 1].set_ylabel("Poisson deviance / ndof")
        axes[irow, 1].set_title(f"{region}: template goodness")

        for ax in axes[irow]:
            ax.set_xlabel(
                r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
            )
            ax.grid(alpha=0.18)
        #endfor
    #endfor

    axes[0, 0].legend(fontsize=7.4, frameon=False)
    axes[0, 1].legend(fontsize=7.4, frameon=False)

    fig.suptitle(
        f"{period.label}: diagnostic alternatives to the nominal Stage-II drivers\n"
        r"$\theta_{\gamma\gamma}$-only explicitly tests separation without "
        r"$M_X^2(ep)$",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_theta_gg_alternative.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)


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
    DIAGNOSTIC ONLY: representative 1D pTmiss projections of three ACTUAL
    nominal shared Stage-II fit rows.

    This function does NOT perform separate/example fits.  The nominal fit is
    performed in every Stage-II energy bin with BOTH production 2D drivers:

        M_X^2(ep) x M_probe^2
        M_X^2(ep) x pTmiss

    The canvas merely projects the pTmiss driver for three representative
    energy bins (0.40-0.50, 0.80-1.00, 1.50-2.00 GeV) so the component shapes
    can be inspected visually.  It is available only in --plot-mode full.
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
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=(
                    "This figure includes Axes that are not compatible "
                    "with tight_layout.*"
                ),
                category=UserWarning,
            )
            fig.tight_layout(rect=rect)
        #endwith
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
            ax.set_yscale("linear")
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
        ax.set_yscale("linear")
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
        ax.set_yscale("linear")
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
    axes[0, 1].set_yscale("linear")
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
    axes[1, 1].set_yscale("linear")
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
    axes[1, 1].set_yscale("linear")
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
        write_rows_csv(
            final_reference_input_rows(stage3_rows),
            Path(stage3_outroot) / "final_reference_inputs.csv",
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

    write_analysis_flow_guide(
        Path(summary_outroot).parent,
        Path(stage2_outroot),
        Path(stage3_outroot),
    )

    make_summary_dashboard(
        shared_rows,
        disc_rows,
        closure_rows,
        stage3_rows,
        ft3_summary_rows,
        Path(summary_outroot),
    )


# -----------------------------------------------------------------------------
# Preflight
# -----------------------------------------------------------------------------

def preflight(
    periods: Sequence[PeriodConfig],
    include_stage2: bool,
    include_stage3: bool = False,
) -> None:
    """Validate required files and exact data-event identifiers before processing."""
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
            "The following required input files do not exist:\n  "
            + "\n  ".join(missing)
        )
    #endif

    if include_stage2:
        detector_missing: List[str] = []
        for p in periods:
            with uproot.open(p.epgamma_data) as root_file:
                found_name, tree = find_tree(root_file, None)
                if "detector2" not in set(tree.keys()):
                    detector_missing.append(
                        f"{p.label}: epgamma data tree '{found_name}' missing detector2"
                    )
                #endif
            #endwith
        #endfor
        if detector_missing:
            raise KeyError(
                "Data-derived photon angular acceptance requires detector2:\n  "
                + "\n  ".join(detector_missing)
            )
        #endif
    #endif

    if include_stage3:
        log(
            "Preflight efficiency stage: verifying exact data-event identifiers "
            "(runnum, evnum) before expensive processing."
        )
        branch_missing: List[str] = []

        for p in periods:
            for label, path in (
                ("epgamma data", p.epgamma_data),
                ("eppi0 data", p.eppi0_data),
            ):
                with uproot.open(path) as root_file:
                    found_name, tree = find_tree(root_file, None)
                    available = set(tree.keys())
                    absent = [
                        b for b in ("runnum", "evnum")
                        if b not in available
                    ]
                    if absent:
                        branch_missing.append(
                            f"{p.label}: {label}: tree '{found_name}' missing "
                            + ", ".join(absent)
                        )
                    else:
                        log(
                            f"Preflight efficiency OK: {p.label}: {label}: "
                            "runnum+evnum available."
                        )
                    #endif
                #endwith
            #endfor
        #endfor

        if branch_missing:
            raise KeyError(
                "Exact data matching cannot run because:\n  "
                + "\n  ".join(branch_missing)
            )
        #endif
    #endif


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
        default=0.4,
        help=(
            "Minimum reconstructed epgamma tag energy in GeV. "
            "Default: 0.4. The nSidis production workflow intentionally "
            "keeps reconstructed tag photons down to the skim threshold."
        ),
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
            "component by --parent-component-tol. Default: 0.4."
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

    parser.add_argument(
        "--skip-grand-diagnostics",
        action="store_true",
        help=(
            "Skip the broad Stage-1 grand diagnostic canvases in a full "
            "production run. By default they are produced because they are "
            "part of the visual validation record."
        ),
    )
    parser.add_argument(
        "--grand-diagnostics-max-entries",
        type=int,
        default=4000000,
        help=(
            "Maximum entries per epgamma sample used only for the broad grand "
            "diagnostic canvases during a full production run. Default: "
            "4,000,000. Use 0 for all entries. This does not truncate the "
            "actual Stage-1/2/3 production analysis."
        ),
    )

    parser.add_argument(
        "--stage1-grand-diagnostics",
        action="store_true",
        help=(
            "When combined with --stage1-only, also run the broad all-variable "
            "grand-diagnostic suite. This is opt-in because it reads dozens "
            "of extra branches and is much heavier in memory/I/O."
        ),
    )

    parser.add_argument(
        "--stage1-only",
        action="store_true",
        help=(
            "Run only the nSidis Stage-1 shape diagnostics and exit. "
            "This reads nSidis data, aaogen epgamma, and dvcsgen epgamma, "
            "makes the standard FT/FD shape canvases plus the grand "
            "all-variable diagnostic suite, and skips eppi0 association, "
            "template fitting, and efficiency extraction."
        ),
    )

    parser.add_argument(
        "--nsidis-study",
        action="store_true",
        help=(
            "Run the isolated loose-nSidis overlap/exclusivity study and exit "
            "before the normal Stage-II/III chain."
        ),
    )
    parser.add_argument(
        "--nsidis-epgamma",
        default=None,
        help=(
            "Override nSidis epgamma ROOT file. Requires exactly one --period."
        ),
    )
    parser.add_argument(
        "--nsidis-eppi0",
        default=None,
        help=(
            "Override nSidis eppi0 ROOT file. Requires exactly one --period."
        ),
    )
    parser.add_argument(
        "--nsidis-probe-m2-scan",
        default="0.02,0.04,0.06,0.08,0.10,0.15,0.20,0.30,0.60",
        help=(
            "Comma-separated |M_X^2(epgamma_tag)| support cuts scanned in the "
            "nSidis overlap study. pTmiss is deliberately left uncut."
        ),
    )
    parser.add_argument(
        "--nsidis-pilot-fit",
        action="store_true",
        help=(
            "Run the nominal nSidis M_X^2(ep) x pTmiss composition fit and "
            "association-first photon-efficiency extraction."
        ),
    )
    parser.add_argument(
        "--nsidis-driver-study",
        action="store_true",
        help=(
            "Run the dedicated denominator-variable study on the nSidis sample. "
            "Compares pTmiss only, theta_gamma_gamma only, |Delta_phi| "
            "only, theta_gamma_gamma x pTmiss, |Delta_phi| x pTmiss, "
            "M_X^2(ep) x pTmiss, and M_X^2(ep) x theta_gamma_gamma on one "
            "common event population. "
            "When enabled with --nsidis-study, the worker exits after this "
            "study and skips the reconstructed-probe efficiency calculation."
        ),
    )
    parser.add_argument(
        "--nsidis-pilot-energy-max",
        type=float,
        default=9.5,
        help=(
            "Maximum predicted probe energy for the optional nSidis pilot fit. "
            "Default: 9.5 GeV."
        ),
    )
    parser.add_argument(
        "--nsidis-eppi0-core-nsigma-scan",
        default="2.5,3.0,3.5,4.0,5.0,6.0",
        help=(
            "Comma-separated Gaussian-core n-sigma values for the nSidis "
            "eppi0 exclusivity scan. Default: 2,2.5,3,3.5,4."
        ),
    )
    parser.add_argument(
        "--nsidis-eppi0-momentum-quantile-scan",
        default="0.95,0.975,0.99,0.995,0.999",
        help=(
            "Comma-separated old-wagon pmiss/pTmiss core quantiles scanned "
            "after the Mx2 core requirement. Default: 0.95,0.975,0.99,0.995,0.999."
        ),
    )
    parser.add_argument(
        "--nsidis-eppi0-target-wagon-retention",
        type=float,
        default=0.95,
        help=(
            "Minimum old-wagon e_p>2 retention used to choose the diagnostic "
            "recommended eppi0 exclusive-core point. Default: 0.95."
        ),
    )
    parser.add_argument(
        "--nsidis-pilot-probe-m2-values",
        default="0.06,0.08,0.10",
        help=(
            "Comma-separated mass-shell support values tested in the nSidis "
            "production-model pilot. Default: 0.06,0.08,0.10 GeV^2."
        ),
    )
    parser.add_argument(
        "--nsidis-pilot-probe-m2-max",
        type=float,
        default=0.08,
        help=(
            "Central mass-shell support used for the detailed nSidis pilot fit "
            "projection canvases. Default: 0.08 GeV^2."
        ),
    )
    parser.add_argument(
        "--diagnostics",
        choices=("selection", "full"),
        default="selection",
        help=(
            "Diagnostic workload. 'selection' (default) runs only the fits and "
            "plots needed for the current discriminator-selection study and skips "
            "expensive closure/mixed/coarse-FT/profile/control diagnostics. "
            "'full' restores the complete development diagnostic suite."
        ),
    )

    # Stage-II denominator-composition controls.
    
    parser.add_argument(
        "--stage2-output-dir",
        default="output/photon_efficiency/stage2",
        help="Stage-II output directory. Default: output/photon_efficiency/stage2.",
    )
    parser.add_argument(
        "--pt-category-split",
        type=float,
        default=0.10,
        help=(
            "SIDIS pT split (GeV) for the two-category composition cross-check. "
            "pT below the split is DVCS-rich; pT at/above the split is pi0-rich. "
            "Both categories are retained and fit simultaneously. Default: 0.10."
        ),
    )

    parser.add_argument(
        "--den-fit-mm2-min",
        type=float,
        default=-0.231,
        help=(
            "Lower M_X^2(ep) selection/fit edge (GeV^2). "
            "Default: -0.231 (Valerii photon-efficiency cut)."
        ),
    )
    parser.add_argument(
        "--den-fit-mm2-max",
        type=float,
        default=0.309,
        help=(
            "Fallback upper M_X^2(epX) selection edge (GeV^2) used only "
            "when a CLASDIS-derived FD pi0/eta valley cannot be determined. "
            "In normal production the period-specific upper edge is derived "
            "from the FD CLASDIS distribution and then imposed identically "
            "in FT and FD. Default fallback: 0.309."
        ),
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
        default=120,
        help=(
            "Number of predicted-probe M^2 bins in each Stage-II template fit. "
            "Default: 120."
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
        default=100,
        help=(
            "pTmiss bins for the optional 2D discriminator fit. Default: 100."
        ),
    )
    parser.add_argument(
        "--disc-theta-max",
        type=float,
        default=4.0,
        help=(
            "Upper theta_gamma_gamma edge (deg) used by the optional "
            "M_X^2(ep) x theta_gamma_gamma discriminator. Default: 4."
        ),
    )
    parser.add_argument(
        "--disc-theta-bins",
        type=int,
        default=100,
        help=(
            "theta_gamma_gamma bins for the optional 2D discriminator fit. "
            "Default: 100."
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




# -----------------------------------------------------------------------------
# Loose nSidis transition study
# -----------------------------------------------------------------------------

def resolve_nsidis_paths(
    period: PeriodConfig,
    args_dict: Dict[str, object],
) -> Tuple[str, str]:
    """Resolve period defaults, optionally overridden from the command line."""
    epg_override = args_dict.get("nsidis_epgamma")
    epi_override = args_dict.get("nsidis_eppi0")

    epg_path = (
        str(epg_override)
        if epg_override
        else period.nsidis_epgamma_data
    )
    epi_path = (
        str(epi_override)
        if epi_override
        else period.nsidis_eppi0_data
    )

    if not epg_path or not epi_path:
        raise FileNotFoundError(
            f"{period.label}: nSidis files are not configured. "
            "Supply --nsidis-epgamma and --nsidis-eppi0."
        )
    #endif
    return str(epg_path), str(epi_path)


def electron_momentum_from_p3(p3: np.ndarray) -> np.ndarray:
    """Magnitude of the reconstructed electron three-momentum."""
    p3 = np.asarray(p3, dtype=float)
    return np.sqrt(np.einsum("ij,ij->i", p3, p3))


def _structured_event_keys(
    raw: Dict[str, np.ndarray],
    mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Collision-free structured runnum/evnum event keys."""
    if "runnum" not in raw or "evnum" not in raw:
        raise KeyError("runnum and evnum are required for nSidis overlap checks.")
    #endif

    run = np.asarray(raw["runnum"], dtype=np.int64)
    ev = np.asarray(raw["evnum"], dtype=np.int64)
    if run.shape != ev.shape:
        raise ValueError("runnum and evnum arrays have different shapes.")
    #endif

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != run.shape:
            raise ValueError("event-key mask has incompatible shape.")
        #endif
        run = run[mask]
        ev = ev[mask]
    #endif

    keys = np.empty(
        run.size,
        dtype=np.dtype([("runnum", "<i8"), ("evnum", "<i8")]),
    )
    keys["runnum"] = run
    keys["evnum"] = ev
    return keys


def event_overlap_masks(
    reference_raw: Dict[str, np.ndarray],
    candidate_raw: Dict[str, np.ndarray],
    reference_mask: Optional[np.ndarray] = None,
    candidate_mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, Dict[str, object]]:
    """
    Exact event-level overlap.

    The returned mask has candidate-entry length and marks entries whose
    (runnum,evnum) occurs in the selected reference-event population.

    For the nSidis transition, reference_mask is used to impose e_p > 2 GeV on
    the old wagon before computing recall, matching the nSidis parent skim.
    """
    ref_keys = _structured_event_keys(reference_raw, reference_mask)
    cand_all = _structured_event_keys(candidate_raw, None)
    cand_keys = _structured_event_keys(candidate_raw, candidate_mask)

    ref_unique = np.unique(ref_keys)
    cand_unique = np.unique(cand_keys)
    common = np.intersect1d(ref_unique, cand_unique, assume_unique=True)

    candidate_overlap = np.isin(cand_all, ref_unique)
    if candidate_mask is not None:
        candidate_overlap &= np.asarray(candidate_mask, dtype=bool)
    #endif

    summary = {
        "reference_entries_after_parent_cut": int(ref_keys.size),
        "candidate_entries_after_parent_cut": int(cand_keys.size),
        "reference_unique_events_after_parent_cut": int(ref_unique.size),
        "candidate_unique_events_after_parent_cut": int(cand_unique.size),
        "common_unique_events": int(common.size),
        "reference_unique_event_recall": (
            float(common.size) / float(ref_unique.size)
            if ref_unique.size else float("nan")
        ),
        "candidate_unique_event_overlap_fraction": (
            float(common.size) / float(cand_unique.size)
            if cand_unique.size else float("nan")
        ),
        "candidate_entry_overlap_fraction": (
            float(np.count_nonzero(candidate_overlap))
            / float(np.count_nonzero(candidate_mask))
            if candidate_mask is not None
            and np.count_nonzero(candidate_mask) > 0
            else (
                float(np.count_nonzero(candidate_overlap))
                / float(candidate_overlap.size)
                if candidate_overlap.size else float("nan")
            )
        ),
    }
    return candidate_overlap, summary


def build_eppi0_exclusivity_features(
    period: PeriodConfig,
    sample: EPPi0Sample,
) -> Dict[str, np.ndarray]:
    """
    Recompute e p -> e p pi0 exclusivity from measured four-vectors.

    This is intentionally independent of the upstream skim cuts and therefore
    provides a common old-wagon / loose-nSidis comparison.
    """
    e3 = sample.electron_p3
    p3 = sample.proton_p3
    pi3 = sample.pi0_p3

    e_E = np.sqrt(np.einsum("ij,ij->i", e3, e3) + M_E * M_E)
    p_E = np.sqrt(np.einsum("ij,ij->i", p3, p3) + M_P * M_P)
    pi_E = np.sqrt(
        np.einsum("ij,ij->i", pi3, pi3)
        + np.asarray(sample.pi0_mass, dtype=float) ** 2
    )

    beam_p = math.sqrt(max(0.0, period.beam_energy**2 - M_E**2))
    miss_E = period.beam_energy + M_P - e_E - p_E - pi_E

    px = -e3[:, 0] - p3[:, 0] - pi3[:, 0]
    py = -e3[:, 1] - p3[:, 1] - pi3[:, 1]
    pz = beam_p - e3[:, 2] - p3[:, 2] - pi3[:, 2]

    p2 = px * px + py * py + pz * pz
    pmiss = np.sqrt(p2)
    ptmiss = np.sqrt(px * px + py * py)
    mx2 = miss_E * miss_E - p2

    out = {
        "electron_p": np.asarray(
            electron_momentum_from_p3(e3), dtype=np.float32
        ),
        "Emiss": np.asarray(miss_E, dtype=np.float32),
        "pmiss": np.asarray(pmiss, dtype=np.float32),
        "pTmiss": np.asarray(ptmiss, dtype=np.float32),
        "Mx2": np.asarray(mx2, dtype=np.float32),
        "Mgg": np.asarray(sample.pi0_mass, dtype=np.float32),
    }
    if "runnum" in sample.raw:
        out["runnum"] = np.asarray(sample.raw["runnum"])
    #endif
    if "evnum" in sample.raw:
        out["evnum"] = np.asarray(sample.raw["evnum"])
    #endif
    return out


def _normalized_hist(
    values: np.ndarray,
    bins: np.ndarray,
    mask: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if mask is not None:
        good &= np.asarray(mask, dtype=bool)
    #endif
    hist, edges = np.histogram(values[good], bins=bins)
    hist = hist.astype(float)
    total = float(np.sum(hist))
    if total > 0.0:
        hist /= total
    #endif
    return hist, edges


def robust_gaussian_core(
    values: np.ndarray,
    *,
    broad_abs_max: float,
    iterations: int = 6,
    clip_sigma: float = 2.5,
) -> Tuple[float, float, int]:
    """
    Estimate the central near-zero exclusive peak while rejecting long tails.

    The first scale is MAD-based, then the mean/RMS are iterated inside a
    symmetric clip. This is used only to characterize the old-wagon exclusive
    core, not to fit the inclusive nSidis tails.
    """
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v) & (np.abs(v) < float(broad_abs_max))]
    if v.size < 100:
        return float("nan"), float("nan"), int(v.size)
    #endif

    center = float(np.median(v))
    mad = float(np.median(np.abs(v - center)))
    sigma = 1.4826 * mad
    if not np.isfinite(sigma) or sigma <= 0.0:
        sigma = float(np.std(v))
    #endif

    for _ in range(max(1, int(iterations))):
        keep = np.abs(v - center) < float(clip_sigma) * max(sigma, 1.0e-8)
        vv = v[keep]
        if vv.size < 100:
            break
        #endif
        new_center = float(np.mean(vv))
        new_sigma = float(np.std(vv, ddof=1))
        if not np.isfinite(new_sigma) or new_sigma <= 0.0:
            break
        #endif
        center, sigma = new_center, new_sigma
    #endfor
    return center, sigma, int(vv.size if "vv" in locals() else v.size)


def derive_eppi0_core_model(
    wagon_f: Dict[str, np.ndarray],
    parent_mask: np.ndarray,
) -> Dict[str, float]:
    """
    Characterize the sharp M_X^2(ep pi0) exclusive peak in the old wagon.

    E_miss is intentionally NOT Gaussianized and is not used as a cut.  Its
    broad/asymmetric shape is retained only as a diagnostic.  This avoids the
    previous failure in which an E_miss Gaussian-core requirement limited the
    maximum old-wagon retention to ~71%.
    """
    parent_mask = np.asarray(parent_mask, dtype=bool)

    mu_m, sig_m, n_m = robust_gaussian_core(
        wagon_f["Mx2"][parent_mask],
        broad_abs_max=0.20,
        iterations=8,
        clip_sigma=2.5,
    )
    if not np.isfinite(sig_m) or sig_m <= 0.0:
        raise RuntimeError("Could not determine old-wagon Mx2 exclusive core.")
    #endif

    # Quantiles of E_miss are recorded only to document its asymmetric/non-core
    # behavior; they never enter the selection mask.
    emiss = np.asarray(wagon_f["Emiss"][parent_mask], dtype=float)
    emiss = emiss[np.isfinite(emiss)]
    if emiss.size:
        eq = np.quantile(emiss, [0.01, 0.16, 0.50, 0.84, 0.99])
    else:
        eq = np.full(5, np.nan)
    #endif

    return {
        "Mx2_center_GeV2": float(mu_m),
        "Mx2_sigma_GeV2": float(sig_m),
        "Mx2_core_fit_count": int(n_m),
        "Emiss_used_in_selection": 0,
        "Emiss_q01_GeV": float(eq[0]),
        "Emiss_q16_GeV": float(eq[1]),
        "Emiss_q50_GeV": float(eq[2]),
        "Emiss_q84_GeV": float(eq[3]),
        "Emiss_q99_GeV": float(eq[4]),
    }


def eppi0_core_mask(
    features: Dict[str, np.ndarray],
    model: Dict[str, float],
    nsigma: float,
    pmiss_max: float,
    ptmiss_max: float,
    parent_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Exclusive-candidate mask.

    Selection variables:
      * M_X^2(ep pi0): symmetric robust-core window;
      * |p_miss(ep pi0)|: upper threshold;
      * pTmiss(ep pi0): upper threshold.

    E_miss is diagnostic-only.
    """
    mask = (
        np.isfinite(features["Mx2"])
        & np.isfinite(features["pmiss"])
        & np.isfinite(features["pTmiss"])
        & (
            np.abs(features["Mx2"] - model["Mx2_center_GeV2"])
            <= float(nsigma) * model["Mx2_sigma_GeV2"]
        )
        & (features["pmiss"] <= float(pmiss_max))
        & (features["pTmiss"] <= float(ptmiss_max))
    )
    if parent_mask is not None:
        mask &= np.asarray(parent_mask, dtype=bool)
    #endif
    return mask


def build_eppi0_core_scan(
    period: PeriodConfig,
    wagon_f: Dict[str, np.ndarray],
    nsidis_f: Dict[str, np.ndarray],
    nsidis_event_overlap: np.ndarray,
    wagon_parent_mask: np.ndarray,
    nsidis_parent_mask: np.ndarray,
    model: Dict[str, float],
    nsigma_values: Sequence[float],
    momentum_quantiles: Sequence[float],
    target_wagon_retention: float,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    """
    Scan an Mx2-core + missing-momentum exclusivity selection.

    For each Mx2 n-sigma window, the pmiss and pTmiss thresholds are taken from
    old-wagon quantiles AFTER the Mx2 requirement.  The chosen working point is
    the most rejecting point that ACTUALLY reaches target_wagon_retention.

    If no point reaches the target, the returned status explicitly says so; it
    is never labeled "recommended".
    """
    rows: List[Dict[str, object]] = []
    wagon_parent_mask = np.asarray(wagon_parent_mask, dtype=bool)
    nsidis_parent_mask = np.asarray(nsidis_parent_mask, dtype=bool)

    n_wagon_parent = int(np.count_nonzero(wagon_parent_mask))
    n_ns_nonwagon_parent = int(
        np.count_nonzero(nsidis_parent_mask & ~nsidis_event_overlap)
    )
    n_ns_overlap_parent = int(
        np.count_nonzero(nsidis_parent_mask & nsidis_event_overlap)
    )

    for nsigma in nsigma_values:
        mx2_core = (
            wagon_parent_mask
            & (
                np.abs(wagon_f["Mx2"] - model["Mx2_center_GeV2"])
                <= float(nsigma) * model["Mx2_sigma_GeV2"]
            )
        )

        pvals = np.asarray(wagon_f["pmiss"][mx2_core], dtype=float)
        ptvals = np.asarray(wagon_f["pTmiss"][mx2_core], dtype=float)
        pvals = pvals[np.isfinite(pvals)]
        ptvals = ptvals[np.isfinite(ptvals)]
        if pvals.size < 100 or ptvals.size < 100:
            continue
        #endif

        for q in momentum_quantiles:
            q = float(q)
            pmiss_max = float(np.quantile(pvals, q))
            ptmiss_max = float(np.quantile(ptvals, q))

            wm = eppi0_core_mask(
                wagon_f, model, nsigma, pmiss_max, ptmiss_max,
                wagon_parent_mask,
            )
            nm = eppi0_core_mask(
                nsidis_f, model, nsigma, pmiss_max, ptmiss_max,
                nsidis_parent_mask,
            )
            no = nm & nsidis_event_overlap
            nx = nm & ~nsidis_event_overlap

            wc = int(np.count_nonzero(wm))
            nc = int(np.count_nonzero(nm))
            oc = int(np.count_nonzero(no))
            xc = int(np.count_nonzero(nx))

            xF2_lo, xF2_hi, _ = NOMINAL_THREE_RANGES["xF2"]
            xF2_values_data = _nominal_three_values(
                data_f,
                "xF2",
            )
            xF2_selected = xF2_values_data[masks["data"]]
            xF2_underflow_count = int(
                np.count_nonzero(
                    np.isfinite(xF2_selected)
                    & (xF2_selected < xF2_lo)
                )
            )
            xF2_overflow_count = int(
                np.count_nonzero(
                    np.isfinite(xF2_selected)
                    & (xF2_selected >= xF2_hi)
                )
            )

            rows.append({
                "period": period.key,
                "label": period.label,
                "nsigma_Mx2": float(nsigma),
                "momentum_quantile": q,
                "Mx2_low_GeV2": (
                    model["Mx2_center_GeV2"]
                    - float(nsigma) * model["Mx2_sigma_GeV2"]
                ),
                "Mx2_high_GeV2": (
                    model["Mx2_center_GeV2"]
                    + float(nsigma) * model["Mx2_sigma_GeV2"]
                ),
                "pmiss_max_GeV": pmiss_max,
                "pTmiss_max_GeV": ptmiss_max,
                "Emiss_cut_applied": 0,
                "wagon_epgt2_retention": (
                    wc / n_wagon_parent if n_wagon_parent else float("nan")
                ),
                "nsidis_overlap_event_entry_retention": (
                    oc / n_ns_overlap_parent
                    if n_ns_overlap_parent else float("nan")
                ),
                "nsidis_nonwagon_event_entry_retention": (
                    xc / n_ns_nonwagon_parent
                    if n_ns_nonwagon_parent else float("nan")
                ),
                "wagon_surviving_count": wc,
                "nsidis_surviving_count": nc,
                "nsidis_overlap_event_surviving_count": oc,
                "nsidis_nonwagon_event_surviving_count": xc,
                "nsidis_surviving_overlap_event_fraction": (
                    oc / nc if nc else float("nan")
                ),
            })
        #endfor
    #endfor

    if not rows:
        raise RuntimeError("eppi0 core scan produced no valid working points.")
    #endif

    eligible = [
        row for row in rows
        if np.isfinite(row["wagon_epgt2_retention"])
        and float(row["wagon_epgt2_retention"])
        >= float(target_wagon_retention)
    ]

    if eligible:
        chosen = min(
            eligible,
            key=lambda row: (
                float(row["nsidis_nonwagon_event_entry_retention"]),
                -float(row["wagon_epgt2_retention"]),
            ),
        )
        status = "recommended_diagnostic_working_point"
        meets_target = 1
    else:
        chosen = max(
            rows,
            key=lambda row: float(row["wagon_epgt2_retention"]),
        )
        status = "NO_RECOMMENDATION_target_retention_not_reached"
        meets_target = 0
    #endif

    chosen = dict(chosen)
    chosen["target_wagon_retention"] = float(target_wagon_retention)
    chosen["meets_target_retention"] = int(meets_target)
    chosen["selection_status"] = status
    return rows, chosen


def make_nsidis_epgamma_overlap_canvas(
    period: PeriodConfig,
    wagon_f: Dict[str, np.ndarray],
    nsidis_f: Dict[str, np.ndarray],
    nsidis_event_overlap: np.ndarray,
    wagon_parent_mask: np.ndarray,
    nsidis_parent_mask: np.ndarray,
    outdir: Path,
) -> None:
    """Old-vs-new epgamma comparison below 2 GeV with matched e_p>2 support."""
    wagon_mask = (
        np.asarray(wagon_parent_mask, dtype=bool)
        & wagon_f["valid_tag"]
        & (wagon_f["pred_probe_energy"] >= 0.40)
        & (wagon_f["pred_probe_energy"] < 2.0)
    )
    ns_mask = (
        np.asarray(nsidis_parent_mask, dtype=bool)
        & nsidis_f["valid_tag"]
        & (nsidis_f["pred_probe_energy"] >= 0.40)
        & (nsidis_f["pred_probe_energy"] < 2.0)
    )

    specs = (
        (
            "ep_missing_mass2",
            np.linspace(-0.30, 0.40, 141),
            r"$M_X^2(ep)$ (GeV$^2$)",
        ),
        (
            "pred_probe_mass2",
            np.linspace(-0.40, 0.40, 161),
            r"$M_X^2(ep\gamma_{\rm tag})$ (GeV$^2$)",
        ),
        (
            "stored_pTmiss",
            np.linspace(0.0, 1.0, 121),
            r"$p_{T,\rm miss}(ep\gamma)$ (GeV)",
        ),
        (
            "pred_probe_energy",
            np.linspace(0.40, 2.00, 81),
            r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)",
        ),
    )

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    for ax, (key, bins, xlabel) in zip(axes.ravel(), specs):
        h_old, edges = _normalized_hist(
            wagon_f[key], bins, wagon_mask
        )
        h_overlap, _ = _normalized_hist(
            nsidis_f[key],
            bins,
            ns_mask & nsidis_event_overlap,
        )
        h_extra, _ = _normalized_hist(
            nsidis_f[key],
            bins,
            ns_mask & ~nsidis_event_overlap,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        ax.step(
            centers, h_old, where="mid", linewidth=1.3,
            label=r"original wagon, $e_p>2$ GeV",
        )
        ax.step(
            centers, h_overlap, where="mid", linewidth=1.1,
            label="nSidis: event also in wagon",
        )
        ax.step(
            centers, h_extra, where="mid", linewidth=1.0,
            label="nSidis: event absent from wagon",
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel("unit-normalized entries")
        ax.grid(alpha=0.18)
    #endfor

    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"{period.label}: old wagon versus loose nSidis epgammaX, "
        r"matched $e_p>2$ GeV support and "
        r"$0.4<E_{\rm probe}^{pred}<2$ GeV",
        fontsize=13.0,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_epgamma_overlap_below2.png",
        rect=(0, 0, 1, 0.95),
    )
    plt.close(fig)


def make_nsidis_overlap_vs_energy_canvas(
    period: PeriodConfig,
    nsidis_f: Dict[str, np.ndarray],
    nsidis_event_overlap: np.ndarray,
    nsidis_parent_mask: np.ndarray,
    outdir: Path,
    tag_max: float,
) -> None:
    valid = (
        np.asarray(nsidis_parent_mask, dtype=bool)
        & nsidis_f["valid_tag"]
        & (nsidis_f["pred_probe_energy"] >= 0.40)
        & (nsidis_f["pred_probe_energy"] < tag_max)
    )

    edges = np.asarray(
        [0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 2.00,
         2.50, 3.00, 4.00, 5.00, 6.50, 8.00, tag_max],
        dtype=float,
    )
    edges = np.unique(
        edges[(edges >= 0.40) & (edges <= tag_max)]
    )
    if edges.size < 2:
        return
    #endif
    if edges[-1] < tag_max:
        edges = np.append(edges, tag_max)
    #endif

    total, _ = np.histogram(
        nsidis_f["pred_probe_energy"][valid],
        bins=edges,
    )
    overlap, _ = np.histogram(
        nsidis_f["pred_probe_energy"][
            valid & nsidis_event_overlap
        ],
        bins=edges,
    )
    frac = np.divide(
        overlap.astype(float),
        total.astype(float),
        out=np.full(total.shape, np.nan, dtype=float),
        where=total > 0,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])

    fig, ax = plt.subplots(figsize=(9.5, 5.6))
    ax.plot(centers, frac, marker="o", linewidth=1.2)
    ax.axvline(2.0, linestyle="--", linewidth=1.0)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
    ax.set_ylabel("nSidis entries whose event occurs in old wagon")
    ax.set_title(
        f"{period.label}: old-wagon event overlap versus predicted probe energy"
    )
    ax.grid(alpha=0.18)
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_wagon_overlap_vs_energy.png",
    )
    plt.close(fig)


def make_nsidis_eppi0_core_canvas(
    period: PeriodConfig,
    wagon_f: Dict[str, np.ndarray],
    nsidis_f: Dict[str, np.ndarray],
    nsidis_event_overlap: np.ndarray,
    wagon_parent_mask: np.ndarray,
    nsidis_parent_mask: np.ndarray,
    model: Dict[str, float],
    recommended: Dict[str, object],
    outdir: Path,
) -> None:
    """
    Show the old-wagon exclusive population against the loose nSidis sample.

    E_miss is displayed precisely so its broad/asymmetric behavior is visible,
    but NO E_miss cut is drawn or applied.
    """
    specs = (
        (
            "Mx2",
            np.linspace(-0.20, 0.20, 181),
            r"$M_X^2(ep\pi^0)$ (GeV$^2$)",
            (
                float(recommended["Mx2_low_GeV2"]),
                float(recommended["Mx2_high_GeV2"]),
            ),
            "selection variable",
        ),
        (
            "pTmiss",
            np.linspace(0.0, 0.80, 161),
            r"$p_{T,\rm miss}(ep\pi^0)$ (GeV)",
            (float(recommended["pTmiss_max_GeV"]),),
            "selection variable",
        ),
        (
            "pmiss",
            np.linspace(0.0, 1.5, 181),
            r"$|\vec p_{\rm miss}(ep\pi^0)|$ (GeV)",
            (float(recommended["pmiss_max_GeV"]),),
            "selection variable",
        ),
        (
            "Emiss",
            np.linspace(-1.5, 2.5, 201),
            r"$E_{\rm miss}(ep\pi^0)$ (GeV)",
            (),
            "DIAGNOSTIC ONLY — no cut applied",
        ),
    )

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    for ax, (key, bins, xlabel, cut_lines, subtitle) in zip(
        axes.ravel(), specs
    ):
        h_old, edges = _normalized_hist(
            wagon_f[key], bins, wagon_parent_mask
        )
        h_overlap, _ = _normalized_hist(
            nsidis_f[key],
            bins,
            nsidis_parent_mask & nsidis_event_overlap,
        )
        h_extra, _ = _normalized_hist(
            nsidis_f[key],
            bins,
            nsidis_parent_mask & ~nsidis_event_overlap,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        ax.step(
            centers, h_old, where="mid", linewidth=1.3,
            label=r"old wagon, $e_p>2$ GeV",
        )
        ax.step(
            centers, h_overlap, where="mid", linewidth=1.1,
            label="nSidis: event also in wagon",
        )
        ax.step(
            centers, h_extra, where="mid", linewidth=1.0,
            label="nSidis: event absent from wagon",
        )
        for value in cut_lines:
            ax.axvline(value, linestyle="--", linewidth=0.9)
        #endfor
        ax.set_xlabel(xlabel)
        ax.set_ylabel("unit-normalized entries")
        ax.set_title(subtitle, fontsize=10)
        ax.grid(alpha=0.18)
    #endfor

    axes[0, 0].legend(fontsize=8, frameon=False)
    meets = int(recommended.get("meets_target_retention", 0)) == 1
    status_text = (
        "recommended diagnostic point"
        if meets else
        "best scanned point — target retention NOT reached"
    )
    fig.suptitle(
        f"{period.label}: eppi0X exclusivity study — {status_text}\n"
        f"old-wagon retention = "
        f"{100.0*float(recommended['wagon_epgt2_retention']):.1f}%; "
        r"$E_{\rm miss}$ is shown but is not used",
        fontsize=12.5,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_eppi0_exclusive_core.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)



def make_nsidis_eppi0_core_scan_canvas(
    period: PeriodConfig,
    rows: List[Dict[str, object]],
    recommended: Dict[str, object],
    outdir: Path,
) -> None:
    if not rows:
        return
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))

    quantiles = sorted({
        float(row["momentum_quantile"]) for row in rows
    })
    for q in quantiles:
        rr = sorted(
            [
                row for row in rows
                if abs(float(row["momentum_quantile"]) - q) < 1.0e-12
            ],
            key=lambda row: float(row["nsigma_Mx2"]),
        )
        axes[0].plot(
            [float(row["wagon_epgt2_retention"]) for row in rr],
            [
                float(row["nsidis_nonwagon_event_entry_retention"])
                for row in rr
            ],
            marker="o",
            linewidth=1.0,
            label=f"momentum q={q:.3f}",
        )
    #endfor

    axes[0].scatter(
        [float(recommended["wagon_epgt2_retention"])],
        [float(recommended["nsidis_nonwagon_event_entry_retention"])],
        marker="*",
        s=110,
        zorder=8,
        label="recommended",
    )
    axes[0].set_xlabel(r"old-wagon $e_p>2$ retention")
    axes[0].set_ylabel("retention of nSidis entries from non-wagon events")
    axes[0].set_title("signal retention versus loose-background retention")
    axes[0].grid(alpha=0.18)
    axes[0].legend(fontsize=7.5, frameon=False)

    rr = sorted(rows, key=lambda row: (
        float(row["nsigma_Mx2"]), float(row["momentum_quantile"])
    ))
    x = np.arange(len(rr))
    axes[1].plot(
        x,
        [float(row["wagon_epgt2_retention"]) for row in rr],
        marker="o",
        linewidth=0.9,
        label="old-wagon retention",
    )
    axes[1].plot(
        x,
        [
            float(row["nsidis_nonwagon_event_entry_retention"])
            for row in rr
        ],
        marker="o",
        linewidth=0.9,
        label="non-wagon nSidis retention",
    )
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_xlabel("core-scan working point index")
    axes[1].set_ylabel("retention")
    axes[1].set_title("all scanned working points")
    axes[1].grid(alpha=0.18)
    axes[1].legend(fontsize=8, frameon=False)

    fig.suptitle(
        f"{period.label}: eppi0X Mx2-core + missing-momentum cut scan",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_eppi0_core_scan.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)


def build_nsidis_probe_mass2_scan(
    period: PeriodConfig,
    wagon_f: Dict[str, np.ndarray],
    nsidis_f: Dict[str, np.ndarray],
    nsidis_event_overlap: np.ndarray,
    wagon_parent_mask: np.ndarray,
    nsidis_parent_mask: np.ndarray,
    mm2_min: float,
    mm2_max: float,
    scan_values: Sequence[float],
) -> List[Dict[str, object]]:
    """
    Scan |M_X^2(epgamma_tag)| while deliberately leaving pTmiss uncut.
    """
    rows: List[Dict[str, object]] = []

    for range_name, elo, ehi in (
        ("overlap_0p4_2", 0.40, 2.00),
        ("extension_2_9p5", 2.00, 9.50),
    ):
        wagon_base = (
            np.asarray(wagon_parent_mask, dtype=bool)
            & wagon_f["valid_tag"]
            & (wagon_f["pred_probe_energy"] >= elo)
            & (wagon_f["pred_probe_energy"] < ehi)
            & (wagon_f["ep_missing_mass2"] >= mm2_min)
            & (wagon_f["ep_missing_mass2"] < mm2_max)
        )
        ns_base = (
            np.asarray(nsidis_parent_mask, dtype=bool)
            & nsidis_f["valid_tag"]
            & (nsidis_f["pred_probe_energy"] >= elo)
            & (nsidis_f["pred_probe_energy"] < ehi)
            & (nsidis_f["ep_missing_mass2"] >= mm2_min)
            & (nsidis_f["ep_missing_mass2"] < mm2_max)
        )

        wb = int(np.count_nonzero(wagon_base))
        nb = int(np.count_nonzero(ns_base))
        ob = int(np.count_nonzero(ns_base & nsidis_event_overlap))
        xb = int(np.count_nonzero(ns_base & ~nsidis_event_overlap))

        for cut in scan_values:
            cut = float(cut)
            wm = wagon_base & (
                np.abs(wagon_f["pred_probe_mass2"]) < cut
            )
            nm = ns_base & (
                np.abs(nsidis_f["pred_probe_mass2"]) < cut
            )
            om = nm & nsidis_event_overlap
            xm = nm & ~nsidis_event_overlap

            wc = int(np.count_nonzero(wm))
            nc = int(np.count_nonzero(nm))
            oc = int(np.count_nonzero(om))
            xc = int(np.count_nonzero(xm))

            rows.append({
                "period": period.key,
                "label": period.label,
                "energy_range": range_name,
                "energy_low_GeV": elo,
                "energy_high_GeV": ehi,
                "abs_pred_probe_m2_max_GeV2": cut,
                "wagon_retention": (
                    wc / wb if wb else float("nan")
                ),
                "nsidis_total_retention": (
                    nc / nb if nb else float("nan")
                ),
                "nsidis_wagon_event_retention": (
                    oc / ob if ob else float("nan")
                ),
                "nsidis_nonwagon_event_retention": (
                    xc / xb if xb else float("nan")
                ),
                "surviving_nsidis_wagon_event_fraction": (
                    oc / nc if nc else float("nan")
                ),
                "wagon_surviving_count": wc,
                "nsidis_surviving_count": nc,
                "nsidis_wagon_event_surviving_count": oc,
                "nsidis_nonwagon_event_surviving_count": xc,
            })
        #endfor
    #endfor

    return rows


def make_nsidis_probe_mass2_scan_canvas(
    period: PeriodConfig,
    rows: List[Dict[str, object]],
    outdir: Path,
) -> None:
    rr = sorted(
        [
            row for row in rows
            if row["energy_range"] == "overlap_0p4_2"
        ],
        key=lambda row: float(row["abs_pred_probe_m2_max_GeV2"]),
    )
    if not rr:
        return
    #endif

    x = [
        float(row["abs_pred_probe_m2_max_GeV2"])
        for row in rr
    ]

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2))
    axes[0].plot(
        x,
        [float(row["wagon_retention"]) for row in rr],
        marker="o",
        label=r"old wagon, $e_p>2$",
    )
    axes[0].plot(
        x,
        [
            float(row["nsidis_wagon_event_retention"])
            for row in rr
        ],
        marker="o",
        label="nSidis entries from wagon events",
    )
    axes[0].plot(
        x,
        [
            float(row["nsidis_nonwagon_event_retention"])
            for row in rr
        ],
        marker="o",
        label="nSidis entries from non-wagon events",
    )
    axes[0].set_ylim(-0.03, 1.03)
    axes[0].set_xlabel(
        r"$|M_X^2(ep\gamma_{\rm tag})|_{\max}$ (GeV$^2$)"
    )
    axes[0].set_ylabel("retention")
    axes[0].set_title(
        r"$0.4<E_{\rm probe}^{pred}<2$ GeV; "
        r"$M_X^2(ep)$ fit window already applied"
    )
    axes[0].grid(alpha=0.18)
    axes[0].legend(fontsize=7.5, frameon=False)

    axes[1].plot(
        x,
        [
            float(row["surviving_nsidis_wagon_event_fraction"])
            for row in rr
        ],
        marker="o",
    )
    axes[1].set_ylim(-0.03, 1.03)
    axes[1].set_xlabel(
        r"$|M_X^2(ep\gamma_{\rm tag})|_{\max}$ (GeV$^2$)"
    )
    axes[1].set_ylabel("wagon-event fraction of surviving nSidis entries")
    axes[1].set_title("exclusive-sample enrichment below 2 GeV")
    axes[1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: nSidis mass-shell support scan "
        "(pTmiss deliberately uncut)",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_probe_mass2_scan.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)



def profile_pi0_fraction_uncertainty_one_driver(
    hd: np.ndarray,
    hp: np.ndarray,
    hv: np.ndarray,
    fit: SharedMorphedFitResult,
    *,
    discriminator: str,
    nuisance_shift_prior: float,
    nuisance_sigma_prior: float,
    max_shift_bins: float,
    max_sigma_bins: float,
    delta_nll: float = 0.5,
) -> Dict[str, float]:
    """
    One-sigma profile-likelihood interval for f_pi0 in the one-driver pilot.

    All four template-morph nuisance parameters are re-profiled at each fixed
    f_pi0.  Since the production pilot uses only M_X^2(ep) x pTmiss,
    Delta NLL = 0.5 is the usual 68% one-parameter interval.
    """
    if (
        not fit.success
        or not np.isfinite(fit.pi0_fraction)
        or fit.nuisance is None
    ):
        return {
            "pi0_fraction_err_low": float("nan"),
            "pi0_fraction_err_high": float("nan"),
            "pi0_fraction_profile_low": float("nan"),
            "pi0_fraction_profile_high": float("nan"),
            "pi0_fraction_profile_status": "fit_unavailable",
        }
    #endif

    data = np.asarray(hd, dtype=float).ravel()
    raw_pi0 = np.asarray(hp, dtype=float)
    raw_dvcs = np.asarray(hv, dtype=float)
    ndata = float(np.sum(data))
    if ndata <= 0.0:
        return {
            "pi0_fraction_err_low": float("nan"),
            "pi0_fraction_err_high": float("nan"),
            "pi0_fraction_profile_low": float("nan"),
            "pi0_fraction_profile_high": float("nan"),
            "pi0_fraction_profile_status": "empty_data",
        }
    #endif

    keys = (
        f"{discriminator}_pi0_shift_bins",
        f"{discriminator}_pi0_sigma_bins",
        f"{discriminator}_dvcs_shift_bins",
        f"{discriminator}_dvcs_sigma_bins",
    )
    x_best = np.asarray(
        [float(fit.nuisance.get(key, 0.0)) for key in keys],
        dtype=float,
    )
    bounds = (
        (-float(max_shift_bins), float(max_shift_bins)),
        (0.0, float(max_sigma_bins)),
        (-float(max_shift_bins), float(max_shift_bins)),
        (0.0, float(max_sigma_bins)),
    )

    def objective_fixed_fraction(fraction: float, x: np.ndarray) -> float:
        f = float(np.clip(fraction, 3.0e-4, 1.0 - 3.0e-4))
        tp = morph_template_second_axis(
            raw_pi0, float(x[0]), float(x[1])
        ).ravel()
        td = morph_template_second_axis(
            raw_dvcs, float(x[2]), float(x[3])
        ).ravel()
        mu = np.clip(
            ndata * (f * tp + (1.0 - f) * td),
            1.0e-12,
            None,
        )
        penalty = (
            0.5 * (float(x[0]) / nuisance_shift_prior) ** 2
            + 0.5 * (float(x[1]) / nuisance_sigma_prior) ** 2
            + 0.5 * (float(x[2]) / nuisance_shift_prior) ** 2
            + 0.5 * (float(x[3]) / nuisance_sigma_prior) ** 2
        )
        return float(np.sum(mu - data * np.log(mu)) + penalty)
    #enddef

    cache: Dict[float, Tuple[float, np.ndarray]] = {}

    def profile_objective(fraction: float) -> float:
        f = float(np.clip(fraction, 3.0e-4, 1.0 - 3.0e-4))
        key = round(f, 10)
        if key in cache:
            return cache[key][0]
        #endif

        x0 = x_best
        if cache:
            nearest = min(cache, key=lambda kk: abs(kk - key))
            x0 = cache[nearest][1]
        #endif

        result = minimize(
            lambda xx: objective_fixed_fraction(f, xx),
            np.asarray(x0, dtype=float),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 180, "ftol": 2.0e-9, "maxls": 30},
        )
        if result.success and np.isfinite(result.fun):
            value = float(result.fun)
            state = np.asarray(result.x, dtype=float)
        else:
            value = objective_fixed_fraction(f, np.asarray(x0, dtype=float))
            state = np.asarray(x0, dtype=float)
        #endif

        cache[key] = (value, state)
        return value
    #enddef

    f_best = float(np.clip(fit.pi0_fraction, 3.0e-4, 1.0 - 3.0e-4))
    best_value = profile_objective(f_best)
    target = best_value + float(delta_nll)

    def root_fn(f: float) -> float:
        return profile_objective(f) - target
    #enddef

    f_min = 3.0e-4
    f_max = 1.0 - 3.0e-4
    low_status = "ok"
    high_status = "ok"

    if root_fn(f_min) <= 0.0:
        f_low = f_min
        low_status = "lower_interval_hits_boundary"
    else:
        f_low = float(brentq(
            root_fn, f_min, f_best,
            xtol=2.0e-5, rtol=2.0e-5, maxiter=60,
        ))
    #endif

    if root_fn(f_max) <= 0.0:
        f_high = f_max
        high_status = "upper_interval_hits_boundary"
    else:
        f_high = float(brentq(
            root_fn, f_best, f_max,
            xtol=2.0e-5, rtol=2.0e-5, maxiter=60,
        ))
    #endif

    return {
        "pi0_fraction_err_low": max(0.0, f_best - f_low),
        "pi0_fraction_err_high": max(0.0, f_high - f_best),
        "pi0_fraction_profile_low": f_low,
        "pi0_fraction_profile_high": f_high,
        "pi0_fraction_profile_status": (
            f"{low_status};{high_status}"
        ),
    }




NOMINAL_THREE_VARIABLES = ("Delta_phi", "pTmiss", "xF2")
NOMINAL_THREE_RANGES = {
    "Delta_phi": (math.pi - 0.40, math.pi + 0.40, 80),
    "pTmiss": (0.0, 1.0, 80),
    "xF2": (0.00, 1.20, 96),
}


def _nominal_three_values(feat, variable):
    if variable == "Delta_phi":
        return np.asarray(feat["stored_delta_phi_rad"], dtype=float)
    #endif
    if variable == "pTmiss":
        return np.asarray(feat["stored_pTmiss"], dtype=float)
    #endif
    if variable == "xF2":
        return np.asarray(feat["stored_xF2"], dtype=float)
    #endif
    raise KeyError(variable)


def _nominal_three_hist_values(
    feat,
    variable,
):
    """
    Values used for the nominal 1D composition histograms.

    No clipping or under/overflow folding is performed for xF2.  The raw ROOT
    branch values are histogrammed directly.  Values outside the displayed
    histogram interval are therefore not collapsed into artificial edge-bin
    spikes.
    """
    return _nominal_three_values(
        feat,
        variable,
    )




def _fixed_morph_fraction_error(histograms, fit):
    if (
        not fit.success
        or fit.nuisance is None
        or not np.isfinite(fit.pi0_fraction)
    ):
        return float("nan"), float("nan")
    #endif

    names = tuple(histograms.keys())

    def obj(f):
        f = float(np.clip(f, 3e-4, 1-3e-4))
        total = 0.0
        for name in names:
            hd, hp, hv = histograms[name]
            tp = morph_template_second_axis(
                hp,
                fit.nuisance.get(f"{name}_pi0_shift_bins", 0.0),
                fit.nuisance.get(f"{name}_pi0_sigma_bins", 0.0),
            ).ravel()
            td = morph_template_second_axis(
                hv,
                fit.nuisance.get(f"{name}_dvcs_shift_bins", 0.0),
                fit.nuisance.get(f"{name}_dvcs_sigma_bins", 0.0),
            ).ravel()
            d = np.asarray(hd, dtype=float).ravel()
            nd = float(np.sum(d))
            mu = np.clip(nd*(f*tp + (1-f)*td), 1e-12, None)
            total += float(np.sum(mu - d*np.log(mu))) / len(names)
        #endfor
        return total

    f0 = float(fit.pi0_fraction)
    target = obj(f0) + 0.5
    grid = np.linspace(3e-4, 1-3e-4, 501)
    vals = np.asarray([obj(x) for x in grid])
    allowed = grid[vals <= target]
    if allowed.size == 0:
        return float("nan"), float("nan")
    #endif
    return (
        max(0.0, f0 - float(np.min(allowed))),
        max(0.0, float(np.max(allowed)) - f0),
    )



@dataclass
class TwoCategoryCompositionResult:
    success: bool
    message: str
    pi0_fraction: float = float("nan")
    err_low: float = float("nan")
    err_high: float = float("nan")
    pi0_high_eff: float = float("nan")
    dvcs_high_eff: float = float("nan")
    data_high_fraction: float = float("nan")
    nll: float = float("nan")


def _binomial_nll(k: int, n: int, p: float) -> float:
    """Binomial NLL, dropping only the parameter-independent combinatorial term."""
    if n <= 0:
        return 0.0
    #endif
    p = float(np.clip(p, 1.0e-10, 1.0 - 1.0e-10))
    return -(
        float(k) * math.log(p)
        + float(n - k) * math.log1p(-p)
    )


def fit_two_category_composition(
    data_low: int,
    data_high: int,
    pi0_low: int,
    pi0_high: int,
    dvcs_low: int,
    dvcs_high: int,
) -> TwoCategoryCompositionResult:
    """
    Fit one pi0 fraction from low/high SIDIS-pT categories.

    The data category probability is
        q_high = f_pi0*eps_pi0_high + (1-f_pi0)*eps_dvcs_high.

    Finite aaogen and dvcsgen statistics are included by treating the two MC
    category efficiencies as nuisance parameters constrained by their own
    binomial likelihoods.
    """
    from scipy.optimize import minimize, minimize_scalar

    n_data = int(data_low + data_high)
    n_pi0 = int(pi0_low + pi0_high)
    n_dvcs = int(dvcs_low + dvcs_high)

    if min(n_data, n_pi0, n_dvcs) <= 0:
        return TwoCategoryCompositionResult(
            False,
            "empty data/template category total",
        )
    #endif

    obs_data = float(data_high) / float(n_data)
    obs_pi0 = float(pi0_high) / float(n_pi0)
    obs_dvcs = float(dvcs_high) / float(n_dvcs)

    if abs(obs_pi0 - obs_dvcs) < 1.0e-6:
        return TwoCategoryCompositionResult(
            False,
            "pi0 and DVCS category efficiencies are indistinguishable",
            pi0_high_eff=obs_pi0,
            dvcs_high_eff=obs_dvcs,
            data_high_fraction=obs_data,
        )
    #endif

    f_seed = float(np.clip(
        (obs_data - obs_dvcs) / (obs_pi0 - obs_dvcs),
        1.0e-5,
        1.0 - 1.0e-5,
    ))

    def full_nll(x):
        f, eps_pi0, eps_dvcs = map(float, x)
        q = f * eps_pi0 + (1.0 - f) * eps_dvcs
        return (
            _binomial_nll(data_high, n_data, q)
            + _binomial_nll(pi0_high, n_pi0, eps_pi0)
            + _binomial_nll(dvcs_high, n_dvcs, eps_dvcs)
        )

    result = minimize(
        full_nll,
        x0=np.asarray([f_seed, obs_pi0, obs_dvcs], dtype=float),
        method="L-BFGS-B",
        bounds=(
            (1.0e-8, 1.0 - 1.0e-8),
            (1.0e-8, 1.0 - 1.0e-8),
            (1.0e-8, 1.0 - 1.0e-8),
        ),
        options={"maxiter": 500, "ftol": 1.0e-12},
    )

    if not result.success or not np.all(np.isfinite(result.x)):
        return TwoCategoryCompositionResult(
            False,
            f"optimizer failed: {result.message}",
            pi0_high_eff=obs_pi0,
            dvcs_high_eff=obs_dvcs,
            data_high_fraction=obs_data,
        )
    #endif

    fhat, eps_pi0_hat, eps_dvcs_hat = map(float, result.x)
    nll_min = float(result.fun)

    def profile_nll(fixed_f: float) -> float:
        fixed_f = float(np.clip(fixed_f, 1.0e-8, 1.0 - 1.0e-8))

        def nuisance_nll(y):
            eps_pi0, eps_dvcs = map(float, y)
            q = fixed_f * eps_pi0 + (1.0 - fixed_f) * eps_dvcs
            return (
                _binomial_nll(data_high, n_data, q)
                + _binomial_nll(pi0_high, n_pi0, eps_pi0)
                + _binomial_nll(dvcs_high, n_dvcs, eps_dvcs)
            )

        rr = minimize(
            nuisance_nll,
            x0=np.asarray([eps_pi0_hat, eps_dvcs_hat], dtype=float),
            method="L-BFGS-B",
            bounds=(
                (1.0e-8, 1.0 - 1.0e-8),
                (1.0e-8, 1.0 - 1.0e-8),
            ),
            options={"maxiter": 300, "ftol": 1.0e-12},
        )
        return float(rr.fun) if rr.success else float("inf")

    target = nll_min + 0.5
    grid = np.linspace(1.0e-5, 1.0 - 1.0e-5, 201)
    prof = np.asarray([profile_nll(f) for f in grid], dtype=float)
    allowed = grid[np.isfinite(prof) & (prof <= target)]

    if allowed.size:
        flo = float(np.min(allowed))
        fhi = float(np.max(allowed))

        if flo < fhat:
            rr = minimize_scalar(
                lambda f: abs(profile_nll(f) - target),
                bounds=(max(1.0e-8, flo - 0.01), fhat),
                method="bounded",
                options={"xatol": 1.0e-5},
            )
            if rr.success:
                flo = float(rr.x)
            #endif
        #endif

        if fhi > fhat:
            rr = minimize_scalar(
                lambda f: abs(profile_nll(f) - target),
                bounds=(fhat, min(1.0 - 1.0e-8, fhi + 0.01)),
                method="bounded",
                options={"xatol": 1.0e-5},
            )
            if rr.success:
                fhi = float(rr.x)
            #endif
        #endif

        err_low = max(0.0, fhat - flo)
        err_high = max(0.0, fhi - fhat)
    else:
        err_low = err_high = float("nan")
    #endif

    return TwoCategoryCompositionResult(
        True,
        "ok",
        pi0_fraction=fhat,
        err_low=err_low,
        err_high=err_high,
        pi0_high_eff=eps_pi0_hat,
        dvcs_high_eff=eps_dvcs_hat,
        data_high_fraction=obs_data,
        nll=nll_min,
    )



def run_nsidis_three_variable_nominal_fits(
    period,
    data_f,
    pi0_f,
    dvcs_f,
    *,
    ft_theta_max,
    max_probe_energy,
    mm2_min,
    mm2_max,
    probe_m2_max,
    min_data_count,
    min_template_count,
    nuisance_shift_prior,
    nuisance_sigma_prior,
    max_shift_bins,
    max_sigma_bins,
    source_label,
    pt_category_split=0.10,
    parent_mask=None,
):
    rows = []

    for name, feat in (
        ("data", data_f),
        ("pi0", pi0_f),
        ("dvcs", dvcs_f),
    ):
        for key in (
            "stored_delta_phi_rad",
            "stored_pTmiss",
            "stored_xF2",
            "electron_p",
        ):
            if key not in feat:
                raise RuntimeError(
                    f"{period.label} {source_label} {name}: missing {key}"
                )
            #endif
        #endfor
    #endfor

    for region in ("FT", "FD_all"):
        edges = nsidis_energy_edges_for_region(region, max_probe_energy)
        for ib in range(len(edges)-1):
            elo, ehi = float(edges[ib]), float(edges[ib+1])
            masks = {}
            for name, feat in (
                ("data", data_f),
                ("pi0", pi0_f),
                ("dvcs", dvcs_f),
            ):
                m = stage2_fit_mask(
                    feat, region, ft_theta_max, elo, ehi,
                    mm2_min, mm2_max, probe_m2_max,
                )
                m &= np.isfinite(feat["electron_p"]) & (feat["electron_p"] > 2.0)
                for var in NOMINAL_THREE_VARIABLES:
                    lo, hi, _ = NOMINAL_THREE_RANGES[var]
                    v = _nominal_three_values(feat, var)

                    if var == "xF2":
                        # xF2 is a discriminator, not an event-selection cut.
                        # Keep every finite raw branch value. The histogram
                        # displays only its configured range; no clipping or
                        # edge-bin folding is applied.
                        m &= np.isfinite(v)
                    else:
                        m &= (
                            np.isfinite(v)
                            & (v >= lo)
                            & (v < hi)
                        )
                    #endif
                #endfor
                masks[name] = m
            #endfor

            if parent_mask is not None:
                pm = np.asarray(parent_mask, dtype=bool)
                if pm.shape != masks["data"].shape:
                    raise ValueError("parent_mask shape mismatch")
                #endif
                masks["data"] &= pm
            #endif

            hists = {}
            for var in NOMINAL_THREE_VARIABLES:
                lo, hi, nb = NOMINAL_THREE_RANGES[var]
                be = np.linspace(lo, hi, nb+1)
                hh = []
                for feat, m in (
                    (data_f, masks["data"]),
                    (pi0_f, masks["pi0"]),
                    (dvcs_f, masks["dvcs"]),
                ):
                    h, _ = np.histogram(
                        _nominal_three_hist_values(
                            feat,
                            var,
                        )[m],
                        bins=be,
                    )
                    hh.append(h.astype(float)[None, :])
                #endfor
                hists[var] = tuple(hh)
            #endfor

            counts = {k: int(np.count_nonzero(v)) for k, v in masks.items()}
            enough = (
                counts["data"] >= min_data_count
                and counts["pi0"] >= min_template_count
                and counts["dvcs"] >= min_template_count
            )

            if enough:
                fit = fit_shared_morphed_composition(
                    hists, {},
                    nuisance_shift_prior, nuisance_sigma_prior,
                    max_shift_bins, max_sigma_bins,
                    driver_names=NOMINAL_THREE_VARIABLES,
                )
            else:
                fit = SharedMorphedFitResult(
                    False,
                    f"insufficient stats data={counts['data']} pi0={counts['pi0']} dvcs={counts['dvcs']}",
                )
            #endif

            singles = {}
            if enough:
                for var in NOMINAL_THREE_VARIABLES:
                    sf = fit_shared_morphed_composition(
                        {var: hists[var]}, {},
                        nuisance_shift_prior, nuisance_sigma_prior,
                        max_shift_bins, max_sigma_bins,
                        driver_names=(var,),
                    )
                    singles[var] = float(sf.pi0_fraction) if sf.success else float("nan")
                #endfor
            #endif

            f0 = float(fit.pi0_fraction) if fit.success else float("nan")
            usable_single_fractions = {
                key: value
                for key, value in singles.items()
                if (
                    np.isfinite(value)
                    and value > 1.0e-3
                    and value < 1.0 - 1.0e-3
                )
            }
            vals = np.asarray(
                list(usable_single_fractions.values()),
                dtype=float,
            )
            if np.isfinite(f0) and vals.size:
                diffs = vals - f0
                spread_max = float(np.max(np.abs(diffs)))
                spread_rms = float(np.sqrt(np.mean(diffs*diffs)))
            else:
                spread_max = spread_rms = float("nan")
            #endif

            err_lo, err_hi = _fixed_morph_fraction_error(hists, fit)

            pt_category = TwoCategoryCompositionResult(
                False,
                "SIDIS pT branch unavailable",
            )
            pt_category_counts = {
                "data_low": 0,
                "data_high": 0,
                "pi0_low": 0,
                "pi0_high": 0,
                "dvcs_low": 0,
                "dvcs_high": 0,
            }

            if all(
                "stored_sidis_pT" in feat
                for feat in (data_f, pi0_f, dvcs_f)
            ):
                for sample_name, feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                ):
                    pt_values = np.asarray(
                        feat["stored_sidis_pT"],
                        dtype=float,
                    )
                    base = (
                        masks[sample_name]
                        & np.isfinite(pt_values)
                    )
                    low = base & (
                        pt_values < float(pt_category_split)
                    )
                    high = base & (
                        pt_values >= float(pt_category_split)
                    )
                    pt_category_counts[f"{sample_name}_low"] = int(
                        np.count_nonzero(low)
                    )
                    pt_category_counts[f"{sample_name}_high"] = int(
                        np.count_nonzero(high)
                    )
                #endfor

                pt_category = fit_two_category_composition(
                    pt_category_counts["data_low"],
                    pt_category_counts["data_high"],
                    pt_category_counts["pi0_low"],
                    pt_category_counts["pi0_high"],
                    pt_category_counts["dvcs_low"],
                    pt_category_counts["dvcs_high"],
                )
            #endif

            # Raw xF2 range diagnostics. These counts are informational only:
            # xF2 is not used as an event-selection cut. The histogram itself
            # uses the visible 0.0--1.2 interval without clipping/folding.
            xF2_lo, xF2_hi, _ = NOMINAL_THREE_RANGES["xF2"]
            xF2_values_data = _nominal_three_values(
                data_f,
                "xF2",
            )
            xF2_selected = xF2_values_data[
                masks["data"]
            ]
            xF2_underflow_count = int(
                np.count_nonzero(
                    np.isfinite(xF2_selected)
                    & (xF2_selected < xF2_lo)
                )
            )
            xF2_overflow_count = int(
                np.count_nonzero(
                    np.isfinite(xF2_selected)
                    & (xF2_selected >= xF2_hi)
                )
            )

            rows.append({
                "period": period.key,
                "label": period.label,
                "source": source_label,
                "region": region,
                "energy_low_GeV": elo,
                "energy_high_GeV": ehi,
                "energy_center_GeV": 0.5*(elo+ehi),
                "mean_probe_energy_GeV": masked_finite_mean(
                    data_f["pred_probe_energy"], masks["data"]
                ),
                "probe_m2_support_max_GeV2": float(probe_m2_max),
                "fit_success": int(fit.success),
                "fit_message": fit.message,
                "fit_model": "simultaneous_Delta_phi_pTmiss_xF2",
                "xF2_selection_cut": "none",
                "xF2_histogram_range": "0.00_to_1.20",
                "xF2_data_underflow_count": xF2_underflow_count,
                "xF2_data_overflow_count": xF2_overflow_count,
                "xF2_histogram_overflow_policy": "none_raw_branch_values",
                "pi0_fraction": f0,
                "dvcs_fraction": 1.0-f0 if np.isfinite(f0) else float("nan"),
                "pi0_fraction_err_low": err_lo,
                "pi0_fraction_err_high": err_hi,
                "pi0_fraction_single_Delta_phi": singles.get("Delta_phi", float("nan")),
                "pi0_fraction_single_pTmiss": singles.get("pTmiss", float("nan")),
                "pi0_fraction_single_xF2": singles.get("xF2", float("nan")),
                "pi0_fraction_single_Delta_phi_usable": int(
                    "Delta_phi" in usable_single_fractions
                ),
                "pi0_fraction_single_pTmiss_usable": int(
                    "pTmiss" in usable_single_fractions
                ),
                "pi0_fraction_single_xF2_usable": int(
                    "xF2" in usable_single_fractions
                ),
                "pi0_fraction_model_max_abs_shift": spread_max,
                "pi0_fraction_model_rms_shift": spread_rms,
                "pt_category_split_GeV": float(pt_category_split),
                "pt_category_fit_success": int(pt_category.success),
                "pt_category_fit_message": pt_category.message,
                "pt_category_pi0_fraction": float(pt_category.pi0_fraction),
                "pt_category_pi0_fraction_err_low": float(pt_category.err_low),
                "pt_category_pi0_fraction_err_high": float(pt_category.err_high),
                "pt_category_pi0_high_eff": float(pt_category.pi0_high_eff),
                "pt_category_dvcs_high_eff": float(pt_category.dvcs_high_eff),
                "pt_category_data_high_fraction": float(
                    pt_category.data_high_fraction
                ),
                "pt_category_template_separation": (
                    float(
                        pt_category.pi0_high_eff
                        - pt_category.dvcs_high_eff
                    )
                    if (
                        np.isfinite(pt_category.pi0_high_eff)
                        and np.isfinite(pt_category.dvcs_high_eff)
                    )
                    else float("nan")
                ),
                **{
                    f"pt_category_{key}_count": int(value)
                    for key, value in pt_category_counts.items()
                },
                **(
                    {
                        key: float(value)
                        for key, value in fit.nuisance.items()
                    }
                    if fit.success and fit.nuisance is not None
                    else {}
                ),
                "data_candidate_count": counts["data"],
                "aaogen_candidate_count": counts["pi0"],
                "dvcsgen_candidate_count": counts["dvcs"],
                "poisson_deviance": float(fit.poisson_deviance) if fit.success else float("nan"),
                "ndof": int(fit.ndof) if fit.success else 0,
                "deviance_per_ndof": (
                    float(fit.poisson_deviance/fit.ndof)
                    if fit.success and fit.ndof else float("nan")
                ),
                "status": "nominal" if fit.success else "fit_failed",
            })
        #endfor
    #endfor

    return rows


def run_nsidis_ptmiss_pilot_fits(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    *,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    support_values: Sequence[float],
    mm2_bins: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    min_data_count: int,
    min_template_count: int,
    nuisance_shift_prior: float,
    nuisance_sigma_prior: float,
    max_shift_bins: float,
    max_sigma_bins: float,
    source_label: str,
    parent_mask: Optional[np.ndarray] = None,
    uncertainty_support: Optional[float] = None,
) -> List[Dict[str, object]]:
    """
    Diagnostic production-model pilot using only
        M_X^2(ep) x pTmiss.

    Each mass-shell support value is fit independently in each region/energy bin.
    """
    rows: List[Dict[str, object]] = []
    disc = "mx2_ep_x_pTmiss"

    for support in support_values:
        support = float(support)

        for region in ("FT", "FD_all"):
            edges = nsidis_energy_edges_for_region(
                region, max_probe_energy
            )
            for ib in range(len(edges) - 1):
                elo = float(edges[ib])
                ehi = float(edges[ib + 1])

                masks = {
                    name: stage2_fit_mask(
                        feat,
                        region,
                        ft_theta_max,
                        elo,
                        ehi,
                        mm2_min,
                        mm2_max,
                        support,
                    )
                    for name, feat in (
                        ("data", data_f),
                        ("pi0", pi0_f),
                        ("dvcs", dvcs_f),
                    )
                }
                # The loose nSidis parent skim requires e_p > 2 GeV.
                # Apply the same parent support to data AND both templates.
                # parent_mask is retained only to reproduce the exact selected
                # data/wagon entry population used by overlap diagnostics.
                for _name, _feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                ):
                    if "electron_p" not in _feat:
                        raise KeyError(
                            f"{source_label}: {_name} feature store lacks electron_p"
                        )
                    #endif
                    masks[_name] &= (
                        np.isfinite(_feat["electron_p"])
                        & (_feat["electron_p"] > 2.0)
                    )
                #endfor

                if parent_mask is not None:
                    pmask = np.asarray(parent_mask, dtype=bool)
                    if pmask.shape != masks["data"].shape:
                        raise ValueError(
                            "pilot parent_mask has incompatible data shape"
                        )
                    #endif
                    masks["data"] &= pmask
                #endif

                hists = {}
                for name, feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                ):
                    hists[name] = histogram_for_discriminator(
                        disc,
                        feat,
                        masks[name],
                        mm2_min=mm2_min,
                        mm2_max=mm2_max,
                        probe_m2_max=support,
                        mm2_bins_2d=mm2_bins,
                        probe_m2_bins_2d=48,
                        bins_1d=90,
                        ptmiss_max=ptmiss_max,
                        ptmiss_bins=ptmiss_bins,
                        theta_max=4.0,
                        theta_bins=80,
                    )
                #endfor

                row = {
                    "period": period.key,
                    "label": period.label,
                    "source": source_label,
                    "region": region,
                    "energy_low_GeV": elo,
                    "energy_high_GeV": ehi,
                    "energy_center_GeV": 0.5 * (elo + ehi),
                    "mean_probe_energy_GeV": masked_finite_mean(
                        data_f["pred_probe_energy"], masks["data"]
                    ),
                    "data_mean_probe_energy_GeV": masked_finite_mean(
                        data_f["pred_probe_energy"], masks["data"]
                    ),
                    "aaogen_mean_probe_energy_GeV": masked_finite_mean(
                        pi0_f["pred_probe_energy"], masks["pi0"]
                    ),
                    "dvcsgen_mean_probe_energy_GeV": masked_finite_mean(
                        dvcs_f["pred_probe_energy"], masks["dvcs"]
                    ),
                    "probe_m2_support_max_GeV2": support,
                    "data_candidate_count": int(np.sum(hists["data"])),
                    "aaogen_candidate_count": int(np.sum(hists["pi0"])),
                    "dvcsgen_candidate_count": int(np.sum(hists["dvcs"])),
                    "status": "diagnostic_only_not_efficiency",
                    "fit_model": STAGE2_PRODUCTION_MODEL,
                }

                if (
                    np.sum(hists["data"]) < min_data_count
                    or np.sum(hists["pi0"]) < min_template_count
                    or np.sum(hists["dvcs"]) < min_template_count
                ):
                    row.update({
                        "fit_success": 0,
                        "pi0_fraction": float("nan"),
                        "deviance_per_ndof": float("nan"),
                        "fit_message": "insufficient statistics",
                    })
                    rows.append(row)
                    continue
                #endif

                fit = fit_shared_morphed_composition(
                    {
                        disc: (
                            hists["data"],
                            hists["pi0"],
                            hists["dvcs"],
                        )
                    },
                    {},
                    nuisance_shift_prior,
                    nuisance_sigma_prior,
                    max_shift_bins,
                    max_sigma_bins,
                    driver_names=(disc,),
                )

                row.update({
                    "fit_success": int(fit.success),
                    "fit_message": fit.message,
                    "pi0_fraction": float(fit.pi0_fraction),
                    "dvcs_fraction": (
                        1.0 - float(fit.pi0_fraction)
                        if np.isfinite(fit.pi0_fraction)
                        else float("nan")
                    ),
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

                # The legacy fit.poisson_deviance/ndof is retained in the CSV
                # for backward comparison, but no longer treated as the primary
                # goodness metric.  Compute quality measures that correspond
                # directly to the actual fitted 2D model and its projections.
                row["legacy_sparse_2d_deviance_per_ndof"] = (
                    float(fit.poisson_deviance / fit.ndof)
                    if fit.ndof else float("nan")
                )
                row.update(
                    fit_quality_for_morphed_2d_model(
                        hists["data"],
                        hists["pi0"],
                        hists["dvcs"],
                        fit,
                        disc,
                    )
                )

                if (
                    uncertainty_support is not None
                    and abs(
                        support - float(uncertainty_support)
                    ) < 1.0e-12
                    and fit.success
                ):
                    row.update(
                        profile_pi0_fraction_uncertainty_one_driver(
                            hists["data"],
                            hists["pi0"],
                            hists["dvcs"],
                            fit,
                            discriminator=disc,
                            nuisance_shift_prior=nuisance_shift_prior,
                            nuisance_sigma_prior=nuisance_sigma_prior,
                            max_shift_bins=max_shift_bins,
                            max_sigma_bins=max_sigma_bins,
                        )
                    )
                else:
                    row.update({
                        "pi0_fraction_err_low": float("nan"),
                        "pi0_fraction_err_high": float("nan"),
                        "pi0_fraction_profile_low": float("nan"),
                        "pi0_fraction_profile_high": float("nan"),
                        "pi0_fraction_profile_status": (
                            "not_evaluated_noncentral_support"
                        ),
                    })
                #endif

                rows.append(row)
            #endfor
        #endfor
    #endfor

    return rows


def make_nsidis_pilot_summary_canvas(
    period: PeriodConfig,
    wagon_rows: List[Dict[str, object]],
    nsidis_rows: List[Dict[str, object]],
    central_support: float,
    outdir: Path,
) -> None:
    """
    Summary of overlap agreement, high-energy extension, support stability, and
    projection-level fit quality.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.8))

    for irow, region in enumerate(("FT", "FD_all")):
        ax_frac = axes[irow, 0]
        ax_dev = axes[irow, 1]

        supports = sorted({
            float(row["probe_m2_support_max_GeV2"])
            for row in nsidis_rows
        })
        for support in supports:
            rr = sorted(
                [
                    row for row in nsidis_rows
                    if row["region"] == region
                    and int(row.get("fit_success", 0)) == 1
                    and abs(
                        float(row["probe_m2_support_max_GeV2"]) - support
                    ) < 1.0e-12
                ],
                key=row_energy_coordinate,
            )
            if not rr:
                continue
            #endif

            x = [row_energy_coordinate(row) for row in rr]
            ax_frac.plot(
                x,
                [float(row["pi0_fraction"]) for row in rr],
                marker="o",
                linewidth=1.0,
                label=rf"nSidis $|M^2_{{probe}}|<{support:.2f}$",
            )
            ax_dev.plot(
                x,
                [float(row.get("quality_ptmiss_projection_deviance_per_active_bin", float("nan"))) for row in rr],
                marker="o",
                linewidth=1.0,
                label=rf"nSidis {support:.2f}",
            )
        #endfor

        wr = sorted(
            [
                row for row in wagon_rows
                if row["region"] == region
                and int(row.get("fit_success", 0)) == 1
                and abs(
                    float(row["probe_m2_support_max_GeV2"])
                    - float(central_support)
                ) < 1.0e-12
            ],
            key=row_energy_coordinate,
        )
        if wr:
            x = [float(row["energy_center_GeV"]) for row in wr]
            ax_frac.plot(
                x,
                [float(row["pi0_fraction"]) for row in wr],
                marker="s",
                linestyle="--",
                linewidth=1.2,
                label="old wagon, central support",
            )
            ax_dev.plot(
                x,
                [float(row.get("quality_ptmiss_projection_deviance_per_active_bin", float("nan"))) for row in wr],
                marker="s",
                linestyle="--",
                linewidth=1.2,
                label="old wagon",
            )
        #endif

        ax_frac.axvline(2.0, linestyle="--", linewidth=0.9)
        ax_dev.axvline(2.0, linestyle="--", linewidth=0.9)
        ax_frac.set_ylim(0.0, 1.02)
        ax_frac.set_ylabel(r"$f_{\pi^0}$")
        ax_dev.set_ylabel(r"$p_T$ projection deviance / active bin")
        ax_frac.set_title(f"{region}: composition")
        ax_dev.set_title(f"{region}: visible pT-projection goodness")
        for ax in (ax_frac, ax_dev):
            ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
            ax.grid(alpha=0.18)
        #endfor
    #endfor

    axes[0, 0].legend(fontsize=7.2, frameon=False)
    axes[0, 1].legend(fontsize=7.2, frameon=False)

    fig.suptitle(
        f"{period.label}: nSidis pilot with production "
        r"$M_X^2(ep)\otimes p_{T,\rm miss}$ model",
        fontsize=13.0,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_pilot_summary.png",
        rect=(0, 0, 1, 0.95),
    )
    plt.close(fig)


def make_nsidis_pilot_fit_projection_canvases(
    period: PeriodConfig,
    rows: List[Dict[str, object]],
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    outdir: Path,
    *,
    central_support: float,
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    mm2_bins: int,
    ptmiss_max: float,
    ptmiss_bins: int,
) -> None:
    """
    Show the ACTUAL central-support nSidis pilot fit in every energy bin.

    Layout is optimized for inspection:
      - up to four ENERGY BINS across each figure row;
      - within each energy-bin tile, M_X^2(ep) is the upper half and pTmiss
        is the lower half;
      - both panels are projections of the SAME fitted 2D model.
    """
    disc = "mx2_ep_x_pTmiss"
    row_map = {
        (
            str(row["region"]),
            round(float(row["energy_low_GeV"]), 10),
            round(float(row["energy_high_GeV"]), 10),
        ): row
        for row in rows
        if int(row.get("fit_success", 0)) == 1
        and abs(
            float(row["probe_m2_support_max_GeV2"])
            - float(central_support)
        ) < 1.0e-12
    }

    for region in ("FT", "FD_all"):
        edges = nsidis_energy_edges_for_region(region, 9.5)
        nbins = len(edges) - 1
        ncols = min(4, max(1, nbins))
        ntile_rows = int(math.ceil(nbins / ncols))

        fig = plt.figure(
            figsize=(4.8 * ncols, 5.8 * ntile_rows)
        )
        outer = fig.add_gridspec(
            ntile_rows,
            ncols,
            wspace=0.25,
            hspace=0.42,
        )

        first_axes = None

        for ib in range(nbins):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])
            tile_row = ib // ncols
            tile_col = ib % ncols
            inner = outer[tile_row, tile_col].subgridspec(
                2, 1,
                hspace=0.18,
            )
            ax_mx = fig.add_subplot(inner[0, 0])
            ax_pt = fig.add_subplot(inner[1, 0])
            if first_axes is None:
                first_axes = ax_mx
            #endif

            row = row_map.get(
                (region, round(elo, 10), round(ehi, 10))
            )
            if row is None:
                for ax in (ax_mx, ax_pt):
                    ax.text(
                        0.5, 0.5, "no successful fit",
                        transform=ax.transAxes,
                        ha="center", va="center",
                    )
                    ax.set_axis_off()
                #endfor
                continue
            #endif

            masks = {
                name: stage2_fit_mask(
                    feat,
                    region,
                    ft_theta_max,
                    elo,
                    ehi,
                    mm2_min,
                    mm2_max,
                    central_support,
                )
                for name, feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                )
            }
            for _name, _feat in (
                ("data", data_f),
                ("pi0", pi0_f),
                ("dvcs", dvcs_f),
            ):
                masks[_name] &= (
                    np.isfinite(_feat["electron_p"])
                    & (_feat["electron_p"] > 2.0)
                )
            #endfor

            hists = {}
            for name, feat in (
                ("data", data_f),
                ("pi0", pi0_f),
                ("dvcs", dvcs_f),
            ):
                hists[name] = histogram_for_discriminator(
                    disc,
                    feat,
                    masks[name],
                    mm2_min=mm2_min,
                    mm2_max=mm2_max,
                    probe_m2_max=central_support,
                    mm2_bins_2d=mm2_bins,
                    probe_m2_bins_2d=48,
                    bins_1d=90,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=4.0,
                    theta_bins=80,
                )
            #endfor

            hd = hists["data"]
            hp = hists["pi0"]
            hv = hists["dvcs"]

            tp = morph_template_second_axis(
                hp,
                float(row.get(f"{disc}_pi0_shift_bins", 0.0)),
                float(row.get(f"{disc}_pi0_sigma_bins", 0.0)),
            )
            td = morph_template_second_axis(
                hv,
                float(row.get(f"{disc}_dvcs_shift_bins", 0.0)),
                float(row.get(f"{disc}_dvcs_sigma_bins", 0.0)),
            )

            fpi0 = float(row["pi0_fraction"])
            ndata = float(np.sum(hd))
            rawp = normalized_template(hp).reshape(hp.shape)
            rawd = normalized_template(hv).reshape(hv.shape)

            for axis_name, ax in (
                ("mx2", ax_mx),
                ("ptmiss", ax_pt),
            ):
                if axis_name == "mx2":
                    data_proj = np.sum(hd, axis=1)
                    pshape = np.sum(tp, axis=1)
                    dshape = np.sum(td, axis=1)
                    rawp_shape = np.sum(rawp, axis=1)
                    rawd_shape = np.sum(rawd, axis=1)
                    plot_edges = np.linspace(
                        mm2_min, mm2_max, mm2_bins + 1
                    )
                    xlabel = r"$M_X^2(ep)$ (GeV$^2$)"
                    qkey = (
                        "quality_mx2_projection_deviance_per_active_bin"
                    )
                    ckey = (
                        "quality_mx2_projection_pearson_chi2_per_bin_mu_ge5"
                    )
                else:
                    data_proj = np.sum(hd, axis=0)
                    pshape = np.sum(tp, axis=0)
                    dshape = np.sum(td, axis=0)
                    rawp_shape = np.sum(rawp, axis=0)
                    rawd_shape = np.sum(rawd, axis=0)
                    plot_edges = np.linspace(
                        0.0, ptmiss_max, ptmiss_bins + 1
                    )
                    xlabel = r"$p_{T,\rm miss}$ (GeV)"
                    qkey = (
                        "quality_ptmiss_projection_deviance_per_active_bin"
                    )
                    ckey = (
                        "quality_ptmiss_projection_pearson_chi2_per_bin_mu_ge5"
                    )
                #endif

                for arr in (
                    pshape,
                    dshape,
                    rawp_shape,
                    rawd_shape,
                ):
                    arr /= max(float(np.sum(arr)), 1.0e-30)
                #endfor

                pcomp = ndata * fpi0 * pshape
                dcomp = ndata * (1.0 - fpi0) * dshape
                total = pcomp + dcomp
                rawp_comp = ndata * fpi0 * rawp_shape
                rawd_comp = ndata * (1.0 - fpi0) * rawd_shape

                centers = 0.5 * (
                    plot_edges[:-1] + plot_edges[1:]
                )
                ax.errorbar(
                    centers,
                    data_proj,
                    yerr=np.sqrt(np.maximum(data_proj, 1.0)),
                    fmt="o",
                    ms=1.8,
                    linewidth=0.55,
                    color=SAMPLE_COLORS["data"],
                    label="data",
                    zorder=6,
                )
                ax.step(
                    centers,
                    rawd_comp,
                    where="mid",
                    color=SAMPLE_COLORS["dvcs_mc"],
                    linewidth=0.55,
                    linestyle="--",
                    alpha=0.5,
                    label="BH/DVCS pre-morph",
                )
                ax.step(
                    centers,
                    rawp_comp,
                    where="mid",
                    color=SAMPLE_COLORS["pi0_mc"],
                    linewidth=0.55,
                    linestyle="--",
                    alpha=0.5,
                    label=r"$\pi^0$ pre-morph",
                )
                ax.step(
                    centers,
                    dcomp,
                    where="mid",
                    color=SAMPLE_COLORS["dvcs_mc"],
                    linewidth=1.0,
                    label="BH/DVCS fitted",
                )
                ax.step(
                    centers,
                    pcomp,
                    where="mid",
                    color=SAMPLE_COLORS["pi0_mc"],
                    linewidth=1.0,
                    label=r"$\pi^0$ fitted",
                )
                ax.step(
                    centers,
                    total,
                    where="mid",
                    color=SAMPLE_COLORS["fit"],
                    linewidth=1.3,
                    label="total fit",
                )

                if axis_name == "ptmiss":
                    if region == "FT":
                        ft_plot_max = min(
                            float(ptmiss_max),
                            max(
                                0.30,
                                1.12 * ehi
                                * math.sin(math.radians(5.0)),
                            ),
                        )
                        ax.set_xlim(0.0, ft_plot_max)
                    else:
                        ax.set_xlim(0.0, ptmiss_max)
                    #endif
                #endif

                qval = float(row.get(qkey, float("nan")))
                cval = float(row.get(ckey, float("nan")))
                ax.text(
                    0.02,
                    0.95,
                    rf"$D/N={qval:.2f}$; "
                    rf"$\chi_P^2/N={cval:.2f}$",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=7.0,
                )
                ax.set_xlabel(xlabel, fontsize=8)
                ax.set_ylabel("entries / bin", fontsize=7.5)
                ax.tick_params(labelsize=7)
                ax.grid(alpha=0.15)
            #endfor

            ax_mx.set_title(
                f"{elo:.2f}-{ehi:.2f} GeV; "
                rf"$f_{{\pi^0}}={fpi0:.3f}$",
                fontsize=9.0,
            )
        #endfor

        # Disable empty tiles.
        for itile in range(nbins, ntile_rows * ncols):
            tile_row = itile // ncols
            tile_col = itile % ncols
            ax = fig.add_subplot(outer[tile_row, tile_col])
            ax.set_axis_off()
        #endfor

        if first_axes is not None:
            handles, labels = first_axes.get_legend_handles_labels()
            if handles:
                fig.legend(
                    handles,
                    labels,
                    loc="upper center",
                    ncol=6,
                    fontsize=7.5,
                    frameon=False,
                    bbox_to_anchor=(0.5, 0.965),
                )
            #endif
        #endif

        fig.suptitle(
            f"{period.label}: ACTUAL nSidis pilot {region} fits at "
            rf"$|M_X^2(ep\gamma_{{tag}})|<{central_support:.2f}$ GeV$^2$"
            "\nEach energy-bin tile: upper = "
            r"$M_X^2(ep)$, lower = $p_{T,\rm miss}$; "
            r"same shared 2D $M_X^2(ep)\otimes p_T$ fit",
            fontsize=12.0,
            y=0.995,
        )
        safe_finalize_figure(
            fig,
            outdir
            / f"canvas_nsidis_pilot_fit_projections_{region.lower()}.png",
            rect=(0, 0.035, 1, 0.94),
        )
        plt.close(fig)
    #endfor



def make_nsidis_photon_acceptance_canvas(
    period: PeriodConfig,
    ns_epg: EPGammaSample,
    ns_epg_f: Dict[str, np.ndarray],
    ns_parent_mask: np.ndarray,
    acceptance: PhotonAngularAcceptance,
    outdir: Path,
) -> None:
    """
    Direct diagnostic for the beam-hole hypothesis.

    Left: actually reconstructed photon theta by detector code, with exact
    observed FT/FD boundaries.
    Right: high-energy predicted-FT pTmiss before versus after requiring the
    predicted probe to lie inside the measured FT angular support.
    """
    det = np.asarray(ns_epg.raw["detector2"], dtype=np.int16)
    theta_reco = photon_theta_deg_from_epgamma(ns_epg)
    parent = np.asarray(ns_parent_mask, dtype=bool)
    good = parent & np.isfinite(theta_reco)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.3))

    bins = np.linspace(0.0, min(45.0, acceptance.fd_theta_max_deg + 3.0), 181)
    for code, label in (
        (PHOTON_DETECTOR_FT, "reconstructed FT (detector=0)"),
        (PHOTON_DETECTOR_FD, "reconstructed FD (detector=1)"),
    ):
        m = good & (det == code)
        axes[0].hist(
            theta_reco[m],
            bins=bins,
            histtype="step",
            density=True,
            linewidth=1.2,
            label=label,
        )
    #endfor
    for x in (
        acceptance.ft_theta_min_deg,
        acceptance.ft_theta_max_deg,
        acceptance.fd_theta_min_deg,
        acceptance.fd_theta_max_deg,
    ):
        axes[0].axvline(x, linestyle="--", linewidth=0.9)
    #endfor
    axes[0].set_yscale("linear")
    axes[0].set_xlabel(r"reconstructed photon $\theta_\gamma$ (deg)")
    axes[0].set_ylabel("density")
    axes[0].set_title(
        "fixed photon angular acceptance\n"
        f"FT=[{acceptance.ft_theta_min_deg:.3f},"
        f"{acceptance.ft_theta_max_deg:.3f}] deg; "
        f"FD=[{acceptance.fd_theta_min_deg:.3f},"
        f"{acceptance.fd_theta_max_deg:.3f}] deg"
    )
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].grid(alpha=0.18)

    # What the old theta<=5.5 definition admitted versus the physical FT range.
    base = (
        parent
        & ns_epg_f["valid_tag"]
        & (ns_epg_f["pred_probe_energy"] >= 2.0)
        & (ns_epg_f["pred_probe_energy"] < 9.5)
        & np.isfinite(ns_epg_f["pred_probe_theta_deg"])
        & np.isfinite(ns_epg_f["stored_pTmiss"])
    )
    legacy_ft = base & (
        ns_epg_f["pred_probe_theta_deg"] <= 5.5
    )
    physical_ft = base & stage2_region_mask(
        ns_epg_f, "FT", 5.5
    )

    ptbins = np.linspace(0.0, 0.90, 121)
    axes[1].hist(
        ns_epg_f["stored_pTmiss"][legacy_ft],
        bins=ptbins,
        histtype="step",
        density=True,
        linewidth=1.2,
        label=r"legacy predicted FT: $\theta^{pred}\leq5.5^\circ$",
    )
    axes[1].hist(
        ns_epg_f["stored_pTmiss"][physical_ft],
        bins=ptbins,
        histtype="step",
        density=True,
        linewidth=1.2,
        label="inside fixed FT theta acceptance",
    )
    axes[1].set_xlabel(r"$p_{T,\rm miss}(ep\gamma)$ (GeV)")
    axes[1].set_ylabel("density")
    axes[1].set_title(
        r"$2<E_{\rm probe}^{pred}<9.5$ GeV: beam-hole population check"
    )
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].grid(alpha=0.18)

    fig.suptitle(
        f"{period.label}: fixed FT/FD photon angular acceptance",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_photon_angular_acceptance.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)



def nsidis_central_tag_mask(
    feat: Dict[str, np.ndarray],
    *,
    ft_theta_max: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    max_probe_energy: float,
) -> np.ndarray:
    """
    Union of the actual FT and FD_all denominator support at the central
    mass-shell choice.  stage2_fit_mask supplies fixed detector acceptance and
    theta_ep > 5 deg; e_p > 2 GeV is added explicitly for nSidis compatibility.
    """
    out = np.zeros(len(feat["pred_probe_energy"]), dtype=bool)
    for region in ("FT", "FD_all"):
        out |= stage2_fit_mask(
            feat,
            region,
            ft_theta_max,
            0.40,
            max_probe_energy,
            mm2_min,
            mm2_max,
            probe_m2_max,
        )
    #endfor

    if "electron_p" not in feat:
        raise KeyError("electron_p missing from denominator feature store.")
    #endif
    out &= (
        np.isfinite(feat["electron_p"])
        & (feat["electron_p"] > 2.0)
    )
    return out


def _fit_row_map_for_support(
    rows: Sequence[Dict[str, object]],
    support: float,
) -> Dict[Tuple[str, float, float], Dict[str, object]]:
    out: Dict[Tuple[str, float, float], Dict[str, object]] = {}
    for row in rows:
        if int(row.get("fit_success", 0)) != 1:
            continue
        #endif
        if abs(
            float(row.get("probe_m2_support_max_GeV2", np.nan))
            - float(support)
        ) > 1.0e-10:
            continue
        #endif
        out[
            (
                str(row["region"]),
                round(float(row["energy_low_GeV"]), 9),
                round(float(row["energy_high_GeV"]), 9),
            )
        ] = row
    #endfor
    return out


def build_nsidis_data_efficiency_rows(
    period: PeriodConfig,
    fit_rows: Sequence[Dict[str, object]],
    selected_feat: Dict[str, np.ndarray],
    association: DataAssociationResult,
    *,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
    source_label: str,
    numerator_purity_rows: Optional[
        Sequence[Dict[str, object]]
    ] = None,
) -> List[Dict[str, object]]:
    """
    Association-first DATA efficiency pilot.

        epsilon_data =
          N(reconstructed companion probe)
          --------------------------------
          f_pi0 * N(selected epgamma tags)

    No MC division is performed here yet.
    """
    fit_map = _fit_row_map_for_support(fit_rows, probe_m2_max)
    purity_map = {
        (
            str(row["region"]),
            round(float(row["energy_low_GeV"]), 9),
            round(float(row["energy_high_GeV"]), 9),
        ): row
        for row in (numerator_purity_rows or [])
    }

    stages = (
        "same_event",
        "positive_remainder",
        "tag_mass_shell",
        "probe_energy",
        "probe_pred_consistent",
    )
    for stage in stages:
        if stage not in association.stage_lookup:
            raise KeyError(
                f"{source_label}: missing association stage '{stage}'."
            )
        #endif
    #endfor

    rows: List[Dict[str, object]] = []
    for region in ("FT", "FD_all"):
        edges = nsidis_energy_edges_for_region(
            region, max_probe_energy
        )
        for ib in range(len(edges) - 1):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])

            m = stage2_fit_mask(
                selected_feat,
                region,
                ft_theta_max,
                elo,
                ehi,
                mm2_min,
                mm2_max,
                probe_m2_max,
            )
            if "electron_p" in selected_feat:
                m &= (
                    np.isfinite(selected_feat["electron_p"])
                    & (selected_feat["electron_p"] > 2.0)
                )
            #endif

            n_tags = int(np.count_nonzero(m))
            counts = {
                stage: int(
                    np.count_nonzero(m & association.stage_lookup[stage])
                )
                for stage in stages
            }

            fit = fit_map.get(
                (region, round(elo, 9), round(ehi, 9))
            )
            fpi0 = (
                float(fit["pi0_fraction"])
                if fit is not None
                else float("nan")
            )
            fpi0_err_low = (
                float(fit.get("pi0_fraction_err_low", np.nan))
                if fit is not None else float("nan")
            )
            fpi0_err_high = (
                float(fit.get("pi0_fraction_err_high", np.nan))
                if fit is not None else float("nan")
            )
            fpi0_model_spread = (
                float(
                    fit.get(
                        "pi0_fraction_model_max_abs_shift",
                        np.nan,
                    )
                )
                if fit is not None
                else float("nan")
            )
            data_den = (
                fpi0 * float(n_tags)
                if np.isfinite(fpi0) and fpi0 > 0.0
                else float("nan")
            )

            n_mass_shell = counts["probe_energy"]
            n_final_raw = counts["probe_pred_consistent"]

            purity_row = purity_map.get(
                (region, round(elo, 9), round(ehi, 9))
            )
            numerator_purity = (
                float(purity_row.get("purity_linear", np.nan))
                if purity_row is not None
                else float("nan")
            )
            numerator_purity_err_low = (
                float(
                    purity_row.get(
                        "purity_linear_err_low", np.nan
                    )
                )
                if purity_row is not None
                else float("nan")
            )
            numerator_purity_err_high = (
                float(
                    purity_row.get(
                        "purity_linear_err_high", np.nan
                    )
                )
                if purity_row is not None
                else float("nan")
            )
            numerator_purity_model_shift = (
                float(
                    purity_row.get(
                        "background_model_abs_shift", np.nan
                    )
                )
                if purity_row is not None
                else float("nan")
            )

            n_final = (
                float(n_final_raw) * numerator_purity
                if np.isfinite(numerator_purity)
                else float("nan")
            )

            # The mass-shell-only diagnostic remains unpurified because the
            # purity fit is defined after the FINAL predicted-probe association.
            eff_mass = (
                float(n_mass_shell) / data_den
                if np.isfinite(data_den) and data_den > 0.0
                else float("nan")
            )
            eff_final = (
                float(n_final) / data_den
                if np.isfinite(data_den) and data_den > 0.0
                else float("nan")
            )

            if (
                np.isfinite(fpi0)
                and fpi0 > 0.0
                and np.isfinite(fpi0_model_spread)
                and fpi0_model_spread >= 0.0
                and n_tags > 0
            ):
                f_model_low = max(
                    1.0e-3,
                    fpi0 - fpi0_model_spread,
                )
                f_model_high = min(
                    1.0 - 1.0e-3,
                    fpi0 + fpi0_model_spread,
                )
                eff_at_f_low = (
                    float(n_final)
                    / (f_model_low * float(n_tags))
                )
                eff_at_f_high = (
                    float(n_final)
                    / (f_model_high * float(n_tags))
                )
                eff_model_min = min(
                    eff_final, eff_at_f_low, eff_at_f_high
                )
                eff_model_max = max(
                    eff_final, eff_at_f_low, eff_at_f_high
                )
                eff_model_sys_low = eff_final - eff_model_min
                eff_model_sys_high = eff_model_max - eff_final
                eff_model_sys = max(
                    eff_model_sys_low,
                    eff_model_sys_high,
                )
            else:
                f_model_low = f_model_high = float("nan")
                eff_model_min = eff_model_max = float("nan")
                eff_model_sys_low = eff_model_sys_high = float("nan")
                eff_model_sys = float("nan")
            #endif

            # Counting fluctuation belongs to the raw associated count, then
            # scales by the fitted purity.
            raw_counting_err = provisional_ratio_error(
                float(n_final_raw), data_den
            )
            counting_err = (
                abs(numerator_purity) * raw_counting_err
                if np.isfinite(numerator_purity)
                and np.isfinite(raw_counting_err)
                else float("nan")
            )
            pi0_eff_err_low = float("nan")
            pi0_eff_err_high = float("nan")
            if (
                np.isfinite(eff_final)
                and np.isfinite(fpi0)
                and fpi0 > 0.0
                and np.isfinite(fpi0_err_low)
                and np.isfinite(fpi0_err_high)
                and n_tags > 0
            ):
                f_hi = fpi0 + fpi0_err_high
                f_lo = max(3.0e-4, fpi0 - fpi0_err_low)
                eff_at_f_hi = float(n_final) / (
                    f_hi * float(n_tags)
                )
                eff_at_f_lo = float(n_final) / (
                    f_lo * float(n_tags)
                )
                pi0_eff_err_low = max(
                    0.0, eff_final - eff_at_f_hi
                )
                pi0_eff_err_high = max(
                    0.0, eff_at_f_lo - eff_final
                )
            #endif

            numerator_purity_eff_err_low = (
                abs(eff_final)
                * numerator_purity_err_low
                / numerator_purity
                if np.isfinite(eff_final)
                and np.isfinite(numerator_purity)
                and numerator_purity > 0.0
                and np.isfinite(numerator_purity_err_low)
                else float("nan")
            )
            numerator_purity_eff_err_high = (
                abs(eff_final)
                * numerator_purity_err_high
                / numerator_purity
                if np.isfinite(eff_final)
                and np.isfinite(numerator_purity)
                and numerator_purity > 0.0
                and np.isfinite(numerator_purity_err_high)
                else float("nan")
            )
            numerator_purity_model_eff_shift = (
                abs(eff_final)
                * numerator_purity_model_shift
                / numerator_purity
                if np.isfinite(eff_final)
                and np.isfinite(numerator_purity)
                and numerator_purity > 0.0
                and np.isfinite(numerator_purity_model_shift)
                else float("nan")
            )

            def _quadrature(*values: float) -> float:
                finite = [
                    float(value)
                    for value in values
                    if np.isfinite(value)
                ]
                return (
                    math.sqrt(sum(value * value for value in finite))
                    if finite else float("nan")
                )
            #enddef

            eff_err_low = _quadrature(
                counting_err,
                pi0_eff_err_low,
                numerator_purity_eff_err_low,
            )
            eff_err_high = _quadrature(
                counting_err,
                pi0_eff_err_high,
                numerator_purity_eff_err_high,
            )

            status = "ok"
            if purity_row is None or not np.isfinite(numerator_purity):
                status = "no_successful_numerator_purity_fit"
            elif fit is None:
                status = "no_successful_denominator_fit"
            elif n_tags <= 0:
                status = "no_selected_tags"
            elif counts["same_event"] <= 0:
                status = "no_same_event_eppi0_candidate"
            elif n_final <= 0:
                status = "no_final_reconstructed_probe"
            elif np.isfinite(eff_final) and eff_final > 1.05:
                status = "efficiency_above_unity_check_numerator_purity"
            #endif

            row = {
                "period": period.key,
                "label": period.label,
                "source": source_label,
                "region": region,
                "energy_low_GeV": elo,
                "energy_high_GeV": ehi,
                "energy_center_GeV": 0.5 * (elo + ehi),
                "mean_probe_energy_GeV": masked_finite_mean(
                    selected_feat["pred_probe_energy"], m
                ),
                "probe_m2_support_max_GeV2": float(probe_m2_max),
                "denominator_fit_success": int(fit is not None),
                "fitted_pi0_fraction": fpi0,
                "fitted_pi0_fraction_err_low": fpi0_err_low,
                "fitted_pi0_fraction_err_high": fpi0_err_high,
                "selected_epgamma_tags": n_tags,
                "fitted_pi0_denominator": data_den,
                "reconstructed_probe_mass_shell": n_mass_shell,
                "reconstructed_probe_final_raw": n_final_raw,
                "numerator_pi0_purity": numerator_purity,
                "numerator_pi0_purity_err_low": (
                    numerator_purity_err_low
                ),
                "numerator_pi0_purity_err_high": (
                    numerator_purity_err_high
                ),
                "numerator_pi0_purity_background_model_abs_shift": (
                    numerator_purity_model_shift
                ),
                "reconstructed_probe_final_purity_corrected": n_final,
                "data_efficiency_mass_shell_only": eff_mass,
                "data_efficiency_final": eff_final,
                "composition_model_fraction_low": f_model_low,
                "composition_model_fraction_high": f_model_high,
                "data_efficiency_composition_model_alt_min": eff_model_min,
                "data_efficiency_composition_model_alt_max": eff_model_max,
                "data_efficiency_composition_model_sys_low": eff_model_sys_low,
                "data_efficiency_composition_model_sys_high": eff_model_sys_high,
                "data_efficiency_composition_model_sys": eff_model_sys,
                "counting_error_final": counting_err,
                "pi0_profile_efficiency_error_low": pi0_eff_err_low,
                "pi0_profile_efficiency_error_high": pi0_eff_err_high,
                "numerator_purity_efficiency_error_low": (
                    numerator_purity_eff_err_low
                ),
                "numerator_purity_efficiency_error_high": (
                    numerator_purity_eff_err_high
                ),
                "numerator_purity_background_model_efficiency_shift": (
                    numerator_purity_model_eff_shift
                ),
                "data_efficiency_error_low": eff_err_low,
                "data_efficiency_error_high": eff_err_high,
                "status": status,
                "interpretation": (
                    "data-only association efficiency corrected by "
                    "final-associated numerator pi0 purity; not yet divided "
                    "by aaogen MC efficiency"
                ),
            }
            previous = n_tags
            for stage in stages:
                count = counts[stage]
                row[f"cutflow_{stage}"] = count
                row[f"cutflow_{stage}_fraction_of_tags"] = (
                    count / n_tags if n_tags else float("nan")
                )
                row[f"cutflow_{stage}_conditional_fraction"] = (
                    count / previous if previous else float("nan")
                )
                previous = count
            #endfor
            rows.append(row)
        #endfor
    #endfor
    return rows


def make_nsidis_data_efficiency_canvas(
    period: PeriodConfig,
    wagon_rows: Sequence[Dict[str, object]],
    nsidis_rows: Sequence[Dict[str, object]],
    outdir: Path,
) -> None:
    """Old-wagon overlap comparison plus nSidis high-energy DATA extension."""
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 8.8))

    for irow, region in enumerate(("FT", "FD_all")):
        ax_eff = axes[irow, 0]
        ax_flow = axes[irow, 1]

        nr = sorted(
            [
                r for r in nsidis_rows
                if r["region"] == region
                and np.isfinite(
                    float(r.get("data_efficiency_final", np.nan))
                )
            ],
            key=row_energy_coordinate,
        )
        wr = sorted(
            [
                r for r in wagon_rows
                if r["region"] == region
                and np.isfinite(
                    float(r.get("data_efficiency_final", np.nan))
                )
            ],
            key=row_energy_coordinate,
        )

        if nr:
            x = np.asarray(
                [row_energy_coordinate(r) for r in nr],
                dtype=float,
            )
            ax_eff.plot(
                x,
                [float(r["data_efficiency_final"]) for r in nr],
                marker="o",
                linewidth=1.25,
                label="nSidis final association",
            )
            ax_eff.plot(
                x,
                [float(r["data_efficiency_mass_shell_only"]) for r in nr],
                marker=".",
                linestyle="--",
                linewidth=0.9,
                label="nSidis mass-shell only",
            )
            for stage, label in (
                ("same_event", "same event"),
                ("probe_energy", "mass-shell + threshold"),
                ("probe_pred_consistent", "final pred-consistent"),
            ):
                ax_flow.plot(
                    x,
                    [
                        float(r[f"cutflow_{stage}_fraction_of_tags"])
                        for r in nr
                    ],
                    marker="o",
                    linewidth=1.0,
                    label=label,
                )
            #endfor
        #endif

        if wr:
            xw = np.asarray(
                [row_energy_coordinate(r) for r in wr],
                dtype=float,
            )
            ax_eff.plot(
                xw,
                [float(r["data_efficiency_final"]) for r in wr],
                marker="s",
                linestyle="--",
                linewidth=1.2,
                label="old wagon reference",
            )
        #endif

        for ax in (ax_eff, ax_flow):
            ax.axvline(2.0, linestyle=":", linewidth=1.0)
            ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
            ax.grid(alpha=0.18)
        #endfor
        ax_eff.axhline(1.0, linestyle=":", linewidth=0.9)
        ax_eff.set_ylabel(r"preliminary $\epsilon_{\rm data}$")
        ax_eff.set_title("FT" if region == "FT" else "FD all")
        ax_flow.set_ylabel("fraction of selected epgamma tags")
        ax_flow.set_ylim(bottom=0.0)
        ax_flow.set_title("association cut flow")
    #endfor

    axes[0, 0].legend(fontsize=8, frameon=False)
    axes[0, 1].legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"{period.label}: association-first reconstructed-probe pilot "
        "(DATA ONLY; not yet data/MC)",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_data_efficiency_pilot.png",
        rect=(0, 0, 1, 0.95),
    )
    plt.close(fig)


def make_nsidis_associated_eppi0_exclusivity_canvas(
    period: PeriodConfig,
    eppi0_features: Dict[str, np.ndarray],
    association: DataAssociationResult,
    outdir: Path,
) -> Dict[str, object]:
    """
    Inspect eppi0 exclusivity only AFTER exact tag/probe association.

    This determines whether association itself supplies sufficient numerator
    purity before we impose any additional reconstructed-side exclusivity cut.
    """
    best = association.best_match
    if len(best.eppi0_index) == 0:
        return {
            "best_associations": 0,
            "best_pred_consistent_associations": 0,
        }
    #endif

    passes = np.asarray(
        association.best_diagnostics.get(
            "passes_pred_consistency",
            np.zeros(len(best.eppi0_index), dtype=bool),
        ),
        dtype=bool,
    )
    best_idx = np.asarray(best.eppi0_index, dtype=np.int64)
    final_idx = best_idx[passes]
    best_unique = np.unique(best_idx)
    final_unique = np.unique(final_idx)

    specs = (
        ("Mx2", np.linspace(-0.08, 0.08, 161),
         r"$M_X^2(ep\pi^0)$ (GeV$^2$)"),
        ("pTmiss", np.linspace(0.0, 0.80, 161),
         r"$p_{T,\rm miss}(ep\pi^0)$ (GeV)"),
        ("pmiss", np.linspace(0.0, 1.5, 161),
         r"$|\vec p_{\rm miss}(ep\pi^0)|$ (GeV)"),
        ("Emiss", np.linspace(-1.5, 2.5, 181),
         r"$E_{\rm miss}(ep\pi^0)$ (GeV)"),
    )

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0))
    n_all = len(eppi0_features["Mx2"])
    idx_all = np.arange(n_all)
    best_mask = np.isin(idx_all, best_unique)
    final_mask = np.isin(idx_all, final_unique)

    for ax, (key, bins, xlabel) in zip(axes.ravel(), specs):
        h_all, edges = _normalized_hist(eppi0_features[key], bins)
        h_best, _ = _normalized_hist(
            eppi0_features[key], bins, best_mask
        )
        h_final, _ = _normalized_hist(
            eppi0_features[key], bins, final_mask
        )
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax.step(
            centers, h_all, where="mid", linewidth=0.9,
            label="all nSidis eppi0X",
        )
        ax.step(
            centers, h_best, where="mid", linewidth=1.1,
            label="best mass-shell association",
        )
        ax.step(
            centers, h_final, where="mid", linewidth=1.3,
            label="best + pred-consistent",
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel("unit-normalized entries")
        ax.grid(alpha=0.18)
    #endfor

    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle(
        f"{period.label}: eppi0X exclusivity AFTER tag/probe association",
        fontsize=13,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_associated_eppi0_exclusivity.png",
        rect=(0, 0, 1, 0.95),
    )
    plt.close(fig)

    stress_rows: List[Dict[str, object]] = []
    if final_unique.size:
        mx = np.asarray(eppi0_features["Mx2"][final_unique], dtype=float)
        pt = np.asarray(
            eppi0_features["pTmiss"][final_unique], dtype=float
        )
        for mxmax in (0.02, 0.04, 0.06, 0.10):
            for ptmax in (0.20, 0.30, 0.50):
                good = (
                    np.isfinite(mx)
                    & np.isfinite(pt)
                    & (np.abs(mx) < mxmax)
                    & (pt < ptmax)
                )
                stress_rows.append({
                    "period": period.key,
                    "label": period.label,
                    "abs_Mx2_max_GeV2": mxmax,
                    "pTmiss_max_GeV": ptmax,
                    "unique_final_associated_pi0": int(final_unique.size),
                    "surviving_count": int(np.count_nonzero(good)),
                    "surviving_fraction": float(np.mean(good)),
                    "status": "diagnostic_only_no_cut_applied",
                })
            #endfor
        #endfor
    #endif

    if stress_rows:
        write_rows_csv(
            stress_rows,
            outdir / "nsidis_associated_eppi0_exclusivity_stress.csv",
        )
    #endif

    return {
        "best_associations": int(len(best_idx)),
        "best_pred_consistent_associations": int(np.count_nonzero(passes)),
        "unique_best_eppi0_candidates": int(best_unique.size),
        "unique_final_eppi0_candidates": int(final_unique.size),
        "stress_rows": int(len(stress_rows)),
    }



def build_mc_probe_stage_lookup(
    n_epgamma: int,
    pair_np: Dict[str, np.ndarray],
    *,
    mgg_min: float,
    mgg_max: float,
    remainder_mass2_max: float,
    reco_probe_energy_min: float,
    probe_angle_max_deg: float,
    probe_frac_energy_max: float,
) -> Dict[str, np.ndarray]:
    """
    Convert the internal aaogen tag/probe association into one boolean lookup
    per original aaogen epgamma entry.

    The DATA association stages and MC stages are deliberately aligned:

      same_event              -> successful e/p parent match in MC
      positive_remainder      -> physical P_pi0 - k_tag photon remainder
      tag_mass_shell          -> pi0 mass window + photon-like remainder
      probe_energy            -> reconstructed companion above threshold
      probe_pred_consistent   -> loose predicted/reconstructed agreement

    The pi0 mass window is explicit in MC because the reconstructed MC eppi0
    source need not be assumed to have exactly the same upstream mass window
    as the nSidis data tree.
    """
    stages = (
        "same_event",
        "positive_remainder",
        "tag_mass_shell",
        "probe_energy",
        "probe_pred_consistent",
    )
    lookup = {
        stage: np.zeros(int(n_epgamma), dtype=bool)
        for stage in stages
    }
    if not pair_np:
        return lookup
    #endif

    idx = np.asarray(pair_np["epg_index"], dtype=np.int64)
    valid_index = (idx >= 0) & (idx < int(n_epgamma))
    if not np.any(valid_index):
        return lookup
    #endif

    idx = idx[valid_index]

    def arr(name: str) -> np.ndarray:
        return np.asarray(pair_np[name])[valid_index]
    #enddef

    positive = (
        np.isfinite(arr("reco_probe_energy"))
        & np.isfinite(arr("reco_probe_p"))
        & (arr("reco_probe_energy") > 0.0)
        & (arr("reco_probe_p") > 0.0)
    )
    mass_shell = (
        positive
        & np.isfinite(arr("pi0_mass"))
        & (arr("pi0_mass") >= float(mgg_min))
        & (arr("pi0_mass") <= float(mgg_max))
        & np.isfinite(arr("reco_probe_mass2"))
        & (
            np.abs(arr("reco_probe_mass2"))
            <= float(remainder_mass2_max)
        )
    )
    above_threshold = (
        mass_shell
        & (
            arr("reco_probe_energy")
            >= float(reco_probe_energy_min)
        )
    )
    pred_consistent = (
        above_threshold
        & np.isfinite(arr("pred_probe_energy"))
        & (arr("pred_probe_energy") > 0.0)
        & np.isfinite(arr("probe_opening_residual_deg"))
        & (
            arr("probe_opening_residual_deg")
            <= float(probe_angle_max_deg)
        )
        & np.isfinite(arr("probe_delta_E_over_E"))
        & (
            np.abs(arr("probe_delta_E_over_E"))
            <= float(probe_frac_energy_max)
        )
    )

    lookup["same_event"][np.unique(idx)] = True
    for name, pair_mask in (
        ("positive_remainder", positive),
        ("tag_mass_shell", mass_shell),
        ("probe_energy", above_threshold),
        ("probe_pred_consistent", pred_consistent),
    ):
        lookup[name][np.unique(idx[pair_mask])] = True
    #endfor

    return lookup


def binomial_error_from_counts(
    passed: int,
    total: int,
) -> float:
    """Simple binomial counting uncertainty for an MC reconstruction efficiency."""
    if total <= 0:
        return float("nan")
    #endif
    p = float(passed) / float(total)
    if not np.isfinite(p):
        return float("nan")
    #endif
    return math.sqrt(max(0.0, p * (1.0 - p) / float(total)))


def build_aaogen_mc_efficiency_rows(
    period: PeriodConfig,
    pi0_f: Dict[str, np.ndarray],
    mc_stage_lookup: Dict[str, np.ndarray],
    *,
    ft_theta_max: float,
    max_probe_energy: float,
    mm2_min: float,
    mm2_max: float,
    probe_m2_max: float,
) -> List[Dict[str, object]]:
    """
    Measure aaogen reconstructed-companion efficiency with the exact same
    denominator support used for the data composition fit.

    Unlike data, no fitted f_pi0 factor is needed: every selected aaogen epgamma
    tag is a true exclusive-pi0 denominator tag.
    """
    stages = (
        "same_event",
        "positive_remainder",
        "tag_mass_shell",
        "probe_energy",
        "probe_pred_consistent",
    )
    for stage in stages:
        if stage not in mc_stage_lookup:
            raise KeyError(f"MC stage lookup is missing '{stage}'.")
        #endif
    #endfor

    rows: List[Dict[str, object]] = []
    for region in ("FT", "FD_all"):
        edges = nsidis_energy_edges_for_region(
            region, max_probe_energy
        )
        for ib in range(len(edges) - 1):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])

            mask = stage2_fit_mask(
                pi0_f,
                region,
                ft_theta_max,
                elo,
                ehi,
                mm2_min,
                mm2_max,
                probe_m2_max,
            )
            # Match the nSidis parent support exactly.
            if "electron_p" in pi0_f:
                mask &= (
                    np.isfinite(pi0_f["electron_p"])
                    & (pi0_f["electron_p"] > 2.0)
                )
            #endif

            n_tags = int(np.count_nonzero(mask))
            counts = {
                stage: int(
                    np.count_nonzero(mask & mc_stage_lookup[stage])
                )
                for stage in stages
            }
            n_mass = counts["probe_energy"]
            n_final = counts["probe_pred_consistent"]

            eff_mass = (
                n_mass / n_tags if n_tags else float("nan")
            )
            eff_final = (
                n_final / n_tags if n_tags else float("nan")
            )

            rows.append({
                "period": period.key,
                "label": period.label,
                "source": "aaogen",
                "region": region,
                "energy_low_GeV": elo,
                "energy_high_GeV": ehi,
                "energy_center_GeV": 0.5 * (elo + ehi),
                "mean_probe_energy_GeV": masked_finite_mean(
                    pi0_f["pred_probe_energy"], mask
                ),
                "probe_m2_support_max_GeV2": float(probe_m2_max),
                "selected_true_pi0_tags": n_tags,
                "reconstructed_probe_mass_shell": n_mass,
                "reconstructed_probe_final": n_final,
                "mc_efficiency_mass_shell_only": eff_mass,
                "mc_efficiency_final": eff_final,
                "mc_efficiency_final_stat_error": (
                    binomial_error_from_counts(n_final, n_tags)
                ),
                "status": (
                    "ok" if n_tags > 0
                    else "no_selected_true_pi0_tags"
                ),
                **{
                    f"cutflow_{stage}": counts[stage]
                    for stage in stages
                },
            })
        #endfor
    #endfor
    return rows


def combine_data_mc_scale_factor_rows(
    data_rows: Sequence[Dict[str, object]],
    mc_rows: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    """
    Form SF_gamma = epsilon_data / epsilon_MC for matching period/region/E bins.

    The data error remains provisional because the fitted-pi0 denominator
    uncertainty/covariance is not yet propagated.  MC binomial counting error is
    propagated together with the provisional data counting error.
    """
    mc_map = {
        (
            str(r["period"]),
            str(r["region"]),
            round(float(r["energy_low_GeV"]), 9),
            round(float(r["energy_high_GeV"]), 9),
        ): r
        for r in mc_rows
    }

    out: List[Dict[str, object]] = []
    for d in data_rows:
        key = (
            str(d["period"]),
            str(d["region"]),
            round(float(d["energy_low_GeV"]), 9),
            round(float(d["energy_high_GeV"]), 9),
        )
        m = mc_map.get(key)

        ed = float(d.get("data_efficiency_final", np.nan))
        sed_low = float(
            d.get("data_efficiency_error_low", np.nan)
        )
        sed_high = float(
            d.get("data_efficiency_error_high", np.nan)
        )
        em = (
            float(m.get("mc_efficiency_final", np.nan))
            if m is not None else float("nan")
        )
        sem = (
            float(m.get("mc_efficiency_final_stat_error", np.nan))
            if m is not None else float("nan")
        )

        sf = (
            ed / em
            if np.isfinite(ed) and np.isfinite(em) and em > 0.0
            else float("nan")
        )

        def propagated_sf_error(
            data_error: float,
            mc_error: float,
        ) -> float:
            rel2 = 0.0
            nrel = 0
            if (
                np.isfinite(data_error)
                and np.isfinite(ed)
                and ed > 0.0
            ):
                rel2 += (data_error / ed) ** 2
                nrel += 1
            #endif
            if (
                np.isfinite(mc_error)
                and np.isfinite(em)
                and em > 0.0
            ):
                rel2 += (mc_error / em) ** 2
                nrel += 1
            #endif
            return (
                abs(sf) * math.sqrt(rel2)
                if np.isfinite(sf) and nrel > 0
                else float("nan")
            )
        #enddef

        sf_err_low = propagated_sf_error(sed_low, sem)
        sf_err_high = propagated_sf_error(sed_high, sem)

        data_comp_low = float(
            d.get("data_efficiency_composition_model_sys_low", np.nan)
        )
        data_comp_high = float(
            d.get("data_efficiency_composition_model_sys_high", np.nan)
        )
        sf_comp_low = (
            data_comp_low / em
            if np.isfinite(data_comp_low) and np.isfinite(em) and em > 0.0
            else float("nan")
        )
        sf_comp_high = (
            data_comp_high / em
            if np.isfinite(data_comp_high) and np.isfinite(em) and em > 0.0
            else float("nan")
        )
        sf_comp = (
            max(sf_comp_low, sf_comp_high)
            if np.isfinite(sf_comp_low) and np.isfinite(sf_comp_high)
            else float("nan")
        )

        numerator_purity_sf_model_shift = (
            abs(sf)
            * float(
                d.get(
                    "numerator_pi0_purity_background_model_abs_shift",
                    np.nan,
                )
            )
            / float(d.get("numerator_pi0_purity", np.nan))
            if np.isfinite(sf)
            and np.isfinite(
                float(d.get("numerator_pi0_purity", np.nan))
            )
            and float(d.get("numerator_pi0_purity", np.nan)) > 0.0
            and np.isfinite(
                float(
                    d.get(
                        "numerator_pi0_purity_background_model_abs_shift",
                        np.nan,
                    )
                )
            )
            else float("nan")
        )


        out.append({
            "period": d["period"],
            "label": d["label"],
            "source": d["source"],
            "region": d["region"],
            "energy_low_GeV": d["energy_low_GeV"],
            "energy_high_GeV": d["energy_high_GeV"],
            "energy_center_GeV": d["energy_center_GeV"],
            "mean_probe_energy_GeV": d.get(
                "mean_probe_energy_GeV", d["energy_center_GeV"]
            ),
            "mc_mean_probe_energy_GeV": (
                m.get("mean_probe_energy_GeV", np.nan)
                if m is not None else float("nan")
            ),
            "probe_m2_support_max_GeV2": d[
                "probe_m2_support_max_GeV2"
            ],
            "fitted_pi0_fraction_data": d["fitted_pi0_fraction"],
            "numerator_pi0_purity": d.get(
                "numerator_pi0_purity", np.nan
            ),
            "numerator_pi0_purity_err_low": d.get(
                "numerator_pi0_purity_err_low", np.nan
            ),
            "numerator_pi0_purity_err_high": d.get(
                "numerator_pi0_purity_err_high", np.nan
            ),
            "numerator_pi0_purity_background_model_abs_shift": d.get(
                "numerator_pi0_purity_background_model_abs_shift", np.nan
            ),
            "data_efficiency_final": ed,
            "data_efficiency_error_low": sed_low,
            "data_efficiency_error_high": sed_high,
            "mc_efficiency_final": em,
            "mc_efficiency_stat_error": sem,
            "photon_efficiency_scale_factor": sf,
            "scale_factor_stat_error_low": sf_err_low,
            "scale_factor_stat_error_high": sf_err_high,
            "scale_factor_composition_model_sys_low": sf_comp_low,
            "scale_factor_composition_model_sys_high": sf_comp_high,
            "scale_factor_composition_model_sys": sf_comp,
            "scale_factor_numerator_purity_background_model_abs_shift": (
                numerator_purity_sf_model_shift
            ),
            "data_status": d["status"],
            "mc_status": (
                m["status"] if m is not None else "missing_mc_row"
            ),
            "status": (
                "PRELIMINARY_stat_error_includes_denominator_and_numerator_pi0_fits"
                if np.isfinite(sf)
                else "invalid_scale_factor"
            ),
        })
    #endfor
    return out


def make_nsidis_scale_factor_canvas(
    period: PeriodConfig,
    sf_rows: Sequence[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Per-period efficiency and data/MC scale-factor summary.

    Neighboring energy bins are connected by lines. Statistical uncertainties
    are thin error bars. Composition-model systematics are shaded rectangles
    spanning the full energy-bin width.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.4))

    for irow, region in enumerate(("FT", "FD_all")):
        rr = sorted(
            [
                row for row in sf_rows
                if row["region"] == region
                and np.isfinite(float(row.get("photon_efficiency_scale_factor", np.nan)))
            ],
            key=row_energy_coordinate,
        )
        if not rr:
            continue
        #endif

        x = np.asarray([row_energy_coordinate(row) for row in rr], dtype=float)
        elo = np.asarray([float(row["energy_low_GeV"]) for row in rr], dtype=float)
        ehi = np.asarray([float(row["energy_high_GeV"]) for row in rr], dtype=float)
        ed = np.asarray([float(row["data_efficiency_final"]) for row in rr], dtype=float)
        em = np.asarray([float(row["mc_efficiency_final"]) for row in rr], dtype=float)
        sf = np.asarray([float(row["photon_efficiency_scale_factor"]) for row in rr], dtype=float)

        edlo = np.asarray([float(row.get("data_efficiency_error_low", np.nan)) for row in rr])
        edhi = np.asarray([float(row.get("data_efficiency_error_high", np.nan)) for row in rr])
        emerr = np.asarray([float(row.get("mc_efficiency_stat_error", np.nan)) for row in rr])
        sflo = np.asarray([float(row.get("scale_factor_stat_error_low", np.nan)) for row in rr])
        sfhi = np.asarray([float(row.get("scale_factor_stat_error_high", np.nan)) for row in rr])
        syslo = np.asarray([float(row.get("scale_factor_composition_model_sys_low", np.nan)) for row in rr])
        syshi = np.asarray([float(row.get("scale_factor_composition_model_sys_high", np.nan)) for row in rr])

        edsyslo = syslo * em
        edsyshi = syshi * em

        ax_e = axes[irow, 0]
        ax_s = axes[irow, 1]

        for j in range(len(rr)):
            if np.isfinite(ed[j]) and np.isfinite(edsyslo[j]) and np.isfinite(edsyshi[j]):
                ax_e.fill_between(
                    [elo[j], ehi[j]],
                    [ed[j]-edsyslo[j]]*2,
                    [ed[j]+edsyshi[j]]*2,
                    color="tab:orange", alpha=0.16, linewidth=0.0,
                    label="composition-model systematic" if j == 0 else None,
                )
            #endif
            if np.isfinite(sf[j]) and np.isfinite(syslo[j]) and np.isfinite(syshi[j]):
                ax_s.fill_between(
                    [elo[j], ehi[j]],
                    [sf[j]-syslo[j]]*2,
                    [sf[j]+syshi[j]]*2,
                    color="tab:orange", alpha=0.16, linewidth=0.0,
                    label="composition-model systematic" if j == 0 else None,
                )
            #endif
        #endfor

        ax_e.errorbar(
            x, ed, yerr=np.vstack((edlo, edhi)),
            marker="o", linestyle="-", linewidth=1.2, capsize=2,
            color="black", label="data", zorder=5,
        )
        ax_e.errorbar(
            x, em, yerr=emerr,
            marker="s", linestyle="-", linewidth=1.2, capsize=2,
            label="aaogen MC", zorder=5,
        )
        ax_s.errorbar(
            x, sf, yerr=np.vstack((sflo, sfhi)),
            marker="o", linestyle="-", linewidth=1.2, capsize=2,
            color="black", label="statistical", zorder=5,
        )

        for ax in (ax_e, ax_s):
            ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
            ax.grid(alpha=0.18)
        #endfor

        ax_e.set_ylabel("reconstructed-probe efficiency")
        ax_e.set_title("FT" if region == "FT" else "FD all")
        ax_e.legend(fontsize=7.5, frameon=False)

        ax_s.axhline(1.0, linestyle="--", linewidth=1.0)
        ax_s.set_ylabel(
            r"$SF_{\gamma}=\epsilon_{\mathrm{data}}/\epsilon_{\mathrm{MC}}$"
        )
        ax_s.set_title("data/MC scale factor")
        ax_s.legend(fontsize=7.5, frameon=False)
    #endfor

    fig.suptitle(
        f"{period.label}: photon efficiency and data/MC scale factor\\n"
        "thin bars = statistical; orange boxes = denominator-composition systematic",
        fontsize=12.5, y=0.985,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_data_mc_scale_factor.png",
        rect=(0.0, 0.0, 1.0, 0.92),
    )
    plt.close(fig)


def read_csv_rows_simple(path: Path) -> List[Dict[str, str]]:
    """Read a CSV into string-valued dictionaries for aggregate plotting."""
    if not path.exists():
        return []
    #endif
    with path.open("r", newline="") as f:
        return list(csv.DictReader(f))
    #endwith


def make_nsidis_period_comparison(
    nsroot: Path,
    periods: Sequence[PeriodConfig],
) -> None:
    """
    Final all-period scale-factor comparison: FT and FD-all, full energy range.
    No composition-systematic boxes are drawn because all five periods are
    overlaid on the same canvas.
    """
    period_rows = {}
    period_labels = {}

    for period in periods:
        path = (
            nsroot / "stage3_efficiency" / period.key
            / "nsidis_photon_efficiency_scale_factors.csv"
        )
        rows = read_csv_rows_simple(path)
        if rows:
            period_rows[period.key] = rows
            period_labels[period.key] = period.label
        #endif
    #endfor

    if len(period_rows) < 2:
        return
    #endif

    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.6), squeeze=False)
    axes = axes[0]
    combined_rows = []

    for ax, region in zip(axes, ("FT", "FD_all")):
        for period_key, rows in period_rows.items():
            rr = []
            for row in rows:
                if row.get("region") != region:
                    continue
                #endif
                try:
                    x = row_energy_coordinate(row)
                    y = float(row["photon_efficiency_scale_factor"])
                    el = float(row.get("scale_factor_stat_error_low", "nan"))
                    eh = float(row.get("scale_factor_stat_error_high", "nan"))
                except (KeyError, TypeError, ValueError):
                    continue
                #endtry
                if np.isfinite(x) and np.isfinite(y):
                    rr.append((x, y, el, eh, row))
                #endif
            #endfor

            rr.sort(key=lambda item: item[0])
            if not rr:
                continue
            #endif

            x = np.asarray([item[0] for item in rr], dtype=float)
            y = np.asarray([item[1] for item in rr], dtype=float)
            el = np.asarray([item[2] for item in rr], dtype=float)
            eh = np.asarray([item[3] for item in rr], dtype=float)

            ax.errorbar(
                x, y, yerr=np.vstack((el, eh)),
                marker="o", linestyle="-", linewidth=1.0, capsize=2,
                label=period_labels[period_key],
            )

            for _, _, _, _, raw in rr:
                combined_rows.append(
                    {"period": period_key, "label": period_labels[period_key], **raw}
                )
            #endfor
        #endfor

        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        ax.set_ylabel(
            r"$SF_{\gamma}=\epsilon_{\mathrm{data}}/\epsilon_{\mathrm{MC}}$"
        )
        ax.set_title("FT" if region == "FT" else "FD all")
        ax.grid(alpha=0.18)
        ax.legend(fontsize=8, frameon=False)
    #endfor

    fig.suptitle("Photon-efficiency period comparison", fontsize=13, y=0.985)
    safe_finalize_figure(
        fig,
        nsroot / "canvas_nsidis_period_scale_factor_comparison.png",
        rect=(0.0, 0.0, 1.0, 0.93),
    )
    plt.close(fig)

    write_rows_csv(
        combined_rows,
        nsroot / "nsidis_period_scale_factor_comparison.csv",
    )

    make_nsidis_period_efficiency_summary(nsroot, periods, period_rows)


def make_nsidis_period_efficiency_summary(
    nsroot: Path,
    periods: Sequence[PeriodConfig],
    period_rows: Dict[str, List[Dict[str, str]]],
) -> None:
    """
    Parent-level 3x2 summary of the DATA/MC EFFICIENCY RATIO only.

    One subplot per run period; the sixth panel is intentionally blank.
    Each period shows the FT and FD-all photon-efficiency scale factors,

        SF_gamma = epsilon_data / epsilon_MC,

    versus predicted probe-photon energy.

    Because this is another multi-period summary canvas, composition-model
    systematic boxes are intentionally omitted.  Only statistical error bars
    are shown, with lines connecting neighboring energy bins.
    """
    fig, axes = plt.subplots(
        3,
        2,
        figsize=(13.8, 13.0),
        squeeze=False,
    )
    flat = axes.ravel()

    for ip, period in enumerate(periods):
        ax = flat[ip]
        rows = period_rows.get(
            period.key,
            [],
        )
        if not rows:
            ax.set_axis_off()
            continue
        #endif

        for region, label, marker in (
            ("FT", "FT", "o"),
            ("FD_all", "FD", "s"),
        ):
            rr = sorted(
                [
                    row
                    for row in rows
                    if row.get("region") == region
                ],
                key=row_energy_coordinate,
            )
            if not rr:
                continue
            #endif

            x = np.asarray(
                [
                    row_energy_coordinate(row)
                    for row in rr
                ],
                dtype=float,
            )
            sf = np.asarray(
                [
                    float(
                        row[
                            "photon_efficiency_scale_factor"
                        ]
                    )
                    for row in rr
                ],
                dtype=float,
            )
            sf_low = np.asarray(
                [
                    float(
                        row.get(
                            "scale_factor_stat_error_low",
                            np.nan,
                        )
                    )
                    for row in rr
                ],
                dtype=float,
            )
            sf_high = np.asarray(
                [
                    float(
                        row.get(
                            "scale_factor_stat_error_high",
                            np.nan,
                        )
                    )
                    for row in rr
                ],
                dtype=float,
            )

            finite = np.isfinite(x) & np.isfinite(sf)
            if not np.any(finite):
                continue
            #endif

            elo = np.asarray(
                [
                    float(row["energy_low_GeV"])
                    for row in rr
                ],
                dtype=float,
            )
            ehi = np.asarray(
                [
                    float(row["energy_high_GeV"])
                    for row in rr
                ],
                dtype=float,
            )
            sys_low = np.asarray(
                [
                    float(
                        row.get(
                            "scale_factor_composition_model_sys_low",
                            np.nan,
                        )
                    )
                    for row in rr
                ],
                dtype=float,
            )
            sys_high = np.asarray(
                [
                    float(
                        row.get(
                            "scale_factor_composition_model_sys_high",
                            np.nan,
                        )
                    )
                    for row in rr
                ],
                dtype=float,
            )

            # Composition-model systematic as a shaded box spanning the
            # full energy-bin width.  This summary has one period per panel,
            # so FT/FD systematic boxes remain readable.
            first_sys_label = True
            for j in range(len(rr)):
                if (
                    finite[j]
                    and np.isfinite(sys_low[j])
                    and np.isfinite(sys_high[j])
                ):
                    ax.fill_between(
                        [elo[j], ehi[j]],
                        [
                            sf[j] - sys_low[j],
                            sf[j] - sys_low[j],
                        ],
                        [
                            sf[j] + sys_high[j],
                            sf[j] + sys_high[j],
                        ],
                        alpha=0.14,
                        linewidth=0.0,
                        label=(
                            "composition-model systematic"
                            if first_sys_label
                            else None
                        ),
                    )
                    first_sys_label = False
                #endif
            #endfor

            ax.errorbar(
                x[finite],
                sf[finite],
                yerr=np.vstack(
                    (
                        sf_low[finite],
                        sf_high[finite],
                    )
                ),
                marker=marker,
                linestyle="-",
                linewidth=1.2,
                capsize=2,
                label=label,
                zorder=5,
            )
        #endfor

        ax.axhline(
            1.0,
            linestyle="--",
            linewidth=1.0,
        )
        ax.set_title(period.label)
        ax.set_xlabel(
            r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
        )
        ax.set_ylabel(
            r"$SF_{\gamma}="
            r"\epsilon_{\mathrm{data}}/"
            r"\epsilon_{\mathrm{MC}}$"
        )
        ax.grid(alpha=0.18)
        ax.legend(
            fontsize=8,
            frameon=False,
        )
    #endfor

    for ip in range(
        len(periods),
        len(flat),
    ):
        flat[ip].set_axis_off()
    #endfor

    fig.suptitle(
        "Photon-efficiency data/MC ratio by run period",
        fontsize=13,
        y=0.985,
    )
    safe_finalize_figure(
        fig,
        nsroot
        / "canvas_nsidis_period_scale_factor_by_period.png",
        rect=(0.0, 0.0, 1.0, 0.96),
    )
    plt.close(fig)


def direct_args_dict_keys(function) -> set[str]:
    """Return literal keys used as args_dict['key'] inside a function."""
    source = textwrap.dedent(inspect.getsource(function))
    tree = ast.parse(source)
    keys: set[str] = set()

    for node in ast.walk(tree):
        if not (
            isinstance(node, ast.Subscript)
            and isinstance(node.value, ast.Name)
            and node.value.id == "args_dict"
        ):
            continue
        #endif
        if (
            isinstance(node.slice, ast.Constant)
            and isinstance(node.slice.value, str)
        ):
            keys.add(node.slice.value)
        #endif
    #endfor
    return keys


def module_argparse_destinations() -> set[str]:
    """
    Extract every argparse destination from the entire current source file.

    The parser is assembled across several helper functions, so scanning only
    main() is insufficient.
    """
    source = Path(__file__).read_text()
    tree = ast.parse(source)
    dests: set[str] = set()

    for node in ast.walk(tree):
        if not (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr == "add_argument"
        ):
            continue
        #endif

        options: List[str] = []
        explicit_dest = None

        for arg in node.args:
            if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                options.append(arg.value)
            #endif
        #endfor

        for kw in node.keywords:
            if (
                kw.arg == "dest"
                and isinstance(kw.value, ast.Constant)
                and isinstance(kw.value.value, str)
            ):
                explicit_dest = kw.value.value
            #endif
        #endfor

        if explicit_dest is not None:
            dests.add(explicit_dest)
            continue
        #endif

        long_options = [
            option for option in options if option.startswith("--")
        ]
        if long_options:
            dests.add(
                long_options[0][2:].replace("-", "_")
            )
        #endif
    #endfor

    return dests


def validate_nsidis_argument_contract(
    args_dict: Dict[str, object],
) -> None:
    """
    Fail before expensive I/O if the nSidis worker contains a stale/misspelled
    direct args_dict key.

    This is intentionally executed before preflight_nsidis_study().
    """
    required = direct_args_dict_keys(process_nsidis_study_period)
    parser_dests = module_argparse_destinations()

    missing_from_parser = sorted(required.difference(parser_dests))
    if missing_from_parser:
        raise RuntimeError(
            "nSidis worker/argparse contract failure BEFORE processing: "
            "worker references destination(s) that do not exist in the parser: "
            + ", ".join(missing_from_parser)
        )
    #endif

    missing_from_runtime = sorted(required.difference(args_dict))
    if missing_from_runtime:
        raise RuntimeError(
            "nSidis runtime argument-contract failure BEFORE processing: "
            "vars(args) is missing required worker key(s): "
            + ", ".join(missing_from_runtime)
        )
    #endif




# -----------------------------------------------------------------------------
# Dedicated nSidis denominator-driver study
# -----------------------------------------------------------------------------

NSIDIS_DRIVER_STUDY_MODELS: Tuple[
    Tuple[str, str, str], ...
] = (
    (
        "ptmiss_only",
        r"$p_{T,\mathrm{miss}}$ only",
        "pTmiss",
    ),
    (
        "theta_gg_only",
        r"$\theta_{\gamma\gamma}$ only",
        "theta_gg",
    ),
    (
        "delta_phi_only",
        r"$|\Delta\phi|$ only",
        "delta_phi",
    ),
    (
        "theta_gg_x_ptmiss",
        r"$\theta_{\gamma\gamma}\otimes p_{T,\mathrm{miss}}$",
        "theta_gg_x_pTmiss",
    ),
    (
        "delta_phi_x_ptmiss",
        r"$|\Delta\phi|\otimes p_{T,\mathrm{miss}}$",
        "delta_phi_x_pTmiss",
    ),
    (
        "mx2_ep_x_ptmiss",
        r"$M_X^2(ep)\otimes p_{T,\mathrm{miss}}$",
        "mx2_ep_x_pTmiss",
    ),
    (
        "mx2_ep_x_theta_gg",
        r"$M_X^2(ep)\otimes\theta_{\gamma\gamma}$",
        "mx2_ep_x_theta_gg",
    ),
)


def nsidis_driver_study_common_mask(
    feat: Dict[str, np.ndarray],
    region_name: str,
    ft_theta_max: float,
    elo: float,
    ehi: float,
    probe_m2_max: float,
    mm2_min: float,
    mm2_max: float,
    ptmiss_max: float,
    theta_max: float,
    parent_mask: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Common-population selection used by every candidate discriminator.

    This is intentionally stricter than "candidate-specific available data":
    every method is evaluated on the SAME events, so differences in fitted
    f_pi0 genuinely reflect the driver/model and not a different event sample.

    Emiss2 is deliberately not included as a candidate.  For epgamma,
    the missing energy after the measured e,p,gamma is exactly the same
    kinematic quantity used here as E_probe^pred, i.e. the binning variable.
    Fitting it inside an E_probe^pred bin would therefore be partly circular
    and unusually sensitive to generator spectral shapes.
    """
    n = len(feat["pred_probe_energy"])
    theta_ep = feat.get("theta_ep_deg")
    electron_p = feat.get("electron_p")
    pt = feat.get("stored_pTmiss")
    th = feat.get("stored_theta_gamma_gamma")
    dphi = feat.get("stored_delta_phi_abs_rad")

    mask = (
        np.asarray(feat["valid_tag"], dtype=bool)
        & stage2_region_mask(feat, region_name, ft_theta_max)
        & np.isfinite(feat["pred_probe_energy"])
        & (feat["pred_probe_energy"] >= elo)
        & (feat["pred_probe_energy"] < ehi)
        & np.isfinite(feat["pred_probe_mass2"])
        & np.isfinite(feat["ep_missing_mass2"])
        & (feat["ep_missing_mass2"] >= mm2_min)
        & (feat["ep_missing_mass2"] < mm2_max)
    )

    if theta_ep is not None:
        mask &= (
            np.isfinite(theta_ep)
            & (theta_ep > THETA_EP_MIN_DEG)
        )
    #endif

    if electron_p is not None:
        mask &= (
            np.isfinite(electron_p)
            & (electron_p > 2.0)
        )
    #endif

    if pt is None or th is None or dphi is None:
        return np.zeros(n, dtype=bool)
    #endif

    mask &= (
        np.isfinite(pt)
        & (pt >= 0.0)
        & (pt < ptmiss_max)
        & np.isfinite(th)
        & (th >= 0.0)
        & (th < theta_max)
        & np.isfinite(dphi)
        & (dphi >= 0.0)
        & (dphi <= math.pi)
    )

    if parent_mask is not None:
        mask &= np.asarray(parent_mask, dtype=bool)
    #endif

    return mask


def nsidis_driver_histogram(
    driver: str,
    feat: Dict[str, np.ndarray],
    mask: np.ndarray,
    *,
    mm2_min: float,
    mm2_max: float,
    mm2_bins: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
) -> np.ndarray:
    """
    Build a 2D histogram for every candidate, including 1D candidates.

    1D candidates are represented as shape (1,N) so the existing template
    morphing/profiling machinery can operate on their physical axis.
    """
    mx2 = np.asarray(feat["ep_missing_mass2"], dtype=float)
    pt = np.asarray(feat["stored_pTmiss"], dtype=float)
    th = np.asarray(feat["stored_theta_gamma_gamma"], dtype=float)
    dphi = np.asarray(
        feat["stored_delta_phi_abs_rad"],
        dtype=float,
    )

    if driver == "pTmiss":
        h, _ = np.histogram(
            pt[mask],
            bins=np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
        )
        return h.astype(float)[None, :]
    #endif

    if driver == "theta_gg":
        h, _ = np.histogram(
            th[mask],
            bins=np.linspace(0.0, theta_max, theta_bins + 1),
        )
        return h.astype(float)[None, :]
    #endif

    if driver == "delta_phi":
        h, _ = np.histogram(
            dphi[mask],
            bins=np.linspace(0.0, math.pi, 65),
        )
        return h.astype(float)[None, :]
    #endif

    if driver == "theta_gg_x_pTmiss":
        # pTmiss is intentionally the SECOND axis, so the existing nuisance
        # morphing shifts/smears the same pTmiss coordinate used in the current
        # Mx2 x pTmiss production model.
        h, _, _ = np.histogram2d(
            th[mask],
            pt[mask],
            bins=(
                np.linspace(0.0, theta_max, theta_bins + 1),
                np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            ),
        )
        return h.astype(float)
    #endif

    if driver == "delta_phi_x_pTmiss":
        # Keep pTmiss on the SECOND axis for the same morph-nuisance convention
        # as the other pTmiss-based candidates.
        h, _, _ = np.histogram2d(
            dphi[mask],
            pt[mask],
            bins=(
                np.linspace(0.0, math.pi, 65),
                np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            ),
        )
        return h.astype(float)
    #endif

    if driver == "mx2_ep_x_pTmiss":
        h, _, _ = np.histogram2d(
            mx2[mask],
            pt[mask],
            bins=(
                np.linspace(mm2_min, mm2_max, mm2_bins + 1),
                np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            ),
        )
        return h.astype(float)
    #endif

    if driver == "mx2_ep_x_theta_gg":
        h, _, _ = np.histogram2d(
            mx2[mask],
            th[mask],
            bins=(
                np.linspace(mm2_min, mm2_max, mm2_bins + 1),
                np.linspace(0.0, theta_max, theta_bins + 1),
            ),
        )
        return h.astype(float)
    #endif

    raise ValueError(f"Unknown nSidis driver-study discriminator '{driver}'")


def nsidis_driver_model_components(
    data_hist: np.ndarray,
    pi0_hist: np.ndarray,
    dvcs_hist: np.ndarray,
    fit: SharedMorphedFitResult,
    driver: str,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return fitted pi0, DVCS, and total model arrays in event counts."""
    if not fit.success or not fit.nuisance:
        shape = np.asarray(data_hist).shape
        nan = np.full(shape, np.nan, dtype=float)
        return nan, nan, nan
    #endif

    hp = morph_template_second_axis(
        np.asarray(pi0_hist, dtype=float),
        float(fit.nuisance.get(f"{driver}_pi0_shift_bins", 0.0)),
        float(fit.nuisance.get(f"{driver}_pi0_sigma_bins", 0.0)),
    )
    hv = morph_template_second_axis(
        np.asarray(dvcs_hist, dtype=float),
        float(fit.nuisance.get(f"{driver}_dvcs_shift_bins", 0.0)),
        float(fit.nuisance.get(f"{driver}_dvcs_sigma_bins", 0.0)),
    )

    ndata = float(np.sum(data_hist))
    fpi0 = float(fit.pi0_fraction)
    pi_comp = ndata * fpi0 * hp
    dv_comp = ndata * (1.0 - fpi0) * hv
    total = pi_comp + dv_comp
    return pi_comp, dv_comp, total


def run_nsidis_driver_study(
    period: PeriodConfig,
    data_f: Dict[str, np.ndarray],
    pi0_f: Dict[str, np.ndarray],
    dvcs_f: Dict[str, np.ndarray],
    *,
    parent_mask: Optional[np.ndarray],
    ft_theta_max: float,
    max_probe_energy: float,
    probe_m2_max: float,
    mm2_min: float,
    mm2_max: float,
    mm2_bins: int,
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
) -> Tuple[List[Dict[str, object]], Dict[Tuple[str, float, float, str], Dict[str, object]]]:
    """
    Compare five physically motivated denominator-composition drivers.

    All candidates use the same selected data/MC population in each
    region/energy bin.  The nominal production fit is NOT changed by this
    diagnostic study.
    """
    rows: List[Dict[str, object]] = []
    detail: Dict[
        Tuple[str, float, float, str],
        Dict[str, object],
    ] = {}

    for region in ("FT", "FD_all"):
        edges = nsidis_energy_edges_for_region(
            region, max_probe_energy
        )
        for ib in range(len(edges) - 1):
            elo = float(edges[ib])
            ehi = float(edges[ib + 1])

            masks = {
                "data": nsidis_driver_study_common_mask(
                    data_f, region, ft_theta_max, elo, ehi,
                    probe_m2_max, mm2_min, mm2_max,
                    ptmiss_max, theta_max, parent_mask,
                ),
                "pi0": nsidis_driver_study_common_mask(
                    pi0_f, region, ft_theta_max, elo, ehi,
                    probe_m2_max, mm2_min, mm2_max,
                    ptmiss_max, theta_max, None,
                ),
                "dvcs": nsidis_driver_study_common_mask(
                    dvcs_f, region, ft_theta_max, elo, ehi,
                    probe_m2_max, mm2_min, mm2_max,
                    ptmiss_max, theta_max, None,
                ),
            }

            mean_e = masked_finite_mean(
                data_f["pred_probe_energy"], masks["data"]
            )
            n_data = int(np.count_nonzero(masks["data"]))
            n_pi0 = int(np.count_nonzero(masks["pi0"]))
            n_dvcs = int(np.count_nonzero(masks["dvcs"]))

            for model_name, model_label, driver in NSIDIS_DRIVER_STUDY_MODELS:
                row: Dict[str, object] = {
                    "period": period.key,
                    "label": period.label,
                    "region": region,
                    "energy_low_GeV": elo,
                    "energy_high_GeV": ehi,
                    "mean_probe_energy_GeV": mean_e,
                    "model": model_name,
                    "model_label": model_label,
                    "driver": driver,
                    "common_population_data_count": n_data,
                    "common_population_pi0_count": n_pi0,
                    "common_population_dvcs_count": n_dvcs,
                }

                hd = nsidis_driver_histogram(
                    driver, data_f, masks["data"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    mm2_bins=mm2_bins,
                    ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max, theta_bins=theta_bins,
                )
                hp = nsidis_driver_histogram(
                    driver, pi0_f, masks["pi0"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    mm2_bins=mm2_bins,
                    ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max, theta_bins=theta_bins,
                )
                hv = nsidis_driver_histogram(
                    driver, dvcs_f, masks["dvcs"],
                    mm2_min=mm2_min, mm2_max=mm2_max,
                    mm2_bins=mm2_bins,
                    ptmiss_max=ptmiss_max, ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max, theta_bins=theta_bins,
                )

                if (
                    np.sum(hd) < min_data_count
                    or np.sum(hp) < min_template_count
                    or np.sum(hv) < min_template_count
                ):
                    row.update({
                        "fit_success": 0,
                        "fit_message": "insufficient common-population statistics",
                        "pi0_fraction": float("nan"),
                        "pi0_fraction_err_low": float("nan"),
                        "pi0_fraction_err_high": float("nan"),
                        "quality_full_deviance_per_active_bin": float("nan"),
                        "quality_axis0_deviance_per_active_bin": float("nan"),
                        "quality_axis1_deviance_per_active_bin": float("nan"),
                    })
                    rows.append(row)
                    continue
                #endif

                fit = fit_shared_morphed_composition(
                    {driver: (hd, hp, hv)},
                    {},
                    nuisance_shift_prior,
                    nuisance_sigma_prior,
                    max_shift_bins,
                    max_sigma_bins,
                    driver_names=(driver,),
                )

                row.update({
                    "fit_success": int(fit.success),
                    "fit_message": fit.message,
                    "pi0_fraction": float(fit.pi0_fraction),
                    "dvcs_fraction": (
                        1.0 - float(fit.pi0_fraction)
                        if np.isfinite(fit.pi0_fraction)
                        else float("nan")
                    ),
                    "fit_poisson_deviance": float(fit.poisson_deviance),
                    "fit_ndof_legacy": int(fit.ndof),
                })

                if fit.success:
                    unc = profile_pi0_fraction_uncertainty_one_driver(
                        hd, hp, hv, fit,
                        discriminator=driver,
                        nuisance_shift_prior=nuisance_shift_prior,
                        nuisance_sigma_prior=nuisance_sigma_prior,
                        max_shift_bins=max_shift_bins,
                        max_sigma_bins=max_sigma_bins,
                    )
                    row.update(unc)

                    pi_comp, dv_comp, total = (
                        nsidis_driver_model_components(
                            hd, hp, hv, fit, driver
                        )
                    )
                    q_full = binned_model_quality(hd, total)
                    q0 = binned_model_quality(
                        np.sum(hd, axis=1),
                        np.sum(total, axis=1),
                    )
                    q1 = binned_model_quality(
                        np.sum(hd, axis=0),
                        np.sum(total, axis=0),
                    )

                    row.update({
                        "quality_full_deviance_per_active_bin": float(
                            q_full["deviance_per_active_bin"]
                        ),
                        "quality_full_pearson_chi2_per_bin_mu_ge5": float(
                            q_full["pearson_chi2_per_bin_mu_ge5"]
                        ),
                        "quality_axis0_deviance_per_active_bin": float(
                            q0["deviance_per_active_bin"]
                        ),
                        "quality_axis0_pearson_chi2_per_bin_mu_ge5": float(
                            q0["pearson_chi2_per_bin_mu_ge5"]
                        ),
                        "quality_axis1_deviance_per_active_bin": float(
                            q1["deviance_per_active_bin"]
                        ),
                        "quality_axis1_pearson_chi2_per_bin_mu_ge5": float(
                            q1["pearson_chi2_per_bin_mu_ge5"]
                        ),
                    })
                    if fit.nuisance:
                        row.update(fit.nuisance)
                    #endif

                    detail[
                        (
                            region,
                            round(elo, 9),
                            round(ehi, 9),
                            model_name,
                        )
                    ] = {
                        "row": row,
                        "data_hist": hd,
                        "pi0_hist": hp,
                        "dvcs_hist": hv,
                        "pi0_component": pi_comp,
                        "dvcs_component": dv_comp,
                        "total_model": total,
                        "driver": driver,
                    }
                else:
                    row.update({
                        "pi0_fraction_err_low": float("nan"),
                        "pi0_fraction_err_high": float("nan"),
                        "quality_full_deviance_per_active_bin": float("nan"),
                        "quality_axis0_deviance_per_active_bin": float("nan"),
                        "quality_axis1_deviance_per_active_bin": float("nan"),
                    })
                #endif

                rows.append(row)
            #endfor
        #endfor
    #endfor

    return rows, detail


def make_nsidis_driver_study_summary_canvas(
    period: PeriodConfig,
    rows: Sequence[Dict[str, object]],
    outdir: Path,
) -> None:
    """
    Compact 2x2 variable-selection summary.

    Left column: fitted f_pi0.
    Right column: descriptive full-histogram D/active-bin.
    Top: FT. Bottom: FD_all.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16.0, 10.5))

    for irow, region in enumerate(("FT", "FD_all")):
        marker_cycle = ("o", "s", "^", "D", "v", "P", "X")
        for imodel, (model_name, model_label, _driver) in enumerate(
            NSIDIS_DRIVER_STUDY_MODELS
        ):
            rr = sorted(
                [
                    row for row in rows
                    if row.get("region") == region
                    and row.get("model") == model_name
                    and int(row.get("fit_success", 0)) == 1
                ],
                key=row_energy_coordinate,
            )
            if not rr:
                continue
            #endif

            x = np.asarray(
                [row_energy_coordinate(row) for row in rr],
                dtype=float,
            )
            y = np.asarray(
                [float(row["pi0_fraction"]) for row in rr],
                dtype=float,
            )
            el = np.asarray(
                [
                    float(row.get("pi0_fraction_err_low", np.nan))
                    for row in rr
                ],
                dtype=float,
            )
            eh = np.asarray(
                [
                    float(row.get("pi0_fraction_err_high", np.nan))
                    for row in rr
                ],
                dtype=float,
            )
            q = np.asarray(
                [
                    float(row.get(
                        "quality_full_deviance_per_active_bin",
                        np.nan,
                    ))
                    for row in rr
                ],
                dtype=float,
            )

            axes[irow, 0].errorbar(
                x,
                y,
                yerr=np.vstack((el, eh)),
                marker=marker_cycle[imodel % len(marker_cycle)],
                linewidth=1.0,
                capsize=2,
                label=model_label,
            )
            axes[irow, 1].plot(
                x,
                q,
                marker=marker_cycle[imodel % len(marker_cycle)],
                linewidth=1.0,
                label=model_label,
            )
        #endfor

        goodness_values = []
        for line in axes[irow, 1].get_lines():
            yy = np.asarray(line.get_ydata(), dtype=float)
            goodness_values.extend(
                [float(v) for v in yy if np.isfinite(v) and v > 0.0]
            )
        #endfor
        if goodness_values:
            gmin = min(goodness_values)
            gmax = max(goodness_values)
            if gmax / max(gmin, 1.0e-12) > 8.0:
                axes[irow, 1].set_yscale("linear")
            #endif
        #endif

        axes[irow, 0].set_ylim(0.0, 1.0)
        axes[irow, 0].set_ylabel(r"fitted $f_{\pi^0}$")
        axes[irow, 1].set_ylabel(
            "Poisson deviance / active bin"
        )

        title = "FT" if region == "FT" else "FD all"
        axes[irow, 0].set_title(
            f"{title}: composition stability"
        )
        axes[irow, 1].set_title(
            f"{title}: absolute template goodness"
        )

        for ax in axes[irow, :]:
            ax.set_xlabel(
                r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
            )
            ax.grid(alpha=0.18)
        #endfor
    #endfor

    axes[0, 0].legend(fontsize=8, frameon=False, ncol=2)
    axes[0, 1].legend(fontsize=8, frameon=False, ncol=2)

    fig.suptitle(
        f"{period.label}: dedicated denominator-driver study\n"
        "all methods fit the same common population; "
        r"$|\Delta\phi|$ is folded to $[0,\pi]$; "
        r"$E_{\mathrm{miss}}$ excluded because it is "
        r"$E_{\mathrm{probe}}^{\mathrm{pred}}$",
        fontsize=12.5,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_nsidis_driver_study_summary.png",
        rect=(0, 0, 1, 0.94),
    )
    plt.close(fig)


def _driver_study_axis_edges(
    driver: str,
    *,
    mm2_min: float,
    mm2_max: float,
    mm2_bins: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
) -> Tuple[np.ndarray, np.ndarray, str, str]:
    if driver == "pTmiss":
        return (
            np.asarray([0.0, 1.0]),
            np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            "",
            r"$p_{T,\mathrm{miss}}$ (GeV)",
        )
    #endif
    if driver == "theta_gg":
        return (
            np.asarray([0.0, 1.0]),
            np.linspace(0.0, theta_max, theta_bins + 1),
            "",
            r"$\theta_{\gamma\gamma}$ (deg)",
        )
    #endif
    if driver == "delta_phi":
        return (
            np.asarray([0.0, 1.0]),
            np.linspace(0.0, math.pi, 65),
            "",
            r"$|\Delta\phi|$ (rad)",
        )
    #endif
    if driver == "theta_gg_x_pTmiss":
        return (
            np.linspace(0.0, theta_max, theta_bins + 1),
            np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            r"$\theta_{\gamma\gamma}$ (deg)",
            r"$p_{T,\mathrm{miss}}$ (GeV)",
        )
    #endif
    if driver == "delta_phi_x_pTmiss":
        return (
            np.linspace(0.0, math.pi, 65),
            np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            r"$|\Delta\phi|$ (rad)",
            r"$p_{T,\mathrm{miss}}$ (GeV)",
        )
    #endif
    if driver == "mx2_ep_x_pTmiss":
        return (
            np.linspace(mm2_min, mm2_max, mm2_bins + 1),
            np.linspace(0.0, ptmiss_max, ptmiss_bins + 1),
            r"$M_X^2(ep)$ (GeV$^2$)",
            r"$p_{T,\mathrm{miss}}$ (GeV)",
        )
    #endif
    if driver == "mx2_ep_x_theta_gg":
        return (
            np.linspace(mm2_min, mm2_max, mm2_bins + 1),
            np.linspace(0.0, theta_max, theta_bins + 1),
            r"$M_X^2(ep)$ (GeV$^2$)",
            r"$\theta_{\gamma\gamma}$ (deg)",
        )
    #endif
    raise ValueError(driver)


def make_nsidis_driver_study_focus_canvas(
    period: PeriodConfig,
    detail: Dict[
        Tuple[str, float, float, str],
        Dict[str, object],
    ],
    outdir: Path,
    *,
    focus_low: float,
    focus_high: float,
    mm2_min: float,
    mm2_max: float,
    mm2_bins: int,
    ptmiss_max: float,
    ptmiss_bins: int,
    theta_max: float,
    theta_bins: int,
) -> None:
    """
    Readable visual inspection of the currently suspicious 3.0-4.5 GeV bin.

    Write one canvas per detector region rather than squeezing FT and FD into
    four columns.  Rows are candidate models; columns are the two physical
    projections of a 2D fit.  For a 1D candidate the right panel is left blank.
    """
    nrows = len(NSIDIS_DRIVER_STUDY_MODELS)

    for region in ("FT", "FD_all"):
        fig, axes = plt.subplots(
            nrows,
            2,
            figsize=(13.5, 3.0 * nrows),
            squeeze=False,
        )

        for irow, (model_name, model_label, driver) in enumerate(
            NSIDIS_DRIVER_STUDY_MODELS
        ):
            entry = detail.get(
                (
                    region,
                    round(float(focus_low), 9),
                    round(float(focus_high), 9),
                    model_name,
                )
            )
            ax0 = axes[irow, 0]
            ax1 = axes[irow, 1]

            if entry is None:
                ax0.text(
                    0.5, 0.5, "fit unavailable",
                    transform=ax0.transAxes,
                    ha="center", va="center",
                )
                ax1.axis("off")
                continue
            #endif

            hd = np.asarray(entry["data_hist"], dtype=float)
            hp = np.asarray(entry["pi0_component"], dtype=float)
            hv = np.asarray(entry["dvcs_component"], dtype=float)
            ht = np.asarray(entry["total_model"], dtype=float)
            row = entry["row"]

            xedges, yedges, xlabel0, xlabel1 = (
                _driver_study_axis_edges(
                    driver,
                    mm2_min=mm2_min,
                    mm2_max=mm2_max,
                    mm2_bins=mm2_bins,
                    ptmiss_max=ptmiss_max,
                    ptmiss_bins=ptmiss_bins,
                    theta_max=theta_max,
                    theta_bins=theta_bins,
                )
            )

            def draw_projection(
                ax,
                edges: np.ndarray,
                data: np.ndarray,
                dvcs: np.ndarray,
                pi0: np.ndarray,
                total: np.ndarray,
                xlabel: str,
            ) -> None:
                centers = 0.5 * (edges[:-1] + edges[1:])
                ax.errorbar(
                    centers,
                    data,
                    yerr=np.sqrt(np.clip(data, 0.0, None)),
                    fmt="o",
                    ms=2.6,
                    linewidth=0.7,
                    label="data",
                )
                ax.step(
                    edges[:-1], dvcs,
                    where="post", linewidth=1.0, label="BH/DVCS"
                )
                ax.step(
                    edges[:-1], pi0,
                    where="post", linewidth=1.0, label=r"$\pi^0$"
                )
                ax.step(
                    edges[:-1], total,
                    where="post", linewidth=1.3, label="total"
                )
                ax.set_xlabel(xlabel)
                ax.set_ylabel("entries / bin")
                ax.grid(alpha=0.15)
            #enddef

            if hd.shape[0] == 1:
                draw_projection(
                    ax0,
                    yedges,
                    hd[0],
                    hv[0],
                    hp[0],
                    ht[0],
                    xlabel1,
                )
                ax1.axis("off")
            else:
                draw_projection(
                    ax0,
                    xedges,
                    np.sum(hd, axis=1),
                    np.sum(hv, axis=1),
                    np.sum(hp, axis=1),
                    np.sum(ht, axis=1),
                    xlabel0,
                )
                draw_projection(
                    ax1,
                    yedges,
                    np.sum(hd, axis=0),
                    np.sum(hv, axis=0),
                    np.sum(hp, axis=0),
                    np.sum(ht, axis=0),
                    xlabel1,
                )
            #endif

            fpi0 = float(row.get("pi0_fraction", np.nan))
            elo = float(row.get("pi0_fraction_err_low", np.nan))
            ehi = float(row.get("pi0_fraction_err_high", np.nan))
            q = float(
                row.get(
                    "quality_full_deviance_per_active_bin",
                    np.nan,
                )
            )
            region_label = "FT" if region == "FT" else "FD all"
            title = (
                f"{region_label}: {model_label}; "
                rf"$f_{{\pi^0}}={fpi0:.3f}"
            )
            if np.isfinite(elo) and np.isfinite(ehi):
                title += rf"^{{+{ehi:.3f}}}_{{-{elo:.3f}}}"
            #endif
            title += rf"$; D/active={q:.2f}"
            ax0.set_title(title, fontsize=9.5)
            if ax1.axison:
                ax1.set_title("same fit: second-axis projection", fontsize=9)
            #endif
        #endfor

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                fontsize=8,
                frameon=False,
                ncol=4,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.985),
            )
        #endif

        region_slug = "ft" if region == "FT" else "fd_all"
        region_label = "FT" if region == "FT" else "FD all"
        fig.suptitle(
            f"{period.label}: denominator-driver visual check — {region_label}\n"
            f"{focus_low:.1f}-{focus_high:.1f} GeV; "
            "same common population for every candidate",
            fontsize=13,
            y=0.998,
        )
        safe_finalize_figure(
            fig,
            outdir
            / f"canvas_nsidis_driver_study_focus_3p0_4p5_{region_slug}.png",
            rect=(0, 0, 1, 0.975),
        )
        plt.close(fig)
    #endfor



def make_cross_period_nsidis_driver_study_canvas(
    nsroot: Path,
    periods: Sequence[PeriodConfig],
) -> None:
    """
    Cross-period stability test for candidate denominator methods.

    Plot |Delta f_pi0| between the first two available periods as a function
    of actual mean E_probe^pred.  The purpose is not to choose a method merely
    because it forces periods to agree; this is a robustness diagnostic used
    together with absolute template goodness.
    """
    period_rows: Dict[str, List[Dict[str, str]]] = {}
    labels: Dict[str, str] = {}

    for period in periods:
        path = (
            nsroot / period.key
            / "nsidis_driver_study.csv"
        )
        if not path.exists():
            continue
        #endif
        rows = read_csv_rows_simple(path)
        if rows:
            period_rows[period.key] = rows
            labels[period.key] = period.label
        #endif
    #endfor

    available = [
        period.key for period in periods
        if period.key in period_rows
    ]
    if len(available) < 2:
        return
    #endif

    p0, p1 = available[:2]
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    for ax, region in zip(axes, ("FT", "FD_all")):
        for model_name, model_label, _driver in NSIDIS_DRIVER_STUDY_MODELS:
            map0 = {}
            for row in period_rows[p0]:
                if (
                    row.get("region") == region
                    and row.get("model") == model_name
                    and row.get("fit_success") == "1"
                ):
                    key = (
                        round(float(row["energy_low_GeV"]), 9),
                        round(float(row["energy_high_GeV"]), 9),
                    )
                    map0[key] = row
                #endif
            #endfor

            x = []
            y = []
            for row in period_rows[p1]:
                if not (
                    row.get("region") == region
                    and row.get("model") == model_name
                    and row.get("fit_success") == "1"
                ):
                    continue
                #endif
                key = (
                    round(float(row["energy_low_GeV"]), 9),
                    round(float(row["energy_high_GeV"]), 9),
                )
                r0 = map0.get(key)
                if r0 is None:
                    continue
                #endif
                f0 = float(r0["pi0_fraction"])
                f1 = float(row["pi0_fraction"])
                e0 = row_energy_coordinate(r0)
                e1 = row_energy_coordinate(row)
                if np.isfinite(f0) and np.isfinite(f1):
                    x.append(0.5 * (e0 + e1))
                    y.append(abs(f0 - f1))
                #endif
            #endfor

            if x:
                ax.plot(
                    x, y,
                    marker="o",
                    linewidth=1.0,
                    label=model_label,
                )
            #endif
        #endfor

        ax.set_xlabel(
            r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
        )
        ax.set_ylabel(
            rf"$|f_{{\pi^0}}^{{{labels[p0]}}}"
            rf"-f_{{\pi^0}}^{{{labels[p1]}}}|$"
        )
        ax.set_title("FT" if region == "FT" else "FD all")
        ax.grid(alpha=0.18)
        ax.legend(fontsize=7.5, frameon=False)
    #endfor

    fig.suptitle(
        "nSidis denominator-driver cross-period stability\n"
        "use together with absolute goodness; smaller period difference alone "
        "does not define the preferred model",
        fontsize=12.5,
    )
    safe_finalize_figure(
        fig,
        nsroot / "canvas_nsidis_driver_study_cross_period.png",
        rect=(0, 0, 1, 0.93),
    )
    plt.close(fig)




# -----------------------------------------------------------------------------
# Associated-numerator pi0-purity diagnostic
# -----------------------------------------------------------------------------

def _normalize_template_1d(counts: np.ndarray) -> np.ndarray:
    arr = np.asarray(counts, dtype=float)
    total = float(np.sum(arr))
    if not np.isfinite(total) or total <= 0.0:
        return np.zeros_like(arr, dtype=float)
    #endif
    return arr / total


def _numerator_background_shape(
    edges: np.ndarray,
    slope: float,
    curvature: float = 0.0,
) -> np.ndarray:
    """
    Positive smooth background over the already-skimmed 0.11-0.16 GeV window.

    u is centered and scaled to roughly [-0.5,0.5].  The nominal model has
    curvature=0 (linear); a quadratic alternative is used only as a diagnostic
    model-dependence check.
    """
    centers = 0.5 * (edges[:-1] + edges[1:])
    width = max(float(edges[-1] - edges[0]), 1.0e-12)
    u = (centers - 0.5 * (edges[0] + edges[-1])) / width
    shape = 1.0 + float(slope) * u + float(curvature) * u * u
    shape = np.clip(shape, 1.0e-10, None)
    return _normalize_template_1d(shape)


def _numerator_mass_nll(
    params: np.ndarray,
    observed: np.ndarray,
    signal_shape: np.ndarray,
    edges: np.ndarray,
    background_order: int,
    fixed_purity: Optional[float] = None,
) -> float:
    obs = np.asarray(observed, dtype=float)
    nobs = float(np.sum(obs))
    if nobs <= 0.0:
        return float("inf")
    #endif

    order = int(background_order)
    if fixed_purity is None:
        purity = float(params[0])
        nuisance = np.asarray(params[1:], dtype=float)
    else:
        purity = float(fixed_purity)
        nuisance = np.asarray(params, dtype=float)
    #endif

    if not (0.0 <= purity <= 1.0):
        return float("inf")
    #endif

    if order <= 0:
        slope = 0.0
        curvature = 0.0
    elif order == 1:
        slope = float(nuisance[0]) if len(nuisance) else 0.0
        curvature = 0.0
    else:
        raise ValueError(
            "Associated-numerator production study supports only "
            "constant (order 0) and linear (order 1) backgrounds."
        )
    #endif

    bkg = _numerator_background_shape(edges, slope, curvature)
    model_shape = (
        purity * np.asarray(signal_shape, dtype=float)
        + (1.0 - purity) * bkg
    )
    mu = np.clip(nobs * model_shape, 1.0e-12, None)
    return float(np.sum(mu - obs * np.log(mu)))


def fit_associated_numerator_mass(
    data_counts: np.ndarray,
    aaogen_counts: np.ndarray,
    edges: np.ndarray,
    *,
    background_order: int = 1,
) -> Dict[str, float]:
    """
    Fit the FINAL associated-data M_gg spectrum with an associated-aaogen
    pi0 signal template plus a smooth combinatorial background.

    background_order:
      0 -> constant background
      1 -> linear background (NOMINAL)

    The 68% purity interval is a one-dimensional profile-likelihood interval
    using Delta(-ln L)=0.5.  The result now feeds the data numerator efficiency.
    """
    data_counts = np.asarray(data_counts, dtype=float)
    aaogen_counts = np.asarray(aaogen_counts, dtype=float)
    edges = np.asarray(edges, dtype=float)

    ndata = float(np.sum(data_counts))
    nmc = float(np.sum(aaogen_counts))
    if ndata < 20.0 or nmc < 20.0:
        return {
            "fit_success": 0,
            "fit_message": "insufficient associated mass statistics",
            "purity": float("nan"),
            "purity_err_low": float("nan"),
            "purity_err_high": float("nan"),
            "background_slope": float("nan"),
            "nll": float("nan"),
        }
    #endif

    signal_shape = _normalize_template_1d(aaogen_counts)
    if not np.any(signal_shape > 0.0):
        return {
            "fit_success": 0,
            "fit_message": "empty aaogen signal template",
            "purity": float("nan"),
            "purity_err_low": float("nan"),
            "purity_err_high": float("nan"),
            "background_slope": float("nan"),
            "nll": float("nan"),
        }
    #endif

    order = int(background_order)
    if order <= 0:
        x0 = np.asarray([0.90], dtype=float)
        bounds = ((0.0, 1.0),)
    elif order == 1:
        x0 = np.asarray([0.90, 0.0], dtype=float)
        bounds = ((0.0, 1.0), (-1.9, 1.9))
    else:
        raise ValueError(
            "Use background_order=0 (constant) or 1 (linear)."
        )
    #endif

    result = minimize(
        _numerator_mass_nll,
        x0,
        args=(
            data_counts,
            signal_shape,
            edges,
            order,
            None,
        ),
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 1000, "ftol": 1.0e-12},
    )

    if not result.success or not np.isfinite(result.fun):
        return {
            "fit_success": 0,
            "fit_message": str(result.message),
            "purity": float("nan"),
            "purity_err_low": float("nan"),
            "purity_err_high": float("nan"),
            "background_slope": float("nan"),
            "nll": float("nan"),
        }
    #endif

    purity = float(result.x[0])
    slope = float(result.x[1]) if order == 1 else 0.0
    nll_min = float(result.fun)

    def profile_nll(p: float) -> float:
        if order <= 0:
            return _numerator_mass_nll(
                np.asarray([], dtype=float),
                data_counts,
                signal_shape,
                edges,
                order,
                float(p),
            )
        #endif

        prof = minimize(
            _numerator_mass_nll,
            np.asarray([slope], dtype=float),
            args=(
                data_counts,
                signal_shape,
                edges,
                order,
                float(p),
            ),
            method="L-BFGS-B",
            bounds=((-1.9, 1.9),),
            options={"maxiter": 500, "ftol": 1.0e-11},
        )
        return float(prof.fun) if prof.success else float("inf")
    #enddef

    target = nll_min + 0.5

    def root_function(p: float) -> float:
        return profile_nll(float(p)) - target
    #enddef

    p_low = 0.0
    if purity > 0.0:
        try:
            f0 = root_function(0.0)
            fp = root_function(purity)
            if np.isfinite(f0) and f0 > 0.0 and fp <= 0.0:
                p_low = float(brentq(root_function, 0.0, purity))
            #endif
        except Exception:
            p_low = 0.0
        #endtry
    #endif

    p_high = 1.0
    if purity < 1.0:
        try:
            fp = root_function(purity)
            f1 = root_function(1.0)
            if np.isfinite(f1) and f1 > 0.0 and fp <= 0.0:
                p_high = float(brentq(root_function, purity, 1.0))
            #endif
        except Exception:
            p_high = 1.0
        #endtry
    #endif

    return {
        "fit_success": 1,
        "fit_message": "ok",
        "purity": purity,
        "purity_err_low": float(max(0.0, purity - p_low)),
        "purity_err_high": float(max(0.0, p_high - purity)),
        "background_slope": slope,
        "nll": nll_min,
    }



def _aaogen_final_association_mask(
    pair_np: Dict[str, np.ndarray],
    *,
    mgg_min: float,
    mgg_max: float,
    remainder_mass2_max: float,
    reco_probe_energy_min: float,
    probe_angle_max_deg: float,
    probe_frac_energy_max: float,
) -> np.ndarray:
    """Pair-level version of the final aaogen association definition."""
    if not pair_np:
        return np.zeros(0, dtype=bool)
    #endif

    pi0_mass = np.asarray(pair_np["pi0_mass"], dtype=float)
    reco_E = np.asarray(pair_np["reco_probe_energy"], dtype=float)
    reco_p = np.asarray(pair_np["reco_probe_p"], dtype=float)
    reco_m2 = np.asarray(pair_np["reco_probe_mass2"], dtype=float)
    pred_E = np.asarray(pair_np["pred_probe_energy"], dtype=float)
    opening = np.asarray(
        pair_np["probe_opening_residual_deg"], dtype=float
    )
    frac_E = np.asarray(pair_np["probe_delta_E_over_E"], dtype=float)

    return (
        np.isfinite(pi0_mass)
        & (pi0_mass >= float(mgg_min))
        & (pi0_mass <= float(mgg_max))
        & np.isfinite(reco_E)
        & np.isfinite(reco_p)
        & (reco_E >= float(reco_probe_energy_min))
        & (reco_p > 0.0)
        & np.isfinite(reco_m2)
        & (np.abs(reco_m2) <= float(remainder_mass2_max))
        & np.isfinite(pred_E)
        & (pred_E > 0.0)
        & np.isfinite(opening)
        & (opening <= float(probe_angle_max_deg))
        & np.isfinite(frac_E)
        & (np.abs(frac_E) <= float(probe_frac_energy_max))
    )



def association_stage_rows(
    period: PeriodConfig,
    label: str,
    assoc: DataAssociationResult,
) -> List[Dict[str, object]]:
    """One compact table of exact-association survival fractions."""
    counters = assoc.counters
    n_tags = int(counters.get("epgamma_entries", 0))
    stages = (
        ("same_event", "same_event"),
        ("positive_remainder", "positive_remainder"),
        ("tag_mass_shell", "tag_mass_shell"),
        ("probe_energy", "probe_energy"),
        ("probe_pred_consistent", "probe_pred_consistent"),
    )
    rows: List[Dict[str, object]] = []
    previous = n_tags
    for stage_name, counter_key in stages:
        count = int(counters.get(counter_key, 0))
        rows.append(
            {
                "period": period.key,
                "label": period.label,
                "sample": label,
                "stage": stage_name,
                "count": count,
                "fraction_of_tags": (
                    float(count) / float(n_tags)
                    if n_tags > 0 else float("nan")
                ),
                "conditional_survival": (
                    float(count) / float(previous)
                    if previous > 0 else float("nan")
                ),
            }
        )
        previous = count
    #endfor
    return rows


def log_association_summary(
    period: PeriodConfig,
    label: str,
    assoc: DataAssociationResult,
) -> None:
    c = assoc.counters
    log(
        f"{period.label}: {label} association stages: "
        f"tags={int(c.get('epgamma_entries', 0)):,}, "
        f"same-event={int(c.get('same_event', 0)):,}, "
        f"positive={int(c.get('positive_remainder', 0)):,}, "
        f"mass-shell={int(c.get('tag_mass_shell', 0)):,}, "
        f"Eprobe={int(c.get('probe_energy', 0)):,}, "
        f"pred-consistent={int(c.get('probe_pred_consistent', 0)):,}."
    )



def build_associated_numerator_purity_rows(
    period: PeriodConfig,
    data_assoc: DataAssociationResult,
    mc_pair_np: Dict[str, np.ndarray],
    pi0_features: Dict[str, np.ndarray],
    *,
    central_pi0_tag_mask: np.ndarray,
    max_probe_energy: float,
    mgg_min: float,
    mgg_max: float,
    remainder_mass2_max: float,
    reco_probe_energy_min: float,
    probe_angle_max_deg: float,
    probe_frac_energy_max: float,
    mass_bins: int = 40,
) -> Tuple[
    List[Dict[str, object]],
    Dict[Tuple[str, float, float], Dict[str, np.ndarray]],
]:
    """
    Extract associated-numerator pi0 purity versus probe energy for FT/FD_all.

    Data:
      use the FINAL best same-event eppi0 association and require
      passes_pred_consistency.

    Signal template:
      use reconstructed aaogen eppi0 pairs passing the analogous final
      association definition and whose aaogen epgamma tag belongs to the same
      central denominator support.

    No purity correction is applied to the efficiency in this version.
    """
    best = data_assoc.best_diagnostics
    required_data = (
        "pi0_mass",
        "pred_probe_energy",
        "pred_probe_theta_deg",
        "passes_pred_consistency",
    )
    missing = [name for name in required_data if name not in best]
    if missing:
        raise RuntimeError(
            f"{period.label}: numerator-purity data association is missing "
            f"diagnostics: {', '.join(missing)}"
        )
    #endif

    data_mass = np.asarray(best["pi0_mass"], dtype=float)
    data_E = np.asarray(best["pred_probe_energy"], dtype=float)
    data_theta = np.asarray(best["pred_probe_theta_deg"], dtype=float)
    data_final = np.asarray(
        best["passes_pred_consistency"], dtype=bool
    )
    data_association_empty = (len(data_mass) == 0)

    mc_final = _aaogen_final_association_mask(
        mc_pair_np,
        mgg_min=mgg_min,
        mgg_max=mgg_max,
        remainder_mass2_max=remainder_mass2_max,
        reco_probe_energy_min=reco_probe_energy_min,
        probe_angle_max_deg=probe_angle_max_deg,
        probe_frac_energy_max=probe_frac_energy_max,
    )

    mc_epg_idx = np.asarray(mc_pair_np["epg_index"], dtype=np.int64)
    valid_mc_idx = (
        (mc_epg_idx >= 0)
        & (mc_epg_idx < len(central_pi0_tag_mask))
    )
    mc_final &= valid_mc_idx
    mc_final &= np.asarray(
        central_pi0_tag_mask[mc_epg_idx.clip(
            0, max(0, len(central_pi0_tag_mask) - 1)
        )],
        dtype=bool,
    )

    mc_mass = np.asarray(mc_pair_np["pi0_mass"], dtype=float)
    mc_E = np.asarray(mc_pair_np["pred_probe_energy"], dtype=float)
    mc_theta = np.asarray(
        mc_pair_np["pred_probe_theta_deg"], dtype=float
    )

    mass_edges = np.linspace(
        float(mgg_min),
        float(mgg_max),
        int(mass_bins) + 1,
    )

    rows: List[Dict[str, object]] = []
    detail: Dict[
        Tuple[str, float, float],
        Dict[str, np.ndarray],
    ] = {}

    for region in ("FT", "FD_all"):
        energy_edges = nsidis_energy_edges_for_region(
            region, max_probe_energy
        )
        for ib in range(len(energy_edges) - 1):
            elo = float(energy_edges[ib])
            ehi = float(energy_edges[ib + 1])

            dmask = (
                data_final
                & np.isfinite(data_mass)
                & np.isfinite(data_E)
                & np.isfinite(data_theta)
                & (data_mass >= mgg_min)
                & (data_mass <= mgg_max)
                & (data_E >= elo)
                & (data_E < ehi)
            )
            mmask = (
                mc_final
                & np.isfinite(mc_mass)
                & np.isfinite(mc_E)
                & np.isfinite(mc_theta)
                & (mc_mass >= mgg_min)
                & (mc_mass <= mgg_max)
                & (mc_E >= elo)
                & (mc_E < ehi)
            )

            if region == "FT":
                dmask &= (data_theta >= 2.4) & (data_theta < 5.0)
                mmask &= (mc_theta >= 2.4) & (mc_theta < 5.0)
            else:
                dmask &= (data_theta >= 6.0) & (data_theta < 35.0)
                mmask &= (mc_theta >= 6.0) & (mc_theta < 35.0)
            #endif

            data_hist, _ = np.histogram(
                data_mass[dmask], bins=mass_edges
            )
            mc_hist, _ = np.histogram(
                mc_mass[mmask], bins=mass_edges
            )

            fit_linear = fit_associated_numerator_mass(
                data_hist,
                mc_hist,
                mass_edges,
                background_order=1,
            )
            fit_constant = fit_associated_numerator_mass(
                data_hist,
                mc_hist,
                mass_edges,
                background_order=0,
            )

            mean_E = (
                float(np.mean(data_E[dmask]))
                if np.any(dmask)
                else float("nan")
            )

            p_lin = float(fit_linear["purity"])
            p_const = float(fit_constant["purity"])
            model_shift = (
                abs(p_lin - p_const)
                if np.isfinite(p_lin) and np.isfinite(p_const)
                else float("nan")
            )

            row: Dict[str, object] = {
                "period": period.key,
                "label": period.label,
                "region": region,
                "energy_low_GeV": elo,
                "energy_high_GeV": ehi,
                "mean_probe_energy_GeV": mean_E,
                "associated_data_count": int(np.sum(data_hist)),
                "associated_aaogen_count": int(np.sum(mc_hist)),
                "purity_linear": p_lin,
                "purity_linear_err_low": float(
                    fit_linear["purity_err_low"]
                ),
                "purity_linear_err_high": float(
                    fit_linear["purity_err_high"]
                ),
                "purity_constant": p_const,
                "purity_constant_err_low": float(
                    fit_constant["purity_err_low"]
                ),
                "purity_constant_err_high": float(
                    fit_constant["purity_err_high"]
                ),
                "background_model_abs_shift": model_shift,
                "linear_fit_success": int(
                    fit_linear["fit_success"]
                ),
                "constant_fit_success": int(
                    fit_constant["fit_success"]
                ),
                "association_status": (
                    "no_final_associated_candidates"
                    if data_association_empty
                    else "ok"
                ),
                "nominal_background_model": "linear",
                "purity_applied_to_efficiency": 1,
            }
            rows.append(row)

            if int(fit_linear["fit_success"]) == 1:
                signal_shape = _normalize_template_1d(mc_hist)
                bkg_shape = _numerator_background_shape(
                    mass_edges,
                    float(fit_linear["background_slope"]),
                    0.0,
                )
                ndata = float(np.sum(data_hist))
                sig_comp = ndata * p_lin * signal_shape
                bkg_comp = ndata * (1.0 - p_lin) * bkg_shape
                detail[(region, elo, ehi)] = {
                    "mass_edges": mass_edges,
                    "data_hist": data_hist.astype(float),
                    "aaogen_hist": mc_hist.astype(float),
                    "signal_component": sig_comp,
                    "background_component": bkg_comp,
                    "total_model": sig_comp + bkg_comp,
                }
            #endif
        #endfor
    #endfor

    return rows, detail


def make_associated_numerator_purity_summary_canvas(
    period: PeriodConfig,
    rows: Sequence[Dict[str, object]],
    outdir: Path,
) -> None:
    """One compact FT/FD summary of numerator purity versus energy."""
    fig, axes = plt.subplots(
        1, 2,
        figsize=(12.5, 4.8),
        sharey=True,
    )

    for ax, region in zip(axes, ("FT", "FD_all")):
        rr = sorted(
            [row for row in rows if row["region"] == region],
            key=row_energy_coordinate,
        )
        if rr:
            x = np.asarray(
                [row_energy_coordinate(row) for row in rr],
                dtype=float,
            )
            y = np.asarray(
                [float(row["purity_linear"]) for row in rr],
                dtype=float,
            )
            el = np.asarray(
                [
                    float(row["purity_linear_err_low"])
                    for row in rr
                ],
                dtype=float,
            )
            eh = np.asarray(
                [
                    float(row["purity_linear_err_high"])
                    for row in rr
                ],
                dtype=float,
            )
            yq = np.asarray(
                [float(row["purity_constant"]) for row in rr],
                dtype=float,
            )

            ax.errorbar(
                x,
                y,
                yerr=np.vstack((el, eh)),
                marker="o",
                linestyle="none",
                linewidth=1.1,
                capsize=2,
                label="aaogen template + linear background",
            )
            ax.plot(
                x,
                yq,
                marker="s",
                linestyle="none",
                linewidth=1.0,
                label="constant-background alternative",
            )
        #endif

        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(
            r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)"
        )
        ax.set_title("FT" if region == "FT" else "FD all")
        ax.grid(alpha=0.18)
        ax.legend(fontsize=8, frameon=False)
    #endfor

    axes[0].set_ylabel(
        r"final-associated numerator $\pi^0$ purity"
    )
    fig.suptitle(
        f"{period.label}: numerator $M_{{\\gamma\\gamma}}$ purity diagnostic\n"
        "linear-background purity is applied to data efficiency; "
        "constant-background result is retained as a model variation",
        fontsize=12.5,
    )
    safe_finalize_figure(
        fig,
        outdir / "canvas_numerator_pi0_purity_summary.png",
        rect=(0, 0, 1, 0.91),
    )
    plt.close(fig)


def make_associated_numerator_mass_fit_canvases(
    period: PeriodConfig,
    rows: Sequence[Dict[str, object]],
    detail: Dict[
        Tuple[str, float, float],
        Dict[str, np.ndarray],
    ],
    outdir: Path,
) -> None:
    """
    One readable all-energy M_gg fit canvas per detector region.

    Eight current energy bins fit naturally in a 4x2 canvas; if the binning
    changes, the layout grows automatically.
    """
    for region in ("FT", "FD_all"):
        rr = sorted(
            [row for row in rows if row["region"] == region],
            key=lambda r: float(r["energy_low_GeV"]),
        )
        if not rr:
            continue
        #endif

        ncols = min(4, max(1, len(rr)))
        nrows = int(math.ceil(len(rr) / ncols))
        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(4.7 * ncols, 3.7 * nrows + 0.8),
            squeeze=False,
        )

        for ip, row in enumerate(rr):
            ax = axes[ip // ncols, ip % ncols]
            key = (
                region,
                float(row["energy_low_GeV"]),
                float(row["energy_high_GeV"]),
            )
            entry = detail.get(key)
            if entry is None:
                ax.text(
                    0.5,
                    0.5,
                    "fit unavailable",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                )
                ax.set_axis_off()
                continue
            #endif

            edges = entry["mass_edges"]
            centers = 0.5 * (edges[:-1] + edges[1:])
            data_hist = entry["data_hist"]

            ax.errorbar(
                centers,
                data_hist,
                yerr=np.sqrt(np.clip(data_hist, 0.0, None)),
                fmt="o",
                ms=3,
                linewidth=0.7,
                label="associated data",
            )
            ax.step(
                edges[:-1],
                entry["signal_component"],
                where="post",
                linewidth=1.0,
                label=r"fitted $\pi^0$ signal",
            )
            ax.step(
                edges[:-1],
                entry["background_component"],
                where="post",
                linewidth=1.0,
                label="smooth background",
            )
            ax.step(
                edges[:-1],
                entry["total_model"],
                where="post",
                linewidth=1.3,
                label="total fit",
            )

            p = float(row["purity_linear"])
            elo = float(row["energy_low_GeV"])
            ehi = float(row["energy_high_GeV"])
            ax.set_title(
                f"{elo:.2f}-{ehi:.2f} GeV; "
                rf"$P_{{\pi^0}}={p:.3f}$",
                fontsize=9.5,
            )
            ax.set_xlabel(r"$M_{\gamma\gamma}$ (GeV)")
            ax.set_ylabel("entries / bin")
            ax.grid(alpha=0.16)
        #endfor

        for ip in range(len(rr), nrows * ncols):
            axes[ip // ncols, ip % ncols].set_axis_off()
        #endfor

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.915),
                ncol=4,
                frameon=False,
                fontsize=8,
            )
        #endif

        region_label = "FT" if region == "FT" else "FD all"
        region_slug = "ft" if region == "FT" else "fd_all"
        fig.suptitle(
            f"{period.label}: final-associated numerator "
            f"$M_{{\\gamma\\gamma}}$ fits — {region_label}\n"
            "aaogen signal template + linear combinatorial background "
            "(constant background is the alternative)",
            fontsize=12.5,
            y=0.985,
        )
        safe_finalize_figure(
            fig,
            outdir
            / f"canvas_numerator_pi0_mass_fits_{region_slug}.png",
            rect=(0.0, 0.0, 1.0, 0.84),
        )
        plt.close(fig)
    #endfor



def preflight_nsidis_study(
    periods: Sequence[PeriodConfig],
    args_dict: Dict[str, object],
) -> None:
    """Fail before the long run if paths or required event IDs are missing."""
    missing = []
    schema_problems = []

    for period in periods:
        ns_epg, ns_epi = resolve_nsidis_paths(period, args_dict)

        if bool(args_dict.get("stage1_only", False)):
            checks = [
                ("nSidis epgamma", ns_epg, EPG_REQUIRED),
                ("aaogen epgamma", period.epgamma_pi0_mc, EPG_REQUIRED),
                ("dvcsgen epgamma", period.epgamma_dvcs_mc, EPG_REQUIRED),
            ]
        else:
            checks = [
                ("old epgamma", period.epgamma_data, EPG_REQUIRED),
                ("old eppi0", period.eppi0_data, EPPIO_REQUIRED),
                ("nSidis epgamma", ns_epg, EPG_REQUIRED),
                ("nSidis eppi0", ns_epi, EPPIO_REQUIRED),
            ]
            if bool(args_dict.get("nsidis_pilot_fit", False)):
                checks.extend([
                    ("aaogen epgamma", period.epgamma_pi0_mc, EPG_REQUIRED),
                    ("dvcsgen epgamma", period.epgamma_dvcs_mc, EPG_REQUIRED),
                ])
            #endif
        #endif

        if period.clasdis_epgamma_mc is not None:
            clasdis_path = Path(period.clasdis_epgamma_mc)
            if clasdis_path.exists():
                checks.append(
                    (
                        "clasdis epgammaX",
                        str(clasdis_path),
                        EPG_REQUIRED,
                    )
                )
            #endif
        #endif

        for label, path, required in checks:
            if not Path(path).exists():
                missing.append(f"{period.label}: {label}: {path}")
                continue
            #endif

            needs_event_ids = (
                not bool(args_dict.get("stage1_only", False))
                and label in (
                    "old epgamma",
                    "old eppi0",
                    "nSidis epgamma",
                    "nSidis eppi0",
                )
            )

            with uproot.open(path) as root_file:
                found, tree = find_tree(
                    root_file,
                    args_dict.get("tree"),
                )
                available = set(tree.keys())
                needed = set(required)
                if needs_event_ids:
                    needed |= {"runnum", "evnum"}
                #endif
                if label == "nSidis epgamma":
                    needed |= {"detector2"}
                #endif
                absent = sorted(needed - available)
                if absent:
                    schema_problems.append(
                        f"{period.label}: {label}: tree '{found}' missing "
                        + ", ".join(absent)
                    )
                #endif
            #endwith
        #endfor
    #endfor

    if missing:
        raise FileNotFoundError(
            "Missing nSidis-study input file(s):\n  "
            + "\n  ".join(missing)
        )
    #endif
    if schema_problems:
        raise KeyError(
            "nSidis-study branch preflight failed:\n  "
            + "\n  ".join(schema_problems)
        )
    #endif



def derive_clasdis_epx_valley_cut(
    period: PeriodConfig,
    clasdis_f: Dict[str, np.ndarray],
    *,
    ft_theta_max: float,
    max_probe_energy: float,
    fallback_cut: float,
) -> float:
    """
    Determine the period-specific upper M_X^2(epX) cut from CLASDIS FD events.

    The relevant CLASDIS FD spectrum has a missing-pi0 peak at low M_X^2 and
    a missing-eta peak near m_eta^2.  The desired cut is NOT the eta-peak
    location; it is the local minimum in the trough separating the two peaks.

    To make that definition explicit and robust, we:
      1. build the same minimally-selected FD CLASDIS spectrum used by the
         diagnostic canvas;
      2. identify the pi0 peak in [-0.02, 0.10] GeV^2;
      3. identify the eta peak in [0.24, 0.38] GeV^2;
      4. search ONLY the physical inter-peak trough [0.10, 0.24] GeV^2 and
         take the minimum of the smoothed histogram there.

    The resulting single upper edge is imposed identically in FD and FT.
    ``fallback_cut`` is used only if there are too few CLASDIS events or the
    histogram cannot be formed at all.
    """
    mm2 = np.asarray(clasdis_f["ep_missing_mass2"], dtype=float)
    electron_p = np.asarray(clasdis_f["electron_p"], dtype=float)
    theta_egamma = np.asarray(clasdis_f["theta_egamma_deg"], dtype=float)
    eprobe = np.asarray(clasdis_f["pred_probe_energy"], dtype=float)

    base = (
        np.asarray(clasdis_f["valid_tag"], dtype=bool)
        & stage2_region_mask(clasdis_f, "FD_all", ft_theta_max)
        & np.isfinite(electron_p)
        & (electron_p > 2.0)
        & np.isfinite(theta_egamma)
        & (theta_egamma > THETA_EGAMMA_MIN_DEG)
        & np.isfinite(eprobe)
        & (eprobe >= 0.40)
        & (eprobe < max_probe_energy)
        & np.isfinite(mm2)
        & (mm2 >= -0.10)
        & (mm2 < 0.60)
    )

    values = mm2[base]
    if values.size < 500:
        log(
            f"{period.label}: WARNING: only {values.size:,} CLASDIS FD events "
            "available for the M_X^2(epX) pi0/eta valley determination; "
            f"using fallback upper cut {fallback_cut:.4f} GeV^2."
        )
        return float(fallback_cut)
    #endif

    # 2 MeV^2-wide bins over the displayed range, followed by modest Gaussian
    # smoothing.  The cut is selected from the smoothed spectrum but reported
    # at the corresponding physical histogram-bin center.
    edges = np.linspace(-0.10, 0.60, 351)
    hist, _ = np.histogram(values, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    smooth = gaussian_filter1d(hist.astype(float), sigma=3.0, mode="nearest")

    if not np.any(smooth > 0):
        log(
            f"{period.label}: WARNING: empty CLASDIS FD M_X^2(epX) histogram; "
            f"using fallback upper cut {fallback_cut:.4f} GeV^2."
        )
        return float(fallback_cut)
    #endif

    pi0_window = (centers >= -0.02) & (centers <= 0.10)
    eta_window = (centers >= 0.24) & (centers <= 0.38)
    valley_window = (centers >= 0.10) & (centers <= 0.24)

    pi0_indices = np.flatnonzero(pi0_window)
    eta_indices = np.flatnonzero(eta_window)
    valley_indices = np.flatnonzero(valley_window)

    i_pi0 = int(pi0_indices[np.argmax(smooth[pi0_indices])])
    i_eta = int(eta_indices[np.argmax(smooth[eta_indices])])

    # Restrict the trough search to bins that are simultaneously inside the
    # physical 0.10--0.24 GeV^2 valley window and strictly between the two
    # observed peak locations.  This prevents the eta maximum itself from ever
    # being returned as the cut.
    between = valley_indices[
        (valley_indices > i_pi0) & (valley_indices < i_eta)
    ]
    if between.size == 0:
        # This should be extremely unusual.  Fall back to the minimum in the
        # physical valley window itself rather than to the eta-peak-scale
        # nominal bound.
        between = valley_indices
    #endif

    i_valley = int(between[np.argmin(smooth[between])])
    cut = float(centers[i_valley])

    log(
        f"{period.label}: CLASDIS FD M_X^2(epX): pi0 peak ~= "
        f"{centers[i_pi0]:.4f} GeV^2, eta peak ~= {centers[i_eta]:.4f} GeV^2, "
        f"inter-peak minimum = {cut:.4f} GeV^2. "
        "Using that minimum as the upper M_X^2(epX) cut in BOTH FD and FT."
    )
    return cut


def make_nsidis_shape_comparison_canvases(
    period,
    data_f,
    pi0_f,
    dvcs_f,
    outdir,
    *,
    clasdis_f=None,
    ft_theta_max,
    max_probe_energy,
    mm2_min,
    mm2_max,
    probe_m2_max,
):
    """
    2x8 denominator-shape comparison for FT and FD_all.

    Row 1 is the minimal common support.  Row 2 adds only the accepted
    M_X^2(epX) exclusivity window.  The upper edge is determined period by
    period from the minimum between the missing-pi0 and missing-eta peaks in
    the FD CLASDIS distribution, then imposed identically in FD and FT.

    MM^2(e gamma X), MM^2(ep gamma X), Delta phi, and Angle(gamma,X) remain purely
    diagnostic while the Valerii presentation labels/cut definitions are being
    treated as unreliable.
    """
    specs = (
        (
            "Mx2_ep",
            r"$MM^2(epX)$",
            r"$MM^2(epX)$ (GeV$^2$)",
            -0.20,
            0.60,
            140,
        ),
        (
            "Mx2_egamma",
            r"$MM^2(e\gamma X)$",
            r"$MM^2(e\gamma X)$ (GeV$^2$)",
            1.00,
            5.00,
            140,
        ),
        (
            "Mx2_epgamma",
            r"$MM^2(ep\gamma X)$",
            r"$MM^2(ep\gamma X)$ (GeV$^2$)",
            -0.20,
            0.30,
            100,
        ),
        (
            "Delta_phi_pgamma",
            r"$\Delta\phi(p,\gamma)$",
            r"$\Delta\phi(p,\gamma)$ (deg)",
            -25.0,
            25.0,
            120,
        ),
        (
            "Angle_gammaX",
            r"$\mathrm{Angle}(\gamma,X)$",
            r"$\mathrm{Angle}(\gamma,X)$ (deg)",
            0.0,
            30.0,
            120,
        ),
        (
            "pTmiss",
            r"$p_{T,\rm miss}$",
            r"$p_{T,\rm miss}$ (GeV)",
            0.0,
            1.0,
            100,
        ),
        (
            "xF2",
            r"$x_{F,2}(\gamma)$",
            r"$x_{F,2}(\gamma)$",
            0.00,
            1.20,
            120,
        ),
        (
            "Egamma",
            r"$E_{\gamma,\rm tag}$",
            r"$E_{\gamma,\rm tag}$ (GeV)",
            0.40,
            max_probe_energy,
            100,
        ),
    )

    def vals(feat, key):
        if key == "Mx2_ep":
            return np.asarray(feat["ep_missing_mass2"], dtype=float)
        #endif
        if key == "Mx2_egamma":
            return np.asarray(feat["egamma_missing_mass2"], dtype=float)
        #endif
        if key == "Mx2_epgamma":
            return np.asarray(feat["epgamma_missing_mass2"], dtype=float)
        #endif
        if key == "Delta_phi_pgamma":
            return np.asarray(feat["delta_phi_pgamma_deg"], dtype=float)
        #endif
        if key == "Angle_gammaX":
            return np.asarray(feat["theta_gammaX_deg"], dtype=float)
        #endif
        if key == "pTmiss":
            return np.asarray(feat["stored_pTmiss"], dtype=float)
        #endif
        if key == "xF2":
            return np.asarray(feat["stored_xF2"], dtype=float)
        #endif
        return np.asarray(feat["tag_energy"], dtype=float)

    samples = [
        ("data", data_f, "data", "black"),
        ("pi0", pi0_f, r"exclusive $\pi^0$ (AAO)", "tab:red"),
    ]
    if clasdis_f is not None:
        samples.append(("clasdis", clasdis_f, "SIDIS (CLASDIS)", "tab:orange"))
    #endif
    samples.append(("dvcs", dvcs_f, "BH/DVCS (DVCSgen)", "tab:blue"))

    row_labels = (
        "minimal",
        r"$+\ M_X^2(epX)$ cut",
    )

    for region in ("FT", "FD_all"):
        fig, axes = plt.subplots(2, 8, figsize=(28.0, 8.6))
        log(
            f"{period.label} {region}: denominator support -> "
            f"M_X^2(epX) window [{mm2_min:.4f}, {mm2_max:.4f}) GeV^2"
        )

        for sample_key, feat, label, color in samples:
            minimal = (
                np.asarray(feat["valid_tag"], dtype=bool)
                & stage2_region_mask(feat, region, ft_theta_max)
                & np.isfinite(feat["electron_p"])
                & (feat["electron_p"] > 2.0)
                & np.isfinite(feat["theta_egamma_deg"])
                & (feat["theta_egamma_deg"] > THETA_EGAMMA_MIN_DEG)
                & np.isfinite(feat["pred_probe_energy"])
                & (feat["pred_probe_energy"] >= 0.40)
                & (feat["pred_probe_energy"] < max_probe_energy)
            )
            cut_ep = (
                np.isfinite(feat["ep_missing_mass2"])
                & (feat["ep_missing_mass2"] >= mm2_min)
                & (feat["ep_missing_mass2"] < mm2_max)
            )
            masks = (minimal, minimal & cut_ep)

            counts = [int(np.count_nonzero(m)) for m in masks]
            frac = 100.0 * counts[1] / counts[0] if counts[0] > 0 else 0.0
            log(
                f"  {label}: minimal={counts[0]:,} -> epX cut={counts[1]:,} "
                f"({frac:.2f}% retained)"
            )

            for irow, mask in enumerate(masks):
                for icol, spec in enumerate(specs):
                    key, _, _, lo, hi, nb = spec
                    v = vals(feat, key)
                    use = mask & np.isfinite(v) & (v >= lo) & (v < hi)
                    h, edges = np.histogram(v[use], bins=nb, range=(lo, hi))
                    h = h.astype(float)
                    if np.sum(h) > 0:
                        h /= np.sum(h)
                    #endif
                    centers = 0.5 * (edges[:-1] + edges[1:])
                    axes[irow, icol].step(
                        centers,
                        h,
                        where="mid",
                        linewidth=1.20,
                        color=color,
                        label=label,
                    )
                #endfor
            #endfor
        #endfor

        for irow in range(2):
            for icol, spec in enumerate(specs):
                ax = axes[irow, icol]
                ax.set_xlim(spec[3], spec[4])
                ax.set_xlabel(spec[2])
                ax.grid(alpha=0.16)
                if spec[0] == "Mx2_ep":
                    ax.axvline(mm2_min, linestyle="--", linewidth=0.9, color="0.35")
                    ax.axvline(mm2_max, linestyle="--", linewidth=0.9, color="0.35")
                #endif
                if irow == 0:
                    ax.set_title(spec[1])
                #endif
            #endfor
            axes[irow, 0].set_ylabel(row_labels[irow] + "\nunit-normalized")
        #endfor

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.947),
                ncol=min(4, len(handles)),
                frameon=False,
                fontsize=8.5,
            )
        #endif
        fig.suptitle(
            f"{period.label}: {region} denominator shape comparison\n"
            r"minimal: $p_e>2$ GeV, $\theta_{e\gamma}>5^\circ$; "
            rf"nominal adds ${mm2_min:.3f}<MM^2(epX)<{mm2_max:.3f}$ GeV$^2$ "
            r"(upper edge = FD CLASDIS inter-peak minimum between $\pi^0$ and $\eta$; same cut in FT and FD)",
            fontsize=11.2,
        )
        safe_finalize_figure(
            fig,
            outdir / f"canvas_shape_comparison_{period.key}_{region.lower()}.png",
            rect=(0, 0, 1, 0.91),
        )
        plt.close(fig)
    #endfor



def make_nsidis_three_variable_fit_canvases(
    period,
    rows,
    data_f,
    pi0_f,
    dvcs_f,
    outdir,
    *,
    ft_theta_max,
    mm2_min,
    mm2_max,
    probe_m2_max,
):
    """
    Explicit nominal fit display.

    Columns are successful predicted-probe energy bins; rows are
    Delta_phi, pTmiss, and xF2. Failed bins are omitted.
    """
    for region in ("FT", "FD_all"):
        rr = sorted(
            [
                row
                for row in rows
                if (
                    row["region"] == region
                    and int(row.get("fit_success", 0)) == 1
                )
            ],
            key=lambda row: float(row["energy_low_GeV"]),
        )
        if not rr:
            continue
        #endif

        fig, axes = plt.subplots(
            3,
            len(rr),
            figsize=(max(15.0, 3.55 * len(rr)), 11.2),
            squeeze=False,
        )

        for icol, row in enumerate(rr):
            elo = float(row["energy_low_GeV"])
            ehi = float(row["energy_high_GeV"])

            masks = {}
            for name, feat in (
                ("data", data_f),
                ("pi0", pi0_f),
                ("dvcs", dvcs_f),
            ):
                mask = stage2_fit_mask(
                    feat,
                    region,
                    ft_theta_max,
                    elo,
                    ehi,
                    mm2_min,
                    mm2_max,
                    probe_m2_max,
                )
                mask &= (
                    np.isfinite(feat["electron_p"])
                    & (feat["electron_p"] > 2.0)
                )
                for var in NOMINAL_THREE_VARIABLES:
                    lo, hi, _ = NOMINAL_THREE_RANGES[var]
                    values = _nominal_three_values(
                        feat,
                        var,
                    )
                    if var == "xF2":
                        # Match the production fit: no xF2 selection cut.
                        mask &= np.isfinite(values)
                    else:
                        mask &= (
                            np.isfinite(values)
                            & (values >= lo)
                            & (values < hi)
                        )
                    #endif
                #endfor
                masks[name] = mask
            #endfor

            fpi0 = float(row["pi0_fraction"])

            for irow, var in enumerate(
                NOMINAL_THREE_VARIABLES
            ):
                lo, hi, nbins = NOMINAL_THREE_RANGES[var]
                edges = np.linspace(lo, hi, nbins + 1)
                centers = 0.5 * (
                    edges[:-1] + edges[1:]
                )

                hist = {}
                for name, feat in (
                    ("data", data_f),
                    ("pi0", pi0_f),
                    ("dvcs", dvcs_f),
                ):
                    h, _ = np.histogram(
                        _nominal_three_hist_values(
                            feat,
                            var,
                        )[masks[name]],
                        bins=edges,
                    )
                    hist[name] = h.astype(float)
                #endfor

                hp = hist["pi0"][None, :]
                hv = hist["dvcs"][None, :]

                pshift = float(
                    row.get(
                        f"{var}_pi0_shift_bins", 0.0
                    )
                )
                psigma = float(
                    row.get(
                        f"{var}_pi0_sigma_bins", 0.0
                    )
                )
                dshift = float(
                    row.get(
                        f"{var}_dvcs_shift_bins", 0.0
                    )
                )
                dsigma = float(
                    row.get(
                        f"{var}_dvcs_sigma_bins", 0.0
                    )
                )

                raw_pi0 = normalized_template(
                    hp
                ).reshape(hp.shape).ravel()
                raw_dvcs = normalized_template(
                    hv
                ).reshape(hv.shape).ravel()
                morph_pi0 = morph_template_second_axis(
                    hp, pshift, psigma
                ).ravel()
                morph_dvcs = morph_template_second_axis(
                    hv, dshift, dsigma
                ).ravel()

                ndata = float(np.sum(hist["data"]))
                raw_pi0 *= ndata * fpi0
                raw_dvcs *= ndata * (1.0 - fpi0)
                morph_pi0 *= ndata * fpi0
                morph_dvcs *= (
                    ndata * (1.0 - fpi0)
                )
                total = morph_pi0 + morph_dvcs

                ax = axes[irow, icol]
                data_counts = np.asarray(
                    hist["data"],
                    dtype=float,
                )
                populated = data_counts > 0.0

                ax.errorbar(
                    centers[populated],
                    data_counts[populated],
                    yerr=np.sqrt(
                        data_counts[populated]
                    ),
                    fmt="o",
                    ms=2.0,
                    linewidth=0.5,
                    color="black",
                    label="data",
                    zorder=6,
                )
                ax.step(
                    centers,
                    raw_pi0,
                    where="mid",
                    color="tab:red",
                    linestyle="--",
                    alpha=0.5,
                    linewidth=0.8,
                    label=r"$\pi^0$ pre-morph",
                )
                ax.step(
                    centers,
                    raw_dvcs,
                    where="mid",
                    color="tab:blue",
                    linestyle="--",
                    alpha=0.5,
                    linewidth=0.8,
                    label="DVCS pre-morph",
                )
                ax.step(
                    centers,
                    morph_pi0,
                    where="mid",
                    color="tab:red",
                    linewidth=1.2,
                    label=r"$\pi^0$ morphed",
                )
                ax.step(
                    centers,
                    morph_dvcs,
                    where="mid",
                    color="tab:blue",
                    linewidth=1.2,
                    label="DVCS morphed",
                )
                ax.step(
                    centers,
                    total,
                    where="mid",
                    color="tab:green",
                    linewidth=1.5,
                    label="total fit",
                )

                # Linear y scale is substantially clearer for these
                # component fits.  It also avoids log-axis pathologies in
                # sparsely populated tails.
                ax.set_yscale("linear")
                ax.set_xlim(lo, hi)

                ymax = max(
                    float(np.max(data_counts))
                    if data_counts.size
                    else 0.0,
                    float(np.max(total))
                    if total.size
                    else 0.0,
                    float(np.max(morph_pi0))
                    if morph_pi0.size
                    else 0.0,
                    float(np.max(morph_dvcs))
                    if morph_dvcs.size
                    else 0.0,
                )
                ax.set_ylim(
                    0.0,
                    1.12 * ymax
                    if ymax > 0.0
                    else 1.0,
                )
                ax.grid(alpha=0.16)

                xlabel = (
                    r"$\Delta\phi$ (rad)"
                    if var == "Delta_phi"
                    else (
                        r"$p_{T,\rm miss}$ (GeV)"
                        if var == "pTmiss"
                        else r"$x_{F,2}(\gamma)$"
                    )
                )
                ax.set_xlabel(xlabel)

                if irow == 0:
                    ax.set_title(
                        f"{elo:.2f}-{ehi:.2f} GeV\n"
                        rf"$f_{{\pi^0}}^{{shared}}="
                        rf"{fpi0:.3f}$",
                        fontsize=8.8,
                    )
                #endif

                single = float(
                    row.get(
                        f"pi0_fraction_single_{var}",
                        np.nan,
                    )
                )
                usable = int(
                    row.get(
                        f"pi0_fraction_single_{var}_usable",
                        1,
                    )
                )
                boundary_text = (
                    ""
                    if usable
                    else " (boundary; excluded)"
                )
                ax.text(
                    0.02,
                    0.04,
                    (
                        rf"single $f_{{\pi^0}}="
                        rf"{single:.3f}$"
                        + boundary_text
                        + "\n"
                        + rf"$\pi^0$: shift={pshift:.2f}, "
                        rf"$\sigma$={psigma:.2f} bins"
                        + "\n"
                        + rf"DVCS: shift={dshift:.2f}, "
                        rf"$\sigma$={dsigma:.2f} bins"
                    ),
                    transform=ax.transAxes,
                    fontsize=6.5,
                    va="bottom",
                    ha="left",
                )
            #endfor
        #endfor

        handles, labels = (
            axes[0, 0].get_legend_handles_labels()
        )
        fig.legend(
            handles,
            labels,
            loc="upper center",
            ncol=6,
            fontsize=7.5,
            frameon=False,
            bbox_to_anchor=(0.5, 0.915),
        )
        fig.suptitle(
            f"{period.label}: {region} nominal denominator fits\n"
            r"shared $f_{\pi^0}$ from "
            r"$\Delta\phi$, $p_{T,\rm miss}$, "
            r"$x_{F,2}(\gamma)$; "
            rf"${mm2_min:.2f}<M_X^2(ep)<"
            rf"{mm2_max:.2f}$ GeV$^2$, "
            rf"$|M_X^2(ep\gamma)|<"
            rf"{probe_m2_max:.2f}$ GeV$^2$, "
            r"$p_e>2$ GeV, "
            r"$\theta_{e\gamma}>5^\circ$; "
            r"$x_{F,2}$ is not selection-cut",
            fontsize=11.0,
            y=0.997,
        )
        safe_finalize_figure(
            fig,
            outdir
            / (
                "canvas_nominal_three_variable_fits_"
                f"{period.key}_{region.lower()}.png"
            ),
            rect=(0.0, 0.0, 1.0, 0.84),
        )
        plt.close(fig)
    #endfor



def make_nsidis_pi0_fraction_energy_canvas(
    period: PeriodConfig,
    rows: Sequence[Dict[str, object]],
    outdir: Path,
) -> None:
    """Show shared and single-variable fitted pi0 fractions versus energy."""
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.6), squeeze=False)
    axes = axes[0]

    for ax, region in zip(axes, ("FT", "FD_all")):
        rr = sorted(
            [
                row for row in rows
                if row.get("region") == region
                and int(row.get("fit_success", 0)) == 1
            ],
            key=lambda row: float(row["energy_low_GeV"]),
        )
        if not rr:
            ax.set_axis_off()
            continue
        #endif

        x = np.asarray([row_energy_coordinate(row) for row in rr], dtype=float)
        shared = np.asarray([float(row["pi0_fraction"]) for row in rr], dtype=float)
        ax.plot(
            x, shared, marker="o", linewidth=1.8, color="black",
            label="shared 3-variable fit", zorder=5,
        )

        ptcat = np.asarray(
            [
                float(row.get("pt_category_pi0_fraction", np.nan))
                for row in rr
            ],
            dtype=float,
        )
        ptcat_lo = np.asarray(
            [
                float(row.get("pt_category_pi0_fraction_err_low", np.nan))
                for row in rr
            ],
            dtype=float,
        )
        ptcat_hi = np.asarray(
            [
                float(row.get("pt_category_pi0_fraction_err_high", np.nan))
                for row in rr
            ],
            dtype=float,
        )
        ptcat_ok = np.asarray(
            [
                int(row.get("pt_category_fit_success", 0)) == 1
                for row in rr
            ],
            dtype=bool,
        )
        good_ptcat = np.isfinite(ptcat) & ptcat_ok
        if np.any(good_ptcat):
            ax.errorbar(
                x[good_ptcat],
                ptcat[good_ptcat],
                yerr=np.vstack(
                    (
                        ptcat_lo[good_ptcat],
                        ptcat_hi[good_ptcat],
                    )
                ),
                marker="P",
                linestyle="--",
                linewidth=1.2,
                capsize=2,
                label=r"two-category $p_T$ fit",
                zorder=4,
            )
        #endif

        specs = (
            ("Delta_phi", r"$\Delta\phi$ only", "tab:blue", "s"),
            ("pTmiss", r"$p_{T,\rm miss}$ only", "tab:orange", "^"),
            ("xF2", r"$x_{F,2}(\gamma)$ only", "tab:green", "D"),
        )
        for var, label, color, marker in specs:
            y = np.asarray(
                [float(row.get(f"pi0_fraction_single_{var}", np.nan)) for row in rr],
                dtype=float,
            )
            usable = np.asarray(
                [
                    int(row.get(f"pi0_fraction_single_{var}_usable", 1)) == 1
                    for row in rr
                ],
                dtype=bool,
            )
            good = np.isfinite(x) & np.isfinite(y) & usable
            if np.any(good):
                ax.plot(
                    x[good], y[good], marker=marker, linewidth=1.0,
                    color=color, label=label,
                )
            #endif
            bad = np.isfinite(x) & np.isfinite(y) & (~usable)
            if np.any(bad):
                ax.plot(
                    x[bad], y[bad], linestyle="none", marker=marker,
                    markersize=6.0, markerfacecolor="none",
                    markeredgecolor=color,
                )
            #endif
        #endfor

        ax.set_ylim(0.0, 1.05)
        ax.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        ax.set_ylabel(r"$f_{\pi^0}$")
        ax.set_title("FT" if region == "FT" else "FD all")
        ax.grid(alpha=0.18)
    #endfor

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.91),
        ncol=4, frameon=False, fontsize=8.5,
    )
    fig.suptitle(
        f"{period.label}: Stage-2 fitted $\\pi^0$ fraction versus energy",
        fontsize=13.0, y=0.985,
    )
    safe_finalize_figure(
        fig,
        outdir / f"canvas_pi0_fraction_vs_energy_{period.key}.png",
        rect=(0.0, 0.0, 1.0, 0.82),
    )
    plt.close(fig)





def make_nsidis_pt_category_canvas(
    period: PeriodConfig,
    rows: Sequence[Dict[str, object]],
    outdir: Path,
) -> None:
    """Show the two-category SIDIS-pT composition cross-check."""
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 9.2), squeeze=False)

    for icol, region in enumerate(("FT", "FD_all")):
        rr = sorted(
            [row for row in rows if row.get("region") == region],
            key=lambda row: float(row["energy_low_GeV"]),
        )
        if not rr:
            axes[0, icol].set_axis_off()
            axes[1, icol].set_axis_off()
            continue
        #endif

        x = np.asarray([row_energy_coordinate(row) for row in rr], dtype=float)
        shared = np.asarray(
            [float(row.get("pi0_fraction", np.nan)) for row in rr],
            dtype=float,
        )
        ptcat = np.asarray(
            [float(row.get("pt_category_pi0_fraction", np.nan)) for row in rr],
            dtype=float,
        )
        ptcat_lo = np.asarray(
            [float(row.get("pt_category_pi0_fraction_err_low", np.nan)) for row in rr],
            dtype=float,
        )
        ptcat_hi = np.asarray(
            [float(row.get("pt_category_pi0_fraction_err_high", np.nan)) for row in rr],
            dtype=float,
        )
        ok = np.asarray(
            [int(row.get("pt_category_fit_success", 0)) == 1 for row in rr],
            dtype=bool,
        )

        top = axes[0, icol]
        top.plot(
            x, shared,
            marker="o", linewidth=1.4, color="black",
            label="shared 3-variable fit",
        )
        good = ok & np.isfinite(ptcat)
        if np.any(good):
            top.errorbar(
                x[good],
                ptcat[good],
                yerr=np.vstack((ptcat_lo[good], ptcat_hi[good])),
                marker="P",
                linestyle="--",
                linewidth=1.2,
                capsize=2,
                label=r"two-category $p_T$ fit",
            )
        #endif
        top.set_ylim(0.0, 1.05)
        top.set_ylabel(r"$f_{\pi^0}$")
        top.set_title("FT" if region == "FT" else "FD all")
        top.grid(alpha=0.18)
        top.legend(fontsize=8, frameon=False)

        bottom = axes[1, icol]
        pi0_high = np.asarray(
            [float(row.get("pt_category_pi0_high_eff", np.nan)) for row in rr],
            dtype=float,
        )
        dvcs_high = np.asarray(
            [float(row.get("pt_category_dvcs_high_eff", np.nan)) for row in rr],
            dtype=float,
        )
        data_high = np.asarray(
            [float(row.get("pt_category_data_high_fraction", np.nan)) for row in rr],
            dtype=float,
        )

        bottom.plot(x, pi0_high, marker="o", linewidth=1.2, label=r"$\pi^0$ MC")
        bottom.plot(x, dvcs_high, marker="s", linewidth=1.2, label="BH/DVCS MC")
        bottom.plot(
            x, data_high, marker="^", linewidth=1.2, color="black", label="data"
        )
        bottom.set_ylim(0.0, 1.05)
        bottom.set_xlabel(r"$E_{\mathrm{probe}}^{\mathrm{pred}}$ (GeV)")
        bottom.set_ylabel(r"fraction with $p_T \geq p_T^{\rm split}$")
        bottom.grid(alpha=0.18)
        bottom.legend(fontsize=8, frameon=False)
    #endfor

    split_values = {
        float(row.get("pt_category_split_GeV", np.nan))
        for row in rows
        if np.isfinite(float(row.get("pt_category_split_GeV", np.nan)))
    }
    split = next(iter(split_values)) if len(split_values) == 1 else float("nan")

    fig.suptitle(
        f"{period.label}: two-category SIDIS-$p_T$ composition test\n"
        + (
            rf"$p_T^{{split}}={split:.3g}$ GeV; "
            if np.isfinite(split)
            else ""
        )
        + "finite MC statistics included in likelihood",
        fontsize=12.5,
        y=0.985,
    )
    safe_finalize_figure(
        fig,
        outdir / f"canvas_pt_category_composition_{period.key}.png",
        rect=(0.0, 0.0, 1.0, 0.93),
    )
    plt.close(fig)




def run_grand_stage1_diagnostics_for_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
    production_root: Path,
) -> None:
    """
    Produce the broad all-variable diagnostic suite for a normal production run.

    This deliberately uses an independent, bounded epgamma read.  The grand
    diagnostics are visual cross-checks only; keeping their I/O separate means
    they cannot alter the production event population or force dozens of
    extra branches to remain resident during Stage 2/3.

    By default only the first 4M entries of each data/template tree are used.
    Override with --grand-diagnostics-max-entries; 0 means all entries.
    """
    if bool(args_dict.get("skip_grand_diagnostics", False)):
        return
    #endif

    diag_max_entries = int(
        args_dict.get(
            "grand_diagnostics_max_entries",
            4000000,
        )
    )
    tree_name = args_dict.get("tree")
    angle_mode = str(args_dict["angles"])
    tag_min = float(args_dict["tag_min"])
    tag_max = float(args_dict["tag_max"])
    ft_theta_max = float(args_dict["ft_theta_max"])
    max_probe_energy = float(
        args_dict["nsidis_pilot_energy_max"]
    )
    mm2_min = float(args_dict["den_fit_mm2_min"])
    mm2_max = float(args_dict["den_fit_mm2_max"])
    probe_m2_max = float(
        args_dict["nsidis_pilot_probe_m2_max"]
    )

    ns_epg_path, _ = resolve_nsidis_paths(
        period,
        args_dict,
    )
    outdir = (
        production_root
        / "stage1_grand_diagnostics"
        / period.key
    )
    outdir.mkdir(parents=True, exist_ok=True)

    grand_optional = tuple(
        dict.fromkeys(
            grand_stage1_optional_branches()
            + EPG_OPTIONAL_PI0_MC
            + EPG_OPTIONAL_DVCS_MC
            + EPG_OPTIONAL_DATA
        )
    )

    input_specs = (
        ("data", ns_epg_path),
        ("pi0", period.epgamma_pi0_mc),
        ("dvcs", period.epgamma_dvcs_mc),
    )

    arrays_by_sample = {}
    samples_by_name = {}
    features_by_sample = {}

    log(
        f"{period.label}: producing grand Stage-1 diagnostics "
        f"(max entries/sample = "
        f"{'all' if diag_max_entries == 0 else f'{diag_max_entries:,}'})."
    )

    for sample_name, path in input_specs:
        arrays, found_tree, total = read_branches(
            path,
            EPG_REQUIRED,
            grand_optional,
            tree_name,
            diag_max_entries,
        )
        sample = extract_epgamma(
            arrays,
            angle_mode,
        )
        feat = build_epgamma_denominator_features(
            period,
            sample,
            tag_min,
            tag_max,
            ft_theta_max,
        )
        feat["electron_p"] = np.asarray(
            electron_momentum_from_p3(
                sample.electron_p3
            ),
            dtype=np.float32,
        )

        arrays_by_sample[sample_name] = arrays
        samples_by_name[sample_name] = sample
        features_by_sample[sample_name] = feat

        log(
            f"{period.label}: grand diagnostics {sample_name} "
            f"loaded {len(sample.tag_energy):,}/{total:,} entries "
            f"from tree '{found_tree}'."
        )
    #endfor

    photon_acceptance = infer_photon_angular_acceptance(
        samples_by_name["data"],
        source=(
            f"{period.label} nSidis epgamma data; "
            "grand Stage-1 diagnostics"
        ),
    )
    for feat in features_by_sample.values():
        attach_photon_angular_acceptance(
            feat,
            photon_acceptance,
        )
    #endfor

    make_grand_stage1_diagnostics(
        period,
        arrays_by_sample,
        samples_by_name,
        features_by_sample,
        outdir,
        ft_theta_max=ft_theta_max,
        max_probe_energy=max_probe_energy,
        mm2_min=mm2_min,
        mm2_max=mm2_max,
        probe_m2_max=probe_m2_max,
    )

    # Drop the large diagnostic-only arrays immediately.
    arrays_by_sample.clear()
    samples_by_name.clear()
    features_by_sample.clear()

    log(
        f"{period.label}: grand Stage-1 diagnostics complete: {outdir}"
    )



def process_nsidis_stage1_only_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
    output_root: str,
) -> Dict[str, object]:
    """
    Fast Stage-1 shape-comparison path.

    The normal --stage1-only mode reads only the branches required for the
    standard data/AAO/CLASDIS/DVCS shape comparison. The broad all-variable
    grand diagnostics are opt-in via --stage1-grand-diagnostics.
    """
    t0 = time.perf_counter()
    production_root = Path(output_root)
    stage1_outdir = production_root / "stage1_shape_comparison"
    stage1_outdir.mkdir(parents=True, exist_ok=True)

    max_entries = int(args_dict["max_entries"])
    tree_name = args_dict.get("tree")
    angle_mode = str(args_dict["angles"])
    tag_min = float(args_dict["tag_min"])
    tag_max = float(args_dict["tag_max"])
    ft_theta_max = float(args_dict["ft_theta_max"])
    max_probe_energy = float(args_dict["nsidis_pilot_energy_max"])
    mm2_min = float(args_dict["den_fit_mm2_min"])
    mm2_max = float(args_dict["den_fit_mm2_max"])
    probe_m2_max = float(args_dict["nsidis_pilot_probe_m2_max"])

    ns_epg_path, _ = resolve_nsidis_paths(period, args_dict)

    input_specs = [
        ("data", ns_epg_path, EPG_OPTIONAL_DATA),
        ("pi0", period.epgamma_pi0_mc, EPG_OPTIONAL_PI0_MC),
    ]

    if period.clasdis_epgamma_mc is not None:
        clasdis_path = Path(period.clasdis_epgamma_mc)
        if clasdis_path.exists():
            input_specs.append(
                (
                    "clasdis",
                    str(clasdis_path),
                    EPG_OPTIONAL_CLASDIS_MC,
                )
            )
        else:
            log(
                f"{period.label}: CLASDIS epgammaX configured but not present; "
                f"omitting it from the shape comparison: {clasdis_path}"
            )
        #endif
    #endif

    input_specs.append(
        ("dvcs", period.epgamma_dvcs_mc, EPG_OPTIONAL_DVCS_MC)
    )

    samples_by_name: Dict[str, EPGammaSample] = {}
    features_by_sample: Dict[str, Dict[str, np.ndarray]] = {}

    for sample_name, path, base_optional in input_specs:
        progress_label = f"{period.label}: Stage-1 {sample_name}"
        log(f"{period.label}: Stage-1-only reading {sample_name} epgamma.")

        arrays, found_tree, total = read_branches_with_progress(
            path,
            EPG_REQUIRED,
            base_optional,
            tree_name,
            max_entries,
            progress_label=progress_label,
        )

        sample = extract_epgamma(arrays, angle_mode)
        feat = build_epgamma_denominator_features(
            period,
            sample,
            tag_min,
            tag_max,
            ft_theta_max,
        )
        feat["electron_p"] = np.asarray(
            electron_momentum_from_p3(sample.electron_p3),
            dtype=np.float32,
        )

        # Raw arrays are not needed by the standard shape comparison.
        del arrays

        samples_by_name[sample_name] = sample
        features_by_sample[sample_name] = feat

        log(
            f"{period.label}: {sample_name} tree '{found_tree}', "
            f"loaded {len(sample.tag_energy):,}/{total:,} entries; "
            "feature construction complete."
        )
    #endfor

    photon_acceptance = infer_photon_angular_acceptance(
        samples_by_name["data"],
        source=f"{period.label} nSidis epgamma data; Stage-1-only",
    )
    for feat in features_by_sample.values():
        attach_photon_angular_acceptance(feat, photon_acceptance)
    #endfor

    log(
        f"{period.label}: Stage-1-only fixed photon angular acceptance: "
        "FT 2.4-5.0 deg; FD 6.0-35.0 deg."
    )

    # IMPORTANT: Stage-1-only used to pass the command-line fallback
    # --den-fit-mm2-max straight into the diagnostic canvas.  That meant the
    # CLASDIS FD pi0/eta-valley finder was never called in --stage1-only mode,
    # so the plots misleadingly continued to show the old 0.309 GeV^2 edge.
    # Derive the period-specific upper edge here as well, before making either
    # the FT or FD canvas.  The same FD-derived number is then imposed in both
    # detector regions.
    fallback_mm2_max = float(mm2_max)
    clasdis_feat = features_by_sample.get("clasdis")
    if clasdis_feat is not None:
        mm2_max = derive_clasdis_epx_valley_cut(
            period,
            clasdis_feat,
            ft_theta_max=ft_theta_max,
            max_probe_energy=max_probe_energy,
            fallback_cut=fallback_mm2_max,
        )
    else:
        mm2_max = fallback_mm2_max
        log(
            f"{period.label}: WARNING: Stage-1-only CLASDIS sample is "
            "unavailable; using fallback M_X^2(epX) upper edge "
            f"{mm2_max:.4f} GeV^2 in both FD and FT."
        )
    #endif

    log(
        f"{period.label}: Stage-1-only FINAL M_X^2(epX) window = "
        f"[{mm2_min:.4f}, {mm2_max:.4f}) GeV^2; upper edge derived from "
        "FD CLASDIS and applied identically to FD and FT."
    )

    make_nsidis_shape_comparison_canvases(
        period,
        features_by_sample["data"],
        features_by_sample["pi0"],
        features_by_sample["dvcs"],
        stage1_outdir,
        clasdis_f=features_by_sample.get("clasdis"),
        ft_theta_max=ft_theta_max,
        max_probe_energy=max_probe_energy,
        mm2_min=mm2_min,
        mm2_max=mm2_max,
        probe_m2_max=probe_m2_max,
    )

    samples_by_name.clear()
    features_by_sample.clear()

    if bool(args_dict.get("stage1_grand_diagnostics", False)):
        log(
            f"{period.label}: starting optional grand Stage-1 diagnostics."
        )
        run_grand_stage1_diagnostics_for_period(
            period,
            args_dict,
            production_root,
        )
    else:
        log(
            f"{period.label}: grand Stage-1 diagnostics skipped "
            "(add --stage1-grand-diagnostics to enable them)."
        )
    #endif

    elapsed = float(time.perf_counter() - t0)
    log(
        f"{period.label}: Stage-1-only diagnostics complete "
        f"in {elapsed:.1f} s."
    )
    return {
        "period": period.key,
        "label": period.label,
        "stage1_only": True,
        "wall_time_s": elapsed,
        "status": "stage1_only_complete",
    }


def process_nsidis_study_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
    output_root: str,
) -> Dict[str, object]:
    """
    Default nSidis production path.

    The optional pilot now forms the first preliminary data/aaogen photon-efficiency scale factor after validating denominator composition and association-first numerator purity. The fitted-pi0 denominator uncertainty is not yet propagated.
    """
    if bool(args_dict.get("stage1_only", False)):
        return process_nsidis_stage1_only_period(
            period,
            args_dict,
            output_root,
        )
    #endif

    t0 = time.perf_counter()
    production_root = Path(output_root)
    stage1_outdir = (
        production_root / "stage1_shape_comparison"
    )
    stage2_outdir = (
        production_root / "stage2_composition_fits"
    )
    stage3_outdir = (
        production_root / "stage3_efficiency" / period.key
    )
    for directory in (
        stage1_outdir,
        stage2_outdir,
        stage3_outdir,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    #endfor
    outdir = stage3_outdir

    max_entries = int(args_dict["max_entries"])
    tree_name = args_dict.get("tree")
    angle_mode = str(args_dict["angles"])
    tag_min = float(args_dict["tag_min"])
    tag_max = float(args_dict["tag_max"])
    ft_theta_max = float(args_dict["ft_theta_max"])

    ns_epg_path, ns_epi_path = resolve_nsidis_paths(
        period, args_dict
    )

    loaded = {}
    load_specs = (
        (
            "wagon_epgamma",
            period.epgamma_data,
            EPG_REQUIRED,
            EPG_OPTIONAL_DATA,
            extract_epgamma,
        ),
        (
            "nsidis_epgamma",
            ns_epg_path,
            EPG_REQUIRED,
            EPG_OPTIONAL_DATA,
            extract_epgamma,
        ),
        (
            "wagon_eppi0",
            period.eppi0_data,
            EPPIO_REQUIRED,
            EPPIO_OPTIONAL_DATA,
            extract_eppi0,
        ),
        (
            "nsidis_eppi0",
            ns_epi_path,
            EPPIO_REQUIRED,
            EPPIO_OPTIONAL_DATA,
            extract_eppi0,
        ),
    )

    # The current production configuration intentionally uses the original
    # eppi0 file for both labels.  Reuse the already extracted object when
    # paths are identical instead of paying ROOT I/O + vector construction
    # twice.
    loaded_by_path = {}
    for key, path, required, optional, extractor in load_specs:
        canonical_path = str(Path(path).expanduser().resolve())
        cache_key = (
            canonical_path,
            tuple(required),
            tuple(optional),
            extractor.__name__,
            str(angle_mode),
            int(max_entries),
        )
        if cache_key in loaded_by_path:
            loaded[key] = loaded_by_path[cache_key]
            log(
                f"{period.label}: reusing {key} from already-loaded identical "
                f"ROOT source ({Path(path).name})."
            )
            continue
        #endif

        read_t0 = time.perf_counter()
        log(f"{period.label}: reading {key}.")
        arrays, found_tree, total = read_branches(
            path,
            required,
            optional,
            tree_name,
            max_entries,
        )
        sample = extractor(arrays, angle_mode)
        del arrays
        loaded[key] = sample
        loaded_by_path[cache_key] = sample
        log(
            f"{period.label}: {key} tree '{found_tree}', "
            f"loaded {len(sample.electron_p3):,}/{total:,} entries "
            f"in {time.perf_counter() - read_t0:.1f} s."
        )
    #endfor

    wagon_epg = loaded["wagon_epgamma"]
    ns_epg = loaded["nsidis_epgamma"]
    wagon_epi = loaded["wagon_eppi0"]
    ns_epi = loaded["nsidis_eppi0"]

    # Match the upstream nSidis parent support exactly for overlap checks.
    # Compute these large vector norms once.  They are reused both by the
    # parent masks and the feature stores below.
    wagon_epg_electron_p = np.asarray(
        electron_momentum_from_p3(wagon_epg.electron_p3),
        dtype=np.float32,
    )
    ns_epg_electron_p = np.asarray(
        electron_momentum_from_p3(ns_epg.electron_p3),
        dtype=np.float32,
    )
    wagon_epi_electron_p = np.asarray(
        electron_momentum_from_p3(wagon_epi.electron_p3),
        dtype=np.float32,
    )
    if ns_epi is wagon_epi:
        ns_epi_electron_p = wagon_epi_electron_p
    else:
        ns_epi_electron_p = np.asarray(
            electron_momentum_from_p3(ns_epi.electron_p3),
            dtype=np.float32,
        )
    #endif

    wagon_epg_parent = wagon_epg_electron_p > 2.0
    ns_epg_parent = ns_epg_electron_p > 2.0
    wagon_epi_parent = wagon_epi_electron_p > 2.0
    ns_epi_parent = ns_epi_electron_p > 2.0

    # Use the enormous loose nSidis epgamma sample itself to measure the
    # reconstructed-photon angular support. This is the geometry the predicted
    # probe must fall inside to count as potentially reconstructable.
    feature_t0 = time.perf_counter()
    photon_acceptance = infer_photon_angular_acceptance(
        ns_epg,
        source=f"{period.label} nSidis epgamma data",
        parent_mask=ns_epg_parent,
    )
    log(
        f"{period.label}: fixed photon angular acceptance: "
        "FT 2.4-5.0 deg; FD 6.0-35.0 deg."
    )

    epg_overlap, epg_overlap_summary = event_overlap_masks(
        wagon_epg.raw,
        ns_epg.raw,
        reference_mask=wagon_epg_parent,
        candidate_mask=ns_epg_parent,
    )
    epi_overlap, epi_overlap_summary = event_overlap_masks(
        wagon_epi.raw,
        ns_epi.raw,
        reference_mask=wagon_epi_parent,
        candidate_mask=ns_epi_parent,
    )

    wagon_epg_f = build_epgamma_denominator_features(
        period, wagon_epg, tag_min, tag_max, ft_theta_max
    )
    ns_epg_f = build_epgamma_denominator_features(
        period, ns_epg, tag_min, tag_max, ft_theta_max
    )
    for _feat in (wagon_epg_f, ns_epg_f):
        attach_photon_angular_acceptance(_feat, photon_acceptance)
    #endfor
    wagon_epg_f["electron_p"] = wagon_epg_electron_p
    ns_epg_f["electron_p"] = ns_epg_electron_p

    wagon_epi_f = build_eppi0_exclusivity_features(
        period, wagon_epi
    )
    if ns_epi is wagon_epi:
        ns_epi_f = wagon_epi_f
    else:
        ns_epi_f = build_eppi0_exclusivity_features(
            period, ns_epi
        )
    #endif

    # Retire the old global eppi0 "95% wagon-core" optimization.
    # The full old eppi0 wagon contains long exclusivity tails, so maximizing
    # retention of that entire tree is not the correct numerator-purity target.
    # Numerator exclusivity is evaluated below AFTER exact tag/probe association.
    log(
        f"{period.label}: initial data feature construction/overlap setup "
        f"completed in {time.perf_counter() - feature_t0:.1f} s."
    )

    pilot_wagon_rows: List[Dict[str, object]] = []
    pilot_nsidis_rows: List[Dict[str, object]] = []
    log(
        f"{period.label}: compact production diagnostics enabled: "
        "fit projections, data/MC scale factor, and post-association "
        "eppi0 exclusivity only."
    )

    if (
        bool(args_dict.get("nsidis_pilot_fit", False))
        or bool(args_dict.get("nsidis_driver_study", False))
    ):
        if bool(args_dict.get("nsidis_driver_study", False)):
            log(
                f"{period.label}: dedicated nSidis denominator-driver study "
                "enabled; reading aaogen and dvcsgen epgamma templates."
            )
        else:
            log(
                f"{period.label}: nSidis production-model pilot enabled; "
                "reading aaogen and dvcsgen epgamma templates."
            )
        #endif
        pi_arrays, _, _ = read_branches(
            period.epgamma_pi0_mc,
            EPG_REQUIRED,
            EPG_OPTIONAL_PI0_MC,
            tree_name,
            max_entries,
        )
        pi_epg = extract_epgamma(pi_arrays, angle_mode)
        del pi_arrays

        dv_arrays, _, _ = read_branches(
            period.epgamma_dvcs_mc,
            EPG_REQUIRED,
            EPG_OPTIONAL_DVCS_MC,
            tree_name,
            max_entries,
        )
        dv_epg = extract_epgamma(dv_arrays, angle_mode)
        del dv_arrays

        pi0_f = build_epgamma_denominator_features(
            period, pi_epg, tag_min, tag_max, ft_theta_max
        )
        dvcs_f = build_epgamma_denominator_features(
            period, dv_epg, tag_min, tag_max, ft_theta_max
        )
        for _feat in (pi0_f, dvcs_f):
            attach_photon_angular_acceptance(_feat, photon_acceptance)
        #endfor
        pi0_f["electron_p"] = np.asarray(
            electron_momentum_from_p3(pi_epg.electron_p3),
            dtype=np.float32,
        )
        dvcs_f["electron_p"] = np.asarray(
            electron_momentum_from_p3(dv_epg.electron_p3),
            dtype=np.float32,
        )

        # Production/refactor invariant: all quantities needed by the restored
        # shape-comparison and nominal 3-variable fit canvases must exist before
        # any downstream plotting/fitting starts.
        required_production_features = (
            "tag_energy",
            "pred_probe_energy",
            "pred_probe_mass2",
            "stored_delta_phi_rad",
            "stored_pTmiss",
            "stored_Emiss2",
            "electron_p",
            "valid_tag",
        )
        for sample_name, feat in (
            ("nSidis data", ns_epg_f),
            ("aaogen epgamma", pi0_f),
            ("dvcsgen epgamma", dvcs_f),
        ):
            missing = [
                key
                for key in required_production_features
                if key not in feat
            ]
            if missing:
                raise RuntimeError(
                    f"{period.label}: {sample_name} feature store is "
                    "missing production quantities: "
                    + ", ".join(missing)
                )
            #endif
        #endfor

        # Keep pi_epg resident for the aaogen reconstructed-companion
        # association below. dvcsgen is no longer needed after its features
        # have been built.
        del dv_epg

        support_values = parse_float_edges(
            str(args_dict["nsidis_pilot_probe_m2_values"]),
            "--nsidis-pilot-probe-m2-values",
        )
        central_support = float(
            args_dict["nsidis_pilot_probe_m2_max"]
        )
        if all(
            abs(float(x) - central_support) > 1.0e-12
            for x in support_values
        ):
            support_values = sorted(
                list(support_values) + [central_support]
            )
        #endif

        if bool(args_dict.get("nsidis_driver_study", False)):
            log(
                f"{period.label}: running dedicated denominator-driver study "
                "on a common data/aaogen/dvcsgen population."
            )
            for _sample_label, _feat in (
                ("nSidis data", ns_epg_f),
                ("aaogen pi0 MC", pi0_f),
                ("dvcsgen BH/DVCS MC", dvcs_f),
            ):
                if "stored_delta_phi_abs_rad" not in _feat:
                    raise RuntimeError(
                        f"{period.label}: dedicated driver study requires "
                        f"Delta_phi in {_sample_label}, but the branch was "
                        "not loaded/found. This check occurs before fitting."
                    )
                #endif
            #endfor

            driver_rows, driver_detail = run_nsidis_driver_study(
                period,
                ns_epg_f,
                pi0_f,
                dvcs_f,
                parent_mask=ns_epg_parent,
                ft_theta_max=ft_theta_max,
                max_probe_energy=float(
                    args_dict["nsidis_pilot_energy_max"]
                ),
                probe_m2_max=central_support,
                mm2_min=float(args_dict["den_fit_mm2_min"]),
                mm2_max=float(args_dict["den_fit_mm2_max"]),
                mm2_bins=int(args_dict["den_fit_mm2_bins"]),
                ptmiss_max=float(args_dict["disc_ptmiss_max"]),
                ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
                theta_max=float(args_dict["disc_theta_max"]),
                theta_bins=int(args_dict["disc_theta_bins"]),
                min_data_count=int(args_dict["den_min_data_count"]),
                min_template_count=int(
                    args_dict["den_min_template_count"]
                ),
                nuisance_shift_prior=float(
                    args_dict["morph_shift_prior_bins"]
                ),
                nuisance_sigma_prior=float(
                    args_dict["morph_sigma_prior_bins"]
                ),
                max_shift_bins=float(
                    args_dict["morph_max_shift_bins"]
                ),
                max_sigma_bins=float(
                    args_dict["morph_max_sigma_bins"]
                ),
            )
            write_rows_csv(
                driver_rows,
                outdir / "nsidis_driver_study.csv",
            )
            make_nsidis_driver_study_summary_canvas(
                period,
                driver_rows,
                outdir,
            )
            make_nsidis_driver_study_focus_canvas(
                period,
                driver_detail,
                outdir,
                focus_low=3.0,
                focus_high=4.5,
                mm2_min=float(args_dict["den_fit_mm2_min"]),
                mm2_max=float(args_dict["den_fit_mm2_max"]),
                mm2_bins=int(args_dict["den_fit_mm2_bins"]),
                ptmiss_max=float(args_dict["disc_ptmiss_max"]),
                ptmiss_bins=int(args_dict["disc_ptmiss_bins"]),
                theta_max=float(args_dict["disc_theta_max"]),
                theta_bins=int(args_dict["disc_theta_bins"]),
            )

            summary = {
                "period": period.key,
                "label": period.label,
                "epgamma_overlap": epg_overlap_summary,
                "eppi0_overlap": epi_overlap_summary,
                "eppi0_numerator_strategy": (
                    "driver-study-only; numerator efficiency not evaluated"
                ),
                "nsidis_pilot_fit_enabled": False,
                "nsidis_driver_study_enabled": True,
                "driver_study_rows": len(driver_rows),
                "wall_time_s": float(time.perf_counter() - t0),
            }
            log(
                f"{period.label}: dedicated denominator-driver study complete "
                f"in {summary['wall_time_s']:.1f} s."
            )
            return summary
        #endif

        clasdis_f = None
        if period.clasdis_epgamma_mc is not None:
            clasdis_path = Path(period.clasdis_epgamma_mc)
            if clasdis_path.exists():
                log(
                    f"{period.label}: reading CLASDIS epgammaX for Stage-1 "
                    "shape comparison."
                )
                clasdis_arrays, clasdis_tree, clasdis_total = read_branches(
                    str(clasdis_path),
                    EPG_REQUIRED,
                    EPG_OPTIONAL_CLASDIS_MC,
                    tree_name,
                    max_entries,
                )
                clasdis_epg = extract_epgamma(
                    clasdis_arrays,
                    angle_mode,
                )
                del clasdis_arrays
                clasdis_f = build_epgamma_denominator_features(
                    period,
                    clasdis_epg,
                    tag_min,
                    tag_max,
                    ft_theta_max,
                )
                attach_photon_angular_acceptance(
                    clasdis_f,
                    photon_acceptance,
                )
                clasdis_f["electron_p"] = np.asarray(
                    electron_momentum_from_p3(
                        clasdis_epg.electron_p3
                    ),
                    dtype=np.float32,
                )
                log(
                    f"{period.label}: CLASDIS tree '{clasdis_tree}', "
                    f"loaded {len(clasdis_epg.tag_energy):,}/"
                    f"{clasdis_total:,} entries."
                )
            else:
                log(
                    f"{period.label}: CLASDIS epgammaX is configured but not "
                    f"present yet; Stage-1 comparison will omit it: "
                    f"{clasdis_path}"
                )
            #endif
        #endif

        # Derive the nominal upper M_X^2(epX) edge from the FD CLASDIS
        # pi0/eta valley for this run period.  Use that same edge in both FD
        # and FT throughout all downstream denominator fits and efficiencies.
        fallback_mm2_max = float(args_dict["den_fit_mm2_max"])
        if clasdis_f is not None:
            derived_mm2_max = derive_clasdis_epx_valley_cut(
                period,
                clasdis_f,
                ft_theta_max=ft_theta_max,
                max_probe_energy=float(args_dict["nsidis_pilot_energy_max"]),
                fallback_cut=fallback_mm2_max,
            )
        else:
            derived_mm2_max = fallback_mm2_max
            log(
                f"{period.label}: WARNING: CLASDIS is unavailable, so the "
                f"fallback M_X^2(epX) upper cut {derived_mm2_max:.4f} GeV^2 "
                "will be used in both FD and FT."
            )
        #endif
        args_dict["den_fit_mm2_max"] = float(derived_mm2_max)

        common = {
            "period": period,
            "pi0_f": pi0_f,
            "dvcs_f": dvcs_f,
            "ft_theta_max": ft_theta_max,
            "mm2_min": float(args_dict["den_fit_mm2_min"]),
            "mm2_max": float(args_dict["den_fit_mm2_max"]),
            "support_values": support_values,
            "mm2_bins": int(args_dict["den_fit_mm2_bins"]),
            "ptmiss_max": float(args_dict["disc_ptmiss_max"]),
            "ptmiss_bins": int(args_dict["disc_ptmiss_bins"]),
            "min_data_count": int(args_dict["den_min_data_count"]),
            "min_template_count": int(
                args_dict["den_min_template_count"]
            ),
            "nuisance_shift_prior": float(
                args_dict["morph_shift_prior_bins"]
            ),
            "nuisance_sigma_prior": float(
                args_dict["morph_sigma_prior_bins"]
            ),
            "max_shift_bins": float(
                args_dict["morph_max_shift_bins"]
            ),
            "max_sigma_bins": float(
                args_dict["morph_max_sigma_bins"]
            ),
            "uncertainty_support": central_support,
        }

        make_nsidis_shape_comparison_canvases(
            period, ns_epg_f, pi0_f, dvcs_f, stage1_outdir,
            clasdis_f=clasdis_f,
            ft_theta_max=ft_theta_max,
            max_probe_energy=float(args_dict["nsidis_pilot_energy_max"]),
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
        )

        if clasdis_f is not None:
            del clasdis_f
            del clasdis_epg
        #endif

        # Keep the full visual diagnostic suite in ordinary production runs.
        # It is intentionally generated from a bounded, independent read so
        # the exploratory branch inventory does not bloat Stage-2/3 memory.
        run_grand_stage1_diagnostics_for_period(
            period,
            args_dict,
            production_root,
        )

        pilot_nsidis_rows = run_nsidis_three_variable_nominal_fits(
            period, ns_epg_f, pi0_f, dvcs_f,
            ft_theta_max=ft_theta_max,
            max_probe_energy=float(args_dict["nsidis_pilot_energy_max"]),
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            min_data_count=int(args_dict["den_min_data_count"]),
            min_template_count=int(args_dict["den_min_template_count"]),
            nuisance_shift_prior=float(args_dict["morph_shift_prior_bins"]),
            nuisance_sigma_prior=float(args_dict["morph_sigma_prior_bins"]),
            max_shift_bins=float(args_dict["morph_max_shift_bins"]),
            max_sigma_bins=float(args_dict["morph_max_sigma_bins"]),
            source_label="nsidis",
            pt_category_split=float(args_dict["pt_category_split"]),
            parent_mask=ns_epg_parent,
        )
        pilot_wagon_rows = run_nsidis_three_variable_nominal_fits(
            period, wagon_epg_f, pi0_f, dvcs_f,
            ft_theta_max=ft_theta_max,
            max_probe_energy=min(2.0, float(args_dict["nsidis_pilot_energy_max"])),
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            min_data_count=int(args_dict["den_min_data_count"]),
            min_template_count=int(args_dict["den_min_template_count"]),
            nuisance_shift_prior=float(args_dict["morph_shift_prior_bins"]),
            nuisance_sigma_prior=float(args_dict["morph_sigma_prior_bins"]),
            max_shift_bins=float(args_dict["morph_max_shift_bins"]),
            max_sigma_bins=float(args_dict["morph_max_sigma_bins"]),
            source_label="wagon_epgt2_reference",
            parent_mask=wagon_epg_parent,
        )

        write_rows_csv(
            pilot_nsidis_rows,
            stage2_outdir / f"nsidis_nominal_composition_fits_{period.key}.csv",
        )
        failed_production_bins = [
            row
            for row in pilot_nsidis_rows
            if int(row.get("fit_success", 0)) != 1
        ]
        if failed_production_bins:
            log(
                f"{period.label}: WARNING: "
                f"{len(failed_production_bins)} nominal composition "
                "bin(s) failed; they are omitted from fit canvases "
                "and cannot produce final scale factors."
            )
        #endif
        write_rows_csv(
            pilot_wagon_rows,
            stage2_outdir / f"wagon_nominal_composition_fits_{period.key}.csv",
        )

        make_nsidis_three_variable_fit_canvases(
            period, pilot_nsidis_rows,
            ns_epg_f, pi0_f, dvcs_f, stage2_outdir,
            ft_theta_max=ft_theta_max,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
        )

        make_nsidis_pi0_fraction_energy_canvas(
            period,
            pilot_nsidis_rows,
            stage2_outdir,
        )
        make_nsidis_pt_category_canvas(
            period,
            pilot_nsidis_rows,
            stage2_outdir,
        )
        log(
            f"{period.label}: wrote visible two-category pT composition "
            f"canvas to {stage2_outdir / ('canvas_pt_category_composition_' + period.key + '.png')}."
        )


        # --------------------------------------------------------------
        # aaogen MC reconstructed-companion efficiency.
        # --------------------------------------------------------------
        log(
            f"{period.label}: reading reconstructed aaogen eppi0 MC for "
            "the data/MC efficiency ratio."
        )
        mc_epi_arrays, _, _ = read_branches(
            period.eppi0_pi0_mc,
            EPPIO_REQUIRED,
            EPPIO_OPTIONAL_PI0_MC,
            tree_name,
            max_entries,
        )
        mc_epi = extract_eppi0(mc_epi_arrays, angle_mode)
        del mc_epi_arrays

        log(
            f"{period.label}: matching aaogen epgamma/eppi0 parents for "
            "the MC reconstructed-probe efficiency."
        )
        mc_matches = match_parent_kinematics(
            pi_epg,
            mc_epi,
            tag_min=tag_min,
            tag_max=tag_max,
            component_tolerance=float(
                args_dict["parent_component_tol"]
            ),
            nearest_distance_max=float(
                args_dict["parent_distance_max"]
            ),
            kdtree_workers=int(args_dict["kdtree_workers"]),
            query_chunk_size=int(
                args_dict["kdtree_query_chunk"]
            ),
        )
        mc_pair_np, mc_pair_counters = build_stage1_arrays(
            period,
            pi_epg,
            mc_epi,
            mc_matches,
        )
        mc_stage_lookup = build_mc_probe_stage_lookup(
            len(pi_epg.tag_energy),
            mc_pair_np,
            mgg_min=float(args_dict["assoc_mgg_min"]),
            mgg_max=float(args_dict["assoc_mgg_max"]),
            remainder_mass2_max=float(
                args_dict["stage3_tag_remainder_m2_max"]
            ),
            reco_probe_energy_min=float(
                args_dict["assoc_probe_energy_min"]
            ),
            probe_angle_max_deg=float(
                args_dict["stage3_probe_angle_max_deg"]
            ),
            probe_frac_energy_max=float(
                args_dict["stage3_probe_frac_energy_max"]
            ),
        )

        mc_eff_rows = build_aaogen_mc_efficiency_rows(
            period,
            pi0_f,
            mc_stage_lookup,
            ft_theta_max=ft_theta_max,
            max_probe_energy=float(
                args_dict["nsidis_pilot_energy_max"]
            ),
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
        )
        write_rows_csv(
            mc_eff_rows,
            outdir / "aaogen_mc_efficiency.csv",
        )
        log(
            f"{period.label}: aaogen MC efficiency built from "
            f"{len(mc_matches.epg_index):,} parent matches."
        )

        # Association-first numerator pilot.
        log(
            f"{period.label}: building association-first reconstructed-probe "
            "numerator pilot."
        )

        ns_mask = nsidis_central_tag_mask(
            ns_epg_f,
            ft_theta_max=ft_theta_max,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            max_probe_energy=float(args_dict["nsidis_pilot_energy_max"]),
        )
        wagon_mask = nsidis_central_tag_mask(
            wagon_epg_f,
            ft_theta_max=ft_theta_max,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            max_probe_energy=min(
                2.0, float(args_dict["nsidis_pilot_energy_max"])
            ),
        )

        ns_idx = np.flatnonzero(ns_mask)
        wagon_idx = np.flatnonzero(wagon_mask)

        ns_epg_selected = subset_epgamma_sample(ns_epg, ns_idx)
        wagon_epg_selected = subset_epgamma_sample(
            wagon_epg, wagon_idx
        )
        ns_feat_selected = subset_feature_store(
            ns_epg_f,
            ns_idx,
            len(ns_epg.tag_energy),
        )
        wagon_feat_selected = subset_feature_store(
            wagon_epg_f,
            wagon_idx,
            len(wagon_epg.tag_energy),
        )

        # Critical invariant: selected feature stores must retain the fixed
        # FT/FD acceptance metadata.  Fail here, before the expensive exact
        # eppi0 association, if a future refactor drops it.
        for selected_name, selected_feat in (
            ("nSidis", ns_feat_selected),
            ("wagon", wagon_feat_selected),
        ):
            ft_min, ft_max, fd_min, fd_max = feature_angular_acceptance(
                selected_feat,
                ft_theta_max,
            )
            expected = (2.4, 5.0, 6.0, 35.0)
            observed = (ft_min, ft_max, fd_min, fd_max)
            if not np.allclose(observed, expected, rtol=0.0, atol=1.0e-12):
                raise RuntimeError(
                    f"{period.label}: {selected_name} selected feature store "
                    f"lost fixed photon angular acceptance; "
                    f"observed={observed}, expected={expected}."
                )
            #endif
        #endfor

        log(
            f"{period.label}: exact eppi0 association on "
            f"{len(ns_idx):,} nSidis denominator-supported tags and "
            f"{len(wagon_idx):,} old-wagon tags."
        )

        ns_assoc = match_data_event_candidates(
            period,
            ns_epg_selected,
            ns_epi,
            tag_remainder_m2_max=float(
                args_dict["stage3_tag_remainder_m2_max"]
            ),
            reco_probe_energy_min=float(
                args_dict["assoc_probe_energy_min"]
            ),
            probe_angle_max_deg=float(
                args_dict["stage3_probe_angle_max_deg"]
            ),
            probe_frac_energy_max=float(
                args_dict["stage3_probe_frac_energy_max"]
            ),
        )
        wagon_assoc = match_data_event_candidates(
            period,
            wagon_epg_selected,
            wagon_epi,
            tag_remainder_m2_max=float(
                args_dict["stage3_tag_remainder_m2_max"]
            ),
            reco_probe_energy_min=float(
                args_dict["assoc_probe_energy_min"]
            ),
            probe_angle_max_deg=float(
                args_dict["stage3_probe_angle_max_deg"]
            ),
            probe_frac_energy_max=float(
                args_dict["stage3_probe_frac_energy_max"]
            ),
        )

        log_association_summary(period, "nSidis", ns_assoc)
        log_association_summary(period, "old-wagon", wagon_assoc)
        association_rows = (
            association_stage_rows(period, "nSidis", ns_assoc)
            + association_stage_rows(period, "old-wagon", wagon_assoc)
        )
        write_rows_csv(
            association_rows,
            outdir / "association_stage_summary.csv",
        )

        if int(
            ns_assoc.counters.get("probe_energy", 0)
        ) == 0:
            log(
                f"{period.label}: WARNING: no nSidis candidate survives the "
                "photon-like remainder + reconstructed-probe-energy stage. "
                "Skipping numerator-purity/SF construction for this period; "
                "association_stage_summary.csv contains the failure location."
            )

            summary = {
                "period": period.key,
                "label": period.label,
                "status": "invalid_nsidis_association_no_probe_energy_candidates",
                "wall_time_s": float(time.perf_counter() - t0),
                "association_counters": ns_assoc.counters,
                "wagon_association_counters": wagon_assoc.counters,
            }
            with (outdir / "nsidis_period_summary.json").open(
                "w", encoding="utf-8"
            ) as fout:
                json.dump(summary, fout, indent=2, default=json_default)
            #endif
            return summary
        #endif

        # --------------------------------------------------------------
        # Associated-numerator pi0 purity diagnostic.
        #
        # The eppi0X skim itself is restricted to 0.11 < M_gg < 0.16 GeV.
        # Fit the FINAL associated candidates in that available window with
        # a reconstructed-aaogen signal template + smooth combinatorial
        # background.  The result is deliberately NOT applied to epsilon_data
        # or SF_gamma in this version; inspect/validate it first.
        # --------------------------------------------------------------
        pi0_central_tag_mask = nsidis_central_tag_mask(
            pi0_f,
            ft_theta_max=ft_theta_max,
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            max_probe_energy=float(
                args_dict["nsidis_pilot_energy_max"]
            ),
        )

        log(
            f"{period.label}: fitting final-associated numerator "
            "M_gammagamma purity."
        )
        numerator_purity_rows, numerator_purity_detail = (
            build_associated_numerator_purity_rows(
                period,
                ns_assoc,
                mc_pair_np,
                pi0_f,
                central_pi0_tag_mask=pi0_central_tag_mask,
                max_probe_energy=float(
                    args_dict["nsidis_pilot_energy_max"]
                ),
                mgg_min=float(args_dict["assoc_mgg_min"]),
                mgg_max=float(args_dict["assoc_mgg_max"]),
                remainder_mass2_max=float(
                    args_dict["stage3_tag_remainder_m2_max"]
                ),
                reco_probe_energy_min=float(
                    args_dict["assoc_probe_energy_min"]
                ),
                probe_angle_max_deg=float(
                    args_dict["stage3_probe_angle_max_deg"]
                ),
                probe_frac_energy_max=float(
                    args_dict["stage3_probe_frac_energy_max"]
                ),
            )
        )
        write_rows_csv(
            numerator_purity_rows,
            outdir / "nsidis_numerator_pi0_purity.csv",
        )
        make_associated_numerator_purity_summary_canvas(
            period,
            numerator_purity_rows,
            outdir,
        )
        make_associated_numerator_mass_fit_canvases(
            period,
            numerator_purity_rows,
            numerator_purity_detail,
            outdir,
        )
        log(
            f"{period.label}: numerator-purity diagnostic complete; "
            "linear-background numerator purity has been applied to "
            "epsilon_data and SF_gamma."
        )

        ns_eff_rows = build_nsidis_data_efficiency_rows(
            period,
            pilot_nsidis_rows,
            ns_feat_selected,
            ns_assoc,
            ft_theta_max=ft_theta_max,
            max_probe_energy=float(args_dict["nsidis_pilot_energy_max"]),
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            source_label="nsidis",
            numerator_purity_rows=numerator_purity_rows,
        )
        wagon_eff_rows = build_nsidis_data_efficiency_rows(
            period,
            pilot_wagon_rows,
            wagon_feat_selected,
            wagon_assoc,
            ft_theta_max=ft_theta_max,
            max_probe_energy=min(
                2.0, float(args_dict["nsidis_pilot_energy_max"])
            ),
            mm2_min=float(args_dict["den_fit_mm2_min"]),
            mm2_max=float(args_dict["den_fit_mm2_max"]),
            probe_m2_max=central_support,
            source_label="wagon_epgt2_reference",
        )

        write_rows_csv(
            ns_eff_rows,
            outdir / "nsidis_data_efficiency_pilot.csv",
        )
        write_rows_csv(
            wagon_eff_rows,
            outdir / "wagon_data_efficiency_reference.csv",
        )
        sf_rows = combine_data_mc_scale_factor_rows(
            ns_eff_rows,
            mc_eff_rows,
        )
        wagon_sf_rows = combine_data_mc_scale_factor_rows(
            wagon_eff_rows,
            mc_eff_rows,
        )
        write_rows_csv(
            sf_rows,
            outdir / "nsidis_photon_efficiency_scale_factors.csv",
        )
        write_rows_csv(
            wagon_sf_rows,
            outdir / "wagon_photon_efficiency_scale_factor_reference.csv",
        )
        make_nsidis_scale_factor_canvas(
            period,
            sf_rows,
            outdir,
        )

        assoc_excl_summary = {
            "status": "diagnostic_removed_not_used_in_production"
        }

        with (
            outdir / "nsidis_numerator_association_summary.json"
        ).open("w") as f:
            json.dump(
                {
                    "period": period.key,
                    "label": period.label,
                    "central_probe_m2_support_GeV2": central_support,
                    "selected_nsidis_epgamma_tags": int(len(ns_idx)),
                    "selected_wagon_epgamma_tags": int(len(wagon_idx)),
                    "nsidis_association_counters": ns_assoc.counters,
                    "wagon_association_counters": wagon_assoc.counters,
                    "associated_eppi0_exclusivity": assoc_excl_summary,
                    "aaogen_mc_association": {
                        "parent_matches": int(
                            len(mc_matches.epg_index)
                        ),
                        "pair_counters": mc_pair_counters,
                        "selected_mc_efficiency_rows": int(
                            len(mc_eff_rows)
                        ),
                    },
                    "numerator_pi0_purity": {
                        "status": (
                            "applied to epsilon_data and SF_gamma"
                        ),
                        "signal_model": (
                            "final-associated reconstructed aaogen "
                            "M_gammagamma template"
                        ),
                        "nominal_background_model": "linear",
                        "alternative_background_model": "constant",
                        "available_mass_window_GeV": [
                            float(args_dict["assoc_mgg_min"]),
                            float(args_dict["assoc_mgg_max"]),
                        ],
                        "rows": int(len(numerator_purity_rows)),
                    },
                    "status": (
                        "preliminary association-first data/MC photon "
                        "efficiency scale factor formed; profiled fitted-pi0 "
                        "statistical uncertainty propagated"
                    ),
                },
                f,
                indent=2,
                allow_nan=True,
            )
        #endwith
    #endif

    summary = {
        "period": period.key,
        "label": period.label,
        "beam_energy_GeV": period.beam_energy,
        "nsidis_epgamma_path": ns_epg_path,
        "nsidis_eppi0_path": ns_epi_path,
        "overlap_parent_requirement": "electron momentum > 2 GeV",
        "theta_ep_requirement": f"theta_ep > {THETA_EP_MIN_DEG:g} deg",
        "photon_angular_acceptance": photon_acceptance.as_dict(),
        "epgamma_overlap": epg_overlap_summary,
        "eppi0_overlap": epi_overlap_summary,
        "eppi0_numerator_strategy": (
            "association-first; no global eppi0 exclusivity cut"
        ),
        "nsidis_pilot_fit_enabled": bool(
            args_dict.get("nsidis_pilot_fit", False)
        ),
        "production_denominator_driver": (
            "simultaneous Delta_phi + pTmiss + xF2"
        ),
        "wall_time_s": float(time.perf_counter() - t0),
        "status": (
            "nSidis denominator + association-first DATA/aaogen MC "
            "photon-efficiency scale-factor pilot; fitted-pi0 denominator "
            "uncertainty propagated; numerator-purity correction diagnostic "
            "only and not yet applied"
        ),
    }

    with (outdir / "nsidis_study_summary.json").open("w") as f:
        json.dump(summary, f, indent=2, allow_nan=True)
    #endwith

    log(
        f"{period.label}: nSidis study complete in "
        f"{summary['wall_time_s']:.1f} s."
    )
    return summary


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
    diagnostics_mode = str(args_dict.get("diagnostics", "selection"))
    full_diagnostics = diagnostics_mode == "full"

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
    candidate_model_rows: List[Dict[str, object]] = []
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

        # Derive the actual reconstructed-photon angular support directly from
        # this period's real epgamma data.  detector2=0 is FT, detector2=1 is FD.
        # These exact observed extrema define where a predicted probe could
        # physically have been reconstructed.
        photon_acceptance = infer_photon_angular_acceptance(
            data_epg,
            source=f"{period.label} epgamma data",
            )
        log(
            f"{period.label}: fixed photon angular acceptance: "
            "FT 2.4-5.0 deg; FD 6.0-35.0 deg."
        )

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
        for _feat in (pi0_f, data_f, dvcs_f):
            attach_photon_angular_acceptance(_feat, photon_acceptance)
        #endfor
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

        pi0_control: Dict[str, object] = {}
        if full_diagnostics:
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
            write_rows_csv(
                control_rows,
                stage2_dir / "pi0_control_validation.csv",
            )
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
                f"{period.label}: reconstructed-eppi0 control validation written."
            )
        else:
            log(
                f"{period.label}: selection diagnostics mode: skipping "
                "reconstructed-eppi0 control validation."
            )
        #endif
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
            run_closure=full_diagnostics,
        )
        write_rows_csv(shared_rows, stage2_dir / "denominator_shared_morphed_fits.csv")
        write_rows_csv(closure_rows, stage2_dir / "template_mixture_closure.csv")
        if not compact_plot_enabled(args_dict):
            make_stage2_template_mixture_closure_canvas(
                closure_rows,
                stage2_dir,
            )
        #endif

        # The discriminator decision is now settled: production uses only
        # M_X^2(ep) x pTmiss. The former theta-driver comparison remains
        # available only in --diagnostics full for archival/debug purposes.
        if full_diagnostics:
            # Always perform the focused discriminator-selection study.
            candidate_model_rows = run_stage2_candidate_model_study(
                period,
                data_f,
                pi0_f,
                dvcs_f,
                ft_theta_max,
                max_probe_energy,
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
            )
            write_rows_csv(
                candidate_model_rows,
                stage2_dir / "candidate_model_selection.csv",
            )
            make_candidate_model_summary_canvas(
                period,
                candidate_model_rows,
                stage2_dir,
            )
            make_selection_driver_fit_canvases(
                period,
                data_f,
                pi0_f,
                dvcs_f,
                stage2_dir,
                ft_theta_max,
                max_probe_energy,
                float(args_dict["den_fit_mm2_min"]),
                float(args_dict["den_fit_mm2_max"]),
                float(args_dict["den_fit_probe_m2_max"]),
                int(args_dict["den_fit_mm2_bins"]),
                float(args_dict["disc_ptmiss_max"]),
                int(args_dict["disc_ptmiss_bins"]),
                float(args_dict["disc_theta_max"]),
                int(args_dict["disc_theta_bins"]),
                int(args_dict["den_min_data_count"]),
                int(args_dict["den_min_template_count"]),
            )

        else:
            log(
                f"{period.label}: production discriminator fixed to "
                "M_X^2(ep) x pTmiss; skipped archived driver-comparison study."
            )
        #endif

        # The old nominal-driver canvases are redundant with the new selection
        # canvases because they hid Mx2. Keep them only in full/debug plotting.
        if not compact_plot_enabled(args_dict):
            make_actual_nominal_stage2_driver_canvases(
                period,
                shared_rows,
                data_f,
                pi0_f,
                dvcs_f,
                stage2_dir,
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
            )
        #endif

        if full_diagnostics:
            make_theta_gg_alternative_canvas(
                period,
                stage2_rows,
                stage2_dir,
            )
        #endif

        stage2_spread_rows = discriminator_spread_rows(stage2_rows)
        write_rows_csv(
            stage2_spread_rows,
            stage2_dir / "denominator_discriminator_spread.csv",
        )

        if full_diagnostics:
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
            log(
                f"{period.label}: wrote explicit morph nuisance profiles "
                f"({len(profile_rows)} rows)."
            )

            # Mixed-event accidental/wrong-photon diagnostics were removed.
            # The denominator study now uses only real data plus the AAO,
            # CLASDIS, and DVCSgen samples.

            if not compact_plot_enabled(args_dict):
                make_morph_profile_canvas(
                    period,
                    profile_rows,
                    stage2_dir,
                )
                make_ft_coarse_composition_canvas(
                    period,
                    shared_rows,
                    ft_coarse_rows,
                    ft_coarse_closure_rows,
                    stage2_dir,
                )
                make_stage2_canvases(
                    period,
                    data_f=data_f,
                    pi0_f=pi0_f,
                    dvcs_f=dvcs_f,
                    fit_rows=stage2_rows,
                    spread_rows=stage2_spread_rows,
                    mixed_rows=[],
                    outdir=stage2_dir,
                    ft_theta_max=ft_theta_max,
                    max_probe_energy=max_probe_energy,
                    mm2_min=float(args_dict["den_fit_mm2_min"]),
                    mm2_max=float(args_dict["den_fit_mm2_max"]),
                    probe_m2_max=float(args_dict["den_fit_probe_m2_max"]),
                )
            #endif
        else:
            log(
                f"{period.label}: selection diagnostics mode: skipped "
                "closure/coarse-FT/morph-profile diagnostics."
            )
        #endif

        stage2_summary = summarize_stage2(
            period,
            fit_rows=stage2_rows,
            spread_rows=stage2_spread_rows,
            mixed_rows=[],
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
            write_rows_csv(
                final_reference_input_rows(stage3_rows),
                stage3_dir / "final_reference_inputs.csv",
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

    wall_time_s = float(time.perf_counter() - t0)
    mc_association_summary["wall_time_s"] = wall_time_s

    log(
        f"{period.label}: total worker wall time = {wall_time_s:.1f} s."
    )

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




def final_reference_input_rows(
    stage3_rows: Sequence[Dict[str, object]],
) -> List[Dict[str, object]]:
    """
    Return a compact provenance table containing only quantities that directly
    enter the preliminary reference scale factor.

    For each period/region/energy bin:

        data denominator = fitted_pi0_fraction * data_tag_candidates
        epsilon_data     = data_clean_reconstructed_probe / data denominator
        epsilon_MC       = mc_clean_reconstructed_probe / mc_tag_candidates
        SF_gamma         = epsilon_data / epsilon_MC

    This table is intended to answer "what actually fed the final point?"
    without requiring inspection of the much larger development CSVs.
    """
    fields = (
        "period",
        "label",
        "region",
        "energy_low_GeV",
        "energy_high_GeV",
        "fitted_pi0_fraction",
        "stage2_deviance_per_ndof",
        "data_tag_candidates",
        "data_fitted_pi0_denominator",
        "data_clean_reconstructed_probe",
        "data_efficiency",
        "mc_tag_candidates",
        "mc_clean_reconstructed_probe",
        "mc_efficiency",
        "data_over_mc_scale_factor",
        "scale_factor_counting_error_provisional",
        "status",
    )

    out: List[Dict[str, object]] = []
    for row in stage3_rows:
        out.append({field: row.get(field, "") for field in fields})
    #endfor
    return out


def write_analysis_flow_guide(
    output_root: Path,
    stage2_outroot: Path,
    stage3_outroot: Path,
) -> None:
    """Write a short guide separating nominal inputs from diagnostics."""
    output_root.mkdir(parents=True, exist_ok=True)

    edges = [float(x) for x in stage2_energy_edges(2.0)]
    edge_text = ", ".join(f"{x:g}" for x in edges)

    guide = f"""PHOTON EFFICIENCY ANALYSIS FLOW
================================

NOMINAL quantities that feed the preliminary efficiency scale factor
--------------------------------------------------------------------
Energy-bin edges (GeV):
  {edge_text}

For EACH period, region, and energy bin:

1. Nominal denominator-composition fit
   File:
     {stage2_outroot}/denominator_shared_morphed_fits.csv

   The fit is a shared two-component aaogen-pi0 + dvcsgen-BH/DVCS fit using
   BOTH production 2D discriminator spaces simultaneously:
     - M_X^2(ep) x M_probe^2
     - M_X^2(ep) x pTmiss

   IMPORTANT: these two drivers share ONE fitted f_pi0.  That shared f_pi0 from
   the exact period/region/energy row is what Stage III uses.

   v045 is explicitly a discriminator-selection study. The existing nominal
   shared fit is retained ONLY so Stage III remains backward-compatible if it is
   requested. It is not being promoted as the final denominator model here.

   The high-value selection canvases are:
     - canvas_selection_fit_projections_ft.png
     - canvas_selection_fit_projections_fd_all.png
       These show M_X^2(ep) explicitly as well as pTmiss/theta_gamma_gamma.
     - canvas_candidate_model_summary.png
       Compares M_X^2 x pTmiss, M_X^2 x theta_gg, and a shared pTmiss+theta fit.
     - canvas_theta_gg_alternative.png

   In particular, theta_gg_1d tests theta_gamma_gamma without using M_X^2(ep).
   Neither theta_gg diagnostic currently changes the nominal f_pi0.

2. Data efficiency
     N_data_den = f_pi0 * N_data_tag_candidates
     epsilon_data = N_data_clean_reconstructed_probe / N_data_den

3. MC efficiency
     epsilon_MC =
       N_MC_clean_reconstructed_probe / N_MC_tag_candidates

   The MC numerator uses the internal aaogen epgamma<->eppi0 parent association.
   It is not a separate user-visible analysis stage.

4. Relative correction
     SF_gamma = epsilon_data / epsilon_MC

Compact table of EXACT quantities entering every displayed reference point:
  {stage3_outroot}/final_reference_inputs.csv

Full Stage-III numerical table:
  {stage3_outroot}/reference_efficiency_scale_factors.csv


DIAGNOSTIC ONLY -- these do NOT feed the nominal reference SF
--------------------------------------------------------------
- individual discriminator fits / discriminator spread;
- template-mixture closure studies;
- coarse FT 0.40-1.20 / 1.20-2.00 GeV studies;
- FT three-component wrong-photon diagnostic;
- morph nuisance scans;
- representative FT/FD pTmiss fit-overlay canvases;
- all --plot-mode full development canvases.

The coarse FT bins are diagnostics only.  The preliminary reference efficiency
continues to use the fine Stage-II energy bins listed above.
"""
    (output_root / "README_analysis_flow.txt").write_text(guide)


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



def merge_period_result(
    result: Dict[str, object],
    accum: Dict[str, List[Dict[str, object]]],
) -> None:
    """
    Merge one process_period() payload into aggregate containers.

    Serial and parallel execution intentionally share this exact function so
    result-key refactors cannot diverge between the two code paths.
    """
    required_keys = {
        "mc_association_summary",
        "stage2_summary",
        "stage2_rows",
        "stage2_spread_rows",
        "mixed_diag_rows",
        "profile_rows",
        "shared_rows",
        "closure_rows",
        "ft_coarse_rows",
        "ft_coarse_closure_rows",
        "ft_coarse_three_component_rows",
        "ft_coarse_three_component_closure_rows",
        "ft_coarse_three_component_summary_rows",
        "control_rows",
        "stage3_rows",
        "stage3_summary",
    }
    missing = sorted(required_keys - set(result))
    if missing:
        raise KeyError(
            "process_period() result payload is missing required key(s): "
            + ", ".join(missing)
        )
    #endif

    accum["mc_association_summaries"].append(
        result["mc_association_summary"]
    )

    if result["stage2_summary"] is not None:
        accum["stage2_summaries"].append(result["stage2_summary"])
        accum["all_stage2_rows"].extend(result["stage2_rows"])
        accum["all_stage2_spread_rows"].extend(
            result["stage2_spread_rows"]
        )
        accum["all_mixed_diag_rows"].extend(result["mixed_diag_rows"])
        accum["all_profile_rows"].extend(result["profile_rows"])
        accum["all_shared_rows"].extend(result["shared_rows"])
        accum["all_closure_rows"].extend(result["closure_rows"])
        accum["all_ft_coarse_rows"].extend(result["ft_coarse_rows"])
        accum["all_ft_coarse_closure_rows"].extend(
            result["ft_coarse_closure_rows"]
        )
        accum["all_ft_coarse_three_component_rows"].extend(
            result["ft_coarse_three_component_rows"]
        )
        accum["all_ft_coarse_three_component_closure_rows"].extend(
            result["ft_coarse_three_component_closure_rows"]
        )
        accum["all_ft_coarse_three_component_summary_rows"].extend(
            result["ft_coarse_three_component_summary_rows"]
        )
        accum["all_control_rows"].extend(result["control_rows"])
    #endif

    accum["all_stage3_rows"].extend(result["stage3_rows"])
    if result["stage3_summary"] is not None:
        accum["stage3_summaries"].append(result["stage3_summary"])
    #endif


def run_period_result_merge_self_test() -> None:
    """Fast regression test for the serial/parallel aggregation schema."""
    fake_result = {
        "mc_association_summary": {"period": "selftest"},
        "stage2_summary": {"period": "selftest"},
        "stage2_rows": [{"period": "selftest"}],
        "stage2_spread_rows": [],
        "mixed_diag_rows": [],
        "profile_rows": [],
        "shared_rows": [],
        "closure_rows": [],
        "ft_coarse_rows": [],
        "ft_coarse_closure_rows": [],
        "ft_coarse_three_component_rows": [],
        "ft_coarse_three_component_closure_rows": [],
        "ft_coarse_three_component_summary_rows": [],
        "control_rows": [],
        "stage3_rows": [{"period": "selftest"}],
        "stage3_summary": {"period": "selftest"},
    }
    fake_accum = {
        "mc_association_summaries": [],
        "stage2_summaries": [],
        "all_stage2_rows": [],
        "all_stage2_spread_rows": [],
        "all_mixed_diag_rows": [],
        "all_shared_rows": [],
        "all_closure_rows": [],
        "all_ft_coarse_rows": [],
        "all_ft_coarse_closure_rows": [],
        "all_ft_coarse_three_component_rows": [],
        "all_ft_coarse_three_component_closure_rows": [],
        "all_ft_coarse_three_component_summary_rows": [],
        "all_profile_rows": [],
        "all_control_rows": [],
        "all_stage3_rows": [],
        "stage3_summaries": [],
    }

    merge_period_result(fake_result, fake_accum)

    if len(fake_accum["mc_association_summaries"]) != 1:
        raise RuntimeError("Period-result self-test failed MC summary merge.")
    #endif
    if len(fake_accum["stage2_summaries"]) != 1:
        raise RuntimeError("Period-result self-test failed Stage-II summary merge.")
    #endif
    if len(fake_accum["all_stage2_rows"]) != 1:
        raise RuntimeError("Period-result self-test failed Stage-II row merge.")
    #endif
    if len(fake_accum["all_stage3_rows"]) != 1:
        raise RuntimeError("Period-result self-test failed Stage-III row merge.")
    #endif
    if len(fake_accum["stage3_summaries"]) != 1:
        raise RuntimeError("Period-result self-test failed Stage-III summary merge.")
    #endif


def main() -> int:
    args = parse_args()

    if args.max_entries < 0:
        raise ValueError("--max-entries must be >= 0.")
    #endif
    if args.grand_diagnostics_max_entries < 0:
        raise ValueError("--grand-diagnostics-max-entries must be >= 0.")
    #endif
    if args.pt_category_split <= 0.0:
        raise ValueError("--pt-category-split must be > 0.")
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

    # The nSidis path is now the DEFAULT production workflow.  The old
    # --nsidis-study flag is retained only for backwards-compatible command
    # lines; users no longer need to specify it.
    args.nsidis_study = True

    if args.nsidis_study:
        log(
            "nSidis production tag support: "
            f"{args.tag_min:.3g} <= E_gamma,tag < {args.tag_max:.3g} GeV."
        )
        # Carry the nSidis denominator extraction all the way through the
        # reconstructed-probe efficiency and epsilon_data/epsilon_MC scale
        # factor by default, including E_probe > 2 GeV.
        if (
            not args.nsidis_driver_study
            and not args.stage1_only
        ):
            args.nsidis_pilot_fit = True
        #endif

        if not selected:
            raise ValueError("No periods selected.")
        #endif
        if (
            (args.nsidis_epgamma is not None or args.nsidis_eppi0 is not None)
            and len(selected) != 1
        ):
            raise ValueError(
                "--nsidis-epgamma/--nsidis-eppi0 overrides require exactly one "
                "--period selection."
            )
        #endif
        if not (2.0 <= args.nsidis_pilot_energy_max <= 10.0):
            raise ValueError(
                "--nsidis-pilot-energy-max must lie in [2,10] GeV."
            )
        #endif
        if args.nsidis_pilot_probe_m2_max <= 0.0:
            raise ValueError(
                "--nsidis-pilot-probe-m2-max must be > 0."
            )
        #endif
        if not (
            0.50 <= args.nsidis_eppi0_target_wagon_retention < 1.0
        ):
            raise ValueError(
                "--nsidis-eppi0-target-wagon-retention must be in [0.5,1)."
            )
        #endif

        log(
            "Default production path: nSidis; "
            f"E_probe^pred endpoint = "
            f"{args.nsidis_pilot_energy_max:.3g} GeV. "
            "Only the old wagon overlap reference is capped at 2 GeV."
        )

        ns_args_dict = vars(args).copy()

        # Cheap refactor guard before opening any ROOT file.
        validate_nsidis_argument_contract(ns_args_dict)

        preflight_nsidis_study(selected, ns_args_dict)

        nsroot = Path(args.output_dir)
        nsroot.mkdir(parents=True, exist_ok=True)

        n_nsidis_workers = min(
            max(1, int(args.workers)),
            len(selected),
        )
        log(
            f"nSidis period-level parallelism: "
            f"{n_nsidis_workers} process(es)."
        )
        if (
            args.stage1_only
            and args.stage1_grand_diagnostics
            and n_nsidis_workers > 2
        ):
            log(
                "WARNING: grand Stage-1 diagnostics read dozens of extra "
                "branches; use --workers 2 if node memory is limited."
            )
        #endif

        summaries_by_period: Dict[str, Dict[str, object]] = {}

        if n_nsidis_workers == 1:
            for period in selected:
                summaries_by_period[period.key] = (
                    process_nsidis_study_period(
                        period,
                        ns_args_dict,
                        str(Path(args.output_dir)),
                    )
                )
            #endfor
        else:
            with ProcessPoolExecutor(
                max_workers=n_nsidis_workers
            ) as executor:
                future_to_period = {
                    executor.submit(
                        process_nsidis_study_period,
                        period,
                        ns_args_dict,
                        str(Path(args.output_dir)),
                    ): period
                    for period in selected
                }

                for future in as_completed(future_to_period):
                    period = future_to_period[future]
                    try:
                        summaries_by_period[period.key] = future.result()
                    except Exception as exc:
                        raise RuntimeError(
                            "Parallel nSidis photon-efficiency processing "
                            f"failed for {period.label}: {exc}"
                        ) from exc
                    #endtry
                #endfor
            #endwith
        #endif

        if args.stage1_only:
            stage1_summary_rows = []
            for period in selected:
                summary = summaries_by_period[period.key]
                stage1_summary_rows.append(
                    {
                        "period": period.key,
                        "label": period.label,
                        "stage1_only": 1,
                        "wall_time_s": float(
                            summary.get(
                                "wall_time_s",
                                float("nan"),
                            )
                        ),
                        "status": summary.get(
                            "status",
                            "stage1_only_complete",
                        ),
                    }
                )
            #endfor

            write_rows_csv(
                stage1_summary_rows,
                nsroot / "stage1_run_summary.csv",
            )
            log(
                "Done. Stage-1-only shape comparisons are in "
                f"{nsroot / 'stage1_shape_comparison'}."
            )
            if args.stage1_grand_diagnostics:
                log(
                    "Grand Stage-1 diagnostics are in "
                    f"{nsroot / 'stage1_grand_diagnostics'} "
                    "(one subdirectory per run period)."
                )
            #endif
            return 0
        #endif

        summary_rows = []
        for period in selected:
            summary = summaries_by_period[period.key]
            summary_rows.append({
                "period": period.key,
                "label": period.label,
                "epgamma_reference_event_recall": (
                    summary.get("epgamma_overlap", {}).get(
                        "reference_unique_event_recall",
                        float("nan"),
                    )
                    if isinstance(
                        summary.get("epgamma_overlap", {}), dict
                    )
                    else float("nan")
                ),
                "eppi0_reference_event_recall": (
                    summary.get("eppi0_overlap", {}).get(
                        "reference_unique_event_recall",
                        float("nan"),
                    )
                    if isinstance(
                        summary.get("eppi0_overlap", {}), dict
                    )
                    else float("nan")
                ),
                "eppi0_numerator_strategy": (
                    summary["eppi0_numerator_strategy"]
                ),
                "pilot_fit_enabled": int(
                    summary.get("nsidis_pilot_fit_enabled", False)
                ),
                "driver_study_enabled": int(
                    summary.get("nsidis_driver_study_enabled", False)
                ),
                "wall_time_s": summary["wall_time_s"],
            })
        #endfor

        write_rows_csv(
            summary_rows,
            nsroot / "photon_efficiency_run_summary.csv",
        )
        if args.nsidis_driver_study:
            make_cross_period_nsidis_driver_study_canvas(
                nsroot,
                selected,
            )
        elif args.nsidis_pilot_fit:
            make_nsidis_period_comparison(
                nsroot,
                selected,
            )
        #endif
        log(f"Done. Photon-efficiency outputs are in {nsroot}.")
        return 0
    #endif

    outroot = Path(args.output_dir)
    outroot.mkdir(parents=True, exist_ok=True)

    stage2_outroot = Path(args.stage2_output_dir)
    stage2_outroot.mkdir(parents=True, exist_ok=True)
    stage3_outroot = Path(args.stage3_output_dir)
    if not args.skip_stage3:
        stage3_outroot.mkdir(parents=True, exist_ok=True)
    #endif
    write_analysis_flow_guide(
        outroot,
        stage2_outroot,
        stage3_outroot,
    )

    run_internal_self_tests(outroot)
    run_period_result_merge_self_test()
    log("Internal period-result aggregation self-test passed.")

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
    log(
        f"Diagnostic workload: {args.diagnostics}. "
        + (
            "Focused discriminator-selection study only."
            if args.diagnostics == "selection"
            else "Full development diagnostics enabled."
        )
    )
    nominal_edges = stage2_energy_edges(
        min(2.0, float(args.den_probe_energy_max))
    )
    log(
        "Nominal reference probe-energy bin edges (GeV): "
        + ", ".join(f"{x:g}" for x in nominal_edges)
        + ". Coarse FT bins are diagnostic-only."
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
                "Period-specific FT/FD theta support is derived from reconstructed "
                "real photons using detector2 (0=FT, 1=FD); predicted probes outside "
                "those observed theta intervals are excluded. FD sectors retain exact "
                "wrapped phi intervals [330,30), [30,90), [90,150), [150,210), "
                "[210,270), [270,330)."
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
            "diagnostics_mode": str(args.diagnostics),
            "memory_model": (
                "sample-specific ROOT branch reads; slim sample objects; float32 "
                "persistent Stage-II histogram features; prompt release of large "
                "temporary/sample objects"
            ),
            "fit_acceleration": (
                "shared-morph fits use cached morphed templates and an exact "
                "convex score-equation solve for the profiled pi0 fraction; "
                "physics objective and bounds are unchanged"
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
                "The mixed-event wrong-photon stress-test machinery has been removed; Stage-II template diagnostics use only physical same-event samples."
            ),
        },
    }

    args_dict = vars(args).copy()

    accum: Dict[str, List[Dict[str, object]]] = {
        "mc_association_summaries": [],
        "stage2_summaries": [],
        "all_stage2_rows": [],
        "all_stage2_spread_rows": [],
        "all_mixed_diag_rows": [],
        "all_shared_rows": [],
        "all_closure_rows": [],
        "all_ft_coarse_rows": [],
        "all_ft_coarse_closure_rows": [],
        "all_ft_coarse_three_component_rows": [],
        "all_ft_coarse_three_component_closure_rows": [],
        "all_ft_coarse_three_component_summary_rows": [],
        "all_profile_rows": [],
        "all_control_rows": [],
        "all_stage3_rows": [],
        "stage3_summaries": [],
    }

    if n_processes == 1:
        for period in selected:
            result = process_period(
                period,
                args_dict,
                str(outroot),
                str(stage2_outroot),
            )
            merge_period_result(result, accum)
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
                    merge_period_result(result, accum)
                except Exception as exc:
                    for other in future_to_period:
                        other.cancel()
                    #endfor
                    raise RuntimeError(
                        f"Parallel photon-efficiency processing failed for "
                        f"{period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    mc_association_summaries = accum["mc_association_summaries"]
    stage2_summaries = accum["stage2_summaries"]
    all_stage2_rows = accum["all_stage2_rows"]
    all_stage2_spread_rows = accum["all_stage2_spread_rows"]
    all_mixed_diag_rows = accum["all_mixed_diag_rows"]
    all_shared_rows = accum["all_shared_rows"]
    all_closure_rows = accum["all_closure_rows"]
    all_ft_coarse_rows = accum["all_ft_coarse_rows"]
    all_ft_coarse_closure_rows = accum["all_ft_coarse_closure_rows"]
    all_ft_coarse_three_component_rows = (
        accum["all_ft_coarse_three_component_rows"]
    )
    all_ft_coarse_three_component_closure_rows = (
        accum["all_ft_coarse_three_component_closure_rows"]
    )
    all_ft_coarse_three_component_summary_rows = (
        accum["all_ft_coarse_three_component_summary_rows"]
    )
    all_profile_rows = accum["all_profile_rows"]
    all_control_rows = accum["all_control_rows"]
    all_stage3_rows = accum["all_stage3_rows"]
    stage3_summaries = accum["stage3_summaries"]

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
            write_rows_csv(
                final_reference_input_rows(all_stage3_rows),
                stage3_outroot / "final_reference_inputs.csv",
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
        #endif

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

    log(f"Done. Summary plots are in {outroot / 'summary'}.")
    log(f"Denominator-composition outputs are in {stage2_outroot}.")
    if not args.skip_stage3:
        log(f"Efficiency-reference outputs are in {stage3_outroot}.")
    #endif
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
