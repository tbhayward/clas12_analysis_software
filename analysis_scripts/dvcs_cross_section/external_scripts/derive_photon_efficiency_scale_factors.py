#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors.py

CLAS12 RGA photon-efficiency study -- clean refactor.

CURRENT IMPLEMENTED SCOPE
=========================
Stage 1: Shape Comparison

This script intentionally does ONE thing only:

    Compare six epgamma/event-shape observables in:
      1. loose nSidis epgamma DATA,
      2. dvcsgen epgamma MC,
      3. aaogen epgamma MC,

first after a genuinely minimal common reconstructed-tag selection, and then
again after one explicit exclusivity requirement:

    |M_X^2(epgamma)| < 0.075 GeV^2.

There are NO fits, NO eppi0 input files, NO efficiency calculation, NO
numerator association, NO bootstrap, and NO systematic extraction.

Minimal common selection
------------------------
    e_p > 2 GeV
    theta_ep > 5 deg
    2 <= E_gamma(tag) < 9.5 GeV

Importantly, Stage 1 does NOT require:
    - positive missing/probe energy,
    - E_probe^pred > 0.4 GeV,
    - predicted-probe FT/FD angular acceptance.

Those requirements belong to a later tag/probe-denominator stage.  Keeping
them out here preserves the genuine exclusive DVCS/BH population near
E_miss = 0 and makes the raw data/MC shape comparison interpretable.

Plotted observables
-------------------
    M_X^2(epgamma)    computed from P_beam + P_target - P_e - P_p - P_gamma
    M_X^2(ep)         computed from P_beam + P_target - P_e - P_p
    Delta_phi         stored branch
    pTmiss            stored branch
    theta_gamma_gamma stored branch, displayed in degrees
    Emiss2            stored branch (despite the name, this is E_miss)

Each period produces ONE 2x6 canvas:
    top row    = minimal selection only
    bottom row = same population after |M_X^2(epgamma)| < 0.075 GeV^2

The M_X^2(epgamma) top-row panel shows dashed vertical lines at +/-0.075
GeV^2.  The output directory is flat: at most one PNG per requested period.

Only one compact JSON file is written for the whole invocation:
    stage1_summary.json

Default output:
    output/photon_efficiency/stage1_shape_comparison/
"""

from __future__ import annotations

import argparse
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


# =============================================================================
# Fixed physics definitions
# =============================================================================

PROTON_MASS_GEV = 0.9382720813

ELECTRON_P_MIN_GEV = 2.0
THETA_EP_MIN_DEG = 5.0

TAG_E_MIN_GEV = 2.0
TAG_E_MAX_GEV = 9.5
MX2_EPGAMMA_ABS_MAX_GEV2 = 0.075

PROBE_E_MIN_GEV = 0.4
PROBE_E_MAX_GEV = 9.5

FT_THETA_MIN_DEG = 2.4
FT_THETA_MAX_DEG = 5.0

FD_THETA_MIN_DEG = 6.0
FD_THETA_MAX_DEG = 35.0

DEFAULT_TREE_NAME = "PhysicsEvents"
DEFAULT_OUTPUT_DIR = "output/photon_efficiency/stage1_shape_comparison"
DEFAULT_CHUNK_SIZE = 500_000


# =============================================================================
# Plot definitions
# =============================================================================

@dataclass(frozen=True)
class PlotVariable:
    branch: str
    title: str
    xlabel: str
    low: float
    high: float
    bins: int


PLOT_VARIABLES: Tuple[PlotVariable, ...] = (
    PlotVariable(
        "Mx2_epgamma",
        r"$M_X^2(ep\gamma)$",
        r"$M_X^2(ep\gamma)$ (GeV$^2$)",
        -0.20,
        0.50,
        100,
    ),
    PlotVariable(
        "Mx2_ep",
        r"$M_X^2(ep)$",
        r"$M_X^2(ep)$ (GeV$^2$)",
        -0.20,
        0.50,
        100,
    ),
    PlotVariable(
        "Delta_phi",
        r"$\Delta\phi$",
        r"$\Delta\phi$ (rad)",
        math.pi / 2.0,
        math.pi,
        100,
    ),
    PlotVariable(
        "pTmiss",
        r"$p_{T,\mathrm{miss}}$",
        r"$p_{T,\mathrm{miss}}$ (GeV)",
        0.0,
        1.20,
        100,
    ),
    PlotVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$",
        r"$\theta_{\gamma\gamma}$ (deg)",
        0.0,
        4.0,
        80,
    ),
    PlotVariable(
        "Emiss2",
        r"$E_{\mathrm{miss}}$",
        r"$E_{\mathrm{miss}}$ (GeV)",
        -2.0,
        4.0,
        100,
    ),
)

COMPUTED_PLOT_KEYS = frozenset(("Mx2_epgamma", "Mx2_ep"))

STORED_SHAPE_BRANCHES: Tuple[str, ...] = tuple(
    variable.branch
    for variable in PLOT_VARIABLES
    if variable.branch not in COMPUTED_PLOT_KEYS
)



# =============================================================================
# Run-period configuration
# =============================================================================

@dataclass(frozen=True)
class PeriodConfig:
    key: str
    label: str
    beam_energy_GeV: float
    nsidis_epgamma_data: str
    dvcsgen_epgamma_mc: str
    aaogen_epgamma_mc: str


PERIODS: Tuple[PeriodConfig, ...] = (
    PeriodConfig(
        key="fa18_inb",
        label="Fa18 Inb",
        beam_energy_GeV=10.604,
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_fa18_inb_epgamma.root"
        ),
        dvcsgen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "dvcsgen_rga_fa18_inb_epgamma_0.40GeV.root"
        ),
        aaogen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "bkg_rga_fa18_inb_epgamma_0.40GeV.root"
        ),
    ),
    PeriodConfig(
        key="fa18_out",
        label="Fa18 Out",
        beam_energy_GeV=10.604,
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_fa18_out_epgamma.root"
        ),
        dvcsgen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "dvcsgen_rga_fa18_out_epgamma_0.40GeV.root"
        ),
        aaogen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "bkg_rga_fa18_out_epgamma_0.40GeV.root"
        ),
    ),
    PeriodConfig(
        key="sp18_inb",
        label="Sp18 Inb",
        beam_energy_GeV=10.594,
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_sp18_inb_epgamma.root"
        ),
        dvcsgen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "dvcsgen_rga_sp18_inb_epgamma_0.40GeV.root"
        ),
        aaogen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "bkg_rga_sp18_inb_epgamma_0.40GeV.root"
        ),
    ),
    PeriodConfig(
        key="sp18_out",
        label="Sp18 Out",
        beam_energy_GeV=10.594,
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_sp18_out_epgamma.root"
        ),
        dvcsgen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "dvcsgen_rga_sp18_out_epgamma_0.40GeV.root"
        ),
        aaogen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "bkg_rga_sp18_out_epgamma_0.40GeV.root"
        ),
    ),
    PeriodConfig(
        key="sp19_inb",
        label="Sp19 Inb",
        beam_energy_GeV=10.200,
        nsidis_epgamma_data=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/"
            "dvcs/efficiency_study/nSidis_rga_sp19_inb_epgamma.root"
        ),
        dvcsgen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "dvcsgen_rga_sp19_inb_epgamma_0.40GeV.root"
        ),
        aaogen_epgamma_mc=(
            "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/"
            "dvcsgen/dvcsgen_files_greater_than_0.40GeV/"
            "bkg_rga_sp19_inb_epgamma_0.40GeV.root"
        ),
    ),
)

PERIOD_BY_KEY: Dict[str, PeriodConfig] = {
    period.key: period for period in PERIODS
}


# =============================================================================
# ROOT branches
# =============================================================================

KINEMATIC_BRANCHES: Tuple[str, ...] = (
    "e_p",
    "e_theta",
    "e_phi",
    "p1_p",
    "p1_theta",
    "p1_phi",
    "p2_p",
    "p2_theta",
    "p2_phi",
)

SHAPE_BRANCHES: Tuple[str, ...] = STORED_SHAPE_BRANCHES

REQUIRED_BRANCHES: Tuple[str, ...] = (
    *KINEMATIC_BRANCHES,
    *SHAPE_BRANCHES,
)


# =============================================================================
# Small utilities
# =============================================================================

def log(message: str) -> None:
    stamp = time.strftime("%H:%M:%S")
    print(f"[{stamp}] {message}", flush=True)


def cartesian_from_spherical(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    sin_theta = np.sin(theta)
    px = momentum * sin_theta * np.cos(phi)
    py = momentum * sin_theta * np.sin(phi)
    pz = momentum * np.cos(theta)
    return px, py, pz


def infer_angle_unit(
    e_theta: np.ndarray,
    e_phi: np.ndarray,
) -> str:
    """
    Infer whether stored particle angles are radians or degrees.

    CLAS12 analysis trees used here are normally radians.  The check exists
    only to make the new script fail safely if a differently formatted file is
    supplied.
    """
    finite_theta = np.asarray(e_theta, dtype=float)
    finite_theta = finite_theta[np.isfinite(finite_theta)]

    finite_phi = np.asarray(e_phi, dtype=float)
    finite_phi = finite_phi[np.isfinite(finite_phi)]

    if len(finite_theta) == 0 or len(finite_phi) == 0:
        return "rad"
    #endif

    theta_max = float(np.nanmax(np.abs(finite_theta)))
    phi_max = float(np.nanmax(np.abs(finite_phi)))

    if theta_max <= math.pi + 0.25 and phi_max <= 2.0 * math.pi + 0.5:
        return "rad"
    #endif

    if theta_max <= 180.0 + 1.0 and phi_max <= 360.0 + 1.0:
        return "deg"
    #endif

    raise RuntimeError(
        "Could not infer particle-angle units from e_theta/e_phi: "
        f"max |e_theta|={theta_max:g}, max |e_phi|={phi_max:g}."
    )


def to_radians(values: np.ndarray, angle_unit: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if angle_unit == "deg":
        return np.radians(array)
    #endif
    return array


def electron_proton_opening_angle_deg(
    e_px: np.ndarray,
    e_py: np.ndarray,
    e_pz: np.ndarray,
    p_px: np.ndarray,
    p_py: np.ndarray,
    p_pz: np.ndarray,
) -> np.ndarray:
    dot = e_px * p_px + e_py * p_py + e_pz * p_pz
    e_norm = np.sqrt(e_px * e_px + e_py * e_py + e_pz * e_pz)
    p_norm = np.sqrt(p_px * p_px + p_py * p_py + p_pz * p_pz)

    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = dot / (e_norm * p_norm)
    #endwith

    cosine = np.clip(cosine, -1.0, 1.0)
    return np.degrees(np.arccos(cosine))


def stage1_missing_kinematics(
    arrays: Dict[str, np.ndarray],
    beam_energy_GeV: float,
    angle_unit: str,
) -> Dict[str, np.ndarray]:
    """
    Compute the Stage-1 missing-system quantities from measured e, p, gamma.

    P_X(ep) =
        P_beam + P_target - P_e - P_p

    P_X(epgamma) =
        P_beam + P_target - P_e - P_p - P_gamma_tag

    Therefore:

        M_X^2(ep)      = P_X(ep)^2

        M_X^2(epgamma) = P_X(epgamma)^2

    In the tag/probe interpretation, P_X(epgamma) is also the predicted
    probe-photon four-vector.  Thus M_X^2(epgamma) is the same quantity that
    the legacy script called M_probe^2.

    The reconstructed tag photon is massless:
        E_gamma = |p_gamma|

    The electron mass is neglected:
        E_e = |p_e|
    """
    e_p = np.asarray(arrays["e_p"], dtype=float)
    p_p = np.asarray(arrays["p1_p"], dtype=float)
    g_p = np.asarray(arrays["p2_p"], dtype=float)

    e_theta = to_radians(arrays["e_theta"], angle_unit)
    e_phi = to_radians(arrays["e_phi"], angle_unit)

    p_theta = to_radians(arrays["p1_theta"], angle_unit)
    p_phi = to_radians(arrays["p1_phi"], angle_unit)

    g_theta = to_radians(arrays["p2_theta"], angle_unit)
    g_phi = to_radians(arrays["p2_phi"], angle_unit)

    e_px, e_py, e_pz = cartesian_from_spherical(
        e_p,
        e_theta,
        e_phi,
    )
    p_px, p_py, p_pz = cartesian_from_spherical(
        p_p,
        p_theta,
        p_phi,
    )
    g_px, g_py, g_pz = cartesian_from_spherical(
        g_p,
        g_theta,
        g_phi,
    )

    proton_E = np.sqrt(
        np.maximum(
            p_p * p_p + PROTON_MASS_GEV * PROTON_MASS_GEV,
            0.0,
        )
    )

    # Missing system after e and p.
    mx_ep_E = (
        beam_energy_GeV
        + PROTON_MASS_GEV
        - e_p
        - proton_E
    )
    mx_ep_px = -e_px - p_px
    mx_ep_py = -e_py - p_py
    mx_ep_pz = beam_energy_GeV - e_pz - p_pz

    mx2_ep = (
        mx_ep_E * mx_ep_E
        - mx_ep_px * mx_ep_px
        - mx_ep_py * mx_ep_py
        - mx_ep_pz * mx_ep_pz
    )

    # Missing system after e, p, and reconstructed gamma_tag.
    pred_E = mx_ep_E - g_p
    pred_px = mx_ep_px - g_px
    pred_py = mx_ep_py - g_py
    pred_pz = mx_ep_pz - g_pz

    pred_p2 = (
        pred_px * pred_px
        + pred_py * pred_py
        + pred_pz * pred_pz
    )
    pred_p = np.sqrt(np.maximum(pred_p2, 0.0))

    mx2_epgamma = pred_E * pred_E - pred_p2

    with np.errstate(divide="ignore", invalid="ignore"):
        cos_theta = pred_pz / pred_p
    #endwith

    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    pred_theta_deg = np.degrees(np.arccos(cos_theta))

    return {
        "pred_probe_energy": pred_E,
        "pred_probe_theta_deg": pred_theta_deg,
        "Mx2_epgamma": mx2_epgamma,
        "Mx2_ep": mx2_ep,
    }



def common_support_mask(
    arrays: Dict[str, np.ndarray],
    beam_energy_GeV: float,
    angle_unit: str,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Apply the genuinely minimal Stage-1 common selection only.

    The missing/probe four-vector is NOT used to select events here.
    In particular there is no probe-energy or predicted-probe angular cut.
    """
    e_p = np.asarray(arrays["e_p"], dtype=float)
    p_p = np.asarray(arrays["p1_p"], dtype=float)
    g_E = np.asarray(arrays["p2_p"], dtype=float)

    e_theta = to_radians(arrays["e_theta"], angle_unit)
    e_phi = to_radians(arrays["e_phi"], angle_unit)
    p_theta = to_radians(arrays["p1_theta"], angle_unit)
    p_phi = to_radians(arrays["p1_phi"], angle_unit)

    e_px, e_py, e_pz = cartesian_from_spherical(
        e_p, e_theta, e_phi
    )
    p_px, p_py, p_pz = cartesian_from_spherical(
        p_p, p_theta, p_phi
    )

    theta_ep_deg = electron_proton_opening_angle_deg(
        e_px,
        e_py,
        e_pz,
        p_px,
        p_py,
        p_pz,
    )

    finite = (
        np.isfinite(e_p)
        & np.isfinite(p_p)
        & np.isfinite(g_E)
        & np.isfinite(theta_ep_deg)
    )

    electron = finite & (e_p > ELECTRON_P_MIN_GEV)
    opening = electron & (theta_ep_deg > THETA_EP_MIN_DEG)
    accepted = (
        opening
        & (g_E >= TAG_E_MIN_GEV)
        & (g_E < TAG_E_MAX_GEV)
    )

    counters = {
        "input": int(len(e_p)),
        "finite": int(np.count_nonzero(finite)),
        "electron_p_gt_2": int(np.count_nonzero(electron)),
        "theta_ep_gt_5deg": int(np.count_nonzero(opening)),
        "tag_2_to_9p5_GeV": int(np.count_nonzero(accepted)),
    }

    return accepted, counters


# =============================================================================
# Histogram accumulation
# =============================================================================

@dataclass
class HistogramResult:
    sample_key: str
    sample_label: str
    file_path: str
    entries_read: int
    selection_counts: Dict[str, int]
    counts_before: Dict[str, np.ndarray]
    counts_after: Dict[str, np.ndarray]
    angle_unit: str


def empty_histograms() -> Dict[str, np.ndarray]:
    return {
        variable.branch: np.zeros(variable.bins, dtype=np.int64)
        for variable in PLOT_VARIABLES
    }


def accumulate_shape_histograms(
    file_path: str,
    sample_key: str,
    sample_label: str,
    beam_energy_GeV: float,
    tree_name: str,
    max_entries: int,
    chunk_size: int,
) -> HistogramResult:
    """
    Stream the minimal Stage-1 branch set once.

    Both the pre-exclusivity and post-exclusivity histograms are filled during
    the same pass, so the second plot row costs essentially no extra I/O.
    """
    counts_before = empty_histograms()
    counts_after = empty_histograms()

    totals = {
        "input": 0,
        "finite": 0,
        "electron_p_gt_2": 0,
        "theta_ep_gt_5deg": 0,
        "tag_2_to_9p5_GeV": 0,
        "after_mx2_epgamma_cut": 0,
    }

    entries_read = 0
    angle_unit: Optional[str] = None

    with uproot.open(file_path) as root_file:
        tree = root_file[tree_name]
        tree_entries = int(tree.num_entries)
        stop = (
            tree_entries
            if max_entries <= 0
            else min(tree_entries, int(max_entries))
        )

        for arrays_raw in tree.iterate(
            expressions=list(REQUIRED_BRANCHES),
            entry_start=0,
            entry_stop=stop,
            step_size=int(chunk_size),
            library="np",
        ):
            arrays = {
                name: np.asarray(value)
                for name, value in arrays_raw.items()
            }

            n = len(arrays["e_p"])
            if n == 0:
                continue
            #endif

            if angle_unit is None:
                angle_unit = infer_angle_unit(
                    arrays["e_theta"],
                    arrays["e_phi"],
                )
            #endif

            minimal, chunk_counts = common_support_mask(
                arrays,
                beam_energy_GeV,
                angle_unit,
            )
            missing = stage1_missing_kinematics(
                arrays,
                beam_energy_GeV,
                angle_unit,
            )

            mx2_epgamma = np.asarray(
                missing["Mx2_epgamma"],
                dtype=float,
            )
            after = (
                minimal
                & np.isfinite(mx2_epgamma)
                & (np.abs(mx2_epgamma) < MX2_EPGAMMA_ABS_MAX_GEV2)
            )

            entries_read += n
            for key in (
                "input",
                "finite",
                "electron_p_gt_2",
                "theta_ep_gt_5deg",
                "tag_2_to_9p5_GeV",
            ):
                totals[key] += int(chunk_counts[key])
            #endfor
            totals["after_mx2_epgamma_cut"] += int(np.count_nonzero(after))

            for variable in PLOT_VARIABLES:
                if variable.branch in COMPUTED_PLOT_KEYS:
                    values = np.asarray(
                        missing[variable.branch],
                        dtype=float,
                    )
                else:
                    values = np.asarray(
                        arrays[variable.branch],
                        dtype=float,
                    )
                #endif

                # Stored angular branches follow the same angular convention
                # as the event tree.  Display theta_gamma_gamma in degrees.
                if variable.branch == "theta_gamma_gamma":
                    if angle_unit == "rad":
                        values = np.degrees(values)
                    #endif
                #endif

                finite_values = np.isfinite(values)

                before_values = values[minimal & finite_values]
                after_values = values[after & finite_values]

                hist_before, _ = np.histogram(
                    before_values,
                    bins=variable.bins,
                    range=(variable.low, variable.high),
                )
                hist_after, _ = np.histogram(
                    after_values,
                    bins=variable.bins,
                    range=(variable.low, variable.high),
                )

                counts_before[variable.branch] += hist_before.astype(np.int64)
                counts_after[variable.branch] += hist_after.astype(np.int64)
            #endfor
        #endfor

    if angle_unit is None:
        angle_unit = "unknown"
    #endif

    return HistogramResult(
        sample_key=sample_key,
        sample_label=sample_label,
        file_path=file_path,
        entries_read=entries_read,
        selection_counts=totals,
        counts_before=counts_before,
        counts_after=counts_after,
        angle_unit=angle_unit,
    )


# =============================================================================
# Preflight -- intentionally complete before any large ROOT read
# =============================================================================

def inspect_root_file(
    file_path: str,
    tree_name: str,
) -> Dict[str, object]:
    path = Path(file_path)

    if not path.exists():
        raise FileNotFoundError(f"Required ROOT file does not exist: {file_path}")
    #endif

    try:
        with uproot.open(file_path) as root_file:
            if tree_name not in root_file:
                available = [str(key).split(";")[0] for key in root_file.keys()]
                raise RuntimeError(
                    f"{file_path}: tree '{tree_name}' is missing. "
                    f"Available keys: {available}"
                )
            #endif

            tree = root_file[tree_name]
            available_branches = set(str(key) for key in tree.keys())
            missing = [
                branch
                for branch in REQUIRED_BRANCHES
                if branch not in available_branches
            ]

            if missing:
                raise RuntimeError(
                    f"{file_path}: tree '{tree_name}' is missing required "
                    f"Stage-1 branches: {', '.join(missing)}"
                )
            #endif

            return {
                "path": file_path,
                "tree": tree_name,
                "entries": int(tree.num_entries),
                "required_branches_present": True,
            }
    except Exception as exc:
        raise RuntimeError(
            f"Preflight failed for {file_path}: {exc}"
        ) from exc
    #endtry


def preflight_periods(
    periods: Sequence[PeriodConfig],
    tree_name: str,
) -> Dict[str, Dict[str, Dict[str, object]]]:
    """
    Check EVERY requested file/tree/branch before any expensive processing.

    This exists specifically to prevent a 20-minute run from failing because a
    late period is missing a branch or input file.
    """
    report: Dict[str, Dict[str, Dict[str, object]]] = {}

    log("Preflight: checking all requested Stage-1 ROOT inputs.")

    for period in periods:
        report[period.key] = {}

        samples = (
            (
                "data",
                period.nsidis_epgamma_data,
            ),
            (
                "dvcsgen",
                period.dvcsgen_epgamma_mc,
            ),
            (
                "aaogen",
                period.aaogen_epgamma_mc,
            ),
        )

        for sample_key, file_path in samples:
            info = inspect_root_file(file_path, tree_name)
            report[period.key][sample_key] = info

            log(
                f"Preflight OK: {period.label} {sample_key}: "
                f"{info['entries']:,} entries."
            )
        #endfor
    #endfor

    log("Preflight complete: all requested files, trees, and branches are valid.")
    return report


# =============================================================================
# Plotting
# =============================================================================

def normalized_histogram(
    counts: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    values = np.asarray(counts, dtype=float)
    total = float(np.sum(values))

    if total <= 0.0:
        return (
            np.zeros_like(values, dtype=float),
            np.zeros_like(values, dtype=float),
        )
    #endif

    normalized = values / total
    errors = np.sqrt(values) / total
    return normalized, errors


def draw_period_canvas(
    period: PeriodConfig,
    results: Dict[str, HistogramResult],
    outdir: Path,
) -> Path:
    """
    Draw one compact 2x6 canvas:
      top    = minimal selection only
      bottom = after |M_X^2(epgamma)| < 0.075 GeV^2
    """
    fig, axes = plt.subplots(
        2,
        len(PLOT_VARIABLES),
        figsize=(24.0, 8.0),
        squeeze=False,
    )

    data = results["data"]
    dvcs = results["dvcsgen"]
    aaogen = results["aaogen"]

    for row, (row_name, count_attr) in enumerate(
        (
            ("Minimal selection", "counts_before"),
            (
                rf"After $|M_X^2(ep\gamma)|<{MX2_EPGAMMA_ABS_MAX_GEV2:.3f}$ GeV$^2$",
                "counts_after",
            ),
        )
    ):
        for col, variable in enumerate(PLOT_VARIABLES):
            ax = axes[row, col]
            edges = np.linspace(
                variable.low,
                variable.high,
                variable.bins + 1,
            )
            centers = 0.5 * (edges[:-1] + edges[1:])

            data_counts = getattr(data, count_attr)[variable.branch]
            dvcs_counts = getattr(dvcs, count_attr)[variable.branch]
            aaogen_counts = getattr(aaogen, count_attr)[variable.branch]

            data_y, data_err = normalized_histogram(data_counts)
            dvcs_y, _ = normalized_histogram(dvcs_counts)
            aaogen_y, _ = normalized_histogram(aaogen_counts)

            ax.step(
                edges[:-1],
                dvcs_y,
                where="post",
                linewidth=1.35,
                color="tab:blue",
                label="DVCS MC",
            )
            ax.step(
                edges[:-1],
                aaogen_y,
                where="post",
                linewidth=1.35,
                color="red",
                label=r"$e\pi^0$ MC as $ep\gamma$",
            )
            ax.errorbar(
                centers,
                data_y,
                yerr=data_err,
                fmt="o",
                markersize=2.2,
                linewidth=0.5,
                capsize=0,
                color="black",
                label=r"$e'p'\gamma$ data",
                zorder=5,
            )

            ax.set_xlim(variable.low, variable.high)

            if row == 0 and variable.branch == "Mx2_epgamma":
                ax.axvline(
                    -MX2_EPGAMMA_ABS_MAX_GEV2,
                    linestyle="--",
                    linewidth=1.0,
                    color="0.35",
                    label=(
                        rf"$|M_X^2(ep\gamma)|"
                        rf"<{MX2_EPGAMMA_ABS_MAX_GEV2:.3f}$ GeV$^2$"
                    ),
                )
                ax.axvline(
                    +MX2_EPGAMMA_ABS_MAX_GEV2,
                    linestyle="--",
                    linewidth=1.0,
                    color="0.35",
                )
            #endif

            ax.set_xlabel(variable.xlabel)
            if col == 0:
                ax.set_ylabel(
                    f"{row_name}\nunit-normalized entries / bin"
                )
            #endif
            ax.grid(alpha=0.18)
        #endfor
    #endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()

    n_data_before = data.selection_counts["tag_2_to_9p5_GeV"]
    n_dvcs_before = dvcs.selection_counts["tag_2_to_9p5_GeV"]
    n_aaogen_before = aaogen.selection_counts["tag_2_to_9p5_GeV"]

    n_data_after = data.selection_counts["after_mx2_epgamma_cut"]
    n_dvcs_after = dvcs.selection_counts["after_mx2_epgamma_cut"]
    n_aaogen_after = aaogen.selection_counts["after_mx2_epgamma_cut"]

    fig.suptitle(
        f"Stage 1 shape comparison: {period.label}\n"
        r"minimal selection: $E_e>2$ GeV, $\theta_{ep}>5^\circ$, "
        r"$2\leq E_{\gamma,\mathrm{tag}}<9.5$ GeV"
        "\n"
        f"before cut: data={n_data_before:,}, DVCS MC={n_dvcs_before:,}, "
        rf"$e\pi^0$ MC={n_aaogen_before:,}; "
        f"after cut: data={n_data_after:,}, DVCS MC={n_dvcs_after:,}, "
        rf"$e\pi^0$ MC={n_aaogen_after:,}",
        fontsize=13,
        y=0.985,
    )

    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.885),
        ncol=4,
        frameon=False,
        fontsize=9.5,
    )

    # Deliberately reserve much less dead space than the previous version.
    fig.subplots_adjust(
        left=0.055,
        right=0.992,
        bottom=0.085,
        top=0.825,
        wspace=0.30,
        hspace=0.38,
    )

    outpath = outdir / f"stage1_shape_comparison_{period.key}.png"
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
    return outpath


# =============================================================================
# Per-period execution
# =============================================================================

def process_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
) -> Dict[str, object]:
    t0 = time.perf_counter()

    outdir = Path(str(args_dict["output_dir"]))
    outdir.mkdir(parents=True, exist_ok=True)

    sample_specs = (
        ("data", r"$e'p'\gamma$ data", period.nsidis_epgamma_data),
        ("dvcsgen", "DVCS MC", period.dvcsgen_epgamma_mc),
        ("aaogen", r"$e\pi^0$ MC as $ep\gamma$", period.aaogen_epgamma_mc),
    )

    results: Dict[str, HistogramResult] = {}

    for sample_key, sample_label, file_path in sample_specs:
        log(f"{period.label}: streaming {sample_key} Stage-1 shapes.")

        result = accumulate_shape_histograms(
            file_path=file_path,
            sample_key=sample_key,
            sample_label=sample_label,
            beam_energy_GeV=period.beam_energy_GeV,
            tree_name=str(args_dict["tree_name"]),
            max_entries=int(args_dict["max_entries"]),
            chunk_size=int(args_dict["chunk_size"]),
        )
        results[sample_key] = result

        log(
            f"{period.label}: {sample_key}: read {result.entries_read:,}; "
            f"minimal={result.selection_counts['tag_2_to_9p5_GeV']:,}; "
            f"after Mx2 cut={result.selection_counts['after_mx2_epgamma_cut']:,}."
        )
    #endfor

    canvas = draw_period_canvas(period, results, outdir)

    summary = {
        "period": period.key,
        "label": period.label,
        "canvas": str(canvas),
        "samples": {
            key: {
                "entries_read": result.entries_read,
                "minimal_selected": result.selection_counts["tag_2_to_9p5_GeV"],
                "after_mx2_epgamma_cut": result.selection_counts["after_mx2_epgamma_cut"],
            }
            for key, result in results.items()
        },
        "wall_time_s": float(time.perf_counter() - t0),
    }

    log(
        f"{period.label}: Stage 1 complete in "
        f"{summary['wall_time_s']:.1f} s."
    )
    return summary


# =============================================================================
# CLI
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "CLAS12 RGA photon-efficiency clean refactor: "
            "Stage 1 shape comparison only."
        )
    )

    parser.add_argument(
        "--period",
        action="append",
        choices=tuple(PERIOD_BY_KEY.keys()),
        help=(
            "Run period to process. Repeat for multiple periods. "
            "If omitted, all five configured RGA periods are processed."
        ),
    )

    parser.add_argument(
        "--max-entries",
        type=int,
        default=0,
        help=(
            "Maximum entries read from EACH ROOT file. "
            "0 means all entries. Default: 0."
        ),
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=2,
        help=(
            "Independent run periods processed in parallel. "
            "Default: 2 to avoid I/O contention."
        ),
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=DEFAULT_CHUNK_SIZE,
        help=(
            "uproot streaming chunk size. "
            f"Default: {DEFAULT_CHUNK_SIZE:,} entries."
        ),
    )

    parser.add_argument(
        "--tree-name",
        default=DEFAULT_TREE_NAME,
        help=f"ROOT tree name. Default: {DEFAULT_TREE_NAME}.",
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Stage-1 output directory. "
            f"Default: {DEFAULT_OUTPUT_DIR}."
        ),
    )

    return parser


def selected_periods(
    requested: Optional[Sequence[str]],
) -> List[PeriodConfig]:
    if not requested:
        return list(PERIODS)
    #endif

    ordered: List[PeriodConfig] = []
    seen = set()

    for key in requested:
        if key in seen:
            continue
        #endif
        ordered.append(PERIOD_BY_KEY[key])
        seen.add(key)
    #endfor

    return ordered


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    periods = selected_periods(args.period)

    if args.max_entries < 0:
        parser.error("--max-entries must be >= 0.")
    #endif

    if args.workers < 1:
        parser.error("--workers must be >= 1.")
    #endif

    if args.chunk_size < 1:
        parser.error("--chunk-size must be >= 1.")
    #endif

    Path(args.output_dir).mkdir(
        parents=True,
        exist_ok=True,
    )

    log("Stage 1 only: six epgamma shape comparisons; no fits or efficiency calculation.")
    log("Minimal selection: e_p>2 GeV, theta_ep>5 deg, 2<=E_tag<9.5 GeV; no probe-energy or probe-angle cut.")

    # Complete preflight BEFORE any large read.
    preflight_periods(periods, args.tree_name)

    args_dict: Dict[str, object] = {
        "max_entries": args.max_entries,
        "workers": args.workers,
        "chunk_size": args.chunk_size,
        "tree_name": args.tree_name,
        "output_dir": args.output_dir,
    }

    workers = min(
        int(args.workers),
        len(periods),
    )

    summaries: Dict[str, Dict[str, object]] = {}

    if workers == 1:
        for period in periods:
            summaries[period.key] = process_period(
                period,
                args_dict,
            )
        #endfor
    else:
        log(
            f"Period-level parallelism: {workers} process(es)."
        )

        with ProcessPoolExecutor(
            max_workers=workers
        ) as executor:
            futures = {
                executor.submit(
                    process_period,
                    period,
                    args_dict,
                ): period
                for period in periods
            }

            for future in as_completed(futures):
                period = futures[future]
                try:
                    summaries[period.key] = future.result()
                except Exception as exc:
                    # This should now be rare because the complete ROOT
                    # schema preflight happens before workers launch.
                    raise RuntimeError(
                        f"Stage 1 failed for {period.label}: {exc}"
                    ) from exc
                #endtry
            #endfor
        #endwith
    #endif

    compact_summary = {
        "stage": "stage1_shape_comparison",
        "selection": {
            "electron_p_min_GeV": ELECTRON_P_MIN_GEV,
            "theta_ep_min_deg": THETA_EP_MIN_DEG,
            "tag_energy_GeV": [TAG_E_MIN_GEV, TAG_E_MAX_GEV],
            "mx2_epgamma_abs_max_GeV2": MX2_EPGAMMA_ABS_MAX_GEV2,
        },
        "periods": [
            summaries[period.key]
            for period in periods
            if period.key in summaries
        ],
    }

    with (
        Path(args.output_dir) / "stage1_summary.json"
    ).open("w", encoding="utf-8") as stream:
        json.dump(compact_summary, stream, indent=2)
    #endwith

    log(
        "Done. Stage-1 outputs are in "
        f"{Path(args.output_dir)}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
