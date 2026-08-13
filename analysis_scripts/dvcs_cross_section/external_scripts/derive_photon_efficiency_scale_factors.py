#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors.py

CLAS12 RGA photon-efficiency study -- clean refactor.

CURRENT IMPLEMENTED SCOPE
=========================
Stage 1: Shape Comparison
Stage 2: Integrated DVCS/pi0 Template Fits

This script intentionally does ONE thing only:

    Compare five epgamma/event-shape observables in:
      1. loose nSidis epgamma DATA,
      2. dvcsgen epgamma MC,
      3. aaogen epgamma MC,

first after a genuinely minimal common reconstructed-tag selection, and then
again after one explicit exclusivity requirement:

    |M_X^2(epgamma)| < 0.075 GeV^2.

Stage 2 fits the post-exclusivity Delta_phi, pTmiss, and Emiss shapes as
DVCS+pi0 template mixtures with one shared response morph per variable. There is still NO eppi0 numerator efficiency,
NO bootstrap, and NO final systematic extraction.

Stage-1 interpretation note
---------------------------
This clean Stage 1 intentionally retains genuine exclusive BH/DVCS events
near E_miss = 0 and M_X^2(epgamma) = 0. Earlier legacy-development canvases
required a positive predicted-probe energy and predicted-probe detector
acceptance, which removed much of that exclusive DVCS peak. Therefore the
raw M_X^2(epgamma) distribution here is expected to be much more sharply
peaked at zero.

Minimal common selection
------------------------
    e_p > 2 GeV
    theta_egamma > 5 deg
    0.4 <= E_gamma(tag) < 9.5 GeV

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
    Delta_phi         stored branch
    pTmiss            stored branch
    Emiss2            stored branch (despite the name, this is E_miss)
    E_gamma,tag       p2_p branch (reconstructed/tag photon energy)

Each period produces TWO 2x5 canvases, split by detected tag photon:
    FT canvas: detector2 == 0
    FD canvas: detector2 == 1

In each canvas:
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

try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.optimize import minimize, minimize_scalar
except Exception:
    gaussian_filter1d = None
    minimize = None
    minimize_scalar = None


# =============================================================================
# Fixed physics definitions
# =============================================================================

PROTON_MASS_GEV = 0.9382720813

ELECTRON_P_MIN_GEV = 2.0
THETA_EGAMMA_MIN_DEG = 5.0

TAG_E_MIN_GEV = 0.4
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

TAG_REGIONS: Tuple[Tuple[str, str, int], ...] = (
    ("FT", "FT tag", 0),
    ("FD", "FD tag", 1),
)


# =============================================================================
# Stage 2: integrated exclusive-sample composition fit
# =============================================================================

STAGE2_VARIABLE_KEYS: Tuple[str, ...] = (
    "Delta_phi",
    "pTmiss",
    "Emiss2",
)

STAGE2_OUTPUT_DIRNAME = "stage2_integrated_fits"

# Morph priors and hard bounds are intentionally modest.  Their role is to
# absorb small data/MC resolution and centering differences without allowing
# either template to morph into the other physics component.
STAGE2_MORPH_CONFIG = {
    "Delta_phi": {
        "kind": "additive",
        "shift_bound": 0.060,       # rad
        "smear_bound": 0.080,       # rad
        "shift_prior": 0.020,       # rad
        "smear_prior": 0.040,       # rad
    },
    "pTmiss": {
        "kind": "positive_log",
        "shift_bound": 0.25,        # log(x+eps) shift
        "smear_bound": 0.35,        # log-space Gaussian sigma
        "shift_prior": 0.10,
        "smear_prior": 0.15,
    },
    "Emiss2": {
        "kind": "asymmetric_additive",
        "shift_bound": 0.35,        # GeV
        "smear_bound": 0.50,        # GeV, independently left/right
        "shift_prior": 0.15,        # GeV
        "smear_prior": 0.20,        # GeV
    },
}


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
        0.30,
        100,
    ),
    PlotVariable(
        "Delta_phi",
        r"$\Delta\phi$",
        r"$\Delta\phi$ (rad)",
        math.pi - 0.40,
        math.pi + 0.40,
        100,
    ),
    PlotVariable(
        "pTmiss",
        r"$p_{T,\mathrm{miss}}$",
        r"$p_{T,\mathrm{miss}}$ (GeV)",
        0.0,
        1.0,
        100,
    ),
    PlotVariable(
        "Emiss2",
        r"$E_{\mathrm{miss}}$",
        r"$E_{\mathrm{miss}}$ (GeV)",
        -0.10,
        4.0,
        100,
    ),
    PlotVariable(
        "Egamma_tag",
        r"$E_{\gamma,\mathrm{tag}}$",
        r"$E_{\gamma,\mathrm{tag}}$ (GeV)",
        0.4,
        9.5,
        100,
    ),
)

COMPUTED_PLOT_KEYS = frozenset(("Mx2_epgamma", "Mx2_ep"))
ALIASED_PLOT_KEYS = frozenset(("Egamma_tag",))

STORED_SHAPE_BRANCHES: Tuple[str, ...] = tuple(
    variable.branch
    for variable in PLOT_VARIABLES
    if (
        variable.branch not in COMPUTED_PLOT_KEYS
        and variable.branch not in ALIASED_PLOT_KEYS
    )
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
    "detector2",
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


def opening_angle_deg(
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

    Requirements:
        e_p > 2 GeV
        opening angle(e, gamma_tag) > 5 deg
        0.4 <= E_gamma_tag < 9.5 GeV

    The e-gamma opening-angle requirement suppresses radiative events.

    No missing/probe-energy requirement and no predicted-probe angular
    acceptance requirement are imposed in Stage 1.
    """
    e_p = np.asarray(arrays["e_p"], dtype=float)
    g_p = np.asarray(arrays["p2_p"], dtype=float)

    e_theta = to_radians(arrays["e_theta"], angle_unit)
    e_phi = to_radians(arrays["e_phi"], angle_unit)
    g_theta = to_radians(arrays["p2_theta"], angle_unit)
    g_phi = to_radians(arrays["p2_phi"], angle_unit)

    e_px, e_py, e_pz = cartesian_from_spherical(
        e_p,
        e_theta,
        e_phi,
    )
    g_px, g_py, g_pz = cartesian_from_spherical(
        g_p,
        g_theta,
        g_phi,
    )

    theta_egamma_deg = opening_angle_deg(
        e_px,
        e_py,
        e_pz,
        g_px,
        g_py,
        g_pz,
    )

    finite = (
        np.isfinite(e_p)
        & np.isfinite(g_p)
        & np.isfinite(theta_egamma_deg)
    )

    electron = finite & (e_p > ELECTRON_P_MIN_GEV)
    nonradiative = electron & (
        theta_egamma_deg > THETA_EGAMMA_MIN_DEG
    )
    accepted = (
        nonradiative
        & (g_p >= TAG_E_MIN_GEV)
        & (g_p < TAG_E_MAX_GEV)
    )

    counters = {
        "input": int(len(e_p)),
        "finite": int(np.count_nonzero(finite)),
        "electron_p_gt_2": int(np.count_nonzero(electron)),
        "theta_egamma_gt_5deg": int(np.count_nonzero(nonradiative)),
        "tag_0p4_to_9p5_GeV": int(np.count_nonzero(accepted)),
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
    counts_before: Dict[str, Dict[str, np.ndarray]]
    counts_after: Dict[str, Dict[str, np.ndarray]]
    stage2_counts: Dict[str, Dict[str, np.ndarray]]
    angle_unit: str


def empty_histograms() -> Dict[str, np.ndarray]:
    return {
        variable.branch: np.zeros(variable.bins, dtype=np.int64)
        for variable in PLOT_VARIABLES
    }


def empty_region_histograms() -> Dict[str, Dict[str, np.ndarray]]:
    return {
        region_key: empty_histograms()
        for region_key, _, _ in TAG_REGIONS
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

    The selected reconstructed tag photon is split by detector2:
        detector2 == 0 -> FT tag
        detector2 == 1 -> FD tag

    Both tag regions and both pre-/post-exclusivity rows are accumulated in
    this same pass over the ROOT tree.
    """
    counts_before = empty_region_histograms()
    counts_after = empty_region_histograms()
    stage2_counts = empty_region_histograms()

    totals = {
        "input": 0,
        "finite": 0,
        "electron_p_gt_2": 0,
        "theta_egamma_gt_5deg": 0,
        "tag_0p4_to_9p5_GeV": 0,
    }
    for region_key, _, _ in TAG_REGIONS:
        totals[f"{region_key}_minimal"] = 0
        totals[f"{region_key}_after_mx2_epgamma_cut"] = 0
        totals[f"{region_key}_stage2_common"] = 0
    #endfor

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

            detector2 = np.asarray(arrays["detector2"], dtype=int)
            mx2_epgamma = np.asarray(
                missing["Mx2_epgamma"],
                dtype=float,
            )

            entries_read += n
            for key in (
                "input",
                "finite",
                "electron_p_gt_2",
                "theta_egamma_gt_5deg",
                "tag_0p4_to_9p5_GeV",
            ):
                totals[key] += int(chunk_counts[key])
            #endfor

            for region_key, _, detector_value in TAG_REGIONS:
                region_minimal = minimal & (detector2 == detector_value)
                region_after = (
                    region_minimal
                    & np.isfinite(mx2_epgamma)
                    & (
                        np.abs(mx2_epgamma)
                        < MX2_EPGAMMA_ABS_MAX_GEV2
                    )
                )

                totals[f"{region_key}_minimal"] += int(
                    np.count_nonzero(region_minimal)
                )
                totals[
                    f"{region_key}_after_mx2_epgamma_cut"
                ] += int(np.count_nonzero(region_after))

                # Stage 2 uses exactly the SAME events in Delta_phi, pTmiss,
                # and Emiss.  This prevents a variable-range difference from
                # masquerading as a composition difference.
                stage2_common = region_after.copy()
                for stage2_key in STAGE2_VARIABLE_KEYS:
                    stage2_variable = next(
                        variable
                        for variable in PLOT_VARIABLES
                        if variable.branch == stage2_key
                    )
                    stage2_values = np.asarray(
                        arrays[stage2_key],
                        dtype=float,
                    )
                    stage2_common &= (
                        np.isfinite(stage2_values)
                        & (stage2_values >= stage2_variable.low)
                        & (stage2_values < stage2_variable.high)
                    )
                #endfor
                totals[f"{region_key}_stage2_common"] += int(
                    np.count_nonzero(stage2_common)
                )

                for variable in PLOT_VARIABLES:
                    if variable.branch in COMPUTED_PLOT_KEYS:
                        values = np.asarray(
                            missing[variable.branch],
                            dtype=float,
                        )
                    elif variable.branch == "Egamma_tag":
                        # The reconstructed/tag photon energy is p2_p in the
                        # epgamma tree. Photons are massless, so E_gamma=p.
                        values = np.asarray(
                            arrays["p2_p"],
                            dtype=float,
                        )
                    else:
                        # IMPORTANT: stored branches are plotted exactly as
                        # stored in the ROOT tree.
                        values = np.asarray(
                            arrays[variable.branch],
                            dtype=float,
                        )
                    #endif

                    finite_values = np.isfinite(values)

                    before_values = values[
                        region_minimal & finite_values
                    ]
                    after_values = values[
                        region_after & finite_values
                    ]

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

                    counts_before[region_key][
                        variable.branch
                    ] += hist_before.astype(np.int64)
                    counts_after[region_key][
                        variable.branch
                    ] += hist_after.astype(np.int64)

                    if variable.branch in STAGE2_VARIABLE_KEYS:
                        stage2_values = values[
                            stage2_common & finite_values
                        ]
                        hist_stage2, _ = np.histogram(
                            stage2_values,
                            bins=variable.bins,
                            range=(variable.low, variable.high),
                        )
                        stage2_counts[region_key][
                            variable.branch
                        ] += hist_stage2.astype(np.int64)
                    #endif
                #endfor
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
        stage2_counts=stage2_counts,
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
# Stage 2 fit machinery
# =============================================================================

@dataclass
class Stage2VariableFit:
    success: bool
    branch: str
    f_pi0: float = math.nan
    objective: float = math.nan
    deviance: float = math.nan
    ndf: int = 0
    data_counts: Optional[np.ndarray] = None
    model_counts: Optional[np.ndarray] = None
    dvcs_component_counts: Optional[np.ndarray] = None
    pi0_component_counts: Optional[np.ndarray] = None
    morphed_dvcs_shape: Optional[np.ndarray] = None
    morphed_pi0_shape: Optional[np.ndarray] = None
    shared_nuisance: Optional[np.ndarray] = None
    message: str = ""


@dataclass
class Stage2SharedFit:
    success: bool
    f_pi0: float = math.nan
    f_pi0_err: float = math.nan
    objective: float = math.nan
    deviance: float = math.nan
    ndf: int = 0
    variable_results: Optional[Dict[str, Stage2VariableFit]] = None
    individual_fractions: Optional[Dict[str, float]] = None
    message: str = ""


def stage2_variable(branch: str) -> PlotVariable:
    return next(
        variable for variable in PLOT_VARIABLES
        if variable.branch == branch
    )


def normalized_shape(counts: np.ndarray) -> Optional[np.ndarray]:
    shape = np.asarray(counts, dtype=float)
    total = float(np.sum(shape))
    if total <= 0.0 or not math.isfinite(total):
        return None
    #endif
    return np.clip(shape, 0.0, None) / total


def poisson_deviance(
    observed: np.ndarray,
    expected: np.ndarray,
) -> float:
    observed = np.asarray(observed, dtype=float)
    expected = np.clip(np.asarray(expected, dtype=float), 1.0e-12, None)
    positive = observed > 0.0
    terms = expected - observed
    terms[positive] += observed[positive] * np.log(
        observed[positive] / expected[positive]
    )
    return 2.0 * float(np.sum(terms))


def stage2_bin_centers(variable: PlotVariable) -> np.ndarray:
    edges = np.linspace(variable.low, variable.high, variable.bins + 1)
    return 0.5 * (edges[:-1] + edges[1:])


def transform_additive_stage2(
    base_shape: np.ndarray,
    variable: PlotVariable,
    shift: float,
    sigma: float,
) -> Optional[np.ndarray]:
    centers = stage2_bin_centers(variable)
    bin_width = (variable.high - variable.low) / variable.bins

    # Positive shift moves the template toward larger x.
    shifted = np.interp(
        centers - shift,
        centers,
        base_shape,
        left=0.0,
        right=0.0,
    )
    sigma_bins = max(float(sigma) / bin_width, 0.0)
    if sigma_bins > 1.0e-8:
        shifted = gaussian_filter1d(
            shifted,
            sigma=sigma_bins,
            mode="constant",
            cval=0.0,
            truncate=5.0,
        )
    #endif
    shifted = np.clip(shifted, 0.0, None)
    total = float(np.sum(shifted))
    return shifted / total if total > 0.0 else None


def transform_positive_log_stage2(
    base_shape: np.ndarray,
    variable: PlotVariable,
    log_shift: float,
    log_sigma: float,
) -> Optional[np.ndarray]:
    """Morph a positive-definite shape in log(x+epsilon)."""
    centers = stage2_bin_centers(variable)
    bin_width = (variable.high - variable.low) / variable.bins
    epsilon = max(0.5 * bin_width, 1.0e-4)

    source_log = np.log(centers + epsilon) + float(log_shift)
    target_log = np.log(centers + epsilon)

    sigma = max(float(log_sigma), 1.0e-6)
    if log_sigma <= 1.0e-8:
        transformed = np.interp(
            target_log,
            source_log,
            base_shape,
            left=0.0,
            right=0.0,
        )
    else:
        delta = target_log[:, None] - source_log[None, :]
        kernel = np.exp(-0.5 * (delta / sigma) ** 2)
        kernel_sum = np.sum(kernel, axis=0)
        valid = kernel_sum > 0.0
        kernel[:, valid] /= kernel_sum[valid][None, :]
        transformed = kernel @ base_shape
    #endif

    transformed = np.clip(transformed, 0.0, None)
    total = float(np.sum(transformed))
    return transformed / total if total > 0.0 else None


def transform_asymmetric_additive_stage2(
    base_shape: np.ndarray,
    variable: PlotVariable,
    shift: float,
    sigma_left: float,
    sigma_right: float,
) -> Optional[np.ndarray]:
    """
    Additive shift with a split-Gaussian response.

    This is used for E_miss, whose detector/data-MC mismatch can have a
    different width on the negative and positive sides.
    """
    centers = stage2_bin_centers(variable)
    source_means = centers + float(shift)

    delta = centers[:, None] - source_means[None, :]
    widths = np.where(
        delta < 0.0,
        max(float(sigma_left), 1.0e-6),
        max(float(sigma_right), 1.0e-6),
    )
    kernel = np.exp(-0.5 * (delta / widths) ** 2)
    kernel_sum = np.sum(kernel, axis=0)
    valid = kernel_sum > 0.0
    kernel[:, valid] /= kernel_sum[valid][None, :]
    transformed = kernel @ base_shape

    transformed = np.clip(transformed, 0.0, None)
    total = float(np.sum(transformed))
    return transformed / total if total > 0.0 else None


def stage2_nuisance_size(branch: str) -> int:
    return 3 if branch == "Emiss2" else 2


def stage2_single_bounds(branch: str) -> List[Tuple[float, float]]:
    config = STAGE2_MORPH_CONFIG[branch]
    if branch == "Emiss2":
        return [
            (-config["shift_bound"], config["shift_bound"]),
            (0.0, config["smear_bound"]),
            (0.0, config["smear_bound"]),
        ]
    #endif
    return [
        (-config["shift_bound"], config["shift_bound"]),
        (0.0, config["smear_bound"]),
    ]


def stage2_single_start(branch: str) -> np.ndarray:
    if branch == "Emiss2":
        return np.asarray([0.0, 0.05, 0.08], dtype=float)
    #endif
    if branch == "pTmiss":
        return np.asarray([0.0, 0.04], dtype=float)
    #endif
    return np.asarray([0.0, 0.01], dtype=float)


def stage2_transform(
    base_shape: np.ndarray,
    branch: str,
    nuisance: np.ndarray,
) -> Optional[np.ndarray]:
    variable = stage2_variable(branch)
    if branch == "pTmiss":
        return transform_positive_log_stage2(
            base_shape,
            variable,
            float(nuisance[0]),
            float(nuisance[1]),
        )
    #endif
    if branch == "Emiss2":
        return transform_asymmetric_additive_stage2(
            base_shape,
            variable,
            float(nuisance[0]),
            float(nuisance[1]),
            float(nuisance[2]),
        )
    #endif
    return transform_additive_stage2(
        base_shape,
        variable,
        float(nuisance[0]),
        float(nuisance[1]),
    )


def stage2_single_penalty(
    branch: str,
    nuisance: np.ndarray,
) -> float:
    config = STAGE2_MORPH_CONFIG[branch]
    shift = float(nuisance[0])
    value = 0.5 * (shift / config["shift_prior"]) ** 2
    value += 0.5 * (
        float(nuisance[1]) / config["smear_prior"]
    ) ** 2
    if branch == "Emiss2":
        value += 0.5 * (
            float(nuisance[2]) / config["smear_prior"]
        ) ** 2
    #endif
    return value


def stage2_build_model(
    data_counts: np.ndarray,
    dvcs_shape: np.ndarray,
    pi0_shape: np.ndarray,
    branch: str,
    fraction: float,
    nuisance: np.ndarray,
) -> Optional[
    Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]
]:
    nuisance = np.asarray(nuisance, dtype=float)
    if nuisance.size != stage2_nuisance_size(branch):
        return None
    #endif

    # The same detector-response morph is applied to both templates because
    # both are reconstructed e'p'gamma final states.  This removes a large
    # composition/morph degeneracy while preserving the distinct raw physics
    # shapes of dvcsgen and aaogen.
    dvcs_morphed = stage2_transform(
        dvcs_shape,
        branch,
        nuisance,
    )
    pi0_morphed = stage2_transform(
        pi0_shape,
        branch,
        nuisance,
    )
    if dvcs_morphed is None or pi0_morphed is None:
        return None
    #endif

    mixture = (
        (1.0 - fraction) * dvcs_morphed
        + fraction * pi0_morphed
    )
    mixture_sum = float(np.sum(mixture))
    data_total = float(np.sum(data_counts))
    if mixture_sum <= 0.0 or data_total <= 0.0:
        return None
    #endif

    model = data_total * mixture / mixture_sum
    dvcs_component = (
        data_total * (1.0 - fraction) * dvcs_morphed / mixture_sum
    )
    pi0_component = (
        data_total * fraction * pi0_morphed / mixture_sum
    )
    return (
        model,
        dvcs_component,
        pi0_component,
        dvcs_morphed,
        pi0_morphed,
    )


def stage2_fit_shared(
    data_hists: Dict[str, np.ndarray],
    dvcs_hists: Dict[str, np.ndarray],
    pi0_hists: Dict[str, np.ndarray],
    branches: Sequence[str] = STAGE2_VARIABLE_KEYS,
) -> Stage2SharedFit:
    """
    Fit one common pi0 fraction to several projections of the SAME events.

    The fraction is shared across all variables.  Each variable has one
    common data-vs-MC response morph that is applied identically to the DVCS
    and pi0 templates because both are reconstructed e'p'gamma final states.
    The alternating profile procedure mirrors the robust strategy used in the
    exclusivity-selection code, but with this stronger physical constraint.
    """
    if minimize is None or minimize_scalar is None or gaussian_filter1d is None:
        return Stage2SharedFit(
            False,
            message="scipy optimization/ndimage is unavailable",
        )
    #endif

    prepared = {}
    for branch in branches:
        data = np.asarray(data_hists[branch], dtype=float)
        dvcs_shape = normalized_shape(dvcs_hists[branch])
        pi0_shape = normalized_shape(pi0_hists[branch])
        if (
            float(np.sum(data)) <= 0.0
            or dvcs_shape is None
            or pi0_shape is None
        ):
            return Stage2SharedFit(
                False,
                message=f"empty data/template for {branch}",
            )
        #endif
        prepared[branch] = (data, dvcs_shape, pi0_shape)
    #endfor

    def nuisance_bounds(branch: str):
        return stage2_single_bounds(branch)

    def nuisance_start(branch: str):
        return stage2_single_start(branch).copy()

    def variable_objective(
        branch: str,
        fraction: float,
        nuisance: np.ndarray,
        include_penalty: bool,
    ) -> float:
        data, dvcs_shape, pi0_shape = prepared[branch]
        built = stage2_build_model(
            data,
            dvcs_shape,
            pi0_shape,
            branch,
            fraction,
            nuisance,
        )
        if built is None:
            return 1.0e100
        #endif
        value = 0.5 * poisson_deviance(data, built[0])
        if include_penalty:
            value += stage2_single_penalty(
                branch,
                nuisance,
            )
        #endif
        return float(value)

    best_objective = math.inf
    best_fraction = math.nan
    best_nuisances = None

    for initial_fraction in (0.10, 0.30, 0.60, 0.85):
        fraction = float(initial_fraction)
        nuisances = {
            branch: nuisance_start(branch)
            for branch in branches
        }
        previous = math.inf

        for _ in range(12):
            for branch in branches:
                result = minimize(
                    lambda values, b=branch, f=fraction:
                        variable_objective(
                            b,
                            f,
                            values,
                            True,
                        ),
                    nuisances[branch],
                    method="L-BFGS-B",
                    bounds=nuisance_bounds(branch),
                    options={"maxiter": 300, "ftol": 1.0e-10},
                )
                if result.success and np.all(np.isfinite(result.x)):
                    nuisances[branch] = np.asarray(
                        result.x,
                        dtype=float,
                    )
                #endif
            #endfor

            def fraction_objective(candidate: float) -> float:
                return sum(
                    variable_objective(
                        branch,
                        candidate,
                        nuisances[branch],
                        False,
                    )
                    for branch in branches
                )

            f_result = minimize_scalar(
                fraction_objective,
                bounds=(0.0, 1.0),
                method="bounded",
                options={"xatol": 1.0e-5, "maxiter": 200},
            )
            if f_result.success and math.isfinite(float(f_result.x)):
                fraction = float(f_result.x)
            #endif

            current = fraction_objective(fraction)
            if abs(previous - current) <= 1.0e-8 * max(1.0, current):
                break
            #endif
            previous = current
        #endfor

        # Final objective includes the nuisance priors.
        objective = sum(
            variable_objective(
                branch,
                fraction,
                nuisances[branch],
                True,
            )
            for branch in branches
        )
        if objective < best_objective:
            best_objective = float(objective)
            best_fraction = float(fraction)
            best_nuisances = {
                branch: values.copy()
                for branch, values in nuisances.items()
            }
        #endif
    #endfor

    if best_nuisances is None or not math.isfinite(best_fraction):
        return Stage2SharedFit(False, message="shared fit did not converge")
    #endif

    variable_results: Dict[str, Stage2VariableFit] = {}
    total_deviance = 0.0
    total_bins = 0
    nuisance_parameters = 0

    for branch in branches:
        data, dvcs_shape, pi0_shape = prepared[branch]
        nuisance = best_nuisances[branch]
        built = stage2_build_model(
            data,
            dvcs_shape,
            pi0_shape,
            branch,
            best_fraction,
            nuisance,
        )
        if built is None:
            return Stage2SharedFit(
                False,
                message=f"failed to build final {branch} model",
            )
        #endif
        deviance = poisson_deviance(data, built[0])
        total_deviance += deviance
        total_bins += len(data)
        nuisance_parameters += len(nuisance)

        variable_results[branch] = Stage2VariableFit(
            success=True,
            branch=branch,
            f_pi0=best_fraction,
            objective=variable_objective(
                branch,
                best_fraction,
                nuisance,
                True,
            ),
            deviance=deviance,
            ndf=max(1, len(data) - len(nuisance) - 1),
            data_counts=data.copy(),
            model_counts=built[0].copy(),
            dvcs_component_counts=built[1].copy(),
            pi0_component_counts=built[2].copy(),
            morphed_dvcs_shape=built[3].copy(),
            morphed_pi0_shape=built[4].copy(),
            shared_nuisance=nuisance.copy(),
            message="shared-fraction, shared-response integrated fit",
        )
    #endfor

    # Independent one-variable fits are intentionally retained as a diagnostic
    # of model dependence. They do NOT define the nominal shared result.
    individual_fractions = {}
    for branch in branches:
        single = stage2_fit_shared(
            data_hists,
            dvcs_hists,
            pi0_hists,
            branches=(branch,),
        ) if len(branches) > 1 else None
        if single is not None and single.success:
            individual_fractions[branch] = single.f_pi0
        elif len(branches) == 1:
            # Avoid recursion for the single-variable call.
            individual_fractions[branch] = best_fraction
        #endif
    #endfor

    # Curvature estimate of the integrated statistical uncertainty with
    # nuisance parameters held at their final profiled values. This is a
    # development diagnostic; bootstrap/profile uncertainties come later.
    f_error = math.nan
    step = 1.0e-3
    if step < best_fraction < 1.0 - step:
        def fixed_nuisance_objective(candidate: float) -> float:
            return sum(
                variable_objective(
                    branch,
                    candidate,
                    best_nuisances[branch],
                    False,
                )
                for branch in branches
            )
        left = fixed_nuisance_objective(best_fraction - step)
        center = fixed_nuisance_objective(best_fraction)
        right = fixed_nuisance_objective(best_fraction + step)
        curvature = (left - 2.0 * center + right) / (step ** 2)
        if curvature > 0.0:
            f_error = 1.0 / math.sqrt(curvature)
        #endif
    #endif

    return Stage2SharedFit(
        success=True,
        f_pi0=best_fraction,
        f_pi0_err=f_error,
        objective=best_objective,
        deviance=total_deviance,
        ndf=max(
            1,
            total_bins - nuisance_parameters - 1,
        ),
        variable_results=variable_results,
        individual_fractions=individual_fractions,
        message="integrated common-population shared-fraction/shared-response fit",
    )


def run_stage2_internal_self_test() -> None:
    """Catch Stage-2 refactor errors before any expensive ROOT I/O."""
    if minimize is None or minimize_scalar is None or gaussian_filter1d is None:
        raise RuntimeError(
            "Stage 2 requires scipy.optimize and scipy.ndimage."
        )
    #endif

    rng = np.random.default_rng(20260813)
    test_hists = {"data": {}, "dvcs": {}, "pi0": {}}
    true_fraction = 0.42

    for branch in STAGE2_VARIABLE_KEYS:
        variable = stage2_variable(branch)
        x = stage2_bin_centers(variable)

        if branch == "Delta_phi":
            dvcs = np.exp(-0.5 * ((x - math.pi) / 0.045) ** 2)
            pi0 = np.exp(-0.5 * ((x - math.pi) / 0.18) ** 2)
        elif branch == "pTmiss":
            dvcs = np.exp(-x / 0.08)
            pi0 = np.exp(-x / 0.32)
        else:
            dvcs = np.exp(-0.5 * (x / 0.22) ** 2)
            pi0 = np.where(
                x < 0.2,
                np.exp(-0.5 * ((x - 0.15) / 0.45) ** 2),
                np.exp(-(x - 0.2) / 1.2),
            )
        #endif

        dvcs /= np.sum(dvcs)
        pi0 /= np.sum(pi0)
        model = 6000.0 * (
            (1.0 - true_fraction) * dvcs
            + true_fraction * pi0
        )
        data = rng.poisson(model)

        test_hists["data"][branch] = data
        test_hists["dvcs"][branch] = 200000.0 * dvcs
        test_hists["pi0"][branch] = 200000.0 * pi0
    #endfor

    result = stage2_fit_shared(
        test_hists["data"],
        test_hists["dvcs"],
        test_hists["pi0"],
    )
    if (
        not result.success
        or not math.isfinite(result.f_pi0)
        or abs(result.f_pi0 - true_fraction) > 0.12
    ):
        raise RuntimeError(
            "Stage-2 internal fit self-test failed: "
            f"success={result.success}, f_pi0={result.f_pi0}."
        )
    #endif


def format_stage2_nuisance(
    branch: str,
    nuisance: Optional[np.ndarray],
) -> str:
    if nuisance is None:
        return "shared morph unavailable"
    #endif

    values = np.asarray(nuisance, dtype=float)
    if branch == "Delta_phi":
        return (
            rf"shared morph: $\Delta={values[0]:+.4f}$ rad, "
            rf"$\sigma={values[1]:.4f}$ rad"
        )
    #endif
    if branch == "pTmiss":
        return (
            rf"shared log-morph: $\Delta_{{\log}}={values[0]:+.4f}$, "
            rf"$\sigma_{{\log}}={values[1]:.4f}$"
        )
    #endif
    if branch == "Emiss2":
        return (
            rf"shared morph: $\Delta={values[0]:+.4f}$ GeV, "
            rf"$\sigma_L={values[1]:.4f}$ GeV, "
            rf"$\sigma_R={values[2]:.4f}$ GeV"
        )
    #endif
    return "shared morph"
#enddef


def draw_stage2_integrated_canvas(
    period: PeriodConfig,
    region_key: str,
    region_label: str,
    fit: Stage2SharedFit,
    outdir: Path,
    common_counts: Dict[str, int],
) -> Optional[Path]:
    if not fit.success or fit.variable_results is None:
        return None
    #endif

    fig, axes = plt.subplots(
        2,
        len(STAGE2_VARIABLE_KEYS),
        figsize=(15.5, 8.1),
        squeeze=False,
        gridspec_kw={"height_ratios": [3.0, 1.0]},
    )

    for col, branch in enumerate(STAGE2_VARIABLE_KEYS):
        variable = stage2_variable(branch)
        result = fit.variable_results[branch]
        edges = np.linspace(variable.low, variable.high, variable.bins + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])

        ax = axes[0, col]
        data = result.data_counts
        model = np.clip(result.model_counts, 1.0e-12, None)

        # Raw MC shapes, normalized to the data total, are thin dashed curves.
        data_total = float(np.sum(data))
        raw_dvcs = normalized_shape(
            common_counts["dvcs_hists"][branch]
        )
        raw_pi0 = normalized_shape(
            common_counts["pi0_hists"][branch]
        )
        if raw_dvcs is not None:
            ax.step(
                edges[:-1],
                data_total * raw_dvcs,
                where="post",
                linestyle="--",
                linewidth=0.8,
                alpha=0.55,
                color="tab:blue",
                label="raw DVCS shape",
            )
        #endif
        if raw_pi0 is not None:
            ax.step(
                edges[:-1],
                data_total * raw_pi0,
                where="post",
                linestyle="--",
                linewidth=0.8,
                alpha=0.55,
                color="red",
                label=r"raw $e\pi^0$ shape",
            )
        #endif

        ax.step(
            edges[:-1],
            result.dvcs_component_counts,
            where="post",
            linewidth=1.35,
            color="tab:blue",
            label="fitted DVCS component",
        )
        ax.step(
            edges[:-1],
            result.pi0_component_counts,
            where="post",
            linewidth=1.35,
            color="red",
            label=r"fitted $e\pi^0$ component",
        )
        ax.step(
            edges[:-1],
            model,
            where="post",
            linewidth=1.6,
            color="green",
            label="total fit",
        )
        ax.errorbar(
            centers,
            data,
            yerr=np.sqrt(np.maximum(data, 1.0)),
            fmt="o",
            markersize=2.6,
            linewidth=0.6,
            color="black",
            label="data",
            zorder=5,
        )

        independent = (
            fit.individual_fractions.get(branch, math.nan)
            if fit.individual_fractions
            else math.nan
        )
        nuisance_text = format_stage2_nuisance(
            branch,
            result.shared_nuisance,
        )
        ax.set_title(
            f"{variable.title}\n"
            rf"shared $f_{{\pi^0}}={fit.f_pi0:.3f}$; "
            rf"single-var $f_{{\pi^0}}={independent:.3f}$; "
            rf"$D/ndf={result.deviance:.1f}/{result.ndf}$"
            "\n"
            f"{nuisance_text}",
            fontsize=9.4,
        )
        ax.set_ylabel("entries / bin")
        ax.set_xlim(variable.low, variable.high)
        ax.grid(alpha=0.18)

        pull_ax = axes[1, col]
        pull = (data - model) / np.sqrt(model)
        pull_ax.axhline(0.0, color="0.35", linewidth=0.9)
        pull_ax.axhline(+2.0, color="0.6", linewidth=0.7, linestyle="--")
        pull_ax.axhline(-2.0, color="0.6", linewidth=0.7, linestyle="--")
        pull_ax.plot(
            centers,
            pull,
            "o",
            markersize=2.3,
            color="black",
        )
        pull_ax.set_ylim(-6.0, 6.0)
        pull_ax.set_xlim(variable.low, variable.high)
        pull_ax.set_xlabel(variable.xlabel)
        pull_ax.set_ylabel("pull")
        pull_ax.grid(alpha=0.18)
    #endfor

    handles, labels = axes[0, 0].get_legend_handles_labels()

    fpi0_text = rf"$f_{{\pi^0}}={fit.f_pi0:.4f}"
    if math.isfinite(fit.f_pi0_err):
        fpi0_text += rf"\pm{fit.f_pi0_err:.4f}"
    #endif
    fpi0_text += "$"

    fig.suptitle(
        f"Stage 2 integrated fit: {period.label} — {region_label}\n"
        rf"$|M_X^2(ep\gamma)|<{MX2_EPGAMMA_ABS_MAX_GEV2:.3f}$ GeV$^2$; "
        f"common entries: data={common_counts['data_n']:,}, "
        f"DVCS MC={common_counts['dvcs_n']:,}, "
        rf"$e\pi^0$ MC={common_counts['pi0_n']:,}; "
        f"shared {fpi0_text}",
        fontsize=12.5,
        y=0.987,
    )

    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.885),
        ncol=6,
        frameon=False,
        fontsize=8.5,
    )

    fig.subplots_adjust(
        left=0.065,
        right=0.992,
        bottom=0.085,
        top=0.795,
        wspace=0.24,
        hspace=0.25,
    )

    outpath = (
        outdir
        / f"stage2_integrated_fits_{period.key}_{region_key.lower()}.png"
    )
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
    return outpath

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
    region_key: str,
    region_label: str,
) -> Path:
    """
    Draw one compact 2x5 canvas for one reconstructed TAG-photon detector:
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
                rf"After $|M_X^2(ep\gamma)|"
                rf"<{MX2_EPGAMMA_ABS_MAX_GEV2:.3f}$ GeV$^2$",
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

            data_counts = getattr(data, count_attr)[
                region_key
            ][variable.branch]
            dvcs_counts = getattr(dvcs, count_attr)[
                region_key
            ][variable.branch]
            aaogen_counts = getattr(aaogen, count_attr)[
                region_key
            ][variable.branch]

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

            positive_data = data_y > 0.0
            ax.errorbar(
                centers[positive_data],
                data_y[positive_data],
                yerr=data_err[positive_data],
                fmt="o",
                markersize=2.2,
                linewidth=0.5,
                capsize=0,
                color="black",
                label=r"$e'p'\gamma$ data",
                zorder=5,
            )

            ax.set_xlim(variable.low, variable.high)

            # Every Stage-1 shape is displayed on a logarithmic y axis.
            # Zero-content data bins are omitted above; zero-content step
            # bins are naturally ignored by the log transform.
            ax.set_yscale("log")

            positive = np.concatenate(
                (
                    data_y[data_y > 0.0],
                    dvcs_y[dvcs_y > 0.0],
                    aaogen_y[aaogen_y > 0.0],
                )
            )
            if len(positive) > 0:
                ymin = max(
                    float(np.min(positive)) * 0.5,
                    1.0e-6,
                )
                ymax = float(np.max(positive)) * 2.0
                ax.set_ylim(ymin, ymax)
            #endif

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

    n_data_before = data.selection_counts[
        f"{region_key}_minimal"
    ]
    n_dvcs_before = dvcs.selection_counts[
        f"{region_key}_minimal"
    ]
    n_aaogen_before = aaogen.selection_counts[
        f"{region_key}_minimal"
    ]

    n_data_after = data.selection_counts[
        f"{region_key}_after_mx2_epgamma_cut"
    ]
    n_dvcs_after = dvcs.selection_counts[
        f"{region_key}_after_mx2_epgamma_cut"
    ]
    n_aaogen_after = aaogen.selection_counts[
        f"{region_key}_after_mx2_epgamma_cut"
    ]

    fig.suptitle(
        f"Stage 1 shape comparison: {period.label} — {region_label}\n"
        r"minimal selection: $E_e>2$ GeV, "
        r"$\theta_{e\gamma}>5^\circ$, "
        r"$0.4\leq E_{\gamma,\mathrm{tag}}<9.5$ GeV"
        "\n"
        f"before cut: data={n_data_before:,}, "
        f"DVCS MC={n_dvcs_before:,}, "
        rf"$e\pi^0$ MC={n_aaogen_before:,}; "
        f"after cut: data={n_data_after:,}, "
        f"DVCS MC={n_dvcs_after:,}, "
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

    fig.subplots_adjust(
        left=0.055,
        right=0.992,
        bottom=0.085,
        top=0.825,
        wspace=0.30,
        hspace=0.38,
    )

    outpath = (
        outdir
        / f"stage1_shape_comparison_{period.key}_{region_key.lower()}.png"
    )
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
            f"FT={result.selection_counts['FT_minimal']:,}"
            f"->{result.selection_counts['FT_after_mx2_epgamma_cut']:,}; "
            f"FD={result.selection_counts['FD_minimal']:,}"
            f"->{result.selection_counts['FD_after_mx2_epgamma_cut']:,}."
        )
    #endfor

    canvases = {}
    for region_key, region_label, _ in TAG_REGIONS:
        canvases[region_key] = str(
            draw_period_canvas(
                period,
                results,
                outdir,
                region_key,
                region_label,
            )
        )
    #endfor

    # Stage 2: integrated composition fit.  Reuse the Stage-2 histograms
    # accumulated in the same ROOT streaming pass above.
    stage2_outdir = outdir.parent / STAGE2_OUTPUT_DIRNAME
    stage2_outdir.mkdir(parents=True, exist_ok=True)
    stage2_results = {}

    for region_key, region_label, _ in TAG_REGIONS:
        data_hists = {
            branch: results["data"].stage2_counts[region_key][branch]
            for branch in STAGE2_VARIABLE_KEYS
        }
        dvcs_hists = {
            branch: results["dvcsgen"].stage2_counts[region_key][branch]
            for branch in STAGE2_VARIABLE_KEYS
        }
        pi0_hists = {
            branch: results["aaogen"].stage2_counts[region_key][branch]
            for branch in STAGE2_VARIABLE_KEYS
        }

        fit = stage2_fit_shared(
            data_hists,
            dvcs_hists,
            pi0_hists,
        )

        common = {
            "data_n": results["data"].selection_counts[
                f"{region_key}_stage2_common"
            ],
            "dvcs_n": results["dvcsgen"].selection_counts[
                f"{region_key}_stage2_common"
            ],
            "pi0_n": results["aaogen"].selection_counts[
                f"{region_key}_stage2_common"
            ],
            "dvcs_hists": dvcs_hists,
            "pi0_hists": pi0_hists,
        }

        canvas = draw_stage2_integrated_canvas(
            period,
            region_key,
            region_label,
            fit,
            stage2_outdir,
            common,
        )

        stage2_results[region_key] = {
            "success": bool(fit.success),
            "f_pi0": float(fit.f_pi0),
            "f_pi0_err_curvature": float(fit.f_pi0_err),
            "deviance": float(fit.deviance),
            "ndf": int(fit.ndf),
            "individual_fractions": (
                {
                    key: float(value)
                    for key, value in fit.individual_fractions.items()
                }
                if fit.individual_fractions
                else {}
            ),
            "shared_morph_nuisance": (
                {
                    branch: [
                        float(value)
                        for value in variable_result.shared_nuisance
                    ]
                    for branch, variable_result
                    in fit.variable_results.items()
                    if variable_result.shared_nuisance is not None
                }
                if fit.variable_results
                else {}
            ),
            "common_population": {
                "data": int(common["data_n"]),
                "dvcs_mc": int(common["dvcs_n"]),
                "pi0_mc": int(common["pi0_n"]),
            },
            "canvas": str(canvas) if canvas is not None else None,
            "message": fit.message,
        }

        log(
            f"{period.label} {region_label}: Stage 2 "
            f"{'OK' if fit.success else 'FAILED'}; "
            f"f_pi0={fit.f_pi0:.4f}"
            + (
                f" +/- {fit.f_pi0_err:.4f}"
                if math.isfinite(fit.f_pi0_err)
                else ""
            )
        )
    #endfor

    summary = {
        "period": period.key,
        "label": period.label,
        "canvases": canvases,
        "stage2_integrated_fits": stage2_results,
        "samples": {
            key: {
                "entries_read": result.entries_read,
                "FT_minimal": result.selection_counts["FT_minimal"],
                "FT_after_mx2_cut": result.selection_counts[
                    "FT_after_mx2_epgamma_cut"
                ],
                "FD_minimal": result.selection_counts["FD_minimal"],
                "FD_after_mx2_cut": result.selection_counts[
                    "FD_after_mx2_epgamma_cut"
                ],
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
            "Stage 1 shape comparison + Stage 2 integrated composition fit."
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

    log("Stages 1+2: shape comparison followed by integrated DVCS/pi0 composition fits.")
    log("Minimal selection: e_p>2 GeV, theta_egamma>5 deg, 0.4<=E_tag<9.5 GeV; no probe-energy or probe-angle cut.")

    # Exercise the full Stage-2 optimizer before any large ROOT read. This is
    # specifically intended to catch refactor/runtime errors early.
    run_stage2_internal_self_test()
    log("Internal Stage-2 fit self-test passed.")

    # Complete ROOT schema preflight BEFORE any large read.
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
            "theta_egamma_min_deg": THETA_EGAMMA_MIN_DEG,
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

    stage2_summary = {
        "stage": "stage2_integrated_fits",
        "selection": {
            "stage1_minimal": {
                "electron_p_min_GeV": ELECTRON_P_MIN_GEV,
                "theta_egamma_min_deg": THETA_EGAMMA_MIN_DEG,
                "tag_energy_GeV": [TAG_E_MIN_GEV, TAG_E_MAX_GEV],
            },
            "exclusivity": {
                "mx2_epgamma_abs_max_GeV2": MX2_EPGAMMA_ABS_MAX_GEV2,
            },
            "fraction_variables": list(STAGE2_VARIABLE_KEYS),
            "common_population_across_fraction_variables": True,
        },
        "periods": [
            {
                "period": period.key,
                "results": summaries[period.key].get(
                    "stage2_integrated_fits",
                    {},
                ),
            }
            for period in periods
            if period.key in summaries
        ],
    }

    stage2_outdir = Path(args.output_dir).parent / STAGE2_OUTPUT_DIRNAME
    stage2_outdir.mkdir(parents=True, exist_ok=True)
    with (
        stage2_outdir / "stage2_summary.json"
    ).open("w", encoding="utf-8") as stream:
        json.dump(stage2_summary, stream, indent=2)
    #endwith

    log(
        "Done. Stage-1 outputs are in "
        f"{Path(args.output_dir)}; Stage-2 outputs are in "
        f"{stage2_outdir}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
