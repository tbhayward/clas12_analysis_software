#!/usr/bin/env python3
"""
derive_photon_efficiency_scale_factors.py

CLAS12 RGA photon-efficiency study -- clean refactor.

CURRENT IMPLEMENTED SCOPE
=========================
Stage 1: Shape Comparison

This script intentionally does ONE thing only:

    Compare the shapes of six stored epgamma observables in:
      1. loose nSidis epgamma DATA,
      2. dvcsgen epgamma MC,
      3. aaogen epgamma MC,

after applying the same minimal tag/probe support definition to all three
samples.

There are NO fits, NO eppi0 input files, NO efficiency calculation, NO
numerator association, NO bootstrap, and NO systematic extraction in this
version.

The reconstructed epgamma photon is the tag.  The missing four-vector

    P_probe^pred = P_beam + P_target - P_e - P_p - P_gamma(tag)

defines the predicted probe photon used only for the common support selection:

    e_p > 2 GeV
    theta_ep > 5 deg
    2 <= E_tag < 9.5 GeV
    0.4 <= E_probe^pred < 9.5 GeV
    predicted probe in:
        FT: 2.4 <= theta_probe^pred < 5 deg
        OR
        FD: 6 <= theta_probe^pred < 35 deg

No exclusivity cut is imposed.

The six SHAPE variables are read directly from the epgamma ROOT trees:

    Delta_phi
    pTmiss
    theta
    theta_gamma_gamma
    Emiss2       (despite the branch name, this is E_miss, not E_miss^2)
    z

One 2x3 unit-normalized shape-comparison canvas is written per run period.

Default output:
    output/photon_efficiency/stage1_shape_comparison/

Typical commands:
    # Fa18 Inb only
    python derive_photon_efficiency_scale_factors.py --period fa18_inb

    # all configured periods
    python derive_photon_efficiency_scale_factors.py --workers 2

Development convention:
    --max-entries 0 means read all entries.
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
        "Delta_phi",
        r"$\Delta\phi$",
        r"$\Delta\phi$ (rad)",
        2.80,
        3.45,
        80,
    ),
    PlotVariable(
        "pTmiss",
        r"$p_{T,\mathrm{miss}}$",
        r"$p_{T,\mathrm{miss}}$ (GeV)",
        0.0,
        1.0,
        80,
    ),
    PlotVariable(
        "theta",
        r"$\theta$",
        r"$\theta$ (rad)",
        0.0,
        math.pi,
        80,
    ),
    PlotVariable(
        "theta_gamma_gamma",
        r"$\theta_{\gamma\gamma}$",
        r"$\theta_{\gamma\gamma}$ (rad)",
        0.0,
        3.0,
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
    PlotVariable(
        "z",
        r"$z$",
        r"$z$",
        0.0,
        1.6,
        80,
    ),
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

SHAPE_BRANCHES: Tuple[str, ...] = tuple(
    variable.branch for variable in PLOT_VARIABLES
)

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


def predicted_probe_kinematics(
    arrays: Dict[str, np.ndarray],
    beam_energy_GeV: float,
    angle_unit: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return (E_probe^pred, theta_probe^pred_deg).

    Electron energy is computed with the electron-mass-negligible
    approximation E_e = |p_e|, which is more than sufficient on this scale.
    The reconstructed photon is massless, so E_tag = p2_p.
    """
    e_p = np.asarray(arrays["e_p"], dtype=float)
    p_p = np.asarray(arrays["p1_p"], dtype=float)
    g_E = np.asarray(arrays["p2_p"], dtype=float)

    e_theta = to_radians(arrays["e_theta"], angle_unit)
    e_phi = to_radians(arrays["e_phi"], angle_unit)

    p_theta = to_radians(arrays["p1_theta"], angle_unit)
    p_phi = to_radians(arrays["p1_phi"], angle_unit)

    g_theta = to_radians(arrays["p2_theta"], angle_unit)
    g_phi = to_radians(arrays["p2_phi"], angle_unit)

    e_px, e_py, e_pz = cartesian_from_spherical(
        e_p, e_theta, e_phi
    )
    p_px, p_py, p_pz = cartesian_from_spherical(
        p_p, p_theta, p_phi
    )
    g_px, g_py, g_pz = cartesian_from_spherical(
        g_E, g_theta, g_phi
    )

    proton_E = np.sqrt(
        np.maximum(p_p * p_p + PROTON_MASS_GEV * PROTON_MASS_GEV, 0.0)
    )

    pred_E = (
        beam_energy_GeV
        + PROTON_MASS_GEV
        - e_p
        - proton_E
        - g_E
    )

    pred_px = -e_px - p_px - g_px
    pred_py = -e_py - p_py - g_py
    pred_pz = beam_energy_GeV - e_pz - p_pz - g_pz

    pred_p = np.sqrt(
        pred_px * pred_px + pred_py * pred_py + pred_pz * pred_pz
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        cosine = pred_pz / pred_p
    #endwith
    cosine = np.clip(cosine, -1.0, 1.0)
    pred_theta_deg = np.degrees(np.arccos(cosine))

    return pred_E, pred_theta_deg


def common_support_mask(
    arrays: Dict[str, np.ndarray],
    beam_energy_GeV: float,
    angle_unit: str,
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Apply only the already-established tag/probe support definition.

    No exclusivity selection is made in Stage 1.
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

    pred_E, pred_theta_deg = predicted_probe_kinematics(
        arrays,
        beam_energy_GeV,
        angle_unit,
    )

    finite = (
        np.isfinite(e_p)
        & np.isfinite(p_p)
        & np.isfinite(g_E)
        & np.isfinite(theta_ep_deg)
        & np.isfinite(pred_E)
        & np.isfinite(pred_theta_deg)
    )

    electron = finite & (e_p > ELECTRON_P_MIN_GEV)

    opening = electron & (theta_ep_deg > THETA_EP_MIN_DEG)

    tag = (
        opening
        & (g_E >= TAG_E_MIN_GEV)
        & (g_E < TAG_E_MAX_GEV)
    )

    probe_energy = (
        tag
        & (pred_E >= PROBE_E_MIN_GEV)
        & (pred_E < PROBE_E_MAX_GEV)
    )

    probe_ft = (
        probe_energy
        & (pred_theta_deg >= FT_THETA_MIN_DEG)
        & (pred_theta_deg < FT_THETA_MAX_DEG)
    )

    probe_fd = (
        probe_energy
        & (pred_theta_deg >= FD_THETA_MIN_DEG)
        & (pred_theta_deg < FD_THETA_MAX_DEG)
    )

    accepted = probe_ft | probe_fd

    counters = {
        "input": int(len(e_p)),
        "finite": int(np.count_nonzero(finite)),
        "electron_p_gt_2": int(np.count_nonzero(electron)),
        "theta_ep_gt_5deg": int(np.count_nonzero(opening)),
        "tag_2_to_9p5_GeV": int(np.count_nonzero(tag)),
        "probe_0p4_to_9p5_GeV": int(np.count_nonzero(probe_energy)),
        "predicted_probe_FT": int(np.count_nonzero(probe_ft)),
        "predicted_probe_FD": int(np.count_nonzero(probe_fd)),
        "accepted_FT_or_FD": int(np.count_nonzero(accepted)),
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
    counts: Dict[str, np.ndarray]
    underflow: Dict[str, int]
    overflow: Dict[str, int]
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
    Stream only the 15 branches required by Stage 1.

    No multi-million-entry event arrays are retained.  This is intentionally
    histogram-only so the runtime and memory footprint remain small.
    """
    counts = empty_histograms()
    underflow = {
        variable.branch: 0 for variable in PLOT_VARIABLES
    }
    overflow = {
        variable.branch: 0 for variable in PLOT_VARIABLES
    }

    totals = {
        "input": 0,
        "finite": 0,
        "electron_p_gt_2": 0,
        "theta_ep_gt_5deg": 0,
        "tag_2_to_9p5_GeV": 0,
        "probe_0p4_to_9p5_GeV": 0,
        "predicted_probe_FT": 0,
        "predicted_probe_FD": 0,
        "accepted_FT_or_FD": 0,
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

            selection, chunk_counts = common_support_mask(
                arrays,
                beam_energy_GeV,
                angle_unit,
            )

            entries_read += n
            for key in totals:
                totals[key] += int(chunk_counts[key])
            #endfor

            for variable in PLOT_VARIABLES:
                values = np.asarray(
                    arrays[variable.branch],
                    dtype=float,
                )
                selected_values = values[
                    selection & np.isfinite(values)
                ]

                hist, _ = np.histogram(
                    selected_values,
                    bins=variable.bins,
                    range=(variable.low, variable.high),
                )
                counts[variable.branch] += hist.astype(np.int64)

                underflow[variable.branch] += int(
                    np.count_nonzero(selected_values < variable.low)
                )
                overflow[variable.branch] += int(
                    np.count_nonzero(selected_values >= variable.high)
                )
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
        counts=counts,
        underflow=underflow,
        overflow=overflow,
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
    Match the supplied exclusivity-shape aesthetic:
      - blue step = dvcsgen
      - red step = aaogen
      - black points = data
      - unit-normalized entries/bin
      - one shared legend above the 2x3 canvas
    """
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(16.0, 9.2),
    )
    axes_flat = axes.ravel()

    data = results["data"]
    dvcs = results["dvcsgen"]
    aaogen = results["aaogen"]

    for ax, variable in zip(axes_flat, PLOT_VARIABLES):
        edges = np.linspace(
            variable.low,
            variable.high,
            variable.bins + 1,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        data_y, data_err = normalized_histogram(
            data.counts[variable.branch]
        )
        dvcs_y, _ = normalized_histogram(
            dvcs.counts[variable.branch]
        )
        aaogen_y, _ = normalized_histogram(
            aaogen.counts[variable.branch]
        )

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
            markersize=2.4,
            linewidth=0.55,
            capsize=0,
            color="black",
            label=r"$e'p'\gamma$ data",
            zorder=5,
        )

        ax.set_xlim(variable.low, variable.high)
        ax.set_xlabel(variable.xlabel)
        ax.set_ylabel("unit-normalized entries / bin")
        ax.grid(alpha=0.18)
    #endfor

    handles, labels = axes_flat[0].get_legend_handles_labels()

    n_data = data.selection_counts["accepted_FT_or_FD"]
    n_dvcs = dvcs.selection_counts["accepted_FT_or_FD"]
    n_aaogen = aaogen.selection_counts["accepted_FT_or_FD"]

    fig.suptitle(
        f"Stage 1 shape comparison: {period.label}\n"
        "minimal tag/probe support only; no exclusivity cuts\n"
        f"selected entries: data={n_data:,}, "
        f"DVCS MC={n_dvcs:,}, "
        rf"$e\pi^0$ MC as $ep\gamma$={n_aaogen:,}",
        fontsize=14,
        y=0.985,
    )

    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.905),
        ncol=3,
        frameon=False,
        fontsize=10,
    )

    fig.tight_layout(
        rect=(0.02, 0.03, 0.99, 0.84)
    )

    outpath = outdir / f"stage1_shape_comparison_{period.key}.png"
    fig.savefig(outpath, dpi=180)
    plt.close(fig)
    return outpath


# =============================================================================
# Per-period execution
# =============================================================================

def period_output_dir(
    base_output_dir: str,
    period: PeriodConfig,
) -> Path:
    return Path(base_output_dir) / period.key


def process_period(
    period: PeriodConfig,
    args_dict: Dict[str, object],
) -> Dict[str, object]:
    t0 = time.perf_counter()

    outdir = period_output_dir(
        str(args_dict["output_dir"]),
        period,
    )
    outdir.mkdir(parents=True, exist_ok=True)

    sample_specs = (
        (
            "data",
            r"$e'p'\gamma$ data",
            period.nsidis_epgamma_data,
        ),
        (
            "dvcsgen",
            "DVCS MC",
            period.dvcsgen_epgamma_mc,
        ),
        (
            "aaogen",
            r"$e\pi^0$ MC as $ep\gamma$",
            period.aaogen_epgamma_mc,
        ),
    )

    results: Dict[str, HistogramResult] = {}

    for sample_key, sample_label, file_path in sample_specs:
        log(
            f"{period.label}: streaming {sample_key} Stage-1 shapes."
        )

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
            f"{period.label}: {sample_key}: "
            f"read {result.entries_read:,}; "
            f"selected {result.selection_counts['accepted_FT_or_FD']:,}."
        )
    #endfor

    canvas = draw_period_canvas(
        period,
        results,
        outdir,
    )

    summary = {
        "stage": "stage1_shape_comparison",
        "period": asdict(period),
        "selection": {
            "electron_p_min_GeV": ELECTRON_P_MIN_GEV,
            "theta_ep_min_deg": THETA_EP_MIN_DEG,
            "tag_energy_GeV": [
                TAG_E_MIN_GEV,
                TAG_E_MAX_GEV,
            ],
            "probe_energy_GeV": [
                PROBE_E_MIN_GEV,
                PROBE_E_MAX_GEV,
            ],
            "FT_theta_deg": [
                FT_THETA_MIN_DEG,
                FT_THETA_MAX_DEG,
            ],
            "FD_theta_deg": [
                FD_THETA_MIN_DEG,
                FD_THETA_MAX_DEG,
            ],
            "exclusivity_cuts": "none",
        },
        "samples": {},
        "output_canvas": str(canvas),
        "wall_time_s": float(time.perf_counter() - t0),
    }

    for key, result in results.items():
        summary["samples"][key] = {
            "file": result.file_path,
            "entries_read": result.entries_read,
            "angle_unit": result.angle_unit,
            "selection_counts": result.selection_counts,
            "plot_underflow": result.underflow,
            "plot_overflow": result.overflow,
        }
    #endfor

    with (outdir / "stage1_summary.json").open(
        "w",
        encoding="utf-8",
    ) as stream:
        json.dump(summary, stream, indent=2)
    #endwith

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

    log(
        "Stage 1 only: six epgamma shape comparisons; "
        "no fits or efficiency calculation."
    )
    log(
        "Minimal support: "
        "e_p>2 GeV, theta_ep>5 deg, "
        "2<=E_tag<9.5 GeV, 0.4<=E_probe^pred<9.5 GeV, "
        "probe in FT or FD angular acceptance."
    )

    # Complete preflight BEFORE any large read.
    preflight_report = preflight_periods(
        periods,
        args.tree_name,
    )

    with (
        Path(args.output_dir) / "stage1_preflight.json"
    ).open("w", encoding="utf-8") as stream:
        json.dump(preflight_report, stream, indent=2)
    #endwith

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

    ordered_summary = {
        "stage": "stage1_shape_comparison",
        "periods": [
            summaries[period.key]
            for period in periods
            if period.key in summaries
        ],
    }

    with (
        Path(args.output_dir) / "stage1_all_periods_summary.json"
    ).open("w", encoding="utf-8") as stream:
        json.dump(ordered_summary, stream, indent=2)
    #endwith

    log(
        "Done. Stage-1 outputs are in "
        f"{Path(args.output_dir)}."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
#endif
