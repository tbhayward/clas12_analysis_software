#!/usr/bin/env python3
"""
derive_piplus_momentum_corrections_v3.py

First diagnostic pass for RGC pi+ momentum corrections using inclusive
e pi+ X data and the missing-neutron peak in

    Mx2(e pi+) = (p_beam + p_target - p_e' - p_pi+)^2.

This stage does NOT yet derive or apply a pi+ momentum correction. It performs
the prerequisite diagnostics in two reconstruction stages:

    1. no momentum corrections;
    2. finalized electron momentum correction applied event by event.

The purpose is to establish:

    * the neutron-peak fit model;
    * available statistics in FD and CD pion samples;
    * sensible theta binning;
    * period, sector, detector-region, and run stability;
    * the size of the electron-correction effect before fitting pi+ corrections.

FD and CD pions are treated separately.

FD:
    * six sectors from pion phi;
    * equal-population pion-theta cells inside each period and sector.

CD:
    * integrated over azimuth at this stage;
    * equal-population pion-theta cells inside each period.

Three signal models are attempted in each Mx2 cell:

    * Gaussian + quadratic background;
    * Crystal Ball + quadratic background;
    * common-mean double Gaussian + quadratic background.

The lowest-AICc successful model is selected. Fit-quality flags are advisory.
All individual fits include a pull panel.

Expected electron-correction JSON:
    output/electron_diagnostics/json/electron_correction_models.json

The script accepts branch aliases and prints the resolved mapping before
processing. Defaults assume common RGC naming conventions but may be overridden
from the command line.

Primary outputs:
    output/piplus_diagnostics/
        csv/
        json/
        plots/
            uncorrected/
            electron_corrected/
            comparison/
            event_counts/

Dependencies:
    numpy
    pandas
    scipy
    matplotlib
    uproot
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import json
import math
import shutil
import sys
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import uproot


# -----------------------------------------------------------------------------
# Physics constants
# -----------------------------------------------------------------------------

PROTON_MASS_GEV = 0.9382720813
NEUTRON_MASS_GEV = 0.9395654133
NEUTRON_MASS2_GEV2 = NEUTRON_MASS_GEV**2
PION_MASS_GEV = 0.13957039
ELECTRON_MASS_GEV = 0.00051099895
RAD2DEG = 180.0 / math.pi


DEFAULT_INPUTS = {
    "su22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/"
        "rgc_su22_inb_NH3_epi+.root"
    ),
    "fa22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/"
        "rgc_fa22_inb_NH3_epi+.root"
    ),
    "sp23": Path(
        "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/"
        "rgc_sp23_inb_NH3_epi+.root"
    ),
}


@dataclass(frozen=True)
class CalibrationPeriod:
    key: str
    label: str
    source_key: str
    run_min: int
    run_max: int


CALIBRATION_PERIODS = (
    CalibrationPeriod(
        key="su22",
        label="Su22",
        source_key="su22",
        run_min=16042,
        run_max=16788,
    ),
    CalibrationPeriod(
        key="fa22_solenoid_minus",
        label="Fa22, solenoid -1",
        source_key="fa22",
        run_min=16843,
        run_max=17183,
    ),
    CalibrationPeriod(
        key="fa22_solenoid_plus",
        label="Fa22, solenoid +1",
        source_key="fa22",
        run_min=17185,
        run_max=17408,
    ),
    CalibrationPeriod(
        key="sp23",
        label="Sp23",
        source_key="sp23",
        run_min=17477,
        run_max=17811,
    ),
)


BRANCH_ALIASES = {
    "runnum": ("runnum", "run", "run_number"),
    "evnum": ("evnum", "event", "event_number"),
    "e_p": ("e_p", "electron_p"),
    "e_theta": ("e_theta", "electron_theta"),
    "e_phi": ("e_phi", "electron_phi"),
    "pip_p": (
        "pip_p",
        "pi_p",
        "piplus_p",
        "pion_p",
        "h_p",
    ),
    "pip_theta": (
        "pip_theta",
        "pi_theta",
        "piplus_theta",
        "pion_theta",
        "h_theta",
    ),
    "pip_phi": (
        "pip_phi",
        "pi_phi",
        "piplus_phi",
        "pion_phi",
        "h_phi",
    ),
    "pip_detector": (
        "pip_detector",
        "pi_detector",
        "piplus_detector",
        "pion_detector",
        "detector",
    ),
}


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

def beam_energy_from_run(runnum: np.ndarray) -> np.ndarray:
    """Assign the nominal beam energy event by event from run number."""
    runnum = np.asarray(runnum, dtype=int)
    energy = np.full(runnum.shape, np.nan, dtype=float)

    energy[(runnum >= 16042) & (runnum <= 17065)] = 10.5473
    energy[(runnum >= 17067) & (runnum <= 17716)] = 10.5563
    energy[(runnum >= 17717) & (runnum <= 17811)] = 10.5593
    return energy


def wrap_phi_deg(phi_deg: np.ndarray | float) -> np.ndarray:
    """Wrap azimuth to [0, 360) degrees."""
    return np.mod(np.asarray(phi_deg, dtype=float), 360.0)


def fd_sector_from_phi_rad(phi_rad: np.ndarray) -> np.ndarray:
    """Assign the six standard CLAS12 FD sectors from pion azimuth."""
    phi_deg = wrap_phi_deg(np.asarray(phi_rad, dtype=float) * RAD2DEG)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)

    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def fd_local_phi_deg(
    phi_rad: np.ndarray,
    sector: np.ndarray,
) -> np.ndarray:
    """Rotate each FD sector onto sector 1 and return local phi on [0, 60)."""
    phi_deg = wrap_phi_deg(np.asarray(phi_rad, dtype=float) * RAD2DEG)
    sector = np.asarray(sector, dtype=int)
    lower_edge = np.full(phi_deg.shape, np.nan, dtype=float)

    lower_edge[sector == 1] = 330.0
    lower_edge[sector == 2] = 30.0
    lower_edge[sector == 3] = 90.0
    lower_edge[sector == 4] = 150.0
    lower_edge[sector == 5] = 210.0
    lower_edge[sector == 6] = 270.0

    local_phi = np.mod(phi_deg - lower_edge, 360.0)
    local_phi[(local_phi < 0.0) | (local_phi >= 60.0)] = np.nan
    return local_phi


def find_tree(file_path: Path, requested_tree: str) -> str:
    """Resolve the requested ROOT tree, accepting a cycle suffix."""
    with uproot.open(file_path) as root_file:
        keys = root_file.keys()
        key_names = {str(key).split(";")[0]: str(key) for key in keys}

        if requested_tree in key_names:
            return key_names[requested_tree]
        # endif

        if requested_tree in keys:
            return requested_tree
        # endif

        tree_like = []
        for key in keys:
            try:
                obj = root_file[key]
                if hasattr(obj, "arrays"):
                    tree_like.append(str(key))
                # endif
            except Exception:
                continue
            # endtry
        # endfor

    if len(tree_like) == 1:
        return tree_like[0]
    # endif

    raise KeyError(
        f"Tree '{requested_tree}' not found in {file_path}. "
        f"Tree-like objects: {tree_like}"
    )


def parse_branch_override(text: str) -> tuple[str, str]:
    """Parse logical_name=ROOT_branch."""
    if "=" not in text:
        raise argparse.ArgumentTypeError(
            "Branch overrides must have the form logical_name=ROOT_branch."
        )
    # endif

    logical, physical = (piece.strip() for piece in text.split("=", 1))
    if logical not in BRANCH_ALIASES:
        raise argparse.ArgumentTypeError(
            f"Unknown logical branch name '{logical}'. "
            f"Allowed: {sorted(BRANCH_ALIASES)}"
        )
    # endif
    if not physical:
        raise argparse.ArgumentTypeError("ROOT branch name cannot be empty.")
    # endif
    return logical, physical


def resolve_branches(
    file_path: Path,
    tree_name: str,
    overrides: dict[str, str],
) -> dict[str, str]:
    """Resolve required logical branches against available ROOT branches."""
    with uproot.open(file_path) as root_file:
        available = set(root_file[tree_name].keys())
    # endwith

    resolved: dict[str, str] = {}
    missing: list[str] = []

    for logical_name, aliases in BRANCH_ALIASES.items():
        if logical_name in overrides:
            candidate = overrides[logical_name]
            if candidate not in available:
                raise KeyError(
                    f"Requested branch override {logical_name}={candidate} "
                    f"is absent from {file_path}:{tree_name}."
                )
            # endif
            resolved[logical_name] = candidate
            continue
        # endif

        match = next((alias for alias in aliases if alias in available), None)
        if match is None:
            missing.append(logical_name)
        else:
            resolved[logical_name] = match
        # endif
    # endfor

    if missing:
        raise KeyError(
            f"Could not resolve required logical branches {missing} in "
            f"{file_path}:{tree_name}. Available branches include "
            f"{sorted(available)[:80]}"
        )
    # endif

    return resolved


def sanitize_json_value(value: Any) -> Any:
    """Convert NumPy values and non-finite floats for JSON serialization."""
    if isinstance(value, dict):
        return {
            str(key): sanitize_json_value(item)
            for key, item in value.items()
        }
    # endif
    if isinstance(value, (list, tuple)):
        return [sanitize_json_value(item) for item in value]
    # endif
    if isinstance(value, np.ndarray):
        return [sanitize_json_value(item) for item in value.tolist()]
    # endif
    if isinstance(value, np.generic):
        return sanitize_json_value(value.item())
    # endif
    if isinstance(value, float) and not math.isfinite(value):
        return None
    # endif
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write indented JSON safely."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as output:
        json.dump(
            sanitize_json_value(payload),
            output,
            indent=2,
            sort_keys=True,
        )
        output.write("\n")
    # endwith


# -----------------------------------------------------------------------------
# Electron correction handling
# -----------------------------------------------------------------------------

def load_electron_models(path: Path) -> dict[tuple[str, int], dict[str, Any]]:
    """Load selected electron correction models keyed by period and sector."""
    if not path.is_file():
        raise FileNotFoundError(
            f"Electron correction JSON does not exist: {path}"
        )
    # endif

    with path.open("r", encoding="utf-8") as input_file:
        payload = json.load(input_file)
    # endwith

    records = payload.get("models", [])
    model_map: dict[tuple[str, int], dict[str, Any]] = {}

    for record in records:
        if not record.get("success", False):
            continue
        # endif

        period_key = str(record["period_key"])
        sector = int(record["sector"])
        coefficients = np.asarray(
            record.get("selected_coefficients_ascending", []),
            dtype=float,
        )
        if coefficients.size == 0:
            continue
        # endif

        model_map[(period_key, sector)] = record
    # endfor

    expected = {
        (period.key, sector)
        for period in CALIBRATION_PERIODS
        for sector in range(1, 7)
    }
    missing = sorted(expected - set(model_map))
    if missing:
        raise RuntimeError(
            "Electron correction JSON is missing successful models for "
            f"{missing}"
        )
    # endif

    return model_map


def evaluate_electron_model(
    theta_deg: np.ndarray,
    model_record: dict[str, Any],
) -> np.ndarray:
    """Evaluate one selected electron model with boundary clipping."""
    theta = np.asarray(theta_deg, dtype=float)
    theta_clipped = np.clip(
        theta,
        float(model_record["theta_valid_min_deg"]),
        float(model_record["theta_valid_max_deg"]),
    )
    coefficients = np.asarray(
        model_record["selected_coefficients_ascending"],
        dtype=float,
    )

    correction = np.zeros_like(theta_clipped, dtype=float)
    for power, coefficient in enumerate(coefficients):
        correction += coefficient * theta_clipped**power
    # endfor
    return correction


def apply_electron_correction(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    model_map: dict[tuple[str, int], dict[str, Any]],
) -> pd.DataFrame:
    """Apply finalized electron correction while leaving pion momentum fixed."""
    corrected = frame.copy()
    e_sector = fd_sector_from_phi_rad(
        corrected["e_phi_rad"].to_numpy()
    )
    corrected["e_sector"] = e_sector
    corrected["electron_delta_p_over_p"] = np.nan

    for sector in range(1, 7):
        sector_mask = corrected["e_sector"] == sector
        model = model_map[(period.key, sector)]
        corrected.loc[
            sector_mask,
            "electron_delta_p_over_p",
        ] = evaluate_electron_model(
            corrected.loc[sector_mask, "e_theta_deg"].to_numpy(),
            model,
        )
    # endfor

    invalid = corrected["electron_delta_p_over_p"].isna()
    if invalid.any():
        corrected = corrected.loc[~invalid].copy()
    # endif

    corrected["e_p_electron_corrected_gev"] = (
        corrected["e_p_gev"]
        * (1.0 + corrected["electron_delta_p_over_p"])
    )
    return corrected


# -----------------------------------------------------------------------------
# Four-vectors and Mx2
# -----------------------------------------------------------------------------

def momentum_components(
    momentum: np.ndarray,
    theta_rad: np.ndarray,
    phi_rad: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert spherical momentum coordinates to Cartesian components."""
    transverse = momentum * np.sin(theta_rad)
    px = transverse * np.cos(phi_rad)
    py = transverse * np.sin(phi_rad)
    pz = momentum * np.cos(theta_rad)
    return px, py, pz


def calculate_mx2_epi(
    beam_energy_gev: np.ndarray,
    e_p_gev: np.ndarray,
    e_theta_rad: np.ndarray,
    e_phi_rad: np.ndarray,
    pip_p_gev: np.ndarray,
    pip_theta_rad: np.ndarray,
    pip_phi_rad: np.ndarray,
) -> np.ndarray:
    """Calculate Mx2(e pi+) with target proton at rest."""
    beam_energy = np.asarray(beam_energy_gev, dtype=float)

    e_px, e_py, e_pz = momentum_components(
        e_p_gev,
        e_theta_rad,
        e_phi_rad,
    )
    pip_px, pip_py, pip_pz = momentum_components(
        pip_p_gev,
        pip_theta_rad,
        pip_phi_rad,
    )

    e_energy = np.sqrt(e_p_gev**2 + ELECTRON_MASS_GEV**2)
    pip_energy = np.sqrt(pip_p_gev**2 + PION_MASS_GEV**2)

    missing_energy = beam_energy + PROTON_MASS_GEV - e_energy - pip_energy
    missing_px = -e_px - pip_px
    missing_py = -e_py - pip_py
    missing_pz = beam_energy - e_pz - pip_pz

    return (
        missing_energy**2
        - missing_px**2
        - missing_py**2
        - missing_pz**2
    )


# -----------------------------------------------------------------------------
# Input reading
# -----------------------------------------------------------------------------

def read_period_events(
    file_path: Path,
    tree_name: str,
    branch_map: dict[str, str],
    period: CalibrationPeriod,
    step_size: str,
    mx2_preselection_min: float,
    mx2_preselection_max: float,
    e_theta_min_deg: float,
    e_theta_max_deg: float,
    pip_theta_min_deg: float,
    pip_theta_max_deg: float,
) -> pd.DataFrame:
    """Read one period in chunks and retain finite e pi+ events."""
    physical_branches = list(dict.fromkeys(branch_map.values()))
    frames: list[pd.DataFrame] = []

    source = f"{file_path}:{tree_name}"
    for arrays in uproot.iterate(
        source,
        expressions=physical_branches,
        step_size=step_size,
        library="np",
    ):
        values = {
            logical: np.asarray(arrays[physical])
            for logical, physical in branch_map.items()
        }

        runnum = values["runnum"].astype(np.int64)
        period_mask = (
            (runnum >= period.run_min)
            & (runnum <= period.run_max)
        )
        if not np.any(period_mask):
            continue
        # endif

        runnum = runnum[period_mask]
        evnum = values["evnum"][period_mask].astype(np.int64)
        e_p = values["e_p"][period_mask].astype(float)
        e_theta = values["e_theta"][period_mask].astype(float)
        e_phi = values["e_phi"][period_mask].astype(float)
        pip_p = values["pip_p"][period_mask].astype(float)
        pip_theta = values["pip_theta"][period_mask].astype(float)
        pip_phi = values["pip_phi"][period_mask].astype(float)
        pip_detector = values["pip_detector"][period_mask].astype(int)

        e_theta_deg = e_theta * RAD2DEG
        pip_theta_deg = pip_theta * RAD2DEG
        beam_energy = beam_energy_from_run(runnum)

        finite = (
            np.isfinite(e_p)
            & np.isfinite(e_theta)
            & np.isfinite(e_phi)
            & np.isfinite(pip_p)
            & np.isfinite(pip_theta)
            & np.isfinite(pip_phi)
            & np.isfinite(beam_energy)
        )
        selected = (
            finite
            & (e_p > 0.0)
            & (pip_p > 0.0)
            & (e_theta_deg >= e_theta_min_deg)
            & (e_theta_deg <= e_theta_max_deg)
            & (pip_theta_deg >= pip_theta_min_deg)
            & (pip_theta_deg <= pip_theta_max_deg)
            & np.isin(pip_detector, [1, 2])
        )

        if not np.any(selected):
            continue
        # endif

        runnum = runnum[selected]
        evnum = evnum[selected]
        e_p = e_p[selected]
        e_theta = e_theta[selected]
        e_phi = e_phi[selected]
        pip_p = pip_p[selected]
        pip_theta = pip_theta[selected]
        pip_phi = pip_phi[selected]
        pip_detector = pip_detector[selected]
        e_theta_deg = e_theta_deg[selected]
        pip_theta_deg = pip_theta_deg[selected]
        beam_energy = beam_energy[selected]

        mx2_uncorrected = calculate_mx2_epi(
            beam_energy,
            e_p,
            e_theta,
            e_phi,
            pip_p,
            pip_theta,
            pip_phi,
        )

        mx2_mask = (
            np.isfinite(mx2_uncorrected)
            & (mx2_uncorrected >= mx2_preselection_min)
            & (mx2_uncorrected <= mx2_preselection_max)
        )
        if not np.any(mx2_mask):
            continue
        # endif

        pion_sector = fd_sector_from_phi_rad(pip_phi)
        pion_local_phi = fd_local_phi_deg(pip_phi, pion_sector)

        frame = pd.DataFrame(
            {
                "runnum": runnum[mx2_mask],
                "evnum": evnum[mx2_mask],
                "beam_energy_gev": beam_energy[mx2_mask],
                "e_p_gev": e_p[mx2_mask],
                "e_theta_rad": e_theta[mx2_mask],
                "e_theta_deg": e_theta_deg[mx2_mask],
                "e_phi_rad": e_phi[mx2_mask],
                "pip_p_gev": pip_p[mx2_mask],
                "pip_theta_rad": pip_theta[mx2_mask],
                "pip_theta_deg": pip_theta_deg[mx2_mask],
                "pip_phi_rad": pip_phi[mx2_mask],
                "pip_detector": pip_detector[mx2_mask],
                "pion_region": np.where(
                    pip_detector[mx2_mask] == 1,
                    "FD",
                    "CD",
                ),
                "pion_sector": np.where(
                    pip_detector[mx2_mask] == 1,
                    pion_sector[mx2_mask],
                    0,
                ),
                "pion_local_phi_deg": np.where(
                    pip_detector[mx2_mask] == 1,
                    pion_local_phi[mx2_mask],
                    np.nan,
                ),
                "mx2_uncorrected_gev2": mx2_uncorrected[mx2_mask],
            }
        )
        frames.append(frame)
    # endfor

    if not frames:
        return pd.DataFrame()
    # endif

    return pd.concat(frames, ignore_index=True)


# -----------------------------------------------------------------------------
# Peak models
# -----------------------------------------------------------------------------

def quadratic_background(
    x: np.ndarray,
    c0: float,
    c1: float,
    c2: float,
) -> np.ndarray:
    """Quadratic background around the neutron mass squared."""
    dx = x - NEUTRON_MASS2_GEV2
    return c0 + c1 * dx + c2 * dx**2


def gaussian_signal(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma: float,
) -> np.ndarray:
    """Gaussian signal specified by peak amplitude."""
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def gaussian_plus_quadratic(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma: float,
    c0: float,
    c1: float,
    c2: float,
) -> np.ndarray:
    return (
        gaussian_signal(x, amplitude, mean, sigma)
        + quadratic_background(x, c0, c1, c2)
    )


def crystal_ball_signal(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma: float,
    alpha: float,
    n_tail: float,
) -> np.ndarray:
    """
    Left-tailed Crystal Ball signal.

    The function is continuous at t=-alpha and uses peak amplitude rather than
    unit normalization. This is sufficient for histogram-shape fitting.
    """
    sigma = max(float(sigma), 1.0e-9)
    alpha_abs = max(abs(float(alpha)), 1.0e-6)
    n_tail = max(float(n_tail), 1.01)
    t = (x - mean) / sigma

    gaussian = np.exp(-0.5 * t**2)
    a = (n_tail / alpha_abs) ** n_tail * np.exp(
        -0.5 * alpha_abs**2
    )
    b = n_tail / alpha_abs - alpha_abs
    tail = a * np.power(np.maximum(b - t, 1.0e-12), -n_tail)
    return amplitude * np.where(t > -alpha_abs, gaussian, tail)


def crystal_ball_plus_quadratic(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma: float,
    alpha: float,
    n_tail: float,
    c0: float,
    c1: float,
    c2: float,
) -> np.ndarray:
    return (
        crystal_ball_signal(
            x,
            amplitude,
            mean,
            sigma,
            alpha,
            n_tail,
        )
        + quadratic_background(x, c0, c1, c2)
    )


def double_gaussian_signal(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma_core: float,
    broad_fraction: float,
    sigma_ratio: float,
) -> np.ndarray:
    """Common-mean core plus broad Gaussian component."""
    broad_fraction = np.clip(broad_fraction, 0.0, 1.0)
    sigma_ratio = max(sigma_ratio, 1.0)
    core = gaussian_signal(
        x,
        amplitude * (1.0 - broad_fraction),
        mean,
        sigma_core,
    )
    broad = gaussian_signal(
        x,
        amplitude * broad_fraction,
        mean,
        sigma_core * sigma_ratio,
    )
    return core + broad


def double_gaussian_plus_quadratic(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma_core: float,
    broad_fraction: float,
    sigma_ratio: float,
    c0: float,
    c1: float,
    c2: float,
) -> np.ndarray:
    return (
        double_gaussian_signal(
            x,
            amplitude,
            mean,
            sigma_core,
            broad_fraction,
            sigma_ratio,
        )
        + quadratic_background(x, c0, c1, c2)
    )


MODEL_DEFINITIONS = {
    "gaussian": {
        "function": gaussian_plus_quadratic,
        "parameter_count": 6,
    },
    "crystal_ball": {
        "function": crystal_ball_plus_quadratic,
        "parameter_count": 8,
    },
    "double_gaussian": {
        "function": double_gaussian_plus_quadratic,
        "parameter_count": 8,
    },
}


def model_initialization(
    model_name: str,
    fit_x: np.ndarray,
    fit_y: np.ndarray,
    fit_range: tuple[float, float],
) -> tuple[list[float], tuple[list[float], list[float]]]:
    """Return initial values and bounds for one fit model."""
    baseline = max(float(np.percentile(fit_y, 20.0)), 0.0)
    amplitude = max(float(np.max(fit_y) - baseline), 1.0)

    peak_window = (
        (fit_x >= max(fit_range[0], NEUTRON_MASS2_GEV2 - 0.16))
        & (fit_x <= min(fit_range[1], NEUTRON_MASS2_GEV2 + 0.16))
    )
    if np.any(peak_window):
        peak_x = fit_x[peak_window]
        peak_y = fit_y[peak_window]
        mean_guess = float(peak_x[np.argmax(peak_y)])
    else:
        mean_guess = NEUTRON_MASS2_GEV2
    # endif

    mean_low = max(fit_range[0], NEUTRON_MASS2_GEV2 - 0.18)
    mean_high = min(fit_range[1], NEUTRON_MASS2_GEV2 + 0.18)
    mean_guess = float(np.clip(mean_guess, mean_low, mean_high))

    background_scale = max(float(np.max(fit_y)), 1.0)
    common_background_initial = [baseline, 0.0, 0.0]
    common_background_low = [
        -0.5 * background_scale,
        -20.0 * background_scale,
        -100.0 * background_scale,
    ]
    common_background_high = [
        2.0 * background_scale,
        20.0 * background_scale,
        100.0 * background_scale,
    ]

    if model_name == "gaussian":
        initial = [
            amplitude,
            mean_guess,
            0.045,
            *common_background_initial,
        ]
        lower = [
            0.0,
            mean_low,
            0.008,
            *common_background_low,
        ]
        upper = [
            5.0 * background_scale,
            mean_high,
            0.18,
            *common_background_high,
        ]
    elif model_name == "crystal_ball":
        initial = [
            amplitude,
            mean_guess,
            0.040,
            1.5,
            3.0,
            *common_background_initial,
        ]
        lower = [
            0.0,
            mean_low,
            0.008,
            0.5,
            1.05,
            *common_background_low,
        ]
        upper = [
            5.0 * background_scale,
            mean_high,
            0.18,
            5.0,
            30.0,
            *common_background_high,
        ]
    elif model_name == "double_gaussian":
        initial = [
            amplitude,
            mean_guess,
            0.030,
            0.25,
            2.5,
            *common_background_initial,
        ]
        lower = [
            0.0,
            mean_low,
            0.006,
            0.0,
            1.1,
            *common_background_low,
        ]
        upper = [
            5.0 * background_scale,
            mean_high,
            0.12,
            0.8,
            6.0,
            *common_background_high,
        ]
    else:
        raise ValueError(f"Unknown model name: {model_name}")
    # endif

    return initial, (lower, upper)


def signal_component(
    model_name: str,
    x: np.ndarray,
    parameters: np.ndarray,
) -> np.ndarray:
    """Evaluate only the selected signal component."""
    if model_name == "gaussian":
        return gaussian_signal(x, *parameters[:3])
    # endif
    if model_name == "crystal_ball":
        return crystal_ball_signal(x, *parameters[:5])
    # endif
    if model_name == "double_gaussian":
        return double_gaussian_signal(x, *parameters[:5])
    # endif
    raise ValueError(model_name)


def background_component(
    model_name: str,
    x: np.ndarray,
    parameters: np.ndarray,
) -> np.ndarray:
    """Evaluate only the quadratic background component."""
    if model_name == "gaussian":
        return quadratic_background(x, *parameters[3:6])
    # endif
    return quadratic_background(x, *parameters[5:8])


def fit_one_model(
    model_name: str,
    centers: np.ndarray,
    counts: np.ndarray,
    fit_mask: np.ndarray,
    fit_range: tuple[float, float],
) -> dict[str, Any]:
    """Fit one candidate signal model."""
    fit_x = centers[fit_mask]
    fit_y = counts[fit_mask]
    fit_error = np.sqrt(np.maximum(fit_y, 1.0))

    model_function = MODEL_DEFINITIONS[model_name]["function"]
    initial, bounds = model_initialization(
        model_name,
        fit_x,
        fit_y,
        fit_range,
    )

    try:
        parameters, covariance = curve_fit(
            model_function,
            fit_x,
            fit_y,
            p0=initial,
            bounds=bounds,
            sigma=fit_error,
            absolute_sigma=True,
            maxfev=100000,
        )
    except (
        RuntimeError,
        ValueError,
        np.linalg.LinAlgError,
        FloatingPointError,
    ) as exc:
        return {
            "model": model_name,
            "success": False,
            "status": f"fit_failed:{type(exc).__name__}",
        }
    # endtry

    model_at_centers = model_function(fit_x, *parameters)
    residual = fit_y - model_at_centers
    chi2 = float(np.sum((residual / fit_error) ** 2))
    parameter_count = len(parameters)
    n_points = len(fit_x)
    ndf = n_points - parameter_count
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan

    if n_points > parameter_count + 1:
        aicc = (
            chi2
            + 2.0 * parameter_count
            + (
                2.0
                * parameter_count
                * (parameter_count + 1)
                / (n_points - parameter_count - 1)
            )
        )
    else:
        aicc = math.inf
    # endif

    errors = np.sqrt(
        np.maximum(np.diag(covariance), 0.0)
    )
    mean = float(parameters[1])
    mean_error = float(errors[1])

    mean_low = bounds[0][1]
    mean_high = bounds[1][1]
    mean_at_bound = (
        abs(mean - mean_low) < 0.002
        or abs(mean - mean_high) < 0.002
    )

    signal_at_centers = signal_component(
        model_name,
        fit_x,
        parameters,
    )
    signal_peak = float(np.max(signal_at_centers))
    background_at_peak = float(
        background_component(
            model_name,
            np.asarray([mean]),
            parameters,
        )[0]
    )
    peak_significance_proxy = (
        signal_peak / math.sqrt(max(signal_peak + abs(background_at_peak), 1.0))
    )

    review_reasons = []
    if not np.isfinite(chi2_ndf) or chi2_ndf > 3.0:
        review_reasons.append("poor_chi2")
    # endif
    if mean_error > 0.025:
        review_reasons.append("large_mean_error")
    # endif
    if mean_at_bound:
        review_reasons.append("mean_at_bound")
    # endif
    if peak_significance_proxy < 5.0:
        review_reasons.append("weak_peak")
    # endif

    return {
        "model": model_name,
        "success": True,
        "status": (
            "success"
            if not review_reasons
            else "success_flagged_for_review"
        ),
        "review_reasons": review_reasons,
        "parameters": parameters.tolist(),
        "covariance": covariance.tolist(),
        "mean_gev2": mean,
        "mean_error_gev2": mean_error,
        "chi2": chi2,
        "ndf": int(ndf),
        "chi2_ndf": chi2_ndf,
        "aicc": float(aicc),
        "peak_significance_proxy": peak_significance_proxy,
        "fit_x": fit_x,
        "fit_y": fit_y,
        "fit_error": fit_error,
        "fit_model": model_at_centers,
        "pulls": residual / fit_error,
    }


def fit_mx2_peak(
    values: np.ndarray,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Fit the missing-neutron peak with all candidate models."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    counts, edges = np.histogram(
        values,
        bins=histogram_bins,
        range=histogram_range,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    fit_mask = (
        (centers >= fit_range[0])
        & (centers <= fit_range[1])
    )

    empty_diagnostics = {
        "counts": counts,
        "edges": edges,
        "centers": centers,
        "fit_mask": fit_mask,
        "fit_x": np.asarray([], dtype=float),
        "fit_y": np.asarray([], dtype=float),
        "fit_error": np.asarray([], dtype=float),
        "fit_model": np.asarray([], dtype=float),
        "fit_signal": np.asarray([], dtype=float),
        "fit_background": np.asarray([], dtype=float),
        "pulls": np.asarray([], dtype=float),
    }

    if values.size < minimum_events:
        return (
            {
                "success": False,
                "status": "insufficient_events",
                "n_events": int(values.size),
                "selected_model": None,
                "candidate_results": {},
            },
            empty_diagnostics,
        )
    # endif

    candidate_results = {
        model_name: fit_one_model(
            model_name,
            centers,
            counts,
            fit_mask,
            fit_range,
        )
        for model_name in MODEL_DEFINITIONS
    }
    successful = [
        result
        for result in candidate_results.values()
        if result.get("success", False)
    ]

    if not successful:
        return (
            {
                "success": False,
                "status": "all_models_failed",
                "n_events": int(values.size),
                "selected_model": None,
                "candidate_results": candidate_results,
            },
            empty_diagnostics,
        )
    # endif

    selected = min(
        successful,
        key=lambda result: (
            result["aicc"]
            if np.isfinite(result["aicc"])
            else math.inf
        ),
    )
    model_name = selected["model"]
    parameters = np.asarray(selected["parameters"], dtype=float)
    fit_x = selected["fit_x"]

    diagnostics = {
        "counts": counts,
        "edges": edges,
        "centers": centers,
        "fit_mask": fit_mask,
        "fit_x": fit_x,
        "fit_y": selected["fit_y"],
        "fit_error": selected["fit_error"],
        "fit_model": selected["fit_model"],
        "fit_signal": signal_component(
            model_name,
            fit_x,
            parameters,
        ),
        "fit_background": background_component(
            model_name,
            fit_x,
            parameters,
        ),
        "pulls": selected["pulls"],
    }

    summary = {
        "success": True,
        "status": selected["status"],
        "review_reasons": selected["review_reasons"],
        "n_events": int(values.size),
        "selected_model": model_name,
        "mean_gev2": selected["mean_gev2"],
        "mean_error_gev2": selected["mean_error_gev2"],
        "chi2": selected["chi2"],
        "ndf": selected["ndf"],
        "chi2_ndf": selected["chi2_ndf"],
        "aicc": selected["aicc"],
        "peak_significance_proxy": selected["peak_significance_proxy"],
        "candidate_results": {
            name: {
                key: value
                for key, value in result.items()
                if key not in {
                    "fit_x",
                    "fit_y",
                    "fit_error",
                    "fit_model",
                    "pulls",
                }
            }
            for name, result in candidate_results.items()
        },
    }
    return summary, diagnostics


# -----------------------------------------------------------------------------
# Calibration cells
# -----------------------------------------------------------------------------

def equal_population_edges(
    values: np.ndarray,
    bin_count: int,
) -> np.ndarray:
    """Build strictly increasing equal-population edges."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    if values.size == 0:
        return np.linspace(0.0, 1.0, bin_count + 1)
    # endif

    edges = np.quantile(
        values,
        np.linspace(0.0, 1.0, bin_count + 1),
    ).astype(float)

    for index in range(1, len(edges)):
        if edges[index] <= edges[index - 1]:
            edges[index] = edges[index - 1] + 1.0e-6
        # endif
    # endfor
    return edges


def cell_record(
    period: CalibrationPeriod,
    stage: str,
    region: str,
    sector: int,
    cell_type: str,
    theta_bin_index: int | None,
    theta_low_deg: float,
    theta_high_deg: float,
    cell_frame: pd.DataFrame,
    fit_result: dict[str, Any],
) -> dict[str, Any]:
    """Build one flat fit record."""
    return {
        "period_key": period.key,
        "period_label": period.label,
        "stage": stage,
        "pion_region": region,
        "pion_sector": sector,
        "cell_type": cell_type,
        "theta_bin_index": theta_bin_index,
        "theta_low_deg": theta_low_deg,
        "theta_high_deg": theta_high_deg,
        "mean_theta_deg": (
            float(cell_frame["pip_theta_deg"].mean())
            if not cell_frame.empty
            else math.nan
        ),
        "mean_pion_momentum_gev": (
            float(cell_frame["pip_p_gev"].mean())
            if not cell_frame.empty
            else math.nan
        ),
        "n_events": int(len(cell_frame)),
        "success": bool(fit_result["success"]),
        "status": fit_result["status"],
        "review_reasons": "|".join(
            fit_result.get("review_reasons", [])
        ),
        "selected_model": fit_result.get("selected_model"),
        "peak_mean_gev2": fit_result.get("mean_gev2", math.nan),
        "peak_mean_error_gev2": fit_result.get(
            "mean_error_gev2",
            math.nan,
        ),
        "peak_residual_gev2": (
            fit_result.get("mean_gev2", math.nan)
            - NEUTRON_MASS2_GEV2
        ),
        "chi2": fit_result.get("chi2", math.nan),
        "ndf": fit_result.get("ndf", 0),
        "chi2_ndf": fit_result.get("chi2_ndf", math.nan),
        "aicc": fit_result.get("aicc", math.nan),
        "peak_significance_proxy": fit_result.get(
            "peak_significance_proxy",
            math.nan,
        ),
        "candidate_results_json": json.dumps(
            sanitize_json_value(
                fit_result.get("candidate_results", {})
            ),
            sort_keys=True,
        ),
    }


def stage_column(stage: str) -> str:
    """Return the Mx2 column corresponding to a reconstruction stage."""
    if stage == "uncorrected":
        return "mx2_uncorrected_gev2"
    # endif
    if stage == "electron_corrected":
        return "mx2_electron_corrected_gev2"
    # endif
    raise ValueError(stage)


def fit_region_cells(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    stage: str,
    region: str,
    fd_theta_bins: int,
    cd_theta_bins: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    individual_fit_directory: Path,
    save_individual_fits: bool,
) -> list[dict[str, Any]]:
    """Fit integrated and theta-binned cells for one detector region."""
    mx2_column = stage_column(stage)
    region_frame = frame.loc[frame["pion_region"] == region].copy()
    records: list[dict[str, Any]] = []

    sectors = range(1, 7) if region == "FD" else [0]
    theta_bin_count = (
        fd_theta_bins if region == "FD" else cd_theta_bins
    )

    for sector in sectors:
        if region == "FD":
            sample = region_frame.loc[
                region_frame["pion_sector"] == sector
            ].copy()
        else:
            sample = region_frame.copy()
        # endif

        theta_edges = equal_population_edges(
            sample["pip_theta_deg"].to_numpy()
            if not sample.empty
            else np.asarray([], dtype=float),
            theta_bin_count,
        )

        integrated_fit, integrated_diagnostics = fit_mx2_peak(
            sample[mx2_column].to_numpy()
            if not sample.empty
            else np.asarray([], dtype=float),
            histogram_range,
            fit_range,
            histogram_bins,
            minimum_events_integrated,
        )
        records.append(
            cell_record(
                period,
                stage,
                region,
                sector,
                "integrated",
                None,
                float(theta_edges[0]),
                float(theta_edges[-1]),
                sample,
                integrated_fit,
            )
        )

        if save_individual_fits:
            save_individual_peak_fit(
                individual_fit_directory
                / (
                    f"{period.key}_{region.lower()}_sector{sector}_"
                    f"integrated_{stage}.png"
                ),
                integrated_diagnostics,
                integrated_fit,
                (
                    f"{period.label}, {region}"
                    + (
                        f", sector {sector}"
                        if region == "FD"
                        else ""
                    )
                    + f", integrated, {stage.replace('_', ' ')}"
                ),
            )
        # endif

        for theta_index in range(len(theta_edges) - 1):
            theta_low = float(theta_edges[theta_index])
            theta_high = float(theta_edges[theta_index + 1])

            if theta_index == len(theta_edges) - 2:
                mask = (
                    (sample["pip_theta_deg"] >= theta_low)
                    & (sample["pip_theta_deg"] <= theta_high)
                )
            else:
                mask = (
                    (sample["pip_theta_deg"] >= theta_low)
                    & (sample["pip_theta_deg"] < theta_high)
                )
            # endif

            cell = sample.loc[mask].copy()
            fit_result, diagnostics = fit_mx2_peak(
                cell[mx2_column].to_numpy(),
                histogram_range,
                fit_range,
                histogram_bins,
                minimum_events_cell,
            )
            records.append(
                cell_record(
                    period,
                    stage,
                    region,
                    sector,
                    "theta",
                    theta_index,
                    theta_low,
                    theta_high,
                    cell,
                    fit_result,
                )
            )

            if save_individual_fits:
                save_individual_peak_fit(
                    individual_fit_directory
                    / (
                        f"{period.key}_{region.lower()}_sector{sector}_"
                        f"theta{theta_index:02d}_{stage}.png"
                    ),
                    diagnostics,
                    fit_result,
                    (
                        f"{period.label}, {region}"
                        + (
                            f", sector {sector}"
                            if region == "FD"
                            else ""
                        )
                        + f", theta cell {theta_index + 1}/"
                        f"{len(theta_edges) - 1}, "
                        f"{stage.replace('_', ' ')}"
                    ),
                )
            # endif
        # endfor
    # endfor

    return records


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def save_individual_peak_fit(
    output_path: Path,
    diagnostics: dict[str, np.ndarray],
    fit_result: dict[str, Any],
    title: str,
) -> None:
    """Save one Mx2 fit with a lower pull panel."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, (axis, pull_axis) = plt.subplots(
        2,
        1,
        figsize=(7.2, 6.4),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0]},
        constrained_layout=True,
    )

    centers = diagnostics["centers"]
    counts = diagnostics["counts"]
    errors = np.sqrt(np.maximum(counts, 1.0))

    if centers.size:
        axis.errorbar(
            centers,
            counts,
            yerr=errors,
            fmt=".",
            markersize=3,
            linewidth=0.7,
            label="Data",
        )
    # endif

    if fit_result["success"]:
        fit_x = diagnostics["fit_x"]
        axis.plot(
            fit_x,
            diagnostics["fit_model"],
            linewidth=1.5,
            label="Total fit",
        )
        axis.plot(
            fit_x,
            diagnostics["fit_background"],
            linestyle="--",
            linewidth=1.2,
            label="Background",
        )
        axis.plot(
            fit_x,
            diagnostics["fit_signal"]
            + diagnostics["fit_background"],
            linestyle=":",
            linewidth=1.2,
            label="Signal + background",
        )
        axis.axvline(
            fit_result["mean_gev2"],
            linestyle="--",
            linewidth=1.0,
            label=(
                rf"$\mu_{{Mx^2}}$={fit_result['mean_gev2']:.4f}"
            ),
        )
        pull_axis.plot(
            fit_x,
            diagnostics["pulls"],
            ".",
            markersize=4,
        )
    # endif

    axis.axvline(
        NEUTRON_MASS2_GEV2,
        color="black",
        linestyle="-.",
        linewidth=1.0,
        label=r"$m_n^2$",
    )
    axis.set_ylabel("Counts")
    axis.set_title(
        title
        + "\n"
        + (
            f"{fit_result['selected_model']}, "
            f"{fit_result['status']}, "
            f"chi2/ndf={fit_result.get('chi2_ndf', math.nan):.2f}"
            if fit_result["success"]
            else fit_result["status"]
        )
    )
    axis.legend(fontsize=8)

    pull_axis.axhline(0.0, color="black", linewidth=0.8)
    pull_axis.axhline(+3.0, color="black", linestyle="--", linewidth=0.7)
    pull_axis.axhline(-3.0, color="black", linestyle="--", linewidth=0.7)
    pull_axis.set_ylim(-5.5, 5.5)
    pull_axis.set_ylabel("Pull")
    pull_axis.set_xlabel(r"$M_X^2(e\pi^+)$ (GeV$^2$)")

    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_integrated_region_plot(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    stage: str,
    region: str,
    fit_records: pd.DataFrame,
    output_path: Path,
    histogram_range: tuple[float, float],
    histogram_bins: int,
) -> None:
    """Save integrated Mx2 distributions for FD sectors or the CD sample."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mx2_column = stage_column(stage)

    if region == "FD":
        figure, axes = plt.subplots(
            2,
            3,
            figsize=(12.8, 7.7),
            constrained_layout=True,
        )
        axes_flat = axes.ravel()
        sectors = list(range(1, 7))
    else:
        figure, axis = plt.subplots(
            figsize=(8.2, 5.7),
            constrained_layout=True,
        )
        axes_flat = [axis]
        sectors = [0]
    # endif

    for axis, sector in zip(axes_flat, sectors):
        if region == "FD":
            sample = frame.loc[
                (frame["pion_region"] == "FD")
                & (frame["pion_sector"] == sector),
                mx2_column,
            ].to_numpy()
        else:
            sample = frame.loc[
                frame["pion_region"] == "CD",
                mx2_column,
            ].to_numpy()
        # endif

        counts, edges = np.histogram(
            sample,
            bins=histogram_bins,
            range=histogram_range,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])
        axis.step(centers, counts, where="mid", linewidth=1.0)
        axis.axvline(
            NEUTRON_MASS2_GEV2,
            color="black",
            linestyle="--",
            linewidth=1.0,
        )

        record = fit_records.loc[
            (fit_records["period_key"] == period.key)
            & (fit_records["stage"] == stage)
            & (fit_records["pion_region"] == region)
            & (fit_records["pion_sector"] == sector)
            & (fit_records["cell_type"] == "integrated")
        ]
        if not record.empty:
            row = record.iloc[0]
            axis.set_title(
                (
                    f"Sector {sector}: "
                    if region == "FD"
                    else "CD integrated: "
                )
                + (
                    rf"$\mu_{{Mx^2}}$={row['peak_mean_gev2']:.4f}"
                    if row["success"]
                    else row["status"]
                )
            )
        # endif

        axis.set_xlim(histogram_range)
        axis.margins(x=0.0, y=0.04)
        axis.set_xlabel(r"$M_X^2(e\pi^+)$ (GeV$^2$)")
        axis.set_ylabel("Counts")
    # endfor

    figure.suptitle(
        f"{period.label}: {region} pion missing-neutron distributions "
        f"({stage.replace('_', ' ')})"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_stage_overlay_plot(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    region: str,
    output_path: Path,
    histogram_range: tuple[float, float],
    histogram_bins: int,
) -> None:
    """Overlay uncorrected and electron-corrected Mx2 distributions."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if region == "FD":
        figure, axes = plt.subplots(
            2,
            3,
            figsize=(12.8, 7.7),
            constrained_layout=True,
        )
        axes_flat = axes.ravel()
        sectors = list(range(1, 7))
    else:
        figure, axis = plt.subplots(
            figsize=(8.2, 5.7),
            constrained_layout=True,
        )
        axes_flat = [axis]
        sectors = [0]
    # endif

    for axis, sector in zip(axes_flat, sectors):
        if region == "FD":
            sample = frame.loc[
                (frame["pion_region"] == "FD")
                & (frame["pion_sector"] == sector)
            ]
        else:
            sample = frame.loc[frame["pion_region"] == "CD"]
        # endif

        for column, label in (
            ("mx2_uncorrected_gev2", "No corrections"),
            (
                "mx2_electron_corrected_gev2",
                "Electron correction applied",
            ),
        ):
            counts, edges = np.histogram(
                sample[column].to_numpy(),
                bins=histogram_bins,
                range=histogram_range,
            )
            centers = 0.5 * (edges[:-1] + edges[1:])
            axis.step(
                centers,
                counts,
                where="mid",
                linewidth=1.1,
                label=label,
            )
        # endfor

        axis.axvline(
            NEUTRON_MASS2_GEV2,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label=r"$m_n^2$" if sector in (0, 1) else None,
        )
        axis.set_xlim(histogram_range)
        axis.margins(x=0.0, y=0.04)
        axis.set_title(
            f"Sector {sector}" if region == "FD" else "CD integrated"
        )
        axis.set_xlabel(r"$M_X^2(e\pi^+)$ (GeV$^2$)")
        axis.set_ylabel("Counts")
        if sector in (0, 1):
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(
        f"{period.label}: {region} pion effect of electron correction"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_peak_vs_theta_comparison(
    period: CalibrationPeriod,
    region: str,
    fit_records: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot absolute fitted neutron-peak position before and after e correction."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if region == "FD":
        figure, axes = plt.subplots(
            2,
            3,
            figsize=(13.2, 8.0),
            sharey=True,
            constrained_layout=True,
        )
        axes_flat = axes.ravel()
        sectors = list(range(1, 7))
    else:
        figure, axis = plt.subplots(
            figsize=(8.2, 5.7),
            constrained_layout=True,
        )
        axes_flat = [axis]
        sectors = [0]
    # endif

    for axis, sector in zip(axes_flat, sectors):
        for stage, marker, label in (
            ("uncorrected", "o", "No corrections"),
            (
                "electron_corrected",
                "s",
                "Electron correction applied",
            ),
        ):
            cells = fit_records.loc[
                (fit_records["period_key"] == period.key)
                & (fit_records["stage"] == stage)
                & (fit_records["pion_region"] == region)
                & (fit_records["pion_sector"] == sector)
                & (fit_records["cell_type"] == "theta")
                & (fit_records["success"])
            ].sort_values("mean_theta_deg")

            axis.errorbar(
                cells["mean_theta_deg"],
                cells["peak_mean_gev2"],
                yerr=cells["peak_mean_error_gev2"],
                fmt=marker,
                label=label,
            )
        # endfor

        axis.axhline(
            NEUTRON_MASS2_GEV2,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label=r"$m_n^2$" if sector in (0, 1) else None,
        )
        axis.set_title(
            f"Sector {sector}" if region == "FD" else "CD integrated"
        )
        axis.set_xlabel(r"Mean pion $\theta$ (deg)")
        axis.set_ylabel(r"Fitted $\mu_{Mx^2}$ (GeV$^2$)")
        if sector in (0, 1):
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(
        f"{period.label}: {region} missing-neutron peak versus pion theta"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_event_count_plot(
    period: CalibrationPeriod,
    region: str,
    fit_records: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot event counts in the proposed theta calibration cells."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    cells = fit_records.loc[
        (fit_records["period_key"] == period.key)
        & (fit_records["stage"] == "uncorrected")
        & (fit_records["pion_region"] == region)
        & (fit_records["cell_type"] == "theta")
    ].copy()

    figure, axis = plt.subplots(
        figsize=(10.5, 5.8),
        constrained_layout=True,
    )

    if region == "FD":
        for sector in range(1, 7):
            sector_cells = cells.loc[
                cells["pion_sector"] == sector
            ].sort_values("theta_bin_index")
            axis.plot(
                sector_cells["mean_theta_deg"],
                sector_cells["n_events"],
                "o-",
                label=f"Sector {sector}",
            )
        # endfor
    else:
        cells = cells.sort_values("theta_bin_index")
        axis.plot(
            cells["mean_theta_deg"],
            cells["n_events"],
            "o-",
            label="CD",
        )
    # endif

    axis.set_xlabel(r"Mean pion $\theta$ (deg)")
    axis.set_ylabel("Selected events")
    axis.set_yscale("log")
    axis.set_title(
        f"{period.label}: {region} proposed theta-cell populations"
    )
    axis.legend()
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_run_stability(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    stage: str,
    region: str,
    output_path: Path,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_run: int,
    run_bin_width: int,
) -> list[dict[str, Any]]:
    """Fit Mx2 peak versus run for FD sectors or the integrated CD sample."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mx2_column = stage_column(stage)
    region_frame = frame.loc[frame["pion_region"] == region].copy()

    figure, axis = plt.subplots(
        figsize=(12.2, 6.4),
        constrained_layout=True,
    )
    records: list[dict[str, Any]] = []
    all_runs: list[float] = []

    sectors = range(1, 7) if region == "FD" else [0]
    for sector in sectors:
        if region == "FD":
            sample = region_frame.loc[
                region_frame["pion_sector"] == sector
            ].copy()
        else:
            sample = region_frame.copy()
        # endif

        if sample.empty:
            continue
        # endif

        run_min = int(sample["runnum"].min())
        run_max = int(sample["runnum"].max())
        run_starts = np.arange(
            run_min,
            run_max + 1,
            run_bin_width,
            dtype=int,
        )

        run_centers = []
        peak_values = []
        peak_errors = []

        for run_start in run_starts:
            run_stop = run_start + run_bin_width
            run_sample = sample.loc[
                (sample["runnum"] >= run_start)
                & (sample["runnum"] < run_stop)
            ]
            if run_sample.empty:
                continue
            # endif

            fit_result, _ = fit_mx2_peak(
                run_sample[mx2_column].to_numpy(),
                histogram_range,
                fit_range,
                histogram_bins,
                minimum_events_run,
            )

            record = {
                "period_key": period.key,
                "period_label": period.label,
                "stage": stage,
                "pion_region": region,
                "pion_sector": sector,
                "run_start": int(run_start),
                "run_stop_exclusive": int(run_stop),
                "mean_run": float(run_sample["runnum"].mean()),
                "n_events": int(len(run_sample)),
                "success": bool(fit_result["success"]),
                "status": fit_result["status"],
                "selected_model": fit_result.get("selected_model"),
                "peak_mean_gev2": fit_result.get(
                    "mean_gev2",
                    math.nan,
                ),
                "peak_mean_error_gev2": fit_result.get(
                    "mean_error_gev2",
                    math.nan,
                ),
            }
            records.append(record)

            if fit_result["success"]:
                run_centers.append(record["mean_run"])
                peak_values.append(record["peak_mean_gev2"])
                peak_errors.append(record["peak_mean_error_gev2"])
            # endif
        # endfor

        if run_centers:
            all_runs.extend(run_centers)
            axis.errorbar(
                run_centers,
                peak_values,
                yerr=peak_errors,
                fmt="o",
                markersize=3,
                linewidth=0.8,
                label=(
                    f"Sector {sector}"
                    if region == "FD"
                    else "CD"
                ),
            )
        # endif
    # endfor

    axis.axhline(
        NEUTRON_MASS2_GEV2,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=r"$m_n^2$",
    )
    if all_runs:
        minimum = min(all_runs)
        maximum = max(all_runs)
        padding = max(2.0, 0.012 * max(maximum - minimum, 1.0))
        axis.set_xlim(minimum - padding, maximum + padding)
    # endif
    axis.set_xlabel("Run number")
    axis.set_ylabel(r"Fitted $\mu_{Mx^2}$ (GeV$^2$)")
    axis.set_title(
        f"{period.label}: {region} neutron-peak stability versus run "
        f"({stage.replace('_', ' ')})"
    )
    axis.legend(fontsize=8, ncol=2)
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)
    return records


# -----------------------------------------------------------------------------
# One-period worker
# -----------------------------------------------------------------------------

def process_period(
    period: CalibrationPeriod,
    file_path: Path,
    tree_name_requested: str,
    branch_overrides: dict[str, str],
    electron_model_map: dict[tuple[str, int], dict[str, Any]],
    step_size: str,
    mx2_preselection_min: float,
    mx2_preselection_max: float,
    e_theta_min_deg: float,
    e_theta_max_deg: float,
    pip_theta_min_deg: float,
    pip_theta_max_deg: float,
    fd_theta_bins: int,
    cd_theta_bins: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    minimum_events_run: int,
    run_bin_width: int,
    save_individual_fits: bool,
    skip_run_stability: bool,
    plot_directories: dict[str, Path],
) -> dict[str, Any]:
    """Process one RGC period through uncorrected and electron-corrected stages."""
    tree_name = find_tree(file_path, tree_name_requested)
    branch_map = resolve_branches(
        file_path,
        tree_name,
        branch_overrides,
    )

    print(
        f"[worker:{period.key}] {file_path}:{tree_name}; "
        f"branches={branch_map}",
        flush=True,
    )

    frame = read_period_events(
        file_path,
        tree_name,
        branch_map,
        period,
        step_size,
        mx2_preselection_min,
        mx2_preselection_max,
        e_theta_min_deg,
        e_theta_max_deg,
        pip_theta_min_deg,
        pip_theta_max_deg,
    )
    if frame.empty:
        return {
            "period_key": period.key,
            "success": False,
            "status": "no_selected_events",
            "fit_records": [],
            "run_records": [],
            "branch_map": branch_map,
            "n_events": 0,
        }
    # endif

    frame = apply_electron_correction(
        frame,
        period,
        electron_model_map,
    )
    frame["mx2_electron_corrected_gev2"] = calculate_mx2_epi(
        frame["beam_energy_gev"].to_numpy(),
        frame["e_p_electron_corrected_gev"].to_numpy(),
        frame["e_theta_rad"].to_numpy(),
        frame["e_phi_rad"].to_numpy(),
        frame["pip_p_gev"].to_numpy(),
        frame["pip_theta_rad"].to_numpy(),
        frame["pip_phi_rad"].to_numpy(),
    )

    fit_records: list[dict[str, Any]] = []
    run_records: list[dict[str, Any]] = []

    for stage in ("uncorrected", "electron_corrected"):
        for region in ("FD", "CD"):
            fit_records.extend(
                fit_region_cells(
                    frame,
                    period,
                    stage,
                    region,
                    fd_theta_bins,
                    cd_theta_bins,
                    histogram_range,
                    fit_range,
                    histogram_bins,
                    minimum_events_integrated,
                    minimum_events_cell,
                    plot_directories[
                        f"{stage}_individual_fits"
                    ],
                    save_individual_fits,
                )
            )
        # endfor
    # endfor

    fit_frame = pd.DataFrame(fit_records)

    for stage in ("uncorrected", "electron_corrected"):
        for region in ("FD", "CD"):
            save_integrated_region_plot(
                frame,
                period,
                stage,
                region,
                fit_frame,
                plot_directories[f"{stage}_integrated"]
                / f"{period.key}_{region.lower()}_integrated.png",
                histogram_range,
                histogram_bins,
            )
        # endfor
    # endfor

    for region in ("FD", "CD"):
        save_stage_overlay_plot(
            frame,
            period,
            region,
            plot_directories["comparison_overlays"]
            / f"{period.key}_{region.lower()}_electron_effect.png",
            histogram_range,
            histogram_bins,
        )
        save_peak_vs_theta_comparison(
            period,
            region,
            fit_frame,
            plot_directories["comparison_theta"]
            / f"{period.key}_{region.lower()}_peak_vs_theta.png",
        )
        save_event_count_plot(
            period,
            region,
            fit_frame,
            plot_directories["event_counts"]
            / f"{period.key}_{region.lower()}_theta_cell_counts.png",
        )
    # endfor

    if not skip_run_stability:
        for stage in ("uncorrected", "electron_corrected"):
            for region in ("FD", "CD"):
                run_records.extend(
                    save_run_stability(
                        frame,
                        period,
                        stage,
                        region,
                        plot_directories[f"{stage}_run_stability"]
                        / (
                            f"{period.key}_{region.lower()}_"
                            f"run_stability.png"
                        ),
                        histogram_range,
                        fit_range,
                        histogram_bins,
                        minimum_events_run,
                        run_bin_width,
                    )
                )
            # endfor
        # endfor
    # endif

    return {
        "period_key": period.key,
        "success": True,
        "status": "success",
        "fit_records": fit_records,
        "run_records": run_records,
        "branch_map": branch_map,
        "n_events": int(len(frame)),
        "fd_events": int((frame["pion_region"] == "FD").sum()),
        "cd_events": int((frame["pion_region"] == "CD").sum()),
    }


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------

def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose the e pi+ X missing-neutron peak before and after "
            "application of finalized electron momentum corrections."
        )
    )

    parser.add_argument("--su22-file", type=Path, default=DEFAULT_INPUTS["su22"])
    parser.add_argument("--fa22-file", type=Path, default=DEFAULT_INPUTS["fa22"])
    parser.add_argument("--sp23-file", type=Path, default=DEFAULT_INPUTS["sp23"])
    parser.add_argument("--tree-name", default="PhysicsEvents")
    parser.add_argument(
        "--branch",
        action="append",
        type=parse_branch_override,
        default=[],
        metavar="LOGICAL=ROOT_BRANCH",
        help=(
            "Override one branch mapping. Repeat as needed. Logical names: "
            + ", ".join(sorted(BRANCH_ALIASES))
        ),
    )
    parser.add_argument(
        "--electron-correction-json",
        type=Path,
        default=Path(
            "output/electron_diagnostics/json/"
            "electron_correction_models.json"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/piplus_diagnostics"),
    )
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--step-size", default="250 MB")

    parser.add_argument("--e-theta-min", type=float, default=5.0)
    parser.add_argument("--e-theta-max", type=float, default=35.0)
    parser.add_argument("--pip-theta-min", type=float, default=5.0)
    parser.add_argument("--pip-theta-max", type=float, default=140.0)

    parser.add_argument("--fd-theta-bins", type=int, default=8)
    parser.add_argument("--cd-theta-bins", type=int, default=5)

    parser.add_argument("--mx2-preselection-min", type=float, default=-0.5)
    parser.add_argument("--mx2-preselection-max", type=float, default=2.0)
    parser.add_argument("--mx2-hist-min", type=float, default=0.20)
    parser.add_argument("--mx2-hist-max", type=float, default=1.55)
    parser.add_argument("--mx2-fit-min", type=float, default=0.55)
    parser.add_argument("--mx2-fit-max", type=float, default=1.20)
    parser.add_argument("--mx2-bins", type=int, default=90)

    parser.add_argument(
        "--minimum-events-integrated",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "--minimum-events-cell",
        type=int,
        default=1200,
    )
    parser.add_argument(
        "--minimum-events-run",
        type=int,
        default=1000,
    )
    parser.add_argument("--run-bin-width", type=int, default=1)

    parser.add_argument(
        "--skip-individual-fits",
        action="store_true",
    )
    parser.add_argument(
        "--skip-run-stability",
        action="store_true",
    )
    return parser


def validate_arguments(arguments: argparse.Namespace) -> None:
    """Validate numerical command-line arguments."""
    if arguments.workers < 1:
        raise ValueError("--workers must be positive.")
    # endif
    if arguments.fd_theta_bins < 2:
        raise ValueError("--fd-theta-bins must be at least 2.")
    # endif
    if arguments.cd_theta_bins < 2:
        raise ValueError("--cd-theta-bins must be at least 2.")
    # endif
    if not (
        arguments.mx2_preselection_min
        < arguments.mx2_preselection_max
    ):
        raise ValueError("Invalid Mx2 preselection range.")
    # endif
    if not arguments.mx2_hist_min < arguments.mx2_hist_max:
        raise ValueError("Invalid Mx2 histogram range.")
    # endif
    if not arguments.mx2_fit_min < arguments.mx2_fit_max:
        raise ValueError("Invalid Mx2 fit range.")
    # endif
    if not (
        arguments.mx2_hist_min
        <= arguments.mx2_fit_min
        < arguments.mx2_fit_max
        <= arguments.mx2_hist_max
    ):
        raise ValueError("Mx2 fit range must lie inside histogram range.")
    # endif
    if arguments.mx2_bins < 30:
        raise ValueError("--mx2-bins must be at least 30.")
    # endif
    if arguments.run_bin_width < 1:
        raise ValueError("--run-bin-width must be positive.")
    # endif


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> int:
    parser = build_argument_parser()
    arguments = parser.parse_args()

    try:
        validate_arguments(arguments)
    except ValueError as exc:
        parser.error(str(exc))
    # endtry

    input_files = {
        "su22": arguments.su22_file,
        "fa22": arguments.fa22_file,
        "sp23": arguments.sp23_file,
    }
    for source_key, file_path in input_files.items():
        if not file_path.is_file():
            print(
                f"FATAL: {source_key} input does not exist: {file_path}",
                file=sys.stderr,
            )
            return 2
        # endif
    # endfor

    branch_overrides = dict(arguments.branch)
    electron_model_map = load_electron_models(
        arguments.electron_correction_json
    )

    output_directory = arguments.output_dir
    csv_directory = output_directory / "csv"
    json_directory = output_directory / "json"
    plot_root = output_directory / "plots"

    for directory in (
        csv_directory,
        json_directory,
        plot_root,
    ):
        if directory.exists():
            shutil.rmtree(directory)
        # endif
        directory.mkdir(parents=True, exist_ok=True)
    # endfor

    plot_directories = {
        "uncorrected_individual_fits": (
            plot_root / "uncorrected" / "individual_fits"
        ),
        "uncorrected_integrated": (
            plot_root / "uncorrected" / "integrated_mx2"
        ),
        "uncorrected_run_stability": (
            plot_root / "uncorrected" / "run_stability"
        ),
        "electron_corrected_individual_fits": (
            plot_root / "electron_corrected" / "individual_fits"
        ),
        "electron_corrected_integrated": (
            plot_root / "electron_corrected" / "integrated_mx2"
        ),
        "electron_corrected_run_stability": (
            plot_root / "electron_corrected" / "run_stability"
        ),
        "comparison_overlays": (
            plot_root / "comparison" / "integrated_overlays"
        ),
        "comparison_theta": (
            plot_root / "comparison" / "peak_vs_theta"
        ),
        "event_counts": plot_root / "event_counts",
    }
    for directory in plot_directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    # endfor

    histogram_range = (
        arguments.mx2_hist_min,
        arguments.mx2_hist_max,
    )
    fit_range = (
        arguments.mx2_fit_min,
        arguments.mx2_fit_max,
    )

    worker_count = min(arguments.workers, len(CALIBRATION_PERIODS))
    all_fit_records: list[dict[str, Any]] = []
    all_run_records: list[dict[str, Any]] = []
    period_summaries: list[dict[str, Any]] = []

    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        future_map = {}
        for period in CALIBRATION_PERIODS:
            future = executor.submit(
                process_period,
                period,
                input_files[period.source_key],
                arguments.tree_name,
                branch_overrides,
                electron_model_map,
                arguments.step_size,
                arguments.mx2_preselection_min,
                arguments.mx2_preselection_max,
                arguments.e_theta_min,
                arguments.e_theta_max,
                arguments.pip_theta_min,
                arguments.pip_theta_max,
                arguments.fd_theta_bins,
                arguments.cd_theta_bins,
                histogram_range,
                fit_range,
                arguments.mx2_bins,
                arguments.minimum_events_integrated,
                arguments.minimum_events_cell,
                arguments.minimum_events_run,
                arguments.run_bin_width,
                not arguments.skip_individual_fits,
                arguments.skip_run_stability,
                plot_directories,
            )
            future_map[future] = period
        # endfor

        for future in as_completed(future_map):
            period = future_map[future]
            try:
                result = future.result()
            except Exception as exc:
                print(
                    f"FATAL: {period.label} worker failed: "
                    f"{type(exc).__name__}: {exc}",
                    file=sys.stderr,
                )
                raise
            # endtry

            all_fit_records.extend(result["fit_records"])
            all_run_records.extend(result["run_records"])
            period_summaries.append(
                {
                    key: value
                    for key, value in result.items()
                    if key not in {"fit_records", "run_records"}
                }
            )

            if result["success"]:
                print(
                    f"[complete] {period.label}: "
                    f"{result['n_events']:,} events "
                    f"(FD={result['fd_events']:,}, "
                    f"CD={result['cd_events']:,})"
                )
            else:
                print(
                    f"WARNING: {period.label}: {result['status']}",
                    file=sys.stderr,
                )
            # endif
        # endfor
    # endwith

    fit_frame = pd.DataFrame(all_fit_records)
    run_frame = pd.DataFrame(all_run_records)
    summary_frame = pd.DataFrame(period_summaries)

    fit_csv = csv_directory / "mx2_fit_results.csv"
    run_csv = csv_directory / "mx2_run_stability.csv"
    summary_csv = csv_directory / "period_summary.csv"

    fit_frame.to_csv(fit_csv, index=False, float_format="%.12g")
    run_frame.to_csv(run_csv, index=False, float_format="%.12g")
    summary_frame.to_csv(summary_csv, index=False, float_format="%.12g")

    write_json(
        json_directory / "diagnostic_summary.json",
        {
            "metadata": {
                "neutron_mass_gev": NEUTRON_MASS_GEV,
                "neutron_mass2_gev2": NEUTRON_MASS2_GEV2,
                "electron_correction_json": str(
                    arguments.electron_correction_json
                ),
                "fit_models": list(MODEL_DEFINITIONS),
                "selected_model_policy": (
                    "minimum AICc among successful Gaussian, Crystal Ball, "
                    "and common-mean double-Gaussian signal models"
                ),
                "pi_correction_status": (
                    "not derived or applied in this diagnostic version"
                ),
                "fd_theta_bins": arguments.fd_theta_bins,
                "cd_theta_bins": arguments.cd_theta_bins,
                "mx2_histogram_range_gev2": histogram_range,
                "mx2_fit_range_gev2": fit_range,
            },
            "periods": period_summaries,
        },
    )

    print("")
    print("First-pass pi+ momentum-correction diagnostics complete.")
    print(f"Output directory: {output_directory}")
    print(f"Fit records: {len(fit_frame):,}")
    print(f"Run-stability records: {len(run_frame):,}")
    if not fit_frame.empty:
        print(
            "Successful fits: "
            f"{int(fit_frame['success'].sum()):,}/"
            f"{len(fit_frame):,}"
        )
        selected_counts = (
            fit_frame.loc[fit_frame["success"], "selected_model"]
            .value_counts()
        )
        print(
            "Selected fit models: "
            + ", ".join(
                f"{name}={count}"
                for name, count in selected_counts.items()
            )
        )
    # endif
    print("")
    print(
        "This version stops after comparing no correction with electron-only "
        "correction. Inspect these outputs before deriving pi+ coefficients."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
