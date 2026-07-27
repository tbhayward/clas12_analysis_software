#!/usr/bin/env python3
"""
derive_piplus_momentum_corrections_v12.py

Derive and validate Forward Detector pi+ momentum corrections for the RGC
exclusive-pi+ analysis using inclusive e pi+ X data and the missing-neutron
peak

    Mx2(e pi+) = (p_beam + p_target - p_e' - p_pi+)^2.

Only detector == 1 pion candidates enter this analysis. The finalized electron
correction is applied first. For each calibration period and each of the six
Forward Detector sectors, the script:

    1. forms eight equal-population pion-theta cells;
    2. fits every Mx2 spectrum with a Gaussian signal plus quadratic background;
    3. varies the pion momentum magnitude while holding its direction fixed;
    4. solves for the fractional scale correction that places the fitted
       neutron peak at m_n^2;
    5. fits and selects the smooth theta-dependent correction model;
    6. applies the correction and produces closure diagnostics.

Forward Detector pion calibration, correction models, plots, and output
metadata are intentionally excluded.
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
    "x": ("x", "xB", "x_b", "xb"),
    "pip_p": (
        "p_p",
        "pip_p",
        "pi_p",
        "piplus_p",
        "pion_p",
        "h_p",
    ),
    "pip_theta": (
        "p_theta",
        "pip_theta",
        "pi_theta",
        "piplus_theta",
        "pion_theta",
        "h_theta",
    ),
    "pip_phi": (
        "p_phi",
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
    """Read one period in chunks and retain finite FD e pi+ events."""
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
        x_b = values["x"][period_mask].astype(float)
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
            & np.isfinite(x_b)
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
            & (pip_detector == 1)
        )

        if not np.any(selected):
            continue
        # endif

        runnum = runnum[selected]
        evnum = evnum[selected]
        e_p = e_p[selected]
        e_theta = e_theta[selected]
        e_phi = e_phi[selected]
        x_b = x_b[selected]
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

        fd_sector = fd_sector_from_phi_rad(pip_phi)
        fd_local_phi = fd_local_phi_deg(pip_phi, fd_sector)

        frame = pd.DataFrame(
            {
                "runnum": runnum[mx2_mask],
                "evnum": evnum[mx2_mask],
                "beam_energy_gev": beam_energy[mx2_mask],
                "e_p_gev": e_p[mx2_mask],
                "e_theta_rad": e_theta[mx2_mask],
                "e_theta_deg": e_theta_deg[mx2_mask],
                "e_phi_rad": e_phi[mx2_mask],
                "x_b": x_b[mx2_mask],
                "pip_p_gev": pip_p[mx2_mask],
                "pip_theta_rad": pip_theta[mx2_mask],
                "pip_theta_deg": pip_theta_deg[mx2_mask],
                "pip_phi_rad": pip_phi[mx2_mask],
                "pip_detector": pip_detector[mx2_mask],
                "pion_region": "FD",
                "pion_sector": fd_sector[mx2_mask],
                "pion_local_phi_deg": fd_local_phi[mx2_mask],
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
    sigma_gev2 = float(abs(parameters[2]))
    sigma_error_gev2 = float(errors[2])

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
        "sigma_gev2": sigma_gev2,
        "sigma_error_gev2": sigma_error_gev2,
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
    """Fit one Mx2 spectrum with Gaussian signal plus quadratic background."""
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
                "selected_model": "gaussian",
                "candidate_results": {},
            },
            empty_diagnostics,
        )

    selected = fit_one_model(
        "gaussian",
        centers,
        counts,
        fit_mask,
        fit_range,
    )
    if not selected.get("success", False):
        return (
            {
                "success": False,
                "status": selected.get("status", "fit_failed"),
                "n_events": int(values.size),
                "selected_model": "gaussian",
                "candidate_results": {"gaussian": selected},
            },
            empty_diagnostics,
        )

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
        "fit_signal": gaussian_signal(fit_x, *parameters[:3]),
        "fit_background": quadratic_background(fit_x, *parameters[3:6]),
        "pulls": selected["pulls"],
    }
    summary = {
        "success": True,
        "status": selected["status"],
        "review_reasons": selected["review_reasons"],
        "n_events": int(values.size),
        "selected_model": "gaussian",
        "mean_gev2": selected["mean_gev2"],
        "mean_error_gev2": selected["mean_error_gev2"],
        "sigma_gev2": selected["sigma_gev2"],
        "sigma_error_gev2": selected["sigma_error_gev2"],
        "chi2": selected["chi2"],
        "ndf": selected["ndf"],
        "chi2_ndf": selected["chi2_ndf"],
        "aicc": selected["aicc"],
        "peak_significance_proxy": selected["peak_significance_proxy"],
        "candidate_results": {
            "gaussian": {
                key: value
                for key, value in selected.items()
                if key not in {
                    "fit_x",
                    "fit_y",
                    "fit_error",
                    "fit_model",
                    "pulls",
                }
            }
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
        "peak_sigma_gev2": fit_result.get("sigma_gev2", math.nan),
        "peak_sigma_error_gev2": fit_result.get(
            "sigma_error_gev2",
            math.nan,
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
    mapping = {
        "uncorrected": "mx2_uncorrected_gev2",
        "electron_corrected": "mx2_electron_corrected_gev2",
        "electron_pion_corrected": "mx2_electron_pion_corrected_gev2",
    }
    if stage not in mapping:
        raise ValueError(stage)
    return mapping[stage]


def fit_region_cells(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    stage: str,
    fd_theta_bins: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    individual_fit_directory: Path,
    save_individual_fits: bool,
) -> list[dict[str, Any]]:
    """Fit integrated and theta-binned cells for the six FD sectors."""
    mx2_column = stage_column(stage)
    region = "FD"
    region_frame = frame
    records: list[dict[str, Any]] = []
    theta_bin_count = fd_theta_bins

    for sector in range(1, 7):
        sample = region_frame.loc[
            region_frame["pion_sector"] == sector
        ].copy()

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
                        f"{period.label}, {region}, sector {sector}"
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
# Pi+ correction extraction and smooth models
# -----------------------------------------------------------------------------

def scaled_cell_mx2(
    cell: pd.DataFrame,
    fractional_scale: float,
) -> np.ndarray:
    """Recompute Mx2 after scaling only the pion momentum magnitude."""
    return calculate_mx2_epi(
        cell["beam_energy_gev"].to_numpy(),
        cell["e_p_electron_corrected_gev"].to_numpy(),
        cell["e_theta_rad"].to_numpy(),
        cell["e_phi_rad"].to_numpy(),
        cell["pip_p_gev"].to_numpy() * (1.0 + fractional_scale),
        cell["pip_theta_rad"].to_numpy(),
        cell["pip_phi_rad"].to_numpy(),
    )


def solve_cell_pion_scale(
    cell: pd.DataFrame,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
    derivative_step: float,
    maximum_abs_correction: float,
    maximum_iterations: int = 3,
) -> dict[str, Any]:
    """
    Numerically solve for the fractional pion momentum correction in one cell.

    A Newton update is formed from finite-difference peak fits at scale +/- step.
    The peak is refitted at each iterate.  The solution is clipped to the
    configured physical range.
    """
    scale = 0.0
    history: list[dict[str, Any]] = []

    for iteration in range(maximum_iterations):
        central_values = scaled_cell_mx2(cell, scale)
        central_fit, _ = fit_mx2_peak(
            central_values,
            histogram_range,
            fit_range,
            histogram_bins,
            minimum_events,
        )
        if not central_fit["success"]:
            return {
                "success": False,
                "status": f"central_fit_failed:{central_fit['status']}",
                "fractional_correction": math.nan,
                "fractional_correction_error": math.nan,
                "history": history,
            }

        plus_fit, _ = fit_mx2_peak(
            scaled_cell_mx2(cell, scale + derivative_step),
            histogram_range,
            fit_range,
            histogram_bins,
            minimum_events,
        )
        minus_fit, _ = fit_mx2_peak(
            scaled_cell_mx2(cell, scale - derivative_step),
            histogram_range,
            fit_range,
            histogram_bins,
            minimum_events,
        )
        if not plus_fit["success"] or not minus_fit["success"]:
            return {
                "success": False,
                "status": "derivative_fit_failed",
                "fractional_correction": math.nan,
                "fractional_correction_error": math.nan,
                "history": history,
            }

        derivative = (
            plus_fit["mean_gev2"] - minus_fit["mean_gev2"]
        ) / (2.0 * derivative_step)
        if not np.isfinite(derivative) or abs(derivative) < 1.0e-5:
            return {
                "success": False,
                "status": "invalid_mx2_scale_derivative",
                "fractional_correction": math.nan,
                "fractional_correction_error": math.nan,
                "history": history,
            }

        residual = central_fit["mean_gev2"] - NEUTRON_MASS2_GEV2
        update = -residual / derivative
        next_scale = float(
            np.clip(
                scale + update,
                -maximum_abs_correction,
                maximum_abs_correction,
            )
        )
        history.append(
            {
                "iteration": iteration,
                "scale": scale,
                "peak_mean_gev2": central_fit["mean_gev2"],
                "peak_mean_error_gev2": central_fit["mean_error_gev2"],
                "residual_gev2": residual,
                "dmean_dscale_gev2": derivative,
                "newton_update": update,
                "next_scale": next_scale,
            }
        )
        scale = next_scale
        if abs(residual) < 5.0e-4 or abs(update) < 1.0e-5:
            break

    final_fit, _ = fit_mx2_peak(
        scaled_cell_mx2(cell, scale),
        histogram_range,
        fit_range,
        histogram_bins,
        minimum_events,
    )
    if not final_fit["success"]:
        return {
            "success": False,
            "status": f"final_fit_failed:{final_fit['status']}",
            "fractional_correction": math.nan,
            "fractional_correction_error": math.nan,
            "history": history,
        }

    derivative = history[-1]["dmean_dscale_gev2"]
    correction_error = (
        final_fit["mean_error_gev2"] / abs(derivative)
        if np.isfinite(derivative) and derivative != 0.0
        else math.nan
    )
    at_limit = abs(scale) >= 0.999 * maximum_abs_correction
    return {
        "success": True,
        "status": "success_at_limit" if at_limit else "success",
        "fractional_correction": scale,
        "fractional_correction_error": correction_error,
        "final_peak_mean_gev2": final_fit["mean_gev2"],
        "final_peak_mean_error_gev2": final_fit["mean_error_gev2"],
        "final_residual_gev2": final_fit["mean_gev2"] - NEUTRON_MASS2_GEV2,
        "dmean_dscale_gev2": derivative,
        "history": history,
    }


def polynomial_design(theta_deg: np.ndarray, order: int) -> np.ndarray:
    """Design matrix using powers of centered theta for numerical stability."""
    theta_deg = np.asarray(theta_deg, dtype=float)
    return np.column_stack([theta_deg**power for power in range(order + 1)])


def fit_one_correction_polynomial(
    theta_centered_deg: np.ndarray,
    correction: np.ndarray,
    correction_error: np.ndarray,
    order: int,
) -> dict[str, Any]:
    """Weighted least-squares correction polynomial."""
    design = polynomial_design(theta_centered_deg, order)
    sigma = np.asarray(correction_error, dtype=float)
    fallback = np.nanmedian(sigma[np.isfinite(sigma) & (sigma > 0.0)])
    if not np.isfinite(fallback):
        fallback = 0.002
    sigma = np.where(np.isfinite(sigma) & (sigma > 0.0), sigma, fallback)
    weighted_design = design / sigma[:, None]
    weighted_y = correction / sigma
    coefficients, _, _, _ = np.linalg.lstsq(
        weighted_design,
        weighted_y,
        rcond=None,
    )
    prediction = design @ coefficients
    residual = correction - prediction
    chi2 = float(np.sum((residual / sigma) ** 2))
    n_points = len(correction)
    parameter_count = order + 1
    ndf = n_points - parameter_count
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan
    if n_points > parameter_count + 1:
        aicc = (
            chi2
            + 2.0 * parameter_count
            + 2.0 * parameter_count * (parameter_count + 1)
            / (n_points - parameter_count - 1)
        )
    else:
        aicc = math.inf
    covariance = np.linalg.pinv(weighted_design.T @ weighted_design)
    return {
        "order": order,
        "coefficients": coefficients.tolist(),
        "covariance": covariance.tolist(),
        "chi2": chi2,
        "ndf": int(ndf),
        "chi2_ndf": chi2_ndf,
        "aicc": aicc,
        "adequate": bool(np.isfinite(chi2_ndf) and chi2_ndf <= 2.0),
    }


def select_correction_polynomial(
    theta_deg: np.ndarray,
    correction: np.ndarray,
    correction_error: np.ndarray,
    aicc_promotion_threshold: float,
) -> dict[str, Any]:
    """Select the lowest adequate polynomial, requiring material AICc gain."""
    theta_deg = np.asarray(theta_deg, dtype=float)
    theta_reference = float(np.mean(theta_deg))
    centered = theta_deg - theta_reference
    candidates = [
        fit_one_correction_polynomial(
            centered,
            correction,
            correction_error,
            order,
        )
        for order in (1, 2, 3)
    ]

    selected = candidates[0]
    for candidate in candidates[1:]:
        improvement = selected["aicc"] - candidate["aicc"]
        if (
            (not selected["adequate"] and candidate["adequate"])
            or improvement >= aicc_promotion_threshold
        ):
            selected = candidate
        if selected["adequate"]:
            # Once adequate, higher orders still require the stated AICc gain.
            continue

    if not any(candidate["adequate"] for candidate in candidates):
        selected = min(candidates, key=lambda item: item["aicc"])

    return {
        "theta_reference_deg": theta_reference,
        "selected_order": selected["order"],
        "coefficients": selected["coefficients"],
        "covariance": selected["covariance"],
        "chi2": selected["chi2"],
        "ndf": selected["ndf"],
        "chi2_ndf": selected["chi2_ndf"],
        "aicc": selected["aicc"],
        "adequate": selected["adequate"],
        "candidate_models": candidates,
    }


def evaluate_pion_model(
    theta_deg: np.ndarray,
    model: dict[str, Any],
) -> np.ndarray:
    """Evaluate a selected correction polynomial with boundary clamping."""
    theta = np.asarray(theta_deg, dtype=float)
    theta_clamped = np.clip(
        theta,
        float(model["theta_valid_min_deg"]),
        float(model["theta_valid_max_deg"]),
    )
    centered = theta_clamped - float(model["theta_reference_deg"])
    coefficients = np.asarray(model["coefficients"], dtype=float)
    result = np.zeros_like(centered, dtype=float)
    for power, coefficient in enumerate(coefficients):
        result += coefficient * centered**power
    return result


def derive_fd_pion_models(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    fd_theta_bins: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_cell: int,
    derivative_step: float,
    maximum_abs_correction: float,
    aicc_promotion_threshold: float,
) -> tuple[list[dict[str, Any]], dict[int, dict[str, Any]]]:
    """Extract cell corrections and smooth models for all six FD sectors."""
    records: list[dict[str, Any]] = []
    models: dict[int, dict[str, Any]] = {}
    fd = frame.loc[frame["pion_region"] == "FD"].copy()

    for sector in range(1, 7):
        sample = fd.loc[fd["pion_sector"] == sector].copy()
        edges = equal_population_edges(
            sample["pip_theta_deg"].to_numpy(),
            fd_theta_bins,
        )
        sector_records: list[dict[str, Any]] = []

        for theta_index in range(fd_theta_bins):
            low = float(edges[theta_index])
            high = float(edges[theta_index + 1])
            if theta_index == fd_theta_bins - 1:
                mask = (
                    (sample["pip_theta_deg"] >= low)
                    & (sample["pip_theta_deg"] <= high)
                )
            else:
                mask = (
                    (sample["pip_theta_deg"] >= low)
                    & (sample["pip_theta_deg"] < high)
                )
            cell = sample.loc[mask].copy()
            solution = solve_cell_pion_scale(
                cell,
                histogram_range,
                fit_range,
                histogram_bins,
                minimum_events_cell,
                derivative_step,
                maximum_abs_correction,
            )
            record = {
                "period_key": period.key,
                "period_label": period.label,
                "pion_region": "FD",
                "pion_sector": sector,
                "theta_bin_index": theta_index,
                "theta_low_deg": low,
                "theta_high_deg": high,
                "mean_theta_deg": float(cell["pip_theta_deg"].mean()),
                "mean_pion_momentum_gev": float(cell["pip_p_gev"].mean()),
                "n_events": int(len(cell)),
                "success": bool(solution["success"]),
                "status": solution["status"],
                "fractional_correction": solution.get(
                    "fractional_correction", math.nan
                ),
                "correction_percent": 100.0 * solution.get(
                    "fractional_correction", math.nan
                ),
                "fractional_correction_error": solution.get(
                    "fractional_correction_error", math.nan
                ),
                "correction_error_percent": 100.0 * solution.get(
                    "fractional_correction_error", math.nan
                ),
                "final_peak_mean_gev2": solution.get(
                    "final_peak_mean_gev2", math.nan
                ),
                "final_peak_mean_error_gev2": solution.get(
                    "final_peak_mean_error_gev2", math.nan
                ),
                "final_residual_gev2": solution.get(
                    "final_residual_gev2", math.nan
                ),
                "dmean_dscale_gev2": solution.get(
                    "dmean_dscale_gev2", math.nan
                ),
                "solver_history_json": json.dumps(
                    sanitize_json_value(solution.get("history", [])),
                    sort_keys=True,
                ),
            }
            records.append(record)
            sector_records.append(record)

        successful = [
            record for record in sector_records
            if record["success"]
            and np.isfinite(record["fractional_correction"])
        ]
        if len(successful) < 4:
            models[sector] = {
                "success": False,
                "status": "insufficient_successful_cells",
                "period_key": period.key,
                "pion_sector": sector,
            }
            continue

        theta = np.asarray(
            [record["mean_theta_deg"] for record in successful],
            dtype=float,
        )
        correction = np.asarray(
            [record["fractional_correction"] for record in successful],
            dtype=float,
        )
        correction_error = np.asarray(
            [
                record["fractional_correction_error"]
                for record in successful
            ],
            dtype=float,
        )
        model = select_correction_polynomial(
            theta,
            correction,
            correction_error,
            aicc_promotion_threshold,
        )
        model.update(
            {
                "success": True,
                "status": "success",
                "period_key": period.key,
                "period_label": period.label,
                "pion_sector": sector,
                "theta_valid_min_deg": float(np.min(theta)),
                "theta_valid_max_deg": float(np.max(theta)),
                "n_calibration_cells": len(successful),
            }
        )
        models[sector] = model

    return records, models


def apply_fd_pion_correction(
    frame: pd.DataFrame,
    models: dict[int, dict[str, Any]],
) -> pd.DataFrame:
    """Apply smooth FD correction models event by event."""
    result = frame.copy()
    correction = np.zeros(len(result), dtype=float)
    fd_mask = result["pion_region"].to_numpy() == "FD"
    sectors = result["pion_sector"].to_numpy()

    for sector in range(1, 7):
        mask = fd_mask & (sectors == sector)
        model = models.get(sector, {})
        if not np.any(mask):
            continue
        if not model.get("success", False):
            raise RuntimeError(
                f"No valid pion correction model for FD sector {sector}"
            )
        correction[mask] = evaluate_pion_model(
            result.loc[mask, "pip_theta_deg"].to_numpy(),
            model,
        )

    result["pion_fractional_correction"] = correction
    result["pip_p_electron_pion_corrected_gev"] = (
        result["pip_p_gev"].to_numpy() * (1.0 + correction)
    )
    result["mx2_electron_pion_corrected_gev2"] = calculate_mx2_epi(
        result["beam_energy_gev"].to_numpy(),
        result["e_p_electron_corrected_gev"].to_numpy(),
        result["e_theta_rad"].to_numpy(),
        result["e_phi_rad"].to_numpy(),
        result["pip_p_electron_pion_corrected_gev"].to_numpy(),
        result["pip_theta_rad"].to_numpy(),
        result["pip_phi_rad"].to_numpy(),
    )
    return result


def save_period_correction_model_plot(
    period: CalibrationPeriod,
    correction_records: pd.DataFrame,
    models: dict[int, dict[str, Any]],
    output_path: Path,
) -> None:
    """Plot all six FD pion correction models on one standardized canvas."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    period_cells = correction_records.loc[
        (correction_records["period_key"] == period.key)
        & correction_records["success"]
    ].copy()

    # Establish one common vertical range for every sector in the period.
    plotted_percent_values: list[float] = []
    if not period_cells.empty:
        plotted_percent_values.extend(
            period_cells["correction_percent"]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .tolist()
        )
        error = (
            period_cells["correction_error_percent"]
            .replace([np.inf, -np.inf], np.nan)
            .fillna(0.0)
            .to_numpy()
        )
        value = period_cells["correction_percent"].to_numpy()
        finite = np.isfinite(value) & np.isfinite(error)
        plotted_percent_values.extend((value[finite] - error[finite]).tolist())
        plotted_percent_values.extend((value[finite] + error[finite]).tolist())

    for sector in range(1, 7):
        model = models.get(sector, {})
        if not model.get("success", False):
            continue
        theta = np.linspace(
            float(model["theta_valid_min_deg"]),
            float(model["theta_valid_max_deg"]),
            300,
        )
        plotted_percent_values.extend(
            (100.0 * evaluate_pion_model(theta, model)).tolist()
        )

    finite_values = np.asarray(plotted_percent_values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size:
        y_min = float(np.min(finite_values))
        y_max = float(np.max(finite_values))
        span = max(y_max - y_min, 0.25)
        padding = 0.12 * span
        common_y_limits = (y_min - padding, y_max + padding)
    else:
        common_y_limits = (-2.0, 2.0)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharey=True,
        constrained_layout=True,
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        cells = period_cells.loc[
            period_cells["pion_sector"] == sector
        ].sort_values("mean_theta_deg")
        model = models.get(sector, {})

        axis.errorbar(
            cells["mean_theta_deg"],
            cells["correction_percent"],
            yerr=cells["correction_error_percent"],
            fmt="o",
            label="Numerical cell correction",
        )

        if model.get("success", False):
            theta = np.linspace(
                float(model["theta_valid_min_deg"]),
                float(model["theta_valid_max_deg"]),
                300,
            )
            axis.plot(
                theta,
                100.0 * evaluate_pion_model(theta, model),
                linewidth=1.5,
                label=f"Selected order {model['selected_order']}",
            )
            axis.axvline(
                float(model["theta_valid_min_deg"]),
                linestyle=":",
                linewidth=0.9,
            )
            axis.axvline(
                float(model["theta_valid_max_deg"]),
                linestyle=":",
                linewidth=0.9,
            )
            axis.set_title(
                f"Sector {sector}: order {model['selected_order']}, "
                f"$\\chi^2$/ndf={model['chi2_ndf']:.2f}"
            )
        else:
            axis.set_title(f"Sector {sector}: no valid model")

        axis.axhline(
            0.0,
            color="black",
            linestyle="--",
            linewidth=0.9,
        )
        axis.set_ylim(*common_y_limits)
        axis.set_xlabel(r"Mean pion $\theta$ (deg)")
        axis.set_ylabel(r"$\Delta p_{\pi^+}/p_{\pi^+}$ (%)")
        if sector == 1:
            axis.legend(fontsize=8)

    figure.suptitle(f"{period.label}: FD pion momentum-correction models")
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


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
    """Save integrated Mx2 distributions for all six FD sectors."""
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
        sample = frame.loc[
            frame["pion_sector"] == sector,
            mx2_column,
        ].to_numpy()

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
                f"Sector {sector}: "
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




def fit_mx2_local_summary_peak(
    values: np.ndarray,
    histogram_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """
    Fit the integrated missing-neutron peak in a local signal window.

    The ordinary calibration fits span the broader analysis fit range because
    they must remain uniform across many smaller cells. For the final
    all-sector resolution summary, the strongly rising continuum above about
    1.1 GeV^2 can bias a broad quadratic-background fit. This diagnostic fit
    therefore uses a local Gaussian-plus-linear-background model over
    0.68--1.08 GeV^2. It is used only for the final summary plots and never
    enters correction extraction or model selection.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]

    counts, edges = np.histogram(
        values,
        bins=histogram_bins,
        range=histogram_range,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    local_low = max(histogram_range[0], 0.68)
    local_high = min(histogram_range[1], 1.08)
    fit_mask = (centers >= local_low) & (centers <= local_high)
    fit_x = centers[fit_mask]
    fit_y = counts[fit_mask]
    fit_error = np.sqrt(np.maximum(fit_y, 1.0))

    empty = {
        "counts": counts,
        "edges": edges,
        "centers": centers,
        "fit_x": np.asarray([], dtype=float),
        "fit_model": np.asarray([], dtype=float),
        "fit_signal": np.asarray([], dtype=float),
        "fit_background": np.asarray([], dtype=float),
    }

    if values.size < minimum_events or fit_x.size < 8:
        return {
            "success": False,
            "status": "insufficient_events",
        }, empty

    def model(
        x: np.ndarray,
        amplitude: float,
        mean: float,
        sigma: float,
        c0: float,
        c1: float,
    ) -> np.ndarray:
        return (
            gaussian_signal(x, amplitude, mean, sigma)
            + c0
            + c1 * (x - NEUTRON_MASS2_GEV2)
        )

    sideband = (
        (fit_x < NEUTRON_MASS2_GEV2 - 0.12)
        | (fit_x > NEUTRON_MASS2_GEV2 + 0.12)
    )
    if np.count_nonzero(sideband) >= 2:
        background_guess = np.polyfit(
            fit_x[sideband] - NEUTRON_MASS2_GEV2,
            fit_y[sideband],
            1,
        )
        c1_guess = float(background_guess[0])
        c0_guess = float(background_guess[1])
    else:
        c0_guess = float(np.percentile(fit_y, 20.0))
        c1_guess = 0.0

    peak_window = (
        (fit_x >= NEUTRON_MASS2_GEV2 - 0.12)
        & (fit_x <= NEUTRON_MASS2_GEV2 + 0.12)
    )
    peak_index = np.argmax(
        np.where(peak_window, fit_y, -np.inf)
    )
    mean_guess = float(fit_x[peak_index])
    amplitude_guess = max(
        float(fit_y[peak_index] - c0_guess),
        1.0,
    )

    lower = [
        0.0,
        NEUTRON_MASS2_GEV2 - 0.12,
        0.025,
        -np.inf,
        -np.inf,
    ]
    upper = [
        np.inf,
        NEUTRON_MASS2_GEV2 + 0.12,
        0.18,
        np.inf,
        np.inf,
    ]

    try:
        parameters, covariance = curve_fit(
            model,
            fit_x,
            fit_y,
            p0=[
                amplitude_guess,
                mean_guess,
                0.075,
                c0_guess,
                c1_guess,
            ],
            bounds=(lower, upper),
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
            "success": False,
            "status": f"fit_failed:{type(exc).__name__}",
        }, empty

    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    fit_model = model(fit_x, *parameters)
    fit_signal = gaussian_signal(fit_x, *parameters[:3])
    fit_background = (
        parameters[3]
        + parameters[4] * (fit_x - NEUTRON_MASS2_GEV2)
    )
    chi2 = float(
        np.sum(((fit_y - fit_model) / fit_error) ** 2)
    )
    ndf = int(fit_x.size - len(parameters))

    return {
        "success": True,
        "status": "success",
        "mean_gev2": float(parameters[1]),
        "mean_error_gev2": float(errors[1]),
        "sigma_gev2": float(abs(parameters[2])),
        "sigma_error_gev2": float(errors[2]),
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2 / ndf if ndf > 0 else math.nan,
        "fit_range_low_gev2": local_low,
        "fit_range_high_gev2": local_high,
    }, {
        "counts": counts,
        "edges": edges,
        "centers": centers,
        "fit_x": fit_x,
        "fit_model": fit_model,
        "fit_signal": fit_signal,
        "fit_background": fit_background,
    }


def build_integrated_missing_neutron_summary_payload(
    frame: pd.DataFrame,
    histogram_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
) -> dict[str, Any]:
    """Build the compact payload used by individual and combined summaries."""
    sample = frame.loc[frame["pion_region"] == "FD"].copy()
    stage_definitions = (
        (
            "uncorrected",
            "No corrections",
            "mx2_uncorrected_gev2",
            "tab:blue",
        ),
        (
            "electron_corrected",
            "Electron corrected",
            "mx2_electron_corrected_gev2",
            "tab:orange",
        ),
        (
            "electron_pion_corrected",
            r"Electron + $\pi^+$ corrected",
            "mx2_electron_pion_corrected_gev2",
            "tab:green",
        ),
    )

    stages: list[dict[str, Any]] = []
    for key, label, column, color in stage_definitions:
        fit_result, diagnostics = fit_mx2_local_summary_peak(
            sample[column].to_numpy(),
            histogram_range=histogram_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events,
        )
        stages.append({
            "key": key,
            "label": label,
            "color": color,
            "fit_result": fit_result,
            "centers": diagnostics["centers"].tolist(),
            "counts": diagnostics["counts"].tolist(),
            "fit_x": diagnostics["fit_x"].tolist(),
            "fit_model": diagnostics["fit_model"].tolist(),
        })

    before = stages[0]["fit_result"]
    after = stages[-1]["fit_result"]
    if (
        before.get("success", False)
        and after.get("success", False)
        and before.get("sigma_gev2", 0.0) > 0.0
    ):
        improvement = (
            100.0
            * (
                before["sigma_gev2"]
                - after["sigma_gev2"]
            )
            / before["sigma_gev2"]
        )
    else:
        improvement = math.nan

    return {
        "histogram_range": list(histogram_range),
        "stages": stages,
        "resolution_improvement_percent": improvement,
    }


def plot_integrated_missing_neutron_summary_axis(
    axis: plt.Axes,
    period: CalibrationPeriod,
    payload: dict[str, Any],
    show_title: bool = True,
) -> None:
    """Draw one period's integrated summary on a supplied axis."""
    improvement = float(
        payload.get("resolution_improvement_percent", math.nan)
    )

    for stage in payload["stages"]:
        fit_result = stage["fit_result"]
        label = stage["label"]
        if fit_result.get("success", False):
            label += (
                rf": $\mu={fit_result['mean_gev2']:.5f}$ GeV$^2$, "
                rf"$\sigma={fit_result['sigma_gev2']:.5f}$ GeV$^2$"
            )
            if (
                stage["key"] == "electron_pion_corrected"
                and np.isfinite(improvement)
            ):
                label += (
                    rf", improvement $={improvement:.1f}\%$"
                )

        axis.step(
            np.asarray(stage["centers"], dtype=float),
            np.asarray(stage["counts"], dtype=float),
            where="mid",
            linewidth=1.25,
            color=stage["color"],
            label=label,
        )

        fit_x = np.asarray(stage["fit_x"], dtype=float)
        fit_model = np.asarray(stage["fit_model"], dtype=float)
        if fit_x.size:
            axis.plot(
                fit_x,
                fit_model,
                linewidth=1.7,
                color=stage["color"],
            )

    axis.axvline(
        NEUTRON_MASS2_GEV2,
        color="black",
        linestyle=":",
        linewidth=1.1,
        label=r"$m_n^2$",
    )
    axis.set_xlim(*payload["histogram_range"])
    axis.margins(x=0.0, y=0.05)
    axis.set_xlabel(r"$M_X^2(e\pi^+)$ (GeV$^2$)")
    axis.set_ylabel("Counts")
    if show_title:
        axis.set_title(period.label)
    axis.legend(fontsize=7.2)


def save_fully_integrated_missing_neutron_summary(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    output_path: Path,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
) -> dict[str, Any]:
    """
    Save one period's all-sector integrated summary and return its payload.

    The summary fit is deliberately local to the neutron peak so that the
    high-Mx2 continuum does not distort the fitted Gaussian width.
    """
    del fit_range  # The dedicated local summary fit defines its own window.
    output_path.parent.mkdir(parents=True, exist_ok=True)

    payload = build_integrated_missing_neutron_summary_payload(
        frame=frame,
        histogram_range=histogram_range,
        histogram_bins=histogram_bins,
        minimum_events=minimum_events,
    )

    figure, axis = plt.subplots(
        figsize=(9.4, 6.3),
        constrained_layout=True,
    )
    plot_integrated_missing_neutron_summary_axis(
        axis,
        period,
        payload,
        show_title=False,
    )
    axis.set_title(
        f"{period.label}: all-sector integrated missing-neutron peak"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)
    return payload


def save_combined_integrated_missing_neutron_summary(
    payloads: dict[str, dict[str, Any]],
    output_path: Path,
) -> None:
    """Save the requested 2x2 all-period integrated closure canvas."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(16.0, 11.5),
        constrained_layout=True,
    )

    for period, axis in zip(CALIBRATION_PERIODS, axes.ravel()):
        payload = payloads.get(period.key)
        if payload is None:
            axis.text(
                0.5,
                0.5,
                "Summary unavailable",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
            axis.set_title(period.label)
            continue

        plot_integrated_missing_neutron_summary_axis(
            axis,
            period,
            payload,
            show_title=True,
        )

    figure.suptitle(
        "All-sector integrated missing-neutron peak closure"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_three_stage_integrated_plot(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    fit_records: pd.DataFrame,
    output_path: Path,
    histogram_range: tuple[float, float],
    histogram_bins: int,
) -> None:
    """Overlay all three integrated FD Mx2 stages on one six-sector canvas."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharex=True,
        constrained_layout=True,
    )

    stage_definitions = (
        (
            "uncorrected",
            "mx2_uncorrected_gev2",
            "No corrections",
        ),
        (
            "electron_corrected",
            "mx2_electron_corrected_gev2",
            "Electron corrected",
        ),
        (
            "electron_pion_corrected",
            "mx2_electron_pion_corrected_gev2",
            r"Electron + $\pi^+$ corrected",
        ),
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        sample = frame.loc[
            (frame["pion_region"] == "FD")
            & (frame["pion_sector"] == sector)
        ]

        for stage, column, base_label in stage_definitions:
            counts, edges = np.histogram(
                sample[column].to_numpy(),
                bins=histogram_bins,
                range=histogram_range,
            )
            centers = 0.5 * (edges[:-1] + edges[1:])

            record = fit_records.loc[
                (fit_records["period_key"] == period.key)
                & (fit_records["stage"] == stage)
                & (fit_records["pion_region"] == "FD")
                & (fit_records["pion_sector"] == sector)
                & (fit_records["cell_type"] == "integrated")
                & fit_records["success"]
            ]
            if not record.empty:
                mean_value = float(record.iloc[0]["peak_mean_gev2"])
                label = base_label + rf" ($\mu$={mean_value:.4f})"
            else:
                label = base_label

            axis.step(
                centers,
                counts,
                where="mid",
                linewidth=1.15,
                label=label,
            )

        axis.axvline(
            NEUTRON_MASS2_GEV2,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label=r"$m_n^2$",
        )
        axis.set_xlim(*histogram_range)
        axis.margins(x=0.0, y=0.04)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel(r"$M_X^2(e\pi^+)$ (GeV$^2$)")
        axis.set_ylabel("Counts")
        axis.legend(fontsize=7.5)

    figure.suptitle(
        f"{period.label}: FD integrated missing-neutron distributions"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def save_theta_residual_summary(
    period: CalibrationPeriod,
    fit_records: pd.DataFrame,
    output_path: Path,
) -> None:
    """Compare absolute fitted neutron-peak positions for all three stages."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharey=True,
        constrained_layout=True,
    )

    stage_definitions = (
        ("uncorrected", "o", "No corrections"),
        (
            "electron_corrected",
            "s",
            "Electron corrected",
        ),
        (
            "electron_pion_corrected",
            "^",
            r"Electron + $\pi^+$ corrected",
        ),
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        for stage, marker, label in stage_definitions:
            cells = fit_records.loc[
                (fit_records["period_key"] == period.key)
                & (fit_records["stage"] == stage)
                & (fit_records["pion_region"] == "FD")
                & (fit_records["pion_sector"] == sector)
                & (fit_records["cell_type"] == "theta")
                & fit_records["success"]
            ].sort_values("mean_theta_deg")

            axis.errorbar(
                cells["mean_theta_deg"],
                cells["peak_mean_gev2"],
                yerr=cells["peak_mean_error_gev2"],
                fmt=marker,
                label=label,
            )

        axis.axhline(
            NEUTRON_MASS2_GEV2,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label=r"$m_n^2$" if sector == 1 else None,
        )
        axis.set_ylim(0.8, 1.0)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel(r"Mean pion $\theta$ (deg)")
        axis.set_ylabel(r"Fitted $\mu_{M_X^2}$ (GeV$^2$)")
        if sector == 1:
            axis.legend(fontsize=8)

    figure.suptitle(
        f"{period.label}: pion momentum-correction closure versus theta"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)



def save_theta_width_summary(
    period: CalibrationPeriod,
    fit_records: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Compare fitted Gaussian missing-neutron peak widths for all three stages.

    This is a closure diagnostic only and does not affect the pion correction
    extraction or model selection.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharey=True,
        constrained_layout=True,
    )

    stage_definitions = (
        ("uncorrected", "o", "No corrections"),
        ("electron_corrected", "s", "Electron corrected"),
        (
            "electron_pion_corrected",
            "^",
            r"Electron + $\pi^+$ corrected",
        ),
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        for stage, marker, label in stage_definitions:
            cells = fit_records.loc[
                (fit_records["period_key"] == period.key)
                & (fit_records["stage"] == stage)
                & (fit_records["pion_region"] == "FD")
                & (fit_records["pion_sector"] == sector)
                & (fit_records["cell_type"] == "theta")
                & fit_records["success"]
            ].sort_values("mean_theta_deg")

            axis.errorbar(
                cells["mean_theta_deg"],
                cells["peak_sigma_gev2"],
                yerr=cells["peak_sigma_error_gev2"],
                fmt=marker,
                label=label,
            )

        axis.set_title(f"Sector {sector}")
        axis.set_xlabel(r"Mean pion $\theta$ (deg)")
        axis.set_ylabel(
            r"Fitted $\sigma_{M_X^2}$ (GeV$^2$)"
        )
        if sector == 1:
            axis.legend(fontsize=8)

    figure.suptitle(
        f"{period.label}: missing-neutron peak width versus pion theta"
    )
    figure.savefig(
        output_path,
        dpi=180,
        bbox_inches="tight",
        pad_inches=0.06,
    )
    plt.close(figure)


def _diagnostic_equal_width_edges(values: np.ndarray, bin_count: int) -> np.ndarray:
    """Return equal-width edges spanning the finite observed range."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.linspace(0.0, 1.0, bin_count + 1)
    low = float(np.min(finite))
    high = float(np.max(finite))
    if not high > low:
        high = low + 1.0e-6
    return np.linspace(low, high, bin_count + 1)


def save_kinematic_residual_summary(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    variable_column: str,
    variable_label: str,
    binning: str,
    bin_count: int,
    output_path: Path,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_cell: int,
) -> None:
    """
    Compare fitted neutron-peak positions versus a diagnostic variable.

    The diagnostic fits are fully independent of correction extraction and
    model selection. Pion local phi uses three equal-population bins in each
    pion sector. Momentum and electron-theta diagnostics use eight equal-width
    bins spanning the observed minimum and maximum in each pion sector. The
    same event boundaries are used for all three reconstruction stages.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(
        2, 3, figsize=(15.0, 9.0), sharey=True, constrained_layout=True
    )
    stage_definitions = (
        ("uncorrected", "o", "No corrections"),
        ("electron_corrected", "s", "Electron corrected"),
        ("electron_pion_corrected", "^", r"Electron + $\pi^+$ corrected"),
    )

    fd = frame.loc[frame["pion_region"] == "FD"].copy()
    for sector, axis in zip(range(1, 7), axes.ravel()):
        sample = fd.loc[fd["pion_sector"] == sector].copy()
        values = sample[variable_column].to_numpy(dtype=float)
        if binning == "equal_population":
            edges = equal_population_edges(values, bin_count)
        elif binning == "equal_width":
            edges = _diagnostic_equal_width_edges(values, bin_count)
        else:
            raise ValueError(f"Unsupported diagnostic binning: {binning}")

        for stage, marker, label in stage_definitions:
            mx2_column = stage_column(stage)
            x_values, y_values, y_errors = [], [], []
            for index in range(len(edges) - 1):
                low = float(edges[index])
                high = float(edges[index + 1])
                if index == len(edges) - 2:
                    mask = (
                        (sample[variable_column] >= low)
                        & (sample[variable_column] <= high)
                    )
                else:
                    mask = (
                        (sample[variable_column] >= low)
                        & (sample[variable_column] < high)
                    )
                cell = sample.loc[mask]
                fit_result, _ = fit_mx2_peak(
                    cell[mx2_column].to_numpy(),
                    histogram_range,
                    fit_range,
                    histogram_bins,
                    minimum_events_cell,
                )
                if fit_result["success"]:
                    x_values.append(float(cell[variable_column].mean()))
                    y_values.append(fit_result["mean_gev2"])
                    y_errors.append(fit_result["mean_error_gev2"])

            axis.errorbar(
                x_values, y_values, yerr=y_errors, fmt=marker, label=label
            )

        axis.axhline(
            NEUTRON_MASS2_GEV2, color="black", linestyle="--", linewidth=1.0,
            label=r"$m_n^2$" if sector == 1 else None,
        )
        axis.set_ylim(0.8, 1.0)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel(variable_label)
        axis.set_ylabel(r"Fitted $\mu_{M_X^2}$ (GeV$^2$)")
        if sector == 1:
            axis.legend(fontsize=8)

    figure.suptitle(
        f"{period.label}: pion momentum-correction closure versus {variable_label}"
    )
    figure.savefig(
        output_path, dpi=180, bbox_inches="tight", pad_inches=0.06
    )
    plt.close(figure)


def save_kinematic_width_summary(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    variable_column: str,
    variable_label: str,
    binning: str,
    bin_count: int,
    output_path: Path,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_cell: int,
) -> None:
    """Compare fitted missing-neutron widths versus one diagnostic variable."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(
        2, 3, figsize=(15.0, 9.0), sharey=True, constrained_layout=True
    )
    stage_definitions = (
        ("uncorrected", "o", "No corrections"),
        ("electron_corrected", "s", "Electron corrected"),
        ("electron_pion_corrected", "^", r"Electron + $\pi^+$ corrected"),
    )

    fd = frame.loc[frame["pion_region"] == "FD"].copy()
    for sector, axis in zip(range(1, 7), axes.ravel()):
        sample = fd.loc[fd["pion_sector"] == sector].copy()
        values = sample[variable_column].to_numpy(dtype=float)
        if binning == "equal_population":
            edges = equal_population_edges(values, bin_count)
        elif binning == "equal_width":
            edges = _diagnostic_equal_width_edges(values, bin_count)
        else:
            raise ValueError(f"Unsupported diagnostic binning: {binning}")

        for stage, marker, label in stage_definitions:
            mx2_column = stage_column(stage)
            x_values, y_values, y_errors = [], [], []
            for index in range(len(edges) - 1):
                low = float(edges[index])
                high = float(edges[index + 1])
                if index == len(edges) - 2:
                    mask = (
                        (sample[variable_column] >= low)
                        & (sample[variable_column] <= high)
                    )
                else:
                    mask = (
                        (sample[variable_column] >= low)
                        & (sample[variable_column] < high)
                    )
                cell = sample.loc[mask]
                fit_result, _ = fit_mx2_peak(
                    cell[mx2_column].to_numpy(),
                    histogram_range,
                    fit_range,
                    histogram_bins,
                    minimum_events_cell,
                )
                if fit_result["success"]:
                    x_values.append(float(cell[variable_column].mean()))
                    y_values.append(fit_result["sigma_gev2"])
                    y_errors.append(fit_result["sigma_error_gev2"])

            axis.errorbar(
                x_values, y_values, yerr=y_errors, fmt=marker, label=label
            )

        axis.set_title(f"Sector {sector}")
        axis.set_xlabel(variable_label)
        axis.set_ylabel(r"Fitted $\sigma_{M_X^2}$ (GeV$^2$)")
        if sector == 1:
            axis.legend(fontsize=8)

    figure.suptitle(
        f"{period.label}: missing-neutron resolution versus {variable_label}"
    )
    figure.savefig(
        output_path, dpi=180, bbox_inches="tight", pad_inches=0.06
    )
    plt.close(figure)



def save_integrated_xb_peak_and_width_summary(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    output_path: Path,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_cell: int,
    bin_count: int = 8,
) -> None:
    """
    Make a sector-integrated 1x2 missing-neutron closure canvas versus x_B.

    All six FD pion sectors are combined before fitting in each x_B bin.
    The left panel shows the fitted Gaussian mean and the right panel shows
    the fitted Gaussian width for the uncorrected, electron-corrected, and
    electron-plus-pion-corrected stages.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sample = frame.loc[frame["pion_region"] == "FD"].copy()
    finite_x = sample["x_b"].to_numpy(dtype=float)
    finite_x = finite_x[np.isfinite(finite_x)]
    edges = _diagnostic_equal_width_edges(finite_x, bin_count)

    stages = (
        (
            "No corrections",
            "mx2_uncorrected_gev2",
            "o",
        ),
        (
            "Electron corrected",
            "mx2_electron_corrected_gev2",
            "s",
        ),
        (
            r"Electron + $\pi^+$ corrected",
            "mx2_electron_pion_corrected_gev2",
            "^",
        ),
    )

    figure, (mean_axis, width_axis) = plt.subplots(
        1,
        2,
        figsize=(13.0, 5.2),
        constrained_layout=True,
    )

    for label, mx2_column, marker in stages:
        x_values: list[float] = []
        mean_values: list[float] = []
        mean_errors: list[float] = []
        sigma_values: list[float] = []
        sigma_errors: list[float] = []

        for index in range(len(edges) - 1):
            low = float(edges[index])
            high = float(edges[index + 1])
            if index == len(edges) - 2:
                mask = (sample["x_b"] >= low) & (sample["x_b"] <= high)
            else:
                mask = (sample["x_b"] >= low) & (sample["x_b"] < high)

            cell = sample.loc[mask]
            fit_result, _ = fit_mx2_peak(
                cell[mx2_column].to_numpy(),
                histogram_range,
                fit_range,
                histogram_bins,
                minimum_events_cell,
            )

            if fit_result["success"]:
                x_values.append(float(cell["x_b"].mean()))
                mean_values.append(float(fit_result["mean_gev2"]))
                mean_errors.append(float(fit_result["mean_error_gev2"]))
                sigma_values.append(float(fit_result["sigma_gev2"]))
                sigma_errors.append(float(fit_result["sigma_error_gev2"]))

        mean_axis.errorbar(
            x_values,
            mean_values,
            yerr=mean_errors,
            fmt=marker,
            label=label,
        )
        width_axis.errorbar(
            x_values,
            sigma_values,
            yerr=sigma_errors,
            fmt=marker,
            label=label,
        )

    mean_axis.axhline(
        NEUTRON_MASS2_GEV2,
        color="black",
        linestyle="--",
        linewidth=1.0,
        label=r"$m_n^2$",
    )
    mean_axis.set_xlabel(r"Mean $x_B$")
    mean_axis.set_ylabel(r"Fitted $\mu_{M_X^2}$ (GeV$^2$)")
    mean_axis.set_title("Missing-neutron peak position")
    mean_axis.legend(fontsize=9)

    width_axis.set_xlabel(r"Mean $x_B$")
    width_axis.set_ylabel(r"Fitted $\sigma_{M_X^2}$ (GeV$^2$)")
    width_axis.set_title("Missing-neutron peak resolution")
    width_axis.legend(fontsize=9)

    figure.suptitle(
        f"{period.label}: sector-integrated missing-neutron closure versus $x_B$"
    )
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
    """Fit Mx2 peak versus run for the six FD sectors."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mx2_column = stage_column(stage)
    region_frame = frame.loc[frame["pion_region"] == region].copy()

    figure, axis = plt.subplots(
        figsize=(12.2, 6.4),
        constrained_layout=True,
    )
    records: list[dict[str, Any]] = []
    all_runs: list[float] = []

    sectors = range(1, 7)
    for sector in sectors:
        sample = region_frame.loc[
            region_frame["pion_sector"] == sector
        ].copy()

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
            plotted = axis.errorbar(
                run_centers,
                peak_values,
                yerr=peak_errors,
                fmt="o",
                markersize=3,
                linewidth=0.8,
                label=f"Sector {sector}",
            )
            errors_array = np.asarray(peak_errors, dtype=float)
            values_array = np.asarray(peak_values, dtype=float)
            valid = (
                np.isfinite(values_array)
                & np.isfinite(errors_array)
                & (errors_array > 0.0)
            )
            if np.any(valid):
                weights = 1.0 / errors_array[valid]**2
                weighted_mean = float(
                    np.sum(weights * values_array[valid]) / np.sum(weights)
                )
                line_color = plotted.lines[0].get_color()
                axis.axhline(
                    weighted_mean,
                    color=line_color,
                    linestyle=":",
                    linewidth=1.0,
                    alpha=0.9,
                )
                for record in records:
                    if (
                        record["pion_sector"] == sector
                        and record["stage"] == stage
                    ):
                        record["sector_weighted_mean_gev2"] = weighted_mean
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
    axis.set_ylim(0.8, 1.0)
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





def apply_fd_pion_corrections(
    frame: pd.DataFrame,
    fd_models: dict[int, dict[str, Any]],
) -> pd.DataFrame:
    """Apply the derived FD pion correction model in each sector."""
    result = frame.copy()
    correction = np.zeros(len(result), dtype=float)
    sectors = result["pion_sector"].to_numpy()

    for sector in range(1, 7):
        mask = sectors == sector
        if not np.any(mask):
            continue
        model = fd_models.get(sector, {})
        if not model.get("success", False):
            raise RuntimeError(
                f"No valid pion correction model for FD sector {sector}"
            )
        correction[mask] = evaluate_pion_model(
            result.loc[mask, "pip_theta_deg"].to_numpy(),
            model,
        )

    result["pion_fractional_correction"] = correction
    result["pip_p_electron_pion_corrected_gev"] = (
        result["pip_p_gev"].to_numpy() * (1.0 + correction)
    )
    result["mx2_electron_pion_corrected_gev2"] = calculate_mx2_epi(
        result["beam_energy_gev"].to_numpy(),
        result["e_p_electron_corrected_gev"].to_numpy(),
        result["e_theta_rad"].to_numpy(),
        result["e_phi_rad"].to_numpy(),
        result["pip_p_electron_pion_corrected_gev"].to_numpy(),
        result["pip_theta_rad"].to_numpy(),
        result["pip_phi_rad"].to_numpy(),
    )
    return result

def process_period(
    period: CalibrationPeriod,
    file_path: Path,
    requested_tree: str,
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
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    minimum_events_run: int,
    run_bin_width: int,
    derivative_step: float,
    maximum_abs_correction: float,
    aicc_promotion_threshold: float,
    save_individual_fits: bool,
    skip_run_stability: bool,
    plot_directories: dict[str, Path],
) -> dict[str, Any]:
    """Process one calibration period using FD pion events only."""
    tree_name = find_tree(file_path, requested_tree)
    branch_map = resolve_branches(file_path, tree_name, branch_overrides)
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
            "status": "no_selected_fd_events",
            "fit_records": [],
            "run_records": [],
            "correction_records": [],
            "correction_models": {},
            "branch_map": branch_map,
            "n_events": 0,
            "fd_events": 0,
            "integrated_summary_payload": None,
        }

    frame = apply_electron_correction(frame, period, electron_model_map)
    frame["mx2_electron_corrected_gev2"] = calculate_mx2_epi(
        frame["beam_energy_gev"].to_numpy(),
        frame["e_p_electron_corrected_gev"].to_numpy(),
        frame["e_theta_rad"].to_numpy(),
        frame["e_phi_rad"].to_numpy(),
        frame["pip_p_gev"].to_numpy(),
        frame["pip_theta_rad"].to_numpy(),
        frame["pip_phi_rad"].to_numpy(),
    )

    correction_records, fd_models = derive_fd_pion_models(
        frame,
        period,
        fd_theta_bins,
        histogram_range,
        fit_range,
        histogram_bins,
        minimum_events_cell,
        derivative_step,
        maximum_abs_correction,
        aicc_promotion_threshold,
    )

    invalid_fd = [
        sector
        for sector, model in fd_models.items()
        if not model.get("success", False)
    ]
    if invalid_fd:
        raise RuntimeError(
            f"{period.label}: invalid FD correction models in sectors "
            f"{invalid_fd}"
        )

    frame = apply_fd_pion_corrections(frame, fd_models)
    correction_frame = pd.DataFrame(correction_records)

    save_period_correction_model_plot(
        period,
        correction_frame,
        fd_models,
        plot_directories["correction_models"]
        / f"{period.key}_fd_pion_correction_models.png",
    )

    fit_records: list[dict[str, Any]] = []
    run_records: list[dict[str, Any]] = []
    stages = (
        "uncorrected",
        "electron_corrected",
        "electron_pion_corrected",
    )

    for stage in stages:
        fit_records.extend(
            fit_region_cells(
                frame,
                period,
                stage,
                fd_theta_bins,
                histogram_range,
                fit_range,
                histogram_bins,
                minimum_events_integrated,
                minimum_events_cell,
                plot_directories[f"{stage}_individual_fits"],
                save_individual_fits,
            )
        )

    fit_frame = pd.DataFrame(fit_records)

    save_three_stage_integrated_plot(
        frame,
        period,
        fit_frame,
        plot_directories["integrated_mx2"]
        / f"{period.key}_fd_integrated_three_stage.png",
        histogram_range,
        histogram_bins,
    )
    integrated_summary_payload = (
        save_fully_integrated_missing_neutron_summary(
            frame=frame,
            period=period,
            output_path=plot_directories["final_summary"]
            / (
                f"{period.key}_all_sector_integrated_"
                "missing_neutron_summary.png"
            ),
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events_integrated,
        )
    )
    save_theta_residual_summary(
        period,
        fit_frame,
        plot_directories["theta_residuals"]
        / f"{period.key}_fd_theta_residuals.png",
    )
    save_theta_width_summary(
        period,
        fit_frame,
        plot_directories["theta_residuals"]
        / f"{period.key}_fd_theta_widths.png",
    )

    for variable_column, variable_label, binning, bin_count, key, suffix in (
        (
            "pip_p_gev",
            r"Mean pion $p$ (GeV)",
            "equal_width",
            8,
            "pion_momentum_residuals",
            "pion_momentum_residuals",
        ),
        (
            "pion_local_phi_deg",
            r"Mean pion local $\phi$ (deg)",
            "equal_population",
            3,
            "pion_phi_residuals",
            "pion_phi_residuals",
        ),
        (
            "e_p_gev",
            r"Mean electron $p$ (GeV)",
            "equal_width",
            8,
            "electron_momentum_residuals",
            "electron_momentum_residuals",
        ),
        (
            "e_theta_deg",
            r"Mean electron $\theta$ (deg)",
            "equal_width",
            8,
            "electron_theta_residuals",
            "electron_theta_residuals",
        ),
    ):
        save_kinematic_residual_summary(
            frame,
            period,
            variable_column,
            variable_label,
            binning,
            bin_count,
            plot_directories[key] / f"{period.key}_fd_{suffix}.png",
            histogram_range,
            fit_range,
            histogram_bins,
            minimum_events_cell,
        )

    save_integrated_xb_peak_and_width_summary(
        frame,
        period,
        plot_directories["x_residuals"]
        / f"{period.key}_fd_sector_integrated_xb_peak_and_width.png",
        histogram_range,
        fit_range,
        histogram_bins,
        minimum_events_cell,
        8,
    )

    if not skip_run_stability:
        for stage in stages:
            run_records.extend(
                save_run_stability(
                    frame,
                    period,
                    stage,
                    "FD",
                    plot_directories[f"{stage}_run_stability"]
                    / f"{period.key}_fd_run_stability.png",
                    histogram_range,
                    fit_range,
                    histogram_bins,
                    minimum_events_run,
                    run_bin_width,
                )
            )

    model_payload = {
        "FD": {
            str(sector): sanitize_json_value(model)
            for sector, model in fd_models.items()
        }
    }
    fd_events = int(len(frame))
    return {
        "period_key": period.key,
        "success": True,
        "status": "success",
        "fit_records": fit_records,
        "run_records": run_records,
        "correction_records": correction_records,
        "correction_models": model_payload,
        "branch_map": branch_map,
        "n_events": fd_events,
        "fd_events": fd_events,
        "integrated_summary_payload": integrated_summary_payload,
    }


# -----------------------------------------------------------------------------
# Command line
# -----------------------------------------------------------------------------


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Derive FD pi+ momentum corrections from the e pi+ X "
            "missing-neutron peak and perform full closure."
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
    parser.add_argument("--pip-theta-max", type=float, default=45.0)
    parser.add_argument("--fd-theta-bins", type=int, default=8)

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
        "--pion-scale-derivative-step",
        type=float,
        default=0.0025,
        help="Finite-difference step in fractional pion momentum.",
    )
    parser.add_argument(
        "--maximum-abs-pion-correction",
        type=float,
        default=0.10,
        help="Maximum allowed absolute fractional pion correction.",
    )
    parser.add_argument(
        "--aicc-promotion-threshold",
        type=float,
        default=2.0,
        help="Required AICc improvement before promoting polynomial order.",
    )
    parser.add_argument("--skip-individual-fits", action="store_true")
    parser.add_argument("--skip-run-stability", action="store_true")
    return parser


def validate_arguments(arguments: argparse.Namespace) -> None:
    if arguments.workers < 1:
        raise ValueError("--workers must be positive.")
    if arguments.fd_theta_bins < 4:
        raise ValueError("--fd-theta-bins must be at least 4.")
    if not arguments.mx2_preselection_min < arguments.mx2_preselection_max:
        raise ValueError("Invalid Mx2 preselection range.")
    if not arguments.mx2_hist_min < arguments.mx2_hist_max:
        raise ValueError("Invalid Mx2 histogram range.")
    if not arguments.mx2_fit_min < arguments.mx2_fit_max:
        raise ValueError("Invalid Mx2 fit range.")
    if not (
        arguments.mx2_hist_min
        <= arguments.mx2_fit_min
        < arguments.mx2_fit_max
        <= arguments.mx2_hist_max
    ):
        raise ValueError("Mx2 fit range must lie inside histogram range.")
    if arguments.mx2_bins < 30:
        raise ValueError("--mx2-bins must be at least 30.")
    if arguments.run_bin_width < 1:
        raise ValueError("--run-bin-width must be positive.")
    if not 0.0 < arguments.pion_scale_derivative_step < 0.02:
        raise ValueError("Invalid --pion-scale-derivative-step.")
    if not 0.0 < arguments.maximum_abs_pion_correction <= 0.25:
        raise ValueError("Invalid --maximum-abs-pion-correction.")
    if arguments.aicc_promotion_threshold < 0.0:
        raise ValueError("--aicc-promotion-threshold cannot be negative.")


def main() -> int:
    parser = build_argument_parser()
    arguments = parser.parse_args()
    try:
        validate_arguments(arguments)
    except ValueError as exc:
        parser.error(str(exc))

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

    if not arguments.electron_correction_json.is_file():
        print(
            "FATAL: electron correction JSON does not exist: "
            f"{arguments.electron_correction_json}",
            file=sys.stderr,
        )
        return 2

    branch_overrides = dict(arguments.branch)
    electron_model_map = load_electron_models(
        arguments.electron_correction_json
    )

    output_directory = arguments.output_dir
    csv_directory = output_directory / "csv"
    json_directory = output_directory / "json"
    plot_root = output_directory / "plots"
    for directory in (csv_directory, json_directory, plot_root):
        if directory.exists():
            shutil.rmtree(directory)
        directory.mkdir(parents=True, exist_ok=True)

    plot_directories = {
        "uncorrected_individual_fits":
            plot_root / "uncorrected" / "individual_fits",
        "uncorrected_run_stability":
            plot_root / "uncorrected" / "run_stability",
        "electron_corrected_individual_fits":
            plot_root / "electron_corrected" / "individual_fits",
        "electron_corrected_run_stability":
            plot_root / "electron_corrected" / "run_stability",
        "electron_pion_corrected_individual_fits":
            plot_root / "electron_pion_corrected" / "individual_fits",
        "electron_pion_corrected_run_stability":
            plot_root / "electron_pion_corrected" / "run_stability",
        "integrated_mx2":
            plot_root / "closure" / "integrated_mx2",
        "final_summary":
            plot_root / "closure" / "final_summary",
        "theta_residuals":
            plot_root / "closure" / "theta_residuals",
        "pion_momentum_residuals":
            plot_root / "closure" / "pion_momentum_residuals",
        "pion_phi_residuals":
            plot_root / "closure" / "pion_phi_residuals",
        "electron_momentum_residuals":
            plot_root / "closure" / "electron_momentum_residuals",
        "electron_theta_residuals":
            plot_root / "closure" / "electron_theta_residuals",
        "x_residuals":
            plot_root / "closure" / "x_residuals",
        "correction_models":
            plot_root / "correction_models",
    }
    for directory in plot_directories.values():
        directory.mkdir(parents=True, exist_ok=True)

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
    all_correction_records: list[dict[str, Any]] = []
    all_models: dict[str, Any] = {}
    period_summaries: list[dict[str, Any]] = []
    integrated_summary_payloads: dict[str, dict[str, Any]] = {}

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
                histogram_range,
                fit_range,
                arguments.mx2_bins,
                arguments.minimum_events_integrated,
                arguments.minimum_events_cell,
                arguments.minimum_events_run,
                arguments.run_bin_width,
                arguments.pion_scale_derivative_step,
                arguments.maximum_abs_pion_correction,
                arguments.aicc_promotion_threshold,
                not arguments.skip_individual_fits,
                arguments.skip_run_stability,
                plot_directories,
            )
            future_map[future] = period

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

            all_fit_records.extend(result["fit_records"])
            all_run_records.extend(result["run_records"])
            all_correction_records.extend(result["correction_records"])
            all_models[period.key] = result["correction_models"]
            if result.get("integrated_summary_payload") is not None:
                integrated_summary_payloads[period.key] = result[
                    "integrated_summary_payload"
                ]
            period_summaries.append(
                {
                    key: value
                    for key, value in result.items()
                    if key not in {
                        "fit_records",
                        "run_records",
                        "correction_records",
                        "correction_models",
                        "integrated_summary_payload",
                    }
                }
            )
            if result["success"]:
                print(
                    f"[complete] {period.label}: "
                    f"{result['fd_events']:,} FD events"
                )
            else:
                print(
                    f"WARNING: {period.label}: {result['status']}",
                    file=sys.stderr,
                )

    fit_frame = pd.DataFrame(all_fit_records)
    run_frame = pd.DataFrame(all_run_records)
    correction_frame = pd.DataFrame(all_correction_records)
    summary_frame = pd.DataFrame(period_summaries)

    save_combined_integrated_missing_neutron_summary(
        integrated_summary_payloads,
        plot_directories["final_summary"]
        / "all_periods_all_sector_integrated_missing_neutron_summary.png",
    )

    fit_frame.to_csv(
        csv_directory / "mx2_fit_results.csv",
        index=False,
        float_format="%.12g",
    )
    run_frame.to_csv(
        csv_directory / "mx2_run_stability.csv",
        index=False,
        float_format="%.12g",
    )
    correction_frame.to_csv(
        csv_directory / "pion_correction_cells.csv",
        index=False,
        float_format="%.12g",
    )
    summary_frame.to_csv(
        csv_directory / "period_summary.csv",
        index=False,
        float_format="%.12g",
    )

    model_payload = {
        "metadata": {
            "script": "derive_piplus_momentum_corrections_v12.py",
            "neutron_mass_gev": NEUTRON_MASS_GEV,
            "neutron_mass2_gev2": NEUTRON_MASS2_GEV2,
            "fit_model": "Gaussian signal plus quadratic background",
            "detector_scope": "Forward Detector pions only (detector == 1)",
            "electron_correction_json": str(
                arguments.electron_correction_json
            ),
            "fd_theta_bins": arguments.fd_theta_bins,
            "theta_extrapolation": (
                "clamp to minimum/maximum calibration-cell mean theta"
            ),
            "mx2_histogram_range_gev2": histogram_range,
            "mx2_fit_range_gev2": fit_range,
            "pion_scale_derivative_step": (
                arguments.pion_scale_derivative_step
            ),
            "maximum_abs_pion_correction": (
                arguments.maximum_abs_pion_correction
            ),
            "aicc_promotion_threshold": (
                arguments.aicc_promotion_threshold
            ),
        },
        "periods": all_models,
    }
    write_json(
        json_directory / "pion_correction_models.json",
        model_payload,
    )
    write_json(
        json_directory / "analysis_summary.json",
        {
            "metadata": model_payload["metadata"],
            "period_summaries": period_summaries,
        },
    )

    print("")
    print("FD pi+ momentum-correction extraction and closure complete.")
    print(f"Output directory: {output_directory}")
    print(f"Calibration cells: {len(correction_frame):,}")
    print(f"Mx2 fit records: {len(fit_frame):,}")
    print(f"Run-stability records: {len(run_frame):,}")
    print(
        "Correction JSON: "
        f"{json_directory / 'pion_correction_models.json'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
