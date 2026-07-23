#!/usr/bin/env python3
"""
apply_epi_momentum_corrections_v5.py

Rewrite an existing e pi+ X ROOT tree after selecting Forward Detector pion
events (detector == 1) and applying the finalized RGC electron and FD pi+
momentum corrections. The electron and pion directions are held fixed; only
their momentum magnitudes are corrected.

The output tree contains the requested reduced branch set:

    runnum, evnum, helicity, detector, beam_pol, target_pol,
    e_p, e_theta, e_phi, vz_e,
    p_p, p_theta, p_phi, vz_p,
    isrTheta, isrPhi,
    Q2, W, Mx2, x, y, t, tmin, tprime, phi,
    DepA, DepB, DepC, DepV, DepW

Definitions after correction:

    q       = k - k'
    Q2      = -q^2
    W^2     = (P + q)^2
    Mx2     = (P + q - p_pi)^2
    x       = Q2 / (2 M_p nu)
    y       = nu / E_beam
    t       = (q - p_pi)^2

`tmin` is the exact forward two-body value for

    gamma* p -> pi+ n,

using the charged-pion and neutron masses. This is the mesonic definition,
rather than the recoil-hadron/DVCS expression used in the original one-hadron
tree producer.

`phi` is the Trento pion azimuth on [0, 2 pi). `isrTheta` and `isrPhi` are
copied without modification.

Corrected events are required to satisfy 0.1 < x < 0.6.
By default, they must also satisfy Mx2 < 1.5 GeV^2.
Use --no-mx2-skim to retain every successfully reconstructed event.

Momentum-model conventions reproduced from the extraction scripts:

    electron:
        delta p / p = polynomial(theta_e in degrees)
        theta is clamped to the model validity interval;

    FD pi+:
        delta p / p = polynomial(theta_pi - theta_reference)
        theta is clamped to the model validity interval.

For pion detector != 1, no pion calibration exists; p_p is retained unchanged.
The electron correction is still applied.

Example:
    python3 apply_epi_momentum_corrections_v5.py \
        input.root output_corrected.root

Example without the default Mx2 skim:
    python3 apply_epi_momentum_corrections_v5.py \
        input.root output_corrected.root --no-mx2-skim

Dependencies:
    numpy
    awkward
    uproot
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import awkward as ak
import numpy as np
import uproot


# =============================================================================
# Physics constants
# =============================================================================

ELECTRON_MASS_GEV = 0.00051099895
PROTON_MASS_GEV = 0.9382720813
NEUTRON_MASS_GEV = 0.9395654133
PION_MASS_GEV = 0.13957039
RAD2DEG = 180.0 / math.pi
TWO_PI = 2.0 * math.pi

DEFAULT_ELECTRON_JSON = Path(
    "output/electron_diagnostics/json/electron_correction_models.json"
)
DEFAULT_PION_JSON = Path(
    "output/piplus_diagnostics/json/pion_correction_models.json"
)

OUTPUT_BRANCH_DTYPES: dict[str, np.dtype] = {
    "runnum": np.dtype(np.int32),
    "evnum": np.dtype(np.int32),
    "helicity": np.dtype(np.int32),
    "detector": np.dtype(np.int32),
    "beam_pol": np.dtype(np.float64),
    "target_pol": np.dtype(np.float64),
    "e_p": np.dtype(np.float64),
    "e_theta": np.dtype(np.float64),
    "e_phi": np.dtype(np.float64),
    "vz_e": np.dtype(np.float64),
    "p_p": np.dtype(np.float64),
    "p_theta": np.dtype(np.float64),
    "p_phi": np.dtype(np.float64),
    "vz_p": np.dtype(np.float64),
    "isrTheta": np.dtype(np.float64),
    "isrPhi": np.dtype(np.float64),
    "Q2": np.dtype(np.float64),
    "W": np.dtype(np.float64),
    "Mx2": np.dtype(np.float64),
    "x": np.dtype(np.float64),
    "y": np.dtype(np.float64),
    "t": np.dtype(np.float64),
    "tmin": np.dtype(np.float64),
    "tprime": np.dtype(np.float64),
    "phi": np.dtype(np.float64),
    "DepA": np.dtype(np.float64),
    "DepB": np.dtype(np.float64),
    "DepC": np.dtype(np.float64),
    "DepV": np.dtype(np.float64),
    "DepW": np.dtype(np.float64),
}

REQUIRED_INPUT_BRANCHES = tuple(
    name for name in OUTPUT_BRANCH_DTYPES.keys()
    if name != "tprime"
)
COPIED_INPUT_BRANCHES = (
    "runnum",
    "evnum",
    "helicity",
    "detector",
    "beam_pol",
    "target_pol",
    "e_theta",
    "e_phi",
    "vz_e",
    "p_theta",
    "p_phi",
    "vz_p",
    "isrTheta",
    "isrPhi",
)
RECALCULATED_INPUTS = ("e_p", "p_p")


# =============================================================================
# Run periods and calibration records
# =============================================================================

@dataclass(frozen=True)
class CalibrationPeriod:
    key: str
    label: str
    run_min: int
    run_max: int


CALIBRATION_PERIODS = (
    CalibrationPeriod("su22", "Su22", 16042, 16788),
    CalibrationPeriod(
        "fa22_solenoid_minus",
        "Fa22, solenoid -1",
        16843,
        17183,
    ),
    CalibrationPeriod(
        "fa22_solenoid_plus",
        "Fa22, solenoid +1",
        17185,
        17408,
    ),
    CalibrationPeriod("sp23", "Sp23", 17477, 17811),
)


def period_key_from_run(runnum: np.ndarray) -> np.ndarray:
    """Map each event to the correction-period key."""
    run = np.asarray(runnum, dtype=np.int64)
    result = np.full(run.shape, "", dtype="U32")
    for period in CALIBRATION_PERIODS:
        mask = (run >= period.run_min) & (run <= period.run_max)
        result[mask] = period.key
    return result


def beam_energy_from_run(runnum: np.ndarray) -> np.ndarray:
    """
    Assign the nominal RGC beam energy event by event.

        16042-17065: 10.5473 GeV
        17067-17716: 10.5563 GeV
        17717-17811: 10.5593 GeV
    """
    run = np.asarray(runnum, dtype=np.int64)
    energy = np.full(run.shape, np.nan, dtype=np.float64)
    energy[(run >= 16042) & (run <= 17065)] = 10.5473
    energy[(run >= 17067) & (run <= 17716)] = 10.5563
    energy[(run >= 17717) & (run <= 17811)] = 10.5593
    return energy


# =============================================================================
# ROOT and model-file handling
# =============================================================================

def find_tree(file_path: Path, requested_tree: str) -> str:
    """Resolve a ROOT tree name, accepting a cycle suffix."""
    with uproot.open(file_path) as root_file:
        direct_names = {
            str(key).split(";")[0]: str(key)
            for key in root_file.keys()
        }
        if requested_tree in direct_names:
            return direct_names[requested_tree]

        candidates: list[str] = []
        for key in root_file.keys():
            try:
                obj = root_file[key]
                if hasattr(obj, "arrays") and hasattr(obj, "num_entries"):
                    candidates.append(str(key))
            except Exception:
                continue

    if len(candidates) == 1:
        return candidates[0]

    raise KeyError(
        f"Tree '{requested_tree}' was not found in {file_path}. "
        f"Tree-like objects: {candidates}"
    )


def validate_input_branches(file_path: Path, tree_name: str) -> None:
    """Require the branches needed for correction and reduced-tree output."""
    with uproot.open(file_path) as root_file:
        available = set(root_file[tree_name].keys())

    missing = sorted(set(REQUIRED_INPUT_BRANCHES) - available)
    if missing:
        raise KeyError(
            "Input tree is missing required branches: "
            + ", ".join(missing)
        )


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as stream:
        payload = json.load(stream)
    if not isinstance(payload, dict):
        raise TypeError(f"Top level of {path} is not a JSON object.")
    return payload


def load_electron_models(
    path: Path,
) -> dict[tuple[str, int], dict[str, Any]]:
    """
    Load the selected electron models written by
    derive_electron_momentum_corrections_v19.py.
    """
    payload = load_json(path)
    records = payload.get("models")
    if not isinstance(records, list):
        raise KeyError(
            f"{path} does not contain the expected top-level 'models' list."
        )

    models: dict[tuple[str, int], dict[str, Any]] = {}
    for record in records:
        if not isinstance(record, dict) or not record.get("success", False):
            continue
        key = (str(record["period_key"]), int(record["sector"]))
        coefficients = record.get("selected_coefficients_ascending", [])
        if not coefficients:
            continue
        models[key] = record

    expected = {
        (period.key, sector)
        for period in CALIBRATION_PERIODS
        for sector in range(1, 7)
    }
    missing = sorted(expected - set(models))
    if missing:
        raise RuntimeError(
            "Electron model JSON is missing successful period/sector models: "
            f"{missing}"
        )
    return models


def _normalise_pion_period_models(
    period_key: str,
    raw_period: Any,
) -> Iterable[tuple[int, dict[str, Any]]]:
    """
    Yield the successful FD pion models for one calibration period.

    Supported JSON layouts include:

      periods[period][sector] = model
      periods[period]["FD"][sector] = model
      periods[period] = [model, ...]

    Region containers such as "CD" are ignored because the finalized pion
    momentum correction is defined only for FD pions.
    """

    def parse_sector(value: Any) -> int | None:
        try:
            sector = int(value)
        except (TypeError, ValueError):
            return None
        return sector if 1 <= sector <= 6 else None

    if isinstance(raw_period, list):
        for model in raw_period:
            if not isinstance(model, dict):
                continue

            region = str(
                model.get("pion_region", model.get("region", "FD"))
            ).upper()
            if region != "FD":
                continue

            sector = parse_sector(
                model.get("pion_sector", model.get("sector"))
            )
            if sector is not None:
                yield sector, model
        return

    if not isinstance(raw_period, dict):
        raise TypeError(
            f"Unexpected pion-model container for period {period_key}: "
            f"{type(raw_period).__name__}"
        )

    # Newer correction JSONs may explicitly group models by detector region.
    # Select only FD; never attempt to interpret "CD" as a sector number.
    for fd_key in ("FD", "fd", "Forward Detector", "forward_detector"):
        if fd_key in raw_period:
            yield from _normalise_pion_period_models(
                period_key,
                raw_period[fd_key],
            )
            return

    # A dictionary can itself be a single model.
    direct_sector = parse_sector(
        raw_period.get("pion_sector", raw_period.get("sector"))
    )
    if direct_sector is not None and (
        "coefficients" in raw_period or "success" in raw_period
    ):
        region = str(
            raw_period.get("pion_region", raw_period.get("region", "FD"))
        ).upper()
        if region == "FD":
            yield direct_sector, raw_period
        return

    for raw_sector, model in raw_period.items():
        # Explicit non-FD region wrappers are not calibration models.
        if str(raw_sector).upper() in {
            "CD",
            "CENTRAL DETECTOR",
            "CENTRAL_DETECTOR",
            "FT",
            "FORWARD TAGGER",
            "FORWARD_TAGGER",
        }:
            continue

        if not isinstance(model, (dict, list)):
            continue

        if isinstance(model, dict):
            region = str(
                model.get("pion_region", model.get("region", "FD"))
            ).upper()
            if region != "FD":
                continue

            sector = parse_sector(
                model.get("pion_sector", model.get("sector", raw_sector))
            )
            if sector is not None and (
                "coefficients" in model or "success" in model
            ):
                yield sector, model
                continue

        # Tolerate an additional harmless nesting level.
        yield from _normalise_pion_period_models(period_key, model)


def load_pion_models(
    path: Path,
) -> dict[tuple[str, int], dict[str, Any]]:
    """
    Load the selected FD pion models written by
    derive_piplus_momentum_corrections_v9.py.
    """
    payload = load_json(path)
    periods = payload.get("periods")
    if not isinstance(periods, dict):
        raise KeyError(
            f"{path} does not contain the expected top-level 'periods' object."
        )

    models: dict[tuple[str, int], dict[str, Any]] = {}
    for period_key, raw_period in periods.items():
        for sector, model in _normalise_pion_period_models(
            str(period_key), raw_period
        ):
            if not model.get("success", False):
                continue
            coefficients = model.get("coefficients", [])
            if not coefficients:
                continue
            models[(str(period_key), sector)] = model

    expected = {
        (period.key, sector)
        for period in CALIBRATION_PERIODS
        for sector in range(1, 7)
    }
    missing = sorted(expected - set(models))
    if missing:
        raise RuntimeError(
            "Pion model JSON is missing successful FD period/sector models: "
            f"{missing}"
        )
    return models


# =============================================================================
# Momentum corrections
# =============================================================================

def wrap_phi_rad(phi_rad: np.ndarray) -> np.ndarray:
    return np.mod(np.asarray(phi_rad, dtype=np.float64), TWO_PI)


def fd_sector_from_phi_rad(phi_rad: np.ndarray) -> np.ndarray:
    """Assign the six standard CLAS12 FD sectors."""
    phi_deg = np.mod(
        np.asarray(phi_rad, dtype=np.float64) * RAD2DEG,
        360.0,
    )
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)
    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def evaluate_electron_model(
    theta_deg: np.ndarray,
    model: dict[str, Any],
) -> np.ndarray:
    """Evaluate delta p/p with theta boundary clamping."""
    theta = np.asarray(theta_deg, dtype=np.float64)
    theta_clamped = np.clip(
        theta,
        float(model["theta_valid_min_deg"]),
        float(model["theta_valid_max_deg"]),
    )
    coefficients = np.asarray(
        model["selected_coefficients_ascending"],
        dtype=np.float64,
    )
    correction = np.zeros_like(theta_clamped)
    for power, coefficient in enumerate(coefficients):
        correction += coefficient * theta_clamped**power
    return correction


def evaluate_pion_model(
    theta_deg: np.ndarray,
    model: dict[str, Any],
) -> np.ndarray:
    """Evaluate delta p/p versus centered pion theta with clamping."""
    theta = np.asarray(theta_deg, dtype=np.float64)
    theta_clamped = np.clip(
        theta,
        float(model["theta_valid_min_deg"]),
        float(model["theta_valid_max_deg"]),
    )
    centered = theta_clamped - float(model["theta_reference_deg"])
    coefficients = np.asarray(model["coefficients"], dtype=np.float64)
    correction = np.zeros_like(centered)
    for power, coefficient in enumerate(coefficients):
        correction += coefficient * centered**power
    return correction


def apply_momentum_corrections(
    runnum: np.ndarray,
    detector: np.ndarray,
    e_p: np.ndarray,
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    p_p: np.ndarray,
    p_theta: np.ndarray,
    p_phi: np.ndarray,
    electron_models: dict[tuple[str, int], dict[str, Any]],
    pion_models: dict[tuple[str, int], dict[str, Any]],
) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """
    Correct the electron in all events and the pion only for detector == 1.

    Electron and pion polar/azimuthal directions are unchanged.
    """
    run = np.asarray(runnum, dtype=np.int64)
    detector = np.asarray(detector, dtype=np.int32)
    period_keys = period_key_from_run(run)

    e_sector = fd_sector_from_phi_rad(e_phi)
    pion_sector = fd_sector_from_phi_rad(p_phi)

    e_fraction = np.full(run.shape, np.nan, dtype=np.float64)
    pion_fraction = np.zeros(run.shape, dtype=np.float64)

    for period in CALIBRATION_PERIODS:
        period_mask = period_keys == period.key

        for sector in range(1, 7):
            e_mask = period_mask & (e_sector == sector)
            if np.any(e_mask):
                e_fraction[e_mask] = evaluate_electron_model(
                    np.asarray(e_theta[e_mask]) * RAD2DEG,
                    electron_models[(period.key, sector)],
                )

            pion_mask = (
                period_mask
                & (detector == 1)
                & (pion_sector == sector)
            )
            if np.any(pion_mask):
                pion_fraction[pion_mask] = evaluate_pion_model(
                    np.asarray(p_theta[pion_mask]) * RAD2DEG,
                    pion_models[(period.key, sector)],
                )

    invalid_period = period_keys == ""
    invalid_e_sector = (e_sector < 1) | (e_sector > 6)
    invalid_e_model = ~np.isfinite(e_fraction)

    if np.any(invalid_period):
        example_runs = np.unique(run[invalid_period])[:20].tolist()
        raise RuntimeError(
            "Runs outside the calibrated RGC periods were encountered. "
            f"Examples: {example_runs}"
        )
    if np.any(invalid_e_sector):
        raise RuntimeError(
            f"{int(np.count_nonzero(invalid_e_sector))} electrons could not "
            "be assigned to an FD sector."
        )
    if np.any(invalid_e_model):
        raise RuntimeError(
            f"{int(np.count_nonzero(invalid_e_model))} events did not receive "
            "a finite electron correction."
        )

    corrected_e_p = np.asarray(e_p, dtype=np.float64) * (1.0 + e_fraction)
    corrected_p_p = np.asarray(p_p, dtype=np.float64) * (1.0 + pion_fraction)

    counters = {
        "electron_corrected": int(len(run)),
        "fd_pion_corrected": int(len(run)),
    }
    return corrected_e_p, corrected_p_p, counters


# =============================================================================
# Four-vector and kinematic calculations
# =============================================================================

def momentum_components(
    momentum: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    p = np.asarray(momentum, dtype=np.float64)
    theta = np.asarray(theta, dtype=np.float64)
    phi = np.asarray(phi, dtype=np.float64)
    transverse = p * np.sin(theta)
    return (
        transverse * np.cos(phi),
        transverse * np.sin(phi),
        p * np.cos(theta),
    )


def kallen(x: np.ndarray, y: float, z: float) -> np.ndarray:
    return x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z)


def calculate_mesonic_tmin(
    q2: np.ndarray,
    w: np.ndarray,
) -> np.ndarray:
    """
    Exact forward t for gamma* p -> pi+ n.

    In the gamma*-p center-of-mass frame,

        t = -Q2 + m_pi^2
            - 2 (q0* E_pi* - |q*| |p_pi*| cos(theta*)).

    tmin is the forward value at cos(theta*) = +1 and is therefore the value
    closest to zero in the physical two-body region.
    """
    q2 = np.asarray(q2, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    w2 = w * w

    threshold = PION_MASS_GEV + NEUTRON_MASS_GEV
    valid = (
        np.isfinite(q2)
        & np.isfinite(w)
        & (q2 > 0.0)
        & (w >= threshold)
        & (w2 > 0.0)
    )

    result = np.full(w.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return result

    wv = w[valid]
    w2v = w2[valid]
    q2v = q2[valid]

    q0_star = (
        w2v - PROTON_MASS_GEV**2 - q2v
    ) / (2.0 * wv)
    q_star_sq = np.maximum(q0_star * q0_star + q2v, 0.0)
    q_star = np.sqrt(q_star_sq)

    epi_star = (
        w2v + PION_MASS_GEV**2 - NEUTRON_MASS_GEV**2
    ) / (2.0 * wv)
    ppi_star_sq = np.maximum(
        epi_star * epi_star - PION_MASS_GEV**2,
        0.0,
    )
    ppi_star = np.sqrt(ppi_star_sq)

    result[valid] = (
        -q2v
        + PION_MASS_GEV**2
        - 2.0 * (q0_star * epi_star - q_star * ppi_star)
    )
    return result


def calculate_trento_phi(
    q_vec: np.ndarray,
    electron_vec: np.ndarray,
    pion_vec: np.ndarray,
) -> np.ndarray:
    """
    Calculate the Trento pion azimuth on [0, 2 pi).

    The sign convention reproduces the construction in TwoParticles.java:

        lepton-plane normal = q_hat x e_hat
        hadron-plane normal = q_hat x p_pi
        sign ~ (e x p_pi) dot q_hat

    A boost along q does not alter these plane orientations, so the equivalent
    lab-frame construction is used directly.
    """
    q_norm = np.linalg.norm(q_vec, axis=1)
    qhat = np.divide(
        q_vec,
        q_norm[:, None],
        out=np.full_like(q_vec, np.nan),
        where=q_norm[:, None] > 0.0,
    )

    lepton_normal = np.cross(qhat, electron_vec)
    hadron_normal = np.cross(qhat, pion_vec)
    lepton_norm = np.linalg.norm(lepton_normal, axis=1)
    hadron_norm = np.linalg.norm(hadron_normal, axis=1)

    valid = (
        np.isfinite(qhat).all(axis=1)
        & (lepton_norm > 0.0)
        & (hadron_norm > 0.0)
    )

    result = np.full(q_norm.shape, np.nan, dtype=np.float64)
    if not np.any(valid):
        return result

    n_l = lepton_normal[valid] / lepton_norm[valid, None]
    n_h = hadron_normal[valid] / hadron_norm[valid, None]

    cosine = np.clip(np.einsum("ij,ij->i", n_l, n_h), -1.0, 1.0)
    unsigned = np.arccos(cosine)
    sine_sign = np.einsum(
        "ij,ij->i",
        np.cross(electron_vec[valid], pion_vec[valid]),
        qhat[valid],
    )
    result[valid] = np.where(
        sine_sign < 0.0,
        TWO_PI - unsigned,
        unsigned,
    )
    return result


def calculate_kinematics(
    beam_energy: np.ndarray,
    corrected_e_p: np.ndarray,
    e_theta: np.ndarray,
    e_phi: np.ndarray,
    corrected_p_p: np.ndarray,
    p_theta: np.ndarray,
    p_phi: np.ndarray,
) -> dict[str, np.ndarray]:
    """Recalculate every requested derived branch from corrected four-vectors."""
    beam_energy = np.asarray(beam_energy, dtype=np.float64)
    corrected_e_p = np.asarray(corrected_e_p, dtype=np.float64)
    corrected_p_p = np.asarray(corrected_p_p, dtype=np.float64)

    e_px, e_py, e_pz = momentum_components(
        corrected_e_p, e_theta, e_phi
    )
    pi_px, pi_py, pi_pz = momentum_components(
        corrected_p_p, p_theta, p_phi
    )

    electron_energy = np.sqrt(
        corrected_e_p**2 + ELECTRON_MASS_GEV**2
    )
    pion_energy = np.sqrt(corrected_p_p**2 + PION_MASS_GEV**2)
    beam_momentum = np.sqrt(
        np.maximum(beam_energy**2 - ELECTRON_MASS_GEV**2, 0.0)
    )

    q_energy = beam_energy - electron_energy
    q_px = -e_px
    q_py = -e_py
    q_pz = beam_momentum - e_pz
    q_vec = np.column_stack((q_px, q_py, q_pz))
    electron_vec = np.column_stack((e_px, e_py, e_pz))
    pion_vec = np.column_stack((pi_px, pi_py, pi_pz))

    q_three_sq = q_px**2 + q_py**2 + q_pz**2
    q2 = q_three_sq - q_energy**2

    w2 = (
        (PROTON_MASS_GEV + q_energy) ** 2
        - q_three_sq
    )
    w = np.sqrt(np.maximum(w2, 0.0))

    missing_energy = PROTON_MASS_GEV + q_energy - pion_energy
    missing_px = q_px - pi_px
    missing_py = q_py - pi_py
    missing_pz = q_pz - pi_pz
    mx2 = (
        missing_energy**2
        - missing_px**2
        - missing_py**2
        - missing_pz**2
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        x = q2 / (2.0 * PROTON_MASS_GEV * q_energy)
        y = q_energy / beam_energy

    delta_energy = q_energy - pion_energy
    delta_px = q_px - pi_px
    delta_py = q_py - pi_py
    delta_pz = q_pz - pi_pz
    t = (
        delta_energy**2
        - delta_px**2
        - delta_py**2
        - delta_pz**2
    )
    tmin = calculate_mesonic_tmin(q2, w)
    tprime = t - tmin

    with np.errstate(divide="ignore", invalid="ignore"):
        gamma = 2.0 * PROTON_MASS_GEV * x / np.sqrt(q2)
        gamma2 = gamma * gamma
        common = 1.0 - y - y * y * gamma2 / 4.0

        dep_a = (
            1.0 - y + y * y / 2.0 + y * y * gamma2 / 4.0
        ) / (1.0 + gamma2)
        dep_b = common / (1.0 + gamma2)
        dep_c = y * (1.0 - y / 2.0) / np.sqrt(1.0 + gamma2)

        common_sqrt = np.sqrt(np.maximum(common, 0.0))
        dep_v = (
            (2.0 - y) * common_sqrt / (1.0 + gamma2)
        )
        dep_w = (
            y * common_sqrt / np.sqrt(1.0 + gamma2)
        )

    trento_phi = calculate_trento_phi(q_vec, electron_vec, pion_vec)

    return {
        "Q2": q2,
        "W": w,
        "Mx2": mx2,
        "x": x,
        "y": y,
        "t": t,
        "tmin": tmin,
        "tprime": tprime,
        "phi": trento_phi,
        "DepA": dep_a,
        "DepB": dep_b,
        "DepC": dep_c,
        "DepV": dep_v,
        "DepW": dep_w,
    }


# =============================================================================
# Chunk processing and output
# =============================================================================

def finite_output_mask(output: dict[str, np.ndarray]) -> np.ndarray:
    """Reject events for which any corrected floating branch is non-finite."""
    floating_names = [
        name
        for name, dtype in OUTPUT_BRANCH_DTYPES.items()
        if np.issubdtype(dtype, np.floating)
    ]
    mask = np.ones(len(output["runnum"]), dtype=bool)
    for name in floating_names:
        mask &= np.isfinite(output[name])
    return mask


def cast_and_select(
    output: dict[str, np.ndarray],
    mask: np.ndarray,
) -> dict[str, np.ndarray]:
    return {
        name: np.asarray(output[name][mask], dtype=dtype)
        for name, dtype in OUTPUT_BRANCH_DTYPES.items()
    }


def process_chunk(
    arrays: ak.Array,
    electron_models: dict[tuple[str, int], dict[str, Any]],
    pion_models: dict[tuple[str, int], dict[str, Any]],
    apply_mx2_skim: bool,
    mx2_max: float,
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    raw_all = {
        name: ak.to_numpy(arrays[name])
        for name in REQUIRED_INPUT_BRANCHES
    }

    detector_mask = np.asarray(raw_all["detector"]) == 1
    detector_rejected = int(np.count_nonzero(~detector_mask))
    raw = {
        name: np.asarray(values)[detector_mask]
        for name, values in raw_all.items()
    }

    if len(raw["runnum"]) == 0:
        return (
            {
                name: np.asarray([], dtype=dtype)
                for name, dtype in OUTPUT_BRANCH_DTYPES.items()
            },
            {
                "input": int(len(raw_all["runnum"])),
                "written": 0,
                "detector_rejected": detector_rejected,
                "nonfinite_rejected": 0,
                "mx2_rejected": 0,
                "electron_corrected": 0,
                "fd_pion_corrected": 0,
            },
        )

    corrected_e_p, corrected_p_p, correction_counts = (
        apply_momentum_corrections(
            raw["runnum"],
            raw["detector"],
            raw["e_p"],
            raw["e_theta"],
            raw["e_phi"],
            raw["p_p"],
            raw["p_theta"],
            raw["p_phi"],
            electron_models,
            pion_models,
        )
    )

    beam_energy = beam_energy_from_run(raw["runnum"])
    if np.any(~np.isfinite(beam_energy)):
        bad_runs = np.unique(
            np.asarray(raw["runnum"])[~np.isfinite(beam_energy)]
        )[:20].tolist()
        raise RuntimeError(
            f"Could not assign beam energy to runs {bad_runs}"
        )

    derived = calculate_kinematics(
        beam_energy,
        corrected_e_p,
        raw["e_theta"],
        raw["e_phi"],
        corrected_p_p,
        raw["p_theta"],
        raw["p_phi"],
    )

    output: dict[str, np.ndarray] = {}
    for name in COPIED_INPUT_BRANCHES:
        output[name] = np.asarray(raw[name])

    output["e_phi"] = wrap_phi_rad(output["e_phi"])
    output["p_phi"] = wrap_phi_rad(output["p_phi"])
    output["e_p"] = corrected_e_p
    output["p_p"] = corrected_p_p
    output.update(derived)

    valid = finite_output_mask(output)
    finite_rejected = int(np.count_nonzero(~valid))

    x_mask = (output["x"] > 0.1) & (output["x"] < 0.6)
    x_rejected = int(np.count_nonzero(valid & ~x_mask))
    valid &= x_mask

    if apply_mx2_skim:
        skim_mask = output["Mx2"] < mx2_max
        skim_rejected = int(np.count_nonzero(valid & ~skim_mask))
        valid &= skim_mask
    else:
        skim_rejected = 0

    selected = cast_and_select(output, valid)
    counters = {
        "input": int(len(raw_all["runnum"])),
        "written": int(np.count_nonzero(valid)),
        "detector_rejected": detector_rejected,
        "nonfinite_rejected": finite_rejected,
        "x_rejected": x_rejected,
        "mx2_rejected": skim_rejected,
        **correction_counts,
    }
    return selected, counters


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select detector == 1 pions and apply RGC electron and FD pi+ "
            "momentum corrections while rewriting corrected kinematics."
        )
    )
    parser.add_argument("input_root", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("--tree-name", default="PhysicsEvents")
    parser.add_argument(
        "--output-tree-name",
        default=None,
        help="Default: use --tree-name.",
    )
    parser.add_argument(
        "--electron-correction-json",
        type=Path,
        default=DEFAULT_ELECTRON_JSON,
    )
    parser.add_argument(
        "--pion-correction-json",
        type=Path,
        default=DEFAULT_PION_JSON,
    )
    parser.add_argument(
        "--step-size",
        default="250 MB",
        help="uproot iteration step size; default: 250 MB.",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=4,
        help="ZLIB compression level in [0, 9]; default: 4.",
    )
    parser.add_argument(
        "--mx2-max",
        type=float,
        default=1.5,
        help="Default corrected-Mx2 skim upper limit in GeV^2.",
    )
    parser.add_argument(
        "--no-mx2-skim",
        action="store_true",
        help="Disable the default corrected Mx2 < 1.5 GeV^2 skim.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Permit replacement of an existing output file.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    if not args.input_root.is_file():
        raise FileNotFoundError(args.input_root)
    if args.input_root.resolve() == args.output_root.resolve():
        raise ValueError("Input and output ROOT paths must differ.")
    if args.output_root.exists() and not args.overwrite:
        raise FileExistsError(
            f"{args.output_root} already exists. Use --overwrite to replace it."
        )
    if not 0 <= args.compression_level <= 9:
        raise ValueError("--compression-level must be in [0, 9].")
    if not np.isfinite(args.mx2_max):
        raise ValueError("--mx2-max must be finite.")

    input_tree_name = find_tree(args.input_root, args.tree_name)
    output_tree_name = args.output_tree_name or args.tree_name
    validate_input_branches(args.input_root, input_tree_name)

    electron_models = load_electron_models(args.electron_correction_json)
    pion_models = load_pion_models(args.pion_correction_json)

    args.output_root.parent.mkdir(parents=True, exist_ok=True)
    if args.output_root.exists():
        args.output_root.unlink()

    totals = {
        "input": 0,
        "written": 0,
        "detector_rejected": 0,
        "nonfinite_rejected": 0,
        "x_rejected": 0,
        "mx2_rejected": 0,
        "electron_corrected": 0,
        "fd_pion_corrected": 0,
    }

    source = f"{args.input_root}:{input_tree_name}"
    expressions = list(REQUIRED_INPUT_BRANCHES)

    with uproot.recreate(
        args.output_root,
        compression=uproot.ZLIB(args.compression_level),
    ) as output_file:
        output_tree = output_file.mktree(
            output_tree_name,
            OUTPUT_BRANCH_DTYPES,
            title=(
                "RGC e pi+ X data with electron/FD-pion momentum corrections "
                "and recalculated mesonic kinematics including tprime = t - tmin"
            ),
        )

        for chunk_index, arrays in enumerate(
            uproot.iterate(
                source,
                expressions=expressions,
                library="ak",
                step_size=args.step_size,
            ),
            start=1,
        ):
            selected, counters = process_chunk(
                arrays,
                electron_models,
                pion_models,
                apply_mx2_skim=not args.no_mx2_skim,
                mx2_max=args.mx2_max,
            )
            if counters["written"] > 0:
                output_tree.extend(selected)

            for key in totals:
                totals[key] += counters[key]

            print(
                f"[chunk {chunk_index}] "
                f"input={counters['input']:,}, "
                f"written={counters['written']:,}, "
                f"nonfinite={counters['nonfinite_rejected']:,}, "
                f"x skim={counters['x_rejected']:,}, "
                f"Mx2 skim={counters['mx2_rejected']:,}"
            )

    print("")
    print("Momentum-corrected e pi+ X tree complete.")
    print(f"Input:  {args.input_root}:{input_tree_name}")
    print(f"Output: {args.output_root}:{output_tree_name}")
    print(f"Input events:                 {totals['input']:,}")
    print(f"Non-FD pion events rejected:  {totals['detector_rejected']:,}")
    print(f"Written events:               {totals['written']:,}")
    print(f"Non-finite rejected:          {totals['nonfinite_rejected']:,}")
    print(
        f"x <= 0.1 or x >= 0.6 rejected: "
        f"{totals['x_rejected']:,}"
    )
    if args.no_mx2_skim:
        print("Corrected-Mx2 skim:           disabled")
    else:
        print(
            f"Corrected Mx2 >= {args.mx2_max:g} rejected: "
            f"{totals['mx2_rejected']:,}"
        )
    print(f"Electron-corrected events:    {totals['electron_corrected']:,}")
    print(f"FD-pion-corrected events:     {totals['fd_pion_corrected']:,}")
    print(
        f"Non-FD pion momentum unchanged: "
        f"{totals['non_fd_pion_unchanged']:,}"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(
            f"FATAL: {type(exc).__name__}: {exc}",
            file=sys.stderr,
        )
        raise
