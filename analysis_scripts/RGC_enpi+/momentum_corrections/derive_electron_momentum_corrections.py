#!/usr/bin/env python3
"""
derive_electron_momentum_corrections_v1.py

Diagnostic and provisional calibration extraction for RGC electron momentum
corrections using inclusive NH3 eX data.

The input samples are not assumed to be pure elastic samples. In every
calibration cell, the script fits the reconstructed W distribution with a
Gaussian elastic core plus a quadratic background. The fitted elastic-peak
position is then translated into the momentum-scale correction required to
place the peak at the proton mass.

Calibration configurations:
    Su22                 runs 16042-16788, solenoid configuration from data period
    Fa22, solenoid -1    runs 16843-17183
    Fa22, solenoid +1    runs 17185-17408
    Sp23                 runs 17477-17811

All considered data are torus inbending. Electrons are assigned to the six
CLAS12 Forward Detector sectors using their azimuthal angle. The sector-local
azimuth is represented on [0, 60) degrees by rotating every sector onto sector 1.

Primary outputs:
    output/electron_diagnostics_v1/fit_results.csv
    output/electron_diagnostics_v1/electron_correction_cells_v1.json
    output/electron_diagnostics_v1/electron_correction_models_v1.json
    output/electron_diagnostics_v1/plots/*.png

The JSON output is provisional. It records the fitted calibration cells and a
simple weighted polynomial model, but it is not yet intended to be consumed by
the production analysis until the diagnostic plots and model choices have been
reviewed.

Dependencies:
    numpy
    pandas
    scipy
    matplotlib
    uproot
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit

# -----------------------------------------------------------------------------
# Physics and analysis constants
# -----------------------------------------------------------------------------

PROTON_MASS_GEV = 0.9382720813
RAD2DEG = 180.0 / math.pi

DEFAULT_INPUTS = {
    "su22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/elastic/"
        "rgc_su22_inb_eX_Wlessthan2.root"
    ),
    "fa22": Path(
        "/work/clas12/thayward/CLAS12_exclusive/elastic/"
        "rgc_fa22_inb_eX_Wlessthan2.root"
    ),
    "sp23": Path(
        "/work/clas12/thayward/CLAS12_exclusive/elastic/"
        "rgc_sp23_inb_eX_Wlessthan2.root"
    ),
}

REQUIRED_BRANCHES = (
    "fiducial_status",
    "runnum",
    "evnum",
    "e_p",
    "e_theta",
    "e_phi",
    "vz_e",
    "W",
)

OPTIONAL_BRANCHES = (
    "Q2",
    "x",
    "y",
    "num_pos",
    "num_neg",
    "num_neutral",
)


@dataclass(frozen=True)
class CalibrationPeriod:
    key: str
    label: str
    source_key: str
    run_min: int
    run_max: int
    beam_energy_gev: float | None
    torus_polarity: int
    solenoid_polarity: int | None


CALIBRATION_PERIODS = (
    CalibrationPeriod(
        key="su22",
        label="Su22",
        source_key="su22",
        run_min=16042,
        run_max=16788,
        beam_energy_gev=None,
        torus_polarity=-1,
        solenoid_polarity=None,
    ),
    CalibrationPeriod(
        key="fa22_solenoid_minus",
        label="Fa22, solenoid -1",
        source_key="fa22",
        run_min=16843,
        run_max=17183,
        beam_energy_gev=None,
        torus_polarity=-1,
        solenoid_polarity=-1,
    ),
    CalibrationPeriod(
        key="fa22_solenoid_plus",
        label="Fa22, solenoid +1",
        source_key="fa22",
        run_min=17185,
        run_max=17408,
        beam_energy_gev=None,
        torus_polarity=-1,
        solenoid_polarity=+1,
    ),
    CalibrationPeriod(
        key="sp23",
        label="Sp23",
        source_key="sp23",
        run_min=17477,
        run_max=17811,
        beam_energy_gev=None,
        torus_polarity=-1,
        solenoid_polarity=None,
    ),
)


# -----------------------------------------------------------------------------
# General utilities
# -----------------------------------------------------------------------------

def parse_edges(text: str, name: str) -> np.ndarray:
    """Parse a comma-separated list of monotonically increasing bin edges."""
    try:
        edges = np.asarray([float(value.strip()) for value in text.split(",")])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{name} must be a comma-separated list of numbers."
        ) from exc
    # endtry

    if edges.size < 2:
        raise argparse.ArgumentTypeError(f"{name} must contain at least two edges.")
    # endif

    if not np.all(np.diff(edges) > 0.0):
        raise argparse.ArgumentTypeError(
            f"{name} edges must be strictly increasing."
        )
    # endif

    return edges


def wrap_phi_deg(phi_deg: np.ndarray | float) -> np.ndarray:
    """Wrap azimuthal angles to [0, 360) degrees."""
    return np.mod(np.asarray(phi_deg, dtype=float), 360.0)


def fd_sector_from_phi_rad(phi_rad: np.ndarray) -> np.ndarray:
    """
    Assign CLAS12 Forward Detector sectors from azimuthal angle.

    Sector definitions:
        1: [330, 360) U [0, 30)
        2: [30, 90)
        3: [90, 150)
        4: [150, 210)
        5: [210, 270)
        6: [270, 330)
    """
    phi_deg = wrap_phi_deg(np.asarray(phi_rad, dtype=float) * RAD2DEG)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)

    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6

    return sector


def fd_local_phi_deg(phi_rad: np.ndarray, sector: np.ndarray) -> np.ndarray:
    """
    Rotate each FD sector onto sector 1 and return local phi on [0, 60).

    The sector lower edges are:
        sector 1: 330 deg
        sector 2:  30 deg
        sector 3:  90 deg
        sector 4: 150 deg
        sector 5: 210 deg
        sector 6: 270 deg
    """
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


def beam_energy_from_run(runnum: np.ndarray) -> np.ndarray:
    """
    Return the nominal RGC beam energy from run number.

    This broader mapping follows the run ranges supplied by the analysis:
        16042-17065: 10.5473 GeV
        17067-17716: 10.5563 GeV
        17717-17811: 10.5593 GeV

    The calibration-period records below use the requested period-level
    energies. This helper is retained for validation and future event-level use.
    """
    runnum = np.asarray(runnum, dtype=int)
    energy = np.full(runnum.shape, np.nan, dtype=float)

    energy[(runnum >= 16042) & (runnum <= 17065)] = 10.5473
    energy[(runnum >= 17067) & (runnum <= 17716)] = 10.5563
    energy[(runnum >= 17717) & (runnum <= 17811)] = 10.5593

    return energy


def elastic_scattered_momentum(theta_rad: np.ndarray, beam_energy_gev: float) -> np.ndarray:
    """
    Expected scattered-electron momentum for elastic ep scattering.

    The electron mass is neglected:
        p'_elastic = M_p E_b / [M_p + E_b (1 - cos(theta_e))]
    """
    theta_rad = np.asarray(theta_rad, dtype=float)
    denominator = PROTON_MASS_GEV + beam_energy_gev * (1.0 - np.cos(theta_rad))
    return PROTON_MASS_GEV * beam_energy_gev / denominator


def w_from_electron(
    e_p_gev: np.ndarray,
    e_theta_rad: np.ndarray,
    beam_energy_gev: float,
) -> np.ndarray:
    """
    Recalculate W from the measured electron momentum and polar angle.

    The electron mass is neglected:
        W^2 = M_p^2 + 2 M_p (E_b - E') - Q^2
        Q^2 = 2 E_b E' (1 - cos(theta_e))
    """
    e_p_gev = np.asarray(e_p_gev, dtype=float)
    e_theta_rad = np.asarray(e_theta_rad, dtype=float)

    q2 = 2.0 * beam_energy_gev * e_p_gev * (1.0 - np.cos(e_theta_rad))
    w2 = (
        PROTON_MASS_GEV**2
        + 2.0 * PROTON_MASS_GEV * (beam_energy_gev - e_p_gev)
        - q2
    )
    return np.sqrt(np.clip(w2, a_min=0.0, a_max=None))


def momentum_at_w(
    w_gev: float,
    theta_rad: float,
    beam_energy_gev: float,
) -> float:
    """
    Electron momentum corresponding to a specified W at fixed beam energy and angle.

        p' = [M_p^2 + 2 M_p E_b - W^2]
             / [2 (M_p + E_b (1 - cos(theta_e)))]
    """
    numerator = (
        PROTON_MASS_GEV**2
        + 2.0 * PROTON_MASS_GEV * beam_energy_gev
        - w_gev**2
    )
    denominator = 2.0 * (
        PROTON_MASS_GEV + beam_energy_gev * (1.0 - math.cos(theta_rad))
    )
    return numerator / denominator


# -----------------------------------------------------------------------------
# W-peak fitting
# -----------------------------------------------------------------------------

def gaussian_plus_quadratic(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    sigma: float,
    b0: float,
    b1: float,
    b2: float,
) -> np.ndarray:
    """Gaussian elastic peak plus a quadratic background."""
    sigma_safe = np.maximum(np.abs(sigma), 1.0e-6)
    gaussian = amplitude * np.exp(-0.5 * ((x - mean) / sigma_safe) ** 2)
    x_centered = x - PROTON_MASS_GEV
    background = b0 + b1 * x_centered + b2 * x_centered**2
    return gaussian + background


@dataclass
class PeakFitResult:
    success: bool
    status: str
    n_events: int
    n_histogram_entries: int
    peak_mean_gev: float
    peak_mean_error_gev: float
    peak_sigma_gev: float
    peak_sigma_error_gev: float
    amplitude: float
    chi2: float
    ndf: int
    chi2_ndf: float
    fit_range_low_gev: float
    fit_range_high_gev: float


def failed_peak_fit(
    status: str,
    n_events: int,
    n_histogram_entries: int = 0,
    fit_range: tuple[float, float] = (math.nan, math.nan),
) -> PeakFitResult:
    """Construct a failed fit record with NaN numerical fields."""
    return PeakFitResult(
        success=False,
        status=status,
        n_events=n_events,
        n_histogram_entries=n_histogram_entries,
        peak_mean_gev=math.nan,
        peak_mean_error_gev=math.nan,
        peak_sigma_gev=math.nan,
        peak_sigma_error_gev=math.nan,
        amplitude=math.nan,
        chi2=math.nan,
        ndf=0,
        chi2_ndf=math.nan,
        fit_range_low_gev=fit_range[0],
        fit_range_high_gev=fit_range[1],
    )


def fit_w_peak(
    w_values: np.ndarray,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
) -> tuple[PeakFitResult, dict[str, np.ndarray]]:
    """
    Fit the elastic peak in a W distribution.

    The fit uses Poisson-like uncertainties sqrt(max(N, 1)) on histogram bins.
    A robust mode estimate initializes the Gaussian mean. The Gaussian width is
    bounded to remain inside a physically sensible core range.
    """
    finite_w = np.asarray(w_values, dtype=float)
    finite_w = finite_w[np.isfinite(finite_w)]
    n_events = int(finite_w.size)

    histogram_counts, histogram_edges = np.histogram(
        finite_w,
        bins=histogram_bins,
        range=histogram_range,
    )
    histogram_centers = 0.5 * (histogram_edges[:-1] + histogram_edges[1:])

    diagnostics = {
        "counts": histogram_counts,
        "edges": histogram_edges,
        "centers": histogram_centers,
        "fit_x": np.asarray([], dtype=float),
        "fit_y": np.asarray([], dtype=float),
    }

    if n_events < minimum_events:
        return (
            failed_peak_fit(
                status="insufficient_events",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
            ),
            diagnostics,
        )
    # endif

    fit_mask = (
        (histogram_centers >= fit_range[0])
        & (histogram_centers <= fit_range[1])
    )
    x_fit = histogram_centers[fit_mask]
    y_fit = histogram_counts[fit_mask].astype(float)

    if x_fit.size < 10 or np.count_nonzero(y_fit) < 8:
        return (
            failed_peak_fit(
                status="insufficient_populated_bins",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
            ),
            diagnostics,
        )
    # endif

    peak_search_mask = (
        (x_fit >= max(fit_range[0], PROTON_MASS_GEV - 0.10))
        & (x_fit <= min(fit_range[1], PROTON_MASS_GEV + 0.10))
    )
    if not np.any(peak_search_mask):
        return (
            failed_peak_fit(
                status="empty_peak_search_window",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
            ),
            diagnostics,
        )
    # endif

    peak_search_indices = np.flatnonzero(peak_search_mask)
    peak_index = peak_search_indices[np.argmax(y_fit[peak_search_mask])]
    initial_mean = float(x_fit[peak_index])

    edge_mask = (
        (x_fit < initial_mean - 0.06)
        | (x_fit > initial_mean + 0.06)
    )
    if np.count_nonzero(edge_mask) >= 3:
        background_initial = float(np.median(y_fit[edge_mask]))
    else:
        background_initial = float(np.percentile(y_fit, 20.0))
    # endif

    initial_amplitude = max(float(y_fit[peak_index]) - background_initial, 1.0)
    initial_sigma = 0.025

    initial_parameters = (
        initial_amplitude,
        initial_mean,
        initial_sigma,
        max(background_initial, 0.0),
        0.0,
        0.0,
    )

    lower_bounds = (
        0.0,
        PROTON_MASS_GEV - 0.10,
        0.005,
        -np.inf,
        -np.inf,
        -np.inf,
    )
    upper_bounds = (
        np.inf,
        PROTON_MASS_GEV + 0.10,
        0.10,
        np.inf,
        np.inf,
        np.inf,
    )

    sigma_y = np.sqrt(np.maximum(y_fit, 1.0))

    try:
        parameters, covariance = curve_fit(
            gaussian_plus_quadratic,
            x_fit,
            y_fit,
            p0=initial_parameters,
            sigma=sigma_y,
            absolute_sigma=True,
            bounds=(lower_bounds, upper_bounds),
            maxfev=50_000,
        )
    except (RuntimeError, ValueError, FloatingPointError) as exc:
        return (
            failed_peak_fit(
                status=f"fit_failed:{type(exc).__name__}",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
            ),
            diagnostics,
        )
    # endtry

    parameter_errors = np.sqrt(np.clip(np.diag(covariance), 0.0, None))
    fitted_y = gaussian_plus_quadratic(x_fit, *parameters)
    chi2 = float(np.sum(((y_fit - fitted_y) / sigma_y) ** 2))
    ndf = int(x_fit.size - len(parameters))
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan

    diagnostics["fit_x"] = np.linspace(fit_range[0], fit_range[1], 800)
    diagnostics["fit_y"] = gaussian_plus_quadratic(
        diagnostics["fit_x"],
        *parameters,
    )

    result = PeakFitResult(
        success=True,
        status="success",
        n_events=n_events,
        n_histogram_entries=int(histogram_counts.sum()),
        peak_mean_gev=float(parameters[1]),
        peak_mean_error_gev=float(parameter_errors[1]),
        peak_sigma_gev=float(abs(parameters[2])),
        peak_sigma_error_gev=float(parameter_errors[2]),
        amplitude=float(parameters[0]),
        chi2=chi2,
        ndf=ndf,
        chi2_ndf=chi2_ndf,
        fit_range_low_gev=fit_range[0],
        fit_range_high_gev=fit_range[1],
    )
    return result, diagnostics


# -----------------------------------------------------------------------------
# ROOT input
# -----------------------------------------------------------------------------

def find_tree(file_path: Path, requested_tree_name: str) -> str:
    """Find the requested tree or fall back to the only TTree in the file."""
    with uproot.open(file_path) as root_file:
        if requested_tree_name in root_file:
            return requested_tree_name
        # endif

        tree_candidates: list[str] = []
        for key, class_name in root_file.classnames().items():
            if class_name.startswith("TTree"):
                tree_candidates.append(key.split(";")[0])
            # endif
        # endfor

    unique_candidates = sorted(set(tree_candidates))
    if len(unique_candidates) == 1:
        return unique_candidates[0]
    # endif

    raise RuntimeError(
        f"Could not identify tree '{requested_tree_name}' in {file_path}. "
        f"TTree candidates: {unique_candidates}"
    )


def validate_branches(file_path: Path, tree_name: str) -> tuple[str, ...]:
    """Validate required branches and return all requested available branches."""
    with uproot.open(file_path) as root_file:
        tree = root_file[tree_name]
        available = set(tree.keys())
    # endwith

    missing = [branch for branch in REQUIRED_BRANCHES if branch not in available]
    if missing:
        raise RuntimeError(
            f"{file_path}:{tree_name} is missing required branches: {missing}"
        )
    # endif

    selected = list(REQUIRED_BRANCHES)
    selected.extend(branch for branch in OPTIONAL_BRANCHES if branch in available)
    return tuple(selected)


def read_selected_events(
    file_path: Path,
    tree_name: str,
    branches: tuple[str, ...],
    period: CalibrationPeriod,
    step_size: str,
    require_fiducial_status: int | None,
    vertex_min_cm: float,
    vertex_max_cm: float,
    theta_min_deg: float,
    theta_max_deg: float,
    w_preselection_min_gev: float,
    w_preselection_max_gev: float,
) -> pd.DataFrame:
    """
    Read and select events for one calibration period.

    The ROOT file is processed in chunks. Only selected NumPy arrays are retained.
    """
    selected_chunks: list[pd.DataFrame] = []
    source_expression = f"{file_path}:{tree_name}"

    for arrays in uproot.iterate(
        source_expression,
        expressions=list(branches),
        library="np",
        step_size=step_size,
    ):
        runnum = arrays["runnum"].astype(np.int32, copy=False)
        e_p = arrays["e_p"].astype(float, copy=False)
        e_theta = arrays["e_theta"].astype(float, copy=False)
        e_phi = arrays["e_phi"].astype(float, copy=False)
        vz_e = arrays["vz_e"].astype(float, copy=False)
        w_branch = arrays["W"].astype(float, copy=False)

        mask = (
            (runnum >= period.run_min)
            & (runnum <= period.run_max)
            & np.isfinite(e_p)
            & np.isfinite(e_theta)
            & np.isfinite(e_phi)
            & np.isfinite(vz_e)
            & np.isfinite(w_branch)
            & (e_p > 0.0)
            & (e_theta * RAD2DEG >= theta_min_deg)
            & (e_theta * RAD2DEG < theta_max_deg)
            & (vz_e >= vertex_min_cm)
            & (vz_e <= vertex_max_cm)
            & (w_branch >= w_preselection_min_gev)
            & (w_branch <= w_preselection_max_gev)
        )

        if require_fiducial_status is not None:
            mask &= arrays["fiducial_status"] == require_fiducial_status
        # endif

        if not np.any(mask):
            continue
        # endif

        sector = fd_sector_from_phi_rad(e_phi[mask])
        local_phi = fd_local_phi_deg(e_phi[mask], sector)
        valid_sector = (
            (sector >= 1)
            & (sector <= 6)
            & np.isfinite(local_phi)
        )

        if not np.any(valid_sector):
            continue
        # endif

        selected_indices = np.flatnonzero(mask)[valid_sector]
        selected_sector = sector[valid_sector]
        selected_local_phi = local_phi[valid_sector]

        frame_data: dict[str, Any] = {
            "runnum": runnum[selected_indices],
            "evnum": arrays["evnum"][selected_indices].astype(np.int64, copy=False),
            "fiducial_status": arrays["fiducial_status"][selected_indices].astype(
                np.int16,
                copy=False,
            ),
            "e_p": e_p[selected_indices],
            "e_theta_rad": e_theta[selected_indices],
            "e_theta_deg": e_theta[selected_indices] * RAD2DEG,
            "e_phi_rad": e_phi[selected_indices],
            "e_phi_deg": wrap_phi_deg(e_phi[selected_indices] * RAD2DEG),
            "vz_e": vz_e[selected_indices],
            "w_branch": w_branch[selected_indices],
            "sector": selected_sector,
            "local_phi_deg": selected_local_phi,
        }

        for branch in OPTIONAL_BRANCHES:
            if branch in arrays:
                frame_data[branch] = arrays[branch][selected_indices]
            # endif
        # endfor

        chunk_frame = pd.DataFrame(frame_data)
        chunk_frame["beam_energy_gev"] = beam_energy_from_run(
            chunk_frame["runnum"].to_numpy()
        )

        valid_beam_energy = np.isfinite(
            chunk_frame["beam_energy_gev"].to_numpy()
        )
        if not np.all(valid_beam_energy):
            invalid_runs = sorted(
                set(
                    chunk_frame.loc[
                        ~valid_beam_energy,
                        "runnum",
                    ].astype(int)
                )
            )
            raise RuntimeError(
                f"Undefined beam energy for runs in {period.label}: "
                f"{invalid_runs}"
            )
        # endif

        theta_values = chunk_frame["e_theta_rad"].to_numpy()
        momentum_values = chunk_frame["e_p"].to_numpy()
        beam_energy_values = chunk_frame["beam_energy_gev"].to_numpy()

        denominator = (
            PROTON_MASS_GEV
            + beam_energy_values * (1.0 - np.cos(theta_values))
        )
        chunk_frame["p_elastic_gev"] = (
            PROTON_MASS_GEV * beam_energy_values / denominator
        )
        chunk_frame["raw_delta_p_over_p"] = (
            chunk_frame["p_elastic_gev"] - chunk_frame["e_p"]
        ) / chunk_frame["e_p"]

        q2_values = (
            2.0
            * beam_energy_values
            * momentum_values
            * (1.0 - np.cos(theta_values))
        )
        w2_values = (
            PROTON_MASS_GEV**2
            + 2.0
            * PROTON_MASS_GEV
            * (beam_energy_values - momentum_values)
            - q2_values
        )
        chunk_frame["w_recalculated"] = np.sqrt(
            np.clip(w2_values, a_min=0.0, a_max=None)
        )

        selected_chunks.append(chunk_frame)
    # endfor

    if not selected_chunks:
        return pd.DataFrame()
    # endif

    return pd.concat(selected_chunks, ignore_index=True)


# -----------------------------------------------------------------------------
# Calibration-cell construction
# -----------------------------------------------------------------------------

def correction_from_peak(
    peak_w_gev: float,
    peak_w_error_gev: float,
    mean_theta_rad: float,
    mean_e_p_gev: float,
    beam_energy_gev: float,
) -> tuple[float, float, float]:
    """
    Translate a fitted W peak into an electron momentum correction.

    At the mean polar angle of the calibration cell:
        p_at_peak  = momentum corresponding to the fitted W peak
        p_elastic  = momentum corresponding to W = M_p

    The relative correction is:
        delta_p_over_p = (p_elastic - p_at_peak) / mean_measured_p

    The uncertainty is propagated numerically from the W-peak uncertainty.
    """
    p_at_peak = momentum_at_w(
        peak_w_gev,
        mean_theta_rad,
        beam_energy_gev,
    )
    p_elastic = elastic_scattered_momentum(
        np.asarray([mean_theta_rad]),
        beam_energy_gev,
    )[0]
    delta_p_over_p = (p_elastic - p_at_peak) / mean_e_p_gev

    if not np.isfinite(peak_w_error_gev) or peak_w_error_gev <= 0.0:
        return float(delta_p_over_p), math.nan, float(p_at_peak)
    # endif

    peak_high = peak_w_gev + peak_w_error_gev
    peak_low = peak_w_gev - peak_w_error_gev

    correction_high = (
        p_elastic
        - momentum_at_w(peak_high, mean_theta_rad, beam_energy_gev)
    ) / mean_e_p_gev
    correction_low = (
        p_elastic
        - momentum_at_w(peak_low, mean_theta_rad, beam_energy_gev)
    ) / mean_e_p_gev

    correction_error = 0.5 * abs(correction_high - correction_low)
    return (
        float(delta_p_over_p),
        float(correction_error),
        float(p_at_peak),
    )


def make_fit_record(
    period: CalibrationPeriod,
    sector: int,
    cell_type: str,
    theta_bin_index: int | None,
    theta_low_deg: float,
    theta_high_deg: float,
    local_phi_bin_index: int | None,
    local_phi_low_deg: float,
    local_phi_high_deg: float,
    cell_frame: pd.DataFrame,
    fit_result: PeakFitResult,
) -> dict[str, Any]:
    """Combine event-cell metadata with the W-peak fit result."""
    if cell_frame.empty:
        mean_theta_rad = math.nan
        mean_theta_deg = math.nan
        mean_local_phi_deg = math.nan
        mean_e_p_gev = math.nan
        mean_run = math.nan
    else:
        mean_theta_rad = float(cell_frame["e_theta_rad"].mean())
        mean_theta_deg = float(cell_frame["e_theta_deg"].mean())
        mean_local_phi_deg = float(cell_frame["local_phi_deg"].mean())
        mean_e_p_gev = float(cell_frame["e_p"].mean())
        mean_run = float(cell_frame["runnum"].mean())
    # endif

    correction = math.nan
    correction_error = math.nan
    p_at_peak = math.nan

    if (
        fit_result.success
        and np.isfinite(mean_theta_rad)
        and np.isfinite(mean_e_p_gev)
        and mean_e_p_gev > 0.0
    ):
        correction, correction_error, p_at_peak = correction_from_peak(
            peak_w_gev=fit_result.peak_mean_gev,
            peak_w_error_gev=fit_result.peak_mean_error_gev,
            mean_theta_rad=mean_theta_rad,
            mean_e_p_gev=mean_e_p_gev,
            beam_energy_gev=float(cell_frame["beam_energy_gev"].mean()),
        )
    # endif

    record = {
        "period_key": period.key,
        "period_label": period.label,
        "source_key": period.source_key,
        "run_min": period.run_min,
        "run_max": period.run_max,
        "beam_energy_gev": (
            float(cell_frame["beam_energy_gev"].mean())
            if not cell_frame.empty
            else math.nan
        ),
        "torus_polarity": period.torus_polarity,
        "solenoid_polarity": period.solenoid_polarity,
        "sector": sector,
        "cell_type": cell_type,
        "theta_bin_index": theta_bin_index,
        "theta_low_deg": theta_low_deg,
        "theta_high_deg": theta_high_deg,
        "local_phi_bin_index": local_phi_bin_index,
        "local_phi_low_deg": local_phi_low_deg,
        "local_phi_high_deg": local_phi_high_deg,
        "mean_run": mean_run,
        "mean_theta_rad": mean_theta_rad,
        "mean_theta_deg": mean_theta_deg,
        "mean_local_phi_deg": mean_local_phi_deg,
        "mean_e_p_gev": mean_e_p_gev,
        "p_at_fitted_peak_gev": p_at_peak,
        "delta_p_over_p": correction,
        "delta_p_over_p_error": correction_error,
    }
    record.update(asdict(fit_result))
    return record


def fit_calibration_cells(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    theta_edges_deg: np.ndarray,
    local_phi_edges_deg: np.ndarray,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    save_individual_fits: bool,
    individual_fit_directory: Path,
) -> list[dict[str, Any]]:
    """Fit integrated, theta-binned, and theta/local-phi-binned W spectra."""
    records: list[dict[str, Any]] = []

    for sector in range(1, 7):
        sector_frame = frame.loc[frame["sector"] == sector].copy()

        integrated_fit, integrated_diagnostics = fit_w_peak(
            sector_frame["w_recalculated"].to_numpy(),
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events_integrated,
        )
        records.append(
            make_fit_record(
                period=period,
                sector=sector,
                cell_type="sector_integrated",
                theta_bin_index=None,
                theta_low_deg=float(theta_edges_deg[0]),
                theta_high_deg=float(theta_edges_deg[-1]),
                local_phi_bin_index=None,
                local_phi_low_deg=float(local_phi_edges_deg[0]),
                local_phi_high_deg=float(local_phi_edges_deg[-1]),
                cell_frame=sector_frame,
                fit_result=integrated_fit,
            )
        )

        if save_individual_fits:
            save_peak_fit_plot(
                output_path=individual_fit_directory
                / f"{period.key}_sector{sector}_integrated.png",
                diagnostics=integrated_diagnostics,
                fit_result=integrated_fit,
                title=f"{period.label}, sector {sector}, integrated",
            )
        # endif

        for theta_index in range(len(theta_edges_deg) - 1):
            theta_low = float(theta_edges_deg[theta_index])
            theta_high = float(theta_edges_deg[theta_index + 1])
            theta_frame = sector_frame.loc[
                (sector_frame["e_theta_deg"] >= theta_low)
                & (sector_frame["e_theta_deg"] < theta_high)
            ].copy()

            theta_fit, theta_diagnostics = fit_w_peak(
                theta_frame["w_recalculated"].to_numpy(),
                histogram_range=histogram_range,
                fit_range=fit_range,
                histogram_bins=histogram_bins,
                minimum_events=minimum_events_cell,
            )
            records.append(
                make_fit_record(
                    period=period,
                    sector=sector,
                    cell_type="theta",
                    theta_bin_index=theta_index,
                    theta_low_deg=theta_low,
                    theta_high_deg=theta_high,
                    local_phi_bin_index=None,
                    local_phi_low_deg=float(local_phi_edges_deg[0]),
                    local_phi_high_deg=float(local_phi_edges_deg[-1]),
                    cell_frame=theta_frame,
                    fit_result=theta_fit,
                )
            )

            if save_individual_fits:
                save_peak_fit_plot(
                    output_path=individual_fit_directory
                    / (
                        f"{period.key}_sector{sector}_theta"
                        f"{theta_index:02d}.png"
                    ),
                    diagnostics=theta_diagnostics,
                    fit_result=theta_fit,
                    title=(
                        f"{period.label}, sector {sector}, "
                        f"{theta_low:.1f} <= theta_e < {theta_high:.1f} deg"
                    ),
                )
            # endif

            for phi_index in range(len(local_phi_edges_deg) - 1):
                phi_low = float(local_phi_edges_deg[phi_index])
                phi_high = float(local_phi_edges_deg[phi_index + 1])
                cell_frame = theta_frame.loc[
                    (theta_frame["local_phi_deg"] >= phi_low)
                    & (theta_frame["local_phi_deg"] < phi_high)
                ].copy()

                cell_fit, cell_diagnostics = fit_w_peak(
                    cell_frame["w_recalculated"].to_numpy(),
                    histogram_range=histogram_range,
                    fit_range=fit_range,
                    histogram_bins=histogram_bins,
                    minimum_events=minimum_events_cell,
                )
                records.append(
                    make_fit_record(
                        period=period,
                        sector=sector,
                        cell_type="theta_local_phi",
                        theta_bin_index=theta_index,
                        theta_low_deg=theta_low,
                        theta_high_deg=theta_high,
                        local_phi_bin_index=phi_index,
                        local_phi_low_deg=phi_low,
                        local_phi_high_deg=phi_high,
                        cell_frame=cell_frame,
                        fit_result=cell_fit,
                    )
                )

                if save_individual_fits:
                    save_peak_fit_plot(
                        output_path=individual_fit_directory
                        / (
                            f"{period.key}_sector{sector}_theta"
                            f"{theta_index:02d}_phi{phi_index:02d}.png"
                        ),
                        diagnostics=cell_diagnostics,
                        fit_result=cell_fit,
                        title=(
                            f"{period.label}, sector {sector}, "
                            f"{theta_low:.1f} <= theta_e < {theta_high:.1f} deg, "
                            f"{phi_low:.1f} <= phi_local < {phi_high:.1f} deg"
                        ),
                    )
                # endif
            # endfor
        # endfor
    # endfor

    return records


# -----------------------------------------------------------------------------
# Provisional correction-surface model
# -----------------------------------------------------------------------------

MODEL_TERMS = (
    "constant",
    "p",
    "p2",
    "phi_centered",
    "phi_centered2",
    "p_times_phi_centered",
)


def correction_design_matrix(
    momentum_gev: np.ndarray,
    local_phi_deg: np.ndarray,
) -> np.ndarray:
    """
    Construct a controlled six-parameter correction model.

        delta_p/p =
            c0
          + c1 p
          + c2 p^2
          + c3 phi_c
          + c4 phi_c^2
          + c5 p phi_c

    where phi_c = phi_local - 30 degrees.
    """
    p = np.asarray(momentum_gev, dtype=float)
    phi_centered = np.asarray(local_phi_deg, dtype=float) - 30.0

    return np.column_stack(
        (
            np.ones_like(p),
            p,
            p**2,
            phi_centered,
            phi_centered**2,
            p * phi_centered,
        )
    )


def fit_provisional_correction_models(
    fit_frame: pd.DataFrame,
) -> dict[str, Any]:
    """
    Fit one provisional polynomial correction model per period and FD sector.

    Only successful theta/local-phi calibration cells are used. This model is
    deliberately modest and should be judged from the residual diagnostics
    before being promoted to production.
    """
    model_database: dict[str, Any] = {
        "schema_version": 1,
        "status": "provisional_diagnostic_model",
        "definition": "p_corr = p_meas * (1 + delta_p_over_p)",
        "model_terms": list(MODEL_TERMS),
        "phi_definition": "phi_centered_deg = phi_local_deg - 30",
        "periods": {},
    }

    for period in CALIBRATION_PERIODS:
        period_record: dict[str, Any] = {
            "metadata": asdict(period),
            "sectors": {},
        }

        for sector in range(1, 7):
            cells = fit_frame.loc[
                (fit_frame["period_key"] == period.key)
                & (fit_frame["sector"] == sector)
                & (fit_frame["cell_type"] == "theta_local_phi")
                & (fit_frame["success"])
                & np.isfinite(fit_frame["delta_p_over_p"])
                & np.isfinite(fit_frame["delta_p_over_p_error"])
                & (fit_frame["delta_p_over_p_error"] > 0.0)
            ].copy()

            sector_record: dict[str, Any] = {
                "status": "not_fitted",
                "n_cells": int(len(cells)),
                "coefficients": {},
                "covariance": [],
                "chi2": math.nan,
                "ndf": 0,
                "chi2_ndf": math.nan,
                "validity": {},
            }

            if len(cells) >= len(MODEL_TERMS) + 2:
                design = correction_design_matrix(
                    cells["mean_e_p_gev"].to_numpy(),
                    cells["mean_local_phi_deg"].to_numpy(),
                )
                observations = cells["delta_p_over_p"].to_numpy()
                uncertainties = cells["delta_p_over_p_error"].to_numpy()

                weights = 1.0 / uncertainties**2
                normal_matrix = design.T @ (weights[:, None] * design)
                right_hand_side = design.T @ (weights * observations)

                try:
                    covariance = np.linalg.inv(normal_matrix)
                    coefficients = covariance @ right_hand_side
                except np.linalg.LinAlgError:
                    covariance = np.linalg.pinv(normal_matrix)
                    coefficients = covariance @ right_hand_side
                # endtry

                predictions = design @ coefficients
                chi2 = float(
                    np.sum(((observations - predictions) / uncertainties) ** 2)
                )
                ndf = int(len(observations) - len(coefficients))

                sector_record = {
                    "status": "success",
                    "n_cells": int(len(cells)),
                    "coefficients": {
                        term: float(value)
                        for term, value in zip(MODEL_TERMS, coefficients)
                    },
                    "covariance": covariance.tolist(),
                    "chi2": chi2,
                    "ndf": ndf,
                    "chi2_ndf": chi2 / ndf if ndf > 0 else math.nan,
                    "validity": {
                        "p_gev": [
                            float(cells["mean_e_p_gev"].min()),
                            float(cells["mean_e_p_gev"].max()),
                        ],
                        "theta_deg": [
                            float(cells["mean_theta_deg"].min()),
                            float(cells["mean_theta_deg"].max()),
                        ],
                        "local_phi_deg": [
                            float(cells["mean_local_phi_deg"].min()),
                            float(cells["mean_local_phi_deg"].max()),
                        ],
                    },
                }
            # endif

            period_record["sectors"][str(sector)] = sector_record
        # endfor

        model_database["periods"][period.key] = period_record
    # endfor

    return model_database


def evaluate_provisional_model(
    model_record: dict[str, Any],
    momentum_gev: np.ndarray,
    local_phi_deg: np.ndarray,
) -> np.ndarray:
    """Evaluate one provisional sector correction model."""
    coefficients = np.asarray(
        [model_record["coefficients"][term] for term in MODEL_TERMS],
        dtype=float,
    )
    design = correction_design_matrix(momentum_gev, local_phi_deg)
    return design @ coefficients


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def save_peak_fit_plot(
    output_path: Path,
    diagnostics: dict[str, np.ndarray],
    fit_result: PeakFitResult,
    title: str,
) -> None:
    """Save one W-distribution fit diagnostic."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axis = plt.subplots(figsize=(8.0, 5.5))
    centers = diagnostics["centers"]
    counts = diagnostics["counts"]
    errors = np.sqrt(np.maximum(counts, 1.0))

    axis.errorbar(
        centers,
        counts,
        yerr=errors,
        fmt=".",
        markersize=3,
        linewidth=0.8,
        label="Data",
    )

    if fit_result.success:
        axis.plot(
            diagnostics["fit_x"],
            diagnostics["fit_y"],
            linewidth=1.5,
            label="Gaussian + quadratic fit",
        )
        annotation = (
            f"mu_W = {fit_result.peak_mean_gev:.6f} +/- "
            f"{fit_result.peak_mean_error_gev:.6f} GeV\n"
            f"sigma_W = {fit_result.peak_sigma_gev:.6f} GeV\n"
            f"chi2/ndf = {fit_result.chi2_ndf:.2f}\n"
            f"N = {fit_result.n_events}"
        )
    else:
        annotation = (
            f"Fit status: {fit_result.status}\n"
            f"N = {fit_result.n_events}"
        )
    # endif

    axis.axvline(PROTON_MASS_GEV, linestyle="--", linewidth=1.0, label="Proton mass")
    axis.set_xlabel("W [GeV]")
    axis.set_ylabel("Counts")
    axis.set_title(title)
    axis.text(
        0.03,
        0.97,
        annotation,
        transform=axis.transAxes,
        ha="left",
        va="top",
    )
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_period_integrated_w_plot(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    output_path: Path,
    histogram_range: tuple[float, float],
    histogram_bins: int,
) -> None:
    """Save a 2x3 sector canvas of integrated W spectra and fitted peak positions."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 9.0), sharex=True)
    axes_flat = axes.ravel()

    for sector in range(1, 7):
        axis = axes_flat[sector - 1]
        sector_values = frame.loc[
            frame["sector"] == sector,
            "w_recalculated",
        ].to_numpy()

        axis.hist(
            sector_values,
            bins=histogram_bins,
            range=histogram_range,
            histtype="step",
        )
        axis.axvline(PROTON_MASS_GEV, linestyle="--", linewidth=1.0)

        row = fit_frame.loc[
            (fit_frame["period_key"] == period.key)
            & (fit_frame["sector"] == sector)
            & (fit_frame["cell_type"] == "sector_integrated")
        ]

        if len(row) == 1 and bool(row.iloc[0]["success"]):
            peak = float(row.iloc[0]["peak_mean_gev"])
            axis.axvline(peak, linestyle="-", linewidth=1.0)
            axis.set_title(
                f"Sector {sector}: mu_W={peak:.5f} GeV"
            )
        else:
            axis.set_title(f"Sector {sector}: fit unavailable")
        # endif

        axis.set_xlabel("W [GeV]")
        axis.set_ylabel("Counts")
    # endfor

    figure.suptitle(f"{period.label}: inclusive eX elastic-peak diagnostics")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_peak_vs_theta_plot(
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the fitted W peak versus electron polar angle for all sectors."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 9.0), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    for sector in range(1, 7):
        axis = axes_flat[sector - 1]
        cells = fit_frame.loc[
            (fit_frame["period_key"] == period.key)
            & (fit_frame["sector"] == sector)
            & (fit_frame["cell_type"] == "theta")
            & (fit_frame["success"])
        ].sort_values("mean_theta_deg")

        axis.errorbar(
            cells["mean_theta_deg"],
            cells["peak_mean_gev"],
            yerr=cells["peak_mean_error_gev"],
            fmt="o",
            markersize=4,
            linewidth=1.0,
        )
        axis.axhline(PROTON_MASS_GEV, linestyle="--", linewidth=1.0)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Mean electron theta [deg]")
        axis.set_ylabel("Fitted W peak [GeV]")
    # endfor

    figure.suptitle(f"{period.label}: fitted elastic W peak versus electron theta")
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_correction_vs_theta_plot(
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot provisional delta-p/p values versus electron polar angle."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 9.0), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    for sector in range(1, 7):
        axis = axes_flat[sector - 1]
        cells = fit_frame.loc[
            (fit_frame["period_key"] == period.key)
            & (fit_frame["sector"] == sector)
            & (fit_frame["cell_type"] == "theta")
            & (fit_frame["success"])
            & np.isfinite(fit_frame["delta_p_over_p"])
        ].sort_values("mean_theta_deg")

        axis.errorbar(
            cells["mean_theta_deg"],
            100.0 * cells["delta_p_over_p"],
            yerr=100.0 * cells["delta_p_over_p_error"],
            fmt="o",
            markersize=4,
            linewidth=1.0,
        )
        axis.axhline(0.0, linestyle="--", linewidth=1.0)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Mean electron theta [deg]")
        axis.set_ylabel("Required momentum correction [%]")
    # endfor

    figure.suptitle(
        f"{period.label}: correction inferred from fitted elastic W peak"
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_correction_maps(
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    theta_edges_deg: np.ndarray,
    local_phi_edges_deg: np.ndarray,
    output_path: Path,
) -> None:
    """Save theta-versus-local-phi correction maps for all six sectors."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(2, 3, figsize=(16.0, 9.5), sharex=True, sharey=True)
    axes_flat = axes.ravel()

    valid_period_cells = fit_frame.loc[
        (fit_frame["period_key"] == period.key)
        & (fit_frame["cell_type"] == "theta_local_phi")
        & (fit_frame["success"])
        & np.isfinite(fit_frame["delta_p_over_p"])
    ]

    if valid_period_cells.empty:
        for axis in axes_flat:
            axis.text(0.5, 0.5, "No successful cells", ha="center", va="center")
        # endfor
        figure.suptitle(f"{period.label}: no correction maps available")
        figure.tight_layout()
        figure.savefig(output_path, dpi=180)
        plt.close(figure)
        return
    # endif

    correction_percent = 100.0 * valid_period_cells["delta_p_over_p"]
    global_limit = float(np.nanpercentile(np.abs(correction_percent), 98.0))
    global_limit = max(global_limit, 0.05)

    last_image = None
    for sector in range(1, 7):
        axis = axes_flat[sector - 1]
        grid = np.full(
            (len(theta_edges_deg) - 1, len(local_phi_edges_deg) - 1),
            np.nan,
            dtype=float,
        )

        cells = valid_period_cells.loc[
            valid_period_cells["sector"] == sector
        ]

        for _, row in cells.iterrows():
            theta_index = int(row["theta_bin_index"])
            phi_index = int(row["local_phi_bin_index"])
            grid[theta_index, phi_index] = 100.0 * float(row["delta_p_over_p"])
        # endfor

        last_image = axis.pcolormesh(
            local_phi_edges_deg,
            theta_edges_deg,
            grid,
            shading="auto",
            vmin=-global_limit,
            vmax=global_limit,
        )
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Local phi [deg]")
        axis.set_ylabel("Electron theta [deg]")
    # endfor

    if last_image is not None:
        colorbar = figure.colorbar(last_image, ax=axes_flat.tolist(), shrink=0.85)
        colorbar.set_label("Required momentum correction [%]")
    # endif

    figure.suptitle(
        f"{period.label}: correction inferred from fitted W peaks"
    )
    figure.subplots_adjust(
        left=0.07,
        right=0.88,
        bottom=0.08,
        top=0.92,
        wspace=0.22,
        hspace=0.25,
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_run_stability_plot(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    output_path: Path,
    run_bin_width: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
) -> None:
    """
    Fit the sector-integrated W peak in run-number bins and plot run stability.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    run_edges = np.arange(
        period.run_min,
        period.run_max + run_bin_width + 1,
        run_bin_width,
        dtype=int,
    )

    figure, axis = plt.subplots(figsize=(11.0, 6.0))

    for sector in range(1, 7):
        run_centers: list[float] = []
        peak_positions: list[float] = []
        peak_errors: list[float] = []

        sector_frame = frame.loc[frame["sector"] == sector]

        for run_index in range(len(run_edges) - 1):
            run_low = int(run_edges[run_index])
            run_high = int(run_edges[run_index + 1])
            run_frame = sector_frame.loc[
                (sector_frame["runnum"] >= run_low)
                & (sector_frame["runnum"] < run_high)
            ]

            fit_result, _ = fit_w_peak(
                run_frame["w_recalculated"].to_numpy(),
                histogram_range=histogram_range,
                fit_range=fit_range,
                histogram_bins=histogram_bins,
                minimum_events=minimum_events,
            )

            if fit_result.success:
                run_centers.append(float(run_frame["runnum"].mean()))
                peak_positions.append(fit_result.peak_mean_gev)
                peak_errors.append(fit_result.peak_mean_error_gev)
            # endif
        # endfor

        if run_centers:
            axis.errorbar(
                run_centers,
                peak_positions,
                yerr=peak_errors,
                fmt="o",
                markersize=3,
                linewidth=0.8,
                label=f"Sector {sector}",
            )
        # endif
    # endfor

    axis.axhline(PROTON_MASS_GEV, linestyle="--", linewidth=1.0)
    axis.set_xlabel("Run number")
    axis.set_ylabel("Fitted W peak [GeV]")
    axis.set_title(f"{period.label}: elastic-peak stability versus run")
    axis.legend(ncol=2)
    figure.tight_layout()
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_model_residual_plots(
    fit_frame: pd.DataFrame,
    model_database: dict[str, Any],
    output_directory: Path,
) -> None:
    """Plot calibration-cell residuals after applying each provisional model."""
    output_directory.mkdir(parents=True, exist_ok=True)

    for period in CALIBRATION_PERIODS:
        figure, axes = plt.subplots(
            2,
            3,
            figsize=(15.0, 9.0),
            sharex=True,
            sharey=True,
        )
        axes_flat = axes.ravel()

        for sector in range(1, 7):
            axis = axes_flat[sector - 1]
            model_record = model_database["periods"][period.key]["sectors"][
                str(sector)
            ]

            cells = fit_frame.loc[
                (fit_frame["period_key"] == period.key)
                & (fit_frame["sector"] == sector)
                & (fit_frame["cell_type"] == "theta_local_phi")
                & (fit_frame["success"])
                & np.isfinite(fit_frame["delta_p_over_p"])
            ].copy()

            if model_record["status"] != "success" or cells.empty:
                axis.text(0.5, 0.5, "Model unavailable", ha="center", va="center")
                axis.set_title(f"Sector {sector}")
                continue
            # endif

            prediction = evaluate_provisional_model(
                model_record,
                cells["mean_e_p_gev"].to_numpy(),
                cells["mean_local_phi_deg"].to_numpy(),
            )
            residual_percent = 100.0 * (
                cells["delta_p_over_p"].to_numpy() - prediction
            )

            axis.scatter(
                cells["mean_theta_deg"],
                residual_percent,
                s=18,
            )
            axis.axhline(0.0, linestyle="--", linewidth=1.0)
            axis.set_title(
                f"Sector {sector}, chi2/ndf={model_record['chi2_ndf']:.2f}"
            )
            axis.set_xlabel("Mean electron theta [deg]")
            axis.set_ylabel("Cell - model correction [%]")
        # endfor

        figure.suptitle(
            f"{period.label}: residuals of provisional correction model"
        )
        figure.tight_layout()
        figure.savefig(
            output_directory / f"{period.key}_model_residuals.png",
            dpi=180,
        )
        plt.close(figure)
    # endfor


# -----------------------------------------------------------------------------
# Serialization
# -----------------------------------------------------------------------------

def dataframe_to_cell_database(
    fit_frame: pd.DataFrame,
    arguments: argparse.Namespace,
) -> dict[str, Any]:
    """Convert the fit table into a nested JSON-compatible database."""
    database: dict[str, Any] = {
        "schema_version": 1,
        "status": "provisional_calibration_cells",
        "physics": {
            "particle": "electron",
            "beam_energy_definition": {
                "16042-17065": 10.5473,
                "17067-17716": 10.5563,
                "17717-17811": 10.5593
            },
            "reaction": "inclusive NH3 eX with fitted elastic W peak",
            "target_peak_gev": PROTON_MASS_GEV,
            "correction_definition": (
                "p_corr = p_meas * (1 + delta_p_over_p)"
            ),
            "direction_treatment": (
                "electron theta and phi are unchanged; only momentum magnitude "
                "is corrected"
            ),
        },
        "selection": {
            "fiducial_status": arguments.require_fiducial_status,
            "vertex_min_cm": arguments.vertex_min,
            "vertex_max_cm": arguments.vertex_max,
            "theta_min_deg": arguments.theta_min,
            "theta_max_deg": arguments.theta_max,
            "w_preselection_min_gev": arguments.w_preselection_min,
            "w_preselection_max_gev": arguments.w_preselection_max,
        },
        "fit_configuration": {
            "histogram_range_gev": [
                arguments.w_hist_min,
                arguments.w_hist_max,
            ],
            "fit_range_gev": [
                arguments.w_fit_min,
                arguments.w_fit_max,
            ],
            "histogram_bins": arguments.w_bins,
            "minimum_events_integrated": arguments.minimum_events_integrated,
            "minimum_events_cell": arguments.minimum_events_cell,
            "theta_edges_deg": arguments.theta_edges.tolist(),
            "local_phi_edges_deg": arguments.local_phi_edges.tolist(),
            "fit_model": "Gaussian elastic core + quadratic background",
        },
        "periods": {},
    }

    for period in CALIBRATION_PERIODS:
        period_entries = fit_frame.loc[
            fit_frame["period_key"] == period.key
        ].copy()

        database["periods"][period.key] = {
            "metadata": asdict(period),
            "cells": period_entries.replace({np.nan: None}).to_dict(
                orient="records"
            ),
        }
    # endfor

    return database


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write JSON while allowing NaN-like values to be represented as null."""
    path.parent.mkdir(parents=True, exist_ok=True)

    def clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {key: clean(item) for key, item in value.items()}
        # endif
        if isinstance(value, list):
            return [clean(item) for item in value]
        # endif
        if isinstance(value, tuple):
            return [clean(item) for item in value]
        # endif
        if isinstance(value, np.generic):
            return clean(value.item())
        # endif
        if isinstance(value, float) and not math.isfinite(value):
            return None
        # endif
        return value

    with path.open("w", encoding="utf-8") as output_file:
        json.dump(clean(payload), output_file, indent=2, sort_keys=False)
        output_file.write("\n")
    # endwith


# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------

def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Fit inclusive eX elastic W peaks and derive provisional RGC "
            "electron momentum corrections."
        )
    )

    parser.add_argument(
        "--su22-file",
        type=Path,
        default=DEFAULT_INPUTS["su22"],
        help=f"Su22 input ROOT file. Default: {DEFAULT_INPUTS['su22']}",
    )
    parser.add_argument(
        "--fa22-file",
        type=Path,
        default=DEFAULT_INPUTS["fa22"],
        help=f"Fa22 input ROOT file. Default: {DEFAULT_INPUTS['fa22']}",
    )
    parser.add_argument(
        "--sp23-file",
        type=Path,
        default=DEFAULT_INPUTS["sp23"],
        help=f"Sp23 input ROOT file. Default: {DEFAULT_INPUTS['sp23']}",
    )
    parser.add_argument(
        "--tree-name",
        default="PhysicsEvents",
        help="Input TTree name. Default: PhysicsEvents",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/electron_diagnostics_v1"),
        help="Output directory. Default: output/electron_diagnostics_v1",
    )
    parser.add_argument(
        "--step-size",
        default="250 MB",
        help="uproot iteration step size. Default: 250 MB",
    )

    parser.add_argument(
        "--require-fiducial-status",
        type=int,
        default=1,
        help=(
            "Required fiducial_status value. Use a negative value to disable "
            "this requirement. Default: 1"
        ),
    )
    parser.add_argument(
        "--vertex-min",
        type=float,
        default=-10.0,
        help="Minimum electron vertex z [cm]. Default: -10",
    )
    parser.add_argument(
        "--vertex-max",
        type=float,
        default=5.0,
        help="Maximum electron vertex z [cm]. Default: 5",
    )
    parser.add_argument(
        "--theta-min",
        type=float,
        default=5.0,
        help="Minimum electron theta [deg]. Default: 5",
    )
    parser.add_argument(
        "--theta-max",
        type=float,
        default=35.0,
        help="Maximum electron theta [deg]. Default: 35",
    )
    parser.add_argument(
        "--theta-bins",
        default="5,7,9,11,13,15,17,19,21,24,27,31,35",
        help=(
            "Comma-separated electron-theta bin edges [deg]. "
            "Default: 5,7,9,11,13,15,17,19,21,24,27,31,35"
        ),
    )
    parser.add_argument(
        "--local-phi-bins",
        default="0,10,20,30,40,50,60",
        help=(
            "Comma-separated FD local-phi bin edges [deg]. "
            "Default: 0,10,20,30,40,50,60"
        ),
    )

    parser.add_argument(
        "--w-preselection-min",
        type=float,
        default=0.65,
        help="Minimum input W retained [GeV]. Default: 0.65",
    )
    parser.add_argument(
        "--w-preselection-max",
        type=float,
        default=1.45,
        help="Maximum input W retained [GeV]. Default: 1.45",
    )
    parser.add_argument(
        "--w-hist-min",
        type=float,
        default=0.70,
        help="Minimum W histogram edge [GeV]. Default: 0.70",
    )
    parser.add_argument(
        "--w-hist-max",
        type=float,
        default=1.35,
        help="Maximum W histogram edge [GeV]. Default: 1.35",
    )
    parser.add_argument(
        "--w-fit-min",
        type=float,
        default=0.82,
        help="Minimum W fit edge [GeV]. Default: 0.82",
    )
    parser.add_argument(
        "--w-fit-max",
        type=float,
        default=1.12,
        help="Maximum W fit edge [GeV]. Default: 1.12",
    )
    parser.add_argument(
        "--w-bins",
        type=int,
        default=180,
        help="Number of W histogram bins. Default: 180",
    )
    parser.add_argument(
        "--minimum-events-integrated",
        type=int,
        default=2000,
        help="Minimum events for an integrated sector fit. Default: 2000",
    )
    parser.add_argument(
        "--minimum-events-cell",
        type=int,
        default=500,
        help="Minimum events for a theta or theta/local-phi fit. Default: 500",
    )
    parser.add_argument(
        "--run-bin-width",
        type=int,
        default=25,
        help="Run-number width used for stability plots. Default: 25",
    )
    parser.add_argument(
        "--save-individual-fits",
        action="store_true",
        help="Save every individual calibration-cell W fit.",
    )
    parser.add_argument(
        "--skip-run-stability",
        action="store_true",
        help="Skip the run-number stability plots.",
    )

    return parser


def validate_arguments(arguments: argparse.Namespace) -> None:
    """Validate command-line arguments and derive parsed edge arrays."""
    arguments.theta_edges = parse_edges(arguments.theta_bins, "--theta-bins")
    arguments.local_phi_edges = parse_edges(
        arguments.local_phi_bins,
        "--local-phi-bins",
    )

    if arguments.require_fiducial_status < 0:
        arguments.require_fiducial_status = None
    # endif

    if arguments.vertex_min >= arguments.vertex_max:
        raise ValueError("--vertex-min must be smaller than --vertex-max.")
    # endif

    if arguments.theta_min >= arguments.theta_max:
        raise ValueError("--theta-min must be smaller than --theta-max.")
    # endif

    if arguments.w_preselection_min >= arguments.w_preselection_max:
        raise ValueError(
            "--w-preselection-min must be smaller than --w-preselection-max."
        )
    # endif

    if arguments.w_hist_min >= arguments.w_hist_max:
        raise ValueError("--w-hist-min must be smaller than --w-hist-max.")
    # endif

    if arguments.w_fit_min >= arguments.w_fit_max:
        raise ValueError("--w-fit-min must be smaller than --w-fit-max.")
    # endif

    if (
        arguments.w_fit_min < arguments.w_hist_min
        or arguments.w_fit_max > arguments.w_hist_max
    ):
        raise ValueError("The W fit range must lie inside the histogram range.")
    # endif

    if arguments.w_bins < 20:
        raise ValueError("--w-bins must be at least 20.")
    # endif

    if arguments.minimum_events_integrated <= 0:
        raise ValueError("--minimum-events-integrated must be positive.")
    # endif

    if arguments.minimum_events_cell <= 0:
        raise ValueError("--minimum-events-cell must be positive.")
    # endif

    if arguments.run_bin_width <= 0:
        raise ValueError("--run-bin-width must be positive.")
    # endif


def main() -> int:
    """Run the complete diagnostic and provisional calibration workflow."""
    parser = build_argument_parser()
    arguments = parser.parse_args()

    try:
        validate_arguments(arguments)
    except (ValueError, argparse.ArgumentTypeError) as exc:
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
                f"FATAL: input file for {source_key} does not exist: {file_path}",
                file=sys.stderr,
            )
            return 2
        # endif
    # endfor

    output_directory = arguments.output_dir
    plot_directory = output_directory / "plots"
    individual_fit_directory = plot_directory / "individual_fits"
    output_directory.mkdir(parents=True, exist_ok=True)
    plot_directory.mkdir(parents=True, exist_ok=True)

    loaded_frames: dict[str, pd.DataFrame] = {}
    all_fit_records: list[dict[str, Any]] = []

    histogram_range = (arguments.w_hist_min, arguments.w_hist_max)
    fit_range = (arguments.w_fit_min, arguments.w_fit_max)

    for period in CALIBRATION_PERIODS:
        file_path = input_files[period.source_key]
        tree_name = find_tree(file_path, arguments.tree_name)
        branches = validate_branches(file_path, tree_name)

        print(
            f"[read] {period.label}: {file_path}:{tree_name}, "
            f"runs {period.run_min}-{period.run_max}"
        )

        period_frame = read_selected_events(
            file_path=file_path,
            tree_name=tree_name,
            branches=branches,
            period=period,
            step_size=arguments.step_size,
            require_fiducial_status=arguments.require_fiducial_status,
            vertex_min_cm=arguments.vertex_min,
            vertex_max_cm=arguments.vertex_max,
            theta_min_deg=arguments.theta_min,
            theta_max_deg=arguments.theta_max,
            w_preselection_min_gev=arguments.w_preselection_min,
            w_preselection_max_gev=arguments.w_preselection_max,
        )

        if period_frame.empty:
            print(
                f"WARNING: no selected events for {period.label}.",
                file=sys.stderr,
            )
            loaded_frames[period.key] = period_frame
            continue
        # endif

        loaded_frames[period.key] = period_frame
        print(
            f"[read] {period.label}: retained {len(period_frame):,} events; "
            f"runs {period_frame['runnum'].min()}-"
            f"{period_frame['runnum'].max()}"
        )

        period_records = fit_calibration_cells(
            frame=period_frame,
            period=period,
            theta_edges_deg=arguments.theta_edges,
            local_phi_edges_deg=arguments.local_phi_edges,
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=arguments.w_bins,
            minimum_events_integrated=arguments.minimum_events_integrated,
            minimum_events_cell=arguments.minimum_events_cell,
            save_individual_fits=arguments.save_individual_fits,
            individual_fit_directory=individual_fit_directory,
        )
        all_fit_records.extend(period_records)
    # endfor

    if not all_fit_records:
        print("FATAL: no calibration fits were attempted.", file=sys.stderr)
        return 3
    # endif

    fit_frame = pd.DataFrame(all_fit_records)
    fit_csv_path = output_directory / "fit_results_v1.csv"
    fit_frame.to_csv(fit_csv_path, index=False, float_format="%.12g")
    print(f"[write] {fit_csv_path}")

    cell_database = dataframe_to_cell_database(fit_frame, arguments)
    cell_json_path = output_directory / "electron_correction_cells_v1.json"
    write_json(cell_json_path, cell_database)
    print(f"[write] {cell_json_path}")

    model_database = fit_provisional_correction_models(fit_frame)
    model_json_path = output_directory / "electron_correction_models_v1.json"
    write_json(model_json_path, model_database)
    print(f"[write] {model_json_path}")

    for period in CALIBRATION_PERIODS:
        period_frame = loaded_frames.get(period.key, pd.DataFrame())
        if period_frame.empty:
            continue
        # endif

        save_period_integrated_w_plot(
            frame=period_frame,
            period=period,
            fit_frame=fit_frame,
            output_path=plot_directory
            / f"{period.key}_integrated_w_by_sector.png",
            histogram_range=histogram_range,
            histogram_bins=arguments.w_bins,
        )
        save_peak_vs_theta_plot(
            period=period,
            fit_frame=fit_frame,
            output_path=plot_directory
            / f"{period.key}_w_peak_vs_theta.png",
        )
        save_correction_vs_theta_plot(
            period=period,
            fit_frame=fit_frame,
            output_path=plot_directory
            / f"{period.key}_correction_vs_theta.png",
        )
        save_correction_maps(
            period=period,
            fit_frame=fit_frame,
            theta_edges_deg=arguments.theta_edges,
            local_phi_edges_deg=arguments.local_phi_edges,
            output_path=plot_directory
            / f"{period.key}_correction_theta_local_phi_maps.png",
        )

        if not arguments.skip_run_stability:
            save_run_stability_plot(
                frame=period_frame,
                period=period,
                output_path=plot_directory
                / f"{period.key}_w_peak_vs_run.png",
                run_bin_width=arguments.run_bin_width,
                histogram_range=histogram_range,
                fit_range=fit_range,
                histogram_bins=arguments.w_bins,
                minimum_events=arguments.minimum_events_cell,
            )
        # endif
    # endfor

    save_model_residual_plots(
        fit_frame=fit_frame,
        model_database=model_database,
        output_directory=plot_directory,
    )

    successful_cells = fit_frame.loc[
        (fit_frame["cell_type"] == "theta_local_phi")
        & (fit_frame["success"])
    ]
    attempted_cells = fit_frame.loc[
        fit_frame["cell_type"] == "theta_local_phi"
    ]

    print("")
    print("Electron momentum-correction diagnostic complete.")
    print(
        f"Successful theta/local-phi fits: {len(successful_cells):,} / "
        f"{len(attempted_cells):,}"
    )
    print(f"Output directory: {output_directory}")
    print("")
    print(
        "The polynomial JSON is provisional. Review the fitted W spectra, "
        "correction maps, run stability, and model residuals before using it "
        "in production."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
# endif
