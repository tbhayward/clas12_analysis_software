#!/usr/bin/env python3
"""
derive_electron_momentum_corrections_v14.py

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
    output/electron_diagnostics/csv/fit_results_uncorrected.csv
    output/electron_diagnostics/csv/fit_results_corrected.csv
    output/electron_diagnostics/csv/correction_model_comparison.csv
    output/electron_diagnostics/json/electron_correction_cells.json
    output/electron_diagnostics/json/electron_correction_models.json
    output/electron_diagnostics/json/closure_summary.json
    output/electron_diagnostics/plots/<plot_category>/*.png

The workflow derives a smooth theta-dependent correction independently for
each calibration period and FD sector, applies it event by event, recomputes W,
and performs a closure test using exactly the same theta-bin boundaries.

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
import json
import math
import shutil
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
    "runnum",
    "evnum",
    "e_p",
    "e_theta",
    "e_phi",
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
    quality_pass: bool
    quality_flags: str
    fit_mode: str
    sigma_fixed: bool
    mean_constraint_center_gev: float
    mean_constraint_half_width_gev: float
    n_events: int
    n_histogram_entries: int
    peak_mean_gev: float
    peak_mean_error_gev: float
    peak_sigma_gev: float
    peak_sigma_error_gev: float
    amplitude: float
    amplitude_error: float
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
    fit_mode: str = "unknown",
    sigma_fixed: bool = False,
    mean_constraint_center_gev: float = math.nan,
    mean_constraint_half_width_gev: float = math.nan,
) -> PeakFitResult:
    """Construct a failed fit record with NaN numerical fields."""
    return PeakFitResult(
        success=False,
        status=status,
        quality_pass=False,
        quality_flags=status,
        fit_mode=fit_mode,
        sigma_fixed=sigma_fixed,
        mean_constraint_center_gev=mean_constraint_center_gev,
        mean_constraint_half_width_gev=mean_constraint_half_width_gev,
        n_events=n_events,
        n_histogram_entries=n_histogram_entries,
        peak_mean_gev=math.nan,
        peak_mean_error_gev=math.nan,
        peak_sigma_gev=math.nan,
        peak_sigma_error_gev=math.nan,
        amplitude=math.nan,
        amplitude_error=math.nan,
        chi2=math.nan,
        ndf=0,
        chi2_ndf=math.nan,
        fit_range_low_gev=fit_range[0],
        fit_range_high_gev=fit_range[1],
    )


def gaussian_fixed_sigma_plus_quadratic(
    x: np.ndarray,
    amplitude: float,
    mean: float,
    b0: float,
    b1: float,
    b2: float,
    fixed_sigma_gev: float,
) -> np.ndarray:
    """Gaussian with fixed width plus a quadratic continuum."""
    gaussian = amplitude * np.exp(
        -0.5 * ((x - mean) / fixed_sigma_gev) ** 2
    )
    x_centered = x - PROTON_MASS_GEV
    background = b0 + b1 * x_centered + b2 * x_centered**2
    return gaussian + background


def fit_w_peak(
    w_values: np.ndarray,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
    fixed_sigma_gev: float | None = None,
    mean_center_gev: float | None = None,
    mean_half_window_gev: float = 0.040,
) -> tuple[PeakFitResult, dict[str, np.ndarray]]:
    """
    Fit the elastic structure in an inclusive W spectrum.

    Sector-integrated spectra are fit with a floating Gaussian width. The four
    sector-specific theta cells are then fit with sigma_W fixed to the
    sector-integrated value and mu_W restricted to the integrated mean plus or
    minus mean_half_window_gev. Automatic checks remain advisory.
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
        "fit_mask": np.zeros_like(histogram_centers, dtype=bool),
        "fit_centers": np.asarray([], dtype=float),
        "fit_counts": np.asarray([], dtype=float),
        "fit_errors": np.asarray([], dtype=float),
        "fit_model_at_centers": np.asarray([], dtype=float),
        "pulls": np.asarray([], dtype=float),
        "fit_x": np.asarray([], dtype=float),
        "fit_y": np.asarray([], dtype=float),
        "fit_signal_y": np.asarray([], dtype=float),
        "fit_background_y": np.asarray([], dtype=float),
    }

    fit_mode = (
        "cell_fixed_sigma"
        if fixed_sigma_gev is not None
        else "sector_integrated_free_sigma"
    )
    sigma_fixed = fixed_sigma_gev is not None
    constraint_center = (
        float(mean_center_gev)
        if mean_center_gev is not None
        else PROTON_MASS_GEV
    )

    if n_events < minimum_events:
        return (
            failed_peak_fit(
                status="skipped_insufficient_events",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
                fit_mode=fit_mode,
                sigma_fixed=sigma_fixed,
                mean_constraint_center_gev=constraint_center,
                mean_constraint_half_width_gev=mean_half_window_gev,
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
    sigma_y = np.sqrt(np.maximum(y_fit, 1.0))

    diagnostics["fit_mask"] = fit_mask
    diagnostics["fit_centers"] = x_fit
    diagnostics["fit_counts"] = y_fit
    diagnostics["fit_errors"] = sigma_y

    minimum_fit_points = 12 if sigma_fixed else 14
    if x_fit.size < minimum_fit_points or np.count_nonzero(y_fit) < 10:
        return (
            failed_peak_fit(
                status="skipped_insufficient_populated_bins",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
                fit_mode=fit_mode,
                sigma_fixed=sigma_fixed,
                mean_constraint_center_gev=constraint_center,
                mean_constraint_half_width_gev=mean_half_window_gev,
            ),
            diagnostics,
        )
    # endif

    if fixed_sigma_gev is not None:
        if not np.isfinite(fixed_sigma_gev) or fixed_sigma_gev <= 0.0:
            return (
                failed_peak_fit(
                    status="invalid_fixed_sigma",
                    n_events=n_events,
                    n_histogram_entries=int(histogram_counts.sum()),
                    fit_range=fit_range,
                    fit_mode=fit_mode,
                    sigma_fixed=True,
                    mean_constraint_center_gev=constraint_center,
                    mean_constraint_half_width_gev=mean_half_window_gev,
                ),
                diagnostics,
            )
        # endif
        mean_lower = max(
            fit_range[0],
            constraint_center - mean_half_window_gev,
        )
        mean_upper = min(
            fit_range[1],
            constraint_center + mean_half_window_gev,
        )
    else:
        mean_lower = max(fit_range[0], PROTON_MASS_GEV - 0.080)
        mean_upper = min(fit_range[1], PROTON_MASS_GEV + 0.080)
    # endif

    peak_search_mask = (x_fit >= mean_lower) & (x_fit <= mean_upper)
    sideband_exclusion = max(
        0.045,
        2.0 * fixed_sigma_gev if fixed_sigma_gev is not None else 0.055,
    )
    sideband_mask = (
        (x_fit < constraint_center - sideband_exclusion)
        | (x_fit > constraint_center + sideband_exclusion)
    )

    background_coefficients: np.ndarray | None = None
    if np.count_nonzero(sideband_mask) >= 8:
        try:
            background_coefficients = np.polyfit(
                x_fit[sideband_mask] - PROTON_MASS_GEV,
                y_fit[sideband_mask],
                deg=2,
                w=1.0 / sigma_y[sideband_mask],
            )
            background_seed = np.polyval(
                background_coefficients,
                x_fit - PROTON_MASS_GEV,
            )
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            background_seed = np.full_like(
                y_fit,
                float(np.median(y_fit[sideband_mask])),
            )
        # endtry
    else:
        background_seed = np.full_like(
            y_fit,
            float(np.percentile(y_fit, 25.0)),
        )
    # endif

    residual_seed = y_fit - background_seed
    search_indices = np.flatnonzero(peak_search_mask)
    if search_indices.size == 0:
        return (
            failed_peak_fit(
                status="empty_peak_search_window",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
                fit_mode=fit_mode,
                sigma_fixed=sigma_fixed,
                mean_constraint_center_gev=constraint_center,
                mean_constraint_half_width_gev=mean_half_window_gev,
            ),
            diagnostics,
        )
    # endif

    peak_index = search_indices[np.argmax(residual_seed[peak_search_mask])]
    initial_mean = float(np.clip(x_fit[peak_index], mean_lower, mean_upper))
    initial_amplitude = max(float(residual_seed[peak_index]), 1.0)

    if background_coefficients is not None:
        initial_b2 = float(background_coefficients[0])
        initial_b1 = float(background_coefficients[1])
        initial_b0 = float(background_coefficients[2])
    else:
        initial_b0 = max(float(np.median(background_seed)), 0.0)
        initial_b1 = 0.0
        initial_b2 = 0.0
    # endif

    try:
        if fixed_sigma_gev is None:
            initial_sigma = 0.025
            initial_parameters = (
                initial_amplitude,
                initial_mean,
                initial_sigma,
                initial_b0,
                initial_b1,
                initial_b2,
            )
            lower_bounds = (
                0.0,
                mean_lower,
                0.010,
                -np.inf,
                -np.inf,
                -np.inf,
            )
            upper_bounds = (
                np.inf,
                mean_upper,
                0.060,
                np.inf,
                np.inf,
                np.inf,
            )
            parameters, covariance = curve_fit(
                gaussian_plus_quadratic,
                x_fit,
                y_fit,
                p0=initial_parameters,
                sigma=sigma_y,
                absolute_sigma=True,
                bounds=(lower_bounds, upper_bounds),
                maxfev=100_000,
            )
            amplitude, mean, fitted_sigma, b0, b1, b2 = parameters
            parameter_errors = np.sqrt(
                np.clip(np.diag(covariance), 0.0, None)
            )
            amplitude_error = float(parameter_errors[0])
            mean_error = float(parameter_errors[1])
            sigma_error = float(parameter_errors[2])
            model_at_centers = gaussian_plus_quadratic(x_fit, *parameters)
            n_parameters = 6
            sigma_lower = lower_bounds[2]
            sigma_upper = upper_bounds[2]
        else:
            def fixed_sigma_model(
                x: np.ndarray,
                amplitude: float,
                mean: float,
                b0: float,
                b1: float,
                b2: float,
            ) -> np.ndarray:
                return gaussian_fixed_sigma_plus_quadratic(
                    x,
                    amplitude,
                    mean,
                    b0,
                    b1,
                    b2,
                    fixed_sigma_gev,
                )

            initial_parameters = (
                initial_amplitude,
                initial_mean,
                initial_b0,
                initial_b1,
                initial_b2,
            )
            lower_bounds = (
                0.0,
                mean_lower,
                -np.inf,
                -np.inf,
                -np.inf,
            )
            upper_bounds = (
                np.inf,
                mean_upper,
                np.inf,
                np.inf,
                np.inf,
            )
            parameters, covariance = curve_fit(
                fixed_sigma_model,
                x_fit,
                y_fit,
                p0=initial_parameters,
                sigma=sigma_y,
                absolute_sigma=True,
                bounds=(lower_bounds, upper_bounds),
                maxfev=100_000,
            )
            amplitude, mean, b0, b1, b2 = parameters
            fitted_sigma = float(fixed_sigma_gev)
            parameter_errors = np.sqrt(
                np.clip(np.diag(covariance), 0.0, None)
            )
            amplitude_error = float(parameter_errors[0])
            mean_error = float(parameter_errors[1])
            sigma_error = 0.0
            model_at_centers = fixed_sigma_model(x_fit, *parameters)
            n_parameters = 5
            sigma_lower = math.nan
            sigma_upper = math.nan
        # endif
    except (RuntimeError, ValueError, FloatingPointError) as exc:
        return (
            failed_peak_fit(
                status=f"fit_failed:{type(exc).__name__}",
                n_events=n_events,
                n_histogram_entries=int(histogram_counts.sum()),
                fit_range=fit_range,
                fit_mode=fit_mode,
                sigma_fixed=sigma_fixed,
                mean_constraint_center_gev=constraint_center,
                mean_constraint_half_width_gev=mean_half_window_gev,
            ),
            diagnostics,
        )
    # endtry

    chi2 = float(np.sum(((y_fit - model_at_centers) / sigma_y) ** 2))
    ndf = int(x_fit.size - n_parameters)
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan
    pulls = (y_fit - model_at_centers) / sigma_y

    dense_x = np.linspace(fit_range[0], fit_range[1], 800)
    centered_dense_x = dense_x - PROTON_MASS_GEV
    signal_y = amplitude * np.exp(
        -0.5 * ((dense_x - mean) / fitted_sigma) ** 2
    )
    background_y = b0 + b1 * centered_dense_x + b2 * centered_dense_x**2

    diagnostics["fit_model_at_centers"] = model_at_centers
    diagnostics["pulls"] = pulls
    diagnostics["fit_x"] = dense_x
    diagnostics["fit_y"] = signal_y + background_y
    diagnostics["fit_signal_y"] = signal_y
    diagnostics["fit_background_y"] = background_y

    quality_problems: list[str] = []
    mean_margin = min(mean - mean_lower, mean_upper - mean)
    if mean_margin < 0.002:
        quality_problems.append("mean_at_bound")
    # endif

    if fixed_sigma_gev is None:
        sigma_margin = min(
            fitted_sigma - sigma_lower,
            sigma_upper - fitted_sigma,
        )
        if sigma_margin < 0.002:
            quality_problems.append("sigma_at_bound")
        # endif
    # endif

    if not np.isfinite(mean_error) or mean_error > 0.020:
        quality_problems.append("poor_mean_precision")
    # endif
    if not np.isfinite(chi2_ndf) or chi2_ndf > 5.0:
        quality_problems.append("poor_chi2")
    # endif
    if (
        not np.isfinite(amplitude_error)
        or amplitude_error <= 0.0
        or amplitude / amplitude_error < 3.0
    ):
        quality_problems.append("weak_peak")
    # endif

    quality_pass = not quality_problems
    quality_flags = ",".join(quality_problems)
    status = "success" if quality_pass else "success_flagged_for_review"

    result = PeakFitResult(
        success=True,
        status=status,
        quality_pass=quality_pass,
        quality_flags=quality_flags,
        fit_mode=fit_mode,
        sigma_fixed=sigma_fixed,
        mean_constraint_center_gev=constraint_center,
        mean_constraint_half_width_gev=mean_half_window_gev,
        n_events=n_events,
        n_histogram_entries=int(histogram_counts.sum()),
        peak_mean_gev=float(mean),
        peak_mean_error_gev=mean_error,
        peak_sigma_gev=float(abs(fitted_sigma)),
        peak_sigma_error_gev=sigma_error,
        amplitude=float(amplitude),
        amplitude_error=amplitude_error,
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
        w_branch = arrays["W"].astype(float, copy=False)

        mask = (
            (runnum >= period.run_min)
            & (runnum <= period.run_max)
            & np.isfinite(e_p)
            & np.isfinite(e_theta)
            & np.isfinite(e_phi)
            & np.isfinite(w_branch)
            & (e_p > 0.0)
            & (e_theta * RAD2DEG >= theta_min_deg)
            & (e_theta * RAD2DEG < theta_max_deg)
            & (w_branch >= w_preselection_min_gev)
            & (w_branch <= w_preselection_max_gev)
        )

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
            "e_p": e_p[selected_indices],
            "e_theta_rad": e_theta[selected_indices],
            "e_theta_deg": e_theta[selected_indices] * RAD2DEG,
            "e_phi_rad": e_phi[selected_indices],
            "e_phi_deg": wrap_phi_deg(e_phi[selected_indices] * RAD2DEG),
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


def build_theta_edges_with_widest_bin_split(
    theta_values_deg: np.ndarray,
    base_bin_count: int,
) -> tuple[np.ndarray, int]:
    """
    Build equal-population theta edges and split the angularly widest base bin.

    First, base_bin_count quantile bins are constructed independently for one
    FD sector. The base interval with the largest angular width is then divided
    at the median theta value of the events inside that interval. The final
    number of bins is therefore base_bin_count + 1.

    Returns
    -------
    final_edges_deg
        Strictly increasing theta-bin edges.
    split_base_bin_index
        Zero-based index of the base bin that was subdivided.
    """
    theta_values = np.asarray(theta_values_deg, dtype=float)
    theta_values = theta_values[np.isfinite(theta_values)]

    if theta_values.size == 0:
        return np.linspace(0.0, 1.0, base_bin_count + 2), -1
    # endif

    quantiles = np.linspace(0.0, 1.0, base_bin_count + 1)
    base_edges = np.quantile(theta_values, quantiles).astype(float)

    # Guard against tied quantiles. The correction is far below detector
    # angular resolution and only ensures valid interval construction.
    for edge_index in range(1, len(base_edges)):
        if base_edges[edge_index] <= base_edges[edge_index - 1]:
            base_edges[edge_index] = base_edges[edge_index - 1] + 1.0e-6
        # endif
    # endfor

    base_widths = np.diff(base_edges)
    split_base_bin_index = int(np.argmax(base_widths))

    split_low = float(base_edges[split_base_bin_index])
    split_high = float(base_edges[split_base_bin_index + 1])

    if split_base_bin_index == base_bin_count - 1:
        in_split_bin = (
            (theta_values >= split_low)
            & (theta_values <= split_high)
        )
    else:
        in_split_bin = (
            (theta_values >= split_low)
            & (theta_values < split_high)
        )
    # endif

    split_values = theta_values[in_split_bin]
    if split_values.size >= 2:
        split_edge = float(np.median(split_values))
    else:
        split_edge = 0.5 * (split_low + split_high)
    # endif

    # Ensure that the inserted edge remains strictly interior.
    epsilon = max(1.0e-6, 1.0e-6 * abs(split_high - split_low))
    split_edge = float(
        np.clip(
            split_edge,
            split_low + epsilon,
            split_high - epsilon,
        )
    )

    final_edges = np.insert(
        base_edges,
        split_base_bin_index + 1,
        split_edge,
    )

    for edge_index in range(1, len(final_edges)):
        if final_edges[edge_index] <= final_edges[edge_index - 1]:
            final_edges[edge_index] = final_edges[edge_index - 1] + 1.0e-6
        # endif
    # endfor

    return final_edges, split_base_bin_index


def fit_calibration_cells(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    theta_bin_count: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    save_individual_fits: bool,
    individual_fit_directory: Path,
) -> list[dict[str, Any]]:
    """
    Fit one sector-integrated W spectrum and theta-dependent calibration cells.

    For each FD sector:

      1. Build theta_bin_count equal-population base bins.
      2. Identify the base bin with the largest angular width.
      3. Split that widest base bin at its event median.
      4. Fit the resulting theta_bin_count + 1 cells.

    The sector-integrated fit determines sigma_W. The same sigma_W is fixed in
    every theta-cell fit, while mu_W is constrained to lie within 40 MeV of the
    integrated-sector value. All theta cells are integrated over local phi.
    """
    records: list[dict[str, Any]] = []

    for sector in range(1, 7):
        sector_frame = frame.loc[frame["sector"] == sector].copy()
        sector_frame = sector_frame.sort_values("e_theta_deg")

        theta_edges_deg, split_base_bin_index = (
            build_theta_edges_with_widest_bin_split(
                sector_frame["e_theta_deg"].to_numpy()
                if not sector_frame.empty
                else np.asarray([], dtype=float),
                theta_bin_count,
            )
        )
        final_theta_bin_count = len(theta_edges_deg) - 1

        integrated_fit, integrated_diagnostics = fit_w_peak(
            sector_frame["w_recalculated"].to_numpy(),
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events_integrated,
        )
        integrated_record = make_fit_record(
            period=period,
            sector=sector,
            cell_type="sector_integrated",
            theta_bin_index=None,
            theta_low_deg=float(theta_edges_deg[0]),
            theta_high_deg=float(theta_edges_deg[-1]),
            local_phi_bin_index=None,
            local_phi_low_deg=0.0,
            local_phi_high_deg=60.0,
            cell_frame=sector_frame,
            fit_result=integrated_fit,
        )
        integrated_record["theta_base_bin_count"] = theta_bin_count
        integrated_record["theta_final_bin_count"] = final_theta_bin_count
        integrated_record["split_base_bin_index"] = split_base_bin_index
        records.append(integrated_record)

        if save_individual_fits:
            save_peak_fit_plot(
                output_path=individual_fit_directory
                / f"{period.key}_sector{sector}_integrated.png",
                diagnostics=integrated_diagnostics,
                fit_result=integrated_fit,
                title=f"{period.label}, sector {sector}, integrated",
            )
        # endif

        if not integrated_fit.success:
            for theta_index in range(final_theta_bin_count):
                theta_low = float(theta_edges_deg[theta_index])
                theta_high = float(theta_edges_deg[theta_index + 1])
                failed_fit = failed_peak_fit(
                    status="skipped_integrated_fit_failed",
                    n_events=0,
                    fit_range=fit_range,
                    fit_mode="cell_fixed_sigma",
                    sigma_fixed=True,
                )
                failed_record = make_fit_record(
                    period=period,
                    sector=sector,
                    cell_type="theta",
                    theta_bin_index=theta_index,
                    theta_low_deg=theta_low,
                    theta_high_deg=theta_high,
                    local_phi_bin_index=None,
                    local_phi_low_deg=0.0,
                    local_phi_high_deg=60.0,
                    cell_frame=pd.DataFrame(),
                    fit_result=failed_fit,
                )
                failed_record["theta_base_bin_count"] = theta_bin_count
                failed_record["theta_final_bin_count"] = final_theta_bin_count
                failed_record["split_base_bin_index"] = split_base_bin_index
                failed_record["is_split_sub_bin"] = (
                    theta_index == split_base_bin_index
                    or theta_index == split_base_bin_index + 1
                )
                records.append(failed_record)
            # endfor
            continue
        # endif

        for theta_index in range(final_theta_bin_count):
            theta_low = float(theta_edges_deg[theta_index])
            theta_high = float(theta_edges_deg[theta_index + 1])

            if theta_index == final_theta_bin_count - 1:
                theta_mask = (
                    (sector_frame["e_theta_deg"] >= theta_low)
                    & (sector_frame["e_theta_deg"] <= theta_high)
                )
            else:
                theta_mask = (
                    (sector_frame["e_theta_deg"] >= theta_low)
                    & (sector_frame["e_theta_deg"] < theta_high)
                )
            # endif

            theta_frame = sector_frame.loc[theta_mask].copy()
            theta_fit, theta_diagnostics = fit_w_peak(
                theta_frame["w_recalculated"].to_numpy(),
                histogram_range=histogram_range,
                fit_range=fit_range,
                histogram_bins=histogram_bins,
                minimum_events=minimum_events_cell,
                fixed_sigma_gev=integrated_fit.peak_sigma_gev,
                mean_center_gev=integrated_fit.peak_mean_gev,
                mean_half_window_gev=0.040,
            )

            theta_record = make_fit_record(
                period=period,
                sector=sector,
                cell_type="theta",
                theta_bin_index=theta_index,
                theta_low_deg=theta_low,
                theta_high_deg=theta_high,
                local_phi_bin_index=None,
                local_phi_low_deg=0.0,
                local_phi_high_deg=60.0,
                cell_frame=theta_frame,
                fit_result=theta_fit,
            )
            theta_record["theta_base_bin_count"] = theta_bin_count
            theta_record["theta_final_bin_count"] = final_theta_bin_count
            theta_record["split_base_bin_index"] = split_base_bin_index
            theta_record["is_split_sub_bin"] = (
                theta_index == split_base_bin_index
                or theta_index == split_base_bin_index + 1
            )
            records.append(theta_record)

            if save_individual_fits:
                split_note = (
                    ", split widest-bin sub-cell"
                    if theta_record["is_split_sub_bin"]
                    else ""
                )
                save_peak_fit_plot(
                    output_path=individual_fit_directory
                    / (
                        f"{period.key}_sector{sector}_theta"
                        f"{theta_index:02d}.png"
                    ),
                    diagnostics=theta_diagnostics,
                    fit_result=theta_fit,
                    title=(
                        f"{period.label}, sector {sector}, theta cell "
                        f"{theta_index + 1}/{final_theta_bin_count}"
                        f"{split_note}: "
                        f"{theta_low:.3f} <= theta_e < {theta_high:.3f} deg"
                    ),
                )
            # endif
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
# Theta-model selection, application, and closure
# -----------------------------------------------------------------------------

def weighted_polynomial_fit(
    theta_deg: np.ndarray,
    correction: np.ndarray,
    correction_error: np.ndarray,
    degree: int,
    error_floor: float = 2.0e-4,
) -> dict[str, Any]:
    """
    Fit delta-p/p as a polynomial in theta.

    An uncertainty floor of 2e-4 (0.02 percent absolute in delta-p/p) prevents
    extremely small formal peak-fit errors from giving individual calibration
    cells unrealistically large leverage.
    """
    x = np.asarray(theta_deg, dtype=float)
    y = np.asarray(correction, dtype=float)
    yerr = np.asarray(correction_error, dtype=float)

    valid = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(valid) < degree + 2:
        return {
            "success": False,
            "degree": degree,
            "status": "insufficient_points",
            "coefficients_ascending": [],
            "covariance": [],
            "chi2": math.nan,
            "ndf": 0,
            "chi2_ndf": math.nan,
            "aicc": math.nan,
            "n_points": int(np.count_nonzero(valid)),
        }
    # endif

    x = x[valid]
    y = y[valid]
    yerr = yerr[valid]
    effective_error = np.where(
        np.isfinite(yerr) & (yerr > 0.0),
        np.maximum(yerr, error_floor),
        error_floor,
    )
    weights = 1.0 / effective_error

    try:
        coefficients_descending, covariance_descending = np.polyfit(
            x,
            y,
            deg=degree,
            w=weights,
            cov=True,
        )
    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
        return {
            "success": False,
            "degree": degree,
            "status": f"fit_failed:{type(exc).__name__}",
            "coefficients_ascending": [],
            "covariance": [],
            "chi2": math.nan,
            "ndf": 0,
            "chi2_ndf": math.nan,
            "aicc": math.nan,
            "n_points": int(x.size),
        }
    # endtry

    prediction = np.polyval(coefficients_descending, x)
    chi2 = float(np.sum(((y - prediction) / effective_error) ** 2))
    parameter_count = degree + 1
    ndf = int(x.size - parameter_count)
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan

    if x.size > parameter_count + 1:
        aicc = (
            chi2
            + 2.0 * parameter_count
            + (
                2.0
                * parameter_count
                * (parameter_count + 1)
                / (x.size - parameter_count - 1)
            )
        )
    else:
        aicc = math.inf
    # endif

    coefficients_ascending = coefficients_descending[::-1]
    covariance_ascending = covariance_descending[::-1, ::-1]

    return {
        "success": True,
        "degree": degree,
        "status": "success",
        "coefficients_ascending": coefficients_ascending.tolist(),
        "covariance": covariance_ascending.tolist(),
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
        "aicc": float(aicc),
        "n_points": int(x.size),
        "effective_error_floor": error_floor,
    }


def select_theta_correction_models(
    period: CalibrationPeriod,
    fit_records: list[dict[str, Any]],
    linear_chi2_ndf_max: float = 2.5,
    quadratic_aicc_improvement_min: float = 2.0,
) -> list[dict[str, Any]]:
    """
    Fit linear and quadratic correction models and select the lowest adequate
    order independently for every FD sector.

    The linear model is selected whenever chi2/ndf <= 2.5. If it is inadequate,
    the quadratic model is selected when available. A quadratic may also replace
    a marginal linear fit when its AICc improves by at least 2.
    """
    fit_frame = pd.DataFrame(fit_records)
    model_records: list[dict[str, Any]] = []

    for sector in range(1, 7):
        cells = fit_frame.loc[
            (fit_frame["sector"] == sector)
            & (fit_frame["cell_type"] == "theta")
            & (fit_frame["success"])
            & np.isfinite(fit_frame["mean_theta_deg"])
            & np.isfinite(fit_frame["delta_p_over_p"])
        ].sort_values("mean_theta_deg")

        linear = weighted_polynomial_fit(
            cells["mean_theta_deg"].to_numpy(),
            cells["delta_p_over_p"].to_numpy(),
            cells["delta_p_over_p_error"].to_numpy(),
            degree=1,
        )
        quadratic = weighted_polynomial_fit(
            cells["mean_theta_deg"].to_numpy(),
            cells["delta_p_over_p"].to_numpy(),
            cells["delta_p_over_p_error"].to_numpy(),
            degree=2,
        )

        selected_name = "none"
        selected = None
        selection_reason = "no_successful_model"

        if linear["success"]:
            linear_adequate = (
                np.isfinite(linear["chi2_ndf"])
                and linear["chi2_ndf"] <= linear_chi2_ndf_max
            )

            quadratic_materially_better = (
                quadratic["success"]
                and np.isfinite(linear["aicc"])
                and np.isfinite(quadratic["aicc"])
                and quadratic["aicc"]
                <= linear["aicc"] - quadratic_aicc_improvement_min
            )

            if linear_adequate and not quadratic_materially_better:
                selected_name = "linear"
                selected = linear
                selection_reason = "lowest_order_adequate"
            elif quadratic["success"]:
                selected_name = "quadratic"
                selected = quadratic
                selection_reason = (
                    "quadratic_aicc_improvement"
                    if quadratic_materially_better
                    else "linear_inadequate"
                )
            else:
                selected_name = "linear"
                selected = linear
                selection_reason = "quadratic_unavailable"
            # endif
        elif quadratic["success"]:
            selected_name = "quadratic"
            selected = quadratic
            selection_reason = "linear_unavailable"
        # endif

        theta_min = (
            float(cells["theta_low_deg"].min())
            if not cells.empty
            else math.nan
        )
        theta_max = (
            float(cells["theta_high_deg"].max())
            if not cells.empty
            else math.nan
        )

        model_records.append(
            {
                "period_key": period.key,
                "period_label": period.label,
                "sector": sector,
                "success": selected is not None,
                "selected_model": selected_name,
                "selected_degree": (
                    int(selected["degree"])
                    if selected is not None
                    else None
                ),
                "selection_reason": selection_reason,
                "theta_valid_min_deg": theta_min,
                "theta_valid_max_deg": theta_max,
                "selected_coefficients_ascending": (
                    selected["coefficients_ascending"]
                    if selected is not None
                    else []
                ),
                "selected_covariance": (
                    selected["covariance"]
                    if selected is not None
                    else []
                ),
                "selected_chi2": (
                    selected["chi2"] if selected is not None else math.nan
                ),
                "selected_ndf": (
                    selected["ndf"] if selected is not None else 0
                ),
                "selected_chi2_ndf": (
                    selected["chi2_ndf"]
                    if selected is not None
                    else math.nan
                ),
                "linear": linear,
                "quadratic": quadratic,
                "n_theta_cells": int(len(cells)),
            }
        )
    # endfor

    return model_records


def evaluate_selected_theta_model(
    theta_deg: np.ndarray,
    model_record: dict[str, Any],
) -> np.ndarray:
    """Evaluate the selected polynomial with theta clipped to its valid range."""
    theta = np.asarray(theta_deg, dtype=float)
    coefficients = np.asarray(
        model_record["selected_coefficients_ascending"],
        dtype=float,
    )

    if coefficients.size == 0:
        return np.full_like(theta, np.nan, dtype=float)
    # endif

    theta_clipped = np.clip(
        theta,
        model_record["theta_valid_min_deg"],
        model_record["theta_valid_max_deg"],
    )

    result = np.zeros_like(theta_clipped, dtype=float)
    for power, coefficient in enumerate(coefficients):
        result += coefficient * theta_clipped**power
    # endfor
    return result


def apply_theta_corrections(
    frame: pd.DataFrame,
    model_records: list[dict[str, Any]],
) -> pd.DataFrame:
    """
    Apply the selected correction model event by event and recompute W.

    Sign convention:
        p_corrected = p_measured * (1 + delta_p_over_p)
    """
    corrected = frame.copy()
    corrected["model_delta_p_over_p"] = np.nan

    model_by_sector = {
        int(record["sector"]): record
        for record in model_records
        if record["success"]
    }

    for sector, model_record in model_by_sector.items():
        sector_mask = corrected["sector"] == sector
        corrected.loc[
            sector_mask,
            "model_delta_p_over_p",
        ] = evaluate_selected_theta_model(
            corrected.loc[sector_mask, "e_theta_deg"].to_numpy(),
            model_record,
        )
    # endfor

    if corrected["model_delta_p_over_p"].isna().any():
        missing_sectors = sorted(
            corrected.loc[
                corrected["model_delta_p_over_p"].isna(),
                "sector",
            ].unique()
        )
        raise RuntimeError(
            "No selected theta correction model for sectors "
            f"{missing_sectors}"
        )
    # endif

    corrected["e_p_uncorrected_gev"] = corrected["e_p"]
    corrected["e_p"] = corrected["e_p"] * (
        1.0 + corrected["model_delta_p_over_p"]
    )
    corrected["e_p_corrected_gev"] = corrected["e_p"]

    theta = corrected["e_theta_rad"].to_numpy()
    momentum = corrected["e_p"].to_numpy()
    beam_energy = corrected["beam_energy_gev"].to_numpy()

    q2 = (
        2.0
        * beam_energy
        * momentum
        * (1.0 - np.cos(theta))
    )
    w2 = (
        PROTON_MASS_GEV**2
        + 2.0 * PROTON_MASS_GEV * (beam_energy - momentum)
        - q2
    )
    corrected["w_corrected"] = np.sqrt(
        np.clip(w2, a_min=0.0, a_max=None)
    )
    corrected["w_recalculated"] = corrected["w_corrected"]
    return corrected


def theta_edges_from_uncorrected_records(
    fit_records: list[dict[str, Any]],
    sector: int,
) -> np.ndarray:
    """Recover the exact original theta boundaries for one sector."""
    sector_cells = pd.DataFrame(fit_records)
    sector_cells = sector_cells.loc[
        (sector_cells["sector"] == sector)
        & (sector_cells["cell_type"] == "theta")
    ].sort_values("theta_bin_index")

    if sector_cells.empty:
        return np.asarray([], dtype=float)
    # endif

    lows = sector_cells["theta_low_deg"].to_numpy(dtype=float)
    final_high = float(sector_cells["theta_high_deg"].iloc[-1])
    return np.concatenate([lows, [final_high]])


def fit_corrected_closure_cells(
    corrected_frame: pd.DataFrame,
    period: CalibrationPeriod,
    uncorrected_records: list[dict[str, Any]],
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    save_individual_fits: bool,
    individual_fit_directory: Path,
) -> list[dict[str, Any]]:
    """
    Refit corrected W using the exact theta edges from the uncorrected sample.
    """
    records: list[dict[str, Any]] = []

    for sector in range(1, 7):
        sector_frame = corrected_frame.loc[
            corrected_frame["sector"] == sector
        ].copy()
        theta_edges = theta_edges_from_uncorrected_records(
            uncorrected_records,
            sector,
        )

        integrated_fit, integrated_diagnostics = fit_w_peak(
            sector_frame["w_corrected"].to_numpy(),
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events_integrated,
        )
        integrated_record = make_fit_record(
            period=period,
            sector=sector,
            cell_type="sector_integrated",
            theta_bin_index=None,
            theta_low_deg=float(theta_edges[0]),
            theta_high_deg=float(theta_edges[-1]),
            local_phi_bin_index=None,
            local_phi_low_deg=0.0,
            local_phi_high_deg=60.0,
            cell_frame=sector_frame,
            fit_result=integrated_fit,
        )
        integrated_record["analysis_stage"] = "corrected"
        records.append(integrated_record)

        if save_individual_fits:
            save_peak_fit_plot(
                output_path=individual_fit_directory
                / f"{period.key}_sector{sector}_integrated_corrected.png",
                diagnostics=integrated_diagnostics,
                fit_result=integrated_fit,
                title=(
                    f"{period.label}, sector {sector}, corrected integrated W"
                ),
            )
        # endif

        for theta_index in range(len(theta_edges) - 1):
            theta_low = float(theta_edges[theta_index])
            theta_high = float(theta_edges[theta_index + 1])
            if theta_index == len(theta_edges) - 2:
                mask = (
                    (sector_frame["e_theta_deg"] >= theta_low)
                    & (sector_frame["e_theta_deg"] <= theta_high)
                )
            else:
                mask = (
                    (sector_frame["e_theta_deg"] >= theta_low)
                    & (sector_frame["e_theta_deg"] < theta_high)
                )
            # endif

            cell_frame = sector_frame.loc[mask].copy()
            if integrated_fit.success:
                cell_fit, diagnostics = fit_w_peak(
                    cell_frame["w_corrected"].to_numpy(),
                    histogram_range=histogram_range,
                    fit_range=fit_range,
                    histogram_bins=histogram_bins,
                    minimum_events=minimum_events_cell,
                    fixed_sigma_gev=integrated_fit.peak_sigma_gev,
                    mean_center_gev=PROTON_MASS_GEV,
                    mean_half_window_gev=0.040,
                )
            else:
                cell_fit = failed_peak_fit(
                    status="skipped_corrected_integrated_fit_failed",
                    n_events=int(len(cell_frame)),
                    fit_range=fit_range,
                    fit_mode="cell_fixed_sigma",
                    sigma_fixed=True,
                    mean_constraint_center_gev=PROTON_MASS_GEV,
                    mean_constraint_half_width_gev=0.040,
                )
                diagnostics = {
                    "counts": np.asarray([], dtype=float),
                    "edges": np.asarray([], dtype=float),
                    "centers": np.asarray([], dtype=float),
                    "fit_mask": np.asarray([], dtype=bool),
                    "fit_centers": np.asarray([], dtype=float),
                    "fit_counts": np.asarray([], dtype=float),
                    "fit_errors": np.asarray([], dtype=float),
                    "fit_model_at_centers": np.asarray([], dtype=float),
                    "pulls": np.asarray([], dtype=float),
                    "fit_x": np.asarray([], dtype=float),
                    "fit_y": np.asarray([], dtype=float),
                    "fit_signal_y": np.asarray([], dtype=float),
                    "fit_background_y": np.asarray([], dtype=float),
                }
            # endif

            record = make_fit_record(
                period=period,
                sector=sector,
                cell_type="theta",
                theta_bin_index=theta_index,
                theta_low_deg=theta_low,
                theta_high_deg=theta_high,
                local_phi_bin_index=None,
                local_phi_low_deg=0.0,
                local_phi_high_deg=60.0,
                cell_frame=cell_frame,
                fit_result=cell_fit,
            )
            record["analysis_stage"] = "corrected"
            records.append(record)

            if save_individual_fits:
                save_peak_fit_plot(
                    output_path=individual_fit_directory
                    / (
                        f"{period.key}_sector{sector}_theta"
                        f"{theta_index:02d}_corrected.png"
                    ),
                    diagnostics=diagnostics,
                    fit_result=cell_fit,
                    title=(
                        f"{period.label}, sector {sector}, corrected theta "
                        f"cell {theta_index + 1}/{len(theta_edges) - 1}"
                    ),
                )
            # endif
        # endfor
    # endfor

    return records


def save_model_selection_plot(
    period: CalibrationPeriod,
    fit_records: list[dict[str, Any]],
    model_records: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Plot calibration cells with linear, quadratic, and selected models."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fit_frame = pd.DataFrame(fit_records)
    model_by_sector = {
        int(record["sector"]): record for record in model_records
    }

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        constrained_layout=True,
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        cells = fit_frame.loc[
            (fit_frame["sector"] == sector)
            & (fit_frame["cell_type"] == "theta")
            & (fit_frame["success"])
        ].sort_values("mean_theta_deg")
        model = model_by_sector[sector]

        axis.errorbar(
            cells["mean_theta_deg"],
            100.0 * cells["delta_p_over_p"],
            yerr=100.0 * cells["delta_p_over_p_error"],
            fmt="o",
            label="Fitted cells",
        )

        if model["success"]:
            theta_grid = np.linspace(
                model["theta_valid_min_deg"],
                model["theta_valid_max_deg"],
                300,
            )

            for model_name, linestyle in (
                ("linear", "--"),
                ("quadratic", ":"),
            ):
                candidate = model[model_name]
                if not candidate["success"]:
                    continue
                # endif
                coefficients = np.asarray(
                    candidate["coefficients_ascending"],
                    dtype=float,
                )
                values = np.zeros_like(theta_grid)
                for power, coefficient in enumerate(coefficients):
                    values += coefficient * theta_grid**power
                # endfor
                axis.plot(
                    theta_grid,
                    100.0 * values,
                    linestyle=linestyle,
                    linewidth=1.2,
                    label=(
                        f"{model_name.capitalize()} "
                        f"($\\chi^2$/ndf={candidate['chi2_ndf']:.2f})"
                    ),
                )
            # endfor

            selected_values = evaluate_selected_theta_model(
                theta_grid,
                model,
            )
            axis.plot(
                theta_grid,
                100.0 * selected_values,
                linewidth=2.0,
                label=f"Selected: {model['selected_model']}",
            )
        # endif

        axis.axhline(0.0, color="black", linestyle="--", linewidth=0.8)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Mean electron theta (deg)")
        axis.set_ylabel("Momentum correction (%)")
        if sector == 1:
            axis.legend(fontsize=7)
        # endif
    # endfor

    figure.suptitle(
        f"{period.label}: linear and quadratic theta-correction models"
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_before_after_theta_plot(
    period: CalibrationPeriod,
    uncorrected_records: list[dict[str, Any]],
    corrected_records: list[dict[str, Any]],
    output_path: Path,
) -> None:
    """Compare fitted W-peak residuals before and after correction."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    before = pd.DataFrame(uncorrected_records)
    after = pd.DataFrame(corrected_records)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharey=True,
        constrained_layout=True,
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        before_cells = before.loc[
            (before["sector"] == sector)
            & (before["cell_type"] == "theta")
            & (before["success"])
        ].sort_values("mean_theta_deg")
        after_cells = after.loc[
            (after["sector"] == sector)
            & (after["cell_type"] == "theta")
            & (after["success"])
        ].sort_values("mean_theta_deg")

        axis.errorbar(
            before_cells["mean_theta_deg"],
            1000.0 * (
                before_cells["peak_mean_gev"] - PROTON_MASS_GEV
            ),
            yerr=1000.0 * before_cells["peak_mean_error_gev"],
            fmt="o",
            label="Before correction",
        )
        axis.errorbar(
            after_cells["mean_theta_deg"],
            1000.0 * (
                after_cells["peak_mean_gev"] - PROTON_MASS_GEV
            ),
            yerr=1000.0 * after_cells["peak_mean_error_gev"],
            fmt="s",
            label="After correction",
        )
        axis.axhline(0.0, color="black", linestyle="--", linewidth=1.0)
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Mean electron theta (deg)")
        axis.set_ylabel(r"$\mu_W-m_p$ (MeV)")
        if sector == 1:
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(
        f"{period.label}: electron momentum-correction closure versus theta"
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_before_after_integrated_w_plot(
    uncorrected_frame: pd.DataFrame,
    corrected_frame: pd.DataFrame,
    period: CalibrationPeriod,
    output_path: Path,
    histogram_range: tuple[float, float],
    histogram_bins: int,
) -> None:
    """Overlay sector-integrated W before and after correction."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        constrained_layout=True,
    )

    for sector, axis in zip(range(1, 7), axes.ravel()):
        before = uncorrected_frame.loc[
            uncorrected_frame["sector"] == sector,
            "w_recalculated",
        ].to_numpy()
        after = corrected_frame.loc[
            corrected_frame["sector"] == sector,
            "w_corrected",
        ].to_numpy()

        counts_before, edges = np.histogram(
            before,
            bins=histogram_bins,
            range=histogram_range,
        )
        counts_after, _ = np.histogram(
            after,
            bins=histogram_bins,
            range=histogram_range,
        )
        centers = 0.5 * (edges[:-1] + edges[1:])

        axis.step(
            centers,
            counts_before,
            where="mid",
            linewidth=1.2,
            label="Before correction",
        )
        axis.step(
            centers,
            counts_after,
            where="mid",
            linewidth=1.2,
            label="After correction",
        )
        axis.axvline(
            PROTON_MASS_GEV,
            color="black",
            linestyle="--",
            linewidth=1.0,
            label=r"$m_p$" if sector == 1 else None,
        )
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("W (GeV)")
        axis.set_ylabel("Counts")
        if sector == 1:
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(f"{period.label}: integrated W before and after correction")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def build_closure_summary(
    period: CalibrationPeriod,
    uncorrected_records: list[dict[str, Any]],
    corrected_records: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Calculate before/after closure metrics for every sector."""
    before = pd.DataFrame(uncorrected_records)
    after = pd.DataFrame(corrected_records)
    summary: list[dict[str, Any]] = []

    for sector in range(1, 7):
        before_cells = before.loc[
            (before["sector"] == sector)
            & (before["cell_type"] == "theta")
            & (before["success"])
        ]
        after_cells = after.loc[
            (after["sector"] == sector)
            & (after["cell_type"] == "theta")
            & (after["success"])
        ]

        before_residual = (
            before_cells["peak_mean_gev"].to_numpy() - PROTON_MASS_GEV
        )
        after_residual = (
            after_cells["peak_mean_gev"].to_numpy() - PROTON_MASS_GEV
        )

        summary.append(
            {
                "period_key": period.key,
                "period_label": period.label,
                "sector": sector,
                "n_before_cells": int(len(before_cells)),
                "n_after_cells": int(len(after_cells)),
                "before_mean_abs_residual_mev": (
                    1000.0 * float(np.mean(np.abs(before_residual)))
                    if before_residual.size
                    else math.nan
                ),
                "after_mean_abs_residual_mev": (
                    1000.0 * float(np.mean(np.abs(after_residual)))
                    if after_residual.size
                    else math.nan
                ),
                "before_rms_residual_mev": (
                    1000.0 * float(np.sqrt(np.mean(before_residual**2)))
                    if before_residual.size
                    else math.nan
                ),
                "after_rms_residual_mev": (
                    1000.0 * float(np.sqrt(np.mean(after_residual**2)))
                    if after_residual.size
                    else math.nan
                ),
            }
        )
    # endfor

    return summary



# -----------------------------------------------------------------------------
# Parallel period worker
# -----------------------------------------------------------------------------

def process_calibration_period(
    period: CalibrationPeriod,
    file_path: Path,
    tree_name_requested: str,
    step_size: str,
    theta_min_deg: float,
    theta_max_deg: float,
    w_preselection_min_gev: float,
    w_preselection_max_gev: float,
    theta_bin_count: int,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events_integrated: int,
    minimum_events_cell: int,
    save_individual_fits: bool,
    plot_directories: dict[str, Path],
    skip_run_stability: bool,
    run_bin_width: int,
) -> dict[str, Any]:
    """
    Process one calibration period completely inside one worker.

    Event-level data remain local to the worker. Only compact fit, model, and
    closure records are returned to the parent process.
    """
    tree_name = find_tree(file_path, tree_name_requested)
    branches = validate_branches(file_path, tree_name)

    print(
        f"[worker:{period.key}] reading {file_path}:{tree_name}, "
        f"runs {period.run_min}-{period.run_max}",
        flush=True,
    )

    period_frame = read_selected_events(
        file_path=file_path,
        tree_name=tree_name,
        branches=branches,
        period=period,
        step_size=step_size,
        theta_min_deg=theta_min_deg,
        theta_max_deg=theta_max_deg,
        w_preselection_min_gev=w_preselection_min_gev,
        w_preselection_max_gev=w_preselection_max_gev,
    )

    if period_frame.empty:
        return {
            "period_key": period.key,
            "period_label": period.label,
            "success": False,
            "status": "no_selected_events",
            "n_events": 0,
            "run_min_selected": None,
            "run_max_selected": None,
            "uncorrected_fit_records": [],
            "corrected_fit_records": [],
            "model_records": [],
            "closure_records": [],
        }
    # endif

    uncorrected_records = fit_calibration_cells(
        frame=period_frame,
        period=period,
        theta_bin_count=theta_bin_count,
        histogram_range=histogram_range,
        fit_range=fit_range,
        histogram_bins=histogram_bins,
        minimum_events_integrated=minimum_events_integrated,
        minimum_events_cell=minimum_events_cell,
        save_individual_fits=save_individual_fits,
        individual_fit_directory=plot_directories[
            "uncorrected_individual_fits"
        ],
    )
    for record in uncorrected_records:
        record["analysis_stage"] = "uncorrected"
    # endfor

    model_records = select_theta_correction_models(
        period,
        uncorrected_records,
    )
    unsuccessful_models = [
        record["sector"]
        for record in model_records
        if not record["success"]
    ]
    if unsuccessful_models:
        raise RuntimeError(
            f"Unable to select correction models for sectors "
            f"{unsuccessful_models} in {period.label}"
        )
    # endif

    corrected_frame = apply_theta_corrections(
        period_frame,
        model_records,
    )
    corrected_records = fit_corrected_closure_cells(
        corrected_frame=corrected_frame,
        period=period,
        uncorrected_records=uncorrected_records,
        histogram_range=histogram_range,
        fit_range=fit_range,
        histogram_bins=histogram_bins,
        minimum_events_integrated=minimum_events_integrated,
        minimum_events_cell=minimum_events_cell,
        save_individual_fits=save_individual_fits,
        individual_fit_directory=plot_directories[
            "corrected_individual_fits"
        ],
    )
    closure_records = build_closure_summary(
        period,
        uncorrected_records,
        corrected_records,
    )

    uncorrected_fit_frame = pd.DataFrame(uncorrected_records)
    corrected_fit_frame = pd.DataFrame(corrected_records)

    save_period_integrated_w_plot(
        frame=period_frame,
        period=period,
        fit_frame=uncorrected_fit_frame,
        output_path=plot_directories["uncorrected_integrated_w"]
        / f"{period.key}_integrated_w_by_sector.png",
        histogram_range=histogram_range,
        fit_range=fit_range,
        histogram_bins=histogram_bins,
        minimum_events=minimum_events_integrated,
        w_column="w_recalculated",
        stage_label="before correction",
    )
    save_period_integrated_w_plot(
        frame=corrected_frame,
        period=period,
        fit_frame=corrected_fit_frame,
        output_path=plot_directories["corrected_integrated_w"]
        / f"{period.key}_integrated_w_by_sector_corrected.png",
        histogram_range=histogram_range,
        fit_range=fit_range,
        histogram_bins=histogram_bins,
        minimum_events=minimum_events_integrated,
        w_column="w_corrected",
        stage_label="after correction",
    )

    save_model_selection_plot(
        period,
        uncorrected_records,
        model_records,
        plot_directories["model_selection"]
        / f"{period.key}_theta_model_selection.png",
    )
    save_before_after_theta_plot(
        period,
        uncorrected_records,
        corrected_records,
        plot_directories["closure_theta"]
        / f"{period.key}_theta_closure.png",
    )
    save_before_after_integrated_w_plot(
        period_frame,
        corrected_frame,
        period,
        plot_directories["before_after_w"]
        / f"{period.key}_integrated_w_before_after.png",
        histogram_range,
        histogram_bins,
    )

    if not skip_run_stability:
        save_run_stability_plot(
            frame=period_frame,
            period=period,
            output_path=plot_directories["uncorrected_run_stability"]
            / f"{period.key}_w_peak_vs_run.png",
            run_bin_width=run_bin_width,
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events_cell,
            w_column="w_recalculated",
            stage_label="before correction",
        )
        save_run_stability_plot(
            frame=corrected_frame,
            period=period,
            output_path=plot_directories["corrected_run_stability"]
            / f"{period.key}_w_peak_vs_run_corrected.png",
            run_bin_width=run_bin_width,
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events_cell,
            w_column="w_corrected",
            stage_label="after correction",
        )
    # endif

    uncorrected_theta = uncorrected_fit_frame.loc[
        uncorrected_fit_frame["cell_type"] == "theta"
    ]
    corrected_theta = corrected_fit_frame.loc[
        corrected_fit_frame["cell_type"] == "theta"
    ]

    print(
        f"[worker:{period.key}] retained {len(period_frame):,} events; "
        f"uncorrected theta fits "
        f"{int(uncorrected_theta['success'].sum())}/"
        f"{len(uncorrected_theta)}, corrected closure fits "
        f"{int(corrected_theta['success'].sum())}/"
        f"{len(corrected_theta)}",
        flush=True,
    )

    return {
        "period_key": period.key,
        "period_label": period.label,
        "success": True,
        "status": "success",
        "n_events": int(len(period_frame)),
        "run_min_selected": int(period_frame["runnum"].min()),
        "run_max_selected": int(period_frame["runnum"].max()),
        "uncorrected_fit_records": uncorrected_records,
        "corrected_fit_records": corrected_records,
        "model_records": model_records,
        "closure_records": closure_records,
    }


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------


def save_peak_fit_plot(
    output_path: Path,
    diagnostics: dict[str, np.ndarray],
    fit_result: PeakFitResult,
    title: str,
) -> None:
    """Save a W-fit overlay and a pull panel for manual review."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, (axis, pull_axis) = plt.subplots(
        2,
        1,
        figsize=(8.2, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.05},
    )
    centers = diagnostics["centers"]
    counts = diagnostics["counts"]
    errors = np.sqrt(np.maximum(counts, 1.0))

    axis.errorbar(
        centers,
        counts,
        yerr=errors,
        fmt=".",
        markersize=4,
        linewidth=0.8,
        label="Data",
    )

    if diagnostics["fit_x"].size:
        axis.plot(
            diagnostics["fit_x"],
            diagnostics["fit_y"],
            linewidth=1.6,
            label="Total fit",
        )
        axis.plot(
            diagnostics["fit_x"],
            diagnostics["fit_background_y"],
            linestyle="--",
            linewidth=1.2,
            label="Quadratic background",
        )
        axis.plot(
            diagnostics["fit_x"],
            diagnostics["fit_signal_y"]
            + diagnostics["fit_background_y"],
            linestyle=":",
            linewidth=1.2,
            label="Elastic component + background",
        )

        fit_centers = diagnostics["fit_centers"]
        pulls = diagnostics["pulls"]
        pull_axis.axhline(0.0, color="black", linewidth=1.0)
        pull_axis.axhline(3.0, color="black", linestyle=":", linewidth=0.8)
        pull_axis.axhline(-3.0, color="black", linestyle=":", linewidth=0.8)
        pull_axis.plot(
            fit_centers,
            pulls,
            marker="o",
            linestyle="none",
            markersize=3,
        )
        pull_axis.set_ylim(-6.0, 6.0)

        sigma_description = (
            "fixed from sector-integrated fit"
            if fit_result.sigma_fixed
            else "floating"
        )
        quality_description = (
            "automatic checks: pass"
            if fit_result.quality_pass
            else f"review flags: {fit_result.quality_flags}"
        )
        annotation = (
            f"status: {fit_result.status}\n"
            f"$\\mu_W$ = {fit_result.peak_mean_gev:.6f} +/- "
            f"{fit_result.peak_mean_error_gev:.6f} GeV\n"
            f"$\\sigma_W$ = {fit_result.peak_sigma_gev:.6f} GeV "
            f"({sigma_description})\n"
            f"$\\chi^2$/ndf = {fit_result.chi2_ndf:.2f}\n"
            f"N = {fit_result.n_events:,}\n"
            f"{quality_description}"
        )
    else:
        pull_axis.axhline(0.0, color="black", linewidth=1.0)
        pull_axis.text(
            0.5,
            0.5,
            "Fit unavailable",
            transform=pull_axis.transAxes,
            ha="center",
            va="center",
        )
        annotation = (
            f"Fit status: {fit_result.status}\n"
            f"N = {fit_result.n_events:,}"
        )
    # endif

    axis.axvline(
        PROTON_MASS_GEV,
        linestyle="--",
        linewidth=1.0,
        color="black",
        label=r"$m_p$",
    )
    axis.set_ylabel("Counts")
    axis.set_title(title)
    axis.text(
        0.03,
        0.97,
        annotation,
        transform=axis.transAxes,
        ha="left",
        va="top",
        fontsize=8,
    )
    axis.legend(fontsize=8)

    pull_axis.set_xlabel("W (GeV)")
    pull_axis.set_ylabel("Pull")
    figure.subplots_adjust(
        left=0.11,
        right=0.97,
        bottom=0.10,
        top=0.92,
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)



def save_period_integrated_w_plot(
    frame: pd.DataFrame,
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    output_path: Path,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
    w_column: str = "w_recalculated",
    stage_label: str = "uncorrected",
) -> None:
    """Save sector-integrated W spectra with the fitted model overlaid."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(2, 3, figsize=(15.0, 9.0), sharex=True)
    axes_flat = axes.ravel()

    for sector in range(1, 7):
        axis = axes_flat[sector - 1]
        sector_values = frame.loc[
            frame["sector"] == sector,
            w_column,
        ].to_numpy()

        fit_result, diagnostics = fit_w_peak(
            sector_values,
            histogram_range=histogram_range,
            fit_range=fit_range,
            histogram_bins=histogram_bins,
            minimum_events=minimum_events,
        )

        centers = diagnostics["centers"]
        counts = diagnostics["counts"]
        errors = np.sqrt(np.maximum(counts, 1.0))
        axis.errorbar(
            centers,
            counts,
            yerr=errors,
            fmt=".",
            markersize=2,
            linewidth=0.6,
            label="Data",
        )

        if diagnostics["fit_x"].size:
            axis.plot(
                diagnostics["fit_x"],
                diagnostics["fit_y"],
                linewidth=1.4,
                label="Total fit",
            )
            axis.plot(
                diagnostics["fit_x"],
                diagnostics["fit_background_y"],
                linestyle="--",
                linewidth=1.1,
                label="Background",
            )
            axis.plot(
                diagnostics["fit_x"],
                diagnostics["fit_signal_y"]
                + diagnostics["fit_background_y"],
                linestyle=":",
                linewidth=1.1,
                label="Elastic + background",
            )
        # endif

        axis.axvline(PROTON_MASS_GEV, linestyle="--", linewidth=1.0)
        if np.isfinite(fit_result.peak_mean_gev):
            axis.set_title(
                f"Sector {sector}: $\\mu_W$={fit_result.peak_mean_gev:.5f} GeV\n"
                f"{fit_result.status}, chi2/ndf={fit_result.chi2_ndf:.2f}"
            )
        else:
            axis.set_title(f"Sector {sector}: {fit_result.status}")
        # endif

        axis.set_xlabel("W (GeV)")
        axis.set_ylabel("Counts")
        if sector == 1:
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(
        f"{period.label}: inclusive eX elastic-peak fits ({stage_label})"
    )
    figure.savefig(output_path, dpi=180)
    plt.close(figure)



def save_peak_vs_theta_plot(
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot the fitted W peak versus electron polar angle for all sectors."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes_flat = axes.ravel()

    for sector in range(1, 7):
        axis = axes_flat[sector - 1]
        cells = fit_frame.loc[
            (fit_frame["period_key"] == period.key)
            & (fit_frame["sector"] == sector)
            & (fit_frame["cell_type"] == "theta")
            & (fit_frame["success"])
        ].sort_values("mean_theta_deg")

        if cells.empty:
            axis.text(
                0.5,
                0.5,
                "No converged fits",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
        else:
            clean_cells = cells.loc[cells["quality_pass"]]
            flagged_cells = cells.loc[~cells["quality_pass"]]

            if not clean_cells.empty:
                axis.errorbar(
                    clean_cells["mean_theta_deg"],
                    clean_cells["peak_mean_gev"],
                    yerr=clean_cells["peak_mean_error_gev"],
                    fmt="o",
                    markersize=5,
                    linewidth=1.0,
                    label="Passes automatic checks",
                )
            # endif

            if not flagged_cells.empty:
                axis.errorbar(
                    flagged_cells["mean_theta_deg"],
                    flagged_cells["peak_mean_gev"],
                    yerr=flagged_cells["peak_mean_error_gev"],
                    fmt="x",
                    markersize=6,
                    linewidth=1.0,
                    label="Flagged for manual review",
                )
            # endif
        # endif

        axis.axhline(
            PROTON_MASS_GEV,
            linestyle="--",
            linewidth=1.0,
            color="black",
        )
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Mean electron theta (deg)")
        axis.set_ylabel("Fitted W peak (GeV)")
        if sector == 1 and not cells.empty:
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(f"{period.label}: fitted elastic W peak versus electron theta")
    figure.savefig(output_path, dpi=180)
    plt.close(figure)


def save_correction_vs_theta_plot(
    period: CalibrationPeriod,
    fit_frame: pd.DataFrame,
    output_path: Path,
) -> None:
    """Plot provisional delta-p/p values versus electron polar angle."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    figure, axes = plt.subplots(
        2,
        3,
        figsize=(15.0, 9.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
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

        if cells.empty:
            axis.text(
                0.5,
                0.5,
                "No converged fits",
                ha="center",
                va="center",
                transform=axis.transAxes,
            )
        else:
            clean_cells = cells.loc[cells["quality_pass"]]
            flagged_cells = cells.loc[~cells["quality_pass"]]

            if not clean_cells.empty:
                axis.errorbar(
                    clean_cells["mean_theta_deg"],
                    100.0 * clean_cells["delta_p_over_p"],
                    yerr=100.0 * clean_cells["delta_p_over_p_error"],
                    fmt="o",
                    markersize=5,
                    linewidth=1.0,
                    label="Passes automatic checks",
                )
            # endif

            if not flagged_cells.empty:
                axis.errorbar(
                    flagged_cells["mean_theta_deg"],
                    100.0 * flagged_cells["delta_p_over_p"],
                    yerr=100.0 * flagged_cells["delta_p_over_p_error"],
                    fmt="x",
                    markersize=6,
                    linewidth=1.0,
                    label="Flagged for manual review",
                )
            # endif
        # endif

        axis.axhline(0.0, linestyle="--", linewidth=1.0, color="black")
        axis.set_title(f"Sector {sector}")
        axis.set_xlabel("Mean electron theta (deg)")
        axis.set_ylabel("Required momentum correction (%)")
        if sector == 1 and not cells.empty:
            axis.legend(fontsize=8)
        # endif
    # endfor

    figure.suptitle(
        f"{period.label}: correction inferred from fitted elastic W peak"
    )
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
        axis.set_xlabel("Local phi (deg)")
        axis.set_ylabel("Electron theta (deg)")
    # endfor

    if last_image is not None:
        colorbar = figure.colorbar(
            last_image,
            ax=axes_flat.tolist(),
            location="right",
            fraction=0.025,
            pad=0.045,
            shrink=0.88,
        )
        colorbar.set_label("Required momentum correction (%)")
    # endif

    figure.suptitle(
        f"{period.label}: correction inferred from fitted W peaks"
    )
    figure.subplots_adjust(
        left=0.07,
        right=0.82,
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
    w_column: str = "w_recalculated",
    stage_label: str = "uncorrected",
) -> None:
    """
    Fit W for every individual run by default.

    When run_bin_width is greater than one, consecutive run-number intervals are
    grouped. The x coordinate is the event-weighted mean run number. A dashed
    line in each sector's plotting color marks the arithmetic mean of that
    sector's successfully fitted run-level peak positions.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(11.0, 6.0))

    available_runs = np.sort(frame["runnum"].unique().astype(int))

    for sector in range(1, 7):
        run_centers: list[float] = []
        peak_positions: list[float] = []
        peak_errors: list[float] = []
        sector_frame = frame.loc[frame["sector"] == sector]

        if run_bin_width == 1:
            groups = [(run, run + 1) for run in available_runs]
        else:
            first_run = int(available_runs.min())
            last_run = int(available_runs.max())
            edges = np.arange(
                first_run,
                last_run + run_bin_width + 1,
                run_bin_width,
                dtype=int,
            )
            groups = list(zip(edges[:-1], edges[1:]))
        # endif

        for run_low, run_high in groups:
            run_frame = sector_frame.loc[
                (sector_frame["runnum"] >= run_low)
                & (sector_frame["runnum"] < run_high)
            ]

            fit_result, _ = fit_w_peak(
                run_frame[w_column].to_numpy(),
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
            errorbar_container = axis.errorbar(
                run_centers,
                peak_positions,
                yerr=peak_errors,
                fmt="o",
                markersize=3,
                linewidth=0.8,
                label=f"Sector {sector}",
            )
            sector_color = errorbar_container.lines[0].get_color()
            sector_mean_peak = float(np.mean(peak_positions))
            axis.axhline(
                sector_mean_peak,
                linestyle="--",
                linewidth=1.0,
                color=sector_color,
                alpha=0.9,
                label=(
                    rf"Sector {sector} mean "
                    rf"$\mu_W$={sector_mean_peak:.4f} GeV"
                ),
            )
        # endif
    # endfor

    axis.axhline(
        PROTON_MASS_GEV,
        linestyle="--",
        linewidth=1.0,
        color="black",
        label=r"$m_p$",
    )
    axis.set_xlabel("Run number")
    axis.set_ylabel("Fitted W peak (GeV)")
    point_definition = (
        "individual runs"
        if run_bin_width == 1
        else f"{run_bin_width}-run-number groups"
    )
    axis.set_title(
        f"{period.label}: elastic-peak stability versus run "
        f"({point_definition}, {stage_label})"
    )
    axis.legend(ncol=2)
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
            axis.set_xlabel("Mean electron theta (deg)")
            axis.set_ylabel("Residual: fitted cell - smooth model (%)")
        # endfor

        figure.suptitle(
            f"{period.label}: fitted-cell residuals relative to smooth model"
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
            "theta_base_bin_count_per_sector": arguments.theta_bin_count,
            "theta_final_bin_count_per_sector": arguments.theta_bin_count + 1,
            "theta_edges_definition": (
                "equal-population quantile base bins derived independently in "
                "each FD sector; the angularly widest base bin is split at the "
                "median theta of events inside that base bin"
            ),
            "local_phi_treatment": "integrated over 0-60 degrees",
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
        default=Path("output/electron_diagnostics"),
        help="Output directory. Default: output/electron_diagnostics",
    )
    parser.add_argument(
        "--step-size",
        default="250 MB",
        help="uproot iteration step size. Default: 250 MB",
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help=(
            "Number of calibration periods processed simultaneously. "
            "The maximum useful value is 4. Default: 4"
        ),
    )

    parser.add_argument(
        "--theta-min",
        type=float,
        default=5.0,
        help="Minimum electron theta (deg). Default: 5",
    )
    parser.add_argument(
        "--theta-max",
        type=float,
        default=25.0,
        help="Maximum electron theta (deg). Default: 25",
    )
    parser.add_argument(
        "--theta-bin-count",
        type=int,
        default=5,
        help=(
            "Number of equal-population base theta bins derived independently "
            "in each FD sector. The angularly widest base bin is split once, "
            "so the final number of fitted cells is this value plus one. "
            "Default: 5"
        ),
    )

    parser.add_argument(
        "--w-preselection-min",
        type=float,
        default=0.65,
        help="Minimum input W retained (GeV). Default: 0.65",
    )
    parser.add_argument(
        "--w-preselection-max",
        type=float,
        default=1.45,
        help="Maximum input W retained (GeV). Default: 1.45",
    )
    parser.add_argument(
        "--w-hist-min",
        type=float,
        default=0.70,
        help="Minimum W histogram edge (GeV). Default: 0.70",
    )
    parser.add_argument(
        "--w-hist-max",
        type=float,
        default=1.35,
        help="Maximum W histogram edge (GeV). Default: 1.35",
    )
    parser.add_argument(
        "--w-fit-min",
        type=float,
        default=0.82,
        help="Minimum W fit edge (GeV). Default: 0.82",
    )
    parser.add_argument(
        "--w-fit-max",
        type=float,
        default=1.12,
        help="Maximum W fit edge (GeV). Default: 1.12",
    )
    parser.add_argument(
        "--w-bins",
        type=int,
        default=60,
        help="Number of W histogram bins. Default: 60",
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
        default=2000,
        help="Minimum events for a sector-specific theta-bin fit. Default: 2000",
    )
    parser.add_argument(
        "--run-bin-width",
        type=int,
        default=1,
        help=(
            "Number of consecutive run numbers grouped into each stability "
            "point. The default of 1 fits every individual run."
        ),
    )
    parser.add_argument(
        "--skip-individual-fits",
        action="store_true",
        help=(
            "Do not save individual sector-integrated and theta-bin W-fit "
            "overlays. They are saved by default for manual review."
        ),
    )
    parser.add_argument(
        "--skip-run-stability",
        action="store_true",
        help="Skip the run-number stability plots.",
    )

    return parser


def validate_arguments(arguments: argparse.Namespace) -> None:
    """Validate command-line arguments and derive parsed edge arrays."""
    if arguments.theta_min >= arguments.theta_max:
        raise ValueError("--theta-min must be smaller than --theta-max.")
    # endif

    if arguments.theta_bin_count <= 0:
        raise ValueError("--theta-bin-count must be positive.")
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

    if arguments.workers <= 0:
        raise ValueError("--workers must be positive.")
    # endif

    if arguments.workers > 4:
        raise ValueError(
            "--workers cannot exceed 4 because there are only four "
            "independent calibration periods."
        )
    # endif


def main() -> int:
    """Run extraction, model selection, correction application, and closure."""
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
                f"FATAL: input file for {source_key} does not exist: "
                f"{file_path}",
                file=sys.stderr,
            )
            return 2
        # endif
    # endfor

    output_directory = arguments.output_dir
    csv_directory = output_directory / "csv"
    json_directory = output_directory / "json"
    plot_root = output_directory / "plots"

    # Rebuild generated products on every run so discontinued or renamed plots
    # cannot survive from an older script version.
    for generated_directory in (
        csv_directory,
        json_directory,
        plot_root,
    ):
        if generated_directory.exists():
            shutil.rmtree(generated_directory)
        # endif
        generated_directory.mkdir(parents=True, exist_ok=True)
    # endfor

    plot_directories = {
        "uncorrected_individual_fits": (
            plot_root / "uncorrected" / "individual_fits"
        ),
        "uncorrected_integrated_w": (
            plot_root / "uncorrected" / "integrated_w"
        ),
        "uncorrected_theta_trends": (
            plot_root / "uncorrected" / "theta_trends"
        ),
        "uncorrected_run_stability": (
            plot_root / "uncorrected" / "run_stability"
        ),
        "model_selection": plot_root / "models",
        "corrected_individual_fits": (
            plot_root / "closure" / "individual_fits"
        ),
        "corrected_integrated_w": (
            plot_root / "closure" / "integrated_w"
        ),
        "closure_theta": plot_root / "closure" / "theta_residuals",
        "before_after_w": plot_root / "closure" / "before_after_w",
        "corrected_run_stability": (
            plot_root / "closure" / "run_stability"
        ),
    }
    for directory in plot_directories.values():
        directory.mkdir(parents=True, exist_ok=True)
    # endfor

    all_uncorrected_records: list[dict[str, Any]] = []
    all_corrected_records: list[dict[str, Any]] = []
    all_model_records: list[dict[str, Any]] = []
    all_closure_records: list[dict[str, Any]] = []
    period_results: dict[str, dict[str, Any]] = {}

    histogram_range = (arguments.w_hist_min, arguments.w_hist_max)
    fit_range = (arguments.w_fit_min, arguments.w_fit_max)

    worker_count = min(arguments.workers, len(CALIBRATION_PERIODS))
    print(
        f"[parallel] processing {len(CALIBRATION_PERIODS)} calibration "
        f"periods with {worker_count} workers"
    )

    with ProcessPoolExecutor(max_workers=worker_count) as executor:
        future_to_period = {}

        for period in CALIBRATION_PERIODS:
            future = executor.submit(
                process_calibration_period,
                period,
                input_files[period.source_key],
                arguments.tree_name,
                arguments.step_size,
                arguments.theta_min,
                arguments.theta_max,
                arguments.w_preselection_min,
                arguments.w_preselection_max,
                arguments.theta_bin_count,
                histogram_range,
                fit_range,
                arguments.w_bins,
                arguments.minimum_events_integrated,
                arguments.minimum_events_cell,
                not arguments.skip_individual_fits,
                plot_directories,
                arguments.skip_run_stability,
                arguments.run_bin_width,
            )
            future_to_period[future] = period
        # endfor

        for future in as_completed(future_to_period):
            period = future_to_period[future]
            try:
                result = future.result()
            except Exception as exc:
                print(
                    f"FATAL: worker failed for {period.label}: "
                    f"{type(exc).__name__}: {exc}",
                    file=sys.stderr,
                )
                raise
            # endtry

            period_results[period.key] = result
            all_uncorrected_records.extend(
                result["uncorrected_fit_records"]
            )
            all_corrected_records.extend(
                result["corrected_fit_records"]
            )
            all_model_records.extend(result["model_records"])
            all_closure_records.extend(result["closure_records"])

            if result["success"]:
                print(
                    f"[parallel] {period.label}: retained "
                    f"{result['n_events']:,} events; runs "
                    f"{result['run_min_selected']}-"
                    f"{result['run_max_selected']}"
                )
            else:
                print(
                    f"WARNING: {period.label}: {result['status']}",
                    file=sys.stderr,
                )
            # endif
        # endfor
    # endwith

    if not all_uncorrected_records:
        print("FATAL: no calibration fits were attempted.", file=sys.stderr)
        return 3
    # endif

    uncorrected_frame = pd.DataFrame(all_uncorrected_records)
    corrected_frame = pd.DataFrame(all_corrected_records)
    model_frame = pd.json_normalize(all_model_records)
    closure_frame = pd.DataFrame(all_closure_records)

    uncorrected_csv = csv_directory / "fit_results_uncorrected.csv"
    corrected_csv = csv_directory / "fit_results_corrected.csv"
    models_csv = csv_directory / "correction_model_comparison.csv"
    closure_csv = csv_directory / "closure_summary.csv"

    uncorrected_frame.to_csv(
        uncorrected_csv,
        index=False,
        float_format="%.12g",
    )
    corrected_frame.to_csv(
        corrected_csv,
        index=False,
        float_format="%.12g",
    )
    model_frame.to_csv(models_csv, index=False, float_format="%.12g")
    closure_frame.to_csv(closure_csv, index=False, float_format="%.12g")

    for path in (
        uncorrected_csv,
        corrected_csv,
        models_csv,
        closure_csv,
    ):
        print(f"[write] {path}")
    # endfor

    cell_database = dataframe_to_cell_database(
        uncorrected_frame,
        arguments,
    )
    cells_json = json_directory / "electron_correction_cells.json"
    models_json = json_directory / "electron_correction_models.json"
    closure_json = json_directory / "closure_summary.json"

    write_json(cells_json, cell_database)
    write_json(
        models_json,
        {
            "metadata": {
                "model_selection": (
                    "select linear when chi2/ndf <= 2.5 unless quadratic "
                    "improves AICc by at least 2; otherwise select quadratic"
                ),
                "application_sign_convention": (
                    "p_corrected = p_measured * (1 + delta_p_over_p)"
                ),
                "theta_extrapolation_policy": (
                    "theta is clipped to the fitted sector validity range"
                ),
            },
            "models": all_model_records,
        },
    )
    write_json(
        closure_json,
        {
            "metadata": {
                "proton_mass_gev": PROTON_MASS_GEV,
                "theta_edges": (
                    "exact uncorrected extraction edges reused after correction"
                ),
            },
            "closure": all_closure_records,
        },
    )

    for path in (cells_json, models_json, closure_json):
        print(f"[write] {path}")
    # endfor

    for period in CALIBRATION_PERIODS:
        save_peak_vs_theta_plot(
            period=period,
            fit_frame=uncorrected_frame,
            output_path=plot_directories["uncorrected_theta_trends"]
            / f"{period.key}_w_peak_vs_theta.png",
        )
        save_correction_vs_theta_plot(
            period=period,
            fit_frame=uncorrected_frame,
            output_path=plot_directories["uncorrected_theta_trends"]
            / f"{period.key}_correction_vs_theta.png",
        )
        save_peak_vs_theta_plot(
            period=period,
            fit_frame=corrected_frame,
            output_path=plot_directories["closure_theta"]
            / f"{period.key}_corrected_w_peak_vs_theta.png",
        )
        save_correction_vs_theta_plot(
            period=period,
            fit_frame=corrected_frame,
            output_path=plot_directories["closure_theta"]
            / f"{period.key}_residual_correction_vs_theta.png",
        )
    # endfor

    uncorrected_theta = uncorrected_frame.loc[
        uncorrected_frame["cell_type"] == "theta"
    ]
    corrected_theta = corrected_frame.loc[
        corrected_frame["cell_type"] == "theta"
    ]
    selected_model_counts = pd.Series(
        [record["selected_model"] for record in all_model_records]
    ).value_counts()

    print("")
    print("Electron momentum-correction extraction and closure complete.")
    print(
        f"Uncorrected converged theta fits: "
        f"{int(uncorrected_theta['success'].sum()):,} / "
        f"{len(uncorrected_theta):,}"
    )
    print(
        f"Corrected closure theta fits: "
        f"{int(corrected_theta['success'].sum()):,} / "
        f"{len(corrected_theta):,}"
    )
    print(
        "Selected models: "
        + ", ".join(
            f"{name}={count}"
            for name, count in selected_model_counts.items()
        )
    )
    print(f"Output directory: {output_directory}")
    print(f"Period workers used: {worker_count}")
    print("")
    print(
        "The script stops after theta-dependent correction closure. "
        "No local-phi correction is derived or applied."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
# endif
