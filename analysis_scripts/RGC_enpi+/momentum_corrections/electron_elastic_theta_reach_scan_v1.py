#!/usr/bin/env python3
"""
electron_elastic_theta_reach_scan_v2.py

Determine the maximum defensible electron-theta reach of the inclusive elastic
eX calibration sample.

Relative to v1, this version:
  * reproduces the production fitting strategy:
      - sector-integrated Gaussian + quadratic fit with free sigma;
      - theta-cell fits with sigma fixed to the integrated-sector value;
      - theta-cell mean constrained around the integrated-sector centroid;
  * uses deliberately wider, adaptive high-theta bins;
  * performs fit-model stability tests:
      - nominal quadratic background and nominal fit range;
      - linear background;
      - narrower fit range;
      - wider fit range;
      - integrated sigma shifted by +/- its fitted uncertainty;
  * calculates signal yield and background under +/-2 sigma;
  * applies a provisional constant momentum correction in each theta bin,
    recomputes W event by event, and refits the corrected spectrum;
  * classifies a bin using nominal fit quality, fit-variation stability, and
    corrected-W closure;
  * saves every fit overlay at and above a configurable theta threshold;
  * reports the highest stable high-theta interval rather than requiring every
    lower fixed-width bin to pass.

This remains a diagnostic study. It does not overwrite the production
electron correction.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import math
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit


PROTON_MASS_GEV = 0.9382720813

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


@dataclass(frozen=True)
class Period:
    key: str
    label: str
    source_key: str
    run_min: int
    run_max: int


PERIODS = (
    Period("su22", "Su22", "su22", 16042, 16788),
    Period("fa22_solenoid_minus", "Fa22, solenoid -1", "fa22", 16843, 17183),
    Period("fa22_solenoid_plus", "Fa22, solenoid +1", "fa22", 17185, 17408),
    Period("sp23", "Sp23", "sp23", 17477, 17811),
)


def beam_energy_from_run(runnum: np.ndarray) -> np.ndarray:
    runnum = np.asarray(runnum, dtype=int)
    energy = np.full(runnum.shape, np.nan, dtype=float)
    energy[(runnum >= 16042) & (runnum <= 17065)] = 10.5473
    energy[(runnum >= 17067) & (runnum <= 17716)] = 10.5563
    energy[(runnum >= 17717) & (runnum <= 17811)] = 10.5593
    return energy


def fd_sector_from_phi(phi_rad: np.ndarray) -> np.ndarray:
    phi_deg = np.mod(np.degrees(phi_rad), 360.0)
    sector = np.full(phi_deg.shape, -1, dtype=np.int16)
    sector[(phi_deg >= 330.0) | (phi_deg < 30.0)] = 1
    sector[(phi_deg >= 30.0) & (phi_deg < 90.0)] = 2
    sector[(phi_deg >= 90.0) & (phi_deg < 150.0)] = 3
    sector[(phi_deg >= 150.0) & (phi_deg < 210.0)] = 4
    sector[(phi_deg >= 210.0) & (phi_deg < 270.0)] = 5
    sector[(phi_deg >= 270.0) & (phi_deg < 330.0)] = 6
    return sector


def momentum_at_w(w_gev: float, theta_rad: float, beam_energy_gev: float) -> float:
    numerator = (
        PROTON_MASS_GEV**2
        + 2.0 * PROTON_MASS_GEV * beam_energy_gev
        - w_gev**2
    )
    denominator = 2.0 * (
        PROTON_MASS_GEV
        + beam_energy_gev * (1.0 - math.cos(theta_rad))
    )
    return numerator / denominator


def elastic_scattered_momentum(theta_rad: float, beam_energy_gev: float) -> float:
    return (
        PROTON_MASS_GEV * beam_energy_gev
        / (
            PROTON_MASS_GEV
            + beam_energy_gev * (1.0 - math.cos(theta_rad))
        )
    )


def w_from_electron(
    momentum_gev: np.ndarray,
    theta_rad: np.ndarray,
    beam_energy_gev: np.ndarray,
) -> np.ndarray:
    q2 = (
        2.0
        * beam_energy_gev
        * momentum_gev
        * (1.0 - np.cos(theta_rad))
    )
    w2 = (
        PROTON_MASS_GEV**2
        + 2.0 * PROTON_MASS_GEV * (beam_energy_gev - momentum_gev)
        - q2
    )
    return np.sqrt(np.clip(w2, 0.0, None))


def parse_edges(text: str) -> np.ndarray:
    edges = np.asarray([float(value.strip()) for value in text.split(",")])
    if len(edges) < 2 or not np.all(np.diff(edges) > 0.0):
        raise argparse.ArgumentTypeError(
            "Theta edges must be a strictly increasing comma-separated list."
        )
    return edges


def gaussian_quadratic(x, amplitude, mean, sigma, b0, b1, b2):
    dx = x - PROTON_MASS_GEV
    return (
        amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        + b0 + b1 * dx + b2 * dx**2
    )


def gaussian_linear(x, amplitude, mean, sigma, b0, b1):
    dx = x - PROTON_MASS_GEV
    return (
        amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        + b0 + b1 * dx
    )


def fixed_sigma_quadratic(x, amplitude, mean, b0, b1, b2, sigma):
    dx = x - PROTON_MASS_GEV
    return (
        amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        + b0 + b1 * dx + b2 * dx**2
    )


def fixed_sigma_linear(x, amplitude, mean, b0, b1, sigma):
    dx = x - PROTON_MASS_GEV
    return (
        amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
        + b0 + b1 * dx
    )


def histogram_arrays(values, histogram_range, histogram_bins, fit_range):
    counts, edges = np.histogram(
        values,
        bins=histogram_bins,
        range=histogram_range,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (centers >= fit_range[0]) & (centers <= fit_range[1])
    x = centers[mask]
    y = counts[mask].astype(float)
    ey = np.sqrt(np.maximum(y, 1.0))
    return counts, edges, centers, mask, x, y, ey


def fit_peak(
    values: np.ndarray,
    histogram_range: tuple[float, float],
    fit_range: tuple[float, float],
    histogram_bins: int,
    minimum_events: int,
    fixed_sigma: float | None,
    mean_center: float | None,
    mean_half_window: float,
    background_order: int = 2,
) -> tuple[dict, dict]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    n_events = int(finite.size)

    counts, edges, centers, fit_mask, x, y, ey = histogram_arrays(
        finite, histogram_range, histogram_bins, fit_range
    )
    diagnostics = {
        "counts": counts,
        "edges": edges,
        "centers": centers,
        "fit_mask": fit_mask,
        "x": x,
        "y": y,
        "ey": ey,
        "model": np.asarray([]),
        "signal": np.asarray([]),
        "background": np.asarray([]),
    }

    result = {
        "success": False,
        "status": "not_fitted",
        "n_events": n_events,
        "peak_mean_gev": math.nan,
        "peak_mean_error_gev": math.nan,
        "peak_sigma_gev": math.nan,
        "peak_sigma_error_gev": math.nan,
        "amplitude": math.nan,
        "amplitude_error": math.nan,
        "chi2_ndf": math.nan,
        "signal_yield_2sigma": math.nan,
        "background_yield_2sigma": math.nan,
        "signal_to_background_2sigma": math.nan,
    }

    if n_events < minimum_events:
        result["status"] = "insufficient_events"
        return result, diagnostics

    minimum_points = 12 if fixed_sigma is not None else 14
    parameter_count = (
        5 if fixed_sigma is not None and background_order == 2 else
        4 if fixed_sigma is not None else
        6 if background_order == 2 else 5
    )
    if x.size < minimum_points or np.count_nonzero(y) < 10:
        result["status"] = "insufficient_populated_bins"
        return result, diagnostics

    center = (
        float(mean_center)
        if mean_center is not None
        else PROTON_MASS_GEV
    )
    if fixed_sigma is None:
        mean_low = max(fit_range[0], PROTON_MASS_GEV - 0.080)
        mean_high = min(fit_range[1], PROTON_MASS_GEV + 0.080)
        sideband_exclusion = 0.055
    else:
        if not np.isfinite(fixed_sigma) or fixed_sigma <= 0.0:
            result["status"] = "invalid_fixed_sigma"
            return result, diagnostics
        mean_low = max(fit_range[0], center - mean_half_window)
        mean_high = min(fit_range[1], center + mean_half_window)
        sideband_exclusion = max(0.045, 2.0 * fixed_sigma)

    peak_mask = (x >= mean_low) & (x <= mean_high)
    sideband_mask = (
        (x < center - sideband_exclusion)
        | (x > center + sideband_exclusion)
    )

    degree = background_order
    if np.count_nonzero(sideband_mask) >= degree + 5:
        try:
            coefficients = np.polyfit(
                x[sideband_mask] - PROTON_MASS_GEV,
                y[sideband_mask],
                deg=degree,
                w=1.0 / ey[sideband_mask],
            )
            background_seed = np.polyval(
                coefficients,
                x - PROTON_MASS_GEV,
            )
        except Exception:
            coefficients = None
            background_seed = np.full_like(y, np.percentile(y, 25.0))
    else:
        coefficients = None
        background_seed = np.full_like(y, np.percentile(y, 25.0))

    indices = np.flatnonzero(peak_mask)
    if indices.size == 0:
        result["status"] = "empty_peak_window"
        return result, diagnostics

    residual = y - background_seed
    peak_index = indices[np.argmax(residual[peak_mask])]
    amplitude0 = max(float(residual[peak_index]), 1.0)
    mean0 = float(np.clip(x[peak_index], mean_low, mean_high))

    if background_order == 2:
        if coefficients is not None:
            b2, b1, b0 = map(float, coefficients)
        else:
            b0, b1, b2 = float(np.percentile(y, 25.0)), 0.0, 0.0
    else:
        if coefficients is not None:
            b1, b0 = map(float, coefficients)
        else:
            b0, b1 = float(np.percentile(y, 25.0)), 0.0

    try:
        if fixed_sigma is None and background_order == 2:
            p0 = (amplitude0, mean0, 0.025, b0, b1, b2)
            bounds = (
                (0.0, mean_low, 0.010, -np.inf, -np.inf, -np.inf),
                (np.inf, mean_high, 0.080, np.inf, np.inf, np.inf),
            )
            pars, cov = curve_fit(
                gaussian_quadratic, x, y, p0=p0, sigma=ey,
                absolute_sigma=True, bounds=bounds, maxfev=100000
            )
            amplitude, mean, sigma, b0, b1, b2 = pars
            model = gaussian_quadratic(x, *pars)
            signal = amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
            background = b0 + b1 * (x - PROTON_MASS_GEV) + b2 * (
                x - PROTON_MASS_GEV
            )**2
            amplitude_error, mean_error, sigma_error = np.sqrt(
                np.clip(np.diag(cov), 0.0, None)
            )[:3]
        elif fixed_sigma is None:
            p0 = (amplitude0, mean0, 0.025, b0, b1)
            bounds = (
                (0.0, mean_low, 0.010, -np.inf, -np.inf),
                (np.inf, mean_high, 0.080, np.inf, np.inf),
            )
            pars, cov = curve_fit(
                gaussian_linear, x, y, p0=p0, sigma=ey,
                absolute_sigma=True, bounds=bounds, maxfev=100000
            )
            amplitude, mean, sigma, b0, b1 = pars
            model = gaussian_linear(x, *pars)
            signal = amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
            background = b0 + b1 * (x - PROTON_MASS_GEV)
            amplitude_error, mean_error, sigma_error = np.sqrt(
                np.clip(np.diag(cov), 0.0, None)
            )[:3]
        elif background_order == 2:
            def model_function(x_values, amplitude, mean, b0, b1, b2):
                return fixed_sigma_quadratic(
                    x_values, amplitude, mean, b0, b1, b2, fixed_sigma
                )
            p0 = (amplitude0, mean0, b0, b1, b2)
            bounds = (
                (0.0, mean_low, -np.inf, -np.inf, -np.inf),
                (np.inf, mean_high, np.inf, np.inf, np.inf),
            )
            pars, cov = curve_fit(
                model_function, x, y, p0=p0, sigma=ey,
                absolute_sigma=True, bounds=bounds, maxfev=100000
            )
            amplitude, mean, b0, b1, b2 = pars
            sigma = float(fixed_sigma)
            model = model_function(x, *pars)
            signal = amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
            background = b0 + b1 * (x - PROTON_MASS_GEV) + b2 * (
                x - PROTON_MASS_GEV
            )**2
            amplitude_error, mean_error = np.sqrt(
                np.clip(np.diag(cov), 0.0, None)
            )[:2]
            sigma_error = 0.0
        else:
            def model_function(x_values, amplitude, mean, b0, b1):
                return fixed_sigma_linear(
                    x_values, amplitude, mean, b0, b1, fixed_sigma
                )
            p0 = (amplitude0, mean0, b0, b1)
            bounds = (
                (0.0, mean_low, -np.inf, -np.inf),
                (np.inf, mean_high, np.inf, np.inf),
            )
            pars, cov = curve_fit(
                model_function, x, y, p0=p0, sigma=ey,
                absolute_sigma=True, bounds=bounds, maxfev=100000
            )
            amplitude, mean, b0, b1 = pars
            sigma = float(fixed_sigma)
            model = model_function(x, *pars)
            signal = amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)
            background = b0 + b1 * (x - PROTON_MASS_GEV)
            amplitude_error, mean_error = np.sqrt(
                np.clip(np.diag(cov), 0.0, None)
            )[:2]
            sigma_error = 0.0
    except Exception as exc:
        result["status"] = f"fit_failed:{type(exc).__name__}"
        return result, diagnostics

    chi2 = float(np.sum(((y - model) / ey) ** 2))
    ndf = int(len(x) - parameter_count)
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan

    in_window = (x >= mean - 2.0 * sigma) & (x <= mean + 2.0 * sigma)
    signal_yield = float(np.sum(np.clip(signal[in_window], 0.0, None)))
    background_yield = float(np.sum(np.clip(background[in_window], 0.0, None)))
    ratio = (
        signal_yield / background_yield
        if background_yield > 0.0 else math.inf
    )

    diagnostics["model"] = model
    diagnostics["signal"] = signal
    diagnostics["background"] = background

    result.update({
        "success": True,
        "status": "success",
        "peak_mean_gev": float(mean),
        "peak_mean_error_gev": float(mean_error),
        "peak_sigma_gev": float(abs(sigma)),
        "peak_sigma_error_gev": float(sigma_error),
        "amplitude": float(amplitude),
        "amplitude_error": float(amplitude_error),
        "chi2_ndf": chi2_ndf,
        "signal_yield_2sigma": signal_yield,
        "background_yield_2sigma": background_yield,
        "signal_to_background_2sigma": ratio,
    })
    return result, diagnostics


def save_fit_overlay(path: Path, diagnostics: dict, fit_result: dict, title: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    counts = diagnostics["counts"]
    edges = diagnostics["edges"]
    centers = diagnostics["centers"]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.step(centers, counts, where="mid", label="Data")
    if fit_result["success"]:
        x = diagnostics["x"]
        ax.plot(x, diagnostics["model"], label="Total fit")
        ax.plot(x, diagnostics["signal"], linestyle="--", label="Elastic component")
        ax.plot(
            x, diagnostics["background"], linestyle=":", label="Background"
        )
        annotation = (
            rf"$\mu_W={fit_result['peak_mean_gev']:.5f}"
            rf"\pm{fit_result['peak_mean_error_gev']:.5f}$ GeV"
            "\n"
            rf"$\sigma_W={fit_result['peak_sigma_gev']:.5f}$ GeV"
            "\n"
            rf"$\chi^2/\mathrm{{ndf}}={fit_result['chi2_ndf']:.2f}$"
            "\n"
            rf"$S/B(\pm2\sigma)={fit_result['signal_to_background_2sigma']:.2f}$"
        )
    else:
        annotation = fit_result["status"]

    ax.text(
        0.03, 0.97, annotation,
        transform=ax.transAxes, ha="left", va="top"
    )
    ax.axvline(PROTON_MASS_GEV, linestyle="--", label="Proton mass")
    ax.set_xlabel(r"$W$ (GeV)")
    ax.set_ylabel("Counts")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def load_period(period: Period, inputs: dict[str, Path], args) -> pd.DataFrame:
    source = f"{inputs[period.source_key]}:{args.tree_name}"
    chunks = []

    for arrays in uproot.iterate(
        source,
        expressions=["runnum", "e_p", "e_theta", "e_phi", "W"],
        library="np",
        step_size=args.step_size,
    ):
        run = np.asarray(arrays["runnum"], dtype=np.int64)
        momentum = np.asarray(arrays["e_p"], dtype=float)
        theta_rad = np.asarray(arrays["e_theta"], dtype=float)
        phi_rad = np.asarray(arrays["e_phi"], dtype=float)
        input_w = np.asarray(arrays["W"], dtype=float)
        beam_energy = beam_energy_from_run(run)

        theta_deg = np.degrees(theta_rad)
        mask = (
            (run >= period.run_min) & (run <= period.run_max)
            & np.isfinite(momentum)
            & np.isfinite(theta_rad)
            & np.isfinite(phi_rad)
            & np.isfinite(input_w)
            & np.isfinite(beam_energy)
            & (momentum > 0.0)
            & (theta_deg >= args.theta_edges[0])
            & (theta_deg < args.theta_edges[-1])
            & (input_w >= args.w_preselection_min)
            & (input_w <= args.w_preselection_max)
        )
        if not np.any(mask):
            continue

        sector = fd_sector_from_phi(phi_rad[mask])
        good = (sector >= 1) & (sector <= 6)
        if not np.any(good):
            continue

        selected_momentum = momentum[mask][good]
        selected_theta = theta_rad[mask][good]
        selected_energy = beam_energy[mask][good]
        recalculated_w = w_from_electron(
            selected_momentum, selected_theta, selected_energy
        )

        chunks.append(pd.DataFrame({
            "runnum": run[mask][good],
            "sector": sector[good],
            "e_p": selected_momentum,
            "e_theta_rad": selected_theta,
            "e_theta_deg": np.degrees(selected_theta),
            "beam_energy_gev": selected_energy,
            "w": recalculated_w,
        }))

    if not chunks:
        raise RuntimeError(f"{period.label}: no selected events")
    return pd.concat(chunks, ignore_index=True)


def variation_definitions(args, integrated_sigma, integrated_sigma_error):
    sigma_error = (
        integrated_sigma_error
        if np.isfinite(integrated_sigma_error) and integrated_sigma_error > 0
        else 0.10 * integrated_sigma
    )
    return [
        ("nominal", 2, (args.w_fit_min, args.w_fit_max), integrated_sigma),
        ("linear_background", 1, (args.w_fit_min, args.w_fit_max), integrated_sigma),
        (
            "narrow_fit_range", 2,
            (args.w_fit_min + args.fit_range_shift, args.w_fit_max - args.fit_range_shift),
            integrated_sigma
        ),
        (
            "wide_fit_range", 2,
            (args.w_fit_min - args.fit_range_shift, args.w_fit_max + args.fit_range_shift),
            integrated_sigma
        ),
        (
            "sigma_minus", 2, (args.w_fit_min, args.w_fit_max),
            max(args.minimum_fixed_sigma, integrated_sigma - sigma_error)
        ),
        (
            "sigma_plus", 2, (args.w_fit_min, args.w_fit_max),
            integrated_sigma + sigma_error
        ),
    ]


def scan_period(period: Period, inputs: dict[str, Path], args) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = load_period(period, inputs, args)
    nominal_records = []
    variation_records = []

    for sector in range(1, 7):
        sector_frame = frame.loc[frame["sector"] == sector].copy()

        integrated_fit, integrated_diag = fit_peak(
            sector_frame["w"].to_numpy(),
            histogram_range=(args.w_hist_min, args.w_hist_max),
            fit_range=(args.w_fit_min, args.w_fit_max),
            histogram_bins=args.w_bins,
            minimum_events=args.minimum_events_integrated,
            fixed_sigma=None,
            mean_center=None,
            mean_half_window=args.mean_half_window,
            background_order=2,
        )

        save_fit_overlay(
            args.output_dir / "plots" / "integrated"
            / f"{period.key}_sector{sector}.png",
            integrated_diag,
            integrated_fit,
            f"{period.label}, sector {sector}, integrated",
        )

        for theta_index, (theta_low, theta_high) in enumerate(
            zip(args.theta_edges[:-1], args.theta_edges[1:])
        ):
            theta_mask = (
                (sector_frame["e_theta_deg"] >= theta_low)
                & (
                    (sector_frame["e_theta_deg"] <= theta_high)
                    if theta_index == len(args.theta_edges) - 2
                    else (sector_frame["e_theta_deg"] < theta_high)
                )
            )
            cell = sector_frame.loc[theta_mask].copy()

            base = {
                "period_key": period.key,
                "period": period.label,
                "sector": sector,
                "theta_bin_index": theta_index,
                "theta_low_deg": float(theta_low),
                "theta_high_deg": float(theta_high),
                "theta_center_deg": float(0.5 * (theta_low + theta_high)),
                "n_events": int(len(cell)),
                "integrated_fit_success": bool(integrated_fit["success"]),
                "integrated_peak_mean_gev": integrated_fit["peak_mean_gev"],
                "integrated_peak_sigma_gev": integrated_fit["peak_sigma_gev"],
                "integrated_peak_sigma_error_gev": integrated_fit[
                    "peak_sigma_error_gev"
                ],
            }

            if not integrated_fit["success"]:
                record = dict(base)
                record.update({
                    "fit_success": False,
                    "fit_status": "integrated_fit_failed",
                })
                nominal_records.append(record)
                continue

            variants = variation_definitions(
                args,
                integrated_fit["peak_sigma_gev"],
                integrated_fit["peak_sigma_error_gev"],
            )
            nominal_fit = None
            nominal_diag = None

            for variation_name, background_order, fit_range, fixed_sigma in variants:
                fit_result, diagnostics = fit_peak(
                    cell["w"].to_numpy(),
                    histogram_range=(args.w_hist_min, args.w_hist_max),
                    fit_range=fit_range,
                    histogram_bins=args.w_bins,
                    minimum_events=args.minimum_events_cell,
                    fixed_sigma=fixed_sigma,
                    mean_center=integrated_fit["peak_mean_gev"],
                    mean_half_window=args.mean_half_window,
                    background_order=background_order,
                )
                variation_record = dict(base)
                variation_record.update({
                    "variation": variation_name,
                    "background_order": background_order,
                    "fit_range_low_gev": fit_range[0],
                    "fit_range_high_gev": fit_range[1],
                    "fixed_sigma_gev": fixed_sigma,
                    **{f"fit_{key}": value for key, value in fit_result.items()},
                })
                variation_records.append(variation_record)

                if variation_name == "nominal":
                    nominal_fit = fit_result
                    nominal_diag = diagnostics

            assert nominal_fit is not None
            record = dict(base)
            record.update({
                "fit_success": nominal_fit["success"],
                "fit_status": nominal_fit["status"],
                "peak_mean_gev": nominal_fit["peak_mean_gev"],
                "peak_mean_error_gev": nominal_fit["peak_mean_error_gev"],
                "peak_sigma_gev": nominal_fit["peak_sigma_gev"],
                "amplitude": nominal_fit["amplitude"],
                "amplitude_error": nominal_fit["amplitude_error"],
                "chi2_ndf": nominal_fit["chi2_ndf"],
                "signal_yield_2sigma": nominal_fit["signal_yield_2sigma"],
                "background_yield_2sigma": nominal_fit["background_yield_2sigma"],
                "signal_to_background_2sigma": nominal_fit[
                    "signal_to_background_2sigma"
                ],
            })

            if nominal_fit["success"]:
                mean_theta = float(cell["e_theta_rad"].mean())
                mean_momentum = float(cell["e_p"].mean())
                mean_energy = float(cell["beam_energy_gev"].mean())

                momentum_at_peak = momentum_at_w(
                    nominal_fit["peak_mean_gev"], mean_theta, mean_energy
                )
                elastic_momentum = elastic_scattered_momentum(
                    mean_theta, mean_energy
                )
                correction = (
                    elastic_momentum - momentum_at_peak
                ) / mean_momentum

                corrected_momentum = cell["e_p"].to_numpy() * (1.0 + correction)
                corrected_w = w_from_electron(
                    corrected_momentum,
                    cell["e_theta_rad"].to_numpy(),
                    cell["beam_energy_gev"].to_numpy(),
                )
                closure_fit, closure_diag = fit_peak(
                    corrected_w,
                    histogram_range=(args.w_hist_min, args.w_hist_max),
                    fit_range=(args.w_fit_min, args.w_fit_max),
                    histogram_bins=args.w_bins,
                    minimum_events=args.minimum_events_cell,
                    fixed_sigma=integrated_fit["peak_sigma_gev"],
                    mean_center=PROTON_MASS_GEV,
                    mean_half_window=args.mean_half_window,
                    background_order=2,
                )
                record.update({
                    "provisional_delta_p_over_p": correction,
                    "closure_fit_success": closure_fit["success"],
                    "closure_peak_mean_gev": closure_fit["peak_mean_gev"],
                    "closure_peak_mean_error_gev": closure_fit[
                        "peak_mean_error_gev"
                    ],
                    "closure_residual_gev": (
                        closure_fit["peak_mean_gev"] - PROTON_MASS_GEV
                        if closure_fit["success"] else math.nan
                    ),
                    "closure_chi2_ndf": closure_fit["chi2_ndf"],
                })

                if theta_low >= args.save_overlays_above:
                    safe_low = str(theta_low).replace(".", "p")
                    safe_high = str(theta_high).replace(".", "p")
                    base_path = (
                        args.output_dir / "plots" / "high_theta_overlays"
                        / period.key / f"sector{sector}"
                    )
                    save_fit_overlay(
                        base_path / f"theta_{safe_low}_{safe_high}_uncorrected.png",
                        nominal_diag,
                        nominal_fit,
                        (
                            f"{period.label}, sector {sector}, "
                            f"{theta_low:g}-{theta_high:g} deg, uncorrected"
                        ),
                    )
                    save_fit_overlay(
                        base_path / f"theta_{safe_low}_{safe_high}_closure.png",
                        closure_diag,
                        closure_fit,
                        (
                            f"{period.label}, sector {sector}, "
                            f"{theta_low:g}-{theta_high:g} deg, provisional closure"
                        ),
                    )
            else:
                record.update({
                    "provisional_delta_p_over_p": math.nan,
                    "closure_fit_success": False,
                    "closure_peak_mean_gev": math.nan,
                    "closure_peak_mean_error_gev": math.nan,
                    "closure_residual_gev": math.nan,
                    "closure_chi2_ndf": math.nan,
                })

            nominal_records.append(record)

    return pd.DataFrame(nominal_records), pd.DataFrame(variation_records)


def evaluate_stability(
    nominal: pd.DataFrame,
    variations: pd.DataFrame,
    args,
) -> pd.DataFrame:
    result = nominal.copy()

    successful_variations = variations.loc[
        variations["fit_success"].fillna(False)
    ].copy()

    grouped = successful_variations.groupby(
        ["period_key", "sector", "theta_bin_index"]
    )
    summary = grouped.agg(
        variation_success_count=("variation", "count"),
        variation_mean_min_gev=("fit_peak_mean_gev", "min"),
        variation_mean_max_gev=("fit_peak_mean_gev", "max"),
        variation_mean_std_gev=("fit_peak_mean_gev", "std"),
    ).reset_index()
    summary["variation_mean_span_gev"] = (
        summary["variation_mean_max_gev"]
        - summary["variation_mean_min_gev"]
    )

    result = result.merge(
        summary,
        how="left",
        on=["period_key", "sector", "theta_bin_index"],
    )
    result["variation_success_count"] = (
        result["variation_success_count"].fillna(0).astype(int)
    )

    amplitude_significance = (
        result["amplitude"] / result["amplitude_error"]
    )
    result["amplitude_significance"] = amplitude_significance

    result["nominal_quality_pass"] = (
        result["fit_success"].fillna(False)
        & (result["n_events"] >= args.minimum_events_cell)
        & (result["peak_mean_error_gev"] <= args.max_mean_error)
        & (amplitude_significance >= args.minimum_significance)
        & (result["chi2_ndf"] <= args.max_chi2_ndf)
        & (
            result["signal_to_background_2sigma"]
            >= args.minimum_signal_to_background
        )
    )

    result["variation_stability_pass"] = (
        (result["variation_success_count"] >= args.minimum_successful_variations)
        & (
            result["variation_mean_span_gev"]
            <= args.max_variation_mean_span
        )
    )

    result["closure_pass"] = (
        result["closure_fit_success"].fillna(False)
        & (
            np.abs(result["closure_residual_gev"])
            <= np.maximum(
                args.absolute_closure_tolerance,
                args.closure_sigma_multiplier
                * result["closure_peak_mean_error_gev"],
            )
        )
        & (result["closure_chi2_ndf"] <= args.max_chi2_ndf)
    )

    result["usable"] = (
        result["nominal_quality_pass"]
        & result["variation_stability_pass"]
        & result["closure_pass"]
    )

    reasons = []
    for row in result.itertuples():
        failed = []
        if not row.nominal_quality_pass:
            failed.append("nominal_quality")
        if not row.variation_stability_pass:
            failed.append("fit_variations")
        if not row.closure_pass:
            failed.append("closure")
        reasons.append("pass" if not failed else ",".join(failed))
    result["usability_status"] = reasons
    return result


def highest_stable_interval(group: pd.DataFrame, args) -> dict:
    ordered = group.sort_values("theta_low_deg").reset_index(drop=True)
    candidates = ordered.loc[
        ordered["theta_low_deg"] >= args.reach_search_min_theta
    ].copy()

    best_low = math.nan
    best_high = math.nan
    best_count = 0

    current_low = None
    current_high = None
    current_count = 0

    for row in candidates.itertuples():
        if bool(row.usable):
            if current_low is None:
                current_low = row.theta_low_deg
                current_count = 1
            else:
                current_count += 1
            current_high = row.theta_high_deg
        else:
            if current_count > 0:
                candidate_is_better = (
                    best_high is None
                    or not np.isfinite(float(best_high))
                    or (
                        current_high is not None
                        and float(current_high) > float(best_high)
                    )
                )
                if candidate_is_better:
                    best_low, best_high, best_count = (
                        float(current_low),
                        float(current_high),
                        int(current_count),
                    )
            current_low = None
            current_high = None
            current_count = 0

    if current_count > 0:
        candidate_is_better = (
            best_high is None
            or not np.isfinite(float(best_high))
            or (
                current_high is not None
                and float(current_high) > float(best_high)
            )
        )
        if candidate_is_better:
            best_low, best_high, best_count = (
                float(current_low),
                float(current_high),
                int(current_count),
            )

    highest_any = (
        float(candidates.loc[candidates["usable"], "theta_high_deg"].max())
        if candidates["usable"].any() else math.nan
    )
    return {
        "highest_stable_interval_low_deg": best_low,
        "highest_stable_interval_high_deg": best_high,
        "stable_interval_bin_count": best_count,
        "highest_any_usable_theta_deg": highest_any,
        "usable_high_theta_bins": int(candidates["usable"].sum()),
        "total_high_theta_bins": int(len(candidates)),
    }


def make_summary_plots(df: pd.DataFrame, output_dir: Path):
    plot_dir = output_dir / "plots" / "summary"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for period, period_df in df.groupby("period", sort=False):
        safe = re.sub(r"[^a-zA-Z0-9]+", "_", period).strip("_").lower()

        metrics = [
            ("n_events", "Selected events", True),
            ("peak_mean_error_gev", r"$\sigma(\mu_W)$ (GeV)", False),
            ("signal_to_background_2sigma", r"$S/B$ within $\pm2\sigma$", False),
            ("variation_mean_span_gev", "Fit-variation centroid span (GeV)", False),
            ("closure_residual_gev", r"Corrected $\mu_W-M_p$ (GeV)", False),
        ]

        for column, ylabel, log_y in metrics:
            fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
            axes = axes.ravel()
            for sector in range(1, 7):
                ax = axes[sector - 1]
                sdf = period_df.loc[period_df["sector"] == sector]
                ax.plot(sdf["theta_center_deg"], sdf[column], marker="o")
                good = sdf.loc[sdf["usable"]]
                ax.scatter(
                    good["theta_center_deg"], good[column],
                    marker="s", label="Usable"
                )
                if log_y:
                    positive = sdf[column] > 0
                    if positive.any():
                        ax.set_yscale("log")
                ax.set_title(f"Sector {sector}")
                ax.set_xlabel(r"$\theta_e$ (deg)")
                ax.set_ylabel(ylabel)
                ax.grid(True, alpha=0.25)
            fig.suptitle(f"{period}: {ylabel}")
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(plot_dir / f"{safe}_{column}.png", dpi=180)
            plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--su22-file", type=Path, default=DEFAULT_INPUTS["su22"])
    parser.add_argument("--fa22-file", type=Path, default=DEFAULT_INPUTS["fa22"])
    parser.add_argument("--sp23-file", type=Path, default=DEFAULT_INPUTS["sp23"])
    parser.add_argument("--tree-name", default="PhysicsEvents")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/electron_theta_reach_scan_v2"),
    )
    parser.add_argument("--step-size", default="250 MB")
    parser.add_argument("--workers", type=int, default=4)

    parser.add_argument(
        "--theta-edges",
        type=parse_edges,
        default=parse_edges(
            "7,8,9,10,11,12,13,14,15,16,17.5,19,21,23,25"
        ),
        help=(
            "Adaptive theta-bin edges in degrees. Default: "
            "7,8,9,10,11,12,13,14,15,16,17.5,19,21,23,25"
        ),
    )
    parser.add_argument(
        "--save-overlays-above",
        type=float,
        default=13.0,
    )
    parser.add_argument(
        "--reach-search-min-theta",
        type=float,
        default=13.0,
    )

    parser.add_argument("--w-preselection-min", type=float, default=0.65)
    parser.add_argument("--w-preselection-max", type=float, default=1.45)
    parser.add_argument("--w-hist-min", type=float, default=0.70)
    parser.add_argument("--w-hist-max", type=float, default=1.35)
    parser.add_argument("--w-fit-min", type=float, default=0.82)
    parser.add_argument("--w-fit-max", type=float, default=1.12)
    parser.add_argument("--w-bins", type=int, default=60)
    parser.add_argument("--fit-range-shift", type=float, default=0.025)
    parser.add_argument("--mean-half-window", type=float, default=0.040)
    parser.add_argument("--minimum-fixed-sigma", type=float, default=0.010)

    parser.add_argument("--minimum-events-integrated", type=int, default=2000)
    parser.add_argument("--minimum-events-cell", type=int, default=500)
    parser.add_argument("--max-mean-error", type=float, default=0.010)
    parser.add_argument("--minimum-significance", type=float, default=3.0)
    parser.add_argument("--max-chi2-ndf", type=float, default=5.0)
    parser.add_argument(
        "--minimum-signal-to-background",
        type=float,
        default=0.10,
    )
    parser.add_argument(
        "--minimum-successful-variations",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--max-variation-mean-span",
        type=float,
        default=0.015,
    )
    parser.add_argument(
        "--absolute-closure-tolerance",
        type=float,
        default=0.005,
    )
    parser.add_argument(
        "--closure-sigma-multiplier",
        type=float,
        default=2.0,
    )
    return parser


def validate_args(args):
    for path in (args.su22_file, args.fa22_file, args.sp23_file):
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.w_fit_min + args.fit_range_shift >= args.w_fit_max - args.fit_range_shift:
        raise ValueError("The narrow fit-range variation is empty.")
    if args.minimum_successful_variations < 1 or args.minimum_successful_variations > 6:
        raise ValueError("--minimum-successful-variations must be in [1, 6].")


def main():
    args = build_parser().parse_args()
    validate_args(args)

    inputs = {
        "su22": args.su22_file,
        "fa22": args.fa22_file,
        "sp23": args.sp23_file,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_dir = args.output_dir / "csv"
    csv_dir.mkdir(parents=True, exist_ok=True)

    nominal_frames = []
    variation_frames = []

    with ProcessPoolExecutor(max_workers=max(1, min(args.workers, 4))) as pool:
        future_map = {
            pool.submit(scan_period, period, inputs, args): period
            for period in PERIODS
        }
        for future in as_completed(future_map):
            period = future_map[future]
            nominal, variations = future.result()
            nominal_frames.append(nominal)
            variation_frames.append(variations)
            print(f"[theta reach v2] Completed {period.label}")

    nominal = pd.concat(nominal_frames, ignore_index=True)
    variations = pd.concat(variation_frames, ignore_index=True)
    evaluated = evaluate_stability(nominal, variations, args)

    evaluated.to_csv(csv_dir / "theta_bin_fit_scan.csv", index=False)
    variations.to_csv(csv_dir / "theta_bin_fit_variations.csv", index=False)

    reach_records = []
    for (period, sector), group in evaluated.groupby(
        ["period", "sector"], sort=False
    ):
        record = {"period": period, "sector": int(sector)}
        record.update(highest_stable_interval(group, args))
        reach_records.append(record)

    reach = pd.DataFrame(reach_records)
    reach.to_csv(csv_dir / "theta_reach_summary.csv", index=False)

    make_summary_plots(evaluated, args.output_dir)

    print("\nStable high-theta reach by period and sector")
    print(reach.to_string(index=False))
    print(f"\nOutputs written under {args.output_dir}")


if __name__ == "__main__":
    main()
