#!/usr/bin/env python3
"""
extract_structure_function_ratios_yield_check.py

Independent yield-level external check for the RGC exclusive e p -> e' n pi+
asymmetry extraction.

This script deliberately does NOT import or apply the production dilution
factor.  Instead it:

  1. divides each of the 24 (xB, -t') bins into nine phi bins;
  2. extracts Gaussian exclusive-signal areas separately for the four NH3
     (beam helicity, target-polarization sign) states;
  3. extracts beam-helicity-separated exclusive-signal areas from C and CH2;
  4. combines Su22, Fa22 and Sp23 into one yield-level data set, obtains one
     common carbon material-normalization coefficient from the broad
     0.00 <= Mx2 < 0.40 GeV2 control region, and subtracts beam-helicity-
     separated carbon rates from each NH3 spin-state rate;
  5. fits the resulting four hydrogen yield/rate distributions simultaneously
     to the five polarized structure-function ratios lu1, ul1, ul2, ll0 and
     ll1.  The unpolarized u1 and u2 modulations are fixed to zero in this
     deliberately simplified external check.

The NH3 Mx2 spectra are fit simultaneously across the four spin states in each
(period, kinematic bin, phi bin) cell.  The peak centroid and width are fixed
from the established spin-integrated nominal Mx2 selection: because the nominal
cut is mu +/- 2sigma, mu=(low+high)/2 and sigma=(high-low)/4.  The four signal
normalizations remain independent while the smooth background shape is shared.

Run from:
    RGC_enpi+/asymmetry_extraction/

Typical command:
    python extract_structure_function_ratios_yield_check_polarized_only_v9.py

Outputs:
    output/asymmetry_extraction/yield_check/

Important statistical note
--------------------------
This first implementation propagates the per-spectrum Gaussian-area covariance
through the beam-helicity-separated carbon subtraction and then performs a weighted
least-squares physics fit.  It is intended as the first external-check
implementation.  The next statistical upgrade should bootstrap the complete
carbon-normalized Gaussian extraction so that shared auxiliary-target correlations
between the four hydrogen spin states can be carried exactly.
"""

from __future__ import annotations

import argparse
import gc
import csv
import json
import math
import os
import sys

# Each Gaussian fit is already a separate process-level task.  Prevent BLAS/OMP
# libraries from creating additional threads inside each worker and accidentally
# oversubscribing an ifarm node.  Users may still override these before launch.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from scipy.optimize import curve_fit, least_squares

# Reuse definitions of the production binning, cuts and depolarization factors,
# but NOT its likelihood or dilution-factor values.
import extract_structure_function_ratios as nominal


PERIODS = nominal.PERIODS
TARGETS = ("NH3", "C", "CH2")
AUX_TARGETS = ("C", "CH2")
CARBON_CONTROL_MIN_GEV2 = 0.0
CARBON_CONTROL_MAX_GEV2 = 0.40
PHYSICS_PARAMETERS = ("lu1", "ul1", "ul2", "ll0", "ll1")
POLARIZED_PARAMETERS = PHYSICS_PARAMETERS
BEAM_POLARIZATION = nominal.BEAM_POLARIZATION
XB_BINS = nominal.XB_BINS
TP_BINS = nominal.MINUS_TPRIME_BINS_GEV2
N_KIN_BINS = nominal.NUMBER_OF_BINS
N_PHI_BINS = 9
MAX_WORKERS = 8
PHI_EDGES = np.linspace(0.0, 2.0 * math.pi, N_PHI_BINS + 1)
PHI_CENTERS = 0.5 * (PHI_EDGES[:-1] + PHI_EDGES[1:])

DEFAULT_OUTPUT_DIR = Path("output/asymmetry_extraction/yield_check")
DEFAULT_RUN_INFO = Path("clas12_run_info.csv")
DEFAULT_CUT_JSON = nominal.DEFAULT_CUT_JSON
DEFAULT_TREE_NAME = "PhysicsEvents"
DEFAULT_CHUNK_SIZE = "250 MB"

PAPER_VERSIONS_DIR = nominal.PAPER_VERSIONS_DIR
DEFAULT_INPUTS = {
    (period, target): PAPER_VERSIONS_DIR
    / f"rgc_{period}_inb_{target}_epi+_mom_corrections.root"
    for period in PERIODS
    for target in TARGETS
}

BRANCH_ALIASES = {
    "runnum": nominal.BRANCH_ALIASES["runnum"],
    "helicity": nominal.BRANCH_ALIASES["helicity"],
    "xB": nominal.BRANCH_ALIASES["xB"],
    "tprime": nominal.BRANCH_ALIASES["tprime"],
    "Mx2": nominal.BRANCH_ALIASES["Mx2"],
    "phi": nominal.BRANCH_ALIASES["phi"],
    "DepA": nominal.BRANCH_ALIASES["DepA"],
    "DepB": nominal.BRANCH_ALIASES["DepB"],
    "DepC": nominal.BRANCH_ALIASES["DepC"],
    "DepV": nominal.BRANCH_ALIASES["DepV"],
    "DepW": nominal.BRANCH_ALIASES["DepW"],
}


@dataclass(frozen=True)
class RunInfo:
    period: str
    target: str
    run: int
    q_plus: float
    q_minus: float
    target_polarization: float


@dataclass
class Spectrum:
    period: str
    target: str
    kin_bin: int
    phi_bin: int
    helicity: int
    target_sign: int
    mx2: np.ndarray
    dep_ratios: dict[str, float]


@dataclass
class AreaFit:
    area: float
    area_error: float
    background_area: float
    chi2_ndf: float
    n_events: int
    status: str
    mu: float
    sigma: float
    fit_low: float
    fit_high: float
    hist_x: np.ndarray
    hist_y: np.ndarray
    hist_err: np.ndarray
    model_y: np.ndarray
    signal_y: np.ndarray
    background_y: np.ndarray


def ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_branch(tree: Any, aliases: Iterable[str]) -> str:
    names = {str(key).split(";")[0]: str(key).split(";")[0] for key in tree.keys()}
    lower = {name.lower(): name for name in names}
    for alias in aliases:
        if alias in names:
            return names[alias]
        # endif
        if alias.lower() in lower:
            return lower[alias.lower()]
        # endif
    # endfor
    raise KeyError(f"Could not resolve any of {tuple(aliases)}; branches={sorted(names)}")


def resolve_tree(root_file: Any, requested: str, path: Path) -> Any:
    try:
        obj = root_file[requested]
        if hasattr(obj, "arrays"):
            return obj
        # endif
    except Exception:
        pass
    # endtry
    candidates = []
    for key in root_file.keys():
        try:
            obj = root_file[key]
            if hasattr(obj, "arrays"):
                candidates.append(obj)
            # endif
        except Exception:
            continue
        # endtry
    # endfor
    if len(candidates) == 1:
        return candidates[0]
    # endif
    raise RuntimeError(f"Could not uniquely resolve tree {requested!r} in {path}")


def parse_run_info(path: Path) -> dict[tuple[str, str, int], RunInfo]:
    period_lookup = {"Su22": "su22", "Fa22": "fa22", "Sp23": "sp23"}
    records: dict[tuple[str, str, int], RunInfo] = {}
    current: tuple[str, str] | None = None
    for line_number, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        # endif
        if line.startswith("#"):
            current = None
            words = line[1:].strip().split()
            if len(words) >= 3 and words[0] == "RGC":
                period = period_lookup.get(words[1])
                target = words[2]
                if period in PERIODS and target in TARGETS:
                    current = (period, target)
                # endif
            # endif
            continue
        # endif
        if current is None:
            continue
        # endif
        fields = [item.strip() for item in line.split(",")]
        if len(fields) < 4:
            raise RuntimeError(f"Malformed run-info row {line_number}: {raw}")
        # endif
        run = int(fields[0])
        q_plus = float(fields[2])
        q_minus = float(fields[3])
        pt = float(fields[4]) if current[1] == "NH3" and len(fields) > 4 else 0.0
        records[(current[0], current[1], run)] = RunInfo(
            current[0], current[1], run, q_plus, q_minus, pt
        )
    # endfor
    return records


def input_override(text: str) -> tuple[tuple[str, str], Path]:
    left, right = text.split("=", 1)
    period, target = left.split(":", 1)
    period = period.lower().strip()
    target_map = {target.lower(): target for target in TARGETS}
    if period not in PERIODS or target.lower() not in target_map:
        raise ValueError(f"Invalid --input {text!r}; use period:target=/path/file.root")
    # endif
    return (period, target_map[target.lower()]), Path(right).expanduser()


def gaussian_density(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    return np.exp(-0.5 * ((x - mu) / sigma) ** 2) / (math.sqrt(2.0 * math.pi) * sigma)


def spectrum_model(
    x: np.ndarray,
    area: float,
    b0: float,
    b1: float,
    b2: float,
    mu: float,
    sigma: float,
    bin_width: float,
) -> np.ndarray:
    z = x - mu
    density = area * gaussian_density(x, mu, sigma) + b0 + b1 * z + b2 * z * z
    return bin_width * density


def fit_signal_area(values: np.ndarray, mu: float, sigma: float) -> AreaFit:
    fit_low = mu - 5.0 * sigma
    fit_high = mu + 5.0 * sigma
    selected = values[np.isfinite(values) & (values >= fit_low) & (values <= fit_high)]
    n_events = int(selected.size)
    n_hist = max(24, min(50, int(round(math.sqrt(max(n_events, 1)) * 1.8))))
    edges = np.linspace(fit_low, fit_high, n_hist + 1)
    y, _ = np.histogram(selected, bins=edges)
    x = 0.5 * (edges[:-1] + edges[1:])
    width = float(edges[1] - edges[0])
    err = np.sqrt(np.maximum(y.astype(float), 1.0))

    if n_events < 20:
        zeros = np.zeros_like(x, dtype=float)
        return AreaFit(np.nan, np.nan, np.nan, np.nan, n_events, "low_statistics",
                       mu, sigma, fit_low, fit_high, x, y.astype(float), err,
                       zeros, zeros, zeros)
    # endif

    side = np.abs(x - mu) > 2.5 * sigma
    b0_guess = float(np.median(y[side]) / width) if np.any(side) else float(np.median(y) / width)
    area_guess = max(float(np.sum(y)) - b0_guess * (fit_high - fit_low), 1.0)
    p0 = np.array([area_guess, max(b0_guess, 0.0), 0.0, 0.0])
    lower = np.array([0.0, 0.0, -np.inf, -np.inf])
    upper = np.array([max(10.0 * n_events, 10.0), np.inf, np.inf, np.inf])

    def model(xx: np.ndarray, area: float, b0: float, b1: float, b2: float) -> np.ndarray:
        return spectrum_model(xx, area, b0, b1, b2, mu, sigma, width)

    try:
        pars, cov = curve_fit(
            model, x, y.astype(float), p0=p0, sigma=err, absolute_sigma=True,
            bounds=(lower, upper), maxfev=30000,
        )
        area_error = math.sqrt(max(float(cov[0, 0]), 0.0)) if np.all(np.isfinite(cov)) else np.nan
        model_y = model(x, *pars)
        signal_y = width * pars[0] * gaussian_density(x, mu, sigma)
        background_y = model_y - signal_y
        chi2 = float(np.sum(((y - model_y) / err) ** 2))
        ndf = max(len(y) - len(pars), 1)
        background_area = float(np.sum(background_y))
        status = "ok" if np.isfinite(area_error) and area_error > 0.0 else "bad_covariance"
        return AreaFit(float(pars[0]), area_error, background_area, chi2 / ndf,
                       n_events, status, mu, sigma, fit_low, fit_high,
                       x, y.astype(float), err, model_y, signal_y, background_y)
    except Exception as exc:
        zeros = np.zeros_like(x, dtype=float)
        return AreaFit(np.nan, np.nan, np.nan, np.nan, n_events,
                       f"fit_failed:{type(exc).__name__}", mu, sigma,
                       fit_low, fit_high, x, y.astype(float), err,
                       zeros, zeros, zeros)
    # endtry


def fit_four_state_signal_areas(
    values_by_state: Mapping[tuple[int, int], np.ndarray], mu: float, sigma: float,
) -> tuple[dict[tuple[int, int], AreaFit], np.ndarray]:
    """Simultaneously fit the four NH3 (h,s) spectra in one phi cell.

    The four Gaussian areas and background normalizations are independent.
    The *shape* of the smooth quadratic background is shared across the four
    states.  This lets the combined statistics constrain the nuisance shape
    without forcing any spin-state signal normalization to agree.  mu and
    sigma remain fixed from the spin-integrated nominal Mx2 selection.

    Returns the four AreaFit objects and their 4x4 Gaussian-area covariance in
    STATE_ORDER.  The latter is propagated into the final hydrogen-rate GLS.
    """
    state_order = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    fit_low, fit_high = mu - 5.0 * sigma, mu + 5.0 * sigma
    selected = {st: np.asarray(values_by_state.get(st, np.empty(0)), dtype=float) for st in state_order}
    selected = {st: v[np.isfinite(v) & (v >= fit_low) & (v <= fit_high)] for st, v in selected.items()}
    total_events = sum(len(v) for v in selected.values())
    n_hist = max(24, min(50, int(round(math.sqrt(max(total_events / 4.0, 1.0)) * 1.8))))
    edges = np.linspace(fit_low, fit_high, n_hist + 1)
    x = 0.5 * (edges[:-1] + edges[1:])
    width = float(edges[1] - edges[0])
    ys = {st: np.histogram(selected[st], bins=edges)[0].astype(float) for st in state_order}
    errs = {st: np.sqrt(np.maximum(ys[st], 1.0)) for st in state_order}

    if total_events < 80 or any(len(selected[st]) < 10 for st in state_order):
        out = {}
        for st in state_order:
            z = np.zeros_like(x)
            out[st] = AreaFit(np.nan, np.nan, np.nan, np.nan, len(selected[st]),
                              "low_statistics_simultaneous", mu, sigma, fit_low, fit_high,
                              x, ys[st], errs[st], z, z, z)
        # endfor
        return out, np.full((4, 4), np.nan)
    # endif

    area0, b00 = [], []
    for st in state_order:
        side = np.abs(x - mu) > 2.5 * sigma
        b = float(np.median(ys[st][side]) / width) if np.any(side) else float(np.median(ys[st]) / width)
        b = max(b, 1.0e-6)
        b00.append(b)
        area0.append(max(float(np.sum(ys[st])) - b * (fit_high - fit_low), 1.0))
    # endfor
    p0 = np.r_[area0, b00, 0.0, 0.0]
    lower = np.r_[np.zeros(4), np.zeros(4), -50.0, -500.0]
    upper = np.r_[np.full(4, max(10.0 * total_events, 10.0)), np.full(4, np.inf), 50.0, 500.0]

    def components(pars: np.ndarray, st_index: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        area = pars[st_index]
        b0 = pars[4 + st_index]
        c1, c2 = pars[8], pars[9]
        z = x - mu
        signal = width * area * gaussian_density(x, mu, sigma)
        background = width * b0 * (1.0 + c1 * z + c2 * z * z)
        return signal + background, signal, background

    def residuals(pars: np.ndarray) -> np.ndarray:
        chunks = []
        for i, st in enumerate(state_order):
            model, _, _ = components(pars, i)
            chunks.append((ys[st] - model) / errs[st])
        # endfor
        return np.concatenate(chunks)

    try:
        result = least_squares(residuals, p0, bounds=(lower, upper), max_nfev=30000)
        try:
            full_cov = np.linalg.inv(result.jac.T @ result.jac)
        except np.linalg.LinAlgError:
            full_cov = np.linalg.pinv(result.jac.T @ result.jac)
        # endtry
        area_cov = full_cov[:4, :4]
        out = {}
        for i, st in enumerate(state_order):
            model, signal, background = components(result.x, i)
            chi2 = float(np.sum(((ys[st] - model) / errs[st]) ** 2))
            ndf = max(len(x) - 3, 1)
            err = math.sqrt(max(float(area_cov[i, i]), 0.0))
            out[st] = AreaFit(float(result.x[i]), err, float(np.sum(background)), chi2 / ndf,
                              len(selected[st]), "ok" if np.isfinite(err) and err > 0 else "bad_covariance",
                              mu, sigma, fit_low, fit_high, x, ys[st], errs[st], model, signal, background)
        # endfor
        return out, area_cov
    except Exception as exc:
        out = {}
        for st in state_order:
            z = np.zeros_like(x)
            out[st] = AreaFit(np.nan, np.nan, np.nan, np.nan, len(selected[st]),
                              f"fit_failed:{type(exc).__name__}", mu, sigma, fit_low, fit_high,
                              x, ys[st], errs[st], z, z, z)
        # endfor
        return out, np.full((4, 4), np.nan)
    # endtry


def fit_task_worker(task: tuple[Any, ...]) -> tuple[Any, Any, Any]:
    """Worker dispatcher: simultaneous NH3 quartet or one auxiliary spectrum."""
    kind = task[0]
    if kind == "NH3":
        _, key, values_by_state, mu, sigma = task
        fits, cov = fit_four_state_signal_areas(values_by_state, mu, sigma)
        return key, fits, cov
    # endif
    _, key, values, mu, sigma = task
    return key, fit_signal_area(values, mu, sigma), None



def carbon_subtracted_hydrogen_rate(
    nh3_area: float, nh3_error: float, nh3_charge: float,
    carbon_area: float, carbon_error: float, carbon_charge: float,
    alpha: float, alpha_error: float,
) -> tuple[float, float]:
    """Return the combined-period carbon-subtracted free-H signal rate.

    One common material-normalization coefficient alpha is used for all three
    periods and both beam helicities.  The carbon *yield* remains explicitly
    beam-helicity separated, so a genuine carbon LU modulation is retained
    without allowing statistical fluctuations in the material scale itself to
    manufacture a beam-spin asymmetry.
    """
    vals = (nh3_area, nh3_error, nh3_charge, carbon_area, carbon_error,
            carbon_charge, alpha, alpha_error)
    if not all(np.isfinite(v) for v in vals):
        return np.nan, np.nan
    # endif
    if nh3_charge <= 0.0 or carbon_charge <= 0.0 or nh3_error <= 0.0 or carbon_error <= 0.0:
        return np.nan, np.nan
    # endif
    r_nh3 = nh3_area / nh3_charge
    r_c = carbon_area / carbon_charge
    rate = r_nh3 - alpha * r_c
    variance = (nh3_error / nh3_charge) ** 2
    variance += (alpha * carbon_error / carbon_charge) ** 2
    variance += (r_c * alpha_error) ** 2
    return rate, math.sqrt(max(variance, 0.0))


def determine_carbon_normalization(
    events: pd.DataFrame,
    charge_map: Mapping[tuple[str, str, int, int], tuple[float, float]],
) -> pd.DataFrame:
    """Determine one common C->NH3 material scale from all periods/helicities."""
    control = events[(events.Mx2 >= CARBON_CONTROL_MIN_GEV2) &
                     (events.Mx2 < CARBON_CONTROL_MAX_GEV2)]
    n_a = int(np.count_nonzero(control.target == "NH3"))
    n_c = int(np.count_nonzero(control.target == "C"))
    q_a = sum(charge_map[(period, "NH3", h, s)][0]
              for period in PERIODS for h in (-1, 1) for s in (-1, 1))
    q_c = sum(charge_map[(period, "C", h, 0)][0]
              for period in PERIODS for h in (-1, 1))
    if n_a > 0 and n_c > 0 and q_a > 0.0 and q_c > 0.0:
        r_a = n_a / q_a
        r_c = n_c / q_c
        alpha = r_a / r_c if r_c > 0.0 else np.nan
        alpha_error = abs(alpha) * math.sqrt(1.0 / n_a + 1.0 / n_c)
    else:
        alpha = np.nan
        alpha_error = np.nan
    # endif
    return pd.DataFrame([dict(
        nh3_control_events=n_a, carbon_control_events=n_c,
        nh3_control_charge=q_a, carbon_control_charge=q_c,
        alpha=alpha, alpha_error=alpha_error,
    )])


def load_events(
    inputs: Mapping[tuple[str, str], Path],
    tree_name: str,
    chunk_size: str,
    run_info: Mapping[tuple[str, str, int], RunInfo],
) -> pd.DataFrame:
    pieces: list[pd.DataFrame] = []
    for period in PERIODS:
        for target in TARGETS:
            path = inputs[(period, target)].expanduser().resolve()
            if not path.is_file():
                raise FileNotFoundError(f"Missing input ROOT file: {path}")
            # endif
            print(f"[input] {period}/{target}: {path}", flush=True)
            with uproot.open(path) as root_file:
                tree = resolve_tree(root_file, tree_name, path)
                branches = {name: resolve_branch(tree, aliases) for name, aliases in BRANCH_ALIASES.items()}
                needed = list(dict.fromkeys(branches.values()))
                for arrays in tree.iterate(needed, step_size=chunk_size, library="np"):
                    run = np.asarray(arrays[branches["runnum"]], dtype=np.int64)
                    hel = np.asarray(arrays[branches["helicity"]], dtype=np.int8)
                    xb = np.asarray(arrays[branches["xB"]], dtype=float)
                    tp_raw = np.asarray(arrays[branches["tprime"]], dtype=float)
                    # Match the production extractor exactly: the stored tprime
                    # branch is signed and the analysis bins in -tprime.
                    minus_tp = -tp_raw
                    mx2 = np.asarray(arrays[branches["Mx2"]], dtype=float)
                    phi_raw = np.asarray(arrays[branches["phi"]], dtype=float)
                    phi = nominal.angle_to_radians(phi_raw, branches["phi"])
                    phi = np.mod(phi, 2.0 * math.pi)
                    _, _, kin = nominal.bin_indices(xb, minus_tp)
                    phi_bin = np.searchsorted(PHI_EDGES, phi, side="right") - 1
                    phi_bin[phi_bin == N_PHI_BINS] = N_PHI_BINS - 1

                    target_sign = np.zeros(run.shape, dtype=np.int8)
                    valid_run = np.ones(run.shape, dtype=bool)
                    for unique_run in np.unique(run):
                        rec = run_info.get((period, target, int(unique_run)))
                        mask = run == unique_run
                        if rec is None:
                            valid_run[mask] = False
                            continue
                        # endif
                        q = rec.q_plus if np.any(hel[mask] == 1) else rec.q_minus
                        if rec.q_plus <= 0.0 and rec.q_minus <= 0.0:
                            valid_run[mask] = False
                        # endif
                        if target == "NH3":
                            target_sign[mask] = 1 if rec.target_polarization > 0.0 else -1
                        # endif
                    # endfor

                    valid = (
                        valid_run & np.isin(hel, (-1, 1)) & (kin >= 1)
                        & (phi_bin >= 0) & (phi_bin < N_PHI_BINS)
                        & np.isfinite(mx2) & np.isfinite(phi)
                    )
                    selected_run = run[valid]
                    event_pt = np.zeros(np.count_nonzero(valid), dtype=float)
                    if target == "NH3":
                        for unique_run in np.unique(selected_run):
                            rec = run_info.get((period, target, int(unique_run)))
                            if rec is not None:
                                event_pt[selected_run == unique_run] = rec.target_polarization
                            # endif
                        # endfor
                    # endif
                    data = {
                        "period": np.full(np.count_nonzero(valid), period),
                        "target": np.full(np.count_nonzero(valid), target),
                        "run": selected_run, "helicity": hel[valid],
                        "target_sign": target_sign[valid], "kin_bin": kin[valid],
                        "phi_bin": phi_bin[valid], "phi": phi[valid], "Mx2": mx2[valid],
                        "event_pt": event_pt,
                        "event_pb": np.full(np.count_nonzero(valid), BEAM_POLARIZATION[period]),
                    }
                    depA = np.asarray(arrays[branches["DepA"]], dtype=float)[valid]
                    for key in ("DepB", "DepC", "DepV", "DepW"):
                        dep = np.asarray(arrays[branches[key]], dtype=float)[valid]
                        data[f"r{key[-1]}"] = np.divide(dep, depA, out=np.full_like(dep, np.nan), where=np.abs(depA) > 0)
                    # endfor
                    pieces.append(pd.DataFrame(data))
                # endfor
            # endwith
    # endfor
    return pd.concat(pieces, ignore_index=True)


def charge_for_state(
    period: str,
    target: str,
    helicity: int,
    target_sign: int,
    run_info: Mapping[tuple[str, str, int], RunInfo],
    loaded_runs: set[int],
) -> tuple[float, float]:
    qsum = 0.0
    qpt = 0.0
    for run in loaded_runs:
        rec = run_info.get((period, target, int(run)))
        if rec is None:
            continue
        # endif
        if target == "NH3" and (1 if rec.target_polarization > 0.0 else -1) != target_sign:
            continue
        # endif
        q = rec.q_plus if helicity == 1 else rec.q_minus
        qsum += q
        if target == "NH3":
            qpt += q * rec.target_polarization
        # endif
    # endfor
    return qsum, qpt


def physics_shape_from_moments(row: Any, theta: np.ndarray) -> float:
    """Polarized-only cell model using averages of the exact event-level products.

    This deliberately avoids factorizing, e.g.
      <Pb Pt rW cos(phi)> -> <Pb Pt><rW>cos(<phi>).
    The moments are measured from NH3 events inside the nominal Mx2 window.
    """
    lu1, ul1, ul2, ll0, ll1 = theta
    h = int(row.helicity)
    return (
        1.0
        + h * row.m_lu1 * lu1
        + row.m_ul1 * ul1
        + row.m_ul2 * ul2
        + h * row.m_ll0 * ll0
        + h * row.m_ll1 * ll1
    )


def accumulate_exact_moments(
    accumulator: dict[tuple[int, int, int, int], dict[str, float]],
    subset: pd.DataFrame | None,
    cut: Any,
    kin_bin: int,
    phi_bin: int,
    helicity: int,
    target_sign: int,
) -> None:
    """Accumulate exact NH3 event-product moments in one pass.

    Only events inside the nominal Mx2 window enter.  Periods are accumulated
    before division, reproducing the previous all-period event-weighted mean
    without repeatedly scanning the full multi-million-row DataFrame.
    """
    if subset is None or len(subset) == 0:
        return
    # endif
    selected = subset[
        (subset.Mx2 >= cut.low_gev2) & (subset.Mx2 <= cut.high_gev2)
    ]
    if len(selected) == 0:
        return
    # endif

    phi = selected.phi.to_numpy(dtype=float, copy=False)
    pb = selected.event_pb.to_numpy(dtype=float, copy=False)
    pt = selected.event_pt.to_numpy(dtype=float, copy=False)
    rB = selected.rB.to_numpy(dtype=float, copy=False)
    rC = selected.rC.to_numpy(dtype=float, copy=False)
    rV = selected.rV.to_numpy(dtype=float, copy=False)
    rW = selected.rW.to_numpy(dtype=float, copy=False)

    values = {
        "m_lu1": pb * rW * np.sin(phi),
        "m_ul1": pt * rV * np.sin(phi),
        "m_ul2": pt * rB * np.sin(2.0 * phi),
        "m_ll0": pb * pt * rC,
        "m_ll1": pb * pt * rW * np.cos(phi),
    }
    key = (kin_bin, phi_bin, helicity, target_sign)
    acc = accumulator.setdefault(
        key,
        {f"{name}_sum": 0.0 for name in values}
        | {f"{name}_count": 0 for name in values},
    )
    for name, value in values.items():
        finite = np.isfinite(value)
        acc[f"{name}_sum"] += float(np.sum(value[finite]))
        acc[f"{name}_count"] += int(np.count_nonzero(finite))
    # endfor


def finalize_exact_moments(
    accumulator: Mapping[tuple[int, int, int, int], Mapping[str, float]],
) -> dict[tuple[int, int, int, int], dict[str, float]]:
    result: dict[tuple[int, int, int, int], dict[str, float]] = {}
    for key, acc in accumulator.items():
        result[key] = {}
        for name in ("m_lu1", "m_ul1", "m_ul2", "m_ll0", "m_ll1"):
            count = int(acc[f"{name}_count"])
            result[key][name] = (
                float(acc[f"{name}_sum"]) / count if count > 0 else np.nan
            )
        # endfor
    # endfor
    return result


def build_rate_covariance(frame: pd.DataFrame, nh3_cov: Mapping[int, np.ndarray]) -> np.ndarray:
    """Full statistical covariance for the carbon-subtracted hydrogen rates.

    Contributions:
      * simultaneous-NH3 Gaussian-area covariance within each phi cell;
      * one shared carbon Gaussian yield for s=+/- at fixed (phi,h);
      * one global alpha nuisance shared by every rate in the kinematic bin.
    """
    rows = list(frame.itertuples(index=False))
    n = len(rows)
    V = np.zeros((n, n), dtype=float)
    state_order = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    state_index = {st: i for i, st in enumerate(state_order)}
    for i, ri in enumerate(rows):
        for j, rj in enumerate(rows):
            cov = 0.0
            # NH3 simultaneous-fit covariance exists only within one phi cell.
            if int(ri.phi_bin) == int(rj.phi_bin):
                C = nh3_cov.get(int(ri.phi_bin))
                if C is not None and np.shape(C) == (4, 4) and ri.charge > 0 and rj.charge > 0:
                    ii = state_index[(int(ri.helicity), int(ri.target_sign))]
                    jj = state_index[(int(rj.helicity), int(rj.target_sign))]
                    if np.isfinite(C[ii, jj]):
                        cov += float(C[ii, jj]) / (ri.charge * rj.charge)
                    # endif
                # endif
            # endif
            # The same helicity-separated carbon area is subtracted from both
            # target-spin states at fixed phi.  It is therefore fully shared.
            if int(ri.phi_bin) == int(rj.phi_bin) and int(ri.helicity) == int(rj.helicity):
                if np.isfinite(ri.carbon_rate_variance):
                    cov += (ri.carbon_alpha ** 2) * ri.carbon_rate_variance
                # endif
            # endif
            # alpha is one common global nuisance, hence correlated across all
            # phi/helicity/target-spin cells.
            if np.isfinite(ri.carbon_rate) and np.isfinite(rj.carbon_rate):
                cov += ri.carbon_rate * rj.carbon_rate * (ri.carbon_alpha_error ** 2)
            # endif
            V[i, j] = cov
        # endfor
    # endfor
    return 0.5 * (V + V.T)


def fit_physics_bin(frame: pd.DataFrame, nh3_cov: Mapping[int, np.ndarray]) -> tuple[dict[str, float], np.ndarray, float, int]:
    good = frame[np.isfinite(frame["hydrogen_rate"])].copy().reset_index(drop=True)
    if len(good) < 18:
        return {name: np.nan for name in PHYSICS_PARAMETERS}, np.full((5, 5), np.nan), np.nan, 0
    # endif

    V = build_rate_covariance(good, nh3_cov)
    finite_diag = np.isfinite(np.diag(V)) & (np.diag(V) > 0.0)
    if not np.all(finite_diag):
        good = good.loc[finite_diag].reset_index(drop=True)
        V = build_rate_covariance(good, nh3_cov)
    # endif
    if len(good) < 18:
        return {name: np.nan for name in PHYSICS_PARAMETERS}, np.full((5, 5), np.nan), np.nan, 0
    # endif

    # Stable Cholesky whitening.  Tiny numerical jitter is allowed only at the
    # 1e-12 scale of the median diagonal; it is not an added physics error.
    scale = float(np.nanmedian(np.diag(V)))
    jitter = max(scale, 1.0e-30) * 1.0e-12
    try:
        L = np.linalg.cholesky(V + jitter * np.eye(len(V)))
    except np.linalg.LinAlgError:
        eigval, eigvec = np.linalg.eigh(V)
        floor = max(float(np.max(eigval)) * 1.0e-12, 1.0e-30)
        V = (eigvec * np.maximum(eigval, floor)) @ eigvec.T
        L = np.linalg.cholesky(V)
    # endtry

    def raw_residuals(pars: np.ndarray) -> np.ndarray:
        theta = pars[:5]
        norm = math.exp(pars[5])
        result = []
        for row in good.itertuples(index=False):
            shape = physics_shape_from_moments(row, theta)
            result.append(row.hydrogen_rate - norm * shape)
        # endfor
        return np.asarray(result, dtype=float)

    def residuals(pars: np.ndarray) -> np.ndarray:
        return np.linalg.solve(L, raw_residuals(pars))

    positive = good.loc[good.hydrogen_rate > 0.0, "hydrogen_rate"]
    norm0 = max(float(np.nanmedian(positive)) if len(positive) else 1.0, 1.0e-12)
    x0 = np.r_[np.zeros(5), math.log(norm0)]
    lower = np.r_[np.full(5, -1.5), -40.0]
    upper = np.r_[np.full(5, 1.5), 40.0]
    result = least_squares(residuals, x0, bounds=(lower, upper), max_nfev=20000)
    resid = residuals(result.x)
    ndf = max(len(resid) - len(result.x), 1)
    chi2_ndf = float(np.dot(resid, resid) / ndf)
    try:
        covariance_full = np.linalg.inv(result.jac.T @ result.jac)
    except np.linalg.LinAlgError:
        covariance_full = np.linalg.pinv(result.jac.T @ result.jac)
    # endtry
    # IMPORTANT: do not multiply by chi2/ndf.  V contains absolute propagated
    # statistical uncertainties; rescaling by goodness-of-fit would incorrectly
    # inflate the polarized errors when nuisance/unpolarized modeling is imperfect.
    covariance = covariance_full[:5, :5]
    values = {name: float(value) for name, value in zip(PHYSICS_PARAMETERS, result.x[:5])}
    return values, covariance, chi2_ndf, len(good)


def plot_spectrum(path: Path, fit: AreaFit, title: str) -> None:
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.errorbar(fit.hist_x, fit.hist_y, yerr=fit.hist_err, fmt="o", ms=3, label="Data")
    ax.plot(fit.hist_x, fit.model_y, label="Gaussian + background")
    ax.plot(fit.hist_x, fit.signal_y, "--", label="Gaussian signal")
    ax.plot(fit.hist_x, fit.background_y, ":", label="Background")
    ax.set_xlabel(r"$M_X^2$ (GeV$^2$)")
    ax.set_ylabel("Counts / bin")
    ax.set_title(title)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_four_state_rates(path: Path, subset: pd.DataFrame, kin_bin: int) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), sharex=True)
    states = [(1, 1), (-1, 1), (1, -1), (-1, -1)]
    for ax, (h, s) in zip(axes.flat, states):
        part = subset[(subset.helicity == h) & (subset.target_sign == s)].sort_values("phi_bin")
        if not part.empty:
            ax.errorbar(np.degrees(part.phi_center), part.hydrogen_rate,
                        yerr=part.hydrogen_rate_error, fmt="o-", ms=4)
        # endif
        ax.set_title(f"h={h:+d}, target sign={s:+d}")
        ax.set_ylabel("Extracted H signal rate (counts / charge)")
        ax.grid(alpha=0.25)
    # endfor
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\phi$ (deg)")
    # endfor
    fig.suptitle(f"Combined RGC direct-yield free-H rates: kinematic bin {kin_bin}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def load_nominal_results() -> pd.DataFrame | None:
    candidates = [
        Path("output/asymmetry_extraction/nominal_zero_uu/tables/structure_function_ratios.csv"),
        Path("output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv"),
    ]
    for path in candidates:
        if path.is_file():
            frame = pd.read_csv(path)
            print(f"[comparison] loaded nominal results: {path}", flush=True)
            return frame
        # endif
    # endfor
    print(f"[comparison] no nominal table found in: {candidates}; skipping overlays", flush=True)
    return None


def plot_by_xb(output_dir: Path, results: pd.DataFrame, nominal_results: pd.DataFrame | None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    nt = len(TP_BINS)
    for x_index, (x_low, x_high) in enumerate(XB_BINS):
        bins = np.arange(x_index * nt + 1, (x_index + 1) * nt + 1)
        subset = results[results.kin_bin.isin(bins)].copy().sort_values("kin_bin")
        subset["t_center"] = [0.5 * (TP_BINS[(int(b)-1) % nt][0] + TP_BINS[(int(b)-1) % nt][1]) for b in subset.kin_bin]
        fig, axes = _grouped_axes()
        for ip, parameter in enumerate(PHYSICS_PARAMETERS):
            ax = axes[parameter]
            ax.errorbar(subset.t_center, subset[parameter], yerr=subset[f"{parameter}_error"],
                        fmt="o", capsize=2, label="Direct-yield check")
            if nominal_results is not None:
                nsub = nominal_results[nominal_results.bin_number.isin(bins)].copy().sort_values("bin_number")
                if len(nsub):
                    if "mean_minus_tprime_gev2" in nsub.columns:
                        nx = nsub.mean_minus_tprime_gev2.to_numpy(float)
                    else:
                        nx = np.asarray([0.5 * (TP_BINS[(int(b)-1) % nt][0] + TP_BINS[(int(b)-1) % nt][1]) for b in nsub.bin_number])
                    errcol = f"{parameter}_stat" if f"{parameter}_stat" in nsub.columns else f"{parameter}_error"
                    ax.errorbar(nx, nsub[parameter], yerr=nsub[errcol], fmt="s", capsize=2, label="Nominal")
                # endif
            # endif
            ax.axhline(0.0, lw=0.8)
            ax.set_ylabel(nominal.PARAMETER_LABELS.get(parameter, parameter))
            ax.grid(alpha=0.25)
            ax.set_xlabel(r"$-t^\prime$ (GeV$^2$)")
        # endfor
        handles, labels = axes["lu1"].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(0.97, 0.04))
        # endif
        fig.suptitle(rf"${x_low:.2f} \leq x_B < {x_high:.2f}$: direct-yield external check", y=0.995)
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        fig.savefig(output_dir / f"yield_check_xB_{x_index+1}.png", dpi=170)
        plt.close(fig)
    # endfor


def write_nominal_comparison(path: Path, results: pd.DataFrame, nominal_results: pd.DataFrame | None) -> None:
    """Write raw nominal-vs-yield differences for the five polarized terms.

    The two extractions use the same underlying events, so their statistical
    errors are correlated.  No quadrature pull or comparison chi2 is formed
    without an explicit cross-method covariance.
    """
    if nominal_results is None:
        return
    # endif
    comp = results.copy()
    n = nominal_results.copy().rename(columns={"bin_number": "kin_bin"})
    keep = ["kin_bin"]
    for p in PHYSICS_PARAMETERS:
        if p in n.columns:
            keep.append(p)
        if f"{p}_stat" in n.columns:
            keep.append(f"{p}_stat")
        elif f"{p}_error" in n.columns:
            keep.append(f"{p}_error")
        # endif
    # endfor
    n = n[keep].rename(columns={c: f"nominal_{c}" for c in keep if c != "kin_bin"})
    comp = comp.merge(n, on="kin_bin", how="left")
    for p in PHYSICS_PARAMETERS:
        if f"nominal_{p}" in comp.columns:
            comp[f"delta_{p}"] = comp[p] - comp[f"nominal_{p}"]
        # endif
    # endfor
    comp.to_csv(path, index=False)

    summary = {
        "included_parameters": list(PHYSICS_PARAMETERS),
        "comparison": "raw nominal-vs-yield differences only",
        "statistical_significance_quoted": False,
        "reason": "The nominal and direct-yield estimators use the same underlying events and therefore have nonzero cross-method statistical covariance.",
        "required_for_pull_or_chi2": "Common-replica/bootstrap estimate of Cov(X_yield, X_nominal).",
    }
    path.with_name("polarized_comparison_summary.json").write_text(json.dumps(summary, indent=2))

def _grouped_axes() -> tuple[plt.Figure, dict[str, plt.Axes]]:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=False)
    mapping = {
        "lu1": axes[0, 0], "ul1": axes[0, 1], "ul2": axes[0, 2],
        "ll0": axes[1, 0], "ll1": axes[1, 1],
    }
    axes[1, 2].axis("off")
    return fig, mapping

def fit_all_exclusive_nh3_bsa(
    frame: pd.DataFrame,
    nh3_cov: Mapping[int, np.ndarray],
) -> tuple[float, float, float, int]:
    """Yield-level BSA from all exclusive NH3 signal events.

    No carbon subtraction and no dilution factor are used.  The other four
    polarized harmonics are floated as nuisance amplitudes so unequal target-
    polarization exposure cannot leak trivially into LU.  Only lu1 is reported.
    """
    good = frame[np.isfinite(frame["nh3_rate"])].copy().reset_index(drop=True)
    if len(good) < 18:
        return np.nan, np.nan, np.nan, 0
    # endif
    rows = list(good.itertuples(index=False))
    V = np.zeros((len(rows), len(rows)), dtype=float)
    state_order = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
    state_index = {st: i for i, st in enumerate(state_order)}
    for i, ri in enumerate(rows):
        for j, rj in enumerate(rows):
            if int(ri.phi_bin) != int(rj.phi_bin):
                continue
            # endif
            C = nh3_cov.get(int(ri.phi_bin))
            if C is None or np.shape(C) != (4, 4):
                continue
            # endif
            ii = state_index[(int(ri.helicity), int(ri.target_sign))]
            jj = state_index[(int(rj.helicity), int(rj.target_sign))]
            V[i, j] = C[ii, jj] / (ri.charge * rj.charge)
        # endfor
    # endfor
    scale = max(float(np.nanmedian(np.diag(V))), 1.0e-30)
    try:
        L = np.linalg.cholesky(V + 1.0e-12 * scale * np.eye(len(V)))
    except np.linalg.LinAlgError:
        eigval, eigvec = np.linalg.eigh(V)
        floor = max(float(np.max(eigval)) * 1.0e-12, 1.0e-30)
        V = (eigvec * np.maximum(eigval, floor)) @ eigvec.T
        L = np.linalg.cholesky(V)
    # endtry

    def residuals(pars: np.ndarray) -> np.ndarray:
        theta = pars[:5]
        norm = math.exp(pars[5])
        raw = np.asarray([
            row.nh3_rate - norm * physics_shape_from_moments(row, theta)
            for row in rows
        ], dtype=float)
        return np.linalg.solve(L, raw)

    positive = good.loc[good.nh3_rate > 0.0, "nh3_rate"]
    norm0 = max(float(np.nanmedian(positive)) if len(positive) else 1.0, 1.0e-12)
    x0 = np.r_[np.zeros(5), math.log(norm0)]
    result = least_squares(
        residuals, x0,
        bounds=(np.r_[np.full(5, -1.5), -40.0], np.r_[np.full(5, 1.5), 40.0]),
        max_nfev=20000,
    )
    resid = residuals(result.x)
    ndf = max(len(resid) - len(result.x), 1)
    try:
        cov = np.linalg.inv(result.jac.T @ result.jac)
    except np.linalg.LinAlgError:
        cov = np.linalg.pinv(result.jac.T @ result.jac)
    # endtry
    return float(result.x[0]), float(math.sqrt(max(cov[0, 0], 0.0))), float(resid @ resid / ndf), len(resid)



def plot_summary(path: Path, results: pd.DataFrame) -> None:
    fig, axes = _grouped_axes()
    for parameter in PHYSICS_PARAMETERS:
        ax = axes[parameter]
        ax.errorbar(results.kin_bin, results[parameter], yerr=results[f"{parameter}_error"], fmt="o", ms=4)
        ax.axhline(0.0, lw=0.8)
        ax.set_title(nominal.PARAMETER_LABELS.get(parameter, parameter))
        ax.set_xlabel("Kinematic bin")
        ax.grid(alpha=0.2)
    # endfor
    fig.suptitle("Direct-yield external extraction: combined Su22 + Fa22 + Sp23")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--run-info", type=Path, default=DEFAULT_RUN_INFO)
    parser.add_argument("--cut-json", type=Path, default=DEFAULT_CUT_JSON)
    parser.add_argument("--tree", default=DEFAULT_TREE_NAME)
    parser.add_argument("--chunk-size", default=DEFAULT_CHUNK_SIZE)
    parser.add_argument("--input", action="append", default=[], help="period:target=/path/file.root")
    parser.add_argument("--workers", type=int, default=min(MAX_WORKERS, os.cpu_count() or 1),
                        help="Gaussian-fit worker processes (default: up to 8; use 1 for serial).")
    parser.add_argument("--max-spectrum-plots", type=int, default=240,
                        help="Maximum individual Mx2 fit plots (0 means all).")
    args = parser.parse_args()
    if args.workers < 1 or args.workers > MAX_WORKERS:
        parser.error(f"--workers must be between 1 and {MAX_WORKERS}")
    # endif

    out = ensure(args.output_dir)
    tables = ensure(out / "tables")
    diagnostics = ensure(out / "diagnostics")
    spectrum_plots = ensure(diagnostics / "mx2_fits")
    rate_plots = ensure(diagnostics / "hydrogen_rates")
    physics_plots = ensure(out / "physics")

    inputs = dict(DEFAULT_INPUTS)
    for text in args.input:
        key, path = input_override(text)
        inputs[key] = path
    # endfor

    run_info = parse_run_info(args.run_info)
    cuts = nominal.load_channel_cuts(args.cut_json, "nominal")
    events = load_events(inputs, args.tree, args.chunk_size, run_info)
    print(f"[input] selected binned events loaded: {len(events):,}", flush=True)
    print(
        "[performance] exact harmonic moments will be accumulated during spectrum grouping; "
        "the full event DataFrame will be released before Gaussian-fit workers start.",
        flush=True,
    )

    loaded_runs: dict[tuple[str, str], set[int]] = {
        (period, target): set(events[(events.period == period) & (events.target == target)].run.unique().astype(int))
        for period in PERIODS for target in TARGETS
    }

    # Charges are state-specific and use only runs actually present in each ROOT file.
    charge_rows = []
    charge_map: dict[tuple[str, str, int, int], tuple[float, float]] = {}
    for period in PERIODS:
        for target in TARGETS:
            signs = (-1, 1) if target == "NH3" else (0,)
            for h in (-1, 1):
                for s in signs:
                    q, qpt = charge_for_state(period, target, h, s, run_info, loaded_runs[(period, target)])
                    charge_map[(period, target, h, s)] = (q, qpt)
                    charge_rows.append(dict(period=period, target=target, helicity=h,
                                            target_sign=s, charge=q, charge_times_pt=qpt))
                # endfor
            # endfor
        # endfor
    # endfor
    pd.DataFrame(charge_rows).to_csv(tables / "state_charges.csv", index=False)

    carbon_norms = determine_carbon_normalization(events, charge_map)
    carbon_norms.to_csv(tables / "carbon_normalization.csv", index=False)
    alpha = float(carbon_norms.iloc[0].alpha)
    alpha_error = float(carbon_norms.iloc[0].alpha_error)
    print("[carbon normalization] one common all-period/all-helicity scale:", flush=True)
    print(carbon_norms.to_string(index=False), flush=True)

    area_rows: list[dict[str, Any]] = []
    area_lookup: dict[tuple[str, str, int, int, int, int], AreaFit] = {}
    plot_count = 0

    # Build fit inputs in the parent process.  NH3 is fit as one simultaneous
    # four-state (h,s) quartet per period/kinematic/phi cell.  C and CH2 remain
    # independent helicity-separated spectra.
    fit_tasks: list[tuple[Any, ...]] = []
    fit_metadata: dict[Any, dict[str, Any]] = {}
    exact_moment_accumulator: dict[tuple[int, int, int, int], dict[str, float]] = {}
    grouped = {
        key: group for key, group in events.groupby(
            ["period", "target", "kin_bin", "phi_bin", "helicity", "target_sign"],
            sort=False, observed=True,
        )
    }

    for period in PERIODS:
        for kin_bin in range(1, N_KIN_BINS + 1):
            cut = cuts[(period, kin_bin)]
            mu = 0.5 * (cut.low_gev2 + cut.high_gev2)
            sigma = 0.25 * (cut.high_gev2 - cut.low_gev2)
            for phi_bin in range(N_PHI_BINS):
                # Simultaneous NH3 quartet.
                values_by_state = {}
                for h in (-1, 1):
                    for ss in (-1, 1):
                        key = (period, "NH3", kin_bin, phi_bin, h, ss)
                        subset = grouped.get(key)
                        values_by_state[(h, ss)] = (subset.Mx2.to_numpy(dtype=float, copy=True)
                                                     if subset is not None else np.empty(0, dtype=float))
                        dep = ({name: float(np.nanmean(subset[name])) if len(subset) else np.nan
                                for name in ("rB", "rC", "rV", "rW")} if subset is not None
                               else {name: np.nan for name in ("rB", "rC", "rV", "rW")})
                        fit_metadata[key] = dict(period=period, target="NH3", kin_bin=kin_bin,
                                                 phi_bin=phi_bin, phi_center=PHI_CENTERS[phi_bin],
                                                 helicity=h, target_sign=ss, mu=mu, sigma=sigma, **dep)
                        accumulate_exact_moments(
                            exact_moment_accumulator, subset, cut,
                            kin_bin, phi_bin, h, ss,
                        )
                    # endfor
                # endfor
                fit_tasks.append(("NH3", (period, kin_bin, phi_bin), values_by_state, mu, sigma))

                # Auxiliary spectra stay helicity separated.
                for target in ("C", "CH2"):
                    for h in (-1, 1):
                        key = (period, target, kin_bin, phi_bin, h, 0)
                        subset = grouped.get(key)
                        values = subset.Mx2.to_numpy(dtype=float, copy=True) if subset is not None else np.empty(0, dtype=float)
                        dep = ({name: float(np.nanmean(subset[name])) if len(subset) else np.nan
                                for name in ("rB", "rC", "rV", "rW")} if subset is not None
                               else {name: np.nan for name in ("rB", "rC", "rV", "rW")})
                        fit_metadata[key] = dict(period=period, target=target, kin_bin=kin_bin,
                                                 phi_bin=phi_bin, phi_center=PHI_CENTERS[phi_bin],
                                                 helicity=h, target_sign=0, mu=mu, sigma=sigma, **dep)
                        fit_tasks.append(("AUX", key, values, mu, sigma))
                    # endfor
                # endfor
            # endfor
        # endfor
    # endfor

    exact_moment_lookup = finalize_exact_moments(exact_moment_accumulator)
    del exact_moment_accumulator
    del grouped
    del events
    gc.collect()

    print(
        f"[fits] fit tasks: {len(fit_tasks):,} (NH3 quartets + auxiliary spectra); "
        f"workers: {args.workers}",
        flush=True,
    )

    nh3_area_cov_period: dict[tuple[str, int, int], np.ndarray] = {}
    completed_fits = 0
    progress_step = max(1, len(fit_tasks) // 20)

    def consume_fit_result(result: tuple[Any, Any, np.ndarray]) -> None:
        nonlocal plot_count, completed_fits
        task_key, payload, cov = result
        completed_fits += 1
        if isinstance(payload, dict):
            period, kin_bin, phi_bin = task_key
            nh3_area_cov_period[(period, kin_bin, phi_bin)] = cov
            for (h, ss), fit in payload.items():
                key = (period, "NH3", kin_bin, phi_bin, h, ss)
                area_lookup[key] = fit
                meta = fit_metadata[key]
                area_rows.append(dict(**meta, area=fit.area, area_error=fit.area_error,
                                      background_area=fit.background_area, chi2_ndf=fit.chi2_ndf,
                                      n_events=fit.n_events, status=fit.status))
                if (args.max_spectrum_plots == 0 or plot_count < args.max_spectrum_plots) and fit.status == "ok":
                    plot_spectrum(
                        spectrum_plots / f"{period}_NH3_bin{kin_bin:02d}_phi{phi_bin:02d}_h{h:+d}_s{ss:+d}.png",
                        fit,
                        f"{period} NH3; bin {kin_bin}; phi bin {phi_bin}; h={h:+d}; s={ss:+d}",
                    )
                    plot_count += 1
                # endif
            # endfor
        else:
            key, fit = task_key, payload
            area_lookup[key] = fit
            meta = fit_metadata[key]
            area_rows.append(dict(**meta, area=fit.area, area_error=fit.area_error,
                                  background_area=fit.background_area, chi2_ndf=fit.chi2_ndf,
                                  n_events=fit.n_events, status=fit.status))
            if (args.max_spectrum_plots == 0 or plot_count < args.max_spectrum_plots) and fit.status == "ok":
                period, target, kin_bin, phi_bin, h, ss = key
                plot_spectrum(
                    spectrum_plots / f"{period}_{target}_bin{kin_bin:02d}_phi{phi_bin:02d}_h{h:+d}_s{ss:+d}.png",
                    fit,
                    f"{period} {target}; bin {kin_bin}; phi bin {phi_bin}; h={h:+d}",
                )
                plot_count += 1
            # endif
        # endif
        if completed_fits % progress_step == 0 or completed_fits == len(fit_tasks):
            print(
                f"[fits] progress: {completed_fits:,}/{len(fit_tasks):,} "
                f"({100.0 * completed_fits / len(fit_tasks):.1f}%)",
                flush=True,
            )
        # endif

    if args.workers == 1:
        for task in fit_tasks:
            consume_fit_result(fit_task_worker(task))
        # endfor
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(fit_task_worker, task) for task in fit_tasks]
            for future in as_completed(futures):
                consume_fit_result(future.result())
            # endfor
        # endwith
    # endif


    del fit_tasks
    gc.collect()
    area_frame = pd.DataFrame(area_rows)
    area_frame.to_csv(tables / "gaussian_signal_areas.csv", index=False)

    # Combine Su22, Fa22 and Sp23 before the physics extraction.  Gaussian
    # areas remain period-specific because their fixed Mx2 peak shapes/cuts are
    # period-specific, but areas, charges and polarization moments are summed
    # here to form one RGC data set.
    rate_rows = []
    # Combined-period NH3 area covariance, indexed by (kin_bin, phi_bin).
    # Periods are statistically independent, so their covariance matrices add.
    nh3_cov_combined: dict[tuple[int, int], np.ndarray] = {}
    for kin_bin in range(1, N_KIN_BINS + 1):
        for phi_bin in range(N_PHI_BINS):
            mats = [nh3_area_cov_period.get((period, kin_bin, phi_bin)) for period in PERIODS]
            mats = [m for m in mats if m is not None and np.shape(m) == (4, 4) and np.all(np.isfinite(m))]
            nh3_cov_combined[(kin_bin, phi_bin)] = (sum(mats, np.zeros((4, 4))) if mats
                                                        else np.full((4, 4), np.nan))
        # endfor
    # endfor

    for kin_bin in range(1, N_KIN_BINS + 1):
        for phi_bin in range(N_PHI_BINS):
            for h in (-1, 1):
                cfits = [area_lookup[(period, "C", kin_bin, phi_bin, h, 0)] for period in PERIODS]
                cfinite = [f for f in cfits if np.isfinite(f.area) and np.isfinite(f.area_error)]
                carbon_area = sum(f.area for f in cfinite) if cfinite else np.nan
                carbon_error = math.sqrt(sum(f.area_error ** 2 for f in cfinite)) if cfinite else np.nan
                qC = sum(charge_map[(period, "C", h, 0)][0] for period in PERIODS)
                for s in (-1, 1):
                    nfits = [area_lookup[(period, "NH3", kin_bin, phi_bin, h, s)] for period in PERIODS]
                    nfinite = [f for f in nfits if np.isfinite(f.area) and np.isfinite(f.area_error)]
                    nh3_area = sum(f.area for f in nfinite) if nfinite else np.nan
                    nh3_error = math.sqrt(sum(f.area_error ** 2 for f in nfinite)) if nfinite else np.nan
                    qA = sum(charge_map[(period, "NH3", h, s)][0] for period in PERIODS)
                    qpt = sum(charge_map[(period, "NH3", h, s)][1] for period in PERIODS)
                    qpb = sum(charge_map[(period, "NH3", h, s)][0] * BEAM_POLARIZATION[period] for period in PERIODS)
                    qpbpt = sum(charge_map[(period, "NH3", h, s)][1] * BEAM_POLARIZATION[period] for period in PERIODS)
                    h_rate, _ = carbon_subtracted_hydrogen_rate(
                        nh3_area, nh3_error, qA, carbon_area, carbon_error, qC, alpha, alpha_error
                    )
                    carbon_rate = carbon_area / qC if np.isfinite(carbon_area) and qC > 0 else np.nan
                    carbon_rate_variance = ((carbon_error / qC) ** 2
                                            if np.isfinite(carbon_error) and qC > 0 else np.nan)
                    h_error = np.nan  # filled from the full covariance after the table is assembled
                    pt_eff = qpt / qA if qA > 0.0 else np.nan
                    pb_eff = qpb / qA if qA > 0.0 else np.nan
                    pbpt_eff = qpbpt / qA if qA > 0.0 else np.nan
                    nh3_rows = area_frame[
                        (area_frame.target == "NH3") & (area_frame.kin_bin == kin_bin)
                        & (area_frame.phi_bin == phi_bin) & (area_frame.helicity == h)
                        & (area_frame.target_sign == s)
                    ]
                    # Event-weighted depolarization ratios over all three periods.
                    dep = {}
                    for name in ("rB", "rC", "rV", "rW"):
                        vals = nh3_rows[name].to_numpy(dtype=float)
                        weights = np.maximum(nh3_rows.n_events.to_numpy(dtype=float), 1.0)
                        finite_dep = np.isfinite(vals) & np.isfinite(weights)
                        dep[name] = (float(np.average(vals[finite_dep], weights=weights[finite_dep]))
                                     if np.any(finite_dep) else np.nan)
                    # endfor
                    moments = exact_moment_lookup.get(
                        (kin_bin, phi_bin, h, s),
                        {name: np.nan for name in ("m_lu1", "m_ul1", "m_ul2", "m_ll0", "m_ll1")},
                    )
                    rate_rows.append(dict(
                        kin_bin=kin_bin, phi_bin=phi_bin, phi_center=PHI_CENTERS[phi_bin],
                        helicity=h, target_sign=s, charge=qA, pt_effective=pt_eff,
                        pb_effective=pb_eff, pbpt_effective=pbpt_eff, hydrogen_rate=h_rate,
                        hydrogen_rate_error=h_error, nh3_signal_area=nh3_area,
                        nh3_signal_area_error=nh3_error, carbon_signal_area=carbon_area,
                        carbon_signal_area_error=carbon_error, carbon_rate=carbon_rate,
                        carbon_rate_variance=carbon_rate_variance, carbon_alpha=alpha,
                        carbon_alpha_error=alpha_error, **dep, **moments,
                    ))
                # endfor
            # endfor
        # endfor
    # endfor

    rates = pd.DataFrame(rate_rows)
    # Populate display/CSV point errors from the diagonal of the same full
    # covariance matrix used by the generalized least-squares physics fit.
    for kin_bin in range(1, N_KIN_BINS + 1):
        idx = rates.index[rates.kin_bin == kin_bin]
        sub = rates.loc[idx].reset_index(drop=True)
        phi_cov = {pb: nh3_cov_combined[(kin_bin, pb)] for pb in range(N_PHI_BINS)}
        V = build_rate_covariance(sub, phi_cov)
        rates.loc[idx, "hydrogen_rate_error"] = np.sqrt(np.maximum(np.diag(V), 0.0))
    # endfor
    rates["nh3_rate"] = np.divide(
        rates["nh3_signal_area"].to_numpy(dtype=float),
        rates["charge"].to_numpy(dtype=float),
        out=np.full(len(rates), np.nan, dtype=float),
        where=rates["charge"].to_numpy(dtype=float) > 0.0,
    )
    rates.to_csv(tables / "carbon_subtracted_hydrogen_rates.csv", index=False)

    # Apples-to-apples BSA check: use every exclusive NH3 signal event, with no
    # carbon subtraction and no dilution factor.  This is the yield-level
    # analogue of the nominal undiluted LU term.
    all_exclusive_bsa_rows = []
    for kin_bin in range(1, N_KIN_BINS + 1):
        sub = rates[rates.kin_bin == kin_bin].copy()
        phi_cov = {pb: nh3_cov_combined[(kin_bin, pb)] for pb in range(N_PHI_BINS)}
        lu, lu_err, chi2_ndf, npoints = fit_all_exclusive_nh3_bsa(sub, phi_cov)
        all_exclusive_bsa_rows.append({
            "kin_bin": kin_bin, "lu1_all_exclusive_nh3": lu,
            "lu1_all_exclusive_nh3_error": lu_err,
            "chi2_ndf": chi2_ndf, "npoints": npoints,
        })
    # endfor
    pd.DataFrame(all_exclusive_bsa_rows).to_csv(
        tables / "all_exclusive_nh3_bsa.csv", index=False
    )

    result_rows = []
    covariance_payload: dict[str, Any] = {}
    for kin_bin in range(1, N_KIN_BINS + 1):
        subset = rates[rates.kin_bin == kin_bin]
        plot_four_state_rates(rate_plots / f"hydrogen_rates_bin{kin_bin:02d}.png", subset, kin_bin)
        phi_cov = {pb: nh3_cov_combined[(kin_bin, pb)] for pb in range(N_PHI_BINS)}
        values, covariance, chi2_ndf, npoints = fit_physics_bin(subset, phi_cov)
        row: dict[str, Any] = {"kin_bin": kin_bin, "fit_chi2_ndf_polarized_model": chi2_ndf,
                               "chi2_ndf": chi2_ndf, "npoints": npoints}
        for index, name in enumerate(PHYSICS_PARAMETERS):
            row[name] = values[name]
            row[f"{name}_error"] = math.sqrt(max(float(covariance[index, index]), 0.0)) if np.isfinite(covariance[index, index]) else np.nan
        # endfor
        result_rows.append(row)
        covariance_payload[str(kin_bin)] = covariance.tolist()
    # endfor

    results = pd.DataFrame(result_rows)
    results.to_csv(tables / "yield_check_structure_function_ratios.csv", index=False)
    (tables / "yield_check_covariances.json").write_text(json.dumps(covariance_payload, indent=2))
    plot_summary(physics_plots / "yield_check_structure_function_ratios.png", results)
    nominal_results = load_nominal_results()
    plot_by_xb(physics_plots / "by_xB", results, nominal_results)
    write_nominal_comparison(tables / "yield_check_vs_nominal.csv", results, nominal_results)

    # Useful stability/quality figures.
    fig, ax = plt.subplots(figsize=(10, 6))
    ok = area_frame[np.isfinite(area_frame.chi2_ndf)]
    ax.hist(ok.chi2_ndf, bins=np.linspace(0, min(10, max(2, float(ok.chi2_ndf.quantile(0.99)))), 60))
    ax.set_xlabel(r"Gaussian-area fit $\chi^2$/ndf")
    ax.set_ylabel("Spectra")
    ax.set_title("Mx2 fit quality")
    fig.tight_layout(); fig.savefig(diagnostics / "mx2_fit_chi2_distribution.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6))
    finite = rates[np.isfinite(rates.hydrogen_rate) & np.isfinite(rates.hydrogen_rate_error)]
    pullscale = np.abs(finite.hydrogen_rate) / finite.hydrogen_rate_error
    ax.hist(pullscale, bins=60)
    ax.set_xlabel(r"$|R_H|/\delta R_H$")
    ax.set_ylabel("Spin/phi cells")
    ax.set_title("Carbon-subtracted hydrogen-rate statistical significance")
    fig.tight_layout(); fig.savefig(diagnostics / "hydrogen_rate_significance.png", dpi=160); plt.close(fig)

    manifest = {
        "method": "combined-period direct-yield carbon-subtracted hydrogen-rate external check",
        "n_phi_bins": N_PHI_BINS,
        "phi_edges_rad": PHI_EDGES.tolist(),
        "does_not_use_dilution_factor": True,
        "auxiliary_targets_beam_helicity_separated": True,
        "periods_combined_before_physics_fit": True,
        "carbon_material_scale_common_to_periods_and_helicities": True,
        "nh3_beam_and_target_spin_separated": True,
        "signal_shape": "mu/sigma fixed from nominal mu +/- 2 sigma channel-selection cut",
        "statistics_treatment": "Full GLS covariance with exact event-product harmonic moments inside each phi cell; simultaneous-NH3 area covariance + shared carbon-yield covariance + global alpha nuisance covariance; no chi2/ndf rescaling.",
        "nh3_mx2_fit": "simultaneous four-state fit with independent signal/background normalizations and shared quadratic background shape",
        "fit_and_comparison_parameters": list(POLARIZED_PARAMETERS),
        "unpolarized_parameters": "u1 and u2 are fixed to zero and are not fit in this polarized-only external check.",
        "v7_correlation_anecdote": "In the preceding seven-parameter study, mean correlations of u1/u2 with the polarized coefficients were generally only at the few-percent level; this supports using the simplified polarized-only fit as an external check.",
        "nominal_comparison_statistics": "Raw differences only; no pull or comparison chi2 is quoted because the two methods share underlying events and their cross-method covariance has not been estimated.",
        "outputs": {
            "areas": str(tables / "gaussian_signal_areas.csv"),
            "rates": str(tables / "carbon_subtracted_hydrogen_rates.csv"),
            "carbon_normalization": str(tables / "carbon_normalization.csv"),
            "physics": str(tables / "yield_check_structure_function_ratios.csv"),
            "all_exclusive_bsa": str(tables / "all_exclusive_nh3_bsa.csv"),
        },
    }
    (out / "yield_check_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[done] yield-check outputs written to {out.resolve()}")


if __name__ == "__main__":
    main()
# endif
