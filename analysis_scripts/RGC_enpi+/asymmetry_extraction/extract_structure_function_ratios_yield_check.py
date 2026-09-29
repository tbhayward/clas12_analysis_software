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
  3. extracts beam-helicity-separated exclusive-signal areas from C, CH2, He
     and empty-target data;
  4. evaluates the direct five-target Method-1 algebra as a FREE-HYDROGEN RATE
     rather than first constructing f = N_H/N_NH3;
  5. fits the resulting four hydrogen yield/rate distributions simultaneously
     to the same seven longitudinal structure-function ratios used by the
     nominal analysis.

The first implementation intentionally uses independent Gaussian-area fits in
each state.  The peak centroid and width are fixed from the established
spin-integrated nominal Mx2 selection: because the nominal cut is mu +/- 2sigma,
mu=(low+high)/2 and sigma=(high-low)/4.  A later version can replace these
independent fits with one simultaneous Mx2 fit sharing signal/background shape
parameters across spin states.

Run from:
    RGC_enpi+/asymmetry_extraction/

Typical command:
    python extract_structure_function_ratios_yield_check_parallel_v3.py

Outputs:
    output/asymmetry_extraction/yield_check/

Important statistical note
--------------------------
This first implementation propagates the per-spectrum Gaussian-area covariance
through the direct five-target expression and then performs a weighted
least-squares physics fit.  It is intended as the first external-check
implementation.  The next statistical upgrade should bootstrap the complete
five-target Gaussian extraction so that shared auxiliary-target correlations
between the four hydrogen spin states are carried exactly.
"""

from __future__ import annotations

import argparse
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
from concurrent.futures import ProcessPoolExecutor
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
TARGETS = ("NH3", "C", "CH2", "He", "ET")
AUX_TARGETS = ("C", "CH2", "He", "ET")
PHYSICS_PARAMETERS = nominal.PHYSICS_PARAMETERS
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


def fit_signal_area_worker(task: tuple[Any, np.ndarray, float, float]) -> tuple[Any, AreaFit]:
    """Process-pool worker for one independent Mx2 Gaussian-area fit.

    Only the compact Mx2 array and fixed peak parameters are passed to the
    worker.  No pandas frame, ROOT handle or Matplotlib object crosses the
    process boundary.
    """
    key, values, mu, sigma = task
    return key, fit_signal_area(values, mu, sigma)



def direct_method1_hydrogen_rate(
    areas: Mapping[str, float],
    charges: Mapping[str, float],
) -> float:
    """Return the free-H signal rate directly from Method-1 Eq. (10).

    Eq. (10) returns f.  Algebraically multiplying Eq. (10) by N_A/Q_A and
    cancelling N_A gives this expression.  Thus this routine never constructs
    or consumes a dilution factor.

    Charges can be absolute charges rather than fractions: the common charge
    normalization cancels from the expression.
    """
    nA, nC = float(areas["NH3"]), float(areas["C"])
    nCH, nMT, nf = float(areas["CH2"]), float(areas["He"]), float(areas["ET"])
    qA, qC = float(charges["NH3"]), float(charges["C"])
    qCH, qHe, qf = float(charges["CH2"]), float(charges["He"]), float(charges["ET"])
    if min(qA, qC, qCH, qHe, qf) <= 0.0:
        return np.nan
    # endif
    first = -nMT * qA + nA * qHe
    second = (
        -0.579353 * nMT * qC * qCH * qf
        + (nf * qC * qCH - 3.50431 * nCH * qC * qf + 3.08366 * nC * qCH * qf) * qHe
    )
    bracket = (
        35.88 * nMT * qC * qCH * qf
        - nf * qC * qCH * qHe
        - 43.3586 * nCH * qC * qf * qHe
        + 8.47866 * nC * qCH * qf * qHe
    )
    denominator = qA * qHe * bracket
    if denominator == 0.0:
        return np.nan
    # endif
    return 12.3729 * first * second / denominator


def propagate_direct_rate_error(
    areas: Mapping[str, float],
    errors: Mapping[str, float],
    charges: Mapping[str, float],
) -> float:
    variance = 0.0
    for target in TARGETS:
        error = float(errors[target])
        if not np.isfinite(error) or error <= 0.0:
            return np.nan
        # endif
        value = float(areas[target])
        step = max(1.0e-4 * max(abs(value), 1.0), 1.0e-3 * error)
        plus = dict(areas)
        minus = dict(areas)
        plus[target] = value + step
        minus[target] = max(value - step, 0.0)
        actual_step = plus[target] - minus[target]
        if actual_step <= 0.0:
            return np.nan
        # endif
        derivative = (
            direct_method1_hydrogen_rate(plus, charges)
            - direct_method1_hydrogen_rate(minus, charges)
        ) / actual_step
        variance += derivative * derivative * error * error
    # endfor
    return math.sqrt(max(variance, 0.0))


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
                    data = {
                        "period": np.full(np.count_nonzero(valid), period),
                        "target": np.full(np.count_nonzero(valid), target),
                        "run": run[valid], "helicity": hel[valid],
                        "target_sign": target_sign[valid], "kin_bin": kin[valid],
                        "phi_bin": phi_bin[valid], "phi": phi[valid], "Mx2": mx2[valid],
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


def physics_shape(phi: float, h: int, pb: float, pt_eff: float,
                  rB: float, rC: float, rV: float, rW: float,
                  theta: np.ndarray) -> float:
    u1, u2, lu1, ul1, ul2, ll0, ll1 = theta
    return (
        1.0 + rV * u1 * math.cos(phi) + rB * u2 * math.cos(2.0 * phi)
        + h * pb * rW * lu1 * math.sin(phi)
        + pt_eff * (rV * ul1 * math.sin(phi) + rB * ul2 * math.sin(2.0 * phi))
        + h * pb * pt_eff * (rC * ll0 + rW * ll1 * math.cos(phi))
    )


def fit_physics_bin(frame: pd.DataFrame) -> tuple[dict[str, float], np.ndarray, float, int]:
    good = frame[
        np.isfinite(frame["hydrogen_rate"]) & np.isfinite(frame["hydrogen_rate_error"])
        & (frame["hydrogen_rate_error"] > 0.0)
    ].copy()
    if len(good) < 18:
        return {name: np.nan for name in PHYSICS_PARAMETERS}, np.full((7, 7), np.nan), np.nan, 0
    # endif

    # One free normalization per period absorbs luminosity-independent acceptance
    # and the overall unpolarized rate.  Phi-dependent acceptance is assumed to
    # cancel between spin states, as in the nominal extraction.
    period_indices = {period: index for index, period in enumerate(PERIODS)}

    def residuals(pars: np.ndarray) -> np.ndarray:
        theta = pars[:7]
        norms = np.exp(pars[7:])
        result = []
        for row in good.itertuples(index=False):
            shape = physics_shape(
                row.phi_center, int(row.helicity), BEAM_POLARIZATION[row.period],
                row.pt_effective, row.rB, row.rC, row.rV, row.rW, theta,
            )
            prediction = norms[period_indices[row.period]] * shape
            result.append((row.hydrogen_rate - prediction) / row.hydrogen_rate_error)
        # endfor
        return np.asarray(result)

    initial_rate = max(float(np.nanmedian(good["hydrogen_rate"])), 1.0e-9)
    x0 = np.r_[np.zeros(7), np.log(np.full(3, initial_rate))]
    lower = np.r_[np.full(7, -2.5), np.full(3, -30.0)]
    upper = np.r_[np.full(7, 2.5), np.full(3, 30.0)]
    result = least_squares(residuals, x0, bounds=(lower, upper), max_nfev=50000)
    r = residuals(result.x)
    ndf = max(len(r) - len(result.x), 1)
    chi2_ndf = float(np.sum(r * r) / ndf)
    try:
        covariance_all = np.linalg.inv(result.jac.T @ result.jac) * chi2_ndf
        covariance = covariance_all[:7, :7]
    except np.linalg.LinAlgError:
        covariance = np.full((7, 7), np.nan)
    # endtry
    values = {name: float(value) for name, value in zip(PHYSICS_PARAMETERS, result.x[:7])}
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
        state = subset[(subset.helicity == h) & (subset.target_sign == s)]
        for period in PERIODS:
            part = state[state.period == period].sort_values("phi_bin")
            if part.empty:
                continue
            # endif
            ax.errorbar(np.degrees(part.phi_center), part.hydrogen_rate,
                        yerr=part.hydrogen_rate_error, fmt="o-", ms=3, label=period)
        # endfor
        ax.set_title(f"h={h:+d}, target sign={s:+d}")
        ax.set_ylabel("Extracted H signal rate (counts / charge)")
        ax.grid(alpha=0.25)
    # endfor
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\phi$ (deg)")
    # endfor
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(f"Yield-check extracted free-H rates: kinematic bin {kin_bin}")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_summary(path: Path, results: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 4, figsize=(18, 12), sharex=True)
    axes = axes.flat
    for index, parameter in enumerate(PHYSICS_PARAMETERS):
        ax = axes[index]
        ax.errorbar(results.kin_bin, results[parameter], yerr=results[f"{parameter}_error"], fmt="o", ms=4)
        ax.axhline(0.0, lw=0.8)
        ax.set_title(nominal.PARAMETER_LABELS.get(parameter, parameter))
        ax.set_xlabel("Kinematic bin")
        ax.grid(alpha=0.2)
    # endfor
    for index in range(len(PHYSICS_PARAMETERS), 12):
        axes[index].axis("off")
    # endfor
    fig.suptitle("Direct-yield external extraction")
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

    area_rows: list[dict[str, Any]] = []
    area_lookup: dict[tuple[str, str, int, int, int, int], AreaFit] = {}
    plot_count = 0

    # Build all independent fit inputs in the parent process.  This keeps ROOT
    # I/O and pandas filtering out of workers and sends only the Mx2 arrays
    # required by curve_fit.  The expensive numerical fits then run in parallel.
    fit_tasks: list[tuple[Any, np.ndarray, float, float]] = []
    fit_metadata: dict[Any, dict[str, Any]] = {}
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
                for target in TARGETS:
                    signs = (-1, 1) if target == "NH3" else (0,)
                    for h in (-1, 1):
                        for s in signs:
                            key = (period, target, kin_bin, phi_bin, h, s)
                            subset = grouped.get(key)
                            if subset is None:
                                values = np.empty(0, dtype=float)
                                dep = {name: np.nan for name in ("rB", "rC", "rV", "rW")}
                            else:
                                values = subset.Mx2.to_numpy(dtype=float, copy=True)
                                dep = {name: float(np.nanmean(subset[name])) if len(subset) else np.nan
                                       for name in ("rB", "rC", "rV", "rW")}
                            # endif
                            fit_tasks.append((key, values, mu, sigma))
                            fit_metadata[key] = dict(
                                period=period, target=target, kin_bin=kin_bin, phi_bin=phi_bin,
                                phi_center=PHI_CENTERS[phi_bin], helicity=h, target_sign=s,
                                mu=mu, sigma=sigma, **dep,
                            )
                        # endfor
                    # endfor
                # endfor
            # endfor
        # endfor
    # endfor

    print(f"[fits] Gaussian spectra: {len(fit_tasks):,}; workers: {args.workers}", flush=True)
    if args.workers == 1:
        fit_results = map(fit_signal_area_worker, fit_tasks)
    else:
        executor = ProcessPoolExecutor(max_workers=args.workers)
        # chunksize amortizes process-IPC overhead without making tasks so coarse
        # that one worker can dominate the tail of the calculation.
        chunksize = max(1, len(fit_tasks) // (args.workers * 32))
        fit_results = executor.map(fit_signal_area_worker, fit_tasks, chunksize=chunksize)
    # endif

    try:
        for key, fit in fit_results:
            area_lookup[key] = fit
            meta = fit_metadata[key]
            area_rows.append(dict(
                **meta, area=fit.area, area_error=fit.area_error,
                background_area=fit.background_area, chi2_ndf=fit.chi2_ndf,
                n_events=fit.n_events, status=fit.status,
            ))
            if (args.max_spectrum_plots == 0 or plot_count < args.max_spectrum_plots) and fit.status == "ok":
                period, target, kin_bin, phi_bin, h, s = key
                plot_spectrum(
                    spectrum_plots / f"{period}_{target}_bin{kin_bin:02d}_phi{phi_bin:02d}_h{h:+d}_s{s:+d}.png",
                    fit, f"{period} {target}; bin {kin_bin}; phi bin {phi_bin}; h={h:+d}; s={s:+d}"
                )
                plot_count += 1
            # endif
        # endfor
    finally:
        if args.workers != 1:
            executor.shutdown(wait=True)
        # endif
    # endtry

    area_frame = pd.DataFrame(area_rows)
    area_frame.to_csv(tables / "gaussian_signal_areas.csv", index=False)

    rate_rows = []
    for period in PERIODS:
        for kin_bin in range(1, N_KIN_BINS + 1):
            for phi_bin in range(N_PHI_BINS):
                for h in (-1, 1):
                    # Auxiliary targets retain beam-helicity separation.
                    aux_areas = {}
                    aux_errors = {}
                    aux_charges = {}
                    for target in AUX_TARGETS:
                        fit = area_lookup[(period, target, kin_bin, phi_bin, h, 0)]
                        aux_areas[target] = fit.area
                        aux_errors[target] = fit.area_error
                        aux_charges[target] = charge_map[(period, target, h, 0)][0]
                    # endfor
                    for s in (-1, 1):
                        nh3_fit = area_lookup[(period, "NH3", kin_bin, phi_bin, h, s)]
                        areas = {"NH3": nh3_fit.area, **aux_areas}
                        errors = {"NH3": nh3_fit.area_error, **aux_errors}
                        qA, qpt = charge_map[(period, "NH3", h, s)]
                        charges = {"NH3": qA, **aux_charges}
                        h_rate = direct_method1_hydrogen_rate(areas, charges)
                        h_error = propagate_direct_rate_error(areas, errors, charges)
                        pt_eff = qpt / qA if qA > 0.0 else np.nan
                        nh3_rows = area_frame[
                            (area_frame.period == period) & (area_frame.target == "NH3")
                            & (area_frame.kin_bin == kin_bin) & (area_frame.phi_bin == phi_bin)
                            & (area_frame.helicity == h) & (area_frame.target_sign == s)
                        ]
                        dep = {name: float(nh3_rows.iloc[0][name]) if len(nh3_rows) else np.nan
                               for name in ("rB", "rC", "rV", "rW")}
                        rate_rows.append(dict(
                            period=period, kin_bin=kin_bin, phi_bin=phi_bin,
                            phi_center=PHI_CENTERS[phi_bin], helicity=h, target_sign=s,
                            charge=qA, pt_effective=pt_eff, hydrogen_rate=h_rate,
                            hydrogen_rate_error=h_error, nh3_signal_area=nh3_fit.area,
                            nh3_signal_area_error=nh3_fit.area_error, **dep,
                        ))
                    # endfor
                # endfor
            # endfor
        # endfor
    # endfor

    rates = pd.DataFrame(rate_rows)
    rates.to_csv(tables / "direct_hydrogen_rates.csv", index=False)

    result_rows = []
    covariance_payload: dict[str, Any] = {}
    for kin_bin in range(1, N_KIN_BINS + 1):
        subset = rates[rates.kin_bin == kin_bin]
        plot_four_state_rates(rate_plots / f"hydrogen_rates_bin{kin_bin:02d}.png", subset, kin_bin)
        values, covariance, chi2_ndf, npoints = fit_physics_bin(subset)
        row: dict[str, Any] = {"kin_bin": kin_bin, "chi2_ndf": chi2_ndf, "npoints": npoints}
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
    ax.set_title("Direct hydrogen-rate statistical significance")
    fig.tight_layout(); fig.savefig(diagnostics / "hydrogen_rate_significance.png", dpi=160); plt.close(fig)

    manifest = {
        "method": "direct five-target hydrogen-rate external check",
        "n_phi_bins": N_PHI_BINS,
        "phi_edges_rad": PHI_EDGES.tolist(),
        "does_not_use_dilution_factor": True,
        "auxiliary_targets_beam_helicity_separated": True,
        "nh3_beam_and_target_spin_separated": True,
        "signal_shape": "mu/sigma fixed from nominal mu +/- 2 sigma channel-selection cut",
        "statistics_warning": "First implementation uses local Gaussian-fit covariance; full bootstrap shared-target covariance is the planned upgrade.",
        "outputs": {
            "areas": str(tables / "gaussian_signal_areas.csv"),
            "rates": str(tables / "direct_hydrogen_rates.csv"),
            "physics": str(tables / "yield_check_structure_function_ratios.csv"),
        },
    }
    (out / "yield_check_manifest.json").write_text(json.dumps(manifest, indent=2))
    print(f"[done] yield-check outputs written to {out.resolve()}")


if __name__ == "__main__":
    main()
# endif
