#!/usr/bin/env python3
"""
validate_carbon_normalization_v7.py

Dedicated validation and stress-test program for the RGC exclusive e pi+
carbon-normalization dilution-factor method.

Run from

    RGC_enpi+/dilution_factor/

beside the production script

    determine_dilution_factor.py

The program deliberately imports the production dilution-factor module and
uses its ROOT loading, branch resolution, kinematic bins, final Mx2 cuts,
charge handling, and direct five-target dilution-factor expression.  This
prevents the validation from silently drifting away from the production
calculation.

The diagnostics answer four questions:

  1. How sensitive is Method 1 to both lower and upper control-window edges?
  2. Does a period-wide carbon scale suffice, or are xB-row, -t' column, or
     bin-local scales required?
  3. Does scaled carbon close between disjoint low-Mx2 subregions that are not
     expected to contain the exclusive neutron signal?
  4. What carbon transfer factor would be required in the signal window for
     Method 1 to reproduce Method 2?
  5. Is that required transfer factor consistent with one common correction,
     a period-dependent correction, or a kinematic dependence?
  6. Do cumulative low-side background integrals from scaled carbon and the
     direct five-target method agree without the unstable narrow-bin ratios?
  7. Which auxiliary-target normalization most strongly controls Method 2,
     and what coherent shift would be required to reproduce Method 1?
  8. Can a CH2-minus-carbon empirical hydrogen template describe the NH3
     spectrum, and does an independent template fit prefer a reduced carbon
     transfer factor beneath the neutron region?

No high-side NH3/carbon closure requirement is imposed.  Above the neutron
region, NH3 can contain real hydrogen-origin events from channels other than
exclusive ep -> e n pi+, so disagreement there is not interpreted as failure
of the carbon normalization.  Such regions are plotted only as descriptive
shape diagnostics when requested.

Statistical treatment
---------------------
Control-window scans use disjoint atomic Mx2 intervals constructed from every
requested lower and upper edge.  Their counts are independently Poisson
resampled and recombined in every replica.  This preserves covariance among
all overlapping control-window choices.  Final signal counts are also
Poisson-resampled, and every carbon scale and Method-1 dilution factor is
recomputed in each replica.  Method 2 is recalculated from the same target
replicas in the signal window, preserving Method-1/Method-2 covariance.

The program uses no more than seven workers.  It does not alter production
cuts or production dilution-factor outputs.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from collections import defaultdict
from dataclasses import dataclass
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

MAXIMUM_WORKERS = 7
DEFAULT_REPLICAS = 10000
DEFAULT_SEED = 27072026
DEFAULT_OUTPUT = Path("output/carbon_normalization_validation")
DEFAULT_MODULE = Path(__file__).resolve().with_name("determine_dilution_factor.py")

# The default production-like Method-1 normalization window is [0.00, 0.40).
# Command-line overrides may replace either edge. These scans explicitly vary
# both edges, and duplicate entries are removed during canonicalization.
DEFAULT_CONTROL_WINDOWS: tuple[tuple[float, float], ...] = (
    (0.00, 0.25),
    (0.00, 0.30),
    (0.00, 0.35),
    (0.00, 0.40),
    (0.00, 0.45),
    (0.00, 0.50),
    (-0.20, 0.40),
    (-0.10, 0.40),
    (0.05, 0.40),
    (0.10, 0.40),
)

# Disjoint low-Mx2 closure regions.  These are deliberately below the neutron
# peak.  The program tests transfer among them but does not assume any one is
# automatically signal-free; the results are diagnostics.
DEFAULT_LOW_REGIONS: tuple[tuple[float, float], ...] = (
    (-0.20, 0.00),
    (0.00, 0.20),
    (0.20, 0.40),
)

DEFAULT_SPECTRUM_RANGE = (-0.20, 1.50)
DEFAULT_SPECTRUM_BINS = 80
DEFAULT_CUMULATIVE_UPPER_MIN = 0.15
DEFAULT_CUMULATIVE_UPPER_MAX = 0.65
DEFAULT_CUMULATIVE_UPPER_STEP = 0.025
PERIOD_GROUP = "period"
XB_GROUP = "xB-row"
TPRIME_GROUP = "tprime-column"
BIN_GROUP = "bin-local"
GROUPINGS = (PERIOD_GROUP, XB_GROUP, TPRIME_GROUP, BIN_GROUP)
SENSITIVITY_FRACTION = 0.01
REQUIRED_SCALE_MIN = 0.20
REQUIRED_SCALE_MAX = 2.50
REQUIRED_SCALE_GRID = 241
TEMPLATE_FIT_MIN = -0.20
TEMPLATE_FIT_MAX = 1.35
TEMPLATE_POLYNOMIAL_ORDER = 2
NOMINAL_R_CN = 6.0 / 7.0
DEFAULT_R_CN_MIN = 0.50
DEFAULT_R_CN_MAX = 1.10
DEFAULT_R_CN_POINTS = 121
DEFAULT_NOMINAL_CONTROL_WINDOW = (0.0, 0.4)
NOMINAL_CONTROL_WINDOW = DEFAULT_NOMINAL_CONTROL_WINDOW
R_CN_BISECTION_ITERATIONS = 48


@dataclass(frozen=True)
class ValidationCounts:
    """Observed counts needed for one period and target."""

    # shape: (24 bins, n atomic intervals)
    atomic_control: np.ndarray
    # shape: (24 bins, 3 cut variations)
    signal: np.ndarray
    # shape: (24 bins, n low regions)
    low_regions: np.ndarray
    # shape: (24 bins, n spectrum bins)
    spectrum: np.ndarray


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, payload: Any) -> None:
    ensure_dir(path.parent)
    path.write_text(json.dumps(payload, indent=2, allow_nan=False), encoding="utf-8")


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def import_production_module(path: Path):
    if not path.is_file():
        raise FileNotFoundError(
            f"Production dilution-factor module not found: {path}. "
            "Place this validation script beside determine_dilution_factor.py "
            "or pass --dilution-module."
        )
    spec = importlib.util.spec_from_file_location("rgc_dilution_production", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import production module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_window(text: str) -> tuple[float, float]:
    pieces = text.replace(",", ":").split(":")
    if len(pieces) != 2:
        raise argparse.ArgumentTypeError("Window must be MIN:MAX")
    low, high = map(float, pieces)
    if not low < high:
        raise argparse.ArgumentTypeError("Window requires MIN < MAX")
    return low, high


def cumulative_upper_edges(minimum: float, maximum: float, step: float) -> list[float]:
    if step <= 0:
        raise ValueError("Cumulative upper-edge step must be positive")
    if not minimum < maximum:
        raise ValueError("Cumulative upper-edge minimum must be below maximum")
    count = int(round((maximum - minimum) / step))
    values = [minimum + i * step for i in range(count + 1)]
    if values[-1] < maximum - 1e-10:
        values.append(maximum)
    return [round(value, 10) for value in values]


def canonical_windows(
    values: list[tuple[float, float]] | None,
    cumulative_edges: list[float],
) -> list[tuple[float, float]]:
    source = list(DEFAULT_CONTROL_WINDOWS if not values else values)
    source.extend((0.0, upper) for upper in cumulative_edges)
    unique = sorted({(round(a, 10), round(b, 10)) for a, b in source})
    if NOMINAL_CONTROL_WINDOW not in unique:
        unique.append(NOMINAL_CONTROL_WINDOW)
        unique.sort()
    return unique


def structured_windows(windows: list[tuple[float, float]]) -> list[tuple[float, float]]:
    requested = set(DEFAULT_CONTROL_WINDOWS)
    return [window for window in windows if window in requested]


def combined_indices(number_of_bins: int, n_t: int) -> tuple[np.ndarray, np.ndarray]:
    indices = np.arange(number_of_bins, dtype=int)
    return indices // n_t, indices % n_t


def group_ids(grouping: str, number_of_bins: int, n_t: int) -> np.ndarray:
    xb_index, t_index = combined_indices(number_of_bins, n_t)
    if grouping == PERIOD_GROUP:
        return np.zeros(number_of_bins, dtype=int)
    if grouping == XB_GROUP:
        return xb_index
    if grouping == TPRIME_GROUP:
        return t_index
    if grouping == BIN_GROUP:
        return np.arange(number_of_bins, dtype=int)
    raise ValueError(f"Unknown grouping: {grouping}")


def group_count(grouping: str, n_x: int, n_t: int) -> int:
    return {
        PERIOD_GROUP: 1,
        XB_GROUP: n_x,
        TPRIME_GROUP: n_t,
        BIN_GROUP: n_x * n_t,
    }[grouping]


def ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    output = np.full(np.broadcast_shapes(numerator.shape, denominator.shape), np.nan)
    np.divide(numerator, denominator, out=output, where=denominator != 0)
    return output


def bin_mask(dataset: dict[str, np.ndarray], x_range: tuple[float, float], t_range: tuple[float, float]) -> np.ndarray:
    return (
        (dataset["xB"] >= x_range[0])
        & (dataset["xB"] < x_range[1])
        & (dataset["minus_tprime_gev2"] >= t_range[0])
        & (dataset["minus_tprime_gev2"] < t_range[1])
    )


def build_atomic_edges(windows: list[tuple[float, float]]) -> np.ndarray:
    edges = sorted({value for window in windows for value in window})
    if len(edges) < 2:
        raise RuntimeError("At least two unique control-window edges are required")
    return np.asarray(edges, dtype=float)


def window_interval_mask(
    atomic_edges: np.ndarray,
    window: tuple[float, float],
) -> np.ndarray:
    lows = atomic_edges[:-1]
    highs = atomic_edges[1:]
    return (lows >= window[0] - 1e-12) & (highs <= window[1] + 1e-12)


def count_dataset(
    dataset: dict[str, np.ndarray],
    cuts: dict[tuple[str, int, int], dict[str, tuple[float, float]]],
    period: str,
    xb_bins: tuple[tuple[float, float], ...],
    t_bins: tuple[tuple[float, float], ...],
    atomic_edges: np.ndarray,
    low_regions: tuple[tuple[float, float], ...],
    spectrum_edges: np.ndarray,
) -> ValidationCounts:
    """Vectorized counting for one target sample.

    Events are assigned to the fixed xB and -t' bins once with searchsorted.
    Two-dimensional bincount histograms then replace the former repeated
    24-bin masking and histogram loop.  Only the bin-dependent signal windows
    require a short loop over the 24 analysis bins.
    """
    n_x = len(xb_bins)
    n_t = len(t_bins)
    n_bins = n_x * n_t
    x_edges = np.asarray([xb_bins[0][0], *[value[1] for value in xb_bins]], dtype=float)
    t_edges = np.asarray([t_bins[0][0], *[value[1] for value in t_bins]], dtype=float)

    x = np.asarray(dataset["xB"], dtype=float)
    t = np.asarray(dataset["minus_tprime_gev2"], dtype=float)
    mx2 = np.asarray(dataset["mx2_gev2"], dtype=float)
    ix = np.searchsorted(x_edges, x, side="right") - 1
    it = np.searchsorted(t_edges, t, side="right") - 1
    valid = (ix >= 0) & (ix < n_x) & (it >= 0) & (it < n_t) & np.isfinite(mx2)
    combined = (ix[valid] * n_t + it[valid]).astype(np.int64, copy=False)
    values = mx2[valid]

    def histogram_by_analysis_bin(edges: np.ndarray) -> np.ndarray:
        interval = np.searchsorted(edges, values, side="right") - 1
        good = (interval >= 0) & (interval < len(edges) - 1)
        flat = combined[good] * (len(edges) - 1) + interval[good]
        return np.bincount(
            flat, minlength=n_bins * (len(edges) - 1)
        ).reshape(n_bins, len(edges) - 1).astype(np.int64, copy=False)

    atomic = histogram_by_analysis_bin(atomic_edges)
    spectrum = histogram_by_analysis_bin(spectrum_edges)
    low = np.stack([
        np.bincount(
            combined[(values >= lo) & (values < hi)], minlength=n_bins
        )
        for lo, hi in low_regions
    ], axis=1).astype(np.int64, copy=False)

    signal = np.zeros((n_bins, 3), dtype=np.int64)
    order = ("tight", "nominal", "loose")
    for bin_index in range(n_bins):
        bin_values = values[combined == bin_index]
        ix_bin, it_bin = divmod(bin_index, n_t)
        cut = cuts[(period, ix_bin, it_bin)]
        lows = np.asarray([cut[name][0] for name in order])
        highs = np.asarray([cut[name][1] for name in order])
        signal[bin_index] = np.count_nonzero(
            (bin_values[:, np.newaxis] >= lows)
            & (bin_values[:, np.newaxis] < highs),
            axis=0,
        )

    return ValidationCounts(atomic, signal, low, spectrum)

def count_worker(args: tuple[Any, ...]) -> tuple[tuple[str, str], ValidationCounts]:
    (
        key,
        dataset,
        cuts,
        xb_bins,
        t_bins,
        atomic_edges,
        low_regions,
        spectrum_edges,
    ) = args
    period, _target = key
    return key, count_dataset(
        dataset, cuts, period, xb_bins, t_bins, atomic_edges, low_regions, spectrum_edges
    )


def build_all_counts(
    loaded: dict[tuple[str, str], Any],
    cuts: dict[tuple[str, int, int], Any],
    module: Any,
    atomic_edges: np.ndarray,
    low_regions: tuple[tuple[float, float], ...],
    spectrum_edges: np.ndarray,
    workers: int,
) -> dict[tuple[str, str], ValidationCounts]:
    plain_cuts = {
        key: {name: tuple(entry.interval(name)) for name in ("tight", "nominal", "loose")}
        for key, entry in cuts.items()
    }
    jobs = [
        (
            key,
            {
                "xB": np.asarray(dataset.xB),
                "minus_tprime_gev2": np.asarray(dataset.minus_tprime_gev2),
                "mx2_gev2": np.asarray(dataset.mx2_gev2),
            },
            plain_cuts,
            module.XB_BINS,
            module.MINUS_TPRIME_BINS_GEV2,
            atomic_edges,
            low_regions,
            spectrum_edges,
        )
        for key, dataset in loaded.items()
    ]
    output: dict[tuple[str, str], ValidationCounts] = {}
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as executor:
        futures = [executor.submit(count_worker, job) for job in jobs]
        for completed, future in enumerate(as_completed(futures), start=1):
            key, value = future.result()
            output[key] = value
            print(f"Counted {completed:2d}/{len(jobs)}: {key[0]} {key[1]}")
    return output


def control_projection_matrix(
    atomic_edges: np.ndarray,
    windows: list[tuple[float, float]],
) -> np.ndarray:
    """Return an atomic-interval to control-window incidence matrix."""
    return np.stack(
        [window_interval_mask(atomic_edges, window) for window in windows], axis=1
    ).astype(np.int8, copy=False)


def controls_for_windows(
    atomic: np.ndarray,
    atomic_edges: np.ndarray,
    windows: list[tuple[float, float]],
) -> np.ndarray:
    projection = control_projection_matrix(atomic_edges, windows)
    return np.einsum("...i,iw->...w", atomic, projection, optimize=True)

def scales_for_grouping(
    nh3_control: np.ndarray,
    carbon_control: np.ndarray,
    grouping: str,
    n_x: int,
    n_t: int,
) -> np.ndarray:
    """Return one scale per bin and window, broadcast from the chosen groups."""
    n_bins = n_x * n_t
    gids = group_ids(grouping, n_bins, n_t)
    scales = np.full_like(nh3_control, np.nan, dtype=float)
    for gid in range(group_count(grouping, n_x, n_t)):
        members = gids == gid
        n = np.sum(nh3_control[members], axis=0)
        c = np.sum(carbon_control[members], axis=0)
        scales[members] = ratio(n, c)
    return scales


def method1_values(
    nh3_signal: np.ndarray,
    carbon_signal: np.ndarray,
    scale_by_bin: np.ndarray,
) -> np.ndarray:
    # signal shape (..., 24, 3); scale shape (..., 24, n_windows)
    numerator = nh3_signal[..., :, np.newaxis, :] - (
        scale_by_bin[..., :, :, np.newaxis] * carbon_signal[..., :, np.newaxis, :]
    )
    return ratio(numerator, nh3_signal[..., :, np.newaxis, :])


def method2_values(
    signal_counts: dict[str, np.ndarray],
    period: str,
    module: Any,
    charge_fractions: dict[str, dict[str, float]],
) -> np.ndarray:
    return module.method2_equation10_from_counts(signal_counts, charge_fractions[period])



def method2_with_rcn(
    signal_counts: dict[str, np.ndarray],
    charge_fractions_period: dict[str, float],
    module: Any,
    r_cn: float | np.ndarray,
    nominal_r_cn: float = NOMINAL_R_CN,
) -> np.ndarray:
    """Evaluate Method 2 for an adjustable carbon-to-nitrogen relation.

    The vetted production expression is normalized at

        sigma_C = nominal_r_cn * sigma_N.

    Holding the measured pure-carbon yield fixed, changing the assumed ratio to
    ``r_cn`` changes the nitrogen cross section inferred from that carbon proxy
    by ``nominal_r_cn / r_cn``.  The production expression is therefore
    reevaluated with an effective carbon-proxy yield

        N_C(effective) = N_C * nominal_r_cn / r_cn.

    At ``r_cn == nominal_r_cn`` this reproduces the production Method-2 value
    exactly.  This is a controlled one-parameter study of the fixed C--N
    relation; it is not a Faraday-cup or target-yield renormalization.
    """
    r = np.asarray(r_cn, dtype=float)
    scale = ratio(nominal_r_cn, r)
    shifted = {target: np.asarray(values, dtype=float) for target, values in signal_counts.items()}
    shifted = dict(shifted)
    shifted["C"] = shifted["C"] * scale
    return module.method2_equation10_from_counts(shifted, charge_fractions_period)


def solve_required_rcn(
    signal_counts: dict[str, np.ndarray],
    target_method1: np.ndarray,
    charge_fractions_period: dict[str, float],
    module: Any,
    r_min: float,
    r_max: float,
    nominal_r_cn: float = NOMINAL_R_CN,
) -> np.ndarray:
    """Solve f2(r_CN) = f1 with vectorized bisection.

    Elements without a finite sign-changing bracket are returned as NaN.
    """
    target = np.asarray(target_method1, dtype=float)
    low = np.full(target.shape, float(r_min), dtype=float)
    high = np.full(target.shape, float(r_max), dtype=float)
    f_low = method2_with_rcn(signal_counts, charge_fractions_period, module, low, nominal_r_cn) - target
    f_high = method2_with_rcn(signal_counts, charge_fractions_period, module, high, nominal_r_cn) - target
    valid = np.isfinite(target) & np.isfinite(f_low) & np.isfinite(f_high) & (f_low * f_high <= 0.0)
    for _ in range(R_CN_BISECTION_ITERATIONS):
        mid = 0.5 * (low + high)
        f_mid = method2_with_rcn(signal_counts, charge_fractions_period, module, mid, nominal_r_cn) - target
        same_side = np.signbit(f_mid) == np.signbit(f_low)
        update_low = valid & same_side
        update_high = valid & ~same_side
        low = np.where(update_low, mid, low)
        f_low = np.where(update_low, f_mid, f_low)
        high = np.where(update_high, mid, high)
        f_high = np.where(update_high, f_mid, f_high)
    answer = 0.5 * (low + high)
    return np.where(valid, answer, np.nan)

def summarize(samples: np.ndarray) -> dict[str, np.ndarray]:
    finite = np.isfinite(samples)
    valid = np.mean(finite, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.nanmean(samples, axis=0)
        std = np.nanstd(samples, axis=0, ddof=1)
        p16 = np.nanpercentile(samples, 16, axis=0)
        p50 = np.nanpercentile(samples, 50, axis=0)
        p84 = np.nanpercentile(samples, 84, axis=0)
    return {"mean": mean, "std": std, "p16": p16, "p50": p50, "p84": p84, "valid": valid}


def validation_bootstrap_worker(
    period: str,
    observed: dict[str, dict[str, np.ndarray]],
    windows: list[tuple[float, float]],
    atomic_edges: np.ndarray,
    groupings: tuple[str, ...],
    n_x: int,
    n_t: int,
    charge_fractions_period: dict[str, float],
    replicas: int,
    seed: int,
    module_path: str,
    chunk_index: int,
    nominal_control_window: tuple[float, float],
) -> dict[str, Any]:
    """Generate one independent bootstrap chunk for one period."""
    module = import_production_module(Path(module_path))
    rng = np.random.default_rng(seed)
    n_bins = n_x * n_t
    n_windows = len(windows)
    method1 = {
        grouping: np.full((replicas, n_bins, n_windows, 3), np.nan)
        for grouping in groupings
    }
    method2 = np.full((replicas, n_bins, 3), np.nan)
    alpha_required = np.full((replicas, n_bins, 3), np.nan)
    transfer_ratio = np.full((replicas, n_bins, 3), np.nan)
    kappa_period = np.full((replicas, 3), np.nan)
    required_rcn = np.full((replicas, n_bins, 3), np.nan)
    required_rcn_period = np.full((replicas, 3), np.nan)
    alpha = {
        grouping: np.full((replicas, n_bins, n_windows), np.nan)
        for grouping in groupings
    }

    projection = control_projection_matrix(atomic_edges, windows)
    nominal_window_index = windows.index(nominal_control_window)
    grouping_metadata = {}
    for grouping in groupings:
        gids = group_ids(grouping, n_bins, n_t)
        n_groups = group_count(grouping, n_x, n_t)
        membership = np.equal.outer(np.arange(n_groups), gids).astype(np.int8)
        grouping_metadata[grouping] = (gids, membership)

    # Keep transient Poisson arrays bounded while using vectorized linear
    # algebra within each batch.
    batch_size = min(512, replicas)
    for start in range(0, replicas, batch_size):
        stop = min(replicas, start + batch_size)
        batch = stop - start
        control_by_target: dict[str, np.ndarray] = {}
        signal_by_target: dict[str, np.ndarray] = {}
        for target, arrays in observed.items():
            atomic_rep = rng.poisson(
                arrays["atomic"][np.newaxis, ...],
                size=(batch,) + arrays["atomic"].shape,
            )
            signal_by_target[target] = rng.poisson(
                arrays["signal"][np.newaxis, ...],
                size=(batch,) + arrays["signal"].shape,
            )
            control_by_target[target] = np.einsum(
                "bni,iw->bnw", atomic_rep, projection, optimize=True
            )

        for grouping in groupings:
            gids, membership = grouping_metadata[grouping]
            nh3_group = np.einsum(
                "gn,bnw->bgw", membership, control_by_target["NH3"], optimize=True
            )
            carbon_group = np.einsum(
                "gn,bnw->bgw", membership, control_by_target["C"], optimize=True
            )
            scale_by_group = ratio(nh3_group, carbon_group)
            scale = scale_by_group[:, gids, :]
            alpha[grouping][start:stop] = scale
            method1[grouping][start:stop] = method1_values(
                signal_by_target["NH3"], signal_by_target["C"], scale
            )

        m2_batch = module.method2_equation10_from_counts(
            signal_by_target, charge_fractions_period
        )
        method2[start:stop] = m2_batch
        required = ratio(
            signal_by_target["NH3"] * (1.0 - m2_batch),
            signal_by_target["C"],
        )
        alpha_required[start:stop] = required
        nominal_period_alpha = alpha[PERIOD_GROUP][
            start:stop, :, nominal_window_index
        ]
        transfer_ratio[start:stop] = ratio(
            required, nominal_period_alpha[..., np.newaxis]
        )
        predicted_background = (
            nominal_period_alpha[..., np.newaxis] * signal_by_target["C"]
        )
        required_background = signal_by_target["NH3"] * (1.0 - m2_batch)
        kappa_period[start:stop] = ratio(
            np.sum(required_background, axis=1),
            np.sum(predicted_background, axis=1),
        )
        nominal_m1 = method1[PERIOD_GROUP][
            start:stop, :, nominal_window_index, :
        ]
        required_rcn[start:stop] = solve_required_rcn(
            signal_by_target,
            nominal_m1,
            charge_fractions_period,
            module,
            DEFAULT_R_CN_MIN,
            DEFAULT_R_CN_MAX,
        )
        integrated_counts = {
            target: np.sum(values, axis=1)
            for target, values in signal_by_target.items()
        }
        integrated_m1 = ratio(
            np.sum(signal_by_target["NH3"], axis=1)
            - np.sum(
                nominal_period_alpha[..., np.newaxis] * signal_by_target["C"],
                axis=1,
            ),
            np.sum(signal_by_target["NH3"], axis=1),
        )
        required_rcn_period[start:stop] = solve_required_rcn(
            integrated_counts,
            integrated_m1,
            charge_fractions_period,
            module,
            DEFAULT_R_CN_MIN,
            DEFAULT_R_CN_MAX,
        )

    return {
        "period": period,
        "chunk_index": chunk_index,
        "method1": method1,
        "method2": method2,
        "alpha": alpha,
        "alpha_required": alpha_required,
        "transfer_ratio": transfer_ratio,
        "kappa_period": kappa_period,
        "required_rcn": required_rcn,
        "required_rcn_period": required_rcn_period,
    }


def _replica_chunks(replicas: int, jobs: int) -> list[int]:
    jobs = max(1, min(jobs, replicas))
    base, remainder = divmod(replicas, jobs)
    return [base + (index < remainder) for index in range(jobs)]


def run_bootstrap(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    module_path: Path,
    windows: list[tuple[float, float]],
    atomic_edges: np.ndarray,
    charge_fractions: dict[str, dict[str, float]],
    replicas: int,
    seed: int,
    workers: int,
) -> dict[str, Any]:
    """Parallel bootstrap across period/replica chunks, capped at seven workers."""
    n_x = len(module.XB_BINS)
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    workers = max(1, min(MAXIMUM_WORKERS, workers))

    # Use enough chunks to occupy every available worker while keeping each
    # returned object reasonably large.  Replica order has no statistical
    # significance, but chunks are reassembled deterministically.
    chunks_per_period = max(1, math.ceil(workers / len(module.PERIODS)))
    chunk_sizes = _replica_chunks(replicas, chunks_per_period)
    seed_sequence = np.random.SeedSequence(seed)
    child_seeds = iter(seed_sequence.spawn(len(module.PERIODS) * chunks_per_period))

    futures = {}
    with ProcessPoolExecutor(max_workers=workers) as executor:
        for period in module.PERIODS:
            observed = {
                target: {
                    "atomic": counts[(period, target)].atomic_control,
                    "signal": counts[(period, target)].signal,
                }
                for target in module.TARGETS
            }
            for chunk_index, chunk_size in enumerate(chunk_sizes):
                child_seed = int(next(child_seeds).generate_state(1, dtype=np.uint64)[0])
                future = executor.submit(
                    validation_bootstrap_worker,
                    period, observed, windows, atomic_edges, GROUPINGS,
                    n_x, n_t, charge_fractions[period], chunk_size, child_seed,
                    str(module_path), chunk_index, NOMINAL_CONTROL_WINDOW,
                )
                futures[future] = (period, chunk_index, chunk_size)

        chunks: dict[str, list[dict[str, Any]]] = defaultdict(list)
        for future in as_completed(futures):
            period, chunk_index, chunk_size = futures[future]
            chunks[period].append(future.result())
            print(
                f"Completed bootstrap chunk {chunk_index + 1}/{chunks_per_period} "
                f"for {period} ({chunk_size:,} replicas)."
            )

    output: dict[str, Any] = {}
    for period in module.PERIODS:
        ordered = sorted(chunks[period], key=lambda item: item["chunk_index"])
        output[period] = {
            "method1": {
                grouping: np.concatenate(
                    [item["method1"][grouping] for item in ordered], axis=0
                )
                for grouping in GROUPINGS
            },
            "method2": np.concatenate([item["method2"] for item in ordered], axis=0),
            "alpha": {
                grouping: np.concatenate(
                    [item["alpha"][grouping] for item in ordered], axis=0
                )
                for grouping in GROUPINGS
            },
            "alpha_required": np.concatenate(
                [item["alpha_required"] for item in ordered], axis=0
            ),
            "transfer_ratio": np.concatenate(
                [item["transfer_ratio"] for item in ordered], axis=0
            ),
            "kappa_period": np.concatenate(
                [item["kappa_period"] for item in ordered], axis=0
            ),
            "required_rcn": np.concatenate(
                [item["required_rcn"] for item in ordered], axis=0
            ),
            "required_rcn_period": np.concatenate(
                [item["required_rcn_period"] for item in ordered], axis=0
            ),
        }
        print(f"Assembled {replicas:,} validation replicas for {period}.")
    return output

def central_estimators(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    windows: list[tuple[float, float]],
    atomic_edges: np.ndarray,
    charge_fractions: dict[str, dict[str, float]],
) -> dict[str, Any]:
    output: dict[str, Any] = {}
    n_x = len(module.XB_BINS)
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    for period in module.PERIODS:
        controls = {
            target: controls_for_windows(counts[(period, target)].atomic_control, atomic_edges, windows)
            for target in module.TARGETS
        }
        signals = {target: counts[(period, target)].signal.astype(float) for target in module.TARGETS}
        result: dict[str, Any] = {"method1": {}, "alpha": {}}
        for grouping in GROUPINGS:
            scale = scales_for_grouping(controls["NH3"], controls["C"], grouping, n_x, n_t)
            result["alpha"][grouping] = scale
            result["method1"][grouping] = method1_values(signals["NH3"], signals["C"], scale)
        result["method2"] = method2_values(signals, period, module, charge_fractions)
        result["alpha_required"] = ratio(
            signals["NH3"] * (1.0 - result["method2"]),
            signals["C"],
        )
        nominal_window_index = windows.index(NOMINAL_CONTROL_WINDOW)
        nominal_period_alpha = result["alpha"][PERIOD_GROUP][:, nominal_window_index]
        result["transfer_ratio"] = ratio(
            result["alpha_required"],
            nominal_period_alpha[:, np.newaxis],
        )
        predicted_nominal_background = nominal_period_alpha[:, np.newaxis] * signals["C"]
        required_background = signals["NH3"] * (1.0 - result["method2"])
        result["kappa_period"] = ratio(
            np.sum(required_background, axis=0),
            np.sum(predicted_nominal_background, axis=0),
        )
        nominal_m1 = result["method1"][PERIOD_GROUP][:, nominal_window_index, :]
        result["required_rcn"] = solve_required_rcn(
            signals, nominal_m1, charge_fractions[period], module,
            DEFAULT_R_CN_MIN, DEFAULT_R_CN_MAX,
        )
        integrated_counts = {target: np.sum(values, axis=0) for target, values in signals.items()}
        integrated_m1 = ratio(
            np.sum(signals["NH3"], axis=0)
            - np.sum(nominal_period_alpha[:, np.newaxis] * signals["C"], axis=0),
            np.sum(signals["NH3"], axis=0),
        )
        result["integrated_method1"] = integrated_m1
        result["required_rcn_period"] = solve_required_rcn(
            integrated_counts, integrated_m1, charge_fractions[period], module,
            DEFAULT_R_CN_MIN, DEFAULT_R_CN_MAX,
        )
        output[period] = result
    return output


def flatten_bin_labels(module: Any) -> list[str]:
    labels = []
    for ix, xb in enumerate(module.XB_BINS):
        for it, tr in enumerate(module.MINUS_TPRIME_BINS_GEV2):
            labels.append(f"{ix+1}:{it+1}")
    return labels


def window_label(window: tuple[float, float]) -> str:
    return f"[{window[0]:.2f}, {window[1]:.2f})"


def plot_window_scan_summary(
    output: Path,
    module: Any,
    windows: list[tuple[float, float]],
    central: dict[str, Any],
    bootstrap: dict[str, Any],
) -> None:
    nominal_window = windows.index(NOMINAL_CONTROL_WINDOW)
    nominal_cut = list(module.CUT_VARIATIONS).index("nominal")
    display_windows = structured_windows(windows)
    display_indices = [windows.index(window) for window in display_windows]
    x = np.arange(len(display_windows))
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(14, 11), sharex=True, sharey=True)
    for ax, period in zip(axes, module.PERIODS):
        all_values = central[period]["method1"][PERIOD_GROUP][:, :, nominal_cut]
        values = all_values[:, display_indices]
        reference = all_values[:, nominal_window][:, np.newaxis]
        shift = values - reference
        all_replicas = bootstrap[period]["method1"][PERIOD_GROUP][:, :, :, nominal_cut]
        replicas = all_replicas[:, :, display_indices]
        replica_shift = replicas - all_replicas[:, :, nominal_window][:, :, np.newaxis]
        mean_shift = np.nanmean(shift, axis=0)
        max_shift = np.nanmax(np.abs(shift), axis=0)
        mean_rep = np.nanmean(replica_shift, axis=1)
        mean_err = np.nanstd(mean_rep, axis=0, ddof=1)
        ax.errorbar(x, mean_shift, yerr=mean_err, fmt="o", label="Mean bin shift")
        ax.plot(x, max_shift, "s", label="Maximum |bin shift|")
        ax.axhline(0.0, linewidth=1)
        ax.set_ylabel(r"$\Delta f_1$")
        ax.set_title(module.PERIOD_LABELS[period])
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xticks(x)
    axes[-1].set_xticklabels([window_label(w) for w in display_windows], rotation=45, ha="right")
    axes[-1].set_xlabel(r"Carbon normalization window in $M_X^2$ (GeV$^2$)")
    fig.suptitle("Method-1 sensitivity to both carbon-control boundaries")
    fig.tight_layout()
    fig.savefig(output / "control_window_scan_summary.png", dpi=180)
    plt.close(fig)


def plot_grouping_comparison(
    output: Path,
    module: Any,
    windows: list[tuple[float, float]],
    central: dict[str, Any],
    bootstrap: dict[str, Any],
) -> None:
    w = windows.index(NOMINAL_CONTROL_WINDOW)
    cut = list(module.CUT_VARIATIONS).index("nominal")
    x = np.arange(1, module.NUMBER_OF_BINS + 1)
    offsets = np.linspace(-0.24, 0.24, len(GROUPINGS) + 1)
    for period in module.PERIODS:
        fig, ax = plt.subplots(figsize=(16, 6))
        for offset, grouping in zip(offsets, GROUPINGS):
            y = central[period]["method1"][grouping][:, w, cut]
            yerr = np.nanstd(bootstrap[period]["method1"][grouping][:, :, w, cut], axis=0, ddof=1)
            ax.errorbar(x + offset, y, yerr=yerr, fmt="o", capsize=2, label=grouping)
        y2 = central[period]["method2"][:, cut]
        y2err = np.nanstd(bootstrap[period]["method2"][:, :, cut], axis=0, ddof=1)
        ax.errorbar(x + offsets[-1], y2, yerr=y2err, fmt="o", capsize=2, label="direct five-target")
        ax.set_ylim(0.1, 0.6)
        ax.set_xticks(x)
        ax.set_xlabel("Combined kinematic-bin number")
        ax.set_ylabel("Dilution factor")
        ax.set_title(f"{module.PERIOD_LABELS[period]} — normalization granularity, nominal cuts")
        ax.grid(alpha=0.25)
        ax.legend(ncol=3)
        fig.tight_layout()
        fig.savefig(output / f"normalization_grouping_{period}.png", dpi=180)
        plt.close(fig)


def plot_method_difference_summary(
    output: Path,
    module: Any,
    windows: list[tuple[float, float]],
    central: dict[str, Any],
    bootstrap: dict[str, Any],
) -> None:
    w = windows.index(NOMINAL_CONTROL_WINDOW)
    cut = list(module.CUT_VARIATIONS).index("nominal")
    x = np.arange(1, module.NUMBER_OF_BINS + 1)
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(16, 11), sharex=True, sharey=True)
    for ax, period in zip(axes, module.PERIODS):
        for grouping in GROUPINGS:
            delta = central[period]["method1"][grouping][:, w, cut] - central[period]["method2"][:, cut]
            rep_delta = (
                bootstrap[period]["method1"][grouping][:, :, w, cut]
                - bootstrap[period]["method2"][:, :, cut]
            )
            err = np.nanstd(rep_delta, axis=0, ddof=1)
            ax.errorbar(x, delta, yerr=err, fmt="o", capsize=2, label=grouping)
        ax.axhline(0.0, linewidth=1)
        ax.set_ylabel(r"$f_1-f_2$")
        ax.set_title(module.PERIOD_LABELS[period])
        ax.grid(alpha=0.25)
        ax.legend(ncol=4)
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(x)
    fig.suptitle("Does finer carbon normalization reduce the Method-1/Method-2 difference?")
    fig.tight_layout()
    fig.savefig(output / "method_difference_by_normalization_granularity.png", dpi=180)
    plt.close(fig)


def low_region_closure_table(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    low_regions: tuple[tuple[float, float], ...],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    for period in module.PERIODS:
        nh3 = counts[(period, "NH3")].low_regions.astype(float)
        carbon = counts[(period, "C")].low_regions.astype(float)
        for normalization_index, normalization_region in enumerate(low_regions):
            alpha = ratio(np.sum(nh3[:, normalization_index]), np.sum(carbon[:, normalization_index]))
            for validation_index, validation_region in enumerate(low_regions):
                if validation_index == normalization_index:
                    continue
                predicted = alpha * carbon[:, validation_index]
                observed = nh3[:, validation_index]
                closure = ratio(predicted, observed)
                approximate_variance = ratio(predicted, observed**2) + ratio(predicted**2, observed**3)
                for bin_index in range(module.NUMBER_OF_BINS):
                    ix, it = divmod(bin_index, n_t)
                    rows.append({
                        "period": period,
                        "bin_number": bin_index + 1,
                        "x_index": ix,
                        "tprime_index": it,
                        "normalization_region_min": normalization_region[0],
                        "normalization_region_max": normalization_region[1],
                        "validation_region_min": validation_region[0],
                        "validation_region_max": validation_region[1],
                        "alpha_period": alpha,
                        "nh3_validation_count": observed[bin_index],
                        "carbon_validation_count": carbon[bin_index, validation_index],
                        "scaled_carbon_over_nh3": closure[bin_index],
                        "approx_stat_uncertainty": math.sqrt(approximate_variance[bin_index]) if np.isfinite(approximate_variance[bin_index]) and approximate_variance[bin_index] >= 0 else math.nan,
                    })
    return pd.DataFrame(rows)


def plot_low_region_closure(output: Path, table: pd.DataFrame, module: Any) -> None:
    combinations = table[[
        "normalization_region_min", "normalization_region_max",
        "validation_region_min", "validation_region_max"
    ]].drop_duplicates().to_dict("records")
    for period in module.PERIODS:
        fig, axes = plt.subplots(len(combinations), 1, figsize=(15, 3.2 * len(combinations)), sharex=True)
        axes = np.atleast_1d(axes)
        for ax, combo in zip(axes, combinations):
            selected = table[(table.period == period)]
            for key, value in combo.items():
                selected = selected[np.isclose(selected[key], value)]
            x = selected.bin_number.to_numpy()
            y = selected.scaled_carbon_over_nh3.to_numpy()
            e = selected.approx_stat_uncertainty.to_numpy()
            ax.errorbar(x, y, yerr=e, fmt="o", capsize=2)
            ax.axhline(1.0, linewidth=1)
            ax.set_ylabel("Scaled C / NH$_3$")
            ax.set_title(
                f"Normalize [{combo['normalization_region_min']:.2f}, {combo['normalization_region_max']:.2f}); "
                f"validate [{combo['validation_region_min']:.2f}, {combo['validation_region_max']:.2f})"
            )
            ax.grid(alpha=0.25)
        axes[-1].set_xlabel("Combined kinematic-bin number")
        axes[-1].set_xticks(np.arange(1, module.NUMBER_OF_BINS + 1))
        fig.suptitle(f"{module.PERIOD_LABELS[period]} — disjoint low-$M_X^2$ transfer closure")
        fig.tight_layout()
        fig.savefig(output / f"low_region_closure_{period}.png", dpi=180)
        plt.close(fig)


def build_spectrum_comparison(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    charge_fractions: dict[str, dict[str, float]],
    spectrum_edges: np.ndarray,
    nominal_window: tuple[float, float],
) -> pd.DataFrame:
    centers = 0.5 * (spectrum_edges[:-1] + spectrum_edges[1:])
    rows: list[dict[str, Any]] = []
    for period in module.PERIODS:
        nh3_spectrum = counts[(period, "NH3")].spectrum.astype(float)
        carbon_spectrum = counts[(period, "C")].spectrum.astype(float)
        control_nh3 = np.sum(nh3_spectrum[:, (centers >= nominal_window[0]) & (centers < nominal_window[1])])
        control_c = np.sum(carbon_spectrum[:, (centers >= nominal_window[0]) & (centers < nominal_window[1])])
        alpha = float(ratio(control_nh3, control_c))

        # Apply the exact direct five-target expression independently in every
        # Mx2 histogram bin.  This is a descriptive background-shape diagnostic;
        # low-statistics bins are retained and flagged rather than clipped.
        target_counts = {
            target: counts[(period, target)].spectrum.astype(float)
            for target in module.TARGETS
        }
        f2_spectrum = module.method2_equation10_from_counts(target_counts, charge_fractions[period])
        five_target_background = nh3_spectrum * (1.0 - f2_spectrum)
        scaled_carbon = alpha * carbon_spectrum

        summed_nh3 = np.sum(nh3_spectrum, axis=0)
        summed_scaled_c = np.sum(scaled_carbon, axis=0)
        summed_five = np.sum(five_target_background, axis=0)
        for i, center in enumerate(centers):
            rows.append({
                "period": period,
                "mx2_low": spectrum_edges[i],
                "mx2_high": spectrum_edges[i + 1],
                "mx2_center": center,
                "nh3_count": summed_nh3[i],
                "scaled_carbon_background": summed_scaled_c[i],
                "five_target_background": summed_five[i],
                "method1_dilution_factor": float(1.0 - ratio(summed_scaled_c[i], summed_nh3[i])),
                "method2_dilution_factor": float(1.0 - ratio(summed_five[i], summed_nh3[i])),
                "scaled_carbon_over_five_target": float(ratio(summed_scaled_c[i], summed_five[i])),
                "period_carbon_scale": alpha,
                "normalization_window_min": nominal_window[0],
                "normalization_window_max": nominal_window[1],
            })
    return pd.DataFrame(rows)


def plot_spectrum_comparison(
    output: Path,
    table: pd.DataFrame,
    module: Any,
    nominal_window: tuple[float, float],
) -> None:
    """Compare Method-1 and Method-2 backgrounds and dilution versus Mx2.

    The displayed lower edge is the upper edge of the selected Method-1
    normalization window, while the upper edge remains fixed at 1.4 GeV^2.
    This keeps the comparison focused entirely above the normalization region.
    """
    plot_min = float(nominal_window[1])
    plot_max = 1.4
    if not np.isfinite(plot_min) or plot_min >= plot_max:
        raise ValueError(
            "The Method-1 normalization-window upper edge must be finite and "
            f"below {plot_max:.2f} GeV^2 to make the Mx2 comparison plot; "
            f"received {plot_min:.6g} GeV^2."
        )
    for period in module.PERIODS:
        selected = table[table.period == period].sort_values("mx2_center")
        fig, axes = plt.subplots(
            2, 1, figsize=(13, 8), sharex=True,
            gridspec_kw={"height_ratios": [2.0, 1.2]},
        )
        centers = selected.mx2_center.to_numpy()

        axes[0].step(centers, selected.nh3_count, where="mid", label="NH$_3$ data")
        axes[0].step(
            centers, selected.scaled_carbon_background, where="mid",
            label="Method 1: scaled-carbon background",
        )
        axes[0].step(
            centers, selected.five_target_background, where="mid",
            label="Method 2: direct five-target background",
        )
        axes[0].set_ylabel("Counts")
        axes[0].set_yscale("symlog", linthresh=1.0)
        axes[0].grid(alpha=0.25)
        axes[0].legend()
        axes[0].set_title(
            f"{module.PERIOD_LABELS[period]} — period-integrated Method-1/Method-2 comparison"
        )

        axes[1].plot(
            centers, selected.method1_dilution_factor, "o-",
            markersize=3, label="Method 1",
        )
        axes[1].plot(
            centers, selected.method2_dilution_factor, "s-",
            markersize=3, label="Method 2",
        )
        axes[1].axhline(0.0, linewidth=0.8)
        axes[1].set_ylabel("Dilution factor")
        axes[1].grid(alpha=0.25)
        axes[1].legend()

        axes[1].set_xlabel(r"$M_X^2$ (GeV$^2$)")

        for ax in axes:
            # The selected Method-1 normalization window may lie partly or
            # entirely outside this focused comparison range.  Shade only the
            # visible overlap.
            visible_low = max(nominal_window[0], plot_min)
            visible_high = min(nominal_window[1], plot_max)
            if visible_high > visible_low:
                ax.axvspan(visible_low, visible_high, alpha=0.12, label=None)
            ax.set_xlim(plot_min, plot_max)

        axes[0].text(
            0.01, 0.03,
            rf"Method-1 normalization window: ${nominal_window[0]:.2f} \leq M_X^2 < {nominal_window[1]:.2f}$ GeV$^2$",
            transform=axes[0].transAxes, va="bottom", ha="left",
        )
        fig.tight_layout()
        fig.savefig(output / f"method_comparison_vs_mx2_{period}.png", dpi=180)
        plt.close(fig)



def plot_cumulative_window_scan(
    output: Path,
    module: Any,
    windows: list[tuple[float, float]],
    cumulative_edges: list[float],
    central: dict[str, Any],
    bootstrap: dict[str, Any],
) -> None:
    cut = list(module.CUT_VARIATIONS).index("nominal")
    indices = [windows.index((0.0, upper)) for upper in cumulative_edges]
    fig, axes = plt.subplots(len(module.PERIODS), 3, figsize=(17, 11), sharex=True, sharey="col")
    for row, period in enumerate(module.PERIODS):
        alpha = central[period]["alpha"][PERIOD_GROUP][:, indices]
        f1 = central[period]["method1"][PERIOD_GROUP][:, indices, cut]
        f2 = central[period]["method2"][:, cut]
        alpha_mean = np.nanmean(alpha, axis=0)
        f1_mean = np.nanmean(f1, axis=0)
        delta_mean = np.nanmean(f1 - f2[:, np.newaxis], axis=0)

        alpha_rep = bootstrap[period]["alpha"][PERIOD_GROUP][:, :, indices]
        f1_rep = bootstrap[period]["method1"][PERIOD_GROUP][:, :, indices, cut]
        f2_rep = bootstrap[period]["method2"][:, :, cut]
        alpha_err = np.nanstd(np.nanmean(alpha_rep, axis=1), axis=0, ddof=1)
        f1_err = np.nanstd(np.nanmean(f1_rep, axis=1), axis=0, ddof=1)
        delta_err = np.nanstd(np.nanmean(f1_rep - f2_rep[:, :, np.newaxis], axis=1), axis=0, ddof=1)

        axes[row, 0].errorbar(cumulative_edges, alpha_mean, yerr=alpha_err, fmt="o-")
        axes[row, 1].errorbar(cumulative_edges, f1_mean, yerr=f1_err, fmt="o-")
        axes[row, 2].errorbar(cumulative_edges, delta_mean, yerr=delta_err, fmt="o-")
        axes[row, 2].axhline(0.0, linewidth=1)
        for ax in axes[row]:
            ax.axvline(0.4, linewidth=1, linestyle="--")
            ax.grid(alpha=0.25)
        axes[row, 0].set_ylabel(module.PERIOD_LABELS[period])

    axes[0, 0].set_title(r"Period carbon scale $\alpha_C(0,u)$")
    axes[0, 1].set_title(r"Mean Method-1 dilution $\langle f_1\rangle$")
    axes[0, 2].set_title(r"Mean difference $\langle f_1-f_2\rangle$")
    for ax in axes[-1]:
        ax.set_xlabel(r"Upper control edge $u$ in $M_X^2$ (GeV$^2$)")
    fig.suptitle("Dense cumulative carbon-normalization scan with lower edge fixed at zero")
    fig.tight_layout()
    fig.savefig(output / "cumulative_upper_edge_scan.png", dpi=180)
    plt.close(fig)


def build_required_transfer_table(
    module: Any,
    central: dict[str, Any],
    bootstrap: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    for period in module.PERIODS:
        for bin_index in range(module.NUMBER_OF_BINS):
            ix, it = divmod(bin_index, n_t)
            for cut_index, cut in enumerate(module.CUT_VARIATIONS):
                required_rep = bootstrap[period]["alpha_required"][:, bin_index, cut_index]
                ratio_rep = bootstrap[period]["transfer_ratio"][:, bin_index, cut_index]
                rows.append({
                    "period": period,
                    "bin_number": bin_index + 1,
                    "x_index": ix,
                    "tprime_index": it,
                    "cut": cut,
                    "alpha_required": central[period]["alpha_required"][bin_index, cut_index],
                    "alpha_required_stat_uncertainty": np.nanstd(required_rep, ddof=1),
                    "required_over_control_alpha": central[period]["transfer_ratio"][bin_index, cut_index],
                    "required_over_control_stat_uncertainty": np.nanstd(ratio_rep, ddof=1),
                    "valid_required_fraction": np.mean(np.isfinite(required_rep)),
                    "valid_ratio_fraction": np.mean(np.isfinite(ratio_rep)),
                })
        for cut_index, cut in enumerate(module.CUT_VARIATIONS):
            kappa_rep = bootstrap[period]["kappa_period"][:, cut_index]
            rows.append({
                "period": period,
                "bin_number": 0,
                "x_index": -1,
                "tprime_index": -1,
                "cut": cut,
                "alpha_required": np.nan,
                "alpha_required_stat_uncertainty": np.nan,
                "required_over_control_alpha": central[period]["kappa_period"][cut_index],
                "required_over_control_stat_uncertainty": np.nanstd(kappa_rep, ddof=1),
                "valid_required_fraction": np.nan,
                "valid_ratio_fraction": np.mean(np.isfinite(kappa_rep)),
            })
    return pd.DataFrame(rows)


def plot_required_transfer_factor(
    output: Path,
    table: pd.DataFrame,
    module: Any,
) -> None:
    selected = table[(table.cut == "nominal") & (table.bin_number > 0)]
    x = np.arange(1, module.NUMBER_OF_BINS + 1)
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(16, 11), sharex=True, sharey=True)
    for ax, period in zip(axes, module.PERIODS):
        rows = selected[selected.period == period].sort_values("bin_number")
        ax.errorbar(
            x,
            rows.required_over_control_alpha,
            yerr=rows.required_over_control_stat_uncertainty,
            fmt="o",
            capsize=2,
            label=r"$\alpha_C^{\rm required}/\alpha_C^{\rm control}$",
        )
        summary = table[
            (table.period == period) & (table.cut == "nominal") & (table.bin_number == 0)
        ].iloc[0]
        ax.axhline(summary.required_over_control_alpha, linewidth=1.5, label="Period-integrated transfer factor")
        ax.axhline(1.0, linewidth=1, linestyle="--")
        ax.set_ylabel(r"$R_\alpha$")
        ax.set_title(
            f"{module.PERIOD_LABELS[period]}: integrated $\\kappa$ = "
            f"{summary.required_over_control_alpha:.4f} ± "
            f"{summary.required_over_control_stat_uncertainty:.4f}"
        )
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(x)
    fig.suptitle("Signal-window carbon transfer required for Method 1 to reproduce Method 2")
    fig.tight_layout()
    fig.savefig(output / "required_carbon_transfer_factor.png", dpi=180)
    plt.close(fig)


def plot_required_transfer_kinematics(
    output: Path,
    table: pd.DataFrame,
    module: Any,
) -> None:
    selected = table[(table.cut == "nominal") & (table.bin_number > 0)]
    fig, axes = plt.subplots(len(module.PERIODS), 2, figsize=(15, 11), sharey="col")
    for row, period in enumerate(module.PERIODS):
        rows = selected[selected.period == period]
        by_x = rows.groupby("x_index").agg(
            value=("required_over_control_alpha", "mean"),
            error=("required_over_control_stat_uncertainty", lambda x: float(np.sqrt(np.nansum(np.asarray(x) ** 2)) / max(1, np.count_nonzero(np.isfinite(x))))),
        )
        by_t = rows.groupby("tprime_index").agg(
            value=("required_over_control_alpha", "mean"),
            error=("required_over_control_stat_uncertainty", lambda x: float(np.sqrt(np.nansum(np.asarray(x) ** 2)) / max(1, np.count_nonzero(np.isfinite(x))))),
        )
        axes[row, 0].errorbar(by_x.index + 1, by_x.value, yerr=by_x.error, fmt="o-")
        axes[row, 1].errorbar(by_t.index + 1, by_t.value, yerr=by_t.error, fmt="o-")
        for ax in axes[row]:
            ax.axhline(1.0, linewidth=1, linestyle="--")
            ax.grid(alpha=0.25)
        axes[row, 0].set_ylabel(module.PERIOD_LABELS[period])
    axes[0, 0].set_title(r"Mean $R_\alpha$ by $x_B$ row")
    axes[0, 1].set_title(r"Mean $R_\alpha$ by $-t'$ column")
    axes[-1, 0].set_xlabel(r"$x_B$ row index")
    axes[-1, 1].set_xlabel(r"$-t'$ column index")
    fig.suptitle("Kinematic structure of the required carbon transfer factor")
    fig.tight_layout()
    fig.savefig(output / "required_carbon_transfer_kinematic_dependence.png", dpi=180)
    plt.close(fig)



def weighted_linear_fit(x: np.ndarray, y: np.ndarray, sigma: np.ndarray) -> dict[str, float]:
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    x = x[finite]
    y = y[finite]
    sigma = sigma[finite]
    if x.size < 2:
        return {"intercept": math.nan, "slope": math.nan, "intercept_uncertainty": math.nan,
                "slope_uncertainty": math.nan, "chi2": math.nan, "ndf": 0}
    design = np.column_stack([np.ones_like(x), x - np.mean(x)])
    weight = 1.0 / sigma**2
    normal = design.T @ (weight[:, np.newaxis] * design)
    try:
        covariance = np.linalg.inv(normal)
    except np.linalg.LinAlgError:
        return {"intercept": math.nan, "slope": math.nan, "intercept_uncertainty": math.nan,
                "slope_uncertainty": math.nan, "chi2": math.nan, "ndf": 0}
    parameters = covariance @ (design.T @ (weight * y))
    residual = y - design @ parameters
    chi2 = float(np.sum((residual / sigma) ** 2))
    return {
        "intercept": float(parameters[0]),
        "slope": float(parameters[1]),
        "intercept_uncertainty": float(math.sqrt(max(covariance[0, 0], 0.0))),
        "slope_uncertainty": float(math.sqrt(max(covariance[1, 1], 0.0))),
        "chi2": chi2,
        "ndf": int(x.size - 2),
    }


def build_transfer_fit_summary(table: pd.DataFrame, module: Any) -> pd.DataFrame:
    selected = table[(table.cut == "nominal") & (table.bin_number > 0)].copy()
    rows: list[dict[str, Any]] = []
    for period in module.PERIODS:
        period_rows = selected[selected.period == period]
        y = period_rows.required_over_control_alpha.to_numpy(float)
        sigma = period_rows.required_over_control_stat_uncertainty.to_numpy(float)
        finite = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
        weights = np.zeros_like(sigma)
        weights[finite] = 1.0 / sigma[finite] ** 2
        constant = float(np.sum(weights * np.nan_to_num(y)) / np.sum(weights)) if np.sum(weights) > 0 else math.nan
        constant_error = float(math.sqrt(1.0 / np.sum(weights))) if np.sum(weights) > 0 else math.nan
        chi2 = float(np.sum(((y[finite] - constant) / sigma[finite]) ** 2)) if np.isfinite(constant) else math.nan
        xfit = weighted_linear_fit(period_rows.x_index.to_numpy(float), y, sigma)
        tfit = weighted_linear_fit(period_rows.tprime_index.to_numpy(float), y, sigma)
        rows.append({
            "scope": period,
            "constant_transfer_factor": constant,
            "constant_uncertainty": constant_error,
            "constant_chi2": chi2,
            "constant_ndf": int(np.count_nonzero(finite) - 1),
            "x_index_slope": xfit["slope"],
            "x_index_slope_uncertainty": xfit["slope_uncertainty"],
            "x_fit_chi2": xfit["chi2"],
            "x_fit_ndf": xfit["ndf"],
            "tprime_index_slope": tfit["slope"],
            "tprime_index_slope_uncertainty": tfit["slope_uncertainty"],
            "tprime_fit_chi2": tfit["chi2"],
            "tprime_fit_ndf": tfit["ndf"],
        })
    all_y = selected.required_over_control_alpha.to_numpy(float)
    all_sigma = selected.required_over_control_stat_uncertainty.to_numpy(float)
    finite = np.isfinite(all_y) & np.isfinite(all_sigma) & (all_sigma > 0)
    weights = np.zeros_like(all_sigma)
    weights[finite] = 1.0 / all_sigma[finite] ** 2
    constant = float(np.sum(weights * np.nan_to_num(all_y)) / np.sum(weights)) if np.sum(weights) > 0 else math.nan
    error = float(math.sqrt(1.0 / np.sum(weights))) if np.sum(weights) > 0 else math.nan
    chi2 = float(np.sum(((all_y[finite] - constant) / all_sigma[finite]) ** 2)) if np.isfinite(constant) else math.nan
    rows.append({
        "scope": "all-periods",
        "constant_transfer_factor": constant,
        "constant_uncertainty": error,
        "constant_chi2": chi2,
        "constant_ndf": int(np.count_nonzero(finite) - 1),
        "x_index_slope": math.nan,
        "x_index_slope_uncertainty": math.nan,
        "x_fit_chi2": math.nan,
        "x_fit_ndf": 0,
        "tprime_index_slope": math.nan,
        "tprime_index_slope_uncertainty": math.nan,
        "tprime_fit_chi2": math.nan,
        "tprime_fit_ndf": 0,
    })
    return pd.DataFrame(rows)


def build_cumulative_background_comparison(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    charge_fractions: dict[str, dict[str, float]],
    spectrum_edges: np.ndarray,
    cumulative_edges: list[float],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    centers = 0.5 * (spectrum_edges[:-1] + spectrum_edges[1:])
    nominal_mask = (centers >= 0.0) & (centers < 0.4)
    for period in module.PERIODS:
        period_spectra = {
            target: np.sum(counts[(period, target)].spectrum.astype(float), axis=0)
            for target in module.TARGETS
        }
        alpha = float(ratio(
            np.sum(period_spectra["NH3"][nominal_mask]),
            np.sum(period_spectra["C"][nominal_mask]),
        ))
        for upper in cumulative_edges:
            mask = (centers >= 0.0) & (centers < upper)
            cumulative = {
                target: np.asarray(np.sum(values[mask]), dtype=float)
                for target, values in period_spectra.items()
            }
            f2 = float(module.method2_equation10_from_counts(cumulative, charge_fractions[period]))
            five_background = cumulative["NH3"] * (1.0 - f2)
            carbon_background = alpha * cumulative["C"]
            rows.append({
                "period": period,
                "upper_edge": upper,
                "nh3_cumulative": cumulative["NH3"],
                "scaled_carbon_cumulative_background": carbon_background,
                "five_target_cumulative_background": five_background,
                "scaled_carbon_over_five_target": float(ratio(carbon_background, five_background)),
                "method2_cumulative_dilution": f2,
                "nominal_control_alpha": alpha,
            })
    return pd.DataFrame(rows)


def plot_cumulative_background_comparison(output: Path, table: pd.DataFrame, module: Any) -> None:
    fig, axes = plt.subplots(len(module.PERIODS), 2, figsize=(15, 11), sharex=True, sharey="col")
    for row, period in enumerate(module.PERIODS):
        selected = table[table.period == period]
        axes[row, 0].plot(selected.upper_edge, selected.scaled_carbon_cumulative_background, "o-", label="Scaled carbon")
        axes[row, 0].plot(selected.upper_edge, selected.five_target_cumulative_background, "s-", label="Direct five-target")
        axes[row, 1].plot(selected.upper_edge, selected.scaled_carbon_over_five_target, "o-")
        axes[row, 1].axhline(1.0, linewidth=1, linestyle="--")
        for ax in axes[row]:
            ax.axvline(0.4, linewidth=1, linestyle="--")
            ax.grid(alpha=0.25)
        axes[row, 0].set_ylabel(module.PERIOD_LABELS[period])
        axes[row, 0].legend()
    axes[0, 0].set_title("Cumulative inferred background from 0 to upper edge")
    axes[0, 1].set_title("Scaled carbon / direct five-target background")
    axes[-1, 0].set_xlabel(r"Upper edge in $M_X^2$ (GeV$^2$)")
    axes[-1, 1].set_xlabel(r"Upper edge in $M_X^2$ (GeV$^2$)")
    fig.suptitle("Stable cumulative background comparison; narrow-bin nonlinear ratios removed")
    fig.tight_layout()
    fig.savefig(output / "cumulative_background_comparison.png", dpi=180)
    plt.close(fig)



def _method2_with_target_scale(
    signal_counts: dict[str, np.ndarray],
    target: str,
    scale: float,
    charge_fractions_period: dict[str, float],
    module: Any,
) -> np.ndarray:
    shifted = {
        name: np.asarray(values, dtype=float) * (scale if name == target else 1.0)
        for name, values in signal_counts.items()
    }
    return module.method2_equation10_from_counts(shifted, charge_fractions_period)


def build_target_sensitivity_table(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    charge_fractions: dict[str, dict[str, float]],
    central: dict[str, Any],
    windows: list[tuple[float, float]],
) -> pd.DataFrame:
    """Finite-difference sensitivity of Method 2 to coherent target normalization.

    S_T = d f2 / d ln N_T is evaluated with symmetric +/-1% changes.  The table
    also reports the approximate percent target-normalization shift needed to
    move Method 2 to Method 1 under the local linear approximation.
    """
    rows: list[dict[str, Any]] = []
    nominal_window_index = windows.index(NOMINAL_CONTROL_WINDOW)
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    eps = SENSITIVITY_FRACTION
    for period in module.PERIODS:
        signals = {
            target: counts[(period, target)].signal.astype(float)
            for target in module.TARGETS
        }
        f2 = central[period]["method2"]
        f1 = central[period]["method1"][PERIOD_GROUP][:, nominal_window_index, :]
        for target in module.TARGETS:
            plus = _method2_with_target_scale(
                signals, target, 1.0 + eps, charge_fractions[period], module
            )
            minus = _method2_with_target_scale(
                signals, target, 1.0 - eps, charge_fractions[period], module
            )
            sensitivity = (plus - minus) / (2.0 * eps)
            required_linear_fraction = ratio(f1 - f2, sensitivity)
            for bin_index in range(module.NUMBER_OF_BINS):
                ix, it = divmod(bin_index, n_t)
                for cut_index, cut in enumerate(module.CUT_VARIATIONS):
                    rows.append({
                        "period": period,
                        "bin_number": bin_index + 1,
                        "x_index": ix,
                        "tprime_index": it,
                        "cut": cut,
                        "target": target,
                        "method1_value": f1[bin_index, cut_index],
                        "method2_value": f2[bin_index, cut_index],
                        "d_method2_d_ln_target_yield": sensitivity[bin_index, cut_index],
                        "method2_shift_for_plus_1_percent": plus[bin_index, cut_index] - f2[bin_index, cut_index],
                        "linear_required_fractional_target_shift": required_linear_fraction[bin_index, cut_index],
                        "linear_required_percent_target_shift": 100.0 * required_linear_fraction[bin_index, cut_index],
                    })
    return pd.DataFrame(rows)


def _solve_required_scale_one(
    signal_counts: dict[str, np.ndarray],
    target: str,
    desired: float,
    charge_fractions_period: dict[str, float],
    module: Any,
) -> float:
    if not np.isfinite(desired):
        return math.nan
    grid = np.linspace(REQUIRED_SCALE_MIN, REQUIRED_SCALE_MAX, REQUIRED_SCALE_GRID)
    values = []
    for scale in grid:
        f2 = float(_method2_with_target_scale(
            signal_counts, target, float(scale), charge_fractions_period, module
        ))
        values.append(f2 - desired)
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    for index in range(len(grid) - 1):
        if not finite[index] or not finite[index + 1]:
            continue
        if values[index] == 0:
            return float(grid[index])
        if values[index] * values[index + 1] > 0:
            continue
        lo, hi = float(grid[index]), float(grid[index + 1])
        flo, fhi = float(values[index]), float(values[index + 1])
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            fmid = float(_method2_with_target_scale(
                signal_counts, target, mid, charge_fractions_period, module
            )) - desired
            if not np.isfinite(fmid):
                break
            if abs(fmid) < 1e-12:
                return mid
            if flo * fmid <= 0:
                hi, fhi = mid, fmid
            else:
                lo, flo = mid, fmid
        return 0.5 * (lo + hi)
    return math.nan


def build_required_target_scale_table(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    charge_fractions: dict[str, dict[str, float]],
    central: dict[str, Any],
    windows: list[tuple[float, float]],
) -> pd.DataFrame:
    """Solve for a coherent target-yield factor that would force f2=f1."""
    rows: list[dict[str, Any]] = []
    nominal_window_index = windows.index(NOMINAL_CONTROL_WINDOW)
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    for period in module.PERIODS:
        full = {
            target: counts[(period, target)].signal.astype(float)
            for target in module.TARGETS
        }
        f1 = central[period]["method1"][PERIOD_GROUP][:, nominal_window_index, :]
        for target in module.TARGETS:
            for bin_index in range(module.NUMBER_OF_BINS):
                ix, it = divmod(bin_index, n_t)
                for cut_index, cut in enumerate(module.CUT_VARIATIONS):
                    scalar_counts = {
                        name: np.asarray(values[bin_index, cut_index], dtype=float)
                        for name, values in full.items()
                    }
                    solved = _solve_required_scale_one(
                        scalar_counts, target, float(f1[bin_index, cut_index]),
                        charge_fractions[period], module
                    )
                    rows.append({
                        "period": period,
                        "bin_number": bin_index + 1,
                        "x_index": ix,
                        "tprime_index": it,
                        "cut": cut,
                        "target": target,
                        "required_yield_scale": solved,
                        "required_percent_shift": 100.0 * (solved - 1.0) if np.isfinite(solved) else math.nan,
                        "solution_in_scan_range": bool(np.isfinite(solved)),
                    })
            # Yield-weighted period-integrated solution for each cut.
            for cut_index, cut in enumerate(module.CUT_VARIATIONS):
                scalar_counts = {
                    name: np.asarray(np.sum(values[:, cut_index]), dtype=float)
                    for name, values in full.items()
                }
                nh3 = scalar_counts["NH3"]
                carbon = scalar_counts["C"]
                alpha = central[period]["alpha"][PERIOD_GROUP][0, nominal_window_index]
                desired = float((nh3 - alpha * carbon) / nh3) if nh3 else math.nan
                solved = _solve_required_scale_one(
                    scalar_counts, target, desired, charge_fractions[period], module
                )
                rows.append({
                    "period": period,
                    "bin_number": 0,
                    "x_index": -1,
                    "tprime_index": -1,
                    "cut": cut,
                    "target": target,
                    "required_yield_scale": solved,
                    "required_percent_shift": 100.0 * (solved - 1.0) if np.isfinite(solved) else math.nan,
                    "solution_in_scan_range": bool(np.isfinite(solved)),
                })
    return pd.DataFrame(rows)


def plot_target_sensitivity(output: Path, table: pd.DataFrame, module: Any) -> None:
    selected = table[table.cut == "nominal"].copy()
    summary = selected.groupby(["period", "target"], as_index=False).agg(
        sensitivity=("d_method2_d_ln_target_yield", "mean"),
        rms=("d_method2_d_ln_target_yield", lambda x: float(np.sqrt(np.nanmean(np.asarray(x) ** 2)))),
    )
    targets = list(module.TARGETS)
    fig, axes = plt.subplots(1, len(module.PERIODS), figsize=(16, 5), sharey=True)
    all_values = summary.sensitivity.to_numpy(dtype=float)
    limit = 1.15 * np.nanmax(np.abs(all_values)) if np.any(np.isfinite(all_values)) else 1.0
    for ax, period in zip(axes, module.PERIODS):
        local = summary[summary.period == period].set_index("target").reindex(targets)
        ax.bar(np.arange(len(targets)), local.sensitivity)
        ax.axhline(0.0, linewidth=1)
        ax.set_xticks(np.arange(len(targets)), targets, rotation=35, ha="right")
        ax.set_ylim(-limit, limit)
        ax.set_title(module.PERIOD_LABELS[period])
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel(r"$\partial f_2/\partial\ln N_T$")
    fig.suptitle("Method-2 sensitivity to coherent target-yield normalization")
    fig.tight_layout()
    fig.savefig(output / "method2_target_normalization_sensitivity.png", dpi=180)
    plt.close(fig)


def plot_required_target_scales(output: Path, table: pd.DataFrame, module: Any) -> None:
    selected = table[(table.cut == "nominal") & (table.bin_number == 0)].copy()
    targets = list(module.TARGETS)
    fig, ax = plt.subplots(figsize=(12, 6))
    width = 0.24
    x = np.arange(len(targets))
    for index, period in enumerate(module.PERIODS):
        local = selected[selected.period == period].set_index("target").reindex(targets)
        ax.bar(x + (index - 1) * width, local.required_percent_shift, width, label=module.PERIOD_LABELS[period])
    ax.axhline(0.0, linewidth=1)
    ax.set_xticks(x, targets)
    ax.set_ylabel("Required coherent target-yield shift (%)")
    ax.set_title("Target normalization change required to force Method 2 to Method 1")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "required_target_normalization_shifts.png", dpi=180)
    plt.close(fig)


def _weighted_linear_solution(design: np.ndarray, y: np.ndarray, variance: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, int]:
    good = np.all(np.isfinite(design), axis=1) & np.isfinite(y) & np.isfinite(variance) & (variance > 0)
    x = design[good]
    yy = y[good]
    vv = variance[good]
    if len(yy) <= x.shape[1]:
        return np.full(x.shape[1], np.nan), np.full((x.shape[1], x.shape[1]), np.nan), math.nan, 0
    w = 1.0 / vv
    normal = x.T @ (w[:, None] * x)
    rhs = x.T @ (w * yy)
    try:
        covariance = np.linalg.inv(normal)
        parameters = covariance @ rhs
    except np.linalg.LinAlgError:
        covariance = np.linalg.pinv(normal)
        parameters = covariance @ rhs
    residual = yy - x @ parameters
    chi2 = float(np.sum(w * residual * residual))
    return parameters, covariance, chi2, len(yy) - x.shape[1]


def build_empirical_template_fit(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    spectrum_edges: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit period-integrated NH3 spectra with C and CH2-C empirical templates.

    beta is measured from CH2/C in the nominal low-side control region, making
    H_aux = CH2 - beta*C have zero integral there.  The absolute normalization
    of H_aux is arbitrary; its fitted coefficient absorbs it.  Four linear
    fit variants probe carbon fixing/freeing and the polynomial boundary.
    """
    centers = 0.5 * (spectrum_edges[:-1] + spectrum_edges[1:])
    fit_mask = (centers >= TEMPLATE_FIT_MIN) & (centers < TEMPLATE_FIT_MAX)
    control_mask = (centers >= 0.0) & (centers < 0.4)
    above = centers >= 0.4
    hinge = np.clip(centers - 0.4, 0.0, None)
    rows: list[dict[str, Any]] = []
    curves: list[dict[str, Any]] = []
    for period in module.PERIODS:
        spectra = {
            target: np.sum(counts[(period, target)].spectrum.astype(float), axis=0)
            for target in module.TARGETS
        }
        alpha_control = float(ratio(np.sum(spectra["NH3"][control_mask]), np.sum(spectra["C"][control_mask])))
        beta_ch2_c = float(ratio(np.sum(spectra["CH2"][control_mask]), np.sum(spectra["C"][control_mask])))
        hydrogen = spectra["CH2"] - beta_ch2_c * spectra["C"]
        y = spectra["NH3"]
        variance = np.maximum(y, 1.0)
        polynomial_anchored = np.column_stack([hinge ** order for order in range(1, TEMPLATE_POLYNOMIAL_ORDER + 1)])
        polynomial_free = np.column_stack([above.astype(float), polynomial_anchored])
        variants = {
            "fixed_C_plus_H": (False, np.empty((len(centers), 0))),
            "fixed_C_plus_H_plus_anchored_poly": (False, polynomial_anchored),
            "free_C_plus_H_plus_anchored_poly": (True, polynomial_anchored),
            "free_C_plus_H_plus_free_intercept_poly": (True, polynomial_free),
        }
        for variant, (free_carbon, poly) in variants.items():
            fixed = np.zeros_like(y)
            columns = []
            names = []
            if free_carbon:
                columns.append(spectra["C"])
                names.append("carbon_scale")
            else:
                fixed = alpha_control * spectra["C"]
            columns.append(hydrogen)
            names.append("hydrogen_template_scale")
            for index in range(poly.shape[1]):
                columns.append(poly[:, index])
                names.append(f"polynomial_coefficient_{index}")
            design = np.column_stack(columns)
            parameters, covariance, chi2, ndf = _weighted_linear_solution(
                design[fit_mask], (y - fixed)[fit_mask], variance[fit_mask]
            )
            prediction = fixed + design @ parameters
            parameter_map = {name: parameters[i] for i, name in enumerate(names)}
            carbon_scale = parameter_map.get("carbon_scale", alpha_control)
            carbon_scale_unc = (
                math.sqrt(max(covariance[names.index("carbon_scale"), names.index("carbon_scale")], 0.0))
                if "carbon_scale" in names and np.all(np.isfinite(covariance)) else 0.0
            )
            rows.append({
                "period": period,
                "variant": variant,
                "alpha_control": alpha_control,
                "beta_ch2_over_c_control": beta_ch2_c,
                "fitted_carbon_scale": carbon_scale,
                "fitted_carbon_scale_uncertainty": carbon_scale_unc,
                "fitted_over_control_carbon_scale": float(ratio(carbon_scale, alpha_control)),
                "hydrogen_template_scale": parameter_map.get("hydrogen_template_scale", math.nan),
                "chi2": chi2,
                "ndf": ndf,
                "chi2_per_ndf": chi2 / ndf if ndf > 0 else math.nan,
            })
            for index, center in enumerate(centers):
                curves.append({
                    "period": period,
                    "variant": variant,
                    "mx2_center": center,
                    "nh3": y[index],
                    "carbon_component": carbon_scale * spectra["C"][index],
                    "hydrogen_template": hydrogen[index],
                    "prediction": prediction[index],
                    "residual": y[index] - prediction[index],
                })
    return pd.DataFrame(rows), pd.DataFrame(curves)


def plot_empirical_template_fits(output: Path, summary: pd.DataFrame, curves: pd.DataFrame, module: Any) -> None:
    preferred = "free_C_plus_H_plus_free_intercept_poly"
    fig, axes = plt.subplots(len(module.PERIODS), 2, figsize=(15, 11), sharex=True, sharey="col")
    y_max = float(np.nanmax(curves.nh3)) * 1.08
    residual_limit = float(np.nanpercentile(np.abs(curves.residual), 99)) * 1.2
    for row, period in enumerate(module.PERIODS):
        local = curves[(curves.period == period) & (curves.variant == preferred)]
        axes[row, 0].step(local.mx2_center, local.nh3, where="mid", label="NH3")
        axes[row, 0].plot(local.mx2_center, local.prediction, label="Template fit")
        axes[row, 0].plot(local.mx2_center, local.carbon_component, label="Carbon component")
        axes[row, 1].axhline(0.0, linewidth=1)
        axes[row, 1].plot(local.mx2_center, local.residual)
        for ax in axes[row]:
            ax.axvline(0.4, linewidth=1, linestyle="--")
            ax.grid(alpha=0.25)
        axes[row, 0].set_ylim(0, y_max)
        axes[row, 1].set_ylim(-residual_limit, residual_limit)
        axes[row, 0].set_ylabel(module.PERIOD_LABELS[period])
        axes[row, 0].legend()
    axes[0, 0].set_title("NH3 empirical-template decomposition")
    axes[0, 1].set_title("NH3 - fitted model")
    axes[-1, 0].set_xlabel(r"$M_X^2$ (GeV$^2$)")
    axes[-1, 1].set_xlabel(r"$M_X^2$ (GeV$^2$)")
    fig.suptitle("CH2-minus-carbon hydrogen template and carbon-transfer fit")
    fig.tight_layout()
    fig.savefig(output / "empirical_hydrogen_template_fit.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12, 6))
    variants = list(summary.variant.drop_duplicates())
    x = np.arange(len(variants))
    width = 0.24
    for index, period in enumerate(module.PERIODS):
        local = summary[summary.period == period].set_index("variant").reindex(variants)
        ax.bar(x + (index - 1) * width, local.fitted_over_control_carbon_scale, width, label=module.PERIOD_LABELS[period])
    ax.axhline(1.0, linewidth=1, linestyle="--")
    ax.set_xticks(x, variants, rotation=20, ha="right")
    ax.set_ylabel("Fitted carbon scale / control-region scale")
    ax.set_title("Independent empirical-template carbon transfer factors")
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "empirical_template_carbon_transfer_summary.png", dpi=180)
    plt.close(fig)



def build_required_rcn_table(module: Any, central: dict[str, Any], bootstrap: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    for period in module.PERIODS:
        for cut_index, cut in enumerate(module.CUT_VARIATIONS):
            period_samples = bootstrap[period]["required_rcn_period"][:, cut_index]
            rows.append({
                "period": period,
                "bin_number": 0,
                "x_index": -1,
                "tprime_index": -1,
                "cut": cut,
                "required_r_cn": central[period]["required_rcn_period"][cut_index],
                "required_r_cn_stat_uncertainty": np.nanstd(period_samples, ddof=1),
                "valid_replica_fraction": np.mean(np.isfinite(period_samples)),
                "nominal_r_cn": NOMINAL_R_CN,
                "required_over_nominal_r_cn": float(ratio(central[period]["required_rcn_period"][cut_index], NOMINAL_R_CN)),
                "scope": "period-integrated",
            })
            for bin_index in range(module.NUMBER_OF_BINS):
                ix, it = divmod(bin_index, n_t)
                samples = bootstrap[period]["required_rcn"][:, bin_index, cut_index]
                value = central[period]["required_rcn"][bin_index, cut_index]
                rows.append({
                    "period": period,
                    "bin_number": bin_index + 1,
                    "x_index": ix,
                    "tprime_index": it,
                    "cut": cut,
                    "required_r_cn": value,
                    "required_r_cn_stat_uncertainty": np.nanstd(samples, ddof=1),
                    "valid_replica_fraction": np.mean(np.isfinite(samples)),
                    "nominal_r_cn": NOMINAL_R_CN,
                    "required_over_nominal_r_cn": float(ratio(value, NOMINAL_R_CN)),
                    "scope": "bin",
                })
    return pd.DataFrame(rows)


def build_rcn_fit_summary(table: pd.DataFrame, module: Any) -> pd.DataFrame:
    selected = table[(table.cut == "nominal") & (table.bin_number > 0)].copy()
    rows: list[dict[str, Any]] = []
    for scope, local in [(period, selected[selected.period == period]) for period in module.PERIODS] + [("all-periods", selected)]:
        y = local.required_r_cn.to_numpy(float)
        sigma = local.required_r_cn_stat_uncertainty.to_numpy(float)
        finite = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
        w = np.where(finite, 1.0 / sigma**2, 0.0)
        sw = np.sum(w)
        constant = float(np.sum(w * np.nan_to_num(y)) / sw) if sw > 0 else math.nan
        error = float(math.sqrt(1.0 / sw)) if sw > 0 else math.nan
        chi2 = float(np.sum(((y[finite] - constant) / sigma[finite]) ** 2)) if np.isfinite(constant) else math.nan
        record = {
            "scope": scope,
            "constant_r_cn": constant,
            "constant_uncertainty": error,
            "constant_chi2": chi2,
            "constant_ndf": int(np.count_nonzero(finite) - 1),
            "nominal_r_cn": NOMINAL_R_CN,
            "constant_over_nominal": float(ratio(constant, NOMINAL_R_CN)),
        }
        if scope != "all-periods":
            xfit = weighted_linear_fit(local.x_index.to_numpy(float), y, sigma)
            tfit = weighted_linear_fit(local.tprime_index.to_numpy(float), y, sigma)
            record.update({
                "x_index_slope": xfit["slope"],
                "x_index_slope_uncertainty": xfit["slope_uncertainty"],
                "x_fit_chi2": xfit["chi2"],
                "x_fit_ndf": xfit["ndf"],
                "tprime_index_slope": tfit["slope"],
                "tprime_index_slope_uncertainty": tfit["slope_uncertainty"],
                "tprime_fit_chi2": tfit["chi2"],
                "tprime_fit_ndf": tfit["ndf"],
            })
        rows.append(record)
    return pd.DataFrame(rows)


def _replica_weighted_constant(samples: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Fit one constant to every bootstrap replica using fixed bin weights."""
    values = np.asarray(samples, dtype=float)
    base_weights = np.asarray(weights, dtype=float)
    finite = np.isfinite(values) & np.isfinite(base_weights)[np.newaxis, :] & (base_weights[np.newaxis, :] > 0)
    effective_weights = np.where(finite, base_weights[np.newaxis, :], 0.0)
    denominator = np.sum(effective_weights, axis=1)
    numerator = np.sum(effective_weights * np.where(finite, values, 0.0), axis=1)
    return np.divide(
        numerator,
        denominator,
        out=np.full(values.shape[0], np.nan, dtype=float),
        where=denominator > 0,
    )


def _replica_weighted_linear(
    x: np.ndarray,
    samples: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit intercept and slope to every replica with fixed inverse-variance weights."""
    x_values = np.asarray(x, dtype=float)
    values = np.asarray(samples, dtype=float)
    base_weights = np.asarray(weights, dtype=float)
    finite = (
        np.isfinite(values)
        & np.isfinite(x_values)[np.newaxis, :]
        & np.isfinite(base_weights)[np.newaxis, :]
        & (base_weights[np.newaxis, :] > 0)
    )
    w = np.where(finite, base_weights[np.newaxis, :], 0.0)
    y = np.where(finite, values, 0.0)
    x_row = x_values[np.newaxis, :]

    sw = np.sum(w, axis=1)
    sx = np.sum(w * x_row, axis=1)
    sy = np.sum(w * y, axis=1)
    sxx = np.sum(w * x_row * x_row, axis=1)
    sxy = np.sum(w * x_row * y, axis=1)
    determinant = sw * sxx - sx * sx

    intercept = np.full(values.shape[0], np.nan, dtype=float)
    slope = np.full(values.shape[0], np.nan, dtype=float)
    valid = determinant > 0
    intercept[valid] = (
        sxx[valid] * sy[valid] - sx[valid] * sxy[valid]
    ) / determinant[valid]
    slope[valid] = (
        sw[valid] * sxy[valid] - sx[valid] * sy[valid]
    ) / determinant[valid]
    return intercept, slope


def _distribution_summary(samples: np.ndarray) -> dict[str, float]:
    values = np.asarray(samples, dtype=float)
    finite = np.isfinite(values)
    if np.count_nonzero(finite) < 2:
        return {
            "bootstrap_mean": math.nan,
            "bootstrap_std": math.nan,
            "bootstrap_p16": math.nan,
            "bootstrap_median": math.nan,
            "bootstrap_p84": math.nan,
            "valid_replica_fraction": float(np.mean(finite)),
        }
    selected = values[finite]
    return {
        "bootstrap_mean": float(np.mean(selected)),
        "bootstrap_std": float(np.std(selected, ddof=1)),
        "bootstrap_p16": float(np.percentile(selected, 16)),
        "bootstrap_median": float(np.percentile(selected, 50)),
        "bootstrap_p84": float(np.percentile(selected, 84)),
        "valid_replica_fraction": float(np.mean(finite)),
    }


def build_rcn_replica_fit_summary(
    table: pd.DataFrame,
    bootstrap: dict[str, Any],
    module: Any,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Fit the required r_CN replica by replica.

    Fixed inverse-variance weights are obtained from the ensemble standard
    deviations in required_r_cn.csv.  Applying those same weights to every
    correlated replica preserves the bin-to-bin covariance in the fitted
    constant and slopes.  The spread of the replica-level fit parameters is
    therefore the quoted statistical uncertainty.
    """
    nominal_cut = list(module.CUT_VARIATIONS).index("nominal")
    selected = table[(table.cut == "nominal") & (table.bin_number > 0)].copy()
    summary_rows: list[dict[str, Any]] = []
    replica_rows: list[dict[str, Any]] = []
    period_payload: dict[str, dict[str, Any]] = {}

    for period in module.PERIODS:
        local = selected[selected.period == period].sort_values("bin_number")
        central_y = local.required_r_cn.to_numpy(float)
        sigma = local.required_r_cn_stat_uncertainty.to_numpy(float)
        weights = np.where(np.isfinite(sigma) & (sigma > 0), 1.0 / sigma**2, 0.0)
        samples = np.asarray(
            bootstrap[period]["required_rcn"][:, :, nominal_cut],
            dtype=float,
        )
        x_index = local.x_index.to_numpy(float)
        t_index = local.tprime_index.to_numpy(float)

        constant_samples = _replica_weighted_constant(samples, weights)
        x_intercept_samples, x_slope_samples = _replica_weighted_linear(
            x_index, samples, weights
        )
        t_intercept_samples, t_slope_samples = _replica_weighted_linear(
            t_index, samples, weights
        )

        central_constant = _replica_weighted_constant(
            central_y[np.newaxis, :], weights
        )[0]
        central_x = weighted_linear_fit(x_index, central_y, sigma)
        central_t = weighted_linear_fit(t_index, central_y, sigma)

        constant_stats = _distribution_summary(constant_samples)
        x_slope_stats = _distribution_summary(x_slope_samples)
        t_slope_stats = _distribution_summary(t_slope_samples)
        x_intercept_stats = _distribution_summary(x_intercept_samples)
        t_intercept_stats = _distribution_summary(t_intercept_samples)

        finite = np.isfinite(central_y) & np.isfinite(sigma) & (sigma > 0)
        constant_chi2 = (
            float(np.sum(((central_y[finite] - central_constant) / sigma[finite]) ** 2))
            if np.isfinite(central_constant)
            else math.nan
        )
        summary_rows.append({
            "scope": period,
            "fit_type": "period",
            "constant_r_cn": central_constant,
            "constant_stat_uncertainty": constant_stats["bootstrap_std"],
            "constant_bootstrap_mean": constant_stats["bootstrap_mean"],
            "constant_bootstrap_median": constant_stats["bootstrap_median"],
            "constant_bootstrap_p16": constant_stats["bootstrap_p16"],
            "constant_bootstrap_p84": constant_stats["bootstrap_p84"],
            "constant_valid_replica_fraction": constant_stats["valid_replica_fraction"],
            "constant_chi2_using_marginal_errors": constant_chi2,
            "constant_ndf": int(np.count_nonzero(finite) - 1),
            "nominal_r_cn": NOMINAL_R_CN,
            "constant_over_nominal": float(ratio(central_constant, NOMINAL_R_CN)),
            "x_index_intercept": central_x["intercept"],
            "x_index_intercept_stat_uncertainty": x_intercept_stats["bootstrap_std"],
            "x_index_slope": central_x["slope"],
            "x_index_slope_stat_uncertainty": x_slope_stats["bootstrap_std"],
            "x_index_slope_bootstrap_p16": x_slope_stats["bootstrap_p16"],
            "x_index_slope_bootstrap_p84": x_slope_stats["bootstrap_p84"],
            "x_index_slope_valid_replica_fraction": x_slope_stats["valid_replica_fraction"],
            "x_fit_chi2_using_marginal_errors": central_x["chi2"],
            "x_fit_ndf": central_x["ndf"],
            "tprime_index_intercept": central_t["intercept"],
            "tprime_index_intercept_stat_uncertainty": t_intercept_stats["bootstrap_std"],
            "tprime_index_slope": central_t["slope"],
            "tprime_index_slope_stat_uncertainty": t_slope_stats["bootstrap_std"],
            "tprime_index_slope_bootstrap_p16": t_slope_stats["bootstrap_p16"],
            "tprime_index_slope_bootstrap_p84": t_slope_stats["bootstrap_p84"],
            "tprime_index_slope_valid_replica_fraction": t_slope_stats["valid_replica_fraction"],
            "tprime_fit_chi2_using_marginal_errors": central_t["chi2"],
            "tprime_fit_ndf": central_t["ndf"],
        })

        period_payload[period] = {
            "central": central_y,
            "sigma": sigma,
            "weights": weights,
            "samples": samples,
            "x_index": x_index,
            "t_index": t_index,
            "constant_samples": constant_samples,
            "x_slope_samples": x_slope_samples,
            "t_slope_samples": t_slope_samples,
        }
        for replica_index in range(samples.shape[0]):
            replica_rows.append({
                "scope": period,
                "replica": replica_index,
                "constant_r_cn": constant_samples[replica_index],
                "x_index_intercept": x_intercept_samples[replica_index],
                "x_index_slope": x_slope_samples[replica_index],
                "tprime_index_intercept": t_intercept_samples[replica_index],
                "tprime_index_slope": t_slope_samples[replica_index],
            })

    # Pair replica indices across periods.  The period bootstraps are
    # statistically independent, so this constructs the proper joint
    # statistical ensemble for a common all-period fit.
    replica_count = min(
        payload["samples"].shape[0] for payload in period_payload.values()
    )
    all_central = np.concatenate(
        [period_payload[period]["central"] for period in module.PERIODS]
    )
    all_sigma = np.concatenate(
        [period_payload[period]["sigma"] for period in module.PERIODS]
    )
    all_weights = np.where(
        np.isfinite(all_sigma) & (all_sigma > 0),
        1.0 / all_sigma**2,
        0.0,
    )
    all_samples = np.concatenate(
        [
            period_payload[period]["samples"][:replica_count]
            for period in module.PERIODS
        ],
        axis=1,
    )
    all_constant_samples = _replica_weighted_constant(all_samples, all_weights)
    all_central_constant = _replica_weighted_constant(
        all_central[np.newaxis, :], all_weights
    )[0]
    all_stats = _distribution_summary(all_constant_samples)
    all_finite = np.isfinite(all_central) & np.isfinite(all_sigma) & (all_sigma > 0)
    all_chi2 = (
        float(np.sum(
            ((all_central[all_finite] - all_central_constant) / all_sigma[all_finite]) ** 2
        ))
        if np.isfinite(all_central_constant)
        else math.nan
    )
    summary_rows.append({
        "scope": "all-periods",
        "fit_type": "common",
        "constant_r_cn": all_central_constant,
        "constant_stat_uncertainty": all_stats["bootstrap_std"],
        "constant_bootstrap_mean": all_stats["bootstrap_mean"],
        "constant_bootstrap_median": all_stats["bootstrap_median"],
        "constant_bootstrap_p16": all_stats["bootstrap_p16"],
        "constant_bootstrap_p84": all_stats["bootstrap_p84"],
        "constant_valid_replica_fraction": all_stats["valid_replica_fraction"],
        "constant_chi2_using_marginal_errors": all_chi2,
        "constant_ndf": int(np.count_nonzero(all_finite) - 1),
        "nominal_r_cn": NOMINAL_R_CN,
        "constant_over_nominal": float(ratio(all_central_constant, NOMINAL_R_CN)),
        "x_index_intercept": math.nan,
        "x_index_intercept_stat_uncertainty": math.nan,
        "x_index_slope": math.nan,
        "x_index_slope_stat_uncertainty": math.nan,
        "x_index_slope_bootstrap_p16": math.nan,
        "x_index_slope_bootstrap_p84": math.nan,
        "x_index_slope_valid_replica_fraction": math.nan,
        "x_fit_chi2_using_marginal_errors": math.nan,
        "x_fit_ndf": 0,
        "tprime_index_intercept": math.nan,
        "tprime_index_intercept_stat_uncertainty": math.nan,
        "tprime_index_slope": math.nan,
        "tprime_index_slope_stat_uncertainty": math.nan,
        "tprime_index_slope_bootstrap_p16": math.nan,
        "tprime_index_slope_bootstrap_p84": math.nan,
        "tprime_index_slope_valid_replica_fraction": math.nan,
        "tprime_fit_chi2_using_marginal_errors": math.nan,
        "tprime_fit_ndf": 0,
    })
    for replica_index in range(replica_count):
        replica_rows.append({
            "scope": "all-periods",
            "replica": replica_index,
            "constant_r_cn": all_constant_samples[replica_index],
            "x_index_intercept": math.nan,
            "x_index_slope": math.nan,
            "tprime_index_intercept": math.nan,
            "tprime_index_slope": math.nan,
        })

    return pd.DataFrame(summary_rows), pd.DataFrame(replica_rows)


def plot_rcn_replica_fit_summary(
    output: Path,
    summary: pd.DataFrame,
    module: Any,
) -> None:
    """Plot covariance-aware constant and slope results."""
    period_rows = summary[summary.scope.isin(module.PERIODS)].set_index("scope")
    common = summary[summary.scope == "all-periods"].iloc[0]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    x = np.arange(len(module.PERIODS))
    labels = [module.PERIOD_LABELS[period] for period in module.PERIODS]

    constants = period_rows.loc[list(module.PERIODS), "constant_r_cn"].to_numpy(float)
    constant_errors = period_rows.loc[
        list(module.PERIODS), "constant_stat_uncertainty"
    ].to_numpy(float)
    axes[0].errorbar(x, constants, yerr=constant_errors, fmt="o", label="Period fits")
    axes[0].axhline(common.constant_r_cn, linewidth=1, label="Common fit")
    axes[0].axhspan(
        common.constant_r_cn - common.constant_stat_uncertainty,
        common.constant_r_cn + common.constant_stat_uncertainty,
        alpha=0.15,
    )
    axes[0].axhline(NOMINAL_R_CN, linewidth=1, linestyle="--", label=r"Nominal $6/7$")
    axes[0].set_xticks(x, labels, rotation=20)
    axes[0].set_ylabel(r"Required $r_{CN}$")
    axes[0].set_title("Replica-level constant fits")
    axes[0].legend()

    x_slopes = period_rows.loc[list(module.PERIODS), "x_index_slope"].to_numpy(float)
    x_errors = period_rows.loc[
        list(module.PERIODS), "x_index_slope_stat_uncertainty"
    ].to_numpy(float)
    axes[1].errorbar(x, x_slopes, yerr=x_errors, fmt="o")
    axes[1].axhline(0.0, linewidth=1, linestyle="--")
    axes[1].set_xticks(x, labels, rotation=20)
    axes[1].set_ylabel(r"Slope in $r_{CN}$ per $x_B$ row")
    axes[1].set_title(r"$x_B$-row dependence")

    t_slopes = period_rows.loc[
        list(module.PERIODS), "tprime_index_slope"
    ].to_numpy(float)
    t_errors = period_rows.loc[
        list(module.PERIODS), "tprime_index_slope_stat_uncertainty"
    ].to_numpy(float)
    axes[2].errorbar(x, t_slopes, yerr=t_errors, fmt="o")
    axes[2].axhline(0.0, linewidth=1, linestyle="--")
    axes[2].set_xticks(x, labels, rotation=20)
    axes[2].set_ylabel(r"Slope in $r_{CN}$ per $-t'$ column")
    axes[2].set_title(r"$-t'$-column dependence")

    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle(r"Replica-by-replica fits of the required C--N relation")
    fig.tight_layout()
    fig.savefig(output / "required_r_cn_replica_fit_summary.png", dpi=180)
    plt.close(fig)


def build_rcn_scan_table(
    counts: dict[tuple[str, str], ValidationCounts],
    module: Any,
    charge_fractions: dict[str, dict[str, float]],
    central: dict[str, Any],
    r_min: float,
    r_max: float,
    points: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    grid = np.linspace(r_min, r_max, points)
    nominal_cut = list(module.CUT_VARIATIONS).index("nominal")
    for period in module.PERIODS:
        integrated = {
            target: np.sum(counts[(period, target)].signal.astype(float), axis=0)
            for target in module.TARGETS
        }
        target_f1 = central[period]["integrated_method1"][nominal_cut]
        for value in grid:
            f2 = float(method2_with_rcn(
                {target: values[nominal_cut] for target, values in integrated.items()},
                charge_fractions[period], module, value,
            ))
            rows.append({
                "period": period,
                "r_cn": value,
                "method2": f2,
                "method1_reference": target_f1,
                "method2_minus_method1": f2 - target_f1,
            })
    return pd.DataFrame(rows)


def plot_required_rcn(output: Path, table: pd.DataFrame, module: Any) -> None:
    nominal = table[(table.cut == "nominal") & (table.bin_number > 0)]
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(14, 11), sharex=True, sharey=True)
    labels = flatten_bin_labels(module)
    for ax, period in zip(axes, module.PERIODS):
        local = nominal[nominal.period == period].sort_values("bin_number")
        ax.errorbar(local.bin_number, local.required_r_cn, yerr=local.required_r_cn_stat_uncertainty, fmt="o")
        integrated = table[(table.period == period) & (table.cut == "nominal") & (table.bin_number == 0)].iloc[0]
        ax.axhline(NOMINAL_R_CN, linewidth=1, linestyle="--", label=r"Nominal $6/7$")
        ax.axhline(integrated.required_r_cn, linewidth=1, label="Period-integrated required value")
        ax.set_ylabel(r"Required $r_{CN}$")
        ax.set_title(module.PERIOD_LABELS[period])
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel(r"Kinematic bin ($x_B$ row:$-t'$ column)")
    axes[-1].set_xticks(np.arange(1, module.NUMBER_OF_BINS + 1), labels, rotation=90)
    fig.suptitle(r"Carbon-to-nitrogen ratio required for $f_2(r_{CN})=f_1$")
    fig.tight_layout()
    fig.savefig(output / "required_r_cn_by_bin.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(len(module.PERIODS), 2, figsize=(14, 11), sharex="col", sharey=True)
    for row, period in enumerate(module.PERIODS):
        local = nominal[nominal.period == period]
        for column, key in enumerate(("x_index", "tprime_index")):
            grouped = local.groupby(key, as_index=False).apply(
                lambda frame: pd.Series({
                    "value": np.average(frame.required_r_cn, weights=1.0 / frame.required_r_cn_stat_uncertainty**2),
                    "error": math.sqrt(1.0 / np.sum(1.0 / frame.required_r_cn_stat_uncertainty**2)),
                }), include_groups=False
            ).reset_index(drop=True)
            axes[row, column].errorbar(grouped[key] + 1, grouped.value, yerr=grouped.error, fmt="o-")
            axes[row, column].axhline(NOMINAL_R_CN, linewidth=1, linestyle="--")
            axes[row, column].grid(alpha=0.25)
        axes[row, 0].set_ylabel(module.PERIOD_LABELS[period] + "\n" + r"Required $r_{CN}$")
    axes[0, 0].set_title(r"Dependence on $x_B$ row")
    axes[0, 1].set_title(r"Dependence on $-t'$ column")
    axes[-1, 0].set_xlabel(r"$x_B$ row index")
    axes[-1, 1].set_xlabel(r"$-t'$ column index")
    fig.suptitle(r"Kinematic dependence of the required C--N relation")
    fig.tight_layout()
    fig.savefig(output / "required_r_cn_kinematic_dependence.png", dpi=180)
    plt.close(fig)


def plot_rcn_scan(output: Path, scan: pd.DataFrame, module: Any) -> None:
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(11, 10), sharex=True, sharey=True)
    for ax, period in zip(axes, module.PERIODS):
        local = scan[scan.period == period]
        ax.plot(local.r_cn, local.method2, label=r"Method 2, $f_2(r_{CN})$")
        ax.plot(local.r_cn, local.method1_reference, linestyle="--", label="Method 1 reference")
        ax.axvline(NOMINAL_R_CN, linewidth=1, linestyle=":", label=r"Nominal $r_{CN}=6/7$")
        ax.set_ylabel("Dilution factor")
        ax.set_title(module.PERIOD_LABELS[period])
        ax.grid(alpha=0.25)
        ax.legend()
    axes[-1].set_xlabel(r"Assumed $r_{CN}=\sigma_C/\sigma_N$")
    fig.suptitle(r"Period-integrated Method 2 response to the C--N relation")
    fig.tight_layout()
    fig.savefig(output / "method2_r_cn_scan.png", dpi=180)
    plt.close(fig)

def build_flat_results(
    module: Any,
    windows: list[tuple[float, float]],
    central: dict[str, Any],
    bootstrap: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    for period in module.PERIODS:
        m2_summary = summarize(bootstrap[period]["method2"])
        m1_summaries = {
            grouping: summarize(bootstrap[period]["method1"][grouping])
            for grouping in GROUPINGS
        }
        for bin_index in range(module.NUMBER_OF_BINS):
            ix, it = divmod(bin_index, n_t)
            for cut_index, cut in enumerate(module.CUT_VARIATIONS):
                for window_index, window in enumerate(windows):
                    for grouping in GROUPINGS:
                        m1 = central[period]["method1"][grouping][bin_index, window_index, cut_index]
                        m2 = central[period]["method2"][bin_index, cut_index]
                        rep_delta = (
                            bootstrap[period]["method1"][grouping][:, bin_index, window_index, cut_index]
                            - bootstrap[period]["method2"][:, bin_index, cut_index]
                        )
                        rows.append({
                            "period": period,
                            "bin_number": bin_index + 1,
                            "x_index": ix,
                            "tprime_index": it,
                            "cut": cut,
                            "grouping": grouping,
                            "control_min": window[0],
                            "control_max": window[1],
                            "method1_value": m1,
                            "method1_stat_uncertainty": m1_summaries[grouping]["std"][bin_index, window_index, cut_index],
                            "method2_value": m2,
                            "method2_stat_uncertainty": m2_summary["std"][bin_index, cut_index],
                            "method1_minus_method2": m1 - m2,
                            "difference_stat_uncertainty": np.nanstd(rep_delta, ddof=1),
                            "alpha": central[period]["alpha"][grouping][bin_index, window_index],
                            "alpha_stat_uncertainty": np.nanstd(bootstrap[period]["alpha"][grouping][:, bin_index, window_index], ddof=1),
                            "alpha_required_signal": central[period]["alpha_required"][bin_index, cut_index],
                            "alpha_required_signal_stat_uncertainty": np.nanstd(bootstrap[period]["alpha_required"][:, bin_index, cut_index], ddof=1),
                            "required_over_nominal_period_alpha": central[period]["transfer_ratio"][bin_index, cut_index],
                            "required_over_nominal_period_alpha_stat_uncertainty": np.nanstd(bootstrap[period]["transfer_ratio"][:, bin_index, cut_index], ddof=1),
                        })
    return pd.DataFrame(rows)


def write_summary_json(
    path: Path,
    module: Any,
    windows: list[tuple[float, float]],
    flat: pd.DataFrame,
) -> None:
    nominal = flat[
        (flat.cut == "nominal")
        & np.isclose(flat.control_min, 0.0)
        & np.isclose(flat.control_max, 0.4)
    ]
    payload: dict[str, Any] = {
        "purpose": "Carbon-normalization validation; diagnostic only",
        "control_windows": windows,
        "method1_nominal_normalization_window": list(NOMINAL_CONTROL_WINDOW),
        "high_side_closure_policy": (
            "No high-side NH3/carbon closure requirement is imposed because "
            "NH3 can contain genuine hydrogen-origin nonexclusive channels."
        ),
        "periods": {},
    }
    for period in module.PERIODS:
        payload["periods"][period] = {}
        for grouping in GROUPINGS:
            selected = nominal[(nominal.period == period) & (nominal.grouping == grouping)]
            payload["periods"][period][grouping] = {
                "mean_method1": selected.method1_value.mean(),
                "mean_method2": selected.method2_value.mean(),
                "mean_method1_minus_method2": selected.method1_minus_method2.mean(),
                "rms_method1_minus_method2": float(np.sqrt(np.mean(selected.method1_minus_method2.to_numpy() ** 2))),
                "maximum_abs_method1_minus_method2": selected.method1_minus_method2.abs().max(),
            }
    write_json(path, json_safe(payload))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--dilution-module", type=Path, default=DEFAULT_MODULE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--tree", default=None, help="Override production module's default ROOT tree name")
    parser.add_argument("--workers", type=int, default=7)
    parser.add_argument("--replicas", type=int, default=DEFAULT_REPLICAS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--cut-json", type=Path, default=None)
    parser.add_argument("--charge-fractions-json", type=Path, default=None)
    parser.add_argument("--control-window", action="append", type=parse_window, default=None,
                        help="Control window MIN:MAX; repeat to replace the default scan")
    parser.add_argument(
        "--normalization-window-min", type=float,
        default=DEFAULT_NOMINAL_CONTROL_WINDOW[0],
        help="Lower edge of the production-like Method-1 carbon normalization window",
    )
    parser.add_argument(
        "--normalization-window-max", type=float,
        default=DEFAULT_NOMINAL_CONTROL_WINDOW[1],
        help="Upper edge of the production-like Method-1 carbon normalization window",
    )
    parser.add_argument("--spectrum-min", type=float, default=DEFAULT_SPECTRUM_RANGE[0])
    parser.add_argument("--spectrum-max", type=float, default=DEFAULT_SPECTRUM_RANGE[1])
    parser.add_argument("--spectrum-bins", type=int, default=DEFAULT_SPECTRUM_BINS)
    parser.add_argument("--cumulative-upper-min", type=float, default=DEFAULT_CUMULATIVE_UPPER_MIN)
    parser.add_argument("--cumulative-upper-max", type=float, default=DEFAULT_CUMULATIVE_UPPER_MAX)
    parser.add_argument("--cumulative-upper-step", type=float, default=DEFAULT_CUMULATIVE_UPPER_STEP)
    parser.add_argument("--input", action="append", default=[],
                        help="Override standard input as period:target=/path.root")
    parser.add_argument("--template-fit-min", type=float, default=TEMPLATE_FIT_MIN)
    parser.add_argument("--template-fit-max", type=float, default=TEMPLATE_FIT_MAX)
    parser.add_argument("--r-cn-min", type=float, default=DEFAULT_R_CN_MIN)
    parser.add_argument("--r-cn-max", type=float, default=DEFAULT_R_CN_MAX)
    parser.add_argument("--r-cn-points", type=int, default=DEFAULT_R_CN_POINTS)
    parser.add_argument("--run-isr", action="store_true",
                        help="Repeat the complete validation using ISR inputs and ISR cut JSON")
    parser.add_argument("--isr-cut-json", type=Path, default=None)
    parser.add_argument("--isr-input", action="append", default=[],
                        help="Override ISR input as period:target=/path.root")
    return parser


def resolve_inputs(module: Any, overrides: list[str], isr: bool) -> dict[str, dict[str, Path]]:
    base = module.DEFAULT_ISR_INPUTS if isr else module.DEFAULT_INPUTS
    output = {period: {target: Path(path) for target, path in targets.items()} for period, targets in base.items()}
    for text in overrides:
        period, target, path = module.parse_input_override(text)
        output[period][target] = path
    return output


def resolve_cut_json(module: Any, args: argparse.Namespace, isr: bool) -> Path:
    if isr:
        if args.isr_cut_json is not None:
            return args.isr_cut_json
        return module.resolve_isr_cut_json(module.DEFAULT_CHANNEL_SELECTION_MANIFEST, None)
    return args.cut_json if args.cut_json is not None else module.DEFAULT_CUT_JSON


def run_one_variant(
    label: str,
    module: Any,
    module_path: Path,
    input_paths: dict[str, dict[str, Path]],
    cut_json: Path,
    charge_fractions: dict[str, dict[str, float]],
    args: argparse.Namespace,
    windows: list[tuple[float, float]],
    cumulative_edges: list[float],
) -> None:
    root = args.output_dir / label
    tables = root / "tables"
    plots = root / "plots"
    metadata = root / "metadata"
    for path in (tables, plots, metadata):
        ensure_dir(path)

    tree_name = args.tree or module.DEFAULT_TREE_NAME
    workers = max(1, min(MAXIMUM_WORKERS, int(args.workers)))
    print(f"\n[{label}] Loading inputs with {workers} workers...")
    loaded = module.load_all_datasets(input_paths, tree_name, workers)
    cuts = module.load_exclusivity_cuts(cut_json)

    atomic_edges = build_atomic_edges(windows)
    spectrum_edges = np.linspace(args.spectrum_min, args.spectrum_max, args.spectrum_bins + 1)
    counts = build_all_counts(
        loaded, cuts, module, atomic_edges, DEFAULT_LOW_REGIONS, spectrum_edges, workers
    )
    central = central_estimators(counts, module, windows, atomic_edges, charge_fractions)
    bootstrap = run_bootstrap(
        counts, module, module_path, windows, atomic_edges, charge_fractions,
        args.replicas, args.seed + (1000003 if label == "isr" else 0), workers
    )

    flat = build_flat_results(module, windows, central, bootstrap)
    flat.to_csv(tables / "carbon_normalization_validation.csv", index=False)

    closure = low_region_closure_table(counts, module, DEFAULT_LOW_REGIONS)
    closure.to_csv(tables / "low_mx2_transfer_closure.csv", index=False)

    transfer = build_required_transfer_table(module, central, bootstrap)
    transfer.to_csv(tables / "required_carbon_transfer_factor.csv", index=False)
    transfer_fit = build_transfer_fit_summary(transfer, module)
    transfer_fit.to_csv(tables / "required_carbon_transfer_fit_summary.csv", index=False)

    required_rcn = build_required_rcn_table(module, central, bootstrap)
    required_rcn.to_csv(tables / "required_r_cn.csv", index=False)
    rcn_fit = build_rcn_fit_summary(required_rcn, module)
    rcn_fit.to_csv(tables / "required_r_cn_fit_summary_legacy_independent_bins.csv", index=False)
    rcn_replica_fit, rcn_replica_parameters = build_rcn_replica_fit_summary(
        required_rcn, bootstrap, module
    )
    rcn_replica_fit.to_csv(
        tables / "required_r_cn_fit_summary.csv", index=False
    )
    rcn_replica_parameters.to_csv(
        tables / "required_r_cn_replica_fit_parameters.csv", index=False
    )
    rcn_scan = build_rcn_scan_table(
        counts, module, charge_fractions, central,
        args.r_cn_min, args.r_cn_max, args.r_cn_points,
    )
    rcn_scan.to_csv(tables / "method2_r_cn_scan.csv", index=False)

    spectrum_comparison = build_spectrum_comparison(
        counts, module, charge_fractions, spectrum_edges, NOMINAL_CONTROL_WINDOW
    )
    spectrum_comparison.to_csv(
        tables / "method_comparison_vs_mx2.csv", index=False
    )

    cumulative_background = build_cumulative_background_comparison(
        counts, module, charge_fractions, spectrum_edges, cumulative_edges
    )
    cumulative_background.to_csv(tables / "cumulative_background_comparison.csv", index=False)

    target_sensitivity = build_target_sensitivity_table(
        counts, module, charge_fractions, central, windows
    )
    target_sensitivity.to_csv(tables / "method2_target_normalization_sensitivity.csv", index=False)

    required_target_scales = build_required_target_scale_table(
        counts, module, charge_fractions, central, windows
    )
    required_target_scales.to_csv(tables / "required_target_normalization_shifts.csv", index=False)

    template_summary, template_curves = build_empirical_template_fit(
        counts, module, spectrum_edges
    )
    template_summary.to_csv(tables / "empirical_hydrogen_template_fit_summary.csv", index=False)
    template_curves.to_csv(tables / "empirical_hydrogen_template_fit_curves.csv", index=False)

    write_summary_json(metadata / "summary.json", module, windows, flat)
    write_json(metadata / "configuration.json", json_safe({
        "variant": label,
        "production_module": str(module_path.resolve()),
        "cut_json": str(cut_json.resolve()),
        "tree": tree_name,
        "workers": workers,
        "replicas": args.replicas,
        "seed": args.seed,
        "performance": {
            "vectorized_event_binning": True,
            "vectorized_control_window_projection": True,
            "bootstrap_replica_chunk_parallelism": True,
            "maximum_workers": MAXIMUM_WORKERS,
        },
        "control_windows": windows,
        "method1_nominal_normalization_window": list(NOMINAL_CONTROL_WINDOW),
        "low_mx2_closure_regions": DEFAULT_LOW_REGIONS,
        "spectrum_range": [args.spectrum_min, args.spectrum_max],
        "spectrum_bins": args.spectrum_bins,
        "cumulative_upper_edges": cumulative_edges,
        "method2_target_sensitivity_fraction": SENSITIVITY_FRACTION,
        "required_target_scale_scan": [REQUIRED_SCALE_MIN, REQUIRED_SCALE_MAX, REQUIRED_SCALE_GRID],
        "carbon_to_nitrogen_study": {
            "nominal_r_cn": NOMINAL_R_CN,
            "scan_range": [args.r_cn_min, args.r_cn_max],
            "scan_points": args.r_cn_points,
            "definition": "r_cn = sigma_C / sigma_N",
            "implementation": "effective pure-carbon proxy N_C * (nominal_r_cn / r_cn), exact at r_cn=6/7",
            "fit_uncertainties": "required r_cn constant and slope fits are repeated replica by replica with fixed inverse-marginal-variance weights; the replica spread is quoted",
            "legacy_independent_bin_fit_table": "required_r_cn_fit_summary_legacy_independent_bins.csv",
            "thermal_contraction_variations": False,
        },
        "empirical_template_fit": {
            "fit_range": [TEMPLATE_FIT_MIN, TEMPLATE_FIT_MAX],
            "control_region": list(NOMINAL_CONTROL_WINDOW),
            "polynomial_order": TEMPLATE_POLYNOMIAL_ORDER,
            "hydrogen_template": "CH2 - beta*C, with beta from the control region",
        },
        "inputs": {p: {t: str(v) for t, v in d.items()} for p, d in input_paths.items()},
    }))

    plot_window_scan_summary(plots, module, windows, central, bootstrap)
    plot_cumulative_window_scan(plots, module, windows, cumulative_edges, central, bootstrap)
    plot_grouping_comparison(plots, module, windows, central, bootstrap)
    plot_method_difference_summary(plots, module, windows, central, bootstrap)
    plot_low_region_closure(plots, closure, module)
    plot_required_transfer_factor(plots, transfer, module)
    plot_required_transfer_kinematics(plots, transfer, module)
    plot_required_rcn(plots, required_rcn, module)
    plot_rcn_replica_fit_summary(plots, rcn_replica_fit, module)
    plot_rcn_scan(plots, rcn_scan, module)
    plot_spectrum_comparison(
        plots, spectrum_comparison, module, NOMINAL_CONTROL_WINDOW
    )
    plot_cumulative_background_comparison(plots, cumulative_background, module)
    plot_target_sensitivity(plots, target_sensitivity, module)
    plot_required_target_scales(plots, required_target_scales, module)
    plot_empirical_template_fits(plots, template_summary, template_curves, module)

    print(f"[{label}] Validation products written to {root.resolve()}")


def main() -> int:
    args = build_parser().parse_args()
    if args.replicas < 2:
        raise ValueError("--replicas must be at least 2")
    if args.spectrum_bins < 10:
        raise ValueError("--spectrum-bins must be at least 10")
    if not args.spectrum_min < args.spectrum_max:
        raise ValueError("--spectrum-min must be below --spectrum-max")
    if not args.template_fit_min < args.template_fit_max:
        raise ValueError("--template-fit-min must be below --template-fit-max")
    global TEMPLATE_FIT_MIN, TEMPLATE_FIT_MAX, DEFAULT_R_CN_MIN, DEFAULT_R_CN_MAX, NOMINAL_CONTROL_WINDOW
    if not args.normalization_window_min < args.normalization_window_max:
        raise ValueError(
            "--normalization-window-min must be below --normalization-window-max"
        )
    NOMINAL_CONTROL_WINDOW = (
        round(float(args.normalization_window_min), 10),
        round(float(args.normalization_window_max), 10),
    )
    TEMPLATE_FIT_MIN = float(args.template_fit_min)
    TEMPLATE_FIT_MAX = float(args.template_fit_max)
    if not 0.0 < args.r_cn_min < args.r_cn_max:
        raise ValueError("Require 0 < --r-cn-min < --r-cn-max")
    if args.r_cn_points < 3:
        raise ValueError("--r-cn-points must be at least 3")
    DEFAULT_R_CN_MIN = float(args.r_cn_min)
    DEFAULT_R_CN_MAX = float(args.r_cn_max)

    module_path = args.dilution_module.resolve()
    module = import_production_module(module_path)
    cumulative_edges = cumulative_upper_edges(
        args.cumulative_upper_min, args.cumulative_upper_max, args.cumulative_upper_step
    )
    windows = canonical_windows(args.control_window, cumulative_edges)
    charge_fractions = module.load_charge_fractions(args.charge_fractions_json)
    module.validate_fraction_table(charge_fractions)

    standard_inputs = resolve_inputs(module, args.input, isr=False)
    standard_cuts = resolve_cut_json(module, args, isr=False)
    run_one_variant(
        "nominal", module, module_path, standard_inputs, standard_cuts,
        charge_fractions, args, windows, cumulative_edges
    )

    if args.run_isr:
        isr_inputs = resolve_inputs(module, args.isr_input, isr=True)
        isr_cuts = resolve_cut_json(module, args, isr=True)
        run_one_variant(
            "isr", module, module_path, isr_inputs, isr_cuts,
            charge_fractions, args, windows, cumulative_edges
        )

    print("\nCarbon-normalization validation complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
