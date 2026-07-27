#!/usr/bin/env python3
"""
validate_carbon_normalization_v2.py

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

# The nominal window remains [0.00, 0.40).  These scans explicitly vary both
# edges.  Duplicate nominal entries are removed during canonicalization.
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

DEFAULT_SPECTRUM_RANGE = (-0.20, 1.40)
DEFAULT_SPECTRUM_BINS = 80
DEFAULT_CUMULATIVE_UPPER_MIN = 0.15
DEFAULT_CUMULATIVE_UPPER_MAX = 0.65
DEFAULT_CUMULATIVE_UPPER_STEP = 0.025
PERIOD_GROUP = "period"
XB_GROUP = "xB-row"
TPRIME_GROUP = "tprime-column"
BIN_GROUP = "bin-local"
GROUPINGS = (PERIOD_GROUP, XB_GROUP, TPRIME_GROUP, BIN_GROUP)


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
    if (0.0, 0.4) not in unique:
        unique.append((0.0, 0.4))
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
    n_bins = len(xb_bins) * len(t_bins)
    atomic = np.zeros((n_bins, len(atomic_edges) - 1), dtype=np.int64)
    signal = np.zeros((n_bins, 3), dtype=np.int64)
    low = np.zeros((n_bins, len(low_regions)), dtype=np.int64)
    spectrum = np.zeros((n_bins, len(spectrum_edges) - 1), dtype=np.int64)

    index = 0
    for ix, xb_range in enumerate(xb_bins):
        for it, t_range in enumerate(t_bins):
            base = bin_mask(dataset, xb_range, t_range)
            values = dataset["mx2_gev2"][base]
            atomic[index], _ = np.histogram(values, bins=atomic_edges)
            spectrum[index], _ = np.histogram(values, bins=spectrum_edges)
            for region_index, (lo, hi) in enumerate(low_regions):
                low[index, region_index] = int(np.count_nonzero((values >= lo) & (values < hi)))

            cut = cuts[(period, ix, it)]
            for cut_index, name in enumerate(("tight", "nominal", "loose")):
                lo, hi = cut[name]
                signal[index, cut_index] = int(np.count_nonzero((values >= lo) & (values < hi)))
            index += 1

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


def controls_for_windows(
    atomic: np.ndarray,
    atomic_edges: np.ndarray,
    windows: list[tuple[float, float]],
) -> np.ndarray:
    values = []
    for window in windows:
        selected = window_interval_mask(atomic_edges, window)
        values.append(np.sum(atomic[..., selected], axis=-1))
    return np.stack(values, axis=-1)


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
) -> dict[str, Any]:
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
    alpha = {
        grouping: np.full((replicas, n_bins, n_windows), np.nan)
        for grouping in groupings
    }
    batch_size = min(250, replicas)

    interval_masks = [window_interval_mask(atomic_edges, window) for window in windows]

    for start in range(0, replicas, batch_size):
        stop = min(replicas, start + batch_size)
        batch = stop - start
        control_by_target: dict[str, np.ndarray] = {}
        signal_by_target: dict[str, np.ndarray] = {}
        for target, arrays in observed.items():
            atomic_rep = rng.poisson(arrays["atomic"][np.newaxis, ...], size=(batch,) + arrays["atomic"].shape)
            signal_rep = rng.poisson(arrays["signal"][np.newaxis, ...], size=(batch,) + arrays["signal"].shape)
            control_by_target[target] = np.stack(
                [np.sum(atomic_rep[..., mask], axis=-1) for mask in interval_masks], axis=-1
            )
            signal_by_target[target] = signal_rep

        for grouping in groupings:
            scale = np.full((batch, n_bins, n_windows), np.nan)
            gids = group_ids(grouping, n_bins, n_t)
            for gid in range(group_count(grouping, n_x, n_t)):
                members = gids == gid
                nh3_sum = np.sum(control_by_target["NH3"][:, members, :], axis=1)
                c_sum = np.sum(control_by_target["C"][:, members, :], axis=1)
                scale[:, members, :] = ratio(nh3_sum, c_sum)[:, np.newaxis, :]
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
        nominal_window_index = windows.index((0.0, 0.4))
        nominal_period_alpha = alpha[PERIOD_GROUP][start:stop, :, nominal_window_index]
        transfer_ratio[start:stop] = ratio(required, nominal_period_alpha[..., np.newaxis])

        predicted_nominal_background = (
            nominal_period_alpha[..., np.newaxis] * signal_by_target["C"]
        )
        required_background = signal_by_target["NH3"] * (1.0 - m2_batch)
        kappa_period[start:stop] = ratio(
            np.sum(required_background, axis=1),
            np.sum(predicted_nominal_background, axis=1),
        )

    return {
        "period": period,
        "method1": method1,
        "method2": method2,
        "alpha": alpha,
        "alpha_required": alpha_required,
        "transfer_ratio": transfer_ratio,
        "kappa_period": kappa_period,
    }


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
    output: dict[str, Any] = {}
    n_x = len(module.XB_BINS)
    n_t = len(module.MINUS_TPRIME_BINS_GEV2)
    with ProcessPoolExecutor(max_workers=min(workers, len(module.PERIODS))) as executor:
        futures = {}
        for period_index, period in enumerate(module.PERIODS):
            observed = {
                target: {
                    "atomic": counts[(period, target)].atomic_control,
                    "signal": counts[(period, target)].signal,
                }
                for target in module.TARGETS
            }
            future = executor.submit(
                validation_bootstrap_worker,
                period,
                observed,
                windows,
                atomic_edges,
                GROUPINGS,
                n_x,
                n_t,
                charge_fractions[period],
                replicas,
                seed + 100003 * period_index,
                str(module_path),
            )
            futures[future] = period
        for future in as_completed(futures):
            period = futures[future]
            output[period] = future.result()
            print(f"Completed {replicas:,} validation replicas for {period}.")
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
        nominal_window_index = windows.index((0.0, 0.4))
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
    nominal_window = windows.index((0.0, 0.4))
    nominal_cut = list(module.CUT_VARIATIONS).index("nominal")
    display_windows = structured_windows(windows)
    display_indices = [windows.index(window) for window in display_windows]
    x = np.arange(len(display_windows))
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(14, 11), sharex=True)
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
    w = windows.index((0.0, 0.4))
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
    w = windows.index((0.0, 0.4))
    cut = list(module.CUT_VARIATIONS).index("nominal")
    x = np.arange(1, module.NUMBER_OF_BINS + 1)
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(16, 11), sharex=True)
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
                "scaled_carbon_over_five_target": float(ratio(summed_scaled_c[i], summed_five[i])),
                "period_carbon_scale": alpha,
            })
    return pd.DataFrame(rows)


def plot_spectrum_comparison(output: Path, table: pd.DataFrame, module: Any) -> None:
    for period in module.PERIODS:
        selected = table[table.period == period]
        fig, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
        centers = selected.mx2_center.to_numpy()
        axes[0].step(centers, selected.nh3_count, where="mid", label="NH$_3$ data")
        axes[0].step(centers, selected.scaled_carbon_background, where="mid", label="Scaled-carbon background")
        axes[0].step(centers, selected.five_target_background, where="mid", label="Direct five-target background")
        axes[0].set_ylabel("Counts")
        axes[0].set_yscale("symlog", linthresh=1.0)
        axes[0].grid(alpha=0.25)
        axes[0].legend()
        axes[0].set_title(f"{module.PERIOD_LABELS[period]} — period-integrated background-shape comparison")
        axes[1].plot(centers, selected.scaled_carbon_over_five_target, "o")
        axes[1].axhline(1.0, linewidth=1)
        axes[1].set_ylabel("Scaled C / five-target")
        axes[1].set_xlabel(r"$M_X^2$ (GeV$^2$)")
        axes[1].grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(output / f"background_shape_comparison_{period}.png", dpi=180)
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
    fig, axes = plt.subplots(len(module.PERIODS), 3, figsize=(17, 11), sharex=True)
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
    fig, axes = plt.subplots(len(module.PERIODS), 1, figsize=(16, 11), sharex=True)
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
    fig, axes = plt.subplots(len(module.PERIODS), 2, figsize=(15, 11))
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
    fig, axes = plt.subplots(len(module.PERIODS), 2, figsize=(15, 11), sharex=True)
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
    parser.add_argument("--spectrum-min", type=float, default=DEFAULT_SPECTRUM_RANGE[0])
    parser.add_argument("--spectrum-max", type=float, default=DEFAULT_SPECTRUM_RANGE[1])
    parser.add_argument("--spectrum-bins", type=int, default=DEFAULT_SPECTRUM_BINS)
    parser.add_argument("--cumulative-upper-min", type=float, default=DEFAULT_CUMULATIVE_UPPER_MIN)
    parser.add_argument("--cumulative-upper-max", type=float, default=DEFAULT_CUMULATIVE_UPPER_MAX)
    parser.add_argument("--cumulative-upper-step", type=float, default=DEFAULT_CUMULATIVE_UPPER_STEP)
    parser.add_argument("--input", action="append", default=[],
                        help="Override standard input as period:target=/path.root")
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

    cumulative_background = build_cumulative_background_comparison(
        counts, module, charge_fractions, spectrum_edges, cumulative_edges
    )
    cumulative_background.to_csv(tables / "cumulative_background_comparison.csv", index=False)

    write_summary_json(metadata / "summary.json", module, windows, flat)
    write_json(metadata / "configuration.json", json_safe({
        "variant": label,
        "production_module": str(module_path.resolve()),
        "cut_json": str(cut_json.resolve()),
        "tree": tree_name,
        "workers": workers,
        "replicas": args.replicas,
        "seed": args.seed,
        "control_windows": windows,
        "low_mx2_closure_regions": DEFAULT_LOW_REGIONS,
        "spectrum_range": [args.spectrum_min, args.spectrum_max],
        "spectrum_bins": args.spectrum_bins,
        "cumulative_upper_edges": cumulative_edges,
        "inputs": {p: {t: str(v) for t, v in d.items()} for p, d in input_paths.items()},
    }))

    plot_window_scan_summary(plots, module, windows, central, bootstrap)
    plot_cumulative_window_scan(plots, module, windows, cumulative_edges, central, bootstrap)
    plot_grouping_comparison(plots, module, windows, central, bootstrap)
    plot_method_difference_summary(plots, module, windows, central, bootstrap)
    plot_low_region_closure(plots, closure, module)
    plot_required_transfer_factor(plots, transfer, module)
    plot_required_transfer_kinematics(plots, transfer, module)
    plot_cumulative_background_comparison(plots, cumulative_background, module)

    print(f"[{label}] Validation products written to {root.resolve()}")


def main() -> int:
    args = build_parser().parse_args()
    if args.replicas < 2:
        raise ValueError("--replicas must be at least 2")
    if args.spectrum_bins < 10:
        raise ValueError("--spectrum-bins must be at least 10")
    if not args.spectrum_min < args.spectrum_max:
        raise ValueError("--spectrum-min must be below --spectrum-max")

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
