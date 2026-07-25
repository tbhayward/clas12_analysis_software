#!/usr/bin/env python3
"""
determine_dilution_factors_v1.py

Determine the RGC exclusive e pi+ dilution factor in the fixed 4 xB by 6
(-tprime) bins using three complementary methods:

  Method 1
    Carbon subtraction with exactly the period-wide raw-count normalization
    used by channel_selection_mx2_fits_v23.py:

        alpha_C = N_NH3(0.20 <= Mx2 < 0.40) /
                  N_C  (0.20 <= Mx2 < 0.40)

        f = (N_NH3 - alpha_C N_C) / N_NH3.

  Method 2
    Direct auxiliary-target dilution factor from Eq. (10) of the standalone
    RGC dilution-factor note, using charge-normalized rates from NH3, C, CH2,
    He-bath (MT), and empty-target/foils-only (F) data.

  Method 3
    Packing fraction from Eq. (11), followed by the nonlinear dilution-factor
    expression in Eq. (14).  A per-bin packing fraction is retained as a
    diagnostic.  The production-style Method-3 calculation uses a period-wide
    packing fraction extracted from all 24 kinematic bins for the same Mx2 cut
    variation, because packing fraction is a target-cell property rather than
    a kinematic observable.

The program is intended to be run from

    RGC_enpi+/dilution_factor/

while the channel-selection outputs are in the sibling directory

    RGC_enpi+/channel_selection/

The default exclusivity-cut JSON is therefore

    ../channel_selection/output/channel_selection_mx2_fit_stability/
        final_carbon_assisted_cuts/tables/
        final_carbon_assisted_mx2_cuts_v23.json

All ROOT inputs default to the finalized momentum-corrected paper_versions
files.  The program uses at most seven worker processes.  It writes complete
JSON and CSV products, compact downstream JSON, statistical covariance and
correlation matrices, and diagnostic plots.

Statistical treatment
---------------------
Counts are modeled as a Poisson point process.  For each period, target, and
kinematic bin, events are classified by their complete membership pattern in
four possibly overlapping selections:

    control, tight, nominal, loose.

The 16 disjoint membership-pattern counts are independently Poisson-resampled.
Summing those replicas reconstructs the selected counts while preserving the
correct statistical covariance for nested or overlapping windows.  The shared
Method-1 carbon scale, the shared period packing fraction, and all three cut
variations are recalculated in every replica.  This yields statistically
consistent uncertainties and full 24-bin covariance matrices.

No dilution factor or packing fraction is clipped to a physical interval.
Out-of-range values are retained and flagged as diagnostics.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


# =============================================================================
# Fixed analysis definitions
# =============================================================================

PERIODS: tuple[str, ...] = ("su22", "fa22", "sp23")
TARGETS: tuple[str, ...] = ("NH3", "C", "CH2", "He", "ET")
CUT_VARIATIONS: tuple[str, ...] = ("tight", "nominal", "loose")
SELECTION_NAMES: tuple[str, ...] = ("control", "tight", "nominal", "loose")

PERIOD_LABELS: dict[str, str] = {
    "su22": "Su22",
    "fa22": "Fa22",
    "sp23": "Sp23",
}

TARGET_LABELS: dict[str, str] = {
    "NH3": "NH$_3$",
    "C": "C",
    "CH2": "CH$_2$",
    "He": "He bath",
    "ET": "Empty target",
}

XB_BINS: tuple[tuple[float, float], ...] = (
    (0.10, 0.25),
    (0.25, 0.35),
    (0.35, 0.45),
    (0.45, 0.60),
)

MINUS_TPRIME_BINS_GEV2: tuple[tuple[float, float], ...] = (
    (0.05, 0.25),
    (0.25, 0.45),
    (0.45, 0.65),
    (0.65, 0.85),
    (0.85, 1.05),
    (1.05, 1.25),
)

NUMBER_OF_BINS = len(XB_BINS) * len(MINUS_TPRIME_BINS_GEV2)
MAXIMUM_WORKERS = 7

DEFAULT_TREE_NAME = "PhysicsEvents"
DEFAULT_CONTROL_MIN_GEV2 = 0.20
DEFAULT_CONTROL_MAX_GEV2 = 0.40
DEFAULT_REPLICAS = 10000
DEFAULT_SEED = 7302026

DEFAULT_OUTPUT_DIR = Path("output/dilution_factor_determination")
DEFAULT_CUT_JSON = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "final_carbon_assisted_cuts/tables/"
    "final_carbon_assisted_mx2_cuts_v23.json"
)

PAPER_VERSIONS_DIR = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
    "paper_versions"
)

DEFAULT_INPUTS: dict[str, dict[str, Path]] = {
    period: {
        target: PAPER_VERSIONS_DIR
        / f"rgc_{period}_inb_{target}_epi+_mom_corrections.root"
        for target in TARGETS
    }
    for period in PERIODS
}

# Relative accumulated-charge fractions supplied for this analysis.  They are
# sufficient because the auxiliary-target equations are homogeneous in the
# charge-normalized rates; the common total charge cancels.
DEFAULT_CHARGE_FRACTIONS: dict[str, dict[str, float]] = {
    "su22": {
        "NH3": 0.7624,
        "C": 0.0747,
        "CH2": 0.0390,
        "He": 0.0249,
        "ET": 0.0989,
    },
    "fa22": {
        "NH3": 0.5839,
        "C": 0.1992,
        "CH2": 0.1802,
        "He": 0.0301,
        "ET": 0.0066,
    },
    "sp23": {
        "NH3": 0.4543,
        "C": 0.1509,
        "CH2": 0.1966,
        "He": 0.1205,
        "ET": 0.0777,
    },
}

BRANCH_ALIASES: dict[str, tuple[str, ...]] = {
    "xB": ("xB", "x", "xb", "x_b"),
    "tprime": ("tprime", "t_prime", "tp", "tPrime"),
    "Mx2": (
        "Mx2",
        "mx2",
        "Mx2_epi",
        "Mx2_epip",
        "missing_mass_squared",
        "missing_mass2",
    ),
}

# Numerical coefficients from Eqs. (10), (11), and (14) for 15NH3.
# In the note's notation:
#   A  -> NH3
#   C  -> carbon
#   CH -> CH2
#   MT -> cell/foils immersed in liquid helium
#   F  -> foils-only / empty target
EQ10_EQ11_EQ14_COEFFICIENTS: dict[str, float] = {
    "ch2_carbon_hydrogen_combination": 0.880502,
    "he_hydrogen_combination": 0.16604,
    "empty_hydrogen_combination": -0.285536,
    "ch2_carbon_nuclear_combination": 0.195667,
    "he_nuclear_combination": -0.82744,
    "empty_nuclear_combination": 0.023106,
    "packing_prefactor": 0.50734,
}


# =============================================================================
# Data containers
# =============================================================================

@dataclass(frozen=True)
class DatasetSpec:
    period: str
    target: str
    file_path: str
    tree_name: str


@dataclass(frozen=True)
class LoadedDataset:
    period: str
    target: str
    file_path: str
    tree_name: str
    x_branch: str
    tprime_branch: str
    mx2_branch: str
    xB: np.ndarray
    minus_tprime_gev2: np.ndarray
    mx2_gev2: np.ndarray


@dataclass(frozen=True)
class CutEntry:
    period: str
    x_index: int
    t_index: int
    xB_min: float
    xB_max: float
    minus_tprime_min_gev2: float
    minus_tprime_max_gev2: float
    mu_gev2: float
    mu_error_gev2: float
    sigma_gev2: float
    sigma_error_gev2: float
    classification: str
    tight_min_gev2: float
    tight_max_gev2: float
    nominal_min_gev2: float
    nominal_max_gev2: float
    loose_min_gev2: float
    loose_max_gev2: float

    def interval(self, variation: str) -> tuple[float, float]:
        if variation == "tight":
            return self.tight_min_gev2, self.tight_max_gev2
        # endif
        if variation == "nominal":
            return self.nominal_min_gev2, self.nominal_max_gev2
        # endif
        if variation == "loose":
            return self.loose_min_gev2, self.loose_max_gev2
        # endif
        raise KeyError(f"Unknown cut variation: {variation}")


@dataclass(frozen=True)
class CountPayload:
    """Observed counts and disjoint membership-pattern counts."""

    period: str
    target: str
    x_index: int
    t_index: int
    bin_number: int
    total_kinematic_count: int
    selected_counts: tuple[int, int, int, int]
    pattern_counts: tuple[int, ...]


# =============================================================================
# General utilities
# =============================================================================

def ensure_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    # endif
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    # endif
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    # endif
    if isinstance(value, (np.integer,)):
        return int(value)
    # endif
    if isinstance(value, (np.floating,)):
        value = float(value)
    # endif
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    # endif
    if isinstance(value, Path):
        return str(value)
    # endif
    return value


def write_json(path: Path, payload: Any) -> None:
    ensure_directory(path.parent)
    path.write_text(
        json.dumps(json_safe(payload), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
        # endfor
    # endwith
    return digest.hexdigest()


def combined_bin_number(x_index: int, t_index: int) -> int:
    return x_index * len(MINUS_TPRIME_BINS_GEV2) + t_index + 1


def parse_input_override(text: str) -> tuple[str, str, Path]:
    try:
        left, path_text = text.split("=", 1)
        period, target = left.split(":", 1)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "Input overrides must have the form PERIOD:TARGET=/path/file.root"
        ) from exc
    # endtry
    period = period.strip().lower()
    target_lookup = {item.lower(): item for item in TARGETS}
    if period not in PERIODS:
        raise argparse.ArgumentTypeError(
            f"Unknown period {period!r}; expected one of {PERIODS}."
        )
    # endif
    if target.strip().lower() not in target_lookup:
        raise argparse.ArgumentTypeError(
            f"Unknown target {target!r}; expected one of {TARGETS}."
        )
    # endif
    return period, target_lookup[target.strip().lower()], Path(path_text).expanduser()


def validate_fraction_table(fractions: dict[str, dict[str, float]]) -> None:
    for period in PERIODS:
        if period not in fractions:
            raise RuntimeError(f"Missing charge fractions for {period}.")
        # endif
        for target in TARGETS:
            if target not in fractions[period]:
                raise RuntimeError(
                    f"Missing charge fraction for {period}/{target}."
                )
            # endif
            value = float(fractions[period][target])
            if not math.isfinite(value) or value <= 0.0:
                raise RuntimeError(
                    f"Invalid charge fraction for {period}/{target}: {value}"
                )
            # endif
        # endfor
        total = sum(float(fractions[period][target]) for target in TARGETS)
        if abs(total - 1.0) > 5.0e-4:
            raise RuntimeError(
                f"Charge fractions for {period} sum to {total:.8f}, not 1."
            )
        # endif
    # endfor


def load_charge_fractions(path: Path | None) -> dict[str, dict[str, float]]:
    if path is None:
        fractions = {
            period: dict(values)
            for period, values in DEFAULT_CHARGE_FRACTIONS.items()
        }
    else:
        payload = json.loads(path.read_text(encoding="utf-8"))
        source = payload.get("charge_fractions", payload)
        fractions = {
            str(period).lower(): {
                str(target): float(value)
                for target, value in values.items()
            }
            for period, values in source.items()
        }
    # endif
    validate_fraction_table(fractions)
    return fractions


def resolve_branch(tree: uproot.behaviors.TTree.TTree, aliases: Iterable[str]) -> str:
    available = set(str(key).split(";")[0] for key in tree.keys())
    for alias in aliases:
        if alias in available:
            return alias
        # endif
    # endfor
    raise RuntimeError(
        f"None of the branch aliases {tuple(aliases)} exists. "
        f"Available branches include: {sorted(available)[:80]}"
    )


# =============================================================================
# Input loading and cut parsing
# =============================================================================

def load_one_dataset(spec: DatasetSpec) -> LoadedDataset:
    path = Path(spec.file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing ROOT input: {path}")
    # endif

    with uproot.open(path) as root_file:
        if spec.tree_name not in root_file:
            available = [str(key).split(";")[0] for key in root_file.keys()]
            raise RuntimeError(
                f"Tree {spec.tree_name!r} not found in {path}. "
                f"Available keys: {available}"
            )
        # endif
        tree = root_file[spec.tree_name]
        x_branch = resolve_branch(tree, BRANCH_ALIASES["xB"])
        tprime_branch = resolve_branch(tree, BRANCH_ALIASES["tprime"])
        mx2_branch = resolve_branch(tree, BRANCH_ALIASES["Mx2"])
        arrays = tree.arrays(
            [x_branch, tprime_branch, mx2_branch],
            library="np",
        )
    # endwith

    xB = np.asarray(arrays[x_branch], dtype=np.float64)
    minus_tprime = -np.asarray(arrays[tprime_branch], dtype=np.float64)
    mx2 = np.asarray(arrays[mx2_branch], dtype=np.float64)
    finite = np.isfinite(xB) & np.isfinite(minus_tprime) & np.isfinite(mx2)

    return LoadedDataset(
        period=spec.period,
        target=spec.target,
        file_path=str(path.resolve()),
        tree_name=spec.tree_name,
        x_branch=x_branch,
        tprime_branch=tprime_branch,
        mx2_branch=mx2_branch,
        xB=xB[finite],
        minus_tprime_gev2=minus_tprime[finite],
        mx2_gev2=mx2[finite],
    )


def load_all_datasets(
    input_paths: dict[str, dict[str, Path]],
    tree_name: str,
    workers: int,
) -> dict[tuple[str, str], LoadedDataset]:
    specs = [
        DatasetSpec(period, target, str(input_paths[period][target]), tree_name)
        for period in PERIODS
        for target in TARGETS
    ]
    loaded: dict[tuple[str, str], LoadedDataset] = {}

    with ProcessPoolExecutor(max_workers=workers) as executor:
        future_to_spec = {
            executor.submit(load_one_dataset, spec): spec
            for spec in specs
        }
        for completed, future in enumerate(as_completed(future_to_spec), start=1):
            spec = future_to_spec[future]
            dataset = future.result()
            loaded[(dataset.period, dataset.target)] = dataset
            print(
                f"Loaded {completed:2d}/{len(specs)}: "
                f"{PERIOD_LABELS[dataset.period]} {dataset.target} "
                f"({dataset.mx2_gev2.size:,} finite events)"
            )
        # endfor
    # endwith

    return loaded


def load_exclusivity_cuts(path: Path) -> dict[tuple[str, int, int], CutEntry]:
    if not path.is_file():
        raise FileNotFoundError(
            "Exclusivity-cut JSON was not found.  The default assumes this "
            "script is run from RGC_enpi+/dilution_factor and channel selection "
            f"was run in RGC_enpi+/channel_selection.  Missing path: {path}"
        )
    # endif

    payload = json.loads(path.read_text(encoding="utf-8"))
    periods_payload = payload.get("periods")
    if not isinstance(periods_payload, dict):
        raise RuntimeError(
            f"Cut JSON {path} does not contain the expected 'periods' object."
        )
    # endif

    cuts: dict[tuple[str, int, int], CutEntry] = {}
    for period in PERIODS:
        entries = periods_payload.get(period)
        if not isinstance(entries, list):
            raise RuntimeError(f"Cut JSON has no list for period {period}.")
        # endif
        for item in entries:
            tight = item["tight"]
            nominal = item["nominal"]
            loose = item["loose"]
            entry = CutEntry(
                period=period,
                x_index=int(item["x_index"]),
                t_index=int(item["t_index"]),
                xB_min=float(item["xB_min"]),
                xB_max=float(item["xB_max"]),
                minus_tprime_min_gev2=float(item["minus_tprime_min_gev2"]),
                minus_tprime_max_gev2=float(item["minus_tprime_max_gev2"]),
                mu_gev2=float(item["mu_gev2"]),
                mu_error_gev2=float(item["mu_error_gev2"]),
                sigma_gev2=float(item["sigma_gev2"]),
                sigma_error_gev2=float(item["sigma_error_gev2"]),
                classification=str(item.get("classification", "unknown")),
                tight_min_gev2=float(tight[0]),
                tight_max_gev2=float(tight[1]),
                nominal_min_gev2=float(nominal[0]),
                nominal_max_gev2=float(nominal[1]),
                loose_min_gev2=float(loose[0]),
                loose_max_gev2=float(loose[1]),
            )
            cuts[(period, entry.x_index, entry.t_index)] = entry
        # endfor
    # endfor

    expected = {
        (period, x_index, t_index)
        for period in PERIODS
        for x_index in range(len(XB_BINS))
        for t_index in range(len(MINUS_TPRIME_BINS_GEV2))
    }
    missing = sorted(expected - set(cuts))
    extra = sorted(set(cuts) - expected)
    if missing or extra:
        raise RuntimeError(
            f"Cut-table coverage mismatch. Missing={missing}; extra={extra}."
        )
    # endif

    for key, entry in cuts.items():
        expected_x = XB_BINS[entry.x_index]
        expected_t = MINUS_TPRIME_BINS_GEV2[entry.t_index]
        values_match = (
            np.allclose([entry.xB_min, entry.xB_max], expected_x)
            and np.allclose(
                [entry.minus_tprime_min_gev2, entry.minus_tprime_max_gev2],
                expected_t,
            )
        )
        if not values_match:
            raise RuntimeError(
                f"Cut-table bin edges for {key} do not match the fixed analysis binning."
            )
        # endif
        for variation in CUT_VARIATIONS:
            low, high = entry.interval(variation)
            if not (math.isfinite(low) and math.isfinite(high) and low < high):
                raise RuntimeError(f"Invalid {variation} interval for {key}: {(low, high)}")
            # endif
        # endfor
    # endfor

    return cuts


# =============================================================================
# Counting and Poisson-membership representation
# =============================================================================

def membership_pattern_counts(masks: list[np.ndarray]) -> tuple[int, ...]:
    if len(masks) != len(SELECTION_NAMES):
        raise ValueError("Expected exactly four selection masks.")
    # endif
    if not masks:
        return tuple()
    # endif
    pattern = np.zeros(masks[0].size, dtype=np.uint8)
    for bit, mask in enumerate(masks):
        pattern |= np.asarray(mask, dtype=np.uint8) << bit
    # endfor
    counts = np.bincount(pattern, minlength=1 << len(masks))
    return tuple(int(value) for value in counts)


def selection_counts_from_patterns(pattern_counts: np.ndarray) -> np.ndarray:
    """Convert [...,16] membership-pattern counts to [...,4] selections."""
    bits = np.arange(16, dtype=np.uint8)
    membership = np.stack(
        [((bits >> bit) & 1).astype(np.int64) for bit in range(4)],
        axis=1,
    )
    return np.asarray(pattern_counts) @ membership


def count_one_dataset(
    dataset: LoadedDataset,
    cuts: dict[tuple[str, int, int], CutEntry],
    control_min_gev2: float,
    control_max_gev2: float,
) -> list[CountPayload]:
    rows: list[CountPayload] = []
    for x_index, (x_min, x_max) in enumerate(XB_BINS):
        x_mask = (dataset.xB >= x_min) & (dataset.xB < x_max)
        for t_index, (t_min, t_max) in enumerate(MINUS_TPRIME_BINS_GEV2):
            kinematic_mask = (
                x_mask
                & (dataset.minus_tprime_gev2 >= t_min)
                & (dataset.minus_tprime_gev2 < t_max)
            )
            values = dataset.mx2_gev2[kinematic_mask]
            entry = cuts[(dataset.period, x_index, t_index)]
            masks = [
                (values >= control_min_gev2) & (values < control_max_gev2),
            ]
            for variation in CUT_VARIATIONS:
                low, high = entry.interval(variation)
                masks.append((values >= low) & (values < high))
            # endfor
            patterns = membership_pattern_counts(masks)
            selected = selection_counts_from_patterns(
                np.asarray(patterns, dtype=np.int64)
            )
            rows.append(
                CountPayload(
                    period=dataset.period,
                    target=dataset.target,
                    x_index=x_index,
                    t_index=t_index,
                    bin_number=combined_bin_number(x_index, t_index),
                    total_kinematic_count=int(values.size),
                    selected_counts=tuple(int(value) for value in selected),
                    pattern_counts=patterns,
                )
            )
        # endfor
    # endfor
    return rows


def build_count_arrays(
    loaded: dict[tuple[str, str], LoadedDataset],
    cuts: dict[tuple[str, int, int], CutEntry],
    control_min_gev2: float,
    control_max_gev2: float,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, Any]]]:
    """
    Return arrays with shape:

      pattern_counts[period,target,bin,16]
      observed_counts[period,target,bin,4]
    """
    patterns = np.zeros(
        (len(PERIODS), len(TARGETS), NUMBER_OF_BINS, 16),
        dtype=np.int64,
    )
    observed = np.zeros(
        (len(PERIODS), len(TARGETS), NUMBER_OF_BINS, 4),
        dtype=np.int64,
    )
    rows: list[dict[str, Any]] = []

    for p_index, period in enumerate(PERIODS):
        for target_index, target in enumerate(TARGETS):
            payloads = count_one_dataset(
                loaded[(period, target)],
                cuts,
                control_min_gev2,
                control_max_gev2,
            )
            for payload in payloads:
                b = payload.bin_number - 1
                patterns[p_index, target_index, b, :] = payload.pattern_counts
                observed[p_index, target_index, b, :] = payload.selected_counts
                row = {
                    "period": period,
                    "target": target,
                    "x_index": payload.x_index,
                    "t_index": payload.t_index,
                    "bin_number": payload.bin_number,
                    "xB_min": XB_BINS[payload.x_index][0],
                    "xB_max": XB_BINS[payload.x_index][1],
                    "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[payload.t_index][0],
                    "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[payload.t_index][1],
                    "total_kinematic_count": payload.total_kinematic_count,
                }
                for selection_index, selection in enumerate(SELECTION_NAMES):
                    row[f"count_{selection}"] = payload.selected_counts[selection_index]
                # endfor
                rows.append(row)
            # endfor
        # endfor
    # endfor

    return patterns, observed, rows


# =============================================================================
# Physics estimators
# =============================================================================

def safe_divide(numerator: np.ndarray | float, denominator: np.ndarray | float) -> np.ndarray:
    numerator_array = np.asarray(numerator, dtype=np.float64)
    denominator_array = np.asarray(denominator, dtype=np.float64)
    result_shape = np.broadcast_shapes(numerator_array.shape, denominator_array.shape)
    result = np.full(result_shape, np.nan, dtype=np.float64)
    np.divide(
        numerator_array,
        denominator_array,
        out=result,
        where=np.isfinite(denominator_array) & (denominator_array != 0.0),
    )
    return result


def charge_normalized_rates(
    selected_counts: np.ndarray,
    period: str,
    charge_fractions: dict[str, dict[str, float]],
) -> dict[str, np.ndarray]:
    return {
        target: np.asarray(selected_counts[target], dtype=np.float64)
        / float(charge_fractions[period][target])
        for target in TARGETS
    }


def method1_from_counts(
    counts_by_target: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """
    Exact channel-selection carbon normalization.

    Input arrays have shape [...,bin,4], with selection index 0=control and
    1/2/3=tight/nominal/loose.
    """
    nh3 = np.asarray(counts_by_target["NH3"], dtype=np.float64)
    carbon = np.asarray(counts_by_target["C"], dtype=np.float64)
    alpha = safe_divide(np.sum(nh3[..., 0], axis=-1), np.sum(carbon[..., 0], axis=-1))
    # alpha has the leading replica shape, if any. Broadcast over bins/cuts.
    expanded_alpha = alpha[..., np.newaxis, np.newaxis]
    dilution = 1.0 - safe_divide(
        expanded_alpha * carbon[..., 1:4],
        nh3[..., 1:4],
    )
    return dilution, alpha


def auxiliary_combinations(rates: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    c = EQ10_EQ11_EQ14_COEFFICIENTS
    hydrogen_like = (
        c["ch2_carbon_hydrogen_combination"] * rates["CH2"] * rates["C"]
        + c["he_hydrogen_combination"] * rates["He"]
        + c["empty_hydrogen_combination"] * rates["ET"]
    )
    nuclear_like = (
        c["ch2_carbon_nuclear_combination"] * rates["CH2"] * rates["C"]
        + c["he_nuclear_combination"] * rates["He"]
        + c["empty_nuclear_combination"] * rates["ET"]
    )
    return hydrogen_like, nuclear_like


def method2_equation10(rates: dict[str, np.ndarray]) -> np.ndarray:
    hydrogen_like, nuclear_like = auxiliary_combinations(rates)
    numerator = (rates["NH3"] - rates["He"]) * hydrogen_like
    denominator = rates["NH3"] * nuclear_like
    return safe_divide(numerator, denominator)


def packing_fraction_equation11(rates: dict[str, np.ndarray]) -> np.ndarray:
    _, nuclear_like = auxiliary_combinations(rates)
    prefactor = EQ10_EQ11_EQ14_COEFFICIENTS["packing_prefactor"]
    return prefactor * safe_divide(rates["NH3"] * rates["He"], nuclear_like)


def integrated_packing_fraction_equation11(rates: dict[str, np.ndarray]) -> np.ndarray:
    """Integrate rates over all 24 bins before applying Eq. (11)."""
    integrated = {
        target: np.sum(value, axis=-2)
        for target, value in rates.items()
    }
    return packing_fraction_equation11(integrated)


def method3_equation14(
    rates: dict[str, np.ndarray],
    period_packing_fraction: np.ndarray,
) -> np.ndarray:
    hydrogen_like, nuclear_like = auxiliary_combinations(rates)
    prefactor = EQ10_EQ11_EQ14_COEFFICIENTS["packing_prefactor"]
    # period_packing_fraction has shape [...,3 cuts]. Insert the bin axis.
    pf = np.asarray(period_packing_fraction, dtype=np.float64)[..., np.newaxis, :]
    denominator = nuclear_like + safe_divide(prefactor * rates["He"], pf)
    return safe_divide(hydrogen_like, denominator)


def observed_estimators_for_period(
    observed_period: np.ndarray,
    period: str,
    charge_fractions: dict[str, dict[str, float]],
) -> dict[str, Any]:
    counts = {
        target: observed_period[target_index]
        for target_index, target in enumerate(TARGETS)
    }
    method1, alpha = method1_from_counts(counts)
    selected_counts = {
        target: values[:, 1:4]
        for target, values in counts.items()
    }
    rates = charge_normalized_rates(selected_counts, period, charge_fractions)
    method2 = method2_equation10(rates)
    pf_bin = packing_fraction_equation11(rates)
    pf_period = integrated_packing_fraction_equation11(rates)
    method3 = method3_equation14(rates, pf_period)
    return {
        "method1": method1,
        "method1_alpha": alpha,
        "method2": method2,
        "packing_fraction_bin": pf_bin,
        "packing_fraction_period": pf_period,
        "method3": method3,
        "rates": rates,
    }


# =============================================================================
# Statistical replicas
# =============================================================================

def bootstrap_period_worker(
    period_index: int,
    pattern_counts_period: np.ndarray,
    charge_fractions_period: dict[str, float],
    replicas: int,
    seed: int,
) -> dict[str, Any]:
    """Generate correlated Poisson replicas for one run period in bounded-memory batches."""
    rng = np.random.default_rng(seed)
    batch_size = min(250, replicas)

    method1_output = np.full((replicas, NUMBER_OF_BINS, 3), np.nan, dtype=np.float64)
    method2_output = np.full_like(method1_output, np.nan)
    method3_output = np.full_like(method1_output, np.nan)
    pf_bin_output = np.full_like(method1_output, np.nan)
    alpha_output = np.full(replicas, np.nan, dtype=np.float64)
    pf_period_output = np.full((replicas, 3), np.nan, dtype=np.float64)

    for start in range(0, replicas, batch_size):
        stop = min(start + batch_size, replicas)
        number_in_batch = stop - start
        # Shape: batch,target,bin,16.  Batching avoids a multi-gigabyte
        # allocation when the default 10,000 replicas are requested.
        pattern_replicas = rng.poisson(
            lam=pattern_counts_period[np.newaxis, ...],
            size=(number_in_batch,) + pattern_counts_period.shape,
        )
        selected = selection_counts_from_patterns(pattern_replicas)
        counts = {
            target: selected[:, target_index, :, :]
            for target_index, target in enumerate(TARGETS)
        }

        method1, alpha = method1_from_counts(counts)
        selected_for_auxiliary = {
            target: value[..., 1:4]
            for target, value in counts.items()
        }
        rates = {
            target: selected_for_auxiliary[target] / charge_fractions_period[target]
            for target in TARGETS
        }
        method2 = method2_equation10(rates)
        pf_bin = packing_fraction_equation11(rates)
        pf_period = integrated_packing_fraction_equation11(rates)
        method3 = method3_equation14(rates, pf_period)

        method1_output[start:stop] = method1
        method2_output[start:stop] = method2
        method3_output[start:stop] = method3
        pf_bin_output[start:stop] = pf_bin
        alpha_output[start:stop] = alpha
        pf_period_output[start:stop] = pf_period
    # endfor

    return {
        "period_index": period_index,
        "method1": method1_output,
        "method1_alpha": alpha_output,
        "method2": method2_output,
        "packing_fraction_bin": pf_bin_output,
        "packing_fraction_period": pf_period_output,
        "method3": method3_output,
    }


def run_bootstrap(
    pattern_counts: np.ndarray,
    charge_fractions: dict[str, dict[str, float]],
    replicas: int,
    seed: int,
    workers: int,
) -> dict[str, dict[str, np.ndarray]]:
    output: dict[str, dict[str, np.ndarray]] = {}
    with ProcessPoolExecutor(max_workers=min(workers, len(PERIODS))) as executor:
        futures = {
            executor.submit(
                bootstrap_period_worker,
                period_index,
                pattern_counts[period_index],
                charge_fractions[period],
                replicas,
                seed + 100003 * period_index,
            ): period
            for period_index, period in enumerate(PERIODS)
        }
        for future in as_completed(futures):
            period = futures[future]
            result = future.result()
            output[period] = {
                key: value
                for key, value in result.items()
                if key != "period_index"
            }
            print(f"Completed {replicas:,} statistical replicas for {PERIOD_LABELS[period]}.")
        # endfor
    # endwith
    return output


def summarize_replicas(
    central: np.ndarray,
    replicas: np.ndarray,
) -> dict[str, np.ndarray]:
    finite = np.isfinite(replicas)
    valid_fraction = np.mean(finite, axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        standard_deviation = np.nanstd(replicas, axis=0, ddof=1)
        p16 = np.nanpercentile(replicas, 16.0, axis=0)
        p50 = np.nanpercentile(replicas, 50.0, axis=0)
        p84 = np.nanpercentile(replicas, 84.0, axis=0)
    # endwith
    return {
        "central": central,
        "stat_uncertainty": standard_deviation,
        "p16": p16,
        "p50": p50,
        "p84": p84,
        "valid_replica_fraction": valid_fraction,
    }


def pairwise_covariance(samples: np.ndarray) -> np.ndarray:
    """Pairwise-complete covariance for samples shaped [replica,bin]."""
    n_bins = samples.shape[1]
    covariance = np.full((n_bins, n_bins), np.nan, dtype=np.float64)
    for i in range(n_bins):
        for j in range(i, n_bins):
            valid = np.isfinite(samples[:, i]) & np.isfinite(samples[:, j])
            if np.count_nonzero(valid) >= 2:
                value = np.cov(samples[valid, i], samples[valid, j], ddof=1)[0, 1]
                covariance[i, j] = value
                covariance[j, i] = value
            # endif
        # endfor
    # endfor
    return covariance


def covariance_to_correlation(covariance: np.ndarray) -> np.ndarray:
    sigma = np.sqrt(np.diag(covariance))
    denominator = np.outer(sigma, sigma)
    return safe_divide(covariance, denominator)


# =============================================================================
# Output assembly
# =============================================================================

def quality_flags(value: float, valid_fraction: float) -> list[str]:
    flags: list[str] = []
    if not math.isfinite(value):
        flags.append("nonfinite_central_value")
    else:
        if value < 0.0:
            flags.append("below_zero")
        # endif
        if value > 1.0:
            flags.append("above_one")
        # endif
    # endif
    if not math.isfinite(valid_fraction) or valid_fraction < 0.99:
        flags.append("less_than_99_percent_valid_replicas")
    # endif
    return flags


def scalar_record(summary: dict[str, np.ndarray], index: tuple[int, ...]) -> dict[str, Any]:
    value = float(summary["central"][index])
    valid_fraction = float(summary["valid_replica_fraction"][index])
    return {
        "value": value,
        "stat_uncertainty": float(summary["stat_uncertainty"][index]),
        "bootstrap_p16": float(summary["p16"][index]),
        "bootstrap_p50": float(summary["p50"][index]),
        "bootstrap_p84": float(summary["p84"][index]),
        "valid_replica_fraction": valid_fraction,
        "flags": quality_flags(value, valid_fraction),
    }


def build_output_tables(
    observed: np.ndarray,
    central_results: dict[str, dict[str, Any]],
    bootstrap_results: dict[str, dict[str, np.ndarray]],
    charge_fractions: dict[str, dict[str, float]],
    cuts: dict[tuple[str, int, int], CutEntry],
) -> tuple[pd.DataFrame, dict[str, Any], dict[str, Any], dict[str, Any]]:
    flat_rows: list[dict[str, Any]] = []
    master_periods: dict[str, Any] = {}
    compact_periods: dict[str, Any] = {}
    summaries_by_period: dict[str, Any] = {}

    for p_index, period in enumerate(PERIODS):
        central = central_results[period]
        replicas = bootstrap_results[period]
        summaries = {
            "method1": summarize_replicas(central["method1"], replicas["method1"]),
            "method2": summarize_replicas(central["method2"], replicas["method2"]),
            "method3": summarize_replicas(central["method3"], replicas["method3"]),
            "packing_fraction_bin": summarize_replicas(
                central["packing_fraction_bin"], replicas["packing_fraction_bin"]
            ),
            "packing_fraction_period": summarize_replicas(
                central["packing_fraction_period"], replicas["packing_fraction_period"]
            ),
        }
        alpha_summary = summarize_replicas(
            np.asarray(central["method1_alpha"]),
            np.asarray(replicas["method1_alpha"]),
        )
        summaries["method1_alpha"] = alpha_summary
        summaries_by_period[period] = summaries

        period_bins: list[dict[str, Any]] = []
        compact_bins: list[dict[str, Any]] = []
        for b in range(NUMBER_OF_BINS):
            x_index = b // len(MINUS_TPRIME_BINS_GEV2)
            t_index = b % len(MINUS_TPRIME_BINS_GEV2)
            cut_entry = cuts[(period, x_index, t_index)]
            cut_payload: dict[str, Any] = {}
            compact_cut_payload: dict[str, Any] = {}
            for cut_index, variation in enumerate(CUT_VARIATIONS):
                selection_index = cut_index + 1
                counts_payload = {
                    target: int(observed[p_index, target_index, b, selection_index])
                    for target_index, target in enumerate(TARGETS)
                }
                rates_payload = {
                    target: float(central["rates"][target][b, cut_index])
                    for target in TARGETS
                }
                method1_record = scalar_record(summaries["method1"], (b, cut_index))
                method2_record = scalar_record(summaries["method2"], (b, cut_index))
                method3_record = scalar_record(summaries["method3"], (b, cut_index))
                pf_bin_record = scalar_record(
                    summaries["packing_fraction_bin"], (b, cut_index)
                )
                pf_period_record = scalar_record(
                    summaries["packing_fraction_period"], (cut_index,)
                )
                low, high = cut_entry.interval(variation)
                cut_payload[variation] = {
                    "mx2_min_gev2": low,
                    "mx2_max_gev2": high,
                    "counts": counts_payload,
                    "charge_normalized_rates_relative_units": rates_payload,
                    "method1_carbon_subtraction": method1_record,
                    "method2_equation10": method2_record,
                    "packing_fraction_equation11_bin_diagnostic": pf_bin_record,
                    "packing_fraction_equation11_period_integrated": pf_period_record,
                    "method3_equation14": method3_record,
                }
                compact_cut_payload[variation] = {
                    "method1": method1_record,
                    "method2": method2_record,
                    "method3": method3_record,
                    "packing_fraction_used_by_method3": pf_period_record,
                }

                row: dict[str, Any] = {
                    "period": period,
                    "period_label": PERIOD_LABELS[period],
                    "x_index": x_index,
                    "t_index": t_index,
                    "bin_number": b + 1,
                    "xB_min": XB_BINS[x_index][0],
                    "xB_max": XB_BINS[x_index][1],
                    "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[t_index][0],
                    "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[t_index][1],
                    "cut_variation": variation,
                    "mx2_min_gev2": low,
                    "mx2_max_gev2": high,
                    "method1_carbon_scale_period": float(alpha_summary["central"]),
                    "method1_carbon_scale_stat_uncertainty": float(
                        alpha_summary["stat_uncertainty"]
                    ),
                }
                for target in TARGETS:
                    row[f"count_{target}"] = counts_payload[target]
                    row[f"relative_charge_fraction_{target}"] = charge_fractions[period][target]
                    row[f"rate_{target}_relative_units"] = rates_payload[target]
                # endfor
                for label, record in (
                    ("method1", method1_record),
                    ("method2", method2_record),
                    ("method3", method3_record),
                    ("packing_fraction_bin", pf_bin_record),
                    ("packing_fraction_period", pf_period_record),
                ):
                    row[f"{label}_value"] = record["value"]
                    row[f"{label}_stat_uncertainty"] = record["stat_uncertainty"]
                    row[f"{label}_bootstrap_p16"] = record["bootstrap_p16"]
                    row[f"{label}_bootstrap_p50"] = record["bootstrap_p50"]
                    row[f"{label}_bootstrap_p84"] = record["bootstrap_p84"]
                    row[f"{label}_valid_replica_fraction"] = record[
                        "valid_replica_fraction"
                    ]
                    row[f"{label}_flags"] = ";".join(record["flags"])
                # endfor
                flat_rows.append(row)
            # endfor

            bin_record = {
                "x_index": x_index,
                "t_index": t_index,
                "bin_number": b + 1,
                "xB_min": XB_BINS[x_index][0],
                "xB_max": XB_BINS[x_index][1],
                "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[t_index][0],
                "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[t_index][1],
                "exclusivity_fit": {
                    "mu_gev2": cut_entry.mu_gev2,
                    "mu_error_gev2": cut_entry.mu_error_gev2,
                    "sigma_gev2": cut_entry.sigma_gev2,
                    "sigma_error_gev2": cut_entry.sigma_error_gev2,
                    "classification": cut_entry.classification,
                },
                "cuts": cut_payload,
            }
            period_bins.append(bin_record)
            compact_bins.append(
                {
                    "x_index": x_index,
                    "t_index": t_index,
                    "bin_number": b + 1,
                    "xB_min": XB_BINS[x_index][0],
                    "xB_max": XB_BINS[x_index][1],
                    "minus_tprime_min_gev2": MINUS_TPRIME_BINS_GEV2[t_index][0],
                    "minus_tprime_max_gev2": MINUS_TPRIME_BINS_GEV2[t_index][1],
                    "cuts": compact_cut_payload,
                }
            )
        # endfor

        master_periods[period] = {
            "charge_fractions": charge_fractions[period],
            "method1_period_carbon_scale": scalar_record(alpha_summary, ()),
            "period_packing_fraction_by_cut": {
                variation: scalar_record(
                    summaries["packing_fraction_period"], (cut_index,)
                )
                for cut_index, variation in enumerate(CUT_VARIATIONS)
            },
            "bins": period_bins,
        }
        compact_periods[period] = {
            "charge_fractions": charge_fractions[period],
            "recommended_nominal_method": "method2_equation10_pending_result_review",
            "period_packing_fraction_by_cut": master_periods[period][
                "period_packing_fraction_by_cut"
            ],
            "bins": compact_bins,
        }
    # endfor

    frame = pd.DataFrame(flat_rows)
    return frame, master_periods, compact_periods, summaries_by_period


def write_covariance_products(
    output_dir: Path,
    bootstrap_results: dict[str, dict[str, np.ndarray]],
) -> dict[str, Any]:
    ensure_directory(output_dir)
    manifest: dict[str, Any] = {}
    for period in PERIODS:
        manifest[period] = {}
        for method_key in ("method1", "method2", "method3", "packing_fraction_bin"):
            manifest[period][method_key] = {}
            samples = bootstrap_results[period][method_key]
            for cut_index, variation in enumerate(CUT_VARIATIONS):
                cut_samples = samples[:, :, cut_index]
                covariance = pairwise_covariance(cut_samples)
                correlation = covariance_to_correlation(covariance)
                stem = f"{method_key}_{period}_{variation}_v1"
                cov_path = output_dir / f"{stem}_covariance.json"
                corr_path = output_dir / f"{stem}_correlation.json"
                write_json(
                    cov_path,
                    {
                        "period": period,
                        "quantity": method_key,
                        "cut_variation": variation,
                        "bin_order": list(range(1, NUMBER_OF_BINS + 1)),
                        "matrix": covariance,
                    },
                )
                write_json(
                    corr_path,
                    {
                        "period": period,
                        "quantity": method_key,
                        "cut_variation": variation,
                        "bin_order": list(range(1, NUMBER_OF_BINS + 1)),
                        "matrix": correlation,
                    },
                )
                manifest[period][method_key][variation] = {
                    "covariance": str(cov_path),
                    "correlation": str(corr_path),
                }
            # endfor
        # endfor
    # endfor
    return manifest


# =============================================================================
# Plotting
# =============================================================================

def reshape_kinematic(values: np.ndarray) -> np.ndarray:
    return np.asarray(values).reshape(len(XB_BINS), len(MINUS_TPRIME_BINS_GEV2))


def annotate_heatmap(ax: plt.Axes, matrix: np.ndarray, fmt: str = ".3f") -> None:
    for x_index in range(matrix.shape[0]):
        for t_index in range(matrix.shape[1]):
            value = matrix[x_index, t_index]
            text = "nan" if not np.isfinite(value) else format(value, fmt)
            ax.text(t_index, x_index, text, ha="center", va="center", fontsize=8)
        # endfor
    # endfor


def configure_heatmap_axes(ax: plt.Axes) -> None:
    ax.set_xticks(range(len(MINUS_TPRIME_BINS_GEV2)))
    ax.set_xticklabels(
        [f"{low:.2f}–{high:.2f}" for low, high in MINUS_TPRIME_BINS_GEV2],
        rotation=35,
        ha="right",
    )
    ax.set_yticks(range(len(XB_BINS)))
    ax.set_yticklabels([f"{low:.2f}–{high:.2f}" for low, high in XB_BINS])
    ax.set_xlabel(r"$-t^{\prime}$ (GeV$^2$)")
    ax.set_ylabel("$x_B$")


def plot_method_heatmaps(
    output_dir: Path,
    period: str,
    summaries: dict[str, Any],
) -> list[str]:
    ensure_directory(output_dir)
    paths: list[str] = []
    for method_key, title in (
        ("method1", "Method 1: carbon subtraction"),
        ("method2", "Method 2: Eq. (10)"),
        ("method3", "Method 3: Eq. (14)"),
        ("packing_fraction_bin", "Packing fraction: Eq. (11) bin diagnostic"),
    ):
        for cut_index, variation in enumerate(CUT_VARIATIONS):
            central = reshape_kinematic(summaries[method_key]["central"][:, cut_index])
            uncertainty = reshape_kinematic(
                summaries[method_key]["stat_uncertainty"][:, cut_index]
            )
            fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
            image = axes[0].imshow(central, aspect="auto", origin="upper")
            annotate_heatmap(axes[0], central)
            configure_heatmap_axes(axes[0])
            axes[0].set_title("Central value")
            fig.colorbar(image, ax=axes[0])

            image = axes[1].imshow(uncertainty, aspect="auto", origin="upper")
            annotate_heatmap(axes[1], uncertainty)
            configure_heatmap_axes(axes[1])
            axes[1].set_title("Statistical uncertainty")
            fig.colorbar(image, ax=axes[1])

            fig.suptitle(f"{PERIOD_LABELS[period]} — {title} — {variation}")
            fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
            path = output_dir / f"{method_key}_{period}_{variation}_heatmaps_v1.png"
            fig.savefig(path, dpi=180)
            plt.close(fig)
            paths.append(str(path))
        # endfor
    # endfor
    return paths


def plot_three_method_comparison(
    output_dir: Path,
    period: str,
    summaries: dict[str, Any],
) -> str:
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    offsets = {"method1": -0.18, "method2": 0.0, "method3": 0.18}
    labels = {"method1": "Method 1", "method2": "Method 2, Eq. (10)", "method3": "Method 3, Eq. (14)"}

    fig, axes = plt.subplots(3, 1, figsize=(17, 14), sharex=True)
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        ax = axes[cut_index]
        for method_key in ("method1", "method2", "method3"):
            ax.errorbar(
                bins + offsets[method_key],
                summaries[method_key]["central"][:, cut_index],
                yerr=summaries[method_key]["stat_uncertainty"][:, cut_index],
                marker="o",
                linestyle="none",
                markersize=4,
                capsize=2,
                label=labels[method_key],
            )
        # endfor
        ax.axhline(0.0, linestyle=":", linewidth=1.0)
        ax.axhline(1.0, linestyle=":", linewidth=1.0)
        ax.set_ylabel("Dilution factor")
        ax.set_title(variation.capitalize())
        ax.grid(alpha=0.25)
        ax.legend(ncol=3)
    # endfor
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(bins)
    fig.suptitle(f"{PERIOD_LABELS[period]} dilution-factor method comparison")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    path = output_dir / f"three_method_comparison_{period}_v1.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_packing_fraction_summary(
    output_dir: Path,
    period: str,
    summaries: dict[str, Any],
) -> str:
    ensure_directory(output_dir)
    bins = np.arange(1, NUMBER_OF_BINS + 1)
    fig, axes = plt.subplots(3, 1, figsize=(17, 14), sharex=True)
    for cut_index, variation in enumerate(CUT_VARIATIONS):
        ax = axes[cut_index]
        central = summaries["packing_fraction_bin"]["central"][:, cut_index]
        uncertainty = summaries["packing_fraction_bin"]["stat_uncertainty"][:, cut_index]
        period_value = summaries["packing_fraction_period"]["central"][cut_index]
        period_uncertainty = summaries["packing_fraction_period"]["stat_uncertainty"][cut_index]
        ax.errorbar(
            bins,
            central,
            yerr=uncertainty,
            marker="o",
            linestyle="none",
            markersize=4,
            capsize=2,
            label="Per-bin Eq. (11)",
        )
        ax.axhline(period_value, linewidth=1.2, label="Period-integrated Eq. (11)")
        ax.axhspan(
            period_value - period_uncertainty,
            period_value + period_uncertainty,
            alpha=0.15,
        )
        ax.set_ylabel("Packing fraction")
        ax.set_title(variation.capitalize())
        ax.grid(alpha=0.25)
        ax.legend()
    # endfor
    axes[-1].set_xlabel("Combined kinematic-bin number")
    axes[-1].set_xticks(bins)
    fig.suptitle(f"{PERIOD_LABELS[period]} packing-fraction diagnostics")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    path = output_dir / f"packing_fraction_summary_{period}_v1.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def plot_method1_control_summary(
    output_dir: Path,
    central_results: dict[str, dict[str, Any]],
    summaries_by_period: dict[str, Any],
) -> str:
    ensure_directory(output_dir)
    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(PERIODS))
    values = np.array([central_results[p]["method1_alpha"] for p in PERIODS], dtype=float)
    errors = np.array(
        [summaries_by_period[p]["method1_alpha"]["stat_uncertainty"] for p in PERIODS],
        dtype=float,
    )
    ax.errorbar(x, values, yerr=errors, marker="o", linestyle="none", capsize=3)
    ax.set_xticks(x)
    ax.set_xticklabels([PERIOD_LABELS[p] for p in PERIODS])
    ax.set_ylabel("Raw-count NH$_3$/C normalization")
    ax.set_title(r"Method-1 period-wide carbon normalization" "\n" r"$0.20 \leq M_x^2 < 0.40$ GeV$^2$")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    path = output_dir / "method1_period_carbon_scales_v1.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return str(path)


def write_plots(
    output_dir: Path,
    central_results: dict[str, dict[str, Any]],
    summaries_by_period: dict[str, Any],
) -> list[str]:
    paths: list[str] = []
    paths.append(
        plot_method1_control_summary(output_dir, central_results, summaries_by_period)
    )
    for period in PERIODS:
        period_dir = output_dir / period
        paths.extend(plot_method_heatmaps(period_dir, period, summaries_by_period[period]))
        paths.append(
            plot_three_method_comparison(period_dir, period, summaries_by_period[period])
        )
        paths.append(
            plot_packing_fraction_summary(period_dir, period, summaries_by_period[period])
        )
    # endfor
    return paths


# =============================================================================
# Command-line interface and main program
# =============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Determine RGC enpi+ dilution factors with three methods."
    )
    parser.add_argument(
        "--tree",
        default=DEFAULT_TREE_NAME,
        help=f"ROOT tree name (default: {DEFAULT_TREE_NAME}).",
    )
    parser.add_argument(
        "--exclusivity-json",
        type=Path,
        default=DEFAULT_CUT_JSON,
        help=(
            "Final carbon-assisted Mx2-cut JSON.  Default is the sibling "
            "channel_selection output path."
        ),
    )
    parser.add_argument(
        "--input",
        action="append",
        default=[],
        metavar="PERIOD:TARGET=PATH",
        help=(
            "Override one ROOT input. May be repeated. Example: "
            "--input su22:NH3=/path/file.root"
        ),
    )
    parser.add_argument(
        "--charge-fractions-json",
        type=Path,
        default=None,
        help="Optional JSON overriding the built-in relative charge fractions.",
    )
    parser.add_argument(
        "--control-min",
        type=float,
        default=DEFAULT_CONTROL_MIN_GEV2,
        help="Method-1 carbon control-window minimum (GeV^2).",
    )
    parser.add_argument(
        "--control-max",
        type=float,
        default=DEFAULT_CONTROL_MAX_GEV2,
        help="Method-1 carbon control-window maximum (GeV^2).",
    )
    parser.add_argument(
        "--replicas",
        type=int,
        default=DEFAULT_REPLICAS,
        help=f"Poisson replicas per period (default: {DEFAULT_REPLICAS}).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Base random seed (default: {DEFAULT_SEED}).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=MAXIMUM_WORKERS,
        help="Worker processes; hard-capped at 7 (default: 7).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Stable output directory (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Write numerical products without rendering plots.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    if not (math.isfinite(args.control_min) and math.isfinite(args.control_max)):
        raise RuntimeError("Control-window endpoints must be finite.")
    # endif
    if args.control_min >= args.control_max:
        raise RuntimeError("Control-window minimum must be below the maximum.")
    # endif
    if args.replicas < 2:
        raise RuntimeError("At least two replicas are required.")
    # endif

    workers = max(1, min(int(args.workers), MAXIMUM_WORKERS, os.cpu_count() or 1))
    output_dir = args.output_dir.resolve()
    tables_dir = output_dir / "tables"
    covariance_dir = output_dir / "covariance"
    plots_dir = output_dir / "plots"
    for directory in (output_dir, tables_dir, covariance_dir, plots_dir):
        ensure_directory(directory)
    # endfor

    input_paths = {
        period: dict(targets)
        for period, targets in DEFAULT_INPUTS.items()
    }
    for override_text in args.input:
        period, target, path = parse_input_override(override_text)
        input_paths[period][target] = path
    # endfor

    charge_fractions = load_charge_fractions(args.charge_fractions_json)
    cut_json = args.exclusivity_json.resolve()
    cuts = load_exclusivity_cuts(cut_json)

    print("=" * 78)
    print("RGC exclusive enpi+ dilution-factor determination")
    print("=" * 78)
    print(f"Working directory:       {Path.cwd()}")
    print(f"Exclusivity JSON:        {cut_json}")
    print(f"Output directory:        {output_dir}")
    print(f"Tree:                    {args.tree}")
    print(f"Workers:                 {workers} (hard maximum {MAXIMUM_WORKERS})")
    print(f"Poisson replicas/period: {args.replicas:,}")
    print(
        f"Method-1 control region: [{args.control_min:.6g}, "
        f"{args.control_max:.6g}) GeV^2"
    )
    print()

    loaded = load_all_datasets(input_paths, args.tree, workers)
    pattern_counts, observed_counts, count_rows = build_count_arrays(
        loaded,
        cuts,
        args.control_min,
        args.control_max,
    )
    count_frame = pd.DataFrame(count_rows)
    count_path = tables_dir / "target_counts_all_selections_v1.csv"
    count_frame.to_csv(count_path, index=False)

    central_results = {
        period: observed_estimators_for_period(
            observed_counts[p_index],
            period,
            charge_fractions,
        )
        for p_index, period in enumerate(PERIODS)
    }

    bootstrap_results = run_bootstrap(
        pattern_counts,
        charge_fractions,
        args.replicas,
        args.seed,
        workers,
    )

    flat_frame, master_periods, compact_periods, summaries_by_period = build_output_tables(
        observed_counts,
        central_results,
        bootstrap_results,
        charge_fractions,
        cuts,
    )
    flat_csv_path = tables_dir / "dilution_factors_all_methods_v1.csv"
    flat_frame.to_csv(flat_csv_path, index=False)

    covariance_manifest = write_covariance_products(
        covariance_dir,
        bootstrap_results,
    )

    plot_paths: list[str] = []
    if not args.skip_plots:
        plot_paths = write_plots(plots_dir, central_results, summaries_by_period)
    # endif

    provenance = {
        "script": Path(__file__).name,
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "working_directory": str(Path.cwd()),
        "tree_name": args.tree,
        "maximum_worker_cap": MAXIMUM_WORKERS,
        "workers_used": workers,
        "bootstrap_replicas_per_period": args.replicas,
        "bootstrap_seed": args.seed,
        "exclusivity_json": str(cut_json),
        "exclusivity_json_sha256": sha256_file(cut_json),
        "control_region_gev2": [args.control_min, args.control_max],
        "root_inputs": {
            period: {
                target: {
                    "path": loaded[(period, target)].file_path,
                    "branches": {
                        "xB": loaded[(period, target)].x_branch,
                        "tprime": loaded[(period, target)].tprime_branch,
                        "Mx2": loaded[(period, target)].mx2_branch,
                    },
                }
                for target in TARGETS
            }
            for period in PERIODS
        },
        "charge_fractions": charge_fractions,
        "equation_coefficients": EQ10_EQ11_EQ14_COEFFICIENTS,
        "method1_definition": (
            "Exact channel-selection raw-count normalization pooled over all "
            "24 bins in the control window; no additional charge normalization."
        ),
        "method2_definition": "Standalone dilution-factor note Eq. (10).",
        "method3_definition": (
            "Per-bin Eq. (11) diagnostics; period-integrated Eq. (11) packing "
            "fraction inserted into nonlinear Eq. (14)."
        ),
        "statistical_model": (
            "Independent Poisson replicas of 16 disjoint membership-pattern "
            "counts for control/tight/nominal/loose selections."
        ),
    }

    master_payload = {
        "schema_version": 1,
        "analysis": "RGC exclusive enpi+ dilution-factor determination",
        "status": (
            "All three methods are retained. Eq. (10) is the provisional "
            "nominal candidate pending inspection of the results."
        ),
        "binning": {
            "xB": [list(interval) for interval in XB_BINS],
            "minus_tprime_gev2": [list(interval) for interval in MINUS_TPRIME_BINS_GEV2],
            "combined_bin_order": (
                "x_index major, t_index minor, one-based bin_number"
            ),
        },
        "cut_variations": list(CUT_VARIATIONS),
        "provenance": provenance,
        "periods": master_periods,
        "covariance_manifest": covariance_manifest,
        "plot_paths": plot_paths,
    }
    master_json_path = output_dir / "dilution_factors_v1.json"
    write_json(master_json_path, master_payload)

    compact_payload = {
        "schema_version": 1,
        "analysis": "RGC exclusive enpi+ dilution factors for downstream asymmetries",
        "provisional_nominal_method": "method2_equation10",
        "note": (
            "Loose, nominal, and tight values are all carried downstream. "
            "No exclusivity-cut systematic is assigned at this stage."
        ),
        "binning": master_payload["binning"],
        "periods": compact_periods,
        "source_master_json": str(master_json_path),
        "source_exclusivity_json": str(cut_json),
        "source_exclusivity_json_sha256": provenance["exclusivity_json_sha256"],
    }
    compact_json_path = output_dir / "dilution_factors_production_v1.json"
    write_json(compact_json_path, compact_payload)

    configuration_path = output_dir / "configuration_v1.json"
    write_json(configuration_path, provenance)

    print()
    print("Dilution-factor calculation complete.")
    print(f"  Master JSON:       {master_json_path}")
    print(f"  Downstream JSON:   {compact_json_path}")
    print(f"  Flat CSV:          {flat_csv_path}")
    print(f"  Count CSV:         {count_path}")
    print(f"  Covariance dir:    {covariance_dir}")
    if not args.skip_plots:
        print(f"  Plot directory:    {plots_dir}")
    # endif
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted by user.", file=sys.stderr)
        raise SystemExit(130)
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
