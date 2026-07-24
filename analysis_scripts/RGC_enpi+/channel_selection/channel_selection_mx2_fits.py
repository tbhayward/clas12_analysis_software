#!/usr/bin/env python3
"""
channel_selection_mx2_fits_v5.py

Deterministic pass-1 fit-stability study of the missing-neutron Mx2 peak used for the RGC
exclusive e pi+ channel selection.

For each of the three RGC run periods, before and after the finalized momentum
corrections, the script:

  1. reads xB, tprime, and Mx2 from the PhysicsEvents ROOT tree;
  2. divides the sample into the fixed 4 xB by 6 (-tprime) bins;
  3. fits each Mx2 distribution with a Gaussian signal and four deterministic
     background hypotheses:

         linear polynomial,
         quadratic polynomial,
         cubic polynomial,
         quartic polynomial;

     The script also recommends the lowest polynomial order from quadratic
     through quartic that gives an acceptable signal-region fit, with the
     signal region defined from the fitted Gaussian as mu +/- 2 sigma;

  4. extracts the Gaussian mean, width, and signal yield;
  5. records fit quality, covariance information, model differences, and
     quality-control flags;
  6. writes compact JSON, flat CSV, LaTeX tables, and diagnostic plots.

This is intentionally the deterministic first pass. It does not perform the
Poisson-replica study, Carbon subtraction, or the Mx2-window asymmetry
systematic. Those are later stages.

The branch tprime is negated explicitly before applying the requested -tprime
binning.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import uproot


# =============================================================================
# Physics constants and fixed analysis binning
# =============================================================================

NEUTRON_MASS_GEV = 0.9395654133
NEUTRON_MASS2_GEV2 = NEUTRON_MASS_GEV**2

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

DEFAULT_INPUTS: dict[str, dict[str, Path]] = {
    "su22": {
        "before": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "enpi+/rgc_su22_inb_NH3_epi+_2.root"
        ),
        "after": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "paper_versions/rgc_su22_inb_NH3_epi+_mom_corrections.root"
        ),
    },
    "fa22": {
        "before": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "enpi+/rgc_fa22_inb_NH3_epi+_2.root"
        ),
        "after": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "paper_versions/rgc_fa22_inb_NH3_epi+_mom_corrections.root"
        ),
    },
    "sp23": {
        "before": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "enpi+/rgc_sp23_inb_NH3_epi+_2.root"
        ),
        "after": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "paper_versions/rgc_sp23_inb_NH3_epi+_mom_corrections.root"
        ),
    },
}

PERIOD_LABELS = {
    "su22": "Su22",
    "fa22": "Fa22",
    "sp23": "Sp23",
}

STAGE_LABELS = {
    "before": "Before momentum corrections",
    "after": "After momentum corrections",
}

BACKGROUND_LABELS = {
    "linear": "Linear background",
    "quadratic": "Quadratic background",
    "cubic": "Cubic background",
    "quartic": "Quartic background",
}

BACKGROUND_MODELS: tuple[str, ...] = (
    "linear",
    "quadratic",
    "cubic",
    "quartic",
)

ALTERNATIVE_BACKGROUND_MODELS: tuple[str, ...] = tuple(
    model for model in BACKGROUND_MODELS if model != "quadratic"
)

NOMINAL_BACKGROUND_MODEL = "quadratic"

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


# =============================================================================
# Data containers
# =============================================================================

@dataclass(frozen=True)
class InputDataset:
    period: str
    stage: str
    file_path: str
    tree_name: str
    x_branch: str
    tprime_branch: str
    mx2_branch: str


@dataclass(frozen=True)
class FitJob:
    period: str
    stage: str
    x_index: int
    t_index: int
    x_min: float
    x_max: float
    minus_tprime_min_gev2: float
    minus_tprime_max_gev2: float
    values: np.ndarray
    histogram_min_gev2: float
    histogram_max_gev2: float
    histogram_bins: int
    fit_min_gev2: float
    fit_max_gev2: float
    minimum_events: int
    corrected_mean_max_gev2: float


# =============================================================================
# General helpers
# =============================================================================

def json_safe(value: Any) -> Any:
    """Recursively convert NumPy objects and non-finite floats for JSON."""
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
    if isinstance(value, (np.floating, float)):
        numeric = float(value)
        return numeric if math.isfinite(numeric) else None
    # endif
    if isinstance(value, Path):
        return str(value)
    # endif
    return value


def ensure_directory(path: Path) -> None:
    """Create a directory and all missing parents."""
    path.mkdir(parents=True, exist_ok=True)


def bin_identifier(x_index: int, t_index: int) -> str:
    """Return a stable one-based bin identifier."""
    return f"x{x_index + 1}_t{t_index + 1}"


def combined_bin_number(x_index: int, t_index: int) -> int:
    """Map a 4 x 6 bin coordinate to a one-based sequential bin number."""
    return x_index * len(MINUS_TPRIME_BINS_GEV2) + t_index + 1


def parse_dataset_override(text: str) -> tuple[str, str, Path]:
    """Parse PERIOD:STAGE=/path/to/file.root."""
    if "=" not in text or ":" not in text.split("=", 1)[0]:
        raise argparse.ArgumentTypeError(
            "Dataset overrides must have the form PERIOD:STAGE=/path/file.root."
        )
    # endif

    key_text, path_text = text.split("=", 1)
    period, stage = (piece.strip().lower() for piece in key_text.split(":", 1))
    if period not in DEFAULT_INPUTS:
        raise argparse.ArgumentTypeError(
            f"Unknown period '{period}'. Expected one of {sorted(DEFAULT_INPUTS)}."
        )
    # endif
    if stage not in {"before", "after"}:
        raise argparse.ArgumentTypeError(
            f"Unknown stage '{stage}'. Expected 'before' or 'after'."
        )
    # endif
    return period, stage, Path(path_text.strip()).expanduser()


def parse_branch_override(text: str) -> tuple[str, str]:
    """Parse LOGICAL=ROOT_BRANCH."""
    if "=" not in text:
        raise argparse.ArgumentTypeError(
            "Branch overrides must have the form LOGICAL=ROOT_BRANCH."
        )
    # endif

    logical, physical = (piece.strip() for piece in text.split("=", 1))
    if logical not in BRANCH_ALIASES:
        raise argparse.ArgumentTypeError(
            f"Unknown logical branch '{logical}'. Expected one of "
            f"{sorted(BRANCH_ALIASES)}."
        )
    # endif
    if not physical:
        raise argparse.ArgumentTypeError("ROOT branch name cannot be empty.")
    # endif
    return logical, physical


def find_tree(file_path: Path, requested_tree: str) -> str:
    """Resolve a ROOT tree name, including ROOT cycle suffixes."""
    with uproot.open(file_path) as root_file:
        key_map = {str(key).split(";")[0]: str(key) for key in root_file.keys()}
        if requested_tree in key_map:
            return key_map[requested_tree]
        # endif

        tree_like: list[str] = []
        for key in root_file.keys():
            try:
                obj = root_file[key]
                if hasattr(obj, "arrays"):
                    tree_like.append(str(key))
                # endif
            except Exception:
                continue
            # endtry
        # endfor
    # endwith

    if len(tree_like) == 1:
        return tree_like[0]
    # endif

    raise KeyError(
        f"Tree '{requested_tree}' was not found in {file_path}. "
        f"Tree-like objects: {tree_like}"
    )


def resolve_branches(
    file_path: Path,
    tree_name: str,
    branch_overrides: dict[str, str],
) -> dict[str, str]:
    """Resolve required logical branches to physical ROOT branch names."""
    with uproot.open(file_path) as root_file:
        tree = root_file[tree_name]
        available = {str(key).split(";")[0] for key in tree.keys()}
    # endwith

    resolved: dict[str, str] = {}
    for logical, aliases in BRANCH_ALIASES.items():
        if logical in branch_overrides:
            physical = branch_overrides[logical]
            if physical not in available:
                raise KeyError(
                    f"Requested branch override {logical}={physical} is absent "
                    f"from {file_path}:{tree_name}."
                )
            # endif
            resolved[logical] = physical
            continue
        # endif

        match = next((alias for alias in aliases if alias in available), None)
        if match is None:
            raise KeyError(
                f"Could not resolve logical branch '{logical}' in "
                f"{file_path}:{tree_name}. Tried {aliases}."
            )
        # endif
        resolved[logical] = match
    # endfor
    return resolved


def preflight_inputs(
    input_paths: dict[str, dict[str, Path]],
    requested_tree: str,
    branch_overrides: dict[str, str],
) -> list[InputDataset]:
    """Validate all six inputs and resolve tree and branch names."""
    datasets: list[InputDataset] = []
    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            file_path = input_paths[period][stage]
            if not file_path.is_file():
                raise FileNotFoundError(
                    f"Missing {period} {stage} input file: {file_path}"
                )
            # endif

            tree_name = find_tree(file_path, requested_tree)
            branches = resolve_branches(
                file_path=file_path,
                tree_name=tree_name,
                branch_overrides=branch_overrides,
            )
            datasets.append(
                InputDataset(
                    period=period,
                    stage=stage,
                    file_path=str(file_path),
                    tree_name=tree_name,
                    x_branch=branches["xB"],
                    tprime_branch=branches["tprime"],
                    mx2_branch=branches["Mx2"],
                )
            )
        # endfor
    # endfor
    return datasets


def load_dataset_arrays(
    dataset: InputDataset,
    step_size: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Read and minimally preselect xB, -tprime, and Mx2 from one ROOT file."""
    x_chunks: list[np.ndarray] = []
    minus_tprime_chunks: list[np.ndarray] = []
    mx2_chunks: list[np.ndarray] = []

    entries_read = 0
    finite_entries = 0
    kinematic_entries = 0

    expressions = [
        dataset.x_branch,
        dataset.tprime_branch,
        dataset.mx2_branch,
    ]

    source = f"{dataset.file_path}:{dataset.tree_name}"
    for arrays in uproot.iterate(
        source,
        expressions=expressions,
        library="np",
        step_size=step_size,
    ):
        x_b = np.asarray(arrays[dataset.x_branch], dtype=float)
        tprime = np.asarray(arrays[dataset.tprime_branch], dtype=float)
        mx2 = np.asarray(arrays[dataset.mx2_branch], dtype=float)
        minus_tprime = -tprime

        entries_read += int(x_b.size)
        finite_mask = (
            np.isfinite(x_b)
            & np.isfinite(minus_tprime)
            & np.isfinite(mx2)
        )
        finite_entries += int(np.count_nonzero(finite_mask))

        kinematic_mask = (
            finite_mask
            & (x_b >= XB_BINS[0][0])
            & (x_b < XB_BINS[-1][1])
            & (minus_tprime >= MINUS_TPRIME_BINS_GEV2[0][0])
            & (minus_tprime < MINUS_TPRIME_BINS_GEV2[-1][1])
        )
        kinematic_entries += int(np.count_nonzero(kinematic_mask))

        if not np.any(kinematic_mask):
            continue
        # endif

        x_chunks.append(x_b[kinematic_mask])
        minus_tprime_chunks.append(minus_tprime[kinematic_mask])
        mx2_chunks.append(mx2[kinematic_mask])
    # endfor

    if x_chunks:
        x_all = np.concatenate(x_chunks)
        minus_tprime_all = np.concatenate(minus_tprime_chunks)
        mx2_all = np.concatenate(mx2_chunks)
    else:
        x_all = np.asarray([], dtype=float)
        minus_tprime_all = np.asarray([], dtype=float)
        mx2_all = np.asarray([], dtype=float)
    # endif

    counts = {
        "entries_read": entries_read,
        "finite_entries": finite_entries,
        "kinematic_entries": kinematic_entries,
    }
    return x_all, minus_tprime_all, mx2_all, counts


# =============================================================================
# Signal and background models
# =============================================================================

def gaussian_signal_from_yield(
    x: np.ndarray,
    signal_yield: float,
    mean_gev2: float,
    sigma_gev2: float,
    histogram_bin_width_gev2: float,
) -> np.ndarray:
    """Gaussian counts per histogram bin, parameterized by total area."""
    sigma_safe = max(float(sigma_gev2), 1.0e-12)
    normalization = (
        float(signal_yield)
        * float(histogram_bin_width_gev2)
        / (sigma_safe * math.sqrt(2.0 * math.pi))
    )
    return normalization * np.exp(
        -0.5 * ((np.asarray(x, dtype=float) - mean_gev2) / sigma_safe) ** 2
    )


def linear_background(
    x: np.ndarray,
    c0: float,
    c1: float,
) -> np.ndarray:
    """Linear background centered at the neutron mass squared."""
    dx = np.asarray(x, dtype=float) - NEUTRON_MASS2_GEV2
    return c0 + c1 * dx


def quadratic_background(
    x: np.ndarray,
    c0: float,
    c1: float,
    c2: float,
) -> np.ndarray:
    """Quadratic background centered at the neutron mass squared."""
    dx = np.asarray(x, dtype=float) - NEUTRON_MASS2_GEV2
    return c0 + c1 * dx + c2 * dx**2


def cubic_background(
    x: np.ndarray,
    c0: float,
    c1: float,
    c2: float,
    c3: float,
) -> np.ndarray:
    """Cubic background centered at the neutron mass squared."""
    dx = np.asarray(x, dtype=float) - NEUTRON_MASS2_GEV2
    return c0 + c1 * dx + c2 * dx**2 + c3 * dx**3


def quartic_background(
    x: np.ndarray,
    c0: float,
    c1: float,
    c2: float,
    c3: float,
    c4: float,
) -> np.ndarray:
    """Quartic background centered at the neutron mass squared."""
    dx = np.asarray(x, dtype=float) - NEUTRON_MASS2_GEV2
    return c0 + c1 * dx + c2 * dx**2 + c3 * dx**3 + c4 * dx**4



def build_total_model(
    background_model: str,
    histogram_bin_width_gev2: float,
) -> Callable[..., np.ndarray]:
    """Create one Gaussian-plus-background model with a fixed bin width."""
    if background_model == "linear":
        def model(
            x: np.ndarray,
            signal_yield: float,
            mean_gev2: float,
            sigma_gev2: float,
            c0: float,
            c1: float,
        ) -> np.ndarray:
            return gaussian_signal_from_yield(
                x,
                signal_yield,
                mean_gev2,
                sigma_gev2,
                histogram_bin_width_gev2,
            ) + linear_background(x, c0, c1)
        # enddef
        return model
    # endif

    if background_model == "quadratic":
        def model(
            x: np.ndarray,
            signal_yield: float,
            mean_gev2: float,
            sigma_gev2: float,
            c0: float,
            c1: float,
            c2: float,
        ) -> np.ndarray:
            return gaussian_signal_from_yield(
                x,
                signal_yield,
                mean_gev2,
                sigma_gev2,
                histogram_bin_width_gev2,
            ) + quadratic_background(x, c0, c1, c2)
        # enddef
        return model
    # endif

    if background_model == "cubic":
        def model(
            x: np.ndarray,
            signal_yield: float,
            mean_gev2: float,
            sigma_gev2: float,
            c0: float,
            c1: float,
            c2: float,
            c3: float,
        ) -> np.ndarray:
            return gaussian_signal_from_yield(
                x,
                signal_yield,
                mean_gev2,
                sigma_gev2,
                histogram_bin_width_gev2,
            ) + cubic_background(x, c0, c1, c2, c3)
        # enddef
        return model
    # endif

    if background_model == "quartic":
        def model(
            x: np.ndarray,
            signal_yield: float,
            mean_gev2: float,
            sigma_gev2: float,
            c0: float,
            c1: float,
            c2: float,
            c3: float,
            c4: float,
        ) -> np.ndarray:
            return gaussian_signal_from_yield(
                x,
                signal_yield,
                mean_gev2,
                sigma_gev2,
                histogram_bin_width_gev2,
            ) + quartic_background(x, c0, c1, c2, c3, c4)
        # enddef
        return model
    # endif



    raise ValueError(f"Unknown background model: {background_model}")


def evaluate_background(
    background_model: str,
    x: np.ndarray,
    parameters: np.ndarray,
) -> np.ndarray:
    """Evaluate only the fitted background component."""
    if background_model == "linear":
        return linear_background(x, *parameters[3:5])
    # endif
    if background_model == "quadratic":
        return quadratic_background(x, *parameters[3:6])
    # endif
    if background_model == "cubic":
        return cubic_background(x, *parameters[3:7])
    # endif
    if background_model == "quartic":
        return quartic_background(x, *parameters[3:8])
    # endif
    raise ValueError(background_model)


def initial_values_and_bounds(
    background_model: str,
    fit_x: np.ndarray,
    fit_y: np.ndarray,
    histogram_bin_width_gev2: float,
    fit_min_gev2: float,
    fit_max_gev2: float,
    stage: str,
    corrected_mean_max_gev2: float,
) -> tuple[list[float], tuple[list[float], list[float]]]:
    """Construct stable initial values and bounded parameter ranges."""
    if fit_x.size == 0:
        raise ValueError("Cannot initialize a fit with no fit-range bins.")
    # endif

    baseline = max(float(np.percentile(fit_y, 20.0)), 0.0)
    peak_window = (
        (fit_x >= max(fit_min_gev2, NEUTRON_MASS2_GEV2 - 0.16))
        & (fit_x <= min(fit_max_gev2, NEUTRON_MASS2_GEV2 + 0.16))
    )

    if np.any(peak_window):
        local_y = fit_y[peak_window]
        peak_height = max(float(np.max(local_y) - baseline), 1.0)
    else:
        peak_height = max(float(np.max(fit_y) - baseline), 1.0)
    # endif

    # Initialize every Gaussian at the physical neutron mass squared rather
    # than at the largest local histogram bin. The latter can be captured by
    # a broad continuum structure in the difficult low-tprime bins.
    mean_guess = NEUTRON_MASS2_GEV2

    mean_low = max(fit_min_gev2, NEUTRON_MASS2_GEV2 - 0.20)
    if stage == "after":
        mean_high = min(fit_max_gev2, corrected_mean_max_gev2)
    else:
        # Before-correction fits are diagnostic and can have genuinely shifted
        # peaks, so retain the wider historical upper bound.
        mean_high = min(fit_max_gev2, NEUTRON_MASS2_GEV2 + 0.20)
    # endif

    mean_guess = float(
        np.clip(
            mean_guess,
            mean_low + 1.0e-4,
            mean_high - 1.0e-4,
        )
    )

    sigma_guess = 0.075
    yield_guess = (
        peak_height
        * sigma_guess
        * math.sqrt(2.0 * math.pi)
        / histogram_bin_width_gev2
    )
    total_fit_counts = max(float(np.sum(fit_y)), 1.0)
    yield_upper = max(3.0 * total_fit_counts, 10.0 * yield_guess, 100.0)
    background_scale = max(float(np.max(fit_y)), 1.0)

    common_initial = [yield_guess, mean_guess, sigma_guess]
    common_lower = [0.0, mean_low, 0.010]
    common_upper = [yield_upper, mean_high, 0.180]

    if background_model == "linear":
        initial = [*common_initial, baseline, 0.0]
        lower = [
            *common_lower,
            -0.50 * background_scale,
            -30.0 * background_scale,
        ]
        upper = [
            *common_upper,
            3.00 * background_scale,
            30.0 * background_scale,
        ]
    elif background_model == "quadratic":
        initial = [*common_initial, baseline, 0.0, 0.0]
        lower = [
            *common_lower,
            -0.50 * background_scale,
            -30.0 * background_scale,
            -150.0 * background_scale,
        ]
        upper = [
            *common_upper,
            3.00 * background_scale,
            30.0 * background_scale,
            150.0 * background_scale,
        ]
    elif background_model == "cubic":
        initial = [*common_initial, baseline, 0.0, 0.0, 0.0]
        lower = [
            *common_lower,
            -0.50 * background_scale,
            -30.0 * background_scale,
            -150.0 * background_scale,
            -750.0 * background_scale,
        ]
        upper = [
            *common_upper,
            3.00 * background_scale,
            30.0 * background_scale,
            150.0 * background_scale,
            750.0 * background_scale,
        ]
    elif background_model == "quartic":
        initial = [*common_initial, baseline, 0.0, 0.0, 0.0, 0.0]
        lower = [
            *common_lower,
            -0.50 * background_scale,
            -30.0 * background_scale,
            -150.0 * background_scale,
            -750.0 * background_scale,
            -4000.0 * background_scale,
        ]
        upper = [
            *common_upper,
            3.00 * background_scale,
            30.0 * background_scale,
            150.0 * background_scale,
            750.0 * background_scale,
            4000.0 * background_scale,
        ]
    else:
        raise ValueError(background_model)
    # endif

    return initial, (lower, upper)


def parameter_at_bound(
    value: float,
    lower: float,
    upper: float,
    relative_tolerance: float = 2.0e-3,
    absolute_tolerance: float = 2.0e-4,
) -> bool:
    """Check whether a fitted parameter is numerically close to either bound."""
    scale = max(abs(lower), abs(upper), abs(value), 1.0)
    tolerance = max(absolute_tolerance, relative_tolerance * scale)
    return abs(value - lower) <= tolerance or abs(value - upper) <= tolerance


def fit_background_model(
    background_model: str,
    centers: np.ndarray,
    counts: np.ndarray,
    fit_mask: np.ndarray,
    histogram_bin_width_gev2: float,
    fit_min_gev2: float,
    fit_max_gev2: float,
    stage: str,
    corrected_mean_max_gev2: float,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Fit one histogram with a Gaussian and one background family."""
    fit_x = np.asarray(centers[fit_mask], dtype=float)
    fit_y = np.asarray(counts[fit_mask], dtype=float)
    fit_error = np.sqrt(np.maximum(fit_y, 1.0))

    model_function = build_total_model(
        background_model=background_model,
        histogram_bin_width_gev2=histogram_bin_width_gev2,
    )

    try:
        initial, bounds = initial_values_and_bounds(
            background_model=background_model,
            fit_x=fit_x,
            fit_y=fit_y,
            histogram_bin_width_gev2=histogram_bin_width_gev2,
            fit_min_gev2=fit_min_gev2,
            fit_max_gev2=fit_max_gev2,
            stage=stage,
            corrected_mean_max_gev2=corrected_mean_max_gev2,
        )
        parameters, covariance = curve_fit(
            model_function,
            fit_x,
            fit_y,
            p0=initial,
            bounds=bounds,
            sigma=fit_error,
            absolute_sigma=True,
            maxfev=200000,
        )
    except (
        RuntimeError,
        ValueError,
        np.linalg.LinAlgError,
        FloatingPointError,
        OverflowError,
    ) as exc:
        result = {
            "background_model": background_model,
            "success": False,
            "status": f"fit_failed:{type(exc).__name__}",
            "review_reasons": ["fit_failed"],
        }
        diagnostics = {
            "fit_x": fit_x,
            "fit_y": fit_y,
            "fit_error": fit_error,
            "fit_total": np.asarray([], dtype=float),
            "fit_signal": np.asarray([], dtype=float),
            "fit_background": np.asarray([], dtype=float),
            "pulls": np.asarray([], dtype=float),
        }
        return result, diagnostics
    # endtry

    parameters = np.asarray(parameters, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    diagonal = np.diag(covariance)
    errors = np.sqrt(np.where(diagonal >= 0.0, diagonal, np.nan))

    fit_total = model_function(fit_x, *parameters)
    fit_signal = gaussian_signal_from_yield(
        fit_x,
        parameters[0],
        parameters[1],
        parameters[2],
        histogram_bin_width_gev2,
    )
    fit_background = evaluate_background(
        background_model,
        fit_x,
        parameters,
    )
    residuals = fit_y - fit_total
    pulls = residuals / fit_error

    chi2 = float(np.sum(pulls**2))
    number_of_parameters = int(parameters.size)
    number_of_points = int(fit_x.size)
    ndf = number_of_points - number_of_parameters
    chi2_ndf = chi2 / ndf if ndf > 0 else math.nan

    fitted_mean = float(parameters[1])
    fitted_sigma = float(abs(parameters[2]))
    signal_region_mask = (
        (fit_x >= fitted_mean - 2.0 * fitted_sigma)
        & (fit_x <= fitted_mean + 2.0 * fitted_sigma)
    )
    signal_region_chi2 = float(np.sum(pulls[signal_region_mask] ** 2))
    signal_region_npoints = int(np.count_nonzero(signal_region_mask))

    # This is an intentionally local diagnostic. The background parameters are
    # constrained by the complete fit range, so a mathematically exact local
    # number of degrees of freedom is not uniquely defined. Subtracting the
    # three Gaussian parameters gives a transparent and reproducible
    # approximation for comparing models in the signal region.
    signal_region_ndf_approx = max(signal_region_npoints - 3, 1)
    signal_region_chi2_ndf_approx = (
        signal_region_chi2 / signal_region_ndf_approx
    )
    signal_region_chi2_per_bin = (
        signal_region_chi2 / signal_region_npoints
        if signal_region_npoints > 0
        else math.nan
    )

    if number_of_points > number_of_parameters + 1:
        aicc = (
            chi2
            + 2.0 * number_of_parameters
            + (
                2.0
                * number_of_parameters
                * (number_of_parameters + 1)
                / (number_of_points - number_of_parameters - 1)
            )
        )
    else:
        aicc = math.inf
    # endif

    signal_yield = float(parameters[0])
    signal_yield_error = float(errors[0])
    mean_gev2 = float(parameters[1])
    mean_error_gev2 = float(errors[1])
    sigma_gev2 = float(abs(parameters[2]))
    sigma_error_gev2 = float(errors[2])

    dense_x = np.linspace(fit_min_gev2, fit_max_gev2, 500)
    dense_background = evaluate_background(
        background_model,
        dense_x,
        parameters,
    )
    background_minimum = float(np.min(dense_background))

    fitted_mean = float(parameters[1])
    fitted_sigma = float(abs(parameters[2]))
    dense_signal_mask_2sigma = (
        (dense_x >= fitted_mean - 2.0 * fitted_sigma)
        & (dense_x <= fitted_mean + 2.0 * fitted_sigma)
    )
    dense_signal_mask_3sigma = (
        (dense_x >= fitted_mean - 3.0 * fitted_sigma)
        & (dense_x <= fitted_mean + 3.0 * fitted_sigma)
    )
    background_minimum_2sigma = (
        float(np.min(dense_background[dense_signal_mask_2sigma]))
        if np.any(dense_signal_mask_2sigma)
        else math.nan
    )
    background_minimum_3sigma = (
        float(np.min(dense_background[dense_signal_mask_3sigma]))
        if np.any(dense_signal_mask_3sigma)
        else math.nan
    )

    peak_signal = float(
        gaussian_signal_from_yield(
            np.asarray([mean_gev2]),
            signal_yield,
            mean_gev2,
            sigma_gev2,
            histogram_bin_width_gev2,
        )[0]
    )
    background_at_peak = float(
        evaluate_background(
            background_model,
            np.asarray([mean_gev2]),
            parameters,
        )[0]
    )
    characteristic_count_scale = max(
        float(np.max(fit_y)) if fit_y.size else 0.0,
        abs(background_at_peak),
        1.0,
    )
    background_negativity_tolerance = max(
        1.0,
        0.01 * characteristic_count_scale,
    )
    peak_significance_proxy = peak_signal / math.sqrt(
        max(peak_signal + abs(background_at_peak), 1.0)
    )

    review_reasons: list[str] = []
    if not np.all(np.isfinite(parameters)):
        review_reasons.append("nonfinite_parameters")
    # endif
    if not np.all(np.isfinite(covariance)):
        review_reasons.append("nonfinite_covariance")
    # endif
    if (
        not np.isfinite(signal_region_chi2_ndf_approx)
        or signal_region_chi2_ndf_approx > 3.0
    ):
        review_reasons.append("poor_signal_region_chi2")
    # endif
    if not np.isfinite(chi2_ndf) or chi2_ndf > 5.0:
        review_reasons.append("poor_global_chi2_diagnostic")
    # endif
    if not np.isfinite(mean_error_gev2) or mean_error_gev2 > 0.025:
        review_reasons.append("large_mean_error")
    # endif
    if not np.isfinite(sigma_error_gev2) or sigma_error_gev2 > 0.025:
        review_reasons.append("large_sigma_error")
    # endif
    if peak_significance_proxy < 5.0:
        review_reasons.append("weak_peak")
    # endif
    if background_model in {"linear", "quadratic", "cubic", "quartic"}:
        if (
            np.isfinite(background_minimum_3sigma)
            and background_minimum_3sigma < -background_negativity_tolerance
        ):
            review_reasons.append("negative_background_in_3sigma_region")
        # endif
        if background_minimum < -background_negativity_tolerance:
            review_reasons.append("negative_background_global_diagnostic")
        # endif
    # endif

    bound_names = ["signal_yield", "mean", "sigma"]
    for index, name in enumerate(bound_names):
        if parameter_at_bound(
            float(parameters[index]),
            float(bounds[0][index]),
            float(bounds[1][index]),
        ):
            review_reasons.append(f"{name}_at_bound")
        # endif
    # endfor

    result = {
        "background_model": background_model,
        "success": True,
        "status": (
            "success" if not review_reasons else "success_flagged_for_review"
        ),
        "review_reasons": sorted(set(review_reasons)),
        "signal_yield": signal_yield,
        "signal_yield_error": signal_yield_error,
        "mean_gev2": mean_gev2,
        "mean_error_gev2": mean_error_gev2,
        "sigma_gev2": sigma_gev2,
        "sigma_error_gev2": sigma_error_gev2,
        "chi2": chi2,
        "ndf": int(ndf),
        "chi2_ndf": chi2_ndf,
        "signal_region_definition": "fitted mean +/- 2 fitted sigma",
        "signal_region_chi2": signal_region_chi2,
        "signal_region_npoints": signal_region_npoints,
        "signal_region_ndf_approx": signal_region_ndf_approx,
        "signal_region_chi2_ndf_approx": signal_region_chi2_ndf_approx,
        "signal_region_chi2_per_bin": signal_region_chi2_per_bin,
        "aicc": float(aicc),
        "peak_significance_proxy": peak_significance_proxy,
        "background_minimum_in_fit_range": background_minimum,
        "background_minimum_in_2sigma_region": background_minimum_2sigma,
        "background_minimum_in_3sigma_region": background_minimum_3sigma,
        "background_negativity_tolerance": background_negativity_tolerance,
        "parameters": parameters.tolist(),
        "parameter_errors": errors.tolist(),
        "covariance": covariance.tolist(),
    }
    diagnostics = {
        "fit_x": fit_x,
        "fit_y": fit_y,
        "fit_error": fit_error,
        "fit_total": fit_total,
        "fit_signal": fit_signal,
        "fit_background": fit_background,
        "pulls": pulls,
    }
    return result, diagnostics


def fit_one_job(job: FitJob) -> dict[str, Any]:
    """Fit four deterministic polynomial hypotheses for one kinematic bin."""
    values = np.asarray(job.values, dtype=float)
    values = values[np.isfinite(values)]

    counts, edges = np.histogram(
        values,
        bins=job.histogram_bins,
        range=(job.histogram_min_gev2, job.histogram_max_gev2),
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    histogram_bin_width_gev2 = float(edges[1] - edges[0])
    fit_mask = (
        (centers >= job.fit_min_gev2)
        & (centers <= job.fit_max_gev2)
    )

    base = {
        "period": job.period,
        "stage": job.stage,
        "period_label": PERIOD_LABELS[job.period],
        "stage_label": STAGE_LABELS[job.stage],
        "x_index": job.x_index,
        "t_index": job.t_index,
        "bin_id": bin_identifier(job.x_index, job.t_index),
        "bin_number": combined_bin_number(job.x_index, job.t_index),
        "xB_min": job.x_min,
        "xB_max": job.x_max,
        "minus_tprime_min_gev2": job.minus_tprime_min_gev2,
        "minus_tprime_max_gev2": job.minus_tprime_max_gev2,
        "number_of_events": int(values.size),
        "histogram_entries": int(np.sum(counts)),
        "histogram_bin_width_gev2": histogram_bin_width_gev2,
        "counts": counts.tolist(),
        "edges": edges.tolist(),
        "models": {},
    }

    if values.size < job.minimum_events:
        base["status"] = "insufficient_events"
        base["models"] = {
            model_name: {
                "background_model": model_name,
                "success": False,
                "status": "insufficient_events",
                "review_reasons": ["insufficient_events"],
            }
            for model_name in BACKGROUND_MODELS
        }
        return base
    # endif

    for model_name in BACKGROUND_MODELS:
        fit_result, _ = fit_background_model(
            background_model=model_name,
            centers=centers,
            counts=counts,
            fit_mask=fit_mask,
            histogram_bin_width_gev2=histogram_bin_width_gev2,
            fit_min_gev2=job.fit_min_gev2,
            fit_max_gev2=job.fit_max_gev2,
            stage=job.stage,
            corrected_mean_max_gev2=job.corrected_mean_max_gev2,
        )
        base["models"][model_name] = fit_result
    # endfor

    polynomial_candidates = ("quadratic", "cubic", "quartic")
    hard_rejection_reasons = {
        "nonfinite_parameters",
        "nonfinite_covariance",
        "signal_yield_at_bound",
        "mean_at_bound",
        "sigma_at_bound",
        "negative_background_in_3sigma_region",
    }

    viable_candidates: list[str] = []
    passing_candidates: list[str] = []

    for model_name in polynomial_candidates:
        candidate = base["models"][model_name]
        reasons = set(candidate.get("review_reasons", []))
        local_score = candidate.get(
            "signal_region_chi2_ndf_approx",
            math.nan,
        )
        is_viable = (
            candidate.get("success", False)
            and np.isfinite(local_score)
            and not bool(reasons & hard_rejection_reasons)
        )
        if is_viable:
            viable_candidates.append(model_name)
            if local_score <= 2.0:
                passing_candidates.append(model_name)
            # endif
        # endif
    # endfor

    if passing_candidates:
        recommended_model = passing_candidates[0]
        recommendation_reason = (
            "lowest viable polynomial order with signal-region "
            "chi2/ndf_approx <= 2.0"
        )
    elif viable_candidates:
        recommended_model = min(
            viable_candidates,
            key=lambda model_name: (
                base["models"][model_name]["signal_region_chi2_ndf_approx"],
                base["models"][model_name]["aicc"],
            ),
        )
        recommendation_reason = (
            "no viable polynomial reached signal-region chi2/ndf_approx <= 2.0; "
            "selected the viable polynomial with the smallest local score"
        )
    else:
        fallback_candidates = [
            model_name
            for model_name in polynomial_candidates
            if (
                base["models"][model_name].get("success", False)
                and np.isfinite(
                    base["models"][model_name].get(
                        "signal_region_chi2_ndf_approx",
                        math.nan,
                    )
                )
            )
        ]
        if fallback_candidates:
            recommended_model = min(
                fallback_candidates,
                key=lambda model_name: (
                    base["models"][model_name]["signal_region_chi2_ndf_approx"],
                    base["models"][model_name]["aicc"],
                ),
            )
            recommendation_reason = (
                "no polynomial passed the hard viability criteria; selected "
                "the successful polynomial with the smallest local score"
            )
        else:
            recommended_model = None
            recommendation_reason = (
                "no successful polynomial fit with a finite signal-region score"
            )
        # endif
    # endif

    base["recommended_background_model"] = recommended_model
    base["recommended_background_order"] = {
        "quadratic": 2,
        "cubic": 3,
        "quartic": 4,
    }.get(recommended_model)
    base["recommendation_reason"] = recommendation_reason

    if recommended_model is None:
        base["status"] = "fit_failed"
        base["fit_accepted"] = False
        base["fit_acceptance_reason"] = "no recommended polynomial fit"
        return base
    # endif

    nominal = base["models"][recommended_model]
    recommended_reasons = set(nominal.get("review_reasons", []))
    recommended_local_score = nominal.get(
        "signal_region_chi2_ndf_approx",
        math.nan,
    )
    fit_accepted = (
        nominal.get("success", False)
        and np.isfinite(recommended_local_score)
        and recommended_local_score <= 3.0
        and not bool(recommended_reasons & hard_rejection_reasons)
    )
    base["fit_accepted"] = bool(fit_accepted)
    if fit_accepted:
        base["status"] = "accepted"
        base["fit_acceptance_reason"] = (
            "recommended fit passes the hard viability criteria and has "
            "signal-region chi2/ndf_approx <= 3.0"
        )
    else:
        base["status"] = "recommended_but_unresolved"
        base["fit_acceptance_reason"] = (
            "recommended fit is the best available extraction but does not "
            "pass the final acceptance criteria"
        )
    # endif

    if nominal.get("success", False):
        base["recommended_selection_windows_gev2"] = {
            "1.5_sigma": {
                "minimum": nominal["mean_gev2"] - 1.5 * nominal["sigma_gev2"],
                "maximum": nominal["mean_gev2"] + 1.5 * nominal["sigma_gev2"],
            },
            "2.0_sigma": {
                "minimum": nominal["mean_gev2"] - 2.0 * nominal["sigma_gev2"],
                "maximum": nominal["mean_gev2"] + 2.0 * nominal["sigma_gev2"],
            },
            "3.0_sigma": {
                "minimum": nominal["mean_gev2"] - 3.0 * nominal["sigma_gev2"],
                "maximum": nominal["mean_gev2"] + 3.0 * nominal["sigma_gev2"],
            },
        }

        model_differences: dict[str, Any] = {}
        for alternative in BACKGROUND_MODELS:
            if alternative == recommended_model:
                continue
            # endif
            alt = base["models"][alternative]
            if alt.get("success", False):
                model_differences[alternative] = {
                    "delta_mean_gev2": alt["mean_gev2"] - nominal["mean_gev2"],
                    "delta_sigma_gev2": alt["sigma_gev2"] - nominal["sigma_gev2"],
                    "delta_signal_yield": (
                        alt["signal_yield"] - nominal["signal_yield"]
                    ),
                    "relative_signal_yield_difference": (
                        (alt["signal_yield"] - nominal["signal_yield"])
                        / nominal["signal_yield"]
                        if nominal["signal_yield"] != 0.0
                        else math.nan
                    ),
                    "delta_aicc": alt["aicc"] - nominal["aicc"],
                }
            else:
                model_differences[alternative] = {
                    "status": alt.get("status", "fit_failed")
                }
            # endif
        # endfor
        base["model_differences_from_recommended"] = model_differences
    # endif

    return base


# =============================================================================
# Job construction and execution
# =============================================================================

def build_fit_jobs(
    datasets: list[InputDataset],
    histogram_min_gev2: float,
    histogram_max_gev2: float,
    histogram_bins: int,
    fit_min_gev2: float,
    fit_max_gev2: float,
    minimum_events: int,
    corrected_mean_max_gev2: float,
    step_size: str,
) -> tuple[list[FitJob], dict[str, Any]]:
    """Load each dataset and construct all 144 deterministic fit jobs."""
    jobs: list[FitJob] = []
    loading_metadata: dict[str, Any] = {}

    for dataset in datasets:
        print(
            f"Reading {PERIOD_LABELS[dataset.period]} {dataset.stage}: "
            f"{dataset.file_path}",
            flush=True,
        )
        x_b, minus_tprime, mx2, counts = load_dataset_arrays(
            dataset=dataset,
            step_size=step_size,
        )

        metadata_key = f"{dataset.period}:{dataset.stage}"
        loading_metadata[metadata_key] = {
            "dataset": asdict(dataset),
            "counts": counts,
        }

        for x_index, (x_min, x_max) in enumerate(XB_BINS):
            x_mask = (x_b >= x_min) & (x_b < x_max)
            for t_index, (t_min, t_max) in enumerate(MINUS_TPRIME_BINS_GEV2):
                mask = (
                    x_mask
                    & (minus_tprime >= t_min)
                    & (minus_tprime < t_max)
                )
                jobs.append(
                    FitJob(
                        period=dataset.period,
                        stage=dataset.stage,
                        x_index=x_index,
                        t_index=t_index,
                        x_min=x_min,
                        x_max=x_max,
                        minus_tprime_min_gev2=t_min,
                        minus_tprime_max_gev2=t_max,
                        values=np.asarray(mx2[mask], dtype=float),
                        histogram_min_gev2=histogram_min_gev2,
                        histogram_max_gev2=histogram_max_gev2,
                        histogram_bins=histogram_bins,
                        fit_min_gev2=fit_min_gev2,
                        fit_max_gev2=fit_max_gev2,
                        minimum_events=minimum_events,
                        corrected_mean_max_gev2=corrected_mean_max_gev2,
                    )
                )
            # endfor
        # endfor
    # endfor

    return jobs, loading_metadata


def execute_fit_jobs(jobs: list[FitJob], workers: int) -> list[dict[str, Any]]:
    """Run all deterministic fits with at most seven worker processes."""
    actual_workers = max(1, min(int(workers), 7, os.cpu_count() or 1))
    print("Running channel_selection_mx2_fits_v5.py", flush=True)
    print(
        f"Running {len(jobs)} kinematic-bin jobs with {actual_workers} worker(s).",
        flush=True,
    )

    if actual_workers == 1:
        results = []
        for index, job in enumerate(jobs, start=1):
            results.append(fit_one_job(job))
            print(f"Completed fit job {index}/{len(jobs)}", flush=True)
        # endfor
        return results
    # endif

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=actual_workers) as executor:
        future_map = {
            executor.submit(fit_one_job, job): job
            for job in jobs
        }
        completed = 0
        for future in as_completed(future_map):
            job = future_map[future]
            try:
                result = future.result()
            except Exception as exc:
                raise RuntimeError(
                    f"Fit worker failed for {job.period} {job.stage} "
                    f"x index {job.x_index}, t index {job.t_index}."
                ) from exc
            # endtry
            results.append(result)
            completed += 1
            print(f"Completed fit job {completed}/{len(jobs)}", flush=True)
        # endfor
    # endwith

    results.sort(
        key=lambda item: (
            item["period"],
            item["stage"],
            item["x_index"],
            item["t_index"],
        )
    )
    return results


# =============================================================================
# Result flattening and output tables
# =============================================================================

def flatten_results(results: list[dict[str, Any]]) -> pd.DataFrame:
    """Create one flat row per period-stage-bin-background fit."""
    rows: list[dict[str, Any]] = []
    for item in results:
        common = {
            "period": item["period"],
            "period_label": item["period_label"],
            "stage": item["stage"],
            "stage_label": item["stage_label"],
            "bin_id": item["bin_id"],
            "bin_number": item["bin_number"],
            "x_index": item["x_index"],
            "t_index": item["t_index"],
            "xB_min": item["xB_min"],
            "xB_max": item["xB_max"],
            "minus_tprime_min_gev2": item["minus_tprime_min_gev2"],
            "minus_tprime_max_gev2": item["minus_tprime_max_gev2"],
            "number_of_events": item["number_of_events"],
            "histogram_entries": item["histogram_entries"],
            "histogram_bin_width_gev2": item["histogram_bin_width_gev2"],
            "recommended_background_model": item.get(
                "recommended_background_model",
                NOMINAL_BACKGROUND_MODEL,
            ),
            "recommended_background_order": item.get(
                "recommended_background_order",
                math.nan,
            ),
            "recommendation_reason": item.get("recommendation_reason", ""),
            "fit_accepted": bool(item.get("fit_accepted", False)),
            "fit_acceptance_reason": item.get(
                "fit_acceptance_reason",
                "",
            ),
        }
        for model_name in BACKGROUND_MODELS:
            fit = item["models"][model_name]
            row = {
                **common,
                "background_model": model_name,
                "background_label": BACKGROUND_LABELS[model_name],
                "is_nominal_quadratic": model_name == NOMINAL_BACKGROUND_MODEL,
                "is_recommended": model_name
                == item.get(
                    "recommended_background_model",
                    NOMINAL_BACKGROUND_MODEL,
                ),
                "success": fit.get("success", False),
                "status": fit.get("status", "unknown"),
                "review_reasons": ";".join(fit.get("review_reasons", [])),
                "signal_yield": fit.get("signal_yield", math.nan),
                "signal_yield_error": fit.get("signal_yield_error", math.nan),
                "mean_gev2": fit.get("mean_gev2", math.nan),
                "mean_error_gev2": fit.get("mean_error_gev2", math.nan),
                "sigma_gev2": fit.get("sigma_gev2", math.nan),
                "sigma_error_gev2": fit.get("sigma_error_gev2", math.nan),
                "chi2": fit.get("chi2", math.nan),
                "ndf": fit.get("ndf", math.nan),
                "chi2_ndf": fit.get("chi2_ndf", math.nan),
                "signal_region_chi2": fit.get(
                    "signal_region_chi2",
                    math.nan,
                ),
                "signal_region_npoints": fit.get(
                    "signal_region_npoints",
                    math.nan,
                ),
                "signal_region_ndf_approx": fit.get(
                    "signal_region_ndf_approx",
                    math.nan,
                ),
                "signal_region_chi2_ndf_approx": fit.get(
                    "signal_region_chi2_ndf_approx",
                    math.nan,
                ),
                "signal_region_chi2_per_bin": fit.get(
                    "signal_region_chi2_per_bin",
                    math.nan,
                ),
                "aicc": fit.get("aicc", math.nan),
                "peak_significance_proxy": fit.get(
                    "peak_significance_proxy",
                    math.nan,
                ),
                "background_minimum_in_fit_range": fit.get(
                    "background_minimum_in_fit_range",
                    math.nan,
                ),
                "background_minimum_in_2sigma_region": fit.get(
                    "background_minimum_in_2sigma_region",
                    math.nan,
                ),
                "background_minimum_in_3sigma_region": fit.get(
                    "background_minimum_in_3sigma_region",
                    math.nan,
                ),
                "background_negativity_tolerance": fit.get(
                    "background_negativity_tolerance",
                    math.nan,
                ),
            }

            differences = item.get(
                "model_differences_from_recommended",
                {},
            ).get(model_name, {})
            is_recommended = model_name == item.get(
                "recommended_background_model",
                NOMINAL_BACKGROUND_MODEL,
            )
            row["delta_mean_from_recommended_gev2"] = differences.get(
                "delta_mean_gev2",
                0.0 if is_recommended else math.nan,
            )
            row["delta_sigma_from_recommended_gev2"] = differences.get(
                "delta_sigma_gev2",
                0.0 if is_recommended else math.nan,
            )
            row["delta_signal_yield_from_recommended"] = differences.get(
                "delta_signal_yield",
                0.0 if is_recommended else math.nan,
            )
            row["relative_signal_yield_difference_from_recommended"] = differences.get(
                "relative_signal_yield_difference",
                0.0 if is_recommended else math.nan,
            )
            row["delta_aicc_from_recommended"] = differences.get(
                "delta_aicc",
                0.0 if is_recommended else math.nan,
            )

            windows = item.get("recommended_selection_windows_gev2", {})
            for window_name in ("1.5_sigma", "2.0_sigma", "3.0_sigma"):
                row[f"window_{window_name}_min_gev2"] = windows.get(
                    window_name,
                    {},
                ).get("minimum", math.nan)
                row[f"window_{window_name}_max_gev2"] = windows.get(
                    window_name,
                    {},
                ).get("maximum", math.nan)
            # endfor
            rows.append(row)
        # endfor
    # endfor

    frame = pd.DataFrame(rows)
    frame.sort_values(
        ["period", "stage", "bin_number", "background_model"],
        inplace=True,
        ignore_index=True,
    )
    return frame


def write_latex_nominal_tables(frame: pd.DataFrame, output_path: Path) -> None:
    """Write corrected adaptively selected mu/sigma tables."""
    selected = frame[
        (frame["stage"] == "after")
        & frame["is_recommended"]
    ].copy()

    lines: list[str] = []
    lines.append("% Auto-generated by channel_selection_mx2_fits_v5.py")
    lines.append("% Corrected Gaussian fits with adaptively selected polynomial backgrounds.")
    lines.append("")

    for period in ("su22", "fa22", "sp23"):
        period_frame = selected[selected["period"] == period].copy()
        period_frame.sort_values("bin_number", inplace=True)
        lines.extend(
            [
                "\\begin{table}[htbp]",
                "\\centering",
                "\\small",
                "\\caption{Missing-neutron peak fit results for "
                f"{PERIOD_LABELS[period]} after momentum corrections. "
                "Each bin uses the adaptively recommended quadratic, cubic, "
                "or quartic background.}",
                f"\\label{{tab:mx2_peak_{period}}}",
                "\\begin{tabular}{cccccccc}",
                "\\hline",
                "$x_{B,\\min}$ & $x_{B,\\max}$ & $(-t')_{\\min}$ & "
                "$(-t')_{\\max}$ & $\\mu$ & $\\delta\\mu$ & "
                "$\\sigma$ & $\\delta\\sigma$ \\\\",
                " &  & (GeV$^2$) & (GeV$^2$) & (GeV$^2$) & "
                "(GeV$^2$) & (GeV$^2$) & (GeV$^2$) \\\\",
                "\\hline",
            ]
        )

        for _, row in period_frame.iterrows():
            if bool(row["success"]):
                values = (
                    f"{row['xB_min']:.4f} & {row['xB_max']:.4f} & "
                    f"{row['minus_tprime_min_gev2']:.4f} & "
                    f"{row['minus_tprime_max_gev2']:.4f} & "
                    f"{row['mean_gev2']:.4f} & {row['mean_error_gev2']:.4f} & "
                    f"{row['sigma_gev2']:.4f} & {row['sigma_error_gev2']:.4f} "
                    "\\\\"
                )
            else:
                values = (
                    f"{row['xB_min']:.4f} & {row['xB_max']:.4f} & "
                    f"{row['minus_tprime_min_gev2']:.4f} & "
                    f"{row['minus_tprime_max_gev2']:.4f} & "
                    "-- & -- & -- & -- \\\\"
                )
            # endif
            lines.append(values)
        # endfor

        lines.extend(
            [
                "\\hline",
                "\\end{tabular}",
                "\\end{table}",
                "",
            ]
        )
    # endfor

    output_path.write_text("\n".join(lines), encoding="utf-8")


# =============================================================================
# Plotting helpers
# =============================================================================

def result_lookup(
    results: list[dict[str, Any]],
) -> dict[tuple[str, str, int, int], dict[str, Any]]:
    """Index result dictionaries by period, stage, x index, and t index."""
    return {
        (
            item["period"],
            item["stage"],
            item["x_index"],
            item["t_index"],
        ): item
        for item in results
    }


def evaluate_model_dense(
    item: dict[str, Any],
    model_name: str,
    x_dense: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Evaluate total, signal, and background curves for plotting."""
    fit = item["models"][model_name]
    if not fit.get("success", False):
        return None
    # endif

    parameters = np.asarray(fit["parameters"], dtype=float)
    bin_width = float(item["histogram_bin_width_gev2"])
    model_function = build_total_model(model_name, bin_width)
    total = model_function(x_dense, *parameters)
    signal = gaussian_signal_from_yield(
        x_dense,
        parameters[0],
        parameters[1],
        parameters[2],
        bin_width,
    )
    background = evaluate_background(model_name, x_dense, parameters)
    return total, signal, background


def plot_spectrum_canvases(
    results: list[dict[str, Any]],
    output_dir: Path,
    fit_min_gev2: float,
    fit_max_gev2: float,
) -> None:
    """Make one 4 x 6 adaptively selected polynomial-fit canvas per dataset."""
    lookup = result_lookup(results)
    ensure_directory(output_dir)

    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            fig, axes = plt.subplots(4, 6, figsize=(20, 13), squeeze=False)
            for x_index, (x_min, x_max) in enumerate(XB_BINS):
                for t_index, (t_min, t_max) in enumerate(MINUS_TPRIME_BINS_GEV2):
                    ax = axes[x_index][t_index]
                    item = lookup[(period, stage, x_index, t_index)]
                    counts = np.asarray(item["counts"], dtype=float)
                    edges = np.asarray(item["edges"], dtype=float)
                    centers = 0.5 * (edges[:-1] + edges[1:])
                    errors = np.sqrt(np.maximum(counts, 1.0))

                    ax.errorbar(
                        centers,
                        counts,
                        yerr=errors,
                        fmt=".",
                        markersize=2.6,
                        linewidth=0.6,
                        capsize=0.0,
                        label="Data",
                    )
                    ax.axvline(
                        fit_min_gev2,
                        linewidth=0.7,
                        linestyle=":",
                        alpha=0.55,
                    )
                    ax.axvline(
                        fit_max_gev2,
                        linewidth=0.7,
                        linestyle=":",
                        alpha=0.55,
                    )

                    x_dense = np.linspace(fit_min_gev2, fit_max_gev2, 500)
                    selected_model = item.get("recommended_background_model")
                    if selected_model is None:
                        selected_model = NOMINAL_BACKGROUND_MODEL
                    # endif
                    evaluated = evaluate_model_dense(
                        item,
                        selected_model,
                        x_dense,
                    )
                    nominal = item["models"][selected_model]
                    if evaluated is not None:
                        total, signal, background = evaluated
                        ax.plot(x_dense, total, linewidth=1.2, label="Total fit")
                        ax.plot(
                            x_dense,
                            signal,
                            linewidth=1.0,
                            linestyle="--",
                            label="Gaussian",
                        )
                        ax.plot(
                            x_dense,
                            background,
                            linewidth=1.0,
                            linestyle=":",
                            label="Background",
                        )

                        mean = nominal["mean_gev2"]
                        sigma = nominal["sigma_gev2"]
                        ax.axvline(mean - 2.0 * sigma, linewidth=0.8, linestyle="--")
                        ax.axvline(mean + 2.0 * sigma, linewidth=0.8, linestyle="--")
                        annotation = (
                            f"$\\mu$={mean:.4f}\n"
                            f"$\\sigma$={sigma:.4f}\n"
                            f"model={selected_model}\n"
                            f"$\\chi^2$/ndf$_{{\\pm2\\sigma,approx}}$="
                            f"{nominal['signal_region_chi2_ndf_approx']:.2f}\n"
                            f"{'ACCEPTED' if item.get('fit_accepted', False) else 'UNRESOLVED'}"
                        )
                    else:
                        annotation = nominal.get("status", "fit failed")
                    # endif

                    ax.text(
                        0.03,
                        0.96,
                        annotation,
                        transform=ax.transAxes,
                        ha="left",
                        va="top",
                        fontsize=7,
                    )
                    ax.set_title(
                        f"{x_min:.2f} < $x_B$ < {x_max:.2f}\n"
                        f"{t_min:.2f} < $-t'$ < {t_max:.2f}",
                        fontsize=8,
                    )
                    ax.grid(True, alpha=0.35)
                    ax.tick_params(labelsize=7)
                    if x_index == 3:
                        ax.set_xlabel("$M_x^2$ (GeV$^2$)", fontsize=8)
                    # endif
                    if t_index == 0:
                        ax.set_ylabel("Counts", fontsize=8)
                    # endif
                # endfor
            # endfor

            handles, labels = axes[0][0].get_legend_handles_labels()
            if handles:
                fig.legend(
                    handles,
                    labels,
                    loc="upper center",
                    bbox_to_anchor=(0.5, 0.955),
                    ncol=4,
                    fontsize=9,
                    frameon=True,
                )
            # endif
            fig.suptitle(
                f"{PERIOD_LABELS[period]}: {STAGE_LABELS[stage]}\n"
                "Gaussian signal plus adaptively selected polynomial background",
                fontsize=15,
                y=0.995,
            )
            fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.900))
            output_path = output_dir / f"mx2_fits_{period}_{stage}.png"
            fig.savefig(output_path, dpi=180)
            plt.close(fig)
        # endfor
    # endfor


def plot_background_model_canvases(
    results: list[dict[str, Any]],
    output_dir: Path,
    fit_min_gev2: float,
    fit_max_gev2: float,
) -> None:
    """Overlay all three fitted total models for corrected data."""
    lookup = result_lookup(results)
    ensure_directory(output_dir)

    for period in ("su22", "fa22", "sp23"):
        fig, axes = plt.subplots(4, 6, figsize=(20, 13), squeeze=False)
        for x_index, (x_min, x_max) in enumerate(XB_BINS):
            for t_index, (t_min, t_max) in enumerate(MINUS_TPRIME_BINS_GEV2):
                ax = axes[x_index][t_index]
                item = lookup[(period, "after", x_index, t_index)]
                counts = np.asarray(item["counts"], dtype=float)
                edges = np.asarray(item["edges"], dtype=float)
                centers = 0.5 * (edges[:-1] + edges[1:])
                errors = np.sqrt(np.maximum(counts, 1.0))
                ax.errorbar(
                    centers,
                    counts,
                    yerr=errors,
                    fmt=".",
                    markersize=2.4,
                    linewidth=0.6,
                    capsize=0.0,
                    label="Data",
                )
                ax.axvline(
                    fit_min_gev2,
                    linewidth=0.7,
                    linestyle=":",
                    alpha=0.55,
                )
                ax.axvline(
                    fit_max_gev2,
                    linewidth=0.7,
                    linestyle=":",
                    alpha=0.55,
                )

                x_dense = np.linspace(fit_min_gev2, fit_max_gev2, 500)
                for model_name in BACKGROUND_MODELS:
                    evaluated = evaluate_model_dense(item, model_name, x_dense)
                    if evaluated is None:
                        continue
                    # endif
                    total, _, _ = evaluated
                    ax.plot(
                        x_dense,
                        total,
                        linewidth=1.1,
                        label=BACKGROUND_LABELS[model_name],
                    )
                # endfor

                ax.set_title(
                    f"{x_min:.2f} < $x_B$ < {x_max:.2f}\n"
                    f"{t_min:.2f} < $-t'$ < {t_max:.2f}",
                    fontsize=8,
                )
                ax.grid(True, alpha=0.35)
                ax.tick_params(labelsize=7)
                if x_index == 3:
                    ax.set_xlabel("$M_x^2$ (GeV$^2$)", fontsize=8)
                # endif
                if t_index == 0:
                    ax.set_ylabel("Counts", fontsize=8)
                # endif
            # endfor
        # endfor

        handles, labels = axes[0][0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.955),
                ncol=4,
                fontsize=9,
                frameon=True,
            )
        # endif
        fig.suptitle(
            f"{PERIOD_LABELS[period]} after momentum corrections: "
            "background-model comparison",
            fontsize=15,
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.900))
        fig.savefig(
            output_dir / f"mx2_background_models_{period}_after.png",
            dpi=180,
        )
        plt.close(fig)
    # endfor



def padded_limits(values: np.ndarray, fractional_padding: float = 0.08) -> tuple[float, float]:
    """Return finite padded limits for a family of matching summary plots."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (0.0, 1.0)
    # endif
    low = float(np.min(finite))
    high = float(np.max(finite))
    span = high - low
    if span <= 0.0:
        span = max(abs(low), 1.0) * 0.10
    # endif
    padding = fractional_padding * span
    return low - padding, high + padding


def symmetric_limits(values: np.ndarray, fractional_padding: float = 0.10) -> tuple[float, float]:
    """Return symmetric limits about zero for matching variation plots."""
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return (-1.0, 1.0)
    # endif
    maximum = float(np.max(np.abs(finite)))
    if maximum <= 0.0:
        maximum = 1.0
    # endif
    maximum *= 1.0 + fractional_padding
    return -maximum, maximum

def plot_before_after_summary(frame: pd.DataFrame, output_dir: Path) -> None:
    """Plot recommended mu and sigma before versus after momentum corrections."""
    ensure_directory(output_dir)
    recommended = frame[frame["is_recommended"] & frame["success"]].copy()

    mu_low, mu_high = padded_limits(
        np.concatenate(
            [
                recommended["mean_gev2"].to_numpy(dtype=float),
                np.asarray([NEUTRON_MASS2_GEV2], dtype=float),
            ]
        ),
        fractional_padding=0.08,
    )
    sigma_low, sigma_high = padded_limits(
        recommended["sigma_gev2"].to_numpy(dtype=float),
        fractional_padding=0.10,
    )
    sigma_low = max(0.0, sigma_low)

    for period in ("su22", "fa22", "sp23"):
        period_frame = recommended[recommended["period"] == period]

        for quantity, error, ylabel, title_word, limits, filename_prefix in (
            (
                "mean_gev2",
                "mean_error_gev2",
                "Gaussian mean $\\mu$ (GeV$^2$)",
                "mean",
                (mu_low, mu_high),
                "mu",
            ),
            (
                "sigma_gev2",
                "sigma_error_gev2",
                "Gaussian width $\\sigma$ (GeV$^2$)",
                "width",
                (sigma_low, sigma_high),
                "sigma",
            ),
        ):
            fig, ax = plt.subplots(figsize=(12, 5.5))
            for stage in ("before", "after"):
                selected = period_frame[
                    period_frame["stage"] == stage
                ].sort_values("bin_number")
                accepted = selected[selected["fit_accepted"]]
                unresolved = selected[~selected["fit_accepted"]]

                ax.errorbar(
                    accepted["bin_number"],
                    accepted[quantity],
                    yerr=accepted[error],
                    fmt="o",
                    markersize=4,
                    capsize=2,
                    label=f"{STAGE_LABELS[stage]} accepted",
                )
                ax.errorbar(
                    unresolved["bin_number"],
                    unresolved[quantity],
                    yerr=unresolved[error],
                    fmt="o",
                    markerfacecolor="none",
                    markersize=5,
                    capsize=2,
                    label=f"{STAGE_LABELS[stage]} unresolved",
                )
            # endfor

            if quantity == "mean_gev2":
                ax.axhline(
                    NEUTRON_MASS2_GEV2,
                    linewidth=1.0,
                    linestyle="--",
                    label="$m_n^2$",
                )
            # endif

            ax.set_xlabel("Kinematic bin")
            ax.set_ylabel(ylabel)
            ax.set_title(
                f"{PERIOD_LABELS[period]} missing-neutron peak {title_word}"
            )
            ax.set_xlim(0.5, 24.5)
            ax.set_ylim(*limits)
            ax.set_xticks(np.arange(1, 25))
            ax.grid(True, alpha=0.4)
            ax.legend(ncol=2)
            fig.tight_layout()
            fig.savefig(
                output_dir / f"{filename_prefix}_before_after_{period}.png",
                dpi=180,
            )
            plt.close(fig)
        # endfor
    # endfor


def plot_corrected_period_summary(frame: pd.DataFrame, output_dir: Path) -> None:
    """Compare corrected recommended mu and sigma among the run periods."""
    ensure_directory(output_dir)
    selected = frame[
        (frame["stage"] == "after")
        & frame["is_recommended"]
        & frame["success"]
    ].copy()

    mu_low, mu_high = padded_limits(
        np.concatenate(
            [
                selected["mean_gev2"].to_numpy(dtype=float),
                np.asarray([NEUTRON_MASS2_GEV2], dtype=float),
            ]
        ),
        fractional_padding=0.08,
    )
    sigma_low, sigma_high = padded_limits(
        selected["sigma_gev2"].to_numpy(dtype=float),
        fractional_padding=0.10,
    )
    sigma_low = max(0.0, sigma_low)

    for quantity, error, ylabel, title_word, limits, filename_prefix in (
        (
            "mean_gev2",
            "mean_error_gev2",
            "Gaussian mean $\\mu$ (GeV$^2$)",
            "mean",
            (mu_low, mu_high),
            "mu",
        ),
        (
            "sigma_gev2",
            "sigma_error_gev2",
            "Gaussian width $\\sigma$ (GeV$^2$)",
            "width",
            (sigma_low, sigma_high),
            "sigma",
        ),
    ):
        fig, ax = plt.subplots(figsize=(12, 5.5))
        for period in ("su22", "fa22", "sp23"):
            period_frame = selected[
                selected["period"] == period
            ].sort_values("bin_number")
            accepted = period_frame[period_frame["fit_accepted"]]
            unresolved = period_frame[~period_frame["fit_accepted"]]

            ax.errorbar(
                accepted["bin_number"],
                accepted[quantity],
                yerr=accepted[error],
                fmt="o",
                markersize=4,
                capsize=2,
                label=f"{PERIOD_LABELS[period]} accepted",
            )
            ax.errorbar(
                unresolved["bin_number"],
                unresolved[quantity],
                yerr=unresolved[error],
                fmt="o",
                markerfacecolor="none",
                markersize=5,
                capsize=2,
                label=f"{PERIOD_LABELS[period]} unresolved",
            )
        # endfor

        if quantity == "mean_gev2":
            ax.axhline(
                NEUTRON_MASS2_GEV2,
                linewidth=1.0,
                linestyle="--",
                label="$m_n^2$",
            )
        # endif

        ax.set_xlabel("Kinematic bin")
        ax.set_ylabel(ylabel)
        ax.set_title(
            f"Corrected missing-neutron peak {title_word} by run period"
        )
        ax.set_xlim(0.5, 24.5)
        ax.set_ylim(*limits)
        ax.set_xticks(np.arange(1, 25))
        ax.grid(True, alpha=0.4)
        ax.legend(ncol=2)
        fig.tight_layout()
        fig.savefig(
            output_dir / f"{filename_prefix}_corrected_period_comparison.png",
            dpi=180,
        )
        plt.close(fig)
    # endfor


def plot_model_variations(frame: pd.DataFrame, output_dir: Path) -> None:
    """Plot corrected model shifts relative to the adaptively recommended model."""
    ensure_directory(output_dir)
    selected = frame[
        (frame["stage"] == "after")
        & frame["success"]
    ].copy()

    yield_low, yield_high = symmetric_limits(
        100.0
        * selected[
            "relative_signal_yield_difference_from_recommended"
        ].to_numpy(dtype=float),
        fractional_padding=0.10,
    )
    mu_delta_low, mu_delta_high = symmetric_limits(
        selected["delta_mean_from_recommended_gev2"].to_numpy(dtype=float),
        fractional_padding=0.10,
    )
    sigma_delta_low, sigma_delta_high = symmetric_limits(
        selected["delta_sigma_from_recommended_gev2"].to_numpy(dtype=float),
        fractional_padding=0.10,
    )

    for period in ("su22", "fa22", "sp23"):
        period_frame = selected[selected["period"] == period]

        fig, ax = plt.subplots(figsize=(12, 5.5))
        for model_name in BACKGROUND_MODELS:
            model_frame = period_frame[
                period_frame["background_model"] == model_name
            ].sort_values("bin_number")
            ax.plot(
                model_frame["bin_number"],
                100.0
                * model_frame[
                    "relative_signal_yield_difference_from_recommended"
                ],
                marker="o",
                markersize=4,
                linewidth=1.0,
                label=BACKGROUND_LABELS[model_name],
            )
        # endfor
        ax.axhline(0.0, linewidth=1.0, linestyle="--")
        ax.set_xlabel("Kinematic bin")
        ax.set_ylabel("Signal-yield difference from recommended model (%)")
        ax.set_ylim(yield_low, yield_high)
        ax.set_title(
            f"{PERIOD_LABELS[period]} corrected background-model dependence\n"
            "relative to the adaptively recommended model"
        )
        ax.set_xlim(0.5, 24.5)
        ax.set_xticks(np.arange(1, 25))
        ax.grid(True, alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            output_dir / f"signal_yield_model_variation_{period}.png",
            dpi=180,
        )
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(12, 5.5))
        for model_name in BACKGROUND_MODELS:
            model_frame = period_frame[
                period_frame["background_model"] == model_name
            ].sort_values("bin_number")
            ax.plot(
                model_frame["bin_number"],
                model_frame["delta_mean_from_recommended_gev2"],
                marker="o",
                markersize=4,
                linewidth=1.0,
                label=BACKGROUND_LABELS[model_name],
            )
        # endfor
        ax.axhline(0.0, linewidth=1.0, linestyle="--")
        ax.set_xlabel("Kinematic bin")
        ax.set_ylabel("$\\mu_{model}-\\mu_{quadratic}$ (GeV$^2$)")
        ax.set_ylim(mu_delta_low, mu_delta_high)
        ax.set_title(
            f"{PERIOD_LABELS[period]} corrected mean dependence on background"
        )
        ax.set_xlim(0.5, 24.5)
        ax.set_xticks(np.arange(1, 25))
        ax.grid(True, alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            output_dir / f"mu_model_variation_{period}.png",
            dpi=180,
        )
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(12, 5.5))
        for model_name in BACKGROUND_MODELS:
            model_frame = period_frame[
                period_frame["background_model"] == model_name
            ].sort_values("bin_number")
            ax.plot(
                model_frame["bin_number"],
                model_frame["delta_sigma_from_recommended_gev2"],
                marker="o",
                markersize=4,
                linewidth=1.0,
                label=BACKGROUND_LABELS[model_name],
            )
        # endfor
        ax.axhline(0.0, linewidth=1.0, linestyle="--")
        ax.set_xlabel("Kinematic bin")
        ax.set_ylabel("$\\sigma_{model}-\\sigma_{quadratic}$ (GeV$^2$)")
        ax.set_ylim(sigma_delta_low, sigma_delta_high)
        ax.set_title(
            f"{PERIOD_LABELS[period]} corrected width dependence on background"
        )
        ax.set_xlim(0.5, 24.5)
        ax.set_xticks(np.arange(1, 25))
        ax.grid(True, alpha=0.4)
        ax.legend()
        fig.tight_layout()
        fig.savefig(
            output_dir / f"sigma_model_variation_{period}.png",
            dpi=180,
        )
        plt.close(fig)
    # endfor



def plot_problem_bin_pull_diagnostics(
    results: list[dict[str, Any]],
    output_dir: Path,
    fit_min_gev2: float,
    fit_max_gev2: float,
) -> None:
    """Plot corrected spectra and pulls for difficult bins 1, 7, 13, and 19."""
    ensure_directory(output_dir)
    lookup = result_lookup(results)
    problem_bins = (1, 7, 13, 19)

    for period in ("su22", "fa22", "sp23"):
        fig, axes = plt.subplots(
            len(problem_bins),
            2,
            figsize=(13, 13),
            squeeze=False,
        )

        for row_index, bin_number in enumerate(problem_bins):
            x_index = (bin_number - 1) // len(MINUS_TPRIME_BINS_GEV2)
            t_index = (bin_number - 1) % len(MINUS_TPRIME_BINS_GEV2)
            item = lookup[(period, "after", x_index, t_index)]
            model_name = item.get("recommended_background_model")
            spectrum_ax, pull_ax = axes[row_index]

            if model_name is None:
                spectrum_ax.text(
                    0.5,
                    0.5,
                    "No recommended fit",
                    ha="center",
                    va="center",
                    transform=spectrum_ax.transAxes,
                )
                pull_ax.axis("off")
                continue
            # endif

            fit = item["models"][model_name]
            counts = np.asarray(item["counts"], dtype=float)
            edges = np.asarray(item["edges"], dtype=float)
            centers = 0.5 * (edges[:-1] + edges[1:])
            errors = np.sqrt(np.maximum(counts, 1.0))

            spectrum_ax.errorbar(
                centers,
                counts,
                yerr=errors,
                fmt=".",
                markersize=3,
                linewidth=0.6,
                label="Data",
            )

            x_dense = np.linspace(fit_min_gev2, fit_max_gev2, 600)
            evaluated = evaluate_model_dense(item, model_name, x_dense)
            if evaluated is not None:
                total, signal, background = evaluated
                spectrum_ax.plot(x_dense, total, linewidth=1.2, label="Total fit")
                spectrum_ax.plot(
                    x_dense,
                    signal,
                    linewidth=1.0,
                    linestyle="--",
                    label="Gaussian",
                )
                spectrum_ax.plot(
                    x_dense,
                    background,
                    linewidth=1.0,
                    linestyle=":",
                    label="Background",
                )
            # endif

            mean = float(fit["mean_gev2"])
            sigma = float(fit["sigma_gev2"])
            spectrum_ax.axvline(mean - 2.0 * sigma, linestyle="--", linewidth=0.8)
            spectrum_ax.axvline(mean + 2.0 * sigma, linestyle="--", linewidth=0.8)

            fit_mask = (
                (centers >= fit_min_gev2)
                & (centers <= fit_max_gev2)
            )
            fit_centers = centers[fit_mask]
            fit_counts = counts[fit_mask]
            fit_errors = errors[fit_mask]
            model = build_total_model(
                model_name,
                float(item["histogram_bin_width_gev2"]),
            )
            parameters = np.asarray(fit["parameters"], dtype=float)
            pulls = (
                fit_counts - model(fit_centers, *parameters)
            ) / fit_errors

            pull_ax.axhline(0.0, linewidth=0.9)
            pull_ax.axhline(2.0, linewidth=0.7, linestyle=":")
            pull_ax.axhline(-2.0, linewidth=0.7, linestyle=":")
            pull_ax.plot(
                fit_centers,
                pulls,
                marker="o",
                markersize=3,
                linewidth=0.8,
            )
            pull_ax.axvspan(
                mean - 2.0 * sigma,
                mean + 2.0 * sigma,
                alpha=0.10,
            )
            pull_ax.set_ylim(-8.0, 8.0)

            status = "accepted" if item.get("fit_accepted", False) else "unresolved"
            spectrum_ax.set_title(
                f"Bin {bin_number}: {model_name}, {status}\\n"
                f"$\\\\chi^2$/ndf$_{{\\\\pm2\\\\sigma,approx}}$="
                f"{fit['signal_region_chi2_ndf_approx']:.2f}",
                fontsize=9,
            )
            pull_ax.set_title("Pulls over full fit range", fontsize=9)
            spectrum_ax.set_ylabel("Counts")
            pull_ax.set_ylabel("(Data - fit) / uncertainty")
            spectrum_ax.grid(True, alpha=0.35)
            pull_ax.grid(True, alpha=0.35)

            if row_index == len(problem_bins) - 1:
                spectrum_ax.set_xlabel("$M_x^2$ (GeV$^2$)")
                pull_ax.set_xlabel("$M_x^2$ (GeV$^2$)")
            # endif
        # endfor

        handles, labels = axes[0][0].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                labels,
                loc="upper center",
                ncol=4,
                fontsize=9,
            )
        # endif
        fig.suptitle(
            f"{PERIOD_LABELS[period]} corrected difficult-bin diagnostics",
            y=0.995,
        )
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.965))
        fig.savefig(
            output_dir / f"problem_bin_pulls_{period}_after.png",
            dpi=180,
        )
        plt.close(fig)
    # endfor


def write_fit_status_summary(frame: pd.DataFrame, output_path: Path) -> None:
    """Write a concise text report of fit success and flagged bins."""
    lines: list[str] = []
    lines.append("Missing-neutron Mx2 deterministic fit summary")
    lines.append("=" * 52)
    lines.append("")

    grouped = frame.groupby(["period", "stage", "background_model"], sort=True)
    for (period, stage, model_name), group in grouped:
        successful = int(group["success"].sum())
        flagged = int(
            np.count_nonzero(group["status"] == "success_flagged_for_review")
        )
        failed = int(len(group) - successful)
        lines.append(
            f"{PERIOD_LABELS[period]:4s} | {stage:6s} | "
            f"{model_name:11s} : successful={successful:2d}/24, "
            f"flagged={flagged:2d}, failed={failed:2d}"
        )
    # endfor

    lines.append("")
    recommended = frame[frame["is_recommended"]].copy()
    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            group = recommended[
                (recommended["period"] == period)
                & (recommended["stage"] == stage)
            ]
            accepted = int(group["fit_accepted"].sum())
            unresolved = int(len(group) - accepted)
            lines.append(
                f"{PERIOD_LABELS[period]:4s} | {stage:6s} | recommended fits: "
                f"accepted={accepted:2d}/24, unresolved={unresolved:2d}"
            )
        # endfor
    # endfor

    lines.append("")
    lines.append("Flagged or failed fits")
    lines.append("-" * 52)
    problem_rows = frame[
        (~frame["success"])
        | (frame["status"] == "success_flagged_for_review")
    ]
    if problem_rows.empty:
        lines.append("None")
    else:
        for _, row in problem_rows.iterrows():
            lines.append(
                f"{row['period']} {row['stage']} bin {int(row['bin_number']):02d} "
                f"{row['background_model']}: {row['status']} "
                f"[{row['review_reasons']}]"
            )
        # endfor
    # endif

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# =============================================================================
# JSON assembly
# =============================================================================

def build_compact_json(
    results: list[dict[str, Any]],
    datasets: list[InputDataset],
    loading_metadata: dict[str, Any],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """Build the production-oriented deterministic fit JSON."""
    payload: dict[str, Any] = {
        "schema_version": 1,
        "analysis": "RGC exclusive enpi+ missing-neutron peak fits",
        "script": Path(__file__).name,
        "tree_name_requested": args.tree,
        "physics": {
            "neutron_mass_gev": NEUTRON_MASS_GEV,
            "neutron_mass_squared_gev2": NEUTRON_MASS2_GEV2,
            "minus_tprime_definition": "-tprime",
        },
        "binning": {
            "xB": [list(bounds) for bounds in XB_BINS],
            "minus_tprime_gev2": [
                list(bounds) for bounds in MINUS_TPRIME_BINS_GEV2
            ],
        },
        "histogram": {
            "minimum_gev2": args.hist_min,
            "maximum_gev2": args.hist_max,
            "number_of_bins": args.hist_bins,
        },
        "fit": {
            "minimum_gev2": args.fit_min,
            "maximum_gev2": args.fit_max,
            "minimum_events": args.minimum_events,
            "signal_model": "Gaussian parameterized by integrated yield",
            "reference_background_model": NOMINAL_BACKGROUND_MODEL,
            "production_fit_per_bin": "recommended_background_model",
            "background_models": list(BACKGROUND_MODELS),
            "adaptive_polynomial_selection": {
                "candidate_orders": [2, 3, 4],
                "preferred_threshold_signal_region_chi2_ndf_approx": 2.0,
                "final_acceptance_threshold_signal_region_chi2_ndf_approx": 3.0,
                "background_positivity_region": "fitted mean +/- 3 fitted sigma",
                "background_negativity_tolerance": (
                    "max(1 count, 1 percent of the characteristic count scale)"
                ),
                "policy": (
                    "Choose the lowest viable polynomial order satisfying the "
                    "preferred local threshold. If none passes, choose the viable "
                    "order with the smallest local score. Recommendation and final "
                    "acceptance are stored separately."
                ),
            },
            "signal_region_fit_quality": {
                "definition": "fitted mean +/- 2 fitted sigma",
                "ndf_convention": (
                    "number of histogram bins in the signal region minus the "
                    "three Gaussian parameters; this is a local diagnostic, not "
                    "a strict independent goodness-of-fit probability"
                ),
                "global_chi2_ndf_also_stored": True,
            },
            "replicas_performed": False,
        },
        "execution": {
            "requested_workers": args.workers,
            "maximum_allowed_workers": 7,
            "actual_workers": max(
                1,
                min(args.workers, 7, os.cpu_count() or 1),
            ),
        },
        "datasets": [asdict(dataset) for dataset in datasets],
        "loading_metadata": loading_metadata,
        "periods": {},
    }

    lookup = result_lookup(results)
    for period in ("su22", "fa22", "sp23"):
        payload["periods"][period] = {
            "label": PERIOD_LABELS[period],
            "before": {"bins": {}},
            "after": {"bins": {}},
        }
        for stage in ("before", "after"):
            for x_index, _ in enumerate(XB_BINS):
                for t_index, _ in enumerate(MINUS_TPRIME_BINS_GEV2):
                    item = lookup[(period, stage, x_index, t_index)]
                    compact_item = {
                        key: value
                        for key, value in item.items()
                        if key not in {"counts", "edges"}
                    }
                    payload["periods"][period][stage]["bins"][
                        item["bin_id"]
                    ] = compact_item
                # endfor
            # endfor
        # endfor
    # endfor

    return payload


# =============================================================================
# Command line and main
# =============================================================================

def build_argument_parser() -> argparse.ArgumentParser:
    """Construct the command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Fit the RGC e pi+ missing-neutron Mx2 peak in 24 kinematic bins "
            "for three run periods before and after momentum corrections."
        )
    )
    parser.add_argument(
        "--tree",
        default="PhysicsEvents",
        help="Requested ROOT tree name (default: PhysicsEvents).",
    )
    parser.add_argument(
        "--input",
        action="append",
        type=parse_dataset_override,
        default=[],
        metavar="PERIOD:STAGE=FILE",
        help=(
            "Override one default input. PERIOD is su22, fa22, or sp23; "
            "STAGE is before or after. May be repeated."
        ),
    )
    parser.add_argument(
        "--branch",
        action="append",
        type=parse_branch_override,
        default=[],
        metavar="LOGICAL=ROOT_BRANCH",
        help=(
            "Override branch resolution for xB, tprime, or Mx2. "
            "May be repeated."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("output/channel_selection_mx2_fit_stability"),
        help=(
            "Stable output directory (default: "
            "output/channel_selection_mx2_fit_stability)."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=7,
        help="Number of worker processes; hard-capped at 7 (default: 7).",
    )
    parser.add_argument(
        "--step-size",
        default="250 MB",
        help="uproot iteration step size (default: 250 MB).",
    )
    parser.add_argument(
        "--hist-min",
        type=float,
        default=0.20,
        help="Mx2 histogram minimum in GeV^2 (default: 0.20).",
    )
    parser.add_argument(
        "--hist-max",
        type=float,
        default=1.55,
        help="Mx2 histogram maximum in GeV^2 (default: 1.55).",
    )
    parser.add_argument(
        "--hist-bins",
        type=int,
        default=60,
        help="Number of Mx2 histogram bins (default: 60).",
    )
    parser.add_argument(
        "--fit-min",
        type=float,
        default=0.50,
        help="Mx2 total-fit minimum in GeV^2 (default: 0.50).",
    )
    parser.add_argument(
        "--fit-max",
        type=float,
        default=1.35,
        help="Mx2 total-fit maximum in GeV^2 (default: 1.35).",
    )
    parser.add_argument(
        "--corrected-mean-max",
        type=float,
        default=0.910,
        help=(
            "Upper bound on the corrected-data Gaussian mean in GeV^2 "
            "(default: 0.910). Before-correction fits retain a wider bound."
        ),
    )
    parser.add_argument(
        "--minimum-events",
        type=int,
        default=150,
        help="Minimum events required in a kinematic bin (default: 150).",
    )
    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate numerical command-line options."""
    if args.workers < 1:
        raise ValueError("--workers must be at least 1.")
    # endif
    if args.hist_bins < 20:
        raise ValueError("--hist-bins must be at least 20.")
    # endif
    if args.hist_min >= args.hist_max:
        raise ValueError("--hist-min must be smaller than --hist-max.")
    # endif
    if args.fit_min >= args.fit_max:
        raise ValueError("--fit-min must be smaller than --fit-max.")
    # endif
    if args.fit_min < args.hist_min or args.fit_max > args.hist_max:
        raise ValueError("The fit range must lie inside the histogram range.")
    # endif
    if args.minimum_events < 1:
        raise ValueError("--minimum-events must be positive.")
    # endif


def main() -> int:
    """Run the deterministic pass-1 fit-stability study."""
    parser = build_argument_parser()
    args = parser.parse_args()

    try:
        validate_arguments(args)

        input_paths = {
            period: {
                stage: Path(path)
                for stage, path in stage_map.items()
            }
            for period, stage_map in DEFAULT_INPUTS.items()
        }
        for period, stage, path in args.input:
            input_paths[period][stage] = path
        # endfor

        branch_overrides = dict(args.branch)
        datasets = preflight_inputs(
            input_paths=input_paths,
            requested_tree=args.tree,
            branch_overrides=branch_overrides,
        )

        output_dir = args.output_dir.expanduser().resolve()
        table_dir = output_dir / "tables"
        plot_dir = output_dir / "plots"
        spectrum_dir = plot_dir / "spectra"
        model_canvas_dir = plot_dir / "background_models"
        summary_dir = plot_dir / "summaries"
        pull_diagnostic_dir = plot_dir / "problem_bin_diagnostics"
        for path in (
            output_dir,
            table_dir,
            plot_dir,
            spectrum_dir,
            model_canvas_dir,
            summary_dir,
            pull_diagnostic_dir,
        ):
            ensure_directory(path)
        # endfor

        jobs, loading_metadata = build_fit_jobs(
            datasets=datasets,
            histogram_min_gev2=args.hist_min,
            histogram_max_gev2=args.hist_max,
            histogram_bins=args.hist_bins,
            fit_min_gev2=args.fit_min,
            fit_max_gev2=args.fit_max,
            minimum_events=args.minimum_events,
            corrected_mean_max_gev2=args.corrected_mean_max,
            step_size=args.step_size,
        )
        results = execute_fit_jobs(jobs=jobs, workers=args.workers)
        frame = flatten_results(results)

        csv_path = table_dir / "mx2_peak_fit_results_v5.csv"
        json_path = table_dir / "mx2_peak_fit_results_v5.json"
        latex_path = table_dir / "mx2_peak_nominal_corrected_tables_v5.tex"
        status_path = table_dir / "mx2_peak_fit_status_v5.txt"

        frame.to_csv(csv_path, index=False)
        compact_json = build_compact_json(
            results=results,
            datasets=datasets,
            loading_metadata=loading_metadata,
            args=args,
        )
        json_path.write_text(
            json.dumps(json_safe(compact_json), indent=2, sort_keys=False) + "\n",
            encoding="utf-8",
        )
        write_latex_nominal_tables(frame, latex_path)
        write_fit_status_summary(frame, status_path)

        plot_spectrum_canvases(
            results=results,
            output_dir=spectrum_dir,
            fit_min_gev2=args.fit_min,
            fit_max_gev2=args.fit_max,
        )
        plot_background_model_canvases(
            results=results,
            output_dir=model_canvas_dir,
            fit_min_gev2=args.fit_min,
            fit_max_gev2=args.fit_max,
        )
        plot_before_after_summary(frame, summary_dir)
        plot_corrected_period_summary(frame, summary_dir)
        plot_model_variations(frame, summary_dir)
        plot_problem_bin_pull_diagnostics(
            results=results,
            output_dir=pull_diagnostic_dir,
            fit_min_gev2=args.fit_min,
            fit_max_gev2=args.fit_max,
        )

        recommended_after = frame[
            (frame["stage"] == "after")
            & frame["is_recommended"]
        ]
        successful_recommended_after = int(
            recommended_after["success"].sum()
        )
        accepted_recommended_after = int(
            recommended_after["fit_accepted"].sum()
        )
        unresolved_recommended_after = int(
            len(recommended_after) - accepted_recommended_after
        )

        print("")
        print("Deterministic Mx2 fit study complete.")
        print(f"Output directory: {output_dir}")
        print(f"Flat fit table:    {csv_path}")
        print(f"Compact JSON:      {json_path}")
        print(f"LaTeX tables:      {latex_path}")
        print(f"Status report:     {status_path}")
        print(
            "Corrected recommended fits: "
            f"{successful_recommended_after}/72 successful; "
            f"{accepted_recommended_after} accepted; "
            f"{unresolved_recommended_after} unresolved."
        )
        return 0
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        return 1
    # endtry


if __name__ == "__main__":
    raise SystemExit(main())
# endif
