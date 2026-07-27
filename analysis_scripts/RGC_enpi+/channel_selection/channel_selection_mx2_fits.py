#!/usr/bin/env python3
"""
channel_selection_mx2_fits_v30.py

Derive the RGC exclusive e pi+ missing-neutron Mx2 cuts used by downstream
analysis. The production cuts are always obtained from the ordinary non-ISR
NH3 and carbon samples. ISR samples are fitted only as a retained diagnostic
cross-check and never enter the exported production cut table. The complete
nominal and ISR analysis trees are retained under output/nominal and output/isr;
no fit plots, polynomial-method plots, carbon-subtraction diagnostics, tables,
or fit-status products are discarded after completion.

The physics fit remains the established carbon-assisted corrected-data fit:
for each of the fixed 4 xB by 6 (-tprime) bins, Su22, Fa22, and Sp23 share one
corrected Gaussian mean while retaining independent Gaussian widths, yields,
and residual-background parameters. The exported tight, nominal, and loose
windows are respectively mu +/- 2 sigma, mu +/- 3 sigma, and mu +/- 4 sigma
for each run period.

The public output is intentionally compact:

  output/channel_selection_mx2_fit_stability/
    nominal/
      complete nominal fit tables, plots, carbon diagnostics, and cut products
    isr/
      complete ISR fit tables, plots, carbon diagnostics, and cut products
    final_carbon_assisted_cuts/
      tables/final_carbon_assisted_mx2_cuts.json
      tables/final_carbon_assisted_mx2_cuts.csv
      exports/final_carbon_assisted_mx2_cuts.py
      exports/final_carbon_assisted_mx2_cuts.h
    diagnostics/
      plots/nominal_vs_isr_mu_sigma_v27.png
      plots/nominal_vs_isr_parameter_differences_v27.png
      tables/nominal_vs_isr_fit_comparison_v27.csv
      tables/nominal_vs_isr_fit_comparison_v27.json
    channel_selection_manifest.json

A v24 JSON compatibility alias is also written at the exact default path
currently consumed by determine_dilution_factor.py. Intermediate fit products
are created in a temporary work directory and removed after the final exports
and nominal-versus-ISR diagnostics have been written.

The branch tprime is negated explicitly before applying the requested
-minus-tprime binning.
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
import shutil
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, least_squares
from scipy.ndimage import gaussian_filter1d
import uproot


# =============================================================================
# Physics constants and fixed analysis binning
# =============================================================================

NEUTRON_MASS_GEV = 0.9395654133
NEUTRON_MASS2_GEV2 = NEUTRON_MASS_GEV**2
CARBON_RESIDUAL_THRESHOLD_GEV2 = 0.40

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

DEFAULT_CARBON_INPUTS: dict[str, dict[str, Path]] = {
    "su22": {
        "before": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "enpi+/rgc_su22_inb_C_epi+_2.root"
        ),
        "after": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "paper_versions/rgc_su22_inb_C_epi+_mom_corrections.root"
        ),
    },
    "fa22": {
        "before": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "enpi+/rgc_fa22_inb_C_epi+_2.root"
        ),
        "after": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "paper_versions/rgc_fa22_inb_C_epi+_mom_corrections.root"
        ),
    },
    "sp23": {
        "before": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "enpi+/rgc_sp23_inb_C_epi+_2.root"
        ),
        "after": Path(
            "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/"
            "paper_versions/rgc_sp23_inb_C_epi+_mom_corrections.root"
        ),
    },
}


DEFAULT_ISR_INPUTS: dict[str, dict[str, Path]] = {
    "su22": {
        "before": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+_ISR_2.root"),
        "after": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions/rgc_su22_inb_NH3_epi+_ISR_externalISR_mom_corrections.root"),
    },
    "fa22": {
        "before": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_NH3_epi+_ISR_2.root"),
        "after": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions/rgc_fa22_inb_NH3_epi+_ISR_externalISR_mom_corrections.root"),
    },
    "sp23": {
        "before": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_NH3_epi+_ISR_2.root"),
        "after": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions/rgc_sp23_inb_NH3_epi+_ISR_externalISR_mom_corrections.root"),
    },
}

DEFAULT_ISR_CARBON_INPUTS: dict[str, dict[str, Path]] = {
    "su22": {
        "before": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_C_epi+_ISR_2.root"),
        "after": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions/rgc_su22_inb_C_epi+_ISR_externalISR_mom_corrections.root"),
    },
    "fa22": {
        "before": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_fa22_inb_C_epi+_ISR_2.root"),
        "after": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions/rgc_fa22_inb_C_epi+_ISR_externalISR_mom_corrections.root"),
    },
    "sp23": {
        "before": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_sp23_inb_C_epi+_ISR_2.root"),
        "after": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions/rgc_sp23_inb_C_epi+_ISR_externalISR_mom_corrections.root"),
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

@dataclass(frozen=True)
class JointFitConfig:
    """Configuration for corrected-data simultaneous fits."""

    aicc_improvement_required: float
    local_improvement_fraction: float
    accepted_local_score_max: float
    marginal_local_score_max: float
    sideband_exclusion_min_gev2: float
    sideband_exclusion_max_gev2: float
    range_study_maxima_gev2: tuple[float, ...]
    enable_neighbor_mean_fallback: bool
    replicas: int
    replica_seed: int



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


def parse_float_list(value: str) -> tuple[float, ...]:
    """Parse a comma-separated list of finite floating-point values."""
    pieces = [piece.strip() for piece in value.split(",") if piece.strip()]
    if not pieces:
        raise argparse.ArgumentTypeError(
            "Expected at least one comma-separated floating-point value."
        )
    # endif
    try:
        values = tuple(float(piece) for piece in pieces)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Could not parse floating-point list: {value}"
        ) from exc
    # endtry
    if not all(np.isfinite(item) for item in values):
        raise argparse.ArgumentTypeError("All values must be finite.")
    # endif
    return values


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




def polynomial_parameter_count(background_model: str) -> int:
    """Return the number of polynomial coefficients."""
    return {
        "linear": 2,
        "quadratic": 3,
        "cubic": 4,
        "quartic": 5,
    }[background_model]


def sideband_background_initial(
    background_model: str,
    fit_x: np.ndarray,
    fit_y: np.ndarray,
    exclusion_min_gev2: float = 0.74,
    exclusion_max_gev2: float = 1.04,
) -> list[float]:
    """Estimate polynomial coefficients from lower and upper sidebands.

    Coefficients are returned in ascending powers of
    dx = Mx2 - neutron_mass_squared, matching the background functions.
    """
    fit_x = np.asarray(fit_x, dtype=float)
    fit_y = np.asarray(fit_y, dtype=float)
    coefficient_count = polynomial_parameter_count(background_model)
    degree = coefficient_count - 1

    sideband_mask = (
        (fit_x < exclusion_min_gev2)
        | (fit_x > exclusion_max_gev2)
    )
    sideband_mask &= np.isfinite(fit_x) & np.isfinite(fit_y)

    if np.count_nonzero(sideband_mask) < degree + 2:
        baseline = max(float(np.percentile(fit_y, 20.0)), 0.0)
        return [baseline, *([0.0] * degree)]
    # endif

    dx = fit_x[sideband_mask] - NEUTRON_MASS2_GEV2
    y = fit_y[sideband_mask]
    weights = 1.0 / np.sqrt(np.maximum(y, 1.0))

    try:
        descending = np.polyfit(
            dx,
            y,
            deg=degree,
            w=weights,
        )
        ascending = descending[::-1]
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        baseline = max(float(np.percentile(fit_y, 20.0)), 0.0)
        return [baseline, *([0.0] * degree)]
    # endtry

    if not np.all(np.isfinite(ascending)):
        baseline = max(float(np.percentile(fit_y, 20.0)), 0.0)
        return [baseline, *([0.0] * degree)]
    # endif

    return [float(value) for value in ascending]


def polynomial_bounds(
    background_model: str,
    background_scale: float,
) -> tuple[list[float], list[float]]:
    """Return broad deterministic bounds for polynomial coefficients."""
    scale = max(float(background_scale), 1.0)
    lower_by_order = {
        "linear": [-0.50 * scale, -30.0 * scale],
        "quadratic": [-0.50 * scale, -30.0 * scale, -150.0 * scale],
        "cubic": [
            -0.50 * scale,
            -30.0 * scale,
            -150.0 * scale,
            -750.0 * scale,
        ],
        "quartic": [
            -0.50 * scale,
            -30.0 * scale,
            -150.0 * scale,
            -750.0 * scale,
            -4000.0 * scale,
        ],
    }
    upper_by_order = {
        "linear": [3.00 * scale, 30.0 * scale],
        "quadratic": [3.00 * scale, 30.0 * scale, 150.0 * scale],
        "cubic": [
            3.00 * scale,
            30.0 * scale,
            150.0 * scale,
            750.0 * scale,
        ],
        "quartic": [
            3.00 * scale,
            30.0 * scale,
            150.0 * scale,
            750.0 * scale,
            4000.0 * scale,
        ],
    }
    return lower_by_order[background_model], upper_by_order[background_model]


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
    background_initial = sideband_background_initial(
        background_model=background_model,
        fit_x=fit_x,
        fit_y=fit_y,
    )
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
        initial = [*common_initial, *background_initial]
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
        initial = [*common_initial, *background_initial]
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
        initial = [*common_initial, *background_initial]
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
        initial = [*common_initial, *background_initial]
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
# Joint corrected-data fits
# =============================================================================

def histogram_for_job(
    job: FitJob,
    fit_max_override_gev2: float | None = None,
) -> dict[str, Any]:
    """Build a deterministic histogram and fit mask for one job."""
    values = np.asarray(job.values, dtype=float)
    values = values[np.isfinite(values)]
    counts, edges = np.histogram(
        values,
        bins=job.histogram_bins,
        range=(job.histogram_min_gev2, job.histogram_max_gev2),
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    fit_max = (
        float(fit_max_override_gev2)
        if fit_max_override_gev2 is not None
        else float(job.fit_max_gev2)
    )
    fit_mask = (
        (centers >= job.fit_min_gev2)
        & (centers <= fit_max)
    )
    return {
        "values": values,
        "counts": np.asarray(counts, dtype=float),
        "edges": np.asarray(edges, dtype=float),
        "centers": np.asarray(centers, dtype=float),
        "fit_mask": np.asarray(fit_mask, dtype=bool),
        "fit_min_gev2": float(job.fit_min_gev2),
        "fit_max_gev2": fit_max,
        "histogram_bin_width_gev2": float(edges[1] - edges[0]),
    }


def joint_parameter_layout(
    background_model: str,
    periods: tuple[str, ...],
    fixed_mean_gev2: float | None,
) -> dict[str, Any]:
    """Describe the deterministic joint-parameter vector."""
    coefficient_count = polynomial_parameter_count(background_model)
    cursor = 0
    mean_index: int | None = None
    if fixed_mean_gev2 is None:
        mean_index = cursor
        cursor += 1
    # endif

    period_slices: dict[str, dict[str, Any]] = {}
    for period in periods:
        yield_index = cursor
        sigma_index = cursor + 1
        coefficient_slice = slice(
            cursor + 2,
            cursor + 2 + coefficient_count,
        )
        period_slices[period] = {
            "yield_index": yield_index,
            "sigma_index": sigma_index,
            "coefficient_slice": coefficient_slice,
        }
        cursor += 2 + coefficient_count
    # endfor

    return {
        "mean_index": mean_index,
        "period_slices": period_slices,
        "number_of_parameters": cursor,
        "coefficient_count": coefficient_count,
    }


def evaluate_polynomial_coefficients(
    background_model: str,
    x: np.ndarray,
    coefficients: np.ndarray,
) -> np.ndarray:
    """Evaluate one polynomial from ascending coefficients."""
    if background_model == "linear":
        return linear_background(x, *coefficients)
    # endif
    if background_model == "quadratic":
        return quadratic_background(x, *coefficients)
    # endif
    if background_model == "cubic":
        return cubic_background(x, *coefficients)
    # endif
    if background_model == "quartic":
        return quartic_background(x, *coefficients)
    # endif
    raise ValueError(background_model)


def joint_initial_values_and_bounds(
    background_model: str,
    jobs_by_period: dict[str, FitJob],
    histograms: dict[str, dict[str, Any]],
    periods: tuple[str, ...],
    corrected_mean_max_gev2: float,
    fixed_mean_gev2: float | None,
    sideband_exclusion_min_gev2: float,
    sideband_exclusion_max_gev2: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Build sideband-initialized joint parameters and bounds."""
    layout = joint_parameter_layout(
        background_model=background_model,
        periods=periods,
        fixed_mean_gev2=fixed_mean_gev2,
    )

    initial: list[float] = []
    lower: list[float] = []
    upper: list[float] = []

    if fixed_mean_gev2 is None:
        mean_low = max(
            max(job.fit_min_gev2 for job in jobs_by_period.values()),
            NEUTRON_MASS2_GEV2 - 0.08,
        )
        mean_high = min(
            min(job.fit_max_gev2 for job in jobs_by_period.values()),
            corrected_mean_max_gev2,
        )
        initial.append(
            float(
                np.clip(
                    NEUTRON_MASS2_GEV2,
                    mean_low + 1.0e-4,
                    mean_high - 1.0e-4,
                )
            )
        )
        lower.append(mean_low)
        upper.append(mean_high)
    # endif

    for period in periods:
        histogram = histograms[period]
        fit_mask = histogram["fit_mask"]
        fit_x = histogram["centers"][fit_mask]
        fit_y = histogram["counts"][fit_mask]
        bin_width = histogram["histogram_bin_width_gev2"]

        baseline = max(float(np.percentile(fit_y, 20.0)), 0.0)
        peak_mask = (
            (fit_x >= NEUTRON_MASS2_GEV2 - 0.16)
            & (fit_x <= NEUTRON_MASS2_GEV2 + 0.16)
        )
        peak_height = max(
            float(np.max(fit_y[peak_mask]) - baseline)
            if np.any(peak_mask)
            else float(np.max(fit_y) - baseline),
            1.0,
        )
        sigma_guess = 0.075
        yield_guess = (
            peak_height
            * sigma_guess
            * math.sqrt(2.0 * math.pi)
            / bin_width
        )
        total_counts = max(float(np.sum(fit_y)), 1.0)
        yield_upper = max(3.0 * total_counts, 10.0 * yield_guess, 100.0)
        background_scale = max(float(np.max(fit_y)), 1.0)
        background_initial = sideband_background_initial(
            background_model=background_model,
            fit_x=fit_x,
            fit_y=fit_y,
            exclusion_min_gev2=sideband_exclusion_min_gev2,
            exclusion_max_gev2=sideband_exclusion_max_gev2,
        )
        background_lower, background_upper = polynomial_bounds(
            background_model=background_model,
            background_scale=background_scale,
        )

        initial.extend([yield_guess, sigma_guess, *background_initial])
        lower.extend([0.0, 0.010, *background_lower])
        upper.extend([yield_upper, 0.180, *background_upper])
    # endfor

    return (
        np.asarray(initial, dtype=float),
        np.asarray(lower, dtype=float),
        np.asarray(upper, dtype=float),
        layout,
    )


def fit_joint_background_model(
    background_model: str,
    jobs: list[FitJob],
    config: JointFitConfig,
    fit_max_override_gev2: float | None = None,
    fixed_mean_gev2: float | None = None,
    counts_override: dict[str, np.ndarray] | None = None,
) -> dict[str, Any]:
    """Fit corrected Su22/Fa22/Sp23 spectra with one shared Gaussian mean."""
    jobs_by_period = {job.period: job for job in jobs}
    periods = tuple(
        period for period in ("su22", "fa22", "sp23")
        if period in jobs_by_period
    )
    if len(periods) != 3:
        return {
            "background_model": background_model,
            "success": False,
            "status": "joint_fit_requires_three_periods",
            "review_reasons": ["missing_period"],
            "period_results": {},
        }
    # endif

    histograms = {
        period: histogram_for_job(
            jobs_by_period[period],
            fit_max_override_gev2=fit_max_override_gev2,
        )
        for period in periods
    }
    if counts_override is not None:
        for period in periods:
            histograms[period]["counts"] = np.asarray(
                counts_override[period],
                dtype=float,
            )
        # endfor
    # endif

    if any(
        histograms[period]["values"].size
        < jobs_by_period[period].minimum_events
        for period in periods
    ):
        return {
            "background_model": background_model,
            "success": False,
            "status": "insufficient_events",
            "review_reasons": ["insufficient_events"],
            "period_results": {},
        }
    # endif

    try:
        initial, lower, upper, layout = joint_initial_values_and_bounds(
            background_model=background_model,
            jobs_by_period=jobs_by_period,
            histograms=histograms,
            periods=periods,
            corrected_mean_max_gev2=min(
                job.corrected_mean_max_gev2 for job in jobs
            ),
            fixed_mean_gev2=fixed_mean_gev2,
            sideband_exclusion_min_gev2=(
                config.sideband_exclusion_min_gev2
            ),
            sideband_exclusion_max_gev2=(
                config.sideband_exclusion_max_gev2
            ),
        )

        def residual_function(parameters: np.ndarray) -> np.ndarray:
            if fixed_mean_gev2 is None:
                mean_gev2 = float(
                    parameters[layout["mean_index"]]
                )
            else:
                mean_gev2 = float(fixed_mean_gev2)
            # endif

            residual_blocks: list[np.ndarray] = []
            for period in periods:
                period_layout = layout["period_slices"][period]
                histogram = histograms[period]
                mask = histogram["fit_mask"]
                fit_x = histogram["centers"][mask]
                fit_y = histogram["counts"][mask]
                fit_error = np.sqrt(np.maximum(fit_y, 1.0))
                signal_yield = float(
                    parameters[period_layout["yield_index"]]
                )
                sigma_gev2 = float(
                    parameters[period_layout["sigma_index"]]
                )
                coefficients = parameters[
                    period_layout["coefficient_slice"]
                ]
                predicted = gaussian_signal_from_yield(
                    fit_x,
                    signal_yield,
                    mean_gev2,
                    sigma_gev2,
                    histogram["histogram_bin_width_gev2"],
                ) + evaluate_polynomial_coefficients(
                    background_model,
                    fit_x,
                    coefficients,
                )
                residual_blocks.append((fit_y - predicted) / fit_error)
            # endfor
            return np.concatenate(residual_blocks)
        # enddef

        optimization = least_squares(
            residual_function,
            x0=initial,
            bounds=(lower, upper),
            max_nfev=250000,
            method="trf",
            x_scale="jac",
        )
        if not optimization.success:
            raise RuntimeError(optimization.message)
        # endif

        parameters = np.asarray(optimization.x, dtype=float)
        pulls_all = residual_function(parameters)
        chi2 = float(np.sum(pulls_all**2))
        number_of_points = int(pulls_all.size)
        number_of_parameters = int(parameters.size)
        ndf = number_of_points - number_of_parameters
        chi2_ndf = chi2 / ndf if ndf > 0 else math.nan

        jacobian = np.asarray(optimization.jac, dtype=float)
        covariance = np.linalg.pinv(jacobian.T @ jacobian)
        errors = np.sqrt(
            np.where(np.diag(covariance) >= 0.0, np.diag(covariance), np.nan)
        )

        if number_of_points > number_of_parameters + 1:
            aicc = (
                chi2
                + 2.0 * number_of_parameters
                + (
                    2.0
                    * number_of_parameters
                    * (number_of_parameters + 1)
                    / (
                        number_of_points
                        - number_of_parameters
                        - 1
                    )
                )
            )
        else:
            aicc = math.inf
        # endif
    except (
        RuntimeError,
        ValueError,
        np.linalg.LinAlgError,
        FloatingPointError,
        OverflowError,
    ) as exc:
        return {
            "background_model": background_model,
            "success": False,
            "status": f"joint_fit_failed:{type(exc).__name__}",
            "review_reasons": ["fit_failed"],
            "period_results": {},
        }
    # endtry

    if fixed_mean_gev2 is None:
        mean_index = layout["mean_index"]
        mean_gev2 = float(parameters[mean_index])
        mean_error_gev2 = float(errors[mean_index])
    else:
        mean_gev2 = float(fixed_mean_gev2)
        mean_error_gev2 = 0.0
    # endif

    joint_review_reasons: list[str] = []
    if not np.all(np.isfinite(parameters)):
        joint_review_reasons.append("nonfinite_parameters")
    # endif
    if not np.all(np.isfinite(covariance)):
        joint_review_reasons.append("nonfinite_covariance")
    # endif
    if fixed_mean_gev2 is None and parameter_at_bound(
        mean_gev2,
        float(lower[layout["mean_index"]]),
        float(upper[layout["mean_index"]]),
    ):
        joint_review_reasons.append("mean_at_bound")
    # endif
    if not np.isfinite(mean_error_gev2) or mean_error_gev2 > 0.025:
        joint_review_reasons.append("large_mean_error")
    # endif

    period_results: dict[str, dict[str, Any]] = {}
    local_chi2_total = 0.0
    local_npoints_total = 0

    for period in periods:
        period_layout = layout["period_slices"][period]
        histogram = histograms[period]
        mask = histogram["fit_mask"]
        fit_x = histogram["centers"][mask]
        fit_y = histogram["counts"][mask]
        fit_error = np.sqrt(np.maximum(fit_y, 1.0))
        signal_yield = float(
            parameters[period_layout["yield_index"]]
        )
        sigma_gev2 = float(
            parameters[period_layout["sigma_index"]]
        )
        coefficients = np.asarray(
            parameters[period_layout["coefficient_slice"]],
            dtype=float,
        )
        signal_yield_error = float(
            errors[period_layout["yield_index"]]
        )
        sigma_error_gev2 = float(
            errors[period_layout["sigma_index"]]
        )

        fit_signal = gaussian_signal_from_yield(
            fit_x,
            signal_yield,
            mean_gev2,
            sigma_gev2,
            histogram["histogram_bin_width_gev2"],
        )
        fit_background = evaluate_polynomial_coefficients(
            background_model,
            fit_x,
            coefficients,
        )
        fit_total = fit_signal + fit_background
        pulls = (fit_y - fit_total) / fit_error
        period_chi2 = float(np.sum(pulls**2))
        period_npoints = int(fit_x.size)
        period_parameter_count_approx = (
            2 + polynomial_parameter_count(background_model)
            + 1.0 / len(periods)
        )
        period_ndf_approx = max(
            int(round(period_npoints - period_parameter_count_approx)),
            1,
        )
        period_chi2_ndf_approx = period_chi2 / period_ndf_approx

        local_mask = (
            (fit_x >= mean_gev2 - 2.0 * sigma_gev2)
            & (fit_x <= mean_gev2 + 2.0 * sigma_gev2)
        )
        local_chi2 = float(np.sum(pulls[local_mask] ** 2))
        local_npoints = int(np.count_nonzero(local_mask))
        local_score = (
            local_chi2 / local_npoints
            if local_npoints > 0
            else math.nan
        )
        local_ndf_approx = max(local_npoints - 2, 1)
        local_chi2_ndf_approx = local_chi2 / local_ndf_approx
        local_chi2_total += local_chi2
        local_npoints_total += local_npoints

        dense_x = np.linspace(
            histogram["fit_min_gev2"],
            histogram["fit_max_gev2"],
            2000,
        )
        dense_background = evaluate_polynomial_coefficients(
            background_model,
            dense_x,
            coefficients,
        )
        region_2sigma = (
            (dense_x >= mean_gev2 - 2.0 * sigma_gev2)
            & (dense_x <= mean_gev2 + 2.0 * sigma_gev2)
        )
        region_3sigma = (
            (dense_x >= mean_gev2 - 3.0 * sigma_gev2)
            & (dense_x <= mean_gev2 + 3.0 * sigma_gev2)
        )
        background_minimum = float(np.min(dense_background))
        background_minimum_2sigma = (
            float(np.min(dense_background[region_2sigma]))
            if np.any(region_2sigma)
            else math.nan
        )
        background_minimum_3sigma = (
            float(np.min(dense_background[region_3sigma]))
            if np.any(region_3sigma)
            else math.nan
        )
        background_scale = max(float(np.max(fit_y)), 1.0)
        negativity_tolerance = max(1.0, 0.01 * background_scale)

        period_review_reasons = list(joint_review_reasons)
        if parameter_at_bound(
            signal_yield,
            float(lower[period_layout["yield_index"]]),
            float(upper[period_layout["yield_index"]]),
        ):
            period_review_reasons.append("signal_yield_at_bound")
        # endif
        if parameter_at_bound(
            sigma_gev2,
            float(lower[period_layout["sigma_index"]]),
            float(upper[period_layout["sigma_index"]]),
        ):
            period_review_reasons.append("sigma_at_bound")
        # endif
        if (
            np.isfinite(background_minimum_3sigma)
            and background_minimum_3sigma < -negativity_tolerance
        ):
            period_review_reasons.append(
                "negative_background_in_3sigma_region"
            )
        # endif
        if background_minimum < -negativity_tolerance:
            period_review_reasons.append(
                "negative_background_global_diagnostic"
            )
        # endif
        if not np.isfinite(sigma_error_gev2) or sigma_error_gev2 > 0.025:
            period_review_reasons.append("large_sigma_error")
        # endif
        if not np.isfinite(local_score) or local_score >= 3.0:
            period_review_reasons.append("poor_signal_region_chi2")
        # endif

        signal_plus_background = np.maximum(fit_total, 1.0)
        peak_significance_proxy = (
            signal_yield
            / math.sqrt(float(np.sum(signal_plus_background)))
            if signal_yield > 0.0
            else 0.0
        )
        if peak_significance_proxy < 5.0:
            period_review_reasons.append("weak_peak")
        # endif

        local_parameters = np.asarray(
            [
                signal_yield,
                mean_gev2,
                sigma_gev2,
                *coefficients.tolist(),
            ],
            dtype=float,
        )
        local_errors = np.asarray(
            [
                signal_yield_error,
                mean_error_gev2,
                sigma_error_gev2,
                *errors[
                    period_layout["coefficient_slice"]
                ].tolist(),
            ],
            dtype=float,
        )

        period_results[period] = {
            "background_model": background_model,
            "success": True,
            "status": (
                "success"
                if not period_review_reasons
                else "success_flagged_for_review"
            ),
            "review_reasons": sorted(set(period_review_reasons)),
            "signal_yield": signal_yield,
            "signal_yield_error": signal_yield_error,
            "mean_gev2": mean_gev2,
            "mean_error_gev2": mean_error_gev2,
            "sigma_gev2": sigma_gev2,
            "sigma_error_gev2": sigma_error_gev2,
            "chi2": period_chi2,
            "ndf": period_ndf_approx,
            "chi2_ndf": period_chi2_ndf_approx,
            "signal_region_definition": "shared mean +/- 2 period sigma",
            "signal_region_chi2": local_chi2,
            "signal_region_npoints": local_npoints,
            "signal_region_ndf_approx": local_ndf_approx,
            "signal_region_chi2_ndf_approx": local_chi2_ndf_approx,
            "signal_region_chi2_per_bin": local_score,
            "aicc": float(aicc),
            "joint_aicc": float(aicc),
            "joint_chi2": chi2,
            "joint_ndf": int(ndf),
            "joint_chi2_ndf": chi2_ndf,
            "peak_significance_proxy": peak_significance_proxy,
            "background_minimum_in_fit_range": background_minimum,
            "background_minimum_in_2sigma_region": (
                background_minimum_2sigma
            ),
            "background_minimum_in_3sigma_region": (
                background_minimum_3sigma
            ),
            "background_negativity_tolerance": negativity_tolerance,
            "parameters": local_parameters.tolist(),
            "parameter_errors": local_errors.tolist(),
            "covariance": covariance.tolist(),
            "shared_mean": True,
            "shared_mean_fixed": fixed_mean_gev2 is not None,
            "parameter_source": (
                "fixed_neighbor_mean"
                if fixed_mean_gev2 is not None
                else "joint_fit"
            ),
        }
    # endfor

    aggregate_local_score = (
        local_chi2_total / local_npoints_total
        if local_npoints_total > 0
        else math.nan
    )
    maximum_period_local_score = max(
        (
            result["signal_region_chi2_per_bin"]
            for result in period_results.values()
            if np.isfinite(result["signal_region_chi2_per_bin"])
        ),
        default=math.nan,
    )

    hard_rejection_reasons = {
        "nonfinite_parameters",
        "nonfinite_covariance",
        "signal_yield_at_bound",
        "mean_at_bound",
        "sigma_at_bound",
        "negative_background_in_3sigma_region",
    }
    viable = all(
        result["success"]
        and not bool(
            set(result["review_reasons"])
            & hard_rejection_reasons
        )
        for result in period_results.values()
    )

    return {
        "background_model": background_model,
        "success": True,
        "status": "success" if viable else "success_flagged_for_review",
        "review_reasons": sorted(set(joint_review_reasons)),
        "period_results": period_results,
        "periods": list(periods),
        "shared_mean_gev2": mean_gev2,
        "shared_mean_error_gev2": mean_error_gev2,
        "shared_mean_fixed": fixed_mean_gev2 is not None,
        "joint_chi2": chi2,
        "joint_ndf": int(ndf),
        "joint_chi2_ndf": chi2_ndf,
        "aicc": float(aicc),
        "aggregate_local_score": aggregate_local_score,
        "maximum_period_local_score": maximum_period_local_score,
        "viable": bool(viable),
        "fit_max_gev2": (
            fit_max_override_gev2
            if fit_max_override_gev2 is not None
            else jobs[0].fit_max_gev2
        ),
        "parameters": parameters.tolist(),
        "parameter_errors": errors.tolist(),
        "covariance": covariance.tolist(),
    }


def choose_joint_background_model(
    models: dict[str, dict[str, Any]],
    config: JointFitConfig,
) -> tuple[str | None, str]:
    """Choose the lowest justified corrected-data polynomial order."""
    candidates = ("quadratic", "cubic", "quartic")
    viable = [
        model_name
        for model_name in candidates
        if (
            models[model_name].get("success", False)
            and models[model_name].get("viable", False)
            and np.isfinite(
                models[model_name].get(
                    "aggregate_local_score",
                    math.nan,
                )
            )
            and np.isfinite(models[model_name].get("aicc", math.nan))
        )
    ]
    if not viable:
        successful = [
            model_name
            for model_name in candidates
            if (
                models[model_name].get("success", False)
                and np.isfinite(
                    models[model_name].get(
                        "aggregate_local_score",
                        math.nan,
                    )
                )
            )
        ]
        if not successful:
            return None, "no successful corrected joint polynomial fit"
        # endif
        selected = min(
            successful,
            key=lambda name: (
                models[name]["maximum_period_local_score"],
                models[name]["aicc"],
            ),
        )
        return (
            selected,
            "no model passed hard viability checks; selected the successful "
            "model with the smallest worst-period local score",
        )
    # endif

    selected = viable[0]
    selection_steps = [f"started from {selected}"]
    selected_index = candidates.index(selected)

    for next_name in candidates[selected_index + 1:]:
        if next_name not in viable:
            continue
        # endif
        current = models[selected]
        challenger = models[next_name]
        delta_aicc = challenger["aicc"] - current["aicc"]
        current_local = current["aggregate_local_score"]
        challenger_local = challenger["aggregate_local_score"]
        improvement_fraction = (
            (current_local - challenger_local) / current_local
            if current_local > 0.0
            else 0.0
        )
        earned_complexity = (
            delta_aicc <= -config.aicc_improvement_required
            and improvement_fraction
            >= config.local_improvement_fraction
        )
        rescue_unresolved = (
            current["maximum_period_local_score"]
            >= config.marginal_local_score_max
            and challenger["maximum_period_local_score"]
            < current["maximum_period_local_score"]
            and delta_aicc <= -config.aicc_improvement_required
        )

        if earned_complexity or rescue_unresolved:
            selected = next_name
            selection_steps.append(
                f"promoted to {next_name}: delta AICc={delta_aicc:.2f}, "
                f"aggregate local improvement={100.0 * improvement_fraction:.1f}%"
            )
        else:
            selection_steps.append(
                f"stopped before {next_name}: delta AICc={delta_aicc:.2f}, "
                f"aggregate local improvement={100.0 * improvement_fraction:.1f}%"
            )
            break
        # endif
    # endfor

    return selected, "; ".join(selection_steps)


def classify_joint_fit(
    selected_joint: dict[str, Any],
    config: JointFitConfig,
) -> tuple[str, bool, str]:
    """Classify one corrected joint fit."""
    if not selected_joint.get("success", False):
        return "unresolved", False, "selected joint fit failed"
    # endif
    if not selected_joint.get("viable", False):
        return (
            "unresolved",
            False,
            "selected joint fit failed hard viability criteria",
        )
    # endif

    worst_score = selected_joint.get(
        "maximum_period_local_score",
        math.nan,
    )
    if not np.isfinite(worst_score):
        return "unresolved", False, "nonfinite local fit-quality score"
    # endif
    if worst_score < config.accepted_local_score_max:
        return (
            "accepted",
            True,
            "all periods have local chi2 per signal-region bin below "
            f"{config.accepted_local_score_max:.2f}",
        )
    # endif
    if worst_score < config.marginal_local_score_max:
        return (
            "marginal",
            True,
            "hard viability criteria pass and the worst-period local score "
            f"is below {config.marginal_local_score_max:.2f}",
        )
    # endif
    return (
        "unresolved",
        False,
        "worst-period local score exceeds the marginal threshold",
    )


def base_result_from_job(
    job: FitJob,
    histogram: dict[str, Any],
) -> dict[str, Any]:
    """Create the common per-period result record."""
    return {
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
        "number_of_events": int(histogram["values"].size),
        "histogram_entries": int(np.sum(histogram["counts"])),
        "histogram_bin_width_gev2": (
            histogram["histogram_bin_width_gev2"]
        ),
        "counts": histogram["counts"].astype(int).tolist(),
        "edges": histogram["edges"].tolist(),
        "models": {},
    }


def assemble_joint_period_results(
    jobs: list[FitJob],
    joint_models: dict[str, dict[str, Any]],
    selected_model: str | None,
    recommendation_reason: str,
    config: JointFitConfig,
    fit_range_stability: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Convert one joint corrected fit into three legacy-compatible records."""
    jobs_by_period = {job.period: job for job in jobs}
    histograms = {
        period: histogram_for_job(job)
        for period, job in jobs_by_period.items()
    }

    if selected_model is None:
        output: list[dict[str, Any]] = []
        for period, job in jobs_by_period.items():
            base = base_result_from_job(job, histograms[period])
            base["models"] = {
                model_name: joint_models[model_name].get(
                    "period_results",
                    {},
                ).get(
                    period,
                    {
                        "background_model": model_name,
                        "success": False,
                        "status": joint_models[model_name].get(
                            "status",
                            "fit_failed",
                        ),
                        "review_reasons": joint_models[model_name].get(
                            "review_reasons",
                            ["fit_failed"],
                        ),
                    },
                )
                for model_name in BACKGROUND_MODELS
            }
            base["recommended_background_model"] = None
            base["recommended_background_order"] = None
            base["recommendation_reason"] = recommendation_reason
            base["fit_quality_class"] = "unresolved"
            base["fit_accepted"] = False
            base["status"] = "unresolved"
            base["fit_acceptance_reason"] = (
                "no corrected joint model could be selected"
            )
            output.append(base)
        # endfor
        return output
    # endif

    selected_joint = joint_models[selected_model]
    quality_class, fit_accepted, acceptance_reason = classify_joint_fit(
        selected_joint=selected_joint,
        config=config,
    )

    output = []
    for period, job in jobs_by_period.items():
        base = base_result_from_job(job, histograms[period])
        base["models"] = {}
        for model_name in BACKGROUND_MODELS:
            joint_model = joint_models[model_name]
            base["models"][model_name] = joint_model.get(
                "period_results",
                {},
            ).get(
                period,
                {
                    "background_model": model_name,
                    "success": False,
                    "status": joint_model.get("status", "fit_failed"),
                    "review_reasons": joint_model.get(
                        "review_reasons",
                        ["fit_failed"],
                    ),
                },
            )
        # endfor

        base["recommended_background_model"] = selected_model
        base["recommended_background_order"] = {
            "quadratic": 2,
            "cubic": 3,
            "quartic": 4,
        }[selected_model]
        base["recommendation_reason"] = recommendation_reason
        base["fit_quality_class"] = quality_class
        base["fit_accepted"] = fit_accepted
        base["status"] = quality_class
        base["fit_acceptance_reason"] = acceptance_reason
        base["joint_fit"] = {
            "shared_mean_gev2": selected_joint["shared_mean_gev2"],
            "shared_mean_error_gev2": (
                selected_joint["shared_mean_error_gev2"]
            ),
            "shared_mean_fixed": selected_joint["shared_mean_fixed"],
            "background_model": selected_model,
            "joint_chi2": selected_joint["joint_chi2"],
            "joint_ndf": selected_joint["joint_ndf"],
            "joint_chi2_ndf": selected_joint["joint_chi2_ndf"],
            "joint_aicc": selected_joint["aicc"],
            "aggregate_local_score": (
                selected_joint["aggregate_local_score"]
            ),
            "maximum_period_local_score": (
                selected_joint["maximum_period_local_score"]
            ),
            "parameter_source": base["models"][selected_model].get(
                "parameter_source",
                "joint_fit",
            ),
        }
        if fit_range_stability is not None:
            base["fit_range_stability"] = fit_range_stability
        # endif

        nominal = base["models"][selected_model]
        if nominal.get("success", False):
            base["recommended_selection_windows_gev2"] = {
                label: {
                    "minimum": (
                        nominal["mean_gev2"]
                        - multiplier * nominal["sigma_gev2"]
                    ),
                    "maximum": (
                        nominal["mean_gev2"]
                        + multiplier * nominal["sigma_gev2"]
                    ),
                }
                for label, multiplier in (
                    ("1.5_sigma", 1.5),
                    ("2.0_sigma", 2.0),
                    ("3.0_sigma", 3.0),
                )
            }

            differences: dict[str, Any] = {}
            for alternative in BACKGROUND_MODELS:
                if alternative == selected_model:
                    continue
                # endif
                alt = base["models"][alternative]
                if alt.get("success", False):
                    differences[alternative] = {
                        "delta_mean_gev2": (
                            alt["mean_gev2"] - nominal["mean_gev2"]
                        ),
                        "delta_sigma_gev2": (
                            alt["sigma_gev2"] - nominal["sigma_gev2"]
                        ),
                        "delta_signal_yield": (
                            alt["signal_yield"]
                            - nominal["signal_yield"]
                        ),
                        "relative_signal_yield_difference": (
                            (
                                alt["signal_yield"]
                                - nominal["signal_yield"]
                            )
                            / nominal["signal_yield"]
                            if nominal["signal_yield"] != 0.0
                            else math.nan
                        ),
                        "delta_aicc": (
                            alt["aicc"] - nominal["aicc"]
                        ),
                    }
                else:
                    differences[alternative] = {
                        "status": alt.get("status", "fit_failed")
                    }
                # endif
            # endfor
            base["model_differences_from_recommended"] = differences
        # endif

        output.append(base)
    # endfor

    return output


def fit_one_joint_corrected_bin(
    jobs: list[FitJob],
    config: JointFitConfig,
) -> list[dict[str, Any]]:
    """Fit one corrected kinematic bin across all three periods."""
    joint_models = {
        model_name: fit_joint_background_model(
            background_model=model_name,
            jobs=jobs,
            config=config,
        )
        for model_name in BACKGROUND_MODELS
    }
    selected_model, reason = choose_joint_background_model(
        models=joint_models,
        config=config,
    )

    range_stability: dict[str, Any] = {}
    if selected_model is not None:
        nominal = joint_models[selected_model]
        for fit_max in config.range_study_maxima_gev2:
            varied = fit_joint_background_model(
                background_model=selected_model,
                jobs=jobs,
                config=config,
                fit_max_override_gev2=fit_max,
            )
            if varied.get("success", False):
                range_stability[f"{fit_max:.2f}"] = {
                    "fit_max_gev2": fit_max,
                    "shared_mean_gev2": varied["shared_mean_gev2"],
                    "delta_shared_mean_gev2": (
                        varied["shared_mean_gev2"]
                        - nominal["shared_mean_gev2"]
                    ),
                    "periods": {
                        period: {
                            "sigma_gev2": varied["period_results"][
                                period
                            ]["sigma_gev2"],
                            "delta_sigma_gev2": (
                                varied["period_results"][period][
                                    "sigma_gev2"
                                ]
                                - nominal["period_results"][period][
                                    "sigma_gev2"
                                ]
                            ),
                            "signal_yield": varied["period_results"][
                                period
                            ]["signal_yield"],
                            "relative_signal_yield_change": (
                                (
                                    varied["period_results"][period][
                                        "signal_yield"
                                    ]
                                    - nominal["period_results"][period][
                                        "signal_yield"
                                    ]
                                )
                                / nominal["period_results"][period][
                                    "signal_yield"
                                ]
                                if nominal["period_results"][period][
                                    "signal_yield"
                                ] != 0.0
                                else math.nan
                            ),
                        }
                        for period in ("su22", "fa22", "sp23")
                    },
                }
            else:
                range_stability[f"{fit_max:.2f}"] = {
                    "fit_max_gev2": fit_max,
                    "status": varied.get("status", "fit_failed"),
                }
            # endif
        # endfor
    # endif

    return assemble_joint_period_results(
        jobs=jobs,
        joint_models=joint_models,
        selected_model=selected_model,
        recommendation_reason=reason,
        config=config,
        fit_range_stability=range_stability,
    )


def neighboring_shared_mean(
    x_index: int,
    t_index: int,
    corrected_results: list[dict[str, Any]],
) -> float | None:
    """Interpolate a shared mean from accepted neighboring t-prime bins."""
    unique_by_bin: dict[tuple[int, int], dict[str, Any]] = {}
    for item in corrected_results:
        key = (item["x_index"], item["t_index"])
        unique_by_bin.setdefault(key, item)
    # endfor

    neighbors: list[tuple[int, float, float]] = []
    for (candidate_x, candidate_t), item in unique_by_bin.items():
        if candidate_x != x_index:
            continue
        # endif
        if item.get("fit_quality_class") == "unresolved":
            continue
        # endif
        joint = item.get("joint_fit", {})
        mean = joint.get("shared_mean_gev2", math.nan)
        error = joint.get("shared_mean_error_gev2", math.nan)
        if np.isfinite(mean):
            neighbors.append(
                (
                    candidate_t,
                    float(mean),
                    float(error) if np.isfinite(error) else 0.010,
                )
            )
        # endif
    # endfor

    if not neighbors:
        return None
    # endif

    neighbors.sort(key=lambda value: value[0])
    lower = [value for value in neighbors if value[0] < t_index]
    upper = [value for value in neighbors if value[0] > t_index]

    if lower and upper:
        t0, mean0, _ = lower[-1]
        t1, mean1, _ = upper[0]
        fraction = (t_index - t0) / (t1 - t0)
        return float(mean0 + fraction * (mean1 - mean0))
    # endif

    nearest = min(neighbors, key=lambda value: abs(value[0] - t_index))
    return float(nearest[1])


def apply_neighbor_mean_fallbacks(
    corrected_results: list[dict[str, Any]],
    corrected_groups: list[list[FitJob]],
    config: JointFitConfig,
) -> list[dict[str, Any]]:
    """Refit unresolved corrected bins with a neighboring shared mean fixed."""
    if not config.enable_neighbor_mean_fallback:
        return corrected_results
    # endif

    results_by_key: dict[tuple[str, int, int], dict[str, Any]] = {
        (item["period"], item["x_index"], item["t_index"]): item
        for item in corrected_results
    }

    for jobs in corrected_groups:
        x_index = jobs[0].x_index
        t_index = jobs[0].t_index
        representative = results_by_key[
            ("su22", x_index, t_index)
        ]
        if representative.get("fit_quality_class") != "unresolved":
            continue
        # endif

        fixed_mean = neighboring_shared_mean(
            x_index=x_index,
            t_index=t_index,
            corrected_results=list(results_by_key.values()),
        )
        if fixed_mean is None:
            continue
        # endif

        selected_model = representative.get(
            "recommended_background_model"
        )
        if selected_model not in {"quadratic", "cubic", "quartic"}:
            selected_model = "quartic"
        # endif

        fixed_joint = fit_joint_background_model(
            background_model=selected_model,
            jobs=jobs,
            config=config,
            fixed_mean_gev2=fixed_mean,
        )
        if not fixed_joint.get("success", False):
            continue
        # endif

        quality_class, accepted, reason = classify_joint_fit(
            selected_joint=fixed_joint,
            config=config,
        )
        if quality_class == "unresolved":
            continue
        # endif

        original_models = {
            model_name: {
                "background_model": model_name,
                "success": False,
                "status": "not_refit_in_neighbor_fallback",
                "review_reasons": ["not_refit_in_neighbor_fallback"],
                "period_results": {},
            }
            for model_name in BACKGROUND_MODELS
        }
        original_models[selected_model] = fixed_joint
        replacement = assemble_joint_period_results(
            jobs=jobs,
            joint_models=original_models,
            selected_model=selected_model,
            recommendation_reason=(
                "original free-mean joint fit was unresolved; refit with "
                f"neighbor-interpolated shared mean fixed at {fixed_mean:.6f} GeV^2"
            ),
            config=config,
            fit_range_stability=representative.get(
                "fit_range_stability",
                {},
            ),
        )
        for item in replacement:
            item["fit_quality_class"] = quality_class
            item["fit_accepted"] = accepted
            item["status"] = quality_class
            item["fit_acceptance_reason"] = reason
            item["fallback_applied"] = True
            results_by_key[
                (item["period"], item["x_index"], item["t_index"])
            ] = item
        # endfor
    # endfor

    return sorted(
        results_by_key.values(),
        key=lambda item: (
            item["period"],
            item["stage"],
            item["x_index"],
            item["t_index"],
        ),
    )


def run_joint_replicas(
    corrected_results: list[dict[str, Any]],
    corrected_groups: list[list[FitJob]],
    config: JointFitConfig,
) -> None:
    """Run optional fixed-model Poisson replicas in place."""
    if config.replicas <= 0:
        return
    # endif

    rng = np.random.default_rng(config.replica_seed)
    results_by_key = {
        (item["period"], item["x_index"], item["t_index"]): item
        for item in corrected_results
    }

    for group_index, jobs in enumerate(corrected_groups, start=1):
        representative = results_by_key[
            ("su22", jobs[0].x_index, jobs[0].t_index)
        ]
        model_name = representative.get("recommended_background_model")
        if model_name not in {"quadratic", "cubic", "quartic"}:
            continue
        # endif

        fixed_mean = None
        parameter_source = representative.get(
            "joint_fit",
            {},
        ).get("parameter_source", "joint_fit")
        if parameter_source == "fixed_neighbor_mean":
            fixed_mean = representative["joint_fit"]["shared_mean_gev2"]
        # endif

        histograms = {
            job.period: histogram_for_job(job)
            for job in jobs
        }
        nominal_period_results = {
            period: results_by_key[
                (period, jobs[0].x_index, jobs[0].t_index)
            ]["models"][model_name]
            for period in ("su22", "fa22", "sp23")
        }

        expected_counts: dict[str, np.ndarray] = {}
        for period, histogram in histograms.items():
            nominal = nominal_period_results[period]
            parameters = np.asarray(nominal["parameters"], dtype=float)
            centers = histogram["centers"]
            expected = (
                gaussian_signal_from_yield(
                    centers,
                    parameters[0],
                    parameters[1],
                    parameters[2],
                    histogram["histogram_bin_width_gev2"],
                )
                + evaluate_background(
                    model_name,
                    centers,
                    parameters,
                )
            )
            expected_counts[period] = np.maximum(expected, 0.0)
        # endfor

        replica_values = {
            period: {
                "mean": [],
                "sigma": [],
                "yield": [],
            }
            for period in ("su22", "fa22", "sp23")
        }
        successful = 0

        for _ in range(config.replicas):
            fluctuated = {
                period: rng.poisson(expected_counts[period]).astype(float)
                for period in ("su22", "fa22", "sp23")
            }
            replica = fit_joint_background_model(
                background_model=model_name,
                jobs=jobs,
                config=config,
                fixed_mean_gev2=fixed_mean,
                counts_override=fluctuated,
            )
            if not replica.get("success", False):
                continue
            # endif
            successful += 1
            for period in ("su22", "fa22", "sp23"):
                period_fit = replica["period_results"][period]
                replica_values[period]["mean"].append(
                    period_fit["mean_gev2"]
                )
                replica_values[period]["sigma"].append(
                    period_fit["sigma_gev2"]
                )
                replica_values[period]["yield"].append(
                    period_fit["signal_yield"]
                )
            # endfor
        # endfor

        for period in ("su22", "fa22", "sp23"):
            item = results_by_key[
                (period, jobs[0].x_index, jobs[0].t_index)
            ]
            summary: dict[str, Any] = {
                "requested": config.replicas,
                "successful": successful,
                "model_order_fixed": True,
                "background_model": model_name,
            }
            for quantity in ("mean", "sigma", "yield"):
                values = np.asarray(
                    replica_values[period][quantity],
                    dtype=float,
                )
                summary[quantity] = {
                    "mean": (
                        float(np.mean(values))
                        if values.size > 0
                        else math.nan
                    ),
                    "standard_deviation": (
                        float(np.std(values, ddof=1))
                        if values.size > 1
                        else math.nan
                    ),
                    "p16": (
                        float(np.percentile(values, 16.0))
                        if values.size > 0
                        else math.nan
                    ),
                    "median": (
                        float(np.percentile(values, 50.0))
                        if values.size > 0
                        else math.nan
                    ),
                    "p84": (
                        float(np.percentile(values, 84.0))
                        if values.size > 0
                        else math.nan
                    ),
                }
            # endfor
            item["replica_summary"] = summary
        # endfor

        print(
            f"Completed replica bin {group_index}/{len(corrected_groups)} "
            f"({successful}/{config.replicas} successful)",
            flush=True,
        )
    # endfor


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


def execute_fit_jobs(
    jobs: list[FitJob],
    workers: int,
    joint_config: JointFitConfig,
) -> list[dict[str, Any]]:
    """Fit before-correction jobs independently and corrected bins jointly."""
    actual_workers = max(1, min(int(workers), 7, os.cpu_count() or 1))
    print("Running channel_selection_mx2_fits_v30.py", flush=True)

    before_jobs = [job for job in jobs if job.stage == "before"]
    after_jobs = [job for job in jobs if job.stage == "after"]

    corrected_group_map: dict[tuple[int, int], list[FitJob]] = {}
    for job in after_jobs:
        corrected_group_map.setdefault(
            (job.x_index, job.t_index),
            [],
        ).append(job)
    # endfor
    corrected_groups = [
        sorted(group, key=lambda item: item.period)
        for _, group in sorted(corrected_group_map.items())
    ]

    print(
        f"Running {len(before_jobs)} independent before-correction jobs and "
        f"{len(corrected_groups)} three-period corrected joint jobs with "
        f"{actual_workers} worker(s).",
        flush=True,
    )

    before_results: list[dict[str, Any]] = []
    if actual_workers == 1:
        for index, job in enumerate(before_jobs, start=1):
            before_results.append(fit_one_job(job))
            print(
                f"Completed before-correction fit {index}/{len(before_jobs)}",
                flush=True,
            )
        # endfor
    else:
        with ProcessPoolExecutor(max_workers=actual_workers) as executor:
            future_map = {
                executor.submit(fit_one_job, job): job
                for job in before_jobs
            }
            completed = 0
            for future in as_completed(future_map):
                job = future_map[future]
                try:
                    before_results.append(future.result())
                except Exception as exc:
                    raise RuntimeError(
                        f"Before-correction worker failed for {job.period}, "
                        f"x index {job.x_index}, t index {job.t_index}."
                    ) from exc
                # endtry
                completed += 1
                print(
                    f"Completed before-correction fit "
                    f"{completed}/{len(before_jobs)}",
                    flush=True,
                )
            # endfor
        # endwith
    # endif

    corrected_results: list[dict[str, Any]] = []
    if actual_workers == 1:
        for index, group in enumerate(corrected_groups, start=1):
            corrected_results.extend(
                fit_one_joint_corrected_bin(group, joint_config)
            )
            print(
                f"Completed corrected joint fit "
                f"{index}/{len(corrected_groups)}",
                flush=True,
            )
        # endfor
    else:
        with ProcessPoolExecutor(max_workers=actual_workers) as executor:
            future_map = {
                executor.submit(
                    fit_one_joint_corrected_bin,
                    group,
                    joint_config,
                ): group
                for group in corrected_groups
            }
            completed = 0
            for future in as_completed(future_map):
                group = future_map[future]
                try:
                    corrected_results.extend(future.result())
                except Exception as exc:
                    raise RuntimeError(
                        "Corrected joint worker failed for "
                        f"x index {group[0].x_index}, "
                        f"t index {group[0].t_index}."
                    ) from exc
                # endtry
                completed += 1
                print(
                    f"Completed corrected joint fit "
                    f"{completed}/{len(corrected_groups)}",
                    flush=True,
                )
            # endfor
        # endwith
    # endif

    corrected_results = apply_neighbor_mean_fallbacks(
        corrected_results=corrected_results,
        corrected_groups=corrected_groups,
        config=joint_config,
    )
    run_joint_replicas(
        corrected_results=corrected_results,
        corrected_groups=corrected_groups,
        config=joint_config,
    )

    results = [*before_results, *corrected_results]
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
            "fit_quality_class": item.get(
                "fit_quality_class",
                item.get("status", ""),
            ),
            "shared_mean_fit": bool(
                item.get("joint_fit", {}).get("shared_mean_gev2") is not None
            ),
            "shared_mean_gev2": item.get(
                "joint_fit",
                {},
            ).get("shared_mean_gev2", math.nan),
            "shared_mean_error_gev2": item.get(
                "joint_fit",
                {},
            ).get("shared_mean_error_gev2", math.nan),
            "joint_chi2_ndf": item.get(
                "joint_fit",
                {},
            ).get("joint_chi2_ndf", math.nan),
            "joint_aicc": item.get(
                "joint_fit",
                {},
            ).get("joint_aicc", math.nan),
            "joint_aggregate_local_score": item.get(
                "joint_fit",
                {},
            ).get("aggregate_local_score", math.nan),
            "joint_maximum_period_local_score": item.get(
                "joint_fit",
                {},
            ).get("maximum_period_local_score", math.nan),
            "parameter_source": item.get(
                "joint_fit",
                {},
            ).get("parameter_source", "independent_fit"),
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
    lines.append("% Auto-generated by channel_selection_mx2_fits_v30.py")
    lines.append("% Corrected simultaneous three-period Gaussian fits with a shared mean and adaptively selected polynomial backgrounds.")
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
# Phase-1 carbon-template diagnostics
# =============================================================================

def fit_job_index(
    jobs: list[FitJob],
) -> dict[tuple[str, str, int, int], FitJob]:
    """Index fit jobs by period, stage, xB bin, and -tprime bin."""
    return {
        (job.period, job.stage, job.x_index, job.t_index): job
        for job in jobs
    }


def carbon_group_key(job: FitJob, grouping: str) -> tuple[Any, ...]:
    """Return the requested carbon-normalization group key."""
    if grouping == "tprime-column":
        return (job.period, job.stage, job.t_index)
    # endif
    if grouping == "xB-row":
        return (job.period, job.stage, job.x_index)
    # endif
    if grouping == "bin":
        return (job.period, job.stage, job.x_index, job.t_index)
    # endif
    if grouping == "period":
        return (job.period, job.stage)
    # endif
    raise ValueError(f"Unknown carbon grouping: {grouping}")


def carbon_group_label(key: tuple[Any, ...], grouping: str) -> str:
    """Build a human-readable carbon-normalization group label."""
    prefix = f"{PERIOD_LABELS[key[0]]}, {STAGE_LABELS[key[1]]}"
    if grouping == "tprime-column":
        low, high = MINUS_TPRIME_BINS_GEV2[int(key[2])]
        return f"{prefix}, {low:.2f} < -t' < {high:.2f} GeV^2"
    # endif
    if grouping == "xB-row":
        low, high = XB_BINS[int(key[2])]
        return f"{prefix}, {low:.2f} < xB < {high:.2f}"
    # endif
    if grouping == "bin":
        return (
            f"{prefix}, "
            f"{bin_identifier(int(key[2]), int(key[3]))}"
        )
    # endif
    return prefix


def count_in_interval(
    values: np.ndarray,
    minimum: float,
    maximum: float,
) -> int:
    """Count finite values in a half-open interval."""
    array = np.asarray(values, dtype=float)
    return int(
        np.count_nonzero(
            np.isfinite(array)
            & (array >= minimum)
            & (array < maximum)
        )
    )


def ratio_with_poisson_error(
    numerator: int,
    denominator: int,
) -> tuple[float, float]:
    """Calculate a count ratio and independent-Poisson uncertainty."""
    if denominator <= 0:
        return math.nan, math.nan
    # endif
    value = numerator / denominator
    if numerator <= 0:
        return value, math.nan
    # endif
    uncertainty = value * math.sqrt(
        1.0 / numerator + 1.0 / denominator
    )
    return value, uncertainty


def build_carbon_phase1_diagnostics(
    nh3_jobs: list[FitJob],
    carbon_jobs: list[FitJob],
    grouping: str,
    control_min: float,
    control_max: float,
    scan_boundaries: tuple[float, ...],
    hist_min: float,
    hist_max: float,
    hist_bins: int,
    validation_low: tuple[float, float],
    validation_high: tuple[float, float],
) -> dict[str, Any]:
    """Construct Phase-1 low-Mx2 carbon-normalization diagnostics."""
    nh3 = fit_job_index(nh3_jobs)
    carbon = fit_job_index(carbon_jobs)
    if set(nh3) != set(carbon):
        raise RuntimeError(
            "NH3 and carbon job keys differ; the two samples must use the "
            "same periods, correction stages, and kinematic binning."
        )
    # endif

    grouped: dict[tuple[Any, ...], list[tuple[FitJob, FitJob]]] = {}
    for key, nh3_job in nh3.items():
        grouped.setdefault(
            carbon_group_key(nh3_job, grouping),
            [],
        ).append((nh3_job, carbon[key]))
    # endfor

    group_rows: list[dict[str, Any]] = []
    scan_rows: list[dict[str, Any]] = []
    group_info: dict[tuple[Any, ...], dict[str, Any]] = {}

    for key, pairs in sorted(grouped.items()):
        nh3_control = sum(
            count_in_interval(a.values, control_min, control_max)
            for a, _ in pairs
        )
        carbon_control = sum(
            count_in_interval(b.values, control_min, control_max)
            for _, b in pairs
        )
        scale, scale_error = ratio_with_poisson_error(
            nh3_control,
            carbon_control,
        )
        row = {
            "grouping": grouping,
            "group_key": "|".join(str(item) for item in key),
            "group_label": carbon_group_label(key, grouping),
            "period": key[0],
            "stage": key[1],
            "x_index": (
                int(key[2]) if grouping == "xB-row" else math.nan
            ),
            "t_index": (
                int(key[2])
                if grouping == "tprime-column"
                else math.nan
            ),
            "control_min_gev2": control_min,
            "control_max_gev2": control_max,
            "nh3_control_count": nh3_control,
            "carbon_control_count": carbon_control,
            "carbon_scale": scale,
            "carbon_scale_uncertainty": scale_error,
            "relative_scale_uncertainty": (
                scale_error / scale
                if np.isfinite(scale_error)
                and np.isfinite(scale)
                and scale != 0.0
                else math.nan
            ),
            "number_of_bins_pooled": len(pairs),
        }
        group_rows.append(row)
        group_info[key] = row

        for boundary in sorted(set(scan_boundaries)):
            n_count = sum(
                count_in_interval(a.values, control_min, boundary)
                for a, _ in pairs
            )
            c_count = sum(
                count_in_interval(b.values, control_min, boundary)
                for _, b in pairs
            )
            value, error = ratio_with_poisson_error(n_count, c_count)
            scan_rows.append(
                {
                    "grouping": grouping,
                    "group_key": row["group_key"],
                    "group_label": row["group_label"],
                    "period": key[0],
                    "stage": key[1],
                    "control_min_gev2": control_min,
                    "control_max_gev2": boundary,
                    "nh3_control_count": n_count,
                    "carbon_control_count": c_count,
                    "carbon_scale": value,
                    "carbon_scale_uncertainty": error,
                }
            )
        # endfor
    # endfor

    edges = np.linspace(hist_min, hist_max, hist_bins + 1)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_rows: list[dict[str, Any]] = []
    payloads: dict[tuple[str, str, int, int], dict[str, Any]] = {}

    for key in sorted(nh3):
        nh3_job = nh3[key]
        carbon_job = carbon[key]
        group_key = carbon_group_key(nh3_job, grouping)
        info = group_info[group_key]
        scale = float(info["carbon_scale"])
        scale_error = float(info["carbon_scale_uncertainty"])

        nh3_counts, _ = np.histogram(nh3_job.values, bins=edges)
        carbon_counts, _ = np.histogram(carbon_job.values, bins=edges)
        nh3_counts = nh3_counts.astype(float)
        carbon_counts = carbon_counts.astype(float)

        scaled = scale * carbon_counts
        scaled_error = np.sqrt(
            np.maximum(scale**2 * carbon_counts, 0.0)
            + np.maximum((scale_error * carbon_counts) ** 2, 0.0)
        )
        subtraction = nh3_counts - scaled
        subtraction_error = np.sqrt(
            np.maximum(nh3_counts, 0.0) + scaled_error**2
        )

        ratio = np.full_like(nh3_counts, np.nan)
        ratio_error = np.full_like(nh3_counts, np.nan)
        valid = (nh3_counts > 0.0) & (carbon_counts > 0.0) & (scale > 0.0)
        ratio[valid] = nh3_counts[valid] / scaled[valid]
        ratio_error[valid] = ratio[valid] * np.sqrt(
            1.0 / nh3_counts[valid]
            + 1.0 / carbon_counts[valid]
            + (
                (scale_error / scale) ** 2
                if np.isfinite(scale_error)
                else 0.0
            )
        )

        local_n = count_in_interval(
            nh3_job.values,
            control_min,
            control_max,
        )
        local_c = count_in_interval(
            carbon_job.values,
            control_min,
            control_max,
        )
        local_scale, local_error = ratio_with_poisson_error(
            local_n,
            local_c,
        )

        def closure(low: float, high: float) -> dict[str, float]:
            mask = (centers >= low) & (centers < high)
            observed = float(np.sum(nh3_counts[mask]))
            expected = float(np.sum(scaled[mask]))
            uncertainty = math.sqrt(
                float(np.sum(nh3_counts[mask]))
                + float(np.sum(scaled_error[mask] ** 2))
            )
            return {
                "observed": observed,
                "scaled_carbon": expected,
                "ratio": (
                    observed / expected
                    if expected > 0.0
                    else math.nan
                ),
                "pull": (
                    (observed - expected) / uncertainty
                    if uncertainty > 0.0
                    else math.nan
                ),
            }
        # enddef

        low_closure = closure(*validation_low)
        high_closure = closure(*validation_high)

        record = {
            "period": nh3_job.period,
            "stage": nh3_job.stage,
            "x_index": nh3_job.x_index,
            "t_index": nh3_job.t_index,
            "bin_number": combined_bin_number(
                nh3_job.x_index,
                nh3_job.t_index,
            ),
            "bin_id": bin_identifier(
                nh3_job.x_index,
                nh3_job.t_index,
            ),
            "xB_min": nh3_job.x_min,
            "xB_max": nh3_job.x_max,
            "minus_tprime_min_gev2": nh3_job.minus_tprime_min_gev2,
            "minus_tprime_max_gev2": nh3_job.minus_tprime_max_gev2,
            "grouping": grouping,
            "group_key": info["group_key"],
            "control_min_gev2": control_min,
            "control_max_gev2": control_max,
            "group_nh3_control_count": info["nh3_control_count"],
            "group_carbon_control_count": info["carbon_control_count"],
            "group_carbon_scale": scale,
            "group_carbon_scale_uncertainty": scale_error,
            "bin_nh3_control_count": local_n,
            "bin_carbon_control_count": local_c,
            "bin_local_carbon_scale": local_scale,
            "bin_local_carbon_scale_uncertainty": local_error,
            "local_over_group_scale": (
                local_scale / scale
                if np.isfinite(local_scale)
                and np.isfinite(scale)
                and scale != 0.0
                else math.nan
            ),
            "low_validation_min_gev2": validation_low[0],
            "low_validation_max_gev2": validation_low[1],
            "low_validation_ratio": low_closure["ratio"],
            "low_validation_pull": low_closure["pull"],
            "high_validation_min_gev2": validation_high[0],
            "high_validation_max_gev2": validation_high[1],
            "high_validation_ratio": high_closure["ratio"],
            "high_validation_pull": high_closure["pull"],
            "negative_subtracted_bins": int(
                np.count_nonzero(subtraction < 0.0)
            ),
            "subtracted_histogram_integral": float(
                np.sum(subtraction)
            ),
        }
        bin_rows.append(record)
        payloads[key] = {
            "record": record,
            "centers": centers,
            "edges": edges,
            "nh3": nh3_counts,
            "nh3_error": np.sqrt(np.maximum(nh3_counts, 0.0)),
            "scaled_carbon": scaled,
            "scaled_carbon_error": scaled_error,
            "subtraction": subtraction,
            "subtraction_error": subtraction_error,
            "ratio": ratio,
            "ratio_error": ratio_error,
        }
    # endfor

    return {
        "configuration": {
            "phase": 1,
            "carbon_used_in_production_fit": False,
            "production_fit_changed": False,
            "grouping": grouping,
            "control_region_gev2": [control_min, control_max],
            "boundary_scan_gev2": list(sorted(set(scan_boundaries))),
            "validation_low_region_gev2": list(validation_low),
            "validation_high_region_gev2": list(validation_high),
            "normalization_is_count_based": True,
            "normalization_is_not_fitted": True,
            "diagnostic_subtraction_only": True,
            "summary_products": [
                "nominal normalization versus grouping coordinate",
                "normalization versus control-region endpoint",
                "local normalization divided by pooled normalization",
                "negative-bin fraction after diagnostic subtraction",
                "low- and high-side closure pull distributions",
                "focused corrected subtraction for bins 1, 7, 13, and 19",
            ],
        },
        "group_rows": group_rows,
        "scan_rows": scan_rows,
        "bin_rows": bin_rows,
        "payloads": payloads,
    }


def draw_carbon_phase1_canvas(
    diagnostics: dict[str, Any],
    period: str,
    stage: str,
    mode: str,
    output_path: Path,
) -> None:
    """Draw one 4x6 Phase-1 carbon diagnostic canvas."""
    payloads = diagnostics["payloads"]
    fig, axes = plt.subplots(
        len(XB_BINS),
        len(MINUS_TPRIME_BINS_GEV2),
        figsize=(24, 15),
        sharex=True,
        sharey=(mode == "ratio"),
    )

    for x_index in range(len(XB_BINS)):
        for t_index in range(len(MINUS_TPRIME_BINS_GEV2)):
            ax = axes[x_index, t_index]
            payload = payloads[(period, stage, x_index, t_index)]
            record = payload["record"]
            x = payload["centers"]

            if mode == "overlay":
                ax.errorbar(
                    x,
                    payload["nh3"],
                    yerr=payload["nh3_error"],
                    fmt="o",
                    markersize=2.0,
                    linewidth=0.7,
                    label="NH$_3$",
                )
                ax.step(
                    x,
                    payload["scaled_carbon"],
                    where="mid",
                    linewidth=1.3,
                    label="Scaled carbon",
                )
                ax.fill_between(
                    x,
                    np.maximum(
                        payload["scaled_carbon"]
                        - payload["scaled_carbon_error"],
                        0.0,
                    ),
                    payload["scaled_carbon"]
                    + payload["scaled_carbon_error"],
                    step="mid",
                    alpha=0.20,
                )
                title_extra = (
                    f"$\\alpha_C$={record['group_carbon_scale']:.3g}"
                    f"$\\pm${record['group_carbon_scale_uncertainty']:.2g}"
                )
            elif mode == "subtraction":
                ax.errorbar(
                    x,
                    payload["subtraction"],
                    yerr=payload["subtraction_error"],
                    fmt="o",
                    markersize=2.0,
                    linewidth=0.7,
                )
                ax.axhline(0.0, linewidth=0.8)
                title_extra = (
                    f"{record['negative_subtracted_bins']} negative bins"
                )
            elif mode == "ratio":
                valid = (
                    np.isfinite(payload["ratio"])
                    & np.isfinite(payload["ratio_error"])
                )
                ax.errorbar(
                    x[valid],
                    payload["ratio"][valid],
                    yerr=payload["ratio_error"][valid],
                    fmt="o",
                    markersize=2.0,
                    linewidth=0.7,
                )
                ax.axhline(1.0, linestyle="--", linewidth=0.9)
                ax.set_ylim(0.0, 5.0)
                title_extra = (
                    f"L={record['low_validation_ratio']:.2g}, "
                    f"H={record['high_validation_ratio']:.2g}"
                )
            else:
                raise ValueError(mode)
            # endif

            ax.axvspan(
                record["control_min_gev2"],
                record["control_max_gev2"],
                alpha=0.12,
            )
            ax.axvline(
                NEUTRON_MASS2_GEV2,
                linestyle="--",
                linewidth=0.8,
            )
            ax.set_title(
                f"Bin {record['bin_number']}: {title_extra}",
                fontsize=9,
            )
            ax.grid(alpha=0.22)
            if x_index == len(XB_BINS) - 1:
                ax.set_xlabel("$M_x^2$ (GeV$^2$)")
            # endif
            if t_index == 0:
                x_low, x_high = XB_BINS[x_index]
                ylabel = {
                    "overlay": "Counts",
                    "subtraction": "NH$_3$ - scaled C",
                    "ratio": "NH$_3$ / scaled C",
                }[mode]
                ax.set_ylabel(
                    f"{x_low:.2f} < $x_B$ < {x_high:.2f}\n{ylabel}"
                )
            # endif
        # endfor
    # endfor

    if mode == "overlay":
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.995),
            ncol=2,
            frameon=True,
        )
    # endif

    subtitle = {
        "overlay": "NH$_3$ and low-$M_x^2$ normalized carbon",
        "subtraction": (
            "Diagnostic NH$_3$ minus scaled carbon "
            "(not used in production fits)"
        ),
        "ratio": "Carbon-template closure ratio",
    }[mode]
    fig.suptitle(
        f"{PERIOD_LABELS[period]}: {STAGE_LABELS[stage]}\n"
        f"Phase-1 {subtitle}",
        fontsize=15,
        y=0.945 if mode == "overlay" else 0.985,
    )
    fig.tight_layout(
        rect=(0.0, 0.0, 1.0, 0.895 if mode == "overlay" else 0.945)
    )
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_carbon_scale_scans(
    diagnostics: dict[str, Any],
    output_dir: Path,
) -> list[Path]:
    """Plot normalization scale versus the low-Mx2 upper boundary."""
    frame = pd.DataFrame(diagnostics["scan_rows"])
    outputs: list[Path] = []
    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            subset = frame[
                (frame["period"] == period)
                & (frame["stage"] == stage)
            ]
            if subset.empty:
                continue
            # endif
            fig, ax = plt.subplots(figsize=(11, 7))
            for _, group in subset.groupby("group_key", sort=True):
                group = group.sort_values("control_max_gev2")
                ax.errorbar(
                    group["control_max_gev2"],
                    group["carbon_scale"],
                    yerr=group["carbon_scale_uncertainty"],
                    marker="o",
                    linewidth=1.1,
                    markersize=4,
                    label=group["group_label"].iloc[0],
                )
            # endfor
            ax.set_xlabel(
                "Low-$M_x^2$ control-region upper edge (GeV$^2$)"
            )
            ax.set_ylabel("NH$_3$ / carbon normalization")
            ax.set_title(
                f"{PERIOD_LABELS[period]}: {STAGE_LABELS[stage]}\n"
                "Phase-1 carbon-normalization stability"
            )
            ax.grid(alpha=0.25)
            ax.legend(fontsize=7, ncol=2)
            fig.tight_layout()
            path = output_dir / (
                f"carbon_phase1_scale_scan_{period}_{stage}_v27.png"
            )
            fig.savefig(path, dpi=180)
            plt.close(fig)
            outputs.append(path)
        # endfor
    # endfor
    return outputs



def plot_carbon_group_scales(
    diagnostics: dict[str, Any],
    output_dir: Path,
) -> list[Path]:
    """Plot nominal carbon normalizations versus the grouping coordinate."""
    frame = pd.DataFrame(diagnostics["group_rows"])
    outputs: list[Path] = []
    grouping = diagnostics["configuration"]["grouping"]

    for stage in ("before", "after"):
        subset = frame[frame["stage"] == stage].copy()
        if subset.empty:
            continue
        # endif

        fig, ax = plt.subplots(figsize=(11, 7))
        for period in ("su22", "fa22", "sp23"):
            period_frame = subset[subset["period"] == period].copy()
            if period_frame.empty:
                continue
            # endif
            if grouping == "tprime-column":
                period_frame = period_frame.sort_values("t_index")
                x_values = period_frame["t_index"].to_numpy(dtype=float) + 1.0
            elif grouping == "xB-row":
                period_frame = period_frame.sort_values("x_index")
                x_values = period_frame["x_index"].to_numpy(dtype=float) + 1.0
            else:
                period_frame = period_frame.reset_index(drop=True)
                x_values = np.arange(1, len(period_frame) + 1, dtype=float)
            # endif
            ax.errorbar(
                x_values,
                period_frame["carbon_scale"],
                yerr=period_frame["carbon_scale_uncertainty"],
                marker="o",
                linewidth=1.2,
                markersize=5,
                label=PERIOD_LABELS[period],
            )
        # endfor

        if grouping == "tprime-column":
            ax.set_xlabel("$-t'$ column")
            ax.set_xticks(np.arange(1, len(MINUS_TPRIME_BINS_GEV2) + 1))
            labels = [
                f"{low:.2f}–{high:.2f}"
                for low, high in MINUS_TPRIME_BINS_GEV2
            ]
            ax.set_xticklabels(labels, rotation=25, ha="right")
        elif grouping == "xB-row":
            ax.set_xlabel("$x_B$ row")
            ax.set_xticks(np.arange(1, len(XB_BINS) + 1))
            labels = [f"{low:.2f}–{high:.2f}" for low, high in XB_BINS]
            ax.set_xticklabels(labels)
        else:
            ax.set_xlabel("Normalization group index")
        # endif

        ax.set_ylabel("NH$_3$ / carbon normalization")
        ax.set_title(
            f"{STAGE_LABELS[stage]}\n"
            "Nominal Phase-1 carbon normalization"
        )
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        path = output_dir / f"carbon_phase1_group_scales_{stage}_v27.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)
    # endfor
    return outputs


def plot_carbon_local_scale_closure(
    diagnostics: dict[str, Any],
    output_dir: Path,
) -> list[Path]:
    """Plot each bin's local control-region scale divided by its pooled scale."""
    frame = pd.DataFrame(diagnostics["bin_rows"])
    outputs: list[Path] = []

    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            subset = frame[
                (frame["period"] == period)
                & (frame["stage"] == stage)
            ].sort_values("bin_number")
            if subset.empty:
                continue
            # endif
            fig, ax = plt.subplots(figsize=(12, 7))
            for t_index in range(len(MINUS_TPRIME_BINS_GEV2)):
                column = subset[subset["t_index"] == t_index]
                ax.plot(
                    column["bin_number"],
                    column["local_over_group_scale"],
                    marker="o",
                    linewidth=1.0,
                    markersize=4,
                    label=f"$-t'$ column {t_index + 1}",
                )
            # endfor
            ax.axhline(1.0, linestyle="--", linewidth=1.0)
            for boundary in (6.5, 12.5, 18.5):
                ax.axvline(boundary, linewidth=0.7, alpha=0.4)
            # endfor
            ax.set_xlabel("Kinematic bin number")
            ax.set_ylabel("Local scale / pooled-group scale")
            ax.set_title(
                f"{PERIOD_LABELS[period]}: {STAGE_LABELS[stage]}\n"
                "Low-$M_x^2$ normalization closure by kinematic bin"
            )
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8, ncol=2)
            fig.tight_layout()
            path = output_dir / (
                f"carbon_phase1_local_over_group_"
                f"{period}_{stage}_v27.png"
            )
            fig.savefig(path, dpi=180)
            plt.close(fig)
            outputs.append(path)
        # endfor
    # endfor
    return outputs


def plot_carbon_negative_bin_fraction(
    diagnostics: dict[str, Any],
    output_dir: Path,
) -> list[Path]:
    """Plot the fraction of histogram bins negative after diagnostic subtraction."""
    frame = pd.DataFrame(diagnostics["bin_rows"])
    histogram_bins = int(
        len(next(iter(diagnostics["payloads"].values()))["centers"])
    )
    outputs: list[Path] = []

    for stage in ("before", "after"):
        fig, ax = plt.subplots(figsize=(12, 7))
        for period in ("su22", "fa22", "sp23"):
            subset = frame[
                (frame["period"] == period)
                & (frame["stage"] == stage)
            ].sort_values("bin_number")
            if subset.empty:
                continue
            # endif
            fraction = subset["negative_subtracted_bins"] / histogram_bins
            ax.plot(
                subset["bin_number"],
                fraction,
                marker="o",
                linewidth=1.1,
                markersize=4,
                label=PERIOD_LABELS[period],
            )
        # endfor
        for boundary in (6.5, 12.5, 18.5):
            ax.axvline(boundary, linewidth=0.7, alpha=0.4)
        # endfor
        ax.set_xlabel("Kinematic bin number")
        ax.set_ylabel("Fraction of negative subtracted histogram bins")
        ax.set_ylim(bottom=0.0)
        ax.set_title(
            f"{STAGE_LABELS[stage]}\n"
            "Diagnostic NH$_3$ minus scaled-carbon negativity"
        )
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        path = output_dir / (
            f"carbon_phase1_negative_bin_fraction_{stage}_v27.png"
        )
        fig.savefig(path, dpi=180)
        plt.close(fig)
        outputs.append(path)
    # endfor
    return outputs


def plot_carbon_closure_pulls(
    diagnostics: dict[str, Any],
    output_dir: Path,
) -> list[Path]:
    """Plot low- and high-side carbon-template closure pulls."""
    frame = pd.DataFrame(diagnostics["bin_rows"])
    outputs: list[Path] = []

    for stage in ("before", "after"):
        subset = frame[frame["stage"] == stage]
        if subset.empty:
            continue
        # endif
        for region_name, column, label in (
            (
                "low",
                "low_validation_pull",
                "Low-side validation-region pull",
            ),
            (
                "high",
                "high_validation_pull",
                "High-side validation-region pull",
            ),
        ):
            fig, ax = plt.subplots(figsize=(10, 7))
            finite_all: list[float] = []
            for period in ("su22", "fa22", "sp23"):
                values = subset.loc[
                    subset["period"] == period,
                    column,
                ].to_numpy(dtype=float)
                values = values[np.isfinite(values)]
                finite_all.extend(values.tolist())
                if values.size > 0:
                    ax.hist(
                        values,
                        bins=np.linspace(-8.0, 8.0, 33),
                        histtype="step",
                        linewidth=1.5,
                        label=PERIOD_LABELS[period],
                    )
                # endif
            # endfor
            ax.axvline(0.0, linestyle="--", linewidth=1.0)
            ax.axvline(-3.0, linestyle=":", linewidth=0.9)
            ax.axvline(3.0, linestyle=":", linewidth=0.9)
            ax.set_xlabel(label)
            ax.set_ylabel("Kinematic bins")
            title = f"{STAGE_LABELS[stage]}: carbon-template closure pulls"
            if finite_all:
                values_all = np.asarray(finite_all, dtype=float)
                title += (
                    f"\nmean={np.mean(values_all):.2f}, "
                    f"RMS={np.std(values_all):.2f}"
                )
            # endif
            ax.set_title(title)
            ax.grid(alpha=0.25)
            ax.legend()
            fig.tight_layout()
            path = output_dir / (
                f"carbon_phase1_{region_name}_closure_pulls_"
                f"{stage}_v27.png"
            )
            fig.savefig(path, dpi=180)
            plt.close(fig)
            outputs.append(path)
        # endfor
    # endfor
    return outputs


def plot_problematic_carbon_subtractions(
    diagnostics: dict[str, Any],
    output_dir: Path,
) -> Path:
    """Draw focused corrected-data subtractions for bins 1, 7, 13, and 19."""
    problematic = ((0, 0), (1, 0), (2, 0), (3, 0))
    fig, axes = plt.subplots(
        3,
        4,
        figsize=(18, 12),
        sharex=True,
        sharey=False,
    )
    for row, period in enumerate(("su22", "fa22", "sp23")):
        for column, (x_index, t_index) in enumerate(problematic):
            ax = axes[row, column]
            payload = diagnostics["payloads"][
                (period, "after", x_index, t_index)
            ]
            record = payload["record"]
            ax.errorbar(
                payload["centers"],
                payload["subtraction"],
                yerr=payload["subtraction_error"],
                fmt="o",
                markersize=2.5,
                linewidth=0.8,
            )
            ax.axhline(0.0, linewidth=0.8)
            ax.axvline(
                NEUTRON_MASS2_GEV2,
                linestyle="--",
                linewidth=0.9,
            )
            ax.axvspan(
                record["control_min_gev2"],
                record["control_max_gev2"],
                alpha=0.12,
            )
            ax.grid(alpha=0.25)
            if row == 0:
                ax.set_title(f"Bin {record['bin_number']}")
            # endif
            if column == 0:
                ax.set_ylabel(
                    f"{PERIOD_LABELS[period]}\nNH$_3$ - scaled C"
                )
            # endif
            if row == 2:
                ax.set_xlabel("$M_x^2$ (GeV$^2$)")
            # endif
        # endfor
    # endfor
    fig.suptitle(
        "Corrected-data Phase-1 carbon subtraction in the four "
        "problematic lowest-$-t'$ bins\n"
        "Diagnostic only; not used in production fits",
        fontsize=15,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.955))
    path = output_dir / "carbon_phase1_problematic_bins_after_v27.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def write_carbon_phase1_products(
    nominal_diagnostics: dict[str, Any],
    column_diagnostics: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    """Write streamlined and organized Phase-1 carbon products."""
    summary_dir = output_dir / "summary"
    overlays_dir = output_dir / "overlays"
    subtraction_dir = output_dir / "subtractions"
    ratios_dir = output_dir / "ratios"
    tables_dir = output_dir / "tables"
    for directory in (
        output_dir,
        summary_dir,
        overlays_dir,
        subtraction_dir,
        ratios_dir,
        tables_dir,
    ):
        ensure_directory(directory)
    # endfor

    nominal_group_csv = (
        tables_dir / "carbon_phase1_period_scales_v27.csv"
    )
    column_group_csv = (
        tables_dir / "carbon_phase1_tprime_column_scales_v27.csv"
    )
    nominal_scan_csv = (
        tables_dir / "carbon_phase1_period_scale_scan_v27.csv"
    )
    nominal_bin_csv = (
        tables_dir / "carbon_phase1_bin_diagnostics_v27.csv"
    )
    json_path = tables_dir / "carbon_phase1_diagnostics_v27.json"

    pd.DataFrame(nominal_diagnostics["group_rows"]).to_csv(
        nominal_group_csv,
        index=False,
    )
    pd.DataFrame(column_diagnostics["group_rows"]).to_csv(
        column_group_csv,
        index=False,
    )
    pd.DataFrame(nominal_diagnostics["scan_rows"]).to_csv(
        nominal_scan_csv,
        index=False,
    )
    pd.DataFrame(nominal_diagnostics["bin_rows"]).to_csv(
        nominal_bin_csv,
        index=False,
    )

    json_payload = {
        "nominal_period_wide": {
            "configuration": nominal_diagnostics["configuration"],
            "group_scales": nominal_diagnostics["group_rows"],
            "scale_scan": nominal_diagnostics["scan_rows"],
            "bin_diagnostics": nominal_diagnostics["bin_rows"],
        },
        "tprime_column_cross_check": {
            "configuration": column_diagnostics["configuration"],
            "group_scales": column_diagnostics["group_rows"],
        },
    }
    json_path.write_text(
        json.dumps(json_safe(json_payload), indent=2) + "\n",
        encoding="utf-8",
    )

    plots: list[Path] = []
    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            destinations = {
                "overlay": overlays_dir,
                "subtraction": subtraction_dir,
                "ratio": ratios_dir,
            }
            for mode, directory in destinations.items():
                path = directory / (
                    f"carbon_phase1_{mode}_{period}_{stage}_v27.png"
                )
                draw_carbon_phase1_canvas(
                    nominal_diagnostics,
                    period,
                    stage,
                    mode,
                    path,
                )
                plots.append(path)
            # endfor
        # endfor
    # endfor

    # Retain the important -tprime stability demonstration as a cross-check.
    plots.extend(
        plot_carbon_group_scales(column_diagnostics, summary_dir)
    )
    plots.extend(
        plot_carbon_scale_scans(nominal_diagnostics, summary_dir)
    )
    plots.append(
        plot_problematic_carbon_subtractions(
            nominal_diagnostics,
            summary_dir,
        )
    )

    return {
        "period_scale_csv": nominal_group_csv,
        "tprime_scale_csv": column_group_csv,
        "scan_csv": nominal_scan_csv,
        "bin_csv": nominal_bin_csv,
        "json": json_path,
        "plots": plots,
        "summary_dir": summary_dir,
    }




# =============================================================================
# Carbon-template display helpers
# =============================================================================

def smooth_template_for_display(
    values: np.ndarray,
    sigma_bins: float = 1.25,
) -> np.ndarray:
    """
    Smooth a finite-statistics template for plotting only.

    The statistical fit always uses the original unsmoothed histogram bins.
    This helper is used solely to draw a visually continuous total-model curve.
    """
    array = np.asarray(values, dtype=float)
    if array.size < 3 or sigma_bins <= 0.0:
        return array.copy()
    # endif
    smoothed = gaussian_filter1d(
        array,
        sigma=sigma_bins,
        mode="nearest",
        truncate=3.0,
    )
    original_sum = float(np.sum(array))
    smoothed_sum = float(np.sum(smoothed))
    if smoothed_sum > 0.0 and np.isfinite(original_sum):
        smoothed *= original_sum / smoothed_sum
    # endif
    return smoothed


# =============================================================================
# Candidate carbon-template fits
# =============================================================================

def carbon_residual_background(
    x: np.ndarray,
    coefficients: np.ndarray,
    threshold_gev2: float = CARBON_RESIDUAL_THRESHOLD_GEV2,
) -> np.ndarray:
    """
    Evaluate a residual background that turns on continuously at the threshold.

    For Mx2 below threshold_gev2 the residual is identically zero. Above the
    threshold it is expanded in positive powers of

        dx = Mx2 - threshold_gev2,

    with no constant term. Therefore the residual approaches zero from the
    right and is continuous at the threshold:

        P1: a1 * dx
        P2: a1 * dx + a2 * dx^2
    """
    x_values = np.asarray(x, dtype=float)
    dx = np.maximum(x_values - float(threshold_gev2), 0.0)
    answer = np.zeros_like(dx, dtype=float)
    for power, coefficient in enumerate(coefficients, start=1):
        answer += float(coefficient) * dx**power
    # endfor
    return answer


def carbon_fit_aicc(chi2: float, n_points: int, n_parameters: int) -> float:
    """Return the small-sample corrected Akaike information criterion."""
    aic = float(chi2) + 2.0 * float(n_parameters)
    denominator = n_points - n_parameters - 1
    if denominator <= 0:
        return math.inf
    # endif
    return aic + (
        2.0 * n_parameters * (n_parameters + 1.0) / denominator
    )


def fit_one_carbon_candidate(
    nh3_jobs: list[FitJob],
    carbon_jobs: list[FitJob],
    scale_rows: dict[tuple[str, str], dict[str, float]],
    stage: str,
    x_index: int,
    t_index: int,
    fit_min_gev2: float,
    fit_max_gev2: float,
    residual_order: int,
    alpha_constraint_scale: float = 1.0,
) -> dict[str, Any]:
    """
    Fit one bin with shared mu, constrained alphas, and a thresholded residual.

    The residual P1/P2 contribution is exactly zero below
    CARBON_RESIDUAL_THRESHOLD_GEV2 and is expanded in positive powers of
    (Mx2-threshold) above it, with no constant term. It is therefore
    continuous and approaches zero from the right.
    """
    if residual_order not in (1, 2):
        raise ValueError('residual_order must be 1 or 2')
    # endif
    nh3_lookup = fit_job_index(nh3_jobs)
    carbon_lookup = fit_job_index(carbon_jobs)
    periods = ('su22', 'fa22', 'sp23')
    prepared: dict[str, dict[str, Any]] = {}
    for period in periods:
        key=(period,stage,x_index,t_index)
        nj=nh3_lookup[key]; cj=carbon_lookup[key]
        edges=np.linspace(nj.histogram_min_gev2,nj.histogram_max_gev2,nj.histogram_bins+1)
        centers=0.5*(edges[:-1]+edges[1:]); width=float(edges[1]-edges[0])
        n,_=np.histogram(nj.values,bins=edges); c,_=np.histogram(cj.values,bins=edges)
        n=n.astype(float); c=c.astype(float)
        mask=(centers>=fit_min_gev2)&(centers<=fit_max_gev2)
        row=scale_rows[(period,stage)]
        alpha0=float(row['carbon_scale']); alpha_err=float(row['carbon_scale_uncertainty'])*alpha_constraint_scale
        prepared[period]={
            'centers':centers,'width':width,'nh3':n,'carbon':c,'fit_mask':mask,
            'alpha0':alpha0,'alpha_error':max(alpha_err,1.0e-12),
        }
    # endfor

    ncoef=residual_order
    # shared mu; then alpha, yield, sigma, thresholded residual coefficients
    # for each period. There is no residual constant term.
    initial=[NEUTRON_MASS2_GEV2]; lower=[0.82]; upper=[0.95]
    for period in periods:
        item=prepared[period]; mask=item['fit_mask']
        residual=np.maximum(item['nh3'][mask]-item['alpha0']*item['carbon'][mask],0.0)
        peak=(item['centers'][mask]>=0.72)&(item['centers'][mask]<=1.04)
        y0=max(float(np.sum(residual[peak])),10.0)
        side=~peak
        # The thresholded residual has no constant term. Start the linear
        # coefficient from the above-threshold sideband scale and initialize
        # any quadratic coefficient to zero.
        above_threshold = (
            item['centers'][mask] >= CARBON_RESIDUAL_THRESHOLD_GEV2
        )
        if np.any(side & above_threshold):
            side_x = (
                item['centers'][mask][side & above_threshold]
                - CARBON_RESIDUAL_THRESHOLD_GEV2
            )
            side_y = residual[side & above_threshold]
            denominator = float(np.sum(side_x**2))
            slope0 = (
                float(np.sum(side_x * side_y) / denominator)
                if denominator > 0.0
                else 0.0
            )
        else:
            slope0 = 0.0
        # endif
        residual_initial = [slope0] + [0.0] * (residual_order - 1)
        initial.extend(
            [item['alpha0'], y0, 0.055] + residual_initial
        )
        lower.extend([0.0, 0.0, 0.015] + [-1.0e7] * ncoef)
        upper.extend(
            [max(10.0 * item['alpha0'], 1.0), 1.0e10, 0.20]
            + [1.0e7] * ncoef
        )
    # endfor

    stride=3+ncoef
    def unpack(p: np.ndarray, idx: int) -> tuple[float,float,float,np.ndarray]:
        off=1+idx*stride
        return float(p[off]),float(p[off+1]),float(p[off+2]),np.asarray(p[off+3:off+3+ncoef],dtype=float)
    # enddef

    def residual_function(p: np.ndarray) -> np.ndarray:
        mu=float(p[0]); blocks=[]
        for idx,period in enumerate(periods):
            item=prepared[period]; mask=item['fit_mask']
            alpha,yield_,sigma,coeff=unpack(p,idx)
            signal=gaussian_signal_from_yield(item['centers'],yield_,mu,sigma,item['width'])
            background=carbon_residual_background(item['centers'],coeff)
            prediction=alpha*item['carbon']+signal+background
            # NH3 and finite-carbon counting variance. Alpha is handled by its correlated prior.
            variance=np.maximum(item['nh3']+alpha**2*item['carbon'],1.0)
            blocks.append((item['nh3'][mask]-prediction[mask])/np.sqrt(variance[mask]))
            blocks.append(np.asarray([(alpha-item['alpha0'])/item['alpha_error']]))
        # endfor
        return np.concatenate(blocks)
    # enddef

    opt=least_squares(residual_function,np.asarray(initial),bounds=(np.asarray(lower),np.asarray(upper)),max_nfev=50000,xtol=1e-11,ftol=1e-11,gtol=1e-11)
    p=np.asarray(opt.x,dtype=float); pulls=residual_function(p)
    chi2=float(np.sum(pulls**2)); npts=int(pulls.size); npar=int(p.size); ndf=max(npts-npar,1)
    aicc=carbon_fit_aicc(chi2,npts,npar)
    cov=np.full((npar,npar),np.nan)
    if opt.jac.size and npts>npar:
        try: cov=np.linalg.pinv(opt.jac.T@opt.jac)*(chi2/ndf)
        except np.linalg.LinAlgError: pass
    # endif
    mu=float(p[0]); muerr=math.sqrt(max(float(cov[0,0]),0.0)) if np.isfinite(cov[0,0]) else math.nan
    records=[]; payloads={}; local_scores=[]
    for idx,period in enumerate(periods):
        item=prepared[period]; alpha,yield_,sigma,coeff=unpack(p,idx)
        signal=gaussian_signal_from_yield(item['centers'],yield_,mu,sigma,item['width'])
        residual_bg=carbon_residual_background(item['centers'],coeff)
        carbon_component=alpha*item['carbon']; total=carbon_component+signal+residual_bg
        variance=np.maximum(item['nh3']+alpha**2*item['carbon'],1.0)
        mask=item['fit_mask']; local=(item['nh3'][mask]-total[mask])/np.sqrt(variance[mask])
        local_chi2=float(np.sum(local**2)); local_ndf=max(int(np.count_nonzero(mask))-(3+ncoef),1)
        local_score=local_chi2/local_ndf; local_scores.append(local_score)
        off=1+idx*stride
        diag=[cov[off+j,off+j] for j in range(stride)]
        errs=[math.sqrt(max(float(v),0.0)) if np.isfinite(v) else math.nan for v in diag]
        significance=yield_/errs[1] if np.isfinite(errs[1]) and errs[1]>0 else math.nan
        rec={
            'period':period,'stage':stage,'x_index':x_index,'t_index':t_index,
            'bin_number':combined_bin_number(x_index,t_index),'bin_id':bin_identifier(x_index,t_index),
            'success':bool(opt.success),'message':str(opt.message),'residual_order':residual_order,
            'model':(
                f'Gaussian + constrained carbon + continuous-threshold '
                f'P{residual_order}'
            ),
            'shared_mean_gev2':mu,'shared_mean_error_gev2':muerr,
            'alpha_fit':alpha,'alpha_fit_error':errs[0],'alpha_control':item['alpha0'],
            'alpha_control_error':item['alpha_error'],'alpha_pull':(alpha-item['alpha0'])/item['alpha_error'],
            'signal_yield':yield_,'signal_yield_error':errs[1],'signal_significance':significance,
            'sigma_gev2':sigma,'sigma_error_gev2':errs[2],
            'residual_coefficients':coeff.tolist(),'residual_coefficient_errors':errs[3:],
            'local_chi2_ndf':local_score,'joint_chi2':chi2,'joint_ndf':ndf,'joint_chi2_ndf':chi2/ndf,
            'aicc':aicc,'fit_min_gev2':fit_min_gev2,'fit_max_gev2':fit_max_gev2,
            'residual_threshold_gev2':CARBON_RESIDUAL_THRESHOLD_GEV2,
            'residual_is_zero_below_threshold':True,
            'residual_is_continuous_at_threshold':True,
            'residual_has_constant_term':False,
            'production_fit_changed':False,
        }
        records.append(rec)
        payloads[period]={'record':rec,'centers':item['centers'],'nh3':item['nh3'],'error':np.sqrt(variance),
                          'carbon':carbon_component,'signal':signal,'residual':residual_bg,'total':total,'fit_mask':mask}
    # endfor
    return {'stage':stage,'x_index':x_index,'t_index':t_index,'success':bool(opt.success),
            'residual_order':residual_order,'shared_mean_gev2':mu,'shared_mean_error_gev2':muerr,
            'joint_chi2':chi2,'joint_ndf':ndf,'joint_chi2_ndf':chi2/ndf,'aicc':aicc,
            'maximum_local_chi2_ndf':max(local_scores),'period_records':records,'payloads':payloads}


def classify_carbon_candidate(selected: dict[str, Any], alternative: dict[str, Any]) -> str:
    """Classify peak identifiability and model stability."""
    records=selected['period_records']
    min_sig=min((r['signal_significance'] for r in records if np.isfinite(r['signal_significance'])),default=-math.inf)
    mean=selected['shared_mean_gev2']; meanerr=selected['shared_mean_error_gev2']
    delta_model=abs(mean-alternative['shared_mean_gev2'])
    max_local=selected['maximum_local_chi2_ndf']
    widths=[r['sigma_gev2'] for r in records]
    if (selected['success'] and min_sig>=5.0 and max_local<=2.0 and 0.84<=mean<=0.925
            and all(0.02<w<0.14 for w in widths) and delta_model<=max(0.008,2.0*meanerr)):
        return 'resolved'
    # endif
    if (selected['success'] and min_sig>=3.0 and max_local<=3.0 and 0.82<=mean<=0.95
            and all(0.015<=w<=0.18 for w in widths) and delta_model<=0.020):
        return 'marginal'
    # endif
    return 'unresolved'


def run_carbon_template_fits(
    nh3_jobs: list[FitJob],
    carbon_jobs: list[FitJob],
    period_diagnostics: dict[str, Any],
    fit_min_gev2: float,
    fit_max_gev2: float,
    control_min_gev2: float,
    control_endpoints_gev2: tuple[float, ...],
) -> list[dict[str, Any]]:
    """Run P1/P2 candidates plus fit-range and control-endpoint stability checks."""
    nominal_scales={(r['period'],r['stage']):r for r in period_diagnostics['group_rows']}
    results=[]
    for stage in ('before','after'):
        for xi in range(len(XB_BINS)):
            for ti in range(len(MINUS_TPRIME_BINS_GEV2)):
                p1=fit_one_carbon_candidate(nh3_jobs,carbon_jobs,nominal_scales,stage,xi,ti,fit_min_gev2,fit_max_gev2,1)
                p2=fit_one_carbon_candidate(nh3_jobs,carbon_jobs,nominal_scales,stage,xi,ti,fit_min_gev2,fit_max_gev2,2)
                # P2 promotion requires meaningful AICc improvement and stable peak parameters.
                p2_improvement=p1['aicc']-p2['aicc']
                mean_shift=abs(p2['shared_mean_gev2']-p1['shared_mean_gev2'])
                selected=p2 if (p2_improvement>=10.0 and mean_shift<=0.015 and p2['maximum_local_chi2_ndf']<p1['maximum_local_chi2_ndf']) else p1
                alternative=p1 if selected is p2 else p2
                selected['selected_model']=selected['residual_order']
                selected['p1_aicc']=p1['aicc']; selected['p2_aicc']=p2['aicc']; selected['p2_aicc_improvement']=p2_improvement
                selected['model_mean_difference_gev2']=mean_shift
                selected['classification']=classify_carbon_candidate(selected,alternative)
                selected['candidate_summaries']={
                    'P1':{k:p1[k] for k in ('aicc','joint_chi2_ndf','shared_mean_gev2','shared_mean_error_gev2','maximum_local_chi2_ndf')},
                    'P2':{k:p2[k] for k in ('aicc','joint_chi2_ndf','shared_mean_gev2','shared_mean_error_gev2','maximum_local_chi2_ndf')},
                }
                # Deterministic fit-range variations.
                variations=[]
                for lo,hi in ((fit_min_gev2,1.30),(fit_min_gev2,1.40),(0.45,fit_max_gev2)):
                    if lo>=hi or abs(lo-fit_min_gev2)<1e-9 and abs(hi-fit_max_gev2)<1e-9: continue
                    v=fit_one_carbon_candidate(nh3_jobs,carbon_jobs,nominal_scales,stage,xi,ti,lo,hi,selected['residual_order'])
                    variations.append({'fit_min_gev2':lo,'fit_max_gev2':hi,'shared_mean_gev2':v['shared_mean_gev2'],
                                       'delta_mean_gev2':v['shared_mean_gev2']-selected['shared_mean_gev2'],'joint_chi2_ndf':v['joint_chi2_ndf']})
                # Control endpoint variations: rebuild period-wide count scales.
                endpoint_vars=[]
                for endpoint in sorted(set(control_endpoints_gev2)):
                    diag=build_carbon_phase1_diagnostics(nh3_jobs,carbon_jobs,'period',control_min_gev2,endpoint,(endpoint,),
                                                         nh3_jobs[0].histogram_min_gev2,nh3_jobs[0].histogram_max_gev2,nh3_jobs[0].histogram_bins,
                                                         (0.40,0.65),(1.15,1.35))
                    scales={(r['period'],r['stage']):r for r in diag['group_rows']}
                    v=fit_one_carbon_candidate(nh3_jobs,carbon_jobs,scales,stage,xi,ti,fit_min_gev2,fit_max_gev2,selected['residual_order'])
                    endpoint_vars.append({'control_max_gev2':endpoint,'shared_mean_gev2':v['shared_mean_gev2'],
                                          'delta_mean_gev2':v['shared_mean_gev2']-selected['shared_mean_gev2'],'joint_chi2_ndf':v['joint_chi2_ndf']})
                selected['fit_range_variations']=variations; selected['control_endpoint_variations']=endpoint_vars
                selected['maximum_abs_fit_range_delta_mean_gev2']=max((abs(v['delta_mean_gev2']) for v in variations),default=0.0)
                selected['maximum_abs_control_delta_mean_gev2']=max((abs(v['delta_mean_gev2']) for v in endpoint_vars),default=0.0)
                for r in selected['period_records']:
                    r.update({'selected_model':selected['selected_model'],'classification':selected['classification'],
                              'p1_aicc':selected['p1_aicc'],'p2_aicc':selected['p2_aicc'],
                              'p2_aicc_improvement':selected['p2_aicc_improvement'],
                              'model_mean_difference_gev2':selected['model_mean_difference_gev2'],
                              'maximum_abs_fit_range_delta_mean_gev2':selected['maximum_abs_fit_range_delta_mean_gev2'],
                              'maximum_abs_control_delta_mean_gev2':selected['maximum_abs_control_delta_mean_gev2']})
                # endfor
                results.append(selected)
            # endfor
        # endfor
    # endfor
    return results


def plot_carbon_template_fit_canvases(
    results: list[dict[str, Any]],
    output_dir: Path,
) -> list[Path]:
    """Draw selected carbon-template fits with clean display components."""
    ensure_directory(output_dir)
    lookup = {
        (result["stage"], result["x_index"], result["t_index"]): result
        for result in results
    }
    outputs: list[Path] = []

    for period in ("su22", "fa22", "sp23"):
        for stage in ("before", "after"):
            fig, axes = plt.subplots(
                len(XB_BINS),
                len(MINUS_TPRIME_BINS_GEV2),
                figsize=(24, 15),
                sharex=True,
                sharey=False,
            )
            for x_index in range(len(XB_BINS)):
                for t_index in range(len(MINUS_TPRIME_BINS_GEV2)):
                    ax = axes[x_index, t_index]
                    joint = lookup[(stage, x_index, t_index)]
                    payload = joint["payloads"][period]
                    record = payload["record"]
                    centers = np.asarray(payload["centers"], dtype=float)
                    carbon = np.asarray(payload["carbon"], dtype=float)
                    signal = np.asarray(payload["signal"], dtype=float)
                    residual = np.asarray(payload["residual"], dtype=float)

                    dense_x = np.linspace(
                        float(np.min(centers)),
                        float(np.max(centers)),
                        800,
                    )
                    dense_residual = carbon_residual_background(
                        dense_x,
                        np.asarray(
                            record["residual_coefficients"],
                            dtype=float,
                        ),
                    )
                    display_carbon = smooth_template_for_display(carbon)
                    dense_carbon = np.interp(
                        dense_x,
                        centers,
                        display_carbon,
                    )
                    dense_signal = np.interp(
                        dense_x,
                        centers,
                        signal,
                    )
                    dense_total = (
                        dense_carbon
                        + dense_signal
                        + dense_residual
                    )

                    ax.errorbar(
                        centers,
                        payload["nh3"],
                        yerr=payload["error"],
                        fmt="o",
                        markersize=2.0,
                        linewidth=0.7,
                        label="NH$_3$",
                        zorder=5,
                    )
                    ax.step(
                        centers,
                        carbon,
                        where="mid",
                        linewidth=1.0,
                        label="Constrained carbon",
                        zorder=2,
                    )
                    ax.plot(
                        dense_x,
                        dense_signal,
                        linewidth=1.2,
                        label="Gaussian signal",
                        zorder=3,
                    )
                    ax.plot(
                        dense_x,
                        dense_residual,
                        linestyle="--",
                        linewidth=1.0,
                        label=(f"Residual $P_{joint['selected_model']}$ "       f"($M_x^2\\geq{CARBON_RESIDUAL_THRESHOLD_GEV2:.2f}$)"),
                        zorder=3,
                    )
                    ax.plot(
                        dense_x,
                        dense_total,
                        linewidth=1.5,
                        label="Total model (display-smoothed)",
                        zorder=4,
                    )
                    ax.axvline(
                        CARBON_RESIDUAL_THRESHOLD_GEV2,
                        linestyle="--",
                        linewidth=0.7,
                        alpha=0.7,
                    )
                    ax.axvline(
                        NEUTRON_MASS2_GEV2,
                        linestyle=":",
                        linewidth=0.8,
                    )
                    ax.set_title(
                        f"Bin {record['bin_number']}: "
                        f"{joint['classification']}, "
                        f"$\\mu$={record['shared_mean_gev2']:.4f}, "
                        f"$\\chi^2_\\nu$={record['local_chi2_ndf']:.2f}",
                        fontsize=8,
                    )
                    ax.grid(alpha=0.22)
                    if x_index == len(XB_BINS) - 1:
                        ax.set_xlabel("$M_x^2$ (GeV$^2$)")
                    # endif
                    if t_index == 0:
                        x_low, x_high = XB_BINS[x_index]
                        ax.set_ylabel(
                            f"{x_low:.2f} < $x_B$ < {x_high:.2f}\nCounts"
                        )
                    # endif
                # endfor
            # endfor

            handles, labels = axes[0, 0].get_legend_handles_labels()
            fig.suptitle(
                f"{PERIOD_LABELS[period]}: {STAGE_LABELS[stage]}\n"
                "Candidate Gaussian + constrained carbon + "
                "residual-polynomial fits",
                fontsize=15,
                y=0.94,
            )
            fig.legend(
                handles,
                labels,
                loc="upper center",
                bbox_to_anchor=(0.5, 0.885),
                ncol=5,
                frameon=True,
            )
            fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.82))
            path = output_dir / (
                f"carbon_template_fits_{period}_{stage}_v27.png"
            )
            fig.savefig(path, dpi=180)
            plt.close(fig)
            outputs.append(path)
        # endfor
    # endfor
    return outputs


def build_carbon_comparison_table(results: list[dict[str, Any]], production_frame: pd.DataFrame) -> pd.DataFrame:
    """Compare selected carbon fits directly with existing polynomial fits."""
    rows=[]
    for joint in results:
        stage=joint['stage']; xi=joint['x_index']; ti=joint['t_index']
        prod=production_frame[(production_frame['stage']==stage)&(production_frame['x_index']==xi)&(production_frame['t_index']==ti)&production_frame['is_recommended']]
        prod_mu=float(prod['shared_mean_gev2'].dropna().iloc[0]) if 'shared_mean_gev2' in prod and not prod['shared_mean_gev2'].dropna().empty else math.nan
        rows.append({'stage':stage,'x_index':xi,'t_index':ti,'bin_number':combined_bin_number(xi,ti),'bin_id':bin_identifier(xi,ti),
                     'carbon_classification':joint['classification'],'carbon_selected_residual_order':joint['selected_model'],
                     'carbon_shared_mean_gev2':joint['shared_mean_gev2'],'carbon_shared_mean_error_gev2':joint['shared_mean_error_gev2'],
                     'carbon_neutron_mass_deviation_gev2':joint['shared_mean_gev2']-NEUTRON_MASS2_GEV2,
                     'carbon_joint_chi2_ndf':joint['joint_chi2_ndf'],'carbon_maximum_local_chi2_ndf':joint['maximum_local_chi2_ndf'],
                     'carbon_model_mean_difference_gev2':joint['model_mean_difference_gev2'],
                     'carbon_fit_range_stability_gev2':joint['maximum_abs_fit_range_delta_mean_gev2'],
                     'carbon_control_endpoint_stability_gev2':joint['maximum_abs_control_delta_mean_gev2'],
                     'polynomial_shared_mean_gev2':prod_mu,'polynomial_neutron_mass_deviation_gev2':prod_mu-NEUTRON_MASS2_GEV2 if np.isfinite(prod_mu) else math.nan,
                     'carbon_minus_polynomial_mean_gev2':joint['shared_mean_gev2']-prod_mu if np.isfinite(prod_mu) else math.nan})
    # endfor
    return pd.DataFrame(rows)


def write_carbon_template_fit_products(results: list[dict[str, Any]], production_frame: pd.DataFrame, table_dir: Path, plot_dir: Path) -> dict[str, Any]:
    """Write candidate fit records, stability details, comparison, and plots."""
    ensure_directory(table_dir); ensure_directory(plot_dir)
    rows=[r for j in results for r in j['period_records']]
    csv_path=table_dir/'carbon_template_fit_results_v27.csv'; json_path=table_dir/'carbon_template_fit_results_v27.json'
    comparison_path=table_dir/'carbon_vs_polynomial_comparison_v27.csv'
    pd.DataFrame(rows).to_csv(csv_path,index=False)
    comparison=build_carbon_comparison_table(results,production_frame); comparison.to_csv(comparison_path,index=False)
    compact=[]
    for j in results:
        compact.append({k:v for k,v in j.items() if k!='payloads'})
    # endfor
    json_path.write_text(json.dumps(json_safe({'description':'Candidate constrained-carbon fits; original production fits remain unchanged.','results':compact}),indent=2)+'\n',encoding='utf-8')
    plots=plot_carbon_template_fit_canvases(results,plot_dir)
    return {'csv':csv_path,'json':json_path,'comparison_csv':comparison_path,'plots':plots}



# =============================================================================
# Final carbon-assisted cut products and method comparisons
# =============================================================================

def _selected_corrected_carbon_results(
    results: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Return corrected selected carbon-assisted results ordered by bin."""
    selected = [
        item
        for item in results
        if item["stage"] == "after"
    ]
    selected.sort(
        key=lambda item: combined_bin_number(
            item["x_index"],
            item["t_index"],
        )
    )
    return selected


def build_final_carbon_cut_table(
    results: list[dict[str, Any]],
) -> pd.DataFrame:
    """
    Build period-specific tight, nominal, and loose Mx2 windows.

    Each kinematic bin retains the shared fitted mean from the simultaneous
    corrected-data fit. Su22, Fa22, and Sp23 retain their independently fitted
    Gaussian widths and therefore receive independent 2-, 3-, and 4-sigma
    windows.
    """
    rows: list[dict[str, Any]] = []
    periods = ("su22", "fa22", "sp23")

    for joint in _selected_corrected_carbon_results(results):
        period_records = {
            record["period"]: record
            for record in joint["period_records"]
        }
        mean = float(joint["shared_mean_gev2"])
        mean_error = float(joint["shared_mean_error_gev2"])

        row: dict[str, Any] = {
            "x_index": int(joint["x_index"]),
            "t_index": int(joint["t_index"]),
            "bin_number": combined_bin_number(
                joint["x_index"],
                joint["t_index"],
            ),
            "bin_id": bin_identifier(
                joint["x_index"],
                joint["t_index"],
            ),
            "xB_min": XB_BINS[joint["x_index"]][0],
            "xB_max": XB_BINS[joint["x_index"]][1],
            "minus_tprime_min_gev2": (
                MINUS_TPRIME_BINS_GEV2[joint["t_index"]][0]
            ),
            "minus_tprime_max_gev2": (
                MINUS_TPRIME_BINS_GEV2[joint["t_index"]][1]
            ),
            "shared_mean_gev2": mean,
            "shared_mean_error_gev2": mean_error,
            "selected_residual_order": int(joint["selected_model"]),
            "classification": str(joint["classification"]),
            "joint_chi2_ndf": float(joint["joint_chi2_ndf"]),
            "maximum_local_chi2_ndf": float(
                joint["maximum_local_chi2_ndf"]
            ),
            "fit_range_stability_gev2": float(
                joint["maximum_abs_fit_range_delta_mean_gev2"]
            ),
            "control_endpoint_stability_gev2": float(
                joint["maximum_abs_control_delta_mean_gev2"]
            ),
            "model_mean_difference_gev2": float(
                joint["model_mean_difference_gev2"]
            ),
        }

        sigma_values: list[float] = []
        for period in periods:
            record = period_records[period]
            sigma = float(record["sigma_gev2"])
            sigma_error = float(record["sigma_error_gev2"])
            sigma_values.append(sigma)
            row[f"sigma_{period}_gev2"] = sigma
            row[f"sigma_{period}_error_gev2"] = sigma_error

            for label, multiple in (
                ("tight", 2.0),
                ("nominal", 3.0),
                ("loose", 4.0),
            ):
                row[f"{label}_{period}_sigma_multiple"] = multiple
                row[f"{label}_{period}_min_gev2"] = (
                    mean - multiple * sigma
                )
                row[f"{label}_{period}_max_gev2"] = (
                    mean + multiple * sigma
                )
            # endfor
        # endfor

        row["mean_sigma_for_summary_gev2"] = float(
            np.mean(sigma_values)
        )
        rows.append(row)
    # endfor

    frame = pd.DataFrame(rows)
    frame.sort_values("bin_number", inplace=True, ignore_index=True)
    return frame

def write_python_cut_lookup(
    cuts: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write an importable period-first Python lookup for downstream stages."""
    periods = ("su22", "fa22", "sp23")
    lines = [
        '"""Auto-generated period-specific carbon-assisted Mx2 cuts."""',
        "",
        'CUT_VERSION = "v21"',
        'CUT_UNITS = "GeV^2"',
        "SIGMA_MULTIPLES = {\"tight\": 2.0, \"nominal\": 3.0, \"loose\": 4.0}",
        "",
        "MX2_CUTS_GEV2 = {",
    ]
    for period in periods:
        lines.append(f'    "{period}": {{')
        for row in cuts.itertuples(index=False):
            lines.extend(
                [
                    f"        ({row.x_index}, {row.t_index}): {{",
                    f'            "tight": ({getattr(row, f"tight_{period}_min_gev2"):.12g}, '
                    f'{getattr(row, f"tight_{period}_max_gev2"):.12g}),',
                    f'            "nominal": ({getattr(row, f"nominal_{period}_min_gev2"):.12g}, '
                    f'{getattr(row, f"nominal_{period}_max_gev2"):.12g}),',
                    f'            "loose": ({getattr(row, f"loose_{period}_min_gev2"):.12g}, '
                    f'{getattr(row, f"loose_{period}_max_gev2"):.12g}),',
                    f'            "mu": {row.shared_mean_gev2:.12g},',
                    f'            "mu_error": {row.shared_mean_error_gev2:.12g},',
                    f'            "sigma": {getattr(row, f"sigma_{period}_gev2"):.12g},',
                    f'            "sigma_error": {getattr(row, f"sigma_{period}_error_gev2"):.12g},',
                    f'            "classification": "{row.classification}",',
                    "        },",
                ]
            )
        # endfor
        lines.append("    },")
    # endfor
    lines.extend(
        [
            "}",
            "",
            "def get_mx2_cut(period, x_index, t_index, variation=\"nominal\"):",
            "    return MX2_CUTS_GEV2[period][(x_index, t_index)][variation]",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")

def write_cpp_cut_lookup(
    cuts: pd.DataFrame,
    output_path: Path,
) -> None:
    """Write a C++ header with period-specific cuts and a lookup helper."""
    periods = ("su22", "fa22", "sp23")
    enum_names = {"su22": "Su22", "fa22": "Fa22", "sp23": "Sp23"}
    lines = [
        "#pragma once",
        "",
        "#include <array>",
        "#include <stdexcept>",
        "",
        "namespace rgc_enpiplus {",
        "",
        "enum class RunPeriod { Su22, Fa22, Sp23 };",
        "enum class CutVariation { Tight, Nominal, Loose };",
        "",
        "struct Mx2Window { double low; double high; };",
        "",
        "struct Mx2CutEntry {",
        "    RunPeriod period;",
        "    int x_index;",
        "    int t_index;",
        "    Mx2Window tight;",
        "    Mx2Window nominal;",
        "    Mx2Window loose;",
        "    double mu;",
        "    double mu_error;",
        "    double sigma;",
        "    double sigma_error;",
        "};",
        "",
        f"inline constexpr std::array<Mx2CutEntry, {len(cuts) * len(periods)}> ",
        "kMx2CutsGeV2 = {{",
    ]
    for period in periods:
        for row in cuts.itertuples(index=False):
            lines.append(
                "    Mx2CutEntry{"
                f"RunPeriod::{enum_names[period]}, {row.x_index}, {row.t_index}, "
                f"Mx2Window{{{getattr(row, f'tight_{period}_min_gev2'):.12g}, "
                f"{getattr(row, f'tight_{period}_max_gev2'):.12g}}}, "
                f"Mx2Window{{{getattr(row, f'nominal_{period}_min_gev2'):.12g}, "
                f"{getattr(row, f'nominal_{period}_max_gev2'):.12g}}}, "
                f"Mx2Window{{{getattr(row, f'loose_{period}_min_gev2'):.12g}, "
                f"{getattr(row, f'loose_{period}_max_gev2'):.12g}}}, "
                f"{row.shared_mean_gev2:.12g}, {row.shared_mean_error_gev2:.12g}, "
                f"{getattr(row, f'sigma_{period}_gev2'):.12g}, "
                f"{getattr(row, f'sigma_{period}_error_gev2'):.12g}"
                "},"
            )
        # endfor
    # endfor
    lines.extend(
        [
            "};",
            "",
            "inline Mx2Window GetMissingNeutronCut(",
            "    RunPeriod period, int x_index, int t_index,",
            "    CutVariation variation = CutVariation::Nominal) {",
            "    for (const auto& entry : kMx2CutsGeV2) {",
            "        if (entry.period == period && entry.x_index == x_index &&",
            "            entry.t_index == t_index) {",
            "            if (variation == CutVariation::Tight) return entry.tight;",
            "            if (variation == CutVariation::Loose) return entry.loose;",
            "            return entry.nominal;",
            "        }",
            "    }",
            '    throw std::out_of_range("Missing Mx2 cut entry");',
            "}",
            "",
            "}  // namespace rgc_enpiplus",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")

def plot_final_cut_canvases(
    results: list[dict[str, Any]],
    cuts: pd.DataFrame,
    output_dir: Path,
) -> list[Path]:
    """
    Plot corrected carbon-assisted fits with clearly separated components.

    The fitted carbon template remains a binned step histogram. The Gaussian
    and residual polynomial are smooth curves. A lightly smoothed copy of the
    carbon template is used only to draw the total-model curve; the statistical
    fit always uses the original unsmoothed carbon histogram bins.
    """
    ensure_directory(output_dir)
    cut_lookup = {
        (int(row.x_index), int(row.t_index)): row
        for row in cuts.itertuples(index=False)
    }
    result_lookup_ = {
        (item["x_index"], item["t_index"]): item
        for item in _selected_corrected_carbon_results(results)
    }
    outputs: list[Path] = []

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    for period in ("su22", "fa22", "sp23"):
        fig, axes = plt.subplots(
            len(XB_BINS),
            len(MINUS_TPRIME_BINS_GEV2),
            figsize=(24, 15),
            sharex=True,
            sharey=False,
        )
        for x_index in range(len(XB_BINS)):
            for t_index in range(len(MINUS_TPRIME_BINS_GEV2)):
                ax = axes[x_index, t_index]
                joint = result_lookup_[(x_index, t_index)]
                payload = joint["payloads"][period]
                cut = cut_lookup[(x_index, t_index)]
                centers = np.asarray(payload["centers"], dtype=float)

                # Finer x-grid for smooth display curves only.
                dense_x = np.linspace(
                    float(np.min(centers)),
                    float(np.max(centers)),
                    800,
                )
                gaussian_dense = np.interp(
                    dense_x,
                    centers,
                    np.asarray(payload["signal"], dtype=float),
                )
                residual_dense = carbon_residual_background(
                    dense_x,
                    np.asarray(
                        payload["record"]["residual_coefficients"],
                        dtype=float,
                    ),
                )
                raw_carbon = np.asarray(payload["carbon"], dtype=float)
                display_carbon = smooth_template_for_display(raw_carbon)
                carbon_dense = np.interp(
                    dense_x,
                    centers,
                    display_carbon,
                )
                total_dense = carbon_dense + gaussian_dense + residual_dense

                ax.errorbar(
                    centers,
                    payload["nh3"],
                    yerr=payload["error"],
                    fmt="o",
                    markersize=2.0,
                    linewidth=0.7,
                    label="NH$_3$",
                    zorder=5,
                )
                ax.step(
                    centers,
                    payload["carbon"],
                    where="mid",
                    linewidth=1.0,
                    label="Constrained carbon",
                    zorder=2,
                )
                ax.plot(
                    dense_x,
                    gaussian_dense,
                    linewidth=1.2,
                    label="Gaussian signal",
                    zorder=3,
                )
                ax.plot(
                    dense_x,
                    residual_dense,
                    linewidth=1.0,
                    linestyle="--",
                    label=(f"Residual $P_{joint['selected_model']}$ "       f"($M_x^2\\geq{CARBON_RESIDUAL_THRESHOLD_GEV2:.2f}$)"),
                    zorder=3,
                )
                ax.plot(
                    dense_x,
                    total_dense,
                    linewidth=1.5,
                    label="Total model (display-smoothed)",
                    zorder=4,
                )

                ax.axvline(
                    CARBON_RESIDUAL_THRESHOLD_GEV2,
                    linestyle="--",
                    linewidth=0.7,
                    alpha=0.7,
                    zorder=1,
                )

                # Shaded nominal selection and lighter systematic boundaries.
                ax.axvspan(
                    getattr(cut, f"nominal_{period}_min_gev2"),
                    getattr(cut, f"nominal_{period}_max_gev2"),
                    alpha=0.12,
                    zorder=0,
                )
                ax.axvline(
                    getattr(cut, f"nominal_{period}_min_gev2"),
                    linestyle="--",
                    linewidth=1.0,
                    zorder=1,
                )
                ax.axvline(
                    getattr(cut, f"nominal_{period}_max_gev2"),
                    linestyle="--",
                    linewidth=1.0,
                    zorder=1,
                )
                for boundary in (
                    getattr(cut, f"tight_{period}_min_gev2"),
                    getattr(cut, f"tight_{period}_max_gev2"),
                ):
                    ax.axvline(
                        boundary,
                        linestyle=":",
                        linewidth=0.8,
                        alpha=0.8,
                        zorder=1,
                    )
                # endfor
                for boundary in (
                    getattr(cut, f"loose_{period}_min_gev2"),
                    getattr(cut, f"loose_{period}_max_gev2"),
                ):
                    ax.axvline(
                        boundary,
                        linestyle="-.",
                        linewidth=0.8,
                        alpha=0.8,
                        zorder=1,
                    )
                # endfor

                period_sigma = getattr(
                    cut,
                    f"sigma_{period}_gev2",
                )
                period_sigma_error = getattr(
                    cut,
                    f"sigma_{period}_error_gev2",
                )
                ax.set_title(
                    f"Bin {cut.bin_number}: {cut.classification}, "
                    f"$\\mu$={cut.shared_mean_gev2:.4f}, "
                    f"$\\sigma_{{{PERIOD_LABELS[period]}}}$="
                    f"{period_sigma:.4f}$\\pm${period_sigma_error:.4f}",
                    fontsize=8,
                )
                ax.grid(alpha=0.22)
                if x_index == len(XB_BINS) - 1:
                    ax.set_xlabel("$M_x^2$ (GeV$^2$)")
                # endif
                if t_index == 0:
                    low_x, high_x = XB_BINS[x_index]
                    ax.set_ylabel(
                        f"{low_x:.2f} < $x_B$ < {high_x:.2f}\nCounts"
                    )
                # endif
            # endfor
        # endfor

        handles = [
            Line2D([0], [0], marker="o", linestyle="none", label="NH$_3$"),
            Line2D([0], [0], linestyle="-", label="Constrained carbon"),
            Line2D([0], [0], linestyle="-", label="Gaussian signal"),
            Line2D([0], [0], linestyle="--", label="Residual polynomial"),
            Line2D([0], [0], linestyle="-", linewidth=1.5, label="Total model (display-smoothed)"),
            Patch(alpha=0.12, label="Nominal (3$\\sigma$)"),
            Line2D([0], [0], linestyle=":", label="Tight (2$\\sigma$)"),
            Line2D([0], [0], linestyle="-.", label="Loose (4$\\sigma$)"),
        ]
        labels = [handle.get_label() for handle in handles]

        fig.suptitle(
            f"{PERIOD_LABELS[period]}: After momentum corrections\n"
            "Final carbon-assisted missing-neutron cut definitions",
            fontsize=15,
            y=0.955,
        )
        fig.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.895),
            ncol=8,
            frameon=True,
        )
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.84))
        output_path = output_dir / (
            f"final_carbon_assisted_cuts_{period}_v27.png"
        )
        fig.savefig(output_path, dpi=180)
        plt.close(fig)
        outputs.append(output_path)
    # endfor

    return outputs


def plot_carbon_vs_polynomial_mu_sigma(
    carbon_results: list[dict[str, Any]],
    production_frame: pd.DataFrame,
    output_path: Path,
) -> Path:
    """
    Compare fitted mu and sigma values from carbon and polynomial methods.

    Top: shared fitted mu values.
    Middle: delta mu = mu_carbon - mu_polynomial.
    Bottom: period-specific fitted sigma values.
    """
    corrected_carbon = _selected_corrected_carbon_results(carbon_results)
    carbon_by_bin = {
        combined_bin_number(item["x_index"], item["t_index"]): item
        for item in corrected_carbon
    }
    polynomial = production_frame[
        (production_frame["stage"] == "after")
        & production_frame["is_recommended"]
    ].copy()

    bins = np.arange(
        1,
        len(XB_BINS) * len(MINUS_TPRIME_BINS_GEV2) + 1,
    )
    carbon_mu: list[float] = []
    carbon_mu_error: list[float] = []
    polynomial_mu: list[float] = []
    polynomial_mu_error: list[float] = []
    classifications: list[str] = []

    carbon_sigma: dict[str, list[float]] = {
        period: [] for period in ("su22", "fa22", "sp23")
    }
    polynomial_sigma: dict[str, list[float]] = {
        period: [] for period in ("su22", "fa22", "sp23")
    }

    for bin_number in bins:
        joint = carbon_by_bin[int(bin_number)]
        carbon_mu.append(float(joint["shared_mean_gev2"]))
        carbon_mu_error.append(float(joint["shared_mean_error_gev2"]))
        classifications.append(str(joint["classification"]))

        polynomial_bin = polynomial[
            polynomial["bin_number"] == bin_number
        ]
        shared_values = polynomial_bin["shared_mean_gev2"].dropna()
        shared_errors = polynomial_bin[
            "shared_mean_error_gev2"
        ].dropna()
        polynomial_mu.append(
            float(shared_values.iloc[0])
            if not shared_values.empty
            else math.nan
        )
        polynomial_mu_error.append(
            float(shared_errors.iloc[0])
            if not shared_errors.empty
            else math.nan
        )

        carbon_period_records = {
            record["period"]: record
            for record in joint["period_records"]
        }
        for period in ("su22", "fa22", "sp23"):
            carbon_sigma[period].append(
                float(carbon_period_records[period]["sigma_gev2"])
            )
            polynomial_row = polynomial_bin[
                polynomial_bin["period"] == period
            ]
            polynomial_sigma[period].append(
                float(polynomial_row["sigma_gev2"].iloc[0])
                if not polynomial_row.empty
                else math.nan
            )
        # endfor
    # endfor

    carbon_mu_array = np.asarray(carbon_mu, dtype=float)
    polynomial_mu_array = np.asarray(polynomial_mu, dtype=float)
    delta_mu = carbon_mu_array - polynomial_mu_array
    resolved_mask = np.asarray(
        [value == "resolved" for value in classifications],
        dtype=bool,
    )
    marginal_mask = ~resolved_mask

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(16, 14),
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 0.8, 1.25]},
    )

    ax_mu = axes[0]
    ax_mu.errorbar(
        bins,
        carbon_mu_array,
        yerr=np.asarray(carbon_mu_error, dtype=float),
        marker="o",
        linewidth=1.2,
        markersize=4,
        label="Carbon-assisted",
    )
    ax_mu.errorbar(
        bins,
        polynomial_mu_array,
        yerr=np.asarray(polynomial_mu_error, dtype=float),
        marker="s",
        linewidth=1.0,
        markersize=3.5,
        linestyle="--",
        label="Polynomial-only",
    )
    ax_mu.scatter(
        bins[resolved_mask],
        carbon_mu_array[resolved_mask],
        marker="o",
        s=42,
        facecolors="none",
        label="Resolved",
        zorder=6,
    )
    ax_mu.scatter(
        bins[marginal_mask],
        carbon_mu_array[marginal_mask],
        marker="D",
        s=34,
        facecolors="none",
        label="Marginal",
        zorder=6,
    )
    ax_mu.axhline(
        NEUTRON_MASS2_GEV2,
        linestyle=":",
        linewidth=1.0,
        label="$m_n^2$",
    )
    ax_mu.set_ylabel("Shared fitted $\\mu$ (GeV$^2$)")
    ax_mu.set_title(
        "Carbon-assisted versus polynomial-only corrected fits"
    )
    ax_mu.grid(alpha=0.25)
    ax_mu.legend(ncol=5, fontsize=9)

    ax_delta = axes[1]
    ax_delta.axhline(0.0, linestyle=":", linewidth=1.0)
    ax_delta.scatter(
        bins[resolved_mask],
        delta_mu[resolved_mask],
        marker="o",
        s=34,
        label="Resolved",
    )
    ax_delta.scatter(
        bins[marginal_mask],
        delta_mu[marginal_mask],
        marker="D",
        s=34,
        label="Marginal",
    )
    ax_delta.plot(
        bins,
        delta_mu,
        linewidth=0.8,
        alpha=0.7,
    )
    ax_delta.set_ylabel(
        "$\\Delta\\mu=\\mu_{\\mathrm{C}}-"
        "\\mu_{\\mathrm{poly}}$ (GeV$^2$)"
    )
    ax_delta.grid(alpha=0.25)
    ax_delta.legend(ncol=2, fontsize=9)

    ax_sigma = axes[2]
    for period in ("su22", "fa22", "sp23"):
        ax_sigma.plot(
            bins,
            carbon_sigma[period],
            marker="o",
            linewidth=1.1,
            markersize=3.5,
            label=f"{PERIOD_LABELS[period]}, carbon",
        )
        ax_sigma.plot(
            bins,
            polynomial_sigma[period],
            marker="s",
            linestyle="--",
            linewidth=1.0,
            markersize=3.0,
            label=f"{PERIOD_LABELS[period]}, polynomial",
        )
    # endfor
    ax_sigma.set_xlabel("Combined kinematic-bin number")
    ax_sigma.set_ylabel("Fitted $\\sigma$ (GeV$^2$)")
    ax_sigma.set_xticks(bins)
    ax_sigma.grid(alpha=0.25)
    ax_sigma.legend(ncol=3, fontsize=9)

    fig.tight_layout()
    ensure_directory(output_path.parent)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_master_cut_summary(
    cuts: pd.DataFrame,
    output_path: Path,
) -> Path:
    """
    Plot a compact two-panel summary of the final period-specific cuts.

    The upper bands use the mean of the three fitted period widths only as a
    compact visual summary. The exported and applied cuts remain strictly
    period specific. The lower panel shows each fitted width with its fit
    uncertainty.
    """
    frame = cuts.sort_values("bin_number").copy()
    bins = frame["bin_number"].to_numpy(dtype=int)
    mean = frame["shared_mean_gev2"].to_numpy(dtype=float)
    mean_error = frame["shared_mean_error_gev2"].to_numpy(dtype=float)
    sigma_summary = frame["mean_sigma_for_summary_gev2"].to_numpy(dtype=float)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 10),
        sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0]},
    )

    ax_bounds = axes[0]
    # Draw widest first so the nested bands remain visible.
    ax_bounds.fill_between(
        bins,
        mean - 4.0 * sigma_summary,
        mean + 4.0 * sigma_summary,
        alpha=0.10,
        label="Loose (4$\\sigma$), mean period width",
    )
    ax_bounds.fill_between(
        bins,
        mean - 3.0 * sigma_summary,
        mean + 3.0 * sigma_summary,
        alpha=0.13,
        label="Nominal (3$\\sigma$), mean period width",
    )
    ax_bounds.fill_between(
        bins,
        mean - 2.0 * sigma_summary,
        mean + 2.0 * sigma_summary,
        alpha=0.17,
        label="Tight (2$\\sigma$), mean period width",
    )
    ax_bounds.errorbar(
        bins,
        mean,
        yerr=mean_error,
        marker="o",
        linewidth=1.2,
        markersize=4,
        capsize=2,
        label="Shared fitted $\\mu$",
        zorder=5,
    )
    ax_bounds.axhline(
        NEUTRON_MASS2_GEV2,
        linestyle=":",
        linewidth=1.0,
        label="$m_n^2$",
    )
    ax_bounds.set_ylabel("$M_x^2$ position (GeV$^2$)")
    ax_bounds.set_title(
        "Final carbon-assisted missing-neutron cut summary"
    )
    ax_bounds.grid(alpha=0.25)
    ax_bounds.legend(ncol=3, fontsize=9)

    ax_width = axes[1]
    period_offsets = {
        "su22": -0.12,
        "fa22": 0.00,
        "sp23": 0.12,
    }
    for period in ("su22", "fa22", "sp23"):
        x_values = bins.astype(float) + period_offsets[period]
        ax_width.errorbar(
            x_values,
            frame[f"sigma_{period}_gev2"],
            yerr=frame[f"sigma_{period}_error_gev2"],
            marker="o",
            linestyle="none",
            linewidth=1.0,
            elinewidth=1.0,
            markersize=4,
            capsize=3,
            capthick=1.0,
            label=f"{PERIOD_LABELS[period]} $\\sigma$",
            zorder=3,
        )
    # endfor
    ax_width.set_xlabel("Combined kinematic-bin number")
    ax_width.set_ylabel("Fitted $\\sigma$ (GeV$^2$)")
    ax_width.set_xticks(bins)
    ax_width.grid(alpha=0.25)
    ax_width.legend(ncol=3, fontsize=9)

    fig.tight_layout()
    ensure_directory(output_path.parent)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path

def write_final_carbon_cut_products(
    results: list[dict[str, Any]],
    production_frame: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    """Write final cut tables, lookup files, plots, and method comparison."""
    tables_dir = output_dir / "tables"
    plots_dir = output_dir / "plots"
    exports_dir = output_dir / "exports"
    for directory in (
        output_dir,
        tables_dir,
        plots_dir,
        exports_dir,
    ):
        ensure_directory(directory)
    # endfor

    cuts = build_final_carbon_cut_table(results)
    csv_path = tables_dir / "final_carbon_assisted_mx2_cuts.csv"
    json_path = tables_dir / "final_carbon_assisted_mx2_cuts.json"
    python_path = exports_dir / "final_carbon_assisted_mx2_cuts.py"
    cpp_path = exports_dir / "final_carbon_assisted_mx2_cuts.h"
    metadata_path = exports_dir / "final_carbon_assisted_mx2_cuts_metadata.json"
    comparison_plot = (
        plots_dir / "carbon_vs_polynomial_mu_sigma_v27.png"
    )
    master_summary_plot = (
        plots_dir / "final_carbon_cut_master_summary_v27.png"
    )

    cuts.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(
            json_safe(
                {
                    "method": (
                        "Shared carbon-assisted mean with independent period-specific "
                        "Gaussian widths and cut windows. The residual "
                        "P1/P2 background is exactly zero below Mx2=0.40 "
                        "GeV^2."
                    ),
                    "tight_sigma_multiple": 2.0,
                    "nominal_sigma_multiple": 3.0,
                    "loose_sigma_multiple": 4.0,
                    "periods": {
                        period: [
                            {
                                "x_index": int(row.x_index),
                                "t_index": int(row.t_index),
                                "bin_number": int(row.bin_number),
                                "bin_id": row.bin_id,
                                "xB_min": float(row.xB_min),
                                "xB_max": float(row.xB_max),
                                "minus_tprime_min_gev2": float(
                                    row.minus_tprime_min_gev2
                                ),
                                "minus_tprime_max_gev2": float(
                                    row.minus_tprime_max_gev2
                                ),
                                "mu_gev2": float(row.shared_mean_gev2),
                                "mu_error_gev2": float(
                                    row.shared_mean_error_gev2
                                ),
                                "sigma_gev2": float(
                                    getattr(row, f"sigma_{period}_gev2")
                                ),
                                "sigma_error_gev2": float(
                                    getattr(
                                        row,
                                        f"sigma_{period}_error_gev2",
                                    )
                                ),
                                "tight": [
                                    float(getattr(row, f"tight_{period}_min_gev2")),
                                    float(getattr(row, f"tight_{period}_max_gev2")),
                                ],
                                "nominal": [
                                    float(getattr(row, f"nominal_{period}_min_gev2")),
                                    float(getattr(row, f"nominal_{period}_max_gev2")),
                                ],
                                "loose": [
                                    float(getattr(row, f"loose_{period}_min_gev2")),
                                    float(getattr(row, f"loose_{period}_max_gev2")),
                                ],
                                "classification": row.classification,
                            }
                            for row in cuts.itertuples(index=False)
                        ]
                        for period in ("su22", "fa22", "sp23")
                    },
                    "flat_rows": cuts.to_dict(orient="records"),
                }
            ),
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    write_python_cut_lookup(cuts, python_path)
    write_cpp_cut_lookup(cuts, cpp_path)
    metadata_path.write_text(
        json.dumps(
            {
                "version": "v29",
                "output_filename_policy": "stable_version_independent",
                "units": "GeV^2",
                "fit_method": (
                    "Shared corrected-data mean; constrained period-wide "
                    "carbon template; adaptive thresholded P1/P2 residual that is "
                    "zero below Mx2=0.40 GeV^2; independent "
                    "period-specific Gaussian widths."
                ),
                "shared_mean_across_periods": True,
                "period_specific_widths": True,
                "sigma_multiples": {
                    "tight": 2.0,
                    "nominal": 3.0,
                    "loose": 4.0,
                },
                "carbon_control_region_gev2": [0.00, 0.40],
                "carbon_residual_background": {
                    "threshold_gev2": 0.40,
                    "below_threshold": "identically zero",
                    "at_threshold": "continuous and exactly zero",
                    "above_threshold": (
                        "adaptive P1 or P2 in powers of "
                        "(Mx2-threshold), with no constant term"
                    ),
                },
                "available_periods": ["su22", "fa22", "sp23"],
                "lookup_files": {
                    "python": python_path.name,
                    "cpp": cpp_path.name,
                    "csv": csv_path.name,
                    "json": json_path.name,
                },
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    cut_plots = plot_final_cut_canvases(
        results,
        cuts,
        plots_dir,
    )
    plot_carbon_vs_polynomial_mu_sigma(
        carbon_results=results,
        production_frame=production_frame,
        output_path=comparison_plot,
    )
    plot_master_cut_summary(
        cuts=cuts,
        output_path=master_summary_plot,
    )

    return {
        "csv": csv_path,
        "json": json_path,
        "python_lookup": python_path,
        "cpp_lookup": cpp_path,
        "metadata": metadata_path,
        "comparison_plot": comparison_plot,
        "master_summary_plot": master_summary_plot,
        "cut_plots": cut_plots,
    }




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
                            f"{item.get('fit_quality_class', 'accepted').upper()}"
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
            fit_strategy_label = (
                "Joint shared-mean Gaussian signal plus adaptively selected "
                "polynomial background"
                if stage == "after"
                else (
                    "Independent Gaussian signal plus adaptively selected "
                    "polynomial background"
                )
            )
            fig.suptitle(
                f"{PERIOD_LABELS[period]}: {STAGE_LABELS[stage]}\n"
                f"{fit_strategy_label}",
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
                f"Bin {bin_number}: {model_name}, {status}\n"
                f"$\\chi^2$/ndf$_{{\\pm 2\\sigma,\\mathrm{{approx}}}}$="
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
            "corrected_fit_strategy": {
                "method": (
                    "simultaneous Su22/Fa22/Sp23 fit in each kinematic bin"
                ),
                "shared_parameter": "Gaussian mean",
                "period_specific_parameters": (
                    "signal yield, Gaussian sigma, and background coefficients"
                ),
                "initial_mean_gev2": NEUTRON_MASS2_GEV2,
                "corrected_mean_upper_guardrail_gev2": (
                    args.corrected_mean_max
                ),
                "neighbor_mean_fallback_enabled": (
                    not args.disable_neighbor_mean_fallback
                ),
            },
            "background_initialization": {
                "method": "weighted polynomial sideband fit",
                "excluded_region_gev2": [
                    args.sideband_exclusion_min,
                    args.sideband_exclusion_max,
                ],
                "coefficients_remain_free_in_full_fit": True,
            },
            "adaptive_polynomial_selection": {
                "candidate_orders": [2, 3, 4],
                "minimum_aicc_improvement": (
                    args.joint_aicc_improvement
                ),
                "minimum_fractional_local_improvement": (
                    args.joint_local_improvement_fraction
                ),
                "accepted_worst_period_local_score_max": (
                    args.accepted_local_score_max
                ),
                "marginal_worst_period_local_score_max": (
                    args.marginal_local_score_max
                ),
                "background_positivity_region": (
                    "shared mean +/- 3 period-specific sigma"
                ),
                "background_negativity_tolerance": (
                    "max(1 count, 1 percent of the characteristic count scale)"
                ),
                "policy": (
                    "Start with the lowest viable order and increase complexity "
                    "only when the AICc and aggregate local-score improvements "
                    "meet the configured thresholds. Stop at quartic."
                ),
            },
            "signal_region_fit_quality": {
                "definition": (
                    "shared mean +/- 2 period-specific fitted sigma"
                ),
                "primary_local_score": (
                    "chi2 divided by the number of signal-region histogram bins"
                ),
                "global_joint_chi2_ndf_also_stored": True,
            },
            "fit_range_study_maxima_gev2": [1.30, 1.35, 1.40],
            "replicas_performed": args.replicas > 0,
            "replicas_per_corrected_bin": args.replicas,
            "replica_model_selection_repeated": False,
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
# Nominal-versus-ISR comparison products
# =============================================================================

def load_final_cut_frame(path: Path) -> pd.DataFrame:
    """Load and validate one final carbon-assisted cut table."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing final cut table: {path}")
    # endif
    frame = pd.read_csv(path)
    required = {
        "bin_number", "x_index", "t_index",
        "shared_mean_gev2", "shared_mean_error_gev2",
        "sigma_su22_gev2", "sigma_su22_error_gev2",
        "sigma_fa22_gev2", "sigma_fa22_error_gev2",
        "sigma_sp23_gev2", "sigma_sp23_error_gev2",
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise RuntimeError(f"Final cut table {path} is missing columns: {missing}")
    # endif
    return frame.sort_values("bin_number").reset_index(drop=True)


def plot_nominal_vs_isr_mu_sigma(
    nominal: pd.DataFrame,
    isr: pd.DataFrame,
    output_path: Path,
) -> Path:
    """Compare corrected shared means and widths for nominal and ISR."""
    ensure_directory(output_path.parent)
    bins = nominal["bin_number"].to_numpy(dtype=int)
    fig, axes = plt.subplots(2, 1, figsize=(16, 11), sharex=True)

    axes[0].errorbar(
        bins - 0.08,
        nominal["shared_mean_gev2"],
        yerr=nominal["shared_mean_error_gev2"],
        marker="o", linestyle="none", capsize=2, label="Nominal",
    )
    axes[0].errorbar(
        bins + 0.08,
        isr["shared_mean_gev2"],
        yerr=isr["shared_mean_error_gev2"],
        marker="s", linestyle="none", capsize=2, label="ISR",
    )
    axes[0].axhline(NEUTRON_MASS2_GEV2, linestyle=":", linewidth=1.0, label="$m_n^2$")
    axes[0].set_ylabel("Shared fitted $\\mu$ (GeV$^2$)")
    axes[0].set_title("Nominal versus ISR carbon-assisted corrected fits")
    axes[0].grid(alpha=0.25)
    axes[0].legend(ncol=3)

    offsets = {"su22": -0.18, "fa22": 0.0, "sp23": 0.18}
    for period in ("su22", "fa22", "sp23"):
        axes[1].errorbar(
            bins + offsets[period] - 0.045,
            nominal[f"sigma_{period}_gev2"],
            yerr=nominal[f"sigma_{period}_error_gev2"],
            marker="o", linestyle="none", capsize=2,
            label=f"{PERIOD_LABELS[period]} nominal",
        )
        axes[1].errorbar(
            bins + offsets[period] + 0.045,
            isr[f"sigma_{period}_gev2"],
            yerr=isr[f"sigma_{period}_error_gev2"],
            marker="s", linestyle="none", capsize=2,
            label=f"{PERIOD_LABELS[period]} ISR",
        )
    # endfor
    axes[1].set_xlabel("Combined kinematic-bin number")
    axes[1].set_ylabel("Fitted $\\sigma$ (GeV$^2$)")
    axes[1].set_xticks(bins)
    axes[1].grid(alpha=0.25)
    axes[1].legend(ncol=3, fontsize=9)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_nominal_vs_isr_parameter_differences(
    nominal: pd.DataFrame,
    isr: pd.DataFrame,
    output_path: Path,
) -> Path:
    """Plot ISR-minus-nominal mean and width differences."""
    ensure_directory(output_path.parent)
    bins = nominal["bin_number"].to_numpy(dtype=int)
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    delta_mu = isr["shared_mean_gev2"] - nominal["shared_mean_gev2"]
    delta_mu_error = np.sqrt(
        isr["shared_mean_error_gev2"] ** 2
        + nominal["shared_mean_error_gev2"] ** 2
    )
    axes[0].errorbar(
        bins, delta_mu, yerr=delta_mu_error,
        marker="o", linestyle="none", capsize=2,
    )
    axes[0].axhline(0.0, linestyle=":", linewidth=1.0)
    axes[0].set_ylabel("$\\mu_{ISR}-\\mu_{nominal}$ (GeV$^2$)")
    axes[0].set_title("ISR-minus-nominal fitted-parameter differences")
    axes[0].grid(alpha=0.25)

    offsets = {"su22": -0.12, "fa22": 0.0, "sp23": 0.12}
    for period in ("su22", "fa22", "sp23"):
        delta_sigma = (
            isr[f"sigma_{period}_gev2"]
            - nominal[f"sigma_{period}_gev2"]
        )
        delta_sigma_error = np.sqrt(
            isr[f"sigma_{period}_error_gev2"] ** 2
            + nominal[f"sigma_{period}_error_gev2"] ** 2
        )
        axes[1].errorbar(
            bins + offsets[period], delta_sigma,
            yerr=delta_sigma_error, marker="o",
            linestyle="none", capsize=2,
            label=PERIOD_LABELS[period],
        )
    # endfor
    axes[1].axhline(0.0, linestyle=":", linewidth=1.0)
    axes[1].set_xlabel("Combined kinematic-bin number")
    axes[1].set_ylabel("$\\sigma_{ISR}-\\sigma_{nominal}$ (GeV$^2$)")
    axes[1].set_xticks(bins)
    axes[1].grid(alpha=0.25)
    axes[1].legend(ncol=3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def write_nominal_isr_comparison_products(
    nominal_result: dict[str, Any],
    isr_result: dict[str, Any],
    output_dir: Path,
) -> dict[str, Any]:
    """Write third-location nominal-versus-ISR summaries and tables."""
    plots_dir = output_dir / "plots"
    tables_dir = output_dir / "tables"
    ensure_directory(plots_dir)
    ensure_directory(tables_dir)

    nominal = load_final_cut_frame(Path(nominal_result["final_cut_csv"]))
    isr = load_final_cut_frame(Path(isr_result["final_cut_csv"]))
    merged = nominal.merge(
        isr,
        on=["bin_number", "x_index", "t_index"],
        suffixes=("_nominal", "_isr"),
        validate="one_to_one",
    )
    merged["delta_mu_isr_minus_nominal_gev2"] = (
        merged["shared_mean_gev2_isr"]
        - merged["shared_mean_gev2_nominal"]
    )
    for period in ("su22", "fa22", "sp23"):
        merged[f"delta_sigma_{period}_isr_minus_nominal_gev2"] = (
            merged[f"sigma_{period}_gev2_isr"]
            - merged[f"sigma_{period}_gev2_nominal"]
        )
    # endfor

    csv_path = tables_dir / "nominal_vs_isr_fit_comparison_v27.csv"
    json_path = tables_dir / "nominal_vs_isr_fit_comparison_v27.json"
    mu_sigma_plot = plots_dir / "nominal_vs_isr_mu_sigma_v27.png"
    difference_plot = plots_dir / "nominal_vs_isr_parameter_differences_v27.png"

    merged.to_csv(csv_path, index=False)
    json_path.write_text(
        json.dumps(
            json_safe({
                "schema_version": 1,
                "nominal_final_cut_json": nominal_result["final_cut_json"],
                "isr_final_cut_json": isr_result["final_cut_json"],
                "rows": merged.to_dict(orient="records"),
            }),
            indent=2,
        ) + "\n",
        encoding="utf-8",
    )
    plot_nominal_vs_isr_mu_sigma(nominal, isr, mu_sigma_plot)
    plot_nominal_vs_isr_parameter_differences(nominal, isr, difference_plot)
    return {
        "csv": csv_path,
        "json": json_path,
        "mu_sigma_plot": mu_sigma_plot,
        "difference_plot": difference_plot,
    }




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
        "--carbon-input",
        action="append",
        type=parse_dataset_override,
        default=[],
        metavar="PERIOD:STAGE=FILE",
        help=(
            "Override one default carbon input for the additive Phase-1 "
            "diagnostic. May be repeated."
        ),
    )
    parser.add_argument(
        "--isr-input",
        action="append",
        type=parse_dataset_override,
        default=[],
        metavar="PERIOD:STAGE=FILE",
        help="Override one default ISR NH3 input. May be repeated.",
    )
    parser.add_argument(
        "--isr-carbon-input",
        action="append",
        type=parse_dataset_override,
        default=[],
        metavar="PERIOD:STAGE=FILE",
        help=(
            "Override one default ISR carbon input. Uncorrected defaults use "
            "*_ISR_2.root and corrected defaults use "
            "*_ISR_externalISR_mom_corrections.root."
        ),
    )
    parser.add_argument(
        "--disable-isr",
        action="store_true",
        help="Disable the parallel ISR analysis.",
    )
    parser.add_argument(
        "--disable-carbon-phase1",
        action="store_true",
        help=(
            "Disable the additive Phase-1 carbon diagnostic. The original "
            "fit workflow remains unchanged."
        ),
    )
    parser.add_argument(
        "--carbon-control-min",
        type=float,
        default=0.0,
        help="Low-Mx2 control-region minimum in GeV^2 (default: 0.0).",
    )
    parser.add_argument(
        "--carbon-control-max",
        type=float,
        default=0.40,
        help="Low-Mx2 control-region maximum in GeV^2 (default: 0.40).",
    )
    parser.add_argument(
        "--carbon-control-boundaries",
        type=parse_float_list,
        default=(0.30, 0.35, 0.40, 0.45, 0.50),
        metavar="V1,V2,...",
        help=(
            "Comma-separated control-region upper edges for the scale "
            "stability scan (default: 0.30,0.35,0.40,0.45,0.50)."
        ),
    )
    parser.add_argument(
        "--carbon-normalization-grouping",
        choices=("tprime-column", "xB-row", "bin", "period"),
        default="period",
        help=(
            "Nominal carbon normalization grouping. The default pools all "
            "24 kinematic bins within each period and correction stage. "
            "The -tprime-column result is always retained as a stability "
            "cross-check."
        ),
    )
    parser.add_argument(
        "--carbon-validation-low-min",
        type=float,
        default=0.40,
        help="Independent low-side validation minimum (default: 0.40).",
    )
    parser.add_argument(
        "--carbon-validation-low-max",
        type=float,
        default=0.65,
        help="Independent low-side validation maximum (default: 0.65).",
    )
    parser.add_argument(
        "--carbon-validation-high-min",
        type=float,
        default=1.15,
        help="Independent high-side validation minimum (default: 1.15).",
    )
    parser.add_argument(
        "--carbon-validation-high-max",
        type=float,
        default=1.35,
        help="Independent high-side validation maximum (default: 1.35).",
    )
    parser.add_argument(
        "--disable-carbon-template-fits",
        action="store_true",
        help=(
            "Disable the additive exploratory Gaussian + fixed scaled-carbon "
            "+ linear-residual fits. These fits never replace the original "
            "production fits."
        ),
    )
    parser.add_argument(
        "--carbon-template-fit-min",
        type=float,
        default=0.40,
        help=(
            "Minimum Mx2 used by the exploratory carbon-template fits in "
            "GeV^2 (default: 0.40)."
        ),
    )
    parser.add_argument(
        "--carbon-template-fit-max",
        type=float,
        default=1.35,
        help=(
            "Maximum Mx2 used by the exploratory carbon-template fits in "
            "GeV^2 (default: 1.35)."
        ),
    )
    parser.add_argument(
        "--carbon-template-control-endpoints",
        type=parse_float_list,
        default=(0.35, 0.40, 0.45),
        metavar="V1,V2,...",
        help=(
            "Control-region upper endpoints used for candidate-fit stability "
            "checks (default: 0.35,0.40,0.45)."
        ),
    )
    parser.add_argument(
        "--regular-fit-bin-factor",
        type=int,
        default=1,
        help=(
            "For the non-carbon polynomial fits only, reduce the requested "
            "histogram-bin count by this integer factor. This test revision "
            "defaults to 1, i.e. the standard binning. Carbon diagnostics and "
            "carbon-assisted fits retain the full --hist-bins value."
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
        default=0.0,
        help="Mx2 histogram minimum in GeV^2 (default: 0.0).",
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
        default=0.930,
        help=(
            "Upper bound on the corrected-data Gaussian mean in GeV^2 "
            "(default: 0.930). This is a diagnostic guardrail rather than a "
            "physics constraint; before-correction fits retain a wider bound."
        ),
    )
    parser.add_argument(
        "--minimum-events",
        type=int,
        default=150,
        help="Minimum events required in a kinematic bin (default: 150).",
    )
    parser.add_argument(
        "--joint-aicc-improvement",
        type=float,
        default=10.0,
        help=(
            "Minimum AICc decrease required before increasing the corrected "
            "background order (default: 10.0)."
        ),
    )
    parser.add_argument(
        "--joint-local-improvement-fraction",
        type=float,
        default=0.20,
        help=(
            "Minimum fractional improvement in the aggregate corrected local "
            "score before increasing background order (default: 0.20)."
        ),
    )
    parser.add_argument(
        "--accepted-local-score-max",
        type=float,
        default=2.0,
        help=(
            "Worst-period local chi2-per-bin threshold for accepted corrected "
            "fits (default: 2.0)."
        ),
    )
    parser.add_argument(
        "--marginal-local-score-max",
        type=float,
        default=3.0,
        help=(
            "Worst-period local chi2-per-bin threshold for marginal corrected "
            "fits (default: 3.0)."
        ),
    )
    parser.add_argument(
        "--sideband-exclusion-min",
        type=float,
        default=0.74,
        help=(
            "Lower edge of the provisional neutron region excluded during "
            "background initialization in GeV^2 (default: 0.74)."
        ),
    )
    parser.add_argument(
        "--sideband-exclusion-max",
        type=float,
        default=1.04,
        help=(
            "Upper edge of the provisional neutron region excluded during "
            "background initialization in GeV^2 (default: 1.04)."
        ),
    )
    parser.add_argument(
        "--disable-neighbor-mean-fallback",
        action="store_true",
        help=(
            "Disable the deterministic neighboring-bin fixed-mean fallback "
            "for unresolved corrected bins."
        ),
    )
    parser.add_argument(
        "--replicas",
        type=int,
        default=0,
        help=(
            "Optional Poisson replicas per corrected bin after model freezing "
            "(default: 0). Use 100 for development and 1000 or more for final."
        ),
    )
    parser.add_argument(
        "--replica-seed",
        type=int,
        default=20260724,
        help="Random seed for optional Poisson replicas.",
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
    if args.joint_aicc_improvement < 0.0:
        raise ValueError("--joint-aicc-improvement must be nonnegative.")
    # endif
    if not 0.0 <= args.joint_local_improvement_fraction <= 1.0:
        raise ValueError(
            "--joint-local-improvement-fraction must lie between 0 and 1."
        )
    # endif
    if args.accepted_local_score_max <= 0.0:
        raise ValueError("--accepted-local-score-max must be positive.")
    # endif
    if args.marginal_local_score_max <= args.accepted_local_score_max:
        raise ValueError(
            "--marginal-local-score-max must exceed "
            "--accepted-local-score-max."
        )
    # endif
    if args.sideband_exclusion_min >= args.sideband_exclusion_max:
        raise ValueError(
            "--sideband-exclusion-min must be below "
            "--sideband-exclusion-max."
        )
    # endif
    if args.replicas < 0:
        raise ValueError("--replicas must be nonnegative.")
    # endif
    if args.regular_fit_bin_factor < 1:
        raise ValueError("--regular-fit-bin-factor must be at least 1.")
    # endif
    if args.hist_bins // args.regular_fit_bin_factor < 10:
        raise ValueError(
            "The reduced regular-fit histogram must contain at least 10 bins."
        )
    # endif
    if args.carbon_control_min >= args.carbon_control_max:
        raise ValueError(
            "--carbon-control-min must be below --carbon-control-max."
        )
    # endif
    if (
        args.carbon_control_min < args.hist_min
        or args.carbon_control_max > args.hist_max
    ):
        raise ValueError(
            "The carbon control region must lie inside the histogram range."
        )
    # endif
    if any(
        value <= args.carbon_control_min or value > args.hist_max
        for value in args.carbon_control_boundaries
    ):
        raise ValueError(
            "Carbon scan boundaries must exceed the control minimum and "
            "remain inside the histogram range."
        )
    # endif
    if (
        args.carbon_validation_low_min
        >= args.carbon_validation_low_max
        or args.carbon_validation_high_min
        >= args.carbon_validation_high_max
    ):
        raise ValueError("Carbon validation-region bounds are invalid.")
    # endif
    if (
        args.carbon_template_fit_min
        >= args.carbon_template_fit_max
    ):
        raise ValueError(
            "--carbon-template-fit-min must be below "
            "--carbon-template-fit-max."
        )
    # endif
    if (
        args.carbon_template_fit_min < args.hist_min
        or args.carbon_template_fit_max > args.hist_max
    ):
        raise ValueError(
            "The carbon-template fit range must lie inside the histogram range."
        )
    # endif


def run_analysis_variant(
    args: argparse.Namespace,
    sample_variant: str,
    input_paths: dict[str, dict[str, Path]],
    carbon_input_paths: dict[str, dict[str, Path]],
    output_dir: Path,
) -> dict[str, Any]:
    """Run the complete existing workflow for one sample variant."""
    branch_overrides = dict(args.branch)
    datasets = preflight_inputs(
        input_paths=input_paths,
        requested_tree=args.tree,
        branch_overrides=branch_overrides,
    )

    carbon_datasets: list[InputDataset] = []
    if not args.disable_carbon_phase1:
        carbon_datasets = preflight_inputs(
            input_paths=carbon_input_paths,
            requested_tree=args.tree,
            branch_overrides=branch_overrides,
        )
    # endif

    output_dir = output_dir.expanduser().resolve()
    table_dir = output_dir / "tables"
    plot_dir = output_dir / "plots"
    spectrum_dir = plot_dir / "spectra"
    model_canvas_dir = plot_dir / "background_models"
    summary_dir = plot_dir / "summaries"
    pull_diagnostic_dir = plot_dir / "problem_bin_diagnostics"
    carbon_phase1_dir = output_dir / "carbon_phase1"
    carbon_template_fit_dir = output_dir / "carbon_template_fits"
    carbon_template_fit_table_dir = carbon_template_fit_dir / "tables"
    carbon_template_fit_plot_dir = carbon_template_fit_dir / "plots"
    final_carbon_cut_dir = output_dir / "final_carbon_assisted_cuts"
    for path in (
        output_dir,
        table_dir,
        plot_dir,
        spectrum_dir,
        model_canvas_dir,
        summary_dir,
        pull_diagnostic_dir,
        carbon_phase1_dir,
        carbon_template_fit_dir,
        carbon_template_fit_table_dir,
        carbon_template_fit_plot_dir,
        final_carbon_cut_dir,
    ):
        ensure_directory(path)
    # endfor

    regular_fit_histogram_bins = (
        args.hist_bins // args.regular_fit_bin_factor
    )
    jobs, loading_metadata = build_fit_jobs(
        datasets=datasets,
        histogram_min_gev2=args.hist_min,
        histogram_max_gev2=args.hist_max,
        histogram_bins=regular_fit_histogram_bins,
        fit_min_gev2=args.fit_min,
        fit_max_gev2=args.fit_max,
        minimum_events=args.minimum_events,
        corrected_mean_max_gev2=args.corrected_mean_max,
        step_size=args.step_size,
    )

    carbon_nh3_jobs: list[FitJob] = []
    carbon_jobs: list[FitJob] = []
    carbon_loading_metadata: dict[str, Any] = {}
    if not args.disable_carbon_phase1:
        carbon_nh3_jobs, _ = build_fit_jobs(
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
        carbon_jobs, carbon_loading_metadata = build_fit_jobs(
            datasets=carbon_datasets,
            histogram_min_gev2=args.hist_min,
            histogram_max_gev2=args.hist_max,
            histogram_bins=args.hist_bins,
            fit_min_gev2=args.fit_min,
            fit_max_gev2=args.fit_max,
            minimum_events=1,
            corrected_mean_max_gev2=args.corrected_mean_max,
            step_size=args.step_size,
        )
    # endif

    joint_config = JointFitConfig(
        aicc_improvement_required=args.joint_aicc_improvement,
        local_improvement_fraction=(
            args.joint_local_improvement_fraction
        ),
        accepted_local_score_max=args.accepted_local_score_max,
        marginal_local_score_max=args.marginal_local_score_max,
        sideband_exclusion_min_gev2=args.sideband_exclusion_min,
        sideband_exclusion_max_gev2=args.sideband_exclusion_max,
        range_study_maxima_gev2=(1.30, 1.35, 1.40),
        enable_neighbor_mean_fallback=(
            not args.disable_neighbor_mean_fallback
        ),
        replicas=args.replicas,
        replica_seed=args.replica_seed,
    )
    results = execute_fit_jobs(
        jobs=jobs,
        workers=args.workers,
        joint_config=joint_config,
    )
    frame = flatten_results(results)

    csv_path = table_dir / "mx2_peak_fit_results_v27.csv"
    json_path = table_dir / "mx2_peak_fit_results_v27.json"
    latex_path = table_dir / "mx2_peak_nominal_corrected_tables_v27.tex"
    status_path = table_dir / "mx2_peak_fit_status_v27.txt"

    frame.to_csv(csv_path, index=False)
    compact_json = build_compact_json(
        results=results,
        datasets=datasets,
        loading_metadata=loading_metadata,
        args=args,
    )
    compact_json["sample_variant"] = sample_variant
    compact_json["input_paths"] = {
        period: {stage: str(path) for stage, path in stage_map.items()}
        for period, stage_map in input_paths.items()
    }
    compact_json["carbon_input_paths"] = {
        period: {stage: str(path) for stage, path in stage_map.items()}
        for period, stage_map in carbon_input_paths.items()
    }
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

    carbon_products: dict[str, Any] | None = None
    carbon_template_products: dict[str, Any] | None = None
    final_carbon_cut_products: dict[str, Any] | None = None
    if not args.disable_carbon_phase1:
        # Nominal: one normalization per period and correction stage.
        carbon_period_diagnostics = build_carbon_phase1_diagnostics(
            nh3_jobs=carbon_nh3_jobs,
            carbon_jobs=carbon_jobs,
            grouping="period",
            control_min=args.carbon_control_min,
            control_max=args.carbon_control_max,
            scan_boundaries=args.carbon_control_boundaries,
            hist_min=args.hist_min,
            hist_max=args.hist_max,
            hist_bins=args.hist_bins,
            validation_low=(
                args.carbon_validation_low_min,
                args.carbon_validation_low_max,
            ),
            validation_high=(
                args.carbon_validation_high_min,
                args.carbon_validation_high_max,
            ),
        )
        # Cross-check retained for demonstrating stability versus -tprime.
        carbon_column_diagnostics = build_carbon_phase1_diagnostics(
            nh3_jobs=carbon_nh3_jobs,
            carbon_jobs=carbon_jobs,
            grouping="tprime-column",
            control_min=args.carbon_control_min,
            control_max=args.carbon_control_max,
            scan_boundaries=args.carbon_control_boundaries,
            hist_min=args.hist_min,
            hist_max=args.hist_max,
            hist_bins=args.hist_bins,
            validation_low=(
                args.carbon_validation_low_min,
                args.carbon_validation_low_max,
            ),
            validation_high=(
                args.carbon_validation_high_min,
                args.carbon_validation_high_max,
            ),
        )
        metadata = {
            "carbon_loading_metadata": carbon_loading_metadata,
            "carbon_inputs": {
                period: {
                    stage: str(path)
                    for stage, path in stage_map.items()
                }
                for period, stage_map in (
                    carbon_input_paths or {}
                ).items()
            },
        }
        carbon_period_diagnostics["configuration"].update(metadata)
        carbon_column_diagnostics["configuration"].update(metadata)

        carbon_products = write_carbon_phase1_products(
            nominal_diagnostics=carbon_period_diagnostics,
            column_diagnostics=carbon_column_diagnostics,
            output_dir=carbon_phase1_dir,
        )

        if not args.disable_carbon_template_fits:
            carbon_template_results = run_carbon_template_fits(
                nh3_jobs=jobs,
                carbon_jobs=carbon_jobs,
                period_diagnostics=carbon_period_diagnostics,
                fit_min_gev2=args.carbon_template_fit_min,
                fit_max_gev2=args.carbon_template_fit_max,
                control_min_gev2=args.carbon_control_min,
                control_endpoints_gev2=(
                    args.carbon_template_control_endpoints
                ),
            )
            carbon_template_products = (
                write_carbon_template_fit_products(
                    results=carbon_template_results,
                    production_frame=frame,
                    table_dir=carbon_template_fit_table_dir,
                    plot_dir=carbon_template_fit_plot_dir,
                )
            )
            final_carbon_cut_products = (
                write_final_carbon_cut_products(
                    results=carbon_template_results,
                    production_frame=frame,
                    output_dir=final_carbon_cut_dir,
                )
            )
        # endif
    # endif

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
    print(f"{sample_variant.upper()} deterministic Mx2 fit study complete.")
    print(f"Output directory: {output_dir}")
    print(f"Flat fit table:    {csv_path}")
    print(f"Compact JSON:      {json_path}")
    print(f"LaTeX tables:      {latex_path}")
    print(f"Status report:     {status_path}")
    print(
        "Regular-fit bins:   "
        f"{regular_fit_histogram_bins} "
        f"(carbon binning: {args.hist_bins})"
    )
    if carbon_products is not None:
        print(f"Carbon Phase-1:    {carbon_products['json']}")
        print(
            "Carbon normalization: period-wide nominal; "
            "-tprime-column stability cross-check; "
            f"{args.carbon_control_min:.3f} <= Mx2 < "
            f"{args.carbon_control_max:.3f} GeV^2"
        )
    # endif
    if carbon_template_products is not None:
        print(
            "Carbon template fits: "
            f"{carbon_template_products['json']}"
        )
    # endif
    if final_carbon_cut_products is not None:
        print(
            "Final carbon cuts:  "
            f"{final_carbon_cut_products['csv']}"
        )
        print(
            "Method comparison:  "
            f"{final_carbon_cut_products['comparison_plot']}"
        )
        print(
            "Cut summary:        "
            f"{final_carbon_cut_products['master_summary_plot']}"
        )
    # endif
    print(
        "Corrected recommended fits: "
        f"{successful_recommended_after}/72 successful; "
        f"{accepted_recommended_after} accepted; "
        f"{unresolved_recommended_after} unresolved."
    )

    variant_manifest_path = output_dir / "analysis_variant_manifest.json"
    variant_manifest = {
        "schema_version": 1,
        "sample_variant": sample_variant,
        "input_paths": {
            period: {stage: str(path) for stage, path in stage_map.items()}
            for period, stage_map in input_paths.items()
        },
        "carbon_input_paths": {
            period: {stage: str(path) for stage, path in stage_map.items()}
            for period, stage_map in carbon_input_paths.items()
        },
        "fit_table_csv": str(csv_path),
        "compact_fit_json": str(json_path),
        "final_cut_csv": (
            str(final_carbon_cut_products["csv"])
            if final_carbon_cut_products is not None else None
        ),
        "final_cut_json": (
            str(final_carbon_cut_products["json"])
            if final_carbon_cut_products is not None else None
        ),
        "python_cut_lookup": (
            str(final_carbon_cut_products["python_lookup"])
            if final_carbon_cut_products is not None else None
        ),
        "cpp_cut_lookup": (
            str(final_carbon_cut_products["cpp_lookup"])
            if final_carbon_cut_products is not None else None
        ),
        "final_cut_metadata": (
            str(final_carbon_cut_products["metadata"])
            if final_carbon_cut_products is not None else None
        ),
    }
    variant_manifest_path.write_text(
        json.dumps(json_safe(variant_manifest), indent=2) + "\n",
        encoding="utf-8",
    )
    return {
        "sample_variant": sample_variant,
        "output_dir": str(output_dir),
        "fit_table_csv": str(csv_path),
        "compact_fit_json": str(json_path),
        "variant_manifest": str(variant_manifest_path),
        "final_cut_csv": variant_manifest["final_cut_csv"],
        "final_cut_json": variant_manifest["final_cut_json"],
        "python_cut_lookup": variant_manifest["python_cut_lookup"],
        "cpp_cut_lookup": variant_manifest["cpp_cut_lookup"],
        "final_cut_metadata": variant_manifest["final_cut_metadata"],
    }



def clone_input_map(
    defaults: dict[str, dict[str, Path]],
) -> dict[str, dict[str, Path]]:
    return {
        period: {stage: Path(path) for stage, path in stage_map.items()}
        for period, stage_map in defaults.items()
    }


def apply_input_overrides(
    input_paths: dict[str, dict[str, Path]],
    overrides: list[tuple[str, str, Path]],
) -> None:
    for period, stage, path in overrides:
        input_paths[period][stage] = path
    # endfor


def copy_nominal_production_products(
    nominal_result: dict[str, Any],
    root_output: Path,
) -> dict[str, str]:
    """
    Publish stable, version-independent nominal cut products.

    Downstream analyses consume these canonical paths directly. Script version
    numbers are intentionally absent from all production filenames.
    """
    destination_root = root_output / "final_carbon_assisted_cuts"
    tables_dir = destination_root / "tables"
    exports_dir = destination_root / "exports"
    ensure_directory(tables_dir)
    ensure_directory(exports_dir)

    source_csv = Path(nominal_result["final_cut_csv"]).resolve()
    source_json = Path(nominal_result["final_cut_json"]).resolve()
    source_python = Path(nominal_result["python_cut_lookup"]).resolve()
    source_cpp = Path(nominal_result["cpp_cut_lookup"]).resolve()
    source_metadata = Path(
        nominal_result["final_cut_metadata"]
    ).resolve()

    canonical_csv = (
        tables_dir / "final_carbon_assisted_mx2_cuts.csv"
    )
    canonical_json = (
        tables_dir / "final_carbon_assisted_mx2_cuts.json"
    )
    canonical_python = (
        exports_dir / "final_carbon_assisted_mx2_cuts.py"
    )
    canonical_cpp = (
        exports_dir / "final_carbon_assisted_mx2_cuts.h"
    )
    canonical_metadata = (
        exports_dir / "final_carbon_assisted_mx2_cuts_metadata.json"
    )

    for source, destination in (
        (source_csv, canonical_csv),
        (source_json, canonical_json),
        (source_python, canonical_python),
        (source_cpp, canonical_cpp),
        (source_metadata, canonical_metadata),
    ):
        destination_resolved = destination.resolve()
        if source != destination_resolved:
            shutil.copy2(source, destination_resolved)
        # endif
    # endfor

    return {
        "csv": str(canonical_csv),
        "json": str(canonical_json),
        "python_lookup": str(canonical_python),
        "cpp_lookup": str(canonical_cpp),
        "metadata": str(canonical_metadata),
    }


def main() -> int:
    """Derive production cuts from nominal data and retain ISR only as a diagnostic."""
    parser = build_argument_parser()
    args = parser.parse_args()

    try:
        validate_arguments(args)
        root_output = args.output_dir.expanduser().resolve()
        nominal_work = root_output / "nominal"
        isr_work = root_output / "isr"
        diagnostics_output = root_output / "diagnostics"

        ensure_directory(root_output)
        for analysis_directory in (
            nominal_work,
            isr_work,
            diagnostics_output,
        ):
            if analysis_directory.exists():
                shutil.rmtree(analysis_directory)
            # endif
        # endfor

        nominal_inputs = clone_input_map(DEFAULT_INPUTS)
        nominal_carbon_inputs = clone_input_map(DEFAULT_CARBON_INPUTS)
        apply_input_overrides(nominal_inputs, args.input)
        apply_input_overrides(nominal_carbon_inputs, args.carbon_input)

        nominal_result = run_analysis_variant(
            args,
            "nominal",
            nominal_inputs,
            nominal_carbon_inputs,
            nominal_work,
        )
        if nominal_result["final_cut_json"] is None:
            raise RuntimeError(
                "Nominal carbon-assisted cuts were not produced. Do not use "
                "--disable-carbon-phase1 or --disable-carbon-template-fits for "
                "the production channel-selection run."
            )
        # endif

        production_products = copy_nominal_production_products(
            nominal_result,
            root_output,
        )

        comparison_products: dict[str, Any] | None = None
        if not args.disable_isr:
            isr_inputs = clone_input_map(DEFAULT_ISR_INPUTS)
            isr_carbon_inputs = clone_input_map(DEFAULT_ISR_CARBON_INPUTS)
            apply_input_overrides(isr_inputs, args.isr_input)
            apply_input_overrides(isr_carbon_inputs, args.isr_carbon_input)
            isr_result = run_analysis_variant(
                args,
                "isr",
                isr_inputs,
                isr_carbon_inputs,
                isr_work,
            )
            if isr_result["final_cut_csv"] is None:
                raise RuntimeError(
                    "ISR diagnostic fits did not produce a final cut table."
                )
            # endif
            comparison_products = write_nominal_isr_comparison_products(
                nominal_result,
                isr_result,
                diagnostics_output,
            )
        # endif

        manifest_path = root_output / "channel_selection_manifest.json"
        manifest = {
            "schema_version": 2,
            "production_policy": (
                "All downstream selections use the non-ISR nominal "
                "carbon-assisted cuts. ISR fits are diagnostic only."
            ),
            "production_products": production_products,
            "retained_analysis_products": {
                "nominal": str(nominal_work),
                "isr": (
                    str(isr_work) if not args.disable_isr else None
                ),
                "nominal_isr_comparisons": (
                    str(diagnostics_output)
                    if not args.disable_isr else None
                ),
            },
            "diagnostics": (
                {key: str(value) for key, value in comparison_products.items()}
                if comparison_products is not None else None
            ),
            "dilution_factor_contract": {
                "compatible": True,
                "filenames_are_version_independent": True,
                "default_json_expected_by_determine_dilution_factor": (
                    production_products["json"]
                ),
                "required_json_schema": (
                    "periods.<period>[] with x_index, t_index, bin edges, "
                    "mu_gev2, mu_error_gev2, sigma_gev2, "
                    "sigma_error_gev2, tight, nominal, loose, and "
                    "classification"
                ),
            },
        }
        ensure_directory(root_output)
        manifest_path.write_text(
            json.dumps(json_safe(manifest), indent=2) + "\n",
            encoding="utf-8",
        )

        print("")
        print("Nominal production channel-selection analysis complete.")
        print(f"Production cut JSON: {production_products['json']}")
        print(
            "Stable downstream cut JSON: "
            f"{production_products['json']}"
        )
        print(f"Complete nominal products: {nominal_work}")
        if not args.disable_isr:
            print(f"Complete ISR products:     {isr_work}")
        # endif
        if comparison_products is not None:
            print(
                "Nominal/ISR diagnostic plot: "
                f"{comparison_products['mu_sigma_plot']}"
            )
            print(
                "Nominal/ISR difference plot: "
                f"{comparison_products['difference_plot']}"
            )
        # endif
        print(f"Manifest: {manifest_path}")
        return 0
    except Exception as exc:
        print(f"FATAL ERROR: {exc}", file=sys.stderr)
        return 1
    # endtry


if __name__ == "__main__":
    raise SystemExit(main())
# endif

