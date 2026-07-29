#!/usr/bin/env python3
"""
plot_calibration_v14.py

Fast calibration diagnostics for the CLAS12 DVCS calibration trees.

Current modules
---------------
1. Forward Tagger photon occupancy
   - One 2x2 summary per run period.
   - Data before/after the FT fiducial cut on the top row.
   - Reconstructed MC before/after on the bottom row.

2. Electromagnetic-calorimeter strip matching
   - Separate output directories for PCal, ECin, and ECout.
   - Separate electron and photon summaries.
   - Six sectors by three strip coordinates: lu, lv, and lw.
   - One-dimensional data/MC overlays use 4.5 cm strip-width bins and
     logarithmic vertical axes.
   - Original-style two-dimensional lu-lv, lu-lw, and lw-lv occupancy
     maps are also produced with logarithmic color scales. Raw counts are
     used, with independent data and MC scales.
   - Existing fiducial exclusion intervals are shaded but are NOT applied to
     the plotted samples, so data/GEMC mismodeling remains visible.
   - The period-dependent PCal sector-2 lv exclusions are overlaid only for the applicable RGA periods.
   - PCal-electron mean momentum and mean scattering-angle maps are made
     versus lu-lv, lu-lw, and lv-lw before and after the cuts.
   - Sector-resolved one-dimensional mean momentum and mean theta curves
     are also made versus lu, lv, and lw.

Performance model
-----------------
One worker process handles one complete run period.  Each ROOT file is read
once, and the FT and calorimeter histograms are accumulated in the same chunk
loop.  At most five period workers are used.

Default RGA input
-----------------
/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/dvcs_calibration/

Default RGC input
-----------------
/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/

Default outputs
---------------
output/calibration/forward_tagger/
output/calibration/calorimeter/before_cuts/pcal/
output/calibration/calorimeter/before_cuts/ecin/
output/calibration/calorimeter/before_cuts/ecout/
output/calibration/calorimeter/after_cuts/pcal/
output/calibration/calorimeter/after_cuts/ecin/
output/calibration/calorimeter/after_cuts/ecout/

Examples
--------
All five RGA periods:

  python3 external_scripts/plot_calibration_v14.py

RGC data:

  python3 external_scripts/plot_calibration_v14.py --rgc

One period, limited test:

  python3 external_scripts/plot_calibration_v14.py \
      --period rga_sp18_inb \
      --max-events 5000000 \
      --workers 1

Only calorimeter plots:

  python3 external_scripts/plot_calibration_v14.py --skip-ft

Only FT plots:

  python3 external_scripts/plot_calibration_v14.py --skip-calorimeter
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Circle
import numpy as np

try:
    import uproot
except ImportError as exc:
    raise SystemExit(
        "ERROR: uproot is required. Install it with "
        "`python3 -m pip install uproot`."
    ) from exc


# =============================================================================
# Global configuration
# =============================================================================

TREE_NAME_DEFAULT = "PhysicsEvents"

RGA_INPUT_DIR_DEFAULT = Path(
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/dvcs_calibration"
)
RGC_INPUT_DIR_DEFAULT = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration"
)
OUTPUT_BASE_DEFAULT = Path("output/calibration")

PHOTON_PID = 22
ELECTRON_PID = 11
INVALID_SENTINEL = -9999.0
MAX_WORKERS_HARD_LIMIT = 5

PARTICLE_KEYS = ("electron", "photon")
PARTICLE_PIDS = np.asarray((ELECTRON_PID, PHOTON_PID), dtype=np.int32)

LAYER_KEYS = ("pcal", "ecin", "ecout")
LAYER_NUMBERS = (1, 4, 7)
LAYER_LABELS = {
    "pcal": "PCal",
    "ecin": "ECin",
    "ecout": "ECout",
}

COORD_KEYS = ("lu", "lv", "lw")
COORD_LABELS = {
    "lu": r"$lu$ (cm)",
    "lv": r"$lv$ (cm)",
    "lw": r"$lw$ (cm)",
}

PAIR_DEFINITIONS: tuple[tuple[str, str, str], ...] = (
    ("lu_lv", "lu", "lv"),
    ("lu_lw", "lu", "lw"),
    ("lw_lv", "lw", "lv"),
)

# PCal electron kinematic maps use all three coordinate-pair projections.
# Values are indices into PAIR_DEFINITIONS.
KINEMATIC_PAIR_INDICES: tuple[int, ...] = (0, 1, 2)

KINEMATIC_KEYS: tuple[str, ...] = ("momentum", "theta")
KINEMATIC_LABELS = {
    "momentum": "Mean electron momentum (GeV)",
    "theta": "Mean electron theta (deg)",
}
KINEMATIC_FILENAME_KEYS = {
    "momentum": "mean_momentum",
    "theta": "mean_theta",
}


@dataclass(frozen=True)
class Dataset:
    key: str
    label: str
    data_file: str
    mc_file: str | None = None


@dataclass(frozen=True)
class FTFiducialConfig:
    inner_radius_cm: float = 8.5
    outer_radius_cm: float = 15.5
    holes: tuple[tuple[float, float, float], ...] = (
        (-8.42, 9.89, 1.60),
        (-9.89, -5.33, 1.60),
        (-6.15, -13.00, 2.30),
        (3.70, -6.50, 2.00),
    )


@dataclass(frozen=True)
class HistogramConfig:
    ft_bins: int = 100
    ft_xy_min_cm: float = -20.0
    ft_xy_max_cm: float = 20.0
    cal_bin_width_cm: float = 4.5
    cal_min_cm: float = 0.0
    cal_max_cm: float = 450.0

    @property
    def ft_edges(self) -> np.ndarray:
        return np.linspace(
            self.ft_xy_min_cm,
            self.ft_xy_max_cm,
            self.ft_bins + 1,
            dtype=np.float64,
        )

    @property
    def cal_edges(self) -> np.ndarray:
        span = self.cal_max_cm - self.cal_min_cm
        bins_float = span / self.cal_bin_width_cm
        bins = int(round(bins_float))
        if not np.isclose(bins_float, bins, rtol=0.0, atol=1.0e-10):
            raise ValueError(
                "Calorimeter range must be an integer multiple of "
                "the requested strip width."
            )
        #endif
        return self.cal_min_cm + self.cal_bin_width_cm * np.arange(
            bins + 1,
            dtype=np.float64,
        )

    @property
    def cal_bins(self) -> int:
        return self.cal_edges.size - 1


@dataclass
class SampleResult:
    ft_before: np.ndarray
    ft_after: np.ndarray
    cal_counts_before: np.ndarray
    cal_counts_after: np.ndarray
    cal_pair_counts_before: np.ndarray
    cal_pair_counts_after: np.ndarray
    pcal_electron_kin_counts_before: np.ndarray
    pcal_electron_kin_counts_after: np.ndarray
    pcal_electron_kin_sums_before: np.ndarray
    pcal_electron_kin_sums_after: np.ndarray
    pcal_electron_kin_1d_counts_before: np.ndarray
    pcal_electron_kin_1d_counts_after: np.ndarray
    pcal_electron_kin_1d_sums_before: np.ndarray
    pcal_electron_kin_1d_sums_after: np.ndarray
    rows_read: int
    photon_rows: int
    valid_ft_photons: int
    fiducial_ft_photons: int
    elapsed_seconds: float


@dataclass
class PeriodResult:
    dataset: Dataset
    data: SampleResult
    mc: SampleResult | None
    worker_pid: int
    elapsed_seconds: float


RGA_DATASETS: tuple[Dataset, ...] = (
    Dataset(
        "rga_sp18_inb",
        "RGA Sp18 Inb",
        "DVCSWagon_rga_sp18_inb.root",
        "dvcsgen_rga_sp18_inb.root",
    ),
    Dataset(
        "rga_sp18_out",
        "RGA Sp18 Out",
        "DVCSWagon_rga_sp18_out.root",
        "dvcsgen_rga_sp18_out.root",
    ),
    Dataset(
        "rga_fa18_inb",
        "RGA Fa18 Inb",
        "DVCSWagon_rga_fa18_inb.root",
        "dvcsgen_rga_fa18_inb.root",
    ),
    Dataset(
        "rga_fa18_out",
        "RGA Fa18 Out",
        "DVCSWagon_rga_fa18_out.root",
        "dvcsgen_rga_fa18_out.root",
    ),
    Dataset(
        "rga_sp19_inb",
        "RGA Sp19 Inb",
        "DVCSWagon_rga_sp19_inb.root",
        "dvcsgen_rga_sp19_inb.root",
    ),
)

RGC_DATASETS: tuple[Dataset, ...] = (
    Dataset(
        "rgc_su22_inb",
        "RGC Su22 Inb",
        "rgc_su22_inb_NH3_epi+X_calibration.root",
    ),
    Dataset(
        "rgc_fa22_inb",
        "RGC Fa22 Inb",
        "rgc_fa22_inb_NH3_epi+X_calibration.root",
    ),
    Dataset(
        "rgc_sp23_inb",
        "RGC Sp23 Inb",
        "rgc_sp23_inb_NH3_epi+X_calibration.root",
    ),
)


# =============================================================================
# Existing calorimeter exclusions
# =============================================================================

# Each tuple is:
#   (layer key, sector, coordinate key, lower bound, upper bound,
#    applicable dataset keys)
#
# The generic PCal edge requirement is handled separately below.  These
# intervals describe only the localized dead-channel exclusions used for
# data/MC matching.
ALL_RGA_DATASET_KEYS: frozenset[str] = frozenset(
    dataset.key for dataset in RGA_DATASETS
)
SP18_DATASET_KEYS: frozenset[str] = frozenset(
    ("rga_sp18_inb", "rga_sp18_out")
)
FA18_DATASET_KEYS: frozenset[str] = frozenset(
    ("rga_fa18_inb", "rga_fa18_out")
)
SP19_DATASET_KEYS: frozenset[str] = frozenset(("rga_sp19_inb",))

SP18_AND_SP19_DATASET_KEYS: frozenset[str] = (
    SP18_DATASET_KEYS | SP19_DATASET_KEYS
)
FA18_AND_SP19_DATASET_KEYS: frozenset[str] = (
    FA18_DATASET_KEYS | SP19_DATASET_KEYS
)

DEAD_STRIP_EXCLUSIONS: tuple[
    tuple[str, int, str, float, float, frozenset[str]], ...
] = (
    # PCal exclusions common to all five RGA periods.
    ("pcal", 1, "lw", 72.0, 94.5, ALL_RGA_DATASET_KEYS),
    ("pcal", 1, "lw", 220.5, 234.0, ALL_RGA_DATASET_KEYS),
    ("pcal", 3, "lv", 346.5, 378.0, ALL_RGA_DATASET_KEYS),
    ("pcal", 4, "lv", 229.5, 243.0, ALL_RGA_DATASET_KEYS),
    ("pcal", 6, "lw", 166.5, 193.5, ALL_RGA_DATASET_KEYS),

    # Period-dependent PCal exclusions.
    ("pcal", 1, "lw", 328.5, 342.0, SP18_DATASET_KEYS),
    ("pcal", 2, "lv", 31.5, 49.5, SP18_AND_SP19_DATASET_KEYS),
    ("pcal", 2, "lv", 99.0, 117.0, FA18_AND_SP19_DATASET_KEYS),

    # ECin exclusions common to all five RGA periods.
    ("ecin", 1, "lv", 67.5, 94.5, ALL_RGA_DATASET_KEYS),
    ("ecin", 5, "lv", 0.0, 23.5, ALL_RGA_DATASET_KEYS),

    # ECout exclusions common to all five RGA periods.
    ("ecout", 1, "lv", 0.0, 67.5, ALL_RGA_DATASET_KEYS),
    ("ecout", 5, "lu", 198.0, 220.5, ALL_RGA_DATASET_KEYS),

    # Spring-2018-specific ECout exclusion.
    ("ecout", 2, "lw", 306.0, 328.5, SP18_DATASET_KEYS),
    ("ecout", 3, "lu", 360.0, 373.5, SP18_DATASET_KEYS),
    ("ecout", 4, "lw", 171.0, 189.0, SP18_DATASET_KEYS),
)


def pcal_edge_threshold_cm(strictness: int) -> float:
    if strictness == 1:
        return 9.0
    #endif
    if strictness == 2:
        return 13.5
    #endif
    if strictness == 3:
        return 18.0
    #endif
    raise ValueError(f"Unsupported calorimeter strictness: {strictness}")


def exclusion_intervals(
    dataset_key: str,
    layer_key: str,
    sector: int,
    coordinate_key: str,
    strictness: int,
    is_rgc: bool,
) -> list[tuple[float, float, str]]:
    """
    Return currently implemented exclusion intervals for one panel.

    The generic PCal edge cut is included for lv and lw.  The RGA dead-strip
    removals are included only for strictness >= 2, matching the Java logic.
    """

    intervals: list[tuple[float, float, str]] = []

    if layer_key == "pcal" and coordinate_key in ("lv", "lw"):
        intervals.append(
            (0.0, pcal_edge_threshold_cm(strictness), "PCal edge")
        )
    #endif

    if strictness < 2 or is_rgc:
        return intervals
    #endif

    for (
        cut_layer,
        cut_sector,
        cut_coordinate,
        lower,
        upper,
        applicable_dataset_keys,
    ) in DEAD_STRIP_EXCLUSIONS:
        if cut_layer != layer_key:
            continue
        #endif
        if cut_sector != sector:
            continue
        #endif
        if cut_coordinate != coordinate_key:
            continue
        #endif
        if dataset_key not in applicable_dataset_keys:
            continue
        #endif

        intervals.append((lower, upper, "RGA period-specific exclusion"))
    #endfor

    return intervals


def calorimeter_fiducial_mask(
    particle_mask: np.ndarray,
    sectors: np.ndarray,
    coordinate_arrays: dict[tuple[str, str], np.ndarray],
    dataset_key: str,
    strictness: int,
    is_rgc: bool,
) -> np.ndarray:
    """
    Apply the complete particle-level calorimeter fiducial decision.

    This reproduces the corrected Java logic:
      * PCal lv and lw must exceed the strictness-dependent edge threshold.
      * For strictness >= 2, applicable RGA dead-strip intervals are removed.
      * PCal sector-2 exclusions are selected explicitly by run period.

    The result is calculated once per particle species and ROOT chunk, then
    reused for every layer, sector, coordinate, and 2D coordinate pair.
    """

    accepted = particle_mask.copy()
    threshold_cm = pcal_edge_threshold_cm(strictness)

    pcal_lv = coordinate_arrays[("pcal", "lv")]
    pcal_lw = coordinate_arrays[("pcal", "lw")]

    accepted &= np.isfinite(pcal_lv)
    accepted &= np.isfinite(pcal_lw)
    accepted &= pcal_lv != INVALID_SENTINEL
    accepted &= pcal_lw != INVALID_SENTINEL
    accepted &= pcal_lv >= threshold_cm
    accepted &= pcal_lw >= threshold_cm

    if strictness < 2 or is_rgc:
        return accepted
    #endif

    for (
        layer_key,
        sector,
        coordinate_key,
        lower,
        upper,
        applicable_dataset_keys,
    ) in DEAD_STRIP_EXCLUSIONS:
        if dataset_key not in applicable_dataset_keys:
            continue
        #endif

        values = coordinate_arrays[(layer_key, coordinate_key)]
        in_interval = (
            (sectors == sector)
            & np.isfinite(values)
            & (values != INVALID_SENTINEL)
            & (values > lower)
            & (values < upper)
        )
        accepted &= ~in_interval
    #endfor

    return accepted


# =============================================================================
# Command line
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create FT occupancy and calorimeter strip-matching diagnostics "
            "for the configured RGA or RGC calibration trees."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--rga",
        action="store_true",
        help="Process the five RGA periods; this is the default.",
    )
    mode.add_argument(
        "--rgc",
        action="store_true",
        help="Process the three RGC data periods.",
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Override the selected mode's input directory.",
    )
    parser.add_argument(
        "--output-base",
        type=Path,
        default=OUTPUT_BASE_DEFAULT,
        help="Base directory under which calibration modules are written.",
    )
    parser.add_argument(
        "--tree",
        default=TREE_NAME_DEFAULT,
        help="ROOT TTree name.",
    )
    parser.add_argument(
        "--period",
        action="append",
        default=[],
        metavar="KEY",
        help="Select one period key; may be repeated.",
    )
    parser.add_argument(
        "--list-periods",
        action="store_true",
        help="List valid period keys and exit.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=MAX_WORKERS_HARD_LIMIT,
        help="Requested period workers; always capped at five.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Rows per uproot chunk.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional per-file row limit for tests.",
    )
    parser.add_argument(
        "--skip-ft",
        action="store_true",
        help="Do not create FT plots or read FT branches.",
    )
    parser.add_argument(
        "--skip-calorimeter",
        action="store_true",
        help="Do not create calorimeter plots or read calorimeter branches.",
    )
    parser.add_argument(
        "--ft-bins",
        type=int,
        default=100,
        help="FT bins per coordinate.",
    )
    parser.add_argument(
        "--ft-xy-min",
        type=float,
        default=-20.0,
        help="Minimum FT coordinate in cm.",
    )
    parser.add_argument(
        "--ft-xy-max",
        type=float,
        default=20.0,
        help="Maximum FT coordinate in cm.",
    )
    parser.add_argument(
        "--cal-bin-width",
        type=float,
        default=4.5,
        help="Calorimeter strip-coordinate bin width in cm.",
    )
    parser.add_argument(
        "--cal-min",
        type=float,
        default=0.0,
        help="Minimum calorimeter coordinate in cm.",
    )
    parser.add_argument(
        "--cal-max",
        type=float,
        default=450.0,
        help="Maximum calorimeter coordinate in cm.",
    )
    parser.add_argument(
        "--calorimeter-strictness",
        type=int,
        choices=(1, 2, 3),
        default=2,
        help=(
            "Strictness whose existing PCal edge and dead-strip exclusions "
            "are overlaid. Histograms themselves remain uncut."
        ),
    )
    parser.add_argument(
        "--kinematic-min-bin-count",
        type=int,
        default=10,
        help=(
            "Minimum entries required to display a PCal electron mean-"
            "kinematics bin."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="PNG resolution.",
    )

    args = parser.parse_args()

    if args.skip_ft and args.skip_calorimeter:
        parser.error("--skip-ft and --skip-calorimeter cannot both be set")
    #endif
    if args.workers <= 0:
        parser.error("--workers must be positive")
    #endif
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")
    #endif
    if args.max_events is not None and args.max_events <= 0:
        parser.error("--max-events must be positive when supplied")
    #endif
    if args.ft_bins <= 0:
        parser.error("--ft-bins must be positive")
    #endif
    if args.cal_bin_width <= 0.0:
        parser.error("--cal-bin-width must be positive")
    #endif
    if args.ft_xy_max <= args.ft_xy_min:
        parser.error("--ft-xy-max must exceed --ft-xy-min")
    #endif
    if args.cal_max <= args.cal_min:
        parser.error("--cal-max must exceed --cal-min")
    #endif
    if args.kinematic_min_bin_count <= 0:
        parser.error("--kinematic-min-bin-count must be positive")
    #endif
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    #endif

    return args


def choose_datasets(
    args: argparse.Namespace,
) -> tuple[tuple[Dataset, ...], Path, str]:
    is_rgc = bool(args.rgc)
    datasets = RGC_DATASETS if is_rgc else RGA_DATASETS
    input_dir = (
        args.input_dir
        if args.input_dir is not None
        else (RGC_INPUT_DIR_DEFAULT if is_rgc else RGA_INPUT_DIR_DEFAULT)
    )
    mode_label = "RGC" if is_rgc else "RGA"

    if args.period:
        by_key = {dataset.key: dataset for dataset in datasets}
        unknown = sorted(set(args.period) - set(by_key))
        if unknown:
            raise ValueError(
                f"Unknown {mode_label} period key(s): {', '.join(unknown)}. "
                f"Valid keys: {', '.join(by_key)}"
            )
        #endif

        requested = set(args.period)
        datasets = tuple(
            dataset for dataset in datasets if dataset.key in requested
        )
    #endif

    return datasets, input_dir, mode_label


# =============================================================================
# ROOT reading and accumulation
# =============================================================================


def required_branches(
    do_ft: bool,
    do_calorimeter: bool,
) -> tuple[str, ...]:
    branches: list[str] = ["particle_pid"]

    if do_ft:
        branches.extend(("ft_x", "ft_y"))
    #endif

    if do_calorimeter:
        branches.extend(("p", "theta", "cal_sector"))
        for layer_number in LAYER_NUMBERS:
            for coordinate_key in COORD_KEYS:
                branches.append(f"cal_{coordinate_key}_{layer_number}")
            #endfor
        #endfor
    #endif

    return tuple(branches)


def validate_root_tree(
    file_path: Path,
    tree_name: str,
    branches: tuple[str, ...],
) -> int:
    if not file_path.is_file():
        raise FileNotFoundError(f"Input ROOT file does not exist: {file_path}")
    #endif

    with uproot.open(file_path) as root_file:
        if tree_name not in root_file:
            keys = ", ".join(str(key) for key in root_file.keys())
            raise KeyError(
                f"Tree '{tree_name}' is absent from {file_path}. "
                f"Available keys: {keys}"
            )
        #endif

        tree = root_file[tree_name]
        available = set(tree.keys())
        missing = [branch for branch in branches if branch not in available]
        if missing:
            raise KeyError(
                f"Tree '{tree_name}' in {file_path} is missing: "
                f"{', '.join(missing)}"
            )
        #endif

        return int(tree.num_entries)


def ft_fiducial_mask(
    x_cm: np.ndarray,
    y_cm: np.ndarray,
    config: FTFiducialConfig,
) -> np.ndarray:
    radius_squared = x_cm * x_cm + y_cm * y_cm
    accepted = (
        (radius_squared >= config.inner_radius_cm**2)
        & (radius_squared <= config.outer_radius_cm**2)
    )

    for center_x_cm, center_y_cm, hole_radius_cm in config.holes:
        delta_x = x_cm - center_x_cm
        delta_y = y_cm - center_y_cm
        accepted &= (
            delta_x * delta_x + delta_y * delta_y
            >= hole_radius_cm**2
        )
    #endfor

    return accepted


def empty_cal_counts(histogram_config: HistogramConfig) -> np.ndarray:
    # Axes:
    #   particle, layer, sector, coordinate, histogram bin
    return np.zeros(
        (
            len(PARTICLE_KEYS),
            len(LAYER_KEYS),
            6,
            len(COORD_KEYS),
            histogram_config.cal_bins,
        ),
        dtype=np.uint64,
    )


def empty_cal_pair_counts(
    histogram_config: HistogramConfig,
) -> np.ndarray:
    # Axes:
    #   particle, layer, sector, coordinate pair, x bin, y bin
    return np.zeros(
        (
            len(PARTICLE_KEYS),
            len(LAYER_KEYS),
            6,
            len(PAIR_DEFINITIONS),
            histogram_config.cal_bins,
            histogram_config.cal_bins,
        ),
        dtype=np.uint32,
    )


def empty_pcal_electron_kinematic_counts(
    histogram_config: HistogramConfig,
) -> np.ndarray:
    # Axes:
    #   sector, requested coordinate pair, x bin, y bin
    return np.zeros(
        (
            6,
            len(KINEMATIC_PAIR_INDICES),
            histogram_config.cal_bins,
            histogram_config.cal_bins,
        ),
        dtype=np.uint32,
    )


def empty_pcal_electron_kinematic_sums(
    histogram_config: HistogramConfig,
) -> np.ndarray:
    # Axes:
    #   kinematic quantity, sector, requested coordinate pair, x bin, y bin
    return np.zeros(
        (
            len(KINEMATIC_KEYS),
            6,
            len(KINEMATIC_PAIR_INDICES),
            histogram_config.cal_bins,
            histogram_config.cal_bins,
        ),
        dtype=np.float64,
    )


def empty_pcal_electron_kinematic_1d_counts(
    histogram_config: HistogramConfig,
) -> np.ndarray:
    # Axes:
    #   sector, coordinate, histogram bin
    return np.zeros(
        (
            6,
            len(COORD_KEYS),
            histogram_config.cal_bins,
        ),
        dtype=np.uint32,
    )


def empty_pcal_electron_kinematic_1d_sums(
    histogram_config: HistogramConfig,
) -> np.ndarray:
    # Axes:
    #   kinematic quantity, sector, coordinate, histogram bin
    return np.zeros(
        (
            len(KINEMATIC_KEYS),
            6,
            len(COORD_KEYS),
            histogram_config.cal_bins,
        ),
        dtype=np.float64,
    )


def accumulate_sample(
    file_path: Path,
    tree_name: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    chunk_size: int,
    max_events: int | None,
    do_ft: bool,
    do_calorimeter: bool,
    dataset_key: str,
    calorimeter_strictness: int,
    is_rgc: bool,
) -> SampleResult:
    start_time = time.perf_counter()
    branches = required_branches(do_ft, do_calorimeter)
    validate_root_tree(file_path, tree_name, branches)

    ft_shape = (histogram_config.ft_bins, histogram_config.ft_bins)
    ft_before = np.zeros(ft_shape, dtype=np.uint64)
    ft_after = np.zeros(ft_shape, dtype=np.uint64)
    cal_counts_before = empty_cal_counts(histogram_config)
    cal_counts_after = empty_cal_counts(histogram_config)
    cal_pair_counts_before = empty_cal_pair_counts(histogram_config)
    cal_pair_counts_after = empty_cal_pair_counts(histogram_config)
    pcal_electron_kin_counts_before = (
        empty_pcal_electron_kinematic_counts(histogram_config)
    )
    pcal_electron_kin_counts_after = (
        empty_pcal_electron_kinematic_counts(histogram_config)
    )
    pcal_electron_kin_sums_before = (
        empty_pcal_electron_kinematic_sums(histogram_config)
    )
    pcal_electron_kin_sums_after = (
        empty_pcal_electron_kinematic_sums(histogram_config)
    )
    pcal_electron_kin_1d_counts_before = (
        empty_pcal_electron_kinematic_1d_counts(histogram_config)
    )
    pcal_electron_kin_1d_counts_after = (
        empty_pcal_electron_kinematic_1d_counts(histogram_config)
    )
    pcal_electron_kin_1d_sums_before = (
        empty_pcal_electron_kinematic_1d_sums(histogram_config)
    )
    pcal_electron_kin_1d_sums_after = (
        empty_pcal_electron_kinematic_1d_sums(histogram_config)
    )

    ft_edges = histogram_config.ft_edges
    cal_edges = histogram_config.cal_edges

    rows_read = 0
    photon_rows = 0
    valid_ft_photons = 0
    fiducial_ft_photons = 0

    iterator = uproot.iterate(
        f"{file_path}:{tree_name}",
        expressions=list(branches),
        step_size=chunk_size,
        entry_stop=max_events,
        library="np",
    )

    for arrays in iterator:
        pid = np.asarray(arrays["particle_pid"], dtype=np.int32)
        rows_read += int(pid.size)

        if do_ft:
            photon_mask = pid == PHOTON_PID
            photon_rows += int(np.count_nonzero(photon_mask))

            ft_x = np.asarray(arrays["ft_x"], dtype=np.float64)
            ft_y = np.asarray(arrays["ft_y"], dtype=np.float64)

            valid_ft_mask = (
                photon_mask
                & np.isfinite(ft_x)
                & np.isfinite(ft_y)
                & (ft_x != INVALID_SENTINEL)
                & (ft_y != INVALID_SENTINEL)
            )

            if np.any(valid_ft_mask):
                valid_x = ft_x[valid_ft_mask]
                valid_y = ft_y[valid_ft_mask]
                valid_ft_photons += int(valid_x.size)

                accepted = ft_fiducial_mask(
                    valid_x,
                    valid_y,
                    fiducial_config,
                )
                fiducial_ft_photons += int(np.count_nonzero(accepted))

                chunk_before, _, _ = np.histogram2d(
                    valid_x,
                    valid_y,
                    bins=(ft_edges, ft_edges),
                )
                ft_before += chunk_before.astype(np.uint64, copy=False)

                if np.any(accepted):
                    chunk_after, _, _ = np.histogram2d(
                        valid_x[accepted],
                        valid_y[accepted],
                        bins=(ft_edges, ft_edges),
                    )
                    ft_after += chunk_after.astype(np.uint64, copy=False)
                #endif
            #endif
        #endif

        if do_calorimeter:
            sectors = np.asarray(arrays["cal_sector"], dtype=np.int16)
            electron_momentum = np.asarray(arrays["p"], dtype=np.float64)
            electron_theta = np.asarray(arrays["theta"], dtype=np.float64)
            valid_sector = (sectors >= 1) & (sectors <= 6)

            # Read each calorimeter coordinate array exactly once per chunk.
            coordinate_arrays = {
                (layer_key, coordinate_key): np.asarray(
                    arrays[f"cal_{coordinate_key}_{layer_number}"],
                    dtype=np.float64,
                )
                for layer_key, layer_number in zip(
                    LAYER_KEYS,
                    LAYER_NUMBERS,
                )
                for coordinate_key in COORD_KEYS
            }

            for particle_index, particle_pid in enumerate(PARTICLE_PIDS):
                particle_mask = valid_sector & (pid == particle_pid)
                if not np.any(particle_mask):
                    continue
                #endif

                after_cut_mask = calorimeter_fiducial_mask(
                    particle_mask=particle_mask,
                    sectors=sectors,
                    coordinate_arrays=coordinate_arrays,
                    dataset_key=dataset_key,
                    strictness=calorimeter_strictness,
                    is_rgc=is_rgc,
                )

                for layer_index, layer_key in enumerate(LAYER_KEYS):
                    layer_arrays = {
                        coordinate_key: coordinate_arrays[
                            (layer_key, coordinate_key)
                        ]
                        for coordinate_key in COORD_KEYS
                    }

                    for sector in range(1, 7):
                        sector_mask = particle_mask & (sectors == sector)
                        sector_after_mask = (
                            after_cut_mask & (sectors == sector)
                        )

                        if not np.any(sector_mask):
                            continue
                        #endif

                        for coordinate_index, coordinate_key in enumerate(COORD_KEYS):
                            values = layer_arrays[coordinate_key]
                            valid_coordinate = (
                                np.isfinite(values)
                                & (values != INVALID_SENTINEL)
                                & (values >= histogram_config.cal_min_cm)
                                & (values <= histogram_config.cal_max_cm)
                            )

                            before_mask = sector_mask & valid_coordinate
                            if np.any(before_mask):
                                counts, _ = np.histogram(
                                    values[before_mask],
                                    bins=cal_edges,
                                )
                                cal_counts_before[
                                    particle_index,
                                    layer_index,
                                    sector - 1,
                                    coordinate_index,
                                ] += counts.astype(np.uint64, copy=False)
                            #endif

                            after_mask = sector_after_mask & valid_coordinate
                            if np.any(after_mask):
                                counts, _ = np.histogram(
                                    values[after_mask],
                                    bins=cal_edges,
                                )
                                cal_counts_after[
                                    particle_index,
                                    layer_index,
                                    sector - 1,
                                    coordinate_index,
                                ] += counts.astype(np.uint64, copy=False)
                            #endif
                        #endfor

                        for pair_index, (_, x_key, y_key) in enumerate(PAIR_DEFINITIONS):
                            x_values = layer_arrays[x_key]
                            y_values = layer_arrays[y_key]
                            valid_pair = (
                                np.isfinite(x_values)
                                & np.isfinite(y_values)
                                & (x_values != INVALID_SENTINEL)
                                & (y_values != INVALID_SENTINEL)
                                & (x_values >= histogram_config.cal_min_cm)
                                & (x_values <= histogram_config.cal_max_cm)
                                & (y_values >= histogram_config.cal_min_cm)
                                & (y_values <= histogram_config.cal_max_cm)
                            )

                            before_pair_mask = sector_mask & valid_pair
                            if np.any(before_pair_mask):
                                pair_counts, _, _ = np.histogram2d(
                                    x_values[before_pair_mask],
                                    y_values[before_pair_mask],
                                    bins=(cal_edges, cal_edges),
                                )
                                cal_pair_counts_before[
                                    particle_index,
                                    layer_index,
                                    sector - 1,
                                    pair_index,
                                ] += pair_counts.astype(np.uint32, copy=False)
                            #endif

                            after_pair_mask = sector_after_mask & valid_pair
                            if np.any(after_pair_mask):
                                pair_counts, _, _ = np.histogram2d(
                                    x_values[after_pair_mask],
                                    y_values[after_pair_mask],
                                    bins=(cal_edges, cal_edges),
                                )
                                cal_pair_counts_after[
                                    particle_index,
                                    layer_index,
                                    sector - 1,
                                    pair_index,
                                ] += pair_counts.astype(np.uint32, copy=False)
                            #endif
                        #endfor

                        # Mean electron momentum and theta maps are needed only
                        # for PCal and the lu-lv / lu-lw coordinate pairs.
                        if particle_pid == ELECTRON_PID and layer_key == "pcal":
                            valid_momentum = (
                                np.isfinite(electron_momentum)
                                & (electron_momentum > 0.0)
                            )
                            valid_theta = (
                                np.isfinite(electron_theta)
                                & (electron_theta >= 0.0)
                            )

                            for coordinate_index, coordinate_key in enumerate(
                                COORD_KEYS
                            ):
                                coordinate_values = layer_arrays[
                                    coordinate_key
                                ]
                                valid_coordinate_1d = (
                                    np.isfinite(coordinate_values)
                                    & (
                                        coordinate_values
                                        != INVALID_SENTINEL
                                    )
                                    & (
                                        coordinate_values
                                        >= histogram_config.cal_min_cm
                                    )
                                    & (
                                        coordinate_values
                                        <= histogram_config.cal_max_cm
                                    )
                                )
                                valid_kinematics_1d = (
                                    valid_coordinate_1d
                                    & valid_momentum
                                    & valid_theta
                                )

                                before_1d_mask = (
                                    sector_mask & valid_kinematics_1d
                                )
                                if np.any(before_1d_mask):
                                    counts_1d, _ = np.histogram(
                                        coordinate_values[before_1d_mask],
                                        bins=cal_edges,
                                    )
                                    momentum_sum_1d, _ = np.histogram(
                                        coordinate_values[before_1d_mask],
                                        bins=cal_edges,
                                        weights=electron_momentum[
                                            before_1d_mask
                                        ],
                                    )
                                    theta_sum_1d, _ = np.histogram(
                                        coordinate_values[before_1d_mask],
                                        bins=cal_edges,
                                        weights=electron_theta[
                                            before_1d_mask
                                        ],
                                    )

                                    pcal_electron_kin_1d_counts_before[
                                        sector - 1,
                                        coordinate_index,
                                    ] += counts_1d.astype(
                                        np.uint32,
                                        copy=False,
                                    )
                                    pcal_electron_kin_1d_sums_before[
                                        0,
                                        sector - 1,
                                        coordinate_index,
                                    ] += momentum_sum_1d
                                    pcal_electron_kin_1d_sums_before[
                                        1,
                                        sector - 1,
                                        coordinate_index,
                                    ] += theta_sum_1d
                                #endif

                                after_1d_mask = (
                                    sector_after_mask
                                    & valid_kinematics_1d
                                )
                                if np.any(after_1d_mask):
                                    counts_1d, _ = np.histogram(
                                        coordinate_values[after_1d_mask],
                                        bins=cal_edges,
                                    )
                                    momentum_sum_1d, _ = np.histogram(
                                        coordinate_values[after_1d_mask],
                                        bins=cal_edges,
                                        weights=electron_momentum[
                                            after_1d_mask
                                        ],
                                    )
                                    theta_sum_1d, _ = np.histogram(
                                        coordinate_values[after_1d_mask],
                                        bins=cal_edges,
                                        weights=electron_theta[
                                            after_1d_mask
                                        ],
                                    )

                                    pcal_electron_kin_1d_counts_after[
                                        sector - 1,
                                        coordinate_index,
                                    ] += counts_1d.astype(
                                        np.uint32,
                                        copy=False,
                                    )
                                    pcal_electron_kin_1d_sums_after[
                                        0,
                                        sector - 1,
                                        coordinate_index,
                                    ] += momentum_sum_1d
                                    pcal_electron_kin_1d_sums_after[
                                        1,
                                        sector - 1,
                                        coordinate_index,
                                    ] += theta_sum_1d
                                #endif
                            #endfor

                            for kin_pair_slot, pair_index in enumerate(
                                KINEMATIC_PAIR_INDICES
                            ):
                                _, x_key, y_key = PAIR_DEFINITIONS[pair_index]
                                x_values = layer_arrays[x_key]
                                y_values = layer_arrays[y_key]

                                valid_coordinates = (
                                    np.isfinite(x_values)
                                    & np.isfinite(y_values)
                                    & (x_values != INVALID_SENTINEL)
                                    & (y_values != INVALID_SENTINEL)
                                    & (
                                        x_values
                                        >= histogram_config.cal_min_cm
                                    )
                                    & (
                                        x_values
                                        <= histogram_config.cal_max_cm
                                    )
                                    & (
                                        y_values
                                        >= histogram_config.cal_min_cm
                                    )
                                    & (
                                        y_values
                                        <= histogram_config.cal_max_cm
                                    )
                                )
                                valid_kinematics = (
                                    valid_coordinates
                                    & valid_momentum
                                    & valid_theta
                                )

                                before_kin_mask = (
                                    sector_mask & valid_kinematics
                                )
                                if np.any(before_kin_mask):
                                    kin_counts, _, _ = np.histogram2d(
                                        x_values[before_kin_mask],
                                        y_values[before_kin_mask],
                                        bins=(cal_edges, cal_edges),
                                    )
                                    momentum_sum, _, _ = np.histogram2d(
                                        x_values[before_kin_mask],
                                        y_values[before_kin_mask],
                                        bins=(cal_edges, cal_edges),
                                        weights=electron_momentum[
                                            before_kin_mask
                                        ],
                                    )
                                    theta_sum, _, _ = np.histogram2d(
                                        x_values[before_kin_mask],
                                        y_values[before_kin_mask],
                                        bins=(cal_edges, cal_edges),
                                        weights=electron_theta[
                                            before_kin_mask
                                        ],
                                    )

                                    pcal_electron_kin_counts_before[
                                        sector - 1,
                                        kin_pair_slot,
                                    ] += kin_counts.astype(
                                        np.uint32,
                                        copy=False,
                                    )
                                    pcal_electron_kin_sums_before[
                                        0,
                                        sector - 1,
                                        kin_pair_slot,
                                    ] += momentum_sum
                                    pcal_electron_kin_sums_before[
                                        1,
                                        sector - 1,
                                        kin_pair_slot,
                                    ] += theta_sum
                                #endif

                                after_kin_mask = (
                                    sector_after_mask & valid_kinematics
                                )
                                if np.any(after_kin_mask):
                                    kin_counts, _, _ = np.histogram2d(
                                        x_values[after_kin_mask],
                                        y_values[after_kin_mask],
                                        bins=(cal_edges, cal_edges),
                                    )
                                    momentum_sum, _, _ = np.histogram2d(
                                        x_values[after_kin_mask],
                                        y_values[after_kin_mask],
                                        bins=(cal_edges, cal_edges),
                                        weights=electron_momentum[
                                            after_kin_mask
                                        ],
                                    )
                                    theta_sum, _, _ = np.histogram2d(
                                        x_values[after_kin_mask],
                                        y_values[after_kin_mask],
                                        bins=(cal_edges, cal_edges),
                                        weights=electron_theta[
                                            after_kin_mask
                                        ],
                                    )

                                    pcal_electron_kin_counts_after[
                                        sector - 1,
                                        kin_pair_slot,
                                    ] += kin_counts.astype(
                                        np.uint32,
                                        copy=False,
                                    )
                                    pcal_electron_kin_sums_after[
                                        0,
                                        sector - 1,
                                        kin_pair_slot,
                                    ] += momentum_sum
                                    pcal_electron_kin_sums_after[
                                        1,
                                        sector - 1,
                                        kin_pair_slot,
                                    ] += theta_sum
                                #endif
                            #endfor
                        #endif
                    #endfor
                #endfor
            #endfor
        #endif
    #endfor

    return SampleResult(
        ft_before=ft_before,
        ft_after=ft_after,
        cal_counts_before=cal_counts_before,
        cal_counts_after=cal_counts_after,
        cal_pair_counts_before=cal_pair_counts_before,
        cal_pair_counts_after=cal_pair_counts_after,
        pcal_electron_kin_counts_before=(
            pcal_electron_kin_counts_before
        ),
        pcal_electron_kin_counts_after=(
            pcal_electron_kin_counts_after
        ),
        pcal_electron_kin_sums_before=pcal_electron_kin_sums_before,
        pcal_electron_kin_sums_after=pcal_electron_kin_sums_after,
        pcal_electron_kin_1d_counts_before=(
            pcal_electron_kin_1d_counts_before
        ),
        pcal_electron_kin_1d_counts_after=(
            pcal_electron_kin_1d_counts_after
        ),
        pcal_electron_kin_1d_sums_before=(
            pcal_electron_kin_1d_sums_before
        ),
        pcal_electron_kin_1d_sums_after=(
            pcal_electron_kin_1d_sums_after
        ),
        rows_read=rows_read,
        photon_rows=photon_rows,
        valid_ft_photons=valid_ft_photons,
        fiducial_ft_photons=fiducial_ft_photons,
        elapsed_seconds=time.perf_counter() - start_time,
    )


def process_period(
    dataset: Dataset,
    input_dir: Path,
    tree_name: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    chunk_size: int,
    max_events: int | None,
    do_ft: bool,
    do_calorimeter: bool,
    calorimeter_strictness: int,
    is_rgc: bool,
) -> PeriodResult:
    start_time = time.perf_counter()

    data = accumulate_sample(
        input_dir / dataset.data_file,
        tree_name,
        histogram_config,
        fiducial_config,
        chunk_size,
        max_events,
        do_ft,
        do_calorimeter,
        dataset.key,
        calorimeter_strictness,
        is_rgc,
    )

    mc = None
    if dataset.mc_file is not None:
        mc = accumulate_sample(
            input_dir / dataset.mc_file,
            tree_name,
            histogram_config,
            fiducial_config,
            chunk_size,
            max_events,
            do_ft,
            do_calorimeter,
            dataset.key,
            calorimeter_strictness,
            is_rgc,
        )
    #endif

    return PeriodResult(
        dataset=dataset,
        data=data,
        mc=mc,
        worker_pid=os.getpid(),
        elapsed_seconds=time.perf_counter() - start_time,
    )


# =============================================================================
# FT plotting
# =============================================================================


def draw_ft_boundaries(
    axis: plt.Axes,
    config: FTFiducialConfig,
) -> None:
    for radius_cm in (config.inner_radius_cm, config.outer_radius_cm):
        axis.add_patch(
            Circle(
                (0.0, 0.0),
                radius_cm,
                fill=False,
                linewidth=1.5,
                linestyle="--",
                edgecolor="black",
            )
        )
    #endfor

    for center_x_cm, center_y_cm, radius_cm in config.holes:
        axis.add_patch(
            Circle(
                (center_x_cm, center_y_cm),
                radius_cm,
                fill=False,
                linewidth=1.5,
                edgecolor="black",
            )
        )
    #endfor


def configure_ft_axis(
    axis: plt.Axes,
    config: HistogramConfig,
) -> None:
    axis.set_xlabel(r"FT $x$ (cm)")
    axis.set_ylabel(r"FT $y$ (cm)")
    axis.set_xlim(config.ft_xy_min_cm, config.ft_xy_max_cm)
    axis.set_ylim(config.ft_xy_min_cm, config.ft_xy_max_cm)
    axis.set_aspect("equal", adjustable="box")


def shared_log_norm(
    first: np.ndarray,
    second: np.ndarray,
) -> LogNorm | None:
    maximum = max(float(np.max(first)), float(np.max(second)))
    if maximum <= 0.0:
        return None
    #endif
    return LogNorm(vmin=1.0, vmax=max(1.0, maximum))


def draw_ft_panel(
    axis: plt.Axes,
    histogram: np.ndarray,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    norm: LogNorm | None,
):
    extent = (
        histogram_config.ft_xy_min_cm,
        histogram_config.ft_xy_max_cm,
        histogram_config.ft_xy_min_cm,
        histogram_config.ft_xy_max_cm,
    )
    image = axis.imshow(
        histogram.T,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        cmap="viridis",
        norm=norm,
        rasterized=True,
    )
    axis.set_title(title)
    configure_ft_axis(axis, histogram_config)
    draw_ft_boundaries(axis, fiducial_config)
    return image


def ft_efficiency_text(result: SampleResult) -> str:
    if result.valid_ft_photons <= 0:
        return "0 / 0 valid FT photons"
    #endif

    efficiency = (
        100.0
        * result.fiducial_ft_photons
        / result.valid_ft_photons
    )
    return (
        f"{result.fiducial_ft_photons:,} / "
        f"{result.valid_ft_photons:,} valid FT photons "
        f"({efficiency:.2f}%)"
    )


def draw_missing_mc_ft_panel(
    axis: plt.Axes,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
) -> None:
    axis.set_title(title)
    configure_ft_axis(axis, histogram_config)
    draw_ft_boundaries(axis, fiducial_config)
    axis.text(
        0.5,
        0.5,
        "MC not available\nfor RGC mode",
        transform=axis.transAxes,
        ha="center",
        va="center",
        fontsize=14,
    )


def save_ft_summary(
    result: PeriodResult,
    output_base: Path,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    dpi: int,
) -> Path:
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.5, 11.5),
        constrained_layout=True,
    )

    data_norm = shared_log_norm(
        result.data.ft_before,
        result.data.ft_after,
    )
    data_image = draw_ft_panel(
        axes[0, 0],
        result.data.ft_before,
        "Data before FT fiducial cut",
        histogram_config,
        fiducial_config,
        data_norm,
    )
    draw_ft_panel(
        axes[0, 1],
        result.data.ft_after,
        "Data after FT fiducial cut",
        histogram_config,
        fiducial_config,
        data_norm,
    )
    colorbar = figure.colorbar(
        data_image,
        ax=axes[0, :],
        pad=0.015,
        fraction=0.035,
    )
    colorbar.set_label("Photon entries per bin")

    if result.mc is not None:
        mc_norm = shared_log_norm(
            result.mc.ft_before,
            result.mc.ft_after,
        )
        mc_image = draw_ft_panel(
            axes[1, 0],
            result.mc.ft_before,
            "MC before FT fiducial cut",
            histogram_config,
            fiducial_config,
            mc_norm,
        )
        draw_ft_panel(
            axes[1, 1],
            result.mc.ft_after,
            "MC after FT fiducial cut",
            histogram_config,
            fiducial_config,
            mc_norm,
        )
        colorbar = figure.colorbar(
            mc_image,
            ax=axes[1, :],
            pad=0.015,
            fraction=0.035,
        )
        colorbar.set_label("Photon entries per bin")
    else:
        draw_missing_mc_ft_panel(
            axes[1, 0],
            "MC before FT fiducial cut",
            histogram_config,
            fiducial_config,
        )
        draw_missing_mc_ft_panel(
            axes[1, 1],
            "MC after FT fiducial cut",
            histogram_config,
            fiducial_config,
        )
    #endif

    figure.suptitle(
        f"{result.dataset.label}: Forward Tagger photon occupancy",
        fontsize=17,
    )

    data_caption = f"Data: {ft_efficiency_text(result.data)}"
    mc_caption = (
        f"MC: {ft_efficiency_text(result.mc)}"
        if result.mc is not None
        else "MC: not available"
    )
    figure.text(
        0.5,
        0.005,
        f"{data_caption}    |    {mc_caption}",
        ha="center",
        va="bottom",
        fontsize=10,
    )

    output_dir = output_base / "forward_tagger"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{result.dataset.key}_ft_summary.png"
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return output_path


# =============================================================================
# Calorimeter strip-matching plotting
# =============================================================================


def normalized_histogram(
    counts: np.ndarray,
    bin_width_cm: float,
) -> np.ndarray:
    total = float(np.sum(counts))
    if total <= 0.0:
        return np.zeros_like(counts, dtype=np.float64)
    #endif
    return counts.astype(np.float64) / (total * bin_width_cm)


def draw_exclusions(
    axis: plt.Axes,
    intervals: list[tuple[float, float, str]],
) -> None:
    for lower, upper, _ in intervals:
        axis.axvspan(
            lower,
            upper,
            alpha=0.18,
            hatch="//",
            edgecolor="black",
            linewidth=0.8,
        )
    #endfor


def positive_minimum(*arrays: np.ndarray) -> float | None:
    minima: list[float] = []
    for array in arrays:
        positive = np.asarray(array)[np.asarray(array) > 0.0]
        if positive.size > 0:
            minima.append(float(np.min(positive)))
        #endif
    #endfor
    return min(minima) if minima else None


def save_calorimeter_matching_summary(
    result: PeriodResult,
    output_base: Path,
    histogram_config: HistogramConfig,
    layer_index: int,
    particle_index: int,
    strictness: int,
    is_rgc: bool,
    cut_stage: str,
    dpi: int,
) -> Path:
    if cut_stage not in ("before_cuts", "after_cuts"):
        raise ValueError(f"Unknown calorimeter cut stage: {cut_stage}")
    #endif

    layer_key = LAYER_KEYS[layer_index]
    layer_label = LAYER_LABELS[layer_key]
    particle_key = PARTICLE_KEYS[particle_index]
    particle_label = particle_key.capitalize()

    edges = histogram_config.cal_edges
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width_cm = histogram_config.cal_bin_width_cm

    figure, axes = plt.subplots(
        6,
        3,
        figsize=(17.0, 20.0),
        sharex=True,
        constrained_layout=True,
    )

    legend_handles = None
    any_exclusion = False

    for sector_index in range(6):
        sector = sector_index + 1

        for coordinate_index, coordinate_key in enumerate(COORD_KEYS):
            axis = axes[sector_index, coordinate_index]

            data_counts_array = (
                result.data.cal_counts_before
                if cut_stage == "before_cuts"
                else result.data.cal_counts_after
            )
            data_counts = data_counts_array[
                particle_index,
                layer_index,
                sector_index,
                coordinate_index,
            ]
            data_density = normalized_histogram(data_counts, bin_width_cm)

            data_line = axis.step(
                centers,
                data_density,
                where="mid",
                linewidth=1.35,
                label="Data",
            )[0]

            mc_line = None
            mc_density = np.zeros_like(data_density)
            if result.mc is not None:
                mc_counts_array = (
                    result.mc.cal_counts_before
                    if cut_stage == "before_cuts"
                    else result.mc.cal_counts_after
                )
                mc_counts = mc_counts_array[
                    particle_index,
                    layer_index,
                    sector_index,
                    coordinate_index,
                ]
                mc_density = normalized_histogram(mc_counts, bin_width_cm)
                mc_line = axis.step(
                    centers,
                    mc_density,
                    where="mid",
                    linewidth=1.25,
                    label="MC",
                )[0]
            #endif

            intervals = exclusion_intervals(
                dataset_key=result.dataset.key,
                layer_key=layer_key,
                sector=sector,
                coordinate_key=coordinate_key,
                strictness=strictness,
                is_rgc=is_rgc,
            )
            if intervals:
                any_exclusion = True
                draw_exclusions(axis, intervals)
            #endif

            if sector_index == 0:
                axis.set_title(COORD_LABELS[coordinate_key], fontsize=13)
            #endif

            if coordinate_index == 0:
                axis.set_ylabel(
                    f"Sector {sector}\nNormalized entries / cm",
                    fontsize=10,
                )
            #endif

            axis.set_xlim(
                histogram_config.cal_min_cm,
                histogram_config.cal_max_cm,
            )
            axis.set_yscale("log")
            minimum = positive_minimum(data_density, mc_density)
            maximum = max(float(np.max(data_density)), float(np.max(mc_density)))
            if minimum is not None and maximum > 0.0:
                axis.set_ylim(
                    max(minimum * 0.5, 1.0e-8),
                    maximum * 2.0,
                )
            #endif
            axis.grid(alpha=0.18, which="both")

            data_entries = int(np.sum(data_counts))
            if result.mc is not None:
                mc_entries = int(
                    np.sum(
                        mc_counts_array[
                            particle_index,
                            layer_index,
                            sector_index,
                            coordinate_index,
                        ]
                    )
                )
                entry_text = f"D {data_entries:,}\nMC {mc_entries:,}"
            else:
                entry_text = f"D {data_entries:,}"
            #endif

            axis.text(
                0.985,
                0.95,
                entry_text,
                transform=axis.transAxes,
                ha="right",
                va="top",
                fontsize=7.5,
            )

            if legend_handles is None:
                legend_handles = (
                    [data_line, mc_line]
                    if mc_line is not None
                    else [data_line]
                )
            #endif
        #endfor
    #endfor

    for axis in axes[-1, :]:
        axis.set_xlabel("Strip coordinate (cm)")
    #endfor

    cut_description = (
        "fiducial exclusions are displayed but not applied"
        if cut_stage == "before_cuts"
        else "complete calorimeter fiducial cuts are applied"
    )
    subtitle = (
        "4.5 cm bins; logarithmic scale; independent unit-area "
        f"normalization; {cut_description}"
    )
    if not any_exclusion:
        subtitle = (
            "4.5 cm bins; logarithmic scale; independent unit-area "
            f"normalization; {cut_description}"
        )
    #endif

    figure.suptitle(
        f"{result.dataset.label}: {layer_label} {particle_label} "
        f"data/MC strip matching\n{subtitle}",
        fontsize=17,
    )

    if legend_handles is not None:
        legend_handles = [
            handle for handle in legend_handles if handle is not None
        ]
        legend_labels = [handle.get_label() for handle in legend_handles]
        figure.legend(
            legend_handles,
            legend_labels,
            loc="upper right",
            bbox_to_anchor=(0.985, 0.985),
            frameon=True,
        )
    #endif

    figure.text(
        0.012,
        0.005,
        (
            f"Overlay strictness = {strictness}. "
            + (
                "Raw valid calorimeter-coordinate entries are plotted; "
                "fiducial exclusions are not imposed."
                if cut_stage == "before_cuts"
                else "The complete particle-level calorimeter fiducial "
                "selection is imposed."
            )
        ),
        ha="left",
        va="bottom",
        fontsize=9,
    )

    output_dir = output_base / "calorimeter" / cut_stage / layer_key
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (
        f"{result.dataset.key}_{layer_key}_{particle_key}_strip_matching.png"
    )
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return output_path


def row_log_norm(
    histograms: list[np.ndarray],
) -> LogNorm | None:
    """
    Build one logarithmic normalization for one complete data or MC row.

    Zero-count bins are ignored. The resulting normalization is shared across
    all six sectors in that row so sector-to-sector occupancy differences
    remain visible, while data and MC use independent scales.
    """

    positive_minima: list[float] = []
    maxima: list[float] = []

    for histogram in histograms:
        positive = histogram[histogram > 0]
        if positive.size == 0:
            continue
        #endif

        positive_minima.append(float(np.min(positive)))
        maxima.append(float(np.max(positive)))
    #endfor

    if not maxima:
        return None
    #endif

    minimum = min(positive_minima)
    maximum = max(maxima)

    return LogNorm(
        vmin=max(1.0, minimum),
        vmax=max(maximum, minimum),
    )


def draw_2d_exclusions(
    axis: plt.Axes,
    x_intervals: list[tuple[float, float, str]],
    y_intervals: list[tuple[float, float, str]],
) -> None:
    for lower, upper, _ in x_intervals:
        axis.axvspan(
            lower,
            upper,
            alpha=0.12,
            hatch="//",
            edgecolor="white",
            linewidth=0.6,
        )
    #endfor

    for lower, upper, _ in y_intervals:
        axis.axhspan(
            lower,
            upper,
            alpha=0.12,
            hatch="\\\\",
            edgecolor="white",
            linewidth=0.6,
        )
    #endfor


def save_calorimeter_2d_matching_summary(
    result: PeriodResult,
    output_base: Path,
    histogram_config: HistogramConfig,
    layer_index: int,
    particle_index: int,
    pair_index: int,
    strictness: int,
    is_rgc: bool,
    cut_stage: str,
    dpi: int,
) -> Path:
    if cut_stage not in ("before_cuts", "after_cuts"):
        raise ValueError(f"Unknown calorimeter cut stage: {cut_stage}")
    #endif

    layer_key = LAYER_KEYS[layer_index]
    layer_label = LAYER_LABELS[layer_key]
    particle_key = PARTICLE_KEYS[particle_index]
    particle_label = particle_key.capitalize()
    pair_key, x_key, y_key = PAIR_DEFINITIONS[pair_index]

    has_mc = result.mc is not None
    number_of_rows = 2 if has_mc else 1

    # Manual spacing is used here rather than constrained_layout because the
    # latter leaves excessive vertical and horizontal gaps for a 2x6 canvas.
    figure, axes = plt.subplots(
        number_of_rows,
        6,
        figsize=(19.5, 6.9 if has_mc else 3.9),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    figure.subplots_adjust(
        left=0.050,
        right=0.925,
        bottom=0.105,
        top=0.835,
        wspace=0.075,
        hspace=0.090,
    )

    data_pair_counts = (
        result.data.cal_pair_counts_before
        if cut_stage == "before_cuts"
        else result.data.cal_pair_counts_after
    )
    data_histograms = [
        data_pair_counts[
            particle_index,
            layer_index,
            sector_index,
            pair_index,
        ]
        for sector_index in range(6)
    ]
    data_norm = row_log_norm(data_histograms)

    mc_histograms: list[np.ndarray] = []
    mc_norm = None
    if has_mc:
        assert result.mc is not None
        mc_pair_counts = (
            result.mc.cal_pair_counts_before
            if cut_stage == "before_cuts"
            else result.mc.cal_pair_counts_after
        )
        mc_histograms = [
            mc_pair_counts[
                particle_index,
                layer_index,
                sector_index,
                pair_index,
            ]
            for sector_index in range(6)
        ]
        mc_norm = row_log_norm(mc_histograms)
    #endif

    extent = (
        histogram_config.cal_min_cm,
        histogram_config.cal_max_cm,
        histogram_config.cal_min_cm,
        histogram_config.cal_max_cm,
    )

    first_data_image = None
    first_mc_image = None

    for sector_index in range(6):
        sector = sector_index + 1

        x_intervals = exclusion_intervals(
            result.dataset.key,
            layer_key,
            sector,
            x_key,
            strictness,
            is_rgc,
        )
        y_intervals = exclusion_intervals(
            result.dataset.key,
            layer_key,
            sector,
            y_key,
            strictness,
            is_rgc,
        )

        data_axis = axes[0, sector_index]
        data_image = data_axis.imshow(
            data_histograms[sector_index].T,
            origin="lower",
            extent=extent,
            interpolation="nearest",
            aspect="equal",
            cmap="viridis",
            norm=data_norm,
            rasterized=True,
        )
        if first_data_image is None:
            first_data_image = data_image
        #endif

        data_axis.set_title(f"Sector {sector}", fontsize=11, pad=4)
        if not has_mc:
            data_axis.set_xlabel(COORD_LABELS[x_key], fontsize=9)
        else:
            data_axis.tick_params(labelbottom=False)
        #endif
        if sector_index == 0:
            data_axis.set_ylabel(
                f"Data\n{COORD_LABELS[y_key]}",
                fontsize=9,
            )
        #endif
        data_axis.tick_params(labelsize=8)
        draw_2d_exclusions(data_axis, x_intervals, y_intervals)

        if has_mc:
            assert result.mc is not None
            mc_axis = axes[1, sector_index]
            mc_image = mc_axis.imshow(
                mc_histograms[sector_index].T,
                origin="lower",
                extent=extent,
                interpolation="nearest",
                aspect="equal",
                cmap="viridis",
                norm=mc_norm,
                rasterized=True,
            )
            if first_mc_image is None:
                first_mc_image = mc_image
            #endif

            mc_axis.set_xlabel(COORD_LABELS[x_key], fontsize=9)
            if sector_index == 0:
                mc_axis.set_ylabel(
                    f"MC\n{COORD_LABELS[y_key]}",
                    fontsize=9,
                )
            #endif
            mc_axis.tick_params(labelsize=8)
            draw_2d_exclusions(mc_axis, x_intervals, y_intervals)
        #endif
    #endfor

    # Separate colorbars preserve the very different absolute statistics in
    # data and MC. Each bar applies to all six sectors in its corresponding row.
    if first_data_image is not None:
        data_colorbar_axis = figure.add_axes(
            [0.942, 0.535 if has_mc else 0.175, 0.012, 0.300 if has_mc else 0.625]
        )
        data_colorbar = figure.colorbar(
            first_data_image,
            cax=data_colorbar_axis,
        )
        data_colorbar.set_label(
            "Data entries per 4.5 cm × 4.5 cm bin",
            fontsize=9,
        )
        data_colorbar.ax.tick_params(labelsize=8)
    #endif

    if has_mc and first_mc_image is not None:
        mc_colorbar_axis = figure.add_axes(
            [0.942, 0.105, 0.012, 0.300]
        )
        mc_colorbar = figure.colorbar(
            first_mc_image,
            cax=mc_colorbar_axis,
        )
        mc_colorbar.set_label(
            "MC entries per 4.5 cm × 4.5 cm bin",
            fontsize=9,
        )
        mc_colorbar.ax.tick_params(labelsize=8)
    #endif

    figure.suptitle(
        f"{result.dataset.label}: {layer_label} {particle_label} "
        f"{x_key}-{y_key} occupancy\n"
        "Raw entries; logarithmic color scales are independent for data and MC; "
        + (
            "shaded exclusions are displayed but not applied"
            if cut_stage == "before_cuts"
            else "complete calorimeter fiducial cuts are applied"
        ),
        fontsize=14,
        y=0.965,
    )

    output_dir = output_base / "calorimeter" / cut_stage / layer_key
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (
        f"{result.dataset.key}_{layer_key}_{particle_key}_"
        f"{pair_key}_2d_matching.png"
    )
    figure.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0.08,
    )
    plt.close(figure)
    return output_path


def mean_map(
    sums: np.ndarray,
    counts: np.ndarray,
    minimum_count: int,
) -> np.ma.MaskedArray:
    valid = counts >= minimum_count
    means = np.zeros_like(sums, dtype=np.float64)
    np.divide(
        sums,
        counts,
        out=means,
        where=valid,
    )
    return np.ma.array(means, mask=~valid)


def finite_mean_limits(
    maps: list[np.ma.MaskedArray],
) -> tuple[float, float] | None:
    values: list[np.ndarray] = []

    for mean_values in maps:
        compressed = mean_values.compressed()
        compressed = compressed[np.isfinite(compressed)]
        if compressed.size > 0:
            values.append(compressed)
        #endif
    #endfor

    if not values:
        return None
    #endif

    combined = np.concatenate(values)
    minimum = float(np.min(combined))
    maximum = float(np.max(combined))

    if np.isclose(minimum, maximum):
        padding = max(abs(minimum) * 0.01, 1.0e-6)
        return minimum - padding, maximum + padding
    #endif

    return minimum, maximum


def save_pcal_electron_kinematic_summary(
    result: PeriodResult,
    output_base: Path,
    histogram_config: HistogramConfig,
    pair_slot: int,
    kinematic_index: int,
    strictness: int,
    is_rgc: bool,
    cut_stage: str,
    minimum_bin_count: int,
    dpi: int,
) -> Path:
    if cut_stage not in ("before_cuts", "after_cuts"):
        raise ValueError(f"Unknown calorimeter cut stage: {cut_stage}")
    #endif

    pair_index = KINEMATIC_PAIR_INDICES[pair_slot]
    pair_key, x_key, y_key = PAIR_DEFINITIONS[pair_index]
    kinematic_key = KINEMATIC_KEYS[kinematic_index]
    kinematic_label = KINEMATIC_LABELS[kinematic_key]

    has_mc = result.mc is not None
    number_of_rows = 2 if has_mc else 1

    figure, axes = plt.subplots(
        number_of_rows,
        6,
        figsize=(19.5, 6.9 if has_mc else 3.9),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    figure.subplots_adjust(
        left=0.050,
        right=0.925,
        bottom=0.105,
        top=0.835,
        wspace=0.075,
        hspace=0.090,
    )

    data_counts = (
        result.data.pcal_electron_kin_counts_before
        if cut_stage == "before_cuts"
        else result.data.pcal_electron_kin_counts_after
    )
    data_sums = (
        result.data.pcal_electron_kin_sums_before
        if cut_stage == "before_cuts"
        else result.data.pcal_electron_kin_sums_after
    )

    data_maps = [
        mean_map(
            data_sums[kinematic_index, sector_index, pair_slot],
            data_counts[sector_index, pair_slot],
            minimum_bin_count,
        )
        for sector_index in range(6)
    ]

    mc_maps: list[np.ma.MaskedArray] = []
    if has_mc:
        assert result.mc is not None
        mc_counts = (
            result.mc.pcal_electron_kin_counts_before
            if cut_stage == "before_cuts"
            else result.mc.pcal_electron_kin_counts_after
        )
        mc_sums = (
            result.mc.pcal_electron_kin_sums_before
            if cut_stage == "before_cuts"
            else result.mc.pcal_electron_kin_sums_after
        )
        mc_maps = [
            mean_map(
                mc_sums[kinematic_index, sector_index, pair_slot],
                mc_counts[sector_index, pair_slot],
                minimum_bin_count,
            )
            for sector_index in range(6)
        ]
    #endif

    limits = finite_mean_limits(data_maps + mc_maps)
    if limits is None:
        vmin = None
        vmax = None
    else:
        vmin, vmax = limits
    #endif

    extent = (
        histogram_config.cal_min_cm,
        histogram_config.cal_max_cm,
        histogram_config.cal_min_cm,
        histogram_config.cal_max_cm,
    )

    colormap = plt.get_cmap("viridis").copy()
    colormap.set_bad("white")

    first_image = None

    for sector_index in range(6):
        sector = sector_index + 1

        x_intervals = exclusion_intervals(
            result.dataset.key,
            "pcal",
            sector,
            x_key,
            strictness,
            is_rgc,
        )
        y_intervals = exclusion_intervals(
            result.dataset.key,
            "pcal",
            sector,
            y_key,
            strictness,
            is_rgc,
        )

        data_axis = axes[0, sector_index]
        data_image = data_axis.imshow(
            data_maps[sector_index].T,
            origin="lower",
            extent=extent,
            interpolation="nearest",
            aspect="equal",
            cmap=colormap,
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )
        if first_image is None:
            first_image = data_image
        #endif

        data_axis.set_title(f"Sector {sector}", fontsize=11, pad=4)
        if has_mc:
            data_axis.tick_params(labelbottom=False)
        else:
            data_axis.set_xlabel(COORD_LABELS[x_key], fontsize=9)
        #endif
        if sector_index == 0:
            data_axis.set_ylabel(
                f"Data\n{COORD_LABELS[y_key]}",
                fontsize=9,
            )
        #endif
        data_axis.tick_params(labelsize=8)
        draw_2d_exclusions(data_axis, x_intervals, y_intervals)

        if has_mc:
            mc_axis = axes[1, sector_index]
            mc_axis.imshow(
                mc_maps[sector_index].T,
                origin="lower",
                extent=extent,
                interpolation="nearest",
                aspect="equal",
                cmap=colormap,
                vmin=vmin,
                vmax=vmax,
                rasterized=True,
            )
            mc_axis.set_xlabel(COORD_LABELS[x_key], fontsize=9)
            if sector_index == 0:
                mc_axis.set_ylabel(
                    f"MC\n{COORD_LABELS[y_key]}",
                    fontsize=9,
                )
            #endif
            mc_axis.tick_params(labelsize=8)
            draw_2d_exclusions(mc_axis, x_intervals, y_intervals)
        #endif
    #endfor

    if first_image is not None:
        colorbar_axis = figure.add_axes(
            [0.942, 0.105, 0.012, 0.730]
        )
        colorbar = figure.colorbar(
            first_image,
            cax=colorbar_axis,
        )
        colorbar.set_label(kinematic_label, fontsize=9)
        colorbar.ax.tick_params(labelsize=8)
    #endif

    cut_description = (
        "fiducial exclusions are displayed but not applied"
        if cut_stage == "before_cuts"
        else "complete calorimeter fiducial cuts are applied"
    )

    figure.suptitle(
        f"{result.dataset.label}: PCal electron {kinematic_key} "
        f"versus {x_key}-{y_key}\n"
        f"Common data/MC scale; minimum {minimum_bin_count} entries per bin; "
        f"{cut_description}",
        fontsize=14,
        y=0.965,
    )

    output_dir = (
        output_base
        / "calorimeter"
        / cut_stage
        / "pcal"
        / "kinematics"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (
        f"{result.dataset.key}_pcal_electron_{pair_key}_"
        f"{KINEMATIC_FILENAME_KEYS[kinematic_key]}.png"
    )
    figure.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0.08,
    )
    plt.close(figure)
    return output_path


def mean_1d_values(
    sums: np.ndarray,
    counts: np.ndarray,
    minimum_count: int,
) -> np.ma.MaskedArray:
    valid = counts >= minimum_count
    means = np.zeros_like(sums, dtype=np.float64)
    np.divide(
        sums,
        counts,
        out=means,
        where=valid,
    )
    return np.ma.array(means, mask=~valid)


def save_pcal_electron_kinematic_1d_summary(
    result: PeriodResult,
    output_base: Path,
    histogram_config: HistogramConfig,
    coordinate_index: int,
    kinematic_index: int,
    strictness: int,
    is_rgc: bool,
    cut_stage: str,
    minimum_bin_count: int,
    dpi: int,
) -> Path:
    if cut_stage not in ("before_cuts", "after_cuts"):
        raise ValueError(f"Unknown calorimeter cut stage: {cut_stage}")
    #endif

    coordinate_key = COORD_KEYS[coordinate_index]
    kinematic_key = KINEMATIC_KEYS[kinematic_index]
    kinematic_label = KINEMATIC_LABELS[kinematic_key]

    data_counts = (
        result.data.pcal_electron_kin_1d_counts_before
        if cut_stage == "before_cuts"
        else result.data.pcal_electron_kin_1d_counts_after
    )
    data_sums = (
        result.data.pcal_electron_kin_1d_sums_before
        if cut_stage == "before_cuts"
        else result.data.pcal_electron_kin_1d_sums_after
    )

    mc_counts = None
    mc_sums = None
    if result.mc is not None:
        mc_counts = (
            result.mc.pcal_electron_kin_1d_counts_before
            if cut_stage == "before_cuts"
            else result.mc.pcal_electron_kin_1d_counts_after
        )
        mc_sums = (
            result.mc.pcal_electron_kin_1d_sums_before
            if cut_stage == "before_cuts"
            else result.mc.pcal_electron_kin_1d_sums_after
        )
    #endif

    edges = histogram_config.cal_edges
    centers = 0.5 * (edges[:-1] + edges[1:])

    figure, axes = plt.subplots(
        3,
        2,
        figsize=(13.0, 11.0),
        sharex=True,
        constrained_layout=True,
    )
    axes_flat = axes.ravel()

    for sector_index, axis in enumerate(axes_flat):
        sector = sector_index + 1

        data_mean = mean_1d_values(
            data_sums[
                kinematic_index,
                sector_index,
                coordinate_index,
            ],
            data_counts[sector_index, coordinate_index],
            minimum_bin_count,
        )
        axis.plot(
            centers,
            data_mean,
            marker="o",
            markersize=2.5,
            linewidth=1.1,
            label="Data",
        )

        if result.mc is not None:
            assert mc_counts is not None
            assert mc_sums is not None
            mc_mean = mean_1d_values(
                mc_sums[
                    kinematic_index,
                    sector_index,
                    coordinate_index,
                ],
                mc_counts[sector_index, coordinate_index],
                minimum_bin_count,
            )
            axis.plot(
                centers,
                mc_mean,
                marker="o",
                markersize=2.5,
                linewidth=1.1,
                label="MC",
            )
        #endif

        intervals = exclusion_intervals(
            dataset_key=result.dataset.key,
            layer_key="pcal",
            sector=sector,
            coordinate_key=coordinate_key,
            strictness=strictness,
            is_rgc=is_rgc,
        )
        draw_exclusions(axis, intervals)

        axis.set_title(f"Sector {sector}")
        axis.set_xlim(
            histogram_config.cal_min_cm,
            histogram_config.cal_max_cm,
        )
        axis.set_ylabel(kinematic_label)
        axis.grid(alpha=0.20)

        if sector_index >= 4:
            axis.set_xlabel(COORD_LABELS[coordinate_key])
        #endif

        if sector_index == 0:
            axis.legend()
        #endif
    #endfor

    cut_description = (
        "fiducial exclusions are displayed but not applied"
        if cut_stage == "before_cuts"
        else "complete calorimeter fiducial cuts are applied"
    )

    figure.suptitle(
        f"{result.dataset.label}: PCal electron {kinematic_key} "
        f"versus {coordinate_key}\n"
        f"4.5 cm strip bins; minimum {minimum_bin_count} entries per bin; "
        f"{cut_description}",
        fontsize=15,
    )

    output_dir = (
        output_base
        / "calorimeter"
        / cut_stage
        / "pcal"
        / "kinematics"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / (
        f"{result.dataset.key}_pcal_electron_{coordinate_key}_"
        f"{KINEMATIC_FILENAME_KEYS[kinematic_key]}_1d.png"
    )
    figure.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        pad_inches=0.08,
    )
    plt.close(figure)
    return output_path


# =============================================================================
# Main
# =============================================================================


def format_sample_report(
    label: str,
    sample: SampleResult,
) -> str:
    return (
        f"{label}: rows={sample.rows_read:,}, "
        f"PID-22 rows={sample.photon_rows:,}, "
        f"valid FT photons={sample.valid_ft_photons:,}, "
        f"FT fiducial={sample.fiducial_ft_photons:,}, "
        f"time={sample.elapsed_seconds:.1f} s"
    )


def main() -> int:
    args = parse_arguments()

    try:
        datasets, input_dir, mode_label = choose_datasets(args)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    #endtry

    if args.list_periods:
        print(f"{mode_label} period keys:")
        for dataset in datasets:
            print(f"  {dataset.key:<16} {dataset.label}")
        #endfor
        return 0
    #endif

    if not datasets:
        print("ERROR: no run periods selected.", file=sys.stderr)
        return 2
    #endif

    do_ft = not args.skip_ft
    do_calorimeter = not args.skip_calorimeter
    is_rgc = mode_label == "RGC"

    histogram_config = HistogramConfig(
        ft_bins=args.ft_bins,
        ft_xy_min_cm=args.ft_xy_min,
        ft_xy_max_cm=args.ft_xy_max,
        cal_bin_width_cm=args.cal_bin_width,
        cal_min_cm=args.cal_min,
        cal_max_cm=args.cal_max,
    )
    fiducial_config = FTFiducialConfig()

    effective_workers = min(
        MAX_WORKERS_HARD_LIMIT,
        args.workers,
        len(datasets),
    )

    print("=" * 80)
    print("CLAS12 calibration diagnostics")
    print("=" * 80)
    print(f"Mode:                   {mode_label}")
    print(f"Input directory:        {input_dir}")
    print(f"Output base:            {args.output_base}")
    print(f"Tree:                   {args.tree}")
    print(f"Periods:                {len(datasets)}")
    print(f"Workers:                {effective_workers}")
    print(f"Chunk size:             {args.chunk_size:,}")
    print(f"FT enabled:             {do_ft}")
    print(f"Calorimeter enabled:    {do_calorimeter}")
    print(f"Calorimeter bin width:  {args.cal_bin_width:.3f} cm")
    print(f"Calorimeter strictness: {args.calorimeter_strictness}")
    print(
        f"Kinematic min count:    {args.kinematic_min_bin_count}"
    )
    if args.max_events is not None:
        print(f"Per-file row limit:     {args.max_events:,}")
    #endif
    print("=" * 80)

    overall_start = time.perf_counter()
    completed: dict[str, PeriodResult] = {}
    failures: dict[str, str] = {}

    with ProcessPoolExecutor(max_workers=effective_workers) as executor:
        future_to_dataset = {
            executor.submit(
                process_period,
                dataset,
                input_dir,
                args.tree,
                histogram_config,
                fiducial_config,
                args.chunk_size,
                args.max_events,
                do_ft,
                do_calorimeter,
                args.calorimeter_strictness,
                is_rgc,
            ): dataset
            for dataset in datasets
        }

        for future in as_completed(future_to_dataset):
            dataset = future_to_dataset[future]

            try:
                result = future.result()
            except Exception as exc:
                failures[dataset.key] = (
                    f"{type(exc).__name__}: {exc}\n"
                    f"{traceback.format_exc()}"
                )
                print(
                    f"[FAILED] {dataset.label}: "
                    f"{type(exc).__name__}: {exc}",
                    file=sys.stderr,
                    flush=True,
                )
                continue
            #endtry

            completed[dataset.key] = result
            print(
                f"[READ] {dataset.label} completed by PID "
                f"{result.worker_pid} in {result.elapsed_seconds:.1f} s",
                flush=True,
            )
            print(
                f"       {format_sample_report('data', result.data)}",
                flush=True,
            )
            if result.mc is not None:
                print(
                    f"       {format_sample_report('MC', result.mc)}",
                    flush=True,
                )
            #endif
        #endfor
    #endwith

    # Plot in configured period order for deterministic filenames and logs.
    for dataset in datasets:
        result = completed.get(dataset.key)
        if result is None:
            continue
        #endif

        if do_ft:
            output_path = save_ft_summary(
                result,
                args.output_base,
                histogram_config,
                fiducial_config,
                args.dpi,
            )
            print(f"[PLOT] {output_path}", flush=True)
        #endif

        if do_calorimeter:
            for cut_stage in ("before_cuts", "after_cuts"):
                for coordinate_index in range(len(COORD_KEYS)):
                    for kinematic_index in range(len(KINEMATIC_KEYS)):
                        output_path = (
                            save_pcal_electron_kinematic_1d_summary(
                                result=result,
                                output_base=args.output_base,
                                histogram_config=histogram_config,
                                coordinate_index=coordinate_index,
                                kinematic_index=kinematic_index,
                                strictness=args.calorimeter_strictness,
                                is_rgc=is_rgc,
                                cut_stage=cut_stage,
                                minimum_bin_count=(
                                    args.kinematic_min_bin_count
                                ),
                                dpi=args.dpi,
                            )
                        )
                        print(f"[PLOT] {output_path}", flush=True)
                    #endfor
                #endfor

                for pair_slot in range(len(KINEMATIC_PAIR_INDICES)):
                    for kinematic_index in range(len(KINEMATIC_KEYS)):
                        output_path = save_pcal_electron_kinematic_summary(
                            result=result,
                            output_base=args.output_base,
                            histogram_config=histogram_config,
                            pair_slot=pair_slot,
                            kinematic_index=kinematic_index,
                            strictness=args.calorimeter_strictness,
                            is_rgc=is_rgc,
                            cut_stage=cut_stage,
                            minimum_bin_count=(
                                args.kinematic_min_bin_count
                            ),
                            dpi=args.dpi,
                        )
                        print(f"[PLOT] {output_path}", flush=True)
                    #endfor
                #endfor

                for layer_index in range(len(LAYER_KEYS)):
                    for particle_index in range(len(PARTICLE_KEYS)):
                        output_path = save_calorimeter_matching_summary(
                            result=result,
                            output_base=args.output_base,
                            histogram_config=histogram_config,
                            layer_index=layer_index,
                            particle_index=particle_index,
                            strictness=args.calorimeter_strictness,
                            is_rgc=is_rgc,
                            cut_stage=cut_stage,
                            dpi=args.dpi,
                        )
                        print(f"[PLOT] {output_path}", flush=True)

                        for pair_index in range(len(PAIR_DEFINITIONS)):
                            output_path = save_calorimeter_2d_matching_summary(
                                result=result,
                                output_base=args.output_base,
                                histogram_config=histogram_config,
                                layer_index=layer_index,
                                particle_index=particle_index,
                                pair_index=pair_index,
                                strictness=args.calorimeter_strictness,
                                is_rgc=is_rgc,
                                cut_stage=cut_stage,
                                dpi=args.dpi,
                            )
                            print(f"[PLOT] {output_path}", flush=True)
                        #endfor
                    #endfor
                #endfor
            #endfor
        #endif
    #endfor

    elapsed = time.perf_counter() - overall_start
    print("=" * 80)
    print(
        f"Completed {len(completed)} of {len(datasets)} period(s) "
        f"in {elapsed:.1f} s."
    )

    if failures:
        print(f"{len(failures)} period(s) failed:", file=sys.stderr)
        for dataset in datasets:
            failure = failures.get(dataset.key)
            if failure is None:
                continue
            #endif
            print(
                f"  {dataset.label}: {failure.splitlines()[0]}",
                file=sys.stderr,
            )
        #endfor
        return 1
    #endif

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
