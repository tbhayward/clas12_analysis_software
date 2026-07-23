#!/usr/bin/env python3
"""
plot_calibration_ft_v2.py

Fast, streamlined Forward Tagger photon-occupancy diagnostics for the CLAS12
calibration ROOT trees.

For each run period, the script creates exactly one 2x2 summary figure:

    top left:     data before the FT fiducial cut
    top right:    data after the FT fiducial cut
    bottom left:  reconstructed MC before the FT fiducial cut
    bottom right: reconstructed MC after the FT fiducial cut

The five RGA run periods are processed by default.  The --rgc override instead
processes the three RGC data sets.  Since no corresponding RGC MC files are
configured, the lower two panels are marked "MC not available."

Run periods are parallelized with at most five worker processes.  Each worker
handles one complete period, reading its data and MC files sequentially.  This
keeps ROOT-file access simple while using all available period-level
parallelism.

Only these branches are read:

    particle_pid
    ft_x
    ft_y

Default input locations
-----------------------
RGA:
  /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/dvcs_calibration/

RGC:
  /work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/

Default output location
-----------------------
  output/calibration/forward_tagger/

Examples
--------
Process all five RGA periods with up to five workers:

  python3 external_scripts/plot_calibration_ft_v2.py

Process the three RGC periods:

  python3 external_scripts/plot_calibration_ft_v2.py --rgc

Process selected periods only:

  python3 external_scripts/plot_calibration_ft_v2.py \
      --period rga_sp18_inb \
      --period rga_fa18_inb

Limit processing for a quick test:

  python3 external_scripts/plot_calibration_ft_v2.py \
      --max-events 5000000 \
      --workers 2
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
from typing import Sequence

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
# Configuration
# =============================================================================

TREE_NAME_DEFAULT = "PhysicsEvents"

RGA_INPUT_DIR_DEFAULT = Path(
    "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/dvcs_calibration"
)
RGC_INPUT_DIR_DEFAULT = Path(
    "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration"
)
OUTPUT_DIR_DEFAULT = Path("output/calibration/forward_tagger")

FT_BRANCHES = ("particle_pid", "ft_x", "ft_y")
PHOTON_PID = 22
INVALID_SENTINEL = -9999.0
MAX_WORKERS_HARD_LIMIT = 5


@dataclass(frozen=True)
class Dataset:
    """Input files and display metadata for one run period."""

    key: str
    label: str
    data_file: str
    mc_file: str | None = None


@dataclass(frozen=True)
class FTFiducialConfig:
    """FT annulus and excluded circular regions."""

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
    """FT coordinate histogram binning."""

    bins: int = 100
    xy_min_cm: float = -20.0
    xy_max_cm: float = 20.0

    @property
    def edges(self) -> np.ndarray:
        return np.linspace(
            self.xy_min_cm,
            self.xy_max_cm,
            self.bins + 1,
            dtype=np.float64,
        )


@dataclass
class OccupancyResult:
    """Accumulated before/after occupancy and counters for one ROOT file."""

    before: np.ndarray
    after: np.ndarray
    rows_read: int
    photon_rows: int
    valid_ft_photons: int
    fiducial_ft_photons: int
    elapsed_seconds: float


@dataclass
class PeriodResult:
    """Data and optional MC results for one run period."""

    dataset: Dataset
    data: OccupancyResult
    mc: OccupancyResult | None
    worker_pid: int
    elapsed_seconds: float


RGA_DATASETS: tuple[Dataset, ...] = (
    Dataset(
        key="rga_sp18_inb",
        label="RGA Sp18 Inb",
        data_file="DVCSWagon_rga_sp18_inb.root",
        mc_file="dvcsgen_rga_sp18_inb.root",
    ),
    Dataset(
        key="rga_sp18_out",
        label="RGA Sp18 Out",
        data_file="DVCSWagon_rga_sp18_out.root",
        mc_file="dvcsgen_rga_sp18_out.root",
    ),
    Dataset(
        key="rga_fa18_inb",
        label="RGA Fa18 Inb",
        data_file="DVCSWagon_rga_fa18_inb.root",
        mc_file="dvcsgen_rga_fa18_inb.root",
    ),
    Dataset(
        key="rga_fa18_out",
        label="RGA Fa18 Out",
        data_file="DVCSWagon_rga_fa18_out.root",
        mc_file="dvcsgen_rga_fa18_out.root",
    ),
    Dataset(
        key="rga_sp19_inb",
        label="RGA Sp19 Inb",
        data_file="DVCSWagon_rga_sp19_inb.root",
        mc_file="dvcsgen_rga_sp19_inb.root",
    ),
)

RGC_DATASETS: tuple[Dataset, ...] = (
    Dataset(
        key="rgc_su22_inb",
        label="RGC Su22 Inb",
        data_file="rgc_su22_inb_NH3_epi+X_calibration.root",
    ),
    Dataset(
        key="rgc_fa22_inb",
        label="RGC Fa22 Inb",
        data_file="rgc_fa22_inb_NH3_epi+X_calibration.root",
    ),
    Dataset(
        key="rgc_sp23_inb",
        label="RGC Sp23 Inb",
        data_file="rgc_sp23_inb_NH3_epi+X_calibration.root",
    ),
)


# =============================================================================
# Command-line interface
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create one four-panel FT occupancy summary for each configured "
            "RGA or RGC run period."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--rga",
        action="store_true",
        help="Process the five RGA data/MC periods; this is the default.",
    )
    mode.add_argument(
        "--rgc",
        action="store_true",
        help="Process the three RGC data periods instead.",
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Override the input directory for the selected mode.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR_DEFAULT,
        help="Single directory receiving all FT summary figures.",
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
        help=(
            "Process only this period key. May be repeated. "
            "Use --list-periods to display valid keys."
        ),
    )
    parser.add_argument(
        "--list-periods",
        action="store_true",
        help="Print valid period keys for the selected mode and exit.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=MAX_WORKERS_HARD_LIMIT,
        help=(
            "Requested number of period-level worker processes. The effective "
            "value is capped at five and at the number of selected periods."
        ),
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="ROOT rows read per uproot chunk.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Optional maximum rows read from each data or MC file.",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=100,
        help="Number of bins along each FT coordinate axis.",
    )
    parser.add_argument(
        "--xy-min",
        type=float,
        default=-20.0,
        help="Minimum plotted FT coordinate in cm.",
    )
    parser.add_argument(
        "--xy-max",
        type=float,
        default=20.0,
        help="Maximum plotted FT coordinate in cm.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="PNG resolution.",
    )

    args = parser.parse_args()

    if args.workers <= 0:
        parser.error("--workers must be positive")
    #endif
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")
    #endif
    if args.max_events is not None and args.max_events <= 0:
        parser.error("--max-events must be positive when supplied")
    #endif
    if args.bins <= 0:
        parser.error("--bins must be positive")
    #endif
    if args.xy_max <= args.xy_min:
        parser.error("--xy-max must exceed --xy-min")
    #endif
    if args.dpi <= 0:
        parser.error("--dpi must be positive")
    #endif

    return args


def choose_datasets(
    args: argparse.Namespace,
) -> tuple[tuple[Dataset, ...], Path, str]:
    use_rgc = bool(args.rgc)
    datasets = RGC_DATASETS if use_rgc else RGA_DATASETS
    default_input_dir = RGC_INPUT_DIR_DEFAULT if use_rgc else RGA_INPUT_DIR_DEFAULT
    input_dir = args.input_dir if args.input_dir is not None else default_input_dir
    mode_label = "RGC" if use_rgc else "RGA"

    if args.period:
        by_key = {dataset.key: dataset for dataset in datasets}
        unknown = sorted(set(args.period) - set(by_key))
        if unknown:
            valid = ", ".join(by_key)
            raise ValueError(
                f"Unknown {mode_label} period key(s): {', '.join(unknown)}. "
                f"Valid keys: {valid}"
            )
        #endif

        requested = set(args.period)
        datasets = tuple(
            dataset for dataset in datasets if dataset.key in requested
        )
    #endif

    return datasets, input_dir, mode_label


# =============================================================================
# Fast occupancy accumulation
# =============================================================================


def validate_root_tree(file_path: Path, tree_name: str) -> int:
    """Validate one ROOT input and return its tree entry count."""

    if not file_path.is_file():
        raise FileNotFoundError(f"Input ROOT file does not exist: {file_path}")
    #endif

    try:
        with uproot.open(file_path) as root_file:
            if tree_name not in root_file:
                keys = ", ".join(str(key) for key in root_file.keys())
                raise KeyError(
                    f"Tree '{tree_name}' is absent from {file_path}. "
                    f"Available top-level keys: {keys}"
                )
            #endif

            tree = root_file[tree_name]
            available_branches = set(tree.keys())
            missing = [
                branch for branch in FT_BRANCHES
                if branch not in available_branches
            ]
            if missing:
                raise KeyError(
                    f"Tree '{tree_name}' in {file_path} is missing required "
                    f"branch(es): {', '.join(missing)}"
                )
            #endif

            return int(tree.num_entries)
    except OSError as exc:
        raise OSError(f"Could not open ROOT file {file_path}: {exc}") from exc


def ft_fiducial_mask(
    x_cm: np.ndarray,
    y_cm: np.ndarray,
    config: FTFiducialConfig,
) -> np.ndarray:
    """
    Return the FT annulus-plus-hole mask.

    Squared distances are used throughout to avoid unnecessary square roots.
    """

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


def accumulate_occupancy(
    file_path: Path,
    tree_name: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    chunk_size: int,
    max_events: int | None,
) -> OccupancyResult:
    """Read one ROOT file and accumulate only the required occupancy maps."""

    start_time = time.perf_counter()
    validate_root_tree(file_path, tree_name)

    edges = histogram_config.edges
    shape = (histogram_config.bins, histogram_config.bins)
    before = np.zeros(shape, dtype=np.uint64)
    after = np.zeros(shape, dtype=np.uint64)

    rows_read = 0
    photon_rows = 0
    valid_ft_photons = 0
    fiducial_ft_photons = 0

    iterator = uproot.iterate(
        f"{file_path}:{tree_name}",
        expressions=list(FT_BRANCHES),
        step_size=chunk_size,
        entry_stop=max_events,
        library="np",
    )

    for arrays in iterator:
        pid = np.asarray(arrays["particle_pid"])
        x_cm = np.asarray(arrays["ft_x"], dtype=np.float64)
        y_cm = np.asarray(arrays["ft_y"], dtype=np.float64)

        rows_read += int(pid.size)

        photon_mask = pid == PHOTON_PID
        photon_rows += int(np.count_nonzero(photon_mask))

        valid_mask = (
            photon_mask
            & np.isfinite(x_cm)
            & np.isfinite(y_cm)
            & (x_cm != INVALID_SENTINEL)
            & (y_cm != INVALID_SENTINEL)
        )

        if not np.any(valid_mask):
            continue
        #endif

        valid_x = x_cm[valid_mask]
        valid_y = y_cm[valid_mask]
        valid_ft_photons += int(valid_x.size)

        accepted_mask = ft_fiducial_mask(
            valid_x,
            valid_y,
            fiducial_config,
        )
        fiducial_ft_photons += int(np.count_nonzero(accepted_mask))

        chunk_before, _, _ = np.histogram2d(
            valid_x,
            valid_y,
            bins=(edges, edges),
        )
        before += chunk_before.astype(np.uint64, copy=False)

        if np.any(accepted_mask):
            chunk_after, _, _ = np.histogram2d(
                valid_x[accepted_mask],
                valid_y[accepted_mask],
                bins=(edges, edges),
            )
            after += chunk_after.astype(np.uint64, copy=False)
        #endif
    #endfor

    return OccupancyResult(
        before=before,
        after=after,
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
) -> PeriodResult:
    """Worker entry point: process one complete run period."""

    start_time = time.perf_counter()

    data_result = accumulate_occupancy(
        file_path=input_dir / dataset.data_file,
        tree_name=tree_name,
        histogram_config=histogram_config,
        fiducial_config=fiducial_config,
        chunk_size=chunk_size,
        max_events=max_events,
    )

    mc_result = None
    if dataset.mc_file is not None:
        mc_result = accumulate_occupancy(
            file_path=input_dir / dataset.mc_file,
            tree_name=tree_name,
            histogram_config=histogram_config,
            fiducial_config=fiducial_config,
            chunk_size=chunk_size,
            max_events=max_events,
        )
    #endif

    return PeriodResult(
        dataset=dataset,
        data=data_result,
        mc=mc_result,
        worker_pid=os.getpid(),
        elapsed_seconds=time.perf_counter() - start_time,
    )


# =============================================================================
# Plotting
# =============================================================================


def draw_ft_boundaries(
    axis: plt.Axes,
    config: FTFiducialConfig,
) -> None:
    """Overlay the complete FT fiducial geometry."""

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


def configure_xy_axis(
    axis: plt.Axes,
    histogram_config: HistogramConfig,
) -> None:
    axis.set_xlabel(r"FT $x$ (cm)")
    axis.set_ylabel(r"FT $y$ (cm)")
    axis.set_xlim(
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
    )
    axis.set_ylim(
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
    )
    axis.set_aspect("equal", adjustable="box")


def shared_log_norm(
    first: np.ndarray,
    second: np.ndarray,
) -> LogNorm | None:
    """
    Construct one logarithmic normalization for a before/after row.

    Sharing the normalization makes the occupancy loss visually meaningful.
    """

    maximum = max(float(np.max(first)), float(np.max(second)))
    if maximum <= 0.0:
        return None
    #endif

    return LogNorm(vmin=1.0, vmax=max(1.0, maximum))


def draw_occupancy_panel(
    axis: plt.Axes,
    histogram: np.ndarray,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    norm: LogNorm | None,
):
    extent = (
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
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
    configure_xy_axis(axis, histogram_config)
    draw_ft_boundaries(axis, fiducial_config)
    return image


def efficiency_text(result: OccupancyResult) -> str:
    if result.valid_ft_photons <= 0:
        return "0 / 0 valid FT photons"
    #endif

    efficiency_percent = (
        100.0
        * result.fiducial_ft_photons
        / result.valid_ft_photons
    )
    return (
        f"{result.fiducial_ft_photons:,} / "
        f"{result.valid_ft_photons:,} valid FT photons "
        f"({efficiency_percent:.2f}%)"
    )


def draw_missing_mc_panel(
    axis: plt.Axes,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
) -> None:
    axis.set_title(title)
    configure_xy_axis(axis, histogram_config)
    draw_ft_boundaries(axis, fiducial_config)
    axis.text(
        0.5,
        0.5,
        "MC not available\nfor this RGC mode",
        transform=axis.transAxes,
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=15,
    )


def save_period_summary(
    result: PeriodResult,
    output_dir: Path,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    dpi: int,
) -> Path:
    """Create the only plot emitted for one run period."""

    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.5, 11.5),
        constrained_layout=True,
    )

    data_norm = shared_log_norm(
        result.data.before,
        result.data.after,
    )

    data_before_image = draw_occupancy_panel(
        axis=axes[0, 0],
        histogram=result.data.before,
        title="Data before FT fiducial cut",
        histogram_config=histogram_config,
        fiducial_config=fiducial_config,
        norm=data_norm,
    )
    draw_occupancy_panel(
        axis=axes[0, 1],
        histogram=result.data.after,
        title="Data after FT fiducial cut",
        histogram_config=histogram_config,
        fiducial_config=fiducial_config,
        norm=data_norm,
    )

    data_colorbar = figure.colorbar(
        data_before_image,
        ax=axes[0, :],
        pad=0.015,
        fraction=0.035,
    )
    data_colorbar.set_label("Photon entries per bin")

    if result.mc is not None:
        mc_norm = shared_log_norm(
            result.mc.before,
            result.mc.after,
        )

        mc_before_image = draw_occupancy_panel(
            axis=axes[1, 0],
            histogram=result.mc.before,
            title="MC before FT fiducial cut",
            histogram_config=histogram_config,
            fiducial_config=fiducial_config,
            norm=mc_norm,
        )
        draw_occupancy_panel(
            axis=axes[1, 1],
            histogram=result.mc.after,
            title="MC after FT fiducial cut",
            histogram_config=histogram_config,
            fiducial_config=fiducial_config,
            norm=mc_norm,
        )

        mc_colorbar = figure.colorbar(
            mc_before_image,
            ax=axes[1, :],
            pad=0.015,
            fraction=0.035,
        )
        mc_colorbar.set_label("Photon entries per bin")
    else:
        draw_missing_mc_panel(
            axis=axes[1, 0],
            title="MC before FT fiducial cut",
            histogram_config=histogram_config,
            fiducial_config=fiducial_config,
        )
        draw_missing_mc_panel(
            axis=axes[1, 1],
            title="MC after FT fiducial cut",
            histogram_config=histogram_config,
            fiducial_config=fiducial_config,
        )
    #endif

    figure.suptitle(
        f"{result.dataset.label}: Forward Tagger photon occupancy",
        fontsize=17,
    )

    data_caption = f"Data: {efficiency_text(result.data)}"
    if result.mc is not None:
        mc_caption = f"MC: {efficiency_text(result.mc)}"
    else:
        mc_caption = "MC: not available"
    #endif

    figure.text(
        0.5,
        0.006,
        f"{data_caption}    |    {mc_caption}",
        horizontalalignment="center",
        verticalalignment="bottom",
        fontsize=10,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{result.dataset.key}_ft_summary.png"
    figure.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(figure)

    return output_path


# =============================================================================
# Main execution
# =============================================================================


def format_sample_report(
    sample_label: str,
    result: OccupancyResult,
) -> str:
    return (
        f"{sample_label}: rows={result.rows_read:,}, "
        f"PID-22 rows={result.photon_rows:,}, "
        f"valid FT photons={result.valid_ft_photons:,}, "
        f"fiducial={result.fiducial_ft_photons:,}, "
        f"time={result.elapsed_seconds:.1f} s"
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
        print("ERROR: no run periods were selected.", file=sys.stderr)
        return 2
    #endif

    histogram_config = HistogramConfig(
        bins=args.bins,
        xy_min_cm=args.xy_min,
        xy_max_cm=args.xy_max,
    )
    fiducial_config = FTFiducialConfig()

    effective_workers = min(
        MAX_WORKERS_HARD_LIMIT,
        args.workers,
        len(datasets),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print("Forward Tagger occupancy diagnostics")
    print("=" * 78)
    print(f"Mode:             {mode_label}")
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Tree:             {args.tree}")
    print(f"Periods:          {len(datasets)}")
    print(f"Workers:          {effective_workers}")
    print(f"Chunk size:       {args.chunk_size:,} rows")
    if args.max_events is not None:
        print(f"Per-file limit:   {args.max_events:,} rows")
    #endif
    print("=" * 78)

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
                    f"[FAILED] {dataset.label}: {type(exc).__name__}: {exc}",
                    file=sys.stderr,
                    flush=True,
                )
                continue
            #endtry

            completed[dataset.key] = result
            print(
                f"[READ] {dataset.label} completed by PID {result.worker_pid} "
                f"in {result.elapsed_seconds:.1f} s",
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

    # Plot in configured run-period order for deterministic output and logs.
    for dataset in datasets:
        result = completed.get(dataset.key)
        if result is None:
            continue
        #endif

        output_path = save_period_summary(
            result=result,
            output_dir=args.output_dir,
            histogram_config=histogram_config,
            fiducial_config=fiducial_config,
            dpi=args.dpi,
        )
        print(f"[PLOT] {output_path}", flush=True)
    #endfor

    elapsed = time.perf_counter() - overall_start
    print("=" * 78)
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
            first_line = failure.splitlines()[0]
            print(
                f"  {dataset.label}: {first_line}",
                file=sys.stderr,
            )
        #endfor
        return 1
    #endif

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
