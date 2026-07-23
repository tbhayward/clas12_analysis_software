#!/usr/bin/env python3
"""
plot_calibration_ft_v1.py

First-stage calibration diagnostics for the CLAS12 Forward Tagger (FT)
calorimeter fiducial selection.

The script is intentionally organized so that later detector-calibration and
fiducial-cut studies can be added as independent analysis modules.  Version 1
contains only the FT photon diagnostics adapted from the legacy calibration.cpp:

  1. FT photon hit occupancy before the fiducial selection.
  2. FT photon hit occupancy after the fiducial selection.
  3. Mean FT energy in each (x_FT, y_FT) bin.
  4. Mean FT energy with low-response bins highlighted.
  5. A four-panel overview containing all of the above.

By default, all five RGA run periods are processed.  The --rgc override instead
processes the three RGC enpi+ calibration data sets.  RGA data and DVCS MC are
both plotted; the currently available RGC inputs contain data only.

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
  python3 external_scripts/plot_calibration_ft_v1.py

  python3 external_scripts/plot_calibration_ft_v1.py --rgc

  python3 external_scripts/plot_calibration_ft_v1.py \
      --max-events 5000000 --chunk-size 250000

  python3 external_scripts/plot_calibration_ft_v1.py \
      --period rga_sp18_inb --period rga_fa18_inb

Dependencies
------------
  numpy, matplotlib, uproot
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

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
        "ERROR: uproot is required. Install it with `python3 -m pip install uproot`."
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

FT_BRANCHES = ("particle_pid", "ft_x", "ft_y", "ft_energy")
INVALID_SENTINEL = -9999.0
PHOTON_PID = 22


@dataclass(frozen=True)
class Dataset:
    """Input files and display metadata for one run period."""

    key: str
    label: str
    data_file: str
    mc_file: str | None = None


@dataclass(frozen=True)
class FTFiducialConfig:
    """Definition of the FT annulus and excluded circular regions."""

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
    """Binning and display range for all FT x-y histograms."""

    bins: int = 100
    xy_min_cm: float = -20.0
    xy_max_cm: float = 20.0

    @property
    def edges(self) -> np.ndarray:
        return np.linspace(self.xy_min_cm, self.xy_max_cm, self.bins + 1)


@dataclass
class FTHistograms:
    """Accumulated FT histograms and event counters for one input sample."""

    occupancy_all: np.ndarray
    occupancy_pass: np.ndarray
    energy_sum: np.ndarray
    energy_count: np.ndarray
    rows_read: int = 0
    photon_rows: int = 0
    valid_ft_photons: int = 0
    fiducial_ft_photons: int = 0

    @classmethod
    def empty(cls, bins: int) -> "FTHistograms":
        shape = (bins, bins)
        return cls(
            occupancy_all=np.zeros(shape, dtype=np.float64),
            occupancy_pass=np.zeros(shape, dtype=np.float64),
            energy_sum=np.zeros(shape, dtype=np.float64),
            energy_count=np.zeros(shape, dtype=np.float64),
        )

    @property
    def mean_energy(self) -> np.ndarray:
        result = np.full_like(self.energy_sum, np.nan, dtype=np.float64)
        np.divide(
            self.energy_sum,
            self.energy_count,
            out=result,
            where=self.energy_count > 0.0,
        )
        return result


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
# Command-line handling
# =============================================================================


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create FT calorimeter occupancy, mean-energy, and fiducial-cut "
            "diagnostics for all configured RGA periods or, with --rgc, RGC periods."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--rga",
        action="store_true",
        help="Process the five RGA DVCS data/MC pairs (the default mode).",
    )
    mode.add_argument(
        "--rgc",
        action="store_true",
        help="Process the three RGC enpi+ calibration data files instead.",
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=None,
        help="Override the default input directory for the selected mode.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=OUTPUT_DIR_DEFAULT,
        help="Top-level output directory for this FT study.",
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
            "Process only the requested period key. May be repeated. "
            "Use --list-periods to display valid keys."
        ),
    )
    parser.add_argument(
        "--list-periods",
        action="store_true",
        help="Print the period keys for the selected mode and exit.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500_000,
        help="Number of tree rows read per uproot chunk.",
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
        help="Minimum plotted x_FT and y_FT coordinate in cm.",
    )
    parser.add_argument(
        "--xy-max",
        type=float,
        default=20.0,
        help="Maximum plotted x_FT and y_FT coordinate in cm.",
    )
    parser.add_argument(
        "--minimum-energy-bin-count",
        type=int,
        default=1,
        help=(
            "Minimum photons required in an x-y bin before that bin enters the "
            "mean-energy statistics and low-response classification. A value of "
            "1 reproduces the legacy behavior."
        ),
    )
    parser.add_argument(
        "--low-response-sigma",
        type=float,
        default=1.0,
        help=(
            "Flag populated x-y bins whose mean energy lies this many standard "
            "deviations below the unweighted mean over populated bins."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=180,
        help="PNG resolution.",
    )

    args = parser.parse_args()

    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")
    if args.max_events is not None and args.max_events <= 0:
        parser.error("--max-events must be positive when supplied")
    if args.bins <= 0:
        parser.error("--bins must be positive")
    if args.xy_max <= args.xy_min:
        parser.error("--xy-max must be larger than --xy-min")
    if args.minimum_energy_bin_count <= 0:
        parser.error("--minimum-energy-bin-count must be positive")
    if args.low_response_sigma < 0.0:
        parser.error("--low-response-sigma cannot be negative")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")

    return args


# =============================================================================
# Input validation and histogram accumulation
# =============================================================================


def choose_datasets(args: argparse.Namespace) -> tuple[tuple[Dataset, ...], Path, str]:
    use_rgc = bool(args.rgc)
    datasets = RGC_DATASETS if use_rgc else RGA_DATASETS
    default_input_dir = RGC_INPUT_DIR_DEFAULT if use_rgc else RGA_INPUT_DIR_DEFAULT
    input_dir = args.input_dir if args.input_dir is not None else default_input_dir
    mode_label = "RGC" if use_rgc else "RGA"

    if args.period:
        available = {dataset.key: dataset for dataset in datasets}
        unknown = sorted(set(args.period) - set(available))
        if unknown:
            valid = ", ".join(available)
            raise ValueError(
                f"Unknown {mode_label} period key(s): {', '.join(unknown)}. "
                f"Valid keys: {valid}"
            )
        #endif
        requested = set(args.period)
        datasets = tuple(dataset for dataset in datasets if dataset.key in requested)
    #endif

    return datasets, input_dir, mode_label


def validate_root_tree(file_path: Path, tree_name: str) -> int:
    if not file_path.is_file():
        raise FileNotFoundError(f"Input ROOT file does not exist: {file_path}")
    #endif

    try:
        with uproot.open(file_path) as root_file:
            if tree_name not in root_file:
                available = ", ".join(str(key) for key in root_file.keys())
                raise KeyError(
                    f"Tree '{tree_name}' is absent from {file_path}. "
                    f"Available top-level keys: {available}"
                )
            #endif

            tree = root_file[tree_name]
            branch_names = set(tree.keys())
            missing = [name for name in FT_BRANCHES if name not in branch_names]
            if missing:
                raise KeyError(
                    f"Tree '{tree_name}' in {file_path} is missing required FT "
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
    radius_cm = np.hypot(x_cm, y_cm)
    accepted = (
        (radius_cm >= config.inner_radius_cm)
        & (radius_cm <= config.outer_radius_cm)
    )

    for center_x_cm, center_y_cm, hole_radius_cm in config.holes:
        distance_cm = np.hypot(x_cm - center_x_cm, y_cm - center_y_cm)
        accepted &= distance_cm >= hole_radius_cm
    #endfor

    return accepted


def accumulate_ft_histograms(
    file_path: Path,
    tree_name: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    chunk_size: int,
    max_events: int | None,
) -> FTHistograms:
    histograms = FTHistograms.empty(histogram_config.bins)
    edges = histogram_config.edges

    entry_stop = max_events
    iterator = uproot.iterate(
        f"{file_path}:{tree_name}",
        expressions=list(FT_BRANCHES),
        step_size=chunk_size,
        entry_stop=entry_stop,
        library="np",
    )

    for arrays in iterator:
        pid = np.asarray(arrays["particle_pid"])
        x_cm = np.asarray(arrays["ft_x"], dtype=np.float64)
        y_cm = np.asarray(arrays["ft_y"], dtype=np.float64)
        energy_gev = np.asarray(arrays["ft_energy"], dtype=np.float64)

        histograms.rows_read += int(pid.size)

        photon = pid == PHOTON_PID
        histograms.photon_rows += int(np.count_nonzero(photon))

        valid = (
            photon
            & np.isfinite(x_cm)
            & np.isfinite(y_cm)
            & np.isfinite(energy_gev)
            & (x_cm != INVALID_SENTINEL)
            & (y_cm != INVALID_SENTINEL)
            & (energy_gev != INVALID_SENTINEL)
        )

        if not np.any(valid):
            continue
        #endif

        valid_x = x_cm[valid]
        valid_y = y_cm[valid]
        valid_energy = energy_gev[valid]
        histograms.valid_ft_photons += int(valid_x.size)

        accepted = ft_fiducial_mask(valid_x, valid_y, fiducial_config)
        histograms.fiducial_ft_photons += int(np.count_nonzero(accepted))

        occupancy_all, _, _ = np.histogram2d(
            valid_x,
            valid_y,
            bins=(edges, edges),
        )
        occupancy_pass, _, _ = np.histogram2d(
            valid_x[accepted],
            valid_y[accepted],
            bins=(edges, edges),
        )
        energy_sum, _, _ = np.histogram2d(
            valid_x,
            valid_y,
            bins=(edges, edges),
            weights=valid_energy,
        )
        energy_count, _, _ = np.histogram2d(
            valid_x,
            valid_y,
            bins=(edges, edges),
        )

        histograms.occupancy_all += occupancy_all
        histograms.occupancy_pass += occupancy_pass
        histograms.energy_sum += energy_sum
        histograms.energy_count += energy_count
    #endfor

    return histograms


# =============================================================================
# Plotting
# =============================================================================


def draw_ft_boundaries(axis: plt.Axes, config: FTFiducialConfig) -> None:
    axis.add_patch(
        Circle(
            (0.0, 0.0),
            config.inner_radius_cm,
            fill=False,
            linewidth=1.6,
            linestyle="--",
            edgecolor="black",
        )
    )
    axis.add_patch(
        Circle(
            (0.0, 0.0),
            config.outer_radius_cm,
            fill=False,
            linewidth=1.6,
            linestyle="--",
            edgecolor="black",
        )
    )

    for center_x_cm, center_y_cm, radius_cm in config.holes:
        axis.add_patch(
            Circle(
                (center_x_cm, center_y_cm),
                radius_cm,
                fill=False,
                linewidth=1.6,
                edgecolor="black",
            )
        )
    #endfor


def configure_xy_axis(axis: plt.Axes, histogram_config: HistogramConfig) -> None:
    axis.set_xlabel(r"FT $x$ (cm)")
    axis.set_ylabel(r"FT $y$ (cm)")
    axis.set_xlim(histogram_config.xy_min_cm, histogram_config.xy_max_cm)
    axis.set_ylim(histogram_config.xy_min_cm, histogram_config.xy_max_cm)
    axis.set_aspect("equal", adjustable="box")


def nonempty_lognorm(histogram: np.ndarray) -> LogNorm | None:
    positive = histogram[histogram > 0.0]
    if positive.size == 0:
        return None
    #endif
    return LogNorm(vmin=max(1.0, float(np.min(positive))), vmax=float(np.max(positive)))


def draw_occupancy(
    figure: plt.Figure,
    axis: plt.Axes,
    histogram: np.ndarray,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
) -> None:
    extent = (
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
    )
    norm = nonempty_lognorm(histogram)
    image = axis.imshow(
        histogram.T,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        norm=norm,
        cmap="viridis",
    )
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("Photon entries per bin")
    axis.set_title(title)
    configure_xy_axis(axis, histogram_config)
    draw_ft_boundaries(axis, fiducial_config)


def energy_bin_statistics(
    histograms: FTHistograms,
    minimum_bin_count: int,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    mean_energy = histograms.mean_energy
    eligible = (
        np.isfinite(mean_energy)
        & (histograms.energy_count >= float(minimum_bin_count))
    )

    values = mean_energy[eligible]
    if values.size == 0:
        return mean_energy, eligible, float("nan"), float("nan")
    #endif

    global_mean = float(np.mean(values))
    global_std = float(np.std(values, ddof=0))
    return mean_energy, eligible, global_mean, global_std


def draw_mean_energy(
    figure: plt.Figure,
    axis: plt.Axes,
    histograms: FTHistograms,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    minimum_bin_count: int,
    low_response_sigma: float,
    highlight_low_response: bool,
) -> tuple[float, float, int]:
    mean_energy, eligible, global_mean, global_std = energy_bin_statistics(
        histograms,
        minimum_bin_count,
    )

    displayed = np.ma.masked_where(~eligible, mean_energy)
    extent = (
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
        histogram_config.xy_min_cm,
        histogram_config.xy_max_cm,
    )
    image = axis.imshow(
        displayed.T,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        cmap="viridis",
    )
    colorbar = figure.colorbar(image, ax=axis, pad=0.02)
    colorbar.set_label("Mean FT energy (GeV)")

    low_response = np.zeros_like(eligible, dtype=bool)
    if np.isfinite(global_mean) and np.isfinite(global_std):
        threshold = global_mean - low_response_sigma * global_std
        low_response = eligible & (mean_energy < threshold)

        if highlight_low_response and np.any(low_response):
            low_overlay = np.ma.masked_where(~low_response, np.ones_like(mean_energy))
            axis.imshow(
                low_overlay.T,
                origin="lower",
                extent=extent,
                interpolation="nearest",
                cmap="Reds",
                alpha=0.72,
                vmin=0.0,
                vmax=1.0,
            )
        #endif
    #endif

    axis.set_title(title)
    configure_xy_axis(axis, histogram_config)
    draw_ft_boundaries(axis, fiducial_config)

    if np.isfinite(global_mean):
        annotation = (
            f"Populated-bin mean = {global_mean:.3f} GeV\n"
            f"Populated-bin std. dev. = {global_std:.3f} GeV\n"
            f"Eligible bins = {np.count_nonzero(eligible):,}"
        )
        if highlight_low_response:
            annotation += f"\nFlagged bins = {np.count_nonzero(low_response):,}"
        #endif
        axis.text(
            0.02,
            0.98,
            annotation,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "none"},
        )
    #endif

    return global_mean, global_std, int(np.count_nonzero(low_response))


def save_single_occupancy_plot(
    output_path: Path,
    histogram: np.ndarray,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    dpi: int,
) -> None:
    figure, axis = plt.subplots(figsize=(7.4, 6.4), constrained_layout=True)
    draw_occupancy(
        figure,
        axis,
        histogram,
        title,
        histogram_config,
        fiducial_config,
    )
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def save_single_energy_plot(
    output_path: Path,
    histograms: FTHistograms,
    title: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    minimum_bin_count: int,
    low_response_sigma: float,
    highlight_low_response: bool,
    dpi: int,
) -> tuple[float, float, int]:
    figure, axis = plt.subplots(figsize=(7.4, 6.4), constrained_layout=True)
    statistics = draw_mean_energy(
        figure,
        axis,
        histograms,
        title,
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        highlight_low_response,
    )
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)
    return statistics


def save_overview_plot(
    output_path: Path,
    histograms: FTHistograms,
    period_label: str,
    sample_label: str,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    minimum_bin_count: int,
    low_response_sigma: float,
    dpi: int,
) -> None:
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(13.0, 11.0),
        constrained_layout=True,
    )
    figure.suptitle(f"{period_label}: Forward Tagger photon diagnostics ({sample_label})")

    draw_occupancy(
        figure,
        axes[0, 0],
        histograms.occupancy_all,
        "FT hit occupancy before fiducial cut",
        histogram_config,
        fiducial_config,
    )
    draw_occupancy(
        figure,
        axes[0, 1],
        histograms.occupancy_pass,
        "FT hit occupancy after fiducial cut",
        histogram_config,
        fiducial_config,
    )
    draw_mean_energy(
        figure,
        axes[1, 0],
        histograms,
        "Mean FT energy by hit position",
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        False,
    )
    draw_mean_energy(
        figure,
        axes[1, 1],
        histograms,
        (
            "Low-response bins: mean energy below "
            f"global mean - {low_response_sigma:g} sigma"
        ),
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        True,
    )

    acceptance = (
        histograms.fiducial_ft_photons / histograms.valid_ft_photons
        if histograms.valid_ft_photons > 0
        else float("nan")
    )
    figure.text(
        0.5,
        0.005,
        (
            f"Rows read: {histograms.rows_read:,}    |    "
            f"PID 22 rows: {histograms.photon_rows:,}    |    "
            f"Valid FT photons: {histograms.valid_ft_photons:,}    |    "
            f"Fiducial-pass photons: {histograms.fiducial_ft_photons:,} "
            f"({100.0 * acceptance:.2f}%)"
        ),
        ha="center",
        va="bottom",
        fontsize=10,
    )
    figure.savefig(output_path, dpi=dpi)
    plt.close(figure)


def write_summary_file(
    output_path: Path,
    dataset: Dataset,
    sample_label: str,
    input_path: Path,
    total_tree_entries: int,
    histograms: FTHistograms,
    global_mean_energy: float,
    global_std_energy: float,
    low_response_bins: int,
    minimum_bin_count: int,
    low_response_sigma: float,
    fiducial_config: FTFiducialConfig,
) -> None:
    acceptance = (
        histograms.fiducial_ft_photons / histograms.valid_ft_photons
        if histograms.valid_ft_photons > 0
        else float("nan")
    )

    lines = [
        f"Period: {dataset.label}",
        f"Sample: {sample_label}",
        f"Input file: {input_path}",
        f"Tree entries: {total_tree_entries}",
        f"Rows read: {histograms.rows_read}",
        f"PID 22 rows: {histograms.photon_rows}",
        f"Valid FT photons: {histograms.valid_ft_photons}",
        f"Fiducial-pass FT photons: {histograms.fiducial_ft_photons}",
        f"Fiducial-pass fraction: {acceptance:.10f}",
        f"Mean over eligible mean-energy bins (GeV): {global_mean_energy:.10g}",
        f"Std. dev. over eligible mean-energy bins (GeV): {global_std_energy:.10g}",
        f"Low-response bins: {low_response_bins}",
        f"Minimum energy-bin count: {minimum_bin_count}",
        f"Low-response threshold (sigma): {low_response_sigma}",
        f"FT inner radius (cm): {fiducial_config.inner_radius_cm}",
        f"FT outer radius (cm): {fiducial_config.outer_radius_cm}",
        "Excluded holes (x_cm, y_cm, radius_cm):",
    ]
    lines.extend(
        f"  {x_cm:.6g}, {y_cm:.6g}, {radius_cm:.6g}"
        for x_cm, y_cm, radius_cm in fiducial_config.holes
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def make_sample_plots(
    dataset: Dataset,
    sample_label: str,
    input_path: Path,
    total_tree_entries: int,
    histograms: FTHistograms,
    period_output_dir: Path,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    minimum_bin_count: int,
    low_response_sigma: float,
    dpi: int,
) -> None:
    sample_key = sample_label.lower()
    sample_output_dir = period_output_dir / sample_key
    sample_output_dir.mkdir(parents=True, exist_ok=True)

    prefix = f"{dataset.key}_{sample_key}_ft"

    save_single_occupancy_plot(
        sample_output_dir / f"{prefix}_hit_occupancy_before_cut.png",
        histograms.occupancy_all,
        f"{dataset.label}: FT photon hit occupancy ({sample_label}, before cut)",
        histogram_config,
        fiducial_config,
        dpi,
    )
    save_single_occupancy_plot(
        sample_output_dir / f"{prefix}_hit_occupancy_after_cut.png",
        histograms.occupancy_pass,
        f"{dataset.label}: FT photon hit occupancy ({sample_label}, after cut)",
        histogram_config,
        fiducial_config,
        dpi,
    )
    global_mean, global_std, _ = save_single_energy_plot(
        sample_output_dir / f"{prefix}_mean_energy.png",
        histograms,
        f"{dataset.label}: mean FT energy ({sample_label})",
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        False,
        dpi,
    )
    _, _, low_response_bins = save_single_energy_plot(
        sample_output_dir / f"{prefix}_mean_energy_low_response.png",
        histograms,
        (
            f"{dataset.label}: low-response FT bins ({sample_label}; "
            f"below mean - {low_response_sigma:g} sigma)"
        ),
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        True,
        dpi,
    )
    save_overview_plot(
        sample_output_dir / f"{prefix}_overview.png",
        histograms,
        dataset.label,
        sample_label,
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        dpi,
    )
    write_summary_file(
        sample_output_dir / f"{prefix}_summary.txt",
        dataset,
        sample_label,
        input_path,
        total_tree_entries,
        histograms,
        global_mean,
        global_std,
        low_response_bins,
        minimum_bin_count,
        low_response_sigma,
        fiducial_config,
    )


# =============================================================================
# Main driver
# =============================================================================


def process_sample(
    dataset: Dataset,
    sample_label: str,
    input_path: Path,
    tree_name: str,
    output_dir: Path,
    histogram_config: HistogramConfig,
    fiducial_config: FTFiducialConfig,
    chunk_size: int,
    max_events: int | None,
    minimum_bin_count: int,
    low_response_sigma: float,
    dpi: int,
) -> None:
    total_entries = validate_root_tree(input_path, tree_name)
    rows_to_read = min(total_entries, max_events) if max_events is not None else total_entries

    print(
        f"  [{sample_label}] {input_path.name}: "
        f"reading {rows_to_read:,} of {total_entries:,} tree rows"
    )

    histograms = accumulate_ft_histograms(
        input_path,
        tree_name,
        histogram_config,
        fiducial_config,
        chunk_size,
        max_events,
    )

    period_output_dir = output_dir / dataset.key
    make_sample_plots(
        dataset,
        sample_label,
        input_path,
        total_entries,
        histograms,
        period_output_dir,
        histogram_config,
        fiducial_config,
        minimum_bin_count,
        low_response_sigma,
        dpi,
    )

    fraction = (
        100.0 * histograms.fiducial_ft_photons / histograms.valid_ft_photons
        if histograms.valid_ft_photons > 0
        else float("nan")
    )
    print(
        f"    valid FT photons = {histograms.valid_ft_photons:,}; "
        f"fiducial pass = {histograms.fiducial_ft_photons:,} ({fraction:.2f}%)"
    )


def main() -> int:
    args = parse_arguments()

    try:
        datasets, input_dir, mode_label = choose_datasets(args)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    if args.list_periods:
        print(f"Available {mode_label} period keys:")
        for dataset in datasets:
            print(f"  {dataset.key:<16} {dataset.label}")
        #endfor
        return 0
    #endif

    histogram_config = HistogramConfig(
        bins=args.bins,
        xy_min_cm=args.xy_min,
        xy_max_cm=args.xy_max,
    )
    fiducial_config = FTFiducialConfig()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"FT calibration mode: {mode_label}")
    print(f"Input directory: {input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Tree: {args.tree}")
    print(f"Periods: {', '.join(dataset.label for dataset in datasets)}")

    failures: list[str] = []

    for dataset in datasets:
        print(f"\nProcessing {dataset.label}")

        data_path = input_dir / dataset.data_file
        try:
            process_sample(
                dataset=dataset,
                sample_label="Data",
                input_path=data_path,
                tree_name=args.tree,
                output_dir=args.output_dir,
                histogram_config=histogram_config,
                fiducial_config=fiducial_config,
                chunk_size=args.chunk_size,
                max_events=args.max_events,
                minimum_bin_count=args.minimum_energy_bin_count,
                low_response_sigma=args.low_response_sigma,
                dpi=args.dpi,
            )
        except Exception as exc:
            message = f"{dataset.label} data: {exc}"
            failures.append(message)
            print(f"  ERROR: {message}", file=sys.stderr)
        #endtry

        if dataset.mc_file is not None:
            mc_path = input_dir / dataset.mc_file
            try:
                process_sample(
                    dataset=dataset,
                    sample_label="MC",
                    input_path=mc_path,
                    tree_name=args.tree,
                    output_dir=args.output_dir,
                    histogram_config=histogram_config,
                    fiducial_config=fiducial_config,
                    chunk_size=args.chunk_size,
                    max_events=args.max_events,
                    minimum_bin_count=args.minimum_energy_bin_count,
                    low_response_sigma=args.low_response_sigma,
                    dpi=args.dpi,
                )
            except Exception as exc:
                message = f"{dataset.label} MC: {exc}"
                failures.append(message)
                print(f"  ERROR: {message}", file=sys.stderr)
            #endtry
        #endif
    #endfor

    if failures:
        print("\nFT calibration completed with failures:", file=sys.stderr)
        for failure in failures:
            print(f"  - {failure}", file=sys.stderr)
        #endfor
        return 1
    #endif

    print("\nFT calibration plots completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
#endif
