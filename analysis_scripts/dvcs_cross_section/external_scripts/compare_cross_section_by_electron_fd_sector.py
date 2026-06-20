#!/usr/bin/env python3
"""
compare_cross_section_by_sector.py

Compare integrated DVCS cross sections from separate sector-resolved CSV files.

This is a generic sector-comparison wrapper around the shared plotting engine in:

  compare_cross_section_by_topology.py

It can be used for electron, photon, or proton sector studies.

Examples:

  Electron FD sectors:

    python3 compare_cross_section_by_sector.py \
      elec_sec1.csv elec_sec2.csv elec_sec3.csv \
      elec_sec4.csv elec_sec5.csv elec_sec6.csv \
      --particle electron

  Photon FD sectors:

    python3 compare_cross_section_by_sector.py \
      photon_sec1.csv photon_sec2.csv photon_sec3.csv \
      photon_sec4.csv photon_sec5.csv photon_sec6.csv \
      --particle photon

  Proton sectors:

    python3 compare_cross_section_by_sector.py \
      proton_sec1.csv proton_sec2.csv proton_sec3.csv \
      proton_sec4.csv proton_sec5.csv proton_sec6.csv \
      --particle proton

  Proton CD sectors, if only three sector files exist:

    python3 compare_cross_section_by_sector.py \
      proton_cd_sec1.csv proton_cd_sec2.csv proton_cd_sec3.csv \
      --particle proton \
      --sector-system "CD sector"

Each input CSV should have the same normal pass-2 cross-section column names,
for example:

  normed cross sections, ep->epg, exp, Fa18 Inb, unpol

The script makes one 2x3 canvas per projection:

  xB, Q2, |t|, theta_e, theta_p, theta_gamma

The phi projection is intentionally omitted by the shared plotting engine.

Each canvas has one subplot per run period. Each subplot has two panels:

  top:    absolute integrated cross sections
  bottom: each sector divided by the arithmetic average of available sectors
          in that same run-period panel and projected bin.

Run-period panel order:

  top row:    Sp18 Inb, Sp18 Out, empty
  bottom row: Fa18 Inb, Fa18 Out, Sp19 Inb

A diagnostic chi2/ndf is printed in each run-period panel.

This wrapper performs an explicit input-file preflight check before calling the
shared plotting engine. That avoids misleading pandas EmptyDataError traces when
one sector file is missing, zero bytes, unreadable, or not parseable.

By default, this wrapper also forces all ratio panels in the shared plotting
engine to use a fixed y-axis range:

  0.0 <= ratio <= 2.0

This can be changed with:

  --ratio-ymin VALUE
  --ratio-ymax VALUE
"""

import argparse
import os
import sys
import time
from typing import List, Tuple

import pandas as pd

import compare_cross_section_by_topology as topology_engine


DEFAULT_PARTICLE = "electron"
DEFAULT_SECTOR_SYSTEM = "FD sector"

DEFAULT_RATIO_YMIN = 0.0
DEFAULT_RATIO_YMAX = 2.0


def sanitize_for_path(text: str) -> str:
    """
    Convert a short human-readable label into a safe filename/directory component.
    """

    cleaned = text.strip().lower()
    cleaned = cleaned.replace("/", "_")
    cleaned = cleaned.replace("\\", "_")
    cleaned = cleaned.replace("-", "_")
    cleaned = cleaned.replace(" ", "_")

    while "__" in cleaned:
        cleaned = cleaned.replace("__", "_")

    return cleaned.strip("_")


def default_labels_for_n_files(n_files: int) -> List[str]:
    return [f"S{i}" for i in range(1, n_files + 1)]


def default_comparison_name(particle: str, sector_system: str) -> str:
    return f"{particle} {sector_system}"


def default_output_prefix(particle: str, sector_system: str) -> str:
    return sanitize_for_path(default_comparison_name(particle, sector_system))


def default_output_dir(particle: str, sector_system: str) -> str:
    return f"output/{default_output_prefix(particle, sector_system)}_comparison"


def log(msg: str) -> None:
    """
    Use the shared topology-engine logger when available.
    """

    if hasattr(topology_engine, "log"):
        topology_engine.log(msg)
    else:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_first_bytes(path: str, nbytes: int = 512) -> bytes:
    """
    Read the first bytes of a file for diagnostics.
    """

    try:
        with open(path, "rb") as fin:
            return fin.read(nbytes)
    except Exception as exc:
        return f"<could not read first bytes: {exc}>".encode("utf-8", errors="replace")


def file_diagnostic_string(path: str) -> str:
    """
    Return a compact diagnostic string for one path.
    """

    exists = os.path.exists(path)

    if not exists:
        return (
            f"path={path!r}\n"
            f"  exists=False"
        )

    try:
        size = os.path.getsize(path)
    except Exception as exc:
        return (
            f"path={path!r}\n"
            f"  exists=True\n"
            f"  size=<could not stat: {exc}>"
        )

    first = read_first_bytes(path, 256)

    return (
        f"path={path!r}\n"
        f"  exists=True\n"
        f"  size={size} bytes\n"
        f"  first_256_bytes={first!r}"
    )


def preflight_check_one_csv(path: str, label: str, template: str) -> None:
    """
    Check that one CSV exists, is non-empty, and can be parsed by pandas.

    This intentionally reads only the first few rows. The full read is still
    performed later by compare_cross_section_by_topology.read_inputs().
    """

    log(f"Preflight input {label}: {path!r}")

    if not os.path.exists(path):
        raise FileNotFoundError(
            f"Input CSV for {label} does not exist.\n"
            f"{file_diagnostic_string(path)}"
        )

    size = os.path.getsize(path)

    if size <= 0:
        raise RuntimeError(
            f"Input CSV for {label} is zero bytes.\n"
            f"{file_diagnostic_string(path)}"
        )

    try:
        preview = pd.read_csv(path, low_memory=False, nrows=3)
    except pd.errors.EmptyDataError as exc:
        raise RuntimeError(
            f"Input CSV for {label} exists but pandas found no parseable columns.\n"
            f"{file_diagnostic_string(path)}"
        ) from exc
    except Exception as exc:
        raise RuntimeError(
            f"Input CSV for {label} exists but failed a pandas preview read.\n"
            f"{type(exc).__name__}: {exc}\n"
            f"{file_diagnostic_string(path)}"
        ) from exc

    if preview.shape[1] <= 0:
        raise RuntimeError(
            f"Input CSV for {label} preview read produced zero columns.\n"
            f"{file_diagnostic_string(path)}"
        )

    log(
        f"Preflight OK for {label}: "
        f"{preview.shape[1]} columns visible in preview, "
        f"{size / (1024.0 * 1024.0):.2f} MB"
    )

    expected_any = [
        template.format(period="Fa18 Inb"),
        template.format(period="Fa18 Out"),
        template.format(period="Sp19 Inb"),
        template.format(period="Sp18 Inb"),
        template.format(period="Sp18 Out"),
    ]

    available = set(str(c) for c in preview.columns)
    found = [c for c in expected_any if c in available]

    if not found:
        print()
        print(f"[compare_cross_section_by_sector] WARNING for {label}:")
        print("  File is parseable, but none of the expected cross-section columns")
        print("  from --xs-template were found in the preview header.")
        print()
        print("  Expected one of:")
        for c in expected_any:
            print(f"    {c}")
        print()
        print("  First 25 available columns:")
        for c in list(preview.columns[:25]):
            print(f"    {c}")
        print()
        print("  Continuing anyway; the shared read_inputs() routine will perform")
        print("  the final required-column checks.")
        print()


def preflight_check_all_csvs(paths: List[str], labels: List[str], template: str) -> None:
    """
    Check all inputs before the shared plotting engine tries to read them.
    """

    log("============================================================")
    log("Preflight CSV checks")
    log("============================================================")

    if len(paths) != len(labels):
        raise ValueError(
            f"Internal error: paths length ({len(paths)}) does not match "
            f"labels length ({len(labels)})."
        )

    for path, label in zip(paths, labels):
        preflight_check_one_csv(
            path=path,
            label=label,
            template=template,
        )

    log("All preflight CSV checks passed.")


def force_shared_engine_ratio_ylim(ratio_ymin: float, ratio_ymax: float) -> None:
    """
    Force the shared compare_cross_section_by_topology.py plotting engine to use
    a fixed ratio-panel y-axis range.

    The shared engine computes ratio_ylim through dynamic_ratio_ylim_from_values().
    This wrapper replaces that function at runtime with a fixed-range version.
    """

    if ratio_ymax <= ratio_ymin:
        raise ValueError(
            f"Invalid ratio y-axis range: ymin={ratio_ymin}, ymax={ratio_ymax}"
        )

    def fixed_ratio_ylim_from_values(_values) -> Tuple[float, float]:
        return float(ratio_ymin), float(ratio_ymax)

    if not hasattr(topology_engine, "dynamic_ratio_ylim_from_values"):
        raise RuntimeError(
            "The shared plotting engine does not define "
            "dynamic_ratio_ylim_from_values(), so this wrapper cannot force "
            "the ratio-panel y-axis range. Update compare_cross_section_by_topology.py "
            "or remove this wrapper-level patch."
        )

    topology_engine.dynamic_ratio_ylim_from_values = fixed_ratio_ylim_from_values

    log(
        "Forced shared plotting-engine ratio y-axis range: "
        f"{ratio_ymin:g} to {ratio_ymax:g}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare sector-resolved DVCS cross sections from separate CSV files. "
            "Use --particle electron/photon/proton to control labels, titles, "
            "and output defaults."
        )
    )

    parser.add_argument(
        "csv_files",
        nargs="+",
        help=(
            "Sector CSV files. Usually six files ordered S1 S2 S3 S4 S5 S6. "
            "Three files are also allowed for studies with only three sectors."
        ),
    )

    parser.add_argument(
        "--particle",
        choices=["electron", "photon", "proton"],
        default=DEFAULT_PARTICLE,
        help=(
            "Which particle's sectors are being compared. This controls the "
            "default output directory, output prefix, and figure-title text. "
            f"Default: {DEFAULT_PARTICLE}"
        ),
    )

    parser.add_argument(
        "--sector-system",
        default=DEFAULT_SECTOR_SYSTEM,
        help=(
            "Detector-sector system text used in default titles/output names. "
            "Examples: 'FD sector', 'CD sector'. "
            f"Default: '{DEFAULT_SECTOR_SYSTEM}'"
        ),
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help=(
            "Labels for the input CSVs. If omitted, labels are generated as "
            "S1, S2, ..., SN."
        ),
    )

    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Output directory. If omitted, this is derived from --particle and "
            "--sector-system, e.g. output/photon_fd_sector_comparison."
        ),
    )

    parser.add_argument(
        "--output-prefix",
        default=None,
        help=(
            "Output filename prefix. If omitted, this is derived from --particle "
            "and --sector-system, e.g. photon_fd_sector."
        ),
    )

    parser.add_argument(
        "--comparison-name",
        default=None,
        help=(
            "Text used in figure titles. If omitted, this is derived from "
            "--particle and --sector-system, e.g. 'photon FD sector'."
        ),
    )

    parser.add_argument(
        "--xs-template",
        default=topology_engine.DEFAULT_XS_TEMPLATE,
        help=(
            "Cross-section column template with {period}. "
            f"Default: {topology_engine.DEFAULT_XS_TEMPLATE!r}"
        ),
    )

    parser.add_argument(
        "--no-width-weighting",
        action="store_true",
        help="Use raw sums instead of bin-width weighted integrations.",
    )

    parser.add_argument(
        "--phi-degrees",
        action="store_true",
        help="Use phi bin widths in degrees instead of radians.",
    )

    parser.add_argument(
        "--theta-bins",
        type=int,
        default=topology_engine.DEFAULT_THETA_BINS,
        help=f"Number of derived theta bins. Default: {topology_engine.DEFAULT_THETA_BINS}",
    )

    parser.add_argument(
        "--theta-binning-period",
        default=topology_engine.DEFAULT_THETA_BINNING_PERIOD,
        help=(
            "Theta binning reference period. "
            f"Default: {topology_engine.DEFAULT_THETA_BINNING_PERIOD!r}"
        ),
    )

    parser.add_argument(
        "--include-bin-to-bin-sys",
        action="store_true",
        help="Draw outer stat+fractional-bin-to-bin-systematic bars.",
    )

    parser.add_argument(
        "--bin-to-bin-sys-frac",
        type=float,
        default=0.10,
        help="Fractional bin-to-bin systematic uncertainty. Default: 0.10",
    )

    parser.add_argument(
        "--ratio-ymin",
        type=float,
        default=DEFAULT_RATIO_YMIN,
        help=f"Fixed minimum y value for all ratio panels. Default: {DEFAULT_RATIO_YMIN}",
    )

    parser.add_argument(
        "--ratio-ymax",
        type=float,
        default=DEFAULT_RATIO_YMAX,
        help=f"Fixed maximum y value for all ratio panels. Default: {DEFAULT_RATIO_YMAX}",
    )

    parser.add_argument(
        "--skip-preflight",
        action="store_true",
        help=(
            "Skip the explicit input-file preflight checks and call the shared "
            "read_inputs() routine directly."
        ),
    )

    return parser.parse_args()


def finalize_runtime_options(args: argparse.Namespace) -> argparse.Namespace:
    n_files = len(args.csv_files)

    if n_files < 2:
        raise ValueError("At least two sector CSV files are required.")

    if args.labels is None:
        args.labels = default_labels_for_n_files(n_files)

    if len(args.labels) != n_files:
        raise ValueError(
            f"Number of --labels entries ({len(args.labels)}) must match "
            f"number of CSV files ({n_files})."
        )

    if args.output_prefix is None:
        args.output_prefix = default_output_prefix(
            particle=args.particle,
            sector_system=args.sector_system,
        )

    if args.output_dir is None:
        args.output_dir = default_output_dir(
            particle=args.particle,
            sector_system=args.sector_system,
        )

    if args.comparison_name is None:
        args.comparison_name = default_comparison_name(
            particle=args.particle,
            sector_system=args.sector_system,
        )

    if "{period}" not in args.xs_template:
        raise ValueError("--xs-template must contain {period}")

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive")

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative")

    if args.ratio_ymax <= args.ratio_ymin:
        raise ValueError(
            f"--ratio-ymax must be greater than --ratio-ymin. "
            f"Got ymin={args.ratio_ymin}, ymax={args.ratio_ymax}."
        )

    return args


def print_runtime_summary(args: argparse.Namespace) -> None:
    log("============================================================")
    log("compare_cross_section_by_sector.py")
    log("============================================================")
    log(f"Particle: {args.particle}")
    log(f"Sector system: {args.sector_system}")
    log(f"Comparison name: {args.comparison_name}")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {args.labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"Output prefix: {args.output_prefix}")
    log(f"XS template: {args.xs_template}")
    log(f"Theta bins: {args.theta_bins}")
    log(f"Theta binning period: {args.theta_binning_period}")
    log(f"Width weighting: {'disabled' if args.no_width_weighting else 'enabled'}")
    log(f"Phi bin widths: {'degrees' if args.phi_degrees else 'radians'}")
    log(f"Bin-to-bin sys bars: {'enabled' if args.include_bin_to_bin_sys else 'disabled'}")
    log(f"Bin-to-bin sys frac: {args.bin_to_bin_sys_frac}")
    log(f"Ratio y-axis range: {args.ratio_ymin:g} to {args.ratio_ymax:g}")
    log(f"Preflight checks: {'disabled' if args.skip_preflight else 'enabled'}")
    log("phi dependence plots: disabled")


def read_inputs_with_diagnostics(args: argparse.Namespace):
    """
    Call the shared read_inputs(), but make pandas EmptyDataError messages less
    ambiguous if something still fails after preflight.
    """

    try:
        return topology_engine.read_inputs(
            paths=args.csv_files,
            labels=args.labels,
            template=args.xs_template,
        )
    except pd.errors.EmptyDataError as exc:
        print()
        print("[compare_cross_section_by_sector] ERROR:")
        print("  The shared read_inputs() routine hit pandas EmptyDataError.")
        print("  Reprinting file diagnostics for every input:")
        print()

        for path, label in zip(args.csv_files, args.labels):
            print(f"--- {label} ---")
            print(file_diagnostic_string(path))
            print()

        raise RuntimeError(
            "Shared read_inputs() failed with EmptyDataError. "
            "Check the diagnostics above for the file with zero size, no header, "
            "or unreadable contents."
        ) from exc


def main() -> None:
    t0 = time.perf_counter()

    try:
        args = parse_args()
        args = finalize_runtime_options(args)

        print_runtime_summary(args)

        force_shared_engine_ratio_ylim(
            ratio_ymin=args.ratio_ymin,
            ratio_ymax=args.ratio_ymax,
        )

        if not args.skip_preflight:
            preflight_check_all_csvs(
                paths=args.csv_files,
                labels=args.labels,
                template=args.xs_template,
            )

        dfs = read_inputs_with_diagnostics(args)

        dfs = topology_engine.prepare_dataframes(
            dfs=dfs,
            theta_binning_period=args.theta_binning_period,
            theta_bins=args.theta_bins,
            phi_degrees=args.phi_degrees,
        )

        topology_engine.make_all_projection_canvases(
            dfs=dfs,
            output_dir=args.output_dir,
            output_prefix=args.output_prefix,
            template=args.xs_template,
            no_width_weighting=args.no_width_weighting,
            include_bin_to_bin_sys=args.include_bin_to_bin_sys,
            frac=args.bin_to_bin_sys_frac,
            comparison_name=args.comparison_name,
        )

        log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
        log("Done.")

    except Exception as exc:
        print()
        print("[compare_cross_section_by_sector] FATAL:")
        print(f"  {type(exc).__name__}: {exc}")
        print()
        sys.exit(1)


if __name__ == "__main__":
    main()