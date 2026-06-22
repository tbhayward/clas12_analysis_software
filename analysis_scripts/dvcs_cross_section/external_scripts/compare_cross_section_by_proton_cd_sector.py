#!/usr/bin/env python3
"""
compare_cross_section_by_proton_cd_sector.py

Compare proton-CD-sector DVCS pass-2 quantities from separate CSV files.

Expected usage:

  python3 compare_cross_section_by_proton_cd_sector.py \
    cd_sec1.csv cd_sec2.csv cd_sec3.csv

The original cross-section comparison functionality is preserved. In addition,
this wrapper now asks the shared plotting engine in compare_cross_section_by_topology.py
to produce step-by-step diagnostic plots for raw yields, current factors,
normalized raw yields, generated/reconstructed MC yields, pi0 contamination,
signal yields, acceptance, and acceptance-corrected yields.

All ratio panels are standardized to 0 <= ratio <= 2 by default.
"""

import argparse
import os
import sys
import time
from typing import List, Tuple

import pandas as pd

import compare_cross_section_by_topology as topology_engine


DEFAULT_OUTPUT_DIR = "output/proton_cd_sector_comparison"
DEFAULT_OUTPUT_PREFIX = "proton_cd_sector"
DEFAULT_LABELS = ["CD S1", "CD S2", "CD S3"]
DEFAULT_COMPARISON_NAME = "proton CD sector"

DEFAULT_RATIO_YMIN = 0.0
DEFAULT_RATIO_YMAX = 2.0


def log(msg: str) -> None:
    if hasattr(topology_engine, "log"):
        topology_engine.log(msg)
    else:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_first_bytes(path: str, nbytes: int = 512) -> bytes:
    try:
        with open(path, "rb") as fin:
            return fin.read(nbytes)
    except Exception as exc:
        return f"<could not read first bytes: {exc}>".encode("utf-8", errors="replace")


def file_diagnostic_string(path: str) -> str:
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
        print(f"[compare_cross_section_by_proton_cd_sector] WARNING for {label}:")
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
    if ratio_ymax <= ratio_ymin:
        raise ValueError(
            f"Invalid ratio y-axis range: ymin={ratio_ymin}, ymax={ratio_ymax}"
        )

    def fixed_ratio_ylim_from_values(_values) -> Tuple[float, float]:
        return float(ratio_ymin), float(ratio_ymax)

    topology_engine.dynamic_ratio_ylim_from_values = fixed_ratio_ylim_from_values

    if hasattr(topology_engine, "DEFAULT_RATIO_YMIN"):
        topology_engine.DEFAULT_RATIO_YMIN = float(ratio_ymin)

    if hasattr(topology_engine, "DEFAULT_RATIO_YMAX"):
        topology_engine.DEFAULT_RATIO_YMAX = float(ratio_ymax)

    log(
        "Forced shared plotting-engine ratio y-axis range: "
        f"{ratio_ymin:g} to {ratio_ymax:g}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare proton-CD-sector DVCS pass-2 quantities from separate CSV files."
    )

    parser.add_argument(
        "csv_files",
        nargs="+",
        help="Proton-CD-sector CSV files. Usually three files ordered CD S1, CD S2, CD S3.",
    )

    parser.add_argument(
        "--labels",
        nargs="+",
        default=None,
        help=(
            "Labels for the input CSVs. If omitted and three files are given, "
            "the default is CD S1, CD S2, CD S3. Otherwise labels are generated."
        ),
    )

    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help=f"Output directory. Default: {DEFAULT_OUTPUT_DIR}",
    )

    parser.add_argument(
        "--output-prefix",
        default=DEFAULT_OUTPUT_PREFIX,
        help=f"Output filename prefix. Default: {DEFAULT_OUTPUT_PREFIX}",
    )

    parser.add_argument(
        "--comparison-name",
        default=DEFAULT_COMPARISON_NAME,
        help=f"Text used in figure titles. Default: {DEFAULT_COMPARISON_NAME!r}",
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
        help=(
            "Use raw sums instead of bin-width weighted integrations for "
            "width-weighted metrics. Yield diagnostics are always raw sums."
        ),
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
        help="Skip the explicit input-file preflight checks.",
    )

    parser.add_argument(
        "--no-step-diagnostics",
        action="store_true",
        help="Disable expanded step-by-step diagnostic canvases.",
    )

    parser.add_argument(
        "--diagnostic-metrics",
        nargs="+",
        default=None,
        help=(
            "Optional subset of diagnostic metric tags to plot. "
            "Use output/.../step_diagnostics/*_step_diagnostics_index.txt for the tag list."
        ),
    )

    return parser.parse_args()


def finalize_runtime_options(args: argparse.Namespace) -> argparse.Namespace:
    n_files = len(args.csv_files)

    if n_files < 2:
        raise ValueError("At least two sector CSV files are required.")

    if args.labels is None:
        if n_files == len(DEFAULT_LABELS):
            args.labels = list(DEFAULT_LABELS)
        else:
            args.labels = [f"CD S{i}" for i in range(1, n_files + 1)]

    if len(args.labels) != n_files:
        raise ValueError(
            f"Number of --labels entries ({len(args.labels)}) must match "
            f"number of CSV files ({n_files})."
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
    log("compare_cross_section_by_proton_cd_sector.py")
    log("============================================================")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {args.labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"Output prefix: {args.output_prefix}")
    log(f"Comparison name: {args.comparison_name}")
    log(f"XS template: {args.xs_template}")
    log(f"Theta bins: {args.theta_bins}")
    log(f"Theta binning period: {args.theta_binning_period}")
    log(f"Width weighting: {'disabled' if args.no_width_weighting else 'enabled'}")
    log(f"Phi bin widths: {'degrees' if args.phi_degrees else 'radians'}")
    log(f"Bin-to-bin sys bars: {'enabled' if args.include_bin_to_bin_sys else 'disabled'}")
    log(f"Bin-to-bin sys frac: {args.bin_to_bin_sys_frac}")
    log(f"Ratio y-axis range: {args.ratio_ymin:g} to {args.ratio_ymax:g}")
    log(f"Preflight checks: {'disabled' if args.skip_preflight else 'enabled'}")
    log(f"Step diagnostics: {'disabled' if args.no_step_diagnostics else 'enabled'}")
    if args.diagnostic_metrics:
        log(f"Diagnostic metric subset: {args.diagnostic_metrics}")
    log("phi dependence plots: disabled")


def read_inputs_with_diagnostics(args: argparse.Namespace):
    try:
        return topology_engine.read_inputs(
            paths=args.csv_files,
            labels=args.labels,
            template=args.xs_template,
        )
    except pd.errors.EmptyDataError as exc:
        print()
        print("[compare_cross_section_by_proton_cd_sector] ERROR:")
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

        # Original/final cross-section comparison canvases, same filenames as before.
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

        # Expanded upstream/downstream pass-2 diagnostic canvases.
        if not args.no_step_diagnostics:
            topology_engine.make_all_step_diagnostic_canvases(
                dfs=dfs,
                output_dir=args.output_dir,
                output_prefix=args.output_prefix,
                no_width_weighting=args.no_width_weighting,
                include_bin_to_bin_sys=args.include_bin_to_bin_sys,
                frac=args.bin_to_bin_sys_frac,
                comparison_name=args.comparison_name,
                metric_tags=args.diagnostic_metrics,
            )

        log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
        log("Done.")

    except Exception as exc:
        print()
        print("[compare_cross_section_by_proton_cd_sector] FATAL:")
        print(f"  {type(exc).__name__}: {exc}")
        print()
        sys.exit(1)


if __name__ == "__main__":
    main()
