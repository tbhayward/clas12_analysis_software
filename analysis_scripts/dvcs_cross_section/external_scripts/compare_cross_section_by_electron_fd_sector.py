#!/usr/bin/env python3
"""
compare_cross_section_by_sector.py

Compare integrated DVCS cross sections from separate sector-resolved CSV files.

This is a generic sector-comparison wrapper around the shared plotting engine in:

  compare_cross_section_by_topology.py

It can be used for electron, photon, or proton sector studies.

Examples:

  Electron FD sectors:

    python compare_cross_section_by_sector.py \
      elec_sec1.csv elec_sec2.csv elec_sec3.csv \
      elec_sec4.csv elec_sec5.csv elec_sec6.csv \
      --particle electron

  Photon FD sectors:

    python compare_cross_section_by_sector.py \
      photon_sec1.csv photon_sec2.csv photon_sec3.csv \
      photon_sec4.csv photon_sec5.csv photon_sec6.csv \
      --particle photon

  Proton sectors:

    python compare_cross_section_by_sector.py \
      proton_sec1.csv proton_sec2.csv proton_sec3.csv \
      proton_sec4.csv proton_sec5.csv proton_sec6.csv \
      --particle proton

  Proton CD sectors, if only three sector files exist:

    python compare_cross_section_by_sector.py \
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

Defaults depend on --particle:

  --particle electron:
    output dir    = output/electron_fd_sector_comparison
    output prefix = electron_fd_sector
    title text    = electron FD sector

  --particle photon:
    output dir    = output/photon_fd_sector_comparison
    output prefix = photon_fd_sector
    title text    = photon FD sector

  --particle proton:
    output dir    = output/proton_fd_sector_comparison
    output prefix = proton_fd_sector
    title text    = proton FD sector

Use --sector-system to override the detector-sector text, for example:

  --sector-system "CD sector"

Use --comparison-name to override the full figure-title comparison text.
"""

import argparse
import time
from typing import List

from compare_cross_section_by_topology import (
    DEFAULT_THETA_BINS,
    DEFAULT_THETA_BINNING_PERIOD,
    DEFAULT_XS_TEMPLATE,
    log,
    read_inputs,
    prepare_dataframes,
    make_all_projection_canvases,
)


DEFAULT_PARTICLE = "electron"
DEFAULT_SECTOR_SYSTEM = "FD sector"


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
    # endwhile

    return cleaned.strip("_")


def default_labels_for_n_files(n_files: int) -> List[str]:
    return [f"S{i}" for i in range(1, n_files + 1)]


def default_comparison_name(particle: str, sector_system: str) -> str:
    return f"{particle} {sector_system}"


def default_output_prefix(particle: str, sector_system: str) -> str:
    return sanitize_for_path(default_comparison_name(particle, sector_system))


def default_output_dir(particle: str, sector_system: str) -> str:
    return f"output/{default_output_prefix(particle, sector_system)}_comparison"


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
        default=DEFAULT_XS_TEMPLATE,
        help=(
            "Cross-section column template with {period}. "
            f"Default: {DEFAULT_XS_TEMPLATE!r}"
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
        default=DEFAULT_THETA_BINS,
        help=f"Number of derived theta bins. Default: {DEFAULT_THETA_BINS}",
    )

    parser.add_argument(
        "--theta-binning-period",
        default=DEFAULT_THETA_BINNING_PERIOD,
        help=(
            "Theta binning reference period. "
            f"Default: {DEFAULT_THETA_BINNING_PERIOD!r}"
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

    return parser.parse_args()


def finalize_runtime_options(args: argparse.Namespace) -> argparse.Namespace:
    n_files = len(args.csv_files)

    if n_files < 2:
        raise ValueError("At least two sector CSV files are required.")
    # endif

    if args.labels is None:
        args.labels = default_labels_for_n_files(n_files)
    # endif

    if len(args.labels) != n_files:
        raise ValueError(
            f"Number of --labels entries ({len(args.labels)}) must match "
            f"number of CSV files ({n_files})."
        )
    # endif

    if args.output_prefix is None:
        args.output_prefix = default_output_prefix(
            particle=args.particle,
            sector_system=args.sector_system,
        )
    # endif

    if args.output_dir is None:
        args.output_dir = default_output_dir(
            particle=args.particle,
            sector_system=args.sector_system,
        )
    # endif

    if args.comparison_name is None:
        args.comparison_name = default_comparison_name(
            particle=args.particle,
            sector_system=args.sector_system,
        )
    # endif

    return args


def main() -> None:
    t0 = time.perf_counter()

    args = parse_args()
    args = finalize_runtime_options(args)

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
    log("phi dependence plots: disabled")

    if "{period}" not in args.xs_template:
        raise ValueError("--xs-template must contain {period}")
    # endif

    if args.theta_bins <= 0:
        raise ValueError("--theta-bins must be positive")
    # endif

    if args.bin_to_bin_sys_frac < 0.0:
        raise ValueError("--bin-to-bin-sys-frac must be non-negative")
    # endif

    dfs = read_inputs(
        paths=args.csv_files,
        labels=args.labels,
        template=args.xs_template,
    )

    dfs = prepare_dataframes(
        dfs=dfs,
        theta_binning_period=args.theta_binning_period,
        theta_bins=args.theta_bins,
        phi_degrees=args.phi_degrees,
    )

    make_all_projection_canvases(
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


if __name__ == "__main__":
    main()
# endif