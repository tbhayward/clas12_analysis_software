#!/usr/bin/env python3
"""
compare_cross_section_by_proton_cd_sector.py

Compare integrated DVCS cross sections from separate proton-CD-sector CSV files.

Expected usage:

  python compare_cross_section_by_proton_cd_sector.py \
    cd_sec1.csv cd_sec2.csv cd_sec3.csv

Each input CSV should have the same normal pass-2 cross-section column names,
for example:

  normed cross sections, ep->epg, exp, Fa18 Inb, unpol

The script makes one 2x3 canvas per projection:

  xB, Q2, |t|, phi, theta_e, theta_p, theta_gamma

Each canvas has one subplot per run period. Each subplot has two panels:

  top:    absolute integrated cross sections
  bottom: each sector divided by the arithmetic average of available sectors
          in that same run-period panel and projected bin.

A diagnostic chi2/ndf is printed in each run-period panel.

Outputs by default:

  output/proton_cd_sector_comparison/proton_cd_sector_xB_comparison.png
  output/proton_cd_sector_comparison/proton_cd_sector_Q2_comparison.png
  output/proton_cd_sector_comparison/proton_cd_sector_t_comparison.png
  output/proton_cd_sector_comparison/proton_cd_sector_phi_comparison.png
  output/proton_cd_sector_comparison/proton_cd_sector_e_theta_comparison.png
  output/proton_cd_sector_comparison/proton_cd_sector_p_theta_comparison.png
  output/proton_cd_sector_comparison/proton_cd_sector_g_theta_comparison.png
"""

import argparse
import time

from compare_cross_section_by_topology import (
    DEFAULT_THETA_BINS,
    DEFAULT_THETA_BINNING_PERIOD,
    DEFAULT_XS_TEMPLATE,
    log,
    read_inputs,
    prepare_dataframes,
    make_all_projection_canvases,
)


DEFAULT_OUTPUT_DIR = "output/proton_cd_sector_comparison"
DEFAULT_OUTPUT_PREFIX = "proton_cd_sector"
DEFAULT_LABELS = ["CD S1", "CD S2", "CD S3"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare proton-CD-sector DVCS cross sections from separate CSVs."
    )

    parser.add_argument(
        "csv_files",
        nargs=3,
        help="CSV files in order: CD S1 CD S2 CD S3",
    )

    parser.add_argument(
        "--labels",
        nargs=3,
        default=DEFAULT_LABELS,
        help="Labels for the three CSVs. Default: 'CD S1' 'CD S2' 'CD S3'",
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
        "--xs-template",
        default=DEFAULT_XS_TEMPLATE,
        help=f"Cross-section column template with {{period}}. Default: {DEFAULT_XS_TEMPLATE!r}",
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
        help=f"Theta binning reference period. Default: {DEFAULT_THETA_BINNING_PERIOD!r}",
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
        "--comparison-name",
        default="proton CD sector",
        help="Text used in figure titles. Default: proton CD sector",
    )

    return parser.parse_args()


def main() -> None:
    t0 = time.perf_counter()
    args = parse_args()

    log("============================================================")
    log("compare_cross_section_by_proton_cd_sector.py")
    log("============================================================")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {args.labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"XS template: {args.xs_template}")
    log(f"Theta bins: {args.theta_bins}")
    log(f"Theta binning period: {args.theta_binning_period}")

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