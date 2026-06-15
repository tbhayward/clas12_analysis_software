#!/usr/bin/env python3
"""
compare_cross_sections_by_proton_cd_sector.py

Compare integrated DVCS cross sections from separate proton-CD-sector CSV files.

Expected usage:

  python compare_cross_sections_by_proton_cd_sector.py S1.csv S2.csv S3.csv

Each input CSV should have the same normal pass-2 cross-section column names,
for example:

  normed cross sections, ep->epg, exp, Fa18 Inb, unpol

Outputs one 2x3 canvas per projection:

  xB, Q2, |t|, phi, theta_e, theta_p, theta_gamma
"""

from compare_cross_section_by_topology import (
    DEFAULT_THETA_BINS,
    DEFAULT_THETA_BINNING_PERIOD,
    DEFAULT_XS_TEMPLATE,
    PROJECTION_ORDER,
    Timer,
    log,
    read_inputs,
    prepare_dataframes,
    make_projection_canvas,
)

import argparse
import time


DEFAULT_OUTPUT_DIR = "output/proton_cd_sector_comparison"
DEFAULT_OUTPUT_PREFIX = "proton_cd_sector"
DEFAULT_LABELS = ["CD S1", "CD S2", "CD S3"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare proton-CD-sector DVCS cross sections from separate CSVs.")
    p.add_argument("csv_files", nargs=3, help="CSV files in order: CD S1 CD S2 CD S3")
    p.add_argument("--labels", nargs=3, default=DEFAULT_LABELS, help="Labels for the three CSVs.")
    p.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    p.add_argument("--xs-template", default=DEFAULT_XS_TEMPLATE, help="Cross-section column template with {period}.")
    p.add_argument("--no-width-weighting", action="store_true")
    p.add_argument("--phi-degrees", action="store_true")
    p.add_argument("--theta-bins", type=int, default=DEFAULT_THETA_BINS)
    p.add_argument("--theta-binning-period", default=DEFAULT_THETA_BINNING_PERIOD)
    p.add_argument("--include-bin-to-bin-sys", action="store_true")
    p.add_argument("--bin-to-bin-sys-frac", type=float, default=0.10)
    return p.parse_args()


def main() -> None:
    t0 = time.perf_counter()
    args = parse_args()

    log("============================================================")
    log("compare_cross_sections_by_proton_cd_sector.py")
    log("============================================================")
    log(f"CSV files: {args.csv_files}")
    log(f"Labels: {args.labels}")
    log(f"Output dir: {args.output_dir}")
    log(f"XS template: {args.xs_template}")

    dfs = read_inputs(args.csv_files, args.labels, args.xs_template)
    dfs = prepare_dataframes(dfs, args.theta_binning_period, args.theta_bins, args.phi_degrees)

    for projection in PROJECTION_ORDER:
        with Timer(f"plotting {projection}"):
            make_projection_canvas(
                dfs=dfs,
                projection=projection,
                output_dir=args.output_dir,
                output_prefix=args.output_prefix,
                template=args.xs_template,
                no_width_weighting=args.no_width_weighting,
                include_bin_to_bin_sys=args.include_bin_to_bin_sys,
                frac=args.bin_to_bin_sys_frac,
            )

    log(f"TOTAL RUNTIME: {time.perf_counter() - t0:.3f} s")
    log("Done.")


if __name__ == "__main__":
    main()
# endif