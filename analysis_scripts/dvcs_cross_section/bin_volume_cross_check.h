#ifndef BIN_VOLUME_CROSS_CHECK_H
#define BIN_VOLUME_CROSS_CHECK_H

#include <string>

/**
 * Compare bin volumes between Hayward (pass-2) and Lee (pass-1).
 *
 * Reads directly from CSVs:
 *   - Lee CSV column:      "bin_volume"
 *   - Hayward CSV column:  "bin_volume, 10.6 GeV"
 *
 * The Hayward column is stored as a three-tuple "(value, stat_err, syst_err)".
 * We use the central value and the statistical uncertainty (second entry).
 *
 * For each valid bin, we compare bin_volume(Lee) and bin_volume(Hayward),
 * organized by xB, Q2, and -t. For each (xB, Q2, -t) cell, we take all rows
 * matching those ranges and use their provided phiavg values as x-coordinates.
 * We do NOT rebin phi.
 *
 * Produces comparison plots (overlay and ratio) organized by xB bins with
 * panels for each (Q2, t) combination.
 *
 * Output filenames (per xB index ix):
 *   bin_volume_counts_xB_<ix>.png
 *   bin_volume_ratio_xB_<ix>.png
 *
 * @param lee_csv_path      Path to Lee's pass-1 CSV (e.g. imports/all_bin_v3.csv)
 * @param hayward_csv_path  Path to Hayward's pass-2 CSV (e.g. output/csvs/dvcs_pass2_analysis.csv)
 * @param output_base_dir   Directory for output plots (e.g. output/cross_check/lee/bin_volume)
 */
void plot_bin_volume_cross_checks(const std::string& lee_csv_path,
                                  const std::string& hayward_csv_path,
                                  const std::string& output_base_dir);

#endif // BIN_VOLUME_CROSS_CHECK_H