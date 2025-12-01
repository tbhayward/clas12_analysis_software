#ifndef BIN_CENTERING_CROSS_CHECK_H
#define BIN_CENTERING_CROSS_CHECK_H

#include <string>

/**
 * Compare bin-centering corrections between Hayward (pass-2) and Lee (pass-1).
 *
 * Reads directly from CSVs:
 *   - Lee CSV column:       "Fbin"
 *   - Hayward CSV column:   "bin_volume, 10.6 GeV"
 *
 * This is organized exactly like the radiative correction cross-check:
 * we make overlay and ratio plots as a function of phiavg for each (Q2, t)
 * panel within each xB bin.
 *
 * @param lee_csv_path      Path to Lee's pass-1 CSV (e.g. imports/all_bin_v3.csv)
 * @param hayward_csv_path  Path to Hayward's pass-2 CSV (e.g. output/csvs/dvcs_pass2_analysis.csv)
 * @param output_base_dir   Directory for output plots
 *                          (e.g. output/cross_check/lee/bin_centering)
 */
void plot_bin_centering_cross_checks(const std::string& lee_csv_path,
                                     const std::string& hayward_csv_path,
                                     const std::string& output_base_dir);

#endif // BIN_CENTERING_CROSS_CHECK_H