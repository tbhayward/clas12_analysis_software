#ifndef RAD_CORRECTION_CROSS_CHECK_H
#define RAD_CORRECTION_CROSS_CHECK_H

#include <string>

/**
 * Compare radiative correction factors (Frad) between Hayward (pass-2) and Lee (pass-1).
 *
 * Reads directly from CSVs:
 *   - Lee CSV column: "Frad"
 *   - Hayward CSV column: "Frad, 10.6 GeV"
 *
 * Since radiative corrections are beam-energy dependent (not polarity dependent),
 * there is no inbending/outbending split - just one comparison for 10.6 GeV.
 *
 * Produces comparison plots (overlay and ratio) organized by xB bins with
 * panels for each (Q2, t) combination.
 *
 * @param lee_csv_path      Path to Lee's pass-1 CSV (e.g. imports/all_bin_v3.csv)
 * @param hayward_csv_path  Path to Hayward's pass-2 CSV (e.g. output/csvs/dvcs_pass2_analysis.csv)
 * @param output_base_dir   Directory for output plots (e.g. output/cross_check/lee/rad_corrections)
 */
void plot_rad_correction_cross_checks(const std::string& lee_csv_path,
                                      const std::string& hayward_csv_path,
                                      const std::string& output_base_dir);

#endif // RAD_CORRECTION_CROSS_CHECK_H