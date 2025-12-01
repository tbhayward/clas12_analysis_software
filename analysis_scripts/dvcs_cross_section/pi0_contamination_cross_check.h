#ifndef PI0_CONTAMINATION_CROSS_CHECK_H
#define PI0_CONTAMINATION_CROSS_CHECK_H

#include <string>

/**
 * Compare pi0 contamination ratios between Hayward (pass-2) and Lee (pass-1).
 *
 * Reads directly from CSVs:
 *   - Lee CSV columns: "contamination ratio, inbending", "contamination ratio, outbending"
 *   - Hayward CSV columns: "contamination ratio, Fa18 Inb", "contamination ratio, Fa18 Out"
 *     (these are three-tuples "(value, stat_err, syst_err)" - we extract value and stat_err)
 *
 * Produces comparison plots (counts overlay and ratio) for Fa18 Inb and Fa18 Out,
 * organized by xB bins with panels for each (Q2, t) combination.
 *
 * @param lee_csv_path      Path to Lee's pass-1 CSV (e.g. imports/all_bin_v3.csv)
 * @param hayward_csv_path  Path to Hayward's pass-2 CSV (e.g. output/csvs/dvcs_pass2_analysis.csv)
 * @param output_base_dir   Directory for output plots (e.g. output/cross_check/lee/pi0_contamination)
 */
void plot_pi0_contam_cross_checks(const std::string& lee_csv_path,
                                  const std::string& hayward_csv_path,
                                  const std::string& output_base_dir);

#endif // PI0_CONTAMINATION_CROSS_CHECK_H