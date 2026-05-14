#ifndef CROSS_SECTION_CROSS_CHECK_H
#define CROSS_SECTION_CROSS_CHECK_H

#include <string>

/**
 * Compare experimental cross sections between Hayward (pass-2) and Lee (pass-1).
 *
 * Reads directly from CSVs:
 *   - Lee CSV column:     "cross sections, ep->epg, exp"
 *   - Hayward CSV column: "cross sections, ep->epg, exp, Fa18, unpol"
 *
 * The Hayward column is stored as a three-tuple "(value, stat_err, syst_err)".
 * We use the central value and the statistical uncertainty (second entry).
 *
 * For each valid bin, we compare cross_section(Lee) and cross_section(Hayward),
 * organized by xB, Q^{2}, and -t. For each (xB, Q^{2}, -t) cell, we take all
 * rows matching those ranges and use their provided phiavg values as the
 * x-coordinates. We do NOT rebin phi.
 *
 * Produces comparison plots (overlay and ratio) organized by xB bins with
 * panels for each (Q^{2}, t) combination. The counts canvases overlay
 * the BH theory prediction as a green dashed curve.
 *
 * Output filenames (per xB index ix):
 *   cross_section_counts_xB_<ix>.png
 *   cross_section_ratio_xB_<ix>.png
 *
 * @param lee_csv_path      Path to Lee's pass-1 CSV (e.g. imports/all_bin_v3.csv)
 * @param hayward_csv_path  Path to Hayward's pass-2 CSV (e.g. output/csvs/dvcs_pass2_analysis.csv)
 * @param output_base_dir   Directory for output plots (e.g. output/cross_check/lee/cross_sections)
 * @param theory_json_root  Root directory containing theory JSONs, e.g. output/jsons/cross_sections
 */
void plot_cross_section_cross_checks(const std::string& lee_csv_path,
                                     const std::string& hayward_csv_path,
                                     const std::string& output_base_dir,
                                     const std::string& theory_json_root = "output/jsons/cross_sections");

#endif // CROSS_SECTION_CROSS_CHECK_H