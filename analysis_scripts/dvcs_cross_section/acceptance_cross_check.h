#ifndef ACCEPTANCE_CROSS_CHECK_H
#define ACCEPTANCE_CROSS_CHECK_H

#include <string>

/**
 * Compare acceptances between Hayward (pass-2) and Lee (pass-1).
 *
 * Reads directly from CSVs:
 *   - Lee column:     "acceptance, ep->epg, sim(KM15)"
 *   - Hayward columns:
 *        "acceptance, Fa18 Inb"
 *        "acceptance, Fa18 Out"
 *
 * The Hayward columns are stored as three-tuples "(value, stat_err, syst_err)".
 * We use the central values and the statistical uncertainties (second entries).
 *
 * For each valid bin, we compare acceptance(Lee) with acceptance(Hayward Inb/Out),
 * organized by xB, Q^{2}, and -t. For each (xB, Q^{2}, -t) cell, we take all
 * rows matching those ranges and use their provided phiavg values as the
 * x-coordinates. We do NOT rebin phi.
 *
 * For each xB bin ix, we produce:
 *   - acceptance_counts_xB_<ix>.png
 *       Overlays: Lee, Fa18 Inb, Fa18 Out
 *   - acceptance_ratio_xB_<ix>.png
 *       Ratios: Fa18 Inb / Lee and Fa18 Out / Lee, plus y = 1 line
 *
 * @param lee_csv_path      Path to Lee's pass-1 CSV (e.g. imports/all_bin_v3.csv)
 * @param hayward_csv_path  Path to Hayward's pass-2 CSV (e.g. output/csvs/dvcs_pass2_analysis.csv)
 * @param output_base_dir   Directory for output plots (e.g. output/cross_check/lee/acceptance)
 */
void plot_acceptance_cross_checks(const std::string& lee_csv_path,
                                  const std::string& hayward_csv_path,
                                  const std::string& output_base_dir);

#endif // ACCEPTANCE_CROSS_CHECK_H