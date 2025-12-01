#ifndef UNFOLDED_YIELDS_CROSS_CHECK_H
#define UNFOLDED_YIELDS_CROSS_CHECK_H

#include <string>

// Compare Lee/pass-1 acceptance-corrected yields to
// Hayward/pass-2 unfolded (acceptance-corrected) yields and
// make a ROOT plot.
//
// - lee_csv     : Lee/pass-1 CSV (e.g. imports/all_bin_v3.csv)
// - hayward_csv : Hayward/pass-2 CSV (dvcs_pass2_analysis.csv)
// - out_dir     : output directory for plots (e.g. output/cross_check/lee/unfolding)
void plot_unfolded_yields_cross_checks(const std::string &lee_csv,
                                       const std::string &hayward_csv,
                                       const std::string &out_dir);

#endif // UNFOLDED_YIELDS_CROSS_CHECK_H