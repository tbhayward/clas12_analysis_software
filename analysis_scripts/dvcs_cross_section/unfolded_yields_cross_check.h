#ifndef UNFOLDED_YIELDS_CROSS_CHECK_H
#define UNFOLDED_YIELDS_CROSS_CHECK_H

#include <string>

// Compare unfolded acceptance-corrected yields between Lee/pass-1 CSV
// and Hayward/pass-2 CSV (Fa18, unpolarized). Produces a single 1D plot.
void plot_unfolded_yields_cross_checks(const std::string &lee_csv,
                                       const std::string &hayward_csv,
                                       const std::string &out_dir);

#endif // UNFOLDED_YIELDS_CROSS_CHECK_H