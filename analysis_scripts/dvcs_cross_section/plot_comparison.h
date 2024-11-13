#ifndef PLOT_COMPARISON_H
#define PLOT_COMPARISON_H

#include <string>
#include <vector>
#include "bin_boundaries.h"

// Adjusted function signature to match plot_comparison.cpp
void plot_comparison(
    const std::string& output_dir,
    const std::vector<BinBoundary>& bin_boundaries,
    const std::string& csv_fa18_inb,
    const std::string& csv_fa18_out
);

#endif // PLOT_COMPARISON_H