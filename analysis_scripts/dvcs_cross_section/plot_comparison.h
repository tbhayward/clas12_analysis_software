#ifndef PLOT_COMPARISON_H
#define PLOT_COMPARISON_H

#include <string>
#include <vector>
#include <map>
#include "unfolding_data.h"
#include "bin_boundaries.h"

void plot_comparison(
    const std::string& output_dir,
    const std::vector<BinBoundary>& bin_boundaries,
    const std::string& previous_data_csv,
    const std::map<std::string, std::vector<UnfoldingData>>& all_unfolding_data
);

#endif // PLOT_COMPARISON_H