#ifndef PLOT_YIELD_COMPARISON_H
#define PLOT_YIELD_COMPARISON_H

#include <string>
#include <vector>
#include "all_bin_data.h"  // Include the struct for bin data

// Function declaration for plotting yield comparison
void plot_yield_comparison(const std::string& output_dir, int xB_bin, const std::vector<AllBinData>& all_bin_data);

#endif  // PLOT_YIELD_COMPARISON_H