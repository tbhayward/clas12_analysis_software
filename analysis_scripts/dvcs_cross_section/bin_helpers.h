#ifndef BIN_HELPERS_H
#define BIN_HELPERS_H

#include <string>
#include <vector>
#include "bin_boundaries.h"  // Include your bin boundary definition

// Function to clean up the bin label
std::string clean_bin_label(const std::string& label);

// Function to precompute relevant bins based on xB_bin
std::vector<int> precompute_relevant_bins(int xB_bin, const std::vector<BinBoundary>& bin_boundaries);

#endif // BIN_HELPERS_H