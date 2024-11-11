#ifndef BIN_HELPERS_H
#define BIN_HELPERS_H

#include <vector>
#include "bin_boundaries.h"  // Include your bin boundary definition

// Function to precompute relevant bins based on xB_bin
std::vector<int> precompute_relevant_bins(int xB_bin, const std::vector<BinBoundary>& bin_boundaries);

// Find the next perfect square greater than or equal to the given number
int next_perfect_square(int n);

#endif // BIN_HELPERS_H