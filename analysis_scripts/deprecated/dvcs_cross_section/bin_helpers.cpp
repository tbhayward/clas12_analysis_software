#include "bin_helpers.h"
#include <iostream>
#include <set>
#include <vector>
#include <cmath>  
#include <algorithm>

// Function to precompute relevant bins based on xB_bin
std::vector<int> precompute_relevant_bins(int xB_bin, const std::vector<BinBoundary>& bin_boundaries) {
    std::vector<int> relevant_bins;
    // std::cout << "Entering precompute_relevant_bins with xB_bin = " << xB_bin << std::endl;

    // Extract unique xB ranges from bin_boundaries
    std::set<std::pair<double, double>> xB_ranges_set;
    for (const auto& bin : bin_boundaries) {
        xB_ranges_set.insert(std::make_pair(bin.xB_low, bin.xB_high));
    }

    // Convert set to vector for indexing
    std::vector<std::pair<double, double>> xB_ranges(xB_ranges_set.begin(), xB_ranges_set.end());

    // Sort the xB_ranges vector
    std::sort(xB_ranges.begin(), xB_ranges.end());

    // Check if xB_bin index is within range
    if (xB_bin < 0 || xB_bin >= static_cast<int>(xB_ranges.size())) {
        std::cerr << "Error: xB_bin index out of range: " << xB_bin << std::endl;
        return relevant_bins;
    }

    // Get the xB range corresponding to xB_bin
    double xB_low = xB_ranges[xB_bin].first;
    double xB_high = xB_ranges[xB_bin].second;
    std::cout << "xB_bin range: [" << xB_low << ", " << xB_high << "]" << std::endl;

    // Select bins where bin.xB_low and bin.xB_high match xB_low and xB_high
    for (size_t bin_idx = 0; bin_idx < bin_boundaries.size(); ++bin_idx) {
        const auto& bin = bin_boundaries[bin_idx];

        if (bin.xB_low == xB_low && bin.xB_high == xB_high) {
            relevant_bins.push_back(bin_idx);
            // std::cout << "Bin_idx " << bin_idx << " added to relevant_bins." << std::endl;
        }
    }

    // std::cout << "Total relevant bins found: " << relevant_bins.size() << std::endl;
    return relevant_bins;
}

// Find the next perfect square greater than or equal to the given number
int next_perfect_square(int n) {
    int square_root = std::ceil(std::sqrt(n));
    return square_root * square_root;
}