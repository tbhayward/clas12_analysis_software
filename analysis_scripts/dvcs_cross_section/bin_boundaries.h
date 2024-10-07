#ifndef BIN_BOUNDARIES_H
#define BIN_BOUNDARIES_H

#include <vector>
#include <string>

struct BinBoundary {
    std::string bin_label;  // Add this field to hold the bin label
    double xB_low, xB_high;
    double Q2_low, Q2_high;
    double t_low, t_high;
};

// Function to read bin boundaries from CSV
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename);

// Function to map xB, Q2, and t to a bin index
int map_to_bin(const std::vector<BinBoundary>& bin_boundaries, double xB, double Q2, double t);

#endif