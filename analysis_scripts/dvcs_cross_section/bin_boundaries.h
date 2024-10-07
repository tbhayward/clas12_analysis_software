#ifndef BIN_BOUNDARIES_H
#define BIN_BOUNDARIES_H

#include <string>
#include <vector>

// Structure to store bin boundaries
struct BinBoundary {
    double xB_low, xB_high;
    double Q2_low, Q2_high;
    double t_low, t_high;
};

// Function declarations
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename);

// Function to map xB, Q2, and t to a bin
int map_to_bin(const std::vector<BinBoundary>& bin_boundaries, double xB, double Q2, double t);

#endif