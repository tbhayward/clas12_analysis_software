#ifndef BIN_BOUNDARIES_H
#define BIN_BOUNDARIES_H

#include <vector>
#include <string>

struct BinBoundary {
    int bin_index;         // Bin Index
    std::string bin_name;  // Bin Name
    double xB_low;         // xBmin
    double xB_high;        // xBmax
    double Q2_low;         // Q2min
    double Q2_high;        // Q2max
    double t_low;          // t_abs_min
    double t_high;         // t_abs_max
    double phi_low;        // phimin
    double phi_high;       // phimax
    // Add other members if necessary
};

// Function to read bin boundaries from CSV
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename);

// Function to map xB, Q2, t, and phi to a bin index (if needed)
int map_to_bin(const std::vector<BinBoundary>& bin_boundaries, double xB, double Q2, double t, double phi);

#endif