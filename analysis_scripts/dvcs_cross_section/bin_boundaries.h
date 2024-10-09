#ifndef BIN_BOUNDARIES_H
#define BIN_BOUNDARIES_H

#include <vector>
#include <string>

struct BinBoundary {
    int bin_number;          // Bin number
    std::string bin_label;   // Bin label for identification
    double xB_low, xB_high;  // xB min and max
    double Q2_low, Q2_high;  // Q2 min and max
    double t_low, t_high;    // |t| min and max
    double xB_avg;           // Average xB
    double Q2_avg;           // Average Q2
    double t_avg;            // Average |t|
};

// Function to read bin boundaries from CSV
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename);

// Function to map xB, Q2, and t to a bin index
int map_to_bin(const std::vector<BinBoundary>& bin_boundaries, double xB, double Q2, double t);

#endif