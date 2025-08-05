// bin_boundaries.h

#ifndef BIN_BOUNDARIES_H
#define BIN_BOUNDARIES_H

#include <string>
#include <vector>

struct BinBoundary {
    int unnamed;        // First column, unnamed int (global bin number)
    int bin_name;       // Second column, bin name (int, corresponds to 3D bin)
    double xB_low;      // xBmin
    double xB_high;     // xBmax
    double xB_avg;      // xBavg
    double Q2_low;      // Q2min
    double Q2_high;     // Q2max
    double Q2_avg;      // Q2avg
    double t_low;       // t_abs_min
    double t_high;      // t_abs_max
    double t_avg;       // t_abs_avg
    double phi_low;     // phimin
    double phi_high;    // phimax
    double phi_avg;     // phiavg
    // Add other members if necessary
};

// Function to read bin boundaries from CSV
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename);

#endif // BIN_BOUNDARIES_H