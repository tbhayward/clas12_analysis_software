#ifndef ALL_BIN_DATA_H
#define ALL_BIN_DATA_H

#include <vector>
#include <string>

// Struct to hold bin data (bin number, xB, Q2, t, phi values)
struct AllBinData {
    int bin_number;  
    double xB_min;
    double xB_max;
    double xB_avg;
    double Q2_min;
    double Q2_max;
    double Q2_avg;
    double t_min;
    double t_max;
    double t_avg;
    double phi_min;
    double phi_max;
    double phi_avg;
};

// Function to read the bin data from CSV
std::vector<AllBinData> read_bin_data(const std::string& filename);

// Function to print the bin data for debugging
void print_bin_data(const std::vector<AllBinData>& bin_data);

#endif  // ALL_BIN_DATA_H