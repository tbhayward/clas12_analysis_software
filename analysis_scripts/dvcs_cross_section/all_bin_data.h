#ifndef ALL_BIN_DATA_H
#define ALL_BIN_DATA_H

#include <vector>
#include <string>

// Struct to hold only the bin number for now
struct AllBinData {
    int bin_number;    // Bin number (to track the row)
    std::string bin_name;  // Bin name from the second column
};

// Function to read and print the bin names from the CSV file
std::vector<AllBinData> read_bin_names(const std::string& filename);

// Function to print the bin names for debugging
void print_bin_names(const std::vector<AllBinData>& all_bin_data);

#endif  // ALL_BIN_DATA_H