#ifndef ALL_BIN_DATA_H
#define ALL_BIN_DATA_H

#include <vector>
#include <string>

// Function to read only the bin names from the CSV
std::vector<std::string> read_bin_names(const std::string& filename);

// Function to print the bin names for debugging
void print_bin_names(const std::vector<std::string>& bin_names);

#endif  // ALL_BIN_DATA_H