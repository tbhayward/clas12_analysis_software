#include "all_bin_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Function to read only the bin name (second column) from the CSV file
std::vector<AllBinData> read_bin_names(const std::string& filename) {
    std::vector<AllBinData> all_bin_data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return all_bin_data;
    }

    std::string line;
    int bin_number = 0;

    // Skip the first row (header)
    std::getline(file, line);

    // Read data from the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        AllBinData bin;

        // Increment bin number to track row
        bin.bin_number = ++bin_number;

        // Get the second column (bin name)
        // Skip the first column
        std::getline(ss, token, '\t');  // First column (skipped)

        // Read the second column (bin name)
        std::getline(ss, token, '\t');  // Second column (bin name)
        bin.bin_name = token;

        // Add this bin entry to the list
        all_bin_data.push_back(bin);
    }

    file.close();
    return all_bin_data;
}

// Function to print the bin names for debugging
void print_bin_names(const std::vector<AllBinData>& all_bin_data) {
    for (const auto& bin : all_bin_data) {
        std::cout << "Bin: " << bin.bin_number << ", Bin Name: " << bin.bin_name << std::endl;
    }
}