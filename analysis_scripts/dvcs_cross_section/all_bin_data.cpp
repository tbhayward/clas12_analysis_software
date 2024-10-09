#include "all_bin_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Function to read the CSV and extract only the bin name
std::vector<std::string> read_bin_names(const std::string& filename) {
    std::vector<std::string> bin_names;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bin_names;
    }

    std::string line;
    // Skip the first line (header)
    std::getline(file, line);

    // Read data from the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;

        // Skip the first column (we don't need it)
        std::getline(ss, token, ',');

        // Second column: bin name (read it)
        std::getline(ss, token, ',');
        bin_names.push_back(token);  // Add bin name to the vector
    }

    file.close();
    return bin_names;
}

// Function to print the bin names for debugging
void print_bin_names(const std::vector<std::string>& bin_names) {
    for (const auto& name : bin_names) {
        std::cout << "Bin Name: " << name << std::endl;
    }
}