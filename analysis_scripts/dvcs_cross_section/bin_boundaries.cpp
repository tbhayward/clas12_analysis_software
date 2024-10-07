#include "bin_boundaries.h"
#include <fstream>
#include <sstream>
#include <iostream>

// Function definition to read bin boundaries from CSV
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename) {
    std::vector<BinBoundary> bin_boundaries;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bin_boundaries;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        BinBoundary bin;

        // Assuming CSV format: xB_low,xB_high,Q2_low,Q2_high,t_low,t_high
        std::getline(ss, token, ','); bin.xB_low = std::stod(token);
        std::getline(ss, token, ','); bin.xB_high = std::stod(token);
        std::getline(ss, token, ','); bin.Q2_low = std::stod(token);
        std::getline(ss, token, ','); bin.Q2_high = std::stod(token);
        std::getline(ss, token, ','); bin.t_low = std::stod(token);
        std::getline(ss, token, ','); bin.t_high = std::stod(token);

        bin_boundaries.push_back(bin);
    }

    file.close();
    return bin_boundaries;
}

// Function to map xB, Q2, and t to a bin index
int map_to_bin(const std::vector<BinBoundary>& bin_boundaries, double xB, double Q2, double t) {
    for (size_t i = 0; i < bin_boundaries.size(); ++i) {
        const BinBoundary& bin = bin_boundaries[i];

        if (xB >= bin.xB_low && xB < bin.xB_high &&
            Q2 >= bin.Q2_low && Q2 < bin.Q2_high &&
            t >= bin.t_low && t < bin.t_high) {
            return static_cast<int>(i);
        }
    }
    return -1;  // Return -1 if the bin is not found
}