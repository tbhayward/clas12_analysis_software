#include "bin_boundaries.h"
#include <fstream>
#include <sstream>
#include <iostream>

// Function definition to read bin boundaries from tab-delimited CSV
std::vector<BinBoundary> read_bin_boundaries(const std::string& filename) {
    std::vector<BinBoundary> bin_boundaries;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bin_boundaries;
    }

    std::string line;
    int line_num = 0;

    // Skip first two lines (header and descriptive line)
    while (line_num < 2 && std::getline(file, line)) {
        ++line_num;
    }

    // Now read the actual data from the third line onward
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        BinBoundary bin;

        // First column: bin number (skip it)
        std::getline(ss, token, '\t');  // Skip bin number

        // Second column: bin name (read it)
        std::getline(ss, token, '\t');
        bin.bin_label = token;
        // std::cout << "Processing bin: " << bin.bin_label << std::endl;

        try {
            // xB_low (min)
            std::getline(ss, token, '\t');
            bin.xB_low = std::stod(token);

            // xB_high (max)
            std::getline(ss, token, '\t');
            bin.xB_high = std::stod(token);

            // Skip average xB
            std::getline(ss, token, '\t');

            // Q2_low (min)
            std::getline(ss, token, '\t');
            bin.Q2_low = std::stod(token);

            // Q2_high (max)
            std::getline(ss, token, '\t');
            bin.Q2_high = std::stod(token);

            // Skip average Q2
            std::getline(ss, token, '\t');

            // t_low (min)
            std::getline(ss, token, '\t');
            bin.t_low = std::stod(token);

            // t_high (max)
            std::getline(ss, token, '\t');
            bin.t_high = std::stod(token);

            // Skip average t
            std::getline(ss, token, '\t');

        } catch (const std::invalid_argument& e) {
            std::cerr << "Invalid argument while parsing token: '" << token << "' in bin " << bin.bin_label << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << "Out of range error while parsing token: '" << token << "' in bin " << bin.bin_label << std::endl;
        }

        bin_boundaries.push_back(bin);
    }

    // // Output test to verify
    // for (const auto& bin : bin_boundaries) {
    //     std::cout << "Bin: " << bin.bin_label
    //               << " xB [" << bin.xB_low << ", " << bin.xB_high << "]"
    //               << " Q2 [" << bin.Q2_low << ", " << bin.Q2_high << "]"
    //               << " t [" << bin.t_low << ", " << bin.t_high << "]"
    //               << std::endl;
    // }

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