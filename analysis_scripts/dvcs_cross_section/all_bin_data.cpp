#include "all_bin_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>  // for rounding precision
#include <vector>
#include <string>

// Function to read the bin data from CSV
std::vector<AllBinData> read_bin_data(const std::string& filename) {
    std::vector<AllBinData> bin_data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return bin_data;
    }

    std::string line;
    int line_num = 0;

    // Skip the first line (header)
    std::getline(file, line);

    // Read data from the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        AllBinData bin;

        // First column: bin name (convert it to an int and store as bin_number)
        std::getline(ss, token, ',');
        bin.bin_number = std::stoi(token);  // Convert string to int

        // xB_min, xB_max, xB_avg
        std::getline(ss, token, ',');
        bin.xB_min = std::stod(token);
        std::getline(ss, token, ',');
        bin.xB_max = std::stod(token);
        std::getline(ss, token, ',');
        bin.xB_avg = std::stod(token);

        // Q2_min, Q2_max, Q2_avg
        std::getline(ss, token, ',');
        bin.Q2_min = std::stod(token);
        std::getline(ss, token, ',');
        bin.Q2_max = std::stod(token);
        std::getline(ss, token, ',');
        bin.Q2_avg = std::stod(token);

        // t_min, t_max, t_avg
        std::getline(ss, token, ',');
        bin.t_min = std::stod(token);
        std::getline(ss, token, ',');
        bin.t_max = std::stod(token);
        std::getline(ss, token, ',');
        bin.t_avg = std::stod(token);

        // phi_min, phi_max, phi_avg
        std::getline(ss, token, ',');
        bin.phi_min = std::stod(token);
        std::getline(ss, token, ',');
        bin.phi_max = std::stod(token);
        std::getline(ss, token, ',');
        bin.phi_avg = std::stod(token);

        // Round all values to four decimal places
        bin.xB_min = std::round(bin.xB_min * 10000) / 10000;
        bin.xB_max = std::round(bin.xB_max * 10000) / 10000;
        bin.xB_avg = std::round(bin.xB_avg * 10000) / 10000;

        bin.Q2_min = std::round(bin.Q2_min * 10000) / 10000;
        bin.Q2_max = std::round(bin.Q2_max * 10000) / 10000;
        bin.Q2_avg = std::round(bin.Q2_avg * 10000) / 10000;

        bin.t_min = std::round(bin.t_min * 10000) / 10000;
        bin.t_max = std::round(bin.t_max * 10000) / 10000;
        bin.t_avg = std::round(bin.t_avg * 10000) / 10000;

        bin.phi_min = std::round(bin.phi_min * 10000) / 10000;
        bin.phi_max = std::round(bin.phi_max * 10000) / 10000;
        bin.phi_avg = std::round(bin.phi_avg * 10000) / 10000;

        // Add bin to vector
        bin_data.push_back(bin);
    }

    file.close();
    return bin_data;
}

// Function to print the bin data for debugging
void print_bin_data(const std::vector<AllBinData>& bin_data) {
    for (const auto& bin : bin_data) {
        std::cout << "Bin Number: " << bin.bin_number
                  << ", xB [" << bin.xB_min << ", " << bin.xB_max << "], Avg: " << bin.xB_avg
                  << ", Q2 [" << bin.Q2_min << ", " << bin.Q2_max << "], Avg: " << bin.Q2_avg
                  << ", t [" << bin.t_min << ", " << bin.t_max << "], Avg: " << bin.t_avg
                  << ", Phi [" << bin.phi_min << ", " << bin.phi_max << "], Avg: " << bin.phi_avg
                  << std::endl;
    }
}