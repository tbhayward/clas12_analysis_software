#include "all_bin_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>  // For std::round

// Function to read all_bin_v3.csv
std::vector<AllBinData> read_bin_data(const std::string& filename) {
    std::vector<AllBinData> all_bin_data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return all_bin_data;
    }

    std::string line;

    // Skip the first line (header)
    std::getline(file, line);

    // Read data from the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        AllBinData bin;

        // First column: bin number
        std::getline(ss, token, ',');
        bin.bin_number = std::stoi(token);

        // Second column: bin name (convert it to an integer)
        std::getline(ss, token, ',');
        bin.bin_name = std::stoi(token);

        // xB_min (3rd column)
        std::getline(ss, token, ',');
        bin.xB_min = std::stod(token);

        // xB_max (4th column)
        std::getline(ss, token, ',');
        bin.xB_max = std::stod(token);

        // xB_avg (5th column)
        std::getline(ss, token, ',');
        bin.xB_avg = std::stod(token);

        // Q2_min (6th column)
        std::getline(ss, token, ',');
        bin.Q2_min = std::stod(token);

        // Q2_max (7th column)
        std::getline(ss, token, ',');
        bin.Q2_max = std::stod(token);

        // Q2_avg (8th column)
        std::getline(ss, token, ',');
        bin.Q2_avg = std::stod(token);

        // t_min (9th column)
        std::getline(ss, token, ',');
        bin.t_min = std::stod(token);

        // t_max (10th column)
        std::getline(ss, token, ',');
        bin.t_max = std::stod(token);

        // t_avg (11th column)
        std::getline(ss, token, ',');
        bin.t_avg = std::stod(token);

        // phi_min (12th column)
        std::getline(ss, token, ',');
        bin.phi_min = std::stod(token);

        // phi_max (13th column)
        std::getline(ss, token, ',');
        bin.phi_max = std::stod(token);

        // phi_avg (14th column)
        std::getline(ss, token, ',');
        bin.phi_avg = std::stod(token);

        // Skip the next few columns until we reach the yields (starting from 22nd column)
        for (int i = 0; i < 7; ++i) {
            std::getline(ss, token, ',');
        }

        // Read the yield data (22nd to 33rd columns)
        std::getline(ss, token, ',');
        bin.yield_epg_FD_FD_inb = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_epg_CD_FD_inb = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_epg_CD_FT_inb = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_epg_FD_FD_out = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_epg_CD_FD_out = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_epg_CD_FT_out = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_eppi0_FD_FD_inb = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_eppi0_CD_FD_inb = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_eppi0_CD_FT_inb = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_eppi0_FD_FD_out = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_eppi0_CD_FD_out = std::stoi(token);

        std::getline(ss, token, ',');
        bin.yield_eppi0_CD_FT_out = std::stoi(token);

        // Round all relevant values to 4 decimal places
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

        // Add the bin to the list
        all_bin_data.push_back(bin);
    }

    file.close();
    return all_bin_data;
}

// Function to print the bin data for debugging
void print_bin_data(const std::vector<AllBinData>& all_bin_data) {
    for (const auto& bin : all_bin_data) {
        std::cout << "Bin Number: " << bin.bin_number
                  << ", Bin Name: " << bin.bin_name
                  << ", xB Range: [" << bin.xB_min << ", " << bin.xB_max << "]"
                  << ", xB Avg: " << bin.xB_avg
                  << ", Q2 Range: [" << bin.Q2_min << ", " << bin.Q2_max << "]"
                  << ", Q2 Avg: " << bin.Q2_avg
                  << ", t Range: [" << bin.t_min << ", " << bin.t_max << "]"
                  << ", t Avg: " << bin.t_avg
                  << ", Phi Range: [" << bin.phi_min << ", " << bin.phi_max << "]"
                  << ", Phi Avg: " << bin.phi_avg
                  << ", Yields EPG (Inb): " << bin.yield_epg_FD_FD_inb << ", " << bin.yield_epg_CD_FD_inb << ", " << bin.yield_epg_CD_FT_inb
                  << ", Yields EPG (Out): " << bin.yield_epg_FD_FD_out << ", " << bin.yield_epg_CD_FD_out << ", " << bin.yield_epg_CD_FT_out
                  << ", Yields EPPI0 (Inb): " << bin.yield_eppi0_FD_FD_inb << ", " << bin.yield_eppi0_CD_FD_inb << ", " << bin.yield_eppi0_CD_FT_inb
                  << ", Yields EPPI0 (Out): " << bin.yield_eppi0_FD_FD_out << ", " << bin.yield_eppi0_CD_FD_out << ", " << bin.yield_eppi0_CD_FT_out
                  << std::endl;
    }
}