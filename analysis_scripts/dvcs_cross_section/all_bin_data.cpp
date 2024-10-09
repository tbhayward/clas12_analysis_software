#include "all_bin_data.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

// Function to read all_bin_v3.csv
std::vector<AllBinData> read_all_bin_v3(const std::string& filename) {
    std::vector<AllBinData> all_bin_data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return all_bin_data;
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

        // Skip the first column (bin name)
        std::getline(ss, token, '\t');

        // Second column: bin number (int)
        std::getline(ss, token, '\t');
        bin.bin_number = std::stoi(token);

        // Skip columns until phi min (12th column)
        for (int i = 0; i < 10; ++i) {
            std::getline(ss, token, '\t');
        }

        // Phi min (12th column)
        std::getline(ss, token, '\t');
        bin.phi_min = std::stoi(token);

        // Phi max (13th column)
        std::getline(ss, token, '\t');
        bin.phi_max = std::stoi(token);

        // Phi avg (14th column)
        std::getline(ss, token, '\t');
        bin.phi_avg = std::stod(token);

        // Skip to yield columns (starting from 22nd column)
        for (int i = 0; i < 7; ++i) {
            std::getline(ss, token, '\t');
        }

        // Yields for epg (22nd to 27th columns)
        std::getline(ss, token, '\t');
        bin.yield_epg_FD_FD_inb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_epg_CD_FD_inb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_epg_CD_FT_inb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_epg_FD_FD_outb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_epg_CD_FD_outb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_epg_CD_FT_outb = std::stoi(token);

        // Yields for eppi0 (28th to 33rd columns)
        std::getline(ss, token, '\t');
        bin.yield_eppi0_FD_FD_inb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_eppi0_CD_FD_inb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_eppi0_CD_FT_inb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_eppi0_FD_FD_outb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_eppi0_CD_FD_outb = std::stoi(token);
        std::getline(ss, token, '\t');
        bin.yield_eppi0_CD_FT_outb = std::stoi(token);

        // Add to the list of bins
        all_bin_data.push_back(bin);
    }

    file.close();
    return all_bin_data;
}

// Function to print the data line by line for debugging
void print_all_bin_data(const std::vector<AllBinData>& all_bin_data) {
    for (const auto& bin : all_bin_data) {
        std::cout << "Bin: " << bin.bin_number
                  << ", Phi Range: [" << bin.phi_min << ", " << bin.phi_max << "]"
                  << ", Phi Avg: " << bin.phi_avg
                  << ", Yields EPG (Inb): " << bin.yield_epg_FD_FD_inb << ", " << bin.yield_epg_CD_FD_inb << ", " << bin.yield_epg_CD_FT_inb
                  << ", Yields EPG (Outb): " << bin.yield_epg_FD_FD_outb << ", " << bin.yield_epg_CD_FD_outb << ", " << bin.yield_epg_CD_FT_outb
                  << ", Yields EPPI0 (Inb): " << bin.yield_eppi0_FD_FD_inb << ", " << bin.yield_eppi0_CD_FD_inb << ", " << bin.yield_eppi0_CD_FT_inb
                  << ", Yields EPPI0 (Outb): " << bin.yield_eppi0_FD_FD_outb << ", " << bin.yield_eppi0_CD_FD_outb << ", " << bin.yield_eppi0_CD_FT_outb
                  << std::endl;
    }
}