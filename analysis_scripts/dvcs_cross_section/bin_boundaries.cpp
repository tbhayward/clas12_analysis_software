// read_bin_boundaries.cpp

#include "bin_boundaries.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::vector<BinBoundary> read_bin_boundaries(const std::string& filename) {
    std::vector<BinBoundary> bin_boundaries;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open binning file: " << filename << std::endl;
        return bin_boundaries;
    }

    std::string line;
    // Read header line (the first line with column names)
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string item;
        BinBoundary bin;

        // Read unnamed (first column)
        if (!std::getline(ss, item, ',')) continue;
        bin.unnamed = std::stoi(item);

        // Read bin name (second column)
        if (!std::getline(ss, item, ',')) continue;
        bin.bin_name = std::stoi(item);

        // Read xBmin (xB_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.xB_low = std::stod(item);

        // Read xBmax (xB_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.xB_high = std::stod(item);

        // Read xBavg
        if (!std::getline(ss, item, ',')) continue;
        bin.xB_avg = std::stod(item);

        // Read Q2min (Q2_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.Q2_low = std::stod(item);

        // Read Q2max (Q2_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.Q2_high = std::stod(item);

        // Read Q2avg
        if (!std::getline(ss, item, ',')) continue;
        bin.Q2_avg = std::stod(item);

        // Read t_abs_min (t_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.t_low = std::stod(item);

        // Read t_abs_max (t_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.t_high = std::stod(item);

        // Read t_abs_avg
        if (!std::getline(ss, item, ',')) continue;
        bin.t_avg = std::stod(item);

        // Read phimin (phi_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.phi_low = std::stod(item);

        // Read phimax (phi_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.phi_high = std::stod(item);

        // Read phiavg
        if (!std::getline(ss, item, ',')) continue;
        bin.phi_avg = std::stod(item);

        // Add the bin to the vector
        bin_boundaries.push_back(bin);
    }

    infile.close();
    return bin_boundaries;
}