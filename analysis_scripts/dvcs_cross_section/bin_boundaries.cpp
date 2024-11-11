// read_bin_boundaries.cpp

#include "bin_boundaries.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

std::vector<BinBoundary> read_bin_boundaries(const std::string& filename) {
    std::vector<BinBoundary> bin_boundaries;
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Unable to open binning file: " << filename << std::endl;
        return bin_boundaries;
    }

    std::string line;
    // Skip header line if present
    std::getline(infile, line); // Adjust this if your CSV does not have a header

    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        std::string item;
        BinBoundary bin;

        // Read and discard the first two columns (index and bin name)
        std::getline(ss, item, ','); // Skip index
        std::getline(ss, item, ','); // Skip bin name

        // Now read xBmin (xB_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.xB_low = std::stod(item);

        // xBmax (xB_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.xB_high = std::stod(item);

        // Skip xBavg
        if (!std::getline(ss, item, ',')) continue;

        // Q2min (Q2_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.Q2_low = std::stod(item);

        // Q2max (Q2_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.Q2_high = std::stod(item);

        // Skip Q2avg
        if (!std::getline(ss, item, ',')) continue;

        // t_abs_min (t_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.t_low = std::stod(item);

        // t_abs_max (t_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.t_high = std::stod(item);

        // Skip t_abs_avg
        if (!std::getline(ss, item, ',')) continue;

        // phimin (phi_low)
        if (!std::getline(ss, item, ',')) continue;
        bin.phi_low = std::stod(item);

        // phimax (phi_high)
        if (!std::getline(ss, item, ',')) continue;
        bin.phi_high = std::stod(item);

        // Optionally, skip phiavg
        if (!std::getline(ss, item, ',')) continue;

        // Convert phi values to degrees if necessary
        // If phi values are already in degrees, you can skip this
        // bin.phi_low *= RAD_TO_DEG;
        // bin.phi_high *= RAD_TO_DEG;

        bin_boundaries.push_back(bin);
    }

    infile.close();
    return bin_boundaries;
}