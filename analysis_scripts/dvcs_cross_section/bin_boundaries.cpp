#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <iostream>

struct BinBoundaries {
    std::vector<double> xB_bins;
    std::vector<double> Q2_bins;
    std::vector<double> t_bins;
    std::vector<double> phi_bins;
};

BinBoundaries read_bin_boundaries(const std::string& filename) {
    std::ifstream file(filename);
    std::string line;
    BinBoundaries boundaries;
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return boundaries;
    }

    while (getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::vector<double> current_bin;
        while (getline(ss, field, ',')) {
            current_bin.push_back(stod(field));
        }

        // Assume first row defines the boundaries
        boundaries.xB_bins.push_back(current_bin[0]);
        boundaries.Q2_bins.push_back(current_bin[1]);
        boundaries.t_bins.push_back(current_bin[2]);
        boundaries.phi_bins.push_back(current_bin[3]);
    }

    return boundaries;
}

// Function to assign an event to a bin based on its xB, Q2, and t values
std::tuple<int, int, int> assign_bin(double xB, double Q2, double t, const BinBoundaries& boundaries) {
    int xB_bin = -1, Q2_bin = -1, t_bin = -1;
    
    // Find the correct xB bin
    for (size_t i = 0; i < boundaries.xB_bins.size() - 1; ++i) {
        if (xB >= boundaries.xB_bins[i] && xB < boundaries.xB_bins[i + 1]) {
            xB_bin = i;
            break;
        }
    }

    // Find the correct Q2 bin
    for (size_t i = 0; i < boundaries.Q2_bins.size() - 1; ++i) {
        if (Q2 >= boundaries.Q2_bins[i] && Q2 < boundaries.Q2_bins[i + 1]) {
            Q2_bin = i;
            break;
        }
    }

    // Find the correct t bin
    for (size_t i = 0; i < boundaries.t_bins.size() - 1; ++i) {
        if (t >= boundaries.t_bins[i] && t < boundaries.t_bins[i + 1]) {
            t_bin = i;
            break;
        }
    }

    return std::make_tuple(xB_bin, Q2_bin, t_bin);
}