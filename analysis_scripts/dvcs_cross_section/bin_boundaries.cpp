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

void create_directories(const std::string& base_output_dir) {
    // Define the necessary subdirectories
    std::string exclusivity_dir = base_output_dir + "/exclusivity_plots";
    std::string dvcs_dir = exclusivity_dir + "/dvcs";
    std::string eppi0_dir = exclusivity_dir + "/eppi0";
    std::string pi0_mass_dir = base_output_dir + "/pi0_mass";
    std::string comparison_dir = base_output_dir + "/data_mc_comparison";
    std::string dvcs_comparison_dir = comparison_dir + "/dvcs";
    std::string eppi0_comparison_dir = comparison_dir + "/eppi0";

    // Create the directories if they don't exist
    std::vector<std::string> dirs = {base_output_dir, exclusivity_dir, dvcs_dir, eppi0_dir, pi0_mass_dir, comparison_dir, dvcs_comparison_dir, eppi0_comparison_dir};
    for (const auto& dir : dirs) {
        if (!fs::exists(dir)) {
            if (fs::create_directory(dir)) {
                std::cout << "Created directory: " << dir << std::endl;
            } else {
                std::cerr << "Error: Failed to create directory: " << dir << std::endl;
            }
        }
    }
}