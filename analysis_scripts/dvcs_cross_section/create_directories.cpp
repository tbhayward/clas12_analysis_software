#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>

// Use the filesystem library
namespace fs = std::filesystem;

void create_directories(const std::string& base_output_dir) {
    // Define the necessary subdirectories
    std::string exclusivity_dir = base_output_dir + "/exclusivity_plots";
    std::string dvcs_dir = exclusivity_dir + "/dvcs";
    std::string eppi0_dir = exclusivity_dir + "/eppi0";
    std::string pi0_mass_dir = base_output_dir + "/pi0_mass";
    std::string data_mc_comparison_dir = base_output_dir + "/data_mc_comparison/dvcs";

    // Additional directories for the unfolding plots
    std::string unfolded_dir = base_output_dir + "/unfolded";
    std::string unfolded_dvcs_dir = unfolded_dir + "/dvcs";
    std::string unfolded_eppi0_dir = unfolded_dir + "/eppi0";

    // Check and create the directories if they don't exist
    if (!fs::exists(base_output_dir)) {
        if (fs::create_directories(base_output_dir)) {
            std::cout << "Created directory: " << base_output_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create base directory: " << base_output_dir << std::endl;
        }
    }

    if (!fs::exists(exclusivity_dir)) {
        if (fs::create_directories(exclusivity_dir)) {
            std::cout << "Created directory: " << exclusivity_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create exclusivity_plots directory: " << exclusivity_dir << std::endl;
        }
    }

    if (!fs::exists(dvcs_dir)) {
        if (fs::create_directories(dvcs_dir)) {
            std::cout << "Created directory: " << dvcs_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create dvcs directory: " << dvcs_dir << std::endl;
        }
    }

    if (!fs::exists(eppi0_dir)) {
        if (fs::create_directories(eppi0_dir)) {
            std::cout << "Created directory: " << eppi0_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create eppi0 directory: " << eppi0_dir << std::endl;
        }
    }

    if (!fs::exists(pi0_mass_dir)) {
        if (fs::create_directories(pi0_mass_dir)) {
            std::cout << "Created directory: " << pi0_mass_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create pi0_mass directory: " << pi0_mass_dir << std::endl;
        }
    }

    if (!fs::exists(data_mc_comparison_dir)) {
        if (fs::create_directories(data_mc_comparison_dir)) {
            std::cout << "Created directory: " << data_mc_comparison_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create data_mc_comparison/dvcs directory: " << data_mc_comparison_dir << std::endl;
        }
    }

    // Create directories for unfolded plots
    if (!fs::exists(unfolded_dir)) {
        if (fs::create_directories(unfolded_dir)) {
            std::cout << "Created directory: " << unfolded_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create unfolded directory: " << unfolded_dir << std::endl;
        }
    }

    if (!fs::exists(unfolded_dvcs_dir)) {
        if (fs::create_directories(unfolded_dvcs_dir)) {
            std::cout << "Created directory: " << unfolded_dvcs_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create unfolded/dvcs directory: " << unfolded_dvcs_dir << std::endl;
        }
    }

    if (!fs::exists(unfolded_eppi0_dir)) {
        if (fs::create_directories(unfolded_eppi0_dir)) {
            std::cout << "Created directory: " << unfolded_eppi0_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create unfolded/eppi0 directory: " << unfolded_eppi0_dir << std::endl;
        }
    }
}