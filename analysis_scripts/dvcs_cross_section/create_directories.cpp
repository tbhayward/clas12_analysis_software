#include "create_directories.h"
#include <filesystem>
#include <iostream>

// Use the filesystem library
namespace fs = std::filesystem;

void create_directories(const std::string& base_output_dir) {
    // Define the necessary subdirectories
    std::string exclusivity_dir = base_output_dir + "/exclusivity_plots";
    std::string dvcs_dir = exclusivity_dir + "/dvcs";
    std::string eppi0_dir = exclusivity_dir + "/eppi0";

    // Check and create the directories if they don't exist
    if (!fs::exists(base_output_dir)) {
        if (fs::create_directory(base_output_dir)) {
            std::cout << "Created directory: " << base_output_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create base directory: " << base_output_dir << std::endl;
        }
    }

    if (!fs::exists(exclusivity_dir)) {
        if (fs::create_directory(exclusivity_dir)) {
            std::cout << "Created directory: " << exclusivity_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create exclusivity_plots directory: " << exclusivity_dir << std::endl;
        }
    }

    if (!fs::exists(dvcs_dir)) {
        if (fs::create_directory(dvcs_dir)) {
            std::cout << "Created directory: " << dvcs_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create dvcs directory: " << dvcs_dir << std::endl;
        }
    }

    if (!fs::exists(eppi0_dir)) {
        if (fs::create_directory(eppi0_dir)) {
            std::cout << "Created directory: " << eppi0_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create eppi0 directory: " << eppi0_dir << std::endl;
        }
    }
}