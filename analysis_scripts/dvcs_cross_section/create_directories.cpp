// create_directories.cpp

#include <filesystem>
#include <iostream>
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
    std::string data_mc_comparison_dir = base_output_dir + "/data_mc_comparison";
    std::string dvcs_comparison_dir = data_mc_comparison_dir + "/dvcs";
    std::string eppi0_comparison_dir = data_mc_comparison_dir + "/eppi0";

    // Additional directories for the unfolding plots
    std::string unfolded_dir = base_output_dir + "/unfolded";
    std::string unfolded_dvcs_dir = unfolded_dir + "/dvcs";
    std::string unfolded_eppi0_dir = unfolded_dir + "/eppi0";

    // Directories for specific periods: Fa18Inb, Fa18Out, Sp19Inb (for dvcs)
    std::string fa18inb_yields_dir = unfolded_dvcs_dir + "/Fa18Inb/yields";
    std::string fa18inb_acceptances_dir = unfolded_dvcs_dir + "/Fa18Inb/acceptances";
    std::string fa18out_yields_dir = unfolded_dvcs_dir + "/Fa18Out/yields";
    std::string fa18out_acceptances_dir = unfolded_dvcs_dir + "/Fa18Out/acceptances";
    std::string sp19inb_yields_dir = unfolded_dvcs_dir + "/Sp19Inb/yields";
    std::string sp19inb_acceptances_dir = unfolded_dvcs_dir + "/Sp19Inb/acceptances";

    // eppi0 directories for each period
    std::string fa18inb_eppi0_yields_dir = unfolded_eppi0_dir + "/Fa18Inb/yields";
    std::string fa18inb_eppi0_acceptances_dir = unfolded_eppi0_dir + "/Fa18Inb/acceptances";
    std::string fa18out_eppi0_yields_dir = unfolded_eppi0_dir + "/Fa18Out/yields";
    std::string fa18out_eppi0_acceptances_dir = unfolded_eppi0_dir + "/Fa18Out/acceptances";
    std::string sp19inb_eppi0_yields_dir = unfolded_eppi0_dir + "/Sp19Inb/yields";
    std::string sp19inb_eppi0_acceptances_dir = unfolded_eppi0_dir + "/Sp19Inb/acceptances";

    // Directory for contamination plots
    std::string contamination_plots_dir = base_output_dir + "/contamination_plots";

    // Check and create the directories if they don't exist
    std::vector<std::string> directories = {
        exclusivity_dir, dvcs_dir, eppi0_dir, pi0_mass_dir,
        dvcs_comparison_dir, eppi0_comparison_dir,
        unfolded_dvcs_dir, unfolded_eppi0_dir,
        fa18inb_yields_dir, fa18inb_acceptances_dir, fa18out_yields_dir, fa18out_acceptances_dir, sp19inb_yields_dir, sp19inb_acceptances_dir,
        fa18inb_eppi0_yields_dir, fa18inb_eppi0_acceptances_dir, fa18out_eppi0_yields_dir, fa18out_eppi0_acceptances_dir, sp19inb_eppi0_yields_dir, sp19inb_eppi0_acceptances_dir,
        contamination_plots_dir  // Added contamination_plots_dir
    };

    // Loop through the directories and create them if necessary
    for (const auto& dir : directories) {
        if (!fs::exists(dir)) {
            if (fs::create_directories(dir)) {
                std::cout << "Created directory: " << dir << std::endl;
            } else {
                std::cerr << "Error: Failed to create directory: " << dir << std::endl;
            }
        }
    }
}