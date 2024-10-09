#include <filesystem>
#include <iostream>
#include <string>

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

    // Directories for Yields and Acceptance plots
    std::string unfolded_yields_dir = unfolded_dvcs_dir + "/yields";
    std::string unfolded_acceptances_dir = unfolded_dvcs_dir + "/acceptances";

    // Directories for specific periods: Fa18Inb, Fa18Out, Sp19Inb
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

    // Create new directories for yields and acceptance plots
    if (!fs::exists(unfolded_yields_dir)) {
        if (fs::create_directories(unfolded_yields_dir)) {
            std::cout << "Created directory: " << unfolded_yields_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create yields directory: " << unfolded_yields_dir << std::endl;
        }
    }

    if (!fs::exists(unfolded_acceptances_dir)) {
        if (fs::create_directories(unfolded_acceptances_dir)) {
            std::cout << "Created directory: " << unfolded_acceptances_dir << std::endl;
        } else {
            std::cerr << "Error: Failed to create acceptances directory: " << unfolded_acceptances_dir << std::endl;
        }
    }

    // Create the specific period subdirectories for dvcs and eppi0 (Fa18Inb, Fa18Out, Sp19Inb)
    std::vector<std::string> period_dirs = {
        fa18inb_yields_dir, fa18inb_acceptances_dir,
        fa18out_yields_dir, fa18out_acceptances_dir,
        sp19inb_yields_dir, sp19inb_acceptances_dir,
        fa18inb_eppi0_yields_dir, fa18inb_eppi0_acceptances_dir,
        fa18out_eppi0_yields_dir, fa18out_eppi0_acceptances_dir,
        sp19inb_eppi0_yields_dir, sp19inb_eppi0_acceptances_dir
    };

    for (const auto& dir : period_dirs) {
        if (!fs::exists(dir)) {
            if (fs::create_directories(dir)) {
                std::cout << "Created directory: " << dir << std::endl;
            } else {
                std::cerr << "Error: Failed to create directory: " << dir << std::endl;
            }
        }
    }
}