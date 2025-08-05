#include "make_dirs.h"
#include <iostream>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

void makeOutputDirs() {
    // Base output folders
    std::vector<std::string> bases = {
        "output/jsons", "output/exclusivity_plots", "output/mean_kinematic_plots"
    };
    for (auto& dir : bases) {
        if (!fs::exists(dir)) {
            std::cout << "Creating directory: " << dir << std::endl;
            fs::create_directories(dir);
        }
    }

    // contamination_plots subdirectories
    std::string contBase = "output/contamination_plots";
    if (!fs::exists(contBase)) fs::create_directories(contBase);
    std::vector<std::string> contSubs = {
        "sp18_inb", "sp18_out", "fa18_inb", "fa18_out", "sp19_inb"
    };
    for (auto& sub : contSubs) {
        std::string path = contBase + "/" + sub;
        if (!fs::exists(path)) {
            std::cout << "Creating directory: " << path << std::endl;
            fs::create_directories(path);
        }
    }

    // bsa_plots subdirectories
    std::string bsaBase = "output/bsa_plots";
    if (!fs::exists(bsaBase)) fs::create_directories(bsaBase);
    std::vector<std::string> bsaSubs = {
        "sp18_inb", "sp18_out", "fa18_inb_supplemental", "fa18_inb", "fa18_out", "10.6GeV", 
        "sp19_inb"
    };
    for (auto& sub : bsaSubs) {
        std::string path = bsaBase + "/" + sub;
        if (!fs::exists(path)) {
            std::cout << "Creating directory: " << path << std::endl;
            fs::create_directories(path);
        }
    }

    // acceptance subdirectories
    std::string accBase = "output/acceptance";
    if (!fs::exists(accBase)) fs::create_directories(accBase);
    std::vector<std::string> accSubs = {
        "fa18_inb", "fa18_out", "sp19_inb"
    };
    for (auto& sub : accSubs) {
        std::string path = accBase + "/" + sub;
        if (!fs::exists(path)) {
            std::cout << "Creating directory: " << path << std::endl;
            fs::create_directories(path);
        }
    }
}