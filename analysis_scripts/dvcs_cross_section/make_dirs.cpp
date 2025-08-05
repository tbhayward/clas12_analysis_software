#include "make_dirs.h"
#include <iostream>
#include <filesystem>
#include <vector>
namespace fs = std::filesystem;

void makeOutputDirs() {
    std::vector<std::string> bases = {
        "output/jsons","output/exclusivity_plots",
        "output/mean_kinematic_plots"
    };
    for (auto& dir : bases)
        if (!fs::exists(dir)) {
            std::cout << "Creating directory: " << dir << "
";
            fs::create_directories(dir);
        }

    std::string contBase = "output/contamination_plots";
    if (!fs::exists(contBase)) fs::create_directories(contBase);
    for (auto& sub : {"sp18_inb","sp18_out","fa18_inb","fa18_out","sp19_inb"}) {
        std::string path = contBase + "/" + sub;
        if (!fs::exists(path)) {
            std::cout << "Creating directory: " << path << "
";
            fs::create_directories(path);
        }
    }

    std::string bsaBase = "output/bsa_plots";
    if (!fs::exists(bsaBase)) fs::create_directories(bsaBase);
    for (auto& sub : {"sp18_inb","sp18_out","fa18_inb_supplemental",
                      "sp18_inb","fa18_out","10.6GeV","sp19_inb"}) {
        std::string path = bsaBase + "/" + sub;
        if (!fs::exists(path)) {
            std::cout << "Creating directory: " << path << "
";
            fs::create_directories(path);
        }
    }

    std::string accBase = "output/acceptance";
    if (!fs::exists(accBase)) fs::create_directories(accBase);
    for (auto& sub : {"fa18_inb","fa18_out","sp19_inb"}) {
        std::string path = accBase + "/" + sub;
        if (!fs::exists(path)) {
            std::cout << "Creating directory: " << path << "
";
            fs::create_directories(path);
        }
    }
}