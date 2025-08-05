#include "make_dirs.h"
#include <iostream>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

void makeOutputDirs() {
    // Base output folders
    std::vector<std::string> bases = {
        "output/jsons","output/exclusivity_plots",
        "output/mean_kinematic_plots"
    };
    std::vector<std::string> created;
    for (auto& dir : bases) {
        if (!fs::exists(dir)) {
            fs::create_directories(dir);
            created.push_back(dir);
        }
    }
    if (!created.empty()) {
        std::cout << "[Created] ";
        for (size_t i = 0; i < created.size(); ++i) {
            std::cout << created[i] << (i+1<created.size()?",":"");
        }
        std::cout << "
";
    }

    // contamination_plots subdirectories
    std::string contBase = "output/contamination_plots";
    if (!fs::exists(contBase)) fs::create_directories(contBase);
    created.clear();
    std::vector<std::string> contSubs = {"sp18_inb","sp18_out","fa18_inb","fa18_out","sp19_inb"};
    for (auto& sub : contSubs) {
        std::string path = contBase + "/" + sub;
        if (!fs::exists(path)) {
            fs::create_directories(path);
            created.push_back(path);
        }
    }
    if (!created.empty()) {
        std::cout << "[Created] ";
        for (size_t i = 0; i < created.size(); ++i) {
            std::cout << created[i] << (i+1<created.size()?",":"");
        }
        std::cout << "
";
    }

    // bsa_plots subdirectories
    std::string bsaBase = "output/bsa_plots";
    if (!fs::exists(bsaBase)) fs::create_directories(bsaBase);
    created.clear();
    std::vector<std::string> bsaSubs = {"sp18_inb","sp18_out","fa18_inb_supplemental",
                                        "fa18_inb","fa18_out","10.6GeV","sp19_inb"};
    for (auto& sub : bsaSubs) {
        std::string path = bsaBase + "/" + sub;
        if (!fs::exists(path)) {
            fs::create_directories(path);
            created.push_back(path);
        }
    }
    if (!created.empty()) {
        std::cout << "[Created] ";
        for (size_t i = 0; i < created.size(); ++i) {
            std::cout << created[i] << (i+1<created.size()?",":"");
        }
        std::cout << "
";
    }

    // acceptance subdirectories
    std::string accBase = "output/acceptance";
    if (!fs::exists(accBase)) fs::create_directories(accBase);
    created.clear();
    std::vector<std::string> accSubs = {"fa18_inb","fa18_out","sp19_inb"};
    for (auto& sub : accSubs) {
        std::string path = accBase + "/" + sub;
        if (!fs::exists(path)) {
            fs::create_directories(path);
            created.push_back(path);
        }
    }
    if (!created.empty()) {
        std::cout << "[Created] ";
        for (size_t i = 0; i < created.size(); ++i) {
            std::cout << created[i] << (i+1<created.size()?",":"");
        }
        std::cout << "
";
    }
}