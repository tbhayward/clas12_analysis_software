// make_dirs.cpp
#include "make_dirs.h"
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

namespace fs = std::filesystem;

void makeOutputDirs() {
    struct Category { std::string base; std::vector<std::string> subs; };
    std::vector<Category> cats = {
        // jsons now has a subdir for the per-period/per-topology JSONs
        {"output/jsons", {"individual_cuts"}},
        {"output/exclusivity_plots", {}},
        {"output/mean_kinematic_plots", {}},
        {"output/bsa_plots",
            {"sp18_inb","sp18_out","fa18_inb_supplemental",
             "sp18_inb","sp18_out","10.6GeV","sp19_inb"}},
        {"output/acceptance",
            {"fa18_inb","fa18_out","sp19_inb"}}
    };

    for (const auto& cat : cats) {
        std::vector<std::string> created;

        // create base directory
        if (!fs::exists(cat.base)) {
            fs::create_directories(cat.base);
            created.push_back(cat.base);
        }
        // create subdirectories
        for (const auto& sub : cat.subs) {
            std::string path = cat.base + "/" + sub;
            if (!fs::exists(path)) {
                fs::create_directories(path);
                created.push_back(path);
            }
        }
        // print all created paths on one line
        if (!created.empty()) {
            std::cout << "[Created] ";
            for (size_t i = 0; i < created.size(); ++i) {
                std::cout << created[i];
                if (i + 1 < created.size()) std::cout << ",";
            }
            std::cout << std::endl;
        }
    }
}