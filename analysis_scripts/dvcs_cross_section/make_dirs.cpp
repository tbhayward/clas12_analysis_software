#include "make_dirs.h"
#include <iostream>
#include <filesystem>
#include <vector>
#include <string>

namespace fs = std::filesystem;

void makeOutputDirs() {
    struct Category { std::string base;
                       std::vector<std::string> subs; };
    std::vector<Category> cats = {
        {"output/jsons", {}},
        {"output/exclusivity_plots", {}},
        {"output/mean_kinematic_plots", {}},
        {"output/contamination_plots",
            {"sp18_inb","sp18_out","fa18_inb","fa18_out","sp19_inb"}},
        {"output/bsa_plots",
            {"sp18_inb","sp18_out","fa18_inb_supplemental",
             "sp18_inb","sp18_out","10.6GeV","sp19_inb"}},
        {"output/acceptance",
            {"fa18_inb","fa18_out","sp19_inb"}}
    };

    for (const auto& cat : cats) {
        std::vector<std::string> created;
        if (!fs::exists(cat.base)) {
            fs::create_directories(cat.base);
            created.push_back(cat.base);
        }
        for (const auto& sub : cat.subs) {
            std::string path = cat.base + "/" + sub;
            if (!fs::exists(path)) {
                fs::create_directories(path);
                created.push_back(path);
            }
        }
        if (!created.empty()) {
            std::cout << "[Created] ";
            for (size_t i = 0; i < created.size(); ++i) {
                std::cout << created[i];
                if (i + 1 < created.size()) std::cout << ",";
            }
            std::cout << "
";
        }
    }
}