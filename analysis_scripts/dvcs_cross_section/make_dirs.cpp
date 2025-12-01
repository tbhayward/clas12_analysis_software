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
        {"output/csvs"},
        {"output/jsons"},
        {"output/exclusivity_plots", {}},
        {"output/total_counts_plots",
            {"Fa18_Inb", "Fa18_Out", "Fa18_Inb_Supp", "Sp18_Inb",
                "Sp18_Out", "Sp19_Inb", "Fa18", "Sp18", "10.6 GeV"
            }
        },
        {"output/contamination_plots",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb",
                "Sp18_Out", "Sp19_Inb"
            }
        },
        {"output/signal_yield_plots",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb",
                "Sp18_Out", "Sp19_Inb"
            }
        },
        {"output/acceptance",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb",
                "Sp18_Out", "Sp19_Inb"
            }
        },
        {"output/unfolded_yields",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb",
                "Sp18_Out", "Sp19_Inb", "Fa18", "Sp18",
                "10.6_GeV"
            }
        },
        {"output/radiative_correction_plots",
            {"10.60","10.2"}},
        {"output/bin_volume",
            {"10.60","10.2"}},
        {"output/bin_centering_plots",
            {"10.60","10.2"}},
        {"output/cross_sections",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb",
                "Sp18_Out", "Sp19_Inb", "Fa18", "Sp18",
                "10.6_GeV"
            }
        },
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