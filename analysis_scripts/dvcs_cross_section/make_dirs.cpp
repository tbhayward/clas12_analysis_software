// make_dirs.cpp
#include "make_dirs.h"

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

void makeOutputDirs() {
    struct Category {
        const char* base;
        std::vector<const char*> subs;
    };

    // Construct this immutable directory specification only once. The function
    // is idempotent, so a later call remains safe even though main.cpp now calls
    // it only once at startup.
    static const std::vector<Category> categories = {
        {"output/csvs", {}},
        {"output/jsons", {}},
        {"output/exclusivity_plots", {}},
        {"output/total_counts_plots",
            {"Fa18_Inb", "Fa18_Out", "Fa18_Inb_Supp", "Sp18_Inb",
             "Sp18_Out", "Sp19_Inb", "Fa18", "Sp18", "10.6 GeV"}},
        {"output/contamination_plots",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb", "Sp18_Out", "Sp19_Inb"}},
        {"output/signal_yield_plots",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb", "Sp18_Out", "Sp19_Inb"}},
        {"output/acceptance",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb", "Sp18_Out", "Sp19_Inb"}},
        {"output/unfolded_yields",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb", "Sp18_Out", "Sp19_Inb",
             "Fa18", "Sp18", "10.6_GeV"}},
        {"output/radiative_correction_plots", {"10.60", "10.2"}},
        {"output/bin_volume", {"10.60", "10.2"}},
        {"output/bin_centering_plots", {"10.60", "10.2"}},
        {"output/cross_sections",
            {"Fa18_Inb", "Fa18_Out", "Sp18_Inb", "Sp18_Out", "Sp19_Inb",
             "Fa18", "Sp18", "10.6_GeV"}}
    };

    for (const Category& category : categories) {
        std::vector<fs::path> created;
        created.reserve(1 + category.subs.size());

        const fs::path base(category.base);

        // create_directories() already checks whether the path exists. Avoiding
        // a separate exists() call removes one filesystem query per path.
        if (fs::create_directories(base)) {
            created.push_back(base);
        }

        for (const char* sub : category.subs) {
            const fs::path path = base / sub;
            if (fs::create_directories(path)) {
                created.push_back(path);
            }
        }

        // Preserve the previous grouped creation-message behavior.
        if (!created.empty()) {
            std::cout << "[Created] ";
            for (std::size_t i = 0; i < created.size(); ++i) {
                std::cout << created[i].string();
                if (i + 1 < created.size()) {
                    std::cout << ',';
                }
            }
            std::cout << '\n';
        }
    }
}
