// setup_directories.cpp
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <string>
#include <vector>

void setup_directories() {
    std::vector<std::string> directories = {
        "output",
        "output/results",
        "output/binned_plots",
        "output/correlation_plots",
        "output/integrated_plots",
        "output/individual_chi2_fits",
        "output/structure_function",
        "output/misid_plots",
        "output/dependence_plots"
    };

    for (const auto& dir : directories) {
        struct stat info;
        if (stat(dir.c_str(), &info) != 0) {
            std::cout << "Directory " << dir << " does not exist. Creating it." << std::endl;
            if (mkdir(dir.c_str(), 0755) != 0) {
                std::cerr << "Error creating directory " << dir << std::endl;
            }
        } else if (!(info.st_mode & S_IFDIR)) {
            std::cerr << dir << " is not a directory!" << std::endl;
        }
    }
}