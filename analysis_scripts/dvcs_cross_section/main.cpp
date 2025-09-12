#include "make_dirs.h"
#include "load_trees.h"
#include "exclusivity_cuts.h"
#include "load_binning_scheme.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    // Create necessary output directories
    makeOutputDirs();
    std::cout << "Output directories ready." << std::endl;

    // --- Load binning scheme ---
    const std::string csv_file_path = "imports/integrated_bin_v2.csv";
    auto binning_scheme = load_binning_scheme(csv_file_path);
    std::cout << "Loaded binning scheme: " << binning_scheme.size() << " bins" << std::endl;

    // Containers for different tree categories
    std::map<std::string, TTree*> dataTrees;        // DVCS data
    std::map<std::string, TTree*> genMcTrees;       // DVCS generated MC (unused here)
    std::map<std::string, TTree*> recMcTrees;       // DVCS reconstructed MC
    std::map<std::string, TTree*> eppi0DataTrees;   // eppi0 data
    std::map<std::string, TTree*> eppi0GenMcTrees;  // eppi0 generated MC (unused here)
    std::map<std::string, TTree*> eppi0RecMcTrees;  // eppi0 reconstructed MC

    // Load all trees from files
    loadTrees(dataTrees, genMcTrees, recMcTrees, eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees);
    std::cout << "All trees loaded successfully." << std::endl;

    // Run exclusivity cut extraction (single-threaded for stability)
    runAllExclusivityCuts(
        dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
        "output/jsons", "output/exclusivity_plots", 1
    );

    std::cout << "Exclusivity-cut stage finished." << std::endl;
    return 0;
}