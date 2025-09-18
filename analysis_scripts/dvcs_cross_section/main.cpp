#include "make_dirs.h"
#include "load_trees.h"
#include "exclusivity_cuts.h"
#include "load_binning_scheme.h"
#include "bin_means.h"
#include "total_counts.h"
#include "pi0_contamination.h"   // <-- needed for compute_pi0_contamination_helicity
#include "bsa.h"                 // <-- needed for compute_and_plot_bsa_helicity

#include <filesystem>            // <-- needed for std::filesystem::path
#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    // Create necessary output directories
    makeOutputDirs();
    std::cout << "Output directories ready." << std::endl;

    // Root of output tree (used by several stages)
    const std::string output_root = "output";

    // --- Load binning scheme ---
    const std::string csv_file_path = "imports/integrated_bin_v2.csv";
    auto binning_scheme = load_binning_scheme(csv_file_path);
    std::cout << "Loaded binning scheme: " << binning_scheme.size() << " bins" << std::endl;

    // Containers for different tree categories
    std::map<std::string, TTree*> dataTrees;        // DVCS data
    std::map<std::string, TTree*> genMcTrees;       // DVCS generated MC 
    std::map<std::string, TTree*> recMcTrees;       // DVCS reconstructed MC
    std::map<std::string, TTree*> eppi0DataTrees;   // eppi0 data
    std::map<std::string, TTree*> eppi0GenMcTrees;  // eppi0 generated MC 
    std::map<std::string, TTree*> eppi0RecMcTrees;  // eppi0 reconstructed MC

    // Load all trees from files
    loadTrees(dataTrees, genMcTrees, recMcTrees, eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees);
    std::cout << "All trees loaded successfully." << std::endl;

    // // Run exclusivity cut extraction (single-threaded for stability)
    // runAllExclusivityCuts(
    //     dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
    //     "output/jsons", "output/exclusivity_plots", 1
    // );
    // std::cout << "Exclusivity-cut stage finished." << std::endl;

    // --------- Global bin-averaged kinematics ----------
    std::vector<std::string> dvcs_periods = {
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "DVCS_Sp18_out", "DVCS_Sp18_inb", "DVCS_Fa18_inb_supp"
    };
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)"};
    const std::string analysis_type = "dvcs";
    const std::string output_json_means = "output/jsons/bin_means_global.json";

    calculate_bin_means(dvcs_periods, topologies, analysis_type, binning_scheme, output_json_means, 
        dataTrees);

    // --------- Total counts after exclusivity cuts (by helicity) ----------
    const std::string cuts_json_path   = "output/jsons/combined_cuts.json"; 
    // produced by exclusivity_cuts
    const std::string output_counts_js = "output/jsons/total_counts.json";

    compute_total_counts(dvcs_periods, topologies, binning_scheme, dataTrees, cuts_json_path, 
        output_counts_js);

    // Helicity-resolved π0 contamination
    // NOTE: pass the OUTPUT ROOT ("output") so the implementation writes:
    //   - per-period JSONs to output/jsons/contamination/
    //   - combined JSON to output/jsons/
    //   - plots to output/contamination_plots/...
    compute_pi0_contamination_helicity(
        dvcs_periods,
        topologies,
        binning_scheme,
        dataTrees,
        eppi0DataTrees,
        eppi0RecMcTrees,   // reco MC (keys "*_rec_mc")
        eppi0RecMcTrees,   // bkg MC  (keys "*_bkg")
        cuts_json_path,
        output_root        // <<< was "output/contamination"; must be the root "output"
    );

    // Beam-Spin Asymmetry: reads total_counts.json and contamination JSONs,
    // writes per-period BSA fits to output/jsons/BSA_fits/BSA_fits_<period>.json,
    // writes all-periods file to output/jsons/BSA_fits_all_periods.json,
    // writes 10.6 GeV combined to output/jsons/BSA_fits_combined_10p6.json,
    // and plots to output/bsa_plots/<runTag>/...
    namespace fs = std::filesystem;
    const std::string contamination_dir = (fs::path(output_root) / "jsons" / "contamination").string();
    compute_and_plot_bsa_helicity(
        dvcs_periods,
        topologies,
        binning_scheme,
        dataTrees,            // DVCS trees (for beam_pol extraction)
        output_counts_js,     // total_counts.json path
        contamination_dir,    // directory with contamination_<period>.json files
        output_root
    );

    std::cout << "All done." << std::endl;
    return 0;
}