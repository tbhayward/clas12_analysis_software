#include "make_dirs.h"
#include "load_trees.h"
#include "exclusivity_cuts.h"
#include "load_binning_scheme.h"
#include "bin_means.h"
#include "total_counts.h"
#include "pi0_contamination.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    makeOutputDirs();
    std::cout << "Output directories ready." << std::endl;

    const std::string csv_file_path = "imports/integrated_bin_v2.csv";
    auto binning_scheme = load_binning_scheme(csv_file_path);
    std::cout << "Loaded binning scheme: " << binning_scheme.size() << " bins" << std::endl;

    std::map<std::string, TTree*> dataTrees;
    std::map<std::string, TTree*> genMcTrees;
    std::map<std::string, TTree*> recMcTrees;
    std::map<std::string, TTree*> eppi0DataTrees;
    std::map<std::string, TTree*> eppi0GenMcTrees;
    std::map<std::string, TTree*> eppi0RecMcTrees;
    // NOTE: load_trees currently loads eppi0_bkg into eppi0RecMcTrees by tag "*_bkg".
    loadTrees(dataTrees, genMcTrees, recMcTrees, eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees);
    std::cout << "All trees loaded successfully." << std::endl;

    // // runAllExclusivityCuts(...); // (kept commented during dev)

    std::vector<std::string> dvcs_periods = {
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "DVCS_Sp18_out", "DVCS_Sp18_inb"
        // Fa18_inb_supp handled by COPY step
    };
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)"};

    const std::string output_json_means = "output/jsons/bin_means_global.json";
    calculate_bin_means(dvcs_periods, topologies, "dvcs", binning_scheme, output_json_means, dataTrees);

    const std::string cuts_json_path   = "output/jsons/combined_cuts.json";
    const std::string output_counts_js = "output/jsons/total_counts.json";
    compute_total_counts(dvcs_periods, topologies, binning_scheme, dataTrees, cuts_json_path, output_counts_js);

    // --- NEW: helicity-resolved π0 contamination (no plotting yet) ---
    // eppi0 reco MC container holds both "*_rec_mc" and "*_bkg" based on your loader.
    // We pass it twice (once as reco MC, once as bkg), but lookup keys differ ("*_rec_mc" vs "*_bkg").
    compute_pi0_contamination_helicity(
        dvcs_periods,
        topologies,
        binning_scheme,
        dataTrees,
        eppi0DataTrees,
        eppi0RecMcTrees,   // reco MC (keys "*_rec_mc")
        eppi0RecMcTrees,   // bkg MC  (keys "*_bkg")
        cuts_json_path,
        "output/contamination"
    );

    std::cout << "All done." << std::endl;
    return 0;
}