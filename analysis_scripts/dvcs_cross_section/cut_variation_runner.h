#ifndef CUT_VARIATION_RUNNER_H
#define CUT_VARIATION_RUNNER_H

#include "global_cuts.h"
#include <map>
#include <string>

class TTree;

struct AutomaticCutVariationOptions {
    bool enabled = false;

    // Retained for source compatibility. The Python exclusivity stage controls
    // its own diagnostic output and does not use the legacy C++ extraction flag.
    bool make_exclusivity_extraction_plots = false;

    bool make_final_diagnostic_plots = true;
    bool use_pass1_tight_instability_rule = true;
    double tight_relative_difference_threshold = 0.50;
    int max_workers = 7;

    std::string nominal_csv = "output/csvs/dvcs_pass2_analysis.csv";
    std::string output_dir = "output/cut_variation_systematics";

    // Production Python exclusivity integration.
    std::string python_executable = "python3";
    std::string python_script_path =
        "plot_exclusivity_data_dvcs_pi0_mc.py";
    std::string production_cuts_dir = "output/jsons";

    double tight_containment = 0.90;
    double nominal_containment = 0.95;
    double loose_containment = 0.98;
};

bool run_automatic_cut_variation_systematics(
    const AutomaticCutVariationOptions& options,
    const GlobalCutConfig& nominal_global_cuts,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::map<std::string, TTree*>& currentStudyGenMcTrees,
    const std::map<std::string, TTree*>& currentStudyRecMcTrees,
    bool use_nobkg_dvcs_mc_for_acceptance,
    bool use_epg_mc_current_factor_for_eppi0_bkg);

#endif
