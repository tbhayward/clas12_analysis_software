#ifndef CUT_VARIATION_RUNNER_H
#define CUT_VARIATION_RUNNER_H

#include "global_cuts.h"
#include <map>
#include <string>

class TTree;

struct AutomaticCutVariationOptions {
    bool enabled = false;
    bool make_exclusivity_extraction_plots = false;
    bool make_final_diagnostic_plots = true;
    int max_workers = 7;
    std::string nominal_csv = "output/csvs/dvcs_pass2_analysis.csv";
    std::string output_dir = "output/cut_variation_systematics";
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
