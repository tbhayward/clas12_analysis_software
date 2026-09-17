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

    // Detector-normalization mode must follow the nominal production run.
    // When eppi0 is enabled, each cut variation re-derives its own sequential
    // theta_p * p_p efficiency map using that variation's global/fiducial and
    // exclusivity cuts, then applies that variation-specific map to reconstructed
    // DVCS MC.  When disabled, variations restore the Krishna proton-efficiency
    // correction exactly as the nominal --no-eppi0-normalization mode does.
    bool use_eppi0_production_normalization = true;
    std::string eppi0_charge_csv_path = "imports/integrated_luminosity/global.csv";
    std::string eppi0_normalization_json_path = "imports/eppi0_aao_normalization_inputs.json";
    std::string current_response_model_json =
        "output/dvcs_current_dependence/calibration/current_response_model.json";

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
