#ifndef PYTHON_EXCLUSIVITY_RUNNER_H
#define PYTHON_EXCLUSIVITY_RUNNER_H

#include <string>

struct PythonExclusivityOptions {
    bool enabled = true;
    bool force_rerun = true;

    std::string python_executable = "python3";
    std::string script_path = "plot_exclusivity_data_dvcs_pi0_mc.py";
    std::string global_cuts_json = "output/jsons/global_cuts_config.json";
    std::string output_directory = "output/exclusivity_fit";
    std::string install_directory = "output/jsons";

    int workers = 7;

    double tight_containment = 0.90;
    double nominal_containment = 0.95;
    double loose_containment = 0.98;
};

// Run the production Python exclusivity fit, validate its three containment
// JSONs and atomically install the canonical files consumed downstream.
//
// Installed files:
//   combined_cuts_90.json
//   combined_cuts_95.json
//   combined_cuts_98.json
//   combined_cuts.json       (identical to the nominal 95% result)
//
// Returns true on success. All failures are reported to stderr.
bool run_python_exclusivity_analysis(
    const PythonExclusivityOptions& options);

#endif
