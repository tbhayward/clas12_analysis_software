#ifndef PYTHON_EXCLUSIVITY_RUNNER_H
#define PYTHON_EXCLUSIVITY_RUNNER_H

#include <string>

struct PythonExclusivityOptions {
    bool enabled = true;
    bool force_rerun = true;

    std::string python_executable = "python3";
    std::string script_path = "plot_exclusivity_data_dvcs_pi0_mc.py";
    std::string global_cuts_json = "output/jsons/global_cuts_config.json";

    // The Python stage writes its complete diagnostic output beneath this
    // directory. Its cut JSON files are expected in <output_directory>/jsons.
    std::string output_directory = "output/exclusivity_fit";

    // Validated cut JSON files are copied here for the downstream C++ stages.
    std::string install_directory = "output/jsons";

    int workers = 7;

    double tight_containment = 0.90;
    double nominal_containment = 0.95;
    double loose_containment = 0.98;
};

// Runs the production Python exclusivity analysis, validates its three
// containment JSON files, and installs the canonical downstream JSONs.
//
// Installed files:
//   combined_cuts_90.json
//   combined_cuts_95.json
//   combined_cuts_98.json
//   combined_cuts.json      (an exact copy of combined_cuts_95.json)
//
// Returns false immediately on any configuration, execution, validation, or
// installation failure.
bool run_python_exclusivity_analysis(
    const PythonExclusivityOptions& options = PythonExclusivityOptions{}
);

#endif
