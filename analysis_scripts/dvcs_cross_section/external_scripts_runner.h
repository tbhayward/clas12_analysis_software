#ifndef EXTERNAL_SCRIPTS_RUNNER_H
#define EXTERNAL_SCRIPTS_RUNNER_H

#include <string>

struct ExternalScriptOptions {
    std::string scripts_directory = "external_scripts";
    std::string python_executable = "python";

    bool include_bin_to_bin_systematics = true;
    bool use_simple_clas6_cross_check = true;
};

// Runs the external Python cross-section studies sequentially.
// Returns false immediately if any script cannot be launched or exits nonzero.
bool run_external_cross_section_scripts(
    const std::string& analysis_csv,
    const ExternalScriptOptions& options = ExternalScriptOptions{}
);

#endif
