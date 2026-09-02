#ifndef EXTERNAL_SCRIPTS_RUNNER_H
#define EXTERNAL_SCRIPTS_RUNNER_H

#include <string>

struct ExternalScriptOptions {
    std::string scripts_directory = "external_scripts";
    std::string python_executable = "python";

    // Authoritative published CLAS12 pass-1 cross-section table.  Legacy
    // all_bin_v3.csv remains the bin-geometry/correction source elsewhere.
    std::string published_pass1_cross_section_table =
        "imports/clasdb_E214M1.txt";

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
