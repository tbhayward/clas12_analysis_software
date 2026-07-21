#ifndef SYSTEMATICS_RUNNER_H
#define SYSTEMATICS_RUNNER_H

#include <string>

struct SystematicsRunnerOptions {
    std::string executable = "./main_systematics";
    std::string pass1_systematics_csv = "imports/pass1_systematic_summary.csv";
};

// Run the standalone CSV-only systematics executable and wait for it to finish.
// Returns true only when the executable exits normally with status zero.
bool run_main_systematics(
    const std::string& analysis_csv,
    const SystematicsRunnerOptions& options = SystematicsRunnerOptions{}
);

#endif
