// pi0_corrected_counts.h
#pragma once
#include <string>
#include <vector>

// Use the project's existing Binning struct
#include "load_binning_scheme.h"

/**
 * Compute DVCS counts corrected for π0 contamination, helicity-resolved and per φ-bin.
 *
 * Inputs:
 *  - periods:               e.g. {"DVCS_Fa18_inb", ...}
 *  - binning_scheme:        from load_binning_scheme(...) (only used for meta / validation)
 *  - total_counts_json:     path to "output/jsons/total_counts.json"
 *  - contamination_dir:     directory containing "contamination_<period>.json" files
 *  - out_root_dir:          pass "output" (writer handles subdirs)
 *
 * Outputs:
 *  - Per-period JSONs:      output/jsons/pi0_corrected_counts_<period>.json
 *  - Combined JSON:         output/jsons/pi0_corrected_counts_all_periods.json
 *
 * Notes:
 *  - Uncertainty propagation for N_corr = N * (1 - c) assumes
 *      Var(N_corr) = (1 - c)^2 * Var(N) + N^2 * Var(c),
 *    with Var(N)=N (Poisson) and Var(c) from contamination JSON.
 *  - Negative corrected counts are clamped at 0 (and σ left as computed).
 */
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& out_root_dir
);