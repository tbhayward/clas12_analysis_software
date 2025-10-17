#pragma once
#include <string>
#include <vector>
#include "load_binning_scheme.h"

/**
 * Compute DVCS counts corrected for π0 contamination, helicity-resolved and per φ-bin,
 * with uncertainty propagation:
 *   N_corr,h = N_h * (1 - c_h)
 *   Var(N_corr,h) = (1 - c_h)^2 * Var(N_h) + (N_h)^2 * Var(c_h) with Var(N_h)=N_h (Poisson).
 *
 * Inputs:
 *  - periods:               e.g. {"DVCS_Fa18_inb", ...}
 *  - binning_scheme:        from load_binning_scheme(...)
 *  - total_counts_json:     path to output/jsons/total_counts.json (has "groups", including combined)
 *  - contamination_dir:     directory with per-period files contamination_<period>.json
 *  - contamination_combined: file output/jsons/pi0_contamination_combined.json (has "periods" incl. combined groups)
 *  - out_root_dir:          pass "output"
 *
 * Writes:
 *  - Per-group JSONs:       output/jsons/pi0_corrected_counts_<GROUP>.json
 *                           (GROUP is each period plus Spring2018, Fall2018, 10.6_GeV if present)
 *  - Master combined JSON:  output/jsons/pi0_corrected_counts_all_groups.json
 *  - Plots:                 output/pi0_corrected_plots/<GROUP>/plot_pi0corr_<GROUP>_xB_<ix>.png
 */
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& contamination_combined,
    const std::string& out_root_dir
);