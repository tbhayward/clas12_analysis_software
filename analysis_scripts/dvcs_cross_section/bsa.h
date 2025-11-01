#ifndef BSA_H
#define BSA_H

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "load_binning_scheme.h" // provides struct Binning

class TTree;

/**
 * Compute helicity-resolved Beam-Spin Asymmetry (BSA)
 * using pi0-corrected helicity counts instead of raw total_counts.
 *
 * Inputs:
 *   - pi0_corrected_counts_json_path : JSON file containing all groups,
 *     e.g. "output/jsons/pi0_corrected_counts_all_groups.json"
 *   - dvcsDataTrees                  : map from runTag to DVCS TTree* (for beam polarization)
 *   - binning_scheme                 : xB, Q2, t bin definitions
 *
 * Outputs:
 *   - JSONs:   output/jsons/BSA_fits/BSA_fits_<period>.json
 *              output/jsons/BSA_fits_all_periods.json
 *              output/jsons/BSA_fits_combined_10.6.json
 *   - Plots:   output/bsa_plots/<runTag>/plot_bsa_<period>_xB_<ix>.png
 *              output/bsa_plots/10.6_combined/...
 */
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& pi0_corrected_counts_json_path,
    const std::string& out_root_dir);

#endif // BSA_H