#pragma once

#include <map>
#include <string>
#include <vector>
#include <TTree.h>

#include "load_binning_scheme.h"  

// Computes and plots BSA using:
// - total_counts.json (helicity-resolved counts by bin/group)
// - contamination_<period>.json files (helicity-resolved pi0 contamination)
// - beam_pol from DVCS trees for polarization scaling
// Writes per-period fit JSONs under output/jsons/BSA_fits/,
// an all-periods JSON under output/jsons/BSA_fits_all_periods.json,
// and the 10.6 GeV combined JSON under output/jsons/BSA_fits_combined_10p6.json.
// Also writes plots under output/bsa_plots/<runTag>/...
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& total_counts_json_path,
    const std::string& contamination_jsons_dir,
    const std::string& out_root_dir
);