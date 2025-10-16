#pragma once
#include <map>
#include <string>
#include <vector>
class TTree;

#include "load_binning_scheme.h" // Binning

void compute_total_counts(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,  // "output/jsons/combined_cuts.json"
    const std::string& out_json_path        // "output/jsons/total_counts.json"
);