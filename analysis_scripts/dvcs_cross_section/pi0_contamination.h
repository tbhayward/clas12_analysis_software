#pragma once
#include <map>
#include <string>
#include <vector>
class TTree;
#include "load_binning_scheme.h"

// Compute helicity-resolved π0 contamination and write JSONs/plots.
// NOTE: requires both reconstructed eπ0 MC and eπ0 background MC maps.
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,  // <-- ensure this param exists
    const std::string& combined_cuts_json_path,
    const std::string& out_root_dir);