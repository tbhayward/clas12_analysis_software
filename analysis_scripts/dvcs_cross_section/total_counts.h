#ifndef TOTAL_COUNTS_H
#define TOTAL_COUNTS_H

#include <map>
#include <string>
#include <vector>
#include <TTree.h>
#include "load_binning_scheme.h"

// Updated declaration with 7 parameters
void compute_total_counts(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_json_path,
    const std::string& out_root_dir   // e.g. "output"
);

#endif