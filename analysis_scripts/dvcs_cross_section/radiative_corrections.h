// radiative_corrections.h
#ifndef RADIATIVE_CORRECTIONS_H
#define RADIATIVE_CORRECTIONS_H

#include <map>
#include <string>
#include <vector>

class TTree;

// Use the canonical Binning definition
#include "load_binning_scheme.h"

// Compute radiative corrections using GENERATED trees only:
//   RC(phi) = (N_gen_rad(cell,phi)/N_gen_rad_total) / (N_gen_born(cell,phi)/N_gen_born_total)
//
// Notes:
// - 'topologies' kept for API compatibility (ignored in implementation).
// - recMcTrees_* and combined_cuts_json_path are also ignored.
void compute_radiative_corrections(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees_norad,
    const std::map<std::string, TTree*>& recMcTrees_norad,
    const std::map<std::string, TTree*>& genMcTrees_rad,
    const std::map<std::string, TTree*>& recMcTrees_rad,
    const std::string& combined_cuts_json_path,
    const std::string& out_root_dir
);

#endif // RADIATIVE_CORRECTIONS_H