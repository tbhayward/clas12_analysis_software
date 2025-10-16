#ifndef BSA_H
#define BSA_H

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "load_binning_scheme.h" // provides struct Binning

class TTree;

// Compute helicity-resolved BSA including per-bin polarization and π0 contamination,
// write JSONs, and plot per-period grids plus a 10.6 GeV combined set.
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& total_counts_json_path,
    const std::string& contamination_dir,
    const std::string& out_root_dir);

#endif // BSA_H