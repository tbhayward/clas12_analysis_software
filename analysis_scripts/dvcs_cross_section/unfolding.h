#ifndef UNFOLDING_H
#define UNFOLDING_H

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "load_binning_scheme.h" // struct Binning { double xBmin,xBmax,Q2min,Q2max,tmin,tmax; }

class TTree;

// Compute unfolded yields N/A per φ-bin for each kinematic bin, helicity-resolved.
// Inputs:
//  - periods:                     e.g. {"DVCS_Fa18_inb", ...}  (skip Fa18_inb_supp in your caller)
//  - binning_scheme:              result of load_binning_scheme(...)
//  - total_counts_json_path:      "output/jsons/total_counts.json"
//  - out_root_dir:                "output"  (acceptance JSONs expected at output/jsons/acceptance_<period>.json)
//
// Writes per-period JSONs to:     output/jsons/unfolded_<period>.json
// Makes plots to:                 output/unfolding/<runTag>/plot_unfolded_<period>_xB_<ix>.png
void compute_and_plot_unfolding(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json_path,
    const std::string& out_root_dir);



#endif // UNFOLDING_H