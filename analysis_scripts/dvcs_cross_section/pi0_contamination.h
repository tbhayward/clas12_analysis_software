#ifndef PI0_CONTAMINATION_H
#define PI0_CONTAMINATION_H

#include <TTree.h>
#include <map>
#include <string>
#include <vector>

#include "load_binning_scheme.h" // Binning struct

// Toggle: also emit a contamination JSON for Fa18_inb_supp by copying Fa18_inb
// Default true (as requested)
extern bool COPY_CONTAM_TO_FA18_INB_SUPP;

// Compute helicity-resolved pi0 contamination for the requested DVCS periods.
// - periods: e.g. {"DVCS_Fa18_inb", "DVCS_Fa18_out", ...}
// - topologies: {"(FD,FD)","(CD,FD)","(CD,FT)"}; an event is accepted if it matches any
// - binning_scheme: integrated_bin_v2.csv already loaded
// - trees: all categories already loaded via load_trees.cpp
//   * DVCS data trees accessed with keys like "fa18_inb", "sp18_out", etc. (periodToRunTagKey)
//   * eppi0 data trees accessed with keys like "fa18_inb_eppi0", ...
//   * eppi0 background MC (mis-ID to DVCS) trees accessed with keys like "fa18_inb_bkg"
//   * eppi0 reco MC (for N_pi0_reco) with keys like "fa18_inb_rec_mc"
// - combined_cuts_json: output/jsons/combined_cuts.json (from exclusivity_cuts)
// - bin_means_json: output/jsons/bin_means_global.json (for future plotting; not required here)
// - out_dir: directory where contamination_*.json will be written (e.g. "output/contamination")
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,   // reco MC for π0
    const std::map<std::string, TTree*>& eppi0BkgTrees,     // mis-ID background MC to DVCS
    const std::string& combined_cuts_json,
    const std::string& out_dir
);

#endif // PI0_CONTAMINATION_H