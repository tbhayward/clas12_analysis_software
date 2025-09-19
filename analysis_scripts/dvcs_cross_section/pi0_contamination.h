#ifndef PI0_CONTAMINATION_H
#define PI0_CONTAMINATION_H

#include <TTree.h>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "load_binning_scheme.h" // Binning

// -------------------- User toggles --------------------
// Copy Fa18_inb contamination JSON to Fa18_inb_supp (default: true)
extern bool COPY_CONTAM_TO_FA18_INB_SUPP;

// Make ROOT plots for contamination after computing (default: true)
extern bool ENABLE_PI0_CONTAMINATION_PLOTS;

// -------------------- Public API --------------------
// Compute helicity-resolved pi0 contamination and write JSONs:
//   - Per-period files: output/jsons/contamination/contamination_<period>.json
//   - Combined file:    output/jsons/pi0_contamination_combined.json
// Plots (if enabled) go to: output/contamination_plots/<runTag>/
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,                 // e.g. {"DVCS_Fa18_inb", ...}
    const std::vector<std::string>& topologies,              // {"(FD,FD)","(CD,FD)","(CD,FT)"}
    const std::vector<Binning>& binning_scheme,              // from load_binning_scheme(...)
    const std::map<std::string, TTree*>& dvcsDataTrees,      // keys: "fa18_inb", "sp18_out", ...
    const std::map<std::string, TTree*>& eppi0DataTrees,     // keys: "<tag>_eppi0"
    const std::map<std::string, TTree*>& eppi0RecMcTrees,    // keys: "<tag>_rec_mc"
    const std::map<std::string, TTree*>& eppi0BkgTrees,      // keys: "<tag>_bkg"
    const std::string& combined_cuts_json,                   // "output/jsons/combined_cuts.json"
    const std::string& out_root_dir                          // pass "output" (writer handles subdirs)
);

// -------------------- Optional: plot-from-JSON helper --------------------
// Re-generate canvases later from a saved JSON without recomputing.
void plot_pi0_contamination_from_json(
    const std::string& period,                               // "DVCS_Fa18_inb"
    const std::vector<Binning>& binning_scheme,
    const std::string& contamination_json_path,              // "output/jsons/contamination/contamination_<period>.json"
    const std::string& out_dir_plots                         // "output/contamination_plots/<runTag>"
);

#endif // PI0_CONTAMINATION_H