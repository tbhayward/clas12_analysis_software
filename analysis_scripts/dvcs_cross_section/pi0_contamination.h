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
//   contamination_<period>.json per DVCS period (and optional copy to Fa18_inb_supp).
// NOTE: Plotting (ROOT canvases) is produced automatically at the end if
//       ENABLE_PI0_CONTAMINATION_PLOTS is true.
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,                 // e.g. {"DVCS_Fa18_inb", ...}
    const std::vector<std::string>& topologies,              // {"(FD,FD)","(CD,FD)","(CD,FT)"}
    const std::vector<Binning>& binning_scheme,              // from load_binning_scheme(...)
    const std::map<std::string, TTree*>& dvcsDataTrees,      // keys: "fa18_inb", "sp18_out", ...
    const std::map<std::string, TTree*>& eppi0DataTrees,     // keys: "<tag>_eppi0"
    const std::map<std::string, TTree*>& eppi0RecMcTrees,    // keys: "<tag>_rec_mc" (and "<tag>_bkg" are also in here)
    const std::map<std::string, TTree*>& eppi0BkgTrees,      // keys: "<tag>_bkg" (same map as above is fine)
    const std::string& combined_cuts_json,                   // "output/jsons/combined_cuts.json"
    const std::string& out_dir                               // e.g. "output/contamination"
);

// -------------------- Optional: plot-from-JSON helper --------------------
// If you’d like to (re)generate canvases later from saved JSON files (without recomputing),
// you can call this function separately (not used by default).
void plot_pi0_contamination_from_json(
    const std::string& period,                               // "DVCS_Fa18_inb"
    const std::vector<Binning>& binning_scheme,
    const std::string& contamination_json_path,              // "output/contamination/contamination_<period>.json"
    const std::string& out_dir_plots                         // "output/contamination/plots"
);

#endif // PI0_CONTAMINATION_H