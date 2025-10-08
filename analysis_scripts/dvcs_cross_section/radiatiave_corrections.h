#ifndef RADIATIVE_CORRECTIONS_H
#define RADIATIVE_CORRECTIONS_H

#include <map>
#include <string>
#include <vector>

class TTree;
struct Binning;

// Compute radiative corrections R_C per (xB, Q2, |t|, phi) bin for each period,
// using reconstructed MC with and without radiative effects, after applying
// MC 3σ exclusivity cuts and the standard global kinematic cuts + topology filters.
//
// Inputs:
//  - periods:     e.g. {"DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb", ...}
//  - topologies:  {"(FD,FD)","(CD,FD)","(CD,FT)"}
//  - binning:     from load_binning_scheme(...)
//  - recMcTrees:  non-radiative reconstructed DVCS MC, keys: "sp18_inb_rec", "fa18_inb_rec", ...
//  - radRecMcTrees: radiative reconstructed DVCS MC, keys: "sp18_inb_rec_rad", ...
//  - cuts_json_path: "output/jsons/combined_cuts.json" (produced by exclusivity_cuts)
//  - out_root_dir:   "output"
//
// Outputs:
//  - JSON per period: output/jsons/radcorr_<period>.json
//  - Combined summary: output/jsons/radcorr_all_periods.json
//  - Plots per xB-slice: output/radiative_correction_plots/<runTag>/plot_radcorr_<period>_xB_<ix>.png
void compute_radiative_corrections(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::map<std::string, TTree*>& radRecMcTrees,
    const std::string& cuts_json_path,
    const std::string& out_root_dir
);

#endif // RADIATIVE_CORRECTIONS_H