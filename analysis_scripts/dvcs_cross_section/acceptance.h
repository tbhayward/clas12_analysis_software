#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <map>
#include <string>
#include <vector>
class TTree;

#include "load_binning_scheme.h" // Binning { xBmin,xBmax,Q2min,Q2max,tmin,tmax }

/**
 * Compute acceptance A(phi) = N_rec(phi) / N_gen(phi) per (xB, Q2, |t|) bin
 * for each period, using DVCS MC:
 *   - generated MC:   <runTag>_gen        (no MC exclusivity cuts)
 *   - reconstructed:  <runTag>_rec        (apply MC exclusivity cuts)
 *
 * Input:
 *   periods:   {"DVCS_Sp18_inb", "DVCS_Sp18_out", "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb"}
 *   topologies: vector of strings {"(FD,FD)","(CD,FD)","(CD,FT)"} to keep
 *   binning:   result of load_binning_scheme(...)
 *   genMcTrees: map tag->TTree*, with tags sp18_inb_gen, ...
 *   recMcTrees: map tag->TTree*, with tags sp18_inb_rec, ...
 *   cuts_json_path: path to output/jsons/combined_cuts.json (used for thresholds; falls back to defaults)
 *   out_root_dir:   "output"
 *
 * Output:
 *   - JSON per period: output/jsons/acceptance_<period>.json
 *   - Plots per period: output/acceptance/<runTag>/plot_acceptance_<period>_xB_<ix>.png
 */
void compute_and_plot_acceptance(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& cuts_json_path,
    const std::string& out_root_dir);

#endif // ACCEPTANCE_H