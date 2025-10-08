#ifndef ACCEPTANCE_H
#define ACCEPTANCE_H

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "load_binning_scheme.h" // struct Binning { double xBmin,xBmax,Q2min,Q2max,tmin,tmax; }

class TTree;

// Compute acceptance A(φ) = N_reco_pass / N_gen in each (xB,Q²,|t|,φ) bin,
// per run period. Reconstructed MC uses the *MC-side* 3σ exclusivity cuts
// from output/jsons/combined_cuts.json; generated MC is uncut.
//
// Uncertainty: **binomial**
//   sigma(A) = sqrt( A * (1 - A) / Ngen )
//
// Inputs:
//  - periods:        {"DVCS_Sp18_inb","DVCS_Sp18_out","DVCS_Fa18_inb","DVCS_Fa18_out","DVCS_Sp19_inb"}
//                    (skip supplemental externally)
//  - topologies:     {"(FD,FD)","(CD,FD)","(CD,FT)"}  (applied to reconstructed MC)
//  - binning:        from load_binning_scheme(...)
//  - genMcTrees:     keys "<tag>_gen"
//  - recMcTrees:     keys "<tag>_rec"
//  - cuts_json_path: "output/jsons/combined_cuts.json"
//  - out_root_dir:   "output"
void compute_and_plot_acceptance(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& cuts_json_path,
    const std::string& out_root_dir);

#endif // ACCEPTANCE_H