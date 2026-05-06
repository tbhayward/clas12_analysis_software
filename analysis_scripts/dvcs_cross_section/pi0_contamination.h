#ifndef PI0_CONTAMINATION_H
#define PI0_CONTAMINATION_H

#include <map>
#include <string>

class TTree;

// Computes the pi0 contamination ratio using only quantities already stored in
// dvcs_pass2_analysis.csv. The TTree maps are retained in the signature for
// compatibility with the existing main.cpp call, but this implementation does
// not loop over ROOT trees.
//
// For each CSV row and period:
//
//   c = N_mis * N_pi0_data / (N_pi0_rec_mc * N_dvcs_data)
//
// where the terms are read from existing CSV columns:
//
//   N_dvcs_data   = sum_topo normalized raw yield, ep->epg,        topo, exp, period, unpol
//   N_pi0_data    = sum_topo normalized raw yield, ep->eppi0,      topo, exp, period, unpol
//   N_pi0_rec_mc  = sum_topo reconstructed current corrected yield, ep->eppi0,      topo, mc, period
//   N_mis         = sum_topo reconstructed current corrected yield, ep->eppi0->epg, topo, mc, period
//
// The output is written to:
//
//   contamination ratio, <period>
//
// as a triplet:
//
//   (value,stat,sys)
//
// with sys currently set to 0.
bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*> &dvcsDataTrees,
    const std::map<std::string, TTree*> &eppi0DataTrees,
    const std::map<std::string, TTree*> &eppi0RecMcTrees,
    const std::map<std::string, TTree*> &eppi0BkgTrees,
    const std::string &combined_cuts_json,
    const std::string &csv_main,
    const std::string &output_root_dir,
    int max_workers);

#endif // PI0_CONTAMINATION_H
