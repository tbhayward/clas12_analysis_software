#ifndef PI0_CONTAMINATION_H
#define PI0_CONTAMINATION_H

#include <map>
#include <string>

class TTree;

// Computes per-bin pi0 contamination ratio:
//   c_i = Npi0mis * (Npi0_exp / Npi0_sim) / NDVCS_exp
//
// Where each N is counted in the same (xB,Q2,|t|,phi) binning (CSV rows).
// This implementation SUMS over topologies (FD_FD, CD_FD, CD_FT) because the
// existing CSV schema provides only one column per period:
//
//   "contamination ratio, <PeriodDisplay>"
//
// Cut-block selection is explicit:
//   - NDVCS_exp and Npi0mis use DVCS_* cut blocks
//   - Npi0_exp and Npi0sim use eppi0_* cut blocks
//
bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*> &dvcsDataTrees,
    const std::map<std::string, TTree*> &eppi0DataTrees,
    const std::map<std::string, TTree*> &eppi0RecMcTrees,
    const std::map<std::string, TTree*> &eppi0BkgTrees,
    const std::string &combined_cuts_json,
    const std::string &csv_main,
    const std::string &output_root_dir,
    int max_workers);

#endif