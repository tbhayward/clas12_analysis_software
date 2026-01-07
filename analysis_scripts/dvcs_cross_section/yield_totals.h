// yield_totals.h
#ifndef YIELD_TOTALS_H
#define YIELD_TOTALS_H

#include <map>
#include <string>

class TTree;

// Compute total reconstructed yields after exclusivity cuts.
//  - dvcsDataTrees: keys like "DVCS_Fa18_inb", "DVCS_Sp18_out", etc.
//  - dvcsRecMcTrees: keys like "DVCS_Fa18_inb_rec", ...
//  - combined_cuts_json: path to combined_cuts.json produced by runAllExclusivityCuts
//  - output_txt: path to a text file where a summary will be written.
//
// Returns true on success, false on any fatal problem (missing cuts, etc).
bool compute_yield_totals(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& output_txt);

#endif // YIELD_TOTALS_H