#ifndef BRANCH_DATA_MC_COMPARISON_H
#define BRANCH_DATA_MC_COMPARISON_H

#include <map>
#include <string>

class TTree;

void runAllBranchDataMcComparisons(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& outPlotDir
);

#endif