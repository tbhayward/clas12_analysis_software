#ifndef BRANCH_DATA_MC_COMPARISON_H
#define BRANCH_DATA_MC_COMPARISON_H

#include <map>
#include <string>

class TTree;

bool runAllBranchDataMcComparisons(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& outPlotDir,
    int max_workers = 5
);

#endif // BRANCH_DATA_MC_COMPARISON_H
