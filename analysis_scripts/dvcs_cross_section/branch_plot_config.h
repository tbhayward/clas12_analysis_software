#ifndef BRANCH_PLOT_CONFIG_H
#define BRANCH_PLOT_CONFIG_H

#include <map>
#include <string>

struct BranchHistConfig {
    int nbins;
    double xlow;
    double xhigh;
    bool enabled;
};

const std::map<std::string, BranchHistConfig>& getBranchPlotConfigs();
bool hasBranchPlotConfig(const std::string& branch_name);
const BranchHistConfig* findBranchPlotConfig(const std::string& branch_name);

#endif