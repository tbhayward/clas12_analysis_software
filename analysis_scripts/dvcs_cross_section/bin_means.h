#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include "load_binning_scheme.h"
#include <TTree.h>

#include <map>
#include <string>
#include <tuple>
#include <vector>

// Compute GLOBAL (combined across periods/topologies) bin-averaged kinematics.
// - dvcs_periods: e.g. {"DVCS_Fa18_inb","DVCS_Fa18_out","DVCS_Sp19_inb","DVCS_Sp18_out","DVCS_Fa18_inb_supp"}
// - topologies:   e.g. {"(FD,FD)","(CD,FD)","(CD,FT)"}
// - analysis_type: currently unused (kept for parity with Python)
// - binning_scheme: vector of Binning{xBmin,xBmax,Q2min,Q2max,tmin,tmax}
// - output_json: path (e.g. "output/jsons/bin_means_global.json")
// - dvcsDataTrees: map<tag,TTree*> where tags are like "fa18_inb","fa18_out","sp19_inb","sp18_out","fa18_inb_supp"
void calculate_bin_means(const std::vector<std::string>& dvcs_periods,
                         const std::vector<std::string>& topologies,
                         const std::string& analysis_type,
                         const std::vector<Binning>& binning_scheme,
                         const std::string& output_json,
                         const std::map<std::string, TTree*>& dvcsDataTrees);

#endif // BIN_MEANS_H