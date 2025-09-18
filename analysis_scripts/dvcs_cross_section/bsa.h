#pragma once
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

class TTree;

struct Binning {
    double xBmin, xBmax;
    double Q2min, Q2max;
    double tmin,  tmax;
};

void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,                 // e.g. {"DVCS_Sp18_inb", ...}
    const std::vector<std::string>& topologies,              // {"(FD,FD)", "(CD,FD)", "(CD,FT)"}
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,      // keyed by runTag: "sp18_inb", ...
    const std::string& total_counts_json_path,               // output/jsons/total_counts/total_counts.json (from your total_counts.cpp)
    const std::string& contamination_dir,                    // output/jsons/contamination/
    const std::string& out_root_dir                          // usually "output"
);