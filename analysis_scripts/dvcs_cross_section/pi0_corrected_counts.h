// pi0_contamination.h
#pragma once

#include <string>
#include <vector>

struct Binning;
class TTree;

// contamination APIs (already implemented in pi0_contamination.cpp)
void compute_pi0_contamination_helicity(
  const std::vector<std::string>& periods,
  const std::vector<std::string>& topologies,
  const std::vector<Binning>& binning_scheme,
  const std::map<std::string, TTree*>& dvcsDataTrees,
  const std::map<std::string, TTree*>& eppi0DataTrees,
  const std::map<std::string, TTree*>& eppi0RecMcTrees,
  const std::map<std::string, TTree*>& eppi0BkgTrees,
  const std::string& combined_cuts_json,
  const std::string& out_root_dir
);

void plot_pi0_contamination_from_json(
  const std::string& period,
  const std::vector<Binning>& binning_scheme,
  const std::string& contamination_json_path,
  const std::string& out_dir_plots
);

// single (extern) definitions — defined only once in pi0_contamination.cpp
extern bool COPY_CONTAM_TO_FA18_INB_SUPP;
extern bool ENABLE_PI0_CONTAMINATION_PLOTS;
extern bool VERBOSE_CONTAM_DEBUG;