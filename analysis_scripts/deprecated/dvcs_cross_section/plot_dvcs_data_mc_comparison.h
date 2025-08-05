#ifndef PLOT_DVCS_DATA_MC_COMPARISON_H
#define PLOT_DVCS_DATA_MC_COMPARISON_H

#include <string>
#include <vector>
#include <TTreeReader.h>
#include "bin_boundaries.h"

void plot_dvcs_data_mc_comparison(const std::string& output_dir, 
                                  const std::string& analysisType, 
                                  const std::string& dataset, 
                                  int xB_bin,
                                  const std::vector<BinBoundary>& bin_boundaries, 
                                  TTreeReader& data_reader, 
                                  TTreeReader& mc_gen_reader, 
                                  TTreeReader& mc_rec_reader);

#endif