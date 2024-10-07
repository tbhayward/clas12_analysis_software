#ifndef PLOT_DVCS_DATA_MC_COMPARISON_H
#define PLOT_DVCS_DATA_MC_COMPARISON_H

#include <TTreeReader.h>
#include <string>

// Function to plot the normalized data and MC comparisons for DVCS
void plot_dvcs_data_mc_comparison(TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader, 
                                  const std::string& output_dir, int xB_bin);

#endif