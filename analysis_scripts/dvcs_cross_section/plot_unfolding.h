// plot_unfolding.h
#ifndef PLOT_UNFOLDING_H
#define PLOT_UNFOLDING_H

#include <TH1D.h>
#include <TCanvas.h>
#include <TTreeReader.h>
#include <cmath>
#include <string>
#include <vector>

#include "bin_boundaries.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include "unfolding_data.h"  // Include the new header

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Function prototype for plot_unfolding
std::vector<UnfoldingData> plot_unfolding(const std::string& output_dir, 
                                          const std::string& analysisType, 
                                          int xB_bin,
                                          const std::vector<BinBoundary>& bin_boundaries, 
                                          std::vector<TTreeReader>& data_readers,  
                                          std::vector<TTreeReader>& mc_gen_readers,  
                                          std::vector<TTreeReader>& mc_rec_readers, 
                                          std::vector<TTreeReader>& eppi0_readers,  
                                          std::vector<TTreeReader>& mc_gen_aaogen_readers,  
                                          std::vector<TTreeReader>& mc_rec_aaogen_readers);

#endif // PLOT_UNFOLDING_H