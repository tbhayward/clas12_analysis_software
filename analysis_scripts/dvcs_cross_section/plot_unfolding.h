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

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Function prototype for plot_unfolding
void plot_unfolding(const std::string& output_dir, 
                    const std::string& analysisType, 
                    int xB_bin,
                    const std::vector<BinBoundary>& bin_boundaries, 
                    const std::vector<TTreeReader>& data_readers,  // Changed to vector of TTreeReader
                    const std::vector<TTreeReader>& mc_gen_readers,  // Changed to vector of TTreeReader
                    const std::vector<TTreeReader>& mc_rec_readers);  // Changed to vector of TTreeReader

#endif // PLOT_UNFOLDING_H