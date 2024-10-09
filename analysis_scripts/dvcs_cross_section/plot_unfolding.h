#ifndef PLOT_UNFOLDING_H
#define PLOT_UNFOLDING_H

#include <string>
#include <vector>
#include "bin_boundaries.h"
#include <TTreeReader.h>

// Function to plot unfolded distributions and acceptance
void plot_unfolding(const std::string& output_dir, 
                    const std::string& analysisType, 
                    int xB_bin,
                    const std::vector<BinBoundary>& bin_boundaries, 
                    TTreeReader& data_reader);

#endif  // PLOT_UNFOLDING_H