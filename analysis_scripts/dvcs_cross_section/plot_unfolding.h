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
#include "unfolding_data.h"  // Include the UnfoldingData structure

// Function prototype for plot_unfolding
std::map<std::string, std::vector<UnfoldingData>> plot_unfolding(
    const std::string& output_dir,
    int xB_bin,
    const std::vector<BinBoundary>& bin_boundaries,
    std::vector<TTreeReader>& data_readers,
    std::vector<TTreeReader>& mc_gen_readers,
    std::vector<TTreeReader>& mc_rec_readers,
    std::vector<TTreeReader>& eppi0_readers,
    std::vector<TTreeReader>& mc_gen_aaogen_readers,
    std::vector<TTreeReader>& mc_rec_aaogen_readers
);

#endif // PLOT_UNFOLDING_H