#ifndef CALCULATE_CONTAMINATION_H
#define CALCULATE_CONTAMINATION_H

#include <string>
#include <vector>
#include <map>
#include "unfolding_data.h"
#include "bin_boundaries.h"
#include <TTreeReader.h>

// Function declaration
void calculate_contamination(
    const std::string& base_output_dir,
    int xB_bin,
    const std::vector<BinBoundary>& bin_boundaries,
    std::vector<TTreeReader>& data_readers,
    std::vector<TTreeReader>& eppi0_readers,
    std::vector<TTreeReader>& mc_rec_aaogen_readers,
    std::vector<TTreeReader>& mc_rec_eppi0_bkg_readers,
    std::map<std::string, std::vector<UnfoldingData>>& bin_data
);

#endif // CALCULATE_CONTAMINATION_H