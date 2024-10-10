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

// Struct to hold the unfolding data for each phi bin within an xB-Q2-t bin
struct UnfoldingData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;

    // Phi bin data (min and max for each bin)
    std::vector<double> phi_min;  // Lower boundary for each phi bin
    std::vector<double> phi_max;  // Upper boundary for each phi bin

    // 2D Vectors to store raw yields, acceptance, and unfolded yields for each period and topology
    std::vector<int> raw_yields_Fa18Inb;  // Topology x Phi bins for Fa18Inb
    std::vector<int> raw_yields_Fa18Out;  // Topology x Phi bins for Fa18Out
    std::vector<int> raw_yields_Sp19Inb;  // Topology x Phi bins for Sp19Inb

    std::vector<double> acceptance_Fa18Inb;  // Acceptance for Fa18Inb
    std::vector<double> acceptance_Fa18Out;  // Acceptance for Fa18Out
    std::vector<double> acceptance_Sp19Inb;  // Acceptance for Sp19Inb

    std::vector<double> unfolded_yields_Fa18Inb;  // Unfolded yields for Fa18Inb
    std::vector<double> unfolded_yields_Fa18Out;  // Unfolded yields for Fa18Out
    std::vector<double> unfolded_yields_Sp19Inb;  // Unfolded yields for Sp19Inb
};

// Function prototype for plot_unfolding
std::vector<UnfoldingData> plot_unfolding(const std::string& output_dir, 
                                          const std::string& analysisType, 
                                          int xB_bin,
                                          const std::vector<BinBoundary>& bin_boundaries, 
                                          std::vector<TTreeReader>& data_readers,  // Pass by reference
                                          std::vector<TTreeReader>& mc_gen_readers,  // Pass by reference
                                          std::vector<TTreeReader>& mc_rec_readers);  // Pass by reference

#endif // PLOT_UNFOLDING_H