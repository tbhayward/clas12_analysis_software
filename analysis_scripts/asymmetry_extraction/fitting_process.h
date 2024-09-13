// fitting_process.h
#ifndef FITTING_PROCESS_H
#define FITTING_PROCESS_H

#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <Rtypes.h>  // For ROOT types like Int_t, Double_t
#include "TH1D.h"
#include "TH2D.h" // Include for 2D histograms
#include "TF1.h"
#include "TF2.h" // Include for 2D fit functions

// single hadron
void negLogLikelihood_single_hadron(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void performMLMFits_single_hadron(const char* output_file, const char* kinematic_file, const std::string& prefix, const std::vector<std::pair<double, double>>& dilutionFactors);
void plotHistogramAndFit_single_hadron(TH1D* histogram, TF1* fitFunction, int binIndex, int asymmetryIndex, const std::string& prefix);
TH1D* createHistogramForBin_single_hadron(const char* histName, int binIndex, const std::string& prefix, int asymmetry_index, const std::vector<std::pair<double, double>>& dilutionFactors);
void performChi2Fits_single_hadron(const char* output_file, const char* kinematic_file, const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index, const std::vector<std::pair<double, double>>& dilutionFactors);

#endif // FITTING_PROCESS_H