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

// inclusive
std::tuple<int, int, int, int, double, double, double> getInclusiveCounts(int binIndex, 
  const std::string& prefix);
void calculate_inclusive(const char* output_file, const char* kinematic_file,
  const std::string& prefix, int asymmetry_index);

// single hadron
void negLogLikelihood_single_hadron(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void performMLMFits_single_hadron(const char* output_file, const char* kinematic_file, const std::string& prefix);
void plotHistogramAndFit_single_hadron(TH1D* histogram, TF1* fitFunction, int binIndex, int asymmetryIndex, const std::string& prefix);
TH1D* createHistogramForBin_single_hadron(const char* histName, int binIndex, const std::string& prefix, int asymmetry_index);
void performChi2Fits_single_hadron(const char* output_file, const char* kinematic_file, const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index);


// b2b dihadron (dSIDIS)
void negLogLikelihood_b2b_dihadron(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void performMLMFits_b2b_dihadron(const char* output_file, const char* kinematic_file, const std::string& prefix);
TH2D* createHistogramForBin_b2b_dihadron(const char* histName, int binIndex, const std::string& prefix, int asymmetry_index);
void performChi2Fits_b2b_dihadron(const char* output_file, const char* kinematic_file, const std::string& prefix, int asymmetry_index);


// dvcs
void negLogLikelihood_dvcs(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void performMLMFits_dvcs(const char* output_file, const char* kinematic_file, const std::string& prefix);
void plotHistogramAndFit_dvcs(TH1D* histogram, TF1* fitFunction, int binIndex, int asymmetryIndex, const std::string& prefix);
TH1D* createHistogramForBin_dvcs(const char* histName, int binIndex, const std::string& prefix, int asymmetry_index);
void performChi2Fits_dvcs(const char* output_file, const char* kinematic_file, const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index);


// eppi0
void plotHistogramAndFit_eppi0(TH1D* histogram, TF1* fitFunction, int binIndex, int asymmetryIndex, const std::string& prefix);
TH1D* createHistogramForBin_eppi0(const char* histName, int binIndex, const std::string& prefix, int asymmetry_index);
void performChi2Fits_eppi0(const char* output_file, const char* kinematic_file, const char* kinematicPlot_file, const std::string& prefix, int asymmetry_index);

// === GeneralExclusive (simultaneous) ===
void performChi2Fits_GeneralExclusive(const char* outFile, const char* kinFile,
  const char* kinPlotFile, const std::string& prefix);

// Build the three phi-binned asymmetry histograms for a given bin.
void createHistogramForBin_GeneralExclusive(const std::string& prefix, int binIndex,
  TH1D*& hBSA, TH1D*& hTSA, TH1D*& hDSA);

// Plot hists + best-fit curves (3 pads) and save a PNG per bin.
void plotHistogramAndFit_GeneralExclusive(TH1D* hBSA, TH1D* hTSA, TH1D* hDSA, const double pars[7],
  // dep means for overlay curves
  double depA, double depB, double depC, double depV, double depW,
  const std::string& prefix, int binIndex, const char* kinPlotFile);

#endif // FITTING_PROCESS_H