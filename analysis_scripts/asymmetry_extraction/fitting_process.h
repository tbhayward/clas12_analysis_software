// fitting_process.h
#ifndef FITTING_PROCESS_H
#define FITTING_PROCESS_H

void negLogLikelihood_single_hadron(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void performMLMFits_single_hadron(const char* output_file, const char* kinematic_file, const std::string& prefix);
void plotHistogramAndFit_single_hadron(TH1D* histogram, TF1* fitFunction, int binIndex, int asymmetryIndex, const std::string& prefix);
TH1D* createHistogramForBin_single_hadron(const char* histName, int binIndex, const std::string& prefix, int asymmetry_index);
void performChi2Fits_single_hadron(const char* output_file, const char* kinematic_file, const std::string& prefix, int asymmetry_index);

#endif // FITTING_PROCESS_H
