#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
// tbhayward libraries
#include "common_vars.h"  // Include the common header
#include "load_bins_from_csv.h"
#include "load_run_info_from_csv.h"
#include "dilution_factor.h"
#include "calculate_dilution_factors.h"
#include "asymmetry_fits.h"
#include "BaseKinematicCuts.h"
#include "KinematicCuts.h"
#include "InclusiveKinematicCuts.h"
#include "SingleHadronKinematicCuts.h"
#include "B2BDihadronKinematicCuts.h"
#include "DihadronKinematicCuts.h"
#include "dvcsKinematicCuts.h"
#include "formatLabelName.h"
#include "readChi2Fits.h"
#include "histConfigs.h"
#include "charge_accumulation.h"
#include "plot_data.h"
#include "modifyTree.h"
#include "fitting_process.h" // Include your header file

// Fractional charge values
const double xA = 0.71297;
const double xC = 0.08399;
const double xCH = 0.03600;
const double xHe = 0.07535;
const double xf = 0.09168;

std::vector<std::pair<double, double>> calculate_dilution_factors() {

    // Load ROOT files and trees
    TFile* nh3File = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/cphi2024_epX/data/processed_data/dilution/rgc_su22_inb_NH3_epX.root");
    TFile* cFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/cphi2024_epX/data/processed_data/dilution/rgc_su22_inb_C_epX.root");
    TFile* chFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/cphi2024_epX/data/processed_data/dilution/rgc_su22_inb_CH2_epX.root");
    TFile* heFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/cphi2024_epX/data/processed_data/dilution/rgc_su22_inb_He_bath_epX.root");
    TFile* emptyFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/cphi2024_epX/data/processed_data/dilution/rgc_su22_inb_ET_epX.root");

    TTree* nh3 = (TTree*)nh3File->Get("PhysicsEvents");
    TTree* c = (TTree*)cFile->Get("PhysicsEvents");
    TTree* ch = (TTree*)chFile->Get("PhysicsEvents");
    TTree* he = (TTree*)heFile->Get("PhysicsEvents");
    TTree* empty = (TTree*)emptyFile->Get("PhysicsEvents");

    // Create local TTreeReader objects for each tree
    TTreeReader nh3Reader(nh3);
    TTreeReader cReader(c);
    TTreeReader chReader(ch);
    TTreeReader heReader(he);
    TTreeReader emptyReader(empty);

    // Pointers for local kinematic cuts, dynamically allocated based on the channel
    BaseKinematicCuts* nh3Cuts = nullptr;
    BaseKinematicCuts* cCuts = nullptr;
    BaseKinematicCuts* chCuts = nullptr;
    BaseKinematicCuts* heCuts = nullptr;
    BaseKinematicCuts* emptyCuts = nullptr;

    // Allocate the appropriate kinematic cuts based on the channel
    switch (channel) {
        case 0:
            nh3Cuts = new InclusiveKinematicCuts(nh3Reader);
            cCuts = new InclusiveKinematicCuts(cReader);
            chCuts = new InclusiveKinematicCuts(chReader);
            heCuts = new InclusiveKinematicCuts(heReader);
            emptyCuts = new InclusiveKinematicCuts(emptyReader);
            break;
        case 1:
            nh3Cuts = new SingleHadronKinematicCuts(nh3Reader);
            cCuts = new SingleHadronKinematicCuts(cReader);
            chCuts = new SingleHadronKinematicCuts(chReader);
            heCuts = new SingleHadronKinematicCuts(heReader);
            emptyCuts = new SingleHadronKinematicCuts(emptyReader);
            break;
        case 2:
            nh3Cuts = new B2BDihadronKinematicCuts(nh3Reader);
            cCuts = new B2BDihadronKinematicCuts(cReader);
            chCuts = new B2BDihadronKinematicCuts(chReader);
            heCuts = new B2BDihadronKinematicCuts(heReader);
            emptyCuts = new B2BDihadronKinematicCuts(emptyReader);
            break;
        case 3:
            nh3Cuts = new DihadronKinematicCuts(nh3Reader);
            cCuts = new DihadronKinematicCuts(cReader);
            chCuts = new DihadronKinematicCuts(chReader);
            heCuts = new DihadronKinematicCuts(heReader);
            emptyCuts = new DihadronKinematicCuts(emptyReader);
            break;
        case 4:
            nh3Cuts = new dvcsKinematicCuts(nh3Reader);
            cCuts = new dvcsKinematicCuts(cReader);
            chCuts = new dvcsKinematicCuts(chReader);
            heCuts = new dvcsKinematicCuts(heReader);
            emptyCuts = new dvcsKinematicCuts(emptyReader);
            break;
        default:
            std::cerr << "Invalid channel specified." << std::endl;
            return {}; // Return an empty vector to indicate failure
    }

    std::vector<std::pair<double, double>> dilutionResults;
    TGraphErrors* gr_dilution = new TGraphErrors();

    // Loop over each bin
    for (size_t binIndex = 0; binIndex < allBins[currentFits].size() - 1; ++binIndex) {
        double varMin = allBins[currentFits][binIndex];
        double varMax = allBins[currentFits][binIndex + 1];

        // Create histograms for each target type
        TH1D *h_nh3 = new TH1D("h_nh3", "", 1, varMin, varMax);
        TH1D *h_c = new TH1D("h_c", "", 1, varMin, varMax);
        TH1D *h_ch = new TH1D("h_ch", "", 1, varMin, varMax);
        TH1D *h_he = new TH1D("h_he", "", 1, varMin, varMax);
        TH1D *h_empty = new TH1D("h_empty", "", 1, varMin, varMax);

        double sumCurrentVariable = 0.0;
        int count = 0;

        // Helper function to fill histograms based on kinematic cuts and track mean
        auto fill_histogram = [&](TTreeReader& reader, TH1D* hist, BaseKinematicCuts* cuts) {
            TTreeReaderValue<double> currentVariable(reader, propertyNames[currentFits].c_str());

            while (reader.Next()) {
                bool passedKinematicCuts = cuts->applyCuts(currentFits, false);
                if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
                    hist->Fill(*currentVariable);
                    sumCurrentVariable += *currentVariable;
                    ++count;
                }
            }
        };

        fill_histogram(nh3Reader, h_nh3, nh3Cuts);
        fill_histogram(cReader, h_c, cCuts);
        fill_histogram(chReader, h_ch, chCuts);
        fill_histogram(heReader, h_he, heCuts);
        fill_histogram(emptyReader, h_empty, emptyCuts);

        // Calculate the mean value of currentVariable in this bin
        double meanCurrentVariable = (count > 0) ? (sumCurrentVariable / count) : (varMin + varMax) / 2.0;

        // Retrieve bin contents
        double nA = h_nh3->GetBinContent(1);
        double nC = h_c->GetBinContent(1);
        double nCH = h_ch->GetBinContent(1);
        double nMT = h_he->GetBinContent(1);
        double nf = h_empty->GetBinContent(1);

        // Calculate dilution factor for this bin
        double dilution = (23.0 * (-nMT * xA + nA * xHe) * 
                           (-0.511667 * nMT * xC * xCH * xf + 
                            (1.0 * nf * xC * xCH - 
                             3.41833 * nCH * xC * xf + 
                             2.93 * nC * xCH * xf) * xHe)) / 
                          (nA * xHe * 
                           (62.6461 * nMT * xC * xCH * xf + 
                            1.0 * nf * xC * xCH * xHe - 
                            78.6217 * nCH * xC * xf * xHe + 
                            14.9756 * nC * xCH * xf * xHe));
        // Calculate error for this bin
        double error = TMath::Sqrt((TMath::Power(nA, 2) / (xA * xA)) + 
                                   (TMath::Power(nC, 2) / (xC * xC)) + 
                                   (TMath::Power(nCH, 2) / (xCH * xCH)) + 
                                   (TMath::Power(nMT, 2) / (xHe * xHe)) + 
                                   (TMath::Power(nf, 2) / (xf * xf)));
        std::cout << dilution << " " << error << std::endl;

        // Add the dilution factor and error to the TGraphErrors
        gr_dilution->SetPoint(binIndex, meanCurrentVariable, dilution);
        gr_dilution->SetPointError(binIndex, 0, error);

        // Store the original dilution and error for now
        dilutionResults.emplace_back(dilution, error);

        // Clean up histograms
        delete h_nh3;
        delete h_c;
        delete h_ch;
        delete h_he;
        delete h_empty;
    }

    // Fit the TGraphErrors to a cubic polynomial
    TF1* fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x^2 + [3]*x^3", allBins[currentFits].front(), allBins[currentFits].back());
    gr_dilution->Fit(fitFunc, "Q");

    // Calculate the chi2/ndf
    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    double scaleFactor = TMath::Sqrt(chi2 / ndf);

    // Scale the errors by sqrt(chi2/ndf)
    for (size_t binIndex = 0; binIndex < dilutionResults.size(); ++binIndex) {
        dilutionResults[binIndex].second *= scaleFactor;  // Scale the error
    }

    // Clean up
    delete gr_dilution;
    delete fitFunc;

    // Clean up dynamically allocated kinematic cuts
    delete nh3Cuts;
    delete cCuts;
    delete chCuts;
    delete heCuts;
    delete emptyCuts;

    // Close and delete ROOT files
    nh3File->Close(); delete nh3File;
    cFile->Close(); delete cFile;
    chFile->Close(); delete chFile;
    heFile->Close(); delete heFile;
    emptyFile->Close(); delete emptyFile;

    return dilutionResults;
}