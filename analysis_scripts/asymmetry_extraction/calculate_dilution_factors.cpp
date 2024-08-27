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

double calculate_dilution_error(double nA, double nC, double nCH, double nMT, double nf) {
    double term1 = 0.734694 * nA * nC * pow((nA - nMT), 2) * 
                   pow((- 0.554976 * nCH - 0.0373109 * nf + 0.592287 * nMT), 2);

    double term2 = 0.456108 * nA * pow((nA - nMT), 2) * 
                   pow((0.704358 * nC + 0.295642 * nf - nMT), 2) * nMT;

    double term3 = 0.0281176 * nA * nf * pow((nA - nMT), 2) * 
                   pow((0.190722 * nC - 1.19072 * nCH + nMT), 2);

    double term4 = 0.0135714 * pow((nC - 1.16667 * nCH + 0.341297 * nf - 0.17463 * nMT), 2) * 
                   pow(nMT, 2) * pow((nC - 5.25 * nCH + 0.0667755 * nf + 4.18322 * nMT), 2);

    double term5 = 0.022405 * nA * nMT * 
                   pow((0.778288 * pow(nC, 2) + 4.76701 * pow(nCH, 2) - 1.45518 * nCH * nf + 0.0177374 * pow(nf, 2) + 
                        nC * (-4.99401 * nCH + 0.317598 * nf - 0.271825 * nMT) + 
                        nA * (3.39167 * nC - 4.51192 * nCH + 1.12025 * nf) + 
                        1.42708 * nCH * nMT - 0.0181513 * nf * nMT - 0.568553 * pow(nMT, 2)), 2);

    double denominator = pow(nA, 3) * 
                         pow((0.135913 * nC - 0.713541 * nCH + 0.00907563 * nf + 0.568553 * nMT), 4);

    double sigma_df = 0.713541 * sqrt((term1 + term2 + term3 + term4 + term5 )/ denominator);

    return sigma_df;
}

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
            reader.Restart();
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
        double error = calculate_dilution_error(nA, nC, nCH, nMT, nf);

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
    gr_dilution->Fit(fitFunc, "QN");

    // Calculate the chi2/ndf
    double chi2 = fitFunc->GetChisquare();
    int ndf = fitFunc->GetNDF();
    double scalingFactor = std::sqrt(chi2 / ndf);

    // // Loop over all points in gr_dilution and scale the errors
    // int nPoints = gr_dilution->GetN();
    // for (int i = 0; i < nPoints; ++i) {
    //     double x, y;
    //     gr_dilution->GetPoint(i, x, y);
    //     double ey = gr_dilution->GetErrorY(i);
    //     gr_dilution->SetPointError(i, gr_dilution->GetErrorX(i), ey * scalingFactor);
    // }

    // Now plot the TGraphErrors
    // Get the prefix from propertyNames
    std::string prefix = propertyNames[currentFits];

    // Retrieve the HistConfig for the current prefix
    HistConfig config = histConfigs.find(prefix) != histConfigs.end() ? histConfigs[prefix] : HistConfig{100, 0, 1};

    // Create a canvas for plotting
    TCanvas* canvas = new TCanvas("c_dilution", "Dilution Factor Plot", 800, 600);
    canvas->SetLeftMargin(0.15);
    canvas->SetBottomMargin(0.15);

    // Configure the graph
    gr_dilution->SetTitle("");
    gr_dilution->GetXaxis()->SetTitle(formatLabelName(prefix).c_str());
    gr_dilution->GetXaxis()->SetLimits(config.xMin, config.xMax);
    gr_dilution->GetXaxis()->SetTitleSize(0.05);
    gr_dilution->GetYaxis()->SetTitle("D_{f}");
    gr_dilution->GetYaxis()->SetTitleSize(0.05);
    gr_dilution->GetYaxis()->SetTitleOffset(1.6);
    gr_dilution->GetYaxis()->SetRangeUser(0.0, 0.4);

    // Draw the graph with error bars
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->SetMarkerColor(kBlack);
    gr_dilution->Draw("AP");

    // Save the canvas to a file
    std::string outputDir = "output/dilution_factor_plots/";
    std::string outputFileName = outputDir + "df_" + binNames[currentFits] + "_" + prefix + ".png";
    canvas->SaveAs(outputFileName.c_str());

    // Clean up
    delete canvas;
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