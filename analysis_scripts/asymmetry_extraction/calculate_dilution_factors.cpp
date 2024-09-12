#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
// tbhayward libraries
#include "common_vars.h"
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
#include "fitting_process.h"

// Fractional charge values for Total
const double xAtotal = 0.72032;
const double xCtotal = 0.08184;
const double xCHtotal = 0.03508;
const double xHetotal = 0.07343;
const double xftotal = 0.08934;

// Fractional charge values for Period 1
const double xAperiod_1 = 0.16639;
const double xCperiod_1 = 0.24392;
const double xCHperiod_1 = 0.10457;
const double xHeperiod_1 = 0.21885;
const double xfperiod_1 = 0.26627;

// Fractional charge values for Period 2
const double xAperiod_2 = 0.08956;
const double xCperiod_2 = 0.26640;
const double xCHperiod_2 = 0.11420;
const double xHeperiod_2 = 0.23902;
const double xfperiod_2 = 0.29081;

// Fractional charge values for Period 3
const double xAperiod_3 = 0.13839;
const double xCperiod_3 = 0.25211;
const double xCHperiod_3 = 0.10808;
const double xHeperiod_3 = 0.22620;
const double xfperiod_3 = 0.27522;

// Fractional charge values for Period 4
const double xAperiod_4 = 0.26162;
const double xCperiod_4 = 0.21605;
const double xCHperiod_4 = 0.09262;
const double xHeperiod_4 = 0.19385;
const double xfperiod_4 = 0.23585;

// Fractional charge values for Period 5
const double xAperiod_5 = 0.19162;
const double xCperiod_5 = 0.23654;
const double xCHperiod_5 = 0.10140;
const double xHeperiod_5 = 0.21223;
const double xfperiod_5 = 0.25821;

// Fractional charge values for Period 6
const double xAperiod_6 = 0.27241;
const double xCperiod_6 = 0.21289;
const double xCHperiod_6 = 0.09127;
const double xHeperiod_6 = 0.19102;
const double xfperiod_6 = 0.23241;

// Fractional charge values for Period 7
const double xAperiod_7 = 0.18628;
const double xCperiod_7 = 0.23810;
const double xCHperiod_7 = 0.10207;
const double xHeperiod_7 = 0.21363;
const double xfperiod_7 = 0.25992;

// Fractional charge values for Period 8
const double xAperiod_8 = 0.40136;
const double xCperiod_8 = 0.17517;
const double xCHperiod_8 = 0.07509;
const double xHeperiod_8 = 0.15717;
const double xfperiod_8 = 0.19122;

// Fractional charge values for Period 9
const double xAperiod_9 = 0.13814;
const double xCperiod_9 = 0.25218;
const double xCHperiod_9 = 0.10811;
const double xHeperiod_9 = 0.22627;
const double xfperiod_9 = 0.27530;

double calculate_dilution_error(double nA, double nC, double nCH, double nMT, double nf, 
                                double xA, double xC, double xCH, double xHe, double xf) {
    // First part of the expression
    double term1 = 3988.9 * nA * nf * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) * pow(xHe, 2) *
                  pow((-1.0 * nMT * xA + nA * xHe), 2) *
                  pow((1.0 * nMT * xC * xCH - 1.19072 * nCH * xC * xf + 0.190722 * nC * xCH * xHe), 2);

    // Second part of the expression
    double term2 = 64705.7 * nA * nCH * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) * pow(xHe, 2) *
                  pow((-1.0 * nMT * xA + nA * xHe), 2) *
                  pow((1.0 * nMT * xC * xf - 0.295642 * nf * xC * xHe - 0.704359 * nC * xf * xHe), 2);

    // Third part of the expression
    double term3 = 8.5849 * nA * nC * pow(xCH, 2) * pow(xf, 2) * pow(xHe, 2) *
                  pow((-1.0 * nMT * xA + nA * xHe), 2) *
                  pow((65.2613 * nMT * xC * xCH * xf + (-4.11113 * nf * xC * xCH - 61.1502 * nCH * xC * xf) * xHe), 2);

    // Fourth part of the expression
    double term4 = 1027.46 * pow(nMT, 2) * pow(xA, 2) *
                  pow((1.0 * nMT * xC * xCH * xf + (-1.9544 * nf * xC * xCH + 6.68077 * nCH * xC * xf - 5.72638 * nC * xCH * xf) * xHe), 2) *
                  pow((1.0 * nMT * xC * xCH * xf + (0.0159627 * nf * xC * xCH - 1.25501 * nCH * xC * xf + 0.239051 * nC * xCH * xf) * xHe), 2);

    // Fifth part of the expression
    double term5 = 0.261803 * nA * nMT *
                  pow((62.6461 * pow(nMT, 2) * xA * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) +
                       nMT * xA * xC * xCH * xf * (2.0 * nf * xC * xCH - 157.243 * nCH * xC * xf + 29.9512 * nC * xCH * xf) * xHe +
                       (-1.9544 * pow(nf, 2) * xA * pow(xC, 2) * pow(xCH, 2) +
                        nf * xC * xCH * (160.339 * nCH * xA * xC - 34.9946 * nC * xA * xCH - 123.435 * nA * xC * xCH) * xf +
                        (-525.254 * pow(nCH, 2) * xA * pow(xC, 2) +
                         nCH * xC * (550.266 * nC * xA + 497.146 * nA * xC) * xCH +
                         nC * (-85.756 * nC * xA - 373.711 * nA * xC) * pow(xCH, 2)) * xf)) *
                      pow(xHe, 2),2);

    // Denominator of the expression
    double denominator = pow(nA, 3) * pow(xHe, 2) * pow((62.6461 * nMT * xC * xCH * xf + 1.0 * nf * xC * xCH * xHe - 78.6217 * nCH * xC * xf * xHe + 14.9756 * nC * xCH * xf * xHe), 4);

    // Final error calculation
    double sigma_df = 23.0 * sqrt((term1 + term2 + term3 + term4 + term5) / denominator);

    return sigma_df;
}

std::pair<double, double> calculate_dilution_and_error(double nA, double nC, double nCH, double nMT, double nf, 
                                                       double xA, double xC, double xCH, double xHe, double xf) {
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
    
    double error = calculate_dilution_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf);
    
    return std::make_pair(dilution, error);
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

    // Read helicity and target polarization
    TTreeReaderValue<int> runnum(nh3Reader, "runnum");
    TTreeReaderValue<int> helicity(nh3Reader, "helicity");
    TTreeReaderValue<double> target_pol(nh3Reader, "target_pol");

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
        TGraphErrors* gr_dilution[10];  // Array of pointers to TGraphErrors
        for (int i = 0; i < 10; ++i) {
            gr_dilution[i] = new TGraphErrors();  // Create a new TGraphErrors for each index
        }

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
            auto fill_histogram = [&](TTreeReader& reader, TH1D* hist, BaseKinematicCuts* cuts, 
                bool is_nh3) {
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

            // Call fill_histogram for each target type
            fill_histogram(nh3Reader, h_nh3, nh3Cuts, true);  // NH3 data
            fill_histogram(cReader, h_c,  cCuts, false);                      // Carbon data
            fill_histogram(chReader, h_ch, chCuts, false);                   // CH2 data
            fill_histogram(heReader, h_he,  heCuts, false);                   // Helium data
            fill_histogram(emptyReader, h_empty,  emptyCuts, false);  

            // Calculate the mean value of currentVariable in this bin
            double meanCurrentVariable = (count > 0) ? (sumCurrentVariable / count) : (varMin + varMax) / 2.0;

            // Retrieve bin contents
            double nA = h_nh3->GetBinContent(1);
            double nC = h_c->GetBinContent(1);
            double nCH = h_ch->GetBinContent(1);
            double nMT = h_he->GetBinContent(1);
            double nf = h_empty->GetBinContent(1);

            /// Calculate dilution factors for the general case
            auto [dilution, error] = calculate_dilution_and_error(nA, nC, nCH, nMT, nf, xAtotal, xCtotal, xCHtotal, xHetotal, xftotal);

            // Add the dilution factor and error to the TGraphErrors
            gr_dilution[0]->SetPoint(binIndex, meanCurrentVariable, dilution);
            gr_dilution[0]->SetPointError(binIndex, 0, error);

            // Store the original dilution and error for now
            dilutionResults.emplace_back(dilution, error);

            // Clean up histograms
            delete h_nh3;
            delete h_c;
            delete h_ch;
            delete h_he;
            delete h_empty;
        }

        // Plot the original dilution factor
        std::string prefix = propertyNames[currentFits];
        HistConfig config = histConfigs.find(prefix) != histConfigs.end() ? histConfigs[prefix] : HistConfig{100, 0, 1};

        TCanvas* canvas = new TCanvas("c_dilution", "Dilution Factor Plot", 800, 600);
        canvas->SetLeftMargin(0.15);
        canvas->SetBottomMargin(0.15);

        gr_dilution[0]->SetTitle("");
        gr_dilution[0]->GetXaxis()->SetTitle(formatLabelName(prefix).c_str());
        gr_dilution[0]->GetXaxis()->SetLimits(config.xMin, config.xMax);
        gr_dilution[0]->GetXaxis()->SetTitleSize(0.05);
        gr_dilution[0]->GetYaxis()->SetTitle("D_{f}");
        gr_dilution[0]->GetYaxis()->SetTitleSize(0.05);
        gr_dilution[0]->GetYaxis()->SetTitleOffset(1.6);
        gr_dilution[0]->GetYaxis()->SetRangeUser(0.0, 0.4);
        gr_dilution[0]->SetMarkerStyle(20);
        gr_dilution[0]->SetMarkerColor(kBlack);
        gr_dilution[0]->Draw("AP");

        std::string outputDir = "output/dilution_factor_plots/";
        std::string outputFileName = outputDir + "df_" + binNames[currentFits] + "_" + prefix + ".png";
        canvas->SaveAs(outputFileName.c_str());


        // Clean up
        delete canvas;
        delete gr_dilution;
        delete fitFunc;

        delete nh3Cuts;
        delete cCuts;
        delete chCuts;
        delete heCuts;
        delete emptyCuts;

        nh3File->Close(); delete nh3File;
        cFile->Close(); delete cFile;
        chFile->Close(); delete chFile;
        heFile->Close(); delete heFile;
        emptyFile->Close(); delete emptyFile;

        return dilutionResults;
    }