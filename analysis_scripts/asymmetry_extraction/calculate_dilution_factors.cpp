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
    // Load ROOT files and trees (as before)
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

    // Create local TTreeReader objects (as before)
    TTreeReader nh3Reader(nh3);
    TTreeReader cReader(c);
    TTreeReader chReader(ch);
    TTreeReader heReader(he);
    TTreeReader emptyReader(empty);

    TTreeReaderValue<int> helicity(nh3Reader, "helicity");
    TTreeReaderValue<int> runnum(nh3Reader, "runnum");

    // Declare output directory
    std::string outputDir = "output/dilution_factor_plots/";

    // Read kinematic cuts for the specific channel
    BaseKinematicCuts* nh3Cuts = nullptr;
    BaseKinematicCuts* cCuts = nullptr;
    BaseKinematicCuts* chCuts = nullptr;
    BaseKinematicCuts* heCuts = nullptr;
    BaseKinematicCuts* emptyCuts = nullptr;

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
            return {};
    }

    std::vector<std::pair<double, double>> dilutionResults;
    TGraphErrors* gr_dilution_total = new TGraphErrors();
    std::vector<TGraphErrors*> gr_dilution_periods(9);
    for (int i = 0; i < 9; i++) {
        gr_dilution_periods[i] = new TGraphErrors();
    }

    // NH3 periods defined as a pair of run numbers (start, end)
    std::vector<std::pair<int, int>> nh3_periods = {
        {16137, 16148}, {16156, 16178}, {16211, 16228}, {16231, 16260},
        {16318, 16333}, {16335, 16357}, {16709, 16720}, {16721, 16766}, {16767, 16772}
    };

    // Loop over each bin
    for (size_t binIndex = 0; binIndex < allBins[currentFits].size() - 1; ++binIndex) {
        double varMin = allBins[currentFits][binIndex];
        double varMax = allBins[currentFits][binIndex + 1];

        // Create histograms for total and each period
        TH1D* h_nh3_total = new TH1D("h_nh3_total", "", 1, varMin, varMax);
        std::vector<TH1D*> h_nh3_periods(9);
        for (int i = 0; i < 9; i++) {
            h_nh3_periods[i] = new TH1D(Form("h_nh3_period_%d", i + 1), "", 1, varMin, varMax);
        }

        double sumCurrentVariable = 0.0;
        int count = 0;

        // Fill histograms for NH3 total
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

        fill_histogram(nh3Reader, h_nh3_total, nh3Cuts);

        // Calculate the mean value of `currentVariable` in this bin
        double meanCurrentVariable = (count > 0) ? (sumCurrentVariable / count) : (varMin + varMax) / 2.0;

        // Fill period histograms based on run number
        while (nh3Reader.Next()) {
            bool passedKinematicCuts = nh3Cuts->applyCuts(currentFits, false);
            if (!passedKinematicCuts) continue;

            int run = *runnum;
            for (int i = 0; i < nh3_periods.size(); ++i) {
                if (run >= nh3_periods[i].first && run <= nh3_periods[i].second) {
                    h_nh3_periods[i]->Fill(varMin);  // Fill based on current bin
                    break;
                }
            }
        }

        // Retrieve bin contents for total
        double nA_total = h_nh3_total->GetBinContent(1);
        double nC = cReader.GetTree()->GetEntries();
        double nCH = chReader.GetTree()->GetEntries();
        double nMT = heReader.GetTree()->GetEntries();
        double nf = emptyReader.GetTree()->GetEntries();

        // Calculate dilution factors for the total case
        auto [dilution_total, error_total] = calculate_dilution_and_error(nA_total, nC, nCH, nMT, nf, xAtotal, xCtotal, xCHtotal, xHetotal, xftotal);
        gr_dilution_total->SetPoint(binIndex, meanCurrentVariable, dilution_total); // Plot at the mean value
        gr_dilution_total->SetPointError(binIndex, 0, error_total);

        for (size_t binIndex = 0; binIndex < allBins[currentFits].size() - 1; ++binIndex) {
            double varMin = allBins[currentFits][binIndex];
            double varMax = allBins[currentFits][binIndex + 1];

            // Create histograms for total and each period
            TH1D* h_nh3_total = new TH1D("h_nh3_total", "", 1, varMin, varMax);
            TH1D* h_c = new TH1D("h_c", "", 1, varMin, varMax);
            TH1D* h_ch = new TH1D("h_ch", "", 1, varMin, varMax);
            TH1D* h_he = new TH1D("h_he", "", 1, varMin, varMax);
            TH1D* h_empty = new TH1D("h_empty", "", 1, varMin, varMax);

            double sumCurrentVariable = 0.0;
            int count = 0;

            // Fill histograms for NH3 total
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

            // Fill the histograms for each target type
            fill_histogram(nh3Reader, h_nh3_total, nh3Cuts);
            fill_histogram(cReader, h_c, cCuts);
            fill_histogram(chReader, h_ch, chCuts);
            fill_histogram(heReader, h_he, heCuts);
            fill_histogram(emptyReader, h_empty, emptyCuts);

            // Calculate the mean value of `currentVariable` in this bin
            double meanCurrentVariable = (count > 0) ? (sumCurrentVariable / count) : (varMin + varMax) / 2.0;

            // Retrieve bin contents for total and each target type
            double nA_total = h_nh3_total->GetBinContent(1);
            double nC = h_c->GetBinContent(1);
            double nCH = h_ch->GetBinContent(1);
            double nMT = h_he->GetBinContent(1);
            double nf = h_empty->GetBinContent(1);

            // Calculate dilution factors for the total case
            auto [dilution_total, error_total] = calculate_dilution_and_error(
                nA_total, nC, nCH, nMT, nf, xAtotal, xCtotal, xCHtotal, xHetotal, xftotal
            );
            gr_dilution_total->SetPoint(binIndex, meanCurrentVariable, dilution_total); // Plot at the mean value
            gr_dilution_total->SetPointError(binIndex, 0, error_total);

            for (int i = 0; i < 9; i++) {
                double nA_period = h_nh3_total->GetBinContent(1);  // Period-specific bin content

                // Use the correct fractional charge values for each period (as described before)
                double xA_period, xC_period, xCH_period, xHe_period, xf_period;
                switch (i) {
                    case 0: xA_period = xAperiod_1; xC_period = xCperiod_1; xCH_period = xCHperiod_1; xHe_period = xHeperiod_1; xf_period = xfperiod_1; break;
                    case 1: xA_period = xAperiod_2; xC_period = xCperiod_2; xCH_period = xCHperiod_2; xHe_period = xHeperiod_2; xf_period = xfperiod_2; break;
                    case 2: xA_period = xAperiod_3; xC_period = xCperiod_3; xCH_period = xCHperiod_3; xHe_period = xHeperiod_3; xf_period = xfperiod_3; break;
                    case 3: xA_period = xAperiod_4; xC_period = xCperiod_4; xCH_period = xCHperiod_4; xHe_period = xHeperiod_4; xf_period = xfperiod_4; break;
                    case 4: xA_period = xAperiod_5; xC_period = xCperiod_5; xCH_period = xCHperiod_5; xHe_period = xHeperiod_5; xf_period = xfperiod_5; break;
                    case 5: xA_period = xAperiod_6; xC_period = xCperiod_6; xCH_period = xCHperiod_6; xHe_period = xHeperiod_6; xf_period = xfperiod_6; break;
                    case 6: xA_period = xAperiod_7; xC_period = xCperiod_7; xCH_period = xCHperiod_7; xHe_period = xHeperiod_7; xf_period = xfperiod_7; break;
                    case 7: xA_period = xAperiod_8; xC_period = xCperiod_8; xCH_period = xCHperiod_8; xHe_period = xHeperiod_8; xf_period = xfperiod_8; break;
                    case 8: xA_period = xAperiod_9; xC_period = xCperiod_9; xCH_period = xCHperiod_9; xHe_period = xHeperiod_9; xf_period = xfperiod_9; break;
                }

                // Calculate dilution factor for each period
                auto [dilution_period, error_period] = calculate_dilution_and_error(
                    nA_period, nC, nCH, nMT, nf, 
                    xA_period, xC_period, xCH_period, xHe_period, xf_period
                );
                gr_dilution_periods[i]->SetPoint(binIndex, meanCurrentVariable + 0.005 * (i + 1), dilution_period);
                gr_dilution_periods[i]->SetPointError(binIndex, 0, error_period);
            }

            // Clean up histograms
            delete h_nh3_total;
            delete h_c;
            delete h_ch;
            delete h_he;
            delete h_empty;
        }

        // Print values for each bin
        std::cout << "df_" << propertyNames[currentFits] << "_0: {" << meanCurrentVariable << ", " << dilution_total << ", " << error_total << "}" << std::endl;
        for (int i = 0; i < 9; i++) {
            std::cout << "df_" << propertyNames[currentFits] << "_" << (i + 1) << ": {"
                      << meanCurrentVariable + 0.0025 * (i + 1) << ", " << gr_dilution_periods[i]->GetY()[binIndex] << ", "
                      << gr_dilution_periods[i]->GetEY()[binIndex] << "}" << std::endl;
        }

        // Clean up histograms
        delete h_nh3_total;
        for (int i = 0; i < 9; i++) {
            delete h_nh3_periods[i];
        }
    }

    // Plot and save total dilution factor
    std::string prefix = propertyNames[currentFits];
    HistConfig config = histConfigs.find(prefix) != histConfigs.end() ? histConfigs[prefix] : HistConfig{100, 0, 1};

    TCanvas* canvas1 = new TCanvas("c_dilution_total", "Total Dilution Factor", 800, 600);
    canvas1->SetLeftMargin(0.15);
    canvas1->SetBottomMargin(0.15);
    gr_dilution_total->SetTitle("");
    gr_dilution_total->GetXaxis()->SetTitle(formatLabelName(propertyNames[currentFits]).c_str());
    gr_dilution_total->GetXaxis()->SetLimits(config.xMin, config.xMax);
    gr_dilution_total->GetXaxis()->SetTitleSize(0.05);
    gr_dilution_total->GetYaxis()->SetTitle("D_{f}");
    gr_dilution_total->GetYaxis()->SetTitleSize(0.05);
    gr_dilution_total->GetYaxis()->SetTitleOffset(1.6);
    gr_dilution_total->GetYaxis()->SetRangeUser(0.0, 0.4);
    gr_dilution_total->SetMarkerStyle(20);
    gr_dilution_total->SetMarkerColor(kBlack);
    gr_dilution_total->Draw("AP");

    canvas1->SaveAs((outputDir + "df_total_" + propertyNames[currentFits] + ".png").c_str());

    // Plot and save period-based dilution factor
    TCanvas* canvas2 = new TCanvas("c_dilution_periods", "Dilution Factor by Periods", 800, 600);
    canvas2->SetLeftMargin(0.15);
    canvas2->SetBottomMargin(0.15);
    gr_dilution_periods[0]->SetTitle("");
    gr_dilution_periods[0]->GetXaxis()->SetTitle(formatLabelName(propertyNames[currentFits]).c_str());
    gr_dilution_periods[0]->GetXaxis()->SetLimits(config.xMin, config.xMax);
    gr_dilution_periods[0]->GetXaxis()->SetTitleSize(0.05);
    gr_dilution_periods[0]->GetYaxis()->SetTitle("D_{f}");
    gr_dilution_periods[0]->GetYaxis()->SetTitleSize(0.05);
    gr_dilution_periods[0]->GetYaxis()->SetTitleOffset(1.6);
    gr_dilution_periods[0]->GetYaxis()->SetRangeUser(0.0, 0.4);

    int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet, kYellow, kPink};
    int markers[] = {21, 22, 23, 24, 25, 26, 27, 28, 29};
    for (int i = 0; i < 9; i++) {
        gr_dilution_periods[i]->SetMarkerStyle(markers[i]);
        gr_dilution_periods[i]->SetMarkerColor(colors[i]);
        if (i == 0) {
            gr_dilution_periods[i]->Draw("AP");
        } else {
            gr_dilution_periods[i]->Draw("P SAME");
        }
    }

    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    for (int i = 0; i < 9; i++) {
        legend->AddEntry(gr_dilution_periods[i], Form("Period %d", i + 1), "p");
    }
    legend->Draw();

    canvas2->SaveAs((outputDir + "df_by_periods_" + propertyNames[currentFits] + ".png").c_str());

    // Clean up
    delete canvas1;
    delete canvas2;
    delete gr_dilution_total;
    for (int i = 0; i < 9; i++) {
        delete gr_dilution_periods[i];
    }

    return dilutionResults;
}