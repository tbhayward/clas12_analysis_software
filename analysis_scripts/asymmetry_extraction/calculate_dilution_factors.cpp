#include <TGraphErrors.h>
#include <TF1.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <cmath>
#include <vector>
#include <numeric>
#include <TLatex.h>
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
const double xAtotal = 0.76103;
const double xCtotal = 0.06223;
const double xCHtotal = 0.03246;
const double xHetotal = 0.06509;
const double xftotal = 0.07919;

// Fractional charge values for Period 1
const double xAperiod_1 = 0.17657;
const double xCperiod_1 = 0.21442;
const double xCHperiod_1 = 0.11185;
const double xHeperiod_1 = 0.22428;
const double xfperiod_1 = 0.27288;

// Fractional charge values for Period 2
const double xAperiod_2 = 0.09580;
const double xCperiod_2 = 0.23545;
const double xCHperiod_2 = 0.12282;
const double xHeperiod_2 = 0.24628;
const double xfperiod_2 = 0.29965;

// Fractional charge values for Period 3
const double xAperiod_3 = 0.15864;
const double xCperiod_3 = 0.21909;
const double xCHperiod_3 = 0.11428;
const double xHeperiod_3 = 0.22917;
const double xfperiod_3 = 0.27882;

// Fractional charge values for Period 4
const double xAperiod_4 = 0.25309;
const double xCperiod_4 = 0.19449;
const double xCHperiod_4 = 0.10145;
const double xHeperiod_4 = 0.20344;
const double xfperiod_4 = 0.24752;

// Fractional charge values for Period 5
const double xAperiod_5 = 0.20247;
const double xCperiod_5 = 0.20767;
const double xCHperiod_5 = 0.10833;
const double xHeperiod_5 = 0.21723;
const double xfperiod_5 = 0.26430;

// Fractional charge values for Period 6
const double xAperiod_6 = 0.28764;
const double xCperiod_6 = 0.18550;
const double xCHperiod_6 = 0.09676;
const double xHeperiod_6 = 0.19403;
const double xfperiod_6 = 0.23607;

// Fractional charge values for Period 7
const double xAperiod_7 = 0.19732;
const double xCperiod_7 = 0.20902;
const double xCHperiod_7 = 0.10903;
const double xHeperiod_7 = 0.21863;
const double xfperiod_7 = 0.26601;

// Fractional charge values for Period 8
const double xAperiod_8 = 0.40967;
const double xCperiod_8 = 0.15372;
const double xCHperiod_8 = 0.08019;
const double xHeperiod_8 = 0.16079;
const double xfperiod_8 = 0.19563;

// Fractional charge values for Period 9
const double xAperiod_9 = 0.14680;
const double xCperiod_9 = 0.22217;
const double xCHperiod_9 = 0.11589;
const double xHeperiod_9 = 0.23239;
const double xfperiod_9 = 0.28274;

// NH3 periods defined as a pair of run numbers (start, end)
std::vector<std::pair<int, int>> nh3_periods = {
    {16137, 16148}, {16156, 16178}, {16211, 16228}, {16231, 16260},
    {16318, 16333}, {16335, 16357}, {16709, 16720}, {16721, 16766}, {16767, 16772}
};

// Function to calculate standard deviation of a vector
double calculate_standard_deviation(const std::vector<double>& values) {
    if (values.size() < 2) return 0.0;
    double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
    double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0,
                                       std::plus<double>(), [&](double a, double b) { return (b - mean) * (b - mean); });
    return std::sqrt(sq_sum / (values.size() - 1));
}

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
    // double dilution = (23.0 * (-nMT * xA + nA * xHe) * 
    //                    (-0.511667 * nMT * xC * xCH * xf + 
    //                     (1.0 * nf * xC * xCH - 
    //                      3.41833 * nCH * xC * xf + 
    //                      2.93 * nC * xCH * xf) * xHe)) / 
    //                   (nA * xHe * 
    //                    (62.6461 * nMT * xC * xCH * xf + 
    //                     1.0 * nf * xC * xCH * xHe - 
    //                     78.6217 * nCH * xC * xf * xHe + 
    //                     14.9756 * nC * xCH * xf * xHe));

    // double packing_fraction = (0.565345 * xC * xCH * xf * (nMT * xA - 1.0 * nA * xHe)) /
    //                       (xA * (1.0 * nMT * xC * xCH * xf +
    //                              0.0159627 * nf * xC * xCH * xHe -
    //                              1.25501 * nCH * xC * xf * xHe +
    //                              0.23905 * nC * xCH * xf * xHe));

    double dilution = (0.699832)*(nA/xA-nMT/xHe)/(1.25055*nCH/xCH-0.23688*nC/xC-0.013668*nf/xf-nMT/xMT);
    
    double error = calculate_dilution_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf);
    
    return std::make_pair(dilution, error);
}

std::vector<std::pair<double, double>> calculate_dilution_factors() {

    // Load ROOT files and trees
    TFile* nh3File = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epX/rgc_su22_inb_NH3_epX.root");
    TFile* cFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epX/rgc_su22_inb_C_epX.root");
    TFile* chFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epX/dilution_factor/rgc_su22_inb_CH2_epX.root");
    TFile* heFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epX/dilution_factor/rgc_su22_inb_He_epX.root");
    TFile* emptyFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epX/dilution_factor/rgc_su22_inb_ET_epX.root");

    // TFile* nh3File = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+X/rgc_su22_inb_NH3_epi+X.root");
    // TFile* cFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+X/rgc_su22_inb_C_epi+X.root");
    // TFile* chFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+X/dilution_factor/rgc_su22_inb_CH2_epi+X.root");
    // TFile* heFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+X/dilution_factor/rgc_su22_inb_He_epi+X.root");
    // TFile* emptyFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+X/dilution_factor/rgc_su22_inb_ET_epi+X.root");

    // TFile* nh3File = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi-X/rgc_su22_inb_NH3_epi-X.root");
    // TFile* cFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi-X/rgc_su22_inb_C_epi-X.root");
    // TFile* chFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi-X/dilution_factor/rgc_su22_inb_CH2_epi-X.root");
    // TFile* heFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi-X/dilution_factor/rgc_su22_inb_He_epi-X.root");
    // TFile* emptyFile = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi-X/dilution_factor/rgc_su22_inb_ET_epi-X.root");

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
        TH1D* h_nh3[10];  // Array of histograms for each NH3 period
        // Initialize the histograms
        for (int i = 0; i < 10; ++i) {
            h_nh3[i] = new TH1D(Form("h_nh3_%d", i), "", 1, varMin, varMax);
        }
        TH1D *h_c = new TH1D("h_c", "", 1, varMin, varMax);
        TH1D *h_ch = new TH1D("h_ch", "", 1, varMin, varMax);
        TH1D *h_he = new TH1D("h_he", "", 1, varMin, varMax);
        TH1D *h_empty = new TH1D("h_empty", "", 1, varMin, varMax);

        double sumCurrentVariable = 0.0;
        int count = 0;

        // Helper function to fill histograms based on kinematic cuts and track mean
        auto fill_histogram = [&](TTreeReader& reader, TH1D* hist, BaseKinematicCuts* cuts, 
            bool is_nh3, int min_run = 0, int max_run = 1000000) {
            TTreeReaderValue<double> currentVariable(reader, propertyNames[currentFits].c_str());

            while (reader.Next()) {
                bool passedKinematicCuts = cuts->applyCuts(currentFits, false);
                if (is_nh3) {
                    if (*runnum < min_run || *runnum > max_run) continue;
                }
                if (*currentVariable >= varMin && *currentVariable < varMax && passedKinematicCuts) {
                    hist->Fill(*currentVariable);
                    sumCurrentVariable += *currentVariable;
                    ++count;

                }
            }
            reader.Restart();
        };

        // Call fill_histogram for each target type
        fill_histogram(nh3Reader, h_nh3[0], nh3Cuts, true, 0, 1000000);  // NH3 data
        fill_histogram(cReader, h_c,  cCuts, false, 0, 1000000);  // Carbon data
        fill_histogram(chReader, h_ch, chCuts, false, 0, 1000000); // CH2 data
        fill_histogram(heReader, h_he,  heCuts, false, 0, 1000000);  // Helium data
        fill_histogram(emptyReader, h_empty,  emptyCuts, false, 0, 1000000);  // Empty target data

        // Calculate the mean value of currentVariable in this bin
        double meanCurrentVariable = (count > 0) ? (sumCurrentVariable / count) : (varMin + varMax) / 2.0;

        // Retrieve bin contents
        double nA = h_nh3[0]->GetBinContent(1);
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

        // // Loop over each NH3 period to fill histograms and calculate dilution factors
        // for (int i = 0; i < nh3_periods.size(); ++i) {
        //     std::cout << "On period " << (i+1) << "." << std::endl;
        //     int min_run = nh3_periods[i].first;
        //     int max_run = nh3_periods[i].second;

        //     // Reset the variables for this NH3 period
        //     sumCurrentVariable = 0.0;
        //     count = 0;

        //     // Fill the histogram for the NH3 data only for the current period
        //     fill_histogram(nh3Reader, h_nh3[i+1], nh3Cuts, true, min_run, max_run);

        //     // Calculate the mean value of currentVariable in this bin
        //     double meanCurrentVariable = (count > 0) ? (sumCurrentVariable / count) : (varMin + varMax) / 2.0;

        //     // Retrieve bin contents for other target types (already filled outside loop)
        //     double nA = h_nh3[i+1]->GetBinContent(1);

        //     // Define variables to hold the fractional charge values
        //     double xA, xC, xCH, xHe, xf;

        //     // Use switch statement to assign the correct fractional charge values for each period
        //     switch (i + 1) {
        //         case 1:
        //             xA = xAperiod_1; xC = xCperiod_1; xCH = xCHperiod_1; xHe = xHeperiod_1; xf = xfperiod_1;
        //             break;
        //         case 2:
        //             xA = xAperiod_2; xC = xCperiod_2; xCH = xCHperiod_2; xHe = xHeperiod_2; xf = xfperiod_2;
        //             break;
        //         case 3:
        //             xA = xAperiod_3; xC = xCperiod_3; xCH = xCHperiod_3; xHe = xHeperiod_3; xf = xfperiod_3;
        //             break;
        //         case 4:
        //             xA = xAperiod_4; xC = xCperiod_4; xCH = xCHperiod_4; xHe = xHeperiod_4; xf = xfperiod_4;
        //             break;
        //         case 5:
        //             xA = xAperiod_5; xC = xCperiod_5; xCH = xCHperiod_5; xHe = xHeperiod_5; xf = xfperiod_5;
        //             break;
        //         case 6:
        //             xA = xAperiod_6; xC = xCperiod_6; xCH = xCHperiod_6; xHe = xHeperiod_6; xf = xfperiod_6;
        //             break;
        //         case 7:
        //             xA = xAperiod_7; xC = xCperiod_7; xCH = xCHperiod_7; xHe = xHeperiod_7; xf = xfperiod_7;
        //             break;
        //         case 8:
        //             xA = xAperiod_8; xC = xCperiod_8; xCH = xCHperiod_8; xHe = xHeperiod_8; xf = xfperiod_8;
        //             break;
        //         case 9:
        //             xA = xAperiod_9; xC = xCperiod_9; xCH = xCHperiod_9; xHe = xHeperiod_9; xf = xfperiod_9;
        //             break;
        //         default:
        //             // Use total values for index 0 or if something goes wrong
        //             xA = xAtotal; xC = xCtotal; xCH = xCHtotal; xHe = xHetotal; xf = xftotal;
        //             break;
        //     }

        //     // Calculate dilution factors for the current NH3 period
        //     auto [dilution, error] = calculate_dilution_and_error(nA, nC, nCH, nMT, nf, xA, xC, xCH, xHe, xf);

        //     // Add the dilution factor and error to the corresponding TGraphErrors
        //     gr_dilution[i + 1]->SetPoint(binIndex, meanCurrentVariable, dilution);
        //     gr_dilution[i + 1]->SetPointError(binIndex, 0, error);

        // }

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
    gr_dilution[0]->GetYaxis()->SetRangeUser(0.15, 0.35);
    gr_dilution[0]->SetMarkerStyle(20);
    gr_dilution[0]->SetMarkerColor(kBlack);
    gr_dilution[0]->Draw("AP");

    std::string outputDir = "output/dilution_factor_plots/";
    std::string outputFileName = outputDir + "df_" + binNames[currentFits] + "_" + prefix + ".pdf";
    canvas->SaveAs(outputFileName.c_str());


    //

    // // Create a new canvas for the period-specific dilution factors
    // TCanvas* canvas_periods = new TCanvas("c_dilution_periods", "Dilution Factor Plot by Period", 800, 600);
    // canvas_periods->SetLeftMargin(0.15);
    // canvas_periods->SetBottomMargin(0.15);

    // // Create a legend in the top right, increasing its size and adding a border
    // TLegend* legend = new TLegend(0.7, 0.5, 0.90, 0.90);  // Make the legend larger
    // legend->SetTextSize(0.04);  // Increase text size
    // legend->SetBorderSize(1);   // Add black border around the legend
    // legend->SetLineColor(kBlack);  // Set border color to black

    // // Define different marker styles and colors for each period
    // int markerStyles[9] = {20, 21, 22, 23, 24, 25, 26, 27, 28};  // Different marker styles
    // int markerColors[9] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet, kYellow+2, kPink+7};  // Different colors

    // // Vector to store the y-values for standard deviation calculation
    // std::vector<double> y_values;

    // for (int i = 0; i < 9; ++i) {
    //     gr_dilution[i + 1]->SetTitle("");  // No title
    //     gr_dilution[i + 1]->GetXaxis()->SetTitle(formatLabelName(prefix).c_str());
    //     gr_dilution[i + 1]->GetXaxis()->SetLimits(config.xMin, config.xMax*1.2);
    //     gr_dilution[i + 1]->GetXaxis()->SetTitleSize(0.05);
    //     gr_dilution[i + 1]->GetYaxis()->SetTitle("D_{f}");
    //     gr_dilution[i + 1]->GetYaxis()->SetTitleSize(0.05);
    //     gr_dilution[i + 1]->GetYaxis()->SetTitleOffset(1.6);
    //     gr_dilution[i + 1]->GetYaxis()->SetRangeUser(0.15, 0.30);
        
    //     // Set marker style and color
    //     gr_dilution[i + 1]->SetMarkerStyle(markerStyles[i]);
    //     gr_dilution[i + 1]->SetMarkerColor(markerColors[i]);
        
    //     // Add to legend
    //     legend->AddEntry(gr_dilution[i + 1], Form("Period %d", i + 1), "p");

    //     // Collect y-values for standard deviation calculation
    //     if (binNames[currentFits] == "integrated") {
    //         for (int j = 0; j < gr_dilution[i + 1]->GetN(); ++j) {
    //             double x, y;
    //             gr_dilution[i + 1]->GetPoint(j, x, y);  // Gets the value of the point, not the error
    //             y_values.push_back(y);  // Collect y-values
    //         }
    //     }
        
    //     // Draw on canvas (use "P same" for subsequent plots to overlay them)
    //     if (i == 0) {
    //         gr_dilution[i + 1]->Draw("AP");  // Draw the first graph
    //     } else {
    //         gr_dilution[i + 1]->Draw("P same");  // Overlay the rest
    //     }
    // }

    // // Draw the legend
    // legend->Draw();

    // // If prefix is "integrated", calculate the standard deviation and display it
    // if (binNames[currentFits] == "integrated") {
    //     double stddev = calculate_standard_deviation(y_values);

    //     // Create a TLatex object to draw the standard deviation on the canvas
    //     TLatex latex;
    //     latex.SetNDC();  // Set to normalized device coordinates
    //     latex.SetTextSize(0.04);  // Set text size
    //     latex.DrawLatex(0.23, 0.65, Form("Std Dev: %.4f", stddev));  // Display stddev in the top left
    // }

    // // Save the canvas
    // std::string outputFileNamePeriods = outputDir + "df_periods_" + binNames[currentFits] + "_" + prefix + ".pdf";
    // canvas_periods->SaveAs(outputFileNamePeriods.c_str());

    return dilutionResults;
}