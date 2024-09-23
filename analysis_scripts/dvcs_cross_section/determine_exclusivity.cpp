#include "determine_exclusivity.h"
#include "histConfigs.h"
#include "formatLabelName.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>
#include <cmath>  // for radian conversion

void determine_exclusivity(TTreeReader& dataReader, TTreeReader& mcReader, const std::string& outputDir) {
    // Set up global style options (remove stat boxes)
    gStyle->SetOptStat(0);

    // Create a 2x4 canvas
    TCanvas* canvas = new TCanvas("exclusivity_canvas", "Exclusivity Plots", 1600, 800);
    canvas->Divide(4, 2);  // 2 rows, 4 columns

    // Vector of variables for plotting
    std::vector<std::string> variables = {"open_angle_ep2", "Mx2_2", "theta_gamma_gamma", "placeholder", "Emiss2", "Mx2", "Mx2_1", "pTmiss"};

    // Normalization condition: theta between 14 and 18 degrees in radians
    const double theta_min = 14.0 * TMath::Pi() / 180.0;  // 14 degrees in radians
    const double theta_max = 18.0 * TMath::Pi() / 180.0;  // 18 degrees in radians

    // Readers for e_theta variable for normalization
    TTreeReaderValue<double> eTheta_data(dataReader, "e_theta");
    TTreeReaderValue<double> eTheta_mc(mcReader, "e_theta");

    int total_electrons_data = 0;
    int total_electrons_mc = 0;

    // Count total electrons in data (within 14 to 18 degrees)
    while (dataReader.Next()) {
        if (*eTheta_data >= theta_min && *eTheta_data <= theta_max) {
            ++total_electrons_data;
        }
    }
    dataReader.Restart();  // Restart the reader after counting

    // Count total electrons in MC (within 14 to 18 degrees)
    while (mcReader.Next()) {
        if (*eTheta_mc >= theta_min && *eTheta_mc <= theta_max) {
            ++total_electrons_mc;
        }
    }
    mcReader.Restart();  // Restart the reader after counting

    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i] == "placeholder") {
            continue;  // Skip placeholder for now
        }

        // Retrieve histogram bin settings
        HistConfig config = histConfigs[variables[i]];

        // Create histogram names
        std::string hist_data_name = "data_" + variables[i];
        std::string hist_mc_name = "mc_" + variables[i];

        // Create histograms for data and MC
        TH1D* hist_data = new TH1D(hist_data_name.c_str(), formatLabelName(variables[i]).c_str(), config.bins, config.min, config.max);
        TH1D* hist_mc = new TH1D(hist_mc_name.c_str(), formatLabelName(variables[i]).c_str(), config.bins, config.min, config.max);

        // Create readers for the variable
        TTreeReaderValue<double> dataVar(dataReader, variables[i].c_str());
        TTreeReaderValue<double> mcVar(mcReader, variables[i].c_str());

        // Fill data histograms
        while (dataReader.Next()) {
            hist_data->Fill(*dataVar);
        }
        dataReader.Restart();  // Restart reader for the next variable

        // Fill MC histograms
        while (mcReader.Next()) {
            hist_mc->Fill(*mcVar);
        }
        mcReader.Restart();  // Restart reader for the next variable

        // Normalize histograms based on the total electron count
        hist_data->Scale(1.0 / total_electrons_data);
        hist_mc->Scale(1.0 / total_electrons_mc);

        // Draw the histograms
        canvas->cd(i + 1);
        hist_data->SetLineColor(kBlue);
        hist_mc->SetLineColor(kRed);
        hist_data->SetXTitle(formatLabelName(variables[i]).c_str());
        hist_data->SetYTitle("Normalized counts");
        hist_data->Draw("HIST");
        hist_mc->Draw("HIST SAME");

        // Add a legend with the count information
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_data, ("Data (" + std::to_string(hist_data->GetEntries()) + ")").c_str(), "l");
        legend->AddEntry(hist_mc, ("MC (" + std::to_string(hist_mc->GetEntries()) + ")").c_str(), "l");
        legend->Draw();
    }

    // Save the canvas
    canvas->SaveAs((outputDir + "/exclusivity_plots_rga_fa18_inb.png").c_str());

    // Clean up
    delete canvas;
}