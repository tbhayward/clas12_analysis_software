#include "determine_exclusivity.h"
#include "histConfigs.h"
#include "formatLabelName.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>

void determine_exclusivity(TTreeReader& dataReader, TTreeReader& mcReader, const std::string& outputDir) {
    // Create a 2x4 canvas
    TCanvas* canvas = new TCanvas("exclusivity_canvas", "Exclusivity Plots", 1600, 800);
    canvas->Divide(4, 2);  // 2 rows, 4 columns

    std::vector<std::string> variables = {"open_angle_ep2", "Mx2_2", "theta_gamma_gamma", "placeholder", "Emiss2", "Mx2", "Mx2_1", "pTmiss"};

    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i] == "placeholder") {
            continue; // Skip placeholder for now
        }

        // Retrieve histogram bin settings
        HistConfig config = histConfigs[variables[i]];

        // Create histogram names
        std::string hist_data_name = "data_" + variables[i];
        std::string hist_mc_name = "mc_" + variables[i];

        // Create histograms for data and MC
        TH1D* hist_data = new TH1D(hist_data_name.c_str(), formatLabelName(variables[i]).c_str(), config.bins, config.min, config.max);
        TH1D* hist_mc = new TH1D(hist_mc_name.c_str(), formatLabelName(variables[i]).c_str(), config.bins, config.min, config.max);

        // Fill histograms using TTreeReader values
        TTreeReaderValue<double> dataVar(dataReader, variables[i].c_str());
        TTreeReaderValue<double> mcVar(mcReader, variables[i].c_str());

        while (dataReader.Next()) {
            hist_data->Fill(*dataVar);
        }
        while (mcReader.Next()) {
            hist_mc->Fill(*mcVar);
        }

        // Normalize histograms
        hist_data->Scale(1.0 / hist_data->Integral());
        hist_mc->Scale(1.0 / hist_mc->Integral());

        // Draw the histograms
        canvas->cd(i + 1);
        hist_data->SetLineColor(kBlue);
        hist_mc->SetLineColor(kRed);
        hist_data->Draw("HIST");
        hist_mc->Draw("HIST SAME");

        // Add a legend
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