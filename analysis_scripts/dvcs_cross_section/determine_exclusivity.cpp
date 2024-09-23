#include "determine_exclusivity.h"
#include "histConfigs.h"
#include "formatLabelName.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void determine_exclusivity(TTreeReader& dataReader, TTreeReader& mcReader, const std::string& output_dir) {
    // Set up the canvas with 2 rows and 4 columns
    TCanvas* canvas = new TCanvas("canvas", "Exclusivity Variables", 1600, 800);
    canvas->Divide(4, 2);  // 4 columns, 2 rows

    // Disable stat boxes
    gStyle->SetOptStat(0);

    // Variables to plot
    std::vector<std::string> variables = {"open_angle_ep2", "Mx2_2", "theta_gamma_gamma", "placeholder", 
                                          "Emiss2", "Mx2", "Mx2_1", "pTmiss"};

    // Create TTreeReaderValue objects for each variable for both data and MC
    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i] == "placeholder") continue;  // Skip placeholder variable for now

        // Restart the readers to avoid errors
        dataReader.Restart();
        mcReader.Restart();

        // Set up TTreeReaderValues for data and MC
        TTreeReaderValue<Double_t> dataValue(dataReader, variables[i].c_str());
        TTreeReaderValue<Double_t> mcValue(mcReader, variables[i].c_str());

        // Create histograms for both data and MC
        std::string hist_data_name = "hist_data_" + variables[i];
        std::string hist_mc_name = "hist_mc_" + variables[i];
        const HistConfig& config = histConfigs[variables[i]];

        TH1D* hist_data = new TH1D(hist_data_name.c_str(), hist_data_name.c_str(), config.nBins, config.xMin, config.xMax);
        TH1D* hist_mc = new TH1D(hist_mc_name.c_str(), hist_mc_name.c_str(), config.nBins, config.xMin, config.xMax);

        // Fill data histograms
        while (dataReader.Next()) {
            hist_data->Fill(*dataValue);
        }

        // Fill MC histograms
        while (mcReader.Next()) {
            hist_mc->Fill(*mcValue);
        }

        // Normalize histograms
        hist_data->Scale(1.0 / hist_data->Integral());
        hist_mc->Scale(1.0 / hist_mc->Integral());

        // Select the canvas pad
        canvas->cd(i + 1);

        // Draw histograms
        hist_data->SetLineColor(kBlue);
        hist_mc->SetLineColor(kRed);
        hist_data->Draw("HIST");
        hist_mc->Draw("HIST SAME");

        // Set axis labels using formatLabelName
        hist_data->GetXaxis()->SetTitle(formatLabelName(variables[i]).c_str());
        hist_data->GetYaxis()->SetTitle("normalized counts");

        // Create legend
        int data_counts = hist_data->GetEntries();
        int mc_counts = hist_mc->GetEntries();
        TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
        legend->AddEntry(hist_data, ("Data, " + std::to_string(data_counts)).c_str(), "l");
        legend->AddEntry(hist_mc, ("MC, " + std::to_string(mc_counts)).c_str(), "l");
        legend->Draw();
    }

    // Save the canvas to the output directory
    std::string output_file = output_dir + "/exclusivity_plots_rga_fa18_inb.png";
    canvas->SaveAs(output_file.c_str());

    std::cout << "Plots saved to " << output_file << std::endl;

    // Clean up
    delete canvas;
}