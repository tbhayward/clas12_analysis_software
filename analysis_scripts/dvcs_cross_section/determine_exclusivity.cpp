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
    canvas->Divide(2, 4);

    std::vector<std::string> variables = {"open_angle_ep2", "Mx2_2", "theta_gamma_gamma", "placeholder", "Emiss2", "Mx2", "Mx2_1", "pTmiss"};

    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i] == "placeholder") {
            continue; // Skip placeholder for now
        }

        // Retrieve histogram bin settings
        HistConfig config = histConfigs[variables[i]];

        // Create histograms for data and MC
        TH1D* hData = new TH1D(("data_" + variables[i]).c_str(), formatLabelName(variables[i]).c_str(), config.bins, config.min, config.max);
        TH1D* hMC = new TH1D(("mc_" + variables[i]).c_str(), formatLabelName(variables[i]).c_str(), config.bins, config.min, config.max);

        // Fill histograms using TTreeReader values
        TTreeReaderValue<double> dataVar(dataReader, variables[i].c_str());
        TTreeReaderValue<double> mcVar(mcReader, variables[i].c_str());

        while (dataReader.Next()) {
            hData->Fill(*dataVar);
        }
        while (mcReader.Next()) {
            hMC->Fill(*mcVar);
        }

        // Draw the histograms
        canvas->cd(i + 1);
        hData->SetLineColor(kBlue);
        hMC->SetLineColor(kRed);
        hData->Draw("HIST");
        hMC->Draw("HIST SAME");

        // Add a legend
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(hData, "Data", "l");
        legend->AddEntry(hMC, "MC", "l");
        legend->Draw();
    }

    // Save the canvas
    canvas->SaveAs((outputDir + "/exclusivity_plots_rga_fa18_inb.png").c_str());

    // Clean up
    delete canvas;
}