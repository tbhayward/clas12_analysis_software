#include "determine_exclusivity.h"
#include "histConfigs.h"
#include "formatLabelName.h"
#include "kinematic_cuts.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>
#include <cmath>  // for radian conversion

void determine_exclusivity(TTreeReader& dataReader, TTreeReader& mcReader, const std::string& outputDir, const std::string& plotTitle) {
    // Set up global style options (remove stat boxes)
    gStyle->SetOptStat(0);

    // Increase font size for axis labels and title
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetLegendTextSize(0.04);

    // Create a 2x4 canvas for original plots
    TCanvas* canvas = new TCanvas("exclusivity_canvas", "Exclusivity Plots", 1600, 800);
    canvas->Divide(4, 2);  // 2 rows, 4 columns

    // Create a 2x4 canvas for plots with "Loose Cuts"
    TCanvas* canvas_loose_cuts = new TCanvas("exclusivity_canvas_loose_cuts", "Exclusivity Plots with Loose Cuts", 1600, 800);
    canvas_loose_cuts->Divide(4, 2);

    // Vector of variables for plotting
    std::vector<std::string> variables = {"open_angle_ep2", "Mx2_2", "theta_gamma_gamma", "placeholder", "Emiss2", "Mx2", "Mx2_1", "pTmiss"};

    // Readers for e_theta and other relevant variables for cuts
    TTreeReaderValue<double> t_data(dataReader, "t");
    TTreeReaderValue<double> t_mc(mcReader, "t");
    TTreeReaderValue<double> open_angle_ep2_data(dataReader, "open_angle_ep2");
    TTreeReaderValue<double> open_angle_ep2_mc(mcReader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_data(dataReader, "Emiss2");
    TTreeReaderValue<double> Emiss2_mc(mcReader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(dataReader, "Mx2_1");
    TTreeReaderValue<double> Mx2_1_mc(mcReader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(dataReader, "pTmiss");
    TTreeReaderValue<double> pTmiss_mc(mcReader, "pTmiss");

    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i] == "placeholder") {
            continue;  // Skip placeholder for now
        }

        // Retrieve histogram bin settings
        HistConfig config = histConfigs[variables[i]];

        // Create histogram names
        std::string hist_data_name = "data_" + variables[i];
        std::string hist_mc_name = "mc_" + variables[i];
        std::string hist_data_name_loose = "data_loose_" + variables[i];
        std::string hist_mc_name_loose = "mc_loose_" + variables[i];

        // Create histograms for data and MC
        TH1D* hist_data = new TH1D(hist_data_name.c_str(), "", config.bins, config.min, config.max);
        TH1D* hist_mc = new TH1D(hist_mc_name.c_str(), "", config.bins, config.min, config.max);
        TH1D* hist_data_loose = new TH1D(hist_data_name_loose.c_str(), "", config.bins, config.min, config.max);
        TH1D* hist_mc_loose = new TH1D(hist_mc_name_loose.c_str(), "", config.bins, config.min, config.max);

        // Create readers for the variable
        TTreeReaderValue<double> dataVar(dataReader, variables[i].c_str());
        TTreeReaderValue<double> mcVar(mcReader, variables[i].c_str());

        // Fill data histograms
        while (dataReader.Next()) {
            hist_data->Fill(*dataVar);
            // Apply kinematic cuts and fill the "Loose Cuts" histograms if they pass
            if (apply_kinematic_cuts(*t_data, *open_angle_ep2_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {
                hist_data_loose->Fill(*dataVar);
            }
        }
        dataReader.Restart();  // Restart reader for the next variable

        // Fill MC histograms
        while (mcReader.Next()) {
            hist_mc->Fill(*mcVar);
            // Apply kinematic cuts and fill the "Loose Cuts" histograms if they pass
            if (apply_kinematic_cuts(*t_mc, *open_angle_ep2_mc, *Emiss2_mc, *Mx2_1_mc, *pTmiss_mc)) {
                hist_mc_loose->Fill(*mcVar);
            }
        }
        mcReader.Restart();  // Restart reader for the next variable

        // Normalize histograms based on their integrals
        if (hist_data->Integral() != 0) { hist_data->Scale(1.0 / hist_data->Integral()); }
        if (hist_mc->Integral() != 0) { hist_mc->Scale(1.0 / hist_mc->Integral()); }
        if (hist_data_loose->Integral() != 0) { hist_data_loose->Scale(1.0 / hist_data_loose->Integral()); }
        if (hist_mc_loose->Integral() != 0) { hist_mc_loose->Scale(1.0 / hist_mc_loose->Integral()); }

        // Get the maximum value for setting y-axis range
        double max_data = hist_data->GetMaximum();
        double max_mc = hist_mc->GetMaximum();
        double max_data_loose = hist_data_loose->GetMaximum();
        double max_mc_loose = hist_mc_loose->GetMaximum();
        double y_max = 1.35 * std::max({max_data, max_mc});
        double y_max_loose = 1.35 * std::max({max_data_loose, max_mc_loose});

        // Use const char* for the title
        const char* originalTitle = plotTitle.c_str();
        // Set the title for the histograms before drawing them
        hist_data->SetTitle(originalTitle);
        // Set the title for the histograms before drawing them
        std::string looseCutsTitle = plotTitle + " ; Loose Cuts";
        hist_data->SetTitle(plotTitle.c_str());
        hist_data_loose->SetTitle(looseCutsTitle.c_str());  // Set loose cuts title properly


        // Draw the histograms for original plots
        canvas->cd(i + 1);
        gPad->SetLeftMargin(0.15);  // Add left padding
        gPad->SetBottomMargin(0.15);  // Add bottom padding
        hist_data->SetLineColor(kBlue);
        hist_mc->SetLineColor(kRed);
        hist_data->SetXTitle(formatLabelName(variables[i]).c_str());
        hist_data->SetYTitle("Normalized counts");
        hist_data->GetYaxis()->SetRangeUser(0, y_max);
        hist_data->Draw("HIST");
        hist_mc->Draw("HIST SAME");

        // Add a legend with the count information (integer format)
        TLegend* legend = new TLegend(0.375, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_data, ("Data (" + std::to_string(static_cast<int>(hist_data->GetEntries())) + " counts)").c_str(), "l");
        legend->AddEntry(hist_mc, ("MC (" + std::to_string(static_cast<int>(hist_mc->GetEntries())) + " counts)").c_str(), "l");
        legend->Draw();

        // Draw the histograms for "Loose Cuts" plots
        canvas_loose_cuts->cd(i + 1);
        gPad->SetLeftMargin(0.15);  // Add left padding
        gPad->SetBottomMargin(0.15);  // Add bottom padding
        hist_data_loose->SetLineColor(kBlue);
        hist_mc_loose->SetLineColor(kRed);
        hist_data_loose->SetXTitle(formatLabelName(variables[i]).c_str());
        hist_data_loose->SetYTitle("Normalized counts");
        hist_data_loose->GetYaxis()->SetRangeUser(0, y_max_loose);
        hist_data_loose->Draw("HIST");
        hist_mc_loose->Draw("HIST SAME");

        // Add a legend with the count information for "Loose Cuts" (integer format)
        TLegend* legend_loose = new TLegend(0.25, 0.7, 0.9, 0.9);
        legend_loose->AddEntry(hist_data_loose, ("Data (" + std::to_string(static_cast<int>(hist_data_loose->GetEntries())) + " counts; Loose Cuts)").c_str(), "l");
        legend_loose->AddEntry(hist_mc_loose, ("MC (" + std::to_string(static_cast<int>(hist_mc_loose->GetEntries())) + " counts; Loose Cuts)").c_str(), "l");
        legend_loose->Draw();
    }

    // Save the canvases with updated names to reflect the plotTitle input
    std::string cleanTitle = plotTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');  // Replace spaces with underscores
    canvas->SaveAs((outputDir + "/exclusivity_plots_" + cleanTitle + ".png").c_str());
    canvas_loose_cuts->SaveAs((outputDir + "/exclusivity_plots_" + cleanTitle + "_loose_cuts.png").c_str());

    // Clean up
    delete canvas;
    delete canvas_loose_cuts;
}