#include "determine_exclusivity.h"
#include "histConfigs.h"
#include "formatLabelName.h"
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

    // Create histograms for each detector configuration (FD,FD), (CD,FD), (CD,FT)
    std::map<std::string, std::vector<TH1D*>> histograms_data;
    std::map<std::string, std::vector<TH1D*>> histograms_mc;

    std::vector<std::string> variables = {"open_angle_ep2", "Mx2_2", "theta_gamma_gamma", "Emiss2", "Mx2", "Mx2_1", "pTmiss"};
    std::vector<std::string> configurations = {"FD_FD", "CD_FD", "CD_FT"};
    
    for (const auto& config : configurations) {
        for (const auto& var : variables) {
            // Initialize histograms for each configuration and variable
            histograms_data[config].push_back(new TH1D(("data_" + config + "_" + var).c_str(), "", histConfigs[var].bins, histConfigs[var].min, histConfigs[var].max));
            histograms_mc[config].push_back(new TH1D(("mc_" + config + "_" + var).c_str(), "", histConfigs[var].bins, histConfigs[var].min, histConfigs[var].max));
        }
    }

    // Normalization condition: theta between 14 and 18 degrees in radians
    const double theta_min = 14.0 * 3.14159 / 180.0;  // Convert 14 degrees to radians
    const double theta_max = 18.0 * 3.14159 / 180.0;  // Convert 18 degrees to radians

    // Readers for e_theta variable (detector1 and detector2 are commented out for now)
    TTreeReaderValue<double> eTheta_data(dataReader, "e_theta");
    TTreeReaderValue<double> eTheta_mc(mcReader, "e_theta");

    // Commenting out detector1 and detector2 logic temporarily
    /*
    TTreeReaderValue<int> detector1_data(dataReader, "detector1");
    TTreeReaderValue<int> detector2_data(dataReader, "detector2");
    TTreeReaderValue<int> detector1_mc(mcReader, "detector1");
    TTreeReaderValue<int> detector2_mc(mcReader, "detector2");
    */

    // Total electron counts for normalization
    std::map<std::string, int> total_electrons_data = {{"FD_FD", 0}, {"CD_FD", 0}, {"CD_FT", 0}};
    std::map<std::string, int> total_electrons_mc = {{"FD_FD", 0}, {"CD_FD", 0}, {"CD_FT", 0}};

    // Fill histograms for all configurations (for now, we fill all three configurations for every event)
    while (dataReader.Next()) {
        for (const auto& config_key : configurations) {
            // Normalization for electrons in the specified theta range
            if (*eTheta_data >= theta_min && *eTheta_data <= theta_max) {
                ++total_electrons_data[config_key];
            }

            // Fill histograms for variables
            for (size_t i = 0; i < variables.size(); ++i) {
                TTreeReaderValue<double> var_reader(dataReader, variables[i].c_str());
                histograms_data[config_key][i]->Fill(*var_reader);
            }
        }
    }
    dataReader.Restart();  // Restart reader for MC

    while (mcReader.Next()) {
        for (const auto& config_key : configurations) {
            // Normalization for electrons in the specified theta range
            if (*eTheta_mc >= theta_min && *eTheta_mc <= theta_max) {
                ++total_electrons_mc[config_key];
            }

            // Fill histograms for variables
            for (size_t i = 0; i < variables.size(); ++i) {
                TTreeReaderValue<double> var_reader(mcReader, variables[i].c_str());
                histograms_mc[config_key][i]->Fill(*var_reader);
            }
        }
    }
    mcReader.Restart();  // Restart reader

    // Now process each configuration and create plots
    for (const auto& config : configurations) {
        TCanvas* canvas = new TCanvas(("canvas_" + config).c_str(), ("Exclusivity Plots: " + config).c_str(), 1600, 800);
        canvas->Divide(4, 2);  // 2 rows, 4 columns

        for (size_t i = 0; i < variables.size(); ++i) {
            canvas->cd(i + 1);
            gPad->SetLeftMargin(0.15);  // Add left padding
            gPad->SetBottomMargin(0.15);  // Add bottom padding
            
            // Normalize histograms based on total electrons
            if (total_electrons_data[config] > 0)
                histograms_data[config][i]->Scale(1.0 / total_electrons_data[config]);
            if (total_electrons_mc[config] > 0)
                histograms_mc[config][i]->Scale(1.0 / total_electrons_mc[config]);

            // Set the plotting styles for points with error bars (no horizontal errors)
            histograms_data[config][i]->SetMarkerStyle(20);
            histograms_data[config][i]->SetMarkerColor(kBlue);
            histograms_mc[config][i]->SetMarkerStyle(24);
            histograms_mc[config][i]->SetMarkerColor(kRed);

            // Set axis labels
            histograms_data[config][i]->SetXTitle(formatLabelName(variables[i]).c_str());
            histograms_data[config][i]->SetYTitle("Normalized counts");

            // Set Y-axis range
            double max_data = histograms_data[config][i]->GetMaximum();
            double max_mc = histograms_mc[config][i]->GetMaximum();
            double y_max = 1.35 * std::max(max_data, max_mc);  // Set y-axis range from 0 to 1.35 * max

            histograms_data[config][i]->GetYaxis()->SetRangeUser(0, y_max);

            // Draw histograms with error bars
            histograms_data[config][i]->Draw("E1");  // Points with error bars
            histograms_mc[config][i]->Draw("E1 SAME");

            // Add a legend with blue/red text for Data/MC
            TLegend* legend = new TLegend(0.375, 0.7, 0.9, 0.9);
            legend->AddEntry(histograms_data[config][i], ("#color[4]{Data} (" + std::to_string(static_cast<int>(histograms_data[config][i]->GetEntries())) + " counts)").c_str(), "p");
            legend->AddEntry(histograms_mc[config][i], ("#color[2]{MC} (" + std::to_string(static_cast<int>(histograms_mc[config][i]->GetEntries())) + " counts)").c_str(), "p");
            legend->Draw();

            // Add the title for each subplot
            histograms_data[config][i]->SetTitle((plotTitle + "; " + config).c_str());
        }

        // Save the canvas for each configuration
        canvas->SaveAs((outputDir + "/exclusivity_plots_" + config + ".png").c_str());

        // Clean up
        delete canvas;
    }

    // Clean up histograms
    for (const auto& config : configurations) {
        for (auto hist : histograms_data[config]) delete hist;
        for (auto hist : histograms_mc[config]) delete hist;
    }
}