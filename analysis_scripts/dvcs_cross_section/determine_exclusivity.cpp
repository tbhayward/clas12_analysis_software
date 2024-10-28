#include "determine_exclusivity.h"
#include "histConfigs.h"
#include "formatLabelName.h"
#include "kinematic_cuts.h"
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>
#include <TF1.h>
#include <filesystem>
#include <cmath>  // for radian conversion

// Define the asymmetric Gaussian function
Double_t asymmetricGaussian(Double_t* x, Double_t* par) {
    // Parameters:
    // par[0] = Normalization
    // par[1] = Mean (mu)
    // par[2] = Sigma_left (sigma for x < mu)
    // par[3] = Sigma_right (sigma for x >= mu)
    Double_t xx = x[0];
    Double_t mu = par[1];
    Double_t sigma = (xx < mu) ? par[2] : par[3];
    Double_t norm = par[0];
    return norm * exp(-0.5 * ((xx - mu) / sigma) * ((xx - mu) / sigma));
}

void determine_exclusivity(const std::string& analysisType, const std::string& topology, TTreeReader& dataReader, TTreeReader& mcReader, const std::string& outputDir, const std::string& plotTitle) {
    // Set up global style options (remove stat boxes)
    gStyle->SetOptStat(0);

    // Increase font size for axis labels and title
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetLabelSize(0.04, "XY");
    gStyle->SetLegendTextSize(0.04);

    // Create the correct subdirectory based on the analysis type (dvcs or eppi0)
    std::string analysis_dir = (analysisType == "dvcs") ? "dvcs" : "eppi0";
    std::string final_output_dir = outputDir + "/exclusivity_plots/" + analysis_dir;

    // Check and create the necessary directory if it doesn't exist
    if (!std::filesystem::exists(final_output_dir)) {
        std::filesystem::create_directories(final_output_dir);
    }

    // Create a 2x4 canvas for original plots
    TCanvas* canvas = new TCanvas("exclusivity_canvas", "Exclusivity Plots", 1600, 800);
    canvas->Divide(4, 2);  // 2 rows, 4 columns

    // Create a 2x4 canvas for plots with "Loose Cuts"
    TCanvas* canvas_loose_cuts = new TCanvas("exclusivity_canvas_loose_cuts", "Exclusivity Plots with Loose Cuts", 1600, 800);
    canvas_loose_cuts->Divide(4, 2);

    // Vector of variables for plotting (switch between "dvcs" and "eppi0")
    std::vector<std::string> variables;
    if (analysisType == "dvcs") {
        variables = {"open_angle_ep2", "theta_gamma_gamma", "pTmiss", "xF", "Emiss2", "Mx2", "Mx2_1", "Mx2_2"};
    } else if (analysisType == "eppi0") {
        variables = {"open_angle_ep2", "theta_pi0_pi0", "pTmiss", "xF", "Emiss2", "Mx2", "Mx2_1", "Mx2_2"};
    } else {
        throw std::runtime_error("Invalid analysis type! Must be 'dvcs' or 'eppi0'");
    }

    // Variables to perform fits on
    std::vector<std::string> variables_to_fit_gaussian = {"xF", "Emiss2"};
    std::vector<std::string> variables_to_fit_asymmetric_gaussian = {"pTmiss", "Mx2", "Mx2_1", "Mx2_2"};

    // Readers for relevant variables for cuts
    TTreeReaderValue<double> t_data(dataReader, "t1");
    TTreeReaderValue<double> t_mc(mcReader, "t1");
    TTreeReaderValue<double> open_angle_ep2_data(dataReader, "open_angle_ep2");
    TTreeReaderValue<double> open_angle_ep2_mc(mcReader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_data(dataReader, "Emiss2");
    TTreeReaderValue<double> Emiss2_mc(mcReader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(dataReader, "Mx2_1");
    TTreeReaderValue<double> Mx2_1_mc(mcReader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(dataReader, "pTmiss");
    TTreeReaderValue<double> pTmiss_mc(mcReader, "pTmiss");

    // Readers for detector status variables (detector1 and detector2)
    TTreeReaderValue<int> detector1_data(dataReader, "detector1");
    TTreeReaderValue<int> detector2_data(dataReader, "detector2");
    TTreeReaderValue<int> detector1_mc(mcReader, "detector1");
    TTreeReaderValue<int> detector2_mc(mcReader, "detector2");

    // Conditional initialization of "theta_neutral_neutral" based on analysis type
    TTreeReaderValue<double>* theta_neutral_neutral_data;
    TTreeReaderValue<double>* theta_neutral_neutral_mc;
    if (analysisType == "dvcs") {
        theta_neutral_neutral_data = new TTreeReaderValue<double>(dataReader, "theta_gamma_gamma");
        theta_neutral_neutral_mc = new TTreeReaderValue<double>(mcReader, "theta_gamma_gamma");
    } else if (analysisType == "eppi0") {
        theta_neutral_neutral_data = new TTreeReaderValue<double>(dataReader, "theta_pi0_pi0");
        theta_neutral_neutral_mc = new TTreeReaderValue<double>(mcReader, "theta_pi0_pi0");
    } 

    for (size_t i = 0; i < variables.size(); ++i) {

        // Retrieve histogram bin settings
        HistConfig config = histConfigs[variables[i]];

        // Create histogram names
        std::string hist_data_name = "data_" + variables[i];
        std::string hist_mc_name = "mc_" + variables[i];
        std::string hist_data_name_loose = "data_loose_" + variables[i];
        std::string hist_mc_name_loose = "mc_loose_" + variables[i];

        // Create histograms for data and MC
        TH1D* hist_data = new TH1D(hist_data_name.c_str(), "", config.bins, config.min, config.max);
        hist_data->SetDirectory(0);  // Prevent ROOT from storing histograms globally
        TH1D* hist_mc = new TH1D(hist_mc_name.c_str(), "", config.bins, config.min, config.max);
        hist_mc->SetDirectory(0);
        TH1D* hist_data_loose = new TH1D(hist_data_name_loose.c_str(), "", config.bins, config.min, config.max);
        hist_data_loose->SetDirectory(0);
        TH1D* hist_mc_loose = new TH1D(hist_mc_name_loose.c_str(), "", config.bins, config.min, config.max);
        hist_mc_loose->SetDirectory(0);

        // Set marker styles and colors for data and MC histograms
        hist_data->SetMarkerStyle(20);
        hist_data->SetMarkerColor(kBlue);
        hist_data->SetLineColor(kBlue);
        hist_data->SetMarkerSize(0.5);  // Make points smaller

        hist_mc->SetMarkerStyle(21);
        hist_mc->SetMarkerColor(kRed);
        hist_mc->SetLineColor(kRed);
        hist_mc->SetMarkerSize(0.5);

        hist_data_loose->SetMarkerStyle(20);
        hist_data_loose->SetMarkerColor(kBlue);
        hist_data_loose->SetLineColor(kBlue);
        hist_data_loose->SetMarkerSize(0.5);

        hist_mc_loose->SetMarkerStyle(21);
        hist_mc_loose->SetMarkerColor(kRed);
        hist_mc_loose->SetLineColor(kRed);
        hist_mc_loose->SetMarkerSize(0.5);

        // Create readers for the variable
        TTreeReaderValue<double> dataVar(dataReader, variables[i].c_str());
        TTreeReaderValue<double> mcVar(mcReader, variables[i].c_str());

        // Check detector status and fill histograms based on topology
        while (dataReader.Next()) {
            if ((topology == "(FD,FD)" && *detector1_data == 1 && *detector2_data == 1) ||
                (topology == "(CD,FD)" && *detector1_data == 2 && *detector2_data == 1) ||
                (topology == "(CD,FT)" && *detector1_data == 2 && *detector2_data == 0)) {

                hist_data->Fill(*dataVar);            
                if (apply_kinematic_cuts(*t_data, *open_angle_ep2_data, **theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {
                    hist_data_loose->Fill(*dataVar);
                }
            }
        }
        dataReader.Restart();  // Restart reader for the next variable

        while (mcReader.Next()) {
            if ((topology == "(FD,FD)" && *detector1_mc == 1 && *detector2_mc == 1) ||
                (topology == "(CD,FD)" && *detector1_mc == 2 && *detector2_mc == 1) ||
                (topology == "(CD,FT)" && *detector1_mc == 2 && *detector2_mc == 0)) {

                hist_mc->Fill(*mcVar);
                if (apply_kinematic_cuts(*t_mc, *open_angle_ep2_mc, **theta_neutral_neutral_mc, *Emiss2_mc, *Mx2_1_mc, *pTmiss_mc)) {
                    hist_mc_loose->Fill(*mcVar);
                }
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
        double y_max = 1.75 * std::max({max_data, max_mc});
        double y_max_loose = 1.75 * std::max({max_data_loose, max_mc_loose});

        hist_data->SetTitle(plotTitle.c_str());
        hist_data_loose->SetTitle(plotTitle.c_str());  

        // Draw the histograms for original plots
        canvas->cd(i + 1);
        gPad->SetLeftMargin(0.15);  // Add left padding
        gPad->SetBottomMargin(0.15);  // Add bottom padding
        hist_data->SetXTitle(formatLabelName(variables[i], analysisType).c_str());  // Pass analysisType to formatLabelName
        hist_data->SetYTitle("Normalized counts");
        hist_data->GetYaxis()->SetRangeUser(0, y_max);
        hist_data->Draw("E1");
        hist_mc->Draw("E1 SAME");

        // Add a legend with the count information (integer format)
        TLegend* legend = new TLegend(0.4, 0.7, 0.9, 0.9);
        legend->AddEntry(hist_data, ("Data (" + std::to_string(static_cast<int>(hist_data->GetEntries())) + " events)").c_str(), "lep");
        legend->AddEntry(hist_mc, ("MC (" + std::to_string(static_cast<int>(hist_mc->GetEntries())) + " events)").c_str(), "lep");
        legend->SetTextSize(0.03);  // Set smaller font size for legend
        legend->Draw();

        // Draw the histograms for "Loose Cuts" plots
        canvas_loose_cuts->cd(i + 1);
        gPad->SetLeftMargin(0.15);  // Add left padding
        gPad->SetBottomMargin(0.15);  // Add bottom padding
        hist_data_loose->SetXTitle(formatLabelName(variables[i], analysisType).c_str());  // Pass analysisType to formatLabelName
        hist_data_loose->SetYTitle("Normalized counts");
        hist_data_loose->GetYaxis()->SetRangeUser(0, y_max_loose);
        hist_data_loose->Draw("E1");
        hist_mc_loose->Draw("E1 SAME");

        // Check if the variable is in the list of variables to fit
        if (std::find(variables_to_fit_gaussian.begin(), variables_to_fit_gaussian.end(), variables[i]) != variables_to_fit_gaussian.end() ||
            std::find(variables_to_fit_asymmetric_gaussian.begin(), variables_to_fit_asymmetric_gaussian.end(), variables[i]) != variables_to_fit_asymmetric_gaussian.end()) {

            TF1* fit_data;
            TF1* fit_mc;

            // Determine the fit function based on the variable
            if (std::find(variables_to_fit_asymmetric_gaussian.begin(), variables_to_fit_asymmetric_gaussian.end(), variables[i]) != variables_to_fit_asymmetric_gaussian.end()) {
                // Use Asymmetric Gaussian function
                fit_data = new TF1("fit_data", asymmetricGaussian, hist_data_loose->GetXaxis()->GetXmin(), hist_data_loose->GetXaxis()->GetXmax(), 5);
                fit_mc = new TF1("fit_mc", asymmetricGaussian, hist_mc_loose->GetXaxis()->GetXmin(), hist_mc_loose->GetXaxis()->GetXmax(), 5);

                // Set initial parameter guesses
                double peak_data = hist_data_loose->GetMaximum();
                double peak_mc = hist_mc_loose->GetMaximum();
                double mean_guess = hist_data_loose->GetMean();

                fit_data->SetParameters(peak_data, mean_guess, 0.01, 0.01);
                fit_mc->SetParameters(peak_mc, mean_guess, 0.01, 0.01);

                // Optional: Set parameter limits if needed
                fit_data->SetParLimits(2, 0.05, 1);  // sigma_left
                fit_data->SetParLimits(3, 0.05, 1);  // sigma_right
                // Similar for fit_mc
            } else {
                // Use Gaussian function
                fit_data = new TF1("fit_data", "gaus", hist_data_loose->GetXaxis()->GetXmin(), hist_data_loose->GetXaxis()->GetXmax());
                fit_mc = new TF1("fit_mc", "gaus", hist_mc_loose->GetXaxis()->GetXmin(), hist_mc_loose->GetXaxis()->GetXmax());
            }

            // Fit data
            hist_data_loose->Fit(fit_data, "RQ0");
            fit_data->SetLineColor(kBlue);
            fit_data->SetLineStyle(2);  // Dashed line
            fit_data->SetLineWidth(2);  // Slightly thicker line
            fit_data->Draw("SAME");

            // Fit MC
            hist_mc_loose->Fit(fit_mc, "RQ0");
            fit_mc->SetLineColor(kRed);
            fit_mc->SetLineStyle(2);
            fit_mc->SetLineWidth(2);
            fit_mc->Draw("SAME");

            // Get mu and sigma from fits
            double mu_data, sigma_data, mu_mc, sigma_mc;

            if (std::find(variables_to_fit_asymmetric_gaussian.begin(), variables_to_fit_asymmetric_gaussian.end(), variables[i]) != variables_to_fit_asymmetric_gaussian.end()) {
                // For Asymmetric Gaussian
                mu_data = fit_data->GetParameter(1);
                double sigma_left_data = fit_data->GetParameter(2);
                double sigma_right_data = fit_data->GetParameter(3);

                mu_mc = fit_mc->GetParameter(1);
                double sigma_left_mc = fit_mc->GetParameter(2);
                double sigma_right_mc = fit_mc->GetParameter(3);

                // For simplicity, we can report the average sigma
                sigma_data = (sigma_left_data + sigma_right_data) / 2.0;
                sigma_mc = (sigma_left_mc + sigma_right_mc) / 2.0;
            } else {
                // Gaussian parameters: [0]=norm, [1]=mean, [2]=sigma
                mu_data = fit_data->GetParameter(1);
                sigma_data = fit_data->GetParameter(2);
                mu_mc = fit_mc->GetParameter(1);
                sigma_mc = fit_mc->GetParameter(2);
            }

            // Add a legend with mu and sigma from fits
            TLegend* legend_loose = new TLegend(0.275, 0.6, 0.9, 0.9);
            legend_loose->AddEntry(hist_data_loose, ("Data (" + std::to_string(static_cast<int>(hist_data_loose->GetEntries())) + " events; Loose Cuts)").c_str(), "lep");
            legend_loose->AddEntry(fit_data, Form("Data Fit: #mu=%.3f, #sigma=%.3f", mu_data, sigma_data), "l");
            legend_loose->AddEntry(hist_mc_loose, ("MC (" + std::to_string(static_cast<int>(hist_mc_loose->GetEntries())) + " events; Loose Cuts)").c_str(), "lep");
            legend_loose->AddEntry(fit_mc, Form("MC Fit: #mu=%.3f, #sigma=%.3f", mu_mc, sigma_mc), "l");
            legend_loose->SetTextSize(0.03);  // Set smaller font size for "Loose Cuts" legend
            legend_loose->Draw();

            // Note: Do not delete fit functions here to keep them on the canvas
        } else {
            // Add a legend without fits
            TLegend* legend_loose = new TLegend(0.275, 0.7, 0.9, 0.9);
            legend_loose->AddEntry(hist_data_loose, ("Data (" + std::to_string(static_cast<int>(hist_data_loose->GetEntries())) + " events; Loose Cuts)").c_str(), "lep");
            legend_loose->AddEntry(hist_mc_loose, ("MC (" + std::to_string(static_cast<int>(hist_mc_loose->GetEntries())) + " events; Loose Cuts)").c_str(), "lep");
            legend_loose->SetTextSize(0.03);  // Set smaller font size for "Loose Cuts" legend
            legend_loose->Draw();
        }
    }

    // Save the canvases with updated names to reflect the plotTitle and topology input
    std::string cleanTitle = plotTitle;
    std::replace(cleanTitle.begin(), cleanTitle.end(), ' ', '_');  // Replace spaces with underscores
    canvas->SaveAs((final_output_dir + "/exclusivity_plots_" + cleanTitle + "_" + topology + "_0_cuts.png").c_str());
    canvas_loose_cuts->SaveAs((final_output_dir + "/exclusivity_plots_" + cleanTitle + "_" + topology + "_1_cuts.png").c_str());

    // Clean up
    delete canvas; delete canvas_loose_cuts;
    delete theta_neutral_neutral_data; delete theta_neutral_neutral_mc;
}