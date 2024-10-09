#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>
#include <cmath>
#include <string>
#include <vector>
#include "bin_boundaries.h"
#include <algorithm>
#include <cctype>
#include "kinematic_cuts.h"
#include "bin_helpers.h"

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

void plot_unfolding(const std::string& output_dir, 
                    const std::string& analysisType, 
                    int xB_bin,
                    const std::vector<BinBoundary>& bin_boundaries, 
                    TTreeReader& data_reader) {

    // List of topologies to loop over
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)"};

    for (const auto& topology : topologies) {
        // Precompute the relevant bins for the xB_bin
        std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);

        int n_Q2t_bins = relevant_bins.size();
        std::cout << "Current xB_bin = " << xB_bin << ", Topology = " << topology << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

        // Adjust canvas size
        const int base_canvas_width = 1200;
        const int base_canvas_height = 800;
        int canvas_width = static_cast<int>(1.5 * base_canvas_width);
        int canvas_height = static_cast<int>(1.5 * base_canvas_height);

        // For certain xB bins (3 and 4), we use a fixed 6x5 layout, otherwise use next perfect square logic
        int n_subplots;
        int n_columns;
        int n_rows;

        if (xB_bin == 3 || xB_bin == 4) {
            n_subplots = 30;
            n_columns = 5;
            n_rows = 6;
        } else {
            n_subplots = next_perfect_square(n_Q2t_bins);
            n_columns = std::sqrt(n_subplots);
            n_rows = std::ceil(static_cast<double>(n_Q2t_bins) / n_columns);
        }

        TCanvas* canvas = new TCanvas("c1", "Unfolded Distributions", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        gStyle->SetOptStat(0);

        std::vector<TH1D*> h_data_histograms(n_Q2t_bins);

        // Create histograms only for the relevant bins
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Update title to include the analysisType, topology, and averages
            std::string title = Form("%s, %s: x_{B} avg: %.2f, Q^{2} avg: %.2f, -t avg: %.2f", 
                                     analysisType.c_str(), topology.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            h_data_histograms[idx] = new TH1D(Form("h_data_%d", idx), title.c_str(), 24, 0, 360);

            // Increase axis label and title font size
            h_data_histograms[idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms[idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms[idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms[idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms[idx]->GetXaxis()->SetTitle("#phi");
            h_data_histograms[idx]->GetYaxis()->SetTitle("Counts");
        }

        // Restart the reader before looping over data
        data_reader.Restart();

        // Readers for necessary branches
        TTreeReaderValue<double> phi_data(data_reader, "phi");
        TTreeReaderValue<double> xB_data(data_reader, "x");
        TTreeReaderValue<double> Q2_data(data_reader, "Q2");
        TTreeReaderValue<double> t1_data(data_reader, "t1");
        TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
        TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
        TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
        TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");

        // Import detector information for topology check
        TTreeReaderValue<int> detector1_data(data_reader, "detector1");
        TTreeReaderValue<int> detector2_data(data_reader, "detector2");

        // Handle theta_neutral_neutral based on analysis type (dvcs or eppi0)
        TTreeReaderValue<double>* theta_neutral_neutral_data;
        if (analysisType == "dvcs") {
            theta_neutral_neutral_data = new TTreeReaderValue<double>(data_reader, "theta_gamma_gamma");
        } else if (analysisType == "eppi0") {
            theta_neutral_neutral_data = new TTreeReaderValue<double>(data_reader, "theta_pi0_pi0");
        }

        // Data histograms filling with topology check
        std::cout << "Started data " << std::endl;
        while (data_reader.Next()) {
            double phi_deg = *phi_data * RAD_TO_DEG;

            for (int idx = 0; idx < n_Q2t_bins; ++idx) {
                const auto& bin = bin_boundaries[relevant_bins[idx]];

                // Apply bin cuts, kinematic cuts, and topology conditions
                if ((*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                    *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                    std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high) &&
                    ((topology == "(FD,FD)" && *detector1_data == 1 && *detector2_data == 1) ||
                     (topology == "(CD,FD)" && *detector1_data == 2 && *detector2_data == 1) ||
                     (topology == "(CD,FT)" && *detector1_data == 2 && *detector2_data == 0)) &&
                    apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, **theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

                    // Fill the data histogram
                    h_data_histograms[idx]->Fill(phi_deg);
                    break;
                }
            }
        }

        // Normalize histograms and plot in each subplot
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            TPad* pad = (TPad*)canvas->cd(idx + 1);  // Get the current pad (subplot)
        
            // Set the margins for the subplot
            pad->SetLeftMargin(0.15);   // Add space to the left
            pad->SetBottomMargin(0.15); // Add space below

            TH1D* h_data = h_data_histograms[idx];

            // Find the maximum bin content for plotting
            double max_value = h_data->GetMaximum();

            // Set the maximum for the histograms
            h_data->SetMaximum(1.35 * max_value);

            // Draw histograms
            h_data->SetMarkerColor(kBlue);
            h_data->SetMarkerStyle(20);
            h_data->SetLineColor(kBlue);
            h_data->Draw("E1");

            // Add legend
            TLegend* legend = new TLegend(0.575, 0.45, 0.9, 0.75);
            legend->AddEntry(h_data, "Data", "lep");
            legend->SetTextSize(0.04);
            legend->Draw();
        }

        // Save canvas to the output directory
        std::string channel_dir = (analysisType == "dvcs") ? "/dvcs" : "/eppi0";
        std::string filename = output_dir + "/unfolded" + channel_dir + "/unfolded_" + analysisType + "_" + topology + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        canvas->SaveAs(filename.c_str());

        // Clean up histograms and canvas
        for (auto& h : h_data_histograms) delete h;
        delete canvas;
        delete theta_neutral_neutral_data;
    }
}