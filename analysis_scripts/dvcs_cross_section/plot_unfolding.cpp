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

    // List of topologies
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)"};

    // Precompute the relevant bins for the xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);
    int n_Q2t_bins = relevant_bins.size();
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of rows and columns for the canvases
    int n_columns = 0;
    int n_rows = 0;

    if (xB_bin == 3 || xB_bin == 4) {
        n_columns = 5;
        n_rows = 6;
    } else {
        int n_subplots = next_perfect_square(n_Q2t_bins);
        n_columns = std::sqrt(n_subplots);
        n_rows = std::ceil(static_cast<double>(n_Q2t_bins) / n_columns);
    }

    // Create three sets of histograms, one for each topology
    std::vector<std::vector<TH1D*>> h_data_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));

    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Title with <xB>, <Q²>, and <-t> notation
            std::string title = Form("%s, %s: <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     analysisType.c_str(), topologies[topo_idx].c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            h_data_histograms[topo_idx][idx] = new TH1D(Form("h_data_%zu_%d", topo_idx, idx), title.c_str(), 24, 0, 360);

            // Set axis label and title font size
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetTitle("#phi");
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitle("Counts");
        }
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

    // Fill histograms based on topology in a single loop over data
    std::cout << "Started data " << std::endl;
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Apply bin cuts and kinematic cuts
            if ((*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high)) {

                // Fill histograms for the appropriate topology
                if (*detector1_data == 1 && *detector2_data == 1) {  // (FD,FD)
                    h_data_histograms[0][idx]->Fill(phi_deg);
                } else if (*detector1_data == 2 && *detector2_data == 1) {  // (CD,FD)
                    h_data_histograms[1][idx]->Fill(phi_deg);
                } else if (*detector1_data == 2 && *detector2_data == 0) {  // (CD,FT)
                    h_data_histograms[2][idx]->Fill(phi_deg);
                }
            }
        }
    }

    // Plot each topology on its own canvas
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        TCanvas* canvas = new TCanvas(Form("c_%zu", topo_idx), "Unfolded Distributions", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        gStyle->SetOptStat(0);

        // Plot histograms for this topology
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            TPad* pad = (TPad*)canvas->cd(idx + 1);  // Get the current pad (subplot)
        
            // Set margins
            pad->SetLeftMargin(0.15);
            pad->SetBottomMargin(0.15);

            TH1D* h_data = h_data_histograms[topo_idx][idx];

            // Find the maximum value for scaling
            double max_value = h_data->GetMaximum();
            h_data->SetMaximum(1.35 * max_value);

            // Draw the histogram
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

        // Save the canvas
        std::string channel_dir = (analysisType == "dvcs") ? "/dvcs" : "/eppi0";
        std::string filename = output_dir + "/unfolded" + channel_dir + "/unfolded_" + analysisType + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        canvas->SaveAs(filename.c_str());

        // Clean up the canvas
        delete canvas;
    }

    // Clean up histograms and memory
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_data_histograms[topo_idx]) {
            delete h;
        }
    }
    delete theta_neutral_neutral_data;
}