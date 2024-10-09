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

    // Precompute the relevant bins for the xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);

    int n_Q2t_bins = relevant_bins.size();
    std::cout << "Current xB_bin = " << xB_bin << ", Number of Q2t bins: " << n_Q2t_bins << std::endl;

    // Create a vector of histograms for each topology
    std::vector<std::vector<TH1D*>> h_data_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));

    // Create histograms for each topology and relevant bin
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        const std::string& topology = topologies[topo_idx];
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Update title to include the analysisType, topology, and averages
            std::string title = Form("%s, %s: x_{B} avg: %.2f, Q^{2} avg: %.2f, -t avg: %.2f", 
                                     analysisType.c_str(), topology.c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            h_data_histograms[topo_idx][idx] = new TH1D(Form("h_data_%zu_%d", topo_idx, idx), title.c_str(), 24, 0, 360);

            // Increase axis label and title font size
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

    // Data histograms filling with topology check
    std::cout << "Started data " << std::endl;
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Apply bin cuts and kinematic cuts
            if ((*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high)) {

                // Loop over topologies and fill corresponding histograms
                for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
                    const std::string& topology = topologies[topo_idx];

                    if (((topology == "(FD,FD)" && *detector1_data == 1 && *detector2_data == 1) ||
                         (topology == "(CD,FD)" && *detector1_data == 2 && *detector2_data == 1) ||
                         (topology == "(CD,FT)" && *detector1_data == 2 && *detector2_data == 0) ) &&
                        apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, **theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {
                        
                        // Fill the corresponding topology histogram
                        h_data_histograms[topo_idx][idx]->Fill(phi_deg);
                    }
                }
            }
        }
    }

    // Create canvases for each topology
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        // Create a new canvas for the current topology
        TCanvas* canvas = new TCanvas(Form("c_%zu", topo_idx), "Unfolded Distributions", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        // Normalize histograms and plot in each subplot
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            TPad* pad = (TPad*)canvas->cd(idx + 1);  // Get the current pad (subplot)
        
            // Set the margins for the subplot
            pad->SetLeftMargin(0.15);   // Add space to the left
            pad->SetBottomMargin(0.15); // Add space below

            TH1D* h_data = h_data_histograms[topo_idx][idx];

            // Find the maximum bin content for plotting
            double max_value = h_data->GetMaximum();

            // Set the maximum for the histograms
            h_data->SetMaximum(1.35 * max_value);

            // Draw histograms
            h_data->SetMarkerColor(kBlue);  // Color for this topology
            h_data->SetMarkerStyle(20);
            h_data->SetLineColor(kBlue);
            h_data->Draw("E1");

            // Add legend
            TLegend* legend = new TLegend(0.575, 0.45, 0.9, 0.75);
            legend->AddEntry(h_data, Form("Data %s", topologies[topo_idx].c_str()), "lep");
            legend->SetTextSize(0.04);
            legend->Draw();
        }

        // Save canvas to the output directory
        std::string channel_dir = (analysisType == "dvcs") ? "/dvcs" : "/eppi0";
        std::string filename = output_dir + "/unfolded" + channel_dir + "/unfolded_" + analysisType + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        canvas->SaveAs(filename.c_str());

        // Clean up the canvas
        delete canvas;
    }

    // Clean up histograms
    for (auto& topo_histograms : h_data_histograms) {
        for (auto& h : topo_histograms) {
            delete h;
        }
    }
    delete theta_neutral_neutral_data;
}