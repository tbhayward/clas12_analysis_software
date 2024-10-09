#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
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
                    TTreeReader& data_reader, 
                    TTreeReader& mc_gen_reader, 
                    TTreeReader& mc_rec_reader) {

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

    // Create histograms for data, rec MC, gen MC, and acceptance
    std::vector<std::vector<TH1D*>> h_data_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_mc_gen_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_mc_rec_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<std::vector<TH1D*>> h_acceptance_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));

    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            std::string title = Form("%s, %s: <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     analysisType.c_str(), topologies[topo_idx].c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            h_data_histograms[topo_idx][idx] = new TH1D(Form("h_data_%zu_%d", topo_idx, idx), title.c_str(), 24, 0, 360);
            h_mc_gen_histograms[topo_idx][idx] = new TH1D(Form("h_mc_gen_%zu_%d", topo_idx, idx), "MC Gen", 24, 0, 360);
            h_mc_rec_histograms[topo_idx][idx] = new TH1D(Form("h_mc_rec_%zu_%d", topo_idx, idx), "MC Rec", 24, 0, 360);
            h_acceptance_histograms[topo_idx][idx] = new TH1D(Form("h_acceptance_%zu_%d", topo_idx, idx), "Acceptance", 24, 0, 360);

            // Set axis label and title font size
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetTitle("#phi");
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
        }
    }

    // Restart the reader before looping over data
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Readers for necessary branches (assuming same variables for MC)
    TTreeReaderValue<double> phi_data(data_reader, "phi");
    TTreeReaderValue<double> xB_data(data_reader, "x");
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");

    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");

    // Fill histograms for data and MC
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high)) {

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

    // Fill histograms for MC generated and reconstructed
    while (mc_gen_reader.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];
            h_mc_gen_histograms[0][idx]->Fill(phi_mc_gen_deg);
        }
    }

    while (mc_rec_reader.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];
            h_mc_rec_histograms[0][idx]->Fill(phi_mc_rec_deg);
        }
    }

    // Compute and plot acceptance histograms
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            h_acceptance_histograms[topo_idx][idx]->Divide(h_mc_rec_histograms[topo_idx][idx], h_mc_gen_histograms[topo_idx][idx], 1, 1, "B");
            h_acceptance_histograms[topo_idx][idx]->GetYaxis()->SetTitle("Acceptance");
        }
    }

    // Plot and save histograms for yields and acceptances
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        TCanvas* canvas_yield = new TCanvas(Form("c_yield_%zu", topo_idx), "Unfolded Yields", canvas_width, canvas_height);
        TCanvas* canvas_acceptance = new TCanvas(Form("c_acceptance_%zu", topo_idx), "Acceptance", canvas_width, canvas_height);
        canvas_yield->Divide(n_columns, n_rows);
        canvas_acceptance->Divide(n_columns, n_rows);

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            // Yield plot
            TPad* pad_yield = (TPad*)canvas_yield->cd(idx + 1);
            pad_yield->SetLeftMargin(0.15);
            pad_yield->SetBottomMargin(0.15);
            h_data_histograms[topo_idx][idx]->Draw("E1");

            // Acceptance plot
            TPad* pad_acceptance = (TPad*)canvas_acceptance->cd(idx + 1);
            pad_acceptance->SetLeftMargin(0.15);
            pad_acceptance->SetBottomMargin(0.15);
            h_acceptance_histograms[topo_idx][idx]->Draw("E1");
        }

        // Save the canvases
        std::string channel_dir = (analysisType == "dvcs") ? "/dvcs" : "/eppi0";
        std::string filename_yield = output_dir + "/unfolded" + channel_dir + "/yields/yields_" + analysisType + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        std::string filename_acceptance = output_dir + "/unfolded" + channel_dir + "/acceptances/acceptances_" + analysisType + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";

        canvas_yield->SaveAs(filename_yield.c_str());
        canvas_acceptance->SaveAs(filename_acceptance.c_str());

        delete canvas_yield;
        delete canvas_acceptance;
    }

    // Clean up histograms
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_data_histograms[topo_idx]) delete h;
        for (auto& h : h_mc_gen_histograms[topo_idx]) delete h;
        for (auto& h : h_mc_rec_histograms[topo_idx]) delete h;
        for (auto& h : h_acceptance_histograms[topo_idx]) delete h;
    }
}