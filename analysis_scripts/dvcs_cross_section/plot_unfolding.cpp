#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>
#include <string>
#include <vector>
#include <filesystem>
#include "bin_boundaries.h"
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

    // Set style to remove stat boxes
    gStyle->SetOptStat(0);

    // Ensure the necessary directories are created
    create_directories(output_dir);

    // List of topologies and combined option
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)", "combined"};
    std::string channel_dir = (analysisType == "dvcs") ? "/dvcs" : "/eppi0";

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
    int n_columns = (xB_bin == 3 || xB_bin == 4) ? 5 : std::sqrt(next_perfect_square(n_Q2t_bins));
    int n_rows = std::ceil(static_cast<double>(n_Q2t_bins) / n_columns);

    // Create histograms for data, rec MC, gen MC, and acceptance (for combined only)
    std::vector<std::vector<TH1D*>> h_data_histograms(topologies.size(), std::vector<TH1D*>(n_Q2t_bins));
    std::vector<TH1D*> h_mc_gen_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_mc_rec_histograms(n_Q2t_bins);
    std::vector<TH1D*> h_acceptance_histograms(n_Q2t_bins);

    // Initialize histograms
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            // Create title string with the channel and bin information
            std::string title = Form("%s, %s: <x_{B}>: %.2f, <Q^{2}>: %.2f, <-t>: %.2f", 
                                     analysisType.c_str(), topologies[topo_idx].c_str(),
                                     bin.xB_avg, bin.Q2_avg, std::abs(bin.t_avg));

            // Create histogram for the data yields with the appropriate title
            h_data_histograms[topo_idx][idx] = new TH1D(Form("h_data_%zu_%d", topo_idx, idx), title.c_str(), 24, 0, 360);

            // Set axis labels and format for data histograms
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetLabelSize(0.05);
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetLabelSize(0.05);
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetTitleSize(0.06);
            h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitleSize(0.06);
            h_data_histograms[topo_idx][idx]->GetXaxis()->SetTitle("#phi");

            // Set Y-axis title based on whether it's a "Raw Yield" or "Unfolded Yield"
            if (topo_idx == 3) {
                h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitle("Unfolded Yield");
            } else {
                h_data_histograms[topo_idx][idx]->GetYaxis()->SetTitle("Raw Yield");
            }

            // Set markers and draw vertical error bars without horizontal error bars
            h_data_histograms[topo_idx][idx]->SetMarkerStyle(20);
            h_data_histograms[topo_idx][idx]->SetMarkerSize(1.2);
            h_data_histograms[topo_idx][idx]->SetDrawOption("P E1");  // 'P' for points, 'E1' for vertical error bars only

            // Create histograms for MC and acceptance if it's the combined histogram
            if (topo_idx == 3) {
                h_mc_gen_histograms[idx] = new TH1D(Form("h_mc_gen_combined_%d", idx), title.c_str(), 24, 0, 360);
                h_mc_rec_histograms[idx] = new TH1D(Form("h_mc_rec_combined_%d", idx), title.c_str(), 24, 0, 360);
                h_acceptance_histograms[idx] = new TH1D(Form("h_acceptance_combined_%d", idx), title.c_str(), 24, 0, 360);

                // Set axis labels and format for acceptance histograms
                h_acceptance_histograms[idx]->GetXaxis()->SetLabelSize(0.05);
                h_acceptance_histograms[idx]->GetYaxis()->SetLabelSize(0.05);
                h_acceptance_histograms[idx]->GetXaxis()->SetTitleSize(0.06);
                h_acceptance_histograms[idx]->GetYaxis()->SetTitleSize(0.06);
                h_acceptance_histograms[idx]->GetXaxis()->SetTitle("#phi");
                h_acceptance_histograms[idx]->GetYaxis()->SetTitle("Acceptance");

                // Set markers and draw vertical error bars without horizontal error bars for acceptance histograms
                h_acceptance_histograms[idx]->SetMarkerStyle(20);
                h_acceptance_histograms[idx]->SetMarkerSize(1.2);
                h_acceptance_histograms[idx]->SetDrawOption("P E1");
            }
        }
    }

    // Readers for necessary branches in all datasets (data, mc_gen, mc_rec)
    TTreeReaderValue<int> detector1_data(data_reader, "detector1");
    TTreeReaderValue<int> detector2_data(data_reader, "detector2");
    TTreeReaderValue<double> phi_data(data_reader, "phi");
    TTreeReaderValue<double> xB_data(data_reader, "x");
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");

    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
    TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
    TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
    TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");

    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");
    TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
    TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
    TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_rec(mc_rec_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec(mc_rec_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec(mc_rec_reader, "pTmiss");

    // Handle theta_neutral_neutral based on analysis type (dvcs or eppi0)
    TTreeReaderValue<double>* theta_neutral_neutral_data;
    TTreeReaderValue<double>* theta_neutral_neutral_mc_gen;
    TTreeReaderValue<double>* theta_neutral_neutral_mc_rec;
    if (analysisType == "dvcs") {
        theta_neutral_neutral_data = new TTreeReaderValue<double>(data_reader, "theta_gamma_gamma");
        theta_neutral_neutral_mc_gen = new TTreeReaderValue<double>(mc_gen_reader, "theta_gamma_gamma");
        theta_neutral_neutral_mc_rec = new TTreeReaderValue<double>(mc_rec_reader, "theta_gamma_gamma");
    } else if (analysisType == "eppi0") {
        theta_neutral_neutral_data = new TTreeReaderValue<double>(data_reader, "theta_pi0_pi0");
        theta_neutral_neutral_mc_gen = new TTreeReaderValue<double>(mc_gen_reader, "theta_pi0_pi0");
        theta_neutral_neutral_mc_rec = new TTreeReaderValue<double>(mc_rec_reader, "theta_pi0_pi0");
    }

    // Print before starting loops
    std::cout << "starting data" << std::endl;
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, **theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

                // Fill data histograms based on detector topologies
                if (*detector1_data == 1 && *detector2_data == 1) {  // (FD,FD)
                    h_data_histograms[0][idx]->Fill(phi_deg);
                } else if (*detector1_data == 2 && *detector2_data == 1) {  // (CD,FD)
                    h_data_histograms[1][idx]->Fill(phi_deg);
                } else if (*detector1_data == 2 && *detector2_data == 0) {  // (CD,FT)
                    h_data_histograms[2][idx]->Fill(phi_deg);
                }

                // Also fill the combined histogram
                h_data_histograms[3][idx]->Fill(phi_deg);
            }
        }
    }

    std::cout << "starting gen mc" << std::endl;
    while (mc_gen_reader.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_gen >= bin.xB_low && *xB_mc_gen <= bin.xB_high &&
                *Q2_mc_gen >= bin.Q2_low && *Q2_mc_gen <= bin.Q2_high &&
                std::abs(*t1_mc_gen) >= bin.t_low && std::abs(*t1_mc_gen) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, **theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen)) {

                // Fill generated MC histogram (combined only)
                h_mc_gen_histograms[idx]->Fill(phi_mc_gen_deg);
            }
        }
    }

    std::cout << "starting rec mc" << std::endl;
    while (mc_rec_reader.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            const auto& bin = bin_boundaries[relevant_bins[idx]];

            if ((*xB_mc_rec >= bin.xB_low && *xB_mc_rec <= bin.xB_high &&
                *Q2_mc_rec >= bin.Q2_low && *Q2_mc_rec <= bin.Q2_high &&
                std::abs(*t1_mc_rec) >= bin.t_low && std::abs(*t1_mc_rec) <= bin.t_high) &&
                apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, **theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec)) {

                // Fill reconstructed MC histogram (combined only)
                h_mc_rec_histograms[idx]->Fill(phi_mc_rec_deg);
            }
        }
    }

    // Compute acceptance for combined
    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        if (h_mc_gen_histograms[idx]->Integral() > 0) {
            h_acceptance_histograms[idx]->Divide(h_mc_rec_histograms[idx], h_mc_gen_histograms[idx], 1, 1, "B");
        }
    }

    // Plot and save histograms (data divided by acceptance for combined)
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        TCanvas* canvas_yield = new TCanvas(Form("c_yield_%zu", topo_idx), "Yields", canvas_width, canvas_height);
        canvas_yield->Divide(n_columns, n_rows);

        for (int idx = 0; idx < n_Q2t_bins; ++idx) {
            TPad* pad_yield = (TPad*)canvas_yield->cd(idx + 1);
            pad_yield->SetLeftMargin(0.15);
            pad_yield->SetBottomMargin(0.15);

            // If this is the combined topology, we need to unfold by dividing data by acceptance
            if (topo_idx == 3) {  // Combined histograms
                for (int bin = 1; bin <= h_data_histograms[topo_idx][idx]->GetNbinsX(); ++bin) {
                    double acceptance_value = h_acceptance_histograms[idx]->GetBinContent(bin);
                    if (acceptance_value > 0) {
                        double data_value = h_data_histograms[topo_idx][idx]->GetBinContent(bin);
                        double unfolded_value = data_value / acceptance_value;
                        h_data_histograms[topo_idx][idx]->SetBinContent(bin, unfolded_value);

                        // Calculate and set appropriate error for unfolded values
                        double data_error = h_data_histograms[topo_idx][idx]->GetBinError(bin);
                        double unfolded_error = data_error / acceptance_value;
                        h_data_histograms[topo_idx][idx]->SetBinError(bin, unfolded_error);
                    } else {
                        // Set the unfolded value to 0 if the acceptance is 0
                        h_data_histograms[topo_idx][idx]->SetBinContent(bin, 0);
                        h_data_histograms[topo_idx][idx]->SetBinError(bin, 0);
                    }
                }
            }

            // Draw the histogram
            h_data_histograms[topo_idx][idx]->Draw("E1");
        }

        std::string filename_yield = output_dir + "/unfolded" + channel_dir + "/yields/yields_" + analysisType + "_" + topologies[topo_idx] + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
        canvas_yield->SaveAs(filename_yield.c_str());

        delete canvas_yield;
    }

    // Plot and save acceptance for combined only
    TCanvas* canvas_acceptance = new TCanvas("c_acceptance_combined", "Acceptance Combined", canvas_width, canvas_height);
    canvas_acceptance->Divide(n_columns, n_rows);

    for (int idx = 0; idx < n_Q2t_bins; ++idx) {
        TPad* pad_acceptance = (TPad*)canvas_acceptance->cd(idx + 1);
        pad_acceptance->SetLeftMargin(0.20);
        pad_acceptance->SetBottomMargin(0.15);
        h_acceptance_histograms[idx]->Draw("E1");
    }

    std::string filename_acceptance = output_dir + "/unfolded" + channel_dir + "/acceptances/acceptances_combined_" + analysisType + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    canvas_acceptance->SaveAs(filename_acceptance.c_str());

    delete canvas_acceptance;

    // Clean up histograms
    for (size_t topo_idx = 0; topo_idx < topologies.size(); ++topo_idx) {
        for (auto& h : h_data_histograms[topo_idx]) delete h;
    }

    for (auto& h : h_mc_gen_histograms) delete h;
    for (auto& h : h_mc_rec_histograms) delete h;
    for (auto& h : h_acceptance_histograms) delete h;

    delete theta_neutral_neutral_data;
    delete theta_neutral_neutral_mc_gen;
    delete theta_neutral_neutral_mc_rec;

    // Reset the readers after each iteration
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();
}