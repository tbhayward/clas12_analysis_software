// plot_dvcs_data_mc_comparison.cpp

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
#include <map>

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Helper function to find the closest lower power of 10 for a given value
double closest_lower_power_of_ten(double value) {
    if (value <= 0) return 0.1;  // Safety check for non-positive values
    return std::pow(10, std::floor(std::log10(value)));
}

void plot_dvcs_data_mc_comparison(const std::string& output_dir, 
                                  const std::string& analysisType, 
                                  const std::string& dataset, 
                                  int xB_bin,
                                  const std::vector<BinBoundary>& bin_boundaries, 
                                  TTreeReader& data_reader, 
                                  TTreeReader& mc_gen_reader, 
                                  TTreeReader& mc_rec_reader) {

    // Precompute the relevant bins for the xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);

    // Group the bins by (xB, Q2, t)
    std::map<std::tuple<double, double, double>, std::vector<int>> bin_groups;

    for (int idx : relevant_bins) {
        const auto& bin = bin_boundaries[idx];
        auto key = std::make_tuple(bin.xB_low, bin.Q2_low, bin.t_low);
        bin_groups[key].push_back(idx);
    }

    // Number of (xB, Q2, t) bins
    int n_bins = bin_groups.size();

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of columns and rows for the canvas
    int n_columns = std::sqrt(n_bins);
    int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

    TCanvas* canvas = new TCanvas("c1", "Data vs MC", canvas_width, canvas_height);
    canvas->Divide(n_columns, n_rows);

    gStyle->SetOptStat(0);

    // Histograms per (xB, Q2, t) bin
    std::vector<TH1D*> h_data_histograms;
    std::vector<TH1D*> h_mc_gen_histograms;
    std::vector<TH1D*> h_mc_rec_histograms;

    int pad_idx = 1;  // Pad index for canvas

    // Map to store the index mapping from (xB, Q2, t) keys to histogram indices
    std::map<std::tuple<double, double, double>, int> bin_index_map;

    // Create histograms for each (xB, Q2, t) bin with phi bins from bin boundaries
    for (const auto& group : bin_groups) {
        const auto& key = group.first;
        const auto& idx_list = group.second;

        // Collect phi bin edges for this (xB, Q2, t) bin
        std::vector<double> phi_bin_edges;

        for (int idx : idx_list) {
            const auto& bin = bin_boundaries[idx];
            phi_bin_edges.push_back(bin.phi_low);
        }
        // Add the last upper edge
        phi_bin_edges.push_back(bin_boundaries[idx_list.back()].phi_high);

        // Remove duplicate phi edges and sort them
        std::sort(phi_bin_edges.begin(), phi_bin_edges.end());
        phi_bin_edges.erase(std::unique(phi_bin_edges.begin(), phi_bin_edges.end()), phi_bin_edges.end());

        // Create histograms
        std::string title = Form("%s, %s: x_{B} [%.2f, %.2f], Q^{2} [%.2f, %.2f], -t [%.2f, %.2f]", 
                                 analysisType.c_str(), 
                                 dataset.c_str(),  
                                 std::get<0>(key), std::get<0>(key) + bin_boundaries[idx_list[0]].xB_high - bin_boundaries[idx_list[0]].xB_low,
                                 std::get<1>(key), std::get<1>(key) + bin_boundaries[idx_list[0]].Q2_high - bin_boundaries[idx_list[0]].Q2_low,
                                 std::get<2>(key), std::get<2>(key) + bin_boundaries[idx_list[0]].t_high - bin_boundaries[idx_list[0]].t_low);

        // Convert phi edges to degrees
        std::vector<double> phi_edges_deg;
        for (double edge : phi_bin_edges) {
            phi_edges_deg.push_back(edge);
        }

        // Create histograms
        TH1D* h_data = new TH1D(Form("h_data_%d", pad_idx - 1), title.c_str(), phi_edges_deg.size() - 1, &phi_edges_deg[0]);
        TH1D* h_mc_gen = new TH1D(Form("h_mc_gen_%d", pad_idx - 1), "gen", phi_edges_deg.size() - 1, &phi_edges_deg[0]);
        TH1D* h_mc_rec = new TH1D(Form("h_mc_rec_%d", pad_idx - 1), "rec", phi_edges_deg.size() - 1, &phi_edges_deg[0]);

        // Axis labels and styles
        h_data->GetXaxis()->SetTitle("#phi [deg]");
        h_data->GetYaxis()->SetTitle("Counts");
        h_data->GetXaxis()->SetLabelSize(0.05);
        h_data->GetYaxis()->SetLabelSize(0.05);
        h_data->GetXaxis()->SetTitleSize(0.06);
        h_data->GetYaxis()->SetTitleSize(0.06);

        h_mc_gen->GetXaxis()->SetLabelSize(0.05);
        h_mc_gen->GetYaxis()->SetLabelSize(0.05);
        h_mc_gen->GetXaxis()->SetTitleSize(0.06);
        h_mc_gen->GetYaxis()->SetTitleSize(0.06);

        h_mc_rec->GetXaxis()->SetLabelSize(0.05);
        h_mc_rec->GetYaxis()->SetLabelSize(0.05);
        h_mc_rec->GetXaxis()->SetTitleSize(0.06);
        h_mc_rec->GetYaxis()->SetTitleSize(0.06);

        h_data_histograms.push_back(h_data);
        h_mc_gen_histograms.push_back(h_mc_gen);
        h_mc_rec_histograms.push_back(h_mc_rec);

        bin_index_map[key] = pad_idx - 1;  // Map the key to histogram index

        ++pad_idx;
    }

    // Readers for data and MC
    TTreeReaderValue<double> phi_data(data_reader, "phi2");
    TTreeReaderValue<double> xB_data(data_reader, "x");
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");

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

    // Restart the readers before looping over data
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Fill data histograms
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;  
        double xB = *xB_data;
        double Q2 = *Q2_data;
        double t_abs = std::abs(*t1_data);

        // Find the (xB, Q2, t) bin
        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            const auto& bin_example = bin_boundaries[idx_list[0]];  // Use first bin as representative

            if (xB >= bin_example.xB_low && xB <= bin_example.xB_high &&
                Q2 >= bin_example.Q2_low && Q2 <= bin_example.Q2_high &&
                t_abs >= bin_example.t_low && t_abs <= bin_example.t_high &&
                apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, **theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

                int hist_idx = bin_index_map[key];
                TH1D* h_data = h_data_histograms[hist_idx];

                // Find the phi bin within this (xB, Q2, t) bin
                for (int idx : idx_list) {
                    const auto& bin = bin_boundaries[idx];
                    if (phi_deg >= bin.phi_low && phi_deg < bin.phi_high) {
                        h_data->Fill(phi_deg);
                        break;
                    }
                }
                break;  // Exit the loop once the event is assigned
            }
        }
    }

    // Fill mc_gen histograms
    while (mc_gen_reader.Next()) {
        double phi_deg = *phi_mc_gen * RAD_TO_DEG; 
        double xB = *xB_mc_gen;
        double Q2 = *Q2_mc_gen;
        double t_abs = std::abs(*t1_mc_gen);

        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            const auto& bin_example = bin_boundaries[idx_list[0]];

            if (xB >= bin_example.xB_low && xB <= bin_example.xB_high &&
                Q2 >= bin_example.Q2_low && Q2 <= bin_example.Q2_high &&
                t_abs >= bin_example.t_low && t_abs <= bin_example.t_high &&
                apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, **theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen)) {

                int hist_idx = bin_index_map[key];
                TH1D* h_mc_gen = h_mc_gen_histograms[hist_idx];

                // Find the phi bin within this (xB, Q2, t) bin
                for (int idx : idx_list) {
                    const auto& bin = bin_boundaries[idx];
                    if (phi_deg >= bin.phi_low && phi_deg < bin.phi_high) {
                        h_mc_gen->Fill(phi_deg);
                        break;
                    }
                }
                break;
            }
        }
    }

    // Fill mc_rec histograms
    while (mc_rec_reader.Next()) {
        double phi_deg = *phi_mc_rec * RAD_TO_DEG; 
        double xB = *xB_mc_rec;
        double Q2 = *Q2_mc_rec;
        double t_abs = std::abs(*t1_mc_rec);

        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            const auto& bin_example = bin_boundaries[idx_list[0]];

            if (xB >= bin_example.xB_low && xB <= bin_example.xB_high &&
                Q2 >= bin_example.Q2_low && Q2 <= bin_example.Q2_high &&
                t_abs >= bin_example.t_low && t_abs <= bin_example.t_high &&
                apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, **theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec)) {

                int hist_idx = bin_index_map[key];
                TH1D* h_mc_rec = h_mc_rec_histograms[hist_idx];

                // Find the phi bin within this (xB, Q2, t) bin
                for (int idx : idx_list) {
                    const auto& bin = bin_boundaries[idx];
                    if (phi_deg >= bin.phi_low && phi_deg < bin.phi_high) {
                        h_mc_rec->Fill(phi_deg);
                        break;
                    }
                }
                break;
            }
        }
    }

    // Normalize histograms and plot
    pad_idx = 1;
    for (const auto& group : bin_groups) {
        const auto& key = group.first;

        int hist_idx = bin_index_map[key];

        canvas->cd(pad_idx);

        TPad* pad = (TPad*)gPad;

        // Adjust margins if necessary
        pad->SetLeftMargin(0.15);
        pad->SetBottomMargin(0.15);

        TH1D* h_data = h_data_histograms[hist_idx];
        TH1D* h_mc_gen = h_mc_gen_histograms[hist_idx];
        TH1D* h_mc_rec = h_mc_rec_histograms[hist_idx];

        // Scaling and plotting as before
        double integral_data = h_data->Integral();
        double integral_mc_rec = h_mc_rec->Integral();
        double integral_mc_gen = h_mc_gen->Integral();

        // Scale reconstructed MC to have the same integral as data
        if (integral_mc_rec > 0 && integral_data > 0) {
            h_mc_rec->Scale(integral_data / integral_mc_rec);
        }

        // Now scale generated MC using the ratio of its original integral to the original reconstructed MC integral
        if (integral_mc_gen > 0 && integral_mc_rec > 0) {
            h_mc_gen->Scale((integral_mc_gen / integral_mc_rec) * (integral_data / integral_mc_rec));
        }

        // Set styles
        h_data->SetMarkerColor(kBlue);
        h_data->SetMarkerStyle(20);
        h_data->SetLineColor(kBlue);

        h_mc_gen->SetLineColor(kBlack);
        h_mc_gen->SetLineStyle(2);

        h_mc_rec->SetMarkerColor(kRed);
        h_mc_rec->SetMarkerStyle(22);
        h_mc_rec->SetLineColor(kRed);

        // Draw histograms
        h_data->Draw("E1");
        h_mc_rec->Draw("E1 SAME");
        h_mc_gen->Draw("HIST SAME");

        // Add legend
        TLegend* legend = new TLegend(0.575, 0.45, 0.9, 0.75);
        legend->AddEntry(h_data, "Data", "lep");
        legend->AddEntry(h_mc_rec, "Reconstructed MC", "lep");
        legend->AddEntry(h_mc_gen, "Generated MC", "l");
        legend->SetTextSize(0.04);
        legend->Draw();

        ++pad_idx;
    }

    // Save canvas
    std::string filename = output_dir + "/phi_data_mc_comparison_" + analysisType + "_" + dataset + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    canvas->SaveAs(filename.c_str());

    // Clean up
    for (auto& h : h_data_histograms) delete h;
    for (auto& h : h_mc_gen_histograms) delete h;
    for (auto& h : h_mc_rec_histograms) delete h;
    delete canvas;

    delete theta_neutral_neutral_data;
    delete theta_neutral_neutral_mc_gen;
    delete theta_neutral_neutral_mc_rec;
}