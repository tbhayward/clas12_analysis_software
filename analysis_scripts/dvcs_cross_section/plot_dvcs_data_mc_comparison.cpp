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
#include <set>

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

    // Collect all unique phi bin edges from bin boundaries
    std::set<double> global_phi_edges_set;

    for (const auto& idx : relevant_bins) {
        const auto& bin = bin_boundaries[idx];
        global_phi_edges_set.insert(bin.phi_low);
        global_phi_edges_set.insert(bin.phi_high);
    }

    // Convert set to sorted vector
    std::vector<double> global_phi_edges(global_phi_edges_set.begin(), global_phi_edges_set.end());
    std::sort(global_phi_edges.begin(), global_phi_edges.end());

    // Ensure we have at least two edges to create histograms
    if (global_phi_edges.size() < 2) {
        std::cerr << "Error: Insufficient global phi bin edges. Cannot create histograms." << std::endl;
        return;
    }

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of columns and rows for the canvas
    int n_bins = bin_groups.size();
    int n_columns = std::ceil(std::sqrt(n_bins));
    int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

    TCanvas* canvas = new TCanvas("c1", "Data vs MC", canvas_width, canvas_height);
    canvas->Divide(n_columns, n_rows);

    gStyle->SetOptStat(0);

    // Histograms per (xB, Q2, t) bin
    std::map<std::tuple<double, double, double>, TH1D*> h_data_histograms;
    std::map<std::tuple<double, double, double>, TH1D*> h_mc_gen_histograms;
    std::map<std::tuple<double, double, double>, TH1D*> h_mc_rec_histograms;

    int pad_idx = 1;  // Pad index for canvas

    // Restart the readers before declaring TTreeReaderValue objects
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Handle theta_neutral_neutral based on analysis type (dvcs or eppi0)
    std::string theta_variable_name;
    if (analysisType == "dvcs") {
        theta_variable_name = "theta_gamma_gamma";
    } else if (analysisType == "eppi0") {
        theta_variable_name = "theta_pi0_pi0";
    } else {
        std::cerr << "Error: Unknown analysisType '" << analysisType << "'" << std::endl;
        return;
    }

    // Readers for data
    TTreeReaderValue<double> phi_data(data_reader, "phi2");
    TTreeReaderValue<double> xB_data(data_reader, "x");
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");
    TTreeReaderValue<double> theta_neutral_neutral_data(data_reader, theta_variable_name.c_str());

    // Readers for MC-generated data
    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi2");
    TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
    TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
    TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen(mc_gen_reader, theta_variable_name.c_str());

    // Readers for MC-reconstructed data
    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi2");
    TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
    TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
    TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_rec(mc_rec_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec(mc_rec_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec(mc_rec_reader, "pTmiss");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec(mc_rec_reader, theta_variable_name.c_str());

    // Create histograms for each (xB, Q2, t) bin using global phi bin edges
    for (const auto& group : bin_groups) {
        const auto& key = group.first;

        // Create histograms with global phi bin edges
        std::string title = Form("%s, %s: x_{B} [%.2f, %.2f], Q^{2} [%.2f, %.2f], -t [%.2f, %.2f]", 
                                 analysisType.c_str(), 
                                 dataset.c_str(),  
                                 std::get<0>(key), bin_boundaries[group.second[0]].xB_high,
                                 std::get<1>(key), bin_boundaries[group.second[0]].Q2_high,
                                 std::get<2>(key), bin_boundaries[group.second[0]].t_high);

        TH1D* h_data = new TH1D(Form("h_data_%d", pad_idx - 1), title.c_str(), global_phi_edges.size() - 1, &global_phi_edges[0]);
        TH1D* h_mc_gen = new TH1D(Form("h_mc_gen_%d", pad_idx - 1), "gen", global_phi_edges.size() - 1, &global_phi_edges[0]);
        TH1D* h_mc_rec = new TH1D(Form("h_mc_rec_%d", pad_idx - 1), "rec", global_phi_edges.size() - 1, &global_phi_edges[0]);

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

        h_data_histograms[key] = h_data;
        h_mc_gen_histograms[key] = h_mc_gen;
        h_mc_rec_histograms[key] = h_mc_rec;

        ++pad_idx;
    }

    // Fill data histograms
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;  
        double xB = *xB_data;
        double Q2 = *Q2_data;
        double t_abs = std::abs(*t1_data);

        // Adjust phi_deg to be within [0, 360)
        phi_deg = std::fmod(phi_deg + 360.0, 360.0);

        // Find the (xB, Q2, t) bin
        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            const auto& bin_example = bin_boundaries[idx_list[0]];  // Use first bin as representative

            if (xB >= bin_example.xB_low && xB <= bin_example.xB_high &&
                Q2 >= bin_example.Q2_low && Q2 <= bin_example.Q2_high &&
                t_abs >= bin_example.t_low && t_abs <= bin_example.t_high &&
                apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {

                TH1D* h_data = h_data_histograms[key];
                if (h_data) {
                    h_data->Fill(phi_deg);
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

        // Adjust phi_deg to be within [0, 360)
        phi_deg = std::fmod(phi_deg + 360.0, 360.0);

        for (const auto& group : bin_groups) {
            const auto& key = group.first;

            const auto& bin_example = bin_boundaries[group.second[0]];

            if (xB >= bin_example.xB_low && xB <= bin_example.xB_high &&
                Q2 >= bin_example.Q2_low && Q2 <= bin_example.Q2_high &&
                t_abs >= bin_example.t_low && t_abs <= bin_example.t_high &&
                apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, *theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen)) {

                TH1D* h_mc_gen = h_mc_gen_histograms[key];
                if (h_mc_gen) {
                    h_mc_gen->Fill(phi_deg);
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

        // Adjust phi_deg to be within [0, 360)
        phi_deg = std::fmod(phi_deg + 360.0, 360.0);

        for (const auto& group : bin_groups) {
            const auto& key = group.first;

            const auto& bin_example = bin_boundaries[group.second[0]];

            if (xB >= bin_example.xB_low && xB <= bin_example.xB_high &&
                Q2 >= bin_example.Q2_low && Q2 <= bin_example.Q2_high &&
                t_abs >= bin_example.t_low && t_abs <= bin_example.t_high &&
                apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, *theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec)) {

                TH1D* h_mc_rec = h_mc_rec_histograms[key];
                if (h_mc_rec) {
                    h_mc_rec->Fill(phi_deg);
                }
                break;
            }
        }
    }

    // Normalize histograms and plot
    pad_idx = 1;
    for (const auto& group : bin_groups) {
        const auto& key = group.first;

        canvas->cd(pad_idx);

        TPad* pad = (TPad*)gPad;

        // Adjust margins if necessary
        pad->SetLeftMargin(0.15);
        pad->SetBottomMargin(0.15);

        TH1D* h_data = h_data_histograms[key];
        TH1D* h_mc_gen = h_mc_gen_histograms[key];
        TH1D* h_mc_rec = h_mc_rec_histograms[key];

        if (!h_data || !h_mc_gen || !h_mc_rec) {
            std::cerr << "Error: Histograms not found for bin group. Skipping plotting for this bin group." << std::endl;
            ++pad_idx;
            continue;
        }

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
    for (auto& entry : h_data_histograms) delete entry.second;
    for (auto& entry : h_mc_gen_histograms) delete entry.second;
    for (auto& entry : h_mc_rec_histograms) delete entry.second;
    delete canvas;
}