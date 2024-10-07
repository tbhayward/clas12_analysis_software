#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>  // For adding text in each subplot
#include <cmath>     // For conversion from radians to degrees
#include <string>
#include <vector>
#include "bin_boundaries.h"  // Include the header where BinBoundary is defined

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Function to calculate the number of Q²-t bins for the current xB bin
int count_Q2t_bins_for_xB(int xB_bin, const std::vector<BinBoundary>& bin_boundaries) {
    int n_Q2t_bins = 0;
    for (const auto& bin : bin_boundaries) {
        if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
            n_Q2t_bins++;
        }
    }
    return n_Q2t_bins;
}

// Plot function for DVCS data/MC comparison
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, const std::vector<BinBoundary>& bin_boundaries, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {

    // Count the number of Q²-t bins for the current xB bin
    int n_Q2t_bins = count_Q2t_bins_for_xB(xB_bin, bin_boundaries);

    // Calculate the number of subplots and make the canvas as square as possible
    int n_subplots = std::pow(std::ceil(std::sqrt(n_Q2t_bins)), 2);  // Perfect square of subplots
    int n_columns = std::ceil(std::sqrt(n_subplots));  // Rows and columns for a square layout
    int n_rows = n_columns;

    // Create canvas with dynamic subdivision
    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 1200, 800);
    canvas->Divide(n_columns, n_rows);  // Dynamically divide the canvas

    // Disable stat boxes globally
    gStyle->SetOptStat(0);

    // Create histograms for each Q²-t bin
    std::vector<TH1D*> h_data_histograms;
    std::vector<TH1D*> h_mc_gen_histograms;
    std::vector<TH1D*> h_mc_rec_histograms;

    // Loop through the bin boundaries and create histograms for each subplot
    int histogram_idx = 0;
    for (const auto& bin : bin_boundaries) {
        if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
            h_data_histograms.push_back(new TH1D(Form("h_data_%d", histogram_idx), "Reconstructed Data", 24, 0, 360));
            h_mc_gen_histograms.push_back(new TH1D(Form("h_mc_gen_%d", histogram_idx), "Generated MC", 24, 0, 360));
            h_mc_rec_histograms.push_back(new TH1D(Form("h_mc_rec_%d", histogram_idx), "Reconstructed MC", 24, 0, 360));
            histogram_idx++;
        }
    }

    // Restart the readers before looping over data
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Reinitialize readers for each tree
    TTreeReaderValue<double> phi_data(data_reader, "phi");
    TTreeReaderValue<double> xB_data(data_reader, "x");  // Use 'x' for xB
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");  // Use t1 instead of t

    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");

    // Fill histograms with data and MC according to cuts
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;  // Convert phi from radians to degrees

        histogram_idx = 0;
        for (const auto& bin : bin_boundaries) {
            if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
                if (*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                    *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                    std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high) {
                    h_data_histograms[histogram_idx]->Fill(phi_deg);
                }
                histogram_idx++;
            }
        }
    }

    while (mc_gen_reader.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;

        histogram_idx = 0;
        for (const auto& bin : bin_boundaries) {
            if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
                h_mc_gen_histograms[histogram_idx]->Fill(phi_mc_gen_deg);
                histogram_idx++;
            }
        }
    }

    while (mc_rec_reader.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;

        histogram_idx = 0;
        for (const auto& bin : bin_boundaries) {
            if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
                h_mc_rec_histograms[histogram_idx]->Fill(phi_mc_rec_deg);
                histogram_idx++;
            }
        }
    }

    // Normalize histograms and plot in each subplot
    histogram_idx = 0;
    for (int subplot_idx = 1; subplot_idx <= n_Q2t_bins; ++subplot_idx) {
        canvas->cd(subplot_idx);

        TH1D* h_data = h_data_histograms[histogram_idx];
        TH1D* h_mc_gen = h_mc_gen_histograms[histogram_idx];
        TH1D* h_mc_rec = h_mc_rec_histograms[histogram_idx];

        // Normalize histograms if they are not empty
        if (h_data->Integral() > 0) h_data->Scale(1.0 / h_data->Integral());
        if (h_mc_gen->Integral() > 0) h_mc_gen->Scale(1.0 / h_mc_gen->Integral());
        if (h_mc_rec->Integral() > 0) h_mc_rec->Scale(1.0 / h_mc_rec->Integral());

        // Find the maximum bin content across all histograms
        double max_value = std::max({h_data->GetMaximum(), h_mc_gen->GetMaximum(), h_mc_rec->GetMaximum()});

        // Set the y-axis range from 0 to 1.2 times the maximum value
        h_data->SetMaximum(1.2 * max_value);
        h_data->SetMinimum(0);

        // Set colors and styles
        h_data->SetMarkerColor(kBlue);
        h_data->SetMarkerStyle(20);  // Data points with error bars
        h_data->SetLineColor(kBlue);

        h_mc_gen->SetLineColor(kBlack);
        h_mc_gen->SetLineStyle(2);  // Black line for generated MC

        h_mc_rec->SetMarkerColor(kRed);
        h_mc_rec->SetMarkerStyle(22);  // Red points for reconstructed MC
        h_mc_rec->SetLineColor(kRed);

        // Draw histograms on the same canvas for each subplot
        h_data->Draw("E1");           // Data with error bars
        h_mc_rec->Draw("E1 SAME");    // Reconstructed MC with error bars
        h_mc_gen->Draw("HIST SAME");  // Generated MC as a line

        // Add kinematic constraints as text in each subplot
        TLatex latex;
        latex.SetTextSize(0.03);
        latex.DrawLatexNDC(0.15, 0.85, Form("x_{B}: %.2f - %.2f", bin_boundaries[histogram_idx].xB_low, bin_boundaries[histogram_idx].xB_high));
        latex.DrawLatexNDC(0.15, 0.8, Form("Q^{2}: %.2f - %.2f", bin_boundaries[histogram_idx].Q2_low, bin_boundaries[histogram_idx].Q2_high));
        latex.DrawLatexNDC(0.15, 0.75, Form("|t|: %.2f - %.2f", std::abs(bin_boundaries[histogram_idx].t_low), std::abs(bin_boundaries[histogram_idx].t_high)));

        // Add legend to every subplot
        TLegend* legend = new TLegend(0.5, 0.6, 0.9, 0.9);
        legend->AddEntry(h_data, "Data", "lep");
        legend->AddEntry(h_mc_rec, "Reconstructed MC", "lep");
        legend->AddEntry(h_mc_gen, "Generated MC", "l");
        legend->SetTextSize(0.04); 
        legend->Draw();

        histogram_idx++;
    }

    // Save canvas to the output directory
    std::string filename = output_dir + "/phi_data_mc_comparison_xB_bin_" + std::to_string(xB_bin) + ".png";
    canvas->SaveAs(filename.c_str());

    // Clean up histograms and canvas
    for (auto& h : h_data_histograms) delete h;
    for (auto& h : h_mc_gen_histograms) delete h;
    for (auto& h : h_mc_rec_histograms) delete h;
    delete canvas;
}