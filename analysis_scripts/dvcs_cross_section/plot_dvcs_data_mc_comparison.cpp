#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>  // For conversion from radians to degrees
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

// Simplified Plot function for DVCS data/MC comparison
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

    // Create histograms for data and MC
    TH1D* h_data = new TH1D("h_data", "Reconstructed Data", 24, 0, 360);  // 24 bins for phi (15-degree bins)
    TH1D* h_mc_gen = new TH1D("h_mc_gen", "Generated MC", 24, 0, 360);
    TH1D* h_mc_rec = new TH1D("h_mc_rec", "Reconstructed MC", 24, 0, 360);

    // Restart the readers before looping over data
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Reinitialize readers for each tree
    TTreeReaderValue<double> phi_data(data_reader, "phi");
    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");

    // Read data from trees and fill histograms
    while (data_reader.Next()) {
        double phi_deg = *phi_data * RAD_TO_DEG;  // Convert phi from radians to degrees
        h_data->Fill(phi_deg);
    }

    while (mc_gen_reader.Next()) {
        double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;  // Convert phi from radians to degrees
        h_mc_gen->Fill(phi_mc_gen_deg);
    }

    while (mc_rec_reader.Next()) {
        double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;  // Convert phi from radians to degrees
        h_mc_rec->Fill(phi_mc_rec_deg);
    }

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

    // Loop through each subplot and draw the histograms
    for (int i = 1; i <= n_subplots; ++i) {
        canvas->cd(i);  // Go to the next subplot

        // Draw histograms on the same canvas for each subplot
        h_data->Draw("E1");           // Data with error bars
        h_mc_gen->Draw("HIST SAME");  // Generated MC as a line
        h_mc_rec->Draw("E1 SAME");    // Reconstructed MC with error bars

        // Set axis labels in the first subplot
        if (i == 1) {
            h_data->GetXaxis()->SetTitle("#phi");
            h_data->GetYaxis()->SetTitle("Normalized Counts");
        }
    }

    // Create and add a legend to the first subplot
    canvas->cd(1);  // Ensure we're in the first subplot
    TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend->AddEntry(h_data, "Data", "lep");
    legend->AddEntry(h_mc_gen, "Generated MC", "l");
    legend->AddEntry(h_mc_rec, "Reconstructed MC", "lep");
    legend->Draw();

    // Save canvas to the output directory
    std::string filename = output_dir + "/phi_data_mc_comparison_xB_bin_" + std::to_string(xB_bin) + ".png";
    canvas->SaveAs(filename.c_str());

    // Clean up
    delete h_data;
    delete h_mc_gen;
    delete h_mc_rec;
    delete legend;
    delete canvas;
}