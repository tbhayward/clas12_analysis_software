#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "bin_boundaries.h"  // Include the header where BinBoundary is defined
#include <cmath>  // For conversion from radians to degrees and std::abs

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Plot function for DVCS data/MC comparison with dynamic subplot layout
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, const std::vector<BinBoundary>& bin_boundaries, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {

    // Get the number of Q²-t bins for the current xB bin
    int n_Q2t_bins = 0;
    for (const auto& bin : bin_boundaries) {
        if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
            n_Q2t_bins++;
        }
    }

    // Create canvas with dynamic number of subplots
    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 1200, 800);
    canvas->Divide((n_Q2t_bins + 2) / 3, 3);  // Create enough subplots

    // Disable stat boxes globally
    gStyle->SetOptStat(0);

    // Loop over the first Q²-t bin for debugging (single subplot)
    int subplot_idx = 1;
    // for (int i = 0; i < n_Q2t_bins; ++i) {
    for (int i = 0; i < 1; ++i) {  // Single bin for now
        const BinBoundary& bin = bin_boundaries[i];

        // Ensure the bin corresponds to the current xB_bin
        if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {
            // Create histograms for the current Q²-t bin
            TH1D* h_data = new TH1D("h_data", "Reconstructed Data", 24, 0, 360);  // 24 bins for phi (15-degree bins)
            TH1D* h_mc_gen = new TH1D("h_mc_gen", "Generated MC", 24, 0, 360);
            TH1D* h_mc_rec = new TH1D("h_mc_rec", "Reconstructed MC", 24, 0, 360);

            // Restart the readers before looping over data
            data_reader.Restart();
            mc_gen_reader.Restart();
            mc_rec_reader.Restart();

            // Reinitialize readers before looping over the data
            TTreeReaderValue<double> phi_data(data_reader, "phi");
            TTreeReaderValue<double> Q2_data(data_reader, "Q2");
            TTreeReaderValue<double> t1_data(data_reader, "t1");  // Now using t1 instead of t
            TTreeReaderValue<double> xB_data(data_reader, "x");  // xB is the x branch

            TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
            TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
            TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");  // Now using t1 instead of t
            TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");

            TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");
            TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
            TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");  // Now using t1 instead of t
            TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");

            // Read data from trees and fill histograms with checks for xB, Q², and absolute value of t1
            while (data_reader.Next()) {
                // Convert phi from radians to degrees
                double phi_deg = *phi_data * RAD_TO_DEG;

                if (*xB_data >= bin.xB_low && *xB_data <= bin.xB_high &&
                    *Q2_data >= bin.Q2_low && *Q2_data <= bin.Q2_high &&
                    std::abs(*t1_data) >= bin.t_low && std::abs(*t1_data) <= bin.t_high) {  // Use absolute value of t1
                    h_data->Fill(phi_deg);
                }
            }
            while (mc_gen_reader.Next()) {
                double phi_mc_gen_deg = *phi_mc_gen * RAD_TO_DEG;

                if (*xB_mc_gen >= bin.xB_low && *xB_mc_gen <= bin.xB_high &&
                    *Q2_mc_gen >= bin.Q2_low && *Q2_mc_gen <= bin.Q2_high &&
                    std::abs(*t1_mc_gen) >= bin.t_low && std::abs(*t1_mc_gen) <= bin.t_high) {  // Use absolute value of t1
                    h_mc_gen->Fill(phi_mc_gen_deg);
                }
            }
            while (mc_rec_reader.Next()) {
                double phi_mc_rec_deg = *phi_mc_rec * RAD_TO_DEG;

                if (*xB_mc_rec >= bin.xB_low && *xB_mc_rec <= bin.xB_high &&
                    *Q2_mc_rec >= bin.Q2_low && *Q2_mc_rec <= bin.Q2_high &&
                    std::abs(*t1_mc_rec) >= bin.t_low && std::abs(*t1_mc_rec) <= bin.t_high) {  // Use absolute value of t1
                    h_mc_rec->Fill(phi_mc_rec_deg);
                }
            }

            // Check if h_data and h_mc_rec are both empty, and skip this subplot if so
            if (h_data->Integral() == 0 && h_mc_rec->Integral() == 0) {
                delete h_data;
                delete h_mc_gen;
                delete h_mc_rec;
                continue;
            }

            // Normalize the histograms to their integrals
            if (h_data->Integral() > 0) h_data->Scale(1.0 / h_data->Integral());
            if (h_mc_gen->Integral() > 0) h_mc_gen->Scale(1.0 / h_mc_gen->Integral());
            if (h_mc_rec->Integral() > 0) h_mc_rec->Scale(1.0 / h_mc_rec->Integral());

            // Set colors and styles for each histogram
            h_data->SetMarkerColor(kBlue);
            h_data->SetMarkerStyle(20);  // Data points with error bars
            h_data->SetLineColor(kBlue);

            h_mc_gen->SetLineColor(kRed);
            h_mc_gen->SetLineStyle(2);  // Dotted line for generated MC

            h_mc_rec->SetLineColor(kGreen);
            h_mc_rec->SetLineStyle(3);  // Dashed line for reconstructed MC

            // Select the next available subplot and activate the pad
            canvas->cd(subplot_idx);
            TPad* pad = (TPad*)canvas->cd(subplot_idx);
            pad->SetLogy();  // Optional: Enable log scale for better visibility

            // Ensure histograms are drawn in the correct order
            h_data->Draw("E1");           // Data with error bars
            h_mc_gen->Draw("HIST SAME");  // Generated MC as a line
            h_mc_rec->Draw("HIST SAME");  // Reconstructed MC as a line

            // Add legend in the first subplot
            if (subplot_idx == 1) {
                TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.9);
                legend->AddEntry(h_data, "Data", "lep");
                legend->AddEntry(h_mc_gen, "Generated MC", "l");
                legend->AddEntry(h_mc_rec, "Reconstructed MC", "l");
                legend->Draw();
            }

            // Move to the next subplot
            subplot_idx++;

            // Clean up histograms after drawing
            delete h_data;
            delete h_mc_gen;
            delete h_mc_rec;
        }
    }

    // Save canvas to the output directory
    std::string filename = output_dir + "/xB_bin_" + std::to_string(xB_bin) + ".png";
    canvas->SaveAs(filename.c_str());

    // Clean up canvas
    delete canvas;
}