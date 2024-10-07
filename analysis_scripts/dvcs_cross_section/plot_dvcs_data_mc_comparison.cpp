#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include "bin_boundaries.h"  // Include the header where BinBoundary is defined

// Plot function for DVCS data/MC comparison with dynamic subplot layout
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, const std::vector<BinBoundary>& bin_boundaries, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader);
    std::cout << "We entered plot_dvcs..." << std::endl;

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

    // Loop over all Q²-t bins for this xB bin
    int subplot_idx = 1;
    for (size_t i = 0; i < bin_boundaries.size(); ++i) {
        const BinBoundary& bin = bin_boundaries[i];

        // Check if the bin belongs to the current xB_bin
        if (bin.xB_low == bin_boundaries[xB_bin].xB_low && bin.xB_high == bin_boundaries[xB_bin].xB_high) {

            // Create histograms for the current Q²-t bin
            TH1D* h_data = new TH1D("h_data", "Reconstructed Data", 24, 0, 360);  // 24 bins for phi
            TH1D* h_mc_gen = new TH1D("h_mc_gen", "Generated MC", 24, 0, 360);
            TH1D* h_mc_rec = new TH1D("h_mc_rec", "Reconstructed MC", 24, 0, 360);

            // Read data from trees and fill histograms
            TTreeReaderValue<double> phi_data(data_reader, "phi");
            while (data_reader.Next()) {
                h_data->Fill(*phi_data);
            }

            TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi");
            while (mc_gen_reader.Next()) {
                h_mc_gen->Fill(*phi_mc_gen);
            }

            TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi");
            while (mc_rec_reader.Next()) {
                h_mc_rec->Fill(*phi_mc_rec);
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

            // Select the next available subplot
            canvas->cd(subplot_idx++);
            h_data->Draw("E1");         // Data with error bars
            h_mc_gen->Draw("HIST SAME");  // Generated MC as a line
            h_mc_rec->Draw("HIST SAME");  // Reconstructed MC as a line

            // Add legend in the first subplot
            if (subplot_idx == 2) {
                TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.9);
                legend->AddEntry(h_data, "Data", "lep");
                legend->AddEntry(h_mc_gen, "Generated MC", "l");
                legend->AddEntry(h_mc_rec, "Reconstructed MC", "l");
                legend->Draw();
            }

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