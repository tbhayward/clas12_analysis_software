#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>  // For conversion from radians to degrees
#include <string>
#include <vector>

// Include the header where BinBoundary is defined
#include "bin_boundaries.h"  

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Simplified Plot function for DVCS data/MC comparison
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, const std::vector<BinBoundary>& bin_boundaries, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {

    // The parameters xB_bin and bin_boundaries are currently unused in this simplified version.
    (void)xB_bin;  // Silence unused variable warning for now
    (void)bin_boundaries;

    // Create a single canvas (no subdivision)
    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 1200, 800);

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

    // Check if histograms are filled
    std::cout << "Data: " << h_data->Integral() << ", MC Gen: " << h_mc_gen->Integral() << ", MC Rec: " << h_mc_rec->Integral() << std::endl;

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

    h_mc_gen->SetMarkerColor(kRed);
    h_mc_gen->SetMarkerStyle(21);  // Generated MC points
    h_mc_gen->SetLineColor(kRed);  // Set line color in case of connecting dots later

    h_mc_rec->SetMarkerColor(kGreen);
    h_mc_rec->SetMarkerStyle(22);  // Reconstructed MC points
    h_mc_rec->SetLineColor(kGreen);  // Set line color in case of connecting dots later

    // Draw histograms on the same canvas
    h_data->Draw("E1 P");   // Draw data with error bars as points
    h_mc_gen->Draw("Hist SAME");  // Draw generated MC as points
    h_mc_rec->Draw("Hist SAME");  // Draw reconstructed MC as points

    // Set axis labels
    h_data->GetXaxis()->SetTitle("#phi");
    h_data->GetYaxis()->SetTitle("Normalized Counts");

    // Create and add a legend
    TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend->AddEntry(h_data, "Data", "lep");
    legend->AddEntry(h_mc_gen, "Generated MC", "lep");
    legend->AddEntry(h_mc_rec, "Reconstructed MC", "lep");
    legend->Draw();

    // Save canvas to the output directory
    std::string filename = output_dir + "/phi_data_mc_comparison.png";
    canvas->SaveAs(filename.c_str());

    // Clean up
    delete h_data;
    delete h_mc_gen;
    delete h_mc_rec;
    delete legend;
    delete canvas;
}