#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cmath>  // For conversion from radians to degrees

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Simplified Plot function for DVCS data/MC comparison without bin boundary checks
void plot_dvcs_data_mc_comparison(const std::string& output_dir, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {

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

    // Set colors and styles
    h_data->SetMarkerColor(kBlue);
    h_data->SetMarkerStyle(20);  // Data points with error bars
    h_data->SetLineColor(kBlue);

    h_mc_gen->SetLineColor(kRed);
    h_mc_gen->SetLineStyle(2);  // Dotted line for generated MC

    h_mc_rec->SetLineColor(kGreen);
    h_mc_rec->SetLineStyle(3);  // Dashed line for reconstructed MC

    // Draw histograms on the same canvas
    h_data->Draw("E1");
    h_mc_gen->Draw("HIST SAME");
    h_mc_rec->Draw("HIST SAME");

    // Create and add a legend
    TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend->AddEntry(h_data, "Data", "lep");
    legend->AddEntry(h_mc_gen, "Generated MC", "l");
    legend->AddEntry(h_mc_rec, "Reconstructed MC", "l");
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