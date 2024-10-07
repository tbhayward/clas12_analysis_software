#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>

// Plot function for DVCS data/MC comparison
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {
    // Create canvas and histograms for plotting
    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 800, 600);
    TH1D* h_data = new TH1D("h_data", "Reconstructed Data", 10, 0, 360);
    TH1D* h_mc_gen = new TH1D("h_mc_gen", "Generated MC", 10, 0, 360);
    TH1D* h_mc_rec = new TH1D("h_mc_rec", "Reconstructed MC", 10, 0, 360);

    // Loop over the tree and fill histograms
    TTreeReaderValue<double> phi(data_reader, "phi");
    while (data_reader.Next()) {
        h_data->Fill(*phi);
    }
    
    TTreeReaderValue<double> phi_gen(mc_gen_reader, "phi");
    while (mc_gen_reader.Next()) {
        h_mc_gen->Fill(*phi_gen);
    }
    
    TTreeReaderValue<double> phi_rec(mc_rec_reader, "phi");
    while (mc_rec_reader.Next()) {
        h_mc_rec->Fill(*phi_rec);
    }

    // Normalize the histograms to their integral
    if (h_data->Integral() > 0) h_data->Scale(1.0 / h_data->Integral());
    if (h_mc_gen->Integral() > 0) h_mc_gen->Scale(1.0 / h_mc_gen->Integral());
    if (h_mc_rec->Integral() > 0) h_mc_rec->Scale(1.0 / h_mc_rec->Integral());

    // Set colors and styles
    h_data->SetLineColor(kBlue);
    h_mc_gen->SetLineColor(kRed);
    h_mc_rec->SetLineColor(kGreen);

    // Draw histograms
    h_data->Draw("hist");
    h_mc_gen->Draw("hist same");
    h_mc_rec->Draw("hist same");

    // Add legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h_data, "Data", "l");
    legend->AddEntry(h_mc_gen, "Generated MC", "l");
    legend->AddEntry(h_mc_rec, "Reconstructed MC", "l");
    legend->Draw();

    // Save canvas
    canvas->SaveAs((output_dir + "/xB_bin_" + std::to_string(xB_bin) + ".png").c_str());

    // Clean up
    delete h_data;
    delete h_mc_gen;
    delete h_mc_rec;
    delete canvas;
}