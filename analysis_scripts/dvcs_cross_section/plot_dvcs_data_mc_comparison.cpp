#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

// Plot function for DVCS data/MC comparison
void plot_dvcs_data_mc_comparison(const std::string& output_dir, int xB_bin, TTreeReader& data_reader, TTreeReader& mc_gen_reader, TTreeReader& mc_rec_reader) {
    // Create a canvas and histograms for plotting with 24 phi bins (15 degrees each)
    TCanvas* canvas = new TCanvas("c1", "Data vs MC", 800, 600);
    TH1D* h_data = new TH1D("h_data", "Reconstructed Data", 24, 0, 360);  // 24 bins for phi
    TH1D* h_mc_gen = new TH1D("h_mc_gen", "Generated MC", 24, 0, 360);
    TH1D* h_mc_rec = new TH1D("h_mc_rec", "Reconstructed MC", 24, 0, 360);

    std::cout << "We entered plot_dvcs..." << std::endl;

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

    // Draw histograms as points with error bars for data, lines for MC
    h_data->Draw("E1");         // Data with error bars
    h_mc_gen->Draw("HIST SAME");  // Generated MC as a line
    h_mc_rec->Draw("HIST SAME");  // Reconstructed MC as a line

    // Add legend
    TLegend* legend = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend->AddEntry(h_data, "Data", "lep");
    legend->AddEntry(h_mc_gen, "Generated MC", "l");
    legend->AddEntry(h_mc_rec, "Reconstructed MC", "l");
    legend->Draw();

    // Save canvas to the output directory
    std::string filename = output_dir + "/xB_bin_" + std::to_string(xB_bin) + ".png";
    canvas->SaveAs(filename.c_str());

    // Clean up
    delete h_data;
    delete h_mc_gen;
    delete h_mc_rec;
    delete legend;
    delete canvas;
}