#include "plot_pi0_mass.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>
#include <filesystem>
#include <iostream>

void plot_pi0_mass(TTreeReader& dataReader1, TTreeReader& dataReader2, TTreeReader& dataReader3,
                   TTreeReader& mcReader1, TTreeReader& mcReader2, TTreeReader& mcReader3,
                   const std::string& outputDir) {
    // Set up global style options (remove stat boxes)
    gStyle->SetOptStat(0);

    // Check and create the necessary directory if it doesn't exist
    std::string pi0_mass_dir = outputDir + "/pi0_mass";
    if (!std::filesystem::exists(pi0_mass_dir)) {
        std::filesystem::create_directories(pi0_mass_dir);
    }

    // Create a canvas with 1x3 layout (wider aspect ratio)
    TCanvas* canvas = new TCanvas("pi0_mass_canvas", "Pi0 Mass Plots", 2400, 800);  // Wider canvas
    canvas->Divide(3, 1);  // 1 row, 3 columns

    // Create histograms for the data and MC
    TH1D* hist_data1 = new TH1D("hist_data1", "Fa18 Inb", 100, 0.11, 0.16);
    TH1D* hist_mc1 = new TH1D("hist_mc1", "Fa18 Inb", 100, 0.11, 0.16);
    TH1D* hist_data2 = new TH1D("hist_data2", "Fa18 Out", 100, 0.11, 0.16);
    TH1D* hist_mc2 = new TH1D("hist_mc2", "Fa18 Out", 100, 0.11, 0.16);
    TH1D* hist_data3 = new TH1D("hist_data3", "Sp19 Inb", 100, 0.11, 0.16);
    TH1D* hist_mc3 = new TH1D("hist_mc3", "Sp19 Inb", 100, 0.11, 0.16);

    // Set directories to 0 to avoid automatic ROOT management
    hist_data1->SetDirectory(0);
    hist_mc1->SetDirectory(0);
    hist_data2->SetDirectory(0);
    hist_mc2->SetDirectory(0);
    hist_data3->SetDirectory(0);
    hist_mc3->SetDirectory(0);

    // Readers for Mh variable
    TTreeReaderValue<double> Mh_data1(dataReader1, "Mh");
    TTreeReaderValue<double> Mh_mc1(mcReader1, "Mh");
    TTreeReaderValue<double> Mh_data2(dataReader2, "Mh");
    TTreeReaderValue<double> Mh_mc2(mcReader2, "Mh");
    TTreeReaderValue<double> Mh_data3(dataReader3, "Mh");
    TTreeReaderValue<double> Mh_mc3(mcReader3, "Mh");

    // Fill histograms for each data reader
    while (dataReader1.Next()) {
        if (*Mh_data1 >= 0.11 && *Mh_data1 <= 0.16) {
            hist_data1->Fill(*Mh_data1);
        }
    }
    while (mcReader1.Next()) {
        if (*Mh_mc1 >= 0.11 && *Mh_mc1 <= 0.16) {
            hist_mc1->Fill(*Mh_mc1);
        }
    }

    while (dataReader2.Next()) {
        if (*Mh_data2 >= 0.11 && *Mh_data2 <= 0.16) {
            hist_data2->Fill(*Mh_data2);
        }
    }
    while (mcReader2.Next()) {
        if (*Mh_mc2 >= 0.11 && *Mh_mc2 <= 0.16) {
            hist_mc2->Fill(*Mh_mc2);
        }
    }

    while (dataReader3.Next()) {
        if (*Mh_data3 >= 0.11 && *Mh_data3 <= 0.16) {
            hist_data3->Fill(*Mh_data3);
        }
    }
    while (mcReader3.Next()) {
        if (*Mh_mc3 >= 0.11 && *Mh_mc3 <= 0.16) {
            hist_mc3->Fill(*Mh_mc3);
        }
    }

    // Normalize the histograms based on their integrals
    if (hist_data1->Integral() != 0) hist_data1->Scale(1.0 / hist_data1->Integral());
    if (hist_mc1->Integral() != 0) hist_mc1->Scale(1.0 / hist_mc1->Integral());
    if (hist_data2->Integral() != 0) hist_data2->Scale(1.0 / hist_data2->Integral());
    if (hist_mc2->Integral() != 0) hist_mc2->Scale(1.0 / hist_mc2->Integral());
    if (hist_data3->Integral() != 0) hist_data3->Scale(1.0 / hist_data3->Integral());
    if (hist_mc3->Integral() != 0) hist_mc3->Scale(1.0 / hist_mc3->Integral());

    // Draw the histograms on the canvas
    canvas->cd(1);
    hist_data1->SetLineColor(kBlue);
    hist_mc1->SetLineColor(kRed);
    hist_data1->SetXTitle("M_{#gamma#gamma} (GeV)");
    hist_data1->Draw("HIST");
    hist_mc1->Draw("HIST SAME");

    // Add legend for the first plot
    TLegend* legend1 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend1->AddEntry(hist_data1, "Data", "l");
    legend1->AddEntry(hist_mc1, "MC", "l");
    legend1->Draw();

    canvas->cd(2);
    hist_data2->SetLineColor(kBlue);
    hist_mc2->SetLineColor(kRed);
    hist_data2->SetXTitle("M_{#gamma#gamma} (GeV)");
    hist_data2->Draw("HIST");
    hist_mc2->Draw("HIST SAME");

    // Add legend for the second plot
    TLegend* legend2 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend2->AddEntry(hist_data2, "Data", "l");
    legend2->AddEntry(hist_mc2, "MC", "l");
    legend2->Draw();

    canvas->cd(3);
    hist_data3->SetLineColor(kBlue);
    hist_mc3->SetLineColor(kRed);
    hist_data3->SetXTitle("M_{#gamma#gamma} (GeV)");
    hist_data3->Draw("HIST");
    hist_mc3->Draw("HIST SAME");

    // Add legend for the third plot
    TLegend* legend3 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend3->AddEntry(hist_data3, "Data", "l");
    legend3->AddEntry(hist_mc3, "MC", "l");
    legend3->Draw();

    // Save the canvas
    canvas->SaveAs((pi0_mass_dir + "/pi0_mass_comparison.png").c_str());

    // Clean up memory
    delete hist_data1;
    delete hist_mc1;
    delete hist_data2;
    delete hist_mc2;
    delete hist_data3;
    delete hist_mc3;
    delete canvas;
}