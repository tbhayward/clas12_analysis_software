#include "plot_pi0_mass.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <TF1.h>
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
    TTreeReaderValue<double> Mh_data1(dataReader1, "Mh_gammagamma");
    TTreeReaderValue<double> Mh_mc1(mcReader1, "Mh_gammagamma");
    TTreeReaderValue<double> Mh_data2(dataReader2, "Mh_gammagamma");
    TTreeReaderValue<double> Mh_mc2(mcReader2, "Mh_gammagamma");
    TTreeReaderValue<double> Mh_data3(dataReader3, "Mh_gammagamma");
    TTreeReaderValue<double> Mh_mc3(mcReader3, "Mh_gammagamma");

    // Fill histograms for each data reader
    while (dataReader1.Next()) {
        hist_data1->Fill(*Mh_data1);
    }
    while (mcReader1.Next()) {
        hist_mc1->Fill(*Mh_mc1);
    }

    while (dataReader2.Next()) {
        hist_data2->Fill(*Mh_data2);
    }
    while (mcReader2.Next()) {
        hist_mc2->Fill(*Mh_mc2);
    }

    while (dataReader3.Next()) {
        hist_data3->Fill(*Mh_data3);
    }
    while (mcReader3.Next()) {
        hist_mc3->Fill(*Mh_mc3);
    }

    // Normalize the histograms based on their integrals
    if (hist_data1->Integral() != 0) hist_data1->Scale(1.0 / hist_data1->Integral());
    if (hist_mc1->Integral() != 0) hist_mc1->Scale(1.0 / hist_mc1->Integral());
    if (hist_data2->Integral() != 0) hist_data2->Scale(1.0 / hist_data2->Integral());
    if (hist_mc2->Integral() != 0) hist_mc2->Scale(1.0 / hist_mc2->Integral());
    if (hist_data3->Integral() != 0) hist_data3->Scale(1.0 / hist_data3->Integral());
    if (hist_mc3->Integral() != 0) hist_mc3->Scale(1.0 / hist_mc3->Integral());

    // Create dynamic lines for each plot
    TLine* pi0_mass_line1 = new TLine(0.135, 0, 0.135, hist_data1->GetMaximum() * 1.1);
    TLine* pi0_mass_line2 = new TLine(0.135, 0, 0.135, hist_data2->GetMaximum() * 1.1);
    TLine* pi0_mass_line3 = new TLine(0.135, 0, 0.135, hist_data3->GetMaximum() * 1.1);
    pi0_mass_line1->SetLineColor(kGray + 2);
    pi0_mass_line1->SetLineStyle(7);
    pi0_mass_line2->SetLineColor(kGray + 2);
    pi0_mass_line2->SetLineStyle(7);
    pi0_mass_line3->SetLineColor(kGray + 2);
    pi0_mass_line3->SetLineStyle(7);

    // Fit histograms with Gaussian + constant model
    TF1* fit_data1 = new TF1("fit_data1", "gaus(0) + pol0(3)", 0.11, 0.16);
    TF1* fit_mc1 = new TF1("fit_mc1", "gaus(0) + pol0(3)", 0.11, 0.16);
    TF1* fit_data2 = new TF1("fit_data2", "gaus(0) + pol0(3)", 0.11, 0.16);
    TF1* fit_mc2 = new TF1("fit_mc2", "gaus(0) + pol0(3)", 0.11, 0.16);
    TF1* fit_data3 = new TF1("fit_data3", "gaus(0) + pol0(3)", 0.11, 0.16);
    TF1* fit_mc3 = new TF1("fit_mc3", "gaus(0) + pol0(3)", 0.11, 0.16);

    // Set initial parameter guesses
    fit_data1->SetParameters(0.02, 0.135, 0.01);
    fit_mc1->SetParameters(0.02, 0.135, 0.01);
    fit_data2->SetParameters(0.02, 0.135, 0.01);
    fit_mc2->SetParameters(0.02, 0.135, 0.01);
    fit_data3->SetParameters(0.02, 0.135, 0.01);
    fit_mc3->SetParameters(0.02, 0.135, 0.01);

    // Perform fits
    canvas->cd(1);
    hist_data1->Fit(fit_data1, "R");
    hist_mc1->Fit(fit_mc1, "R");
    pi0_mass_line1->Draw("SAME");

    canvas->cd(2);
    hist_data2->Fit(fit_data2, "R");
    hist_mc2->Fit(fit_mc2, "R");
    pi0_mass_line2->Draw("SAME");

    canvas->cd(3);
    hist_data3->Fit(fit_data3, "R");
    hist_mc3->Fit(fit_mc3, "R");
    pi0_mass_line3->Draw("SAME");

    // Get fit parameters
    double mu_data1 = fit_data1->GetParameter(1);
    double sigma_data1 = fit_data1->GetParameter(2);
    double mu_mc1 = fit_mc1->GetParameter(1);
    double sigma_mc1 = fit_mc1->GetParameter(2);

    double mu_data2 = fit_data2->GetParameter(1);
    double sigma_data2 = fit_data2->GetParameter(2);
    double mu_mc2 = fit_mc2->GetParameter(1);
    double sigma_mc2 = fit_mc2->GetParameter(2);

    double mu_data3 = fit_data3->GetParameter(1);
    double sigma_data3 = fit_data3->GetParameter(2);
    double mu_mc3 = fit_mc3->GetParameter(1);
    double sigma_mc3 = fit_mc3->GetParameter(2);

    // Create legends for the three subplots
    TLegend* legend1 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend1->AddEntry(hist_data1, Form("#color[4]{Data},  #mu = %.3f,  #sigma = %.3f", mu_data1, sigma_data1), "p");
    legend1->AddEntry(hist_mc1, Form("#color[2]{MC},  #mu = %.3f,  #sigma = %.3f", mu_mc1, sigma_mc1), "p");
    legend1->Draw();

    TLegend* legend2 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend2->AddEntry(hist_data2, Form("#color[4]{Data},  #mu = %.3f,  #sigma = %.3f", mu_data2, sigma_data2), "p");
    legend2->AddEntry(hist_mc2, Form("#color[2]{MC},  #mu = %.3f,  #sigma = %.3f", mu_mc2, sigma_mc2), "p");
    legend2->Draw();

    TLegend* legend3 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend3->AddEntry(hist_data3, Form("#color[4]{Data},  #mu = %.3f,  #sigma = %.3f", mu_data3, sigma_data3), "p");
    legend3->AddEntry(hist_mc3, Form("#color[2]{MC},  #mu = %.3f,  #sigma = %.3f", mu_mc3, sigma_mc3), "p");
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
    delete pi0_mass_line1;
    delete pi0_mass_line2;
    delete pi0_mass_line3;
    delete canvas;
}