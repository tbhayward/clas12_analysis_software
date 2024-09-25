#include "plot_pi0_mass.h"
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReaderValue.h>
#include <TLine.h>
#include <filesystem>
#include <iostream>
#include <TF1.h>

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

    // Create a dashed gray line at the pi0 mass (0.135 GeV)
    TLine* pi0_mass_line = new TLine(0.135, 0, 0.135, 0.25);
    pi0_mass_line->SetLineColor(kGray + 2);
    pi0_mass_line->SetLineStyle(7);  // Dashed line

    // Determine the maximum value for y-axis scaling
    double y_max1 = 1.1 * std::max(hist_data1->GetMaximum(), hist_mc1->GetMaximum());
    double y_max2 = 1.1 * std::max(hist_data2->GetMaximum(), hist_mc2->GetMaximum());
    double y_max3 = 1.1 * std::max(hist_data3->GetMaximum(), hist_mc3->GetMaximum());

    // Draw the histograms on the canvas as points with error bars
    canvas->cd(1);
    hist_data1->SetLineColor(kBlue);
    hist_data1->SetMarkerColor(kBlue);
    hist_data1->SetMarkerStyle(20);  // Points with error bars
    hist_mc1->SetLineColor(kRed);
    hist_mc1->SetMarkerColor(kRed);
    hist_mc1->SetMarkerStyle(24);
    hist_data1->SetXTitle("M_{#gamma#gamma} (GeV)");
    hist_data1->GetYaxis()->SetRangeUser(0, y_max1);

    // Draw the data and MC histograms
    hist_data1->Draw("E1");
    hist_mc1->Draw("E1 SAME");

    // Create and draw the vertical line at pi0 mass
    pi0_mass_line->Draw("SAME");

    // Create a Gaussian function for fitting the data histogram within the range 0.12 to 0.15 GeV
    TF1* gausFit = new TF1("gausFit", "gaus", 0.12, 0.15);
    gausFit->SetLineColor(kBlue);  // Set the line color to blue
    gausFit->SetLineWidth(2);      // Set the line width for better visibility

    // Fit the data histogram with the Gaussian function
    hist_data1->Fit(gausFit, "R");  // "R" ensures the fit is within the specified range
    double mu = gausFit->GetParameter(1);      // Mean (μ)
    double sigma = gausFit->GetParameter(2);   // Sigma (σ)

    // Create a Gaussian function for fitting the MC histogram within the same range
    TF1* gausFitMC = new TF1("gausFitMC", "gaus", 0.12, 0.15);
    gausFitMC->SetLineColor(kRed);  // Set the line color to red
    gausFitMC->SetLineWidth(2);     // Set the line width for better visibility

    // Fit the MC histogram with the Gaussian function
    hist_mc1->Fit(gausFitMC, "R");  // "R" ensures the fit is within the specified range
    double muMC = gausFitMC->GetParameter(1);      // Mean (μ) for MC
    double sigmaMC = gausFitMC->GetParameter(2);   // Sigma (σ) for MC

    // Add legend for the first plot with colored text
    TLegend* legend1 = new TLegend(0.65, 0.7, 0.9, 0.9);  // Adjusted position to accommodate more entries
    legend1->SetTextSize(0.04);  // Increase the text size slightly

    char dataLegendEntry[200];
    sprintf(dataLegendEntry, "#color[4]{Data (#mu = %.4f, #sigma = %.4f)}", mu, sigma);
    legend1->AddEntry(hist_data1, dataLegendEntry, "p");

    char mcLegendEntry[200];
    sprintf(mcLegendEntry, "#color[2]{MC (#mu = %.4f, #sigma = %.4f)}", muMC, sigmaMC);
    legend1->AddEntry(hist_mc1, mcLegendEntry, "p");

    legend1->Draw();

    canvas->cd(2);
    hist_data2->SetLineColor(kBlue);
    hist_data2->SetMarkerColor(kBlue);
    hist_data2->SetMarkerStyle(20);
    hist_mc2->SetLineColor(kRed);
    hist_mc2->SetMarkerColor(kRed);
    hist_mc2->SetMarkerStyle(24);
    hist_data2->SetXTitle("M_{#gamma#gamma} (GeV)");
    hist_data2->GetYaxis()->SetRangeUser(0, y_max2);
    hist_data2->Draw("E1");
    hist_mc2->Draw("E1 SAME");
    pi0_mass_line->Draw("SAME");

    // Add legend for the second plot with colored text
    TLegend* legend2 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend2->AddEntry(hist_data2, "#color[4]{Data}", "p");
    legend2->AddEntry(hist_mc2, "#color[2]{MC}", "p");
    legend2->Draw();

    canvas->cd(3);
    hist_data3->SetLineColor(kBlue);
    hist_data3->SetMarkerColor(kBlue);
    hist_data3->SetMarkerStyle(20);
    hist_mc3->SetLineColor(kRed);
    hist_mc3->SetMarkerColor(kRed);
    hist_mc3->SetMarkerStyle(24);
    hist_data3->SetXTitle("M_{#gamma#gamma} (GeV)");
    hist_data3->GetYaxis()->SetRangeUser(0, y_max3);
    hist_data3->Draw("E1");
    hist_mc3->Draw("E1 SAME");
    pi0_mass_line->Draw("SAME");

    // Add legend for the third plot with colored text
    TLegend* legend3 = new TLegend(0.7, 0.75, 0.9, 0.9);
    legend3->AddEntry(hist_data3, "#color[4]{Data}", "p");
    legend3->AddEntry(hist_mc3, "#color[2]{MC}", "p");
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
    delete pi0_mass_line;
    delete canvas;
}