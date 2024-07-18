#include <iostream>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <sstream>

void dilution_factors_epX(const char* nh3_file, const char* c_file) {
    // Open the ROOT files
    TFile *nh3 = TFile::Open(nh3_file);
    TFile *carbon = TFile::Open(c_file);
    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    // Get the PhysicsEvents trees
    TTree *tree_nh3;
    TTree *tree_carbon;
    nh3->GetObject("PhysicsEvents", tree_nh3);
    carbon->GetObject("PhysicsEvents", tree_carbon);
    if (!tree_nh3 || !tree_carbon) {
        std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
        nh3->Close();
        carbon->Close();
        return;
    }

    // Create histograms for xF2
    // TH1D *h_xF2_nh3 = new TH1D("h_xF2_nh3", "x_{F2} Distribution; x_{F2}; Counts", 100, -2.5, 1);
    // TH1D *h_xF2_carbon = new TH1D("h_xF2_carbon", "x_{F2} Distribution; x_{F2}; Counts", 100, -2.5, 1);
    TH1D *h_xF2_nh3 = new TH1D("h_xF2_nh3", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);
    TH1D *h_xF2_carbon = new TH1D("h_xF2_carbon", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);

    // Fill the histograms
    // tree_nh3->Draw("xF2>>h_xF2_nh3");
    // tree_carbon->Draw("xF2>>h_xF2_carbon");
    tree_nh3->Draw("Mx>>h_xF2_nh3");
    tree_carbon->Draw("Mx>>h_xF2_carbon");

    // Create canvas and divide it into four panels
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(2, 2);

    // First panel: plot xF2 histograms
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_xF2_nh3->SetLineColor(kBlue);
    h_xF2_carbon->SetLineColor(kRed);
    h_xF2_nh3->Draw();
    h_xF2_carbon->Draw("SAME");

    // Add legend
    TLegend *leg = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg->AddEntry(h_xF2_nh3, "NH_{3}", "l");
    leg->AddEntry(h_xF2_carbon, "C", "l");
    leg->Draw();

    // Remove statboxes
    h_xF2_nh3->SetStats(0);
    h_xF2_carbon->SetStats(0);

    // Second panel: ratio of NH3 to Carbon counts
    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_ratio = new TGraphErrors();
    for (int i = 1; i <= h_xF2_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xF2_nh3->GetBinContent(i);
        double c_counts = h_xF2_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio->SetPoint(i - 1, h_xF2_nh3->GetBinCenter(i), ratio);
            gr_ratio->SetPointError(i - 1, 0, error);
        }
    }
    // Set y-axis range from 5 to 15
    gr_ratio->GetYaxis()->SetRangeUser(9, 15);

    // gr_ratio->SetTitle("NH_{3} to Carbon Ratio; x_{F2}; Ratio");
    gr_ratio->SetTitle("NH_{3} to Carbon Ratio; M_{x} (GeV); Ratio");
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");

    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const = new TF1("fit_const", "[0]", -2, -0.5);
    gr_ratio->Fit(fit_const, "R");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");

    // Add the dotted-dashed line from -0.5 to 3
    TF1 *dotted_line = new TF1("dotted_line", "[0]", -0.5, 3);
    dotted_line->SetParameter(0, fit_const->GetParameter(0)); // Use the same constant value as the fit
    dotted_line->SetLineColor(kRed);
    dotted_line->SetLineStyle(7); // Set line style to dotted-dashed
    dotted_line->Draw("SAME");

    // Add fit constant value and uncertainty
    double fit_value = fit_const->GetParameter(0);
    double fit_error = fit_const->GetParError(0);

    // Retrieve chi2 and NDF
    double chi2 = fit_const->GetChisquare();
    int ndf = fit_const->GetNDF();
    double chi2_ndf = chi2 / ndf;

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.85, Form("Fit Const, s = %.3f #pm %.3f", fit_value, fit_error));
    latex.DrawLatex(0.20, 0.80, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2_ndf));

    // Save the canvas
    c1->SaveAs("dilution_factors.pdf");
    // Clean up
    nh3->Close();
    carbon->Close();
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file>" << std::endl;
        return 1;
    }

    dilution_factors_epX(argv[1], argv[2]);
    return 0;
}