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

void dilution_factors(const char* nh3_file, const char* c_file) {
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
    TH1D *h_xF2_nh3 = new TH1D("h_xF2_nh3", "xF2 Distribution; xF2; Counts", 100, -2.5, 1);
    TH1D *h_xF2_carbon = new TH1D("h_xF2_carbon", "xF2 Distribution; xF2; Counts", 100, -2.5, 1);

    // Fill the histograms
    tree_nh3->Draw("xF2>>h_xF2_nh3");
    tree_carbon->Draw("xF2>>h_xF2_carbon");

    // Create canvas and divide it into four panels
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(2, 2);

    // First panel: plot xF2 histograms
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    h_xF2_nh3->SetLineColor(kBlue);
    h_xF2_carbon->SetLineColor(kRed);
    h_xF2_nh3->Draw();
    h_xF2_carbon->Draw("SAME");

    // Add legend
    TLegend *leg = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg->AddEntry(h_xF2_nh3, "NH3", "l");
    leg->AddEntry(h_xF2_carbon, "Carbon", "l");
    leg->Draw();

    // Second panel: ratio of NH3 to Carbon counts
    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_ratio = new TGraphErrors();
    TH1D *h_ratio = new TH1D("h_ratio", "NH3/Carbon Ratio; xF2; Ratio", 100, -2.5, 0);
    for (int i = 1; i <= h_xF2_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xF2_nh3->GetBinContent(i);
        double c_counts = h_xF2_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio->SetPoint(i - 1, h_xF2_nh3->GetBinCenter(i), ratio);
            gr_ratio->SetPointError(i - 1, 0, error);
            h_ratio->SetBinContent(i, ratio);
            h_ratio->SetBinError(i, error);
        }
    }
    gr_ratio->SetTitle("NH3 to Carbon Ratio; xF2; Ratio");
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");

    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const = new TF1("fit_const", "[0]", -2.5, -1);
    h_ratio->Fit(fit_const, "R");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");

    // Third panel: pTpT histograms scaled by the fit constant
    c1->cd(3);
    gPad->SetLeftMargin(0.15);
    TH1D *h_pTpT_nh3 = new TH1D("h_pTpT_nh3", "pTpT Distribution; pTpT; Counts", 100, 0, 1);
    TH1D *h_pTpT_carbon = new TH1D("h_pTpT_carbon", "pTpT Distribution; pTpT; Counts", 100, 0, 1);
    tree_nh3->Draw("pTpT>>h_pTpT_nh3");
    tree_carbon->Draw("pTpT>>h_pTpT_carbon");

    double scale_factor = fit_const->GetParameter(0);
    h_pTpT_carbon->Scale(scale_factor);

    h_pTpT_nh3->SetLineColor(kBlue);
    h_pTpT_carbon->SetLineColor(kRed);
    h_pTpT_nh3->Draw();
    h_pTpT_carbon->Draw("SAME");

    // Add legend
    TLegend *leg_pTpT = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg_pTpT->AddEntry(h_pTpT_nh3, "NH3", "l");
    leg_pTpT->AddEntry(h_pTpT_carbon, "Carbon (scaled)", "l");
    leg_pTpT->Draw();

    // Fourth panel: (NH3 - Carbon) / NH3 with fit to a third-degree polynomial
    c1->cd(4);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= h_xF2_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xF2_nh3->GetBinContent(i);
        double c_counts = h_xF2_carbon->GetBinContent(i) * scale_factor;
        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;
            double error = std::sqrt(nh3_counts + std::pow(c_counts * scale_factor, 2)) / nh3_counts;
            gr_dilution->SetPoint(i - 1, h_xF2_nh3->GetBinCenter(i), dilution);
            gr_dilution->SetPointError(i - 1, 0, error);
        }
    }
    gr_dilution->SetTitle("Dilution Factor; xF2; (NH3 - Carbon) / NH3");
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");

    // Fit to a third-degree polynomial
    TF1 *fit_poly = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2 + [3]*x^3", -2.5, 1);
    gr_dilution->Fit(fit_poly, "R");
    fit_poly->SetLineColor(kRed);
    fit_poly->Draw("SAME");

    // Add fit parameters box
    TPaveText *pt = new TPaveText(0.7, 0.7, 0.9, 0.9, "brNDC");
    pt->SetBorderSize(1);
    pt->SetFillStyle(0);
    pt->AddText(Form("p0 = %.3f", fit_poly->GetParameter(0)));
    pt->AddText(Form("p1 = %.3f", fit_poly->GetParameter(1)));
    pt->AddText(Form("p2 = %.3f", fit_poly->GetParameter(2)));
    pt->AddText(Form("p3 = %.3f", fit_poly->GetParameter(3)));
    pt->Draw();

    // Save the canvas
    c1->SaveAs("dilution_factors_updated.png");

    // Clean up
    nh3->Close();
    carbon->Close();
    delete h_xF2_nh3;
    delete h_xF2_carbon;
    delete h_pTpT_nh3;
    delete h_pTpT_carbon;
    delete gr_ratio;
    delete fit_const;
    delete gr_dilution;
    delete fit_poly;
    delete stats;
    delete c1;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file>" << std::endl;
        return 1;
    }

    dilution_factors(argv[1], argv[2]);
    return 0;
}