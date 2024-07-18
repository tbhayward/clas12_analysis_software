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
    TFile* nh3 = TFile::Open(nh3_file);
    TFile* carbon = TFile::Open(c_file);
    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    // Get the PhysicsEvents trees
    TTree* tree_nh3;
    TTree* tree_carbon;
    nh3->GetObject("PhysicsEvents", tree_nh3);
    carbon->GetObject("PhysicsEvents", tree_carbon);
    if (!tree_nh3 || !tree_carbon) {
        std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
        nh3->Close();
        carbon->Close();
        return;
    }

    // Create canvas and divide it into four panels
    TCanvas* c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(2, 2);

    //~~~~~~~~~~~~~~~~~~~~//
    // First set of plots: Mx
    //~~~~~~~~~~~~~~~~~~~~//
    // Create histograms for Mx
    TH1D* h_Mx_nh3 = new TH1D("h_Mx_nh3", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);
    TH1D* h_Mx_carbon = new TH1D("h_Mx_carbon", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);

    // Fill the histograms
    tree_nh3->Draw("Mx>>h_Mx_nh3");
    tree_carbon->Draw("Mx>>h_Mx_carbon");

    // First panel: plot Mx histograms
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_Mx_nh3->SetLineColor(kBlue);
    h_Mx_carbon->SetLineColor(kRed);
    h_Mx_nh3->Draw();
    h_Mx_carbon->Draw("SAME");

    // Add legend
    TLegend* leg_Mx = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg_Mx->AddEntry(h_Mx_nh3, "NH_{3}", "l");
    leg_Mx->AddEntry(h_Mx_carbon, "C", "l");
    leg_Mx->Draw();
    // Remove statboxes
    h_Mx_nh3->SetStats(0);
    h_Mx_carbon->SetStats(0);

    // Second panel: ratio of NH3 to Carbon counts for Mx
    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    TGraphErrors* gr_ratio_Mx = new TGraphErrors();
    for (int i = 1; i <= h_Mx_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_Mx_nh3->GetBinContent(i);
        double c_counts = h_Mx_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio_Mx->SetPoint(i - 1, h_Mx_nh3->GetBinCenter(i), ratio);
            gr_ratio_Mx->SetPointError(i - 1, 0, error);
        }
    }
    gr_ratio_Mx->GetYaxis()->SetRangeUser(9, 15);
    gr_ratio_Mx->SetTitle("NH_{3} to Carbon Ratio; M_{x} (GeV); Ratio");
    gr_ratio_Mx->SetMarkerStyle(20);
    gr_ratio_Mx->Draw("AP");

    // Fit the data from -2 to -0.5 to a constant for Mx
    TF1* fit_const_Mx = new TF1("fit_const_Mx", "[0]", -2, -0.5);
    gr_ratio_Mx->Fit(fit_const_Mx, "R");
    fit_const_Mx->SetLineColor(kRed);
    fit_const_Mx->Draw("SAME");
    // Add the dotted-dashed line from -0.5 to 3 for Mx
    TF1* dotted_line_Mx = new TF1("dotted_line_Mx", "[0]", -0.5, 3);
    dotted_line_Mx->SetParameter(0, fit_const_Mx->GetParameter(0)); // Use the same constant value as the fit
    dotted_line_Mx->SetLineColor(kRed);
    dotted_line_Mx->SetLineStyle(7); // Set line style to dotted-dashed
    dotted_line_Mx->Draw("SAME");

    // Add fit constant value and uncertainty for Mx
    double fit_value_Mx = fit_const_Mx->GetParameter(0);
    double fit_error_Mx = fit_const_Mx->GetParError(0);
    // Retrieve chi2 and NDF for Mx
    double chi2_Mx = fit_const_Mx->GetChisquare();
    int ndf_Mx = fit_const_Mx->GetNDF();
    double chi2_ndf_Mx = chi2_Mx / ndf_Mx;
    TLatex latex_Mx;
    latex_Mx.SetNDC();
    latex_Mx.SetTextSize(0.04);
    latex_Mx.DrawLatex(0.20, 0.85, Form("Fit Const, s = %.3f #pm %.3f", fit_value_Mx, fit_error_Mx));
    latex_Mx.DrawLatex(0.20, 0.80, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_Mx, ndf_Mx, chi2_ndf_Mx));

    //~~~~~~~~~~~~~~~~~~~~//
    // Second set of plots: xF
    //~~~~~~~~~~~~~~~~~~~~//
    // Create histograms for xF
    TH1D* h_xF_nh3 = new TH1D("h_xF_nh3", "x_{F} Distribution; x_{F}; Counts", 100, -2.5, 1);
    TH1D* h_xF_carbon = new TH1D("h_xF_carbon", "x_{F} Distribution; x_{F}; Counts", 100, -2.5, 1);

    // Fill the histograms
    tree_nh3->Draw("xF>>h_xF_nh3");
    tree_carbon->Draw("xF>>h_xF_carbon");

    // Third panel: plot xF histograms
    c1->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_xF_nh3->SetLineColor(kBlue);
    h_xF_carbon->SetLineColor(kRed);
    h_xF_nh3->Draw();
    h_xF_carbon->Draw("SAME");

    // Add legend
    TLegend* leg_xF = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg_xF->AddEntry(h_xF_nh3, "NH_{3}", "l");
    leg_xF->AddEntry(h_xF_carbon, "C", "l");
    leg_xF->Draw();
    // Remove statboxes
    h_xF_nh3->SetStats(0);
    h_xF_carbon->SetStats(0);

    // Fourth panel: ratio of NH3 to Carbon counts for xF
    c1->cd(4);
    gPad->SetLeftMargin(0.15);
    TGraphErrors* gr_ratio_xF = new TGraphErrors();
    for (int i = 1; i <= h_xF_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xF_nh3->GetBinContent(i);
        double c_counts = h_xF_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio_xF->SetPoint(i - 1, h_xF_nh3->GetBinCenter(i), ratio);
            gr_ratio_xF->SetPointError(i - 1, 0, error);
        }
    }
    gr_ratio_xF->GetYaxis()->SetRangeUser(9, 15);
    gr_ratio_xF->SetTitle("NH_{3} to Carbon Ratio; x_{F}; Ratio");
    gr_ratio_xF->SetMarkerStyle(20);
    gr_ratio_xF->Draw("AP");

    // Fit the data from -2 to -1.25 to a constant for xF
    TF1* fit_const_xF = new TF1("fit_const_xF", "[0]", -2, -1.25);
    gr_ratio_xF->Fit(fit_const_xF, "R");
    fit_const_xF->SetLineColor(kRed);
    fit_const_xF->Draw("SAME");
    // Add the dotted-dashed line from -1.25 to 1 for xF
    TF1* dotted_line_xF = new TF1("dotted_line_xF", "[0]", -1.25, 1);
    dotted_line_xF->SetParameter(0, fit_const_xF->GetParameter(0)); // Use the same constant value as the fit
    dotted_line_xF->SetLineColor(kRed);
    dotted_line_xF->SetLineStyle(7); // Set line style to dotted-dashed
    dotted_line_xF->Draw("SAME");

    // Add fit constant value and uncertainty for xF
    double fit_value_xF = fit_const_xF->GetParameter(0);
    double fit_error_xF = fit_const_xF->GetParError(0);
    // Retrieve chi2 and NDF for xF
    double chi2_xF = fit_const_xF->GetChisquare();
    int ndf_xF = fit_const_xF->GetNDF();
    double chi2_ndf_xF = chi2_xF / ndf_xF;
    TLatex latex_xF;
    latex_xF.SetNDC();
    latex_xF.SetTextSize(0.04);
    latex_xF.DrawLatex(0.20, 0.85, Form("Fit Const, s = %.3f #pm %.3f", fit_value_xF, fit_error_xF));
    latex_xF.DrawLatex(0.20, 0.80, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_xF, ndf_xF, chi2_ndf_xF));

    // Save the canvas
    c1->SaveAs("dilution_factors.pdf");

    // Clean up
    nh3->Close();
    carbon->Close();
    }

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " " << std::endl;
        return 1;
    }
    dilution_factors_epX(argv[1], argv[2]);
    return 0;
}