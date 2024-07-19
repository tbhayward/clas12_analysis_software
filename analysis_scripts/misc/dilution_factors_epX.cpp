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
#include <utility> 
#include "Math/MinimizerOptions.h"

std::pair<double, double> scale_normalization(const char* nh3_file, const char* c_file) {
    // Open the ROOT files
    TFile *nh3 = TFile::Open(nh3_file);
    TFile *carbon = TFile::Open(c_file);
    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return std::make_pair(-1, -1);
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
        return std::make_pair(-1, -1);
    }

    // Create histograms for Mx and xF
    TH1D *h_Mx_nh3 = 
        new TH1D("h_Mx_nh3", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);
    TH1D *h_Mx_carbon = 
        new TH1D("h_Mx_carbon", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);
    // Fill the histograms
    tree_nh3->Draw("Mx>>h_Mx_nh3");
    tree_carbon->Draw("Mx>>h_Mx_carbon");
    //
    TH1D *h_xF_nh3 = 
        new TH1D("h_xF_nh3", "x_{F} Distribution; x_{F} (GeV); Counts", 100, -2, 1);
    TH1D *h_xF_carbon = 
        new TH1D("h_xF_carbon", "x_{F} Distribution; x_{F} (GeV); Counts", 100, -2, 1);
    // Fill the histograms
    tree_nh3->Draw("xF>>h_xF_nh3");
    tree_carbon->Draw("xF>>h_xF_carbon");

    // Create canvas and divide it into four panels
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(2, 2);

    // First panel: plot Mx histograms
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_Mx_nh3->SetLineColor(kBlue);
    h_Mx_carbon->SetLineColor(kRed);
    h_Mx_nh3->Draw();
    h_Mx_carbon->Draw("SAME");

    // Add legend
    TLegend *leg = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg->AddEntry(h_Mx_nh3, "NH_{3}", "l");
    leg->AddEntry(h_Mx_carbon, "C", "l");
    leg->Draw();

    // Remove statboxes
    h_Mx_nh3->SetStats(0);
    h_Mx_carbon->SetStats(0);

    // Second panel: ratio of NH3 to Carbon counts
    c1->cd(2);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_ratio = new TGraphErrors();
    for (int i = 1; i <= h_Mx_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_Mx_nh3->GetBinContent(i);
        double c_counts = h_Mx_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio->SetPoint(i - 1, h_Mx_nh3->GetBinCenter(i), ratio);
            gr_ratio->SetPointError(i - 1, 0, error);
        }
    }
    // Set y-axis range from 5 to 15
    gr_ratio->GetYaxis()->SetRangeUser(9, 15);
    // Set x-axis range
    gr_ratio->GetXaxis()->SetLimits(-2, 3);

    // gr_ratio->SetTitle("NH_{3} to Carbon Ratio; x_{F2}; Ratio");
    gr_ratio->SetTitle("NH_{3} to Carbon Ratio; M_{x} (GeV); Ratio");
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");

    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const = new TF1("fit_const", "[0]", -2, -0.5);
    gr_ratio->Fit(fit_const, "R");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");

    // Add the dotted-dashed line from -0.5 to 3
    TF1 *dotted_line = new TF1("dotted_line", "[0]", -0.5, 3);
    dotted_line->SetParameter(0, fit_const->GetParameter(0)); 
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

    // Third panel: plot xF histograms
    c1->cd(3);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_xF_nh3->SetLineColor(kBlue);
    h_xF_carbon->SetLineColor(kRed);
    h_xF_nh3->Draw();
    h_xF_carbon->Draw("SAME");

    // Add legend
    TLegend *leg_xF = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg_xF->AddEntry(h_xF_nh3, "NH_{3}", "l");
    leg_xF->AddEntry(h_xF_carbon, "C", "l");
    leg_xF->Draw();

    // Remove statboxes
    h_xF_nh3->SetStats(0);
    h_xF_carbon->SetStats(0);

    // Fourth panel: ratio of NH3 to Carbon counts
    c1->cd(4);
    gPad->SetLeftMargin(0.15);
    TGraphErrors *gr_ratio_xF = new TGraphErrors();
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
    // Set y-axis range from 5 to 15
    gr_ratio_xF->GetYaxis()->SetRangeUser(9, 15);
    // Set x-axis range
    gr_ratio_xF->GetXaxis()->SetLimits(-2, 1);

    gr_ratio_xF->SetTitle("NH_{3} to Carbon Ratio; x_{F} (GeV); Ratio");
    gr_ratio_xF->SetMarkerStyle(20);
    gr_ratio_xF->Draw("AP");

    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const_xF = new TF1("fit_const", "[0]", -2, -1.0);
    gr_ratio_xF->Fit(fit_const_xF, "R");
    fit_const_xF->SetLineColor(kRed);
    fit_const_xF->Draw("SAME");

    // Add the dotted-dashed line from -0.5 to 3
    TF1 *dotted_line_xF = new TF1("dotted_line_xF", "[0]", -1, 1);
    dotted_line_xF->SetParameter(0, fit_const_xF->GetParameter(0));
    dotted_line_xF->SetLineColor(kRed);
    dotted_line_xF->SetLineStyle(7); // Set line style to dotted-dashed
    dotted_line_xF->Draw("SAME");

    // Add fit constant value and uncertainty
    double fit_value_xF = fit_const_xF->GetParameter(0);
    double fit_error_xF = fit_const_xF->GetParError(0);

    // Retrieve chi2 and NDF
    double chi2_xF = fit_const_xF->GetChisquare();
    int ndf_xF = fit_const_xF->GetNDF();
    double chi2_ndf_xF = chi2_xF / ndf;

    TLatex latex_xF;
    latex_xF.SetNDC();
    latex_xF.SetTextSize(0.04);
    latex_xF.DrawLatex(0.20, 0.85, 
        Form("Fit Const, s = %.3f #pm %.3f", fit_value_xF, fit_error_xF));
    latex_xF.DrawLatex(0.20, 0.80, 
        Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_xF, ndf_xF, chi2_ndf_xF));

    // Save the canvas
    c1->SaveAs("output/scale_constant.pdf");
    // Clean up
    nh3->Close();
    carbon->Close();
    delete c1;

    // Return the fit value and error as a pair
    return std::make_pair(fit_value, fit_error);
}

double one_dimensional(const char* nh3_file, const char* c_file, 
    std::pair<double, double> fit_constant) {
    double scale_factor = fit_constant.first;
    double scale_error = fit_constant.second;

    // Open the ROOT files
    TFile *nh3 = TFile::Open(nh3_file);
    TFile *carbon = TFile::Open(c_file);
    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return 0;
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
        return 0;
    }

    // Create canvas and divide it into four panels
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1600, 1200);
    c1->Divide(3, 2);

    // Third panel: pT histograms scaled by the fit constant with propagated errors
    c1->cd(1);
    gPad->SetLeftMargin(0.15);
    TH1D *h_pT_nh3 = 
        new TH1D("h_pT_nh3", "P_{T} Distribution; P_{T} (GeV); Counts", 100, 0, 1.0);
    TH1D *h_pT_carbon = 
        new TH1D("h_pT_carbon", "P_{T} Distribution; P_{T} (GeV); Counts", 100, 0, 1.0);
    tree_nh3->Draw("pT>>h_pT_nh3");
    tree_carbon->Draw("pT>>h_pT_carbon");
    TH1D *h_pT_carbon_scaled = (TH1D*)h_pT_carbon->Clone("h_pT_carbon_scaled");
    h_pT_carbon_scaled->SetTitle("P_{T} Distribution; P_{T} (GeV); Counts (Scaled)");

    for (int i = 1; i <= h_pT_carbon->GetNbinsX(); ++i) {
        double bin_content = h_pT_carbon->GetBinContent(i);
        double bin_error = h_pT_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_pT_carbon_scaled->SetBinContent(i, new_content);
        h_pT_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= h_pT_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_pT_nh3->GetBinContent(i);
        double nh3_error = h_pT_nh3->GetBinError(i);
        double c_counts = h_pT_carbon_scaled->GetBinContent(i);
        double c_error = h_pT_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution->SetPoint(i - 1, h_pT_nh3->GetBinCenter(i), dilution);
            gr_dilution->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution->SetTitle("; P_{T} (GeV); D_f = (NH3 - s*C) / NH3");
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");
    gr_dilution->GetXaxis()->SetRangeUser(0, 1);
    gr_dilution->GetYaxis()->SetRangeUser(0.05, 0.15);

    // Fit to a third-degree polynomial
    TF1 *fit_poly = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2 + [3]*x^3", 0, 1.0);
    gr_dilution->Fit(fit_poly, "R");
    fit_poly->SetLineColor(kRed);
    fit_poly->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0 = fit_poly->GetParameter(0);
    double p0_err = fit_poly->GetParError(0);
    double p1 = fit_poly->GetParameter(1);
    double p1_err = fit_poly->GetParError(1);
    double p2 = fit_poly->GetParameter(2);
    double p2_err = fit_poly->GetParError(2);
    double p3 = fit_poly->GetParameter(3);
    double p3_err = fit_poly->GetParError(3);

    // Retrieve chi2 and NDF
    double chi2 = fit_poly->GetChisquare();
    int ndf = fit_poly->GetNDF();
    double chi2_ndf = chi2 / ndf;

    // Add fit parameters box
    TPaveText *pt = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    pt->SetBorderSize(1);
    pt->SetFillStyle(1001); // Solid fill style
    pt->SetFillColor(kWhite); // White background
    pt->AddText(Form("p0 = %.3f +/- %.3f", p0, p0_err));
    pt->AddText(Form("p1 = %.3f +/- %.3f", p1, p1_err));
    pt->AddText(Form("p2 = %.3f +/- %.3f", p2, p2_err));
    pt->AddText(Form("p3 = %.3f +/- %.3f", p3, p3_err));
    pt->Draw();

    // Add chi2/ndf in the top left
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2_ndf));

    std::cout << "if (prefix == \"PT\") { return " << p0 << 
        "+" << p1 << "*currentVariable" << p2 << "std::pow(currentVariable,2)" <<
        p3 << "*std::pow(currentVariable,3);" << std::endl;

    // Save the canvas
    c1->SaveAs("output/one_dimensional.pdf");
    // Clean up
    nh3->Close();
    carbon->Close();

    return 0;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file>" << std::endl;
        return 1;
    }

    std::pair<double, double> fit_constant = scale_normalization(argv[1], argv[2]);
    one_dimensional(argv[1], argv[2], fit_constant);
}