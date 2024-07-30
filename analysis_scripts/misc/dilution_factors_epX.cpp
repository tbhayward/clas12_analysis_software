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
#include <TError.h> // Include the header for error handling

std::string reformatRange(const std::string &range) {
    std::string formatted;
    std::stringstream ss(range);
    std::string token;
    std::string lower_bound, upper_bound;
    std::string variable;

    while (std::getline(ss, token, ' ')) {
        if (token.find(">") != std::string::npos && token.find("<") == std::string::npos) {
            lower_bound = token.substr(0, token.find(">"));
            variable = token.substr(token.find(">") + 1);
        } else if (token.find("<") != std::string::npos && token.find(">") == std::string::npos) {
            upper_bound = token.substr(token.find("<") + 1);
        }
    }

    formatted = lower_bound + " < " + variable + " < " + upper_bound;
    return formatted;
}

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
        new TH1D("h_Mx_nh3", "M_{x} Distribution; M_{x} (GeV); Counts", 50, -2, 3);
    TH1D *h_Mx_carbon = 
        new TH1D("h_Mx_carbon", "M_{x} Distribution; M_{x} (GeV); Counts", 50, -2, 3);
    // Fill the histograms
    tree_nh3->Draw("Mx>>h_Mx_nh3");
    tree_carbon->Draw("Mx>>h_Mx_carbon");
    //
    TH1D *h_xF_nh3 = 
        new TH1D("h_xF_nh3", "x_{F} Distribution; x_{F} (GeV); Counts", 50, -2, 1);
    TH1D *h_xF_carbon = 
        new TH1D("h_xF_carbon", "x_{F} Distribution; x_{F} (GeV); Counts", 50, -2, 1);
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
    // gr_ratio->GetYaxis()->SetRangeUser(1, 2);
    // Set x-axis range
    gr_ratio->GetXaxis()->SetLimits(-2, 3);

    // gr_ratio->SetTitle("NH_{3} to Carbon Ratio; x_{F2}; Ratio");
    gr_ratio->SetTitle("NH_{3} to Carbon Ratio; M_{x} (GeV); Ratio");
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");

    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const = new TF1("fit_const", "[0]", -2, -0.0);
    gr_ratio->Fit(fit_const, "RQ");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");

    // Add the dotted-dashed line from -0.5 to 3
    TF1 *dotted_line = new TF1("dotted_line", "[0]", -0.0, 3);
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

    gr_ratio_xF->SetTitle("NH_{3} to Carbon Ratio; x_{F}; Ratio");
    gr_ratio_xF->SetMarkerStyle(20);
    gr_ratio_xF->Draw("AP");

    // Fit the data from -2.5 to -1 to a constant
    TF1 *fit_const_xF = new TF1("fit_const", "[0]", -2, -1.0);
    gr_ratio_xF->Fit(fit_const_xF, "RQ");
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

    //pT histograms scaled by the fit constant with propagated errors
    c1->cd(4);
    gPad->SetLeftMargin(0.15);
    TH1D *h_pT_nh3 = 
        new TH1D("h_pT_nh3", "P_{T} Distribution; P_{T} (GeV); Counts", 50, 0, 1.0);
    TH1D *h_pT_carbon = 
        new TH1D("h_pT_carbon", "P_{T} Distribution; P_{T} (GeV); Counts", 50, 0, 1.0);
    tree_nh3->Draw("pT>>h_pT_nh3","Mx > 1.4");
    tree_carbon->Draw("pT>>h_pT_carbon","Mx > 1.4");
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
    gr_dilution->SetTitle("; P_{T} (GeV); D_{f} = (NH3 - s*C) / NH3");
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");
    gr_dilution->GetXaxis()->SetRangeUser(0, 1);
    gr_dilution->GetYaxis()->SetRangeUser(0.00, 0.20);

    // Fit to a third-degree polynomial
    TF1 *fit_poly = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2", 0, 1.0);
    gr_dilution->Fit(fit_poly, "RQ");
    fit_poly->SetLineColor(kRed);
    fit_poly->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0 = fit_poly->GetParameter(0);
    double p0_err = fit_poly->GetParError(0);
    double p1 = fit_poly->GetParameter(1);
    double p1_err = fit_poly->GetParError(1);
    double p2 = fit_poly->GetParameter(2);
    double p2_err = fit_poly->GetParError(2);

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
    pt->Draw();

    // Add chi2/ndf in the top left
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2_ndf));

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    c1->cd(1);
    gPad->SetLeftMargin(0.15);

    // Third panel: x histograms scaled by the fit constant with propagated errors
    TH1D *h_x_nh3 = 
        new TH1D("h_x_nh3", "x_{B} Distribution; x_{B} (GeV); Counts", 50, 0.06, 0.6);
    TH1D *h_x_carbon = 
        new TH1D("h_x_carbon", "x_{B} Distribution; x_{B} (GeV); Counts", 50, 0.06, 0.6);
    tree_nh3->Draw("x>>h_x_nh3","Mx > 1.4");
    tree_carbon->Draw("x>>h_x_carbon","Mx > 1.4");
    TH1D *h_x_carbon_scaled = (TH1D*)h_x_carbon->Clone("h_x_carbon_scaled");
    h_x_carbon_scaled->SetTitle("x_{B} Distribution; x_{B}; Counts (Scaled)");

    for (int i = 1; i <= h_x_carbon->GetNbinsX(); ++i) {
        double bin_content = h_x_carbon->GetBinContent(i);
        double bin_error = h_x_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_x_carbon_scaled->SetBinContent(i, new_content);
        h_x_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution_x = new TGraphErrors();
    for (int i = 1; i <= h_x_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_x_nh3->GetBinContent(i);
        double nh3_error = h_x_nh3->GetBinError(i);
        double c_counts = h_x_carbon_scaled->GetBinContent(i);
        double c_error = h_x_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution_x->SetPoint(i - 1, h_x_nh3->GetBinCenter(i), dilution);
            gr_dilution_x->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution_x->SetTitle("; x_{B}; D_{f} = (NH3 - s*C) / NH3");
    gr_dilution_x->SetMarkerStyle(20);
    gr_dilution_x->Draw("AP");
    gr_dilution_x->GetXaxis()->SetRangeUser(0, 0.6);
    gr_dilution_x->GetYaxis()->SetRangeUser(0.00, 0.20);

    // Fit to a third-degree polynomial
    TF1 *fit_poly_x = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2", 0.06, 0.6);
    gr_dilution_x->Fit(fit_poly_x, "RQ");
    fit_poly_x->SetLineColor(kRed);
    fit_poly_x->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0_x = fit_poly_x->GetParameter(0);
    double p0_err_x = fit_poly_x->GetParError(0);
    double p1_x = fit_poly_x->GetParameter(1);
    double p1_err_x = fit_poly_x->GetParError(1);
    double p2_x = fit_poly_x->GetParameter(2);
    double p2_err_x = fit_poly_x->GetParError(2);

    // Retrieve chi2 and NDF
    double chi2_x = fit_poly_x->GetChisquare();
    int ndf_x = fit_poly_x->GetNDF();
    double chi2_ndf_x = chi2_x / ndf;

    // Add fit parameters box
    TPaveText *x = new TPaveText(0.15, 0.7, 0.55, 0.9, "brNDC");
    x->SetBorderSize(1);
    x->SetFillStyle(1001); // Solid fill style
    x->SetFillColor(kWhite); // White background
    x->AddText(Form("p0 = %.3f +/- %.3f", p0_x, p0_err_x));
    x->AddText(Form("p1 = %.3f +/- %.3f", p1_x, p1_err_x));
    x->AddText(Form("p2 = %.3f +/- %.3f", p2_x, p2_err_x));
    x->Draw();

    // Add chi2/ndf in the top left
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, 
        Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_x, ndf_x, chi2_ndf_x));


    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    c1->cd(2);
    gPad->SetLeftMargin(0.15);

    // Third panel: z histograms scaled by the fit constant with propagated errors
    TH1D *h_z_nh3 = 
        new TH1D("h_z_nh3", "z Distribution; z (GeV); Counts", 50, 0.1, 0.75);
    TH1D *h_z_carbon = 
        new TH1D("h_z_carbon", "z Distribution; z (GeV); Counts", 50, 0.1, 0.75);
    tree_nh3->Draw("z>>h_z_nh3","Mx > 1.4");
    tree_carbon->Draw("z>>h_z_carbon","Mx > 1.4");
    TH1D *h_z_carbon_scaled = (TH1D*)h_z_carbon->Clone("h_z_carbon_scaled");
    h_z_carbon_scaled->SetTitle("z Distribution; z; Counts (Scaled)");

    for (int i = 1; i <= h_z_carbon->GetNbinsX(); ++i) {
        double bin_content = h_z_carbon->GetBinContent(i);
        double bin_error = h_z_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_z_carbon_scaled->SetBinContent(i, new_content);
        h_z_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution_z = new TGraphErrors();
    for (int i = 1; i <= h_z_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_z_nh3->GetBinContent(i);
        double nh3_error = h_z_nh3->GetBinError(i);
        double c_counts = h_z_carbon_scaled->GetBinContent(i);
        double c_error = h_z_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution_z->SetPoint(i - 1, h_z_nh3->GetBinCenter(i), dilution);
            gr_dilution_z->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution_z->SetTitle("; z; D_{f} = (NH3 - s*C) / NH3");
    gr_dilution_z->SetMarkerStyle(20);
    gr_dilution_z->Draw("AP");
    gr_dilution_z->GetXaxis()->SetRangeUser(0.1, 0.75);
    gr_dilution_z->GetYaxis()->SetRangeUser(0.00, 0.20);

    // Fit to a third-degree polynomial
    TF1 *fit_poly_z = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2 + [3]*x^3", 0.1, 0.75);
    gr_dilution_z->Fit(fit_poly_z, "RQ");
    fit_poly_z->SetLineColor(kRed);
    fit_poly_z->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0_z = fit_poly_z->GetParameter(0);
    double p0_err_z = fit_poly_z->GetParError(0);
    double p1_z = fit_poly_z->GetParameter(1);
    double p1_err_z = fit_poly_z->GetParError(1);
    double p2_z = fit_poly_z->GetParameter(2);
    double p2_err_z = fit_poly_z->GetParError(2);
    double p3_z = fit_poly_z->GetParameter(3);
    double p3_err_z = fit_poly_z->GetParError(3);

    // Retrieve chi2 and NDF
    double chi2_z = fit_poly_z->GetChisquare();
    int ndf_z = fit_poly_z->GetNDF();
    double chi2_ndf_z = chi2_z / ndf;

    // Add fit parameters box
    TPaveText *z = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    z->SetBorderSize(1);
    z->SetFillStyle(1001); // Solid fill style
    z->SetFillColor(kWhite); // White background
    z->AddText(Form("p0 = %.3f +/- %.3f", p0_z, p0_err_z));
    z->AddText(Form("p1 = %.3f +/- %.3f", p1_z, p1_err_z));
    z->AddText(Form("p2 = %.3f +/- %.3f", p2_z, p2_err_z));
    z->AddText(Form("p3 = %.3f +/- %.3f", p3_z, p3_err_z));
    z->Draw();

    // Add chi2/ndf in the top left
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.30, 0.15, 
        Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_z, ndf_z, chi2_ndf_z));

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    c1->cd(5);
    gPad->SetLeftMargin(0.15);

    // Third panel: x_{F} histograms scaled by the fit constant with propagated errors
    TH1D *h_xF_nh3 = 
        new TH1D("h_xF_nh3", "x_{F} Distribution; x_{F} (GeV); Counts", 50, -0.8, 0.5);
    TH1D *h_xF_carbon = 
        new TH1D("h_xF_carbon", "x_{F} Distribution; x_{F} (GeV); Counts", 50, -0.8, 0.5);
    tree_nh3->Draw("xF>>h_xF_nh3","Mx > 1.4");
    tree_carbon->Draw("xF>>h_xF_carbon","Mx > 1.4");
    TH1D *h_xF_carbon_scaled = (TH1D*)h_xF_carbon->Clone("h_xF_carbon_scaled");
    h_xF_carbon_scaled->SetTitle("x_{F} Distribution; x_{F}; Counts (Scaled)");

    for (int i = 1; i <= h_xF_carbon->GetNbinsX(); ++i) {
        double bin_content = h_xF_carbon->GetBinContent(i);
        double bin_error = h_xF_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_xF_carbon_scaled->SetBinContent(i, new_content);
        h_xF_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution_xF = new TGraphErrors();
    for (int i = 1; i <= h_xF_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_xF_nh3->GetBinContent(i);
        double nh3_error = h_xF_nh3->GetBinError(i);
        double c_counts = h_xF_carbon_scaled->GetBinContent(i);
        double c_error = h_xF_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution_xF->SetPoint(i - 1, h_xF_nh3->GetBinCenter(i), dilution);
            gr_dilution_xF->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution_xF->SetTitle("; x_{F}; D_{f} = (NH3 - s*C) / NH3");
    gr_dilution_xF->SetMarkerStyle(20);
    gr_dilution_xF->Draw("AP");
    gr_dilution_xF->GetXaxis()->SetRangeUser(-0.8, 0.5);
    gr_dilution_xF->GetYaxis()->SetRangeUser(0.00, 0.20);

    // Fit to a third-degree polynomial
    TF1 *fit_poly_xF = new TF1("fit_poly", "[0] + [1]*x + [2]*x^2 + [3]*x^3", -0.8, 0.5);
    gr_dilution_xF->Fit(fit_poly_xF, "RQ");
    fit_poly_xF->SetLineColor(kRed);
    fit_poly_xF->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0_xF = fit_poly_xF->GetParameter(0);
    double p0_err_xF = fit_poly_xF->GetParError(0);
    double p1_xF = fit_poly_xF->GetParameter(1);
    double p1_err_xF = fit_poly_xF->GetParError(1);
    double p2_xF = fit_poly_xF->GetParameter(2);
    double p2_err_xF = fit_poly_xF->GetParError(2);
    double p3_xF = fit_poly_xF->GetParameter(3);
    double p3_err_xF = fit_poly_xF->GetParError(3);

    // Retrieve chi2 and NDF
    double chi2_xF = fit_poly_xF->GetChisquare();
    int ndf_xF = fit_poly_xF->GetNDF();
    double chi2_ndf_xF = chi2_xF / ndf;

    // Add fit parameters box
    TPaveText *xF = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    xF->SetBorderSize(1);
    xF->SetFillStyle(1001); // Solid fill style
    xF->SetFillColor(kWhite); // White background
    xF->AddText(Form("p0 = %.3f +/- %.3f", p0_xF, p0_err_xF));
    xF->AddText(Form("p1 = %.3f +/- %.3f", p1_xF, p1_err_xF));
    xF->AddText(Form("p2 = %.3f +/- %.3f", p2_xF, p2_err_xF));
    xF->AddText(Form("p3 = %.3f +/- %.3f", p3_xF, p3_err_xF));
    xF->Draw();

    // Add chi2/ndf in the top left
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.30, 0.15, 
        Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_xF, ndf_xF, chi2_ndf_xF));

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    c1->cd(3);
    gPad->SetLeftMargin(0.15);

    // Third panel: zeta histograms scaled by the fit constant with propagated errors
    TH1D *h_zeta_nh3 = 
        new TH1D("h_zeta_nh3", "#zeta Distribution; #zeta (GeV); Counts", 50, 0.3, 0.8);
    TH1D *h_zeta_carbon = 
        new TH1D("h_zeta_carbon", "z Distribution; z (GeV); Counts", 50, 0.3, 0.8);
    tree_nh3->Draw("zeta>>h_zeta_nh3","Mx > 1.4");
    tree_carbon->Draw("zeta>>h_zeta_carbon","Mx > 1.4");
    TH1D *h_zeta_carbon_scaled = (TH1D*)h_zeta_carbon->Clone("h_zeta_carbon_scaled");
    h_zeta_carbon_scaled->SetTitle("#zeta Distribution; #zeta; Counts (Scaled)");

    for (int i = 1; i <= h_zeta_carbon->GetNbinsX(); ++i) {
        double bin_content = h_zeta_carbon->GetBinContent(i);
        double bin_error = h_zeta_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_zeta_carbon_scaled->SetBinContent(i, new_content);
        h_zeta_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution_zeta = new TGraphErrors();
    for (int i = 1; i <= h_zeta_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_zeta_nh3->GetBinContent(i);
        double nh3_error = h_zeta_nh3->GetBinError(i);
        double c_counts = h_zeta_carbon_scaled->GetBinContent(i);
        double c_error = h_zeta_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution_zeta->SetPoint(i - 1, h_zeta_nh3->GetBinCenter(i), dilution);
            gr_dilution_zeta->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution_zeta->SetTitle("; #zeta; D_{f} = (NH3 - s*C) / NH3");
    gr_dilution_zeta->SetMarkerStyle(20);
    gr_dilution_zeta->Draw("AP");
    gr_dilution_zeta->GetXaxis()->SetRangeUser(0.3, 0.8);
    gr_dilution_zeta->GetYaxis()->SetRangeUser(0.00, 0.20);

    // Fit to a third-degree polynomial
    TF1 *fit_poly_zeta = new TF1("fit_poly", "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0.3, 0.8);
    gr_dilution_zeta->Fit(fit_poly_zeta, "RQ");
    fit_poly_zeta->SetLineColor(kRed);
    fit_poly_zeta->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0_zeta = fit_poly_zeta->GetParameter(0);
    double p0_err_zeta = fit_poly_zeta->GetParError(0);
    double p1_zeta = fit_poly_zeta->GetParameter(1);
    double p1_err_zeta = fit_poly_zeta->GetParError(1);
    double p2_zeta = fit_poly_zeta->GetParameter(2);
    double p2_err_zeta = fit_poly_zeta->GetParError(2);
    double p3_zeta = fit_poly_zeta->GetParameter(3);
    double p3_err_zeta = fit_poly_zeta->GetParError(3);
    double p4_zeta = fit_poly_zeta->GetParameter(4);
    double p4_err_zeta = fit_poly_zeta->GetParError(4);

    // Retrieve chi2 and NDF
    double chi2_zeta = fit_poly_zeta->GetChisquare();
    int ndf_zeta = fit_poly_zeta->GetNDF();
    double chi2_ndf_zeta = chi2_zeta / ndf;

    // Add fit parameters box
    TPaveText *zeta = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    zeta->SetBorderSize(1);
    zeta->SetFillStyle(1001); // Solid fill style
    zeta->SetFillColor(kWhite); // White background
    zeta->AddText(Form("p0 = %.3f +/- %.3f", p0_zeta, p0_err_zeta));
    zeta->AddText(Form("p1 = %.3f +/- %.3f", p1_zeta, p1_err_zeta));
    zeta->AddText(Form("p2 = %.3f +/- %.3f", p2_zeta, p2_err_zeta));
    zeta->AddText(Form("p3 = %.3f +/- %.3f", p3_zeta, p3_err_zeta));
    zeta->AddText(Form("p4 = %.3f +/- %.3f", p4_zeta, p4_err_zeta));
    zeta->Draw();

    // Add chi2/ndf in the top left
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, 
        Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_zeta, ndf_zeta, chi2_ndf_zeta));

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    c1->cd(6);
    gPad->SetLeftMargin(0.15);

    // Third panel: x histograms scaled by the fit constant with propagated errors
    TH1D *h_Mx_nh3 = 
        new TH1D("h_Mx_nh3", "M_x Distribution; M_x (GeV); Counts", 50, 0.0, 3);
    TH1D *h_Mx_carbon = 
        new TH1D("h_Mx_carbon", "M_x Distribution; M_x (GeV); Counts", 50, 0.0, 3);
    tree_nh3->Draw("Mx>>h_Mx_nh3");
    tree_carbon->Draw("Mx>>h_Mx_carbon");
    TH1D *h_Mx_carbon_scaled = (TH1D*)h_Mx_carbon->Clone("h_Mx_carbon_scaled");
    h_Mx_carbon_scaled->SetTitle("M_{x} (GeV) Distribution; M_{x} (GeV); Counts (Scaled)");

    for (int i = 1; i <= h_Mx_carbon->GetNbinsX(); ++i) {
        double bin_content = h_Mx_carbon->GetBinContent(i);
        double bin_error = h_Mx_carbon->GetBinError(i);

        double new_content = bin_content * scale_factor;
        double new_error = 
            new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (scale_error / scale_factor) * (scale_error / scale_factor));
        h_Mx_carbon_scaled->SetBinContent(i, new_content);
        h_Mx_carbon_scaled->SetBinError(i, new_error);
    }

    TGraphErrors *gr_dilution_Mx = new TGraphErrors();
    for (int i = 1; i <= h_Mx_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_Mx_nh3->GetBinContent(i);
        double nh3_error = h_Mx_nh3->GetBinError(i);
        double c_counts = h_Mx_carbon_scaled->GetBinContent(i);
        double c_error = h_Mx_carbon_scaled->GetBinError(i);

        if (nh3_counts > 0) {
            double dilution = (nh3_counts - c_counts) / nh3_counts;

            // Propagate the error
            double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));

            gr_dilution_Mx->SetPoint(i - 1, h_Mx_nh3->GetBinCenter(i), dilution);
            gr_dilution_Mx->SetPointError(i - 1, 0, dilution_error);
        }
    }
    gr_dilution_Mx->SetTitle("; M_{x} (GeV); D_{f} = (NH3 - s*C) / NH3");
    gr_dilution_Mx->SetMarkerStyle(20);
    gr_dilution_Mx->Draw("AP");
    gr_dilution_Mx->GetXaxis()->SetRangeUser(0.0, 3.0);
    gr_dilution_Mx->GetYaxis()->SetRangeUser(0.00, 0.20);

    // Fit to a third-degree polynomial
    // TF1 *fit_poly_Mx = new TF1("fit_poly", 
    //     "[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4", 0.0, 3.0);
    TF1 *fit_poly_Mx = new TF1("fit_gauss_sum", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[4])/[5])^2) + [6]*exp(-0.5*((x-[7])/[8])^2)", 0.0, 3.0);
    fit_poly_Mx->SetParameters(0.037, 0.770, 0.056, 0.026, 1.15, 0.29, 0.108, 1.689, 1.813);
    gr_dilution_Mx->Fit(fit_poly_Mx, "RQ");
    fit_poly_Mx->SetLineColor(kRed);
    fit_poly_Mx->Draw("SAME");

    // Retrieve fit parameters and their errors
    double p0_Mx = fit_poly_Mx->GetParameter(0);
    double p0_err_Mx = fit_poly_Mx->GetParError(0);
    double p1_Mx = fit_poly_Mx->GetParameter(1);
    double p1_err_Mx = fit_poly_Mx->GetParError(1);
    double p2_Mx = fit_poly_Mx->GetParameter(2);
    double p2_err_Mx = fit_poly_Mx->GetParError(2);
    double p3_Mx = fit_poly_Mx->GetParameter(3);
    double p3_err_Mx = fit_poly_Mx->GetParError(3);
    double p4_Mx = fit_poly_Mx->GetParameter(4);
    double p4_err_Mx = fit_poly_Mx->GetParError(4);
    double p5_Mx = fit_poly_Mx->GetParameter(5);
    double p5_err_Mx = fit_poly_Mx->GetParError(5);
    double p6_Mx = fit_poly_Mx->GetParameter(6);
    double p6_err_Mx = fit_poly_Mx->GetParError(6);
    double p7_Mx = fit_poly_Mx->GetParameter(7);
    double p7_err_Mx = fit_poly_Mx->GetParError(7);
    double p8_Mx = fit_poly_Mx->GetParameter(8);
    double p8_err_Mx = fit_poly_Mx->GetParError(8);

    // Retrieve chi2 and NDF
    double chi2_Mx = fit_poly_Mx->GetChisquare();
    int ndf_Mx = fit_poly_Mx->GetNDF();
    double chi2_ndf_Mx = chi2_Mx / ndf;

    // Add fit parameters box
    TPaveText *Mx = new TPaveText(0.5, 0.6, 0.9, 0.9, "brNDC");
    Mx->SetBorderSize(1);
    Mx->SetFillStyle(1001); // Solid fill style
    Mx->SetFillColor(kWhite); // White background
    Mx->AddText(Form("amp_{1} = %.3f +/- %.3f", p0_Mx, p0_err_Mx));
    Mx->AddText(Form("#mu_{1} = %.3f +/- %.3f", p1_Mx, p1_err_Mx));
    Mx->AddText(Form("#sigma_{1} = %.3f +/- %.3f", p2_Mx, p2_err_Mx));
    Mx->AddText(Form("amp_{2} = %.3f +/- %.3f", p3_Mx, p3_err_Mx));
    Mx->AddText(Form("#mu_{2} = %.3f +/- %.3f", p4_Mx, p4_err_Mx));
    Mx->AddText(Form("#sigma_{2} = %.3f +/- %.3f", p5_Mx, p5_err_Mx));
    Mx->AddText(Form("amp_{3} = %.3f +/- %.3f", p6_Mx, p6_err_Mx));
    Mx->AddText(Form("#mu_{3} = %.3f +/- %.3f", p7_Mx, p7_err_Mx));
    Mx->AddText(Form("#sigma_{3} = %.3f +/- %.3f", p8_Mx, p8_err_Mx));
    Mx->Draw();

    // Add chi2/ndf in the top left
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, 
        Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2_Mx, ndf_Mx, chi2_ndf_Mx));

    // Save the canvas
    c1->SaveAs("output/one_dimensional.pdf");
    // Clean up
    nh3->Close();
    carbon->Close();
    delete c1;

    std::cout << std::endl << std::endl << std::endl;
    std::cout << "if (prefix == \"x\") { return " << p0_x << 
        "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" <<
        std::endl;

    std::cout << "if (prefix == \"z\") { return " << p0_z << 
        "+" << p1_z << "*currentVariable+" << p2_z << "*std::pow(currentVariable,2)+" <<
        p3_z << "*std::pow(currentVariable,3); }" << std::endl;

    std::cout << "if (prefix == \"zeta\") { return " << p0_zeta << 
        "+" << p1_zeta << "*currentVariable+" << p2_zeta << "*std::pow(currentVariable,2)+" <<
        p3_zeta << "*std::pow(currentVariable,3)+" << p4_zeta << 
        "*std::pow(currentVariable,4); }" << std::endl;

    std::cout << "if (prefix == \"PT\") { return " << p0 << 
        "+" << p1 << "*currentVariable+" << p2 << "*std::pow(currentVariable,2); }" << std::endl;

    std::cout << "if (prefix == \"xF\") { return " << p0_xF << 
        "+" << p1_xF << "*currentVariable+" << p2_xF << "*std::pow(currentVariable,2)+" <<
        p3_xF << "*std::pow(currentVariable,3); }" << std::endl;

    std::cout << "if (prefix == \"Mx\") { return " 
          << p0_Mx << "*std::exp(-0.5*std::pow((currentVariable-" << p1_Mx << ")/" << p2_Mx << ",2)) + " 
          << p3_Mx << "*std::exp(-0.5*std::pow((currentVariable-" << p4_Mx << ")/" << p5_Mx << ",2)) + " 
          << p6_Mx << "*std::exp(-0.5*std::pow((currentVariable-" << p7_Mx << ")/" << p8_Mx << ",2)); }" 
          << std::endl;

    return 0;
}

double multi_dimensional(const char* nh3_file, const char* c_file, std::pair<double, double> fit_constant) {
    double scale_factor = fit_constant.first;
    double scale_error = fit_constant.second;

    // Suppress warnings
    // gErrorIgnoreLevel = kError;

    // Open the ROOT files
    TFile *nh3 = TFile::Open(nh3_file);
    TFile *carbon = TFile::Open(c_file);
    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        return 0;
    }

    for (int k = 0; k < 4; ++k) {

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

    std::string canvasName = "c1_" + std::to_string(k);
    TCanvas *c1 = new TCanvas(canvasName.c_str(), "Dilution Factor Analysis", 1600, 2000);
    if (k==0) {
        c1->Divide(3, 5);
    } else if (k==1) {
        c1->Divide(4, 5);
    } else {
        c1->Divide(5, 5);
    }

    std::string y_title;
    std::string y_range;
    switch (k) {
        case 0:
            y_title = "0.30<y<0.45";
            y_range = "0.30 < y && y < 0.45";
            break;
        case 1:
            y_title = "0.45<y<0.55";
            y_range = "0.45 < y && y < 0.55";
            break;
        case 2:
            y_title = "0.55<y<0.65";
            y_range = "0.55 < y && y < 0.65";
            break;
        case 3:
            y_title = "0.65<y<0.75";
            y_range = "0.65 < y && y < 0.75";
            break;
    }

    int max_Q2_bin = 5;
    if (k==0) {
        max_Q2_bin = 3;
    } else if (k==1) {
        max_Q2_bin = 4;
    }
    for (int j = 0; j < max_Q2_bin; ++j) {
        std::string Q2_range;
        std::string Q2y_prefix;
        std::string Q2_title;
        switch (j) {
            case 0:
                Q2_range = "1.00<Q2 && Q2<2.00";
                Q2_title = "1<Q^{2}<2";
                switch (k) {
                    case 0: 
                        Q2y_prefix = "Q2y4";
                        break;
                    case 1: 
                        Q2y_prefix = "Q2y3";
                        break;
                    case 2: 
                        Q2y_prefix = "Q2y2";
                        break;
                    case 3: 
                        Q2y_prefix = "Q2y1";
                        break;
                }
                break;
            case 1:
                Q2_range = "2.00<Q2 && Q2<3.00";
                Q2_title = "2<Q^{2}<3";
                switch (k) {
                    case 0: 
                        Q2y_prefix = "Q2y8";
                        break;
                    case 1: 
                        Q2y_prefix = "Q2y7";
                        break;
                    case 2: 
                        Q2y_prefix = "Q2y6";
                        break;
                    case 3: 
                        Q2y_prefix = "Q2y5";
                        break;
                }
                break;
            case 2:
                Q2_range = "3.00<Q2 && Q2<4.00";
                Q2_title = "3<Q^{2}<4";
                switch (k) {
                    case 0: 
                        Q2y_prefix = "Q2y12";
                        break;
                    case 1: 
                        Q2y_prefix = "Q2y11";
                        break;
                    case 2: 
                        Q2y_prefix = "Q2y10";
                        break;
                    case 3: 
                        Q2y_prefix = "Q2y9";
                        break;
                }
                break;
            case 3:
                Q2_range = "4.00<Q2 && Q2<5.00";
                Q2_title = "4<Q^{2}<5";
                switch (k) {
                    case 1: 
                        Q2y_prefix = "Q2y15";
                        break;
                    case 2: 
                        Q2y_prefix = "Q2y14";
                        break;
                    case 3: 
                        Q2y_prefix = "Q2y13";
                        break;
                }
                break;
            case 4:
                Q2_range = "5.00<Q2 && Q2<7.00";
                Q2_title = "5<Q^{2}<7";
                switch (k) {
                    case 2: 
                        Q2y_prefix = "Q2y17";
                        break;
                    case 3: 
                        Q2y_prefix = "Q2y16";
                        break;
                }
                break;
        }
    for (int i = 0; i < 5; ++i) {
        std::string z_range;
        std::string z_prefix;
        std::string z_title;
        switch (i) {
            case 0:
                z_range = "0.10<z && z<0.25";
                z_title = "0.10<z<0.25";
                z_prefix = "z1";
                break;
            case 1:
                z_range = "0.25<z && z<0.35";
                z_title = "0.25<z<0.35";
                z_prefix = "z2";
                break;
            case 2:
                z_range = "0.35<z && z<0.45";
                z_title = "0.35<z<0.45";
                z_prefix = "z3";
                break;
            case 3:
                z_range = "0.45<z && z<0.55";
                z_title = "0.45<z<0.55";
                z_prefix = "z4";
                break;
            case 4:
                z_range = "0.55<z && z<0.75";
                z_title = "0.55<z<0.75";
                z_prefix = "z5";
                break;
        }

        std::string cuts = "Mx>1.4 && "+Q2_range+" && "+y_range+" && "+z_range;
        c1->cd(5*j+(i+1)); // Pads are numbered from 1 to 25
        gPad->SetLeftMargin(0.15);

        // Create unique names for histograms and graphs
        std::string h_pT_nh3_name = "h_pT_nh3_" + std::to_string(k) + std::to_string(j) + std::to_string(i);
        std::string h_pT_carbon_name = "h_pT_carbon_" + std::to_string(k) + std::to_string(j) + std::to_string(i);
        std::string h_pT_carbon_scaled_name = "h_pT_carbon_scaled_" + std::to_string(k) + std::to_string(j) + std::to_string(i);
        std::string gr_dilution_name = "gr_dilution_" + std::to_string(k) + std::to_string(j) + std::to_string(i);
        std::string fit_poly_name = "fit_poly_" + std::to_string(k) + std::to_string(j) + std::to_string(i);
        std::string pave_text_name = "pave_text_" + std::to_string(k) + std::to_string(j) + std::to_string(i);

        // Create histograms
        TH1D *h_pT_nh3 = new 
            TH1D(h_pT_nh3_name.c_str(), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
        TH1D *h_pT_carbon = new 
            TH1D(h_pT_carbon_name.c_str(), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);

        // Draw histograms
        tree_nh3->Draw(("pT>>" + h_pT_nh3_name).c_str(), cuts.c_str());
        tree_carbon->Draw(("pT>>" + h_pT_carbon_name).c_str(), cuts.c_str());

        // Clone and scale carbon histogram
        TH1D *h_pT_carbon_scaled = (TH1D*)h_pT_carbon->Clone(h_pT_carbon_scaled_name.c_str());
        h_pT_carbon_scaled->SetTitle("P_{T} Distribution; P_{T} (GeV); Counts (Scaled)");

        for (int j = 1; j <= h_pT_carbon->GetNbinsX(); ++j) {
            double bin_content = h_pT_carbon->GetBinContent(j);
            double bin_error = h_pT_carbon->GetBinError(j);

            double new_content = bin_content * scale_factor;
            double new_error = new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + (scale_error / scale_factor) * (scale_error / scale_factor));
            h_pT_carbon_scaled->SetBinContent(j, new_content);
            h_pT_carbon_scaled->SetBinError(j, new_error);
        }

        // Create TGraphErrors for dilution factor
        TGraphErrors *gr_dilution = new TGraphErrors();
        gr_dilution->SetName(gr_dilution_name.c_str());

        for (int j = 1; j <= h_pT_nh3->GetNbinsX(); ++j) {
            double nh3_counts = h_pT_nh3->GetBinContent(j);
            double nh3_error = h_pT_nh3->GetBinError(j);
            double c_counts = h_pT_carbon_scaled->GetBinContent(j);
            double c_error = h_pT_carbon_scaled->GetBinError(j);

            if (nh3_counts > 0) {
                double dilution = (nh3_counts - c_counts) / nh3_counts;
                double dilution_error = std::sqrt(std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) + std::pow(1.0 / nh3_counts * c_error, 2));
                gr_dilution->SetPoint(j - 1, h_pT_nh3->GetBinCenter(j), dilution);
                gr_dilution->SetPointError(j - 1, 0, dilution_error);
            }
        }

        // Use the reformatted strings in the title
        std::string title = Q2_title + " , " + y_title + " , " + z_title;
        gr_dilution->SetTitle((title + "; P_{T} (GeV); D_{f} = (NH3 - s*C) / NH3").c_str());
        gr_dilution->SetTitle((title + "; P_{T} (GeV); D_{f} = (NH3 - s*C) / NH3").c_str());
        gr_dilution->SetMarkerStyle(20);
        // Draw the graph and set axis ranges
        gr_dilution->Draw("AP");
        gr_dilution->GetXaxis()->SetLimits(0, 1);
        gr_dilution->GetXaxis()->SetRangeUser(0, 1);
        gr_dilution->GetYaxis()->SetRangeUser(0.00, 0.30);

        // Increase the size of axis labels and titles
        gr_dilution->GetXaxis()->SetTitleSize(0.05);  // Increase title size
        gr_dilution->GetYaxis()->SetTitleSize(0.05);  // Increase title size
        gr_dilution->GetXaxis()->SetLabelSize(0.04);  // Increase label size
        gr_dilution->GetYaxis()->SetLabelSize(0.04);  // Increase label size

        // Fit to a polynomial
        TF1 *fit_poly = new TF1(fit_poly_name.c_str(), "[0]", 0, 1.0);
        // TF1 *fit_poly = new TF1(fit_poly_name.c_str(), "[0] + [1]*x + [2]*x^2", 0, 1.0);
        gr_dilution->Fit(fit_poly, "RQ");
        fit_poly->SetLineColor(kRed);
        fit_poly->Draw("SAME");

        // Retrieve fit parameters and their errors
        double p0 = fit_poly->GetParameter(0);
        double p0_err = fit_poly->GetParError(0);
        // double p1 = fit_poly->GetParameter(1);
        // double p1_err = fit_poly->GetParError(1);
        // double p2 = fit_poly->GetParameter(2);
        // double p2_err = fit_poly->GetParError(2);

        // Retrieve chi2 and NDF
        double chi2 = fit_poly->GetChisquare();
        int ndf = fit_poly->GetNDF();
        double chi2_ndf = chi2 / ndf;

        // Add fit parameters box
        TPaveText *pt = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
        pt->SetName(pave_text_name.c_str());
        pt->SetBorderSize(1);
        pt->SetFillStyle(1001); // Solid fill style
        pt->SetFillColor(kWhite); // White background
        pt->AddText(Form("p0 = %.3f +/- %.3f", p0, p0_err));
        // pt->AddText(Form("p1 = %.3f +/- %.3f", p1, p1_err));
        // pt->AddText(Form("p2 = %.3f +/- %.3f", p2, p2_err));
        pt->Draw();

        // Add chi2/ndf in the top left
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2_ndf));

        // Print the fit formula
        std::cout << "if (prefix == \"" << Q2y_prefix << z_prefix << "\") { return " << p0 << "; }" << std::endl << std::endl;
        // std::cout << "if (prefix == \"" << z_prefix << "\") { return " << p0 << "+" << p1 << "*currentVariable+" << p2 << "*std::pow(currentVariable,2); }" << std::endl;
    }
    }

    // Save the canvas
    std::string file = "output/" + std::to_string(k) + ".pdf";
    c1->SaveAs(file.c_str());
    }

    // Clean up
    nh3->Close();
    carbon->Close();
    // delete c1;

    return 0;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file>" << std::endl;
        return 1;
    }

    std::pair<double, double> fit_constant = scale_normalization(argv[1], argv[2]);
    one_dimensional(argv[1], argv[2], fit_constant);
    multi_dimensional(argv[1], argv[2], fit_constant);
}