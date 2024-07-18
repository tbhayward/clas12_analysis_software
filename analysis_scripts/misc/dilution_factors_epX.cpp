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

void plot_histograms_and_ratios(TTree* tree_nh3, TTree* tree_carbon, const char* branch_name, const char* title, 
                                double hist_min, double hist_max, double fit_min, double fit_max, 
                                double plot_y_min, double plot_y_max, TCanvas* canvas, int pad_hist, int pad_ratio) {
    // Create histograms
    std::string hist_name_nh3 = std::string(branch_name) + "_nh3";
    std::string hist_name_carbon = std::string(branch_name) + "_carbon";
    TH1D* h_nh3 = new TH1D(hist_name_nh3.c_str(), (title + std::string(" Distribution; ") + title + "; Counts").c_str(), 100, hist_min, hist_max);
    TH1D* h_carbon = new TH1D(hist_name_carbon.c_str(), (title + std::string(" Distribution; ") + title + "; Counts").c_str(), 100, hist_min, hist_max);

    // Fill the histograms
    tree_nh3->Draw(Form("%s>>%s", branch_name, hist_name_nh3.c_str()));
    tree_carbon->Draw(Form("%s>>%s", branch_name, hist_name_carbon.c_str()));

    // First panel: plot histograms
    canvas->cd(pad_hist);
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy(); // Log scale to better see differences
    h_nh3->SetLineColor(kBlue);
    h_carbon->SetLineColor(kRed);
    h_nh3->Draw();
    h_carbon->Draw("SAME");

    // Add legend
    TLegend* leg = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg->AddEntry(h_nh3, "NH_{3}", "l");
    leg->AddEntry(h_carbon, "C", "l");
    leg->Draw();

    // Remove statboxes
    h_nh3->SetStats(0);
    h_carbon->SetStats(0);

    // Second panel: ratio of NH3 to Carbon counts
    canvas->cd(pad_ratio);
    gPad->SetLeftMargin(0.15);
    TGraphErrors* gr_ratio = new TGraphErrors();
    for (int i = 1; i <= h_nh3->GetNbinsX(); ++i) {
        double nh3_counts = h_nh3->GetBinContent(i);
        double c_counts = h_carbon->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio->SetPoint(i - 1, h_nh3->GetBinCenter(i), ratio);
            gr_ratio->SetPointError(i - 1, 0, error);
        }
    }
    gr_ratio->GetYaxis()->SetRangeUser(plot_y_min, plot_y_max);
    gr_ratio->SetTitle((std::string("NH_{3} to Carbon Ratio; ") + title + "; Ratio").c_str());
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");

    // Fit the data to a constant
    TF1* fit_const = new TF1("fit_const", "[0]", fit_min, fit_max);
    gr_ratio->Fit(fit_const, "R");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");

    // Add the dotted-dashed line for the remaining range
    TF1* dotted_line = new TF1("dotted_line", "[0]", fit_max, hist_max);
    dotted_line->SetParameter(0, fit_const->GetParameter(0));
    dotted_line->SetLineColor(kRed);
    dotted_line->SetLineStyle(7);
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
}

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

    // Plot Mx histograms and ratios
    plot_histograms_and_ratios(tree_nh3, tree_carbon, "Mx", "M_{x} (GeV)", -2, 3, -2, -0.5, 9, 15, c1, 1, 2);

    // Plot xF histograms and ratios
    plot_histograms_and_ratios(tree_nh3, tree_carbon, "xF", "x_{F}", -2.5, 1, -2, -1.25, 5, 15, c1, 3, 4);

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