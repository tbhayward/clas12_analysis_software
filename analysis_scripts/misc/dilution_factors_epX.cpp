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
#include <TLine.h>
#include <sstream>

// Function to create and fill histograms
void createAndFillHistograms(TTree* tree, const char* branchName, TH1D* hist) {
    tree->Draw(Form("%s>>%s", branchName, hist->GetName()));
}

// Function to calculate ratios and fit
TGraphErrors* calculateRatiosAndFit(TH1D* h1, TH1D* h2, double fitStart, double fitEnd, double plotMin, double plotMax, TF1* &fit_const, TF1* &dotted_line) {
    TGraphErrors* gr_ratio = new TGraphErrors();
    for (int i = 1; i <= h1->GetNbinsX(); ++i) {
        double nh3_counts = h1->GetBinContent(i);
        double c_counts = h2->GetBinContent(i);
        if (c_counts > 0) {
            double ratio = nh3_counts / c_counts;
            double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
            gr_ratio->SetPoint(i - 1, h1->GetBinCenter(i), ratio);
            gr_ratio->SetPointError(i - 1, 0, error);
        }
    }
    gr_ratio->GetYaxis()->SetRangeUser(plotMin, plotMax);

    fit_const = new TF1("fit_const", "[0]", fitStart, fitEnd);
    gr_ratio->Fit(fit_const, "R");

    dotted_line = new TF1("dotted_line", "[0]", fitEnd, h1->GetXaxis()->GetXmax());
    dotted_line->SetParameter(0, fit_const->GetParameter(0));
    dotted_line->SetLineColor(kRed);
    dotted_line->SetLineStyle(7);

    return gr_ratio;
}

// Function to plot histograms and fits
void plotHistogramsAndFits(TCanvas* canvas, TH1D* h1, TH1D* h2, TGraphErrors* gr_ratio, TF1* fit_const, TF1* dotted_line, const char* title) {
    canvas->cd();
    gPad->SetLeftMargin(0.15);
    gPad->SetLogy();
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);
    h1->Draw();
    h2->Draw("SAME");

    TLegend* leg = new TLegend(0.75, 0.8, 0.9, 0.9);
    leg->AddEntry(h1, "NH_{3}", "l");
    leg->AddEntry(h2, "C", "l");
    leg->Draw();

    h1->SetStats(0);
    h2->SetStats(0);

    canvas->cd();
    gPad->SetLeftMargin(0.15);
    gr_ratio->SetTitle(title);
    gr_ratio->SetMarkerStyle(20);
    gr_ratio->Draw("AP");
    fit_const->SetLineColor(kRed);
    fit_const->Draw("SAME");
    dotted_line->Draw("SAME");
}

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

    // First Canvas
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(2, 2);

    // Histograms for Mx
    TH1D *h_Mx_nh3 = new TH1D("h_Mx_nh3", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);
    TH1D *h_Mx_carbon = new TH1D("h_Mx_carbon", "M_{x} Distribution; M_{x} (GeV); Counts", 100, -2, 3);
    createAndFillHistograms(tree_nh3, "Mx", h_Mx_nh3);
    createAndFillHistograms(tree_carbon, "Mx", h_Mx_carbon);

    // Fit and plot for Mx
    TF1 *fit_const_Mx, *dotted_line_Mx;
    TGraphErrors* gr_ratio_Mx = calculateRatiosAndFit(h_Mx_nh3, h_Mx_carbon, -2, -0.5, 9, 15, fit_const_Mx, dotted_line_Mx);
    plotHistogramsAndFits(c1->cd(1), h_Mx_nh3, h_Mx_carbon, gr_ratio_Mx, fit_const_Mx, dotted_line_Mx, "NH_{3} to Carbon Ratio; M_{x} (GeV); Ratio");

    // Histograms for xF
    TH1D *h_xF_nh3 = new TH1D("h_xF_nh3", "x_{F} Distribution; x_{F}; Counts", 100, -2.5, 1);
    TH1D *h_xF_carbon = new TH1D("h_xF_carbon", "x_{F} Distribution; x_{F}; Counts", 100, -2.5, 1);
    createAndFillHistograms(tree_nh3, "xF", h_xF_nh3);
    createAndFillHistograms(tree_carbon, "xF", h_xF_carbon);

    // Fit and plot for xF
    TF1 *fit_const_xF, *dotted_line_xF;
    TGraphErrors* gr_ratio_xF = calculateRatiosAndFit(h_xF_nh3, h_xF_carbon, -2, -0.5, 5, 15, fit_const_xF, dotted_line_xF);
    plotHistogramsAndFits(c1->cd(2), h_xF_nh3, h_xF_carbon, gr_ratio_xF, fit_const_xF, dotted_line_xF, "NH_{3} to Carbon Ratio; x_{F}; Ratio");

    // Second Canvas
    TCanvas *c2 = new TCanvas("c2", "Dilution Factor Analysis", 1200, 1600);
    c2->Divide(2, 4);

    // Functions to create and fill pT and other histograms
    auto createAndFillHistogramsForPT = [&](TTree* tree, const char* branchName, TH1D* hist) {
        tree->Draw(Form("%s>>%s", branchName, hist->GetName()));
    };

    auto calculateRatiosAndFitForPT = [&](TH1D* h1, TH1D* h2, double fitStart, double fitEnd, double plotMin, double plotMax, TF1* &fit_const, TF1* &dotted_line) {
        TGraphErrors* gr_ratio = new TGraphErrors();
        for (int i = 1; i <= h1->GetNbinsX(); ++i) {
            double nh3_counts = h1->GetBinContent(i);
            double c_counts = h2->GetBinContent(i);
            if (c_counts > 0) {
                double ratio = nh3_counts / c_counts;
                double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
                gr_ratio->SetPoint(i - 1, h1->GetBinCenter(i), ratio);
                gr_ratio->SetPointError(i - 1, 0, error);
            }
        }
        gr_ratio->GetYaxis()->SetRangeUser(plotMin, plotMax);

        fit_const = new TF1("fit_const", "[0]", fitStart, fitEnd);
        gr_ratio->Fit(fit_const, “R”);
        dotted_line = new TF1("dotted_line", "[0]", fitEnd, h1->GetXaxis()->GetXmax());
        dotted_line->SetParameter(0, fit_const->GetParameter(0));
        dotted_line->SetLineColor(kRed);
        dotted_line->SetLineStyle(7);

        return gr_ratio;
    };

    auto plotHistogramsAndFitsForPT = [&](TCanvas* canvas, TH1D* h1, TH1D* h2, TGraphErrors* gr_ratio, TF1* fit_const, TF1* dotted_line, const char* title) {
        canvas->cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetLogy();
        h1->SetLineColor(kBlue);
        h2->SetLineColor(kRed);
        h1->Draw();
        h2->Draw("SAME");

        TLegend* leg = new TLegend(0.75, 0.8, 0.9, 0.9);
        leg->AddEntry(h1, "NH_{3}", "l");
        leg->AddEntry(h2, "C", "l");
        leg->Draw();

        h1->SetStats(0);
        h2->SetStats(0);

        canvas->cd();
        gPad->SetLeftMargin(0.15);
        gr_ratio->SetTitle(title);
        gr_ratio->SetMarkerStyle(20);
        gr_ratio->Draw("AP");
        fit_const->SetLineColor(kRed);
        fit_const->Draw("SAME");
        dotted_line->Draw("SAME");
    };

    // Histograms for pT
    TH1D *h_pT_nh3 = new TH1D("h_pT_nh3", "P_{T} Distribution; P_{T} (GeV); Counts", 100, 0, 1.0);
    TH1D *h_pT_carbon = new TH1D("h_pT_carbon", "P_{T} Distribution; P_{T} (GeV); Counts", 100, 0, 1.0);
    createAndFillHistogramsForPT(tree_nh3, "pT", h_pT_nh3);
    createAndFillHistogramsForPT(tree_carbon, "pT", h_pT_carbon);

    // Fit and plot for pT
    TF1 *fit_const_pT, *dotted_line_pT;
    TGraphErrors* gr_ratio_pT = calculateRatiosAndFitForPT(h_pT_nh3, h_pT_carbon, -2, -0.5, 5, 15, fit_const_pT, dotted_line_pT);
    plotHistogramsAndFitsForPT(c2->cd(1), h_pT_nh3, h_pT_carbon, gr_ratio_pT, fit_const_pT, dotted_line_pT, "NH_{3} to Carbon Ratio; P_{T} (GeV); Ratio");

    // Add black divider lines to canvas
    c2->cd();
    TLine *line = new TLine();
    line->SetLineColor(kBlack);
    line->SetLineWidth(2);
    line->DrawLineNDC(0.5, 0, 0.5, 1);
    line->DrawLineNDC(0, 0.5, 1, 0.5);

    // Save the canvases
    c1->SaveAs("dilution_factors_Mx_xF.pdf");
    c2->SaveAs("dilution_factors_pT.pdf");

    // Clean up
    nh3->Close();
    carbon->Close();
    // Clean up objects
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file>" << std::endl;
        return 1;
    }
    dilution_factors_epX(argv[1], argv[2]);
    return 0;
}