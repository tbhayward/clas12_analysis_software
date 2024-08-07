#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TF1.h>
#include <TROOT.h>
#include <iostream>
#include <cmath>
#include <utility>

// Constants for the dilution factor calculation
const double L_C = 1.5;
const double L = 5.86;
const double L_CH = 3;
const double rho_He = 0.145 / 4;
const double rho_C = 2.0 / 12;
const double rho_CH = 1.0 / 14;
const double rho_A = 0.92 / 17;

// Fractional charge values
const double xA = 0.23895;
const double xC = 0.08751;
const double xCH = 0.09614;
const double xHe = 0.33822;
const double xf = 0.23918;

double calculate_dilution_factor(double nA, double nC, double nCH, double nMT, double nf) {
    return (23.0 * (-nMT * xA + nA * xHe) * 
            (-0.511667 * nMT * xC * xCH * xf + 
             (1.0 * nf * xC * xCH - 
              3.41833 * nCH * xC * xf + 
              2.93 * nC * xCH * xf) * xHe)) / 
           (nA * xHe * 
            (62.6461 * nMT * xC * xCH * xf + 
             1.0 * nf * xC * xCH * xHe - 
             78.6217 * nCH * xC * xf * xHe + 
             14.9756 * nC * xCH * xf * xHe));
}

void plot_dilution_factor(const char* variable_name, const char* x_title, double x_min, double x_max, int n_bins, 
                          TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad) {
    canvas->cd(pad);
    gPad->SetLeftMargin(0.15);

    TH1D *h_nh3 = new TH1D(Form("h_%s_nh3", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_c = new TH1D(Form("h_%s_c", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_ch = new TH1D(Form("h_%s_ch", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_he = new TH1D(Form("h_%s_he", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_empty = new TH1D(Form("h_%s_empty", variable_name), "", n_bins, x_min, x_max);

    nh3->Draw(Form("%s>>h_%s_nh3", variable_name, variable_name));
    c->Draw(Form("%s>>h_%s_c", variable_name, variable_name));
    ch->Draw(Form("%s>>h_%s_ch", variable_name, variable_name));
    he->Draw(Form("%s>>h_%s_he", variable_name, variable_name));
    empty->Draw(Form("%s>>h_%s_empty", variable_name, variable_name));

    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= n_bins; ++i) {
        double nA = h_nh3->GetBinContent(i);
        double nC = h_c->GetBinContent(i);
        double nCH = h_ch->GetBinContent(i);
        double nMT = h_he->GetBinContent(i);
        double nf = h_empty->GetBinContent(i);

        double dilution = calculate_dilution_factor(nA, nC, nCH, nMT, nf);
        double error = 0.01; // Placeholder error

        gr_dilution->SetPoint(i - 1, h_nh3->GetBinCenter(i), dilution);
        gr_dilution->SetPointError(i - 1, 0, error);
    }

    gr_dilution->SetTitle(Form(";%s;D_{f}", x_title));
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");
    gr_dilution->GetXaxis()->SetRangeUser(x_min, x_max);
    gr_dilution->GetYaxis()->SetRangeUser(0.10, 0.30);

    // Fit and plot
    TF1 *fit_func = new TF1("fit_func", "[0] + [1]*x + [2]*x^2", x_min, x_max);
    gr_dilution->Fit(fit_func, "RQ");
    fit_func->SetLineColor(kRed);
    fit_func->Draw("SAME");

    // Print fit parameters and chi2/ndf
    double chi2 = fit_func->GetChisquare();
    int ndf = fit_func->GetNDF();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2 / ndf));

    // Add fit parameters box
    TPaveText *pt = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
    pt->SetBorderSize(1);
    pt->SetFillStyle(1001);
    pt->SetFillColor(kWhite);
    pt->AddText(Form("p0 = %.3f +/- %.3f", fit_func->GetParameter(0), fit_func->GetParError(0)));
    pt->AddText(Form("p1 = %.3f +/- %.3f", fit_func->GetParameter(1), fit_func->GetParError(1)));
    pt->AddText(Form("p2 = %.3f +/- %.3f", fit_func->GetParameter(2), fit_func->GetParError(2)));
    pt->Draw();

    // Clean up histograms
    delete h_nh3;
    delete h_c;
    delete h_ch;
    delete h_he;
    delete h_empty;
}

std::pair<TF1*, TGraphErrors*> fit_and_plot_dilution(const char* variable_name, const char* x_title, double x_min, double x_max, int n_bins, 
                          TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad) {
    // Call the plotting function
    plot_dilution_factor(variable_name, x_title, x_min, x_max, n_bins, nh3, c, ch, he, empty, canvas, pad);
    
    // Return the fit function and graph
    TF1* fit_func = (TF1*)gROOT->FindObject("fit_func");
    TGraphErrors* gr_dilution = (TGraphErrors*)gROOT->FindObject("gr_dilution");
    
    return std::make_pair(fit_func, gr_dilution);
}

void one_dimensional(TFile* nh3_file, TFile* c_file, TFile* ch_file, TFile* he_file, TFile* empty_file) {
    // Get the PhysicsEvents trees
    TTree* nh3 = (TTree*)nh3_file->Get("PhysicsEvents");
    TTree* c = (TTree*)c_file->Get("PhysicsEvents");
    TTree* ch = (TTree*)ch_file->Get("PhysicsEvents");
    TTree* he = (TTree*)he_file->Get("PhysicsEvents");
    TTree* empty = (TTree*)empty_file->Get("PhysicsEvents");

    // Create a canvas and divide it into 1 row and 3 columns
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1600, 600);
    c1->Divide(3, 1);

    // Fit and plot for x-Bjorken
    auto fit_x = fit_and_plot_dilution("x", "x_{B} (GeV)", 0.06, 0.6, 50, nh3, c, ch, he, empty, c1, 1);

    // Fit and plot for transverse momentum
    auto fit_pT = fit_and_plot_dilution("pT", "P_{T} (GeV)", 0, 1.0, 50, nh3, c, ch, he, empty, c1, 2);

    // Fit and plot for x-Feynman
    auto fit_xF = fit_and_plot_dilution("xF", "x_{F} (GeV)", -0.8, 0.5, 50, nh3, c, ch, he, empty, c1, 3);

    // Save the canvas as a PNG file
    c1->SaveAs("output/one_dimensional.png");

    // Print the fit functions for each variable
    std::cout << "Fit for x_Bjorken: " << fit_x.first->GetExpFormula("p") << std::endl;
    std::cout << "Fit for P_T: " << fit_pT.first->GetExpFormula("p") << std::endl;
    std::cout << "Fit for x_Feynman: " << fit_xF.first->GetExpFormula("p") << std::endl;

    // Clean up
    delete c1;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <NH3 ROOT file> <Carbon ROOT file> <CH2 ROOT file> <Helium ROOT file> <Empty ROOT file>" << std::endl;
        return 1;
    }

    // Open the ROOT files
    TFile *nh3 = TFile::Open(argv[1]);
    TFile *c = TFile::Open(argv[2]);
    TFile *ch = TFile::Open(argv[3]);
    TFile *he = TFile::Open(argv[4]);
    TFile *empty = TFile::Open(argv[5]);

    // Check if files opened successfully
    if (!nh3 || nh3->IsZombie() || !c || c->IsZombie() || !ch || ch->IsZombie() || !he || he->IsZombie() || !empty || empty->IsZombie()) {
        std::cerr << "Error opening one or more files!" << std::endl;
        if (nh3) nh3->Close();
        if (c) c->Close();
        if (ch) ch->Close();
        if (he) he->Close();
        if (empty) empty->Close();
        return 2;
    }

    // Call the one-dimensional function
    one_dimensional(nh3, c, ch, he, empty);

    // Safely close the ROOT files
    nh3->Close();
    c->Close();
    ch->Close();
    he->Close();
    empty->Close();

    return 0;
}