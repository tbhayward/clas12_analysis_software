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
                          TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad, 
                          std::map<std::string, std::string>& fit_results) {
    // Create histograms for each target
    TH1D* h_nh3 = new TH1D(Form("h_%s_nh3", variable_name), "", n_bins, x_min, x_max);
    TH1D* h_c = new TH1D(Form("h_%s_c", variable_name), "", n_bins, x_min, x_max);
    TH1D* h_ch = new TH1D(Form("h_%s_ch", variable_name), "", n_bins, x_min, x_max);
    TH1D* h_he = new TH1D(Form("h_%s_he", variable_name), "", n_bins, x_min, x_max);
    TH1D* h_empty = new TH1D(Form("h_%s_empty", variable_name), "", n_bins, x_min, x_max);

    // Fill histograms
    nh3->Draw(Form("%s>>h_%s_nh3", variable_name, variable_name), "Mx > 1.4");
    c->Draw(Form("%s>>h_%s_c", variable_name, variable_name), "Mx > 1.4");
    ch->Draw(Form("%s>>h_%s_ch", variable_name, variable_name), "Mx > 1.4");
    he->Draw(Form("%s>>h_%s_he", variable_name, variable_name), "Mx > 1.4");
    empty->Draw(Form("%s>>h_%s_empty", variable_name, variable_name), "Mx > 1.4");

    // Calculate dilution factor for each bin
    TGraphErrors* gr_dilution = new TGraphErrors();
    for (int i = 1; i <= h_nh3->GetNbinsX(); ++i) {
        double nA = h_nh3->GetBinContent(i);
        double nC = h_c->GetBinContent(i);
        double nCH = h_ch->GetBinContent(i);
        double nMT = h_he->GetBinContent(i);
        double nf = h_empty->GetBinContent(i);

        double dilution = calculate_dilution_factor(nA, nC, nCH, nMT, nf);
        gr_dilution->SetPoint(i - 1, h_nh3->GetBinCenter(i), dilution);
        gr_dilution->SetPointError(i - 1, 0, 0.01); // Placeholder for the error
    }

    // Draw on the canvas
    canvas->cd(pad);
    gPad->SetLeftMargin(0.15);
    gr_dilution->SetTitle(Form(";%s;D_{f}", x_title));
    gr_dilution->SetMarkerStyle(20);
    gr_dilution->Draw("AP");
    gr_dilution->GetXaxis()->SetRangeUser(x_min, x_max);
    gr_dilution->GetYaxis()->SetRangeUser(0.10, 0.30);

    // Fit to a fourth-degree polynomial
    TF1* fit_func = new TF1("fit_func", "[0] + [1]*x + [2]*x^2 + [3]*x^3 + [4]*x^4", x_min, x_max);
    gr_dilution->Fit(fit_func, "RQ");
    fit_func->SetLineColor(kRed);
    fit_func->Draw("SAME");

    // Retrieve fit parameters
    double p0 = fit_func->GetParameter(0);
    double p1 = fit_func->GetParameter(1);
    double p2 = fit_func->GetParameter(2);
    double p3 = fit_func->GetParameter(3);
    double p4 = fit_func->GetParameter(4);

    // Store the fit result in the map
    fit_results[variable_name] = Form("%f+%f*currentVariable+%f*std::pow(currentVariable,2)+%f*std::pow(currentVariable,3)+%f*std::pow(currentVariable,4)", 
                                      p0, p1, p2, p3, p4);

    // Optional: Add chi2/ndf in the top left
    double chi2 = fit_func->GetChisquare();
    int ndf = fit_func->GetNDF();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2 / ndf));

    // Clean up
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

void one_dimensional(TFile* nh3, TFile* c, TFile* ch, TFile* he, TFile* empty) {
    // Create a canvas and divide it into 1 row and 3 columns
    TCanvas* c1 = new TCanvas("c1", "Dilution Factor Analysis", 1600, 600);
    c1->Divide(3, 1);

    // Create a map to store the fit results
    std::map<std::string, std::string> fit_results;

    // Plot for x-Bjorken
    plot_dilution_factor("x", "x_{B} (GeV)", 0.06, 0.6, 50, nh3, c, ch, he, empty, c1, 1, fit_results);

    // Plot for transverse momentum
    plot_dilution_factor("pT", "P_{T} (GeV)", 0, 1.0, 50, nh3, c, ch, he, empty, c1, 2, fit_results);

    // Plot for x-Feynman
    plot_dilution_factor("xF", "x_{F} (GeV)", -0.8, 0.5, 50, nh3, c, ch, he, empty, c1, 3, fit_results);

    // Save the canvas as a PNG file
    c1->SaveAs("output/one_dimensional_placeholder.png");

    // Print the fit results in the desired format
    std::cout << "\n\n";
    if (fit_results.find("x") != fit_results.end()) {
        std::cout << "if (prefix == \"x\") { return " << fit_results["x"] << "; }\n";
    }
    if (fit_results.find("pT") != fit_results.end()) {
        std::cout << "if (prefix == \"pT\") { return " << fit_results["pT"] << "; }\n";
    }
    if (fit_results.find("xF") != fit_results.end()) {
        std::cout << "if (prefix == \"xF\") { return " << fit_results["xF"] << "; }\n";
    }

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