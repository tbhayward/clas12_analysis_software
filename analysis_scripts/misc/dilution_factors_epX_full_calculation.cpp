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
const double xA = 0.71526;
const double xC = 0.07111;
const double xCH = 0.03049;
const double xHe = 0.10728;
const double xf = 0.07586;


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

// double calculate_dilution_error(double nA, double nC, double nCH, double nMT, double nf) {
//     double term1 = 3988.9 * nA * nf * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) * pow(xHe, 2) * 
//                    pow(-1.0 * nMT * xA + nA * xHe, 2) * 
//                    pow(1.0 * nMT * xC * xCH - 1.19072 * nCH * xC * xHe + 0.190722 * nC * xCH * xHe, 2);
    
//     double term2 = 64705.8 * nA * nCH * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) * pow(xHe, 2) * 
//                    pow(-1.0 * nMT * xA + nA * xHe, 2) * 
//                    pow(1.0 * nMT * xC * xf - 0.295642 * nf * xC * xHe - 0.704358 * nC * xf * xHe, 2);
    
//     double term3 = 36563.4 * nA * nC * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) * pow(xHe, 2) * 
//                    pow(-1.0 * nMT * xA + nA * xHe, 2) * 
//                    pow(1.0 * nMT * xCH * xf - 0.0629946 * nf * xCH * xHe - 0.937005 * nCH * xf * xHe, 2);
    
//     double term4 = 1027.45 * pow(nMT, 2) * pow(xA, 2) * 
//                    pow(1.0 * nMT * xC * xCH * xf + 
//                        (-1.9544 * nf * xC * xCH + 6.68078 * nCH * xC * xf - 5.72638 * nC * xCH * xf) * xHe, 2) * 
//                    pow(1.0 * nMT * xC * xCH * xf + 
//                        (0.0159627 * nf * xC * xCH - 1.25501 * nCH * xC * xf + 0.23905 * nC * xCH * xf) * xHe, 2);
    
//     double term5 = 0.261803 * nA * nMT * 
//                    pow(62.6461 * pow(nMT, 2) * xA * pow(xC, 2) * pow(xCH, 2) * pow(xf, 2) + 
//                        nMT * xC * xCH * xf * 
//                        (2.0 * nf * xA * xC * xCH + 
//                         (-157.243 * nCH * xA * xC + 29.9511 * nC * xA * xCH + 
//                          7.10543e-15 * nA * xC * xCH) * xf) * xHe + 
//                        (-1.9544 * pow(nf, 2) * xA * pow(xC, 2) * pow(xCH, 2) + 
//                         nf * xC * xCH * 
//                         (160.339 * nCH * xA * xC - 34.9946 * nC * xA * xCH - 
//                          123.435 * nA * xC * xCH) * xf + 
//                         (-525.254 * pow(nCH, 2) * xA * pow(xC, 2) + 
//                          nCH * xC * 
//                          (550.266 * nC * xA + 497.147 * nA * xC) * xCH + 
//                          nC * 
//                          (-85.7558 * nC * xA - 373.711 * nA * xC) * pow(xCH, 2)) * pow(xf, 2)) * pow(xHe, 2), 2);
    
//     double denominator = pow(nA, 3) * pow(xHe, 2) * 
//                          pow(62.6461 * nMT * xC * xCH * xf + 1.0 * nf * xC * xCH * xHe - 
//                              78.6217 * nCH * xC * xf * xHe + 
//                              14.9756 * nC * xCH * xf * xHe, 4);
    
//     double sigma_df = 23.0 * sqrt(term1 + term2 + term3 + term4 + term5) / denominator;
    
//     return sigma_df;
// }

double calculate_new_dilution_error(double nA, double nC, double nCH, double nMT, double nf) {
    double term1 = 0.734694 * nA * nC * pow((nA - nMT), 2) * 
                   pow((2.77556e-17 * nC - 0.554976 * nCH - 0.0373109 * nf + 0.592287 * nMT), 2);

    double term2 = 0.456108 * nA * pow((nA - nMT), 2) * 
                   pow((0.704358 * nC + 0.295642 * nf - nMT), 2) * nMT;

    double term3 = 0.0281176 * nA * nf * pow((nA - nMT), 2) * 
                   pow((0.190722 * nC - 1.19072 * nCH + nMT), 2);

    double term4 = 0.0135714 * pow((nC - 1.16667 * nCH + 0.341297 * nf - 0.17463 * nMT), 2) * 
                   pow(nMT, 2) * pow((nC - 5.25 * nCH + 0.0667755 * nf + 4.18322 * nMT), 2);

    double term5 = 0.022405 * nA * nMT * 
                   pow((0.778288 * pow(nC, 2) + 4.76701 * pow(nCH, 2) - 1.45518 * nCH * nf + 0.0177374 * pow(nf, 2) + 
                        nC * (-4.99401 * nCH + 0.317598 * nf - 0.271825 * nMT) + 
                        nA * (3.39167 * nC - 4.51192 * nCH + 1.12025 * nf - 1.11022e-16 * nMT) + 
                        1.42708 * nCH * nMT - 0.0181513 * nf * nMT - 0.568553 * pow(nMT, 2)), 2);

    double denominator = pow(nA, 3) * 
                         pow((0.135913 * nC - 0.713541 * nCH + 0.00907563 * nf + 0.568553 * nMT), 4);

    double sigma_df = 0.713541 * sqrt(term1 + term2 + term3 + term4 + term5) / denominator;


    std::cout << term1 << " " << term2 << " " << term3 << " " << term4 << " " << term5 << " " << denominator << std::endl;
    std::cout << pow(term1 + term2 + term3 + term4 + term5,0.5) << " " << denominator << std::endl << std::endl;
    return sigma_df;
    // return term1;
}

double calculate_simple_error(double nh3_counts, double nh3_error, double c_counts, double c_error) {
    // Propagate the error using the simplified method
    double dilution_error = std::sqrt(
                std::pow((c_counts / (nh3_counts * nh3_counts)) * nh3_error, 2) +
                std::pow(1.0 / nh3_counts * c_error, 2));
    return dilution_error;
}

void plot_dilution_factor(const char* variable_name, const char* x_title, double x_min, double x_max, int n_bins, 
                          TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad, bool skip_fit = false) {
    canvas->cd(pad);
    gPad->SetLeftMargin(0.15);

    // Create histograms for data
    TH1D *h_nh3 = new TH1D(Form("h_%s_nh3", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_c = new TH1D(Form("h_%s_c", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_ch = new TH1D(Form("h_%s_ch", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_he = new TH1D(Form("h_%s_he", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_empty = new TH1D(Form("h_%s_empty", variable_name), "", n_bins, x_min, x_max);

    // Fill histograms with data
    nh3->Draw(Form("%s>>h_%s_nh3", variable_name, variable_name));
    c->Draw(Form("%s>>h_%s_c", variable_name, variable_name));
    ch->Draw(Form("%s>>h_%s_ch", variable_name, variable_name));
    he->Draw(Form("%s>>h_%s_he", variable_name, variable_name));
    empty->Draw(Form("%s>>h_%s_empty", variable_name, variable_name));

    // Scale carbon counts and update their errors
    // double s = 1;       // scale factor for carbon counts
    double s = 11.306;  // Uncomment when using full statistics
    // double s_error = 0.110;  // uncertainty in the scale factor
    double s_error = 0.0;  // uncertainty in the scale factor
    TH1D *h_c_scaled = (TH1D*)h_c->Clone(Form("h_%s_c_scaled", variable_name));
    for (int i = 1; i <= h_c->GetNbinsX(); ++i) {
        double bin_content = h_c->GetBinContent(i);
        double bin_error = h_c->GetBinError(i);

        double new_content = bin_content * s;
        double new_error = new_content * std::sqrt((bin_error / bin_content) * (bin_error / bin_content) + 
            (s_error / s) * (s_error / s));
        h_c_scaled->SetBinContent(i, new_content);
        h_c_scaled->SetBinError(i, new_error);
    }

    // Calculate dilution factor and its error
    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= n_bins; ++i) {
        double nA = h_nh3->GetBinContent(i);
        double nA_error = h_nh3->GetBinError(i);
        double nC = h_c->GetBinContent(i);
        double nC_scaled = h_c_scaled->GetBinContent(i);
        double nC_scaled_error = h_c_scaled->GetBinError(i);

        double nCH = h_ch->GetBinContent(i);
        double nMT = h_he->GetBinContent(i);
        double nf = h_empty->GetBinContent(i);

        double dilution = calculate_dilution_factor(nA, nC, nCH, nMT, nf);
        // double error = calculate_dilution_error(nA, nC, nCH, nMT, nf);
        double error = calculate_new_dilution_error(nA/1000000, nC/1000000, nCH/1000000, nMT/1000000, nf/1000000);
        std::cout << error << std::endl;
        // double error = calculate_simple_error(nA, nA_error, nC_scaled, nC_scaled_error);

        // For integrated plot, set the point at the center of the plot range
        double x_position = skip_fit ? (x_min + x_max) / 2 : h_nh3->GetBinCenter(i);

        gr_dilution->SetPoint(i - 1, x_position, dilution);
        gr_dilution->SetPointError(i - 1, 0, error);
    }

    gr_dilution->SetTitle(Form(";%s;D_{f}", x_title));
    gr_dilution->SetMarkerStyle(20);

    // Customizing the integrated plot
    if (skip_fit) {
        gr_dilution->GetXaxis()->SetRangeUser(0, 1);  // Set x-axis range manually
        gr_dilution->GetXaxis()->SetLabelSize(0);     // Remove x-axis labels
        gr_dilution->GetXaxis()->SetTickLength(0);    // Remove x-axis ticks
        gr_dilution->SetTitle(";Integrated;D_{f}");
    }
    
    gr_dilution->Draw("AP");
    gr_dilution->GetYaxis()->SetRangeUser(0.10, 0.30);

    // Fit and plot (skip fit for the integrated version)
    if (!skip_fit) {
        TF1 *fit_func = new TF1("fit_func", "[0] + [1]*x + [2]*x^2", x_min, x_max);
        gr_dilution->Fit(fit_func, "RQ");
        fit_func->SetLineColor(kRed);
        fit_func->Draw("SAME");

        // Print fit parameters and chi2/ndf
        double chi2 = fit_func->GetChisquare();
        int ndf = fit_func->GetNDF();
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035); // Decrease the font size
        latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2 / ndf));

        // Add fit parameters box
        TPaveText *pt = new TPaveText(0.55, 0.7, 0.9, 0.9, "brNDC");
        pt->SetBorderSize(1);
        pt->SetFillStyle(1001);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.035); // Decrease the font size
        pt->AddText(Form("p0 = %.3f +/- %.3f", fit_func->GetParameter(0), fit_func->GetParError(0)));
        pt->AddText(Form("p1 = %.3f +/- %.3f", fit_func->GetParameter(1), fit_func->GetParError(1)));
        pt->AddText(Form("p2 = %.3f +/- %.3f", fit_func->GetParameter(2), fit_func->GetParError(2)));
        pt->Draw();
    } else {
        // For integrated plot, display the value in the top right corner
        TPaveText *pt = new TPaveText(0.55, 0.7, 0.9, 0.9, "brNDC");
        pt->SetBorderSize(1);
        pt->SetFillStyle(1001);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.035); // Decrease the font size
        pt->AddText(Form("p0 = %.3f +/- %.3f", gr_dilution->GetY()[0], gr_dilution->GetErrorY(0)));
        pt->Draw();
    }

    // Clean up histograms
    delete h_nh3;
    delete h_c;
    delete h_ch;
    delete h_he;
    delete h_empty;
    delete h_c_scaled;
}

std::pair<TF1*, TGraphErrors*> fit_and_plot_dilution(const char* variable_name, const char* x_title, double x_min, double x_max, int n_bins,
TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad, bool skip_fit = false) {
    // Call the plotting function
    plot_dilution_factor(variable_name, x_title, x_min, x_max, n_bins, nh3, c, ch, he, empty, canvas, pad, skip_fit);
    // Return the fit function and graph
    TF1* fit_func = nullptr;
    TGraphErrors* gr_dilution = (TGraphErrors*)gPad->GetPrimitive("gr_dilution");

    if (!skip_fit) {
        fit_func = (TF1*)gPad->GetPrimitive("fit_func");
    }

    return std::make_pair(fit_func, gr_dilution);
}

void one_dimensional(TFile* nh3_file, TFile* c_file, TFile* ch_file, TFile* he_file, TFile* empty_file) {
    // Get the PhysicsEvents trees
    TTree* nh3 = (TTree*)nh3_file->Get("PhysicsEvents");
    TTree* c = (TTree*)c_file->Get("PhysicsEvents");
    TTree* ch = (TTree*)ch_file->Get("PhysicsEvents");
    TTree* he = (TTree*)he_file->Get("PhysicsEvents");
    TTree* empty = (TTree*)empty_file->Get("PhysicsEvents");
    // Create a canvas and divide it into 2 rows and 2 columns
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1600, 1200);
    c1->Divide(2, 2);

    // Integrated version (single bin)
    auto fit_integrated = fit_and_plot_dilution("x", "", 0.0, 1.0, 1, nh3, c, ch, he, empty, c1, 1, true);

    // Fit and plot for x-Bjorken
    auto fit_x = fit_and_plot_dilution("x", "x_{B} (GeV)", 0.06, 0.6, 25, nh3, c, ch, he, empty, c1, 2);

    // Fit and plot for transverse momentum
    auto fit_pT = fit_and_plot_dilution("pT", "P_{T} (GeV)", 0, 1.0, 25, nh3, c, ch, he, empty, c1, 3);
    // Fit and plot for x-Feynman
    auto fit_xF = fit_and_plot_dilution("xF", "x_{F} (GeV)", -0.8, 0.5, 25, nh3, c, ch, he, empty, c1, 4);

    // Save the canvas as a PNG file
    c1->SaveAs("output/one_dimensional.png");

    // Prepare to print the fit functions for each variable
    std::cout << std::endl << std::endl << std::endl;

    if (fit_x.first) {
        double p0_x = fit_x.first->GetParameter(0);
        double p1_x = fit_x.first->GetParameter(1);
        double p2_x = fit_x.first->GetParameter(2);
        std::cout << "if (prefix == \"x\") { return " << p0_x << 
            "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" << std::endl;
    }

    if (fit_pT.first) {
        double p0_PT = fit_pT.first->GetParameter(0);
        double p1_PT = fit_pT.first->GetParameter(1);
        double p2_PT = fit_pT.first->GetParameter(2);
        std::cout << "if (prefix == \"PT\") { return " << p0_PT << 
            "+" << p1_PT << "*currentVariable+" << p2_PT << "*std::pow(currentVariable,2); }" << std::endl;
    }

    if (fit_xF.first) {
        double p0_xF = fit_xF.first->GetParameter(0);
        double p1_xF = fit_xF.first->GetParameter(1);
        double p2_xF = fit_xF.first->GetParameter(2);
        std::cout << "if (prefix == \"xF\") { return " << p0_xF << 
            "+" << p1_xF << "*currentVariable+" << p2_xF << "*std::pow(currentVariable,2); }" << std::endl;
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