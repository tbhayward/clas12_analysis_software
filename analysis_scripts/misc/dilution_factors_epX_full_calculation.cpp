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

double calculate_dilution_error(double nA, double nC, double nCH, double nMT, double nf) {
    double term1 = 0.734694 * nA * nC * pow((nA - nMT), 2) * 
                   pow((- 0.554976 * nCH - 0.0373109 * nf + 0.592287 * nMT), 2);

    double term2 = 0.456108 * nA * pow((nA - nMT), 2) * 
                   pow((0.704358 * nC + 0.295642 * nf - nMT), 2) * nMT;

    double term3 = 0.0281176 * nA * nf * pow((nA - nMT), 2) * 
                   pow((0.190722 * nC - 1.19072 * nCH + nMT), 2);

    double term4 = 0.0135714 * pow((nC - 1.16667 * nCH + 0.341297 * nf - 0.17463 * nMT), 2) * 
                   pow(nMT, 2) * pow((nC - 5.25 * nCH + 0.0667755 * nf + 4.18322 * nMT), 2);

    double term5 = 0.022405 * nA * nMT * 
                   pow((0.778288 * pow(nC, 2) + 4.76701 * pow(nCH, 2) - 1.45518 * nCH * nf + 0.0177374 * pow(nf, 2) + 
                        nC * (-4.99401 * nCH + 0.317598 * nf - 0.271825 * nMT) + 
                        nA * (3.39167 * nC - 4.51192 * nCH + 1.12025 * nf) + 
                        1.42708 * nCH * nMT - 0.0181513 * nf * nMT - 0.568553 * pow(nMT, 2)), 2);

    double denominator = pow(nA, 3) * 
                         pow((0.135913 * nC - 0.713541 * nCH + 0.00907563 * nf + 0.568553 * nMT), 4);

    double sigma_df = 0.713541 * sqrt((term1 + term2 + term3 + term4 + term5 )/ denominator);

    return sigma_df;
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

    // Calculate dilution factor and its error
    TGraphErrors *gr_dilution = new TGraphErrors();
    for (int i = 1; i <= n_bins; ++i) {
        double nA = h_nh3->GetBinContent(i);
        double nC = h_c->GetBinContent(i);

        double nCH = h_ch->GetBinContent(i);
        double nMT = h_he->GetBinContent(i);
        double nf = h_empty->GetBinContent(i);

        double dilution = calculate_dilution_factor(nA, nC, nCH, nMT, nf);
        double error = calculate_dilution_error(nA/xA, nC/xC, nCH/xCH, nMT/xHe, nf/xf);

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

    double chi2_scale_factor = 1.0;

    // Fit and plot (skip fit for the integrated version)
    if (!skip_fit) {
        TF1 *fit_func = new TF1("fit_func", "[0] + [1]*x + [2]*x^2 + [3]*x^3", x_min, x_max);
        gr_dilution->Fit(fit_func, "RQ");
        fit_func->SetLineColor(kRed);

        // Calculate chi2/ndf scaling factor
        double chi2 = fit_func->GetChisquare();
        int ndf = fit_func->GetNDF();
        chi2_scale_factor = std::sqrt(chi2 / ndf);
        
        // Rescale the errors
        for (int i = 0; i < gr_dilution->GetN(); ++i) {
            double x, y;
            gr_dilution->GetPoint(i, x, y);
            gr_dilution->SetPointError(i, 0, gr_dilution->GetErrorY(i) * chi2_scale_factor);
        }

        // Refit with scaled errors
        gr_dilution->Fit(fit_func, "RQ");
        fit_func->Draw("SAME");

        // Print fit parameters and chi2/ndf
        chi2 = fit_func->GetChisquare();
        ndf = fit_func->GetNDF();
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
        pt->AddText(Form("p3 = %.3f +/- %.3f", fit_func->GetParameter(3), fit_func->GetParError(3)));
        pt->Draw();
    } else {
        // For integrated plot, scale errors by the average scale factor from other fits
        double avg_scale_factor = chi2_scale_factor;
        for (int i = 0; i < gr_dilution->GetN(); ++i) {
            double x, y;
            gr_dilution->GetPoint(i, x, y);
            gr_dilution->SetPointError(i, 0, gr_dilution->GetErrorY(i) * avg_scale_factor);
        }

        // For integrated plot, display the value in the top right corner
        TPaveText *pt = new TPaveText(0.55, 0.7, 0.9, 0.9, "brNDC");
        pt->SetBorderSize(1);
        pt->SetFillStyle(1001);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.035); // Decrease the font size
        pt->AddText(Form("p0 = %.4f +/- %.4f", gr_dilution->GetY()[0], gr_dilution->GetErrorY(0)));
        pt->Draw();
    }

    // Clean up histograms
    delete h_nh3;
    delete h_c;
    delete h_ch;
    delete h_he;
    delete h_empty;
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

double multi_dimensional(const char* nh3_file, const char* c_file, const char* ch_file, const char* he_file, const char* empty_file) {
    // Open the ROOT files
    TFile *nh3 = TFile::Open(nh3_file);
    TFile *carbon = TFile::Open(c_file);
    TFile *ch = TFile::Open(ch_file);
    TFile *he = TFile::Open(he_file);
    TFile *empty = TFile::Open(empty_file);

    if (!nh3 || nh3->IsZombie() || !carbon || carbon->IsZombie() || !ch || ch->IsZombie() || !he || he->IsZombie() || !empty || empty->IsZombie()) {
        std::cerr << "Error opening files!" << std::endl;
        if (nh3) nh3->Close();
        if (carbon) carbon->Close();
        if (ch) ch->Close();
        if (he) he->Close();
        if (empty) empty->Close();
        return 0;
    }

    for (int k = 0; k < 4; ++k) {
        // Get the PhysicsEvents trees
        TTree *tree_nh3;
        TTree *tree_carbon;
        TTree *tree_ch;
        TTree *tree_he;
        TTree *tree_empty;
        nh3->GetObject("PhysicsEvents", tree_nh3);
        carbon->GetObject("PhysicsEvents", tree_carbon);
        ch->GetObject("PhysicsEvents", tree_ch);
        he->GetObject("PhysicsEvents", tree_he);
        empty->GetObject("PhysicsEvents", tree_empty);

        if (!tree_nh3 || !tree_carbon || !tree_ch || !tree_he || !tree_empty) {
            std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
            nh3->Close();
            carbon->Close();
            ch->Close();
            he->Close();
            empty->Close();
            return 0;
        }

        std::string canvasName = "c1_" + std::to_string(k);
        TCanvas *c1 = new TCanvas(canvasName.c_str(), "Dilution Factor Analysis", 1600, 2000);
        if (k == 0) {
            c1->Divide(3, 5);
        } else if (k == 1) {
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
        if (k == 0) {
            max_Q2_bin = 3;
        } else if (k == 1) {
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

                std::string cuts = "Mx>1.4 && " + Q2_range + " && " + y_range + " && " + z_range;
                c1->cd(5 * j + (i + 1)); // Pads are numbered from 1 to 25
                gPad->SetLeftMargin(0.15);

                // Create histograms for different targets
                TH1D *h_pT_nh3 = new TH1D(Form("h_pT_nh3_%d%d%d", k, j, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
                TH1D *h_pT_c = new TH1D(Form("h_pT_c_%d%d%d", k, j, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
                TH1D *h_pT_ch = new TH1D(Form("h_pT_ch_%d%d%d", k, j, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
                TH1D *h_pT_he = new TH1D(Form("h_pT_he_%d%d%d", k, j, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
                TH1D *h_pT_empty = new TH1D(Form("h_pT_empty_%d%d%d", k, j, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
                // Draw histograms
                tree_nh3->Draw(Form("pT>>h_pT_nh3_%d%d%d", k, j, i), cuts.c_str());
                tree_carbon->Draw(Form("pT>>h_pT_c_%d%d%d", k, j, i), cuts.c_str());
                tree_ch->Draw(Form("pT>>h_pT_ch_%d%d%d", k, j, i), cuts.c_str());
                tree_he->Draw(Form("pT>>h_pT_he_%d%d%d", k, j, i), cuts.c_str());
                tree_empty->Draw(Form("pT>>h_pT_empty_%d%d%d", k, j, i), cuts.c_str());

                // Here you would include the steps for calculating and fitting the dilution factor,
                // similar to the previous method with all five target types. Since this is just the 
                // minimal working function, I'm skipping that part for now.
                
                // Cleanup the histograms after each iteration
                delete h_pT_nh3;
                delete h_pT_c;
                delete h_pT_ch;
                delete h_pT_he;
                delete h_pT_empty;
            }
        }
        // Save the canvas after all pads are filled
        c1->SaveAs(Form("output/multidimensional_ybin_%d.png", k));
        delete c1;
    }

    // Close the files
    nh3->Close();
    carbon->Close();
    ch->Close();
    he->Close();
    empty->Close();

    return 1;  // Return some meaningful value
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
    // one_dimensional(nh3, c, ch, he, empty);
    multi_dimensional(nh3, c, ch, he, empty);

    // Safely close the ROOT files
    nh3->Close();
    c->Close();
    ch->Close();
    he->Close();
    empty->Close();

    return 0;
}