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
#include <vector>
#include <TLegend.h>

// Constants for the dilution factor calculation
const double L_C = 1.5;
const double L = 5.86;
const double L_CH = 3;
const double rho_He = 0.145 / 4;
const double rho_C = 2.0 / 12;
const double rho_CH = 1.0 / 14;
const double rho_A = 0.92 / 17;

// Fractional charge values
const double xA = 0.71297;
const double xC = 0.08399;
const double xCH = 0.03600;
const double xHe = 0.07535;
const double xf = 0.09168;

// Total accumulated charge
const double nc_A = 3599616.7128999997; 
const double nc_C = 424018.35203;
const double nc_CH = 181775.30390000003;
const double nc_He = 380446.68830000004;
const double nc_ET = 462880.79;

// Declare vectors to store the dynamically allocated objects
std::vector<TGraphErrors*> dilution_graphs;
std::vector<TF1*> fit_functions;

// Function to create and normalize histograms
TH1D* create_normalized_histogram(TTree* tree, const char* branch, const char* hist_name, double total_charge, double min_range, double max_range) {
    TH1D* hist = new TH1D(hist_name, "", 100, min_range, max_range);
    tree->Draw(Form("%s>>%s", branch, hist_name));
    hist->Scale(1.0 / total_charge);
    return hist;
}

// Function to create a TGraphErrors from a histogram
TGraphErrors* create_tgrapherrors(TH1D* hist, int color, int marker_style) {
    int n_bins = hist->GetNbinsX();
    TGraphErrors* graph = new TGraphErrors(n_bins);
    graph->SetMarkerColor(color);
    graph->SetMarkerStyle(marker_style);
    for (int i = 1; i <= n_bins; ++i) {
        double x = hist->GetBinCenter(i);
        double y = hist->GetBinContent(i);
        double ex = 0;
        double ey = hist->GetBinError(i);
        graph->SetPoint(i - 1, x, y);
        graph->SetPointError(i - 1, ex, ey);
    }
    return graph;
}

void plot_dilution_kinematics(TFile* nh3, TFile* carbon, TFile* ch, TFile* he, TFile* empty) {
    // Branches to plot
    std::vector<std::string> branches = {"e_p", "e_theta*180.0/3.14159", "p_p", "p_theta*180.0/3.14159", "Q2", "x", "pT", "xF"};
    std::vector<std::string> x_labels = {
        "e_{p} (GeV)", "e_{#theta} (degrees)", "p_{p} (GeV)", "p_{#theta} (degrees)", 
        "Q^{2} (GeV^{2})", "x_{B}", "P_{T} (GeV)", "x_{F}"
    };
    std::vector<std::pair<double, double>> ranges = {
        {2, 8}, {0, 45}, {0, 4}, {0, 60}, {1, 8}, {0, 0.6}, {0, 1.2}, {-1, 1}
    };

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
        return;
    }

    // Create a 2x4 canvas
    TCanvas* c1 = new TCanvas("c1", "Dilution Kinematics", 1600, 800);
    c1->Divide(4, 2);

    // Loop over the branches and plot them
    for (size_t i = 0; i < branches.size(); ++i) {
        c1->cd(i + 1)->SetLeftMargin(0.15); // Add padding to the left

        // Create histograms for each target type
        TH1D* h_nh3 = create_normalized_histogram(tree_nh3, branches[i].c_str(), Form("h_nh3_%zu", i), nc_A, ranges[i].first, ranges[i].second);
        TH1D* h_carbon = create_normalized_histogram(tree_carbon, branches[i].c_str(), Form("h_carbon_%zu", i), nc_C, ranges[i].first, ranges[i].second);
        TH1D* h_ch = create_normalized_histogram(tree_ch, branches[i].c_str(), Form("h_ch_%zu", i), nc_CH, ranges[i].first, ranges[i].second);
        TH1D* h_he = create_normalized_histogram(tree_he, branches[i].c_str(), Form("h_he_%zu", i), nc_He, ranges[i].first, ranges[i].second);
        TH1D* h_empty = create_normalized_histogram(tree_empty, branches[i].c_str(), Form("h_empty_%zu", i), nc_ET, ranges[i].first, ranges[i].second);

        // Create TGraphErrors for each histogram
        TGraphErrors* gr_nh3 = create_tgrapherrors(h_nh3, kRed, 20);
        TGraphErrors* gr_carbon = create_tgrapherrors(h_carbon, kBlue, 21);
        TGraphErrors* gr_ch = create_tgrapherrors(h_ch, kGreen, 22);
        TGraphErrors* gr_he = create_tgrapherrors(h_he, kMagenta, 23);
        TGraphErrors* gr_empty = create_tgrapherrors(h_empty, kCyan, 24);

        // Determine the y-axis range
        double max_y = std::max({h_nh3->GetMaximum(), h_carbon->GetMaximum(), h_ch->GetMaximum(), h_he->GetMaximum(), h_empty->GetMaximum()});
        gr_nh3->GetYaxis()->SetRangeUser(0, max_y * 1.1);

        // Draw the graphs
        gr_nh3->Draw("AP");
        gr_carbon->Draw("P SAME");
        gr_ch->Draw("P SAME");
        gr_he->Draw("P SAME");
        gr_empty->Draw("P SAME");

        // Set axis labels
        gr_nh3->GetXaxis()->SetTitle(x_labels[i].c_str());
        gr_nh3->GetYaxis()->SetTitle("Counts/nC");

        // Remove the default title
        gr_nh3->SetTitle("");

        // Create legend
        TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(gr_nh3, "NH3", "P");
        legend->AddEntry(gr_carbon, "C", "P");
        legend->AddEntry(gr_ch, "CH", "P");
        legend->AddEntry(gr_he, "He", "P");
        legend->AddEntry(gr_empty, "Empty", "P");
        legend->Draw();

        // Clean up histograms
        delete h_nh3;
        delete h_carbon;
        delete h_ch;
        delete h_he;
        delete h_empty;
    }

    // Add a global title to the canvas and center it
    c1->cd();
    TLatex title;
    title.SetNDC();
    title.SetTextSize(0.04);
    title.SetTextAlign(22);  // Center alignment
    title.DrawLatex(0.5, 0.97, "epX, Q^{2} > 1 GeV^{2}, W > 2, y < 0.75");

    // Save the canvas
    c1->SaveAs("output/dilution_kinematics.png");

    // Clean up
    delete c1;
}

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
                          TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad, bool skip_fit = false, bool isMx = false) {
    canvas->cd(pad);
    gPad->SetLeftMargin(0.15);

    // Define the base cuts for vz
    // std::string vz_cuts = "-10 < vz_e && vz_e < 1 && -10 < vz_p && vz_p < 1 && runnum != 16317 && runnum != 16742";
    std::string vz_cuts = "x > 0.06 && x < 0.60 && pT > 0 && pT < 1.2 && xF > -1 && xF < 1 && -10 < vz_e && vz_e < 1 && -10 < vz_p && vz_p < 1 && runnum != 16317 && runnum != 16742";
    // std::string vz_cuts = "x > 0.14 && x < 0.21 && pT > 0 && pT < 1.2 && xF > -1 && xF < 1 && -10 < vz_e && vz_e < 1 && -10 < vz_p && vz_p < 1 && runnum != 16317 && runnum != 16742";
    // std::string vz_cuts = "z > 0.19 && z < 0.30 && pT > 0 && pT < 1.2 && xF > -1 && xF < 1 && -10 < vz_e && vz_e < 1 && -10 < vz_p && vz_p < 1 && runnum != 16317 && runnum != 16742";


    // Define the combined cuts based on the value of isMx
    std::string combined_cuts;
    if (isMx) {
        // Apply vz cuts and Mx > 0 if isMx is true
        combined_cuts = "Mx > 0 && " + vz_cuts;
    } else {
        // Apply both Mx > 1.35 and vz cuts if isMx is false
        combined_cuts = "Mx > 1.35 && " + vz_cuts;
        // combined_cuts = "Mx > 0.55 && " + vz_cuts;
        // combined_cuts = "Mx < 1.35 && Mx > 0 && " + vz_cuts;
    }

    // Define the combined cuts based on the value of isMx
    std::string combined_cuts_all;
    combined_cuts_all = "Mx> 0.55 && " + vz_cuts;

    // Create histograms for data using the appropriate cuts
    TH1D *h_nh3 = new TH1D(Form("h_%s_nh3", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_c = new TH1D(Form("h_%s_c", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_ch = new TH1D(Form("h_%s_ch", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_he = new TH1D(Form("h_%s_he", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_empty = new TH1D(Form("h_%s_empty", variable_name), "", n_bins, x_min, x_max);

    // Create histograms for data using the appropriate cuts
    TH1D *h_nh3_all = new TH1D(Form("h_%s_nh3_all", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_c_all = new TH1D(Form("h_%s_c_all", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_ch_all = new TH1D(Form("h_%s_ch_all", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_he_all = new TH1D(Form("h_%s_he_all", variable_name), "", n_bins, x_min, x_max);
    TH1D *h_empty_all = new TH1D(Form("h_%s_empty_all", variable_name), "", n_bins, x_min, x_max);

    // Draw the histograms with the appropriate cuts
    nh3->Draw(Form("%s>>h_%s_nh3", variable_name, variable_name), combined_cuts.c_str());
    c->Draw(Form("%s>>h_%s_c", variable_name, variable_name), combined_cuts.c_str());
    ch->Draw(Form("%s>>h_%s_ch", variable_name, variable_name), combined_cuts.c_str());
    he->Draw(Form("%s>>h_%s_he", variable_name, variable_name), combined_cuts.c_str());
    empty->Draw(Form("%s>>h_%s_empty", variable_name, variable_name), combined_cuts.c_str());

    // Draw the histograms with the appropriate cuts
    nh3->Draw(Form("%s>>h_%s_nh3_all", variable_name, variable_name), combined_cuts_all.c_str());
    c->Draw(Form("%s>>h_%s_c_all", variable_name, variable_name), combined_cuts_all.c_str());
    ch->Draw(Form("%s>>h_%s_ch_all", variable_name, variable_name), combined_cuts_all.c_str());
    he->Draw(Form("%s>>h_%s_he_all", variable_name, variable_name), combined_cuts_all.c_str());
    empty->Draw(Form("%s>>h_%s_empty_all", variable_name, variable_name), combined_cuts_all.c_str());

    // Calculate dilution factor and its error
    TGraphErrors *gr_dilution = new TGraphErrors();
    TGraphErrors *gr_dilution_all = new TGraphErrors();
    for (int i = 1; i <= n_bins; ++i) {
        double nA = h_nh3->GetBinContent(i);
        double nC = h_c->GetBinContent(i);
        double nCH = h_ch->GetBinContent(i);
        double nMT = h_he->GetBinContent(i);
        double nf = h_empty->GetBinContent(i);

        double nA_all = h_nh3_all->GetBinContent(i);
        double nC_all = h_c_all->GetBinContent(i);
        double nCH_all = h_ch_all->GetBinContent(i);
        double nMT_all = h_he_all->GetBinContent(i);
        double nf_all = h_empty_all->GetBinContent(i);

        double dilution = calculate_dilution_factor(nA, nC, nCH, nMT, nf);
        double error = calculate_dilution_error(nA/xA, nC/xC, nCH/xCH, nMT/xHe, nf/xf);

        double dilution_all = calculate_dilution_factor(nA_all, nC_all, nCH_all, nMT_all, nf_all);
        double error_all = calculate_dilution_error(nA_all/xA, nC_all/xC, nCH_all/xCH, nMT_all/xHe, nf_all/xf);

        // For integrated plot, set the point at the center of the plot range
        double x_position = skip_fit ? (x_min + x_max) / 2 : h_nh3->GetBinCenter(i);

        gr_dilution->SetPoint(i - 1, x_position, dilution);
        gr_dilution->SetPointError(i - 1, 0, error);

        gr_dilution_all->SetPoint(i - 1, x_position, dilution_all);
        gr_dilution_all->SetPointError(i - 1, 0, error_all);

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
    gr_dilution->GetYaxis()->SetRangeUser(0.10, 0.40);

    double chi2_scale_factor = 1.0;
    double chi2_scale_factor_all = 1.0;
    double chi2_scale_factor_exclusive = 1.0;

    // Fit and plot (skip fit for the integrated version)
    if (!skip_fit) {
        TF1 *fit_func;
        TF1 *fit_func_all;
        TF1 *fit_func_exclusive;
        if (isMx) {
            // Two Gaussians + Quadratic Polynomial Background
            fit_func = new TF1("fit_func",
                "[0]*exp(-0.5*((x-[1])/[2])^2) + "  // Gaussian 1
                "[3]*exp(-0.5*((x-[4])/[5])^2) + "  // Gaussian 2
                "[6] + [7]*x + [8]*x^2 + [9]*x^3",            // Quadratic Polynomial
                x_min, x_max);

            // Initial guesses
            fit_func->SetParameters(0.05, 0.135, 0.02, 0.5, 0.770, 0.1, 0.1, 0.2, 0.0, 0.0);

            // // Set parameter limits for Gaussians
            fit_func->SetParLimits(0, 0.0, 0.25); // Amplitude 1 must be positive
            fit_func->SetParLimits(1, 0.135 - 0.015, 0.135 + 0.015); // pi0 mass limits in GeV
            // fit_func->SetParLimits(2, 0, 0.3); // pi0 sigma limits in GeV

            fit_func->SetParLimits(3, 0.0, 0.25); // Amplitude 2 must be positive
            fit_func->SetParLimits(4, 0.770 - 0.015, 0.770 + 0.015); // rho0 mass limits in GeV
            // fit_func->SetParLimits(5, 0, 0.15); // rho0 sigma limits in GeV

            // The coefficients of the polynomial are left unconstrained for now
        } else {
            // Use a cubic polynomial fit for other variables
            fit_func = new TF1("fit_func", "[0] + [1]*x + [2]*x^2 + [3]*x^3", x_min, x_max);
            fit_func_all = new TF1("fit_func_all", "[0] + [1]*x + [2]*x^2 + [3]*x^3", x_min, x_max);
            fit_func_exclusive = new TF1("fit_func_exclusive", "[0] + [1]*x + [2]*x^2 + [3]*x^3", x_min, x_max);
        }

        gr_dilution->Fit(fit_func, "RQ");
        fit_func->SetLineColor(kBlack);
        fit_func->SetLineStyle(2); // Dashed line

        if (!isMx) {
            gr_dilution_all->Fit(fit_func_all, "RQ");
            fit_func_all->SetLineColor(kRed);
            fit_func_all->SetLineStyle(2); // Dashed line
            // fit_func_all->SetLineWidth(1); // Set thinner line
        }

        // Calculate chi2/ndf scaling factor
        double chi2 = fit_func->GetChisquare();
        int ndf = fit_func->GetNDF();
        chi2_scale_factor = std::sqrt(chi2 / ndf);

        double chi2_all; int ndf_all;

        if (!isMx) {
            // Calculate chi2/ndf scaling factor
            chi2_all = fit_func_all->GetChisquare();
            ndf_all = fit_func_all->GetNDF();
            chi2_scale_factor_all = std::sqrt(chi2_all / ndf_all);
        }
        
        // Rescale the errors
        for (int i = 0; i < gr_dilution->GetN(); ++i) {
            double x, y;
            gr_dilution->GetPoint(i, x, y);
            gr_dilution->SetPointError(i, 0, gr_dilution->GetErrorY(i) * chi2_scale_factor);
        }

        if (!isMx) {
            for (int i = 0; i < gr_dilution_all->GetN(); ++i) {
                double x, y;
                gr_dilution_all->GetPoint(i, x, y);
                gr_dilution_all->SetPointError(i, 0, gr_dilution_all->GetErrorY(i) * chi2_scale_factor_all);
            }
        }

        // Refit with scaled errors
        gr_dilution->Fit(fit_func, "RQ");
        fit_func->Draw("SAME");

        if (!isMx) {
            gr_dilution_all->Fit(fit_func, "RQ");
            fit_func_all->Draw("SAME");

            // gr_dilution_exclusive->Fit(fit_func, "RQ");
            // fit_func_exclusive->Draw("SAME");
        }

        // Commented out chi2/ndf printing
        /*
        // Print chi2/ndf information
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035); // Decrease the font size
        latex.DrawLatex(0.20, 0.15, Form("#chi^{2}/NDF = %.2f / %d = %.2f", chi2, ndf, chi2 / ndf));
        */

        if (!isMx) {
            // Add red note "0 < M_{x} < 1.35 (GeV)"
            // TLatex latex_red;
            // latex_red.SetNDC();
            // latex_red.SetTextSize(0.035); // Adjust the font size if needed
            // latex_red.SetTextColor(kRed); // Set the text color to red
            // latex_red.DrawLatex(0.2, 0.85, "0.55 < M_{x} (GeV) < 1.35");

            // Add blue note "0 < M_{x} (GeV)" just below the red one
            TLatex latex_blue;
            latex_blue.SetNDC();
            latex_blue.SetTextSize(0.035); // Adjust the font size if needed
            latex_blue.SetTextColor(kRed); // Set the text color to blue
            latex_blue.DrawLatex(0.2, 0.85, "M_{x} > 0.55 GeV");
        }

        // Add fit parameters box
        double box_x1 = (isMx) ? 0.45 : 0.55;
        double box_y1 = (isMx) ? 0.50 : 0.7; // Slightly lower start position for Mx plot
        double box_y2 = (isMx) ? 0.9 : 0.9; // Increase vertical size more for Mx plot, but within plot limits
        TPaveText *pt = new TPaveText(box_x1, box_y1, 0.9, box_y2, "brNDC");
        pt->SetBorderSize(1);
        pt->SetFillStyle(1001);
        pt->SetFillColor(kWhite);
        pt->SetTextSize(0.035); // Decrease the font size

        // Add text for the parameters based on the fit function used
        if (isMx) {
            for (int p = 0; p < 6; p += 3) { // Adjust the loop to iterate over two Gaussians
                pt->AddText(Form("Amp%d = %.3f +/- %.3f", p/3+1, fit_func->GetParameter(p), fit_func->GetParError(p)));
                
                if (p == 0) {
                    pt->AddText(Form("#mu_{1} = %.3f +/- %.3f", fit_func->GetParameter(p+1), fit_func->GetParError(p+1)));
                    pt->AddText(Form("#sigma_{1} = %.3f +/- %.3f", fit_func->GetParameter(p+2), fit_func->GetParError(p+2)));
                } else if (p == 3) {
                    pt->AddText(Form("#rho^{0} mass (GeV) = %.3f +/- %.3f", fit_func->GetParameter(p+1), fit_func->GetParError(p+1)));
                    pt->AddText(Form("#sigma#rho^{0} mass (GeV) = %.3f +/- %.3f", fit_func->GetParameter(p+2), fit_func->GetParError(p+2)));
                }
            }
            
            // Add the polynomial coefficients
            pt->AddText(Form("Const = %.3f +/- %.3f", fit_func->GetParameter(6), fit_func->GetParError(6)));
            pt->AddText(Form("Linear = %.3f +/- %.3f", fit_func->GetParameter(7), fit_func->GetParError(7)));
            pt->AddText(Form("Quadratic = %.3f +/- %.3f", fit_func->GetParameter(8), fit_func->GetParError(8)));
        } else {
            for (int p = 0; p < fit_func->GetNpar(); ++p) {
                pt->AddText(Form("p%d = %.3f +/- %.3f", p, fit_func->GetParameter(p), fit_func->GetParError(p)));
            }
        }
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
    // Add a title to the plot with phase space parameters
    TLatex title;
    title.SetNDC();
    if (isMx) {
        title.SetTextSize(0.045);
        title.DrawLatex(0.15, 0.92, "Q^{2} > 1.0 GeV^{2}, W > 2 GeV, y < 0.75, M_{x} > 0 GeV");
        title.SetTextSize(0.035);
        // title.DrawLatex(0.15, 0.92, "Q^{2} > 1.0 GeV^{2}, W > 2 GeV, y < 0.75, M_{x} > 0 GeV, 0.14 < x_{B} < 0.21");
        // title.DrawLatex(0.15, 0.92, "Q^{2} > 1.0 GeV^{2}, W > 2 GeV, y < 0.75, M_{x} > 0 GeV, 0.19 < z < 0.30");
    } else {
        title.SetTextSize(0.045);
        title.DrawLatex(0.15, 0.92, "Q^{2} > 1.0 GeV^{2}, W > 2 GeV, y < 0.75, M_{x} > 1.35 GeV");
        // title.SetTextSize(0.035);
        // title.DrawLatex(0.15, 0.92, "Q^{2} > 1.0 GeV^{2}, W > 2 GeV, y < 0.75, M_{x} > 1.35 GeV, 0.14 < x_{B} < 0.21");
        // title.DrawLatex(0.15, 0.92, "Q^{2} > 1.0 GeV^{2}, W > 2 GeV, y < 0.75, M_{x} > 1.35 GeV, 0.19 < z < 0.30");
    }

    // Clean up histograms
    delete h_nh3;
    delete h_c;
    delete h_ch;
    delete h_he;
    delete h_empty;

    delete h_nh3_all;
    delete h_c_all;
    delete h_ch_all;
    delete h_he_all;
    delete h_empty_all;
}

std::pair<TF1*, TGraphErrors*> fit_and_plot_dilution(const char* variable_name, const char* x_title, double x_min, double x_max, int n_bins,
TTree* nh3, TTree* c, TTree* ch, TTree* he, TTree* empty, TCanvas* canvas, int pad, bool skip_fit = false, bool isMx = false) {
    // Call the plotting function
    plot_dilution_factor(variable_name, x_title, x_min, x_max, n_bins, nh3, c, ch, he, empty, canvas, pad, skip_fit, isMx);
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
    // Create a canvas and divide it into 3 rows and 3 columns
    TCanvas *c1 = new TCanvas("c1", "Dilution Factor Analysis", 1600, 1200);
    c1->Divide(3, 3);

    // Prepare to print the fit functions for each variable
    std::cout << std::endl << std::endl;

    // Integrated version (single bin)
    auto fit_integrated = fit_and_plot_dilution("x", "", 0.0, 1.0, 1, nh3, c, ch, he, empty, c1, 1, true, false);

    // Fit and plot for Q2
    auto fit_Q2 = fit_and_plot_dilution("Q2", "Q^{2} (GeV^{2})", 1, 9, 25, nh3, c, ch, he, empty, c1, 2, false, false);
    if (fit_Q2.first) {
        double p0_x = fit_Q2.first->GetParameter(0);
        double p1_x = fit_Q2.first->GetParameter(1);
        double p2_x = fit_Q2.first->GetParameter(2);
        std::cout << "if (prefix == \"Q2\") { return " << p0_x << 
            "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for x-Bjorken
    auto fit_x = fit_and_plot_dilution("x", "x_{B} (GeV)", 0.06, 0.6, 25, nh3, c, ch, he, empty, c1, 3, false, false);
    if (fit_x.first) {
        double p0_x = fit_x.first->GetParameter(0);
        double p1_x = fit_x.first->GetParameter(1);
        double p2_x = fit_x.first->GetParameter(2);
        std::cout << "if (prefix == \"x\") { return " << p0_x << 
            "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for y
    auto fit_y = fit_and_plot_dilution("y", "y", 0.3, 0.75, 25, nh3, c, ch, he, empty, c1, 4, false, false);
    if (fit_y.first) {
        double p0_x = fit_y.first->GetParameter(0);
        double p1_x = fit_y.first->GetParameter(1);
        double p2_x = fit_y.first->GetParameter(2);
        std::cout << "if (prefix == \"y\") { return " << p0_x << 
            "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for z
    auto fit_z = fit_and_plot_dilution("z", "z", 0.06, 0.8, 25, nh3, c, ch, he, empty, c1, 5, false, false);
    if (fit_z.first) {
        double p0_x = fit_z.first->GetParameter(0);
        double p1_x = fit_z.first->GetParameter(1);
        double p2_x = fit_z.first->GetParameter(2);
        std::cout << "if (prefix == \"z\") { return " << p0_x << 
            "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for zeta
    auto fit_zeta = fit_and_plot_dilution("zeta", "#zeta", 0.3, 0.7, 25, nh3, c, ch, he, empty, c1, 6, false, false);
    if (fit_zeta.first) {
        double p0_x = fit_zeta.first->GetParameter(0);
        double p1_x = fit_zeta.first->GetParameter(1);
        double p2_x = fit_zeta.first->GetParameter(2);
        std::cout << "if (prefix == \"zeta\") { return " << p0_x << 
            "+" << p1_x << "*currentVariable+" << p2_x << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for transverse momentum
    auto fit_pT = fit_and_plot_dilution("pT", "P_{T} (GeV)", 0, 1.0, 25, nh3, c, ch, he, empty, c1, 7, false, false);
    if (fit_pT.first) {
        double p0_PT = fit_pT.first->GetParameter(0);
        double p1_PT = fit_pT.first->GetParameter(1);
        double p2_PT = fit_pT.first->GetParameter(2);
        std::cout << "if (prefix == \"PT\") { return " << p0_PT << 
            "+" << p1_PT << "*currentVariable+" << p2_PT << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for x-Feynman
    auto fit_xF = fit_and_plot_dilution("xF", "x_{F}", -0.8, 0.5, 25, nh3, c, ch, he, empty, c1, 8, false, false);
    if (fit_xF.first) {
        double p0_xF = fit_xF.first->GetParameter(0);
        double p1_xF = fit_xF.first->GetParameter(1);
        double p2_xF = fit_xF.first->GetParameter(2);
        std::cout << "if (prefix == \"xF\") { return " << p0_xF <<
        "+" << p1_xF << "*currentVariable+" << p2_xF << "*std::pow(currentVariable,2); }" << std::endl;
    }

    // Fit and plot for Mx
    auto fit_Mx = fit_and_plot_dilution("Mx", "M_{x} (GeV)", 0 , 2.75, 50, nh3, c, ch, he, empty, c1, 9, false, true);
    if (fit_Mx.first) {
        double amp1 = fit_Mx.first->GetParameter(0);
        double mean1 = fit_Mx.first->GetParameter(1);
        double sigma1 = fit_Mx.first->GetParameter(2);

        double amp2 = fit_Mx.first->GetParameter(3);
        double mean2 = fit_Mx.first->GetParameter(4);
        double sigma2 = fit_Mx.first->GetParameter(5);

        double constTerm = fit_Mx.first->GetParameter(6);
        double linearTerm = fit_Mx.first->GetParameter(7);
        double quadTerm = fit_Mx.first->GetParameter(8);

        std::cout << "if (prefix == \"Mx\") {"
                  << " return " << amp1 << "*exp(-0.5*std::pow((currentVariable - " << mean1 
                  << ") / " << sigma1 << ", 2)) + "
                  << amp2 << "*exp(-0.5*std::pow((currentVariable - " << mean2 
                  << ") / " << sigma2 << ", 2)) + "
                  << constTerm << " + "
                  << linearTerm << "*currentVariable + "
                  << quadTerm << "*std::pow(currentVariable, 2); }"
                  << std::endl;
    }

    // Save the canvas as a PNG file
    c1->SaveAs("output/one_dimensional.png");

    // Clean up
    delete c1;
}

std::vector<TH1D*> create_and_draw_histograms(TTree* tree_nh3, TTree* tree_carbon, TTree* tree_ch, TTree* tree_he, TTree* tree_empty, const std::string& cuts, int k, int i) {
    // Create histograms for different targets
    TH1D *h_pT_nh3 = new TH1D(Form("h_pT_nh3_%d%d", k, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
    TH1D *h_pT_c = new TH1D(Form("h_pT_c_%d%d", k, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
    TH1D *h_pT_ch = new TH1D(Form("h_pT_ch_%d%d", k, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
    TH1D *h_pT_he = new TH1D(Form("h_pT_he_%d%d", k, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);
    TH1D *h_pT_empty = new TH1D(Form("h_pT_empty_%d%d", k, i), "P_{T} Distribution; P_{T} (GeV); Counts", 9, 0, 1.0);

    // Define the additional cuts
    std::string additional_cuts = "x > 0.06 && x < 0.60 && pT > 0 && pT < 1.2 && xF > -1 && xF < 1 && Mx > 1.35 && -10 < vz_e && vz_e < 1 && -10 < vz_p && vz_p < 1 && runnum != 16317 && runnum != 16742";
    // std::string additional_cuts = "x > 0.06 && x < 0.60 && pT > 0 && pT < 1.2 && xF > -1 && xF < 1 && Mx > 0.55 && -10 < vz_e && vz_e < 1 && -10 < vz_p && vz_p < 1 && runnum != 16317 && runnum != 16742";

    // Combine the existing cuts with the additional cuts
    std::string combined_cuts = cuts + " && " + additional_cuts;

    // Draw histograms with the combined cuts
    tree_nh3->Draw(Form("pT>>h_pT_nh3_%d%d", k, i), combined_cuts.c_str());
    tree_carbon->Draw(Form("pT>>h_pT_c_%d%d", k, i), combined_cuts.c_str());
    tree_ch->Draw(Form("pT>>h_pT_ch_%d%d", k, i), combined_cuts.c_str());
    tree_he->Draw(Form("pT>>h_pT_he_%d%d", k, i), combined_cuts.c_str());
    tree_empty->Draw(Form("pT>>h_pT_empty_%d%d", k, i), combined_cuts.c_str());

    // Store histograms in a vector
    std::vector<TH1D*> histograms = {h_pT_nh3, h_pT_c, h_pT_ch, h_pT_he, h_pT_empty};
    return histograms;
}

double multi_dimensional(TFile* nh3, TFile* carbon, TFile* ch, TFile* he, TFile* empty) {
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
        return 0;
    }

    // Q2-x bin information
    struct BinInfo {
        std::string x_title;
        std::string x_range;
        std::string q2_range;
    };

    std::vector<BinInfo> bins = {
        {"x < 0.1", "x < 0.1", "Q2 > 0"},
        {"0.1 < x < 0.14", "x > 0.1 && x < 0.14 && Q2 < 1.50", "Q2 < 1.50"},
        {"0.1 < x < 0.14", "x > 0.1 && x < 0.14 && Q2 >= 1.50 && Q2 < 1.70", "1.50 <= Q2 < 1.70"},
        {"0.1 < x < 0.14", "x > 0.1 && x < 0.14 && Q2 >= 1.70", "Q2 >= 1.70"},
        {"0.14 < x < 0.21", "x > 0.14 && x < 0.21 && Q2 < 1.50", "Q2 < 1.50"},
        {"0.14 < x < 0.21", "x > 0.14 && x < 0.21 && Q2 >= 1.50 && Q2 < 1.70", "1.50 <= Q2 < 1.70"},
        {"0.14 < x < 0.21", "x > 0.14 && x < 0.21 && Q2 >= 1.70 && Q2 < 2.00", "1.70 <= Q2 < 2.00"},
        {"0.14 < x < 0.21", "x > 0.14 && x < 0.21 && Q2 >= 2.00", "Q2 >= 2.00"},
        {"0.21 < x < 0.30", "x > 0.21 && x < 0.30 && Q2 < 2.20", "Q2 < 2.20"},
        {"0.21 < x < 0.30", "x > 0.21 && x < 0.30 && Q2 >= 2.20 && Q2 < 2.60", "2.20 <= Q2 < 2.60"},
        {"0.21 < x < 0.30", "x > 0.21 && x < 0.30 && Q2 >= 2.60", "Q2 >= 2.60"},
        {"0.30 < x < 0.42", "x > 0.30 && x < 0.42 && Q2 < 3.20", "Q2 < 3.20"},
        {"0.30 < x < 0.42", "x > 0.30 && x < 0.42 && Q2 >= 3.20", "Q2 >= 3.20"},
        {"x > 0.42", "x >= 0.42", "Q2 > 0"}
    };

    // Define z bins
    std::vector<std::string> z_ranges = {
        "0 < z && z <= 0.19",
        "0.19 < z && z <= 0.30",
        "0.30 < z && z <= 0.42",
        "0.42 < z && z <= 1.00"
    };

    std::vector<std::string> z_titles = {
        "0 < z <= 0.19",
        "0.19 < z <= 0.30",
        "0.30 < z <= 0.42",
        "0.42 < z <= 1.00"
    };

    std::vector<std::string> z_prefixes = {
        "z1", "z2", "z3", "z4"
    };

    // Create canvases based on grouped x bins
    std::vector<TCanvas*> canvases;
    int padIndex = 1;

    for (size_t k = 0; k < bins.size(); ++k) {
        if (k == 0 || k == 1 || k == 4 || k == 8 || k == 11 || k == 13) {
            // Determine the number of rows needed
            int nRows = 1;
            if (k == 4) nRows = 4;
            else if (k == 1 || k == 8) nRows = 3;
            else if (k == 11) nRows = 2;

            // Create a new canvas
            TCanvas* c1 = new TCanvas(Form("c1_%zu", canvases.size()), "Dilution Factor Analysis", 1600, 400 * nRows);
            c1->Divide(4, nRows);  // 4 columns, dynamic rows
            canvases.push_back(c1);
            padIndex = 1;
        }

        // Loop over z bins
        for (size_t i = 0; i < z_ranges.size(); ++i) {
            std::string cuts = bins[k].x_range + " && " + z_ranges[i];
            canvases.back()->cd(padIndex++);
            gPad->SetLeftMargin(0.15);

            // Call the function to create and draw histograms
            std::vector<TH1D*> histograms = create_and_draw_histograms(tree_nh3, tree_carbon, tree_ch, tree_he, tree_empty, cuts, k, i);
            
            // Access the histograms using the vector
            TH1D *h_pT_nh3 = histograms[0];
            TH1D *h_pT_c = histograms[1];
            TH1D *h_pT_ch = histograms[2];
            TH1D *h_pT_he = histograms[3];
            TH1D *h_pT_empty = histograms[4];

            // Inside loop after creating the histograms
            int n_bins = h_pT_nh3->GetNbinsX();
            TGraphErrors *gr_dilution = new TGraphErrors(n_bins);

            for (int bin = 1; bin <= n_bins; ++bin) {
                // Get bin contents for each target type
                double nA = h_pT_nh3->GetBinContent(bin);
                double nC = h_pT_c->GetBinContent(bin);
                double nCH = h_pT_ch->GetBinContent(bin);
                double nMT = h_pT_he->GetBinContent(bin);
                double nf = h_pT_empty->GetBinContent(bin);
                // Calculate the dilution factor
                double dilution = calculate_dilution_factor(nA, nC, nCH, nMT, nf);
                double dilution_error = calculate_dilution_error(nA / xA, nC / xC, nCH / xCH, nMT / xHe, nf / xf);

                // Get the bin center
                double x_position = h_pT_nh3->GetBinCenter(bin);

                // Set the dilution factor point and error in the TGraphErrors
                gr_dilution->SetPoint(bin - 1, x_position, dilution);
                gr_dilution->SetPointError(bin - 1, 0, dilution_error);
            }

            // Use the reformatted strings in the title
            std::string title = Form("Q2x%d (%s), %s, %s", k + 1, bins[k].q2_range.c_str(), bins[k].x_title.c_str(), z_titles[i].c_str());
            gr_dilution->SetTitle((title + "; P_{T} (GeV); D_{f}").c_str());
            gr_dilution->SetMarkerStyle(20);

            // Draw the TGraphErrors on the canvas
            gr_dilution->Draw("AP");
            gr_dilution->GetXaxis()->SetLimits(0, 1);
            gr_dilution->GetYaxis()->SetRangeUser(0.00, 0.50);  // Set the y-axis range from 0.0 to 0.5
            // Increase the size of axis labels and titles
            gr_dilution->GetXaxis()->SetTitleSize(0.05);  // Increase title size
            gr_dilution->GetYaxis()->SetTitleSize(0.05);  // Increase title size
            gr_dilution->GetXaxis()->SetLabelSize(0.04);  // Increase label size
            gr_dilution->GetYaxis()->SetLabelSize(0.04);  // Increase label size
            // Fit the dilution factor to a constant function
            TF1 *fit_func = new TF1("fit_func", "[0]", 0, 1.0);  // Constant fit
            gr_dilution->Fit(fit_func, "RQ");
            fit_func->SetLineColor(kRed);
            fit_func->Draw("SAME");
            // Retrieve fit parameters and chi-squared
            double p0 = fit_func->GetParameter(0);
            double p0_err = fit_func->GetParError(0);
            double chi2 = fit_func->GetChisquare();
            int ndf = fit_func->GetNDF();
            double chi2_ndf = chi2 / ndf;

            // Calculate chi2/ndf scaling factor
            double chi2_scale_factor = std::sqrt(chi2_ndf);

            // Rescale the errors
            for (int bin = 0; bin < gr_dilution->GetN(); ++bin) {
                double x, y;
                gr_dilution->GetPoint(bin, x, y);
                gr_dilution->SetPointError(bin, 0, gr_dilution->GetErrorY(bin) * chi2_scale_factor);
            }

            // Refit with scaled errors
            gr_dilution->Fit(fit_func, "RQ");
            fit_func->Draw("SAME");

            // Retrieve updated chi2 and NDF
            chi2 = fit_func->GetChisquare();
            ndf = fit_func->GetNDF();
            chi2_ndf = chi2 / ndf;

            // Retrieve fit parameters and chi-squared
            p0 = fit_func->GetParameter(0);
            p0_err = fit_func->GetParError(0);

            // Add fit parameters and chi-squared box
            TPaveText *pt = new TPaveText(0.5, 0.7, 0.9, 0.9, "brNDC");
            pt->SetBorderSize(1);
            pt->SetFillStyle(1001);  // Solid fill style
            pt->SetFillColor(kWhite);  // White background
            pt->AddText(Form("p0 = %.3f +/- %.3f", p0, p0_err));
            pt->Draw();

            // Generate the return statement with random Gaussian variation
            std::cout << "if (prefix == \"" << "Q2x" << k + 1 << z_prefixes[i] << "\") {"
                      << " double sigma = " << p0_err << ";"
                      << " return " << p0 << " + rand_gen.Gaus(0, sigma); }" << std::endl << std::endl;

            // Store the objects in vectors for later cleanup
            dilution_graphs.push_back(gr_dilution);
            fit_functions.push_back(fit_func);
        }
    }

    // Save the canvases after all processing is complete
    for (size_t i = 0; i < canvases.size(); ++i) {
        canvases[i]->SaveAs(Form("output/multidimensional_xbin_%zu.png", i));
        delete canvases[i];
    }

    // Clean up the histograms and fit functions
    for (auto graph : dilution_graphs) delete graph;
    for (auto func : fit_functions) delete func;

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

    // Call the plot_dilution_kinematics function
    // plot_dilution_kinematics(nh3, c, ch, he, empty);
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