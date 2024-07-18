#include <iostream>
#include <string>
#include <cmath>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLatex.h>

// Function to create and fill histograms
void createAndFillHistograms(TTree* tree_nh3, TTree* tree_carbon, const char* branch_name,
                             TH1D*& h_nh3, TH1D*& h_carbon, double hist_min, double hist_max) {
    std::string hist_name_nh3 = std::string(branch_name) + "_nh3";
    std::string hist_name_carbon = std::string(branch_name) + "_carbon";
    h_nh3 = new TH1D(hist_name_nh3.c_str(), (std::string(branch_name) + " Distribution; " + branch_name + "; Counts").c_str(), 100, hist_min, hist_max);
    h_carbon = new TH1D(hist_name_carbon.c_str(), (std::string(branch_name) + " Distribution; " + branch_name + "; Counts").c_str(), 100, hist_min, hist_max);

    tree_nh3->Draw(Form("%s>>%s", branch_name, hist_name_nh3.c_str()));
    tree_carbon->Draw(Form("%s>>%s", branch_name, hist_name_carbon.c_str()));
}

// Function to plot histograms and ratios, and return the fitted constant for Mx
double plotHistogramsAndRatios(TTree* tree_nh3, TTree* tree_carbon) {
    // Create canvas and divide it into six panels
    TCanvas* c1 = new TCanvas("c1", "Dilution Factor Analysis", 1200, 1200);
    c1->Divide(3, 2);

    double fitted_const_Mx = 0;

    // Define the variables and their properties
    const char* branches[] = {"Mx", "xF", "zeta"};
    const char* titles[] = {"M_{x} (GeV)", "x_{F}", "#zeta"};
    double hist_mins[] = {-2, -2, 0};
    double hist_maxs[] = {3, 1, 1.5};
    double fit_mins[] = {-2, -2, 1.0};
    double fit_maxs[] = {-0.5, -1.25, 1.5};
    double plot_y_mins[] = {9, 5, 5};
    double plot_y_maxs[] = {15, 15, 15};

    for (int i = 0; i < 3; ++i) {
        TH1D* h_nh3 = nullptr;
        TH1D* h_carbon = nullptr;
        createAndFillHistograms(tree_nh3, tree_carbon, branches[i], h_nh3, h_carbon, hist_mins[i], hist_maxs[i]);

        // Plot histograms
        c1->cd(2 * i + 1);
        gPad->SetLeftMargin(0.15);
        gPad->SetLogy();
        h_nh3->SetLineColor(kBlue);
        h_carbon->SetLineColor(kRed);
        h_nh3->Draw();
        h_carbon->Draw("SAME");

        TLegend* leg = new TLegend(0.75, 0.8, 0.9, 0.9);
        leg->AddEntry(h_nh3, "NH_{3}", "l");
        leg->AddEntry(h_carbon, "C", "l");
        leg->Draw();

        h_nh3->SetStats(0);
        h_carbon->SetStats(0);

        // Plot ratios
        c1->cd(2 * i + 2);
        gPad->SetLeftMargin(0.15);
        TGraphErrors* gr_ratio = new TGraphErrors();
        for (int j = 1; j <= h_nh3->GetNbinsX(); ++j) {
            double nh3_counts = h_nh3->GetBinContent(j);
            double c_counts = h_carbon->GetBinContent(j);
            if (c_counts > 0) {
                double ratio = nh3_counts / c_counts;
                double error = ratio * std::sqrt(1 / nh3_counts + 1 / c_counts);
                gr_ratio->SetPoint(j - 1, h_nh3->GetBinCenter(j), ratio);
                gr_ratio->SetPointError(j - 1, 0, error);
            }
        }
        gr_ratio->GetYaxis()->SetRangeUser(plot_y_mins[i], plot_y_maxs[i]);
        gr_ratio->SetTitle((std::string(titles[i]) + " Ratio; " + titles[i] + "; Ratio").c_str());
        gr_ratio->SetMarkerStyle(20);
        gr_ratio->Draw("AP");

        // Fit the data to a constant
        TF1* fit_const = new TF1("fit_const", "[0]", fit_mins[i], fit_maxs[i]);
        gr_ratio->Fit(fit_const, "R");
        fit_const->SetLineColor(kRed);
        fit_const->Draw("SAME");

        // Add the dotted-dashed line for the remaining range
        TF1* dotted_line = new TF1("dotted_line", "[0]", fit_maxs[i], hist_maxs[i]);
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

        if (i == 0) {
            fitted_const_Mx = fit_value;
        }

        delete h_nh3;
        delete h_carbon;
        delete gr_ratio;
        delete fit_const;
        delete dotted_line;
        delete leg;
    }

    // Save the canvas
    c1->SaveAs("dilution_factors.pdf");

    // Clean up
    delete c1;

    return fitted_const_Mx;
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

    // Plot histograms and ratios, and get the fitted constant for Mx
    double fitted_const_Mx = plotHistogramsAndRatios(tree_nh3, tree_carbon);

    // Output the fitted constant for Mx
    std::cout << "Fitted constant for Mx: " << fitted_const_Mx << std::endl;

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