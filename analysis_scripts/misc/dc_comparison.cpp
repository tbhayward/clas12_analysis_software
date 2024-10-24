#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printf("Usage: ./dc_comparison file1.root file2.root\n");
        return 1;
    }

    const char* file1_name = argv[1];
    const char* file2_name = argv[2];

    // Set style options
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Open the ROOT files
    TFile* file1 = TFile::Open(file1_name);
    TFile* file2 = TFile::Open(file2_name);

    if (!file1 || !file2 || file1->IsZombie() || file2->IsZombie()) {
        printf("Error opening files.\n");
        return 1;
    }

    // Retrieve the trees
    TTree* tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree* tree2 = (TTree*)file2->Get("PhysicsEvents");

    if (!tree1 || !tree2) {
        printf("Error retrieving trees from files.\n");
        return 1;
    }

    // ----------------------------------------
    // Plot 1: Mx2 Comparison (Updated)
    // ----------------------------------------

    // Define histograms for Mx2
    TH1F* h1 = new TH1F("h1", "", 33, 0.6, 1.25);
    TH1F* h2 = new TH1F("h2", "", 33, 0.6, 1.25);

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);

    // Define cuts
    double e_theta_cut = 5.0 * TMath::Pi() / 180.0; // Convert degrees to radians
    TString cuts = Form("fiducial_status==3 && detector1==1 && detector2==1 && e_theta>%f", e_theta_cut);

    // Fill histograms for Mx2
    tree1->Draw("Mx2>>h1", cuts);
    tree2->Draw("Mx2>>h2", cuts);

    // Fit histograms with Gaussian plus linear polynomial
    TF1* fitFunc1 = new TF1("fitFunc1", "gaus(0)+pol1(3)", 0.6, 1.25);
    TF1* fitFunc2 = new TF1("fitFunc2", "gaus(0)+pol1(3)", 0.6, 1.25);

    double proton_mass_sq = 0.938 * 0.938; // Proton mass squared in GeV^2

    // Set initial parameters: [Amplitude, Mean, Sigma, p0, p1]
    fitFunc1->SetParameters(h1->GetMaximum(), proton_mass_sq, 0.01, 1, 0);
    fitFunc2->SetParameters(h2->GetMaximum(), proton_mass_sq, 0.01, 1, 0);

    h1->Fit(fitFunc1, "R");
    h2->Fit(fitFunc2, "R");

    // Set line style and color for the fit functions
    fitFunc1->SetLineColor(kRed);
    fitFunc1->SetLineStyle(2); // Dashed line
    fitFunc1->SetLineWidth(2);

    fitFunc2->SetLineColor(kBlue);
    fitFunc2->SetLineStyle(2); // Dashed line
    fitFunc2->SetLineWidth(2);

    // Create canvas for Mx2 plot
    TCanvas* c1 = new TCanvas("c1", "Mx2 Comparison", 800, 600);
    c1->SetMargin(0.12, 0.05, 0.12, 0.05);

    // Set axis labels for Mx2 plot
    h1->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetLabelSize(0.04);
    h1->GetXaxis()->SetTitleOffset(1.0);
    h1->GetYaxis()->SetTitleOffset(1.2);

    // Adjust y-axis range for Mx2 plot
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    h1->SetMaximum(1.2 * TMath::Max(max1, max2));

    // Draw histograms for Mx2 plot
    h1->Draw("E");
    h2->Draw("E SAME");

    // Draw the fitted functions on top of the histograms
    fitFunc1->Draw("SAME");
    fitFunc2->Draw("SAME");

    // Create a legend for Mx2 plot
    TLegend* legend = new TLegend(0.50, 0.8, 0.95, 0.95);
    legend->SetBorderSize(1);  // Solid border
    legend->SetLineColor(kBlack);  // Black border line
    legend->SetFillColor(kWhite);  // White background
    legend->SetTextSize(0.02);  // Smaller text

    // Retrieve fit parameters for Mx2 plot
    double mu1 = fitFunc1->GetParameter(1);
    double sigma1 = fitFunc1->GetParameter(2);
    double counts1 = h1->GetEntries();

    double mu2 = fitFunc2->GetParameter(1);
    double sigma2 = fitFunc2->GetParameter(2);
    double counts2 = h2->GetEntries();

    // Add entries to legend for Mx2 plot
    legend->AddEntry(h1, Form("cj 10.1.0: N=%.0f, #mu=%.3f, #sigma=%.3f", counts1, mu1, sigma1), "l");
    legend->AddEntry(h2, Form("DCBetaTimeWalk: N=%.0f, #mu=%.3f, #sigma=%.3f", counts2, mu2, sigma2), "l");
    legend->Draw();

    // Save the Mx2 plot
    c1->SaveAs("/home/thayward/dc_study/dc_comparison.png");

    // ----------------------------------------
    // Plot 2: vz_e Comparison (Updated)
    // ----------------------------------------

    // Define histograms for vz_e
    TH1F* h1_vz = new TH1F("h1_vz", "", 100, 23, 30);
    TH1F* h2_vz = new TH1F("h2_vz", "", 100, 23, 30);

    h1_vz->SetLineColor(kRed);
    h2_vz->SetLineColor(kBlue);
    h1_vz->SetLineWidth(2);
    h2_vz->SetLineWidth(2);

    // Fill histograms for vz_e
    tree1->Draw("vz_e>>h1_vz", cuts);
    tree2->Draw("vz_e>>h2_vz", cuts);

    // Fit histograms with Gaussian plus linear polynomial
    TF1* fitFunc1_vz = new TF1("fitFunc1_vz", "gaus(0)+pol1(3)", 23, 30);
    TF1* fitFunc2_vz = new TF1("fitFunc2_vz", "gaus(0)+pol1(3)", 23, 30);

    // Set initial parameters for vz_e fit: [Amplitude, Mean, Sigma, p0, p1]
    fitFunc1_vz->SetParameters(h1_vz->GetMaximum(), h1_vz->GetMean(), h1_vz->GetRMS(), 1, 0);
    fitFunc2_vz->SetParameters(h2_vz->GetMaximum(), h2_vz->GetMean(), h2_vz->GetRMS(), 1, 0);

    h1_vz->Fit(fitFunc1_vz, "R");
    h2_vz->Fit(fitFunc2_vz, "R");

    // Set line style and color for the fit functions
    fitFunc1_vz->SetLineColor(kRed);
    fitFunc1_vz->SetLineStyle(2); // Dashed line
    fitFunc1_vz->SetLineWidth(2);

    fitFunc2_vz->SetLineColor(kBlue);
    fitFunc2_vz->SetLineStyle(2); // Dashed line
    fitFunc2_vz->SetLineWidth(2);

    // Create canvas for vz_e plot
    TCanvas* c2 = new TCanvas("c2", "vz_e Comparison", 800, 600);
    c2->SetMargin(0.12, 0.05, 0.12, 0.05);

    // Set axis labels for vz_e plot
    h1_vz->GetXaxis()->SetTitle("v_{z}^{e} (cm)");
    h1_vz->GetYaxis()->SetTitle("Counts");
    h1_vz->GetXaxis()->SetTitleSize(0.05);
    h1_vz->GetYaxis()->SetTitleSize(0.05);
    h1_vz->GetXaxis()->SetLabelSize(0.04);
    h1_vz->GetYaxis()->SetLabelSize(0.04);
    h1_vz->GetXaxis()->SetTitleOffset(1.0);
    h1_vz->GetYaxis()->SetTitleOffset(1.2);

    // Adjust y-axis range for vz_e plot
    double max1_vz = h1_vz->GetMaximum();
    double max2_vz = h2_vz->GetMaximum();
    h1_vz->SetMaximum(1.2 * TMath::Max(max1_vz, max2_vz));

    // Draw histograms for vz_e plot
    h1_vz->Draw("E");
    h2_vz->Draw("E SAME");

    // Draw the fitted functions on top of the histograms
    fitFunc1_vz->Draw("SAME");
    fitFunc2_vz->Draw("SAME");

    // Create a legend for vz_e plot
    TLegend* legend_vz = new TLegend(0.50, 0.8, 0.95, 0.95);
    legend_vz->SetBorderSize(1);  // Solid border
    legend_vz->SetLineColor(kBlack);  // Black border line
    legend_vz->SetFillColor(kWhite);  // White background
    legend_vz->SetTextSize(0.02);  // Smaller text

    // Retrieve fit parameters for vz_e plot
    double mu1_vz = fitFunc1_vz->GetParameter(1);
    double sigma1_vz = fitFunc1_vz->GetParameter(2);
    double counts1_vz = h1_vz->GetEntries();

    double mu2_vz = fitFunc2_vz->GetParameter(1);
    double sigma2_vz = fitFunc2_vz->GetParameter(2);
    double counts2_vz = h2_vz->GetEntries();

    // Add entries to legend for vz_e plot
    legend_vz->AddEntry(h1_vz, Form("cj 10.1.0: N=%.0f, #mu=%.2f, #sigma=%.2f", counts1_vz, mu1_vz, sigma1_vz), "l");
    legend_vz->AddEntry(h2_vz, Form("DCBetaTimeWalk: N=%.0f, #mu=%.2f, #sigma=%.2f", counts2_vz, mu2_vz, sigma2_vz), "l");
    legend_vz->Draw();

    // Save the vz_e plot
    c2->SaveAs("/home/thayward/dc_study/dc_comparison_vze.png");

    // ----------------------------------------
    // Expanded Analyses: Subplots for e_theta bins
    // ----------------------------------------

    // Define e_theta bins in degrees and convert to radians
    const int nThetaBins = 5;
    double thetaBinsDeg[nThetaBins+1] = {5, 10, 15, 20, 25, 30};
    double thetaBinsRad[nThetaBins+1];
    for (int i = 0; i <= nThetaBins; ++i) {
        thetaBinsRad[i] = thetaBinsDeg[i] * TMath::Pi() / 180.0;
    }

    // Create canvases for subplots
    TCanvas* c3 = new TCanvas("c3", "Mx2 Comparison by e_theta bins", 1200, 800);
    c3->Divide(3, 2); // 2x3 canvas

    TCanvas* c4 = new TCanvas("c4", "vz_e Comparison by e_theta bins", 1200, 800);
    c4->Divide(3, 2); // 2x3 canvas

    // Vectors to hold histograms and fit functions for cleanup
    std::vector<TH1F*> h1_theta_v, h2_theta_v;
    std::vector<TF1*> fit1_theta_v, fit2_theta_v;
    std::vector<TLegend*> leg_theta_v;

    std::vector<TH1F*> h1_vz_theta_v, h2_vz_theta_v;
    std::vector<TF1*> fit1_vz_theta_v, fit2_vz_theta_v;
    std::vector<TLegend*> leg_vz_theta_v;

    // Loop over e_theta bins
    for (int i = 0; i < nThetaBins; ++i) {
        // Define cuts for this theta bin
        double thetaMin = thetaBinsRad[i];
        double thetaMax = thetaBinsRad[i+1];
        TString thetaCut = Form("fiducial_status==3 && detector1==1 && detector2==1 && e_theta>%f && e_theta<=%f", thetaMin, thetaMax);

        // Histograms for Mx2
        TH1F* h1_theta = new TH1F(Form("h1_theta_%d", i), "", 33, 0.6, 1.25);
        TH1F* h2_theta = new TH1F(Form("h2_theta_%d", i), "", 33, 0.6, 1.25);
        h1_theta->SetLineColor(kRed);
        h2_theta->SetLineColor(kBlue);
        h1_theta->SetLineWidth(2);
        h2_theta->SetLineWidth(2);

        // Fill histograms
        tree1->Draw(Form("Mx2>>h1_theta_%d", i), thetaCut);
        tree2->Draw(Form("Mx2>>h2_theta_%d", i), thetaCut);

        // Fit histograms
        TF1* fit1_theta = new TF1(Form("fit1_theta_%d", i), "gaus(0)+pol1(3)", 0.6, 1.25);
        TF1* fit2_theta = new TF1(Form("fit2_theta_%d", i), "gaus(0)+pol1(3)", 0.6, 1.25);
        fit1_theta->SetParameters(h1_theta->GetMaximum(), proton_mass_sq, 0.01, 1, 0);
        fit2_theta->SetParameters(h2_theta->GetMaximum(), proton_mass_sq, 0.01, 1, 0);

        h1_theta->Fit(fit1_theta, "R");
        h2_theta->Fit(fit2_theta, "R");

        fit1_theta->SetLineColor(kRed);
        fit1_theta->SetLineStyle(2);
        fit1_theta->SetLineWidth(2);

        fit2_theta->SetLineColor(kBlue);
        fit2_theta->SetLineStyle(2);
        fit2_theta->SetLineWidth(2);

        // Store histograms and fits for cleanup
        h1_theta_v.push_back(h1_theta);
        h2_theta_v.push_back(h2_theta);
        fit1_theta_v.push_back(fit1_theta);
        fit2_theta_v.push_back(fit2_theta);

        // Draw on canvas c3
        c3->cd(i+1);
        h1_theta->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h1_theta->GetYaxis()->SetTitle("Counts");
        h1_theta->SetTitle(Form("%.0f#circ < #theta_{e} < %.0f#circ", thetaBinsDeg[i], thetaBinsDeg[i+1]));

        // Adjust y-axis range
        double max_h1 = h1_theta->GetMaximum();
        double max_h2 = h2_theta->GetMaximum();
        h1_theta->SetMaximum(1.2 * TMath::Max(max_h1, max_h2));

        h1_theta->Draw("E");
        h2_theta->Draw("E SAME");
        fit1_theta->Draw("SAME");
        fit2_theta->Draw("SAME");

        // Legend for subplot
        TLegend* leg_theta = new TLegend(0.50, 0.75, 0.95, 0.95);
        leg_theta->SetBorderSize(1);
        leg_theta->SetLineColor(kBlack);
        leg_theta->SetFillColor(kWhite);
        leg_theta->SetTextSize(0.02);

        double mu1_theta = fit1_theta->GetParameter(1);
        double counts1_theta = h1_theta->GetEntries();

        double mu2_theta = fit2_theta->GetParameter(1);
        double counts2_theta = h2_theta->GetEntries();

        leg_theta->AddEntry(h1_theta, Form("cj 10.1.0: N=%.0f, #mu=%.3f", counts1_theta, mu1_theta), "l");
        leg_theta->AddEntry(h2_theta, Form("DCBetaTimeWalk: N=%.0f, #mu=%.3f", counts2_theta, mu2_theta), "l");
        leg_theta->Draw();

        leg_theta_v.push_back(leg_theta);

        // Histograms for vz_e
        TH1F* h1_vz_theta = new TH1F(Form("h1_vz_theta_%d", i), "", 100, 23, 30);
        TH1F* h2_vz_theta = new TH1F(Form("h2_vz_theta_%d", i), "", 100, 23, 30);
        h1_vz_theta->SetLineColor(kRed);
        h2_vz_theta->SetLineColor(kBlue);
        h1_vz_theta->SetLineWidth(2);
        h2_vz_theta->SetLineWidth(2);

        // Fill histograms
        tree1->Draw(Form("vz_e>>h1_vz_theta_%d", i), thetaCut);
        tree2->Draw(Form("vz_e>>h2_vz_theta_%d", i), thetaCut);

        // Fit histograms
        TF1* fit1_vz_theta = new TF1(Form("fit1_vz_theta_%d", i), "gaus(0)+pol1(3)", 23, 30);
        TF1* fit2_vz_theta = new TF1(Form("fit2_vz_theta_%d", i), "gaus(0)+pol1(3)", 23, 30);
        fit1_vz_theta->SetParameters(h1_vz_theta->GetMaximum(), h1_vz_theta->GetMean(), h1_vz_theta->GetRMS(), 1, 0);
        fit2_vz_theta->SetParameters(h2_vz_theta->GetMaximum(), h2_vz_theta->GetMean(), h2_vz_theta->GetRMS(), 1, 0);

        h1_vz_theta->Fit(fit1_vz_theta, "R");
        h2_vz_theta->Fit(fit2_vz_theta, "R");

        fit1_vz_theta->SetLineColor(kRed);
        fit1_vz_theta->SetLineStyle(2);
        fit1_vz_theta->SetLineWidth(2);

        fit2_vz_theta->SetLineColor(kBlue);
        fit2_vz_theta->SetLineStyle(2);
        fit2_vz_theta->SetLineWidth(2);

        // Store histograms and fits for cleanup
        h1_vz_theta_v.push_back(h1_vz_theta);
        h2_vz_theta_v.push_back(h2_vz_theta);
        fit1_vz_theta_v.push_back(fit1_vz_theta);
        fit2_vz_theta_v.push_back(fit2_vz_theta);

        // Draw on canvas c4
        c4->cd(i+1);
        h1_vz_theta->GetXaxis()->SetTitle("v_{z}^{e} (cm)");
        h1_vz_theta->GetYaxis()->SetTitle("Counts");
        h1_vz_theta->SetTitle(Form("%.0f#circ < #theta_{e} < %.0f#circ", thetaBinsDeg[i], thetaBinsDeg[i+1]));

        // Adjust y-axis range
        double max_h1_vz = h1_vz_theta->GetMaximum();
        double max_h2_vz = h2_vz_theta->GetMaximum();
        h1_vz_theta->SetMaximum(1.2 * TMath::Max(max_h1_vz, max_h2_vz));

        h1_vz_theta->Draw("E");
        h2_vz_theta->Draw("E SAME");
        fit1_vz_theta->Draw("SAME");
        fit2_vz_theta->Draw("SAME");

        // Legend for subplot
        TLegend* leg_vz_theta = new TLegend(0.50, 0.75, 0.95, 0.95);
        leg_vz_theta->SetBorderSize(1);
        leg_vz_theta->SetLineColor(kBlack);
        leg_vz_theta->SetFillColor(kWhite);
        leg_vz_theta->SetTextSize(0.02);

        double mu1_vz_theta = fit1_vz_theta->GetParameter(1);
        double counts1_vz_theta = h1_vz_theta->GetEntries();

        double mu2_vz_theta = fit2_vz_theta->GetParameter(1);
        double counts2_vz_theta = h2_vz_theta->GetEntries();

        leg_vz_theta->AddEntry(h1_vz_theta, Form("cj 10.1.0: N=%.0f, #mu=%.2f", counts1_vz_theta, mu1_vz_theta), "l");
        leg_vz_theta->AddEntry(h2_vz_theta, Form("DCBetaTimeWalk: N=%.0f, #mu=%.2f", counts2_vz_theta, mu2_vz_theta), "l");
        leg_vz_theta->Draw();

        leg_vz_theta_v.push_back(leg_vz_theta);
    }

    // Save the canvases with subplots
    c3->SaveAs("/home/thayward/dc_study/dc_comparison_theta_Mx2.png");
    c4->SaveAs("/home/thayward/dc_study/dc_comparison_theta_vze.png");

    // Clean up
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete h1;
    delete h2;
    delete fitFunc1;
    delete fitFunc2;
    delete legend;
    delete h1_vz;
    delete h2_vz;
    delete fitFunc1_vz;
    delete fitFunc2_vz;
    delete legend_vz;

    // Clean up histograms and fits for subplots
    for (size_t i = 0; i < h1_theta_v.size(); ++i) {
        delete h1_theta_v[i];
        delete h2_theta_v[i];
        delete fit1_theta_v[i];
        delete fit2_theta_v[i];
        delete leg_theta_v[i];

        delete h1_vz_theta_v[i];
        delete h2_vz_theta_v[i];
        delete fit1_vz_theta_v[i];
        delete fit2_vz_theta_v[i];
        delete leg_vz_theta_v[i];
    }

    file1->Close();
    file2->Close();

    return 0;
}