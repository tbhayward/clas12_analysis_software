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
    // Plot 1: Mx2 Comparison (Same as before)
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
    TF1* fitFunc1 = new TF1("fitFunc1", "gaus(0)+pol1(3)", 0.6, 1.4);
    TF1* fitFunc2 = new TF1("fitFunc2", "gaus(0)+pol1(3)", 0.6, 1.4);

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
    c1->SaveAs("/home/thayward/dc_comparison.png");

    // ----------------------------------------
    // Plot 2: e_vz Comparison
    // ----------------------------------------

    // Define histograms for e_vz
    TH1F* h1_vz = new TH1F("h1_vz", "", 100, -10, 5);
    TH1F* h2_vz = new TH1F("h2_vz", "", 100, -10, 5);

    h1_vz->SetLineColor(kRed);
    h2_vz->SetLineColor(kBlue);
    h1_vz->SetLineWidth(2);
    h2_vz->SetLineWidth(2);

    // Fill histograms for e_vz
    tree1->Draw("vz_e>>h1_vz", cuts);
    tree2->Draw("vz_e>>h2_vz", cuts);

    // Fit histograms with Gaussian
    TF1* fitFunc1_vz = new TF1("fitFunc1_vz", "gaus", -15, 10);
    TF1* fitFunc2_vz = new TF1("fitFunc2_vz", "gaus", -15, 10);

    // Set initial parameters for e_vz fit: [Amplitude, Mean, Sigma]
    fitFunc1_vz->SetParameters(h1_vz->GetMaximum(), h1_vz->GetMean(), h1_vz->GetRMS());
    fitFunc2_vz->SetParameters(h2_vz->GetMaximum(), h2_vz->GetMean(), h2_vz->GetRMS());

    h1_vz->Fit(fitFunc1_vz, "R");
    h2_vz->Fit(fitFunc2_vz, "R");

    // Set line style and color for the fit functions
    fitFunc1_vz->SetLineColor(kRed);
    fitFunc1_vz->SetLineStyle(2); // Dashed line
    fitFunc1_vz->SetLineWidth(2);

    fitFunc2_vz->SetLineColor(kBlue);
    fitFunc2_vz->SetLineStyle(2); // Dashed line
    fitFunc2_vz->SetLineWidth(2);

    // Create canvas for e_vz plot
    TCanvas* c2 = new TCanvas("c2", "e_vz Comparison", 800, 600);
    c2->SetMargin(0.12, 0.05, 0.12, 0.05);

    // Set axis labels for e_vz plot
    h1_vz->GetXaxis()->SetTitle("v_{z}^{e}");
    h1_vz->GetYaxis()->SetTitle("Counts");
    h1_vz->GetXaxis()->SetTitleSize(0.05);
    h1_vz->GetYaxis()->SetTitleSize(0.05);
    h1_vz->GetXaxis()->SetLabelSize(0.04);
    h1_vz->GetYaxis()->SetLabelSize(0.04);
    h1_vz->GetXaxis()->SetTitleOffset(1.0);
    h1_vz->GetYaxis()->SetTitleOffset(1.2);

    // Adjust y-axis range for e_vz plot
    double max1_vz = h1_vz->GetMaximum();
    double max2_vz = h2_vz->GetMaximum();
    h1_vz->SetMaximum(1.2 * TMath::Max(max1_vz, max2_vz));

    // Draw histograms for e_vz plot
    h1_vz->Draw("E");
    h2_vz->Draw("E SAME");

    // Draw the fitted functions on top of the histograms
    fitFunc1_vz->Draw("SAME");
    fitFunc2_vz->Draw("SAME");

    // Create a legend for e_vz plot
    TLegend* legend_vz = new TLegend(0.50, 0.8, 0.95, 0.95);
    legend_vz->SetBorderSize(1);  // Solid border
    legend_vz->SetLineColor(kBlack);  // Black border line
    legend_vz->SetFillColor(kWhite);  // White background
    legend_vz->SetTextSize(0.02);  // Smaller text

    // Retrieve fit parameters for e_vz plot
    double mu1_vz = fitFunc1_vz->GetParameter(1);
    double sigma1_vz = fitFunc1_vz->GetParameter(2);
    double counts1_vz = h1_vz->GetEntries();

    double mu2_vz = fitFunc2_vz->GetParameter(1);
    double sigma2_vz = fitFunc2_vz->GetParameter(2);
    double counts2_vz = h2_vz->GetEntries();

    // Add entries to legend for e_vz plot
    legend_vz->AddEntry(h1_vz, Form("cj 10.1.0: N=%.0f, #mu=%.2f, #sigma=%.2f", counts1_vz, mu1_vz, sigma1_vz), "l");
    legend_vz->AddEntry(h2_vz, Form("DCBetaTimeWalk: N=%.0f, #mu=%.2f, #sigma=%.2f", counts2_vz, mu2_vz, sigma2_vz), "l");
    legend_vz->Draw();

    // Save the e_vz plot
    c2->SaveAs("/home/thayward/dc_comparison_vze.png");

    // Clean up
    delete c1;
    delete c2;
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
    file1->Close();
    file2->Close();

    return 0;
}