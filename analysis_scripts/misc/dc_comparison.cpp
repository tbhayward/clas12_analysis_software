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

    // **Decrease the number of bins by a factor of two**
    // Define histograms with 25 bins (half of 50) and starting at 0.6
    TH1F* h1 = new TH1F("h1", "", 25, 0.6, 1.2);
    TH1F* h2 = new TH1F("h2", "", 25, 0.6, 1.2);

    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h1->SetLineWidth(2);
    h2->SetLineWidth(2);

    // Define cuts
    double e_theta_cut = 5.0 * TMath::Pi() / 180.0; // Convert degrees to radians
    TString cuts = Form("fiducial_status==3 && detector1==1 && detector2==1 && e_theta>%f", e_theta_cut);

    // Fill histograms
    tree1->Draw("Mx2>>h1", cuts);
    tree2->Draw("Mx2>>h2", cuts);

    // **Adjust the fit to be Gaussian plus linear polynomial**
    // Fit histograms with Gaussian plus linear polynomial within the new range
    TF1* fitFunc1 = new TF1("fitFunc1", "gaus(0)+pol1(3)", 0.6, 1.4);
    TF1* fitFunc2 = new TF1("fitFunc2", "gaus(0)+pol1(3)", 0.6, 1.4);

    double proton_mass_sq = 0.938 * 0.938; // Proton mass squared in GeV^2

    // Set initial parameters: [Amplitude, Mean, Sigma, p0, p1]
    fitFunc1->SetParameters(h1->GetMaximum(), proton_mass_sq, 0.01, 1, 0);
    fitFunc2->SetParameters(h2->GetMaximum(), proton_mass_sq, 0.01, 1, 0);

    h1->Fit(fitFunc1, "R");
    h2->Fit(fitFunc2, "R");

    // Add the fitted functions as colored dashed lines
    // Set line style and color for the fit functions
    fitFunc1->SetLineColor(kRed);
    fitFunc1->SetLineStyle(2); // Dashed line
    fitFunc1->SetLineWidth(2);

    fitFunc2->SetLineColor(kBlue);
    fitFunc2->SetLineStyle(2); // Dashed line
    fitFunc2->SetLineWidth(2);

    // Create canvas
    TCanvas* c1 = new TCanvas("c1", "Mx2 Comparison", 800, 600);
    c1->SetMargin(0.12, 0.05, 0.12, 0.05);

    // Set axis labels
    h1->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h1->GetYaxis()->SetTitle("Counts");
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);
    h1->GetXaxis()->SetLabelSize(0.04);
    h1->GetYaxis()->SetLabelSize(0.04);
    h1->GetXaxis()->SetTitleOffset(1.0);
    h1->GetYaxis()->SetTitleOffset(1.2);

    // Adjust y-axis range
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    h1->SetMaximum(1.2 * TMath::Max(max1, max2));

    // Draw histograms
    h1->Draw("E");
    h2->Draw("E SAME");

    // Draw the fitted functions on top of the histograms
    fitFunc1->Draw("SAME");
    fitFunc2->Draw("SAME");

    // **Adjust the legend**
    // Create a smaller legend in the top right corner with a solid black border
    TLegend* legend = new TLegend(0.50, 0.8, 0.95, 0.95);
    legend->SetBorderSize(1);  // Solid border
    legend->SetLineColor(kBlack);  // Black border line
    legend->SetFillColor(kWhite);  // White background
    legend->SetTextSize(0.02);  // Much smaller text

    // Retrieve fit parameters
    double mu1 = fitFunc1->GetParameter(1);
    double sigma1 = fitFunc1->GetParameter(2);
    double counts1 = h1->GetEntries();

    double mu2 = fitFunc2->GetParameter(1);
    double sigma2 = fitFunc2->GetParameter(2);
    double counts2 = h2->GetEntries();

    // Add entries to legend
    legend->AddEntry(h1, Form("cj 10.1.0: N=%.0f, #mu=%.3f, #sigma=%.3f", counts1, mu1, sigma1), "l");
    legend->AddEntry(h2, Form("DCBetaTimeWalk: N=%.0f, #mu=%.3f, #sigma=%.3f", counts2, mu2, sigma2), "l");
    legend->Draw();

    // Save the canvas
    c1->SaveAs("/home/thayward/dc_comparison.png");

    // Clean up
    delete c1;
    delete h1;
    delete h2;
    delete fitFunc1;
    delete fitFunc2;
    delete legend;
    file1->Close();
    file2->Close();

    return 0;
}