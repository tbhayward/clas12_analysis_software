#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TMath.h>
#include <TF1.h>
#include <iostream>

void plot_dvcs_energy_loss_validation(const char* file1, const char* file2, const char* titleSuffix) {
    // Open the ROOT files
    TFile *f1 = new TFile(file1);
    TFile *f2 = new TFile(file2);

    // Access the "PhysicsEvents" trees
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Prepare variables for reading the branches
    Double_t p1_theta1, p1_theta2;
    Double_t Mxgammasquared_1, Mxgammasquared_2;
    Double_t eta2_1, eta2_2;
    Double_t t1_1, t1_2;
    Double_t theta_gamma_gamma_1, theta_gamma_gamma_2;
    Double_t Emiss2_1, Emiss2_2;
    Double_t pTmiss_1, pTmiss_2;

    // Set branch addresses
    tree1->SetBranchAddress("p1_theta", &p1_theta1);
    tree1->SetBranchAddress("Mxgammasquared", &Mxgammasquared_1);
    tree1->SetBranchAddress("eta2", &eta2_1);
    tree1->SetBranchAddress("t1", &t1_1);
    tree1->SetBranchAddress("theta_gamma_gamma", &theta_gamma_gamma_1);
    tree1->SetBranchAddress("Emiss2", &Emiss2_1);
    tree1->SetBranchAddress("pTmiss", &pTmiss_1);

    tree2->SetBranchAddress("p1_theta", &p1_theta2);
    tree2->SetBranchAddress("Mxgammasquared", &Mxgammasquared_2);
    tree2->SetBranchAddress("eta2", &eta2_2);
    tree2->SetBranchAddress("t1", &t1_2);
    tree2->SetBranchAddress("theta_gamma_gamma", &theta_gamma_gamma_2);
    tree2->SetBranchAddress("Emiss2", &Emiss2_2);
    tree2->SetBranchAddress("pTmiss", &pTmiss_2);

    // Create canvas and divide it into 3x4 subplots
    TCanvas *c1 = new TCanvas("c1", "DVCS Energy Loss Validation", 1200, 900);
    c1->Divide(4, 3);

    // Set the theta bins and corresponding histogram ranges
    const int nBins = 11;
    Double_t thetaBins[nBins + 1] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 65}; // 5 to 65 degrees

    TH1D *h1[nBins];
    TH1D *h2[nBins];

    for (int i = 0; i < nBins; ++i) {
        c1->cd(i + 1);

        // Create histograms for each theta bin
        h1[i] = new TH1D(Form("h1_%d", i), Form("M_{xp}^{2} GeV^{2} for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 100, -0.5, 0.5);
        h2[i] = new TH1D(Form("h2_%d", i), Form("M_{xp}^{2} GeV^{2} for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 100, -0.5, 0.5);

        // Fill the histograms with the cuts applied
        Long64_t nEntries1 = tree1->GetEntries();
        for (Long64_t j = 0; j < nEntries1; ++j) {
            tree1->GetEntry(j);
            Double_t thetaDeg1 = p1_theta1 * (180.0 / TMath::Pi()); // Convert to degrees
            if (thetaDeg1 >= thetaBins[i] && thetaDeg1 < thetaBins[i + 1] &&
                eta2_1 < 0 && t1_1 > -2 && theta_gamma_gamma_1 < 0.6 &&
                Emiss2_1 < 0.5 && pTmiss_1 < 0.125) {
                h1[i]->Fill(Mxgammasquared_1);
            }
        }

        Long64_t nEntries2 = tree2->GetEntries();
        for (Long64_t j = 0; j < nEntries2; ++j) {
            tree2->GetEntry(j);
            Double_t thetaDeg2 = p1_theta2 * (180.0 / TMath::Pi()); // Convert to degrees
            if (thetaDeg2 >= thetaBins[i] && thetaDeg2 < thetaBins[i + 1] &&
                eta2_2 < 0 && t1_2 > -2 && theta_gamma_gamma_2 < 0.6 &&
                Emiss2_2 < 0.5 && pTmiss_2 < 0.125) {
                h2[i]->Fill(Mxgammasquared_2);
            }
        }

        // Draw histograms as points with error bars
        h1[i]->SetMarkerStyle(20);
        h1[i]->SetMarkerColor(kBlack);
        h1[i]->Draw("E");
        h2[i]->SetMarkerStyle(21);
        h2[i]->SetMarkerColor(kRed);
        h2[i]->Draw("E SAME");

        // Fit histograms to Gaussian
        TF1 *fit1 = new TF1(Form("fit1_%d", i), "gaus", -2, 2);
        TF1 *fit2 = new TF1(Form("fit2_%d", i), "gaus", -2, 2);
        h1[i]->Fit(fit1, "Q");
        h2[i]->Fit(fit2, "Q");

        // Draw fitted functions
        fit1->SetLineColor(kBlack);
        fit1->Draw("SAME");
        fit2->SetLineColor(kRed);
        fit2->Draw("SAME");

        // Get fit parameters
        Double_t mu1 = fit1->GetParameter(1);
        Double_t sigma1 = fit1->GetParameter(2);
        Double_t mu2 = fit2->GetParameter(1);
        Double_t sigma2 = fit2->GetParameter(2);

        // Add legend with mu and sigma values
        TLegend *legend = new TLegend(0.6, 0.75, 0.9, 0.9);
        legend->AddEntry(h1[i], Form("Uncorrected: #mu=%.3f, #sigma=%.3f", mu1, sigma1), "lep");
        legend->AddEntry(h2[i], Form("Corrected: #mu=%.3f, #sigma=%.3f", mu2, sigma2), "lep");
        legend->Draw();

        // Label the axes
        h1[i]->GetXaxis()->SetTitle("M_{xp}^{2} GeV^{2}");
        h1[i]->GetYaxis()->SetTitle("Counts");
    }

    // Save the canvas as a PDF
    c1->SaveAs("output/dvcs_energy_loss_validation.pdf");

    // Clean up
    delete c1;
    delete f1;
    delete f2;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <titleSuffix>" << std::endl;
        return 1;
    }

    plot_dvcs_energy_loss_validation(argv[1], argv[2], argv[3]);
    return 0;
}