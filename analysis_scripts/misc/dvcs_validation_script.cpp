#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TString.h>
#include <TMath.h>
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
    Double_t thetaBins[nBins + 1] = {5, 11, 17, 23, 29, 35, 41, 47, 53, 59, 65, 70}; // 5 to 70 degrees

    TH1D *h1[nBins];
    TH1D *h2[nBins];

    for (int i = 0; i < nBins; ++i) {
        c1->cd(i + 1);

        // Create histograms for each theta bin
        h1[i] = new TH1D(Form("h1_%d", i), Form("M_{xp}^{2} GeV^{2} for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 100, -2, 2);
        h2[i] = new TH1D(Form("h2_%d", i), Form("M_{xp}^{2} GeV^{2} for #theta [%.0f, %.0f] %s", thetaBins[i], thetaBins[i + 1], titleSuffix), 100, -2, 2);

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

        // Draw histograms on the same pad
        h1[i]->SetLineColor(kBlack);
        h1[i]->Draw("HIST");
        h2[i]->SetLineColor(kRed);
        h2[i]->Draw("HIST SAME");

        // Remove the stat box
        h1[i]->SetStats(0);
        h2[i]->SetStats(0);

        // Add legend
        TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
        legend->AddEntry(h1[i], "Uncorrected", "l");
        legend->AddEntry(h2[i], "Corrected", "l");
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