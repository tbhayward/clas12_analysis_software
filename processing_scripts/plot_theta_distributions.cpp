#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>

void plot_theta_distributions() {
    // Open the ROOT file and get the tree
    TFile *file = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/cphi2024_rho0/data/new_sep6/rgc_su22_inb_NH3_eppi+pi-X.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: PhysicsEvents tree not found!" << std::endl;
        file->Close();
        return;
    }

    // Create histograms for symmetric and antisymmetric conditions
    TH1F *h_symmetric = new TH1F("h_symmetric", "Symmetric", 100, 0, TMath::Pi() / 2);
    TH1F *h_antisymmetric = new TH1F("h_antisymmetric", "Antisymmetric", 100, 0, TMath::Pi() / 2);

    // Set axis labels
    h_symmetric->GetXaxis()->SetTitle("#theta");
    h_symmetric->GetYaxis()->SetTitle("Counts");
    h_antisymmetric->GetXaxis()->SetTitle("#theta");
    h_antisymmetric->GetYaxis()->SetTitle("Counts");

    // Apply the selection cuts and fill the histograms
    tree->Draw("theta >> h_symmetric", "abs(p2_p - p3_p) < 1", "goff");
    tree->Draw("theta >> h_antisymmetric", "abs(p2_p - p3_p) > 1", "goff");

    // Create a canvas with 1x2 subplots
    TCanvas *c = new TCanvas("c", "Theta Distributions", 1200, 600);
    c->Divide(2, 1);

    // Draw the histograms
    c->cd(1);
    h_symmetric->Draw();
    c->cd(2);
    h_antisymmetric->Draw();

    // Save the canvas as an image
    c->SaveAs("/home/thayward/theta_distributions.png");

    // Clean up
    file->Close();
    delete c;
    delete h_symmetric;
    delete h_antisymmetric;
}