#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>

void plot_Q2_vs_epsilon(const char* inputFile) {
    // Open the input ROOT file
    TFile *file = TFile::Open(inputFile);
    if (!file || file->IsZombie()) {
        printf("Error opening file %s\n", inputFile);
        return;
    }

    // Load the PhysicsEvents tree
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        printf("Error: Tree 'PhysicsEvents' not found in file %s\n", inputFile);
        return;
    }

    // Set up the variables
    double Q2, DepA, DepB, x;
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("DepA", &DepA);
    tree->SetBranchAddress("DepB", &DepB);
    tree->SetBranchAddress("x", &x);

    // Create a 2D histogram
    TH2D *hist = new TH2D("Q2_vs_epsilon", "Q^{2} vs #epsilon;Q^{2} (GeV^{2});#epsilon = DepB / DepA", 
                          100, 0, 10, 100, 0, 2);

    // Loop over the events and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (x >= 0.14 && x <= 0.18 && DepA != 0) {
            double epsilon = DepB / DepA;
            hist->Fill(Q2, epsilon);
        }
    }

    // Draw and save the histogram
    TCanvas *c = new TCanvas("c", "Q^{2} vs #epsilon", 800, 600);
    gStyle->SetOptStat(0);
    hist->Draw("COLZ");
    c->SaveAs("output/Q2vsepsilon.png");

    // Clean up
    delete c;
    delete hist;
    file->Close();
    delete file;
}