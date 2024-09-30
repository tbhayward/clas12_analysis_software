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
    TH2D *hist = new TH2D("Q2_vs_epsilon", "", 100, 0, 6, 100, 0, 1);
    hist->SetTitle("Q^{2} vs #epsilon for x_{B} ~ 0.16; Q^{2} (GeV^{2}); #epsilon = DepB / DepA");

    // Loop over the events and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < 1000000; ++i) {
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

// Add a main function for standalone compilation
int main(int argc, char** argv) {
    if (argc != 2) {
        printf("Usage: %s input_file.root\n", argv[0]);
        return 1;
    }
    const char* inputFile = argv[1];
    plot_Q2_vs_epsilon(inputFile);
    return 0;
}