#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>

void analyzeROOTFile(const char* fileName) {
    // Open the ROOT file
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return;
    }

    // Load the tree
    TTree* tree;
    file->GetObject("PhysicsEvents", tree);
    if (!tree) {
        std::cerr << "Tree 'PhysicsEvents' not found in file: " << fileName << std::endl;
        return;
    }

    // Define a variable to hold the branch value
    double Mx;
    tree->SetBranchAddress("Mx", &Mx);

    // Create a histogram for Mx
    TH1D* hMx = new TH1D("hMx", "Mx Distribution;Mx (GeV/c^2);Entries", 125, 0.6, 1.1);

    // Loop over the tree entries and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        hMx->Fill(Mx);
    }

    // Define a fitting function: Gaussian + second degree polynomial
    TF1* fitFunc = new TF1("fitFunc", "gaus(0)+pol2(3)", 0, 1.5);
    // Gaussian initial parameters: constant, mean (neutron mass), sigma (0.01)
    fitFunc->SetParameters(500, 0.939, 0.01);
    // Polynomial initial parameters can remain as the default

    // Perform the fit
    hMx->Fit(fitFunc, "R");

    // Plot the histogram and the fitted function
    TCanvas* canvas = new TCanvas("canvas", "Mx Fit", 800, 600);
    hMx->SetLineColor(kBlack);
    hMx->Draw();
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    // Save the canvas to an image file
    canvas->SaveAs("output.png");

    // Print out the mean and sigma values of the Gaussian component
    std::cout << "Mean = " << fitFunc->GetParameter(1) << ", Sigma = " << fitFunc->GetParameter(2) << std::endl;

    // Cleanup
    delete file; // Automatically deletes associated objects like TTree, TH1D, etc.
    delete canvas;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <ROOT file>" << std::endl;
        return 1;
    }
    analyzeROOTFile(argv[1]);
    return 0;
}
