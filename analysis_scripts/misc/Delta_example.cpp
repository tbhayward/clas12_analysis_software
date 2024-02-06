#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
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

    // Define variables to hold the branch values
    double mc_p_p, p_p;
    tree->SetBranchAddress("mc_p_p", &mc_p_p);
    tree->SetBranchAddress("p_p", &p_p);

    // Create a histogram for Delta p
    TH1F* hDeltaP = new TH1F("hDeltaP", "#Delta p (GeV);#Delta p;Entries", 100, -0.3, 0.3);

    // Loop over the tree entries and fill the histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        float deltaP = mc_p_p - p_p;
        hDeltaP->Fill(deltaP);
    }

    // Fit the histogram to a Gaussian
    TF1* gaussFit = new TF1("gaussFit", "gaus", -0.3, 0.3);
    hDeltaP->Fit(gaussFit, "RQ"); // "RQ" option for Range and Quiet mode

    // Plot the histogram and the fitted function
    TCanvas* canvas = new TCanvas("canvas", "Delta p Analysis", 800, 600);
    hDeltaP->SetLineColor(kBlack);
    hDeltaP->Draw();
    gaussFit->SetLineColor(kRed);
    gaussFit->Draw("same");

    // Save the canvas to an image file
    canvas->SaveAs("output.png");

    // Print out the mean and sigma values
    std::cout << "Mean = " << gaussFit->GetParameter(1) << ", Sigma = " << gaussFit->GetParameter(2) << std::endl;

    // Cleanup
    delete file; // Automatically deletes associated objects like TTree, TH1F, etc.
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
