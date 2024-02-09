#include <TFile.h>
#include <TTree.h>
#include <iostream>

int main(int argc, char *argv[]) {
    // Check for correct number of arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <inputFile.root> <outputFile.root>" << std::endl;
        return 1;
    }

    // Open the input file
    TFile *inputFile = new TFile(argv[1], "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error opening input file!" << std::endl;
        return 1;
    }

    // Get the tree from the input file
    TTree *inputTree = nullptr;
    inputFile->GetObject("PhysicsEvents", inputTree);
    if (!inputTree) {
        std::cerr << "Error: Tree 'PhysicsEvents' not found in input file!" << std::endl;
        inputFile->Close();
        return 1;
    }

    // Create the output file and a clone of the tree to work on
    TFile *outputFile = new TFile(argv[2], "RECREATE");
    TTree *outputTree = inputTree->CloneTree(0); // Clone structure, but no entries yet

    // Set branch addresses for input tree
    double Mx, Mx1, Mx2, z1, xF1, xF2;
    inputTree->SetBranchAddress("Mx", &Mx);
    inputTree->SetBranchAddress("Mx1", &Mx1);
    inputTree->SetBranchAddress("Mx2", &Mx2);
    inputTree->SetBranchAddress("z1", &z1);
    inputTree->SetBranchAddress("xF1", &xF1);
    inputTree->SetBranchAddress("xF2", &xF2);

    // Loop over events and apply cuts
    Long64_t nEntries = inputTree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        inputTree->GetEntry(i);

        if (Mx > 0.95 && Mx1 > 1.8 && Mx2 > 1.4 && z1 > 0.2 && xF1 > 0.0 && xF2 < 0.0) {
            outputTree->Fill(); // Fill the output tree with the current entry if it passes the cuts
        }
    }

    // Write the output tree to the output file
    outputFile->cd();
    outputTree->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "Processing complete. Output saved to " << argv[2] << std::endl;

    return 0;
}
