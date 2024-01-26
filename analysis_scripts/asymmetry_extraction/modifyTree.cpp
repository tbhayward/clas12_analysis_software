// modifyTree.cpp
#include "modifyTree.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>

void modifyTree(const char* inputFileName, const char* outputFileName) {
    // Open the original file
    TFile* inputFile = new TFile(inputFileName, "READ");
    if (!inputFile || !inputFile->IsOpen()) {
        std::cerr << "Error opening input file." << std::endl;
        return;
    }

    // Get the tree from the file
    TTree* tree = (TTree*)inputFile->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error getting tree from input file." << std::endl;
        inputFile->Close();
        return;
    }

    // Clone the tree into a new file
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    TTree* clonedTree = tree->CloneTree(0); // Clone without copying the entries

    // Create or modify the "runnum" branch
    Int_t runnumValue = 11;
    TBranch* runnumBranch = tree->GetBranch("runnum");
    TBranch* newRunnumBranch = nullptr;
    if (runnumBranch) {
        newRunnumBranch = clonedTree->Branch("runnum", &runnumValue, "runnum/I");
    } else {
        newRunnumBranch = clonedTree->Branch("runnum", &runnumValue, "runnum/I");
    }

    // Create or modify the "target_pol" branch
    Double_t targetPolValue = 0.0;
    TBranch* targetPolBranch = tree->GetBranch("target_pol");
    TBranch* newTargetPolBranch = nullptr;
    if (targetPolBranch) {
        newTargetPolBranch = clonedTree->Branch("target_pol", &targetPolValue, "target_pol/D");
    } else {
        newTargetPolBranch = clonedTree->Branch("target_pol", &targetPolValue, "target_pol/D");
    }

    // Copy the entries from the original tree to the cloned tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i); // Load the original entry
        runnumValue = 11;  // Modify "runnum"
        targetPolValue = 0.0; // Modify "target_pol"
        clonedTree->Fill(); // Fill the cloned tree with the modified entry
    }

    // Write and close the new file
    clonedTree->Write();
    outputFile->Close();
    inputFile->Close();

    std::cout << "New file with modified tree created: " << outputFileName << std::endl;
}
