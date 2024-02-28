#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <string>

void filterEvents(const char* inputFile) {
    // Generate the output file name
    std::string outputFile = std::string(inputFile).substr(0, std::string(inputFile).rfind(".")) + "_skimmed.root";

    // Open the input ROOT file
    TFile* inFile = TFile::Open(inputFile, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error opening file: " << inputFile << std::endl;
        return;
    }

    // Load the PhysicsEvents tree
    TTree* tree;
    inFile->GetObject("PhysicsEvents", tree);
    if (!tree) {
        std::cerr << "Error: 'PhysicsEvents' tree not found in file: " << inputFile << std::endl;
        inFile->Close();
        return;
    }

    // Define variables to hold the branch data
    double p_p, Q2, z, xF, Mx, y;

    // Set branch addresses
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("xF", &xF);
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("y", &y);

    // Create a new file and clone the tree structure
    TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
    TTree* outTree = tree->CloneTree(0); // Clone structure but don't copy the data yet

    // Filtering logic
    Long64_t nentries = tree->GetEntries();
    Long64_t nselected = 0;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        if (p_p > 1.25 && Q2 > 2.0 && z > 0.15 && xF > 0 && Mx > 1.5 && y > 0.30 && y < 0.75) {
            outTree->Fill(); // Copy this entry to the output tree
            nselected++;
        }
    }

    // for (Long64_t i = 0; i < nentries; i++) {
    //     tree->GetEntry(i);
    //     if (Mx > 1.4 && y < 0.75) {
    //         outTree->Fill(); // Copy this entry to the output tree
    //         nselected++;
    //     }
    // }

    // Save the skimmed tree and close files
    outTree->AutoSave();
    outFile->Close();
    inFile->Close();

    // Print the results
    std::cout << "Initial number of events: " << nentries << std::endl;
    std::cout << "Number of events after filtering: " << nselected << std::endl;
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input ROOT file>" << std::endl;
        return 1;
    }

    filterEvents(argv[1]);

    return 0;
}
