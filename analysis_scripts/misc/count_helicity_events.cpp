#include <iostream>
#include <TFile.h>
#include <TTree.h>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <ROOT file path>" << std::endl;
        return 1;
    }

    const char* fileName = argv[1];
    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return 1;
    }

    TTree* tree = nullptr;
    file->GetObject("PhysicsEvents", tree);
    if (!tree) {
        std::cerr << "Error: 'PhysicsEvents' tree not found in file: " << fileName << std::endl;
        file->Close();
        return 1;
    }

    int helicity = 0;
    tree->SetBranchAddress("helicity", &helicity);

    int countHelicity1 = 0, countHelicity0 = 0, countHelicityMinus1 = 0;
    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t iEntry = 0; iEntry < nEntries; ++iEntry) {
        tree->GetEntry(iEntry);
        
        if (helicity == 1) ++countHelicity1;
        else if (helicity == 0) ++countHelicity0;
        else if (helicity == -1) ++countHelicityMinus1;
    }

    std::cout << "Number of events with helicity = 1: " << countHelicity1 << std::endl;
    std::cout << "Number of events with helicity = 0: " << countHelicity0 << std::endl;
    std::cout << "Number of events with helicity = -1: " << countHelicityMinus1 << std::endl;

    file->Close();
    return 0;
}
