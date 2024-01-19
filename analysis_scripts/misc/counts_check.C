#include <iostream>
#include <map>
#include <TFile.h>
#include <TTree.h>

void counts_check() {
    // Open the ROOT file
    TFile *file = TFile::Open("/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_eX.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Access the tree
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error accessing tree 'PhysicsEvents'" << std::endl;
        return;
    }

    // Set up variables to read branches
    Int_t runnum, helicity;
    Double_t DepolarizationA, DepolarizationW;
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("helicity", &helicity);
    tree->SetBranchAddress("DepolarizationA", &DepolarizationA);
    tree->SetBranchAddress("DepolarizationW", &DepolarizationW);

    // Counters for npp, npm, nmp, nmm
    int npp = 0, npm = 0, nmp = 0, nmm = 0;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Skip run 16297
        if (runnum == 16297) continue;

        // Determine target polarization
        bool isPositiveTarget = (runnum == 16320 || runnum == 16327);
        bool isNegativeTarget = (runnum == 16346 || runnum == 16353);

        // Increment counters based on helicity and target polarization
        if (helicity > 0) {  // Positive helicity
            if (isPositiveTarget) npp++;
            if (isNegativeTarget) npm++;
        } else if (helicity < 0) {  // Negative helicity
            if (isPositiveTarget) nmp++;
            if (isNegativeTarget) nmm++;
        }
    }

    // Print results
    std::cout << "Counts:" << std::endl;
    std::cout << "npp: " << npp << std::endl;
    std::cout << "npm: " << npm << std::endl;
    std::cout << "nmp: " << nmp << std::endl;
    std::cout << "nmm: " << nmm << std::endl;

    // Close the file
    file->Close();
}
