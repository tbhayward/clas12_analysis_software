#include <iostream>
#include <map>
#include <TFile.h>
#include <TTree.h>

void counts_check() {
    // Open the ROOT file
    TFile *file = TFile::Open("/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epX.root");
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
    Double_t x, DepolarizationA, DepolarizationW;
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("helicity", &helicity);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("DepA", &DepolarizationA);
    tree->SetBranchAddress("DepW", &DepolarizationW);

    // Counters and sums for calculations
    int npp = 0, npm = 0, nmp = 0, nmm = 0;
    double sum_x = 0, sum_DepW = 0, sum_DepA = 0;
    int count_x = 0, count_Dep = 0;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Skip run 16297
        if (runnum == 16297) continue;

        // Add to sums and counts for mean calculations
        sum_x += x;
        count_x++;
        sum_DepW += DepolarizationW;
        sum_DepA += DepolarizationA;
        count_Dep++;

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

    // Calculate overall means
    double mean_x = (count_x > 0) ? sum_x / count_x : 0;
    double mean_DepW_over_DepA = (count_Dep > 0 && sum_DepA != 0) ? sum_DepW / sum_DepA : 0;

    // Print results
    std::cout << "Counts:" << std::endl;
    std::cout << "npp: " << npp << std::endl;
    std::cout << "npm: " << npm << std::endl;
    std::cout << "nmp: " << nmp << std::endl;
    std::cout << "nmm: " << nmm << std::endl;
    std::cout << "\nMean x for runs excluding 16297: " << mean_x << std::endl;
    std::cout << "Mean DepolarizationW over DepolarizationA: " << mean_DepW_over_DepA << std::endl;

    // Close the file
    file->Close();
}
