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
    Int_t runnum;
    Float_t x;
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("x", &x);

    // Maps to store counts and sums for each runnum
    std::map<int, int> counts;
    std::map<int, double> sum_x;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Skip run 16297 for x calculations
        if (runnum != 16297) {
            sum_x[runnum] += x;
        }

        // Count entries per runnum
        counts[runnum]++;
    }

    // Calculate means and print results
    std::cout << "Runnum\tEntries\tMean x (excluding run 16297)" << std::endl;
    for (const auto &count : counts) {
        double mean_x = (count.first != 16297) ? sum_x[count.first] / count.second : 0;
        std::cout << count.first << "\t" << count.second << "\t" << mean_x << std::endl;
    }

    // Close the file
    file->Close();
}

