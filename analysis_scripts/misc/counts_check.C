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
    Double_t x, W, A;  // Changed x to Double_t, added W and A
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("x", &x);  // Corrected type for x
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("A", &A);

    // Variables for overall mean calculations
    double sum_x = 0, sum_W = 0, sum_A = 0;
    int count_x = 0, count_WA = 0;

    // Map to store counts for each runnum
    std::map<int, int> counts;

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Count entries per runnum
        counts[runnum]++;

        // For x calculations (excluding run 16297)
        if (runnum != 16297) {
            sum_x += x;
            count_x++;
        }

        // For W and A calculations
        sum_W += W;
        sum_A += A;
        count_WA++;
    }

    // Calculate overall means
    double mean_x = (count_x > 0) ? sum_x / count_x : 0;
    double mean_W = (count_WA > 0) ? sum_W / count_WA : 0;
    double mean_A = (count_WA > 0) ? sum_A / count_WA : 0;
    double mean_W_over_A = (mean_A > 0) ? mean_W / mean_A : 0;

    // Print results
    std::cout << "Entries per Runnum:" << std::endl;
    for (const auto &count : counts) {
        std::cout << "Runnum " << count.first << ": " << count.second << " entries" << std::endl;
    }

    std::cout << "\nOverall mean x (excluding run 16297): " << mean_x << std::endl;
    std::cout << "Overall mean W over mean A: " << mean_W_over_A << std::endl;

    // Close the file
    file->Close();
}
