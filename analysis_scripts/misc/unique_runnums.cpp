#include <iostream>
#include <set>
#include "TFile.h"
#include "TTree.h"

int main() {
    // Open the ROOT file
    TFile *file = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epi+pX/rgc_su22_inb_NH3_epi+pX.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Get the tree from the file
    TTree *tree;
    file->GetObject("PhysicsData", tree);
    if (!tree) {
        std::cerr << "Error: PhysicsData tree not found!" << std::endl;
        file->Close();
        return 1;
    }

    // Set the branch address for runnum
    Int_t runnum;
    tree->SetBranchAddress("runnum", &runnum);

    // Use a set to store unique run numbers
    std::set<Int_t> unique_runnums;

    // Loop over the entries in the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        unique_runnums.insert(runnum);
    }

    // Print the unique run numbers
    std::cout << "Unique run numbers:" << std::endl;
    for (const auto &num : unique_runnums) {
        std::cout << num << std::endl;
    }

    // Print the total number of unique run numbers
    std::cout << "Total number of unique run numbers: " << unique_runnums.size() << std::endl;

    // Clean up
    file->Close();

    return 0;
}