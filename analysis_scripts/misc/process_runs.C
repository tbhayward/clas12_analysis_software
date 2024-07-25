#include <iostream>
#include <set>
#include <TFile.h>
#include <TTree.h>

void process_runs() {
    // Open the ROOT file
    TFile *file = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/epX/rgc_su22_inb_C_epX.root");

    // Check if file was opened successfully
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file." << std::endl;
        return;
    }

    // Get the TTree object named "PhysicsEvents"
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error accessing PhysicsEvents tree." << std::endl;
        file->Close();
        return;
    }

    // Variable to hold the values of runnum
    Int_t runnum;
    
    // Set branch address
    tree->SetBranchAddress("runnum", &runnum);

    // Set to store unique runnum entries
    std::set<Int_t> unique_runs;

    // Loop over all entries in the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        std::cout << runnum << " " << std::endl;
        unique_runs.insert(runnum);
    }

    // Print the total number of unique entries
    std::cout << "Total number of unique runnum entries: " << unique_runs.size() << std::endl;

    // Print each unique entry
    for (const auto& run : unique_runs) {
        std::cout << run << std::endl;
    }

    // Close the file
    file->Close();
}