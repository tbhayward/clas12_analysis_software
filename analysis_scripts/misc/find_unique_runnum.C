#include <iostream>
#include <set>
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

void find_unique_runnum(TString file1, TString file2) {
    // Open the first ROOT file
    TFile *f1 = TFile::Open(file1);
    if (!f1 || !f1->IsOpen()) {
        std::cerr << "Error opening file: " << file1 << std::endl;
        return;
    }

    // Open the second ROOT file
    TFile *f2 = TFile::Open(file2);
    if (!f2 || !f2->IsOpen()) {
        std::cerr << "Error opening file: " << file2 << std::endl;
        return;
    }

    // Load the "PhysicsEvents" tree from both files
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    if (!tree1 || !tree2) {
        std::cerr << "Error: Tree 'PhysicsEvents' not found in one or both files." << std::endl;
        return;
    }

    // Set up a set to store unique runnum values
    std::set<int> uniqueRunnums;

    // Variable to hold the "runnum" value
    int runnum;

    // Process the first tree
    tree1->SetBranchAddress("runnum", &runnum);
    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; i++) {
        tree1->GetEntry(i);
        uniqueRunnums.insert(runnum);
    }

    // Process the second tree
    tree2->SetBranchAddress("runnum", &runnum);
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; i++) {
        tree2->GetEntry(i);
        uniqueRunnums.insert(runnum);
    }

    // Print out the unique runnum values
    std::cout << "Unique runnum values:" << std::endl;
    for (const auto& r : uniqueRunnums) {
        std::cout << r << std::endl;
    }

    // Clean up
    f1->Close();
    f2->Close();
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: root -l -b -q 'find_unique_runnum.C(\"file1.root\", \"file2.root\")'" << std::endl;
        return 1;
    }

    find_unique_runnum(argv[1], argv[2]);

    return 0;
}