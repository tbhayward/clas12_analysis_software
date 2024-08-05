#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <set>

void find_missing_events(const char* file1, const char* file2, const char* prefix) {
    // Open the ROOT files
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);

    if (!f1 || !f2 || f1->IsZombie() || f2->IsZombie()) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }

    // Get the PhysicsEvents tree
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    if (!tree1 || !tree2) {
        std::cerr << "Error: Could not find the PhysicsEvents tree in one of the files" << std::endl;
        return;
    }

    // Set up the evnum branch
    int evnum1, evnum2;
    tree1->SetBranchAddress("evnum", &evnum1);
    tree2->SetBranchAddress("evnum", &evnum2);

    // Store evnum from the first file in a set
    std::set<int> evnum_set;
    Long64_t nentries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nentries1; ++i) {
        tree1->GetEntry(i);
        evnum_set.insert(evnum1);
    }

    // Open the output ROOT file and text file
    TString output_root_file = TString(prefix) + ".root";
    TString output_txt_file = TString(prefix) + ".txt";
    TFile *fout = new TFile(output_root_file, "RECREATE");
    TTree *new_tree = tree2->CloneTree(0); // Clone structure of tree2 but with 0 entries
    std::ofstream txt_out(output_txt_file);

    // Loop over the second tree and find events not in the first file
    Long64_t nentries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nentries2; ++i) {
        tree2->GetEntry(i);
        if (evnum_set.find(evnum2) == evnum_set.end()) {
            new_tree->Fill(); // Fill the tree with the missing event
            txt_out << evnum2 << std::endl; // Write the event number to the text file
        }
    }

    // Write and close the output ROOT file and text file
    new_tree->Write();
    fout->Close();
    txt_out.close();

    // Clean up
    delete f1;
    delete f2;
    delete fout;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " file1.root file2.root prefix" << std::endl;
        return 1;
    }

    find_missing_events(argv[1], argv[2], argv[3]);
    return 0;
}