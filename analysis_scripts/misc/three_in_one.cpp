// combine_root_trees.cpp
// Combine the "PhysicsEvents" trees from three ROOT files into one output file.
//
// Usage:
//   ./combine_root_trees file1.root file2.root file3.root output.root
//
// Build (csh-friendly):
//   g++ -O2 -std=c++17 combine_root_trees.cpp `root-config --cflags --libs` -o combine_root_trees

#include <iostream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

static bool file_has_tree(const std::string& path, const std::string& treename) {
    TFile f(path.c_str(), "READ");
    if (f.IsZombie()) {
        std::cerr << "[ERROR] Cannot open file: " << path << "\n";
        return false;
    }
    TTree* t = dynamic_cast<TTree*>(f.Get(treename.c_str()));
    if (!t) {
        std::cerr << "[ERROR] Tree \"" << treename << "\" not found in file: " << path << "\n";
        f.Close();
        return false;
    }
    f.Close();
    return true;
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " file1.root file2.root file3.root output.root\n";
        return 1;
    }

    const std::string in1 = argv[1];
    const std::string in2 = argv[2];
    const std::string in3 = argv[3];
    const std::string out = argv[4];
    const std::string treename = "PhysicsEvents";

    // Basic sanity checks on inputs
    if (!file_has_tree(in1, treename)) return 2;
    if (!file_has_tree(in2, treename)) return 2;
    if (!file_has_tree(in3, treename)) return 2;

    // Create a chain and add inputs
    TChain chain(treename.c_str());
    chain.Add(in1.c_str());
    chain.Add(in2.c_str());
    chain.Add(in3.c_str());

    if (chain.GetNtrees() != 3) {
        std::cerr << "[ERROR] Expected to attach 3 input files, but chain has " << chain.GetNtrees() << "\n";
        return 3;
    }

    std::cout << "[INFO] Total entries across inputs: " << chain.GetEntries() << "\n";

    // Create output file and clone the chain into it
    TFile fout(out.c_str(), "RECREATE");
    if (fout.IsZombie()) {
        std::cerr << "[ERROR] Cannot create output file: " << out << "\n";
        return 4;
    }

    fout.cd();
    // Use "fast" option when possible (fast cloning), fallback handled by ROOT internally
    TTree* merged = chain.CloneTree(-1, "fast");
    if (!merged) {
        std::cerr << "[ERROR] CloneTree returned null. Aborting.\n";
        fout.Close();
        return 5;
    }

    merged->SetName(treename.c_str()); // ensure the tree is named consistently
    merged->Write("", TObject::kOverwrite);
    fout.Write();
    fout.Close();

    std::cout << "[OK] Wrote combined tree \"" << treename << "\" with "
              << merged->GetEntries() << " entries to " << out << "\n";
    return 0;
}