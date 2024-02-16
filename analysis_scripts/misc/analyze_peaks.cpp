#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <set>
#include <vector>

void analyzePions() {
    // Open the ROOT file
    TFile* file = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/mc/epi+X/rga_fa18_inb_clasdis_50nA_epi+X.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    
    // Access the tree
    TTree* tree = dynamic_cast<TTree*>(file->Get("PhysicsEvents"));
    if (!tree) {
        std::cerr << "Error: Tree 'PhysicsEvents' not found" << std::endl;
        file->Close();
        return;
    }

    // Define variables to hold branch data
    Double_t Q2, W, y, xF, Mx, z, p_p, pT;
    Int_t mc_p1_parent;

    // Set branch addresses
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("xF", &xF);
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("pT", &pT); // Assuming pT is a branch; adjust as needed
    tree->SetBranchAddress("mc_p1_parent", &mc_p1_parent);

    // Container for unique mc_p1_parent values
    std::set<int> uniqueParents;

    // Loop through the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Apply the kinematic conditions
        if (Q2 > 1 && W > 2 && y < 0.75 && xF > 0 && Mx > 1.5 && z > 0.2 && p_p > 1.2 &&
            y > 0.65 && y < 0.75 && z > 0.24 && z < 0.29 && Q2 > 2.0 && Q2 < 2.5) {
            // Store unique mc_p1_parent
            uniqueParents.insert(mc_p1_parent);
        }
    }

    // Print unique mc_p1_parent values
    std::cout << "Unique mc_p1_parent values in the filtered set:" << std::endl;
    for (int parent : uniqueParents) {
        std::cout << parent << std::endl;
    }

    // Define pT bins (adjust if necessary)
    std::vector<double> pTBins = {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00};

    // Cleanup
    file->Close();
    delete file;
}

int main() {
    analyzePions();
    return 0;
}
