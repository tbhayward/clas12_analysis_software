#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <set>
#include <vector>
#include <map> // Include map header

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
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("mc_p1_parent", &mc_p1_parent);

    // Containers for unique mc_p1_parent values and event counts
    std::map<int, int> parentEventCounts;
    int totalEventsMeetingCriteria = 0;

    // Loop through the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Apply the kinematic conditions
        if (Q2 > 1 && W > 2 && y < 0.75 && xF > 0 && Mx > 1.5 && z > 0.2 && p_p > 1.2 &&
            y > 0.65 && y < 0.75 && z > 0.24 && z < 0.29 && Q2 > 2.0 && Q2 < 2.5) {
            // Increment count for this mc_p1_parent
            parentEventCounts[mc_p1_parent]++;
            totalEventsMeetingCriteria++;
        }
    }

    // Print unique mc_p1_parent values and their corresponding percentages
    std::cout << "Percentage of events for each unique mc_p1_parent:" << std::endl;
    for (const auto& pair : parentEventCounts) {
        double percentage = 100.0 * pair.second / totalEventsMeetingCriteria;
        std::cout << "mc_p1_parent = " << pair.first << ": " << percentage << "%" << std::endl;
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
