#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>

// Function to determine the Q2-y bin based on given Q2 and y values.
// This function is a simplified version for Q2-y bin 1, adjust for full range as needed.
int DetermineQ2yBin(float Q2, float y) {
    if (Q2 > 2.000 && Q2 < 2.423) {
        if (y > 0.650 && y < 0.750) {
            return 1; // Q2-y bin 1
        }
    }
    return -1; // Not in Q2-y bin 1
}

// Main function
int main() {
    // Open the ROOT file
    TFile* file = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_inb_epi+X_skimmed.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open the file." << std::endl;
        return -1;
    }

    // Get the TTree
    TTree* tree;
    file->GetObject("PhysicsEvents", tree); // Replace "tree_name" with the actual tree name
    if (!tree) {
        std::cerr << "Tree not found in the file." << std::endl;
        return -1;
    }

    // Define variables to hold the data from the tree branches
    double Q2, y, phi, pT, z;
    // Assuming these are the branch names, adjust as necessary
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("z", &z);

    // Define the bin edges for z and pT
    float pT_edges[] = {0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0};
    float z_edges[] = {0.15, 0.2, 0.24, 0.29, 0.4, 0.73};
    int num_pT_bins = sizeof(pT_edges)/sizeof(float) - 1;
    int num_z_bins = sizeof(z_edges)/sizeof(float) - 1;

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Q2-y Bin 1 Phi Distributions", 2000, 1200);
    canvas->Divide(num_pT_bins, num_z_bins); // Adjust the division based on the number of bins

    std::cout << std::endl << "Creating histograms." << std::endl;
    // Create histograms for each z-pT bin
    std::vector<TH1F*> histograms;
    for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
        histograms.push_back(new TH1F(Form("hist_%d", i), "Phi Distribution;Phi;Counts", 24, 0, 2*3.14159));
    }

    std::cout << "Looping over data." << std::endl;
    // Loop over the tree and fill histograms
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Determine if the event is in Q2-y bin 1
        if (DetermineQ2yBin(Q2, y) != 1) continue;

        // Determine the corresponding pT and z bins
        int pT_bin = -1, z_bin = -1;
        for (int j = 0; j < num_pT_bins; ++j) {
            if (pT > pT_edges[j] && pT <= pT_edges[j+1]) {
                pT_bin = j;
                break;
            }
        }
        for (int k = num_z_bins - 1; k >= 0; k--) {
            if (z > z_edges[k] && z <= z_edges[k+1]) {
                z_bin = k;
                break;
            }
        }

        // Fill the corresponding histogram if the event is in a valid bin
        if (pT_bin != -1 && z_bin != -1) {
            int histIndex = pT_bin * num_z_bins + (num_z_bins - z_bin - 1);
            histograms[histIndex]->Fill(phi);
            // std::cout << z_bin << " " << pT_bin << " " << (z_bin * num_pT_bins + pT_bin + 1) << std::endl;
        }
    }

    // Plot and save histograms
    for (size_t i = 0; i < histograms.size(); ++i) {
        if (histograms[i]->GetEntries() > 0) {
            canvas->cd(i + 1);
            histograms[i]->DrawNormalized();
        }
    }

    // canvas->cd(0 + 1);
    // histograms[0]->DrawNormalized();

    // Save the canvas
    canvas->SaveAs("output/Q2y1.png");

    // Cleanup
    delete canvas;
    for (auto& hist : histograms) delete hist;
    file->Close();
    delete file;

    return 0;
}
