#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>

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
        histograms.push_back(new TH1F(Form("hist_%d", i), ";#phi;Normalized Counts", 24, 0, 2*3.14159));
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
        for (int k = 0; k < num_z_bins; k++) {
            if (z > z_edges[k] && z <= z_edges[k+1]) {
                z_bin = num_z_bins - k - 1;
                break;
            }
        }
        // Fill the corresponding histogram if the event is in a valid bin
        if (pT_bin != -1 && z_bin != -1) {
            int histIndex = z_bin * num_pT_bins + pT_bin;
            histograms[histIndex]->Fill(phi);
        }
    }


    /* ~~~~~~~~~~~~~~~~~~~~~~ */ 
    // Declare the TLatex object here, before the loop
    TLatex latex;
    latex.SetTextSize(0.05);
    latex.SetNDC();

    int currentQ2yBin = 1;
    // Loop over the histograms to draw them
    for (size_t i = 0; i < histograms.size(); ++i) {
        if (histograms[i]->GetEntries() > 0) {
            // Navigate to the correct pad
            canvas->cd(i + 1);

            // Adjust pad margins to add space around the plots
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.15);
            gPad->SetTopMargin(0.15);
            gPad->SetBottomMargin(0.15);

            // Remove the stat box
            histograms[i]->SetStats(0);

            // Change the line color to a darker blue
            histograms[i]->SetLineColor(kBlue+2);
            histograms[i]->SetLineWidth(2); // Increase line width

            // Increase font size for axis labels
            histograms[i]->GetXaxis()->SetLabelSize(0.04); // Adjust as needed
            histograms[i]->GetYaxis()->SetLabelSize(0.04); // Adjust as needed
            
            // Increase font size for axis titles
            histograms[i]->GetXaxis()->SetTitleSize(0.04); // Adjust as needed
            histograms[i]->GetYaxis()->SetTitleSize(0.04); // Adjust as needed

            // Draw the histogram
            histograms[i]->DrawNormalized("HIST");

            // Display z-pT bin information as 'z-P_{T} bin: histIndex'
            // Note: Adjust the positioning (x, y coordinates) as needed
            latex.DrawLatexNDC(0.5, 0.85, Form("Q2-y bin: %d, z-P_{T} bin: %zu", currentQ2yBin, i));
        }
    }


    // Save the canvas
    canvas->SaveAs("output/Q2y1.png");

    // Cleanup
    delete canvas;
    for (auto& hist : histograms) delete hist;
    file->Close();
    delete file;

    return 0;
}
