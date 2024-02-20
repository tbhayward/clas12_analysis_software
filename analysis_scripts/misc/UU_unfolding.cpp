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
    // Open the ROOT files for data and Monte Carlo
    TFile* fData = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_inb_epi+X_skimmed.root");
    TFile* fMCReco = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_inb_clasdis_50nA_epi+X_skimmed.root");
    if (!fData || fData->IsZombie() || !fMCReco || fMCReco->IsZombie()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return -1;
    }

    // Get the TTrees
    TTree* tData;
    TTree* tMCReco;
    fData->GetObject("PhysicsEvents", tData);
    fMCReco->GetObject("PhysicsEvents", tMCReco);
    if (!tData || !tMCReco) {
        std::cerr << "Tree not found in one or more files." << std::endl;
        return -1;
    }

    // Define variables for both trees
    double Q2Data, yData, phiData, pTData, zData;
    double Q2MC, yMC, phiMC, pTMC, zMC;
    tData->SetBranchAddress("Q2", &Q2Data);
    tData->SetBranchAddress("y", &yData);
    tData->SetBranchAddress("phi", &phiData);
    tData->SetBranchAddress("pT", &pTData);
    tData->SetBranchAddress("z", &zData);
    tMCReco->SetBranchAddress("Q2", &Q2MC);
    tMCReco->SetBranchAddress("y", &yMC);
    tMCReco->SetBranchAddress("phi", &phiMC);
    tMCReco->SetBranchAddress("pT", &pTMC);
    tMCReco->SetBranchAddress("z", &zMC);

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
    std::vector<TH1F*> hData, hMCReco;
    for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
        hData.push_back(new TH1F(Form("hData_%d", i), ";#phi;Normalized Counts", 24, 0, 2*3.14159));
        hMCReco.push_back(new TH1F(Form("hMCReco_%d", i), ";#phi;Normalized Counts", 24, 0, 2*3.14159));
    }

    std::cout << "Looping over data." << std::endl;
    // Fill histograms for data
    Long64_t nDataEntries = tData->GetEntries();
    for (Long64_t i = 0; i < nDataEntries; ++i) {
        tData->GetEntry(i);
        if (DetermineQ2yBin(Q2Data, yData) != 1) continue;

        // Determine the corresponding pT and z bins
        int pT_bin = -1, z_bin = -1;
        for (int j = 0; j < num_pT_bins; ++j) {
            if (pTData > pT_edges[j] && pTData <= pT_edges[j+1]) {
                pT_bin = j;
                break;
            }
        }
        for (int k = 0; k < num_z_bins; k++) {
            if (zData > z_edges[k] && zData <= z_edges[k+1]) {
                z_bin = num_z_bins - k - 1;
                break;
            }
        }
        // Fill the corresponding histogram if the event is in a valid bin
        if (pT_bin != -1 && z_bin != -1) {
            int histIndex = z_bin * num_pT_bins + pT_bin;
            hData[histIndex]->Fill(phiData);
        }
    }

    std::cout << "Looping over reconstructed MC." << std::endl;
    // Fill histograms for rec MC
    Long64_t nMCEntries = hMCReco->GetEntries();
    for (Long64_t i = 0; i < nMCEntries; ++i) {
        hMCReco->GetEntry(i);
        if (DetermineQ2yBin(Q2MC, yMC) != 1) continue;

        // Determine the corresponding pT and z bins
        int pT_bin = -1, z_bin = -1;
        for (int j = 0; j < num_pT_bins; ++j) {
            if (pTMC > pT_edges[j] && pTMC <= pT_edges[j+1]) {
                pT_bin = j;
                break;
            }
        }
        for (int k = 0; k < num_z_bins; k++) {
            if (zMC > z_edges[k] && zMC <= z_edges[k+1]) {
                z_bin = num_z_bins - k - 1;
                break;
            }
        }
        // Fill the corresponding histogram if the event is in a valid bin
        if (pT_bin != -1 && z_bin != -1) {
            int histIndex = z_bin * num_pT_bins + pT_bin;
            hMCReco[histIndex]->Fill(phiMC);
        }
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~ */ 
    // Declare the TLatex object here, before the loop
    TLatex latex;
    latex.SetTextSize(0.08);
    latex.SetNDC();

    int currentQ2yBin = 1;
    // Loop over the histograms to draw them
    for (size_t i = 0; i < hData.size(); ++i) {
        if (hData[i]->GetEntries() > 0 && hMCReco[i]->GetEntries() > 0) {
            // Navigate to the correct pad
            canvas->cd(i + 1);

            // Adjust pad margins to add space around the plots
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.15);
            gPad->SetTopMargin(0.15);
            gPad->SetBottomMargin(0.15);

            // Remove the stat box
            hData[i]->SetStats(0); hMCReco[i]->SetStats(0);

            // Change the line color to a darker blue
            hData[i]->SetLineColor(kBlue+2);
            hData[i]->SetLineWidth(2); // Increase line width
            hMCReco[i]->SetLineColor(kBlue+2);
            hMCReco[i]->SetLineWidth(2); // Increase line width

            // Increase font size for axis labels
            hData[i]->GetXaxis()->SetLabelSize(0.04); // Adjust as needed
            hData[i]->GetYaxis()->SetLabelSize(0.04); // Adjust as needed
            // Increase font size for axis titles
            hData[i]->GetXaxis()->SetTitleSize(0.04); // Adjust as needed
            hData[i]->GetYaxis()->SetTitleSize(0.04); // Adjust as needed

            // Draw the histogram
            hData[i]->DrawNormalized("HIST");
            hMCReco[i]->DrawNormalized("HIST same");

            // Display z-pT bin information as 'z-P_{T} bin: histIndex'
            // Note: Adjust the positioning (x, y coordinates) as needed
            latex.DrawLatexNDC(0.35, 0.88, Form("Q2-y bin: %d, z-P_{T} bin: %zu", currentQ2yBin, i));
        }
    }


    // Save the canvas
    canvas->SaveAs("output/Q2y1.png");

    // Cleanup
    delete canvas;
    for (auto& hist : hData) delete hist;
    fData->Close();
    fMC->Close();
    delete fData,fMCReco;

    return 0;
}
