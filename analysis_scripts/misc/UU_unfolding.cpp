#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>

// Function to determine the Q2-y bin based on given Q2 and y values.
int DetermineQ2yBin(float Q2, float y) {
    // Q2-y Bins 1-4
    if (Q2 > 2.000 && Q2 <= 2.423) {
        if (y > 0.650 && y <= 0.750) return 1;
        if (y > 0.550 && y <= 0.650) return 2;
        if (y > 0.450 && y <= 0.550) return 3;
        if (y > 0.300 && y <= 0.450) return 4;
    }
    // Q2-y Bins 5-8
    else if (Q2 > 2.423 && Q2 <= 2.987) {
        if (y > 0.650 && y <= 0.750) return 5;
        if (y > 0.550 && y <= 0.650) return 6;
        if (y > 0.450 && y <= 0.550) return 7;
        if (y > 0.300 && y <= 0.450) return 8;
    }
    // Q2-y Bins 9-12
    else if (Q2 > 2.987 && Q2 <= 3.974) {
        if (y > 0.650 && y <= 0.750) return 9;
        if (y > 0.550 && y <= 0.650) return 10;
        if (y > 0.450 && y <= 0.550) return 11;
        if (y > 0.350 && y <= 0.450) return 12;
    }
    // Q2-y Bins 13-14
    else if (Q2 > 3.974 && Q2 <= 5.384) {
        if (y > 0.650 && y <= 0.750) return 13;
        if (y > 0.550 && y <= 0.650) return 14;
    }
    // Q2-y Bin 15
    else if (Q2 > 3.974 && Q2 <= 5.948 && y > 0.450 && y <= 0.550) {
        return 15;
    }
    // Q2-y Bin 16
    else if (Q2 > 5.384 && Q2 <= 9.896 && y > 0.650 && y <= 0.750) {
        return 16;
    }
    // Q2-y Bin 17
    else if (Q2 > 5.384 && Q2 <= 7.922 && y > 0.550 && y <= 0.650) {
        return 17;
    }

    return -1; // Not in any defined Q2-y bin
}


// Main function
int main() {
    // Open the ROOT files for data and Monte Carlo
    TFile* fData = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_inb_epi+X_skimmed.root");
    TFile* fMCReco = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_inb_clasdis_50nA_epi+X_skimmed.root");
    TFile* fMCGene = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_inb_clasdis_50nA_gen_epi+X_skimmed.root");
    if (!fData || fData->IsZombie() || !fMCReco || fMCReco->IsZombie() || !fMCGene || fMCGene->IsZombie()) {
        std::cerr << "Failed to open one or more files." << std::endl;
        return -1;
    }

    // Get the TTrees
    TTree* tData;
    TTree* tMCReco;
    TTree* tMCGene;
    fData->GetObject("PhysicsEvents", tData);
    fMCReco->GetObject("PhysicsEvents", tMCReco);
    fMCGene->GetObject("PhysicsEvents", tMCGene);
    if (!tData || !tMCReco || !tMCGene) {
        std::cerr << "Tree not found in one or more files." << std::endl;
        return -1;
    }

    // Define variables for both trees
    double Q2Data, yData, phiData, pTData, zData;
    double Q2MC, yMC, phiMC, pTMC, zMC;
    double Q2Gen, yGen, phiGen, pTGen, zGen;
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
    tMCGene->SetBranchAddress("Q2", &Q2Gen);
    tMCGene->SetBranchAddress("y", &yGen);
    tMCGene->SetBranchAddress("phi", &phiGen);
    tMCGene->SetBranchAddress("pT", &pTGen);
    tMCGene->SetBranchAddress("z", &zGen);

    // Define the bin edges for z and pT
    float pT_edges[] = {0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0};
    float z_edges[] = {0.15, 0.2, 0.24, 0.29, 0.4, 0.73};
    int num_pT_bins = sizeof(pT_edges)/sizeof(float) - 1;
    int num_z_bins = sizeof(z_edges)/sizeof(float) - 1;

    std::cout << std::endl << "Creating histograms." << std::endl;
    // Create histograms for each z-pT bin
    std::vector<std::vector<TH1F*>> hData(17), hMCReco(17), hMCGene(17);
    for (int bin = 0; bin < 17; ++bin) {
        for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
            hData[bin].push_back(new TH1F(Form("hData_bin%d_%d", bin+1, i), ";#phi;Normalized Counts", 24, 0, 2*3.14159));
            hMCReco[bin].push_back(new TH1F(Form("hMCReco_bin%d_%d", bin+1, i), ";#phi;Normalized Counts", 24, 0, 2*3.14159));
            hMCGene[bin].push_back(new TH1F(Form("hMCGene_bin%d_%d", bin+1, i), ";#phi;Normalized Counts", 24, 0, 2*3.14159));
        }
    }


    std::cout << "Looping over data." << std::endl;
    // Fill histograms for data
    Long64_t nDataEntries = tData->GetEntries();
    for (Long64_t i = 0; i < nDataEntries; ++i) {
        tData->GetEntry(i);
        int binIndex = DetermineQ2yBin(Q2Data, yData) - 1; // Adjusted for 0-based indexing
        if (binIndex >= 0) {

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
                hData[binIndex][histIndex]->Fill(phiData);
            }
        }
    }

    std::cout << "Looping over reconstructed MC." << std::endl;
    // Fill histograms for rec MC
    Long64_t nMCEntries = tMCReco->GetEntries();
    for (Long64_t i = 0; i < nMCEntries; ++i) {
        tMCReco->GetEntry(i);
        int binIndex = DetermineQ2yBin(Q2MC, yMC) - 1; // Adjusted for 0-based indexing
        if (binIndex >= 0) {

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
                hMCReco[binIndex][histIndex]->Fill(phiMC);
            }
        }
    }

    std::cout << "Looping over generated MC." << std::endl;
    // Fill histograms for gen MC
    Long64_t nGenEntries = tMCGene->GetEntries();
    for (Long64_t i = 0; i < nGenEntries; ++i) {
        tMCGene->GetEntry(i);
        int binIndex = DetermineQ2yBin(Q2Gen, yGen) - 1; // Adjusted for 0-based indexing
        if (binIndex >= 0) {

            // Determine the corresponding pT and z bins
            int pT_bin = -1, z_bin = -1;
            for (int j = 0; j < num_pT_bins; ++j) {
                if (pTGen > pT_edges[j] && pTGen <= pT_edges[j+1]) {
                    pT_bin = j;
                    break;
                }
            }
            for (int k = 0; k < num_z_bins; k++) {
                if (zGen > z_edges[k] && zGen <= z_edges[k+1]) {
                    z_bin = num_z_bins - k - 1;
                    break;
                }
            }
            // Fill the corresponding histogram if the event is in a valid bin
            if (pT_bin != -1 && z_bin != -1) {
                int histIndex = z_bin * num_pT_bins + pT_bin;
                hMCGene[binIndex][histIndex]->Fill(phiGen);
            }
        }
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~ */ 
    // Declare the TLatex object here, before the loop
    TLatex latex; latex.SetTextSize(0.09); latex.SetNDC();
    for (int bin = 0; bin < 17; ++bin) {
        TCanvas* canvas = new TCanvas(Form("canvas_bin_%d", bin+1), Form("Q2-y Bin %d Phi Distributions", bin+1), 2000, 1200);
        canvas->Divide(num_pT_bins, num_z_bins);
        
        for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
            canvas->cd(i + 1);

            // Adjust pad margins to add space around the plots
            gPad->SetLeftMargin(0.25);
            gPad->SetRightMargin(0.2);
            gPad->SetTopMargin(0.2);
            gPad->SetBottomMargin(0.2);

            // Remove the stat box
            hData[bin][i]->SetStats(0); hMCReco[bin][i]->SetStats(0); hMCGene[bin][i]->SetStats(0);

            // Change the line color to a darker blue
            hData[bin][i]->SetLineColor(kBlue+2);
            hData[bin][i]->SetLineWidth(2); // Increase line width
            hMCReco[bin][i]->SetLineColor(kRed+2);
            hMCReco[bin][i]->SetLineWidth(2); // Increase line width
            hMCGene[bin][i]->SetLineColor(kGreen+2);
            hMCGene[bin][i]->SetLineWidth(2); // Increase line width

            // Increase font size for axis labels
            hData[bin][i]->GetXaxis()->SetLabelSize(0.06); // Adjust as needed
            hData[bin][i]->GetYaxis()->SetLabelSize(0.06); // Adjust as needed
            hMCReco[bin][i]->GetXaxis()->SetLabelSize(0.06); // Adjust as needed
            hMCReco[bin][i]->GetYaxis()->SetLabelSize(0.06); // Adjust as needed
            hMCGene[bin][i]->GetXaxis()->SetLabelSize(0.06); // Adjust as needed
            hMCGene[bin][i]->GetYaxis()->SetLabelSize(0.06); // Adjust as needed
            // Increase font size for axis titles
            hData[bin][i]->GetXaxis()->SetTitleSize(0.07); // Adjust as needed
            hData[bin][i]->GetYaxis()->SetTitleSize(0.07); // Adjust as needed
            hMCReco[bin][i]->GetXaxis()->SetTitleSize(0.07); // Adjust as needed
            hMCReco[bin][i]->GetYaxis()->SetTitleSize(0.07); // Adjust as needed
            hMCGene[bin][i]->GetXaxis()->SetTitleSize(0.07); // Adjust as needed
            hMCGene[bin][i]->GetYaxis()->SetTitleSize(0.07); // Adjust as needed

            if (hMCReco[bin][i]->GetEntries() > 100) {
                hData[bin][i]->DrawNormalized("HIST");
                hMCReco[bin][i]->DrawNormalized("HIST same");
                hMCGene[bin][i]->DrawNormalized("HIST same");

                // Display z-pT bin information as 'z-P_{T} bin: histIndex'
                // Note: Adjust the positioning (x, y coordinates) as needed
                latex.DrawLatexNDC(0.10, 0.86, Form("Q2-y bin: %d, z-P_{T} bin: %zu", (bin+1), (i+1)));
            }
        }
        
        canvas->SaveAs(Form("output/Q2yBin_%d.png", bin+1));
        delete canvas;
    }


    std::vector<std::vector<TH1F*>> hAcceptance(17);
    for (int bin = 0; bin < 17; ++bin) {
        hAcceptance[bin].resize(num_pT_bins * num_z_bins);
        for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
            // Ensure MC reco histogram has entries to avoid division by zero
            if (hMCReco[bin][i]->GetEntries() > 0) {
                hAcceptance[bin][i] = (TH1F*)hMCGene[bin][i]->Clone(Form("hAcceptance_bin%d_%d", bin+1, i));
                hAcceptance[bin][i]->Divide(hMCReco[bin][i]); // Divide generated by reconstructed
            }
        }
    }

    for (int bin = 0; bin < 17; ++bin) {
    TCanvas* unfoldedCanvas = new TCanvas(Form("unfolded_canvas_bin_%d", bin+1), Form("Unfolded Q2-y Bin %d Phi Distributions", bin+1), 2000, 1200);
    unfoldedCanvas->Divide(num_pT_bins, num_z_bins);
    
    for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
        unfoldedCanvas->cd(i + 1);
        
        // Apply acceptance correction if both histograms have entries
        if (hData[bin][i]->GetEntries() > 0 && hAcceptance[bin][i] != nullptr) {
            TH1F* hUnfolded = (TH1F*)hData[bin][i]->Clone(Form("hUnfolded_bin%d_%d", bin+1, i));
            hUnfolded->Multiply(hAcceptance[bin][i]); // Multiply data by acceptance correction
            
            // Plotting style adjustments
            hUnfolded->SetLineColor(kBrown);
            hUnfolded->SetLineWidth(2);
            hUnfolded->GetYaxis()->SetTitle("Counts");
            hUnfolded->GetXaxis()->SetTitleSize(0.07);
            hUnfolded->GetYaxis()->SetTitleSize(0.07);
            hUnfolded->GetXaxis()->SetLabelSize(0.06);
            hUnfolded->GetYaxis()->SetLabelSize(0.06);
            
            hUnfolded->Draw("HIST"); // Draw unfolded histogram
            
            // Label
            latex.DrawLatexNDC(0.10, 0.86, Form("Q2-y bin: %d, z-P_{T} bin: %zu", (bin+1), (i+1)));
        }
        
        // Adjust pad margins as before
        gPad->SetLeftMargin(0.25);
        gPad->SetRightMargin(0.2);
        gPad->SetTopMargin(0.2);
        gPad->SetBottomMargin(0.2);
        }

        // Save the unfolded canvas
        unfoldedCanvas->SaveAs(Form("output/unfolded_Q2yBin_%d.png", bin+1));
        delete unfoldedCanvas;
    }

    // Cleanup
    for (int bin = 0; bin < 17; ++bin) {
        for (auto& hist : hData[bin]) delete hist;
        for (auto& hist : hMCReco[bin]) delete hist;
        for (auto& hist : hMCGene[bin]) delete hist;
    }

    fData->Close();
    fMCReco->Close();
    fMCGene->Close();
    delete fData;
    delete fMCReco;
    delete fMCGene;


    return 0;
}
