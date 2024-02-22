#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <iostream>
#include <TLatex.h>
#include <TGraphErrors.h> 
#include <TF1.h>
#include <fstream> 
#include <vector>
#include <map>

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

// Define bin edges for z for each Q2-y bin
std::map<int, std::vector<float>> zEdges = {
    {1, {0.15, 0.2, 0.24, 0.29, 0.40, 0.73}},
    {2, {0.18, 0.23, 0.26, 0.31, 0.38, 0.50, 0.74}},
    {3, {0.22, 0.28, 0.35, 0.45, 0.60, 0.78}},
    {4, {0.26, 0.32, 0.37, 0.43, 0.50, 0.60, 0.71}},
    {5, {0.15, 0.19, 0.24, 0.29, 0.38, 0.50, 0.73}},
    {6, {0.18, 0.23, 0.30, 0.39, 0.50, 0.78}},
    {7, {0.18, 0.23, 0.30, 0.39, 0.50, 0.78}},
    {8, {0.26, 0.32, 0.36, 0.40, 0.45, 0.53, 0.72}},
    {9, {0.15, 0.20, 0.24, 0.30, 0.38, 0.48, 0.72}},
    {10, {0.18, 0.23, 0.26, 0.32, 0.40, 0.50, 0.72}},
    {11, {0.21, 0.26, 0.32, 0.40, 0.50, 0.70}},
    {12, {0.26, 0.32, 0.40, 0.50, 0.70}},
    {13, {0.15, 0.20, 0.24, 0.30, 0.40, 0.72}},
    {14, {0.18, 0.23, 0.27, 0.33, 0.44, 0.74}},
    {15, {0.21, 0.28, 0.35, 0.47, 0.72}},
    {16, {0.15, 0.20, 0.25, 0.32, 0.41, 0.71}},
    {17, {0.18, 0.23, 0.30, 0.38, 0.48, 0.72}}
};

// Define bin edges for pT for each Q2-y bin
std::map<int, std::vector<float>> pTEdges = {
    {1, {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00}},
    {2, {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00}},
    {3, {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00}},
    {4, {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80}},
    {5, {0.05, 0.22, 0.32, 0.41, 0.51, 0.65, 1.00}},
    {6, {0.05, 0.22, 0.32, 0.41, 0.51, 0.65, 1.00}},
    {7, {0.05, 0.20, 0.30, 0.40, 0.50, 0.65, 1.00}},
    {8, {0.05, 0.20, 0.30, 0.40, 0.52, 0.75}},
    {9, {0.05, 0.22, 0.30, 0.38, 0.46, 0.60, 0.95}},
    {10, {0.05, 0.22, 0.32, 0.41, 0.51, 0.65, 1.00}},
    {11, {0.05, 0.20, 0.31, 0.40, 0.50, 0.64, 0.95}},
    {12, {0.05, 0.22, 0.32, 0.41, 0.51, 0.67}},
    {13, {0.05, 0.23, 0.34, 0.43, 0.55, 0.90}},
    {14, {0.05, 0.23, 0.34, 0.44, 0.55, 0.90}},
    {15, {0.05, 0.23, 0.34, 0.45, 0.58, 0.90}},
    {16, {0.05, 0.24, 0.36, 0.55, 0.80}},
    {17, {0.05, 0.23, 0.36, 0.51, 0.85}}
};


// Function to find bin index given value and bin edges
int findBinIndex(float value, const std::vector<float>& edges) {
    for (size_t i = 0; i < edges.size() - 1; i++) {
        if (value > edges[i] && value <= edges[i + 1]) {
            return i;
        }
    }
    return -1; // Value outside bin edges
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
    double Q2Data, yData, phiData, pTData, zData, DepAData, DepBData, DepVData;
    double Q2MC, yMC, phiMC, pTMC, zMC;
    double Q2Gen, yGen, phiGen, pTGen, zGen;
    tData->SetBranchAddress("Q2", &Q2Data);
    tData->SetBranchAddress("y", &yData);
    tData->SetBranchAddress("phi", &phiData);
    tData->SetBranchAddress("pT", &pTData);
    tData->SetBranchAddress("z", &zData);
    tData->SetBranchAddress("DepA", &DepAData);
    tData->SetBranchAddress("DepB", &DepBData);
    tData->SetBranchAddress("DepV", &DepVData);
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

    struct BinParams {
        double sumDepA = 0, sumDepB = 0, sumDepV = 0, sumPT = 0;
        int count = 0; 
    };
    std::vector<std::vector<BinParams>> allBinParams;
    // After defining zEdges and pTEdges
    allBinParams.resize(zEdges.size()); // Ensure there's a vector for each Q2-y bin
    int num_z_bins[17], num_pT_bins[17];
    for (int i = 0; i <= zEdges.size()-1; ++i) { // Assuming bins are 1-indexed based on your map
        num_z_bins[i] = zEdges[i+1].size() - 1; // Number of z bins for this Q2-y bin
        num_pT_bins[i] = pTEdges[i+1].size() - 1; // Number of pT bins for this Q2-y bin
        int totalbins = num_z_bins[i] * num_pT_bins[i]; // Total number of z-pT bin combinations for this Q2-y bin
        allBinParams[i].resize(totalbins); // Resize the vector for this Q2-y bin to hold all combinations
    }
    std::cout << std::endl << "Creating histograms." << std::endl;
    // Create histograms for each z-pT bin
    std::vector<std::vector<TH1F*>> hData(17), hMCReco(17), hMCGene(17);
    for (int bin = 0; bin < 17; ++bin) {
        // std::cout << bin << " " << num_pT_bins[bin] << " " << num_z_bins[bin] << " " << num_pT_bins[bin] * num_z_bins[bin] << std::endl;
        for (int i = 0; i < num_pT_bins[bin] * num_z_bins[bin]; ++i) {
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
        // std::cout << binIndex << std::endl;

        if (binIndex >= 0) {

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            if (zData < currentZEdges[0] || zData > currentZEdges[num_z_bins[binIndex]]) {
                continue;
            }
            if (pTData < currentPTEdges[0] || pTData > currentPTEdges[num_pT_bins[binIndex]]) {
                continue;
            }

            int z_bin = findBinIndex(zData, currentZEdges);
            int pT_bin = findBinIndex(pTData, currentPTEdges);

            if (pT_bin != -1 && z_bin != -1) {
                int histIndex = z_bin * num_pT_bins[binIndex] + pT_bin;
                // std::cout << "z: " << zData;
                // std::cout << ", z_bin: " << z_bin;
                // std::cout << ", pT: " << pTData;
                // std::cout << ", pT_bin: " << pT_bin;
                // std::cout << ", histIndex " << histIndex << std::endl;
                hData[binIndex][histIndex]->Fill(phiData);
                allBinParams[binIndex][histIndex].sumDepA += DepAData;
                allBinParams[binIndex][histIndex].sumDepB += DepBData;
                allBinParams[binIndex][histIndex].sumDepV += DepVData;
                allBinParams[binIndex][histIndex].sumPT += pTData; 
                allBinParams[binIndex][histIndex].count++;
            }
        }
    }

    std::cout << "Looping over reconstructed MC." << std::endl;
    // Fill histograms for MC
    Long64_t nMCEntries = tMCReco->GetEntries();
    for (Long64_t i = 0; i < nMCEntries; ++i) {
        tMCReco->GetEntry(i);

        int binIndex = DetermineQ2yBin(Q2MC, yMC) - 1; // Adjusted for 0-based indexing

        if (binIndex >= 0) {

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            if (zMC < currentZEdges[0] || zMC > currentZEdges[num_z_bins[binIndex]]) {
                continue;
            }
            if (pTMC < currentPTEdges[0] || pTMC > currentPTEdges[num_pT_bins[binIndex]]) {
                continue;
            }

            int z_bin = findBinIndex(zMC, currentZEdges);
            int pT_bin = findBinIndex(pTMC, currentPTEdges);
            // Fill the corresponding histogram if the event is in a valid bin
            if (pT_bin != -1 && z_bin != -1) {
                int histIndex = z_bin * num_pT_bins[binIndex] + pT_bin;
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

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            if (zGen < currentZEdges[0] || zGen > currentZEdges[num_z_bins[binIndex]]) {
                continue;
            }
            if (pTGen < currentPTEdges[0] || pTGen > currentPTEdges[num_pT_bins[binIndex]]) {
                continue;
            }

            int z_bin = findBinIndex(zGen, currentZEdges);
            int pT_bin = findBinIndex(pTGen, currentPTEdges);
            // Fill the corresponding histogram if the event is in a valid bin
            if (pT_bin != -1 && z_bin != -1) {
                int histIndex = z_bin * num_pT_bins[binIndex] + pT_bin;
                hMCGene[binIndex][histIndex]->Fill(phiGen);
            }
        }
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~ */ 
    // Declare the TLatex object here, before the loop
    TLatex latex; latex.SetTextSize(0.09); latex.SetNDC();
    for (int bin = 0; bin < 17; ++bin) {
        TCanvas* canvas = new TCanvas(Form("canvas_bin_%d", bin + 1), Form("Q2-y Bin %d Phi Distributions", bin + 1), 2000, 1200);
        canvas->Divide(num_pT_bins[bin], num_z_bins[bin]);

        // Adjust loop to iterate over z and pT bins separately
        for (int z_bin = 0; z_bin < num_z_bins[bin]; ++z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int histIndex = z_bin * num_pT_bins[bin] + pT_bin;
                // Calculate pad number with the highest z bins at the top and pT increasing left to right
                int padNumber = (num_z_bins[bin] - z_bin - 1) * num_pT_bins[bin] + pT_bin + 1;
                canvas->cd(padNumber);

                // Adjust pad margins to add space around the plots
                gPad->SetLeftMargin(0.25);
                gPad->SetRightMargin(0.2);
                gPad->SetTopMargin(0.2);
                gPad->SetBottomMargin(0.2);

                // Setup for hData histograms
                TH1F* hDataHist = hData[bin][histIndex];
                hDataHist->SetStats(0); // Remove the stat box
                hDataHist->SetLineColor(kBlue + 2);
                hDataHist->SetLineWidth(2);
                hDataHist->GetXaxis()->SetLabelSize(0.06); 
                hDataHist->GetYaxis()->SetLabelSize(0.06); 
                hDataHist->GetXaxis()->SetTitleSize(0.07);
                hDataHist->GetYaxis()->SetTitleSize(0.07); 

                // Setup for hMCReco histograms
                TH1F* hMCRecoHist = hMCReco[bin][histIndex];
                hMCRecoHist->SetLineColor(kRed + 2);
                hMCRecoHist->SetLineWidth(2);
                hMCRecoHist->GetXaxis()->SetLabelSize(0.06); 
                hMCRecoHist->GetYaxis()->SetLabelSize(0.06); 
                hMCRecoHist->GetXaxis()->SetTitleSize(0.07);
                hMCRecoHist->GetYaxis()->SetTitleSize(0.07);

                // Setup for hMCGene histograms
                TH1F* hMCGeneHist = hMCGene[bin][histIndex];
                hMCGeneHist->SetLineColor(kGreen + 2);
                hMCGeneHist->SetLineWidth(2);
                hMCGeneHist->GetXaxis()->SetLabelSize(0.06); 
                hMCGeneHist->GetYaxis()->SetLabelSize(0.06); 
                hMCGeneHist->GetXaxis()->SetTitleSize(0.07);
                hMCGeneHist->GetYaxis()->SetTitleSize(0.07);

                // Draw histograms if they have sufficient entries
                if (hDataHist->GetEntries() > 500) {
                    hDataHist->DrawNormalized("HIST");
                    if (hMCRecoHist->GetEntries() > 500) {
                        hMCRecoHist->DrawNormalized("HIST same");
                    }
                    if (hMCGeneHist->GetEntries() > 500) {
                        hMCGeneHist->DrawNormalized("HIST same");
                    }

                    // Display z-pT bin information
                    latex.DrawLatexNDC(0.10, 0.86, Form("Q2-y bin: %d, z-P_{T} bin: z%d-pT%d", bin + 1, num_z_bins[bin] - z_bin, pT_bin + 1));
                }
            }
        }

        canvas->SaveAs(Form("output/Q2yBin_%d.png", bin + 1));
        delete canvas;
    }



    std::vector<std::vector<TH1F*>> hAcceptance(17);
    for (int bin = 0; bin < 17; ++bin) {
        hAcceptance[bin].resize(num_pT_bins[bin] * num_z_bins[bin]);
        for (int i = 0; i < num_pT_bins[bin] * num_z_bins[bin]; ++i) {
            // Ensure MC reco histogram has entries to avoid division by zero
            if (hMCReco[bin][i]->GetEntries() > 100) {
                hAcceptance[bin][i] = (TH1F*)hMCGene[bin][i]->Clone(Form("hAcceptance_bin%d_%d", bin+1, i));
                hAcceptance[bin][i]->Divide(hMCReco[bin][i]); // Divide generated by reconstructed
            }
        }
    }

    struct FitParams { double A, B, C, errA, errB, errC, chi2ndf; };
    std::vector<std::vector<FitParams>> allFitParams(zEdges.size());

    for (int bin = 0; bin < zEdges.size(); ++bin) {
        // Corrected indexing for num_z_bins and num_pT_bins arrays
        int totalBins = num_z_bins[bin] * num_pT_bins[bin];
        allFitParams[bin].resize(totalBins);
    }

    for (int bin = 0; bin < 17; ++bin) {
        TCanvas* unfoldedCanvas = new TCanvas(Form("unfolded_canvas_bin_%d", bin + 1), Form("Unfolded Q2-y Bin %d Phi Distributions", bin + 1), 2000, 1200);
        unfoldedCanvas->Divide(num_pT_bins[bin], num_z_bins[bin]);

        for (int z_bin = 0; z_bin < num_z_bins[bin]; ++z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int histIndex = z_bin * num_pT_bins[bin] + pT_bin;
                int padNumber = (num_z_bins[bin] - z_bin - 1) * num_pT_bins[bin] + pT_bin + 1;
                unfoldedCanvas->cd(padNumber);

                if (hData[bin][histIndex]->GetEntries() > 100 && hAcceptance[bin][histIndex] != nullptr) {
                    TH1F* hUnfolded = (TH1F*)hData[bin][histIndex]->Clone(Form("hUnfolded_bin%d_%d", bin + 1, histIndex));
                    hUnfolded->Multiply(hAcceptance[bin][histIndex]); // Apply acceptance correction

                    TF1* fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
                    hUnfolded->Fit(fitFunc, "Q"); // Quiet mode fit

                    // Extracting fit parameters and their errors
                    FitParams params = {
                        fitFunc->GetParameter(0),
                        fitFunc->GetParameter(1),
                        fitFunc->GetParameter(2),
                        fitFunc->GetParError(0),
                        fitFunc->GetParError(1),
                        fitFunc->GetParError(2),
                        fitFunc->GetChisquare() / fitFunc->GetNDF()
                    };
                    allFitParams[bin][histIndex] = params;

                    TGraphErrors* gUnfolded = new TGraphErrors();
                    for (int binX = 1; binX <= hUnfolded->GetNbinsX(); ++binX) {
                        gUnfolded->SetPoint(binX - 1, hUnfolded->GetBinCenter(binX), hUnfolded->GetBinContent(binX));
                        gUnfolded->SetPointError(binX - 1, 0., hUnfolded->GetBinError(binX));
                    }

                    gUnfolded->SetMarkerStyle(20);
                    gUnfolded->SetMarkerColor(kBlack);
                    gUnfolded->SetLineColor(kBlack);
                    gUnfolded->Draw("AP");

                    gUnfolded->GetXaxis()->SetLabelSize(0.06);
                    gUnfolded->GetYaxis()->SetLabelSize(0.06);
                    gUnfolded->GetXaxis()->SetTitleSize(0.07);
                    gUnfolded->GetYaxis()->SetTitleSize(0.07);

                    fitFunc->SetLineColor(kRed);
                    fitFunc->Draw("same");

                    TLatex latexParams;
                    latexParams.SetNDC();
                    latexParams.SetTextSize(0.08);
                    // latexParams.DrawLatex(0.3, 0.375, Form("A = %.2f #pm %.2f", params.A, params.errA));
                    latexParams.DrawLatex(0.3, 0.335, Form("B = %.2f #pm %.3f", params.B, params.errB));
                    latexParams.DrawLatex(0.3, 0.275, Form("C = %.2f #pm %.3f", params.C, params.errC));
                    latexParams.DrawLatex(0.3, 0.215, Form("#chi^{2}/ndf = %.2f", params.chi2ndf));

                    // Adjusting this display to correctly label each bin according to your new structure
                    latex.DrawLatexNDC(0.14, 0.86, Form("Q2-y bin: %d, z-P_{T} bin: z%d-pT%d", bin + 1, num_z_bins[bin] - z_bin, pT_bin + 1));

                    delete hUnfolded; // Clean up
                }
            }
        }

        unfoldedCanvas->SaveAs(Form("output/unfolded_Q2yBin_%d.png", bin + 1));
        delete unfoldedCanvas;
    }

    std::ofstream capobiancoFile("output/capobianco_cross_check.txt");
    for (size_t bin = 0; bin < allFitParams.size(); ++bin) {
        int current_bin = 1;
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                // Calculate the linear index based on z_bin and pT_bin
                int i = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allFitParams[bin][i];
                if (params.A != 0) { // Check if the fit was performed
                    // Print Q2-y bin heading
                    capobiancoFile << "Q2-y Bin " << bin + 1 << std::endl;
                    capobiancoFile << "z-PT bin: " << current_bin
                       << ", A = " << params.A << " +/- " << params.errA
                       << ", B = " << params.B << " +/- " << params.errB
                       << ", C = " << params.C << " +/- " << params.errC
                       << ", chi2/NDF = " << params.chi2ndf << std::endl;
                } else {
                    // If no fit was performed due to insufficient statistics
                    capobiancoFile << "Sub-bin (z" << num_z_bins[bin] - z_bin << "-pT" << pT_bin + 1 << "): No fit performed due to insufficient statistics." << std::endl;
                }
                current_bin++;
            }
        }
    }
    capobiancoFile.close(); // Close the file after writing

    struct StructureFunction {
        double meanPT;
        double value;
        double error;
    };

    std::ofstream structureFile("output/structure_functions.txt");

    // Loop over Q2-y bins
    for (int bin = 0; bin < 17; ++bin) {

        int current_bin = 1;
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& binParams = allBinParams[bin][index];
                if (binParams.count > 0) {
                    double meanPT = binParams.sumPT / binParams.count;
                    // Assuming the existence of an allFitParams similar structure to hold fit parameters for B and C
                    const auto& fitParams = allFitParams[bin][index];
                    double structureB = fitParams.B * binParams.sumDepA / binParams.sumDepV;
                    double structureC = fitParams.C * binParams.sumDepA / binParams.sumDepB;
                    double structureBerr = fitParams.errB * binParams.sumDepA / binParams.sumDepV;
                    double structureCerr = fitParams.errC * binParams.sumDepA / binParams.sumDepB;
                    structureFile << "Q2-y Bin " << bin + 1 << std::endl;
                    structureFile << "z-pT Bin " << current_bin << ": "
                    << "B = {" << meanPT << ", " << structureB << ", " << structureBerr << "}, "
                    << "C = {" << meanPT << ", " << structureC << ", " << structureCerr << "}, " << std::endl;
                    current_bin++;
                }
            }
        }
    }

    structureFile.close();


    fData->Close();
    fMCReco->Close();
    fMCGene->Close();
    delete fData;
    delete fMCReco;
    delete fMCGene;


    return 0;
}
