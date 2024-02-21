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

#include <map>
#include <vector>

// Define bin edges for z for each Q2-y bin
std::map<int, std::vector<float>> zEdges = {
    {1, {0.15, 0.2, 0.24, 0.29, 0.40, 0.73}},
    {2, {0.18, 0.23, 0.26, 0.31, 0.38, 0.50, 0.74}},
    {3, {0.22, 0.28, 0.35, 0.45, 0.60, 0.78}},
    {4, {0.26, 0.32, 0.37, 0.43, 0.50, 0.60}},
    {5, {0.15, 0.19, 0.24, 0.29, 0.38, 0.50, 0.73}},
    {6, {0.18, 0.23, 0.30, 0.39, 0.50, 0.78}},
    {7, {0.18, 0.23, 0.30, 0.39, 0.50, 0.78}},
    {8, {0.26, 0.32, 0.36, 0.40, 0.53, 0.72}},
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
    {5, {0.05, 0.22, 0.22, 0.32, 0.41, 0.51, 0.65, 1.00}},
    {6, {0.05, 0.22, 0.32, 0.41, 0.51, 0.51, 0.65, 1.00}},
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

    // Define the bin edges for z and pT
    float pT_edges[] = {0.05, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 1.0};
    float z_edges[] = {0.15, 0.2, 0.24, 0.29, 0.4, 0.73};
    int num_pT_bins = sizeof(pT_edges)/sizeof(float) - 1;
    int num_z_bins = sizeof(z_edges)/sizeof(float) - 1;

    struct BinParams {
        double sumDepA = 0, sumDepB = 0, sumDepV = 0, sumPT = 0;
        int count = 0; 
    };
    std::vector<std::vector<BinParams>> allBinParams(17, std::vector<BinParams>(num_pT_bins * num_z_bins));


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
        std::cout << binIndex << std::endl;
        if (binIndex >= 0) {

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];
            std::cout << zData << " " < z_bin << " " << pTData << " " << pT_bin << std::endl << std::endl;
            int z_bin = findBinIndex(zData, currentZEdges);
            int pT_bin = findBinIndex(pTData, currentPTEdges);
            // Fill the corresponding histogram if the event is in a valid bin
            if (pT_bin != -1 && z_bin != -1) {
                int histIndex = z_bin * num_pT_bins + pT_bin;
                hData[binIndex][histIndex]->Fill(phiData);
                allBinParams[binIndex][histIndex].sumDepA += DepAData;
                allBinParams[binIndex][histIndex].sumDepB += DepBData;
                allBinParams[binIndex][histIndex].sumDepV += DepVData;
                allBinParams[binIndex][histIndex].sumPT += pTData; // Assuming pTData is your pT variable
                allBinParams[binIndex][histIndex].count++;
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

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            int z_bin = findBinIndex(zMC, currentZEdges);
            int pT_bin = findBinIndex(pTMC, currentPTEdges);
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

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            int z_bin = findBinIndex(zGen, currentZEdges);
            int pT_bin = findBinIndex(pTGen, currentPTEdges);
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
            if (hMCReco[bin][i]->GetEntries() > 100) {
                hAcceptance[bin][i] = (TH1F*)hMCGene[bin][i]->Clone(Form("hAcceptance_bin%d_%d", bin+1, i));
                hAcceptance[bin][i]->Divide(hMCReco[bin][i]); // Divide generated by reconstructed
            }
        }
    }

    struct FitParams { double A, B, C, errA, errB, errC, chi2ndf; };
    std::vector<std::vector<FitParams>> allFitParams(17, std::vector<FitParams>(num_pT_bins * num_z_bins));
    for (int bin = 0; bin < 17; ++bin) {
        TCanvas* unfoldedCanvas = new TCanvas(Form("unfolded_canvas_bin_%d", bin+1), Form("Unfolded Q2-y Bin %d Phi Distributions", bin+1), 2000, 1200);
        unfoldedCanvas->Divide(num_pT_bins, num_z_bins);
        
        for (int i = 0; i < num_pT_bins * num_z_bins; ++i) {
            unfoldedCanvas->cd(i + 1);
            
            // Apply acceptance correction if both histograms have entries
            if (hData[bin][i]->GetEntries() > 100 && hAcceptance[bin][i] != nullptr) {
                TH1F* hUnfolded = (TH1F*)hData[bin][i]->Clone(Form("hUnfolded_bin%d_%d", bin+1, i));
                hUnfolded->Multiply(hAcceptance[bin][i]); // Multiply data by acceptance correction
                // Assuming gUnfolded is already properly filled and drawn
                TF1* fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
                hUnfolded->Fit(fitFunc, "Q"); // "Q" for quiet mode - doesn't print fit results to the console

                // Save fit parameters
                double A = fitFunc->GetParameter(0);
                double B = fitFunc->GetParameter(1);
                double C = fitFunc->GetParameter(2);

                // Optionally, save the errors of the parameters too
                double errA = fitFunc->GetParError(0);
                double errB = fitFunc->GetParError(1);
                double errC = fitFunc->GetParError(2);

                double chi2 = fitFunc->GetChisquare();
                double ndf = fitFunc->GetNDF();
                double chi2ndf = ndf > 0 ? chi2 / ndf : 0;

                // Inside the loop, after the fit
                FitParams params;
                params.A = A;
                params.B = B;
                params.C = C;
                params.errA = errA;
                params.errB = errB;
                params.errC = errC;
                params.chi2ndf = chi2ndf;

                allFitParams[bin][i] = params;

                // Create TGraphErrors from hUnfolded
                TGraphErrors* gUnfolded = new TGraphErrors();
                gUnfolded->SetName(Form("gUnfolded_bin%d_%d", bin+1, i));
                for (int binX = 1; binX <= hUnfolded->GetNbinsX(); ++binX) {
                    double x = hUnfolded->GetBinCenter(binX);
                    double y = hUnfolded->GetBinContent(binX);
                    double yError = hUnfolded->GetBinError(binX);
                    
                    // Set the point and its vertical error
                    gUnfolded->SetPoint(binX-1, x, y);
                    gUnfolded->SetPointError(binX-1, 0., yError); // Set 0 as the x-error
                }
                // Inside the loop, after using hUnfolded
                delete hUnfolded; // Delete the cloned histogram

                // Prepare the canvas and draw the graph
                unfoldedCanvas->cd(i + 1);
                gUnfolded->SetMarkerStyle(20); // Choose style suitable for your data points
                gUnfolded->SetMarkerColor(kBlack);
                gUnfolded->SetLineColor(kBlack);
                gUnfolded->Draw("AP"); // Draw the graph to ensure axes are created

                // Now adjust axis label and title sizes
                gUnfolded->GetXaxis()->SetLabelSize(0.06); // Adjust as needed
                gUnfolded->GetYaxis()->SetLabelSize(0.06); // Adjust as needed
                gUnfolded->GetXaxis()->SetTitleSize(0.07); // Adjust as needed
                gUnfolded->GetYaxis()->SetTitleSize(0.07); // Adjust as needed

                fitFunc->SetLineColor(kRed);
                fitFunc->Draw("same");

                // Display parameters on the plot
                TLatex latexParams;
                latexParams.SetNDC();
                latexParams.SetTextSize(0.04); // Adjust as needed
                latexParams.DrawLatex(0.375, 0.375, Form("A = %.2f #pm %.2f", A, errA));
                latexParams.DrawLatex(0.375, 0.325, Form("B = %.2f #pm %.2f", B, errB));
                latexParams.DrawLatex(0.375, 0.275, Form("C = %.2f #pm %.2f", C, errC));
                
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

    std::ofstream capobiancoFile("output/capobianco_cross_check.txt");
    for (size_t bin = 0; bin < allFitParams.size(); ++bin) {
        for (size_t i = 0; i < allFitParams[bin].size(); ++i) {
            const auto& params = allFitParams[bin][i];
            if (params.A != 0) { // Assuming -1 is the sentinel value for "fit not performed"
                capobiancoFile << "Bin " << bin+1 << ", Sub-bin " << i+1 << ": "
                               << "A = " << params.A << " +/- " << params.errA
                               << ", B = " << params.B << " +/- " << params.errB
                               << ", C = " << params.C << " +/- " << params.errC
                               << ", chi2/NDF = " << params.chi2ndf << std::endl;
            } else {
                capobiancoFile << "Q2-y bin: " << bin+1 << ", z-PT bin: " << i+1 << ": No fit performed due to insufficient statistics." << std::endl;
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
    for (size_t q2yBin = 0; q2yBin < allBinParams.size(); ++q2yBin) {
        // Assuming 5 z
        std::vector<std::vector<StructureFunction>> structureFunctionsB(5);
        std::vector<std::vector<StructureFunction>> structureFunctionsC(5);

        // Loop over all pT bins within each Q2-y bin to fill the structure functions
        for (size_t zpTBin = 0; zpTBin < allBinParams[q2yBin].size(); ++zpTBin) {
            if (allBinParams[q2yBin][zpTBin].count > 0) {
                // Calculate mean values and structure functions
                double meanDepA = allBinParams[q2yBin][zpTBin].sumDepA / allBinParams[q2yBin][zpTBin].count;
                double meanDepB = allBinParams[q2yBin][zpTBin].sumDepB / allBinParams[q2yBin][zpTBin].count;
                double meanDepV = allBinParams[q2yBin][zpTBin].sumDepV / allBinParams[q2yBin][zpTBin].count;
                double meanPT = allBinParams[q2yBin][zpTBin].sumPT / allBinParams[q2yBin][zpTBin].count;

                const auto& params = allFitParams[q2yBin][zpTBin];
                double structureB = params.B * meanDepA / meanDepV;
                double structureC = params.C * meanDepA / meanDepB;

                // Determine which z-bin this pT bin belongs to
                int zBinIndex = zpTBin / (pTEdges[q2yBin].size() - 1); 
                structureFunctionsB[zBinIndex].push_back({meanPT, structureB, params.errB});
                structureFunctionsC[zBinIndex].push_back({meanPT, structureC, params.errC});
            }
        }

        // Now, print the aggregated lists for each Q2-y and z-bin
        structureFile << "Q2-y bin: " << q2yBin + 1 << std::endl;
        for (size_t zBin = 0; zBin < 5; ++zBin) {
            structureFile << "z-bin " << zBin + 1 << " B: {";
            for (const auto& func : structureFunctionsB[zBin]) {
                structureFile << "{" << func.meanPT << ", " << func.value << ", " << func.error << "}, ";
            }
            structureFile << "}" << std::endl;

            structureFile << "z-bin " << zBin + 1 << " C: {";
            for (const auto& func : structureFunctionsC[zBin]) {
                structureFile << "{" << func.meanPT << ", " << func.value << ", " << func.error << "}, ";
            }
            structureFile << "}" << std::endl;
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
