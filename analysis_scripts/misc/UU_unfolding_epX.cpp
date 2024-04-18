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
#include <TPaveText.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TMinuit.h>


// Global pointer to histogram, so it can be accessed in the chi2 function
TH1F *hUnfoldedFilteredGlobal;

// Helper function to join vector of strings with a delimiter
std::string join(const std::vector<std::string>& vec, const std::string& delim) {
    std::string result;
    for (const auto& s : vec) {
        if (!result.empty()) result += delim;
        result += s;
    }
    return result;
}

// Fit function: a trigonometric polynomial
double fitFunction(double x, double *par) {
    return par[0] * (1 + par[1] * cos(x) + par[2] * cos(2*x));
}

// Chi-square function
void chiSquare(int &npar, double *gin, double &f, double *par, int iflag) {
    double chi2 = 0.0;
    int nbins = hUnfoldedFilteredGlobal->GetNbinsX();
    for (int i = 1; i <= nbins; i++) {
        double x = hUnfoldedFilteredGlobal->GetBinCenter(i);
        double meas = hUnfoldedFilteredGlobal->GetBinContent(i);
        double err = hUnfoldedFilteredGlobal->GetBinError(i);
        double fit = fitFunction(x, par);
        
        if (err > 0 && err < 1e10) {
            chi2 += pow((meas - fit) / err, 2);
        }
    }
    f = chi2;
}


// Function to determine the Q2-y bin based on given Q2 and y values.
int DetermineQ2yBin(float Q2, float y) {
    // Q2-y Bins 1-4
    if (Q2 > 1.000 && Q2 <= 2.000) {
        if (y > 0.650 && y <= 0.750) return 1;
        if (y > 0.550 && y <= 0.650) return 2;
        if (y > 0.450 && y <= 0.550) return 3;
        if (y > 0.300 && y <= 0.450) return 4;
    }
    // Q2-y Bins 5-8
    else if (Q2 > 2.000 && Q2 <= 3.000) {
        if (y > 0.650 && y <= 0.750) return 5;
        if (y > 0.550 && y <= 0.650) return 6;
        if (y > 0.450 && y <= 0.550) return 7;
        if (y > 0.300 && y <= 0.450) return 8;
    }
    // Q2-y Bins 9-12
    else if (Q2 > 3.000 && Q2 <= 4.000) {
        if (y > 0.650 && y <= 0.750) return 9;
        if (y > 0.550 && y <= 0.650) return 10;
        if (y > 0.450 && y <= 0.550) return 11;
        if (y > 0.300 && y <= 0.450) return 12;
    }
    // Q2-y Bins 13-15
    else if (Q2 > 4.000 && Q2 <= 5.000) {
        if (y > 0.650 && y <= 0.750) return 13;
        if (y > 0.550 && y <= 0.650) return 14;
        if (y > 0.450 && y <= 0.550) return 15;
    }
    // Q2-y Bins 16-17
    else if (Q2 > 5.000 && Q2 <= 7.000) {
        if (y > 0.650 && y <= 0.750) return 16;
        if (y > 0.550 && y <= 0.650) return 17;
    }

    return -1; // Not in any defined Q2-y bin
}

// Define bin edges for z for each Q2-y bin
std::map<int, std::vector<float>> zEdges = {
    {1, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {2, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {3, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {4, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {5, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {6, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {7, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {8, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {9, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {10, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {11, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {12, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {13, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {14, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {15, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {16, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}},
    {17, {0.10, 0.25, 0.35, 0.45, 0.55, 0.75}}
};

// Define bin edges for pT for each Q2-y bin
std::map<int, std::vector<float>> pTEdges = {
    {1, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {2, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {3, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {4, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {5, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {6, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {7, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {8, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {9, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {10, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {11, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {12, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {13, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {14, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {15, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {16, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}},
    {17, {0.00, 0.21, 0.34, 0.47, 0.60, 1.00}}
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
    // inbending files
    TFile* fData = TFile::Open("/volatile/clas12/thayward/multiplicities/data/epi+X/rga_sp19_inb_epi+X.root");
    TFile* fMCReco = TFile::Open("/volatile/clas12/thayward/multiplicities/data/epi+X/rga_sp19_inb_clasdis_50nAbkg_rec_epi+X.root");
    TFile* fMCGene = TFile::Open("/volatile/clas12/thayward/multiplicities/data/epi+X/rga_sp19_inb_clasdis_50nAbkg_gen_epi+X.root");
    // outbending files
    // TFile* fData = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_out_epX_skimmed.root");
    // TFile* fMCReco = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_out_clasdis_50nA_epX_skimmed.root");
    // TFile* fMCGene = TFile::Open("/volatile/clas12/thayward/UU_validation/root_files/rga_fa18_out_clasdis_50nA_gen_epX_skimmed.root");
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
    double e_phiData, p_phiData, e_phiMC, p_phiMC, e_phiGen, p_phiGen;
    tData->SetBranchAddress("e_phi", &e_phiData);
    tData->SetBranchAddress("p_phi", &p_phiData);
    tMCReco->SetBranchAddress("e_phi", &e_phiMC);
    tMCReco->SetBranchAddress("p_phi", &p_phiMC);
    tMCGene->SetBranchAddress("e_phi", &e_phiGen);
    tMCGene->SetBranchAddress("p_phi", &p_phiGen);
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

    // record depolarization values to turn the asymmetries into structure functions
    // record PT for each event to calculate mean PT of the bin
    struct BinParams {
        double sumDepA = 0, sumDepB = 0, sumDepV = 0, sumPT = 0;
        int count = 0; 
    };
    std::vector<std::vector<BinParams>> allBinParams;
    // After defining zEdges and pTEdges
    allBinParams.resize(zEdges.size()); // Ensure there's a vector for each Q2-y bin
    int num_z_bins[17], num_pT_bins[17];
    for (int i = 0; i <= zEdges.size()-1; ++i) { //  bins are 1-indexed=
        num_z_bins[i] = zEdges[i+1].size() - 1; // Number of z bins for this Q2-y bin
        num_pT_bins[i] = pTEdges[i+1].size() - 1; // Number of pT bins for this Q2-y bin
        int totalbins = num_z_bins[i] * num_pT_bins[i]; // Total number of z-pT bin combinations for this Q2-y bin
        allBinParams[i].resize(totalbins); // Resize the vector for this Q2-y bin to hold all combinations
    }
    std::cout << std::endl << "Creating histograms." << std::endl;
    // Create histograms for each z-pT bin
    std::vector<std::vector<TH1F*>> hData(17), hMCReco(17), hMCGene(17);
    // 17 Q2-y bins in total
    for (int bin = 0; bin < 17; ++bin) { 
        // std::cout << bin << " " << num_pT_bins[bin] << " " << num_z_bins[bin] << " " << num_pT_bins[bin] * num_z_bins[bin] << std::endl;
        for (int i = 0; i < num_pT_bins[bin] * num_z_bins[bin]; ++i) {
            hData[bin].push_back(new TH1F(Form("hData_bin%d_%d", bin+1, i), ";#phi;Normalized Counts", 24, 0, 2*TMath::Pi()));
            hMCReco[bin].push_back(new TH1F(Form("hMCReco_bin%d_%d", bin+1, i), ";#phi;Normalized Counts", 24, 0, 2*TMath::Pi()));
            hMCGene[bin].push_back(new TH1F(Form("hMCGene_bin%d_%d", bin+1, i), ";#phi;Normalized Counts", 24, 0, 2*TMath::Pi()));
        }
    }

    double phi_min = 0;
    double phi_max = 2*3.14159+0.01;
    Long64_t nDataEntries = tData->GetEntries();
    std::cout << "Looping over data. " << nDataEntries << " entries." << std::endl;
    for (Long64_t i = 0; i < nDataEntries; ++i) {
        tData->GetEntry(i);

        int binIndex = DetermineQ2yBin(Q2Data, yData) - 1; // Adjusted for 0-based indexing

        if (binIndex >= 0) {

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            // extra check to make sure data is not outside the bins
            if (zData < currentZEdges[0] || zData > currentZEdges[num_z_bins[binIndex]]) {
                continue;
            }
            if (pTData < currentPTEdges[0] || pTData > currentPTEdges[num_pT_bins[binIndex]]) {
                continue;
            }

            int z_bin = findBinIndex(zData, currentZEdges);
            int pT_bin = findBinIndex(pTData, currentPTEdges);

            // if (phiData<3.14159/3 || phiData>(2*3.14159)-3.14159/3) {
            //     continue;
            // }
            if (pT_bin != -1 && z_bin != -1 && p_phiData > phi_min && p_phiData < phi_max) {
            // if (pT_bin != -1 && z_bin != -1 && p_phiData < phi_min || p_phiData > phi_max) {
                // index of bins is flipped (per Kyungseon's request), first bin is highest z (top left)
                // "binIndex" is the Q2-y bin
                int histIndex = z_bin * num_pT_bins[binIndex] + pT_bin;
                hData[binIndex][histIndex]->Fill(phiData);
                allBinParams[binIndex][histIndex].sumDepA += DepAData;
                allBinParams[binIndex][histIndex].sumDepB += DepBData;
                allBinParams[binIndex][histIndex].sumDepV += DepVData;
                allBinParams[binIndex][histIndex].sumPT += pTData; 
                allBinParams[binIndex][histIndex].count++;
            }
        }
    }

    Long64_t nMCEntries = tMCReco->GetEntries();
    std::cout << "Looping over reconstructed MC. " << nMCEntries << " entries." << std::endl;
    for (Long64_t i = 0; i < nMCEntries; ++i) {
        tMCReco->GetEntry(i);

        int binIndex = DetermineQ2yBin(Q2MC, yMC) - 1; // Adjusted for 0-based indexing

        if (binIndex >= 0) {

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            // extra check to make sure data is not outside the bins
            if (zMC < currentZEdges[0] || zMC > currentZEdges[num_z_bins[binIndex]]) {
                continue;
            }
            if (pTMC < currentPTEdges[0] || pTMC > currentPTEdges[num_pT_bins[binIndex]]) {
                continue;
            }

            int z_bin = findBinIndex(zMC, currentZEdges);
            int pT_bin = findBinIndex(pTMC, currentPTEdges);
            // Fill the corresponding histogram if the event is in a valid bin
            // if (phiMC<3.14159/3 || phiMC>(2*3.14159)-3.14159/3) {
            //     continue;
            // }
            if (pT_bin != -1 && z_bin != -1 && p_phiMC > phi_min && p_phiMC < phi_max) {
            // if (pT_bin != -1 && z_bin != -1 && p_phiMC < phi_min || p_phiMC > phi_max) {
                // index of bins is flipped (per Kyungseon's request), first bin is highest z (top left)
                // "binIndex" is the Q2-y bin
                int histIndex = z_bin * num_pT_bins[binIndex] + pT_bin;
                hMCReco[binIndex][histIndex]->Fill(phiMC);
            }
        }
    }


    
    // Fill histograms for gen MC
    Long64_t nGenEntries = tMCGene->GetEntries();
    std::cout << "Looping over generated MC. " << nGenEntries << " entries." << std::endl;
    for (Long64_t i = 0; i < nGenEntries; ++i) {
        tMCGene->GetEntry(i);
        int binIndex = DetermineQ2yBin(Q2Gen, yGen) - 1; // Adjusted for 0-based indexing

        if (binIndex >= 0) {

            const auto& currentZEdges = zEdges[binIndex+1]; 
            const auto& currentPTEdges = pTEdges[binIndex+1];

            // extra check to make sure data is not outside the bins
            if (zGen < currentZEdges[0] || zGen > currentZEdges[num_z_bins[binIndex]]) {
                continue;
            }
            if (pTGen < currentPTEdges[0] || pTGen > currentPTEdges[num_pT_bins[binIndex]]) {
                continue;
            }

            int z_bin = findBinIndex(zGen, currentZEdges);
            int pT_bin = findBinIndex(pTGen, currentPTEdges);
            // Fill the corresponding histogram if the event is in a valid bin
            // if (phiGen<3.14159/3 || phiGen>(2*3.14159)-3.14159/3) {
            //     continue;
            // }
            if (pT_bin != -1 && z_bin != -1 && p_phiGen > phi_min && p_phiGen < phi_max) {
            // if (pT_bin != -1 && z_bin != -1 && p_phiGen < phi_min || p_phiGen > phi_max) {
                // index of bins is flipped (per Kyungseon's request), first bin is highest z (top left)
                // "binIndex" is the Q2-y bin
                int histIndex = z_bin * num_pT_bins[binIndex] + pT_bin;
                hMCGene[binIndex][histIndex]->Fill(phiGen);
            }
        }
    }

    /* ~~~~~~~~~~~~~~~~~~~~~~ */ 
    // First loop is for plotting the normalized data and reconstructed and generated Monte Carlo plots
    // Declare the TLatex object here, before the loop
    TLatex latex; latex.SetTextSize(0.09); latex.SetNDC();
    for (int bin = 0; bin < 17; ++bin) {
        TCanvas* canvas = new TCanvas(Form("canvas_bin_%d", bin + 1), Form("Q2-y Bin %d Phi Distributions", bin + 1), 2000, 1200);
        canvas->Divide(num_pT_bins[bin], num_z_bins[bin]);
        // create a canvas with subplots for each possible z-PT bin (some will be empty)

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

                // // Set the y-axis scale minimum and maximum
                // hMCGeneHist->SetMinimum(0); // For example, set minimum to 0
                // hMCGeneHist->SetMaximum(hMCGeneHist->GetMaximum() * 1.2); // Set maximum to 120% of the current maximum value for some headroom

                // Draw histograms if they have sufficient entries
                // just picked 100, probably should have done another number
                if (hMCGeneHist->GetEntries() > 400) {
                    hDataHist->DrawNormalized("HIST");
                    hMCRecoHist->DrawNormalized("HIST same");
                    hMCGeneHist->DrawNormalized("HIST same");

                    // Display z-pT bin information
                    latex.DrawLatexNDC(0.10, 0.86, Form("Q2-y bin: %d, z-PT bin: %d", bin + 1, padNumber));
                }
            }
        }

        canvas->SaveAs(Form("output3/yields/Q2yBin_%d.png", bin + 1));
        delete canvas;
    }

    // Calculate the acceptance in each of the bins
    std::vector<std::vector<TH1F*>> hAcceptance(17);
    for (int bin = 0; bin < 17; ++bin) {
        hAcceptance[bin].resize(num_pT_bins[bin] * num_z_bins[bin]);
        for (int i = 0; i < num_pT_bins[bin] * num_z_bins[bin]; ++i) {
            // Ensure MC reco histogram has entries to avoid division by zero
            if (hMCReco[bin][i]->GetEntries() > 0) {
                hAcceptance[bin][i] = (TH1F*)hMCGene[bin][i]->Clone(Form("hAcceptance_bin%d_%d", bin+1, i));
                hAcceptance[bin][i]->Divide(hMCReco[bin][i]); // Divide generated by reconstructed

                // Manually set bin errors based on Poisson statistics
                for (int binX = 1; binX <= hAcceptance[bin][i]->GetNbinsX(); ++binX) {
                    double N1 = hMCGene[bin][i]->GetBinContent(binX);
                    double N2 = hMCReco[bin][i]->GetBinContent(binX);
                    double ratio = hAcceptance[bin][i]->GetBinContent(binX);

                    // Calculate the error on the ratio
                    if (N1 > 0 && N2 > 0) {
                        double error = ratio * sqrt((1 / N1) + (1 / N2));
                        hAcceptance[bin][i]->SetBinError(binX, error);
                    } else {
                        // Set error to 0 if either N1 or N2 is 0 to avoid division by zero
                        hAcceptance[bin][i]->SetBinError(binX, 1e20);
                    }
                }
            }
        }
    }

    // Calculate the acceptance in each of the bins
    std::vector<std::vector<TH1F*>> hAcceptanceInverse(17);
    for (int bin = 0; bin < 17; ++bin) {
        hAcceptanceInverse[bin].resize(num_pT_bins[bin] * num_z_bins[bin]);
        for (int i = 0; i < num_pT_bins[bin] * num_z_bins[bin]; ++i) {
            // Ensure MC reco histogram has entries to avoid division by zero
            if (hMCReco[bin][i]->GetEntries() > 1) {
                hAcceptanceInverse[bin][i] = (TH1F*)hMCReco[bin][i]->Clone(Form("hAcceptanceInverse_bin%d_%d", bin+1, i));
                hAcceptanceInverse[bin][i]->Divide(hMCGene[bin][i]); // Divide reconstructed by generated

                // Manually set bin errors based on Poisson statistics
                for (int binX = 1; binX <= hAcceptanceInverse[bin][i]->GetNbinsX(); ++binX) {
                    double N1 = hMCReco[bin][i]->GetBinContent(binX);
                    double N2 = hMCGene[bin][i]->GetBinContent(binX);
                    double ratio = hAcceptanceInverse[bin][i]->GetBinContent(binX);

                    // Calculate the error on the ratio
                    if (N1 > 0 && N2 > 0) {
                        double error = ratio * sqrt((1 / N1) + (1 / N2));
                        hAcceptanceInverse[bin][i]->SetBinError(binX, error);
                    } else {
                        // Set error to 0 if either N1 or N2 is 0 to avoid division by zero
                        hAcceptanceInverse[bin][i]->SetBinError(binX, 0);
                    }
                }
            }
        }
    }

    // Plot the acceptance function (ratio of gen to rec MC). 
    for (int bin = 0; bin < 17; ++bin) {
        TCanvas* acceptanceCanvas = new TCanvas(Form("acceptance_canvas_bin_%d", bin + 1), Form("Acceptance Q2-y Bin %d", bin + 1), 2000, 1200);
        acceptanceCanvas->Divide(num_pT_bins[bin], num_z_bins[bin]); 

        for (int z_bin = 0; z_bin < num_z_bins[bin]; ++z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int histIndex = z_bin * num_pT_bins[bin] + pT_bin;
                int padNumber = (num_z_bins[bin] - z_bin - 1) * num_pT_bins[bin] + pT_bin + 1;
                acceptanceCanvas->cd(padNumber);

                // Set pad margins
                gPad->SetLeftMargin(0.3);
                gPad->SetRightMargin(0.1);
                gPad->SetTopMargin(0.2);
                gPad->SetBottomMargin(0.2);

                if (hAcceptanceInverse[bin][histIndex] != nullptr) {
                    TH1F* hAcc = hAcceptanceInverse[bin][histIndex];
                    hAcc->SetStats(0);

                    // // Set the y-axis scale minimum and maximum
                    hAcc->SetMinimum(0); 
                    hAcc->SetMaximum(0.1); 
                    
                    hAcc->SetMarkerStyle(20);
                    hAcc->Draw("PE"); // "PE" for drawing error bars with points

                    // Set axis titles and sizes
                    hAcc->GetXaxis()->SetTitle("#phi");
                    hAcc->GetYaxis()->SetTitle("Acceptance");
                    hAcc->GetXaxis()->SetTitleSize(0.07);
                    hAcc->GetYaxis()->SetTitleSize(0.07);
                    hAcc->GetXaxis()->SetTitleOffset(1.2);
                    hAcc->GetYaxis()->SetTitleOffset(1.2);
                    hAcc->GetXaxis()->SetLabelSize(0.06);
                    hAcc->GetYaxis()->SetLabelSize(0.06);

                    // Display z-pT bin information
                    latex.DrawLatexNDC(0.10, 0.86, Form("Q2-y bin: %d, z-PT bin: %d", bin + 1, padNumber));
                }
            }
        }
        acceptanceCanvas->SaveAs(Form("output3/acceptance/acceptance_Q2yBin_%d.png", bin + 1)); // Save each canvas
        delete acceptanceCanvas; // Clean up
    }

    // create struct to hold the values and errors of the fits as well as chi2
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
                // calculate which pad we are in (highest z bins are first)
                unfoldedCanvas->cd(padNumber);

                // Adjust pad margins to add space around the plots
                gPad->SetLeftMargin(0.3);
                gPad->SetRightMargin(0.1);
                gPad->SetTopMargin(0.2);
                gPad->SetBottomMargin(0.2);

                if (hData[bin][histIndex]->GetEntries() > 1 && hMCGene[bin][histIndex]->GetEntries() > 400 && hAcceptance[bin][histIndex] != nullptr) {
                    TH1F* hUnfolded = (TH1F*)hData[bin][histIndex]->Clone(Form("hUnfolded_bin%d_%d", bin + 1, histIndex));
                    hUnfolded->Multiply(hAcceptance[bin][histIndex]); // Apply acceptance correction

                    // calculate and set the errors for each bin
                    for (int binX = 1; binX <= hUnfolded->GetNbinsX(); ++binX) {
                        double a = hData[bin][histIndex]->GetBinContent(binX);
                        double c = hMCReco[bin][histIndex]->GetBinContent(binX);
                        double b = hAcceptance[bin][histIndex]->GetBinContent(binX);
                        double sigma_a = hData[bin][histIndex]->GetBinError(binX);
                        double sigma_b = hAcceptance[bin][histIndex]->GetBinError(binX);
                        hUnfolded->SetBinContent(binX, a*b);

                        // Calculate the uncertainty using error propagation 
                        if (a > 0 && b > 0) { // Check to avoid division by zero
                            double f = a * b; 
                            double sigma_f = f * sqrt(pow(sigma_a / a, 2) + pow(sigma_b / b, 2));
                            hUnfolded->SetBinError(binX, sigma_f);
                        } else if (a == 0 || b == 0) {
                            hUnfolded->SetBinError(binX, 1e20);
                        }

                        if (bin == 0 && padNumber == 2) { std::cout << a << " " << c << " " << b << " " << hUnfolded->GetBinContent(binX) << " " << hUnfolded->GetBinError(binX) << std::endl; }
                    }

                    // TF1* fitFunc = new TF1("fitFunc", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
                    // fitFunc->SetParameters(0, 0, 0);
                    // Threshold for acceptance
                    double acceptanceThreshold = 1/0.000000000000001; // lower number is the percentage threshold
                    // Clone the original histogram to preserve the data
                    TH1F* hUnfoldedFiltered = (TH1F*)hUnfolded->Clone("hUnfoldedFiltered");
                    // Loop over bins and only keep those with acceptance above the threshold
                    for (int binX = 0; binX <= hUnfolded->GetNbinsX(); ++binX) {
                        double acceptance = hAcceptance[bin][histIndex]->GetBinContent(binX);
                        if (acceptance > acceptanceThreshold || acceptance == 0) {
                            // For bins below the threshold, set content and error in the filtered histogram to indicate exclusion
                            hUnfoldedFiltered->SetBinContent(binX, 0);
                            hUnfoldedFiltered->SetBinError(binX, 1e20); // Set a very high error
                        }
                    }

                    TMinuit minuit(3);
                    minuit.SetPrintLevel(-1);
                    minuit.SetErrorDef(1);
                    minuit.SetFCN(chiSquare);

                    minuit.DefineParameter(0, "p0", 0.8*hUnfoldedFiltered->GetMaximum(), 0.00001, 0, 0);
                    minuit.DefineParameter(1, "p1", -0.1, 0.00001, 0, 0);
                    minuit.DefineParameter(2, "p2", 0.0, 0.00001, 0, 0);

                    hUnfoldedFilteredGlobal = hUnfoldedFiltered;

                    double arglist[10];
                    int ierflg = 0;

                    // Migrad with increased iterations and tolerance
                    arglist[0] = 50000;
                    arglist[1] = 0.01;
                    minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

                    if (ierflg == 0) {
                        // Proceed only if MIGRAD converged
                        minuit.mnexcm("HESSE", arglist, 1, ierflg);
                        if (ierflg == 0) {
                            // Proceed to MINOS only if HESSE was successful
                            minuit.mnexcm("MINOS", arglist, 1, ierflg);
                        } else {
                            std::cerr << "HESSE failed to converge." << std::endl;
                        }
                    } else {
                        std::cerr << "MIGRAD failed to converge." << std::endl;
                    }

                    // For calculating asymmetric errors, use Minos
                    // minuit.Minos();

                    double par[3], err[3];
                    for (int i = 0; i < 3; i++) {
                        minuit.GetParameter(i, par[i], err[i]);
                    }

                    // Corrected mnstat call
                    double chi2, dum1, dum2;
                    int nvpar, nparx, icstat;
                    minuit.mnstat(chi2, dum1, dum2, nvpar, nparx, icstat);
                    int num_empty_bins = 0;
                    for (int binX = 1; binX < hUnfoldedFilteredGlobal->GetNbinsX(); ++binX) {
                        if (hUnfoldedFiltered->GetBinError(binX) > 1e10) {
                            num_empty_bins++;
                        }
                    }
                    int ndf = hUnfoldedFilteredGlobal->GetNbinsX() - nvpar - num_empty_bins;

                    FitParams params = {
                        par[0], par[1], par[2],
                        err[0], err[1], err[2],
                        chi2 / ndf
                    };

                    allFitParams[bin][histIndex] = params;

                    TGraphErrors* gUnfolded = new TGraphErrors();

                    for (int binX = 1; binX <= hUnfoldedFiltered->GetNbinsX(); ++binX) {
                        // if (hUnfoldedFiltered->GetBinContent(binX) != 0) {
                            double acceptance = hAcceptance[bin][histIndex]->GetBinContent(binX);
                            if (acceptance > acceptanceThreshold || hUnfoldedFiltered->GetBinError(binX) < 1 || hUnfoldedFiltered->GetBinContent(binX) == 0) {
                                gUnfolded->SetPoint(binX - 1, 9, hUnfoldedFiltered->GetBinContent(binX));
                                gUnfolded->SetPointError(binX - 1, 0., 1e20);
                            } else {
                                gUnfolded->SetPoint(binX - 1, hUnfoldedFiltered->GetBinCenter(binX), hUnfoldedFiltered->GetBinContent(binX));
                                gUnfolded->SetPointError(binX - 1, 0., hUnfoldedFiltered->GetBinError(binX));
                            }
                        // }
                    }

                    // if (bin == 0 && padNumber == 2) {
                    //     for (int binX = 1; binX <= hUnfolded->GetNbinsX(); ++binX) {
                    //         Double_t x, y;
                    //         gUnfolded->GetPoint(binX, x, y); // Get the x and y values of the point at index i
                    //         Double_t yerr = gUnfolded->GetErrorY(binX);
                    //         std::cout << binX << " " << x << " " << y << " " << yerr << std::endl;
                    //     }
                    // }

                    double stepSize = 0.001; // Define a step size for scanning the range
                    double maxVal = -DBL_MAX; // Start with the smallest possible double
                    double minVal = DBL_MAX;  // Start with the largest possible double

                    for (double x = 0; x <= 2 * TMath::Pi(); x += stepSize) {
                        double val = fitFunction(x, par);
                        if (val > maxVal) maxVal = val;
                        if (val < minVal) minVal = val;
                    }
                    // // Find the maximum value of the function between 0 and 2*pi
                    // double maxVal = fitFunc->GetMaximum(0, 2 * TMath::Pi());
                    // double minVal = fitFunc->GetMinimum(0, 2 * TMath::Pi());

                    // Assuming the mean is still best estimated at x = 0
                    // double mean = fitFunc->Eval(0);
                    double mean = fitFunction(TMath::Pi()/2, par);

                    // Now use the maximum value to determine the amplitude
                    // Amplitude is now taken as the maximum deviation from the mean
                    double amplitude = TMath::Max(TMath::Abs(maxVal - mean), TMath::Abs(mean - minVal));

                    // Setting y-axis limits to +/- 2 times the amplitude around the mean
                    double yMin = mean - 2.0 * amplitude;
                    double yMax = mean + 2.0 * amplitude;

                    gUnfolded->SetMinimum(yMin);
                    gUnfolded->SetMaximum(yMax);

                    // Setting x-axis range to 0 to 2pi
                    gUnfolded->GetXaxis()->SetLimits(0, TMath::TwoPi());

                    gUnfolded->GetXaxis()->SetTitle("#phi");
                    gUnfolded->GetYaxis()->SetTitle("Counts");

                    gUnfolded->GetXaxis()->SetTitleSize(0.07);
                    gUnfolded->GetYaxis()->SetTitleSize(0.07);

                    gUnfolded->GetXaxis()->SetTitleOffset(1.2);
                    gUnfolded->GetYaxis()->SetTitleOffset(1.2);

                    gUnfolded->SetMarkerStyle(20);
                    gUnfolded->SetMarkerColor(kBlack);
                    gUnfolded->SetLineColor(kBlack);
                    gUnfolded->Draw("AP");

                    gUnfolded->GetXaxis()->SetLabelSize(0.06);
                    gUnfolded->GetYaxis()->SetLabelSize(0.06);
                    gUnfolded->GetXaxis()->SetTitleSize(0.07);
                    gUnfolded->GetYaxis()->SetTitleSize(0.07);

                    // fitFunc->SetLineColor(kRed);
                    // fitFunc->Draw("same");
                    // For drawing the fit result, you can create a TF1 object with the final parameters:
                    TF1* fitResult = new TF1("fitResult", "[0]*(1 + [1]*cos(x) + [2]*cos(2*x))", 0, 2*TMath::Pi());
                    fitResult->SetParameters(par);
                    fitResult->SetLineColor(kRed);
                    fitResult->Draw("same");

                    TPaveText *pt = new TPaveText(0.075, 0.0725, 0.6, 0.325 , "brNDC");
                    pt->SetBorderSize(1); // Set border size
                    pt->SetLineColor(kBlack); // Set border color
                    pt->SetFillColor(kWhite); // Set solid background color
                    pt->SetTextAlign(12); // Align text left and vertically centered
                    pt->SetTextSize(0.07); // Set text size
                    // Add lines of text
                    pt->AddText(Form("A = %.2f #pm %.1f", params.A, params.errA));
                    pt->AddText(Form("B = %.2f #pm %.3f", params.B, params.errB));
                    pt->AddText(Form("C = %.2f #pm %.3f", params.C, params.errC));
                    pt->AddText(Form("#chi^{2}/ndf = %.2f", params.chi2ndf));

                    // Draw the TPaveText
                    pt->Draw();
                    // Adjusting this display to correctly label each bin according to the new structure
                    latex.DrawLatexNDC(0.18, 0.86, Form("Q2-y bin: %d, z-PT bin: %d", bin + 1, padNumber));

                    delete hUnfolded; // Clean up
                }
            }
        }

        unfoldedCanvas->SaveAs(Form("output3/unfolded/unfolded_Q2yBin_%d.png", bin + 1));
        delete unfoldedCanvas;
    }

    std::ofstream capobiancoFile("output3/capobianco_cross_check.txt");
    for (size_t bin = 0; bin < allFitParams.size(); ++bin) {
        int current_bin = 1;
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                // Calculate the linear index based on z_bin and pT_bin
                int i = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allFitParams[bin][i];
                if (params.A != 0) { // Check if the fit was performed
                    // Print Q2-y bin heading
                    capobiancoFile << "Q2-y Bin " << bin + 1;
                    capobiancoFile << ", z-PT bin: " << current_bin
                       << ", A = " << params.A << " +/- " << params.errA
                       << ", B = " << params.B << " +/- " << params.errB
                       << ", C = " << params.C << " +/- " << params.errC
                       << ", chi2/NDF = " << params.chi2ndf
                       << ", counts = " << hMCReco[bin][i]->GetEntries() << std::endl;
                } else {
                    // If no fit was performed due to insufficient statistics
                    capobiancoFile << "Q2-y Bin " << bin + 1;
                    capobiancoFile << ", z-PT bin: " << current_bin << 
                    ", No fit performed due to insufficient statistics." << std::endl;
                }
                current_bin++;
            }
        }
    }
    capobiancoFile.close(); // Close the file after writing

    // also print out a version for myself where the asymmetries are scaled by the structure functions
    struct StructureFunction {
        double meanPT;
        double value;
        double error;
    };

    std::ofstream structureFile("output3/structure_functions.txt");
    // Loop over Q2-y bins
    for (int bin = 0; bin < allFitParams.size(); ++bin) {
        int current_bin = 1;
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allBinParams[bin][index];
                const auto& fitParams = allFitParams[bin][index];
                double meanPT = params.sumPT / params.count;
                if (fitParams.B != 0) {
                    double structureB = fitParams.B * params.sumDepA / params.sumDepV;
                    double structureC = fitParams.C * params.sumDepA / params.sumDepB;
                    double structureBerr = fitParams.errB * params.sumDepA / params.sumDepV;
                    double structureCerr = fitParams.errC * params.sumDepA / params.sumDepB;
                    structureFile << "Q2-y Bin " << (bin+1);
                    structureFile << ", z-pT Bin " << current_bin << ": "
                    << "B = {" << meanPT << ", " << structureB << ", " << structureBerr << "}, "
                    << "C = {" << meanPT << ", " << structureC << ", " << structureCerr << "}, " << std::endl;
                    current_bin++;
                } else {
                    structureFile << "Q2-y Bin " << bin + 1;
                    structureFile << ", z-pT Bin " << current_bin << ": "
                    << "B = {-1.00, 0.00, 1000.00}, "
                    << "C = {-1.00, 0.00, 1000.00}" << std::endl;
                    current_bin++;
                }
            }
        }
    }
    structureFile.close();


    std::ofstream structureFile2("output3/mathematica.txt");
    // Loop over Q2-y bins
    for (int bin = 0; bin < allFitParams.size(); ++bin) {
        structureFile2 << "sfQ2y" << (bin+1) << "B = {";
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allBinParams[bin][index];
                const auto& fitParams = allFitParams[bin][index];
                double meanPT = params.sumPT / params.count;
                if (fitParams.B != 0) {
                    double structureB = fitParams.B * params.sumDepA / params.sumDepV;
                    double structureBerr = fitParams.errB * params.sumDepA / params.sumDepV;
                    structureFile2 << "{" << meanPT << ", " << structureB << ", " << structureBerr << "}";
                } else {
                    structureFile2 << "{-1.00, 0.00, 1000.00}";
                }
                if (!(z_bin == 0 && pT_bin == num_pT_bins[bin]-1)) {
                    structureFile2 << ",";
                }
            }
        }
        structureFile2 << "};\n\n" << std::endl;
    }
    for (int bin = 0; bin < allFitParams.size(); ++bin) {
        structureFile2 << "sfQ2y" << (bin+1) << "C = {";
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allBinParams[bin][index];
                const auto& fitParams = allFitParams[bin][index];
                double meanPT = params.sumPT / params.count;
                if (fitParams.B != 0) {
                    double structureC = fitParams.C * params.sumDepA / params.sumDepB;
                    double structureCerr = fitParams.errC * params.sumDepA / params.sumDepB;
                    structureFile2 << "{" << meanPT << ", " << structureC << ", " << structureCerr << "}";
                } else {
                    structureFile2 << "{-1.00, 0.00, 1000.00}";
                }
                if (!(z_bin == 0 && pT_bin == num_pT_bins[bin]-1)) {
                    structureFile2 << ",";
                }
            }
        }
        structureFile2 << "};\n\n" << std::endl;
    }


    // Loop over Q2-y bins
    for (int bin = 0; bin < allFitParams.size(); ++bin) {
        structureFile2 << "unfQ2y" << (bin+1) << "B = {";
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allBinParams[bin][index];
                const auto& fitParams = allFitParams[bin][index];
                double meanPT = params.sumPT / params.count;
                if (fitParams.B != 0) {
                    double structureB = fitParams.B ;
                    double structureBerr = fitParams.errB;
                    structureFile2 << "{" << meanPT << ", " << structureB << ", " << structureBerr << "}";
                } else {
                    structureFile2 << "{-1.00, 0.00, 1000.00}";
                }
                if (!(z_bin == 0 && pT_bin == num_pT_bins[bin]-1)) {
                    structureFile2 << ",";
                }
            }
        }
        structureFile2 << "};\n\n" << std::endl;
    }
    for (int bin = 0; bin < allFitParams.size(); ++bin) {
        structureFile2 << "unfQ2y" << (bin+1) << "C = {";
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& params = allBinParams[bin][index];
                const auto& fitParams = allFitParams[bin][index];
                double meanPT = params.sumPT / params.count;
                if (fitParams.B != 0) {
                    double structureC = fitParams.C;
                    double structureCerr = fitParams.errC;
                    structureFile2 << "{" << meanPT << ", " << structureC << ", " << structureCerr << "}";
                } else {
                    structureFile2 << "{-1.00, 0.00, 1000.00}";
                }
                if (!(z_bin == 0 && pT_bin == num_pT_bins[bin]-1)) {
                    structureFile2 << ",";
                }
            }
        }
        structureFile2 << "};\n\n" << std::endl;
    }
    for (int bin = 0; bin < allFitParams.size(); ++bin) {
        structureFile2 << "chi2ndfQ2y" << (bin+1) << " = {";
        // Iterate through z and pT bins in the desired order
        for (int z_bin = num_z_bins[bin] - 1; z_bin >= 0; --z_bin) {
            for (int pT_bin = 0; pT_bin < num_pT_bins[bin]; ++pT_bin) {
                int index = z_bin * num_pT_bins[bin] + pT_bin;
                const auto& fitParams = allFitParams[bin][index];
                if (fitParams.B != 0) {
                    double structureC = fitParams.C;
                    double structureCerr = fitParams.errC;
                    structureFile2 << fitParams.chi2ndf;
                } else {
                    structureFile2 << "-1.00";
                }
                if (!(z_bin == 0 && pT_bin == num_pT_bins[bin]-1)) {
                    structureFile2 << ",";
                }
            }
        }
        structureFile2 << "};\n\n" << std::endl;
    }
    for (int bin = 0; bin < hAcceptanceInverse.size(); ++bin) { 
        structureFile2 << "acceptanceQ2y" << (bin+1) << " = {";
        for (int index = 0; index < hAcceptanceInverse[bin].size(); ++index) {
            structureFile2 << "{";
            if (hAcceptanceInverse[bin][index] != nullptr) { // Move check outside the loop
                // Loop through bins only if histogram exists
                for (int binX = 1; binX <= hAcceptanceInverse[bin][index]->GetNbinsX(); ++binX) { // binX should start from 1 for ROOT histograms
                    structureFile2 << hAcceptanceInverse[bin][index]->GetBinContent(binX);
                    if (binX < hAcceptanceInverse[bin][index]->GetNbinsX()) {
                        structureFile2 << ",";
                    }
                }
            } else {
                // If histogram does not exist, output placeholder for entire histogram
                structureFile2 << "0"; // Adjusted to output just one "0" as a placeholder
            }
            structureFile2 << "}";
            if (index < hAcceptanceInverse[bin].size() - 1) {
                structureFile2 << ",";
            }
        }
        structureFile2 << "};\n\n"; 
    }




    structureFile.close();
    structureFile2.close();

    fData->Close();
    fMCReco->Close();
    fMCGene->Close();
    delete fData;
    delete fMCReco;
    delete fMCGene;

    return 0;
}
