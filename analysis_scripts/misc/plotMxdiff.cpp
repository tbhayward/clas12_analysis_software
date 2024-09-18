#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"

using namespace std;

int main(int argc, char* argv[]) {
    // Check for input arguments
    if (argc != 2) {
        cout << "Usage: " << argv[0] << " input_file.root" << endl;
        return 1;
    }

    string inputFilename = argv[1];

    // Open the ROOT file
    TFile* inputFile = TFile::Open(inputFilename.c_str(), "READ");
    if (!inputFile || !inputFile->IsOpen()) {
        cout << "Error: Cannot open input file " << inputFilename << endl;
        return 1;
    }

    // Get the tree
    TTree* tree = (TTree*)inputFile->Get("PhysicsEvents");
    if (!tree) {
        cout << "Error: Cannot find tree 'PhysicsEvents' in file " << inputFilename << endl;
        return 1;
    }

    // Define variables to read from the tree
    Double_t Mx, Q2, W, y, target_pol;
    Int_t helicity;

    // Set branch addresses
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("helicity", &helicity);
    tree->SetBranchAddress("target_pol", &target_pol);

    // Define bin edges
    const int nBins = 21; // Number of bins
    double binEdges[nBins+1] = {
        0, 0.1, 0.2, 0.3, 0.4, 0.5,
        0.6, 0.7, 0.8, 0.9, 1.0,
        1.1, 1.2, 1.3, 1.4, 1.6,
        1.8, 2.0, 2.2, 2.4, 2.8
    };

    // Prepare arrays to hold results
    vector<double> binCenters(nBins);
    vector<double> delta(nBins, 0);
    vector<double> deltaError(nBins, 0);

    // For counting events
    vector<int> N_same_sign(nBins, 0);
    vector<int> N_opposite_sign(nBins, 0);

    // Get number of entries
    Long64_t nEntries = tree->GetEntries();

    // Loop over tree entries
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Apply selection criteria
        if (Q2 <= 1 || W <= 2 || y >= 0.80) continue;

        // Find the bin index for Mx
        int binIndex = -1;
        for (int j = 0; j < nBins; ++j) {
            if (Mx >= binEdges[j] && Mx < binEdges[j+1]) {
                binIndex = j;
                break;
            }
        }
        if (binIndex == -1) continue; // Mx out of range

        // Check for zero values
        if (helicity == 0 || target_pol == 0.0) continue;

        // Check signs of helicity and target_pol
        if ((helicity > 0 && target_pol > 0.0) || (helicity < 0 && target_pol < 0.0)) {
            // Same sign
            N_same_sign[binIndex]++;
        } else if ((helicity > 0 && target_pol < 0.0) || (helicity < 0 && target_pol > 0.0)) {
            // Opposite sign
            N_opposite_sign[binIndex]++;
        } else {
            // One of them is zero or some unexpected value, ignore
            continue;
        }
    }

    // Calculate delta and errors
    for (int i = 0; i < nBins; ++i) {
        delta[i] = N_same_sign[i] - N_opposite_sign[i];
        deltaError[i] = sqrt(N_same_sign[i] + N_opposite_sign[i]); // Assuming Poisson errors
        binCenters[i] = (binEdges[i] + binEdges[i+1]) / 2.0;
    }

    // Create TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nBins, &binCenters[0], &delta[0], 0, &deltaError[0]);

    graph->SetTitle("");

    // Set axis labels
    graph->GetXaxis()->SetTitle("M_{x} (GeV)");
    graph->GetYaxis()->SetTitle("#Delta (++ + --) - (+- + -+)");

    // Set marker style
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1);
    graph->SetLineWidth(2);

    // Create canvas and draw graph
    TCanvas* c1 = new TCanvas("c1", "Delta vs Mx", 800, 600);
    graph->Draw("AP");

    // Save the canvas as PNG
    system("mkdir -p output");
    c1->SaveAs("output/plotMxdiff.png");

    // Clean up
    delete c1;
    delete graph;
    inputFile->Close();
    delete inputFile;

    return 0;
}