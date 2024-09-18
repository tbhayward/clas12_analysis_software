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
    double binEdges[] = {
        0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
        0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 
        1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 
        1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.20, 2.25, 2.30, 2.35, 2.4, 2.45
    };
    const int nBins = sizeof(binEdges)/sizeof(binEdges[0]) - 1; // nBins = 47

    // Prepare arrays to hold results
    vector<double> binCenters(nBins, 0.0);      // Mean Mx in each bin
    vector<int>    binCounts(nBins, 0);         // Number of events in each bin
    vector<double> delta(nBins, 0.0);
    vector<double> deltaError(nBins, 0.0);

    // For counting events
    vector<int> N_same_sign(nBins, 0);
    vector<int> N_opposite_sign(nBins, 0);

    // For accumulating Mx values to compute mean
    vector<double> totalMx(nBins, 0.0);

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

        // Accumulate Mx values and counts for mean calculation
        totalMx[binIndex] += Mx;
        binCounts[binIndex]++;

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

    // Calculate delta, errors, and bin centers (mean Mx)
    for (int i = 0; i < nBins; ++i) {
        delta[i] = N_same_sign[i] - N_opposite_sign[i];
        deltaError[i] = sqrt(N_same_sign[i] + N_opposite_sign[i]); // Assuming Poisson errors

        if (binCounts[i] > 0) {
            binCenters[i] = totalMx[i] / binCounts[i]; // Mean Mx in this bin
        } else {
            binCenters[i] = (binEdges[i] + binEdges[i+1]) / 2.0; // Default to bin midpoint
        }
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

    // Adjust the left margin to add padding space
    c1->SetLeftMargin(0.20); // Increase the left margin from default (0.1) to 0.15

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