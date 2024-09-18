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
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TMultiGraph.h"

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
    Double_t Mx, Q2, W, y, target_pol, x, t, tmin;
    Int_t helicity;

    // Set branch addresses
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("tmin", &tmin);
    tree->SetBranchAddress("helicity", &helicity);
    tree->SetBranchAddress("target_pol", &target_pol);

    // Define bin edges
    double binEdges[] = {
        0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 
        0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 
        1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 
        1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45
    };
    const int nBins = sizeof(binEdges)/sizeof(binEdges[0]) - 1; // nBins = 49

    // Number of subplots and datasets
    const int nSubplots = 3;
    const int nDatasets = 2;

    // x ranges for each subplot
    double xRanges[nSubplots][2] = {
        {0.16, 0.20},
        {0.28, 0.32},
        {0.35, 0.40}
    };

    // x_B midpoints for each subplot
    double xB[nSubplots] = {0.18, 0.30, 0.375}; // Midpoints of each x range

    // t - tmin cuts for each dataset
    // Dataset 0: t - tmin > -1 (blue circles)
    // Dataset 1: t - tmin < -2 (red squares)

    // Prepare arrays to hold results
    vector<vector<vector<double>>> binCenters(nSubplots, vector<vector<double>>(nDatasets, vector<double>(nBins, 0.0))); // Mean Mx
    vector<vector<vector<int>>> binCounts(nSubplots, vector<vector<int>>(nDatasets, vector<int>(nBins, 0))); // Counts
    vector<vector<vector<double>>> delta(nSubplots, vector<vector<double>>(nDatasets, vector<double>(nBins, 0.0)));
    vector<vector<vector<double>>> deltaError(nSubplots, vector<vector<double>>(nDatasets, vector<double>(nBins, 0.0)));

    vector<vector<vector<int>>> N_same_sign(nSubplots, vector<vector<int>>(nDatasets, vector<int>(nBins, 0)));
    vector<vector<vector<int>>> N_opposite_sign(nSubplots, vector<vector<int>>(nDatasets, vector<int>(nBins, 0)));
    vector<vector<vector<double>>> totalMx(nSubplots, vector<vector<double>>(nDatasets, vector<double>(nBins, 0.0)));

    // Get number of entries
    Long64_t nEntries = tree->GetEntries();

    // Loop over tree entries
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Apply selection criteria
        if (Q2 <= 1 || W <= 2 || y >= 0.80) continue;

        // Loop over subplots
        for (int s = 0; s < nSubplots; ++s) {
            // Check x range for this subplot
            if (x < xRanges[s][0] || x > xRanges[s][1]) continue;

            // For each dataset
            for (int d = 0; d < nDatasets; ++d) {
                // Apply t - tmin cut
                bool t_tmin_condition = false;
                if (d == 0) {
                    t_tmin_condition = (t - tmin) > -1;
                } else if (d == 1) {
                    t_tmin_condition = (t - tmin) < -2;
                }

                if (!t_tmin_condition) continue;

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
                totalMx[s][d][binIndex] += Mx;
                binCounts[s][d][binIndex]++;

                // Check signs of helicity and target_pol
                if ((helicity > 0 && target_pol > 0.0) || (helicity < 0 && target_pol < 0.0)) {
                    // Same sign
                    N_same_sign[s][d][binIndex]++;
                } else if ((helicity > 0 && target_pol < 0.0) || (helicity < 0 && target_pol > 0.0)) {
                    // Opposite sign
                    N_opposite_sign[s][d][binIndex]++;
                }
            } // End of dataset loop
        } // End of subplot loop
    } // End of event loop

    // Calculate delta, errors, and bin centers (mean Mx)
    for (int s = 0; s < nSubplots; ++s) {
        for (int d = 0; d < nDatasets; ++d) {
            for (int i = 0; i < nBins; ++i) {
                delta[s][d][i] = N_same_sign[s][d][i] - N_opposite_sign[s][d][i];
                deltaError[s][d][i] = sqrt(N_same_sign[s][d][i] + N_opposite_sign[s][d][i]); // Assuming Poisson errors

                if (binCounts[s][d][i] > 0) {
                    binCenters[s][d][i] = totalMx[s][d][i] / binCounts[s][d][i]; // Mean Mx in this bin
                } else {
                    binCenters[s][d][i] = (binEdges[i] + binEdges[i+1]) / 2.0; // Default to bin midpoint
                }
            }
        }
    }

    // Create canvas and divide into 1x3 subplots
    TCanvas* c1 = new TCanvas("c1", "Delta vs Mx", 2400, 600); // Width increased to accommodate 3 plots
    c1->Divide(3,1);

    // Colors and markers for datasets
    int colors[nDatasets] = {kBlue, kRed};
    int markers[nDatasets] = {20, 21}; // Circle and square

    // Legends for t - tmin cuts
    const char* t_tmin_labels[nDatasets] = {"t - t_{min} > -1", "t - t_{min} < -2"};

    // Loop over subplots to create and draw graphs
    for (int s = 0; s < nSubplots; ++s) {
        c1->cd(s+1);

        // Adjust the left margin for the first subplot
        if (s == 0) gPad->SetLeftMargin(0.20);
        else gPad->SetLeftMargin(0.10); // Smaller margins for other subplots

        // Adjust bottom margin
        gPad->SetBottomMargin(0.15);

        // Create multigraph to hold both datasets
        TMultiGraph* mg = new TMultiGraph();

        // Create legend
        TLegend* legend = new TLegend(0.15, 0.75, 0.40, 0.90); // Move to the top left
        legend->SetBorderSize(1); // Add black box border
        legend->SetTextSize(0.04);

        // Add entry for x_B (for each subplot)
        legend->AddEntry((TObject*)0, Form("x_{B} = %.3f", xB[s]), "");

        for (int d = 0; d < nDatasets; ++d) {
            // Prepare data vectors
            vector<double> x_vals;
            vector<double> y_vals;
            vector<double> x_errs;
            vector<double> y_errs;

            for (int i = 0; i < nBins; ++i) {
                if (binCounts[s][d][i] > 0) {
                    x_vals.push_back(binCenters[s][d][i]);
                    y_vals.push_back(delta[s][d][i]);
                    x_errs.push_back(0);
                    y_errs.push_back(deltaError[s][d][i]);
                }
            }

            // Create TGraphErrors
            int nPoints = x_vals.size();
            if (nPoints == 0) continue; // No data for this dataset

            TGraphErrors* graph = new TGraphErrors(nPoints, &x_vals[0], &y_vals[0], &x_errs[0], &y_errs[0]);

            graph->SetMarkerStyle(markers[d]);
            graph->SetMarkerSize(1);
            graph->SetLineWidth(1);
            graph->SetMarkerColor(colors[d]);
            graph->SetLineColor(colors[d]);

            mg->Add(graph);

            // Add entry to legend for t - tmin cuts
            legend->AddEntry(graph, t_tmin_labels[d], "p");
        }

        mg->Draw("AP");

        // Set axis labels
        mg->GetXaxis()->SetTitle("M_{x} (GeV)");
        mg->GetYaxis()->SetTitle("#Delta (++ + --) - (+- + -+)");
        mg->GetXaxis()->SetTitleSize(0.05);
        mg->GetYaxis()->SetTitleSize(0.05);
        mg->GetXaxis()->SetTitleOffset(1.0);
        mg->GetYaxis()->SetTitleOffset(1.5);

        // Draw legend
        legend->Draw();

        // Add text for x range
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.15, 0.92, Form("%.2f < x < %.2f", xRanges[s][0], xRanges[s][1]));
    }

    // Save the canvas as PNG
    system("mkdir -p output/epX_plots");
    c1->SaveAs("output/epX_plots/plotMxdiff.png");

    // Clean up
    delete c1;
    inputFile->Close();
    delete inputFile;

    return 0;
}