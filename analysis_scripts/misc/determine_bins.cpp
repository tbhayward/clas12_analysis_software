// determine_bins.cpp

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

std::vector<double> find_bin_edges(std::vector<double>& values, int n_bins, double min_bin_width = 0, double max_bin_width = 0, bool enforce_max_on_last_bin = false) {
    // Ensure values are sorted
    std::sort(values.begin(), values.end());
    size_t N = values.size();
    std::vector<double> bin_edges;
    bin_edges.push_back(values.front()); // Start from min value

    size_t bin_start_idx = 0;
    for (int i = 1; i <= n_bins; ++i) {
        size_t target_idx = i * N / n_bins;
        if (target_idx >= N) target_idx = N - 1;
        double edge_value = values[target_idx];

        // Adjust edge_value to satisfy min_bin_width
        double bin_width = edge_value - bin_edges.back();
        if (min_bin_width > 0 && bin_width < min_bin_width) {
            edge_value = bin_edges.back() + min_bin_width;
            // Find the index of edge_value
            auto it = std::lower_bound(values.begin() + bin_start_idx, values.end(), edge_value);
            target_idx = it - values.begin();
            if (target_idx >= N) {
                target_idx = N - 1;
                edge_value = values[target_idx];
            } else {
                edge_value = values[target_idx];
            }
        }

        // Adjust edge_value to satisfy max_bin_width
        if (max_bin_width > 0 && bin_width > max_bin_width && (!enforce_max_on_last_bin || i < n_bins)) {
            edge_value = bin_edges.back() + max_bin_width;
            // Find the index of edge_value
            auto it = std::upper_bound(values.begin() + bin_start_idx, values.end(), edge_value);
            target_idx = it - values.begin();
            if (target_idx >= N) {
                target_idx = N - 1;
                edge_value = values[target_idx];
            } else {
                edge_value = values[target_idx];
            }
        }

        bin_edges.push_back(edge_value);
        bin_start_idx = target_idx;
    }

    // Ensure the last bin edge is the maximum value
    if (bin_edges.back() < values.back()) {
        bin_edges.back() = values.back();
    }

    return bin_edges;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Error: Please provide a ROOT file as input." << std::endl;
        return 1;
    }

    std::string filename = argv[1];

    // Open the ROOT file
    TFile* file = TFile::Open(filename.c_str(), "READ");
    if (!file || !file->IsOpen()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return 1;
    }

    // Get the TTree
    TTree* tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: Cannot find TTree PhysicsEvents in file " << filename << std::endl;
        file->Close();
        return 1;
    }

    // Set up branches
    double pT, xi, x, Q2, Mx2;
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("xi", &xi);
    tree->SetBranchAddress("x", &x);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("Mx2", &Mx2);

    // Collect values satisfying the cut
    std::vector<double> pT_values;
    std::vector<double> xi_values;
    std::vector<double> x_values;
    std::vector<double> Q2_values;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (Mx2 > 1.8225) {
            pT_values.push_back(pT);
            xi_values.push_back(xi);
            x_values.push_back(x);
            Q2_values.push_back(Q2);
        }
    }

    // Determine pT bins
    std::vector<double> pT_bin_edges = find_bin_edges(pT_values, 2, 0.3);

    // Loop over pT bins
    for (size_t i = 0; i < pT_bin_edges.size() - 1; ++i) {
        double pT_min = pT_bin_edges[i];
        double pT_max = pT_bin_edges[i + 1];

        // Collect xi values in this pT bin
        std::vector<double> xi_in_pT_bin;
        std::vector<size_t> xi_indices;
        for (size_t j = 0; j < pT_values.size(); ++j) {
            if (pT_values[j] >= pT_min && pT_values[j] < pT_max) {
                xi_in_pT_bin.push_back(xi_values[j]);
                xi_indices.push_back(j);
            }
        }

        if (xi_in_pT_bin.empty()) continue;

        // Determine xi bins
        std::vector<double> xi_bin_edges = find_bin_edges(xi_in_pT_bin, 3, 0.2);

        // Loop over xi bins
        for (size_t k = 0; k < xi_bin_edges.size() - 1; ++k) {
            double xi_min = xi_bin_edges[k];
            double xi_max = xi_bin_edges[k + 1];

            // Collect x values in this xi bin
            std::vector<double> x_in_xi_bin;
            std::vector<size_t> x_indices;
            for (size_t l = 0; l < xi_in_pT_bin.size(); ++l) {
                if (xi_in_pT_bin[l] >= xi_min && xi_in_pT_bin[l] < xi_max) {
                    size_t idx = xi_indices[l];
                    x_in_xi_bin.push_back(x_values[idx]);
                    x_indices.push_back(idx);
                }
            }

            if (x_in_xi_bin.empty()) continue;

            // Determine x bins
            std::vector<double> x_bin_edges = find_bin_edges(x_in_xi_bin, 4, 0, 0.12, true);

            // Loop over x bins
            for (size_t m = 0; m < x_bin_edges.size() - 1; ++m) {
                double x_min = x_bin_edges[m];
                double x_max = x_bin_edges[m + 1];

                // Collect Q2 values in this x bin
                std::vector<double> Q2_in_x_bin;
                for (size_t n = 0; n < x_in_xi_bin.size(); ++n) {
                    if (x_in_xi_bin[n] >= x_min && x_in_xi_bin[n] < x_max) {
                        size_t idx = x_indices[n];
                        Q2_in_x_bin.push_back(Q2_values[idx]);
                    }
                }

                if (Q2_in_x_bin.empty()) continue;

                // Determine Q2 bins
                std::vector<double> Q2_bin_edges = find_bin_edges(Q2_in_x_bin, 3);

                // Loop over Q2 bins
                for (size_t o = 0; o < Q2_bin_edges.size() - 1; ++o) {
                    double Q2_min = Q2_bin_edges[o];
                    double Q2_max = Q2_bin_edges[o + 1];

                    // Print out the bin bounds in the required format
                    std::cout << pT_min << " <= *pT && " << pT_max << " > *pT && ";
                    std::cout << xi_min << " <= *xi && " << xi_max << " > *xi && ";
                    std::cout << x_min << " <= *x && " << x_max << " > *x && ";
                    std::cout << "*Q2 > " << Q2_min << " && *Q2 < " << Q2_max << std::endl;
                }
            }
        }
    }

    file->Close();
    return 0;
}