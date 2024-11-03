// determine_bins.cpp

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include <iomanip> // For formatting
#include <sstream> // For string streams

struct BinDefinition {
    int bin_number;
    int pT_bin;
    int xi_bin;
    int x_bin;
    int Q2_bin;
    double pT_min;
    double pT_max;
    double xi_min;
    double xi_max;
    double x_min;
    double x_max;
    double Q2_min;
    double Q2_max;
};

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

    // Prepare to store bin definitions
    std::vector<BinDefinition> bin_definitions;
    int bin_counter = 1; // To assign a unique bin number

    // Loop over pT bins
    for (size_t i = 0; i < pT_bin_edges.size() - 1; ++i) {
        double pT_min = pT_bin_edges[i];
        double pT_max = pT_bin_edges[i + 1];
        int pT_bin = i + 1;

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

        // Determine xi bins (now 2 bins instead of 3)
        std::vector<double> xi_bin_edges = find_bin_edges(xi_in_pT_bin, 2, 0.2);

        // Loop over xi bins
        for (size_t k = 0; k < xi_bin_edges.size() - 1; ++k) {
            double xi_min = xi_bin_edges[k];
            double xi_max = xi_bin_edges[k + 1];
            int xi_bin = k + 1;

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
                int x_bin = m + 1;

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
                    int Q2_bin = o + 1;

                    // Print out the bin bounds in the required format
                    std::cout << pT_min << " <= *pT && " << pT_max << " > *pT && ";
                    std::cout << xi_min << " <= *xi && " << xi_max << " > *xi && ";
                    std::cout << x_min << " <= *x && " << x_max << " > *x && ";
                    std::cout << "*Q2 > " << Q2_min << " && *Q2 < " << Q2_max << std::endl;

                    // Store bin definition
                    BinDefinition bin_def;
                    bin_def.bin_number = bin_counter++;
                    bin_def.pT_bin = pT_bin;
                    bin_def.xi_bin = xi_bin;
                    bin_def.x_bin = x_bin;
                    bin_def.Q2_bin = Q2_bin;
                    bin_def.pT_min = pT_min;
                    bin_def.pT_max = pT_max;
                    bin_def.xi_min = xi_min;
                    bin_def.xi_max = xi_max;
                    bin_def.x_min = x_min;
                    bin_def.x_max = x_max;
                    bin_def.Q2_min = Q2_min;
                    bin_def.Q2_max = Q2_max;
                    bin_definitions.push_back(bin_def);
                }
            }
        }
    }

    file->Close();

    // Print LaTeX table
    std::cout << "\n% LaTeX table of bin definitions\n";
    std::cout << "\\begin{table}[h!]\n";
    std::cout << "\\centering\n";
    std::cout << "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}\n";
    std::cout << "\\hline\n";
    std::cout << "Bin & $p_T$ Bin & $\\xi$ Bin & $x$ Bin & $Q^2$ Bin & $p_T$ Range & $\\xi$ Range & $x$ Range & $Q^2$ Range \\\\ \\hline\n";

    for (size_t i = 0; i < bin_definitions.size(); ++i) {
        BinDefinition& bin = bin_definitions[i];

        std::stringstream ss;
        ss << std::fixed << std::setprecision(4);
        ss << bin.bin_number << " & ";
        ss << bin.pT_bin << " & " << bin.xi_bin << " & " << bin.x_bin << " & " << bin.Q2_bin << " & ";
        ss << "$[" << bin.pT_min << ", " << bin.pT_max << ")$ & ";
        ss << "$[" << bin.xi_min << ", " << bin.xi_max << ")$ & ";
        ss << "$[" << bin.x_min << ", " << bin.x_max << ")$ & ";
        ss << "$[" << bin.Q2_min << ", " << bin.Q2_max << ")$ \\\\ \\hline";
        std::cout << ss.str() << std::endl;
    }

    std::cout << "\\end{tabular}\n";
    std::cout << "\\caption{Bin definitions}\n";
    std::cout << "\\label{tab:bin_definitions}\n";
    std::cout << "\\end{table}\n";

    return 0;
}