#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <set>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLine.h"
#include "TMath.h"

struct HistConfig {
    int bins;
    double x_min;
    double x_max;
};

// Function to format axis labels
std::string formatLabelName(const std::string& original) {
    std::map<std::string, std::string> specialLabels = {
        {"Q2", "Q^{2} (GeV^{2})"},
        {"W", "W (GeV)"},
        // ... [Same as before, include all special labels]
    };
  
    if (specialLabels.find(original) != specialLabels.end()) {
        return specialLabels[original];
    }

    std::string formatted = original;
    size_t pos = 0;
    while ((pos = formatted.find('_', pos)) != std::string::npos) {
        formatted.replace(pos, 1, "_{");
        size_t closing = formatted.find('_', pos + 2);
        if (closing == std::string::npos) {
            closing = formatted.length();
        }
        formatted.insert(closing, "}");
        pos = closing + 1;
    }

    if (formatted.find("theta") != std::string::npos) {
        formatted.replace(formatted.find("theta"), 5, "#theta");
    }

    if (formatted.find("zeta") != std::string::npos) {
        formatted.replace(formatted.find("zeta"), 5, "#zeta");
    }

    if (formatted.find("phi") != std::string::npos) {
        formatted.replace(formatted.find("phi"), 3, "#phi");
    }

    if (formatted.find("eta") != std::string::npos && 
        formatted.find("theta") == std::string::npos && 
        formatted.find("zeta") == std::string::npos) {
        formatted.replace(formatted.find("eta"), 3, "#eta");
    }
  
    return formatted;
}

int main() {
    // Set style to remove stat boxes and increase font sizes
    gStyle->SetOptStat(0);
    gStyle->SetLegendTextSize(0.04);
    gStyle->SetLabelSize(0.05, "XYZ");
    gStyle->SetTitleSize(0.06, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");

    const double pi = TMath::Pi();

    // Open the first ROOT file and get the tree "tree"
    TFile *file1 = TFile::Open("/volatile/clas12/thayward/cross_check_rgc_epX/step1_EB_yields/dilks_files/rgc_su22_inb_NH3_run016346_EB_yields.root");
    if (!file1 || file1->IsZombie()) {
        std::cerr << "Error opening file1!" << std::endl;
        return -1;
    }
    TTree *tree1 = (TTree*)file1->Get("tree");
    if (!tree1) {
        std::cerr << "Error getting tree from file1!" << std::endl;
        return -1;
    }

    // Open the second ROOT file and get the tree "PhysicsEvents"
    TFile *file2 = TFile::Open("/volatile/clas12/thayward/cross_check_rgc_epX/step1_EB_yields/hayward_files/rgc_su22_inb_NH3_run016346_EB_yields.root");
    if (!file2 || file2->IsZombie()) {
        std::cerr << "Error opening file2!" << std::endl;
        return -1;
    }
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents");
    if (!tree2) {
        std::cerr << "Error getting tree from file2!" << std::endl;
        return -1;
    }

    // List of branches to process
    std::vector<std::string> branches = {
        "fiducial_status",
        "runnum",
        "num_pos",
        "num_neg",
        "num_neutral",
        "evnum",
        "helicity",
        "detector",
        "beam_pol",
        "target_pol",
        "e_p",
        "e_theta",
        "e_phi",
        "vz_e",
        "p_p",
        "p_theta",
        "p_phi",
        "vz_p",
        "open_angle",
        "Q2",
        "W",
        "Mx2",
        "x",
        "y",
        "t",
        "tmin",
        "z",
        "xF",
        "pT",
        "xi",
        "eta",
        "phi",
        "DepA",
        "DepB",
        "DepC",
        "DepV",
        "DepW"
    };

    // Histogram configurations
    std::map<std::string, HistConfig> histConfigs = {
        {"fiducial_status", {50, -5, 5}},
        {"runnum", {50, 16135, 16773}},
        {"num_pos", {50, 0, 10}},
        {"num_neg", {50, 0, 10}},
        {"num_neutral", {50, 0, 10}},
        {"evnum", {50, 0, 0}},
        {"helicity", {4, -2, 2}},
        {"detector", {50, 0, 10}},
        {"beam_pol", {50, -1, 1}},
        {"target_pol", {50, -1, 1}},
        {"e_p", {50, 1, 9}},
        {"e_theta", {50, 0, 0.9}},
        {"e_phi", {50, 0, 2 * pi}},
        {"vz_e", {50, -15, 15}},
        {"p_p", {50, 0, 6}},
        {"p_theta", {50, 0, 0.9}},
        {"p_phi", {50, 0, 2 * pi}},
        {"vz_p", {50, -15, 15}},
        {"open_angle", {50, 0, 180}},
        {"Q2", {50, 0, 10}},
        {"W", {50, 2, 4}},
        {"Mx2", {50, -6, 10}},
        {"x", {50, 0, 0.7}},
        {"y", {50, 0.0, 1.00}},
        {"t", {50, -12, 1}},
        {"tmin", {50, -0.5, 0}},
        {"z", {50, 0, 1}},
        {"xF", {50, -3, 1}},
        {"pT", {50, 0, 1.2}},
        {"xi", {50, -1, 2}},
        {"eta", {50, -3, 3}},
        {"phi", {50, 0, 2 * pi}},
        {"DepA", {50, 0, 1}},
        {"DepB", {50, 0, 1}},
        {"DepC", {50, 0, 1}},
        {"DepV", {50, 0, 2}},
        {"DepW", {50, 0, 1}}
    };

    // Branch types
    std::map<std::string, std::string> branchTypes = {
        {"fiducial_status", "I"},
        {"runnum", "I"},
        {"num_pos", "I"},
        {"num_neg", "I"},
        {"num_neutral", "I"},
        {"evnum", "I"},
        {"helicity", "I"},
        {"detector", "I"},
        {"beam_pol", "D"},
        {"target_pol", "D"},
        {"e_p", "D"},
        {"e_theta", "D"},
        {"e_phi", "D"},
        {"vz_e", "D"},
        {"p_p", "D"},
        {"p_theta", "D"},
        {"p_phi", "D"},
        {"vz_p", "D"},
        {"open_angle", "D"},
        {"Q2", "D"},
        {"W", "D"},
        {"Mx2", "D"},
        {"x", "D"},
        {"y", "D"},
        {"t", "D"},
        {"tmin", "D"},
        {"z", "D"},
        {"xF", "D"},
        {"pT", "D"},
        {"xi", "D"},
        {"eta", "D"},
        {"phi", "D"},
        {"DepA", "D"},
        {"DepB", "D"},
        {"DepC", "D"},
        {"DepV", "D"},
        {"DepW", "D"}
    };

    // Get total number of entries
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();

    // Set output directory
    std::string outputDir = "output/unique_branches";
    gSystem->mkdir(outputDir.c_str(), kTRUE);

    for (const auto& branchName : branches) {
        bool exists_in_tree1 = tree1->GetBranch(branchName.c_str()) != nullptr;
        bool exists_in_tree2 = tree2->GetBranch(branchName.c_str()) != nullptr;

        // If the branch exists in both trees, skip it
        if (exists_in_tree1 && exists_in_tree2) {
            continue;
        }

        // Get hist config for this branch
        HistConfig histConfig;
        if (histConfigs.find(branchName) != histConfigs.end()) {
            histConfig = histConfigs[branchName];
        } else {
            // If not specified, use default values
            histConfig = {50, 0, 1};
            std::cerr << "Histogram config for " << branchName << " not found. Using default values." << std::endl;
        }

        // Determine the type of the branch
        std::string branchType = branchTypes[branchName];

        if (exists_in_tree1 && !exists_in_tree2) {
            // Branch exists only in tree1 (Dilks)
            if (branchType == "I") {
                int value = 0;
                tree1->SetBranchAddress(branchName.c_str(), &value);

                TH1D *hist = new TH1D(("hist_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

                for (Long64_t i = 0; i < nEntries1; ++i) {
                    tree1->GetEntry(i);
                    hist->Fill(value);
                }

                // Normalize histogram
                if (hist->Integral() != 0)
                    hist->Scale(1.0 / hist->Integral());

                // Create canvas and draw
                TCanvas *c = new TCanvas(("c_" + branchName).c_str(), "", 800, 600);
                c->SetGrid();
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.15);

                hist->SetTitle("");
                hist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
                hist->GetYaxis()->SetTitle("Normalized Counts");
                hist->Draw("HIST");

                // Save plot
                std::string outputFileName = outputDir + "/" + branchName + "_Dilks.png";
                c->SaveAs(outputFileName.c_str());

                // Clean up
                delete c;
                delete hist;

            } else if (branchType == "D") {
                double value = 0;
                tree1->SetBranchAddress(branchName.c_str(), &value);

                TH1D *hist = new TH1D(("hist_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

                // Adjust phi variables from [-pi, pi] to [0, 2*pi]
                bool isPhiVariable = (branchName == "e_phi" || branchName == "p_phi" || branchName == "phi");

                for (Long64_t i = 0; i < nEntries1; ++i) {
                    tree1->GetEntry(i);
                    if (isPhiVariable) {
                        if (value < 0)
                            value += 2 * pi;
                    }
                    hist->Fill(value);
                }

                // Normalize histogram
                if (hist->Integral() != 0)
                    hist->Scale(1.0 / hist->Integral());

                // Create canvas and draw
                TCanvas *c = new TCanvas(("c_" + branchName).c_str(), "", 800, 600);
                c->SetGrid();
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.15);

                hist->SetTitle("");
                hist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
                hist->GetYaxis()->SetTitle("Normalized Counts");
                hist->Draw("HIST");

                // Save plot
                std::string outputFileName = outputDir + "/" + branchName + "_Dilks.png";
                c->SaveAs(outputFileName.c_str());

                // Clean up
                delete c;
                delete hist;

            } else {
                std::cerr << "Unknown branch type for " << branchName << std::endl;
                continue;
            }
        }

        if (exists_in_tree2 && !exists_in_tree1) {
            // Branch exists only in tree2 (Hayward)
            if (branchType == "I") {
                int value = 0;
                tree2->SetBranchAddress(branchName.c_str(), &value);

                TH1D *hist = new TH1D(("hist_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

                for (Long64_t i = 0; i < nEntries2; ++i) {
                    tree2->GetEntry(i);
                    hist->Fill(value);
                }

                // Normalize histogram
                if (hist->Integral() != 0)
                    hist->Scale(1.0 / hist->Integral());

                // Create canvas and draw
                TCanvas *c = new TCanvas(("c_" + branchName).c_str(), "", 800, 600);
                c->SetGrid();
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.15);

                hist->SetTitle("");
                hist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
                hist->GetYaxis()->SetTitle("Normalized Counts");
                hist->Draw("HIST");

                // Save plot
                std::string outputFileName = outputDir + "/" + branchName + "_Hayward.png";
                c->SaveAs(outputFileName.c_str());

                // Clean up
                delete c;
                delete hist;

            } else if (branchType == "D") {
                double value = 0;
                tree2->SetBranchAddress(branchName.c_str(), &value);

                TH1D *hist = new TH1D(("hist_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

                // Adjust phi variables from [-pi, pi] to [0, 2*pi]
                bool isPhiVariable = (branchName == "e_phi" || branchName == "p_phi" || branchName == "phi");

                for (Long64_t i = 0; i < nEntries2; ++i) {
                    tree2->GetEntry(i);
                    if (isPhiVariable) {
                        if (value < 0)
                            value += 2 * pi;
                    }
                    hist->Fill(value);
                }

                // Normalize histogram
                if (hist->Integral() != 0)
                    hist->Scale(1.0 / hist->Integral());

                // Create canvas and draw
                TCanvas *c = new TCanvas(("c_" + branchName).c_str(), "", 800, 600);
                c->SetGrid();
                c->SetLeftMargin(0.15);
                c->SetBottomMargin(0.15);

                hist->SetTitle("");
                hist->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
                hist->GetYaxis()->SetTitle("Normalized Counts");
                hist->Draw("HIST");

                // Save plot
                std::string outputFileName = outputDir + "/" + branchName + "_Hayward.png";
                c->SaveAs(outputFileName.c_str());

                // Clean up
                delete c;
                delete hist;

            } else {
                std::cerr << "Unknown branch type for " << branchName << std::endl;
                continue;
            }
        }
    }

    // Close files
    file1->Close();
    file2->Close();

    return 0;
}