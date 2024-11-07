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
#include "TGraphErrors.h"
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
        {"Delta_eta", "#Delta#eta"},
        {"Delta_phi", "#Delta#phi"},
        {"Delta_phi12", "#Delta#phi_{12}"},
        {"Delta_phi13", "#Delta#phi_{13}"},
        {"Delta_phi23", "#Delta#phi_{23}"},
        {"eta1", "#eta_{1}"},
        {"eta2", "#eta_{2}"},
        {"eta3", "#eta_{3}"},
        {"eta12", "#eta_{12}"},
        {"eta13", "#eta_{13}"},
        {"eta23", "#eta_{23}"},
        {"pT", "P_{T} (GeV)"},
        {"pT1", "P_{1T} (GeV)"},
        {"pT2", "P_{2T} (GeV)"},
        {"pT3", "P_{3T} (GeV)"},
        {"pT12", "P_{12T} (GeV)"},
        {"pT13", "P_{13T} (GeV)"},
        {"pT23", "P_{23T} (GeV)"},
        {"pTpT", "P_{1T}P_{2T} (GeV^{2})"},
        {"Mh", "M_{h} (GeV)"},
        {"Mh12", "M_{h12} (GeV)"},
        {"Mh13", "M_{h13} (GeV)"},
        {"Mh23", "M_{h23} (GeV)"},
        {"t", "t (GeV^{2})"},
        {"phi1", "#phi_{1}"},
        {"phi2", "#phi_{2}"},
        {"phi3", "#phi_{3}"},
        {"phi12", "#phi_{12}"},
        {"phi13", "#phi_{13}"},
        {"phiR", "#phi_{R}"},
        {"tmin", "t_{min} (GeV^{2})"},
        {"e_p", "e_{p} (GeV)"},
        {"Mx", "M_{x} (GeV)"},
        {"Mx1", "M_{x1} (GeV)"},
        {"Mx2", "M_{x2} (GeV)"},
        {"Mx3", "M_{x3} (GeV)"},
        {"Mx12", "M_{x12} (GeV)"},
        {"Mx13", "M_{x13} (GeV)"},
        {"Mx23", "M_{x23} (GeV)"},
        {"p_p", "p_{p} (GeV)"},
        {"p1_p", "p1_{p} (GeV)"},
        {"p2_p", "p2_{p} (GeV)"},
        {"t1", "t_{1}"},
        {"t2", "t_{2}"},
        {"t3", "t_{3}"},
        {"xF", "x_{F}"},
        {"xF1", "x_{F1}"},
        {"xF2", "x_{F2}"},
        {"x", "x_{B}"},
        {"z1", "z_{1}"},
        {"z2", "z_{2}"},
        {"z3", "z_{3}"},
        {"z12", "z_{12}"},
        {"z13", "z_{13}"},
        {"z23", "z_{23}"},
        {"zeta", "#zeta"},
        {"zeta1", "#zeta_{1}"},
        {"zeta2", "#zeta_{2}"},
        {"zeta3", "#zeta_{3}"},
        {"zeta12", "#zeta_{12}"},
        {"zeta13", "#zeta_{13}"},
        {"zeta23", "#zeta_{23}"},
        {"xi", "#xi"},
        {"xi1", "#xi_{1}"},
        {"xi2", "#xi_{2}"},
        {"xi3", "#xi_{3}"},
        {"xi12", "#xi_{12}"},
        {"xi13", "#xi_{13}"},
        {"xi23", "#xi_{23}"},
        {"vz_e", "v_{z_{e}} (cm)"},
        {"vz_p", "v_{z_{p}} (cm)"},
        {"vz_p1", "v_{z_{p1}} (cm)"},
        {"vz_p2", "v_{z_{p2}} (cm)"},
        {"vz_p3", "v_{z_{p3}} (cm)"},
        {"Emiss2", "E_{miss} (GeV)"},
        {"theta_gamma_gamma", "#theta_{#gamma#gamma}"},
        {"pTmiss", "p_{t miss} (GeV)"},
        {"Mxgammasquared", "M_{e'#gammaX}^{2} (GeV^{2})"}
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

    // Set output directories
    std::string comparisonDir = "output/comparison";
    std::string uniqueDir = "output/unique_events";
    gSystem->mkdir(comparisonDir.c_str(), kTRUE);
    gSystem->mkdir(uniqueDir.c_str(), kTRUE);

    // Part 1: Compare branches that exist in both trees
    for (const auto& branchName : branches) {
        bool exists_in_tree1 = tree1->GetBranch(branchName.c_str()) != nullptr;
        bool exists_in_tree2 = tree2->GetBranch(branchName.c_str()) != nullptr;

        if (!exists_in_tree1 || !exists_in_tree2) {
            std::cerr << "Branch " << branchName << " does not exist in both trees. Skipping comparison." << std::endl;
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

        if (branchType == "I") {
            // Integer type
            int value1 = 0;
            int value2 = 0;

            // Set branch addresses for value1 and value2
            tree1->SetBranchAddress(branchName.c_str(), &value1);
            tree2->SetBranchAddress(branchName.c_str(), &value2);

            // Create histograms
            TH1D *hist1 = new TH1D(("hist1_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);
            TH1D *hist2 = new TH1D(("hist2_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

            // Fill histograms for tree1
            for (Long64_t i = 0; i < nEntries1; ++i) {
                tree1->GetEntry(i);
                hist1->Fill(value1);
            }
            // Fill histograms for tree2
            for (Long64_t i = 0; i < nEntries2; ++i) {
                tree2->GetEntry(i);
                hist2->Fill(value2);
            }

            // Normalize histograms
            if (hist1->Integral() != 0)
                hist1->Scale(1.0 / hist1->Integral());
            if (hist2->Integral() != 0)
                hist2->Scale(1.0 / hist2->Integral());

            // Compute ratio and uncertainties
            std::vector<double> xValues;
            std::vector<double> yValues;
            std::vector<double> xErrors;
            std::vector<double> yErrors;

            for (int bin = 1; bin <= hist1->GetNbinsX(); ++bin) {
                double binCenter = hist1->GetBinCenter(bin);
                double content1 = hist1->GetBinContent(bin);
                double content2 = hist2->GetBinContent(bin);
                double error1 = hist1->GetBinError(bin);
                double error2 = hist2->GetBinError(bin);

                // Skip bins where either content is zero to avoid division by zero
                if (content1 > 0 && content2 > 0) {
                    double ratio = content1 / content2;
                    double ratioError = ratio * sqrt( (error1/content1)*(error1/content1) + (error2/content2)*(error2/content2) );
                    xValues.push_back(binCenter);
                    yValues.push_back(ratio);
                    xErrors.push_back(0);
                    yErrors.push_back(ratioError);
                }
            }

            // Create TGraphErrors
            TGraphErrors *graph = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], &xErrors[0], &yErrors[0]);

            // Create canvas and draw
            TCanvas *c = new TCanvas(("c_" + branchName).c_str(), "", 800, 600);
            c->SetGrid();
            c->SetLeftMargin(0.15);   // Increase left margin
            c->SetBottomMargin(0.15); // Increase bottom margin

            graph->SetTitle("");
            graph->GetYaxis()->SetTitle("Dilks / Hayward");
            graph->GetYaxis()->SetRangeUser(0, 2);
            graph->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(0.8);
            graph->Draw("AP");

            // Add dashed line at y=1
            TLine *line = new TLine(histConfig.x_min, 1, histConfig.x_max, 1);
            line->SetLineColor(kGray+2);
            line->SetLineStyle(2);
            line->Draw();

            // Add legend with total entry counts
            TLegend *legend = new TLegend(0.55, 0.75, 0.9, 0.9); // Adjusted position
            legend->SetTextSize(0.045); // Increased text size
            legend->AddEntry((TObject*)0, ("Dilks Entries: " + std::to_string(nEntries1)).c_str(), "");
            legend->AddEntry((TObject*)0, ("Hayward Entries: " + std::to_string(nEntries2)).c_str(), "");
            legend->Draw();

            // Save plot
            std::string outputFileName = comparisonDir + "/" + branchName + ".png";
            c->SaveAs(outputFileName.c_str());

            // Clean up
            delete c;
            delete hist1;
            delete hist2;
            delete graph;
            delete line;

        } else if (branchType == "D") {
            // Double type
            double value1 = 0;
            double value2 = 0;

            // Set branch addresses for value1 and value2
            tree1->SetBranchAddress(branchName.c_str(), &value1);
            tree2->SetBranchAddress(branchName.c_str(), &value2);

            // Create histograms
            TH1D *hist1 = new TH1D(("hist1_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);
            TH1D *hist2 = new TH1D(("hist2_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

            // Adjust phi variables from [-pi, pi] to [0, 2*pi] for tree1
            bool isPhiVariable = (branchName == "e_phi" || branchName == "p_phi" || branchName == "phi");

            // Fill histograms for tree1
            for (Long64_t i = 0; i < nEntries1; ++i) {
                tree1->GetEntry(i);
                if (isPhiVariable) {
                    if (value1 < 0)
                        value1 += 2 * pi;
                }
                hist1->Fill(value1);
            }
            // Fill histograms for tree2
            for (Long64_t i = 0; i < nEntries2; ++i) {
                tree2->GetEntry(i);
                hist2->Fill(value2);
            }

            // Normalize histograms
            if (hist1->Integral() != 0)
                hist1->Scale(1.0 / hist1->Integral());
            if (hist2->Integral() != 0)
                hist2->Scale(1.0 / hist2->Integral());

            // Compute ratio and uncertainties
            std::vector<double> xValues;
            std::vector<double> yValues;
            std::vector<double> xErrors;
            std::vector<double> yErrors;

            for (int bin = 1; bin <= hist1->GetNbinsX(); ++bin) {
                double binCenter = hist1->GetBinCenter(bin);
                double content1 = hist1->GetBinContent(bin);
                double content2 = hist2->GetBinContent(bin);
                double error1 = hist1->GetBinError(bin);
                double error2 = hist2->GetBinError(bin);

                // Skip bins where either content is zero to avoid division by zero
                if (content1 > 0 && content2 > 0) {
                    double ratio = content1 / content2;
                    double ratioError = ratio * sqrt( (error1/content1)*(error1/content1) + (error2/content2)*(error2/content2) );
                    xValues.push_back(binCenter);
                    yValues.push_back(ratio);
                    xErrors.push_back(0);
                    yErrors.push_back(ratioError);
                }
            }

            // Create TGraphErrors
            TGraphErrors *graph = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], &xErrors[0], &yErrors[0]);

            // Create canvas and draw
            TCanvas *c = new TCanvas(("c_" + branchName).c_str(), "", 800, 600);
            c->SetGrid();
            c->SetLeftMargin(0.15);   // Increase left margin
            c->SetBottomMargin(0.15); // Increase bottom margin

            graph->SetTitle("");
            graph->GetYaxis()->SetTitle("Dilks / Hayward");
            graph->GetYaxis()->SetRangeUser(0.4, 1.6);
            graph->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(0.8);
            graph->Draw("AP");

            // Add legend with total entry counts
            TLegend *legend = new TLegend(0.3, 0.75, 0.9, 0.9); // Adjusted position
            legend->SetTextSize(0.04); // Increased text size
            legend->AddEntry((TObject*)0, ("Dilks Entries: " + std::to_string(nEntries1)).c_str(), "");
            legend->AddEntry((TObject*)0, ("Hayward Entries: " + std::to_string(nEntries2)).c_str(), "");
            legend->Draw();

            // Save plot
            std::string outputFileName = comparisonDir + "/" + branchName + ".png";
            c->SaveAs(outputFileName.c_str());

            // Clean up
            delete c;
            delete hist1;
            delete hist2;
            delete graph;

        } else {
            std::cerr << "Unknown branch type for " << branchName << std::endl;
            continue;
        }
    }

    // Part 2: Find entries that are unique to each tree and plot branch distributions
    std::cout << "\nFinding entries unique to each tree and plotting branch distributions...\n" << std::endl;

    // Read 'evnum' from both trees
    std::set<int> evnum_set1;
    std::set<int> evnum_set2;

    int evnum1 = 0;
    int evnum2 = 0;

    tree1->SetBranchAddress("evnum", &evnum1);
    tree2->SetBranchAddress("evnum", &evnum2);

    // Populate sets of event numbers
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        evnum_set1.insert(evnum1);
    }
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        evnum_set2.insert(evnum2);
    }

    // Find event numbers unique to each tree
    std::vector<int> unique_to_tree1;
    std::vector<int> unique_to_tree2;

    std::set_difference(evnum_set1.begin(), evnum_set1.end(), evnum_set2.begin(), evnum_set2.end(),
                        std::back_inserter(unique_to_tree1));
    std::set_difference(evnum_set2.begin(), evnum_set2.end(), evnum_set1.begin(), evnum_set1.end(),
                        std::back_inserter(unique_to_tree2));

    // Part 2a: Plot branch distributions for events unique to Dilks tree
    if (!unique_to_tree1.empty()) {
        std::cout << "Plotting branch distributions for events unique to Dilks tree..." << std::endl;

        // Create maps to hold branch values
        std::map<std::string, std::vector<double>> branchValues;

        // Initialize branch values
        for (const auto& branchName : branches) {
            branchValues[branchName] = std::vector<double>();
        }

        // Set branch addresses
        std::map<std::string, double> doubleValues;
        std::map<std::string, int> intValues;

        for (const auto& branchName : branches) {
            std::string branchType = branchTypes[branchName];
            if (branchType == "I") {
                intValues[branchName] = 0;
                tree1->SetBranchAddress(branchName.c_str(), &intValues[branchName]);
            } else if (branchType == "D") {
                doubleValues[branchName] = 0;
                tree1->SetBranchAddress(branchName.c_str(), &doubleValues[branchName]);
            }
        }

        // Create a set for faster lookup
        std::set<int> unique_evnums_set(unique_to_tree1.begin(), unique_to_tree1.end());

        // Loop over tree1 and collect branch values for unique events
        for (Long64_t i = 0; i < nEntries1; ++i) {
            tree1->GetEntry(i);
            if (unique_evnums_set.count(evnum1)) {
                for (const auto& branchName : branches) {
                    std::string branchType = branchTypes[branchName];
                    if (branchType == "I") {
                        branchValues[branchName].push_back(intValues[branchName]);
                    } else if (branchType == "D") {
                        branchValues[branchName].push_back(doubleValues[branchName]);
                    }
                }
            }
        }

        // Now, for each branch, create a histogram and plot
        for (const auto& branchName : branches) {
            HistConfig histConfig = histConfigs[branchName];
            std::string branchType = branchTypes[branchName];

            TH1D *hist = new TH1D(("hist_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

            // Fill histogram
            for (const auto& value : branchValues[branchName]) {
                double val = value;
                // Adjust phi variables from [-pi, pi] to [0, 2*pi]
                bool isPhiVariable = (branchName == "e_phi" || branchName == "p_phi" || branchName == "phi");
                if (isPhiVariable && branchType == "D") {
                    if (val < 0)
                        val += 2 * pi;
                }
                hist->Fill(val);
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
            std::string outputFileName = uniqueDir + "/" + branchName + "_Dilks.png";
            c->SaveAs(outputFileName.c_str());

            // Clean up
            delete c;
            delete hist;
        }
    } else {
        std::cout << "No unique events found in Dilks tree." << std::endl;
    }

    // Part 2b: Plot branch distributions for events unique to Hayward tree
    if (!unique_to_tree2.empty()) {
        std::cout << "Plotting branch distributions for events unique to Hayward tree..." << std::endl;

        // Create maps to hold branch values
        std::map<std::string, std::vector<double>> branchValues;

        // Initialize branch values
        for (const auto& branchName : branches) {
            branchValues[branchName] = std::vector<double>();
        }

        // Set branch addresses
        std::map<std::string, double> doubleValues;
        std::map<std::string, int> intValues;

        for (const auto& branchName : branches) {
            std::string branchType = branchTypes[branchName];
            if (branchType == "I") {
                intValues[branchName] = 0;
                tree2->SetBranchAddress(branchName.c_str(), &intValues[branchName]);
            } else if (branchType == "D") {
                doubleValues[branchName] = 0;
                tree2->SetBranchAddress(branchName.c_str(), &doubleValues[branchName]);
            }
        }

        // Create a set for faster lookup
        std::set<int> unique_evnums_set(unique_to_tree2.begin(), unique_to_tree2.end());

        // Loop over tree2 and collect branch values for unique events
        for (Long64_t i = 0; i < nEntries2; ++i) {
            tree2->GetEntry(i);
            if (unique_evnums_set.count(evnum2)) {
                for (const auto& branchName : branches) {
                    std::string branchType = branchTypes[branchName];
                    if (branchType == "I") {
                        branchValues[branchName].push_back(intValues[branchName]);
                    } else if (branchType == "D") {
                        branchValues[branchName].push_back(doubleValues[branchName]);
                    }
                }
            }
        }

        // Now, for each branch, create a histogram and plot
        for (const auto& branchName : branches) {
            HistConfig histConfig = histConfigs[branchName];
            std::string branchType = branchTypes[branchName];

            TH1D *hist = new TH1D(("hist_" + branchName).c_str(), "", histConfig.bins, histConfig.x_min, histConfig.x_max);

            // Fill histogram
            for (const auto& value : branchValues[branchName]) {
                double val = value;
                // Adjust phi variables from [-pi, pi] to [0, 2*pi]
                bool isPhiVariable = (branchName == "e_phi" || branchName == "p_phi" || branchName == "phi");
                if (isPhiVariable && branchType == "D") {
                    if (val < 0)
                        val += 2 * pi;
                }
                hist->Fill(val);
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
            std::string outputFileName = uniqueDir + "/" + branchName + "_Hayward.png";
            c->SaveAs(outputFileName.c_str());

            // Clean up
            delete c;
            delete hist;
        }
    } else {
        std::cout << "No unique events found in Hayward tree." << std::endl;
    }

    // Close files
    file1->Close();
    file2->Close();

    return 0;
}