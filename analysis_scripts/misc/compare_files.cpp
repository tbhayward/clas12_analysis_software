#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <map>
#include <cmath>

using namespace std;

// Created 9/6/23
// Created to compare epi+X and epX nSidis distributions between pass-1 and preliminary pass-2 

struct HistConfig {
    int bins;
    double min;
    double max;
};

std::map<std::string, HistConfig> histConfigs = {
    {"beam_pol", {20, 0.80, 1.00}},
    {"e_p", {200, 2, 9}},
    {"e_phi", {200, 0, 2 * M_PI}},
    {"eta", {200, -1, 3}},
    {"e_theta", {200, 0, 40 * (M_PI / 180.0)}}, // Convert degree to radian
    {"Mx", {200, -4, 3}},
    {"Mx2", {200, -10, 10}},
    {"phi", {200, 0, 2 * M_PI}},
    {"p_p", {200, 0, 6}},
    {"p_phi", {200, 0, 2 * M_PI}},
    {"pT", {200, 0, 1.3}},
    {"p_theta", {200, 0, 90 * (M_PI / 180.0)}}, // Convert degree to radian
    {"Q2", {200, 0, 9}},
    {"t", {200, -10, 1}},
    {"tmin", {200, -0.5, 0}},
    {"vze", {200, -15, 15}},
    {"vzp", {200, -15, 15}},
    {"W", {200, 2, 4}},
    {"x", {200, 0, 1}},
    {"xF", {200, -1, 1}},
    {"y", {200, 0, 1}},
    {"z", {200, 0, 1}},
    {"zeta", {200, 0, 1}}
};

std::string formatBranchName(const std::string& original) {
    std::map<std::string, std::string> specialLabels = {
        {"Q2", "Q^{2} (GeV^{2})"},
        {"W", "W (GeV)"},
        {"pT", "P_{T} (GeV)"},
        {"t", "t (GeV^{2})"},
        {"tmin", "t_{min} (GeV^{2})"},
        {"e_p", "e_{p} (GeV)"},
        {"p_p", "p_{p} (GeV)"}
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

    if (formatted.find("phi") != std::string::npos) {
        formatted.replace(formatted.find("phi"), 3, "#phi");
    }

    if (formatted.find("eta") != std::string::npos && formatted.find("theta") == std::string::npos) {
        formatted.replace(formatted.find("eta"), 3, "#eta");
    }
  
    return formatted;
}

void createHistograms(TTree* tree1, TTree* tree2, const char* outDir) {
    TObjArray* branches1 = tree1->GetListOfBranches();
    TObjArray* branches2 = tree2->GetListOfBranches();

    if (branches1->GetEntries() != branches2->GetEntries()) {
        std::cout << "Number of branches mismatch. Exiting." << std::endl;
        return;
    }

    for (int i = 0; i < branches1->GetEntries(); ++i) {
        const char* branchName = branches1->At(i)->GetName();
        std::string formattedBranchName = formatBranchName(branchName);
        
        HistConfig config = histConfigs[branchName];
        
        // Special handling for theta values to convert from radians to degrees
        double factor = 1.0;
        if (std::string(branchName) == "e_theta" || std::string(branchName) == "p_theta") {
            factor = 180.0 / M_PI;
        }

        TH1F hist1(Form("%s_1", branchName), "", config.bins, config.min * factor, config.max * factor);
        TH1F hist2(Form("%s_2", branchName), "", config.bins, config.min * factor, config.max * factor);

        std::string cutCondition = "";
        if (std::strcmp(branchName, "Mx") != 0 && std::strcmp(branchName, "Mx2") != 0) {
            cutCondition = "Mx > 1.5";
        }

        tree1->Draw(Form("%s * %f >> %s_1", branchName, factor, branchName), cutCondition.c_str());
        tree2->Draw(Form("%s * %f >> %s_2", branchName, factor, branchName), cutCondition.c_str());

        // ... (Your existing code for drawing, setting colors, and stats)

        TPaveText* stats = new TPaveText(0.65, 0.85, 0.85, 0.95, "NDC");
        stats->AddText(Form("pass-1: %d", int(hist1.GetEntries())));
        stats->AddText(Form("pass-2: %d", int(hist2.GetEntries())));
        stats->SetTextAlign(12);
        stats->Draw("same");

        canvas.SaveAs(Form("%s/%s.png", outDir, branchName));
    }
}



void compare_files(std::string root_file1_path, std::string root_file2_path) {

    TFile* file1 = new TFile(root_file1_path.c_str(), "READ");
    TFile* file2 = new TFile(root_file2_path.c_str(), "READ");

    if (!file1->IsOpen() || !file2->IsOpen()) {
        cout << "Error opening ROOT files (is the location correct?). Exiting." << endl;
    }

    TTree* tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree* tree2 = (TTree*)file2->Get("PhysicsEvents");

    if (!tree1 || !tree2) {
        cout << "Error getting trees from ROOT files." << endl;
    }

    createHistograms(tree1, tree2, "output");

    file1->Close();
    file2->Close();
}
