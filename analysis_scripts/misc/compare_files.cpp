// Created 9/6/23
// Created to compare epi+X and epX nSidis distributions between pass-1 and preliminary pass-2 

#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <string>
#include <map>
#include <cstring>
#include <algorithm>
#include <cmath>

using namespace std;

struct HistConfig {
    int bins;
    double min;
    double max;
};

std::map<std::string, HistConfig> histConfigs = {
    {"beam_pol", {20, 0.80, 1.00}},
    {"DepA", {200, 0, 1}},
    {"DepB", {200, 0, 1}},
    {"DepC", {200, 0, 1}},
    {"DepV", {200, 0, 2}},
    {"DepW", {200, 0, 1}},
    {"e_p", {200, 2, 8}},
    {"e_phi", {200, 0, 360}},
    {"eta", {200, -1, 3}},
    {"e_theta", {200, 0, 40}}, // Convert degree to radian
    {"evnum", {200, 0, 0}},
    {"helicity", {2, -2, 2}},
    {"Mx", {200, -4, 3}},
    {"Mx2", {200, -10, 10}},
    {"phi", {200, 0, 2 * M_PI}},
    {"p_p", {200, 0, 6}},
    {"p_phi", {200, 0, 360}},
    {"pT", {200, 0, 1.2}},
    {"p_theta", {200, 0, 60}}, // Convert degree to radian
    {"Q2", {200, 0, 9}},
    {"runnum", {200, 0, 0}},
    {"t", {200, -10, 1}},
    {"tmin", {200, -0.5, 0}},
    {"vz_e", {200, -15, 15}},
    {"vz_p", {200, -15, 15}},
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
        {"Mx", "M_{x} (GeV)"},
        {"Mx2", "M_{x}^{2} (GeV)"},
        {"p_p", "p_{p} (GeV)"},
        {"xF", "x_{F}"},
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

void createHistograms(TTree* tree1, TTree* tree2, 
    std::string data_set_1_name, std::string data_set_2_name, const char* outDir) {
    TObjArray* branches1 = tree1->GetListOfBranches();
    TObjArray* branches2 = tree2->GetListOfBranches();

    if (branches1->GetEntries() != branches2->GetEntries()) {
        std::cout << "Number of branches mismatch. Exiting." << std::endl;
        return;
    }

    for (int i = 0; i < branches1->GetEntries(); ++i) {
        const char* branchName = branches1->At(i)->GetName();
        std::string formattedBranchName = formatBranchName(branchName);
        TCanvas canvas(branchName, "Canvas", 1600, 600);  // Width doubled for side-by-side panels
        canvas.SetLeftMargin(0.25);
        TPad *pad1 = new TPad("pad1", "The pad with the function",0.0,0.0,0.5,1.0,21);
        pad1->SetFillColor(0);  // Set the fill color to white for pad1
        pad1->Draw();
        pad1->cd();  // Set current pad to pad1


        HistConfig config = histConfigs[branchName];
        TH1F hist1(Form("%s_1", branchName), "", config.bins, config.min, config.max);
        TH1F hist2(Form("%s_2", branchName), "", config.bins, config.min, config.max);

        std::string cutCondition = "";
        if (std::strcmp(branchName, "Mx") != 0 && std::strcmp(branchName, "Mx2") != 0) {
            cutCondition = "Mx > 1.5";
        }

        std::string drawCommand1 = Form("%s>>%s_1", branchName, branchName);
        std::string drawCommand2 = Form("%s>>%s_2", branchName, branchName);

        // Convert to degrees if necessary
        if (strcmp(branchName, "e_phi") == 0 || strcmp(branchName, "e_theta") == 0 ||
            strcmp(branchName, "p_phi") == 0 || strcmp(branchName, "p_theta") == 0) {
            drawCommand1 = Form("%s * (180 / TMath::Pi())>>%s_1", branchName, branchName);
            drawCommand2 = Form("%s * (180 / TMath::Pi())>>%s_2", branchName, branchName);
        }

        tree1->Draw(drawCommand1.c_str(), cutCondition.c_str());
        tree2->Draw(drawCommand2.c_str(), cutCondition.c_str());

        hist1.SetLineColor(kRed);
        hist2.SetLineColor(kBlue);

        hist1.Draw(""); 
        hist2.Draw("same"); 

        hist1.SetStats(0);
        hist2.SetStats(0);

        hist1.GetXaxis()->SetTitle(formattedBranchName.c_str());
        hist1.GetYaxis()->SetTitle("Counts");

        double max_value = std::max(hist1.GetMaximum(), hist2.GetMaximum());
        hist1.SetMaximum(max_value * 1.1);
        hist2.SetMaximum(max_value * 1.1);

        TPaveText* stats = new TPaveText(0.65, 0.85, 0.85, 0.95, "NDC");
        stats->SetBorderSize(1);  // Draw a border
        stats->SetFillColor(0);  // Transparent fill
        stats->AddText(Form("%s counts: %d", data_set_1_name.c_str(), int(hist1.GetEntries())));
        stats->AddText(Form("%s counts: %d", data_set_2_name.c_str(), int(hist2.GetEntries())));
        stats->SetTextAlign(12);
        stats->Draw("same");

        // pad with ratio
        canvas.cd();  // Switch back to the main canvas before creating a new pad
        TPad *pad2 = new TPad("pad2", "The pad with the ratio",0.5,0.0,1.0,1.0,21);
        pad2->SetFillColor(0);  // Set the fill color to white for pad2
        pad2->Draw();
        pad2->cd();  // Set current pad to pad2
        
        TH1F ratioHist(Form("%s_ratio", branchName), "", config.bins, config.min, config.max);
        ratioHist.Divide(&hist2, &hist1);
        ratioHist.SetLineColor(kBlack);
        ratioHist.SetMinimum(0.75);  // Set Y-range
        ratioHist.SetMaximum(2.25);  // Set Y-range
        ratioHist.GetXaxis()->SetTitle(formattedBranchName.c_str());
        std::string yAxisTitle = Form("%s/%s counts", data_set_2_name.c_str(), 
            data_set_1_name.c_str());
        ratioHist.GetYaxis()->SetTitle(yAxisTitle.c_str());
        ratioHist.Draw("HIST");

        // Ratio stats box
        TPaveText* ratioStats = new TPaveText(0.65, 0.85, 0.85, 0.95, "NDC");
        ratioHist.SetStats(0);
        ratioStats->SetBorderSize(1);  // Draw a border
        ratioStats->SetFillColor(0);  // Transparent fill
        double overallRatio = (hist2.GetEntries() != 0 && hist1.GetEntries() != 0) ? 
            hist2.GetEntries() / hist1.GetEntries() : 0;
        ratioStats->AddText(Form("Overall Ratio: %.2f", overallRatio));
        ratioStats->SetTextAlign(12);
        ratioStats->Draw("same");

        canvas.SaveAs(Form("%s/%s.png", outDir, branchName));
    }
}

void compare_files(std::string root_file1_path, std::string root_file2_path, 
    std::string data_set_1_name, std::string data_set_2_name) {
    gStyle->SetCanvasColor(0);

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

    createHistograms(tree1, tree2, data_set_1_name, data_set_2_name, "output");

    file1->Close();
    file2->Close();
}
