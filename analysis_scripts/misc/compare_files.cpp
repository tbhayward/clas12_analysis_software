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

using namespace std;

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
        TCanvas canvas(branchName, "Canvas", 800, 600);

        TH1F hist1(Form("%s_1", branchName), "", 200, 0, 0);
        TH1F hist2(Form("%s_2", branchName), "", 200, 0, 0);

        std::string cutCondition = "";
        if (std::strcmp(branchName, "Mx") != 0 && std::strcmp(branchName, "Mx2") != 0) {
            cutCondition = "Mx > 1.5";
        }

        tree1->Draw(Form("%s>>%s_1", branchName, branchName), cutCondition.c_str());
        tree2->Draw(Form("%s>>%s_2", branchName, branchName), cutCondition.c_str());

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

        TPaveText* stats1 = new TPaveText(0.15, 0.85, 0.35, 0.95, "NDC");
        TPaveText* stats2 = new TPaveText(0.65, 0.85, 0.85, 0.95, "NDC");

        stats1->AddText("pass-1");
        stats1->AddText(Form("Counts: %d", int(hist1.GetEntries())));
        stats1->SetTextAlign(12);
        stats1->SetTextColor(kRed);
        stats1->Draw("same");

        stats2->AddText("pass-2");
        stats2->AddText(Form("Counts: %d", int(hist2.GetEntries())));
        stats2->SetTextAlign(12);
        stats2->SetTextColor(kBlue);
        stats2->Draw("same");

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
