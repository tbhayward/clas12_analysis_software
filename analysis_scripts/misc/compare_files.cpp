#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TKey.h>
#include <TList.h>
#include <TPaveText.h>  
#include <iostream>

using namespace std;

// Created 9/6/23
// Created to compare epi+X and epX nSidis distributions between pass-1 and preliminary pass-2 

#include <TCanvas.h>
#include <TFile.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TTree.h>
#include <iostream>
#include <algorithm> // for std::max
#include <string>

std::string formatBranchName(const std::string& original) {
    std::string formatted = original;
    size_t pos = 0;
    while ((pos = formatted.find('_', pos)) != std::string::npos) {
        formatted.replace(pos, 1, "_{");
        size_t closing = formatted.find('_', pos + 2); // find next underscore
        if (closing == std::string::npos) {
            closing = formatted.length(); // if no more underscores, append at the end
        }
        formatted.insert(closing, "}");
        pos = closing + 1;
    }
    
    // Replace specific names with Greek letters, only if they exist
    size_t theta_pos = formatted.find("theta");
    if (theta_pos != std::string::npos) {
        formatted.replace(theta_pos, 5, "#theta");
    }

    size_t phi_pos = formatted.find("phi");
    if (phi_pos != std::string::npos) {
        formatted.replace(phi_pos, 3, "#phi");
    }

    size_t eta_pos = formatted.find("eta");
    if (eta_pos != std::string::npos) {
        formatted.replace(eta_pos, 3, "#eta");
    }

    // Add more replacements as needed
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

        TH1F hist1(Form("%s_1", branchName), "", 100, 0, 0); // Empty title
        TH1F hist2(Form("%s_2", branchName), "", 100, 0, 0); // Empty title

        tree1->Draw(Form("%s>>%s_1", branchName, branchName));
        tree2->Draw(Form("%s>>%s_2", branchName, branchName));

        // Set line colors
        hist1.SetLineColor(kRed);
        hist2.SetLineColor(kBlue);

        // Draw histograms as error bars
        hist1.Draw("E");
        hist2.Draw("Esame");

        // Hide default stats box
        hist1.SetStats(0);
        hist2.SetStats(0);

        // Set axis labels
        hist1.GetXaxis()->SetTitle(formattedBranchName.c_str());
        hist1.GetYaxis()->SetTitle("Counts");

        // Determine the max value between the two histograms
        double max_value = std::max(hist1.GetMaximum(), hist2.GetMaximum());
        hist1.SetMaximum(max_value * 1.1); // Add some margin
        hist2.SetMaximum(max_value * 1.1); // Add some margin

        // Create and draw custom stats boxes
        TPaveText* stats1 = new TPaveText(0.75, 0.75, 0.95, 0.95, "NDC");
        stats1->AddText("pass-1");
        stats1->AddText(Form("Counts: %d", int(hist1.GetEntries())));
        stats1->SetTextAlign(12); // Align text to the left
        stats1->SetTextColor(kRed);
        stats1->SetBorderSize(1); // Add border
        stats1->Draw("same");

        TPaveText* stats2 = new TPaveText(0.55, 0.75, 0.75, 0.95, "NDC");
        stats2->AddText("pass-2");
        stats2->AddText(Form("Counts: %d", int(hist2.GetEntries())));
        stats2->SetTextAlign(12); // Align text to the left
        stats2->SetTextColor(kBlue);
        stats2->SetBorderSize(1); // Add border
        stats2->Draw("same");

        // Save the canvas to a file
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
