#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
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
        // ... [same as before]
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
        // ... [same as before]
        "DepW"
    };

    // Histogram configurations
    std::map<std::string, HistConfig> histConfigs = {
        // ... [same as before]
        {"DepW", {50, 0, 1}}
    };

    // Branch types
    std::map<std::string, std::string> branchTypes = {
        // ... [same as before]
        {"DepW", "D"}
    };

    // Get total number of entries
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();

    // Variables to hold Mx2 values
    double Mx2_value1 = 0;
    double Mx2_value2 = 0;

    // Set branch addresses for Mx2
    tree1->SetBranchAddress("Mx2", &Mx2_value1);
    tree2->SetBranchAddress("Mx2", &Mx2_value2);

    // Variables to count entries with Mx2 > 0
    Long64_t nEntriesMx2Pos1 = 0;
    Long64_t nEntriesMx2Pos2 = 0;

    // First, count the number of entries with Mx2 > 0 in each tree
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        if (Mx2_value1 > 0) {
            ++nEntriesMx2Pos1;
        }
    }
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        if (Mx2_value2 > 0) {
            ++nEntriesMx2Pos2;
        }
    }

    // Now loop over branches
    for (const auto& branchName : branches) {
        // Check if branch exists in both trees
        if (!tree1->GetBranch(branchName.c_str())) {
            std::cerr << "Branch " << branchName << " does not exist in tree1" << std::endl;
            continue;
        }
        if (!tree2->GetBranch(branchName.c_str())) {
            std::cerr << "Branch " << branchName << " does not exist in tree2" << std::endl;
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
            tree1->SetBranchAddress(branchName.c_str(), &value1);
            tree2->SetBranchAddress(branchName.c_str(), &value2);

            // Create histograms
            TH1D *hist1 = new TH1D(("hist1_" + branchName).c_str(), branchName.c_str(), histConfig.bins, histConfig.x_min, histConfig.x_max);
            TH1D *hist2 = new TH1D(("hist2_" + branchName).c_str(), branchName.c_str(), histConfig.bins, histConfig.x_min, histConfig.x_max);

            // Fill histograms for tree1
            for (Long64_t i = 0; i < nEntries1; ++i) {
                tree1->GetEntry(i);
                if (Mx2_value1 > 0) {
                    hist1->Fill(value1);
                }
            }
            // Fill histograms for tree2
            for (Long64_t i = 0; i < nEntries2; ++i) {
                tree2->GetEntry(i);
                if (Mx2_value2 > 0) {
                    hist2->Fill(value2);
                }
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
            TCanvas *c = new TCanvas(("c_" + branchName).c_str(), branchName.c_str(), 800, 600);
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

            // Add legend with updated entry counts
            TLegend *legend = new TLegend(0.55, 0.75, 0.9, 0.9); // Adjusted position
            legend->SetTextSize(0.045); // Increased text size
            legend->AddEntry((TObject*)0, ("Dilks Entries (Mx2 > 0): " + std::to_string(nEntriesMx2Pos1)).c_str(), "");
            legend->AddEntry((TObject*)0, ("Hayward Entries (Mx2 > 0): " + std::to_string(nEntriesMx2Pos2)).c_str(), "");
            legend->Draw();

            // Save plot
            std::string outputDir = "output/comparison";
            gSystem->mkdir(outputDir.c_str(), kTRUE);
            std::string outputFileName = outputDir + "/" + branchName + ".png";
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
            tree1->SetBranchAddress(branchName.c_str(), &value1);
            tree2->SetBranchAddress(branchName.c_str(), &value2);

            // Create histograms
            TH1D *hist1 = new TH1D(("hist1_" + branchName).c_str(), branchName.c_str(), histConfig.bins, histConfig.x_min, histConfig.x_max);
            TH1D *hist2 = new TH1D(("hist2_" + branchName).c_str(), branchName.c_str(), histConfig.bins, histConfig.x_min, histConfig.x_max);

            // Adjust phi variables from [-pi, pi] to [0, 2*pi] for tree1
            bool isPhiVariable = (branchName == "e_phi" || branchName == "p_phi" || branchName == "phi");

            // Fill histograms for tree1
            for (Long64_t i = 0; i < nEntries1; ++i) {
                tree1->GetEntry(i);
                if (Mx2_value1 > 0) {
                    if (isPhiVariable) {
                        if (value1 < 0)
                            value1 += 2 * pi;
                    }
                    hist1->Fill(value1);
                }
            }
            // Fill histograms for tree2
            for (Long64_t i = 0; i < nEntries2; ++i) {
                tree2->GetEntry(i);
                if (Mx2_value2 > 0) {
                    hist2->Fill(value2);
                }
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
            TCanvas *c = new TCanvas(("c_" + branchName).c_str(), branchName.c_str(), 800, 600);
            c->SetGrid();
            c->SetLeftMargin(0.15);   // Increase left margin
            c->SetBottomMargin(0.15); // Increase bottom margin

            graph->SetTitle("");
            graph->GetYaxis()->SetTitle("Dilks / Hayward");
            graph->GetYaxis()->SetRangeUser(0.0, 2);
            graph->GetXaxis()->SetTitle(formatLabelName(branchName).c_str());
            graph->SetMarkerStyle(20);
            graph->SetMarkerSize(0.8);
            graph->Draw("AP");

            // Add legend with updated entry counts
            TLegend *legend = new TLegend(0.3, 0.75, 0.9, 0.9); // Adjusted position
            legend->SetTextSize(0.04); // Increased text size
            legend->AddEntry((TObject*)0, ("Dilks Entries: " + std::to_string(nEntriesMx2Pos1)).c_str(), "");
            legend->AddEntry((TObject*)0, ("Hayward Entries: " + std::to_string(nEntriesMx2Pos2)).c_str(), "");
            legend->Draw();

            // Save plot
            std::string outputDir = "output/comparison";
            gSystem->mkdir(outputDir.c_str(), kTRUE);
            std::string outputFileName = outputDir + "/" + branchName + ".png";
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

    // Close files
    file1->Close();
    file2->Close();

    return 0;
}