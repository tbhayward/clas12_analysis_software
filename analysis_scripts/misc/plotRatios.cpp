#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <algorithm>

void plotRatios(const char* file1, const char* file2, const char* file3, const char* file4) {
    // Open the ROOT files
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TFile *f3 = TFile::Open(file3);
    TFile *f4 = TFile::Open(file4);

    if (!f1 || !f2 || !f3 || !f4 || f1->IsZombie() || f2->IsZombie() || f3->IsZombie() || f4->IsZombie()) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }

    // Get the trees
    TTree *tree1 = (TTree*)f1->Get("tree");
    TTree *tree2 = (TTree*)f2->Get("tree");
    TTree *tree3 = (TTree*)f3->Get("tree");
    TTree *tree4 = (TTree*)f4->Get("tree");

    if (!tree1) {
        std::cerr << "Error: Could not find the tree in one of the files" << std::endl;
        return;
    }

    // Get the list of branches
    TObjArray *branches1 = tree1->GetListOfBranches();
    TObjArray *branches2 = tree2->GetListOfBranches();
    TObjArray *branches3 = tree3->GetListOfBranches();
    TObjArray *branches4 = tree4->GetListOfBranches();

    if (branches1->GetEntries() != branches2->GetEntries() || branches3->GetEntries() != branches4->GetEntries()) {
        std::cerr << "Error: The trees have different number of branches" << std::endl;
        return;
    }

    // Loop over branches
    for (int i = 0; i < branches1->GetEntries(); ++i) {
        TBranch *branch1 = (TBranch*)branches1->At(i);
        TBranch *branch2 = (TBranch*)branches2->At(i);
        TBranch *branch3 = (TBranch*)branches3->At(i);
        TBranch *branch4 = (TBranch*)branches4->At(i);

        TString branchName1 = branch1->GetName();
        TString branchName2 = branch2->GetName();
        TString branchName3 = branch3->GetName();
        TString branchName4 = branch4->GetName();

        if (branchName1 != branchName2 || branchName3 != branchName4) {
            std::cerr << "Error: Branch names do not match: " << branchName1 << " vs " << branchName2 
                      << " or " << branchName3 << " vs " << branchName4 << std::endl;
            continue;
        }

        float value1, value2, value3, value4;

        tree1->SetBranchAddress(branchName1, &value1);
        tree2->SetBranchAddress(branchName2, &value2);
        tree3->SetBranchAddress(branchName3, &value3);
        tree4->SetBranchAddress(branchName4, &value4);

        Long64_t nentries1 = tree1->GetEntries();
        Long64_t nentries2 = tree2->GetEntries();
        Long64_t nentries3 = tree3->GetEntries();
        Long64_t nentries4 = tree4->GetEntries();

        // Create histograms for ratios
        TH1F *h1 = new TH1F(branchName1 + "_ratio1", branchName1 + "_ratio1", 50, 0, 5); // Assuming max ratio range to 5
        TH1F *h2 = new TH1F(branchName3 + "_ratio2", branchName3 + "_ratio2", 50, 0, 5); // Assuming max ratio range to 5

        // Fill histograms for the first set of files
        for (Long64_t j = 0; j < nentries1 && j < nentries2; ++j) {
            tree1->GetEntry(j);
            tree2->GetEntry(j);
            if (value2 != 0) { // Avoid division by zero
                h1->Fill(value1 / value2);
            }
        }

        // Fill histograms for the second set of files
        for (Long64_t j = 0; j < nentries3 && j < nentries4; ++j) {
            tree3->GetEntry(j);
            tree4->GetEntry(j);
            if (value4 != 0) { // Avoid division by zero
                h2->Fill(value3 / value4);
            }
        }

        // Create canvases and draw the histograms
        TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
        h1->GetXaxis()->SetTitle(branchName1);
        h1->GetYaxis()->SetTitle("ratio (NH3/C)");
        h1->Draw();

        TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
        h2->GetXaxis()->SetTitle(branchName3);
        h2->GetYaxis()->SetTitle("ratio (NH3/C)");
        h2->Draw();

        // Save the canvases as PNG files
        TString filename1 = TString::Format("output/%s_ratio1.png", branchName1.Data());
        TString filename2 = TString::Format("output/%s_ratio2.png", branchName3.Data());
        c1->SaveAs(filename1);
        c2->SaveAs(filename2);

        // Clean up
        delete c1;
        delete c2;
        delete h1;
        delete h2;
    }

    f1->Close();
    f2->Close();
    f3->Close();
    f4->Close();
}

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " file1.root file2.root file3.root file4.root" << std::endl;
        return 1;
    }

    plotRatios(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}