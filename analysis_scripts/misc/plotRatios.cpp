#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <stdexcept>

void plotRatios(const char* file1, const char* file2) {
    // Open the ROOT files
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);

    if (!f1 || !f2 || f1->IsZombie() || f2->IsZombie()) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }

    // Get the trees
    TTree *tree1 = (TTree*)f1->Get("tree");
    TTree *tree2 = (TTree*)f2->Get("tree");

    if (!tree1 || !tree2) {
        std::cerr << "Error: Could not find the tree in one of the files" << std::endl;
        return;
    }

    // Get the list of branches
    TObjArray *branches1 = tree1->GetListOfBranches();
    TObjArray *branches2 = tree2->GetListOfBranches();

    if (branches1->GetEntries() != branches2->GetEntries()) {
        std::cerr << "Error: The trees have different number of branches" << std::endl;
        return;
    }

    // Loop over branches
    for (int i = 0; i < branches1->GetEntries(); ++i) {
        TBranch *branch1 = (TBranch*)branches1->At(i);
        TBranch *branch2 = (TBranch*)branches2->At(i);

        TString branchName1 = branch1->GetName();
        TString branchName2 = branch2->GetName();

        if (branchName1 != branchName2) {
            std::cerr << "Error: Branch names do not match: " << branchName1 << " vs " << branchName2 << std::endl;
            continue;
        }

        std::vector<float> values1;
        std::vector<float> values2;

        float value1, value2;

        tree1->SetBranchAddress(branchName1, &value1);
        tree2->SetBranchAddress(branchName2, &value2);

        Long64_t nentries1 = tree1->GetEntries();
        Long64_t nentries2 = tree2->GetEntries();

        if (nentries1 != nentries2) {
            std::cerr << "Error: The trees have different number of entries" << std::endl;
            return;
        }

        for (Long64_t j = 0; j < nentries1; ++j) {
            tree1->GetEntry(j);
            tree2->GetEntry(j);
            if (value2 != 0) { // Avoid division by zero
                values1.push_back(value1 / value2);
            } else {
                values1.push_back(0); // or handle it in some other way
            }
        }

        // Create histogram for ratios
        TH1F *h = new TH1F(branchName1, branchName1, 100, 0, *std::max_element(values1.begin(), values1.end()));
        for (auto ratio : values1) {
            h->Fill(ratio);
        }

        // Create a canvas and draw the histogram
        TCanvas *c = new TCanvas("c", "c", 800, 600);
        h->GetXaxis()->SetTitle(branchName1);
        h->GetYaxis()->SetTitle("ratio (NH3/C)");
        h->Draw();

        // Save the canvas as a PNG file
        TString filename = TString::Format("output/%s.png", branchName1.Data());
        c->SaveAs(filename);

        // Clean up
        delete c;
        delete h;
    }

    f1->Close();
    f2->Close();
}

int main(int argc, char **argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " file1.root file2.root" << std::endl;
        return 1;
    }

    plotRatios(argv[1], argv[2]);
    return 0;
}