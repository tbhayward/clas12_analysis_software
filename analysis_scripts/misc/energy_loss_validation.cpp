#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <iostream>

void compareTrees(const char* file1, const char* file2, const char* output) {
    // Define the momentum bin edges
    std::vector<double> binEdges = {0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0};
    int nBins = binEdges.size() - 1;

    // Open ROOT files and get trees
    TFile* f1 = new TFile(file1);
    TFile* f2 = new TFile(file2);
    TTree* tree1 = (TTree*)f1->Get("PhysicsEvents"); // Replace "tree_name" with your tree name
    TTree* tree2 = (TTree*)f2->Get("PhysicsEvents"); // Same as above

    // Create histograms for each bin
    std::vector<TH1D*> hist1, hist2;
    for (int i = 0; i < nBins; ++i) {
        hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 0, 4));
        hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 0, 4));
    }

    // Set branch addresses
    double p_p, Mx;
    tree1->SetBranchAddress("p_p", &p_p);
    tree1->SetBranchAddress("Mx2", &Mx);
    tree2->SetBranchAddress("p_p", &p_p);
    tree2->SetBranchAddress("Mx2", &Mx);

    // Fill histograms
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        for (int b = 0; b < nBins; ++b) {
            if (p_p >= binEdges[b] && p_p < binEdges[b+1])
                hist1[b]->Fill(Mx);
        }
    }
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        for (int b = 0; b < nBins; ++b) {
            if (p_p >= binEdges[b] && p_p < binEdges[b+1])
                hist2[b]->Fill(Mx);
        }
    }

    // Create a canvas and draw histograms
    TCanvas* c1 = new TCanvas("c1", "Comparison", 1200, 800);
    c1->Divide(TMath::CeilNint(sqrt(nBins)), TMath::CeilNint(sqrt(nBins)));

    for (int i = 0; i < nBins; ++i) {
    c1->cd(i+1);

    // Find the maximum value in both histograms for this bin
    double maxVal = TMath::Max(hist1[i]->GetMaximum(), hist2[i]->GetMaximum());

    // Set the range of y-axis to 0 - 10% more than the max value
    hist1[i]->SetMaximum(maxVal * 1.10); 

    // Set x and y axis labels
    hist1[i]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    hist1[i]->GetYaxis()->SetTitle("Counts");

    // Draw the histograms
    hist1[i]->SetLineColor(kBlack);
    hist2[i]->SetLineColor(kRed);
    hist1[i]->Draw();
    hist2[i]->Draw("same");

    // Create and add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Adjust these coordinates as needed
    legend->AddEntry(hist1[i], "Uncorrected", "l");
    legend->AddEntry(hist2[i], "Corrected", "l");
    legend->Draw();

    // Remove the statistics box
    gStyle->SetOptStat(0);

    // Save the canvas
    c1->SaveAs(output);
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <file1> <file2> <output>" << std::endl;
        return 1;
    }
    compareTrees(argv[1], argv[2], argv[3]);
    return 0;
}
