#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>


void compareTrees(const char* file1, const char* file2, const char* output, double lineValue) {
    // Define the momentum bin edges
    std::vector<double> binEdges = {0,0.6,0.7,0.8,1.0,1.2,1.4,1.8,2.4,3.0};
    int nBins = binEdges.size() - 1;

    // Open ROOT files and get trees
    TFile* f1 = new TFile(file1);
    TFile* f2 = new TFile(file2);
    TTree* tree1 = (TTree*)f1->Get("PhysicsEvents"); // Replace "tree_name" with your tree name
    TTree* tree2 = (TTree*)f2->Get("PhysicsEvents"); // Same as above

    // Create histograms for each bin
    std::vector<TH1D*> hist1, hist2;
    for (int i = 0; i < nBins; ++i) {
        hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 0, 1));
        hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 0, 1));
    }

    // Set branch addresses
    double p_p, Mx2;
    tree1->SetBranchAddress("p_p", &p_p);
    tree1->SetBranchAddress("Mx2", &Mx2);
    tree2->SetBranchAddress("p_p", &p_p);
    tree2->SetBranchAddress("Mx2", &Mx2);

    // Fill histograms
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        for (int b = 0; b < nBins; ++b) {
            if (p_p >= binEdges[b] && p_p < binEdges[b+1])
                hist1[b]->Fill(Mx2);
        }
    }
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        for (int b = 0; b < nBins; ++b) {
            if (p_p >= binEdges[b] && p_p < binEdges[b+1])
                hist2[b]->Fill(Mx2);
        }
    }

    TCanvas* c1 = new TCanvas("c1", "Comparison", 1200, 800);
    c1->Divide(TMath::CeilNint(sqrt(nBins)), TMath::CeilNint(sqrt(nBins)));
    double globalFontSize = 0.04; // You can adjust this value as needed
    for (int i = 0; i < nBins; ++i) {
        c1->cd(i+1);
        gPad->SetLeftMargin(0.15);  

        // Setting the title of the histogram to include the bin range
        char title[100];
        sprintf(title, "%.1f < p (GeV) < %.1f; M_{x}^{2} (GeV^{2}); Counts", 
            binEdges[i], binEdges[i+1]);
        hist1[i]->SetTitle(title);

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

        // Create and draw a vertical line
        TLine *line = new TLine(lineValue, hist1[i]->GetMinimum(), lineValue, maxVal * 1.10);
        line->SetLineStyle(2); // Set the line style to dashed
        line->Draw();

        // Remove the statistics box
        gStyle->SetOptStat(0);

        // Adjust font size and center axis labels
        hist1[i]->GetXaxis()->SetTitleSize(globalFontSize);
        hist1[i]->GetYaxis()->SetTitleSize(globalFontSize);
        hist1[i]->GetXaxis()->SetLabelSize(globalFontSize);
        hist1[i]->GetYaxis()->SetLabelSize(globalFontSize);
        
        hist1[i]->GetXaxis()->CenterTitle();
        hist1[i]->GetYaxis()->CenterTitle();
    }


    // Save the canvas
    c1->SaveAs(output);
}

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cout << "Usage: " << argv[0] << " <file1> <file2> <output> <peak value>" << std::endl;
        return 1;
    }
    compareTrees(argv[1], argv[2], argv[3], atof(argv[4]));
    return 0;
}
