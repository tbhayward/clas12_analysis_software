#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input_root_file" << std::endl;
        return 1;
    }

    // Load the ROOT file and the tree
    TFile *f = TFile::Open(argv[1]);
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file " << argv[1] << std::endl;
        return 1;
    }

    TTree *tree = (TTree*)f->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: Could not find tree PhysicsEvents in file " << argv[1] << std::endl;
        return 1;
    }

    // Set style options
    gStyle->SetOptStat(0); // Remove the stats box

    // We need two distributions of Mx2_2:
    // 1) For 0 < |t2| < 0.1
    // 2) For 0 < |t| - |tmin| < 0.1

    // Create histograms
    TH1F *h1 = new TH1F("h1", "Mx2_2 Distribution", 250, 0, 5);
    TH1F *h2 = new TH1F("h2", "Mx2_2 Distribution", 250, 0, 5);

    // Draw into histograms with cuts
    // For the first distribution: |t2| < 0.1
    tree->Draw("Mx2_2 >> h1", "fabs(t2)<0.1", "goff"); //endfor

    // For the second distribution: 0 < |t| - |tmin| < 0.1
    // This means we require fabs(t)-fabs(tmin) to be between 0 and 0.1
    tree->Draw("Mx2_2 >> h2", "(fabs(t2)-fabs(tmin))>0 && (fabs(t2)-fabs(tmin))<0.1", "goff"); //endfor

    // Set line colors and styles
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);

    // Create a canvas and draw
    TCanvas *c = new TCanvas("c", "tmin_resolution_test", 800, 600);
    h1->GetXaxis()->SetTitle("M_{x}^{2}");
    h1->Draw(); //endfor
    h2->Draw("same"); //endfor

    // Add a legend
    TLegend *leg = new TLegend(0.65,0.75,0.85,0.85);
    leg->AddEntry(h1, "0 < |t2| < 0.1", "l");
    leg->AddEntry(h2, "0 < |t| - |tmin| < 0.1", "l");
    leg->Draw(); //endfor

    // Save output
    c->SaveAs("output/tmin_resolution_test.png"); //endfor

    // Clean up
    f->Close();
    delete f;
    delete c;
    delete h1;
    delete h2;
    delete leg;

    return 0;
}
//endfor