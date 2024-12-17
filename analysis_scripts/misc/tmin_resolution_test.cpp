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
    TH1F *h1 = new TH1F("h1", "Mx2_2 Distribution", 250, 0, 5);
    TH1F *h2 = new TH1F("h2", "Mx2_2 Distribution", 250, 0, 5);

    // Draw into histograms with cuts
    tree->Draw("Mx2_2 >> h1", "fabs(t2)<0.1", "goff"); //endfor
    tree->Draw("Mx2_2 >> h2", "(fabs(t2)-fabs(tmin))>0 && (fabs(t2)-fabs(tmin))<0.1", "goff"); //endfor

    // Set line colors and styles
    h1->SetLineColor(kBlue);
    h2->SetLineColor(kRed);

    // Now create histograms for |t2| and (|t2|-|tmin|) from -1 to 10
    TH1F *h3 = new TH1F("h3", "|t2| Distribution", 100, -1, 8);
    TH1F *h4 = new TH1F("h4", "|t2|-|tmin| Distribution", 100, -1, 8);

    // Fill those histograms
    tree->Draw("fabs(t2) >> h3", "", "goff"); //endfor
    tree->Draw("(fabs(t2)-fabs(tmin)) >> h4", "", "goff"); //endfor

    h3->SetLineColor(kBlue);
    h4->SetLineColor(kGreen+2);

    // Create a canvas with 1x2 subplots
    TCanvas *c = new TCanvas("c", "tmin_resolution_test", 1200, 600);
    c->Divide(2,1); //endfor

    // Left pad: the original Mx2_2 distributions
    c->cd(1); //endfor
    h1->GetXaxis()->SetTitle("M_{x}^{2}");
    h1->Draw(); //endfor
    h2->Draw("same"); //endfor
    TLegend *leg1 = new TLegend(0.65,0.75,0.85,0.85);
    leg1->AddEntry(h1, "0 < |t2| < 0.1", "l");
    leg1->AddEntry(h2, "0 < |t2|-|tmin| < 0.1", "l");
    leg1->Draw(); //endfor

    // Right pad: distributions of |t2| and |t2|-|tmin|
    c->cd(2); //endfor
    h3->GetXaxis()->SetTitle("Value");
    h3->SetTitle(""); // no title for cleaner look
    h3->Draw(); //endfor
    h4->Draw("same"); //endfor
    TLegend *leg2 = new TLegend(0.65,0.75,0.85,0.85);
    leg2->AddEntry(h3, "|t2| distribution", "l");
    leg2->AddEntry(h4, "|t2|-|tmin| distribution", "l");
    leg2->Draw(); //endfor

    // Save output
    c->SaveAs("output/tmin_resolution_test.png"); //endfor

    // Clean up
    f->Close();
    delete f;
    delete c;
    delete h1;
    delete h2;
    delete h3;
    delete h4;
    delete leg1;
    delete leg2;

    return 0;
}
//endfor