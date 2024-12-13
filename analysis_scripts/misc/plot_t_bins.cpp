#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <utility>

int main(int argc, char** argv) {
    // Check command line arguments
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.root" << std::endl;
        return 1; //endif
    } //endif

    const char* infile = argv[1];

    // Open the input ROOT file
    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file " << infile << std::endl;
        return 1; //endif
    } //endif

    // Retrieve the PhysicsEvents TTree
    TTree *tree = (TTree*)f->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: PhysicsEvents tree not found in file." << std::endl;
        f->Close();
        return 1; //endif
    } //endif

    // Define the -t bins
    std::vector<std::pair<double, double>> tBins = {
        {0.0,   0.5},
        {0.5,   1.0},
        {1.0,   2.0},
        {2.0,   3.0},
        {3.0,   6.0},
        {6.0,  20.0}
    };

    // Turn off stats
    gStyle->SetOptStat(0);

    // Create a 3x2 canvas
    TCanvas *c = new TCanvas("c", "Mx2 in t-bins", 1200, 800);
    c->Divide(3,2);

    std::vector<TH1F*> hists;
    hists.reserve(tBins.size());

    // Loop over t-bins and plot Mx2
    for (size_t i = 0; i < tBins.size(); i++) {
        c->cd(i+1);
        double tmin = tBins[i].first;
        double tmax = tBins[i].second;

        // Create histogram for this bin
        TString hname = Form("hMx2_bin%zu", i);
        TH1F *h = new TH1F(hname, "", 100, -1, 6);
        h->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h->GetYaxis()->SetTitle("Normalized Counts");

        // Apply cut on -t
        TString cut = Form("(-t>%f && -t<%f)", tmin, tmax);

        // Draw into histogram
        tree->Draw("Mx2 >> "+hname, cut, "E");

        // Normalize each histogram
        double integral = h->Integral();
        if (integral > 0) {
            h->Scale(1.0/integral);
        } //endif

        // Add label with t-bin range
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.15, 0.85, Form("%.1f < -t < %.1f", tmin, tmax));

        hists.push_back(h);
    } //endfor

    // Create output directory and save the plot
    gSystem->Exec("mkdir -p output");
    c->SaveAs("output/rho_contribution.png");

    f->Close();
    return 0;
} //end of main