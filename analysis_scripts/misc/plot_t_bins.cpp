#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TText.h>
#include <TLatex.h>
#include <iostream>
#include <string>
#include <vector>
#include <TString.h>

void plot_t_bins(const char* infile = "input.root") {

    // Open the file and retrieve the tree
    TFile *f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error opening file " << infile << std::endl;
        return;
    } //endif

    TTree *tree = (TTree*)f->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: PhysicsEvents tree not found in file." << std::endl;
        f->Close();
        return;
    } //endif

    // Bins for -t
    // We'll define the -t bins as requested:
    // [0,0.5], [0.5,1], [1,2], [2,3], [3,6], [6,20]
    std::vector<std::pair<double,double>> tBins = {
        {0.0,   0.5},
        {0.5,   1.0},
        {1.0,   2.0},
        {2.0,   3.0},
        {3.0,   6.0},
        {6.0,  20.0}
    };

    // Turn off stats box globally
    gStyle->SetOptStat(0);

    // Create a canvas with 3x2 pads
    TCanvas *c = new TCanvas("c", "Mx2 in t-bins", 1200, 800);
    c->Divide(3,2);

    // We'll store histograms in a vector to avoid memory issues
    std::vector<TH1F*> hists;
    hists.reserve(tBins.size());

    // Loop over t-bins, draw Mx2 with cuts
    for (size_t i = 0; i < tBins.size(); i++) {
        c->cd(i+1);
        double tmin = tBins[i].first;
        double tmax = tBins[i].second;

        // Create a histogram name/title
        TString hname = Form("hMx2_bin%zu", i);
        TH1F *h = new TH1F(hname, "", 100, -1, 6);
        h->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h->GetYaxis()->SetTitle("Normalized Counts");

        // Cut on -t to be within [tmin, tmax]
        // Condition: t is negative, so -t is positive.
        // The cut string: "(-t>tmin && -t<tmax)"
        TString cut = Form("(-t>%f && -t<%f)", tmin, tmax);

        // Draw Mx2 into the histogram
        tree->Draw("Mx2 >> "+hname, cut, "E");

        // Normalize the histogram to its integral if integral > 0
        double integral = h->Integral();
        if (integral > 0) {
            h->Scale(1.0/integral);
        } //endif

        // Add a label with the t-bin range
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.05);
        latex.DrawLatex(0.15, 0.85, Form("%.1f < -t < %.1f", tmin, tmax));

        hists.push_back(h);
    } //endfor

    // Save the canvas as a PNG
    system("mkdir -p output"); // Make sure the output directory exists
    c->SaveAs("output/rho_contribution.png");

    // Clean up
    f->Close();
} 