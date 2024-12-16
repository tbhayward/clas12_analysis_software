#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TLegend.h>
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

    // Define the original 6 -t bins
    std::vector<std::pair<double,double>> tBins = {
        {0.0,   0.1},
        {0.1,   0.5},
        {0.5,   1.0},
        {1.0,   2.0},
        {2.0,   3.0},
        {3.0,   20.0}
    };

    // Turn off stats
    gStyle->SetOptStat(0);

    // Create a single canvas
    TCanvas *c = new TCanvas("c", "Mx2 in t-bins", 1200, 800);
    // Add some padding on the right so the axis and legend don't get cut off
    c->SetRightMargin(0.1);

    // We'll store histograms
    std::vector<TH1F*> hists;
    hists.reserve(tBins.size());

    // Create a legend in the top-right corner
    TLegend *leg = new TLegend(0.65, 0.7, 0.9, 0.9);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    // Pick colors for the 6 histograms
    Color_t colors[6] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+1};

    bool firstHist = true; // To know when to draw the axis
    for (size_t i = 0; i < tBins.size(); i++) {
        double tmin = tBins[i].first;
        double tmax = tBins[i].second;

        TString hname = Form("hMx2_bin%zu", i);
        // Twice as many bins as before: previously 100, now 200
        TH1F *h = new TH1F(hname, "", 200, -1, 6);
        // Set axis titles
        h->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        h->GetYaxis()->SetTitle("Normalized Counts");

        // Apply cut on -t
        TString cut = Form("(-(t2-tmin)>%f && -(t2-tmin)<%f)", tmin, tmax);
        tree->Draw("Mx2_2 >> "+hname, cut, "goff"); // Fill histogram in memory

        // Normalize the histogram
        double integral = h->Integral();
        if (integral > 0) {
            h->Scale(1.0/integral);
        } //endif

        // Set line style and color
        h->SetLineColor(colors[i % 6]);
        h->SetLineWidth(2);

        // Draw the histogram
        if (firstHist) {
            h->Draw("hist");
            firstHist = false;
        } else {
            h->Draw("hist same");
        } //endif

        // Add to the legend
        leg->AddEntry(h, Form("%.1f < -(t-t_{min}) < %.1f", tmin, tmax), "l");

        hists.push_back(h);
    } //endfor

    // Draw the legend
    leg->Draw();

    // Create output directory and save the plot
    gSystem->Exec("mkdir -p output");
    c->SaveAs("output/rho_contribution.png");

    f->Close();
    return 0;
} //end of main