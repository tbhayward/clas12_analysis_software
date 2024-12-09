#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input_10p5GeV.root input_22GeV.root" << std::endl;
        return 1;
    }

    // Load input files
    TFile *f1 = TFile::Open(argv[1]);
    TFile *f2 = TFile::Open(argv[2]);

    if(!f1 || f1->IsZombie()) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 1;
    }
    if(!f2 || f2->IsZombie()) {
        std::cerr << "Error opening file: " << argv[2] << std::endl;
        return 1;
    }

    // Get trees (assuming they are named "tree")
    TTree *t1 = (TTree*)f1->Get("tree");
    TTree *t2 = (TTree*)f2->Get("tree");

    if(!t1) {
        std::cerr << "Tree named 'tree' not found in " << argv[1] << std::endl;
        return 1;
    }
    if(!t2) {
        std::cerr << "Tree named 'tree' not found in " << argv[2] << std::endl;
        return 1;
    }

    // Create histograms for different pT ranges and files
    // We'll use 50 bins from -1 to 1 for xF
    TH1F *h10_0pT05 = new TH1F("h10_0pT05","",50,-1.0,1.0);
    TH1F *h10_05pT15 = new TH1F("h10_05pT15","",50,-1.0,1.0);
    TH1F *h22_0pT05 = new TH1F("h22_0pT05","",50,-1.0,1.0);
    TH1F *h22_05pT15 = new TH1F("h22_05pT15","",50,-1.0,1.0);

    // Draw into histograms
    t1->Draw("xF>>h10_0pT05","pT>0 && pT<0.5");
    t1->Draw("xF>>h10_05pT15","pT>0.5 && pT<1.5");

    t2->Draw("xF>>h22_0pT05","pT>0 && pT<0.5");
    t2->Draw("xF>>h22_05pT15","pT>0.5 && pT<1.5");

    // Normalize histograms to their integral
    if (h10_0pT05->Integral() > 0) h10_0pT05->Scale(1.0/h10_0pT05->Integral());
    if (h10_05pT15->Integral() > 0) h10_05pT15->Scale(1.0/h10_05pT15->Integral());
    if (h22_0pT05->Integral() > 0) h22_0pT05->Scale(1.0/h22_0pT05->Integral());
    if (h22_05pT15->Integral() > 0) h22_05pT15->Scale(1.0/h22_05pT15->Integral());

    // Set line colors and styles
    // First file (10.5 GeV): red
    h10_0pT05->SetLineColor(kRed);
    h10_05pT15->SetLineColor(kRed);
    // Second file (22 GeV): blue
    h22_0pT05->SetLineColor(kBlue);
    h22_05pT15->SetLineColor(kBlue);

    // Line styles: lower pT bin dashed (2), higher pT bin solid (1)
    h10_0pT05->SetLineStyle(2);
    h10_05pT15->SetLineStyle(1);
    h22_0pT05->SetLineStyle(2);
    h22_05pT15->SetLineStyle(1);

    // Remove stats boxes
    gStyle->SetOptStat(0);
    h10_0pT05->SetStats(0);
    h10_05pT15->SetStats(0);
    h22_0pT05->SetStats(0);
    h22_05pT15->SetStats(0);

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1","c1",800,600);

    // Draw first histogram (set axis labels)
    h10_0pT05->GetXaxis()->SetTitle("x_{F}");
    h10_0pT05->GetYaxis()->SetTitle("normalized counts");
    h10_0pT05->Draw("hist");
    h10_05pT15->Draw("hist same");
    h22_0pT05->Draw("hist same");
    h22_05pT15->Draw("hist same");

    // Add a legend
    TLegend *leg = new TLegend(0.15,0.75,0.5,0.9);
    leg->AddEntry(h10_0pT05,"10.5 GeV, 0 < P_{T} (GeV) < 0.5","l");
    leg->AddEntry(h10_05pT15,"10.5 GeV, 0.5 < P_{T} (GeV) < 1.5","l");
    leg->AddEntry(h22_0pT05,"22 GeV, 0 < P_{T} (GeV) < 0.5","l");
    leg->AddEntry(h22_05pT15,"22 GeV, 0.5 < P_{T} (GeV) < 1.5","l");
    leg->Draw();

    // Create output directory if needed, then save the plot
    // (Assuming the directory "output" already exists or you have permissions)
    c1->SaveAs("output/22gev_comparison.png");

    // Cleanup
    f1->Close();
    f2->Close();

    return 0;
}