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

    // Get trees
    TTree *t1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *t2 = (TTree*)f2->Get("PhysicsEvents");

    if(!t1) {
        std::cerr << "Tree named 'PhysicsEvents' not found in " << argv[1] << std::endl;
        return 1;
    }
    if(!t2) {
        std::cerr << "Tree named 'PhysicsEvents' not found in " << argv[2] << std::endl;
        return 1;
    }

    // We will now do the same analysis twice:
    // Left pad: use pT2 (proton)
    // Right pad: use pT1 (pion)

    // ---------------------------
    // Proton (pT2) histograms (Left pad)
    // ---------------------------
    TH1F *h10_0pT05_pT2 = new TH1F("h10_0pT05_pT2","",50,-1.0,1.0);
    TH1F *h10_05pT15_pT2 = new TH1F("h10_05pT15_pT2","",50,-1.0,1.0);
    TH1F *h22_0pT05_pT2 = new TH1F("h22_0pT05_pT2","",50,-1.0,1.0);
    TH1F *h22_05pT15_pT2 = new TH1F("h22_05pT15_pT2","",50,-1.0,1.0);

    // Draw into histograms for pT2
    t1->Draw("xF>>h10_0pT05_pT2","pT2>0 && pT2<0.5","");
    t1->Draw("xF>>h10_05pT15_pT2","pT2>0.5 && pT2<1.5","");

    t2->Draw("xF>>h22_0pT05_pT2","pT2>0 && pT2<0.5","");
    t2->Draw("xF>>h22_05pT15_pT2","pT2>0.5 && pT2<1.5","");

    // Normalize histograms
    if (h10_0pT05_pT2->Integral() > 0) h10_0pT05_pT2->Scale(1.0/h10_0pT05_pT2->Integral());
    if (h10_05pT15_pT2->Integral() > 0) h10_05pT15_pT2->Scale(1.0/h10_05pT15_pT2->Integral());
    if (h22_0pT05_pT2->Integral() > 0) h22_0pT05_pT2->Scale(1.0/h22_0pT05_pT2->Integral());
    if (h22_05pT15_pT2->Integral() > 0) h22_05pT15_pT2->Scale(1.0/h22_05pT15_pT2->Integral());

    // Determine max for pT2 set
    double max_val_pT2 = 0;
    double temp_max = h10_0pT05_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;
    temp_max = h10_05pT15_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;
    temp_max = h22_0pT05_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;
    temp_max = h22_05pT15_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;

    // Set colors and styles for pT2
    h10_0pT05_pT2->SetLineColor(kRed);   h10_0pT05_pT2->SetLineStyle(2);
    h10_05pT15_pT2->SetLineColor(kRed);  h10_05pT15_pT2->SetLineStyle(1);
    h22_0pT05_pT2->SetLineColor(kBlue);  h22_0pT05_pT2->SetLineStyle(2);
    h22_05pT15_pT2->SetLineColor(kBlue); h22_05pT15_pT2->SetLineStyle(1);

    h10_0pT05_pT2->SetStats(0);
    h10_05pT15_pT2->SetStats(0);
    h22_0pT05_pT2->SetStats(0);
    h22_05pT15_pT2->SetStats(0);

    // ---------------------------
    // Pion (pT1) histograms (Right pad)
    // ---------------------------
    TH1F *h10_0pT05_pT1 = new TH1F("h10_0pT05_pT1","",50,-1.0,1.0);
    TH1F *h10_05pT15_pT1 = new TH1F("h10_05pT15_pT1","",50,-1.0,1.0);
    TH1F *h22_0pT05_pT1 = new TH1F("h22_0pT05_pT1","",50,-1.0,1.0);
    TH1F *h22_05pT15_pT1 = new TH1F("h22_05pT15_pT1","",50,-1.0,1.0);

    // Draw into histograms for pT1
    t1->Draw("xF>>h10_0pT05_pT1","pT1>0 && pT1<0.5","");
    t1->Draw("xF>>h10_05pT15_pT1","pT1>0.5 && pT1<1.5","");

    t2->Draw("xF>>h22_0pT05_pT1","pT1>0 && pT1<0.5","");
    t2->Draw("xF>>h22_05pT15_pT1","pT1>0.5 && pT1<1.5","");

    // Normalize histograms (pT1)
    if (h10_0pT05_pT1->Integral() > 0) h10_0pT05_pT1->Scale(1.0/h10_0pT05_pT1->Integral());
    if (h10_05pT15_pT1->Integral() > 0) h10_05pT15_pT1->Scale(1.0/h10_05pT15_pT1->Integral());
    if (h22_0pT05_pT1->Integral() > 0) h22_0pT05_pT1->Scale(1.0/h22_0pT05_pT1->Integral());
    if (h22_05pT15_pT1->Integral() > 0) h22_05pT15_pT1->Scale(1.0/h22_05pT15_pT1->Integral());

    // Determine max for pT1 set
    double max_val_pT1 = 0;
    temp_max = h10_0pT05_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;
    temp_max = h10_05pT15_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;
    temp_max = h22_0pT05_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;
    temp_max = h22_05pT15_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;

    // Set colors and styles for pT1
    h10_0pT05_pT1->SetLineColor(kRed);   h10_0pT05_pT1->SetLineStyle(2);
    h10_05pT15_pT1->SetLineColor(kRed);  h10_05pT15_pT1->SetLineStyle(1);
    h22_0pT05_pT1->SetLineColor(kBlue);  h22_0pT05_pT1->SetLineStyle(2);
    h22_05pT15_pT1->SetLineColor(kBlue); h22_05pT15_pT1->SetLineStyle(1);

    h10_0pT05_pT1->SetStats(0);
    h10_05pT15_pT1->SetStats(0);
    h22_0pT05_pT1->SetStats(0);
    h22_05pT15_pT1->SetStats(0);

    gStyle->SetOptStat(0);

    // Create a 1x2 canvas
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1); // 1 row, 2 columns

    // Left pad: pT2 (proton)
    c1->cd(1);
    h10_0pT05_pT2->GetXaxis()->SetTitle("x_{F}");
    h10_0pT05_pT2->GetYaxis()->SetTitle("normalized counts");
    h10_0pT05_pT2->SetMaximum(1.25 * max_val_pT2);
    // Set pad title as "proton"
    // One way is to set the title on the histogram:
    h10_0pT05_pT2->SetTitle("proton");
    h10_0pT05_pT2->Draw("hist");
    h10_05pT15_pT2->Draw("hist same");
    h22_0pT05_pT2->Draw("hist same");
    h22_05pT15_pT2->Draw("hist same");

    // Legend on left pad
    {
        TLegend *leg = new TLegend(0.15,0.75,0.5,0.9);
        leg->AddEntry(h10_0pT05_pT2,"10.5 GeV, 0 < P_{T} < 0.5","l");
        leg->AddEntry(h10_05pT15_pT2,"10.5 GeV, 0.5 < P_{T} < 1.5","l");
        leg->AddEntry(h22_0pT05_pT2,"22 GeV, 0 < P_{T} < 0.5","l");
        leg->AddEntry(h22_05pT15_pT2,"22 GeV, 0.5 < P_{T} < 1.5","l");
        leg->Draw();
    }

    // Right pad: pT1 (pion)
    c1->cd(2);
    h10_0pT05_pT1->GetXaxis()->SetTitle("x_{F}");
    h10_0pT05_pT1->GetYaxis()->SetTitle("normalized counts");
    h10_0pT05_pT1->SetMaximum(1.25 * max_val_pT1);
    // Set pad title as "pion"
    h10_0pT05_pT1->SetTitle("pion");
    h10_0pT05_pT1->Draw("hist");
    h10_05pT15_pT1->Draw("hist same");
    h22_0pT05_pT1->Draw("hist same");
    h22_05pT15_pT1->Draw("hist same");

    // Legend on right pad
    {
        TLegend *leg2 = new TLegend(0.15,0.75,0.5,0.9);
        leg2->AddEntry(h10_0pT05_pT1,"10.5 GeV, 0 < P_{T} < 0.5","l");
        leg2->AddEntry(h10_05pT15_pT1,"10.5 GeV, 0.5 < P_{T} < 1.5","l");
        leg2->AddEntry(h22_0pT05_pT1,"22 GeV, 0 < P_{T} < 0.5","l");
        leg2->AddEntry(h22_05pT15_pT1,"22 GeV, 0.5 < P_{T} < 1.5","l");
        leg2->Draw();
    }

    c1->SaveAs("output/22gev_comparison.png");

    // Cleanup
    f1->Close();
    f2->Close();

    return 0;
}