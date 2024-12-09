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

    // ---------------------------
    // Proton (pT2) histograms (Left pad)
    // ---------------------------
    TH1F *h10_0pT05_pT2 = new TH1F("h10_0pT05_pT2","",50,-1.5,1.5);
    TH1F *h10_05pT15_pT2 = new TH1F("h10_05pT15_pT2","",50,-1.5,1.5);
    TH1F *h22_0pT05_pT2 = new TH1F("h22_0pT05_pT2","",50,-1.5,1.5);
    TH1F *h22_05pT15_pT2 = new TH1F("h22_05pT15_pT2","",50,-1.5,1.5);

    // Draw into histograms for pT2
    t1->Draw("xF2>>h10_0pT05_pT2","pT2>0 && pT2<0.5","");
    t1->Draw("xF2>>h10_05pT15_pT2","pT2>0.5 && pT2<1.5","");

    t2->Draw("xF2>>h22_0pT05_pT2","pT2>0 && pT2<0.5","");
    t2->Draw("xF2>>h22_05pT15_pT2","pT2>0.5 && pT2<1.5","");

    // Normalize histograms for pT2
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
    TH1F *h10_0pT05_pT1 = new TH1F("h10_0pT05_pT1","",50,-1.5,1.5);
    TH1F *h10_05pT15_pT1 = new TH1F("h10_05pT15_pT1","",50,-1.5,1.5);
    TH1F *h22_0pT05_pT1 = new TH1F("h22_0pT05_pT1","",50,-1.5,1.5);
    TH1F *h22_05pT15_pT1 = new TH1F("h22_05pT15_pT1","",50,-1.5,1.5);

    // Draw into histograms for pT1
    t1->Draw("xF1>>h10_0pT05_pT1","pT1>0 && pT1<0.5","");
    t1->Draw("xF1>>h10_05pT15_pT1","pT1>0.5 && pT1<1.5","");

    t2->Draw("xF1>>h22_0pT05_pT1","pT1>0 && pT1<0.5","");
    t2->Draw("xF1>>h22_05pT15_pT1","pT1>0.5 && pT1<1.5","");

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

    // Create a 1x2 canvas (original functionality)
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1); // 1 row, 2 columns

    // Left pad: pT2 (proton)
    c1->cd(1);
    gPad->SetLogy(); // Set log scale for left pad
    h10_0pT05_pT2->GetXaxis()->SetTitle("x_{F}");
    h10_0pT05_pT2->GetYaxis()->SetTitle("normalized counts");
    h10_0pT05_pT2->SetMaximum(5 * max_val_pT2);
    h10_0pT05_pT2->SetTitle("proton");
    h10_0pT05_pT2->Draw("hist");
    h10_05pT15_pT2->Draw("hist same");
    h22_0pT05_pT2->Draw("hist same");
    h22_05pT15_pT2->Draw("hist same");

    {
        TLegend *leg = new TLegend(0.7,0.75,0.9,0.9);
        leg->AddEntry(h10_0pT05_pT2,"10.5 GeV, 0 < P_{T} < 0.5","l");
        leg->AddEntry(h10_05pT15_pT2,"10.5 GeV, 0.5 < P_{T} < 1.5","l");
        leg->AddEntry(h22_0pT05_pT2,"22 GeV, 0 < P_{T} < 0.5","l");
        leg->AddEntry(h22_05pT15_pT2,"22 GeV, 0.5 < P_{T} < 1.5","l");
        leg->Draw();
    }

    // Right pad: pT1 (pion)
    c1->cd(2);
    gPad->SetLogy(); // Set log scale for right pad
    h10_0pT05_pT1->GetXaxis()->SetTitle("x_{F}");
    h10_0pT05_pT1->GetYaxis()->SetTitle("normalized counts");
    h10_0pT05_pT1->SetMaximum(5 * max_val_pT1);
    h10_0pT05_pT1->SetTitle("pion");
    h10_0pT05_pT1->Draw("hist");
    h10_05pT15_pT1->Draw("hist same");
    h22_0pT05_pT1->Draw("hist same");
    h22_05pT15_pT1->Draw("hist same");

    {
        TLegend *leg2 = new TLegend(0.7,0.75,0.9,0.9);
        leg2->AddEntry(h10_0pT05_pT1,"10.5 GeV, 0 < P_{T} < 0.5","l");
        leg2->AddEntry(h10_05pT15_pT1,"10.5 GeV, 0.5 < P_{T} < 1.5","l");
        leg2->AddEntry(h22_0pT05_pT1,"22 GeV, 0 < P_{T} < 0.5","l");
        leg2->AddEntry(h22_05pT15_pT1,"22 GeV, 0.5 < P_{T} < 1.5","l");
        leg2->Draw();
    }

    c1->SaveAs("output/22gev_comparison.png");

    // Now create a new canvas for the ratio
    // Ratio = (22 GeV histogram) / (10.5 GeV histogram)
    // We'll do this for pT2 (proton) as an example, giving two lines (low pT and high pT).

    // Create ratio histograms for pT2
    TH1F *ratio_low_pT_pT2 = (TH1F*)h22_0pT05_pT2->Clone("ratio_low_pT_pT2");
    ratio_low_pT_pT2->Divide(h10_0pT05_pT2);

    TH1F *ratio_high_pT_pT2 = (TH1F*)h22_05pT15_pT2->Clone("ratio_high_pT_pT2");
    ratio_high_pT_pT2->Divide(h10_05pT15_pT2);

    // Set styles: black line for high pT, dashed black line for low pT
    ratio_low_pT_pT2->SetLineColor(kBlack);
    ratio_low_pT_pT2->SetLineStyle(2);  // dashed
    ratio_high_pT_pT2->SetLineColor(kBlack);
    ratio_high_pT_pT2->SetLineStyle(1); // solid

    ratio_low_pT_pT2->SetStats(0);
    ratio_high_pT_pT2->SetStats(0);

    // Determine max for ratio
    double max_ratio = 0;
    temp_max = ratio_low_pT_pT2->GetMaximum(); if(temp_max > max_ratio) max_ratio = temp_max;
    temp_max = ratio_high_pT_pT2->GetMaximum(); if(temp_max > max_ratio) max_ratio = temp_max;

    // Create second canvas for ratio
    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    gPad->SetLogy();
    // "with all similar aesthetics" - we keep axis titles, etc.
    ratio_low_pT_pT2->GetXaxis()->SetTitle("x_{F}");
    ratio_low_pT_pT2->GetYaxis()->SetTitle("ratio");
    ratio_low_pT_pT2->SetMaximum(1.25 * max_ratio);

    // Draw ratio histograms
    ratio_low_pT_pT2->Draw("hist");
    ratio_high_pT_pT2->Draw("hist same");

    // Legend: just say which pT bin is which
    {
        TLegend *leg_ratio = new TLegend(0.7,0.75,0.9,0.9);
        leg_ratio->AddEntry(ratio_low_pT_pT2,"0 < P_{T} < 0.5","l");
        leg_ratio->AddEntry(ratio_high_pT_pT2,"0.5 < P_{T} < 1.5","l");
        leg_ratio->Draw();
    }

    c2->SaveAs("output/22gev_ratio.png");

    // Cleanup
    f1->Close();
    f2->Close();

    return 0;
}