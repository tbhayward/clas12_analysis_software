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
    TH1F *h10_0pT06_pT2 = new TH1F("h10_0pT06_pT2","",50,-1.0,1.0);
    TH1F *h10_06pT16_pT2 = new TH1F("h10_06pT16_pT2","",50,-1.0,1.0);
    TH1F *h22_0pT06_pT2 = new TH1F("h22_0pT06_pT2","",50,-1.0,1.0);
    TH1F *h22_06pT16_pT2 = new TH1F("h22_06pT16_pT2","",50,-1.0,1.0);

    // Draw into histograms for pT2
    t1->Draw("xF2>>h10_0pT06_pT2","pT2>0.0 && pT2<0.6","");
    t1->Draw("xF2>>h10_06pT16_pT2","pT2>0.6 && pT2<1.6","");

    t2->Draw("xF2>>h22_0pT06_pT2","pT2>0.0 && pT2<0.6","");
    t2->Draw("xF2>>h22_06pT16_pT2","pT2>0.6 && pT2<1.6","");

    // Normalize histograms for pT2
    if (h10_0pT06_pT2->Integral() > 0) h10_0pT06_pT2->Scale(1.0/h10_0pT06_pT2->Integral());
    if (h10_06pT16_pT2->Integral() > 0) h10_06pT16_pT2->Scale(1.0/h10_06pT16_pT2->Integral());
    if (h22_0pT06_pT2->Integral() > 0) h22_0pT06_pT2->Scale(1.0/h22_0pT06_pT2->Integral());
    if (h22_06pT16_pT2->Integral() > 0) h22_06pT16_pT2->Scale(1.0/h22_06pT16_pT2->Integral());

    // Determine max for pT2 set
    double max_val_pT2 = 0;
    double temp_max = h10_0pT06_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;
    temp_max = h10_06pT16_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;
    temp_max = h22_0pT06_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;
    temp_max = h22_06pT16_pT2->GetMaximum(); if (temp_max > max_val_pT2) max_val_pT2 = temp_max;

    // Set colors and styles for pT2
    h10_0pT06_pT2->SetLineColor(kRed);   h10_0pT06_pT2->SetLineStyle(2);
    h10_06pT16_pT2->SetLineColor(kRed);  h10_06pT16_pT2->SetLineStyle(1);
    h22_0pT06_pT2->SetLineColor(kBlue);  h22_0pT06_pT2->SetLineStyle(2);
    h22_06pT16_pT2->SetLineColor(kBlue); h22_06pT16_pT2->SetLineStyle(1);

    h10_0pT06_pT2->SetStats(0);
    h10_06pT16_pT2->SetStats(0);
    h22_0pT06_pT2->SetStats(0);
    h22_06pT16_pT2->SetStats(0);

    // ---------------------------
    // Pion (pT1) histograms (Right pad)
    // ---------------------------
    TH1F *h10_0pT06_pT1 = new TH1F("h10_0pT06_pT1","",50,-1.0,1.0);
    TH1F *h10_06pT16_pT1 = new TH1F("h10_06pT16_pT1","",50,-1.0,1.0);
    TH1F *h22_0pT06_pT1 = new TH1F("h22_0pT06_pT1","",50,-1.0,1.0);
    TH1F *h22_06pT16_pT1 = new TH1F("h22_06pT16_pT1","",50,-1.0,1.0);

    // Draw into histograms for pT1
    t1->Draw("xF1>>h10_0pT06_pT1","pT1>0.0 && pT1<0.6","");
    t1->Draw("xF1>>h10_06pT16_pT1","pT1>0.6 && pT1<1.6","");

    t2->Draw("xF1>>h22_0pT06_pT1","pT1>0.0 && pT1<0.6","");
    t2->Draw("xF1>>h22_06pT16_pT1","pT1>0.6 && pT1<1.6","");

    // Normalize histograms (pT1)
    if (h10_0pT06_pT1->Integral() > 0) h10_0pT06_pT1->Scale(1.0/h10_0pT06_pT1->Integral());
    if (h10_06pT16_pT1->Integral() > 0) h10_06pT16_pT1->Scale(1.0/h10_06pT16_pT1->Integral());
    if (h22_0pT06_pT1->Integral() > 0) h22_0pT06_pT1->Scale(1.0/h22_0pT06_pT1->Integral());
    if (h22_06pT16_pT1->Integral() > 0) h22_06pT16_pT1->Scale(1.0/h22_06pT16_pT1->Integral());

    // Determine max for pT1 set
    double max_val_pT1 = 0;
    temp_max = h10_0pT06_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;
    temp_max = h10_06pT16_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;
    temp_max = h22_0pT06_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;
    temp_max = h22_06pT16_pT1->GetMaximum(); if (temp_max > max_val_pT1) max_val_pT1 = temp_max;

    // Set colors and styles for pT1
    h10_0pT06_pT1->SetLineColor(kRed);   h10_0pT06_pT1->SetLineStyle(2);
    h10_06pT16_pT1->SetLineColor(kRed);  h10_06pT16_pT1->SetLineStyle(1);
    h22_0pT06_pT1->SetLineColor(kBlue);  h22_0pT06_pT1->SetLineStyle(2);
    h22_06pT16_pT1->SetLineColor(kBlue); h22_06pT16_pT1->SetLineStyle(1);

    h10_0pT06_pT1->SetStats(0);
    h10_06pT16_pT1->SetStats(0);
    h22_0pT06_pT1->SetStats(0);
    h22_06pT16_pT1->SetStats(0);

    gStyle->SetOptStat(0);

    // Create the first canvas (1x2)
    TCanvas *c1 = new TCanvas("c1","c1",1200,600);
    c1->Divide(2,1);

    // Left pad: pT2 (proton)
    c1->cd(1);
    gPad->SetLogy();
    h10_0pT06_pT2->GetXaxis()->SetTitle("x_{F}");
    h10_0pT06_pT2->GetYaxis()->SetTitle("normalized counts");
    h10_0pT06_pT2->SetMaximum(1.25 * max_val_pT2);
    h10_0pT06_pT2->SetTitle("proton");
    h10_0pT06_pT2->Draw("hist");
    h10_06pT16_pT2->Draw("hist same");
    h22_0pT06_pT2->Draw("hist same");
    h22_06pT16_pT2->Draw("hist same");

    TLegend *leg1 = new TLegend(0.55,0.75,0.9,0.9);
    leg1->AddEntry(h10_0pT06_pT2,"10.5 GeV, 0 < P_{T} < 0.6","l");
    leg1->AddEntry(h10_06pT16_pT2,"10.5 GeV, 0.6 < P_{T} < 1.6","l");
    leg1->AddEntry(h22_0pT06_pT2,"22 GeV, 0 < P_{T} < 0.6","l");
    leg1->AddEntry(h22_06pT16_pT2,"22 GeV, 0.6 < P_{T} < 1.6","l");
    leg1->Draw();

    // Right pad: pT1 (pion)
    c1->cd(2);
    gPad->SetLogy();
    h10_0pT06_pT1->GetXaxis()->SetTitle("x_{F}");
    h10_0pT06_pT1->GetYaxis()->SetTitle("normalized counts");
    h10_0pT06_pT1->SetMaximum(1.25 * max_val_pT1);
    h10_0pT06_pT1->SetTitle("pion");
    h10_0pT06_pT1->Draw("hist");
    h10_06pT16_pT1->Draw("hist same");
    h22_0pT06_pT1->Draw("hist same");
    h22_06pT16_pT1->Draw("hist same");

    TLegend *leg2 = new TLegend(0.55,0.75,0.9,0.9);
    leg2->AddEntry(h10_0pT06_pT1,"10.5 GeV, 0 < P_{T} < 0.6","l");
    leg2->AddEntry(h10_06pT16_pT1,"10.5 GeV, 0.6 < P_{T} < 1.6","l");
    leg2->AddEntry(h22_0pT06_pT1,"22 GeV, 0 < P_{T} < 0.6","l");
    leg2->AddEntry(h22_06pT16_pT1,"22 GeV, 0.6 < P_{T} < 1.6","l");
    leg2->Draw();

    c1->SaveAs("output/22gev_comparison.png");

    // Create the second canvas (ratios)
    TCanvas *c2 = new TCanvas("c2","c2",1200,600);
    c2->Divide(2,1);

    // Create and plot ratio histograms for pT2
    TH1F *ratio_low_pT_pT2 = (TH1F*)h22_0pT06_pT2->Clone("ratio_low_pT_pT2");
    ratio_low_pT_pT2->Divide(h10_0pT06_pT2);
    TH1F *ratio_high_pT_pT2 = (TH1F*)h22_06pT16_pT2->Clone("ratio_high_pT_pT2");
    ratio_high_pT_pT2->Divide(h10_06pT16_pT2);

    ratio_low_pT_pT2->SetLineColor(kBlack);
    ratio_low_pT_pT2->SetLineStyle(2);
    ratio_high_pT_pT2->SetLineColor(kBlack);
    ratio_high_pT_pT2->SetLineStyle(1);

    double max_ratio_pT2 = std::max(ratio_low_pT_pT2->GetMaximum(), ratio_high_pT_pT2->GetMaximum()) * 1.25;

    c2->cd(1);
    ratio_low_pT_pT2->GetXaxis()->SetTitle("x_{F}");
    ratio_low_pT_pT2->GetYaxis()->SetTitle("ratio");
    ratio_low_pT_pT2->SetMaximum(max_ratio_pT2);
    ratio_low_pT_pT2->Draw("hist");
    ratio_high_pT_pT2->Draw("hist same");

    TLegend *leg_ratio_pT2 = new TLegend(0.6,0.75,0.9,0.9);
    leg_ratio_pT2->AddEntry(ratio_low_pT_pT2,"0 < P_{T} < 0.6","l");
    leg_ratio_pT2->AddEntry(ratio_high_pT_pT2,"0.6 < P_{T} < 1.6","l");
    leg_ratio_pT2->Draw();

    // Create and plot ratio histograms for pT1
    TH1F *ratio_low_pT_pT1 = (TH1F*)h22_0pT06_pT1->Clone("ratio_low_pT_pT1");
    ratio_low_pT_pT1->Divide(h10_0pT06_pT1);
    TH1F *ratio_high_pT_pT1 = (TH1F*)h22_06pT16_pT1->Clone("ratio_high_pT_pT1");
    ratio_high_pT_pT1->Divide(h10_06pT16_pT1);

    ratio_low_pT_pT1->SetLineColor(kBlack);
    ratio_low_pT_pT1->SetLineStyle(2);
    ratio_high_pT_pT1->SetLineColor(kBlack);
    ratio_high_pT_pT1->SetLineStyle(1);

    double max_ratio_pT1 = std::max(ratio_low_pT_pT1->GetMaximum(), ratio_high_pT_pT1->GetMaximum()) * 1.25;

    c2->cd(2);
    ratio_low_pT_pT1->GetXaxis()->SetTitle("x_{F}");
    ratio_low_pT_pT1->GetYaxis()->SetTitle("ratio");
    ratio_low_pT_pT1->SetMaximum(max_ratio_pT1);
    ratio_low_pT_pT1->Draw("hist");
    ratio_high_pT_pT1->Draw("hist same");

    TLegend *leg_ratio_pT1 = new TLegend(0.6,0.75,0.9,0.9);
    leg_ratio_pT1->AddEntry(ratio_low_pT_pT1,"0 < P_{T} < 0.6","l");
    leg_ratio_pT1->AddEntry(ratio_high_pT_pT1,"0.6 < P_{T} < 1.6","l");
    leg_ratio_pT1->Draw();

    c2->SaveAs("output/22gev_ratios.png");

    // Cleanup
    f1->Close();
    f2->Close();

    return 0;
}