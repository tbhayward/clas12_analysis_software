#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>

void plotRatios(const char* file1, const char* file2, const char* file3, const char* file4) {
    // Open the ROOT files
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TFile *f3 = TFile::Open(file3);
    TFile *f4 = TFile::Open(file4);

    if (!f1 || !f2 || !f3 || !f4 || f1->IsZombie() || f2->IsZombie() || f3->IsZombie() || f4->IsZombie()) {
        std::cerr << "Error opening files" << std::endl;
        return;
    }

    // Get the trees
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");
    TTree *tree3 = (TTree*)f3->Get("PhysicsEvents");
    TTree *tree4 = (TTree*)f4->Get("PhysicsEvents");

    if (!tree1 || !tree2 || !tree3 || !tree4) {
        std::cerr << "Error: Could not find the PhysicsEvents tree in one of the files" << std::endl;
        return;
    }

    const int nbins = 30;
    TH1F *h1_p_p = new TH1F("h1_p_p", "p_p", nbins, 0, 4);
    TH1F *h2_p_p = new TH1F("h2_p_p", "p_p", nbins, 0, 4);
    TH1F *h3_p_p = new TH1F("h3_p_p", "p_p", nbins, 0, 4);
    TH1F *h4_p_p = new TH1F("h4_p_p", "p_p", nbins, 0, 4);

    TH1F *h1_xF = new TH1F("h1_xF", "xF", nbins, -2, 1);
    TH1F *h2_xF = new TH1F("h2_xF", "xF", nbins, -2, 1);
    TH1F *h3_xF = new TH1F("h3_xF", "xF", nbins, -2, 1);
    TH1F *h4_xF = new TH1F("h4_xF", "xF", nbins, -2, 1);

    TH1F *h1_p_theta = new TH1F("h1_p_theta", "p_theta", nbins, 0, 30);
    TH1F *h2_p_theta = new TH1F("h2_p_theta", "p_theta", nbins, 0, 30);
    TH1F *h3_p_theta = new TH1F("h3_p_theta", "p_theta", nbins, 0, 30);
    TH1F *h4_p_theta = new TH1F("h4_p_theta", "p_theta", nbins, 0, 30);

    TH1F *h1_Mx = new TH1F("h1_Mx", "Mx", nbins, -2, 3);
    TH1F *h2_Mx = new TH1F("h2_Mx", "Mx", nbins, -2, 3);
    TH1F *h3_Mx = new TH1F("h3_Mx", "Mx", nbins, -2, 3);
    TH1F *h4_Mx = new TH1F("h4_Mx", "Mx", nbins, -2, 3);

    double p_p1, p_p2, p_p3, p_p4;
    double xF1, xF2, xF3, xF4;
    double p_theta1, p_theta2, p_theta3, p_theta4;
    double Mx1, Mx2, Mx3, Mx4;

    tree1->SetBranchAddress("p_p", &p_p1);
    tree2->SetBranchAddress("p_p", &p_p2);
    tree3->SetBranchAddress("p_p", &p_p3);
    tree4->SetBranchAddress("p_p", &p_p4);

    tree1->SetBranchAddress("xF", &xF1);
    tree2->SetBranchAddress("xF", &xF2);
    tree3->SetBranchAddress("xF", &xF3);
    tree4->SetBranchAddress("xF", &xF4);

    tree1->SetBranchAddress("p_theta", &p_theta1);
    tree2->SetBranchAddress("p_theta", &p_theta2);
    tree3->SetBranchAddress("p_theta", &p_theta3);
    tree4->SetBranchAddress("p_theta", &p_theta4);

    tree1->SetBranchAddress("Mx", &Mx1);
    tree2->SetBranchAddress("Mx", &Mx2);
    tree3->SetBranchAddress("Mx", &Mx3);
    tree4->SetBranchAddress("Mx", &Mx4);

    Long64_t nentries1 = tree1->GetEntries();
    Long64_t nentries2 = tree2->GetEntries();
    Long64_t nentries3 = tree3->GetEntries();
    Long64_t nentries4 = tree4->GetEntries();

    // Fill histograms for the first set of files
    for (Long64_t j = 0; j < nentries1; ++j) {
        tree1->GetEntry(j);
        // if (Mx1 > 1.4) {
            h1_p_p->Fill(p_p1);
            h1_xF->Fill(xF1);
            h1_p_theta->Fill(p_theta1 * 180.0 / 3.14159);
        // }
        h1_Mx->Fill(Mx1); // Always fill Mx histogram
    }

    for (Long64_t j = 0; j < nentries2; ++j) {
        tree2->GetEntry(j);
        // if (Mx2 > 1.4) {
            h2_p_p->Fill(p_p2);
            h2_xF->Fill(xF2);
            h2_p_theta->Fill(p_theta2 * 180.0 / 3.14159);
        // }
        h2_Mx->Fill(Mx2); // Always fill Mx histogram
    }

    // Fill histograms for the second set of files
    for (Long64_t j = 0; j < nentries3; ++j) {
        tree3->GetEntry(j);
        // if (Mx3 > 1.4) {
            h3_p_p->Fill(p_p3);
            h3_xF->Fill(xF3);
            h3_p_theta->Fill(p_theta3 * 180.0 / 3.14159);
        // }
        h3_Mx->Fill(Mx3); // Always fill Mx histogram
    }

    for (Long64_t j = 0; j < nentries4; ++j) {
        tree4->GetEntry(j);
        // if (Mx4 > 1.4) {
            h4_p_p->Fill(p_p4);
            h4_xF->Fill(xF4);
            h4_p_theta->Fill(p_theta4 * 180.0 / 3.14159);
        // }
        h4_Mx->Fill(Mx4); // Always fill Mx histogram
    }

    gStyle->SetOptStat(0);

    // Create ratio histograms
    TH1F *ratio_p_p_1 = (TH1F*)h1_p_p->Clone("ratio_p_p_1");
    ratio_p_p_1->Divide(h2_p_p);
    ratio_p_p_1->GetXaxis()->SetTitle("p_p");

    TH1F *ratio_p_p_2 = (TH1F*)h3_p_p->Clone("ratio_p_p_2");
    ratio_p_p_2->Divide(h4_p_p);

    TH1F *ratio_xF_1 = (TH1F*)h1_xF->Clone("ratio_xF_1");
    ratio_xF_1->Divide(h2_xF);
    ratio_xF_1->GetXaxis()->SetTitle("xF");

    TH1F *ratio_xF_2 = (TH1F*)h3_xF->Clone("ratio_xF_2");
    ratio_xF_2->Divide(h4_xF);

    TH1F *ratio_p_theta_1 = (TH1F*)h1_p_theta->Clone("ratio_p_theta_1");
    ratio_p_theta_1->Divide(h2_p_theta);
    ratio_p_theta_1->GetXaxis()->SetTitle("p_theta (degrees)");

    TH1F *ratio_p_theta_2 = (TH1F*)h3_p_theta->Clone("ratio_p_theta_2");
    ratio_p_theta_2->Divide(h4_p_theta);

    TH1F *ratio_Mx_1 = (TH1F*)h1_Mx->Clone("ratio_Mx_1");
    ratio_Mx_1->Divide(h2_Mx);
    ratio_Mx_1->GetXaxis()->SetTitle("Mx");

    TH1F *ratio_Mx_2 = (TH1F*)h3_Mx->Clone("ratio_Mx_2");
    ratio_Mx_2->Divide(h4_Mx);

    // Save the ratio histograms
    TCanvas *c = new TCanvas("c", "c", 800, 600);

    // Plot ratios for p_p
    c->cd();
    ratio_p_p_1->SetLineColor(kRed);
    ratio_p_p_1->Draw();
    ratio_p_p_2->SetLineColor(kBlue);
    ratio_p_p_2->Draw("SAME");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    ratio_p_p_1->SetMinimum(1.4);
    ratio_p_p_1->SetMaximum(2.0);
    legend->AddEntry(ratio_p_p_1, "preliminary", "l");
    legend->AddEntry(ratio_p_p_2, "pass-1", "l");
    legend->Draw();

    c->SaveAs("output/ratio_p_p.png");

    // Plot ratios for xF
    c->cd();
    ratio_xF_1->SetLineColor(kRed);
    ratio_xF_1->Draw();
    ratio_xF_2->SetLineColor(kBlue);
    ratio_xF_2->Draw("SAME");

    legend->Clear();
    legend->AddEntry(ratio_xF_1, "preliminary", "l");
    legend->AddEntry(ratio_xF_2, "pass-1", "l");
    legend->Draw();

    c->SaveAs("output/ratio_xF.png");

    // Plot ratios for p_theta
    c->cd();
    ratio_p_theta_1->SetLineColor(kRed);
    ratio_p_theta_1->Draw();
    ratio_p_theta_2->SetLineColor(kBlue);
    ratio_p_theta_2->Draw("SAME");

    legend->Clear();
    legend->AddEntry(ratio_p_theta_1, "preliminary", "l");
    legend->AddEntry(ratio_p_theta_2, "pass-1", "l");
    legend->Draw();

    c->SaveAs("output/ratio_p_theta.png");

    // Plot ratios for Mx
    c->cd();
    ratio_Mx_1->SetMinimum(1.4);
    ratio_Mx_1->SetMaximum(2.0);
    ratio_Mx_1->SetLineColor(kRed);
    ratio_Mx_1->Draw();
    ratio_Mx_2->SetLineColor(kBlue);
    ratio_Mx_2->Draw("SAME");

    legend->Clear();
    legend->AddEntry(ratio_Mx_1, "preliminary", "l");
    legend->AddEntry(ratio_Mx_2, "pass-1", "l");
    legend->Draw();

    c->SaveAs("output/ratio_Mx.png");

    // Clean up
    delete c;
    delete h1_p_p;
    delete h2_p_p;
    delete h3_p_p;
    delete h4_p_p;
    delete h1_xF;
    delete h2_xF;
    delete h3_xF;
    delete h4_xF;
    delete h1_p_theta;
    delete h2_p_theta;
    delete h3_p_theta;
    delete h4_p_theta;
    delete ratio_p_p_1;
    delete ratio_p_p_2;
    delete ratio_xF_1;
    delete ratio_xF_2;
    delete ratio_p_theta_1;
    delete ratio_p_theta_2;

    f1->Close();
    f2->Close();
    f3->Close();
    f4->Close();
    }

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " file1.root file2.root file3.root file4.root" << std::endl;
        return 1;
    }
    plotRatios(argv[1], argv[2], argv[3], argv[4]);
    return 0;
}
   