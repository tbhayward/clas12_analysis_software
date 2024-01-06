#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <iostream>

void pair_production_rate(const char* file1, const char* file2, const char* output) {
    // Open ROOT files and get trees
    TFile* f1 = new TFile(file1);
    TFile* f2 = new TFile(file2);
    TTree* tree1 = (TTree*)f1->Get("PhysicsEvents"); // Adjust tree name if necessary
    TTree* tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Define histograms for each variable
    TH1D* h_e_p1 = new TH1D("h_e_p1", "e_p", 100, 2.6, 10);
    TH1D* h_e_p2 = new TH1D("h_e_p2", "e_p", 100, 2.6, 10);
    TH1D* h_Q21 = new TH1D("h_Q21", "Q^2", 100, 0, 10); // Adjust range as needed
    TH1D* h_Q22 = new TH1D("h_Q22", "Q^2", 100, 0, 10);
    TH1D* h_W1 = new TH1D("h_W1", "W", 100, 0, 3); // Adjust range as needed
    TH1D* h_W2 = new TH1D("h_W2", "W", 100, 0, 3);

    // // Set branch addresses
    // double e_p, Q2, W;
    // tree1->SetBranchAddress("e_p", &e_p);
    // tree1->SetBranchAddress("Q2", &Q2);
    // tree1->SetBranchAddress("W", &W);
    // tree2->SetBranchAddress("e_p", &e_p);
    // tree2->SetBranchAddress("Q2", &Q2);
    // tree2->SetBranchAddress("W", &W);

    // // Fill histograms
    // Long64_t nEntries1 = tree1->GetEntries();
    // Long64_t nEntries2 = tree2->GetEntries();
    // for (Long64_t i = 0; i < nEntries1; ++i) {
    //     tree1->GetEntry(i);
    //     h_e_p1->Fill(e_p);
    //     h_Q21->Fill(Q2);
    //     h_W1->Fill(W);
    // }
    // for (Long64_t i = 0; i < nEntries2; ++i) {
    //     tree2->GetEntry(i);
    //     h_e_p2->Fill(e_p);
    //     h_Q22->Fill(Q2);
    //     h_W2->Fill(W);
    // }

    // // Normalize histograms
    // h_e_p1->Scale(1.0 / 378709);
    // h_e_p2->Scale(1.0 / 266653);
    // h_Q21->Scale(1.0 / 378709);
    // h_Q22->Scale(1.0 / 266653);
    // h_W1->Scale(1.0 / 378709);
    // h_W2->Scale(1.0 / 266653);

    // // Create a canvas
    // TCanvas* c1 = new TCanvas("c1", "Pair Production Rate", 1200, 800);
    // c1->Divide(2, 3);

    // // Draw histograms and ratios
    // // e_p
    // c1->cd(1);
    // gPad->SetLogy(1);
    // h_e_p1->SetLineColor(kBlue);
    // h_e_p1->Draw();
    // h_e_p2->SetLineColor(kRed);
    // h_e_p2->Draw("same");

    // c1->cd(2);
    // gPad->SetLogy(1);
    // TH1D* h_ratio_e_p = (TH1D*)h_e_p2->Clone();
    // h_ratio_e_p->Divide(h_e_p1);
    // h_ratio_e_p->Draw();

    // // Q^2
    // c1->cd(3);
    // gPad->SetLogy(1);
    // h_Q21->Draw();
    // h_Q22->Draw("same");

    // c1->cd(4);
    // gPad->SetLogy(1);
    // TH1D* h_ratio_Q2 = (TH1D*)h_Q22->Clone();
    // h_ratio_Q2->Divide(h_Q21);
    // h_ratio_Q2->Draw();

    // // W
    // c1->cd(5);
    // gPad->SetLogy(1);
    // h_W1->Draw();
    // h_W2->Draw("same");

    // c1->cd(6);
    // gPad->SetLogy(1);
    // TH1D* h_ratio_W = (TH1D*)h_W2->Clone();
    // h_ratio_W->Divide(h_W1);
    // h_ratio_W->Draw();

    // // Add legends
    // TLegend* leg1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg1->AddEntry(h_e_p1, "run 5038, e^-", "l");
    // leg1->AddEntry(h_e_p2, "run 5482, e^+", "l");
    // leg1->Draw();

    // TLegend* leg2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg2->AddEntry(h_Q21, "run 5038, e^-", "l");
    // leg2->AddEntry(h_Q22, "run 5482, e^+", "l");
    // leg2->Draw();

    // TLegend* leg3 = new TLegend(0.7, 0.7, 0.9, 0.9);
    // leg3->AddEntry(h_W1, "run 5038, e^-", "l");
    // leg3->AddEntry(h_W2, "run 5482, e^+", "l");
    // leg3->Draw();

    // // Save the canvas
    // c1->SaveAs(output);

    // // Clean up
    // delete f1;
    // delete f2;
    // delete h_e_p1;
    // delete h_e_p2;
    // delete h_Q21;
    // delete h_Q22;
    // delete h_W1;
    // delete h_W2;
    // delete h_ratio_e_p;
    // delete h_ratio_Q2;
    // delete h_ratio_W;
    // delete leg1;
    // delete leg2;
    // delete leg3;
    // delete c1;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <file1> <file2> <output>" << std::endl;
        return 1;
    }
    pair_production_rate(argv[1], argv[2], argv[3]);
    return 0;
}
