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
    TH1D* h_e_p1 = new TH1D("h_e_p1", ";e_{p} (GeV);Normalized Counts", 100, 2.6, 10);
    TH1D* h_e_p2 = new TH1D("h_e_p2", ";e_{p} (GeV);Normalized Counts", 100, 2.6, 10);
    TH1D* h_Q21 = new TH1D("h_Q21", ";Q^{2} (GeV^{2});Normalized Counts", 100, 0, 10); 
    TH1D* h_Q22 = new TH1D("h_Q22", ";Q^{2} (GeV^{2});Normalized Counts", 100, 0, 10);
    TH1D* h_W1 = new TH1D("h_W1", ";W (GeV);Normalized Counts", 100, 0.8, 4); 
    TH1D* h_W2 = new TH1D("h_W2", ";W (GeV);Normalized Counts", 100, 0.8, 4);

    // Set label and title sizes for y-axis for all histograms
    double labelFontSize = 0.05; // Adjust as needed
    double titleFontSize = 0.06; // Adjust as needed

    h_e_p1->GetYaxis()->SetLabelSize(labelFontSize);
    h_e_p1->GetYaxis()->SetTitleSize(titleFontSize);
    h_e_p2->GetYaxis()->SetLabelSize(labelFontSize);
    h_e_p2->GetYaxis()->SetTitleSize(titleFontSize);
    h_Q21->GetYaxis()->SetLabelSize(labelFontSize);
    h_Q21->GetYaxis()->SetTitleSize(titleFontSize);
    h_Q22->GetYaxis()->SetLabelSize(labelFontSize);
    h_Q22->GetYaxis()->SetTitleSize(titleFontSize);
    h_W1->GetYaxis()->SetLabelSize(labelFontSize);
    h_W1->GetYaxis()->SetTitleSize(titleFontSize);
    h_W2->GetYaxis()->SetLabelSize(labelFontSize);
    h_W2->GetYaxis()->SetTitleSize(titleFontSize);

    // Set branch addresses
    double e_p, Q2, W;
    tree1->SetBranchAddress("e_p", &e_p);
    tree1->SetBranchAddress("Q2", &Q2);
    tree1->SetBranchAddress("W", &W);
    tree2->SetBranchAddress("e_p", &e_p);
    tree2->SetBranchAddress("Q2", &Q2);
    tree2->SetBranchAddress("W", &W);

    // Fill histograms
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        h_e_p1->Fill(e_p);
        h_Q21->Fill(Q2);
        h_W1->Fill(W);
    }
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        h_e_p2->Fill(e_p);
        h_Q22->Fill(Q2);
        h_W2->Fill(W);
    }

    // Normalize histograms
    h_e_p1->Scale(1.0 / 378709);
    h_e_p2->Scale(1.0 / 266653);
    h_Q21->Scale(1.0 / 378709);
    h_Q22->Scale(1.0 / 266653);
    h_W1->Scale(1.0 / 378709);
    h_W2->Scale(1.0 / 266653);

    // Define larger font sizes
    double legendFontSize = 0.05; // Adjust as needed
    double leftMargin = 0.15; // Increase left margin
    double bottomMargin = 0.15; // Increase bottom margin

    // Create a canvas
    TCanvas* c1 = new TCanvas("c1", "Pair Production Rate", 1200, 800);
    c1->Divide(2, 3);

    // Draw histograms and ratios
    // e_p
    c1->cd(1);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(leftMargin);
    gPad->SetBottomMargin(bottomMargin);
    h_e_p1->SetLineColor(kBlue);
    h_e_p1->SetLabelSize(labelFontSize);
    h_e_p1->SetTitleSize(titleFontSize);
    h_e_p1->Draw();
    h_e_p2->SetLineColor(kRed);
    h_e_p2->SetLabelSize(labelFontSize);
    h_e_p2->SetTitleSize(titleFontSize);
    h_e_p2->Draw("same");
    h_e_p1->SetMaximum(10); // 100% higher than the maximum
    h_e_p1->SetMinimum(1e-5); // Minimum set to 10e-5
    h_e_p1->SetStats(0); // Remove stat box
    h_e_p2->SetStats(0); // Remove stat box
    h_e_p1->GetXaxis()->SetLabelSize(labelFontSize);
    h_e_p1->GetXaxis()->SetTitleSize(titleFontSize); // Increase x-axis title size
    h_e_p1->GetYaxis()->SetLabelSize(labelFontSize);
    h_e_p1->GetYaxis()->SetTitleSize(titleFontSize); // Increase y-axis title size


    c1->cd(2);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(leftMargin);
    gPad->SetBottomMargin(bottomMargin);
    TH1D* h_ratio_e_p = (TH1D*)h_e_p2->Clone();
    h_ratio_e_p->Divide(h_e_p1);
    h_ratio_e_p->SetLineColor(kBlack);
    h_ratio_e_p->SetTitle(";e_{p} (GeV);Ratio");
    h_ratio_e_p->SetLabelSize(labelFontSize);
    h_ratio_e_p->SetTitleSize(titleFontSize);
    h_ratio_e_p->GetYaxis()->SetLabelSize(labelFontSize);
    h_ratio_e_p->GetYaxis()->SetTitleSize(titleFontSize); // Increase y-axis title size
    h_ratio_e_p->Draw();
    h_ratio_e_p->SetMaximum(10.0); // 100% higher than the maximum
    h_ratio_e_p->SetMinimum(1e-5); // Minimum set to 10e-5
    h_ratio_e_p->SetStats(0); // Remove stat box


    // Q^2
    c1->cd(3);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(leftMargin);
    gPad->SetBottomMargin(bottomMargin);
    h_Q21->Draw();
    h_Q21->SetLineColor(kBlue);
    h_Q21->SetLabelSize(labelFontSize);
    h_Q21->SetTitleSize(titleFontSize);
    h_Q22->Draw("same");
    h_Q22->SetLineColor(kRed);
    h_Q22->SetLabelSize(labelFontSize);
    h_Q22->SetTitleSize(titleFontSize);
    h_Q21->SetMaximum(10.0); // 100% higher than the maximum
    h_Q21->SetMinimum(1e-5); // Minimum set to 10e-5
    h_Q21->SetStats(0); // Remove stat box
    h_Q22->SetStats(0); // Remove stat box
    h_Q21->GetXaxis()->SetLabelSize(labelFontSize);
    h_Q21->GetXaxis()->SetTitleSize(titleFontSize); // Increase x-axis title size
    h_Q22->GetYaxis()->SetLabelSize(labelFontSize);
    h_Q22->GetYaxis()->SetTitleSize(titleFontSize); // Increase y-axis title size

    c1->cd(4);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(leftMargin);
    gPad->SetBottomMargin(bottomMargin);
    TH1D* h_ratio_Q2 = (TH1D*)h_Q22->Clone();
    h_ratio_Q2->Divide(h_Q21);
    h_ratio_Q2->SetTitle(";Q^{2} (GeV^2);Ratio");
    h_ratio_Q2->SetLineColor(kBlack);
    h_ratio_Q2->SetLabelSize(labelFontSize);
    h_ratio_Q2->SetTitleSize(titleFontSize);
    h_ratio_Q2->Draw();
    h_ratio_Q2->SetMaximum(10.0); // 100% higher than the maximum
    h_ratio_Q2->SetMinimum(1e-5); // Minimum set to 10e-5
    h_ratio_Q2->SetStats(0); // Remove stat box
    h_ratio_Q2->GetXaxis()->SetLabelSize(labelFontSize);
    h_ratio_Q2->GetXaxis()->SetTitleSize(titleFontSize); // Increase x-axis title size

    // W
    c1->cd(5);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(leftMargin);
    gPad->SetBottomMargin(bottomMargin);
    h_W1->SetLabelSize(labelFontSize);
    h_W1->SetTitleSize(titleFontSize);
    h_W1->Draw();
    h_W1->SetLineColor(kBlue);
    h_W2->SetLabelSize(labelFontSize);
    h_W2->SetTitleSize(titleFontSize);
    h_W2->Draw("same");
    h_W2->SetLineColor(kRed);
    h_W1->SetMaximum(10.0); // 100% higher than the maximum
    h_W1->SetMinimum(1e-5); // Minimum set to 10e-5
    h_W1->SetStats(0); // Remove stat box
    h_W2->SetStats(0); // Remove stat box
    h_W1->GetXaxis()->SetLabelSize(labelFontSize);
    h_W1->GetXaxis()->SetTitleSize(titleFontSize); // Increase x-axis title size
    h_W2->GetYaxis()->SetLabelSize(labelFontSize);
    h_W2->GetYaxis()->SetTitleSize(titleFontSize); // Increase y-axis title size

    c1->cd(6);
    gPad->SetLogy(1);
    gPad->SetLeftMargin(leftMargin);
    gPad->SetBottomMargin(bottomMargin);
    TH1D* h_ratio_W = (TH1D*)h_W2->Clone();
    h_ratio_W->Divide(h_W1);
    h_ratio_W->SetTitle(";W (GeV);Ratio");
    h_ratio_W->SetLineColor(kBlack);
    h_ratio_W->SetLabelSize(labelFontSize);
    h_ratio_W->SetTitleSize(titleFontSize);
    h_ratio_W->Draw();
    h_ratio_W->SetMaximum(10.0); // 100% higher than the maximum
    h_ratio_W->SetMinimum(1e-5); // Minimum set to 10e-5
    h_ratio_W->SetStats(0); // Remove stat box
    h_ratio_W->GetXaxis()->SetLabelSize(labelFontSize);
    h_ratio_W->GetXaxis()->SetTitleSize(titleFontSize); // Increase x-axis title size

    // Add legends to the left-side plots
    TLegend* leg_e_p = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_e_p->AddEntry(h_e_p1, "run 5038, e^{-}", "l");
    leg_e_p->AddEntry(h_e_p2, "run 5482, e^{+}", "l");
    c1->cd(1);
    leg_e_p->Draw();

    TLegend* leg_Q2 = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_Q2->AddEntry(h_Q21, "run 5038, e^{-}", "l");
    leg_Q2->AddEntry(h_Q22, "run 5482, e^{+}", "l");
    c1->cd(3);
    leg_Q2->Draw();

    TLegend* leg_W = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg_W->AddEntry(h_W1, "run 5038, e^{-}", "l");
    leg_W->AddEntry(h_W2, "run 5482, e^{+}", "l");
    c1->cd(5);
    leg_W->Draw();

    // Save the canvas
    c1->SaveAs(output);

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
