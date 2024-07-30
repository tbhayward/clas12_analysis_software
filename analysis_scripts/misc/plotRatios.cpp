#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TString.h>
#include <TSystem.h>
#include <iostream>
#include <TGraphErrors.h>
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

    TH1F *h1_p_theta = new TH1F("h1_p_theta", "p_theta", nbins, 0, 80);
    TH1F *h2_p_theta = new TH1F("h2_p_theta", "p_theta", nbins, 0, 80);
    TH1F *h3_p_theta = new TH1F("h3_p_theta", "p_theta", nbins, 0, 80);
    TH1F *h4_p_theta = new TH1F("h4_p_theta", "p_theta", nbins, 0, 80);

    TH1F *h1_Mx = new TH1F("h1_Mx", "Mx", nbins, -3, 4);
    TH1F *h2_Mx = new TH1F("h2_Mx", "Mx", nbins, -3, 4);
    TH1F *h3_Mx = new TH1F("h3_Mx", "Mx", nbins, -3, 4);
    TH1F *h4_Mx = new TH1F("h4_Mx", "Mx", nbins, -3, 4);

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
        h1_p_p->Fill(p_p1);
        h1_xF->Fill(xF1);
        h1_p_theta->Fill(p_theta1 * 180.0 / 3.14159);
        h1_Mx->Fill(Mx1); // Always fill Mx histogram
    }

    for (Long64_t j = 0; j < nentries2; ++j) {
        tree2->GetEntry(j);
        h2_p_p->Fill(p_p2);
        h2_xF->Fill(xF2);
        h2_p_theta->Fill(p_theta2 * 180.0 / 3.14159);
        h2_Mx->Fill(Mx2); // Always fill Mx histogram
    }

    // Fill histograms for the second set of files
    for (Long64_t j = 0; j < nentries3; ++j) {
        tree3->GetEntry(j);
        h3_p_p->Fill(p_p3);
        h3_xF->Fill(xF3);
        h3_p_theta->Fill(p_theta3 * 180.0 / 3.14159);
        h3_Mx->Fill(Mx3); // Always fill Mx histogram
    }

    for (Long64_t j = 0; j < nentries4; ++j) {
        tree4->GetEntry(j);
        h4_p_p->Fill(p_p4);
        h4_xF->Fill(xF4);
        h4_p_theta->Fill(p_theta4 * 180.0 / 3.14159);
        h4_Mx->Fill(Mx4); // Always fill Mx histogram
    }

    gStyle->SetOptStat(0);

    // Function to create TGraphErrors from ratio histograms
    auto createRatioGraph = [](TH1F* h1, TH1F* h2) {
        int n = h1->GetNbinsX();
        std::vector<double> x, y, ex, ey;
        for (int i = 1; i <= n; ++i) {
            double bin1 = h1->GetBinContent(i);
            double bin2 = h2->GetBinContent(i);
            double error1 = h1->GetBinError(i);
            double error2 = h2->GetBinError(i);
            if (bin2 != 0) {
                double ratio = bin1 / bin2;
                double error = ratio * std::sqrt((error1 / bin1) * (error1 / bin1) + (error2 / bin2) * (error2 / bin2));
                x.push_back(h1->GetBinCenter(i));
                y.push_back(ratio);
                ex.push_back(0);
                ey.push_back(error);
            }
        }
        return new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
    };

    // Create TGraphErrors for each ratio
    TGraphErrors* graph_p_p_1 = createRatioGraph(h1_p_p, h2_p_p);
    TGraphErrors* graph_p_p_2 = createRatioGraph(h3_p_p, h4_p_p);
    TGraphErrors* graph_xF_1 = createRatioGraph(h1_xF, h2_xF);
    TGraphErrors* graph_xF_2 = createRatioGraph(h3_xF, h4_xF);
    TGraphErrors* graph_p_theta_1 = createRatioGraph(h1_p_theta, h2_p_theta);
    TGraphErrors* graph_p_theta_2 = createRatioGraph(h3_p_theta, h4_p_theta);
    TGraphErrors* graph_Mx_1 = createRatioGraph(h1_Mx, h2_Mx);
    TGraphErrors* graph_Mx_2 = createRatioGraph(h3_Mx, h4_Mx);

    // Function to plot TGraphErrors
    auto plotGraph = [](TGraphErrors* g1, TGraphErrors* g2, const char* xTitle, const char* yTitle, const char* filename, double yMin = 1.0, double yMax = 3.0) {
        TCanvas* c = new TCanvas("c", "c", 800, 600);
        g1->SetLineColor(kRed);
        g1->SetMarkerColor(kRed);
        g1->SetMarkerStyle(20);
        g1->GetXaxis()->SetTitle(xTitle);
        g1->GetYaxis()->SetTitle(yTitle);
        g1->GetYaxis()->SetRangeUser(yMin, yMax);
        g1->Draw("AP");

        g2->SetLineColor(kBlue);
        g2->SetMarkerColor(kBlue);
        g2->SetMarkerStyle(21);
        g2->Draw("P SAME");

        TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
        legend->AddEntry(g1, "preliminary", "lp");
        legend->AddEntry(g2, "pass-1", "lp");
        legend->Draw();

        c->SaveAs(filename);
        delete c;
    };

    // Plot the graphs
    plotGraph(graph_p_p_1, graph_p_p_2, "p_p", "Ratio (NH3/C)", "output/ratio_p_p.png", 1.4, 2.0);
    plotGraph(graph_xF_1, graph_xF_2, "xF", "Ratio (NH3/C)", "output/ratio_xF.png");
    plotGraph(graph_p_theta_1, graph_p_theta_2, "p_theta (degrees) (NH3/C)", "Ratio", "output/ratio_p_theta.png");
    plotGraph(graph_Mx_1, graph_Mx_2, "Mx", "Ratio (NH3/C) ", "output/ratio_Mx.png");

    // Clean up
    delete graph_p_p_1;
    delete graph_p_p_2;
    delete graph_xF_1;
    delete graph_xF_2;
    delete graph_p_theta_1;
    delete graph_p_theta_2;
    delete graph_Mx_1;
    delete graph_Mx_2;

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
    delete h1_Mx;
    delete h2_Mx;
    delete h3_Mx;
    delete h4_Mx;

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