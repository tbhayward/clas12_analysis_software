#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TAxis.h>
#include <iostream>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char* argv[])
{
    if (argc != 4) {
        printf("Usage: %s file1.root file2.root file3.root\n", argv[0]);
        return 1;
    }

    // Set style to make the plot look nice and scientific
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetLabelSize(0.04, "xy");
    gStyle->SetTitleSize(0.05, "xy");

    // Open the ROOT files containing the PhysicsEvents trees
    TFile *file1 = TFile::Open(argv[1]);
    TFile *file2 = TFile::Open(argv[2]);
    TFile *file3 = TFile::Open(argv[3]);

    if (!file1 || !file2 || !file3) {
        printf("Error opening one of the files.\n");
        return 1;
    }

    // Retrieve the PhysicsEvents trees from each file
    TTree *tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents");
    TTree *tree3 = (TTree*)file3->Get("PhysicsEvents");

    if (!tree1 || !tree2 || !tree3) {
        printf("Error retrieving PhysicsEvents tree from one of the files.\n");
        return 1;
    }

    // Create histograms for Mx2 from each tree with the new range
    TH1F *h1 = new TH1F("h1", "", 150, -4, 8);
    TH1F *h2 = new TH1F("h2", "", 150, -4, 8);
    TH1F *h3 = new TH1F("h3", "", 150, -4, 8);

    // Fill the histograms with Mx2 data
    tree1->Draw("Mx2>>h1", "", "goff");
    tree2->Draw("Mx2>>h2", "", "goff");
    tree3->Draw("Mx2>>h3", "", "goff");

    // Normalize the histograms to the specified values
    double norm1 = 0.3725622171315408;
    double norm2 = 0.050721607409;
    double norm3 = 0.3725622171315408;
    // double norm1 = 1;
    // double norm2 = 1;
    // double norm3 = 1;

    h1->Scale(1.0 / norm1);
    h2->Scale(1.0 / norm2);
    h3->Scale(1.0 / norm3);

    // Set styles for the histograms
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h2->SetLineColor(kBlue);
    h2->SetLineWidth(2);
    h3->SetLineColor(kGreen+2);
    h3->SetLineWidth(2);

    // Set axis labels
    h1->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h1->GetYaxis()->SetTitle("Counts/nC");

    // Create a canvas to draw the histograms
    TCanvas *c = new TCanvas("c", "", 800, 600);

    // Set canvas margins for additional padding
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);

    // Adjust the maximum value for better display
    double max1 = h1->GetMaximum();
    double max2 = h2->GetMaximum();
    double max3 = h3->GetMaximum();
    double max = std::max({max1, max2, max3});
    h1->SetMaximum(max * 1.1);

    // Draw the histograms on the same canvas
    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");

    // Add a legend with the specified labels
    TLegend *leg = new TLegend(0.45, 0.7, 0.9, 0.9);
    leg->AddEntry(h1, "H_{2}", "l");
    leg->AddEntry(h2, "C", "l");
    leg->AddEntry(h3, "C, p_{z} = Gauss(0,0.033)", "l");
    leg->SetTextSize(0.04);
    leg->Draw();

    // Create the output directory if it doesn't exist
    struct stat info;
    if (stat("output", &info) != 0) {
        // Directory does not exist; create it
        if (mkdir("output", 0777) != 0) {
            perror("Error creating output directory");
            return 1;
        }
    } else if (!(info.st_mode & S_IFDIR)) {
        // 'output' exists but is not a directory
        fprintf(stderr, "Error: 'output' exists but is not a directory\n");
        return 1;
    }

    // Save the plot to the specified file
    c->SaveAs("output/fermi_motion_Mx2.png");

    // Close the files and clean up
    file1->Close();
    file2->Close();
    file3->Close();

    delete h1;
    delete h2;
    delete h3;
    delete c;
    delete leg;

    return 0;
}