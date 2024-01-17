#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

void vertex_study() {
    // Open the ROOT file
    TFile *file = TFile::Open("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    if (file == 0) {
        // If the file couldn't be opened, print an error and return
        printf("Error: could not open the file.\n");
        return;
    }

    // Get the TTree
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (tree == 0) {
        // If the tree couldn't be found, print an error and return
        printf("Error: could not find the tree.\n");
        return;
    }

    // Create the histogram
    TH1D *h1 = new TH1D("h1", "v_{z} distribution;v_{z};counts", 100, -15, 5);
    h1->SetLineColor(kBlack);

    // Draw the histogram
    tree->Draw("vz_e>>h1");

    // Customize the histogram appearance
    h1->GetXaxis()->CenterTitle(true);
    h1->GetYaxis()->CenterTitle(true);
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);

    // Create a canvas to draw the histogram
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    h1->Draw();

    // Remove the stat box
    h1->SetStats(0);

    // Add a legend
    TLegend *leg = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg->AddEntry(h1, "e^{-}", "l");
    leg->Draw();

    // Save the histogram as an image
    c1->SaveAs("vz_distribution.png");

    // Close the file
    file->Close();
}

