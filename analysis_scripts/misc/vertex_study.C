// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>

void vertex_study() {
    // Import the file and get the tree

    // RGA Fa18 Inb
    TFile *rgafa18inbeXfile = new TFile("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    TTree *rgafa18inbeXtree = (TTree*)rgafa18inbeXfile->Get("PhysicsEvents");
    TFile *rgafa18inbepimXfile = new TFile("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_epi-X.root");
    TTree *rgafa18inbepimXtree = (TTree*)rgafa18inbepimXfile->Get("PhysicsEvents");
    TFile *rgafa18inbekmXfile = new TFile("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_ek-X.root");
    TTree *rgafa18inbekmXtree = (TTree*)rgafa18inbekmXfile->Get("PhysicsEvents");

    // Create a histogram for vz_e
    TH1F *h_rgafa18inbeX = new TH1F("h_rgafa18inbeX", "RGA Fa18 Inb", 100, -15, 10);
    h_rgafa18inbeX->SetLineColor(kBlack);

    // Fill the histogram from the tree
    rgafa18inbeXtree->Draw("vz_e>>h_rgafa18inbeX");

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetBottomMargin(0.15);

    // Style the histogram
    h_rgafa18inbeX->GetXaxis()->SetTitle("v_{z} (cm)");
    h_rgafa18inbeX->GetYaxis()->SetTitle("counts");
    h_rgafa18inbeX->GetXaxis()->CenterTitle();
    h_rgafa18inbeX->GetYaxis()->CenterTitle();
    h_rgafa18inbeX->GetXaxis()->SetTitleSize(0.05);
    h_rgafa18inbeX->GetYaxis()->SetTitleSize(0.05);

    // Draw the histogram on the canvas
    h_rgafa18inbeX->Draw("HIST");

    // Add a legend
    TLegend *legend = new TLegend(0.7 - c1->GetLeftMargin(), 0.9, 0.9 - c1->GetLeftMargin(), 0.95);
    legend->AddEntry(h_rgafa18inbeX, "e^{-}", "l");
    legend->Draw();

    // Remove the stat box
    gStyle->SetOptStat(0);

    // Save the canvas as a PNG file
    c1->SaveAs("output/neg_vertices.png");
}

int main() {
    vertex_study();
    return 0;
}
