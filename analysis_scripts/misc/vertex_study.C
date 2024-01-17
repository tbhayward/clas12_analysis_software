// Include necessary ROOT headers
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>

void vertex_study() {
    // Import the file and get the tree
    TFile *file = new TFile("/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root");
    TTree *tree = (TTree*)file->Get("PhysicsEvents");

    // Create a histogram for vz_e
    TH1F *h1 = new TH1F("h1", "Vertex Study", 100, -15, 5);
    h1->SetLineColor(kBlack);

    // Fill the histogram from the tree
    tree->Draw("vz_e>>h1");

    // Style the histogram
    h1->GetXaxis()->SetTitle("v_{z}");
    h1->GetYaxis()->SetTitle("counts");
    h1->GetXaxis()->CenterTitle();
    h1->GetYaxis()->CenterTitle();
    h1->GetXaxis()->SetTitleSize(0.05);
    h1->GetYaxis()->SetTitleSize(0.05);

    // Add a legend
    TLegend *legend = new TLegend(0.7, 0.9, 0.9, 0.95);
    legend->AddEntry(h1, "e^{-}", "l");
    legend->Draw();

    // Remove the stat box
    gStyle->SetOptStat(0);

    // Draw the histogram
    h1->Draw("HIST");

    // Save the canvas as a PNG file
    c1->SaveAs("output/neg_vertices.png");
}

int main() {
    vertex_study();
    return 0;
}
