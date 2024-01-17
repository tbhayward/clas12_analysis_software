#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>

void vertex_study() {
    gStyle->SetOptStat(0); // Turn off the statistics box

    TCanvas* c1 = new TCanvas("c1", "Vertex Position Study", 800, 600);

    TString file_path = "/volatile/clas12/thayward/vertex_studies/rga/fa18_inb/rga_fa18_inb_eX.root";
    TFile* file = new TFile(file_path);

    TTree* tree = (TTree*)file->Get("PhysicsEvents");

    TH1F* hist = new TH1F("hist_eX", "fa18_inb;v_{z};Counts", 100, -15, 10);
    tree->Draw("vz_e>>hist_eX");

    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->Draw();

    file->Close(); // Close the file after use

    c1->SaveAs("output/single_histogram.png");
}
