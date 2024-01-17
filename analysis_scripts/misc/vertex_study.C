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
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return;
    } else {
        std::cout << "Opened file: " << file_path << std::endl;
    }

    TTree* tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Tree PhysicsEvents not found in file: " << file_path << std::endl;
        file->Close();
        return;
    } else {
        std::cout << "Found tree PhysicsEvents in file." << std::endl;
    }

    TH1F* hist = new TH1F("hist_eX", "fa18_inb;v_{z};Counts", 100, -15, 10);
    tree->Draw("vz_e>>hist_eX");
    
    int nEntries = hist->GetEntries();
    std::cout << "Number of entries in histogram: " << nEntries << std::endl;

    if (nEntries > 0) {
        hist->SetLineColor(kBlack);
        hist->Draw();
    } else {
        std::cout << "Histogram is empty!" << std::endl;
    }

    file->Close(); // Close the file after use

    c1->SaveAs("output/single_histogram.png");
}
