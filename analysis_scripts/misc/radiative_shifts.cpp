#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <iostream>

int main(int argc, char** argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <branch_variable> <x_axis_label>" << std::endl;
        return 1;
    }

    const char* file1_name = argv[1];
    const char* file2_name = argv[2];
    const char* branch_name = argv[3];
    const char* x_axis_label = argv[4];

    // Open the input files
    TFile* file1 = TFile::Open(file1_name);
    TFile* file2 = TFile::Open(file2_name);
    if (!file1 || !file2 || file1->IsZombie() || file2->IsZombie()) {
        std::cerr << "Error: Could not open one or both ROOT files." << std::endl;
        return 1;
    }

    // Load the trees
    TTree* tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree* tree2 = (TTree*)file2->Get("PhysicsEvents");
    if (!tree1 || !tree2) {
        std::cerr << "Error: Could not find tree 'PhysicsEvents' in one or both files." << std::endl;
        return 1;
    }

    // Define histograms
    TH1D* hist1 = new TH1D("hist1", branch_name, 100, -10, 10);
    TH1D* hist2 = new TH1D("hist2", branch_name, 100, -10, 10);

    // Project the branch variable into the histograms
    tree1->Project("hist1", branch_name);
    tree2->Project("hist2", branch_name);

    // Customize histogram appearance
    hist1->SetLineColor(kBlue);
    hist2->SetLineColor(kRed);
    hist1->GetXaxis()->SetTitle(Form("#it{%s}", x_axis_label)); // Set LaTeX-style title
    hist1->GetYaxis()->SetTitle("Counts");
    hist2->GetXaxis()->SetTitle(Form("#it{%s}", x_axis_label));
    hist2->GetYaxis()->SetTitle("Counts");

    // Remove the stat boxes
    hist1->SetStats(0);
    hist2->SetStats(0);

    // Create a canvas and draw the histograms
    TCanvas* canvas = new TCanvas("canvas", "Comparison of Branch Variable", 800, 600);
    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, "File 1", "l");
    legend->AddEntry(hist2, "File 2", "l");
    legend->Draw();

    // Save the plot
    canvas->SaveAs("output/rad_study/test.pdf");

    // Clean up
    delete canvas;
    delete hist1;
    delete hist2;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    return 0;
}