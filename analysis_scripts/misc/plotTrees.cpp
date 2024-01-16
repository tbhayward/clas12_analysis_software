#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>

void plotTrees(const char* file1, const char* file2) {
    // Open the ROOT files and get the trees
    TFile *f1 = TFile::Open(file1);
    TFile *f2 = TFile::Open(file2);
    TTree *tree1 = (TTree*)f1->Get("PhysicsEvents");
    TTree *tree2 = (TTree*)f2->Get("PhysicsEvents");

    // Create histograms
    TH1F *hist1 = new TH1F("hist1", "M_{X} Distribution;M_{X} (GeV);Counts", 100, 0, 4);
    TH1F *hist2 = new TH1F("hist2", "M_{X} Distribution;M_{X} (GeV);Counts", 100, 0, 4);

    // Set histogram colors
    hist1->SetLineColor(kRed);
    hist2->SetLineColor(kBlue);

    // Draw from trees
    tree1->Draw("Mx>>hist1");
    tree2->Draw("Mx>>hist2");

    // Create a canvas to draw histograms
    TCanvas *c = new TCanvas("c", "M_{X} Distributions", 800, 600);
    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");

    // Create and add legend
    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg->AddEntry(hist1, "e#pi^{+}X", "l");
    leg->AddEntry(hist2, "e#pi^{-}X", "l");
    leg->Draw();

    // Save the canvas
    c->SaveAs("Mx_distribution.png");
}

int main(int argc, char **argv) {
    if (argc != 3) {
        printf("Usage: %s <ROOT file 1> <ROOT file 2>\n", argv[0]);
        return 1;
    }
    plotTrees(argv[1], argv[2]);
    return 0;
}
