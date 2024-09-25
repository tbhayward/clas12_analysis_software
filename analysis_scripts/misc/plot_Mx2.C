#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TF1.h>

void plot_Mx2() {
    // Open the ROOT file
    TFile *file = TFile::Open("/Users/tbhayward/Desktop/updated_script_test/rga_fa18_inb_epX.root");
    if (!file || file->IsZombie()) {
        printf("Error: Could not open file.\n");
        return;
    }

    // Load the tree
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        printf("Error: Could not load tree.\n");
        file->Close();
        return;
    }

    // Create a histogram for Mx2 with 200 bins between -1 and 1
    TH1F *h_Mx2 = new TH1F("h_Mx2", "Mx^{2} distribution", 200, -0.5, 0.4);

    // Draw the Mx2 branch into the histogram
    tree->Draw("Mx2 >> h_Mx2");

    // Set axis labels
    h_Mx2->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
    h_Mx2->GetYaxis()->SetTitle("Counts");

    // Set histogram style to data points with vertical error bars
    h_Mx2->SetMarkerStyle(20); // Circle marker
    h_Mx2->SetMarkerSize(0.7);   // Slightly smaller marker size
    h_Mx2->SetMarkerColor(kBlack); // Marker color black
    h_Mx2->SetLineColor(kBlack);   // Line color black

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "Mx^{2} Plot", 800, 600);

    // Remove the stat box
    gStyle->SetOptStat(0);

    // Set additional spacing on the left to avoid cutting off labels
    c1->SetLeftMargin(0.15);

    // Draw the histogram with error bars (E1 option)
    h_Mx2->Draw("E1");

    // Define a fitting function: quadratic + two gaussians
    TF1 *fitFunc = new TF1("fitFunc", "[0] + [1]*x + [2]*x*x + [3]*x*x + [4]*exp(-0.5*((x-[5])/[6])**2))", -1, 1);

    // Set initial parameters for the fit
    fitFunc->SetParameters(100, 0, 0, 0, 50, 0, 0.1); // Initial guesses for params

    // Set parameter names for clarity
    fitFunc->SetParNames("const", "linear", "quad", "cubic", "Gaussian1_amp", "Gaussian1_mean", "Gaussian1_sigma");

    // Perform the fit
    h_Mx2->Fit("fitFunc");

    // Save the plot
    c1->SaveAs("/Users/tbhayward/Desktop/rga_Mx2_plot_with_fit.png");

    // Clean up
    file->Close();
    delete c1;
    delete h_Mx2;
}