#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip> // For std::setprecision
#include <sstream> // For std::stringstream
#include <cmath> // For M_PI

using namespace std;

void analyzePions() {
    // Open the ROOT file
    TFile* file = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/mc/epi+X/rga_fa18_inb_clasdis_50nA_epi+X.root");
    if (!file || file->IsZombie()) {
        cerr << "Error opening file" << endl;
        return;
    }

    // Access the tree
    TTree* tree = dynamic_cast<TTree*>(file->Get("PhysicsEvents"));
    if (!tree) {
        cerr << "Error: Tree 'PhysicsEvents' not found" << endl;
        file->Close();
        return;
    }

    // Define variables to hold branch data
    Double_t p_theta, Q2, W, y, xF, Mx, z, p_p, pT, phi;
    Int_t mc_p1_parent;

    // Set branch addresses
    tree->SetBranchAddress("p_theta", &p_theta);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("xF", &xF);
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("mc_p1_parent", &mc_p1_parent);

    // Define pT bins and create a histogram for theta distributions
    vector<double> pTBins = {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00};
    vector<TH1F*> thetaHistograms;

    for (size_t i = 0; i < pTBins.size() - 1; ++i) {
        stringstream histName, histTitle;
        histName << "thetaHist_" << i;
        histTitle << fixed << setprecision(2) << pTBins[i] << " < P_T (GeV) < " << pTBins[i+1] << ";#theta (degrees);Counts";
        thetaHistograms.push_back(new TH1F(histName.str().c_str(), "", 80, 0, 40)); // 80 bins for better resolution
    }

    // Loop through the tree to fill histograms
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        double thetaDeg = p_theta * 180.0 / M_PI; // Convert to degrees

        if (Q2 > 1 && W > 2 && y < 0.75 && xF > 0 && Mx > 1.5 && z > 0.2 && p_p > 1.2 && y > 0.65 && y < 0.75 && z > 0.24 && z < 0.29 && Q2 > 2.0 && Q2 < 2.5) {
            for (size_t bin = 0; bin < pTBins.size() - 1; ++bin) {
                if (pT > pTBins[bin] && pT <= pTBins[bin + 1]) {
                    thetaHistograms[bin]->Fill(thetaDeg);
                }
            }
        }
    }

    // Create a new canvas for theta distribution
    TCanvas* thetaCanvas = new TCanvas("thetaCanvas", "Combined Theta Distributions", 800, 600);
    thetaCanvas->cd();

    // Create legend
    TLegend* legend = new TLegend(0.1, 0.7, 0.3, 0.9);
    legend->SetTextSize(0.03);

    for (size_t i = 0; i < thetaHistograms.size(); ++i) {
        thetaHistograms[i]->SetStats(0); // Explicitly hide the stats box for each histogram
    }

    // Custom colors for visibility
    int colors[] = {kBlack, kBlue, kRed, kGreen+2, kMagenta, kCyan, kOrange, kGray};

    for (size_t i = 0; i < thetaHistograms.size(); ++i) {
        thetaHistograms[i]->SetLineColor(colors[i % 8]);
        thetaHistograms[i]->SetLineWidth(2);
        stringstream label;
        label << fixed << setprecision(2) << pTBins[i] << " < P_T < " << pTBins[i+1] << " GeV";
        legend->AddEntry(thetaHistograms[i], label.str().c_str(), "l");

        if (i == 0) thetaHistograms[i]->Draw("HIST");
        else thetaHistograms[i]->Draw("HIST SAME");
    }

    legend->Draw();

    // Save the combined theta distribution canvas to a file
    thetaCanvas->SaveAs("output/combined_theta_distribution.png");

    // Cleanup
    file->Close();
}

int main() {
    analyzePions();
    return 0;
}
