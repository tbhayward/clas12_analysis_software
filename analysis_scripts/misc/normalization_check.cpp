#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include <map>
#include <iostream>

struct RunInfo {
    double total_charge;
    double pos_beam_charge;
    double neg_beam_charge;
    double target_pol;
    double target_pol_unc;
};

std::map<int, RunInfo> parseCSV(const std::string& filename) {
    std::map<int, RunInfo> runInfoMap;
    std::ifstream infile(filename);
    std::string line;

    while (std::getline(infile, line)) {
        if (line[0] == '#') continue; // Skip comment lines
        std::istringstream ss(line);
        int runnum;
        RunInfo info;
        char comma;
        ss >> runnum >> comma >> info.total_charge >> comma >> info.pos_beam_charge >> comma >> info.neg_beam_charge >> comma >> info.target_pol >> comma >> info.target_pol_unc;
        runInfoMap[runnum] = info;
    }
    return runInfoMap;
}

void normalization_check(const char* inputFileName) {
    // Open the input ROOT file and get the tree
    TFile* inputFile = TFile::Open(inputFileName);
    TTree* tree = (TTree*)inputFile->Get("PhysicsEvents");

    // Import the CSV file
    std::map<int, RunInfo> runInfoMap = parseCSV("/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv");

    // Set up branches
    int runnum, helicity;
    double Delta_phi;
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("helicity", &helicity);
    tree->SetBranchAddress("Delta_phi", &Delta_phi);

    // Create histograms
    TH1F* h1 = new TH1F("h1", "Helicity +, Target Pol +;#Delta_{#phi};Counts/nC", 200, 0, 2 * TMath::Pi());
    TH1F* h2 = new TH1F("h2", "Helicity +, Target Pol -;#Delta_{#phi};Counts/nC", 200, 0, 2 * TMath::Pi());
    TH1F* h3 = new TH1F("h3", "Helicity -, Target Pol +;#Delta_{#phi};Counts/nC", 200, 0, 2 * TMath::Pi());
    TH1F* h4 = new TH1F("h4", "Helicity -, Target Pol -;#Delta_{#phi};Counts/nC", 200, 0, 2 * TMath::Pi());

    // Accumulate charges for normalization
    double charge_h1 = 0.0, charge_h2 = 0.0, charge_h3 = 0.0, charge_h4 = 0.0;

    // Loop over the tree
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (helicity == 0) continue; // Skip if helicity is 0

        double target_pol = runInfoMap[runnum].target_pol;
        double pos_beam_charge = runInfoMap[runnum].pos_beam_charge;
        double neg_beam_charge = runInfoMap[runnum].neg_beam_charge;

        if (helicity > 0) {
            if (target_pol > 0) {
                h1->Fill(Delta_phi);
                charge_h1 += pos_beam_charge;
            } else {
                h2->Fill(Delta_phi);
                charge_h2 += pos_beam_charge;
            }
        } else {
            if (target_pol > 0) {
                h3->Fill(Delta_phi);
                charge_h3 += neg_beam_charge;
            } else {
                h4->Fill(Delta_phi);
                charge_h4 += neg_beam_charge;
            }
        }
    }

    // Normalize histograms
    if (charge_h1 > 0) h1->Scale(1.0 / charge_h1);
    if (charge_h2 > 0) h2->Scale(1.0 / charge_h2);
    if (charge_h3 > 0) h3->Scale(1.0 / charge_h3);
    if (charge_h4 > 0) h4->Scale(1.0 / charge_h4);

    // Draw histograms
    TCanvas* c = new TCanvas("c", "Normalization Check", 800, 600);
    gStyle->SetOptStat(0);
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);
    h3->SetLineColor(kGreen);
    h4->SetLineColor(kMagenta);

    h1->Draw("HIST");
    h2->Draw("HIST SAME");
    h3->Draw("HIST SAME");
    h4->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.75, 0.75, 0.9, 0.9);
    legend->AddEntry(h1, "Helicity +, Target Pol +", "l");
    legend->AddEntry(h2, "Helicity +, Target Pol -", "l");
    legend->AddEntry(h3, "Helicity -, Target Pol +", "l");
    legend->AddEntry(h4, "Helicity -, Target Pol -", "l");
    legend->Draw();

    c->SaveAs("normalization_check.pdf");

    // Clean up
    delete h1;
    delete h2;
    delete h3;
    delete h4;
    delete legend;
    delete c;
    inputFile->Close();
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input ROOT file>" << std::endl;
        return 1;
    }

    normalization_check(argv[1]);
    return 0;
}