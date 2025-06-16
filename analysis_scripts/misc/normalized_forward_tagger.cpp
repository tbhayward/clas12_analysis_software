/*****************************************************************************
 * File: forward_tagger_comparison.cpp
 *
 * Compile with:
 *   g++ forward_tagger_comparison.cpp -o forward_tagger_comparison \
 *       `root-config --cflags --libs`
 *
 * Run with:
 *   ./forward_tagger_comparison [Nevents] [dataFile] [mcFile]
 *
 * What it does:
 *   1) Creates output/ and output/ft/ if they do not exist.
 *   2) Loops over data and MC, filling photon (PID=22) FT hit‐position histograms.
 *   3) Makes two unnormalized 2D plots (data vs. MC).
 *   4) Normalizes each histogram to unit integral.
 *   5) Computes and plots the ratio (data/MC) in a third canvas.
 *****************************************************************************/

#include <iostream>
#include <string>
#include <cmath>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TString.h"

// ----------------------------------------------------------------------------
// Helper to unify the color scale of three TH2D histograms
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* a, TH2D* b, TH2D* c) {
    double mn =  1e9, mx = -1e9;
    for (auto* h : {a, b, c}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto* h : {a, b, c}) {
        h->SetMinimum(mn);
        h->SetMaximum(mx);
    }
}

int main(int argc, char** argv) {
    // --- Create output directories if needed ---
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/ft", kTRUE);

    // 1) Parse arguments
    Long64_t maxEvents = -1;
    if (argc > 1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents == 0) maxEvents = -1;
    }
    bool useDataFile = (argc > 2);
    std::string dataFile = useDataFile ? argv[2] : "";
    bool useMCFile = (argc > 3);
    std::string mcFile   = useMCFile   ? argv[3]   : "";

    // 2) Set up TChains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if (useDataFile) dataCh.Add(dataFile.c_str());
    else             dataCh.Add("/work/clas12/thayward/.../ft_data.root");
    if (useMCFile)   mcCh.Add(mcFile.c_str());
    else             mcCh.Add("/work/clas12/thayward/.../ft_mc.root");

    // 3) Branch addresses
    Int_t    pid;
    Double_t ft_x, ft_y, ft_energy;
    dataCh.SetBranchAddress("particle_pid", &pid);
    dataCh.SetBranchAddress("ft_x",         &ft_x);
    dataCh.SetBranchAddress("ft_y",         &ft_y);
    dataCh.SetBranchAddress("ft_energy",    &ft_energy);
    mcCh.SetBranchAddress("particle_pid",   &pid);
    mcCh.SetBranchAddress("ft_x",           &ft_x);
    mcCh.SetBranchAddress("ft_y",           &ft_y);
    mcCh.SetBranchAddress("ft_energy",      &ft_energy);

    // 4) Book histograms
    const int NB = 200;
    TH2D* h_data = new TH2D("h_ft_data",
      "FT Hit Position (data); x_{FT} (cm); y_{FT} (cm)",
      NB, -20, 20, NB, -20, 20);
    TH2D* h_mc   = new TH2D("h_ft_mc",
      "FT Hit Position (mc);   x_{FT} (cm); y_{FT} (cm)",
      NB, -20, 20, NB, -20, 20);

    gStyle->SetOptStat(0);

    // 5) Fill DATA
    Long64_t nD = dataCh.GetEntries();
    if (maxEvents > 0 && maxEvents < nD) nD = maxEvents;
    for (Long64_t i = 0; i < nD; ++i) {
        dataCh.GetEntry(i);
        if (pid != 22)          continue;
        if (ft_energy == -9999) continue;
        h_data->Fill(ft_x, ft_y);
    }

    // 6) Fill MC
    Long64_t nM = mcCh.GetEntries();
    if (maxEvents > 0 && maxEvents < nM) nM = maxEvents;
    for (Long64_t i = 0; i < nM; ++i) {
        mcCh.GetEntry(i);
        if (pid != 22)          continue;
        if (ft_energy == -9999) continue;
        h_mc->Fill(ft_x, ft_y);
    }

    // 7) Draw unnormalized plots
    {
        TCanvas c("c_ft_data","FT Data Unnormalized",600,600);
        c.cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        h_data->Draw("COLZ");
        c.SaveAs("output/ft/ft_data.png");
    }
    {
        TCanvas c("c_ft_mc","FT MC Unnormalized",600,600);
        c.cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        h_mc->Draw("COLZ");
        c.SaveAs("output/ft/ft_mc.png");
    }

    // 8) Normalize histograms
    double Idata = h_data->Integral();
    if (Idata > 0) h_data->Scale(1.0/Idata);
    double Imc   = h_mc->Integral();
    if (Imc   > 0) h_mc->Scale(1.0/Imc);

    // 9) Compute ratio
    TH2D* h_ratio = (TH2D*)h_data->Clone("h_ft_ratio");
    h_ratio->SetTitle("FT Hit Position Data/MC; x_{FT} (cm); y_{FT} (cm)");
    h_ratio->Divide(h_mc);

    // 10) Draw ratio with unified color scale
    SetSame2DScale(h_data, h_mc, h_ratio);
    {
        TCanvas c("c_ft_ratio","FT Data/MC Ratio",600,600);
        c.cd();
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        h_ratio->Draw("COLZ");
        c.SaveAs("output/ft/ft_ratio.png");
    }

    return 0;
}