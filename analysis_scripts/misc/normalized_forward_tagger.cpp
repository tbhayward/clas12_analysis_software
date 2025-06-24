/*
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
 *   2) Applies fiducial cuts to forward tagger hits.
 *   3) Loops over data and MC, filling photon (PID=22) FT hit-position histograms.
 *   4) Makes two unnormalized 2D plots (data vs. MC).
 *   5) Normalizes each histogram to unit integral.
 *   6) Computes and plots the ratio (data/MC) as "ft_ratio.png" with mean/stddev in a legend box.
 *   7) Creates a second plot "ft_ratio_outliers.png" showing a four-color map:
 *        • 0.5 ≤ ratio ≤ 2       → blue
 *        • 1/3 ≤ ratio < 0.5 and 2 < ratio ≤ 3   → light red (kOrange)
 *        • 1/5 ≤ ratio < 1/3 and 3 < ratio ≤ 5   → dark red (kRed)
 *        • ratio < 1/5 or ratio > 5  → pink (kPink)
 */

#include <iostream>
#include <string>
#include <cmath>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TLegend.h"

// ----------------------------------------------------------------------------
// Forward Tagger fiducial cut
// ----------------------------------------------------------------------------
bool forward_tagger_fiducial_cut(double x, double y) {
    double r = std::hypot(x, y);
    // radial acceptance
    if (r < 8.5 || r > 15.5) return false;
    // holes defined by {radius, center_x, center_y}
    const double holes[4][3] = {
        {1.60, -8.42,  9.89},
        {1.60, -9.89, -5.33},
        {2.30, -6.15,-13.00},
        {2.00,  3.70, -6.50}
    };
    for (auto& h : holes) {
        double hr = h[0], cx = h[1], cy = h[2];
        if (std::hypot(x - cx, y - cy) < hr) return false;
    }
    return true;
}

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
    // Create output directories if needed
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/ft", kTRUE);

    // 1) Parse arguments
    Long64_t maxEvents = -1;
    if (argc > 1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents == 0) maxEvents = -1;
    }
    bool useDataFile = (argc > 2);
    bool useMCFile   = (argc > 3);
    std::string dataFile = useDataFile ? argv[2] : "";
    std::string mcFile   = useMCFile   ? argv[3] : "";

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
    mcCh .SetBranchAddress("particle_pid",  &pid);
    mcCh .SetBranchAddress("ft_x",          &ft_x);
    mcCh .SetBranchAddress("ft_y",          &ft_y);
    mcCh .SetBranchAddress("ft_energy",     &ft_energy);

    // 4) Book histograms
    const int NB = 200;
    TH2D* h_data = new TH2D("h_ft_data",
      "FT Hit Position (data); x_{FT} (cm); y_{FT} (cm)",
      NB, -20, 20, NB, -20, 20);
    TH2D* h_mc   = new TH2D("h_ft_mc",
      "FT Hit Position (mc);   x_{FT} (cm); y_{FT} (cm)",
      NB, -20, 20, NB, -20, 20);

    gStyle->SetOptStat(0);

    // 5) Fill DATA with fiducial cut
    Long64_t nD = dataCh.GetEntries();
    if (maxEvents > 0 && maxEvents < nD) nD = maxEvents;
    for (Long64_t i = 0; i < nD; ++i) {
        dataCh.GetEntry(i);
        if (pid != 22 || ft_energy == -9999) continue;
        if (!forward_tagger_fiducial_cut(ft_x, ft_y)) continue;
        h_data->Fill(ft_x, ft_y);
    }

    // 6) Fill MC with fiducial cut
    Long64_t nM = mcCh.GetEntries();
    if (maxEvents > 0 && maxEvents < nM) nM = maxEvents;
    for (Long64_t i = 0; i < nM; ++i) {
        mcCh.GetEntry(i);
        if (pid != 22 || ft_energy == -9999) continue;
        if (!forward_tagger_fiducial_cut(ft_x, ft_y)) continue;
        h_mc->Fill(ft_x, ft_y);
    }

    // 7) Draw unnormalized plots
    {
        TCanvas c("c_ft_data","FT Data Unnormalized",600,600);
        c.cd(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15);
        h_data->Draw("COLZ");
        c.SaveAs("output/ft/ft_data.png");
    }
    {
        TCanvas c("c_ft_mc","FT MC Unnormalized",600,600);
        c.cd(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15);
        h_mc->Draw("COLZ");
        c.SaveAs("output/ft/ft_mc.png");
    }

    // 8) Normalize histograms
    double Idata = h_data->Integral(); if (Idata > 0) h_data->Scale(1.0/Idata);
    double Imc   = h_mc->Integral();   if (Imc   > 0) h_mc->Scale(1.0/Imc);

    // 9) Compute ratio
    TH2D* h_ratio = (TH2D*)h_data->Clone("h_ft_ratio");
    h_ratio->SetTitle("FT Hit Position Data/MC; x_{FT} (cm); y_{FT} (cm)");
    h_ratio->Divide(h_mc);

    // compute mean and stddev of non-zero bins
    int cnt = 0; double sum = 0, sum2 = 0;
    for (int ix = 1; ix <= h_ratio->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= h_ratio->GetNbinsY(); ++iy) {
            double v = h_ratio->GetBinContent(ix, iy);
            if (v == 0) continue;
            sum  += v; 
            sum2 += v*v; 
            ++cnt;
        }
    }
    double mean  = cnt ? sum/cnt : 0;
    double sigma = cnt ? std::sqrt(sum2/cnt - mean*mean) : 0;

    // 10) Draw ratio with legend box
    SetSame2DScale(h_data, h_mc, h_ratio);
    {
        TCanvas c("c_ft_ratio","FT Data/MC Ratio",600,600);
        c.cd(); gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15);
        h_ratio->Draw("COLZ");
        TLegend leg(0.6, 0.7, 0.9, 0.9);
        leg.SetFillColor(kWhite);
        leg.SetBorderSize(1);
        leg.SetTextSize(0.03);
        leg.AddEntry((TObject*)0, Form("Mean = %.4f",  mean), "");
        leg.AddEntry((TObject*)0, Form("StdDev = %.4f", sigma), "");
        leg.Draw();
        c.SaveAs("output/ft/ft_ratio.png");
    }

    // 11) Draw outliers map with four-color categories
    {
        TCanvas c("c_ft_outliers","FT Ratio Outliers",600,600);
        c.cd(); 
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.15);
        gPad->SetLogz(0);

        // clone & reset
        TH2D* h_map = (TH2D*)h_ratio->Clone("h_ft_map");
        h_map->Reset();

        int nX = h_ratio->GetNbinsX(), nY = h_ratio->GetNbinsY();
        for (int ix = 1; ix <= nX; ++ix) {
            for (int iy = 1; iy <= nY; ++iy) {
                double v = h_ratio->GetBinContent(ix, iy);
                if (v <= 0) continue;
                int cat = 4; // default: <1/5 or >5 → pink
                if (v >= 0.5 && v <= 2.0) {
                    cat = 1; // blue
                }
                else if ((v >= 1.0/3.0 && v < 0.5) ||
                         (v > 2.0   && v <= 3.0)) {
                    cat = 2; // light red (orange)
                }
                else if ((v >= 1.0/5.0 && v < 1.0/3.0) ||
                         (v > 3.0   && v <= 5.0)) {
                    cat = 3; // dark red
                }
                h_map->SetBinContent(ix, iy, cat);
            }
        }

        // set four-color palette
        Int_t palette[4] = { kBlue, kOrange, kRed, kPink };
        gStyle->SetPalette(4, palette);

        // discrete contours
        h_map->SetContour(4);
        h_map->SetContourLevel(0, 1);
        h_map->SetContourLevel(1, 2);
        h_map->SetContourLevel(2, 3);
        h_map->SetContourLevel(3, 4);

        h_map->Draw("COLZ");
        c.SaveAs("output/ft/ft_ratio_outliers.png");
    }

    return 0;
}