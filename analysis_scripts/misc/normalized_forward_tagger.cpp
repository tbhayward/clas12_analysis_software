/*
 * File: calorimeter_comparison.cpp
 *
 * Compile with:
 *   g++ calorimeter_comparison.cpp -o calorimeter_comparison \
 *       `root-config --cflags --libs`
 *
 * Run with:
 *   ./calorimeter_comparison [Nevents] [dataFile] [mcFile]
 *
 * What it does:
 *   1) Creates output/ and output/cal/ if they do not exist.
 *   2) Loops over data and MC, filling photon (PID=22) and electron (PID=11)
 *      calorimeter hit-position histograms for PCal, ECin, and ECout.
 *   3) Makes unnormalized 2D plots (data vs. MC) with all three layers on one canvas.
 *   4) Normalizes each histogram and computes data/MC ratio.
 *   5) Creates a 3-panel ratio plot with mean (μ) and stddev (σ) annotated per layer.
 *   6) Creates a second 3-panel plot highlighting bins with ratio >2 or <0.5 in red.
 */

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <utility>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TBox.h"
#include "TString.h"

// ----------------------------------------------------------------------------
// Helper to unify color scale for three TH2D histograms
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* a, TH2D* b, TH2D* c) {
    double mn = 1e9, mx = -1e9;
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
    // 1) Create output directories
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/cal",kTRUE);

    // 2) Parse arguments
    Long64_t maxEvents = -1;
    if (argc > 1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents == 0) maxEvents = -1;
    }
    bool useData = (argc > 2);
    bool useMC   = (argc > 3);
    std::string dataFile = useData ? argv[2] : "";
    std::string mcFile   = useMC   ? argv[3] : "";

    // 3) Set up TChains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if (useData) dataCh.Add(dataFile.c_str());
    else         dataCh.Add("/work/clas12/thayward/.../cal_data.root");
    if (useMC)   mcCh.Add(mcFile.c_str());
    else         mcCh.Add("/work/clas12/thayward/.../cal_mc.root");

    // 4) Branch addresses
    Int_t    pid;
    Double_t cal_x_1, cal_y_1;
    Double_t cal_x_4, cal_y_4;
    Double_t cal_x_7, cal_y_7;
    dataCh.SetBranchAddress("particle_pid", &pid);
    dataCh.SetBranchAddress("cal_x_1",      &cal_x_1);
    dataCh.SetBranchAddress("cal_y_1",      &cal_y_1);
    dataCh.SetBranchAddress("cal_x_4",      &cal_x_4);
    dataCh.SetBranchAddress("cal_y_4",      &cal_y_4);
    dataCh.SetBranchAddress("cal_x_7",      &cal_x_7);
    dataCh.SetBranchAddress("cal_y_7",      &cal_y_7);
    mcCh.SetBranchAddress("particle_pid",    &pid);
    mcCh.SetBranchAddress("cal_x_1",         &cal_x_1);
    mcCh.SetBranchAddress("cal_y_1",         &cal_y_1);
    mcCh.SetBranchAddress("cal_x_4",         &cal_x_4);
    mcCh.SetBranchAddress("cal_y_4",         &cal_y_4);
    mcCh.SetBranchAddress("cal_x_7",         &cal_x_7);
    mcCh.SetBranchAddress("cal_y_7",         &cal_y_7);

    // 5) Configuration
    const int NB = 200;
    const double xmin = -450.0, xmax = 450.0;
    const double ymin = -450.0, ymax = 450.0;

    std::vector<std::pair<int,std::string>> species = {{22,"photon"},{11,"electron"}};
    std::vector<std::string> layerLabels = {"PCal","ECin","ECout"};

    for (auto& sp : species) {
        int pidVal = sp.first;
        std::string label = sp.second;

        // Histograms for data, MC, and ratio
        std::vector<TH2D*> hData(3), hMC(3), hRatio(3);
        for (int i = 0; i < 3; ++i) {
            hData[i] = new TH2D(
                Form("h_data_%s_%s", label.c_str(), layerLabels[i].c_str()),
                Form("%s Data %s; x (cm); y (cm)", label.c_str(), layerLabels[i].c_str()),
                NB, xmin, xmax, NB, ymin, ymax);
            hMC[i] = (TH2D*)hData[i]->Clone(
                Form("h_mc_%s_%s", label.c_str(), layerLabels[i].c_str()));
            hMC[i]->SetTitle(
                Form("%s MC %s", label.c_str(), layerLabels[i].c_str()));
        }

        // 6) Fill data
        Long64_t nD = dataCh.GetEntries();
        if (maxEvents > 0 && maxEvents < nD) nD = maxEvents;
        for (Long64_t ev = 0; ev < nD; ++ev) {
            dataCh.GetEntry(ev);
            if (pid != pidVal) continue;
            double xs[3] = {cal_x_1, cal_x_4, cal_x_7};
            double ys[3] = {cal_y_1, cal_y_4, cal_y_7};
            for (int i = 0; i < 3; ++i) {
                if (xs[i] == -9999 || ys[i] == -9999) continue;
                hData[i]->Fill(xs[i], ys[i]);
            }
        }
        // 7) Fill MC
        Long64_t nM = mcCh.GetEntries();
        if (maxEvents > 0 && maxEvents < nM) nM = maxEvents;
        for (Long64_t ev = 0; ev < nM; ++ev) {
            mcCh.GetEntry(ev);
            if (pid != pidVal) continue;
            double xs[3] = {cal_x_1, cal_x_4, cal_x_7};
            double ys[3] = {cal_y_1, cal_y_4, cal_y_7};
            for (int i = 0; i < 3; ++i) {
                if (xs[i] == -9999 || ys[i] == -9999) continue;
                hMC[i]->Fill(xs[i], ys[i]);
            }
        }

        gStyle->SetOptStat(0);

        // 8) Unnormalized plots
        {
            TCanvas c("c_data_uncut", "Data Uncut", 1800, 600);
            c.Divide(3,1);
            SetSame2DScale(hData[0], hData[1], hData[2]);
            for (int i = 0; i < 3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                hData[i]->Draw("COLZ");
            }
            c.SaveAs(Form("output/cal/data_uncut_%s.png", label.c_str()));
        }
        {
            TCanvas c("c_mc_uncut", "MC Uncut", 1800, 600);
            c.Divide(3,1);
            SetSame2DScale(hMC[0], hMC[1], hMC[2]);
            for (int i = 0; i < 3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                hMC[i]->Draw("COLZ");
            }
            c.SaveAs(Form("output/cal/mc_uncut_%s.png", label.c_str()));
        }

        // 9) Normalize
        for (int i = 0; i < 3; ++i) {
            double id = hData[i]->Integral(); if (id > 0) hData[i]->Scale(1.0/id);
            double im = hMC[i]->Integral();   if (im > 0) hMC[i]->Scale(1.0/im);
        }
        // 10) Compute ratio and stats
        std::vector<double> mean(3), sigma(3);
        for (int i = 0; i < 3; ++i) {
            hRatio[i] = (TH2D*)hData[i]->Clone(
                Form("h_ratio_%s_%s", label.c_str(), layerLabels[i].c_str()));
            hRatio[i]->SetTitle(
                Form("%s Data/MC %s", label.c_str(), layerLabels[i].c_str()));
            hRatio[i]->Divide(hMC[i]);
            int cnt = 0; double s = 0, s2 = 0;
            for (int ix = 1; ix <= hRatio[i]->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= hRatio[i]->GetNbinsY(); ++iy) {
                    double v = hRatio[i]->GetBinContent(ix, iy);
                    if (v == 0) continue;
                    s += v; s2 += v*v; cnt++;
                }
            }
            mean[i]  = cnt ? s/cnt : 0;
            sigma[i] = cnt ? sqrt(s2/cnt - mean[i]*mean[i]) : 0;
        }

        // 11) Ratio canvas with stats
        {
            TCanvas c("c_ratio", "Data/MC Ratio", 1800, 600);
            c.Divide(3,1);
            SetSame2DScale(hRatio[0], hRatio[1], hRatio[2]);
            for (int i = 0; i < 3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                hRatio[i]->Draw("COLZ");
                TLatex tex;
                tex.SetNDC(); tex.SetTextAlign(31); tex.SetTextSize(0.03);
                tex.DrawLatex(0.95, 0.95, Form("μ=%.3f", mean[i]));
                tex.DrawLatex(0.95, 0.90, Form("σ=%.3f", sigma[i]));
            }
            c.SaveAs(Form("output/cal/ratio_%s.png", label.c_str()));
        }

        // 12) Outliers canvas: highlight bins with ratio >2 or <0.5
        {
            TCanvas c("c_outliers", "Outliers (r>2 or r<0.5)", 1800, 600);
            c.Divide(3,1);
            SetSame2DScale(hRatio[0], hRatio[1], hRatio[2]);
            for (int i = 0; i < 3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                hRatio[i]->Draw("COLZ");
                for (int ix = 1; ix <= hRatio[i]->GetNbinsX(); ++ix) {
                    for (int iy = 1; iy <= hRatio[i]->GetNbinsY(); ++iy) {
                        double v = hRatio[i]->GetBinContent(ix, iy);
                        if (v == 0) continue;
                        if (v > 2.0 || v < 0.5) {
                            double x1 = hRatio[i]->GetXaxis()->GetBinLowEdge(ix);
                            double x2 = hRatio[i]->GetXaxis()->GetBinUpEdge(ix);
                            double y1 = hRatio[i]->GetYaxis()->GetBinLowEdge(iy);
                            double y2 = hRatio[i]->GetYaxis()->GetBinUpEdge(iy);
                            TBox box(x1, y1, x2, y2);
                            box.SetFillColor(kRed);
                            box.SetFillStyle(1001);
                            box.SetLineColor(kRed);
                            box.Draw("same");
                        }
                    }
                }
            }
            c.SaveAs(Form("output/cal/ratio_%s_outliers.png", label.c_str()));
        }
    }
    return 0;
}