/*
 * File: calorimeter_comparison.cpp
 *
 * Compile with:
 *   g++ calorimeter_comparison.cpp -o calorimeter_comparison `root-config --cflags --libs`
 *
 * Run with:
 *   ./calorimeter_comparison [Nevents] [dataFile] [mcFile]
 *
 * What it does:
 *   1) Creates output/ and output/cal/ if they do not exist.
 *   2) Loops over data and MC, filling photon (PID=22) and electron (PID=11)
 *      calorimeter hit-position histograms for PCal, ECin, and ECout.
 *   3) Makes unnormalized 2D plots (data vs. MC) with three panels (PCal/ECin/ECout).
 *   4) Normalizes each histogram and computes data/MC ratio.
 *   5) Creates a 3-panel ratio plot with mean (μ) and stddev (σ) in a legend box.
 *   6) Creates a second 3-panel "ratio_outliers" map histogram using two colors:
 *      blue where 0.5 ≤ ratio ≤ 2, red otherwise.
 *   7) Adds padding margins so y-axis labels are not clipped.
 */

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "TChain.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TSystem.h"
#include "TLegend.h"

// Helper to unify the color scale of three TH2D histograms
void SetSame2DScale(TH2D* a, TH2D* b, TH2D* c) {
    double mn = 1e9, mx = -1e9;
    for (auto* h : {a,b,c}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto* h : {a,b,c}) {
        h->SetMinimum(mn);
        h->SetMaximum(mx);
    }
}

int main(int argc, char** argv) {
    // Create output directories
    gSystem->mkdir("output",    kTRUE);
    gSystem->mkdir("output/cal",kTRUE);

    // Parse arguments
    Long64_t maxEvents = -1;
    if (argc > 1) { maxEvents = std::stoll(argv[1]); if (maxEvents == 0) maxEvents = -1; }
    bool useData = (argc > 2);
    bool useMC   = (argc > 3);
    std::string dataFile = useData ? argv[2] : "";
    std::string mcFile   = useMC   ? argv[3] : "";

    // Set up TChains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if (useData) dataCh.Add(dataFile.c_str()); else dataCh.Add("/work/clas12/thayward/.../cal_data.root");
    if (useMC)   mcCh.Add(mcFile.c_str());   else mcCh.Add("/work/clas12/thayward/.../cal_mc.root");

    // Branch addresses
    Int_t    pid;
    Double_t cal_x_1,cal_y_1,cal_x_4,cal_y_4,cal_x_7,cal_y_7;
    dataCh.SetBranchAddress("particle_pid", &pid);
    dataCh.SetBranchAddress("cal_x_1",      &cal_x_1);
    dataCh.SetBranchAddress("cal_y_1",      &cal_y_1);
    dataCh.SetBranchAddress("cal_x_4",      &cal_x_4);
    dataCh.SetBranchAddress("cal_y_4",      &cal_y_4);
    dataCh.SetBranchAddress("cal_x_7",      &cal_x_7);
    dataCh.SetBranchAddress("cal_y_7",      &cal_y_7);
    mcCh .SetBranchAddress("particle_pid", &pid);
    mcCh .SetBranchAddress("cal_x_1",      &cal_x_1);
    mcCh .SetBranchAddress("cal_y_1",      &cal_y_1);
    mcCh .SetBranchAddress("cal_x_4",      &cal_x_4);
    mcCh .SetBranchAddress("cal_y_4",      &cal_y_4);
    mcCh .SetBranchAddress("cal_x_7",      &cal_x_7);
    mcCh .SetBranchAddress("cal_y_7",      &cal_y_7);

    // Config
    const int NB = 200;
    const double xmin=-450, xmax=450, ymin=-450, ymax=450;
    std::vector<std::pair<int,std::string>> species = {{22,"photon"},{11,"electron"}};
    std::vector<std::string> layers = {"PCal","ECin","ECout"};

    for (auto& sp : species) {
        int pidVal = sp.first;
        std::string label = sp.second;

        // Default palette for regular plots
        gStyle->SetOptStat(0);
        gStyle->SetPalette(1);

        // Book
        std::vector<TH2D*> hD(3), hM(3), hR(3);
        for (int i=0; i<3; ++i) {
            hD[i] = new TH2D(Form("hD_%s_%s",label.c_str(),layers[i].c_str()),
                Form("%s Data %s; x (cm); y (cm)",label.c_str(),layers[i].c_str()),
                NB,xmin,xmax,NB,ymin,ymax);
            hM[i] = (TH2D*)hD[i]->Clone(Form("hM_%s_%s",label.c_str(),layers[i].c_str()));
        }
        // Fill data
        Long64_t nD = dataCh.GetEntries(); if (maxEvents>0 && maxEvents<nD) nD = maxEvents;
        for (Long64_t ev=0; ev<nD; ++ev) {
            dataCh.GetEntry(ev);
            if (pid != pidVal) continue;
            double xs[3] = {cal_x_1,cal_x_4,cal_x_7};
            double ys[3] = {cal_y_1,cal_y_4,cal_y_7};
            for (int i=0; i<3; ++i) if (xs[i]!=-9999 && ys[i]!=-9999) hD[i]->Fill(xs[i],ys[i]);
        }
        // Fill MC
        Long64_t nM = mcCh.GetEntries(); if (maxEvents>0 && maxEvents<nM) nM = maxEvents;
        for (Long64_t ev=0; ev<nM; ++ev) {
            mcCh.GetEntry(ev);
            if (pid != pidVal) continue;
            double xs[3] = {cal_x_1,cal_x_4,cal_x_7};
            double ys[3] = {cal_y_1,cal_y_4,cal_y_7};
            for (int i=0; i<3; ++i) if (xs[i]!=-9999 && ys[i]!=-9999) hM[i]->Fill(xs[i],ys[i]);
        }
        // Unnormalized Data
        TCanvas c1("c_data_uncut","Data Uncut",1800,600); c1.Divide(3,1);
        SetSame2DScale(hD[0],hD[1],hD[2]);
        for (int i=0; i<3; ++i) {
            c1.cd(i+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15);
            gPad->SetLogz(); hD[i]->Draw("COLZ");
        }
        c1.SaveAs(Form("output/cal/data_uncut_%s.png",label.c_str()));
        // Unnormalized MC
        TCanvas c2("c_mc_uncut","MC Uncut",1800,600); c2.Divide(3,1);
        SetSame2DScale(hM[0],hM[1],hM[2]);
        for (int i=0; i<3; ++i) {
            c2.cd(i+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15);
            gPad->SetLogz(); hM[i]->Draw("COLZ");
        }
        c2.SaveAs(Form("output/cal/mc_uncut_%s.png",label.c_str()));
        // Normalize & ratio
        for (int i=0; i<3; ++i) {
            double id = hD[i]->Integral(); if (id>0) hD[i]->Scale(1.0/id);
            double im = hM[i]->Integral(); if (im>0) hM[i]->Scale(1.0/im);
            hR[i] = (TH2D*)hD[i]->Clone(Form("hR_%s_%s",label.c_str(),layers[i].c_str()));
            hR[i]->Divide(hM[i]);
        }
        // Stats
        std::vector<double> mu(3), sigma(3);
        for (int i=0; i<3; ++i) {
            int cnt=0; double s=0, s2=0;
            for (int ix=1; ix<=NB; ++ix) for (int iy=1; iy<=NB; ++iy) {
                double v=hR[i]->GetBinContent(ix,iy);
                if (v<=0) continue;
                s+=v; s2+=v*v; cnt++;
            }
            mu[i]    = cnt? s/cnt : 0;
            sigma[i] = cnt? sqrt(s2/cnt - mu[i]*mu[i]) : 0;
        }
        // Ratio
        gStyle->SetPalette(1);
        TCanvas c3("c_ratio","Data/MC Ratio",1800,600); c3.Divide(3,1);
        SetSame2DScale(hR[0],hR[1],hR[2]);
        for (int i=0; i<3; ++i) {
            c3.cd(i+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15);
            gPad->SetLogz(); hR[i]->Draw("COLZ");
            TLegend leg(0.6, 0.7, 0.9, 0.9);
            leg.SetFillColor(kWhite); leg.SetBorderSize(1); leg.SetTextSize(0.03);
            leg.AddEntry((TObject*)0,Form("Mean=%.3f",mu[i]),"");
            leg.AddEntry((TObject*)0,Form("StdDev=%.3f",sigma[i]),"");
            leg.Draw();
        }
        c3.SaveAs(Form("output/cal/ratio_%s.png",label.c_str()));
        // Outliers map
        TCanvas c4("c_outliers","Outliers (map)",1800,600); c4.Divide(3,1);
        for (int i=0; i<3; ++i) {
            c4.cd(i+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15); gPad->SetLogz(0);
            TH2D* hMap = (TH2D*)hR[i]->Clone(Form("hMap_%s_%s",label.c_str(),layers[i].c_str()));
            hMap->Reset();
            for (int ix=1; ix<=NB; ++ix) {
                for (int iy=1; iy<=NB; ++iy) {
                    double v = hR[i]->GetBinContent(ix, iy);
                    if (v<=0) continue;
                    int lvl = (v<0.5||v>2.0)? 2 : 1;
                    hMap->SetBinContent(ix, iy, lvl);
                }
            }
            int palette[2] = {kBlue, kRed};
            gStyle->SetPalette(2, palette);
            hMap->SetContour(2);
            hMap->SetContourLevel(0,1);
            hMap->SetContourLevel(1,2);
            hMap->Draw("COLZ");
        }
        c4.SaveAs(Form("output/cal/ratio_%s_outliers.png",label.c_str()));
    }
    return 0;
}
