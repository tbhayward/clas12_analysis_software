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
 *   2) Applies fiducial cuts to calorimeter hits if desired.
 *   3) Loops once over data and once over MC, filling photon (PID=22)
 *      and electron (PID=11) hit-position histograms for PCal, ECin, ECout.
 *   4) Normalizes each histogram and computes data/MC ratio.
 *   5) Draws all regular "COLZ" canvases for both species: Data, MC, Ratio.
 *   6) Switches to two-color palette and draws outlier maps for both species.
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
#include "TLegend.h"

// ----------------------------------------------------------------------------
// Helper to unify the color scale of three TH2D histograms
// ----------------------------------------------------------------------------
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
    // create output directories
    gSystem->mkdir("output", kTRUE);
    gSystem->mkdir("output/cal", kTRUE);

    // parse args
    Long64_t maxEvents = -1;
    if (argc>1) { maxEvents = std::stoll(argv[1]); if(maxEvents==0) maxEvents=-1; }
    bool useData = (argc>2), useMC = (argc>3);
    std::string dataFile = useData? argv[2]:"";
    std::string mcFile   = useMC?   argv[3]:"";

    // setup chains
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if(useData) dataCh.Add(dataFile.c_str()); else dataCh.Add("/work/clas12/.../cal_data.root");
    if(useMC)   mcCh.Add(mcFile.c_str());   else mcCh.Add("/work/clas12/.../cal_mc.root");

    // branches
    Int_t pid;
    Double_t x1,y1,x4,y4,x7,y7;
    dataCh.SetBranchAddress("particle_pid", &pid);
    dataCh.SetBranchAddress("cal_x_1",      &x1);
    dataCh.SetBranchAddress("cal_y_1",      &y1);
    dataCh.SetBranchAddress("cal_x_4",      &x4);
    dataCh.SetBranchAddress("cal_y_4",      &y4);
    dataCh.SetBranchAddress("cal_x_7",      &x7);
    dataCh.SetBranchAddress("cal_y_7",      &y7);
    mcCh .SetBranchAddress("particle_pid", &pid);
    mcCh .SetBranchAddress("cal_x_1",      &x1);
    mcCh .SetBranchAddress("cal_y_1",      &y1);
    mcCh .SetBranchAddress("cal_x_4",      &x4);
    mcCh .SetBranchAddress("cal_y_4",      &y4);
    mcCh .SetBranchAddress("cal_x_7",      &x7);
    mcCh .SetBranchAddress("cal_y_7",      &y7);

    const int NB=200;
    const double xmin=-450, xmax=450, ymin=-450, ymax=450;
    std::vector<std::pair<int,std::string>> species = {{22,"photon"},{11,"electron"}};
    std::vector<std::string> layers = {"PCal","ECin","ECout"};

    // histos: [species][layer]
    std::vector<std::vector<TH2D*>> hData(2), hMC(2), hRatio(2);
    for(int s=0; s<2; ++s) {
        for(int i=0; i<3; ++i) {
            auto lbl = Form("%s_%s", species[s].second.c_str(), layers[i].c_str());
            hData[s].push_back(new TH2D(Form("hD_%s",lbl),Form("%s Data %s; x; y",species[s].second.c_str(),layers[i].c_str()),NB,xmin,xmax,NB,ymin,ymax));
            hMC  [s].push_back((TH2D*)hData[s][i]->Clone(Form("hM_%s",lbl)));
        }
    }

    // fill data
    Long64_t nD = dataCh.GetEntries(); if(maxEvents>0 && maxEvents<nD) nD=maxEvents;
    for(Long64_t ev=0; ev<nD; ++ev) {
        dataCh.GetEntry(ev);
        for(int s=0; s<2; ++s) {
            if(pid != species[s].first) continue;
            double xs[3]={x1,x4,x7}, ys[3]={y1,y4,y7};
            for(int i=0;i<3;i++) if(xs[i]!=-9999 && ys[i]!=-9999)
                hData[s][i]->Fill(xs[i], ys[i]);
        }
    }
    // fill MC
    Long64_t nM = mcCh.GetEntries(); if(maxEvents>0 && maxEvents<nM) nM=maxEvents;
    for(Long64_t ev=0; ev<nM; ++ev) {
        mcCh.GetEntry(ev);
        for(int s=0; s<2; ++s) {
            if(pid != species[s].first) continue;
            double xs[3]={x1,x4,x7}, ys[3]={y1,y4,y7};
            for(int i=0;i<3;i++) if(xs[i]!=-9999 && ys[i]!=-9999)
                hMC[s][i]->Fill(xs[i], ys[i]);
        }
    }

    // normalize & ratio
    std::vector<std::vector<double>> mu(2, std::vector<double>(3)), sigma(2, std::vector<double>(3));
    for(int s=0;s<2;++s) {
        for(int i=0;i<3;++i) {
            double id=hData[s][i]->Integral(); if(id>0) hData[s][i]->Scale(1/id);
            double im=hMC  [s][i]->Integral(); if(im>0) hMC  [s][i]->Scale(1/im);
            hRatio[s].push_back((TH2D*)hData[s][i]->Clone(Form("hR_%s_%s",species[s].second.c_str(),layers[i].c_str())));
            hRatio[s][i]->Divide(hMC[s][i]);
            // stats
            int cnt=0; double sum=0,sum2=0;
            for(int ix=1;ix<=NB;++ix) for(int iy=1;iy<=NB;++iy) {
                double v=hRatio[s][i]->GetBinContent(ix,iy);
                if(v<=0) continue; sum+=v; sum2+=v*v; cnt++;
            }
            mu[s][i]    = cnt? sum/cnt:0;
            sigma[s][i] = cnt? std::sqrt(sum2/cnt - mu[s][i]*mu[s][i]):0;
        }
    }

    gStyle->SetOptStat(0);
    // draw regular for both species
    for(int s=0;s<2;++s) {
        // Data uncut
        TCanvas cD(Form("c_data_%s",species[s].second.c_str()),Form("%s Data Uncut",species[s].second.c_str()),1800,600);
        cD.Divide(3,1);
        SetSame2DScale(hData[s][0],hData[s][1],hData[s][2]);
        for(int i=0;i<3;++i) { cD.cd(i+1); gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz(); hData[s][i]->Draw("COLZ"); }
        cD.SaveAs(Form("output/cal/data_uncut_%s.png",species[s].second.c_str()));
        // MC uncut
        TCanvas cM(Form("c_mc_%s",species[s].second.c_str()),Form("%s MC Uncut",species[s].second.c_str()),1800,600);
        cM.Divide(3,1);
        SetSame2DScale(hMC[s][0],hMC[s][1],hMC[s][2]);
        for(int i=0;i<3;++i) { cM.cd(i+1); gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz(); hMC[s][i]->Draw("COLZ"); }
        cM.SaveAs(Form("output/cal/mc_uncut_%s.png",species[s].second.c_str()));
        // Ratio
        TCanvas cR(Form("c_ratio_%s",species[s].second.c_str()),Form("%s Data/MC Ratio",species[s].second.c_str()),1800,600);
        cR.Divide(3,1);
        SetSame2DScale(hRatio[s][0],hRatio[s][1],hRatio[s][2]);
        for(int i=0;i<3;++i) {
            cR.cd(i+1);
            gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz();
            hRatio[s][i]->Draw("COLZ");
            TLegend leg(.6,.7,.9,.9); leg.SetFillColor(kWhite); leg.SetBorderSize(1); leg.SetTextSize(0.03);
            leg.AddEntry((TObject*)0,Form("Mean=%.3f",mu[s][i]),"");
            leg.AddEntry((TObject*)0,Form("StdDev=%.3f",sigma[s][i]),"");
            leg.Draw();
        }
        cR.SaveAs(Form("output/cal/ratio_%s.png",species[s].second.c_str()));
    }

    // two-color palette for outliers
    int outPal[2]={kBlue,kRed};
    gStyle->SetPalette(2,outPal);
    for(int s=0;s<2;++s) {
        TCanvas cO(Form("c_outliers_%s",species[s].second.c_str()),Form("%s Ratio Outliers",species[s].second.c_str()),1800,600);
        cO.Divide(3,1);
        for(int i=0;i<3;++i) {
            cO.cd(i+1);
            gPad->SetLeftMargin(.15); gPad->SetRightMargin(.15); gPad->SetLogz(0);
            TH2D* map = (TH2D*)hRatio[s][i]->Clone("map"); map->Reset();
            for(int ix=1;ix<=NB;++ix) for(int iy=1;iy<=NB;++iy) {
                double v=hRatio[s][i]->GetBinContent(ix,iy);
                if(v<=0) continue;
                int lvl = (v<0.5||v>2.0)?2:1;
                map->SetBinContent(ix,iy,lvl);
            }
            map->SetContour(2); map->SetContourLevel(0,1); map->SetContourLevel(1,2);
            map->Draw("COLZ");
        }
        cO.SaveAs(Form("output/cal/ratio_%s_outliers.png",species[s].second.c_str()));
    }

    return 0;
}
