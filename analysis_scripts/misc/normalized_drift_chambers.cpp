/*****************************************************************************
 * File: normalized_drift_chambers.cpp
 *
 * Compile with:
 *   g++ normalized_drift_chambers.cpp -std=c++11 -o normalized_drift_chambers \
 *       `root-config --cflags --libs`
 *
 * Run with:
 *   ./normalized_drift_chambers [NeventsSIDISDVCS] [dataFile] [mcFile]
 *
 * This version:
 *   • Loops once over data and once over MC, filling all histograms for both
 *     electrons (PID=11) and protons (PID=2212) in one pass.
 *   • Normalizes each individual histogram.
 *   • Computes data/MC ratio histograms for both species.
 *   • Saves all five canvases per species, suffixed with `_electron` or `_proton`.
 *****************************************************************************/

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include "TChain.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TEllipse.h"
#include "TString.h"

// ----------------------------------------------------------------------------
// Helpers to unify scales
// ----------------------------------------------------------------------------
void SetSame2DScale(TH2D* a, TH2D* b, TH2D* c) {
    double mn =  1e9, mx = -1e9;
    for (auto* h : {a,b,c}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto* h : {a,b,c}) {
        h->SetMinimum(mn);
        h->SetMaximum(mx);
    }
}
void SetSame1DScale(TH1* a, TH1* b, TH1* c) {
    double mn =  1e9, mx = -1e9;
    for (auto* h : {a,b,c}) {
        mn = std::min(mn, h->GetMinimum());
        mx = std::max(mx, h->GetMaximum());
    }
    for (auto* h : {a,b,c}) {
        h->SetMinimum(mn);
        h->SetMaximum(mx);
    }
}

// ----------------------------------------------------------------------------
// Bundle for each species
// ----------------------------------------------------------------------------
struct HistoSet {
    // data
    TH2D *d_r1, *d_r2, *d_r3;
    TH2D *d_r1c,*d_r2c,*d_r3c;
    TH1D *d_e1,  *d_e2,  *d_e3;
    TH2D *d_rp1, *d_rp2, *d_rp3;
    // mc
    TH2D *m_r1, *m_r2, *m_r3;
    TH2D *m_r1c,*m_r2c,*m_r3c;
    TH1D *m_e1,  *m_e2,  *m_e3;
    TH2D *m_rp1, *m_rp2, *m_rp3;
    // ratios
    TH2D *r_r1, *r_r2, *r_r3;
    TH2D *r_r1c,*r_r2c,*r_r3c;
    TH1D *r_e1, *r_e2, *r_e3;
    TH2D *r_rp1,*r_rp2,*r_rp3;
};

int main(int argc, char** argv) {
    //--------------------------------------------------------------------------
    // 1) Parse arguments
    //--------------------------------------------------------------------------
    Long64_t maxEvents = -1;
    if (argc>1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents==0) maxEvents=-1;
    }
    bool useDataFile = (argc>2), useMCFile = (argc>3);
    std::string dataFile = useDataFile? argv[2] : "";
    std::string mcFile   = useMCFile  ? argv[3] : "";

    //--------------------------------------------------------------------------
    // 2) Setup chains
    //--------------------------------------------------------------------------
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if (useDataFile) {
        dataCh.Add(dataFile.c_str());
    } else {
        dataCh.Add("/work/clas12/thayward/CLAS12_SIDIS/.../sidisdvcs_rgc_su22_inb_calibration.root");
        dataCh.Add("/work/clas12/thayward/CLAS12_SIDIS/.../sidisdvcs_rgc_su22_inb_calibration_1.root");
        dataCh.Add("/work/clas12/thayward/CLAS12_SIDIS/.../sidisdvcs_rgc_su22_inb_calibration_2.root");
    }
    if (useMCFile) {
        mcCh.Add(mcFile.c_str());
    } else {
        mcCh.Add("/work/clas12/thayward/CLAS12_SIDIS/.../clasdis_rgc_su22_inb_neutron_calibration.root");
    }

    //--------------------------------------------------------------------------
    // 3) Branch addresses
    //--------------------------------------------------------------------------
    Int_t    pid;
    Double_t x6,y6, x18,y18, x36,y36;
    Double_t e6,e18,e36;
    for (auto* ch : {&dataCh, &mcCh}) {
        ch->SetBranchAddress("particle_pid", &pid);
        ch->SetBranchAddress("traj_x_6",      &x6);
        ch->SetBranchAddress("traj_y_6",      &y6);
        ch->SetBranchAddress("traj_x_18",     &x18);
        ch->SetBranchAddress("traj_y_18",     &y18);
        ch->SetBranchAddress("traj_x_36",     &x36);
        ch->SetBranchAddress("traj_y_36",     &y36);
        ch->SetBranchAddress("traj_edge_6",   &e6);
        ch->SetBranchAddress("traj_edge_18",  &e18);
        ch->SetBranchAddress("traj_edge_36",  &e36);
    }

    //--------------------------------------------------------------------------
    // 4) Constants & prepare histo-sets
    //--------------------------------------------------------------------------
    const double R1=140.0, R2=215.0, R3=290.0;
    const int NB2=200, NB1=125, nbR1=150, nbR2=230, nbR3=400, nbPhi=180;

    // species vector: { PID, label }
    std::vector<std::pair<int,const char*>> species = {
        {11,"electron"}, {2212,"proton"}
    };
    std::map<int,HistoSet> H;

    // create for each species
    for (auto& sp : species) {
        int p=sp.first; const char* lab=sp.second;
        HistoSet &h = H[p];

        // data
        h.d_r1  = new TH2D(Form("d_r1_%s",lab), Form("Data R1 (%s); x; y",lab),
                           NB2,-180,180, NB2,-180,180);
        h.d_r2  = new TH2D(Form("d_r2_%s",lab), Form("Data R2 (%s); x; y",lab),
                           NB2,-280,280, NB2,-280,280);
        h.d_r3  = new TH2D(Form("d_r3_%s",lab), Form("Data R3 (%s); x; y",lab),
                           NB2,-450,450, NB2,-450,450);
        h.d_r1c = (TH2D*)h.d_r1->Clone(Form("d_r1c_%s",lab));
        h.d_r2c = (TH2D*)h.d_r2->Clone(Form("d_r2c_%s",lab));
        h.d_r3c = (TH2D*)h.d_r3->Clone(Form("d_r3c_%s",lab));
        h.d_e1  = new TH1D(Form("d_e1_%s",lab),Form("Data edge R1 (%s); edge (cm)",lab),
                           NB1,0,125);
        h.d_e2  = new TH1D(Form("d_e2_%s",lab),Form("Data edge R2 (%s); edge (cm)",lab),
                           NB1,0,125);
        h.d_e3  = new TH1D(Form("d_e3_%s",lab),Form("Data edge R3 (%s); edge (cm)",lab),
                           NB1,0,125);
        h.d_rp1 = new TH2D(Form("d_rp1_%s",lab),Form("Data R1 ρ-φ (%s); ρ; φ (deg)",lab),
                           nbR1,20,170, nbPhi,0,360);
        h.d_rp2 = new TH2D(Form("d_rp2_%s",lab),Form("Data R2 ρ-φ (%s); ρ; φ (deg)",lab),
                           nbR2,20,250, nbPhi,0,360);
        h.d_rp3 = new TH2D(Form("d_rp3_%s",lab),Form("Data R3 ρ-φ (%s); ρ; φ (deg)",lab),
                           nbR3,20,420, nbPhi,0,360);

        // mc = clones of data shapes
        h.m_r1  = (TH2D*)h.d_r1 ->Clone(Form("m_r1_%s",lab));
        h.m_r2  = (TH2D*)h.d_r2 ->Clone(Form("m_r2_%s",lab));
        h.m_r3  = (TH2D*)h.d_r3 ->Clone(Form("m_r3_%s",lab));
        h.m_r1c = (TH2D*)h.d_r1c->Clone(Form("m_r1c_%s",lab));
        h.m_r2c = (TH2D*)h.d_r2c->Clone(Form("m_r2c_%s",lab));
        h.m_r3c = (TH2D*)h.d_r3c->Clone(Form("m_r3c_%s",lab));
        h.m_e1  = (TH1D*)h.d_e1 ->Clone(Form("m_e1_%s",lab));
        h.m_e2  = (TH1D*)h.d_e2 ->Clone(Form("m_e2_%s",lab));
        h.m_e3  = (TH1D*)h.d_e3 ->Clone(Form("m_e3_%s",lab));
        h.m_rp1 = (TH2D*)h.d_rp1->Clone(Form("m_rp1_%s",lab));
        h.m_rp2 = (TH2D*)h.d_rp2->Clone(Form("m_rp2_%s",lab));
        h.m_rp3 = (TH2D*)h.d_rp3->Clone(Form("m_rp3_%s",lab));

        // ratios: clones of data that we will Divide(mc)
        h.r_r1  = (TH2D*)h.d_r1 ->Clone(Form("r_r1_%s",lab));
        h.r_r2  = (TH2D*)h.d_r2 ->Clone(Form("r_r2_%s",lab));
        h.r_r3  = (TH2D*)h.d_r3 ->Clone(Form("r_r3_%s",lab));
        h.r_r1c = (TH2D*)h.d_r1c->Clone(Form("r_r1c_%s",lab));
        h.r_r2c = (TH2D*)h.d_r2c->Clone(Form("r_r2c_%s",lab));
        h.r_r3c = (TH2D*)h.d_r3c->Clone(Form("r_r3c_%s",lab));
        h.r_e1  = (TH1D*)h.d_e1 ->Clone(Form("r_e1_%s",lab));
        h.r_e2  = (TH1D*)h.d_e2 ->Clone(Form("r_e2_%s",lab));
        h.r_e3  = (TH1D*)h.d_e3 ->Clone(Form("r_e3_%s",lab));
        h.r_rp1 = (TH2D*)h.d_rp1->Clone(Form("r_rp1_%s",lab));
        h.r_rp2 = (TH2D*)h.d_rp2->Clone(Form("r_rp2_%s",lab));
        h.r_rp3 = (TH2D*)h.d_rp3->Clone(Form("r_rp3_%s",lab));
    }

    gStyle->SetOptStat(0);

    //--------------------------------------------------------------------------
    // 5) Fill DATA in one pass
    //--------------------------------------------------------------------------
    Long64_t nD = dataCh.GetEntries();
    if (maxEvents>0 && maxEvents<nD) nD = maxEvents;
    for (Long64_t i=0; i<nD; ++i) {
        dataCh.GetEntry(i);
        // only electron or proton
        if (!H.count(pid)) continue;
        auto& h = H[pid];
        // edges
        if (e6<=3||e18<=3||e36<=10) continue;
        // R1
        if (x6!=-9999) {
            h.d_r1->Fill(x6,y6);
            double r=std::hypot(x6,y6);
            if (r<R1) h.d_r1c->Fill(x6,y6);
            double phi = std::atan2(y6,x6)-M_PI; if(phi<0) phi+=2*M_PI;
            h.d_rp1->Fill(r,phi*180./M_PI);
        }
        // R2
        if (x18!=-9999) {
            h.d_r2->Fill(x18,y18);
            double r=std::hypot(x18,y18);
            if(r<R2) h.d_r2c->Fill(x18,y18);
            double phi=std::atan2(y18,x18)-M_PI; if(phi<0)phi+=2*M_PI;
            h.d_rp2->Fill(r,phi*180./M_PI);
        }
        // R3
        if (x36!=-9999) {
            h.d_r3->Fill(x36,y36);
            double r=std::hypot(x36,y36);
            if(r<R3) h.d_r3c->Fill(x36,y36);
            double phi=std::atan2(y36,x36)-M_PI; if(phi<0)phi+=2*M_PI;
            h.d_rp3->Fill(r,phi*180./M_PI);
        }
        h.d_e1->Fill(e6);
        h.d_e2->Fill(e18);
        h.d_e3->Fill(e36);
    }

    //--------------------------------------------------------------------------
    // 6) Fill MC in one pass
    //--------------------------------------------------------------------------
    Long64_t nM = mcCh.GetEntries();
    if (maxEvents>0 && maxEvents<nM) nM = maxEvents;
    for (Long64_t i=0; i<nM; ++i) {
        mcCh.GetEntry(i);
        if (!H.count(pid)) continue;
        auto& h = H[pid];
        if (e6<=3||e18<=3||e36<=10) continue;
        if (x6!=-9999) {
            h.m_r1->Fill(x6,y6);
            double r=std::hypot(x6,y6);
            if(r<R1) h.m_r1c->Fill(x6,y6);
            double phi=std::atan2(y6,x6)-M_PI; if(phi<0)phi+=2*M_PI;
            h.m_rp1->Fill(r,phi*180./M_PI);
        }
        if (x18!=-9999) {
            h.m_r2->Fill(x18,y18);
            double r=std::hypot(x18,y18);
            if(r<R2) h.m_r2c->Fill(x18,y18);
            double phi=std::atan2(y18,x18)-M_PI; if(phi<0)phi+=2*M_PI;
            h.m_rp2->Fill(r,phi*180./M_PI);
        }
        if (x36!=-9999) {
            h.m_r3->Fill(x36,y36);
            double r=std::hypot(x36,y36);
            if(r<R3) h.m_r3c->Fill(x36,y36);
            double phi=std::atan2(y36,x36)-M_PI; if(phi<0)phi+=2*M_PI;
            h.m_rp3->Fill(r,phi*180./M_PI);
        }
        h.m_e1->Fill(e6);
        h.m_e2->Fill(e18);
        h.m_e3->Fill(e36);
    }

    //--------------------------------------------------------------------------
    // 7) Normalize & form ratio, then draw per species
    //--------------------------------------------------------------------------
    for (auto& sp : species) {
        int p = sp.first; const char* lab = sp.second;
        auto& h = H[p];

        // normalize data
        for (auto* hh : { (TH1*)h.d_r1, (TH1*)h.d_r2, (TH1*)h.d_r3,
                          (TH1*)h.d_r1c,(TH1*)h.d_r2c,(TH1*)h.d_r3c,
                          h.d_e1,h.d_e2,h.d_e3,
                          (TH1*)h.d_rp1,(TH1*)h.d_rp2,(TH1*)h.d_rp3 }) {
            double I=hh->Integral(); if(I>0) hh->Scale(1./I);
        }
        // normalize mc
        for (auto* hh : { (TH1*)h.m_r1, (TH1*)h.m_r2, (TH1*)h.m_r3,
                          (TH1*)h.m_r1c,(TH1*)h.m_r2c,(TH1*)h.m_r3c,
                          h.m_e1,h.m_e2,h.m_e3,
                          (TH1*)h.m_rp1,(TH1*)h.m_rp2,(TH1*)h.m_rp3 }) {
            double I=hh->Integral(); if(I>0) hh->Scale(1./I);
        }

        // form ratios: data / mc
        h.r_r1 ->Divide(h.m_r1 );
        h.r_r2 ->Divide(h.m_r2 );
        h.r_r3 ->Divide(h.m_r3 );
        h.r_r1c->Divide(h.m_r1c);
        h.r_r2c->Divide(h.m_r2c);
        h.r_r3c->Divide(h.m_r3c);
        h.r_e1 ->Divide(h.m_e1 );
        h.r_e2 ->Divide(h.m_e2 );
        h.r_e3 ->Divide(h.m_e3 );
        h.r_rp1->Divide(h.m_rp1);
        h.r_rp2->Divide(h.m_rp2);
        h.r_rp3->Divide(h.m_rp3);

        // unify scales on ratio canvases
        SetSame2DScale(h.r_r1,h.r_r2,h.r_r3);
        SetSame2DScale(h.r_r1c,h.r_r2c,h.r_r3c);
        SetSame1DScale(h.r_e1,h.r_e2,h.r_e3);
        SetSame2DScale(h.r_rp1,h.r_rp2,h.r_rp3);

        // -- 1) uncut --
        {
            TCanvas c(Form("ratio_uncut_%s",lab),Form("%s Ratio Uncut",lab),1800,600);
            c.Divide(3,1);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? h.r_r1 : i==1? h.r_r2 : h.r_r3)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/ratio_uncut_%s.png",lab));
        }
        // -- 2) circle overlay --
        {
            TCanvas c(Form("ratio_circle_%s",lab),Form("%s Ratio Circle",lab),1800,600);
            c.Divide(3,1);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                TH2D* hh = i==0? h.r_r1: i==1? h.r_r2: h.r_r3;
                double R = i==0? R1: i==1? R2: R3;
                hh->Draw("COLZ");
                TEllipse cir(0,0,R,R);
                cir.SetLineColor(kRed);
                cir.SetLineWidth(2);
                cir.SetFillStyle(0);
                cir.Draw("same");
            }
            c.SaveAs(Form("output/normalization/ratio_circle_%s.png",lab));
        }
        // -- 3) cut --
        {
            TCanvas c(Form("ratio_cut_%s",lab),Form("%s Ratio Cut",lab),1800,600);
            c.Divide(3,1);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? h.r_r1c : i==1? h.r_r2c : h.r_r3c)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/ratio_cut_%s.png",lab));
        }
        // -- 4) edges --
        {
            TCanvas c(Form("ratio_edges_%s",lab),Form("%s Ratio Edges",lab),1800,600);
            c.Divide(3,1);
            c.cd(1); h.r_e1->Draw("HIST");
            c.cd(2); h.r_e2->Draw("HIST");
            c.cd(3); h.r_e3->Draw("HIST");
            c.SaveAs(Form("output/normalization/ratio_edges_%s.png",lab));
        }
        // -- 5) rad-phi --
        {
            TCanvas c(Form("ratio_radphi_%s",lab),Form("%s Ratio ρ-φ",lab),1800,600);
            c.Divide(3,1);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? h.r_rp1 : i==1? h.r_rp2 : h.r_rp3)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/ratio_radphi_%s.png",lab));
        }
    }

    return 0;
}