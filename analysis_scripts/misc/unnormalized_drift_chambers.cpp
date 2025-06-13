/*****************************************************************************
 * File: unnormalized_drift_chambers.cpp
 *
 * Compile with:
 *   g++ unnormalized_drift_chambers.cpp -std=c++11 -o unnormalized_drift_chambers \
 *       `root-config --cflags --libs`
 *
 * Run with:
 *   ./unnormalized_drift_chambers [Nevents] [dataFile] [mcFile]
 *
 * This version:
 *   • Loops once over data and once over MC, filling both electron (PID=11)
 *     and proton (PID=2212) histograms in a single pass.
 *   • Normalizes and draws separate canvases for electrons and protons.
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
// Bundle of histograms for one PID
// ----------------------------------------------------------------------------
struct HistoSet {
    TH2D *r1, *r2, *r3;
    TH2D *r1c,*r2c,*r3c;
    TH1D *e1, *e2, *e3;
    TH2D *rp1,*rp2,*rp3;
};

int main(int argc, char** argv) {
    // ------------------------------------------------------------------------
    // 1) Parse arguments
    // ------------------------------------------------------------------------
    Long64_t maxEvents = -1;
    if (argc > 1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents == 0) maxEvents = -1;
    }
    bool useDataFile = (argc > 2);
    std::string dataFile = useDataFile ? argv[2] : "";
    bool useMCFile   = (argc > 3);
    std::string mcFile   = useMCFile   ? argv[3] : "";

    // ------------------------------------------------------------------------
    // 2) Set up TChains for data and MC
    // ------------------------------------------------------------------------
    TChain dataCh("PhysicsEvents"), mcCh("PhysicsEvents");
    if (useDataFile) {
        dataCh.Add(dataFile.c_str());
    } else {
        dataCh.Add("/work/clas12/thayward/.../sidisdvcs_0.root");
        dataCh.Add("/work/clas12/thayward/.../sidisdvcs_1.root");
        dataCh.Add("/work/clas12/thayward/.../sidisdvcs_2.root");
    }
    if (useMCFile) {
        mcCh.Add(mcFile.c_str());
    } else {
        mcCh.Add("/work/clas12/thayward/.../clasdis.root");
    }

    // ------------------------------------------------------------------------
    // 3) Branch addresses
    // ------------------------------------------------------------------------
    Int_t    pid;
    Double_t x6,y6, x18,y18, x36,y36;
    Double_t e6,e18,e36;
    dataCh.SetBranchAddress("particle_pid",&pid);
    dataCh.SetBranchAddress("traj_x_6",&x6);
    dataCh.SetBranchAddress("traj_y_6",&y6);
    dataCh.SetBranchAddress("traj_x_18",&x18);
    dataCh.SetBranchAddress("traj_y_18",&y18);
    dataCh.SetBranchAddress("traj_x_36",&x36);
    dataCh.SetBranchAddress("traj_y_36",&y36);
    dataCh.SetBranchAddress("traj_edge_6",&e6);
    dataCh.SetBranchAddress("traj_edge_18",&e18);
    dataCh.SetBranchAddress("traj_edge_36",&e36);

    mcCh.SetBranchAddress("particle_pid",&pid);
    mcCh.SetBranchAddress("traj_x_6",&x6);
    mcCh.SetBranchAddress("traj_y_6",&y6);
    mcCh.SetBranchAddress("traj_x_18",&x18);
    mcCh.SetBranchAddress("traj_y_18",&y18);
    mcCh.SetBranchAddress("traj_x_36",&x36);
    mcCh.SetBranchAddress("traj_y_36",&y36);
    mcCh.SetBranchAddress("traj_edge_6",&e6);
    mcCh.SetBranchAddress("traj_edge_18",&e18);
    mcCh.SetBranchAddress("traj_edge_36",&e36);

    // ------------------------------------------------------------------------
    // 4) Constants and histogram maps
    // ------------------------------------------------------------------------
    const double R1 = 140.0, R2 = 215.0, R3 = 290.0;
    const int NB2 = 200, NB1 = 125;
    const int nbR1 = 150, nbR2 = 230, nbR3 = 400, nbPhi = 180;

    // We will build histograms for electrons (11) and protons (2212)
    std::vector<std::pair<int,const char*>> species = {
        { 11,   "electron" },
        { 2212, "proton"   }
    };
    std::cout << " HELLO WORLD " << std::endl;

    std::map<int,HistoSet> Hdata, Hmc;

    // Create and clone histograms
    for (auto& sp : species) {
        int pidval   = sp.first;
        const char* label = sp.second;

        // Data histos
        HistoSet &d = Hdata[pidval];
        d.r1  = new TH2D(Form("r1_%s_data", label), Form("Data R1 (%s); x; y", label),
                         NB2, -180, 180, NB2, -180, 180);
        d.r2  = new TH2D(Form("r2_%s_data", label), Form("Data R2 (%s); x; y", label),
                         NB2, -280, 280, NB2, -280, 280);
        d.r3  = new TH2D(Form("r3_%s_data", label), Form("Data R3 (%s); x; y", label),
                         NB2, -450, 450, NB2, -450, 450);

        d.r1c = new TH2D(Form("r1c_%s_data", label), Form("Data R1 cut (%s); x; y", label),
                         NB2, -180, 180, NB2, -180, 180);
        d.r2c = new TH2D(Form("r2c_%s_data", label), Form("Data R2 cut (%s); x; y", label),
                         NB2, -280, 280, NB2, -280, 280);
        d.r3c = new TH2D(Form("r3c_%s_data", label), Form("Data R3 cut (%s); x; y", label),
                         NB2, -450, 450, NB2, -450, 450);

        d.e1  = new TH1D(Form("e1_%s_data", label), Form("Data edge R1 (%s); edge (cm)", label),
                         NB1, 0, 125);
        d.e2  = new TH1D(Form("e2_%s_data", label), Form("Data edge R2 (%s); edge (cm)", label),
                         NB1, 0, 125);
        d.e3  = new TH1D(Form("e3_%s_data", label), Form("Data edge R3 (%s); edge (cm)", label),
                         NB1, 0, 125);

        d.rp1 = new TH2D(Form("rp1_%s_data", label), Form("Data R1 ρ-φ (%s); ρ; φ (deg)", label),
                         nbR1, 20, 170, nbPhi, 0, 360);
        d.rp2 = new TH2D(Form("rp2_%s_data", label), Form("Data R2 ρ-φ (%s); ρ; φ (deg)", label),
                         nbR2, 20, 250, nbPhi, 0, 360);
        d.rp3 = new TH2D(Form("rp3_%s_data", label), Form("Data R3 ρ-φ (%s); ρ; φ (deg)", label),
                         nbR3, 20, 420, nbPhi, 0, 360);

        // MC histos: clone each Data histogram into MC version
        HistoSet &m = Hmc[pidval];
        m.r1  = (TH2D*)d.r1 ->Clone(Form("r1_%s_mc",  label));
        m.r2  = (TH2D*)d.r2 ->Clone(Form("r2_%s_mc",  label));
        m.r3  = (TH2D*)d.r3 ->Clone(Form("r3_%s_mc",  label));
        m.r1c = (TH2D*)d.r1c->Clone(Form("r1c_%s_mc", label));
        m.r2c = (TH2D*)d.r2c->Clone(Form("r2c_%s_mc", label));
        m.r3c = (TH2D*)d.r3c->Clone(Form("r3c_%s_mc", label));
        m.e1  = (TH1D*)d.e1 ->Clone(Form("e1_%s_mc",  label));
        m.e2  = (TH1D*)d.e2 ->Clone(Form("e2_%s_mc",  label));
        m.e3  = (TH1D*)d.e3 ->Clone(Form("e3_%s_mc",  label));
        m.rp1 = (TH2D*)d.rp1->Clone(Form("rp1_%s_mc", label));
        m.rp2 = (TH2D*)d.rp2->Clone(Form("rp2_%s_mc", label));
        m.rp3 = (TH2D*)d.rp3->Clone(Form("rp3_%s_mc", label));
    }

    gStyle->SetOptStat(0);

    // ------------------------------------------------------------------------
    // 5) Fill DATA histograms in one pass
    // ------------------------------------------------------------------------
    Long64_t nD = dataCh.GetEntries();
    if (maxEvents > 0 && maxEvents < nD) nD = maxEvents;
    for (Long64_t i = 0; i < nD; ++i) {
        dataCh.GetEntry(i);

        // only electrons or protons
        if (!Hdata.count(pid)) continue;
        HistoSet &h = Hdata[pid];

        // edge cuts
        if (e6 <= 3 || e18 <= 3 || e36 <= 10) continue;

        // Region 1
        if (x6 != -9999) {
            h.r1->Fill(x6, y6);
            double r = std::hypot(x6, y6);
            if (r < R1) h.r1c->Fill(x6, y6);
            double phi = std::atan2(y6, x6) - M_PI; if (phi < 0) phi += 2*M_PI;
            h.rp1->Fill(r, phi * 180.0 / M_PI);
        }
        // Region 2
        if (x18 != -9999) {
            h.r2->Fill(x18, y18);
            double r = std::hypot(x18, y18);
            if (r < R2) h.r2c->Fill(x18, y18);
            double phi = std::atan2(y18, x18) - M_PI; if (phi < 0) phi += 2*M_PI;
            h.rp2->Fill(r, phi * 180.0 / M_PI);
        }
        // Region 3
        if (x36 != -9999) {
            h.r3->Fill(x36, y36);
            double r = std::hypot(x36, y36);
            if (r < R3) h.r3c->Fill(x36, y36);
            double phi = std::atan2(y36, x36) - M_PI; if (phi < 0) phi += 2*M_PI;
            h.rp3->Fill(r, phi * 180.0 / M_PI);
        }

        h.e1->Fill(e6);
        h.e2->Fill(e18);
        h.e3->Fill(e36);
    }

    // ------------------------------------------------------------------------
    // 6) Fill MC histograms in one pass
    // ------------------------------------------------------------------------
    Long64_t nM = mcCh.GetEntries();
    if (maxEvents > 0 && maxEvents < nM) nM = maxEvents;
    for (Long64_t i = 0; i < nM; ++i) {
        mcCh.GetEntry(i);

        if (!Hmc.count(pid)) continue;
        HistoSet &h = Hmc[pid];

        if (e6 <= 3 || e18 <= 3 || e36 <= 10) continue;

        if (x6 != -9999) {
            h.r1->Fill(x6, y6);
            double r = std::hypot(x6, y6);
            if (r < R1) h.r1c->Fill(x6, y6);
            double phi = std::atan2(y6, x6) - M_PI; if (phi < 0) phi += 2*M_PI;
            h.rp1->Fill(r, phi * 180.0 / M_PI);
        }
        if (x18 != -9999) {
            h.r2->Fill(x18, y18);
            double r = std::hypot(x18, y18);
            if (r < R2) h.r2c->Fill(x18, y18);
            double phi = std::atan2(y18, x18) - M_PI; if (phi < 0) phi += 2*M_PI;
            h.rp2->Fill(r, phi * 180.0 / M_PI);
        }
        if (x36 != -9999) {
            h.r3->Fill(x36, y36);
            double r = std::hypot(x36, y36);
            if (r < R3) h.r3c->Fill(x36, y36);
            double phi = std::atan2(y36, x36) - M_PI; if (phi < 0) phi += 2*M_PI;
            h.rp3->Fill(r, phi * 180.0 / M_PI);
        }

        h.e1->Fill(e6);
        h.e2->Fill(e18);
        h.e3->Fill(e36);
    }

    // ------------------------------------------------------------------------
    // 7) Normalize & draw for each species
    // ------------------------------------------------------------------------
    for (auto& sp : species) {
        int pidval    = sp.first;
        const char* label = sp.second;
        HistoSet &d = Hdata[pidval];
        HistoSet &m = Hmc[   pidval];

        // normalize data
        for (auto* h : std::vector<TH1*>{ (TH1*)d.r1,(TH1*)d.r2,(TH1*)d.r3,
                                          (TH1*)d.r1c,(TH1*)d.r2c,(TH1*)d.r3c,
                                          d.e1,d.e2,d.e3,
                                          (TH1*)d.rp1,(TH1*)d.rp2,(TH1*)d.rp3 }) {
            double I = h->Integral();
            if (I>0) h->Scale(1.0/I);
        }
        // normalize MC
        for (auto* h : std::vector<TH1*>{ (TH1*)m.r1,(TH1*)m.r2,(TH1*)m.r3,
                                          (TH1*)m.r1c,(TH1*)m.r2c,(TH1*)m.r3c,
                                          m.e1,m.e2,m.e3,
                                          (TH1*)m.rp1,(TH1*)m.rp2,(TH1*)m.rp3 }) {
            double I = h->Integral();
            if (I>0) h->Scale(1.0/I);
        }

        // ---- Data: uncut ----
        {
            TCanvas c(Form("data_uncut_%s",label),Form("%s Data Uncut",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(d.r1,d.r2,d.r3);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                (i==0? d.r1 : i==1? d.r2 : d.r3)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/data_uncut_%s.png",label));
        }

        // ---- Data: circle overlay ----
        {
            TCanvas c(Form("data_circle_%s",label),Form("%s Data Circle",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(d.r1,d.r2,d.r3);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                TH2* hh = i==0? d.r1 : i==1? d.r2 : d.r3;
                double R = i==0? R1 : i==1? R2 : R3;
                hh->Draw("COLZ");
                TEllipse cir(0,0,R,R);
                cir.SetLineColor(kRed);
                cir.SetLineWidth(2);
                cir.SetFillStyle(0);
                cir.Draw("same");
            }
            c.SaveAs(Form("output/normalization/data_circle_%s.png",label));
        }

        // ---- Data: cut ----
        {
            TCanvas c(Form("data_cut_%s",label),Form("%s Data Cut",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(d.r1c,d.r2c,d.r3c);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLeftMargin(.15);
                gPad->SetRightMargin(.15);
                gPad->SetLogz();
                (i==0? d.r1c : i==1? d.r2c : d.r3c)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/data_cut_%s.png",label));
        }

        // ---- Data: edges ----
        {
            TCanvas c(Form("data_edges_%s",label),Form("%s Data Edges",label),1800,600);
            c.Divide(3,1);
            SetSame1DScale(d.e1,d.e2,d.e3);
            c.cd(1); d.e1->Draw("HIST");
            c.cd(2); d.e2->Draw("HIST");
            c.cd(3); d.e3->Draw("HIST");
            c.SaveAs(Form("output/normalization/data_edges_%s.png",label));
        }

        // ---- Data: rad-phi ----
        {
            TCanvas c(Form("data_radphi_%s",label),Form("%s Data ρ-φ",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(d.rp1,d.rp2,d.rp3);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? d.rp1 : i==1? d.rp2 : d.rp3)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/data_radphi_%s.png",label));
        }

        // ---- MC: uncut ----
        {
            TCanvas c(Form("mc_uncut_%s",label),Form("%s MC Uncut",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(m.r1,m.r2,m.r3);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? m.r1 : i==1? m.r2 : m.r3)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/mc_uncut_%s.png",label));
        }

        // ---- MC: circle overlay ----
        {
            TCanvas c(Form("mc_circle_%s",label),Form("%s MC Circle",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(m.r1,m.r2,m.r3);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                TH2* hh = i==0? m.r1 : i==1? m.r2 : m.r3;
                double R = i==0? R1 : i==1? R2 : R3;
                hh->Draw("COLZ");
                TEllipse cir(0,0,R,R);
                cir.SetLineColor(kRed);
                cir.SetLineWidth(2);
                cir.SetFillStyle(0);
                cir.Draw("same");
            }
            c.SaveAs(Form("output/normalization/mc_circle_%s.png",label));
        }

        // ---- MC: cut ----
        {
            TCanvas c(Form("mc_cut_%s",label),Form("%s MC Cut",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(m.r1c,m.r2c,m.r3c);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? m.r1c : i==1? m.r2c : m.r3c)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/mc_cut_%s.png",label));
        }

        // ---- MC: edges ----
        {
            TCanvas c(Form("mc_edges_%s",label),Form("%s MC Edges",label),1800,600);
            c.Divide(3,1);
            SetSame1DScale(m.e1,m.e2,m.e3);
            c.cd(1); m.e1->Draw("HIST");
            c.cd(2); m.e2->Draw("HIST");
            c.cd(3); m.e3->Draw("HIST");
            c.SaveAs(Form("output/normalization/mc_edges_%s.png",label));
        }

        // ---- MC: rad-phi ----
        {
            TCanvas c(Form("mc_radphi_%s",label),Form("%s MC ρ-φ",label),1800,600);
            c.Divide(3,1);
            SetSame2DScale(m.rp1,m.rp2,m.rp3);
            for (int i=0; i<3; ++i) {
                c.cd(i+1);
                gPad->SetLogz();
                (i==0? m.rp1 : i==1? m.rp2 : m.rp3)->Draw("COLZ");
            }
            c.SaveAs(Form("output/normalization/mc_radphi_%s.png",label));
        }
    }

    return 0;
}