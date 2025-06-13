/*****************************************************************************
 * File: unnormalized_drift_chambers.cpp
 *
 * Compile with:
 *   g++ unnormalized_drift_chambers.cpp -o unnormalized_drift_chambers \
 *       $(root-config --cflags --libs)
 *
 * Run with:
 *   ./unnormalized_drift_chambers [NeventsSIDISDVCS] [dataFile] [mcFile]
 *
 * This version:
 *   • Fills and normalizes data and MC histograms separately.
 *   • Does NOT compute any data/MC ratios.
 *   • Saves Data and MC canvases in output/normalization/.
 *****************************************************************************/

#include <iostream>
#include <string>
#include <cmath>
#include "TChain.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TPad.h"
#include "TEllipse.h"

// Helper: unify 2D color scale
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

// Helper: unify 1D y-scale
void SetSame1DScale(TH1* a, TH1* b, TH1* c) {
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
    // 1) Parse arguments
    Long64_t maxEvents = -1;
    if (argc > 1) {
        maxEvents = std::stoll(argv[1]);
        if (maxEvents == 0) maxEvents = -1;
    }
    bool useDataFile = (argc > 2);
    std::string dataFile = useDataFile ? argv[2] : "";
    bool useMCFile = (argc > 3);
    std::string mcFile = useMCFile ? argv[3] : "";

    // 2) Setup TChains
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

    // 3) Branch addresses
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
    mcCh  .SetBranchAddress("particle_pid",&pid);
    mcCh  .SetBranchAddress("traj_x_6",&x6);
    mcCh  .SetBranchAddress("traj_y_6",&y6);
    mcCh  .SetBranchAddress("traj_x_18",&x18);
    mcCh  .SetBranchAddress("traj_y_18",&y18);
    mcCh  .SetBranchAddress("traj_x_36",&x36);
    mcCh  .SetBranchAddress("traj_y_36",&y36);
    mcCh  .SetBranchAddress("traj_edge_6",&e6);
    mcCh  .SetBranchAddress("traj_edge_18",&e18);
    mcCh  .SetBranchAddress("traj_edge_36",&e36);

    // 4) Create histograms
    const int NB2 = 200;
    TH2D *d_r1 = new TH2D("d_r1","Data R1;x;y",NB2,-180,180,NB2,-180,180),
         *d_r2 = new TH2D("d_r2","Data R2;x;y",NB2,-280,280,NB2,-280,280),
         *d_r3 = new TH2D("d_r3","Data R3;x;y",NB2,-450,450,NB2,-450,450);
    TH2D *m_r1 = new TH2D("m_r1","MC   R1;x;y",NB2,-180,180,NB2,-180,180),
         *m_r2 = new TH2D("m_r2","MC   R2;x;y",NB2,-280,280,NB2,-280,280),
         *m_r3 = new TH2D("m_r3","MC   R3;x;y",NB2,-450,450,NB2,-450,450);
    TH2D *d_r1c = new TH2D("d_r1c","Data R1 cut;x;y",NB2,-180,180,NB2,-180,180),
         *d_r2c = new TH2D("d_r2c","Data R2 cut;x;y",NB2,-280,280,NB2,-280,280),
         *d_r3c = new TH2D("d_r3c","Data R3 cut;x;y",NB2,-450,450,NB2,-450,450);
    TH2D *m_r1c = new TH2D("m_r1c","MC   R1 cut;x;y",NB2,-180,180,NB2,-180,180),
         *m_r2c = new TH2D("m_r2c","MC   R2 cut;x;y",NB2,-280,280,NB2,-280,280),
         *m_r3c = new TH2D("m_r3c","MC   R3 cut;x;y",NB2,-450,450,NB2,-450,450);
    const double R1=140.0, R2=215.0, R3=290.0;

    const int NB1 = 125;
    TH1D *d_e1 = new TH1D("d_e1","Data edge R1;cm",NB1,0,125),
         *d_e2 = new TH1D("d_e2","Data edge R2;cm",NB1,0,125),
         *d_e3 = new TH1D("d_e3","Data edge R3;cm",NB1,0,125);
    TH1D *m_e1 = new TH1D("m_e1","MC   edge R1;cm",NB1,0,125),
         *m_e2 = new TH1D("m_e2","MC   edge R2;cm",NB1,0,125),
         *m_e3 = new TH1D("m_e3","MC   edge R3;cm",NB1,0,125);

    // rad-phi
    int nbR1=150, nbR2=230, nbR3=400, nbPhi=180;
    TH2D *d_rp1 = new TH2D("d_rp1","Data R1 ρ-φ;ρ;φ",nbR1,20,170,nbPhi,0,360),
         *d_rp2 = new TH2D("d_rp2","Data R2 ρ-φ;ρ;φ",nbR2,20,250,nbPhi,0,360),
         *d_rp3 = new TH2D("d_rp3","Data R3 ρ-φ;ρ;φ",nbR3,20,420,nbPhi,0,360);
    TH2D *m_rp1 = new TH2D("m_rp1","MC   R1 ρ-φ;ρ;φ",nbR1,20,170,nbPhi,0,360),
         *m_rp2 = new TH2D("m_rp2","MC   R2 ρ-φ;ρ;φ",nbR2,20,250,nbPhi,0,360),
         *m_rp3 = new TH2D("m_rp3","MC   R3 ρ-φ;ρ;φ",nbR3,20,420,nbPhi,0,360);

    // 5) Fill Data
    Long64_t nD = dataCh.GetEntries();
    if (maxEvents>0 && maxEvents<nD) nD = maxEvents;
    for (Long64_t i=0; i<nD; ++i) {
        dataCh.GetEntry(i);
        if (pid!=11) continue;
        if (e6<=3 || e18<=3 || e36<=10) continue;
        if (x6!=-9999) {
            d_r1->Fill(x6,y6);
            double r = std::hypot(x6,y6);
            if (r<R1) d_r1c->Fill(x6,y6);
            double φ = std::atan2(y6,x6)-M_PI;
            if (φ<0) φ+=2*M_PI;
            d_rp1->Fill(r,φ*180.0/M_PI);
        }
        if (x18!=-9999) {
            d_r2->Fill(x18,y18);
            double r = std::hypot(x18,y18);
            if (r<R2) d_r2c->Fill(x18,y18);
            double φ = std::atan2(y18,x18)-M_PI;
            if (φ<0) φ+=2*M_PI;
            d_rp2->Fill(r,φ*180.0/M_PI);
        }
        if (x36!=-9999) {
            d_r3->Fill(x36,y36);
            double r = std::hypot(x36,y36);
            if (r<R3) d_r3c->Fill(x36,y36);
            double φ = std::atan2(y36,x36)-M_PI;
            if (φ<0) φ+=2*M_PI;
            d_rp3->Fill(r,φ*180.0/M_PI);
        }
        d_e1->Fill(e6);
        d_e2->Fill(e18);
        d_e3->Fill(e36);
    }

    // 6) Fill MC
    Long64_t nM = mcCh.GetEntries();
    if (maxEvents>0 && maxEvents<nM) nM = maxEvents;
    for (Long64_t i=0; i<nM; ++i) {
        mcCh.GetEntry(i);
        if (pid!=11) continue;
        if (e6<=3 || e18<=3 || e36<=10) continue;
        if (x6!=-9999) {
            m_r1->Fill(x6,y6);
            double r = std::hypot(x6,y6);
            if (r<R1) m_r1c->Fill(x6,y6);
            double φ = std::atan2(y6,x6)-M_PI;
            if (φ<0) φ+=2*M_PI;
            m_rp1->Fill(r,φ*180.0/M_PI);
        }
        if (x18!=-9999) {
            m_r2->Fill(x18,y18);
            double r = std::hypot(x18,y18);
            if (r<R2) m_r2c->Fill(x18,y18);
            double φ = std::atan2(y18,x18)-M_PI;
            if (φ<0) φ+=2*M_PI;
            m_rp2->Fill(r,φ*180.0/M_PI);
        }
        if (x36!=-9999) {
            m_r3->Fill(x36,y36);
            double r = std::hypot(x36,y36);
            if (r<R3) m_r3c->Fill(x36,y36);
            double φ = std::atan2(y36,x36)-M_PI;
            if (φ<0) φ+=2*M_PI;
            m_rp3->Fill(r,φ*180.0/M_PI);
        }
        m_e1->Fill(e6);
        m_e2->Fill(e18);
        m_e3->Fill(e36);
    }

    // 7) Normalize histograms using TH1* arrays
    TH1* dataHists[] = {
        d_r1, d_r2, d_r3,
        d_r1c, d_r2c, d_r3c,
        d_e1, d_e2, d_e3,
        d_rp1, d_rp2, d_rp3
    };
    for (TH1* h : dataHists) {
        double I = h->Integral();
        if (I>0) h->Scale(1.0/I);
    }
    TH1* mcHists[] = {
        m_r1, m_r2, m_r3,
        m_r1c, m_r2c, m_r3c,
        m_e1, m_e2, m_e3,
        m_rp1, m_rp2, m_rp3
    };
    for (TH1* h : mcHists) {
        double I = h->Integral();
        if (I>0) h->Scale(1.0/I);
    }

    gStyle->SetOptStat(0);

    // 8) Draw & save Data canvases
    {
        TCanvas c("data_uncut","Data Uncut",1800,600);
        c.Divide(3,1);
        SetSame2DScale(d_r1,d_r2,d_r3);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.15);
            gPad->SetLogz();
            (i==0? d_r1 : i==1? d_r2 : d_r3)->Draw("COLZ");
        }
        c.SaveAs("output/normalization/data_uncut.png");
    }
    {
        TCanvas c("data_circle","Data Circle Overlay",1800,600);
        c.Divide(3,1);
        SetSame2DScale(d_r1,d_r2,d_r3);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15); gPad->SetLogz();
            TH2* h = (i==0? d_r1 : i==1? d_r2 : d_r3);
            double R = (i==0? R1 : i==1? R2 : R3);
            h->Draw("COLZ");
            TEllipse cir(0,0,R,R);
            cir.SetLineColor(kRed); cir.SetLineWidth(2); cir.SetFillStyle(0);
            cir.Draw("same");
        }
        c.SaveAs("output/normalization/data_circle.png");
    }
    {
        TCanvas c("data_cut","Data Cut",1800,600);
        c.Divide(3,1);
        SetSame2DScale(d_r1c,d_r2c,d_r3c);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.15); gPad->SetLogz();
            (i==0? d_r1c : i==1? d_r2c : d_r3c)->Draw("COLZ");
        }
        c.SaveAs("output/normalization/data_cut.png");
    }
    {
        TCanvas c("data_edges","Data Edges",1800,600);
        c.Divide(3,1);
        SetSame1DScale(d_e1,d_e2,d_e3);
        d_e1->Draw("HIST"); c.cd(2); d_e2->Draw("HIST"); c.cd(3); d_e3->Draw("HIST");
        c.SaveAs("output/normalization/data_edges.png");
    }
    {
        TCanvas c("data_radphi","Data ρ-φ",1800,600);
        c.Divide(3,1);
        SetSame2DScale(d_rp1,d_rp2,d_rp3);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLogz();
            (i==0? d_rp1 : i==1? d_rp2 : d_rp3)->Draw("COLZ");
        }
        c.SaveAs("output/normalization/data_radphi.png");
    }

    // 9) Draw & save MC canvases (mirror)
    {
        TCanvas c("mc_uncut","MC Uncut",1800,600);
        c.Divide(3,1);
        SetSame2DScale(m_r1,m_r2,m_r3);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLogz();
            (i==0? m_r1 : i==1? m_r2 : m_r3)->Draw("COLZ");
        }
        c.SaveAs("output/normalization/mc_uncut.png");
    }
    {
        TCanvas c("mc_circle","MC Circle Overlay",1800,600);
        c.Divide(3,1);
        SetSame2DScale(m_r1,m_r2,m_r3);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLogz();
            TH2* h = (i==0? m_r1 : i==1? m_r2 : m_r3);
            double R = (i==0? R1 : i==1? R2 : R3);
            h->Draw("COLZ");
            TEllipse cir(0,0,R,R);
            cir.SetLineColor(kRed); cir.SetLineWidth(2); cir.SetFillStyle(0);
            cir.Draw("same");
        }
        c.SaveAs("output/normalization/mc_circle.png");
    }
    {
        TCanvas c("mc_cut","MC Cut",1800,600);
        c.Divide(3,1);
        SetSame2DScale(m_r1c,m_r2c,m_r3c);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLogz();
            (i==0? m_r1c : i==1? m_r2c : m_r3c)->Draw("COLZ");
        }
        c.SaveAs("output/normalization/mc_cut.png");
    }
    {
        TCanvas c("mc_edges","MC Edges",1800,600);
        c.Divide(3,1);
        SetSame1DScale(m_e1,m_e2,m_e3);
        m_e1->Draw("HIST"); c.cd(2); m_e2->Draw("HIST"); c.cd(3); m_e3->Draw("HIST");
        c.SaveAs("output/normalization/mc_edges.png");
    }
    {
        TCanvas c("mc_radphi","MC ρ-φ",1800,600);
        c.Divide(3,1);
        SetSame2DScale(m_rp1,m_rp2,m_rp3);
        for (int i=0; i<3; ++i) {
            c.cd(i+1);
            gPad->SetLogz();
            (i==0? m_rp1 : i==1? m_rp2 : m_rp3)->Draw("COLZ");
        }
        c.SaveAs("output/normalization/mc_radphi.png");
    }

    return 0;
}