#!/usr/bin/env python3
"""
First-pass photon-efficiency exclusivity study for CLAS12 RGA Fa18 inbending.

Purpose
-------
Define a clean e'p'gamma_tag -> gamma_probe denominator WITHOUT requiring a
reconstructed probe photon.  This version deliberately stops before extracting
an efficiency.  It shows how each exclusivity requirement changes the predicted
probe phase space, especially the FT/high-energy region.

Default cumulative stages
-------------------------
0 baseline : W>2, standard proton, tag beta+fiducial
1 Mx2(ep)  : -0.231 < Mx2(ep) < 0.309 GeV^2
2 Mx2(eg)  : Mx2(e gamma_tag) > 1.4 GeV^2
3 coplanar : |Delta phi_copl| < 5.7 deg (Trento-plane residual)

angle(gamma_tag,X) and Delta t=t_p-t_gamma are DIAGNOSTICS ONLY in v1.
No reconstructed-probe requirement is made anywhere.

Performance
-----------
Uses ROOT RDataFrame, enables up to 8 ROOT worker threads, disables unused
branches implicitly through RDataFrame column access, books all actions before
execution, and runs the sample graphs together with ROOT.RDF.RunGraphs.
"""

import argparse
import glob
import os
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

BASE = "/work/clas12/thayward/photon_efficiency/ROOT_trees"
TREE = "PhotonEfficiency"

SAMPLES = {
    "data":    ("Data",       "data/fa18_inb/*.root",    ROOT.kBlack),
    "aaogen":  ("AAOgen",     "aaogen/fa18_inb/*.root",  ROOT.kRed + 1),
    "clasdis": ("CLASDIS",    "clasdis/fa18_inb/*.root", ROOT.kBlue + 1),
    "dvcsgen": ("DVCSgen",    "dvcsgen/fa18_inb/*.root", ROOT.kGreen + 2),
}

STAGES = [
    ("baseline", "Baseline", "baseline"),
    ("mx2ep", " + M_{X}^{2}(ep)", "cut_mx2ep"),
    ("mx2eg", " + M_{X}^{2}(e#gamma)", "cut_mx2eg"),
    ("coplanarity", " + coplanarity", "cut_coplanarity"),
]

ROOT.gInterpreter.Declare(r'''
#include <cmath>
#include <algorithm>
#include "Math/Vector4D.h"
#include "TVector3.h"

static double pe2_wrap180(double x) {
    while (x <= -180.0) x += 360.0;
    while (x >   180.0) x -= 360.0;
    return x;
}

static TVector3 pe2_vec(double p, double th_deg, double ph_deg) {
    const double d=M_PI/180.0, th=th_deg*d, ph=ph_deg*d;
    return TVector3(p*std::sin(th)*std::cos(ph),
                    p*std::sin(th)*std::sin(ph),
                    p*std::cos(th));
}

double pe2_angle(double th1,double ph1,double th2,double ph2) {
    TVector3 a=pe2_vec(1.0,th1,ph1), b=pe2_vec(1.0,th2,ph2);
    if (a.Mag2()<=0 || b.Mag2()<=0) return -999.0;
    double c=std::max(-1.0,std::min(1.0,a.Unit().Dot(b.Unit())));
    return std::acos(c)*180.0/M_PI;
}

double pe2_mx2_eg(double Eb,double ep,double eth,double eph,
                  double gp,double gth,double gph) {
    const double me=0.00051099895, mp=0.9382720813;
    const double pb=std::sqrt(std::max(0.0,Eb*Eb-me*me));
    TVector3 e=pe2_vec(ep,eth,eph), g=pe2_vec(gp,gth,gph);
    const double Ee=std::sqrt(ep*ep+me*me);
    const double E=Eb+mp-Ee-gp;
    const double px=-e.X()-g.X(), py=-e.Y()-g.Y(), pz=pb-e.Z()-g.Z();
    return E*E-px*px-py*py-pz*pz;
}

double pe2_trento_phi(double Eb,double ep,double eth,double eph,
                      double hp,double hth,double hph) {
    const double me=0.00051099895;
    const double pb=std::sqrt(std::max(0.0,Eb*Eb-me*me));
    TVector3 k(0,0,pb), kp=pe2_vec(ep,eth,eph), h=pe2_vec(hp,hth,hph);
    TVector3 q=k-kp;
    if (q.Mag2()<=1e-18) return -999.0;
    TVector3 qhat=q.Unit(), nl=q.Cross(k), nh=q.Cross(h);
    if (nl.Mag2()<=1e-18 || nh.Mag2()<=1e-18) return -999.0;
    nl=nl.Unit(); nh=nh.Unit();
    double c=std::max(-1.0,std::min(1.0,nl.Dot(nh)));
    double s=qhat.Dot(nl.Cross(nh));
    return pe2_wrap180(std::atan2(s,c)*180.0/M_PI);
}

double pe2_copl(double Eb,double ep,double eth,double eph,
                double pp,double pth,double pph,
                double gp,double gth,double gph) {
    double a=pe2_trento_phi(Eb,ep,eth,eph,pp,pth,pph);
    double b=pe2_trento_phi(Eb,ep,eth,eph,gp,gth,gph);
    if (a < -900 || b < -900) return -999.0;
    double d=pe2_wrap180(a-b);
    double sign=(d>=0.0)?1.0:-1.0;
    return sign*(180.0-std::fabs(d));
}

double pe2_delta_t(double Eb,double ep,double eth,double eph,
                   double pp,double pth,double pph,
                   double gp,double gth,double gph) {
    const double me=0.00051099895, mp=0.9382720813;
    const double pb=std::sqrt(std::max(0.0,Eb*Eb-me*me));
    ROOT::Math::PxPyPzEVector k(0,0,pb,Eb), target(0,0,0,mp);
    TVector3 ev=pe2_vec(ep,eth,eph), pv=pe2_vec(pp,pth,pph), gv=pe2_vec(gp,gth,gph);
    ROOT::Math::PxPyPzEVector kp(ev.X(),ev.Y(),ev.Z(),std::sqrt(ep*ep+me*me));
    ROOT::Math::PxPyPzEVector pr(pv.X(),pv.Y(),pv.Z(),std::sqrt(pp*pp+mp*mp));
    ROOT::Math::PxPyPzEVector ga(gv.X(),gv.Y(),gv.Z(),gp);
    auto q=k-kp;
    return (target-pr).M2()-(q-ga).M2();
}
''')


def files_for(pattern):
    return sorted(glob.glob(os.path.join(BASE, pattern)))


def make_rdf(files):
    v=ROOT.std.vector("string")()
    for f in files:
        v.push_back(f)
    return ROOT.RDataFrame(TREE, v)


def define_columns(df):
    cols={str(x) for x in df.GetColumnNames()}
    required={
        "W","p_pass_standard","tag_pass_beta","tag_pass_fiducial",
        "Mx_ep","probe_corr_p","probe_corr_theta","probe_corr_phi",
        "beam_energy","e_p","e_theta","e_phi",
        "p_corr_p","p_corr_theta","p_corr_phi",
        "tag_corr_p","tag_corr_theta","tag_corr_phi"
    }
    missing=sorted(required-cols)
    if missing:
        raise RuntimeError("Missing required branches: "+", ".join(missing))

    mx2ep="Mx2_ep" if "Mx2_ep" in cols else "(Mx_ep>=0 ? Mx_ep*Mx_ep : -Mx_ep*Mx_ep)"
    d=(df.Define("v_mx2ep",mx2ep)
         .Define("v_mx2eg","pe2_mx2_eg(beam_energy,e_p,e_theta,e_phi,tag_corr_p,tag_corr_theta,tag_corr_phi)")
         .Define("v_copl","pe2_copl(beam_energy,e_p,e_theta,e_phi,p_corr_p,p_corr_theta,p_corr_phi,tag_corr_p,tag_corr_theta,tag_corr_phi)")
         .Define("v_angle_gX","pe2_angle(tag_corr_theta,tag_corr_phi,probe_corr_theta,probe_corr_phi)")
         .Define("v_dt","pe2_delta_t(beam_energy,e_p,e_theta,e_phi,p_corr_p,p_corr_theta,p_corr_phi,tag_corr_p,tag_corr_theta,tag_corr_phi)")
         .Define("probe_region","probe_corr_theta<=5.5 ? 0 : (probe_corr_theta<=36.0 ? 1 : 2)"))

    # Baseline uses only reconstructed e,p,tag quality; never the reconstructed probe.
    d=d.Define("baseline","W>2.0 && p_pass_standard==1 && tag_pass_beta==1 && tag_pass_fiducial==1 && probe_corr_p>0 && probe_corr_theta>=0")
    d=d.Define("cut_mx2ep","baseline && v_mx2ep>-0.231 && v_mx2ep<0.309")
    d=d.Define("cut_mx2eg","cut_mx2ep && v_mx2eg>1.4")
    d=d.Define("cut_coplanarity","cut_mx2eg && fabs(v_copl)<5.7")
    return d


def book_sample(key, df):
    out={"counts":{},"h1":{},"h2":{}}
    for ik,(stage,label,cut) in enumerate(STAGES):
        node=df.Filter(cut, stage)
        out["counts"][stage]=node.Count()
        for reg,rsel in [("all","true"),("FT","probe_region==0"),("FD","probe_region==1")]:
            n=node.Filter(rsel)
            out["h1"][(stage,reg,"E")]=n.Histo1D((f"hE_{key}_{stage}_{reg}","",106,0,10.6),"probe_corr_p")
            out["h1"][(stage,reg,"theta")]=n.Histo1D((f"hth_{key}_{stage}_{reg}","",90,0,45),"probe_corr_theta")

    final=df.Filter("cut_coplanarity")
    for reg,rsel in [("all","true"),("FT","probe_region==0"),("FD","probe_region==1")]:
        n=final.Filter(rsel)
        out["h1"][("final",reg,"angle")]=n.Histo1D((f"hang_{key}_{reg}","",120,0,30),"v_angle_gX")
        out["h1"][("final",reg,"dt")]=n.Histo1D((f"hdt_{key}_{reg}","",150,-1,4),"v_dt")
        out["h2"][(reg,"angle")]=n.Histo2D((f"h2ang_{key}_{reg}","",106,0,10.6,120,0,30),"probe_corr_p","v_angle_gX")
        out["h2"][(reg,"dt")]=n.Histo2D((f"h2dt_{key}_{reg}","",106,0,10.6,150,-1,4),"probe_corr_p","v_dt")
    return out


def clone(h):
    x=h.GetValue().Clone()
    x.SetDirectory(0)
    return x


def unit(h):
    x=h.Clone(h.GetName()+"_unit")
    x.SetDirectory(0)
    integ=x.Integral("width")
    if integ>0: x.Scale(1.0/integ)
    return x


def style(h,key):
    h.SetLineColor(SAMPLES[key][2]); h.SetMarkerColor(SAMPLES[key][2]); h.SetLineWidth(2)
    if key=="data":
        h.SetMarkerStyle(20); h.SetMarkerSize(0.55)



def prep_hist(results,key,stage,region,var,normalize=True):
    h=clone(results[key]["h1"][(stage,region,var)])
    style(h,key)
    if normalize:
        integ=h.Integral("width")
        if integ>0: h.Scale(1.0/integ)
    return h

def setup_pad(p, right=0.04):
    p.SetLeftMargin(0.14); p.SetRightMargin(right); p.SetBottomMargin(0.14); p.SetTopMargin(0.10)
    p.SetTicks(1,1)

def draw_figure1(results,outdir):
    c=ROOT.TCanvas("c_fig1","",1500,1050); c.Divide(2,2)
    labels=["Baseline","M_{X}^{2}(ep)","M_{X}^{2}(e#gamma)","Coplanarity"]
    keep=[]
    for ipad,key in enumerate(SAMPLES,1):
        p=c.cd(ipad); setup_pad(p)
        vals=[int(results[key]["counts"][st].GetValue()) for st,_,_ in STAGES]
        den=max(1,vals[0])
        g=ROOT.TGraph(len(vals)); keep.append(g)
        for i,v in enumerate(vals): g.SetPoint(i,i,v/den)
        g.SetLineWidth(3); g.SetMarkerStyle(20); g.SetMarkerSize(1.2)
        g.SetLineColor(SAMPLES[key][2]); g.SetMarkerColor(SAMPLES[key][2])
        g.SetTitle(""); g.Draw("ALP")
        g.GetYaxis().SetTitle("Fraction of baseline hypotheses retained")
        g.GetYaxis().SetRangeUser(0,1.08); g.GetXaxis().SetLimits(-0.25,3.25)
        g.GetXaxis().SetNdivisions(4,False)
        for i,lab in enumerate(labels): g.GetXaxis().ChangeLabel(i+1,-1,0.035,-1,-1,-1,lab)
        tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.047)
        tx.DrawLatex(0.17,0.92,SAMPLES[key][0])
        tx.SetTextSize(0.035); tx.DrawLatex(0.17,0.85,f"Baseline N = {vals[0]:,}")
    c.SaveAs(str(Path(outdir)/"01_exclusivity_cutflow.png"))

def draw_stage_overlay_pad(p,results,key,region,title):
    setup_pad(p)
    stage_cols=[ROOT.kGray+2,ROOT.kBlue+1,ROOT.kOrange+7,ROOT.kRed+1]
    hs=[]; ymax=0.0
    for i,(st,lab,cut) in enumerate(STAGES):
        h=clone(results[key]["h1"][(st,region,"E")])
        integ=h.Integral("width")
        if integ>0: h.Scale(1.0/integ)
        h.SetLineColor(stage_cols[i]); h.SetLineWidth(3)
        hs.append(h); ymax=max(ymax,h.GetMaximum())
    for i,h in enumerate(hs):
        h.SetTitle(""); h.GetXaxis().SetTitle("Predicted probe momentum p_{#gamma,probe} (GeV)")
        h.GetYaxis().SetTitle("Unit-normalized density"); h.GetYaxis().SetRangeUser(0,1.25*ymax)
        h.Draw("HIST" if i==0 else "HIST SAME")
    leg=ROOT.TLegend(0.48,0.62,0.91,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.033)
    for h,(_,lab,_) in zip(hs,STAGES): leg.AddEntry(h,lab.strip().lstrip("+").strip(),"l")
    leg.Draw()
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.044); tx.DrawLatex(0.17,0.92,title)
    return hs,leg

def draw_figure2(results,outdir):
    c=ROOT.TCanvas("c_fig2","",1500,1050); c.Divide(2,2); keep=[]
    keep += draw_stage_overlay_pad(c.cd(1),results,"data","FT","DATA — FT predicted probe")
    keep += draw_stage_overlay_pad(c.cd(2),results,"aaogen","FT","AAOgen — FT predicted probe")
    keep += draw_stage_overlay_pad(c.cd(3),results,"data","FD","DATA — FD predicted probe")
    keep += draw_stage_overlay_pad(c.cd(4),results,"aaogen","FD","AAOgen — FD predicted probe")
    c.SaveAs(str(Path(outdir)/"02_probe_momentum_survival.png"))

def draw_final_overlay_pad(p,results,region,var,xlabel,title):
    setup_pad(p); hs=[]; ymax=0
    for key in SAMPLES:
        h=prep_hist(results,key,"coplanarity",region,var,True); hs.append((key,h)); ymax=max(ymax,h.GetMaximum())
    for i,(key,h) in enumerate(hs):
        h.SetTitle(""); h.GetXaxis().SetTitle(xlabel); h.GetYaxis().SetTitle("Unit-normalized density")
        h.GetYaxis().SetRangeUser(0,1.25*ymax)
        opt=("E1" if key=="data" else "HIST")+("" if i==0 else " SAME"); h.Draw(opt)
    leg=ROOT.TLegend(0.60,0.64,0.91,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.034)
    for key,h in hs: leg.AddEntry(h,SAMPLES[key][0],"lep" if key=="data" else "l")
    leg.Draw()
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.044); tx.DrawLatex(0.17,0.92,title)
    return hs,leg

def draw_figure3(results,outdir):
    c=ROOT.TCanvas("c_fig3","",1500,1050); c.Divide(2,2); keep=[]
    keep += draw_final_overlay_pad(c.cd(1),results,"FT","E","Predicted probe momentum p_{#gamma,probe} (GeV)","Final sample — FT momentum")
    keep += draw_final_overlay_pad(c.cd(2),results,"FD","E","Predicted probe momentum p_{#gamma,probe} (GeV)","Final sample — FD momentum")
    keep += draw_final_overlay_pad(c.cd(3),results,"FT","theta","Predicted probe #theta_{#gamma,probe} (deg)","Final sample — FT angle")
    keep += draw_final_overlay_pad(c.cd(4),results,"FD","theta","Predicted probe #theta_{#gamma,probe} (deg)","Final sample — FD angle")
    c.SaveAs(str(Path(outdir)/"03_final_sample_composition.png"))

def draw_2d_pad(p,results,key,var,title):
    setup_pad(p,0.16)
    h=clone(results[key]["h2"][("FT",var)])
    h.SetTitle(""); h.GetXaxis().SetTitle("Predicted probe momentum p_{#gamma,probe} (GeV)")
    if var=="angle": h.GetYaxis().SetTitle("#theta(#gamma_{tag},X) (deg)")
    else: h.GetYaxis().SetTitle("#Delta t=t_{p}-t_{#gamma} (GeV^{2})")
    h.Draw("COLZ")
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043); tx.DrawLatex(0.17,0.92,title)
    return h

def draw_figure4(results,outdir):
    c=ROOT.TCanvas("c_fig4","",1500,1050); c.Divide(2,2); keep=[]
    keep.append(draw_2d_pad(c.cd(1),results,"data","angle","DATA — FT"))
    keep.append(draw_2d_pad(c.cd(2),results,"aaogen","angle","AAOgen — FT"))
    keep.append(draw_2d_pad(c.cd(3),results,"data","dt","DATA — FT"))
    keep.append(draw_2d_pad(c.cd(4),results,"aaogen","dt","AAOgen — FT"))
    c.SaveAs(str(Path(outdir)/"04_next_cut_optimization.png"))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=8,help="ROOT worker threads (1-8; default 8)")
    ap.add_argument("--output",default="output_RGA_study")
    args=ap.parse_args()
    workers=max(1,min(8,args.workers))
    ROOT.EnableImplicitMT(workers)
    Path(args.output).mkdir(parents=True,exist_ok=True)

    print("="*72)
    print("Photon-efficiency RGA study — exclusivity optimization")
    print(f"ROOT worker threads: {workers}")
    print("No reconstructed-probe requirement; angle(gamma,X) and Delta-t are QA only.")
    print("="*72)

    results={}; all_actions=[]
    for key,(label,pat,color) in SAMPLES.items():
        fs=files_for(pat)
        if not fs: raise RuntimeError(f"No ROOT files for {key}: {pat}")
        print(f"{label:10s}: {len(fs):3d} ROOT files")
        df=define_columns(make_rdf(fs))
        results[key]=book_sample(key,df)
        all_actions += list(results[key]["counts"].values())
        all_actions += list(results[key]["h1"].values())
        all_actions += list(results[key]["h2"].values())

    print(f"Executing {len(all_actions)} booked actions ...")
    ROOT.RDF.RunGraphs(all_actions)

    with open(Path(args.output)/"cutflow.txt","w") as f:
        hdr="sample"+"".join(f" {st[0]:>16s}" for st in STAGES)
        print(hdr); f.write(hdr+"\n")
        for key in SAMPLES:
            vals=[int(results[key]["counts"][st[0]].GetValue()) for st in STAGES]
            line=f"{key:10s}"+"".join(f" {v:16d}" for v in vals)
            print(line); f.write(line+"\n")

    draw_figure1(results,args.output)
    draw_figure2(results,args.output)
    draw_figure3(results,args.output)
    draw_figure4(results,args.output)
    print(f"Done. Four summary figures + cutflow.txt written to: {args.output}")

if __name__=="__main__":
    main()
