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


def draw_overlay(results, outdir, stage, region, var, xlabel):
    c=ROOT.TCanvas(f"c_{stage}_{region}_{var}","",950,760)
    c.SetLeftMargin(0.13); c.SetRightMargin(0.04); c.SetBottomMargin(0.13); c.SetTopMargin(0.08)
    hs=[]
    ymax=0
    for key in SAMPLES:
        h=unit(clone(results[key]["h1"][(stage,region,var)])); style(h,key); hs.append((key,h)); ymax=max(ymax,h.GetMaximum())
    first=True
    for key,h in hs:
        h.SetTitle("")
        h.GetXaxis().SetTitle(xlabel); h.GetYaxis().SetTitle("Unit-normalized density")
        h.GetYaxis().SetRangeUser(0,max(1e-12,1.28*ymax))
        opt="E1" if key=="data" else "HIST"
        if not first: opt += " SAME"
        h.Draw(opt); first=False
    leg=ROOT.TLegend(0.58,0.67,0.91,0.89); leg.SetBorderSize(0); leg.SetFillStyle(0)
    for key,h in hs: leg.AddEntry(h,SAMPLES[key][0],"lep" if key=="data" else "l")
    leg.Draw()
    txt=ROOT.TLatex(); txt.SetNDC(); txt.SetTextSize(0.037); txt.DrawLatex(0.15,0.94,f"{region} predicted probe: {dict((a,b) for a,b,_ in STAGES)[stage]}")
    c.SaveAs(str(Path(outdir)/f"probe_{var}_{stage}_{region}.png"))


def draw_2d(results,outdir,key,region,var,ylabel):
    h=clone(results[key]["h2"][(region,var)])
    c=ROOT.TCanvas(f"c2_{key}_{region}_{var}","",980,760)
    c.SetLeftMargin(0.13); c.SetRightMargin(0.16); c.SetBottomMargin(0.13); c.SetTopMargin(0.08)
    h.SetTitle(""); h.GetXaxis().SetTitle("Predicted probe momentum p_{#gamma,probe} (GeV)"); h.GetYaxis().SetTitle(ylabel)
    h.Draw("COLZ")
    txt=ROOT.TLatex(); txt.SetNDC(); txt.SetTextSize(0.037); txt.DrawLatex(0.15,0.94,f"{SAMPLES[key][0]}: after M_{{X}}^{{2}}(ep), M_{{X}}^{{2}}(e#gamma), coplanarity")
    c.SaveAs(str(Path(outdir)/f"{var}_vs_probe_p_{key}_{region}.png"))


def draw_survival(results,outdir,region):
    c=ROOT.TCanvas(f"csurv_{region}","",950,760)
    c.SetLeftMargin(0.13); c.SetRightMargin(0.04); c.SetBottomMargin(0.13); c.SetTopMargin(0.08)
    mg=ROOT.TMultiGraph(); graphs=[]
    for key in SAMPLES:
        vals=[]
        for stage,_,_ in STAGES:
            h=clone(results[key]["h1"][(stage,region,"E")]); vals.append(h.Integral())
        den=vals[0] if vals[0]>0 else 1
        g=ROOT.TGraph(len(vals)); g.SetName(f"g_{key}_{region}")
        for i,v in enumerate(vals): g.SetPoint(i,i,v/den)
        g.SetLineColor(SAMPLES[key][2]); g.SetMarkerColor(SAMPLES[key][2]); g.SetLineWidth(2); g.SetMarkerStyle(20)
        mg.Add(g,"LP"); graphs.append((key,g))
    mg.Draw("A"); mg.SetTitle(""); mg.GetYaxis().SetTitle("Fraction of baseline hypotheses retained"); mg.GetXaxis().SetTitle("")
    mg.GetYaxis().SetRangeUser(0,1.08); mg.GetXaxis().SetLimits(-0.25,len(STAGES)-0.75)
    ax=mg.GetXaxis(); ax.SetNdivisions(len(STAGES),False)
    labels=["Baseline","M_{X}^{2}(ep)","M_{X}^{2}(e#gamma)","Coplanarity"]
    for i,s in enumerate(labels): ax.ChangeLabel(i+1,-1,-1,-1,-1,-1,s)
    leg=ROOT.TLegend(0.63,0.68,0.91,0.89); leg.SetBorderSize(0); leg.SetFillStyle(0)
    for key,g in graphs: leg.AddEntry(g,SAMPLES[key][0],"lp")
    leg.Draw();
    txt=ROOT.TLatex(); txt.SetNDC(); txt.SetTextSize(0.04); txt.DrawLatex(0.15,0.94,f"Exclusivity cut survival: {region} predicted probe")
    c.SaveAs(str(Path(outdir)/f"cut_survival_{region}.png"))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=8,help="ROOT worker threads (1-8; default 8)")
    ap.add_argument("--output",default="output/photon_efficiency_exclusivity_study")
    args=ap.parse_args()
    workers=max(1,min(8,args.workers))
    ROOT.EnableImplicitMT(workers)
    Path(args.output).mkdir(parents=True,exist_ok=True)

    print("="*72)
    print("Photon-efficiency exclusivity study v1")
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
        hdr="sample"+"".join(f" {s[0]:>16s}" for s in STAGES)
        print(hdr); f.write(hdr+"\n")
        for key in SAMPLES:
            vals=[int(results[key]["counts"][s[0]].GetValue()) for s in STAGES]
            line=f"{key:10s}"+"".join(f" {v:16d}" for v in vals)
            print(line); f.write(line+"\n")

    for stage,_,_ in STAGES:
        for region in ("all","FT","FD"):
            draw_overlay(results,args.output,stage,region,"E","Predicted probe momentum p_{#gamma,probe} (GeV)")
            draw_overlay(results,args.output,stage,region,"theta","Predicted probe #theta_{#gamma,probe} (deg)")
    for region in ("all","FT","FD"):
        draw_survival(results,args.output,region)
        for key in SAMPLES:
            draw_2d(results,args.output,key,region,"angle","#angle(#gamma_{tag},X) (deg)")
            draw_2d(results,args.output,key,region,"dt","#Delta t=t_{p}-t_{#gamma} (GeV^{2})")

    print(f"Done. Output: {args.output}")

if __name__=="__main__":
    main()
