#!/usr/bin/env python3
"""
CLAS12 RGA photon-efficiency study.

Stage 1: characterize and optimize the exclusive-pi0 tag sample without ever
requiring the probe photon to be reconstructed.  The script is intentionally
structured so the same observables, detector-region definitions, probe-momentum
bins, and MC weights can be reused in the later component-normalization,
pi0-purity, probe-matching, efficiency, and closure stages.

Input skim caveat:
PhotonEfficiency already has broad skim requirements
  -1 < Mx2(ep) < 2 GeV^2
  -0.25 < Mx2(ep gamma_tag) < 0.25 GeV^2
and positive missing-photon energy.  Offline plots must be interpreted within
that preselected population.
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
    v = ROOT.std.vector("string")()
    for f in files:
        v.push_back(f)
    return ROOT.RDataFrame(TREE, v)


def define_columns(df, sample_key):
    cols = {str(x) for x in df.GetColumnNames()}
    required = {
        "W", "p_pass_standard", "tag_pass_beta", "tag_pass_fiducial",
        "Mx_ep", "probe_corr_p", "probe_corr_theta", "probe_corr_phi",
        "beam_energy", "e_p", "e_theta", "e_phi",
        "p_corr_p", "p_corr_theta", "p_corr_phi",
        "tag_corr_p", "tag_corr_theta", "tag_corr_phi"
    }
    missing = sorted(required - cols)
    if missing:
        raise RuntimeError("Missing required branches: " + ", ".join(missing))

    mx2ep = "Mx2_ep" if "Mx2_ep" in cols else "(Mx_ep>=0 ? Mx_ep*Mx_ep : -Mx_ep*Mx_ep)"
    mx2epg = "Mx2_epg_corr" if "Mx2_epg_corr" in cols else "0.0"
    weight = "1.0" if sample_key == "data" or "mc_weight" not in cols else "mc_weight"

    d = (df.Define("v_mx2ep", mx2ep)
           .Define("v_mx2epg", mx2epg)
           .Define("v_mx2eg", "pe2_mx2_eg(beam_energy,e_p,e_theta,e_phi,tag_corr_p,tag_corr_theta,tag_corr_phi)")
           .Define("v_copl", "pe2_copl(beam_energy,e_p,e_theta,e_phi,p_corr_p,p_corr_theta,p_corr_phi,tag_corr_p,tag_corr_theta,tag_corr_phi)")
           .Define("v_angle_gX", "pe2_angle(tag_corr_theta,tag_corr_phi,probe_corr_theta,probe_corr_phi)")
           .Define("v_dt", "pe2_delta_t(beam_energy,e_p,e_theta,e_phi,p_corr_p,p_corr_theta,p_corr_phi,tag_corr_p,tag_corr_theta,tag_corr_phi)")
           .Define("probe_region", "probe_corr_theta<=5.5 ? 0 : (probe_corr_theta<=36.0 ? 1 : 2)")
           .Define("ana_weight", weight))

    # Important: the production skim has already required:
    #   -1 < Mx2(ep) < 2 GeV^2
    #   -0.25 < Mx2(ep gamma_tag) < 0.25 GeV^2
    # plus a positive-energy missing-photon hypothesis.
    #
    # This offline baseline adds only reconstructed e/p/tag quality.  It NEVER
    # requires a reconstructed probe photon.
    d = d.Define(
        "baseline",
        "W>2.0 && p_pass_standard==1 && tag_pass_beta==1 && tag_pass_fiducial==1 "
        "&& probe_corr_p>0 && probe_corr_theta>=0"
    )

    # Current working selection.  angle(gamma_tag,X), Delta-t, and any tightening
    # of Mx2(ep gamma_tag) remain diagnostics until their probe-phase-space impact
    # has been examined.
    d = d.Define("cut_mx2ep", "baseline && v_mx2ep>-0.231 && v_mx2ep<0.309")
    d = d.Define("cut_mx2eg", "cut_mx2ep && v_mx2eg>1.4")
    d = d.Define("cut_coplanarity", "cut_mx2eg && fabs(v_copl)<5.7")
    return d


def clone(result):
    h = result.GetValue().Clone()
    h.SetDirectory(0)
    return h


def unit(h):
    q = h.Clone(h.GetName() + "_unit")
    q.SetDirectory(0)
    integ = q.Integral("width")
    if integ > 0:
        q.Scale(1.0 / integ)
    return q


def style(h, key):
    h.SetLineColor(SAMPLES[key][2])
    h.SetMarkerColor(SAMPLES[key][2])
    h.SetLineWidth(2)
    if key == "data":
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.45)


def setup_pad(pad, right=0.04):
    pad.SetLeftMargin(0.14)
    pad.SetRightMargin(right)
    pad.SetBottomMargin(0.14)
    pad.SetTopMargin(0.10)
    pad.SetTicks(1, 1)


# Histogram definitions used both now and in the later normalization stage.
# Keep these centralized so the fit can reuse exactly the same observables.
OBS = {
    "mx2ep":   ("v_mx2ep",       90, -0.50, 1.30, "M_{X}^{2}(ep) (GeV^{2})"),
    "mx2epg":  ("v_mx2epg",     100, -0.25, 0.25, "M_{X}^{2}(ep#gamma_{tag}) (GeV^{2})"),
    "mx2eg":   ("v_mx2eg",      100,  0.00, 6.00, "M_{X}^{2}(e#gamma_{tag}) (GeV^{2})"),
    "copl":    ("v_copl",       120, -30.0, 30.0, "#Delta#phi_{copl} (deg)"),
    "angle":   ("v_angle_gX",    90,   0.0, 30.0, "#theta(#gamma_{tag},X) (deg)"),
    "dt":      ("v_dt",         125,  -1.0, 4.00, "#Delta t=t_{p}-t_{#gamma} (GeV^{2})"),
    "tag_p":   ("tag_corr_p",    90,   0.0, 9.0, "Tag photon momentum p_{#gamma,tag} (GeV)"),
    "tag_th":  ("tag_corr_theta",90,   0.0, 45.0, "Tag photon #theta_{#gamma,tag} (deg)"),
    "tag_phi": ("tag_corr_phi", 120, -180., 180., "Tag photon #phi_{#gamma,tag} (deg)"),
    "probe_p": ("probe_corr_p",  90,   0.0, 9.0, "Predicted probe momentum p_{#gamma,probe} (GeV)"),
    "probe_th":("probe_corr_theta",90,  0.0, 45.0, "Predicted probe #theta_{#gamma,probe} (deg)"),
}

# Predicted-probe momentum slices chosen to expose the FT/high-energy issue and
# retained for the later purity fit.  Last bin includes the beam-energy endpoint.
PROBE_P_BINS = [(0.4,2.0), (2.0,4.0), (4.0,6.0), (6.0,8.0), (8.0,9.0)]


def h1(node, name, obs_key, weight=False):
    col, nb, lo, hi, _ = OBS[obs_key]
    model = (name, "", nb, lo, hi)
    return node.Histo1D(model, col, "ana_weight") if weight else node.Histo1D(model, col)


def h2_probe(node, name, obs_key):
    col, nb, lo, hi, _ = OBS[obs_key]
    return node.Histo2D((name, "", 45, 0.0, 9.0, nb, lo, hi), "probe_corr_p", col)


def book_sample(key, df):
    """Book the cut-optimization study in one RDF graph per sample.

    The scan base deliberately excludes coplanarity and theta(gamma_tag,X), so
    neither variable is preselected before we test it.  The existing Mx2(e gamma)
    >1.4 requirement is retained because the previous study already showed that
    it strongly suppresses the DVCS-like epgamma topology while preserving AAOgen.
    """
    out = {"counts": {}, "shape": {}, "corr": {}, "scan": {}, "survival": {}}

    base = df.Filter("baseline", "baseline")
    # Common optimization base: retain lower Mx2(ep) edge, but leave its upper
    # edge free for the scan.  Do not apply coplanarity or angle cuts here.
    optbase = df.Filter("baseline && v_mx2ep>-0.231 && v_mx2eg>1.4", "optimization base")

    out["counts"]["baseline"] = base.Count()
    out["counts"]["optbase"] = optbase.Count()

    # Shapes at the common optimization base.  These are the distributions that
    # determine the three cuts being optimized, not post-cut versions of them.
    for obs in ("mx2ep", "copl", "angle", "probe_p", "probe_th", "tag_p"):
        out["shape"][("optbase", obs)] = h1(optbase, f"h_{key}_opt_{obs}", obs)

    # Correlations needed to diagnose sculpting of the predicted probe phase space.
    for region, rcut in (("FT", "probe_region==0"), ("FD", "probe_region==1")):
        rn = optbase.Filter(rcut)
        out["corr"][(region,"copl")] = h2_probe(rn, f"h2_{key}_{region}_copl", "copl")
        out["corr"][(region,"angle")] = h2_probe(rn, f"h2_{key}_{region}_angle", "angle")

    # Individual scans.  Values intentionally bracket the visually interesting
    # region rather than assuming the previous values were optimal.
    mx_hi_vals = [0.18,0.20,0.22,0.24,0.26,0.28,0.309]
    anti_vals  = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0]
    ang_vals   = [5.0,6.0,7.0,8.0,9.2,10.0,12.0,15.0,30.0]
    out["scan_values"] = {"mxhi":mx_hi_vals, "anti":anti_vals, "angle":ang_vals}

    for i,x in enumerate(mx_hi_vals):
        out["scan"][("mxhi",i)] = optbase.Filter(f"v_mx2ep<{x}").Count()
    # For the anti-coplanarity scan, 0 means no anti-coplanarity veto.
    for i,x in enumerate(anti_vals):
        expr = "1" if x==0 else f"fabs(v_copl)>{x}"
        out["scan"][("anti",i)] = optbase.Filter(expr).Count()
    # 30 deg is effectively the no-angle-cut reference within the plotted range.
    for i,x in enumerate(ang_vals):
        out["scan"][("angle",i)] = optbase.Filter(f"v_angle_gX<{x}").Count()

    # Combined grid: upper Mx2(ep), anti-coplanarity lower bound, and angle max.
    # Keep it deliberately modest (7*9*9=567 counts/sample) so the optimization
    # is comprehensive without exploding into thousands of histograms.
    for im,mxhi in enumerate(mx_hi_vals):
        for ic,anti in enumerate(anti_vals):
            for ia,amax in enumerate(ang_vals):
                terms=[f"v_mx2ep<{mxhi}", f"v_angle_gX<{amax}"]
                if anti>0: terms.append(f"fabs(v_copl)>{anti}")
                out["scan"][("grid",im,ic,ia)] = optbase.Filter(" && ".join(terms)).Count()

    # Probe-energy survival for representative candidate cuts.  This directly
    # tests the concern that anti-coplanarity may bias pi0 energy sharing.
    candidates = {
        "base": "1",
        "mx024": "v_mx2ep<0.24",
        "anti3": "fabs(v_copl)>3.0",
        "anti5": "fabs(v_copl)>5.0",
        "ang92": "v_angle_gX<9.2",
        "combo_loose": "v_mx2ep<0.24 && fabs(v_copl)>3.0 && v_angle_gX<9.2",
        "combo_tight": "v_mx2ep<0.22 && fabs(v_copl)>5.0 && v_angle_gX<9.2",
    }
    out["candidate_exprs"] = candidates
    for region,rcut in (("FT","probe_region==0"),("FD","probe_region==1")):
        rb=optbase.Filter(rcut)
        for cname,expr in candidates.items():
            node=rb.Filter(expr)
            out["survival"][(region,cname)] = h1(node,f"hsurv_{key}_{region}_{cname}","probe_p")
            out["counts"][(region,cname)] = node.Count()

    return out


def graph_from_scan(results, scan_name, sample_key, denom_key="optbase"):
    vals=results[sample_key]["scan_values"][scan_name]
    den=float(results[sample_key]["counts"][denom_key].GetValue())
    g=ROOT.TGraph(len(vals))
    for i,x in enumerate(vals):
        n=float(results[sample_key]["scan"][(scan_name,i)].GetValue())
        g.SetPoint(i,x,n/den if den>0 else 0.0)
    g.SetLineColor(SAMPLES[sample_key][2]); g.SetMarkerColor(SAMPLES[sample_key][2])
    g.SetLineWidth(2); g.SetMarkerStyle(20); g.SetMarkerSize(0.8)
    return g


def draw_scan_pad(pad, results, scan_name, xtitle, title):
    setup_pad(pad)
    graphs=[]
    frame=ROOT.TH1D(f"frame_{scan_name}_{pad.GetNumber()}","",100,
                    min(results["data"]["scan_values"][scan_name]),
                    max(results["data"]["scan_values"][scan_name]))
    frame.SetDirectory(0); frame.SetMinimum(0); frame.SetMaximum(1.05)
    frame.GetXaxis().SetTitle(xtitle); frame.GetYaxis().SetTitle("Fraction of optimization-base sample retained")
    frame.Draw("AXIS")
    for key in SAMPLES:
        g=graph_from_scan(results,scan_name,key); g.Draw("LP SAME"); graphs.append((key,g))
    leg=ROOT.TLegend(0.57,0.62,0.91,0.87); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.031)
    for key,g in graphs: leg.AddEntry(g,SAMPLES[key][0],"lp")
    leg.Draw()
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043); tx.DrawLatex(0.17,0.92,title)
    return frame,graphs,leg


def draw_shape_overlay_opt(pad, results, obs, title, lines=None):
    setup_pad(pad); hs=[]; ymax=0.0
    for key in SAMPLES:
        h=unit(clone(results[key]["shape"][("optbase",obs)])); style(h,key); hs.append((key,h)); ymax=max(ymax,h.GetMaximum())
    for i,(key,h) in enumerate(hs):
        h.SetTitle(""); h.GetXaxis().SetTitle(OBS[obs][4]); h.GetYaxis().SetTitle("Unit-normalized density")
        h.GetYaxis().SetRangeUser(0,max(1e-12,1.28*ymax)); h.Draw(("E1" if key=="data" else "HIST")+(" SAME" if i else ""))
    keep=[]
    for x in (lines or []):
        l=ROOT.TLine(x,0,x,1.12*ymax); l.SetLineStyle(2); l.SetLineWidth(2); l.Draw(); keep.append(l)
    leg=ROOT.TLegend(0.60,0.64,0.91,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.031)
    for key,h in hs: leg.AddEntry(h,SAMPLES[key][0],"lep" if key=="data" else "l")
    leg.Draw(); tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043); tx.DrawLatex(0.17,0.92,title)
    return hs,keep,leg


def draw_figure1(results,outdir):
    c=ROOT.TCanvas("c01","",1500,1050); c.Divide(2,2); keep=[]
    keep += list(draw_shape_overlay_opt(c.cd(1),results,"mx2ep","Optimization base: M_{X}^{2}(ep)",[0.20,0.24,0.309]))
    keep += list(draw_shape_overlay_opt(c.cd(2),results,"copl","Optimization base: tag-proton coplanarity",[-5,-3,3,5]))
    keep += list(draw_shape_overlay_opt(c.cd(3),results,"angle","Optimization base: #theta(#gamma_{tag},X)",[9.2]))
    keep += list(draw_shape_overlay_opt(c.cd(4),results,"probe_p","Optimization base: predicted probe hypothesis"))
    c.SaveAs(str(Path(outdir)/"01_optimization_base_shapes.png"))


def draw_figure2(results,outdir):
    c=ROOT.TCanvas("c02","",1500,1050); c.Divide(2,2); keep=[]
    keep += list(draw_scan_pad(c.cd(1),results,"mxhi","Upper M_{X}^{2}(ep) cut (GeV^{2})","Scan: reject #eta/other meson region"))
    keep += list(draw_scan_pad(c.cd(2),results,"anti","Minimum |#Delta#phi_{copl}| (deg)","Scan: anti-coplanarity DVCS veto"))
    keep += list(draw_scan_pad(c.cd(3),results,"angle","Maximum #theta(#gamma_{tag},X) (deg)","Scan: missing-photon consistency"))
    p=c.cd(4); setup_pad(p); p.Clear(); tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.046)
    tx.DrawLatex(0.12,0.83,"Optimization base")
    tx.SetTextSize(0.035)
    tx.DrawLatex(0.12,0.70,"M_{X}^{2}(ep)>-0.231 GeV^{2}")
    tx.DrawLatex(0.12,0.62,"M_{X}^{2}(e#gamma_{tag})>1.4 GeV^{2}")
    tx.DrawLatex(0.12,0.54,"No coplanarity or #theta(#gamma_{tag},X) cut")
    tx.DrawLatex(0.12,0.40,"Curves are retention, not purity.")
    tx.DrawLatex(0.12,0.32,"Final choice must also preserve probe phase space.")
    c.SaveAs(str(Path(outdir)/"02_cut_retention_scans.png"))


def ratio_hist(num,den,name):
    r=num.Clone(name); r.SetDirectory(0); r.Divide(den); return r


def draw_survival_pad(pad,results,key,region,title):
    setup_pad(pad); base=clone(results[key]["survival"][(region,"base")]); curves=[]
    names=[("mx024","M_{X}^{2}(ep)<0.24"),("anti3","|#Delta#phi|>3^{#circ}"),("anti5","|#Delta#phi|>5^{#circ}"),("ang92","#theta(#gamma_{tag},X)<9.2^{#circ}"),("combo_loose","combined: 0.24, 3^{#circ}, 9.2^{#circ}"),("combo_tight","combined: 0.22, 5^{#circ}, 9.2^{#circ}")]
    cols=[ROOT.kBlue+1,ROOT.kMagenta+1,ROOT.kRed+1,ROOT.kGreen+2,ROOT.kOrange+7,ROOT.kViolet+1]
    frame=ROOT.TH1D(f"frsurv_{key}_{region}","",90,0,9); frame.SetDirectory(0); frame.SetMinimum(0); frame.SetMaximum(1.08)
    frame.GetXaxis().SetTitle("Predicted probe-hypothesis momentum p_{#gamma,probe} (GeV)"); frame.GetYaxis().SetTitle("Retention relative to optimization base"); frame.Draw("AXIS")
    for (cname,label),col in zip(names,cols):
        h=ratio_hist(clone(results[key]["survival"][(region,cname)]),base,f"r_{key}_{region}_{cname}"); h.SetLineColor(col); h.SetLineWidth(2); h.Draw("HIST SAME"); curves.append((label,h))
    leg=ROOT.TLegend(0.45,0.58,0.91,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.027)
    for label,h in curves: leg.AddEntry(h,label,"l")
    leg.Draw(); tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043); tx.DrawLatex(0.17,0.92,title)
    return frame,curves,leg


def draw_figure3(results,outdir):
    # AAOgen panels are the decisive sculpting test; DATA panels show how much of
    # the observed population each candidate removes in the same phase space.
    c=ROOT.TCanvas("c03","",1500,1050); c.Divide(2,2); keep=[]
    keep += list(draw_survival_pad(c.cd(1),results,"aaogen","FT","AAOgen FT: probe-energy sculpting"))
    keep += list(draw_survival_pad(c.cd(2),results,"data","FT","Data FT: probe-hypothesis retention"))
    keep += list(draw_survival_pad(c.cd(3),results,"aaogen","FD","AAOgen FD: probe-energy sculpting"))
    keep += list(draw_survival_pad(c.cd(4),results,"data","FD","Data FD: probe-hypothesis retention"))
    c.SaveAs(str(Path(outdir)/"03_probe_energy_retention.png"))


def draw_2d_pad_opt(pad,results,key,region,obs,title):
    setup_pad(pad,0.16); h=clone(results[key]["corr"][(region,obs)]); h.SetTitle("")
    h.GetXaxis().SetTitle("Predicted probe-hypothesis momentum p_{#gamma,probe} (GeV)"); h.GetXaxis().SetRangeUser(0,9)
    h.GetYaxis().SetTitle(OBS[obs][4]); h.Draw("COLZ")
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043); tx.DrawLatex(0.17,0.92,title)
    return h


def draw_figure4(results,outdir):
    c=ROOT.TCanvas("c04","",1500,1050); c.Divide(2,2); keep=[]
    keep.append(draw_2d_pad_opt(c.cd(1),results,"data","FT","copl","Data FT: coplanarity vs probe hypothesis"))
    keep.append(draw_2d_pad_opt(c.cd(2),results,"aaogen","FT","copl","AAOgen FT: coplanarity vs probe hypothesis"))
    keep.append(draw_2d_pad_opt(c.cd(3),results,"data","FT","angle","Data FT: #theta(#gamma_{tag},X) vs probe hypothesis"))
    keep.append(draw_2d_pad_opt(c.cd(4),results,"aaogen","FT","angle","AAOgen FT: #theta(#gamma_{tag},X) vs probe hypothesis"))
    c.SaveAs(str(Path(outdir)/"04_FT_cut_correlations.png"))


def best_grid_rows(results):
    """Rank cut combinations by MC separation only, without pretending MC counts
    are physically normalized.  Score = geometric mean of AAO retention and the
    two rejection fractions.  It is a diagnostic ordering, not a purity estimate.
    """
    vals=results["aaogen"]["scan_values"]
    den={k:float(results[k]["counts"]["optbase"].GetValue()) for k in ("aaogen","clasdis","dvcsgen")}
    rows=[]
    for im,mxhi in enumerate(vals["mxhi"]):
      for ic,anti in enumerate(vals["anti"]):
       for ia,amax in enumerate(vals["angle"]):
        fr={}
        for k in den:
            n=float(results[k]["scan"][("grid",im,ic,ia)].GetValue()); fr[k]=n/den[k] if den[k]>0 else 0
        # Favor AAO retention and rejection of both backgrounds, but do not call
        # this purity because absolute MC normalizations are not yet established.
        score=(max(fr["aaogen"],0)*max(1-fr["clasdis"],0)*max(1-fr["dvcsgen"],0))**(1/3)
        rows.append((score,mxhi,anti,amax,fr["aaogen"],fr["clasdis"],fr["dvcsgen"]))
    return sorted(rows,reverse=True)


def write_tables(results,outdir):
    p=Path(outdir)
    with open(p/"cut_optimization_scan.txt","w") as f:
        f.write("# Fractions are relative to the common optimization base.\n")
        f.write("# They are retention/rejection diagnostics, NOT purity estimates; MC components are not normalized yet.\n")
        f.write("# score mx2ep_upper anti_copl_min angle_max aaogen_ret clasdis_ret dvcsgen_ret\n")
        for row in best_grid_rows(results)[:80]:
            f.write("%.6f %.3f %.1f %.1f %.6f %.6f %.6f\n"%row)
    with open(p/"candidate_cut_counts.txt","w") as f:
        f.write("# Unweighted counts after the common optimization base and representative candidate cuts.\n")
        f.write("sample region candidate count\n")
        for key in SAMPLES:
            for region in ("FT","FD"):
                for cname in results[key]["candidate_exprs"]:
                    f.write(f"{key} {region} {cname} {int(results[key]['counts'][(region,cname)].GetValue())}\n")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=8,help="ROOT worker threads (1-8; default 8)")
    ap.add_argument("--output",default="output_RGA_study")
    args=ap.parse_args(); workers=max(1,min(8,args.workers)); ROOT.EnableImplicitMT(workers)
    Path(args.output).mkdir(parents=True,exist_ok=True)
    print("="*78)
    print("Photon-efficiency RGA study - exclusive-pi0 cut optimization")
    print(f"ROOT worker threads: {workers}")
    print("No reconstructed-probe requirement.")
    print("Optimization base: Mx2(ep)>-0.231 and Mx2(e gamma_tag)>1.4; no coplanarity/angle cut.")
    print("Scanning Mx2(ep) upper edge, anti-coplanarity veto, and theta(gamma_tag,X) maximum.")
    print("Photon-energy plotting range: 0-9 GeV.")
    print("="*78)
    results={}; actions=[]
    for key,(label,pat,color) in SAMPLES.items():
        fs=files_for(pat)
        if not fs: raise RuntimeError(f"No ROOT files for {key}: {pat}")
        print(f"{label:10s}: {len(fs):3d} ROOT files")
        results[key]=book_sample(key,define_columns(make_rdf(fs),key))
        for group in ("counts","shape","corr","scan","survival"):
            actions += list(results[key][group].values())
    print(f"Executing {len(actions)} booked actions ...")
    ROOT.RDF.RunGraphs(actions)
    write_tables(results,args.output)
    draw_figure1(results,args.output); draw_figure2(results,args.output)
    draw_figure3(results,args.output); draw_figure4(results,args.output)
    print("Top 10 MC-separation scan points (diagnostic only; not purity):")
    print(" score  Mx2hi anti angle  AAOret CLASDISret DVCSret")
    for r in best_grid_rows(results)[:10]:
        print(f" {r[0]:.3f}  {r[1]:.3f}  {r[2]:.1f}  {r[3]:.1f}   {r[4]:.3f}   {r[5]:.3f}      {r[6]:.3f}")
    print(f"Done. Four focused figures + two text tables written to: {args.output}")


if __name__=="__main__":
    main()
