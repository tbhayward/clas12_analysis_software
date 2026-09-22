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
import math
import itertools
import numpy as np

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


# Final-validation candidates.  Keep the names stable: these definitions are
# intended to flow directly into the subsequent component-normalization and
# efficiency stages.
CANDIDATES = {
    "A": "v_mx2ep<0.24 && fabs(v_copl)>2.0 && v_angle_gX<9.2",
    "B": "v_mx2ep<0.24 && fabs(v_copl)>2.0 && v_angle_gX<8.0",
    "C": "v_mx2ep<0.22 && fabs(v_copl)>2.0 && v_angle_gX<8.0",
    "D": "v_mx2ep<0.24 && fabs(v_copl)>3.0 && v_angle_gX<8.0",
}
CANDIDATE_LABELS = {
    "A": "A: 0.24, 2^{#circ}, 9.2^{#circ}",
    "B": "B: 0.24, 2^{#circ}, 8^{#circ}",
    "C": "C: 0.22, 2^{#circ}, 8^{#circ}",
    "D": "D: 0.24, 3^{#circ}, 8^{#circ}",
}
CANDIDATE_COLORS = {"A": ROOT.kBlue+1, "B": ROOT.kGreen+2,
                    "C": ROOT.kMagenta+1, "D": ROOT.kOrange+7}

# Candidate B is the reference only for *validation plots*.  No final cut is
# hard-coded by this choice; all four candidates are counted and compared.
REFERENCE_CANDIDATE = "B"
NOMINAL_CANDIDATE = "B"
FIT_P_RANGE = (0.4, 8.0)

# Control regions are deliberately disjoint in purpose, not used in the
# efficiency denominator.  They are retained now so the later normalization
# fit can be validated against background-enriched DATA regions.
CONTROL_REGIONS = {
    "eta": (
        "v_mx2ep>0.26 && v_mx2ep<0.38 && fabs(v_copl)>2.0 && v_angle_gX<8.0",
        "Competing-meson control: 0.26<M_{X}^{2}(ep)<0.38"
    ),
    "central": (
        "v_mx2ep<0.24 && fabs(v_copl)<2.0 && v_angle_gX<8.0",
        "Central-coplanarity control: |#Delta#phi_{copl}|<2^{#circ}"
    ),
    "highangle": (
        "v_mx2ep<0.24 && fabs(v_copl)>2.0 && v_angle_gX>12.0 && v_angle_gX<25.0",
        "High-angle control: 12^{#circ}<#theta(#gamma_{tag},X)<25^{#circ}"
    ),
}

# Four mutually exclusive regions used by the Stage-2 component-normalization fit.
# Candidate B is the pi0-enriched signal region.  The other three are controls.
FIT_REGIONS = {
    "signal": CANDIDATES[NOMINAL_CANDIDATE],
    "eta": CONTROL_REGIONS["eta"][0],
    "central": CONTROL_REGIONS["central"][0],
    "highangle": CONTROL_REGIONS["highangle"][0],
}
FIT_COMPONENTS = ("aaogen", "clasdis", "dvcsgen")


def book_sample(key, df):
    """Book overview, denominator validation, and normalization inputs in one RDF graph."""
    out = {"counts": {}, "survival": {}, "selected": {}, "slices": {},
           "control": {}, "overview": {}, "fit": {}}

    # Restore the useful early distributions.  These are intentionally booked
    # before the final Candidate-B cuts so a reader can see why each cut exists.
    base = df.Filter("baseline", "baseline")
    out["counts"]["baseline"] = base.Count()
    for obs in ("mx2ep", "mx2eg", "copl", "angle", "mx2epg", "tag_p", "probe_p"):
        out["overview"][("baseline", obs)] = h1(base, f"hov_{key}_base_{obs}", obs)
    n1_mx2ep = base.Filter("v_mx2eg>1.4")
    n1_mx2eg = base.Filter("v_mx2ep>-0.231 && v_mx2ep<0.24")
    n1_copl = base.Filter("v_mx2ep>-0.231 && v_mx2ep<0.24 && v_mx2eg>1.4 && v_angle_gX<8.0")
    n1_angle = base.Filter("v_mx2ep>-0.231 && v_mx2ep<0.24 && v_mx2eg>1.4 && fabs(v_copl)>2.0")
    for tag,node,obs in (("n1_mx2ep",n1_mx2ep,"mx2ep"),("n1_mx2eg",n1_mx2eg,"mx2eg"),
                         ("n1_copl",n1_copl,"copl"),("n1_angle",n1_angle,"angle")):
        out["overview"][(tag,obs)] = h1(node, f"hov_{key}_{tag}_{obs}", obs)

    optbase = df.Filter("baseline && v_mx2ep>-0.231 && v_mx2eg>1.4", "optimization base")
    out["counts"]["optbase"] = optbase.Count()

    for region, rcut in (("FT", "probe_region==0"), ("FD", "probe_region==1")):
        rb = optbase.Filter(rcut)
        out["counts"][(region, "base")] = rb.Count()
        out["survival"][(region, "base")] = h1(rb, f"hs_{key}_{region}_base", "probe_p")

        for cname, expr in CANDIDATES.items():
            node = rb.Filter(expr, f"candidate {cname} {region}")
            out["counts"][(region, cname)] = node.Count()
            out["survival"][(region, cname)] = h1(node, f"hs_{key}_{region}_{cname}", "probe_p")
            for obs in ("mx2epg", "tag_p", "tag_th", "tag_phi", "dt", "probe_p", "probe_th"):
                out["selected"][(region, cname, obs)] = h1(node, f"hsel_{key}_{region}_{cname}_{obs}", obs)
            for ib, (plo, phi) in enumerate(PROBE_P_BINS):
                ps = node.Filter(f"probe_corr_p>={plo} && probe_corr_p<{phi}")
                out["counts"][(region, cname, "pbin", ib)] = ps.Count()
                out["slices"][(region, cname, ib, "mx2epg")] = h1(ps, f"hsl_{key}_{region}_{cname}_{ib}_mx2epg", "mx2epg")
                out["slices"][(region, cname, ib, "dt")] = h1(ps, f"hsl_{key}_{region}_{cname}_{ib}_dt", "dt")
                # Weighted local-purity template.  DATA has ana_weight=1, while
                # MC keeps its native event weight.  Candidate A/C/D are booked
                # now so their spread can become a denominator-selection
                # systematic without another pass over the ROOT trees.
                out["fit"][(region, "local", cname, ib, "mx2epg")] = h1(
                    ps, f"hlocal_{key}_{region}_{cname}_{ib}_mx2epg",
                    "mx2epg", weight=(key!="data"))

        # Stage-2 normalization inputs.  The primary fit uses one shared set of
        # component scale factors across FT+FD.  Separate FT/FD fits are retained
        # as a closure diagnostic, not as the nominal normalization.
        for fr, expr in FIT_REGIONS.items():
            fn = rb.Filter(expr, f"fit region {fr} {region}")
            out["fit"][(region, fr, "count")] = fn.Count()
            out["fit"][(region, fr, "sumw")] = fn.Sum("ana_weight")
            out["fit"][(region, fr, "sumw2")] = fn.Define(
                f"w2_{region}_{fr}_{key}", "ana_weight*ana_weight").Sum(f"w2_{region}_{fr}_{key}")
            out["fit"][(region, fr, "probe_p")] = h1(fn, f"hfit_{key}_{region}_{fr}_p", "probe_p", weight=(key!="data"))
            # Shape reserved for independent post-fit closure; it is NOT used to determine coefficients.
            closure_obs = "mx2epg" if fr=="signal" else ("mx2ep" if fr=="eta" else ("copl" if fr=="central" else "angle"))
            out["fit"][(region, fr, "shape")] = h1(fn, f"hfit_{key}_{region}_{fr}_shape", closure_obs, weight=(key!="data"))
            for ib,(plo,phi) in enumerate(PROBE_P_BINS):
                pn=fn.Filter(f"probe_corr_p>={plo} && probe_corr_p<{phi}")
                out["fit"][(region, fr, ib, "count")] = pn.Count()
                out["fit"][(region, fr, ib, "sumw")] = pn.Sum("ana_weight")

    # Legacy compact control plots, integrated over predicted detector region.
    for cr, (expr, _) in CONTROL_REGIONS.items():
        node = optbase.Filter(expr, f"control {cr}")
        out["counts"][("control", cr)] = node.Count()
        obs = "mx2ep" if cr == "eta" else ("copl" if cr == "central" else "angle")
        out["control"][(cr, obs)] = h1(node, f"hcr_{key}_{cr}_{obs}", obs)
        out["control"][(cr, "probe_p")] = h1(node, f"hcr_{key}_{cr}_probe_p", "probe_p")
    return out


def ratio_hist(num, den, name):
    r = num.Clone(name); r.SetDirectory(0); r.Divide(den); return r


def draw_candidate_retention_pad(pad, results, key, region, title):
    setup_pad(pad)
    base = clone(results[key]["survival"][(region, "base")])
    frame = ROOT.TH1D(f"fr_{key}_{region}", "", 90, 0, 9); frame.SetDirectory(0)
    frame.SetMinimum(0); frame.SetMaximum(1.08)
    frame.GetXaxis().SetTitle("Predicted probe-hypothesis momentum p_{#gamma,probe} (GeV)")
    frame.GetYaxis().SetTitle("Candidate retention / optimization base")
    frame.Draw("AXIS")
    curves = []
    for cname in CANDIDATES:
        h = ratio_hist(clone(results[key]["survival"][(region, cname)]), base,
                       f"rr_{key}_{region}_{cname}")
        h.SetLineColor(CANDIDATE_COLORS[cname]); h.SetLineWidth(2); h.Draw("HIST SAME")
        curves.append((cname, h))
    leg = ROOT.TLegend(0.48, 0.64, 0.91, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.029)
    for cname, h in curves: leg.AddEntry(h, CANDIDATE_LABELS[cname], "l")
    leg.Draw()
    tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043); tx.DrawLatex(0.17, 0.92, title)
    return frame, curves, leg


def draw_overlay(pad, results, getter, obs, title, lines=None):
    setup_pad(pad); hs=[]; ymax=0.0
    for key in SAMPLES:
        h = unit(clone(getter(results[key], key))); style(h, key)
        hs.append((key, h)); ymax=max(ymax, h.GetMaximum())
    for i, (key, h) in enumerate(hs):
        h.SetTitle(""); h.GetXaxis().SetTitle(OBS[obs][4]); h.GetYaxis().SetTitle("Unit-normalized density")
        h.GetYaxis().SetRangeUser(0, max(1e-12, 1.28*ymax))
        h.Draw(("E1" if key=="data" else "HIST") + (" SAME" if i else ""))
    ls=[]
    for x in (lines or []):
        l=ROOT.TLine(x,0,x,1.12*ymax); l.SetLineStyle(2); l.SetLineWidth(2); l.Draw(); ls.append(l)
    leg=ROOT.TLegend(0.60,0.64,0.91,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.031)
    for key,h in hs: leg.AddEntry(h,SAMPLES[key][0],"lep" if key=="data" else "l")
    leg.Draw(); tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.041); tx.DrawLatex(0.17,0.92,title)
    return hs, ls, leg


def draw_figure1(results, outdir):
    c=ROOT.TCanvas("c01","",1500,1050); c.Divide(2,2); keep=[]
    keep += list(draw_candidate_retention_pad(c.cd(1),results,"aaogen","FT","AAOgen FT: candidate stability"))
    keep += list(draw_candidate_retention_pad(c.cd(2),results,"data","FT","Data FT: candidate stability"))
    keep += list(draw_candidate_retention_pad(c.cd(3),results,"aaogen","FD","AAOgen FD: candidate stability"))
    keep += list(draw_candidate_retention_pad(c.cd(4),results,"data","FD","Data FD: candidate stability"))
    c.SaveAs(str(Path(outdir)/"01_candidate_probe_phase_space.png"))


def draw_figure2(results, outdir):
    c=ROOT.TCanvas("c02","",1500,1050); c.Divide(2,2); keep=[]; cand=REFERENCE_CANDIDATE
    keep += list(draw_overlay(c.cd(1),results,lambda r,k:r["selected"][("FT",cand,"tag_p")],"tag_p",f"Candidate {cand}, FT: tag momentum"))
    keep += list(draw_overlay(c.cd(2),results,lambda r,k:r["selected"][("FT",cand,"tag_th")],"tag_th",f"Candidate {cand}, FT: tag polar angle"))
    keep += list(draw_overlay(c.cd(3),results,lambda r,k:r["selected"][("FT",cand,"mx2epg")],"mx2epg",f"Candidate {cand}, FT: M_{{X}}^{{2}}(ep#gamma_{{tag}})"))
    keep += list(draw_overlay(c.cd(4),results,lambda r,k:r["selected"][("FD",cand,"mx2epg")],"mx2epg",f"Candidate {cand}, FD: M_{{X}}^{{2}}(ep#gamma_{{tag}})"))
    c.SaveAs(str(Path(outdir)/"02_reference_candidate_shapes.png"))


def draw_figure3(results, outdir):
    # FT is the difficult detector region and the immediate motivation for this
    # study.  Show the normalization-sensitive loose-skim missing-mass shape in
    # four broad probe-energy intervals.  The 8-9 GeV bin is retained in tables
    # but omitted here because it is expected to be sparse.
    c=ROOT.TCanvas("c03","",1500,1050); c.Divide(2,2); keep=[]; cand=REFERENCE_CANDIDATE
    for ipad, ib in enumerate(range(4), start=1):
        lo,hi=PROBE_P_BINS[ib]
        keep += list(draw_overlay(c.cd(ipad),results,
            lambda r,k,ib=ib:r["slices"][("FT",cand,ib,"mx2epg")],
            "mx2epg",f"Candidate {cand}, FT: {lo:g}<p_{{probe}}<{hi:g} GeV"))
    c.SaveAs(str(Path(outdir)/"03_FT_normalization_shape_vs_probe_energy.png"))


def draw_figure4(results, outdir):
    c=ROOT.TCanvas("c04","",1500,1050); c.Divide(2,2); keep=[]
    keep += list(draw_overlay(c.cd(1),results,lambda r,k:r["control"][("eta","mx2ep")],"mx2ep",CONTROL_REGIONS["eta"][1],[0.26,0.38]))
    central = draw_overlay(c.cd(2),results,lambda r,k:r["control"][("central","copl")],"copl",CONTROL_REGIONS["central"][1],[-2,2])
    for _key,_h in central[0]:
        _h.GetXaxis().SetRangeUser(-3.0,3.0)
    keep += list(central)
    keep += list(draw_overlay(c.cd(3),results,lambda r,k:r["control"][("highangle","angle")],"angle",CONTROL_REGIONS["highangle"][1],[12,25]))
    cand=REFERENCE_CANDIDATE
    keep += list(draw_overlay(c.cd(4),results,lambda r,k:r["selected"][("FT",cand,"dt")],"dt",f"Candidate {cand}, FT: #Delta t validation"))
    c.SaveAs(str(Path(outdir)/"04_background_control_regions.png"))



def draw_selection_overview(results, outdir):
    c=ROOT.TCanvas("coverview","",1500,1050); c.Divide(2,2); keep=[]
    keep += list(draw_overlay(c.cd(1),results,lambda r,k:r["overview"][("n1_mx2ep","mx2ep")],"mx2ep","N-1: M_{X}^{2}(ep)",[-0.231,0.24]))
    keep += list(draw_overlay(c.cd(2),results,lambda r,k:r["overview"][("n1_mx2eg","mx2eg")],"mx2eg","N-1: M_{X}^{2}(e#gamma_{tag})",[1.4]))
    keep += list(draw_overlay(c.cd(3),results,lambda r,k:r["overview"][("n1_copl","copl")],"copl","N-1: tag-proton coplanarity",[-2,2]))
    keep += list(draw_overlay(c.cd(4),results,lambda r,k:r["overview"][("n1_angle","angle")],"angle","N-1: #theta(#gamma_{tag},X)",[8]))
    c.SaveAs(str(Path(outdir)/"01_exclusivity_Nminus1.png"))

    c2=ROOT.TCanvas("coverview2","",1500,1050); c2.Divide(2,2); keep2=[]
    keep2 += list(draw_overlay(c2.cd(1),results,lambda r,k:r["overview"][("baseline","mx2epg")],"mx2epg","Loose-skim M_{X}^{2}(ep#gamma_{tag})"))
    keep2 += list(draw_overlay(c2.cd(2),results,lambda r,k:r["overview"][("baseline","tag_p")],"tag_p","Baseline tag-photon momentum"))
    keep2 += list(draw_overlay(c2.cd(3),results,lambda r,k:r["overview"][("baseline","probe_p")],"probe_p","Baseline predicted-probe momentum"))
    keep2 += list(draw_overlay(c2.cd(4),results,lambda r,k:r["overview"][("baseline","angle")],"angle","Baseline #theta(#gamma_{tag},X)"))
    c2.SaveAs(str(Path(outdir)/"02_baseline_kinematics.png"))


def _fit_nonnegative(A, y, sigma):
    """Tiny non-negative weighted least-squares solver for three components."""
    A=np.asarray(A,float); y=np.asarray(y,float); sigma=np.asarray(sigma,float)
    good=np.isfinite(y)&np.isfinite(sigma)&(sigma>0)&np.all(np.isfinite(A),axis=1)
    A=A[good]; y=y[good]; sigma=sigma[good]
    Aw=A/sigma[:,None]; yw=y/sigma
    best=None
    n=A.shape[1]
    for mask in range(1,1<<n):
        active=[j for j in range(n) if mask&(1<<j)]
        x=np.zeros(n)
        sol,*_=np.linalg.lstsq(Aw[:,active],yw,rcond=None)
        if np.any(sol<0): continue
        x[active]=sol
        resid=(A@x-y)/sigma
        chi2=float(resid@resid)
        if best is None or chi2<best[0]: best=(chi2,x)
    if best is None: return np.zeros(n), float("inf"), 0
    chi2,x=best
    ndf=max(0,len(y)-np.count_nonzero(x>0))
    return x,chi2,ndf


def perform_normalization_fit(results):
    """Fit AAOgen/CLASDIS/DVCSgen scale factors from region yields only.

    Nominal coefficients are shared across FT and FD.  Detector-specific fits are
    diagnostic.  Shape histograms are deliberately excluded from the fit so they
    remain genuine closure tests.
    """
    fits={}
    for fitname,regions in (("shared",("FT","FD")),("FT",("FT",)),("FD",("FD",))):
        rows=[]; y=[]; sig=[]; labels=[]
        for det in regions:
            for fr in FIT_REGIONS:
                yd=float(results["data"]["fit"][(det,fr,"count")].GetValue())
                mc=[float(results[k]["fit"][(det,fr,"sumw")].GetValue()) for k in FIT_COMPONENTS]
                if yd<=0 or sum(mc)<=0: continue
                rows.append(mc); y.append(yd); sig.append(math.sqrt(max(yd,1.0))); labels.append((det,fr))
        x,chi2,ndf=_fit_nonnegative(rows,y,sig)
        fits[fitname]={"scale":dict(zip(FIT_COMPONENTS,x)),"chi2":chi2,"ndf":ndf,"rows":labels,
                       "A":np.asarray(rows,float),"y":np.asarray(y,float),"sigma":np.asarray(sig,float)}
    return fits


def normalization_products(results, fits):
    """Compute post-fit region closure and differential pi0 purity from nominal shared fit."""
    scales=fits["shared"]["scale"]
    rows=[]
    for det in ("FT","FD"):
        for fr in FIT_REGIONS:
            data=float(results["data"]["fit"][(det,fr,"count")].GetValue())
            comp={k:scales[k]*float(results[k]["fit"][(det,fr,"sumw")].GetValue()) for k in FIT_COMPONENTS}
            pred=sum(comp.values())
            rows.append((det,fr,data,comp,pred))
    purity=[]
    for det in ("FT","FD"):
        for ib,(lo,hi) in enumerate(PROBE_P_BINS):
            data=float(results["data"]["fit"][(det,"signal",ib,"count")].GetValue())
            comp={k:scales[k]*float(results[k]["fit"][(det,"signal",ib,"sumw")].GetValue()) for k in FIT_COMPONENTS}
            pred=sum(comp.values()); p=comp["aaogen"]/pred if pred>0 else float("nan")
            purity.append((det,ib,lo,hi,data,comp,pred,p))
    return rows,purity


def draw_normalization_summary(results, fits, outdir):
    scales=fits["shared"]["scale"]
    # Region-yield closure: absolute yields, not unit-normalized shapes.
    c=ROOT.TCanvas("cnorm","",1500,850); c.Divide(2,1); keep=[]
    for ipad,det in enumerate(("FT","FD"),1):
        pad=c.cd(ipad); setup_pad(pad)
        frame=ROOT.TH1D(f"frnorm_{det}","",4,0,4); frame.SetDirectory(0)
        for i,fr in enumerate(FIT_REGIONS,1): frame.GetXaxis().SetBinLabel(i,fr)
        datah=frame.Clone(f"datanorm_{det}"); datah.SetDirectory(0); datah.SetMarkerStyle(20); datah.SetLineColor(ROOT.kBlack)
        predh=frame.Clone(f"prednorm_{det}"); predh.SetDirectory(0); predh.SetLineWidth(3); predh.SetLineColor(ROOT.kRed+1)
        ymax=0
        for i,fr in enumerate(FIT_REGIONS,1):
            d=float(results["data"]["fit"][(det,fr,"count")].GetValue())
            pred=sum(scales[k]*float(results[k]["fit"][(det,fr,"sumw")].GetValue()) for k in FIT_COMPONENTS)
            datah.SetBinContent(i,d); datah.SetBinError(i,math.sqrt(max(d,1)))
            predh.SetBinContent(i,pred); ymax=max(ymax,d,pred)
        frame.SetMinimum(0); frame.SetMaximum(1.25*ymax)
        frame.GetYaxis().SetTitle("Events")
        frame.GetXaxis().SetTitle("Fit/control region")
        frame.Draw("AXIS")
        predh.Draw("HIST SAME"); datah.Draw("E1 SAME")
        leg=ROOT.TLegend(.58,.72,.90,.87); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.AddEntry(datah,"Data","lep"); leg.AddEntry(predh,"Post-fit MC sum","l"); leg.Draw()
        tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(.045); tx.DrawLatex(.17,.92,f"{det}: four-region yield closure")
        q=fits[det]
        tx2=ROOT.TLatex(); tx2.SetNDC(); tx2.SetTextSize(.032)
        tx2.DrawLatex(.17,.855,f"diagnostic fit: #chi^{{2}}/ndf = {q['chi2']:.1f}/{q['ndf']}")
        keep += [frame,datah,predh,leg,tx,tx2]
    c.SaveAs(str(Path(outdir)/"01_four_region_yield_closure.png"))

    # Differential signal-region composition: this is the quantity that will feed the denominator correction.
    c2=ROOT.TCanvas("cpurity","",1500,850); c2.Divide(2,1); keep2=[]
    for ipad,det in enumerate(("FT","FD"),1):
        pad=c2.cd(ipad); setup_pad(pad)
        h=ROOT.TH1D(f"hpur_{det}","",len(PROBE_P_BINS),0,len(PROBE_P_BINS)); h.SetDirectory(0); h.SetMarkerStyle(20); h.SetLineWidth(2)
        for i,(lo,hi) in enumerate(PROBE_P_BINS,1):
            comp={k:scales[k]*float(results[k]["fit"][(det,"signal",i-1,"sumw")].GetValue()) for k in FIT_COMPONENTS}
            tot=sum(comp.values()); h.SetBinContent(i,comp["aaogen"]/tot if tot>0 else 0); h.GetXaxis().SetBinLabel(i,f"{lo:g}-{hi:g}")
        h.SetMinimum(0); h.SetMaximum(1.05); h.GetYaxis().SetTitle("Post-fit #pi^{0} fraction"); h.GetXaxis().SetTitle("Predicted probe momentum (GeV)"); h.Draw("E1")
        tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(.045); tx.DrawLatex(.17,.92,f"{det}: Candidate B composition")
        keep2 += [h,tx]
    c2.SaveAs(str(Path(outdir)/"02_pi0_fraction_vs_probe_momentum.png"))

    # Independent shape closure.  Coefficients came only from four integrated yields.
    for det in ("FT","FD"):
        c3=ROOT.TCanvas(f"cshape_{det}","",1500,1050); c3.Divide(2,2); keep3=[]
        for ipad,fr in enumerate(FIT_REGIONS,1):
            pad=c3.cd(ipad); setup_pad(pad)
            hd=clone(results["data"]["fit"][(det,fr,"shape")]); hd.SetMarkerStyle(20); hd.SetLineColor(ROOT.kBlack)
            hs=[]; total=None
            for k in FIT_COMPONENTS:
                hm=clone(results[k]["fit"][(det,fr,"shape")]); hm.Scale(scales[k]); style(hm,k); hs.append((k,hm))
                if total is None: total=hm.Clone(f"htot_{det}_{fr}"); total.SetDirectory(0)
                else: total.Add(hm)
            ymax=max(hd.GetMaximum(),total.GetMaximum())*1.28 if total else hd.GetMaximum()*1.28
            hd.SetTitle(""); hd.GetYaxis().SetTitle("Events"); hd.GetYaxis().SetRangeUser(0,max(ymax,1)); hd.Draw("E1")
            for k,hm in hs: hm.Draw("HIST SAME")
            total.SetLineColor(ROOT.kMagenta+2); total.SetLineWidth(3); total.Draw("HIST SAME"); hd.Draw("E1 SAME")
            leg=ROOT.TLegend(.57,.61,.91,.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(.028); leg.AddEntry(hd,"Data","lep")
            for k,hm in hs: leg.AddEntry(hm,SAMPLES[k][0],"l")
            leg.AddEntry(total,"Post-fit sum","l"); leg.Draw()
            tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(.041); tx.DrawLatex(.17,.92,f"{det}: {fr} shape closure")
            keep3 += [hd,total,leg,tx]+[x[1] for x in hs]
        c3.SaveAs(str(Path(outdir)/f"03_{det}_independent_shape_closure.png"))



def _hist_arrays(h):
    """Return bin contents and errors without under/overflow."""
    n=h.GetNbinsX()
    y=np.array([h.GetBinContent(i) for i in range(1,n+1)],dtype=float)
    e=np.array([h.GetBinError(i) for i in range(1,n+1)],dtype=float)
    return y,e


def _fit_hist_templates(data_h, template_hists):
    """Non-negative binned template fit.

    The fit is deliberately local: one detector region, one predicted-probe
    momentum bin, and one candidate definition.  DATA statistical uncertainty
    defines the nominal weighting.  MC statistical precision is reported
    separately and can be promoted into the fit uncertainty once closure is
    accepted.
    """
    y,_=_hist_arrays(data_h)
    cols=[]
    for h in template_hists:
        a,_=_hist_arrays(h); cols.append(a)
    A=np.column_stack(cols)
    # Empty DATA bins still carry finite Poisson information; use variance >=1.
    sigma=np.sqrt(np.maximum(y,1.0))
    x,chi2,ndf=_fit_nonnegative(A,y,sigma)

    pred=A@x
    active=np.where(x>0)[0]
    xerr=np.full(len(x),np.nan)
    if len(active):
        Aw=A[:,active]/sigma[:,None]
        try:
            cov=np.linalg.inv(Aw.T@Aw)
        except np.linalg.LinAlgError:
            cov=np.linalg.pinv(Aw.T@Aw)
        errs=np.sqrt(np.maximum(np.diag(cov),0.0))
        for j,e in zip(active,errs): xerr[j]=e

    integrals=np.array([h.Integral() for h in template_hists],dtype=float)
    yields=x*integrals
    yerrs=xerr*integrals
    total=float(np.sum(yields))
    frac=float(yields[0]/total) if total>0 else float("nan")

    # A simple shape-conditioning diagnostic.  Large values warn that two
    # templates are too similar for their individual normalizations to be
    # determined robustly from this distribution alone.
    good=np.sum(A,axis=1)>0
    cond=float("inf")
    if np.count_nonzero(good)>len(template_hists):
        M=A[good]
        norms=np.sqrt(np.sum(M*M,axis=0))
        ok=norms>0
        if np.count_nonzero(ok)>=2:
            cond=float(np.linalg.cond(M[:,ok]/norms[ok]))

    return {
        "scale":dict(zip(FIT_COMPONENTS,x)),
        "scale_err":dict(zip(FIT_COMPONENTS,xerr)),
        "yield":dict(zip(FIT_COMPONENTS,yields)),
        "yield_err":dict(zip(FIT_COMPONENTS,yerrs)),
        "total":total, "pi0_fraction":frac,
        "chi2":chi2, "ndf":ndf, "condition":cond,
        "pred_bins":pred,
    }


def perform_local_purity_fits(results):
    """Fit Candidate-B Mx2(ep gamma_tag) directly in detector/p_probe bins.

    Candidate B is nominal.  A/C/D are fitted in parallel so the same machinery
    is ready for the cut-variation systematic.  The 8-9 GeV bin is retained as
    a diagnostic but is not intended to define an independent final efficiency
    point unless its statistics prove adequate.
    """
    out={}
    for det in ("FT","FD"):
        for cname in CANDIDATES:
            for ib,(lo,hi) in enumerate(PROBE_P_BINS):
                hd=clone(results["data"]["fit"][(det,"local",cname,ib,"mx2epg")])
                hm=[clone(results[k]["fit"][(det,"local",cname,ib,"mx2epg")]) for k in FIT_COMPONENTS]
                out[(det,cname,ib)] = _fit_hist_templates(hd,hm)
    return out


def draw_local_purity_fits(results, localfits, outdir):
    """Presentation-oriented local fits: four useful p bins per detector."""
    p=Path(outdir); p.mkdir(parents=True,exist_ok=True)
    cname=NOMINAL_CANDIDATE

    for det in ("FT","FD"):
        c=ROOT.TCanvas(f"clocal_{det}","",1600,1100); c.Divide(2,2); keep=[]
        for ipad,ib in enumerate(range(4),1):
            lo,hi=PROBE_P_BINS[ib]
            pad=c.cd(ipad); setup_pad(pad)
            hd=clone(results["data"]["fit"][(det,"local",cname,ib,"mx2epg")])
            hd.SetMarkerStyle(20); hd.SetMarkerSize(.45); hd.SetLineColor(ROOT.kBlack)
            fit=localfits[(det,cname,ib)]
            comps=[]; total=None
            for k in FIT_COMPONENTS:
                h=clone(results[k]["fit"][(det,"local",cname,ib,"mx2epg")])
                h.Scale(fit["scale"][k]); style(h,k); comps.append((k,h))
                if total is None:
                    total=h.Clone(f"hloc_tot_{det}_{ib}"); total.SetDirectory(0)
                else:
                    total.Add(h)
            ymax=max(hd.GetMaximum(),total.GetMaximum())*1.30 if total else hd.GetMaximum()*1.30
            hd.SetTitle("")
            hd.GetXaxis().SetTitle(OBS["mx2epg"][4])
            hd.GetYaxis().SetTitle("Events")
            hd.GetYaxis().SetRangeUser(0,max(1.0,ymax))
            hd.Draw("E1")
            for k,h in comps: h.Draw("HIST SAME")
            total.SetLineColor(ROOT.kMagenta+2); total.SetLineWidth(3); total.Draw("HIST SAME")
            hd.Draw("E1 SAME")
            leg=ROOT.TLegend(.55,.59,.91,.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(.027)
            leg.AddEntry(hd,"Data","lep")
            for k,h in comps: leg.AddEntry(h,SAMPLES[k][0],"l")
            leg.AddEntry(total,"Local fit sum","l"); leg.Draw()
            tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(.040)
            tx.DrawLatex(.17,.92,f"{det}, Candidate {cname}: {lo:g}<p_{{probe}}<{hi:g} GeV")
            tx2=ROOT.TLatex(); tx2.SetNDC(); tx2.SetTextSize(.030)
            tx2.DrawLatex(.17,.855,f"f_{{#pi^{{0}}}}={fit['pi0_fraction']:.3f},  #chi^{{2}}/ndf={fit['chi2']:.1f}/{fit['ndf']}")
            keep += [hd,total,leg,tx,tx2]+[x[1] for x in comps]
        c.SaveAs(str(p/f"04_{det}_candidateB_local_mx2epg_fits.png"))

    # Nominal pi0 fraction and inferred pi0 tag yield.  These are the two
    # quantities Stage 3 will need; the latter is the actual efficiency denominator.
    c=ROOT.TCanvas("clocal_summary","",1500,900); c.Divide(2,1); keep=[]
    for ipad,det in enumerate(("FT","FD"),1):
        pad=c.cd(ipad); setup_pad(pad)
        h=ROOT.TH1D(f"hlocpur_{det}","",4,0,4); h.SetDirectory(0); h.SetMarkerStyle(20); h.SetLineWidth(2)
        for i,ib in enumerate(range(4),1):
            lo,hi=PROBE_P_BINS[ib]
            h.GetXaxis().SetBinLabel(i,f"{lo:g}-{hi:g}")
            h.SetBinContent(i,localfits[(det,cname,ib)]["pi0_fraction"])
        h.SetMinimum(0); h.SetMaximum(1.05)
        h.GetXaxis().SetTitle("Predicted probe momentum (GeV)")
        h.GetYaxis().SetTitle("Locally fitted #pi^{0} fraction")
        h.Draw("E1")
        tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(.045); tx.DrawLatex(.17,.92,f"{det}: Candidate B local purity")
        keep += [h,tx]
    c.SaveAs(str(p/"05_local_pi0_fraction_vs_probe_momentum.png"))


def write_local_purity_tables(results, localfits, outdir):
    p=Path(outdir); p.mkdir(parents=True,exist_ok=True)
    with open(p/"local_candidate_template_fits.txt","w") as f:
        f.write("# Local binned Mx2(ep gamma_tag) template fits.\n")
        f.write("# AAOgen is the exclusive-pi0 component.  Candidate B is nominal; A/C/D are cut variations.\n")
        f.write("# condition is a template-similarity diagnostic: very large values mean component separation is ill-conditioned.\n")
        f.write("detector candidate pmin pmax data Npi0 Nclasdis Ndvcs total pi0_fraction chi2 ndf condition aaoscale clasdisscale dvcsscale\n")
        for det in ("FT","FD"):
            for cname in CANDIDATES:
                for ib,(lo,hi) in enumerate(PROBE_P_BINS):
                    q=localfits[(det,cname,ib)]
                    data=float(results["data"]["fit"][(det,"local",cname,ib,"mx2epg")].GetValue().Integral())
                    f.write(
                        f"{det} {cname} {lo:g} {hi:g} {data:.8g} "
                        f"{q['yield']['aaogen']:.8g} {q['yield']['clasdis']:.8g} {q['yield']['dvcsgen']:.8g} "
                        f"{q['total']:.8g} {q['pi0_fraction']:.8g} {q['chi2']:.8g} {q['ndf']} {q['condition']:.8g} "
                        f"{q['scale']['aaogen']:.8g} {q['scale']['clasdis']:.8g} {q['scale']['dvcsgen']:.8g}\n"
                    )

    with open(p/"candidateB_pi0_tag_denominator.txt","w") as f:
        f.write("# Candidate-B quantities intended to flow into Stage 3.\n")
        f.write("# Npi0_tag is the locally fitted AAOgen contribution in DATA's denominator sample.\n")
        f.write("# This is NOT yet a photon efficiency.\n")
        f.write("detector pmin pmax data Npi0_tag pi0_fraction fit_chi2 fit_ndf condition\n")
        for det in ("FT","FD"):
            for ib,(lo,hi) in enumerate(PROBE_P_BINS):
                q=localfits[(det,NOMINAL_CANDIDATE,ib)]
                data=float(results["data"]["fit"][(det,"local",NOMINAL_CANDIDATE,ib,"mx2epg")].GetValue().Integral())
                f.write(f"{det} {lo:g} {hi:g} {data:.8g} {q['yield']['aaogen']:.8g} "
                        f"{q['pi0_fraction']:.8g} {q['chi2']:.8g} {q['ndf']} {q['condition']:.8g}\n")

    with open(p/"local_fit_guardrails.txt","w") as f:
        f.write("1. The old four-region shared fit is retained only as a diagnostic; its chi2 failure means it is not the nominal purity model.\n")
        f.write("2. Nominal purity now comes from Candidate-B Mx2(ep gamma_tag) fits local in detector region and predicted-probe momentum.\n")
        f.write("3. FT and FD are intentionally allowed to have different fitted mixtures; they are different detector/phase-space populations.\n")
        f.write("4. A/C/D are fitted with the same machinery and will provide a denominator-selection systematic if local closure is acceptable.\n")
        f.write("5. Large template condition numbers flag bins where AAOgen/background shapes are too similar to separate reliably.\n")
        f.write("6. Do not start reconstructed-probe efficiency extraction until Candidate-B local fits and control-region transfer checks are accepted.\n")

def write_normalization_tables(results, fits, outdir):
    p=Path(outdir); p.mkdir(parents=True,exist_ok=True)
    rows,purity=normalization_products(results,fits)
    with open(p/"normalization_fit.txt","w") as f:
        f.write("# Four-region component-normalization DIAGNOSTIC. This is no longer the nominal purity model.\n")
        f.write("# Fits use integrated region yields only; plotted shapes are independent closure tests.\n")
        for name in ("shared","FT","FD"):
            q=fits[name]; f.write(f"fit {name} chi2 {q['chi2']:.6g} ndf {q['ndf']}\n")
            for k in FIT_COMPONENTS: f.write(f"  scale {k} {q['scale'][k]:.10g}\n")
        f.write("\n# nominal post-fit region closure\n# detector region data aaogen clasdis dvcsgen total data_over_total\n")
        for det,fr,data,comp,pred in rows:
            f.write(f"{det} {fr} {data:.8g} {comp['aaogen']:.8g} {comp['clasdis']:.8g} {comp['dvcsgen']:.8g} {pred:.8g} {data/pred if pred else float('nan'):.8g}\n")
    with open(p/"pi0_fraction_vs_probe_momentum.txt","w") as f:
        f.write("# Candidate-B composition from the failed shared-fit diagnostic; do not use as nominal purity.\n")
        f.write("# This is a normalization result, NOT yet a photon efficiency.\n")
        f.write("detector pmin pmax data aaogen clasdis dvcsgen total pi0_fraction data_over_total\n")
        for det,ib,lo,hi,data,comp,pred,pur in purity:
            f.write(f"{det} {lo:g} {hi:g} {data:.8g} {comp['aaogen']:.8g} {comp['clasdis']:.8g} {comp['dvcsgen']:.8g} {pred:.8g} {pur:.8g} {data/pred if pred else float('nan'):.8g}\n")
    with open(p/"normalization_interpretation_guardrails.txt","w") as f:
        f.write("1. The shared FT+FD fit is a diagnostic stress test, not the nominal purity extraction.\n")
        f.write("2. FT-only and FD-only scales are diagnostics for detector-dependent closure.\n")
        f.write("3. The fit uses only four-region yields; shape plots are independent validation.\n")
        f.write("4. Do not call pi0_fraction a measured purity unless yield and shape closure are acceptable.\n")
        f.write("5. Do not introduce reconstructed-probe matching until this denominator normalization is accepted.\n")

def write_tables(results, outdir):
    p=Path(outdir)
    with open(p/"candidate_validation_counts.txt","w") as f:
        f.write("# Unweighted denominator-side counts. No reconstructed-probe requirement.\n")
        f.write("# Candidates: A=(0.24,2,9.2), B=(0.24,2,8), C=(0.22,2,8), D=(0.24,3,8).\n")
        f.write("sample region candidate count retention_vs_region_base\n")
        for key in SAMPLES:
            for region in ("FT","FD"):
                den=float(results[key]["counts"][(region,"base")].GetValue())
                for cname in CANDIDATES:
                    n=float(results[key]["counts"][(region,cname)].GetValue())
                    f.write(f"{key} {region} {cname} {int(n)} {n/den if den else 0:.8f}\n")

    with open(p/"candidate_probe_bin_counts.txt","w") as f:
        f.write("# Counts after each candidate in predicted-probe momentum bins.\n")
        f.write("sample region candidate pmin pmax count\n")
        for key in SAMPLES:
            for region in ("FT","FD"):
                for cname in CANDIDATES:
                    for ib,(lo,hi) in enumerate(PROBE_P_BINS):
                        n=int(results[key]["counts"][(region,cname,"pbin",ib)].GetValue())
                        f.write(f"{key} {region} {cname} {lo:g} {hi:g} {n}\n")

    with open(p/"control_region_counts.txt","w") as f:
        f.write("# Background-enriched control regions; these are NOT efficiency denominators.\n")
        f.write("sample control count\n")
        for key in SAMPLES:
            for cr in CONTROL_REGIONS:
                n=int(results[key]["counts"][("control",cr)].GetValue())
                f.write(f"{key} {cr} {n}\n")

    with open(p/"next_stage_contract.txt","w") as f:
        f.write("Final-validation contract for the next normalization/efficiency stage\n")
        f.write("================================================================\n")
        f.write("1. No reconstructed probe is required in any denominator selection here.\n")
        f.write("2. Candidate B is the nominal denominator; A/C/D are retained as cut-variation systematic checks.\n")
        f.write("3. Pi0 purity must be determined differentially in probe phase space; do not use one global purity number.\n")
        f.write("4. Candidate B plus eta/central/highangle are retained as normalization/control diagnostics.\n")
        f.write("5. Nominal purity is determined locally from Candidate-B Mx2(ep gamma_tag), separately in FT/FD and probe-momentum bins.\n")
        f.write("6. Only after normalization/purity validation should reconstructed-probe matching be introduced for the numerator.\n")
        f.write("7. Final correction convention should remain epsilon_data/epsilon_MC; cross-section correction is its inverse.\n")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=8,help="ROOT worker threads (1-8; default 8)")
    ap.add_argument("--output",default="output_RGA_study")
    args=ap.parse_args(); workers=max(1,min(8,args.workers)); ROOT.EnableImplicitMT(workers)

    outroot=Path(args.output)
    overview_dir=outroot/"00_selection_overview"
    validation_dir=outroot/"01_denominator_validation"
    normalization_dir=outroot/"02_component_normalization"
    efficiency_dir=outroot/"03_efficiency_extraction"
    final_dir=outroot/"04_final_corrections"
    for d in (overview_dir,validation_dir,normalization_dir,efficiency_dir,final_dir): d.mkdir(parents=True,exist_ok=True)

    print("="*78)
    print("Photon-efficiency RGA study - denominator validation + component normalization")
    print(f"ROOT worker threads: {workers}")
    print("No reconstructed-probe requirement.")
    print("Nominal denominator candidate: B = Mx2(ep)<0.24, |coplanarity|>2 deg, theta(tag,X)<8 deg.")
    print("Normalization diagnostic: global four-region AAOgen/CLASDIS/DVCSgen stress test.")
    print("Nominal purity extraction: local Candidate-B Mx2(ep gamma_tag) template fits in FT/FD and p_probe bins.")
    print("Photon-energy plotting range: 0-9 GeV.")
    print("="*78)
    results={}; actions=[]
    for key,(label,pat,color) in SAMPLES.items():
        fs=files_for(pat)
        if not fs: raise RuntimeError(f"No ROOT files for {key}: {pat}")
        print(f"{label:10s}: {len(fs):3d} ROOT files")
        results[key]=book_sample(key,define_columns(make_rdf(fs),key))
        for group in ("counts","survival","selected","slices","control","overview","fit"):
            actions += list(results[key][group].values())
    print(f"Executing {len(actions)} booked actions ...")
    ROOT.RDF.RunGraphs(actions)

    draw_selection_overview(results,overview_dir)
    write_tables(results,validation_dir)
    draw_figure1(results,validation_dir); draw_figure2(results,validation_dir)
    draw_figure3(results,validation_dir); draw_figure4(results,validation_dir)

    # First retain the deliberately over-constrained four-region fit as a
    # diagnostic.  Its failure is useful evidence, but it is no longer the
    # nominal purity model.
    fits=perform_normalization_fit(results)
    write_normalization_tables(results,fits,normalization_dir)
    draw_normalization_summary(results,fits,normalization_dir)

    # Nominal denominator-purity extraction: fit the actual Candidate-B
    # Mx2(ep gamma_tag) spectrum locally in detector region and p_probe.
    localfits=perform_local_purity_fits(results)
    write_local_purity_tables(results,localfits,normalization_dir)
    draw_local_purity_fits(results,localfits,normalization_dir)

    print("Candidate retentions relative to each detector-region optimization base:")
    print(" sample   region   A       B       C       D")
    for key in SAMPLES:
        for region in ("FT","FD"):
            den=float(results[key]["counts"][(region,"base")].GetValue())
            vals=[]
            for cname in CANDIDATES:
                n=float(results[key]["counts"][(region,cname)].GetValue()); vals.append(n/den if den else 0.0)
            print(f" {key:8s} {region:>3s}   " + "  ".join(f"{x:.3f}" for x in vals))

    print("\nGlobal four-region normalization diagnostics (not nominal purity):")
    for name in ("shared","FT","FD"):
        q=fits[name]
        scales=" ".join(f"{k}={q['scale'][k]:.4g}" for k in FIT_COMPONENTS)
        print(f"  {name:6s}: chi2/ndf={q['chi2']:.3g}/{q['ndf']}  {scales}")

    print("\nNominal Candidate-B local purity fits:")
    for det in ("FT","FD"):
        vals=[]
        for ib,(lo,hi) in enumerate(PROBE_P_BINS[:4]):
            q=localfits[(det,NOMINAL_CANDIDATE,ib)]
            vals.append(f"{lo:g}-{hi:g}: fpi0={q['pi0_fraction']:.3f}, chi2/ndf={q['chi2']:.1f}/{q['ndf']}, cond={q['condition']:.1f}")
        print(f"  {det}: " + " | ".join(vals))
    print("\nOutput layout:")
    for d in (overview_dir,validation_dir,normalization_dir,efficiency_dir,final_dir): print(f"  {d}/")
    print("Stage 3/4 remain empty: reconstructed-probe matching waits for acceptance of the local Candidate-B purity/closure results.")


if __name__=="__main__":
    main()
