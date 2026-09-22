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


def book_sample(key, df):
    """Book the final cut-validation stage in one RDF graph per sample.

    No reconstructed-probe requirement is introduced.  Everything here is
    denominator-side validation: candidate stability, probe-phase-space
    preservation, normalization-variable shapes, and background control regions.
    """
    out = {"counts": {}, "survival": {}, "selected": {}, "slices": {}, "control": {}}
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

            # These are the observables most useful for the next normalization
            # stage.  Keep them for every candidate even though only B is drawn
            # in the compact summary figures.
            for obs in ("mx2epg", "tag_p", "tag_th", "tag_phi", "dt", "probe_p", "probe_th"):
                out["selected"][(region, cname, obs)] = h1(
                    node, f"hsel_{key}_{region}_{cname}_{obs}", obs
                )

            # Predicted-probe slices: essential because the pi0 fraction must
            # ultimately be differential rather than a single global number.
            for ib, (plo, phi) in enumerate(PROBE_P_BINS):
                ps = node.Filter(f"probe_corr_p>={plo} && probe_corr_p<{phi}")
                out["counts"][(region, cname, "pbin", ib)] = ps.Count()
                out["slices"][(region, cname, ib, "mx2epg")] = h1(
                    ps, f"hsl_{key}_{region}_{cname}_{ib}_mx2epg", "mx2epg"
                )
                out["slices"][(region, cname, ib, "dt")] = h1(
                    ps, f"hsl_{key}_{region}_{cname}_{ib}_dt", "dt"
                )

    # Control regions use the same optimization base and are not restricted to
    # FT/FD: they are intended to constrain/validate component normalization.
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
    keep += list(draw_overlay(c.cd(2),results,lambda r,k:r["control"][("central","copl")],"copl",CONTROL_REGIONS["central"][1],[-2,2]))
    keep += list(draw_overlay(c.cd(3),results,lambda r,k:r["control"][("highangle","angle")],"angle",CONTROL_REGIONS["highangle"][1],[12,25]))
    cand=REFERENCE_CANDIDATE
    keep += list(draw_overlay(c.cd(4),results,lambda r,k:r["selected"][("FT",cand,"dt")],"dt",f"Candidate {cand}, FT: #Delta t validation"))
    c.SaveAs(str(Path(outdir)/"04_background_control_regions.png"))


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
        f.write("2. Candidate definitions A-D are centralized in CANDIDATES and must be reused unchanged for normalization tests.\n")
        f.write("3. Pi0 purity must be determined differentially in probe phase space; do not use one global purity number.\n")
        f.write("4. Control regions eta/central/highangle are validation regions, never denominator signal regions.\n")
        f.write("5. Component normalization must be established before interpreting DATA as a pi0 percentage.\n")
        f.write("6. Only after normalization/purity validation should reconstructed-probe matching be introduced for the numerator.\n")
        f.write("7. Final correction convention should remain epsilon_data/epsilon_MC; cross-section correction is its inverse.\n")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=8,help="ROOT worker threads (1-8; default 8)")
    ap.add_argument("--output",default="output_RGA_study")
    args=ap.parse_args(); workers=max(1,min(8,args.workers)); ROOT.EnableImplicitMT(workers)
    Path(args.output).mkdir(parents=True,exist_ok=True)
    print("="*78)
    print("Photon-efficiency RGA study - final exclusive-pi0 denominator validation")
    print(f"ROOT worker threads: {workers}")
    print("No reconstructed-probe requirement.")
    print("Common base: Mx2(ep)>-0.231 and Mx2(e gamma_tag)>1.4.")
    print("Validating candidates A-D, control regions, and probe-energy dependence.")
    print(f"Reference candidate for compact shape plots: {REFERENCE_CANDIDATE}")
    print("Photon-energy plotting range: 0-9 GeV.")
    print("="*78)
    results={}; actions=[]
    for key,(label,pat,color) in SAMPLES.items():
        fs=files_for(pat)
        if not fs: raise RuntimeError(f"No ROOT files for {key}: {pat}")
        print(f"{label:10s}: {len(fs):3d} ROOT files")
        results[key]=book_sample(key,define_columns(make_rdf(fs),key))
        for group in ("counts","survival","selected","slices","control"):
            actions += list(results[key][group].values())
    print(f"Executing {len(actions)} booked actions ...")
    ROOT.RDF.RunGraphs(actions)
    write_tables(results,args.output)
    draw_figure1(results,args.output); draw_figure2(results,args.output)
    draw_figure3(results,args.output); draw_figure4(results,args.output)

    print("Candidate retentions relative to each detector-region optimization base:")
    print(" sample   region   A       B       C       D")
    for key in SAMPLES:
        for region in ("FT","FD"):
            den=float(results[key]["counts"][(region,"base")].GetValue())
            vals=[]
            for cname in CANDIDATES:
                n=float(results[key]["counts"][(region,cname)].GetValue())
                vals.append(n/den if den else 0.0)
            print(f" {key:8s} {region:>3s}   " + "  ".join(f"{x:.3f}" for x in vals))
    print("Control-region counts written to control_region_counts.txt.")
    print(f"Done. Four focused figures + validation tables written to: {args.output}")


if __name__=="__main__":
    main()
