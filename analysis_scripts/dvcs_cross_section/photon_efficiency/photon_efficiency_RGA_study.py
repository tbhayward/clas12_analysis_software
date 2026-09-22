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
    "tag_p":   ("tag_corr_p",   106,   0.0, 10.6, "Tag photon momentum p_{#gamma,tag} (GeV)"),
    "tag_th":  ("tag_corr_theta",90,   0.0, 45.0, "Tag photon #theta_{#gamma,tag} (deg)"),
    "tag_phi": ("tag_corr_phi", 120, -180., 180., "Tag photon #phi_{#gamma,tag} (deg)"),
    "probe_p": ("probe_corr_p", 106,   0.0, 10.6, "Predicted probe momentum p_{#gamma,probe} (GeV)"),
    "probe_th":("probe_corr_theta",90,  0.0, 45.0, "Predicted probe #theta_{#gamma,probe} (deg)"),
}

# Predicted-probe momentum slices chosen to expose the FT/high-energy issue and
# retained for the later purity fit.  Last bin includes the beam-energy endpoint.
PROBE_P_BINS = [(0.4,2.0), (2.0,4.0), (4.0,6.0), (6.0,8.0), (8.0,10.6)]


def h1(node, name, obs_key, weight=False):
    col, nb, lo, hi, _ = OBS[obs_key]
    model = (name, "", nb, lo, hi)
    return node.Histo1D(model, col, "ana_weight") if weight else node.Histo1D(model, col)


def h2_probe(node, name, obs_key):
    col, nb, lo, hi, _ = OBS[obs_key]
    return node.Histo2D((name, "", 53, 0.0, 10.6, nb, lo, hi), "probe_corr_p", col)


def book_sample(key, df):
    out = {"counts": {}, "shape": {}, "weighted": {}, "slice": {}, "corr": {}}
    base = df.Filter("baseline", "baseline")
    c1 = df.Filter("cut_mx2ep", "Mx2(ep)")
    c2 = df.Filter("cut_mx2eg", "Mx2(e gamma)")
    final = df.Filter("cut_coplanarity", "coplanarity")

    nodes = {"baseline": base, "mx2ep": c1, "mx2eg": c2, "final": final}
    for st, node in nodes.items():
        out["counts"][st] = node.Count()

    # N-1 distributions: all other currently active cuts are applied, but the
    # plotted variable's own cut is omitted.  Mx2(ep gamma_tag), angle(gamma,X),
    # and Delta-t have no active offline cut yet, so they use the final selection.
    nminus1 = {
        "mx2ep": df.Filter("baseline && v_mx2eg>1.4 && fabs(v_copl)<5.7"),
        "mx2eg": df.Filter("baseline && v_mx2ep>-0.231 && v_mx2ep<0.309 && fabs(v_copl)<5.7"),
        "copl":  df.Filter("baseline && v_mx2ep>-0.231 && v_mx2ep<0.309 && v_mx2eg>1.4"),
        "mx2epg": final,
        "angle": final,
        "dt": final,
    }
    for obs, node in nminus1.items():
        out["shape"][("nminus1", obs)] = h1(node, f"h_{key}_n1_{obs}", obs)

    # Tag and predicted-probe kinematics after the current working selection.
    for obs in ("tag_p", "tag_th", "tag_phi", "probe_p", "probe_th"):
        out["shape"][("final", obs)] = h1(final, f"h_{key}_final_{obs}", obs)
        # Weighted copy is intentionally retained now so the next normalization
        # stage does not require redesigning/re-reading the analysis structure.
        out["weighted"][("final", obs)] = h1(final, f"hw_{key}_final_{obs}", obs, True)

    # Strong exclusivity variables in predicted-probe momentum slices.  These
    # are the direct inputs for deciding whether one global pi0 fraction is
    # defensible or whether purity must be fitted versus probe momentum/region.
    for region, rcut in (("FT", "probe_region==0"), ("FD", "probe_region==1")):
        rnode = final.Filter(rcut)
        for obs in ("mx2epg", "angle", "dt"):
            out["corr"][(region, obs)] = h2_probe(rnode, f"h2_{key}_{region}_{obs}", obs)
        for ib, (plo, phi) in enumerate(PROBE_P_BINS):
            sn = rnode.Filter(f"probe_corr_p>={plo} && probe_corr_p<{phi}")
            out["counts"][(region, ib)] = sn.Count()
            for obs in ("mx2epg", "angle"):
                out["slice"][(region, ib, obs)] = h1(sn, f"hs_{key}_{region}_{ib}_{obs}", obs)

    return out


def draw_shape_overlay(pad, results, obs, where, title, cut_lines=None):
    setup_pad(pad)
    hs = []
    ymax = 0.0
    for key in SAMPLES:
        h = unit(clone(results[key]["shape"][(where, obs)]))
        style(h, key)
        hs.append((key, h))
        ymax = max(ymax, h.GetMaximum())

    for i, (key, h) in enumerate(hs):
        h.SetTitle("")
        h.GetXaxis().SetTitle(OBS[obs][4])
        h.GetYaxis().SetTitle("Unit-normalized density")
        h.GetYaxis().SetRangeUser(0.0, max(1e-12, 1.28*ymax))
        opt = "E1" if key == "data" else "HIST"
        if i:
            opt += " SAME"
        h.Draw(opt)

    lines = []
    if cut_lines:
        for x in cut_lines:
            line = ROOT.TLine(x, 0.0, x, 1.12*ymax)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.Draw()
            lines.append(line)

    leg = ROOT.TLegend(0.60, 0.64, 0.91, 0.88)
    leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.032)
    for key, h in hs:
        leg.AddEntry(h, SAMPLES[key][0], "lep" if key == "data" else "l")
    leg.Draw()

    tx = ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043)
    tx.DrawLatex(0.17, 0.92, title)
    return hs, lines, leg


def draw_figure1(results, outdir):
    # The three active offline cuts plus the already-skimmed epgamma missing mass.
    c = ROOT.TCanvas("c01", "", 1500, 1050); c.Divide(2, 2); keep = []
    keep += list(draw_shape_overlay(c.cd(1), results, "mx2ep", "nminus1",
                                    "N-1: M_{X}^{2}(ep)", [-0.231, 0.309]))
    keep += list(draw_shape_overlay(c.cd(2), results, "mx2eg", "nminus1",
                                    "N-1: M_{X}^{2}(e#gamma_{tag})", [1.4]))
    keep += list(draw_shape_overlay(c.cd(3), results, "copl", "nminus1",
                                    "N-1: tag-proton coplanarity", [-5.7, 5.7]))
    keep += list(draw_shape_overlay(c.cd(4), results, "mx2epg", "nminus1",
                                    "After current cuts: M_{X}^{2}(ep#gamma_{tag})", None))
    c.SaveAs(str(Path(outdir)/"01_tag_exclusivity_shapes.png"))


def draw_figure2(results, outdir):
    c = ROOT.TCanvas("c02", "", 1500, 1050); c.Divide(2, 2); keep = []
    keep += list(draw_shape_overlay(c.cd(1), results, "tag_p", "final",
                                    "Tag photon momentum after current cuts"))
    keep += list(draw_shape_overlay(c.cd(2), results, "tag_th", "final",
                                    "Tag photon polar angle after current cuts"))
    keep += list(draw_shape_overlay(c.cd(3), results, "tag_phi", "final",
                                    "Tag photon azimuth after current cuts"))
    keep += list(draw_shape_overlay(c.cd(4), results, "angle", "nminus1",
                                    "#theta(#gamma_{tag},X) after current cuts", [9.2]))
    c.SaveAs(str(Path(outdir)/"02_tag_sample_kinematics.png"))


def draw_figure3(results, outdir):
    c = ROOT.TCanvas("c03", "", 1500, 1050); c.Divide(2, 2); keep = []
    keep += list(draw_shape_overlay(c.cd(1), results, "probe_p", "final",
                                    "Predicted probe momentum"))
    keep += list(draw_shape_overlay(c.cd(2), results, "probe_th", "final",
                                    "Predicted probe polar angle", [5.5]))
    keep += list(draw_shape_overlay(c.cd(3), results, "dt", "nminus1",
                                    "#Delta t diagnostic after current cuts"))
    keep += list(draw_shape_overlay(c.cd(4), results, "mx2epg", "nminus1",
                                    "Loose-skim M_{X}^{2}(ep#gamma_{tag}) window", [-0.25, 0.25]))
    c.SaveAs(str(Path(outdir)/"03_probe_and_remaining_diagnostics.png"))


def draw_slice_panel(pad, results, region, ib, obs):
    setup_pad(pad)
    plo, phi = PROBE_P_BINS[ib]
    hs=[]; ymax=0.0
    for key in SAMPLES:
        h=unit(clone(results[key]["slice"][(region,ib,obs)]))
        style(h,key); hs.append((key,h)); ymax=max(ymax,h.GetMaximum())
    for i,(key,h) in enumerate(hs):
        h.SetTitle(""); h.GetXaxis().SetTitle(OBS[obs][4]); h.GetYaxis().SetTitle("Unit-normalized density")
        h.GetYaxis().SetRangeUser(0,max(1e-12,1.28*ymax))
        opt=("E1" if key=="data" else "HIST") + (" SAME" if i else "")
        h.Draw(opt)
    leg=ROOT.TLegend(0.60,0.64,0.91,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.030)
    for key,h in hs: leg.AddEntry(h,SAMPLES[key][0],"lep" if key=="data" else "l")
    leg.Draw()
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.041)
    tx.DrawLatex(0.17,0.92,f"{region}: {plo:g} < p_{{probe}} < {phi:g} GeV")
    return hs,leg


def draw_figure4(results, outdir):
    # FT is the immediate concern.  Five p_probe slices show whether the tag
    # sample's exclusivity shape changes as we approach the DVCS-like high-E region.
    c=ROOT.TCanvas("c04","",1800,1050); c.Divide(3,2); keep=[]
    for ib in range(len(PROBE_P_BINS)):
        keep += list(draw_slice_panel(c.cd(ib+1),results,"FT",ib,"mx2epg"))
    p=c.cd(6); setup_pad(p); p.Clear()
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.050)
    tx.DrawLatex(0.12,0.83,"FT tag-sample composition diagnostic")
    tx.SetTextSize(0.037)
    tx.DrawLatex(0.12,0.69,"M_{X}^{2}(ep#gamma_{tag}) in predicted-probe")
    tx.DrawLatex(0.12,0.62,"momentum slices after current working cuts.")
    tx.DrawLatex(0.12,0.48,"No #theta(#gamma_{tag},X) or #Delta t cut.")
    tx.DrawLatex(0.12,0.34,"Use this before choosing the purity-fit model.")
    c.SaveAs(str(Path(outdir)/"04_FT_tag_shape_vs_probe_momentum.png"))


def draw_2d_pad(pad, results, key, region, obs):
    setup_pad(pad, 0.16)
    h=clone(results[key]["corr"][(region,obs)])
    h.SetTitle("")
    h.GetXaxis().SetTitle("Predicted probe momentum p_{#gamma,probe} (GeV)")
    h.GetYaxis().SetTitle(OBS[obs][4])
    h.Draw("COLZ")
    tx=ROOT.TLatex(); tx.SetNDC(); tx.SetTextSize(0.043)
    tx.DrawLatex(0.17,0.92,f"{SAMPLES[key][0]} - {region}")
    return h


def draw_figure5(results, outdir):
    # One variable per canvas would be clearer, but we intentionally keep the
    # first pass to one 2x2 summary: all four physical samples, same FT observable.
    c=ROOT.TCanvas("c05","",1500,1050); c.Divide(2,2); keep=[]
    for i,key in enumerate(SAMPLES,1):
        keep.append(draw_2d_pad(c.cd(i),results,key,"FT","angle"))
    c.SaveAs(str(Path(outdir)/"05_FT_angle_gX_vs_probe_momentum.png"))


def write_tables(results, outdir):
    p=Path(outdir)
    with open(p/"cutflow.txt","w") as f:
        f.write("# Counts are unweighted hypothesis counts.\n")
        f.write("# The input skim already imposed -1<Mx2(ep)<2 and -0.25<Mx2(epgamma_tag)<0.25 GeV^2.\n")
        f.write("sample baseline mx2ep mx2eg coplanarity\n")
        for key in SAMPLES:
            vals=[int(results[key]["counts"][st].GetValue()) for st in ("baseline","mx2ep","mx2eg","final")]
            f.write(key+" "+" ".join(str(v) for v in vals)+"\n")

    with open(p/"probe_slice_counts.txt","w") as f:
        f.write("# Unweighted counts after current working cuts; useful for judging fit feasibility.\n")
        f.write("sample region pmin pmax count\n")
        for key in SAMPLES:
            for region in ("FT","FD"):
                for ib,(plo,phi) in enumerate(PROBE_P_BINS):
                    n=int(results[key]["counts"][(region,ib)].GetValue())
                    f.write(f"{key} {region} {plo:g} {phi:g} {n}\n")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--workers",type=int,default=8,help="ROOT worker threads (1-8; default 8)")
    ap.add_argument("--output",default="output_RGA_study")
    args=ap.parse_args()
    workers=max(1,min(8,args.workers))
    ROOT.EnableImplicitMT(workers)
    Path(args.output).mkdir(parents=True,exist_ok=True)

    print("="*76)
    print("Photon-efficiency RGA study - tag-sample characterization")
    print(f"ROOT worker threads: {workers}")
    print("No reconstructed-probe requirement.")
    print("Current offline cuts: Mx2(ep), Mx2(e gamma_tag), coplanarity.")
    print("Mx2(ep gamma_tag), angle(gamma_tag,X), and Delta-t remain diagnostics.")
    print("="*76)

    results={}; actions=[]
    for key,(label,pat,color) in SAMPLES.items():
        fs=files_for(pat)
        if not fs:
            raise RuntimeError(f"No ROOT files for {key}: {pat}")
        print(f"{label:10s}: {len(fs):3d} ROOT files")
        df=define_columns(make_rdf(fs),key)
        results[key]=book_sample(key,df)
        for group in ("counts","shape","weighted","slice","corr"):
            actions += list(results[key][group].values())

    print(f"Executing {len(actions)} booked actions ...")
    ROOT.RDF.RunGraphs(actions)

    write_tables(results,args.output)
    draw_figure1(results,args.output)
    draw_figure2(results,args.output)
    draw_figure3(results,args.output)
    draw_figure4(results,args.output)
    draw_figure5(results,args.output)

    print(f"Done. Five focused figures + two text tables written to: {args.output}")


if __name__=="__main__":
    main()
