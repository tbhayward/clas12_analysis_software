#!/usr/bin/env python3
"""
CLAS12 photon-efficiency analysis -- restart from the raw e'p'gamma1 hypotheses.

Current step
------------
Use every row in the PhotonEfficiency hypothesis tree, requiring W > 2 GeV and angle(e',gamma1) > 8 deg.
Each row is one reconstructed e'p'gamma1 hypothesis produced by the skim.

Exclusivity plots are followed by high- and low-energy normalization studies.\n\nPlot each of the following observables on its own 1x3 canvas:
  1) E_gamma1                     0.4 to 9 GeV
  2) missing energy e'p'gamma1    0 to 9 GeV
  3) Mx2(e'p')                   -0.5 to 1.0 GeV^2
  4) Mx2(e'p'gamma1)             -0.1 to 0.15 GeV^2
  5) Mx2(e'gamma1)                -20 to 20 GeV^2

Pads are cumulative selections, ordered top-left, top-right, bottom-left, bottom-right:
  1) W > 2 GeV and angle(e',gamma1) > 8 deg
  2) additionally Mx2(e'p') < 0.15 GeV^2
  3) additionally -0.05 < Mx2(e'p'gamma1) < 0.05 GeV^2

Samples:
  Data black, DVCSgen green, AAOgen red, CLASDIS blue.

The script retains ROOT implicit multithreading, books all histogram actions
before executing them, automatically discovers completed skim-version-4 ROOT
files, and clears output/ at the beginning of every run.
"""

import argparse
import glob
import os
import shutil
import sys

import ROOT
import numpy as np
from scipy.ndimage import gaussian_filter1d, shift as ndimage_shift
from scipy.optimize import minimize

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

BASE = "/work/clas12/thayward/photon_efficiency/ROOT_trees"
TREE = "PhotonEfficiency"
EVENT_TREE = "PhotonEfficiencyEvents"
EXPECTED_SKIM_VERSION = 4

# Requested draw order / colors.
SAMPLES = [
    ("data", "Data"),
    ("dvcsgen", "DVCSgen"),
    ("aaogen", "AAOgen"),
    ("clasdis", "CLASDIS incl."),
]

COLORS = {
    "data": ROOT.kBlack,
    "dvcsgen": ROOT.kGreen + 2,
    "aaogen": ROOT.kRed + 1,
    "clasdis": ROOT.kBlue + 1,
}


# ---------------------------------------------------------------------------
# C++ helpers callable from RDataFrame
# ---------------------------------------------------------------------------

ROOT.gInterpreter.Declare(r"""
#include <cmath>

double pe_angle_deg(double th1_deg, double ph1_deg, double th2_deg, double ph2_deg) {
    // Saved theta/phi branches are in DEGREES (producer thetaDeg/phiDeg).
    const double d2r = M_PI / 180.0;
    const double th1 = th1_deg*d2r, ph1 = ph1_deg*d2r;
    const double th2 = th2_deg*d2r, ph2 = ph2_deg*d2r;
    const double dot = std::sin(th1)*std::sin(th2)*std::cos(ph1-ph2)
                     + std::cos(th1)*std::cos(th2);
    const double c = std::max(-1.0, std::min(1.0, dot));
    return std::acos(c) * 180.0 / M_PI;
}

double pe_mgg_deg(double p1, double th1_deg, double ph1_deg,
                  double p2, double th2_deg, double ph2_deg) {
    const double d2r = M_PI / 180.0;
    const double th1=th1_deg*d2r, ph1=ph1_deg*d2r;
    const double th2=th2_deg*d2r, ph2=ph2_deg*d2r;
    const double dot=std::sin(th1)*std::sin(th2)*std::cos(ph1-ph2)+std::cos(th1)*std::cos(th2);
    const double c=std::max(-1.0,std::min(1.0,dot));
    const double m2=2.0*p1*p2*(1.0-c);
    return m2>=0.0 ? std::sqrt(m2) : -1.0;
}

double pe_mx2_egamma(double ebeam, double ep, double eth, double eph,
                     double gp, double gth, double gph) {
    const double me = 0.00051099895;
    const double mp = 0.9382720813;
    const double Ee = std::sqrt(ep*ep + me*me);
    const double d2r = M_PI / 180.0;
    eth *= d2r; eph *= d2r; gth *= d2r; gph *= d2r;
    const double ex = ep*std::sin(eth)*std::cos(eph);
    const double ey = ep*std::sin(eth)*std::sin(eph);
    const double ez = ep*std::cos(eth);
    const double gx = gp*std::sin(gth)*std::cos(gph);
    const double gy = gp*std::sin(gth)*std::sin(gph);
    const double gz = gp*std::cos(gth);
    const double Em = ebeam + mp - Ee - gp;
    const double px = -ex - gx;
    const double py = -ey - gy;
    const double pz = ebeam - ez - gz;
    return Em*Em - px*px - py*py - pz*pz;
}

double pe_emiss_epg(double ebeam, double ep, double pp, double gp) {
    const double me = 0.00051099895;
    const double mp = 0.9382720813;

    const double Ee = std::sqrt(ep*ep + me*me);
    const double Ep = std::sqrt(pp*pp + mp*mp);

    // Missing energy for e p -> e' p' gamma1 X:
    // E_miss = E_beam + M_p - E_e' - E_p' - E_gamma1.
    return ebeam + mp - Ee - Ep - gp;
}

ROOT::VecOps::RVec<double> pe_mgg_tag_probe(
        double tag_p, double tag_th, double tag_ph, int tag_index,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<double>& neutral_p,
        const ROOT::VecOps::RVec<double>& neutral_th,
        const ROOT::VecOps::RVec<double>& neutral_ph) {
    ROOT::VecOps::RVec<double> out;
    const auto n = neutral_idx.size();
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (neutral_idx[i] < 0) continue;
        if (neutral_idx[i] == tag_index) continue;
        if (neutral_pid[i] != 22) continue;
        if (!(neutral_p[i] > 0.0)) continue;
        const double dot = std::sin((tag_th)*M_PI/180.0)*std::sin((neutral_th[i])*M_PI/180.0)*std::cos((tag_ph-neutral_ph[i])*M_PI/180.0)
                         + std::cos((tag_th)*M_PI/180.0)*std::cos((neutral_th[i])*M_PI/180.0);
        const double c = std::max(-1.0, std::min(1.0, dot));
        const double m2 = 2.0*tag_p*neutral_p[i]*(1.0-c);
        if (m2 >= 0.0) out.push_back(std::sqrt(m2));
    }
    return out;
}



ROOT::VecOps::RVec<double> pe_mgg_unique_rec_pair(
        double e_th, double e_ph,
        double tag_p, double tag_th, double tag_ph, int tag_index,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<double>& neutral_p,
        const ROOT::VecOps::RVec<double>& neutral_th,
        const ROOT::VecOps::RVec<double>& neutral_ph) {
    // Reconstruct unique REC photon pairs from the hypothesis tree without using X.
    //
    // Every REC photon can appear as a tag on its own hypothesis row.  Requiring
    // probe REC index > tag REC index therefore keeps a given unordered pair only
    // once, rather than filling gamma_i-gamma_j and gamma_j-gamma_i separately.
    //
    // This diagnostic deliberately uses no epgammaX exclusivity cuts and no
    // nearest-to-X information.  Both photons must satisfy the same basic
    // p >= 0.4 GeV and electron-photon opening-angle > 8 deg requirements.
    ROOT::VecOps::RVec<double> out;
    for (size_t i = 0; i < neutral_idx.size(); ++i) {
        if (neutral_idx[i] < 0) continue;
        if (neutral_idx[i] <= tag_index) continue;
        if (neutral_pid[i] != 22) continue;
        if (!(tag_p >= 0.4) || !(neutral_p[i] >= 0.4)) continue;

        const double ce =
            std::sin((e_th)*M_PI/180.0)*std::sin((neutral_th[i])*M_PI/180.0)*std::cos((e_ph-neutral_ph[i])*M_PI/180.0)
            + std::cos((e_th)*M_PI/180.0)*std::cos((neutral_th[i])*M_PI/180.0);
        const double ae =
            std::acos(std::max(-1.0, std::min(1.0, ce))) * 180.0 / M_PI;
        if (!(ae > 8.0)) continue;

        const double dot =
            std::sin((tag_th)*M_PI/180.0)*std::sin((neutral_th[i])*M_PI/180.0)*std::cos((tag_ph-neutral_ph[i])*M_PI/180.0)
            + std::cos((tag_th)*M_PI/180.0)*std::cos((neutral_th[i])*M_PI/180.0);
        const double c = std::max(-1.0, std::min(1.0, dot));
        const double m2 = 2.0*tag_p*neutral_p[i]*(1.0-c);
        if (m2 >= 0.0) out.push_back(std::sqrt(m2));
    }
    return out;
}

ROOT::VecOps::RVec<double> pe_mgg_nearest_pid22(
        double tag_p, double tag_th, double tag_ph, int tag_index,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<double>& neutral_p,
        const ROOT::VecOps::RVec<double>& neutral_th,
        const ROOT::VecOps::RVec<double>& neutral_ph,
        const ROOT::VecOps::RVec<double>& neutral_da) {
    ROOT::VecOps::RVec<double> out;
    int best = -1; double best_da = 1e99;
    for (size_t i=0;i<neutral_idx.size();++i) {
        if (neutral_idx[i] < 0 || neutral_idx[i] == tag_index || neutral_pid[i] != 22) continue;
        if (!(neutral_p[i] > 0.0) || !std::isfinite(neutral_da[i])) continue;
        if (neutral_da[i] < best_da) { best_da=neutral_da[i]; best=(int)i; }
    }
    if (best < 0) return out;
    const double dot = std::sin((tag_th)*M_PI/180.0)*std::sin((neutral_th[best])*M_PI/180.0)*std::cos((tag_ph-neutral_ph[best])*M_PI/180.0)
                     + std::cos((tag_th)*M_PI/180.0)*std::cos((neutral_th[best])*M_PI/180.0);
    const double c=std::max(-1.0,std::min(1.0,dot));
    const double m2=2.0*tag_p*neutral_p[best]*(1.0-c);
    if (m2>=0.0) out.push_back(std::sqrt(m2));
    return out;
}

ROOT::VecOps::RVec<double> pe_mgg_truth_matched_pi0_reco(
        double tag_p, double tag_th, double tag_ph, int tag_index,
        int mc_tag_index, int mc_tag_pid, int mc_tag_parent,
        int mc_probe_index, int mc_probe_pid, int mc_probe_parent,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<double>& neutral_p,
        const ROOT::VecOps::RVec<double>& neutral_th,
        const ROOT::VecOps::RVec<double>& neutral_ph,
        const ROOT::VecOps::RVec<int>& neutral_mc_index,
        const ROOT::VecOps::RVec<int>& neutral_mc_pid) {
    ROOT::VecOps::RVec<double> out;
    // AAOgen sanity check: reconstructed tag and reconstructed partner must map
    // to the two saved generated pi0-photon roles.
    if (mc_tag_index < 0 || mc_probe_index < 0 || mc_tag_index == mc_probe_index) return out;
    if (mc_tag_pid != 22 || mc_probe_pid != 22) return out;
    if (mc_tag_parent != 111 || mc_probe_parent != 111) return out;
    for (size_t i=0;i<neutral_idx.size();++i) {
        if (neutral_idx[i] < 0 || neutral_idx[i] == tag_index || neutral_pid[i] != 22) continue;
        if (neutral_mc_index[i] != mc_probe_index || neutral_mc_pid[i] != 22) continue;
        if (!(neutral_p[i] > 0.0)) continue;
        const double dot=std::sin((tag_th)*M_PI/180.0)*std::sin((neutral_th[i])*M_PI/180.0)*std::cos((tag_ph-neutral_ph[i])*M_PI/180.0)
                        +std::cos((tag_th)*M_PI/180.0)*std::cos((neutral_th[i])*M_PI/180.0);
        const double c=std::max(-1.0,std::min(1.0,dot));
        const double m2=2.0*tag_p*neutral_p[i]*(1.0-c);
        if (m2>=0.0) out.push_back(std::sqrt(m2));
    }
    return out;
}

double pe_mgg_truth_pi0(double p1,double th1,double ph1,int i1,int pid1,int par1,
                         double p2,double th2,double ph2,int i2,int pid2,int par2) {
    if (i1<0 || i2<0 || i1==i2 || pid1!=22 || pid2!=22 || par1!=111 || par2!=111) return -1.0;
    if (!(p1>0.0) || !(p2>0.0)) return -1.0;
    const double dot=std::sin((th1)*M_PI/180.0)*std::sin((th2)*M_PI/180.0)*std::cos((ph1-ph2)*M_PI/180.0)+std::cos((th1)*M_PI/180.0)*std::cos((th2)*M_PI/180.0);
    const double c=std::max(-1.0,std::min(1.0,dot));
    const double m2=2.0*p1*p2*(1.0-c);
    return m2>=0.0 ? std::sqrt(m2) : -1.0;
}
""")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--period", default="fa18_inb")
    parser.add_argument("--output-dir", default="output")
    parser.add_argument(
        "--threads", type=int, default=2,
        help="ROOT implicit-multithreading worker count (default: 2)",
    )
    return parser.parse_args()


def clear_output_directory(path):
    """Delete everything currently inside output/ before making new plots."""
    if os.path.isdir(path):
        for name in os.listdir(path):
            target = os.path.join(path, name)
            if os.path.isdir(target) and not os.path.islink(target):
                shutil.rmtree(target)
            else:
                os.remove(target)
    else:
        os.makedirs(path, exist_ok=True)

    print(f"Cleared output directory: {os.path.abspath(path)}")


def usable_file(path):
    """Accept only nonempty PhotonEfficiency skim-version-4 ROOT files."""
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        return False, "cannot open"

    tree = f.Get(TREE)
    if not tree:
        f.Close()
        return False, f"missing {TREE}"

    nentries = int(tree.GetEntries())
    if nentries == 0:
        f.Close()
        return False, "empty tree"

    if not tree.GetBranch("skim_version"):
        f.Close()
        return False, "missing skim_version"

    tree.GetEntry(0)
    version = int(getattr(tree, "skim_version"))
    f.Close()

    if version != EXPECTED_SKIM_VERSION:
        return False, f"skim_version={version}, expected {EXPECTED_SKIM_VERSION}"

    return True, "ok"


def discover(sample, period):
    pattern = os.path.join(BASE, sample, period, "*_photon_efficiency.root")
    candidates = sorted(glob.glob(pattern))
    good = []

    print(f"\n[{sample}] {pattern}")
    if not candidates:
        print("  no completed ROOT files found")
        return good

    for path in candidates:
        ok, reason = usable_file(path)
        name = os.path.basename(path)
        if ok:
            good.append(path)
            print(f"  USE   {name}")
        else:
            print(f"  SKIP  {name}: {reason}")

    print(f"  -> {len(good)} usable file(s)")
    return good


def make_chain(paths):
    chain = ROOT.TChain(TREE)
    for path in paths:
        chain.Add(path)
    return chain


def make_dataframe(chain):
    """
    Start from every e'p'gamma1 hypothesis saved in PhotonEfficiency.

    The baseline exclusive sample imposes W > 2 GeV and angle(e',gamma1) > 8 deg.
    There is no beta, photon fiducial, additional missing-mass, inferred-probe,
    or gamma2 requirement before the two displayed exclusivity cuts.
    """
    df = ROOT.RDataFrame(chain)

    # The skim currently saves photon corrected quantities as identity copies
    # of the raw reconstructed photon quantities, so tag_corr_p is the saved
    # reconstructed gamma1 energy/momentum used by the analysis.
    df = df.Define("E_gamma1", "tag_corr_p")

    # Use the corrected proton momentum, matching the producer's physics
    # kinematics and the final analysis convention.
    df = df.Define(
        "Emiss_epg",
        "pe_emiss_epg(beam_energy,e_p,p_corr_p,E_gamma1)",
    )

    df = df.Define(
        "Mx2_egamma1",
        "pe_mx2_egamma(beam_energy,e_p,e_theta,e_phi,E_gamma1,tag_corr_theta,tag_corr_phi)",
    )

    df = df.Define(
        "angle_e_gamma1_deg",
        "pe_angle_deg(e_theta,e_phi,tag_corr_theta,tag_corr_phi)",
    )
    return (df
            .Filter("W > 2.0", "W_gt_2_GeV")
            .Filter("angle_e_gamma1_deg > 8.0", "angle_e_gamma1_gt_8_deg"))


def draw_canvases(dfs, output_dir, period):
    """Draw one 1x3 canvas per observable through the cumulative cut sequence."""
    plots = [
        ("E_gamma1", ";E_{#gamma1} (GeV);Unit-normalized entries", 172, 0.4, 9.0, False, "gamma1_energy"),
        ("Emiss_epg", ";E_{miss}(e'p'#gamma1) (GeV);Unit-normalized entries", 180, 0.0, 9.0, True, "missing_energy_epgamma1"),
        ("Mx2_ep", ";M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries", 180, -0.5, 1.0, False, "mx2_ep"),
        ("Mx2_epg_raw", ";M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries", 180, -0.1, 0.15, True, "mx2_epgamma1"),
        ("Mx2_egamma1", ";M^{2}_{X}(e'#gamma1) (GeV^{2});Unit-normalized entries", 150, 0.0, 15.0, False, "mx2_egamma1"),
    ]

    stages = [
        ("", lambda df: df),
        ("M^{2}_{X}(e'p') < 0.13 GeV^{2}",
         lambda df: df.Filter("Mx2_ep < 0.13", "Mx2_ep_lt_0p13")),
        ("-0.05 < M^{2}_{X}(e'p'#gamma1) < 0.05 GeV^{2}",
         lambda df: df.Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", "Mx2_epg_window")),
        ("M^{2}_{X}(e'#gamma1) > 1.0 GeV^{2}",
         lambda df: df.Filter("Mx2_egamma1 > 1.0", "Mx2_egamma1_gt_1p0")),
    ]

    stage_dfs = {}
    for sample, _label in SAMPLES:
        if sample not in dfs:
            continue
        stage_dfs[(sample, 0)] = dfs[sample]
        for istage in range(1, len(stages)):
            stage_dfs[(sample, istage)] = stages[istage][1](stage_dfs[(sample, istage - 1)])

    booked = {}
    actions = []
    unique = str(abs(hash((output_dir, period))))
    for irow, (expr, title, nbins, xmin, xmax, _logy, _slug) in enumerate(plots):
        for icol in range(len(stages)):
            booked[(irow, icol)] = {}
            for sample, label in SAMPLES:
                key = (sample, icol)
                if key not in stage_dfs:
                    continue
                h = stage_dfs[key].Histo1D(
                    (f"h_{irow}_{icol}_{sample}_{unique}", title, nbins, xmin, xmax), expr
                )
                booked[(irow, icol)][sample] = (label, h)
                actions.append(h)

    if actions:
        ROOT.RDF.RunGraphs(actions)

    keep = list(actions)
    written = []
    for irow, (_expr, _title, _nbins, _xmin, _xmax, logy, slug) in enumerate(plots):
        canvas = ROOT.TCanvas(f"c_{slug}_{unique}", "", 2200, 650)
        canvas.Divide(len(stages), 1, 0.002, 0.002)
        keep.append(canvas)

        for icol, (stage_title, _filter) in enumerate(stages):
            pad = canvas.cd(icol + 1)
            pad.SetTicks(1, 1)
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.035)
            pad.SetBottomMargin(0.13)
            pad.SetTopMargin(0.12)
            if logy:
                pad.SetLogy(True)

            legend = ROOT.TLegend(0.55, 0.68, 0.90, 0.88)
            legend.SetBorderSize(0)
            legend.SetFillStyle(0)
            legend.SetTextSize(0.028)
            keep.append(legend)

            histograms = []
            ymax = 0.0
            positive_min = None
            for sample, label in SAMPLES:
                if sample not in booked[(irow, icol)]:
                    continue
                _, result = booked[(irow, icol)][sample]
                hist = result.GetValue()
                hist.SetDirectory(0)
                hist.SetStats(0)
                hist.SetLineColor(COLORS[sample])
                hist.SetLineWidth(3)
                entries = int(round(hist.GetEntries()))
                integral = hist.Integral(1, hist.GetNbinsX())
                if integral > 0.0:
                    hist.Scale(1.0 / integral)
                ymax = max(ymax, hist.GetMaximum())
                if logy:
                    for ibin in range(1, hist.GetNbinsX() + 1):
                        value = hist.GetBinContent(ibin)
                        if value > 0.0 and (positive_min is None or value < positive_min):
                            positive_min = value
                histograms.append((sample, label, hist))
                keep.append(hist)
                legend.AddEntry(hist, f"{label} (N={entries:,})", "l")

            first = True
            for _sample, _label, hist in histograms:
                if logy:
                    ymin = max((positive_min or 1.0e-6) * 0.5, 1.0e-7)
                    hist.SetMinimum(ymin)
                    hist.SetMaximum(5.0 * ymax if ymax > 0.0 else 1.0)
                else:
                    hist.SetMinimum(0.0)
                    hist.SetMaximum(1.25 * ymax if ymax > 0.0 else 1.0)
                hist.GetXaxis().SetTitleSize(0.042)
                hist.GetYaxis().SetTitleSize(0.040)
                hist.GetXaxis().SetLabelSize(0.034)
                hist.GetYaxis().SetLabelSize(0.034)
                hist.GetXaxis().SetTitleOffset(1.05)
                hist.GetYaxis().SetTitleOffset(1.45)
                hist.Draw("HIST" if first else "HIST SAME")
                first = False

            if stage_title:
                title = ROOT.TLatex()
                title.SetNDC(True)
                title.SetTextAlign(22)
                title.SetTextSize(0.046)
                title.DrawLatex(0.52, 0.925, stage_title)
                keep.append(title)
            legend.Draw()

        output_file = os.path.join(output_dir, f"{irow + 1}_{period}_{slug}.png")
        canvas.SaveAs(output_file)
        written.append(output_file)

    return keep, written

def draw_normalization_step1(dfs, output_dir, period):
    """Plot the E_gamma1 > 4 GeV normalization control region after exclusivity cuts."""
    plots = [
        ("E_gamma1", ";E_{#gamma1} (GeV);Unit-normalized entries", 120, 0.4, 9.0, False),
        ("Emiss_epg", ";E_{miss}(e'p'#gamma1) (GeV);Unit-normalized entries", 120, 0.0, 9.0, True),
        ("Mx2_ep", ";M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries", 120, -0.5, 1.0, False),
        ("Mx2_epg_raw", ";M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries", 120, -0.1, 0.15, True),
    ]

    selected = {}
    for sample, _label in SAMPLES:
        if sample not in dfs:
            continue
        selected[sample] = (dfs[sample]
            .Filter("Mx2_ep < 0.13", "norm_Mx2_ep_lt_0p13")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", "norm_Mx2_epg_window")
            .Filter("Mx2_egamma1 > 1.0", "norm_Mx2_egamma1_gt_1p0")
            .Filter("E_gamma1 > 4.0", "norm_Egamma1_gt_4"))

    unique = str(abs(hash((output_dir, period, "normalization_step1"))))
    booked = {}
    actions = []
    for iplot, (expr, title, nbins, xmin, xmax, _logy) in enumerate(plots):
        booked[iplot] = {}
        for sample, label in SAMPLES:
            if sample not in selected:
                continue
            h = selected[sample].Histo1D(
                (f"h_norm1_{iplot}_{sample}_{unique}", title, nbins, xmin, xmax), expr
            )
            booked[iplot][sample] = (label, h)
            actions.append(h)

    if actions:
        ROOT.RDF.RunGraphs(actions)

    canvas = ROOT.TCanvas(f"c_norm1_{unique}", "", 1500, 1100)
    canvas.Divide(2, 2, 0.002, 0.002)
    keep = [canvas] + list(actions)
    canvas.cd()
    canvas_title = ROOT.TLatex()
    canvas_title.SetNDC(True)
    canvas_title.SetTextAlign(22)
    canvas_title.SetTextSize(0.028)
    canvas_title.DrawLatex(0.50, 0.975, "E_{#gamma1} > 4 GeV")
    keep.append(canvas_title)

    for iplot, (_expr, _title, _nbins, _xmin, _xmax, logy) in enumerate(plots):
        pad = canvas.cd(iplot + 1)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.13)
        pad.SetTopMargin(0.11)
        if logy:
            pad.SetLogy(True)

        legend = ROOT.TLegend(0.55, 0.68, 0.90, 0.88)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.030)
        keep.append(legend)

        histograms = []
        ymax = 0.0
        positive_min = None
        for sample, label in SAMPLES:
            if sample not in booked[iplot]:
                continue
            _, result = booked[iplot][sample]
            hist = result.GetValue()
            hist.SetDirectory(0)
            hist.SetStats(0)
            hist.SetLineColor(COLORS[sample])
            hist.SetLineWidth(3)
            entries = int(round(hist.GetEntries()))
            integral = hist.Integral(1, hist.GetNbinsX())
            if integral > 0.0:
                hist.Scale(1.0 / integral)
            ymax = max(ymax, hist.GetMaximum())
            if logy:
                for ibin in range(1, hist.GetNbinsX() + 1):
                    value = hist.GetBinContent(ibin)
                    if value > 0.0 and (positive_min is None or value < positive_min):
                        positive_min = value
            histograms.append(hist)
            keep.append(hist)
            legend.AddEntry(hist, f"{label} (N={entries:,})", "l")

        first = True
        for hist in histograms:
            if logy:
                ymin = max((positive_min or 1.0e-6) * 0.5, 1.0e-7)
                hist.SetMinimum(ymin)
                hist.SetMaximum(5.0 * ymax if ymax > 0.0 else 1.0)
            else:
                hist.SetMinimum(0.0)
                hist.SetMaximum(1.25 * ymax if ymax > 0.0 else 1.0)
            hist.GetXaxis().SetTitleSize(0.045)
            hist.GetYaxis().SetTitleSize(0.043)
            hist.GetXaxis().SetLabelSize(0.036)
            hist.GetYaxis().SetLabelSize(0.036)
            hist.GetXaxis().SetTitleOffset(1.05)
            hist.GetYaxis().SetTitleOffset(1.45)
            hist.Draw("HIST" if first else "HIST SAME")
            first = False
        legend.Draw()

    output_file = os.path.join(output_dir, f"1_{period}_Egamma1_gt_4_GeV.png")
    canvas.SaveAs(output_file)
    return keep, output_file



def fit_two_poisson_templates(data_hist, dvcs_hist, aao_hist):
    """Fit Data_i = A*DVCS_i + B*AAO_i with a binned Poisson likelihood.

    The minimization is a two-parameter Newton iteration with positivity-preserving
    line search.  The returned covariance is the inverse local curvature matrix
    of -log L at the optimum.
    """
    import math

    bins = []
    for ibin in range(1, data_hist.GetNbinsX() + 1):
        d = float(data_hist.GetBinContent(ibin))
        x = float(dvcs_hist.GetBinContent(ibin))
        y = float(aao_hist.GetBinContent(ibin))
        if x > 0.0 or y > 0.0 or d > 0.0:
            bins.append((d, x, y))

    if not bins:
        raise RuntimeError("No populated bins available for the normalization fit")

    sum_d = sum(v[0] for v in bins)
    sum_x = sum(v[1] for v in bins)
    sum_y = sum(v[2] for v in bins)
    scale = sum_d / max(sum_x + sum_y, 1.0)
    A = max(scale, 1.0e-12)
    B = max(scale, 1.0e-12)

    def nll(a, b):
        total = 0.0
        for d, x, y in bins:
            mu = max(a*x + b*y, 1.0e-300)
            total += mu - d*math.log(mu)
        return total

    for _ in range(100):
        gA = gB = hAA = hAB = hBB = 0.0
        for d, x, y in bins:
            mu = max(A*x + B*y, 1.0e-300)
            common = 1.0 - d/mu
            gA += x*common
            gB += y*common
            curv = d/(mu*mu)
            hAA += x*x*curv
            hAB += x*y*curv
            hBB += y*y*curv

        det = hAA*hBB - hAB*hAB
        if det <= 0.0:
            raise RuntimeError("Normalization-fit curvature matrix is singular")
        stepA = ( hBB*gA - hAB*gB) / det
        stepB = (-hAB*gA + hAA*gB) / det

        old = nll(A, B)
        alpha = 1.0
        accepted = False
        while alpha > 1.0e-8:
            trialA = A - alpha*stepA
            trialB = B - alpha*stepB
            if trialA > 0.0 and trialB > 0.0 and nll(trialA, trialB) <= old:
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            break

        rel = max(abs(trialA-A)/max(A,1.0e-12), abs(trialB-B)/max(B,1.0e-12))
        A, B = trialA, trialB
        if rel < 1.0e-10:
            break

    # Recompute local curvature/covariance at the solution.
    hAA = hAB = hBB = 0.0
    for d, x, y in bins:
        mu = max(A*x + B*y, 1.0e-300)
        curv = d/(mu*mu)
        hAA += x*x*curv
        hAB += x*y*curv
        hBB += y*y*curv
    det = hAA*hBB - hAB*hAB
    varA = hBB/det
    varB = hAA/det
    covAB = -hAB/det
    errA = math.sqrt(max(varA, 0.0))
    errB = math.sqrt(max(varB, 0.0))
    corr = covAB/(errA*errB) if errA > 0.0 and errB > 0.0 else 0.0

    # Poisson deviance: useful goodness-of-description diagnostic.
    deviance = 0.0
    used = 0
    for d, x, y in bins:
        mu = A*x + B*y
        if mu <= 0.0:
            continue
        if d > 0.0:
            deviance += 2.0*(mu - d + d*math.log(d/mu))
        else:
            deviance += 2.0*mu
        used += 1
    ndof = max(used - 2, 0)
    return A, B, errA, errB, corr, deviance, ndof


def draw_normalization_fit(dfs, output_dir, period):
    """Fit A*DVCSgen + B*AAOgen to raw Data E_gamma1 counts for E_gamma1 > 4 GeV.

    The same fitted A and B are then applied without refitting to Emiss(epgamma1),
    Mx2(ep), and Mx2(epgamma1).  CLASDIS is deliberately excluded from the fit.
    """
    needed = ("data", "dvcsgen", "aaogen")
    if any(sample not in dfs for sample in needed):
        print("WARNING: normalization fit skipped; Data, DVCSgen, and AAOgen are all required")
        return [], None

    selected = {}
    for sample in needed:
        selected[sample] = (dfs[sample]
            .Filter("Mx2_ep < 0.15", f"fit_{sample}_Mx2_ep_lt_0p15")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", f"fit_{sample}_Mx2_epg_window")
            .Filter("E_gamma1 > 4.0", f"fit_{sample}_Egamma1_gt_4"))

    plots = [
        ("E_gamma1", ";E_{#gamma1} (GeV);Entries", 100, 4.0, 9.0, False, "fit observable"),
        ("Emiss_epg", ";E_{miss}(e'p'#gamma1) (GeV);Entries", 120, 0.0, 9.0, True, "validation"),
        ("Mx2_ep", ";M^{2}_{X}(e'p') (GeV^{2});Entries", 120, -0.5, 1.0, False, "validation"),
        ("Mx2_epg_raw", ";M^{2}_{X}(e'p'#gamma1) (GeV^{2});Entries", 120, -0.1, 0.15, True, "validation"),
    ]

    unique = str(abs(hash((output_dir, period, "normalization_fit"))))
    booked = {}
    actions = []
    for iplot, (expr, title, nbins, xmin, xmax, _logy, _role) in enumerate(plots):
        booked[iplot] = {}
        for sample in needed:
            h = selected[sample].Histo1D(
                (f"h_normfit_{iplot}_{sample}_{unique}", title, nbins, xmin, xmax), expr
            )
            booked[iplot][sample] = h
            actions.append(h)
    ROOT.RDF.RunGraphs(actions)

    raw = {}
    for iplot in range(len(plots)):
        raw[iplot] = {}
        for sample in needed:
            h = booked[iplot][sample].GetValue()
            h.SetDirectory(0)
            raw[iplot][sample] = h

    A, B, errA, errB, corr, dev, ndof = fit_two_poisson_templates(
        raw[0]["data"], raw[0]["dvcsgen"], raw[0]["aaogen"]
    )
    n_dvcs = A * raw[0]["dvcsgen"].Integral(1, raw[0]["dvcsgen"].GetNbinsX())
    n_aao = B * raw[0]["aaogen"].Integral(1, raw[0]["aaogen"].GetNbinsX())
    n_model = n_dvcs + n_aao
    frac_dvcs = n_dvcs/n_model if n_model > 0.0 else 0.0
    frac_aao = n_aao/n_model if n_model > 0.0 else 0.0

    print("\nNormalization fit from raw E_gamma1 counts, E_gamma1 > 4 GeV")
    print("  Model: Data_i = A * DVCSgen_i + B * AAOgen_i")
    print("  CLASDIS is excluded from the fit.")
    print(f"  A = {A:.8g} +/- {errA:.3g}")
    print(f"  B = {B:.8g} +/- {errB:.3g}")
    print(f"  corr(A,B) = {corr:+.4f}")
    print(f"  Poisson deviance / dof = {dev:.2f} / {ndof}" if ndof else f"  Poisson deviance = {dev:.2f}")
    print(f"  fitted DVCSgen yield = {n_dvcs:,.1f} ({100.0*frac_dvcs:.2f}%)")
    print(f"  fitted AAOgen yield  = {n_aao:,.1f} ({100.0*frac_aao:.2f}%)")
    print("  The same A and B are applied to the other three panels without refitting.\n")

    canvas = ROOT.TCanvas(f"c_normfit_{unique}", "", 1500, 1100)
    canvas.Divide(2, 2, 0.002, 0.002)
    keep = [canvas] + actions

    for iplot, (_expr, _title, _nbins, _xmin, _xmax, logy, role) in enumerate(plots):
        pad = canvas.cd(iplot + 1)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.13)
        pad.SetTopMargin(0.12)
        if logy:
            pad.SetLogy(True)

        data = raw[iplot]["data"].Clone(f"h_draw_data_{iplot}_{unique}")
        dvcs = raw[iplot]["dvcsgen"].Clone(f"h_draw_dvcs_{iplot}_{unique}")
        aao = raw[iplot]["aaogen"].Clone(f"h_draw_aao_{iplot}_{unique}")
        total = dvcs.Clone(f"h_draw_total_{iplot}_{unique}")
        dvcs.Scale(A)
        aao.Scale(B)
        total.Reset("ICES")
        total.Add(dvcs)
        total.Add(aao)
        for h in (data, dvcs, aao, total):
            h.SetDirectory(0)
            h.SetStats(0)
        data.SetLineColor(ROOT.kBlack); data.SetLineWidth(3)
        dvcs.SetLineColor(COLORS["dvcsgen"]); dvcs.SetLineWidth(2)
        aao.SetLineColor(COLORS["aaogen"]); aao.SetLineWidth(2)
        total.SetLineColor(ROOT.kMagenta + 2); total.SetLineWidth(4)
        keep.extend([data, dvcs, aao, total])

        ymax = max(data.GetMaximum(), total.GetMaximum(), dvcs.GetMaximum(), aao.GetMaximum())
        if logy:
            positive = []
            for h in (data, total, dvcs, aao):
                positive += [h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1) if h.GetBinContent(i) > 0.0]
            ymin = max((min(positive) if positive else 1.0)*0.5, 0.1)
            data.SetMinimum(ymin)
            data.SetMaximum(max(10.0*ymax, 10.0))
        else:
            data.SetMinimum(0.0)
            data.SetMaximum(1.25*ymax if ymax > 0.0 else 1.0)
        data.GetXaxis().SetTitleSize(0.045)
        data.GetYaxis().SetTitleSize(0.043)
        data.GetXaxis().SetLabelSize(0.036)
        data.GetYaxis().SetLabelSize(0.036)
        data.GetXaxis().SetTitleOffset(1.05)
        data.GetYaxis().SetTitleOffset(1.45)
        data.Draw("E1")
        dvcs.Draw("HIST SAME")
        aao.Draw("HIST SAME")
        total.Draw("HIST SAME")
        data.Draw("E1 SAME")

        leg = ROOT.TLegend(0.47, 0.67, 0.90, 0.88)
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.029)
        leg.AddEntry(data, f"Data (N={int(data.GetEntries()):,})", "lep")
        leg.AddEntry(total, "A#timesDVCSgen + B#timesAAOgen", "l")
        leg.AddEntry(dvcs, f"A#timesDVCSgen (A={A:.4g})", "l")
        leg.AddEntry(aao, f"B#timesAAOgen (B={B:.4g})", "l")
        leg.Draw()
        keep.append(leg)

        label = ROOT.TLatex()
        label.SetNDC(True); label.SetTextAlign(13); label.SetTextSize(0.034)
        label.DrawLatex(0.16, 0.955, "E_{#gamma1} > 4 GeV: " + role)
        keep.append(label)

    output_file = os.path.join(output_dir, f"2_{period}_normalization_fit_Egamma1_gt_4_GeV.png")
    canvas.SaveAs(output_file)
    return keep, output_file



def _th1_to_numpy(hist):
    """Return regular-bin contents from a TH1 as a float numpy array."""
    return np.asarray([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)], dtype=float)


def _morph_1d_counts(counts, shift_bins, sigma_bins):
    """Apply the DVCS-analysis shift+additional-Gaussian-smearing morph to 1D counts.

    The integral is explicitly restored after morphing, so A and B remain pure
    normalization coefficients rather than absorbing morph-induced yield changes.
    """
    x = np.asarray(counts, dtype=float)
    original = float(np.sum(x))
    if original <= 0.0:
        return np.zeros_like(x)
    out = x.copy()
    if sigma_bins > 1.0e-8:
        out = gaussian_filter1d(out, sigma=float(sigma_bins), mode="nearest")
    if abs(shift_bins) > 1.0e-8:
        out = ndimage_shift(out, shift=float(shift_bins), order=1, mode="nearest", prefilter=False)
    out = np.clip(out, 0.0, None)
    total = float(np.sum(out))
    if total <= 0.0:
        return x.copy()
    return out * (original / total)


def fit_two_poisson_templates_common_morph(data_hist, dvcs_hist, aao_hist, seed_A, seed_B):
    """Fit A*DVCS + B*AAO with one common shift/smearing applied to both templates.

    This mirrors the shift+additional-smearing template morph used in the DVCS
    analysis, but deliberately shares the morph between the two components so
    template-shape freedom cannot independently fake a change in composition.
    Shift and smearing are measured in E_gamma1 histogram bins.
    """
    d = _th1_to_numpy(data_hist)
    x0 = _th1_to_numpy(dvcs_hist)
    y0 = _th1_to_numpy(aao_hist)

    def objective(par):
        A, B, shift_bins, sigma_bins = [float(v) for v in par]
        x = _morph_1d_counts(x0, shift_bins, sigma_bins)
        y = _morph_1d_counts(y0, shift_bins, sigma_bins)
        mu = np.clip(A*x + B*y, 1.0e-12, None)
        return float(np.sum(mu - d*np.log(mu)))

    starts = [
        np.asarray([seed_A, seed_B, 0.0, 0.5]),
        np.asarray([seed_A, seed_B, 0.0, 1.5]),
        np.asarray([seed_A, seed_B, -1.0, 1.0]),
        np.asarray([seed_A, seed_B, +1.0, 1.0]),
    ]
    bounds = ((1.0e-12, None), (1.0e-12, None), (-4.0, 4.0), (0.0, 4.0))
    results = [minimize(objective, q, method="L-BFGS-B", bounds=bounds,
                        options={"maxiter": 500, "ftol": 1.0e-12}) for q in starts]
    res = min(results, key=lambda r: float(r.fun))
    A, B, shift_bins, sigma_bins = [float(v) for v in res.x]
    x = _morph_1d_counts(x0, shift_bins, sigma_bins)
    y = _morph_1d_counts(y0, shift_bins, sigma_bins)
    mu = np.clip(A*x + B*y, 1.0e-12, None)

    dev = 0.0
    used = 0
    for di, mui in zip(d, mu):
        if mui <= 0.0:
            continue
        dev += 2.0*(mui - di + di*np.log(di/mui)) if di > 0.0 else 2.0*mui
        used += 1
    ndof = max(used - 4, 0)
    return {
        "success": bool(res.success), "message": str(res.message),
        "A": A, "B": B, "shift_bins": shift_bins, "sigma_bins": sigma_bins,
        "deviance": float(dev), "ndof": int(ndof), "dvcs": x, "aao": y,
    }


def _fill_th1_from_numpy(template_hist, values, name):
    h = template_hist.Clone(name)
    h.Reset("ICES")
    h.SetDirectory(0)
    for i, value in enumerate(values, start=1):
        h.SetBinContent(i, float(value))
    return h


def draw_normalization_morph_comparison(dfs, output_dir, period):
    """Compare nominal and common-morphed E_gamma1 normalization fits.

    The morph is fitted ONLY to E_gamma1.  Its A and B are then applied without
    refitting to the other observables.  The E_gamma1 shift/smearing itself is
    not blindly transferred to those different variables.
    """
    needed = ("data", "dvcsgen", "aaogen")
    if any(sample not in dfs for sample in needed):
        return [], None

    selected = {}
    for sample in needed:
        selected[sample] = (dfs[sample]
            .Filter("Mx2_ep < 0.15", f"morph_{sample}_Mx2_ep_lt_0p15")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", f"morph_{sample}_Mx2_epg_window")
            .Filter("E_gamma1 > 4.0", f"morph_{sample}_Egamma1_gt_4"))

    plots = [
        ("E_gamma1", ";E_{#gamma1} (GeV);Entries", 100, 4.0, 9.0, False, "morph fit observable"),
        ("Emiss_epg", ";E_{miss}(e'p'#gamma1) (GeV);Entries", 120, 0.0, 9.0, True, "validation: A,B only"),
        ("Mx2_ep", ";M^{2}_{X}(e'p') (GeV^{2});Entries", 120, -0.5, 1.0, False, "validation: A,B only"),
        ("Mx2_epg_raw", ";M^{2}_{X}(e'p'#gamma1) (GeV^{2});Entries", 120, -0.1, 0.15, True, "validation: A,B only"),
    ]
    unique = str(abs(hash((output_dir, period, "normalization_morph_comparison"))))
    booked, actions = {}, []
    for ip, (expr, title, nb, lo, hi, _log, _role) in enumerate(plots):
        booked[ip] = {}
        for sample in needed:
            h = selected[sample].Histo1D((f"h_morph_{ip}_{sample}_{unique}", title, nb, lo, hi), expr)
            booked[ip][sample] = h
            actions.append(h)
    ROOT.RDF.RunGraphs(actions)
    raw = {ip: {} for ip in range(len(plots))}
    for ip in raw:
        for sample in needed:
            h = booked[ip][sample].GetValue(); h.SetDirectory(0); raw[ip][sample] = h

    A0, B0, eA0, eB0, corr0, dev0, ndof0 = fit_two_poisson_templates(
        raw[0]["data"], raw[0]["dvcsgen"], raw[0]["aaogen"])
    mf = fit_two_poisson_templates_common_morph(raw[0]["data"], raw[0]["dvcsgen"], raw[0]["aaogen"], A0, B0)
    A1, B1 = mf["A"], mf["B"]
    dA = 100.0*(A1/A0 - 1.0) if A0 else float("nan")
    dB = 100.0*(B1/B0 - 1.0) if B0 else float("nan")
    binw = raw[0]["data"].GetXaxis().GetBinWidth(1)

    print("\nNormalization-template morphing comparison, E_gamma1 > 4 GeV")
    print("  Nominal fit: Data_i = A*DVCSgen_i + B*AAOgen_i")
    print(f"    A_nominal = {A0:.8g} +/- {eA0:.3g}")
    print(f"    B_nominal = {B0:.8g} +/- {eB0:.3g}")
    print(f"    deviance/dof = {dev0:.2f}/{ndof0}")
    print("  Common-morph fit: same shift and additional Gaussian smearing for both templates")
    print(f"    A_morph = {A1:.8g}   change vs nominal = {dA:+.3f}%")
    print(f"    B_morph = {B1:.8g}   change vs nominal = {dB:+.3f}%")
    print(f"    shift = {mf['shift_bins']:+.4f} bins = {mf['shift_bins']*binw:+.5f} GeV")
    print(f"    additional smearing sigma = {mf['sigma_bins']:.4f} bins = {mf['sigma_bins']*binw:.5f} GeV")
    print(f"    deviance/dof = {mf['deviance']:.2f}/{mf['ndof']}")
    print(f"    minimizer success = {mf['success']} ({mf['message']})")
    print("  Validation panels use A_morph and B_morph without refitting and without transferring")
    print("  the E_gamma1 shift/smearing numerically to the other observables.\n")

    canvas = ROOT.TCanvas(f"c_morphcmp_{unique}", "", 1500, 1100)
    canvas.Divide(2, 2, 0.002, 0.002)
    keep = [canvas] + actions
    for ip, (_expr, _title, _nb, _lo, _hi, logy, role) in enumerate(plots):
        pad = canvas.cd(ip+1); pad.SetTicks(1,1); pad.SetLeftMargin(0.14); pad.SetRightMargin(0.04); pad.SetBottomMargin(0.13); pad.SetTopMargin(0.12)
        if logy: pad.SetLogy(True)
        data = raw[ip]["data"].Clone(f"h_mcmp_data_{ip}_{unique}"); data.SetDirectory(0)
        if ip == 0:
            dvcs = _fill_th1_from_numpy(raw[ip]["dvcsgen"], mf["dvcs"], f"h_mcmp_dvcs_{ip}_{unique}")
            aao = _fill_th1_from_numpy(raw[ip]["aaogen"], mf["aao"], f"h_mcmp_aao_{ip}_{unique}")
        else:
            dvcs = raw[ip]["dvcsgen"].Clone(f"h_mcmp_dvcs_{ip}_{unique}"); dvcs.SetDirectory(0)
            aao = raw[ip]["aaogen"].Clone(f"h_mcmp_aao_{ip}_{unique}"); aao.SetDirectory(0)
        dvcs.Scale(A1); aao.Scale(B1)
        total = dvcs.Clone(f"h_mcmp_total_{ip}_{unique}"); total.Add(aao); total.SetDirectory(0)
        for h in (data,dvcs,aao,total): h.SetStats(0)
        data.SetLineColor(ROOT.kBlack); data.SetLineWidth(3)
        dvcs.SetLineColor(COLORS["dvcsgen"]); dvcs.SetLineWidth(2)
        aao.SetLineColor(COLORS["aaogen"]); aao.SetLineWidth(2)
        total.SetLineColor(ROOT.kMagenta+2); total.SetLineWidth(4)
        keep += [data,dvcs,aao,total]
        ymax=max(h.GetMaximum() for h in (data,dvcs,aao,total))
        if logy:
            pos=[h.GetBinContent(i) for h in (data,dvcs,aao,total) for i in range(1,h.GetNbinsX()+1) if h.GetBinContent(i)>0]
            data.SetMinimum(max((min(pos) if pos else 1.0)*0.5,0.1)); data.SetMaximum(max(10*ymax,10.0))
        else:
            data.SetMinimum(0); data.SetMaximum(1.25*ymax if ymax>0 else 1)
        data.GetXaxis().SetTitleSize(0.045); data.GetYaxis().SetTitleSize(0.043); data.GetXaxis().SetLabelSize(0.036); data.GetYaxis().SetLabelSize(0.036); data.GetXaxis().SetTitleOffset(1.05); data.GetYaxis().SetTitleOffset(1.45)
        data.Draw("E1"); dvcs.Draw("HIST SAME"); aao.Draw("HIST SAME"); total.Draw("HIST SAME"); data.Draw("E1 SAME")
        leg=ROOT.TLegend(0.43,0.66,0.90,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.027)
        leg.AddEntry(data,f"Data (N={int(data.GetEntries()):,})","lep")
        leg.AddEntry(total,"Morphed A#timesDVCSgen + B#timesAAOgen" if ip==0 else "A_{morph}#timesDVCSgen + B_{morph}#timesAAOgen","l")
        leg.AddEntry(dvcs,f"DVCSgen (A_{{morph}}={A1:.4g})","l"); leg.AddEntry(aao,f"AAOgen (B_{{morph}}={B1:.4g})","l"); leg.Draw(); keep.append(leg)
        label=ROOT.TLatex(); label.SetNDC(True); label.SetTextAlign(13); label.SetTextSize(0.032); label.DrawLatex(0.16,0.955,"E_{#gamma1} > 4 GeV: "+role); keep.append(label)

    out=os.path.join(output_dir,f"3_{period}_normalization_common_morph_Egamma1_gt_4_GeV.png")
    canvas.SaveAs(out)
    return keep,out

def print_egamma1_survival_scan(dfs):
    """Print exclusive-sample survival versus E_gamma1 threshold for every sample."""
    thresholds = [0.4, 1.0, 2.0, 3.0, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5]

    exclusive = {}
    totals = {}
    scans = {}
    actions = []

    for sample, _label in SAMPLES:
        if sample not in dfs:
            continue
        ex = (dfs[sample]
              .Filter("Mx2_ep < 0.15", f"scan_{sample}_Mx2_ep_lt_0p15")
              .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05",
                      f"scan_{sample}_Mx2_epg_window"))
        exclusive[sample] = ex
        totals[sample] = ex.Count()
        actions.append(totals[sample])
        scans[sample] = {}
        for threshold in thresholds:
            tag = str(threshold).replace('.', 'p')
            handle = ex.Filter(
                f"E_gamma1 > {threshold:.6g}",
                f"scan_{sample}_Egamma1_gt_{tag}",
            ).Count()
            scans[sample][threshold] = handle
            actions.append(handle)

    if actions:
        ROOT.RDF.RunGraphs(actions)

    print("\nE_gamma1 survival scan after the full exclusive selection")
    print("  Denominator = all rows after Mx2(ep) < 0.15 and -0.05 < Mx2(epgamma1) < 0.05")
    print("  The same selected rows feed E_gamma1, Emiss(epgamma1), Mx2(ep), and Mx2(epgamma1),")
    print("  so N and the survival percentage are identical for each of those four distributions.\n")

    header = f"{'E_gamma1 cut':>14s}"
    active = [(sample, label) for sample, label in SAMPLES if sample in scans]
    for _sample, label in active:
        header += f" | {label:^25s}"
    print(header)
    print('-' * len(header))

    for threshold in thresholds:
        row = f"> {threshold:4.1f} GeV   "
        for sample, _label in active:
            total = int(totals[sample].GetValue())
            n = int(scans[sample][threshold].GetValue())
            pct = 100.0 * n / total if total else 0.0
            row += f" | {n:10,d} ({pct:6.2f}%)"
        print(row)




# ---------------------------------------------------------------------------
# v7 normalization workflow
# ---------------------------------------------------------------------------

def _exclusive_df(df, tag):
    """Final epgammaX exclusivity selection used by normalization and probe stages."""
    return (df
            .Filter("Mx2_ep < 0.13", f"{tag}_Mx2_ep_lt_0p13")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", f"{tag}_Mx2_epg_window")
            .Filter("Mx2_egamma1 > 1.0", f"{tag}_Mx2_egamma1_gt_1p0"))


def _book_raw_hists(selected, specs, tag):
    booked, actions = {}, []
    for ip, (expr, title, nb, lo, hi, logy) in enumerate(specs):
        booked[ip] = {}
        for sample, df in selected.items():
            h = df.Histo1D((f"h_{tag}_{ip}_{sample}", title, nb, lo, hi), expr)
            booked[ip][sample] = h
            actions.append(h)
    if actions:
        ROOT.RDF.RunGraphs(actions)
    raw = {ip: {} for ip in range(len(specs))}
    for ip in raw:
        for sample, handle in booked[ip].items():
            h = handle.GetValue(); h.SetDirectory(0); raw[ip][sample] = h
    return raw, actions


def _style_count_panel(pad, data, components, total, labels, logy=False, title_text=None):
    pad.SetTicks(1, 1); pad.SetLeftMargin(0.14); pad.SetRightMargin(0.04)
    pad.SetBottomMargin(0.13); pad.SetTopMargin(0.13)
    if logy: pad.SetLogy(True)
    data.SetStats(0); data.SetLineColor(ROOT.kBlack); data.SetLineWidth(3); data.SetMarkerStyle(20); data.SetMarkerSize(0.55)
    total.SetStats(0); total.SetLineColor(ROOT.kMagenta+2); total.SetLineWidth(4)
    for sample, h in components.items():
        h.SetStats(0); h.SetLineColor(COLORS[sample]); h.SetLineWidth(2)
    allh = [data, total] + list(components.values())
    ymax = max(h.GetMaximum() for h in allh)
    if logy:
        pos = [h.GetBinContent(i) for h in allh for i in range(1, h.GetNbinsX()+1) if h.GetBinContent(i) > 0]
        data.SetMinimum(max((min(pos) if pos else 1.0)*0.5, 0.1)); data.SetMaximum(max(10*ymax, 10.0))
    else:
        data.SetMinimum(0.0); data.SetMaximum(1.25*ymax if ymax > 0 else 1.0)
    data.GetXaxis().SetTitleSize(0.045); data.GetYaxis().SetTitleSize(0.043)
    data.GetXaxis().SetLabelSize(0.036); data.GetYaxis().SetLabelSize(0.036)
    data.GetXaxis().SetTitleOffset(1.05); data.GetYaxis().SetTitleOffset(1.45)
    data.Draw("E1")
    for h in components.values(): h.Draw("HIST SAME")
    total.Draw("HIST SAME"); data.Draw("E1 SAME")
    leg = ROOT.TLegend(0.40, 0.66, 0.89, 0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.027)
    leg.AddEntry(data, f"Data (N={int(data.GetEntries()):,})", "lep")
    leg.AddEntry(total, labels["total"], "l")
    for sample in ("dvcsgen", "aaogen", "clasdis"):
        if sample in components: leg.AddEntry(components[sample], labels[sample], "l")
    leg.Draw()
    txt = None
    if title_text:
        txt = ROOT.TLatex(); txt.SetNDC(True); txt.SetTextAlign(13); txt.SetTextSize(0.034)
        txt.DrawLatex(0.16, 0.965, title_text)
    return [leg] + ([txt] if txt else [])


def draw_high_energy_normalization_v7(dfs, output_dir, period):
    """High-E control region: nominal and common-morph fits in one 2x2 canvas.

    Top row = nominal coefficients; bottom row = common-morph coefficients.
    Left = E_gamma1 fit observable; right = E_miss validation.  Also returns the
    nominal, morph, and arithmetic-mean A/B coefficient sets for the low-E step.
    """
    needed = ("data", "dvcsgen", "aaogen")
    if any(s not in dfs for s in needed): return [], None, None
    selected = {s: _exclusive_df(dfs[s], f"v7hi_{s}").Filter("E_gamma1 > 4.0", f"v7hi_{s}_Egt4") for s in needed}
    specs = [
        ("E_gamma1", ";E_{#gamma1} (GeV);Entries", 100, 4.0, 9.0, False),
        ("Emiss_epg", ";E_{miss}(e'p'#gamma1) (GeV);Entries", 120, 0.0, 9.0, True),
    ]
    unique = str(abs(hash((period, "v7_high_norm"))))
    raw, actions = _book_raw_hists(selected, specs, f"v7hi_{unique}")
    A0,B0,eA,eB,corr,dev0,nd0 = fit_two_poisson_templates(raw[0]["data"],raw[0]["dvcsgen"],raw[0]["aaogen"])
    mf = fit_two_poisson_templates_common_morph(raw[0]["data"],raw[0]["dvcsgen"],raw[0]["aaogen"],A0,B0)
    A1,B1 = mf["A"],mf["B"]
    Am,Bm = 0.5*(A0+A1),0.5*(B0+B1)
    binw = raw[0]["data"].GetXaxis().GetBinWidth(1)
    coeffs = {"nominal":(A0,B0), "morph":(A1,B1), "mean":(Am,Bm)}
    print("\nHigh-E normalization coefficients, E_gamma1 > 4 GeV")
    print(f"  nominal: A={A0:.8g}, B={B0:.8g}, deviance/dof={dev0:.2f}/{nd0}")
    print(f"  morph:   A={A1:.8g}, B={B1:.8g}, shift={mf['shift_bins']*binw:+.5f} GeV, smear={mf['sigma_bins']*binw:.5f} GeV, deviance/dof={mf['deviance']:.2f}/{mf['ndof']}")
    print(f"  mean:    A={Am:.8g}, B={Bm:.8g}")

    canvas=ROOT.TCanvas(f"c_v7hi_{unique}","",1500,1100); canvas.Divide(2,2,0.002,0.002)
    keep=[canvas]+actions
    for row,(A,B,is_morph,rowname) in enumerate(((A0,B0,False,"No morph"),(A1,B1,True,"Common morph"))):
        for ip in range(2):
            pad=canvas.cd(row*2+ip+1)
            data=raw[ip]["data"].Clone(f"d_v7hi_{row}_{ip}_{unique}"); data.SetDirectory(0)
            if is_morph and ip==0:
                dv=_fill_th1_from_numpy(raw[ip]["dvcsgen"],mf["dvcs"],f"dv_v7hi_{row}_{ip}_{unique}")
                aa=_fill_th1_from_numpy(raw[ip]["aaogen"],mf["aao"],f"aa_v7hi_{row}_{ip}_{unique}")
            else:
                dv=raw[ip]["dvcsgen"].Clone(f"dv_v7hi_{row}_{ip}_{unique}"); dv.SetDirectory(0)
                aa=raw[ip]["aaogen"].Clone(f"aa_v7hi_{row}_{ip}_{unique}"); aa.SetDirectory(0)
            dv.Scale(A); aa.Scale(B); total=dv.Clone(f"tot_v7hi_{row}_{ip}_{unique}"); total.Add(aa); total.SetDirectory(0)
            labels={"total":"A#timesDVCSgen + B#timesAAOgen","dvcsgen":f"DVCSgen (A={A:.4g})","aaogen":f"AAOgen (B={B:.4g})"}
            keep += [data,dv,aa,total] + _style_count_panel(pad,data,{"dvcsgen":dv,"aaogen":aa},total,labels,specs[ip][5],f"E_{{#gamma1}} > 4 GeV: {rowname}")
    out=os.path.join(output_dir,f"2_{period}_highE_normalization_nominal_vs_morph.png"); canvas.SaveAs(out)
    return keep,out,coeffs


def draw_low_energy_shapes_v7(dfs, output_dir, period):
    """Unit-normalized low-E companion to the high-E control-region plot."""
    specs=[
        ("E_gamma1",";E_{#gamma1} (GeV);Unit-normalized entries",100,0.4,4.0,False),
        ("Emiss_epg",";E_{miss}(e'p'#gamma1) (GeV);Unit-normalized entries",120,0.0,9.0,True),
        ("Mx2_ep",";M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries",120,-0.5,1.0,False),
        ("Mx2_epg_raw",";M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries",120,-0.1,0.15,True),
    ]
    selected={s:_exclusive_df(dfs[s],f"v7loShape_{s}").Filter("E_gamma1 < 4.0",f"v7loShape_{s}_Elt4") for s,_ in SAMPLES if s in dfs}
    unique=str(abs(hash((period,"v7_low_shapes")))); raw,actions=_book_raw_hists(selected,specs,f"v7los_{unique}")
    canvas=ROOT.TCanvas(f"c_v7los_{unique}","",1500,1100); canvas.Divide(2,2,0.002,0.002); keep=[canvas]+actions
    canvas.cd()
    canvas_title=ROOT.TLatex(); canvas_title.SetNDC(True); canvas_title.SetTextAlign(22); canvas_title.SetTextSize(0.028)
    canvas_title.DrawLatex(0.50,0.975,"E_{#gamma1} < 4 GeV"); keep.append(canvas_title)
    for ip,spec in enumerate(specs):
        pad=canvas.cd(ip+1); pad.SetTicks(1,1); pad.SetLeftMargin(0.14); pad.SetRightMargin(0.04); pad.SetBottomMargin(0.13); pad.SetTopMargin(0.11)
        if spec[5]: pad.SetLogy(True)
        leg=ROOT.TLegend(0.50,0.68,0.88,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.028); keep.append(leg)
        hs=[]; ymax=0.; pmin=None
        for s,label in SAMPLES:
            if s not in raw[ip]: continue
            h=raw[ip][s].Clone(f"norm_v7los_{ip}_{s}_{unique}"); h.SetDirectory(0); h.SetStats(0); h.SetLineColor(COLORS[s]); h.SetLineWidth(3)
            n=int(round(h.GetEntries())); integ=h.Integral(1,h.GetNbinsX());
            if integ>0: h.Scale(1./integ)
            ymax=max(ymax,h.GetMaximum()); hs.append(h); keep.append(h); leg.AddEntry(h,f"{label} (N={n:,})","l")
            if spec[5]:
                vals=[h.GetBinContent(i) for i in range(1,h.GetNbinsX()+1) if h.GetBinContent(i)>0]
                if vals: pmin=min(vals) if pmin is None else min(pmin,min(vals))
        for j,h in enumerate(hs):
            if spec[5]: h.SetMinimum(max((pmin or 1e-6)*0.5,1e-7)); h.SetMaximum(max(5*ymax,1.0))
            else: h.SetMinimum(0); h.SetMaximum(1.25*ymax if ymax else 1)
            h.GetXaxis().SetTitleSize(0.045); h.GetYaxis().SetTitleSize(0.043); h.GetXaxis().SetLabelSize(0.036); h.GetYaxis().SetLabelSize(0.036); h.GetYaxis().SetTitleOffset(1.45)
            h.Draw("HIST" if j==0 else "HIST SAME")
        leg.Draw()
    out=os.path.join(output_dir,f"1b_{period}_Egamma1_lt_4_GeV.png"); canvas.SaveAs(out); return keep,out


def _poisson_deviance_arrays(d,mu,npar):
    d=np.asarray(d,float); mu=np.clip(np.asarray(mu,float),1e-12,None); dev=0.; used=0
    for di,mi in zip(d,mu):
        dev += 2*(mi-di+di*np.log(di/mi)) if di>0 else 2*mi; used+=1
    return float(dev),max(used-npar,0)


def fit_clasdis_fixed_ab(data_hist,dvcs_hist,aao_hist,clas_hist,A,B):
    d=_th1_to_numpy(data_hist); x=_th1_to_numpy(dvcs_hist); y=_th1_to_numpy(aao_hist); z=_th1_to_numpy(clas_hist)
    base=A*x+B*y
    def obj(q):
        C=max(float(q[0]),1e-12); mu=np.clip(base+C*z,1e-12,None); return float(np.sum(mu-d*np.log(mu)))
    seed=max((np.sum(d)-np.sum(base))/max(np.sum(z),1.),1e-8)
    res=minimize(obj,[seed],method="L-BFGS-B",bounds=((1e-12,None),),options={"maxiter":500,"ftol":1e-12})
    C=float(res.x[0]); mu=base+C*z; dev,nd=_poisson_deviance_arrays(d,mu,1)
    return {"C":C,"deviance":dev,"ndof":nd,"success":bool(res.success)}


def fit_clasdis_fixed_ab_common_morph(data_hist,dvcs_hist,aao_hist,clas_hist,A,B):
    d=_th1_to_numpy(data_hist); x0=_th1_to_numpy(dvcs_hist); y0=_th1_to_numpy(aao_hist); z0=_th1_to_numpy(clas_hist)
    seed=fit_clasdis_fixed_ab(data_hist,dvcs_hist,aao_hist,clas_hist,A,B)["C"]
    def obj(q):
        C,sh,sg=[float(v) for v in q]
        x=_morph_1d_counts(x0,sh,sg); y=_morph_1d_counts(y0,sh,sg); z=_morph_1d_counts(z0,sh,sg)
        mu=np.clip(A*x+B*y+C*z,1e-12,None); return float(np.sum(mu-d*np.log(mu)))
    starts=[[seed,0,.5],[seed,0,1.5],[seed,-1,1],[seed,1,1]]
    bounds=((1e-12,None),(-4.,4.),(0.,4.))
    rr=[minimize(obj,q,method="L-BFGS-B",bounds=bounds,options={"maxiter":500,"ftol":1e-12}) for q in starts]; res=min(rr,key=lambda r:float(r.fun))
    C,sh,sg=[float(v) for v in res.x]; x=_morph_1d_counts(x0,sh,sg); y=_morph_1d_counts(y0,sh,sg); z=_morph_1d_counts(z0,sh,sg); mu=A*x+B*y+C*z
    dev,nd=_poisson_deviance_arrays(d,mu,3)
    return {"C":C,"shift_bins":sh,"sigma_bins":sg,"dvcs":x,"aao":y,"clasdis":z,"deviance":dev,"ndof":nd,"success":bool(res.success)}


def draw_low_energy_normalization_v7(dfs, output_dir, period, coeffs):
    """Low-E normalization in the same 2x2 layout as the high-E study.

    High-E A/B are held fixed.  The top row uses the nominal high-E A/B and
    floats only CLASDIS C.  The bottom row uses the morphed high-E A/B and
    floats C together with a common low-E E_gamma1 shift/smearing.  The mean
    A/B case is still calculated and printed as the third normalization case,
    but is not given a separate canvas.
    """
    needed=("data","dvcsgen","aaogen","clasdis")
    if coeffs is None or any(s not in dfs for s in needed): return [],None
    selected={s:_exclusive_df(dfs[s],f"v8lo_{s}").Filter("E_gamma1 < 4.0",f"v8lo_{s}_Elt4") for s in needed}
    specs=[("E_gamma1",";E_{#gamma1} (GeV);Entries",90,0.4,4.0,False),
           ("Emiss_epg",";E_{miss}(e'p'#gamma1) (GeV);Entries",120,0.,9.,True)]
    unique=str(abs(hash((period,"v8_low_norm")))); raw,actions=_book_raw_hists(selected,specs,f"v8lon_{unique}")
    keep=list(actions)

    A0,B0=coeffs["nominal"]
    f0=fit_clasdis_fixed_ab(raw[0]["data"],raw[0]["dvcsgen"],raw[0]["aaogen"],raw[0]["clasdis"],A0,B0)
    A1,B1=coeffs["morph"]
    f1=fit_clasdis_fixed_ab_common_morph(raw[0]["data"],raw[0]["dvcsgen"],raw[0]["aaogen"],raw[0]["clasdis"],A1,B1)
    Am,Bm=coeffs["mean"]
    fm=fit_clasdis_fixed_ab(raw[0]["data"],raw[0]["dvcsgen"],raw[0]["aaogen"],raw[0]["clasdis"],Am,Bm)

    print("\nLow-E normalization, E_gamma1 < 4 GeV; A and B fixed from high-E fits")
    print(f"  nominal: A={A0:.8g}, B={B0:.8g}, C_CLASDIS={f0['C']:.8g}, deviance/dof={f0['deviance']:.2f}/{f0['ndof']}")
    print(f"  morph  : A={A1:.8g}, B={B1:.8g}, C_CLASDIS={f1['C']:.8g}, deviance/dof={f1['deviance']:.2f}/{f1['ndof']}, shift={f1['shift_bins']:+.4f} bins, smear={f1['sigma_bins']:.4f} bins")
    print(f"  mean   : A={Am:.8g}, B={Bm:.8g}, C_CLASDIS={fm['C']:.8g}, deviance/dof={fm['deviance']:.2f}/{fm['ndof']}")

    canvas=ROOT.TCanvas(f"c_v8low_{unique}","",1500,1100); canvas.Divide(2,2,0.002,0.002); keep.append(canvas)
    rows=(("No morph",A0,B0,f0,False),("Common morph",A1,B1,f1,True))
    for row,(rowname,A,B,fit,ism) in enumerate(rows):
        C=fit["C"]
        for ip in range(2):
            pad=canvas.cd(row*2+ip+1)
            data=raw[ip]["data"].Clone(f"d_v8low_{row}_{ip}_{unique}"); data.SetDirectory(0)
            if ism and ip==0:
                dv=_fill_th1_from_numpy(raw[ip]["dvcsgen"],fit["dvcs"],f"dv_v8low_{row}_{ip}_{unique}")
                aa=_fill_th1_from_numpy(raw[ip]["aaogen"],fit["aao"],f"aa_v8low_{row}_{ip}_{unique}")
                cl=_fill_th1_from_numpy(raw[ip]["clasdis"],fit["clasdis"],f"cl_v8low_{row}_{ip}_{unique}")
            else:
                dv=raw[ip]["dvcsgen"].Clone(f"dv_v8low_{row}_{ip}_{unique}"); dv.SetDirectory(0)
                aa=raw[ip]["aaogen"].Clone(f"aa_v8low_{row}_{ip}_{unique}"); aa.SetDirectory(0)
                cl=raw[ip]["clasdis"].Clone(f"cl_v8low_{row}_{ip}_{unique}"); cl.SetDirectory(0)
            dv.Scale(A); aa.Scale(B); cl.Scale(C)
            total=dv.Clone(f"tot_v8low_{row}_{ip}_{unique}"); total.Add(aa); total.Add(cl); total.SetDirectory(0)
            labels={"total":"A#timesDVCS + B#timesAAO + C#timesCLASDIS",
                    "dvcsgen":f"DVCSgen (A={A:.4g})",
                    "aaogen":f"AAOgen (B={B:.4g})",
                    "clasdis":f"CLASDIS (C={C:.4g})"}
            keep += [data,dv,aa,cl,total] + _style_count_panel(
                pad,data,{"dvcsgen":dv,"aaogen":aa,"clasdis":cl},total,labels,
                specs[ip][5],f"E_{{#gamma1}} < 4 GeV: {rowname}")
    out=os.path.join(output_dir,f"3_{period}_lowE_normalization_nominal_vs_morph.png")
    canvas.SaveAs(out)
    low_coeffs = {
        "nominal": (A0, B0, f0["C"]),
        "morph": (A1, B1, f1["C"]),
        "mean": (Am, Bm, fm["C"]),
    }
    return keep,out,low_coeffs


def draw_full_range_denominator(dfs, output_dir, period, coeffs):
    """Final epgammaX denominator check over the complete E_gamma1 range.

    Each row is one coefficient choice (nominal, morph, arithmetic mean).
    The coefficients alone are applied here; no histogram morph is transferred.
    """
    needed=("data","dvcsgen","aaogen","clasdis")
    if coeffs is None or any(x not in dfs for x in needed): return [],None
    selected={x:_exclusive_df(dfs[x],f"finalden_{x}") for x in needed}
    specs=[("E_gamma1",";E_{#gamma1} (GeV);Entries",170,0.4,8.9,False),
           ("Emiss_epg",";E_{miss}(e'p'#gamma1) (GeV);Entries",180,0.0,9.0,True)]
    unique=str(abs(hash((period,"final_denominator"))))
    raw,actions=_book_raw_hists(selected,specs,f"finalden_{unique}")
    canvas=ROOT.TCanvas(f"c_finalden_{unique}","",1500,1500); canvas.Divide(2,3,0.002,0.002)
    keep=[canvas]+actions
    for row,key in enumerate(("nominal","morph","mean")):
        A,B,C=coeffs[key]
        label={"nominal":"No morph coefficients","morph":"Morph coefficients","mean":"Mean coefficients"}[key]
        for ip,spec in enumerate(specs):
            pad=canvas.cd(row*2+ip+1)
            data=raw[ip]["data"].Clone(f"d_finalden_{row}_{ip}_{unique}"); data.SetDirectory(0)
            dv=raw[ip]["dvcsgen"].Clone(f"dv_finalden_{row}_{ip}_{unique}"); dv.SetDirectory(0); dv.Scale(A)
            aa=raw[ip]["aaogen"].Clone(f"aa_finalden_{row}_{ip}_{unique}"); aa.SetDirectory(0); aa.Scale(B)
            cl=raw[ip]["clasdis"].Clone(f"cl_finalden_{row}_{ip}_{unique}"); cl.SetDirectory(0); cl.Scale(C)
            tot=dv.Clone(f"tot_finalden_{row}_{ip}_{unique}"); tot.SetDirectory(0); tot.Add(aa); tot.Add(cl)
            labels={"total":"A#timesDVCS + B#timesAAO + C#timesCLASDIS",
                    "dvcsgen":f"DVCSgen (A={A:.4g})","aaogen":f"AAOgen (B={B:.4g})","clasdis":f"CLASDIS (C={C:.4g})"}
            keep += [data,dv,aa,cl,tot] + _style_count_panel(pad,data,{"dvcsgen":dv,"aaogen":aa,"clasdis":cl},tot,labels,spec[5],label)
    out=os.path.join(output_dir,f"4_{period}_final_epgammaX_denominator.png"); canvas.SaveAs(out)
    return keep,out


def draw_probe_mgg(dfs, output_dir, period, coeffs):
    """M(gamma_tag gamma_probe) after the final epgammaX selection.

    Data remain in measured counts. MC components are scaled by the arithmetic-
    mean A/B/C normalization coefficients obtained in the immediately preceding
    denominator-normalization step. Individual MC components are shown as thin
    dashed lines; their sum is shown as a separate solid curve.
    """
    if coeffs is None:
        print("WARNING: no denominator normalization coefficients; skipping probe Mgg plot.")
        return [], None

    selected = {}
    for sample, df in dfs.items():
        selected[sample] = (_exclusive_df(df, f"probe_{sample}")
            .Define("Mgg_tag_probe",
                "pe_mgg_tag_probe(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,"
                "neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi)"))

    A, B, C = coeffs["mean"]
    scales = {"data": 1.0, "dvcsgen": A, "aaogen": B, "clasdis": C}
    unique = str(abs(hash((period, "probe_mgg_scaled_components"))))
    booked, actions = {}, []

    for sample, df in selected.items():
        h = df.Histo1D(
            (f"h_probe_mgg_{sample}_{unique}",
             ";M_{#gamma_{tag}#gamma_{probe}} (GeV);Normalized tag-probe combinations",
             160, 0.0, 0.8),
            "Mgg_tag_probe")
        booked[sample] = h
        actions.append(h)

    if actions:
        ROOT.RDF.RunGraphs(actions)

    canvas = ROOT.TCanvas(f"c_probe_mgg_{unique}", "", 1050, 800)
    canvas.SetTicks(1, 1)
    canvas.SetLeftMargin(0.13)
    canvas.SetRightMargin(0.04)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)

    keep = [canvas] + actions
    leg = ROOT.TLegend(0.43, 0.62, 0.89, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.026)
    keep.append(leg)

    drawn = {}
    raw_pairs = {}
    for sample, label in SAMPLES:
        if sample not in booked:
            continue
        h = booked[sample].GetValue().Clone(f"hp_probe_{sample}_{unique}")
        h.SetDirectory(0)
        h.SetStats(0)
        raw_pairs[sample] = int(round(h.GetEntries()))

        if sample != "data":
            h.Scale(scales[sample])
            h.SetLineColor(COLORS[sample])
            h.SetLineWidth(2)
            h.SetLineStyle(2)
        else:
            h.SetLineColor(COLORS["data"])
            h.SetLineWidth(3)
            h.SetLineStyle(1)

        drawn[sample] = h
        keep.append(h)

    # Sum the already-scaled MC components. Use a distinct violet/magenta-family
    # ROOT color so it is visually separate from black Data and the component colors.
    combined = None
    for sample in ("dvcsgen", "aaogen", "clasdis"):
        if sample not in drawn:
            continue
        if combined is None:
            combined = drawn[sample].Clone(f"h_probe_combined_mc_{unique}")
            combined.SetDirectory(0)
        else:
            combined.Add(drawn[sample])

    if combined is not None:
        combined.SetStats(0)
        combined.SetLineColor(ROOT.kMagenta + 2)
        combined.SetLineWidth(4)
        combined.SetLineStyle(1)
        keep.append(combined)

    ymax = 0.0
    if "data" in drawn:
        ymax = max(ymax, drawn["data"].GetMaximum())
    if combined is not None:
        ymax = max(ymax, combined.GetMaximum())
    for sample in ("dvcsgen", "aaogen", "clasdis"):
        if sample in drawn:
            ymax = max(ymax, drawn[sample].GetMaximum())

    # Draw Data first to establish axes, then the combined prediction, then thin
    # component overlays, and finally Data again so it remains visually prominent.
    first_hist = drawn.get("data")
    if first_hist is None:
        first_hist = combined
    if first_hist is None and drawn:
        first_hist = next(iter(drawn.values()))

    if first_hist is not None:
        first_hist.SetMinimum(0.0)
        first_hist.SetMaximum(1.22 * ymax if ymax else 1.0)
        first_hist.GetXaxis().SetTitleSize(0.047)
        first_hist.GetYaxis().SetTitleSize(0.043)
        first_hist.GetXaxis().SetLabelSize(0.038)
        first_hist.GetYaxis().SetLabelSize(0.038)
        first_hist.GetYaxis().SetTitleOffset(1.35)
        first_hist.Draw("HIST")

    if combined is not None and combined is not first_hist:
        combined.Draw("HIST SAME")

    for sample in ("dvcsgen", "aaogen", "clasdis"):
        if sample in drawn and drawn[sample] is not first_hist:
            drawn[sample].Draw("HIST SAME")

    if "data" in drawn and drawn["data"] is not first_hist:
        drawn["data"].Draw("HIST SAME")
    elif "data" in drawn:
        drawn["data"].Draw("HIST SAME")

    if "data" in drawn:
        leg.AddEntry(drawn["data"], f"Data (pairs={raw_pairs['data']:,})", "l")
    if combined is not None:
        leg.AddEntry(combined, "Combined MC", "l")
    if "dvcsgen" in drawn:
        leg.AddEntry(drawn["dvcsgen"], f"DVCSgen (#times{A:.4g})", "l")
    if "aaogen" in drawn:
        leg.AddEntry(drawn["aaogen"], f"AAOgen (#times{B:.4g})", "l")
    if "clasdis" in drawn:
        leg.AddEntry(drawn["clasdis"], f"CLASDIS incl. (#times{C:.4g})", "l")
    leg.Draw()

    print("\nProbe Mgg normalization (mean coefficients from preceding denominator step):")
    print("  Data:     1")
    print(f"  DVCSgen:  A = {A:.8g}")
    print(f"  AAOgen:   B = {B:.8g}")
    print(f"  CLASDIS:  C = {C:.8g}")
    print("  Combined MC = A*DVCSgen + B*AAOgen + C*CLASDIS")

    out = os.path.join(output_dir, f"1_{period}_Mgg_tag_probe.png")
    canvas.SaveAs(out)
    return keep, out


def draw_probe_aaogen_pi0_fit(dfs, output_dir, period, coeffs):
    """Compare single- and double-Gaussian pi0 signal models on normalized AAOgen.

    AAOgen is scaled by the same mean B coefficient used in the preceding probe
    comparison.  The nominal comparison uses 0.05 < Mgg < 0.22 GeV and a
    quadratic background.  The double-Gaussian signal has a common mean and two
    widths; both Gaussian components count as pi0 signal.

    A fit-window scan is also performed for both signal models to quantify how
    stable the extracted normalized pi0 yield is against reasonable changes of
    the fitted mass interval.
    """
    if coeffs is None or "aaogen" not in dfs:
        print("WARNING: AAOgen or normalization coefficients unavailable; skipping pi0 fit.")
        return [], None, None

    _, B, _ = coeffs["mean"]

    df = (_exclusive_df(dfs["aaogen"], "probe_pi0fit_aaogen")
          .Define("Mgg_tag_probe_fit",
              "pe_mgg_tag_probe(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,"
              "neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi)"))

    unique = str(abs(hash((period, "probe_aaogen_pi0_fit_models"))))
    hptr = df.Histo1D(
        (f"h_probe_aaogen_pi0fit_{unique}",
         ";M_{#gamma_{tag}#gamma_{probe}} (GeV);Normalized tag-probe combinations",
         160, 0.0, 0.8),
        "Mgg_tag_probe_fit")
    ROOT.RDF.RunGraphs([hptr])

    h = hptr.GetValue().Clone(f"h_probe_aaogen_pi0fit_scaled_{unique}")
    h.SetDirectory(0)
    h.Sumw2()
    h.Scale(B)
    h.SetStats(0)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.65)
    h.SetMarkerColor(COLORS["aaogen"])
    h.SetLineColor(COLORS["aaogen"])
    h.SetLineWidth(1)
    bin_width = h.GetXaxis().GetBinWidth(1)

    def _seed_background(hist, lo, hi):
        blo = hist.GetXaxis().FindBin(lo + 0.005)
        bhi = hist.GetXaxis().FindBin(hi - 0.005)
        return max(0.0, 0.5 * (hist.GetBinContent(blo) + hist.GetBinContent(bhi)))

    def _fit_single(hist, lo, hi, suffix):
        f = ROOT.TF1(f"f_pi0_single_{unique}_{suffix}", "gaus(0)+pol2(3)", lo, hi)
        peak = max(hist.GetBinContent(hist.GetXaxis().FindBin(0.135)), 1.0e-9)
        bg = _seed_background(hist, lo, hi)
        amp = max(peak - bg, 0.25 * peak, 1.0e-9)
        f.SetParameters(amp, 0.135, 0.012, bg, 0.0, 0.0)
        f.SetParLimits(0, 0.0, max(10.0 * peak, 1.0))
        f.SetParLimits(1, 0.115, 0.155)
        f.SetParLimits(2, 0.003, 0.040)
        r = hist.Fit(f, "SQR0", "", lo, hi)

        sig = ROOT.TF1(f"f_pi0_single_sig_{unique}_{suffix}", "gaus", lo, hi)
        sig.SetParameters(f.GetParameter(0), f.GetParameter(1), abs(f.GetParameter(2)))
        yld = sig.Integral(lo, hi) / bin_width
        return f, r, sig, yld

    def _fit_double(hist, lo, hi, suffix):
        # Two Gaussian components share one mean:
        # [0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2) + pol2(5)
        form = ("[0]*exp(-0.5*((x-[1])/[2])^2)"
                "+[3]*exp(-0.5*((x-[1])/[4])^2)"
                "+[5]+[6]*x+[7]*x*x")
        f = ROOT.TF1(f"f_pi0_double_{unique}_{suffix}", form, lo, hi)
        peak = max(hist.GetBinContent(hist.GetXaxis().FindBin(0.135)), 1.0e-9)
        bg = _seed_background(hist, lo, hi)
        sigheight = max(peak - bg, 0.25 * peak, 1.0e-9)
        f.SetParameters(0.75 * sigheight, 0.135, 0.009,
                        0.25 * sigheight, 0.020,
                        bg, 0.0, 0.0)
        f.SetParLimits(0, 0.0, max(10.0 * peak, 1.0))
        f.SetParLimits(1, 0.115, 0.155)
        f.SetParLimits(2, 0.003, 0.020)
        f.SetParLimits(3, 0.0, max(10.0 * peak, 1.0))
        f.SetParLimits(4, 0.010, 0.060)
        r = hist.Fit(f, "SQR0", "", lo, hi)

        sigform = ("[0]*exp(-0.5*((x-[1])/[2])^2)"
                   "+[3]*exp(-0.5*((x-[1])/[4])^2)")
        sig = ROOT.TF1(f"f_pi0_double_sig_{unique}_{suffix}", sigform, lo, hi)
        sig.SetParameters(f.GetParameter(0), f.GetParameter(1), abs(f.GetParameter(2)),
                          f.GetParameter(3), abs(f.GetParameter(4)))
        yld = sig.Integral(lo, hi) / bin_width
        return f, r, sig, yld

    nominal_lo, nominal_hi = 0.050, 0.220
    f_single, r_single, sig_single, y_single = _fit_single(
        h, nominal_lo, nominal_hi, "nominal")
    f_double, r_double, sig_double, y_double = _fit_double(
        h, nominal_lo, nominal_hi, "nominal")

    # Components for display.
    bg_single = ROOT.TF1(f"f_pi0_single_bg_{unique}", "pol2", nominal_lo, nominal_hi)
    for ip in range(3):
        bg_single.SetParameter(ip, f_single.GetParameter(ip + 3))

    bg_double = ROOT.TF1(
        f"f_pi0_double_bg_{unique}", "[0]+[1]*x+[2]*x*x", nominal_lo, nominal_hi)
    for ip in range(3):
        bg_double.SetParameter(ip, f_double.GetParameter(ip + 5))

    # Window stability scan requested explicitly.
    windows = [
        (0.050, 0.200),
        (0.060, 0.210),
        (0.070, 0.220),
        (0.050, 0.220),
    ]
    scan = []
    scan_keep = []
    for iw, (lo, hi) in enumerate(windows):
        fs, rs, ss, ys = _fit_single(h, lo, hi, f"scan_s_{iw}")
        fd, rd, sd, yd = _fit_double(h, lo, hi, f"scan_d_{iw}")
        scan.append({
            "lo": lo, "hi": hi,
            "single_yield": ys,
            "single_chi2": fs.GetChisquare(),
            "single_ndf": fs.GetNDF(),
            "double_yield": yd,
            "double_chi2": fd.GetChisquare(),
            "double_ndf": fd.GetNDF(),
        })
        scan_keep.extend([fs, rs, ss, fd, rd, sd])

    single_yields = [x["single_yield"] for x in scan]
    double_yields = [x["double_yield"] for x in scan]
    single_span = max(single_yields) - min(single_yields)
    double_span = max(double_yields) - min(double_yields)

    canvas = ROOT.TCanvas(f"c_probe_aaogen_pi0fit_{unique}", "", 1500, 720)
    canvas.Divide(2, 1, 0.002, 0.002)
    keep = [canvas, hptr, h, f_single, r_single, sig_single, bg_single,
            f_double, r_double, sig_double, bg_double] + scan_keep

    def _draw_panel(padnum, total, signal, background, title, yld, extra_lines):
        pad = canvas.cd(padnum)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.10)

        hp = h.Clone(f"h_pi0_panel_{unique}_{padnum}")
        hp.SetDirectory(0)
        hp.GetXaxis().SetRangeUser(0.0, 0.30)
        hp.SetMinimum(0.0)
        hp.SetMaximum(1.28 * max(h.GetMaximum(), total.GetMaximum(nominal_lo, nominal_hi)))
        hp.GetXaxis().SetTitleSize(0.050)
        hp.GetYaxis().SetTitleSize(0.045)
        hp.GetXaxis().SetLabelSize(0.040)
        hp.GetYaxis().SetLabelSize(0.040)
        hp.GetYaxis().SetTitleOffset(1.45)
        hp.Draw("E1")

        total.SetLineColor(ROOT.kMagenta + 2)
        total.SetLineWidth(3)
        total.SetLineStyle(1)
        signal.SetLineColor(ROOT.kBlue + 1)
        signal.SetLineWidth(3)
        signal.SetLineStyle(1)
        background.SetLineColor(ROOT.kGray + 2)
        background.SetLineWidth(3)
        background.SetLineStyle(2)
        background.Draw("SAME")
        signal.Draw("SAME")
        total.Draw("SAME")
        hp.Draw("E1 SAME")

        tt = ROOT.TLatex()
        tt.SetNDC(True)
        tt.SetTextAlign(22)
        tt.SetTextSize(0.038)
        tt.DrawLatex(0.50, 0.955, title)

        info = ROOT.TLatex()
        info.SetNDC(True)
        info.SetTextSize(0.031)
        y0 = 0.88
        for il, line in enumerate(extra_lines):
            info.DrawLatex(0.17, y0 - 0.045 * il, line)

        leg = ROOT.TLegend(0.51, 0.64, 0.92, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.027)
        leg.AddEntry(hp, f"AAOgen #times B ({B:.4g})", "lep")
        leg.AddEntry(total, "Signal + quadratic background", "l")
        leg.AddEntry(signal, "#pi^{0} signal", "l")
        leg.AddEntry(background, "Quadratic background", "l")
        leg.Draw()

        keep.extend([hp, tt, info, leg])

    _draw_panel(
        1, f_single, sig_single, bg_single, "Single-Gaussian signal",
        y_single,
        [
            f"#mu = {f_single.GetParameter(1):.4f} GeV",
            f"#sigma = {abs(f_single.GetParameter(2)):.4f} GeV",
            f"N_{{#pi^{{0}}}} = {y_single:,.0f}",
            f"#chi^{{2}}/ndf = {f_single.GetChisquare():.1f}/{f_single.GetNDF()}",
        ])

    _draw_panel(
        2, f_double, sig_double, bg_double, "Common-mean double-Gaussian signal",
        y_double,
        [
            f"#mu = {f_double.GetParameter(1):.4f} GeV",
            f"#sigma_{{core}} = {abs(f_double.GetParameter(2)):.4f} GeV",
            f"#sigma_{{broad}} = {abs(f_double.GetParameter(4)):.4f} GeV",
            f"N_{{#pi^{{0}}}} = {y_double:,.0f}",
            f"#chi^{{2}}/ndf = {f_double.GetChisquare():.1f}/{f_double.GetNDF()}",
        ])

    out = os.path.join(output_dir, f"2_{period}_AAOgen_Mgg_pi0_fit.png")
    canvas.SaveAs(out)

    # A compact second canvas makes the window dependence visible rather than
    # burying it only in console output.
    cscan = ROOT.TCanvas(f"c_probe_aaogen_pi0scan_{unique}", "", 900, 700)
    cscan.SetTicks(1, 1)
    cscan.SetLeftMargin(0.14)
    cscan.SetRightMargin(0.05)
    cscan.SetBottomMargin(0.14)
    cscan.SetTopMargin(0.08)

    nwin = len(scan)
    gs = ROOT.TGraph(nwin)
    gd = ROOT.TGraph(nwin)
    for i, row in enumerate(scan):
        gs.SetPoint(i, i + 1, row["single_yield"])
        gd.SetPoint(i, i + 1, row["double_yield"])

    gs.SetMarkerStyle(20)
    gs.SetMarkerSize(1.2)
    gs.SetLineWidth(2)
    gs.SetLineColor(ROOT.kBlue + 1)
    gs.SetMarkerColor(ROOT.kBlue + 1)
    gd.SetMarkerStyle(21)
    gd.SetMarkerSize(1.2)
    gd.SetLineWidth(2)
    gd.SetLineColor(ROOT.kMagenta + 2)
    gd.SetMarkerColor(ROOT.kMagenta + 2)

    ymin = 0.985 * min(single_yields + double_yields)
    ymax = 1.015 * max(single_yields + double_yields)
    frame = ROOT.TH1D(f"h_pi0_scan_frame_{unique}",
                      ";Fit window;Normalized fitted N_{#pi^{0}}",
                      nwin, 0.5, nwin + 0.5)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(ymin)
    frame.SetMaximum(ymax)
    for i, row in enumerate(scan, start=1):
        frame.GetXaxis().SetBinLabel(i, f"{row['lo']:.2f}-{row['hi']:.2f}")
    frame.GetXaxis().SetTitleSize(0.047)
    frame.GetYaxis().SetTitleSize(0.043)
    frame.GetXaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetLabelSize(0.038)
    frame.GetYaxis().SetTitleOffset(1.45)
    frame.Draw()
    gs.Draw("LP SAME")
    gd.Draw("LP SAME")

    lscan = ROOT.TLegend(0.55, 0.73, 0.90, 0.88)
    lscan.SetBorderSize(0)
    lscan.SetFillStyle(0)
    lscan.SetTextSize(0.030)
    lscan.AddEntry(gs, "Single Gaussian", "lp")
    lscan.AddEntry(gd, "Double Gaussian", "lp")
    lscan.Draw()

    scan_out = os.path.join(output_dir, f"3_{period}_AAOgen_Mgg_pi0_fit_window_scan.png")
    cscan.SaveAs(scan_out)
    keep.extend([cscan, gs, gd, frame, lscan])

    result = {
        "B": B,
        "nominal_lo": nominal_lo,
        "nominal_hi": nominal_hi,
        "single_yield": y_single,
        "single_chi2": f_single.GetChisquare(),
        "single_ndf": f_single.GetNDF(),
        "double_yield": y_double,
        "double_chi2": f_double.GetChisquare(),
        "double_ndf": f_double.GetNDF(),
        "single_window_span": single_span,
        "double_window_span": double_span,
        "scan": scan,
        "window_scan_output": scan_out,
    }

    print("\\nAAOgen pi0 Mgg signal-model comparison (pre-fiducial-volume step):")
    print(f"  AAOgen normalization B = {B:.8g}")
    print(f"  nominal fit range = {nominal_lo:.3f}--{nominal_hi:.3f} GeV")
    print("  single Gaussian + quadratic background:")
    print(f"    N_pi0 = {y_single:.2f}")
    print(f"    chi2/ndf = {f_single.GetChisquare():.2f}/{f_single.GetNDF()}")
    print("  common-mean double Gaussian + quadratic background:")
    print(f"    N_pi0 = {y_double:.2f}")
    print(f"    chi2/ndf = {f_double.GetChisquare():.2f}/{f_double.GetNDF()}")
    print("\\n  Fit-window stability:")
    for row in scan:
        print(f"    {row['lo']:.3f}--{row['hi']:.3f} GeV:"
              f" single N={row['single_yield']:.2f}, chi2/ndf={row['single_chi2']:.1f}/{row['single_ndf']};"
              f" double N={row['double_yield']:.2f}, chi2/ndf={row['double_chi2']:.1f}/{row['double_ndf']}")
    print(f"  single-Gaussian full yield span = {single_span:.2f}"
          f" ({100.0 * single_span / y_single:.3f}% of nominal)")
    print(f"  double-Gaussian full yield span = {double_span:.2f}"
          f" ({100.0 * double_span / y_double:.3f}% of nominal)")
    print("  NOTE: these are normalized reconstructed AAOgen yields before the")
    print("        measurable probe fiducial-volume restriction.")

    return keep, out, result

def _event_grouped_unique_tag_mgg(df, hist_name, title, tolerance=1.0e-10):
    """Build every unique reconstructed-photon pair for each sequential e+p block.

    For MC we deliberately do NOT use evnum as an event identifier.  The skim
    producer writes all tag-photon hypotheses for one selected electron+proton
    hypothesis consecutively.  We therefore delimit blocks by changes in the
    saved electron and proton kinematics (plus source_file_hash), verify their
    within-block stability, deduplicate photons by tag_rec_index, and fill every
    unordered photon pair exactly once.
    """
    cols = ["source_file_hash", "tag_rec_index",
            "e_p", "e_theta", "e_phi", "e_vz",
            "p_rec_index", "p_raw_p", "p_raw_theta", "p_raw_phi",
            "p_corr_p", "p_corr_theta", "p_corr_phi", "p_vz",
            "tag_corr_p", "tag_corr_theta", "tag_corr_phi"]
    arr = df.AsNumpy(cols)
    h = ROOT.TH1D(hist_name, title, 150, 0.0, 0.30)
    h.SetDirectory(0)
    nrows = len(arr["tag_rec_index"])
    if nrows == 0:
        return h, 0, 0, {}

    sh = np.asarray(arr["source_file_hash"], dtype=np.uint64)
    ri = np.asarray(arr["tag_rec_index"], dtype=np.int64)
    pri = np.asarray(arr["p_rec_index"], dtype=np.int64)
    ep_cols = ["e_p", "e_theta", "e_phi", "e_vz",
               "p_raw_p", "p_raw_theta", "p_raw_phi",
               "p_corr_p", "p_corr_theta", "p_corr_phi", "p_vz"]
    kin = np.column_stack([np.asarray(arr[x], dtype=np.float64) for x in ep_cols])
    pp = np.asarray(arr["tag_corr_p"], dtype=np.float64)
    # Producer stores all theta/phi values in degrees. Convert once before trig.
    tt = np.deg2rad(np.asarray(arr["tag_corr_theta"], dtype=np.float64))
    ph = np.deg2rad(np.asarray(arr["tag_corr_phi"], dtype=np.float64))

    # A new block begins whenever the source file, selected proton REC index,
    # or any saved electron/proton kinematic quantity changes beyond tolerance.
    same_source = sh[1:] == sh[:-1]
    same_proton = pri[1:] == pri[:-1]
    delta = np.max(np.abs(kin[1:] - kin[:-1]), axis=1)
    same_kin = delta <= tolerance
    same = same_source & same_proton & same_kin
    starts = np.r_[0, np.nonzero(~same)[0] + 1]
    stops = np.r_[starts[1:], nrows]

    n_blocks = len(starts)
    n_blocks_ge2 = 0
    n_pairs = 0
    n_duplicate_tag_rows = 0
    max_within_block_spread = 0.0
    block_sizes = []

    for lo, hi in zip(starts, stops):
        # Explicitly verify that the whole block, not merely adjacent rows,
        # carries the same e+p hypothesis.
        if hi - lo > 1:
            spread = float(np.max(np.abs(kin[lo:hi] - kin[lo])))
            max_within_block_spread = max(max_within_block_spread, spread)

        # Preserve first appearance of each valid REC photon in this sequential
        # e+p block. Repeated tag_rec_index rows are not separate photons.
        seen = set()
        idx = []
        for j in range(lo, hi):
            rec = int(ri[j])
            if rec < 0:
                continue
            if rec in seen:
                n_duplicate_tag_rows += 1
                continue
            seen.add(rec)
            idx.append(j)
        block_sizes.append(len(idx))
        if len(idx) < 2:
            continue
        n_blocks_ge2 += 1
        for ia in range(len(idx) - 1):
            i = idx[ia]
            js = np.asarray(idx[ia+1:], dtype=np.int64)
            dot = (np.sin(tt[i]) * np.sin(tt[js]) * np.cos(ph[i] - ph[js])
                   + np.cos(tt[i]) * np.cos(tt[js]))
            m2 = 2.0 * pp[i] * pp[js] * (1.0 - np.clip(dot, -1.0, 1.0))
            for mass in np.sqrt(np.maximum(m2, 0.0)):
                h.Fill(float(mass))
                n_pairs += 1

    bs = np.asarray(block_sizes, dtype=np.int64)
    diagnostics = {
        "rows": nrows,
        "blocks": n_blocks,
        "blocks_ge2": n_blocks_ge2,
        "pairs": n_pairs,
        "duplicate_tag_rows": n_duplicate_tag_rows,
        "max_within_block_spread": max_within_block_spread,
        "block_size_mean": float(np.mean(bs)) if len(bs) else 0.0,
        "block_size_max": int(np.max(bs)) if len(bs) else 0,
        "block_size_1": int(np.count_nonzero(bs == 1)),
        "block_size_2": int(np.count_nonzero(bs == 2)),
        "block_size_ge3": int(np.count_nonzero(bs >= 3)),
        "tolerance": tolerance,
    }
    return h, n_blocks, n_blocks_ge2, diagnostics

def draw_probe_mgg_truth_diagnostic(dfs, output_dir, period):
    """AAOgen event-level reconstructed-photon pair audit.

    MC evnum is deliberately not used. Consecutive rows are grouped only while
    their saved electron and selected-proton kinematics remain identical within
    a tight tolerance. Unique tag_rec_index photons are then paired once.

    The retained neutral arrays are NOT used.  Neither inferred X nor MC truth
    enters the pair construction.
    """
    if "aaogen" not in dfs:
        return [], None

    # dfs["aaogen"] already has W > 2 GeV and angle(e',tag) > 8 deg.  Since
    # every unique photon is recovered through a row where that photon itself
    # is the tag, this means every photon entering the collection satisfies the
    # same >8 degree requirement.
    base = dfs["aaogen"]
    excl = _exclusive_df(dfs["aaogen"], "probe_event_grouped_exclusive")

    unique = str(abs(hash((period, "probe_event_grouped_pair_audit"))))
    hbase, nev_base, nev2_base, diag_base = _event_grouped_unique_tag_mgg(
        base, f"h_mgg_event_grouped_base_{unique}",
        ";M_{#gamma#gamma} (GeV);Unit-normalized unique event-level #gamma#gamma pairs")
    hexcl, nev_excl, nev2_excl, diag_excl = _event_grouped_unique_tag_mgg(
        excl, f"h_mgg_event_grouped_excl_{unique}",
        ";M_{#gamma#gamma} (GeV);Unit-normalized unique event-level #gamma#gamma pairs")

    canvas = ROOT.TCanvas(f"c_event_grouped_pair_audit_{unique}", "", 1050, 900)
    canvas.Divide(1, 2)
    keep = [canvas, hbase, hexcl]

    panels = [
        (hbase, "AAOgen: sequential e+p blocks, unique reconstructed photons",
         nev_base, nev2_base),
        (hexcl, "AAOgen: sequential e+p blocks after final ep#gammaX cuts",
         nev_excl, nev2_excl),
    ]

    for ipad, (h, title, nev, nev2) in enumerate(panels, 1):
        pad = canvas.cd(ipad)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.11)

        h.SetStats(0)
        h.SetLineColor(ROOT.kRed + 1)
        h.SetLineWidth(3)
        npairs = int(round(h.GetEntries()))
        integ = h.Integral(1, h.GetNbinsX())
        if integ > 0:
            h.Scale(1.0 / integ)
        h.SetMinimum(0.0)
        h.SetMaximum(1.25 * h.GetMaximum() if h.GetMaximum() > 0 else 1.0)
        h.GetXaxis().SetTitleSize(0.050)
        h.GetYaxis().SetTitleSize(0.045)
        h.GetXaxis().SetLabelSize(0.040)
        h.GetYaxis().SetLabelSize(0.040)
        h.GetYaxis().SetTitleOffset(1.25)
        h.Draw("HIST")

        line = ROOT.TLine(0.1349768, 0.0, 0.1349768, h.GetMaximum())
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw()

        tex = ROOT.TLatex()
        tex.SetNDC(True)
        tex.SetTextFont(42)
        tex.SetTextSize(0.040)
        tex.DrawLatex(0.13, 0.94, title)
        tex.SetTextSize(0.031)
        tex.DrawLatex(0.64, 0.84, f"events = {nev:,}")
        tex.DrawLatex(0.64, 0.79, f"events with #geq2 #gamma = {nev2:,}")
        tex.DrawLatex(0.64, 0.74, f"pairs = {npairs:,}")
        keep += [line, tex]

    out = os.path.join(output_dir, f"2_{period}_AAOgen_event_grouped_REC_pair_audit.png")
    canvas.SaveAs(out)

    print("\nAAOgen sequential e+p reconstructed-pair audit:")
    print("  MC evnum is NOT used for grouping.")
    print("  Construction: consecutive rows remain in one block only while source_file_hash,")
    print("  p_rec_index, and all saved electron/proton kinematics agree within 1e-10.")
    print("  Within each block, photons are deduplicated by tag_rec_index and every unordered pair is filled once.")
    print("  neutral_[0..4], inferred X, nearest-to-X, and MC truth are NOT used.")
    print(f"  Baseline: {nev_base:,} e+p blocks; {nev2_base:,} with >=2 unique photons; {int(hbase.GetEntries()):,} pairs")
    print(f"    block sizes: N1={diag_base['block_size_1']:,}, N2={diag_base['block_size_2']:,}, N>=3={diag_base['block_size_ge3']:,}, max={diag_base['block_size_max']}")
    print(f"    duplicate tag rows removed={diag_base['duplicate_tag_rows']:,}; max within-block e/p spread={diag_base['max_within_block_spread']:.3e}")
    print(f"  Exclusive: {nev_excl:,} e+p blocks; {nev2_excl:,} with >=2 unique photons; {int(hexcl.GetEntries()):,} pairs")
    print(f"    block sizes: N1={diag_excl['block_size_1']:,}, N2={diag_excl['block_size_2']:,}, N>=3={diag_excl['block_size_ge3']:,}, max={diag_excl['block_size_max']}")
    print(f"    duplicate tag rows removed={diag_excl['duplicate_tag_rows']:,}; max within-block e/p spread={diag_excl['max_within_block_spread']:.3e}")
    print("  EXPECTATION: AAOgen should show a visible pi0 peak near 0.135 GeV if the")
    print("  reconstructed tag photons retained in PhotonEfficiency contain both pi0 daughters.")
    return keep, out


def _mgg_from_p_theta_phi(p1, th1_deg, ph1_deg, p2, th2_deg, ph2_deg):
    """Photon-pair mass from saved skim angles, which are in degrees."""
    th1, ph1, th2, ph2 = np.deg2rad([th1_deg, ph1_deg, th2_deg, ph2_deg])
    dot = (np.sin(th1)*np.sin(th2)*np.cos(ph1-ph2) + np.cos(th1)*np.cos(th2))
    c = max(-1.0, min(1.0, float(dot)))
    return np.sqrt(max(0.0, 2.0*float(p1)*float(p2)*(1.0-c)))

def dump_aaogen_truth_rec_events(files, output_dir, period, max_events=20):
    """AAOgen generator-bookkeeping, angle-unit, and truth->REC diagnostics.

    This intentionally does NOT require gen_n_pi0_photon.  That counter comes
    from MC::Lund, whereas the saved gen_gamma_* arrays come from MC::Particle.
    AAOgen can therefore have two genuine generated photons while the LUND pi0
    counters are zero/missing.  We inspect those two truth systems separately.
    """
    if not files:
        return None
    ch = ROOT.TChain(EVENT_TREE)
    for f in files:
        ch.Add(f)
    nentries = int(ch.GetEntries())
    if nentries == 0:
        print("WARNING: AAOgen PhotonEfficiencyEvents tree is empty/missing.")
        return None

    txt_out = os.path.join(output_dir, f"3_{period}_AAOgen_truth_REC_event_dump.txt")
    png_out = os.path.join(output_dir, f"3_{period}_AAOgen_generator_truth_audit.png")
    branch_names = {b.GetName() for b in ch.GetListOfBranches()}
    needed = ["W", "has_mc_particle", "has_mc_lund", "gen_n_photon", "gen_n_pi0_photon",
              "gen_n_pi0", "gen_gamma_total", "gen_gamma_p", "gen_gamma_theta", "gen_gamma_phi",
              "gen_gamma_index", "gen_gamma_n_rec_matches", "gen_gamma_n_rec_pid22",
              "gen_gamma_best_rec_index", "gen_gamma_best_rec_pid", "gen_gamma_best_rec_p",
              "gen_gamma_best_rec_theta", "gen_gamma_best_rec_phi", "gen_gamma_best_rec_delta_alpha"]
    missing = [x for x in needed if x not in branch_names]
    if missing:
        print("\nERROR: PhotonEfficiencyEvents diagnostic is missing required branches:")
        for x in missing:
            print(f"  - {x}")
        return None

    rdf = ROOT.RDataFrame(ch).Filter("W > 2.0", "W > 2 GeV")
    # Fast aggregate bookkeeping.  These actions run in compiled ROOT code.
    h_nph = rdf.Histo1D((f"h_evt_nph_{period}", ";gen_n_photon;Events", 9, -0.5, 8.5), "gen_n_photon")
    h_ngt = rdf.Histo1D((f"h_evt_ngt_{period}", ";gen_gamma_total (MC::Particle);Events", 9, -0.5, 8.5), "gen_gamma_total")
    h_npi = rdf.Histo1D((f"h_evt_npi_{period}", ";gen_n_pi0 (MC::Lund);Events", 6, -0.5, 5.5), "gen_n_pi0")
    h_npig = rdf.Histo1D((f"h_evt_npig_{period}", ";gen_n_pi0_photon (MC::Lund);Events", 9, -0.5, 8.5), "gen_n_pi0_photon")
    n_w = rdf.Count()
    n_lund = rdf.Filter("has_mc_lund == 1").Count()
    n_part = rdf.Filter("has_mc_particle == 1").Count()
    n_two_mcpart = rdf.Filter("gen_gamma_total == 2").Count()
    n_two_lund = rdf.Filter("gen_n_photon == 2").Count()
    n_two_pi0lund = rdf.Filter("gen_n_pi0_photon == 2").Count()
    two = rdf.Filter("gen_gamma_total == 2") \
        .Define("Mgg_gen_first2", "pe_mgg_deg(gen_gamma_p[0],gen_gamma_theta[0],gen_gamma_phi[0],gen_gamma_p[1],gen_gamma_theta[1],gen_gamma_phi[1])")
    h_mgen = two.Histo1D((f"h_evt_mgen_{period}", ";Generated M_{#gamma#gamma} (GeV);Events", 160, 0.0, 0.30), "Mgg_gen_first2")
    actions=[h_nph,h_ngt,h_npi,h_npig,n_w,n_lund,n_part,n_two_mcpart,n_two_lund,n_two_pi0lund,h_mgen]
    ROOT.RDF.RunGraphs(actions)

    # Plot the bookkeeping and the actual MC::Particle two-photon invariant mass.
    c=ROOT.TCanvas(f"c_gen_truth_audit_{period}","",1200,900); c.Divide(2,2)
    hs=[h_ngt.GetValue(),h_nph.GetValue(),h_npig.GetValue(),h_mgen.GetValue()]
    titles=["Saved generated photons from MC::Particle","Stable photons counted from MC::Lund",
            "Direct #pi^{0} daughter photons counted from MC::Lund","MC::Particle events with exactly two generated photons"]
    keep=[c]
    for i,(h,title) in enumerate(zip(hs,titles),1):
        pad=c.cd(i); pad.SetTicks(1,1); pad.SetLeftMargin(0.13); pad.SetBottomMargin(0.13); pad.SetTopMargin(0.11)
        h.SetStats(0); h.SetLineWidth(3); h.SetTitle(title); h.Draw("HIST")
        keep.append(h)
        if i==4:
            line=ROOT.TLine(0.1349768,0.0,0.1349768,max(1.0,h.GetMaximum()*1.05)); line.SetLineStyle(2); line.SetLineWidth(2); line.Draw(); keep.append(line)
    c.SaveAs(png_out)

    # Small direct TTree trace: stop after max_events useful two-photon events.
    # This is intentionally bounded and therefore cannot reproduce the old multi-minute loop.
    examples=[]
    scanned=0
    for i in range(min(nentries, 20000)):
        ch.GetEntry(i); scanned += 1
        if float(ch.W) <= 2.0 or int(ch.gen_gamma_total) != 2:
            continue
        p0,p1=float(ch.gen_gamma_p[0]),float(ch.gen_gamma_p[1])
        t0,t1=float(ch.gen_gamma_theta[0]),float(ch.gen_gamma_theta[1])
        f0,f1=float(ch.gen_gamma_phi[0]),float(ch.gen_gamma_phi[1])
        m_formula=float(_mgg_from_p_theta_phi(p0,t0,f0,p1,t1,f1))
        # Independent Cartesian construction using degree->radian conversion.
        tr0,fr0,tr1,fr1=np.deg2rad([t0,f0,t1,f1])
        v0=np.array([p0*np.sin(tr0)*np.cos(fr0),p0*np.sin(tr0)*np.sin(fr0),p0*np.cos(tr0)])
        v1=np.array([p1*np.sin(tr1)*np.cos(fr1),p1*np.sin(tr1)*np.sin(fr1),p1*np.cos(tr1)])
        m2_cart=(p0+p1)**2-float(np.dot(v0+v1,v0+v1))
        m_cart=float(np.sqrt(max(0.0,m2_cart)))
        r0,r1=int(ch.gen_gamma_best_rec_index[0]),int(ch.gen_gamma_best_rec_index[1])
        pid0,pid1=int(ch.gen_gamma_best_rec_pid[0]),int(ch.gen_gamma_best_rec_pid[1])
        mrec=None
        if r0>=0 and r1>=0 and r0!=r1 and pid0==22 and pid1==22:
            mrec=float(_mgg_from_p_theta_phi(float(ch.gen_gamma_best_rec_p[0]),float(ch.gen_gamma_best_rec_theta[0]),float(ch.gen_gamma_best_rec_phi[0]),
                                             float(ch.gen_gamma_best_rec_p[1]),float(ch.gen_gamma_best_rec_theta[1]),float(ch.gen_gamma_best_rec_phi[1])))
        examples.append(dict(W=float(ch.W),nph=int(ch.gen_n_photon),npi0ph=int(ch.gen_n_pi0_photon),npi0=int(ch.gen_n_pi0),
                             haslund=int(ch.has_mc_lund),m=m_formula,mcart=m_cart,mrec=mrec,
                             g0=(int(ch.gen_gamma_index[0]),p0,t0,f0,int(ch.gen_gamma_n_rec_matches[0]),int(ch.gen_gamma_n_rec_pid22[0]),r0,pid0,float(ch.gen_gamma_best_rec_p[0]),float(ch.gen_gamma_best_rec_theta[0]),float(ch.gen_gamma_best_rec_phi[0])),
                             g1=(int(ch.gen_gamma_index[1]),p1,t1,f1,int(ch.gen_gamma_n_rec_matches[1]),int(ch.gen_gamma_n_rec_pid22[1]),r1,pid1,float(ch.gen_gamma_best_rec_p[1]),float(ch.gen_gamma_best_rec_theta[1]),float(ch.gen_gamma_best_rec_phi[1]))))
        if len(examples)>=max_events:
            break

    nw=int(n_w.GetValue()); nl=int(n_lund.GetValue()); npart=int(n_part.GetValue()); ntwo=int(n_two_mcpart.GetValue())
    with open(txt_out,"w") as f:
        f.write("AAOgen generator bookkeeping + angle-unit + truth->REC audit\n")
        f.write("All saved theta/phi branches are DEGREES: process_photon_efficiency.groovy writes thetaDeg()/phiDeg().\n")
        f.write("W > 2 GeV is enforced for aggregate counts and displayed examples.\n\n")
        f.write(f"PhotonEfficiencyEvents entries: {nentries:,}\nW>2 entries: {nw:,}\n")
        f.write(f"has MC::Particle: {npart:,} ({100*npart/nw if nw else 0:.2f}%)\n")
        f.write(f"has MC::Lund: {nl:,} ({100*nl/nw if nw else 0:.2f}%)\n")
        f.write(f"gen_gamma_total==2 (MC::Particle): {ntwo:,}\n")
        f.write(f"gen_n_photon==2 (MC::Lund stable photons): {int(n_two_lund.GetValue()):,}\n")
        f.write(f"gen_n_pi0_photon==2 (MC::Lund direct pi0 daughters): {int(n_two_pi0lund.GetValue()):,}\n")
        f.write("\nIMPORTANT: gen_n_photon/gen_n_pi0_photon are MC::Lund counters; gen_gamma_total/gen_gamma_* are MC::Particle. They must not be assumed equivalent.\n")
        f.write(f"\nDirect examples: scanned at most {scanned:,} entries to obtain {len(examples)} W>2, gen_gamma_total==2 events.\n")
        for j,x in enumerate(examples):
            f.write(f"\n[{j}] W={x['W']:.5f} has_lund={x['haslund']} gen_n_photon={x['nph']} gen_n_pi0_photon={x['npi0ph']} gen_n_pi0={x['npi0']}\n")
            f.write(f"    Mgg_gen(degree formula)={x['m']:.8f} GeV  Mgg_gen(Cartesian)={x['mcart']:.8f} GeV  delta={abs(x['m']-x['mcart']):.3e} GeV\n")
            f.write(f"    Mgg_bestREC_PID22={x['mrec'] if x['mrec'] is not None else 'NA'}\n")
            f.write(f"    gamma0(index,p,theta_deg,phi_deg,nmatch,npid22,best_rec,best_pid,best_p,best_theta_deg,best_phi_deg)={x['g0']}\n")
            f.write(f"    gamma1(index,p,theta_deg,phi_deg,nmatch,npid22,best_rec,best_pid,best_p,best_theta_deg,best_phi_deg)={x['g1']}\n")

    print("\nAAOgen generator/truth diagnostic:")
    print("  CRITICAL UNIT CHECK: producer stores theta/phi in DEGREES, not radians.")
    print(f"  W>2 event-tree entries: {nw:,}; has MC::Particle={npart:,}; has MC::Lund={nl:,}")
    print(f"  gen_gamma_total==2 (MC::Particle): {ntwo:,}")
    print(f"  gen_n_photon==2 (MC::Lund): {int(n_two_lund.GetValue()):,}")
    print(f"  gen_n_pi0_photon==2 (MC::Lund): {int(n_two_pi0lund.GetValue()):,}")
    print(f"  wrote generator audit -> {png_out}")
    print(f"  wrote bounded event trace -> {txt_out}")
    return txt_out

def main():
    args = parse_args()

    if args.threads < 1:
        print("ERROR: --threads must be >= 1", file=sys.stderr)
        return 1

    ROOT.EnableImplicitMT(args.threads)
    print(f"ROOT implicit multithreading: {args.threads} worker(s)")

    clear_output_directory(args.output_dir)

    dfs = {}
    count_handles = {}

    print(
        f"Looking for completed skim-version-{EXPECTED_SKIM_VERSION} files..."
    )

    for sample, label in SAMPLES:
        files = discover(sample, args.period)
        if not files:
            continue

        chain = make_chain(files)
        print(
            f"  {label}: {len(files)} file(s), "
            f"{chain.GetEntries():,} PhotonEfficiency rows before W cut"
        )

        df = make_dataframe(chain)
        dfs[sample] = df
        count_handles[sample] = df.Count()

    if not dfs:
        print("ERROR: no usable input files.", file=sys.stderr)
        return 1

    # Materialize the W>2 counts in coordinated ROOT event loops.
    ROOT.RDF.RunGraphs(list(count_handles.values()))
    print("\nRows entering the first canvas (W > 2 GeV and angle(e',gamma1) > 8 deg):")
    for sample, label in SAMPLES:
        if sample in count_handles:
            print(f"  {label:<14s} {int(count_handles[sample].GetValue()):,}")

    exclusivity_dir = os.path.join(args.output_dir, "exclusivity_selection")
    os.makedirs(exclusivity_dir, exist_ok=True)
    keep, written = draw_canvases(dfs, exclusivity_dir, args.period)

    normalization_dir = os.path.join(args.output_dir, "normalization")
    os.makedirs(normalization_dir, exist_ok=True)
    norm_keep, norm_output = draw_normalization_step1(dfs, normalization_dir, args.period)
    keep.extend(norm_keep)
    lowshape_keep, lowshape_output = draw_low_energy_shapes_v7(dfs, normalization_dir, args.period)
    keep.extend(lowshape_keep)
    high_keep, high_output, high_coeffs = draw_high_energy_normalization_v7(dfs, normalization_dir, args.period)
    keep.extend(high_keep)
    low_keep, low_output, denominator_coeffs = draw_low_energy_normalization_v7(dfs, normalization_dir, args.period, high_coeffs)
    keep.extend(low_keep)
    final_keep, final_output = draw_full_range_denominator(dfs, normalization_dir, args.period, denominator_coeffs)
    keep.extend(final_keep)

    probe_dir = os.path.join(args.output_dir, "probe")
    os.makedirs(probe_dir, exist_ok=True)
    probe_keep, probe_output = draw_probe_mgg(dfs, probe_dir, args.period, denominator_coeffs)
    keep.extend(probe_keep)
    pi0fit_keep, pi0fit_output, pi0fit_result = draw_probe_aaogen_pi0_fit(
        dfs, probe_dir, args.period, denominator_coeffs)
    keep.extend(pi0fit_keep)
    _ = keep  # Keep ROOT objects alive through SaveAs().

    for output_file in written:
        print(f"\nWrote: {output_file}")
    print(f"\nWrote: {norm_output}")
    print(f"\nWrote: {lowshape_output}")
    if high_output:
        print(f"\nWrote: {high_output}")
    if low_output:
        print(f"\nWrote: {low_output}")
    if final_output:
        print(f"\nWrote: {final_output}")
    if probe_output:
        print(f"\nWrote: {probe_output}")
    if pi0fit_output:
        print(f"\nWrote: {pi0fit_output}")
    if pi0fit_result and pi0fit_result.get("window_scan_output"):
        print(f"\nWrote: {pi0fit_result['window_scan_output']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
