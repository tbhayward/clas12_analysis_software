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
    // Saved theta/phi branches are in radians.
    const double th1 = th1_deg, ph1 = ph1_deg;
    const double th2 = th2_deg, ph2 = ph2_deg;
    const double dot = std::sin(th1)*std::sin(th2)*std::cos(ph1-ph2)
                     + std::cos(th1)*std::cos(th2);
    const double c = std::max(-1.0, std::min(1.0, dot));
    return std::acos(c) * 180.0 / M_PI;
}

double pe_mx2_egamma(double ebeam, double ep, double eth, double eph,
                     double gp, double gth, double gph) {
    const double me = 0.00051099895;
    const double mp = 0.9382720813;
    const double Ee = std::sqrt(ep*ep + me*me);
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
        const double dot = std::sin(tag_th)*std::sin(neutral_th[i])*std::cos(tag_ph-neutral_ph[i])
                         + std::cos(tag_th)*std::cos(neutral_th[i]);
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
    const double dot = std::sin(tag_th)*std::sin(neutral_th[best])*std::cos(tag_ph-neutral_ph[best])
                     + std::cos(tag_th)*std::cos(neutral_th[best]);
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
        const double dot=std::sin(tag_th)*std::sin(neutral_th[i])*std::cos(tag_ph-neutral_ph[i])
                        +std::cos(tag_th)*std::cos(neutral_th[i]);
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
    const double dot=std::sin(th1)*std::sin(th2)*std::cos(ph1-ph2)+std::cos(th1)*std::cos(th2);
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
        ("Mx2_egamma1", ";M^{2}_{X}(e'#gamma1) (GeV^{2});Unit-normalized entries", 160, -20.0, 20.0, False, "mx2_egamma1"),
    ]

    stages = [
        ("", lambda df: df),
        ("M^{2}_{X}(e'p') < 0.15 GeV^{2}",
         lambda df: df.Filter("Mx2_ep < 0.15", "Mx2_ep_lt_0p15")),
        ("-0.05 < M^{2}_{X}(e'p'#gamma1) < 0.05 GeV^{2}",
         lambda df: df.Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", "Mx2_epg_window")),
    ]

    stage_dfs = {}
    for sample, _label in SAMPLES:
        if sample not in dfs:
            continue
        stage_dfs[(sample, 0)] = dfs[sample]
        stage_dfs[(sample, 1)] = stages[1][1](stage_dfs[(sample, 0)])
        stage_dfs[(sample, 2)] = stages[2][1](stage_dfs[(sample, 1)])

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
        canvas = ROOT.TCanvas(f"c_{slug}_{unique}", "", 1800, 650)
        canvas.Divide(3, 1, 0.002, 0.002)
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
            .Filter("Mx2_ep < 0.15", "norm_Mx2_ep_lt_0p15")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", "norm_Mx2_epg_window")
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

    for iplot, (_expr, _title, _nbins, _xmin, _xmax, logy) in enumerate(plots):
        pad = canvas.cd(iplot + 1)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.13)
        pad.SetTopMargin(0.08)
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
    return (df
            .Filter("Mx2_ep < 0.15", f"{tag}_Mx2_ep_lt_0p15")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", f"{tag}_Mx2_epg_window"))


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
    for ip,spec in enumerate(specs):
        pad=canvas.cd(ip+1); pad.SetTicks(1,1); pad.SetLeftMargin(0.14); pad.SetRightMargin(0.04); pad.SetBottomMargin(0.13); pad.SetTopMargin(0.08)
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


def draw_probe_mgg(dfs, output_dir, period):
    """First probe study: M(gamma_tag gamma_probe) for every retained PID-22 partner.

    The tag is the photon defining the selected epgammaX denominator row.  Every
    other retained neutral candidate with PID 22 is treated as a possible probe.
    No tag-probe angular matching or pi0-mass requirement is imposed here.
    """
    selected={}
    for sample,df in dfs.items():
        selected[sample]=(_exclusive_df(df,f"probe_{sample}")
            .Define("Mgg_tag_probe",
                "pe_mgg_tag_probe(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,"
                "neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi)"))
    unique=str(abs(hash((period,"probe_mgg"))))
    booked={}; actions=[]
    for sample,df in selected.items():
        h=df.Histo1D((f"h_probe_mgg_{sample}_{unique}",
                      ";M_{#gamma_{tag}#gamma_{probe}} (GeV);Unit-normalized tag-probe combinations",
                      160,0.0,0.8),"Mgg_tag_probe")
        booked[sample]=h; actions.append(h)
    if actions: ROOT.RDF.RunGraphs(actions)
    canvas=ROOT.TCanvas(f"c_probe_mgg_{unique}","",1000,800); canvas.SetTicks(1,1); canvas.SetLeftMargin(0.13); canvas.SetRightMargin(0.04); canvas.SetBottomMargin(0.13); canvas.SetTopMargin(0.08)
    keep=[canvas]+actions; hs=[]; ymax=0.0
    leg=ROOT.TLegend(0.50,0.66,0.89,0.88); leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.029); keep.append(leg)
    for sample,label in SAMPLES:
        if sample not in booked: continue
        h=booked[sample].GetValue().Clone(f"hp_probe_{sample}_{unique}"); h.SetDirectory(0); h.SetStats(0); h.SetLineColor(COLORS[sample]); h.SetLineWidth(3)
        n=int(round(h.GetEntries())); integ=h.Integral(1,h.GetNbinsX())
        if integ>0: h.Scale(1.0/integ)
        ymax=max(ymax,h.GetMaximum()); hs.append(h); keep.append(h); leg.AddEntry(h,f"{label} (pairs={n:,})","l")
    for i,h in enumerate(hs):
        h.SetMinimum(0.0); h.SetMaximum(1.22*ymax if ymax else 1.0)
        h.GetXaxis().SetTitleSize(0.047); h.GetYaxis().SetTitleSize(0.043); h.GetXaxis().SetLabelSize(0.038); h.GetYaxis().SetLabelSize(0.038); h.GetYaxis().SetTitleOffset(1.35)
        h.Draw("HIST" if i==0 else "HIST SAME")
    leg.Draw()
    out=os.path.join(output_dir,f"1_{period}_Mgg_tag_probe.png"); canvas.SaveAs(out)
    return keep,out


def draw_probe_mgg_truth_diagnostic(dfs, output_dir, period):
    """AAOgen-only four-vector/truth audit for the missing pi0 peak.

    The four panels deliberately answer one question at a time:
      1) does the inclusive reconstructed pairing reproduce the surprising smooth shape?
      2) does choosing only the nearest retained PID-22 candidate recover the peak?
      3) when the saved MC roles say tag/probe are pi0 daughters and the reconstructed
         neutral is matched to the saved generated probe, does reconstructed Mgg peak?
      4) do the saved generated tag/probe four-vectors themselves give m_pi0?
    """
    if "aaogen" not in dfs:
        return [], None
    df=_exclusive_df(dfs["aaogen"],"probe_truth_audit_aaogen")
    df=(df
        .Define("Mgg_audit_all","pe_mgg_tag_probe(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi)")
        .Define("Mgg_audit_nearest","pe_mgg_nearest_pid22(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi,neutral_delta_alpha)")
        .Define("Mgg_audit_truth_reco","pe_mgg_truth_matched_pi0_reco(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,mc_tag_index,mc_tag_pid,mc_tag_parent,mc_probe_index,mc_probe_pid,mc_probe_parent,neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi,neutral_mc_index,neutral_mc_pid)")
        .Define("Mgg_audit_truth_gen","pe_mgg_truth_pi0(mc_tag_p,mc_tag_theta,mc_tag_phi,mc_tag_index,mc_tag_pid,mc_tag_parent,mc_probe_p,mc_probe_theta,mc_probe_phi,mc_probe_index,mc_probe_pid,mc_probe_parent)"))
    unique=str(abs(hash((period,"probe_truth_audit"))))
    specs=[
      ("Mgg_audit_all","All retained reconstructed PID-22 partners"),
      ("Mgg_audit_nearest","Nearest-to-X retained PID-22 partner only"),
      ("Mgg_audit_truth_reco","Reco pair matched to saved #pi^{0} truth tag/probe"),
      ("Mgg_audit_truth_gen","Saved generated #pi^{0} tag/probe four-vectors")]
    booked=[]; actions=[]
    for i,(col,title) in enumerate(specs):
        d=df if i<3 else df.Filter("Mgg_audit_truth_gen >= 0.0","valid generated pi0 pair")
        h=d.Histo1D((f"h_probe_truth_audit_{i}_{unique}",f";M_{{#gamma#gamma}} (GeV);Unit-normalized entries",120,0.0,0.30),col)
        booked.append((h,title)); actions.append(h)
    ROOT.RDF.RunGraphs(actions)
    c=ROOT.TCanvas(f"c_probe_truth_audit_{unique}","",1400,1050); c.Divide(2,2)
    keep=[c]+actions
    for i,(hp,title) in enumerate(booked,1):
        pad=c.cd(i); pad.SetTicks(1,1); pad.SetLeftMargin(0.13); pad.SetRightMargin(0.04); pad.SetBottomMargin(0.13); pad.SetTopMargin(0.12)
        h=hp.GetValue().Clone(f"h_probe_truth_audit_draw_{i}_{unique}"); h.SetDirectory(0); h.SetStats(0); h.SetLineColor(ROOT.kRed+1); h.SetLineWidth(3)
        n=int(round(h.GetEntries())); integ=h.Integral(1,h.GetNbinsX())
        if integ>0: h.Scale(1.0/integ)
        h.SetMinimum(0.0); h.SetMaximum(1.22*h.GetMaximum() if h.GetMaximum()>0 else 1.0)
        h.GetXaxis().SetTitleSize(0.047); h.GetYaxis().SetTitleSize(0.043); h.GetXaxis().SetLabelSize(0.038); h.GetYaxis().SetLabelSize(0.038); h.GetYaxis().SetTitleOffset(1.35)
        h.Draw("HIST")
        line=ROOT.TLine(0.1349768,0.0,0.1349768,h.GetMaximum()); line.SetLineStyle(2); line.SetLineWidth(2); line.Draw()
        tex=ROOT.TLatex(); tex.SetNDC(True); tex.SetTextSize(0.034); tex.SetTextFont(42); tex.DrawLatex(0.14,0.94,title); tex.DrawLatex(0.62,0.86,f"N = {n:,}")
        keep += [h,line,tex]
    out=os.path.join(output_dir,f"2_{period}_AAOgen_Mgg_truth_audit.png"); c.SaveAs(out)
    print("\nAAOgen M(gamma gamma) truth audit:")
    print("  Panel 1: all reconstructed PID-22 partners -- reproduces the original construction.")
    print("  Panel 2: nearest retained PID-22 to inferred X -- tests combinatorial dilution.")
    print("  Panel 3: reconstructed pair whose saved MC identities are the pi0 tag/probe photons -- tests REC four-vectors/matching.")
    print("  Panel 4: same saved generated pi0 tag/probe four-vectors -- must peak at ~0.135 GeV if truth roles are sane.")
    return keep,out

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
    probe_keep, probe_output = draw_probe_mgg(dfs, probe_dir, args.period)
    keep.extend(probe_keep)
    probe_audit_keep, probe_audit_output = draw_probe_mgg_truth_diagnostic(dfs, probe_dir, args.period)
    keep.extend(probe_audit_keep)
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
    if probe_audit_output:
        print(f"\nWrote: {probe_audit_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
