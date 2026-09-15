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

ROOT::VecOps::RVec<double> pe_delta_p_tag_probe(
        double probe_missing_p, int tag_index,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<double>& neutral_p) {
    // Signed momentum residual for each reconstructed candidate probe photon:
    //
    //   Delta p = |p_X| - |p_gamma_candidate|
    //
    // probe_missing_p is the magnitude of the missing three-momentum after
    // subtracting e', p', and the tag photon gamma1.  In an exclusive
    // e'p'gamma1 gamma2 event, X is gamma2, so Delta p should peak near zero
    // for the correct reconstructed partner.
    //
    // Use exactly the same reconstructed partner population as pe_mgg_tag_probe:
    // retained PID-22 neutrals, excluding the tag itself, with positive momentum.
    ROOT::VecOps::RVec<double> out;
    const auto n = neutral_idx.size();
    out.reserve(n);
    if (!(probe_missing_p > 0.0)) return out;
    for (size_t i = 0; i < n; ++i) {
        if (neutral_idx[i] < 0) continue;
        if (neutral_idx[i] == tag_index) continue;
        if (neutral_pid[i] != 22) continue;
        if (!(neutral_p[i] > 0.0)) continue;
        out.push_back(probe_missing_p - neutral_p[i]);
    }
    return out;
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


ROOT::VecOps::RVec<double> pe_mgg_truth_by_matched_mc(
        double tag_p, double tag_th, double tag_ph, int tag_index,
        int mc_tag_index, int mc_tag_pid,
        double mc_tag_p, double mc_tag_th, double mc_tag_ph,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<double>& neutral_p,
        const ROOT::VecOps::RVec<double>& neutral_th,
        const ROOT::VecOps::RVec<double>& neutral_ph,
        const ROOT::VecOps::RVec<int>& neutral_mc_index,
        const ROOT::VecOps::RVec<int>& neutral_mc_pid,
        const ROOT::VecOps::RVec<double>& neutral_mc_p,
        const ROOT::VecOps::RVec<double>& neutral_mc_th,
        const ROOT::VecOps::RVec<double>& neutral_mc_ph) {
    ROOT::VecOps::RVec<double> out;
    if (mc_tag_index < 0 || mc_tag_pid != 22 || !(mc_tag_p > 0.0)) return out;
    for (size_t i=0;i<neutral_idx.size();++i) {
        if (neutral_idx[i] < 0 || neutral_idx[i] == tag_index || neutral_pid[i] != 22) continue;
        if (i >= neutral_mc_index.size() || i >= neutral_mc_pid.size() ||
            i >= neutral_mc_p.size() || i >= neutral_mc_th.size() || i >= neutral_mc_ph.size()) continue;
        if (neutral_mc_index[i] < 0 || neutral_mc_pid[i] != 22) continue;
        if (neutral_mc_index[i] == mc_tag_index || !(neutral_mc_p[i] > 0.0)) continue;

        // Truth test uses the two generated photons actually associated with this
        // reconstructed pair.  No parent-PDG or saved "probe role" assumption.
        const double tdot = std::sin(mc_tag_th*M_PI/180.0)*std::sin(neutral_mc_th[i]*M_PI/180.0)*
                            std::cos((mc_tag_ph-neutral_mc_ph[i])*M_PI/180.0) +
                            std::cos(mc_tag_th*M_PI/180.0)*std::cos(neutral_mc_th[i]*M_PI/180.0);
        const double tc = std::max(-1.0,std::min(1.0,tdot));
        const double tm2 = 2.0*mc_tag_p*neutral_mc_p[i]*(1.0-tc);
        if (!(tm2 >= 0.0)) continue;
        const double tm = std::sqrt(tm2);
        if (!(tm > 0.125 && tm < 0.145)) continue;

        const double rdot = std::sin(tag_th*M_PI/180.0)*std::sin(neutral_th[i]*M_PI/180.0)*
                            std::cos((tag_ph-neutral_ph[i])*M_PI/180.0) +
                            std::cos(tag_th*M_PI/180.0)*std::cos(neutral_th[i]*M_PI/180.0);
        const double rc = std::max(-1.0,std::min(1.0,rdot));
        const double rm2 = 2.0*tag_p*neutral_p[i]*(1.0-rc);
        if (rm2 >= 0.0) out.push_back(std::sqrt(rm2));
    }
    return out;
}

ROOT::VecOps::RVec<int> pe_mgg_truth_match_stage(
        int tag_index, int mc_tag_index, int mc_tag_pid, double mc_tag_p,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<int>& neutral_mc_index,
        const ROOT::VecOps::RVec<int>& neutral_mc_pid,
        const ROOT::VecOps::RVec<double>& neutral_mc_p) {
    ROOT::VecOps::RVec<int> out;
    for (size_t i=0;i<neutral_idx.size();++i) {
        if (neutral_idx[i] < 0 || neutral_idx[i] == tag_index || neutral_pid[i] != 22) continue;
        int stage=1; // reconstructed tag + reconstructed PID22 partner exists
        if (mc_tag_index >= 0) stage=2; else { out.push_back(stage); continue; }
        if (mc_tag_pid == 22) stage=3; else { out.push_back(stage); continue; }
        if (mc_tag_p > 0.0) stage=4; else { out.push_back(stage); continue; }
        if (i < neutral_mc_index.size() && neutral_mc_index[i] >= 0) stage=5; else { out.push_back(stage); continue; }
        if (i < neutral_mc_pid.size() && neutral_mc_pid[i] == 22) stage=6; else { out.push_back(stage); continue; }
        if (i < neutral_mc_p.size() && neutral_mc_p[i] > 0.0) stage=7; else { out.push_back(stage); continue; }
        if (neutral_mc_index[i] != mc_tag_index) stage=8;
        out.push_back(stage);
    }
    return out;
}

ROOT::VecOps::RVec<int> pe_mgg_truth_role_diagnostics(
        int tag_index, int mc_tag_index, int mc_tag_pid, int mc_tag_parent,
        int mc_probe_index, int mc_probe_pid, int mc_probe_parent,
        const ROOT::VecOps::RVec<int>& neutral_idx,
        const ROOT::VecOps::RVec<int>& neutral_pid,
        const ROOT::VecOps::RVec<int>& neutral_mc_index,
        const ROOT::VecOps::RVec<int>& neutral_mc_pid) {
    // One bit-mask per reconstructed tag-partner combination.  This diagnoses
    // the *old* assumptions; it is not used to define truth in the new closure.
    ROOT::VecOps::RVec<int> out;
    for (size_t i=0;i<neutral_idx.size();++i) {
        if (neutral_idx[i] < 0 || neutral_idx[i] == tag_index || neutral_pid[i] != 22) continue;
        int mask=0;
        if (mc_tag_index >= 0) mask |= 1;
        if (mc_tag_pid == 22) mask |= 2;
        if (mc_tag_parent == 111) mask |= 4;
        if (mc_probe_index >= 0) mask |= 8;
        if (mc_probe_pid == 22) mask |= 16;
        if (mc_probe_parent == 111) mask |= 32;
        if (i < neutral_mc_index.size() && neutral_mc_index[i] == mc_probe_index) mask |= 64;
        if (i < neutral_mc_pid.size() && neutral_mc_pid[i] == 22) mask |= 128;
        out.push_back(mask);
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
    bounds = ((1.0e-12, None), (1.0e-12, None), (-5.0, 5.0), (0.0, 4.0))
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



def draw_probe_delta_p(dfs, output_dir, period, coeffs):
    """Plot Delta p = |p_X| - |p_gamma_candidate| for candidate probe photons.

    p_X is the missing three-momentum after e', p', and the tag photon gamma1:
        p_X = p_beam + p_target - p_e' - p_p' - p_gamma1.
    The target has zero three-momentum, so the saved probe_raw_p branch is
    precisely |p_X| from that missing four-vector.

    Candidate gamma_probe objects are exactly the same retained reconstructed
    PID-22 partners used in the Mgg probe plot.  The same final epgammaX
    selection and the same A/B/C normalization coefficients are used.

    Two views are written on one canvas:
      left:  broad -4 < Delta p < 4 GeV
      right: zoom  -1 < Delta p < 1 GeV
    No Delta-p cut is applied.
    """
    if coeffs is None:
        print("WARNING: no denominator normalization coefficients; skipping probe Delta-p plot.")
        return [], None

    selected = {}
    for sample, df in dfs.items():
        selected[sample] = (_exclusive_df(df, f"probe_deltap_{sample}")
            .Define("delta_p_tag_probe",
                "pe_delta_p_tag_probe(probe_raw_p,tag_rec_index,"
                "neutral_idx,neutral_pid,neutral_p)"))

    A, B, C = coeffs["mean"]
    scales = {"data": 1.0, "dvcsgen": A, "aaogen": B, "clasdis": C}
    unique = str(abs(hash((period, "probe_delta_p_scaled_components"))))

    # Book a single broad histogram; the second pad is only an axis zoom of the
    # identical histogram, so both panels contain exactly the same combinations.
    booked, actions = {}, []
    for sample, df in selected.items():
        h = df.Histo1D(
            (f"h_probe_deltap_{sample}_{unique}",
             ";#Delta p = |p_{X}| - |p_{#gamma_{probe}}| (GeV);Normalized tag-probe combinations",
             320, -4.0, 4.0),
            "delta_p_tag_probe")
        booked[sample] = h
        actions.append(h)

    if actions:
        ROOT.RDF.RunGraphs(actions)

    drawn = {}
    raw_pairs = {}
    keep = list(actions)
    for sample, _label in SAMPLES:
        if sample not in booked:
            continue
        h = booked[sample].GetValue().Clone(f"hp_deltap_{sample}_{unique}")
        h.SetDirectory(0)
        h.SetStats(0)
        raw_pairs[sample] = int(round(h.GetEntries()))
        if sample == "data":
            h.SetLineColor(ROOT.kBlack)
            h.SetLineWidth(3)
            h.SetLineStyle(1)
        else:
            h.Scale(scales[sample])
            h.SetLineColor(COLORS[sample])
            h.SetLineWidth(2)
            h.SetLineStyle(2)
        drawn[sample] = h
        keep.append(h)

    combined = None
    for sample in ("dvcsgen", "aaogen", "clasdis"):
        if sample not in drawn:
            continue
        if combined is None:
            combined = drawn[sample].Clone(f"h_probe_deltap_combined_{unique}")
            combined.SetDirectory(0)
        else:
            combined.Add(drawn[sample])
    if combined is not None:
        combined.SetStats(0)
        combined.SetLineColor(ROOT.kMagenta + 2)
        combined.SetLineWidth(4)
        combined.SetLineStyle(1)
        keep.append(combined)

    canvas = ROOT.TCanvas(f"c_probe_deltap_{unique}", "", 1500, 720)
    canvas.Divide(2, 1, 0.002, 0.002)
    keep.append(canvas)

    def draw_panel(pad, xmin, xmax, title):
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.11)

        first = drawn.get("data", combined)
        if first is None and drawn:
            first = next(iter(drawn.values()))
        if first is None:
            return

        all_h = list(drawn.values()) + ([combined] if combined is not None else [])
        ymax = 0.0
        for h in all_h:
            b1 = h.GetXaxis().FindBin(xmin + 1e-9)
            b2 = h.GetXaxis().FindBin(xmax - 1e-9)
            for ib in range(b1, b2 + 1):
                ymax = max(ymax, h.GetBinContent(ib))

        first.GetXaxis().SetRangeUser(xmin, xmax)
        first.SetMinimum(0.0)
        first.SetMaximum(1.20 * ymax if ymax > 0 else 1.0)
        first.GetXaxis().SetTitleSize(0.046)
        first.GetYaxis().SetTitleSize(0.041)
        first.GetXaxis().SetLabelSize(0.036)
        first.GetYaxis().SetLabelSize(0.036)
        first.GetYaxis().SetTitleOffset(1.55)
        first.Draw("HIST")

        if combined is not None and combined is not first:
            combined.GetXaxis().SetRangeUser(xmin, xmax)
            combined.Draw("HIST SAME")
        for sample in ("dvcsgen", "aaogen", "clasdis"):
            if sample in drawn and drawn[sample] is not first:
                drawn[sample].GetXaxis().SetRangeUser(xmin, xmax)
                drawn[sample].Draw("HIST SAME")
        if "data" in drawn:
            drawn["data"].GetXaxis().SetRangeUser(xmin, xmax)
            drawn["data"].Draw("HIST SAME")

        zero = ROOT.TLine(0.0, 0.0, 0.0, first.GetMaximum())
        zero.SetLineColor(ROOT.kGray + 2)
        zero.SetLineStyle(3)
        zero.SetLineWidth(2)
        zero.Draw()
        keep.append(zero)

        leg = ROOT.TLegend(0.53, 0.62, 0.94, 0.87)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.027)
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

        lab = ROOT.TLatex()
        lab.SetNDC(True)
        lab.SetTextAlign(22)
        lab.SetTextSize(0.036)
        lab.DrawLatex(0.53, 0.955, title)
        keep.extend([leg, lab])

    draw_panel(canvas.cd(1), -5.0, 5.0, "Probe momentum residual: broad view")
    draw_panel(canvas.cd(2), -1.0, 1.0, "Probe momentum residual: near-exclusive region")

    out = os.path.join(output_dir, f"5_{period}_delta_p_tag_probe.png")
    canvas.SaveAs(out)

    print("\nProbe Delta-p comparison:")
    print("  Definition: Delta p = |p_X(e'p'gamma1)| - |p_gamma_probe|")
    print("  p_X uses the saved probe_raw_p missing-momentum magnitude.")
    print("  Candidate probes are the same retained PID-22 partners used in Mgg.")
    print("  Same _exclusive_df() selection as probe Mgg; NO Delta-p cut applied.")
    print(f"  MC scales: DVCSgen={A:.8g}, AAOgen={B:.8g}, CLASDIS={C:.8g}")
    print(f"  Wrote broad (-5,5) and zoomed (-1,1) views to {out}")
    return keep, out

def draw_probe_integrated_delta_p_efficiency(dfs, output_dir, period, coeffs):
    """Integrated Delta-p matching study; no p/theta/phi binning."""
    if coeffs is None or "data" not in dfs or "aaogen" not in dfs:
        print("WARNING: data/AAOgen unavailable; skipping integrated Delta-p efficiency.")
        return [], None, None
    A, B, C = coeffs["mean"]
    unique = str(abs(hash((period, "integrated_delta_p_efficiency"))))
    selected = {}
    for sample in ("data", "aaogen"):
        selected[sample] = (_exclusive_df(dfs[sample], f"intdp_{sample}")
            .Define("intdp_vec", "pe_delta_p_tag_probe(probe_raw_p,tag_rec_index,neutral_idx,neutral_pid,neutral_p)"))
    hs, acts = {}, []
    for sample in ("data", "aaogen"):
        h = selected[sample].Histo1D((f"h_intdp_{sample}_{unique}", "", 200, -1.0, 1.0), "intdp_vec")
        hs[sample] = h; acts.append(h)
    ROOT.RDF.RunGraphs(acts)
    hdata = hs["data"].GetValue().Clone(f"h_intdp_data_draw_{unique}"); hdata.SetDirectory(0)
    hmc = hs["aaogen"].GetValue().Clone(f"h_intdp_aaogen_draw_{unique}"); hmc.SetDirectory(0); hmc.Scale(B)

    def fit_peak(h, name):
        ax=h.GetXaxis(); b1=ax.FindBin(-0.45); b2=ax.FindBin(0.45)
        ib=max(range(b1,b2+1), key=lambda b:h.GetBinContent(b)); mode=ax.GetBinCenter(ib)
        # Exclusive peak + smooth combinatorial background.  The old Gaussian-only
        # fit let the broad tails inflate sigma badly.  Fit data and AAOgen independently.
        lo,hi=-1.00,1.00
        f=ROOT.TF1(name,"gaus(0)+pol2(3)",lo,hi)
        amp=max(h.GetBinContent(ib),1.0)
        edge=0.5*(h.GetBinContent(ax.FindBin(lo+0.03))+h.GetBinContent(ax.FindBin(hi-0.03)))
        f.SetParameters(amp,mode,0.12,max(edge,0.0),0.0,0.0)
        f.SetParLimits(0,0.0,max(10.0*amp,1.0e9)); f.SetParLimits(1,-0.35,0.45); f.SetParLimits(2,0.015,0.35)
        h.Fit(f,"QNR")
        mu = float(f.GetParameter(1))
        sigma = abs(float(f.GetParameter(2)))
        return f, mu, sigma

    fdata, mu_data, sig_data = fit_peak(hdata, f"f_intdp_data_{unique}")
    fmc, mu_mc, sig_mc = fit_peak(hmc, f"f_intdp_mc_{unique}")

    # One best candidate per already-selected epgammaX row.  This preserves the
    # denominator population while preventing multiple probe photons on that row
    # from producing multiple numerator counts.
    if not hasattr(draw_probe_integrated_delta_p_efficiency, "_helper_declared"):
        ROOT.gInterpreter.Declare(r'''
        double pe_best_delta_p(const ROOT::VecOps::RVec<double>& v, double mu) {
            if (v.empty()) return 1.0e9;
            double best=v[0], d=std::abs(v[0]-mu);
            for (size_t i=1;i<v.size();++i) {
                const double di=std::abs(v[i]-mu);
                if (di<d) { d=di; best=v[i]; }
            }
            return best;
        }

        ROOT::VecOps::RVec<double> pe_best_delta_p_truth(
            double probe_p,double mu,int tag_index,
            int mc_tag_index,int mc_tag_pid,double mc_tag_p,double mc_tag_th,double mc_tag_ph,
            const ROOT::VecOps::RVec<int>& neutral_idx,
            const ROOT::VecOps::RVec<int>& neutral_pid,
            const ROOT::VecOps::RVec<double>& neutral_p,
            const ROOT::VecOps::RVec<int>& neutral_mc_index,
            const ROOT::VecOps::RVec<int>& neutral_mc_pid,
            const ROOT::VecOps::RVec<double>& neutral_mc_p,
            const ROOT::VecOps::RVec<double>& neutral_mc_th,
            const ROOT::VecOps::RVec<double>& neutral_mc_ph) {
            ROOT::VecOps::RVec<double> out;
            int best=-1; double bestdp=1e9,bestdist=1e99;
            for(size_t i=0;i<neutral_idx.size();++i){
                if(neutral_idx[i]<0 || neutral_idx[i]==tag_index || neutral_pid[i]!=22 || !(neutral_p[i]>0)) continue;
                const double dp=probe_p-neutral_p[i], d=std::abs(dp-mu);
                if(d<bestdist){bestdist=d; bestdp=dp; best=(int)i;}
            }
            if(best<0) return out;
            bool truth=false; const size_t i=(size_t)best;
            if(mc_tag_index>=0 && mc_tag_pid==22 && mc_tag_p>0 &&
               i<neutral_mc_index.size() && i<neutral_mc_pid.size() &&
               i<neutral_mc_p.size() && i<neutral_mc_th.size() && i<neutral_mc_ph.size() &&
               neutral_mc_index[i]>=0 && neutral_mc_pid[i]==22 &&
               neutral_mc_index[i]!=mc_tag_index && neutral_mc_p[i]>0){
                const double dot=std::sin(mc_tag_th*M_PI/180.)*std::sin(neutral_mc_th[i]*M_PI/180.)*
                    std::cos((mc_tag_ph-neutral_mc_ph[i])*M_PI/180.)+
                    std::cos(mc_tag_th*M_PI/180.)*std::cos(neutral_mc_th[i]*M_PI/180.);
                const double c=std::max(-1.,std::min(1.,dot));
                const double m2=2.*mc_tag_p*neutral_mc_p[i]*(1.-c);
                if(m2>=0){const double m=std::sqrt(m2); truth=(m>0.125 && m<0.145);}
            }
            out.push_back(bestdp); out.push_back(truth?1.0:0.0); return out;
        }
        ''')
        draw_probe_integrated_delta_p_efficiency._helper_declared = True
    # Nominal candidate ranking is purely reconstructed and fit-independent:
    # choose the photon with the smallest |Delta p|.  This avoids feeding the
    # unstable Gaussian+polynomial peak position back into candidate selection.
    candidate_mu = 0.0
    ddata=selected["data"].Define("best_dp",f"pe_best_delta_p(intdp_vec,{candidate_mu:.17g})")
    dmc=selected["aaogen"].Define("best_dp",f"pe_best_delta_p(intdp_vec,{candidate_mu:.17g})")
    dmc=dmc.Define("best_dp_truth_info",
        f"pe_best_delta_p_truth(probe_raw_p,{candidate_mu:.17g},tag_rec_index,mc_tag_index,mc_tag_pid,"
        "mc_tag_p,mc_tag_theta,mc_tag_phi,neutral_idx,neutral_pid,neutral_p,"
        "neutral_mc_index,neutral_mc_pid,neutral_mc_p,neutral_mc_theta,neutral_mc_phi)")
    dmc=dmc.Define("best_dp_truth","best_dp_truth_info.size()>1 ? best_dp_truth_info[1] : -1.0")
    den_data,den_mc=ddata.Count(),dmc.Count(); cand_data=ddata.Filter("best_dp < 1e8").Count(); cand_mc=dmc.Filter("best_dp < 1e8").Count()
    hbdp=ddata.Filter("best_dp < 1e8").Histo1D((f"h_bestdp_data_{unique}","",200,-1.0,1.0),"best_dp")
    hbmc=dmc.Filter("best_dp < 1e8").Histo1D((f"h_bestdp_mc_{unique}","",200,-1.0,1.0),"best_dp")
    htruth=dmc.Filter("best_dp < 1e8 && best_dp_truth > 0.5").Histo1D((f"h_bestdp_truth_{unique}","",200,-1.0,1.0),"best_dp")
    hwrong=dmc.Filter("best_dp < 1e8 && best_dp_truth < 0.5").Histo1D((f"h_bestdp_wrong_{unique}","",200,-1.0,1.0),"best_dp")
    ntruth=dmc.Filter("best_dp < 1e8 && best_dp_truth > 0.5").Count()
    nwrong=dmc.Filter("best_dp < 1e8 && best_dp_truth < 0.5").Count()
    ROOT.RDF.RunGraphs([den_data,den_mc,cand_data,cand_mc,hbdp,hbmc,htruth,hwrong,ntruth,nwrong])
    nd=int(den_data.GetValue()); nm=int(den_mc.GetValue())
    # Refit the actual one-best-candidate-per-row spectra.  These are the fits used
    # for the final mu, sigma, background subtraction, and efficiency.
    hdata=hbdp.GetValue().Clone(f"h_bestdp_data_draw_{unique}"); hdata.SetDirectory(0)
    hmc=hbmc.GetValue().Clone(f"h_bestdp_mc_draw_{unique}"); hmc.SetDirectory(0); hmc.Scale(B)
    fdata,mu_data,sig_data=fit_peak(hdata,f"f_intdp_data_best_{unique}")
    fmc,mu_mc,sig_mc=fit_peak(hmc,f"f_intdp_mc_best_{unique}")
    for label, mu_value, sigma_value in (
        ("data", mu_data, sig_data),
        ("AAOgen", mu_mc, sig_mc),
    ):
        # Native-float finite checks without relying on the math module.
        # NaN fails x == x; +/-inf are rejected by the explicit magnitude bound.
        if (mu_value != mu_value or sigma_value != sigma_value
                or abs(mu_value) >= 1.0e100 or abs(sigma_value) >= 1.0e100
                or sigma_value <= 0.0):
            raise RuntimeError(
                f"Integrated Delta-p {label} fit returned invalid parameters: "
                f"mu={mu_value}, sigma={sigma_value}"
            )

    def components(f,prefix):
        sg=ROOT.TF1(prefix+"_sig","gaus",-1.00,1.00); sg.SetParameters(f.GetParameter(0),f.GetParameter(1),f.GetParameter(2))
        bg=ROOT.TF1(prefix+"_bg","pol2",-1.00,1.00); bg.SetParameters(f.GetParameter(3),f.GetParameter(4),f.GetParameter(5))
        return sg,bg
    sigfunc_data,bg_data=components(fdata,f"intdp_data_{unique}")
    sigfunc_mc,bg_mc=components(fmc,f"intdp_mc_{unique}")

    def bgsub_window(h,bg,mu,sig,n,scale):
        # Convert all scalar inputs to native Python numbers before arithmetic.
        # This avoids cppyy/PyROOT operator dispatch (NotImplementedError).
        mu = float(mu)
        sig = abs(float(sig))
        n = float(n)
        scale = float(scale)
        lo = max(-1.0, mu - n*sig)
        hi = min(1.0, mu + n*sig)
        ax = h.GetXaxis()
        b1=ax.FindBin(lo+1e-9); b2=ax.FindBin(hi-1e-9); obs=float(h.Integral(b1,b2))/scale
        bkg=float(bg.Integral(lo,hi)/h.GetBinWidth(1))/scale; signal=max(0.0,obs-bkg)
        return obs,bkg,signal
    results={}
    for n in (1,2,3):
        od,bd,sd=bgsub_window(hdata,bg_data,mu_data,sig_data,n,1.0)
        om,bm,sm=bgsub_window(hmc,bg_mc,mu_mc,sig_mc,n,B)
        ed=sd/nd if nd else float('nan'); em=sm/nm if nm else float('nan')
        results[n]={"Ndata":sd,"Nmc":sm,"obs_data":od,"bkg_data":bd,"obs_mc":om,"bkg_mc":bm,"eff_data":ed,"eff_mc":em,"ratio":ed/em if em>0 else float('nan')}

    c=ROOT.TCanvas(f"c_intdp_eff_{unique}","",1500,720)
    c.Divide(2,1,0.002,0.002)
    keep=[c,hdata,hmc,fdata,fmc,sigfunc_data,bg_data,sigfunc_mc,bg_mc]
    panels = [
        (hdata,fdata,sigfunc_data,bg_data,mu_data,sig_data,
         "Data: integrated #Delta p fit",ROOT.kBlack),
        (hmc,fmc,sigfunc_mc,bg_mc,mu_mc,sig_mc,
         "AAOgen: integrated #Delta p fit",ROOT.kRed+1),
    ]
    for ipad,(h,f,sg,bg,mu,sig,title,col) in enumerate(panels,1):
        pad=c.cd(ipad); pad.SetTicks(1,1); pad.SetLeftMargin(0.14); pad.SetBottomMargin(0.14); pad.SetTopMargin(0.11); pad.SetRightMargin(0.04)
        h.SetStats(0); h.SetLineColor(col); h.SetMarkerColor(col); h.SetMarkerStyle(20); h.SetMarkerSize(0.55)
        h.GetXaxis().SetTitle("#Delta p = |p_{X}| - |p_{#gamma_{probe}}| (GeV)"); h.GetYaxis().SetTitle("Normalized candidate combinations"); h.GetYaxis().SetTitleOffset(1.55); h.Draw("E1")
        f.SetLineColor(ROOT.kMagenta+2); f.SetLineWidth(3); f.Draw("SAME")
        sg.SetLineColor(ROOT.kBlue); sg.SetLineWidth(3); sg.Draw("SAME")
        bg.SetLineColor(ROOT.kGray+2); bg.SetLineStyle(2); bg.SetLineWidth(2); bg.Draw("SAME")
        lab=ROOT.TLatex(); lab.SetNDC(True); lab.SetTextSize(0.031)
        lab.DrawLatex(0.18,0.92,title)
        lab.DrawLatex(0.18,0.865,f"#mu = {mu:+.4f} GeV")
        lab.DrawLatex(0.18,0.815,f"#sigma = {sig:.4f} GeV")
        chi2ndf = float(f.GetChisquare())/float(f.GetNDF()) if f.GetNDF() > 0 else float("nan")
        lab.DrawLatex(0.18,0.765,f"#chi^{{2}}/ndf = {chi2ndf:.2f}")
        y=0.705
        for n in (1,2,3):
            r=results[n]
            txt=(f"{n}#sigma: #epsilon_{{data}}={r['eff_data']:.4f}, C_{{#gamma}}={r['ratio']:.4f}" if ipad==1 else f"{n}#sigma: #epsilon_{{AAO}}={r['eff_mc']:.4f}")
            lab.DrawLatex(0.18,y,txt); y-=0.05
        keep.append(lab)
    out=os.path.join(output_dir,f"6_{period}_integrated_delta_p_efficiency.png"); c.SaveAs(out)

    # Direct AAOgen truth decomposition of exactly the same best-candidate population.
    ht=htruth.GetValue().Clone(f"h_bestdp_truth_draw_{unique}"); ht.SetDirectory(0); ht.Scale(B)
    hw=hwrong.GetValue().Clone(f"h_bestdp_wrong_draw_{unique}"); hw.SetDirectory(0); hw.Scale(B)
    nt,nw=int(ntruth.GetValue()),int(nwrong.GetValue())
    ct=ROOT.TCanvas(f"c_bestdp_truth_{unique}","",1450,760)
    ct.SetTicks(1,1); ct.SetLeftMargin(0.12); ct.SetRightMargin(0.04); ct.SetBottomMargin(0.14); ct.SetTopMargin(0.11)
    hmc.SetStats(0); hmc.SetLineColor(ROOT.kBlack); hmc.SetLineWidth(2); hmc.GetXaxis().SetRangeUser(-1.0,1.0)
    hmc.GetXaxis().SetTitle("#Delta p = |p_{X}| - |p_{#gamma_{probe}}| (GeV)")
    hmc.GetYaxis().SetTitle("Normalized best-candidate rows"); hmc.GetYaxis().SetTitleOffset(1.35); hmc.Draw("HIST")
    ht.SetLineColor(ROOT.kBlue+1); ht.SetLineWidth(3); ht.Draw("HIST SAME")
    hw.SetLineColor(ROOT.kRed+1); hw.SetLineWidth(3); hw.SetLineStyle(2); hw.Draw("HIST SAME")
    lt=ROOT.TLegend(0.56,0.69,0.93,0.87); lt.SetBorderSize(0); lt.SetFillStyle(0); lt.SetTextSize(0.031)
    lt.AddEntry(hmc,"AAOgen all best candidates","l"); lt.AddEntry(ht,f"Truth #pi^{{0}} pair ({nt:,})","l")
    lt.AddEntry(hw,f"Wrong/combinatorial ({nw:,})","l"); lt.Draw()
    tt=ROOT.TLatex(); tt.SetNDC(True); tt.SetTextAlign(22); tt.SetTextSize(0.036)
    tt.DrawLatex(0.52,0.955,"AAOgen best-candidate #Delta p: direct MC truth decomposition")
    ts=ROOT.TLatex(); ts.SetNDC(True); ts.SetTextSize(0.026)
    ts.DrawLatex(0.15,0.88,"Truth: distinct matched MC photons with 0.125 < M_{#gamma#gamma}^{gen} < 0.145 GeV")
    truth_out=os.path.join(output_dir,f"7_{period}_AAOgen_delta_p_truth_decomposition.png"); ct.SaveAs(truth_out)
    keep.extend([ct,ht,hw,lt,tt,ts])
    print("\nAAOgen best-candidate Delta-p truth decomposition:")
    print("  Same best-candidate rule as the efficiency numerator.")
    print("  Truth = distinct matched MC photons with 0.125 < Mgg_gen < 0.145 GeV.")
    print(f"  truth pi0={nt:,}; wrong/combinatorial={nw:,}; truth fraction={nt/(nt+nw):.6f}" if nt+nw else "  no candidates")
    print(f"  Wrote truth decomposition to {truth_out}")

    # ------------------------------------------------------------------
    # Truth-template extraction.
    #
    # The Gaussian+polynomial decomposition above is retained only as a
    # diagnostic.  The nominal efficiency extraction below uses the measured
    # AAOgen truth-matched Delta-p shape as the signal template and the AAOgen
    # wrong/combinatorial best-candidate shape as the background template.
    # ------------------------------------------------------------------
    def hist_arrays(h):
        nb=h.GetNbinsX()
        x=np.asarray([h.GetBinCenter(i) for i in range(1,nb+1)],dtype=float)
        y=np.asarray([h.GetBinContent(i) for i in range(1,nb+1)],dtype=float)
        return x,y

    def poisson_deviance(obs, pred):
        obs=np.asarray(obs,dtype=float)
        pred=np.clip(np.asarray(pred,dtype=float),1.0e-12,None)
        terms=np.where(obs>0.0,2.0*(pred-obs+obs*np.log(obs/pred)),2.0*pred)
        return float(np.sum(terms))

    def normalized_template(v):
        v=np.clip(np.asarray(v,dtype=float),0.0,None)
        s=float(np.sum(v))
        return v/s if s>0.0 else np.zeros_like(v)

    def morph_template_shape(shape, shift_bins=0.0, sigma_bins=0.0):
        # _morph_1d_counts preserves the integral; templates passed here are
        # normalized to unit area, so the returned shape is also unit area.
        return normalized_template(_morph_1d_counts(
            normalized_template(shape),float(shift_bins),float(sigma_bins)))

    xbins,data_counts=hist_arrays(hdata)
    _,all_mc_counts=hist_arrays(hmc)
    _,truth_scaled=hist_arrays(ht)
    _,wrong_scaled=hist_arrays(hw)
    truth_shape=normalized_template(truth_scaled)
    wrong_shape=normalized_template(wrong_scaled)
    bin_width=float(hdata.GetBinWidth(1))

    def fit_template_mixture(obs, signal_shape, background_shape,
                             fit_lo=-1.0, fit_hi=1.0, allow_signal_morph=False):
        obs=np.asarray(obs,dtype=float)
        sig0=normalized_template(signal_shape)
        bg0=normalized_template(background_shape)
        mask=(xbins>=fit_lo-1.0e-12)&(xbins<=fit_hi+1.0e-12)
        if not np.any(mask):
            raise RuntimeError("Template fit has an empty fit range.")

        # Parameters are full-range component yields.  For a restricted fit
        # range, the expected counts in the fitted bins are simply the relevant
        # fraction of each full template.
        total=max(float(np.sum(obs[mask])),1.0)
        starts=[]
        if allow_signal_morph:
            for frac in (0.70,0.85,0.95):
                for sh in (0.0,5.0,10.0):
                    starts.append(np.asarray([frac*total,(1.0-frac)*total,sh,1.0]))
            bounds=[(0.0,5.0*total),(0.0,5.0*total),(-20.0,20.0),(0.0,20.0)]
        else:
            starts=[np.asarray([0.85*total,0.15*total]),
                    np.asarray([0.70*total,0.30*total]),
                    np.asarray([0.95*total,0.05*total])]
            bounds=[(0.0,5.0*total),(0.0,5.0*total)]

        def objective(par):
            ns=float(par[0]); nb=float(par[1])
            if allow_signal_morph:
                sig=morph_template_shape(sig0,float(par[2]),float(par[3]))
            else:
                sig=sig0
            pred=ns*sig+nb*bg0
            p=np.clip(pred[mask],1.0e-12,None)
            d=obs[mask]
            return float(np.sum(p-d*np.log(p)))

        fits=[minimize(objective,s,bounds=bounds,method="L-BFGS-B",
                       options={"maxiter":1000,"ftol":1.0e-12}) for s in starts]
        res=min(fits,key=lambda r:float(r.fun))
        par=np.asarray(res.x,dtype=float)
        sig=morph_template_shape(sig0,par[2],par[3]) if allow_signal_morph else sig0
        pred=par[0]*sig+par[1]*bg0
        dev=poisson_deviance(obs[mask],pred[mask])
        npar=4 if allow_signal_morph else 2
        ndf=max(int(np.count_nonzero(mask))-npar,1)
        return {
            "success":bool(res.success),"status":int(res.status),
            "message":str(res.message),"nsig":float(par[0]),"nbg":float(par[1]),
            "shift_bins":float(par[2]) if allow_signal_morph else 0.0,
            "sigma_bins":float(par[3]) if allow_signal_morph else 0.0,
            "shift_GeV":(float(par[2])*bin_width if allow_signal_morph else 0.0),
            "sigma_GeV":(float(par[3])*bin_width if allow_signal_morph else 0.0),
            "signal_shape":sig,"background_shape":bg0,"prediction":pred,
            "deviance":dev,"ndf":ndf,"deviance_ndf":dev/ndf,
            "fit_lo":float(fit_lo),"fit_hi":float(fit_hi)
        }

    # AAOgen closure is deliberately blind: the target is the total best-candidate
    # spectrum.  Truth labels enter only through the two fixed component shapes.
    closure=fit_template_mixture(all_mc_counts,truth_shape,wrong_shape,
                                 -1.0,1.0,allow_signal_morph=False)
    truth_known=float(np.sum(truth_scaled))
    wrong_known=float(np.sum(wrong_scaled))
    closure_truth_bias=(closure["nsig"]-truth_known)/truth_known if truth_known>0 else float("nan")
    closure_wrong_bias=(closure["nbg"]-wrong_known)/wrong_known if wrong_known>0 else float("nan")

    # Data nominal: allow the AAOgen truth signal template a translation and
    # additional Gaussian smearing.  The wrong/combinatorial template shape is
    # fixed; both component normalizations float.
    data_template=fit_template_mixture(data_counts,truth_shape,wrong_shape,
                                       -1.0,1.0,allow_signal_morph=True)
    data_nomorph=fit_template_mixture(data_counts,truth_shape,wrong_shape,
                                      -1.0,1.0,allow_signal_morph=False)

    # Range stability of the physically motivated template extraction.
    template_ranges=[(-0.50,0.60),(-0.70,0.80),(-1.00,1.00)]
    template_range_results=[]
    for flo,fhi in template_ranges:
        rr=fit_template_mixture(data_counts,truth_shape,wrong_shape,
                                flo,fhi,allow_signal_morph=True)
        template_range_results.append(rr)

    nd_now=int(den_data.GetValue())
    nm_now=int(den_mc.GetValue())
    eff_mc_truth=float(nt)/float(nm_now) if nm_now else float("nan")
    eff_data_template=data_template["nsig"]/float(nd_now) if nd_now else float("nan")
    Cgamma_template=(eff_data_template/eff_mc_truth
                     if eff_mc_truth>0.0 else float("nan"))

    # Utility for drawing a numpy prediction as a ROOT histogram.
    def array_hist(name, arr, color, style=1, width=2):
        hh=hdata.Clone(name); hh.SetDirectory(0); hh.Reset("ICES")
        for ib,val in enumerate(np.asarray(arr,dtype=float),1):
            hh.SetBinContent(ib,float(val)); hh.SetBinError(ib,0.0)
        hh.SetLineColor(color); hh.SetLineStyle(style); hh.SetLineWidth(width)
        hh.SetMarkerStyle(1)
        return hh

    # 8: blind AAOgen closure.
    cclosure=ROOT.TCanvas(f"c_dp_template_closure_{unique}","",1450,780)
    cclosure.SetTicks(1,1); cclosure.SetLeftMargin(0.12); cclosure.SetRightMargin(0.04)
    cclosure.SetBottomMargin(0.14); cclosure.SetTopMargin(0.11)
    hmc.SetStats(0); hmc.SetLineColor(ROOT.kBlack); hmc.SetLineWidth(2)
    hmc.GetXaxis().SetRangeUser(-1.0,1.0)
    hmc.GetXaxis().SetTitle("#Delta p = |p_{X}| - |p_{#gamma_{probe}}| (GeV)")
    hmc.GetYaxis().SetTitle("Normalized best-candidate rows"); hmc.GetYaxis().SetTitleOffset(1.35)
    hmc.Draw("HIST")
    hcl_sig=array_hist(f"h_cl_sig_{unique}",closure["nsig"]*closure["signal_shape"],ROOT.kBlue+1,1,3)
    hcl_bg=array_hist(f"h_cl_bg_{unique}",closure["nbg"]*closure["background_shape"],ROOT.kRed+1,2,3)
    hcl_tot=array_hist(f"h_cl_tot_{unique}",closure["prediction"],ROOT.kMagenta+2,1,3)
    hcl_tot.Draw("HIST SAME"); hcl_sig.Draw("HIST SAME"); hcl_bg.Draw("HIST SAME")
    lcl=ROOT.TLegend(0.57,0.66,0.93,0.87); lcl.SetBorderSize(0); lcl.SetFillStyle(0); lcl.SetTextSize(0.030)
    lcl.AddEntry(hmc,"AAOgen total (truth labels hidden)","l")
    lcl.AddEntry(hcl_tot,"Fitted truth + wrong templates","l")
    lcl.AddEntry(hcl_sig,"Fitted truth-template component","l")
    lcl.AddEntry(hcl_bg,"Fitted wrong-template component","l"); lcl.Draw()
    tx=ROOT.TLatex(); tx.SetNDC(True); tx.SetTextSize(0.029)
    tx.DrawLatex(0.15,0.92,"AAOgen blind template-yield closure")
    tx.DrawLatex(0.15,0.865,f"Known truth yield = {truth_known:.0f}; fitted = {closure['nsig']:.0f}")
    tx.DrawLatex(0.15,0.815,f"Truth-yield bias = {100.0*closure_truth_bias:+.2f}%")
    tx.DrawLatex(0.15,0.765,f"deviance/ndf = {closure['deviance_ndf']:.2f}")
    closure_out=os.path.join(output_dir,f"8_{period}_AAOgen_delta_p_template_closure.png")
    cclosure.SaveAs(closure_out)
    keep.extend([cclosure,hcl_sig,hcl_bg,hcl_tot,lcl,tx])

    # 9: Data truth-template extraction.
    ctemp=ROOT.TCanvas(f"c_dp_template_data_{unique}","",1450,780)
    ctemp.SetTicks(1,1); ctemp.SetLeftMargin(0.12); ctemp.SetRightMargin(0.04)
    ctemp.SetBottomMargin(0.14); ctemp.SetTopMargin(0.11)
    hdata.SetStats(0); hdata.SetLineColor(ROOT.kBlack); hdata.SetLineWidth(2)
    hdata.GetXaxis().SetRangeUser(-1.0,1.0)
    hdata.GetXaxis().SetTitle("#Delta p = |p_{X}| - |p_{#gamma_{probe}}| (GeV)")
    hdata.GetYaxis().SetTitle("Best-candidate rows"); hdata.GetYaxis().SetTitleOffset(1.35)
    hdata.Draw("E1")
    hdt_sig=array_hist(f"h_dt_sig_{unique}",data_template["nsig"]*data_template["signal_shape"],ROOT.kBlue+1,1,3)
    hdt_bg=array_hist(f"h_dt_bg_{unique}",data_template["nbg"]*data_template["background_shape"],ROOT.kRed+1,2,3)
    hdt_tot=array_hist(f"h_dt_tot_{unique}",data_template["prediction"],ROOT.kMagenta+2,1,3)
    hdt_tot.Draw("HIST SAME"); hdt_sig.Draw("HIST SAME"); hdt_bg.Draw("HIST SAME")
    ldt=ROOT.TLegend(0.57,0.66,0.93,0.87); ldt.SetBorderSize(0); ldt.SetFillStyle(0); ldt.SetTextSize(0.030)
    ldt.AddEntry(hdata,"Data","lep"); ldt.AddEntry(hdt_tot,"Template fit","l")
    ldt.AddEntry(hdt_sig,"Morphed AAOgen truth signal","l")
    ldt.AddEntry(hdt_bg,"AAOgen wrong/combinatorial","l"); ldt.Draw()
    td=ROOT.TLatex(); td.SetNDC(True); td.SetTextSize(0.028)
    td.DrawLatex(0.15,0.92,"Data #Delta p truth-template extraction")
    td.DrawLatex(0.15,0.87,f"signal shift = {data_template['shift_GeV']:+.4f} GeV, extra smear = {data_template['sigma_GeV']:.4f} GeV")
    td.DrawLatex(0.15,0.82,f"deviance/ndf = {data_template['deviance_ndf']:.2f}")
    td.DrawLatex(0.15,0.77,f"#epsilon_{{data}} = {eff_data_template:.4f}, #epsilon_{{AAO}}^{{truth}} = {eff_mc_truth:.4f}, C_{{#gamma}} = {Cgamma_template:.4f}")
    data_template_out=os.path.join(output_dir,f"9_{period}_data_delta_p_truth_template_fit.png")
    ctemp.SaveAs(data_template_out)
    keep.extend([ctemp,hdt_sig,hdt_bg,hdt_tot,ldt,td])

    # 10: compact template-method stability plot.
    cts=ROOT.TCanvas(f"c_dp_template_stability_{unique}","",1450,760)
    cts.SetTicks(1,1); cts.SetLeftMargin(0.12); cts.SetRightMargin(0.04)
    cts.SetBottomMargin(0.18); cts.SetTopMargin(0.11)
    fr=ROOT.TH1D(f"h_dp_template_stability_{unique}","",4,0.5,4.5); fr.SetDirectory(0); fr.SetStats(0)
    labels=["no morph [-1,1]","morph [-0.5,0.6]","morph [-0.7,0.8]","morph [-1,1]"]
    vals=[]
    ed0=data_nomorph["nsig"]/float(nd_now) if nd_now else float("nan")
    vals.append(ed0/eff_mc_truth if eff_mc_truth>0 else float("nan"))
    for rr in template_range_results:
        ee=rr["nsig"]/float(nd_now) if nd_now else float("nan")
        vals.append(ee/eff_mc_truth if eff_mc_truth>0 else float("nan"))
    finite=[v for v in vals if np.isfinite(v)]
    if finite:
        fr.SetMinimum(max(0.0,min(finite)-0.08)); fr.SetMaximum(max(finite)+0.08)
    else:
        fr.SetMinimum(0.0); fr.SetMaximum(1.2)
    fr.GetYaxis().SetTitle("C_{#gamma} = #epsilon_{data}/#epsilon_{AAOgen}^{truth}")
    fr.GetYaxis().SetTitleOffset(1.35); fr.GetXaxis().SetLabelSize(0.035)
    for i,lab in enumerate(labels,1): fr.GetXaxis().SetBinLabel(i,lab)
    fr.Draw("AXIS")
    gg=ROOT.TGraph()
    for i,v in enumerate(vals,1):
        if np.isfinite(v): gg.SetPoint(gg.GetN(),float(i),float(v))
    gg.SetMarkerStyle(20); gg.SetMarkerSize(1.4); gg.SetLineWidth(2); gg.Draw("PL SAME")
    tts=ROOT.TLatex(); tts.SetNDC(True); tts.SetTextAlign(22); tts.SetTextSize(0.036)
    tts.DrawLatex(0.53,0.955,"Truth-template #Delta p extraction stability")
    template_stability_out=os.path.join(output_dir,f"10_{period}_delta_p_truth_template_stability.png")
    cts.SaveAs(template_stability_out)
    keep.extend([cts,fr,gg,tts])

    print("\nAAOgen blind Delta-p template closure:")
    print(f"  known truth yield={truth_known:.1f}; fitted truth yield={closure['nsig']:.1f}; bias={100.0*closure_truth_bias:+.3f}%")
    print(f"  known wrong yield={wrong_known:.1f}; fitted wrong yield={closure['nbg']:.1f}; bias={100.0*closure_wrong_bias:+.3f}%")
    print(f"  deviance/ndf={closure['deviance_ndf']:.4f}; fit success={closure['success']} status={closure['status']}")
    print("\nData Delta-p truth-template extraction:")
    print(f"  signal yield={data_template['nsig']:.1f}; background yield={data_template['nbg']:.1f}")
    print(f"  signal shift={data_template['shift_GeV']:+.6f} GeV; extra smear={data_template['sigma_GeV']:.6f} GeV")
    print(f"  deviance/ndf={data_template['deviance_ndf']:.4f}; fit success={data_template['success']} status={data_template['status']}")
    print(f"  epsilon_data={eff_data_template:.6f}; epsilon_AAOgen_truth={eff_mc_truth:.6f}; C_gamma={Cgamma_template:.6f}")
    print("  Range stability:")
    for rr in template_range_results:
        ee=rr["nsig"]/float(nd_now) if nd_now else float("nan")
        cc=ee/eff_mc_truth if eff_mc_truth>0 else float("nan")
        print(f"    [{rr['fit_lo']:+.2f},{rr['fit_hi']:+.2f}] GeV: C_gamma={cc:.6f}, shift={rr['shift_GeV']:+.5f} GeV, smear={rr['sigma_GeV']:.5f} GeV, dev/ndf={rr['deviance_ndf']:.3f}")
    print(f"  Wrote closure to {closure_out}")
    print(f"  Wrote Data template fit to {data_template_out}")
    print(f"  Wrote template stability to {template_stability_out}")

    # ------------------------------------------------------------------
    # Legacy Gaussian+polynomial model/range study retained below ONLY as
    # a diagnostic demonstrating why that decomposition is not nominal.
    # ------------------------------------------------------------------
    # Signal/background model and fit-range stability.
    #
    # IMPORTANT: hdata/hmc are already the one-best-candidate-per-denominator-row
    # spectra.  We do NOT redo candidate selection for each fit variation.  This
    # isolates the uncertainty from the signal/background decomposition itself.
    # ------------------------------------------------------------------
    stability_ranges = [(-0.50,0.60), (-0.70,0.80), (-1.00,1.00)]
    stability_orders = (1,2,3)
    stability = []
    stability_keep = []

    def stability_fit(h, name, order, fit_lo, fit_hi):
        ax=h.GetXaxis()
        sb1=ax.FindBin(max(fit_lo,-0.45)+1e-9)
        sb2=ax.FindBin(min(fit_hi,+0.45)-1e-9)
        ib=max(range(sb1,sb2+1), key=lambda b:h.GetBinContent(b))
        mode=float(ax.GetBinCenter(ib))
        amp=max(float(h.GetBinContent(ib)),1.0)
        edge=0.5*(float(h.GetBinContent(ax.FindBin(fit_lo+0.03)))+
                  float(h.GetBinContent(ax.FindBin(fit_hi-0.03))))
        f=ROOT.TF1(name,f"gaus(0)+pol{order}(3)",fit_lo,fit_hi)
        pars=[amp,mode,0.12,max(edge,0.0)]+[0.0]*order
        f.SetParameters(*pars)
        f.SetParLimits(0,0.0,max(10.0*amp,1.0e9))
        f.SetParLimits(1,-0.35,0.45)
        f.SetParLimits(2,0.015,0.35)
        fitptr=h.Fit(f,"QNRS")
        try:
            status=int(fitptr)
        except Exception:
            status=-999
        mu=float(f.GetParameter(1))
        sigma=abs(float(f.GetParameter(2)))
        ndf=int(f.GetNDF())
        chi2=float(f.GetChisquare())
        valid=(status==0 and mu==mu and sigma==sigma and sigma>0.0 and
               abs(mu)<1e100 and abs(sigma)<1e100 and ndf>0)
        return f,mu,sigma,status,chi2,ndf,valid

    def stability_bg(f, order, name, fit_lo, fit_hi):
        bg=ROOT.TF1(name,f"pol{order}",fit_lo,fit_hi)
        bg.SetParameters(*[float(f.GetParameter(3+i)) for i in range(order+1)])
        return bg

    def stability_window(h,bg,mu,sigma,n,scale,fit_lo,fit_hi):
        # Never extrapolate the polynomial background outside the interval
        # actually used to determine it.
        req_lo=float(mu)-float(n)*float(sigma)
        req_hi=float(mu)+float(n)*float(sigma)
        lo=max(float(fit_lo),req_lo)
        hi=min(float(fit_hi),req_hi)
        if hi<=lo:
            return None
        ax=h.GetXaxis()
        b1=ax.FindBin(lo+1e-9)
        b2=ax.FindBin(hi-1e-9)
        obs=float(h.Integral(b1,b2))/float(scale)
        bkg=float(bg.Integral(lo,hi)/h.GetBinWidth(1))/float(scale)
        signal=max(0.0,obs-bkg)
        coverage=(hi-lo)/(req_hi-req_lo) if req_hi>req_lo else 0.0
        return obs,bkg,signal,coverage

    for order in stability_orders:
        for ir,(fit_lo,fit_hi) in enumerate(stability_ranges):
            fd,mud,sd,std,c2d,ndfd,vd=stability_fit(
                hdata,f"f_stab_data_p{order}_r{ir}_{unique}",order,fit_lo,fit_hi)
            fm,mum,sm,stm,c2m,ndfm,vm=stability_fit(
                hmc,f"f_stab_mc_p{order}_r{ir}_{unique}",order,fit_lo,fit_hi)
            bgd=stability_bg(fd,order,f"bg_stab_data_p{order}_r{ir}_{unique}",fit_lo,fit_hi)
            bgm=stability_bg(fm,order,f"bg_stab_mc_p{order}_r{ir}_{unique}",fit_lo,fit_hi)
            stability_keep.extend([fd,fm,bgd,bgm])
            item={"order":order,"fit_lo":fit_lo,"fit_hi":fit_hi,
                  "data_status":std,"mc_status":stm,
                  "data_chi2ndf":c2d/ndfd if ndfd>0 else float("nan"),
                  "mc_chi2ndf":c2m/ndfm if ndfm>0 else float("nan"),
                  "data_mu":mud,"data_sigma":sd,"mc_mu":mum,"mc_sigma":sm,
                  "valid":bool(vd and vm),"windows":{}}
            if item["valid"]:
                for n in (1,2,3):
                    wd=stability_window(hdata,bgd,mud,sd,n,1.0,fit_lo,fit_hi)
                    wm=stability_window(hmc,bgm,mum,sm,n,B,fit_lo,fit_hi)
                    if wd is None or wm is None:
                        continue
                    od,bd,sgd,covd=wd
                    om,bm,sgm,covm=wm
                    ed=sgd/nd if nd else float("nan")
                    em=sgm/nm if nm else float("nan")
                    item["windows"][n]={
                        "ratio":ed/em if em>0 else float("nan"),
                        "eff_data":ed,"eff_mc":em,
                        "coverage_data":covd,"coverage_mc":covm}
            stability.append(item)

    # One compact summary canvas.  Each x-bin is one background/range choice.
    cstab=ROOT.TCanvas(f"c_intdp_stability_{unique}","",1550,820)
    cstab.SetTicks(1,1); cstab.SetLeftMargin(0.11); cstab.SetRightMargin(0.04)
    cstab.SetBottomMargin(0.23); cstab.SetTopMargin(0.10)
    nchoices=len(stability)
    frame=ROOT.TH1D(f"h_intdp_stability_frame_{unique}","",nchoices,0.5,nchoices+0.5)
    frame.SetDirectory(0); frame.SetStats(0)
    frame.GetYaxis().SetTitle("C_{#gamma} = #epsilon_{data} / #epsilon_{AAOgen}")
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetLabelSize(0.028); frame.GetXaxis().LabelsOption("v")

    ratio_values=[]
    for i,item in enumerate(stability,1):
        frame.GetXaxis().SetBinLabel(
            i,f"pol{item['order']} [{item['fit_lo']:+.1f},{item['fit_hi']:+.1f}]")
        if item["valid"]:
            for n in (1,2,3):
                if n in item["windows"]:
                    v=item["windows"][n]["ratio"]
                    if v==v: ratio_values.append(v)
    if ratio_values:
        ymin=max(0.0,min(ratio_values)-0.08); ymax=max(ratio_values)+0.08
        if ymax-ymin<0.20:
            mid=0.5*(ymin+ymax); ymin=max(0.0,mid-0.10); ymax=mid+0.10
    else:
        ymin,ymax=0.5,1.1
    frame.SetMinimum(ymin); frame.SetMaximum(ymax); frame.Draw("AXIS")

    specs={1:(ROOT.kBlue+1,20),2:(ROOT.kRed+1,21),3:(ROOT.kGreen+2,22)}
    stab_graphs={}
    for n in (1,2,3):
        g=ROOT.TGraph(); g.SetName(f"g_intdp_stability_{n}s_{unique}")
        ip=0
        for i,item in enumerate(stability,1):
            if item["valid"] and n in item["windows"]:
                v=item["windows"][n]["ratio"]
                if v==v:
                    g.SetPoint(ip,float(i),float(v)); ip+=1
        col,marker=specs[n]
        g.SetMarkerColor(col); g.SetLineColor(col); g.SetMarkerStyle(marker)
        g.SetMarkerSize(1.25); g.SetLineWidth(2); g.Draw("PL SAME")
        stab_graphs[n]=g

    legstab=ROOT.TLegend(0.14,0.74,0.34,0.88)
    legstab.SetBorderSize(0); legstab.SetFillStyle(0); legstab.SetTextSize(0.031)
    for n in (1,2,3): legstab.AddEntry(stab_graphs[n],f"{n}#sigma window","lp")
    legstab.Draw()
    titlestab=ROOT.TLatex(); titlestab.SetNDC(True); titlestab.SetTextAlign(22)
    titlestab.SetTextSize(0.036)
    titlestab.DrawLatex(0.53,0.955,
        "Integrated #Delta p efficiency: background-model and fit-range stability")
    notestab=ROOT.TLatex(); notestab.SetNDC(True); notestab.SetTextSize(0.025)
    notestab.DrawLatex(0.39,0.86,
        "Best-candidate population fixed; only signal/background fit changes")

    stability_out=os.path.join(
        output_dir,f"11_{period}_LEGACY_integrated_delta_p_model_range_stability.png")
    cstab.SaveAs(stability_out)
    keep.extend(stability_keep+[cstab,frame,legstab,titlestab,notestab]+list(stab_graphs.values()))

    # ------------------------------------------------------------------
    # Visual inspection grids: one 3x3 canvas for Data and one for AAOgen.
    # Rows are polynomial orders 1,2,3; columns are the three fit ranges.
    # Each pad shows the SAME best-candidate histogram used in the numerical
    # stability scan, together with total fit, Gaussian signal, and polynomial
    # background.  This makes pathological decompositions immediately visible.
    # ------------------------------------------------------------------
    def draw_stability_fit_grid(sample_name, hsource, fit_key_prefix, color, output_name):
        cgrid=ROOT.TCanvas(f"c_intdp_grid_{sample_name}_{unique}","",1800,1500)
        cgrid.Divide(3,3,0.001,0.001)
        grid_keep=[cgrid]
        for ipad,item in enumerate(stability,1):
            pad=cgrid.cd(ipad)
            pad.SetTicks(1,1)
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.035)
            pad.SetBottomMargin(0.14)
            pad.SetTopMargin(0.12)

            # Clone the common best-candidate spectrum so every pad owns a
            # separately styled drawable histogram.
            hh=hsource.Clone(f"h_grid_{sample_name}_{ipad}_{unique}")
            hh.SetDirectory(0)
            hh.SetStats(0)
            hh.SetLineColor(color)
            hh.SetMarkerColor(color)
            hh.SetMarkerStyle(20)
            hh.SetMarkerSize(0.38)
            hh.GetXaxis().SetRangeUser(-1.0,1.0)
            hh.GetXaxis().SetTitle("#Delta p = |p_{X}| - |p_{#gamma_{probe}}| (GeV)")
            hh.GetYaxis().SetTitle("Normalized candidate combinations")
            hh.GetXaxis().SetTitleSize(0.045)
            hh.GetYaxis().SetTitleSize(0.045)
            hh.GetXaxis().SetLabelSize(0.037)
            hh.GetYaxis().SetLabelSize(0.037)
            hh.GetYaxis().SetTitleOffset(1.45)
            hh.Draw("E1")

            order=item["order"]
            fit_lo=item["fit_lo"]
            fit_hi=item["fit_hi"]
            if sample_name=="data":
                ftotal=next(obj for obj in stability_keep
                            if obj.GetName()==f"f_stab_data_p{order}_r{(ipad-1)%3}_{unique}")
                bg=next(obj for obj in stability_keep
                        if obj.GetName()==f"bg_stab_data_p{order}_r{(ipad-1)%3}_{unique}")
                mu=item["data_mu"]; sigma=item["data_sigma"]
                chi2ndf=item["data_chi2ndf"]; status=item["data_status"]
            else:
                ftotal=next(obj for obj in stability_keep
                            if obj.GetName()==f"f_stab_mc_p{order}_r{(ipad-1)%3}_{unique}")
                bg=next(obj for obj in stability_keep
                        if obj.GetName()==f"bg_stab_mc_p{order}_r{(ipad-1)%3}_{unique}")
                mu=item["mc_mu"]; sigma=item["mc_sigma"]
                chi2ndf=item["mc_chi2ndf"]; status=item["mc_status"]

            # Build the Gaussian component directly from this variation's total fit.
            sg=ROOT.TF1(f"sig_grid_{sample_name}_{ipad}_{unique}","gaus",fit_lo,fit_hi)
            sg.SetParameters(float(ftotal.GetParameter(0)),
                             float(ftotal.GetParameter(1)),
                             float(ftotal.GetParameter(2)))

            ftotal.SetLineColor(ROOT.kMagenta+2)
            ftotal.SetLineWidth(3)
            sg.SetLineColor(ROOT.kBlue)
            sg.SetLineWidth(2)
            bg.SetLineColor(ROOT.kGray+2)
            bg.SetLineStyle(2)
            bg.SetLineWidth(2)
            ftotal.Draw("SAME")
            sg.Draw("SAME")
            bg.Draw("SAME")

            # Mark the actual fit boundaries.  This is particularly useful for
            # judging why the fitted background changes as the range is widened.
            ymax=max(float(hh.GetMaximum())*1.08,1.0)
            llo=ROOT.TLine(fit_lo,0.0,fit_lo,ymax)
            lhi=ROOT.TLine(fit_hi,0.0,fit_hi,ymax)
            for line in (llo,lhi):
                line.SetLineColor(ROOT.kGray+1)
                line.SetLineStyle(3)
                line.SetLineWidth(1)
                line.Draw()

            lab=ROOT.TLatex()
            lab.SetNDC(True)
            lab.SetTextSize(0.033)
            lab.DrawLatex(0.17,0.935,
                          f"{sample_name.upper()}  pol{order}, [{fit_lo:+.1f},{fit_hi:+.1f}] GeV")
            lab.DrawLatex(0.17,0.885,
                          f"#mu={mu:+.4f} GeV, #sigma={sigma:.4f} GeV")
            lab.DrawLatex(0.17,0.835,
                          f"#chi^{{2}}/ndf={chi2ndf:.2f}, status={status}")
            if item["valid"] and 2 in item["windows"]:
                lab.DrawLatex(0.17,0.785,
                              f"C_{{#gamma}}(2#sigma)={item['windows'][2]['ratio']:.4f}")

            grid_keep.extend([hh,sg,llo,lhi,lab])

        cgrid.cd()
        cgrid.Update()
        cgrid.SaveAs(output_name)
        return cgrid,grid_keep

    data_grid_out=os.path.join(
        output_dir,f"12_{period}_LEGACY_integrated_delta_p_fit_grid_data.png")
    mc_grid_out=os.path.join(
        output_dir,f"13_{period}_LEGACY_integrated_delta_p_fit_grid_AAOgen.png")
    cgrid_data,grid_data_keep=draw_stability_fit_grid(
        "data",hdata,"data",ROOT.kBlack,data_grid_out)
    cgrid_mc,grid_mc_keep=draw_stability_fit_grid(
        "aaogen",hmc,"mc",ROOT.kRed+1,mc_grid_out)
    keep.extend(grid_data_keep+grid_mc_keep)

    print("\\nIntegrated Delta-p model/range stability:")
    print("  Scan = pol1/pol2/pol3 backgrounds x fit ranges [-0.5,+0.6], [-0.7,+0.8], [-1.0,+1.0] GeV.")
    print("  Best-candidate population is held fixed for every variation.")
    print("  Polynomial background is integrated only inside its fitted range; coverage is printed when an n-sigma window is truncated.")
    print("  model/range | data chi2/ndf | AAO chi2/ndf | Cgamma(1s) Cgamma(2s) Cgamma(3s)")
    for item in stability:
        label=f"pol{item['order']} [{item['fit_lo']:+.2f},{item['fit_hi']:+.2f}]"
        if not item["valid"]:
            print(f"  {label} | FIT FAILED data={item['data_status']} AAO={item['mc_status']}")
            continue
        vals=[]
        for n in (1,2,3):
            if n not in item["windows"]:
                vals.append("NA"); continue
            w=item["windows"][n]
            vals.append(f"{w['ratio']:.6f}[cov {100*w['coverage_data']:.0f}/{100*w['coverage_mc']:.0f}%]")
        print(f"  {label} | {item['data_chi2ndf']:.3f} | {item['mc_chi2ndf']:.3f} | "+" ".join(vals))
    for n in (1,2,3):
        vals=[item["windows"][n]["ratio"] for item in stability
              if item["valid"] and n in item["windows"]
              and item["windows"][n]["ratio"]==item["windows"][n]["ratio"]]
        if vals:
            mean=sum(vals)/len(vals)
            rms=(sum((v-mean)**2 for v in vals)/len(vals))**0.5
            print(f"  {n}sigma stability: mean={mean:.6f}, min={min(vals):.6f}, max={max(vals):.6f}, span={max(vals)-min(vals):.6f}, RMS={rms:.6f}")
    print(f"  Wrote stability summary to {stability_out}")
    results["truth_output"]=truth_out
    results["stability"]=stability
    results["stability_output"]=stability_out
    results["stability_data_grid_output"]=data_grid_out
    results["stability_mc_grid_output"]=mc_grid_out
    print(f"  Wrote Data 3x3 fit grid to {data_grid_out}")
    print(f"  Wrote AAOgen 3x3 fit grid to {mc_grid_out}")

    # Promote the truth-template result to the nominal integrated result.
    results["nominal_method"]="AAOgen truth-template Delta-p extraction"
    results["nominal_eff_data"]=eff_data_template
    results["nominal_eff_mc"]=eff_mc_truth
    results["nominal_Cgamma"]=Cgamma_template
    results["template_closure"]=closure
    results["template_data_fit"]=data_template
    results["template_range_results"]=template_range_results
    results["template_closure_output"]=closure_out
    results["template_data_output"]=data_template_out
    results["template_stability_output"]=template_stability_out

    print("\nIntegrated Delta-p efficiency study (NO kinematic binning):")
    print("  NOMINAL method = AAOgen truth-matched Delta-p signal template + AAOgen wrong/combinatorial template.")
    print("  Candidate ranking = smallest |Delta p|; it does not depend on a fitted Gaussian peak.")
    print("  AAOgen efficiency numerator = directly truth-matched best candidates; no Gaussian window and no background subtraction.")
    print("  Data signal template may shift and receive additional Gaussian smearing; signal/background normalizations float.")
    print("  Legacy Gaussian+polynomial fits are diagnostic only and do NOT define the nominal C_gamma.")
    print(f"  Denominator rows: data={nd_now:,}; AAOgen={nm_now:,}")
    print(f"  Rows with >=1 probe candidate: data={int(cand_data.GetValue()):,}; AAOgen={int(cand_mc.GetValue()):,}")
    print(f"  AAOgen truth-matched best candidates={nt:,}; epsilon_AAOgen_truth={eff_mc_truth:.6f}")
    print(f"  Data fitted signal={data_template['nsig']:.1f}; epsilon_data={eff_data_template:.6f}")
    print(f"  NOMINAL C_gamma={Cgamma_template:.6f}")
    return keep,out,results

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
        pad.SetLeftMargin(0.17)
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
        hp.GetYaxis().SetTitleOffset(1.75)
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
        y0 = 0.84
        for il, line in enumerate(extra_lines):
            info.DrawLatex(0.17, y0 - 0.045 * il, line)

        leg = ROOT.TLegend(0.56, 0.64, 0.95, 0.88)
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
        "single_params": [f_single.GetParameter(i) for i in range(6)],
        "double_params": [f_double.GetParameter(i) for i in range(8)],
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


def draw_probe_aaogen_truth_closure(dfs, output_dir, period, coeffs, fit_result):
    """Truth-label the *same selected AAOgen tag-probe combinations* used by the fit.

    This is the closure that is actually relevant to the fitted pi0 yield:
      - start from the same PhotonEfficiency AAOgen rows as plot 2;
      - apply exactly _exclusive_df(), hence the identical epgammaX cuts;
      - build exactly the same tag + retained-PID22-partner combinations;
      - use saved MC::Particle/REC matching only to label each combination as
        a true pi0 daughter pair or a combinatorial/wrong pair.

    A true pair is defined directly from the two reconstructed photons' own
    REC-to-MC associations: both must match distinct generated photons, and the
    invariant mass of those two matched generated photons must lie in the pi0
    mass window 0.125--0.145 GeV.  Parent-PDG and saved "probe role" fields are
    diagnosed separately but are not required for the truth label.

    The inclusive selected spectrum is therefore decomposed bin-by-bin as
        inclusive = truth-matched pi0 + combinatorial/wrong.
    All three spectra are scaled by the same AAOgen normalization B used in
    plots 1--3.  Closure is evaluated over the nominal fit interval, exactly
    matching the integration interval used for the fitted signal yields.
    """
    if coeffs is None or fit_result is None or "aaogen" not in dfs:
        print("WARNING: AAOgen truth-closure inputs unavailable; skipping.")
        return [], None, None

    _, B, _ = coeffs["mean"]
    lo = float(fit_result["nominal_lo"])
    hi = float(fit_result["nominal_hi"])
    unique = str(abs(hash((period, "aaogen_selected_pair_truth_closure"))))

    # Exact same selected rows and exact same inclusive Mgg construction as plot 2.
    df = (_exclusive_df(dfs["aaogen"], "truthclosure_selected_aaogen")
          .Define("closure_Mgg_all",
              "pe_mgg_tag_probe(tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,"
              "neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi)")
          .Define("closure_Mgg_true",
              "pe_mgg_truth_by_matched_mc("
              "tag_corr_p,tag_corr_theta,tag_corr_phi,tag_rec_index,"
              "mc_tag_index,mc_tag_pid,mc_tag_p,mc_tag_theta,mc_tag_phi,"
              "neutral_idx,neutral_pid,neutral_p,neutral_theta,neutral_phi,"
              "neutral_mc_index,neutral_mc_pid,neutral_mc_p,neutral_mc_theta,neutral_mc_phi)")
          .Define("closure_match_stage",
              "pe_mgg_truth_match_stage(tag_rec_index,mc_tag_index,mc_tag_pid,mc_tag_p,"
              "neutral_idx,neutral_pid,neutral_mc_index,neutral_mc_pid,neutral_mc_p)")
          .Define("closure_role_mask",
              "pe_mgg_truth_role_diagnostics(tag_rec_index,mc_tag_index,mc_tag_pid,mc_tag_parent,"
              "mc_probe_index,mc_probe_pid,mc_probe_parent,neutral_idx,neutral_pid,"
              "neutral_mc_index,neutral_mc_pid)"))

    hall_ptr = df.Histo1D(
        (f"h_closure_all_{unique}",
         ";M_{#gamma_{tag}#gamma_{probe}} (GeV);Normalized tag-probe combinations",
         160, 0.0, 0.8), "closure_Mgg_all")
    htrue_ptr = df.Histo1D(
        (f"h_closure_true_{unique}",
         ";M_{#gamma_{tag}#gamma_{probe}} (GeV);Normalized tag-probe combinations",
         160, 0.0, 0.8), "closure_Mgg_true")
    hstage_ptr = df.Histo1D((f"h_closure_stage_{unique}", ";stage;combinations", 8, 0.5, 8.5), "closure_match_stage")
    hrole_ptr = df.Histo1D((f"h_closure_role_{unique}", ";mask;combinations", 256, -0.5, 255.5), "closure_role_mask")
    ROOT.RDF.RunGraphs([hall_ptr, htrue_ptr, hstage_ptr, hrole_ptr])

    hall = hall_ptr.GetValue().Clone(f"h_closure_all_scaled_{unique}")
    htrue = htrue_ptr.GetValue().Clone(f"h_closure_true_scaled_{unique}")
    hall.SetDirectory(0); htrue.SetDirectory(0)
    hall.Sumw2(); htrue.Sumw2()
    hall.Scale(B); htrue.Scale(B)

    hbg = hall.Clone(f"h_closure_bg_scaled_{unique}")
    hbg.SetDirectory(0)
    hbg.Add(htrue, -1.0)

    # Same bin convention as the fit. The fit yields are integrals over [lo,hi].
    ax = hall.GetXaxis()
    blo = ax.FindBin(lo + 1.0e-9)
    bhi = ax.FindBin(hi - 1.0e-9)
    truth_yield = float(htrue.Integral(blo, bhi))
    inclusive_yield = float(hall.Integral(blo, bhi))
    bg_yield = float(hbg.Integral(blo, bhi))

    fit_single = float(fit_result["single_yield"])
    fit_double = float(fit_result["double_yield"])
    bias_single = 100.0 * (fit_single / truth_yield - 1.0) if truth_yield > 0 else float("nan")
    bias_double = 100.0 * (fit_double / truth_yield - 1.0) if truth_yield > 0 else float("nan")

    # Reconstruct exactly the nominal double-Gaussian fit components from plot 2.
    dp = fit_result["double_params"]
    f_sig = ROOT.TF1(
        f"f_closure_double_signal_{unique}",
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)",
        lo, hi)
    for i in range(5):
        f_sig.SetParameter(i, dp[i])

    f_bg = ROOT.TF1(f"f_closure_double_bg_{unique}", "pol2", lo, hi)
    for i in range(3):
        f_bg.SetParameter(i, dp[5+i])

    # Canvas: left = truth decomposition; right = direct truth-vs-fit closure.
    c = ROOT.TCanvas(f"c_selected_truthclosure_{unique}", "", 1500, 720)
    c.Divide(2, 1, 0.002, 0.002)
    keep = [c, hall_ptr, htrue_ptr, hstage_ptr, hrole_ptr, hall, htrue, hbg, f_sig, f_bg]

    # Left panel.
    p = c.cd(1)
    p.SetTicks(1,1); p.SetLeftMargin(0.16); p.SetRightMargin(0.04)
    p.SetBottomMargin(0.14); p.SetTopMargin(0.11)

    hall.SetStats(0)
    hall.SetLineColor(ROOT.kRed + 1); hall.SetLineWidth(2)
    hall.SetMarkerColor(ROOT.kRed + 1); hall.SetMarkerStyle(20); hall.SetMarkerSize(0.55)
    hall.GetXaxis().SetRangeUser(0.0, 0.30)
    hall.GetXaxis().SetTitleSize(0.048); hall.GetYaxis().SetTitleSize(0.043)
    hall.GetYaxis().SetTitleOffset(1.60)
    hall.SetMaximum(1.18 * max(hall.GetMaximum(), 1.0))
    hall.Draw("E1")

    htrue.SetLineColor(ROOT.kBlue + 1); htrue.SetLineWidth(3)
    htrue.SetMarkerColor(ROOT.kBlue + 1); htrue.SetMarkerStyle(24); htrue.SetMarkerSize(0.50)
    htrue.Draw("E1 SAME")

    hbg.SetLineColor(ROOT.kGray + 2); hbg.SetLineWidth(2); hbg.SetLineStyle(2)
    hbg.Draw("HIST SAME")
    hall.Draw("E1 SAME")

    leg1 = ROOT.TLegend(0.48, 0.68, 0.94, 0.87)
    leg1.SetBorderSize(0); leg1.SetFillStyle(0); leg1.SetTextSize(0.029)
    leg1.AddEntry(hall, "All selected AAOgen combinations", "lep")
    leg1.AddEntry(htrue, "Truth-matched #pi^{0} pairs", "lep")
    leg1.AddEntry(hbg, "Combinatorial / wrong pairs", "l")
    leg1.Draw()

    title1 = ROOT.TLatex(); title1.SetNDC(True); title1.SetTextAlign(22); title1.SetTextSize(0.036)
    title1.DrawLatex(0.54, 0.955, "Same selected AAOgen combinations: truth decomposition")
    keep += [leg1, title1]

    # Right panel: truth signal against the actual fitted signal/background.
    p = c.cd(2)
    p.SetTicks(1,1); p.SetLeftMargin(0.16); p.SetRightMargin(0.04)
    p.SetBottomMargin(0.14); p.SetTopMargin(0.11)

    htrue2 = htrue.Clone(f"h_closure_true_compare_{unique}")
    htrue2.SetDirectory(0); htrue2.SetStats(0)
    htrue2.GetXaxis().SetRangeUser(0.0, 0.30)
    htrue2.GetYaxis().SetTitleOffset(1.60)
    htrue2.SetMaximum(1.25 * max(htrue.GetMaximum(), f_sig.GetMaximum(), 1.0))
    htrue2.Draw("E1")

    f_sig.SetLineColor(ROOT.kMagenta + 2); f_sig.SetLineWidth(3)
    f_bg.SetLineColor(ROOT.kGray + 2); f_bg.SetLineWidth(2); f_bg.SetLineStyle(2)
    f_sig.Draw("SAME")
    f_bg.Draw("SAME")

    leg2 = ROOT.TLegend(0.53, 0.67, 0.94, 0.87)
    leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.028)
    leg2.AddEntry(htrue2, "Truth-matched #pi^{0} pairs", "lep")
    leg2.AddEntry(f_sig, "Double-Gaussian fitted signal", "l")
    leg2.AddEntry(f_bg, "Quadratic fitted background", "l")
    leg2.Draw()

    info = ROOT.TLatex(); info.SetNDC(True); info.SetTextSize(0.027)
    info.DrawLatex(0.19, 0.84, f"{lo:.2f} < M_{{#gamma#gamma}} < {hi:.2f} GeV")
    info.DrawLatex(0.19, 0.79, f"Truth N_{{#pi^{{0}}}} = {truth_yield:,.0f}")
    info.DrawLatex(0.19, 0.74, f"Single-G N_{{#pi^{{0}}}} = {fit_single:,.0f}  ({bias_single:+.2f}%)")
    info.DrawLatex(0.19, 0.69, f"Double-G N_{{#pi^{{0}}}} = {fit_double:,.0f}  ({bias_double:+.2f}%)")

    title2 = ROOT.TLatex(); title2.SetNDC(True); title2.SetTextAlign(22); title2.SetTextSize(0.036)
    title2.DrawLatex(0.54, 0.955, "Fit closure against known AAOgen truth")
    keep += [htrue2, leg2, info, title2]

    out = os.path.join(output_dir, f"4_{period}_AAOgen_Mgg_pi0_truth_closure.png")
    c.SaveAs(out)

    hstage = hstage_ptr.GetValue()
    hrole = hrole_ptr.GetValue()
    stage_exact = {i: int(round(hstage.GetBinContent(i))) for i in range(1,9)}
    # Because each combination stores the highest stage reached, cumulative N(stage>=k)
    # is the sum from k through 8.
    stage_cumulative = {k: sum(stage_exact[i] for i in range(k,9)) for k in range(1,9)}

    role_counts = {"tag_index":0,"tag_pid22":0,"tag_parent111":0,"probe_index":0,
                   "probe_pid22":0,"probe_parent111":0,"partner_eq_probe":0,"partner_pid22":0,
                   "all_old_requirements":0}
    for ib in range(1, hrole.GetNbinsX()+1):
        n=int(round(hrole.GetBinContent(ib)))
        if not n: continue
        mask=int(round(hrole.GetBinCenter(ib)))
        if mask & 1: role_counts["tag_index"] += n
        if mask & 2: role_counts["tag_pid22"] += n
        if mask & 4: role_counts["tag_parent111"] += n
        if mask & 8: role_counts["probe_index"] += n
        if mask & 16: role_counts["probe_pid22"] += n
        if mask & 32: role_counts["probe_parent111"] += n
        if mask & 64: role_counts["partner_eq_probe"] += n
        if mask & 128: role_counts["partner_pid22"] += n
        if mask == 255: role_counts["all_old_requirements"] += n

    result = {
        "truth_yield": truth_yield,
        "inclusive_yield": inclusive_yield,
        "background_yield": bg_yield,
        "single_fit_yield": fit_single,
        "double_fit_yield": fit_double,
        "single_bias_pct": bias_single,
        "double_bias_pct": bias_double,
        "raw_all_pairs": int(round(hall_ptr.GetValue().GetEntries())),
        "raw_truth_pairs": int(round(htrue_ptr.GetValue().GetEntries())),
        "match_stage_cumulative": stage_cumulative,
        "old_role_counts": role_counts,
    }

    diag_out = os.path.join(output_dir, f"4b_{period}_AAOgen_Mgg_truth_matching_diagnostics.txt")
    with open(diag_out, "w") as fout:
        fout.write("AAOgen selected reconstructed-pair truth-matching diagnostics\n")
        fout.write("============================================================\n\n")
        fout.write("Population: exactly the same _exclusive_df() AAOgen rows and tag + retained PID22 partner combinations used in probe plot 2.\n\n")
        labels={1:"reconstructed tag + reconstructed PID22 partner exists",
                2:"tag has a saved MC particle index",
                3:"tag's matched MC particle is a photon (PID 22)",
                4:"tag's matched MC photon has a saved positive momentum",
                5:"partner has a saved MC particle index",
                6:"partner's matched MC particle is a photon (PID 22)",
                7:"partner's matched MC photon has a saved positive momentum",
                8:"tag and partner are matched to two distinct MC particles"}
        fout.write("Direct REC -> MC association chain (cumulative counts):\n")
        for k in range(1,9): fout.write(f"  stage {k}: {stage_cumulative[k]:,}  {labels[k]}\n")
        fout.write("\nOld role/parent assumptions, tested independently (NOT used by new truth definition):\n")
        for key,val in role_counts.items(): fout.write(f"  {key:24s} {val:,}\n")
        fout.write("\nNew truth definition: after stage 8, compute Mgg from the two matched generated-photon four-vectors; call the reconstructed pair a true pi0 pair when 0.125 < generated Mgg < 0.145 GeV.\n")
        fout.write(f"Raw truth-labeled reconstructed combinations: {result['raw_truth_pairs']:,}\n")
        fout.write(f"B-scaled truth yield in {lo:.3f}--{hi:.3f} GeV reconstructed fit interval: {truth_yield:.3f}\n")
        fout.write(f"Single-G fit yield: {fit_single:.3f}; closure bias: {bias_single:+.3f}%\n")
        fout.write(f"Double-G fit yield: {fit_double:.3f}; closure bias: {bias_double:+.3f}%\n")

    print("\nAAOgen selected-pair pi0 truth closure:")
    print("  Population: SAME selected reconstructed AAOgen tag-partner combinations as plot 2.")
    print("  New truth label: match each reconstructed photon independently to its saved MC photon;")
    print("                   require two distinct matched MC photons; compute their generated Mgg;")
    print("                   generated 0.125 < Mgg < 0.145 GeV => true pi0 pair.")
    print("  Parent==111 and saved mc_probe_index are NOT used to define truth.")
    print(f"  nominal reconstructed fit interval = {lo:.3f}--{hi:.3f} GeV")
    print("\n  REC -> MC matching chain (cumulative):")
    stage_labels=["reconstructed pair", "tag has MC index", "tag MC PID=22", "tag MC momentum saved",
                  "partner has MC index", "partner MC PID=22", "partner MC momentum saved", "distinct MC indices"]
    for k,label in enumerate(stage_labels,1):
        print(f"    {k}. {label:29s}: {stage_cumulative[k]:,}")
    print("\n  Diagnostic of the assumptions that produced the previous zero:")
    print(f"    tag mc_parent == 111:                 {role_counts['tag_parent111']:,}")
    print(f"    saved probe MC index exists:          {role_counts['probe_index']:,}")
    print(f"    saved probe MC PID == 22:             {role_counts['probe_pid22']:,}")
    print(f"    saved probe mc_parent == 111:         {role_counts['probe_parent111']:,}")
    print(f"    reconstructed partner == saved probe: {role_counts['partner_eq_probe']:,}")
    print(f"    ALL old requirements simultaneously:  {role_counts['all_old_requirements']:,}")
    print("\n  New direct generated-mass truth result:")
    print(f"    raw selected combinations = {result['raw_all_pairs']:,}")
    print(f"    raw truth-matched pi0 combinations = {result['raw_truth_pairs']:,}")
    print(f"    B-scaled inclusive combinations in fit interval = {inclusive_yield:.2f}")
    print(f"    B-scaled truth pi0 combinations in fit interval = {truth_yield:.2f}")
    print(f"    B-scaled combinatorial/wrong combinations in fit interval = {bg_yield:.2f}")
    print(f"    single-Gaussian fitted signal = {fit_single:.2f}  closure bias = {bias_single:+.3f}%")
    print(f"    double-Gaussian fitted signal = {fit_double:.2f}  closure bias = {bias_double:+.3f}%")
    print(f"  Detailed matching diagnostic -> {diag_out}")

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
    deltap_keep, deltap_output = draw_probe_delta_p(
        dfs, probe_dir, args.period, denominator_coeffs)
    keep.extend(deltap_keep)
    inteff_keep, inteff_output, inteff_result = draw_probe_integrated_delta_p_efficiency(
        dfs, probe_dir, args.period, denominator_coeffs)
    keep.extend(inteff_keep)
    keep.extend(probe_keep)
    pi0fit_keep, pi0fit_output, pi0fit_result = draw_probe_aaogen_pi0_fit(
        dfs, probe_dir, args.period, denominator_coeffs)
    keep.extend(pi0fit_keep)

    truthclosure_keep, truthclosure_output, truthclosure_result = draw_probe_aaogen_truth_closure(
        dfs, probe_dir, args.period, denominator_coeffs, pi0fit_result)
    keep.extend(truthclosure_keep)
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
    if deltap_output:
        print(f"\nWrote: {deltap_output}")
    if inteff_output:
        print(f"\nWrote: {inteff_output}")
    if inteff_result and inteff_result.get("truth_output"):
        print(f"\nWrote: {inteff_result['truth_output']}")
    if inteff_result and inteff_result.get("template_closure_output"):
        print(f"\nWrote: {inteff_result['template_closure_output']}")
    if inteff_result and inteff_result.get("template_data_output"):
        print(f"\nWrote: {inteff_result['template_data_output']}")
    if inteff_result and inteff_result.get("template_stability_output"):
        print(f"\nWrote: {inteff_result['template_stability_output']}")
    if inteff_result and inteff_result.get("stability_output"):
        print(f"\nWrote: {inteff_result['stability_output']}")
    if inteff_result and inteff_result.get("stability_data_grid_output"):
        print(f"\nWrote: {inteff_result['stability_data_grid_output']}")
    if inteff_result and inteff_result.get("stability_mc_grid_output"):
        print(f"\nWrote: {inteff_result['stability_mc_grid_output']}")
    if pi0fit_output:
        print(f"\nWrote: {pi0fit_output}")
    if pi0fit_result and pi0fit_result.get("window_scan_output"):
        print(f"\nWrote: {pi0fit_result['window_scan_output']}")
    if truthclosure_output:
        print(f"\nWrote: {truthclosure_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
