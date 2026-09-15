#!/usr/bin/env python3
"""
CLAS12 photon-efficiency analysis -- restart from the raw e'p'gamma1 hypotheses.

Current step
------------
Use every row in the PhotonEfficiency hypothesis tree, requiring only W > 2 GeV.
Each row is one reconstructed e'p'gamma1 hypothesis produced by the skim.

Plot the following observables as rows:
  1) E_gamma1                     1 to 10 GeV
  2) missing energy e'p'gamma1    1 to 10 GeV
  3) Mx2(e'p')                   -0.5 to 2.0 GeV^2
  4) Mx2(e'p'gamma1)             -0.1 to 0.2 GeV^2

Columns are cumulative selections:
  1) W > 2 GeV only
  2) additionally -0.05 < Mx2(e'p'gamma1) < 0.05 GeV^2
  3) additionally Mx2(e'p') < 0.25 GeV^2

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

double pe_emiss_epg(double ebeam, double ep, double pp, double gp) {
    const double me = 0.00051099895;
    const double mp = 0.9382720813;

    const double Ee = std::sqrt(ep*ep + me*me);
    const double Ep = std::sqrt(pp*pp + mp*mp);

    // Missing energy for e p -> e' p' gamma1 X:
    // E_miss = E_beam + M_p - E_e' - E_p' - E_gamma1.
    return ebeam + mp - Ee - Ep - gp;
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

    No analysis selection is imposed here except W > 2 GeV.  In particular,
    there is no beta, photon fiducial, electron-photon opening-angle,
    missing-mass, inferred-probe, or gamma2 requirement at this stage.
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

    return df.Filter("W > 2.0", "W_gt_2_GeV")


def draw_canvas(dfs, output_file):
    """Draw the four observables through the requested cumulative cut sequence."""
    plots = [
        ("E_gamma1", "#gamma1 energy;E_{#gamma1} (GeV);Unit-normalized entries", 180, 1.0, 10.0),
        ("Emiss_epg", "Missing energy e'p'#gamma1;E_{miss}(e'p'#gamma1) (GeV);Unit-normalized entries", 180, 1.0, 10.0),
        ("Mx2_ep", "Missing mass squared e'p';M^{2}_{X}(e'p') (GeV^{2});Unit-normalized entries", 200, -0.5, 2.0),
        ("Mx2_epg_raw", "Missing mass squared e'p'#gamma1;M^{2}_{X}(e'p'#gamma1) (GeV^{2});Unit-normalized entries", 180, -0.1, 0.2),
    ]

    stages = [
        ("W > 2 GeV", lambda df: df),
        ("-0.05 < M^{2}_{X}(e'p'#gamma1) < 0.05 GeV^{2}",
         lambda df: df.Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", "Mx2_epg_window")),
        ("M^{2}_{X}(e'p') < 0.25 GeV^{2}",
         lambda df: df.Filter("Mx2_ep < 0.25", "Mx2_ep_lt_0p25")),
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
    unique = str(abs(hash(output_file)))
    for icol in range(len(stages)):
        for irow, (expr, title, nbins, xmin, xmax) in enumerate(plots):
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

    canvas = ROOT.TCanvas("c_photon_efficiency_restart", "", 2100, 2100)
    canvas.Divide(3, 4, 0.002, 0.002)
    keep = [canvas] + list(actions)

    for irow, _plot in enumerate(plots):
        for icol, (stage_title, _filter) in enumerate(stages):
            pad_number = irow * 3 + icol + 1
            pad = canvas.cd(pad_number)
            pad.SetTicks(1, 1)
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.035)
            pad.SetBottomMargin(0.13)
            pad.SetTopMargin(0.12)
            if irow == 3:
                pad.SetLogy(True)

            legend = ROOT.TLegend(0.66, 0.67, 0.94, 0.88)
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
                integral = hist.Integral(1, hist.GetNbinsX())
                if integral > 0.0:
                    hist.Scale(1.0 / integral)
                ymax = max(ymax, hist.GetMaximum())
                if irow == 3:
                    for ibin in range(1, hist.GetNbinsX() + 1):
                        value = hist.GetBinContent(ibin)
                        if value > 0.0 and (positive_min is None or value < positive_min):
                            positive_min = value
                histograms.append((sample, label, hist))
                keep.append(hist)
                legend.AddEntry(hist, label, "l")

            first = True
            for _sample, _label, hist in histograms:
                if irow == 3:
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

            title = ROOT.TLatex()
            title.SetNDC(True)
            title.SetTextAlign(22)
            title.SetTextSize(0.033)
            title.DrawLatex(0.52, 0.955, stage_title)
            keep.append(title)
            legend.Draw()

    canvas.SaveAs(output_file)
    return keep

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
    print("\nRows entering the first canvas (W > 2 GeV only):")
    for sample, label in SAMPLES:
        if sample in count_handles:
            print(f"  {label:<14s} {int(count_handles[sample].GetValue()):,}")

    output_file = os.path.join(
        args.output_dir,
        f"1_photon_efficiency_epgamma_hypotheses_{args.period}.png",
    )

    keep = draw_canvas(dfs, output_file)
    _ = keep  # Keep ROOT objects alive through SaveAs().

    print(f"\nWrote: {output_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
