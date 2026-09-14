#!/usr/bin/env python3
"""
Quick QA plots for the CLAS12 photon-efficiency skim.

- Automatically discovers completed ROOT files for data, CLASDIS, AAOgen, DVCSgen.
- Skips empty/corrupt files and files that are not skim_version == 4.
- Uses every completed file currently available.
- Makes one 2x3 canvas with overlaid, unit-area distributions.
- Does NOT apply the final missing-mass cuts, so the distributions remain visible.
- Applies only the standard reconstructed-object QA cuts by default:
      p_pass_standard == 1
      tag_pass_beta == 1
      tag_pass_fiducial == 1

Usage:
    python3 plot_photon_efficiency_quick_qa.py

Optional:
    python3 plot_photon_efficiency_quick_qa.py --period fa18_inb
    python3 plot_photon_efficiency_quick_qa.py --no-object-cuts
    python3 plot_photon_efficiency_quick_qa.py --output my_qa.png
"""

import argparse
import glob
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

BASE = "/work/clas12/thayward/photon_efficiency/ROOT_trees"
TREE = "PhotonEfficiency"
EXPECTED_SKIM_VERSION = 4

SAMPLES = [
    ("data",    "Data"),
    ("clasdis", "CLASDIS (exclusive ep#pi^{0} removed)"),
    ("aaogen",  "AAOgen"),
    ("dvcsgen", "DVCSgen"),
]

# ROOT color palette; deliberately fixed so the same sample always has the same color.
COLORS = {
    "data": ROOT.kBlack,
    "clasdis": ROOT.kBlue + 1,
    "aaogen": ROOT.kRed + 1,
    "dvcsgen": ROOT.kGreen + 2,
}

# (expression, title, nbins, xmin, xmax, vertical cut lines)
PLOTS = [
    (
        "Mx2_ep",
        "M^{2}_{X}(e'p') wide;M^{2}_{X}(e'p') [GeV^{2}];Unit-normalized entries",
        140, -1.0, 2.0, [-0.10, 0.15],
    ),
    (
        "Mx2_ep",
        "M^{2}_{X}(e'p') near photon peak;M^{2}_{X}(e'p') [GeV^{2}];Unit-normalized entries",
        140, -0.30, 0.45, [-0.10, 0.15],
    ),
    (
        "Mx2_epg_raw",
        "Raw M^{2}_{X}(e'p'#gamma_{tag});M^{2}_{X}(e'p'#gamma_{tag}) [GeV^{2}];Unit-normalized entries",
        140, -0.25, 0.25, [-0.10, 0.10],
    ),
    (
        "Mx2_epg_corr",
        "Analysis M^{2}_{X}(e'p'#gamma_{tag});M^{2}_{X}(e'p'#gamma_{tag}) [GeV^{2}];Unit-normalized entries",
        140, -0.25, 0.25, [-0.10, 0.10],
    ),
    (
        "W",
        "DIS invariant mass W;W [GeV];Unit-normalized entries",
        120, 1.5, 4.5, [2.0],
    ),
    (
        "minus_t",
        "-t distribution;-t [GeV^{2}];Unit-normalized entries",
        120, 0.0, 2.5, [],
    ),
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--period", default="fa18_inb")
    p.add_argument(
        "--output",
        default=None,
        help="Output PNG. Default: photon_efficiency_quick_qa_<period>.png",
    )
    p.add_argument(
        "--no-object-cuts",
        action="store_true",
        help="Do not require p_pass_standard, tag_pass_beta, tag_pass_fiducial.",
    )
    return p.parse_args()


def file_is_usable(path):
    """Return (usable, entries, reason)."""
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        return False, 0, "cannot open"

    t = f.Get(TREE)
    if not t:
        f.Close()
        return False, 0, f"missing {TREE}"

    n = int(t.GetEntries())
    if n == 0:
        f.Close()
        return False, 0, "empty tree"

    if not t.GetBranch("skim_version"):
        f.Close()
        return False, n, "missing skim_version"

    # Inspect the first entry only. Every file should contain a single skim version.
    t.GetEntry(0)
    try:
        version = int(getattr(t, "skim_version"))
    except Exception:
        f.Close()
        return False, n, "cannot read skim_version"

    f.Close()

    if version != EXPECTED_SKIM_VERSION:
        return False, n, f"skim_version={version}, expected {EXPECTED_SKIM_VERSION}"

    return True, n, "ok"


def discover_sample(sample, period):
    pattern = os.path.join(BASE, sample, period, "*_photon_efficiency.root")
    candidates = sorted(glob.glob(pattern))

    good = []
    total_entries = 0

    print(f"\n[{sample}] {pattern}")
    if not candidates:
        print("  no completed ROOT files found")
        return good, total_entries

    for path in candidates:
        usable, n, reason = file_is_usable(path)
        name = os.path.basename(path)
        if usable:
            good.append(path)
            total_entries += n
            print(f"  USE   {name:70s} {n:12d} entries")
        else:
            print(f"  SKIP  {name:70s} {reason}")

    print(f"  -> {len(good)} usable file(s), {total_entries} hypothesis rows")
    return good, total_entries


def make_chain(paths):
    chain = ROOT.TChain(TREE)
    for path in paths:
        chain.Add(path)
    return chain


def branch_exists(chain, name):
    branches = chain.GetListOfBranches()
    return bool(branches and branches.FindObject(name))


def main():
    args = parse_args()
    output = args.output or f"photon_efficiency_quick_qa_{args.period}.png"

    chains = {}
    sample_entries = {}

    print(f"Looking for completed skim-version-{EXPECTED_SKIM_VERSION} files...")
    for sample, label in SAMPLES:
        paths, n = discover_sample(sample, args.period)
        if paths:
            chains[sample] = make_chain(paths)
            sample_entries[sample] = n

    if not chains:
        print("\nERROR: no usable ROOT files found.", file=sys.stderr)
        return 1

    print("\nSamples that will be plotted:")
    for sample, label in SAMPLES:
        if sample in chains:
            print(f"  {label}: {len(chains[sample].GetListOfFiles())} file(s), "
                  f"{chains[sample].GetEntries()} entries")

    # These are object-quality cuts only. Deliberately do NOT impose the final
    # Mx2(ep) or Mx2(epg) cuts, because this script is intended to inspect them.
    if args.no_object_cuts:
        base_cut = "1"
        print("\nObject cuts: NONE")
    else:
        base_cut = "p_pass_standard==1 && tag_pass_beta==1 && tag_pass_fiducial==1"
        print(f"\nObject cuts: {base_cut}")

    canvas = ROOT.TCanvas("c_qa", "Photon-efficiency skim QA", 1700, 1050)
    canvas.Divide(3, 2, 0.002, 0.002)

    # Keep Python references alive until SaveAs().
    histograms = []
    legends = []
    lines = []

    for ipad, (expr, title, nbins, xmin, xmax, cut_positions) in enumerate(PLOTS, start=1):
        pad = canvas.cd(ipad)
        pad.SetTicks(1, 1)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.035)
        pad.SetBottomMargin(0.12)
        pad.SetTopMargin(0.10)

        legend = ROOT.TLegend(0.48, 0.68, 0.94, 0.89)
        legend.SetBorderSize(0)
        legend.SetFillStyle(0)
        legend.SetTextSize(0.028)
        legends.append(legend)

        panel_hists = []
        ymax = 0.0

        for sample, label in SAMPLES:
            if sample not in chains:
                continue

            chain = chains[sample]
            if not branch_exists(chain, expr):
                print(f"WARNING: {sample} does not have branch '{expr}'; skipping it in panel {ipad}")
                continue

            hname = f"h_{ipad}_{sample}"
            hist = ROOT.TH1D(hname, title, nbins, xmin, xmax)
            hist.SetDirectory(0)
            hist.SetLineColor(COLORS[sample])
            hist.SetLineWidth(2 if sample != "data" else 3)
            hist.SetStats(0)

            # Unweighted by design: this is a reconstruction/shape QA plot.
            nfilled = chain.Draw(f"{expr}>>{hname}", base_cut, "goff")
            integral = hist.Integral(1, hist.GetNbinsX())

            if integral > 0:
                hist.Scale(1.0 / integral)
            else:
                print(f"WARNING: {sample}, panel {expr}: histogram is empty after cuts")

            ymax = max(ymax, hist.GetMaximum())
            panel_hists.append((sample, label, hist, int(nfilled)))
            histograms.append(hist)

        if not panel_hists:
            frame = ROOT.TH1D(f"empty_{ipad}", title, nbins, xmin, xmax)
            frame.SetDirectory(0)
            frame.SetStats(0)
            frame.Draw()
            histograms.append(frame)
            continue

        # Give enough headroom for legend and avoid autoscale differences.
        first = True
        for sample, label, hist, nfilled in panel_hists:
            hist.SetMaximum(ymax * 1.30 if ymax > 0 else 1.0)
            hist.SetMinimum(0.0)
            hist.GetXaxis().SetTitleSize(0.042)
            hist.GetYaxis().SetTitleSize(0.040)
            hist.GetXaxis().SetLabelSize(0.034)
            hist.GetYaxis().SetLabelSize(0.034)
            hist.GetXaxis().SetTitleOffset(1.18)
            hist.GetYaxis().SetTitleOffset(1.42)

            hist.Draw("HIST" if first else "HIST SAME")
            first = False
            legend.AddEntry(hist, f"{label} (N={nfilled:,})", "l")

        # Draw nominal final-analysis cut positions only as visual guides.
        for x in cut_positions:
            line = ROOT.TLine(x, 0.0, x, ymax * 1.05)
            line.SetLineStyle(2)
            line.SetLineWidth(2)
            line.SetLineColor(ROOT.kGray + 2)
            line.Draw("SAME")
            lines.append(line)

        legend.Draw()

    canvas.cd()
    canvas.SaveAs(output)

    print(f"\nSaved: {os.path.abspath(output)}")
    print("\nNotes:")
    print("  * Histograms are unit-normalized, so this compares SHAPES, not relative rates.")
    print("  * MC generator weights are intentionally not used for this quick reconstruction QA.")
    print("  * Final missing-mass cuts are NOT applied; dashed lines show their nominal locations.")
    print("  * AAOgen/DVCSgen are included automatically as soon as completed ROOT files exist.")
    print("  * Empty QADB-rejected data files (e.g. run 5032) are skipped automatically.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
