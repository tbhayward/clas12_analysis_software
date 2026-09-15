#!/usr/bin/env python3
"""
CLAS12 photon-efficiency analysis -- restart from the raw e'p'gamma1 hypotheses.

Current step
------------
Use every row in the PhotonEfficiency hypothesis tree, requiring W > 2 GeV and angle(e',gamma1) > 8 deg.
Each row is one reconstructed e'p'gamma1 hypothesis produced by the skim.

Plot each of the following observables on its own 2x2 canvas:
  1) E_gamma1                     0.4 to 9 GeV
  2) missing energy e'p'gamma1    0 to 9 GeV
  3) Mx2(e'p')                   -0.5 to 1.0 GeV^2
  4) Mx2(e'p'gamma1)             -0.1 to 0.15 GeV^2
  5) Mx2(e'gamma1)                -20 to 20 GeV^2

Pads are cumulative selections, ordered top-left, top-right, bottom-left, bottom-right:
  1) W > 2 GeV and angle(e',gamma1) > 8 deg
  2) additionally Mx2(e'p') < 0.18 GeV^2
  3) additionally -0.05 < Mx2(e'p'gamma1) < 0.05 GeV^2
  4) additionally Mx2(e'gamma1) > 0 GeV^2

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
        ("M^{2}_{X}(e'p') < 0.18 GeV^{2}",
         lambda df: df.Filter("Mx2_ep < 0.18", "Mx2_ep_lt_0p18")),
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
    """Plot the E_gamma1 > 5 GeV normalization control region after exclusivity cuts."""
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
            .Filter("Mx2_ep < 0.18", "norm_Mx2_ep_lt_0p18")
            .Filter("Mx2_epg_raw > -0.05 && Mx2_epg_raw < 0.05", "norm_Mx2_epg_window")
            .Filter("E_gamma1 > 5.0", "norm_Egamma1_gt_5"))

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

    output_file = os.path.join(output_dir, f"1_{period}_Egamma1_gt_5_GeV.png")
    canvas.SaveAs(output_file)
    return keep, output_file


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
              .Filter("Mx2_ep < 0.18", f"scan_{sample}_Mx2_ep_lt_0p18")
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
    print("  Denominator = all rows after Mx2(ep) < 0.18 and -0.05 < Mx2(epgamma1) < 0.05")
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

    print_egamma1_survival_scan(dfs)

    normalization_dir = os.path.join(args.output_dir, "normalization")
    os.makedirs(normalization_dir, exist_ok=True)
    norm_keep, norm_output = draw_normalization_step1(dfs, normalization_dir, args.period)
    keep.extend(norm_keep)
    _ = keep  # Keep ROOT objects alive through SaveAs().

    for output_file in written:
        print(f"\nWrote: {output_file}")
    print(f"\nWrote: {norm_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
