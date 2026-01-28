#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mx2_mc_mx2_fraction_table.py

Read a ROOT file at runtime, load the PhysicsEvents TTree, select events with:
  -inf <= Mx2   <= 1.07        (for the TABLE only)
  0.10 <= x     <= 0.60
  -1.0 <= tprime <= 0.0

Then compute the percentage distribution of mc_Mx2 values in the bins:
  [-inf, 1.07]
  [1.07, 1.22]
  [1.22, 1.37]
  [1.37, 1.52]
  [1.52, 1.77]
  [1.77, inf]

Also produces an output plot comparing Mx2 and mc_Mx2 histograms with ONLY the
(x, tprime) cuts applied (no Mx2 window), in the range [-1, 4], and saves it to:
  output/Mx2_comparison.png

Usage (tcsh):
  python3 mx2_mc_mx2_fraction_table.py /path/to/file.root
"""

import sys
import os
import math
import ROOT


def die(msg):
    raise SystemExit(f"FATAL: {msg}")


def format_pct(x):
    return f"{100.0 * x:8.3f}%"


def ensure_outdir(path):
    d = os.path.dirname(path)
    if d and (not os.path.isdir(d)):
        os.makedirs(d, exist_ok=True)
    #endif


def bin_label(lo, hi):
    if math.isinf(lo) and (lo < 0.0) and (not math.isinf(hi)):
        return f"[-inf, {hi:.2f}]"
    #endif
    if (not math.isinf(lo)) and math.isinf(hi) and (hi > 0.0):
        return f"[{lo:.2f}, inf]"
    #endif
    if (not math.isinf(lo)) and (not math.isinf(hi)):
        return f"[{lo:.2f}, {hi:.2f}]"
    #endif
    return f"[{lo}, {hi}]"


def in_bin(val, lo, hi, is_last_bin):
    # Half-open [lo, hi) except last bin includes upper edge.
    if not is_last_bin:
        return (val >= lo) and (val < hi)
    #endif
    return (val >= lo) and (val <= hi)
    #endif


def main():
    if len(sys.argv) != 2:
        die("Expected exactly 1 argument: /path/to/input.root")

    inpath = sys.argv[1]

    ROOT.gROOT.SetBatch(True)

    f = ROOT.TFile.Open(inpath, "READ")
    if not f or f.IsZombie():
        die(f"Could not open input ROOT file: {inpath}")

    t = f.Get("PhysicsEvents")
    if not t:
        die("Tree 'PhysicsEvents' not found in file.")

    # Fail-fast branch checks.
    required = ["Mx2", "mc_Mx2", "x", "tprime"]
    missing = []
    for br in required:
        if not (t.GetBranch(br) or t.GetLeaf(br)):
            missing.append(br)
        #endif
    #endfor
    if len(missing) > 0:
        die(f"Missing required branch/leaf(s) in PhysicsEvents: {', '.join(missing)}")

    # --------------------------
    # Selections for TABLE
    # --------------------------
    # "Mx2 between -inf and 1.07" means: Mx2 <= 1.07 (no lower bound).
    mx2_sel_max = 1.07

    x_sel_min = 0.10
    x_sel_max = 0.60

    tprime_sel_min = -1.0
    tprime_sel_max = 0.0

    # --------------------------
    # mc_Mx2 bins (with infinities)
    # --------------------------
    bins = [
        (-math.inf, 1.07),
        (1.07, 1.22),
        (1.22, 1.37),
        (1.37, 1.52),
        (1.52, 1.77),
        (1.77, math.inf),
    ]
    counts = [0] * len(bins)

    total_selected = 0
    total_outside_mc_range = 0  # should remain 0 with infinite endpoints

    n_entries = int(t.GetEntries())
    for i in range(n_entries):
        t.GetEntry(i)

        mx2 = float(getattr(t, "Mx2"))
        if mx2 > mx2_sel_max:
            continue
        #endif

        x = float(getattr(t, "x"))
        if x < x_sel_min or x > x_sel_max:
            continue
        #endif

        tp = float(getattr(t, "tprime"))
        if tp < tprime_sel_min or tp > tprime_sel_max:
            continue
        #endif

        mc_mx2 = float(getattr(t, "mc_Mx2"))
        total_selected += 1

        placed = False
        for k in range(len(bins)):
            lo, hi = bins[k]
            is_last = (k == (len(bins) - 1))
            if in_bin(mc_mx2, lo, hi, is_last):
                counts[k] += 1
                placed = True
                break
            #endif
        #endfor

        if not placed:
            total_outside_mc_range += 1
        #endif
    #endfor

    if total_selected == 0:
        die(
            "No entries passed the TABLE selection: "
            f"Mx2 <= {mx2_sel_max:.2f}, "
            f"{x_sel_min:.2f} <= x <= {x_sel_max:.2f}, "
            f"{tprime_sel_min:.2f} <= tprime <= {tprime_sel_max:.2f}."
        )

    # --------------------------
    # Print table
    # --------------------------
    print("")
    print("------------------------------------------------------------")
    print(f"Input file: {inpath}")
    print("Tree: PhysicsEvents")
    print("TABLE selection:")
    print(f"  Mx2 <= {mx2_sel_max:.2f}")
    print(f"  {x_sel_min:.2f} <= x      <= {x_sel_max:.2f}")
    print(f"  {tprime_sel_min:.2f} <= tprime <= {tprime_sel_max:.2f}")
    print(f"Selected entries (for table): {total_selected}")
    print("------------------------------------------------------------")
    print("")
    print(f"{'mc_Mx2 bin':>20}   {'count':>12}   {'percent':>12}")
    print(f"{'-'*20}   {'-'*12}   {'-'*12}")

    for k in range(len(bins)):
        lo, hi = bins[k]
        label = bin_label(lo, hi)
        frac = counts[k] / float(total_selected)
        print(f"{label:>20}   {counts[k]:12d}   {format_pct(frac):>12}")
    #endfor

    if total_outside_mc_range > 0:
        print("")
        print(
            f"WARNING: {total_outside_mc_range} selected events were not assigned to a bin. "
            "This should not happen with infinite endpoints; check for NaNs in mc_Mx2."
        )
    #endif

    print("")

    # --------------------------
    # Make comparison plot
    #   - apply ONLY x and tprime cuts (no Mx2 window)
    #   - histogram range [-1, 4]
    #   - remove stat box
    #   - legend for both
    # --------------------------
    out_png = "output/Mx2_comparison.png"
    ensure_outdir(out_png)

    ROOT.gStyle.SetOptStat(0)

    nbins = 200
    h_mx2 = ROOT.TH1F("h_mx2", "Mx2 Comparison;M_{x}^{2} (GeV^{2});Counts", nbins, -1.0, 4.0)
    h_mc  = ROOT.TH1F("h_mc_mx2", "Mx2 Comparison;M_{x}^{2} (GeV^{2});Counts", nbins, -1.0, 4.0)

    h_mx2.SetLineWidth(2)
    h_mc.SetLineWidth(2)
    h_mc.SetLineColor(ROOT.kRed)

    n_entries = int(t.GetEntries())
    for i in range(n_entries):
        t.GetEntry(i)

        x = float(getattr(t, "x"))
        if x < x_sel_min or x > x_sel_max:
            continue
        #endif

        tp = float(getattr(t, "tprime"))
        if tp < tprime_sel_min or tp > tprime_sel_max:
            continue
        #endif

        mx2 = float(getattr(t, "Mx2"))
        mc_mx2 = float(getattr(t, "mc_Mx2"))

        h_mx2.Fill(mx2)
        h_mc.Fill(mc_mx2)
    #endfor

    if h_mx2.GetEntries() <= 0:
        die("Histogram for Mx2 has zero entries after (x, tprime) cuts.")
    #endif
    if h_mc.GetEntries() <= 0:
        die("Histogram for mc_Mx2 has zero entries after (x, tprime) cuts.")
    #endif

    c = ROOT.TCanvas("c_mx2_comp", "Mx2 comparison", 900, 700)
    c.SetLeftMargin(0.13)

    ymax = max(h_mx2.GetMaximum(), h_mc.GetMaximum())
    h_mx2.SetMaximum(1.15 * ymax)

    h_mx2.Draw("HIST")
    h_mc.Draw("HIST SAME")

    leg = ROOT.TLegend(0.62, 0.72, 0.89, 0.89)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(1)
    leg.AddEntry(h_mx2, "Mx2", "l")
    leg.AddEntry(h_mc, "mc_Mx2", "l")
    leg.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC(True)
    latex.SetTextSize(0.032)
    latex.DrawLatex(
        0.14,
        0.92,
        f"{x_sel_min:.2f} #leq x #leq {x_sel_max:.2f},  {tprime_sel_min:.1f} #leq t' #leq {tprime_sel_max:.1f}",
    )

    c.SaveAs(out_png)

    print(f"Saved plot: {out_png}")
    print("Done.")


if __name__ == "__main__":
    main()
#endif