#!/usr/bin/env python3

import sys
import os
import ROOT

def normalize_histogram(h):
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    #endif
#enddef


def main():

    if len(sys.argv) != 7:
        print("Usage: python3 compare_mx2.py file1.root file2.root file3.root label1 label2 label3")
        sys.exit(1)
    #endif

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    file3_path = sys.argv[3]

    label1 = sys.argv[4]
    label2 = sys.argv[5]
    label3 = sys.argv[6]

    file1 = ROOT.TFile.Open(file1_path)
    file2 = ROOT.TFile.Open(file2_path)
    file3 = ROOT.TFile.Open(file3_path)

    if not file1 or file1.IsZombie():
        print(f"Error: cannot open {file1_path}")
        sys.exit(1)
    #endif
    if not file2 or file2.IsZombie():
        print(f"Error: cannot open {file2_path}")
        sys.exit(1)
    #endif
    if not file3 or file3.IsZombie():
        print(f"Error: cannot open {file3_path}")
        sys.exit(1)
    #endif

    tree1 = file1.Get("PhysicsEvents")
    tree2 = file2.Get("PhysicsEvents")
    tree3 = file3.Get("PhysicsEvents")

    if not tree1 or not tree2 or not tree3:
        print("Error: could not load PhysicsEvents tree.")
        sys.exit(1)
    #endif

    ROOT.gStyle.SetOptStat(0)

    # -----------------------------
    # Global Mx2 comparison canvas
    # -----------------------------

    nBins = 150
    xMin = 0.0
    xMax = 2.0

    h1 = ROOT.TH1D("h1", "", nBins, xMin, xMax)
    h2 = ROOT.TH1D("h2", "", nBins, xMin, xMax)
    h3 = ROOT.TH1D("h3", "", nBins, xMin, xMax)

    tree1.Draw("Mx2>>h1", "", "goff")
    tree2.Draw("Mx2>>h2", "", "goff")
    tree3.Draw("Mx2>>h3", "", "goff")

    normalize_histogram(h1)
    normalize_histogram(h2)
    normalize_histogram(h3)

    h1.SetLineColor(ROOT.kBlue)
    h2.SetLineColor(ROOT.kOrange+1)
    h3.SetLineColor(ROOT.kGreen+2)

    h1.SetLineWidth(3)
    h2.SetLineWidth(3)
    h3.SetLineWidth(3)

    h1.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
    h1.GetYaxis().SetTitle("Probability density (GeV^{-2})")

    max_val = max(h1.GetMaximum(), h2.GetMaximum(), h3.GetMaximum())
    h1.SetMaximum(1.10 * max_val)

    canvas = ROOT.TCanvas("canvas", "Mx2 Comparison", 900, 650)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.08)

    h1.Draw("hist")
    h2.Draw("hist same")
    h3.Draw("hist same")

    legend = ROOT.TLegend(0.15, 0.72, 0.48, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.04)

    legend.AddEntry(h1, label1, "l")
    legend.AddEntry(h2, label2, "l")
    legend.AddEntry(h3, label3, "l")
    legend.Draw()

    os.makedirs("output", exist_ok=True)

    canvas.SaveAs("output/fermi_motion_problem.png")

    # -------------------------------------
    # 1x3 canvas binned in x (x_B)
    # -------------------------------------

    x_bins = [
        (0.05, 0.20),
        (0.20, 0.40),
        (0.40, 0.60)
    ]

    canvas_binned = ROOT.TCanvas("canvas_binned", "Mx2 xB-binned", 1400, 450)
    canvas_binned.Divide(3,1)

    for i, (xmin_bin, xmax_bin) in enumerate(x_bins):

        canvas_binned.cd(i+1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)

        h1b = ROOT.TH1D(f"h1b_{i}", "", nBins, xMin, xMax)
        h2b = ROOT.TH1D(f"h2b_{i}", "", nBins, xMin, xMax)
        h3b = ROOT.TH1D(f"h3b_{i}", "", nBins, xMin, xMax)

        cut = f"x > {xmin_bin} && x < {xmax_bin}"

        tree1.Draw(f"Mx2>>h1b_{i}", cut, "goff")
        tree2.Draw(f"Mx2>>h2b_{i}", cut, "goff")
        tree3.Draw(f"Mx2>>h3b_{i}", cut, "goff")

        normalize_histogram(h1b)
        normalize_histogram(h2b)
        normalize_histogram(h3b)

        h1b.SetLineColor(ROOT.kBlue)
        h2b.SetLineColor(ROOT.kOrange+1)
        h3b.SetLineColor(ROOT.kGreen+2)

        h1b.SetLineWidth(3)
        h2b.SetLineWidth(3)
        h3b.SetLineWidth(3)

        h1b.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
        h1b.GetYaxis().SetTitle("Probability density (GeV^{-2})")

        max_val_bin = max(h1b.GetMaximum(), h2b.GetMaximum(), h3b.GetMaximum())
        h1b.SetMaximum(1.10 * max_val_bin)

        h1b.Draw("hist")
        h2b.Draw("hist same")
        h3b.Draw("hist same")

        # x_B label annotation
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.DrawLatex(0.18, 0.85, f"{xmin_bin:.2f} < x_{{B}} < {xmax_bin:.2f}")

    #endfor

    canvas_binned.SaveAs("output/fermi_motion_problem_xB_binned.png")

    print("Saved:")
    print("  output/fermi_motion_problem.png")
    print("  output/fermi_motion_problem_xB_binned.png")

    file1.Close()
    file2.Close()
    file3.Close()

#enddef


if __name__ == "__main__":
    main()
#endif