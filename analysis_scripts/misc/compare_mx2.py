#!/usr/bin/env python3

import sys
import os
import ROOT

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

    # Open ROOT files
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
        print("Error: could not load PhysicsEvents tree from one of the files.")
        sys.exit(1)
    #endif

    ROOT.gStyle.SetOptStat(0)

    # Histogram parameters
    nBins = 75
    xMin = 0.0
    xMax = 2.0

    h1 = ROOT.TH1D("h1", "", nBins, xMin, xMax)
    h2 = ROOT.TH1D("h2", "", nBins, xMin, xMax)
    h3 = ROOT.TH1D("h3", "", nBins, xMin, xMax)

    # Fill histograms
    tree1.Draw("Mx2>>h1", "", "goff")
    tree2.Draw("Mx2>>h2", "", "goff")
    tree3.Draw("Mx2>>h3", "", "goff")

    # Normalize to integral
    if h1.Integral() > 0:
        h1.Scale(1.0 / h1.Integral())
    #endif

    if h2.Integral() > 0:
        h2.Scale(1.0 / h2.Integral())
    #endif

    if h3.Integral() > 0:
        h3.Scale(1.0 / h3.Integral())
    #endif

    # Styling (match your aesthetic)
    h1.SetLineColor(ROOT.kBlue)
    h1.SetLineWidth(3)

    h2.SetLineColor(ROOT.kOrange+1)
    h2.SetLineWidth(3)

    h3.SetLineColor(ROOT.kGreen+2)
    h3.SetLineWidth(3)

    h1.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
    h1.GetYaxis().SetTitle("Probability density (GeV^{-2})")

    h1.GetXaxis().CenterTitle()
    h1.GetYaxis().CenterTitle()

    h1.GetXaxis().SetTitleSize(0.045)
    h1.GetYaxis().SetTitleSize(0.045)

    h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetLabelSize(0.04)

    # Determine maximum for overlay scaling
    max_val = max(h1.GetMaximum(), h2.GetMaximum(), h3.GetMaximum())
    h1.SetMaximum(1.10 * max_val)

    # Canvas
    canvas = ROOT.TCanvas("canvas", "Mx2 Comparison", 900, 650)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.08)

    h1.Draw("hist")
    h2.Draw("hist same")
    h3.Draw("hist same")

    # Legend (boxed, top-left)
    legend = ROOT.TLegend(0.15, 0.72, 0.48, 0.88)
    legend.SetBorderSize(1)
    legend.SetFillStyle(1001)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetTextSize(0.04)

    legend.AddEntry(h1, label1, "l")
    legend.AddEntry(h2, label2, "l")
    legend.AddEntry(h3, label3, "l")
    legend.Draw()

    # Ensure output directory exists
    os.makedirs("output", exist_ok=True)

    canvas.SaveAs("output/fermi_motion_problem.png")

    print("Saved plot to output/fermi_motion_problem.png")

    file1.Close()
    file2.Close()
    file3.Close()

#enddef

if __name__ == "__main__":
    main()
#endif