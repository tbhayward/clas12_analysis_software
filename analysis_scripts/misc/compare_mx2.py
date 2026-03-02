#!/usr/bin/env python3

import sys
import os
import ROOT

def main():

    if len(sys.argv) != 5:
        print("Usage: python3 compare_mx2.py file1.root file2.root label1 label2")
        sys.exit(1)
    #endif

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    label1 = sys.argv[3]
    label2 = sys.argv[4]

    # Open ROOT files
    file1 = ROOT.TFile.Open(file1_path)
    if not file1 or file1.IsZombie():
        print(f"Error: cannot open {file1_path}")
        sys.exit(1)
    #endif

    file2 = ROOT.TFile.Open(file2_path)
    if not file2 or file2.IsZombie():
        print(f"Error: cannot open {file2_path}")
        sys.exit(1)
    #endif

    tree1 = file1.Get("PhysicsEvents")
    tree2 = file2.Get("PhysicsEvents")

    if not tree1 or not tree2:
        print("Error: could not load PhysicsEvents tree from one of the files.")
        sys.exit(1)
    #endif

    ROOT.gStyle.SetOptStat(0)

    # Histogram parameters
    nBins = 150
    xMin = 0.0
    xMax = 1.5

    h1 = ROOT.TH1D("h1", "", nBins, xMin, xMax)
    h2 = ROOT.TH1D("h2", "", nBins, xMin, xMax)

    # Fill histograms
    tree1.Draw("Mx2>>h1", "", "goff")
    tree2.Draw("Mx2>>h2", "", "goff")

    # Normalize to integral
    if h1.Integral() > 0:
        h1.Scale(1.0 / h1.Integral())
    #endif

    if h2.Integral() > 0:
        h2.Scale(1.0 / h2.Integral())
    #endif

    # Styling
    h1.SetLineColor(ROOT.kBlack)
    h1.SetLineWidth(2)

    h2.SetLineColor(ROOT.kRed)
    h2.SetLineWidth(2)

    h1.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
    h1.GetYaxis().SetTitle("Normalized counts")

    h1.GetXaxis().CenterTitle()
    h1.GetYaxis().CenterTitle()

    # Determine maximum for proper overlay scaling
    max_val = max(h1.GetMaximum(), h2.GetMaximum())
    h1.SetMaximum(1.15 * max_val)

    # Canvas
    canvas = ROOT.TCanvas("canvas", "Mx2 Comparison", 800, 600)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)

    h1.Draw("hist")
    h2.Draw("hist same")

    # Legend
    legend = ROOT.TLegend(0.60, 0.75, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.AddEntry(h1, label1, "l")
    legend.AddEntry(h2, label2, "l")
    legend.Draw()

    # Ensure output directory exists
    os.makedirs("output", exist_ok=True)

    canvas.SaveAs("output/fermi_motion_problem.png")

    print("Saved plot to output/fermi_motion_problem.png")

    file1.Close()
    file2.Close()

#enddef

if __name__ == "__main__":
    main()
#endif