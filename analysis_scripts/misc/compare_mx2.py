#!/usr/bin/env python3

import sys
import os
import ROOT
import math

def normalize_histogram(h):
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    #endif
#enddef


def style_histograms(h1, h2, h3):
    h1.SetLineColor(ROOT.kBlue)
    h2.SetLineColor(ROOT.kOrange+1)
    h3.SetLineColor(ROOT.kGreen+2)

    h1.SetLineWidth(3)
    h2.SetLineWidth(3)
    h3.SetLineWidth(3)
#enddef


def draw_global(tree1, tree2, tree3, var, nBins, xmin, xmax,
                xTitle, yTitle, outfile,
                label1, label2, label3):

    h1 = ROOT.TH1D(f"h1_{var}", "", nBins, xmin, xmax)
    h2 = ROOT.TH1D(f"h2_{var}", "", nBins, xmin, xmax)
    h3 = ROOT.TH1D(f"h3_{var}", "", nBins, xmin, xmax)

    tree1.Draw(f"{var}>>h1_{var}", "", "goff")
    tree2.Draw(f"{var}>>h2_{var}", "", "goff")
    tree3.Draw(f"{var}>>h3_{var}", "", "goff")

    normalize_histogram(h1)
    normalize_histogram(h2)
    normalize_histogram(h3)

    style_histograms(h1, h2, h3)

    h1.GetXaxis().SetTitle(xTitle)
    h1.GetYaxis().SetTitle(yTitle)

    max_val = max(h1.GetMaximum(), h2.GetMaximum(), h3.GetMaximum())
    h1.SetMaximum(1.10 * max_val)

    canvas = ROOT.TCanvas(f"c_{var}", var, 900, 650)
    canvas.SetLeftMargin(0.13)
    canvas.SetBottomMargin(0.12)

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

    canvas.SaveAs(outfile)
#enddef


def draw_binned(tree1, tree2, tree3, var, nBins, xmin, xmax,
                xTitle, yTitle, outfile,
                label1, label2, label3):

    x_bins = [(0.05,0.20),(0.20,0.40),(0.40,0.60)]

    canvas = ROOT.TCanvas(f"c_{var}_binned", var+"_binned", 1400, 450)
    canvas.Divide(3,1)

    for i,(xmin_bin,xmax_bin) in enumerate(x_bins):

        canvas.cd(i+1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)

        h1 = ROOT.TH1D(f"h1_{var}_{i}", "", nBins, xmin, xmax)
        h2 = ROOT.TH1D(f"h2_{var}_{i}", "", nBins, xmin, xmax)
        h3 = ROOT.TH1D(f"h3_{var}_{i}", "", nBins, xmin, xmax)

        cut = f"x > {xmin_bin} && x < {xmax_bin}"

        tree1.Draw(f"{var}>>h1_{var}_{i}", cut, "goff")
        tree2.Draw(f"{var}>>h2_{var}_{i}", cut, "goff")
        tree3.Draw(f"{var}>>h3_{var}_{i}", cut, "goff")

        normalize_histogram(h1)
        normalize_histogram(h2)
        normalize_histogram(h3)

        style_histograms(h1, h2, h3)

        h1.GetXaxis().SetTitle(xTitle)
        h1.GetYaxis().SetTitle(yTitle)

        max_val = max(h1.GetMaximum(), h2.GetMaximum(), h3.GetMaximum())
        h1.SetMaximum(1.10 * max_val)

        h1.Draw("hist")
        h2.Draw("hist same")
        h3.Draw("hist same")

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.DrawLatex(0.18,0.85,f"{xmin_bin:.2f} < x_{{B}} < {xmax_bin:.2f}")

    #endfor

    canvas.SaveAs(outfile)
#enddef


def main():

    if len(sys.argv) != 7:
        print("Usage: python3 compare_mx2.py file1.root file2.root file3.root label1 label2 label3")
        sys.exit(1)
    #endif

    file1 = ROOT.TFile.Open(sys.argv[1])
    file2 = ROOT.TFile.Open(sys.argv[2])
    file3 = ROOT.TFile.Open(sys.argv[3])

    label1 = sys.argv[4]
    label2 = sys.argv[5]
    label3 = sys.argv[6]

    tree1 = file1.Get("PhysicsEvents")
    tree2 = file2.Get("PhysicsEvents")
    tree3 = file3.Get("PhysicsEvents")

    ROOT.gStyle.SetOptStat(0)
    os.makedirs("output", exist_ok=True)

    # ---------------- Existing plots ----------------

    draw_global(tree1,tree2,tree3,"Mx2",150,0.0,2.0,
                "M_{x}^{2} (GeV^{2})","Probability density (GeV^{-2})",
                "output/fermi_motion_problem.png",
                label1,label2,label3)

    draw_binned(tree1,tree2,tree3,"Mx2",50,0.0,2.0,
                "M_{x}^{2} (GeV^{2})","Probability density (GeV^{-2})",
                "output/fermi_motion_problem_xB_binned.png",
                label1,label2,label3)

    draw_global(tree1,tree2,tree3,"phi",100,0.0,2.0*math.pi,
                "#phi","Probability density",
                "output/phi_comparison.png",
                label1,label2,label3)

    draw_binned(tree1,tree2,tree3,"phi",100,0.0,2.0*math.pi,
                "#phi","Probability density",
                "output/phi_xB_binned.png",
                label1,label2,label3)

    # ---------------- NEW W plots ----------------

    draw_global(tree1,tree2,tree3,"W",120,2.0,6.0,
                "W (GeV)","Probability density (GeV^{-1})",
                "output/W_comparison.png",
                label1,label2,label3)

    draw_binned(tree1,tree2,tree3,"W",80,2.0,6.0,
                "W (GeV)","Probability density (GeV^{-1})",
                "output/W_xB_binned.png",
                label1,label2,label3)

    # ---------------- NEW xF plots ----------------

    draw_global(tree1,tree2,tree3,"xF",120,-1.0,2.0,
                "x_{F}","Probability density",
                "output/xF_comparison.png",
                label1,label2,label3)

    draw_binned(tree1,tree2,tree3,"xF",80,-1.0,2.0,
                "x_{F}","Probability density",
                "output/xF_xB_binned.png",
                label1,label2,label3)

    print("All plots saved in output/")

    file1.Close()
    file2.Close()
    file3.Close()

#enddef


if __name__ == "__main__":
    main()
#endif