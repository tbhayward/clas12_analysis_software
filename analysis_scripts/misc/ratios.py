#!/usr/bin/env python3
import ROOT
import math

# Turn off the default stats box
ROOT.gStyle.SetOptStat(0)

# Open the two ROOT files
file1 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/eppi+pi-/eppi+pi-_cj11.root")
file2 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/eppi+pi-/eppi+pi-_cj13.root")

# Load the PhysicsEvents trees
tree1 = file1.Get("PhysicsEvents")
tree2 = file2.Get("PhysicsEvents")

# Create histograms for e_p from 1 to 4.5 GeV (35 bins of 0.1 GeV each)
h1 = ROOT.TH1F("h1", "e_{p} distribution; e_{p} (GeV); counts", 35, 1.0, 4.5)
h2 = ROOT.TH1F("h2", "e_{p} distribution; e_{p} (GeV); counts", 35, 1.0, 4.5)

# Fill the histograms
tree2.Draw("e_p >> h2")
tree1.Draw("e_p >> h1")

# Style the histograms
h1.SetLineColor(ROOT.kBlack)
h2.SetLineColor(ROOT.kRed)
h1.SetLineWidth(2)
h2.SetLineWidth(2)

# Calculate total counts
total1 = h1.Integral()
total2 = h2.Integral()

# Create a canvas with 1 row, 2 columns
c = ROOT.TCanvas("c", "e_p comparison", 1200, 500)
c.Divide(2, 1)

# --- Left pad: the two histograms ---
c.cd(1)
h1.Draw("HIST")
h2.Draw("HIST SAME")

# Add a legend in the top-right corner with totals
leg = ROOT.TLegend(0.60, 0.60, 0.88, 0.88)
leg.AddEntry(h1, f"cj 11.2.0 (N={int(total1)})", "l")
leg.AddEntry(h2, f"cj 13.0.3 (N={int(total2)})", "l")
leg.Draw()

# --- Right pad: ratio plot with explicit Poisson errors (n2/n1) ---
c.cd(2)

# Create an empty ratio histogram with the same binning
h_ratio = ROOT.TH1F(
    "h_ratio",
    "Ratio; e_{p} (GeV); cj13.0.3 / cj11.2.0",
    35, 1.0, 4.5
)

n_bins = h1.GetNbinsX()

for i in range(1, n_bins + 1):
    n1 = h1.GetBinContent(i)
    n2 = h2.GetBinContent(i)
    if n1 > 0:
        # Poisson errors on each count
        err1 = math.sqrt(n1)
        err2 = math.sqrt(n2)
        # Ratio n2/n1 and propagated error:
        r = n2 / n1
        err = math.sqrt((err2 / n1)**2 + (n2 * err1 / (n1**2))**2)
    else:
        r = 0.0
        err = 0.0

    h_ratio.SetBinContent(i, r)
    h_ratio.SetBinError(i, err)

# Set the y-axis range for ratio
h_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
h_ratio.SetLineColor(ROOT.kBlue)
h_ratio.SetLineWidth(2)

# Draw the ratio with error bars
h_ratio.Draw("E1")

# Save the canvas to a PDF file
c.SaveAs("/u/home/thayward/ep.pdf")