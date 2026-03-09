#!/usr/bin/env python3

import ROOT
import os
from array import array

# ------------------------------------------------------------------
# Input data
# ------------------------------------------------------------------

x  = [0.88, 1.15, 1.29, 1.47, 1.65]
y  = [0.602441964, 0.534944992, 0.202629081, 0.060591614, -0.061152670]
ey = [0.021986625, 0.065386550, 0.082171868, 0.136798488, 0.114795575]

n = len(x)

# Convert to ROOT-compatible arrays
x_arr  = array('d', x)
y_arr  = array('d', y)
ex_arr = array('d', [0.0]*n)
ey_arr = array('d', ey)

# ------------------------------------------------------------------
# Ensure output directory exists
# ------------------------------------------------------------------

os.makedirs("output", exist_ok=True)

# ------------------------------------------------------------------
# Create canvas
# ------------------------------------------------------------------

c = ROOT.TCanvas("c", "c", 800, 600)

# Graph with errors
g = ROOT.TGraphErrors(n, x_arr, y_arr, ex_arr, ey_arr)

g.SetMarkerStyle(20)
g.SetMarkerSize(1.2)
g.SetLineWidth(2)

# Axis labels and title (ROOT LaTeX formatting)
g.SetTitle("0.25 < x_{B} < 0.35; 0.05 < -t' < 0.25")

g.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
g.GetYaxis().SetTitle("F_{LL}/F_{UU}")

g.GetXaxis().CenterTitle()
g.GetYaxis().CenterTitle()

# Draw
g.Draw("AP")

# ------------------------------------------------------------------
# Save figure
# ------------------------------------------------------------------

c.SaveAs("output/quick_Mx2_distribution.png")