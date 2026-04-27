#!/usr/bin/env python3

import os
import sys
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

if len(sys.argv) != 2:
    print("Usage: python plot_2d_phi_vs_theta.py input.root")
    sys.exit(1)
#endif

input_file_name = sys.argv[1]
tree_name = "PhysicsEvents"
output_dir = "output/photon_phi_dependence_study"
output_file_name = os.path.join(output_dir, "2D_phi_vs_theta.pdf")

os.makedirs(output_dir, exist_ok=True)

root_file = ROOT.TFile.Open(input_file_name, "READ")
if not root_file or root_file.IsZombie():
    print("Error: could not open input file: {}".format(input_file_name))
    sys.exit(1)
#endif

tree = root_file.Get(tree_name)
if not tree:
    print("Error: could not find tree '{}' in file: {}".format(tree_name, input_file_name))
    root_file.Close()
    sys.exit(1)
#endif

required_branches = ["e_phi", "e_theta", "open_angle"]
for branch_name in required_branches:
    if not tree.GetBranch(branch_name):
        print("Error: missing required branch '{}' in tree '{}'".format(branch_name, tree_name))
        root_file.Close()
        sys.exit(1)
    #endif
#endfor

hist = ROOT.TH2D(
    "hist_e_phi_vs_e_theta",
    "Electron #phi vs Electron #theta;#phi_{e} (deg);#theta_{e} (deg)",
    180, -180.0, 180.0,
    140, 0.0, 70.0
)

n_entries = tree.GetEntries()

for i in range(n_entries):
    tree.GetEntry(i)

    if tree.open_angle <= 3.0:
        continue
    #endif

    e_phi_deg = tree.e_phi * 180.0 / math.pi
    e_theta_deg = tree.e_theta * 180.0 / math.pi

    hist.Fill(e_phi_deg, e_theta_deg)
#endfor

canvas = ROOT.TCanvas("canvas", "canvas", 900, 700)
canvas.SetLeftMargin(0.13)
canvas.SetRightMargin(0.15)
canvas.SetBottomMargin(0.12)

hist.Draw("COLZ")

canvas.SaveAs(output_file_name)

root_file.Close()

print("Saved: {}".format(output_file_name))