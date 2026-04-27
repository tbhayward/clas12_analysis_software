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
max_entries = 1000000

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

required_branches = [
    "e_phi",
    "e_theta",
    "p_phi",
    "p_theta",
    "open_angle",
]

for branch_name in required_branches:
    if not tree.GetBranch(branch_name):
        print("Error: missing required branch '{}' in tree '{}'".format(branch_name, tree_name))
        root_file.Close()
        sys.exit(1)
    #endif
#endfor

hist_e = ROOT.TH2D(
    "hist_e_phi_vs_e_theta",
    ";#phi_{e} (deg);#theta_{e} (deg)",
    180, 0.0, 360.0,
    140, 0.0, 35.0
)

hist_p = ROOT.TH2D(
    "hist_p_phi_vs_p_theta",
    ";#phi_{#gamma} (deg);#theta_{#gamma} (deg)",
    180, 0.0, 360.0,
    140, 0.0, 35.0
)

def wrap_phi_deg(phi_deg):
    while phi_deg < 0.0:
        phi_deg += 360.0
    #endwhile

    while phi_deg >= 360.0:
        phi_deg -= 360.0
    #endwhile

    return phi_deg
#enddef

n_entries = tree.GetEntries()
n_entries_to_process = min(n_entries, max_entries)

print("Tree entries: {}".format(n_entries))
print("Processing entries: {}".format(n_entries_to_process))

for i in range(n_entries_to_process):
    tree.GetEntry(i)

    if tree.open_angle <= 3.0:
        continue
    #endif

    e_phi_deg = wrap_phi_deg(tree.e_phi * 180.0 / math.pi)
    e_theta_deg = tree.e_theta * 180.0 / math.pi

    p_phi_deg = wrap_phi_deg(tree.p_phi * 180.0 / math.pi)
    p_theta_deg = tree.p_theta * 180.0 / math.pi

    hist_e.Fill(e_phi_deg, e_theta_deg)
    hist_p.Fill(p_phi_deg, p_theta_deg)
#endfor

canvas = ROOT.TCanvas("canvas", "canvas", 1400, 600)
canvas.Divide(2, 1)

canvas.cd(1)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetRightMargin(0.16)
ROOT.gPad.SetBottomMargin(0.13)
hist_e.Draw("COLZ")

canvas.cd(2)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetRightMargin(0.16)
ROOT.gPad.SetBottomMargin(0.13)
hist_p.Draw("COLZ")

canvas.SaveAs(output_file_name)

root_file.Close()

print("Saved: {}".format(output_file_name))