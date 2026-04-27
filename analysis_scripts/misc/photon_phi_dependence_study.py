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

sector_hists = []
sector_colors = [
    ROOT.kBlack,
    ROOT.kRed + 1,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kOrange + 7,
]

for sector_index in range(6):
    sector = sector_index + 1

    hist = ROOT.TH1D(
        "hist_p_theta_sector_{}".format(sector),
        ";#theta_{#gamma} (deg);Counts",
        140, 0.0, 35.0
    )

    hist.SetLineColor(sector_colors[sector_index])
    hist.SetMarkerColor(sector_colors[sector_index])
    hist.SetLineWidth(2)

    sector_hists.append(hist)
#endfor

def wrap_phi_deg(phi_deg):
    while phi_deg < 0.0:
        phi_deg += 360.0
    #endwhile

    while phi_deg >= 360.0:
        phi_deg -= 360.0
    #endwhile

    return phi_deg
#enddef

def get_clas12_sector_from_phi_deg(phi_deg):
    phi = wrap_phi_deg(phi_deg)

    if (phi >= 330.0 and phi < 360.0) or (phi >= 0.0 and phi < 30.0):
        return 1
    #endif

    if phi >= 30.0 and phi < 90.0:
        return 2
    #endif

    if phi >= 90.0 and phi < 150.0:
        return 3
    #endif

    if phi >= 150.0 and phi < 210.0:
        return 4
    #endif

    if phi >= 210.0 and phi < 270.0:
        return 5
    #endif

    if phi >= 270.0 and phi < 330.0:
        return 6
    #endif

    return 0
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

    sector = get_clas12_sector_from_phi_deg(p_phi_deg)
    if sector >= 1 and sector <= 6:
        sector_hists[sector - 1].Fill(p_theta_deg)
    #endif
#endfor

canvas = ROOT.TCanvas("canvas", "canvas", 2100, 600)
canvas.Divide(3, 1)

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

canvas.cd(3)
ROOT.gPad.SetLeftMargin(0.14)
ROOT.gPad.SetRightMargin(0.05)
ROOT.gPad.SetBottomMargin(0.13)

max_sector_count = 0.0
for hist in sector_hists:
    if hist.GetMaximum() > max_sector_count:
        max_sector_count = hist.GetMaximum()
    #endif
#endfor

if max_sector_count <= 0.0:
    max_sector_count = 1.0
#endif

sector_hists[0].SetMaximum(1.20 * max_sector_count)
sector_hists[0].Draw("HIST")

legend = ROOT.TLegend(0.58, 0.60, 0.92, 0.88)
legend.SetBorderSize(1)
legend.SetFillStyle(1001)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextSize(0.035)

legend.AddEntry(sector_hists[0], "Sector 1", "l")

for sector_index in range(1, 6):
    sector_hists[sector_index].Draw("HIST SAME")
    legend.AddEntry(sector_hists[sector_index], "Sector {}".format(sector_index + 1), "l")
#endfor

legend.Draw()

canvas.SaveAs(output_file_name)

root_file.Close()

print("Saved: {}".format(output_file_name))