#!/usr/bin/env python3

import ROOT

# Input / output
input_file = "/scratch/thayward/test.root"
tree_name  = "PhysicsEvents"
output_png = "/u/home/thayward/Mx2.png"

# Open file
f = ROOT.TFile.Open(input_file)
if not f or f.IsZombie():
    raise RuntimeError("Failed to open file: {}".format(input_file))
#endif

# Get tree
tree = f.Get(tree_name)
if not tree:
    raise RuntimeError("Tree '{}' not found in file.".format(tree_name))
#endif

# Create histogram
hist = ROOT.TH1D("h_Mx2", "Mx2;Mx2 (GeV^{2});Counts", 100, 0.0, 4.0)

# Fill histogram
for event in tree:
    hist.Fill(event.Mx2)
#endfor

# Draw and save
canvas = ROOT.TCanvas("c", "c", 800, 600)
hist.Draw()
canvas.SaveAs(output_png)

print("Saved:", output_png)

f.Close()