import ROOT
import os

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# input file
# ------------------------------------------------

file_rga = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_fa18_inb_epi+.root"

# ------------------------------------------------
# open file
# ------------------------------------------------

f_rga = ROOT.TFile.Open(file_rga)
if not f_rga or f_rga.IsZombie():
    raise RuntimeError("Failed to open file: " + file_rga)
#endif

# ------------------------------------------------
# get tree
# ------------------------------------------------

tree_rga = f_rga.Get("PhysicsEvents")

if not tree_rga:
    raise RuntimeError("PhysicsEvents tree missing in RGA file")
#endif

# ------------------------------------------------
# histogram
# ------------------------------------------------

nbins = 200
xmin  = 0.4
xmax  = 1.4

h_rga = ROOT.TH1F("h_rga",";Mx2 (GeV^{2});Normalized Counts",nbins,xmin,xmax)

# ------------------------------------------------
# fill histogram
# ------------------------------------------------

tree_rga.Draw("Mx2>>h_rga","","goff")

# ------------------------------------------------
# normalize to integral
# ------------------------------------------------

integral = h_rga.Integral()

if integral > 0:
    h_rga.Scale(1.0 / integral)
#endif

# ------------------------------------------------
# style
# ------------------------------------------------

h_rga.SetLineColor(ROOT.kBlack)
h_rga.SetLineWidth(2)

# ------------------------------------------------
# set y-axis max
# ------------------------------------------------

h_rga.SetMaximum(1.2 * h_rga.GetMaximum())

# ------------------------------------------------
# canvas
# ------------------------------------------------

c = ROOT.TCanvas("c","Mx2 peak",800,600)

h_rga.Draw("hist")

# ------------------------------------------------
# output
# ------------------------------------------------

os.makedirs("output", exist_ok=True)

c.SaveAs("output/epi+_Mx2_peak.png")