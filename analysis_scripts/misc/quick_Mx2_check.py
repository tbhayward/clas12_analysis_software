import ROOT
import os

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# input files
# ------------------------------------------------

file_rga = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_fa18_inb_epi+.root"
file_rgc = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rgc_su22_inb_NH3_epi+.root"

# ------------------------------------------------
# open files
# ------------------------------------------------

f_rga = ROOT.TFile.Open(file_rga)
if not f_rga or f_rga.IsZombie():
    raise RuntimeError("Failed to open file: " + file_rga)
#endif

f_rgc = ROOT.TFile.Open(file_rgc)
if not f_rgc or f_rgc.IsZombie():
    raise RuntimeError("Failed to open file: " + file_rgc)
#endif

# ------------------------------------------------
# get trees
# ------------------------------------------------

tree_rga = f_rga.Get("PhysicsEvents")
tree_rgc = f_rgc.Get("PhysicsEvents")

if not tree_rga:
    raise RuntimeError("PhysicsEvents tree missing in RGA file")
#endif

if not tree_rgc:
    raise RuntimeError("PhysicsEvents tree missing in RGC file")
#endif

# ------------------------------------------------
# histograms
# ------------------------------------------------

nbins = 200
xmin  = 0.4
xmax  = 1.4

h_rga = ROOT.TH1F("h_rga",";Mx2 (GeV^{2});Normalized Counts",nbins,xmin,xmax)
h_rgc = ROOT.TH1F("h_rgc",";Mx2 (GeV^{2});Normalized Counts",nbins,xmin,xmax)

# ------------------------------------------------
# fill histograms
# ------------------------------------------------

tree_rga.Draw("Mx2>>h_rga","","goff")
tree_rgc.Draw("Mx2>>h_rgc","","goff")

# ------------------------------------------------
# normalize to integral
# ------------------------------------------------

int_rga = h_rga.Integral()
int_rgc = h_rgc.Integral()

if int_rga > 0:
    h_rga.Scale(1.0 / int_rga)
#endif

if int_rgc > 0:
    h_rgc.Scale(1.0 / int_rgc)
#endif

# ------------------------------------------------
# style
# ------------------------------------------------

h_rga.SetLineColor(ROOT.kBlack)
h_rga.SetLineWidth(2)

h_rgc.SetLineColor(ROOT.kRed)
h_rgc.SetLineWidth(2)

# ------------------------------------------------
# determine y-axis max
# ------------------------------------------------

max_val = max(h_rga.GetMaximum(), h_rgc.GetMaximum())
h_rga.SetMaximum(1.2 * max_val)

# ------------------------------------------------
# canvas
# ------------------------------------------------

c = ROOT.TCanvas("c","Mx2 comparison",800,600)

h_rga.Draw("hist")
h_rgc.Draw("hist same")

# ------------------------------------------------
# legend
# ------------------------------------------------

leg = ROOT.TLegend(0.65,0.75,0.88,0.88)
leg.AddEntry(h_rga,"RGA","l")
leg.AddEntry(h_rgc,"RGC","l")
leg.Draw()

# ------------------------------------------------
# output
# ------------------------------------------------

os.makedirs("output", exist_ok=True)

c.SaveAs("output/epi+_Mx2_peak.png")