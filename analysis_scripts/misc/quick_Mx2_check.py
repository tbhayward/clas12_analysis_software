import ROOT
import os

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# input files
# ------------------------------------------------

file_rga = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_fa18_inb_epi+.root"
file_mc  = "/work/clas12/thayward/CLAS12_exclusive/enpi+/mc/clasdis_rga_fa18_inb_epi+X.root"

# ------------------------------------------------
# open files
# ------------------------------------------------

f_rga = ROOT.TFile.Open(file_rga)
if not f_rga or f_rga.IsZombie():
    raise RuntimeError("Failed to open file: " + file_rga)
#endif

f_mc = ROOT.TFile.Open(file_mc)
if not f_mc or f_mc.IsZombie():
    raise RuntimeError("Failed to open file: " + file_mc)
#endif

# ------------------------------------------------
# get trees
# ------------------------------------------------

tree_rga = f_rga.Get("PhysicsEvents")
tree_mc  = f_mc.Get("PhysicsEvents")

if not tree_rga:
    raise RuntimeError("PhysicsEvents tree missing in RGA file")
#endif

if not tree_mc:
    raise RuntimeError("PhysicsEvents tree missing in clasdis file")
#endif

# ------------------------------------------------
# histogram definitions
# ------------------------------------------------

nbins = 200
xmin  = 0.4
xmax  = 1.4

h_rga = ROOT.TH1F("h_rga",";Mx2 (GeV^{2});Normalized Counts",nbins,xmin,xmax)
h_mc  = ROOT.TH1F("h_mc",";Mx2 (GeV^{2});Normalized Counts",nbins,xmin,xmax)

# ------------------------------------------------
# fill histograms
# ------------------------------------------------

tree_rga.Draw("Mx2>>h_rga","","goff")
tree_mc.Draw("Mx2>>h_mc","","goff")

# ------------------------------------------------
# normalize to integral
# ------------------------------------------------

int_rga = h_rga.Integral()
int_mc  = h_mc.Integral()

if int_rga > 0:
    h_rga.Scale(1.0 / int_rga)
#endif

if int_mc > 0:
    h_mc.Scale(1.0 / int_mc)
#endif

# ------------------------------------------------
# style
# ------------------------------------------------

h_rga.SetLineColor(ROOT.kBlack)
h_rga.SetLineWidth(2)

h_mc.SetLineColor(ROOT.kRed)
h_mc.SetLineWidth(2)

# ------------------------------------------------
# determine y-axis range
# ------------------------------------------------

max_val = max(h_rga.GetMaximum(), h_mc.GetMaximum())

h_rga.SetMaximum(1.2 * max_val)
h_rga.SetMinimum(1e-6)

# ------------------------------------------------
# canvas
# ------------------------------------------------

c = ROOT.TCanvas("c","Mx2 comparison",800,600)
c.SetLogy()

h_rga.Draw("hist")
h_mc.Draw("hist same")

# ------------------------------------------------
# legend
# ------------------------------------------------

leg = ROOT.TLegend(0.65,0.75,0.88,0.88)
leg.AddEntry(h_rga,"RGA","l")
leg.AddEntry(h_mc,"clasdis","l")
leg.Draw()

# ------------------------------------------------
# output
# ------------------------------------------------

os.makedirs("output", exist_ok=True)

c.SaveAs("output/epi+_Mx2_peak.png")