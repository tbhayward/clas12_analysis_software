#!/usr/bin/env python3
import ROOT
import math

# Turn off the default stats box
ROOT.gStyle.SetOptStat(0)

# Open the two ROOT files (dipion)
file1 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root")
file2 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root")

# Load the PhysicsEvents trees
tree1 = file1.Get("PhysicsEvents")
tree2 = file2.Get("PhysicsEvents")

# Create histograms for e_p from 1 to 4.5 GeV (35 bins of 0.1 GeV each)
h1 = ROOT.TH1F("h1", "e_{p} distribution; e_{p} (GeV); counts", 35, 1.0, 4.5)
h2 = ROOT.TH1F("h2", "e_{p} distribution; e_{p} (GeV); counts", 35, 1.0, 4.5)

# Create histograms for e_theta in degrees from 5 to 35 (30 bins of 1° each)
h_e1 = ROOT.TH1F("h_e1", "e_{θ} distribution; e_{θ} (deg); counts", 30, 5.0, 35.0)
h_e2 = ROOT.TH1F("h_e2", "e_{θ} distribution; e_{θ} (deg); counts", 30, 5.0, 35.0)

# Fill the e_p histograms
tree1.Draw("e_p >> h1")
tree2.Draw("e_p >> h2")

# Fill the e_theta histograms (convert rad→deg)
tree1.Draw("e_theta*180./TMath::Pi() >> h_e1")
tree2.Draw("e_theta*180./TMath::Pi() >> h_e2")

# Style the histograms
for hist, color in [(h1, ROOT.kBlack), (h_e1, ROOT.kBlue)]:
    hist.SetLineColor(color)
    hist.SetLineWidth(2)
for hist, color in [(h2, ROOT.kRed), (h_e2, ROOT.kGreen+2)]:
    hist.SetLineColor(color)
    hist.SetLineWidth(2)

# Calculate total counts
total1 = h1.Integral()
total2 = h2.Integral()
total_e1 = h_e1.Integral()
total_e2 = h_e2.Integral()

# Create a 2x2 canvas
c = ROOT.TCanvas("c", "comparison", 1200, 1000)
c.Divide(2, 2)

# --- Top-left: e_p histograms ---
c.cd(1)
h1.Draw("HIST")
h2.Draw("HIST SAME")
max_bin = max(h1.GetMaximum(), h2.GetMaximum())
h1.GetYaxis().SetRangeUser(0, 1.2 * max_bin)
leg1 = ROOT.TLegend(0.60, 0.65, 0.88, 0.88)
leg1.AddEntry(h1, f"cj 11.2.0 (N={int(total1)})", "l")
leg1.AddEntry(h2, f"cj 13.0.3 (N={int(total2)})", "l")
leg1.Draw()

# --- Top-right: e_p ratio ---
c.cd(2)
h_ratio = ROOT.TH1F("h_ratio", "Ratio; e_{p} (GeV); cj13.0.3 / cj11.2.0", 35, 1.0, 4.5)
for i in range(1, h1.GetNbinsX()+1):
    n1 = h1.GetBinContent(i)
    n2 = h2.GetBinContent(i)
    if n1 > 0:
        err1 = math.sqrt(n1)
        err2 = math.sqrt(n2)
        r = n2/n1
        err = math.sqrt((err2/n1)**2 + (n2*err1/(n1*n1))**2)
    else:
        r, err = 0.0, 0.0
    h_ratio.SetBinContent(i, r)
    h_ratio.SetBinError(i, err)
h_ratio.SetLineColor(ROOT.kBlue)
h_ratio.SetLineWidth(2)
h_ratio.GetYaxis().SetRangeUser(0.8, 1.2)
h_ratio.Draw("E1")

# --- Bottom-left: e_theta histograms ---
c.cd(3)
h_e1.Draw("HIST")
h_e2.Draw("HIST SAME")
max_bin_e = max(h_e1.GetMaximum(), h_e2.GetMaximum())
h_e1.GetYaxis().SetRangeUser(0, 1.2 * max_bin_e)
leg2 = ROOT.TLegend(0.60, 0.65, 0.88, 0.88)
leg2.AddEntry(h_e1, f"cj 11.2.0 (N={int(total_e1)})", "l")
leg2.AddEntry(h_e2, f"cj 13.0.3 (N={int(total_e2)})", "l")
leg2.Draw()

# --- Bottom-right: e_theta ratio ---
c.cd(4)
h_ratio_t = ROOT.TH1F("h_ratio_t", "Ratio; e_{θ} (deg); cj13.0.3 / cj11.2.0", 30, 5.0, 35.0)
for i in range(1, h_e1.GetNbinsX()+1):
    n1 = h_e1.GetBinContent(i)
    n2 = h_e2.GetBinContent(i)
    if n1 > 0:
        err1 = math.sqrt(n1)
        err2 = math.sqrt(n2)
        r = n2/n1
        err = math.sqrt((err2/n1)**2 + (n2*err1/(n1*n1))**2)
    else:
        r, err = 0.0, 0.0
    h_ratio_t.SetBinContent(i, r)
    h_ratio_t.SetBinError(i, err)
h_ratio_t.SetLineColor(ROOT.kBlue)
h_ratio_t.SetLineWidth(2)
h_ratio_t.GetYaxis().SetRangeUser(0.8, 1.2)
h_ratio_t.Draw("E1")

# Save the canvas
c.SaveAs("/u/home/thayward/ep_etheta_comparison.pdf")