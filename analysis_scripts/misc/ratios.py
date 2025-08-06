#!/usr/bin/env python3
import ROOT
import math

# Turn off the default stats box
ROOT.gStyle.SetOptStat(0)

# Open the two ROOT files (dipion)
file1 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root")
file2 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root")
tree1 = file1.Get("PhysicsEvents")
tree2 = file2.Get("PhysicsEvents")

# Helper to draw ratio with labels
def draw_ratio(h_num, h_den, n_bins, x_min, x_max, xlabel):
    hr = ROOT.TH1F("", "", n_bins, x_min, x_max)
    for i in range(1, n_bins+1):
        n = h_num.GetBinContent(i)
        d = h_den.GetBinContent(i)
        if d > 0:
            err_n = math.sqrt(n)
            err_d = math.sqrt(d)
            r = n/d
            err = math.sqrt((err_n/d)**2 + (n*err_d/(d*d))**2)
        else:
            r, err = 0.0, 0.0
        hr.SetBinContent(i, r)
        hr.SetBinError(i, err)
    hr.SetLineColor(ROOT.kBlue)
    hr.SetLineWidth(2)
    hr.SetTitle(f"Ratio; {xlabel}; cj13.0.3 / cj11.2.0")
    hr.GetYaxis().SetRangeUser(0.8, 1.2)
    return hr

# Create all histograms (70 bins each)
h_ep1  = ROOT.TH1F("h_ep1",  "e_{p}; e_{p} (GeV); counts", 70, 1.0, 4.5)
h_ep2  = ROOT.TH1F("h_ep2",  "e_{p}; e_{p} (GeV); counts", 70, 1.0, 4.5)
h_e1   = ROOT.TH1F("h_e1",   "e_{θ}; e_{θ} (deg); counts", 70, 5.0, 35.0)
h_e2   = ROOT.TH1F("h_e2",   "e_{θ}; e_{θ} (deg); counts", 70, 5.0, 35.0)

h_p1p1 = ROOT.TH1F("h_p1p1", "p1_{p}; p1_{p} (GeV); counts", 70, 1.0, 4.5)
h_p1p2 = ROOT.TH1F("h_p1p2", "p1_{p}; p1_{p} (GeV); counts", 70, 1.0, 4.5)
h_p11  = ROOT.TH1F("h_p11",  "#theta_{π+}; #theta_{π+} (deg); counts", 70, 5.0, 35.0)
h_p12  = ROOT.TH1F("h_p12",  "#theta_{π+}; #theta_{π+} (deg); counts", 70, 5.0, 35.0)

h_p2p1 = ROOT.TH1F("h_p2p1", "p2_{p}; p2_{p} (GeV); counts", 70, 1.0, 4.5)
h_p2p2 = ROOT.TH1F("h_p2p2", "p2_{p}; p2_{p} (GeV); counts", 70, 1.0, 4.5)
h_p21  = ROOT.TH1F("h_p21",  "#theta_{π-}; #theta_{π-} (deg); counts", 70, 5.0, 35.0)
h_p22  = ROOT.TH1F("h_p22",  "#theta_{π-}; #theta_{π-} (deg); counts", 70, 5.0, 35.0)

# Fill histograms
tree1.Draw("e_p >> h_ep1")
tree2.Draw("e_p >> h_ep2")
tree1.Draw("e_theta*180./TMath::Pi() >> h_e1")
tree2.Draw("e_theta*180./TMath::Pi() >> h_e2")

tree1.Draw("p1_p >> h_p1p1")
tree2.Draw("p1_p >> h_p1p2")
tree1.Draw("p1_theta*180./TMath::Pi() >> h_p11")
tree2.Draw("p1_theta*180./TMath::Pi() >> h_p12")

tree1.Draw("p2_p >> h_p2p1")
tree2.Draw("p2_p >> h_p2p2")
tree1.Draw("p2_theta*180./TMath::Pi() >> h_p21")
tree2.Draw("p2_theta*180./TMath::Pi() >> h_p22")

# Style: blue for cj11, red for cj13
for h in (h_ep1, h_e1, h_p1p1, h_p11, h_p2p1, h_p21):
    h.SetLineColor(ROOT.kBlue); h.SetLineWidth(2)
for h in (h_ep2, h_e2, h_p1p2, h_p12, h_p2p2, h_p22):
    h.SetLineColor(ROOT.kRed);  h.SetLineWidth(2)

# Integrals
I = lambda h: int(h.Integral())
tot_ep1, tot_ep2 = I(h_ep1), I(h_ep2)
tot_e1,  tot_e2  = I(h_e1),  I(h_e2)
tot_p1p1, tot_p1p2 = I(h_p1p1), I(h_p1p2)
tot_p11, tot_p12 = I(h_p11), I(h_p12)
tot_p2p1, tot_p2p2 = I(h_p2p1), I(h_p2p2)
tot_p21, tot_p22 = I(h_p21), I(h_p22)

# Precompute maxima
max_ep  = max(h_ep1.GetMaximum(), h_ep2.GetMaximum())
max_e   = max(h_e1.GetMaximum(),  h_e2.GetMaximum())
max_p1p = max(h_p1p1.GetMaximum(),h_p1p2.GetMaximum())
max_p11 = max(h_p11.GetMaximum(), h_p12.GetMaximum())
max_p2p = max(h_p2p1.GetMaximum(),h_p2p2.GetMaximum())
max_p21 = max(h_p21.GetMaximum(), h_p22.GetMaximum())

# Canvas 1: electron
c1 = ROOT.TCanvas("c1","electron",1200,1000); c1.Divide(2,2)
c1.cd(1); ROOT.gPad.SetLeftMargin(0.15)
h_ep1.Draw("HIST"); h_ep2.Draw("HIST SAME")
h_ep1.GetYaxis().SetRangeUser(0,1.2*max_ep)
leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
leg.AddEntry(h_ep1,f"cj11.2.0 (N={tot_ep1})","l")
leg.AddEntry(h_ep2,f"cj13.0.3 (N={tot_ep2})","l"); leg.Draw()
c1.cd(2)
draw_ratio(h_ep2,h_ep1,70,1.0,4.5,"e_{p} (GeV)").Draw("E1")
c1.cd(3); ROOT.gPad.SetLeftMargin(0.15)
h_e1.Draw("HIST"); h_e2.Draw("HIST SAME")
h_e1.GetYaxis().SetRangeUser(0,1.2*max_e)
leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
leg.AddEntry(h_e1,f"cj11.2.0 (N={tot_e1})","l")
leg.AddEntry(h_e2,f"cj13.0.3 (N={tot_e2})","l"); leg.Draw()
c1.cd(4)
draw_ratio(h_e2,h_e1,70,5.0,35.0,"e_{θ} (deg)").Draw("E1")
c1.SaveAs("/u/home/thayward/ep_comparison.pdf")

# --- Canvas 2: pi+ ---
c2 = ROOT.TCanvas("c2", "#pi+ comparison", 1200, 1000)
c2.Divide(2, 2)

c2.cd(1); ROOT.gPad.SetLeftMargin(0.15)
h_p1p1.Draw("HIST"); h_p1p2.Draw("HIST SAME")
h_p1p1.GetYaxis().SetRangeUser(0,1.2*max_p1p)
leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
leg.AddEntry(h_p1p1,f"cj11.2.0 (N={tot_p1p1})","l")
leg.AddEntry(h_p1p2,f"cj13.0.3 (N={tot_p1p2})","l"); leg.Draw()

c2.cd(2)
draw_ratio(h_p1p2,h_p1p1,70,1.0,4.5,"p1_{p} (GeV)").Draw("E1")

c2.cd(3); ROOT.gPad.SetLeftMargin(0.15)
h_p11.Draw("HIST"); h_p12.Draw("HIST SAME")
h_p11.GetYaxis().SetRangeUser(0,1.2*max_p11)
leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
leg.AddEntry(h_p11,f"cj11.2.0 (N={tot_p11})","l")
leg.AddEntry(h_p12,f"cj13.0.3 (N={tot_p12})","l"); leg.Draw()

c2.cd(4)
draw_ratio(h_p12,h_p11,70,5.0,35.0,"#theta_{π+} (deg)").Draw("E1")

c2.SaveAs("/u/home/thayward/p1_comparison.pdf")


# --- Canvas 3: pi- ---
c3 = ROOT.TCanvas("c3", "#pi- comparison", 1200, 1000)
c3.Divide(2, 2)

c3.cd(1); ROOT.gPad.SetLeftMargin(0.15)
h_p2p1.Draw("HIST"); h_p2p2.Draw("HIST SAME")
h_p2p1.GetYaxis().SetRangeUser(0,1.2*max_p2p)
leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
leg.AddEntry(h_p2p1,f"cj11.2.0 (N={tot_p2p1})","l")
leg.AddEntry(h_p2p2,f"cj13.0.3 (N={tot_p2p2})","l"); leg.Draw()

c3.cd(2)
draw_ratio(h_p2p2,h_p2p1,70,1.0,4.5,"p2_{p} (GeV)").Draw("E1")

c3.cd(3); ROOT.gPad.SetLeftMargin(0.15)
h_p21.Draw("HIST"); h_p22.Draw("HIST SAME")
h_p21.GetYaxis().SetRangeUser(0,1.2*max_p21)
leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
leg.AddEntry(h_p21,f"cj11.2.0 (N={tot_p21})","l")
leg.AddEntry(h_p22,f"cj13.0.3 (N={tot_p22})","l"); leg.Draw()

c3.cd(4)
draw_ratio(h_p22,h_p21,70,5.0,35.0,"#theta_{π-} (deg)").Draw("E1")

c3.SaveAs("/u/home/thayward/p2_comparison.pdf")