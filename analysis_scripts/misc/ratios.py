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

# --- Create histograms ---

# e_p histograms: 70 bins from 1 to 4.5 GeV
h_ep1 = ROOT.TH1F("h_ep1", "e_{p} distribution; e_{p} (GeV); counts", 70, 1.0, 4.5)
h_ep2 = ROOT.TH1F("h_ep2", "e_{p} distribution; e_{p} (GeV); counts", 70, 1.0, 4.5)

# e_theta histograms: 70 bins from 5 to 35 degrees
h_e1 = ROOT.TH1F("h_e1", "e_{θ} distribution; e_{θ} (deg); counts", 70, 5.0, 35.0)
h_e2 = ROOT.TH1F("h_e2", "e_{θ} distribution; e_{θ} (deg); counts", 70, 5.0, 35.0)

# p1_theta histograms: 70 bins from 5 to 35 degrees
h_p11 = ROOT.TH1F("h_p11", "#theta_{π+} distribution; #theta_{π+} (deg); counts", 70, 5.0, 35.0)
h_p12 = ROOT.TH1F("h_p12", "#theta_{π+} distribution; #theta_{π+} (deg); counts", 70, 5.0, 35.0)

# p2_theta histograms: 70 bins from 5 to 35 degrees
h_p21 = ROOT.TH1F("h_p21", "#theta_{π-} distribution; #theta_{π-} (deg); counts", 70, 5.0, 35.0)
h_p22 = ROOT.TH1F("h_p22", "#theta_{π-} distribution; #theta_{π-} (deg); counts", 70, 5.0, 35.0)

# --- Fill histograms ---

tree1.Draw("e_p >> h_ep1")
tree2.Draw("e_p >> h_ep2")

tree1.Draw("e_theta*180./TMath::Pi() >> h_e1")
tree2.Draw("e_theta*180./TMath::Pi() >> h_e2")

tree1.Draw("p1_theta*180./TMath::Pi() >> h_p11")
tree2.Draw("p1_theta*180./TMath::Pi() >> h_p12")

tree1.Draw("p2_theta*180./TMath::Pi() >> h_p21")
tree2.Draw("p2_theta*180./TMath::Pi() >> h_p22")

# --- Style histograms: blue for cj11, red for cj13 ---

for hist in (h_ep1, h_e1, h_p11, h_p21):
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(2)

for hist in (h_ep2, h_e2, h_p12, h_p22):
    hist.SetLineColor(ROOT.kRed)
    hist.SetLineWidth(2)

# --- Integrals for legends ---

tot_ep1 = h_ep1.Integral()
tot_ep2 = h_ep2.Integral()
tot_e1  = h_e1.Integral()
tot_e2  = h_e2.Integral()
tot_p11 = h_p11.Integral()
tot_p12 = h_p12.Integral()
tot_p21 = h_p21.Integral()
tot_p22 = h_p22.Integral()

# Function to draw ratio with Poisson errors
def draw_ratio(h_num, h_den, n_bins, x_min, x_max):
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
    hr.GetYaxis().SetRangeUser(0.8, 1.2)
    return hr

# --- Canvas 1: e_p & e_theta ---

c1 = ROOT.TCanvas("c1", "e_p and e_θ comparison", 1200, 1000)
c1.Divide(2, 2)

# Top-left: e_p histograms (add left padding)
c1.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
h_ep1.Draw("HIST")
h_ep2.Draw("HIST SAME")
max_ep = max(h_ep1.GetMaximum(), h_ep2.GetMaximum())
h_ep1.GetYaxis().SetRangeUser(0, 1.2*max_ep)
leg1 = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
leg1.AddEntry(h_ep1, f"cj 11.2.0 (N={int(tot_ep1)})", "l")
leg1.AddEntry(h_ep2, f"cj 13.0.3 (N={int(tot_ep2)})", "l")
leg1.Draw()

# Top-right: e_p ratio
c1.cd(2)
r_ep = draw_ratio(h_ep2, h_ep1, 70, 1.0, 4.5)
r_ep.Draw("E1")

# Bottom-left: e_theta histograms (add left padding)
c1.cd(3)
ROOT.gPad.SetLeftMargin(0.15)
h_e1.Draw("HIST")
h_e2.Draw("HIST SAME")
max_e = max(h_e1.GetMaximum(), h_e2.GetMaximum())
h_e1.GetYaxis().SetRangeUser(0, 1.2*max_e)
leg2 = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
leg2.AddEntry(h_e1, f"cj 11.2.0 (N={int(tot_e1)})", "l")
leg2.AddEntry(h_e2, f"cj 13.0.3 (N={int(tot_e2)})", "l")
leg2.Draw()

# Bottom-right: e_theta ratio
c1.cd(4)
r_e = draw_ratio(h_e2, h_e1, 70, 5.0, 35.0)
r_e.Draw("E1")

c1.SaveAs("/u/home/thayward/ep_etheta_comparison.pdf")

# --- Canvas 2: e_p & p1_theta ---

c2 = ROOT.TCanvas("c2", "e_p and #theta_{π+} comparison", 1200, 1000)
c2.Divide(2, 2)

# Top-left: e_p (add left padding)
c2.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
h_ep1.Draw("HIST")
h_ep2.Draw("HIST SAME")
h_ep1.GetYaxis().SetRangeUser(0, 1.2*max_ep)
leg3 = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
leg3.AddEntry(h_ep1, f"cj 11.2.0 (N={int(tot_ep1)})", "l")
leg3.AddEntry(h_ep2, f"cj 13.0.3 (N={int(tot_ep2)})", "l")
leg3.Draw()

# Top-right: e_p ratio
c2.cd(2)
r_ep.Draw("E1")

# Bottom-left: p1_theta histograms (add left padding)
c2.cd(3)
ROOT.gPad.SetLeftMargin(0.15)
h_p11.Draw("HIST")
h_p12.Draw("HIST SAME")
max_p11 = max(h_p11.GetMaximum(), h_p12.GetMaximum())
h_p11.GetYaxis().SetRangeUser(0, 1.2*max_p11)
leg4 = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
leg4.AddEntry(h_p11, f"cj 11.2.0 (N={int(tot_p11)})", "l")
leg4.AddEntry(h_p12, f"cj 13.0.3 (N={int(tot_p12)})", "l")
leg4.Draw()

# Bottom-right: p1_theta ratio
c2.cd(4)
r_p1 = draw_ratio(h_p12, h_p11, 70, 5.0, 35.0)
r_p1.Draw("E1")

c2.SaveAs("/u/home/thayward/ep_p1theta_comparison.pdf")

# --- Canvas 3: e_p & p2_theta ---

c3 = ROOT.TCanvas("c3", "e_p and #theta_{π-} comparison", 1200, 1000)
c3.Divide(2, 2)

# Top-left: e_p (add left padding)
c3.cd(1)
ROOT.gPad.SetLeftMargin(0.15)
h_ep1.Draw("HIST")
h_ep2.Draw("HIST SAME")
h_ep1.GetYaxis().SetRangeUser(0, 1.2*max_ep)
leg5 = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
leg5.AddEntry(h_ep1, f"cj 11.2.0 (N={int(tot_ep1)})", "l")
leg5.AddEntry(h_ep2, f"cj 13.0.3 (N={int(tot_ep2)})", "l")
leg5.Draw()

# Top-right: e_p ratio
c3.cd(2)
r_ep.Draw("E1")

# Bottom-left: p2_theta histograms (add left padding)
c3.cd(3)
ROOT.gPad.SetLeftMargin(0.15)
h_p21.Draw("HIST")
h_p22.Draw("HIST SAME")
max_p21 = max(h_p21.GetMaximum(), h_p22.GetMaximum())
h_p21.GetYaxis().SetRangeUser(0, 1.2*max_p21)
leg6 = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
leg6.AddEntry(h_p21, f"cj 11.2.0 (N={int(tot_p21)})", "l")
leg6.AddEntry(h_p22, f"cj 13.0.3 (N={int(tot_p22)})", "l")
leg6.Draw()

# Bottom-right: p2_theta ratio
c3.cd(4)
r_p2 = draw_ratio(h_p22, h_p21, 70, 5.0, 35.0)
r_p2.Draw("E1")

c3.SaveAs("/u/home/thayward/ep_p2theta_comparison.pdf")