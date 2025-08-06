#!/usr/bin/env python3
import ROOT
import math

# Turn off the default stats box
ROOT.gStyle.SetOptStat(0)

# Open the two ROOT files (dipion)
file11 = ROOT.TFile.Open(
    "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root")
file13 = ROOT.TFile.Open(
    "/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root")
tree11 = file11.Get("PhysicsEvents")
tree13 = file13.Get("PhysicsEvents")

# Global counter so each ratio hist has a unique name
_ratio_counter = 0

def draw_ratio(h_num, h_den, nbins, xmin, xmax, xlabel):
    """Make a ratio histogram with Poisson errors and unique name."""
    global _ratio_counter
    name = f"ratio_{_ratio_counter}"
    _ratio_counter += 1

    hr = ROOT.TH1F(name, f"Ratio; {xlabel}; ratio", nbins, xmin, xmax)
    for i in range(1, nbins+1):
        n = h_num.GetBinContent(i)
        d = h_den.GetBinContent(i)
        if d > 0:
            err_n = math.sqrt(n)
            err_d = math.sqrt(d)
            r = n / d
            err = math.sqrt((err_n/d)**2 + (n*err_d/(d*d))**2)
        else:
            r, err = 0.0, 0.0
        hr.SetBinContent(i, r)
        hr.SetBinError(i, err)

    hr.SetLineColor(ROOT.kBlue)
    hr.SetLineWidth(2)
    hr.GetYaxis().SetRangeUser(0.8, 1.2)
    return hr

def make_canvas(p_branch, theta_branch, p_label, theta_label,
                canvas_name, out_pdf, p_min, p_max, t_min, t_max):
    """
    Build a 2×2 canvas:
      [0,0] momentum histograms
      [0,1] momentum ratio
      [1,0] theta histograms
      [1,1] theta ratio
    """
    # binning
    p_bins = 70
    t_bins = 70

    # 1) Create histograms
    h_p11 = ROOT.TH1F("h_p11",  f"{p_label}; {p_label} (GeV); counts",
                      p_bins, p_min, p_max)
    h_p13 = ROOT.TH1F("h_p13",  f"{p_label}; {p_label} (GeV); counts",
                      p_bins, p_min, p_max)
    h_t11 = ROOT.TH1F("h_t11",  f"{theta_label}; {theta_label} (deg); counts",
                      t_bins, t_min, t_max)
    h_t13 = ROOT.TH1F("h_t13",  f"{theta_label}; {theta_label} (deg); counts",
                      t_bins, t_min, t_max)

    # 2) Fill
    tree11.Draw(f"{p_branch} >> h_p11")
    tree13.Draw(f"{p_branch} >> h_p13")
    tree11.Draw(f"{theta_branch}*180./TMath::Pi() >> h_t11")
    tree13.Draw(f"{theta_branch}*180./TMath::Pi() >> h_t13")

    # 3) Style
    for h in (h_p11, h_t11):
        h.SetLineColor(ROOT.kBlue)
        h.SetLineWidth(2)
    for h in (h_p13, h_t13):
        h.SetLineColor(ROOT.kRed)
        h.SetLineWidth(2)

    # 4) Legends and maxima
    tot_p11 = int(h_p11.Integral())
    tot_p13 = int(h_p13.Integral())
    tot_t11 = int(h_t11.Integral())
    tot_t13 = int(h_t13.Integral())
    max_p = max(h_p11.GetMaximum(), h_p13.GetMaximum())
    max_t = max(h_t11.GetMaximum(), h_t13.GetMaximum())

    # 5) Build canvas
    c = ROOT.TCanvas(canvas_name, canvas_name, 1200, 1000)
    c.Divide(2,2)

    # -- pad 1: momentum histograms
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    h_p11.Draw("HIST")
    h_p13.Draw("HIST SAME")
    h_p11.GetYaxis().SetRangeUser(0, 1.2*max_p)
    leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
    leg.AddEntry(h_p11, f"cj 11.2.0 (N={tot_p11})", "l")
    leg.AddEntry(h_p13, f"cj 13.0.3 (N={tot_p13})", "l")
    leg.Draw()

    # -- pad 2: momentum ratio
    c.cd(2)
    # ensure pad is active
    ROOT.gPad.SetLeftMargin(0.15)
    r_p = draw_ratio(h_p13, h_p11, p_bins, p_min, p_max,
                     f"{p_label} (GeV)")
    r_p.Draw("E1")

    # -- pad 3: theta histograms
    c.cd(3)
    ROOT.gPad.SetLeftMargin(0.15)
    h_t11.Draw("HIST")
    h_t13.Draw("HIST SAME")
    h_t11.GetYaxis().SetRangeUser(0, 1.2*max_t)
    leg2 = ROOT.TLegend(0.6,0.65,0.88,0.88)
    leg2.AddEntry(h_t11, f"cj 11.2.0 (N={tot_t11})", "l")
    leg2.AddEntry(h_t13, f"cj 13.0.3 (N={tot_t13})", "l")
    leg2.Draw()

    # -- pad 4: theta ratio
    c.cd(4)
    ROOT.gPad.SetLeftMargin(0.15)
    r_t = draw_ratio(h_t13, h_t11, t_bins, t_min, t_max,
                     f"{theta_label} (deg)")
    r_t.Draw("E1")

    c.SaveAs(out_pdf)


# electron
make_canvas(
  "e_p", "e_theta",
  "e_{p}", "#theta_{e}",
  "c_electron", "/u/home/thayward/ep_comparison.pdf",
   p_min=1.0, p_max=4.5,  t_min=5.0, t_max=35.0
)

# pi+
make_canvas(
  "p1_p", "p1_theta",
  "#pi^{+}_{p}", "#theta_{#pi^{+}}",
  "c_pi_plus", "/u/home/thayward/p1_comparison.pdf",
   p_min=1.0, p_max=6.0,  t_min=0.0, t_max=60.0
)

# pi-
make_canvas(
  "p2_p", "p2_theta",
  "#pi^{-}_{p}", "#theta_{#pi^{-}}",
  "c_pi_minus", "/u/home/thayward/p2_comparison.pdf",
   p_min=1.0, p_max=6.0,  t_min=0.0, t_max=60.0
)