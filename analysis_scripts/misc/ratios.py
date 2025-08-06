#!/usr/bin/env python3
import ROOT
import math

# Turn off the default stats box
ROOT.gStyle.SetOptStat(0)

# Open the two ROOT files (dipion)
file11 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root")
file13 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root")
tree11 = file11.Get("PhysicsEvents")
tree13 = file13.Get("PhysicsEvents")

def draw_ratio(h_num, h_den, n_bins, x_min, x_max, xlabel):
    """Return a ratio histogram with Poisson errors and axis labels."""
    hr = ROOT.TH1F("", "", n_bins, x_min, x_max)
    for i in range(1, n_bins+1):
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
    # Title format: "hist title; x-axis title; y-axis title"
    hr.SetTitle(f"Ratio; {xlabel}; ratio")
    hr.GetYaxis().SetRangeUser(0.8, 1.2)
    return hr

def make_canvas(p_branch, theta_branch, p_label, theta_label, canvas_name, out_pdf):
    """
    Create a 2×2 canvas:
      top-left: momentum histograms,
      top-right: momentum ratio,
      bottom-left: theta histograms,
      bottom-right: theta ratio.
    """
    # Settings
    p_bins, p_min, p_max = 70, 1.0, 4.5
    t_bins, t_min, t_max = 70, 5.0, 35.0

    # Histograms
    h_p11 = ROOT.TH1F("h_p11",  f"{p_label}; {p_label} (GeV); counts", p_bins, p_min, p_max)
    h_p13 = ROOT.TH1F("h_p13",  f"{p_label}; {p_label} (GeV); counts", p_bins, p_min, p_max)
    h_t11 = ROOT.TH1F("h_t11",  f"{theta_label}; {theta_label} (deg); counts", t_bins, t_min, t_max)
    h_t13 = ROOT.TH1F("h_t13",  f"{theta_label}; {theta_label} (deg); counts", t_bins, t_min, t_max)

    # Fill
    tree11.Draw(f"{p_branch} >> h_p11")
    tree13.Draw(f"{p_branch} >> h_p13")
    tree11.Draw(f"{theta_branch}*180./TMath::Pi() >> h_t11")
    tree13.Draw(f"{theta_branch}*180./TMath::Pi() >> h_t13")

    # Style
    for h in (h_p11, h_t11):
        h.SetLineColor(ROOT.kBlue)
        h.SetLineWidth(2)
    for h in (h_p13, h_t13):
        h.SetLineColor(ROOT.kRed)
        h.SetLineWidth(2)

    # Integrals
    tot_p11 = int(h_p11.Integral())
    tot_p13 = int(h_p13.Integral())
    tot_t11 = int(h_t11.Integral())
    tot_t13 = int(h_t13.Integral())

    # Maxima
    max_p = max(h_p11.GetMaximum(), h_p13.GetMaximum())
    max_t = max(h_t11.GetMaximum(), h_t13.GetMaximum())

    # Canvas
    c = ROOT.TCanvas(canvas_name, canvas_name, 1200, 1000)
    c.Divide(2, 2)

    # Top-left: p histograms
    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    h_p11.Draw("HIST")
    h_p13.Draw("HIST SAME")
    h_p11.GetYaxis().SetRangeUser(0, 1.2 * max_p)
    leg = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
    leg.AddEntry(h_p11, f"cj 11.2.0 (N={tot_p11})", "l")
    leg.AddEntry(h_p13, f"cj 13.0.3 (N={tot_p13})", "l")
    leg.Draw()

    # Top-right: p ratio
    c.cd(2)
    draw_ratio(h_p13, h_p11, p_bins, p_min, p_max, p_label + " (GeV)").Draw("E1")

    # Bottom-left: theta histograms
    c.cd(3)
    ROOT.gPad.SetLeftMargin(0.15)
    h_t11.Draw("HIST")
    h_t13.Draw("HIST SAME")
    h_t11.GetYaxis().SetRangeUser(0, 1.2 * max_t)
    leg = ROOT.TLegend(0.6, 0.65, 0.88, 0.88)
    leg.AddEntry(h_t11, f"cj 11.2.0 (N={tot_t11})", "l")
    leg.AddEntry(h_t13, f"cj 13.0.3 (N={tot_t13})", "l")
    leg.Draw()

    # Bottom-right: theta ratio
    c.cd(4)
    draw_ratio(h_t13, h_t11, t_bins, t_min, t_max, theta_label + " (deg)").Draw("E1")

    c.SaveAs(out_pdf)


# Electron: use e_p and e_theta with "#theta_e"
make_canvas(
    p_branch="e_p",
    theta_branch="e_theta",
    p_label="e_{p}",
    theta_label="#theta_e",
    canvas_name="c_electron",
    out_pdf="/u/home/thayward/ep_comparison.pdf"
)

# Pi+: use p1_p and p1_theta with "#pi^{+}"
make_canvas(
    p_branch="p1_p",
    theta_branch="p1_theta",
    p_label="#pi^{+}_p",
    theta_label="#theta_{#pi^{+}}",
    canvas_name="c_pi_plus",
    out_pdf="/u/home/thayward/p1_comparison.pdf"
)

# Pi-: use p2_p and p2_theta with "#pi^{-}"
make_canvas(
    p_branch="p2_p",
    theta_branch="p2_theta",
    p_label="#pi^{-}_p",
    theta_label="#theta_{#pi^{-}}",
    canvas_name="c_pi_minus",
    out_pdf="/u/home/thayward/p2_comparison.pdf"
)