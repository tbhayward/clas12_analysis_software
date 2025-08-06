#!/usr/bin/env python3
import ROOT
import math

# Turn off the default stats box
ROOT.gStyle.SetOptStat(0)

# Open the two ROOT files (dipion)
f11 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj11.root")
f13 = ROOT.TFile.Open("/volatile/clas12/thayward/rgk_dc_study/dipion/dipion_cj13.root")
t11 = f11.Get("PhysicsEvents")
t13 = f13.Get("PhysicsEvents")

_ratio_id = 0
def make_ratio(hn, hd, nb, xmin, xmax, xlabel):
    """Unique‐named ratio with Poisson errors."""
    global _ratio_id
    rname = f"r{_ratio_id}"
    _ratio_id += 1
    hr = ROOT.TH1F(rname, f"Ratio; {xlabel}; ratio", nb, xmin, xmax)
    for i in range(1, nb+1):
        n = hn.GetBinContent(i)
        d = hd.GetBinContent(i)
        if d>0:
            e_n, e_d = math.sqrt(n), math.sqrt(d)
            val = n/d
            err = math.sqrt((e_n/d)**2 + (n*e_d/(d*d))**2)
        else:
            val, err = 0.0, 0.0
        hr.SetBinContent(i, val)
        hr.SetBinError(i, err)
    hr.SetLineColor(ROOT.kBlue)
    hr.SetLineWidth(2)
    hr.GetYaxis().SetRangeUser(0.8, 1.2)
    return hr

def make_canvas(p_branch, th_branch, plabel, tlabel,
                cname, outfile, pmin, pmax, tmin, tmax):
    p_bins, t_bins = 70, 70

    # create histos
    hp11 = ROOT.TH1F("hp11", f"{plabel}; {plabel} (GeV); counts", p_bins, pmin, pmax)
    hp13 = ROOT.TH1F("hp13", f"{plabel}; {plabel} (GeV); counts", p_bins, pmin, pmax)
    ht11 = ROOT.TH1F("ht11", f"{tlabel}; {tlabel} (deg); counts", t_bins, tmin, tmax)
    ht13 = ROOT.TH1F("ht13", f"{tlabel}; {tlabel} (deg); counts", t_bins, tmin, tmax)

    # fill
    t11.Draw(f"{p_branch}>>hp11")
    t13.Draw(f"{p_branch}>>hp13")
    t11.Draw(f"{th_branch}*180./TMath::Pi()>>ht11")
    t13.Draw(f"{th_branch}*180./TMath::Pi()>>ht13")

    # style
    for h in (hp11, ht11):
        h.SetLineColor(ROOT.kBlue); h.SetLineWidth(2)
    for h in (hp13, ht13):
        h.SetLineColor(ROOT.kRed);  h.SetLineWidth(2)

    # integrals & maxima
    tot11_p = int(hp11.Integral()); tot13_p = int(hp13.Integral())
    tot11_t = int(ht11.Integral()); tot13_t = int(ht13.Integral())
    max_p = max(hp11.GetMaximum(), hp13.GetMaximum())
    max_t = max(ht11.GetMaximum(), ht13.GetMaximum())

    # canvas
    c = ROOT.TCanvas(cname, cname, 1200, 1000)
    c.Divide(2,2)

    # pad1: momentum histos
    c.cd(1); ROOT.gPad.SetLeftMargin(0.15)
    hp11.Draw("HIST"); hp13.Draw("HIST SAME")
    hp11.GetYaxis().SetRangeUser(0, 1.2*max_p)
    leg = ROOT.TLegend(0.6,0.65,0.88,0.88)
    leg.AddEntry(hp11, f"cj11.2.0 (N={tot11_p})", "l")
    leg.AddEntry(hp13, f"cj13.0.3 (N={tot13_p})", "l")
    leg.Draw()

    # pad2: momentum ratio
    c.cd(2); ROOT.gPad.SetLeftMargin(0.15)
    r_p = make_ratio(hp13, hp11, p_bins, pmin, pmax, f"{plabel} (GeV)")
    r_p.Draw("E1")

    # pad3: theta histos
    c.cd(3); ROOT.gPad.SetLeftMargin(0.15)
    ht11.Draw("HIST"); ht13.Draw("HIST SAME")
    ht11.GetYaxis().SetRangeUser(0, 1.2*max_t)
    leg2 = ROOT.TLegend(0.6,0.65,0.88,0.88)
    leg2.AddEntry(ht11, f"cj11.2.0 (N={tot11_t})", "l")
    leg2.AddEntry(ht13, f"cj13.0.3 (N={tot13_t})", "l")
    leg2.Draw()

    # pad4: theta ratio
    c.cd(4); ROOT.gPad.SetLeftMargin(0.15)
    r_t = make_ratio(ht13, ht11, t_bins, tmin, tmax, f"{tlabel} (deg)")
    r_t.Draw("E1")

    c.SaveAs(outfile)

# electron
make_canvas(
  "e_p","e_theta",
  "e_{p}","#theta_{e}",
  "c_electron","/u/home/thayward/ep_comparison.pdf",
  1.0,4.5, 5.0,35.0
)

# pi+
make_canvas(
  "p1_p","p1_theta",
  "#pi^{+}_{p}","#theta_{#pi^{+}}",
  "c_pi_plus","/u/home/thayward/p1_comparison.pdf",
  0.5,3.5, 0.0,60.0
)

# pi-
make_canvas(
  "p2_p","p2_theta",
  "#pi^{-}_{p}","#theta_{#pi^{-}}",
  "c_pi_minus","/u/home/thayward/p2_comparison.pdf",
  0.5,3.5, 0.0,60.0
)