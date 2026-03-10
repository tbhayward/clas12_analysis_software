import ROOT
import math
import random
import os

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# constants
# ------------------------------------------------

M_p = 0.938272
M_n = 0.939565

beamE = 10.55

# typical DIS scale for example
Q2 = 2.0
nu = 5.0

# ------------------------------------------------
# histograms
# ------------------------------------------------

nbins_pf = 60
nbins_pp = 60

pf_max = 0.4
pp_max = 3.0

h_xB = ROOT.TH2D(
    "h_xB",
    "Shift in xB; p_{F} (GeV); p_{p} (GeV)",
    nbins_pf,0,pf_max,
    nbins_pp,0,pp_max
)

h_t = ROOT.TH2D(
    "h_t",
    "Shift in t; p_{F} (GeV); p_{p} (GeV)",
    nbins_pf,0,pf_max,
    nbins_pp,0,pp_max
)

h_mx2 = ROOT.TH2D(
    "h_mx2",
    "Shift in Mx2; p_{F} (GeV); p_{p} (GeV)",
    nbins_pf,0,pf_max,
    nbins_pp,0,pp_max
)

# ------------------------------------------------
# sampling
# ------------------------------------------------

samples_per_bin = 200

print("Starting scan")

for i_pf in range(nbins_pf):

    pf = (i_pf + 0.5) * pf_max / nbins_pf

    for i_pp in range(nbins_pp):

        pp = (i_pp + 0.5) * pp_max / nbins_pp

        shift_xB = 0.0
        shift_t = 0.0
        shift_mx2 = 0.0

        for s in range(samples_per_bin):

            # random orientation
            costh = random.uniform(-1,1)

            # dot product term
            dot = pf * pp * costh

            # ------------------------------------------------
            # xB shift
            # ------------------------------------------------

            xB = Q2 / (2*M_n*nu)

            delta_xB = xB * (pf / M_n)

            shift_xB += abs(delta_xB)

            # ------------------------------------------------
            # t shift
            # ------------------------------------------------

            delta_t = 2.0 * dot

            shift_t += abs(delta_t)

            # ------------------------------------------------
            # Mx2 shift
            # ------------------------------------------------

            delta_mx2 = -2.0 * dot

            shift_mx2 += abs(delta_mx2)

        #endfor

        shift_xB /= samples_per_bin
        shift_t /= samples_per_bin
        shift_mx2 /= samples_per_bin

        h_xB.SetBinContent(i_pf+1,i_pp+1,shift_xB)
        h_t.SetBinContent(i_pf+1,i_pp+1,shift_t)
        h_mx2.SetBinContent(i_pf+1,i_pp+1,shift_mx2)

    #endfor
#endfor

# ------------------------------------------------
# plotting
# ------------------------------------------------

c = ROOT.TCanvas("c","c",1800,600)
c.Divide(3,1)

ROOT.gStyle.SetOptStat(0)

c.cd(1)
h_xB.Draw("COLZ")

c.cd(2)
h_t.Draw("COLZ")

c.cd(3)
h_mx2.Draw("COLZ")

os.makedirs("output",exist_ok=True)

c.SaveAs("output/fermi_motion_shift_map.png")

print("Saved output/fermi_motion_shift_map.png")