import ROOT
import math
import random
import os

ROOT.gROOT.SetBatch(True)

# --------------------------------
# constants
# --------------------------------

M = 0.938272

beamE = 10.55
Q2 = 2.0
nu = 5.0

xB = Q2/(2*M*nu)

q_mag = math.sqrt(Q2 + nu*nu)

# --------------------------------
# scan ranges
# --------------------------------

pf_min = 0.0
pf_max = 0.3

pp_min = 0.3
pp_max = 1.5

nb_pf = 60
nb_pp = 60

samples = 300

# --------------------------------
# histograms
# --------------------------------

h_xB = ROOT.TH2D(
"h_xB",
"Shift in xB; p_{F} (GeV); p_{p} (GeV)",
nb_pf,pf_min,pf_max,
nb_pp,pp_min,pp_max
)

h_t = ROOT.TH2D(
"h_t",
"Shift in t; p_{F} (GeV); p_{p} (GeV)",
nb_pf,pf_min,pf_max,
nb_pp,pp_min,pp_max
)

h_mx2 = ROOT.TH2D(
"h_mx2",
"Shift in Mx2; p_{F} (GeV); p_{p} (GeV)",
nb_pf,pf_min,pf_max,
nb_pp,pp_min,pp_max
)

print("Scanning phase space")

# --------------------------------
# Monte Carlo sampling
# --------------------------------

for i_pf in range(nb_pf):

    pf = pf_min + (i_pf+0.5)*(pf_max-pf_min)/nb_pf

    for i_pp in range(nb_pp):

        pp = pp_min + (i_pp+0.5)*(pp_max-pp_min)/nb_pp

        sxB = 0
        st  = 0
        smx2 = 0

        for s in range(samples):

            costh = random.uniform(-1,1)

            dot = pf * pp * costh

            # xB shift
            delta_xB = abs(xB * (pf*q_mag*costh)/(M*nu))

            # t shift
            delta_t = abs(2*dot)

            # Mx2 shift
            delta_mx2 = abs(2*dot)

            sxB += delta_xB
            st += delta_t
            smx2 += delta_mx2

        #endfor

        sxB /= samples
        st  /= samples
        smx2 /= samples

        h_xB.SetBinContent(i_pf+1,i_pp+1,sxB)
        h_t.SetBinContent(i_pf+1,i_pp+1,st)
        h_mx2.SetBinContent(i_pf+1,i_pp+1,smx2)

    #endfor
#endfor

# --------------------------------
# determine common color scale
# --------------------------------

max_val = max(
    h_xB.GetMaximum(),
    h_t.GetMaximum(),
    h_mx2.GetMaximum()
)

for h in [h_xB,h_t,h_mx2]:
    h.SetMinimum(0)
    h.SetMaximum(max_val)

# --------------------------------
# plotting
# --------------------------------

ROOT.gStyle.SetOptStat(0)

c = ROOT.TCanvas("c","c",1800,600)
c.Divide(3,1)

c.cd(1)
h_xB.Draw("COLZ")

c.cd(2)
h_t.Draw("COLZ")

c.cd(3)
h_mx2.Draw("COLZ")

os.makedirs("output",exist_ok=True)

c.SaveAs("output/fermi_motion_shift_map.png")

print("Saved output/fermi_motion_shift_map.png")