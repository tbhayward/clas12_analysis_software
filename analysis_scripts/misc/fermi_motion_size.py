import ROOT
import math
import random
import os

ROOT.gROOT.SetBatch(True)

input_file = "/scratch/thayward/fermi_check.root"
tree_name = "PhysicsEvents"

beamE = 10.55

m_e  = 0.000511
m_pi = 0.13957
M_p  = 0.938272
M_n  = 0.939565

# ------------------------------------------------
# PANEL 1 : DATA DENSITY
# ------------------------------------------------

f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

print("Building phase-space density")

for ev in t:

    e_p = ev.e_p
    th  = ev.e_theta
    ph  = ev.e_phi

    ex = e_p*math.sin(th)*math.cos(ph)
    ey = e_p*math.sin(th)*math.sin(ph)
    ez = e_p*math.cos(th)

    Ee = math.sqrt(e_p*e_p + m_e*m_e)

    p1 = ev.p1_p
    th = ev.p1_theta
    ph = ev.p1_phi

    pix = p1*math.sin(th)*math.cos(ph)
    piy = p1*math.sin(th)*math.sin(ph)
    piz = p1*math.cos(th)

    Epi = math.sqrt(p1*p1 + m_pi*m_pi)

    p2 = ev.p2_p
    th = ev.p2_theta
    ph = ev.p2_phi

    ppx = p2*math.sin(th)*math.cos(ph)
    ppy = p2*math.sin(th)*math.sin(ph)
    ppz = p2*math.cos(th)

    qx = -ex
    qy = -ey
    qz = beamE - ez
    qE = beamE - Ee

    mx2 = (qE + M_p - Epi)**2 - ((qx - pix)**2 + (qy - piy)**2 + (qz - piz)**2)

    if mx2 >= 1.07:
        continue
    #endif

    pFx = ppx + pix - qx
    pFy = ppy + piy - qy
    pFz = ppz + piz - qz

    pF = math.sqrt(pFx*pFx + pFy*pFy + pFz*pFz)

    if pF < 0.3 and 0.3 < p2 < 1.5:
        h_density.Fill(pF,p2)
    #endif

#endfor


# ------------------------------------------------
# PANELS 2–4 : THEORY MAPS
# ------------------------------------------------

h_dxB = ROOT.TH2D("|#Delta x_{B}|","|#Delta x_{B}|; |p_{F}| (GeV); p'_{p} (GeV)",60,0,0.3,60,0.3,1.5)
h_dt  = ROOT.TH2D("|#Delta t|","|#Delta t|; |p_{F}| (GeV); p'_{p} (GeV)",60,0,0.3,60,0.3,1.5)
h_dmx2= ROOT.TH2D("|#Delta Mx2|","|#Delta Mx^{2}|; |p_{F}| (GeV); p'_{p} (GeV)",60,0,0.3,60,0.3,1.5)

samples = 200

# typical DIS kinematics
Q2 = 2.0
nu = 5.0
xB = Q2/(2*M_n*nu)
q_mag = math.sqrt(Q2 + nu*nu)

print("Scanning theoretical phase space")

for i in range(60):

    pF = (i+0.5)*0.3/60

    for j in range(60):

        pp = 0.3 + (j+0.5)*(1.2/60)

        dxB = 0
        dt  = 0
        dmx2= 0

        for s in range(samples):

            costh = random.uniform(-1,1)

            dxB  += abs(xB*(pF*q_mag*costh)/(M_n*nu))
            dt   += abs(2*pF*pp*costh)
            dmx2 += abs(2*pF*pp*costh)

        #endfor

        dxB  /= samples
        dt   /= samples
        dmx2 /= samples

        h_dxB.SetBinContent(i+1,j+1,dxB)
        h_dt.SetBinContent(i+1,j+1,dt)
        h_dmx2.SetBinContent(i+1,j+1,dmx2)

#endfor


# ------------------------------------------------
# PLOT
# ------------------------------------------------

ROOT.gStyle.SetOptStat(0)

c = ROOT.TCanvas("c","c",2000,500)
c.Divide(4,1)

c.cd(1)
h_density.Draw("COLZ")

c.cd(2)
h_dxB.Draw("COLZ")

c.cd(3)
h_dt.Draw("COLZ")

c.cd(4)
h_dmx2.Draw("COLZ")

os.makedirs("output",exist_ok=True)

c.SaveAs("output/fermi_motion_bias_map.png")