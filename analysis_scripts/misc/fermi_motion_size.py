import ROOT
import math
import os

ROOT.gROOT.SetBatch(True)

input_file = "/scratch/thayward/fermi_check.root"
tree_name = "PhysicsEvents"

beamE = 10.55

m_e  = 0.000511
m_pi = 0.13957
M_p  = 0.938272
M_n  = 0.939565

f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

# ------------------------------------------
# histograms
# ------------------------------------------

h_xB = ROOT.TH2D(
"h_xB",
"Shift in xB; |p_{F}| (GeV); p_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_t = ROOT.TH2D(
"h_t",
"Shift in t; |p_{F}| (GeV); p_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_mx2 = ROOT.TH2D(
"h_mx2",
"Shift in Mx2; |p_{F}| (GeV); p_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

print("Starting event loop")

for ev in t:

    # --------------------------------
    # electron
    # --------------------------------

    e_p = ev.e_p
    th  = ev.e_theta
    ph  = ev.e_phi

    ex = e_p*math.sin(th)*math.cos(ph)
    ey = e_p*math.sin(th)*math.sin(ph)
    ez = e_p*math.cos(th)

    Ee = math.sqrt(e_p*e_p + m_e*m_e)

    # --------------------------------
    # pion
    # --------------------------------

    p1 = ev.p1_p
    th = ev.p1_theta
    ph = ev.p1_phi

    pix = p1*math.sin(th)*math.cos(ph)
    piy = p1*math.sin(th)*math.sin(ph)
    piz = p1*math.cos(th)

    Epi = math.sqrt(p1*p1 + m_pi*m_pi)

    # --------------------------------
    # proton
    # --------------------------------

    p2 = ev.p2_p
    th = ev.p2_theta
    ph = ev.p2_phi

    ppx = p2*math.sin(th)*math.cos(ph)
    ppy = p2*math.sin(th)*math.sin(ph)
    ppz = p2*math.cos(th)

    Ep = math.sqrt(p2*p2 + M_p*M_p)

    # --------------------------------
    # virtual photon
    # --------------------------------

    qx = -ex
    qy = -ey
    qz = beamE - ez
    qE = beamE - Ee

    # --------------------------------
    # missing mass cut
    # --------------------------------

    mx2 = (qE + M_p - Epi)**2 - ((qx - pix)**2 + (qy - piy)**2 + (qz - piz)**2)

    if mx2 >= 1.07:
        continue
    #endif

    # --------------------------------
    # estimate Fermi momentum
    # --------------------------------

    pFx = ppx + pix - qx
    pFy = ppy + piy - qy
    pFz = ppz + piz - qz

    pF = math.sqrt(pFx*pFx + pFy*pFy + pFz*pFz)

    # --------------------------------
    # observables
    # --------------------------------

    Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz)
    nu = qE

    xB = Q2/(2*M_n*nu)

    xB0 = Q2/(2*M_n*nu)

    dxB = abs(xB - xB0)

    # t definitions

    dE = qE - Epi
    dx = qx - pix
    dy = qy - piy
    dz = qz - piz

    t_mes = dE*dE - (dx*dx + dy*dy + dz*dz)

    t_bar = M_n*M_n + M_p*M_p - 2*M_n*Ep

    dt = abs(t_mes - t_bar)

    # Mx2 shift estimate

    dmx2 = abs(2*(pFx*ppx + pFy*ppy + pFz*ppz))

    # --------------------------------
    # fill histograms
    # --------------------------------

    if pF < 0.3 and 0.3 < p2 < 1.5:

        h_xB.Fill(pF,p2,dxB)
        h_t.Fill(pF,p2,dt)
        h_mx2.Fill(pF,p2,dmx2)

#endfor

# ------------------------------------------
# common color scale
# ------------------------------------------

max_val = max(h_xB.GetMaximum(),h_t.GetMaximum(),h_mx2.GetMaximum())

for h in [h_xB,h_t,h_mx2]:
    h.SetMinimum(0)
    h.SetMaximum(max_val)

# ------------------------------------------
# plotting
# ------------------------------------------

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

c.SaveAs("output/fermi_shift_data_map.png")

print("Saved output/fermi_shift_data_map.png")