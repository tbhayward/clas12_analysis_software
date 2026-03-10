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

# ----------------------------------------
# histograms
# ----------------------------------------

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_dxB = ROOT.TH2D(
"h_dxB",
"|x_{B}^{mes}-x_{B}^{true}|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_dt = ROOT.TH2D(
"h_dt",
"|t_{mes}-t_{true}|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_dmx2 = ROOT.TH2D(
"h_dmx2",
"|Mx2_{mes}-Mx2_{true}|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

print("Looping over events")

for ev in t:

    # ----------------------------------
    # electron
    # ----------------------------------

    e_p = ev.e_p
    th  = ev.e_theta
    ph  = ev.e_phi

    ex = e_p*math.sin(th)*math.cos(ph)
    ey = e_p*math.sin(th)*math.sin(ph)
    ez = e_p*math.cos(th)

    Ee = math.sqrt(e_p*e_p + m_e*m_e)

    # ----------------------------------
    # pion
    # ----------------------------------

    p1 = ev.p1_p
    th = ev.p1_theta
    ph = ev.p1_phi

    pix = p1*math.sin(th)*math.cos(ph)
    piy = p1*math.sin(th)*math.sin(ph)
    piz = p1*math.cos(th)

    Epi = math.sqrt(p1*p1 + m_pi*m_pi)

    # ----------------------------------
    # proton
    # ----------------------------------

    p2 = ev.p2_p
    th = ev.p2_theta
    ph = ev.p2_phi

    ppx = p2*math.sin(th)*math.cos(ph)
    ppy = p2*math.sin(th)*math.sin(ph)
    ppz = p2*math.cos(th)

    Ep = math.sqrt(p2*p2 + M_p*M_p)

    # ----------------------------------
    # virtual photon
    # ----------------------------------

    qx = -ex
    qy = -ey
    qz = beamE - ez
    qE = beamE - Ee

    # ----------------------------------
    # mesonic missing mass cut
    # ----------------------------------

    mx2_mes = (qE + M_p - Epi)**2 - ((qx - pix)**2 + (qy - piy)**2 + (qz - piz)**2)

    if mx2_mes >= 1.07:
        continue
    #endif

    # ----------------------------------
    # reconstructed neutron momentum
    # ----------------------------------

    pnx = ppx + pix - qx
    pny = ppy + piy - qy
    pnz = ppz + piz - qz

    pF = math.sqrt(pnx*pnx + pny*pny + pnz*pnz)

    if pF > 0.3 or p2 < 0.3 or p2 > 1.5:
        continue
    #endif

    En = math.sqrt(pF*pF + M_n*M_n)

    # ----------------------------------
    # Q2 and nu
    # ----------------------------------

    Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz)
    nu = qE

    # ----------------------------------
    # xB
    # ----------------------------------

    xB_mes = Q2/(2*M_n*nu)
    xB_true = Q2/(2*(En*nu - (pnx*qx+pny*qy+pnz*qz)))

    dxB = abs(xB_mes-xB_true)

    # ----------------------------------
    # t
    # ----------------------------------

    dE = qE - Epi
    dx = qx - pix
    dy = qy - piy
    dz = qz - piz

    t_mes = dE*dE - (dx*dx+dy*dy+dz*dz)

    t_true = (Ep-En)**2 - ((ppx-pnx)**2 + (ppy-pny)**2 + (ppz-pnz)**2)

    dt = abs(t_mes-t_true)

    # ----------------------------------
    # Mx2 true
    # ----------------------------------

    mx2_true = (qE + En - Epi)**2 - ((qx + pnx - pix)**2 +
                                    (qy + pny - piy)**2 +
                                    (qz + pnz - piz)**2)

    dmx2 = abs(mx2_mes-mx2_true)

    # ----------------------------------
    # fill histograms
    # ----------------------------------

    h_density.Fill(pF,p2)
    h_dxB.Fill(pF,p2,dxB)
    h_dt.Fill(pF,p2,dt)
    h_dmx2.Fill(pF,p2,dmx2)

#endfor

ROOT.gStyle.SetOptStat(0)

c = ROOT.TCanvas("c","c",1800,1000)
c.Divide(3,2)

c.cd(1)
h_density.Draw("COLZ")

c.cd(2)
h_density.Draw("COLZ")

c.cd(3)
h_density.Draw("COLZ")

c.cd(4)
h_dxB.Draw("COLZ")

c.cd(5)
h_dt.Draw("COLZ")

c.cd(6)
h_dmx2.Draw("COLZ")

os.makedirs("output",exist_ok=True)

c.SaveAs("output/fermi_motion_truth_comparison.png")