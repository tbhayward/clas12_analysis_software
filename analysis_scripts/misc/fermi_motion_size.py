import ROOT
import math
import random
import os

ROOT.gROOT.SetBatch(True)

input_file = "/volatile/clas12/users/mkerr/rgc/p_pim/ND3/fall22_nd3_p_pim.root"
tree_name = "PhysicsEvents"

beamE = 10.55

m_e  = 0.000511
m_pi = 0.13957
M_p  = 0.938272
M_n  = 0.939565

f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

# -------------------------------------------
# PANEL 1 : EVENT DENSITY
# -------------------------------------------

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

# -------------------------------------------
# PANELS 2–4 : THEORY MAPS
# -------------------------------------------

h_dxB = ROOT.TH2D(
"h_dxB",
"|#Delta x_{B}|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_dt = ROOT.TH2D(
"h_dt",
"|#Delta t|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_dmx2 = ROOT.TH2D(
"h_dmx2",
"|#Delta Mx^{2}|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

kin_sample = []

print("Scanning events")

for ev in t:

    # electron
    e_p = ev.e_p
    th  = ev.e_theta
    ph  = ev.e_phi

    ex = e_p*math.sin(th)*math.cos(ph)
    ey = e_p*math.sin(th)*math.sin(ph)
    ez = e_p*math.cos(th)

    Ee = math.sqrt(e_p*e_p + m_e*m_e)

    # proton (p1)
    p1 = ev.p1_p
    th = ev.p1_theta
    ph = ev.p1_phi

    ppx = p1*math.sin(th)*math.cos(ph)
    ppy = p1*math.sin(th)*math.sin(ph)
    ppz = p1*math.cos(th)

    Ep = math.sqrt(p1*p1 + M_p*M_p)

    # pion (p2)
    p2 = ev.p2_p
    th = ev.p2_theta
    ph = ev.p2_phi

    pix = p2*math.sin(th)*math.cos(ph)
    piy = p2*math.sin(th)*math.sin(ph)
    piz = p2*math.cos(th)

    Epi = math.sqrt(p2*p2 + m_pi*m_pi)

    # virtual photon
    qx = -ex
    qy = -ey
    qz = beamE - ez
    qE = beamE - Ee

    # mesonic missing mass
    mx2_mes = (qE + M_p - Epi)**2 - ((qx - pix)**2 + (qy - piy)**2 + (qz - piz)**2)

    if mx2_mes >= 1.07:
        continue

    # neutron momentum from conservation
    pnx = ppx + pix - qx
    pny = ppy + piy - qy
    pnz = ppz + piz - qz

    pF = math.sqrt(pnx*pnx + pny*pny + pnz*pnz)

    if pF < 0.3 and 0.3 < p1 < 1.5:
        h_density.Fill(pF,p1)

    kin_sample.append((qx,qy,qz,qE,ppx,ppy,ppz,pix,piy,piz))

# -------------------------------------------
# THEORY SHIFT MAPS
# -------------------------------------------

print("Building shift maps")

samples = 200

for i in range(60):

    pF = (i+0.5)*0.3/60

    for j in range(60):

        pp = 0.3 + (j+0.5)*(1.2/60)

        dxB = 0
        dt  = 0
        dmx2 = 0

        for s in range(samples):

            qx,qy,qz,qE,ppx,ppy,ppz,pix,piy,piz = random.choice(kin_sample)

            costh = random.uniform(-1,1)

            # random Fermi direction
            pf_vec = (
                pF*costh,
                pF*math.sqrt(1-costh**2),
                0
            )

            pnx,pny,pnz = pf_vec

            En = math.sqrt(pF*pF + M_n*M_n)

            Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz)
            nu = qE

            xB_mes = Q2/(2*M_n*nu)
            xB_true = Q2/(2*(En*nu - (pnx*qx+pny*qy+pnz*qz)))

            dxB += abs(xB_mes-xB_true)

            t_mes = (qE-math.sqrt(pix*pix+piy*piy+piz*piz+m_pi*m_pi))**2 \
                    - ((qx-pix)**2 + (qy-piy)**2 + (qz-piz)**2)

            Ep = math.sqrt(pp*pp + M_p*M_p)

            t_true = (Ep-En)**2 - ((ppx-pnx)**2+(ppy-pny)**2+(ppz-pnz)**2)

            dt += abs(t_mes-t_true)

            mx2_true = (qE+En-math.sqrt(pix*pix+piy*piy+piz*piz+m_pi*m_pi))**2 \
                       - ((qx+pnx-pix)**2+(qy+pny-piy)**2+(qz+pnz-piz)**2)

            dmx2 += abs(mx2_mes-mx2_true)

        dxB /= samples
        dt  /= samples
        dmx2 /= samples

        h_dxB.SetBinContent(i+1,j+1,dxB)
        h_dt.SetBinContent(i+1,j+1,dt)
        h_dmx2.SetBinContent(i+1,j+1,dmx2)

# -------------------------------------------
# DRAW
# -------------------------------------------

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

c.SaveAs("output/fermi_bias_diagnostics.png")