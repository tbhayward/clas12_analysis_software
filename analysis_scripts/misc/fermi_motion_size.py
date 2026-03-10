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

f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

# -----------------------------------------
# histograms (DATA)
# -----------------------------------------

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_t_data = ROOT.TH2D(
"h_t_data",
"|t_{mes}-t_{bar}| (data); |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_mx2_data = ROOT.TH2D(
"h_mx2_data",
"|#Delta Mx^{2}| (data); |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

print("Running data loop")

for ev in t:

    # electron
    e_p = ev.e_p
    th  = ev.e_theta
    ph  = ev.e_phi

    ex = e_p*math.sin(th)*math.cos(ph)
    ey = e_p*math.sin(th)*math.sin(ph)
    ez = e_p*math.cos(th)

    Ee = math.sqrt(e_p*e_p + m_e*m_e)

    # pion
    p1 = ev.p1_p
    th = ev.p1_theta
    ph = ev.p1_phi

    pix = p1*math.sin(th)*math.cos(ph)
    piy = p1*math.sin(th)*math.sin(ph)
    piz = p1*math.cos(th)

    Epi = math.sqrt(p1*p1 + m_pi*m_pi)

    # proton
    p2 = ev.p2_p
    th = ev.p2_theta
    ph = ev.p2_phi

    ppx = p2*math.sin(th)*math.cos(ph)
    ppy = p2*math.sin(th)*math.sin(ph)
    ppz = p2*math.cos(th)

    Ep = math.sqrt(p2*p2 + M_p*M_p)

    # virtual photon
    qx = -ex
    qy = -ey
    qz = beamE - ez
    qE = beamE - Ee

    # exclusivity cut
    mx2 = (qE + M_p - Epi)**2 - ((qx - pix)**2 + (qy - piy)**2 + (qz - piz)**2)

    if mx2 >= 1.07:
        continue
    #endif

    # estimate Fermi momentum
    pFx = ppx + pix - qx
    pFy = ppy + piy - qy
    pFz = ppz + piz - qz

    pF = math.sqrt(pFx*pFx + pFy*pFy + pFz*pFz)

    if pF > 0.3 or p2 < 0.3 or p2 > 1.5:
        continue
    #endif

    h_density.Fill(pF,p2)

    # t definitions
    dE = qE - Epi
    dx = qx - pix
    dy = qy - piy
    dz = qz - piz

    t_mes = dE*dE - (dx*dx + dy*dy + dz*dz)

    t_bar = M_n*M_n + M_p*M_p - 2*M_n*Ep

    dt = abs(t_mes - t_bar)

    h_t_data.Fill(pF,p2,dt)

    # approximate Mx2 shift
    dmx2 = abs(2*(pFx*ppx + pFy*ppy + pFz*ppz))

    h_mx2_data.Fill(pF,p2,dmx2)

#endfor

# -----------------------------------------
# theoretical shift maps
# -----------------------------------------

h_blank = ROOT.TH2D(
"h_blank",
"Phase space reference; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_t_theory = ROOT.TH2D(
"h_t_theory",
"Expected |#Delta t|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

h_mx2_theory = ROOT.TH2D(
"h_mx2_theory",
"Expected |#Delta Mx^{2}|; |p_{F}| (GeV); p'_{p} (GeV)",
60,0,0.3,
60,0.3,1.5
)

samples = 200

print("Scanning theoretical phase space")

for i in range(60):

    pF = (i+0.5)*0.3/60

    for j in range(60):

        pp = 0.3 + (j+0.5)*(1.2/60)

        avg_t = 0
        avg_mx2 = 0

        for s in range(samples):

            costh = random.uniform(-1,1)

            avg_t += abs(2*pF*pp*costh)
            avg_mx2 += abs(2*pF*pp*costh)

        #endfor

        avg_t /= samples
        avg_mx2 /= samples

        h_t_theory.SetBinContent(i+1,j+1,avg_t)
        h_mx2_theory.SetBinContent(i+1,j+1,avg_mx2)

#endfor

# -----------------------------------------
# plotting
# -----------------------------------------

ROOT.gStyle.SetOptStat(0)

c = ROOT.TCanvas("c","c",1800,1000)
c.Divide(3,2)

c.cd(1)
h_density.Draw("COLZ")

c.cd(2)
h_t_data.Draw("COLZ")

c.cd(3)
h_mx2_data.Draw("COLZ")

c.cd(4)
h_blank.Draw("COLZ")

c.cd(5)
h_t_theory.Draw("COLZ")

c.cd(6)
h_mx2_theory.Draw("COLZ")

os.makedirs("output",exist_ok=True)

c.SaveAs("output/fermi_motion_diagnostics.png")

print("Saved output/fermi_motion_diagnostics.png")