import ROOT
import math
import os
import sys

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# runtime event limit
# ------------------------------------------------

max_events = 100000

if len(sys.argv) > 1:
    max_events = int(sys.argv[1])
#endif

print("Processing up to", max_events, "events")

# ------------------------------------------------
# input file
# ------------------------------------------------

input_file = "/volatile/clas12/users/mkerr/rgc/p_pim/ND3/fall22_nd3_p_pim.root"
tree_name = "PhysicsEvents"

beamE = 10.55

m_e  = 0.000511
m_pi = 0.13957
M_p  = 0.938272
M_n  = 0.939565

f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

# ------------------------------------------------
# enable only needed branches
# ------------------------------------------------

t.SetBranchStatus("*",0)

branches = [
"e_p","e_theta","e_phi",
"p1_p","p1_theta","p1_phi",
"p2_p","p2_theta","p2_phi"
]

for b in branches:
    t.SetBranchStatus(b,1)
#endfor

# ------------------------------------------------
# histograms
# ------------------------------------------------

nbins_pf = 60
nbins_pp = 60

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p'_{p} (GeV)",
nbins_pf,0,0.3,
nbins_pp,0.3,1.5
)

h_dxB = ROOT.TH2D(
"h_dxB",
"|#Delta x_{B}|; |p_{F}| (GeV); p'_{p} (GeV)",
nbins_pf,0,0.3,
nbins_pp,0.3,1.5
)

h_dt = ROOT.TH2D(
"h_dt",
"|#Delta t| (GeV^{2}); |p_{F}| (GeV); p'_{p} (GeV)",
nbins_pf,0,0.3,
nbins_pp,0.3,1.5
)

h_dmx2 = ROOT.TH2D(
"h_dmx2",
"|#Delta Mx^{2}| (GeV^{2}); |p_{F}| (GeV); p'_{p} (GeV)",
nbins_pf,0,0.3,
nbins_pp,0.3,1.5
)

# ------------------------------------------------
# scan events
# ------------------------------------------------

print("Scanning events")

sum_K = 0
sum_q = 0
sum_nu = 0
count = 0

event_counter = 0

for ev in t:

    event_counter += 1
    if event_counter > max_events:
        break
    #endif

    # electron
    e_p = ev.e_p
    th = ev.e_theta
    ph = ev.e_phi

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

    nu = qE
    q_mag = math.sqrt(qx*qx + qy*qy + qz*qz)

    # mesonic missing mass
    mx2_mes = (qE + M_p - Epi)**2 - ((qx - pix)**2 + (qy - piy)**2 + (qz - piz)**2)

    if mx2_mes >= 1.07:
        continue
    #endif

    # neutron momentum estimate
    pnx = ppx + pix - qx
    pny = ppy + piy - qy
    pnz = ppz + piz - qz

    pF = math.sqrt(pnx*pnx + pny*pny + pnz*pnz)

    if pF < 0.3 and 0.3 < p1 < 1.5:
        h_density.Fill(pF,p1)
    #endif

    # K vector
    Kx = qx - pix
    Ky = qy - piy
    Kz = qz - piz

    Kmag = math.sqrt(Kx*Kx + Ky*Ky + Kz*Kz)

    sum_K += Kmag
    sum_q += q_mag
    sum_nu += nu

    count += 1

#endfor

print("Events used:",count)

# ------------------------------------------------
# average kinematics
# ------------------------------------------------

Kavg = sum_K/count
qavg = sum_q/count
nuavg = sum_nu/count

print("Average |q-p_pi| =",Kavg)
print("Average |q| =",qavg)

# ------------------------------------------------
# compute analytic shifts
# ------------------------------------------------

for i in range(nbins_pf):

    pF = (i+0.5)*0.3/nbins_pf

    for j in range(nbins_pp):

        pp = 0.3 + (j+0.5)*(1.2/nbins_pp)

        # isotropic expectation
        proj = pF/math.sqrt(3)

        dt = 2*Kavg*proj
        dmx2 = 2*Kavg*proj

        dxB = (pF*qavg)/(M_n*nuavg)

        h_dxB.SetBinContent(i+1,j+1,abs(dxB))
        h_dt.SetBinContent(i+1,j+1,abs(dt))
        h_dmx2.SetBinContent(i+1,j+1,abs(dmx2))

#endfor

# ------------------------------------------------
# plotting
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

c.SaveAs("output/fermi_bias_diagnostics.png")

print("Saved output/fermi_bias_diagnostics.png")