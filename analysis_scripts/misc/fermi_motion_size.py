import ROOT
import math
import os
import sys

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# event limit
# ------------------------------------------------

max_events = 100000

if len(sys.argv) > 1:
    max_events = int(sys.argv[1])
#endif

print("Processing up to",max_events,"events")

# ------------------------------------------------
# input
# ------------------------------------------------

input_file = "/volatile/clas12/users/mkerr/rgc/p_pim/ND3/fall22_nd3_p_pim.root"
tree_name  = "PhysicsEvents"

m_e  = 0.000511
m_pi = 0.13957
M_p  = 0.938272
M_n  = 0.939565

# ------------------------------------------------
# beam energy from run number
# ------------------------------------------------

def get_beam_energy(run):

    Eb = 10.6041

    if 5032 <= run <= 5666: Eb = 10.6041
    elif 2365 <= run <= 2598: Eb = 2.22193
    elif 3030 <= run <= 3106: Eb = 6.42313
    elif 3819 <= run <= 3862: Eb = 6.42313
    elif 3172 <= run <= 3817: Eb = 10.5940
    elif 3863 <= run <= 4326: Eb = 10.5940
    elif 5875 <= run <= 6000: Eb = 6.535
    elif 5674 <= run <= 5870: Eb = 7.546
    elif 6616 <= run <= 6783: Eb = 10.1998
    elif 6120 <= run <= 6399: Eb = 10.5986
    elif 6409 <= run <= 6604: Eb = 10.1998
    elif 11093 <= run <= 11283: Eb = 10.4096
    elif 11284 <= run <= 11300: Eb = 4.17179
    elif 11323 <= run <= 11571: Eb = 10.3894
    elif 16042 <= run <= 17065: Eb = 10.5473
    elif 17067 <= run <= 17716: Eb = 10.5563
    elif 17717 <= run <= 17811: Eb = 10.5593
    elif 19204 <= run <= 19659: Eb = 6.39463
    elif 19662 <= run <= 19893: Eb = 8.47757

    return Eb

# ------------------------------------------------
# open ROOT file
# ------------------------------------------------

f = ROOT.TFile.Open(input_file)
t = f.Get(tree_name)

# ------------------------------------------------
# histograms
# ------------------------------------------------

nb_pf = 60
nb_pp = 60

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.5
)

h_dxB = ROOT.TH2D(
"h_dxB",
"|#Delta x_{B}|; |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.5
)

h_dt = ROOT.TH2D(
"h_dt",
"|#Delta t| (GeV^{2}); |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.5
)

h_dmx2 = ROOT.TH2D(
"h_dmx2",
"|#Delta Mx^{2}| (GeV^{2}); |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.5
)

# ------------------------------------------------
# scan events
# ------------------------------------------------

print("Scanning events")

sum_q = 0
sum_nu = 0
sum_K = 0
sum_xB = 0
count = 0

event_counter = 0

for ev in t:

    event_counter += 1
    if event_counter > max_events:
        break
    #endif

    run = ev.runnum
    Eb = get_beam_energy(run)

    # electron
    e_p = ev.e_p
    th  = ev.e_theta
    ph  = ev.e_phi

    ex = e_p*math.sin(th)*math.cos(ph)
    ey = e_p*math.sin(th)*math.sin(ph)
    ez = e_p*math.cos(th)

    Ee = math.sqrt(e_p*e_p + m_e*m_e)

    # proton
    p1 = ev.p1_p
    th = ev.p1_theta
    ph = ev.p1_phi

    ppx = p1*math.sin(th)*math.cos(ph)
    ppy = p1*math.sin(th)*math.sin(ph)
    ppz = p1*math.cos(th)

    Ep = math.sqrt(p1*p1 + M_p*M_p)

    # pion
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
    qz = Eb - ez
    qE = Eb - Ee

    qmag = math.sqrt(qx*qx+qy*qy+qz*qz)
    nu = qE

    Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz)

    xB = Q2/(2*M_n*nu)

    # missing mass
    mx2 = (qE + M_p - Epi)**2 - ((qx-pix)**2+(qy-piy)**2+(qz-piz)**2)

    if mx2 >= 1.07:
        continue
    #endif

    # neutron momentum
    pnx = ppx + pix - qx
    pny = ppy + piy - qy
    pnz = ppz + piz - qz

    pF = math.sqrt(pnx*pnx + pny*pny + pnz*pnz)

    if pF < 0.3 and 0.3 < p1 < 1.5:
        h_density.Fill(pF,p1)
    #endif

    Kx = qx - pix
    Ky = qy - piy
    Kz = qz - piz

    Kmag = math.sqrt(Kx*Kx+Ky*Ky+Kz*Kz)

    sum_q += qmag
    sum_nu += nu
    sum_K += Kmag
    sum_xB += xB

    count += 1

#endfor

print("Events used:",count)

qavg = sum_q/count
nuavg = sum_nu/count
Kavg = sum_K/count
xBavg = sum_xB/count

print("Average |q| =",qavg)
print("Average |q-p_pi| =",Kavg)

# ------------------------------------------------
# analytic shift maps
# ------------------------------------------------

for i in range(nb_pf):

    pF = (i+0.5)*0.3/nb_pf

    for j in range(nb_pp):

        pp = 0.3 + (j+0.5)*(1.2/nb_pp)

        dxB = xBavg*(pF*qavg)/(math.sqrt(3)*M_n*nuavg)

        dt = (2/math.sqrt(3))*pF*pp

        dmx2 = (2/math.sqrt(3))*pF*Kavg

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

c.SaveAs("output/fermi_shift_maps.png")

print("Saved output/fermi_shift_maps.png")