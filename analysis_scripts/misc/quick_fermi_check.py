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

h_t_mes = ROOT.TH1D("h_t_mes","Mesonic definition; t (GeV^2); Counts",200,-14,0.5)
h_t_bar = ROOT.TH1D("h_t_bar","Baryonic definition; t (GeV^2); Counts",200,-2,0.5)
h_t_diff = ROOT.TH1D("h_t_diff","t_{#pi} - t_{p}; #Delta t (GeV^2); Counts",200,-14,14)

print("Starting event loop")

counter = 0

for ev in t:

    if ev.Mx2 >= 1.07:
        continue
    #endif

    counter += 1

    # electron
    e_p = ev.e_p
    e_theta = ev.e_theta
    e_phi = ev.e_phi

    ex = e_p*math.sin(e_theta)*math.cos(e_phi)
    ey = e_p*math.sin(e_theta)*math.sin(e_phi)
    ez = e_p*math.cos(e_theta)

    E_e = math.sqrt(e_p*e_p + m_e*m_e)

    # pion (p1)
    p1_p = ev.p1_p
    p1_theta = ev.p1_theta
    p1_phi = ev.p1_phi

    px_pi = p1_p*math.sin(p1_theta)*math.cos(p1_phi)
    py_pi = p1_p*math.sin(p1_theta)*math.sin(p1_phi)
    pz_pi = p1_p*math.cos(p1_theta)

    E_pi = math.sqrt(p1_p*p1_p + m_pi*m_pi)

    # virtual photon q = k - k'
    qE = beamE - E_e
    qx = -ex
    qy = -ey
    qz = beamE - ez

    # mesonic definition t = (q - p_pi)^2
    dE = qE - E_pi
    dx = qx - px_pi
    dy = qy - py_pi
    dz = qz - pz_pi

    t_mes = dE*dE - (dx*dx + dy*dy + dz*dz)

    # baryonic definition t = (p_n - p_p)^2
    p2_p = ev.p2_p
    E_p = math.sqrt(p2_p*p2_p + M_p*M_p)

    t_bar = M_n*M_n + M_p*M_p - 2.0*M_n*E_p

    t_diff = t_mes - t_bar

    if counter <= 20:
        print("Event",counter)
        print("electron p,theta,phi:",e_p,e_theta,e_phi)
        print("pion p,theta,phi:",p1_p,p1_theta,p1_phi)
        print("proton p:",p2_p)
        print("t_mes =",t_mes)
        print("t_bar =",t_bar)
        print("t_mes - t_bar =",t_diff)
        print("")
    #endif

    if not math.isnan(t_mes):
        h_t_mes.Fill(t_mes)
    #endif

    if not math.isnan(t_bar):
        h_t_bar.Fill(t_bar)
    #endif

    if not math.isnan(t_diff):
        h_t_diff.Fill(t_diff)
    #endif

#endfor

print("Total events passing cut:",counter)

c = ROOT.TCanvas("c","c",1800,500)
c.Divide(3,1)

c.cd(1)
h_t_mes.Draw()

c.cd(2)
h_t_bar.Draw()

c.cd(3)
h_t_diff.Draw()

os.makedirs("output", exist_ok=True)

c.SaveAs("output/quick_fermi_motion_check.png")