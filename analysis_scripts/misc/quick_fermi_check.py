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

# no-cut histograms (black)
h_t_mes_all  = ROOT.TH1D("h_t_mes_all","Mesonic definition; t (GeV^2); Counts",200,-14,0.5)
h_t_bar_all  = ROOT.TH1D("h_t_bar_all","Baryonic definition; t (GeV^2); Counts",200,-2,0.5)
h_t_diff_all = ROOT.TH1D("h_t_diff_all","t_{#pi}-t_{p}; #Delta t (GeV^2); Counts",200,-14,14)

# cut histograms (red)
h_t_mes_cut  = ROOT.TH1D("h_t_mes_cut","Mesonic definition; t (GeV^2); Counts",200,-14,0.5)
h_t_bar_cut  = ROOT.TH1D("h_t_bar_cut","Baryonic definition; t (GeV^2); Counts",200,-2,0.5)
h_t_diff_cut = ROOT.TH1D("h_t_diff_cut","t_{#pi}-t_{p}; #Delta t (GeV^2); Counts",200,-14,14)

# styling
for h in [h_t_mes_all, h_t_bar_all, h_t_diff_all]:
    h.SetLineColor(ROOT.kBlack)

for h in [h_t_mes_cut, h_t_bar_cut, h_t_diff_cut]:
    h.SetLineColor(ROOT.kRed)

print("Starting event loop")

counter = 0

for ev in t:

    counter += 1

    # electron
    e_p = ev.e_p
    e_theta = ev.e_theta
    e_phi = ev.e_phi

    ex = e_p*math.sin(e_theta)*math.cos(e_phi)
    ey = e_p*math.sin(e_theta)*math.sin(e_phi)
    ez = e_p*math.cos(e_theta)

    E_e = math.sqrt(e_p*e_p + m_e*m_e)

    # pion
    p1_p = ev.p1_p
    p1_theta = ev.p1_theta
    p1_phi = ev.p1_phi

    px_pi = p1_p*math.sin(p1_theta)*math.cos(p1_phi)
    py_pi = p1_p*math.sin(p1_theta)*math.sin(p1_phi)
    pz_pi = p1_p*math.cos(p1_theta)

    E_pi = math.sqrt(p1_p*p1_p + m_pi*m_pi)

    # virtual photon
    qE = beamE - E_e
    qx = -ex
    qy = -ey
    qz = beamE - ez

    # mesonic t
    dE = qE - E_pi
    dx = qx - px_pi
    dy = qy - py_pi
    dz = qz - pz_pi

    t_mes = dE*dE - (dx*dx + dy*dy + dz*dz)

    # baryonic t
    p2_p = ev.p2_p
    E_p = math.sqrt(p2_p*p2_p + M_p*M_p)

    t_bar = M_n*M_n + M_p*M_p - 2.0*M_n*E_p

    t_diff = t_mes - t_bar

    # fill no-cut histograms
    if not math.isnan(t_mes):
        h_t_mes_all.Fill(t_mes)
    #endif

    if not math.isnan(t_bar):
        h_t_bar_all.Fill(t_bar)
    #endif

    if not math.isnan(t_diff):
        h_t_diff_all.Fill(t_diff)
    #endif

    # apply Mx2 cut
    if ev.Mx2 < 1.07:

        if not math.isnan(t_mes):
            h_t_mes_cut.Fill(t_mes)
        #endif

        if not math.isnan(t_bar):
            h_t_bar_cut.Fill(t_bar)
        #endif

        if not math.isnan(t_diff):
            h_t_diff_cut.Fill(t_diff)
        #endif

    #endif

#endfor

print("Total events processed:",counter)

c = ROOT.TCanvas("c","c",1800,500)
c.Divide(3,1)

# Mesonic
c.cd(1)
h_t_mes_all.Draw("hist")
h_t_mes_cut.Draw("hist same")

# Baryonic
c.cd(2)
h_t_bar_all.Draw("hist")
h_t_bar_cut.Draw("hist same")

# Difference
c.cd(3)
h_t_diff_all.Draw("hist")
h_t_diff_cut.Draw("hist same")

os.makedirs("output", exist_ok=True)

c.SaveAs("output/quick_fermi_motion_check.png")