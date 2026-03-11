import ROOT
import math
import os
import sys
from array import array

ROOT.gROOT.SetBatch(True)

# ------------------------------------------------
# runtime event limit
# ------------------------------------------------

max_events = 100000000

if len(sys.argv) > 1:
    max_events = int(sys.argv[1])
#endif

print("Processing up to", max_events, "events")

# ------------------------------------------------
# constants
# ------------------------------------------------

m_e  = 0.000511
m_pi = 0.13957
M_p  = 0.938272
M_n  = 0.939565

HBARC = 0.1973269804

# ------------------------------------------------
# beam energy function
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
# nitrogen momentum table
# ------------------------------------------------

k_vals = [
0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,
1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9
]

N14_vals = [
0.163119,0.183235,0.225102,0.262783,0.276924,
0.263379,0.228534,0.183451,0.137812,0.097787,
0.066275,0.043066,0.026994,0.016454,0.009841,
0.005874,0.003586,0.002312,0.001622,0.001247
]

# ------------------------------------------------
# bin edges (0–0.1 fm⁻¹ etc)
# ------------------------------------------------

edges = []

for k in k_vals:
    edges.append(k * HBARC)
#endfor

edges.append((k_vals[-1] + 0.1) * HBARC)

edge_array = array('d', edges)

nb_pf_pdf = len(edges) - 1


# ------------------------------------------------
# input files
# ------------------------------------------------

# files = [
# "/volatile/clas12/users/mkerr/rgc/p_pim/ND3/fall22_nd3_p_pim.root",
# "/volatile/clas12/users/mkerr/rgc/p_pim/ND3/summer22_nd3_p_pim.root",
# "/volatile/clas12/users/mkerr/rgc/p_pim/ND3/spring23_nd3_p_pim.root"
# ]

files = [
"/volatile/clas12/users/mkerr/rgb/rgb_sidisdvcs_p_pim.root",
]

tree_name = "PhysicsEvents"


# ------------------------------------------------
# histogram binning
# ------------------------------------------------

nb_pf = 20
nb_pp = 20


# ------------------------------------------------
# heatmap histograms
# ------------------------------------------------

h_density = ROOT.TH2D(
"h_density",
"Event density; |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.2
)

h_dxB = ROOT.TH2D(
"h_dxB",
"<|#Delta x_{B}|>; |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.2
)

h_dt = ROOT.TH2D(
"h_dt",
"<|#Delta t|>; |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.2
)

h_dmx2 = ROOT.TH2D(
"h_dmx2",
"<|#Delta Mx^{2}|>; |p_{F}| (GeV); p' (GeV)",
nb_pf,0,0.3,
nb_pp,0.3,1.2
)


# ------------------------------------------------
# spectator momentum histograms
# ------------------------------------------------

h_pf = ROOT.TH1D(
"h_pf",
"Extracted |p_{F}| distribution; |p_{F}| (GeV); normalized PDF",
nb_pf_pdf,
edge_array
)

h_nitrogen = ROOT.TH1D(
"h_nitrogen",
"Nitrogen momentum distribution comparison; |p_{F}| (GeV); normalized PDF",
nb_pf_pdf,
edge_array
)

for i,val in enumerate(N14_vals):
    h_nitrogen.SetBinContent(i+1,val)
#endfor


# ------------------------------------------------
# accumulation arrays
# ------------------------------------------------

sum_dxB = [[0]*nb_pp for _ in range(nb_pf)]
sum_dt  = [[0]*nb_pp for _ in range(nb_pf)]
sum_dmx2= [[0]*nb_pp for _ in range(nb_pf)]
counts  = [[0]*nb_pp for _ in range(nb_pf)]


# ------------------------------------------------
# event loop
# ------------------------------------------------

print("Scanning events")

event_counter = 0

for file in files:

    print("Opening",file)

    f = ROOT.TFile.Open(file)
    t = f.Get(tree_name)

    for ev in t:

        event_counter += 1

        if event_counter > max_events:
            break
        #endif

        run = ev.runnum
        Eb = get_beam_energy(run)

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

        ppx = p1*math.sin(th)*math.cos(ph)
        ppy = p1*math.sin(th)*math.sin(ph)
        ppz = p1*math.cos(th)

        Ep = math.sqrt(p1*p1 + M_p*M_p)

        p2 = ev.p2_p
        th = ev.p2_theta
        ph = ev.p2_phi

        pix = p2*math.sin(th)*math.cos(ph)
        piy = p2*math.sin(th)*math.sin(ph)
        piz = p2*math.cos(th)

        Epi = math.sqrt(p2*p2 + m_pi*m_pi)

        qx = -ex
        qy = -ey
        qz = Eb - ez
        qE = Eb - Ee

        Q2 = -(qE*qE - qx*qx - qy*qy - qz*qz)
        nu = qE

        mx2_mes = (qE + M_p - Epi)**2 - ((qx-pix)**2+(qy-piy)**2+(qz-piz)**2)

        if mx2_mes >= 1.07:
            continue
        #endif

        pnx = ppx + pix - qx
        pny = ppy + piy - qy
        pnz = ppz + piz - qz

        pF = math.sqrt(pnx*pnx + pny*pny + pnz*pnz)

        h_pf.Fill(pF)

    #endfor

#endfor


# ------------------------------------------------
# normalize PDFs
# ------------------------------------------------

if h_pf.Integral() > 0:
    h_pf.Scale(1.0/h_pf.Integral())

if h_nitrogen.Integral() > 0:
    h_nitrogen.Scale(1.0/h_nitrogen.Integral())


# ------------------------------------------------
# PDF comparison plot
# ------------------------------------------------

c2 = ROOT.TCanvas("c2","pF comparison",800,600)

max_val = max(h_pf.GetMaximum(),h_nitrogen.GetMaximum())

h_pf.SetMaximum(1.2*max_val)

h_pf.SetLineColor(ROOT.kBlue)
h_pf.SetLineWidth(3)

h_nitrogen.SetLineColor(ROOT.kRed)
h_nitrogen.SetLineWidth(3)

h_pf.Draw("HIST")
h_nitrogen.Draw("HIST SAME")

mean_pf = h_pf.GetMean()
mean_nitrogen = h_nitrogen.GetMean()

legend = ROOT.TLegend(0.55,0.65,0.85,0.85)
legend.AddEntry(h_pf,f"Extracted p_F (mean = {mean_pf:.3f} GeV)","l")
legend.AddEntry(h_nitrogen,f"N14 distribution (mean = {mean_nitrogen:.3f} GeV)","l")
legend.Draw()

os.makedirs("output",exist_ok=True)

c2.SaveAs("output/fermi_pf_pdf_comparison.png")