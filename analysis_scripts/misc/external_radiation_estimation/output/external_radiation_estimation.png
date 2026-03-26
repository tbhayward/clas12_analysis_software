import ROOT
import os
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# ------------------------------------------------
# input/output
# ------------------------------------------------

input_file = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_fa18_inb_epi+.root"
tree_name  = "PhysicsEvents"
output_dir = "output"
output_png = os.path.join(output_dir, "external_radiation_estimation.png")

os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------
# physics constants
# ------------------------------------------------

M_n = 0.939565
M_n2 = M_n * M_n

# ------------------------------------------------
# binning
# ------------------------------------------------

vz_min = -7.0
vz_max =  2.0
vz_step = 1.0

mx2_min = 0.5
mx2_max = 1.1
mx2_bins = 120

vz_edges = []
current_edge = vz_min
while current_edge <= vz_max + 1.0e-9:
    vz_edges.append(current_edge)
    current_edge += vz_step
#endfor

n_vz_bins = len(vz_edges) - 1

# ------------------------------------------------
# open file and tree
# ------------------------------------------------

f = ROOT.TFile.Open(input_file, "READ")
if not f or f.IsZombie():
    raise RuntimeError("Could not open file: " + input_file)
#endif

tree = f.Get(tree_name)
if not tree:
    raise RuntimeError("Could not find tree: " + tree_name)
#endif

# ------------------------------------------------
# book histograms and accumulators
# ------------------------------------------------

histograms = []
sum_vz = [0.0] * n_vz_bins
count_vz = [0] * n_vz_bins

for i in range(n_vz_bins):
    h = ROOT.TH1D(
        "h_mx2_vzbin_%d" % i,
        "",
        mx2_bins,
        mx2_min,
        mx2_max
    )
    h.Sumw2()
    histograms.append(h)
#endfor

# ------------------------------------------------
# event loop
# ------------------------------------------------

n_entries = tree.GetEntries()

for i_entry in range(n_entries):
    tree.GetEntry(i_entry)

    vz_e = float(tree.vz_e)
    mx2  = float(tree.Mx2)

    if vz_e < vz_min or vz_e >= vz_max:
        continue
    #endif

    if mx2 < mx2_min or mx2 > mx2_max:
        continue
    #endif

    vz_bin = int((vz_e - vz_min) / vz_step)

    if vz_bin < 0 or vz_bin >= n_vz_bins:
        continue
    #endif

    histograms[vz_bin].Fill(mx2)
    sum_vz[vz_bin] += vz_e
    count_vz[vz_bin] += 1
#endfor

# ------------------------------------------------
# normalize histograms
# ------------------------------------------------

for i in range(n_vz_bins):
    integral = histograms[i].Integral()
    if integral > 0.0:
        histograms[i].Scale(1.0 / integral)
    #endif
#endfor

# ------------------------------------------------
# fit each histogram with gaus + polynomial background
# y-error on the right plot will be the fitted sigma
# ------------------------------------------------

x_vals = []
y_vals = []
y_errs = []

fit_functions = []

color_list = [
    ROOT.kBlack,
    ROOT.kRed + 1,
    ROOT.kBlue + 1,
    ROOT.kGreen + 2,
    ROOT.kMagenta + 1,
    ROOT.kCyan + 2,
    ROOT.kOrange + 7,
    ROOT.kViolet + 1,
    ROOT.kTeal + 3
]

for i in range(n_vz_bins):
    h = histograms[i]

    if count_vz[i] < 25:
        continue
    #endif

    if h.GetEntries() < 10:
        continue
    #endif

    mean_vz = sum_vz[i] / float(count_vz[i])

    max_bin = h.GetMaximumBin()
    peak_height = h.GetBinContent(max_bin)

    fit_name = "fit_mx2_vzbin_%d" % i
    fit_func = ROOT.TF1(fit_name, "gaus(0)+pol2(3)", mx2_min, mx2_max)

    fit_func.SetParameter(0, peak_height)
    fit_func.SetParameter(1, M_n2)
    fit_func.SetParameter(2, 0.03)
    fit_func.SetParameter(3, 0.0)
    fit_func.SetParameter(4, 0.0)
    fit_func.SetParameter(5, 0.0)

    fit_func.SetParLimits(0, 0.0, 10.0)
    fit_func.SetParLimits(1, 0.80, 0.95)
    fit_func.SetParLimits(2, 0.005, 0.20)

    h.Fit(fit_func, "RQ0")

    mu = fit_func.GetParameter(1)
    sigma = abs(fit_func.GetParameter(2))

    x_vals.append(mean_vz)
    y_vals.append(mu)
    y_errs.append(sigma)

    h.SetLineColor(color_list[i % len(color_list)])
    h.SetLineWidth(2)

    fit_func.SetLineColor(h.GetLineColor())
    fit_func.SetLineStyle(2)
    fit_functions.append(fit_func)
#endfor

# ------------------------------------------------
# build graph: point is fitted peak position, error bar is sigma
# ------------------------------------------------

n_points = len(x_vals)

graph = ROOT.TGraphErrors(n_points)
for i in range(n_points):
    graph.SetPoint(i, x_vals[i], y_vals[i])
    graph.SetPointError(i, 0.0, y_errs[i])
#endfor

graph.SetTitle("")
graph.SetMarkerStyle(20)
graph.SetMarkerSize(1.2)
graph.SetLineWidth(2)

# ------------------------------------------------
# canvas
# ------------------------------------------------

canvas = ROOT.TCanvas("canvas", "external radiation estimation", 1600, 700)
canvas.Divide(2, 1)

# ------------------------------------------------
# left pad: normalized Mx2 distributions in vz_e bins
# ------------------------------------------------

canvas.cd(1)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetRightMargin(0.04)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetTopMargin(0.08)

frame_left = ROOT.TH1D("frame_left", "", 100, mx2_min, mx2_max)
frame_left.SetMinimum(0.0)

max_y = 0.0
for h in histograms:
    if h.GetMaximum() > max_y:
        max_y = h.GetMaximum()
    #endif
#endfor

frame_left.SetMaximum(1.25 * max_y if max_y > 0.0 else 1.0)
frame_left.GetXaxis().SetTitle("M_{X}^{2} (GeV^{2})")
frame_left.GetYaxis().SetTitle("Normalized counts")
frame_left.GetXaxis().CenterTitle()
frame_left.GetYaxis().CenterTitle()
frame_left.Draw()

legend = ROOT.TLegend(0.52, 0.45, 0.93, 0.92)
legend.SetBorderSize(1)
legend.SetFillStyle(1001)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextSize(0.028)

first_drawn = False
for i in range(n_vz_bins):
    h = histograms[i]

    if count_vz[i] < 25:
        continue
    #endif

    draw_opt = "hist same"
    h.Draw(draw_opt)

    if i < len(fit_functions):
        fit_functions[i].Draw("same")
    #endif

    vz_lo = vz_edges[i]
    vz_hi = vz_edges[i + 1]
    legend.AddEntry(h, "%.1f < v_{z,e} < %.1f (cm)" % (vz_lo, vz_hi), "l")
    first_drawn = True
#endfor

legend.Draw()

latex_left = ROOT.TLatex()
latex_left.SetNDC()
latex_left.SetTextSize(0.040)
latex_left.DrawLatex(0.14, 0.93, "Exclusive #pi^{+}: normalized M_{X}^{2} distributions")

# ------------------------------------------------
# right pad: fitted peak position with sigma error bars vs mean vz_e
# ------------------------------------------------

canvas.cd(2)
ROOT.gPad.SetLeftMargin(0.13)
ROOT.gPad.SetRightMargin(0.04)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetTopMargin(0.08)

frame_right = ROOT.TH1D("frame_right", "", 100, vz_min, vz_max)
frame_right.SetMinimum(0.80)
frame_right.SetMaximum(0.95)
frame_right.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
frame_right.GetYaxis().SetTitle("Fitted neutron peak position M_{X}^{2} (GeV^{2})")
frame_right.GetXaxis().CenterTitle()
frame_right.GetYaxis().CenterTitle()
frame_right.Draw()

line_neutron = ROOT.TLine(vz_min, M_n2, vz_max, M_n2)
line_neutron.SetLineStyle(2)
line_neutron.SetLineWidth(2)
line_neutron.Draw("same")

graph.Draw("P SAME")

legend_right = ROOT.TLegend(0.17, 0.78, 0.62, 0.91)
legend_right.SetBorderSize(1)
legend_right.SetFillStyle(1001)
legend_right.SetFillColor(ROOT.kWhite)
legend_right.SetTextSize(0.032)
legend_right.AddEntry(graph, "Point = fitted mean, error bar = fitted #sigma", "lep")
legend_right.AddEntry(line_neutron, "Expected neutron mass squared", "l")
legend_right.Draw()

latex_right = ROOT.TLatex()
latex_right.SetNDC()
latex_right.SetTextSize(0.040)
latex_right.DrawLatex(0.14, 0.93, "Fitted peak position vs v_{z,e}")

# ------------------------------------------------
# save
# ------------------------------------------------

canvas.SaveAs(output_png)

print("Saved:", output_png)

f.Close()