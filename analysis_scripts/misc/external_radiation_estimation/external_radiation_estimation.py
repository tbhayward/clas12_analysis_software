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
# remove the -7 to -6 bin by starting at -6
# ------------------------------------------------

vz_min = -6.0
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
fit_statuses = []

for i in range(n_vz_bins):
    h = histograms[i]

    if count_vz[i] < 25:
        fit_functions.append(None)
        fit_statuses.append(False)
        continue
    #endif

    if h.GetEntries() < 10:
        fit_functions.append(None)
        fit_statuses.append(False)
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

    fit_result = h.Fit(fit_func, "RQ0S")

    fit_ok = False
    if fit_result and int(fit_result.Status()) == 0:
        fit_ok = True
    #endif

    fit_functions.append(fit_func)
    fit_statuses.append(fit_ok)

    if fit_ok:
        mu = fit_func.GetParameter(1)
        sigma = abs(fit_func.GetParameter(2))

        x_vals.append(mean_vz)
        y_vals.append(mu)
        y_errs.append(sigma)
    #endif
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
# canvas layout
# left side: multi-panel of Mx2 fits
# right side: extracted peak position vs vz_e
# ------------------------------------------------

canvas = ROOT.TCanvas("canvas", "external radiation estimation", 1800, 950)

left_margin_frac = 0.72
left_pad = ROOT.TPad("left_pad", "left_pad", 0.00, 0.00, left_margin_frac, 1.00)
right_pad = ROOT.TPad("right_pad", "right_pad", left_margin_frac, 0.00, 1.00, 1.00)

left_pad.Draw()
right_pad.Draw()

# ------------------------------------------------
# left pad split into multi-panel layout
# ------------------------------------------------

left_pad.cd()
left_pad.Divide(4, 2, 0.001, 0.001)

legend = ROOT.TLegend(0.50, 0.72, 0.88, 0.88)
legend.SetBorderSize(1)
legend.SetFillStyle(1001)
legend.SetFillColor(ROOT.kWhite)
legend.SetTextSize(0.026)

legend_added_hist = False
legend_added_fit = False

for i in range(n_vz_bins):
    left_pad.cd(i + 1)

    ROOT.gPad.SetLeftMargin(0.14)
    ROOT.gPad.SetRightMargin(0.05)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.10)

    h = histograms[i]
    h.SetTitle("")
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.7)
    h.SetLineWidth(1)

    local_max = h.GetMaximum()
    y_max = 1.25 * local_max if local_max > 0.0 else 1.0

    h.SetMinimum(0.0)
    h.SetMaximum(y_max)

    h.GetXaxis().SetTitle("M_{X}^{2} (GeV^{2})")
    h.GetYaxis().SetTitle("Normalized counts")
    h.GetXaxis().CenterTitle()
    h.GetYaxis().CenterTitle()
    h.GetXaxis().SetTitleSize(0.06)
    h.GetYaxis().SetTitleSize(0.06)
    h.GetXaxis().SetLabelSize(0.05)
    h.GetYaxis().SetLabelSize(0.05)
    h.GetYaxis().SetTitleOffset(1.10)

    h.Draw("E1")

    if fit_functions[i]:
        fit_functions[i].SetLineColor(ROOT.kRed + 1)
        fit_functions[i].SetLineWidth(2)
        fit_functions[i].Draw("same")
    #endif

    vz_lo = vz_edges[i]
    vz_hi = vz_edges[i + 1]

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.060)
    latex.DrawLatex(0.18, 0.86, "%.1f < v_{z,e} < %.1f (cm)" % (vz_lo, vz_hi))

    if fit_statuses[i]:
        mu = fit_functions[i].GetParameter(1)
        sigma = abs(fit_functions[i].GetParameter(2))
        latex.SetTextSize(0.050)
        latex.DrawLatex(0.18, 0.78, "#mu = %.4f" % mu)
        latex.DrawLatex(0.18, 0.70, "#sigma = %.4f" % sigma)
    #endif

    if i == 0:
        if not legend_added_hist:
            legend.AddEntry(h, "Data", "lep")
            legend_added_hist = True
        #endif
        if fit_functions[i] and not legend_added_fit:
            legend.AddEntry(fit_functions[i], "Gaussian + pol2", "l")
            legend_added_fit = True
        #endif
    #endif
#endfor

# draw the legend in the last panel for convenience
left_pad.cd(n_vz_bins)
legend.Draw()

# ------------------------------------------------
# right pad: fitted peak position with sigma error bars vs mean vz_e
# ------------------------------------------------

right_pad.cd()
ROOT.gPad.SetLeftMargin(0.16)
ROOT.gPad.SetRightMargin(0.06)
ROOT.gPad.SetBottomMargin(0.12)
ROOT.gPad.SetTopMargin(0.08)

frame_right = ROOT.TH1D("frame_right", "", 100, vz_min, vz_max)
frame_right.SetMinimum(0.7)
frame_right.SetMaximum(1.1)
frame_right.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
frame_right.GetYaxis().SetTitle("Fitted neutron peak position M_{X}^{2} (GeV^{2})")
frame_right.GetXaxis().CenterTitle()
frame_right.GetYaxis().CenterTitle()
frame_right.GetXaxis().SetTitleSize(0.05)
frame_right.GetYaxis().SetTitleSize(0.05)
frame_right.GetXaxis().SetLabelSize(0.04)
frame_right.GetYaxis().SetLabelSize(0.04)
frame_right.GetYaxis().SetTitleOffset(1.35)
frame_right.Draw()

line_neutron = ROOT.TLine(vz_min, M_n2, vz_max, M_n2)
line_neutron.SetLineStyle(2)
line_neutron.SetLineWidth(2)
line_neutron.Draw("same")

graph.Draw("P SAME")

legend_right = ROOT.TLegend(0.18, 0.78, 0.72, 0.90)
legend_right.SetBorderSize(1)
legend_right.SetFillStyle(1001)
legend_right.SetFillColor(ROOT.kWhite)
legend_right.SetTextSize(0.028)
legend_right.AddEntry(graph, "Point = fitted mean, error bar = fitted #sigma", "lep")
legend_right.AddEntry(line_neutron, "Expected neutron mass squared", "l")
legend_right.Draw()

latex_right = ROOT.TLatex()
latex_right.SetNDC()
latex_right.SetTextSize(0.042)
latex_right.DrawLatex(0.18, 0.93, "Fitted peak position vs v_{z,e}")

# ------------------------------------------------
# save
# ------------------------------------------------

canvas.SaveAs(output_png)

print("Saved:", output_png)

f.Close()