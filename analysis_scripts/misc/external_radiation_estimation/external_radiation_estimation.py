import ROOT
import os
import math

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# ------------------------------------------------
# input files
# NOTE:
# You pasted the 2nd and 3rd paths as the same string.
# I am keeping the 3rd one as a separate variable so you can
# replace it with the actual epi+pi- data file path if needed.
# ------------------------------------------------

input_file_data_epi = "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/enpi+/rga_fa18_inb_epi+.root"
input_file_mc_epi   = "/work/clas12/thayward/CLAS12_exclusive/enpi+/mc/rec_clasdis_rga_fa18_inb_epi+X.root"
input_file_data_epipim = "/work/clas12/thayward/CLAS12_exclusive/enpi+/mc/rec_clasdis_rga_fa18_inb_epi+X.root"

tree_name = "PhysicsEvents"
output_dir = "output"
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------
# physics constants
# ------------------------------------------------

M_p = 0.938272
M_n = 0.939565

M_p2 = M_p * M_p
M_n2 = M_n * M_n

# ------------------------------------------------
# binning
# twice as many vz bins as before:
# 0.5 cm bins from -6.5 to 1.5
# ------------------------------------------------

vz_min = -6.5
vz_max =  1.5
vz_step = 0.5

mx2_min = 0.6
mx2_max = 1.2
mx2_bins = 60

vz_edges = []
current_edge = vz_min
while current_edge <= vz_max + 1.0e-9:
    vz_edges.append(current_edge)
    current_edge += vz_step
#endfor

n_vz_bins = len(vz_edges) - 1

# ------------------------------------------------
# helpers
# ------------------------------------------------

def make_safe_tag(text):
    out = text.lower()
    out = out.replace(" ", "_")
    out = out.replace("+", "p")
    out = out.replace("-", "m")
    out = out.replace("/", "_")
    return out
#enddef

def process_file(input_file, dataset_label, missing_particle_label, expected_mass2):
    f = ROOT.TFile.Open(input_file, "READ")
    if not f or f.IsZombie():
        raise RuntimeError("Could not open file: " + input_file)
    #endif

    tree = f.Get(tree_name)
    if not tree:
        raise RuntimeError("Could not find tree: " + tree_name + " in " + input_file)
    #endif

    histograms = []
    sum_vz = [0.0] * n_vz_bins
    count_vz = [0] * n_vz_bins

    for i in range(n_vz_bins):
        h = ROOT.TH1D(
            "h_%s_%d" % (make_safe_tag(dataset_label), i),
            "",
            mx2_bins,
            mx2_min,
            mx2_max
        )
        h.Sumw2()
        histograms.append(h)
    #endfor

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

    for i in range(n_vz_bins):
        integral = histograms[i].Integral()
        if integral > 0.0:
            histograms[i].Scale(1.0 / integral)
        #endif
    #endfor

    x_mu_vals = []
    y_mu_vals = []
    y_mu_errs = []

    x_sigma_vals = []
    y_sigma_vals = []
    y_sigma_errs = []

    fit_functions = []
    fit_statuses = []

    print("")
    print("==============================================================")
    print("Dataset: %s" % dataset_label)
    print("Input:   %s" % input_file)
    print("Expected missing %s peak mass^2 = %.6f GeV^2" % (missing_particle_label, expected_mass2))
    print("mu offsets by vz_e bin:")
    print("vz_low   vz_high   <vz_e>      mu           mu_err       offset=mu-Mexp2")
    print("--------------------------------------------------------------")

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

        fit_name = "fit_%s_%d" % (make_safe_tag(dataset_label), i)
        fit_func = ROOT.TF1(fit_name, "gaus(0)+pol3(3)", mx2_min, mx2_max)

        fit_func.SetParameter(0, peak_height)
        fit_func.SetParameter(1, expected_mass2)
        fit_func.SetParameter(2, 0.03)
        fit_func.SetParameter(3, 0.0)
        fit_func.SetParameter(4, 0.0)
        fit_func.SetParameter(5, 0.0)
        fit_func.SetParameter(6, 0.0)

        fit_func.SetParLimits(0, 0.0, 10.0)
        fit_func.SetParLimits(1, 0.75, 1.00)
        fit_func.SetParLimits(2, 0.003, 0.20)

        fit_result = h.Fit(fit_func, "RQ0S")

        fit_ok = False
        if fit_result and int(fit_result.Status()) == 0:
            fit_ok = True
        #endif

        fit_functions.append(fit_func)
        fit_statuses.append(fit_ok)

        if fit_ok:
            mu        = fit_func.GetParameter(1)
            mu_err    = fit_func.GetParError(1)
            sigma     = abs(fit_func.GetParameter(2))
            sigma_err = fit_func.GetParError(2)
            offset    = mu - expected_mass2

            x_mu_vals.append(mean_vz)
            y_mu_vals.append(mu)
            y_mu_errs.append(mu_err)

            x_sigma_vals.append(mean_vz)
            y_sigma_vals.append(sigma)
            y_sigma_errs.append(sigma_err)

            print("%6.2f   %7.2f   %8.4f   %10.6f   %10.6f   %+12.6f" % (
                vz_edges[i],
                vz_edges[i + 1],
                mean_vz,
                mu,
                mu_err,
                offset
            ))
        #endif
    #endfor

    n_mu_points = len(x_mu_vals)
    graph_mu = ROOT.TGraphErrors(n_mu_points)
    for i in range(n_mu_points):
        graph_mu.SetPoint(i, x_mu_vals[i], y_mu_vals[i])
        graph_mu.SetPointError(i, 0.0, y_mu_errs[i])
    #endfor

    graph_mu.SetTitle("")
    graph_mu.SetMarkerStyle(20)
    graph_mu.SetMarkerSize(1.0)
    graph_mu.SetLineWidth(2)

    n_sigma_points = len(x_sigma_vals)
    graph_sigma = ROOT.TGraphErrors(n_sigma_points)
    for i in range(n_sigma_points):
        graph_sigma.SetPoint(i, x_sigma_vals[i], y_sigma_vals[i])
        graph_sigma.SetPointError(i, 0.0, y_sigma_errs[i])
    #endfor

    graph_sigma.SetTitle("")
    graph_sigma.SetMarkerStyle(20)
    graph_sigma.SetMarkerSize(1.0)
    graph_sigma.SetLineWidth(2)

    cubic_fit = None
    if n_mu_points >= 4:
        cubic_fit = ROOT.TF1(
            "cubic_mu_%s" % make_safe_tag(dataset_label),
            "[0] + [1]*x + [2]*x*x + [3]*x*x*x",
            vz_min,
            vz_max
        )
        graph_mu.Fit(cubic_fit, "Q0")
        cubic_fit.SetLineColor(ROOT.kRed + 1)
        cubic_fit.SetLineStyle(2)
        cubic_fit.SetLineWidth(2)

        c0 = cubic_fit.GetParameter(0)
        c1 = cubic_fit.GetParameter(1)
        c2 = cubic_fit.GetParameter(2)
        c3 = cubic_fit.GetParameter(3)

        print("")
        print("Cubic fit for mu(vz_e):")
        print("mu(vz_e) = %.8f + %.8f*vz_e + %.8f*vz_e^2 + %.8f*vz_e^3" % (c0, c1, c2, c3))
    #endif

    canvas = ROOT.TCanvas(
        "canvas_%s" % make_safe_tag(dataset_label),
        dataset_label,
        2200,
        1300
    )

    left_margin_frac = 0.70
    left_pad = ROOT.TPad("left_pad_%s" % make_safe_tag(dataset_label), "", 0.00, 0.00, left_margin_frac, 1.00)
    right_pad = ROOT.TPad("right_pad_%s" % make_safe_tag(dataset_label), "", left_margin_frac, 0.00, 1.00, 1.00)

    left_pad.Draw()
    right_pad.Draw()

    left_pad.cd()
    left_pad.Divide(4, 4, 0.001, 0.001)

    for i in range(n_vz_bins):
        left_pad.cd(i + 1)

        ROOT.gPad.SetLeftMargin(0.19)
        ROOT.gPad.SetRightMargin(0.05)
        ROOT.gPad.SetBottomMargin(0.14)
        ROOT.gPad.SetTopMargin(0.10)

        h = histograms[i]
        h.SetTitle("")
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.55)
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
        h.GetYaxis().SetTitleSize(0.055)
        h.GetXaxis().SetLabelSize(0.05)
        h.GetYaxis().SetLabelSize(0.05)
        h.GetYaxis().SetTitleOffset(1.55)

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
        latex.SetTextSize(0.050)
        latex.DrawLatex(0.22, 0.86, "%.1f < v_{z,e} < %.1f (cm)" % (vz_lo, vz_hi))

        if fit_statuses[i]:
            mu = fit_functions[i].GetParameter(1)
            sigma = abs(fit_functions[i].GetParameter(2))
            latex.SetTextSize(0.042)
            latex.DrawLatex(0.22, 0.78, "#mu = %.4f" % mu)
            latex.DrawLatex(0.22, 0.70, "#sigma = %.4f" % sigma)
        #endif
    #endfor

    right_pad.cd()

    right_top = ROOT.TPad("right_top_%s" % make_safe_tag(dataset_label), "", 0.00, 0.50, 1.00, 1.00)
    right_bot = ROOT.TPad("right_bot_%s" % make_safe_tag(dataset_label), "", 0.00, 0.00, 1.00, 0.50)

    right_top.Draw()
    right_bot.Draw()

    right_top.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.08)

    frame_right_top = ROOT.TH1D("frame_right_top_%s" % make_safe_tag(dataset_label), "", 100, vz_min, vz_max)
    frame_right_top.SetMinimum(0.7)
    frame_right_top.SetMaximum(1.1)
    frame_right_top.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
    frame_right_top.GetYaxis().SetTitle("Fitted peak position #mu of M_{X}^{2} (GeV^{2})")
    frame_right_top.GetXaxis().CenterTitle()
    frame_right_top.GetYaxis().CenterTitle()
    frame_right_top.GetXaxis().SetTitleSize(0.05)
    frame_right_top.GetYaxis().SetTitleSize(0.05)
    frame_right_top.GetXaxis().SetLabelSize(0.04)
    frame_right_top.GetYaxis().SetLabelSize(0.04)
    frame_right_top.GetYaxis().SetTitleOffset(1.35)
    frame_right_top.Draw()

    line_expected = ROOT.TLine(vz_min, expected_mass2, vz_max, expected_mass2)
    line_expected.SetLineStyle(2)
    line_expected.SetLineWidth(2)
    line_expected.Draw("same")

    graph_mu.Draw("P SAME")

    if cubic_fit:
        cubic_fit.Draw("same")
    #endif

    legend_right_top = ROOT.TLegend(0.16, 0.74, 0.94, 0.90)
    legend_right_top.SetBorderSize(1)
    legend_right_top.SetFillStyle(1001)
    legend_right_top.SetFillColor(ROOT.kWhite)
    legend_right_top.SetTextSize(0.023)
    legend_right_top.AddEntry(graph_mu, "Point = fitted #mu, error bar = fit uncertainty on #mu", "lep")
    legend_right_top.AddEntry(line_expected, "Expected missing %s mass^{2}" % missing_particle_label, "l")
    if cubic_fit:
        legend_right_top.AddEntry(cubic_fit, "Cubic fit to #mu(v_{z,e})", "l")
    #endif
    legend_right_top.Draw()

    latex_right_top = ROOT.TLatex()
    latex_right_top.SetNDC()
    latex_right_top.SetTextSize(0.038)
    latex_right_top.DrawLatex(0.16, 0.93, "%s: fitted #mu vs v_{z,e}" % dataset_label)

    right_bot.cd()
    ROOT.gPad.SetLeftMargin(0.16)
    ROOT.gPad.SetRightMargin(0.06)
    ROOT.gPad.SetBottomMargin(0.14)
    ROOT.gPad.SetTopMargin(0.08)

    sigma_ymax = 0.10
    for val, err in zip(y_sigma_vals, y_sigma_errs):
        if val + err > sigma_ymax:
            sigma_ymax = val + err
        #endif
    #endfor
    sigma_ymax *= 1.25

    frame_right_bot = ROOT.TH1D("frame_right_bot_%s" % make_safe_tag(dataset_label), "", 100, vz_min, vz_max)
    frame_right_bot.SetMinimum(0.0)
    frame_right_bot.SetMaximum(sigma_ymax)
    frame_right_bot.GetXaxis().SetTitle("<v_{z,e}> in bin (cm)")
    frame_right_bot.GetYaxis().SetTitle("Fitted #sigma of M_{X}^{2} peak (GeV^{2})")
    frame_right_bot.GetXaxis().CenterTitle()
    frame_right_bot.GetYaxis().CenterTitle()
    frame_right_bot.GetXaxis().SetTitleSize(0.05)
    frame_right_bot.GetYaxis().SetTitleSize(0.05)
    frame_right_bot.GetXaxis().SetLabelSize(0.04)
    frame_right_bot.GetYaxis().SetLabelSize(0.04)
    frame_right_bot.GetYaxis().SetTitleOffset(1.35)
    frame_right_bot.Draw()

    graph_sigma.Draw("P SAME")

    legend_right_bot = ROOT.TLegend(0.18, 0.82, 0.62, 0.90)
    legend_right_bot.SetBorderSize(1)
    legend_right_bot.SetFillStyle(1001)
    legend_right_bot.SetFillColor(ROOT.kWhite)
    legend_right_bot.SetTextSize(0.026)
    legend_right_bot.AddEntry(graph_sigma, "Point = fitted #sigma, error bar = fit uncertainty", "lep")
    legend_right_bot.Draw()

    latex_right_bot = ROOT.TLatex()
    latex_right_bot.SetNDC()
    latex_right_bot.SetTextSize(0.040)
    latex_right_bot.DrawLatex(0.18, 0.93, "%s: fitted #sigma vs v_{z,e}" % dataset_label)

    output_png = os.path.join(output_dir, "external_radiation_estimation_%s.png" % make_safe_tag(dataset_label))
    canvas.SaveAs(output_png)
    print("Saved: %s" % output_png)

    f.Close()
#enddef

# ------------------------------------------------
# run all three
# ------------------------------------------------

process_file(
    input_file_data_epi,
    "data_epi_plus",
    "neutron",
    M_n2
)

process_file(
    input_file_mc_epi,
    "mc_epi_plus",
    "neutron",
    M_n2
)

process_file(
    input_file_data_epipim,
    "data_epi_plus_pi_minus",
    "proton",
    M_p2
)