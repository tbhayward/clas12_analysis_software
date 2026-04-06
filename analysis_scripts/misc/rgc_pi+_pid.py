import ROOT
import os
import math
import sys

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

# ------------------------------------------------
# configuration
# ------------------------------------------------

OUTPUT_DIR = "output/rgc_pi+_particle_id"
TREE_NAME = "PhysicsEvents"

RUN_FILES = [
    ("Su22", "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/rgc_su22_inb_NH3_epi+X_calibration.root"),
    ("Fa22", "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/rgc_fa22_inb_NH3_epi+X_calibration.root"),
    ("Sp23", "/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/rgc_sp23_inb_NH3_epi+X_calibration.root"),
]

# Assumed detector mapping:
#   CD -> track_sector_5
#   FD -> track_sector_6
DETECTOR_CONFIGS = [
    ("CD", "track_sector_5"),
    ("FD", "track_sector_6"),
]

CD_P_BINS = [
    (0.500, 0.750),
    (0.750, 1.000),
    (1.000, 1.250),
    (1.250, 1.500),
]

FD_P_BINS = [
    (0.500, 1.000),
    (1.000, 1.500),
    (1.500, 2.000),
    (2.000, 2.500),
    (2.500, 3.000),
    (3.000, 3.500),
    (3.500, 4.000),
    (4.000, 5.000),
]

CHI2_XMIN = -5.0
CHI2_XMAX = 5.0
CHI2_NBINS = 200

FIT_XMIN = -2.0
FIT_XMAX = 2.0

BETA_P_XMIN = 0.0
BETA_P_XMAX = 8.0
BETA_P_YMIN = 0.0
BETA_P_YMAX = 1.2
BETA_P_NBINS_X = 240
BETA_P_NBINS_Y = 240

VALID_BETA_P_PIDS = set([211, 321, 2212])

# ------------------------------------------------
# helpers
# ------------------------------------------------

def ensure_output_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path)
    #endif
#endfor

def require_file(path):
    if not os.path.isfile(path):
        raise FileNotFoundError("Required input file not found: " + path)
    #endif
#endfor

def require_tree(root_file, tree_name):
    tree = root_file.Get(tree_name)
    if not tree:
        raise RuntimeError("Could not find tree '" + tree_name + "' in file " + root_file.GetName())
    #endif
    return tree
#endfor

def require_branches(tree, branch_names):
    missing = []
    for branch_name in branch_names:
        if not tree.GetBranch(branch_name):
            missing.append(branch_name)
        #endif
    #endfor

    if len(missing) > 0:
        raise RuntimeError(
            "Missing required branches in tree '" + tree.GetName() + "' from file '" + tree.GetCurrentFile().GetName() + "': "
            + ", ".join(missing)
        )
    #endif
#endfor

def set_branch_statuses(tree, active_branches):
    tree.SetBranchStatus("*", 0)
    for name in active_branches:
        tree.SetBranchStatus(name, 1)
    #endfor
#endfor

def get_p_bins(det_label):
    if det_label == "CD":
        return CD_P_BINS
    #endif
    return FD_P_BINS
#endfor

def get_panel_layout(det_label):
    if det_label == "CD":
        return 3, 2
    #endif
    return 5, 2
#endfor

def get_graph_xmax(det_label):
    if det_label == "CD":
        return 1.5
    #endif
    return 5.0
#endfor

def style_hist(hist, x_title, y_title):
    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)
    hist.GetXaxis().SetTitleSize(0.060)
    hist.GetYaxis().SetTitleSize(0.060)
    hist.GetXaxis().SetLabelSize(0.050)
    hist.GetYaxis().SetLabelSize(0.050)
    hist.GetYaxis().SetTitleOffset(1.05)
    hist.GetXaxis().SetRangeUser(CHI2_XMIN, CHI2_XMAX)
#endfor

def style_frame_hist(hist, x_title, y_title):
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)
    hist.GetXaxis().SetTitleSize(0.060)
    hist.GetYaxis().SetTitleSize(0.060)
    hist.GetXaxis().SetLabelSize(0.050)
    hist.GetYaxis().SetLabelSize(0.050)
    hist.GetYaxis().SetTitleOffset(1.05)
#endfor

def style_hist2d(hist, x_title, y_title, z_title):
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetZaxis().SetTitle(z_title)
    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)
    hist.GetZaxis().CenterTitle(True)
    hist.GetXaxis().SetTitleSize(0.050)
    hist.GetYaxis().SetTitleSize(0.050)
    hist.GetZaxis().SetTitleSize(0.050)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetZaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.10)
    hist.GetZaxis().SetTitleOffset(1.10)
#endfor

def nice_bin_label(pmin, pmax):
    return "%.2f < p < %.2f (GeV)" % (pmin, pmax)
#endfor

def create_histograms(run_label, det_label, p_bins):
    hists = []

    h_int = ROOT.TH1D(
        "h_chi2_integrated_%s_%s" % (run_label, det_label),
        "",
        CHI2_NBINS,
        CHI2_XMIN,
        CHI2_XMAX
    )
    h_int.SetDirectory(0)
    hists.append(h_int)

    for i_bin, _ in enumerate(p_bins):
        hist = ROOT.TH1D(
            "h_chi2_%s_%s_bin_%d" % (run_label, det_label, i_bin),
            "",
            CHI2_NBINS,
            CHI2_XMIN,
            CHI2_XMAX
        )
        hist.SetDirectory(0)
        hists.append(hist)
    #endfor

    return hists
#endfor

def create_beta_histogram(run_label, det_label):
    hist = ROOT.TH2D(
        "h2_beta_p_%s_%s" % (run_label, det_label),
        "",
        BETA_P_NBINS_X,
        BETA_P_XMIN,
        BETA_P_XMAX,
        BETA_P_NBINS_Y,
        BETA_P_YMIN,
        BETA_P_YMAX
    )
    hist.SetDirectory(0)
    return hist
#endfor

def make_title_canvas(name, width, height, title_text, ncols, nrows):
    canvas = ROOT.TCanvas(name, "", width, height)

    title_pad = ROOT.TPad(name + "_title", "", 0.0, 0.92, 1.0, 1.0)
    grid_pad = ROOT.TPad(name + "_grid", "", 0.0, 0.0, 1.0, 0.92)

    title_pad.SetFillColor(0)
    title_pad.SetBorderMode(0)
    title_pad.SetMargin(0.02, 0.02, 0.05, 0.05)
    title_pad.Draw()

    grid_pad.SetFillColor(0)
    grid_pad.SetBorderMode(0)
    grid_pad.SetMargin(0.0, 0.0, 0.0, 0.0)
    grid_pad.Draw()
    grid_pad.Divide(ncols, nrows, 0.001, 0.001)

    title_pad.cd()
    title = ROOT.TLatex()
    title.SetNDC(True)
    title.SetTextAlign(22)
    title.SetTextSize(0.50)
    title.DrawLatex(0.50, 0.50, title_text)

    canvas.cd()

    return canvas, title_pad, grid_pad
#endfor

def determine_fit_range(hist, det_label, i_bin_or_none):
    peak_bin = hist.GetMaximumBin()
    peak_x = hist.GetBinCenter(peak_bin)

    fit_min = FIT_XMIN
    fit_max = FIT_XMAX
    fit_mode = "full"

    # For CD, bins after the second bin are fit only on the left side.
    # With the current CD bins, that means:
    #   [1.00,1.25] and [1.25,1.50]
    if det_label == "CD" and i_bin_or_none is not None and i_bin_or_none >= 2:
        fit_mode = "left-only"
        fit_max = peak_x
        if fit_max > FIT_XMAX:
            fit_max = FIT_XMAX
        #endif
        if fit_max < -0.20:
            fit_max = -0.20
        #endif
    #endif

    if fit_max <= fit_min + 0.20:
        fit_min = FIT_XMIN
        fit_max = FIT_XMAX
        fit_mode = "full-fallback"
    #endif

    return fit_min, fit_max, fit_mode
#endfor

def fit_histogram_gaussian(hist, fit_name, det_label, i_bin_or_none):
    n_entries = int(hist.GetEntries())

    if n_entries < 10 or hist.Integral() <= 0:
        return None
    #endif

    fit_min, fit_max, fit_mode = determine_fit_range(hist, det_label, i_bin_or_none)

    func = ROOT.TF1(fit_name, "gaus", CHI2_XMIN, CHI2_XMAX)

    peak_bin = hist.GetMaximumBin()
    peak_x = hist.GetBinCenter(peak_bin)
    peak_y = hist.GetBinContent(peak_bin)

    initial_sigma = hist.GetRMS()
    if initial_sigma <= 0.0:
        initial_sigma = 0.5
    #endif
    if initial_sigma > 2.0:
        initial_sigma = 2.0
    #endif

    func.SetParameters(peak_y, peak_x, initial_sigma)

    fit_result = hist.Fit(func, "RQSN", "", fit_min, fit_max)

    if int(fit_result) != 0:
        return None
    #endif

    return {
        "func": func,
        "fit_min": fit_min,
        "fit_max": fit_max,
        "fit_mode": fit_mode,
    }
#endfor

def build_function_segments(fit_info, seg_prefix):
    if fit_info is None:
        return []
    #endif

    func = fit_info["func"]
    fit_min = fit_info["fit_min"]
    fit_max = fit_info["fit_max"]

    amp = func.GetParameter(0)
    mean = func.GetParameter(1)
    sigma = func.GetParameter(2)

    segments = []

    if fit_min > CHI2_XMIN:
        left_seg = ROOT.TF1(seg_prefix + "_left", "gaus", CHI2_XMIN, fit_min)
        left_seg.SetParameters(amp, mean, sigma)
        left_seg.SetLineColor(ROOT.kRed + 1)
        left_seg.SetLineStyle(2)
        left_seg.SetLineWidth(2)
        segments.append(left_seg)
    #endif

    center_seg = ROOT.TF1(seg_prefix + "_center", "gaus", fit_min, fit_max)
    center_seg.SetParameters(amp, mean, sigma)
    center_seg.SetLineColor(ROOT.kRed + 1)
    center_seg.SetLineStyle(1)
    center_seg.SetLineWidth(2)
    segments.append(center_seg)

    if fit_max < CHI2_XMAX:
        right_seg = ROOT.TF1(seg_prefix + "_right", "gaus", fit_max, CHI2_XMAX)
        right_seg.SetParameters(amp, mean, sigma)
        right_seg.SetLineColor(ROOT.kRed + 1)
        right_seg.SetLineStyle(2)
        right_seg.SetLineWidth(2)
        segments.append(right_seg)
    #endif

    return segments
#endfor

def draw_info_box(fit_info, hist, fit_result_dict):
    box = ROOT.TPaveText(0.56, 0.68, 0.94, 0.90, "NDC")
    box.SetFillColor(ROOT.kWhite)
    box.SetFillStyle(1001)
    box.SetBorderSize(1)
    box.SetTextAlign(12)
    box.SetTextSize(0.050)

    entries = int(hist.GetEntries())
    box.AddText("N = %d" % entries)

    if fit_info is None:
        box.AddText("fit failed")
    else:
        box.AddText("#mu = %.3f #pm %.3f" % (fit_result_dict["mu"], fit_result_dict["mu_err"]))
        box.AddText("#sigma = %.3f #pm %.3f" % (fit_result_dict["sigma"], fit_result_dict["sigma_err"]))
        box.AddText("fit: [%.2f, %.2f]" % (fit_info["fit_min"], fit_info["fit_max"]))
    #endif

    box.Draw()
    return box
#endfor

def fill_histograms(tree, run_label, det_label, detector_branch, p_bins):
    hists = create_histograms(run_label, det_label, p_bins)

    p_sum = [0.0 for _ in range(len(p_bins))]
    p_count = [0 for _ in range(len(p_bins))]

    n_total = tree.GetEntries()

    print("")
    print("============================================================")
    print("Processing %s %s" % (run_label, det_label))
    print("File: %s" % tree.GetCurrentFile().GetName())
    print("Tree entries: %d" % n_total)
    print("Using detector branch: %s" % detector_branch)
    print("============================================================")

    for i_entry in range(n_total):
        tree.GetEntry(i_entry)

        if i_entry > 0 and i_entry % 1000000 == 0:
            print("  processed %d / %d entries" % (i_entry, n_total))
        #endif

        pid = int(tree.particle_pid)
        if pid != 211:
            continue
        #endif

        detector_value = float(getattr(tree, detector_branch))
        if detector_value == -9999:
            continue
        #endif

        p_val = float(tree.p)
        chi2_val = float(tree.particle_chi2pid)

        hists[0].Fill(chi2_val)

        for i_bin, (pmin, pmax) in enumerate(p_bins):
            in_bin = False

            if i_bin < len(p_bins) - 1:
                if p_val >= pmin and p_val < pmax:
                    in_bin = True
                #endif
            else:
                if p_val >= pmin and p_val <= pmax:
                    in_bin = True
                #endif
            #endif

            if in_bin:
                hists[i_bin + 1].Fill(chi2_val)
                p_sum[i_bin] += p_val
                p_count[i_bin] += 1
                break
            #endif
        #endfor
    #endfor

    p_means = []
    for i_bin in range(len(p_bins)):
        if p_count[i_bin] > 0:
            p_means.append(p_sum[i_bin] / float(p_count[i_bin]))
        else:
            p_means.append(0.5 * (p_bins[i_bin][0] + p_bins[i_bin][1]))
        #endif
    #endfor

    return hists, p_means, p_count
#endfor

def print_fit_summary(run_label, det_label, p_bins, fit_results):
    print("")
    print("------------------------------------------------------------")
    print("Gaussian fit summary for %s %s" % (run_label, det_label))
    print("------------------------------------------------------------")
    print("%-24s  %-10s  %-10s  %-10s  %-10s  %-10s  %-14s" % (
        "Bin", "Entries", "mu", "sigma", "mu err", "sigma err", "fit range"
    ))

    integrated = fit_results["integrated"]
    if integrated["fit_info"] is not None:
        print(
            "%-24s  %-10d  %-10.5f  %-10.5f  %-10.5f  %-10.5f  [%.2f, %.2f]" % (
                "Integrated",
                integrated["entries"],
                integrated["mu"],
                integrated["sigma"],
                integrated["mu_err"],
                integrated["sigma_err"],
                integrated["fit_info"]["fit_min"],
                integrated["fit_info"]["fit_max"],
            )
        )
    else:
        print(
            "%-24s  %-10d  %-10s  %-10s  %-10s  %-10s  %-14s" % (
                "Integrated",
                integrated["entries"],
                "failed",
                "failed",
                "failed",
                "failed",
                "failed"
            )
        )
    #endif

    for i_bin, result in enumerate(fit_results["bins"]):
        label = "%.2f-%.2f GeV" % (p_bins[i_bin][0], p_bins[i_bin][1])
        if result["fit_info"] is not None:
            print(
                "%-24s  %-10d  %-10.5f  %-10.5f  %-10.5f  %-10.5f  [%.2f, %.2f]" % (
                    label,
                    result["entries"],
                    result["mu"],
                    result["sigma"],
                    result["mu_err"],
                    result["sigma_err"],
                    result["fit_info"]["fit_min"],
                    result["fit_info"]["fit_max"],
                )
            )
        else:
            print(
                "%-24s  %-10d  %-10s  %-10s  %-10s  %-10s  %-14s" % (
                    label,
                    result["entries"],
                    "failed",
                    "failed",
                    "failed",
                    "failed",
                    "failed"
                )
            )
        #endif
    #endfor

    print("------------------------------------------------------------")
    print("")
#endfor

def build_graph_from_results(run_label, det_label, bin_results):
    x = []
    y = []
    ex = []
    ey = []

    for result in bin_results:
        if result["fit_info"] is None:
            continue
        #endif

        x.append(result["pmean"])
        y.append(result["mu"])
        ex.append(0.0)
        ey.append(result["sigma"])
    #endfor

    graph = ROOT.TGraphErrors(len(x))
    graph.SetName("g_mu_vs_p_%s_%s" % (run_label, det_label))

    for i in range(len(x)):
        graph.SetPoint(i, x[i], y[i])
        graph.SetPointError(i, ex[i], ey[i])
    #endfor

    graph.SetMarkerStyle(20)
    graph.SetMarkerSize(1.0)
    graph.SetLineWidth(2)

    return graph
#endfor

def determine_graph_y_range(bin_results):
    y_min = 999.0
    y_max = -999.0

    for result in bin_results:
        if result["fit_info"] is None:
            continue
        #endif
        low = result["mu"] - result["sigma"]
        high = result["mu"] + result["sigma"]
        if low < y_min:
            y_min = low
        #endif
        if high > y_max:
            y_max = high
        #endif
    #endfor

    if y_min > y_max:
        y_min = -2.0
        y_max = 2.0
    #endif

    pad = 0.20 * (y_max - y_min)
    if pad <= 0.0:
        pad = 1.0
    #endif

    return y_min - pad, y_max + pad
#endfor

def analyze_one_detector(tree, run_label, det_label, detector_branch):
    p_bins = get_p_bins(det_label)
    hists, p_means, p_counts = fill_histograms(tree, run_label, det_label, detector_branch, p_bins)

    fit_results = {
        "integrated": None,
        "bins": []
    }

    h_int = hists[0]
    int_fit_info = fit_histogram_gaussian(h_int, "f_int_%s_%s" % (run_label, det_label), det_label, None)

    fit_results["integrated"] = {
        "entries": int(h_int.GetEntries()),
        "fit_info": int_fit_info,
        "mu": int_fit_info["func"].GetParameter(1) if int_fit_info else None,
        "sigma": abs(int_fit_info["func"].GetParameter(2)) if int_fit_info else None,
        "mu_err": int_fit_info["func"].GetParError(1) if int_fit_info else None,
        "sigma_err": int_fit_info["func"].GetParError(2) if int_fit_info else None
    }

    for i_bin, (pmin, pmax) in enumerate(p_bins):
        hist = hists[i_bin + 1]
        fit_info = fit_histogram_gaussian(hist, "f_bin_%s_%s_%d" % (run_label, det_label, i_bin), det_label, i_bin)

        result = {
            "pmin": pmin,
            "pmax": pmax,
            "pmean": p_means[i_bin],
            "entries": int(hist.GetEntries()),
            "fit_info": fit_info,
            "mu": fit_info["func"].GetParameter(1) if fit_info else None,
            "sigma": abs(fit_info["func"].GetParameter(2)) if fit_info else None,
            "mu_err": fit_info["func"].GetParError(1) if fit_info else None,
            "sigma_err": fit_info["func"].GetParError(2) if fit_info else None
        }
        fit_results["bins"].append(result)
    #endfor

    print_fit_summary(run_label, det_label, p_bins, fit_results)

    ncols, nrows = get_panel_layout(det_label)

    if det_label == "CD":
        canvas_width = 1800
        canvas_height = 1000
    else:
        canvas_width = 2400
        canvas_height = 1000
    #endif

    global_title = "%s %s, particle_pid = 211" % (run_label, det_label)
    canvas, title_pad, grid_pad = make_title_canvas(
        "c_%s_%s" % (run_label, det_label),
        canvas_width,
        canvas_height,
        global_title,
        ncols,
        nrows
    )

    panel_title = ROOT.TLatex()
    panel_title.SetNDC(True)
    panel_title.SetTextAlign(22)
    panel_title.SetTextSize(0.075)

    objects_to_keep = []
    total_panels = 1 + len(p_bins) + 1

    for ipad in range(1, ncols * nrows + 1):
        pad = grid_pad.cd(ipad)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.18)
        if ipad > total_panels:
            pad.Clear()
        #endif
    #endfor

    grid_pad.cd(1)
    style_hist(h_int, "chi2pid", "Counts")
    h_int.SetTitle("")
    h_int.Draw("hist")
    panel_title.DrawLatex(0.50, 0.95, "Integrated")
    int_segments = build_function_segments(fit_results["integrated"]["fit_info"], "seg_int_%s_%s" % (run_label, det_label))
    for seg in int_segments:
        seg.Draw("same")
        objects_to_keep.append(seg)
    #endfor
    int_box = draw_info_box(fit_results["integrated"]["fit_info"], h_int, fit_results["integrated"])
    objects_to_keep.append(int_box)

    for i_bin, result in enumerate(fit_results["bins"]):
        grid_pad.cd(i_bin + 2)
        hist = hists[i_bin + 1]
        style_hist(hist, "chi2pid", "Counts")
        hist.SetTitle("")
        hist.Draw("hist")
        panel_title.DrawLatex(0.50, 0.95, nice_bin_label(result["pmin"], result["pmax"]))
        segments = build_function_segments(result["fit_info"], "seg_%s_%s_%d" % (run_label, det_label, i_bin))
        for seg in segments:
            seg.Draw("same")
            objects_to_keep.append(seg)
        #endfor
        info_box = draw_info_box(result["fit_info"], hist, result)
        objects_to_keep.append(info_box)
    #endfor

    final_panel = 1 + len(p_bins) + 1
    grid_pad.cd(final_panel)
    graph = build_graph_from_results(run_label, det_label, fit_results["bins"])
    y_min, y_max = determine_graph_y_range(fit_results["bins"])
    x_max_graph = get_graph_xmax(det_label)

    frame = ROOT.TH1D(
        "frame_%s_%s" % (run_label, det_label),
        "",
        100,
        0.0,
        x_max_graph
    )
    frame.SetDirectory(0)
    style_frame_hist(frame, "p (GeV)", "#mu(chi2pid)")
    frame.SetMinimum(y_min)
    frame.SetMaximum(y_max)
    frame.SetTitle("")
    frame.Draw()
    line_zero = ROOT.TLine(0.0, 0.0, x_max_graph, 0.0)
    line_zero.SetLineStyle(1)
    line_zero.SetLineWidth(1)
    line_zero.Draw("same")
    graph.Draw("P SAME")
    panel_title.DrawLatex(0.50, 0.95, "#mu(p) with #sigma as y error")

    objects_to_keep.append(frame)
    objects_to_keep.append(graph)
    objects_to_keep.append(line_zero)

    out_name = os.path.join(
        OUTPUT_DIR,
        "particle_chi2pid_pip_%s_%s.png" % (run_label, det_label)
    )
    canvas.SaveAs(out_name)

    print("Saved plot: %s" % out_name)

    return canvas, hists, graph, objects_to_keep
#endfor

def fill_beta_histogram(tree, run_label, det_label, detector_branch):
    hist = create_beta_histogram(run_label, det_label)

    n_total = tree.GetEntries()
    n_selected = 0

    print("")
    print("============================================================")
    print("Building beta vs p histogram for %s %s" % (run_label, det_label))
    print("File: %s" % tree.GetCurrentFile().GetName())
    print("Tree entries: %d" % n_total)
    print("Using detector branch: %s" % detector_branch)
    print("Included particle_pid: 211, 321, 2212")
    print("============================================================")

    for i_entry in range(n_total):
        tree.GetEntry(i_entry)

        if i_entry > 0 and i_entry % 1000000 == 0:
            print("  processed %d / %d entries" % (i_entry, n_total))
        #endif

        pid = int(tree.particle_pid)
        if pid not in VALID_BETA_P_PIDS:
            continue
        #endif

        detector_value = float(getattr(tree, detector_branch))
        if detector_value == -9999:
            continue
        #endif

        p_val = float(tree.p)
        beta_val = float(tree.particle_beta)

        hist.Fill(p_val, beta_val)
        n_selected += 1
    #endfor

    print("Selected tracks for beta vs p: %d" % n_selected)
    print("")

    return hist, n_selected
#endfor

def make_beta_vs_p_plot(tree, run_label, det_label, detector_branch):
    hist, n_selected = fill_beta_histogram(tree, run_label, det_label, detector_branch)

    title_text = "%s %s, beta vs p for particle_pid = 211, 321, 2212" % (run_label, det_label)
    canvas, title_pad, plot_pad = make_title_canvas(
        "c_beta_p_%s_%s" % (run_label, det_label),
        1200,
        900,
        title_text,
        1,
        1
    )

    plot_pad.cd(1)
    plot_pad.SetLeftMargin(0.12)
    plot_pad.SetRightMargin(0.15)
    plot_pad.SetBottomMargin(0.12)
    plot_pad.SetTopMargin(0.06)

    style_hist2d(hist, "p (GeV)", "beta", "Counts")
    hist.SetTitle("")
    hist.Draw("COLZ")

    info_box = ROOT.TPaveText(0.66, 0.78, 0.92, 0.90, "NDC")
    info_box.SetFillColor(ROOT.kWhite)
    info_box.SetFillStyle(1001)
    info_box.SetBorderSize(1)
    info_box.SetTextAlign(12)
    info_box.SetTextSize(0.035)
    info_box.AddText("Included PID:")
    info_box.AddText("211, 321, 2212")
    info_box.AddText("N = %d" % n_selected)
    info_box.Draw()

    out_name = os.path.join(
        OUTPUT_DIR,
        "beta_vs_p_combined_%s_%s.png" % (run_label, det_label)
    )
    canvas.SaveAs(out_name)

    print("Saved plot: %s" % out_name)

    return canvas, hist, info_box
#endfor

# ------------------------------------------------
# main
# ------------------------------------------------

def main():
    ensure_output_dir(OUTPUT_DIR)

    for _, file_path in RUN_FILES:
        require_file(file_path)
    #endfor

    for run_label, file_path in RUN_FILES:
        root_file = ROOT.TFile.Open(file_path, "READ")
        if not root_file or root_file.IsZombie():
            raise RuntimeError("Could not open ROOT file: " + file_path)
        #endif

        tree = require_tree(root_file, TREE_NAME)

        required_branches = [
            "particle_pid",
            "p",
            "particle_beta",
            "particle_chi2pid",
            "track_sector_5",
            "track_sector_6",
        ]
        require_branches(tree, required_branches)
        set_branch_statuses(tree, required_branches)

        for det_label, detector_branch in DETECTOR_CONFIGS:
            analyze_one_detector(tree, run_label, det_label, detector_branch)
            make_beta_vs_p_plot(tree, run_label, det_label, detector_branch)
        #endfor

        root_file.Close()
    #endfor

    print("")
    print("Done. All plots saved in: %s" % OUTPUT_DIR)
    print("")
#endfor

if __name__ == "__main__":
    main()
#endif