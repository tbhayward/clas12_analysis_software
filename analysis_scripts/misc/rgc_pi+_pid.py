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
#   CD -> track_sector_5, track_chi2_5
#   FD -> track_sector_6, track_chi2_6
DETECTOR_CONFIGS = [
    ("CD", "track_sector_5", "track_chi2_5"),
    ("FD", "track_sector_6", "track_chi2_6"),
]

P_BINS = [
    (0.50, 1.00),
    (1.00, 1.50),
    (1.50, 2.00),
    (2.00, 2.50),
    (2.50, 3.00),
    (3.00, 4.00),
    (4.00, 5.00),
    (5.00, 8.00),
]

CHI2_XMIN = -8.0
CHI2_XMAX = 8.0
CHI2_NBINS = 160

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

def style_hist(hist, x_title, y_title):
    hist.SetLineWidth(2)
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)
    hist.GetXaxis().SetTitleSize(0.055)
    hist.GetYaxis().SetTitleSize(0.055)
    hist.GetXaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetLabelSize(0.045)
    hist.GetYaxis().SetTitleOffset(1.10)
#endfor

def determine_fit_window(hist):
    max_bin = hist.GetMaximumBin()
    peak_x = hist.GetBinCenter(max_bin)
    rms = hist.GetRMS()

    if rms <= 0:
        rms = 1.0
    #endif

    fit_min = peak_x - 1.5 * rms
    fit_max = peak_x + 1.5 * rms

    if fit_min < CHI2_XMIN:
        fit_min = CHI2_XMIN
    #endif

    if fit_max > CHI2_XMAX:
        fit_max = CHI2_XMAX
    #endif

    if fit_max <= fit_min:
        fit_min = CHI2_XMIN
        fit_max = CHI2_XMAX
    #endif

    return fit_min, fit_max
#endfor

def fit_histogram_gaussian(hist, fit_name):
    n_entries = hist.GetEntries()

    if n_entries < 10 or hist.Integral() <= 0:
        return None
    #endif

    fit_min, fit_max = determine_fit_window(hist)

    func = ROOT.TF1(fit_name, "gaus", fit_min, fit_max)

    peak_bin = hist.GetMaximumBin()
    peak_x = hist.GetBinCenter(peak_bin)
    peak_y = hist.GetBinContent(peak_bin)

    initial_sigma = hist.GetRMS()
    if initial_sigma <= 0:
        initial_sigma = 1.0
    #endif

    func.SetParameters(peak_y, peak_x, initial_sigma)

    fit_result = hist.Fit(func, "RQSN")

    if int(fit_result) != 0:
        return None
    #endif

    return func
#endfor

def nice_bin_label(pmin, pmax):
    return "%.2f < p < %.2f (GeV)" % (pmin, pmax)
#endfor

def create_histograms(run_label, det_label):
    hists = []

    h_int = ROOT.TH1D(
        "h_chi2_integrated_%s_%s" % (run_label, det_label),
        "",
        CHI2_NBINS,
        CHI2_XMIN,
        CHI2_XMAX
    )
    hists.append(h_int)

    for i_bin, (pmin, pmax) in enumerate(P_BINS):
        hist = ROOT.TH1D(
            "h_chi2_%s_%s_bin_%d" % (run_label, det_label, i_bin),
            "",
            CHI2_NBINS,
            CHI2_XMIN,
            CHI2_XMAX
        )
        hists.append(hist)
    #endfor

    for hist in hists:
        hist.SetDirectory(0)
    #endfor

    return hists
#endfor

def fill_histograms(tree, run_label, det_label, detector_branch, chi2_branch):
    hists = create_histograms(run_label, det_label)

    p_sum = [0.0 for _ in range(len(P_BINS))]
    p_count = [0 for _ in range(len(P_BINS))]

    n_total = tree.GetEntries()

    print("")
    print("============================================================")
    print("Processing %s %s" % (run_label, det_label))
    print("File: %s" % tree.GetCurrentFile().GetName())
    print("Tree entries: %d" % n_total)
    print("Using detector branch: %s" % detector_branch)
    print("Using chi2 branch:     %s" % chi2_branch)
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
        chi2_val = float(getattr(tree, chi2_branch))

        hists[0].Fill(chi2_val)

        for i_bin, (pmin, pmax) in enumerate(P_BINS):
            if p_val >= pmin and p_val < pmax:
                hists[i_bin + 1].Fill(chi2_val)
                p_sum[i_bin] += p_val
                p_count[i_bin] += 1
                break
            #endif
        #endfor
    #endfor

    p_means = []
    for i_bin in range(len(P_BINS)):
        if p_count[i_bin] > 0:
            p_means.append(p_sum[i_bin] / float(p_count[i_bin]))
        else:
            p_means.append(0.5 * (P_BINS[i_bin][0] + P_BINS[i_bin][1]))
        #endif
    #endfor

    return hists, p_means, p_count
#endfor

def print_fit_summary(run_label, det_label, fit_results):
    print("")
    print("------------------------------------------------------------")
    print("Gaussian fit summary for %s %s" % (run_label, det_label))
    print("------------------------------------------------------------")
    print("%-24s  %-12s  %-12s  %-12s  %-12s  %-12s" % ("Bin", "Entries", "mu", "sigma", "mu err", "sigma err"))

    integrated = fit_results["integrated"]
    if integrated["fit"] is not None:
        print(
            "%-24s  %-12d  %-12.5f  %-12.5f  %-12.5f  %-12.5f" % (
                "Integrated",
                integrated["entries"],
                integrated["mu"],
                integrated["sigma"],
                integrated["mu_err"],
                integrated["sigma_err"]
            )
        )
    else:
        print(
            "%-24s  %-12d  %-12s  %-12s  %-12s  %-12s" % (
                "Integrated",
                integrated["entries"],
                "fit failed",
                "fit failed",
                "fit failed",
                "fit failed"
            )
        )
    #endif

    for result in fit_results["bins"]:
        label = "%.2f-%.2f GeV" % (result["pmin"], result["pmax"])
        if result["fit"] is not None:
            print(
                "%-24s  %-12d  %-12.5f  %-12.5f  %-12.5f  %-12.5f" % (
                    label,
                    result["entries"],
                    result["mu"],
                    result["sigma"],
                    result["mu_err"],
                    result["sigma_err"]
                )
            )
        else:
            print(
                "%-24s  %-12d  %-12s  %-12s  %-12s  %-12s" % (
                    label,
                    result["entries"],
                    "fit failed",
                    "fit failed",
                    "fit failed",
                    "fit failed"
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
        if result["fit"] is None:
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
    graph.GetXaxis().SetTitle("p (GeV)")
    graph.GetYaxis().SetTitle("#mu(track #chi^{2})")
    graph.GetXaxis().CenterTitle(True)
    graph.GetYaxis().CenterTitle(True)
    graph.GetXaxis().SetTitleSize(0.055)
    graph.GetYaxis().SetTitleSize(0.055)
    graph.GetXaxis().SetLabelSize(0.045)
    graph.GetYaxis().SetLabelSize(0.045)
    graph.GetYaxis().SetTitleOffset(1.10)

    return graph
#endfor

def analyze_one_detector(tree, run_label, det_label, detector_branch, chi2_branch):
    hists, p_means, p_counts = fill_histograms(tree, run_label, det_label, detector_branch, chi2_branch)

    fit_results = {
        "integrated": None,
        "bins": []
    }

    h_int = hists[0]
    int_fit = fit_histogram_gaussian(h_int, "f_int_%s_%s" % (run_label, det_label))

    fit_results["integrated"] = {
        "entries": int(h_int.GetEntries()),
        "fit": int_fit,
        "mu": int_fit.GetParameter(1) if int_fit else None,
        "sigma": abs(int_fit.GetParameter(2)) if int_fit else None,
        "mu_err": int_fit.GetParError(1) if int_fit else None,
        "sigma_err": int_fit.GetParError(2) if int_fit else None
    }

    for i_bin, (pmin, pmax) in enumerate(P_BINS):
        hist = hists[i_bin + 1]
        fit = fit_histogram_gaussian(hist, "f_bin_%s_%s_%d" % (run_label, det_label, i_bin))

        result = {
            "pmin": pmin,
            "pmax": pmax,
            "pmean": p_means[i_bin],
            "entries": int(hist.GetEntries()),
            "fit": fit,
            "mu": fit.GetParameter(1) if fit else None,
            "sigma": abs(fit.GetParameter(2)) if fit else None,
            "mu_err": fit.GetParError(1) if fit else None,
            "sigma_err": fit.GetParError(2) if fit else None
        }
        fit_results["bins"].append(result)
    #endfor

    print_fit_summary(run_label, det_label, fit_results)

    canvas = ROOT.TCanvas("c_%s_%s" % (run_label, det_label), "", 2400, 1000)
    canvas.Divide(5, 2)

    title_latex = ROOT.TLatex()
    title_latex.SetNDC(True)
    title_latex.SetTextSize(0.050)

    for ipad in range(1, 11):
        pad = canvas.cd(ipad)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.05)
        pad.SetBottomMargin(0.12)
        pad.SetTopMargin(0.10)
    #endfor

    canvas.cd(1)
    style_hist(h_int, "track #chi^{2}", "Counts")
    h_int.SetTitle("")
    h_int.Draw("hist")
    if int_fit:
        int_fit.SetLineWidth(2)
        int_fit.Draw("same")
    #endif
    title_latex.DrawLatex(0.15, 0.88, "Integrated")

    for i_bin, result in enumerate(fit_results["bins"]):
        canvas.cd(i_bin + 2)
        hist = hists[i_bin + 1]
        style_hist(hist, "track #chi^{2}", "Counts")
        hist.SetTitle("")
        hist.Draw("hist")
        if result["fit"]:
            result["fit"].SetLineWidth(2)
            result["fit"].Draw("same")
        #endif
        label = nice_bin_label(result["pmin"], result["pmax"])
        title_latex.DrawLatex(0.15, 0.88, label)
    #endfor

    canvas.cd(10)
    graph = build_graph_from_results(run_label, det_label, fit_results["bins"])

    frame = ROOT.TH1D(
        "frame_%s_%s" % (run_label, det_label),
        "",
        100,
        0.4,
        8.2
    )
    frame.SetDirectory(0)
    style_hist(frame, "p (GeV)", "#mu(track #chi^{2})")
    frame.SetMinimum(-4.0)
    frame.SetMaximum(4.0)
    frame.SetTitle("")
    frame.Draw()
    graph.Draw("P SAME")
    title_latex.DrawLatex(0.15, 0.88, "#mu(p) with #sigma as y error")

    run_title = ROOT.TLatex()
    run_title.SetNDC(True)
    run_title.SetTextSize(0.030)
    run_title.DrawLatex(0.15, 0.96, "%s %s, particle_pid = 211" % (run_label, det_label))

    out_name = os.path.join(
        OUTPUT_DIR,
        "track_chi2_pip_%s_%s.png" % (run_label, det_label)
    )
    canvas.SaveAs(out_name)

    print("Saved plot: %s" % out_name)

    return canvas, hists, graph
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
            "track_sector_5",
            "track_sector_6",
            "track_chi2_5",
            "track_chi2_6",
        ]
        require_branches(tree, required_branches)

        set_branch_statuses(tree, required_branches)

        for det_label, detector_branch, chi2_branch in DETECTOR_CONFIGS:
            analyze_one_detector(tree, run_label, det_label, detector_branch, chi2_branch)
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