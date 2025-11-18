#!/usr/bin/env python3

import os
import sys
import ROOT


def fill_hist_from_tree(tree, hist, branch_name, max_events=None, label=""):
    """Loop over a TTree and fill a histogram with -t1, with -t<1 cut.

    Parameters
    ----------
    tree : ROOT.TTree
        Input tree.
    hist : ROOT.TH1
        Histogram to fill.
    branch_name : str
        Name of the branch to read (e.g. 't1').
    max_events : int or None
        If not None, cap the number of entries processed at this value.
    label : str
        Label for debug printouts (e.g. 'gen_born').
    """
    # Make sure branch/leaf exists
    leaf = tree.GetLeaf(branch_name)
    if not leaf:
        raise RuntimeError(
            "Leaf/branch '{0}' not found in tree '{1}'".format(
                branch_name, tree.GetName()
            )
        )
    #endif

    n_entries_total = tree.GetEntries()
    if max_events is not None and max_events > 0:
        n_entries = min(n_entries_total, max_events)
    else:
        n_entries = n_entries_total
    #endif

    print("[{0}] Filling histogram '{1}' from branch '{2}' with up to {3} entries "
          "(tree has {4})"
          .format(label, hist.GetName(), branch_name, n_entries, n_entries_total))

    debug_print_limit = 10
    debug_values = []

    vmin = None
    vmax = None
    n_kept = 0

    for i in range(n_entries):
        tree.GetEntry(i)
        # Read the branch value and negate it: we plot -t1
        t_val = float(leaf.GetValue())
        value = -t_val  # this is -t1

        # Apply cut: -t < 1.0 (and also require value >= 0 just in case)
        if value < 0.0 or value >= 1.0:
            continue
        #endif

        hist.Fill(value)
        n_kept += 1

        if len(debug_values) < debug_print_limit:
            debug_values.append(value)
        #endif

        if vmin is None or value < vmin:
            vmin = value
        #endif
        if vmax is None or value > vmax:
            vmax = value
        #endif
    #endfor

    print("[{0}] Kept {1} events after -t<1 cut".format(label, n_kept))
    print("[{0}] First {1} filled -{2} values (after cut): {3}"
          .format(label, len(debug_values), branch_name, debug_values))
    if n_kept > 0:
        print("[{0}] Min(-{1})={2:.6g}, Max(-{1})={3:.6g}"
              .format(label, branch_name, vmin, vmax))
    else:
        print("[{0}] WARNING: no events passed the -t<1 cut".format(label))
    #endif
#enddef


def summarize_hist(hist, label=""):
    """Print a quick summary of a histogram."""
    entries = hist.GetEntries()
    integral = hist.Integral()
    mean = hist.GetMean()
    rms = hist.GetRMS()
    print("[{0}] Hist '{1}': entries={2}, integral={3:.6g}, mean={4:.6g}, rms={5:.6g}"
          .format(label, hist.GetName(), int(entries), integral, mean, rms))
#enddef


def main():
    ROOT.gROOT.SetBatch(True)

    # Optional command-line argument: max_events per tree
    max_events = None
    if len(sys.argv) > 1:
        try:
            max_events = int(sys.argv[1])
        except ValueError:
            raise RuntimeError("First argument must be an integer max_events")
        #endtry
        print("Using max_events = {0} per tree".format(max_events))
    #endif

    # Input file paths
    gen_rad_path  = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp19_inb_10200MeV.root"
    gen_born_path = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_10200MeV.root"

    rec_born_path = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_10200MeV.root"
    rec_rad_path  = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp19_inb_10200MeV.root"

    # Open files
    f_gen_rad  = ROOT.TFile.Open(gen_rad_path, "READ")
    f_gen_born = ROOT.TFile.Open(gen_born_path, "READ")
    f_rec_born = ROOT.TFile.Open(rec_born_path, "READ")
    f_rec_rad  = ROOT.TFile.Open(rec_rad_path, "READ")

    if not f_gen_rad or f_gen_rad.IsZombie():
        raise RuntimeError("Could not open file: {0}".format(gen_rad_path))
    #endif
    if not f_gen_born or f_gen_born.IsZombie():
        raise RuntimeError("Could not open file: {0}".format(gen_born_path))
    #endif
    if not f_rec_born or f_rec_born.IsZombie():
        raise RuntimeError("Could not open file: {0}".format(rec_born_path))
    #endif
    if not f_rec_rad or f_rec_rad.IsZombie():
        raise RuntimeError("Could not open file: {0}".format(rec_rad_path))
    #endif

    # Get trees
    t_gen_rad  = f_gen_rad.Get("PhysicsEvents")
    t_gen_born = f_gen_born.Get("PhysicsEvents")
    t_rec_born = f_rec_born.Get("PhysicsEvents")
    t_rec_rad  = f_rec_rad.Get("PhysicsEvents")

    if not t_gen_rad or not t_gen_born or not t_rec_born or not t_rec_rad:
        raise RuntimeError("One or more 'PhysicsEvents' trees could not be found.")
    #endif

    # Histogram settings: -t in [0,1)
    nbins = 100
    x_min = 0.0
    x_max = 1.0

    # Generated MC histograms
    h_gen_born = ROOT.TH1D("h_gen_born", "Generated MC; -t1 (GeV^{2}); Normalized counts", nbins, x_min, x_max)
    h_gen_rad  = ROOT.TH1D("h_gen_rad",  "Generated MC; -t1 (GeV^{2}); Normalized counts", nbins, x_min, x_max)

    # Reconstructed MC histograms
    h_rec_born = ROOT.TH1D("h_rec_born", "Reconstructed MC; -t1 (GeV^{2}); Normalized counts", nbins, x_min, x_max)
    h_rec_rad  = ROOT.TH1D("h_rec_rad",  "Reconstructed MC; -t1 (GeV^{2}); Normalized counts", nbins, x_min, x_max)

    # Style
    for h in (h_gen_born, h_gen_rad, h_rec_born, h_rec_rad):
        h.SetStats(0)
    #endfor

    h_gen_born.SetLineColor(ROOT.kBlack)
    h_gen_born.SetLineWidth(2)
    h_gen_rad.SetLineColor(ROOT.kRed)
    h_gen_rad.SetLineWidth(2)

    h_rec_born.SetLineColor(ROOT.kBlack)
    h_rec_born.SetLineWidth(2)
    h_rec_rad.SetLineColor(ROOT.kRed)
    h_rec_rad.SetLineWidth(2)

    # Fill histograms (with optional max_events cap), printing debugging info
    fill_hist_from_tree(t_gen_born, h_gen_born, "t1", max_events=max_events, label="gen_born")
    fill_hist_from_tree(t_gen_rad,  h_gen_rad,  "t1", max_events=max_events, label="gen_rad")
    fill_hist_from_tree(t_rec_born, h_rec_born, "t1", max_events=max_events, label="rec_born")
    fill_hist_from_tree(t_rec_rad,  h_rec_rad,  "t1", max_events=max_events, label="rec_rad")

    # Summaries before normalization
    summarize_hist(h_gen_born, "gen_born (before norm)")
    summarize_hist(h_gen_rad,  "gen_rad (before norm)")
    summarize_hist(h_rec_born, "rec_born (before norm)")
    summarize_hist(h_rec_rad,  "rec_rad (before norm)")

    # Normalize to unit area so shapes can be compared
    for h, label in (
        (h_gen_born, "gen_born"),
        (h_gen_rad,  "gen_rad"),
        (h_rec_born, "rec_born"),
        (h_rec_rad,  "rec_rad"),
    ):
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)
        #endif
        summarize_hist(h, label + " (after norm)")
    #endfor

    # Create output directory
    out_dir = "output"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    #endif

    # Canvas: 1x2
    c = ROOT.TCanvas("c_sp19_rad_comparison", "sp19_rad_comparison", 1200, 600)
    c.Divide(2, 1)

    # Left pad: generated MC
    c.cd(1)
    ROOT.gPad.SetGrid()
    h_gen_born.SetTitle("Generated MC")
    # Set y range based on max of the two generated histograms
    gen_max = max(h_gen_born.GetMaximum(), h_gen_rad.GetMaximum())
    if gen_max > 0:
        h_gen_born.SetMinimum(0.0)
        h_gen_born.SetMaximum(1.1 * gen_max)
    #endif
    h_gen_born.Draw("HIST")
    h_gen_rad.Draw("HIST SAME")

    leg1 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.AddEntry(h_gen_born, "Born", "l")
    leg1.AddEntry(h_gen_rad,  "Rad",  "l")
    leg1.Draw()

    # Right pad: reconstructed MC
    c.cd(2)
    ROOT.gPad.SetGrid()
    h_rec_born.SetTitle("Reconstructed MC")
    # Set y range based on max of the two reconstructed histograms
    rec_max = max(h_rec_born.GetMaximum(), h_rec_rad.GetMaximum())
    if rec_max > 0:
        h_rec_born.SetMinimum(0.0)
        h_rec_born.SetMaximum(1.1 * rec_max)
    #endif
    h_rec_born.Draw("HIST")
    h_rec_rad.Draw("HIST SAME")

    leg2 = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.AddEntry(h_rec_born, "Born", "l")
    leg2.AddEntry(h_rec_rad,  "Rad",  "l")
    leg2.Draw()

    # Save output
    out_path = os.path.join(out_dir, "sp19_rad_comparison.png")
    c.SaveAs(out_path)
    print("Saved:", out_path)

    # Close files
    f_gen_rad.Close()
    f_gen_born.Close()
    f_rec_born.Close()
    f_rec_rad.Close()
#enddef


if __name__ == "__main__":
    main()
#endif