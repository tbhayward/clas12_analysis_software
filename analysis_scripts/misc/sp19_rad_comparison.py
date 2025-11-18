#!/usr/bin/env python3

import os
import sys
import ROOT


def fill_three_hists_from_tree(tree,
                               h_x, h_q2, h_t,
                               max_events=None,
                               label=""):
    """
    Fill x, Q2, and -t1 histograms from a TTree, applying -t1 < 1.0 cut.

    Parameters
    ----------
    tree : ROOT.TTree
        Input tree (PhysicsEvents).
    h_x : ROOT.TH1
        Histogram for x.
    h_q2 : ROOT.TH1
        Histogram for Q2.
    h_t : ROOT.TH1
        Histogram for -t1.
    max_events : int or None
        If not None, cap the number of entries processed.
    label : str
        Label for debug printouts.
    """
    leaf_t1 = tree.GetLeaf("t1")
    leaf_x  = tree.GetLeaf("x")
    leaf_q2 = tree.GetLeaf("Q2")

    if not leaf_t1:
        raise RuntimeError("Leaf 't1' not found in tree '{0}'".format(tree.GetName()))
    #endif
    if not leaf_x:
        raise RuntimeError("Leaf 'x' not found in tree '{0}'".format(tree.GetName()))
    #endif
    if not leaf_q2:
        raise RuntimeError("Leaf 'Q2' not found in tree '{0}'".format(tree.GetName()))
    #endif

    n_entries_total = tree.GetEntries()
    if max_events is not None and max_events > 0:
        n_entries = min(n_entries_total, max_events)
    else:
        n_entries = n_entries_total
    #endif

    print("[{0}] Filling hists from tree '{1}' with up to {2} entries (tree has {3})"
          .format(label, tree.GetName(), n_entries, n_entries_total))

    debug_limit = 10
    debug_triplets = []

    n_kept = 0

    for i in range(n_entries):
        tree.GetEntry(i)

        t_raw = float(leaf_t1.GetValue())
        t_minus = -t_raw    # -t1

        # Apply -t < 1.0 cut, and require -t >= 0
        if t_minus < 0.0 or t_minus >= 1.0:
            continue
        #endif

        x_val  = float(leaf_x.GetValue())
        q2_val = float(leaf_q2.GetValue())

        # Fill histograms (ranges are set by histogram binning)
        h_t.Fill(t_minus)
        h_x.Fill(x_val)
        h_q2.Fill(q2_val)
        n_kept += 1

        if len(debug_triplets) < debug_limit:
            debug_triplets.append((x_val, q2_val, t_minus))
        #endif
    #endfor

    print("[{0}] Kept {1} events after -t<1 cut".format(label, n_kept))
    print("[{0}] First {1} (x, Q2, -t1) triplets after cut: {2}"
          .format(label, len(debug_triplets), debug_triplets))
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
    # gen_rad_path  = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_sp19_inb_10200MeV.root"
    # gen_born_path = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_sp19_inb_10200MeV.root"

    # rec_born_path = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_10200MeV.root"
    # rec_rad_path  = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_sp19_inb_10200MeV.root"

    gen_rad_path  = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/gen_dvcsgen_rad_rga_fa18_inb_10604MeV.root"
    gen_born_path = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/gen_dvcsgen_rga_fa18_inb_10604MeV.root"

    rec_born_path = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_10604MeV.root"
    rec_rad_path  = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen_rad/rec_dvcsgen_rad_rga_fa18_inb_10604MeV.root"

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

    # Histogram settings
    nbins_x  = 100
    nbins_q2 = 100
    nbins_t  = 100

    x_min, x_max     = 0.0, 1.0
    q2_min, q2_max   = 1.0, 8.0
    t_min, t_max     = 0.0, 1.0   # -t1 range

    # Generated MC histograms
    h_gen_born_x  = ROOT.TH1D("h_gen_born_x",  "Generated MC; x; Normalized counts", nbins_x,  x_min,  x_max)
    h_gen_rad_x   = ROOT.TH1D("h_gen_rad_x",   "Generated MC; x; Normalized counts", nbins_x,  x_min,  x_max)

    h_gen_born_q2 = ROOT.TH1D("h_gen_born_q2", "Generated MC; Q^{2} (GeV^{2}); Normalized counts", nbins_q2, q2_min, q2_max)
    h_gen_rad_q2  = ROOT.TH1D("h_gen_rad_q2",  "Generated MC; Q^{2} (GeV^{2}); Normalized counts", nbins_q2, q2_min, q2_max)

    h_gen_born_t  = ROOT.TH1D("h_gen_born_t",  "Generated MC; -t1 (GeV^{2}); Normalized counts", nbins_t,  t_min,  t_max)
    h_gen_rad_t   = ROOT.TH1D("h_gen_rad_t",   "Generated MC; -t1 (GeV^{2}); Normalized counts", nbins_t,  t_min,  t_max)

    # Reconstructed MC histograms
    h_rec_born_x  = ROOT.TH1D("h_rec_born_x",  "Reconstructed MC; x; Normalized counts", nbins_x,  x_min,  x_max)
    h_rec_rad_x   = ROOT.TH1D("h_rec_rad_x",   "Reconstructed MC; x; Normalized counts", nbins_x,  x_min,  x_max)

    h_rec_born_q2 = ROOT.TH1D("h_rec_born_q2", "Reconstructed MC; Q^{2} (GeV^{2}); Normalized counts", nbins_q2, q2_min, q2_max)
    h_rec_rad_q2  = ROOT.TH1D("h_rec_rad_q2",  "Reconstructed MC; Q^{2} (GeV^{2}); Normalized counts", nbins_q2, q2_min, q2_max)

    h_rec_born_t  = ROOT.TH1D("h_rec_born_t",  "Reconstructed MC; -t1 (GeV^{2}); Normalized counts", nbins_t,  t_min,  t_max)
    h_rec_rad_t   = ROOT.TH1D("h_rec_rad_t",   "Reconstructed MC; -t1 (GeV^{2}); Normalized counts", nbins_t,  t_min,  t_max)

    # Style
    all_hists = [
        h_gen_born_x, h_gen_rad_x, h_gen_born_q2, h_gen_rad_q2, h_gen_born_t, h_gen_rad_t,
        h_rec_born_x, h_rec_rad_x, h_rec_born_q2, h_rec_rad_q2, h_rec_born_t, h_rec_rad_t
    ]
    for h in all_hists:
        h.SetStats(0)
    #endfor

    # Colors: Born black, Rad red
    for h in (h_gen_born_x, h_gen_born_q2, h_gen_born_t,
              h_rec_born_x, h_rec_born_q2, h_rec_born_t):
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(2)
    #endfor
    for h in (h_gen_rad_x, h_gen_rad_q2, h_gen_rad_t,
              h_rec_rad_x, h_rec_rad_q2, h_rec_rad_t):
        h.SetLineColor(ROOT.kRed)
        h.SetLineWidth(2)
    #endfor

    # Fill histograms from each tree
    fill_three_hists_from_tree(t_gen_born,
                               h_gen_born_x, h_gen_born_q2, h_gen_born_t,
                               max_events=max_events,
                               label="gen_born")
    fill_three_hists_from_tree(t_gen_rad,
                               h_gen_rad_x, h_gen_rad_q2, h_gen_rad_t,
                               max_events=max_events,
                               label="gen_rad")
    fill_three_hists_from_tree(t_rec_born,
                               h_rec_born_x, h_rec_born_q2, h_rec_born_t,
                               max_events=max_events,
                               label="rec_born")
    fill_three_hists_from_tree(t_rec_rad,
                               h_rec_rad_x, h_rec_rad_q2, h_rec_rad_t,
                               max_events=max_events,
                               label="rec_rad")

    # Summaries before normalization
    for h, label in (
        (h_gen_born_x,  "gen_born_x (before norm)"),
        (h_gen_rad_x,   "gen_rad_x (before norm)"),
        (h_gen_born_q2, "gen_born_q2 (before norm)"),
        (h_gen_rad_q2,  "gen_rad_q2 (before norm)"),
        (h_gen_born_t,  "gen_born_t (before norm)"),
        (h_gen_rad_t,   "gen_rad_t (before norm)"),
        (h_rec_born_x,  "rec_born_x (before norm)"),
        (h_rec_rad_x,   "rec_rad_x (before norm)"),
        (h_rec_born_q2, "rec_born_q2 (before norm)"),
        (h_rec_rad_q2,  "rec_rad_q2 (before norm)"),
        (h_rec_born_t,  "rec_born_t (before norm)"),
        (h_rec_rad_t,   "rec_rad_t (before norm)"),
    ):
        summarize_hist(h, label)
    #endfor

    # Normalize each histogram to unit area
    for h, label in (
        (h_gen_born_x,  "gen_born_x"),
        (h_gen_rad_x,   "gen_rad_x"),
        (h_gen_born_q2, "gen_born_q2"),
        (h_gen_rad_q2,  "gen_rad_q2"),
        (h_gen_born_t,  "gen_born_t"),
        (h_gen_rad_t,   "gen_rad_t"),
        (h_rec_born_x,  "rec_born_x"),
        (h_rec_rad_x,   "rec_rad_x"),
        (h_rec_born_q2, "rec_born_q2"),
        (h_rec_rad_q2,  "rec_rad_q2"),
        (h_rec_born_t,  "rec_born_t"),
        (h_rec_rad_t,   "rec_rad_t"),
    ):
        integral = h.Integral()
        if integral > 0.0:
            h.Scale(1.0 / integral)
        #endif
        summarize_hist(h, label + " (after norm)")
    #endfor

    # Create output directory
    out_dir = "output"
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    #endif

    # Canvas: 3 rows (x, Q2, t), 2 columns (gen, rec)
    c = ROOT.TCanvas("c_sp19_rad_comparison", "sp19_rad_comparison", 1200, 900)
    c.Divide(2, 3)

    # ---------- Row 1: x ----------
    # Left: generated
    c.cd(1)
    ROOT.gPad.SetGrid()
    h_gen_born_x.SetTitle("Generated MC")
    max_x_gen = max(h_gen_born_x.GetMaximum(), h_gen_rad_x.GetMaximum())
    if max_x_gen > 0.0:
        h_gen_born_x.SetMinimum(0.0)
        h_gen_born_x.SetMaximum(1.1 * max_x_gen)
    #endif
    h_gen_born_x.Draw("HIST")
    h_gen_rad_x.Draw("HIST SAME")

    leg_x_gen = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg_x_gen.SetBorderSize(0)
    leg_x_gen.SetFillStyle(0)
    leg_x_gen.AddEntry(h_gen_born_x, "Born", "l")
    leg_x_gen.AddEntry(h_gen_rad_x,  "Rad",  "l")
    leg_x_gen.Draw()

    # Right: reconstructed
    c.cd(2)
    ROOT.gPad.SetGrid()
    h_rec_born_x.SetTitle("Reconstructed MC")
    max_x_rec = max(h_rec_born_x.GetMaximum(), h_rec_rad_x.GetMaximum())
    if max_x_rec > 0.0:
        h_rec_born_x.SetMinimum(0.0)
        h_rec_born_x.SetMaximum(1.1 * max_x_rec)
    #endif
    h_rec_born_x.Draw("HIST")
    h_rec_rad_x.Draw("HIST SAME")

    leg_x_rec = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg_x_rec.SetBorderSize(0)
    leg_x_rec.SetFillStyle(0)
    leg_x_rec.AddEntry(h_rec_born_x, "Born", "l")
    leg_x_rec.AddEntry(h_rec_rad_x,  "Rad",  "l")
    leg_x_rec.Draw()

    # ---------- Row 2: Q2 ----------
    # Left: generated
    c.cd(3)
    ROOT.gPad.SetGrid()
    h_gen_born_q2.SetTitle("Generated MC")
    max_q2_gen = max(h_gen_born_q2.GetMaximum(), h_gen_rad_q2.GetMaximum())
    if max_q2_gen > 0.0:
        h_gen_born_q2.SetMinimum(0.0)
        h_gen_born_q2.SetMaximum(1.1 * max_q2_gen)
    #endif
    h_gen_born_q2.Draw("HIST")
    h_gen_rad_q2.Draw("HIST SAME")

    leg_q2_gen = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg_q2_gen.SetBorderSize(0)
    leg_q2_gen.SetFillStyle(0)
    leg_q2_gen.AddEntry(h_gen_born_q2, "Born", "l")
    leg_q2_gen.AddEntry(h_gen_rad_q2,  "Rad",  "l")
    leg_q2_gen.Draw()

    # Right: reconstructed
    c.cd(4)
    ROOT.gPad.SetGrid()
    h_rec_born_q2.SetTitle("Reconstructed MC")
    max_q2_rec = max(h_rec_born_q2.GetMaximum(), h_rec_rad_q2.GetMaximum())
    if max_q2_rec > 0.0:
        h_rec_born_q2.SetMinimum(0.0)
        h_rec_born_q2.SetMaximum(1.1 * max_q2_rec)
    #endif
    h_rec_born_q2.Draw("HIST")
    h_rec_rad_q2.Draw("HIST SAME")

    leg_q2_rec = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg_q2_rec.SetBorderSize(0)
    leg_q2_rec.SetFillStyle(0)
    leg_q2_rec.AddEntry(h_rec_born_q2, "Born", "l")
    leg_q2_rec.AddEntry(h_rec_rad_q2,  "Rad",  "l")
    leg_q2_rec.Draw()

    # ---------- Row 3: -t1 ----------
    # Left: generated
    c.cd(5)
    ROOT.gPad.SetGrid()
    h_gen_born_t.SetTitle("Generated MC")
    max_t_gen = max(h_gen_born_t.GetMaximum(), h_gen_rad_t.GetMaximum())
    if max_t_gen > 0.0:
        h_gen_born_t.SetMinimum(0.0)
        h_gen_born_t.SetMaximum(1.1 * max_t_gen)
    #endif
    h_gen_born_t.Draw("HIST")
    h_gen_rad_t.Draw("HIST SAME")

    leg_t_gen = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg_t_gen.SetBorderSize(0)
    leg_t_gen.SetFillStyle(0)
    leg_t_gen.AddEntry(h_gen_born_t, "Born", "l")
    leg_t_gen.AddEntry(h_gen_rad_t,  "Rad",  "l")
    leg_t_gen.Draw()

    # Right: reconstructed
    c.cd(6)
    ROOT.gPad.SetGrid()
    h_rec_born_t.SetTitle("Reconstructed MC")
    max_t_rec = max(h_rec_born_t.GetMaximum(), h_rec_rad_t.GetMaximum())
    if max_t_rec > 0.0:
        h_rec_born_t.SetMinimum(0.0)
        h_rec_born_t.SetMaximum(1.1 * max_t_rec)
    #endif
    h_rec_born_t.Draw("HIST")
    h_rec_rad_t.Draw("HIST SAME")

    leg_t_rec = ROOT.TLegend(0.60, 0.70, 0.88, 0.88)
    leg_t_rec.SetBorderSize(0)
    leg_t_rec.SetFillStyle(0)
    leg_t_rec.AddEntry(h_rec_born_t, "Born", "l")
    leg_t_rec.AddEntry(h_rec_rad_t,  "Rad",  "l")
    leg_t_rec.Draw()

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