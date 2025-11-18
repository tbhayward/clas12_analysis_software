#!/usr/bin/env python3

import os
import sys
import ROOT


def fill_hist_from_tree(tree, hist, branch_name, max_events=None):
    """Loop over a TTree and fill a histogram with -t1.

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
    """
    # Make sure branch exists
    if not hasattr(tree, branch_name):
        raise RuntimeError("Branch '{0}' not found in tree '{1}'".format(
            branch_name, tree.GetName()
        ))
    #endif

    branch = getattr(tree, branch_name)

    n_entries = tree.GetEntries()
    if max_events is not None and max_events > 0:
        n_entries = min(n_entries, max_events)
    #endif

    for i in range(n_entries):
        tree.GetEntry(i)
        value = -float(branch)  # negate t1
        hist.Fill(value)
    #endfor
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

    # Histogram settings
    nbins = 100
    x_min = 0.0
    x_max = 2.0

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

    # Fill histograms (with optional max_events cap)
    fill_hist_from_tree(t_gen_born, h_gen_born, "t1", max_events=max_events)
    fill_hist_from_tree(t_gen_rad,  h_gen_rad,  "t1", max_events=max_events)
    fill_hist_from_tree(t_rec_born, h_rec_born, "t1", max_events=max_events)
    fill_hist_from_tree(t_rec_rad,  h_rec_rad,  "t1", max_events=max_events)

    # Normalize to unit area so shapes can be compared
    for h in (h_gen_born, h_gen_rad, h_rec_born, h_rec_rad):
        integral = h.Integral()
        if integral > 0:
            h.Scale(1.0 / integral)
        #endif
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