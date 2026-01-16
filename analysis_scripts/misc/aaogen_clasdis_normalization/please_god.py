#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
please_god.py

Usage:
  python please_god.py <file1.root> <file2.root>

Reads TTree "PhysicsEvents" and branch/leaf "Mx2" from two ROOT files and overlays
their histograms on a single canvas.

Output:
  output/this_had_better_work_please_god_let_this_work_im_going_insane.png
"""

import os
import sys

import ROOT


def die(msg):
    print(f"FATAL: {msg}", file=sys.stderr)
    sys.exit(1)


def ensure_outdir(path):
    d = os.path.dirname(path)
    if d and (not os.path.isdir(d)):
        os.makedirs(d, exist_ok=True)
    #endif


def open_root_file(path):
    f = ROOT.TFile.Open(path, "READ")
    if (not f) or f.IsZombie():
        die(f"Could not open ROOT file: {path}")
    #endif
    return f


def get_tree(f, name):
    t = f.Get(name)
    if not t:
        die(f"Could not find TTree '{name}' in file: {f.GetName()}")
    #endif
    return t


def list_names_containing(obj_list, needle_lower):
    out = []
    if not obj_list:
        return out
    #endif
    for i in range(obj_list.GetEntries()):
        o = obj_list.At(i)
        if not o:
            continue
        #endif
        n = o.GetName()
        if needle_lower in n.lower():
            out.append(n)
        #endif
    #endfor
    return out


def require_drawable(tree, varname):
    """
    Ensure 'varname' can be drawn from the tree.
    We check both branches and leaves (some files store leaves without a branch name match).
    """
    br = tree.GetBranch(varname)
    lf = tree.GetLeaf(varname)

    if br or lf:
        return
    #endif

    branch_hits = list_names_containing(tree.GetListOfBranches(), varname.lower())
    leaf_hits   = list_names_containing(tree.GetListOfLeaves(),   varname.lower())

    msg = f"'{varname}' not found as a branch or leaf in tree '{tree.GetName()}'."
    if branch_hits:
        msg += f" Branches containing '{varname}': {', '.join(branch_hits)}."
    else:
        msg += f" Branches containing '{varname}': (none)."
    #endif
    if leaf_hits:
        msg += f" Leaves containing '{varname}': {', '.join(leaf_hits)}."
    else:
        msg += f" Leaves containing '{varname}': (none)."
    #endif

    die(msg)


def print_tree_diagnostics(label, file_path, tree, varname):
    n_entries = int(tree.GetEntries())
    print("")
    print(f"--- Diagnostics: {label} ---")
    print(f"File: {file_path}")
    print(f"Tree: {tree.GetName()}")
    print(f"Entries: {n_entries}")

    if n_entries <= 0:
        print("NOTE: Tree has 0 entries. This will produce an empty histogram.")
        return
    #endif

    # Only attempt min/max if the variable is drawable.
    # TTree::GetMinimum/GetMaximum will return 0 in some failure modes,
    # so we also report whether the branch/leaf exists.
    br = tree.GetBranch(varname)
    lf = tree.GetLeaf(varname)
    print(f"Has branch '{varname}': {'YES' if br else 'NO'}")
    print(f"Has leaf   '{varname}': {'YES' if lf else 'NO'}")

    if br or lf:
        try:
            vmin = float(tree.GetMinimum(varname))
            vmax = float(tree.GetMaximum(varname))
            print(f"{varname} min/max (tree-reported): {vmin} / {vmax}")
        except Exception as e:
            print(f"Could not query {varname} min/max: {e}")
        #endtry
    #endif


def hist_from_tree(tree, hist_name, nbins, xmin, xmax, varname="Mx2"):
    if int(tree.GetEntries()) <= 0:
        die(f"Tree '{tree.GetName()}' has 0 entries; cannot fill histogram '{hist_name}'.")
    #endif

    require_drawable(tree, varname)

    h = ROOT.TH1F(hist_name, "", nbins, xmin, xmax)
    h.Sumw2()
    h.SetDirectory(0)

    # Fill via Draw (fast, avoids manual SetBranchAddress pitfalls).
    # IMPORTANT: Use ">>+" only if you intend to accumulate; here we do NOT.
    n_drawn = tree.Draw(f"{varname}>>{hist_name}", "", "goff")

    # ROOT can return 0 (not < 0) if the expression is invalid or nothing is selected.
    # We treat 0 as a hard error here, and print additional debug below.
    if n_drawn <= 0 or h.GetEntries() <= 0:
        # Print a focused set of debugging hints in the fatal message.
        br = tree.GetBranch(varname)
        lf = tree.GetLeaf(varname)
        die(
            f"Histogram '{hist_name}' has zero entries after Draw.\n"
            f"  tree entries = {int(tree.GetEntries())}\n"
            f"  Draw returned = {int(n_drawn)}\n"
            f"  has branch '{varname}' = {'YES' if br else 'NO'}\n"
            f"  has leaf   '{varname}' = {'YES' if lf else 'NO'}\n"
            f"Most common cause: the file/tree you are pointing at is not the one you think it is "
            f"(e.g., '{varname}' is missing or the tree is empty)."
        )
    #endif

    return h


def normalize_unit_area_in_range(h, xmin, xmax):
    bmin = h.FindBin(xmin + 1e-9)
    bmax = h.FindBin(xmax - 1e-9)
    integral = float(h.Integral(bmin, bmax))
    if integral <= 0.0:
        die("Histogram integral over the plotted range is zero; cannot normalize.")
    #endif
    h.Scale(1.0 / integral)


def main():
    if len(sys.argv) != 3:
        die("Usage: python please_god.py <file1.root> <file2.root>")
    #endif

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]

    out_png = "output/this_had_better_work_please_god_let_this_work_im_going_insane.png"
    ensure_outdir(out_png)

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.SetDefaultSumw2(True)

    f1 = open_root_file(file1_path)
    f2 = open_root_file(file2_path)

    t1 = get_tree(f1, "PhysicsEvents")
    t2 = get_tree(f2, "PhysicsEvents")

    # Print diagnostics up front so you can see immediately if mixed_mc.root is empty/changed.
    print_tree_diagnostics("INPUT 1", file1_path, t1, "Mx2")
    print_tree_diagnostics("INPUT 2", file2_path, t2, "Mx2")

    nbins = 200
    xmin = 0.4
    xmax = 2.0

    h1 = hist_from_tree(t1, "h_mx2_1", nbins, xmin, xmax, varname="Mx2")
    h2 = hist_from_tree(t2, "h_mx2_2", nbins, xmin, xmax, varname="Mx2")

    normalize_unit_area_in_range(h1, xmin, xmax)
    normalize_unit_area_in_range(h2, xmin, xmax)

    h1.SetLineWidth(2)
    h2.SetLineWidth(2)
    h1.SetLineColor(ROOT.kBlack)
    h2.SetLineColor(ROOT.kRed)

    h1.SetTitle("")
    h1.GetXaxis().SetTitle("M_{x}^{2} (GeV^{2})")
    h1.GetYaxis().SetTitle("Normalized counts")
    h1.GetXaxis().CenterTitle(True)
    h1.GetYaxis().CenterTitle(True)

    c = ROOT.TCanvas("c_mx2", "c_mx2", 900, 700)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.12)
    c.SetTopMargin(0.08)

    maxy = max(float(h1.GetMaximum()), float(h2.GetMaximum()))
    h1.SetMaximum(1.15 * maxy)

    h1.Draw("HIST")
    h2.Draw("HIST SAME")

    leg = ROOT.TLegend(0.60, 0.75, 0.93, 0.92)
    leg.SetFillStyle(1001)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetBorderSize(1)
    leg.AddEntry(h1, os.path.basename(file1_path), "l")
    leg.AddEntry(h2, os.path.basename(file2_path), "l")
    leg.Draw()

    c.SaveAs(out_png)

    f1.Close()
    f2.Close()

    print("")
    print(f"Wrote: {out_png}")


if __name__ == "__main__":
    main()
#endif