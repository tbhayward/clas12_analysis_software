#!/usr/bin/env python3

import os
import sys
import ROOT

def main():
    # Input ROOT file and tree
    infile = "/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi+_10_10_10.root"
    treename = "PhysicsEvents"

    # Output
    outdir = "output"
    outfile = os.path.join(outdir, "quick_mx_check.png")

    # Batch mode so no GUI pops up
    ROOT.gROOT.SetBatch(True)

    # Ensure output directory exists
    os.makedirs(outdir, exist_ok=True)

    # Open file
    f = ROOT.TFile.Open(infile)
    if not f or f.IsZombie():
        print("Error: could not open file:", infile)
        return 1
    #endif

    # Get tree
    tree = f.Get(treename)
    if not tree:
        print("Error: could not find tree '{}' in file".format(treename))
        f.Close()
        return 1
    #endif

    # Define histogram: 100 bins from 0.8 to 1.1
    nbins = 100
    xmin = 0.8
    xmax = 1.1
    hist = ROOT.TH1F("hMx2", ";Mx2 (GeV^{2});Counts", nbins, xmin, xmax)

    # Draw Mx2 into the histogram with cuts:
    # - 0.8 <= Mx2 <= 1.1
    # - detector1 == 1
    # - p_theta (rad) converted to degrees is < 40
    draw_expr = "Mx2>>hMx2"
    cut_expr = (
        "Mx2 >= 0.8 && Mx2 <= 1.1 && "
        "detector1 == 1 && "
        "p_theta < 0.69"
    )
    tree.Draw(draw_expr, cut_expr, "goff")

    # Make canvas and draw histogram
    c = ROOT.TCanvas("c", "Mx2 quick check", 800, 600)
    hist.SetLineWidth(2)
    hist.Draw("HIST")

    # Save to file
    c.SaveAs(outfile)
    print("Saved plot to:", outfile)

    # Clean up
    f.Close()
    return 0
#endif

if __name__ == "__main__":
    sys.exit(main())
#endif