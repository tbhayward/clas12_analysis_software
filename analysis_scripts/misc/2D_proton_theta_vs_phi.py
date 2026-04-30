#!/usr/bin/env python3

import os
import sys
import math
import ROOT


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_2D_proton_theta_vs_phi.py input.root")
        sys.exit(1)
    #endif

    input_file_name = sys.argv[1]

    if not os.path.isfile(input_file_name):
        print("ERROR: input file does not exist: {}".format(input_file_name))
        sys.exit(1)
    #endif

    root_file = ROOT.TFile.Open(input_file_name, "READ")
    if not root_file or root_file.IsZombie():
        print("ERROR: could not open ROOT file: {}".format(input_file_name))
        sys.exit(1)
    #endif

    tree = root_file.Get("PhysicsEvents")
    if not tree:
        print("ERROR: could not find tree named PhysicsEvents in file: {}".format(input_file_name))
        root_file.Close()
        sys.exit(1)
    #endif

    os.makedirs("output", exist_ok=True)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.045, "XYZ")
    ROOT.gStyle.SetLabelSize(0.040, "XYZ")
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.16)
    ROOT.gStyle.SetPadBottomMargin(0.13)

    hist = ROOT.TH2D(
        "hist_p1_theta_vs_phi2",
        ";#phi Trento;proton #theta",
        180,
        0.0,
        360.0,
        140,
        0.0,
        70.0
    )

    deg = 180.0 / math.pi

    n_entries = tree.GetEntries()
    for i_entry in range(n_entries):
        tree.GetEntry(i_entry)

        p1_theta_deg = tree.p1_theta * deg
        phi2_deg = tree.phi2 * deg

        while phi2_deg < 0.0:
            phi2_deg += 360.0
        #endwhile

        while phi2_deg >= 360.0:
            phi2_deg -= 360.0
        #endwhile

        if p1_theta_deg < 0.0 or p1_theta_deg > 70.0:
            continue
        #endif

        hist.Fill(phi2_deg, p1_theta_deg)
    #endfor

    canvas = ROOT.TCanvas("canvas", "canvas", 1000, 800)
    canvas.cd()

    hist.GetXaxis().SetTitle("#phi Trento")
    hist.GetYaxis().SetTitle("proton #theta")
    hist.GetZaxis().SetTitle("Counts")

    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)
    hist.GetZaxis().CenterTitle(True)

    hist.Draw("COLZ")

    output_file_name = "output/2D_proton_theta_vs_phi.pdf"
    canvas.SaveAs(output_file_name)

    root_file.Close()

    print("Saved: {}".format(output_file_name))


if __name__ == "__main__":
    main()
#endif