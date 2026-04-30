#!/usr/bin/env python3

import os
import sys
import math
import ROOT


def wrap_phi_deg(phi_deg):
    while phi_deg < 0.0:
        phi_deg += 360.0
    #endwhile

    while phi_deg >= 360.0:
        phi_deg -= 360.0
    #endwhile

    return phi_deg


def setup_root_style():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.045, "XYZ")
    ROOT.gStyle.SetLabelSize(0.040, "XYZ")
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.16)
    ROOT.gStyle.SetPadBottomMargin(0.13)


def save_hist2d(hist, output_file_name):
    canvas = ROOT.TCanvas("canvas_" + hist.GetName(), "canvas_" + hist.GetName(), 1000, 800)
    canvas.cd()

    hist.GetXaxis().CenterTitle(True)
    hist.GetYaxis().CenterTitle(True)
    hist.GetZaxis().CenterTitle(True)

    hist.Draw("COLZ")
    canvas.SaveAs(output_file_name)

    print("Saved: {}".format(output_file_name))


def main():
    if len(sys.argv) != 2:
        print("Usage: python plot_2D_proton_theta_comparisons.py input.root")
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

    setup_root_style()

    hist_theta_vs_phi = ROOT.TH2D(
        "hist_p1_theta_vs_phi2",
        ";#phi Trento;proton #theta",
        180,
        0.0,
        360.0,
        140,
        0.0,
        70.0
    )

    hist_theta_vs_x = ROOT.TH2D(
        "hist_p1_theta_vs_x",
        ";x_{B};proton #theta",
        140,
        0.0,
        0.7,
        140,
        0.0,
        70.0
    )

    hist_theta_vs_Q2 = ROOT.TH2D(
        "hist_p1_theta_vs_Q2",
        ";Q^{2} (GeV^{2});proton #theta",
        120,
        1.0,
        7.0,
        140,
        0.0,
        70.0
    )

    hist_theta_vs_minus_t1 = ROOT.TH2D(
        "hist_p1_theta_vs_minus_t1",
        ";-t (GeV^{2});proton #theta",
        100,
        0.0,
        1.0,
        140,
        0.0,
        70.0
    )

    deg = 180.0 / math.pi

    n_entries = tree.GetEntries()

    for i_entry in range(n_entries):
        tree.GetEntry(i_entry)

        p1_theta_deg = tree.p1_theta * deg
        phi2_deg = wrap_phi_deg(tree.phi2 * deg)
        x_value = tree.x
        Q2_value = tree.Q2
        minus_t1_value = -tree.t1

        if p1_theta_deg < 0.0 or p1_theta_deg > 70.0:
            continue
        #endif

        hist_theta_vs_phi.Fill(phi2_deg, p1_theta_deg)

        if x_value >= 0.0 and x_value <= 0.7:
            hist_theta_vs_x.Fill(x_value, p1_theta_deg)
        #endif

        if Q2_value >= 1.0 and Q2_value <= 7.0:
            hist_theta_vs_Q2.Fill(Q2_value, p1_theta_deg)
        #endif

        if minus_t1_value >= 0.0 and minus_t1_value <= 1.0:
            hist_theta_vs_minus_t1.Fill(minus_t1_value, p1_theta_deg)
        #endif
    #endfor

    save_hist2d(hist_theta_vs_phi, "output/2D_proton_theta_vs_phi.pdf")
    save_hist2d(hist_theta_vs_x, "output/2D_proton_theta_vs_xB.pdf")
    save_hist2d(hist_theta_vs_Q2, "output/2D_proton_theta_vs_Q2.pdf")
    save_hist2d(hist_theta_vs_minus_t1, "output/2D_proton_theta_vs_minus_t.pdf")

    root_file.Close()


if __name__ == "__main__":
    main()
#endif