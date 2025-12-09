#!/usr/bin/env python3

import os
import ROOT

def main():
    # Run in batch mode (no GUI)
    ROOT.gROOT.SetBatch(True)

    input_file = "/volatile/clas12/thayward/clasdis_rga_fa18_inb_ek+k-.root"
    tree_name = "PhysicsEvents"

    # Open file
    f = ROOT.TFile.Open(input_file)
    if not f or f.IsZombie():
        print("Error: could not open file {}".format(input_file))
        return
    #endif

    tree = f.Get(tree_name)
    if not tree:
        print("Error: tree {} not found in file".format(tree_name))
        return
    #endif

    n_entries = tree.GetEntries()
    print("Total entries in tree: {}".format(n_entries))

    # Histogram settings
    nbins = 100
    xmin = 0.3
    xmax = 2.0

    # Total Mx histogram
    h_total = ROOT.TH1F(
        "h_mx_total",
        "Mx for selected events;Mx (GeV);Counts",
        nbins,
        xmin,
        xmax
    )

    # One histogram per unique (mc_p1_parent, mc_p2_parent) pair
    hists_by_parent = {}

    # Colors to cycle through for different parent combinations
    color_list = [
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kCyan + 1,
        ROOT.kViolet + 1,
        ROOT.kAzure + 2,
    ]

    # Track unique parent pairs so we only print each once
    seen_pairs = set()

    # Loop over events
    for entry in tree:
        # Apply kaon PID filter
        if entry.matching_p1_pid != 321:
            continue
        #endif
        if entry.matching_p2_pid != -321:
            continue
        #endif

        # Grab Mx
        mx = float(entry.Mx)

        # Parent IDs
        p1_parent = int(entry.mc_p1_parent)
        p2_parent = int(entry.mc_p2_parent)
        pair = (p1_parent, p2_parent)

        # Print each unique pair once, in the order encountered
        if pair not in seen_pairs:
            print(
                "Found new parent pair: mc_p1_parent = {}, mc_p2_parent = {}".format(
                    p1_parent,
                    p2_parent
                )
            )
            seen_pairs.add(pair)
        #endif

        # Lazily create a histogram for this (parent1, parent2) pair
        if pair not in hists_by_parent:
            idx = len(hists_by_parent)
            hist_name = "h_mx_parent_{}_{}".format(p1_parent, p2_parent)
            hist_title = (
                "Mx;Mx (GeV);Counts"
            )
            h = ROOT.TH1F(hist_name, hist_title, nbins, xmin, xmax)

            color = color_list[idx % len(color_list)]
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetLineWidth(2)

            hists_by_parent[pair] = h
        #endif

        # Fill histograms
        h_total.Fill(mx)
        hists_by_parent[pair].Fill(mx)
    #endfor

    # Make sure output directory exists
    os.makedirs("output", exist_ok=True)

    # Canvas setup
    canvas = ROOT.TCanvas("c_mx", "Mx dikaon parents", 800, 600)
    canvas.SetLeftMargin(0.125)
    canvas.cd()

    # Style for total histogram
    h_total.SetLineColor(ROOT.kBlack)
    h_total.SetLineWidth(2)

    # Draw total first
    h_total.Draw("HIST")

    # Track max for y-axis scaling
    max_y = h_total.GetMaximum()

    # Draw each parent-combination histogram
    for pair, hist in hists_by_parent.items():
        hist.Draw("HIST SAME")
        if hist.GetMaximum() > max_y:
            max_y = hist.GetMaximum()
        #endif
    #endfor

    # Rescale y-axis to fit everything nicely
    h_total.SetMaximum(1.2 * max_y if max_y > 0 else 1.0)

    # Legend
    legend = ROOT.TLegend(0.60, 0.60, 0.88, 0.88)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)

    legend.AddEntry(h_total, "All parents", "l")
    for pair, hist in hists_by_parent.items():
        label = "mc_p1_parent={}, mc_p2_parent={}".format(pair[0], pair[1])
        legend.AddEntry(hist, label, "l")
    #endfor

    legend.Draw()

    # Save canvas
    out_name = "output/mx_dikaon_parents.png"
    canvas.SaveAs(out_name)
    print("Saved canvas to {}".format(out_name))

    f.Close()
#endfor

if __name__ == "__main__":
    main()
#endif