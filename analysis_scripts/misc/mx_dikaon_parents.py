#!/usr/bin/env python3

import os
import ROOT

def parent_label(p1_parent):
    """
    Map mc_p1_parent PDG codes to human-readable legend labels.
    """
    if p1_parent == 91 or p1_parent == 92:
        return "string"
    #endif
    if p1_parent == 313:
        return "K^{*0}"
    #endif
    if p1_parent == 323:
        return "K^{*+}"
    #endif
    if p1_parent == 333:
        return "#phi"
    #endif
    if p1_parent == 2212:
        return "exclusive"
    #endif
    return "mc_p1_parent={}".format(p1_parent)
#endif

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

    # Total Mh histogram
    h_total = ROOT.TH1F(
        "h_mh_total",
        "Mh for selected events;Mh (GeV);Counts",
        nbins,
        xmin,
        xmax
    )

    # One histogram per unique mc_p1_parent
    hists_by_p1_parent = {}

    # Colors to cycle through for different parent IDs
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

    # Track unique mc_p1_parent so we only print each once
    seen_p1_parents = set()

    # Loop over events
    for entry in tree:
        # Apply kaon PID filter
        if entry.matching_p1_pid != 321:
            continue
        #endif
        if entry.matching_p2_pid != -321:
            continue
        #endif

        # Grab Mh
        mh = float(entry.Mh)

        # Parent ID of particle 1
        p1_parent = int(entry.mc_p1_parent)

        # Print each unique mc_p1_parent once, in the order encountered
        if p1_parent not in seen_p1_parents:
            print("Found new mc_p1_parent: {}".format(p1_parent))
            seen_p1_parents.add(p1_parent)
        #endif

        # Lazily create a histogram for this mc_p1_parent
        if p1_parent not in hists_by_p1_parent:
            idx = len(hists_by_p1_parent)
            hist_name = "h_mh_p1parent_{}".format(p1_parent)
            hist_title = "Mh;Mh (GeV);Counts"
            h = ROOT.TH1F(hist_name, hist_title, nbins, xmin, xmax)

            color = color_list[idx % len(color_list)]
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetLineWidth(2)

            hists_by_p1_parent[p1_parent] = h
        #endif

        # Fill histograms
        h_total.Fill(mh)
        hists_by_p1_parent[p1_parent].Fill(mh)
    #endfor

    # Make sure output directory exists
    os.makedirs("output", exist_ok=True)

    # Canvas setup
    canvas = ROOT.TCanvas("c_mh", "Mh dikaon parents (by mc_p1_parent)", 800, 600)
    canvas.SetLeftMargin(0.125)
    canvas.cd()

    # Style for total histogram
    h_total.SetLineColor(ROOT.kBlack)
    h_total.SetLineWidth(2)

    # Draw total first
    h_total.Draw("HIST")

    # Track max for y-axis scaling
    max_y = h_total.GetMaximum()

    # Draw each mc_p1_parent histogram
    for p1_parent, hist in hists_by_p1_parent.items():
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
    for p1_parent, hist in hists_by_p1_parent.items():
        label = parent_label(p1_parent)
        legend.AddEntry(hist, label, "l")
    #endfor

    legend.Draw()

    # Save canvas
    out_name = "output/mh_dikaon_parents.png"
    canvas.SaveAs(out_name)
    print("Saved canvas to {}".format(out_name))

    f.Close()
#endif

if __name__ == "__main__":
    main()
#endif