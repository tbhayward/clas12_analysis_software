#!/usr/bin/env python3

import os
import ROOT

def category_for_parent(p1_parent):
    """
    Map mc_p1_parent PDG codes to the categories we care about.
    Returns a string label or None if we want to ignore this parent.
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
    return None
#endif

def main():
    # Run in batch mode (no GUI)
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)  # remove stat box globally

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
    xmin = 0.8
    xmax = 2.4

    # Total Mh histogram
    h_total = ROOT.TH1F(
        "h_mh_total",
        "k^{+} parent;Mh (GeV);Counts",
        nbins,
        xmin,
        xmax
    )
    h_total.SetStats(False)

    # One histogram per category (string, K*0, K*+, phi, exclusive)
    hists_by_category = {}

    # Colors to cycle through for different categories
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

        # Determine category for this parent
        category = category_for_parent(p1_parent)

        # Always fill total histogram
        h_total.Fill(mh)

        # Skip if this parent is not one of the specified categories
        if category is None:
            continue
        #endif

        # Lazily create a histogram for this category
        if category not in hists_by_category:
            idx = len(hists_by_category)
            hist_name = "h_mh_cat_{}".format(category.replace("^", "").replace("{", "").replace("}", ""))
            hist_title = "k^{+} parent;Mh (GeV);Counts"
            h = ROOT.TH1F(hist_name, hist_title, nbins, xmin, xmax)
            h.SetStats(False)

            color = color_list[idx % len(color_list)]
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetLineWidth(2)

            hists_by_category[category] = h
        #endif

        # Fill category histogram
        hists_by_category[category].Fill(mh)
    #endfor

    # Make sure output directory exists
    os.makedirs("output", exist_ok=True)

    # Canvas setup
    canvas = ROOT.TCanvas("c_mh", "K^{+} parent", 800, 600)
    canvas.SetLeftMargin(0.125)
    canvas.cd()

    # Style for total histogram
    h_total.SetLineColor(ROOT.kBlack)
    h_total.SetLineWidth(2)

    # Draw total first
    h_total.Draw("HIST")

    # Track max for y-axis scaling
    max_y = h_total.GetMaximum()

    # Draw each category histogram
    for category, hist in hists_by_category.items():
        hist.Draw("HIST SAME")
        if hist.GetMaximum() > max_y:
            max_y = hist.GetMaximum()
        #endif
    #endfor

    # Rescale y-axis to fit everything nicely
    h_total.SetMaximum(1.2 * max_y if max_y > 0 else 1.0)

    # Legend with a box around it
    legend = ROOT.TLegend(0.60, 0.60, 0.88, 0.88)
    legend.SetBorderSize(1)       # box border
    legend.SetFillStyle(1001)     # solid fill
    legend.SetFillColor(ROOT.kWhite)

    legend.AddEntry(h_total, "All parents", "l")

    # Fixed order for legend entries
    label_order = ["string", "K^{*0}", "K^{*+}", "#phi", "exclusive"]
    for label in label_order:
        if label in hists_by_category:
            legend.AddEntry(hists_by_category[label], label, "l")
        #endif
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