#!/usr/bin/env python3

import os
import argparse
import ROOT


def main():
    parser = argparse.ArgumentParser(
        description="Make 1x3 2D correlation plots of e_phi, p1_phi, and p2_phi versus phi2, with all angles converted from radians to degrees."
    )
    parser.add_argument(
        "input_root",
        help="Input ROOT file containing the PhysicsEvents tree."
    )
    parser.add_argument(
        "--tree",
        default="PhysicsEvents",
        help="Name of the input tree. Default: PhysicsEvents"
    )
    parser.add_argument(
        "--output",
        default="output/correlation.png",
        help="Output PNG file. Default: output/correlation.png"
    )
    parser.add_argument(
        "--bins-y",
        type=int,
        default=180,
        help="Number of y-axis angle bins. Default: 180"
    )
    parser.add_argument(
        "--bins-phi2",
        type=int,
        default=180,
        help="Number of phi2 bins. Default: 180"
    )
    parser.add_argument(
        "--y-min",
        type=float,
        default=-180.0,
        help="Minimum y-axis angle value in degrees after conversion. Default: -180.0"
    )
    parser.add_argument(
        "--y-max",
        type=float,
        default=180.0,
        help="Maximum y-axis angle value in degrees after conversion. Default: 180.0"
    )
    parser.add_argument(
        "--phi2-min",
        type=float,
        default=-180.0,
        help="Minimum phi2 value in degrees after conversion. Default: -180.0"
    )
    parser.add_argument(
        "--phi2-max",
        type=float,
        default=180.0,
        help="Maximum phi2 value in degrees after conversion. Default: 180.0"
    )

    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)

    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    #endif

    infile = ROOT.TFile.Open(args.input_root, "READ")
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Could not open input ROOT file: {args.input_root}")
    #endif

    tree = infile.Get(args.tree)
    if not tree:
        raise RuntimeError(f"Could not find tree '{args.tree}' in file: {args.input_root}")
    #endif

    required_branches = ["e_phi", "p1_phi", "p2_phi", "phi2"]
    available_branches = {branch.GetName() for branch in tree.GetListOfBranches()}

    for branch_name in required_branches:
        if branch_name not in available_branches:
            raise RuntimeError(
                f"Missing required branch '{branch_name}' in tree '{args.tree}'."
            )
        #endif
    #endfor

    rad_to_deg = "180.0 / TMath::Pi()"

    phi2_deg_expression = f"phi2 * {rad_to_deg}"
    e_phi_deg_expression = f"e_phi * {rad_to_deg}"
    p1_phi_deg_expression = f"p1_phi * {rad_to_deg}"
    p2_phi_deg_expression = f"p2_phi * {rad_to_deg}"

    plots = [
        {
            "hist_name": "h_e_phi_phi2",
            "x_expression": phi2_deg_expression,
            "y_expression": e_phi_deg_expression,
            "title": "e_{#phi} vs #phi_{2};#phi_{2} (deg);e_{#phi} (deg)",
        },
        {
            "hist_name": "h_p1_phi_phi2",
            "x_expression": phi2_deg_expression,
            "y_expression": p1_phi_deg_expression,
            "title": "p1_{#phi} vs #phi_{2};#phi_{2} (deg);p1_{#phi} (deg)",
        },
        {
            "hist_name": "h_p2_phi_phi2",
            "x_expression": phi2_deg_expression,
            "y_expression": p2_phi_deg_expression,
            "title": "p2_{#phi} vs #phi_{2};#phi_{2} (deg);p2_{#phi} (deg)",
        },
    ]

    canvas = ROOT.TCanvas("canvas", "Phi correlations with phi2", 1800, 600)
    canvas.Divide(3, 1)

    histograms = []

    for i, plot in enumerate(plots, start=1):
        canvas.cd(i)

        pad = ROOT.gPad
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.16)
        pad.SetBottomMargin(0.13)
        pad.SetTopMargin(0.10)

        hist = ROOT.TH2D(
            plot["hist_name"],
            plot["title"],
            args.bins_phi2,
            args.phi2_min,
            args.phi2_max,
            args.bins_y,
            args.y_min,
            args.y_max,
        )

        draw_expression = (
            f"({plot['y_expression']}):({plot['x_expression']})>>{plot['hist_name']}"
        )

        cut_expression = (
            f"(({plot['x_expression']}) > {args.phi2_min}) && "
            f"(({plot['x_expression']}) < {args.phi2_max}) && "
            f"(({plot['y_expression']}) > {args.y_min}) && "
            f"(({plot['y_expression']}) < {args.y_max})"
        )

        tree.Draw(draw_expression, cut_expression, "COLZ")

        hist.GetXaxis().SetTitleSize(0.045)
        hist.GetYaxis().SetTitleSize(0.045)
        hist.GetZaxis().SetTitleSize(0.040)

        hist.GetXaxis().SetLabelSize(0.040)
        hist.GetYaxis().SetLabelSize(0.040)
        hist.GetZaxis().SetLabelSize(0.035)

        hist.GetXaxis().SetTitleOffset(1.15)
        hist.GetYaxis().SetTitleOffset(1.35)
        hist.GetZaxis().SetTitleOffset(1.25)

        histograms.append(hist)
    #endfor

    canvas.SaveAs(args.output)

    infile.Close()

    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
#endif