#!/usr/bin/env python3

import os
import argparse
import ROOT


def main():
    parser = argparse.ArgumentParser(
        description="Make 1x3 2D correlation plots with all angular branches converted from radians to degrees."
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
        "--bins-theta",
        type=int,
        default=180,
        help="Number of theta bins. Default: 180"
    )
    parser.add_argument(
        "--bins-phi",
        type=int,
        default=180,
        help="Number of phi2 bins. Default: 180"
    )
    parser.add_argument(
        "--theta-min",
        type=float,
        default=0.0,
        help="Minimum theta value in degrees after conversion. Default: 0.0"
    )
    parser.add_argument(
        "--theta-max",
        type=float,
        default=60.0,
        help="Maximum theta value in degrees after conversion. Default: 60.0"
    )
    parser.add_argument(
        "--phi-min",
        type=float,
        default=-180.0,
        help="Minimum phi2 value in degrees after conversion. Default: -180.0"
    )
    parser.add_argument(
        "--phi-max",
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

    required_branches = ["e_theta", "p1_theta", "p2_theta", "phi2"]
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
    e_theta_deg_expression = f"e_theta * {rad_to_deg}"
    p1_theta_deg_expression = f"p1_theta * {rad_to_deg}"
    p2_theta_deg_expression = f"p2_theta * {rad_to_deg}"

    plots = [
        {
            "hist_name": "h_e_theta_phi2",
            "x_expression": phi2_deg_expression,
            "y_expression": e_theta_deg_expression,
            "title": "e_{#theta} vs #phi_{2};#phi_{2} (deg);e_{#theta} (deg)",
        },
        {
            "hist_name": "h_p1_theta_phi2",
            "x_expression": phi2_deg_expression,
            "y_expression": p1_theta_deg_expression,
            "title": "p1_{#theta} vs #phi_{2};#phi_{2} (deg);p1_{#theta} (deg)",
        },
        {
            "hist_name": "h_p2_theta_phi2",
            "x_expression": phi2_deg_expression,
            "y_expression": p2_theta_deg_expression,
            "title": "p2_{#theta} vs #phi_{2};#phi_{2} (deg);p2_{#theta} (deg)",
        },
    ]

    canvas = ROOT.TCanvas("canvas", "Theta correlations with phi2", 1800, 600)
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
            args.bins_phi,
            args.phi_min,
            args.phi_max,
            args.bins_theta,
            args.theta_min,
            args.theta_max,
        )

        draw_expression = (
            f"({plot['y_expression']}):({plot['x_expression']})>>{plot['hist_name']}"
        )

        cut_expression = (
            f"(({plot['x_expression']}) > {args.phi_min}) && "
            f"(({plot['x_expression']}) < {args.phi_max}) && "
            f"(({plot['y_expression']}) > {args.theta_min}) && "
            f"(({plot['y_expression']}) < {args.theta_max})"
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