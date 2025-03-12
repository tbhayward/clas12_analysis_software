import ROOT
import uproot
import numpy as np
import math
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description='Compare two ROOT trees.')
    parser.add_argument('--data1', required=True, help='First ROOT tree file')
    parser.add_argument('--data2', required=True, help='Second ROOT tree file')
    parser.add_argument('--type', type=int, choices=[1,2], required=True,
                        help='1: Use Mx2_1 and p1_theta, 2: Use Mx2_2 and p2_theta')
    parser.add_argument('--label1', required=True, help='Legend label for data1')
    parser.add_argument('--label2', required=True, help='Legend label for data2')
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadBottomMargin(0.12)

    # Configuration based on type
    if args.type == 1:
        mx_branch = 'Mx2_1'
        theta_branch = 'p1_theta'
        xmin, xmax = -0.3, 0.3
        xlabel = "M_{x1}^{2} (GeV^{2})"
        vline = 0.0
    else:
        mx_branch = 'Mx2_2'
        theta_branch = 'p2_theta'
        xmin, xmax = 0.578, 1.178
        xlabel = "M_{x2}^{2} (GeV^{2})"
        vline = 0.880

    def load_data(file_path):
        try:
            f = uproot.open(file_path)
            tree = f[f.keys()[0]]
            return tree[mx_branch].array(library="np"), tree[theta_branch].array(library="np")
        except Exception as e:
            print(f"Error loading {file_path}: {str(e)}")
            sys.exit(1)

    mx1, theta1 = load_data(args.data1)
    mx2, theta2 = load_data(args.data2)

    output_dir = "output/energy_loss_validation"
    os.makedirs(output_dir, exist_ok=True)

    # Equal 5° bins from 20° to 60°
    angular_bins = [(20,25), (25,30), (30,35), (35,40),
                   (40,45), (45,50), (50,55), (55,60)]
    
    canvas = ROOT.TCanvas("canvas", "Comparison", 2400, 1800)
    canvas.Divide(4, 3)
    canvas.SetMargin(0.08, 0.02, 0.08, 0.1)

    # Storage lists to prevent garbage collection
    all_histos = []
    all_lines = []
    fit_results1 = []
    fit_results2 = []
    bin_centers = []

    # Integrated plot (pad 1)
    canvas.cd(1)
    h1 = ROOT.TH1D("h1_int", "Integrated", 100, xmin, xmax)
    h2 = ROOT.TH1D("h2_int", "Integrated", 100, xmin, xmax)
    all_histos.extend([h1, h2])
    
    for x in mx1: h1.Fill(x)
    for x in mx2: h2.Fill(x)
    
    # Process angular bins (pads 2-9)
    for i, (low, high) in enumerate(angular_bins):
        pad_num = i + 2
        canvas.cd(pad_num)
        
        low_rad = math.radians(low)
        high_rad = math.radians(high)
        
        mask1 = (theta1 >= low_rad) & (theta1 < high_rad)
        mask2 = (theta2 >= low_rad) & (theta2 < high_rad)
        
        # Create and store histograms
        h1 = ROOT.TH1D(f"h1_{i}", f"#theta [{low}, {high})", 100, xmin, xmax)
        h2 = ROOT.TH1D(f"h2_{i}", "", 100, xmin, xmax)
        all_histos.extend([h1, h2])
        
        # Fill histograms
        for x in mx1[mask1]: h1.Fill(x)
        for x in mx2[mask2]: h2.Fill(x)
        
        # Perform and store fits
        fit1_valid = False
        if h1.GetEntries() > 10:  # Minimum entries requirement
            fit1 = h1.Fit("gaus", "SQN")
            fit1_valid = fit1.IsValid() if fit1 else False
        fit2_valid = False
        if h2.GetEntries() > 10:
            fit2 = h2.Fit("gaus", "SQN")
            fit2_valid = fit2.IsValid() if fit2 else False

        f1 = h1.GetFunction("gaus") if fit1_valid else None
        f2 = h2.GetFunction("gaus") if fit2_valid else None
        
        # Store results
        fit_results1.append((f1.GetParameter(1), f1.GetParError(1)) if f1 else (0,0))
        fit_results2.append((f2.GetParameter(1), f2.GetParError(1)) if f2 else (0,0))
        bin_centers.append((low + high)/2)

        # Draw components
        for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetMarkerStyle(20)
            h.GetXaxis().SetTitle(xlabel)
            h.GetXaxis().SetTitleSize(0.06)
            h.Draw("PE")
            
            f = h.GetFunction("gaus")
            if f:
                f.SetLineColor(color)
                f.SetLineStyle(2)
                f.Draw("SAME")

        # Add vertical line
        ymax = max(h1.GetMaximum(), h2.GetMaximum()) * 1.1
        line = ROOT.TLine(vline, 0, vline, ymax)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw()
        all_lines.append(line)

        # Legend with stored fit references
        leg = ROOT.TLegend(0.15, 0.7, 0.55, 0.88)
        leg.SetBorderSize(0)
        if f1:
            leg.AddEntry(h1, f"{args.label1} #mu={f1.GetParameter(1):.3f}", "LP")
        if f2:
            leg.AddEntry(h2, f"{args.label2} #mu={f2.GetParameter(1):.3f}", "LP")
        leg.Draw()

    # Mean plot (pad 12)
    canvas.cd(12)
    gr1 = ROOT.TGraphErrors(len(bin_centers))
    gr2 = ROOT.TGraphErrors(len(bin_centers))
    
    for i in range(len(bin_centers)):
        gr1.SetPoint(i, bin_centers[i], fit_results1[i][0])
        gr1.SetPointError(i, 0, fit_results1[i][1])
        gr2.SetPoint(i, bin_centers[i], fit_results2[i][0])
        gr2.SetPointError(i, 0, fit_results2[i][1])
    
    mg = ROOT.TMultiGraph()
    mg.Add(gr1)
    mg.Add(gr2)
    mg.SetTitle(";#theta (degrees);Mean " + xlabel)
    
    for gr, color in [(gr1, ROOT.kBlack), (gr2, ROOT.kRed)]:
        gr.SetMarkerColor(color)
        gr.SetLineColor(color)
        gr.SetMarkerStyle(20)

    mg.Draw("AP")
    mg.GetXaxis().SetLimits(15, 65)
    mg.SetMinimum(xmin if args.type == 1 else 0.5)
    mg.SetMaximum(xmax if args.type == 1 else 1.2)
    
    leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.85)
    leg.AddEntry(gr1, args.label1, "P")
    leg.AddEntry(gr2, args.label2, "P")
    leg.Draw()

    canvas.SaveAs(f"{output_dir}/comparison_{args.label1}_{args.label2}.png")

if __name__ == "__main__":
    main()