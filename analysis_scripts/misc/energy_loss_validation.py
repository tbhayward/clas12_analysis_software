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

    # Suppress ROOT graphics window
    ROOT.gROOT.SetBatch(True)

    # Set ROOT styles
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)

    # Determine branches and ranges
    if args.type == 1:
        mx_branch = 'Mx2_1'
        theta_branch = 'p1_theta'
        xmin, xmax = -0.3, 0.3
    else:
        mx_branch = 'Mx2_2'
        theta_branch = 'p2_theta' 
        xmin, xmax = 0.578, 1.178

    # Load data
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

    # Create output directory
    os.makedirs('output', exist_ok=True)

    # Define angular bins in degrees
    angular_bins = [(5,11), (11,17), (17,23), (23,29), (29,35),
                    (35,41), (41,47), (47,53), (53,59), (59,65)]

    # Create canvas and pads
    canvas = ROOT.TCanvas("canvas", "Comparison", 2400, 1800)
    canvas.Divide(4, 3)

    fit_results1 = []
    fit_results2 = []
    bin_centers = []

    # Process integrated plot (pad 1)
    canvas.cd(1)
    h1 = ROOT.TH1D("h1_int", "Integrated", 100, xmin, xmax)
    h2 = ROOT.TH1D("h2_int", "Integrated", 100, xmin, xmax)
    
    for x in mx1: h1.Fill(x)
    for x in mx2: h2.Fill(x)
    
    fit1 = h1.Fit("gaus", "SQN")
    fit2 = h2.Fit("gaus", "SQN")
    
    # Configure histograms
    for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20)
        h.GetXaxis().SetRangeUser(xmin, xmax)
    
    h1.Draw("PE")
    h2.Draw("PE SAME")
    
    leg = ROOT.TLegend(0.5, 0.65, 0.9, 0.85)
    if fit1 and fit1.IsValid():
        leg.AddEntry(h1, f"{args.label1} #mu={fit1.Parameter(1):.3f} #sigma={fit1.Parameter(2):.3f}", "p")
    else:
        leg.AddEntry(h1, f"{args.label1} (fit failed)", "p")
    
    if fit2 and fit2.IsValid():
        leg.AddEntry(h2, f"{args.label2} #mu={fit2.Parameter(1):.3f} #sigma={fit2.Parameter(2):.3f}", "p")
    else:
        leg.AddEntry(h2, f"{args.label2} (fit failed)", "p")
    leg.Draw()

    # Process angular bins (pads 2-11)
    for i, (low, high) in enumerate(angular_bins):
        pad_num = i + 2
        canvas.cd(pad_num)
        
        # Convert to radians for data selection
        low_rad = math.radians(low)
        high_rad = math.radians(high)
        
        mask1 = (theta1 >= low_rad) & (theta1 < high_rad)
        mask2 = (theta2 >= low_rad) & (theta2 < high_rad)
        
        # Create histograms with fixed ranges
        h1 = ROOT.TH1D(f"h1_{i}", f"{low}-{high} deg", 100, xmin, xmax)
        h2 = ROOT.TH1D(f"h2_{i}", "", 100, xmin, xmax)
        
        for x in mx1[mask1]: h1.Fill(x)
        for x in mx2[mask2]: h2.Fill(x)
        
        # Perform fits only if there are entries
        fit1_valid = False
        if h1.GetEntries() > 0:
            fit1 = h1.Fit("gaus", "SQN")
            fit1_valid = fit1 and fit1.IsValid() and fit1.Ndf() > 0
        
        fit2_valid = False
        if h2.GetEntries() > 0:
            fit2 = h2.Fit("gaus", "SQN")
            fit2_valid = fit2 and fit2.IsValid() and fit2.Ndf() > 0
        
        # Store fit results
        if fit1_valid:
            fit_results1.append((fit1.Parameter(1), fit1.ParError(1)))
        else:
            fit_results1.append((0, 0))  # Placeholder for failed fits
            
        if fit2_valid:
            fit_results2.append((fit2.Parameter(1), fit2.ParError(1)))
        else:
            fit_results2.append((0, 0))  # Placeholder for failed fits
            
        bin_centers.append((low + high)/2)
        
        # Configure histograms
        h1.SetLineColor(ROOT.kBlack)
        h1.SetMarkerColor(ROOT.kBlack)
        h1.SetMarkerStyle(20)
        h2.SetLineColor(ROOT.kRed)
        h2.SetMarkerColor(ROOT.kRed)
        h2.SetMarkerStyle(20)
        
        h1.Draw("PE")
        h2.Draw("PE SAME")
        
        # Add legend
        leg = ROOT.TLegend(0.5, 0.65, 0.9, 0.85)
        if fit1_valid:
            leg.AddEntry(h1, f"{args.label1} #mu={fit1.Parameter(1):.3f}", "p")
        else:
            leg.AddEntry(h1, f"{args.label1} (no fit)", "p")
            
        if fit2_valid:
            leg.AddEntry(h2, f"{args.label2} #mu={fit2.Parameter(1):.3f}", "p")
        else:
            leg.AddEntry(h2, f"{args.label2} (no fit)", "p")
        leg.Draw()

    # Process mean plot (pad 12)
    canvas.cd(12)
    gr1 = ROOT.TGraphErrors(len(bin_centers))
    gr2 = ROOT.TGraphErrors(len(bin_centers))
    
    for i in range(len(bin_centers)):
        gr1.SetPoint(i, bin_centers[i], fit_results1[i][0])
        gr1.SetPointError(i, 0, fit_results1[i][1])
        gr2.SetPoint(i, bin_centers[i], fit_results2[i][0])
        gr2.SetPointError(i, 0, fit_results2[i][1])
    
    mg = ROOT.TMultiGraph()
    gr1.SetMarkerColor(ROOT.kBlack)
    gr1.SetMarkerStyle(20)
    gr2.SetMarkerColor(ROOT.kRed)
    gr2.SetMarkerStyle(20)
    mg.Add(gr1)
    mg.Add(gr2)
    
    mg.Draw("AP")
    mg.GetXaxis().SetTitle("#theta (degrees)")
    mg.GetYaxis().SetTitle("Mean Mx2 [GeV/c^{2}]")
    mg.GetXaxis().SetLimits(0, 70)
    
    # Set y-axis range based on type
    if args.type == 1:
        mg.SetMinimum(-0.2)
        mg.SetMaximum(0.2)
    else:
        mg.SetMinimum(0.5)
        mg.SetMaximum(1.2)
    
    leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.85)
    leg.AddEntry(gr1, args.label1, "p")
    leg.AddEntry(gr2, args.label2, "p")
    leg.Draw()

    # Save to output directory
    canvas.SaveAs(f"output/comparison_{args.label1}_{args.label2}.png")

if __name__ == "__main__":
    main()