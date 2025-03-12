import numpy as np
import ROOT
import uproot
import math
import argparse
import sys
import os

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description='Compare two ROOT trees.')
    parser.add_argument('--data1', required=True, help='First ROOT tree file')
    parser.add_argument('--data2', required=True, help='Second ROOT tree file')
    parser.add_argument('--type', type=int, choices=[1,2], required=True,
                      help='1: Mx2_1/p1_theta, 2: Mx2_2/p2_theta')
    parser.add_argument('--label1', required=True, help='Legend label for data1')
    parser.add_argument('--label2', required=True, help='Legend label for data2')
    args = parser.parse_args()

    # ROOT configuration
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetPadLeftMargin(0.18)
    ROOT.gStyle.SetPadBottomMargin(0.15)

    # Configuration based on type
    config = {
        1: ('Mx2_1', 'p1_theta', -0.3, 0.3, "M_{x1}^{2} (GeV^{2})", 0.0),
        2: ('Mx2_2', 'p2_theta', 0.578, 1.178, "M_{x2}^{2} (GeV^{2})", 0.880)
    }
    mx_branch, theta_branch, xmin, xmax, xlabel, vline = config[args.type]

    # Data loading
    def load_data(file_path):
        try:
            f = uproot.open(file_path)
            tree = f[f.keys()[0]]
            return (tree[mx_branch].array(library="np"), 
                    tree[theta_branch].array(library="np"))
        except Exception as e:
            print(f"Error loading {file_path}: {str(e)}")
            sys.exit(1)

    mx1, theta1 = load_data(args.data1)
    mx2, theta2 = load_data(args.data2)

    # Output setup
    output_dir = "output/energy_loss_validation"
    os.makedirs(output_dir, exist_ok=True)

    # Angular binning (10 bins from 20-65°)
    angular_bins = [(20 + 4.5*i, 20 + 4.5*(i+1)) for i in range(10)]

    # Create canvas
    canvas = ROOT.TCanvas("canvas", "Comparison", 2400, 1800)
    canvas.Divide(4, 3)

    # Storage containers
    fit_results1 = []
    fit_results2 = []
    bin_centers = []

    # Process integrated plot (pad 1)
    canvas.cd(1)
    h1 = ROOT.TH1D("h1_int", "Integrated", 100, xmin, xmax)
    h2 = ROOT.TH1D("h2_int", "Integrated", 100, xmin, xmax)
    
    for x in mx1: h1.Fill(x)
    for x in mx2: h2.Fill(x)
    
    # Perform fits
    fit1 = h1.Fit("gaus", "SQN")
    fit2 = h2.Fit("gaus", "SQN")
    
    # Configure and draw histograms
    for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20)
        h.GetXaxis().SetTitle(xlabel)
        h.GetYaxis().SetTitle("Counts")
        h.GetXaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleSize(0.06)
        h.Draw("PE")
        
        f = h.GetFunction("gaus")
        if f:
            f.SetLineColor(color)
            f.SetLineStyle(2)
            f.Draw("SAME")

    # Add reference line and legend
    line = ROOT.TLine(vline, 0, vline, h1.GetMaximum()*1.1)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray+2)
    line.Draw()
    
    leg = ROOT.TLegend(0.55, 0.15, 0.85, 0.35)
    leg.SetBorderSize(0)
    if fit1 and fit1.IsValid():
        leg.AddEntry(h1, f"{args.label1}: #mu={fit1.Parameter(1):.3f} #sigma={fit1.Parameter(2):.3f}", "p")
    if fit2 and fit2.IsValid():
        leg.AddEntry(h2, f"{args.label2}: #mu={fit2.Parameter(1):.3f} #sigma={fit2.Parameter(2):.3f}", "p")
    leg.Draw()

    # Process angular bins (pads 2-11)
    for i, (low, high) in enumerate(angular_bins):
        pad_num = i + 2
        canvas.cd(pad_num)
        
        low_rad = math.radians(low)
        high_rad = math.radians(high)
        
        # Create histograms
        h1 = ROOT.TH1D(f"h1_{i}", f"{low:.1f}-{high:.1f}°", 100, xmin, xmax)
        h2 = ROOT.TH1D(f"h2_{i}", "", 100, xmin, xmax)
        
        # Fill histograms
        h1.FillN(len(mx1), mx1, np.where((theta1 >= low_rad) & (theta1 < high_rad), 1, 0))
        h2.FillN(len(mx2), mx2, np.where((theta2 >= low_rad) & (theta2 < high_rad), 1, 0))
        
        # Perform fits
        f1, f2 = None, None
        if h1.GetEntries() > 10:
            fit_result = h1.Fit("gaus", "SQN")
            if fit_result and fit_result.IsValid():
                f1 = h1.GetFunction("gaus")
        if h2.GetEntries() > 10:
            fit_result = h2.Fit("gaus", "SQN")
            if fit_result and fit_result.IsValid():
                f2 = h2.GetFunction("gaus")

        # Store results
        fit_results1.append((f1.GetParameter(1), f1.GetParameter(2)) if f1 else (0, 0))
        fit_results2.append((f2.GetParameter(1), f2.GetParameter(2)) if f2 else (0, 0))
        bin_centers.append((low + high)/2)

        # Draw components
        for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetMarkerStyle(20)
            h.GetXaxis().SetTitle(xlabel)
            h.GetYaxis().SetTitle("Counts")
            h.Draw("PE")
            
            f = h.GetFunction("gaus")
            if f:
                f.SetLineColor(color)
                f.SetLineStyle(2)
                f.Draw("SAME")

        # Add reference line and legend
        line = ROOT.TLine(vline, 0, vline, h1.GetMaximum()*1.1)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw()
        
        leg = ROOT.TLegend(0.55, 0.15, 0.85, 0.35)
        leg.SetBorderSize(0)
        if f1:
            leg.AddEntry(h1, f"{args.label1}: #mu={f1.GetParameter(1):.3f}±{f1.GetParError(1):.3f}", "p")
        if f2:
            leg.AddEntry(h2, f"{args.label2}: #mu={f2.GetParameter(1):.3f}±{f2.GetParError(1):.3f}", "p")
        leg.Draw()

    # Final plot: Mean values (pad 12)
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
        gr.SetLineWidth(2)

    mg.Draw("AP")
    mg.GetXaxis().SetLimits(20, 65)
    
    # Set y-axis range
    y_values = [y for y, _ in fit_results1 + fit_results2 if y != 0]
    if y_values:
        ymin = min(y_values) * 0.9
        ymax = max(y_values) * 1.1
        mg.SetMinimum(ymin)
        mg.SetMaximum(ymax)
    
    leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.85)
    leg.AddEntry(gr1, args.label1, "p")
    leg.AddEntry(gr2, args.label2, "p")
    leg.Draw()

    # Save output
    canvas.SaveAs(f"{output_dir}/comparison_{args.label1}_{args.label2}.png")

if __name__ == "__main__":
    main()