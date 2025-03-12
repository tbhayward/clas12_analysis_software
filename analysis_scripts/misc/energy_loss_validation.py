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
    ROOT.gStyle.SetPadLeftMargin(0.17)
    ROOT.gStyle.SetPadBottomMargin(0.15)

    # Configuration
    if args.type == 1:
        mx_branch = 'Mx2_1'
        theta_branch = 'p1_theta'
        xmin, xmax = -0.3, 0.3
        xlabel = "M_{x1}^{2} (GeV^{2})"
        vline = 0.0
        hline_y = 0.0
    else:
        mx_branch = 'Mx2_2'
        theta_branch = 'p2_theta' 
        xmin, xmax = 0.578, 1.178
        xlabel = "M_{x2}^{2} (GeV^{2})"
        vline = 0.880
        hline_y = 0.880

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

    angular_bins = [(20 + 4.5*i, 20 + 4.5*(i+1)) for i in range(10)]

    canvas = ROOT.TCanvas("canvas", "Comparison", 2400, 1800)
    canvas.Divide(4, 3)
    canvas.SetMargin(0.08, 0.02, 0.08, 0.1)

    all_objects = {
        'hists': [],
        'lines': [],
        'legends': []
    }

    # Integrated plot (pad 1)
    canvas.cd(1)
    h1 = ROOT.TH1D("h1_int", "Integrated", 100, xmin, xmax)
    h2 = ROOT.TH1D("h2_int", "Integrated", 100, xmin, xmax)
    all_objects['hists'].extend([h1, h2])
    
    for x in mx1: h1.Fill(x)
    for x in mx2: h2.Fill(x)
    
    fit_result1 = h1.Fit("gaus", "SQN")
    fit_result2 = h2.Fit("gaus", "SQN")
    
    for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        h.SetMarkerStyle(20)
        h.GetXaxis().SetTitle(xlabel)
        h.GetYaxis().SetTitle("Counts")
        h.GetXaxis().SetTitleSize(0.06)
        h.GetYaxis().SetTitleSize(0.06)
    
    h1.Draw("PE")
    h2.Draw("PE SAME")
    
    if fit_result1 and fit_result1.IsValid():
        f1 = h1.GetFunction("gaus")
        if f1:
            f1.SetLineColor(ROOT.kBlack)
            f1.SetLineStyle(2)
            f1.Draw("SAME")
    
    if fit_result2 and fit_result2.IsValid():
        f2 = h2.GetFunction("gaus")
        if f2:
            f2.SetLineColor(ROOT.kRed)
            f2.SetLineStyle(2)
            f2.Draw("SAME")
    
    # Force pad update before getting ymax
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    ymax = ROOT.gPad.GetUymax()
    line = ROOT.TLine(vline, 0, vline, ymax)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray+2)
    line.Draw()
    all_objects['lines'].append(line)
    
    # Adjusted legend position and size
    leg = ROOT.TLegend(0.2, 0.15, 0.9, 0.32)
    leg.SetBorderSize(1)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextSize(0.035)
    if fit_result1 and fit_result1.IsValid():
        leg.AddEntry(h1, f"{args.label1}: #mu={fit_result1.Parameter(1):.3f}#pm{fit_result1.ParError(1):.3f}", "p")
    else:
        leg.AddEntry(h1, f"{args.label1}: No fit", "p")
    if fit_result2 and fit_result2.IsValid():
        leg.AddEntry(h2, f"{args.label2}: #mu={fit_result2.Parameter(1):.3f}#pm{fit_result2.ParError(1):.3f}", "p")
    else:
        leg.AddEntry(h2, f"{args.label2}: No fit", "p")
    leg.Draw()
    all_objects['legends'].append(leg)

    fit_results1 = []
    fit_results2 = []
    bin_centers = []

    # Angular bins (pads 2-11)
    for i, (low, high) in enumerate(angular_bins):
        pad_num = i + 2
        canvas.cd(pad_num)
        
        low_rad = math.radians(low)
        high_rad = math.radians(high)
        
        title = f"{low:.1f}-{high:.1f}#circ"
        h1 = ROOT.TH1D(f"h1_{i}", title, 100, xmin, xmax)
        h2 = ROOT.TH1D(f"h2_{i}", "", 100, xmin, xmax)
        all_objects['hists'].extend([h1, h2])
        
        mask1 = (theta1 >= low_rad) & (theta1 < high_rad)
        mask2 = (theta2 >= low_rad) & (theta2 < high_rad)
        for x in mx1[mask1]: h1.Fill(x)
        for x in mx2[mask2]: h2.Fill(x)
        
        fit1_valid = False
        fit2_valid = False
        f1, f2 = None, None
        
        if h1.GetEntries() > 10:
            fit_result1 = h1.Fit("gaus", "SQN")
            if fit_result1 and fit_result1.IsValid():
                fit1_valid = True
                f1 = h1.GetFunction("gaus")
        
        if h2.GetEntries() > 10:
            fit_result2 = h2.Fit("gaus", "SQN")
            if fit_result2 and fit_result2.IsValid():
                fit2_valid = True
                f2 = h2.GetFunction("gaus")
        
        if f1:
            f1.SetLineColor(ROOT.kBlack)
            f1.SetLineStyle(2)
            f1.Draw("SAME")
        if f2:
            f2.SetLineColor(ROOT.kRed)
            f2.SetLineStyle(2)
            f2.Draw("SAME")
        
        for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetMarkerStyle(20)
            h.GetXaxis().SetTitle(xlabel)
            h.GetYaxis().SetTitle("Counts")
            h.GetXaxis().SetTitleSize(0.06)
            h.GetYaxis().SetTitleSize(0.06)
        
        h1.Draw("PE")
        h2.Draw("PE SAME")
        
        # Force pad update before getting ymax
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        ymax = ROOT.gPad.GetUymax()
        line = ROOT.TLine(vline, 0, vline, ymax)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw()
        all_objects['lines'].append(line)
        
        # Adjusted legend position
        leg = ROOT.TLegend(0.2, 0.15, 0.9, 0.32)
        leg.SetBorderSize(1)
        leg.SetFillColor(ROOT.kWhite)
        leg.SetTextSize(0.035)
        if fit1_valid:
            leg.AddEntry(h1, f"{args.label1}: #mu={fit_result1.Parameter(1):.3f}#pm{fit_result1.ParError(1):.3f}", "p")
        else:
            leg.AddEntry(h1, f"{args.label1}: No fit", "p")
        if fit2_valid:
            leg.AddEntry(h2, f"{args.label2}: #mu={fit_result2.Parameter(1):.3f}#pm{fit_result2.ParError(1):.3f}", "p")
        else:
            leg.AddEntry(h2, f"{args.label2}: No fit", "p")
        leg.Draw()
        all_objects['legends'].append(leg)
        
        fit_results1.append((fit_result1.Parameter(1), fit_result1.ParError(1)) if fit1_valid else (0,0))
        fit_results2.append((fit_result2.Parameter(1), fit_result2.ParError(1)) if fit2_valid else (0,0))
        bin_centers.append((low + high)/2)

    # Final plot (pad 12)
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
    
    if args.type == 1:
        mg.SetMinimum(-0.2)
        mg.SetMaximum(0.2)
    else:
        mg.SetMinimum(0.5)
        mg.SetMaximum(1.2)
    
    hline = ROOT.TLine(20, hline_y, 65, hline_y)
    hline.SetLineStyle(2)
    hline.SetLineColor(ROOT.kGray+2)
    hline.Draw()
    all_objects['lines'].append(hline)
    
    leg = ROOT.TLegend(0.2, 0.8, 0.5, 0.9)
    leg.SetBorderSize(1)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextSize(0.045)
    leg.AddEntry(gr1, args.label1, "p")
    leg.AddEntry(gr2, args.label2, "p")
    leg.Draw()
    all_objects['legends'].append(leg)

    canvas.SaveAs(f"{output_dir}/comparison_{args.label1}_{args.label2}.png")

if __name__ == "__main__":
    main()