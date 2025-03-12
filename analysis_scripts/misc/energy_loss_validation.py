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
            
            # Check for theta cut branch
            available_branches = tree.keys()
            if 'theta_gamma_gamma' in available_branches:
                theta_cut_branch = 'theta_gamma_gamma'
            elif 'theta_pi0_pi0' in available_branches:
                theta_cut_branch = 'theta_pi0_pi0'
            else:
                raise ValueError(f"No valid theta cut branch found in {file_path}")

            # Load all required data
            data = tree.arrays([mx_branch, theta_branch, 'Mx2', 'pTmiss', theta_cut_branch, 'xF', 'Emiss2'], 
                              library="np")
            
            # Apply kinematic cuts
            mask = (
                (data['pTmiss'] < 0.1) &
                (data[theta_cut_branch] < 0.5) &
                (data['xF'] >= -0.15) &
                (data['xF'] <= 0.15) &
                (data['Mx2'] >= -0.02) &  
                (data['Mx2'] <= 0.02) & 
                (data['Emiss2'] >= -0.4) & 
                (data['Emiss2'] <= 0.4)
            )
            
            return data[mx_branch][mask], data[theta_branch][mask]
            
        except Exception as e:
            print(f"Error loading {file_path}: {str(e)}")
            sys.exit(1)

    mx1, theta1 = load_data(args.data1)
    mx2, theta2 = load_data(args.data2)

    output_dir = "output/energy_loss_validation"
    os.makedirs(output_dir, exist_ok=True)

    # Angular bin configuration
    if args.type == 1:
        # Original binning for type 1
        angular_bins = [(0, 25)]
        middle_bins = [(25 + 4.5*i, 25 + 4.5*(i+1)) for i in range(7)]
        angular_bins += middle_bins
        angular_bins += [(56.5, 60.5), (60.5, 90)]
    else:
        # New binning for type 2: 10 evenly spaced bins from 0-30 degrees
        angular_bins = [(3*i, 3*(i+1)) for i in range(10)]

    canvas = ROOT.TCanvas("canvas", "Comparison", 2400, 1800)
    canvas.Divide(4, 3)
    canvas.SetMargin(0.08, 0.02, 0.08, 0.1)

    all_objects = {
        'hists': [],
        'lines': [],
        'legends': [],
        'funcs': []
    }

    # Integrated plot (pad 1)
    canvas.cd(1)
    h1 = ROOT.TH1D("h1_int", "Integrated", 50, xmin, xmax)
    h2 = ROOT.TH1D("h2_int", "Integrated", 50, xmin, xmax)
    all_objects['hists'].extend([h1, h2])
    
    for x in mx1: h1.Fill(x)
    for x in mx2: h2.Fill(x)
    
    # Define Gaussian+linear fit function
    fit_func1 = ROOT.TF1("gaus_poly1", "gaus(0)+pol1(3)", xmin, xmax)
    fit_func1.SetParameters(h1.GetMaximum(), h1.GetMean(), h1.GetRMS(), 0, 0)
    fit_result1 = h1.Fit(fit_func1, "SQ")
    
    fit_func2 = ROOT.TF1("gaus_poly2", "gaus(0)+pol1(3)", xmin, xmax)
    fit_func2.SetParameters(h2.GetMaximum(), h2.GetMean(), h2.GetRMS(), 0, 0)
    fit_result2 = h2.Fit(fit_func2, "SQ")
    
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
        f1 = h1.GetFunction("gaus_poly1")
        if f1:
            f1.SetLineColor(ROOT.kBlack)
            f1.SetLineStyle(2)
            f1.SetLineWidth(2)
            f1.SetRange(xmin, xmax)
            f1.Draw("SAME")
            all_objects['funcs'].append(f1)

    if fit_result2 and fit_result2.IsValid():
        f2 = h2.GetFunction("gaus_poly2")
        if f2:
            f2.SetLineColor(ROOT.kRed)
            f2.SetLineStyle(2)
            f2.SetLineWidth(2)
            f2.SetRange(xmin, xmax)
            f2.Draw("SAME")
            all_objects['funcs'].append(f2)
    
    ROOT.gPad.Modified()
    ROOT.gPad.Update()
    ymax = ROOT.gPad.GetUymax()
    line = ROOT.TLine(vline, 0, vline, ymax)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray+2)
    line.Draw()
    all_objects['lines'].append(line)
    
    leg = ROOT.TLegend(0.2, 0.15, 0.9, 0.28)
    leg.SetBorderSize(1)
    leg.SetFillColor(ROOT.kWhite)
    leg.SetTextSize(0.035)
    if fit_result1 and fit_result1.IsValid():
        leg.AddEntry(h1, f"{args.label1}: #mu={fit_result1.Parameter(1):.3f}#pm{fit_result1.Parameter(2):.3f}", "p")
    else:
        leg.AddEntry(h1, f"{args.label1}: No fit", "p")
    if fit_result2 and fit_result2.IsValid():
        leg.AddEntry(h2, f"{args.label2}: #mu={fit_result2.Parameter(1):.3f}#pm{fit_result2.Parameter(2):.3f}", "p")
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
        h1 = ROOT.TH1D(f"h1_{i}", title, 50, xmin, xmax)
        h2 = ROOT.TH1D(f"h2_{i}", "", 50, xmin, xmax)
        all_objects['hists'].extend([h1, h2])
        
        mask1 = (theta1 >= low_rad) & (theta1 < high_rad)
        mask2 = (theta2 >= low_rad) & (theta2 < high_rad)
        for x in mx1[mask1]: h1.Fill(x)
        for x in mx2[mask2]: h2.Fill(x)
        
        fit1_valid = False
        fit2_valid = False
        f1, f2 = None, None
        
        if h1.GetEntries() > 10:
            fit_func_name = f"gaus_poly1_bin{i}"
            fit_func = ROOT.TF1(fit_func_name, "gaus(0)+pol1(3)", xmin, xmax)
            fit_func.SetParameters(h1.GetMaximum(), h1.GetMean(), h1.GetRMS(), 0, 0)
            fit_result1 = h1.Fit(fit_func, "SQ")
            if fit_result1 and fit_result1.IsValid():
                fit1_valid = True
                f1 = h1.GetFunction(fit_func_name)
        
        if h2.GetEntries() > 10:
            fit_func_name = f"gaus_poly2_bin{i}"
            fit_func = ROOT.TF1(fit_func_name, "gaus(0)+pol1(3)", xmin, xmax)
            fit_func.SetParameters(h2.GetMaximum(), h2.GetMean(), h2.GetRMS(), 0, 0)
            fit_result2 = h2.Fit(fit_func, "SQ")
            if fit_result2 and fit_result2.IsValid():
                fit2_valid = True
                f2 = h2.GetFunction(fit_func_name)
        
        h1.Draw("PE")
        h2.Draw("PE SAME")

        if f1:
            f1.SetLineColor(ROOT.kBlack)
            f1.SetLineStyle(2)
            f1.SetLineWidth(2)
            f1.Draw("SAME")
            all_objects['funcs'].append(f1)
            
        if f2:
            f2.SetLineColor(ROOT.kRed)
            f2.SetLineStyle(2)
            f2.SetLineWidth(2)
            f2.Draw("SAME")
            all_objects['funcs'].append(f2)
        
        for h, color in [(h1, ROOT.kBlack), (h2, ROOT.kRed)]:
            h.SetLineColor(color)
            h.SetMarkerColor(color)
            h.SetMarkerStyle(20)
            h.GetXaxis().SetTitle(xlabel)
            h.GetYaxis().SetTitle("Counts")
            h.GetXaxis().SetTitleSize(0.06)
            h.GetYaxis().SetTitleSize(0.06)
        
        ROOT.gPad.Modified()
        ROOT.gPad.Update()
        ymax = ROOT.gPad.GetUymax()
        line = ROOT.TLine(vline, 0, vline, ymax)
        line.SetLineStyle(2)
        line.SetLineColor(ROOT.kGray+2)
        line.Draw()
        all_objects['lines'].append(line)
        
        leg = ROOT.TLegend(0.2, 0.15, 0.9, 0.28)
        leg.SetBorderSize(1)
        leg.SetFillColor(ROOT.kWhite)
        leg.SetTextSize(0.035)
        if fit1_valid:
            leg.AddEntry(h1, f"{args.label1}: #mu={fit_result1.Parameter(1):.3f}#pm{fit_result1.Parameter(2):.3f}", "p")
        else:
            leg.AddEntry(h1, f"{args.label1}: No fit", "p")
        if fit2_valid:
            leg.AddEntry(h2, f"{args.label2}: #mu={fit_result2.Parameter(1):.3f}#pm{fit_result2.Parameter(2):.3f}", "p")
        else:
            leg.AddEntry(h2, f"{args.label2}: No fit", "p")
        leg.Draw()
        all_objects['legends'].append(leg)
        
        fit_results1.append((fit_result1.Parameter(1), fit_result1.Parameter(2)) if fit1_valid else (0,0))
        fit_results2.append((fit_result2.Parameter(1), fit_result2.Parameter(2)) if fit2_valid else (0,0))
        bin_centers.append((low + high)/2)

    # Final plot (pad 12)
    canvas.cd(12)
    gr1 = ROOT.TGraphErrors(len(bin_centers))
    gr2 = ROOT.TGraphErrors(len(bin_centers))
    
    for i in range(len(bin_centers)):
        gr1.SetPoint(i, bin_centers[i] - 0.25, fit_results1[i][0])
        gr1.SetPointError(i, 0, fit_results1[i][1])
        gr2.SetPoint(i, bin_centers[i] + 0.25, fit_results2[i][0])
        gr2.SetPointError(i, 0, fit_results2[i][1])
    
    mg = ROOT.TMultiGraph()
    mg.Add(gr1)
    mg.Add(gr2)
    mg.SetTitle(";#theta (degrees);Mean " + xlabel)
    
    gr1.SetMarkerColor(ROOT.kBlack)
    gr1.SetLineColor(ROOT.kBlack)
    gr1.SetMarkerStyle(20)
    gr1.SetMarkerSize(1.2)
    gr1.SetLineWidth(2)

    gr2.SetMarkerColor(ROOT.kRed)
    gr2.SetLineColor(ROOT.kRed)
    gr2.SetMarkerStyle(21)
    gr2.SetMarkerSize(1.2)
    gr2.SetLineWidth(2)

    mg.Draw("AP")
    
    # Set x-axis limits based on type
    if args.type == 1:
        mg.GetXaxis().SetLimits(20, 65)
    else:
        mg.GetXaxis().SetLimits(0, 30)
    
    if args.type == 1:
        mg.SetMinimum(-0.2)
        mg.SetMaximum(0.2)
    else:
        mg.SetMinimum(0.5)
        mg.SetMaximum(1.2)
    
    # Adjust horizontal line based on type
    if args.type == 1:
        hline = ROOT.TLine(20, hline_y, 65, hline_y)
    else:
        hline = ROOT.TLine(0, hline_y, 30, hline_y)
        
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