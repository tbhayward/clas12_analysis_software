// Created 9/6/23
// Created to compare epi+X and epX nSidis distributions between pass-1 and preliminary pass-2 

#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TPaveText.h>
#include <TAxis.h>
#include <iostream>
#include <string>
#include <map>
#include <cstring>
#include <algorithm>
#include <cmath>

using namespace std;

std::pair<std::vector<double>, std::vector<double>> calculateAndPlotALU(
    TTree* tree1, const char* branchName, double min_val, double max_val) {

    cout << "Entering calculateAndPlotALU ";
    std::vector<double> ALU_values;
    std::vector<double> ALU_errors;

    // Create a 2D array to hold N+ and N- for each dynamic bin and phi bin
    std::vector<std::vector<double>> N_pos(6, std::vector<double>(12, 0));  
    // 6 dynamic bins, 12 phi bins
    std::vector<std::vector<double>> N_neg(6, std::vector<double>(12, 0));

    std::vector<double> sum_beam_pol(6, 0.0);
    std::vector<int> count_beam_pol(6, 0);
    std::vector<double> sum_W_over_A(6, 0.0);
    std::vector<int> count_W_over_A(6, 0);

    // Loop through the tree to fill N_pos and N_neg
    int helicity;
    double branch_var, phi, beam_pol, W, A;
    tree1->SetBranchAddress(branchName, &branch_var);
    tree1->SetBranchAddress("phi", &phi);
    tree1->SetBranchAddress("helicity", &helicity);
    tree1->SetBranchAddress("beam_pol", &beam_pol);
    tree1->SetBranchAddress("DepW", &W);
    tree1->SetBranchAddress("DepA", &A);

    for (int entry = 0; entry < tree1->GetEntries(); ++entry) {
        tree1->GetEntry(entry);

        // cout << branchName << " branch_var = " << branch_var << " min_val = " << min_val;
        // cout << " max_val = " << max_val << " phi = " << phi << " " << endl;
        if(branch_var < min_val || branch_var > max_val) continue;  
        // Skip entries out of range
        if(phi < 0 || phi > 2 * TMath::Pi()) continue;  
        // Skip entries out of range

        int dyn_bin = int((branch_var - min_val) / ((max_val - min_val) / 6));
        int phi_bin = int(phi / (2 * TMath::Pi() / 12));

        if(dyn_bin < 0 || dyn_bin >= 6) continue;  // Skip invalid indices
        if(phi_bin < 0 || phi_bin >= 12) continue;  // Skip invalid indices

        if (helicity > 0) {
            N_pos[dyn_bin][phi_bin]++;
        } else if (helicity < 0) {
            N_neg[dyn_bin][phi_bin]++;
        }
        if (helicity != 0) {  // assuming non-zero helicity implies valid polarization
            sum_beam_pol[dyn_bin] += beam_pol;
            count_beam_pol[dyn_bin]++;

            double W_over_A = (A != 0) ? W / A : 0;
            sum_W_over_A[dyn_bin] += W_over_A;
            count_W_over_A[dyn_bin]++;
        }

    }

    for (int dyn_bin = 0; dyn_bin < 6; ++dyn_bin) {
        TF1 fitFunc("fitFunc", "[0]*sin(x)", 0, 2 * TMath::Pi());
        TGraphErrors fitGraph;
        for (int phi_bin = 0; phi_bin < 12; ++phi_bin) {
            double phi_val = phi_bin * (2 * TMath::Pi() / 12);
            double mean_beam_pol = (count_beam_pol[dyn_bin] != 0) ? sum_beam_pol[dyn_bin] 
                / count_beam_pol[dyn_bin] : 1.0;
            double ALU = (1 / mean_beam_pol) * 
                (N_pos[dyn_bin][phi_bin] - N_neg[dyn_bin][phi_bin]) / 
                (N_pos[dyn_bin][phi_bin] + N_neg[dyn_bin][phi_bin]);
            double ALU_error = (2 / mean_beam_pol) * 
                TMath::Sqrt((N_pos[dyn_bin][phi_bin] * N_neg[dyn_bin][phi_bin]) / 
                TMath::Power(N_pos[dyn_bin][phi_bin] + N_neg[dyn_bin][phi_bin], 3));

            fitGraph.SetPoint(phi_bin, phi_val, ALU);
            fitGraph.SetPointError(phi_bin, 0, ALU_error);
        }
        fitGraph.Fit(&fitFunc, "Q");
        double A = fitFunc.GetParameter(0);
        double A_error = fitFunc.GetParError(0);
        double mean_W_over_A = (count_W_over_A[dyn_bin] != 0) ? 
            sum_W_over_A[dyn_bin] / count_W_over_A[dyn_bin] : 1.0;
        ALU_values.push_back(A / mean_W_over_A);
        ALU_errors.push_back(A_error / mean_W_over_A);
        cout << "Made it passed ALUs" << endl;
    }
    
    return std::make_pair(ALU_values, ALU_errors);
}


struct HistConfig {
    int bins;
    double min;
    double max;
};

std::map<std::string, HistConfig> histConfigs = {
    {"beam_pol", {20, 0.80, 1.00}},
    {"DepA", {200, 0, 1}},
    {"DepB", {200, 0, 1}},
    {"DepC", {200, 0, 1}},
    {"DepV", {200, 0, 2}},
    {"DepW", {200, 0, 1}},
    {"e_p", {200, 2, 8}},
    {"e_phi", {200, 0, 360}},
    {"eta", {200, -1, 3}},
    {"e_theta", {200, 0, 40}}, // Convert degree to radian
    {"evnum", {200, 0, 0}},
    {"helicity", {2, -2, 2}},
    {"Mx", {200, -4, 3}},
    {"Mx2", {200, -10, 10}},
    {"phi", {200, 0, 2 * M_PI}},
    {"p_p", {200, 0, 6}},
    {"p_phi", {200, 0, 360}},
    {"pT", {200, 0, 1.2}},
    {"p_theta", {200, 0, 60}}, // Convert degree to radian
    {"Q2", {200, 0, 9}},
    {"runnum", {200, 0, 0}},
    {"t", {200, -10, 1}},
    {"tmin", {200, -0.5, 0}},
    {"vz_e", {200, -15, 15}},
    {"vz_p", {200, -15, 15}},
    {"W", {200, 2, 4}},
    {"x", {200, 0, 0.6}},
    {"xF", {200, -1, 1}},
    {"y", {200, 0, 1}},
    {"z", {200, 0, 1}},
    {"zeta", {200, 0.3, 1}}
};

std::string formatBranchName(const std::string& original) {
    std::map<std::string, std::string> specialLabels = {
        {"Q2", "Q^{2} (GeV^{2})"},
        {"W", "W (GeV)"},
        {"pT", "P_{T} (GeV)"},
        {"t", "t (GeV^{2})"},
        {"tmin", "t_{min} (GeV^{2})"},
        {"e_p", "e_{p} (GeV)"},
        {"Mx", "M_{x} (GeV)"},
        {"Mx2", "M_{x}^{2} (GeV)"},
        {"p_p", "p_{p} (GeV)"},
        {"xF", "x_{F}"},
    };
  
    if (specialLabels.find(original) != specialLabels.end()) {
        return specialLabels[original];
    }

    std::string formatted = original;
    size_t pos = 0;
    while ((pos = formatted.find('_', pos)) != std::string::npos) {
        formatted.replace(pos, 1, "_{");
        size_t closing = formatted.find('_', pos + 2);
        if (closing == std::string::npos) {
            closing = formatted.length();
        }
        formatted.insert(closing, "}");
        pos = closing + 1;
    }

    if (formatted.find("theta") != std::string::npos) {
        formatted.replace(formatted.find("theta"), 5, "#theta");
    }

    if (formatted.find("zeta") != std::string::npos) {
        formatted.replace(formatted.find("zeta"), 5, "#zeta");
    }

    if (formatted.find("phi") != std::string::npos) {
        formatted.replace(formatted.find("phi"), 3, "#phi");
    }

    if (formatted.find("eta") != std::string::npos && 
        formatted.find("theta") == std::string::npos && 
        formatted.find("zeta") == std::string::npos) {
        formatted.replace(formatted.find("eta"), 3, "#eta");
    }
  
    return formatted;
}

void createHistograms(TTree* tree1, TTree* tree2, 
    std::string data_set_1_name, std::string data_set_2_name, const char* outDir) {
    TObjArray* branches1 = tree1->GetListOfBranches();
    TObjArray* branches2 = tree2->GetListOfBranches();

    if (branches1->GetEntries() != branches2->GetEntries()) {
        std::cout << "Number of branches mismatch. Exiting." << std::endl;
        return;
    }

    // for (int i = 0; i < branches1->GetEntries(); ++i) {
    for (int i = 13; i < 17; ++i) {
        cout << "we're starting stuff!" << endl;
        const char* branchName = branches1->At(i)->GetName();
        if (std::strcmp(branchName, "runnum") == 0 || std::strcmp(branchName, "evnum") == 0 || 
            std::strcmp(branchName, "phi") == 0 || std::strcmp(branchName, "beam_pol") == 0 || 
            std::strcmp(branchName, "helicity") == 0
            // || std::strcmp(branchName, "Mx") == 0
            // || std::strcmp(branchName, "Mx2") == 0
            ) {
            continue;
        }

        HistConfig config = histConfigs[branchName];
        TH1F hist1(Form("%s_1", branchName), "", config.bins, config.min, config.max);
        TH1F hist2(Form("%s_2", branchName), "", config.bins, config.min, config.max);

        std::string cutCondition = "";
        if (std::strcmp(branchName, "Mx") != 0 && std::strcmp(branchName, "Mx2") != 0) {
            cutCondition = "Mx > 1.5";
        } 

        cout << "passed the cut condition" << endl;
        // Draw a temporary histogram to get statistics
        TH1F tempHist(Form("temp_%s", branchName), "", 1000, config.min, config.max);  
        // 1000 bins for better statistics
        tree1->Draw(Form("%s>>temp_%s", branchName, branchName), cutCondition.c_str());

        // Find the quantile edges
        int nQuantiles = 6;
        double quantiles[nQuantiles];
        double sum = tempHist.GetEntries();
        for (int i = 1; i <= nQuantiles; ++i) {
            quantiles[i-1] = i * (sum / nQuantiles);
        }
        double edges[nQuantiles + 1];
        tempHist.GetQuantiles(nQuantiles, edges, quantiles);
        cout << "passed the quantiles" << endl;
        std::string formattedBranchName = formatBranchName(branchName);
        TCanvas canvas(branchName, "Canvas", 1600, 600);  // Width doubled for side-by-side panels
        TPad *pad1 = new TPad("pad1", "The pad with the function",0.0,0.0,0.33,1.0,21);
        pad1->SetLeftMargin(0.2); pad1->SetBottomMargin(0.2);
        pad1->SetFillColor(0);  // Set the fill color to white for pad1
        pad1->Draw();
        pad1->cd();  // Set current pad to pad1
        cout << "created pad1" << endl;
        std::string drawCommand1 = Form("%s>>%s_1", branchName, branchName);
        std::string drawCommand2 = Form("%s>>%s_2", branchName, branchName);

        // Convert to degrees if necessary
        if (strcmp(branchName, "e_phi") == 0 || strcmp(branchName, "e_theta") == 0 ||
            strcmp(branchName, "p_phi") == 0 || strcmp(branchName, "p_theta") == 0) {
            drawCommand1 = Form("%s * (180 / TMath::Pi())>>%s_1", branchName, branchName);
            drawCommand2 = Form("%s * (180 / TMath::Pi())>>%s_2", branchName, branchName);
        }

        cout << "before, cutCondition = " << cutCondition.c_str() << endl;
        tree1->Draw(drawCommand1.c_str(), cutCondition.c_str());
        tree2->Draw(drawCommand2.c_str(), cutCondition.c_str());
        cout << "after cutCondition" << endl;

        hist1.SetLineColor(kRed);
        hist2.SetLineColor(kBlue);
        hist1.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
        hist1.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
        hist1.GetXaxis()->SetTitleSize(0.05);  // Increase x-axis title size
        hist1.GetYaxis()->SetTitleSize(0.05);  // Increase y-axis title size
        hist2.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
        hist2.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
        hist2.GetXaxis()->SetTitleSize(0.05);  // Increase x-axis title size
        hist2.GetYaxis()->SetTitleSize(0.05);  // Increase y-axis title size

        hist1.Draw(""); 
        hist2.Draw("same"); 

        hist1.SetStats(0);
        hist2.SetStats(0);

        hist1.GetXaxis()->SetTitle(formattedBranchName.c_str());
        hist1.GetYaxis()->SetTitle("Counts");

        double max_value = std::max(hist1.GetMaximum(), hist2.GetMaximum());
        hist1.SetMaximum(max_value * 1.2);
        hist2.SetMaximum(max_value * 1.2);


        // Create the legend
        TLegend *leg1 = new TLegend(0.25, 0.8, 0.9, 0.9);  // Adjust these coordinates as needed
        leg1->SetBorderSize(1);  // border size
        leg1->SetFillColor(0);  // Transparent fill
        leg1->SetTextSize(0.04);  // text size
        // Add entries for each histogram
        leg1->AddEntry(&hist1, Form("%s counts: %d", data_set_1_name.c_str(), 
            int(hist1.GetEntries())), "l");
        leg1->AddEntry(&hist2, Form("%s counts: %d", data_set_2_name.c_str(), 
            int(hist2.GetEntries())), "l");
        // Draw the legend
        leg1->Draw("same");


        // pad with ratio
        canvas.cd();  // Switch back to the main canvas before creating a new pad
        TPad *pad2 = new TPad("pad2", "The pad with the ratio",0.33,0.0,0.66,1.0,21);
        pad2->SetLeftMargin(0.2); pad2->SetBottomMargin(0.2);
        pad2->SetFillColor(0);  // Set the fill color to white for pad2
        pad2->Draw();
        pad2->cd();  // Set current pad to pad2
        TH1F ratioHist(Form("%s_ratio", branchName), "", config.bins, config.min, config.max);
        ratioHist.Divide(&hist2, &hist1);
        ratioHist.SetLineColor(kBlack);
        ratioHist.SetMinimum(0.75);  // Set Y-range
        ratioHist.SetMaximum(2.50);  // Set Y-range
        ratioHist.GetXaxis()->SetTitle(formattedBranchName.c_str());
        std::string yAxisTitle = "ratio";
        ratioHist.GetYaxis()->SetTitle(yAxisTitle.c_str());
        ratioHist.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
        ratioHist.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
        ratioHist.GetXaxis()->SetTitleSize(0.05);  // Increase x-axis title size
        ratioHist.GetYaxis()->SetTitleSize(0.05);  // Increase y-axis title size
        ratioHist.Draw("HIST");
        ratioHist.SetStats(0);  // Disable the statistical box

        // Create the legend
        TLegend *leg2 = new TLegend(0.45, 0.8, 0.9, 0.9);  // Adjust these values as needed
        leg2->SetBorderSize(1);  // border size
        leg2->SetFillColor(0);  // Transparent fill
        leg2->SetTextSize(0.04);  // text size

        // Add text as legend entries. No associated object, so the last parameter is "".
        double overallRatio = (hist2.GetEntries() != 0 && hist1.GetEntries() != 0) ? 
                    hist2.GetEntries() / hist1.GetEntries() : 0;
        leg2->AddEntry((TObject*)0, Form("Overall Ratio: %.2f", overallRatio), "");

        // Draw the legend
        leg2->Draw("same");


        // Third Panel for ALU calculations and fitting
        canvas.cd();  // Switch back to the main canvas before creating a new pad
        TPad *pad3 = new TPad("pad3", "The pad with ALU", 0.7, 0.0, 1.0, 1.0, 21);
        pad3->SetLeftMargin(0.2); pad3->SetBottomMargin(0.2);
        pad3->SetFillColor(0);  // Set the fill color to white for pad3
        pad3->Draw();
        pad3->cd();  // Set current pad to pad3

        // Get min and max values for the branch
        double min_val = hist1.GetXaxis()->GetXmin();
        double max_val = hist1.GetXaxis()->GetXmax();

        std::pair<std::vector<double>, std::vector<double>> result1, result2;

        result1 = calculateAndPlotALU(tree1, branchName, min_val, max_val);
        result2 = calculateAndPlotALU(tree2, branchName, min_val, max_val);

        TGraphErrors aluGraph1(6), aluGraph2(6);  // We have 6 dynamic bins for each

        double offset = 0.1 * ((max_val - min_val) / 6);  // 10% of bin width
        // Populate aluGraph1 and aluGraph2 using result1 and result2
        // (This part can be put into a loop or function for efficiency)
        for (int dyn_bin = 0; dyn_bin < 6; ++dyn_bin) {
            double bin_center = min_val + (dyn_bin + 0.5) * (max_val - min_val) / 6;
            aluGraph1.SetPoint(dyn_bin, bin_center, result1.first[dyn_bin]);
            aluGraph1.SetPointError(dyn_bin, 0, result1.second[dyn_bin]);
            
            aluGraph2.SetPoint(dyn_bin, bin_center + offset, result2.first[dyn_bin]); // Add offset
            aluGraph2.SetPointError(dyn_bin, 0, result2.second[dyn_bin]);
        }
        

        aluGraph1.SetLineColor(kRed); aluGraph1.SetMarkerColor(kRed);
        aluGraph1.SetMarkerStyle(20);
        aluGraph1.SetMarkerSize(1.1);

        aluGraph2.SetLineColor(kBlue); aluGraph2.SetMarkerColor(kBlue);
        aluGraph2.SetMarkerStyle(21);
        aluGraph2.SetMarkerSize(1.1);

        aluGraph1.Draw("AP");
        aluGraph1.GetYaxis()->SetRangeUser(-0.04, 0.04);
        aluGraph1.GetYaxis()->SetTitle("F_{LU}^{sin#phi} / F_{UU}");
        aluGraph1.GetXaxis()->SetTitle(formattedBranchName.c_str());
        aluGraph1.SetTitle("");  // Remove title
        aluGraph1.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
        aluGraph1.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
        aluGraph1.GetXaxis()->SetTitleSize(0.05);  // Increase x-axis title size
        aluGraph1.GetYaxis()->SetTitleSize(0.05);  // Increase y-axis title size

        aluGraph2.Draw("Psame");
        aluGraph2.SetTitle("");  // Remove title
        aluGraph2.GetXaxis()->SetLabelSize(0.04);  // Increase x-axis label size
        aluGraph2.GetYaxis()->SetLabelSize(0.04);  // Increase y-axis label size
        aluGraph2.GetXaxis()->SetTitleSize(0.05);  // Increase x-axis title size
        aluGraph2.GetYaxis()->SetTitleSize(0.05);  // Increase y-axis title size

        // Create the legend at x1, y1, x2, y2
        TLegend *leg3 = new TLegend(0.63, 0.8, 0.9, 0.9);  // Adjust these values as needed
        // Add entries
        leg3->AddEntry(&aluGraph1, data_set_1_name.c_str(), "p");
        leg3->AddEntry(&aluGraph2, data_set_2_name.c_str(), "p");
        // Set the marker colors to match your graphs
        aluGraph1.SetMarkerColor(kRed);
        aluGraph2.SetMarkerColor(kBlue);
        // Customize the legend
        leg3->SetBorderSize(1);  // border size
        leg3->SetTextSize(0.04);  // text size
        // Draw the legend
        leg3->Draw("same");


        // Save the canvas
        canvas.SaveAs(Form("%s/%s.png", outDir, branchName));

        cout << "We're about to delete stuff!" << endl;

        // Delete or remove from directory all dynamically created objects
        hist1.SetDirectory(0);
        hist2.SetDirectory(0);
        ratioHist.SetDirectory(0);

        delete pad1;
        delete pad2;
        delete pad3;
        delete leg1;
        delete leg2;
        delete leg3;
    }
}

void compare_files(std::string root_file1_path, std::string root_file2_path, 
    std::string data_set_1_name, std::string data_set_2_name) {
    gStyle->SetCanvasColor(0);

    TFile* file1 = new TFile(root_file1_path.c_str(), "READ");
    TFile* file2 = new TFile(root_file2_path.c_str(), "READ");

    if (!file1->IsOpen() || !file2->IsOpen()) {
        cout << "Error opening ROOT files (is the location correct?). Exiting." << endl;
    }

    TTree* tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree* tree2 = (TTree*)file2->Get("PhysicsEvents");

    if (!tree1 || !tree2) {
        cout << "Error getting trees from ROOT files." << endl;
    }

    createHistograms(tree1, tree2, data_set_1_name, data_set_2_name, "output");

    file1->Close(); delete file1;
    file2->Close(); delete file2;
}
