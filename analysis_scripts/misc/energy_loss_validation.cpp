#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TMath.h>
#include <vector>
#include <string>
#include <iostream>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
#include <TF1.h>
#include <TGraphErrors.h>

// Check if the branch exists in the tree
bool BranchExists(TTree* tree, const char* branchName) {
    TBranch* branch = tree->GetBranch(branchName);
    return (branch != nullptr);
}


void compareTrees(const char* file1, const char* file2, const char* output,
    const char* output2, double lineValue) {
    // Define the momentum bin edges
    
    // // electrons
    // std::vector<double> binEdges = {6.0,7,7.25,7.5,7.75,8,8.25,8.5,8.67,9,9.15,9.3,9.5,10.6};

    // // pions and kaon
    // std::vector<double> binEdges = {0.6,1.3,1.8,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.8,3.0,3.2,3.4,3.8,5.0};

    // pions and kaon less
    std::vector<double> binEdges = {0.6,2.6,3.8,5.0};


    // // proton 
    // std::vector<double> binEdges = 
    //     {0,0.6,0.625,0.65,0.675,0.7,0.725,0.75,0.8,1.0,1.1,1.2,1.4,1.8,2.4,3.0};
    

    int nBins = binEdges.size() - 1;

    // Open ROOT files and get trees
    TFile* f1 = new TFile(file1);
    TFile* f2 = new TFile(file2);
    TTree* tree1 = (TTree*)f1->Get("PhysicsEvents"); // Replace "tree_name" with your tree name
    TTree* tree2 = (TTree*)f2->Get("PhysicsEvents"); // Same as above

    // Create histograms for each bin
    std::vector<TH1D*> hist1, hist2;
    for (int i = 0; i < nBins; ++i) {
        // electron proton 
        hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 0.75, 1.1));
        hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 0.75, 1.1));

        // // pi+ 
        // hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 0.6, 1.2));
        // hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 0.6, 1.2));

        // // pi-
        // hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 1.0, 1.8));
        // hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 1.0, 1.8));

        // // k+ 
        // hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 1.0, 1.4));
        // hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 1.0, 1.4));

        // // k- 
        // hist1.push_back(new TH1D(Form("hist1_%d", i), "", 50, 0.5, 2.5));
        // hist2.push_back(new TH1D(Form("hist2_%d", i), "", 50, 0.5, 2.5));

        // // proton rho
        // hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 0.1, 0.9));
        // hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 0.1, 0.9));

        // // // proton f2
        // hist1.push_back(new TH1D(Form("hist1_%d", i), "", 100, 1.2, 1.8));
        // hist2.push_back(new TH1D(Form("hist2_%d", i), "", 100, 1.2, 1.8));
    }

    // Set branch addresses
    double p_p, Mx2;
    // Check if 'p_p' branch exists in tree1
    if (BranchExists(tree1, "p2_p")) {
        tree1->SetBranchAddress("p2_p", &p_p);
        tree1->SetBranchAddress("Mx2", &Mx2);
    } else if (BranchExists(tree1, "p_p")){
        tree1->SetBranchAddress("p_p", &p_p);
        tree1->SetBranchAddress("Mx2", &Mx2);
    } else {
        // Use 'e_p' and 'W' instead
        tree1->SetBranchAddress("e_p", &p_p);
        tree1->SetBranchAddress("W", &Mx2);
    }

     // Check if 'p_p' branch exists in tree2
    if (BranchExists(tree2, "p2_p")) {
        tree2->SetBranchAddress("p2_p", &p_p);
        tree2->SetBranchAddress("Mx2", &Mx2);
    } else if (BranchExists(tree2, "p_p")){
        tree2->SetBranchAddress("p_p", &p_p);
        tree2->SetBranchAddress("Mx2", &Mx2);
    } else {
        // Use 'e_p' and 'W' instead
        tree2->SetBranchAddress("e_p", &p_p);
        tree2->SetBranchAddress("W", &Mx2);
    }

    // Fill histograms
    Long64_t nEntries1 = tree1->GetEntries();
    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        for (int b = 0; b < nBins; ++b) {
            if (p_p >= binEdges[b] && p_p < binEdges[b+1])
                hist1[b]->Fill(Mx2);
        }
    }
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        for (int b = 0; b < nBins; ++b) {
            if (p_p >= binEdges[b] && p_p < binEdges[b+1])
                hist2[b]->Fill(Mx2);
        }
    }

    // Before your loop
    std::vector<double> binCenters, meanValues1, meanErrors1, meanValues2, meanErrors2;

    TCanvas* c1 = new TCanvas("c1", "Comparison", 1200, 800);
    c1->Divide(TMath::CeilNint(sqrt(nBins)), TMath::CeilNint(sqrt(nBins)));
    double globalFontSize = 0.04; // You can adjust this value as needed
    for (int i = 0; i < nBins; ++i) {
        c1->cd(i+1);
        gPad->SetLeftMargin(0.15);  

        // Setting the title of the histogram to include the bin range
        char title[100];
        sprintf(title, "%.3f < p (GeV) < %.3f; M_{x}^{2} (GeV^{2}); Counts", 
            binEdges[i], binEdges[i+1]);
        hist1[i]->SetTitle(title);

        // Find the maximum value in both histograms for this bin
        double maxVal = TMath::Max(hist1[i]->GetMaximum(), hist2[i]->GetMaximum());

        // Set the range of y-axis to 0 - 10% more than the max value
        hist1[i]->SetMaximum(maxVal * 1.4); 

        // Set x and y axis labels
        hist1[i]->GetXaxis()->SetTitle("M_{x}^{2} (GeV^{2})");
        hist1[i]->GetYaxis()->SetTitle("Counts");

        // Draw the histograms
        hist1[i]->SetLineColor(kBlue);
        hist2[i]->SetLineColor(kRed);
        hist1[i]->Draw();
        hist2[i]->Draw("same");

        // Create and draw a vertical line
        TLine *line = new TLine(lineValue, hist1[i]->GetMinimum(), lineValue, maxVal * 1.10);
        line->SetLineStyle(2); // Set the line style to dashed
        line->Draw();

        // Remove the statistics box
        gStyle->SetOptStat(0);

        // Adjust font size and center axis labels
        hist1[i]->GetXaxis()->SetTitleSize(globalFontSize);
        hist1[i]->GetYaxis()->SetTitleSize(globalFontSize);
        hist1[i]->GetXaxis()->SetLabelSize(globalFontSize);
        hist1[i]->GetYaxis()->SetLabelSize(globalFontSize);
        
        hist1[i]->GetXaxis()->CenterTitle();
        hist1[i]->GetYaxis()->CenterTitle();

        // Define the fit function range
        double xMin = hist1[i]->GetXaxis()->GetXmin();
        double xMax = hist1[i]->GetXaxis()->GetXmax();
        // Gaussian + Linear Background Function
        TF1 *fitFunc1 = new TF1(Form("fitFunc1_%d", i), "gaus(0) + pol2(3)", xMin, xMax);
        TF1 *fitFunc2 = new TF1(Form("fitFunc2_%d", i), "gaus(0) + pol2(3)", xMin, xMax);
        // gaus(0): Gaussian part with parameters [0, 1, 2] (amplitude, mean, sigma)
        // pol1(3): Linear background with parameters [3, 4, 5] (constant, slope, x^2)

        // Set initial parameter estimates for the fit function
        fitFunc1->SetParameters(100, lineValue, 0.1, 0, 1); // Example values, adjust as needed
        fitFunc1->SetParLimits(0, 0, 10e6); // amplitude limits
        fitFunc1->SetParLimits(1, 0, lineValue + 0.1); // Limit the mean around lineValue
        fitFunc1->SetParLimits(2, 0.025, 1); // Limit the mean around lineValue

        fitFunc2->SetParameters(100, lineValue, 0.1, 0, 1); // Example values, adjust as needed
        fitFunc2->SetParLimits(0, 0, 10e6); // amplitude limits
        fitFunc2->SetParLimits(1, 0, lineValue + 0.1); // Limit the mean around lineValue

        // Perform the fit on the first histogram
        hist1[i]->Fit(fitFunc1, "R");
        fitFunc1->SetLineColor(kBlue);
        fitFunc1->SetLineStyle(1);
        fitFunc1->Draw("SAME");

        // Perform the fit on the second histogram
        hist2[i]->Fit(fitFunc2, "R+");
        fitFunc2->SetLineColor(kRed);
        fitFunc2->Draw("SAME");

        // Retrieve the mean and its error for the first fit
        double mean1 = fitFunc1->GetParameter(1);
        double meanError1 = fitFunc1->GetParError(1);

        // Retrieve the mean and its error for the second fit
        double mean2 = fitFunc2->GetParameter(1);
        double meanError2 = fitFunc2->GetParError(1);

        // Create and add a legend with fit results
        TLegend* legend = new TLegend(0.15, 0.7, 0.63, 0.9); // Adjust these coordinates as needed
        legend->SetTextSize(0.04); // Set the text size. Adjust as needed.
        char entry1[100], entry2[100];
        sprintf(entry1, "Uncorrected, #mu = %.3f #pm %.3f", mean1, meanError1);
        sprintf(entry2, "Corrected, #mu = %.3f #pm %.3f", mean2, meanError2);
        legend->AddEntry(hist1[i], entry1, "l");
        legend->AddEntry(hist2[i], entry2, "l");
        legend->Draw();

        // Inside your loop, after fitting
        double binCenter = (binEdges[i] + binEdges[i+1]) / 2.0;
        binCenters.push_back(binCenter);
        meanValues1.push_back(mean1); // mean1 from the fit of hist1
        meanErrors1.push_back(meanError1); // meanError1 from the fit of hist1
        meanValues2.push_back(mean2); // mean2 from the fit of hist2
        meanErrors2.push_back(meanError2); // meanError2 from the fit of hist2

    }

    // After the loop
    c1->cd(nBins + 1); // Create a new pad for the final plot

    int nPoints = binCenters.size();

    // Graph for the first histogram
    TGraphErrors* gr1 = new TGraphErrors(nPoints, &binCenters[0], &meanValues1[0], 0, 
        &meanErrors1[0]);
    gr1->SetMarkerColor(kBlue);
    gr1->SetLineColor(kBlue);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.6);
    gr1->SetTitle("Fitted Mean Values; Momentum (GeV); Mean of Gaussian (GeV^{2})");

    // Graph for the second histogram
    TGraphErrors* gr2 = new TGraphErrors(nPoints, &binCenters[0], &meanValues2[0], 0, 
        &meanErrors2[0]);
    gr2->SetMarkerColor(kRed);
    gr2->SetLineColor(kRed);
    gr2->SetMarkerStyle(21);
    gr2->SetMarkerSize(0.6);

    // Draw the graphs
    gr1->Draw("AP");
    gr2->Draw("P SAME");

    // Set the Y-axis range to be lineValue +/- offset
    // double yAxisMin = lineValue - 0.12; // proton, pion
    // double yAxisMax = lineValue + 0.04; // proton, pion
    double yAxisMin = lineValue - 0.04; // kaon
    double yAxisMax = lineValue + 0.12; // kaon
    gr1->GetYaxis()->SetRangeUser(yAxisMin, yAxisMax);
    gr2->GetYaxis()->SetRangeUser(yAxisMin, yAxisMax);

    // Center the axis labels
    gr1->GetXaxis()->CenterTitle();
    gr1->GetYaxis()->CenterTitle();
    gr2->GetXaxis()->CenterTitle();
    gr2->GetYaxis()->CenterTitle();

    // Draw a horizontal line at lineValue
    TF1 *line = new TF1("line", Form("%f", lineValue), binEdges.front(), binEdges.back());
    line->SetLineColor(kBlack);
    line->SetLineStyle(3); // Dashed line
    line->Draw("SAME");

    // Create and add a legend
    TLegend* legend = new TLegend(0.1, 0.7, 0.5, 0.9);
    legend->SetTextSize(0.04);
    legend->AddEntry(gr1, "Uncorrected", "p");
    legend->AddEntry(gr2, "Corrected", "p");
    legend->Draw();

    // Optionally, add more settings for axis labels, title, etc.
    c1->SaveAs(output);

    TCanvas* c2 = new TCanvas("c1", "Comparison", 1200, 800);
    c2->cd(1);
    gr1->SetMarkerColor(kBlue);
    gr1->SetLineColor(kBlue);
    gr1->SetMarkerStyle(20);
    gr1->SetMarkerSize(0.6);
    gr2->SetMarkerColor(kRed);
    gr2->SetLineColor(kRed);
    gr2->SetMarkerStyle(20);
    gr2->SetMarkerSize(0.6);
    gr1->Draw("AP");
    gr2->Draw("P SAME");
    line->Draw("SAME");
    // Move the legend to the bottom right
    TLegend* legend2 = new TLegend(0.7, 0.1, 0.9, 0.3);
    legend2->SetTextSize(0.04);
    legend2->AddEntry(gr1, "Uncorrected", "p");
    legend2->AddEntry(gr2, "Corrected", "p");
    legend2->Draw();
    c2->SaveAs(output2);


}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cout << "Usage: " << argv[0] << " <file1> <file2> <output> <output> <peak value>" << std::endl;
        return 1;
    }
    compareTrees(argv[1], argv[2], argv[3], argv[4], atof(argv[5]));
    return 0;
}
