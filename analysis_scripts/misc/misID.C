#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLatex.h>

// void misIDPlot() {
//     // Open the ROOT file and get the tree
//     TFile *file = new TFile("/scratch/thayward/ratios/epi-X_inb.root");
//     TTree *tree = (TTree*)file->Get("PhysicsEvents");

//     // Create a histogram for the fraction calculation
//     TH1F *hFraction = new TH1F("hFraction", ";p (GeV);% k^{-} #rightarrow #pi^{-}", 12, 0, 7);

//     // Loop over the tree and fill the histogram
//     double p_p;
//     int matching_p1_pid;
//     tree->SetBranchAddress("p_p", &p_p);
//     tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

//     for (int i = 0; i < tree->GetEntries(); ++i) {
//         tree->GetEntry(i);
//         if (matching_p1_pid == -321) {
//             hFraction->Fill(p_p);
//         }
//     }

//     // Normalize the histogram to get the fraction
//     hFraction->Scale(100.0 / tree->GetEntries());

//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph = new TGraphErrors();
//     for (int i = 1; i <= hFraction->GetNbinsX(); ++i) {
//         if (hFraction->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph->SetPoint(i-1, hFraction->GetBinCenter(i), hFraction->GetBinContent(i));
//         graph->SetPointError(i-1, 0, hFraction->GetBinError(i)); // Set horizontal error to 0
//     }

//     // Style the graph with markers
//     graph->SetTitle(";p (GeV);% k^{-} #rightarrow #pi^{-}"); // Setting title and axis labels
//     graph->SetMarkerStyle(20);  // Style 20 is a filled circle
//     graph->SetMarkerSize(1.2);  // Adjust the size as needed

//     // Draw the plot using TGraphErrors
//     TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
//     c1->SetLeftMargin(0.15);
//     graph->Draw("AP"); // "AP" to draw the graph with markers and lines

//     // Set axis styles
//     graph->GetXaxis()->CenterTitle();
//     graph->GetXaxis()->SetTitleSize(0.05);
//     graph->GetYaxis()->CenterTitle();
//     graph->GetYaxis()->SetTitleSize(0.05);
//     graph->GetXaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetTitleOffset(1.5); // Adjust if necessary

//     // Save the canvas as a PNG file
//     c1->SaveAs("output/epi-X_misid.png");
// }

void misIDPlot() {
    // Open the ROOT file and get the tree
    TFile *file = new TFile("/scratch/thayward/ratios/ek-X_inb.root");
    TTree *tree = (TTree*)file->Get("PhysicsEvents");

    // Create a histogram for the fraction calculation
    TH1F *hFraction = new TH1F("hFraction", ";p (GeV);% #pi^{-} #rightarrow k^{-}", 12, 0, 7);

    // Loop over the tree and fill the histogram
    double p_p;
    int matching_p1_pid;
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (matching_p1_pid == -211) {
            hFraction->Fill(p_p);
        }
    }

    // Normalize the histogram to get the fraction
    hFraction->Scale(100.0 / tree->GetEntries());

    // Create a TGraphErrors from the histogram with only vertical error bars
    TGraphErrors *graph = new TGraphErrors();
    for (int i = 1; i <= hFraction->GetNbinsX(); ++i) {
        if (hFraction->GetBinContent(i) == 0) {
            continue;
        }
        graph->SetPoint(i-1, hFraction->GetBinCenter(i), hFraction->GetBinContent(i));
        graph->SetPointError(i-1, 0, hFraction->GetBinError(i)); // Set horizontal error to 0
    }

    // Style the graph with markers
    graph->SetTitle(";p (GeV);% #pi^{-} #rightarrow k^{-}"); // Setting title and axis labels
    graph->SetMarkerStyle(20);  // Style 20 is a filled circle
    graph->SetMarkerSize(1.2);  // Adjust the size as needed

    // Draw the plot using TGraphErrors
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    c1->SetLeftMargin(0.15);
    graph->Draw("AP"); // "AP" to draw the graph with markers and lines

    // Set axis styles
    graph->GetXaxis()->CenterTitle();
    graph->GetXaxis()->SetTitleSize(0.05);
    graph->GetYaxis()->CenterTitle();
    graph->GetYaxis()->SetTitleSize(0.05);
    graph->GetXaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetLabelSize(0.04);
    graph->GetYaxis()->SetTitleOffset(1.5); // Adjust if necessary

    // Save the canvas as a PNG file
    c1->SaveAs("output/ek-X_misid.png");
}

// void misIDPlot() {
//     // Open the ROOT file and get the tree
//     TFile *file = new TFile("/scratch/thayward/ratios/epi+X_inb.root");
//     TTree *tree = (TTree*)file->Get("PhysicsEvents");

//     // Create a histogram for the fraction calculation
//     TH1F *hFraction = new TH1F("hFraction", ";p (GeV);% k^{+} #rightarrow #pi^{+}", 12, 0, 7);
//     TH1F *hFraction2 = new TH1F("hFraction2", ";p (GeV);% p #rightarrow #pi^{+}", 12, 0, 7);

//     // Loop over the tree and fill the histogram
//     double p_p;
//     int matching_p1_pid;
//     tree->SetBranchAddress("p_p", &p_p);
//     tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

//     for (int i = 0; i < tree->GetEntries(); ++i) {
//         tree->GetEntry(i);
//         if (matching_p1_pid == 321) {
//             hFraction->Fill(p_p);
//         }
//         if (matching_p1_pid == 2212) {
//             hFraction2->Fill(p_p);
//         }
//     }

//     // Normalize the histogram to get the fraction
//     hFraction->Scale(100.0 / tree->GetEntries());
//     hFraction2->Scale(100.0 / tree->GetEntries());

//     // After scaling the histograms
//     double maxVal1 = hFraction->GetMaximum();
//     double maxVal2 = hFraction2->GetMaximum();
//     double maxVal = maxVal1 > maxVal2 ? maxVal1 : maxVal2;

//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph = new TGraphErrors();
//     for (int i = 1; i <= hFraction->GetNbinsX(); ++i) {
//         cout << hFraction->GetBinContent(i) << endl;
//         if (hFraction->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph->SetPoint(i-1, hFraction->GetBinCenter(i), hFraction->GetBinContent(i));
//         graph->SetPointError(i-1, 0, hFraction->GetBinError(i)); // Set horizontal error to 0
//     }
//     cout << endl << endl << endl;
//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph2 = new TGraphErrors();
//     for (int i = 1; i <= hFraction2->GetNbinsX(); ++i) {
//         cout << hFraction2->GetBinContent(i) << endl;
//         if (hFraction2->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph2->SetPoint(i-1, hFraction2->GetBinCenter(i), hFraction2->GetBinContent(i));
//         graph2->SetPointError(i-1, 0, hFraction2->GetBinError(i));
//     }

//     // Style the graphs with markers and colors
//     graph->SetMarkerStyle(20);
//     graph->SetMarkerSize(1.2);
//     graph->SetMarkerColor(kRed);  // Red for matching_p1_pid == 321
//     graph->SetLineColor(kRed);

//     graph2->SetMarkerStyle(21);
//     graph2->SetMarkerSize(1.2);
//     graph2->SetMarkerColor(kBlue);  // Blue for matching_p1_pid == 2212
//     graph2->SetLineColor(kBlue);

//     // Style the graph with markers
//     graph->SetTitle(";p (GeV);% h #rightarrow #pi^{+}"); // Setting title and axis labels

//     graph->SetMaximum(maxVal * 1.1); // 10% higher than the max value for better visibility
//     graph->SetMinimum(0); // Assuming you want to start from 0

//     // Draw the plot using TGraphErrors
//     TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
//     c1->SetLeftMargin(0.15);
//     // Create an invisible histogram to set the axis range
//     TH1F *frame = new TH1F("frame", ";p (GeV);% h #rightarrow #pi^{+}", 12, 0, 7);
//     frame->SetMaximum(maxVal * 1.1); // Set the maximum y-value
//     frame->SetMinimum(0); // Set the minimum y-value
//     frame->SetStats(0); // No statistics box

//     // Draw the invisible histogram to set up the axis
//     frame->Draw();

//     // Draw the graphs on the same canvas
//     graph->Draw("P SAME"); // Draw the first graph
//     graph2->Draw("P SAME"); // Draw the second graph

//     // Set axis styles
//     graph->GetXaxis()->CenterTitle();
//     graph->GetXaxis()->SetTitleSize(0.05);
//     graph->GetYaxis()->CenterTitle();
//     graph->GetYaxis()->SetTitleSize(0.05);
//     graph->GetXaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetTitleOffset(1.5); // Adjust if necessary

//     // Adding a legend
//     TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9); // Adjust the position as needed
//     legend->AddEntry(graph, "k^{+} #rightarrow #pi^{+}", "lp");
//     legend->AddEntry(graph2, "p #rightarrow #pi^{+}", "lp");
//     legend->Draw();

//     // Save the canvas as a PNG file
//     c1->SaveAs("output/epi+X_misid.png");
// }

// void misIDPlot() {
//     // Open the ROOT file and get the tree
//     TFile *file = new TFile("/scratch/thayward/ratios/ek+X_inb.root");
//     TTree *tree = (TTree*)file->Get("PhysicsEvents");

//     // Create a histogram for the fraction calculation
//     TH1F *hFraction = new TH1F("hFraction", ";p (GeV);% #pi^{+} #rightarrow k^{+}", 12, 0, 7);
//     TH1F *hFraction2 = new TH1F("hFraction2", ";p (GeV);% p #rightarrow k^{+}", 12, 0, 7);

//     // Loop over the tree and fill the histogram
//     double p_p;
//     int matching_p1_pid;
//     tree->SetBranchAddress("p_p", &p_p);
//     tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

//     for (int i = 0; i < tree->GetEntries(); ++i) {
//         tree->GetEntry(i);
//         if (matching_p1_pid == 211) {
//             hFraction->Fill(p_p);
//         }
//         if (matching_p1_pid == 2212) {
//             hFraction2->Fill(p_p);
//         }
//     }

//     // Normalize the histogram to get the fraction
//     hFraction->Scale(100.0 / tree->GetEntries());
//     hFraction2->Scale(100.0 / tree->GetEntries());

//     // After scaling the histograms
//     double maxVal1 = hFraction->GetMaximum();
//     double maxVal2 = hFraction2->GetMaximum();
//     double maxVal = maxVal1 > maxVal2 ? maxVal1 : maxVal2;

//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph = new TGraphErrors();
//     for (int i = 1; i <= hFraction->GetNbinsX(); ++i) {
//         cout << hFraction->GetBinContent(i) << endl;
//         if (hFraction->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph->SetPoint(i-1, hFraction->GetBinCenter(i), hFraction->GetBinContent(i));
//         graph->SetPointError(i-1, 0, hFraction->GetBinError(i)); // Set horizontal error to 0
//     }
//     cout << endl << endl << endl;
//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph2 = new TGraphErrors();
//     for (int i = 1; i <= hFraction2->GetNbinsX(); ++i) {
//         cout << hFraction2->GetBinContent(i) << endl;
//         if (hFraction2->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph2->SetPoint(i-1, hFraction2->GetBinCenter(i), hFraction2->GetBinContent(i));
//         graph2->SetPointError(i-1, 0, hFraction2->GetBinError(i));
//     }

//     // Style the graphs with markers and colors
//     graph->SetMarkerStyle(20);
//     graph->SetMarkerSize(1.2);
//     graph->SetMarkerColor(kRed);  // Red for matching_p1_pid == 321
//     graph->SetLineColor(kRed);

//     graph2->SetMarkerStyle(21);
//     graph2->SetMarkerSize(1.2);
//     graph2->SetMarkerColor(kBlue);  // Blue for matching_p1_pid == 2212
//     graph2->SetLineColor(kBlue);

//     // Style the graph with markers
//     graph->SetTitle(";p (GeV);% h #rightarrow k^{+}"); // Setting title and axis labels

//     graph->SetMaximum(maxVal * 1.1); // 10% higher than the max value for better visibility
//     graph->SetMinimum(0); // Assuming you want to start from 0

//     // Draw the plot using TGraphErrors
//     TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
//     c1->SetLeftMargin(0.15);
//     // Create an invisible histogram to set the axis range
//     TH1F *frame = new TH1F("frame", ";p (GeV);% h #rightarrow #pi^{+}", 12, 0, 7);
//     frame->SetMaximum(maxVal * 1.1); // Set the maximum y-value
//     frame->SetMinimum(0); // Set the minimum y-value
//     frame->SetStats(0); // No statistics box

//     // Draw the invisible histogram to set up the axis
//     frame->Draw();

//     // Draw the graphs on the same canvas
//     graph->Draw("P SAME"); // Draw the first graph
//     graph2->Draw("P SAME"); // Draw the second graph

//     // Set axis styles
//     graph->GetXaxis()->CenterTitle();
//     graph->GetXaxis()->SetTitleSize(0.05);
//     graph->GetYaxis()->CenterTitle();
//     graph->GetYaxis()->SetTitleSize(0.05);
//     graph->GetXaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetTitleOffset(1.5); // Adjust if necessary

//     // Adding a legend
//     TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9); // Adjust the position as needed
//     legend->AddEntry(graph, "#pi^{+} #rightarrow k^{+}", "lp");
//     legend->AddEntry(graph2, "p #rightarrow k^{+}", "lp");
//     legend->Draw();

//     // Save the canvas as a PNG file
//     c1->SaveAs("output/ek+X_misid.png");
// }

// void misIDPlot() {
//     // Open the ROOT file and get the tree
//     TFile *file = new TFile("/scratch/thayward/ratios/epX_inb.root");
//     TTree *tree = (TTree*)file->Get("PhysicsEvents");

//     // Create a histogram for the fraction calculation
//     TH1F *hFraction = new TH1F("hFraction", ";p (GeV);% #pi^{+} #rightarrow p", 12, 0, 7);
//     TH1F *hFraction2 = new TH1F("hFraction2", ";p (GeV);% k^{+} #rightarrow p", 12, 0, 7);

//     // Loop over the tree and fill the histogram
//     double p_p;
//     int matching_p1_pid;
//     tree->SetBranchAddress("p_p", &p_p);
//     tree->SetBranchAddress("matching_p1_pid", &matching_p1_pid);

//     for (int i = 0; i < tree->GetEntries(); ++i) {
//         tree->GetEntry(i);
//         if (matching_p1_pid == 211) {
//             hFraction->Fill(p_p);
//         }
//         if (matching_p1_pid == 321) {
//             hFraction2->Fill(p_p);
//         }
//     }

//     // Normalize the histogram to get the fraction
//     hFraction->Scale(100.0 / tree->GetEntries());
//     hFraction2->Scale(100.0 / tree->GetEntries());

//     // After scaling the histograms
//     double maxVal1 = hFraction->GetMaximum();
//     double maxVal2 = hFraction2->GetMaximum();
//     double maxVal = maxVal1 > maxVal2 ? maxVal1 : maxVal2;

//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph = new TGraphErrors();
//     for (int i = 1; i <= hFraction->GetNbinsX(); ++i) {
//         cout << hFraction->GetBinContent(i) << endl;
//         if (hFraction->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph->SetPoint(i-1, hFraction->GetBinCenter(i), hFraction->GetBinContent(i));
//         graph->SetPointError(i-1, 0, hFraction->GetBinError(i)); // Set horizontal error to 0
//     }
//     cout << endl << endl << endl;
//     // Create a TGraphErrors from the histogram with only vertical error bars
//     TGraphErrors *graph2 = new TGraphErrors();
//     for (int i = 1; i <= hFraction2->GetNbinsX(); ++i) {
//         cout << hFraction2->GetBinContent(i) << endl;
//         if (hFraction2->GetBinContent(i) == 0) {
//             continue;
//         }
//         graph2->SetPoint(i-1, hFraction2->GetBinCenter(i), hFraction2->GetBinContent(i));
//         graph2->SetPointError(i-1, 0, hFraction2->GetBinError(i));
//     }

//     // Style the graphs with markers and colors
//     graph->SetMarkerStyle(20);
//     graph->SetMarkerSize(1.2);
//     graph->SetMarkerColor(kRed);  // Red for matching_p1_pid == 321
//     graph->SetLineColor(kRed);

//     graph2->SetMarkerStyle(21);
//     graph2->SetMarkerSize(1.2);
//     graph2->SetMarkerColor(kBlue);  // Blue for matching_p1_pid == 2212
//     graph2->SetLineColor(kBlue);

//     // Style the graph with markers
//     graph->SetTitle(";p (GeV);% h #rightarrow p"); // Setting title and axis labels

//     graph->SetMaximum(maxVal * 1.1); // 10% higher than the max value for better visibility
//     graph->SetMinimum(0); // Assuming you want to start from 0

//     // Draw the plot using TGraphErrors
//     TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
//     c1->SetLeftMargin(0.15);
//     // Create an invisible histogram to set the axis range
//     TH1F *frame = new TH1F("frame", ";p (GeV);% h #rightarrow #pi^{+}", 12, 0, 7);
//     frame->SetMaximum(maxVal * 1.1); // Set the maximum y-value
//     frame->SetMinimum(0); // Set the minimum y-value
//     frame->SetStats(0); // No statistics box

//     // Draw the invisible histogram to set up the axis
//     frame->Draw();

//     // Draw the graphs on the same canvas
//     graph->Draw("P SAME"); // Draw the first graph
//     graph2->Draw("P SAME"); // Draw the second graph

//     // Set axis styles
//     graph->GetXaxis()->CenterTitle();
//     graph->GetXaxis()->SetTitleSize(0.05);
//     graph->GetYaxis()->CenterTitle();
//     graph->GetYaxis()->SetTitleSize(0.05);
//     graph->GetXaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetLabelSize(0.04);
//     graph->GetYaxis()->SetTitleOffset(1.5); // Adjust if necessary

//     // Adding a legend
//     TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9); // Adjust the position as needed
//     legend->AddEntry(graph, "#pi^{+} #rightarrow p", "lp");
//     legend->AddEntry(graph2, "k^{+} #rightarrow p", "lp");
//     legend->Draw();

//     // Save the canvas as a PNG file
//     c1->SaveAs("output/epX_misid.png");
// }

// This allows the script to be run standalone from the ROOT interpreter
void misID() {
    misIDPlot();
}
