#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <iostream>
#include <iomanip> // For std::setprecision
#include <sstream> // For std::stringstream

// Using namespace declaration
using namespace std;

void analyzePions() {
    // Open the ROOT file
    TFile* file = TFile::Open("/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/mc/epi+X/rga_fa18_inb_clasdis_50nA_epi+X.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }
    
    // Access the tree
    TTree* tree = dynamic_cast<TTree*>(file->Get("PhysicsEvents"));
    if (!tree) {
        std::cerr << "Error: Tree 'PhysicsEvents' not found" << std::endl;
        file->Close();
        return;
    }

    // Define variables to hold branch data
    Double_t Q2, W, y, xF, Mx, z, p_p, pT, phi;
    Int_t mc_p1_parent;

    // Set branch addresses
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("xF", &xF);
    tree->SetBranchAddress("Mx", &Mx);
    tree->SetBranchAddress("z", &z);
    tree->SetBranchAddress("p_p", &p_p);
    tree->SetBranchAddress("pT", &pT);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("mc_p1_parent", &mc_p1_parent);

    // Containers for unique mc_p1_parent values and event counts
    std::map<int, int> parentEventCounts;
    int totalEventsMeetingCriteria = 0;

    // Define pT bins
    std::vector<double> pTBins = {0.05, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00};

    // Create a vector to hold histograms for each pT bin
    std::vector<TH1F*> histograms;

    // Create histograms for each pT bin
    for (size_t i = 0; i < pTBins.size() - 1; ++i) {
        std::stringstream histName, histTitle;
        histName << "hist" << i;
        histTitle << std::fixed << std::setprecision(2) << pTBins[i] << " < P_T (GeV) < " << pTBins[i+1];
        histograms.push_back(new TH1F(histName.str().c_str(), histTitle.str().c_str(), 75, 0, 2*M_PI));
        histograms.back()->SetLineColor(kBlack); // Make the plot black
    }

    // Map of parent pids to vectors of histograms
    std::map<int, std::vector<TH1F*>> pidHistograms;
    std::vector<int> pids = {92, 213, 223, 113, 310, 323}; // List of specific pids

    // Initialize histograms for each pid and pT bin
    for (int pid : pids) {
        std::vector<TH1F*> pidHists;
        for (size_t i = 0; i < pTBins.size() - 1; ++i) {
            std::stringstream histName, histTitle;
            histName << "hist_" << pid << "_" << i;
            histTitle << "PID " << pid << ": " << std::fixed << std::setprecision(2) << pTBins[i] << " < P_T (GeV) < " << pTBins[i+1];
            pidHists.push_back(new TH1F(histName.str().c_str(), histTitle.str().c_str(), 75, 0, 2*M_PI));
        }
        pidHistograms[pid] = pidHists;
    }


    // Set global style options
    gStyle->SetOptStat(0); // Hide the statistics box
    gStyle->SetTextFont(42); // Set font
    gStyle->SetLabelSize(0.04, "XY"); // Set axis label size
    gStyle->SetTitleSize(0.05, "XY"); // Set axis title size
    gStyle->SetTitleOffset(1.2, "Y"); // Adjust Y axis title position

    // Create a 2x4 canvas
    TCanvas* canvas = new TCanvas("canvas", "Phi Distributions", 1200, 800);
    canvas->Divide(4, 2); // 4 columns, 2 rows


    // Loop through the tree
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        // Apply the kinematic conditions
        if (Q2 > 1 && W > 2 && y < 0.75 && xF > 0 && Mx > 1.5 && z > 0.2 && p_p > 1.2 &&
            y > 0.65 && y < 0.75 && z > 0.24 && z < 0.29 && Q2 > 2.0 && Q2 < 2.5) {
            // Increment count for this mc_p1_parent
            parentEventCounts[mc_p1_parent]++;
            totalEventsMeetingCriteria++;

            for (size_t bin = 0; bin < pTBins.size() - 1; ++bin) {
                if (pT > pTBins[bin] && pT <= pTBins[bin + 1]) {
                    histograms[bin]->Fill(phi);

                    // Check for specific pids and fill their histograms
                    if (pidHistograms.find(mc_p1_parent) != pidHistograms.end()) {
                        // std::cout << mc_p1_parent << endl;
                        pidHistograms[mc_p1_parent][bin]->Fill(phi);
                    }
                }
            }
        }
    }

    // Print unique mc_p1_parent values and their corresponding percentages
    std::cout << "mc_p1_parent values with >1.0% of events:" << std::endl;
    for (const auto& pair : parentEventCounts) {
        double percentage = 100.0 * pair.second / totalEventsMeetingCriteria;
        if (percentage > 1.0) { // Check if percentage is greater than 0.1%
            std::cout << "mc_p1_parent = " << pair.first << ": " << percentage << "%" << std::endl;
        }
    }

    for (int pid : pids) {
        for (TH1F* hist : pidHistograms[pid]) {
            switch (pid) {
                case 92: // Dashed blue lines
                    hist->SetLineColor(kBlue);
                    hist->SetLineStyle(2); // Dashed
                    break;
                case 213: // Dashed dotted dark red lines
                    hist->SetLineColor(kRed+2);
                    hist->SetLineStyle(5); // Dashed-dotted
                    break;
                case 223: // Double dashed dark red lines
                    hist->SetLineColor(kRed+2);
                    hist->SetLineStyle(3); // Double dashed
                    break;
                case 113: // Dotted green lines
                    hist->SetLineColor(kGreen+2);
                    hist->SetLineStyle(4); // Dotted
                    break;
                case 310: // Dotted orange lines
                    hist->SetLineColor(kOrange+2);
                    hist->SetLineStyle(2); // Dashed
                    break;
                case 323: // Dotted green lines
                    hist->SetLineColor(kYellow+2);
                    hist->SetLineStyle(3); // Double dashed
                    break;
            }
        }
    }


    for (size_t i = 0; i < histograms.size(); ++i) {
        canvas->cd(i + 1);
        // Adjust pad margins
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.05);

        histograms[i]->GetXaxis()->SetTitle("#phi");
        histograms[i]->GetYaxis()->SetTitle("Counts");
        histograms[i]->SetMinimum(0); // Set y-axis to start at 0
        histograms[i]->Draw("HIST"); // Draw the main histogram first

        // Draw additional pid histograms
        for (int pid : pids) {
            pidHistograms[pid][i]->Draw("HIST SAME");
        }

        // Adjustments, labels, etc.
    }

    // Move to the empty eighth pad
    canvas->cd(8);

    // Create a legend in the eighth pad
    TLegend* legend = new TLegend(0.1, 0.1, 0.9, 0.9);
    legend->SetTextSize(0.04); // Adjust text size as necessary

    // Add an entry for the main histogram (assuming similar style for all)
    legend->AddEntry(histograms[0], "Main Histogram", "l");

    // Add entries for each specific PID with a description and the corresponding line style
    legend->AddEntry(pidHistograms[92][0], "Direct", "l");
    legend->AddEntry(pidHistograms[213][0], "#rho^{0}", "l");
    legend->AddEntry(pidHistograms[223][0], "#rho^{+1}", "l");
    legend->AddEntry(pidHistograms[113][0], "#omega", "l");
    legend->AddEntry(pidHistograms[310][0], "K_{0}^{s}", "l");
    legend->AddEntry(pidHistograms[323][0], "K^{*+}(892)", "l");

    // You might need to draw a representative line for each style in the legend
    // This is handled automatically by specifying the correct options ("l") and having the histograms styled accordingly

    legend->Draw();


    // Save the canvas to a file
    canvas->SaveAs("output/q2y-1_z3.png");
    // Cleanup
    file->Close();
}

int main() {
    analyzePions();
    return 0;
}
