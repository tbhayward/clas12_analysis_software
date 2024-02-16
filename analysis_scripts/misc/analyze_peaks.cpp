#include <TH1F.h> // For histograms
#include <TCanvas.h> // For canvas
#include <TPad.h> // For pad manipulation
#include <TStyle.h> // For global style settings
#include <TLatex.h> // For text on plots
#include <iostream>


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
        std::string histName = "hist" + std::to_string(i);
        std::string histTitle = std::to_string(pTBins[i]) + " < P_T (GeV) < " + std::to_string(pTBins[i+1]);
        histograms.push_back(new TH1F(histName.c_str(), histTitle.c_str(), 24, 0, 2*M_PI));
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
                    break;
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

    for (size_t i = 0; i < histograms.size(); ++i) {
        canvas->cd(i + 1); // Move to the next pad
        histograms[i]->GetXaxis()->SetTitle("#phi");
        histograms[i]->GetYaxis()->SetTitle("Counts");
        histograms[i]->Draw("HIST"); // Draw histogram
        
        // Adjust pad margins
        gPad->SetLeftMargin(0.12);
        gPad->SetBottomMargin(0.12);
        gPad->SetRightMargin(0.05);

        // Add pT bin label
        TLatex latex;
        latex.SetTextSize(0.05); // Increase font size
        latex.DrawLatexNDC(0.15, 0.85, histograms[i]->GetTitle()); // Position adjusted for visibility
    }

    // Save the canvas to a file
    canvas->SaveAs("output/q2y-1_z3.png");
    // Cleanup
    file->Close();
    delete file;
}

int main() {
    analyzePions();
    return 0;
}
