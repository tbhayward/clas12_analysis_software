#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <iostream>
#include <string>
#include <regex>
#include <unordered_set>
#include <utility>
#include <cmath>

// Helper function to format LaTeX-like input for ROOT titles
std::string formatLatexString(const std::string& input) {
    std::string formatted = input;

    std::regex subscript_pattern("\\_\\{([^}]*)\\}");
    formatted = std::regex_replace(formatted, subscript_pattern, "_{$1}");

    std::regex superscript_pattern("\\^\\{([^}]*)\\}");
    formatted = std::regex_replace(formatted, superscript_pattern, "^{$1}");

    return formatted;
}

// Hash function for std::pair<int, int> to allow use in unordered_set
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ (hash2 << 1); // Combine the two hashes
    }
};

int main(int argc, char** argv) {
    if (argc != 11) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <label1> <label2> <branch_variable> <x_axis_label> <x_min> <x_max> <range_low> <range_high>" << std::endl;
        return 1;
    }

    const char* file1_name = argv[1];
    const char* file2_name = argv[2];
    std::string label1 = argv[3];
    std::string label2 = argv[4];
    const char* branch_name = argv[5];
    std::string x_axis_label = argv[6];
    double x_min = std::stod(argv[7]);
    double x_max = std::stod(argv[8]);
    double range_low = std::stod(argv[9]);
    double range_high = std::stod(argv[10]);

    std::string formatted_label = formatLatexString(x_axis_label);

    TFile* file1 = TFile::Open(file1_name);
    TFile* file2 = TFile::Open(file2_name);
    if (!file1 || !file2 || file1->IsZombie() || file2->IsZombie()) {
        std::cerr << "Error: Could not open one or both ROOT files." << std::endl;
        return 1;
    }

    TTree* tree1 = (TTree*)file1->Get("PhysicsEvents");
    TTree* tree2 = (TTree*)file2->Get("PhysicsEvents");
    if (!tree1 || !tree2) {
        std::cerr << "Error: Could not find tree 'PhysicsEvents' in one or both files." << std::endl;
        return 1;
    }

    // Increase bins by a factor of 1.5 and round to the nearest integer
    int num_bins = static_cast<int>(100 * 1.5);

    TH1D* hist1 = new TH1D("hist1", "", num_bins, x_min, x_max);
    TH1D* hist2 = new TH1D("hist2", "", num_bins, x_min, x_max);

    tree1->Project("hist1", branch_name);
    tree2->Project("hist2", branch_name);

    double integral1 = hist1->Integral();
    if (integral1 > 0) {
        hist1->Scale(1.0 / integral1); // Normalize hist1 to 1
        hist2->Scale(1.0 / integral1); // Normalize hist2 by the same factor as hist1
    }

    hist1->SetLineColor(kBlue);
    hist2->SetLineColor(kRed);
    hist1->GetXaxis()->SetTitle(formatted_label.c_str());
    hist1->GetYaxis()->SetTitle("normalized counts"); // Lowercase title
    hist2->GetXaxis()->SetTitle(formatted_label.c_str());
    hist2->GetYaxis()->SetTitle("normalized counts");

    hist1->SetStats(0);
    hist2->SetStats(0);

    // Step 1: Find (runnum, evnum) pairs within the specified range in file2
    std::unordered_set<std::pair<int, int>, pair_hash> matching_event_pairs;
    double branch_value;
    int runnum, evnum;
    tree2->SetBranchAddress(branch_name, &branch_value);
    tree2->SetBranchAddress("runnum", &runnum);
    tree2->SetBranchAddress("evnum", &evnum);

    Long64_t nEntries2 = tree2->GetEntries();
    for (Long64_t i = 0; i < nEntries2; ++i) {
        tree2->GetEntry(i);
        if (branch_value >= range_low && branch_value <= range_high) {
            matching_event_pairs.emplace(runnum, evnum);
        }
    }

    // Step 2: Record positions of matching events in file1
    TH1D* hist3 = new TH1D("hist3", "", num_bins, x_min, x_max);
    tree1->SetBranchAddress("runnum", &runnum);
    tree1->SetBranchAddress("evnum", &evnum);
    tree1->SetBranchAddress(branch_name, &branch_value);

    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        if (matching_event_pairs.find({runnum, evnum}) != matching_event_pairs.end()) {
            hist3->Fill(branch_value);
        }
    }

    // Normalize hist3 by the integral of hist1
    if (integral1 > 0) {
        hist3->Scale(1.0 / integral1);
    }
    hist3->SetLineColor(kBlack);
    hist3->SetLineStyle(2); // Dashed line for third histogram
    hist3->SetStats(0);

    // Create a canvas with extra left margin padding
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    canvas->SetLeftMargin(0.125); // User-preferred padding

    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");
    hist3->Draw("HIST SAME"); // Draw third histogram on the same canvas

    // Use the new input arguments for legend labels
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, label1.c_str(), "l");
    legend->AddEntry(hist2, label2.c_str(), "l");
    legend->AddEntry(hist3, "Selected Events", "l");
    legend->Draw();

    canvas->SaveAs("output/rad_study/test.pdf");

    delete canvas;
    delete hist1;
    delete hist2;
    delete hist3;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    return 0;
}