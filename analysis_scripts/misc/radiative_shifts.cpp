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

// Helper function to get the file name without path or extension
std::string getFileNameWithoutExtension(const std::string& filePath) {
    size_t lastSlashPos = filePath.find_last_of("/");
    std::string fileName = (lastSlashPos == std::string::npos) ? filePath : filePath.substr(lastSlashPos + 1);
    size_t lastDotPos = fileName.find_last_of(".");
    return (lastDotPos == std::string::npos) ? fileName : fileName.substr(0, lastDotPos);
}

int main(int argc, char** argv) {
    if (argc != 12) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <label1> <label2> <branch_variable> <x_axis_label> <x_min> <x_max> <range_low> <range_high> <ratio_lower_bound>" << std::endl;
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
    double ratio_lower_bound = std::stod(argv[11]);

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

    // Step 2: Record positions of matching events in file1 starting from ratio_lower_bound
    TH1D* hist3 = new TH1D("hist3", "", num_bins, x_min, x_max);
    tree1->SetBranchAddress("runnum", &runnum);
    tree1->SetBranchAddress("evnum", &evnum);
    tree1->SetBranchAddress(branch_name, &branch_value);

    Long64_t nEntries1 = tree1->GetEntries();
    for (Long64_t i = 0; i < nEntries1; ++i) {
        tree1->GetEntry(i);
        if (branch_value >= ratio_lower_bound && matching_event_pairs.find({runnum, evnum}) != matching_event_pairs.end()) {
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

    // Set x-axis range for hist3 in the main plot
    hist3->GetXaxis()->SetRangeUser(ratio_lower_bound, x_max);

    // Calculate the maximum y-axis range
    double max_val = std::max({hist1->GetMaximum(), hist2->GetMaximum(), hist3->GetMaximum()});
    hist1->SetMaximum(1.2 * max_val); // Set y-axis range

    // Create a canvas with extra left margin padding
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    canvas->SetLeftMargin(0.125); // User-preferred padding

    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");
    hist3->Draw("HIST SAME"); // Draw third histogram on the same canvas with restricted range

    // Use the new input arguments for legend labels
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, label1.c_str(), "l");
    legend->AddEntry(hist2, label2.c_str(), "l");
    legend->AddEntry(hist3, "Shifted Rad Events", "l"); // Updated label
    legend->Draw();

    // Construct output filename based on the second file's name and branch variable
    std::string file2_identifier = getFileNameWithoutExtension(file2_name);
    std::string output_filename = "output/rad_study/" + file2_identifier + "_" + branch_name + ".pdf";
    canvas->SaveAs(output_filename.c_str());

    // Create the ratio plot of hist3 / hist1 within the specified bounds
    TH1D* ratio_hist = (TH1D*)hist3->Clone("ratio_hist");
    ratio_hist->Divide(hist1);
    ratio_hist->GetXaxis()->SetRangeUser(ratio_lower_bound, x_max); // Set x-axis range for ratio
    ratio_hist->GetYaxis()->SetTitle("rad events / nominal events");
    ratio_hist->GetXaxis()->SetTitle(formatted_label.c_str());

    // Set the y-axis maximum to 1.35 times the maximum value in the ratio histogram
    double ratio_max_val = ratio_hist->GetMaximum();
    ratio_hist->SetMaximum(1.35 * ratio_max_val);

    TCanvas* ratio_canvas = new TCanvas("ratio_canvas", "", 800, 600);
    ratio_canvas->SetLeftMargin(0.125); // User-preferred padding for ratio plot
    ratio_hist->Draw("HIST");

    // Save the ratio plot with an additional "_ratio" suffix
    std::string ratio_output_filename = "output/rad_study/" + file2_identifier + "_" + branch_name + "_ratio.pdf";
    ratio_canvas->SaveAs(ratio_output_filename.c_str());

    delete canvas;
    delete ratio_canvas;
    delete hist1;
    delete hist2;
    delete hist3;
    delete ratio_hist;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    return 0;
}