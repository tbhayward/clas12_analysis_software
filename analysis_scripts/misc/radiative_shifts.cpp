#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <TTreeReader.h>
#include <iostream>
#include <string>
#include <regex>
#include <unordered_set>
#include <utility>
#include <cmath>

// Global Q² and y range variables
// const double Q2_MIN = -2.0;
// const double Q2_MAX = 12.0;
// const double y_MIN = 0.0;
// const double y_MAX = 1.0;
// const double W_MIN = -2;
// const double W_MAX = 10;
// const double Mx2_MIN = -12;
// const double Mx2_MAX = 15;

const double Q2_MIN = 1.0;
const double Q2_MAX = 12.0;
const double y_MIN = 0.0;
const double y_MAX = 0.80;
const double W_MIN = 2;
const double W_MAX = 10;
const double Mx2_MIN = 1.8225;
const double Mx2_MAX = 15;

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
    if (argc != 12 && argc != 14) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <label1> <label2> <branch_variable> <x_axis_label> <x_min> <x_max> <range_low> <range_high> <ratio_lower_bound> [range_low2 range_high2]" << std::endl;
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

    double range_low2 = 0, range_high2 = 0;
    bool has_second_region = (argc == 14);
    if (has_second_region) {
        range_low2 = std::stod(argv[12]);
        range_high2 = std::stod(argv[13]);
    }

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

    // Create expressions for filtering based on Q² and y ranges
    std::string filter_expr = Form("Q2 >= %f && Q2 <= %f && y >= %f && y <= %f && W >= %f && W <= %f && Mx2 >= %f && Mx2 <= %f", Q2_MIN, Q2_MAX, y_MIN, y_MAX, W_MIN, W_MAX, Mx2_MIN, Mx2_MAX);

    // Create histograms for tree projections with the filtering conditions
    TH1D* hist1 = new TH1D("hist1", "", static_cast<int>(100 * 1.5), x_min, x_max);
    TH1D* hist2 = new TH1D("hist2", "", static_cast<int>(100 * 1.5), x_min, x_max);
    tree1->Project("hist1", branch_name, filter_expr.c_str());
    tree2->Project("hist2", branch_name, filter_expr.c_str());

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

    // Step 1: Set up TTreeReader
    TTreeReader reader(tree2);
    TTreeReaderValue<double> W_reader(reader, "W");
    TTreeReaderValue<double> Q2_reader(reader, "Q2");
    TTreeReaderValue<double> y_reader(reader, "y");
    TTreeReaderValue<double> Mx2_reader(reader, "Mx2");
    TTreeReaderValue<double> branch_reader(reader, branch_name); // Dynamically handle branch_name
    TTreeReaderValue<int> runnum_reader(reader, "runnum");
    TTreeReaderValue<int> evnum_reader(reader, "evnum");

    // Step 2: Iterate over entries in tree2
    std::unordered_set<std::pair<int, int>, pair_hash> matching_event_pairs1, matching_event_pairs2;

    while (reader.Next()) {
        double W_filter = *W_reader;
        double Q2_filter = *Q2_reader;
        double y_filter = *y_reader;
        double Mx2_filter = *Mx2_reader;
        double branch_data = *branch_reader; // This dynamically links to branch_name
        int run = *runnum_reader;
        int event = *evnum_reader;

        // Apply the required cuts
        bool passW = W_filter >= W_MIN && W_filter <= W_MAX;
        bool passQ2 = Q2_filter >= Q2_MIN && Q2_filter <= Q2_MAX;
        bool passY = y_filter >= y_MIN && y_filter <= y_MAX;
        bool passMx2 = Mx2_filter >= Mx2_MIN && Mx2_filter <= Mx2_MAX;

        // Check if branch_data falls within the specified ranges
        if (branch_data >= range_low && branch_data <= range_high) {
            matching_event_pairs1.emplace(run, event);
        }
        if (has_second_region && branch_data >= range_low2 && branch_data <= range_high2) {
            matching_event_pairs2.emplace(run, event);
        }
    }

    // Step 3: Create histograms for file1 based on matching events
    TTreeReader reader1(tree1);
    TTreeReaderValue<double> branch_reader1(reader1, branch_name); // Dynamically handle branch_name
    TTreeReaderValue<int> runnum_reader1(reader1, "runnum");
    TTreeReaderValue<int> evnum_reader1(reader1, "evnum");
    TTreeReaderValue<double> W_reader1(reader1, "W");
    TTreeReaderValue<double> Q2_reader1(reader1, "Q2");
    TTreeReaderValue<double> y_reader1(reader1, "y");

    TH1D* hist3 = new TH1D("hist3", "", static_cast<int>(100 * 1.5), x_min, x_max);
    TH1D* hist4 = has_second_region ? new TH1D("hist4", "", static_cast<int>(100 * 1.5), x_min, x_max) : nullptr;

    while (reader1.Next()) {
        double branch_data1 = *branch_reader1;
        double W1 = *W_reader1;
        double Q21 = *Q2_reader1;
        double y1 = *y_reader1;
        double Mx2_filter = *Mx2_reader;
        int run1 = *runnum_reader1;
        int event1 = *evnum_reader1;

        // Apply cuts for file1
        bool passW = W1 >= W_MIN && W1 <= W_MAX;
        bool passQ2 = Q21 >= Q2_MIN && Q21 <= Q2_MAX;
        bool passY = y1 >= y_MIN && y1 <= y_MAX;
        bool passMx2 = Mx2_filter >= Mx2_MIN && Mx2_filter <= Mx2_MAX;

        if (passW && passQ2 && passY && passMx2 && branch_data1 >= ratio_lower_bound) {
            if (matching_event_pairs1.find({run1, event1}) != matching_event_pairs1.end()) {
                hist3->Fill(branch_data1);
            }
            if (has_second_region && matching_event_pairs2.find({run1, event1}) != matching_event_pairs2.end()) {
                hist4->Fill(branch_data1);
            }
        }
    }

    // Normalize histograms
    if (integral1 > 0) {
        hist3->Scale(1.0 / integral1);
        if (has_second_region) hist4->Scale(1.0 / integral1);
    }

    // Set histogram colors and styles
    hist3->SetLineColor(kCyan + 2); // Distinct color
    hist3->SetLineStyle(2);         // Dashed line
    if (has_second_region) {
        hist4->SetLineColor(kMagenta + 1); // Distinct color
        hist4->SetLineStyle(2);            // Dashed line
    }

    // Calculate the maximum y-axis range
    double max_val = std::max({hist1->GetMaximum(), hist2->GetMaximum(), hist3->GetMaximum()});
    if (has_second_region) max_val = std::max(max_val, hist4->GetMaximum());
    hist1->SetMaximum(1.2 * max_val);

    // Create the main canvas
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    canvas->SetLeftMargin(0.125);

    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");
    hist3->Draw("HIST SAME");
    if (has_second_region) hist4->Draw("HIST SAME");

    // Add legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, label1.c_str(), "l");
    legend->AddEntry(hist2, label2.c_str(), "l");
    legend->AddEntry(hist3, Form("Shifted Rad Events [%g, %g]", range_low, range_high), "l");
    if (has_second_region) {
        legend->AddEntry(hist4, Form("Shifted Rad Events [%g, %g]", range_low2, range_high2), "l");
    }
    legend->Draw();

    // Save the main plot
    std::string file2_identifier = getFileNameWithoutExtension(file2_name);
    std::string output_filename = "output/rad_study/" + file2_identifier + "_" + branch_name + ".pdf";
    canvas->SaveAs(output_filename.c_str());

    // Create the ratio plot
    TCanvas* ratio_canvas = new TCanvas("ratio_canvas", "", 800, 600);
    ratio_canvas->SetLeftMargin(0.125);

    TH1D* ratio_hist3 = (TH1D*)hist3->Clone("ratio_hist3");
    ratio_hist3->Divide(hist1);
    ratio_hist3->SetLineColor(kCyan + 2);
    ratio_hist3->SetLineStyle(2);
    ratio_hist3->GetYaxis()->SetTitle("rad events / nominal events");
    ratio_hist3->GetXaxis()->SetTitle(formatted_label.c_str());
    ratio_hist3->SetStats(0);

    ratio_hist3->Draw("HIST");
    if (has_second_region) {
        TH1D* ratio_hist4 = (TH1D*)hist4->Clone("ratio_hist4");
        ratio_hist4->Divide(hist1);
        ratio_hist4->SetLineColor(kMagenta + 1);
        ratio_hist4->SetLineStyle(2);
        ratio_hist4->SetStats(0);
        ratio_hist4->Draw("HIST SAME");

        legend = new TLegend(0.7, 0.7, 0.9, 0.9);
        legend->AddEntry(ratio_hist3, Form("Ratio [%g, %g]", range_low, range_high), "l");
        legend->AddEntry(ratio_hist4, Form("Ratio [%g, %g]", range_low2, range_high2), "l");
        legend->Draw();
    }

    // Adjust y-axis range for ratio plot
    ratio_hist3->SetMaximum(1);

    // Save the ratio plot
    std::string ratio_output_filename = "output/rad_study/" + file2_identifier + "_" + branch_name + "_ratio.pdf";
    ratio_canvas->SaveAs(ratio_output_filename.c_str());

    // Cleanup
    delete canvas;
    delete ratio_canvas;
    delete hist1;
    delete hist2;
    delete hist3;
    if (has_second_region) delete hist4;
    delete file1;
    delete file2;

    return 0;
}