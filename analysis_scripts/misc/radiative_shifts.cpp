#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <iostream>
#include <string>
#include <regex>
#include <cmath>

std::string formatLatexString(const std::string& input) {
    std::string formatted = input;

    std::regex subscript_pattern("\\_\\{([^}]*)\\}");
    formatted = std::regex_replace(formatted, subscript_pattern, "_{$1}");

    std::regex superscript_pattern("\\^\\{([^}]*)\\}");
    formatted = std::regex_replace(formatted, superscript_pattern, "^{$1}");

    return formatted;
}

int main(int argc, char** argv) {
    if (argc != 9) {
        std::cerr << "Usage: " << argv[0] << " <file1.root> <file2.root> <label1> <label2> <branch_variable> <x_axis_label> <x_min> <x_max>" << std::endl;
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

    // Create a canvas with extra left margin padding
    TCanvas* canvas = new TCanvas("canvas", "", 800, 600);
    canvas->SetLeftMargin(0.075); // Increase left margin for y-axis label padding

    hist1->Draw("HIST");
    hist2->Draw("HIST SAME");

    // Use the new input arguments for legend labels
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(hist1, label1.c_str(), "l");
    legend->AddEntry(hist2, label2.c_str(), "l");
    legend->Draw();

    canvas->SaveAs("output/rad_study/test.pdf");

    delete canvas;
    delete hist1;
    delete hist2;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    return 0;
}