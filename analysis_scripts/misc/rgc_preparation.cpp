#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>

void rgc_preparation(const char* inputFile1, const char* inputFile2, const char* outputFile) {
    // Load input files
    TFile *file1 = new TFile(inputFile1);
    TFile *file2 = new TFile(inputFile2);
    if (!file1->IsOpen() || !file2->IsOpen()) {
        std::cerr << "Error opening files!" << std::endl;
        return;
    }

    TTree *tree1 = (TTree*)file1->Get("PhysicsEvents"); // Replace with your actual tree name
    TTree *tree2 = (TTree*)file2->Get("PhysicsEvents"); // Replace with your actual tree name

    double rga_norm = 159661.55+145813.73;
    double rgc_pos_norm = 19355.9+19392.53+21683.25+21621.178;
    double rgc_neg_norm = 21282.264+21217.414+21303.576+21297.766;
    double rgc_carbon_norm = 8883.014+8834.256;

    // Initialize the canvas for plotting
    TCanvas *c1 = new TCanvas("c1", "Data Analysis", 800, 600);

    // Placeholder for actual data analysis and plotting
    // This is where we will call our fitting function (to be implemented)

    // Save the canvas to a file
    c1->SaveAs(outputFile);

    // Clean up
    delete c1;
    file1->Close();
    file2->Close();
    delete file1;
    delete file2;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <InputFile1> <InputFile2> <OutputFile>" << std::endl;
        return 1;
    }
    AnalyzeData(argv[1], argv[2], argv[3]);
    return 0;
}
