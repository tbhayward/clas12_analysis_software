#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>

const std::string output_dir = "output/rgc_ready_for_cooking_plots/";

// Function to create a histogram for a given dataset
TH1D* createHistogram(TTree* tree, const char* name, const char* title, 
        const char* variable, const char* cut, double norm) {
    TH1D* hist = new TH1D(name, title, 80, -4, 4);
    tree->Draw((std::string(variable) + ">>" + name).c_str(), cut, "goff");
    hist->Scale(1.0 / norm);
    hist->SetStats(kFALSE);
    return hist;
}

void rgc_preparation() {
    // File paths and tree names
    const char* rga_eX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_eX.root";
    const char* rga_epipX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+X.root";
    const char* rga_epX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epX.root";
    const char* rga_epippimX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+pi-X.root";

    const char* rgc_eX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_eX.root";
    const char* rgc_epipX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+X.root";
    const char* rgc_epX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epX.root";
    const char* rgc_epippimX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+pi-X.root";

    // Open files and get trees
    TFile* filesOpened[8];
    TTree* trees[8];
    for (int i = 0; i < 8; i++) {
        filesOpened[i] = new TFile(files[i]);
        trees[i] = (TTree*)filesOpened[i]->Get("PhysicsEvents");
    }

    double rga_H2_norm = 159661.55+145813.73;
    double rgc_pos_NH3_norm = 19355.9+19392.53+21683.25+21621.178;
    double rgc_neg_NH3_norm = 21282.264+21217.414+21303.576+21297.766;
    double rgc_C_norm = 8883.014+8834.256;

    // Example: Initialize a canvas for plotting for one of the datasets
    TCanvas *c1 = new TCanvas("c1", "Data Analysis", 800, 600);

    // More analysis and plotting logic here

    // Example: Save the canvas to a file for one of the datasets
    std::string output_file = output_dir + "analysis_output_example.png";
    c1->SaveAs(output_file.c_str());

    // Clean up
    delete c1;
    rga_eX->Close(); rga_epipX->Close(); rga_epX->Close(); rga_epippimX->Close();
    rgc_eX->Close(); rgc_epipX->Close(); rgc_epX->Close(); rgc_epippimX->Close();
    delete rga_eX; delete rga_epipX; delete rga_epX; delete rga_epippimX;
    delete rgc_eX; delete rgc_epipX; delete rgc_epX; delete rgc_epippimX;
}

int main() {
    rgc_preparation();
    return 0;
}
