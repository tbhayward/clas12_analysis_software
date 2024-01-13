#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TGraph.h>
#include <iostream>

const std::string output_dir = "output/";

void rgc_preparation(const char* inputFile1, const char* inputFile2, const char* outputFile) {
    // File paths and tree names
    const char* rga_eX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_eX.root";
    const char* rga_epipX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+X.root";
    const char* rga_epX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epX.root";
    const char* rga_epippimX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rga_ready_for_calibration_epi+pi-X.root";

    const char* rgc_eX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_eX.root";
    const char* rgc_epipX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+X.root";
    const char* rgc_epX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epX.root";
    const char* rgc_epippimX_file = "/volatile/clas12/thayward/rgc_ready_for_cooking/processed_files/rgc_ready_for_calibration_epi+pi-X.root";

    // Loading files
    TFile* rga_eX = new TFile(rga_eX_file);
    TFile* rga_epipX = new TFile(rga_epipX_file);
    TFile* rga_epX = new TFile(rga_epX_file);
    TFile* rga_epippimX = new TFile(rga_epippimX_file);

    TFile* rgc_eX = new TFile(rgc_eX_file);
    TFile* rgc_epipX = new TFile(rgc_epipX_file);
    TFile* rgc_epX = new TFile(rgc_epX_file);
    TFile* rgc_epippimX = new TFile(rgc_epippimX_file);

    // Extracting trees
    TTree* tree_rga_eX = (TTree*)rga_eX->Get("PhysicsEvents");
    TTree* tree_rga_epipX = (TTree*)rga_epipX->Get("PhysicsEvents");
    TTree* tree_rga_epX = (TTree*)rga_epX->Get("PhysicsEvents");
    TTree* tree_rga_epippimX = (TTree*)rga_epippimX->Get("PhysicsEvents");

    TTree* tree_rgc_eX = (TTree*)rgc_eX->Get("PhysicsEvents");
    TTree* tree_rgc_epipX = (TTree*)rgc_epipX->Get("PhysicsEvents");
    TTree* tree_rgc_epX = (TTree*)rgc_epX->Get("PhysicsEvents");
    TTree* tree_rgc_epippimX = (TTree*)rgc_epippimX->Get("PhysicsEvents");

    double rga_norm = 159661.55+145813.73;
    double rgc_pos_norm = 19355.9+19392.53+21683.25+21621.178;
    double rgc_neg_norm = 21282.264+21217.414+21303.576+21297.766;
    double rgc_carbon_norm = 8883.014+8834.256;

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
