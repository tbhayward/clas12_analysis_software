#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

// ROOT Libraries
#include <TROOT.h>
#include <TApplication.h>

// tbhayward imports
#include "create_directories.h"
#include "determine_exclusivity.h"

// Namespace declaration
using namespace std;
namespace fs = std::filesystem;

// Function to check if a file exists
bool checkFileExists(const string& path) {
    return fs::exists(path);
}

int main(int argc, char* argv[]) {
    std::cout << std::endl << std::endl << std::endl;
    // Start the ROOT application
    TApplication theApp("App", nullptr, nullptr);
    gROOT->SetBatch(kTRUE);  // Set ROOT to batch mode

    // Ensure that the correct number of command-line arguments is provided
    if (argc < 5) {
        cout << "Usage: " << argv[0] << " <dir1> <dir2> <dir3> <dir4>" << endl;
        return 1;
    }

    // Store directories
    string dir1 = argv[1];
    string dir2 = argv[2];
    string dir3 = argv[3];
    string dir4 = argv[4];

    // Define filenames for each directory (3 periods, 6 files per period)
    vector<string> data_filenames = {dir1 + "/rga_fa18_inb_epgamma.root", dir1 + "/rga_fa18_out_epgamma.root", dir1 + "/rga_sp19_inb_epgamma.root"};
    vector<string> eppi0_filenames = {dir3 + "/rga_fa18_inb_eppi0.root", dir3 + "/rga_fa18_out_eppi0.root", dir3 + "/rga_sp19_inb_eppi0.root"};
    vector<string> mc_gen_dvcsgen_filenames = {dir2 + "/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root", dir2 + "/gen_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root", dir2 + "/gen_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"};
    vector<string> mc_rec_dvcsgen_filenames = {dir2 + "/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root", dir2 + "/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root", dir2 + "/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"};
    vector<string> mc_gen_aaogen_filenames = {dir4 + "/gen_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root", dir4 + "/gen_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root", dir4 + "/gen_aaogen_norad_rga_sp19_inb_50nA_10604MeV_eppi0.root"};
    vector<string> mc_rec_aaogen_filenames = {dir4 + "/rec_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root", dir4 + "/rec_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root", dir4 + "/rec_aaogen_norad_rga_sp19_inb_50nA_10604MeV_eppi0.root"};

    // Check if all expected files exist
    for (const auto& file : data_filenames) {
        if (!checkFileExists(file)) {
            cerr << "Error: File " << file << " not found." << endl;
            return 2;
        }
    }
    for (const auto& file : eppi0_filenames) {
        if (!checkFileExists(file)) {
            cerr << "Error: File " << file << " not found." << endl;
            return 2;
        }
    }
    for (const auto& file : mc_gen_dvcsgen_filenames) {
        if (!checkFileExists(file)) {
            cerr << "Error: File " << file << " not found." << endl;
            return 2;
        }
    }
    for (const auto& file : mc_rec_dvcsgen_filenames) {
        if (!checkFileExists(file)) {
            cerr << "Error: File " << file << " not found." << endl;
            return 2;
        }
    }
    for (const auto& file : mc_gen_aaogen_filenames) {
        if (!checkFileExists(file)) {
            cerr << "Error: File " << file << " not found." << endl;
            return 2;
        }
    }
    for (const auto& file : mc_rec_aaogen_filenames) {
        if (!checkFileExists(file)) {
            cerr << "Error: File " << file << " not found." << endl;
            return 2;
        }
    }

    cout << "All required files found. Proceeding to load data and MC files." << endl;

    // Create individual TTreeReader objects
    TFile* data_files[3];
    TFile* eppi0_files[3];
    TFile* mc_gen_dvcsgen_files[3];
    TFile* mc_rec_dvcsgen_files[3];
    TFile* mc_gen_aaogen_files[3];
    TFile* mc_rec_aaogen_files[3];

    TTreeReader data_readers[3];
    TTreeReader eppi0_readers[3];
    TTreeReader mc_gen_dvcsgen_readers[3];
    TTreeReader mc_rec_dvcsgen_readers[3];
    TTreeReader mc_gen_aaogen_readers[3];
    TTreeReader mc_rec_aaogen_readers[3];

    // Load all data and MC files
    for (size_t i = 0; i < 3; ++i) {
        // Data DVCS
        data_files[i] = new TFile(data_filenames[i].c_str(), "READ");
        TTree* data_tree = (TTree*)data_files[i]->Get("PhysicsEvents");
        data_readers[i].SetTree(data_tree);

        // Data eppi0
        eppi0_files[i] = new TFile(eppi0_filenames[i].c_str(), "READ");
        TTree* eppi0_tree = (TTree*)eppi0_files[i]->Get("PhysicsEvents");
        eppi0_readers[i].SetTree(eppi0_tree);

        // MC gen dvcsgen
        mc_gen_dvcsgen_files[i] = new TFile(mc_gen_dvcsgen_filenames[i].c_str(), "READ");
        TTree* mc_gen_dvcsgen_tree = (TTree*)mc_gen_dvcsgen_files[i]->Get("PhysicsEvents");
        mc_gen_dvcsgen_readers[i].SetTree(mc_gen_dvcsgen_tree);

        // MC rec dvcsgen
        mc_rec_dvcsgen_files[i] = new TFile(mc_rec_dvcsgen_filenames[i].c_str(), "READ");
        TTree* mc_rec_dvcsgen_tree = (TTree*)mc_rec_dvcsgen_files[i]->Get("PhysicsEvents");
        mc_rec_dvcsgen_readers[i].SetTree(mc_rec_dvcsgen_tree);

        // MC gen aaogen
        mc_gen_aaogen_files[i] = new TFile(mc_gen_aaogen_filenames[i].c_str(), "READ");
        TTree* mc_gen_aaogen_tree = (TTree*)mc_gen_aaogen_files[i]->Get("PhysicsEvents");
        mc_gen_aaogen_readers[i].SetTree(mc_gen_aaogen_tree);

        // MC rec aaogen
        mc_rec_aaogen_files[i] = new TFile(mc_rec_aaogen_filenames[i].c_str(), "READ");
        TTree* mc_rec_aaogen_tree = (TTree*)mc_rec_aaogen_files[i]->Get("PhysicsEvents");
        mc_rec_aaogen_readers[i].SetTree(mc_rec_aaogen_tree);
    }

    // Create necessary directories before proceeding
    std::string base_output_dir = "output";  // Define the base directory
    create_directories(base_output_dir);     // Create the directories

    cout << "Successfully loaded all data and MC trees." << endl << endl;



    // Call determine_exclusivity and plot variables for DVCS
    determine_exclusivity("dvcs", "(CD,FT)", data_readers[0], mc_rec_dvcsgen_readers[0], "output/exclusivity_plots", "e'p'#gamma (Fa18 Inb; CD,FT)";
    determine_exclusivity("dvcs", "(CD,FT)", data_readers[1], mc_rec_dvcsgen_readers[1], "output/exclusivity_plots", "e'p'#gamma (Fa18 Out; CD,FT)");
    determine_exclusivity("dvcs", "(CD,FT)", data_readers[2], mc_rec_dvcsgen_readers[2], "output/exclusivity_plots", "e'p'#gamma (Sp19 Inb; CD,FT)");

    determine_exclusivity("eppi0", "(CD,FT)", eppi0_readers[0], mc_rec_aaogen_readers[0], "output/exclusivity_plots", "e'p'#pi^{0} (Fa18 Inb; CD,FT)");
    determine_exclusivity("eppi0", "(CD,FT)", eppi0_readers[1], mc_rec_aaogen_readers[1], "output/exclusivity_plots", "e'p'#pi^{0} (Fa18 Out; CD,FT)");
    determine_exclusivity("eppi0", "(CD,FT)", eppi0_readers[2], mc_rec_aaogen_readers[2], "output/exclusivity_plots", "e'p'#pi^{0} (Sp19 Inb; CD,FT)");

    // End program
    cout << "Program complete. Additional functionality to be added later." << endl << endl;

    return 0;
}