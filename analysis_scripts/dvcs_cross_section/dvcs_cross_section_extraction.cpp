// Standard C++ Library Headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>  // C++17 for file existence check
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

// ROOT Libraries
#include <TROOT.h>
#include <TApplication.h>

// Namespace declaration
using namespace std;
namespace fs = std::filesystem;  // Shortened for file system ops

// Global TTreeReader variables
TTreeReader dataReader;
TTreeReader mcReader;

bool checkFileExists(const string& path) {
    return fs::exists(path);
}

int main(int argc, char* argv[]) {
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

    // Define expected files in each directory
    vector<string> dir1_files = {"rga_fa18_inb_epgamma.root", "rga_fa18_out_epgamma.root", "rga_sp19_inb_epgamma.root"};
    vector<string> dir2_files = {
        "gen_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root", "rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root",
        "gen_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root", "rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root",
        "gen_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root", "rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"};
    vector<string> dir3_files = {"rga_fa18_inb_eppi0.root", "rga_fa18_out_eppi0.root", "rga_sp19_inb_eppi0.root"};
    vector<string> dir4_files = {
        "gen_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root", "rec_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root",
        "gen_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root", "rec_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root",
        "gen_aaogen_norad_rga_sp19_inb_50nA_10604MeV_eppi0.root", "rec_aaogen_norad_rga_sp19_inb_50nA_10604MeV_eppi0.root"};

    // Check if all expected files exist in each directory
    for (const auto& file : dir1_files) {
        if (!checkFileExists(dir1 + "/" + file)) {
            cerr << "Error: File " << file << " not found in " << dir1 << endl;
            return 2;
        }
    }
    for (const auto& file : dir2_files) {
        if (!checkFileExists(dir2 + "/" + file)) {
            cerr << "Error: File " << file << " not found in " << dir2 << endl;
            return 2;
        }
    }
    for (const auto& file : dir3_files) {
        if (!checkFileExists(dir3 + "/" + file)) {
            cerr << "Error: File " << file << " not found in " << dir3 << endl;
            return 2;
        }
    }
    for (const auto& file : dir4_files) {
        if (!checkFileExists(dir4 + "/" + file)) {
            cerr << "Error: File " << file << " not found in " << dir4 << endl;
            return 2;
        }
    }

    cout << "All required files found. Proceeding to load data and MC files." << endl;

    // Example: Loading one file (extend as needed for each run period)
    TFile* data_file = new TFile((dir1 + "/rga_fa18_inb_epgamma.root").c_str(), "READ");
    TFile* mc_file = new TFile((dir2 + "/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root").c_str(), "READ");

    if (!data_file->IsOpen() || !mc_file->IsOpen()) {
        cerr << "Error: Unable to open data or MC file." << endl;
        return 3;
    }

    // Load trees
    TTree* data_tree = (TTree*)data_file->Get("PhysicsEvents");
    TTree* mc_tree = (TTree*)mc_file->Get("PhysicsEvents");

    if (!data_tree || !mc_tree) {
        cerr << "Error: PhysicsEvents tree not found in one of the files." << endl;
        return 4;
    }

    // Set up TTreeReaders
    dataReader.SetTree(data_tree);
    mcReader.SetTree(mc_tree);

    cout << "Successfully loaded data and MC trees." << endl;

    // End program for now
    cout << "Program complete. Additional functionality to be added later." << endl;

    return 0;
}