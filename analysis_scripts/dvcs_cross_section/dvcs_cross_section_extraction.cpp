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

// Global 3D vector for TTreeReaders
vector<vector<vector<TTreeReader>>> treeReaders(3, vector<vector<TTreeReader>>(2));

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

    // Define expected files for each directory (Data and MC for each run period)
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

    // 3D vector of TTreeReader* for 3 periods, 2 (data/MC), and respective files
    vector<vector<vector<TTreeReader*>>> treeReaders(3, vector<vector<TTreeReader*>>(2));

    // Load Data Trees (index 0 for data, 1 for MC)
    vector<string> data_filenames = {dir1 + "/rga_fa18_inb_epgamma.root", dir1 + "/rga_fa18_out_epgamma.root", dir1 + "/rga_sp19_inb_epgamma.root"};
    vector<string> eppi0_filenames = {dir3 + "/rga_fa18_inb_eppi0.root", dir3 + "/rga_fa18_out_eppi0.root", dir3 + "/rga_sp19_inb_eppi0.root"};

    for (size_t i = 0; i < 3; ++i) {
        // Load DVCS data
        TFile* data_file = new TFile(data_filenames[i].c_str(), "READ");
        if (!data_file->IsOpen()) {
            cerr << "Error: Unable to open DVCS data file " << data_filenames[i] << endl;
            return 3;
        }
        TTree* data_tree = (TTree*)data_file->Get("PhysicsEvents");
        if (!data_tree) {
            cerr << "Error: PhysicsEvents tree not found in " << data_filenames[i] << endl;
            return 4;
        }
        treeReaders[i][0].push_back(new TTreeReader(data_tree));  // DVCS data

        // Load eppi0 data
        TFile* eppi0_file = new TFile(eppi0_filenames[i].c_str(), "READ");
        if (!eppi0_file->IsOpen()) {
            cerr << "Error: Unable to open eppi0 data file " << eppi0_filenames[i] << endl;
            return 3;
        }
        TTree* eppi0_tree = (TTree*)eppi0_file->Get("PhysicsEvents");
        if (!eppi0_tree) {
            cerr << "Error: PhysicsEvents tree not found in " << eppi0_filenames[i] << endl;
            return 4;
        }
        treeReaders[i][0].push_back(new TTreeReader(eppi0_tree));  // eppi0 data
    }

    // Load MC Trees (index 0 for gen dvcsgen, 1 for rec dvcsgen, 2 for gen aaogen, 3 for rec aaogen)
    vector<string> mc_gen_dvcsgen_filenames = {dir2 + "/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root", dir2 + "/gen_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root", dir2 + "/gen_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"};
    vector<string> mc_rec_dvcsgen_filenames = {dir2 + "/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root", dir2 + "/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root", dir2 + "/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"};
    vector<string> mc_gen_aaogen_filenames = {dir4 + "/gen_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root", dir4 + "/gen_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root", dir4 + "/gen_aaogen_norad_rga_sp19_inb_50nA_10604MeV_eppi0.root"};
    vector<string> mc_rec_aaogen_filenames = {dir4 + "/rec_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root", dir4 + "/rec_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root", dir4 + "/rec_aaogen_norad_rga_sp19_inb_50nA_10604MeV_eppi0.root"};

    for (size_t i = 0; i < 3; ++i) {
        // Load gen dvcsgen
        TFile* gen_dvcsgen_file = new TFile(mc_gen_dvcsgen_filenames[i].c_str(), "READ");
        if (!gen_dvcsgen_file->IsOpen()) {
            cerr << "Error: Unable to open gen dvcsgen MC file " << mc_gen_dvcsgen_filenames[i] << endl;
            return 3;
        }
        TTree* gen_dvcsgen_tree = (TTree*)gen_dvcsgen_file->Get("PhysicsEvents");
        treeReaders[i][1].push_back(new TTreeReader(gen_dvcsgen_tree));

        // Load rec dvcsgen
        TFile* rec_dvcsgen_file = new TFile(mc_rec_dvcsgen_filenames[i].c_str(), "READ");
        if (!rec_dvcsgen_file->IsOpen()) {
            cerr << "Error: Unable to open rec dvcsgen MC file " << mc_rec_dvcsgen_filenames[i] << endl;
            return 3;
        }
        TTree* rec_dvcsgen_tree = (TTree*)rec_dvcsgen_file->Get("PhysicsEvents");
        treeReaders[i][1].push_back(new TTreeReader(rec_dvcsgen_tree));

        // Load gen aaogen
        TFile* gen_aaogen_file = new TFile(mc_gen_aaogen_filenames[i].c_str(), "READ");
        if (!gen_aaogen_file->IsOpen()) {
            cerr << "Error: Unable to open gen aaogen MC file " << mc_gen_aaogen_filenames[i] << endl;
            return 3;
        }
        TTree* gen_aaogen_tree = (TTree*)gen_aaogen_file->Get("PhysicsEvents");
        treeReaders[i][1].push_back(new TTreeReader(gen_aaogen_tree));

        // Load rec aaogen
        TFile* rec_aaogen_file = new TFile(mc_rec_aaogen_filenames[i].c_str(), "READ");
        if (!rec_aaogen_file->IsOpen()) {
            cerr << "Error: Unable to open rec aaogen MC file " << mc_rec_aaogen_filenames[i] << endl;
            return 3;
        }
        TTree* rec_aaogen_tree = (TTree*)rec_aaogen_file->Get("PhysicsEvents");
        treeReaders[i][1].push_back(new TTreeReader(rec_aaogen_tree));
    }

    cout << "Successfully loaded all data and MC trees." << endl;

    // Clean up dynamically allocated TTreeReaders
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 2; ++j) {
            for (TTreeReader* reader : treeReaders[i][j]) {
                delete reader;
            }
        }
    }

    // End program
    cout << "Program complete. Additional functionality to be added later." << endl;

    return 0;
}