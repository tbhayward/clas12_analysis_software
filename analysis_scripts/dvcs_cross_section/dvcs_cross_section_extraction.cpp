#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <filesystem>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>

// ROOT Libraries
#include <TROOT.h>
#include <TApplication.h>

// tbhayward imports
#include "create_directories.h"
#include "plot_pi0_mass.h"
#include "determine_exclusivity.h"
#include "plot_dvcs_data_mc_comparison.h"
#include "bin_boundaries.h" 
#include "all_bin_data.h"
#include "plot_yield_comparison.h"
#include "plot_unfolding.h"
#include "calculate_contamination.h"
#include "unfolding_data.h"
#include "write_csv.h"
#include "plot_comparison.h"
#include "plot_cross_section_comparison.h"
#include "plot_cross_section_run_period_comparison.h"
#include "bin_volume.h"

// Namespace declaration
using namespace std;
namespace fs = std::filesystem;

// Function to check if a file exists
bool checkFileExists(const string& path) {
    return fs::exists(path);
}

void call_determine_exclusivity(std::vector<TTreeReader>& data_readers, std::vector<TTreeReader>& mc_rec_dvcsgen_readers, std::vector<TTreeReader>& eppi0_readers, std::vector<TTreeReader>& mc_rec_aaogen_readers) {
    // DVCS channel calls
    data_readers[0].Restart(); mc_rec_dvcsgen_readers[0].Restart();
    determine_exclusivity("dvcs", "(FD,FD)", data_readers[0], mc_rec_dvcsgen_readers[0], "output/exclusivity_plots", "e'p'#gamma (Fa18 Inb, FD,FD)");
    data_readers[0].Restart(); mc_rec_dvcsgen_readers[0].Restart();
    determine_exclusivity("dvcs", "(CD,FD)", data_readers[0], mc_rec_dvcsgen_readers[0], "output/exclusivity_plots", "e'p'#gamma (Fa18 Inb, CD,FD)");
    data_readers[0].Restart(); mc_rec_dvcsgen_readers[0].Restart();
    determine_exclusivity("dvcs", "(CD,FT)", data_readers[0], mc_rec_dvcsgen_readers[0], "output/exclusivity_plots", "e'p'#gamma (Fa18 Inb, CD,FT)");
    data_readers[0].Restart(); mc_rec_dvcsgen_readers[0].Restart();

    data_readers[1].Restart(); mc_rec_dvcsgen_readers[1].Restart(); 
    determine_exclusivity("dvcs", "(FD,FD)", data_readers[1], mc_rec_dvcsgen_readers[1], "output/exclusivity_plots", "e'p'#gamma (Fa18 Out, FD,FD)");
    data_readers[1].Restart(); mc_rec_dvcsgen_readers[1].Restart();
    determine_exclusivity("dvcs", "(CD,FD)", data_readers[1], mc_rec_dvcsgen_readers[1], "output/exclusivity_plots", "e'p'#gamma (Fa18 Out, CD,FD)");
    data_readers[1].Restart(); mc_rec_dvcsgen_readers[1].Restart();
    determine_exclusivity("dvcs", "(CD,FT)", data_readers[1], mc_rec_dvcsgen_readers[1], "output/exclusivity_plots", "e'p'#gamma (Fa18 Out, CD,FT)");
    data_readers[1].Restart(); mc_rec_dvcsgen_readers[1].Restart();

    data_readers[2].Restart(); mc_rec_dvcsgen_readers[2].Restart();
    determine_exclusivity("dvcs", "(FD,FD)", data_readers[2], mc_rec_dvcsgen_readers[2], "output/exclusivity_plots", "e'p'#gamma (Sp19 Inb, FD,FD)");
    data_readers[2].Restart(); mc_rec_dvcsgen_readers[2].Restart();
    determine_exclusivity("dvcs", "(CD,FD)", data_readers[2], mc_rec_dvcsgen_readers[2], "output/exclusivity_plots", "e'p'#gamma (Sp19 Inb, CD,FD)");
    data_readers[2].Restart(); mc_rec_dvcsgen_readers[2].Restart();   
    determine_exclusivity("dvcs", "(CD,FT)", data_readers[2], mc_rec_dvcsgen_readers[2], "output/exclusivity_plots", "e'p'#gamma (Sp19 Inb, CD,FT)");
    data_readers[2].Restart(); mc_rec_dvcsgen_readers[2].Restart();

    // Eppi0 channel calls
    eppi0_readers[0].Restart(); mc_rec_aaogen_readers[0].Restart();
    determine_exclusivity("eppi0", "(CD,FD)", eppi0_readers[0], mc_rec_aaogen_readers[0], "output/exclusivity_plots", "e'p'#pi^{0} (Fa18 Inb, CD,FD)");
    eppi0_readers[0].Restart(); mc_rec_aaogen_readers[0].Restart();
    determine_exclusivity("eppi0", "(FD,FD)", eppi0_readers[0], mc_rec_aaogen_readers[0], "output/exclusivity_plots", "e'p'#pi^{0} (Fa18 Inb, FD,FD)");
    eppi0_readers[0].Restart(); mc_rec_aaogen_readers[0].Restart();

    eppi0_readers[1].Restart(); mc_rec_aaogen_readers[1].Restart();   
    determine_exclusivity("eppi0", "(CD,FD)", eppi0_readers[1], mc_rec_aaogen_readers[1], "output/exclusivity_plots", "e'p'#pi^{0} (Fa18 Out, CD,FD)");
    eppi0_readers[1].Restart(); mc_rec_aaogen_readers[1].Restart();
    determine_exclusivity("eppi0", "(FD,FD)", eppi0_readers[1], mc_rec_aaogen_readers[1], "output/exclusivity_plots", "e'p'#pi^{0} (Fa18 Out, FD,FD)");
    eppi0_readers[1].Restart(); mc_rec_aaogen_readers[1].Restart();

    eppi0_readers[2].Restart(); mc_rec_aaogen_readers[2].Restart(); 
    determine_exclusivity("eppi0", "(CD,FD)", eppi0_readers[2], mc_rec_aaogen_readers[2], "output/exclusivity_plots", "e'p'#pi^{0} (Sp19 Inb, CD,FD)");
    eppi0_readers[2].Restart(); mc_rec_aaogen_readers[2].Restart();
    determine_exclusivity("eppi0", "(FD,FD)", eppi0_readers[2], mc_rec_aaogen_readers[2], "output/exclusivity_plots", "e'p'#pi^{0} (Sp19 Inb, FD,FD)");
    eppi0_readers[2].Restart(); mc_rec_aaogen_readers[2].Restart();
}

// Function to count the number of unique xB bins
int count_unique_xB_bins(const std::vector<BinBoundary>& bin_boundaries) {
    std::set<std::pair<double, double>> unique_xB_bins;  // Use a set to store unique (xB_low, xB_high) pairs
    
    for (const auto& bin : bin_boundaries) {
        unique_xB_bins.emplace(bin.xB_low, bin.xB_high);  // Insert unique xB bin ranges
    }
    
    return unique_xB_bins.size();  // The number of unique xB bins
}

// Helper function to print bin boundaries for verification
void print_bin_boundaries(const std::vector<BinBoundary>& bin_boundaries) {
    std::cout << "Bin boundaries read from CSV:" << std::endl;
    for (const auto& bin : bin_boundaries) {
        std::cout << "Unnamed: " << bin.unnamed << ", "
                  << "Bin Name: " << bin.bin_name << ", "
                  << "xB: [" << bin.xB_low << ", " << bin.xB_high << "], "
                  << "xB_avg: " << bin.xB_avg << ", "
                  << "Q2: [" << bin.Q2_low << ", " << bin.Q2_high << "], "
                  << "Q2_avg: " << bin.Q2_avg << ", "
                  << "t: [" << bin.t_low << ", " << bin.t_high << "], "
                  << "t_avg: " << bin.t_avg << ", "
                  << "phi: [" << bin.phi_low << ", " << bin.phi_high << "], "
                  << "phi_avg: " << bin.phi_avg << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::cout << std::endl << std::endl << std::endl;
    TApplication theApp("App", nullptr, nullptr);
    gROOT->SetBatch(kTRUE);  // Set ROOT to batch mode

    if (argc < 5) {
        std::cout << "Usage: " << argv[0] << " <dir1> <dir2> <dir3> <dir4>" << std::endl;
        return 1;
    }

    std::string dir1 = argv[1];
    std::string dir2 = argv[2];
    std::string dir3 = argv[3];
    std::string dir4 = argv[4];

    // File paths for the binning CSVs
    // this is Sangbaek's pass-1 file
    std::string binning_file = "imports/all_bin_v3.csv";  // Path to new bin boundaries file
    // Read bin boundaries from CSV
    std::vector<BinBoundary> bin_boundaries = read_bin_boundaries(binning_file);
    if (bin_boundaries.empty()) {
        std::cerr << "Error: No bin boundaries read from file." << std::endl;
        return 1;
    }
    // Print out the bin boundaries for verification
    // print_bin_boundaries(bin_boundaries);
    
    // Calculate the number of unique xB bins
    int num_xB_bins = count_unique_xB_bins(bin_boundaries);

    // Define filenames for each directory (3 periods)
    std::vector<std::string> data_filenames = {
        dir1 + "/rga_fa18_inb_epgamma.root",
        dir1 + "/rga_fa18_out_epgamma.root",
        dir1 + "/rga_sp19_inb_epgamma.root"
    };

    std::vector<std::string> eppi0_filenames = {
        dir3 + "/rga_fa18_inb_eppi0.root",
        dir3 + "/rga_fa18_out_eppi0.root",
        dir3 + "/rga_sp19_inb_eppi0.root"
    };

    std::vector<std::string> mc_gen_dvcsgen_filenames = {
        dir2 + "/gen_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root",
        dir2 + "/gen_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root",
        dir2 + "/gen_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"
    };

    std::vector<std::string> mc_rec_dvcsgen_filenames = {
        dir2 + "/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root",
        dir2 + "/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root",
        dir2 + "/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"
    };

    std::vector<std::string> mc_gen_aaogen_filenames = {
        dir4 + "/gen_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root",
        dir4 + "/gen_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root",
        dir4 + "/gen_aaogen_norad_rga_sp19_inb_50nA_10200MeV_eppi0.root"
    };

    std::vector<std::string> mc_rec_aaogen_filenames = {
        dir4 + "/rec_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root",
        dir4 + "/rec_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root",
        dir4 + "/rec_aaogen_norad_rga_sp19_inb_50nA_10200MeV_eppi0.root"
    };

    std::vector<std::string> mc_rec_eppi0_bkg_filenames = {
        dir4 + "/eppi0_bkg_aaogen_norad_rga_fa18_inb_epgamma.root",
        dir4 + "/eppi0_bkg_aaogen_norad_rga_fa18_out_epgamma.root",
        dir4 + "/eppi0_bkg_aaogen_norad_rga_sp19_inb_epgamma.root"
    };

    // Check if all expected files exist
    for (const auto& file : data_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }
    for (const auto& file : eppi0_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }
    for (const auto& file : mc_gen_dvcsgen_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }
    for (const auto& file : mc_rec_dvcsgen_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }
    for (const auto& file : mc_gen_aaogen_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }
    for (const auto& file : mc_rec_aaogen_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }
    // **Check for the new eppi0 background files**
    for (const auto& file : mc_rec_eppi0_bkg_filenames) {
        if (!checkFileExists(file)) {
            std::cerr << "Error: File " << file << " not found." << std::endl;
            return 2;
        }
    }

    std::cout << "All required files found. Proceeding to load data and MC files." << std::endl;

    // **Create individual TTreeReader objects, including the new ones**
    TFile* data_files[3];
    TFile* eppi0_files[3];
    TFile* mc_gen_dvcsgen_files[3];
    TFile* mc_rec_dvcsgen_files[3];
    TFile* mc_gen_aaogen_files[3];
    TFile* mc_rec_aaogen_files[3];
    // **New TFile pointers for the misidentified eppi0 background files**
    TFile* mc_rec_eppi0_bkg_files[3];

    std::vector<TTreeReader> data_readers(3);
    std::vector<TTreeReader> eppi0_readers(3);
    std::vector<TTreeReader> mc_gen_dvcsgen_readers(3);
    std::vector<TTreeReader> mc_rec_dvcsgen_readers(3);
    std::vector<TTreeReader> mc_gen_aaogen_readers(3);
    std::vector<TTreeReader> mc_rec_aaogen_readers(3);
    // **New TTreeReader objects for the misidentified eppi0 background files**
    std::vector<TTreeReader> mc_rec_eppi0_bkg_readers(3);

    // **Load all data and MC files, including the new ones**
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

        // **MC rec eppi0 background**
        mc_rec_eppi0_bkg_files[i] = new TFile(mc_rec_eppi0_bkg_filenames[i].c_str(), "READ");
        TTree* mc_rec_eppi0_bkg_tree = (TTree*)mc_rec_eppi0_bkg_files[i]->Get("PhysicsEvents");
        mc_rec_eppi0_bkg_readers[i].SetTree(mc_rec_eppi0_bkg_tree);
    }

    // Create necessary directories before proceeding
    std::string base_output_dir = "output";  // Define the base directory
    create_directories(base_output_dir);     // Create the directories

    std::cout << "Successfully loaded all data and MC trees and created output directories." << std::endl << std::endl;
    std::string output_dir = base_output_dir + "/data_mc_comparison/dvcs";  // Define the output directory for plots

    // Call the plotting function for the pi0 mass (optional)
    plot_pi0_mass(eppi0_readers[0], eppi0_readers[1], eppi0_readers[2],
                  mc_rec_aaogen_readers[0], mc_rec_aaogen_readers[1], mc_rec_aaogen_readers[2], "output");

    // // Call the exclusivity plots (optional)
    // call_determine_exclusivity(data_readers, mc_rec_dvcsgen_readers, eppi0_readers, mc_rec_aaogen_readers);

    // // Loop over unique xB bins and call the plotting function for DVCS data/MC comparison
    // for (int xB_bin = 0; xB_bin < num_xB_bins; ++xB_bin) {
    //     plot_dvcs_data_mc_comparison(output_dir, "dvcs", "Fa18_Out", xB_bin, bin_boundaries, data_readers[1], mc_gen_dvcsgen_readers[1], mc_rec_dvcsgen_readers[1]);
    // }

    // for (int xB_bin = 0; xB_bin < num_xB_bins; ++xB_bin) {
    // // for (int xB_bin = 0; xB_bin < 2; ++xB_bin) {
    //     plot_dvcs_data_mc_comparison(output_dir, "eppi0", "Fa18_Inb", xB_bin, bin_boundaries, eppi0_readers[0], mc_gen_aaogen_readers[0], mc_rec_aaogen_readers[0]);
    // } //

    // // Create a map to hold all the unfolding data across bins
    // std::map<std::string, std::vector<UnfoldingData>> all_unfolding_data;

    // // Iterate over the xB bins
    // for (int xB_bin = 0; xB_bin < num_xB_bins; ++xB_bin) {
    // // for (int xB_bin = 0; xB_bin < 1; ++xB_bin) {
    //     // // Call the plot_unfolding function for each xB_bin and get the results
    //     std::map<std::string, std::vector<UnfoldingData>> bin_data = plot_unfolding(base_output_dir,
    //         xB_bin, bin_boundaries, data_readers, mc_gen_dvcsgen_readers, mc_rec_dvcsgen_readers, eppi0_readers,
    //         mc_gen_aaogen_readers, mc_rec_aaogen_readers);

    //     // Calculate the contamination factor and update bin_data
    //     calculate_contamination(base_output_dir, xB_bin, bin_boundaries, data_readers, eppi0_readers,
    //         mc_rec_aaogen_readers, mc_rec_eppi0_bkg_readers, bin_data);

    //     // Merge bin_data into all_unfolding_data
    //     for (const auto& [topology, data_vector] : bin_data) {
    //         all_unfolding_data[topology].insert(all_unfolding_data[topology].end(), data_vector.begin(), data_vector.end());
    //     }
    // }

    // // Debug information to check the contents of all_unfolding_data
    // std::cout << "Debug: Number of unfolding_data entries: " << all_unfolding_data.size() << std::endl;

    // // After collecting all the data, write it to a CSV file
    // write_csv("output/unfolding_data.csv", all_unfolding_data);

    // // Call plot_comparison
    // // plot_comparison(binning_file, "output/unfolding_data.csv");
    // // plot_cross_section_comparison(binning_file, "output/unfolding_data.csv");
    // // plot_cross_section_run_period_comparison("output/unfolding_data.csv");

    std::cout << "Program complete. Additional functionality to be added later." << std::endl << std::endl;

    return 0;
}