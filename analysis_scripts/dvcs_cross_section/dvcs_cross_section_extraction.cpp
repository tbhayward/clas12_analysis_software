#include <iostream>
#include <memory>
#include <vector>
#include <set>
#include <filesystem>
#include <stdexcept>
#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TApplication.h>
#include <TROOT.h>

// Custom headers
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

using namespace std;
namespace fs = filesystem;

// Structure to manage ROOT file resources
struct RootResource {
    unique_ptr<TFile> file;
    TTreeReader reader;
    
    RootResource(const string& file_path, const string& tree_name) {
        file = make_unique<TFile>(file_path.c_str(), "READ");
        if (!file->IsOpen()) {
            throw runtime_error("Failed to open file: " + file_path);
        }
        
        TTree* tree = dynamic_cast<TTree*>(file->Get(tree_name.c_str()));
        if (!tree) {
            throw runtime_error("Missing tree '" + tree_name + "' in: " + file_path);
        }
        
        reader.SetTree(tree);
    }

    // Allow non-const access to reader for Restart()
    TTreeReader& getReader() { return reader; }
};

// Helper function to validate file existence
void validate_file(const string& path) {
    if (!fs::exists(path)) {
        throw runtime_error("Missing required file: " + path);
    }
}

// Exclusivity analysis configuration structure
struct ExclusivityConfig {
    string channel;
    string detector_config;
    string legend;
};

void call_determine_exclusivity(vector<RootResource>& data_resources,
                               vector<RootResource>& mc_dvcs_resources,
                               vector<RootResource>& eppi0_resources,
                               vector<RootResource>& mc_aao_resources) {
    // Configuration for DVCS analysis
    const vector<ExclusivityConfig> dvcs_configs = {
        {"dvcs", "(FD,FD)", "e'p'#gamma (%s, FD,FD)"},
        {"dvcs", "(CD,FD)", "e'p'#gamma (%s, CD,FD)"},
        {"dvcs", "(CD,FT)", "e'p'#gamma (%s, CD,FT)"}
    };

    // Configuration for eppi0 analysis
    const vector<ExclusivityConfig> eppi0_configs = {
        {"eppi0", "(CD,FD)", "e'p'#pi^{0} (%s, CD,FD)"},
        {"eppi0", "(FD,FD)", "e'p'#pi^{0} (%s, FD,FD)"}
    };

    const vector<string> periods = {"Fa18 Inb", "Fa18 Out", "Sp19 Inb"};
    const string output_dir = "output/exclusivity_plots";

    // Process DVCS configurations
    for (size_t i = 0; i < 3; ++i) {
        for (const auto& cfg : dvcs_configs) {
            data_resources[i].getReader().Restart();
            mc_dvcs_resources[i].getReader().Restart();
            
            determine_exclusivity(
                cfg.channel,
                cfg.detector_config,
                data_resources[i].getReader(),
                mc_dvcs_resources[i].getReader(),
                output_dir,
                TString::Format(cfg.legend.c_str(), periods[i].c_str())
            );
        }
    }

    // Process eppi0 configurations
    for (size_t i = 0; i < 3; ++i) {
        for (const auto& cfg : eppi0_configs) {
            eppi0_resources[i].getReader().Restart();
            mc_aao_resources[i].getReader().Restart();
            
            determine_exclusivity(
                cfg.channel,
                cfg.detector_config,
                eppi0_resources[i].getReader(),
                mc_aao_resources[i].getReader(),
                output_dir,
                TString::Format(cfg.legend.c_str(), periods[i].c_str())
            );
        }
    }
}

int main(int argc, char* argv[]) {
    cout << "\n=== CLAS12 DVCS Cross Section Analysis ===\n\n";
    TApplication theApp("App", nullptr, nullptr);
    gROOT->SetBatch(kTRUE);

    try {
        // ---------------------------
        // 1. Command Line Processing
        // ---------------------------
        if (argc < 5) {
            cerr << "Usage: " << argv[0] 
                 << " <data_dir> <dvcs_mc_dir> <eppi0_dir> <aao_mc_dir>\n"
                 << "Missing required directory arguments\n";
            return 1;
        }

        const vector<string> dirs = {argv[1], argv[2], argv[3], argv[4]};
        const string binning_file = "imports/all_bin_v3.csv";

        // ---------------------------
        // 2. File Path Configuration
        // ---------------------------
        const vector<string> periods = {"fa18_inb", "fa18_out", "sp19_inb"};
        const vector<string> file_types = {
            "/rga_%s_epgamma.root",          // Data DVCS
            "/gen_dvcsgen_rga_%s_epgamma.root",  // MC Gen DVCS
            "/rec_dvcsgen_rga_%s_epgamma.root",  // MC Rec DVCS
            "/rga_%s_eppi0.root",            // Data eppi0
            "/gen_aaogen_norad_rga_%s_eppi0.root",  // MC Gen AAO
            "/rec_aaogen_norad_rga_%s_eppi0.root",  // MC Rec AAO
            "/eppi0_bkg_aaogen_norad_rga_%s_epgamma.root"  // MC eppi0 bkg
        };

        // ---------------------------
        // 3. File Loading & Validation
        // ---------------------------
        vector<RootResource> data_resources, mc_dvcs_gen, mc_dvcs_rec,
                            eppi0_data, mc_aao_gen, mc_aao_rec, mc_eppi0_bkg;

        for (const auto& period : periods) {
            // Validate and load resources
            auto load_resource = [&](int dir_idx, const string& pattern) {
                string path = dirs[dir_idx] + TString::Format(pattern.c_str(), period.c_str()).Data();
                validate_file(path);
                return RootResource(path, "PhysicsEvents");
            };

            data_resources.push_back(load_resource(0, file_types[0]));
            mc_dvcs_gen.push_back(load_resource(1, file_types[1]));
            mc_dvcs_rec.push_back(load_resource(1, file_types[2]));
            eppi0_data.push_back(load_resource(2, file_types[3]));
            mc_aao_gen.push_back(load_resource(3, file_types[4]));
            mc_aao_rec.push_back(load_resource(3, file_types[5]));
            mc_eppi0_bkg.push_back(load_resource(3, file_types[6]));
        }

        // ---------------------------
        // 4. Analysis Initialization
        // ---------------------------
        create_directories("output");
        auto bin_boundaries = read_bin_boundaries(binning_file);
        if (bin_boundaries.empty()) {
            throw runtime_error("Failed to load bin boundaries from: " + binning_file);
        }

        // Calculate unique xB bins using set
        set<pair<double, double>> unique_xb;
        for (const auto& bin : bin_boundaries) {
            unique_xb.emplace(bin.xB_low, bin.xB_high);
        }
        const int num_xb_bins = unique_xb.size();

        // ---------------------------
        // 5. Core Analysis Pipeline
        // ---------------------------
        // 5.1 Pi0 Mass Validation
        plot_pi0_mass(eppi0_data[0].getReader(), eppi0_data[1].getReader(), eppi0_data[2].getReader(),
                     mc_aao_rec[0].getReader(), mc_aao_rec[1].getReader(), mc_aao_rec[2].getReader(),
                     "output");

        // 5.2 Exclusivity Cuts
        call_determine_exclusivity(data_resources, mc_dvcs_rec, 
                                  eppi0_data, mc_aao_rec);

        // 5.3 Data/MC Comparisons
        const string comp_dir = "output/data_mc_comparison";
        for (int xb = 0; xb < num_xb_bins; ++xb) {
            // DVCS comparison
            plot_dvcs_data_mc_comparison(
                comp_dir + "/dvcs", "dvcs", "Fa18_Out", xb, bin_boundaries,
                data_resources[1].getReader(), mc_dvcs_gen[1].getReader(), mc_dvcs_rec[1].getReader()
            );

            // eppi0 comparison
            plot_dvcs_data_mc_comparison(
                comp_dir + "/eppi0", "eppi0", "Fa18_Inb", xb, bin_boundaries,
                eppi0_data[0].getReader(), mc_aao_gen[0].getReader(), mc_aao_rec[0].getReader()
            );
        }

        // ---------------------------
        // 6. Unfolding & Contamination
        // ---------------------------
        map<string, vector<UnfoldingData>> unfolding_results;

        for (int xb = 0; xb < num_xb_bins; ++xb) {
            vector<TTreeReader> data_readers = {
                data_resources[0].getReader(),
                data_resources[1].getReader(),
                data_resources[2].getReader()
            };

            vector<TTreeReader> mc_gen_readers = {
                mc_dvcs_gen[0].getReader(),
                mc_dvcs_gen[1].getReader(),
                mc_dvcs_gen[2].getReader()
            };

            vector<TTreeReader> mc_rec_readers = {
                mc_dvcs_rec[0].getReader(),
                mc_dvcs_rec[1].getReader(),
                mc_dvcs_rec[2].getReader()
            };

            vector<TTreeReader> eppi0_readers = {
                eppi0_data[0].getReader(),
                eppi0_data[1].getReader(),
                eppi0_data[2].getReader()
            };

            vector<TTreeReader> mc_aao_gen_readers = {
                mc_aao_gen[0].getReader(),
                mc_aao_gen[1].getReader(),
                mc_aao_gen[2].getReader()
            };

            vector<TTreeReader> mc_aao_rec_readers = {
                mc_aao_rec[0].getReader(),
                mc_aao_rec[1].getReader(),
                mc_aao_rec[2].getReader()
            };

            auto bin_data = plot_unfolding(
                "output", xb, bin_boundaries,
                data_readers, mc_gen_readers, mc_rec_readers,
                eppi0_readers, mc_aao_gen_readers, mc_aao_rec_readers
            );

            calculate_contamination(
                "output", xb, bin_boundaries,
                data_readers,
                eppi0_readers,
                mc_aao_rec_readers,
                {mc_eppi0_bkg[0].getReader(), mc_eppi0_bkg[1].getReader(), mc_eppi0_bkg[2].getReader()},
                bin_data
            );

            // Aggregate results
            for (auto& [key, vec] : bin_data) {
                unfolding_results[key].insert(unfolding_results[key].end(),
                                            vec.begin(), vec.end());
            }
        }

        // ---------------------------
        // 7. Final Output & Results
        // ---------------------------
        write_csv("output/unfolding_data.csv", unfolding_results);
        plot_comparison(binning_file, "output/unfolding_data.csv");
        plot_cross_section_comparison(binning_file, "output/unfolding_data.csv");
        plot_cross_section_run_period_comparison("output/unfolding_data.csv");

        cout << "\n=== Analysis Completed Successfully ===\n";
        return 0;

    } catch (const exception& e) {
        cerr << "\n!!! Analysis Failed: " << e.what() << " !!!\n";
        return 2;
    }
}