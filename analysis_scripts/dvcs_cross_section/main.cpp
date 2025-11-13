#include "initialize_pass2_csv.h"
#include "make_dirs.h"
#include "load_trees.h"
#include "periods.h"
#include "global_cuts.h"
#include "exclusivity_cuts.h"
#include "load_binning_scheme.h"
#include "bin_means.h"
#include "total_counts.h"
#include "pi0_contamination.h"
#include "pi0_corrected_counts.h"
#include "bsa.h"
#include "radiative_corrections.h"
#include <filesystem>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include "acceptance.h"
#include "unfolding.h"
#include "bin_volume.h"
#include "uncorrected_cross_section.h"
#include "rad_corrected_cross_section.h"
#include "model_predictions.h"
#include "bin_centering_corrections.h"
#include "models_vs_data_plots.h"

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    // Create necessary output directories
    makeOutputDirs();
    std::cout << "Output directories ready." << std::endl;

    initialize_pass2_csv("imports/all_bin_v3.csv", "output/csvs/dvcs_pass2_analysis.csv");

    // Root of output tree (used by several stages)
    const std::string output_root = "output";

    // // --- Load binning scheme ---
    // const std::string csv_file_path = "imports/integrated_bin_v2.csv";
    // auto binning_scheme = load_binning_scheme(csv_file_path);
    // std::cout << "Loaded binning scheme: " << binning_scheme.size() << " bins" << std::endl;

    // Containers for different tree categories
    std::map<std::string, TTree*> dataTrees;        // DVCS data
    std::map<std::string, TTree*> genMcTrees;       // DVCS generated MC (no-rad)
    std::map<std::string, TTree*> recMcTrees;       // DVCS reconstructed MC (no-rad)
    std::map<std::string, TTree*> eppi0DataTrees;   // eppi0 data
    std::map<std::string, TTree*> eppi0GenMcTrees;  // eppi0 generated MC
    std::map<std::string, TTree*> eppi0RecMcTrees;  // eppi0 reconstructed MC
    std::map<std::string, TTree*> eppi0BkgTrees;    // eppi0 background MC   (FIX: added)
    std::map<std::string, TTree*> radGenMcTrees;    // NEW: DVCS generated MC (radiative)
    std::map<std::string, TTree*> radRecMcTrees;    // NEW: DVCS reconstructed MC (radiative)

    // Load all trees from files
    loadTrees(dataTrees, genMcTrees, recMcTrees,
              eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees, eppi0BkgTrees,
              radGenMcTrees, radRecMcTrees);
    std::cout << "All trees loaded successfully." << std::endl;

    // Run exclusivity cut extraction 
    // Record the exact global cuts used:
    write_global_cuts_config_json("output/jsons");
    runAllExclusivityCuts(
        dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
        "output/jsons", "output/exclusivity_plots", 1
    );
    std::cout << "Exclusivity-cut stage finished." << std::endl;

    // --------- Global bin-averaged kinematics (CSV update) ----------
    {
        const std::string csv_main   = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string csv_backup = "output/csvs/dvcs_pass2_analysis_backup_bin_means.csv";

        // Make a simple backup before modifying
        try {
            std::filesystem::copy_file(csv_main, csv_backup,
                                       std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to " << csv_backup << "\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: Backup failed (" << e.what() << "). Continuing anyway.\n";
        }

        // dataTrees already built (keys like DVCS_Fa18_inb, ...). Launch with up to 5 workers.
        if (!update_bin_means_csv(csv_main, dataTrees, /*max_workers=*/5)) {
            std::cerr << "[main] ERROR: update_bin_means_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // --------- Raw yields (counts) into CSV + plots ----------
    {
        const std::string csv_main  = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";

        // Make a backup (the function also backs up to ..._total_counts.csv)
        try {
            std::filesystem::copy_file(csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_total_counts.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_total_counts.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup failed (" << e.what() << "). Continuing.\n";
        }

        // update_total_counts_csv() now discovers periods/topologies internally
        if (!update_total_counts_csv(csv_main, dataTrees, cuts_json, output_root,
                /*max_workers=*/5)) {
            std::cerr << "[main] ERROR: update_total_counts_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // --------- pi0 contamination (helicity-averaged; DVCS counts from CSV; eppi0 counts re-counted) ----------
    {
        const std::string csv_main  = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";

        // Make dirs up front (preferred place)
        makeOutputDirs();

        // Make a backup (the function also backs up to ..._total_counts.csv)
        try {
            std::filesystem::copy_file(
                csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_pi0_contamination.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_pi0_contamination.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup failed (" << e.what() << "). Continuing.\n";
        }

        // Parallelize across periods (<=0 means auto thread count)
        const int max_workers = 5;

        if (!compute_pi0_contamination_overall(
                eppi0DataTrees,
                eppi0RecMcTrees,
                eppi0BkgTrees,
                cuts_json,
                csv_main,
                output_root,
                max_workers))
        {
            std::cerr << "[main] ERROR: compute_pi0_contamination_overall failed.\n";
            std::exit(EXIT_FAILURE);
        }
        std::cout << "pi0 contamination stage finished.\n";
    }

    // // // --------- Total counts after exclusivity cuts (by helicity) ----------
    // const std::string cuts_json_path   = "output/jsons/combined_cuts.json"; 
    // // // produced by exclusivity_cuts

    // const std::string output_counts_js = "output/jsons/total_counts.json";
    // compute_total_counts(dvcs_periods, topologies, binning_scheme, dataTrees, cuts_json_path,
    // output_counts_js, output_root); 

    // // Helicity-resolved π0 contamination
    // // NOTE: pass the OUTPUT ROOT ("output") so the implementation writes:
    // //   - per-period JSONs to output/jsons/contamination/
    // //   - combined JSON to output/jsons/
    // //   - plots to output/contamination_plots/...
    // compute_pi0_contamination_helicity(
    //     dvcs_periods,
    //     topologies,
    //     binning_scheme,
    //     dataTrees,
    //     eppi0DataTrees,
    //     eppi0RecMcTrees,   // reco MC
    //     eppi0BkgTrees,     // bkg MC
    //     cuts_json_path,
    //     output_root
    // );

    // // --------- π0-corrected helicity counts (per φ) ----------
    // const std::string total_counts_json        = "output/jsons/total_counts.json";
    // const std::string contamination_dir_counts = "output/jsons/contamination";
    // const std::string contamination_combined   = "output/jsons/pi0_contamination_combined.json";

    // compute_pi0_corrected_counts(
    //     dvcs_periods,
    //     binning_scheme,
    //     total_counts_json,        // from compute_total_counts()
    //     contamination_dir_counts, // per-period contamination_*.json live here
    //     contamination_combined,   // combined groups (Spring2018, Fall2018, 10.6_GeV)
    //     output_root               // "output"
    // );
    // std::cout << "π0-corrected counts stage finished." << std::endl;

    // // --------- Radiative corrections (generated Born/Rad, per bin & φ) ----------
    // // Writes:
    // //   - per-group JSONs: output/jsons/radiative_corrections_group_<energy>.json
    // //   - per-period JSONs: output/jsons/radiative_corrections_<period>.json
    // //   - all-groups file: output/jsons/radiative_corrections_all_groups.json
    // //   - plots (ONLY per beam energy): output/radiative_correction_plots/{10.59,10.60,10.2}/...
    // compute_radiative_corrections(
    //     dvcs_periods,
    //     binning_scheme,
    //     genMcTrees,
    //     radGenMcTrees,
    //     output_root
    // );

    // // --------- Acceptance (reco MC with MC cuts / gen MC, per bin & φ) ----------
    // // Writes per-period JSONs: output/jsons/acceptance_<period>.json
    // // Plots to: output/acceptance/<runTag>/plot_acceptance_<period>_xB_<ix>.png
    // {
    //     std::vector<std::string> acc_periods = {
    //         "DVCS_Sp18_inb", "DVCS_Sp18_out",
    //         "DVCS_Fa18_inb", "DVCS_Fa18_out",
    //         "DVCS_Sp19_inb"
    //     }; // intentionally skipping DVCS_Fa18_inb_supp
    //     compute_and_plot_acceptance(
    //         acc_periods,
    //         topologies,
    //         binning_scheme,
    //         genMcTrees,
    //         recMcTrees,
    //         cuts_json_path,
    //         output_root
    //     );
    // }

    // // --------- Unfolding (counts / acceptance), helicity-resolved ----------
    // // Writes per-period JSONs: output/jsons/unfolded_<period>.json
    // // Plots to: output/unfolding/<runTag>/plot_unfolded_<period>_xB_<ix>.png
    // {
    //     std::vector<std::string> unf_periods = {
    //         "DVCS_Sp18_inb", "DVCS_Sp18_out",
    //         "DVCS_Fa18_inb", "DVCS_Fa18_out",
    //         "DVCS_Sp19_inb"
    //     }; // skip DVCS_Fa18_inb_supp on purpose
    //     // note that we pass the pi0_corrected_counts below (i.e. not the original total_counts)
    //     const std::string total_counts_js = "output/jsons/pi0_corrected_counts_all_groups.json";
    //     compute_and_plot_unfolding(
    //         unf_periods,
    //         binning_scheme,
    //         total_counts_js,
    //         output_root
    //     );
    // }

    // // --------- Bin Volume (generator-based φ coverage), per beam energy ----------
    // // Writes per-energy JSONs: output/jsons/bin_volume_<energy>.json
    // // Plots to: output/bin_volume/<energy>/plot_bin_volume_<energy>_xB_<ix>.png
    // compute_and_plot_bin_volume(
    //     binning_scheme,
    //     genMcTrees,
    //     output_root
    // );

    // // For polarized cross sections (existing functionality, but now with proper luminosity calculation)
    // compute_uncorrected_cross_sections(
    //     binning_scheme,
    //     "output/jsons",                    // bin volume JSON directory
    //     "output/jsons",                    // unfolded counts per helicity
    //     "imports/integrated_luminosity",   // luminosity text files
    //     "output/uncorrected_cross_section" // output dir
    // );

    // // NEW: For unpolarized cross sections using total charge
    // compute_unpolarized_cross_sections(
    //     binning_scheme,
    //     "output/jsons",                    // bin volume JSON directory  
    //     "output/jsons",                    // unfolded counts
    //     "imports/integrated_luminosity",   // luminosity text files
    //     "output/uncorrected_cross_section" // output dir
    // );

    // // Comparison function (uses the new luminosity calculation)
    // compare_unpolarized_cross_sections_sp18out_vs_fa18out(
    //     binning_scheme,
    //     "output/jsons",
    //     "output/jsons", 
    //     "imports/integrated_luminosity",
    //     "output/uncorrected_cross_section"
    // );


    // // Multiply uncorrected dσ/dφ by Born/Rad per-φ correction
    // const std::string unx_dir = "output/jsons";                    // reads: uncorrected_xsec_<E>.json
    // const std::string rc_dir  = "output/jsons";                    // reads: radiative_corrections_group_<E>.json
    // const std::string out_dir = "output/rad_corrected_cross_section"; // writes plots here
    // compute_rad_corrected_cross_sections(binning_scheme, unx_dir, rc_dir, out_dir);


    // // Beam-Spin Asymmetry:
    // // - reads total_counts.json and contamination JSONs
    // // - writes per-period fits to output/jsons/BSA_fits/BSA_fits_<period>.json
    // // - writes all-periods file to output/jsons/BSA_fits_all_periods.json
    // // - writes 10.6 GeV combined to output/jsons/BSA_fits_combined_10p6.json
    // // - plots to output/bsa_plots/<runTag>/...
    // namespace fs = std::filesystem;
    // compute_and_plot_bsa_helicity(
    //     dvcs_periods,                                               
    //     topologies,
    //     binning_scheme,
    //     dataTrees,       
    //     (fs::path(output_root) / "jsons" / "pi0_corrected_counts_all_groups.json").string(),
    //     output_root,  
    //     (fs::path(output_root) / "rad_corrected_cross_section" / "jsons").string() // directory with rad_corrected_xsec_<E>.json
    // );

    // // Example: unpolarized predictions (xB,Q2,t,phi,Ebeam)
    // ModelPaths paths; // leave empty to use env/defaults
    // double xB = 0.11; double Q2 = 1.6; double tpos = 0.20; double phi_deg = 180;
    // double xs_vgg  = vgg_xs(xB, Q2, tpos, phi_deg, 10.604, Helicity::Plus, paths, /*globalfit=*/false);
    // double xs_bh   = vgg_bh_only(xB, Q2, tpos, phi_deg, 10.604, paths, /*globalfit=*/false);
    // double xs_km15 = km15_xs(xB, Q2, tpos, phi_deg, 10.604, Helicity::Plus, paths);
    // std::cout << xs_vgg << " " << xs_km15 << " " << xs_bh << std::endl;
    // xs_vgg  = vgg_xs(xB, Q2, tpos, phi_deg, 10.594, Helicity::Plus, paths, /*globalfit=*/false);
    // xs_bh   = vgg_bh_only(xB, Q2, tpos, phi_deg, 10.594, paths, /*globalfit=*/false);
    // xs_km15 = km15_xs(xB, Q2, tpos, phi_deg, 10.594, Helicity::Plus, paths);
    // std::cout << xs_vgg << " " << xs_km15 << " " << xs_bh << std::endl;

    // // --------- Bin-centering corrections using VGG ONLY (compute only) ----------
    // compute_bin_centered_cross_sections(
    //     binning_scheme,
    //     "output/jsons",           // rad_corrected_xsec_<E>.json input
    //     "output/bin_centering",   // JSONs saved to output/bin_centering/jsons
    //     4,
    //     ModelPaths(),
    //     true,
    //     ModelChoice::Both      // VGGOnly, KM15Only, Both
    // );
    // // --------- Bin-centering plots only (no recompute) ----------
    // plot_bin_centered_cross_sections(
    //     binning_scheme,
    //     "output/jsons",                // where rad_corrected_xsec_<E>.json lives (before BC)
    //     "output/bin_centering/jsons",  // where bin_centered_xsec_<E>.json lives (after BC)
    //     "output/bin_centering/plots"   // where to write the PNGs
    // );

    // // --------- 1) Compute and save model predictions (run once, or when you change models/phi grid) ----------
    // compute_model_predictions(
    //     binning_scheme,
    //     "output/bin_centering/jsons",        // reads bin_centered_xsec_<E>.json to know which bins exist
    //     "output/model_predictions",          // will write jsons/model_predictions_<E>.json here
    //     ModelPaths(),                        // configure as needed
    //     true,                                // vgg_globalfit flag
    //     24                                  // phi sampling density (0..360 by 1 deg)
    // );

    // // --------- 2) Plot using the saved predictions (fast; tweak aesthetics without recomputing) ----------
    // plot_models_vs_bincentered(
    //     binning_scheme,
    //     "output/bin_centering/jsons",        // where bin_centered_xsec_<E>.json live (data points)
    //     "output/model_predictions/jsons",    // where model_predictions_<E>.json live (precomputed curves)
    //     "output/bin_centering",              // PNGs go under output/bin_centering/plots
    //     0                                     // unused (predictions carry their own phi grid)
    // );

    std::cout << "All done." << std::endl;
    return 0;
}