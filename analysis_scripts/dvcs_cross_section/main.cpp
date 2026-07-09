#include "initialize_pass2_csv.h"
#include "make_dirs.h"
#include "load_trees.h"
#include "periods.h"
#include "global_cuts.h"
#include "exclusivity_cuts.h"
#include "q2_xb_cut_coverage.h"
#include "load_binning_scheme.h"
#include "bin_means.h"
#include "total_counts.h"
#include "current_dependence.h"
#include "eppi0_normalization.h"
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
#include "yield_totals.h"
#include "bin_volume.h"
#include "model_predictions.h"
#include "bin_centering_corrections.h"
#include "cross_sections.h"
#include "models_vs_data_plots.h"
#include "overall_normalization_study.h"
#include "propagator_study.h"
#include "norm_cross_sections.h"
#include "pass1_paper_plots.h"
#include "branch_data_mc_comparison.h"
#include "pi0_subtracted_kinematics.h"

int main(int argc, char* argv[]) {
    std::cout << "Starting DVCS analysis..." << std::endl;

    // Create necessary output directories
    makeOutputDirs();
    std::cout << "Output directories ready." << std::endl;

    // -------------------------------------------------------------------------
    // Global event-selection toggles for topology/sector systematics studies.
    //
    // Nominal inclusive setting: leave every enable_* flag below false.
    // The selected configuration is installed once here and then used by every
    // stage through default_global_cuts().
    //
    // Detector topology codes:
    //   FD-FD = proton FD, photon FD = detector1=1, detector2=1
    //   CD-FD = proton CD, photon FD = detector1=2, detector2=1
    //   CD-FT = proton CD, photon FT = detector1=2, detector2=0
    //
    // FD sectors:
    //   1: 330-30 deg, 2: 30-90 deg, 3: 90-150 deg,
    //   4: 150-210 deg, 5: 210-270 deg, 6: 270-330 deg
    // CD sectors:
    //   1: 272.5-32.5 deg, 2: 32.5-150.5 deg, 3: 150.5-272.5 deg
    // -------------------------------------------------------------------------
    GlobalCutConfig global_cfg;

    // Single-topology study. Enable exactly one topology by setting this true
    // and editing required_detector1/required_detector2.
    global_cfg.enable_topology_filter = false;
    global_cfg.required_detector1 = 2;  // 1 FD proton, 2 CD proton
    global_cfg.required_detector2 = 0;  // 0 FT photon, 1 FD photon

    // Electron FD sector study. Electron is always FD.
    global_cfg.enable_electron_fd_sector_filter = false;
    global_cfg.electron_fd_sector = 1;

    // Proton FD sector study. This automatically keeps only FD-FD events.
    global_cfg.enable_proton_fd_sector_filter = false;
    global_cfg.proton_fd_sector = 1;

    // Proton CD sector study. This automatically keeps only CD-FD and CD-FT events.
    global_cfg.enable_proton_cd_sector_filter = true;
    global_cfg.proton_cd_sector = 3;

    // Photon FD sector study. This automatically keeps only CD-FD and FD-FD events.
    global_cfg.enable_photon_fd_sector_filter = false;
    global_cfg.photon_fd_sector = 1;

    // Auxiliary fiducial cuts. Enable this single switch to apply the additional
    // FD-sector separation, particle-angle, and FT-photon momentum cuts
    // analysis-wide. The numerical values are defined in GlobalCutConfig.
    global_cfg.enable_auxiliary_fiducial_cuts = true;

    // Existing optional propagator/ycol mirror cut.
    global_cfg.enable_dvcsgen_ycol_cut = false;
    global_cfg.dvcsgen_ycol_cut = 0.005;

    set_default_global_cuts(global_cfg);

    std::cout << "[main] Global cut analysis tag: "
              << global_cuts_analysis_tag(default_global_cuts()) << std::endl;

    initialize_pass2_csv("imports/all_bin_v3.csv", "output/csvs/dvcs_pass2_analysis.csv");

    // Root of output tree (used by several stages)
    const std::string output_root = "output";

    // Containers for different tree categories
    std::map<std::string, TTree*> dataTrees;                  // DVCS data
    std::map<std::string, TTree*> genMcTrees;                 // DVCS generated MC (no-rad, reference current)
    std::map<std::string, TTree*> recMcTrees;                 // DVCS reconstructed MC (no-rad, reference current)
    std::map<std::string, TTree*> eppi0DataTrees;             // eppi0 data
    std::map<std::string, TTree*> eppi0GenMcTrees;            // eppi0 generated MC
    std::map<std::string, TTree*> eppi0RecMcTrees;            // eppi0 reconstructed MC
    std::map<std::string, TTree*> eppi0BkgTrees;              // eppi0 background MC
    std::map<std::string, TTree*> radGenMcTrees;              // DVCS generated MC (radiative)
    std::map<std::string, TTree*> radRecMcTrees;              // DVCS reconstructed MC (radiative)
    std::map<std::string, TTree*> currentStudyGenMcTrees;     // DVCS generated MC for current-dependence study
    std::map<std::string, TTree*> currentStudyRecMcTrees;     // DVCS reconstructed MC for current-dependence study

    // Load all trees from files
    if (!loadTrees(dataTrees, genMcTrees, recMcTrees,
        eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees,
        eppi0BkgTrees,
        radGenMcTrees, radRecMcTrees,
        currentStudyGenMcTrees, currentStudyRecMcTrees)) {
        std::cerr << "[main] FATAL: loadTrees failed.\n";
        return 1;
    }

    std::cout << "All trees loaded successfully." << std::endl;
    std::cout << "Current-study generated MC trees loaded: "
              << currentStudyGenMcTrees.size() << std::endl;
    std::cout << "Current-study reconstructed MC trees loaded: "
              << currentStudyRecMcTrees.size() << std::endl;

    // Run exclusivity cut extraction 
    // Record the exact global cuts used:
    write_global_cuts_config_json("output/jsons");
    runAllExclusivityCuts(
        dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
        "output/jsons", "output/exclusivity_plots", 5
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

    // --------- Q2 vs xB coverage after global + data 3sigma DVCS cuts ----------
    {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";
        const std::string out_dir = "output/q2_xb_cut_coverage";

        if (!plot_q2_xb_cut_coverage(dataTrees, csv_main, cuts_json, out_dir)) {
            std::cerr << "[main] ERROR: plot_q2_xb_cut_coverage failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // Optional DVCS MC acceptance override:
    //   false -> use the standard production-current/background-overlaid DVCS MC
    //            for ep->epg generated/reconstructed counts and apply the fitted
    //            ep->epg MC current correction.
    //   true  -> use the no-background dvcsgen files for ep->epg MC counts and
    //            write ep->epg MC current-efficiency factors as unity.
    //
    // This affects only ep->epg MC. Data, ep->eppi0 MC, and ep->eppi0->epg
    // background MC are kept in their standard treatment.
    const bool use_nobkg_dvcs_mc_for_acceptance = false;

    // Misidentified pi0-background MC current-correction convention:
    //   true  -> default: correct ep->eppi0->epg MC with the ep->epg MC
    //            current-efficiency factor, because the reconstructed final
    //            state being counted is epgamma.
    //   false -> legacy behavior: correct ep->eppi0->epg MC with the ep->eppi0
    //            MC current-efficiency factor.
    const bool use_epg_mc_current_factor_for_eppi0_bkg = true;

    // --------- Raw yields/counts into CSV + plots ----------
    {
        const std::string csv_main  = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";

        // Make a backup
        try {
            std::filesystem::copy_file(csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_total_counts.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_total_counts.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup failed (" << e.what() << "). Continuing.\n";
        }

        // update_total_counts_csv() fills:
        //   - DVCS data raw yields
        //   - eppi0 data raw yields
        //   - DVCS generated/reconstructed MC yields
        //   - eppi0 generated/reconstructed MC yields
        //   - eppi0-background-as-DVCS reconstructed MC yields
        TotalCountsOptions total_count_opts;
        total_count_opts.use_nobkg_dvcs_mc_counts = use_nobkg_dvcs_mc_for_acceptance;

        if (!update_total_counts_csv(csv_main, dataTrees, eppi0DataTrees,
            genMcTrees, recMcTrees,
            eppi0GenMcTrees, eppi0RecMcTrees,
            eppi0BkgTrees,
            cuts_json,
            output_root,
            /*max_workers=*/5,
            total_count_opts,
            currentStudyGenMcTrees,
            currentStudyRecMcTrees)) {
          std::cerr << "[main] ERROR: update_total_counts_csv failed.\n";
          std::exit(EXIT_FAILURE);
        }
    }

    // --------- Current-dependence correction factors ----------
    {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        CurrentDependenceOptions current_opts;
        current_opts.charge_csv_path = "imports/integrated_luminosity/global.csv";
        current_opts.combined_cuts_json = "output/jsons/combined_cuts.json";
        current_opts.output_dir = "output/dvcs_current_dependence";

        // Existing override:
        //   false -> compute current-efficiency factors normally
        //   true  -> write all current-efficiency factors as (1,0)
        current_opts.override_to_unity = false;

        // Fallback mode if the Fa18/Sp19 columns-3-to-5 mode below is disabled.
        current_opts.use_second_column_charge_for_all_unpolarized = true;

        // New mode:
        //   Sp18 Inb / Sp18 Out:
        //     charge_unpol = column 2
        //
        //   Fa18 Inb / Fa18 Out / Sp19 Inb:
        //     charge_unpol = 1.025 * (column 3 + column 4 + column 5)
        current_opts.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized = true;
        current_opts.columns_3_to_5_charge_sum_scale = 1.025;

        // Default behavior:
        //   true  -> write Sp19 Inb current-efficiency factors using Fa18 Inb
        //            because the Sp19 5 nA luminosity-scan point has suspect
        //            Faraday Cup charge.
        //   false -> use the directly fitted Sp19 Inb current-dependence factor.
        current_opts.use_fa18_inb_current_efficiency_for_sp19_inb = true;

        // Match the total-counts override above. If true, ep->epg MC current
        // factors written to the CSV are forced to (1,0). The scan is still
        // processed diagnostically, and the pi0 treatment is unchanged.
        current_opts.use_nobkg_dvcs_mc_counts = use_nobkg_dvcs_mc_for_acceptance;
        current_opts.use_epg_mc_current_factor_for_eppi0_bkg = use_epg_mc_current_factor_for_eppi0_bkg;

        current_opts.max_workers = 5;

        if (!update_current_dependence_factors_csv(csv_main,
                                                   dataTrees,
                                                   eppi0DataTrees,
                                                   currentStudyGenMcTrees,
                                                   currentStudyRecMcTrees,
                                                   eppi0GenMcTrees,
                                                   eppi0RecMcTrees,
                                                   current_opts)) {
            std::cerr << "[main] ERROR: update_current_dependence_factors_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // // --------- eppi0 AAOGEN data/MC normalization + normalized raw yields ----------
    // {
    //     const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

    //     Eppi0NormalizationOptions norm_opts;
    //     norm_opts.charge_csv_path = "imports/integrated_luminosity/global.csv";
    //     norm_opts.combined_cuts_json = "output/jsons/combined_cuts.json";
    //     norm_opts.normalization_json_path = "imports/eppi0_aao_normalization_inputs.json";
    //     norm_opts.output_dir = "output/data_mc_normalization";
    //     norm_opts.override_to_unity = true;
    //     norm_opts.max_workers = 5;

    //     // Set this to true to disable the eppi0 AAOGEN normalization study.
    //     // norm_opts.override_to_unity = true;

    //     if (!update_eppi0_normalization_csv(csv_main,
    //                                         dataTrees,
    //                                         eppi0DataTrees,
    //                                         eppi0RecMcTrees,
    //                                         norm_opts)) {
    //         std::cerr << "[main] ERROR: update_eppi0_normalization_csv failed.\n";
    //         std::exit(EXIT_FAILURE);
    //     }
    // }

    // --------- pi0 contamination (helicity-averaged; bin-by-bin) ----------
    {
        const std::string csv_main  = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";

        makeOutputDirs();

        try {
            std::filesystem::copy_file(
                csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_pi0_contamination.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_pi0_contamination.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup failed (" << e.what() << "). Continuing.\n";
        }

        const int max_workers = 1;

        Pi0ContaminationOptions pi0_contamination_opts;
        pi0_contamination_opts.use_epg_mc_current_factor_for_eppi0_bkg =
            use_epg_mc_current_factor_for_eppi0_bkg;

        if (!compute_pi0_contamination_overall(
                dataTrees,
                eppi0DataTrees,
                eppi0RecMcTrees,
                eppi0BkgTrees,
                cuts_json,
                csv_main,
                output_root,
                max_workers,
                pi0_contamination_opts))
        {
            std::cerr << "[main] ERROR: compute_pi0_contamination_overall failed.\n";
            std::exit(EXIT_FAILURE);
        }
        std::cout << "pi0 contamination stage finished.\n";
    }

    // --------- Pi0-corrected DVCS signal yields (CSV + plots) ----------
    {
        const std::string csv_main   = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string csv_backup = "output/csvs/dvcs_pass2_analysis_backup_signal_yields.csv";

        // Make a backup before modifying the signal-yield columns
        try {
            std::filesystem::copy_file(csv_main, csv_backup,
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_signal_yields.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup for signal yields failed ("
                      << e.what() << "). Continuing.\n";
        }

        if (!update_pi0_corrected_counts_csv(csv_main, output_root)) {
            std::cerr << "[main] ERROR: update_pi0_corrected_counts_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }


    // --------- Pi0-subtracted direct count-based beam-spin asymmetries ----------
    //
    // This stage must run after:
    //   1. update_total_counts_csv(...),
    //   2. update_current_dependence_factors_csv(...),
    //   3. compute_pi0_contamination_overall(...), and
    //   4. update_pi0_corrected_counts_csv(...).
    //
    // The BSA stage uses helicity-separated measured ep->epgamma and ep->eppi0
    // event counts, plus finalized CSV contamination/normalized-yield columns,
    // to form S+ = G+ - f_pi0 P+ and S- = G- - f_pi0 P-.
    {
        const std::string csv_main  = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";

        try {
            std::filesystem::copy_file(
                csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_bsa.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_bsa.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup for BSA failed ("
                      << e.what() << "). Continuing.\n";
        }

        BSAOptions bsa_opts;
        bsa_opts.csv_path = csv_main;
        bsa_opts.combined_cuts_json = cuts_json;
        bsa_opts.output_root = output_root;
        bsa_opts.beam_pol_sp18_inb = 0.8882;
        bsa_opts.beam_pol_sp18_out = 0.8882;
        bsa_opts.beam_pol_fa18_inb = 0.8592;
        bsa_opts.beam_pol_fa18_out = 0.8922;
        bsa_opts.beam_pol_sp19_inb = 0.8453;
        bsa_opts.enable_pi0_subtraction = true;
        bsa_opts.make_plots = true;
        bsa_opts.max_workers = 1;

        if (!update_bsa_counts_csv(dataTrees, eppi0DataTrees, bsa_opts)) {
            std::cerr << "[main] ERROR: update_bsa_counts_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // // --------- Pi0-subtracted DVCS kinematic DATA/MC shape comparisons ----------
    // {
    //     const std::string csv_main   = "output/csvs/dvcs_pass2_analysis.csv";
    //     const std::string cuts_json  = "output/jsons/combined_cuts.json";
    //     const std::string out_dir    = "output/pi0_subtracted_dvcs_kinematics";

    //     if (!plot_pi0_subtracted_dvcs_kinematics(csv_main,
    //                                              dataTrees,
    //                                              recMcTrees,
    //                                              cuts_json,
    //                                              out_dir,
    //                                              /*max_workers=*/5)) {
    //         std::cerr << "[main] ERROR: plot_pi0_subtracted_dvcs_kinematics failed.\n";
    //         std::exit(EXIT_FAILURE);
    //     }
    // }

    // --------- Yield totals by current ----------
    {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";
        const std::string output_txt = "output/yield_totals/yield_totals_by_current.txt";

        if (!compute_yield_totals(csv_main,
                                  dataTrees, eppi0DataTrees,
                                  cuts_json,
                                  output_txt)) {
            std::cerr << "[main] ERROR: compute_yield_totals failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // // // // --------- Data/MC comparison ----------
    // // // runAllBranchDataMcComparisons(
    // // //     dataTrees,
    // // //     recMcTrees,
    // // //     eppi0DataTrees,
    // // //     eppi0RecMcTrees,
    // // //     "output/jsons/combined_cuts.json",
    // // //     "output"
    // // // );

    // --------- DVCS MC acceptance (CSV + plots) ----------
    {
        const std::string csv_main          = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string combined_cuts_json = "output/jsons/combined_cuts.json";
        const std::string global_cuts_json   = "output/jsons/global_cuts_config.json";
        const std::string csv_backup        = "output/csvs/dvcs_pass2_analysis_backup_acceptance.csv";

        // Make a backup before modifying the acceptance columns
        try {
            std::filesystem::copy_file(csv_main, csv_backup,
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_acceptance.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup for acceptance failed ("
                      << e.what() << "). Continuing.\n";
        }

        if (!update_acceptance_csv(csv_main,
                                   genMcTrees,
                                   recMcTrees,
                                   combined_cuts_json,
                                   global_cuts_json,
                                   output_root)) {
            std::cerr << "[main] ERROR: update_acceptance_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // --------- Unfolding: acceptance-corrected DVCS yields ----------
    {
        const std::string csv_main   = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string csv_backup = "output/csvs/dvcs_pass2_analysis_backup_unfolding.csv";

        // Make a backup before modifying the unfolded-yield columns
        try {
            std::filesystem::copy_file(csv_main, csv_backup,
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_unfolding.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup for unfolding failed ("
                      << e.what() << "). Continuing.\n";
        }

        if (!update_unfolded_yields_csv(csv_main, output_root)) {
            std::cerr << "[main] ERROR: update_unfolded_yields_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // // // --------- Radiative corrections (Frad factors) ----------
    // // {
    // //     const std::string csv_main   = "output/csvs/dvcs_pass2_analysis.csv";
    // //     const std::string csv_backup = "output/csvs/dvcs_pass2_analysis_backup_radcorr.csv";

    // //     // Make a backup before modifying the radiative-correction columns
    // //     try {
    // //         std::filesystem::copy_file(
    // //             csv_main,
    // //             csv_backup,
    // //             std::filesystem::copy_options::overwrite_existing
    // //         );
    // //         std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_radcorr.csv\n";
    // //     } catch (const std::exception& e) {
    // //         std::cerr << "[main] WARNING: backup for radiative corrections failed ("
    // //                   << e.what() << "). Continuing.\n";
    // //     }

    // //     if (!update_radiative_corrections_csv(csv_main, genMcTrees, radGenMcTrees, output_root)) {
    // //         std::cerr << "[main] ERROR: update_radiative_corrections_csv failed.\n";
    // //         std::exit(EXIT_FAILURE);
    // //     }
    // // }

    // --------- Kinematic bin volumes into CSV + plots ----------
    {
        const std::string csv_main     = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string out_root_dir = "output";

        // Make a backup specific to the bin-volume step
        try {
            std::filesystem::copy_file(
                csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_bin_volume.csv",
                std::filesystem::copy_options::overwrite_existing
            );
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: failed to backup CSV for bin_volume: "
                      << e.what() << "\n";
        }

        if (!update_bin_volume_csv(csv_main, out_root_dir)) {
            std::cerr << "[main] ERROR: bin_volume step failed.\n";
            return 1;
        }
    }

    // // // --------- Bin-centering corrections (Fbin) into CSV + debug plots ----------
    // // {
    // //     const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

    // //     // Make a dedicated backup before modifying Fbin columns.
    // //     try {
    // //         std::filesystem::copy_file(
    // //             csv_main,
    // //             "output/csvs/dvcs_pass2_analysis_backup_bin_centering.csv",
    // //             std::filesystem::copy_options::overwrite_existing
    // //         );
    // //         std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_bin_centering.csv\n";
    // //     } catch (const std::exception& ex) {
    // //         std::cerr << "[main] WARNING: failed to create bin-centering backup: "
    // //                   << ex.what() << "\n";
    // //     }

    // //     ModelPaths model_paths;   // use env/defaults for dvcsgen and km15_cli
    // //     const bool vgg_globalfit = false;  // set true if you want --globalfit for VGG
    // //     const int  n_steps       = 3;      // sub-bins per dimension (xB,Q2,t,phi)

    // //     if (!update_bin_centering_corrections_csv(
    // //             csv_main,
    // //             n_steps,
    // //             model_paths,
    // //             vgg_globalfit,
    // //             ModelChoice::Both)) {
    // //         std::cerr << "[main] ERROR: bin-centering corrections failed.\n";
    // //         return 1;
    // //     }

    // //     // Debug plots: Fbin vs phi for 10.6 and 10.2 GeV.
    // //     // Uses Fbin triples and phiavg columns from the updated CSV.
    // //     plot_bin_centering_fbin_vs_phi(
    // //         csv_main,
    // //         "output/bin_centering_plots");
    // // }


    // // {
    // //     const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";
    // //     const std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    // //     const std::string outdir = "output/propagator_study";

    // //     // Build the exact keying scheme required by propagator_study.cpp:
    // //     // required keys: fa18_inb, fa18_out, sp18_inb, sp18_out, sp19_inb
    // //     std::map<std::string, std::vector<TTree*>> dataTreesByPeriod;

    // //     auto require_tree = [&](const std::string& in_key,
    // //                             const std::string& out_key) {
    // //         auto it = dataTrees.find(in_key);
    // //         if (it == dataTrees.end() || it->second == nullptr) {
    // //             std::cerr << "[main] FATAL: missing required data tree key \"" << in_key << "\"\n";
    // //             std::cerr << "[main] Available dataTrees keys are:\n";
    // //             for (const auto& kv : dataTrees) {
    // //                 std::cerr << "  - " << kv.first << "\n";
    // //             } //endfor
    // //             std::exit(EXIT_FAILURE);
    // //         }
    // //         dataTreesByPeriod[out_key].push_back(it->second);
    // //     };

    // //     // These input keys must match exactly whatever loadTrees() uses.
    // //     // Based on your comment earlier ("keys like DVCS_Fa18_inb, ..."), this is the expected mapping:
    // //     require_tree("DVCS_Fa18_inb", "fa18_inb");
    // //     require_tree("DVCS_Fa18_out", "fa18_out");
    // //     require_tree("DVCS_Sp18_inb", "sp18_inb");
    // //     require_tree("DVCS_Sp18_out", "sp18_out");
    // //     require_tree("DVCS_Sp19_inb", "sp19_inb");

    // //     if (!propagator_study::run_propagator_study(csv_main,
    // //                                                 dataTreesByPeriod,
    // //                                                 combined_cuts_json,
    // //                                                 outdir)) {
    // //         std::cerr << "[main] ERROR: propagator study failed\n";
    // //         return 1;
    // //     }
    // // }


    // --------- Cross sections (CSV update + theory JSON + plots) ----------
    {
        const std::string csv_main         = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string theory_json_root = "output/jsons/cross_sections";
        const std::string xs_out_root      = "output/cross_sections";

        // // --------- Theory grids (xs_phi_all.json generation) ----------
        // {
        //     const std::string csv_main         = "output/csvs/dvcs_pass2_analysis.csv";
        //     const std::string theory_json_root = "output/jsons/cross_sections";
        //
        //     if (!regenerate_theory_jsons(csv_main, theory_json_root)) {
        //         std::cerr << "[main] ERROR: regenerate_theory_jsons failed.\n";
        //         return 1;
        //     }
        // }

        LumiBuildOptions lumi_opts;

        // Fallback mode if the Fa18/Sp19 columns-3-to-5 mode below is disabled.
        lumi_opts.use_second_column_charge_for_all_unpolarized = true;

        // New mode:
        //   Sp18 Inb / Sp18 Out:
        //     L_unpol = column 2
        //
        //   Fa18 Inb / Fa18 Out / Sp19 Inb:
        //     L_unpol = 1.025 * (column 3 + column 4 + column 5)
        //
        // Polarized luminosities remain:
        //   pos -> column 3
        //   neg -> column 4
        lumi_opts.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized = true;
        lumi_opts.columns_3_to_5_charge_sum_scale = 1.025;

        LumiMap lumi_map = build_lumi_map(lumi_opts);

        if (!compute_cross_sections(csv_main, lumi_map)) {
            std::cerr << "[main] ERROR: compute_cross_sections failed.\n";
        }

        const std::vector<std::string> labels_to_plot = {
            "Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp",
            "Sp18 Inb", "Sp18 Out", "Sp19 Inb",
            "Fa18", "Sp18", "10.6 GeV"
        };

        for (const auto &label : labels_to_plot) {
            if (!plot_cross_sections_for_label(csv_main, label,
                theory_json_root, xs_out_root)) {
                std::cerr << "[main] WARNING: plot_cross_sections_for_label failed for "
                          << label << "\n";
            }
        }
    }


    // --------- Overall BH-edge normalization study ----------
    {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        OverallNormalizationOptions norm_opts;

        // Production-safe default:
        //   true  -> write norm = 1.00 everywhere and skip BH-edge study
        //   false -> run the full BH-edge normalization study and write fitted norms
        norm_opts.override_to_unity = true;

        norm_opts.use_all_points_within_edge_window = true;
        norm_opts.require_positive_dedge = true;
        norm_opts.max_dedge_for_normalization_deg = 10.0;
        norm_opts.norm_x_axis = OverallNormXAxis::XB;
        norm_opts.output_dir = "output/normalization_study";

        const std::vector<std::string> norm_labels = {
            "Fa18 Inb",
            "Fa18 Out",
            "Sp19 Inb",
            "Sp18 Inb",
            "Sp18 Out",
            "Fa18",
            "Sp18",
            "10.6 GeV"
        };

        for (const std::string& label : norm_labels) {
            if (!update_overall_normalization_study_csv(csv_main,
                                                        label,
                                                        "unpol",
                                                        norm_opts)) {
                std::cerr << "[main] ERROR: update_overall_normalization_study_csv failed for "
                          << label << ".\n";
                std::exit(EXIT_FAILURE);
            }
        }
    }

    // --------- DVCS normalized cross sections (CSV + plots) ----------
    {
        const std::string csv_main           = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string theory_json_root   = "output/jsons/cross_sections";
        const std::string out_norm_xsec_root = "output/normed_cross_sections_plots";

        if (!update_normed_cross_sections_csv(csv_main)) {
            std::cerr << "[main] FATAL: update_normed_cross_sections_csv failed.\n";
            return 1;
        }

        const std::vector<std::string> labels = {
            "Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out", "Sp19 Inb",
            "Fa18", "Sp18", "10.6 GeV"
        };

        for (const auto &lab : labels) {
            if (!plot_normed_cross_sections_for_label(csv_main,
                                                      lab,
                                                      theory_json_root,
                                                      out_norm_xsec_root)) {
                std::cerr << "[main] FATAL: plot_normed_cross_sections_for_label failed for "
                          << lab << "\n";
                return 1;
            }
        }
    }

    std::cout << "All done." << std::endl;
    return 0;
}