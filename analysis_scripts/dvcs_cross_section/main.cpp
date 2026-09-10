#include "initialize_pass2_csv.h"
#include "make_dirs.h"
#include "load_trees.h"
#include "periods.h"
#include "global_cuts.h"
#include "python_exclusivity_runner.h"
#include "load_binning_scheme.h"
#include "bin_means.h"
#include "total_counts.h"
#include "current_dependence.h"
#include "current_systematics.h"
#include "eppi0_normalization.h"
#include "pi0_contamination.h"
#include "pi0_corrected_counts.h"
#include "bsa.h"
#include "radiative_corrections.h"
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <sstream>
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
#include "external_scripts_runner.h"
#include "systematics_runner.h"
#include "cut_variation_runner.h"


namespace {

struct SystematicRunSelection {
    bool cuts = true;
    bool current = true;
    bool acceptance = true;
    bool csv_only = true;
    bool external_studies = true;
};

static std::vector<std::string> split_tokens(const std::string& value) {
    std::vector<std::string> out;
    std::stringstream ss(value);
    std::string token;
    while(std::getline(ss,token,',')){
        token.erase(
            std::remove_if(token.begin(),token.end(),
                           [](unsigned char c){return std::isspace(c);}),
            token.end());
        if(!token.empty()) out.push_back(token);
    }
    return out;
}

static SystematicRunSelection parse_systematic_selection(
    int argc,char* argv[]) {

    SystematicRunSelection sel;
    bool explicitly_set=false;

    for(int i=1;i<argc;++i){
        const std::string arg=argv[i];

        if(arg=="--skip-systematics"){
            sel.cuts=false;
            sel.current=false;
            sel.acceptance=false;
            sel.csv_only=false;
            sel.external_studies=false;
            explicitly_set=true;
            continue;
        }

        if(arg=="--systematics"){
            if(i+1>=argc){
                throw std::runtime_error(
                    "--systematics requires a comma-separated value");
            }

            sel.cuts=false;
            sel.current=false;
            sel.acceptance=false;
            sel.csv_only=false;
            sel.external_studies=false;
            explicitly_set=true;

            for(const std::string& token:split_tokens(argv[++i])){
                if(token=="all"){
                    sel.cuts=sel.current=sel.acceptance=sel.csv_only=sel.external_studies=true;
                } else if(token=="none"){
                    sel.cuts=sel.current=sel.acceptance=sel.csv_only=sel.external_studies=false;
                } else if(token=="cuts" ||
                          token=="exclusivity" ||
                          token=="fiducial"){
                    // The automatic cut runner produces exclusivity and
                    // fiducial variations together because they share the
                    // same cut-dependent extraction machinery.
                    sel.cuts=true;
                } else if(token=="current"){
                    sel.current=true;
                } else if(token=="acceptance" || token=="acceptance-reweighting"){
                    sel.acceptance=true;
                } else if(token=="csv" || token=="post"){
                    sel.csv_only=true;
                } else if(token=="external"){
                    sel.external_studies=true;
                } else {
                    throw std::runtime_error(
                        "unknown --systematics token: "+token);
                }
            }
        }
    }

    if(!explicitly_set){
        // Preserve the historical/default behavior: run everything.
        sel=SystematicRunSelection{};
    }

    return sel;
}

} // namespace

int main(int argc, char* argv[]) {
    bool acceptance_reweighting_only = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--acceptance-reweighting-only") {
            acceptance_reweighting_only = true;
        }
    }
    SystematicRunSelection systematic_selection;
    try {
        systematic_selection = parse_systematic_selection(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "[main] FATAL: " << e.what() << "\n";
        std::cerr << "Usage examples:\n"
                  << "  ./dvcs_analysis\n"
                  << "  ./dvcs_analysis --systematics current\n"
                  << "  ./dvcs_analysis --systematics cuts,current,acceptance\n"
                  << "  ./dvcs_analysis --systematics acceptance,csv\n"
                  << "  ./dvcs_analysis --systematics csv\n"
                  << "  ./dvcs_analysis --skip-systematics\n"
                  << "  ./dvcs_analysis --acceptance-reweighting-only\n";
        return 1;
    }

    std::cout << "Starting DVCS analysis..." << std::endl;
    std::cout << "[main] Systematics selection:"
              << " cuts=" << systematic_selection.cuts
              << " current=" << systematic_selection.current
              << " acceptance=" << systematic_selection.acceptance
              << " csv=" << systematic_selection.csv_only
              << " external=" << systematic_selection.external_studies
              << std::endl;

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

    // Integrated-analysis detector-quality exclusions for Sp18 Out. These are
    // enabled by default in GlobalCutConfig and are automatically suspended if
    // any topology or particle-sector study switch below is enabled.
    global_cfg.enable_sp18_out_sector_quality_cuts = true;

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
    global_cfg.enable_proton_cd_sector_filter = false;
    global_cfg.proton_cd_sector = 1;

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

    // Record the exact global-cut configuration before the Python process is
    // launched. The Python exclusivity stage must apply the same minimally
    // selected event population used by every downstream C++ stage.
    write_global_cuts_config_json("output/jsons", global_cfg);

    // Run the production Python template-fit exclusivity extraction before
    // opening the ROOT trees in the C++ process. The runner validates and
    // installs the 90%, 95%, and 98% cut JSONs under output/jsons, with the
    // nominal 95% result also installed as combined_cuts.json.
    PythonExclusivityOptions exclusivity_opts;
    exclusivity_opts.enabled = true;
    exclusivity_opts.force_rerun = true;
    exclusivity_opts.python_executable = "python3";
    exclusivity_opts.script_path = "plot_exclusivity_data_dvcs_pi0_mc.py";
    exclusivity_opts.global_cuts_json = "output/jsons/global_cuts_config.json";
    exclusivity_opts.output_directory = "output/exclusivity_fit";
    exclusivity_opts.install_directory = "output/jsons";
    exclusivity_opts.workers = 7;
    exclusivity_opts.tight_containment = 0.90;
    exclusivity_opts.nominal_containment = 0.95;
    exclusivity_opts.loose_containment = 0.98;

    if (!acceptance_reweighting_only) {
        if (!run_python_exclusivity_analysis(exclusivity_opts)) {
            std::cerr << "[main] FATAL: Python exclusivity optimization failed.\n";
            return 1;
        }

        std::cout << "[main] Python exclusivity-cut stage finished. "
                  << "Using nominal cuts from output/jsons/combined_cuts.json.\n";

        initialize_pass2_csv("imports/all_bin_v3.csv",
                             "output/csvs/dvcs_pass2_analysis.csv");
    } else {
        std::cout << "[main] Acceptance-reweighting-only mode: reusing existing "
                  << "combined cuts, current calibration, and analysis CSV.\n";
    }

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

    if (acceptance_reweighting_only) {
        AcceptanceReweightingOptions arw;
        arw.combined_cuts_json = "output/jsons/combined_cuts.json";
        arw.current_response_model_json =
            "output/dvcs_current_dependence/calibration/current_response_model.json";
        arw.output_dir = "output/systematics/acceptance_reweighting";
        arw.enable_bh_reweighting = false;
        arw.build_bh_grid_if_missing = false;
        arw.install_candidate_as_production_systematic = false;

        if (!run_acceptance_reweighting_study(
                "output/csvs/dvcs_pass2_analysis.csv",
                dataTrees, genMcTrees, recMcTrees, arw)) {
            std::cerr << "[main] FATAL: acceptance reweighting study failed.\n";
            return 1;
        }
        std::cout << "[main] Acceptance reweighting study finished. "
                  << "No production systematic was replaced.\n";
        return 0;
    }

    // --------- Global bin-averaged kinematics (CSV update) ----------
    {
        const std::string csv_main   = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string csv_backup = "output/csvs/dvcs_pass2_analysis_backup_bin_means.csv";

        try {
            std::filesystem::copy_file(csv_main, csv_backup,
                                       std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to " << csv_backup << "\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: Backup failed (" << e.what() << "). Continuing anyway.\n";
        }

        BinMeansOptions bin_mean_opts;
        bin_mean_opts.make_note_outputs = true;
        bin_mean_opts.note_output_dir = "output/bin_means/analysis_note";

        if (!update_bin_means_csv(csv_main, dataTrees, /*max_workers=*/7, bin_mean_opts)) {
            std::cerr << "[main] ERROR: update_bin_means_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    const bool use_nobkg_dvcs_mc_for_acceptance = false;
    const bool use_epg_mc_current_factor_for_eppi0_bkg = true;

    // --------- Current-response calibration + diagnostics ----------
    {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        CurrentDependenceOptions current_opts;
        current_opts.charge_csv_path = "imports/integrated_luminosity/global.csv";
        current_opts.combined_cuts_json = "output/jsons/combined_cuts.json";
        current_opts.output_dir = "output/dvcs_current_dependence";

        current_opts.override_to_unity = false;
        current_opts.use_second_column_charge_for_all_unpolarized = true;
        current_opts.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized = false;
        current_opts.columns_3_to_5_charge_sum_scale = 1.025;
        current_opts.use_fa18_inb_current_efficiency_for_sp19_inb = true;
        current_opts.use_nobkg_dvcs_mc_counts = use_nobkg_dvcs_mc_for_acceptance;
        current_opts.enable_photon_region_current_diagnostic = true;
        current_opts.enable_eppi0_photon_region_current_diagnostic = true;
        current_opts.enable_region_theta_current_diagnostic = true;
        current_opts.enable_eppi0_region_theta_current_diagnostic = false;
        current_opts.enable_exploratory_kinematic_current_diagnostic = false;
        current_opts.enable_relative_ft_fd_photon_efficiency_diagnostic = false;
        current_opts.enable_fa18_sp19_transfer_shape_diagnostic = true;
        current_opts.sp19_transfer_photon_energy_max_GeV = 10.2;
        current_opts.finalized_production_mode = true;
        current_opts.clean_output_dir_before_run = true;

        // Final production current-response prescription:
        //   - regional FT/S1--S6 DATA response for every period/channel,
        //   - plus one common centered linear electron-angle term for
        //     Sp18 Out ep->epg only,
        //   - no legacy post-binning current correction.
        current_opts.use_sp18_out_e_theta_response_model = true;
        current_opts.use_e_theta_linear_data_current_efficiency = false;
        current_opts.response_model_json = "output/dvcs_current_dependence/calibration/current_response_model.json";
        current_opts.apply_legacy_binned_current_corrections = false;
        current_opts.use_epg_mc_current_factor_for_eppi0_bkg = use_epg_mc_current_factor_for_eppi0_bkg;
        current_opts.max_workers = 7;

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


    // --------- Raw yields + event-level current-corrected yields ----------
    {
        const std::string csv_main  = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string cuts_json = "output/jsons/combined_cuts.json";

        try {
            std::filesystem::copy_file(csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_total_counts.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_total_counts.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup failed (" << e.what() << "). Continuing.\n";
        }

        TotalCountsOptions total_count_opts;
        total_count_opts.use_nobkg_dvcs_mc_counts = use_nobkg_dvcs_mc_for_acceptance;
        // The large per-period/topology diagnostic canvases are redundant with
        // the compact analysis-note yield summaries and cost substantial ROOT
        // drawing/output time in every production run.
        total_count_opts.make_plots = false;
        total_count_opts.make_note_outputs = true;
        total_count_opts.apply_event_level_current_correction = true;
        total_count_opts.current_response_model_json = "output/dvcs_current_dependence/calibration/current_response_model.json";
        total_count_opts.require_sp18_out_epg_e_theta_current_model = true;
        total_count_opts.use_epg_mc_current_factor_for_eppi0_bkg =
            use_epg_mc_current_factor_for_eppi0_bkg;
        total_count_opts.write_current_nuisance_responses = true;
        total_count_opts.current_nuisance_response_csv =
            "output/current_systematics/current_yield_nuisance_responses.csv";

        if (!update_total_counts_csv(csv_main, dataTrees, eppi0DataTrees,
            genMcTrees, recMcTrees,
            eppi0GenMcTrees, eppi0RecMcTrees,
            eppi0BkgTrees,
            cuts_json,
            output_root,
            /*max_workers=*/7,
            total_count_opts,
            currentStudyGenMcTrees,
            currentStudyRecMcTrees)) {
            std::cerr << "[main] ERROR: update_total_counts_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }
    }

    // // --------- eppi0 AAOGEN data/MC normalization + normalized raw yields ----------
    // {
    //     const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";
    //
    //     Eppi0NormalizationOptions norm_opts;
    //     norm_opts.charge_csv_path = "imports/integrated_luminosity/global.csv";
    //     norm_opts.combined_cuts_json = "output/jsons/combined_cuts.json";
    //     norm_opts.normalization_json_path = "imports/eppi0_aao_normalization_inputs.json";
    //     norm_opts.output_dir = "output/data_mc_normalization";
    //     norm_opts.override_to_unity = true;
    //     norm_opts.max_workers = 7;
    //
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

        try {
            std::filesystem::copy_file(
                csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_pi0_contamination.csv",
                std::filesystem::copy_options::overwrite_existing);
            std::cout << "[main] Backed up CSV to dvcs_pass2_analysis_backup_pi0_contamination.csv\n";
        } catch (const std::exception& e) {
            std::cerr << "[main] WARNING: backup failed (" << e.what() << "). Continuing.\n";
        }

        const int max_workers = 7;

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
    //   1. update_current_dependence_factors_csv(...),
    //   2. update_total_counts_csv(...),
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
        bsa_opts.max_workers = 7;

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
    //                                              /*max_workers=*/7)) {
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
        const std::string csv_main           = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string combined_cuts_json = "output/jsons/combined_cuts.json";
        const std::string global_cuts_json   = "output/jsons/global_cuts_config.json";
        const std::string csv_backup         = "output/csvs/dvcs_pass2_analysis_backup_acceptance.csv";

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

    // --------- Radiative-correction analysis-note outputs ----------
    // Frad is model-only and unchanged from pass-1.  Production values are
    // taken from imports/all_bin_v3.csv; the old MC recalculation remains
    // available as a cross-check but is not run in the standard pass-2 chain.
    {
        if (!write_radiative_corrections_analysis_note_outputs(
                "imports/all_bin_v3.csv",
                "output/radiative_corrections")) {
            std::cerr << "[main] ERROR: radiative-correction analysis-note outputs failed.\n";
            return 1;
        }
    }

    // --------- Kinematic bin volumes into CSV + plots ----------
    {
        const std::string csv_main     = "output/csvs/dvcs_pass2_analysis.csv";
        const std::string out_root_dir = "output";

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

    // --------- Pass-2 bin-centering correction ----------
    // Recompute the central KM15 factor using the actual pass-2 mean
    // kinematics.  The fast backend batches the entire CSV through a small
    // fixed number of persistent Gepard/KM15 worker processes.
    {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        try {
            std::filesystem::copy_file(
                csv_main,
                "output/csvs/dvcs_pass2_analysis_backup_bin_centering.csv",
                std::filesystem::copy_options::overwrite_existing);
        } catch (const std::exception& ex) {
            std::cerr << "[main] WARNING: failed to create bin-centering backup: "
                      << ex.what() << "\n";
        }

        ModelPaths model_paths;
        const int quadrature_order = 3;  // 3^4 = 81 KM15 points per populated bin.
        if (!update_bin_centering_corrections_csv(
                csv_main, quadrature_order, model_paths, false, ModelChoice::KM15Only)) {
            std::cerr << "[main] ERROR: bin-centering corrections failed.\n";
            return 1;
        }

        plot_bin_centering_fbin_vs_phi(csv_main, "output/bin_centering_plots");
    }


    // // {
    // //     const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";
    // //     const std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    // //     const std::string outdir = "output/propagator_study";
    //
    // //     std::map<std::string, std::vector<TTree*>> dataTreesByPeriod;
    //
    // //     auto require_tree = [&](const std::string& in_key,
    // //                             const std::string& out_key) {
    // //         auto it = dataTrees.find(in_key);
    // //         if (it == dataTrees.end() || it->second == nullptr) {
    // //             std::cerr << "[main] FATAL: missing required data tree key \"" << in_key << "\"\n";
    // //             std::cerr << "[main] Available dataTrees keys are:\n";
    // //             for (const auto& kv : dataTrees) {
    // //                 std::cerr << "  - " << kv.first << "\n";
    // //             }
    // //             std::exit(EXIT_FAILURE);
    // //         }
    // //         dataTreesByPeriod[out_key].push_back(it->second);
    // //     };
    //
    // //     require_tree("DVCS_Fa18_inb", "fa18_inb");
    // //     require_tree("DVCS_Fa18_out", "fa18_out");
    // //     require_tree("DVCS_Sp18_inb", "sp18_inb");
    // //     require_tree("DVCS_Sp18_out", "sp18_out");
    // //     require_tree("DVCS_Sp19_inb", "sp19_inb");
    //
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
        lumi_opts.charge_csv_path =
            "imports/integrated_luminosity/global.csv";

        LumiMap lumi_map = build_lumi_map(lumi_opts);

        if (!compute_cross_sections(csv_main, lumi_map)) {
            std::cerr << "[main] ERROR: compute_cross_sections failed.\n";
        } else if (!write_cross_section_analysis_note_outputs(csv_main, lumi_map)) {
            std::cerr << "[main] WARNING: cross-section analysis-note output generation failed.\n";
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


    // --------- Current-dependent efficiency systematic ----------
    //
    // Fast post-processing of the signed calibration-nuisance responses written
    // by total_counts.cpp.  This does not rerun event loops.  CSV-only finalization
    // depends on these columns, so regenerate them automatically whenever csv_only
    // is requested, even if the user did not explicitly select current.
    if (systematic_selection.current || systematic_selection.csv_only) {
        CurrentSystematicsOptions current_syst_opts;
        current_syst_opts.nuisance_response_csv =
            "output/current_systematics/current_yield_nuisance_responses.csv";
        current_syst_opts.transfer_shape_csv =
            "output/dvcs_current_dependence/analysis_note/current_systematics/"
            "fa18_sp19_50nA_shape_distance.csv";
        current_syst_opts.output_dir =
            "output/current_systematics/analysis_note";
        current_syst_opts.write_csv_columns = true;

        if (!evaluate_current_dependence_systematics(
                "output/csvs/dvcs_pass2_analysis.csv",
                current_syst_opts)) {
            std::cerr
                << "[main] ERROR: current-dependence systematic failed.\n";
            return 1;
        }
    } else {
        std::cout
            << "[main] Skipping current-dependence systematic by request.\n";
    }

    // --------- Automatic exclusivity/fiducial cut point-to-point systematics ----------
    // One top-level switch controls all four nonnominal selections. The runner
    // clones the completed nominal CSV, recomputes only the cut-dependent chain,
    // applies Barlow B >= 1, writes raw/final absolute uncertainties in nb/GeV^4,
    // and produces bin-by-bin diagnostic canvases.
    if (systematic_selection.cuts) {
        AutomaticCutVariationOptions cut_variation_opts;
        cut_variation_opts.enabled = true;
        cut_variation_opts.make_exclusivity_extraction_plots = false;
        cut_variation_opts.make_final_diagnostic_plots = true;
        cut_variation_opts.use_pass1_tight_instability_rule = true;
        cut_variation_opts.tight_relative_difference_threshold = 0.50;
        cut_variation_opts.max_workers = 7;
        cut_variation_opts.nominal_csv = "output/csvs/dvcs_pass2_analysis.csv";
        cut_variation_opts.output_dir = "output/cut_variation_systematics";

        if (!run_automatic_cut_variation_systematics(
                cut_variation_opts,
                global_cfg,
                dataTrees, genMcTrees, recMcTrees,
                eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees, eppi0BkgTrees,
                currentStudyGenMcTrees, currentStudyRecMcTrees,
                use_nobkg_dvcs_mc_for_acceptance,
                use_epg_mc_current_factor_for_eppi0_bkg)) {
            std::cerr << "[main] ERROR: automatic cut-variation systematics failed.\n";
            return 1;
        }
    }

    else {
        std::cout
            << "[main] Skipping automatic exclusivity/fiducial variations by request.\n";
    }

    // --------- Acceptance model-dependence systematic ----------
    // Derive the pass-2 acceptance-model uncertainty from the nominal ->
    // DATA-reweighted acceptance excursion and the synthetic transfer-closure
    // envelope.  This tree-based stage writes fractional candidate columns;
    // main_systematics subsequently converts them to the absolute production
    // acceptance uncertainties and recomputes the point-to-point totals.
    if (systematic_selection.acceptance) {
        AcceptanceReweightingOptions arw;
        arw.combined_cuts_json = "output/jsons/combined_cuts.json";
        arw.current_response_model_json =
            "output/dvcs_current_dependence/calibration/current_response_model.json";
        arw.output_dir = "output/systematics/acceptance_reweighting";
        arw.enable_bh_reweighting = false;
        arw.build_bh_grid_if_missing = false;
        arw.install_candidate_as_production_systematic = false;

        if (!run_acceptance_reweighting_study(
                "output/csvs/dvcs_pass2_analysis.csv",
                dataTrees, genMcTrees, recMcTrees, arw)) {
            std::cerr << "[main] FATAL: acceptance reweighting study failed.\n";
            return 1;
        }
        std::cout << "[main] Acceptance reweighting candidate written; "
                  << "production assignment will be finalized by main_systematics.\n";
    } else {
        std::cout << "[main] Skipping acceptance-reweighting systematic by request.\n";
    }

    // --------- CSV-only systematic uncertainties ----------
    // Run this after all nominal and cut-variation columns have been written,
    // but before any external study reads the analysis CSV. This guarantees
    // that the legacy/integrated scripts use systematics produced from the
    // current analysis configuration rather than columns left by an older run.
    if (systematic_selection.csv_only) {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        SystematicsRunnerOptions systematics_opts;
        systematics_opts.executable = "./main_systematics";
        systematics_opts.pass1_systematics_csv =
            "imports/pass1_systematic_summary.csv";

        if (!run_main_systematics(csv_main, systematics_opts)) {
            std::cerr << "[main] FATAL: main_systematics failed.\n";
            return 1;
        }
    }

    else {
        std::cout
            << "[main] Skipping CSV-only systematics by request.\n";
    }

    // --------- External integrated and legacy cross-section studies ----------
    if (systematic_selection.external_studies) {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        ExternalScriptOptions external_opts;
        external_opts.scripts_directory = "external_scripts";
        external_opts.python_executable = "python";
        external_opts.published_pass1_cross_section_table =
            "imports/clasdb_E214M1.txt";
        external_opts.include_bin_to_bin_systematics = true;
        external_opts.use_simple_clas6_cross_check = true;

        if (!run_external_cross_section_scripts(csv_main, external_opts)) {
            std::cerr << "[main] FATAL: external cross-section scripts failed.\n";
            return 1;
        }
    }

    else {
        std::cout
            << "[main] Skipping external integrated studies by request.\n";
    }

    std::cout << "All done." << std::endl;
    return 0;
}