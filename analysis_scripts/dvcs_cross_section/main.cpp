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
#include "systematics_runner.h"
#include "cut_variation_runner.h"


namespace {

struct SystematicRunSelection {
    bool cuts = true;
    bool current = true;
    bool acceptance = true;
    bool csv_only = true;
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
            explicitly_set=true;

            for(const std::string& token:split_tokens(argv[++i])){
                if(token=="all"){
                    sel.cuts=sel.current=sel.acceptance=sel.csv_only=true;
                } else if(token=="none"){
                    sel.cuts=sel.current=sel.acceptance=sel.csv_only=false;
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
    bool topology_acceptance_study = false;
    std::string topology_study_csv;
    bool eppi0_normalization_only = false;
    bool eppi0_production_test = false;
    bool use_eppi0_production_normalization = false;
    std::string topology_cli;
    int photon_sector_cli = 0;
    int electron_sector_cli = 0;
    int proton_fd_sector_cli = 0;
    int proton_cd_sector_cli = 0;
    bool disable_sp18_out_sector_quality_cuts = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--acceptance-reweighting-only") {
            acceptance_reweighting_only = true;
        } else if (std::string(argv[i]) == "--topology-acceptance-study") {
            topology_acceptance_study = true;
        } else if (std::string(argv[i]) == "--topology-study-csv") {
            if (i + 1 >= argc) { std::cerr << "[main] FATAL: --topology-study-csv requires a path\n"; return 1; }
            topology_study_csv = argv[++i];
        } else if (std::string(argv[i]) == "--eppi0-normalization-only") {
            eppi0_normalization_only = true;
        } else if (std::string(argv[i]) == "--eppi0-production-test") {
            eppi0_production_test = true;
        } else if (std::string(argv[i]) == "--eppi0-normalization") {
            use_eppi0_production_normalization = true;
        } else if (std::string(argv[i]) == "--no-eppi0-normalization") {
            // Backward-compatible explicit spelling of the production default.
            use_eppi0_production_normalization = false;
        } else if (std::string(argv[i]) == "--disable-sp18-out-sector-quality-cuts") {
            disable_sp18_out_sector_quality_cuts = true;
        } else if (std::string(argv[i]) == "--electron-sector") {
            if (i + 1 >= argc) { std::cerr << "[main] FATAL: --electron-sector requires an integer 1--6\n"; return 1; }
            try { electron_sector_cli = std::stoi(argv[++i]); }
            catch (...) { std::cerr << "[main] FATAL: --electron-sector requires an integer 1--6\n"; return 1; }
            if (electron_sector_cli < 1 || electron_sector_cli > 6) {
                std::cerr << "[main] FATAL: --electron-sector must be in [1,6]\n"; return 1;
            }
        } else if (std::string(argv[i]) == "--proton-fd-sector") {
            if (i + 1 >= argc) { std::cerr << "[main] FATAL: --proton-fd-sector requires an integer 1--6\n"; return 1; }
            try { proton_fd_sector_cli = std::stoi(argv[++i]); }
            catch (...) { std::cerr << "[main] FATAL: --proton-fd-sector requires an integer 1--6\n"; return 1; }
            if (proton_fd_sector_cli < 1 || proton_fd_sector_cli > 6) {
                std::cerr << "[main] FATAL: --proton-fd-sector must be in [1,6]\n"; return 1;
            }
        } else if (std::string(argv[i]) == "--proton-cd-sector") {
            if (i + 1 >= argc) { std::cerr << "[main] FATAL: --proton-cd-sector requires an integer 1--3\n"; return 1; }
            try { proton_cd_sector_cli = std::stoi(argv[++i]); }
            catch (...) { std::cerr << "[main] FATAL: --proton-cd-sector requires an integer 1--3\n"; return 1; }
            if (proton_cd_sector_cli < 1 || proton_cd_sector_cli > 3) {
                std::cerr << "[main] FATAL: --proton-cd-sector must be in [1,3]\n"; return 1;
            }
        } else if (std::string(argv[i]) == "--photon-sector") {
            if (i + 1 >= argc) { std::cerr << "[main] FATAL: --photon-sector requires an integer 1--6\n"; return 1; }
            try { photon_sector_cli = std::stoi(argv[++i]); }
            catch (...) { std::cerr << "[main] FATAL: --photon-sector requires an integer 1--6\n"; return 1; }
            if (photon_sector_cli < 1 || photon_sector_cli > 6) {
                std::cerr << "[main] FATAL: --photon-sector must be in [1,6]\n"; return 1;
            }
        } else if (std::string(argv[i]) == "--topology") {
            if (i + 1 >= argc) {
                std::cerr << "[main] FATAL: --topology requires one of FD-FD, CD-FD, CD-FT\n";
                return 1;
            }
            topology_cli = argv[++i];
            std::transform(topology_cli.begin(), topology_cli.end(), topology_cli.begin(),
                           [](unsigned char c){ return static_cast<char>(std::toupper(c)); });
            std::replace(topology_cli.begin(), topology_cli.end(), '_', '-');
            if (topology_cli != "FD-FD" &&
                topology_cli != "CD-FD" &&
                topology_cli != "CD-FT") {
                std::cerr << "[main] FATAL: invalid --topology '" << topology_cli
                          << "'. Allowed values: FD-FD, CD-FD, CD-FT\n";
                return 1;
            }
        }
    }
    const int n_particle_sector_filters =
        (electron_sector_cli != 0) + (proton_fd_sector_cli != 0) +
        (proton_cd_sector_cli != 0) + (photon_sector_cli != 0);
    if (n_particle_sector_filters > 1) {
        std::cerr << "[main] FATAL: choose only one of --electron-sector, --proton-fd-sector, "
                  << "--proton-cd-sector, or --photon-sector per run.\n";
        return 1;
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
                  << "  ./dvcs_analysis --acceptance-reweighting-only\n"
                  << "  ./dvcs_analysis --eppi0-normalization-only\n"
                  << "  ./dvcs_analysis --eppi0-production-test\n"
                  << "  ./dvcs_analysis --eppi0-normalization      # optional eppi0 detector-normalization study\n"
                  << "  ./dvcs_analysis --no-eppi0-normalization   # explicit spelling of default Krishna normalization\n"
                  << "  ./dvcs_analysis --skip-systematics --disable-sp18-out-sector-quality-cuts\n"
                  << "  ./dvcs_analysis --no-eppi0-normalization --skip-systematics --topology FD-FD\n"
                  << "  ./dvcs_analysis --no-eppi0-normalization --skip-systematics --topology CD-FD\n"
                  << "  ./dvcs_analysis --no-eppi0-normalization --skip-systematics --topology CD-FT\n";
        return 1;
    }

    std::cout << "Starting DVCS analysis..." << std::endl;
    std::cout << "[main] Systematics selection:"
              << " cuts=" << systematic_selection.cuts
              << " current=" << systematic_selection.current
              << " acceptance=" << systematic_selection.acceptance
              << " csv=" << systematic_selection.csv_only
              << std::endl;
    std::cout << "[main] Production detector normalization: "
              << (use_eppi0_production_normalization
                      ? "eppi0 sequential theta*p on reconstructed DVCS MC (Krishna OFF)"
                      : "Krishna/Neupane proton-efficiency correction (eppi0 OFF)")
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

    // Integrated-analysis detector-quality exclusions for Sp18 Out.
    //
    // Nominal production keeps these ON. The diagnostic command-line switch
    // --disable-sp18-out-sector-quality-cuts turns them OFF globally so the
    // inclusive parent, topology parents, and every particle-sector extraction
    // use the same unrestricted Sp18-Out detector population.
    //
    // This is deliberately explicit rather than automatically tied to a sector
    // filter: an ON/OFF study must use one common setting for the parent and all
    // daughter extractions.
    global_cfg.enable_sp18_out_sector_quality_cuts =
        !disable_sp18_out_sector_quality_cuts;

    // Optional single-topology study selected from the command line:
    //   --topology FD-FD
    //   --topology CD-FD
    //   --topology CD-FT
    // With no --topology argument the nominal inclusive configuration is kept.
    global_cfg.enable_topology_filter = !topology_cli.empty();
    if (topology_cli == "FD-FD") {
        global_cfg.required_detector1 = 1;
        global_cfg.required_detector2 = 1;
    } else if (topology_cli == "CD-FD") {
        global_cfg.required_detector1 = 2;
        global_cfg.required_detector2 = 1;
    } else if (topology_cli == "CD-FT") {
        global_cfg.required_detector1 = 2;
        global_cfg.required_detector2 = 0;
    }

    // Electron FD sector study. Electron is always FD.
    global_cfg.enable_electron_fd_sector_filter = (electron_sector_cli != 0);
    global_cfg.electron_fd_sector = (electron_sector_cli != 0 ? electron_sector_cli : 1);

    // Proton FD sector study. This automatically keeps only FD-FD events.
    global_cfg.enable_proton_fd_sector_filter = (proton_fd_sector_cli != 0);
    global_cfg.proton_fd_sector = (proton_fd_sector_cli != 0 ? proton_fd_sector_cli : 1);

    // Proton CD sector study. This automatically keeps only CD-FD and CD-FT events.
    global_cfg.enable_proton_cd_sector_filter = (proton_cd_sector_cli != 0);
    global_cfg.proton_cd_sector = (proton_cd_sector_cli != 0 ? proton_cd_sector_cli : 1);

    // Optional FD-photon sector selection.  For the controlled topology study
    // use this together with --topology CD-FD so the proton remains CD while
    // the photon is restricted to one FD sector.
    if (photon_sector_cli != 0 && topology_cli != "CD-FD" && topology_cli != "FD-FD") {
        std::cerr << "[main] FATAL: --photon-sector requires an FD-photon topology (CD-FD or FD-FD).\n";
        return 1;
    }
    global_cfg.enable_photon_fd_sector_filter = (photon_sector_cli != 0);
    global_cfg.photon_fd_sector = (photon_sector_cli != 0 ? photon_sector_cli : 1);

    // Auxiliary fiducial cuts. Enable this single switch to apply the additional
    // FD-sector separation, particle-angle, and FT-photon momentum cuts
    // analysis-wide. The numerical values are defined in GlobalCutConfig.
    global_cfg.enable_auxiliary_fiducial_cuts = true;

    // Existing optional propagator/ycol mirror cut.
    global_cfg.enable_dvcsgen_ycol_cut = false;
    global_cfg.dvcsgen_ycol_cut = 0.005;

    set_default_global_cuts(global_cfg);

    std::cout << "[main] Topology selection: "
              << (topology_cli.empty() ? "inclusive" : topology_cli)
              << (photon_sector_cli ? (" photon-S" + std::to_string(photon_sector_cli)) : "")
              << (electron_sector_cli ? (" electron-S" + std::to_string(electron_sector_cli)) : "")
              << (proton_fd_sector_cli ? (" proton-FD-S" + std::to_string(proton_fd_sector_cli)) : "")
              << (proton_cd_sector_cli ? (" proton-CD-S" + std::to_string(proton_cd_sector_cli)) : "")
              << std::endl;
    std::cout << "[main] Sp18 Out detector-quality sector exclusions: "
              << (global_cfg.enable_sp18_out_sector_quality_cuts
                      ? "ON (nominal)"
                      : "OFF (diagnostic)")
              << std::endl;
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

    if (topology_acceptance_study) {
        if (topology_cli.empty() || topology_study_csv.empty()) {
            std::cerr << "[main] FATAL: --topology-acceptance-study requires --topology and --topology-study-csv.\n";
            return 1;
        }
        // Re-derive the cuts for this topology so the reconstructed-MC sample
        // exactly matches the topology extraction.  Do NOT initialize/overwrite
        // the production CSV in this diagnostic mode.
        if (!run_python_exclusivity_analysis(exclusivity_opts)) {
            std::cerr << "[main] FATAL: topology-specific exclusivity optimization failed.\n";
            return 1;
        }
        std::cout << "[main] Topology-acceptance study: topology-specific cuts rederived; production CSV untouched.\n";
    } else if (!acceptance_reweighting_only && !eppi0_normalization_only && !eppi0_production_test) {
        if (!run_python_exclusivity_analysis(exclusivity_opts)) {
            std::cerr << "[main] FATAL: Python exclusivity optimization failed.\n";
            return 1;
        }

        std::cout << "[main] Python exclusivity-cut stage finished. "
                  << "Using nominal cuts from output/jsons/combined_cuts.json.\n";

        initialize_pass2_csv("imports/all_bin_v3.csv",
                             "output/csvs/dvcs_pass2_analysis.csv");
    } else if (acceptance_reweighting_only) {
        std::cout << "[main] Acceptance-reweighting-only mode: reusing existing "
                  << "combined cuts, current calibration, and analysis CSV.\n";
    } else {
        std::cout << "[main] eppi0-normalization-only mode: reusing existing "
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

    if (topology_acceptance_study) {
        const std::string topo_tag = topology_cli +
            (photon_sector_cli ? ("-S" + std::to_string(photon_sector_cli)) : "");
        const std::string outdir = "output/topology_acceptance_model/" + topo_tag;
        if (!run_topology_acceptance_model_study(
                topology_study_csv, genMcTrees, recMcTrees, topo_tag, outdir,
                "output/jsons/combined_cuts.json", 8)) {
            std::cerr << "[main] FATAL: topology acceptance-model study failed.\n";
            return 1;
        }
        std::cout << "[main] Topology acceptance-model study finished. No production CSV was modified.\n";
        return 0;
    }

    if (eppi0_normalization_only || eppi0_production_test) {
        Eppi0NormalizationOptions norm_opts;
        norm_opts.charge_csv_path = "imports/integrated_luminosity/global.csv";
        norm_opts.combined_cuts_json = "output/jsons/combined_cuts.json";
        norm_opts.normalization_json_path = "imports/eppi0_aao_normalization_inputs.json";
        norm_opts.output_dir = "output/data_mc_normalization";
        norm_opts.override_to_unity = false;
        norm_opts.clean_output_dir = true;
        norm_opts.current_response_model_json =
            "output/dvcs_current_dependence/calibration/current_response_model.json";
        norm_opts.write_normalized_yields = false;
        norm_opts.validate_production_fits = eppi0_production_test;
        norm_opts.write_summary_csv = true;
        norm_opts.summary_csv_path =
            "output/data_mc_normalization/eppi0_normalization_summary.csv";
        norm_opts.write_accepted_mc_population = true;
        norm_opts.accepted_mc_population_csv_path =
            "output/data_mc_normalization/accepted_aao_population.csv";
        norm_opts.max_workers = 7;

        std::cout << "[main] Running "
                  << (eppi0_production_test ? "eppi0 PRODUCTION FIT TEST" : "eppi0 normalization DIAGNOSTIC")
                  << " only. Krishna/Neupane proton-efficiency weights are not used by this stage, "
                  << "and production normalized yields will not be overwritten.\n";

        if (!update_eppi0_normalization_csv(
                "output/csvs/dvcs_pass2_analysis.csv",
                dataTrees, eppi0DataTrees, eppi0RecMcTrees, recMcTrees, norm_opts)) {
            std::cerr << "[main] FATAL: eppi0 normalization diagnostic failed.\n";
            return 1;
        }

        std::cout << "[main] "
                  << (eppi0_production_test ? "eppi0 production fit test" : "eppi0 normalization diagnostic")
                  << " finished.\n";
        return 0;
    }

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
    // In topology mode the same calibration is rederived on the selected
    // topology. current_dependence.cpp is region-aware: FD-photon topologies
    // require only S1--S6, while CD-FT requires only FT.
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
        // For single-FD-photon-sector extractions the final event-level current
        // calibration is inherited verbatim from the parent CD-FD extraction
        // below.  Do not attempt to promote a sector-only Sp18-Out theta_e fit
        // into the temporary diagnostic model produced in this block.
        const bool any_single_particle_sector_run =
            global_cfg.enable_photon_fd_sector_filter ||
            global_cfg.enable_electron_fd_sector_filter ||
            global_cfg.enable_proton_fd_sector_filter ||
            global_cfg.enable_proton_cd_sector_filter;
        const bool single_fd_photon_sector_run =
            (global_cfg.enable_photon_fd_sector_filter &&
             global_cfg.enable_topology_filter &&
             global_cfg.required_detector2 == 1);
        current_opts.use_sp18_out_e_theta_response_model = !any_single_particle_sector_run;
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

        // -----------------------------------------------------------------
        // Controlled FD-photon-sector studies inherit ONE current-response
        // calibration from the parent CD-FD extraction.  Current efficiency is
        // a detector calibration; it must not be silently re-fit after S1--S6
        // is selected.  In particular, this preserves exactly the parent
        // Sp18-Out regional intercepts + pooled centered theta_e gradient.
        //
        // Workflow:
        //   1) run --topology CD-FD (no --photon-sector) once; this snapshots
        //      its complete event-level model;
        //   2) every CD-FD-Si run may make its own diagnostics above, but before
        //      any corrected event yields are accumulated we restore the exact
        //      parent model.
        //
        // The snapshot lives outside output/dvcs_current_dependence because
        // that directory is intentionally cleaned by each calibration run.
        // -----------------------------------------------------------------
        const std::filesystem::path current_model =
            "output/dvcs_current_dependence/calibration/current_response_model.json";
        const std::filesystem::path parent_dir =
            "output/topology_sector_calibration";
        const std::filesystem::path parent_model =
            parent_dir / "CD-FD_parent_current_response_model.json";

        const bool parent_cd_fd_run =
            (global_cfg.enable_topology_filter &&
             global_cfg.required_detector1 == 2 &&
             global_cfg.required_detector2 == 1 &&
             !global_cfg.enable_photon_fd_sector_filter);
        const bool sector_cd_fd_run =
            (global_cfg.enable_topology_filter &&
             global_cfg.required_detector1 == 2 &&
             global_cfg.required_detector2 == 1 &&
             global_cfg.enable_photon_fd_sector_filter);

        try {
            if (parent_cd_fd_run) {
                std::filesystem::create_directories(parent_dir);
                std::filesystem::copy_file(
                    current_model, parent_model,
                    std::filesystem::copy_options::overwrite_existing);
                std::cout << "[main] Saved parent CD-FD current-response calibration: "
                          << parent_model.string() << "\n";
            } else if (sector_cd_fd_run) {
                if (!std::filesystem::exists(parent_model)) {
                    std::cerr
                        << "[main] FATAL: photon-sector extraction requires the parent "
                        << "CD-FD current-response calibration. Run --topology CD-FD "
                        << "without --photon-sector first. Missing: "
                        << parent_model.string() << "\n";
                    std::exit(EXIT_FAILURE);
                }
                std::filesystem::copy_file(
                    parent_model, current_model,
                    std::filesystem::copy_options::overwrite_existing);
                std::cout
                    << "[main] Photon-sector S" << global_cfg.photon_fd_sector
                    << ": restored exact parent CD-FD current-response calibration "
                    << "before event-level yield accumulation.\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "[main] FATAL: failed to snapshot/restore parent CD-FD "
                      << "current-response calibration: " << e.what() << "\n";
            std::exit(EXIT_FAILURE);
        }

        // Generic single-particle sector diagnostics (electron, proton-FD,
        // proton-CD) inherit the exact inclusive parent current calibration.
        // This prevents the detector-sector selection from redefining the
        // current-efficiency correction being tested.  Run one inclusive
        // parent extraction first to create the protected snapshot.
        const std::filesystem::path inclusive_parent_model =
            parent_dir / "inclusive_parent_current_response_model.json";
        const bool no_particle_sector_filter =
            !global_cfg.enable_photon_fd_sector_filter &&
            !global_cfg.enable_electron_fd_sector_filter &&
            !global_cfg.enable_proton_fd_sector_filter &&
            !global_cfg.enable_proton_cd_sector_filter;
        const bool inclusive_parent_run =
            !global_cfg.enable_topology_filter && no_particle_sector_filter;
        const bool generic_particle_sector_run =
            global_cfg.enable_electron_fd_sector_filter ||
            global_cfg.enable_proton_fd_sector_filter ||
            global_cfg.enable_proton_cd_sector_filter;
        // Topology-resolved extractions are detector-population diagnostics, not
        // independent current-efficiency calibrations.  They must inherit the
        // exact same event-level current-response model as the inclusive parent.
        // This also makes FD-FD/CD-FD/CD-FT directly comparable and avoids a
        // topology selection silently changing the calibration being tested.
        const bool topology_diagnostic_run = global_cfg.enable_topology_filter;
        try {
            if (inclusive_parent_run) {
                std::filesystem::create_directories(parent_dir);
                std::filesystem::copy_file(
                    current_model, inclusive_parent_model,
                    std::filesystem::copy_options::overwrite_existing);
                std::cout << "[main] Saved inclusive parent current-response calibration: "
                          << inclusive_parent_model.string() << "\n";
            } else if (generic_particle_sector_run || topology_diagnostic_run) {
                if (!std::filesystem::exists(inclusive_parent_model)) {
                    std::cerr << "[main] FATAL: detector-subset extraction requires the inclusive "
                              << "parent current-response calibration. Run the same command once "
                              << "without --topology or a sector option first. Missing: "
                              << inclusive_parent_model.string() << "\n";
                    std::exit(EXIT_FAILURE);
                }
                std::filesystem::copy_file(
                    inclusive_parent_model, current_model,
                    std::filesystem::copy_options::overwrite_existing);
                std::cout << "[main] Restored exact inclusive parent current-response calibration "
                          << (topology_diagnostic_run ? "for topology extraction.\n"
                                                     : "for particle-sector extraction.\n");
            }
        } catch (const std::exception& e) {
            std::cerr << "[main] FATAL: failed to snapshot/restore inclusive parent "
                      << "current-response calibration: " << e.what() << "\n";
            std::exit(EXIT_FAILURE);
        }
    }



    // --------- pass-1-style DVpi0P efficiency-map derivation ----------
    // Derive the proton-theta DATA/AAOgen map in six FD sectors + CD, then
    // reweight reconstructed AAOgen by that theta map and fit the residual
    // DATA/MC dependence versus proton momentum.  total_counts consumes the
    // resulting sequential theta*p map for reconstructed DVCS MC only.
    // Experimental DVCS DATA and generated DVCS MC are never weighted by it.
    if (use_eppi0_production_normalization) {
        const std::string csv_main = "output/csvs/dvcs_pass2_analysis.csv";

        Eppi0NormalizationOptions norm_opts;
        norm_opts.charge_csv_path = "imports/integrated_luminosity/global.csv";
        norm_opts.combined_cuts_json = "output/jsons/combined_cuts.json";
        norm_opts.normalization_json_path = "imports/eppi0_aao_normalization_inputs.json";
        norm_opts.current_response_model_json =
            "output/dvcs_current_dependence/calibration/current_response_model.json";
        norm_opts.output_dir = "output/data_mc_normalization";
        norm_opts.override_to_unity = false;
        norm_opts.write_normalized_yields = false;
        norm_opts.write_summary_csv = true;
        norm_opts.write_accepted_mc_population = true;
        norm_opts.max_workers = 7;

        if (!update_eppi0_normalization_csv(csv_main,
                                            dataTrees,
                                            eppi0DataTrees,
                                            eppi0RecMcTrees,
                                            recMcTrees,
                                            norm_opts)) {
            std::cerr << "[main] FATAL: pass-1-style eppi0 normalization failed.\n";
            std::exit(EXIT_FAILURE);
        }
        std::cout << "[main] Completed eppi0 DATA/AAOgen diagnostic; "
                  << "sequential theta and residual-p efficiency maps are ready for reconstructed DVCS MC.\n";
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
        // Mutually exclusive detector-normalization modes.  Production default:
        // eppi0 normalization OFF and Krishna/Neupane proton-efficiency correction
        // ON.  --eppi0-normalization is an explicit diagnostic alternative that
        // applies the sequential eppi0 theta*p efficiency map to reconstructed
        // DVCS MC and disables the Krishna/Neupane correction.
        total_count_opts.apply_neupane_proton_efficiency_correction =
            !use_eppi0_production_normalization;
        total_count_opts.apply_eppi0_efficiency_to_dvcs_rec_mc =
            use_eppi0_production_normalization;
        total_count_opts.eppi0_efficiency_summary_csv =
            "output/data_mc_normalization/eppi0_normalization_summary.csv";
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
    if (global_cfg.enable_topology_filter) {
        // A topology-restricted cross-section run already fixes the photon/proton
        // detector configuration globally.  The BSA module performs its own
        // independent FD-vs-FT photon subdivision and helicity diagnostics, so
        // running it here would create empty forbidden-topology samples and is
        // unrelated to the topology-resolved unpolarized cross-section study.
        std::cout << "[main] Topology mode: skipping BSA stage (including the "
                  << "BSA FD-vs-FT and helicity-scrambling diagnostics).\n";
    } else {
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
        bsa_opts.pi0_leakage_relative_uncertainty = 0.10;
        bsa_opts.make_plots = true;
        bsa_opts.make_photon_topology_study = true;
        bsa_opts.make_helicity_scrambling_study = true;
        bsa_opts.helicity_scramble_replicas = 100;
        bsa_opts.helicity_scramble_seed = 20260915ULL;
        bsa_opts.max_workers = 7;

        if (!update_bsa_counts_csv(dataTrees, eppi0DataTrees, bsa_opts)) {
            std::cerr << "[main] ERROR: update_bsa_counts_csv failed.\n";
            std::exit(EXIT_FAILURE);
        }

        // Validation only: CEBAF's rapid helicity reversal should make the
        // integrated + and - charges equal to high precision. Sp18 is expected
        // to be unavailable because helicity-resolved Faraday-cup charge was
        // not recoverable for that period.
        (void)write_bsa_helicity_charge_balance();
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

    // --------- Bin-centering-correction analysis-note outputs ----------
    // Production Fbin values are inherited directly from the pass-1 analysis
    // (imports/all_bin_v3.csv), exactly like Frad.  This step is read-only:
    // it produces diagnostics for the reused pass-1 factors and does not
    // modify the pass-2 analysis CSV.
    {
        if (!write_bin_centering_analysis_note_outputs(
                "imports/all_bin_v3.csv",
                "output/bin_centering_corrections")) {
            std::cerr
                << "[main] ERROR: bin-centering-correction analysis-note outputs failed.\n";
            return 1;
        }
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
        // Keep every cut variation on the same detector-normalization branch as
        // the nominal run.  In eppi0 mode the runner re-derives a fresh theta*p
        // map under each varied selection; in --no-eppi0-normalization mode it
        // instead restores Krishna for every variation.
        cut_variation_opts.use_eppi0_production_normalization =
            use_eppi0_production_normalization;
        cut_variation_opts.eppi0_charge_csv_path =
            "imports/integrated_luminosity/global.csv";
        cut_variation_opts.eppi0_normalization_json_path =
            "imports/eppi0_aao_normalization_inputs.json";
        cut_variation_opts.current_response_model_json =
            "output/dvcs_current_dependence/calibration/current_response_model.json";

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

    // Hall-A / CLAS6 / world-data comparisons are intentionally not run here.
    // They are maintained as separate external comparison workflows and were
    // removed from the nominal chain to avoid duplicated runtime.

    std::cout << "All done." << std::endl;
    return 0;
}