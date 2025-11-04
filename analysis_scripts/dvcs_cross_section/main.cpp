#include "make_dirs.h"
#include "load_trees.h"
#include "periods.h"
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

    // Root of output tree (used by several stages)
    const std::string output_root = "output";

    // --- Load binning scheme ---
    const std::string csv_file_path = "imports/integrated_bin_v2.csv";
    auto binning_scheme = load_binning_scheme(csv_file_path);
    std::cout << "Loaded binning scheme: " << binning_scheme.size() << " bins" << std::endl;

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
    // FIX: pass eppi0BkgTrees in the correct slot to match load_trees.h
    loadTrees(dataTrees, genMcTrees, recMcTrees,
              eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees, eppi0BkgTrees,
              radGenMcTrees, radRecMcTrees);
    std::cout << "All trees loaded successfully." << std::endl;

    // // Run exclusivity cut extraction (single-threaded for stability)
    // runAllExclusivityCuts(
    //     dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
    //     "output/jsons", "output/exclusivity_plots", 1
    // );
    // std::cout << "Exclusivity-cut stage finished." << std::endl;

    // --------- Global bin-averaged kinematics ----------
    std::vector<std::string> dvcs_periods;
    for (const auto& P : CANONICAL_PERIODS()) {
        dvcs_periods.push_back(P.tree_key);
    }
    std::vector<std::string> topologies = {"(FD,FD)", "(CD,FD)", "(CD,FT)"};
    const std::string analysis_type = "dvcs";
    const std::string output_json_means = "output/jsons/bin_means_global.json";

    // calculate_bin_means(dvcs_periods, topologies, analysis_type, binning_scheme, output_json_means, 
    //     dataTrees);

    // // --------- Total counts after exclusivity cuts (by helicity) ----------
    const std::string cuts_json_path   = "output/jsons/combined_cuts.json"; 
    // // produced by exclusivity_cuts

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

    // compute_uncorrected_cross_sections(
    //     binning_scheme,
    //     "output/jsons",                    // bin volume JSON directory
    //     "output/jsons",               // unfolded counts per helicity
    //     "imports/integrated_luminosity",        // luminosity text files
    //     "output/uncorrected_cross_section"      // output dir
    // );

    // compare_unpolarized_cross_sections_sp18out_vs_fa18out(
    //     binning_scheme,
    //     "output/jsons",                    // bin volume JSON directory
    //     "output/jsons",                    // unfolded counts
    //     "imports/integrated_luminosity",   // luminosity text files (use total column)
    //     "output/uncorrected_cross_section/" // output dir
    // ); // #endif

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

    // --------- Models vs bin-centered data (plots only; no heavy recalculation) ----------
    plot_models_vs_bincentered(
        binning_scheme,
        "output/bin_centering/jsons",  // where bin_centered_xsec_<E>.json live
        "output/bin_centering",        // plots will go under output/bin_centering/plots
        ModelPaths(),                  // your configured paths to models
        true,                          // vgg_globalfit flag
        24                          // phi sampling density (0..360 by 1 deg)
    );

    std::cout << "All done." << std::endl;
    return 0;
}