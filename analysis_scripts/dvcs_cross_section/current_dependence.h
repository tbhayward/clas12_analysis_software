#ifndef CURRENT_DEPENDENCE_H
#define CURRENT_DEPENDENCE_H

#include <map>
#include <string>

class TTree;

struct CurrentDependenceOptions {
    std::string charge_csv_path = "imports/integrated_luminosity/global.csv";
    std::string combined_cuts_json = "output/jsons/combined_cuts.json";
    std::string output_dir = "output/dvcs_current_dependence";

    // Existing override:
    //   true  -> write all current-efficiency factors as (1,0)
    //   false -> compute the current-efficiency factors normally
    bool override_to_unity = false;

    // Charge-column convention for unpolarized data normalization.
    //
    // If true:
    //   Use column 2 of the integrated-luminosity CSV for all periods.
    //
    // If false:
    //   Preserve the older mixed convention:
    //     Spring 2018 -> column 2
    //     Fall 2018 and Spring 2019 -> column 3 + column 4
    //
    // Spring 2018 is always forced to column 2 internally.
    // This option is ignored for Fa18/Sp19 if the scaled columns-3-to-5 mode below is enabled.
    bool use_second_column_charge_for_all_unpolarized = true;

    // Optional Fa18/Sp19 charge-column mode:
    //
    //   true:
    //     For Fa18 Inb, Fa18 Out, and Sp19 Inb only, use
    //
    //       charge_unpol = columns_3_to_5_charge_sum_scale
    //                      * (column 3 + column 4 + column 5)
    //
    //     for unpolarized data charge normalization.
    //
    //     Spring 2018 still uses column 2.
    //
    //   false:
    //     fall back to use_second_column_charge_for_all_unpolarized / legacy mode.
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized = false;

    // Scale applied in the optional Fa18/Sp19 columns-3-to-5 mode.
    double columns_3_to_5_charge_sum_scale = 1.025;

    // New default behavior:
    //
    // If true:
    //   The Sp19 Inb DATA current response is copied from Fa18 Inb instead of
    //   being taken from the Sp19 Inb luminosity scan. The independently
    //   determined Sp19 MC response is retained.
    //
    // Motivation:
    //   The only low-current Sp19 Inb scan point is run 6616 at 5 nA, whose
    //   Faraday Cup charge total is currently suspect. Therefore the Sp19 Inb
    //   fitted slope is not used by default.
    //
    // This DATA replacement is applied to ep->epg and ep->eppi0, including
    // the regional FT/S1--S6 response used by the event-level model. The raw
    // Sp19 Inb DATA scan is still processed and plotted diagnostically.
    bool use_fa18_inb_current_efficiency_for_sp19_inb = true;

    // Optional DVCS MC acceptance override.
    //
    // If true:
    //   The ep->epg MC current-efficiency factors written to the CSV are forced
    //   to unity because the no-background dvcsgen files are being used for
    //   ep->epg generated/reconstructed MC counts upstream in total_counts.cpp.
    //
    // This affects only ep->epg MC. DATA factors and the ep->eppi0 /
    // ep->eppi0->epg treatment are unchanged. The scan is still processed and
    // plotted diagnostically.
    bool use_nobkg_dvcs_mc_counts = false;

    // If true, DATA raw-yield corrections use a period-dependent linear
    // current-efficiency factor in electron polar angle,
    //
    //     f_current(theta_e) = m theta_e + b,
    //
    // obtained from the same kinematic current-efficiency diagnostic machinery.
    // This is applied to both ep->epg and ep->eppi0 DATA normalized raw yields.
    // The CSV current-efficiency-factor columns still store the integrated
    // factors for bookkeeping and diagnostics; only the DATA yield correction
    // uses the theta_e-dependent factor. If a row has no finite theta_e value
    // or a valid fit cannot be built, the code falls back to the integrated
    // current-efficiency factor for that row.
    bool use_e_theta_linear_data_current_efficiency = false;

    // Photon-region current-dependence diagnostic.
    //
    // If true, the ordinary ep->epgamma current-dependence study is repeated
    // after splitting the photon into seven mutually exclusive regions: FT and
    // FD sectors 1--6. The regional calibration is saved outside the analysis
    // CSV and is used to build the event-level response model consumed by
    // total_counts.cpp. Diagnostic plots/tables are written under
    //
    //   output/dvcs_current_dependence/epg/sector_dependence_diagnostic/
    //
    // and can be inspected independently of the production CSV columns.
    bool enable_photon_region_current_diagnostic = true;

    // Repeat the same FT/S1--S6 diagnostic for direct ep->eppi0 data.
    bool enable_eppi0_photon_region_current_diagnostic = true;

    // Within each neutral-particle region, repeat the current-efficiency study
    // versus electron and neutral-particle polar angle. This is diagnostic only
    // and is used to decide whether a residual theta dependence remains after
    // regionization.
    bool enable_region_theta_current_diagnostic = true;

    // Finalized production mode.
    //
    // true:
    //   - retain the regional FT/S1--S6 calibration needed by production,
    //   - retain the DVCS polar-angle diagnostics needed to justify the
    //     Sp18-Out theta_e term,
    //   - skip exploratory scans that do not feed the final response model,
    //   - prune output/dvcs_current_dependence after a successful run to the
    //     analysis-note figures/tables plus the calibration JSON.
    //
    // false restores the broader exploratory diagnostic output.
    bool finalized_production_mode = true;

    // The original inclusive seven-variable (Q2,xB,t,phi,theta_e,theta_p,
    // theta_gamma) current-efficiency diagnostic predates the final regional
    // response model. It requires additional complete DATA/MC tree passes and
    // is not needed by nominal production.
    bool enable_exploratory_kinematic_current_diagnostic = false;

    // The direct-pi0 angular model-selection scan was used to establish that no
    // angular extension is required for the final eppi0 response. It is kept
    // available for dedicated validation jobs, but is not rerun in production.
    bool enable_eppi0_region_theta_current_diagnostic = false;

    // Remove stale current-dependence output before starting a nominal run.
    // This prevents plots from an older exploratory configuration from being
    // mistaken for products of the current calibration.
    bool clean_output_dir_before_run = true;

    // Production DATA response prescription.  When enabled, only Sp18 Out
    // ep->epgamma uses the region-conditioned common linear electron-angle
    // response established by the diagnostics:
    //
    //   s_r(theta_e) = s_r^(0) + a*(theta_e - theta_bar_r).
    //
    // All other DATA periods/channels remain regional-only.  This is an
    // explicit analysis choice; the code does not automatically promote the
    // variable with the best BIC.
    bool use_sp18_out_e_theta_response_model = true;

    // Diagnostic proxy for relative photon-efficiency modeling.  It compares
    // CD_FT to CD_FD using the double ratio
    //
    //   (DATA/MC)_CD,FT / (DATA/MC)_CD,FD
    //
    // after matching electron and proton polar-angle phase space and applying
    // the same current-response corrections used in production.  Because the
    // proton is CD on both sides, charged-particle acceptance effects largely
    // cancel and the remaining ratio primarily tests relative FT/FD photon
    // reconstruction modeling.
    bool enable_relative_ft_fd_photon_efficiency_diagnostic = false;

    // Calibration product consumed by total_counts.cpp. The file contains the
    // regional DATA current-response slopes, regional reconstructed-MC factors,
    // and run->current lookup used for event-level current correction.
    std::string response_model_json = "output/dvcs_current_dependence/calibration/current_response_model.json";

    // When false, this stage determines/writes the calibration and diagnostics
    // but does not divide already-binned yields. Production now applies the
    // response event-by-event inside total_counts.cpp.
    bool apply_legacy_binned_current_corrections = false;

    // Current correction convention for the misidentified pi0 background MC
    // sample ep->eppi0->epg.
    //
    // true  -> use the ep->epg MC current-efficiency factor. This is the
    //          default because the reconstructed final state being counted is
    //          epgamma, even though the generated sample is ep->eppi0.
    // false -> use the ep->eppi0 MC current-efficiency factor, which was the
    //          previous behavior.
    bool use_epg_mc_current_factor_for_eppi0_bkg = true;

    int max_workers = 7;
};

/**
 * update_current_dependence_factors_csv
 *
 * Performs the current-dependence study adapted from dvcs_current_dependence.py.
 *
 * For each channel and period, it determines:
 *
 *   DATA:
 *     weighted_data_rel = charge-weighted fitted data response at the actual
 *                         current mixture divided by fitted zero-current response
 *
 *   MC:
 *     For ep->epg, mc_ref_rel is the fitted MC efficiency at the reference
 *     current divided by fitted zero-current MC efficiency.
 *
 *     For ep->eppi0, no non-production-current MC exists, so the MC current
 *     factor is derived from the data factor scaled by the DVCS MC/data ratio:
 *
 *         eppi0_mc_factor = eppi0_data_factor * (dvcs_mc_factor / dvcs_data_factor)
 *
 * and writes:
 *
 *   current efficiency factor, ep->epg,   exp, <period> = (weighted_data_rel,stat)
 *   current efficiency factor, ep->epg,   mc,  <period> = (mc_ref_rel,stat)
 *   current efficiency factor, ep->eppi0, exp, <period> = (weighted_data_rel,stat)
 *   current efficiency factor, ep->eppi0, mc,  <period> = (mc_ref_rel,stat)
 *
 * If options.use_fa18_inb_current_efficiency_for_sp19_inb is true, the Sp19 Inb
 * DATA response is copied from Fa18 Inb while the Sp19 MC response is retained.
 * The raw Sp19 DATA fit is still produced in the diagnostic output.
 *
 * In the nominal workflow this stage writes a regional response-model JSON and
 * does NOT modify already-binned yields. total_counts.cpp consumes that model
 * and fills the normalized DATA and reconstructed-current-corrected MC columns
 * event-by-event, while preserving the original unit-weight raw-yield columns.
 *
 * options.apply_legacy_binned_current_corrections=true restores the previous
 * post-binning scalar correction behavior for compatibility tests only.
 */
bool update_current_dependence_factors_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& dvcsGenMcTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const CurrentDependenceOptions& options);

#endif // CURRENT_DEPENDENCE_H