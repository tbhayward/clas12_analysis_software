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
    //   The Sp19 Inb current-efficiency factors written to the CSV are copied
    //   from Fa18 Inb instead of being taken from the Sp19 Inb luminosity scan.
    //
    // Motivation:
    //   The only low-current Sp19 Inb scan point is run 6616 at 5 nA, whose
    //   Faraday Cup charge total is currently suspect. Therefore the Sp19 Inb
    //   fitted slope is not used by default.
    //
    // This replacement is applied to:
    //   current efficiency factor, ep->epg,   exp, Sp19 Inb
    //   current efficiency factor, ep->epg,   mc,  Sp19 Inb
    //   current efficiency factor, ep->eppi0, exp, Sp19 Inb
    //   current efficiency factor, ep->eppi0, mc,  Sp19 Inb
    //
    // The raw Sp19 Inb scan is still processed and plotted for diagnostic
    // purposes, but the saved CSV values used downstream are Fa18 Inb values.
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
    bool use_e_theta_linear_data_current_efficiency = true;

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
 *     weighted_data_rel = event-weighted fitted data response at the actual
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
 * factors written to the CSV are copied from Fa18 Inb. The raw Sp19 Inb fit is
 * still produced in the diagnostic output.
 *
 * It also applies the saved MC current-efficiency factors to the reconstructed
 * MC yield columns and writes the current-corrected reconstructed MC columns
 * needed by downstream eppi0 normalization, pi0-contamination, acceptance, and
 * unfolding modules:
 *
 *   reconstructed current corrected yield, ep->epg, mc, <period>
 *   reconstructed current corrected yield, ep->epg, <topology>, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0, <topology>, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0->epg, mc, <period>
 *   reconstructed current corrected yield, ep->eppi0->epg, <topology>, mc, <period>
 *
 * The correction applied is N_rec,current-corrected = N_rec / f_current^MC.
 *
 * If options.override_to_unity is true, no ROOT loops are performed, all
 * current-efficiency factors are written as (1,0), and the current-corrected MC
 * reconstructed-yield columns are copied from the uncorrected reconstructed-yield
 * columns.
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