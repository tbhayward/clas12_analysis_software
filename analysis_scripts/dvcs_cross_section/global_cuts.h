#ifndef GLOBAL_CUTS_H
#define GLOBAL_CUTS_H

#include <string>
#include <vector>

/*
 * Global, universal event-level cuts for DVCS analysis.
 *
 * Base cuts applied everywhere:
 *   (-t1) < t1_abs_max
 *   open_angle_ep2 > open_angle_min_deg (deg)
 *   pTmiss <= pTmiss_max (GeV)
 *
 * Optional (OFF by default): dvcsgen --ycol cut mirror.
 *
 * Empirically (from your diagnostic), dvcsgen generated with:
 *   --ycol 0.005
 * enforces:
 *   P2_pos > 0.005
 * where:
 *   Q2_calc = - (k - kprime)^2
 *   P2_pos  = 2 (kprime . qgamma) / Q2_calc
 *
 * This requires the beam energy and the electron/photon kinematics.
 */

struct GlobalCutConfig {
    double t1_abs_max = 1.0;
    double open_angle_min_deg = 5.0;
    double pTmiss_max = 0.20; // (GeV)

    // dvcsgen ycol cut mirror (OFF by default)
    bool   enable_dvcsgen_ycol_cut = false;
    double dvcsgen_ycol_cut = 0.005; // dimensionless threshold on P2_pos

    // Runs to exclude globally (data quality blacklist).
    std::vector<int> excluded_runs = {
        3867, 5046, 5047, 5051, 5128, 5129, 5130, 5160, 5163, 5165, 5166, 5167, 5168, 5169,
        5180, 5181, 5182, 5183, 5247, 5448, 5495, 5496, 5615, 5567
    };
};

// -----------------------------------------------------------------------------
// P2 enforcement API
// -----------------------------------------------------------------------------

// Load the global cuts JSON. MUST be called before using passes_global_cuts()
// or is_excluded_run() with default cfg (fail-fast otherwise).
bool load_global_cuts_config(const std::string& json_path);

// True only if load_global_cuts_config() has succeeded.
bool global_cuts_loaded();

// Accessor used as the default cfg provider.
// This should FAIL-FAST if global cuts are not loaded.
const GlobalCutConfig& default_global_cuts();

// Optional: print the active config (useful for provenance).
void print_global_cuts_summary();

// -----------------------------------------------------------------------------
// Cut predicates
// -----------------------------------------------------------------------------

// Base version: only the three universal cuts.
// FAIL-FAST if cfg.enable_dvcsgen_ycol_cut is true (because additional inputs are required).
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg = default_global_cuts());

// Extended version (explicit Ebeam): includes optional dvcsgen ycol cut mirror.
// Angles must be in radians (as in your ROOT branches).
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        double Ebeam,
                        double e_p,
                        double e_theta,
                        double e_phi,
                        double p2_p,
                        double p2_theta,
                        double p2_phi,
                        const GlobalCutConfig& cfg = default_global_cuts());

// Extended version (period label): looks up Ebeam deterministically from the run period label.
// Allowed labels (exact):
//   sp18_inb, sp18_out, fa18_inb, fa18_out, sp19_inb
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const std::string& period_label,
                        double e_p,
                        double e_theta,
                        double e_phi,
                        double p2_p,
                        double p2_theta,
                        double p2_phi,
                        const GlobalCutConfig& cfg = default_global_cuts());

// NEW: check if a run is globally blacklisted.
bool is_excluded_run(int runnum,
                     const GlobalCutConfig& cfg = default_global_cuts());

// Convenience: ROOT-style TCut string for base cuts only.
// FAIL-FAST if cfg.enable_dvcsgen_ycol_cut is true (use the C++ cut instead).
std::string global_cuts_tcut(const GlobalCutConfig& cfg = default_global_cuts());

// Persist the configuration once for provenance.
void write_global_cuts_config_json(const std::string& out_json_dir,
                                  const GlobalCutConfig& cfg = default_global_cuts());

#endif // GLOBAL_CUTS_H