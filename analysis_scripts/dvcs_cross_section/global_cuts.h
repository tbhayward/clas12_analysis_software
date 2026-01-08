#ifndef GLOBAL_CUTS_H
#define GLOBAL_CUTS_H

#include <string>
#include <vector>

// global_cuts.h
/*
 * Global, universal event-level cuts for DVCS analysis.
 *
 * Base cuts applied everywhere:
 *   (-t1) < t1_abs_max
 *   open_angle_ep2 > open_angle_min_deg (deg)
 *   pTmiss <= pTmiss_max (GeV)
 *
 * Optional: dvcsgen --ycol cut mirror.
 */

struct GlobalCutConfig {
    double t1_abs_max = 1.0;
    double open_angle_min_deg = 5.0;
    double pTmiss_max = 0.20; // (GeV)

    // dvcsgen ycol cut mirror (ON by default here, per your current config)
    bool   enable_dvcsgen_ycol_cut = true;
    double dvcsgen_ycol_cut = 0.005; // dimensionless threshold

    // Runs to exclude globally (data quality blacklist).
    std::vector<int> excluded_runs = {
        3867, 5046, 5047, 5051, 5128, 5129, 5130, 5160, 5163, 5165, 5166, 5167, 5168, 5169,
        5180, 5181, 5182, 5183, 5247, 5448, 5495, 5496, 5615, 5567
    };
};

const GlobalCutConfig& default_global_cuts();

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

// Run blacklist helper
bool is_excluded_run(int runnum,
                     const GlobalCutConfig& cfg = default_global_cuts());

// ROOT-style TCut string for base cuts only.
// FAIL-FAST if cfg.enable_dvcsgen_ycol_cut is true (use C++ cut instead).
std::string global_cuts_tcut(const GlobalCutConfig& cfg = default_global_cuts());

// Persist the configuration once for provenance.
void write_global_cuts_config_json(const std::string& out_json_dir,
                                   const GlobalCutConfig& cfg = default_global_cuts());

// dvcsgen ycol scalar (P1_pos) helpers
double dvcsgen_ycol_value(double Ebeam,
                          double e_p,
                          double e_theta,
                          double e_phi,
                          double p2_p,
                          double p2_theta,
                          double p2_phi,
                          const GlobalCutConfig& cfg);

double dvcsgen_ycol_value(const std::string& period_label,
                          double e_p,
                          double e_theta,
                          double e_phi,
                          double p2_p,
                          double p2_theta,
                          double p2_phi,
                          const GlobalCutConfig& cfg);

// Convenience wrappers (use default_global_cuts())
double dvcsgen_ycol_value(double Ebeam,
                          double e_p,
                          double e_theta,
                          double e_phi,
                          double p2_p,
                          double p2_theta,
                          double p2_phi);

double dvcsgen_ycol_value(const std::string& period_label,
                          double e_p,
                          double e_theta,
                          double e_phi,
                          double p2_p,
                          double p2_theta,
                          double p2_phi);

#endif // GLOBAL_CUTS_H