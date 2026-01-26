#ifndef GLOBAL_CUTS_H
#define GLOBAL_CUTS_H

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// GlobalCutConfig
//
// This struct defines analysis-wide cuts that all stages should share.
// It is intentionally "fail-fast":
//   - If enable_dvcsgen_ycol_cut is true, callers must use the overloads that
//     provide the required kinematics.
//   - If enable_topology_filter is true, callers must use the overloads that
//     provide detector1 and detector2.
// -----------------------------------------------------------------------------
struct GlobalCutConfig {
    // Baseline global DVCS exclusivity-style cuts
    double t1_abs_max         = 1.0;  // cut is (-t1) < t1_abs_max
    double open_angle_min_deg = 5.0;  // open_angle_ep2_deg > open_angle_min_deg
    double pTmiss_max         = 0.2;  // pTmiss <= pTmiss_max

    // Optional dvcsgen P2_pos (ycol / propagator) cut
    bool   enable_dvcsgen_ycol_cut = false;
    double dvcsgen_ycol_cut        = 0.005;

    // Optional topology filter (detector1/detector2)
    //
    // detector1 is proton region (expected: 1 FD or 2 CD)
    // detector2 is photon region (expected: 0 FT or 1 FD)
    //
    // If enable_topology_filter is true, ALL callers must use the overloads
    // that provide detector1 and detector2; otherwise we throw FATAL.
    bool enable_topology_filter = true;

    // Required detector assignments when topology filtering is enabled.
    // Must be in {0,1,2} when enable_topology_filter == true.
    // 0 = FT, 1 = FD, 2 = CD
    int required_detector1 = 1;  // proton detector (ignored unless enabled)
    int required_detector2 = 1;  // photon detector (ignored unless enabled)

    // Run blacklist (global)
    std::vector<int> excluded_runs;
};

// Default config instance (process-wide)
const GlobalCutConfig& default_global_cuts();

// dvcsgen ycol helper (exposed for diagnostics/plotting)
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

// -----------------------------------------------------------------------------
// passes_global_cuts overloads
//
// NOTE:
//  - If cfg.enable_dvcsgen_ycol_cut is true, you must call an overload that
//    provides the needed kinematics (or we throw FATAL).
//  - If cfg.enable_topology_filter is true, you must call an overload that
//    provides detector1 and detector2 (or we throw FATAL).
// -----------------------------------------------------------------------------

// Legacy overload: no kinematics, no topology.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg);

// Legacy overload: with ycol kinematics, but no topology.
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
                        const GlobalCutConfig& cfg);

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
                        const GlobalCutConfig& cfg);

// NEW: topology-aware overloads (detector1/detector2), no ycol kinematics.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        const GlobalCutConfig& cfg);

// NEW: topology-aware overloads with ycol kinematics.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        double Ebeam,
                        double e_p,
                        double e_theta,
                        double e_phi,
                        double p2_p,
                        double p2_theta,
                        double p2_phi,
                        const GlobalCutConfig& cfg);

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        const std::string& period_label,
                        double e_p,
                        double e_theta,
                        double e_phi,
                        double p2_p,
                        double p2_theta,
                        double p2_phi,
                        const GlobalCutConfig& cfg);

// Run exclusion helper
bool is_excluded_run(int runnum, const GlobalCutConfig& cfg = default_global_cuts());

// ROOT TCut helper (only valid when ycol and topology filters are disabled)
std::string global_cuts_tcut(const GlobalCutConfig& cfg = default_global_cuts());

// JSON writer for config echoing
void write_global_cuts_config_json(const std::string& out_json_dir,
                                  const GlobalCutConfig& cfg = default_global_cuts());

#endif // GLOBAL_CUTS_H