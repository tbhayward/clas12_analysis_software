#ifndef GLOBAL_CUTS_H
#define GLOBAL_CUTS_H

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// GlobalCutConfig
//
// Analysis-wide cuts shared by all event-loop stages.
//
// Detector codes used by the pass2 trees:
//   detector1 = proton detector: 1 FD, 2 CD
//   detector2 = photon detector: 0 FT, 1 FD
//
// Azimuthal branch convention:
//   e_phi, p1_phi, p2_phi are expected in radians.
//   All sector boundaries below are in degrees.
//
// Sector filters are intentionally fail-fast. If any sector filter is enabled,
// callers must use a sector-aware passes_global_cuts(...) overload that provides
// e_phi, p1_phi, and p2_phi. This prevents silently producing an inclusive result
// when the user intended a sector-restricted result.
// -----------------------------------------------------------------------------
struct GlobalCutConfig {
    // Baseline global DVCS exclusivity-style cuts.
    double t1_abs_max         = 1.0;  // cut is (-t1) < t1_abs_max
    double open_angle_min_deg = 5.0;  // open_angle_ep2_deg > open_angle_min_deg
    // The old global pTmiss pre-cut is now disabled by default.
    // pTmiss is handled as a topology/period-dependent quantile exclusivity cut.
    bool   enable_pTmiss_cut = false;
    double pTmiss_max        = 0.2;  // used only if enable_pTmiss_cut=true

    // Optional dvcsgen P2_pos (ycol / propagator) cut.
    bool   enable_dvcsgen_ycol_cut = false;
    double dvcsgen_ycol_cut        = 0.005;

    // Optional auxiliary fiducial cuts. These are stricter analysis-wide
    // particle-angle/topology cuts used for detector-acceptance diagnostics and
    // systematic studies. Angles in this config are degrees; tree branches are
    // still expected in radians.
    bool   enable_auxiliary_fiducial_cuts = false;
    bool   auxiliary_require_distinct_fd_sectors = true;

    double auxiliary_e_theta_min_deg = 8.0;
    double auxiliary_e_theta_max_deg = 25.0;  // applied only in CD-FD and FD-FD

    double auxiliary_fd_proton_theta_min_deg = 8.0;
    double auxiliary_fd_proton_theta_max_deg = 40.0;

    double auxiliary_fd_photon_theta_min_deg = 8.0;  // applied only in FD-FD
    double auxiliary_fd_photon_theta_max_deg = 40.0; // applied to all FD photons

    double auxiliary_cd_proton_theta_min_deg = 40.0;
    double auxiliary_cd_proton_theta_max_deg = 64.23;

    double auxiliary_ft_photon_p_min_GeV = 4.0;

    // Optional single-topology filter.
    // required_detector1/2 set as (p-region, gamma-region):
    //   FD-FD: (1, 1)
    //   CD-FD: (2, 1)
    //   CD-FT: (2, 0)
    bool enable_topology_filter = false;
    int  required_detector1     = 2;
    int  required_detector2     = 0;

    // Optional electron FD-sector filter. Electron is always FD in this analysis.
    // Valid sector: 1..6.
    bool enable_electron_fd_sector_filter = false;
    int  electron_fd_sector               = 1;

    // Optional proton FD-sector filter. Only applies to events where detector1==1.
    // Therefore this naturally selects the FD-FD topology.
    // Valid sector: 1..6.
    bool enable_proton_fd_sector_filter = false;
    int  proton_fd_sector               = 1;

    // Optional proton CD-sector filter. Only applies to events where detector1==2.
    // Therefore this naturally selects CD-FD and CD-FT topologies.
    // Valid sector: 1..3.
    bool enable_proton_cd_sector_filter = false;
    int  proton_cd_sector               = 1;

    // Optional photon FD-sector filter. Only applies to events where detector2==1.
    // Therefore this naturally selects CD-FD and FD-FD topologies.
    // Valid sector: 1..6.
    bool enable_photon_fd_sector_filter = false;
    int  photon_fd_sector               = 1;

    // Run blacklist.
    std::vector<int> excluded_runs;
};

// Process-wide default config.
const GlobalCutConfig& default_global_cuts();
void set_default_global_cuts(const GlobalCutConfig& cfg);

// Convenience diagnostics.
bool global_cuts_require_sector_phi(const GlobalCutConfig& cfg = default_global_cuts());
bool global_cuts_require_auxiliary_kinematics(const GlobalCutConfig& cfg = default_global_cuts());
std::string global_cuts_analysis_tag(const GlobalCutConfig& cfg = default_global_cuts());

// Sector helpers. Input phi is in radians.
int fd_sector_from_phi_rad(double phi_rad);
int cd_sector_from_phi_rad(double phi_rad);

// dvcsgen ycol helper (exposed for diagnostics/plotting).
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

// Legacy overload: no kinematics, no topology, no sectors.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg);

// Legacy overload: with ycol kinematics, but no topology/sectors.
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

// Topology-aware overloads without sector phi.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        const GlobalCutConfig& cfg);

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

// Full auxiliary-aware overloads. These are preferred whenever
// enable_auxiliary_fiducial_cuts may be true because they provide all particle
// momenta/angles needed by that cut set.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        double Ebeam,
                        double e_p,
                        double e_theta,
                        double e_phi,
                        double p1_theta,
                        double p1_phi,
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
                        double p1_theta,
                        double p1_phi,
                        double p2_p,
                        double p2_theta,
                        double p2_phi,
                        const GlobalCutConfig& cfg);

// Full sector-aware overloads. These are the preferred overloads for production
// event loops because they support topology, ycol, and particle-sector filters.
// They intentionally fail if enable_auxiliary_fiducial_cuts=true because the
// auxiliary cut set also requires p1_theta.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        double e_phi,
                        double p1_phi,
                        double p2_phi,
                        const GlobalCutConfig& cfg);

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        double Ebeam,
                        double e_p,
                        double e_theta,
                        double e_phi,
                        double p1_phi,
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
                        double p1_phi,
                        double p2_p,
                        double p2_theta,
                        double p2_phi,
                        const GlobalCutConfig& cfg);

// Run exclusion helper.
bool is_excluded_run(int runnum, const GlobalCutConfig& cfg = default_global_cuts());

// ROOT TCut helper. Only valid when all cuts are representable as a simple TCut.
std::string global_cuts_tcut(const GlobalCutConfig& cfg = default_global_cuts());

// JSON writer for config echoing.
void write_global_cuts_config_json(const std::string& out_json_dir,
                                  const GlobalCutConfig& cfg = default_global_cuts());

#endif // GLOBAL_CUTS_H
