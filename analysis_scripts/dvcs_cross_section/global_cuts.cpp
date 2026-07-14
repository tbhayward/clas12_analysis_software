#include "global_cuts.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <TLorentzVector.h>

static constexpr double RAD2DEG = 180.0 / M_PI;
static GlobalCutConfig g_default_cfg;

const GlobalCutConfig& default_global_cuts() {
    return g_default_cfg;
}

void set_default_global_cuts(const GlobalCutConfig& cfg) {
    g_default_cfg = cfg;
}

static void fatal_global(const std::string& msg) {
    throw std::runtime_error(std::string("[global_cuts] FATAL: ") + msg);
}

static double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) p += 360.0;
    if (p >= 360.0) p = std::nextafter(360.0, 0.0);
    return p;
}

static bool in_wrapped_range(double phi_deg, double lo_deg, double hi_deg) {
    const double p  = wrap_phi_deg(phi_deg);
    const double lo = wrap_phi_deg(lo_deg);
    const double hi = wrap_phi_deg(hi_deg);
    if (hi > lo) return (p >= lo && p < hi);
    return (p >= lo || p < hi);
}

int fd_sector_from_phi_rad(double phi_rad) {
    const double p = wrap_phi_deg(phi_rad * RAD2DEG);
    if (in_wrapped_range(p, 330.0,  30.0)) return 1;
    if (in_wrapped_range(p,  30.0,  90.0)) return 2;
    if (in_wrapped_range(p,  90.0, 150.0)) return 3;
    if (in_wrapped_range(p, 150.0, 210.0)) return 4;
    if (in_wrapped_range(p, 210.0, 270.0)) return 5;
    if (in_wrapped_range(p, 270.0, 330.0)) return 6;
    return -1;
}

int cd_sector_from_phi_rad(double phi_rad) {
    const double p = wrap_phi_deg(phi_rad * RAD2DEG);
    if (in_wrapped_range(p, 272.5,  32.5)) return 1;
    if (in_wrapped_range(p,  32.5, 150.5)) return 2;
    if (in_wrapped_range(p, 150.5, 272.5)) return 3;
    return -1;
}

bool global_cuts_apply_sp18_out_sector_quality_cuts(const GlobalCutConfig& cfg) {
    const bool diagnostic_study_active =
        cfg.enable_topology_filter ||
        cfg.enable_electron_fd_sector_filter ||
        cfg.enable_proton_fd_sector_filter ||
        cfg.enable_proton_cd_sector_filter ||
        cfg.enable_photon_fd_sector_filter;

    return cfg.enable_sp18_out_sector_quality_cuts && !diagnostic_study_active;
}

bool global_cuts_require_sector_phi(const GlobalCutConfig& cfg) {
    return global_cuts_apply_sp18_out_sector_quality_cuts(cfg) ||
           cfg.enable_electron_fd_sector_filter ||
           cfg.enable_proton_fd_sector_filter   ||
           cfg.enable_proton_cd_sector_filter   ||
           cfg.enable_photon_fd_sector_filter   ||
           (cfg.enable_auxiliary_fiducial_cuts &&
            cfg.auxiliary_require_distinct_fd_sectors);
}

bool global_cuts_require_auxiliary_kinematics(const GlobalCutConfig& cfg) {
    return cfg.enable_auxiliary_fiducial_cuts;
}

static void validate_cfg_or_fatal(const GlobalCutConfig& cfg) {
    auto ok_det = [](int d)->bool { return d == 0 || d == 1 || d == 2; };
    auto ok_fd  = [](int s)->bool { return s >= 1 && s <= 6; };
    auto ok_cd  = [](int s)->bool { return s >= 1 && s <= 3; };

    if (cfg.enable_topology_filter) {
        if (!ok_det(cfg.required_detector1) || !ok_det(cfg.required_detector2)) {
            std::ostringstream ss;
            ss << "enable_topology_filter is true, but required_detector1/2 are invalid: ("
               << cfg.required_detector1 << ", " << cfg.required_detector2 << ").";
            fatal_global(ss.str());
        }
    }

    if (cfg.enable_electron_fd_sector_filter && !ok_fd(cfg.electron_fd_sector)) {
        fatal_global("electron_fd_sector must be in [1,6]");
    }
    if (cfg.enable_proton_fd_sector_filter && !ok_fd(cfg.proton_fd_sector)) {
        fatal_global("proton_fd_sector must be in [1,6]");
    }
    if (cfg.enable_proton_cd_sector_filter && !ok_cd(cfg.proton_cd_sector)) {
        fatal_global("proton_cd_sector must be in [1,3]");
    }
    if (cfg.enable_photon_fd_sector_filter && !ok_fd(cfg.photon_fd_sector)) {
        fatal_global("photon_fd_sector must be in [1,6]");
    }

    if (cfg.enable_proton_fd_sector_filter && cfg.enable_proton_cd_sector_filter) {
        fatal_global("proton FD-sector and proton CD-sector filters are mutually exclusive.");
    }

    if (cfg.enable_auxiliary_fiducial_cuts) {
        auto check_range = [](double lo, double hi) {
            return std::isfinite(lo) && std::isfinite(hi) && lo < hi;
        };
        if (!check_range(cfg.auxiliary_e_theta_min_deg, cfg.auxiliary_e_theta_max_deg)) {
            fatal_global("invalid auxiliary electron theta range");
        }
        if (!check_range(cfg.auxiliary_fd_proton_theta_min_deg, cfg.auxiliary_fd_proton_theta_max_deg)) {
            fatal_global("invalid auxiliary FD-proton theta range");
        }
        if (!check_range(cfg.auxiliary_fd_photon_theta_min_deg, cfg.auxiliary_fd_photon_theta_max_deg)) {
            fatal_global("invalid auxiliary FD-photon theta range");
        }
        if (!check_range(cfg.auxiliary_cd_proton_theta_min_deg, cfg.auxiliary_cd_proton_theta_max_deg)) {
            fatal_global("invalid auxiliary CD-proton theta range");
        }
        if (!std::isfinite(cfg.auxiliary_ft_photon_p_min_GeV) ||
            cfg.auxiliary_ft_photon_p_min_GeV < 0.0) {
            fatal_global("invalid auxiliary FT-photon momentum threshold");
        }
    }
}

std::string global_cuts_analysis_tag(const GlobalCutConfig& cfg) {
    validate_cfg_or_fatal(cfg);

    std::vector<std::string> pieces;
    if (cfg.enable_topology_filter) {
        if (cfg.required_detector1 == 1 && cfg.required_detector2 == 1) pieces.push_back("topo_FD_FD");
        else if (cfg.required_detector1 == 2 && cfg.required_detector2 == 1) pieces.push_back("topo_CD_FD");
        else if (cfg.required_detector1 == 2 && cfg.required_detector2 == 0) pieces.push_back("topo_CD_FT");
        else {
            std::ostringstream ss;
            ss << "topo_det" << cfg.required_detector1 << "_det" << cfg.required_detector2;
            pieces.push_back(ss.str());
        }
    }
    if (cfg.enable_electron_fd_sector_filter) pieces.push_back("eFDsec" + std::to_string(cfg.electron_fd_sector));
    if (cfg.enable_proton_fd_sector_filter)   pieces.push_back("pFDsec" + std::to_string(cfg.proton_fd_sector));
    if (cfg.enable_proton_cd_sector_filter)   pieces.push_back("pCDsec" + std::to_string(cfg.proton_cd_sector));
    if (cfg.enable_photon_fd_sector_filter)   pieces.push_back("gFDsec" + std::to_string(cfg.photon_fd_sector));
    if (cfg.enable_auxiliary_fiducial_cuts)   pieces.push_back("auxfid");
    if (global_cuts_apply_sp18_out_sector_quality_cuts(cfg)) pieces.push_back("sp18out_sector_quality");
    if (pieces.empty()) return "nominal";

    std::ostringstream ss;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (i) ss << "__";
        ss << pieces[i];
    }
    return ss.str();
}

static bool passes_topology_filter(int detector1, int detector2, const GlobalCutConfig& cfg) {
    validate_cfg_or_fatal(cfg);
    if (!cfg.enable_topology_filter) return true;
    return detector1 == cfg.required_detector1 && detector2 == cfg.required_detector2;
}

static bool passes_sp18_out_sector_quality_cuts(const std::string& period_label,
                                                  int detector1,
                                                  int detector2,
                                                  double e_phi,
                                                  const GlobalCutConfig& cfg) {
    if (!global_cuts_apply_sp18_out_sector_quality_cuts(cfg)) return true;
    if (period_label != "sp18_out") return true;

    const int e_sector = fd_sector_from_phi_rad(e_phi);
    if (e_sector < 1) return false;

    if (e_sector == 3) return false;
    if (e_sector == 5 && detector2 == 1) return false;
    if (e_sector == 5 && detector1 == 1) return false;

    return true;
}

static bool passes_sector_filters(int detector1,
                                  int detector2,
                                  double e_phi,
                                  double p1_phi,
                                  double p2_phi,
                                  const GlobalCutConfig& cfg) {
    validate_cfg_or_fatal(cfg);

    if (cfg.enable_electron_fd_sector_filter) {
        if (fd_sector_from_phi_rad(e_phi) != cfg.electron_fd_sector) return false;
    }

    if (cfg.enable_proton_fd_sector_filter) {
        if (detector1 != 1) return false;
        if (fd_sector_from_phi_rad(p1_phi) != cfg.proton_fd_sector) return false;
    }

    if (cfg.enable_proton_cd_sector_filter) {
        if (detector1 != 2) return false;
        if (cd_sector_from_phi_rad(p1_phi) != cfg.proton_cd_sector) return false;
    }

    if (cfg.enable_photon_fd_sector_filter) {
        if (detector2 != 1) return false;
        if (fd_sector_from_phi_rad(p2_phi) != cfg.photon_fd_sector) return false;
    }

    return true;
}

static bool passes_auxiliary_fiducial_cuts(int detector1,
                                           int detector2,
                                           double e_theta,
                                           double e_phi,
                                           double p1_theta,
                                           double p1_phi,
                                           double p2_p,
                                           double p2_theta,
                                           double p2_phi,
                                           const GlobalCutConfig& cfg) {
    validate_cfg_or_fatal(cfg);
    if (!cfg.enable_auxiliary_fiducial_cuts) return true;

    const double e_theta_deg  = e_theta  * RAD2DEG;
    const double p1_theta_deg = p1_theta * RAD2DEG;
    const double p2_theta_deg = p2_theta * RAD2DEG;

    if (!(std::isfinite(e_theta_deg) &&
          std::isfinite(p1_theta_deg) &&
          std::isfinite(p2_theta_deg) &&
          std::isfinite(p2_p))) {
        return false;
    }

    if (!(e_theta_deg > cfg.auxiliary_e_theta_min_deg)) return false;

    // Apply the electron upper-angle requirement only in topologies with an FD
    // photon: CD-FD and FD-FD.
    if (detector2 == 1) {
        if (!(e_theta_deg < cfg.auxiliary_e_theta_max_deg)) return false;
    }

    if (detector1 == 1) {
        if (!(p1_theta_deg > cfg.auxiliary_fd_proton_theta_min_deg)) return false;
        if (!(p1_theta_deg < cfg.auxiliary_fd_proton_theta_max_deg)) return false;
    } else if (detector1 == 2) {
        if (!(p1_theta_deg > cfg.auxiliary_cd_proton_theta_min_deg)) return false;
        if (!(p1_theta_deg < cfg.auxiliary_cd_proton_theta_max_deg)) return false;
    } else {
        return false;
    }

    if (detector2 == 1) {
        // The requested FD-photon lower bound is applied only in FD-FD.
        if (detector1 == 1) {
            if (!(p2_theta_deg > cfg.auxiliary_fd_photon_theta_min_deg)) return false;
        }
        if (!(p2_theta_deg < cfg.auxiliary_fd_photon_theta_max_deg)) return false;
    } else if (detector2 == 0) {
        if (!(p2_p > cfg.auxiliary_ft_photon_p_min_GeV)) return false;
    } else {
        return false;
    }

    if (cfg.auxiliary_require_distinct_fd_sectors) {
        std::vector<int> fd_sectors;
        const int e_sector = fd_sector_from_phi_rad(e_phi);
        if (e_sector < 1) return false;
        fd_sectors.push_back(e_sector);

        if (detector1 == 1) {
            const int p_sector = fd_sector_from_phi_rad(p1_phi);
            if (p_sector < 1) return false;
            fd_sectors.push_back(p_sector);
        }

        if (detector2 == 1) {
            const int g_sector = fd_sector_from_phi_rad(p2_phi);
            if (g_sector < 1) return false;
            fd_sectors.push_back(g_sector);
        }

        for (std::size_t i = 0; i < fd_sectors.size(); ++i) {
            for (std::size_t j = i + 1; j < fd_sectors.size(); ++j) {
                if (fd_sectors[i] == fd_sectors[j]) return false;
            }
        }
    }

    return true;
}

static double beam_energy_for_period_label(const std::string& period_label) {
    if (period_label == "sp18_inb") return 10.594;
    if (period_label == "sp18_out") return 10.594;
    if (period_label == "fa18_inb") return 10.604;
    if (period_label == "fa18_out") return 10.604;
    if (period_label == "sp19_inb") return 10.1998;

    std::ostringstream ss;
    ss << "Unknown period_label '" << period_label << "'. Allowed labels are: "
       << "sp18_inb, sp18_out, fa18_inb, fa18_out, sp19_inb.";
    fatal_global(ss.str());
    return -1.0;
}

static TLorentzVector make_beam_electron(double Ebeam_GeV) {
    const double me = 0.00051099895;
    if (Ebeam_GeV <= me) fatal_global("Ebeam <= me (GeV)");
    const double pz = std::sqrt(Ebeam_GeV * Ebeam_GeV - me * me);
    return TLorentzVector(0.0, 0.0, pz, Ebeam_GeV);
}

static TLorentzVector make_p4_from_p_theta_phi(double p,
                                               double theta_rad,
                                               double phi_rad,
                                               double mass_GeV) {
    if (!(p >= 0.0)) fatal_global("make_p4_from_p_theta_phi: p < 0");
    const double px = p * std::sin(theta_rad) * std::cos(phi_rad);
    const double py = p * std::sin(theta_rad) * std::sin(phi_rad);
    const double pz = p * std::cos(theta_rad);
    const double E  = std::sqrt(p * p + mass_GeV * mass_GeV);
    return TLorentzVector(px, py, pz, E);
}

static double compute_P2_pos(double Ebeam,
                             double e_p,
                             double e_theta,
                             double e_phi,
                             double p2_p,
                             double p2_theta,
                             double p2_phi) {
    const double me = 0.00051099895;
    const TLorentzVector k      = make_beam_electron(Ebeam);
    const TLorentzVector kprime = make_p4_from_p_theta_phi(e_p,  e_theta,  e_phi,  me);
    const TLorentzVector qgamma = make_p4_from_p_theta_phi(p2_p, p2_theta, p2_phi, 0.0);
    const TLorentzVector q = k - kprime;
    const double Q2_calc = -q.M2();
    if (!(Q2_calc > 0.0) || !std::isfinite(Q2_calc)) fatal_global("compute_P2_pos: non-physical Q2_calc");
    const double P2_pos = (2.0 * kprime.Dot(qgamma)) / Q2_calc;
    if (!std::isfinite(P2_pos)) fatal_global("compute_P2_pos: non-finite P2_pos");
    return P2_pos;
}

double dvcsgen_ycol_value(double Ebeam,
                          double e_p,
                          double e_theta,
                          double e_phi,
                          double p2_p,
                          double p2_theta,
                          double p2_phi,
                          const GlobalCutConfig& cfg) {
    (void)cfg;
    return compute_P2_pos(Ebeam, e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi);
}

double dvcsgen_ycol_value(const std::string& period_label,
                          double e_p,
                          double e_theta,
                          double e_phi,
                          double p2_p,
                          double p2_theta,
                          double p2_phi,
                          const GlobalCutConfig& cfg) {
    return dvcsgen_ycol_value(beam_energy_for_period_label(period_label),
                              e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi, cfg);
}

static bool passes_basic_cuts(double t1,
                              double open_angle_ep2_deg,
                              double pTmiss,
                              const GlobalCutConfig& cfg) {
    validate_cfg_or_fatal(cfg);
    if ((-t1) >= cfg.t1_abs_max) return false;
    if (open_angle_ep2_deg <= cfg.open_angle_min_deg) return false;
    if (cfg.enable_pTmiss_cut && pTmiss > cfg.pTmiss_max) return false;
    return true;
}

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg) {
    if (cfg.enable_dvcsgen_ycol_cut) fatal_global("dvcsgen ycol cut enabled but called without kinematics");
    if (cfg.enable_topology_filter) fatal_global("topology filter enabled but called without detector1/detector2");
    if (global_cuts_require_sector_phi(cfg)) fatal_global("sector filter enabled but called without particle phi branches");
    return passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg);
}

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
                        const GlobalCutConfig& cfg) {
    if (cfg.enable_topology_filter) fatal_global("topology filter enabled but called without detector1/detector2");
    if (global_cuts_require_sector_phi(cfg)) fatal_global("sector filter enabled but called without full particle phi branches");
    if (!passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg)) return false;
    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(compute_P2_pos(Ebeam, e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi) > cfg.dvcsgen_ycol_cut)) return false;
    }
    return true;
}

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
                        const GlobalCutConfig& cfg) {
    return passes_global_cuts(t1, open_angle_ep2_deg, pTmiss,
                              beam_energy_for_period_label(period_label),
                              e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi, cfg);
}

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        const GlobalCutConfig& cfg) {
    if (cfg.enable_dvcsgen_ycol_cut) fatal_global("dvcsgen ycol cut enabled but called without kinematics");
    if (global_cuts_require_auxiliary_kinematics(cfg)) fatal_global("auxiliary fiducial cuts enabled but called without particle kinematics");
    if (global_cuts_require_sector_phi(cfg)) fatal_global("sector filter enabled but called without particle phi branches");
    if (!passes_topology_filter(detector1, detector2, cfg)) return false;
    return passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg);
}

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        const std::string& period_label,
                        double e_phi,
                        double p1_phi,
                        double p2_phi,
                        const GlobalCutConfig& cfg) {
    if (!passes_sp18_out_sector_quality_cuts(period_label, detector1, detector2, e_phi, cfg)) return false;
    return passes_global_cuts(t1, open_angle_ep2_deg, pTmiss,
                              detector1, detector2,
                              e_phi, p1_phi, p2_phi, cfg);
}

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
                        const GlobalCutConfig& cfg) {
    if (global_cuts_require_sector_phi(cfg)) fatal_global("sector/auxiliary FD-sector cut enabled but called without p1_phi");
    if (global_cuts_require_auxiliary_kinematics(cfg)) fatal_global("auxiliary fiducial cuts enabled but called without p1_theta");
    if (!passes_topology_filter(detector1, detector2, cfg)) return false;
    if (!passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg)) return false;
    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(compute_P2_pos(Ebeam, e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi) > cfg.dvcsgen_ycol_cut)) return false;
    }
    return true;
}

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
                        const GlobalCutConfig& cfg) {
    if (!passes_sp18_out_sector_quality_cuts(period_label, detector1, detector2, e_phi, cfg)) return false;
    return passes_global_cuts(t1, open_angle_ep2_deg, pTmiss,
                              detector1, detector2,
                              beam_energy_for_period_label(period_label),
                              e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi, cfg);
}

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        int detector1,
                        int detector2,
                        double e_phi,
                        double p1_phi,
                        double p2_phi,
                        const GlobalCutConfig& cfg) {
    if (cfg.enable_dvcsgen_ycol_cut) fatal_global("dvcsgen ycol cut enabled but called sector overload without ycol kinematics");
    if (global_cuts_require_auxiliary_kinematics(cfg)) fatal_global("auxiliary fiducial cuts enabled but called without particle theta/momentum branches");
    if (!passes_topology_filter(detector1, detector2, cfg)) return false;
    if (!passes_sector_filters(detector1, detector2, e_phi, p1_phi, p2_phi, cfg)) return false;
    return passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg);
}

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
                        const GlobalCutConfig& cfg) {
    if (!passes_topology_filter(detector1, detector2, cfg)) return false;
    if (!passes_sector_filters(detector1, detector2, e_phi, p1_phi, p2_phi, cfg)) return false;
    if (!passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg)) return false;
    if (!passes_auxiliary_fiducial_cuts(detector1, detector2,
                                        e_theta, e_phi,
                                        p1_theta, p1_phi,
                                        p2_p, p2_theta, p2_phi,
                                        cfg)) return false;
    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(compute_P2_pos(Ebeam, e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi) > cfg.dvcsgen_ycol_cut)) return false;
    }
    return true;
}

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
                        const GlobalCutConfig& cfg) {
    if (!passes_sp18_out_sector_quality_cuts(period_label, detector1, detector2, e_phi, cfg)) return false;
    return passes_global_cuts(t1, open_angle_ep2_deg, pTmiss,
                              detector1, detector2,
                              beam_energy_for_period_label(period_label),
                              e_p, e_theta, e_phi,
                              p1_theta, p1_phi,
                              p2_p, p2_theta, p2_phi,
                              cfg);
}

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
                        const GlobalCutConfig& cfg) {
    if (global_cuts_require_auxiliary_kinematics(cfg)) fatal_global("auxiliary fiducial cuts enabled but called without p1_theta");
    if (!passes_topology_filter(detector1, detector2, cfg)) return false;
    if (!passes_sector_filters(detector1, detector2, e_phi, p1_phi, p2_phi, cfg)) return false;
    if (!passes_basic_cuts(t1, open_angle_ep2_deg, pTmiss, cfg)) return false;
    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(compute_P2_pos(Ebeam, e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi) > cfg.dvcsgen_ycol_cut)) return false;
    }
    return true;
}

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
                        const GlobalCutConfig& cfg) {
    return passes_global_cuts(t1, open_angle_ep2_deg, pTmiss,
                              detector1, detector2,
                              beam_energy_for_period_label(period_label),
                              e_p, e_theta, e_phi, p1_phi,
                              p2_p, p2_theta, p2_phi, cfg);
}

bool is_excluded_run(int runnum, const GlobalCutConfig& cfg) {
    for (int r : cfg.excluded_runs) {
        if (r == runnum) return true;
    }
    return false;
}

std::string global_cuts_tcut(const GlobalCutConfig& cfg) {
    if (cfg.enable_dvcsgen_ycol_cut) fatal_global("global_cuts_tcut: dvcsgen ycol cut enabled; use C++ cut");
    if (cfg.enable_topology_filter) fatal_global("global_cuts_tcut: topology filter enabled; use C++ cut");
    if (global_cuts_require_sector_phi(cfg)) fatal_global("global_cuts_tcut: sector filter enabled; use C++ cut");
    if (global_cuts_require_auxiliary_kinematics(cfg)) fatal_global("global_cuts_tcut: auxiliary fiducial cuts enabled; use C++ cut");

    std::ostringstream ss;
    ss << "(-t1) < " << std::fixed << std::setprecision(3) << cfg.t1_abs_max
       << " && open_angle_ep2 > " << std::fixed << std::setprecision(1) << cfg.open_angle_min_deg;
    if (cfg.enable_pTmiss_cut) {
        ss << " && pTmiss <= " << std::fixed << std::setprecision(2) << cfg.pTmiss_max;
    }
    return ss.str();
}

void write_global_cuts_config_json(const std::string& out_json_dir,
                                  const GlobalCutConfig& cfg) {
    const std::string path = out_json_dir + "/global_cuts_config.json";
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "[global_cuts] ERROR: cannot open " << path << " for writing.\n";
        return;
    }

    ofs << "{\n"
        << "  \"analysis_tag\": \"" << global_cuts_analysis_tag(cfg) << "\",\n"
        << "  \"t1_abs_max\": " << std::setprecision(6) << cfg.t1_abs_max << ",\n"
        << "  \"open_angle_min_deg\": " << std::setprecision(6) << cfg.open_angle_min_deg << ",\n"
        << "  \"enable_pTmiss_cut\": " << (cfg.enable_pTmiss_cut ? "true" : "false") << ",\n"
        << "  \"pTmiss_max\": " << std::setprecision(6) << cfg.pTmiss_max << ",\n"
        << "  \"enable_dvcsgen_ycol_cut\": " << (cfg.enable_dvcsgen_ycol_cut ? "true" : "false") << ",\n"
        << "  \"dvcsgen_ycol_cut\": " << std::setprecision(6) << cfg.dvcsgen_ycol_cut << ",\n"
        << "  \"enable_auxiliary_fiducial_cuts\": " << (cfg.enable_auxiliary_fiducial_cuts ? "true" : "false") << ",\n"
        << "  \"enable_sp18_out_sector_quality_cuts\": " << (cfg.enable_sp18_out_sector_quality_cuts ? "true" : "false") << ",\n"
        << "  \"sp18_out_sector_quality_cuts_active\": " << (global_cuts_apply_sp18_out_sector_quality_cuts(cfg) ? "true" : "false") << ",\n"
        << "  \"auxiliary_require_distinct_fd_sectors\": " << (cfg.auxiliary_require_distinct_fd_sectors ? "true" : "false") << ",\n"
        << "  \"auxiliary_e_theta_min_deg\": " << std::setprecision(6) << cfg.auxiliary_e_theta_min_deg << ",\n"
        << "  \"auxiliary_e_theta_max_deg\": " << std::setprecision(6) << cfg.auxiliary_e_theta_max_deg << ",\n"
        << "  \"auxiliary_fd_proton_theta_min_deg\": " << std::setprecision(6) << cfg.auxiliary_fd_proton_theta_min_deg << ",\n"
        << "  \"auxiliary_fd_proton_theta_max_deg\": " << std::setprecision(6) << cfg.auxiliary_fd_proton_theta_max_deg << ",\n"
        << "  \"auxiliary_fd_photon_theta_min_deg\": " << std::setprecision(6) << cfg.auxiliary_fd_photon_theta_min_deg << ",\n"
        << "  \"auxiliary_fd_photon_theta_max_deg\": " << std::setprecision(6) << cfg.auxiliary_fd_photon_theta_max_deg << ",\n"
        << "  \"auxiliary_cd_proton_theta_min_deg\": " << std::setprecision(6) << cfg.auxiliary_cd_proton_theta_min_deg << ",\n"
        << "  \"auxiliary_cd_proton_theta_max_deg\": " << std::setprecision(6) << cfg.auxiliary_cd_proton_theta_max_deg << ",\n"
        << "  \"auxiliary_ft_photon_p_min_GeV\": " << std::setprecision(6) << cfg.auxiliary_ft_photon_p_min_GeV << ",\n"
        << "  \"enable_topology_filter\": " << (cfg.enable_topology_filter ? "true" : "false") << ",\n"
        << "  \"required_detector1\": " << cfg.required_detector1 << ",\n"
        << "  \"required_detector2\": " << cfg.required_detector2 << ",\n"
        << "  \"enable_electron_fd_sector_filter\": " << (cfg.enable_electron_fd_sector_filter ? "true" : "false") << ",\n"
        << "  \"electron_fd_sector\": " << cfg.electron_fd_sector << ",\n"
        << "  \"enable_proton_fd_sector_filter\": " << (cfg.enable_proton_fd_sector_filter ? "true" : "false") << ",\n"
        << "  \"proton_fd_sector\": " << cfg.proton_fd_sector << ",\n"
        << "  \"enable_proton_cd_sector_filter\": " << (cfg.enable_proton_cd_sector_filter ? "true" : "false") << ",\n"
        << "  \"proton_cd_sector\": " << cfg.proton_cd_sector << ",\n"
        << "  \"enable_photon_fd_sector_filter\": " << (cfg.enable_photon_fd_sector_filter ? "true" : "false") << ",\n"
        << "  \"photon_fd_sector\": " << cfg.photon_fd_sector << ",\n"
        << "  \"fd_sector_boundaries_deg\": {\n"
        << "    \"1\": [330.0, 30.0],\n"
        << "    \"2\": [30.0, 90.0],\n"
        << "    \"3\": [90.0, 150.0],\n"
        << "    \"4\": [150.0, 210.0],\n"
        << "    \"5\": [210.0, 270.0],\n"
        << "    \"6\": [270.0, 330.0]\n"
        << "  },\n"
        << "  \"cd_sector_boundaries_deg\": {\n"
        << "    \"1\": [272.5, 32.5],\n"
        << "    \"2\": [32.5, 150.5],\n"
        << "    \"3\": [150.5, 272.5]\n"
        << "  },\n"
        << "  \"excluded_runs\": [\n";

    for (std::size_t i = 0; i < cfg.excluded_runs.size(); ++i) {
        ofs << "    " << cfg.excluded_runs[i];
        if (i + 1 < cfg.excluded_runs.size()) ofs << ",";
        ofs << "\n";
    }

    ofs << "  ]\n"
        << "}\n";

    std::cout << "[global_cuts] Wrote " << path << "\n";
}
