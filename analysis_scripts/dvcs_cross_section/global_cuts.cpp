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

static GlobalCutConfig g_default_cfg; // default-initialized to values in header

const GlobalCutConfig& default_global_cuts() {
    return g_default_cfg;
}

static void fatal_global(const std::string& msg) {
    throw std::runtime_error(std::string("[global_cuts] FATAL: ") + msg);
}

static double beam_energy_for_period_label(const std::string& period_label) {
    // Deterministic mapping only. No fallbacks, no heuristics.
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
    const double me = 0.00051099895; // (GeV)
    if (Ebeam_GeV <= me) {
        fatal_global("Ebeam <= me (GeV)");
    }

    const double pz = std::sqrt(Ebeam_GeV * Ebeam_GeV - me * me);
    return TLorentzVector(0.0, 0.0, pz, Ebeam_GeV);
}

static TLorentzVector make_p4_from_p_theta_phi(double p,
                                               double theta_rad,
                                               double phi_rad,
                                               double mass_GeV) {
    if (!(p >= 0.0)) {
        fatal_global("make_p4_from_p_theta_phi: p < 0");
    }

    const double px = p * std::sin(theta_rad) * std::cos(phi_rad);
    const double py = p * std::sin(theta_rad) * std::sin(phi_rad);
    const double pz = p * std::cos(theta_rad);
    const double E  = std::sqrt(p * p + mass_GeV * mass_GeV);

    return TLorentzVector(px, py, pz, E);
}

// P2_pos = 2 (kprime . qgamma) / Q2_calc
// where Q2_calc = - (k - kprime)^2 computed from Ebeam.
static double compute_P2_pos(double Ebeam,
                             double e_p,
                             double e_theta,
                             double e_phi,
                             double p2_p,
                             double p2_theta,
                             double p2_phi) {
    const double me = 0.00051099895; // (GeV)

    const TLorentzVector k      = make_beam_electron(Ebeam);
    const TLorentzVector kprime = make_p4_from_p_theta_phi(e_p,  e_theta,  e_phi,  me);
    const TLorentzVector qgamma = make_p4_from_p_theta_phi(p2_p, p2_theta, p2_phi, 0.0);

    const TLorentzVector q = k - kprime;
    const double Q2_calc = -q.M2();

    if (!(Q2_calc > 0.0) || !std::isfinite(Q2_calc)) {
        fatal_global("compute_P2_pos: non-physical Q2_calc");
    }

    const double dot = kprime.Dot(qgamma);
    const double P2_pos = (2.0 * dot) / Q2_calc;

    if (!std::isfinite(P2_pos)) {
        fatal_global("compute_P2_pos: non-finite P2_pos");
    }

    return P2_pos;
}

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg) {
    if (cfg.enable_dvcsgen_ycol_cut) {
        fatal_global("dvcsgen ycol cut enabled but called base passes_global_cuts() without kinematics");
    }

    if ((-t1) >= cfg.t1_abs_max) return false;
    if (open_angle_ep2_deg <= cfg.open_angle_min_deg) return false;
    if (pTmiss > cfg.pTmiss_max) return false;

    return true;
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
    if ((-t1) >= cfg.t1_abs_max) return false;
    if (open_angle_ep2_deg <= cfg.open_angle_min_deg) return false;
    if (pTmiss > cfg.pTmiss_max) return false;

    if (cfg.enable_dvcsgen_ycol_cut) {
        const double P2_pos = compute_P2_pos(Ebeam, e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi);
        if (!(P2_pos > cfg.dvcsgen_ycol_cut)) return false;
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
    const double Ebeam = beam_energy_for_period_label(period_label);

    return passes_global_cuts(t1,
                              open_angle_ep2_deg,
                              pTmiss,
                              Ebeam,
                              e_p,
                              e_theta,
                              e_phi,
                              p2_p,
                              p2_theta,
                              p2_phi,
                              cfg);
}

bool is_excluded_run(int runnum, const GlobalCutConfig& cfg) {
    static std::unordered_set<int> s_runs;
    static bool initialized = false;

    if (!initialized) {
        initialized = true;
        for (int r : cfg.excluded_runs) {
            s_runs.insert(r);
        }
        std::cout << "[global_cuts] excluded_runs size = " << s_runs.size() << std::endl;
    }

    return s_runs.find(runnum) != s_runs.end();
}

std::string global_cuts_tcut(const GlobalCutConfig& cfg) {
    if (cfg.enable_dvcsgen_ycol_cut) {
        fatal_global("global_cuts_tcut: dvcsgen ycol cut enabled; use C++ cut (needs kinematics)");
    }

    std::ostringstream ss;
    ss << "(-t1) < " << std::fixed << std::setprecision(3) << cfg.t1_abs_max
       << " && open_angle_ep2 > " << std::fixed << std::setprecision(1) << cfg.open_angle_min_deg
       << " && pTmiss <= " << std::fixed << std::setprecision(2) << cfg.pTmiss_max;
    return ss.str();
}

void write_global_cuts_config_json(const std::string& out_json_dir,
                                   const GlobalCutConfig& cfg) {
    std::string path = out_json_dir + "/global_cuts_config.json";
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "[global_cuts] ERROR: cannot open " << path << " for writing.\n";
        return;
    }

    ofs << "{\n"
        << "  \"t1_abs_max\": " << std::setprecision(6) << cfg.t1_abs_max << ",\n"
        << "  \"open_angle_min_deg\": " << std::setprecision(6) << cfg.open_angle_min_deg << ",\n"
        << "  \"pTmiss_max\": " << std::setprecision(6) << cfg.pTmiss_max << ",\n"
        << "  \"enable_dvcsgen_ycol_cut\": " << (cfg.enable_dvcsgen_ycol_cut ? "true" : "false") << ",\n"
        << "  \"dvcsgen_ycol_cut\": " << std::setprecision(6) << cfg.dvcsgen_ycol_cut << ",\n"
        << "  \"beam_energy_by_period_label\": {\n"
        << "    \"sp18_inb\": 10.594,\n"
        << "    \"sp18_out\": 10.594,\n"
        << "    \"fa18_inb\": 10.604,\n"
        << "    \"fa18_out\": 10.604,\n"
        << "    \"sp19_inb\": 10.1998\n"
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