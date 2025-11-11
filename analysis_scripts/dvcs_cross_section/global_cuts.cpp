#include "global_cuts.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

static GlobalCutConfig g_default_cfg; // default-initialized to values in header

const GlobalCutConfig& default_global_cuts() {
    return g_default_cfg;
}

bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg) {
    // (-t1) < cfg.t1_abs_max  <=>  reject if (-t1) >= cfg.t1_abs_max
    if ((-t1) >= cfg.t1_abs_max) return false;

    // open_angle_ep2 > cfg.open_angle_min_deg
    if (open_angle_ep2_deg <= cfg.open_angle_min_deg) return false;

    // pTmiss <= cfg.pTmiss_max
    if (pTmiss > cfg.pTmiss_max) return false;

    return true;
}

std::string global_cuts_tcut(const GlobalCutConfig& cfg) {
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
        << "  \"pTmiss_max\": " << std::setprecision(6) << cfg.pTmiss_max << "\n"
        << "}\n";
    std::cout << "[global_cuts] Wrote " << path << "\n";
}