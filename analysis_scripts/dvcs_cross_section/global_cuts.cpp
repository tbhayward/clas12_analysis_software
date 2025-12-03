#include "global_cuts.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_set>

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

// NEW: global helper, using cfg.excluded_runs
bool is_excluded_run(int runnum, const GlobalCutConfig& cfg) {
    // Build a set once per process using the cfg values.
    static std::unordered_set<int> s_runs;
    static bool initialized = false;

    if (!initialized) {
        initialized = true;
        for (int r : cfg.excluded_runs) {
            s_runs.insert(r);
        }
        std::cout << "[global_cuts] excluded_runs size = "
                  << s_runs.size() << std::endl;
    }

    return s_runs.find(runnum) != s_runs.end();
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
        << "  \"pTmiss_max\": " << std::setprecision(6) << cfg.pTmiss_max << ",\n"
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