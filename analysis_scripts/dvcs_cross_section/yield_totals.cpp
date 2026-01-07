// yield_totals.cpp
// ------------------------------------------------------------
// Count reconstructed events after exclusivity cuts, grouped by
// beam current (nA) for DATA and by period for reconstructed MC.
//
// Cuts:
//   - Global kinematic cuts via passes_global_cuts(...) from global_cuts.h,
//     using the SAME pattern as other modules (e.g. bin_means.cpp):
//        * get config via default_global_cuts()
//        * if cfg.enable_dvcsgen_ycol_cut is enabled, call the kinematics-aware
//          overload using the P2_pos mirror inputs:
//              e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi
//          with Ebeam resolved deterministically from the canonical period label
//          (sp18_inb, sp18_out, fa18_inb, fa18_out, sp19_inb).
//        * otherwise call the base overload (t1/open/pTmiss only).
//     NOTE: No "nu" is used anywhere in the P2 implementation.
//
//   - 3-sigma exclusivity cuts from combined_cuts.json, using topology keys
//     "DVCS_<PeriodDir>_<TopoDir>" and the "data"/"mc" blocks.
//     Apply within_3sigma to any variables present (and bound).
//
// Public API is unchanged (no main.cpp changes).
// ------------------------------------------------------------

#include "yield_totals.h"

#include "periods.h"
#include "global_cuts.h"

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TROOT.h> // ROOT::EnableThreadSafety

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using json = nlohmann::json;

// ---------------- small helpers ----------------

[[noreturn]] void fatal(const std::string& msg) {
    std::cerr << "[yield_totals] FATAL: " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static std::string to_lower_nospace(std::string s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (c == ' ' || c == '\t' || c == '_') continue;
        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    return out;
}

// Canonical period directory names (same logic as in total_counts.cpp)
static std::string canonical_period_dir(const std::string& L) {
    const std::string k = to_lower_nospace(L);
    if (k == "fa18inb")      return "Fa18_Inb";
    if (k == "fa18out")      return "Fa18_Out";
    if (k == "fa18inbsupp")  return "Fa18_Inb_Supp";
    if (k == "sp18inb")      return "Sp18_Inb";
    if (k == "sp18out")      return "Sp18_Out";
    if (k == "sp19inb")      return "Sp19_Inb";
    std::string s = L;
    for (char& c : s) if (c == ' ') c = '_';
    return s;
}

static inline std::string period_dir_for_label(const std::string& L) {
    return canonical_period_dir(L);
}

// Canonical period labels for global_cuts (deterministic; fail-fast if unknown)
static std::string canonical_period_label(const std::string& L) {
    const std::string k = to_lower_nospace(L);
    if (k == "sp18inb") return "sp18_inb";
    if (k == "sp18out") return "sp18_out";
    if (k == "fa18inb") return "fa18_inb";
    if (k == "fa18out") return "fa18_out";
    if (k == "sp19inb") return "sp19_inb";
    if (k == "fa18inbsupp") return "fa18_inb"; // only used if someone forgets to skip it
    fatal("Unknown period label '" + L + "' for canonical_period_label().");
    return "";
}

// ---------------- topology ----------------

enum class Topology { FD_FD = 0, CD_FD = 1, CD_FT = 2 };

static inline const char* topo_dir(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

// ---------------- sigma cuts ----------------

struct SigmaCut {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std  = std::numeric_limits<double>::quiet_NaN();
};

using VarCutMap = std::unordered_map<std::string, SigmaCut>;

struct TopoSigma {
    VarCutMap data;
    VarCutMap mc;
};

using TopoCutMap = std::unordered_map<std::string, TopoSigma>;

static TopoCutMap load_sigma_cuts_both(const std::string& path) {
    TopoCutMap out;
    if (path.empty()) fatal("combined_cuts_json path is empty.");

    std::ifstream fin(path);
    if (!fin.is_open()) {
        fatal("Could not open combined cuts JSON: " + path);
    }

    json j;
    fin >> j;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const json& block = it.value();

        TopoSigma ts;

        if (block.contains("data") && block["data"].is_object()) {
            for (auto vit = block["data"].begin(); vit != block["data"].end(); ++vit) {
                const std::string vname = vit.key();
                const json& vs = vit.value();
                SigmaCut sc;
                if (vs.contains("mean") && vs["mean"].is_number()) sc.mean = vs["mean"].get<double>();
                if (vs.contains("std")  && vs["std"].is_number())  sc.std  = vs["std"].get<double>();
                if (std::isfinite(sc.std) && sc.std > 0.0) {
                    ts.data[vname] = sc;
                }
            }
        }

        if (block.contains("mc") && block["mc"].is_object()) {
            for (auto vit = block["mc"].begin(); vit != block["mc"].end(); ++vit) {
                const std::string vname = vit.key();
                const json& vs = vit.value();
                SigmaCut sc;
                if (vs.contains("mean") && vs["mean"].is_number()) sc.mean = vs["mean"].get<double>();
                if (vs.contains("std")  && vs["std"].is_number())  sc.std  = vs["std"].get<double>();
                if (std::isfinite(sc.std) && sc.std > 0.0) {
                    ts.mc[vname] = sc;
                }
            }
        }

        if (!ts.data.empty() || !ts.mc.empty()) {
            out[key] = std::move(ts);
        }
    }

    std::cout << "[yield_totals] Loaded sigma cuts for " << out.size()
              << " topology keys from " << path << std::endl;
    return out;
}

enum class CutMode { TwoSided, UpperOnly };

static CutMode var_mode(const std::string& var) {
    // Match the convention used elsewhere (e.g. bin_means.cpp).
    if (var == "Emiss2" || var == "pTmiss" || var == "theta_gamma_gamma") return CutMode::UpperOnly;
    return CutMode::TwoSided;
}

static inline bool within_3sigma(double val, const SigmaCut& sc, CutMode mode) {
    if (!std::isfinite(val) || !std::isfinite(sc.mean) || !std::isfinite(sc.std) || sc.std <= 0.0) {
        return true;
    }
    if (mode == CutMode::UpperOnly) {
        return (val <= sc.mean + 3.0 * sc.std);
    }
    return (std::fabs(val - sc.mean) <= 3.0 * sc.std);
}

// ---------------- branch binding (bin_means-style) ----------------

static std::mutex g_root_bind_mutex;

struct BranchBinder {
    // Required for DATA current mapping + run blacklist
    int runnum = 0;              bool has_runnum = false;

    // Required for topology-dependent sigma cuts
    int detector1 = 0;           bool has_det1 = false;
    int detector2 = 0;           bool has_det2 = false;

    // Global cuts
    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open = false; // degrees
    double pTmiss = 0.0;         bool has_pT = false;

    // Sigma-cut variables (optional; applied if present AND present in JSON)
    double Emiss2 = 0.0;         bool has_Em2 = false;
    double Mx2 = 0.0;            bool has_Mx2 = false;
    double Mx2_1 = 0.0;          bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;          bool has_Mx2_2 = false;
    double xF = 0.0;             bool has_xF = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gg = false;

    // P2 cut inputs (only required if cfg.enable_dvcsgen_ycol_cut is true)
    double e_p = 0.0;            bool has_e_p = false;
    double e_theta = 0.0;        bool has_e_theta = false;
    double e_phi = 0.0;          bool has_e_phi = false;

    double p2_p = 0.0;           bool has_p2_p = false;
    double p2_theta = 0.0;       bool has_p2_theta = false;
    double p2_phi = 0.0;         bool has_p2_phi = false;

    void bind(TTree* t, bool require_runnum) {
        if (!t) fatal("BranchBinder::bind: null tree.");

        std::lock_guard<std::mutex> lock(g_root_bind_mutex);

        t->SetBranchStatus("*", 0);

        auto enable = [&](const char* n){
            if (t->GetBranch(n)) t->SetBranchStatus(n, 1);
        };

        // Always attempt to enable (guarded) what we might use
        enable("runnum");
        enable("detector1");
        enable("detector2");
        enable("t1");
        enable("open_angle_ep2");
        enable("pTmiss");

        enable("Emiss2");
        enable("Mx2");
        enable("Mx2_1");
        enable("Mx2_2");
        enable("xF");
        enable("theta_gamma_gamma");

        enable("e_p");
        enable("e_theta");
        enable("e_phi");
        enable("p2_p");
        enable("p2_theta");
        enable("p2_phi");

        t->SetCacheSize(0);

        auto bindD = [&](const char* n, double* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bindI = [&](const char* n, int* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bindI("runnum",    &runnum,    has_runnum);
        bindI("detector1", &detector1, has_det1);
        bindI("detector2", &detector2, has_det2);

        bindD("t1",             &t1,             has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open);
        bindD("pTmiss",         &pTmiss,         has_pT);

        bindD("Emiss2",            &Emiss2,            has_Em2);
        bindD("Mx2",               &Mx2,               has_Mx2);
        bindD("Mx2_1",             &Mx2_1,             has_Mx2_1);
        bindD("Mx2_2",             &Mx2_2,             has_Mx2_2);
        bindD("xF",                &xF,                has_xF);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gg);

        bindD("e_p",     &e_p,     has_e_p);
        bindD("e_theta", &e_theta, has_e_theta);
        bindD("e_phi",   &e_phi,   has_e_phi);

        bindD("p2_p",     &p2_p,     has_p2_p);
        bindD("p2_theta", &p2_theta, has_p2_theta);
        bindD("p2_phi",   &p2_phi,   has_p2_phi);

        if (require_runnum && !has_runnum) {
            fatal("Missing required branch 'runnum' in DATA tree.");
        }
        if (!(has_det1 && has_det2)) {
            fatal("Missing required branches detector1/detector2 in DVCS tree (needed for topology).");
        }
        if (!(has_t1 && has_open && has_pT)) {
            fatal("Missing required branches t1/open_angle_ep2/pTmiss in DVCS tree (needed for global cuts).");
        }
    }

    int topo_index() const {
        if (!has_det1 || !has_det2) return -1;
        if (detector1 == 1 && detector2 == 1) return static_cast<int>(Topology::FD_FD);
        if (detector1 == 2 && detector2 == 1) return static_cast<int>(Topology::CD_FD);
        if (detector1 == 2 && detector2 == 0) return static_cast<int>(Topology::CD_FT);
        return -1;
    }
};

// Apply sigma cuts to any variables present in BOTH: (a) JSON map and (b) tree binding.
static bool passes_sigma_cuts(const TopoSigma& ts, bool is_mc, const BranchBinder& b) {
    const VarCutMap& m = is_mc ? ts.mc : ts.data;
    if (m.empty()) return true;

    auto check = [&](const char* vname, bool has_val, double val) -> bool {
        auto it = m.find(vname);
        if (it == m.end()) return true;      // not requested by JSON
        if (!has_val) return true;           // branch not present: skip explicitly
        const SigmaCut& sc = it->second;
        const CutMode mode = var_mode(vname);
        return within_3sigma(val, sc, mode);
    };

    if (!check("Emiss2",            b.has_Em2,      b.Emiss2)) return false;
    if (!check("Mx2",               b.has_Mx2,      b.Mx2)) return false;
    if (!check("Mx2_1",             b.has_Mx2_1,    b.Mx2_1)) return false;
    if (!check("Mx2_2",             b.has_Mx2_2,    b.Mx2_2)) return false;
    if (!check("pTmiss",            b.has_pT,       b.pTmiss)) return false;
    if (!check("xF",                b.has_xF,       b.xF)) return false;
    if (!check("theta_gamma_gamma", b.has_theta_gg, b.theta_gamma_gamma)) return false;

    // If the JSON contains additional variables not listed above, we ignore them here
    // (explicit behavior; no heuristics for unknown names).
    return true;
}

// ---------------- global cuts wrapper (consistent with other modules) ----------------

static inline bool passes_global(bool is_mc,
                                 int runnum,
                                 const std::string& period_label_internal,
                                 const BranchBinder& b,
                                 const GlobalCutConfig& cfg) {
    // Run blacklist only applies to data (explicit behavior).
    if (!is_mc) {
        if (is_excluded_run(runnum, cfg)) return false;
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal(std::string("P2 cut enabled but required branches are missing for period_label '") +
                  period_label_internal + "'. Missing one or more of: e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  period_label_internal,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss, cfg);
}

// ---------------- run -> current maps ----------------

// Fa18 Inb
static int current_fa18_inb(int run, bool& ok) {
    static const std::unordered_map<int, int> m = {
        // 40 nA
        {5335, 40}, {5339, 40}, {5341, 40},
        {5340, 40}, {5342, 40}, {5343, 40}, {5344, 40}, {5345, 40},

        // 45 nA
        {5032, 45}, {5036, 45}, {5038, 45}, {5039, 45}, {5040, 45}, {5041, 45},
        {5043, 45}, {5045, 45}, {5052, 45}, {5053, 45}, {5116, 45}, {5117, 45},
        {5119, 45}, {5120, 45}, {5124, 45}, {5125, 45}, {5126, 45}, {5127, 45},
        {5139, 45}, {5153, 45}, {5158, 45}, {5162, 45}, {5163, 45}, {5164, 45},
        {5181, 45}, {5191, 45}, {5193, 45}, {5195, 45}, {5196, 45}, {5197, 45},
        {5198, 45}, {5199, 45}, {5200, 45}, {5201, 45}, {5202, 45}, {5203, 45},
        {5204, 45}, {5205, 45}, {5206, 45}, {5208, 45}, {5211, 45}, {5212, 45},
        {5215, 45}, {5216, 45}, {5219, 45}, {5220, 45}, {5221, 45}, {5222, 45},
        {5223, 45}, {5230, 45}, {5231, 45}, {5232, 45}, {5233, 45}, {5234, 45},
        {5235, 45}, {5237, 45}, {5238, 45}, {5247, 45}, {5248, 45}, {5249, 45},
        {5252, 45}, {5253, 45}, {5257, 45}, {5258, 45}, {5259, 45}, {5261, 45},
        {5262, 45}, {5303, 45}, {5304, 45}, {5305, 45}, {5306, 45}, {5307, 45},
        {5310, 45}, {5311, 45}, {5315, 45}, {5317, 45}, {5318, 45}, {5319, 45},
        {5320, 45}, {5323, 45}, {5324, 45}, {5333, 45}, {5334, 45}, {5346, 45},
        {5347, 45}, {5349, 45}, {5351, 45}, {5354, 45}, {5355, 45}, {5367, 45},

        // Extra Fa18 Inb 45 nA runs
        {5046, 45}, {5047, 45}, {5051, 45},
        {5128, 45}, {5129, 45}, {5130, 45},
        {5159, 45}, {5160, 45},
        {5165, 45}, {5166, 45}, {5167, 45}, {5168, 45}, {5169, 45},
        {5180, 45}, {5182, 45}, {5183, 45},
        {5190, 45},
        {5239, 45},
        {5336, 45},

        // 50 nA
        {5356, 50}, {5357, 50}, {5358, 50}, {5359, 50}, {5360, 50}, {5361, 50},
        {5362, 50}, {5366, 50},

        // 55 nA
        {5368, 55}, {5369, 55}, {5372, 55}, {5373, 55}, {5374, 55}, {5375, 55},
        {5376, 55}, {5377, 55}, {5378, 55}, {5379, 55}, {5380, 55}, {5381, 55},
        {5382, 55}, {5383, 55}, {5386, 55}, {5390, 55}, {5391, 55}, {5392, 55},
        {5393, 55}, {5398, 55}, {5400, 55}, {5401, 55}, {5403, 55}, {5404, 55},
        {5406, 55}, {5407, 55}
    };
    auto it = m.find(run);
    if (it == m.end()) {
        ok = false;
        return 0;
    }
    ok = true;
    return it->second;
}

// Fa18 Out
static int current_fa18_out(int run, bool& ok) {
    static const std::unordered_map<int, int> m = {
        // 40 nA
        {5423, 40}, {5424, 40}, {5425, 40}, {5426, 40}, {5428, 40}, {5429, 40},
        {5430, 40}, {5432, 40}, {5434, 40}, {5435, 40}, {5436, 40}, {5437, 40},
        {5438, 40}, {5440, 40}, {5441, 40}, {5442, 40}, {5445, 40}, {5447, 40},
        {5449, 40}, {5450, 40}, {5451, 40}, {5452, 40}, {5453, 40}, {5454, 40},
        {5455, 40}, {5460, 40}, {5464, 40}, {5465, 40}, {5466, 40}, {5467, 40},
        {5468, 40}, {5469, 40}, {5470, 40}, {5471, 40}, {5472, 40}, {5473, 40},
        {5474, 40}, {5475, 40}, {5476, 40}, {5478, 40}, {5479, 40}, {5480, 40},
        {5481, 40}, {5482, 40}, {5483, 40}, {5485, 40}, {5486, 40}, {5487, 40},
        {5497, 40}, {5498, 40}, {5499, 40}, {5500, 40}, {5504, 40},

        // Extra Fa18 Out 40 nA runs
        {5448, 40}, {5495, 40}, {5496, 40},

        // 50 nA
        {5507, 50}, {5516, 50}, {5517, 50}, {5518, 50}, {5519, 50}, {5520, 50},
        {5521, 50}, {5522, 50}, {5523, 50}, {5524, 50}, {5525, 50}, {5526, 50},
        {5527, 50}, {5528, 50}, {5530, 50}, {5532, 50}, {5533, 50}, {5534, 50},
        {5535, 50}, {5536, 50}, {5537, 50}, {5538, 50}, {5540, 50}, {5541, 50},
        {5543, 50}, {5544, 50}, {5545, 50}, {5546, 50}, {5547, 50}, {5548, 50},
        {5549, 50}, {5550, 50}, {5551, 50}, {5552, 50}, {5555, 50}, {5556, 50},
        {5557, 50}, {5558, 50}, {5559, 50}, {5562, 50}, {5569, 50}, {5570, 50},
        {5571, 50}, {5572, 50}, {5573, 50}, {5574, 50}, {5577, 50}, {5578, 50},
        {5591, 50}, {5592, 50}, {5594, 50}, {5597, 50}, {5598, 50}, {5600, 50},
        {5601, 50}, {5602, 50}, {5603, 50}, {5604, 50}, {5606, 50}, {5607, 50},
        {5611, 50}, {5612, 50}, {5613, 50}, {5614, 50}, {5616, 50}, {5618, 50},
        {5619, 50}, {5624, 50}, {5625, 50}, {5626, 50}, {5627, 50}, {5628, 50},
        {5629, 50}, {5630, 50}, {5631, 50}, {5632, 50}, {5633, 50}, {5635, 50},
        {5637, 50}, {5638, 50}, {5639, 50}, {5641, 50}, {5643, 50}, {5644, 50},
        {5645, 50}, {5646, 50}, {5647, 50}, {5648, 50}, {5649, 50}, {5650, 50},
        {5651, 50}, {5652, 50}, {5654, 50}, {5655, 50}, {5656, 50}, {5662, 50},
        {5663, 50}, {5664, 50}, {5665, 50}, {5666, 50},
        {5615, 50},

        // Extra Fa18 Out 50 nA runs
        {5505, 50}, {5567, 50}, {5617, 50}, {5621, 50}, {5623, 50},

        // Extra Fa18 Out 5 nA run
        {5610, 5}
    };
    auto it = m.find(run);
    if (it == m.end()) {
        ok = false;
        return 0;
    }
    ok = true;
    return it->second;
}

// Resolve current for a given period label and run number.
static bool resolve_current_for_label(const std::string& period_label,
                                      int runnum,
                                      int& current_nA) {
    const std::string k = to_lower_nospace(period_label);

    if (k == "fa18inb") {
        bool ok = false;
        int cur = current_fa18_inb(runnum, ok);
        if (!ok) return false;
        current_nA = cur;
        return true;
    }

    if (k == "fa18out") {
        bool ok = false;
        int cur = current_fa18_out(runnum, ok);
        if (!ok) return false;
        current_nA = cur;
        return true;
    }

    if (k == "sp18out") {
        if (runnum >= 3211 && runnum <= 3293) {
            current_nA = 30;
            return true;
        }
        if (runnum >= 3867 && runnum <= 3987) {
            current_nA = 45;
            return true;
        }
        return false;
    }

    if (k == "sp18inb") {
        if (runnum == 3418) {
            current_nA = 70;
            return true;
        }
        if (runnum == 3421 || runnum == 3422) {
            current_nA = 35;
            return true;
        }
        if (runnum == 3429) {
            current_nA = 50;
            return true;
        }

        if (runnum >= 3306 && runnum <= 3411) {
            current_nA = 35;
            return true;
        }
        if (runnum >= 3431 && runnum <= 4325) {
            current_nA = 50;
            return true;
        }

        return false;
    }

    if (k == "sp19inb") {
        current_nA = 50;
        return true;
    }

    return false;
}

// ---------------- core worker ----------------

struct Totals {
    std::map<std::string, std::map<int, long long>> data_by_period_current;
    std::map<std::string, long long> data_uncategorized;
    std::map<std::string, std::set<int>> data_uncategorized_runs;
    std::map<std::string, long long> mc_by_period;
};

static Totals compute_totals_internal(const std::map<std::string, TTree*>& dvcsDataTrees,
                                      const std::map<std::string, TTree*>& dvcsRecMcTrees,
                                      const TopoCutMap& sigmaCuts,
                                      const GlobalCutConfig& globalCuts) {
    Totals totals;

    auto find_topo_sigma = [&](const std::string& topo_key) -> const TopoSigma& {
        auto it = sigmaCuts.find(topo_key);
        if (it == sigmaCuts.end()) {
            fatal("Missing sigma cuts for key: " + topo_key +
                  " (did you run runAllExclusivityCuts first?)");
        }
        return it->second;
    };

    for (const auto& P : CANONICAL_PERIODS()) {
        const std::string tree_key = P.tree_key;
        const std::string label    = P.label;

        const std::string period_dir = period_dir_for_label(label);

        if (period_dir == "Fa18_Inb_Supp") {
            std::cout << "[yield_totals] Skipping period " << label
                      << " (Fa18 Inb Supp; no sigma cuts defined)." << std::endl;
            continue;
        }

        const std::string period_label_internal = canonical_period_label(label);

        // ---------------- DATA ----------------
        {
            auto itT = dvcsDataTrees.find(tree_key);
            if (itT != dvcsDataTrees.end() && itT->second) {
                TTree* t = itT->second;
                std::cout << "[yield_totals] Processing DATA for period " << label
                          << " (tree key " << tree_key << ") with "
                          << static_cast<long long>(t->GetEntries()) << " entries\n";

                BranchBinder b;
                b.bind(t, /*require_runnum=*/true);

                const Long64_t N = t->GetEntries();
                long long kept = 0;

                for (Long64_t i = 0; i < N; ++i) {
                    if (t->GetEntry(i) <= 0) continue;

                    const int run = b.runnum;

                    if (!passes_global(/*is_mc=*/false,
                                       run,
                                       period_label_internal,
                                       b,
                                       globalCuts)) {
                        continue;
                    }

                    const int topo_idx = b.topo_index();
                    if (topo_idx < 0 || topo_idx > 2) continue;

                    const std::string topo_key =
                        std::string("DVCS_") + period_dir + "_" +
                        topo_dir(static_cast<Topology>(topo_idx));

                    const TopoSigma& ts = find_topo_sigma(topo_key);

                    if (!passes_sigma_cuts(ts, /*is_mc=*/false, b)) {
                        continue;
                    }

                    int current = 0;
                    if (resolve_current_for_label(label, run, current)) {
                        totals.data_by_period_current[label][current] += 1;
                    } else {
                        totals.data_uncategorized[label] += 1;
                        totals.data_uncategorized_runs[label].insert(run);
                    }

                    kept++;
                }

                std::cout << "[yield_totals]   DATA kept " << kept
                          << " events after cuts for " << label << std::endl;
            }
        }

        // ---------------- MC (reconstructed DVCS) ----------------
        {
            const std::string rec_key = tree_key + "_rec";
            auto itT = dvcsRecMcTrees.find(rec_key);
            if (itT != dvcsRecMcTrees.end() && itT->second) {
                TTree* t = itT->second;
                std::cout << "[yield_totals] Processing MC for period " << label
                          << " (tree key " << rec_key << ") with "
                          << static_cast<long long>(t->GetEntries()) << " entries\n";

                BranchBinder b;
                b.bind(t, /*require_runnum=*/false);

                const Long64_t N = t->GetEntries();
                long long kept = 0;

                for (Long64_t i = 0; i < N; ++i) {
                    if (t->GetEntry(i) <= 0) continue;

                    if (!passes_global(/*is_mc=*/true,
                                       /*runnum=*/0,
                                       period_label_internal,
                                       b,
                                       globalCuts)) {
                        continue;
                    }

                    const int topo_idx = b.topo_index();
                    if (topo_idx < 0 || topo_idx > 2) continue;

                    const std::string topo_key =
                        std::string("DVCS_") + period_dir + "_" +
                        topo_dir(static_cast<Topology>(topo_idx));

                    const TopoSigma& ts = find_topo_sigma(topo_key);

                    if (!passes_sigma_cuts(ts, /*is_mc=*/true, b)) {
                        continue;
                    }

                    totals.mc_by_period[label] += 1;
                    kept++;
                }

                std::cout << "[yield_totals]   MC kept " << kept
                          << " events after cuts for " << label << std::endl;
            }
        }
    }

    return totals;
}

// ---------------- pretty printing ----------------

static void write_totals_to_stream(std::ostream& os, const Totals& totals) {
    os << "================ Yield totals after exclusivity cuts ================\n\n";

    std::map<int, long long> global_by_current;

    os << "DATA totals by period and current (nA):\n";
    for (const auto& pk : totals.data_by_period_current) {
        const std::string& period = pk.first;
        const auto& by_curr = pk.second;

        os << "  Period: " << period << "\n";
        for (const auto& ck : by_curr) {
            int cur = ck.first;
            long long cnt = ck.second;
            os << "    current " << cur << " nA: " << cnt << "\n";
            global_by_current[cur] += cnt;
        }

        auto it_unc = totals.data_uncategorized.find(period);
        if (it_unc != totals.data_uncategorized.end() && it_unc->second > 0) {
            os << "    current (uncategorized): " << it_unc->second << "\n";
            auto it_runs = totals.data_uncategorized_runs.find(period);
            if (it_runs != totals.data_uncategorized_runs.end() &&
                !it_runs->second.empty()) {
                os << "    uncategorized runs:";
                for (int run : it_runs->second) {
                    os << " " << run;
                }
                os << "\n";
            }
        }
        os << "\n";
    }

    for (const auto& uk : totals.data_uncategorized) {
        const std::string& period = uk.first;
        if (totals.data_by_period_current.count(period) == 0 && uk.second > 0) {
            os << "  Period: " << period << "\n";
            os << "    current (uncategorized): " << uk.second << "\n";
            auto it_runs = totals.data_uncategorized_runs.find(period);
            if (it_runs != totals.data_uncategorized_runs.end() &&
                !it_runs->second.empty()) {
                os << "    uncategorized runs:";
                for (int run : it_runs->second) {
                    os << " " << run;
                }
                os << "\n";
            }
            os << "\n";
        }
    }

    os << "GLOBAL data totals by current (nA):\n";
    for (const auto& gk : global_by_current) {
        os << "  current " << gk.first << " nA: " << gk.second << "\n";
    }
    if (global_by_current.empty()) {
        os << "  [no data events counted]\n";
    }
    os << "\n";

    os << "MC totals by period (reconstructed DVCS, no current split):\n";
    for (const auto& mk : totals.mc_by_period) {
        os << "  Period: " << mk.first << "  ->  " << mk.second << " events\n";
    }
    if (totals.mc_by_period.empty()) {
        os << "  [no MC events counted]\n";
    }

    os << "\n=====================================================================\n";
}

} // end anon namespace

// ---------------- public API ----------------

bool compute_yield_totals(const std::map<std::string, TTree*>& dvcsDataTrees,
                          const std::map<std::string, TTree*>& dvcsRecMcTrees,
                          const std::string& combined_cuts_json,
                          const std::string& output_txt) {
    try {
        // Match other modules: ensure ROOT global locks are initialized.
        ROOT::EnableThreadSafety();

        const TopoCutMap sigmaCuts = load_sigma_cuts_both(combined_cuts_json);
        if (sigmaCuts.empty()) {
            std::cerr << "[yield_totals] ERROR: no sigma cuts loaded from "
                      << combined_cuts_json << std::endl;
            return false;
        }

        // Match other modules: use the global cuts config from default_global_cuts().
        const GlobalCutConfig& globalCuts = default_global_cuts();

        Totals totals = compute_totals_internal(dvcsDataTrees, dvcsRecMcTrees, sigmaCuts, globalCuts);

        {
            std::ofstream ofs(output_txt);
            if (!ofs.is_open()) {
                std::cerr << "[yield_totals] ERROR: cannot open output file "
                          << output_txt << " for writing.\n";
                return false;
            }
            write_totals_to_stream(ofs, totals);
            std::cout << "[yield_totals] Wrote summary to " << output_txt << std::endl;
        }

        write_totals_to_stream(std::cout, totals);
        return true;
    } catch (const std::exception& ex) {
        std::cerr << "[yield_totals] EXCEPTION: " << ex.what() << std::endl;
        return false;
    }
}