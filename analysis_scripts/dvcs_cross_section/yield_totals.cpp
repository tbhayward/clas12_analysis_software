// yield_totals.cpp
// ------------------------------------------------------------
// Quick helper to count reconstructed events after exclusivity
// cuts, grouped by beam current (nA) for DATA and by period
// for MC.
//
// Data:
//   - Uses run-by-run current maps:
//       * Fa18 Inb: 40, 45, 50, 55 nA (explicit run lists +
//                    your extra runs mapped below).
//       * Fa18 Out: 5, 40, 50 nA (explicit run lists +
//                    your extra runs mapped below).
//       * Sp18 Out: 30 nA (3211-3293), 45 nA (3867-3987).
//       * Sp18 Inb: 35 nA (3306-3411),
//                   50 nA (3431-4325),
//                   plus explicit overrides for runs
//                   3418 (70 nA), 3421 (35 nA),
//                   3422 (35 nA), 3429 (50 nA).
//       * Sp19 Inb: all runs treated as 50 nA.
//   - Other periods are counted but placed in an "uncategorized"
//     bucket since we do not have a current map.
//
// MC:
//   - Counts per period (no current split).
//
// Cuts:
//   - Global kinematic cuts via passes_global_cuts(...) from
//     global_cuts.h (same as exclusivity_cuts.cpp).
//   - 3 sigma exclusivity cuts using combined_cuts.json entries
//     "data" for data and "mc" for MC, using topology keys
//     "DVCS_<PeriodDir>_<TopoDir>".
//
// Output:
//   - Prints to stdout and writes to a text file (output_txt).
//   - Also prints unique uncategorized run numbers per period.
// ------------------------------------------------------------

#include "yield_totals.h"

#include "periods.h"
#include "global_cuts.h"

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>

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

struct TopologyResolver {
    int detector1 = 0;
    int detector2 = 0;
    bool have1 = false;
    bool have2 = false;

    void enable_and_bind(TTree* t) {
        if (!t) fatal("TopologyResolver: null tree.");
        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("detector1", 1);
        t->SetBranchStatus("detector2", 1);
        have1 = (t->GetBranch("detector1") != nullptr);
        have2 = (t->GetBranch("detector2") != nullptr);
        if (!(have1 && have2)) {
            fatal("Missing detector1/detector2 in DVCS tree.");
        }
        t->SetBranchAddress("detector1", &detector1);
        t->SetBranchAddress("detector2", &detector2);
    }

    int index() const {
        if (detector1 == 1 && detector2 == 1) return static_cast<int>(Topology::FD_FD);
        if (detector1 == 2 && detector2 == 1) return static_cast<int>(Topology::CD_FD);
        if (detector1 == 2 && detector2 == 0) return static_cast<int>(Topology::CD_FT);
        return -1;
    }
};

// ---------------- branch binding ----------------

struct BranchBinding {
    std::string name;
    enum class Kind {
        kDouble, kFloat, kI32, kU32, kI64, kU64, kI16, kU16, kI8, kU8
    } kind = Kind::kDouble;

    double as_double = std::numeric_limits<double>::quiet_NaN();
    float  as_float  = std::numeric_limits<float>::quiet_NaN();
    int as_i32 = 0; unsigned as_u32 = 0;
    long long as_i64 = 0; unsigned long long as_u64 = 0;
    short as_i16 = 0; unsigned short as_u16 = 0;
    signed char as_i8 = 0; unsigned char as_u8 = 0;
};

static inline double bb_as_double(const BranchBinding& bb) {
    using K = BranchBinding::Kind;
    switch (bb.kind) {
        case K::kDouble: return bb.as_double;
        case K::kFloat:  return static_cast<double>(bb.as_float);
        case K::kI32:    return static_cast<double>(bb.as_i32);
        case K::kU32:    return static_cast<double>(bb.as_u32);
        case K::kI64:    return static_cast<double>(bb.as_i64);
        case K::kU64:    return static_cast<double>(bb.as_u64);
        case K::kI16:    return static_cast<double>(bb.as_i16);
        case K::kU16:    return static_cast<double>(bb.as_u16);
        case K::kI8:     return static_cast<double>(bb.as_i8);
        case K::kU8:     return static_cast<double>(bb.as_u8);
    }
    return std::numeric_limits<double>::quiet_NaN();
}

static inline long long bb_as_ll(const BranchBinding& bb) {
    using K = BranchBinding::Kind;
    switch (bb.kind) {
        case K::kDouble: return static_cast<long long>(std::llround(bb.as_double));
        case K::kFloat:  return static_cast<long long>(std::llround(bb.as_float));
        case K::kI32:    return static_cast<long long>(bb.as_i32);
        case K::kU32:    return static_cast<long long>(bb.as_u32);
        case K::kI64:    return static_cast<long long>(bb.as_i64);
        case K::kU64:    return static_cast<long long>(bb.as_u64);
        case K::kI16:    return static_cast<long long>(bb.as_i16);
        case K::kU16:    return static_cast<long long>(bb.as_u16);
        case K::kI8:     return static_cast<long long>(bb.as_i8);
        case K::kU8:     return static_cast<long long>(bb.as_u8);
    }
    return 0;
}

static std::mutex g_root_bind_mutex;

static void bind_one_exact_enable(TTree* t, const std::string& bname, BranchBinding& bb) {
    if (!t) fatal("bind_one_exact_enable: null tree for branch " + bname);
    std::lock_guard<std::mutex> lock(g_root_bind_mutex);
    t->SetBranchStatus(bname.c_str(), 1);
    TBranch* b = t->GetBranch(bname.c_str());
    if (!b) fatal("Branch not found: " + bname);
    TLeaf* leaf = b->GetLeaf(bname.c_str());
    if (!leaf) {
        leaf = static_cast<TLeaf*>(b->GetListOfLeaves()->First());
        if (!leaf) fatal("Branch has no leaves: " + bname);
    }
    const std::string tn = leaf->GetTypeName();
    if (tn == "Double_t" || tn == "double") {
        bb.kind = BranchBinding::Kind::kDouble;
        t->SetBranchAddress(bname.c_str(), &bb.as_double);
    } else if (tn == "Float_t" || tn == "float") {
        bb.kind = BranchBinding::Kind::kFloat;
        t->SetBranchAddress(bname.c_str(), &bb.as_float);
    } else if (tn == "Int_t" || tn == "int") {
        bb.kind = BranchBinding::Kind::kI32;
        t->SetBranchAddress(bname.c_str(), &bb.as_i32);
    } else if (tn == "UInt_t" || tn == "unsigned int") {
        bb.kind = BranchBinding::Kind::kU32;
        t->SetBranchAddress(bname.c_str(), &bb.as_u32);
    } else if (tn == "Long64_t" || tn == "long long") {
        bb.kind = BranchBinding::Kind::kI64;
        t->SetBranchAddress(bname.c_str(), &bb.as_i64);
    } else if (tn == "ULong64_t" || tn == "unsigned long long") {
        bb.kind = BranchBinding::Kind::kU64;
        t->SetBranchAddress(bname.c_str(), &bb.as_u64);
    } else if (tn == "Short_t" || tn == "short") {
        bb.kind = BranchBinding::Kind::kI16;
        t->SetBranchAddress(bname.c_str(), &bb.as_i16);
    } else if (tn == "UShort_t" || tn == "unsigned short") {
        bb.kind = BranchBinding::Kind::kU16;
        t->SetBranchAddress(bname.c_str(), &bb.as_u16);
    } else if (tn == "Char_t" || tn == "char" || tn == "signed char") {
        bb.kind = BranchBinding::Kind::kI8;
        t->SetBranchAddress(bname.c_str(), &bb.as_i8);
    } else if (tn == "UChar_t" || tn == "unsigned char") {
        bb.kind = BranchBinding::Kind::kU8;
        t->SetBranchAddress(bname.c_str(), &bb.as_u8);
    } else {
        fatal("Unsupported leaf type '" + tn + "' for branch '" + bname + "'");
    }
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

static inline bool within_3sigma(double val, const SigmaCut& sc) {
    if (!std::isfinite(val) || !std::isfinite(sc.mean) || !std::isfinite(sc.std) || sc.std <= 0.0) {
        return true;
    }
    return std::fabs(val - sc.mean) <= 3.0 * sc.std;
}

static bool passes_sigma_cuts(
    const TopoSigma& ts,
    bool is_mc,
    double Emiss2,
    double Mx2,
    double Mx2_1,
    double Mx2_2,
    double pTmiss,
    double xF,
    double theta_gg)
{
    const VarCutMap& m = is_mc ? ts.mc : ts.data;
    if (m.empty()) {
        return true;
    }

    for (const auto& kv : m) {
        const std::string& vname = kv.first;
        const SigmaCut& sc = kv.second;
        double val = std::numeric_limits<double>::quiet_NaN();

        if (vname == "Emiss2")                 val = Emiss2;
        else if (vname == "Mx2")               val = Mx2;
        else if (vname == "Mx2_1")             val = Mx2_1;
        else if (vname == "Mx2_2")             val = Mx2_2;
        else if (vname == "pTmiss")            val = pTmiss;
        else if (vname == "xF")                val = xF;
        else if (vname == "theta_gamma_gamma") val = theta_gg;
        else continue;

        if (!within_3sigma(val, sc)) return false;
    }

    return true;
}

// ---------------- run -> current maps ----------------

// Fa18 Inb
static int current_fa18_inb(int run, bool& ok) {
    static const std::unordered_map<int, int> m = {
        // 40 nA
        {5335, 40}, {5339, 40}, {5341, 40},
        {5340, 40}, {5342, 40}, {5343, 40}, {5344, 40}, {5345, 40}, // reclassified

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

        // Your extra Fa18 Inb 45 nA runs
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

        // 55 nA (5356-5407 are 55; 5400 included)
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

        // Your extra Fa18 Out 40 nA runs
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
        {5615, 50}, // redefined as 50 nA

        // Your extra Fa18 Out 50 nA runs
        {5505, 50}, {5567, 50}, {5617, 50}, {5621, 50}, {5623, 50},

        // Your extra Fa18 Out 5 nA run
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
static bool resolve_current_for_label(
    const std::string& period_label,
    int runnum,
    int& current_nA)
{
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
        // RGA Sp18 Out: 30 nA from 3211-3293, 45 nA from 3867-3987.
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
        // Explicit overrides first (your new info):
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

        // RGA Sp18 Inb ranges:
        //  - 35 nA from 3306-3411
        //  - 50 nA from 3431-4325
        if (runnum >= 3306 && runnum <= 3411) {
            current_nA = 35;
            return true;
        }
        if (runnum >= 3431 && runnum <= 4325) {
            current_nA = 50;
            return true;
        }

        // Other Sp18 Inb runs: unknown current (stay uncategorized)
        return false;
    }

    if (k == "sp19inb") {
        // All Sp19 Inb runs are 50 nA.
        current_nA = 50;
        return true;
    }

    // Other periods (Fa18 Inb Supp, etc.) -> uncategorized.
    return false;
}

// ---------------- core worker ----------------

struct Totals {
    std::map<std::string, std::map<int, long long>> data_by_period_current;
    std::map<std::string, long long> data_uncategorized;
    std::map<std::string, std::set<int>> data_uncategorized_runs;
    std::map<std::string, long long> mc_by_period;
};

static Totals compute_totals_internal(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const TopoCutMap& sigmaCuts)
{
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

        // Skip Fa18 Inb Supp entirely (no sigma cuts defined for this period)
        if (period_dir_for_label(label) == "Fa18_Inb_Supp") {
            std::cout << "[yield_totals] Skipping period " << label
                      << " (Fa18 Inb Supp; no sigma cuts defined)." << std::endl;
            continue;
        }

        // ---------------- DATA ----------------
        {
            auto itT = dvcsDataTrees.find(tree_key);
            if (itT != dvcsDataTrees.end() && itT->second) {
                TTree* t = itT->second;
                std::cout << "[yield_totals] Processing DATA for period " << label
                          << " (tree key " << tree_key << ") with "
                          << static_cast<long long>(t->GetEntries()) << " entries\n";

                TopologyResolver topo;
                topo.enable_and_bind(t);

                BranchBinding b_runnum, b_t1, b_open, b_pTmiss;
                BranchBinding b_Emiss2, b_Mx2, b_Mx2_1, b_Mx2_2, b_theta_gg, b_xF;

                bind_one_exact_enable(t, "runnum",            b_runnum);
                bind_one_exact_enable(t, "t1",                b_t1);
                bind_one_exact_enable(t, "open_angle_ep2",    b_open);
                bind_one_exact_enable(t, "pTmiss",            b_pTmiss);
                bind_one_exact_enable(t, "Emiss2",            b_Emiss2);
                bind_one_exact_enable(t, "Mx2",               b_Mx2);
                bind_one_exact_enable(t, "Mx2_1",             b_Mx2_1);
                bind_one_exact_enable(t, "Mx2_2",             b_Mx2_2);
                bind_one_exact_enable(t, "theta_gamma_gamma", b_theta_gg);
                bind_one_exact_enable(t, "xF",                b_xF);

                Long64_t N = t->GetEntries();
                long long kept = 0;

                for (Long64_t i = 0; i < N; ++i) {
                    if (t->GetEntry(i) <= 0) continue;

                    const double t1        = bb_as_double(b_t1);
                    const double open_deg  = bb_as_double(b_open);
                    const double pT        = bb_as_double(b_pTmiss);

                    if (!passes_global_cuts(t1, open_deg, pT)) continue;

                    const int topo_idx = topo.index();
                    if (topo_idx < 0 || topo_idx > 2) continue;

                    const std::string topo_key =
                        std::string("DVCS_") + period_dir_for_label(label) +
                        "_" + topo_dir(static_cast<Topology>(topo_idx));

                    const TopoSigma& ts = find_topo_sigma(topo_key);

                    const double Emiss2 = bb_as_double(b_Emiss2);
                    const double Mx2    = bb_as_double(b_Mx2);
                    const double Mx2_1  = bb_as_double(b_Mx2_1);
                    const double Mx2_2  = bb_as_double(b_Mx2_2);
                    const double xF     = bb_as_double(b_xF);
                    const double theta  = bb_as_double(b_theta_gg);

                    if (!passes_sigma_cuts(ts, /*is_mc=*/false,
                                           Emiss2, Mx2, Mx2_1, Mx2_2,
                                           pT, xF, theta)) {
                        continue;
                    }

                    const long long run_ll = bb_as_ll(b_runnum);
                    const int run = static_cast<int>(run_ll);

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

                TopologyResolver topo;
                topo.enable_and_bind(t);

                BranchBinding b_t1, b_open, b_pTmiss;
                BranchBinding b_Emiss2, b_Mx2, b_Mx2_1, b_Mx2_2, b_theta_gg, b_xF;

                bind_one_exact_enable(t, "t1",                b_t1);
                bind_one_exact_enable(t, "open_angle_ep2",    b_open);
                bind_one_exact_enable(t, "pTmiss",            b_pTmiss);
                bind_one_exact_enable(t, "Emiss2",            b_Emiss2);
                bind_one_exact_enable(t, "Mx2",               b_Mx2);
                bind_one_exact_enable(t, "Mx2_1",             b_Mx2_1);
                bind_one_exact_enable(t, "Mx2_2",             b_Mx2_2);
                bind_one_exact_enable(t, "theta_gamma_gamma", b_theta_gg);
                bind_one_exact_enable(t, "xF",                b_xF);

                Long64_t N = t->GetEntries();
                long long kept = 0;

                for (Long64_t i = 0; i < N; ++i) {
                    if (t->GetEntry(i) <= 0) continue;

                    const double t1        = bb_as_double(b_t1);
                    const double open_deg  = bb_as_double(b_open);
                    const double pT        = bb_as_double(b_pTmiss);

                    if (!passes_global_cuts(t1, open_deg, pT)) continue;

                    const int topo_idx = topo.index();
                    if (topo_idx < 0 || topo_idx > 2) continue;

                    const std::string topo_key =
                        std::string("DVCS_") + period_dir_for_label(label) +
                        "_" + topo_dir(static_cast<Topology>(topo_idx));

                    const TopoSigma& ts = find_topo_sigma(topo_key);

                    const double Emiss2 = bb_as_double(b_Emiss2);
                    const double Mx2    = bb_as_double(b_Mx2);
                    const double Mx2_1  = bb_as_double(b_Mx2_1);
                    const double Mx2_2  = bb_as_double(b_Mx2_2);
                    const double xF     = bb_as_double(b_xF);
                    const double theta  = bb_as_double(b_theta_gg);

                    if (!passes_sigma_cuts(ts, /*is_mc=*/true,
                                           Emiss2, Mx2, Mx2_1, Mx2_2,
                                           pT, xF, theta)) {
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

static void write_totals_to_stream(
    std::ostream& os,
    const Totals& totals)
{
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

bool compute_yield_totals(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& output_txt)
{
    try {
        const TopoCutMap sigmaCuts = load_sigma_cuts_both(combined_cuts_json);
        if (sigmaCuts.empty()) {
            std::cerr << "[yield_totals] ERROR: no sigma cuts loaded from "
                      << combined_cuts_json << std::endl;
            return false;
        }

        Totals totals = compute_totals_internal(dvcsDataTrees, dvcsRecMcTrees, sigmaCuts);

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