#include "total_counts.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ---------------- simple helpers ----------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}
static inline std::string trim(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, b - a + 1);
}

// Period like "DVCS_Fa18_inb" -> key "fa18_inb"; "DVCS_Sp18_out" -> "sp18_out".
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}

static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;

// wrap phi to [0, 2pi)
static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return 0.0;
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}
static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}

// Build unique (min,max) bins from the scheme for a coordinate (x, Q, t).
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}
static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// ------------- topology + universal kinematic cuts -------------
static inline bool passesTopology_simple(int detector1, int detector2, const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return (detector1 == 1 && detector2 == 1);
    if (topoStr == "(CD,FD)") return (detector1 == 2 && detector2 == 1);
    if (topoStr == "(CD,FT)") return (detector1 == 2 && detector2 == 0);
    return false;
}
static inline std::string topoToKey(const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return "FD_FD";
    if (topoStr == "(CD,FD)") return "CD_FD";
    if (topoStr == "(CD,FT)") return "CD_FT";
    return "FD_FD";
}
static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

// ------------- exclusivity 3sigma cuts: data-only stats -------------
struct Stats { double mean = 0.0; double std = 0.0; };
using VarCutMap = std::map<std::string, Stats>; // var -> {mean,std}
using PeriodTopoCuts = std::map<std::string, VarCutMap>; // "DVCS_Sp18_inb_FD_FD" -> data cuts

// We expect JSON written by exclusivity_cuts.cpp -> writeCombinedJson().
// Very small hand-rolled parser for that specific structure.
// Only extracts the "data" section for each key.
static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr << "[total_counts][ERROR] Cannot open cuts JSON: " << path << std::endl; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    // walk the file and find occurrences of "DVCS_..." keys
    size_t pos = 0;
    while (true) {
        size_t keyStart = s.find('"', pos);
        if (keyStart == std::string::npos) break;
        size_t keyEnd = s.find('"', keyStart + 1);
        if (keyEnd == std::string::npos) break;
        std::string key = s.substr(keyStart + 1, keyEnd - keyStart - 1); // e.g., DVCS_Sp18_inb_FD_FD

        // Find "data": { ... } that follows
        size_t dataPos = s.find("\"data\"", keyEnd);
        if (dataPos == std::string::npos) { pos = keyEnd + 1; continue; }
        size_t braceStart = s.find('{', dataPos);
        if (braceStart == std::string::npos) { pos = keyEnd + 1; continue; }

        // Match braces to extract the "data" object
        int depth = 0; size_t i = braceStart;
        for (; i < s.size(); ++i) {
            if (s[i] == '{') depth++;
            else if (s[i] == '}') { depth--; if (depth == 0) { ++i; break; } }
        }
        if (depth != 0) { pos = keyEnd + 1; continue; }
        std::string dataObj = s.substr(braceStart, i - braceStart); // includes braces

        // Parse entries like "Mx2":{"mean":..., "std":...}
        VarCutMap cuts;
        size_t vpos = 0;
        while (true) {
            size_t vKeyS = dataObj.find('"', vpos);
            if (vKeyS == std::string::npos) break;
            size_t vKeyE = dataObj.find('"', vKeyS + 1);
            if (vKeyE == std::string::npos) break;
            std::string var = dataObj.substr(vKeyS + 1, vKeyE - vKeyS - 1); // variable name

            size_t meanPos = dataObj.find("\"mean\"", vKeyE);
            size_t stdPos  = dataObj.find("\"std\"",  vKeyE);
            if (meanPos == std::string::npos || stdPos == std::string::npos) { vpos = vKeyE + 1; continue; }

            auto readNumberAfterColon = [&](size_t from)->double {
                size_t colon = dataObj.find(':', from);
                if (colon == std::string::npos) return 0.0;
                size_t j = colon + 1;
                while (j < dataObj.size() && std::isspace(static_cast<unsigned char>(dataObj[j]))) ++j;
                size_t k = j;
                while (k < dataObj.size() && (std::isdigit(static_cast<unsigned char>(dataObj[k])) || dataObj[k]=='-' || dataObj[k]=='+' || dataObj[k]=='.' || dataObj[k]=='e' || dataObj[k]=='E')) ++k;
                std::string num = dataObj.substr(j, k - j);
                try { return std::stod(num); } catch (...) { return 0.0; }
            };

            double m = readNumberAfterColon(meanPos);
            double sdev = readNumberAfterColon(stdPos);
            cuts[var] = Stats{m, sdev};

            vpos = vKeyE + 1;
        }

        // Only stash keys that look like DVCS_... (we are counting DVCS)
        if (key.rfind("DVCS_", 0) == 0) out[key] = cuts;

        pos = keyEnd + 1;
    }

    if (out.empty()) {
        std::cerr << "[total_counts][WARN] No DVCS cuts found in " << path << std::endl;
        return false;
    }
    return true;
}

static inline bool within3Sigma(double val, const Stats& s) {
    return (val >= s.mean - 3.0 * s.std) && (val <= s.mean + 3.0 * s.std);
}
static bool passes3SigmaCuts(const VarCutMap& cuts, const std::map<std::string,double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ------------- Branch binder (includes helicity and bin vars) -------------
struct BranchBinderTC {
    // topology helpers
    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    // kinematic cuts
    double t1 = 0.0;             bool has_t1 = false;
    double open_angle_ep2 = 0.0; bool has_open_angle_ep2 = false;
    double pTmiss = 0.0;         bool has_pTmiss = false;

    // binning variables
    double x = 0.0;     bool has_x = false;   // Bjorken x
    double Q2 = 0.0;    bool has_Q2 = false;
    double phi2 = 0.0;  bool has_phi2 = false;
    double Delta_phi = 0.0; bool has_Delta_phi = false;

    // exclusivity variables (for 3sigma checks)
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_1 = 0.0;             bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double xF = 0.0;                bool has_xF = false;

    int helicity = 0; bool has_helicity = false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI = [&](const char* n, int* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };
        auto bindD = [&](const char* n, double* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };

        bindI("detector1", &detector1, has_detector1);
        bindI("detector2", &detector2, has_detector2);
        bindI("helicity",  &helicity,  has_helicity);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);

        bindD("Emiss2", &Emiss2, has_Emiss2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("Mx2_1", &Mx2_1, has_Mx2_1);
        bindD("Mx2_2", &Mx2_2, has_Mx2_2);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bindD("xF", &xF, has_xF);
    }

    double phi() const { return has_phi2 ? phi2 : (has_Delta_phi ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    bool readyForCuts() const {
        return has_detector1 && has_detector2 && has_t1 && has_open_angle_ep2 && has_pTmiss && has_helicity;
    }
    bool readyForBinning() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta_phi);
    }

    std::map<std::string,double> valueMapForCuts() const {
        std::map<std::string,double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
        if (has_theta_gamma_gamma) m["theta_gamma_gamma"] = theta_gamma_gamma;
        if (has_pTmiss) m["pTmiss"] = pTmiss;
        if (has_xF) m["xF"] = xF;
        if (has_Emiss2) m["Emiss2"] = Emiss2;
        if (has_Mx2) m["Mx2"] = Mx2;
        if (has_Mx2_1) m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2) m["Mx2_2"] = Mx2_2;
        return m;
    }
};

// ------------- counting core -------------
struct HelCounts { long long plus = 0; long long minus = 0; };
using BinKey = std::tuple<int,int,int,int>; // (ix, iQ2, it, ip)

// group -> binKey -> HelCounts
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

static inline std::string tupleKeyStr(const BinKey& k) {
    int ix, iQ2, it, ip; std::tie(ix, iQ2, it, ip) = k;
    std::ostringstream os; os << "(" << ix << "," << iQ2 << "," << it << "," << ip << ")";
    return os.str();
}

// Map period "DVCS_Sp18_inb" -> pretty code "DVCS_Sp18_inb" (already pretty)
static inline std::string prettyPeriod(const std::string& p) { return p; }

// Map run-tag "sp18_inb" to group name for individual section output
static inline std::string indivGroupName(const std::string& runTag) { return runTag; }

// Combined group names
static inline std::string combFa18() { return "Fa18_inb_out_combined"; }
static inline std::string combFour() { return "Sp18_and_Fa18_inb_out_combined"; }

// Which individual run-tags participate in which combined groups
static inline bool inFa18(const std::string& runTag) { return (runTag == "fa18_inb" || runTag == "fa18_out"); }
static inline bool inFour(const std::string& runTag) { return (runTag == "sp18_inb" || runTag == "sp18_out" || runTag == "fa18_inb" || runTag == "fa18_out"); }

// JSON writer
static void writeCountsJson(const std::string& outPath,
                            const GroupCounts& groups,
                            int nPhi,
                            const std::vector<std::pair<double,double>>& xB_bins,
                            const std::vector<std::pair<double,double>>& Q2_bins,
                            const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(outPath);
    if (!ofs) { std::cerr << "[total_counts][ERROR] Cannot open output JSON: " << outPath << std::endl; return; }
    ofs << std::fixed << std::setprecision(8);

    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"groups\": {\n";

    bool firstGroup = true;
    for (const auto& gkv : groups) {
        if (!firstGroup) ofs << ",\n";
        firstGroup = false;
        ofs << "    \"" << gkv.first << "\": {\n";
        bool firstBin = true;
        for (const auto& bkv : gkv.second) {
            if (!firstBin) ofs << ",\n";
            firstBin = false;
            std::string keyStr = tupleKeyStr(bkv.first);
            const HelCounts& hc = bkv.second;
            long long total = hc.plus + hc.minus;
            ofs << "      \"" << keyStr << "\": {\"helicity\": {\"+1\": " << hc.plus
                << ", \"-1\": " << hc.minus << "}, \"total\": " << total << "}";
        }
        ofs << "\n    }";
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote " << outPath << std::endl;
}

void compute_total_counts(const std::vector<std::string>& periods,
                          const std::vector<std::string>& topologies,
                          const std::vector<Binning>& binning_scheme,
                          const std::map<std::string, TTree*>& dataTrees,
                          const std::string& cuts_json_path,
                          const std::string& output_json_path)
{
    // Build bin edges
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[total_counts][ERROR] Missing xB/Q2/t bins. Aborting." << std::endl; return;
    }

    // Load exclusivity 3sigma cuts (data) from combined JSON
    PeriodTopoCuts cuts;
    bool ok = loadCombinedCuts(cuts_json_path, cuts);
    if (!ok) std::cerr << "[total_counts][WARN] Proceeding, but no cuts were loaded." << std::endl;

    // Set up output container
    GroupCounts outCounts;

    // Loop over each requested period
    for (const std::string& period : periods) {
        const std::string runTag = periodToRunTagKey(period); // e.g., "sp18_inb"
        auto itTree = dataTrees.find(runTag);
        if (itTree == dataTrees.end() || !itTree->second) {
            std::cerr << "[total_counts][WARN] No data tree for '" << period << "' (key '" << runTag << "'). Skipping." << std::endl;
            continue;
        }
        TTree* t = itTree->second;

        // Bind branches once per tree
        BranchBinderTC b;
        b.bind(t);
        if (!b.readyForCuts() || !b.readyForBinning()) {
            std::cerr << "[total_counts][WARN] Missing required branches for '" << period << "'. Skipping." << std::endl;
            continue;
        }

        const Long64_t nent = t->GetEntries();
        // We will aggregate over all topologies but apply the topology-specific cuts.
        // Determine the DVCS key prefix "DVCS_<PrettyRun>" that exclusivity used.
        // Our combined JSON keys look like "DVCS_Sp18_inb_FD_FD" etc.
        // Build "DVCS_" + Capitalized runTag like exclusivity did.
        auto cap = runTag;
        if (!cap.empty()) cap[0] = std::toupper(cap[0]);
        for (size_t i = 0; i + 1 < cap.size(); ++i) if (cap[i] == '_' && i + 1 < cap.size()) cap[i + 1] = std::toupper(cap[i + 1]);
        const std::string periodCode = std::string("DVCS_") + cap;

        // Event loop once; per event apply the appropriate topology set.
        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);

            if (!applyKinematicCuts_simple(b.t1, b.open_angle_ep2, b.pTmiss)) continue;
            if (b.helicity != +1 && b.helicity != -1) continue;

            // Work out which topology this event is in (first match wins)
            std::string usedTopoKey;
            for (const auto& topoStr : topologies) {
                if (passesTopology_simple(b.detector1, b.detector2, topoStr)) { usedTopoKey = topoToKey(topoStr); break; }
            }
            if (usedTopoKey.empty()) continue;

            // Load cuts for this period+topology if available
            VarCutMap topoCuts;
            auto itCut = cuts.find(periodCode + "_" + usedTopoKey);
            if (itCut != cuts.end()) topoCuts = itCut->second; // else empty => no extra cut

            // Apply 3sigma cuts
            if (!topoCuts.empty()) {
                auto vals = b.valueMapForCuts();
                if (!passes3SigmaCuts(topoCuts, vals)) continue;
            }

            // Binning
            double xB = b.x, Q2 = b.Q2, tt = std::fabs(b.t1), phi = b.phi();
            if (!std::isfinite(xB) || !std::isfinite(Q2) || !std::isfinite(tt) || !std::isfinite(phi)) continue;

            int ix  = findBin(xB, xB_bins);
            int iQ2 = findBin(Q2, Q2_bins);
            int itb = findBin(tt,  t_bins);
            int ip  = phiToBin(phi);
            if (ix < 0 || iQ2 < 0 || itb < 0 || ip < 0 || ip >= N_PHI_BINS) continue;

            BinKey key(ix, iQ2, itb, ip);

            // Individual group
            HelCounts& hc_ind = outCounts[indivGroupName(runTag)][key];
            if (b.helicity == +1) hc_ind.plus++; else hc_ind.minus++;

            // Fa18 combined
            if (inFa18(runTag)) {
                HelCounts& hc_fa = outCounts[combFa18()][key];
                if (b.helicity == +1) hc_fa.plus++; else hc_fa.minus++;
            }
            // Four-period combined
            if (inFour(runTag)) {
                HelCounts& hc_4 = outCounts[combFour()][key];
                if (b.helicity == +1) hc_4.plus++; else hc_4.minus++;
            }
        }

        std::cout << "[total_counts] Counted period " << period << " using cuts key '" << (periodCode) << "_*'" << std::endl;
    }

    // Serialize JSON
    writeCountsJson(output_json_path, outCounts, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
}