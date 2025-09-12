#include "bin_means.h"

#include <TTree.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// ----- constants -----
constexpr int N_PHI_BINS = 12;
constexpr double TWO_PI = 2.0 * M_PI;

// ----- helpers -----

static inline std::string to_lower(std::string s) {
    for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Normalize incoming period labels to keys expected in dvcsDataTrees.
// Examples:
//   "DVCS_Fa18_inb"           -> "fa18_inb"
//   "DVCS_Fa18_out"           -> "fa18_out"
//   "DVCS_Sp19_inb"           -> "sp19_inb"
//   "DVCS_Sp18_out"           -> "sp18_out"
//   "DVCS_Fa18_inb_supp"      -> "fa18_inb_supp"
//   "DVCS_Fa18_inb_supplemental" -> "fa18_inb_supp"
static std::string period_to_tag(std::string s) {
    s = to_lower(std::move(s));
    if (s.rfind("dvcs_", 0) == 0) s.erase(0, 5);  // strip "dvcs_" prefix if present

    // normalize separators
    std::replace(s.begin(), s.end(), '-', '_');

    // normalize supplemental spellings
    if (s.find("supplemental") != std::string::npos || s.find("supplement") != std::string::npos)
        s = "fa18_inb_supp";
    else if (s == "fa18_inb_supp") /* ok */;
    // common canonical tags likely in your TTree map:
    // fa18_inb, fa18_out, sp19_inb, sp18_out, sp18_inb, etc.

    return s;
}

// Parse topology label to detector codes.
// FD=1, CD=2, FT=0 (matches your earlier exclusivity code).
static bool topo_match(int det1, int det2, const std::string& topoLabel) {
    if (topoLabel == "(FD,FD)") return (det1 == 1 && det2 == 1);
    if (topoLabel == "(CD,FD)") return (det1 == 2 && det2 == 1);
    if (topoLabel == "(CD,FT)") return (det1 == 2 && det2 == 0);
    return false;
}

// Global kinematic cuts used in your Python (apply_kinematic_cuts analogue):
// 1) open_angle_ep2 > 5.0 deg
// 2) (-t1) <= 1.0  (equivalently |t1| <= 1.0 if t1 is negative)
// 3) pTmiss <= 0.20
static bool apply_global_cuts(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if (std::fabs(t1) > 1.0)   return false;
    if (pTmiss > 0.20)         return false;
    return true;
}

// Ensure phi in [0, 2π) and return phi-bin index in [0, N_PHI_BINS-1].
static int phi_bin(double phi) {
    double x = std::fmod(phi, TWO_PI);
    if (x < 0.0) x += TWO_PI;
    const double w = TWO_PI / static_cast<double>(N_PHI_BINS);
    int b = static_cast<int>(std::floor(x / w));
    if (b < 0) b = 0;
    if (b >= N_PHI_BINS) b = N_PHI_BINS - 1;
    return b;
}

// Build unique bins (sorted by low,high).
template <typename GetPair>
static std::vector<std::pair<double,double>> unique_bins(const std::vector<Binning>& scheme,
                                                         GetPair getter) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) s.insert(getter(b));
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

// Return index of bin [low,high) that contains value; else -1.
static int find_bin(double value, const std::vector<std::pair<double,double>>& edges) {
    for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
        if (value >= edges[i].first && value < edges[i].second) return i;
    }
    return -1;
}

// Minimal branch binder for needed quantities.
struct Binder {
    // topology
    int detector1 = 0, detector2 = 0; bool has_d1 = false, has_d2 = false;
    // kinematics
    double t1 = 0.0, open_angle_ep2 = 0.0, pTmiss = 0.0;
    bool has_t1 = false, has_oa = false, has_pTmiss = false;
    // averaging targets
    double x = 0.0, Q2 = 0.0;
    bool has_x = false, has_Q2 = false;
    double phi2 = 0.0; bool has_phi2 = false;
    double Delta_phi = 0.0; bool has_Delta_phi = false;

    void bind(TTree* t) {
        auto bindI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bindD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bindI("detector1", &detector1, has_d1);
        bindI("detector2", &detector2, has_d2);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_oa);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);

        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);
    }
    bool has_required_for_cuts() const {
        return has_d1 && has_d2 && has_t1 && has_oa && has_pTmiss;
    }
    bool has_required_for_means() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta_phi);
    }
    double phi_value() const {
        return has_phi2 ? phi2 : (has_Delta_phi ? Delta_phi : 0.0);
    }
};

struct Accum {
    double sum_xB  = 0.0;
    double sum_Q2  = 0.0;
    double sum_t   = 0.0;  // |t|
    double sum_phi = 0.0;
    long long count = 0;
};

static std::string key_tuple_str(int ix, int iQ2, int it, int iPhi) {
    std::ostringstream os;
    os << "[" << ix << "," << iQ2 << "," << it << "," << iPhi << "]";
    return os.str();
}

} // namespace

void calculate_bin_means(const std::vector<std::string>& dvcs_periods,
                         const std::vector<std::string>& topologies,
                         const std::string& /*analysis_type*/,
                         const std::vector<Binning>& binning_scheme,
                         const std::string& output_json,
                         const std::map<std::string, TTree*>& dvcsDataTrees)
{
    std::cout << "[bin_means] Combining periods: ";
    for (const auto& p : dvcs_periods) std::cout << p << " ";
    std::cout << "\n[bin_means]   Topologies: ";
    for (const auto& t : topologies) std::cout << t << " ";
    std::cout << std::endl;

    // Unique bin edges for xB, Q2, |t|
    auto xB_bins = unique_bins(binning_scheme, [](const Binning& b){ return std::make_pair(b.xBmin, b.xBmax); });
    auto Q2_bins = unique_bins(binning_scheme, [](const Binning& b){ return std::make_pair(b.Q2min, b.Q2max); });
    auto t_bins  = unique_bins(binning_scheme, [](const Binning& b){ return std::make_pair(b.tmin,  b.tmax ); });

    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[bin_means] ERROR: empty bin edges (xB/Q2/t). Check binning scheme.\n";
        return;
    }

    // Prepare accumulators
    std::map<std::tuple<int,int,int,int>, Accum> accum;

    // Loop over periods and topologies
    for (const auto& per_in : dvcs_periods) {
        const std::string tag = period_to_tag(per_in);

        auto itTree = dvcsDataTrees.find(tag);
        if (itTree == dvcsDataTrees.end() || !itTree->second) {
            std::cerr << "[bin_means] WARNING: no DVCS data tree for tag \"" << tag
                      << "\" (from \"" << per_in << "\"). Skipping.\n";
            continue;
        }
        TTree* tree = itTree->second;

        Binder b;
        b.bind(tree);

        if (!b.has_required_for_cuts() || !b.has_required_for_means()) {
            std::cerr << "[bin_means] WARNING: required branches missing in tree \"" << tag
                      << "\"; needed: detector1/2, t1, open_angle_ep2, pTmiss, x, Q2, phi2 or Delta_phi.\n";
            continue;
        }

        const Long64_t n = tree->GetEntries();

        for (const auto& topo : topologies) {
            for (Long64_t i = 0; i < n; ++i) {
                tree->GetEntry(i);

                // topology filter
                if (!topo_match(b.detector1, b.detector2, topo)) continue;

                // global cuts
                if (!apply_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

                // bin indices
                const double xB  = b.x;
                const double Q2  = b.Q2;
                const double tt  = std::fabs(b.t1);
                const double phi = b.phi_value();

                const int ix  = find_bin(xB, xB_bins);
                const int iQ2 = find_bin(Q2, Q2_bins);
                const int it  = find_bin(tt,  t_bins);
                const int iP  = phi_bin(phi);

                if (ix < 0 || iQ2 < 0 || it < 0 || iP < 0) continue;

                auto key = std::make_tuple(ix, iQ2, it, iP);
                auto& a = accum[key];
                a.sum_xB  += xB;
                a.sum_Q2  += Q2;
                a.sum_t   += tt;
                a.sum_phi += phi;
                a.count   += 1;
            }
        }
    }

    // Write JSON
    std::ofstream ofs(output_json);
    if (!ofs) {
        std::cerr << "[bin_means] ERROR: cannot open output JSON: " << output_json << "\n";
        return;
    }

    ofs << "{\n";
    bool first = true;
    for (const auto& kv : accum) {
        const auto& key = kv.first;
        const auto& a   = kv.second;
        if (a.count <= 0) continue;

        const double xB_avg  = a.sum_xB  / static_cast<double>(a.count);
        const double Q2_avg  = a.sum_Q2  / static_cast<double>(a.count);
        const double t_avg   = a.sum_t   / static_cast<double>(a.count);
        const double phi_avg = a.sum_phi / static_cast<double>(a.count);

        if (!first) ofs << ",\n";
        first = false;

        int ix, iQ2, it, iP;
        std::tie(ix, iQ2, it, iP) = key;

        ofs << "  \"" << key_tuple_str(ix, iQ2, it, iP) << "\": {"
            << "\"xB_avg\":"  << std::setprecision(10) << xB_avg  << ","
            << "\"Q2_avg\":"  << std::setprecision(10) << Q2_avg  << ","
            << "\"t_avg\":"   << std::setprecision(10) << t_avg   << ","
            << "\"phi_avg\":" << std::setprecision(10) << phi_avg
            << "}";
    }
    ofs << "\n}\n";

    std::cout << "[bin_means] Wrote GLOBAL bin-averaged kinematics to " << output_json << std::endl;
}