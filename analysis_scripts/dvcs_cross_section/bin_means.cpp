#include "bin_means.h"

#include <TTree.h>

#include <algorithm>
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

// ---------------- configuration ----------------
static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;

// Simple kinematic cuts, matching exclusivity_cuts.cpp global cuts:
static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

// Topology: (FD,FD)= (1,1), (CD,FD)= (2,1), (CD,FT)= (2,0)
static inline bool passesTopology_simple(int detector1, int detector2, const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return (detector1 == 1 && detector2 == 1);
    if (topoStr == "(CD,FD)") return (detector1 == 2 && detector2 == 1);
    if (topoStr == "(CD,FT)") return (detector1 == 2 && detector2 == 0);
    return false;
}

// ---------------- helpers ----------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}

// Period like "DVCS_Fa18_inb" -> run tag key "fa18_inb"
// "DVCS_Fa18_inb_supp" -> "fa18_inb_supp"
// "DVCS_Sp18_out" -> "sp18_out", etc.
static std::string periodToRunTagKey(const std::string& period) {
    // Keep substring after first underscore
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    std::string rest = period.substr(pos + 1);  // e.g., "Fa18_inb"
    return toLower(rest);                        // "fa18_inb"
}

// wrap angle to [0, 2π)
static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return 0.0;
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    // if it lands exactly on 2π due to numeric, nudge it back
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}

// Given edges [0, 2π], with N_PHI_BINS equal-width bins
static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1; // guard 2π edge
    return idx;
}

// Build unique (min,max) bins from the scheme
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

// Find bin index for value within vector of [low,high) ranges
static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i) {
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    }
    return -1;
}

// ---------------- branch binder ----------------
// Uses doubles/ints like exclusivity code, but checks branch existence before binding.
struct BranchBinderBM {
    // topology helpers
    int    detector1 = 0; bool has_detector1 = false;
    int    detector2 = 0; bool has_detector2 = false;

    // kinematic cuts
    double t1 = 0.0;               bool has_t1 = false;
    double open_angle_ep2 = 0.0;   bool has_open_angle_ep2 = false;
    double pTmiss = 0.0;           bool has_pTmiss = false;

    // binning variables
    double x = 0.0;                bool has_x = false;       // Bjorken x
    double Q2 = 0.0;               bool has_Q2 = false;
    double phi2 = 0.0;             bool has_phi2 = false;
    double Delta_phi = 0.0;        bool has_Delta_phi = false;

    void bind(TTree* t) {
        if (!t) return;

        auto bindI = [&](const char* n, int* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };
        auto bindD = [&](const char* n, double* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; }
        };

        bindI("detector1", &detector1, has_detector1);
        bindI("detector2", &detector2, has_detector2);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bindD("pTmiss", &pTmiss, has_pTmiss);

        // binning vars
        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        // accept either phi2 or Delta_phi
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);
    }

    // Pick phi from phi2 if available, otherwise Delta_phi.
    double phi() const {
        if (has_phi2) return phi2;
        if (has_Delta_phi) return Delta_phi;
        return std::numeric_limits<double>::quiet_NaN();
    }

    bool readyForCuts() const {
        return has_detector1 && has_detector2 && has_t1 && has_open_angle_ep2 && has_pTmiss;
    }

    bool readyForBinning() const {
        return has_x && has_Q2 && (has_phi2 || has_Delta_phi);
    }
};

// ---------------- main routine ----------------
struct Accum {
    double sum_xB  = 0.0;
    double sum_Q2  = 0.0;
    double sum_t   = 0.0; // |t|
    double sum_phi = 0.0;
    long long count = 0;
};

void calculate_bin_means(
    const std::vector<std::string>& dvcs_periods,
    const std::vector<std::string>& topologies,
    const std::string& analysis_type,
    const std::vector<Binning>& binning_scheme,
    const std::string& output_json,
    const std::map<std::string, TTree*>& dataTrees
) {
    (void)analysis_type; // not used, kept for API parity

    std::cout << "[bin_means] Combining periods:";
    for (const auto& p : dvcs_periods) std::cout << " " << p;
    std::cout << " \n[bin_means]   Topologies:";
    for (const auto& t : topologies) std::cout << " " << t;
    std::cout << " \n" << std::endl;

    // Build unique ranges in each dimension
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[bin_means][ERROR] Missing xB/Q2/t bins from scheme. Aborting." << std::endl;
        return;
    }

    // Accumulator map keyed by (i_xB, i_Q2, i_t, i_phi)
    std::map<std::tuple<int,int,int,int>, Accum> acc;

    // Loop over periods & topologies
    for (const auto& period : dvcs_periods) {
        const std::string key = periodToRunTagKey(period); // e.g. "fa18_inb"
        auto itTree = dataTrees.find(key);
        if (itTree == dataTrees.end() || !itTree->second) {
            std::cerr << "[bin_means][WARN] No data tree found for period '" << period
                      << "' (key='" << key << "'). Skipping." << std::endl;
            continue;
        }
        TTree* t = itTree->second;

        // Bind branches once per tree
        BranchBinderBM b;
        b.bind(t);

        if (!b.readyForCuts()) {
            std::cerr << "[bin_means][WARN] Tree for '" << period
                      << "' lacks required cut branches (detector1/2, t1, open_angle_ep2, pTmiss). Skipping." << std::endl;
            continue;
        }
        if (!b.readyForBinning()) {
            std::cerr << "[bin_means][WARN] Tree for '" << period
                      << "' lacks required binning branches (x, Q2, phi2/Delta_phi). Skipping." << std::endl;
            continue;
        }

        const Long64_t nent = t->GetEntries();
        for (const auto& topo : topologies) {
            // event loop (single pass; we still need to read all entries per topology)
            t->ResetBranchAddresses(); // ensure consistent binding each pass
            b.bind(t);

            for (Long64_t i = 0; i < nent; ++i) {
                t->GetEntry(i);

                if (!passesTopology_simple(b.detector1, b.detector2, topo)) continue;
                if (!applyKinematicCuts_simple(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

                // Get values for binning
                double xB  = b.x;
                double Q2  = b.Q2;
                double tt  = std::fabs(b.t1);
                double phi = b.phi();
                if (!std::isfinite(xB) || !std::isfinite(Q2) || !std::isfinite(tt) || !std::isfinite(phi)) continue;

                int ix  = findBin(xB, xB_bins);
                int iQ2 = findBin(Q2, Q2_bins);
                int it  = findBin(tt,  t_bins);
                int ip  = phiToBin(phi);

                if (ix < 0 || iQ2 < 0 || it < 0 || ip < 0 || ip >= N_PHI_BINS) continue;

                auto key4 = std::make_tuple(ix, iQ2, it, ip);
                auto& a = acc[key4];
                a.sum_xB  += xB;
                a.sum_Q2  += Q2;
                a.sum_t   += tt;
                a.sum_phi += wrapToTwoPi(phi);
                a.count   += 1;
            }
        }
    }

    // Compute means, serialize JSON
    std::ofstream ofs(output_json);
    if (!ofs) {
        std::cerr << "[bin_means][ERROR] Cannot open output JSON: " << output_json << std::endl;
        return;
    }

    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";

    bool firstRec = true;
    for (const auto& kv : acc) {
        const auto& key = kv.first;
        const auto& a   = kv.second;
        if (a.count <= 0) continue;

        int ix, iQ2, it, ip;
        std::tie(ix, iQ2, it, ip) = key;

        double xB_avg  = a.sum_xB  / static_cast<double>(a.count);
        double Q2_avg  = a.sum_Q2  / static_cast<double>(a.count);
        double t_avg   = a.sum_t   / static_cast<double>(a.count);
        double phi_avg = a.sum_phi / static_cast<double>(a.count);

        if (!firstRec) ofs << ",\n";
        firstRec = false;
        ofs << "  \"(" << ix << "," << iQ2 << "," << it << "," << ip << ")\": "
            << "{\"xB_avg\":"  << xB_avg
            << ",\"Q2_avg\":"  << Q2_avg
            << ",\"t_avg\":"   << t_avg
            << ",\"phi_avg\":" << phi_avg
            << "}";
    }

    ofs << "\n}\n";
    ofs.close();

    std::cout << "[bin_means] Saved GLOBAL bin-averaged kinematics to " << output_json << std::endl;
}