#include "bin_means.h"

#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cctype>
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

// Simple kinematic cuts (same as exclusivity_cuts global)
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

static std::vector<std::pair<double,double>> uniqueRanges(
    const std::vector<Binning>& scheme, char which) 
{
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return {s.begin(), s.end()};
}

static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// Canonicalize period -> DVCS_<CapitalizedPeriod>
// e.g., "sp18_inb" -> "DVCS_Sp18_inb", "fa18_inb_supp" -> "DVCS_Fa18_inb_supp"
static std::string periodToDVCSKey(const std::string& period) {
    std::string s = period;
    // lower-case everything first
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    // capitalize only the very first character; keep everything after underscores as-is (lower-case)
    if (!s.empty() && std::isalpha(static_cast<unsigned char>(s[0]))) {
        s[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(s[0])));
    }
    return "DVCS_" + s;
}

// ---------------- branch binder ----------------
struct BranchBinderBM {
    int    detector1 = 0; bool has_detector1 = false;
    int    detector2 = 0; bool has_detector2 = false;

    double t1 = 0.0;               bool has_t1 = false;
    double open_angle_ep2 = 0.0;   bool has_open_angle_ep2 = false;
    double pTmiss = 0.0;           bool has_pTmiss = false;

    double x = 0.0;                bool has_x = false;
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

        bindD("x", &x, has_x);
        bindD("Q2", &Q2, has_Q2);
        bindD("phi2", &phi2, has_phi2);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);
    }

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
    double sum_t   = 0.0;
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
    (void)analysis_type;

    std::cout << "[bin_means] Using periods:";
    for (const auto& p : dvcs_periods) std::cout << " " << p;
    std::cout << "\n" << std::endl;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[bin_means][ERROR] Missing xB/Q2/t bins from scheme. Aborting." << std::endl;
        return;
    }

    std::map<std::tuple<int,int,int,int>, Accum> acc;

    // ---------------- Loop over periods ----------------
    for (const auto& period : dvcs_periods) {
        const std::string key = periodToDVCSKey(period); // e.g., "sp18_inb" -> "DVCS_Sp18_inb"

        auto itTree = dataTrees.find(key);
        if (itTree == dataTrees.end() || !itTree->second) {
            std::cerr << "[bin_means][FATAL] Missing data tree for period: '" << key << "'\n";
            std::cerr << "Available keys:";
            for (const auto& kv : dataTrees) std::cerr << " " << kv.first;
            std::cerr << std::endl;
            std::exit(EXIT_FAILURE);
        }

        TTree* t = itTree->second;
        BranchBinderBM b; b.bind(t);

        if (!b.readyForCuts()) {
            std::cerr << "[bin_means][FATAL] Tree for '" << key
                      << "' lacks cut branches (detector1/2,t1,open_angle_ep2,pTmiss)."
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (!b.readyForBinning()) {
            std::cerr << "[bin_means][FATAL] Tree for '" << key
                      << "' lacks binning vars (x,Q2,phi2/Delta_phi)." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        const Long64_t nent = t->GetEntries();
        for (const auto& topo : topologies) {
            for (Long64_t i = 0; i < nent; ++i) {
                t->GetEntry(i);
                if (!passesTopology_simple(b.detector1, b.detector2, topo)) continue;
                if (!applyKinematicCuts_simple(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

                double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;

                int ix=findBin(xB,xB_bins);
                int iQ=findBin(Q2,Q2_bins);
                int it=findBin(tt,t_bins);
                int ip=phiToBin(phi);
                if (ix<0||iQ<0||it<0||ip<0||ip>=N_PHI_BINS) continue;

                auto key4 = std::make_tuple(ix,iQ,it,ip);
                auto& a = acc[key4];
                a.sum_xB  += xB;
                a.sum_Q2  += Q2;
                a.sum_t   += tt;
                a.sum_phi += wrapToTwoPi(phi);
                a.count++;
            }
        }
    }

    // ---------------- Write JSON ----------------
    std::ofstream ofs(output_json);
    if (!ofs) {
        std::cerr << "[bin_means][FATAL] Cannot open output JSON: " << output_json << std::endl;
        std::exit(EXIT_FAILURE);
    }

    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    bool firstRec = true;
    for (const auto& kv : acc) {
        const auto& key = kv.first;
        const auto& a = kv.second;
        if (a.count <= 0) continue;

        int ix, iQ, it, ip;
        std::tie(ix, iQ, it, ip) = key;

        double xB_avg  = a.sum_xB  / a.count;
        double Q2_avg  = a.sum_Q2  / a.count;
        double t_avg   = a.sum_t   / a.count;
        double phi_avg = a.sum_phi / a.count;

        if (!firstRec) ofs << ",\n";
        firstRec = false;
        ofs << "  \"(" << ix << "," << iQ << "," << it << "," << ip << ")\": {"
            << "\"xB_avg\":"  << xB_avg
            << ",\"Q2_avg\":"  << Q2_avg
            << ",\"t_avg\":"   << t_avg
            << ",\"phi_avg\":" << phi_avg
            << "}";
    }
    ofs << "\n}\n";
    ofs.close();

    std::cout << "[bin_means] Saved global bin-averaged kinematics to " << output_json << std::endl;
}