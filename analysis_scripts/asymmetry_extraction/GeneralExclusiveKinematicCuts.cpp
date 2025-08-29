#include "GeneralExclusiveKinematicCuts.h"
#include "common_vars.h"
#include <string>
#include <cmath>
#include "TMath.h"

// Extra headers for CSV-backed Mx2 window helper
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <limits>
#include <cctype>
#include <iostream>

using std::string;

// ====================================================================
// Tunables for dynamic Mx2 window (μ ± Nσ), CSV location, etc.
// ====================================================================
static constexpr double kNSigma = 2.0;                 // half-width in sigmas
static constexpr double kMx2Min = 0.70;                // hard clip low
static constexpr double kMx2Max = 1.10;                // hard clip high
static constexpr double kFallbackLo = 0.81;            // fallback window
static constexpr double kFallbackHi = 1.00;

static const char* kMx2CSVPath =
    "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/Mx2_fit_params.csv";

// xB slices used when the CSV was produced (must match exactly)
struct SliceDef { const char* name; double xa, xb; };
static const SliceDef kSlices[] = {
    {"Low",     0.10, 0.25},
    {"MidLow",  0.25, 0.35},
    {"MidHigh", 0.35, 0.45},
    {"High",    0.45, 0.60},
};
static constexpr int kNumSlices = 4;

// -t edges (ascending) used when the CSV was produced: [0.05, 0.25, ..., 1.25]
static const double kTEdgesPos[] = {0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25};
static constexpr int kNumTEdges = 7; // 6 bins

// ====================================================================
// CSV-backed Mx2 window provider
//   - Reads (slice_name, t_min, mu, sigma, fit_success) per cell
//   - Lookup by event (x, t) → (μ, σ) → [μ±Nσ], clipped to [kMx2Min, kMx2Max]
//   - Fallback to [kFallbackLo, kFallbackHi] if missing/failed/out-of-bounds
// ====================================================================
namespace {

struct Mx2Config {
    double nSigma     = kNSigma;
    double mx2Min     = kMx2Min;
    double mx2Max     = kMx2Max;
    double fallbackLo = kFallbackLo;
    double fallbackHi = kFallbackHi;
};

struct Key {
    std::string slice;
    int tbin;
    bool operator==(const Key& o) const { return slice == o.slice && tbin == o.tbin; }
};

struct KeyHash {
    std::size_t operator()(Key const& k) const noexcept {
        return std::hash<std::string>{}(k.slice) ^ (std::hash<int>{}(k.tbin) << 1);
    }
};

struct Row { double mu{std::numeric_limits<double>::quiet_NaN()};
             double sigma{std::numeric_limits<double>::quiet_NaN()};
             bool ok{false}; };

static std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b-1]))) --b;
    return s.substr(a, b - a);
}

static std::string sliceName(double x) {
    for (const auto& sl : kSlices) {
        if (x > sl.xa && x < sl.xb) return sl.name;
    }
    return "";
}

static int tbinIndex(double tpos) {
    // Return bin i such that tpos in [edge[i], edge[i+1]); -1 if outside
    if (!(tpos >= kTEdgesPos[0] && tpos < kTEdgesPos[kNumTEdges-1])) return -1;
    for (int i = 0; i < kNumTEdges - 1; ++i) {
        if (tpos >= kTEdgesPos[i] && tpos < kTEdgesPos[i+1]) return i;
    }
    return -1;
}

class Mx2WindowProvider {
public:
    explicit Mx2WindowProvider(const std::string& csvPath, Mx2Config cfg = {})
    : cfg_(cfg) { loaded_ = load(csvPath); }

    // Return (usedDynamic?, (lo, hi))
    std::pair<bool, std::pair<double,double>> windowFor(double x, double t) const {
        const std::string sname = sliceName(x);
        if (sname.empty()) {
            return {false, {cfg_.fallbackLo, cfg_.fallbackHi}};
        }
        const double tpos = -t;
        const int tb = tbinIndex(tpos);
        if (tb < 0) {
            return {false, {cfg_.fallbackLo, cfg_.fallbackHi}};
        }
        const auto it = table_.find(Key{sname, tb});
        if (it == table_.end() || !it->second.ok ||
            !(it->second.sigma > 0.0) ||
            !std::isfinite(it->second.mu) ||
            !std::isfinite(it->second.sigma)) {
            return {false, {cfg_.fallbackLo, cfg_.fallbackHi}};
        }
        double lo = it->second.mu - cfg_.nSigma * it->second.sigma;
        double hi = it->second.mu + cfg_.nSigma * it->second.sigma;
        // Clip to hard limits
        lo = std::max(lo, cfg_.mx2Min);
        hi = std::min(hi, cfg_.mx2Max);
        if (hi <= lo) {
            return {false, {cfg_.fallbackLo, cfg_.fallbackHi}};
        }
        return {true, {lo, hi}};
    }

    bool passes(double Mx2, double x, double t) const {
        const auto used_and_win = windowFor(x, t);
        const auto lo = used_and_win.second.first;
        const auto hi = used_and_win.second.second;
        return (Mx2 > lo) && (Mx2 < hi);
    }

    bool loaded() const { return loaded_; }

private:
    bool load(const std::string& path) {
        std::ifstream in(path);
        if (!in) {
            std::cerr << "[Mx2WindowProvider] WARNING: Could not open CSV: " << path
                      << " — will use fallback window.\n";
            return false;
        }
        std::string line;
        // Header
        if (!std::getline(in, line)) {
            std::cerr << "[Mx2WindowProvider] WARNING: Empty CSV: " << path
                      << " — will use fallback window.\n";
            return false;
        }
        // Rows
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string sname, x_min, x_max, t_min, t_max, mu, sigma, n_entries, fit_success;
            std::getline(ss, sname, ',');
            std::getline(ss, x_min, ',');
            std::getline(ss, x_max, ',');
            std::getline(ss, t_min, ',');
            std::getline(ss, t_max, ',');
            std::getline(ss, mu, ',');
            std::getline(ss, sigma, ',');
            std::getline(ss, n_entries, ',');
            std::getline(ss, fit_success, ',');

            const std::string s_trim = trim(sname);
            if (s_trim.empty()) continue;

            // Identify t-bin from t_min
            double tmin_val = std::numeric_limits<double>::quiet_NaN();
            if (!t_min.empty()) {
                try { tmin_val = std::stod(t_min); } catch (...) { tmin_val = std::numeric_limits<double>::quiet_NaN(); }
            }
            const int tb = tbinIndex(tmin_val);
            if (tb < 0) continue;

            Row r;
            try { r.mu    = mu.empty()    ? std::numeric_limits<double>::quiet_NaN() : std::stod(mu); } catch (...) { r.mu = std::numeric_limits<double>::quiet_NaN(); }
            try { r.sigma = sigma.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(sigma); } catch (...) { r.sigma = std::numeric_limits<double>::quiet_NaN(); }
            r.ok = (!fit_success.empty() && trim(fit_success) == "1" &&
                    std::isfinite(r.mu) && std::isfinite(r.sigma) && r.sigma > 0.0);
            table_.emplace(Key{s_trim, tb}, r);
        }
        return true;
    }

    Mx2Config cfg_;
    std::unordered_map<Key, Row, KeyHash> table_;
    bool loaded_{false};
};

} // namespace

// ====================================================================
// Physical masses (GeV)
// ====================================================================
static constexpr double m_e  = 0.000511;  // electron
static constexpr double m_pi = 0.13957;   // charged pion

//================================================================================
// Constructor: grab every branch we’ll need
//================================================================================
GeneralExclusiveKinematicCuts::GeneralExclusiveKinematicCuts(TTreeReader& reader)
    : BaseKinematicCuts(reader),
      runnum       (reader, "runnum"),
      fiducial_status(reader, "fiducial_status"),

      // Electron‐side branches (added e_p, e_theta)
      e_p          (reader, "e_p"),
      e_theta      (reader, "e_theta"),
      e_phi        (reader, "e_phi"),
      vz_e        (reader, "vz_e"),

      // Pion‐side branches (added p_theta)
      p_p          (reader, "p_p"),
      p_theta      (reader, "p_theta"),
      p_phi        (reader, "p_phi"),
      vz_p        (reader, "vz_p"),

      // Standard DIS / hadron variables
      Q2           (reader, "Q2"),
      W            (reader, "W"),
      Mx2          (reader, "Mx2"),
      xF           (reader, "xF"),
      pT           (reader, "pT"),
      y            (reader, "y"),
      x            (reader, "x"),
      xi           (reader, "xi"),
      phi          (reader, "phi"),
      z            (reader, "z"),
      t            (reader, "t"),
      tmin         (reader, "tmin"),
      tprime       (reader, "tprime"),
      target_pol   (reader, "target_pol")
{}

//================================================================================
// beamEnergy(run): (kept as in your code; default now returns 10.604 when
// out of the listed ranges to avoid zero.)
//================================================================================
static double beamEnergy(int run)
{
    if (run >= 6616  && run <= 6783)   return 10.1998;
    if (run >= 16042 && run <= 17065)  return 10.5473;
    if (run >= 17067 && run <= 17724)  return 10.5563;
    if (run >= 17725 && run <= 17811)  return 10.5593;
    return 10.604;
}

//================================================================================
// compute_t_scalar(…): kept for convenience; not used directly below.
//================================================================================
static double compute_t_scalar(int run,
                               double e_p, double e_theta, double e_phi,
                               double p_p, double p_theta, double p_phi)
{
    // 1) beam energy
    double Eb = beamEnergy(run);
    if (Eb <= 0.0) return 1e6; // invalid run → force fail

    // 2) scattered electron 4‐vector
    double E_e = std::sqrt(e_p*e_p + m_e*m_e);
    double sin_e = std::sin(e_theta);
    double cos_e = std::cos(e_theta);
    double ex = e_p * sin_e * std::cos(e_phi);
    double ey = e_p * sin_e * std::sin(e_phi);
    double ez = e_p * cos_e;

    // 3) pion 4‐vector
    double E_pi = std::sqrt(p_p*p_p + m_pi*m_pi);
    double sin_p = std::sin(p_theta);
    double cos_p = std::cos(p_theta);
    double px = p_p * sin_p * std::cos(p_phi);
    double py = p_p * sin_p * std::sin(p_phi);
    double pz = p_p * cos_p;

    // 4) virtual photon q = (Eb – E_e, –ex, –ey, Eb – ez)
    double E_q = Eb - E_e;
    double qx  = -ex;
    double qy  = -ey;
    double qz  = Eb - ez;

    // 5) Δ = q – p_pi
    double dE = E_q - E_pi;
    double dx = qx  - px;
    double dy = qy  - py;
    double dz = qz  - pz;

    // 6) t = (ΔE)^2 – (dx^2 + dy^2 + dz^2)
    return (dE*dE - (dx*dx + dy*dy + dz*dz));
}

//================================================================================
// applyCuts(…)
//   Now uses dynamic Mx2 window (μ ± Nσ) from CSV for the enpi* xB cases.
//   If the lookup fails, it falls back to [0.81, 1.00].
//================================================================================
bool GeneralExclusiveKinematicCuts::applyCuts(int currentFits, bool isMC)
{
    // Build the provider once (static lifetime) — thread-safe since C++11.
    static const Mx2WindowProvider s_mx2prov(kMx2CSVPath, Mx2Config{});

    // Basic naming lookup
    string property = binNames[currentFits];

    // Common DIS-level cuts
    if (*Q2 <  1.0    ) return false;
    if (*W  <  2.0    ) return false;
    if (*y  >  0.80   ) return false;


    // Helper to apply dynamic Mx2 window
    auto PassesDynamicMx2 = [&](double xval, double tval, double mx2val) -> bool {
        // Returns true if Mx2 is inside dynamic (or fallback) window for this (x, -t).
        return s_mx2prov.passes(mx2val, xval, tval);
    };

    // ----------------------------------------------------------------
    // xB-sliced properties: use dynamic Mx2 window from CSV for each event
    // ----------------------------------------------------------------
    if (*Q2 <  1.0    ) return false;
    if (*W  <  2.0    ) return false;
    if (*y  >  0.80   ) return false;
    if (*e_theta < 6 || *e_theta > 30) return false;
    if (*p_theta < 4 || *p_theta > 30) return false;
    if (std::abs(*vz_e + 2.2) > 5.0) return false;
    if (std::abs(*vz_e - *vz_p) > 7.0) return false;
    if (*Mx2  >  0.86 && *Mx2 < 1.02   ) return false;
    if (*x > 0.1 && *x < 0.60) return false;
    

    if (property == "enpi") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (*x > 0.10 && *x < 0.60);
        if (!goodEvent) return false;
        // return PassesDynamicMx2(*x, *t, *Mx2);
        return *Mx2 > 0.86 && *Mx2 < 1;
    }
    if (property == "enpi") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (*x > 0.10 && *x < 0.60);
        if (!goodEvent) return false;
        // return PassesDynamicMx2(*x, *t, *Mx2);
        return *Mx2 > 0.86 && *Mx2 < 1;
    }
    if (property == "enpiLowxB") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (*x > 0.10 && *x < 0.25);
        if (!goodEvent) return false;
        return PassesDynamicMx2(*x, *t, *Mx2);
    }
    if (property == "enpiMidLowxB") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (*x > 0.25 && *x < 0.35);
        if (!goodEvent) return false;
        return PassesDynamicMx2(*x, *t, *Mx2);
    }
    if (property == "enpiMidHighxB") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (*x > 0.35 && *x < 0.45);
        if (!goodEvent) return false;
        return PassesDynamicMx2(*x, *t, *Mx2);
    }
    if (property == "enpiHighxB") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (*x > 0.45 && *x < 0.60);
        if (!goodEvent) return false;
        return PassesDynamicMx2(*x, *t, *Mx2);
    }
    if (property == "enpiHarut1") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (-*tprime > 0.05 && -*tprime < 0.45);
        if (!goodEvent) return false;
        return *Mx2 > 0.86 && *Mx2 < 1;
        // return PassesDynamicMx2(*x, *tprime, *Mx2);
    }
    if (property == "enpiHarut2") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (-*t > 0.45 && -*t < 0.85);
        if (!goodEvent) return false;
        return *Mx2 > 0.86 && *Mx2 < 1;
        // return PassesDynamicMx2(*x, *tprime, *Mx2);
    }
    if (property == "enpiHarut3") {
        bool goodEvent = (*fiducial_status >= 111) &&
                         (-*tprime > 0.85 && -*tprime < 1.225);
        if (!goodEvent) return false;
        return PassesDynamicMx2(*x, *tprime, *Mx2);
    }

    return false;
}