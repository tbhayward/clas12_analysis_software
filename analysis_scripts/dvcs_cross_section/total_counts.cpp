// total_counts.cpp
// ------------------------------------------------------------
// Compute helicity-separated total counts after exclusivity cuts.
// Strict mode: any missing lookup (trees, branches, dirs, files) causes a fatal exit.
// - Reads: DVCS data trees per period (must have at least: helicity, x, Q2, t1, phi2)
//          plus any branches referenced by the combined cuts JSON.
// - Writes:
//     • <out_root_dir>/jsons/total_counts.json             (master, nested "groups")
//     • <out_root_dir>/jsons/total_counts_<label>.json     (flat per-group file; label = fa18_inb, etc.)
// - Produces per-group phi-binned plots under <out_root_dir>/total_counts_plots/<label>/
// ------------------------------------------------------------

#include "total_counts.h"
#include "periods.h"  // canonical PeriodDef {label, tree_key}

#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLatex.h>
#include <TH1.h>
#include <TGraphErrors.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdarg>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

#include <nlohmann/json.hpp>
using nlohmann::json;

namespace {

// ===============================================
// Logging helpers
// ===============================================
static inline void logf(const char* tag, const char* fmt, ...) {
    std::fprintf(stderr, "[%s] ", tag);
    va_list ap; va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    va_end(ap);
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
}

[[noreturn]] static void fatalf(const char* fmt, ...) {
    std::fprintf(stderr, "[total_counts][FATAL] ");
    va_list ap; va_start(ap, fmt);
    std::vfprintf(stderr, fmt, ap);
    va_end(ap);
    std::fprintf(stderr, "\n");
    std::fflush(stderr);
    std::exit(EXIT_FAILURE);
}

// ===============================================
// Constants and small helpers
// ===============================================
static constexpr int N_PHI_BINS = 12;
static constexpr double TWO_PI  = 2.0 * M_PI;

static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}

static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    if (!std::isfinite(w)) return -1;
    const double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}

static std::vector<std::pair<double,double>>
uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return {s.begin(), s.end()};
}

static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < (int)ranges.size(); ++i) {
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    }
    return -1;
}

static inline std::string keyStr(int ix,int iQ,int it,int ip) {
    std::ostringstream os; os<<"("<<ix<<","<<iQ<<","<<it<<","<<ip<<")";
    return os.str();
}

struct HelCounts { long long plus=0, minus=0; };

// Canonical helpers
static inline bool is_canonical_tree_key(const std::string& k) {
    for (const auto& p : CANONICAL_PERIODS()) if (k == p.tree_key) return true;
    return false;
}

static inline const char* tree_key_to_label_or_fatal(const std::string& k) {
    for (const auto& p : CANONICAL_PERIODS()) if (k == p.tree_key) return p.label;
    std::ostringstream oss;
    oss << "Non-canonical tree key '" << k << "'. Expected one of:";
    for (const auto& p : CANONICAL_PERIODS()) oss << " " << p.tree_key;
    fatalf("%s", oss.str().c_str());
    return ""; // unreachable
}

// ===============================================
// Cut model
// ===============================================
struct SigmaCut {
    std::string branch;
    double center = 0.0;
    double sigma  = 0.0;
    double nsigma = 3.0;
    enum Mode { TWO_SIDED, UPPER, LOWER } mode = TWO_SIDED;
};

struct BaseCut {
    std::string branch;
    bool has_min = false;
    bool has_max = false;
    bool has_eq  = false;
    bool has_neq = false;
    double vmin = 0.0;
    double vmax = 0.0;
    long long eq = 0;
    long long neq = 0;
};

struct PeriodCuts {
    std::vector<SigmaCut> sigma_cuts;
    std::vector<BaseCut>  base_cuts;
};

static SigmaCut::Mode parseMode(const std::string& m) {
    if (m == "two_sided") return SigmaCut::TWO_SIDED;
    if (m == "upper")     return SigmaCut::UPPER;
    if (m == "lower")     return SigmaCut::LOWER;
    fatalf("Unknown sigma cut mode: %s", m.c_str());
    return SigmaCut::TWO_SIDED; // unreachable
}

static void parseSigmaCuts(const json& arr, std::vector<SigmaCut>& out) {
    if (!arr.is_array()) fatalf("sigma_cuts must be an array");
    for (const auto& j : arr) {
        SigmaCut c;
        if (!j.contains("branch")) fatalf("sigma_cuts entry missing 'branch'");
        c.branch = j.at("branch").get<std::string>();
        if (!j.contains("center") || !j.contains("sigma"))
            fatalf("sigma_cuts entry for '%s' missing 'center' or 'sigma'", c.branch.c_str());
        c.center  = j.at("center").get<double>();
        c.sigma   = j.at("sigma").get<double>();
        c.nsigma  = j.value("nsigma", 3.0);
        c.mode    = parseMode(j.value("mode", std::string("two_sided")));
        out.push_back(c);
    } //endfor
}

static void parseBaseCuts(const json& arr, std::vector<BaseCut>& out) {
    if (!arr.is_array()) fatalf("base_cuts must be an array");
    for (const auto& j : arr) {
        BaseCut c;
        if (!j.contains("branch")) fatalf("base_cuts entry missing 'branch'");
        c.branch = j.at("branch").get<std::string>();

        if (j.contains("min")) { c.has_min = true; c.vmin = j.at("min").get<double>(); }
        if (j.contains("max")) { c.has_max = true; c.vmax = j.at("max").get<double>(); }

        if (j.contains("eq"))  { c.has_eq = true;  c.eq  = j.at("eq").get<long long>(); }
        if (j.contains("neq")) { c.has_neq = true; c.neq = j.at("neq").get<long long>(); }

        if (!(c.has_min || c.has_max || c.has_eq || c.has_neq)) {
            fatalf("base_cuts entry for '%s' has no condition (min/max/eq/neq)", c.branch.c_str());
        }
        out.push_back(c);
    } //endfor
}

static PeriodCuts loadCombinedCuts(const std::string& json_path, const std::string& period_key) {
    logf("cuts", "Loading combined cuts: %s (period override key: %s)", json_path.c_str(), period_key.c_str());
    std::ifstream ifs(json_path);
    if (!ifs) fatalf("Cannot open combined cuts JSON: %s", json_path.c_str());
    json J; ifs >> J;

    PeriodCuts cuts;

    // global
    if (J.contains("global")) {
        const auto& G = J.at("global");
        if (G.contains("sigma_cuts")) parseSigmaCuts(G.at("sigma_cuts"), cuts.sigma_cuts);
        if (G.contains("base_cuts"))  parseBaseCuts(G.at("base_cuts"), cuts.base_cuts);
    }

    // period override (replace arrays if present)
    if (J.contains("period_overrides")) {
        const auto& PO = J.at("period_overrides");
        auto it = PO.find(period_key);
        if (it != PO.end()) {
            const auto& P = *it;
            if (P.contains("sigma_cuts")) {
                cuts.sigma_cuts.clear();
                parseSigmaCuts(P.at("sigma_cuts"), cuts.sigma_cuts);
            }
            if (P.contains("base_cuts")) {
                cuts.base_cuts.clear();
                parseBaseCuts(P.at("base_cuts"), cuts.base_cuts);
            }
            logf("cuts", "Applied period override for key: %s", period_key.c_str());
        }
    }

    // sanity
    for (const auto& c : cuts.sigma_cuts) {
        if (c.sigma <= 0.0) fatalf("sigma_cuts: non-positive sigma for branch '%s'", c.branch.c_str());
        if (c.nsigma <= 0.0) fatalf("sigma_cuts: non-positive nsigma for branch '%s'", c.branch.c_str());
    }

    // print a summary
    logf("cuts", "Summary: %zu sigma_cuts, %zu base_cuts", cuts.sigma_cuts.size(), cuts.base_cuts.size());
    for (const auto& c : cuts.sigma_cuts) {
        logf("cuts", "  sigma: branch=%s center=%.6g sigma=%.6g nsigma=%.3g mode=%d",
             c.branch.c_str(), c.center, c.sigma, c.nsigma, (int)c.mode);
    } //endfor
    for (const auto& c : cuts.base_cuts) {
        std::ostringstream oss; oss << "  base: branch=" << c.branch;
        if (c.has_min) oss << " min=" << c.vmin;
        if (c.has_max) oss << " max=" << c.vmax;
        if (c.has_eq)  oss << " eq="  << c.eq;
        if (c.has_neq) oss << " neq=" << c.neq;
        logf("cuts", "%s", oss.str().c_str());
    } //endfor

    return cuts;
}

// Collect unique list of branches required by cuts
static std::vector<std::string> requiredCutBranches(const PeriodCuts& cuts) {
    std::set<std::string> s;
    for (const auto& c : cuts.sigma_cuts) s.insert(c.branch);
    for (const auto& c : cuts.base_cuts)  s.insert(c.branch);
    return {s.begin(), s.end()};
}

// ===============================================
// Branch bindings
// ===============================================
struct BranchBinding {
    std::string name;
    double     as_double = std::numeric_limits<double>::quiet_NaN();
    long long  as_ll     = 0;
    bool       is_int    = false; // set at bind time
};

static bool isIntegerLeaf(TLeaf* leaf) {
    if (!leaf) return false;
    const char* t = leaf->GetTypeName();
    return std::string(t) == "Int_t"    || std::string(t) == "UInt_t"   ||
           std::string(t) == "Short_t"  || std::string(t) == "UShort_t" ||
           std::string(t) == "Char_t"   || std::string(t) == "UChar_t"  ||
           std::string(t) == "Long64_t" || std::string(t) == "ULong64_t";
}

static void bindRequiredBranches_STRICT(
    const std::string& period_key,
    TTree* t,
    const std::vector<std::string>& branch_names,
    std::unordered_map<std::string, BranchBinding>& bindings)
{
    static const std::set<std::string> RESERVED = {"helicity","x","Q2","t1","phi2"};
    bindings.clear();
    bindings.reserve(branch_names.size());

    logf("bind", "[%s] Begin binding %zu cut branches", period_key.c_str(), branch_names.size());
    for (const auto& bname : branch_names) {
        if (RESERVED.count(bname)) {
            logf("bind", "[%s]   SKIP reserved branch '%s' (already bound to local var)", period_key.c_str(), bname.c_str());
            continue;
        }

        TBranch* b = t->GetBranch(bname.c_str());
        if (!b) fatalf("[%s] Required branch for cuts missing: '%s'", period_key.c_str(), bname.c_str());

        TLeaf* leaf = b->GetLeaf(bname.c_str());
        if (!leaf) {
            leaf = (TLeaf*)b->GetListOfLeaves()->First();
            if (!leaf) fatalf("[%s] Branch has no leaves: '%s'", period_key.c_str(), bname.c_str());
        }

        auto [it, inserted] = bindings.emplace(bname, BranchBinding{});
        BranchBinding& bb = it->second;
        bb.name   = bname;
        bb.is_int = isIntegerLeaf(leaf);

        if (bb.is_int) {
            t->SetBranchAddress(bname.c_str(), &bb.as_ll);
            logf("bind", "[%s]   INT  '%s' -> as_ll @%p", period_key.c_str(), bname.c_str(), (void*)&bb.as_ll);
        } else {
            t->SetBranchAddress(bname.c_str(), &bb.as_double);
            logf("bind", "[%s]   REAL '%s' -> as_double @%p", period_key.c_str(), bname.c_str(), (void*)&bb.as_double);
        }
    } //endfor
    logf("bind", "[%s] Done binding cut branches", period_key.c_str());
}

// ===============================================
// Cut evaluation
// ===============================================
static inline bool passBaseCuts(const std::vector<BaseCut>& baseCuts,
                                const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : baseCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false; // strict
        const BranchBinding& bb = it->second;

        if (c.has_eq) {
            long long v = bb.is_int ? bb.as_ll : (long long)std::llround(bb.as_double);
            if (v != c.eq) return false;
        }
        if (c.has_neq) {
            long long v = bb.is_int ? bb.as_ll : (long long)std::llround(bb.as_double);
            if (v == c.neq) return false;
        }
        if (c.has_min) {
            double v = bb.is_int ? (double)bb.as_ll : bb.as_double;
            if (!(v >= c.vmin)) return false;
        }
        if (c.has_max) {
            double v = bb.is_int ? (double)bb.as_ll : bb.as_double;
            if (!(v <= c.vmax)) return false;
        }
    } //endfor
    return true;
}

static inline bool passSigmaCuts(const std::vector<SigmaCut>& sigmaCuts,
                                 const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : sigmaCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false; // strict
        const BranchBinding& bb = it->second;
        const double v = bb.is_int ? (double)bb.as_ll : bb.as_double;
        const double lo = c.center - c.nsigma * c.sigma;
        const double hi = c.center + c.nsigma * c.sigma;

        if (c.mode == SigmaCut::TWO_SIDED) {
            if (v < lo || v > hi) return false;
        } else if (c.mode == SigmaCut::UPPER) {
            if (v > hi) return false;
        } else { // LOWER
            if (v < lo) return false;
        }
    } //endfor
    return true;
}

static inline bool passes3SigmaCuts_STRICT(
    const PeriodCuts& cuts,
    const std::unordered_map<std::string, BranchBinding>& B)
{
    return passBaseCuts(cuts.base_cuts, B) && passSigmaCuts(cuts.sigma_cuts, B);
}

// ===============================================
// JSON writers
// ===============================================
static void write_total_counts_group_json(
    const std::string& out_path,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) fatalf("Cannot open for write: %s", out_path.c_str());
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size() << "},\n";
    ofs << "  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n"; first=false;
        const HelCounts& hc = kv.second;
        int ix,iq,itb,ip; std::tie(ix,iq,itb,ip) = kv.first;
        ofs << "    \"" << keyStr(ix,iq,itb,ip) << "\": {"
            << "\"helicity\":{\"+1\":" << hc.plus
            << ",\"-1\":" << hc.minus << "},"
            << "\"total\":" << (hc.plus + hc.minus)
            << "}";
    } //endfor
    ofs << "\n  }\n}\n";
    ofs.close();
    logf("json", "Wrote %s", out_path.c_str());
}

static void write_total_counts_master_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>>& allGroupsByLabel,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) fatalf("Cannot open for write: %s", out_path.c_str());
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\":  " << t_bins.size() << "},\n";
    ofs << "  \"groups\": {\n";

    bool firstG = true;
    for (const auto& gkv : allGroupsByLabel) {
        if (!firstG) ofs << ",\n";
        firstG = false;
        ofs << "    \"" << gkv.first << "\": { \"bins\": {\n";
        bool firstB = true;
        for (const auto& kv : gkv.second) {
            if (!firstB) ofs << ",\n"; firstB=false;
            const HelCounts& hc = kv.second;
            int ix,iq,itb,ip; std::tie(ix,iq,itb,ip)=kv.first;
            ofs << "      \"" << keyStr(ix,iq,itb,ip) << "\": {"
                << "\"helicity\":{\"+1\":" << hc.plus
                << ",\"-1\":" << hc.minus << "},"
                << "\"total\":" << (hc.plus + hc.minus)
                << "}";
        } //endfor
        ofs << "\n    }}";
    } //endfor
    ofs << "\n  }\n}\n";
    ofs.close();
    logf("json", "Wrote master %s", out_path.c_str());
}

// ===============================================
// Plotting
// ===============================================
static std::vector<double> phiCentersDeg() {
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / static_cast<double>(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) v[i] = (i + 0.5) * step;
    return v;
}

static bool pairAlmostEqual(const std::pair<double,double>& a,
                            const std::pair<double,double>& b,
                            double eps = 1e-9) {
    return std::fabs(a.first  - b.first)  < eps &&
           std::fabs(a.second - b.second) < eps;
}

static void plot_group_counts(
    const std::string& group_label,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::string& out_dir)
{
    namespace fs = std::filesystem;
    logf("plot", "[%s] Creating plots in: %s", group_label.c_str(), out_dir.c_str());

    static const std::vector<double> PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
    for (int i = 0; i < N_PHI_BINS; ++i) X[i] = PHI[i];

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) fatalf("[plot][%s] Cannot create directory: %s (%s)",
                   group_label.c_str(), out_dir.c_str(), ec.message().c_str());

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        logf("plot", "[%s] xB bin %d -> (%.6g, %.6g)", group_label.c_str(), ix,
             xB_bins[ix].first, xB_bins[ix].second);

        std::set<std::pair<double,double>> q2set, tset;
        for (const auto& b : binning_scheme) {
            if (std::make_pair(b.xBmin, b.xBmax) == xB_bins[ix]) {
                q2set.emplace(b.Q2min, b.Q2max);
                tset.emplace(b.tmin,  b.tmax);
            }
        }

        std::vector<std::pair<double,double>> Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double>> Ts (tset.begin(),  tset.end());

        if (Q2s.empty() || Ts.empty()) {
            logf("plot", "[%s]   (no Q2 or t bins under this xB range) — skipping", group_label.c_str());
            continue;
        }

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();

        const int W = 260 * ncols + 120;
        const int H = 220 * nrows + 140;

        TCanvas* c = new TCanvas(Form("c_counts_%s_xB%d", group_label.c_str(), ix), "", W, H);
        if (!c) fatalf("[plot][%s] TCanvas allocation failed", group_label.c_str());

        TPad* pTop  = new TPad("pTop",  "pTop",  0.0, 0.94, 1.0, 1.0);
        TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.94);
        pTop->SetFillStyle(0);  pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        if (!gPad) fatalf("[plot][%s] gPad null after pTop->cd()", group_label.c_str());
        {
            TLatex head; head.SetNDC(); head.SetTextSize(0.055); head.SetTextAlign(22);
            head.DrawLatex(0.5, 0.55,
                Form("%s   x_B (%.3g, %.3g)", group_label.c_str(), xB_bins[ix].first, xB_bins[ix].second));
        }

        const int cells = nrows * ncols;

        for (int r = 0; r < nrows; ++r) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int cell = r * ncols + ccol + 1;
                if (cell < 1 || cell > cells) continue;

                pGrid->cd(cell);
                if (!gPad) fatalf("[plot][%s] gPad null after pGrid->cd(cell)", group_label.c_str());

                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.06);

                // Map local (Q2,t) pair to global indices robustly
                int iQ = -1, itb = -1;
                for (int iq = 0; iq < (int)Q2_bins.size(); ++iq) {
                    if (pairAlmostEqual(Q2_bins[iq], Q2s[r])) { iQ = iq; break; }
                }
                for (int it = 0; it < (int)t_bins.size(); ++it) {
                    if (pairAlmostEqual(t_bins[it], Ts[ccol])) { itb = it; break; }
                }

                std::vector<double> Yp(N_PHI_BINS, 0.0), Ym(N_PHI_BINS, 0.0);
                std::vector<double> EYp(N_PHI_BINS, 0.0), EYm(N_PHI_BINS, 0.0);
                double local_max = 0.0;

                if (iQ >= 0 && itb >= 0) {
                    for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                        auto itbl = table.find({ix, iQ, itb, ip});
                        if (itbl == table.end()) continue;
                        Yp[ip]  = (double)itbl->second.plus;
                        Ym[ip]  = (double)itbl->second.minus;
                        EYp[ip] = std::sqrt(std::max(0.0, Yp[ip]));
                        EYm[ip] = std::sqrt(std::max(0.0, Ym[ip]));
                        local_max = std::max(local_max, std::max(Yp[ip], Ym[ip]));
                    }
                }

                double ymax = std::max(1.0, local_max * 1.15);
                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, ymax);
                if (!frame) fatalf("[plot][%s] DrawFrame returned null", group_label.c_str());
                frame->GetXaxis()->SetTitle("phi (deg)");
                frame->GetYaxis()->SetTitle("Counts");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);

                if (local_max >= 500.0) {
                    gPad->SetLogy(1);
                    frame->GetYaxis()->SetRangeUser(0.5, std::max(1.0, local_max * 1.15));
                } else {
                    gPad->SetLogy(0);
                }

                TGraphErrors* gp = new TGraphErrors(N_PHI_BINS, X.data(), Yp.data(), ex.data(), EYp.data());
                TGraphErrors* gm = new TGraphErrors(N_PHI_BINS, X.data(), Ym.data(), ex.data(), EYm.data());
                if (!gp || !gm) fatalf("[plot][%s] TGraphErrors allocation failed", group_label.c_str());

                gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);  gp->SetLineColor(kRed);  gp->SetLineWidth(1);
                gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue); gm->SetLineWidth(1);
                gp->Draw("PE1 SAME");
                gm->Draw("PE1 SAME");

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.045); lab.SetTextAlign(13); // left-top
                lab.DrawLatex(0.12, 0.88,
                    Form("Q2 (%.3g, %.3g)   t (%.3g, %.3g)",
                         Q2s[r].first, Q2s[r].second, Ts[ccol].first, Ts[ccol].second));
            }
        }

        const std::string fpath = out_dir + "/plot_total_counts_" + group_label + "_xB_" + std::to_string(ix) + ".png";
        c->SaveAs(fpath.c_str());

        std::error_code fec;
        bool exists = fs::exists(fpath, fec);
        auto sz = exists ? fs::file_size(fpath, fec) : 0ULL;
        if (!exists || (sz == 0) || fec) {
            fatalf("[plot][%s] Failed to save plot: %s (exists=%d, size=%llu, ec=%s)",
                   group_label.c_str(), fpath.c_str(), (int)exists,
                   (unsigned long long)sz, (fec ? fec.message().c_str() : "ok"));
        }
        logf("plot", "[%s] Saved %s", group_label.c_str(), fpath.c_str());

        c->Close();
        delete c; c = nullptr;
    }

    // proactive cleanup between groups
    if (gROOT) {
        if (auto* lst = gROOT->GetListOfCanvases()) lst->Delete();
    }
    logf("plot", "[%s] Plotting complete", group_label.c_str());
}

// ===============================================
// Main compute function (STRICT)
// ===============================================
} // namespace

void compute_total_counts(
    const std::vector<std::string>& periods,            // MUST be canonical tree keys (e.g. "DVCS_Fa18_inb")
    const std::vector<std::string>& /*topologies*/,     // not used here
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dataTrees,     // keys are DVCS_* names
    const std::string& combined_cuts_json,              // applied here (base + 3-sigma)
    const std::string& out_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // ROOT in batch, cleaner logs
    if (gROOT) gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    // ---- Validate binning axes ----
    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');
    if (N_PHI_BINS <= 0) fatalf("N_PHI_BINS must be > 0");
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty())
        fatalf("Empty binning axes (xB/Q2/t). Check binning_scheme.");

    logf("init", "Using %d phi bins; xB=%zu, Q2=%zu, t=%zu",
         N_PHI_BINS, xB_bins.size(), Q2_bins.size(), t_bins.size());

    // ---- Prepare output dirs ----
    std::error_code ec;
    const fs::path plots_root = fs::path(out_root_dir) / "total_counts_plots";
    const fs::path jsons_dir  = fs::path(out_root_dir) / "jsons";
    fs::create_directories(plots_root, ec);
    if (ec) fatalf("Cannot create plots root: %s (%s)", plots_root.string().c_str(), ec.message().c_str());
    ec.clear();
    fs::create_directories(jsons_dir, ec);
    if (ec) fatalf("Cannot create jsons dir: %s (%s)", jsons_dir.string().c_str(), ec.message().c_str());

    logf("init", "Output dirs OK: plots_root=%s, jsons_dir=%s",
         plots_root.string().c_str(), jsons_dir.string().c_str());

    // ---- Per-period processing ----
    std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>> allGroupsByLabel;

    for (const auto& period_key : periods) {
        logf("period", "=== Begin period key: %s ===", period_key.c_str());

        if (!is_canonical_tree_key(period_key)) {
            std::ostringstream oss;
            oss << "Periods must be canonical tree keys. Got '" << period_key << "'. Expected one of:";
            for (const auto& p : CANONICAL_PERIODS()) oss << " " << p.tree_key;
            fatalf("%s", oss.str().c_str());
        }
        const char* label = tree_key_to_label_or_fatal(period_key);
        const std::string label_str(label);
        logf("period", "[%s] Canonical label: %s", period_key.c_str(), label);

        auto it = dataTrees.find(period_key);
        if (it == dataTrees.end() || !it->second) {
            std::ostringstream avail;
            avail << "{ ";
            bool first = true;
            for (const auto& kv : dataTrees) {
                if (!first) avail << ", ";
                first = false;
                avail << '"' << kv.first << '"';
            }
            avail << " }";
            fatalf("Missing DVCS data tree for key '%s'. Available keys: %s",
                   period_key.c_str(), avail.str().c_str());
        }

        // ---- Load per-period cuts ----
        const PeriodCuts cuts = loadCombinedCuts(combined_cuts_json, period_key);

        TTree* t = it->second;
        if (!t) fatalf("[%s] Null TTree pointer", period_key.c_str());

        // Baseline branches
        int helicity = 0;
        double x = 0.0, Q2 = 0.0, t1 = 0.0, phi2 = std::numeric_limits<double>::quiet_NaN();

        if (!t->GetBranch("helicity") || !t->GetBranch("x") || !t->GetBranch("Q2") || !t->GetBranch("t1"))
            fatalf("[%s] Required branches (helicity,x,Q2,t1) missing", period_key.c_str());
        if (!t->GetBranch("phi2"))
            fatalf("[%s] Required branch 'phi2' missing", period_key.c_str());

        t->SetBranchAddress("helicity", &helicity);
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("Q2", &Q2);
        t->SetBranchAddress("t1", &t1);
        t->SetBranchAddress("phi2", &phi2);

        logf("bind", "[%s] Bound baseline branches helicity,x,Q2,t1,phi2", period_key.c_str());

        // Bind branches required by cuts (skip reserved)
        const auto neededCutBranches = requiredCutBranches(cuts);
        if (!neededCutBranches.empty()) {
            std::ostringstream oss; oss << "[";
            for (size_t i=0;i<neededCutBranches.size();++i){
                if (i) oss<<", ";
                oss << neededCutBranches[i];
            }
            oss << "]";
            logf("bind", "[%s] Cut branches requested: %s", period_key.c_str(), oss.str().c_str());
        } else {
            logf("bind", "[%s] No extra cut branches requested", period_key.c_str());
        }

        std::unordered_map<std::string, BranchBinding> cutBindings;
        bindRequiredBranches_STRICT(period_key, t, neededCutBranches, cutBindings);

        // Count table
        std::map<std::tuple<int,int,int,int>, HelCounts> table;

        const Long64_t nent = t->GetEntries();
        logf("loop", "[%s] Entries: %lld", period_key.c_str(), (long long)nent);

        int badPhi = 0, badBin = 0;
        for (Long64_t i = 0; i < nent; ++i) {
            if ((i % 2500000LL) == 0) logf("loop", "[%s]   at entry %lld", period_key.c_str(), (long long)i);

            t->GetEntry(i);

            if (helicity != +1 && helicity != -1) continue;
            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi2)) continue;

            const int ix  = findBin(x, xB_bins);
            const int iQ  = findBin(Q2, Q2_bins);
            const int itb = findBin(std::fabs(t1), t_bins);
            int ip  = phiToBin(phi2);

            if (ip < 0) { if (badPhi < 5) logf("warn","[%s] bad phi2=%g at i=%lld", period_key.c_str(), phi2, (long long)i); ++badPhi; continue; }
            if (ix < 0 || iQ < 0 || itb < 0) { if (badBin < 5) logf("warn","[%s] out-of-range bin x=%.4g Q2=%.4g |t1|=%.4g at i=%lld", period_key.c_str(), x,Q2,std::fabs(t1),(long long)i); ++badBin; continue; }

            if (!passes3SigmaCuts_STRICT(cuts, cutBindings)) continue;

            auto& hc = table[{ix, iQ, itb, ip}];
            if (helicity == +1) hc.plus++; else hc.minus++;
        } //endfor

        logf("loop", "[%s] Done loop. badPhi=%d badBin=%d bins_filled=%zu",
             period_key.c_str(), badPhi, badBin, table.size());

        // Save plots for this group
        const fs::path plot_dir = fs::path(out_root_dir) / "total_counts_plots" / label_str;
        plot_group_counts(label_str, table, binning_scheme, xB_bins, Q2_bins, t_bins, plot_dir.string());
        logf("period", "[%s] Plots saved", period_key.c_str());

        allGroupsByLabel[label_str] = std::move(table);

        // Extra safety between periods
        if (gROOT) {
            if (auto* lst = gROOT->GetListOfCanvases()) lst->Delete();
        }
        logf("period", "=== End period key: %s ===", period_key.c_str());
    } //endfor periods

    // ---- Write master and per-group JSON ----
    write_total_counts_master_json(out_json_path, allGroupsByLabel, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    for (const auto& gkv : allGroupsByLabel) {
        const std::string fname = (std::filesystem::path(out_root_dir) / "jsons" / ("total_counts_" + gkv.first + ".json")).string();
        write_total_counts_group_json(fname, gkv.second, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
    } //endfor

    logf("done", "All groups done. Master JSON: %s", out_json_path.c_str());
}