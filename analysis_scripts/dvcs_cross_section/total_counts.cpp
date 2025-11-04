// total_counts.cpp
// ------------------------------------------------------------
// Compute helicity-separated total counts after exclusivity cuts.
// STRICT typed SetBranchAddress bindings (no UB) + hardened plotting/teardown.
// ------------------------------------------------------------

#include "total_counts.h"
#include "periods.h"

#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
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
#include <chrono>
#include <thread>

#include <nlohmann/json.hpp>
using nlohmann::json;

namespace {

// ------------------------------------------------------------
// Constants, helpers
// ------------------------------------------------------------
static constexpr int N_PHI_BINS = 12;
static constexpr double TWO_PI  = 2.0 * M_PI;

[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[total_counts][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}
static inline void fence(const std::string& tag) {
    std::cout << "[fence] " << tag << std::endl;
}

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

// Canonical period helpers
static inline bool is_canonical_tree_key(const std::string& k) {
    for (const auto& p : CANONICAL_PERIODS()) if (k == p.tree_key) return true;
    return false;
}
static inline const char* tree_key_to_label_or_fatal(const std::string& k) {
    for (const auto& p : CANONICAL_PERIODS()) if (k == p.tree_key) return p.label;
    std::ostringstream oss;
    oss << "Non-canonical tree key '" << k << "'. Expected one of:";
    for (const auto& p : CANONICAL_PERIODS()) oss << " " << p.tree_key;
    fatal(oss.str());
    return "";
}

// ------------------------------------------------------------
// Cuts model
// ------------------------------------------------------------
struct SigmaCut {
    std::string branch;
    double center = 0.0;
    double sigma  = 0.0;
    double nsigma = 3.0;
    enum Mode { TWO_SIDED, UPPER, LOWER } mode = TWO_SIDED;
};
struct BaseCut {
    std::string branch;
    bool has_min=false, has_max=false, has_eq=false, has_neq=false;
    double vmin=0.0, vmax=0.0;
    long long eq=0, neq=0;
};
struct PeriodCuts {
    std::vector<SigmaCut> sigma_cuts;
    std::vector<BaseCut>  base_cuts;
};

static SigmaCut::Mode parseMode(const std::string& m) {
    if (m == "two_sided") return SigmaCut::TWO_SIDED;
    if (m == "upper")     return SigmaCut::UPPER;
    if (m == "lower")     return SigmaCut::LOWER;
    fatal("Unknown sigma cut mode: " + m);
    return SigmaCut::TWO_SIDED;
}
static void parseSigmaCuts(const json& arr, std::vector<SigmaCut>& out) {
    if (!arr.is_array()) fatal("sigma_cuts must be an array");
    for (const auto& j : arr) {
        SigmaCut c;
        if (!j.contains("branch")) fatal("sigma_cuts entry missing 'branch'");
        c.branch = j.at("branch").get<std::string>();
        if (!j.contains("center") || !j.contains("sigma"))
            fatal("sigma_cuts '" + c.branch + "' missing center or sigma");
        c.center  = j.at("center").get<double>();
        c.sigma   = j.at("sigma").get<double>();
        c.nsigma  = j.value("nsigma", 3.0);
        c.mode    = parseMode(j.value("mode", std::string("two_sided")));
        out.push_back(c);
    }
}
static void parseBaseCuts(const json& arr, std::vector<BaseCut>& out) {
    if (!arr.is_array()) fatal("base_cuts must be an array");
    for (const auto& j : arr) {
        BaseCut c;
        if (!j.contains("branch")) fatal("base_cuts entry missing 'branch'");
        c.branch = j.at("branch").get<std::string>();
        if (j.contains("min")) { c.has_min = true; c.vmin = j.at("min").get<double>(); }
        if (j.contains("max")) { c.has_max = true; c.vmax = j.at("max").get<double>(); }
        if (j.contains("eq"))  { c.has_eq  = true; c.eq   = j.at("eq").get<long long>(); }
        if (j.contains("neq")) { c.has_neq = true; c.neq  = j.at("neq").get<long long>(); }
        if (!(c.has_min || c.has_max || c.has_eq || c.has_neq))
            fatal("base_cuts '" + c.branch + "' has no condition");
        out.push_back(c);
    }
}
static PeriodCuts loadCombinedCuts(const std::string& json_path, const std::string& period_key) {
    std::ifstream ifs(json_path);
    if (!ifs) fatal("Cannot open combined cuts JSON: " + json_path);
    json J; ifs >> J;

    PeriodCuts cuts;
    if (J.contains("global")) {
        const auto& G = J.at("global");
        if (G.contains("sigma_cuts")) parseSigmaCuts(G.at("sigma_cuts"), cuts.sigma_cuts);
        if (G.contains("base_cuts"))  parseBaseCuts (G.at("base_cuts"),  cuts.base_cuts);
    }
    if (J.contains("period_overrides")) {
        const auto& PO = J.at("period_overrides");
        auto it = PO.find(period_key);
        if (it != PO.end()) {
            const auto& P = *it;
            if (P.contains("sigma_cuts")) { cuts.sigma_cuts.clear(); parseSigmaCuts(P.at("sigma_cuts"), cuts.sigma_cuts); }
            if (P.contains("base_cuts"))  { cuts.base_cuts.clear();  parseBaseCuts (P.at("base_cuts"),  cuts.base_cuts);  }
        }
    }
    for (const auto& c : cuts.sigma_cuts) {
        if (c.sigma  <= 0.0) fatal("sigma_cuts: non-positive sigma for '" + c.branch + "'");
        if (c.nsigma <= 0.0) fatal("sigma_cuts: non-positive nsigma for '" + c.branch + "'");
    }
    return cuts;
}
static std::vector<std::string> requiredCutBranches(const PeriodCuts& cuts) {
    std::set<std::string> s;
    for (const auto& c : cuts.sigma_cuts) s.insert(c.branch);
    for (const auto& c : cuts.base_cuts)  s.insert(c.branch);
    return {s.begin(), s.end()};
}

// ------------------------------------------------------------
// Typed branch binding (no UB).
// ------------------------------------------------------------
enum class LeafKind {
    kDouble, kFloat, kInt64, kUInt64, kInt, kUInt, kShort, kUShort, kChar, kUChar, kBool, kOther
};
static LeafKind kindFromTypeName(const char* tn) {
    std::string t = tn ? std::string(tn) : std::string();
    if (t == "Double_t")  return LeafKind::kDouble;
    if (t == "Float_t")   return LeafKind::kFloat;
    if (t == "Long64_t")  return LeafKind::kInt64;
    if (t == "ULong64_t") return LeafKind::kUInt64;
    if (t == "Int_t")     return LeafKind::kInt;
    if (t == "UInt_t")    return LeafKind::kUInt;
    if (t == "Short_t")   return LeafKind::kShort;
    if (t == "UShort_t")  return LeafKind::kUShort;
    if (t == "Char_t")    return LeafKind::kChar;
    if (t == "UChar_t")   return LeafKind::kUChar;
    if (t == "Bool_t")    return LeafKind::kBool;
    return LeafKind::kOther;
}

struct BranchBinding {
    std::string name;
    LeafKind kind = LeafKind::kOther;
    Double_t  vD = 0.0;
    Float_t   vF = 0.0f;
    Long64_t  vI64 = 0;
    ULong64_t vU64 = 0;
    Int_t     vI = 0;
    UInt_t    vU = 0;
    Short_t   vS = 0;
    UShort_t  vUS = 0;
    Char_t    vC = 0;
    UChar_t   vUC = 0;
    Bool_t    vB = 0;

    inline double asDouble() const {
        switch (kind) {
            case LeafKind::kDouble:  return (double)vD;
            case LeafKind::kFloat:   return (double)vF;
            case LeafKind::kInt64:   return (double)vI64;
            case LeafKind::kUInt64:  return (double)vU64;
            case LeafKind::kInt:     return (double)vI;
            case LeafKind::kUInt:    return (double)vU;
            case LeafKind::kShort:   return (double)vS;
            case LeafKind::kUShort:  return (double)vUS;
            case LeafKind::kChar:    return (double)vC;
            case LeafKind::kUChar:   return (double)vUC;
            case LeafKind::kBool:    return vB ? 1.0 : 0.0;
            default:                 return std::numeric_limits<double>::quiet_NaN();
        }
    }
    inline long long asInt64() const {
        switch (kind) {
            case LeafKind::kInt64:   return (long long)vI64;
            case LeafKind::kUInt64:  return (long long)vU64;
            case LeafKind::kInt:     return (long long)vI;
            case LeafKind::kUInt:    return (long long)vU;
            case LeafKind::kShort:   return (long long)vS;
            case LeafKind::kUShort:  return (long long)vUS;
            case LeafKind::kChar:    return (long long)vC;
            case LeafKind::kUChar:   return (long long)vUC;
            case LeafKind::kBool:    return vB ? 1LL : 0LL;
            case LeafKind::kDouble:  return (long long)std::llround(vD);
            case LeafKind::kFloat:   return (long long)std::llround(vF);
            default:                 return 0LL;
        }
    }
};

static void bindOne_STRICT(TTree* t, BranchBinding& bb) {
    TBranch* b = t->GetBranch(bb.name.c_str());
    if (!b) fatal("Required branch missing: '" + bb.name + "'");
    TLeaf* leaf = b->GetLeaf(bb.name.c_str());
    if (!leaf) {
        leaf = (TLeaf*)b->GetListOfLeaves()->First();
        if (!leaf) fatal("Branch has no leaves: '" + bb.name + "'");
    }
    bb.kind = kindFromTypeName(leaf->GetTypeName());

    switch (bb.kind) {
        case LeafKind::kDouble:  t->SetBranchAddress(bb.name.c_str(), &bb.vD);  break;
        case LeafKind::kFloat:   t->SetBranchAddress(bb.name.c_str(), &bb.vF);  break;
        case LeafKind::kInt64:   t->SetBranchAddress(bb.name.c_str(), &bb.vI64);break;
        case LeafKind::kUInt64:  t->SetBranchAddress(bb.name.c_str(), &bb.vU64);break;
        case LeafKind::kInt:     t->SetBranchAddress(bb.name.c_str(), &bb.vI);  break;
        case LeafKind::kUInt:    t->SetBranchAddress(bb.name.c_str(), &bb.vU);  break;
        case LeafKind::kShort:   t->SetBranchAddress(bb.name.c_str(), &bb.vS);  break;
        case LeafKind::kUShort:  t->SetBranchAddress(bb.name.c_str(), &bb.vUS); break;
        case LeafKind::kChar:    t->SetBranchAddress(bb.name.c_str(), &bb.vC);  break;
        case LeafKind::kUChar:   t->SetBranchAddress(bb.name.c_str(), &bb.vUC); break;
        case LeafKind::kBool:    t->SetBranchAddress(bb.name.c_str(), &bb.vB);  break;
        default:
            fatal("Unsupported leaf type for '" + bb.name + "' (type=" + leaf->GetTypeName() + ")");
    }
}
static void bindMany_STRICT(TTree* t,
                            const std::vector<std::string>& names,
                            std::unordered_map<std::string, BranchBinding>& out) {
    out.clear();
    out.reserve(names.size());
    for (const auto& n : names) {
        auto [it, inserted] = out.emplace(n, BranchBinding{});
        BranchBinding& bb = it->second;
        bb.name = n;
        bindOne_STRICT(t, bb);
    }
}

// ------------------------------------------------------------
// Cut evaluation
// ------------------------------------------------------------
static inline bool passBaseCuts(const std::vector<BaseCut>& baseCuts,
                                const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : baseCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false;
        const BranchBinding& bb = it->second;
        if (c.has_eq)  { if (bb.asInt64() == c.eq ? false : true) {} }
        if (c.has_eq)  { if (bb.asInt64() != c.eq)  return false; }
        if (c.has_neq) { if (bb.asInt64() == c.neq) return false; }
        if (c.has_min) { if (bb.asDouble() <  c.vmin) return false; }
        if (c.has_max) { if (bb.asDouble() >  c.vmax) return false; }
    }
    return true;
}
static inline bool passSigmaCuts(const std::vector<SigmaCut>& sigmaCuts,
                                 const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : sigmaCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false;
        const BranchBinding& bb = it->second;
        const double v  = bb.asDouble();
        const double lo = c.center - c.nsigma * c.sigma;
        const double hi = c.center + c.nsigma * c.sigma;
        if (c.mode == SigmaCut::TWO_SIDED) { if (v < lo || v > hi) return false; }
        else if (c.mode == SigmaCut::UPPER) { if (v > hi) return false; }
        else /*LOWER*/                      { if (v < lo) return false; }
    }
    return true;
}
static inline bool passes3SigmaCuts_STRICT(
    const PeriodCuts& cuts,
    const std::unordered_map<std::string, BranchBinding>& B) {
    return passBaseCuts(cuts.base_cuts, B) && passSigmaCuts(cuts.sigma_cuts, B);
}

// ------------------------------------------------------------
// JSON writers
// ------------------------------------------------------------
static void write_total_counts_group_json(
    const std::string& out_path,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) fatal(std::string("Cannot open for write: ") + out_path);
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
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote " << out_path << "\n";
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
    if (!ofs) fatal(std::string("Cannot open for write: ") + out_path);
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
        }
        ofs << "\n    }}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote master " << out_path << "\n";
}

// ------------------------------------------------------------
// Plotting (with hardened teardown)
// ------------------------------------------------------------
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

    static const std::vector<double> PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
    for (int i = 0; i < N_PHI_BINS; ++i) X[i] = PHI[i];

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) fatal("Cannot create directory: " + out_dir + " (" + ec.message() + ")");

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        std::set<std::pair<double,double>> q2set, tset;
        for (const auto& b : binning_scheme) {
            if (std::make_pair(b.xBmin, b.xBmax) == xB_bins[ix]) {
                q2set.emplace(b.Q2min, b.Q2max);
                tset.emplace(b.tmin,  b.tmax);
            }
        }
        std::vector<std::pair<double,double>> Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double>> Ts (tset.begin(),  tset.end());
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();
        if (nrows <= 0 || ncols <= 0) continue;

        const int W = 260 * ncols + 120;
        const int H = 220 * nrows + 140;

        TCanvas* c = new TCanvas(Form("c_counts_%s_xB%d", group_label.c_str(), ix), "", W, H);
        if (!c) fatal("TCanvas allocation failed");

        TPad* pTop  = new TPad("pTop",  "pTop",  0.0, 0.94, 1.0, 1.0);
        TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.94);
        pTop->SetFillStyle(0);  pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(0.055); head.SetTextAlign(22);
        head.DrawLatex(0.5, 0.55,
            Form("%s   x_B (%.3g, %.3g)", group_label.c_str(), xB_bins[ix].first, xB_bins[ix].second));

        const int cells = nrows * ncols;

        for (int r = 0; r < nrows; ++r) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int cell = r * ncols + ccol + 1;
                if (cell < 1 || cell > cells) continue;

                pGrid->cd(cell);
                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.06);

                int iQ = -1, itb = -1;
                for (int iq = 0; iq < (int)Q2_bins.size(); ++iq) if (pairAlmostEqual(Q2_bins[iq], Q2s[r])) { iQ = iq; break; }
                for (int it = 0; it < (int)t_bins.size();  ++it) if (pairAlmostEqual(t_bins[it],  Ts[ccol])) { itb = it; break; }

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
                gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);  gp->SetLineColor(kRed);  gp->SetLineWidth(1);
                gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue); gm->SetLineWidth(1);
                gp->Draw("PE1 SAME"); gm->Draw("PE1 SAME");

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.045); lab.SetTextAlign(13);
                lab.DrawLatex(0.12, 0.88,
                    Form("Q2 (%.3g, %.3g)   t (%.3g, %.3g)",
                         Q2s[r].first, Q2s[r].second, Ts[ccol].first, Ts[ccol].second));
            }
        }

        const std::string fpath = out_dir + "/plot_total_counts_" + group_label + "_xB_" + std::to_string(ix) + ".png";
        c->Modified(); c->Update();
        c->SaveAs(fpath.c_str());
        // give ROOT a moment to flush I/O and global lists
        gSystem->ProcessEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));

        // sanity check file exists
        std::error_code fec;
        bool exists = fs::exists(fpath, fec);
        auto sz = exists ? fs::file_size(fpath, fec) : 0ULL;
        if (!exists || (sz == 0) || fec) {
            std::ostringstream em;
            em << "[total_counts][FATAL] Failed to save plot: " << fpath
               << " (exists=" << exists
               << ", size=" << sz
               << ", ec=" << (fec ? fec.message() : "ok") << ")";
            std::cerr << em.str() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // close politely
        c->Close();
        delete c;
        gSystem->ProcessEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    }
}

// ------------------------------------------------------------
// Main compute
// ------------------------------------------------------------
} // namespace

void compute_total_counts(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& /*topologies*/,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // Force batch mode: prevents post-SaveAs GUI crashes on headless nodes
    if (gROOT) gROOT->SetBatch(kTRUE);

    // Axes
    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');
    if (N_PHI_BINS <= 0) fatal("N_PHI_BINS must be > 0");
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty())
        fatal("Empty binning axes (xB/Q2/t). Check binning_scheme.");

    // Output dirs
    std::error_code ec;
    const fs::path plots_root = fs::path(out_root_dir) / "total_counts_plots";
    const fs::path jsons_dir  = fs::path(out_root_dir) / "jsons";
    if (!fs::create_directories(plots_root, ec) && ec)
        fatal(std::string("Cannot create plots root: ") + plots_root.string() + " (" + ec.message() + ")");
    ec.clear();
    if (!fs::create_directories(jsons_dir, ec) && ec)
        fatal(std::string("Cannot create jsons dir: ") + jsons_dir.string() + " (" + ec.message() + ")");

    std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>> allGroupsByLabel;

    // Per period
    for (const auto& period_key : periods) {
        if (!is_canonical_tree_key(period_key)) {
            std::ostringstream oss;
            oss << "Periods must be canonical tree keys. Got '" << period_key << "'. Expected one of:";
            for (const auto& p : CANONICAL_PERIODS()) oss << " " << p.tree_key;
            fatal(oss.str());
        }
        const char* label = tree_key_to_label_or_fatal(period_key);
        const std::string label_str(label);

        auto it = dataTrees.find(period_key);
        if (it == dataTrees.end() || !it->second) {
            std::ostringstream avail;
            avail << "{ ";
            bool first = true;
            for (const auto& kv : dataTrees) { if (!first) avail << ", "; first=false; avail << '"' << kv.first << '"'; }
            avail << " }";
            fatal(std::string("Missing DVCS data tree for key '") + period_key + "'. Available keys: " + avail.str());
        }
        TTree* t = it->second;

        fence("period.begin " + label_str);

        // Load cuts and gather required cut branches
        const PeriodCuts cuts = loadCombinedCuts(combined_cuts_json, period_key);
        const auto neededCutBranches = requiredCutBranches(cuts);

        // Bind baseline branches with *typed* bindings
        std::unordered_map<std::string, BranchBinding> base;
        bindMany_STRICT(t, {"helicity","x","Q2","t1","phi2"}, base);

        // Bind cut branches (if any)
        std::unordered_map<std::string, BranchBinding> cutBindings;
        if (!neededCutBranches.empty()) bindMany_STRICT(t, neededCutBranches, cutBindings);

        // Count table
        std::map<std::tuple<int,int,int,int>, HelCounts> table;

        const Long64_t nent = t->GetEntries();
        std::cout << "[loop] [" << period_key << "] Entries: " << nent << std::endl;

        for (Long64_t i = 0; i < nent; ++i) {
            if ((i & ((1<<18)-1)) == 0) std::cout << "[loop] [" << period_key << "] at entry " << i << std::endl;

            t->GetEntry(i);

            const long long helicity = base.at("helicity").asInt64();
            if (helicity != +1 && helicity != -1) continue;

            const double x   = base.at("x").asDouble();
            const double Q2  = base.at("Q2").asDouble();
            const double t1  = base.at("t1").asDouble();
            const double phi = base.at("phi2").asDouble();

            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi)) continue;

            const int ix  = findBin(x, xB_bins);
            const int iQ  = findBin(Q2, Q2_bins);
            const int itb = findBin(std::fabs(t1), t_bins);
            const int ip  = phiToBin(phi);
            if (ix < 0 || iQ < 0 || itb < 0 || ip < 0) continue;

            if (!neededCutBranches.empty()) {
                if (!passes3SigmaCuts_STRICT(cuts, cutBindings)) continue;
            }

            auto& hc = table[{ix, iQ, itb, ip}];
            if (helicity == +1) hc.plus++; else hc.minus++;
        }

        fence("period.plot.begin " + label_str);
        // Save plots per period
        {
            const fs::path plot_dir = fs::path(out_root_dir) / "total_counts_plots" / label_str;
            plot_group_counts(label_str, table, binning_scheme, xB_bins, Q2_bins, t_bins, plot_dir.string());
        }
        fence("period.plot.end " + label_str);

        // Explicit emplace (avoid operator[] rebalancing after heavy ROOT activity)
        fence("map.emplace.begin " + label_str);
        allGroupsByLabel.emplace(label_str, std::move(table));
        fence("map.emplace.end " + label_str);

        // small breather for ROOT global bookkeeping
        gSystem->ProcessEvents();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));

        fence("period.end " + label_str);
    }

    fence("master.json.begin");
    write_total_counts_master_json(out_json_path, allGroupsByLabel,
                                   N_PHI_BINS, xB_bins, Q2_bins, t_bins);
    fence("master.json.end");

    fence("groups.json.begin");
    for (const auto& gkv : allGroupsByLabel) {
        const std::string fname = (std::filesystem::path(out_root_dir) / "jsons" / ("total_counts_" + gkv.first + ".json")).string();
        write_total_counts_group_json(fname, gkv.second, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
    }
    fence("groups.json.end");
}