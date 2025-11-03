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
//
// Per-group JSON:
// {
//   "binning_meta": { "phi_bins": N, "xB_bins": nx, "Q2_bins": nQ, "t_bins": nt },
//   "bins": {
//      "(ix,iQ2,it,ip)": { "helicity": { "+1": Np, "-1": Nm }, "total": Np+Nm }, ...
//   }
// }
//
// Master JSON:
// {
//   "binning_meta": {...},
//   "groups": { "<label>": { "bins": {...} }, ... }
// }
// ------------------------------------------------------------

#include "total_counts.h"
#include "periods.h"  // canonical PeriodDef {label, tree_key}

#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TError.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <algorithm>
#include <cctype>
#include <cmath>
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

// ------------------------------------------------------------
// Helpers and constants
// ------------------------------------------------------------
static constexpr int N_PHI_BINS = 12;
static constexpr double TWO_PI  = 2.0 * M_PI;

struct BinningRange { double lo=0.0; double hi=0.0; };

[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[total_counts][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
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
    fatal(oss.str());
    return ""; // unreachable
}

// ------------------------------------------------------------
// Exclusivity cut model
// ------------------------------------------------------------
struct SigmaCut {
    std::string branch;     // branch to read
    double center = 0.0;    // mean value
    double sigma  = 0.0;    // sigma
    double nsigma = 3.0;    // how many sigmas
    enum Mode { TWO_SIDED, UPPER, LOWER } mode = TWO_SIDED;
};

struct BaseCut {
    std::string branch;    // branch to read
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
    fatal("Unknown sigma cut mode: " + m);
    return SigmaCut::TWO_SIDED; // unreachable
}

static void parseSigmaCuts(const json& arr, std::vector<SigmaCut>& out) {
    if (!arr.is_array()) fatal("sigma_cuts must be an array");
    for (const auto& j : arr) {
        SigmaCut c;
        if (!j.contains("branch")) fatal("sigma_cuts entry missing 'branch'");
        c.branch = j.at("branch").get<std::string>();
        if (!j.contains("center") || !j.contains("sigma"))
            fatal("sigma_cuts entry for '" + c.branch + "' missing 'center' or 'sigma'");
        c.center  = j.at("center").get<double>();
        c.sigma   = j.at("sigma").get<double>();
        c.nsigma  = j.value("nsigma", 3.0);
        c.mode    = parseMode(j.value("mode", std::string("two_sided")));
        out.push_back(c);
    } // #endfor
}

static void parseBaseCuts(const json& arr, std::vector<BaseCut>& out) {
    if (!arr.is_array()) fatal("base_cuts must be an array");
    for (const auto& j : arr) {
        BaseCut c;
        if (!j.contains("branch")) fatal("base_cuts entry missing 'branch'");
        c.branch = j.at("branch").get<std::string>();

        if (j.contains("min")) { c.has_min = true; c.vmin = j.at("min").get<double>(); }
        if (j.contains("max")) { c.has_max = true; c.vmax = j.at("max").get<double>(); }

        if (j.contains("eq"))  { c.has_eq = true;  c.eq  = j.at("eq").get<long long>(); }
        if (j.contains("neq")) { c.has_neq = true; c.neq = j.at("neq").get<long long>(); }

        if (!(c.has_min || c.has_max || c.has_eq || c.has_neq)) {
            fatal("base_cuts entry for '" + c.branch + "' has no condition (min/max/eq/neq)");
        }
        out.push_back(c);
    } // #endfor
}

static PeriodCuts loadCombinedCuts(const std::string& json_path, const std::string& period_key) {
    std::ifstream ifs(json_path);
    if (!ifs) fatal("Cannot open combined cuts JSON: " + json_path);
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
        }
    }

    // sanity
    for (const auto& c : cuts.sigma_cuts) {
        if (c.sigma <= 0.0) fatal("sigma_cuts: non-positive sigma for branch '" + c.branch + "'");
        if (c.nsigma <= 0.0) fatal("sigma_cuts: non-positive nsigma for branch '" + c.branch + "'");
    }

    return cuts;
}

// Collect unique list of branches required by cuts
static std::vector<std::string> requiredCutBranches(const PeriodCuts& cuts) {
    std::set<std::string> s;
    for (const auto& c : cuts.sigma_cuts) s.insert(c.branch);
    for (const auto& c : cuts.base_cuts)  s.insert(c.branch);
    return {s.begin(), s.end()};
}

// SetBranchAddress helpers
struct BranchBinding {
    std::string name;
    double     as_double = std::numeric_limits<double>::quiet_NaN();
    long long  as_ll     = 0;
    bool       is_int    = false; // set at bind time
};

static bool isIntegerLeaf(TLeaf* leaf) {
    if (!leaf) return false;
    const char* t = leaf->GetTypeName();
    // ROOT type names; extend if needed
    return std::string(t) == "Int_t"    || std::string(t) == "UInt_t"   ||
           std::string(t) == "Short_t"  || std::string(t) == "UShort_t" ||
           std::string(t) == "Char_t"   || std::string(t) == "UChar_t"  ||
           std::string(t) == "Long64_t" || std::string(t) == "ULong64_t";
}

static void bindRequiredBranches_STRICT(
    TTree* t,
    const std::vector<std::string>& branch_names,
    std::unordered_map<std::string, BranchBinding>& bindings)
{
    bindings.clear();
    bindings.reserve(branch_names.size()); // avoid rehash so stored addresses stay valid during binding

    for (const auto& bname : branch_names) {
        TBranch* b = t->GetBranch(bname.c_str());
        if (!b) fatal("Required branch for cuts missing: '" + bname + "'");

        TLeaf* leaf = b->GetLeaf(bname.c_str());
        if (!leaf) {
            leaf = (TLeaf*)b->GetListOfLeaves()->First();
            if (!leaf) fatal("Branch has no leaves (unexpected): '" + bname + "'");
        }

        auto [it, inserted] = bindings.emplace(bname, BranchBinding{});
        BranchBinding& bb = it->second;
        bb.name   = bname;
        bb.is_int = isIntegerLeaf(leaf);

        if (bb.is_int) {
            t->SetBranchAddress(bname.c_str(), &bb.as_ll);
        } else {
            t->SetBranchAddress(bname.c_str(), &bb.as_double);
        }
    } // #endfor
}

// Evaluate base cuts
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
    } // #endfor
    return true;
}

// Evaluate sigma cuts
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
    } // #endfor
    return true;
}

// Combined exclusivity predicate: base AND sigma
static inline bool passes3SigmaCuts_STRICT(
    const PeriodCuts& cuts,
    const std::unordered_map<std::string, BranchBinding>& B)
{
    return passBaseCuts(cuts.base_cuts, B) && passSigmaCuts(cuts.sigma_cuts, B);
}

// ------------------------------------------------------------
// JSON writers (strict)
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
    } // #endfor
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
        } // #endfor
        ofs << "\n    }}";
    } // #endfor
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote master " << out_path << "\n";
}

// ------------------------------------------------------------
// Plotting (strict save)
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
    const std::string& group_label, // fa18_inb, etc.
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
    if (ec) {
        std::cerr << "[total_counts][FATAL] Cannot create directory: " << out_dir
                  << " (" << ec.message() << ")\n";
        std::exit(EXIT_FAILURE);
    }

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

        // Unique canvas name to avoid ROOT name clashes
        TCanvas* c = new TCanvas(Form("c_counts_%s_xB%d", group_label.c_str(), ix), "", W, H);
        if (!c) { fatal("TCanvas allocation failed"); }
        c->SetBatch(kTRUE); // ensure no GUI side-effects

        TPad* pTop  = new TPad("pTop",  "pTop",  0.0, 0.94, 1.0, 1.0);
        TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.94);
        pTop->SetFillStyle(0);  pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        if (!gPad) fatal("gPad null after pTop->cd()");
        TLatex head; head.SetNDC(); head.SetTextSize(0.055); head.SetTextAlign(22);
        head.DrawLatex(0.5, 0.55,
            Form("%s   x_B (%.3g, %.3g)", group_label.c_str(), xB_bins[ix].first, xB_bins[ix].second));

        const int cells = nrows * ncols;

        for (int r = 0; r < nrows; ++r) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int cell = r * ncols + ccol + 1;
                if (cell < 1 || cell > cells) continue;

                pGrid->cd(cell);
                if (!gPad) fatal("gPad null after pGrid->cd(cell)");

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
                if (!frame) fatal("DrawFrame returned null");
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

                // Points only, with error bars (no connecting lines)
                TGraphErrors* gp = new TGraphErrors(N_PHI_BINS, X.data(), Yp.data(), ex.data(), EYp.data());
                TGraphErrors* gm = new TGraphErrors(N_PHI_BINS, X.data(), Ym.data(), ex.data(), EYm.data());
                if (!gp || !gm) fatal("TGraphErrors allocation failed");

                gp->SetMarkerStyle(24); // red open
                gp->SetMarkerColor(kRed);
                gp->SetLineColor(kRed);
                gp->SetLineWidth(1);
                gp->Draw("PE1 SAME");     // points only

                gm->SetMarkerStyle(20); // blue filled
                gm->SetMarkerColor(kBlue);
                gm->SetLineColor(kBlue);
                gm->SetLineWidth(1);
                gm->Draw("PE1 SAME");     // points only

                // Annotate the Q2 and t ranges for this subplot
                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.045);
                lab.SetTextAlign(13); // left-top
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
            std::ostringstream em;
            em << "[total_counts][FATAL] Failed to save plot: " << fpath
               << " (exists=" << exists
               << ", size=" << sz
               << ", ec=" << (fec ? fec.message() : "ok") << ")";
            std::cerr << em.str() << std::endl;
            std::exit(EXIT_FAILURE);
        }

        // IMPORTANT: do not delete canvases or clear ROOT lists here.
        // c->Close();           // harmless, but optional
        // delete c;             // DO NOT DELETE: ROOT owns it in gROOT->GetListOfCanvases()
    }
}

} // namespace

// ------------------------------------------------------------
// Main compute function (STRICT, canonical names enforced)
// ------------------------------------------------------------
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

    // Run in batch to avoid any GUI resource ownership surprises.
    if (gROOT) gROOT->SetBatch(kTRUE);

    std::cout << "[init] Using " << N_PHI_BINS << " phi bins; "
              << "xB=" << uniqueRanges(binning_scheme,'x').size()
              << ", Q2=" << uniqueRanges(binning_scheme,'Q').size()
              << ", t="  << uniqueRanges(binning_scheme,'t').size() << "\n";

    // ---- Validate binning axes ----
    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');
    if (N_PHI_BINS <= 0) fatal("N_PHI_BINS must be > 0");
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty())
        fatal("Empty binning axes (xB/Q2/t). Check binning_scheme.");

    // ---- Prepare output dirs ----
    std::error_code ec;
    const fs::path plots_root = fs::path(out_root_dir) / "total_counts_plots";
    const fs::path jsons_dir  = fs::path(out_root_dir) / "jsons";
    if (!fs::create_directories(plots_root, ec) && ec)
        fatal(std::string("Cannot create plots root: ") + plots_root.string() + " (" + ec.message() + ")");
    ec.clear();
    if (!fs::create_directories(jsons_dir, ec) && ec)
        fatal(std::string("Cannot create jsons dir: ") + jsons_dir.string() + " (" + ec.message() + ")");

    // ---- Per-period processing ----
    // Use LABELS for group keys and filenames; use TREE KEYS for lookups.
    std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>> allGroupsByLabel;

    for (const auto& period_key : periods) {
        if (!is_canonical_tree_key(period_key)) {
            std::ostringstream oss;
            oss << "Periods must be canonical tree keys. Got '" << period_key << "'. Expected one of:";
            for (const auto& p : CANONICAL_PERIODS()) oss << " " << p.tree_key;
            fatal(oss.str());
        } // #endif
        const char* label = tree_key_to_label_or_fatal(period_key);

        auto it = dataTrees.find(period_key);
        if (it == dataTrees.end() || !it->second) {
            std::ostringstream avail;
            avail << "{ ";
            bool first = true;
            for (const auto& kv : dataTrees) {
                if (!first) avail << ", ";
                first = false;
                avail << '"' << kv.first << '"';
            } // #endfor
            avail << " }";
            fatal(std::string("Missing DVCS data tree for key '") + period_key + "'. Available keys: " + avail.str());
        } // #endif

        // ---- Load per-period cuts (global with optional period override) ----
        const PeriodCuts cuts = loadCombinedCuts(combined_cuts_json, period_key);

        TTree* t = it->second;
        if (!t) fatal("Null TTree pointer for " + period_key);

        // Strict: require baseline branches (for binning and helicity)
        int helicity = 0;
        double x = 0, Q2 = 0, t1 = 0, phi2 = std::numeric_limits<double>::quiet_NaN();

        if (!t->GetBranch("helicity") || !t->GetBranch("x") || !t->GetBranch("Q2") || !t->GetBranch("t1"))
            fatal(std::string("Required branches (helicity,x,Q2,t1) missing in '") + period_key + "'");
        if (!t->GetBranch("phi2"))
            fatal(std::string("Required branch 'phi2' missing in '") + period_key + "'");

        t->SetBranchAddress("helicity", &helicity);
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("Q2", &Q2);
        t->SetBranchAddress("t1", &t1);
        t->SetBranchAddress("phi2", &phi2);

        // Bind all branches required by exclusivity/base cuts
        const auto neededCutBranches = requiredCutBranches(cuts);

        std::unordered_map<std::string, BranchBinding> cutBindings;
        bindRequiredBranches_STRICT(t, neededCutBranches, cutBindings);

        // Count table
        std::map<std::tuple<int,int,int,int>, HelCounts> table;

        const Long64_t nent = t->GetEntries();
        std::cout << "[loop] [" << period_key << "] Entries: " << nent << "\n";
        for (Long64_t i = 0; i < nent; ++i) {
            if ((i & 0x3FFFF) == 0) {
                std::cout << "[loop] [" << period_key << "] at entry " << i << "\n";
            }
            t->GetEntry(i);

            // Basic sanity for helicity
            if (helicity != +1 && helicity != -1) continue;

            // Compute bin indices (note: |t1| into t-bins)
            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi2)) continue;
            const int ix  = findBin(x, xB_bins);
            const int iQ  = findBin(Q2, Q2_bins);
            const int itb = findBin(std::fabs(t1), t_bins);
            const int ip  = phiToBin(phi2);
            if (ix < 0 || iQ < 0 || itb < 0 || ip < 0) continue;

            // Apply exclusivity cuts (base + 3-sigma), strict
            if (!passes3SigmaCuts_STRICT(cuts, cutBindings)) continue;

            // Count
            auto& hc = table[{ix, iQ, itb, ip}];
            if (helicity == +1) hc.plus++; else hc.minus++;
        } // #endfor

        const std::string label_str(label);

        // Store by LABEL for downstream consistency
        allGroupsByLabel[label_str] = table;

        // Plots go to .../total_counts_plots/<label>/
        const fs::path plot_dir = fs::path(out_root_dir) / "total_counts_plots" / label_str;
        std::cout << "[plot] [" << label_str << "] Creating plots in: " << plot_dir.string() << "\n";
        plot_group_counts(label_str, table, binning_scheme, xB_bins, Q2_bins, t_bins, plot_dir.string());
        std::cout << "[period] [" << period_key << "] Plots saved\n";
    } // #endfor

    // ---- Write master JSON ----
    write_total_counts_master_json(out_json_path, allGroupsByLabel, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // ---- Write one flat JSON per group (label) ----
    for (const auto& gkv : allGroupsByLabel) {
        const std::string fname = (std::filesystem::path(out_root_dir) / "jsons" / ("total_counts_" + gkv.first + ".json")).string();
        write_total_counts_group_json(fname, gkv.second, N_PHI_BINS, xB_bins, Q2_bins, t_bins);
    } // #endfor

    // ---- Safe final cleanup of canvases (let ROOT do the deletion exactly once) ----
    if (gROOT) {
        if (auto* lst = gROOT->GetListOfCanvases()) {
            std::cout << "[cleanup] Deleting " << lst->GetSize() << " canvases at end of run\n";
            lst->Delete(); // SAFE HERE: we never manually deleted canvases
        }
    }
}