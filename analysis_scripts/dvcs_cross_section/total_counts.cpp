#include "total_counts.h"
#include "periods.h"                // CANONICAL_PERIODS(), PeriodDef{label, tree_key}
#include "load_binning_scheme.h"    // shared project type (harmless include here)
#include "csv_io.h"                 // same CSV helper used by bin_means

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TROOT.h>
#include <TError.h>
#include <TLeaf.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// If json header is not already pulled in elsewhere:
#include <nlohmann/json.hpp>

namespace {

[[noreturn]] void fatal(const std::string& msg) {
    std::cerr << "[total_counts] FATAL: " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline double wrap_deg_0_360(double deg) {
    if (!std::isfinite(deg)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(deg, 360.0);
    if (w < 0.0) w += 360.0;
    if (w >= 360.0) w = std::nextafter(360.0, 0.0);
    return w;
}

static inline double rad_to_deg(double rad) { return rad * 180.0 / M_PI; }

// ---------------- cuts ----------------
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
    fatal("Unknown sigma cut mode: " + m);
    return SigmaCut::TWO_SIDED;
}

static void parseSigmaCuts(const nlohmann::json& arr, std::vector<SigmaCut>& out) {
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
    }
}

static void parseBaseCuts(const nlohmann::json& arr, std::vector<BaseCut>& out) {
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
            fatal("base_cuts entry for '" + c.branch + "' has no condition");
        out.push_back(c);
    }
}

static PeriodCuts loadCombinedCuts(const std::string& json_path, const std::string& period_key) {
    std::ifstream ifs(json_path);
    if (!ifs) fatal("Cannot open combined cuts JSON: " + json_path);
    nlohmann::json J; ifs >> J;

    PeriodCuts cuts;
    if (J.contains("global")) {
        const auto& G = J.at("global");
        if (G.contains("sigma_cuts")) parseSigmaCuts(G.at("sigma_cuts"), cuts.sigma_cuts);
        if (G.contains("base_cuts"))  parseBaseCuts(G.at("base_cuts"),  cuts.base_cuts);
    }
    if (J.contains("period_overrides")) {
        const auto& PO = J.at("period_overrides");
        auto it = PO.find(period_key);
        if (it != PO.end()) {
            const auto& P = *it;
            if (P.contains("sigma_cuts")) { cuts.sigma_cuts.clear(); parseSigmaCuts(P.at("sigma_cuts"), cuts.sigma_cuts); }
            if (P.contains("base_cuts"))  { cuts.base_cuts.clear();  parseBaseCuts(P.at("base_cuts"),  cuts.base_cuts);  }
        }
    }
    for (const auto& c : cuts.sigma_cuts) {
        if (c.sigma <= 0.0) fatal("sigma_cuts: non-positive sigma for '" + c.branch + "'");
        if (c.nsigma <= 0.0) fatal("sigma_cuts: non-positive nsigma for '" + c.branch + "'");
    }
    return cuts;
}

// -------------- binding helpers --------------
struct BranchBinding {
    std::string name;
    double     as_double = std::numeric_limits<double>::quiet_NaN();
    long long  as_ll     = 0;
    bool       is_int    = false;
};

static bool isIntegerLeaf(TLeaf* leaf) {
    if (!leaf) return false;
    const std::string t = leaf->GetTypeName();
    return t == "Int_t" || t == "UInt_t" || t == "Short_t" || t == "UShort_t" ||
           t == "Char_t" || t == "UChar_t" || t == "Long64_t" || t == "ULong64_t";
}

static void bindRequiredBranches_STRICT(
    TTree* t,
    const std::vector<std::string>& branch_names,
    std::unordered_map<std::string, BranchBinding>& bindings)
{
    bindings.clear();
    bindings.reserve(branch_names.size());
    for (const auto& bname : branch_names) {
        TBranch* b = t->GetBranch(bname.c_str());
        if (!b) fatal("Required branch for cuts missing: '" + bname + "'");
        TLeaf* leaf = b->GetLeaf(bname.c_str());
        if (!leaf) {
            leaf = (TLeaf*)b->GetListOfLeaves()->First();
            if (!leaf) fatal("Branch has no leaves: '" + bname + "'");
        }
        auto [it, ok] = bindings.emplace(bname, BranchBinding{});
        BranchBinding& bb = it->second;
        bb.name   = bname;
        bb.is_int = isIntegerLeaf(leaf);
        if (bb.is_int) t->SetBranchAddress(bname.c_str(), &bb.as_ll);
        else           t->SetBranchAddress(bname.c_str(), &bb.as_double);
    }
}

static inline bool passBaseCuts(const std::vector<BaseCut>& baseCuts,
                                const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : baseCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false;
        const BranchBinding& bb = it->second;

        if (c.has_eq)  { long long v = bb.is_int ? bb.as_ll : (long long)std::llround(bb.as_double); if (v != c.eq)  return false; }
        if (c.has_neq) { long long v = bb.is_int ? bb.as_ll : (long long)std::llround(bb.as_double); if (v == c.neq) return false; }
        if (c.has_min) { double v    = bb.is_int ? (double)bb.as_ll : bb.as_double; if (!(v >= c.vmin)) return false; }
        if (c.has_max) { double v    = bb.is_int ? (double)bb.as_ll : bb.as_double; if (!(v <= c.vmax)) return false; }
    }
    return true;
}

static inline bool passSigmaCuts(const std::vector<SigmaCut>& sigmaCuts,
                                 const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : sigmaCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false;
        const BranchBinding& bb = it->second;
        const double v  = bb.is_int ? (double)bb.as_ll : bb.as_double;
        const double lo = c.center - c.nsigma * c.sigma;
        const double hi = c.center + c.nsigma * c.sigma;
        if (c.mode == SigmaCut::TWO_SIDED) { if (v < lo || v > hi) return false; }
        else if (c.mode == SigmaCut::UPPER) { if (v > hi) return false; }
        else { if (v < lo) return false; }
    }
    return true;
}

static inline bool passes3SigmaCuts_STRICT(
    const PeriodCuts& cuts,
    const std::unordered_map<std::string, BranchBinding>& B)
{
    return passBaseCuts(cuts.base_cuts, B) && passSigmaCuts(cuts.sigma_cuts, B);
}

// -------------- canonical period helpers --------------
static inline bool is_canonical_tree_key(const std::string& k) {
    for (const auto& p : CANONICAL_PERIODS()) if (k == p.tree_key) return true;
    return false;
}
static inline std::string safe_label_for_key(const std::string& k) {
    for (const auto& p : CANONICAL_PERIODS()) if (k == p.tree_key) return std::string(p.label);
    fatal("Non-canonical period key: " + k);
    return std::string();
}

// -------------- topology: fixed set and detector-based resolver --------------
static const char* TOPO_STRS[3] = {"(FD, FD)","(CD, FD)","(CD, FT)"};
static const char* TOPO_DIRS[3] = {"FD_FD","CD_FD","CD_FT"}; // for filenames

struct TopologyResolver {
    int detector1 = 0;
    int detector2 = 0;
    bool have_det1 = false, have_det2 = false;

    void bind(TTree* t) {
        have_det1 = (t->GetBranch("detector1") != nullptr);
        have_det2 = (t->GetBranch("detector2") != nullptr);
        if (!(have_det1 && have_det2)) {
            fatal("Missing detector1/detector2 in DVCS tree.");
        }
        t->SetBranchAddress("detector1", &detector1);
        t->SetBranchAddress("detector2", &detector2);
    }

    // Return 0:(FD,FD), 1:(CD,FD), 2:(CD,FT), or -1 if unknown
    int index() const {
        if (detector1 == 1 && detector2 == 1) return 0; // FD,FD
        if (detector1 == 2 && detector2 == 1) return 1; // CD,FD
        if (detector1 == 2 && detector2 == 0) return 2; // CD,FT
        return -1;
    }
};

// -------------- CSV helpers --------------
struct CsvCols {
    int c_xb_min = -1, c_xb_max = -1;
    int c_q2_min = -1, c_q2_max = -1;
    int c_tab_min = -1, c_tab_max = -1; // |t|
    int c_phi_min = -1, c_phi_max = -1;
};

static inline std::string period_yield_col_base(const std::string& label,
                                                const std::string& topo_str) {
    std::ostringstream os;
    os << "raw yield, ep->epg, " << topo_str << ", exp, " << label << ", ";
    return os.str();
}

static inline std::string combined_yield_col_base(const std::string& group,
                                                  const std::string& topo_str) {
    std::ostringstream os;
    os << "raw yield, ep->epg, " << topo_str << ", exp, " << group << ", ";
    return os.str();
}

static inline std::string topo_dir_name(int topo_idx) {
    return TOPO_DIRS[topo_idx];
}

struct RowCounts { long long pos=0, neg=0; };

static inline double safe_mean(const std::vector<double>& v) {
    double s=0.0; int n=0;
    for (double x : v) if (std::isfinite(x)) { s+=x; ++n; }
    return n>0 ? s/n : std::numeric_limits<double>::quiet_NaN();
}

// Draw one label+topology set of canvases (one canvas per xB range)
// NOTE: output directory is output/total_counts_plots/<Label> (no topology subfolders)
static void draw_group_canvases(
    const std::string& label,
    const std::string& topo_str,
    const std::string& topo_dir,
    const CsvDoc& csv,
    const CsvCols& cols,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // Write into: output/total_counts_plots/<Label>
    const std::string out_dir = (fs::path(out_root_dir) / "total_counts_plots" / label).string();
    if (!fs::exists(out_dir)) {
        fatal("Missing plot output directory: " + out_dir + " (makeOutputDirs must create it)");
    }

    // Per-label means and phi averages
    const std::string c_phiavg_name = std::string("phiavg, ") + label;
    const std::string c_q2avg_name  = std::string("Q2avg, ") + label;
    const std::string c_tabavg_name = std::string("t_abs_avg, ") + label;
    const std::string c_xbavg_name  = std::string("xBavg, ") + label;

    const int c_phiavg = csv.has_col(c_phiavg_name) ? csv.col_index(c_phiavg_name) : -1;
    const int c_q2avg  = csv.has_col(c_q2avg_name)  ? csv.col_index(c_q2avg_name)  : -1;
    const int c_tabavg = csv.has_col(c_tabavg_name) ? csv.col_index(c_tabavg_name) : -1;
    const int c_xbavg  = csv.has_col(c_xbavg_name)  ? csv.col_index(c_xbavg_name)  : -1;

    // Yield columns for this label+topology
    const std::string ybase = period_yield_col_base(label, topo_str);
    const int c_pos = csv.has_col(ybase + "pos")   ? csv.col_index(ybase + "pos")   : -1;
    const int c_neg = csv.has_col(ybase + "neg")   ? csv.col_index(ybase + "neg")   : -1;

    // Unique xB ranges in rows
    std::set<std::pair<double,double>> xb_set;
    for (int r = 0; r < csv.nrows(); ++r) {
        xb_set.emplace(csv.as_double(r, cols.c_xb_min), csv.as_double(r, cols.c_xb_max));
    }

    for (auto xb : xb_set) {
        // collect Q2 ranges for this xB
        std::set<std::pair<double,double>> q2set;
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, cols.c_xb_min);
            const double xbmax = csv.as_double(r, cols.c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 && std::fabs(xbmax - xb.second) < 1e-9) {
                q2set.emplace(csv.as_double(r, cols.c_q2_min), csv.as_double(r, cols.c_q2_max));
            }
        }

        // find t sets for each Q2, then collect rows per (Q2,t)
        std::set<std::pair<double,double>> tset_all;
        for (auto q2r : q2set) {
            for (int r = 0; r < csv.nrows(); ++r) {
                const double xbmin = csv.as_double(r, cols.c_xb_min);
                const double xbmax = csv.as_double(r, cols.c_xb_max);
                const double q2min = csv.as_double(r, cols.c_q2_min);
                const double q2max = csv.as_double(r, cols.c_q2_max);
                if (std::fabs(xbmin - xb.first) < 1e-9 && std::fabs(xbmax - xb.second) < 1e-9 &&
                    std::fabs(q2min - q2r.first) < 1e-9 && std::fabs(q2max - q2r.second) < 1e-9) {
                    tset_all.emplace(csv.as_double(r, cols.c_tab_min), csv.as_double(r, cols.c_tab_max));
                }
            }
        }
        std::vector<std::pair<double,double>> Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double>> Ts (tset_all.begin(), tset_all.end());
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();
        const int W = 260 * ncols + 120;
        const int H = 220 * nrows + 160;

        // Canvas-level <xB> for title
        std::vector<double> xbmeans;
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, cols.c_xb_min);
            const double xbmax = csv.as_double(r, cols.c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 && std::fabs(xbmax - xb.second) < 1e-9) {
                if (c_xbavg >= 0) xbmeans.push_back(csv.as_double(r, c_xbavg));
                else xbmeans.push_back((xb.first + xb.second) * 0.5);
            }
        }
        const double xb_mean_for_title = safe_mean(xbmeans);

        std::string cname = "c_counts_" + label + "_" + topo_dir + "_xB_" + std::to_string((int)std::round(xb.first*1000.0));
        TCanvas* c = new TCanvas(cname.c_str(), "", W, H);
        TPad* pTop  = new TPad("pTop",  "pTop",  0.0, 0.93, 1.0, 1.0);
        TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.93);
        pTop->SetFillStyle(0);  pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(0.055); head.SetTextAlign(22);
        head.DrawLatex(0.5, 0.55,
            Form("%s   <xB>=%.3g   %s", label.c_str(), xb_mean_for_title, topo_str.c_str()));

        for (int rrow = 0; rrow < nrows; ++rrow) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int cell = rrow * ncols + ccol + 1;
                pGrid->cd(cell);
                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.06);

                const auto& qpair = Q2s[rrow];
                const auto& tpair = Ts[ccol];

                // gather rows for this cell
                std::vector<int> rows;
                rows.reserve(24);
                for (int r = 0; r < csv.nrows(); ++r) {
                    const double xbmin = csv.as_double(r, cols.c_xb_min);
                    const double xbmax = csv.as_double(r, cols.c_xb_max);
                    const double q2min = csv.as_double(r, cols.c_q2_min);
                    const double q2max = csv.as_double(r, cols.c_q2_max);
                    const double tmin  = csv.as_double(r, cols.c_tab_min);
                    const double tmax  = csv.as_double(r, cols.c_tab_max);
                    if (std::fabs(xbmin - xb.first)  < 1e-9 && std::fabs(xbmax - xb.second) < 1e-9 &&
                        std::fabs(q2min - qpair.first) < 1e-9 && std::fabs(q2max - qpair.second) < 1e-9 &&
                        std::fabs(tmin  - tpair.first) < 1e-9 && std::fabs(tmax  - tpair.second) < 1e-9) {
                        rows.push_back(r);
                    }
                }
                if (rows.empty()) {
                    TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 1.0);
                    frame->GetXaxis()->SetTitle("phi (deg)");
                    frame->GetYaxis()->SetTitle("Counts");
                    continue;
                }

                std::sort(rows.begin(), rows.end(), [&](int a, int b){
                    return csv.as_double(a, cols.c_phi_min) < csv.as_double(b, cols.c_phi_min);
                });

                std::vector<double> X, Yp, Ym, EX, EYp, EYm, q2means, tmeans;
                X.reserve(rows.size()); Yp.reserve(rows.size()); Ym.reserve(rows.size());
                EX.resize(rows.size(), 0.0);

                for (int r : rows) {
                    const double pmin = csv.as_double(r, cols.c_phi_min);
                    const double pmax = csv.as_double(r, cols.c_phi_max);
                    double xphi = 0.5 * (pmin + pmax);
                    if (c_phiavg >= 0) {
                        const double pav = csv.as_double(r, c_phiavg);
                        if (std::isfinite(pav) && pav > 0.0 && pav < 360.0) xphi = pav;
                    }
                    X.push_back(xphi);

                    double yp=0.0, yn=0.0;
                    if (c_pos >= 0) yp = csv.as_double(r, c_pos);
                    if (c_neg >= 0) yn = csv.as_double(r, c_neg);
                    Yp.push_back(yp);
                    Ym.push_back(yn);
                    EYp.push_back(std::sqrt(std::max(0.0, yp)));
                    EYm.push_back(std::sqrt(std::max(0.0, yn)));

                    if (c_q2avg >= 0)  q2means.push_back(csv.as_double(r, c_q2avg));
                    else               q2means.push_back(0.5 * (qpair.first + qpair.second));
                    if (c_tabavg >= 0) tmeans.push_back(csv.as_double(r, c_tabavg));
                    else               tmeans.push_back(0.5 * (tpair.first + tpair.second));
                }

                double local_max = 1.0;
                for (double v : Yp) local_max = std::max(local_max, v);
                for (double v : Ym) local_max = std::max(local_max, v);

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, std::max(1.0, local_max * 1.15));
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

                TGraphErrors* gp = new TGraphErrors((int)X.size(), X.data(), Yp.data(), EX.data(), EYp.data());
                TGraphErrors* gm = new TGraphErrors((int)X.size(), X.data(), Ym.data(), EX.data(), EYm.data());
                gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);  gp->SetLineColor(kRed);  gp->SetLineWidth(1);
                gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue); gm->SetLineWidth(1);
                gp->Draw("PE1 SAME");
                gm->Draw("PE1 SAME");

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.045); lab.SetTextAlign(13);
                const double q2m = safe_mean(q2means);
                const double tm  = safe_mean(tmeans);
                lab.DrawLatex(0.12, 0.88, Form("Q2=%.3g  |t|=%.3g", q2m, tm));
            }
        }

        const std::string fpath = out_dir + "/plot_total_counts_" + label + "_" + topo_dir
                                  + "_xB_" + std::to_string((int)std::round(xb.first*1000.0)) + ".png";
        c->SaveAs(fpath.c_str());
        delete c;
    }
}

// Combined group tag
static inline std::string group_tag_for_key(const std::string& period_key) {
    if (period_key.find("Fa18") != std::string::npos) return "Fa18";
    if (period_key.find("Sp18") != std::string::npos) return "Sp18";
    return "10.6 GeV"; // includes Sp19 etc.
}

} // namespace

// ---------------- entry ----------------
bool update_total_counts_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_root_dir,
    int max_workers)
{
    namespace fs = std::filesystem;

    // Assemble canonical DVCS period keys once, standardized for the whole pipeline
    std::vector<std::string> periods;
    periods.reserve(CANONICAL_PERIODS().size());
    for (const auto& P : CANONICAL_PERIODS()) {
        periods.push_back(P.tree_key); // e.g. "DVCS_Fa18_inb"
    }

    // Backup
    try {
        const std::string backup = "output/csvs/dvcs_pass2_analysis_backup_total_counts.csv";
        fs::copy_file(csv_path, backup, fs::copy_options::overwrite_existing);
        std::cout << "[total_counts] Backed up CSV to " << backup << "\n";
    } catch (const std::exception& e) {
        std::cerr << "[total_counts] WARNING: backup failed (" << e.what() << "). Continuing.\n";
    }

    // Load CSV
    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to read " << csv_path << "\n";
        return false;
    }

    // Bin-edge columns
    CsvCols cols;
    cols.c_xb_min  = csv.col_index("xBmin");
    cols.c_xb_max  = csv.col_index("xBmax");
    cols.c_q2_min  = csv.col_index("Q2min");
    cols.c_q2_max  = csv.col_index("Q2max");
    if (csv.has_col("t_abs_min") && csv.has_col("t_abs_max")) {
        cols.c_tab_min = csv.col_index("t_abs_min");
        cols.c_tab_max = csv.col_index("t_abs_max");
    } else {
        cols.c_tab_min = csv.col_index("tmin");
        cols.c_tab_max = csv.col_index("tmax");
    }
    cols.c_phi_min = csv.col_index("phimin");
    cols.c_phi_max = csv.col_index("phimax");
    if (cols.c_xb_min < 0 || cols.c_xb_max < 0 || cols.c_q2_min < 0 || cols.c_q2_max < 0
        || cols.c_tab_min < 0 || cols.c_tab_max < 0 || cols.c_phi_min < 0 || cols.c_phi_max < 0) {
        fatal("Missing one or more required bin-edge columns in CSV.");
    }

    // Ensure per-period/per-topology columns exist as we go
    auto ensure_yield_columns = [&](const std::string& label, const std::string& topo_str) {
        const std::string base = period_yield_col_base(label, topo_str);
        csv.ensure_col(base + "unpol");
        csv.ensure_col(base + "pos");
        csv.ensure_col(base + "neg");
    };

    // Thread-safe counts store: label -> topo_str -> (row -> RowCounts)
    using PeriodRowMap = std::unordered_map<int, RowCounts>;
    std::mutex merge_mtx;
    std::map<std::string, std::map<std::string, PeriodRowMap>> counts_by_label_topo;

    // Warm up ROOT
    ROOT::EnableThreadSafety(); (void)TGraphErrors::Class();

    // Process periods in parallel
    #pragma omp parallel for schedule(dynamic) num_threads(max_workers)
    for (int ip = 0; ip < (int)periods.size(); ++ip) {
        const std::string period_key = periods[ip];
        if (!is_canonical_tree_key(period_key)) {
            #pragma omp critical
            { std::cerr << "[total_counts] Non-canonical period: " << period_key << "\n"; }
            continue;
        }
        const std::string label = safe_label_for_key(period_key);

        auto itT = dvcsDataTrees.find(period_key);
        if (itT == dvcsDataTrees.end() || !itT->second) {
            #pragma omp critical
            { std::cerr << "[total_counts] Missing TTree for " << period_key << "\n"; }
            continue;
        }
        TTree* t = itT->second;
        const PeriodCuts cuts = loadCombinedCuts(combined_cuts_json, period_key);

        // Baseline branches (must exist)
        int helicity = 0; double x=0.0, Q2=0.0, t1=0.0, phi2 = std::numeric_limits<double>::quiet_NaN();
        if (!t->GetBranch("helicity") || !t->GetBranch("x") || !t->GetBranch("Q2") || !t->GetBranch("t1") || !t->GetBranch("phi2")) {
            #pragma omp critical
            { std::cerr << "[total_counts] Required branches missing in " << period_key << "\n"; }
            continue;
        }
        t->SetBranchAddress("helicity", &helicity);
        t->SetBranchAddress("x", &x);
        t->SetBranchAddress("Q2", &Q2);
        t->SetBranchAddress("t1", &t1);
        t->SetBranchAddress("phi2", &phi2);

        // Bind cut branches
        std::vector<std::string> need;
        need.reserve(cuts.base_cuts.size() + cuts.sigma_cuts.size());
        for (const auto& c : cuts.base_cuts)  need.push_back(c.branch);
        for (const auto& c : cuts.sigma_cuts) need.push_back(c.branch);
        std::sort(need.begin(), need.end()); need.erase(std::unique(need.begin(), need.end()), need.end());
        std::unordered_map<std::string, BranchBinding> bind;
        bindRequiredBranches_STRICT(t, need, bind);

        // Topology resolver (detector-based only)
        TopologyResolver topo; topo.bind(t);

        // Thread-local map: topo_str -> (row -> RowCounts)
        std::map<std::string, PeriodRowMap> local_by_topo;

        const Long64_t nent = t->GetEntries();
        for (Long64_t i = 0; i < nent; ++i) {
            t->GetEntry(i);
            if (helicity != +1 && helicity != -1) continue;
            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi2)) continue;
            if (!passes3SigmaCuts_STRICT(cuts, bind)) continue;

            const int topo_idx = topo.index();
            if (topo_idx < 0 || topo_idx > 2) continue;
            const std::string topo_str = TOPO_STRS[topo_idx];

            const double xb  = x;
            const double q2  = Q2;
            const double tab = std::fabs(t1);
            const double phi_deg = wrap_deg_0_360(rad_to_deg(phi2));

            // Find the CSV row whose bin edges contain this event
            int found_row = -1;
            for (int r = 0; r < csv.nrows(); ++r) {
                const double xbmin  = csv.as_double(r, cols.c_xb_min);
                const double xbmax  = csv.as_double(r, cols.c_xb_max);
                const double q2min  = csv.as_double(r, cols.c_q2_min);
                const double q2max  = csv.as_double(r, cols.c_q2_max);
                const double tabmin = csv.as_double(r, cols.c_tab_min);
                const double tabmax = csv.as_double(r, cols.c_tab_max);
                const double pmin   = csv.as_double(r, cols.c_phi_min);
                const double pmax   = csv.as_double(r, cols.c_phi_max);
                if (xb  >= xbmin && xb  < xbmax &&
                    q2  >= q2min && q2  < q2max &&
                    tab >= tabmin && tab < tabmax &&
                    phi_deg >= pmin && phi_deg < pmax) {
                    found_row = r; break;
                }
            }
            if (found_row < 0) continue;

            RowCounts& rc = local_by_topo[topo_str][found_row];
            if (helicity == +1) rc.pos++; else rc.neg++;
        }

        // Merge into global
        {
            std::lock_guard<std::mutex> lock(merge_mtx);
            for (const auto& tk : local_by_topo) {
                const std::string topo_str = tk.first;
                // ensure columns for this label+topology
                ensure_yield_columns(label, topo_str);

                auto& tgt = counts_by_label_topo[label][topo_str];
                for (const auto& rkv : tk.second) {
                    tgt[rkv.first].pos += rkv.second.pos;
                    tgt[rkv.first].neg += rkv.second.neg;
                }
            }
        }
    } // omp

    // Write per-period, per-topology counts to CSV
    for (const auto& lblkv : counts_by_label_topo) {
        const std::string& label = lblkv.first;
        for (int topo_idx = 0; topo_idx < 3; ++topo_idx) {
            const std::string topo_str = TOPO_STRS[topo_idx];
            const auto itTopo = lblkv.second.find(topo_str);
            if (itTopo == lblkv.second.end()) continue;

            const std::string base = period_yield_col_base(label, topo_str);
            const int c_unpol = csv.ensure_col(base + "unpol");
            const int c_pos   = csv.ensure_col(base + "pos");
            const int c_neg   = csv.ensure_col(base + "neg");

            for (const auto& rkv : itTopo->second) {
                const int r = rkv.first;
                const long long pos = rkv.second.pos;
                const long long neg = rkv.second.neg;
                const long long unp = pos + neg;
                csv.set(r, c_unpol, (double)unp);
                csv.set(r, c_pos,   (double)pos);
                csv.set(r, c_neg,   (double)neg);
            }
        }
        std::cout << "[total_counts] wrote per-topology counts for " << label << "\n";
    }

    // Build combined groups per topology (always all three)
    std::map<std::string, std::vector<std::string>> group_members = {
        {"Fa18",    {}},
        {"Sp18",    {}},
        {"10.6 GeV",{}}
    };
    for (const auto& p : periods) {
        const std::string g   = group_tag_for_key(p);
        const std::string lbl = safe_label_for_key(p);
        group_members[g].push_back(lbl);
    }

    for (const auto& gkv : group_members) {
        const std::string group = gkv.first;
        const auto& members = gkv.second;

        for (int topo_idx = 0; topo_idx < 3; ++topo_idx) {
            const std::string topo_str = TOPO_STRS[topo_idx];
            const std::string gbase = combined_yield_col_base(group, topo_str);
            const int c_unpol = csv.ensure_col(gbase + "unpol");
            const int c_pos   = csv.ensure_col(gbase + "pos");
            const int c_neg   = csv.ensure_col(gbase + "neg");

            for (int r = 0; r < csv.nrows(); ++r) {
                long long pos_sum=0, neg_sum=0;
                for (const auto& lbl : members) {
                    const std::string pbase = period_yield_col_base(lbl, topo_str);
                    const int p_pos = csv.col_index(pbase + "pos");
                    const int p_neg = csv.col_index(pbase + "neg");
                    if (p_pos >= 0) pos_sum += (long long)std::llround(csv.as_double(r, p_pos));
                    if (p_neg >= 0) neg_sum += (long long)std::llround(csv.as_double(r, p_neg));
                }
                csv.set(r, c_pos,   (double)pos_sum);
                csv.set(r, c_neg,   (double)neg_sum);
                csv.set(r, c_unpol, (double)(pos_sum + neg_sum));
            }
        }
        std::cout << "[total_counts] wrote combined counts for " << group << " (per topology)\n";
    }

    // Save CSV
    if (!csv.save(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to save " << csv_path << "\n";
        return false;
    }
    std::cout << "[total_counts] Updated raw yields in: " << csv_path << "\n";

    // Plots: per period and per group, per topology (fixed 3)
    for (const auto& p : periods) {
        const std::string label = safe_label_for_key(p);
        for (int topo_idx = 0; topo_idx < 3; ++topo_idx) {
            draw_group_canvases(label, TOPO_STRS[topo_idx], topo_dir_name(topo_idx), csv, cols, out_root_dir);
        }
    }
    for (const auto& gkv : group_members) {
        const std::string group = gkv.first;
        for (int topo_idx = 0; topo_idx < 3; ++topo_idx) {
            draw_group_canvases(group, TOPO_STRS[topo_idx], topo_dir_name(topo_idx), csv, cols, out_root_dir);
        }
    }

    return true;
}