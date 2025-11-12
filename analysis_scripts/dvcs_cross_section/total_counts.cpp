// total_counts.cpp — self-contained (no csv_io.h), exact-typed binder, auto-mkdir, hierarchical bin index (supports non-uniform phi bins)

#include "total_counts.h"
#include "periods.h"                // CANONICAL_PERIODS(), PeriodDef{label, tree_key}
#include "load_binning_scheme.h"    // harmless include

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TROOT.h>
#include <TError.h>
#include <TLeaf.h>
#include <TBranch.h>
#include <TTree.h>
#include <TColor.h>
#include <TTreeCache.h>

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
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <nlohmann/json.hpp>

namespace {

// ---------------------- Minimal CSV helper (self-contained) ----------------------
struct CsvDoc {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;

    static std::vector<std::string> split_csv_line(const std::string& line) {
        std::vector<std::string> out;
        std::string cur;
        bool inq = false;
        for (char c : line) {
            if (c == '"') { inq = !inq; continue; }
            if (c == ',' && !inq) { out.push_back(cur); cur.clear(); }
            else cur.push_back(c);
        }
        out.push_back(cur);
        return out;
    }

    static void write_field(std::ostream& os, const std::string& s) {
        bool needq = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
        if (!needq) { os << s; return; }
        os << '"';
        for (char ch : s) {
            if (ch == '"') os << "\"\"";
            else os << ch;
        }
        os << '"';
    }

    bool load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin.is_open()) {
            std::cerr << "[CsvDoc] FATAL: cannot open CSV: " << path << std::endl;
            return false;
        }
        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[CsvDoc] FATAL: empty CSV: " << path << std::endl;
            return false;
        }
        header = split_csv_line(line);
        index.clear();
        for (int i = 0; i < (int)header.size(); ++i) index[header[i]] = i;

        rows.clear();
        while (std::getline(fin, line)) {
            if (!line.empty()) rows.push_back(split_csv_line(line));
        }
        for (auto& r : rows) r.resize(header.size());
        return true;
    }

    bool save(const std::string& path) const {
        std::ofstream fout(path);
        if (!fout.is_open()) {
            std::cerr << "[CsvDoc] FATAL: cannot write CSV: " << path << std::endl;
            return false;
        }
        for (size_t i = 0; i < header.size(); ++i) {
            write_field(fout, header[i]);
            if (i + 1 < header.size()) fout << ',';
        }
        fout << "\n";
        for (const auto& row : rows) {
            for (size_t i = 0; i < row.size(); ++i) {
                write_field(fout, row[i]);
                if (i + 1 < row.size()) fout << ',';
            }
            fout << "\n";
        }
        return true;
    }

    int nrows() const { return (int)rows.size(); }

    bool has_col(const std::string& name) const {
        return index.find(name) != index.end();
    }

    int col_index(const std::string& name) const {
        auto it = index.find(name);
        return (it == index.end()) ? -1 : it->second;
    }

    static double toD(const std::string& s) {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* e = nullptr;
        double v = std::strtod(s.c_str(), &e);
        if (e == s.c_str()) return std::numeric_limits<double>::quiet_NaN();
        return v;
    }

    double as_double(int r, int c) const {
        if (r < 0 || r >= (int)rows.size() || c < 0 || c >= (int)header.size()) return std::numeric_limits<double>::quiet_NaN();
        return toD(rows[r][c]);
    }

    void set(int r, int c, double v) {
        if (r < 0 || r >= (int)rows.size() || c < 0 || c >= (int)header.size()) return;
        std::ostringstream oss; oss << std::setprecision(16) << std::fixed << v;
        rows[r][c] = oss.str();
    }
};

[[noreturn]] void fatal(const std::string& msg) {
    std::cerr << "[total_counts] FATAL: " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline double pi() { return std::acos(-1.0); }
static inline double wrap_deg_0_360(double deg) {
    if (!std::isfinite(deg)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(deg, 360.0);
    if (w < 0.0) w += 360.0;
    if (w >= 360.0) w = std::nextafter(360.0, 0.0);
    return w;
}

static inline double rad_to_deg(double rad) { return rad * 180.0 / pi(); }
static inline bool in_range(double v, double a, double b) { return (v >= a) && (v < b); }

static bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (!std::isfinite(phi_deg) || !std::isfinite(pmin_deg) || !std::isfinite(pmax_deg)) return false;
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg); // wrap window
}

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
        const std::string mode = j.value("mode", std::string("two_sided"));
        c.mode    = parseMode(mode);
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
        if (c.sigma  <= 0.0) fatal("sigma_cuts: non-positive sigma for '" + c.branch + "'");
        if (c.nsigma <= 0.0) fatal("sigma_cuts: non-positive nsigma for '" + c.branch + "'");
    }
    return cuts;
}

// -------------- binding helpers (exact type) --------------
struct BranchBinding {
    std::string name;
    enum class Kind { kDouble, kFloat, kI32, kU32, kI64, kU64, kI16, kU16, kI8, kU8 } kind = Kind::kDouble;

    double              as_double = std::numeric_limits<double>::quiet_NaN();
    float               as_float  = std::numeric_limits<float>::quiet_NaN();
    int                 as_i32    = 0;
    unsigned int        as_u32    = 0;
    long long           as_i64    = 0;
    unsigned long long  as_u64    = 0;
    short               as_i16    = 0;
    unsigned short      as_u16    = 0;
    signed char         as_i8     = 0;
    unsigned char       as_u8     = 0;
};

static inline bool kind_is_int(BranchBinding::Kind k) {
    using K = BranchBinding::Kind;
    return k == K::kI32 || k == K::kU32 || k == K::kI64 || k == K::kU64 ||
           k == K::kI16 || k == K::kU16 || k == K::kI8  || k == K::kU8;
}

static inline double bb_as_double(const BranchBinding& bb) {
    using K = BranchBinding::Kind;
    switch (bb.kind) {
        case K::kDouble: return bb.as_double;
        case K::kFloat:  return (double)bb.as_float;
        case K::kI32:    return (double)bb.as_i32;
        case K::kU32:    return (double)bb.as_u32;
        case K::kI64:    return (double)bb.as_i64;
        case K::kU64:    return (double)bb.as_u64;
        case K::kI16:    return (double)bb.as_i16;
        case K::kU16:    return (double)bb.as_u16;
        case K::kI8:     return (double)bb.as_i8;
        case K::kU8:     return (double)bb.as_u8;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

static inline long long bb_as_ll(const BranchBinding& bb) {
    using K = BranchBinding::Kind;
    switch (bb.kind) {
        case K::kDouble: return (long long)std::llround(bb.as_double);
        case K::kFloat:  return (long long)std::llround(bb.as_float);
        case K::kI32:    return (long long)bb.as_i32;
        case K::kU32:    return (long long)bb.as_u32;
        case K::kI64:    return (long long)bb.as_i64;
        case K::kU64:    return (long long)bb.as_u64;
        case K::kI16:    return (long long)bb.as_i16;
        case K::kU16:    return (long long)bb.as_u16;
        case K::kI8:     return (long long)bb.as_i8;
        case K::kU8:     return (long long)bb.as_u8;
    }
    return 0;
}

static void bind_one_exact(TTree* t, const std::string& bname, BranchBinding& bb) {
    TBranch* b = t->GetBranch(bname.c_str());
    if (!b) return;
    TLeaf* leaf = b->GetLeaf(bname.c_str());
    if (!leaf) {
        leaf = (TLeaf*)b->GetListOfLeaves()->First();
        if (!leaf) fatal("Branch has no leaves: '" + bname + "'");
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

static void bindRequiredBranches_STRICT(
    TTree* t,
    const std::vector<std::string>& branch_names,
    std::unordered_map<std::string, BranchBinding>& bindings)
{
    bindings.clear();
    bindings.reserve(branch_names.size());
    for (const auto& bname : branch_names) {
        TBranch* b = t->GetBranch(bname.c_str());
        if (!b) continue; // non-existing may be optional in cuts
        auto [it, inserted] = bindings.emplace(bname, BranchBinding{});
        BranchBinding& bb = it->second;
        bb.name = bname;
        bind_one_exact(t, bname, bb);
    }
}

static inline bool passBaseCuts(const std::vector<BaseCut>& baseCuts,
                                const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : baseCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false;
        const BranchBinding& bb = it->second;

        if (c.has_eq)  { if (bb_as_ll(bb) != c.eq)  return false; }
        if (c.has_neq) { if (bb_as_ll(bb) == c.neq) return false; }
        if (c.has_min) { if (!(bb_as_double(bb) >= c.vmin)) return false; }
        if (c.has_max) { if (!(bb_as_double(bb) <= c.vmax)) return false; }
    }
    return true;
}

static inline bool passSigmaCuts(const std::vector<SigmaCut>& sigmaCuts,
                                 const std::unordered_map<std::string, BranchBinding>& B) {
    for (const auto& c : sigmaCuts) {
        const auto it = B.find(c.branch);
        if (it == B.end()) return false;
        const BranchBinding& bb = it->second;
        const double v  = bb_as_double(bb);
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
static const char* TOPO_DIRS[3] = {"FD_FD","CD_FD","CD_FT"};

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

    int index() const {
        if (detector1 == 1 && detector2 == 1) return 0; // FD,FD
        if (detector1 == 2 && detector2 == 1) return 1; // CD,FD
        if (detector1 == 2 && detector2 == 0) return 2; // CD,FT
        return -1;
    }
};

// -------------- CSV helpers for counts and plotting --------------
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

static void ensure_plot_dir(const std::string& path) {
    namespace fs = std::filesystem;
    std::error_code ec;
    if (!fs::exists(path)) {
        fs::create_directories(path, ec);
        if (ec) fatal("Could not create plot output directory: " + path + " (" + ec.message() + ")");
    }
}

// Draw one label+topology set of canvases (one canvas per xB range)
static void draw_group_canvases(
    const std::string& label,
    const std::string& topo_str,
    const std::string& topo_dir,
    const CsvDoc& csv,
    const CsvCols& cols,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const std::string out_dir = (fs::path(out_root_dir) / "total_counts_plots" / label).string();
    ensure_plot_dir(out_dir);

    const std::string c_phiavg_name = std::string("phiavg, ") + label;
    const std::string c_q2avg_name  = std::string("Q2avg, ") + label;
    const std::string c_tabavg_name = std::string("t_abs_avg, ") + label;
    const std::string c_xbavg_name  = std::string("xBavg, ") + label;

    const int c_phiavg = csv.has_col(c_phiavg_name) ? csv.col_index(c_phiavg_name) : -1;
    const int c_q2avg  = csv.has_col(c_q2avg_name)  ? csv.col_index(c_q2avg_name)  : -1;
    const int c_tabavg = csv.has_col(c_tabavg_name) ? csv.col_index(c_tabavg_name) : -1;
    const int c_xbavg  = csv.has_col(c_xbavg_name)  ? csv.col_index(c_xbavg_name)  : -1;

    const std::string ybase = period_yield_col_base(label, topo_str);
    const int c_pos = csv.has_col(ybase + "pos")   ? csv.col_index(ybase + "pos")   : -1;
    const int c_neg = csv.has_col(ybase + "neg")   ? csv.col_index(ybase + "neg")   : -1;

    std::set<std::pair<double,double>> xb_set;
    for (int r = 0; r < csv.nrows(); ++r) {
        xb_set.emplace(csv.as_double(r, cols.c_xb_min), csv.as_double(r, cols.c_xb_max));
    }

    for (auto xb : xb_set) {
        std::set<std::pair<double,double>> q2set;
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, cols.c_xb_min);
            const double xbmax = csv.as_double(r, cols.c_xb_max);
            if (std::fabs(xbmin - xb.first) < 1e-9 && std::fabs(xbmax - xb.second) < 1e-9) {
                q2set.emplace(csv.as_double(r, cols.c_q2_min), csv.as_double(r, cols.c_q2_max));
            }
        }

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

                std::vector<int> rows_for_cell;
                rows_for_cell.reserve(24);
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
                        rows_for_cell.push_back(r);
                    }
                }
                if (rows_for_cell.empty()) {
                    TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 1.0);
                    frame->GetXaxis()->SetTitle("phi (deg)");
                    frame->GetYaxis()->SetTitle("Counts");
                    continue;
                }

                std::sort(rows_for_cell.begin(), rows_for_cell.end(), [&](int a, int b){
                    return csv.as_double(a, cols.c_phi_min) < csv.as_double(b, cols.c_phi_min);
                });

                const std::string ybase = period_yield_col_base(label, topo_str);
                const int c_pos = csv.has_col(ybase + "pos") ? csv.col_index(ybase + "pos") : -1;
                const int c_neg = csv.has_col(ybase + "neg") ? csv.col_index(ybase + "neg") : -1;

                std::vector<double> X, Yp, Ym, EX, EYp, EYm, q2means, tmeans;
                X.reserve(rows_for_cell.size());
                EX.assign(rows_for_cell.size(), 0.0);

                for (int r : rows_for_cell) {
                    const double pmin = csv.as_double(r, cols.c_phi_min);
                    const double pmax = csv.as_double(r, cols.c_phi_max);
                    double xphi = 0.5 * (pmin + pmax); // default center
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

// Enable only the branches we will actually read on this thread.
static void enable_needed_branches(TTree* t,
                                   const std::vector<std::string>& baseline,
                                   const std::vector<std::string>& extra)
{
    t->SetBranchStatus("*", 0);
    for (const auto& b : baseline) if (t->GetBranch(b.c_str())) t->SetBranchStatus(b.c_str(), 1);
    for (const auto& b : extra)    if (t->GetBranch(b.c_str())) t->SetBranchStatus(b.c_str(), 1);
}

// --------- Hierarchical bin indexer (xB -> Q2 -> |t| -> vector of phi bins) ----------
struct CsvCols;
struct BinIndex {
    struct PhiBin { double pmin, pmax; int row; };
    struct TLevel {
        std::vector<double> t_edges;              // size Nt+1
        std::vector<std::vector<PhiBin>> phi;     // size Nt, phi[i] is the list for t-bin i
    };
    struct QLevel {
        std::vector<double> q2_edges;             // size Nq+1
        std::vector<TLevel> t_levels;             // size Nq
    };

    std::vector<double> xb_edges;                 // size Nx+1
    std::vector<QLevel> q_levels;                 // size Nx

    static int bin_index(double v, const std::vector<double>& edges) {
        auto it = std::upper_bound(edges.begin(), edges.end(), v);
        int i = int(it - edges.begin()) - 1;
        if (i < 0 || i+1 >= (int)edges.size()) return -1;
        return i;
    }

    static std::vector<double> collect_edges_for_axis(const CsvDoc& csv, int cmin, int cmax) {
        std::set<double> s;
        for (int r = 0; r < csv.nrows(); ++r) {
            s.insert(csv.as_double(r, cmin));
            s.insert(csv.as_double(r, cmax));
        }
        return std::vector<double>(s.begin(), s.end());
    }

    void build(const CsvDoc& csv, const CsvCols& cols) {
        // Global xB edges
        xb_edges = collect_edges_for_axis(csv, cols.c_xb_min, cols.c_xb_max);
        if (xb_edges.size() < 2) fatal("BinIndex: not enough xB edges.");
        const int Nx = (int)xb_edges.size() - 1;
        q_levels.assign(Nx, QLevel{});

        // For each xB-bin, gather Q2 edges from rows that start in this xB-bin
        std::vector<std::set<double>> q2_sets(Nx);
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, cols.c_xb_min);
            int ix = bin_index(xbmin + 1e-12, xb_edges);
            if (ix < 0) continue;
            q2_sets[ix].insert(csv.as_double(r, cols.c_q2_min));
            q2_sets[ix].insert(csv.as_double(r, cols.c_q2_max));
        }
        for (int ix = 0; ix < Nx; ++ix) {
            if (q2_sets[ix].size() < 2) fatal("BinIndex: not enough Q2 edges in an xB slice.");
            q_levels[ix].q2_edges.assign(q2_sets[ix].begin(), q2_sets[ix].end());
            const int Nq = (int)q_levels[ix].q2_edges.size() - 1;
            q_levels[ix].t_levels.assign(Nq, TLevel{});
        }

        // For each (xB,Q2) bin, gather t edges from rows in that slice
        std::vector<std::vector<std::set<double>>> t_sets;
        t_sets.resize(Nx);
        for (int ix = 0; ix < Nx; ++ix) {
            const int Nq = (int)q_levels[ix].t_levels.size();
            t_sets[ix].resize(Nq);
        }
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, cols.c_xb_min);
            const double q2min = csv.as_double(r, cols.c_q2_min);
            int ix = bin_index(xbmin + 1e-12, xb_edges);
            if (ix < 0) continue;
            int iq = bin_index(q2min + 1e-12, q_levels[ix].q2_edges);
            if (iq < 0) continue;
            t_sets[ix][iq].insert(csv.as_double(r, cols.c_tab_min));
            t_sets[ix][iq].insert(csv.as_double(r, cols.c_tab_max));
        }
        for (int ix = 0; ix < Nx; ++ix) {
            int Nq = (int)q_levels[ix].t_levels.size();
            for (int iq = 0; iq < Nq; ++iq) {
                if (t_sets[ix][iq].size() < 2) fatal("BinIndex: not enough |t| edges in an (xB,Q2) slice.");
                q_levels[ix].t_levels[iq].t_edges.assign(t_sets[ix][iq].begin(), t_sets[ix][iq].end());
                int Nt = (int)q_levels[ix].t_levels[iq].t_edges.size() - 1;
                q_levels[ix].t_levels[iq].phi.assign(Nt, std::vector<PhiBin>{});
            }
        }

        // Distribute each row into its (ix,iq,it) cell and store phi interval
        for (int r = 0; r < csv.nrows(); ++r) {
            const double xbmin = csv.as_double(r, cols.c_xb_min);
            const double q2min = csv.as_double(r, cols.c_q2_min);
            const double tmin  = csv.as_double(r, cols.c_tab_min);
            const double pmin  = csv.as_double(r, cols.c_phi_min);
            const double pmax  = csv.as_double(r, cols.c_phi_max);

            int ix = bin_index(xbmin + 1e-12, xb_edges);
            if (ix < 0) continue;
            int iq = bin_index(q2min + 1e-12, q_levels[ix].q2_edges);
            if (iq < 0) continue;
            int it = bin_index(tmin + 1e-12, q_levels[ix].t_levels[iq].t_edges);
            if (it < 0) continue;

            q_levels[ix].t_levels[iq].phi[it].push_back(PhiBin{pmin, pmax, r});
        }

        // Optional: sanity check that each (ix,iq,it) has at least one phi bin
        // Not strictly required if some cells are intentionally empty.
    }

    int find_row(double xb, double q2, double tab, double phi_deg) const {
        int ix = bin_index(xb, xb_edges);
        if (ix < 0) return -1;
        const QLevel& ql = q_levels[ix];

        int iq = bin_index(q2, ql.q2_edges);
        if (iq < 0) return -1;
        const TLevel& tl = ql.t_levels[iq];

        int it = bin_index(tab, tl.t_edges);
        if (it < 0) return -1;
        const std::vector<PhiBin>& pb = tl.phi[it];
        if (pb.empty()) return -1;

        // Non-uniform phi bins: small linear scan over the few phi bins in this cell.
        for (const auto& bin : pb) {
            if (row_accepts_phi(phi_deg, bin.pmin, bin.pmax)) return bin.row;
        }
        return -1;
    }
};

// ---------------- entry ----------------
} // namespace

bool update_total_counts_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_root_dir,
    int max_workers)
{
    namespace fs = std::filesystem;

    std::vector<std::string> periods;
    periods.reserve(CANONICAL_PERIODS().size());
    for (const auto& P : CANONICAL_PERIODS()) {
        periods.push_back(P.tree_key);
    }

    CsvDoc csv;
    if (!csv.load(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to read " << csv_path << "\n";
        return false;
    }

    CsvCols cols;
    cols.c_xb_min  = csv.col_index("xBmin");
    cols.c_xb_max  = csv.col_index("xBmax");
    cols.c_q2_min  = csv.col_index("Q2min");
    cols.c_q2_max  = csv.col_index("Q2max");
    cols.c_tab_min = csv.col_index("t_abs_min");
    cols.c_tab_max = csv.col_index("t_abs_max");
    cols.c_phi_min = csv.col_index("phimin");
    cols.c_phi_max = csv.col_index("phimax");
    if (cols.c_xb_min < 0 || cols.c_xb_max < 0 || cols.c_q2_min < 0 || cols.c_q2_max < 0
        || cols.c_tab_min < 0 || cols.c_tab_max < 0 || cols.c_phi_min < 0 || cols.c_phi_max < 0) {
        fatal("Missing one or more required bin-edge columns in CSV.");
    }

    // Build hierarchical bin index once (supports non-uniform phi bins per cell)
    BinIndex indexer;
    indexer.build(csv, cols);

    auto ensure_yield_columns = [&](const std::string& label, const std::string& topo_str) {
        const std::string base = period_yield_col_base(label, topo_str);
        csv.ensure_col(base + "unpol");
        csv.ensure_col(base + "pos");
        csv.ensure_col(base + "neg");
    };

    using PeriodRowMap = std::unordered_map<int, RowCounts>;
    std::mutex merge_mtx;
    std::map<std::string, std::map<std::string, PeriodRowMap>> counts_by_label_topo;

    ROOT::EnableThreadSafety(); (void)TGraphErrors::Class();

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

        // Build and de-dup 'need'
        std::vector<std::string> need;
        need.reserve(cuts.base_cuts.size() + cuts.sigma_cuts.size() + 8);
        for (const auto& c : cuts.base_cuts)  need.push_back(c.branch);
        for (const auto& c : cuts.sigma_cuts) need.push_back(c.branch);

        // Explicitly bound branches (exclude these from generic binder)
        const char* explicit_binds[] = {"helicity","x","Q2","t1","phi2","detector1","detector2"};
        for (const char* s : explicit_binds) need.push_back(s);

        std::sort(need.begin(), need.end());
        need.erase(std::unique(need.begin(), need.end()), need.end());

        // Enable only what's needed
        enable_needed_branches(t, {}, need);

        // ROOT I/O cache: cache only what we read on this thread
        t->SetCacheSize(128*1024*1024);
        TTreeCache::SetLearnEntries(10);
        t->AddBranchToCache("helicity", true);
        t->AddBranchToCache("x",        true);
        t->AddBranchToCache("Q2",       true);
        t->AddBranchToCache("t1",       true);
        t->AddBranchToCache("phi2",     true);
        t->AddBranchToCache("detector1",true);
        t->AddBranchToCache("detector2",true);

        // Explicit exact-typed bindings for these known branches
        int helicity = 0;
        double x = 0.0, Q2 = 0.0, t1 = 0.0, phi2 = std::numeric_limits<double>::quiet_NaN();
        if (!t->GetBranch("helicity") || !t->GetBranch("x") || !t->GetBranch("Q2")
            || !t->GetBranch("t1") || !t->GetBranch("phi2")) {
            #pragma omp critical
            { std::cerr << "[total_counts] Required branches missing in " << period_key << "\n"; }
            continue;
        }
        t->SetBranchAddress("helicity", &helicity); // Int_t
        t->SetBranchAddress("x",        &x);        // Double_t
        t->SetBranchAddress("Q2",       &Q2);       // Double_t
        t->SetBranchAddress("t1",       &t1);       // Double_t
        t->SetBranchAddress("phi2",     &phi2);     // Double_t (radians)

        TopologyResolver topo;
        topo.bind(t); // binds detector1, detector2 as int

        // Build the list for the generic binder but DROP explicitly bound names
        std::vector<std::string> need_for_binder;
        need_for_binder.reserve(need.size());
        for (const auto& b : need) {
            bool drop = false;
            for (const char* s : explicit_binds) if (b == s) { drop = true; break; }
            if (!drop) need_for_binder.push_back(b);
        }

        std::unordered_map<std::string, BranchBinding> bind;
        bindRequiredBranches_STRICT(t, need_for_binder, bind);
        // Add those branches to cache too
        for (const auto& kv : bind) t->AddBranchToCache(kv.first.c_str(), true);

        std::map<std::string, PeriodRowMap> local_by_topo;

        const Long64_t nent = t->GetEntries();
        std::cout << "[total_counts] " << label << " entries=" << nent
                  << "  (cuts: base=" << cuts.base_cuts.size()
                  << ", sigma=" << cuts.sigma_cuts.size() << ")\n";

        auto t0 = std::chrono::steady_clock::now();
        const Long64_t report_every = 500000;

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

            const int row = indexer.find_row(xb, q2, tab, phi_deg);
            if (row < 0) continue;

            RowCounts& rc = local_by_topo[topo_str][row];
            if (helicity == +1) rc.pos++; else rc.neg++;

            if ((i % report_every) == 0) {
                auto t1clk = std::chrono::steady_clock::now();
                double sec = std::chrono::duration<double>(t1clk - t0).count();
                double rate = (sec > 0.0) ? (i / sec) : 0.0;
                std::cout << "[total_counts] " << label << " "
                          << i << "/" << nent << " ("
                          << std::fixed << std::setprecision(1)
                          << (100.0 * double(i)/double(nent)) << "%), "
                          << std::setprecision(2) << rate/1e6 << " Mev/s\n";
            }
        }

        {
            std::lock_guard<std::mutex> lock(merge_mtx);
            for (const auto& tk : local_by_topo) {
                const std::string topo_str = tk.first;
                ensure_yield_columns(label, topo_str);

                auto& tgt = counts_by_label_topo[label][topo_str];
                for (const auto& rkv : tk.second) {
                    tgt[rkv.first].pos += rkv.second.pos;
                    tgt[rkv.first].neg += rkv.second.neg;
                }
            }
        }
    } // omp parallel loop

    // Write per-period counts into CSV
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

    // Build groups
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

    // Combined per-group counts
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

    if (!csv.save(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to save " << csv_path << "\n";
        return false;
    }
    std::cout << "[total_counts] Updated raw yields in: " << csv_path << "\n";

    // Plots (auto-creates output dirs)
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