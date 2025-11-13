// total_counts.cpp — raw-yield writer + plotting
// Aesthetics tuned (left pad margin, title size, legend), uses Q^{2} in TLatex,
// canonicalized output dirs to match make_dirs.cpp, and strong CSV diagnostics.

#include "total_counts.h"
#include "periods.h"
#include "load_binning_scheme.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TError.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TLeafElement.h>

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
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <nlohmann/json.hpp>

namespace {

// ---------------- small utils ----------------
static inline bool is_debug() { return std::getenv("TOTAL_COUNTS_DEBUG") != nullptr; }
static inline bool trace_matches() { return std::getenv("TOTAL_COUNTS_TRACE_MATCHES") != nullptr; }
static inline bool list_enabled_branches() { return std::getenv("TOTAL_COUNTS_LIST_ENABLED_BRANCHES") != nullptr; }

[[noreturn]] void fatal(const std::string& msg) {
    std::cerr << "[total_counts] FATAL: " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline double pi() { return std::acos(-1.0); }
static inline double rad_to_deg(double r) { return r * 180.0 / pi(); }
static inline double wrap_deg_0_360(double d) {
    if (!std::isfinite(d)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(d, 360.0);
    if (w < 0.0) w += 360.0;
    if (w >= 360.0) w = std::nextafter(360.0, 0.0);
    return w;
}
static inline bool in_range(double v, double a, double b) { return (v >= a) && (v < b); }

static bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (!std::isfinite(phi_deg) || !std::isfinite(pmin_deg) || !std::isfinite(pmax_deg)) return false;
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg); // wrap
}

// ---------------- CSV helper ----------------
struct CsvDoc {
    std::vector<std::string> header;
    std::unordered_map<std::string,int> index;
    std::vector<std::vector<std::string>> rows;

    static std::vector<std::string> split_csv_line(const std::string& line) {
        std::vector<std::string> out; std::string cur; bool inq=false;
        for (char c : line) {
            if (c == '"') { inq = !inq; continue; }
            if (c == ',' && !inq) { out.push_back(cur); cur.clear(); }
            else cur.push_back(c);
        }
        out.push_back(cur); return out;
    }

    static void write_field(std::ostream& os, const std::string& s) {
        bool needq = s.find(',') != std::string::npos || s.find('"') != std::string::npos;
        if (!needq) { os << s; return; }
        os << '"';
        for (char ch : s) os << (ch=='"'? "\"\"" : std::string(1, ch));
        os << '"';
    }

    static double toD(const std::string& s) {
        if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
        char* e=nullptr; double v = std::strtod(s.c_str(), &e);
        return (e==s.c_str()) ? std::numeric_limits<double>::quiet_NaN() : v;
    }

    bool load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin.is_open()) { std::cerr << "[total_counts] cannot open CSV: " << path << "\n"; return false; }
        std::string line;
        if (!std::getline(fin, line)) { std::cerr << "[total_counts] empty CSV: " << path << "\n"; return false; }
        header = split_csv_line(line);
        index.clear();
        for (int i=0;i<(int)header.size();++i) index[header[i]] = i;
        rows.clear();
        while (std::getline(fin, line)) if (!line.empty()) rows.push_back(split_csv_line(line));
        for (auto& r : rows) r.resize(header.size());
        return true;
    }

    bool save_atomic(const std::string& path) const {
        const std::string tmp = path + ".tmp";
        {
            std::ofstream fout(tmp);
            if (!fout.is_open()) { std::cerr << "[total_counts] cannot write CSV tmp: " << tmp << "\n"; return false; }
            // header
            for (size_t i=0;i<header.size();++i) { write_field(fout, header[i]); if (i+1<header.size()) fout<<','; }
            fout << "\n";
            // rows
            for (const auto& row : rows) {
                for (size_t i=0;i<row.size();++i) { write_field(fout, row[i]); if (i+1<row.size()) fout<<','; }
                fout << "\n";
            }
        }
        std::error_code ec;
        std::filesystem::rename(tmp, path, ec);
        if (ec) {
            std::remove(path.c_str());
            std::filesystem::rename(tmp, path, ec);
            if (ec) {
                std::cerr << "[total_counts] ERROR: atomic rename failed (" << ec.message() << ")\n";
                return false;
            }
        }
        return true;
    }

    int nrows() const { return (int)rows.size(); }

    bool has_col(const std::string& name) const { return index.find(name)!=index.end(); }
    int col_index(const std::string& name) const {
        auto it=index.find(name); return (it==index.end())? -1 : it->second;
    }
    int ensure_col(const std::string& name) {
        auto it=index.find(name);
        if (it!=index.end()) return it->second;
        int new_idx=(int)header.size();
        header.push_back(name); index[name]=new_idx;
        for (auto& r:rows) r.resize(header.size());
        if (is_debug()) std::cout << "[total_counts] ensured CSV column: \"" << name << "\" @ " << new_idx << "\n";
        return new_idx;
    }
    double as_double(int r, int c) const {
        if (r<0 || r>=nrows() || c<0 || c>=(int)header.size()) return std::numeric_limits<double>::quiet_NaN();
        return toD(rows[r][c]);
    }
    void set(int r, int c, double v) {
        if (r<0 || r>=nrows() || c<0 || c>=(int)header.size()) return;
        std::ostringstream oss; oss.setf(std::ios::fixed); oss<<std::setprecision(8)<<v;
        rows[r][c]=oss.str();
    }
};

// ---------------- cuts config ----------------
struct SigmaCut { std::string branch; double center=0, sigma=0, nsigma=3; enum Mode{TWO_SIDED,UPPER,LOWER} mode=TWO_SIDED; };
struct BaseCut  { std::string branch; bool has_min=false,has_max=false,has_eq=false,has_neq=false; double vmin=0,vmax=0; long long eq=0,neq=0; };
struct PeriodCuts { std::vector<SigmaCut> sigma_cuts; std::vector<BaseCut> base_cuts; };

static SigmaCut::Mode parseMode(const std::string& m) {
    if (m=="two_sided") return SigmaCut::TWO_SIDED;
    if (m=="upper")     return SigmaCut::UPPER;
    if (m=="lower")     return SigmaCut::LOWER;
    fatal("Unknown sigma cut mode: "+m); return SigmaCut::TWO_SIDED;
}
static void parseSigmaCuts(const nlohmann::json& arr, std::vector<SigmaCut>& out) {
    if (!arr.is_array()) fatal("sigma_cuts must be array");
    for (const auto& j:arr) {
        SigmaCut c; c.branch=j.at("branch").get<std::string>();
        c.center=j.at("center").get<double>(); c.sigma=j.at("sigma").get<double>();
        c.nsigma=j.value("nsigma",3.0); c.mode=parseMode(j.value("mode",std::string("two_sided")));
        out.push_back(c);
    }
}
static void parseBaseCuts(const nlohmann::json& arr, std::vector<BaseCut>& out) {
    if (!arr.is_array()) fatal("base_cuts must be array");
    for (const auto& j:arr) {
        BaseCut c; c.branch=j.at("branch").get<std::string>();
        if (j.contains("min")) { c.has_min=true; c.vmin=j.at("min").get<double>(); }
        if (j.contains("max")) { c.has_max=true; c.vmax=j.at("max").get<double>(); }
        if (j.contains("eq"))  { c.has_eq =true; c.eq  =j.at("eq").get<long long>(); }
        if (j.contains("neq")) { c.has_neq=true; c.neq =j.at("neq").get<long long>(); }
        if (!(c.has_min||c.has_max||c.has_eq||c.has_neq)) fatal("base_cuts entry for '"+c.branch+"' has no condition");
        out.push_back(c);
    }
}
static PeriodCuts loadCombinedCuts(const std::string& json_path, const std::string& period_key) {
    std::ifstream ifs(json_path);
    if (!ifs) fatal("Cannot open combined cuts JSON: "+json_path);
    nlohmann::json J; ifs >> J;
    PeriodCuts cuts;
    if (J.contains("global")) {
        const auto& G=J.at("global");
        if (G.contains("sigma_cuts")) parseSigmaCuts(G.at("sigma_cuts"), cuts.sigma_cuts);
        if (G.contains("base_cuts"))  parseBaseCuts(G.at("base_cuts"),  cuts.base_cuts);
    }
    if (J.contains("period_overrides")) {
        const auto& PO=J.at("period_overrides");
        auto it=PO.find(period_key);
        if (it!=PO.end()) {
            const auto& P=*it;
            if (P.contains("sigma_cuts")) { cuts.sigma_cuts.clear(); parseSigmaCuts(P.at("sigma_cuts"), cuts.sigma_cuts); }
            if (P.contains("base_cuts"))  { cuts.base_cuts.clear();  parseBaseCuts(P.at("base_cuts"),  cuts.base_cuts); }
        }
    }
    for (const auto& c:cuts.sigma_cuts) {
        if (c.sigma<=0.0) fatal("sigma_cuts: non-positive sigma for '"+c.branch+"'");
        if (c.nsigma<=0.0) fatal("sigma_cuts: non-positive nsigma for '"+c.branch+"'");
    }
    return cuts;
}

// ---------------- branch binding ----------------
struct BranchBinding {
    std::string name;
    enum class Kind { kDouble,kFloat,kI32,kU32,kI64,kU64,kI16,kU16,kI8,kU8 } kind = Kind::kDouble;
    double as_double = std::numeric_limits<double>::quiet_NaN();
    float  as_float  = std::numeric_limits<float>::quiet_NaN();
    int as_i32=0; unsigned as_u32=0; long long as_i64=0; unsigned long long as_u64=0;
    short as_i16=0; unsigned short as_u16=0; signed char as_i8=0; unsigned char as_u8=0;
};
static inline double bb_as_double(const BranchBinding& bb){
    using K=BranchBinding::Kind;
    switch(bb.kind){
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
static inline long long bb_as_ll(const BranchBinding& bb){
    using K=BranchBinding::Kind;
    switch(bb.kind){
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
    TBranch* b=t->GetBranch(bname.c_str());
    if (!b) return;
    TLeaf* leaf=b->GetLeaf(bname.c_str());
    if (!leaf) { leaf=(TLeaf*)b->GetListOfLeaves()->First(); if (!leaf) fatal("Branch has no leaves: "+bname); }
    const std::string tn=leaf->GetTypeName();
    if (tn=="Double_t"||tn=="double") { bb.kind=BranchBinding::Kind::kDouble; t->SetBranchAddress(bname.c_str(), &bb.as_double); }
    else if (tn=="Float_t"||tn=="float") { bb.kind=BranchBinding::Kind::kFloat; t->SetBranchAddress(bname.c_str(), &bb.as_float); }
    else if (tn=="Int_t"||tn=="int") { bb.kind=BranchBinding::Kind::kI32; t->SetBranchAddress(bname.c_str(), &bb.as_i32); }
    else if (tn=="UInt_t"||tn=="unsigned int") { bb.kind=BranchBinding::Kind::kU32; t->SetBranchAddress(bname.c_str(), &bb.as_u32); }
    else if (tn=="Long64_t"||tn=="long long") { bb.kind=BranchBinding::Kind::kI64; t->SetBranchAddress(bname.c_str(), &bb.as_i64); }
    else if (tn=="ULong64_t"||tn=="unsigned long long") { bb.kind=BranchBinding::Kind::kU64; t->SetBranchAddress(bname.c_str(), &bb.as_u64); }
    else if (tn=="Short_t"||tn=="short") { bb.kind=BranchBinding::Kind::kI16; t->SetBranchAddress(bname.c_str(), &bb.as_i16); }
    else if (tn=="UShort_t"||tn=="unsigned short") { bb.kind=BranchBinding::Kind::kU16; t->SetBranchAddress(bname.c_str(), &bb.as_u16); }
    else if (tn=="Char_t"||tn=="char"||tn=="signed char") { bb.kind=BranchBinding::Kind::kI8; t->SetBranchAddress(bname.c_str(), &bb.as_i8); }
    else if (tn=="UChar_t"||tn=="unsigned char") { bb.kind=BranchBinding::Kind::kU8; t->SetBranchAddress(bname.c_str(), &bb.as_u8); }
    else { fatal("Unsupported leaf type '"+tn+"' for branch '"+bname+"'"); }
}
static void bindRequiredBranches_STRICT(TTree* t, const std::vector<std::string>& names, std::unordered_map<std::string,BranchBinding>& out) {
    out.clear(); out.reserve(names.size());
    for (const auto& bname : names) {
        TBranch* b=t->GetBranch(bname.c_str());
        if (!b) continue;
        auto [it,ins]=out.emplace(bname, BranchBinding{});
        it->second.name=bname;
        bind_one_exact(t, bname, it->second);
    }
}

// ---------------- topo + labels ----------------
static const char* TOPO_STRS[3]={"(FD, FD)","(CD, FD)","(CD, FT)"};
static const char* TOPO_DIRS[3]={"FD_FD","CD_FD","CD_FT"};

struct TopologyResolver {
    int detector1=0, detector2=0; bool have1=false, have2=false;
    void bind(TTree* t) {
        have1=(t->GetBranch("detector1")!=nullptr);
        have2=(t->GetBranch("detector2")!=nullptr);
        if (!(have1 && have2)) fatal("Missing detector1/detector2 in DVCS tree.");
        t->SetBranchAddress("detector1", &detector1);
        t->SetBranchAddress("detector2", &detector2);
    }
    int index() const {
        if (detector1==1 && detector2==1) return 0;
        if (detector1==2 && detector2==1) return 1;
        if (detector1==2 && detector2==0) return 2;
        return -1;
    }
};

static inline bool is_canonical_tree_key(const std::string& k) {
    for (const auto& p:CANONICAL_PERIODS()) if (k==p.tree_key) return true;
    return false;
}
static inline std::string safe_label_for_key(const std::string& k) {
    for (const auto& p:CANONICAL_PERIODS()) if (k==p.tree_key) return std::string(p.label);
    fatal("Non-canonical period key: "+k); return {};
}
static inline bool is_group_label(const std::string& L) { return (L=="Fa18"||L=="Sp18"||L=="10.6 GeV"); }

// Canonical period directory names — must match make_dirs.cpp
static std::string period_dir_for_label(const std::string& L) {
    if (L=="Fa18 Inb")      return "Fa18_Inb";
    if (L=="Fa18 Out")      return "Fa18_Out";
    if (L=="Fa18 Inb Supp") return "Fa18_Inb_Supp";
    if (L=="Sp18 Inb")      return "Sp18_Inb";
    if (L=="Sp18 Out")      return "Sp18_Out";
    if (L=="Sp19 Inb")      return "Sp19_Inb";
    std::string s=L; for (char& c:s) if (c==' ') c='_';
    return s;
}

static void ensure_dir(const std::string& p) {
    namespace fs=std::filesystem;
    std::error_code ec;
    if (!fs::exists(p)) {
        fs::create_directories(p, ec);
        if (ec) fatal("Could not create output directory: "+p+" ("+ec.message()+")");
        if (is_debug()) std::cout << "[total_counts] created dir: " << p << "\n";
    }
}

// Resolve the base plot directory (canonical only).
static std::string out_root_for_label(const std::string& label, const std::string& out_root_dir) {
    namespace fs=std::filesystem;
    const std::string base = (fs::path(out_root_dir) / "total_counts_plots").string();
    if (is_group_label(label)) {
        // Groups are literal (Fa18, Sp18, 10.6 GeV).
        return (fs::path(base) / label).string();
    }
    return (fs::path(base) / period_dir_for_label(label)).string();
}

struct CsvCols {
    int c_xb_min, c_xb_max, c_q2_min, c_q2_max, c_tab_min, c_tab_max, c_phi_min, c_phi_max;
};

// ---------------- plotting ----------------
struct CellData { std::vector<double> X,Yp,Ym,EX,EYp,EYm,q2means,tmeans; };

static double safe_mean(const std::vector<double>& v){
    double s=0; int n=0; for(double x:v) if (std::isfinite(x)){ s+=x; ++n; }
    return n? s/n : std::numeric_limits<double>::quiet_NaN();
}

static void draw_group_canvases(
    const std::string& label,
    const std::string& topo_str,
    const std::string& topo_dir,
    const CsvDoc& csv,
    const CsvCols& cols,
    const std::string& out_root_dir)
{
    namespace fs=std::filesystem;
    const std::string base_dir = out_root_for_label(label, out_root_dir);
    ensure_dir(base_dir);
    const std::string out_dir = (fs::path(base_dir) / topo_dir).string();
    ensure_dir(out_dir);

    if (is_debug()) std::cout << "[total_counts] plotting to " << out_dir << " (label=" << label << ", topo=" << topo_str << ")\n";

    const int c_phiavg = csv.has_col("phiavg, "+label)? csv.col_index("phiavg, "+label) : -1;
    const int c_q2avg  = csv.has_col("Q2avg, " +label)? csv.col_index("Q2avg, " +label) : -1;
    const int c_tabavg = csv.has_col("t_abs_avg, "+label)? csv.col_index("t_abs_avg, "+label) : -1;
    const int c_xbavg  = csv.has_col("xBavg, "+label)? csv.col_index("xBavg, "+label) : -1;

    // collect unique xB ranges
    std::set<std::pair<double,double>> xb_set;
    for (int r=0;r<csv.nrows();++r) xb_set.emplace(csv.as_double(r,cols.c_xb_min), csv.as_double(r,cols.c_xb_max));

    // aesthetics: tiny bit more left margin; slightly smaller top title; legend bottom-left nudged left
    const double head_size = 0.58;   // was 0.62
    const double lab_size  = 0.075;
    const double tick_size = 0.060;
    const double latex_size= 0.065;

    for (auto xb : xb_set) {
        // Q2 and |t| sets for this xB
        std::set<std::pair<double,double>> q2set, tset_all;
        for (int r=0;r<csv.nrows();++r) {
            const double xbmin=csv.as_double(r, cols.c_xb_min), xbmax=csv.as_double(r, cols.c_xb_max);
            if (std::fabs(xbmin-xb.first)<1e-9 && std::fabs(xbmax-xb.second)<1e-9)
                q2set.emplace(csv.as_double(r, cols.c_q2_min), csv.as_double(r, cols.c_q2_max));
        }
        for (auto q2r : q2set) {
            for (int r=0;r<csv.nrows();++r) {
                const double xbmin=csv.as_double(r, cols.c_xb_min), xbmax=csv.as_double(r, cols.c_xb_max);
                const double q2min=csv.as_double(r, cols.c_q2_min), q2max=csv.as_double(r, cols.c_q2_max);
                if (std::fabs(xbmin-xb.first)<1e-9 && std::fabs(xbmax-xb.second)<1e-9 &&
                    std::fabs(q2min-q2r.first)<1e-9 && std::fabs(q2max-q2r.second)<1e-9) {
                    tset_all.emplace(csv.as_double(r, cols.c_tab_min), csv.as_double(r, cols.c_tab_max));
                }
            }
        }
        std::vector<std::pair<double,double>> Q2s(q2set.begin(), q2set.end());
        std::vector<std::pair<double,double>> Ts (tset_all.begin(), tset_all.end());
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows=(int)Q2s.size(), ncols=(int)Ts.size();
        const int W = 300*ncols + 160;
        const int H = 260*nrows + 240;

        // build cells + y-max
        std::vector<CellData> cells(nrows*ncols);
        double canvas_max = 1.0; bool need_log=false;

        std::vector<double> xbmeans;
        for (int r=0;r<csv.nrows();++r) {
            const double xbmin=csv.as_double(r, cols.c_xb_min), xbmax=csv.as_double(r, cols.c_xb_max);
            if (std::fabs(xbmin-xb.first)<1e-9 && std::fabs(xbmax-xb.second)<1e-9) {
                xbmeans.push_back((c_xbavg>=0)? csv.as_double(r,c_xbavg) : 0.5*(xb.first+xb.second));
            }
        }
        const double xb_mean_for_title = safe_mean(xbmeans);

        for (int rr=0; rr<nrows; ++rr) for (int cc=0; cc<ncols; ++cc) {
            const auto& qpair = Q2s[rr]; const auto& tpair = Ts[cc];
            std::vector<int> rows_for_cell;
            for (int r=0;r<csv.nrows();++r) {
                const double xbmin=csv.as_double(r, cols.c_xb_min), xbmax=csv.as_double(r, cols.c_xb_max);
                const double q2min=csv.as_double(r, cols.c_q2_min), q2max=csv.as_double(r, cols.c_q2_max);
                const double tmin =csv.as_double(r, cols.c_tab_min), tmax =csv.as_double(r, cols.c_tab_max);
                if (std::fabs(xbmin-xb.first)<1e-9 && std::fabs(xbmax-xb.second)<1e-9 &&
                    std::fabs(q2min-qpair.first)<1e-9 && std::fabs(q2max-qpair.second)<1e-9 &&
                    std::fabs(tmin -tpair.first)<1e-9 && std::fabs(tmax -tpair.second)<1e-9)
                    rows_for_cell.push_back(r);
            }
            std::sort(rows_for_cell.begin(), rows_for_cell.end(), [&](int a,int b){
                return csv.as_double(a, cols.c_phi_min) < csv.as_double(b, cols.c_phi_min);
            });

            CellData& C = cells[rr*ncols+cc];
            C.X.reserve(rows_for_cell.size()); C.EX.assign(rows_for_cell.size(),0.0);

            for (int r : rows_for_cell) {
                const double pmin=csv.as_double(r, cols.c_phi_min), pmax=csv.as_double(r, cols.c_phi_max);
                double xphi=0.5*(pmin+pmax);
                if (c_phiavg>=0) {
                    const double pav=csv.as_double(r, c_phiavg);
                    if (std::isfinite(pav) && pav>0.0 && pav<360.0) xphi=pav;
                }
                C.X.push_back(xphi);

                const std::string ybase = std::string("raw yield, ep->epg, ")+topo_str+", exp, "+label+", ";
                const int c_pos = csv.has_col(ybase+"pos")? csv.col_index(ybase+"pos") : -1;
                const int c_neg = csv.has_col(ybase+"neg")? csv.col_index(ybase+"neg") : -1;

                double yp= (c_pos>=0)? csv.as_double(r, c_pos) : std::numeric_limits<double>::quiet_NaN();
                double yn= (c_neg>=0)? csv.as_double(r, c_neg) : std::numeric_limits<double>::quiet_NaN();

                C.Yp.push_back(yp); C.Ym.push_back(yn);
                C.EYp.push_back(std::isfinite(yp)? std::sqrt(std::max(0.0,yp)) : 0.0);
                C.EYm.push_back(std::isfinite(yn)? std::sqrt(std::max(0.0,yn)) : 0.0);

                if (c_q2avg>=0)  C.q2means.push_back(csv.as_double(r, c_q2avg));
                else             C.q2means.push_back(0.5*(qpair.first+qpair.second));
                if (c_tabavg>=0) C.tmeans.push_back(csv.as_double(r, c_tabavg));
                else             C.tmeans.push_back(0.5*(tpair.first+tpair.second));

                if (std::isfinite(yp)) { canvas_max=std::max(canvas_max, yp); if (yp>=500.0) need_log=true; }
                if (std::isfinite(yn)) { canvas_max=std::max(canvas_max, yn); if (yn>=500.0) need_log=true; }
            }
        }

        std::string cname="c_counts_"+period_dir_for_label(label)+"_"+topo_dir+"_xB_"+std::to_string((int)std::round(xb.first*1000.0));
        TCanvas* c=new TCanvas(cname.c_str(), "", W, H);

        TPad* pTop = new TPad("pTop","pTop",0.0,0.84,1.0,1.0);
        TPad* pGrid= new TPad("pGrid","pGrid",0.0,0.00,1.0,0.84);
        pTop->SetFillStyle(0); pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols,nrows,0.0001,0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(head_size); head.SetTextAlign(22);
        head.DrawLatex(0.5,0.50, Form("%s   <xB>=%.3g   %s", label.c_str(), xb_mean_for_title, topo_str.c_str()));

        const double y_lo = need_log ? 0.5 : 0.0;
        const double y_hi = std::max(1.0, canvas_max*1.15);

        for (int rr=0; rr<nrows; ++rr) for (int cc=0; cc<ncols; ++cc) {
            pGrid->cd(rr*ncols+cc+1);
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.205);   // tiny bit more left padding
            gPad->SetRightMargin(0.07);
            if (need_log) gPad->SetLogy(1); else gPad->SetLogy(0);

            TH1* frame=gPad->DrawFrame(0.0, y_lo, 360.0, y_hi);
            frame->GetXaxis()->SetTitle("phi (deg)");
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(lab_size);
            frame->GetYaxis()->SetTitleSize(lab_size);
            frame->GetXaxis()->SetLabelSize(tick_size);
            frame->GetYaxis()->SetLabelSize(tick_size);

            const CellData& C = cells[rr*ncols+cc];

            TGraphErrors* gp=new TGraphErrors((int)C.X.size(), (double*)C.X.data(), (double*)C.Yp.data(), (double*)C.EX.data(), (double*)C.EYp.data());
            TGraphErrors* gm=new TGraphErrors((int)C.X.size(), (double*)C.X.data(), (double*)C.Ym.data(), (double*)C.EX.data(), (double*)C.EYm.data());
            gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);  gp->SetLineColor(kRed);  gp->SetLineWidth(1);
            gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue); gm->SetLineWidth(1);
            gp->Draw("PE1 SAME"); gm->Draw("PE1 SAME");

            TLatex lab; lab.SetNDC(); lab.SetTextSize(latex_size); lab.SetTextAlign(13);
            const double q2m=safe_mean(C.q2means), tm=safe_mean(C.tmeans);
            lab.DrawLatex(0.12, 0.83, Form("Q^{2}=%.3g  |t|=%.3g", q2m, tm));

            // Legend bottom-left moved a bit left for more room
            TLegend* leg=new TLegend(0.58, 0.74, 0.95, 0.93);
            leg->SetBorderSize(1);
            leg->SetLineColor(kBlack);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetTextSize(latex_size);
            leg->AddEntry(gp,"helicity +","pe");
            leg->AddEntry(gm,"helicity -","pe");
            leg->Draw();
        }

        const std::string fpath=(fs::path(out_dir)/("plot_total_counts_"+period_dir_for_label(label)+"_"+topo_dir+"_xB_"+std::to_string((int)std::round(xb.first*1000.0))+".png")).string();
        c->SaveAs(fpath.c_str());
        delete c;
    }
}

// -------------- selection helpers --------------
static inline bool passBaseCuts(const std::vector<BaseCut>& BC,
                                const std::unordered_map<std::string,BranchBinding>& B) {
    for (const auto& c:BC) {
        auto it=B.find(c.branch); if (it==B.end()) return false;
        const BranchBinding& bb=it->second; const double v=bb_as_double(bb);
        if (c.has_eq)  { if (bb_as_ll(bb)!=c.eq) return false; }
        if (c.has_neq) { if (bb_as_ll(bb)==c.neq) return false; }
        if (c.has_min) { if (!(v>=c.vmin)) return false; }
        if (c.has_max) { if (!(v<=c.vmax)) return false; }
    }
    return true;
}
static inline bool passSigmaCuts(const std::vector<SigmaCut>& SC,
                                 const std::unordered_map<std::string,BranchBinding>& B) {
    for (const auto& c:SC) {
        auto it=B.find(c.branch); if (it==B.end()) return false;
        const double v=bb_as_double(it->second);
        const double lo=c.center - c.nsigma*c.sigma, hi=c.center + c.nsigma*c.sigma;
        if (c.mode==SigmaCut::TWO_SIDED) { if (v<lo || v>hi) return false; }
        else if (c.mode==SigmaCut::UPPER) { if (v>hi) return false; }
        else { if (v<lo) return false; }
    }
    return true;
}
static inline bool passesAllCuts(const PeriodCuts& cuts,
                                 const std::unordered_map<std::string,BranchBinding>& B) {
    return passBaseCuts(cuts.base_cuts,B) && passSigmaCuts(cuts.sigma_cuts,B);
}

// A small struct to understand why events fail to write
struct DebugCounts {
    long long total=0, bad_helicity=0, bad_nan=0, cut_fail=0, topo_bad=0, no_row_match=0;
    long long kin_ok_phi_fail=0, phi_ok_kin_fail=0, both_ok=0;
};

static const char* TOPO_DIR(int idx){ return TOPO_DIRS[idx]; }
static const char* TOPO_STR(int idx){ return TOPO_STRS[idx]; }

struct RowCounts { long long pos=0, neg=0; };

} // end anonymous namespace

// ---------------- EXPORTED FUNCTION (global scope) ----------------
// Note the leading '::' to ensure global scope definition, fixing the link error.
bool update_total_counts_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_root_dir,
    int max_workers)
{
    namespace fs=std::filesystem;
    ROOT::EnableThreadSafety();

    // CSV existence + pre-stats
    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    uintmax_t size_before = fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;
    auto mtime_before = fs::exists(csv_path, ec) ? fs::last_write_time(csv_path, ec) : fs::file_time_type{};

    std::cout << "[total_counts] CSV: " << csv_abs << "  (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) return false;

    CsvCols cols;
    cols.c_xb_min  = csv.col_index("xBmin");
    cols.c_xb_max  = csv.col_index("xBmax");
    cols.c_q2_min  = csv.col_index("Q2min");
    cols.c_q2_max  = csv.col_index("Q2max");
    cols.c_tab_min = csv.col_index("t_abs_min");
    cols.c_tab_max = csv.col_index("t_abs_max");
    cols.c_phi_min = csv.col_index("phimin");
    cols.c_phi_max = csv.col_index("phimax");
    if (cols.c_xb_min<0 || cols.c_xb_max<0 || cols.c_q2_min<0 || cols.c_q2_max<0 ||
        cols.c_tab_min<0 || cols.c_tab_max<0 || cols.c_phi_min<0 || cols.c_phi_max<0) {
        fatal("Missing one or more required bin-edge columns in CSV.");
    }

    // periods list (canonical only)
    std::vector<std::string> periods;
    for (const auto& P:CANONICAL_PERIODS()) periods.push_back(P.tree_key);

    // Helper to ensure yield columns exist (for non-supplemental periods only)
    auto ensure_yield_columns = [&](const std::string& label, const std::string& topo_str) {
        const bool write_ok = (label!="Fa18 Inb Supp");
        if (!write_ok) return std::tuple<int,int,int>(-1,-1,-1);
        const std::string base = std::string("raw yield, ep->epg, ")+topo_str+", exp, "+label+", ";
        const int c_unpol = csv.ensure_col(base+"unpol");
        const int c_pos   = csv.ensure_col(base+"pos");
        const int c_neg   = csv.ensure_col(base+"neg");
        return std::make_tuple(c_unpol,c_pos,c_neg);
    };

    using PeriodRowMap = std::unordered_map<int, RowCounts>;
    std::mutex merge_mtx;
    std::map<std::string, std::map<std::string, PeriodRowMap>> counts_by_label_topo;

    auto cadence_for = [](Long64_t N)->Long64_t {
        Long64_t by_pct = (Long64_t)std::max(1.0, std::floor(0.02 * (double)N));
        Long64_t by_abs = 1000000;
        return std::min(by_pct, by_abs);
    };

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(max_workers)
#endif
    for (int ip=0; ip<(int)periods.size(); ++ip) {
        const std::string period_key = periods[ip];
        if (!is_canonical_tree_key(period_key)) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            std::cerr << "[total_counts] Non-canonical period: " << period_key << "\n";
            continue;
        }
        auto itT = dvcsDataTrees.find(period_key);
        if (itT==dvcsDataTrees.end() || !itT->second) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            std::cerr << "[total_counts] Missing TTree for " << period_key << "\n";
            continue;
        }
        TTree* t = itT->second;
        const std::string label = safe_label_for_key(period_key);

        const PeriodCuts cuts = loadCombinedCuts(combined_cuts_json, period_key);

        // Build needed branch list
        std::vector<std::string> need;
        for (const auto& c:cuts.base_cuts)  need.push_back(c.branch);
        for (const auto& c:cuts.sigma_cuts) need.push_back(c.branch);
        const char* explicit_binds[] = {"helicity","x","Q2","t1","phi2","detector1","detector2"};
        for (const char* s:explicit_binds) need.push_back(s);
        std::sort(need.begin(), need.end());
        need.erase(std::unique(need.begin(), need.end()), need.end());

        // Branch status gating
        t->SetBranchStatus("*", 0);
        for (const auto& b:need) if (t->GetBranch(b.c_str())) t->SetBranchStatus(b.c_str(), 1);

        if (list_enabled_branches()) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            {
                std::cout << "[total_counts]["<<label<<"] Enabled branches:\n";
                for (const auto& b:need) {
                    Bool_t st = t->GetBranchStatus(b.c_str());
                    std::cout << "  - " << b << " : " << (st? "ON":"OFF") << "\n";
                }
            }
        }

        // Typed binds
        int helicity=0; double x=0, Q2=0, t1=0, phi2=std::numeric_limits<double>::quiet_NaN();
        if (!t->GetBranch("helicity") || !t->GetBranch("x") || !t->GetBranch("Q2") ||
            !t->GetBranch("t1") || !t->GetBranch("phi2")) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            std::cerr << "[total_counts] Required branches missing in " << period_key << "\n";
            continue;
        }
        t->SetBranchAddress("helicity",&helicity);
        t->SetBranchAddress("x",&x);
        t->SetBranchAddress("Q2",&Q2);
        t->SetBranchAddress("t1",&t1);
        t->SetBranchAddress("phi2",&phi2);

        TopologyResolver topo; topo.bind(t);

        // Generic binds for cut variables
        std::vector<std::string> need_for_binder;
        for (const auto& b:need) {
            bool drop=false; for (const char* s:explicit_binds) if (b==s) { drop=true; break; }
            if (!drop) need_for_binder.push_back(b);
        }
        std::unordered_map<std::string,BranchBinding> bind;
        bindRequiredBranches_STRICT(t, need_for_binder, bind);

        std::map<std::string, PeriodRowMap> local_by_topo;

        const Long64_t N=t->GetEntries();
        const Long64_t cadence=cadence_for(N);
        Long64_t seen=0, kept=0, used=0;
        DebugCounts dbg;

#ifdef _OPENMP
        #pragma omp critical
#endif
        std::cout << "[total_counts] Start " << label << " with " << (long long)N << " entries\n";

        for (Long64_t i=0;i<N;++i) {
            t->GetEntry(i);
            dbg.total++;

            if (helicity!=+1 && helicity!=-1) { dbg.bad_helicity++; continue; }
            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi2)) { dbg.bad_nan++; continue; }
            ++seen;

            // Cuts
            if (!passesAllCuts(cuts, bind)) { dbg.cut_fail++; continue; }
            ++kept;

            const int topo_idx=topo.index();
            if (topo_idx<0 || topo_idx>2) { dbg.topo_bad++; continue; }
            const std::string topo_str = TOPO_STR(topo_idx);

            const double phi_deg = wrap_deg_0_360(rad_to_deg(phi2));
            bool used_any=false;

            // For diagnostics: separate kin-only and phi-only tests
            bool matched_kin_any=false, matched_phi_any=false;

            for (int r=0;r<csv.nrows();++r) {
                const double xbmin=csv.as_double(r, cols.c_xb_min);
                const double xbmax=csv.as_double(r, cols.c_xb_max);
                const double q2min=csv.as_double(r, cols.c_q2_min);
                const double q2max=csv.as_double(r, cols.c_q2_max);
                const double tabmin=csv.as_double(r, cols.c_tab_min);
                const double tabmax=csv.as_double(r, cols.c_tab_max);
                const double pmin=csv.as_double(r, cols.c_phi_min);
                const double pmax=csv.as_double(r, cols.c_phi_max);

                const bool kin_ok = in_range(x, xbmin, xbmax) &&
                                    in_range(Q2, q2min, q2max) &&
                                    in_range(std::fabs(t1), tabmin, tabmax);
                const bool phi_ok = row_accepts_phi(phi_deg, pmin, pmax);

                matched_kin_any |= kin_ok;
                matched_phi_any |= phi_ok;

                if (kin_ok && phi_ok) {
                    RowCounts& rc = local_by_topo[topo_str][r];
                    if (helicity==+1) rc.pos++; else rc.neg++;
                    used_any=true;

                    if (trace_matches()) {
#ifdef _OPENMP
                        #pragma omp critical
#endif
                        std::cout << "[total_counts]["<<label<<"] match: row="<<r
                                  << " x="<<x<<" Q2="<<Q2<<" |t|="<<std::fabs(t1)
                                  << " phi="<<phi_deg<<"\n";
                }
                }
            }

            if (matched_kin_any && !matched_phi_any) dbg.kin_ok_phi_fail++;
            if (!matched_kin_any && matched_phi_any) dbg.phi_ok_kin_fail++;
            if (matched_kin_any && matched_phi_any)  dbg.both_ok++;
            if (!used_any) dbg.no_row_match++;

            if (used_any) ++used;

            if ((i % cadence)==0 && is_debug()) {
                double pct = (N>0)? 100.0*(double)i/(double)N : 100.0;
                std::cout << "[total_counts]["<<label<<"] " << std::fixed << std::setprecision(1)
                          << pct << "%  i="<<(long long)i<<" seen="<<seen<<" kept="<<kept<<" used="<<used<<"\n";
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            std::cout << "[total_counts] Done " << label << "  seen="<<seen<<" kept="<<kept<<" used="<<used<<"\n";
            if (is_debug()) {
                std::cout << "[total_counts]["<<label<<"] debug: total="<<dbg.total
                          << " bad_helicity="<<dbg.bad_helicity
                          << " bad_nan="<<dbg.bad_nan
                          << " cut_fail="<<dbg.cut_fail
                          << " topo_bad="<<dbg.topo_bad
                          << " no_row_match="<<dbg.no_row_match
                          << " kin_ok_phi_fail="<<dbg.kin_ok_phi_fail
                          << " phi_ok_kin_fail="<<dbg.phi_ok_kin_fail
                          << " both_ok="<<dbg.both_ok
                          << "\n";
            }
        }

        // Merge
        {
            std::lock_guard<std::mutex> lock(merge_mtx);
            for (const auto& tk: local_by_topo) {
                const std::string topo_str=tk.first;
                (void)ensure_yield_columns(label, topo_str); // ensure headers (except supplemental)
                auto& tgt = counts_by_label_topo[label][topo_str];
                for (const auto& rkv : tk.second) {
                    tgt[rkv.first].pos += rkv.second.pos;
                    tgt[rkv.first].neg += rkv.second.neg;
                }
            }
        }
    } // omp periods

    // Write period counts to CSV
    size_t grand_cells_written=0;
    for (const auto& lblkv : counts_by_label_topo) {
        const std::string& label = lblkv.first;
        const bool write_ok = (label!="Fa18 Inb Supp");
        for (int topo_idx=0; topo_idx<3; ++topo_idx) {
            const std::string topo_str = TOPO_STR(topo_idx);
            const auto itTopo = lblkv.second.find(topo_str);
            if (itTopo==lblkv.second.end()) {
                std::cout << "[total_counts] NOTE: no rows for " << label << " / " << topo_str << "\n";
                continue;
            }
            const std::string base = std::string("raw yield, ep->epg, ")+topo_str+", exp, "+label+", ";
            int c_unpol = csv.col_index(base+"unpol");
            int c_pos   = csv.col_index(base+"pos");
            int c_neg   = csv.col_index(base+"neg");

            if (!write_ok) {
                std::cout << "[total_counts] skip CSV write for supplemental: " << label << "\n";
                continue;
            }
            if (c_unpol<0 || c_pos<0 || c_neg<0) {
                std::cout << "[total_counts] HEADER MISSING for " << label << " / " << topo_str
                          << " (expected \"" << base << "unpol\" etc.)\n";
                continue;
            }

            size_t cells_written_here=0;
            bool printed_example=false;
            for (const auto& rkv : itTopo->second) {
                const int r = rkv.first;
                const long long pos=rkv.second.pos, neg=rkv.second.neg, unp=pos+neg;
                csv.set(r, c_unpol, (double)unp); ++cells_written_here;
                csv.set(r, c_pos,   (double)pos); ++cells_written_here;
                csv.set(r, c_neg,   (double)neg); ++cells_written_here;

                if (!printed_example && is_debug()) {
                    std::cout << "[total_counts] example write: row="<<r<<"  "
                              << "unpol="<<unp<<" pos="<<pos<<" neg="<<neg
                              << "  ("<<label<<" / "<<topo_str<<")\n";
                    printed_example=true;
                }
            }
            grand_cells_written += cells_written_here;
            std::cout << "[total_counts] wrote " << cells_written_here
                      << " cells for " << label << " / " << topo_str << "\n";
        }
    }

    // Combined groups
    std::map<std::string, std::vector<std::string>> group_members = {
        {"Fa18",    {"Fa18 Inb","Fa18 Out"}},
        {"Sp18",    {"Sp18 Inb","Sp18 Out"}},
        {"10.6 GeV",{"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"}} // adjust if you want Sp19 aggregated too
    };

    for (const auto& gkv : group_members) {
        const std::string& group = gkv.first;
        const auto& members = gkv.second;
        for (int topo_idx=0; topo_idx<3; ++topo_idx) {
            const std::string topo_str=TOPO_STR(topo_idx);
            const std::string gbase = std::string("raw yield, ep->epg, ")+topo_str+", exp, "+group+", ";
            const int c_unpol = csv.ensure_col(gbase+"unpol");
            const int c_pos   = csv.ensure_col(gbase+"pos");
            const int c_neg   = csv.ensure_col(gbase+"neg");

            size_t group_cells_written=0;
            for (int r=0;r<csv.nrows();++r) {
                long long pos_sum=0, neg_sum=0;
                for (const auto& lbl : members) {
                    if (lbl=="Fa18 Inb Supp") continue;
                    const std::string pbase = std::string("raw yield, ep->epg, ")+topo_str+", exp, "+lbl+", ";
                    const int p_pos = csv.col_index(pbase+"pos");
                    const int p_neg = csv.col_index(pbase+"neg");
                    if (p_pos>=0) pos_sum += (long long)std::llround(csv.as_double(r, p_pos));
                    if (p_neg>=0) neg_sum += (long long)std::llround(csv.as_double(r, p_neg));
                }
                csv.set(r, c_pos,   (double)pos_sum); ++group_cells_written;
                csv.set(r, c_neg,   (double)neg_sum); ++group_cells_written;
                csv.set(r, c_unpol, (double)(pos_sum+neg_sum)); ++group_cells_written;
            }
            grand_cells_written += group_cells_written;
            std::cout << "[total_counts] wrote " << group_cells_written
                      << " cells for combined group " << group << " / " << topo_str << "\n";
        }
    }

    if (grand_cells_written==0) {
        std::cout << "[total_counts] WARNING: wrote zero CSV cells. This usually means no rows matched "
                     "the kinematic/phi windows or header names didn’t match.\n";
    }

    // Save CSV
    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to save " << csv_path << "\n";
        return false;
    }

    // Post-save stats
    uintmax_t size_after = fs::file_size(csv_path, ec);
    auto mtime_after = fs::last_write_time(csv_path, ec);
    std::cout << "[total_counts] Updated raw yields in: " << csv_abs
              << "  (cells written: " << grand_cells_written
              << ", size " << size_before << " -> " << size_after << ")\n";

    // Plots (canonical dirs only)
    for (const auto& P:CANONICAL_PERIODS()) {
        const std::string label = P.label;
        for (int topo_idx=0; topo_idx<3; ++topo_idx)
            draw_group_canvases(label, TOPO_STR(topo_idx), TOPO_DIR(topo_idx), csv, cols, out_root_dir);
    }
    for (const std::string group : {"Fa18","Sp18","10.6 GeV"}) {
        for (int topo_idx=0; topo_idx<3; ++topo_idx)
            draw_group_canvases(group, TOPO_STR(topo_idx), TOPO_DIR(topo_idx), csv, cols, out_root_dir);
    }

    return true;
}