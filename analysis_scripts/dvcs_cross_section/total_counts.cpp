// total_counts.cpp — strict raw-yield writer + plotting + strict group sums
// - Reads DVCS TTrees, applies global + exclusivity cuts (bin_means parity)
// - Bins surviving events into CSV rows by (xB, Q2, |t|, phi)
// - Writes ONLY the canonical raw-yield columns produced by initialize_pass2_csv.cpp
// - Per-period raw yields: required columns, no auto-creation
// - Combined-group raw yields (Fa18, Sp18, 10.6 GeV): required columns, sums of member periods
// - Plots per-period and combined groups
//
// Robustness fix: CSV label resolver maps keys like "fa18_inb" to header's "Fa18 Inb"

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

namespace {

static inline bool is_debug() { return std::getenv("TOTAL_COUNTS_DEBUG") != nullptr; }
static inline bool trace_matches() { return std::getenv("TOTAL_COUNTS_TRACE_MATCHES") != nullptr; }
static inline bool trace_labels() { return std::getenv("TOTAL_COUNTS_TRACE_LABELS") != nullptr; }

[[noreturn]] void fatal(const std::string& msg) {
    std::cerr << "[total_counts] FATAL: " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline double PI() { return 3.14159265358979323846; }
static inline double RAD2DEG(double r) { return r * 180.0 / PI(); }

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
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg); // wrap-around
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
            for (size_t i=0;i<header.size();++i) { write_field(fout, header[i]); if (i+1<header.size()) fout<<','; }
            fout << "\n";
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

// ---------------- discover and resolve CSV labels ----------------
static std::string trim_copy(const std::string& s) {
    size_t a=0, b=s.size();
    while (a<b && std::isspace((unsigned char)s[a])) ++a;
    while (b>a && std::isspace((unsigned char)s[b-1])) --b;
    return s.substr(a,b-a);
}
static std::string norm_key(std::string s) {
    std::string out; out.reserve(s.size());
    for (char c : s) {
        if (c==' ' || c=='_' || c=='\t') continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }
    return out;
}

// Builds an index: normalized label -> canonical label as it appears in the CSV header.
// It scans header fields like: "raw yield, ep->epg, (FD, FD), exp, Fa18 Inb, pos" etc.
static std::unordered_map<std::string,std::string> build_label_index(const CsvDoc& csv) {
    std::unordered_map<std::string,std::string> idx;
    for (const auto& h : csv.header) {
        // Quick filter to speed up
        if (h.find("raw yield, ep->epg") == std::string::npos) continue;
        size_t exp_pos = h.find(", exp, ");
        if (exp_pos == std::string::npos) continue;
        size_t last_comma = h.rfind(',');
        if (last_comma == std::string::npos || last_comma <= exp_pos+7) continue;
        std::string label = trim_copy(h.substr(exp_pos + 7, last_comma - (exp_pos + 7)));
        if (label.empty()) continue;
        std::string k = norm_key(label);
        // only first wins; all headers for pos/neg/unpol will map the same label
        if (!k.empty() && idx.find(k)==idx.end()) idx[k] = label;
    }
    if (trace_labels()) {
        std::cerr << "[total_counts] Discovered CSV labels:\n";
        for (auto& kv : idx) std::cerr << "  key=" << kv.first << " -> \"" << kv.second << "\"\n";
    }
    return idx;
}

// Given a desired label (possibly "fa18_inb"), find the exact CSV label ("Fa18 Inb") if present.
static std::string resolve_csv_label(const std::unordered_map<std::string,std::string>& idx,
                                     const std::string& want) {
    auto it = idx.find(norm_key(want));
    return (it==idx.end()) ? want : it->second;
}

// ---------------- period label helpers (for display/dirs only) ----------------
static inline std::string period_dir_for_label(const std::string& L) {
    if (L=="Fa18 Inb")      return "Fa18_Inb";
    if (L=="Fa18 Out")      return "Fa18_Out";
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

static bool is_group_label(const std::string& L) {
    return (L=="Fa18"||L=="Sp18"||L=="10.6 GeV");
}

static std::string out_root_for_label(const std::string& label, const std::string& out_root_dir) {
    namespace fs=std::filesystem;
    const std::string base = (fs::path(out_root_dir) / "total_counts_plots").string();
    if (is_group_label(label)) return (fs::path(base) / label).string();
    return (fs::path(base) / period_dir_for_label(label)).string();
}

// ---------------- topology ----------------
enum class Topology { FD_FD=0, CD_FD=1, CD_FT=2 };

static inline const char* topo_str(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "(FD, FD)";
        case Topology::CD_FD: return "(CD, FD)";
        case Topology::CD_FT: return "(CD, FT)";
    }
    return "(FD, FD)";
}
static inline const char* topo_dir(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

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
        if (detector1==1 && detector2==1) return (int)Topology::FD_FD;
        if (detector1==2 && detector2==1) return (int)Topology::CD_FD;
        if (detector1==2 && detector2==0) return (int)Topology::CD_FT;
        return -1;
    }
};

// ---------------- plotting helpers ----------------
struct CellData { std::vector<double> X,Yp,Ym,EX,EYp,EYm,q2means,tmeans; };
static double safe_mean(const std::vector<double>& v){
    double s=0; int n=0; for(double x:v) if (std::isfinite(x)){ s+=x; ++n; }
    return n? s/n : std::numeric_limits<double>::quiet_NaN();
}

// ---------------- branch binding (type-aware) ----------------
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

static std::mutex g_root_bind_mutex;

static void bind_one_exact(TTree* t, const std::string& bname, BranchBinding& bb) {
    std::lock_guard<std::mutex> lock(g_root_bind_mutex);
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

// ---------------- simple exclusivity policy (bin_means parity) ----------------
// Hard global cuts:
static inline bool passes_global(double open_angle_ep2_deg, double t1, double pTmiss) {
    if (!(open_angle_ep2_deg > 5.0)) return false;
    if (!((-t1) < 1.0)) return false;       // |t| < 1.0
    if (!(pTmiss <= 0.20)) return false;    // (GeV)
    return true;
}

// ---------------- diagnostics and row counts ----------------
struct DebugCounts {
    long long total=0, bad_helicity=0, bad_nan=0, cut_fail=0, topo_bad=0, no_row_match=0;
    long long kin_ok_phi_fail=0, phi_ok_kin_fail=0, both_ok=0;
};
struct RowCounts { long long pos=0, neg=0; };

struct YieldCols { int unpol=-1, pos=-1, neg=-1; };

static YieldCols require_raw_yield_cols(const CsvDoc& csv,
                                        const std::string& csv_label,
                                        const std::string& topo_str_label) {
    const std::string base = std::string("raw yield, ep->epg, ")+topo_str_label+", exp, "+csv_label+", ";
    YieldCols yc;
    yc.unpol = csv.col_index(base+"unpol");
    yc.pos   = csv.col_index(base+"pos");
    yc.neg   = csv.col_index(base+"neg");
    if (yc.unpol<0 || yc.pos<0 || yc.neg<0) {
        fatal(std::string("Missing canonical raw-yield columns for \"")+csv_label+"\" / \""+topo_str_label+"\" "
              "expected \""+base+"unpol|pos|neg\"");
    }
    return yc;
}

struct CombinedGroup { std::string name; std::vector<std::string> members; };
static const std::vector<CombinedGroup>& combined_groups() {
    static const std::vector<CombinedGroup> G = {
        {"Fa18",     {"Fa18 Inb","Fa18 Out"}},
        {"Sp18",     {"Sp18 Inb","Sp18 Out"}},
        {"10.6 GeV", {"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"}}
    };
    return G;
}
static YieldCols require_group_raw_yield_cols(const CsvDoc& csv,
                                              const std::string& csv_group_label,
                                              const std::string& topo_str_label) {
    const std::string base = std::string("raw yield, ep->epg, ")+topo_str_label+", exp, "+csv_group_label+", ";
    YieldCols yc;
    yc.unpol = csv.col_index(base+"unpol");
    yc.pos   = csv.col_index(base+"pos");
    yc.neg   = csv.col_index(base+"neg");
    if (yc.unpol<0 || yc.pos<0 || yc.neg<0) {
        fatal(std::string("Missing canonical group raw-yield columns for \"")+csv_group_label+"\" / \""+topo_str_label+"\" "
              "expected \""+base+"unpol|pos|neg\"");
    }
    return yc;
}

// ---------------- plotting ----------------
struct CsvCols { int c_xb_min, c_xb_max, c_q2_min, c_q2_max, c_tab_min, c_tab_max, c_phi_min, c_phi_max; };

static void draw_group_canvases(
    const std::string& display_label,             // for titles and directories
    const std::string& topo_str_label,
    const std::string& topo_dir_label,
    const CsvDoc& csv,
    const CsvCols& cols,
    const std::unordered_map<std::string,std::string>& label_idx,
    const std::string& out_root_dir)
{
    namespace fs=std::filesystem;
    const std::string base_dir = out_root_for_label(display_label, out_root_dir);
    ensure_dir(base_dir);
    const std::string out_dir = (fs::path(base_dir) / topo_dir_label).string();
    ensure_dir(out_dir);

    // Resolve the actual label as it appears in the CSV header
    const std::string csv_label = resolve_csv_label(label_idx, display_label);
    if (trace_labels() && csv_label != display_label) {
        std::cerr << "[total_counts] plot label remap: \"" << display_label << "\" -> \"" << csv_label << "\"\n";
    }

    const int c_phiavg = csv.has_col("phiavg, "+csv_label)? csv.col_index("phiavg, "+csv_label) : -1;
    const int c_q2avg  = csv.has_col("Q2avg, " +csv_label)? csv.col_index("Q2avg, " +csv_label) : -1;
    const int c_tabavg = csv.has_col("t_abs_avg, "+csv_label)? csv.col_index("t_abs_avg, "+csv_label) : -1;
    const int c_xbavg  = csv.has_col("xBavg, "+csv_label)? csv.col_index("xBavg, "+csv_label) : -1;

    std::set<std::pair<double,double>> xb_set;
    for (int r=0;r<csv.nrows();++r) xb_set.emplace(csv.as_double(r,cols.c_xb_min), csv.as_double(r,cols.c_xb_max));

    const double head_size = 0.58;
    const double latex_size= 0.065;

    for (auto xb : xb_set) {
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

        std::vector<CellData> cells(nrows*ncols);
        double canvas_max = 1.0;

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

                const std::string ybase = std::string("raw yield, ep->epg, ")+topo_str_label+", exp, "+csv_label+", ";
                const int c_pos = csv.has_col(ybase+"pos")? csv.col_index(ybase+"pos") : -1;
                const int c_neg = csv.has_col(ybase+"neg")? csv.col_index(ybase+"neg") : -1;

                double yp= (c_pos>=0)? csv.as_double(r, c_pos) : std::numeric_limits<double>::quiet_NaN();
                double yn= (c_neg>=0)? csv.as_double(r, c_neg) : std::numeric_limits<double>::quiet_NaN();

                C.Yp.push_back(yp); C.Ym.push_back(yn);
                C.EYp.push_back(std::isfinite(yp)? std::sqrt(std::max(0.0,yp)) : 0.0);
                C.EYm.push_back(std::isfinite(yn)? std::sqrt(std::max(0.0,yn)) : 0.0);

                C.q2means.push_back(0.5*(qpair.first+qpair.second));
                C.tmeans.push_back(0.5*(tpair.first+tpair.second));

                if (std::isfinite(yp)) canvas_max=std::max(canvas_max, yp);
                if (std::isfinite(yn)) canvas_max=std::max(canvas_max, yn);
            }
        }

        std::string cname="c_counts_"+period_dir_for_label(display_label)+"_"+topo_dir_label+"_xB_"+std::to_string((int)std::round(xb.first*1000.0));
        TCanvas* c=new TCanvas(cname.c_str(), "", W, H);

        TPad* pTop = new TPad("pTop","pTop",0.0,0.84,1.0,1.0);
        TPad* pGrid= new TPad("pGrid","pGrid",0.0,0.00,1.0,0.84);
        pTop->SetFillStyle(0); pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols,nrows,0.0001,0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(head_size); head.SetTextAlign(22);
        head.DrawLatex(0.5,0.50, Form("%s   <xB>=%.3g   %s", display_label.c_str(), xb_mean_for_title, topo_str_label.c_str()));

        const double y_lo = 0.0;
        const double y_hi = std::max(1.0, canvas_max*1.15);

        for (int rr=0; rr<nrows; ++rr) for (int cc=0; cc<ncols; ++cc) {
            pGrid->cd(rr*ncols+cc+1);
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.125);
            gPad->SetRightMargin(0.07);

            TH1* frame=gPad->DrawFrame(0.0, y_lo, 360.0, y_hi);
            frame->GetXaxis()->SetTitle("phi (deg)");
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.075);
            frame->GetYaxis()->SetTitleSize(0.075);
            frame->GetXaxis()->SetLabelSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.060);

            const CellData& C = cells[rr*ncols+cc];

            TGraphErrors* gp=new TGraphErrors((int)C.X.size(), (double*)C.X.data(), (double*)C.Yp.data(), (double*)C.EX.data(), (double*)C.EYp.data());
            TGraphErrors* gm=new TGraphErrors((int)C.X.size(), (double*)C.X.data(), (double*)C.Ym.data(), (double*)C.EX.data(), (double*)C.EYm.data());
            gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);  gp->SetLineColor(kRed);  gp->SetLineWidth(1);
            gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue); gm->SetLineWidth(1);
            gp->Draw("PE1 SAME"); gm->Draw("PE1 SAME");

            TLatex lab; lab.SetNDC(); lab.SetTextSize(latex_size); lab.SetTextAlign(13);
            const double q2m=safe_mean(C.q2means), tm=safe_mean(C.tmeans);
            lab.DrawLatex(0.12, 0.83, Form("Q^{2}=%.3g  |t|=%.3g", q2m, tm));

            TLegend* leg=new TLegend(0.58, 0.74, 0.95, 0.93);
            leg->SetBorderSize(1);
            leg->SetLineColor(kBlack);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetTextSize(0.065);
            leg->AddEntry(gp,"helicity +","pe");
            leg->AddEntry(gm,"helicity -","pe");
            leg->Draw();
        }

        const std::string fpath=(fs::path(out_dir)/("plot_total_counts_"+period_dir_for_label(display_label)+"_"+topo_dir_label+"_xB_"+std::to_string((int)std::round(xb.first*1000.0))+".png")).string();
        c->SaveAs(fpath.c_str());
        delete c;
    }
}

} // end anon ns

// ---------------- EXPORTED FUNCTION ----------------
bool update_total_counts_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,   // still accepted; global hard-cuts applied irrespective
    const std::string& out_root_dir,
    int max_workers)
{
    namespace fs=std::filesystem;
    ROOT::EnableThreadSafety();

    const std::string csv_abs = fs::absolute(csv_path).string();
    std::error_code ec;
    uintmax_t size_before = fs::exists(csv_path, ec) ? fs::file_size(csv_path, ec) : 0;

    std::cout << "[total_counts] CSV: " << csv_abs << "  (size=" << size_before << ")\n";

    CsvDoc csv;
    if (!csv.load(csv_path)) return false;

    // Build once: header label index for robust lookups
    const auto label_idx = build_label_index(csv);

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

    // Canonical period keys from periods.h
    std::vector<std::string> period_keys;
    for (const auto& P : CANONICAL_PERIODS()) {
        if (dvcsDataTrees.count(P.tree_key) && dvcsDataTrees.at(P.tree_key))
            period_keys.push_back(P.tree_key);
    }
    if (period_keys.empty()) fatal("No DVCS trees available in dvcsDataTrees.");

    using PeriodRowMap = std::unordered_map<int, RowCounts>;
    std::mutex merge_mtx;
    std::map<std::string, std::map<std::string, PeriodRowMap>> counts_by_label_topo; // keyed by *display* label

    auto cadence_for = [](Long64_t N)->Long64_t {
        Long64_t by_pct = (Long64_t)std::max(1.0, std::floor(0.02 * (double)N));
        Long64_t by_abs = 1000000;
        return std::min(by_pct, by_abs);
    };

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic) num_threads(max_workers)
#endif
    for (int ip=0; ip<(int)period_keys.size(); ++ip) {
        const std::string period_key = period_keys[ip];

        auto itT = dvcsDataTrees.find(period_key);
        if (itT==dvcsDataTrees.end() || !itT->second) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            std::cerr << "[total_counts] Missing TTree for " << period_key << "\n";
            continue;
        }
        TTree* t = itT->second;

        // From periods.h: human-readable label (this is the one you expect like "Fa18 Inb")
        std::string display_label;
        for (const auto& P : CANONICAL_PERIODS()) {
            if (P.tree_key == period_key) { display_label = P.label; break; }
        }
        if (display_label.empty()) {
#ifdef _OPENMP
            #pragma omp critical
#endif
            std::cerr << "[total_counts] Unknown period key naming: " << period_key << "\n";
            continue;
        }

        // Bind needed branches
        BranchBinding b_helicity, b_x, b_Q2, b_t1, b_phi2, b_open_angle, b_pTmiss;
        bind_one_exact(t, "helicity",       b_helicity);
        bind_one_exact(t, "x",              b_x);
        bind_one_exact(t, "Q2",             b_Q2);
        bind_one_exact(t, "t1",             b_t1);
        bind_one_exact(t, "phi2",           b_phi2);           // radians
        bind_one_exact(t, "open_angle_ep2", b_open_angle);     // degrees
        bind_one_exact(t, "pTmiss",         b_pTmiss);         // GeV

        TopologyResolver topo; topo.bind(t);

        const Long64_t N=t->GetEntries();
        const Long64_t cadence=cadence_for(N);
        Long64_t seen=0, kept=0, used=0;
        DebugCounts dbg;

#ifdef _OPENMP
        #pragma omp critical
#endif
        std::cout << "[total_counts] Start " << display_label << " with " << (long long)N << " entries\n";

        std::map<std::string, PeriodRowMap> local_by_topo;

        for (Long64_t i=0;i<N;++i) {
            t->GetEntry(i);
            dbg.total++;

            const long long hel = bb_as_ll(b_helicity);
            if (hel!=+1 && hel!=-1) { dbg.bad_helicity++; continue; }

            const double x   = bb_as_double(b_x);
            const double Q2  = bb_as_double(b_Q2);
            const double t1  = bb_as_double(b_t1);
            const double phi2= bb_as_double(b_phi2);
            const double open_angle_deg = bb_as_double(b_open_angle);
            const double pTmiss         = bb_as_double(b_pTmiss);

            if (!std::isfinite(x) || !std::isfinite(Q2) || !std::isfinite(t1) || !std::isfinite(phi2)
                || !std::isfinite(open_angle_deg) || !std::isfinite(pTmiss)) {
                dbg.bad_nan++; continue;
            }
            ++seen;

            // Global hard cuts (bin_means parity)
            if (!passes_global(open_angle_deg, t1, pTmiss)) { dbg.cut_fail++; continue; }
            ++kept;

            const int topo_idx = topo.index();
            if (topo_idx < 0 || topo_idx > 2) { dbg.topo_bad++; continue; }
            const std::string topoS = topo_str((Topology)topo_idx);

            const double phi_deg = wrap_deg_0_360(RAD2DEG(phi2));
            bool used_any=false;
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
                    RowCounts& rc = local_by_topo[topoS][r];
                    if (hel==+1) rc.pos++; else rc.neg++;
                    used_any=true;

                    if (trace_matches()) {
#ifdef _OPENMP
                        #pragma omp critical
#endif
                        std::cout << "[total_counts]["<<display_label<<"] match: row="<<r
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
                std::cout << "[total_counts]["<<display_label<<"] " << std::fixed << std::setprecision(1)
                          << pct << "%  i="<<(long long)i<<" seen="<<seen<<" kept="<<kept<<" used="<<used<<"\n";
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            std::cout << "[total_counts] Done " << display_label << "  seen="<<seen<<" kept="<<kept<<" used="<<used<<"\n";
            if (is_debug()) {
                std::cout << "[total_counts]["<<display_label<<"] debug: total="<<dbg.total
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

        // Merge results
        {
            std::lock_guard<std::mutex> lock(merge_mtx);
            auto& tgt_by_topo = counts_by_label_topo[display_label]; // keep display label as the key
            for (const auto& tk : local_by_topo) {
                const std::string& topoS = tk.first;
                auto& tgt = tgt_by_topo[topoS];
                for (const auto& rkv : tk.second) {
                    tgt[rkv.first].pos += rkv.second.pos;
                    tgt[rkv.first].neg += rkv.second.neg;
                }
            }
        }
    } // omp periods

    // Write per-period raw yields (resolve labels against header)
    size_t grand_cells_written=0;
    for (const auto& lblkv : counts_by_label_topo) {
        const std::string display_label = lblkv.first;
        const std::string csv_label = resolve_csv_label(label_idx, display_label);
        if (trace_labels() && csv_label != display_label) {
            std::cerr << "[total_counts] write label remap: \"" << display_label << "\" -> \"" << csv_label << "\"\n";
        }
        for (int ti=0; ti<3; ++ti) {
            const std::string topoS = topo_str((Topology)ti);
            auto itTopo = lblkv.second.find(topoS);
            if (itTopo==lblkv.second.end()) continue;

            const YieldCols yc = require_raw_yield_cols(csv, csv_label, topoS);

            size_t cells_written_here=0;
            for (const auto& rkv : itTopo->second) {
                const int r = rkv.first;
                const long long pos=rkv.second.pos, neg=rkv.second.neg, unp=pos+neg;
                csv.set(r, yc.unpol, (double)unp); ++cells_written_here;
                csv.set(r, yc.pos,   (double)pos); ++cells_written_here;
                csv.set(r, yc.neg,   (double)neg); ++cells_written_here;
            }
            grand_cells_written += cells_written_here;
            std::cout << "[total_counts] wrote " << cells_written_here
                      << " cells for " << display_label << " / " << topoS << "\n";
        }
    }

    // Strict combined-group raw yields
    for (const auto& G : combined_groups()) {
        const std::string display_group = G.name;
        const std::string csv_group = resolve_csv_label(label_idx, display_group); // usually identical
        for (int ti=0; ti<3; ++ti) {
            const std::string topoS = topo_str((Topology)ti);
            const YieldCols gy = require_group_raw_yield_cols(csv, csv_group, topoS);

            size_t cells_written_here=0;
            for (int r=0; r<csv.nrows(); ++r) {
                long long pos_sum=0, neg_sum=0;
                for (const auto& member_display : G.members) {
                    const std::string member_csv = resolve_csv_label(label_idx, member_display);
                    const std::string pbase = std::string("raw yield, ep->epg, ")+topoS+", exp, "+member_csv+", ";
                    const int p_pos = csv.col_index(pbase+"pos");
                    const int p_neg = csv.col_index(pbase+"neg");
                    if (p_pos>=0) pos_sum += (long long)std::llround(csv.as_double(r, p_pos));
                    if (p_neg>=0) neg_sum += (long long)std::llround(csv.as_double(r, p_neg));
                }
                csv.set(r, gy.pos,   (double)pos_sum); ++cells_written_here;
                csv.set(r, gy.neg,   (double)neg_sum); ++cells_written_here;
                csv.set(r, gy.unpol, (double)(pos_sum+neg_sum)); ++cells_written_here;
            }

            grand_cells_written += cells_written_here;
            std::cout << "[total_counts] wrote " << cells_written_here
                      << " cells for combined group " << display_group << " / " << topoS << "\n";
        }
    }

    if (grand_cells_written==0) {
        std::cout << "[total_counts] WARNING: wrote zero CSV cells (no matches or headers mismatch).\n";
    }

    // Save CSV
    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to save " << csv_path << "\n";
        return false;
    }

    // Plots (per-period + combined groups); resolve labels for header lookups internally
    for (const auto& P:CANONICAL_PERIODS()) {
        const std::string display_label = P.label;
        for (int ti=0; ti<3; ++ti)
            draw_group_canvases(display_label, topo_str((Topology)ti), topo_dir((Topology)ti), csv, cols, label_idx, out_root_dir);
    }
    for (const std::string display_group : {"Fa18","Sp18","10.6 GeV"}) {
        for (int ti=0; ti<3; ++ti)
            draw_group_canvases(display_group, topo_str((Topology)ti), topo_dir((Topology)ti), csv, cols, label_idx, out_root_dir);
    }

    uintmax_t size_after = fs::file_size(csv_path, ec);
    std::cout << "[total_counts] Updated raw yields in: " << csv_abs
              << "  (cells written: " << grand_cells_written
              << ", size " << size_before << " -> " << size_after << ")\n";

    return true;
}