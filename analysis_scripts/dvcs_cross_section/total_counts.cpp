// total_counts.cpp
// Raw-yield writer + plots with legends, correct axis labels, margins, and
// canonicalized output directories. Combined groups are plotted from in-memory
// aggregated period counts with CSV fallback to ensure non-zero curves.

#include <nlohmann/json.hpp>

#include "total_counts.h"
#include "periods.h"
#include "load_binning_scheme.h"
#include "label_aliases.h"

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
#include <cstdlib>
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

// phi-bin acceptance checker (handles wrap-around)
static bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (!std::isfinite(phi_deg) || !std::isfinite(pmin_deg) || !std::isfinite(pmax_deg)) return false;
    if (pmax_deg > pmin_deg) return in_range(phi_deg, pmin_deg, pmax_deg);
    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg);
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
        for (char ch : s) {
            if (ch=='"') os << "\"\"";
            else os << ch;
        }
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

// ---------------- small utils ----------------
using json = nlohmann::json;

struct SigmaCut {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std  = std::numeric_limits<double>::quiet_NaN();
};

using VarCutMap   = std::unordered_map<std::string, SigmaCut>;       // var -> {mean,std}
using TopoCutMap  = std::unordered_map<std::string, VarCutMap>;      // "DVCS_<PeriodDir>_<TopoDir>" -> VarCutMap

static TopoCutMap load_sigma_cuts_data(const std::string& path) {
    TopoCutMap out;
    if (path.empty()) return out;

    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::cerr << "[total_counts] WARNING: could not open combined cuts JSON: " << path << "\n";
        return out;
    }
    json j; fin >> j;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key(); // e.g. "DVCS_Fa18_Inb_CD_FD"
        const json& block = it.value();

        // We are counting DATA here; use the "data" object if present.
        if (!block.contains("data") || !block["data"].is_object()) continue;

        VarCutMap m;
        for (auto vit = block["data"].begin(); vit != block["data"].end(); ++vit) {
            const std::string vname = vit.key(); // Emiss2, Mx2, ...
            const json& vs = vit.value();
            SigmaCut sc;
            if (vs.contains("mean") && vs["mean"].is_number()) sc.mean = vs["mean"].get<double>();
            if (vs.contains("std")  && vs["std"].is_number())  sc.std  = vs["std"].get<double>();
            if (std::isfinite(sc.std) && sc.std > 0.0) {
                m[vname] = sc;
            }
        }
        if (!m.empty()) out[key] = std::move(m);
    }
    std::cout << "[total_counts] Loaded sigma cuts for " << out.size() << " topology keys from " << path << "\n";
    return out;
}

static inline bool within_3sigma(double val, const SigmaCut& sc) {
    if (!std::isfinite(val) || !std::isfinite(sc.mean) || !std::isfinite(sc.std) || sc.std <= 0.0) return true;
    return std::fabs(val - sc.mean) <= 3.0 * sc.std;
}

static std::string to_lower_nospace(std::string s) {
    std::string out; out.reserve(s.size());
    for (char c : s) {
        if (c==' ' || c=='\t' || c=='_') continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }
    return out;
}

// Skip CSV write/lookup for:
//  - supplemental period(s)
//  - combined groups (Fa18, Sp18, 10.6 GeV)
static bool should_skip_csv_for_label(const std::string& display_label) {
    const std::string k = to_lower_nospace(display_label);
    if (k == "fa18inbsupp" || k == "fa18_inb_supp" || k.find("supp") != std::string::npos) return true;
    if (k == "fa18" || k == "sp18" || k == "10.6gev" || k == "10.6") return true;
    return false;
}

static bool is_group_label(const std::string& L) {
    return (L=="Fa18"||L=="Sp18"||L=="10.6 GeV");
}

static std::vector<std::string> group_members(const std::string& group_label) {
    if (group_label=="Fa18")     return {"Fa18 Inb","Fa18 Out"};
    if (group_label=="Sp18")     return {"Sp18 Inb","Sp18 Out"};
    if (group_label=="10.6 GeV") return {"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"};
    return {};
}

// ------------- Alias-aware header resolution -----------
static int find_col_alias(const std::vector<std::string>& header,
                          const std::string& topo_str_label,
                          const std::string& period_display_label,
                          const std::string& helicity_label) {
    for (const auto& topo_try : topology_aliases(topo_str_label)) {
        for (const auto& per_try : period_aliases(period_display_label)) {
            for (const auto& hel_try : helicity_aliases(helicity_label)) {
                const std::string candidate = "raw yield, ep->epg, " + topo_try + ", exp, " + per_try + ", " + hel_try;
                auto it = std::find(header.begin(), header.end(), candidate);
                if (it != header.end()) return int(it - header.begin());
            }
        }
    }
    return -1;
}

struct YieldCols { int unpol=-1, pos=-1, neg=-1; };

static bool try_raw_yield_cols_alias(const CsvDoc& csv,
                                     const std::string& display_period_label,
                                     const std::string& topo_str_label,
                                     YieldCols& out) {
    out.unpol = find_col_alias(csv.header, topo_str_label, display_period_label, "unpol");
    out.pos   = find_col_alias(csv.header, topo_str_label, display_period_label, "pos");
    out.neg   = find_col_alias(csv.header, topo_str_label, display_period_label, "neg");
    return (out.unpol>=0 && out.pos>=0 && out.neg>=0);
}

// ---------------- directory helpers ----------------
static std::string canonical_period_dir(const std::string& L) {
    const std::string k = to_lower_nospace(L);
    if (k=="fa18inb") return "Fa18_Inb";
    if (k=="fa18out") return "Fa18_Out";
    if (k=="fa18inbsupp") return "Fa18_Inb_Supp";
    if (k=="sp18inb") return "Sp18_Inb";
    if (k=="sp18out") return "Sp18_Out";
    if (k=="sp19inb") return "Sp19_Inb";
    std::string s=L; for (char& c:s) if (c==' ') c='_';
    return s;
}

static std::string canonical_group_dir(const std::string& L) {
    if (L=="Fa18") return "Fa18";
    if (L=="Sp18") return "Sp18";
    if (L=="10.6 GeV") return "10.6 GeV";
    return L;
}

static inline std::string period_dir_for_label(const std::string& L) {
    return canonical_period_dir(L);
}

static std::string out_root_for_label(const std::string& label, const std::string& out_root_dir) {
    namespace fs=std::filesystem;
    const std::string base = (fs::path(out_root_dir) / "total_counts_plots").string();
    if (is_group_label(label)) return (fs::path(base) / canonical_group_dir(label)).string();
    return (fs::path(base) / canonical_period_dir(label)).string();
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

// Merge/move any alias-named subdirs into canonical names to stop duplicates.
static void canonicalize_output_dirs(const std::string& total_counts_root) {
    namespace fs=std::filesystem;
    std::error_code ec;

    auto merge_into = [&](const fs::path& src, const fs::path& dst){
        if (!fs::exists(src)) return;
        fs::create_directories(dst, ec);
        for (auto& p : fs::directory_iterator(src)) {
            fs::path to = dst / p.path().filename();
            fs::rename(p.path(), to, ec);
            if (ec) {
                // fallback copy+remove if rename across filesystems
                ec.clear();
                fs::copy(p.path(), to, fs::copy_options::overwrite_existing | fs::copy_options::recursive, ec);
                fs::remove_all(p.path(), ec);
            }
        }
        fs::remove_all(src, ec);
    };

    if (!fs::exists(total_counts_root)) return;

    for (auto& p : fs::directory_iterator(total_counts_root)) {
        if (!p.is_directory()) continue;
        const std::string name = p.path().filename().string();
        // recognize either a period-like or group-like name and rewrite
        std::string canon = canonical_period_dir(name);
        if (is_group_label(name)) canon = canonical_group_dir(name);
        if (canon != name) {
            merge_into(p.path(), fs::path(total_counts_root)/canon);
        }
    }
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
    void enable_and_bind(TTree* t) {
        t->SetBranchStatus("*", 0);
        t->SetBranchStatus("detector1", 1);
        t->SetBranchStatus("detector2", 1);
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

static std::mutex g_root_bind_mutex;

static void bind_one_exact_enable(TTree* t, const std::string& bname, BranchBinding& bb) {
    std::lock_guard<std::mutex> lock(g_root_bind_mutex);
    t->SetBranchStatus(bname.c_str(), 1);
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

// ---------------- simple exclusivity policy ----------------
static inline bool passes_global(double open_angle_ep2_deg, double t1, double pTmiss) {
    if (!(open_angle_ep2_deg > 5.0)) return false;
    if (!((-t1) < 1.0)) return false;
    if (!(pTmiss <= 0.20)) return false;
    return true;
}

// ---------------- diagnostics and row counts ----------------
struct DebugCounts {
    long long total=0, bad_helicity=0, bad_nan=0, cut_fail=0, topo_bad=0, no_row_match=0;
    long long kin_ok_phi_fail=0, phi_ok_kin_fail=0, both_ok=0;
    long long bad_io=0;
};
struct RowCounts { long long pos=0, neg=0; };

// ---------------- Bin-edge pre-cache ----------------
struct BinCache {
    std::vector<double> xbmin, xbmax, q2min, q2max, tmin, tmax, pmin, pmax;
    int nrows() const { return (int)xbmin.size(); }
};

static BinCache precache_bins(const CsvDoc& csv) {
    auto get = [&](const char* name)->int {
        int idx = csv.col_index(name);
        if (idx<0) fatal(std::string("Missing bin-edge column: ")+name);
        return idx;
    };
    const int cxmin = get("xBmin");
    const int cxmax = get("xBmax");
    const int cqmin = get("Q2min");
    const int cqmax = get("Q2max");
    const int ctmin = get("t_abs_min");
    const int ctmax = get("t_abs_max");
    const int cpmin = get("phimin");
    const int cpmax = get("phimax");

    BinCache B;
    const int NR = csv.nrows();
    B.xbmin.resize(NR); B.xbmax.resize(NR);
    B.q2min.resize(NR); B.q2max.resize(NR);
    B.tmin.resize(NR);  B.tmax.resize(NR);
    B.pmin.resize(NR);  B.pmax.resize(NR);
    for (int r=0;r<NR;++r) {
        B.xbmin[r] = csv.as_double(r, cxmin);
        B.xbmax[r] = csv.as_double(r, cxmax);
        B.q2min[r] = csv.as_double(r, cqmin);
        B.q2max[r] = csv.as_double(r, cqmax);
        B.tmin[r]  = csv.as_double(r, ctmin);
        B.tmax[r]  = csv.as_double(r, ctmax);
        B.pmin[r]  = csv.as_double(r, cpmin);
        B.pmax[r]  = csv.as_double(r, cpmax);
    }
    return B;
}

// ---------------- plotting ----------------
struct CsvCols { int c_xb_min, c_xb_max, c_q2_min, c_q2_max, c_tab_min, c_tab_max, c_phi_min, c_phi_max; };

// Sum from CSV across group members (fallback if aggregated map lacks the row)
static void sum_csv_for_group_row(const CsvDoc& csv,
                                  int row_index,
                                  const std::string& topo_str_label,
                                  const std::vector<std::string>& members,
                                  double& sum_pos, double& sum_neg)
{
    sum_pos = 0.0; sum_neg = 0.0;
    bool any=false;
    for (const auto& m : members) {
        int cpos = find_col_alias(csv.header, topo_str_label, m, "pos");
        int cneg = find_col_alias(csv.header, topo_str_label, m, "neg");
        if (cpos>=0 && cneg>=0) {
            double vp = csv.as_double(row_index, cpos);
            double vn = csv.as_double(row_index, cneg);
            if (std::isfinite(vp)) { sum_pos += vp; any=true; }
            if (std::isfinite(vn)) { sum_neg += vn; any=true; }
        }
    }
    if (!any) { sum_pos = 0.0; sum_neg = 0.0; }
}

static void draw_group_canvases(
    const std::string& display_label,
    const std::string& topo_str_label,
    const std::string& topo_dir_label,
    const CsvDoc& csv,
    const CsvCols& cols,
    const std::string& out_root_dir,
    const std::map<std::string, std::unordered_map<int, RowCounts>>* counts_opt
)
{
    namespace fs=std::filesystem;

    // Clean up any alias dirs once per process (cheap no-op if already clean)
    static std::once_flag canon_once;
    std::call_once(canon_once, [&](){
        canonicalize_output_dirs((fs::path(out_root_dir)/"total_counts_plots").string());
    });

    const std::string base_dir = out_root_for_label(display_label, out_root_dir);
    ensure_dir(base_dir);
    const std::string out_dir = (fs::path(base_dir) / topo_dir_label).string();
    ensure_dir(out_dir);

    int c_phiavg=-1, c_q2avg=-1, c_tabavg=-1, c_xbavg=-1;
    {
        c_phiavg = csv.col_index("phiavg, "+display_label);
        c_q2avg  = csv.col_index("Q2avg, " +display_label);
        c_tabavg = csv.col_index("t_abs_avg, "+display_label);
        c_xbavg  = csv.col_index("xBavg, "+display_label);
        if (c_phiavg<0 || c_q2avg<0 || c_tabavg<0 || c_xbavg<0) {
            for (const auto& alt : period_aliases(display_label)) {
                if (c_phiavg<0) { int ci = csv.col_index("phiavg, "+alt); if (ci>=0) c_phiavg=ci; }
                if (c_q2avg <0) { int ci = csv.col_index("Q2avg, " +alt); if (ci>=0) c_q2avg =ci; }
                if (c_tabavg<0) { int ci = csv.col_index("t_abs_avg, "+alt); if (ci>=0) c_tabavg=ci; }
                if (c_xbavg <0) { int ci = csv.col_index("xBavg, " +alt); if (ci>=0) c_xbavg =ci; }
            }
        }
    }

    std::set<std::pair<double,double>> xb_set;
    for (int r=0;r<csv.nrows();++r) xb_set.emplace(csv.as_double(r,cols.c_xb_min), csv.as_double(r,cols.c_xb_max));

    // visual tuning (smaller titles, restore legend box, keep left padding)
    const double head_size = 0.14;   // was 0.28
    const double label_sz  = 0.050;
    const double title_sz  = 0.045;
    const double leg_txt   = 0.050;

    const bool is_group = is_group_label(display_label);
    const bool force_counts_only = should_skip_csv_for_label(display_label);
    const std::vector<std::string> members = is_group ? group_members(display_label) : std::vector<std::string>{};

    // one-time summary of combined totals if applicable
    static std::once_flag once_flag;
    std::call_once(once_flag, [&](){
        if (is_group && counts_opt) {
            long long totp=0, totn=0;
            for (const auto& tk : *counts_opt) {
                for (const auto& rkv : tk.second) { totp += rkv.second.pos; totn += rkv.second.neg; }
            }
            std::cout << "[total_counts] Combined \""<<display_label<<"\" totals: pos="<<totp<<" neg="<<totn<<"\n";
        }
    });

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

            const int c_pos_csv = force_counts_only ? -1 : find_col_alias(csv.header, topo_str_label, display_label, "pos");
            const int c_neg_csv = force_counts_only ? -1 : find_col_alias(csv.header, topo_str_label, display_label, "neg");

            for (int r : rows_for_cell) {
                const double pmin=csv.as_double(r, cols.c_phi_min), pmax=csv.as_double(r, cols.c_phi_max);
                double xphi=0.5*(pmin+pmax);
                if (c_phiavg>=0) {
                    const double pav=csv.as_double(r, c_phiavg);
                    if (std::isfinite(pav) && pav>0.0 && pav<360.0) xphi=pav;
                }
                C.X.push_back(xphi);

                double yp = std::numeric_limits<double>::quiet_NaN();
                double yn = std::numeric_limits<double>::quiet_NaN();

                if (c_pos_csv>=0 && c_neg_csv>=0) {
                    yp = csv.as_double(r, c_pos_csv);
                    yn = csv.as_double(r, c_neg_csv);
                } else if (counts_opt) {
                    auto topo_it = counts_opt->find(topo_str_label);
                    if (topo_it != counts_opt->end()) {
                        auto rc_it = topo_it->second.find(r);
                        if (rc_it != topo_it->second.end()) {
                            yp = (double)rc_it->second.pos;
                            yn = (double)rc_it->second.neg;
                        } else if (is_group && !members.empty()) {
                            // Fallback: sum CSV across the group's member periods for this row
                            double sp=0.0, sn=0.0;
                            sum_csv_for_group_row(csv, r, topo_str_label, members, sp, sn);
                            yp = sp; yn = sn;
                        } else {
                            yp = 0.0; yn = 0.0;
                        }
                    } else if (is_group && !members.empty()) {
                        double sp=0.0, sn=0.0;
                        sum_csv_for_group_row(csv, r, topo_str_label, members, sp, sn);
                        yp = sp; yn = sn;
                    } else {
                        yp = 0.0; yn = 0.0;
                    }
                } else if (is_group && !members.empty()) {
                    double sp=0.0, sn=0.0;
                    sum_csv_for_group_row(csv, r, topo_str_label, members, sp, sn);
                    yp = sp; yn = sn;
                } else {
                    yp = 0.0; yn = 0.0;
                }

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

        // smaller top bar
        TPad* pTop = new TPad("pTop","pTop",0.0,0.86,1.0,1.0);
        TPad* pGrid= new TPad("pGrid","pGrid",0.0,0.00,1.0,0.86);
        pTop->SetFillStyle(0); pTop->Draw();
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols,nrows,0.0001,0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(head_size); head.SetTextAlign(22);
        head.DrawLatex(0.5,0.58, Form("%s   <xB>=%.3g   %s", display_label.c_str(), xb_mean_for_title, topo_str_label.c_str()));

        const double y_lo = 0.0;
        const double y_hi = std::max(1.0, canvas_max*1.15);

        for (int rr=0; rr<nrows; ++rr) for (int cc=0; cc<ncols; ++cc) {
            pGrid->cd(rr*ncols+cc+1);
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.24);     // a touch more room to prevent clipping at top
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.125);   // requested left padding
            gPad->SetRightMargin(0.07);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            TH1* frame=gPad->DrawFrame(0.0, y_lo, 360.0, y_hi);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Counts");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(title_sz);
            frame->GetYaxis()->SetTitleSize(title_sz);
            frame->GetXaxis()->SetLabelSize(label_sz);
            frame->GetYaxis()->SetLabelSize(label_sz);
            frame->GetYaxis()->SetTitleOffset(1.25);   // keep y-title clear of the left edge
            frame->GetXaxis()->SetTitleOffset(1.05);

            const CellData& C = cells[rr*ncols+cc];

            TGraphErrors* gp=new TGraphErrors((int)C.X.size(), (double*)C.X.data(), (double*)C.Yp.data(), (double*)C.EX.data(), (double*)C.EYp.data());
            TGraphErrors* gm=new TGraphErrors((int)C.X.size(), (double*)C.X.data(), (double*)C.Ym.data(), (double*)C.EX.data(), (double*)C.EYm.data());
            gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);  gp->SetLineColor(kRed);  gp->SetLineWidth(1);
            gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue); gm->SetLineWidth(1);
            gp->Draw("PE1 SAME"); gm->Draw("PE1 SAME");

            // per-pad legend, with border box restored
            TLegend* leg = new TLegend(0.60, 0.72, 0.93, 0.92);
            leg->SetBorderSize(1);    
            leg->SetLineColor(kBlack);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetTextSize(leg_txt);
            leg->AddEntry(gp, "+ helicity", "PE");
            leg->AddEntry(gm, "- helicity", "PE");
            leg->Draw();

            // per-pad kinematics label, nudged clear of the legend and top edge
            TLatex lab; lab.SetNDC(); lab.SetTextSize(0.058); lab.SetTextAlign(13);
            const double q2m=safe_mean(C.q2means), tm=safe_mean(C.tmeans);
            lab.DrawLatex(0.16, 0.88, Form("Q^{2}=%.3g  |t|=%.3g", q2m, tm));
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
    const std::string& combined_cuts_json,
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

    BinCache bins = precache_bins(csv);

    CsvCols cols;
    cols.c_xb_min  = csv.col_index("xBmin");
    cols.c_xb_max  = csv.col_index("xBmax");
    cols.c_q2_min  = csv.col_index("Q2min");
    cols.c_q2_max  = csv.col_index("Q2max");
    cols.c_tab_min = csv.col_index("t_abs_min");
    cols.c_tab_max = csv.col_index("t_abs_max");
    cols.c_phi_min = csv.col_index("phimin");
    cols.c_phi_max = csv.col_index("phimax");

    std::vector<std::string> period_keys;
    for (const auto& P : CANONICAL_PERIODS()) {
        if (dvcsDataTrees.count(P.tree_key) && dvcsDataTrees.at(P.tree_key))
            period_keys.push_back(P.tree_key);
    }
    if (period_keys.empty()) fatal("No DVCS trees available in dvcsDataTrees.");

    using PeriodRowMap = std::unordered_map<int, RowCounts>;
    std::mutex merge_mtx;
    std::map<std::string, std::map<std::string, PeriodRowMap>> counts_by_label_topo;

    auto cadence_for = [](Long64_t N)->Long64_t {
        Long64_t by_pct = (Long64_t)std::max(1.0, std::floor(0.02 * (double)N));
        Long64_t by_abs = 1000000;
        return std::min(by_pct, by_abs);
    };

    int threads = max_workers;
#ifdef _OPENMP
    if (threads <= 0) {
        threads = omp_get_max_threads();
    }
    if (threads < 1) threads = 1;
    #pragma omp parallel for schedule(dynamic) num_threads(threads)
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

        TopologyResolver topo; topo.enable_and_bind(t);

        BranchBinding b_helicity, b_x, b_Q2, b_t1, b_phi2, b_open_angle, b_pTmiss;
        bind_one_exact_enable(t, "helicity",       b_helicity);
        bind_one_exact_enable(t, "x",              b_x);
        bind_one_exact_enable(t, "Q2",             b_Q2);
        bind_one_exact_enable(t, "t1",             b_t1);
        bind_one_exact_enable(t, "phi2",           b_phi2);
        bind_one_exact_enable(t, "open_angle_ep2", b_open_angle);
        bind_one_exact_enable(t, "pTmiss",         b_pTmiss);
        BranchBinding b_Emiss2, b_Mx2, b_Mx2_1, b_Mx2_2, b_theta_gg, b_xF;
        bind_one_exact_enable(t, "Emiss2",            b_Emiss2);
        bind_one_exact_enable(t, "Mx2",               b_Mx2);
        bind_one_exact_enable(t, "Mx2_1",             b_Mx2_1);
        bind_one_exact_enable(t, "Mx2_2",             b_Mx2_2);
        bind_one_exact_enable(t, "theta_gamma_gamma", b_theta_gg);
        bind_one_exact_enable(t, "xF",                b_xF);

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
            const Long64_t nb = t->GetEntry(i);
            if (nb < 0) { dbg.bad_io++; continue; }

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

            if (!passes_global(open_angle_deg, t1, pTmiss)) { dbg.cut_fail++; continue; }
            const int topo_idx = topo.index();
            if (topo_idx < 0 || topo_idx > 2) { dbg.topo_bad++; continue; }
            const std::string topoS = topo_str((Topology)topo_idx);

            // Compose key exactly like the JSON uses: "DVCS_<PeriodDir>_<TopoDir>"
            const std::string topo_key = std::string("DVCS_") + period_dir_for_label(display_label) + "_" + topo_dir((Topology)topo_idx);
            auto itTopoCuts = sigmaCuts.find(topo_key);
            if (itTopoCuts != sigmaCuts.end()) {
                // Only check vars that appear in the JSON block; missing vars are ignored.
                auto getv = [&](const std::string& v)->double {
                    if (v=="Emiss2")            return bb_as_double(b_Emiss2);
                    if (v=="Mx2")               return bb_as_double(b_Mx2);
                    if (v=="Mx2_1")             return bb_as_double(b_Mx2_1);
                    if (v=="Mx2_2")             return bb_as_double(b_Mx2_2);
                    if (v=="theta_gamma_gamma") return bb_as_double(b_theta_gg);
                    if (v=="xF")                return bb_as_double(b_xF);
                    if (v=="pTmiss")            return pTmiss; // also enforce 3σ on pTmiss if provided
                    return std::numeric_limits<double>::quiet_NaN();
                };

                bool ok = true;
                for (const auto& kv : itTopoCuts->second) {
                    const std::string& vname = kv.first;
                    const SigmaCut&    sc    = kv.second;
                    const double val = getv(vname);
                    if (!within_3sigma(val, sc)) { ok = false; break; }
                }
                if (!ok) continue; // fail the 3σ exclusivity for this topology
            }
            ++kept;

            const int topo_idx = topo.index();
            if (topo_idx < 0 || topo_idx > 2) { dbg.topo_bad++; continue; }
            const std::string topoS = topo_str((Topology)topo_idx);

            const double phi_deg = wrap_deg_0_360(RAD2DEG(phi2));
            bool used_any=false;
            bool matched_kin_any=false, matched_phi_any=false;

            const int NR = bins.nrows();
            for (int r=0;r<NR;++r) {
                const bool kin_ok = in_range(x,  bins.xbmin[r], bins.xbmax[r]) &&
                                    in_range(Q2, bins.q2min[r], bins.q2max[r]) &&
                                    in_range(std::fabs(t1), bins.tmin[r], bins.tmax[r]);
                const bool phi_ok = row_accepts_phi(phi_deg, bins.pmin[r], bins.pmax[r]);

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
                          << " bad_io="<<dbg.bad_io
                          << "\n";
            }
        }

        {
            std::lock_guard<std::mutex> lock(merge_mtx);
            auto& tgt_by_topo = counts_by_label_topo[display_label];
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

    // Write per-period raw yields if present and not in skip list
    size_t grand_cells_written=0;
    for (const auto& lblkv : counts_by_label_topo) {
        const std::string display_label = lblkv.first;

        if (should_skip_csv_for_label(display_label)) {
            std::cout << "[total_counts] Skipping CSV write for \"" << display_label << "\" (supplemental or combined-group).\n";
            continue;
        }

        for (int ti=0; ti<3; ++ti) {
            const std::string topoS = topo_str((Topology)ti);
            auto itTopo = lblkv.second.find(topoS);
            if (itTopo==lblkv.second.end()) continue;

            YieldCols yc;
            if (!try_raw_yield_cols_alias(csv, display_label, topoS, yc)) {
                std::cout << "[total_counts] NOTE: no CSV headers for \"" << display_label
                          << "\" / \"" << topoS << "\"; skipping write.\n";
                continue;
            }

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

    if (grand_cells_written==0) {
        std::cout << "[total_counts] WARNING: wrote zero CSV cells (headers missing or skip-list only).\n";
    }

    if (!csv.save_atomic(csv_path)) {
        std::cerr << "[total_counts] ERROR: failed to save " << csv_path << "\n";
        return false;
    }

    // -------- Plotting --------
    // (A) Per-period plots
    for (const auto& P:CANONICAL_PERIODS()) {
        const std::string display_label = P.label;
        const auto it_counts = counts_by_label_topo.find(display_label);
        const std::map<std::string, std::unordered_map<int, RowCounts>>* counts_ptr =
            (it_counts!=counts_by_label_topo.end() ? &it_counts->second : nullptr);
        for (int ti=0; ti<3; ++ti) {
            draw_group_canvases(display_label,
                                topo_str((Topology)ti),
                                topo_dir((Topology)ti),
                                csv, cols, out_root_dir,
                                counts_ptr);
        }
    }

    // (B) Combined-group plots: aggregate from period maps
    struct CombinedGroup { std::string name; std::vector<std::string> members; };
    static const std::vector<CombinedGroup> G = {
        {"Fa18",     {"Fa18 Inb","Fa18 Out"}},
        {"Sp18",     {"Sp18 Inb","Sp18 Out"}},
        {"10.6 GeV", {"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"}}
    };

    for (const auto& GG : G) {
        std::map<std::string, std::unordered_map<int, RowCounts>> agg;
        for (const std::string topoS : {std::string("(FD, FD)"), std::string("(CD, FD)"), std::string("(CD, FT)")}) {
            auto& topo_map = agg[topoS];
            for (const std::string& member_display : GG.members) {
                auto it_mem = counts_by_label_topo.find(member_display);
                if (it_mem == counts_by_label_topo.end()) continue;
                auto it_topo = it_mem->second.find(topoS);
                if (it_topo == it_mem->second.end()) continue;
                for (const auto& rkv : it_topo->second) {
                    RowCounts& dst = topo_map[rkv.first];
                    dst.pos += rkv.second.pos;
                    dst.neg += rkv.second.neg;
                }
            }
        }
        for (int ti=0; ti<3; ++ti) {
            draw_group_canvases(GG.name,
                                topo_str((Topology)ti),
                                topo_dir((Topology)ti),
                                csv, cols, out_root_dir,
                                &agg);
        }
    }

    uintmax_t size_after = fs::file_size(csv_path, ec);
    std::cout << "[total_counts] Updated raw yields in: " << csv_abs
              << "  (cells written: " << grand_cells_written
              << ", size " << size_before << " -> " << size_after << ")\n";

    return true;
}