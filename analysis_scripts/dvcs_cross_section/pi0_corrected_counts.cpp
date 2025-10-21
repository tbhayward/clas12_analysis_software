// π0-corrected counts with errors, for periods AND combined groups.
// Uses total_counts.json ("groups": {...}) and:
//  - per-period contamination_<period>.json (each has "bins": {...})
//  - combined pi0_contamination_combined.json ("periods": {Group: {bins:{...}}})
//
// Outputs one JSON per group + a master JSON, and error-bar plots
// (xB in canvas title; Q² and −t in subplot headers).

#include "pi0_corrected_counts.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TString.h>
#include <TColor.h>

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
#include <utility>
#include <vector>

namespace {

constexpr int    N_PHI_BINS = 12;
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)

// ────────── styles ──────────
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetTitleFont(42,"XYZ");
        gStyle->SetLabelFont(42,"XYZ");
        gStyle->SetTextFont(42);
    }
} _style_bootstrap;

// ────────── helpers ──────────
static inline std::string toLower(std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){return std::tolower(c);});
    return s;
}
static inline std::string titleCaseTail(const std::string& tail) {
    // e.g. "Fa18_inb" -> "Fa18_Inb"
    std::string t = tail;
    if (!t.empty()) t[0] = std::toupper(static_cast<unsigned char>(t[0]));
    for (size_t i = 0; i + 1 < t.size(); ++i) {
        if (t[i] == '_' && std::isalpha(static_cast<unsigned char>(t[i+1]))) {
            t[i+1] = std::toupper(static_cast<unsigned char>(t[i+1]));
        }
    }
    return t;
}

// "DVCS_Fa18_inb" -> "fa18_inb"
static inline std::string periodToRunTagKey(const std::string& period){
    auto pos = period.find('_');
    if (pos==std::string::npos || pos+1>=period.size()) return toLower(period);
    return toLower(period.substr(pos+1));
}

// runTag "fa18_inb" -> "DVCS_Fa18_inb"
static inline std::string dvcsPeriodName(const std::string& runTag) {
    std::string cap = runTag;
    if (!cap.empty()) cap[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[0])));
    for (size_t i = 0; i + 1 < cap.size(); ++i) {
        if (cap[i] == '_' && i + 1 < cap.size())
            cap[i+1] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[i+1])));
    }
    return std::string("DVCS_") + cap;
}

static inline bool containsNoCase(const std::string& s, const std::string& pat){
    return toLower(s).find(toLower(pat)) != std::string::npos;
}

static inline std::string combinedKeyFor(const std::string& group) {
    // Map period groups to their combined key
    if (containsNoCase(group, "Fa18")) return "Fall2018";
    if (containsNoCase(group, "Sp18")) return "Spring2018";
    // 10.6_GeV is already a combined key; return as-is if asked for fallback
    if (containsNoCase(group, "10.6") || containsNoCase(group, "10_6") || containsNoCase(group, "10-6")) return "10.6_GeV";
    return "";
}

static inline std::vector<std::pair<double,double>>
uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which=='x') s.emplace(b.xBmin,b.xBmax);
        else if (which=='Q') s.emplace(b.Q2min,b.Q2max);
        else if (which=='t') s.emplace(b.tmin,b.tmax);
    }
    return {s.begin(), s.end()};
}
static inline int findIndex(const std::pair<double,double>& r,
                            const std::vector<std::pair<double,double>>& R){
    for (int i=0;i<(int)R.size();++i) if (R[i]==r) return i;
    return -1;
}
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}
static inline std::string key4s(int ix,int iQ,int it,int ip){
    std::ostringstream os; os<<"("<<ix<<","<<iQ<<","<<it<<","<<ip<<")"; return os.str();
}

// ────────── structures ──────────
struct HelCounts { long long plus=0, minus=0; };
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

struct Contam { double cp=0.0, ep=0.0, cm=0.0, em=0.0; }; // contamination & errors
using ContamTable = std::map<BinKey, Contam>;

struct CorrBin {
    double Np=0.0, ep=0.0;  // corrected +1 and its σ
    double Nm=0.0, em=0.0;  // corrected −1 and its σ
    double Nt=0.0, et=0.0;  // total and its σ (quadrature)
};

// ────────── JSON parsing utilities (tolerant) ──────────
static bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

// total_counts.json → groups -> (key)->helicity
static bool load_total_counts(const std::string& path, GroupCounts& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[pi0corr][ERROR] open "<<path<<" failed\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[pi0corr][ERROR] 'groups' not found in total_counts\n"; return false; }
    size_t br = s.find('{', gpos); if (br==std::string::npos) return false;

    int d=0; size_t i=br;
    for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string groupsObj = s.substr(br, i-br);

    size_t k=0;
    while (true) {
        size_t q1 = groupsObj.find('"', k); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);

        size_t binsS = groupsObj.find("\"bins\"", q2);
        if (binsS==std::string::npos) { k=q2+1; continue; }
        size_t bbr = groupsObj.find('{', binsS); if (bbr==std::string::npos){ k=q2+1; continue; }

        int d2=0; size_t j=bbr;
        for (; j<groupsObj.size(); ++j){ if(groupsObj[j]=='{') d2++; else if(groupsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string binsObj = groupsObj.substr(bbr, j-bbr);

        std::map<BinKey, HelCounts> table; size_t p=0;
        while (true) {
            size_t k1=binsObj.find('"', p); if (k1==std::string::npos) break;
            size_t k2=binsObj.find('"', k1+1); if (k2==std::string::npos) break;
            BinKey bk; if (!parse_tuple_key(binsObj.substr(k1+1,k2-k1-1), bk)){ p=k2+1; continue; }

            size_t vS = binsObj.find('{', k2); if (vS==std::string::npos) break;
            int d3=0; size_t m=vS;
            for (; m<binsObj.size(); ++m){ if(binsObj[m]=='{') d3++; else if(binsObj[m]=='}'){ d3--; if(!d3){ ++m; break;} } }
            std::string v = binsObj.substr(vS, m-vS);

            auto findLL=[&](const char* pat)->long long{
                size_t pos=v.find(pat); if(pos==std::string::npos) return 0;
                pos=v.find(':',pos); if(pos==std::string::npos) return 0;
                size_t a=pos+1; while(a<v.size() && isspace((unsigned char)v[a])) ++a;
                size_t b=a; while(b<v.size() && (isdigit((unsigned char>)v[b])||v[b]=='-')) ++b;
                try { return std::stoll(v.substr(a,b-a)); } catch(...) { return 0; }
            };

            HelCounts hc;
            hc.plus  = findLL("\"+1\"");
            hc.minus = findLL("\"-1\"");
            table[bk] = hc;
            p = m;
        }
        out[gname] = std::move(table);
        k = j;
    }
    return !out.empty();
}

// per-period contamination file
static bool load_contam_period(const std::string& path, ContamTable& out) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t bpos = s.find("\"bins\"");
    if (bpos==std::string::npos) return false;
    size_t br = s.find('{', bpos); if (br==std::string::npos) return false;

    int d=0; size_t i=br;
    for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    size_t p=0; int nb=0;
    while (true) {
        size_t q1=binsObj.find('"', p); if (q1==std::string::npos) break;
        size_t q2=binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        BinKey bk; if (!parse_tuple_key(binsObj.substr(q1+1,q2-q1-1), bk)) { p=q2+1; continue; }

        size_t vS=binsObj.find('{', q2); if (vS==std::string::npos) break;
        int d2=0; size_t j=vS;
        for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj=binsObj.substr(vS, j-vS);

        auto findD=[&](const std::string& pat)->double{
            size_t pos=obj.find(pat); if(pos==std::string::npos) return 0.0;
            pos=obj.find(':',pos); if(pos==std::string::npos) return 0.0;
            size_t a=pos+1, b=a;
            auto isnum=[](char c){return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E';};
            while(b<obj.size() && isnum(obj[b])) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };

        Contam c;
        // preferred nested object
        c.cp = findD("\"contamination\":{\"+1\":{\"value\"");
        c.ep = findD("\"contamination\":{\"+1\":{\"err\"");
        c.cm = findD("\"contamination\":{\"-1\":{\"value\"");
        c.em = findD("\"contamination\":{\"-1\":{\"err\"");
        // fallback to flat
        if (c.cp==0.0 && c.cm==0.0) {
            c.cp = findD("\"+1\":{\"value\"");
            c.ep = findD("\"+1\":{\"err\"");
            c.cm = findD("\"-1\":{\"value\"");
            c.em = findD("\"-1\":{\"err\"");
        }
        out[bk] = c;
        ++nb; p=j;
    }
    return !out.empty();
}

// combined contamination file: periods -> group -> bins
static bool load_contam_group_from_combined(const std::string& combined_path,
                                            const std::string& group,
                                            ContamTable& out)
{
    std::ifstream ifs(combined_path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t P = s.find("\"periods\""); if (P==std::string::npos) return false;
    size_t Pbr = s.find('{', P); if (Pbr==std::string::npos) return false;

    std::string gq = std::string("\"")+group+std::string("\"");
    size_t G = s.find(gq, Pbr); if (G==std::string::npos) return false;

    size_t binsK = s.find("\"bins\"", G); if (binsK==std::string::npos) return false;
    size_t br = s.find('{', binsK); if (br==std::string::npos) return false;

    int d=0; size_t i=br;
    for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    size_t p=0; int nb=0;
    while (true) {
        size_t q1=binsObj.find('"', p); if (q1==std::string::npos) break;
        size_t q2=binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        BinKey bk; if (!parse_tuple_key(binsObj.substr(q1+1,q2-q1-1), bk)) { p=q2+1; continue; }

        size_t vS=binsObj.find('{', q2); if (vS==std::string::npos) break;
        int d2=0; size_t j=vS;
        for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj=binsObj.substr(vS, j-vS);

        auto findD=[&](const std::string& pat)->double{
            size_t pos=obj.find(pat); if(pos==std::string::npos) return 0.0;
            pos=obj.find(':',pos); if(pos==std::string::npos) return 0.0;
            size_t a=pos+1, b=a;
            auto isnum=[](char c){return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E';};
            while(b<obj.size() && isnum(obj[b])) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };

        Contam c;
        c.cp = findD("\"contamination\":{\"+1\":{\"value\"");
        c.ep = findD("\"contamination\":{\"+1\":{\"err\"");
        c.cm = findD("\"contamination\":{\"-1\":{\"value\"");
        c.em = findD("\"contamination\":{\"-1\":{\"err\"");
        if (c.cp==0.0 && c.cm==0.0) {
            c.cp = findD("\"+1\":{\"value\"");
            c.ep = findD("\"+1\":{\"err\"");
            c.cm = findD("\"-1\":{\"value\"");
            c.em = findD("\"-1\":{\"err\"");
        }
        out[bk]=c; ++nb; p=j;
    }
    return !out.empty();
}

// ────────── tolerant contamination path resolver ──────────
static std::string find_contam_json(const std::string& dir, const std::string& period) {
    namespace fs = std::filesystem;

    auto exists_nonempty = [&](const std::string& p)->bool {
        std::error_code ec;
        if (!fs::exists(p, ec)) return false;
        return fs::file_size(p, ec) > 0;
    };

    std::vector<std::string> candidates;

    // 1) exact as given
    candidates.push_back((fs::path(dir) / ("contamination_" + period + ".json")).string());

    // 2) If it starts with DVCS_, try title-casing tail after DVCS_
    if (period.rfind("DVCS_", 0) == 0) {
        std::string tail = period.substr(5);
        candidates.push_back((fs::path(dir) / ("contamination_DVCS_" + titleCaseTail(tail) + ".json")).string());
        candidates.push_back((fs::path(dir) / ("contamination_DVCS_" + toLower(tail) + ".json")).string());
    } else {
        // 3) Not prefixed: try DVCS_<TitleCase(runTag)>
        const std::string dvcs = dvcsPeriodName(period);
        candidates.push_back((fs::path(dir) / ("contamination_" + dvcs + ".json")).string());

        // 4) DVCS_ + TitleCase tail from dvcs
        if (dvcs.rfind("DVCS_", 0) == 0) {
            std::string tail = dvcs.substr(5);
            candidates.push_back((fs::path(dir) / ("contamination_DVCS_" + titleCaseTail(tail) + ".json")).string());
            candidates.push_back((fs::path(dir) / ("contamination_DVCS_" + toLower(tail) + ".json")).string());
        }
    }

    // 5) case-insensitive fallback
    try {
        for (auto& de : fs::directory_iterator(dir)) {
            if (!de.is_regular_file()) continue;
            const std::string name = de.path().filename().string();
            const std::string want = "contamination_" + period + ".json";
            if (toLower(name) == toLower(want)) {
                candidates.push_back(de.path().string());
            }
        }
    } catch (...) {}

    for (const auto& p : candidates) {
        if (!p.empty() && exists_nonempty(p)) return p;
    }
    return ""; // not found
}

// ────────── JSON writers ──────────
static void write_group_json(const std::string& out_path,
                             int nPhi,
                             const std::vector<std::pair<double,double>>& xB_bins,
                             const std::vector<std::pair<double,double>>& Q2_bins,
                             const std::vector<std::pair<double,double>>& t_bins,
                             const std::map<BinKey, CorrBin>& table)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[pi0corr][ERROR] cannot write "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<xB_bins.size()
       <<", \"Q2_bins\": "<<Q2_bins.size()
       <<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
        const auto& b = kv.second;
        ofs<<"    \""<<key4s(ix,iQ,it,ip)<<"\": {"
           <<"\"helicity\":{"
           <<"\"+1\":{\"value\":"<<b.Np<<",\"err\":"<<b.ep<<"},"
           <<"\"-1\":{\"value\":"<<b.Nm<<",\"err\":"<<b.em<<"}"
           <<"},"
           <<"\"total\":{\"value\":"<<b.Nt<<",\"err\":"<<b.et<<"}"
           <<"}";
    }
    ofs<<"\n  }\n}\n";
}

static void write_master_json(const std::string& out_path,
                              int nPhi,
                              const std::vector<std::pair<double,double>>& xB_bins,
                              const std::vector<std::pair<double,double>>& Q2_bins,
                              const std::vector<std::pair<double,double>>& t_bins,
                              const std::map<std::string, std::map<BinKey, CorrBin>>& groups)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[pi0corr][ERROR] cannot write "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<xB_bins.size()
       <<", \"Q2_bins\": "<<Q2_bins.size()
       <<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"groups\": {\n";

    bool firstG=true;
    for (const auto& gkv : groups) {
        if (!firstG) ofs<<",\n"; firstG=false;
        ofs<<"    \""<<gkv.first<<"\": {\n";
        ofs<<"      \"bins\": {\n";
        bool firstB=true;
        for (const auto& kv : gkv.second) {
            if (!firstB) ofs<<",\n"; firstB=false;
            const auto& b = kv.second;
            int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
            ofs<<"        \""<<key4s(ix,iQ,it,ip)<<"\": {"
               <<"\"helicity\":{"
               <<"\"+1\":{\"value\":"<<b.Np<<",\"err\":"<<b.ep<<"},"
               <<"\"-1\":{\"value\":"<<b.Nm<<",\"err\":"<<b.em<<"}"
               <<"},"
               <<"\"total\":{\"value\":"<<b.Nt<<",\"err\":"<<b.et<<"}"
               <<"}";
        }
        ofs<<"\n      }\n"; // bins
        ofs<<"    }";
    }
    ofs<<"\n  }\n}\n";
}

// ────────── plotting ──────────
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin,b.xBmax) == xBrange) {
            qs.emplace(b.Q2min,b.Q2max);
            ts.emplace(b.tmin,b.tmax);
        }
    }
    Q2_list.assign(qs.begin(), qs.end());
    t_list.assign(ts.begin(), ts.end());
}

static void plot_group(
    const std::string& group,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<BinKey, CorrBin>& corrected,
    const std::map<BinKey, HelCounts>& raw_totals,
    const std::string& out_dir_plots)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    const auto PHI = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];
        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)Q2_slice.size();
        const int ncols = (int)t_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;
        std::string cname = Form("c_pi0corr_%s_xB%d", group.c_str(), ix);
        TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), W, H);

        // header pad
        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.92, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
        // grid pad
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.92);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // canvas title with xB
        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.40);
        head.DrawLatex(0.5, 0.55, Form("#pi^{0}-corrected vs raw counts  %s   x_{B} #in [%.3g, %.3g]",
                                       group.c_str(), xb.first, xb.second));

        for (int r = 0; r < nrows; ++r) {
            int iQ=findIndex(Q2_slice[r], Q2_bins); if (iQ<0) continue;
            for (int cc=0; cc<ncols; ++cc) {
                int it=findIndex(t_slice[cc], t_bins); if (it<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.16);
                gPad->SetLeftMargin(0.14);
                gPad->SetRightMargin(0.08);

                std::vector<double> X, Yp_raw, Ym_raw, Yp_corr, Ym_corr, eX, eYp_corr, eYm_corr;
                X.reserve(N_PHI_BINS);
                eX.assign(N_PHI_BINS, 0.0);

                // gather points for THIS slice
                double ymax = 0.0;
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    BinKey bk(ix, iQ, it, ip);
                    auto it_raw  = raw_totals.find(bk);
                    auto it_corr = corrected.find(bk);
                    if (it_raw == raw_totals.end() && it_corr == corrected.end()) continue;

                    // Require at least raw OR corrected to exist to add a point
                    const double Np_raw = (it_raw  != raw_totals.end()) ? (double)it_raw->second.plus  : 0.0;
                    const double Nm_raw = (it_raw  != raw_totals.end()) ? (double)it_raw->second.minus : 0.0;
                    const double Np_cor = (it_corr != corrected.end())  ? it_corr->second.Np : 0.0;
                    const double Nm_cor = (it_corr != corrected.end())  ? it_corr->second.Nm : 0.0;
                    const double ep_cor = (it_corr != corrected.end())  ? it_corr->second.ep : 0.0;
                    const double em_cor = (it_corr != corrected.end())  ? it_corr->second.em : 0.0;

                    // If everything is zero, skip (avoid 0-point graphs)
                    if (Np_raw==0.0 && Nm_raw==0.0 && Np_cor==0.0 && Nm_cor==0.0) continue;

                    X.push_back(PHI[ip]);
                    Yp_raw.push_back(Np_raw);
                    Ym_raw.push_back(Nm_raw);
                    Yp_corr.push_back(Np_cor);
                    Ym_corr.push_back(Nm_cor);
                    eYp_corr.push_back(ep_cor);
                    eYm_corr.push_back(em_cor);
                    ymax = std::max({ymax, Np_raw, Nm_raw, Np_cor+ep_cor, Nm_cor+em_cor});
                }

                if (X.empty()) { // nothing in this panel
                    continue;
                }

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, (ymax>0? ymax*1.25 : 1.0));
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Counts");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);

                // ---- corrected (AFTER) — filled red/blue with errors
                TGraphErrors* grPc = new TGraphErrors((int)X.size(), X.data(), Yp_corr.data(), eX.data(), eYp_corr.data()); // + helicity
                TGraphErrors* grMc = new TGraphErrors((int)X.size(), X.data(), Ym_corr.data(), eX.data(), eYm_corr.data()); // - helicity
                grPc->SetMarkerStyle(20); grPc->SetMarkerSize(1.0); grPc->SetLineWidth(2); grPc->SetMarkerColor(kRed);  grPc->SetLineColor(kRed);
                grMc->SetMarkerStyle(21); grMc->SetMarkerSize(1.0); grMc->SetLineWidth(2); grMc->SetMarkerColor(kBlue); grMc->SetLineColor(kBlue);
                grPc->Draw("P SAME");
                grMc->Draw("P SAME");

                // ---- raw (BEFORE) — open gray, dashed, no errors
                TGraph* grPr = new TGraph((int)X.size(), X.data(), Yp_raw.data()); // + helicity
                TGraph* grMr = new TGraph((int)X.size(), X.data(), Ym_raw.data()); // - helicity
                grPr->SetMarkerStyle(24); grPr->SetMarkerSize(0.9); grPr->SetLineStyle(2); grPr->SetMarkerColor(kGray+2); grPr->SetLineColor(kGray+2);
                grMr->SetMarkerStyle(25); grMr->SetMarkerSize(0.9); grMr->SetLineStyle(2); grMr->SetMarkerColor(kGray+2); grMr->SetLineColor(kGray+2);
                grPr->Draw("P SAME");
                grMr->Draw("P SAME");

                // legend
                TLegend* leg = new TLegend(0.54, 0.68, 0.92, 0.92);
                leg->SetBorderSize(1); leg->SetLineColor(kBlack); leg->SetFillStyle(0); leg->SetTextSize(0.035);
                leg->AddEntry(grPc, "+ helicity (corr.)", "p");
                leg->AddEntry(grMc, "- helicity (corr.)", "p");
                leg->AddEntry(grPr, "+ helicity (raw)",  "p");
                leg->AddEntry(grMr, "- helicity (raw)",  "p");
                leg->Draw();

                // subplot title with Q² and -t
                TLatex sub; sub.SetNDC(); sub.SetTextSize(0.045); sub.SetTextAlign(13);
                sub.DrawLatex(0.14, 0.96,
                    Form("Q^{2} #in [%.3g, %.3g],   -t #in [%.3g, %.3g]",
                         Q2_slice[r].first, Q2_slice[r].second, t_slice[cc].first, t_slice[cc].second));
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/plot_pi0corr_" << group << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ────────── correction core ──────────
static CorrBin make_corrected(const HelCounts& raw, const Contam& c) {
    // N_corr = N_raw * (1 - c),  sigma^2 = (sqrt(N_raw)*(1-c))^2 + (N_raw*sigma_c)^2
    auto one = [](double Nraw, double c, double ec)->std::pair<double,double>{
        if (Nraw <= 0.0) return {0.0, 0.0};
        if (c < 0.0) c = 0.0;
        const double val = Nraw * (1.0 - c);
        const double sig = std::sqrt( std::max(0.0, (std::sqrt(Nraw)*(1.0 - c))*(std::sqrt(Nraw)*(1.0 - c)) + (Nraw*ec)*(Nraw*ec)) );
        return {val, sig};
    };

    CorrBin out;
    std::tie(out.Np, out.ep) = one((double)raw.plus,  c.cp, c.ep);
    std::tie(out.Nm, out.em) = one((double)raw.minus, c.cm, c.em);
    out.Nt = out.Np + out.Nm;
    out.et = std::sqrt(out.ep*out.ep + out.em*out.em);
    return out;
}

} // namespace

// ======================= PUBLIC ENTRY ===========================
void compute_pi0_corrected_counts(
    const std::vector<std::string>& dvcs_periods,           // list like DVCS_Fa18_inb, ...
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,                   // output/jsons/total_counts.json
    const std::string& contamination_dir_counts,            // output/jsons/contamination
    const std::string& contamination_combined,              // output/jsons/pi0_contamination_combined.json
    const std::string& out_root_dir                         // "output"
) {
    namespace fs = std::filesystem;

    // 1) Load total counts groups
    GroupCounts group_counts;
    if (!load_total_counts(total_counts_json, group_counts)) {
        std::cerr << "[pi0corr][ERROR] Could not load total_counts.json; aborting.\n";
        return;
    }

    // 2) Build bin axes
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // 3) Prepare output dirs
    const fs::path json_dir = fs::path(out_root_dir) / "jsons";
    const fs::path plots_dir_root = fs::path(out_root_dir) / "pi0_corrected_plots";
    std::error_code ec;
    fs::create_directories(json_dir, ec);
    fs::create_directories(plots_dir_root, ec);

    // 4) Load contamination per period (tolerant path)
    std::map<std::string, ContamTable> contam_by_group; // key matches group name in total_counts

    auto load_period_contam = [&](const std::string& periodKey){
        // Try to resolve a per-period file. periodKey may be "fa18_inb" or "DVCS_Fa18_inb".
        std::vector<std::string> logical_periods;

        if (periodKey.rfind("DVCS_", 0) == 0) {
            logical_periods.push_back(periodKey);                       // DVCS_Fa18_inb
            logical_periods.push_back(periodToRunTagKey(periodKey));    // fa18_inb
        } else {
            logical_periods.push_back(periodKey);                        // fa18_inb
            logical_periods.push_back(dvcsPeriodName(periodKey));        // DVCS_Fa18_inb
        }

        for (const auto& per : logical_periods) {
            const std::string cand = find_contam_json(contamination_dir_counts, per);
            if (!cand.empty()) {
                ContamTable ct;
                if (load_contam_period(cand, ct)) {
                    contam_by_group[periodKey] = std::move(ct);
                    return true;
                }
            }
        }
        return false;
    };

    // Preload contamination for every group we have raw counts for.
    for (const auto& gkv : group_counts) {
        const std::string& group = gkv.first;

        // Combined groups come from the combined JSON
        if (group == "Spring2018" || group == "Fall2018" || group == "10.6_GeV") {
            ContamTable ct;
            if (load_contam_group_from_combined(contamination_combined, group, ct)) {
                contam_by_group[group] = std::move(ct);
            } else {
                std::cerr << "[pi0corr][WARN] contamination missing for " << group
                          << " (combined file \"" << contamination_combined << "\") — assuming c=0\n";
            }
            continue;
        }

        // Otherwise it's a plain period/runTag/DVCS_* name
        bool ok = load_period_contam(group);
        if (!ok) {
            // Fallback: use combined contamination for period's umbrella group (Fa18→Fall2018, Sp18→Spring2018)
            const std::string combo = combinedKeyFor(group);
            if (!combo.empty()) {
                ContamTable ct;
                if (load_contam_group_from_combined(contamination_combined, combo, ct)) {
                    contam_by_group[group] = std::move(ct);
                    std::cerr << "[pi0corr][INFO] Using combined contamination \"" << combo
                              << "\" as fallback for period \"" << group << "\"\n";
                    continue;
                }
            }
            std::cerr << "[pi0corr][WARN] contamination missing for " << group
                      << " (checked per-period and combined) — assuming c=0\n";
        }
    }

    // 5) Compute corrected counts per group, write one JSON per group, and plot
    std::map<std::string, std::map<BinKey, CorrBin>> all_groups_corrected;

    for (const auto& gkv : group_counts) {
        const std::string& group = gkv.first;
        const auto& raw_table    = gkv.second;

        const auto itC = contam_by_group.find(group);
        const ContamTable* C = (itC != contam_by_group.end()) ? &itC->second : nullptr;

        std::map<BinKey, CorrBin> corr_table;

        // For every bin present in the raw totals, apply contamination (default 0)
        for (const auto& kv : raw_table) {
            const BinKey& bk = kv.first;
            const HelCounts& raw = kv.second;

            Contam cc; // defaults to zeros
            if (C) {
                auto ic = C->find(bk);
                if (ic != C->end()) cc = ic->second;
            }
            CorrBin cb = make_corrected(raw, cc);
            corr_table[bk] = cb;
        }

        // write per-group JSON
        auto sanitize = [](std::string s)->std::string {
            for (char& c : s) if (c=='/' || c==' ' ) c = '_';
            return s;
        };
        const std::string out_group_json =
            (json_dir / ("pi0_corrected_counts_" + sanitize(group) + ".json")).string();

        write_group_json(out_group_json, N_PHI_BINS, xB_bins, Q2_bins, t_bins, corr_table);
        std::cout << "[pi0corr] Wrote " << out_group_json << "\n";

        // plot per-group
        const std::string out_plot_dir = (plots_dir_root / group).string();
        plot_group(group, binning_scheme, xB_bins, Q2_bins, t_bins, corr_table, raw_table, out_plot_dir);

        // stash for master file
        all_groups_corrected[group] = std::move(corr_table);
    }

    // 6) Write master JSON
    const std::string master =
        (json_dir / "pi0_corrected_counts_all_groups.json").string();
    write_master_json(master, N_PHI_BINS, xB_bins, Q2_bins, t_bins, all_groups_corrected);
    std::cout << "[pi0corr] Wrote " << master << "\n";
}