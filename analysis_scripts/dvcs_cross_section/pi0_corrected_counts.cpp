// pi0_corrected_counts.cpp (HARD-WIRED TO total_counts.json GROUPS)
//
// Behavior:
// - Read groups directly from total_counts.json and process ALL of them.
// - For combined groups ("Spring2018","Fall2018","10.6_GeV") use the combined contamination JSON.
// - For all other groups use per-period contamination files named
//     <contamination_dir_counts>/contamination_<group>.json
//   where <group> matches the group key EXACTLY as it appears in total_counts.json.
// - If a specific (ix,iQ2,it,ip) bin is missing in contamination, assume contamination=0.0±0.0.
//
// Inputs:
//   - total_counts.json                 (required; definitive list of groups to process)
//   - contamination_dir_counts/*.json   (per-period files; exact <group> naming)
//   - contamination_combined.json       (for "Spring2018","Fall2018","10.6_GeV")
//
// Outputs:
//   - <out_root_dir>/jsons/pi0_corrected_counts_<group>.json
//   - <out_root_dir>/jsons/pi0_corrected_counts_all_groups.json
//   - <out_root_dir>/pi0_corrected_plots/<group>/plot_pi0corr_<group>_xB_<ix>.png
// ------------------------------------------------------------------------------

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
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

constexpr int N_PHI_BINS = 12;
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)

[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[pi0corr][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

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
                            const std::vector<std::pair<double,double>>& R) {
    for (int i=0;i<(int)R.size();++i) if (R[i]==r) return i;
    return -1;
}

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}

static inline std::string key4s(int ix,int iQ,int it,int ip) {
    std::ostringstream os; os << "(" << ix << "," << iQ << "," << it << "," << ip << ")";
    return os.str();
}

struct HelCounts { long long plus=0, minus=0; };
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

struct Contam { double cp=0.0, ep=0.0, cm=0.0, em=0.0; };
using ContamTable = std::map<BinKey, Contam>;

struct CorrBin {
    double Np=0.0, ep=0.0;
    double Nm=0.0, em=0.0;
    double Nt=0.0, et=0.0;
};

struct BinningMeta {
    int phi_bins=0;
    size_t nx=0, nQ=0, nt=0;
};

// ---- tiny string-based JSON utilities (strict enough for our files) ----

static std::string objForKey(const std::string& s, const std::string& key) {
    size_t p = s.find(key);
    if (p==std::string::npos) fatal("Key '"+key+"' not found in JSON.");
    size_t br = s.find('{', p);
    if (br==std::string::npos) fatal("Malformed JSON after key '"+key+"'");
    int d=0; size_t i=br;
    for (; i<s.size(); ++i) {
        if (s[i]=='{') d++;
        else if (s[i]=='}') { d--; if (!d) { ++i; break; } }
    }
    if (d!=0) fatal("Unbalanced braces in JSON near key '"+key+"'");
    return s.substr(br, i-br);
}

static long long parseIntAfterColon_strict(const std::string& s, size_t cpos, const std::string& ctx) {
    size_t a = cpos + 1;
    while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
    size_t b = a;
    while (b < s.size() && (std::isdigit((unsigned char)s[b]) || s[b]=='-' || s[b]=='+')) ++b;
    try { return std::stoll(s.substr(a,b-a)); }
    catch(...) { fatal("Non-integer value in "+ctx); }
    return 0;
}

// Extract the raw "token" for the JSON value after a colon, and try to parse into double.
// Accepts either bare numbers or quoted numbers. If token is null/NaN/Inf, prints detailed context and fails.
static double parseDoubleAfterColon_diag(const std::string& whole_json,
                                         size_t colon_pos,
                                         const std::string& ctx,
                                         const std::string& file_hint,
                                         const std::string& group_hint,
                                         const std::string& bin_hint) {
    auto snippet_around = [&](size_t pos)->std::string{
        const size_t L = (pos>60? pos-60:0);
        const size_t R = std::min(whole_json.size(), pos+60);
        std::string sn = whole_json.substr(L, R-L);
        for (char& c : sn) if (std::iscntrl((unsigned char)c)) c = ' ';
        return sn;
    };

    // advance to first non-space after colon
    size_t a = colon_pos + 1;
    while (a < whole_json.size() && std::isspace((unsigned char)whole_json[a])) ++a;

    // read raw token
    std::string raw;
    size_t endpos = a;
    if (a < whole_json.size() && whole_json[a] == '"') {
        // quoted token
        ++endpos;
        while (endpos < whole_json.size() && whole_json[endpos] != '"') ++endpos;
        if (endpos >= whole_json.size()) {
            fatal("Unterminated quoted value in "+ctx+" near: "+snippet_around(a));
        }
        raw = whole_json.substr(a+1, endpos - (a+1));
        ++endpos;
    } else {
        // bare token: read until delimiter
        auto isdelim = [](char c)->bool{
            return std::isspace((unsigned char)c) || c==',' || c=='}' || c==']';
        };
        while (endpos < whole_json.size() && !isdelim(whole_json[endpos])) ++endpos;
        raw = whole_json.substr(a, endpos - a);
    }

    // Trim raw
    auto ltrim=[&](std::string& t){ size_t i=0; while (i<t.size() && std::isspace((unsigned char)t[i])) ++i; t.erase(0,i); };
    auto rtrim=[&](std::string& t){ size_t i=t.size(); while (i>0 && std::isspace((unsigned char)t[i-1])) --i; t.erase(i); };
    ltrim(raw); rtrim(raw);

    // Quick rejects we want to report verbosely
    std::string raw_lower = raw; for (char& c : raw_lower) c = (char)std::tolower((unsigned char)c);
    if (raw.empty() || raw_lower=="nan" || raw_lower=="inf" || raw_lower=="infinity"
        || raw_lower=="-inf" || raw_lower=="-infinity" || raw_lower=="null") {
        std::ostringstream em;
        em << "Non-numeric token in " << ctx
           << "  file=" << file_hint
           << "  group=" << group_hint
           << "  bin=" << bin_hint
           << "  raw_token=\"" << raw << "\""
           << "  json_snippet: ..." << snippet_around(a) << "...";
        fatal(em.str());
    }

    // Try std::stod on raw
    try {
        return std::stod(raw);
    } catch(...) {
        std::ostringstream em;
        em << "Non-numeric value in " << ctx
           << "  file=" << file_hint
           << "  group=" << group_hint
           << "  bin=" << bin_hint
           << "  raw_token=\"" << raw << "\""
           << "  json_snippet: ..." << snippet_around(a) << "...";
        fatal(em.str());
    }
    return 0.0; // unreachable
}

static bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

static BinningMeta load_meta(const std::string& s) {
    std::string meta = objForKey(s, "\"binning_meta\"");
    auto findN = [&](const char* k, const char* ctx)->int{
        size_t p = meta.find(k); if (p==std::string::npos) fatal(std::string("binning_meta missing ")+k);
        size_t c = meta.find(':', p); if (c==std::string::npos) fatal(std::string("binning_meta malformed for ")+k);
        return (int)parseIntAfterColon_strict(meta, c, ctx);
    };
    BinningMeta bm;
    bm.phi_bins = findN("\"phi_bins\"", "phi_bins");
    bm.nx       = (size_t)findN("\"xB_bins\"", "xB_bins");
    bm.nQ       = (size_t)findN("\"Q2_bins\"", "Q2_bins");
    bm.nt       = (size_t)findN("\"t_bins\"",  "t_bins");
    if (bm.phi_bins != N_PHI_BINS) fatal("phi_bins mismatch: expected 12");
    return bm;
}

// total_counts.json -> definitive groups + bins
static GroupCounts load_total_counts_STRICT(const std::string& path, BinningMeta& out_meta,
                                            std::vector<std::string>& group_order) {
    std::ifstream ifs(path);
    if (!ifs) fatal(std::string("Cannot open total_counts_json: ")+path);
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    out_meta = load_meta(s);
    std::string groupsObj = objForKey(s, "\"groups\"");

    GroupCounts out;
    size_t k=0;
    while (true) {
        size_t q1 = groupsObj.find('"', k); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) fatal("Malformed group name in total_counts.json");
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);
        group_order.push_back(gname);

        size_t binsS = groupsObj.find("\"bins\"", q2);
        if (binsS==std::string::npos) fatal("Group "+gname+" has no 'bins' object in total_counts.json");
        int d=0; size_t br = groupsObj.find('{', binsS); if (br==std::string::npos) fatal("Malformed 'bins' in group "+gname);
        size_t j=br;
        for (; j<groupsObj.size(); ++j) {
            if(groupsObj[j]=='{') d++;
            else if(groupsObj[j]=='}'){ d--; if(!d){ ++j; break; } }
        }
        std::string binsObj = groupsObj.substr(br, j-br);

        std::map<BinKey, HelCounts> table; size_t p=0; int nbins=0;
        while (true) {
            size_t k1=binsObj.find('"', p); if (k1==std::string::npos) break;
            size_t k2=binsObj.find('"', k1+1); if (k2==std::string::npos) fatal("Malformed bin key in "+gname);
            BinKey bk; if (!parse_tuple_key(binsObj.substr(k1+1,k2-k1-1), bk)) fatal("Bad bin tuple key in "+gname);

            size_t vS = binsObj.find('{', k2); if (vS==std::string::npos) fatal("Missing bin object in "+gname);
            int d2=0; size_t m=vS;
            for (; m<binsObj.size(); ++m) {
                if(binsObj[m]=='{') d2++;
                else if(binsObj[m]=='}'){ d2--; if(!d2){ ++m; break;} }
            }
            std::string v = binsObj.substr(vS, m-vS);

            size_t pPlus = v.find("\"+1\""); if (pPlus==std::string::npos) fatal("Missing +1 in helicity for "+gname);
            pPlus = v.find(':', pPlus); if (pPlus==std::string::npos) fatal("Malformed +1 in "+gname);
            long long Np = parseIntAfterColon_strict(v, pPlus, gname+" (+1)");

            size_t pMinus = v.find("\"-1\""); if (pMinus==std::string::npos) fatal("Missing -1 in helicity for "+gname);
            pMinus = v.find(':', pMinus); if (pMinus==std::string::npos) fatal("Malformed -1 in "+gname);
            long long Nm = parseIntAfterColon_strict(v, pMinus, gname+" (-1)");

            table[bk] = HelCounts{Np,Nm};
            p = m; ++nbins;
        }
        if (nbins==0) fatal("Group "+gname+" has zero bins in total_counts.json");
        out[gname] = std::move(table);
        k = j;
    }
    if (out.empty()) fatal("No groups parsed from total_counts.json");
    return out;
}

// Per-period contamination file (exact <group> name)
static ContamTable load_contam_period_STRICT(const std::string& path,
                                             BinningMeta& out_meta,
                                             const std::string& group) {
    std::ifstream ifs(path);
    if (!ifs) fatal(std::string("Cannot open contamination file: ")+path);
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    out_meta = load_meta(s);
    std::string binsObj = objForKey(s, "\"bins\"");

    // Small helper: extract the {...} block that immediately follows a label like "\"+1\"" or "\"-1\""
    auto extract_block_after_label = [](const std::string& src, const char* label)->std::string {
        size_t p = src.find(label);
        if (p == std::string::npos) return std::string();
        size_t br = src.find('{', p);
        if (br == std::string::npos) return std::string();
        int d = 0; size_t i = br;
        for (; i < src.size(); ++i) {
            if (src[i] == '{') ++d;
            else if (src[i] == '}') { --d; if (!d) { ++i; break; } }
        }
        if (d != 0) return std::string();
        return src.substr(br, i - br);
    };

    // parse a "value" or "err" inside a small block; with diagnostics
    auto parse_field = [&](const std::string& block, const char* key,
                           const std::string& ctx,
                           const std::string& file_hint,
                           const std::string& group_hint,
                           const std::string& bin_hint)->double {
        size_t p = block.find(key);
        if (p == std::string::npos) {
            std::ostringstream em; em << "Missing " << key << " in " << ctx
                                      << "  file=" << file_hint
                                      << "  group=" << group_hint
                                      << "  bin=" << bin_hint;
            fatal(em.str());
        }
        size_t c = block.find(':', p);
        if (c == std::string::npos) {
            std::ostringstream em; em << "Malformed " << key << " in " << ctx
                                      << "  file=" << file_hint
                                      << "  group=" << group_hint
                                      << "  bin=" << bin_hint;
            fatal(em.str());
        }
        return parseDoubleAfterColon_diag(block, c, ctx, file_hint, group_hint, bin_hint);
    };

    ContamTable out; size_t p = 0; int nb = 0;
    while (true) {
        size_t q1 = binsObj.find('"', p); if (q1 == std::string::npos) break;
        size_t q2 = binsObj.find('"', q1+1); if (q2 == std::string::npos) fatal("Malformed bin key in contamination JSON");
        std::string binKeyStr = binsObj.substr(q1+1, q2-q1-1);
        BinKey bk; if (!parse_tuple_key(binKeyStr, bk)) fatal("Bad bin tuple key in contamination JSON");

        size_t vS = binsObj.find('{', q2); if (vS == std::string::npos) fatal("Missing bin object in contamination JSON");
        int d2 = 0; size_t j = vS;
        for (; j < binsObj.size(); ++j) {
            if (binsObj[j] == '{') ++d2;
            else if (binsObj[j] == '}') { --d2; if (!d2) { ++j; break; } }
        }
        std::string obj = binsObj.substr(vS, j - vS);

        // Narrow to the contamination object first
        std::string contamObj = objForKey(obj, "\"contamination\"");
        // Then isolate +1 and -1 blocks
        std::string blockP = extract_block_after_label(contamObj, "\"+1\"");
        std::string blockM = extract_block_after_label(contamObj, "\"-1\"");
        if (blockP.empty() || blockM.empty()) {
            std::ostringstream em; em << "Missing helicity block(s) in contamination object"
                                      << "  file=" << path << "  group=" << group
                                      << "  bin=" << binKeyStr;
            fatal(em.str());
        }

        Contam c;
        c.cp = parse_field(blockP, "\"value\"", "contamination(+1).value", path, group, binKeyStr);
        c.ep = parse_field(blockP, "\"err\"",   "contamination(+1).err",   path, group, binKeyStr);
        c.cm = parse_field(blockM, "\"value\"", "contamination(-1).value", path, group, binKeyStr);
        c.em = parse_field(blockM, "\"err\"",   "contamination(-1).err",   path, group, binKeyStr);

        out[bk] = c; ++nb; p = j;
    }
    if (nb == 0) fatal("No bins parsed in contamination file: " + path);
    return out;
}

// Combined contamination (exact group name inside "periods")
static ContamTable load_contam_group_from_combined_STRICT(const std::string& combined_path,
                                                          const std::string& group,
                                                          BinningMeta& out_meta) {
    std::ifstream ifs(combined_path);
    if (!ifs) fatal(std::string("Cannot open combined contamination: ")+combined_path);
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    out_meta = load_meta(s);

    std::string periodsObj = objForKey(s, "\"periods\"");
    std::string gq = std::string("\"")+group+std::string("\"");
    size_t G = periodsObj.find(gq);
    if (G==std::string::npos) fatal("Group '"+group+"' not found in combined contamination JSON");

    size_t binsK = periodsObj.find("\"bins\"", G);
    if (binsK==std::string::npos) fatal("Group '"+group+"' has no 'bins' in combined contamination JSON");
    int d=0; size_t br = periodsObj.find('{', binsK); if (br==std::string::npos) fatal("Malformed 'bins' for group "+group);
    size_t i=br; for (; i<periodsObj.size(); ++i) {
        if (periodsObj[i]=='{') ++d;
        else if (periodsObj[i]=='}') { --d; if (!d) { ++i; break; } }
    }
    std::string binsObj = periodsObj.substr(br, i-br);

    auto extract_block_after_label = [](const std::string& src, const char* label)->std::string {
        size_t p = src.find(label);
        if (p == std::string::npos) return std::string();
        size_t br = src.find('{', p);
        if (br == std::string::npos) return std::string();
        int d = 0; size_t i = br;
        for (; i < src.size(); ++i) {
            if (src[i] == '{') ++d;
            else if (src[i] == '}') { --d; if (!d) { ++i; break; } }
        }
        if (d != 0) return std::string();
        return src.substr(br, i - br);
    };

    auto parse_field = [&](const std::string& block, const char* key,
                           const std::string& ctx,
                           const std::string& file_hint,
                           const std::string& group_hint,
                           const std::string& bin_hint)->double {
        size_t p = block.find(key);
        if (p == std::string::npos) {
            std::ostringstream em; em << "Missing " << key << " in " << ctx
                                      << "  file=" << file_hint
                                      << "  group=" << group_hint
                                      << "  bin=" << bin_hint;
            fatal(em.str());
        }
        size_t c = block.find(':', p);
        if (c == std::string::npos) {
            std::ostringstream em; em << "Malformed " << key << " in " << ctx
                                      << "  file=" << file_hint
                                      << "  group=" << group_hint
                                      << "  bin=" << bin_hint;
            fatal(em.str());
        }
        return parseDoubleAfterColon_diag(block, c, ctx, file_hint, group_hint, bin_hint);
    };

    ContamTable out; size_t p = 0; int nb = 0;
    while (true) {
        size_t q1 = binsObj.find('"', p); if (q1 == std::string::npos) break;
        size_t q2 = binsObj.find('"', q1+1); if (q2 == std::string::npos) fatal("Malformed bin key in combined contamination");
        std::string binKeyStr = binsObj.substr(q1+1, q2-q1-1);
        BinKey bk; if (!parse_tuple_key(binKeyStr, bk)) fatal("Bad bin tuple in combined contamination");

        size_t vS = binsObj.find('{', q2); if (vS == std::string::npos) fatal("Missing bin object in combined contamination");
        int d2 = 0; size_t j = vS;
        for (; j < binsObj.size(); ++j) {
            if (binsObj[j] == '{') ++d2;
            else if (binsObj[j] == '}') { --d2; if (!d2) { ++j; break; } }
        }
        std::string obj = binsObj.substr(vS, j - vS);

        // Narrow to contamination, then (+1) and (-1) subblocks
        std::string contamObj = objForKey(obj, "\"contamination\"");
        std::string blockP    = extract_block_after_label(contamObj, "\"+1\"");
        std::string blockM    = extract_block_after_label(contamObj, "\"-1\"");
        if (blockP.empty() || blockM.empty()) {
            std::ostringstream em; em << "Missing helicity block(s) in combined contamination object"
                                      << "  file=" << combined_path << "  group=" << group
                                      << "  bin=" << binKeyStr;
            fatal(em.str());
        }

        Contam c;
        c.cp = parse_field(blockP, "\"value\"", "combined contamination(+1).value", combined_path, group, binKeyStr);
        c.ep = parse_field(blockP, "\"err\"",   "combined contamination(+1).err",   combined_path, group, binKeyStr);
        c.cm = parse_field(blockM, "\"value\"", "combined contamination(-1).value", combined_path, group, binKeyStr);
        c.em = parse_field(blockM, "\"err\"",   "combined contamination(-1).err",   combined_path, group, binKeyStr);

        out[bk] = c; ++nb; p = j;
    }
    if (nb == 0) fatal("No bins parsed for group "+group+" in combined contamination JSON");
    return out;
}

// ---- plotting (unchanged) ----

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
    if (ec) fatal("Cannot create directory: "+out_dir_plots+" ("+ec.message()+")");

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

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.92, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.92);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

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
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.08);

                std::vector<double> X, Yp_raw, Ym_raw, Yp_corr, Ym_corr, eX, eYp_corr, eYm_corr;
                eX.assign(N_PHI_BINS, 0.0);

                double ymax = 0.0;
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    BinKey bk(ix, iQ, it, ip);
                    auto it_raw  = raw_totals.find(bk);
                    auto it_corr = corrected.find(bk);

                    const double Np_raw = (it_raw  != raw_totals.end()) ? (double)it_raw->second.plus  : 0.0;
                    const double Nm_raw = (it_raw  != raw_totals.end()) ? (double)it_raw->second.minus : 0.0;
                    const double Np_cor = (it_corr != corrected.end())  ? it_corr->second.Np : 0.0;
                    const double Nm_cor = (it_corr != corrected.end())  ? it_corr->second.Nm : 0.0;
                    const double ep_cor = (it_corr != corrected.end())  ? it_corr->second.ep : 0.0;
                    const double em_cor = (it_corr != corrected.end())  ? it_corr->second.em : 0.0;

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

                if (X.empty()) continue;

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, (ymax>0? ymax*1.25 : 1.0));
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Counts");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);

                TGraphErrors* grPc = new TGraphErrors((int)X.size(), X.data(), Yp_corr.data(), eX.data(), eYp_corr.data());
                TGraphErrors* grMc = new TGraphErrors((int)X.size(), X.data(), Ym_corr.data(), eX.data(), eYm_corr.data());
                grPc->SetMarkerStyle(20); grPc->SetMarkerSize(1.0); grPc->SetLineWidth(2); grPc->SetMarkerColor(kRed);  grPc->SetLineColor(kRed);
                grMc->SetMarkerStyle(21); grMc->SetMarkerSize(1.0); grMc->SetLineWidth(2); grMc->SetMarkerColor(kBlue); grMc->SetLineColor(kBlue);
                grPc->Draw("P SAME");
                grMc->Draw("P SAME");

                TGraph* grPr = new TGraph((int)X.size(), X.data(), Yp_raw.data());
                TGraph* grMr = new TGraph((int)X.size(), X.data(), Ym_raw.data());
                grPr->SetMarkerStyle(24); grPr->SetMarkerSize(0.9); grPr->SetLineStyle(2); grPr->SetMarkerColor(kGray+2); grPr->SetLineColor(kGray+2);
                grMr->SetMarkerStyle(25); grMr->SetMarkerSize(0.9); grMr->SetLineStyle(2); grMr->SetMarkerColor(kGray+2); grMr->SetLineColor(kGray+2);
                grPr->Draw("P SAME");
                grMr->Draw("P SAME");

                TLegend* leg = new TLegend(0.54, 0.68, 0.92, 0.92);
                leg->SetBorderSize(1); leg->SetLineColor(kBlack); leg->SetFillStyle(0); leg->SetTextSize(0.035);
                leg->AddEntry(grPc, "+ helicity (corr.)", "p");
                leg->AddEntry(grMc, "- helicity (corr.)", "p");
                leg->AddEntry(grPr, "+ helicity (raw)",  "p");
                leg->AddEntry(grMr, "- helicity (raw)",  "p");
                leg->Draw();

                TLatex sub; sub.SetNDC(); sub.SetTextSize(0.045); sub.SetTextAlign(13);
                sub.DrawLatex(0.14, 0.96,
                    Form("Q^{2} #in [%.3g, %.3g],   -t #in [%.3g, %.3g]",
                         Q2_slice[r].first, Q2_slice[r].second,
                         t_slice[cc].first, t_slice[cc].second));
            }
        }

        const std::string fpath =
            out_dir_plots + "/plot_pi0corr_" + group + "_xB_" + std::to_string(ix) + ".png";
        c->SaveAs(fpath.c_str());

        std::error_code fec;
        const bool exists = std::filesystem::exists(fpath, fec);
        const auto size   = exists ? std::filesystem::file_size(fpath, fec) : 0ULL;
        if (!exists || size == 0 || fec) {
            delete c;
            std::ostringstream em;
            em << "Failed to save plot: " << fpath
               << " (exists=" << std::boolalpha << exists
               << ", size=" << size
               << ", ec=" << (fec ? fec.message() : "ok") << ")";
            fatal(em.str());
        }
        delete c;
    }
}

static CorrBin make_corrected_STRICT(const HelCounts& raw, const Contam& c,
                                     const std::string& group, const std::string& binKeyStr) {
    auto one = [&](double Nraw, double val, double err, const char* hel)->std::pair<double,double>{
        if (Nraw < 0.0) fatal("Negative raw count in "+group+" bin "+binKeyStr);
        if (val < 0.0) fatal(std::string("Negative contamination for ")+hel+" in "+group+" bin "+binKeyStr);
        if (!std::isfinite(val) || !std::isfinite(err)) fatal("Non-finite contamination in "+group+" bin "+binKeyStr);
        const double Ncorr = Nraw * (1.0 - val);
        const double sigma = std::sqrt(std::max(0.0,
            (std::sqrt(Nraw)*(1.0 - val))*(std::sqrt(Nraw)*(1.0 - val)) + (Nraw*err)*(Nraw*err)));
        return {Ncorr, sigma};
    };

    CorrBin out;
    std::tie(out.Np, out.ep) = one((double)raw.plus,  c.cp, c.ep, "+1");
    std::tie(out.Nm, out.em) = one((double)raw.minus, c.cm, c.em, "-1");
    out.Nt = out.Np + out.Nm;
    out.et = std::sqrt(out.ep*out.ep + out.em*out.em);
    return out;
}

static void write_group_json(const std::string& out_path,
                             int nPhi,
                             size_t nx, size_t nQ, size_t nt,
                             const std::map<BinKey, CorrBin>& table) {
    std::ofstream ofs(out_path);
    if (!ofs) fatal(std::string("Cannot write: ")+out_path);
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<nx
       <<", \"Q2_bins\": "<<nQ
       <<", \"t_bins\": "<<nt<<"},\n";
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
                              size_t nx, size_t nQ, size_t nt,
                              const std::map<std::string, std::map<BinKey, CorrBin>>& groups) {
    std::ofstream ofs(out_path);
    if (!ofs) fatal(std::string("Cannot write: ")+out_path);
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<nx
       <<", \"Q2_bins\": "<<nQ
       <<", \"t_bins\": "<<nt<<"},\n";
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
        ofs<<"\n      }\n";
        ofs<<"    }";
    }
    ofs<<"\n  }\n}\n";
}

} // end anonymous namespace

void compute_pi0_corrected_counts(
    const std::vector<std::string>& /*dvcs_periods (ignored on purpose)*/,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir_counts,
    const std::string& contamination_combined,
    const std::string& out_root_dir
) {
    namespace fs = std::filesystem;

    BinningMeta totals_meta;
    std::vector<std::string> group_order;
    GroupCounts group_counts = load_total_counts_STRICT(total_counts_json, totals_meta, group_order);

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    if (xB_bins.size()!=totals_meta.nx || Q2_bins.size()!=totals_meta.nQ || t_bins.size()!=totals_meta.nt)
        fatal("Binning scheme sizes mismatch total_counts binning_meta.");

    const fs::path json_dir       = fs::path(out_root_dir) / "jsons";
    const fs::path plots_dir_root = fs::path(out_root_dir) / "pi0_corrected_plots";
    std::error_code ec;
    if (!fs::create_directories(json_dir, ec) && ec)
        fatal(std::string("Cannot create JSON output dir: ")+json_dir.string()+" ("+ec.message()+")");
    ec.clear();
    if (!fs::create_directories(plots_dir_root, ec) && ec)
        fatal(std::string("Cannot create plots root: ")+plots_dir_root.string()+" ("+ec.message()+")");

    std::map<std::string, ContamTable> contam_by_group;

    auto load_period_contam = [&](const std::string& group)->void{
        const fs::path f = fs::path(contamination_dir_counts) / ("contamination_" + group + ".json");
        BinningMeta cm;
        std::cout<<"[pi0corr] Reading per-period contamination: "<<f<<"  group="<<group<<std::endl;
        ContamTable ct = load_contam_period_STRICT(f.string(), cm, group);
        if (cm.phi_bins!=totals_meta.phi_bins || cm.nx!=totals_meta.nx || cm.nQ!=totals_meta.nQ || cm.nt!=totals_meta.nt)
            fatal("Binning meta mismatch between contamination("+group+") and total_counts.");
        contam_by_group[group] = std::move(ct);
    };

    auto load_combined_contam = [&](const std::string& group)->void{
        BinningMeta cm;
        std::cout<<"[pi0corr] Reading combined contamination: "<<contamination_combined<<"  group="<<group<<std::endl;
        ContamTable ct = load_contam_group_from_combined_STRICT(contamination_combined, group, cm);
        if (cm.phi_bins!=totals_meta.phi_bins || cm.nx!=totals_meta.nx || cm.nQ!=totals_meta.nQ || cm.nt!=totals_meta.nt)
            fatal("Binning meta mismatch between combined contamination("+group+") and total_counts.");
        contam_by_group[group] = std::move(ct);
    };

    for (const auto& gname : group_order) {
        if (gname == "Spring2018" || gname == "Fall2018" || gname == "10.6_GeV") {
            load_combined_contam(gname);
        } else {
            load_period_contam(gname);
        }
    }

    std::map<std::string, std::map<BinKey, CorrBin>> all_groups_corrected;

    for (const auto& gname : group_order) {
        const auto& raw_table = group_counts.at(gname);

        const auto itC = contam_by_group.find(gname);
        if (itC == contam_by_group.end())
            fatal("No contamination table loaded for group '"+gname+"'");
        const ContamTable& C = itC->second;

        std::map<BinKey, CorrBin> corr_table;

        for (const auto& kv : raw_table) {
            const BinKey& bk = kv.first;
            const HelCounts& raw = kv.second;

            auto ic = C.find(bk);
            Contam c;
            if (ic != C.end()) {
                c = ic->second;
            } else {
                c = Contam{0.0, 0.0, 0.0, 0.0};
            }

            CorrBin cb = make_corrected_STRICT(raw, c, gname,
                                               key4s(std::get<0>(bk),std::get<1>(bk),std::get<2>(bk),std::get<3>(bk)));
            corr_table[bk] = cb;
        }

        auto path_sanitize = [](std::string s)->std::string {
            for (char& c : s) if (c=='/' || c==' ' ) c = '_';
            return s;
        };
        const std::string out_group_json =
            (json_dir / ("pi0_corrected_counts_" + path_sanitize(gname) + ".json")).string();

        write_group_json(out_group_json, N_PHI_BINS, xB_bins.size(), Q2_bins.size(), t_bins.size(), corr_table);
        std::cout << "[pi0corr] Wrote " << out_group_json << "\n";

        const std::string out_plot_dir = (plots_dir_root / gname).string();
        plot_group(gname, binning_scheme, xB_bins, Q2_bins, t_bins, corr_table, raw_table, out_plot_dir);

        all_groups_corrected[gname] = std::move(corr_table);
    }

    const std::string master = (json_dir / "pi0_corrected_counts_all_groups.json").string();
    write_master_json(master, N_PHI_BINS, xB_bins.size(), Q2_bins.size(), t_bins.size(), all_groups_corrected);
    std::cout << "[pi0corr] Wrote " << master << "\n";
}