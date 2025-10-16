// pi0_corrected_counts.cpp — π0-corrected counts (helicity-resolved)

#include "pi0_corrected_counts.h"
#include "load_binning_scheme.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// ---------------- bin helpers (consistent with contamination stage) ----------------
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
static constexpr int N_PHI_BINS = 12;

static inline std::string keyStr(const BinKey& k) {
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os;
    os << "(" << a << "," << b << "," << c << "," << d << ")";
    return os.str();
}

static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin,b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min,b.Q2max);
        else if (which == 't') s.emplace(b.tmin,b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

// ---------------- data containers ----------------
struct HelCounts {
    long long plus = 0;
    long long minus = 0;
};
struct HelContam {
    double c_plus = 0.0;
    double c_plus_err = 0.0;
    double c_minus = 0.0;
    double c_minus_err = 0.0;
};
struct HelCorr {
    double val_plus = 0.0;
    double err_plus = 0.0;
    double val_minus = 0.0;
    double err_minus = 0.0;
};

// Store everything per-bin
struct BinRecord {
    HelCounts N_data;
    HelContam contam;
    HelCorr   N_corr;
};

// ---------------- JSON parsing helpers ----------------
static double parseNumberAfterColon(const std::string& s, size_t pos, double fallback = 0.0) {
    pos = s.find(':', pos);
    if (pos == std::string::npos) return fallback;
    size_t i = pos + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    auto isnum=[&](char c){return (std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E');};
    while (j < s.size() && isnum(s[j])) ++j;
    try { return std::stod(s.substr(i, j - i)); } catch (...) { return fallback; }
}
static long long parseIntAfterColon(const std::string& s, size_t pos, long long fallback = 0) {
    pos = s.find(':', pos);
    if (pos == std::string::npos) return fallback;
    size_t i = pos + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    while (j < s.size() && (std::isdigit((unsigned char)s[j]) || s[j]=='-')) ++j;
    try { return std::stoll(s.substr(i, j - i)); } catch (...) { return fallback; }
}

// ---------------- load total_counts.json ----------------
static std::map<BinKey, HelCounts> load_total_counts_for_period(
    const std::string& total_counts_path,
    const std::string& period)
{
    std::ifstream ifs(total_counts_path);
    if (!ifs) {
        std::cerr << "[pi0_corr][ERROR] Cannot open total_counts: " << total_counts_path << "\n";
        return {};
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    // The key we want (group name)
    const std::string groupKey = "\""+period+"\"";
    size_t pos = s.find(groupKey);
    if (pos == std::string::npos) {
        std::cerr << "[pi0_corr][WARN] total_counts has no group '" << period << "'.\n";
        return {};
    }

    // Find "bins" inside this group
    pos = s.find("\"bins\"", pos);
    if (pos == std::string::npos) return {};
    pos = s.find('{', pos);
    if (pos == std::string::npos) return {};
    int depth = 0;
    size_t start = pos;
    size_t end = start;
    for (; end < s.size(); ++end) {
        if (s[end] == '{') depth++;
        else if (s[end] == '}') {
            depth--;
            if (depth == 0) { ++end; break; }
        }
    }
    std::string binsSection = s.substr(start, end - start);

    std::map<BinKey, HelCounts> result;
    size_t kpos = 0;
    while (true) {
        size_t ks = binsSection.find('"', kpos);
        if (ks == std::string::npos) break;
        size_t ke = binsSection.find('"', ks+1);
        if (ke == std::string::npos) break;
        std::string key = binsSection.substr(ks+1, ke-ks-1);
        int ix,iQ,it,ip;
        if (std::sscanf(key.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4){kpos=ke+1;continue;}

        size_t subObjStart = binsSection.find('{', ke);
        if (subObjStart == std::string::npos) break;
        int d2=0; size_t j=subObjStart;
        for (; j < binsSection.size(); ++j) {
            if (binsSection[j]=='{') d2++;
            else if (binsSection[j]=='}') {
                d2--;
                if (d2==0) {++j;break;}
            }
        }
        std::string obj = binsSection.substr(subObjStart, j-subObjStart);

        HelCounts hc;
        hc.plus  = parseIntAfterColon(obj, obj.find("\"+1\""));
        hc.minus = parseIntAfterColon(obj, obj.find("\"-1\""));
        result[BinKey(ix,iQ,it,ip)] = hc;

        kpos = j;
    }
    std::cout << "[pi0_corr] Loaded total_counts for " << period << " (" << result.size() << " bins)\n";
    return result;
}

// ---------------- load contamination JSON ----------------
static std::map<BinKey, HelContam> load_contamination_for_period(
    const std::string& contamination_dir,
    const std::string& period)
{
    std::string fname = (fs::path(contamination_dir)/("contamination_"+period+".json")).string();
    std::ifstream ifs(fname);
    if (!ifs) {
        std::cerr << "[pi0_corr][WARN] Missing contamination JSON: " << fname << "\n";
        return {};
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    std::map<BinKey, HelContam> out;

    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) return {};
    pos = s.find('{', pos);
    if (pos == std::string::npos) return {};
    int depth = 0; size_t i=pos;
    for (; i<s.size(); ++i) {
        if (s[i]=='{') depth++;
        else if (s[i]=='}') {depth--; if(!depth){++i;break;}}
    }
    std::string binsObj = s.substr(pos,i-pos);

    size_t kpos=0;
    while(true){
        size_t ks=binsObj.find('"',kpos);
        if(ks==std::string::npos)break;
        size_t ke=binsObj.find('"',ks+1);
        if(ke==std::string::npos)break;
        std::string key=binsObj.substr(ks+1,ke-ks-1);
        int ix,iQ,it,ip;
        if(std::sscanf(key.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4){kpos=ke+1;continue;}
        size_t subStart=binsObj.find('{',ke);
        if(subStart==std::string::npos)break;
        int d2=0; size_t j=subStart;
        for(;j<binsObj.size();++j){if(binsObj[j]=='{')d2++;else if(binsObj[j]=='}'){d2--;if(!d2){++j;break;}}}
        std::string obj=binsObj.substr(subStart,j-subStart);
        HelContam hc;
        hc.c_plus     = parseNumberAfterColon(obj, obj.find("\"+1\":{\"value\""));
        hc.c_plus_err = parseNumberAfterColon(obj, obj.find("\"+1\":{\"err\""));
        hc.c_minus    = parseNumberAfterColon(obj, obj.find("\"-1\":{\"value\""));
        hc.c_minus_err= parseNumberAfterColon(obj, obj.find("\"-1\":{\"err\""));
        out[BinKey(ix,iQ,it,ip)] = hc;
        kpos=j;
    }
    std::cout << "[pi0_corr] Loaded contamination for " << period << " (" << out.size() << " bins)\n";
    return out;
}

// ---------------- main correction computation ----------------
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& out_root_dir)
{
    fs::path json_out_dir = fs::path(out_root_dir)/"jsons";
    fs::create_directories(json_out_dir);

    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');

    std::map<std::string, std::map<BinKey, BinRecord>> allPeriods;

    for(const auto& period : periods) {
        std::string usePeriod = period;
        if(period == "DVCS_Fa18_inb_supp")
            usePeriod = "DVCS_Fa18_inb"; // reuse Fa18_inb contamination for supp

        auto counts  = load_total_counts_for_period(total_counts_json, period);
        auto contam  = load_contamination_for_period(contamination_dir, usePeriod);

        std::map<BinKey, BinRecord> bins;

        for(const auto& kv : counts){
            BinKey k = kv.first;
            const HelCounts& nd = kv.second;
            HelContam hc{};
            auto itc = contam.find(k);
            if(itc != contam.end()) hc = itc->second;

            BinRecord br;
            br.N_data = nd;
            br.contam = hc;

            // Compute corrected counts per helicity
            auto correct_one = [&](long long N, double c, double cerr, double& Nc, double& Nc_err){
                Nc = N * (1.0 - c);
                Nc_err = std::sqrt( std::pow(std::sqrt((double)N)*(1.0-c),2) + std::pow(N*cerr,2) );
            };
            correct_one(nd.plus,  hc.c_plus,  hc.c_plus_err,  br.N_corr.val_plus,  br.N_corr.err_plus);
            correct_one(nd.minus, hc.c_minus, hc.c_minus_err, br.N_corr.val_minus, br.N_corr.err_minus);

            if((hc.c_plus>0)||(hc.c_minus>0)){
                std::cout << "[pi0_corr][" << period << "] bin " << keyStr(k)
                          << " contam(+)=" << hc.c_plus << " contam(-)=" << hc.c_minus
                          << " -> corrected: N+(raw="<<nd.plus<<", corr="<<br.N_corr.val_plus<<") "
                          << "N-(raw="<<nd.minus<<", corr="<<br.N_corr.val_minus<<")\n";
            }

            bins[k] = br;
        }

        allPeriods[period] = bins;

        // Write per-period file
        fs::path out_file = json_out_dir / ("pi0_corrected_counts_" + period + ".json");
        std::ofstream ofs(out_file);
        ofs << std::fixed << std::setprecision(8);
        ofs << "{\n  \"binning_meta\": {\"phi_bins\": " << N_PHI_BINS
            << ", \"xB_bins\": " << xB_bins.size()
            << ", \"Q2_bins\": " << Q2_bins.size()
            << ", \"t_bins\": " << t_bins.size() << "},\n";
        ofs << "  \"bins\": {\n";
        bool first=true;
        for(const auto& kv : bins){
            if(!first) ofs << ",\n";
            first=false;
            const BinRecord& br = kv.second;
            ofs << "    \""<<keyStr(kv.first)<<"\": {"
                << "\"N_data\":{\"helicity\":{\"+1\":"<<br.N_data.plus<<",\"-1\":"<<br.N_data.minus<<"},\"total\":"<<(br.N_data.plus+br.N_data.minus)<<"},"
                << "\"contamination\":{\"+1\":{\"value\":"<<br.contam.c_plus<<",\"err\":"<<br.contam.c_plus_err<<"},"
                << "\"-1\":{\"value\":"<<br.contam.c_minus<<",\"err\":"<<br.contam.c_minus_err<<"}},"
                << "\"N_corrected\":{\"helicity\":{\"+1\":{\"value\":"<<br.N_corr.val_plus<<",\"err\":"<<br.N_corr.err_plus<<"},"
                << "\"-1\":{\"value\":"<<br.N_corr.val_minus<<",\"err\":"<<br.N_corr.err_minus<<"}},"
                << "\"total\":{\"value\":"<<(br.N_corr.val_plus+br.N_corr.val_minus)
                << ",\"err\":"<<std::sqrt(br.N_corr.err_plus*br.N_corr.err_plus + br.N_corr.err_minus*br.N_corr.err_minus)<<"}}}";
        }
        ofs << "\n  }\n}\n";
        ofs.close();
        std::cout << "[pi0_corr] Wrote " << out_file << "\n";
    }

    // Combined file (all periods)
    fs::path comb = json_out_dir / "pi0_corrected_counts_all_periods.json";
    std::ofstream ofs(comb);
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n  \"periods\": {\n";
    bool firstP=true;
    for(const auto& pkv : allPeriods){
        if(!firstP) ofs << ",\n";
        firstP=false;
        ofs << "    \""<<pkv.first<<"\": {\n      \"bins\": {\n";
        bool firstB=true;
        for(const auto& kv : pkv.second){
            if(!firstB) ofs << ",\n";
            firstB=false;
            const BinRecord& br = kv.second;
            ofs << "        \""<<keyStr(kv.first)<<"\": {"
                << "\"N_data\":{\"+1\":"<<br.N_data.plus<<",\"-1\":"<<br.N_data.minus<<"},"
                << "\"contam\":{\"+1\":"<<br.contam.c_plus<<",\"-1\":"<<br.contam.c_minus<<"},"
                << "\"N_corr\":{\"+1\":"<<br.N_corr.val_plus<<",\"-1\":"<<br.N_corr.val_minus<<"}"
                << "}";
        }
        ofs << "\n      }\n    }";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_corr] Wrote combined " << comb << "\n";
}