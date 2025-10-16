// pi0_corrected_counts.cpp — π0-corrected counts (helicity-resolved)

#include "pi0_corrected_counts.h"
#include "load_binning_scheme.h"  // Binning

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

// ---------------- bin helpers (keep consistent with contamination stage) ----------------
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)

static constexpr int    N_PHI_BINS = 12;

static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static inline std::string keyStr(const BinKey& k) {
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os; os<<"("<<a<<","<<b<<","<<c<<","<<d<<")";
    return os.str();
}

// ---------------- types for this stage ----------------
struct HelCounts { long long plus=0, minus=0; };          // raw DVCS counts
struct HelContam { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };
struct HelCorr   { double val_plus=0, err_plus=0, val_minus=0, err_minus=0; };

struct BinRecord {
    HelCounts N_data;    // from total_counts.json
    HelContam contam;    // from contamination_<period>.json
    HelCorr   N_corr;    // computed here
};

// ---------------- tiny JSON helpers ----------------
// Parse a signed/decimal/scientific number starting at or after ':' from 's' near position 'pos'.
static double parseNumberAfterColon(const std::string& s, size_t pos, double fallback = 0.0) {
    pos = s.find(':', pos);
    if (pos == std::string::npos) return fallback;
    size_t i = pos + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    auto isnum=[&](char ch){
        return (std::isdigit(static_cast<unsigned char>(ch)) || ch=='-'||ch=='+'||ch=='.'||ch=='e'||ch=='E');
    };
    while (j < s.size() && isnum(s[j])) ++j;
    try { return std::stod(s.substr(i, j - i)); } catch (...) { return fallback; }
}

static long long parseIntAfterColon(const std::string& s, size_t pos, long long fallback = 0) {
    pos = s.find(':', pos);
    if (pos == std::string::npos) return fallback;
    size_t i = pos + 1;
    while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    size_t j = i;
    while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j]=='-')) ++j;
    try { return std::stoll(s.substr(i, j - i)); } catch (...) { return fallback; }
}

// Return substring of the JSON object that corresponds to the value of the given key entry,
// when we already sit inside an object and have its body as 'obj'. Assumes `"key": { ... }` or `"key": ...`.
static bool extractObjectForKey(const std::string& obj, const std::string& key, std::string& out) {
    size_t p = obj.find("\"" + key + "\"");
    if (p == std::string::npos) return false;
    p = obj.find(':', p); if (p == std::string::npos) return false;
    // If it's an object, capture balanced braces; otherwise capture up to next comma/brace
    size_t i = p + 1;
    while (i < obj.size() && std::isspace(static_cast<unsigned char>(obj[i]))) ++i;
    if (i < obj.size() && obj[i] == '{') {
        int depth = 0; size_t j = i;
        for (; j < obj.size(); ++j) {
            if (obj[j] == '{') depth++;
            else if (obj[j] == '}') { depth--; if (depth == 0) { ++j; break; } }
        }
        if (depth != 0) return false;
        out = obj.substr(i, j - i);
        return true;
    } else {
        // primitive value; take until comma or closing brace
        size_t j = i;
        while (j < obj.size() && obj[j] != ',' && obj[j] != '}') ++j;
        out = obj.substr(i, j - i);
        return true;
    }
}

// Iterate entries of a "bins": { "<key>": { ... }, ... } object and call fn(key,obj)
template <typename Fn>
static void iterateBins(const std::string& binsObj, Fn&& fn) {
    size_t kpos = 0;
    while (true) {
        size_t keyS = binsObj.find('"', kpos);
        if (keyS == std::string::npos) break;
        size_t keyE = binsObj.find('"', keyS+1);
        if (keyE == std::string::npos) break;
        std::string key = binsObj.substr(keyS+1, keyE-keyS-1); // e.g. "(ix,iQ2,it,ip)"

        size_t objS = binsObj.find('{', keyE);
        if (objS == std::string::npos) break;
        int depth=0; size_t j=objS;
        for (; j < binsObj.size(); ++j) {
            if (binsObj[j]=='{') depth++;
            else if (binsObj[j]=='}'){ depth--; if(!depth){ ++j; break; } }
        }
        if (depth != 0) break;
        std::string obj = binsObj.substr(objS, j-objS);

        fn(key, obj);
        kpos = j;
    }
}

// ---------------- readers ----------------
// total_counts.json: we only need helicity-resolved DVCS counts per bin.
// We support either of these inside each bin object (be liberal in what we accept):
//   A) "N_data": {"helicity":{"+1": <int>, "-1": <int> }, ...}
//   B) "helicity": { "+1": <int>, "-1": <int> }
//   C) "N_plus": <int>, "N_minus": <int>   (fallback)
static std::map<BinKey, HelCounts> read_total_counts(const std::string& path) {
    std::map<BinKey, HelCounts> out;

    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[pi0_corr][ERROR] Cannot open total counts JSON: " << path << std::endl;
        return out;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) { std::cerr << "[pi0_corr][ERROR] 'bins' not found in total counts.\n"; return out; }
    pos = s.find('{', pos);
    if (pos == std::string::npos) return out;
    int depth = 0; size_t i = pos;
    for (; i < s.size(); ++i) { if (s[i]=='{') depth++; else if (s[i]=='}'){ depth--; if (!depth){ ++i; break; } } }
    std::string binsObj = s.substr(pos, i-pos);

    iterateBins(binsObj, [&](const std::string& key, const std::string& obj){
        int ix, iQ2, itb, ip;
        if (std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix, &iQ2, &itb, &ip) != 4) return;

        long long Np = 0, Nm = 0;

        std::string NdataObj;
        std::string helicityObj;

        // Prefer "N_data"->"helicity"
        if (extractObjectForKey(obj, "N_data", NdataObj)) {
            if (extractObjectForKey(NdataObj, "helicity", helicityObj)) {
                size_t p1 = helicityObj.find("\"+1\"");
                size_t p2 = helicityObj.find("\"-1\"");
                if (p1 != std::string::npos) Np = parseIntAfterColon(helicityObj, p1, 0);
                if (p2 != std::string::npos) Nm = parseIntAfterColon(helicityObj, p2, 0);
            }
        } else if (extractObjectForKey(obj, "helicity", helicityObj)) { // or "helicity" directly
            size_t p1 = helicityObj.find("\"+1\"");
            size_t p2 = helicityObj.find("\"-1\"");
            if (p1 != std::string::npos) Np = parseIntAfterColon(helicityObj, p1, 0);
            if (p2 != std::string::npos) Nm = parseIntAfterColon(helicityObj, p2, 0);
        } else {
            // fallback names
            size_t pP = obj.find("\"N_plus\"");
            size_t pM = obj.find("\"N_minus\"");
            if (pP != std::string::npos) Np = parseIntAfterColon(obj, pP, 0);
            if (pM != std::string::npos) Nm = parseIntAfterColon(obj, pM, 0);
        }

        out[BinKey(ix,iQ2,itb,ip)] = HelCounts{Np,Nm};
    });

    return out;
}

// contamination_<period>.json: we rely on the format written by writeContaminationJson()
//   "contamination": {
//        "+1": { "value": <double>, "err": <double> },
//        "-1": { "value": <double>, "err": <double> }
//   }
static std::map<BinKey, HelContam> read_contamination(const std::string& path) {
    std::map<BinKey, HelContam> out;

    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[pi0_corr][ERROR] Cannot open contamination JSON: " << path << std::endl;
        return out;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    size_t pos = s.find("\"bins\"");
    if (pos == std::string::npos) { std::cerr << "[pi0_corr][ERROR] 'bins' not found in contamination JSON.\n"; return out; }
    pos = s.find('{', pos);
    if (pos == std::string::npos) return out;
    int depth = 0; size_t i = pos;
    for (; i < s.size(); ++i) { if (s[i]=='{') depth++; else if (s[i]=='}'){ depth--; if (!depth){ ++i; break; } } }
    std::string binsObj = s.substr(pos, i-pos);

    iterateBins(binsObj, [&](const std::string& key, const std::string& obj){
        int ix, iQ2, itb, ip;
        if (std::sscanf(key.c_str(), "(%d,%d,%d,%d)", &ix, &iQ2, &itb, &ip) != 4) return;

        std::string contObj;
        if (!extractObjectForKey(obj, "contamination", contObj)) return;

        double cP=0, eP=0, cM=0, eM=0;

        size_t pP = contObj.find("\"+1\"");
        if (pP != std::string::npos) {
            // read value/err inside "+1": { ... }
            size_t brace = contObj.find('{', pP);
            if (brace != std::string::npos) {
                int d=0; size_t j=brace; for(; j<contObj.size(); ++j){ if(contObj[j]=='{') d++; else if(contObj[j]=='}'){ d--; if(!d){ ++j; break; } } }
                std::string plusObj = contObj.substr(brace, j-brace);
                size_t pv = plusObj.find("\"value\"");
                size_t pe = plusObj.find("\"err\"");
                if (pv != std::string::npos) cP = parseNumberAfterColon(plusObj, pv, 0.0);
                if (pe != std::string::npos) eP = parseNumberAfterColon(plusObj, pe, 0.0);
            }
        }
        size_t pM = contObj.find("\"-1\"");
        if (pM != std::string::npos) {
            size_t brace = contObj.find('{', pM);
            if (brace != std::string::npos) {
                int d=0; size_t j=brace; for(; j<contObj.size(); ++j){ if(contObj[j]=='{') d++; else if(contObj[j]=='}'){ d--; if(!d){ ++j; break; } } }
                std::string minusObj = contObj.substr(brace, j-brace);
                size_t pv = minusObj.find("\"value\"");
                size_t pe = minusObj.find("\"err\"");
                if (pv != std::string::npos) cM = parseNumberAfterColon(minusObj, pv, 0.0);
                if (pe != std::string::npos) eM = parseNumberAfterColon(minusObj, pe, 0.0);
            }
        }

        out[BinKey(ix,iQ2,itb,ip)] = HelContam{cP,eP,cM,eM};
    });

    return out;
}

// ---------------- correction + writer ----------------
static HelCorr correct_one_helicity(long long N, double c, double c_err) {
    // N_corr = N * (1 - c)
    const double N_d   = static_cast<double>(std::max<long long>(0, N));
    const double one_mc = std::max(0.0, 1.0 - c);
    const double val = std::max(0.0, N_d * one_mc);

    // Var(N_corr) = (1 - c)^2 * Var(N) + (N)^2 * Var(c), with Var(N)=N (Poisson)
    const double var = (one_mc*one_mc) * N_d + (N_d*N_d) * (c_err*c_err);
    const double err = std::sqrt(std::max(0.0, var));

    HelCorr out;
    out.val_plus = val;  // NOTE: caller places into the right helicity field
    out.err_plus = err;
    return out;
}

static void write_corrected_json(
    const fs::path& out_path,
    const std::map<BinKey, BinRecord>& table,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins
) {
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[pi0_corr][ERROR] Cannot open " << out_path << "\n"; return; }

    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << N_PHI_BINS
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"bins\": {\n";

    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n";
        first=false;

        const auto& br = kv.second;
        const double tot_corr = br.N_corr.val_plus + br.N_corr.val_minus;
        const double tot_err  = std::sqrt(br.N_corr.err_plus*br.N_corr.err_plus
                                        + br.N_corr.err_minus*br.N_corr.err_minus);

        ofs << "    \"" << keyStr(kv.first) << "\": {";

        // original counts
        ofs << "\"N_data\":{\"helicity\":{\"+1\":" << br.N_data.plus
            << ",\"-1\":" << br.N_data.minus
            << "},\"total\":" << (br.N_data.plus + br.N_data.minus) << "}";

        // contamination inputs
        ofs << ",\"contamination\":{"
            << "\"+1\":{\"value\":" << br.contam.c_plus  << ",\"err\":" << br.contam.c_plus_err  << "},"
            << "\"-1\":{\"value\":" << br.contam.c_minus << ",\"err\":" << br.contam.c_minus_err << "}"
            << "}";

        // corrected counts
        ofs << ",\"N_corrected\":{"
            << "\"helicity\":{"
            << "\"+1\":{\"value\":" << br.N_corr.val_plus  << ",\"err\":" << br.N_corr.err_plus  << "},"
            << "\"-1\":{\"value\":" << br.N_corr.val_minus << ",\"err\":" << br.N_corr.err_minus << "}"
            << "},"
            << "\"total\":{\"value\":" << tot_corr << ",\"err\":" << tot_err << "}"
            << "}";

        ofs << "}";
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_corr] Wrote " << out_path.string() << std::endl;
}

static void write_combined_json(
    const fs::path& out_path,
    const std::map<std::string, std::map<BinKey, BinRecord>>& byPeriod,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins
) {
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[pi0_corr][ERROR] Cannot open combined output " << out_path << "\n"; return; }

    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << N_PHI_BINS
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"periods\": {\n";

    bool firstP = true;
    for (const auto& pkv : byPeriod) {
        if (!firstP) ofs << ",\n";
        firstP = false;

        ofs << "    \"" << pkv.first << "\": {\n";
        ofs << "      \"bins\": {\n";

        bool firstB = true;
        for (const auto& kv : pkv.second) {
            if (!firstB) ofs << ",\n";
            firstB = false;

            const auto& br = kv.second;
            const double tot_corr = br.N_corr.val_plus + br.N_corr.val_minus;
            const double tot_err  = std::sqrt(br.N_corr.err_plus*br.N_corr.err_plus
                                            + br.N_corr.err_minus*br.N_corr.err_minus);

            ofs << "        \"" << keyStr(kv.first) << "\": {";
            ofs << "\"N_data\":{\"helicity\":{\"+1\":" << br.N_data.plus
                << ",\"-1\":" << br.N_data.minus << "},\"total\":" << (br.N_data.plus + br.N_data.minus) << "}";
            ofs << ",\"contamination\":{"
                << "\"+1\":{\"value\":" << br.contam.c_plus  << ",\"err\":" << br.contam.c_plus_err  << "},"
                << "\"-1\":{\"value\":" << br.contam.c_minus << ",\"err\":" << br.contam.c_minus_err << "}"
                << "}";
            ofs << ",\"N_corrected\":{"
                << "\"helicity\":{"
                << "\"+1\":{\"value\":" << br.N_corr.val_plus  << ",\"err\":" << br.N_corr.err_plus  << "},"
                << "\"-1\":{\"value\":" << br.N_corr.val_minus << ",\"err\":" << br.N_corr.err_minus << "}"
                << "},"
                << "\"total\":{\"value\":" << tot_corr << ",\"err\":" << tot_err << "}"
                << "}";
            ofs << "}";
        }

        ofs << "\n      }\n"; // bins
        ofs << "    }";
    }

    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_corr] Wrote combined " << out_path.string() << std::endl;
}

// ---------------- main entry ----------------
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>& binning_scheme,
    const std::string& total_counts_json,
    const std::string& contamination_dir,
    const std::string& out_root_dir
) {
    // binning meta (for JSON headers only)
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[pi0_corr][ERROR] Missing binning ranges." << std::endl;
        return;
    }

    if (!fs::exists(total_counts_json)) {
        std::cerr << "[pi0_corr][ERROR] total_counts.json not found: " << total_counts_json << std::endl;
        return;
    }
    if (!fs::exists(contamination_dir)) {
        std::cerr << "[pi0_corr][ERROR] contamination dir not found: " << contamination_dir << std::endl;
        return;
    }

    // Read the global DVCS counts once
    const auto total_counts_map = read_total_counts(total_counts_json);
    if (total_counts_map.empty()) {
        std::cerr << "[pi0_corr][ERROR] No bins parsed from total_counts.json – aborting." << std::endl;
        return;
    }

    // Prepare output dirs
    const fs::path out_root(out_root_dir);
    const fs::path out_json_dir = out_root / "jsons";
    std::error_code ec;
    fs::create_directories(out_json_dir, ec);

    // Keep all per-period tables to write a combined file
    std::map<std::string, std::map<BinKey, BinRecord>> allPeriods;

    // For each period: read contamination JSON and produce corrected counts
    for (const auto& period : periods) {
        const fs::path contam_path = fs::path(contamination_dir) / ("contamination_" + period + ".json");
        if (!fs::exists(contam_path)) {
            std::cerr << "[pi0_corr][WARN] Missing contamination JSON for " << period
                      << " (" << contam_path.string() << "). Skipping period." << std::endl;
            continue;
        }
        auto contam_map = read_contamination(contam_path.string());
        if (contam_map.empty()) {
            std::cerr << "[pi0_corr][WARN] No bins in contamination JSON for " << period << ". Skipping period." << std::endl;
            continue;
        }

        std::map<BinKey, BinRecord> table;

        // Build per-bin records where we have DVCS counts (primary key space)
        for (const auto& kv : total_counts_map) {
            const BinKey& key = kv.first;
            BinRecord rec;
            rec.N_data = kv.second;

            auto itC = contam_map.find(key);
            if (itC != contam_map.end())
                rec.contam = itC->second; // else leaves zeros

            // Compute helicity +1
            {
                const auto tmp = correct_one_helicity(rec.N_data.plus, rec.contam.c_plus, rec.contam.c_plus_err);
                rec.N_corr.val_plus = tmp.val_plus;
                rec.N_corr.err_plus = tmp.err_plus;
            }
            // Compute helicity -1
            {
                const auto tmp = correct_one_helicity(rec.N_data.minus, rec.contam.c_minus, rec.contam.c_minus_err);
                rec.N_corr.val_minus = tmp.val_plus; // using val_plus/err_plus fields from helper
                rec.N_corr.err_minus = tmp.err_plus;
            }

            table[key] = rec;
        }

        // Write this period
        const fs::path out_period = out_json_dir / ("pi0_corrected_counts_" + period + ".json");
        write_corrected_json(out_period, table, xB_bins, Q2_bins, t_bins);

        // Keep for combined
        allPeriods[period] = std::move(table);
    }

    // Combined file
    const fs::path out_combined = out_json_dir / "pi0_corrected_counts_all_periods.json";
    write_combined_json(out_combined, allPeriods, xB_bins, Q2_bins, t_bins);
}