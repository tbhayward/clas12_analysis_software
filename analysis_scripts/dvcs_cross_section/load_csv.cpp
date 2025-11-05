#include "load_csv.h"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Simple logger
static inline void lee_info(const std::string& s)  { std::cout << "[lee] " << s << std::endl; }
static inline void lee_warn(const std::string& s)  { std::cout << "[lee][warn] " << s << std::endl; }

// Quantize doubles to int64 keys at 1e-6 precision to avoid float mismatches.
static inline int64_t q6(double v) {
    return (int64_t) llround(v * 1000000.0);
}

// Composite key for acceptance lookup
struct AccKey {
    int64_t xBmin, xBmax;
    int64_t Q2min, Q2max;
    int64_t tmin,  tmax;
    int64_t phiavg;

    bool operator==(const AccKey& o) const {
        return xBmin==o.xBmin && xBmax==o.xBmax &&
               Q2min==o.Q2min && Q2max==o.Q2max &&
               tmin==o.tmin && tmax==o.tmax &&
               phiavg==o.phiavg;
    }
};

struct AccKeyHash {
    size_t operator()(const AccKey& k) const {
        // Cheap hash mixer
        uint64_t h = 1469598103934665603ull;
        auto mix = [&](int64_t v){
            h ^= (uint64_t) v;
            h *= 1099511628211ull;
        };
        mix(k.xBmin); mix(k.xBmax);
        mix(k.Q2min); mix(k.Q2max);
        mix(k.tmin);  mix(k.tmax);
        mix(k.phiavg);
        return (size_t) h;
    }
};

// Split CSV line into fields (assumes numeric columns have no embedded commas in quotes)
static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur; cur.reserve(line.size());
    bool in_quotes = false;
    for (char c : line) {
        if (c == '"') { in_quotes = !in_quotes; continue; }
        if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

static inline double to_double(const std::string& s) {
    if (s.empty()) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static inline int to_int(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static std::unordered_map<std::string,int> header_index(const std::vector<std::string>& cols) {
    std::unordered_map<std::string,int> m;
    for (int i = 0; i < (int)cols.size(); ++i) m[cols[i]] = i;
    return m;
}

// Helpers: safe col fetch
static inline double getd(const std::vector<std::string>& row, const std::unordered_map<std::string,int>& idx, const char* key) {
    auto it = idx.find(key);
    if (it == idx.end()) return 0.0;
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return 0.0;
    return to_double(row[j]);
}
static inline int geti(const std::vector<std::string>& row, const std::unordered_map<std::string,int>& idx, const char* key) {
    auto it = idx.find(key);
    if (it == idx.end()) return 0;
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return 0;
    return to_int(row[j]);
}

// Build acceptance map from full_acc.csv
static std::unordered_map<AccKey, std::pair<double,double>, AccKeyHash>
build_acceptance_map(const std::string& full_acc_path)
{
    std::unordered_map<AccKey, std::pair<double,double>, AccKeyHash> accmap;

    std::ifstream f(full_acc_path);
    if (!f.is_open()) {
        lee_warn(std::string("Cannot open ") + full_acc_path);
        return accmap;
    }

    std::string header;
    if (!std::getline(f, header)) {
        lee_warn(std::string("Empty file: ") + full_acc_path);
        return accmap;
    }
    auto hcols = split_csv_line(header);
    auto H = header_index(hcols);

    // Required column names in full_acc.csv:
    // xBmin,xBmax,Q2min,Q2max,t1min,t1max,phi_avg_this_point,acceptance_inb,acceptance_outb
    int rown = 0, ok = 0;
    std::string line;
    while (std::getline(f, line)) {
        ++rown;
        if (line.empty()) continue;
        auto cols = split_csv_line(line);

        double xBmin = getd(cols, H, "xBmin");
        double xBmax = getd(cols, H, "xBmax");
        double Q2min = getd(cols, H, "Q2min");
        double Q2max = getd(cols, H, "Q2max");
        double tmin  = getd(cols, H, "t1min");
        double tmax  = getd(cols, H, "t1max");
        double phiav = getd(cols, H, "phi_avg_this_point");

        double ainb  = getd(cols, H, "acceptance_inb");
        double aout  = getd(cols, H, "acceptance_outb");

        AccKey k { q6(xBmin), q6(xBmax), q6(Q2min), q6(Q2max), q6(tmin), q6(tmax), q6(phiav) };
        accmap.emplace(k, std::make_pair(ainb, aout));
        ++ok;
    }
    lee_info("full_acc rows read: " + std::to_string(rown) + " ; stored: " + std::to_string(ok));
    return accmap;
}

std::vector<LeeRow> load_lee_csvs(const std::string& all_bin_v3_path,
                                  const std::string& full_acc_path,
                                  int& matched_acc,
                                  int& unmatched_acc)
{
    matched_acc = 0; unmatched_acc = 0;

    // Build acceptance lookup first
    auto accmap = build_acceptance_map(full_acc_path);

    std::vector<LeeRow> out;

    std::ifstream f(all_bin_v3_path);
    if (!f.is_open()) {
        lee_warn(std::string("Cannot open ") + all_bin_v3_path);
        return out;
    }

    std::string header;
    if (!std::getline(f, header)) {
        lee_warn(std::string("Empty file: ") + all_bin_v3_path);
        return out;
    }
    auto hcols = split_csv_line(header);
    auto H = header_index(hcols);

    // Column names we rely on in all_bin_v3.csv (checked via /mnt/data inspection):
    // Binning
    //  xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phimin,phimax,phiavg,"valid bin"
    // Raw yields
    //  raw yield, ep->epg, (FD, FD), exp, inbending
    //  raw yield, ep->epg, (CD, FD), exp, inbending
    //  raw yield, ep->epg, (CD, FT), exp, inbending
    //  raw yield, ep->epg, (FD, FD), exp, outbending
    //  raw yield, ep->epg, (CD, FD), exp, outbending
    //  raw yield, ep->epg, (CD, FT), exp, outbending
    //
    // Contamination and signal
    //  contamination ratio, inbending
    //  contamination ratio, outbending
    //  signal yield, ep->epg, exp, inbending
    //  signal yield, ep->epg, exp, outbending
    //
    // Acceptance-corrected yield
    //  acceptance corrected yield, ep->epg, exp
    //
    // Factors
    //  Frad, Fbin, bin_volume

    int kept = 0, skipped = 0, rown = 0;
    std::string line;
    while (std::getline(f, line)) {
        ++rown;
        if (line.empty()) continue;
        auto cols = split_csv_line(line);

        int valid = geti(cols, H, "valid bin");
        if (valid != 1) { ++skipped; continue; }

        LeeRow r{};

        r.xBmin  = getd(cols, H, "xBmin");
        r.xBmax  = getd(cols, H, "xBmax");
        r.Q2min  = getd(cols, H, "Q2min");
        r.Q2max  = getd(cols, H, "Q2max");
        r.tmin   = getd(cols, H, "t_abs_min");
        r.tmax   = getd(cols, H, "t_abs_max");
        r.phimin = getd(cols, H, "phimin");
        r.phimax = getd(cols, H, "phimax");
        r.phiavg = getd(cols, H, "phiavg");

        // Raw yields: sum per topology for inb/out
        double inb_fd_fd = getd(cols, H, "raw yield, ep->epg, (FD, FD), exp, inbending");
        double inb_cd_fd = getd(cols, H, "raw yield, ep->epg, (CD, FD), exp, inbending");
        double inb_cd_ft = getd(cols, H, "raw yield, ep->epg, (CD, FT), exp, inbending");
        double out_fd_fd = getd(cols, H, "raw yield, ep->epg, (FD, FD), exp, outbending");
        double out_cd_fd = getd(cols, H, "raw yield, ep->epg, (CD, FD), exp, outbending");
        double out_cd_ft = getd(cols, H, "raw yield, ep->epg, (CD, FT), exp, outbending");

        r.raw_inb_sum  = inb_fd_fd + inb_cd_fd + inb_cd_ft;
        r.raw_out_sum  = out_fd_fd + out_cd_fd + out_cd_ft;
        r.raw_combined = r.raw_inb_sum + r.raw_out_sum;

        r.contam_inb   = getd(cols, H, "contamination ratio, inbending");
        r.contam_out   = getd(cols, H, "contamination ratio, outbending");
        r.signal_inb   = getd(cols, H, "signal yield, ep->epg, exp, inbending");
        r.signal_out   = getd(cols, H, "signal yield, ep->epg, exp, outbending");

        // Acceptance from full_acc via keyed lookup
        AccKey k{ q6(r.xBmin), q6(r.xBmax), q6(r.Q2min), q6(r.Q2max),
                  q6(r.tmin),  q6(r.tmax),  q6(r.phiavg) };

        auto it = accmap.find(k);
        if (it != accmap.end()) {
            r.acc_inb = it->second.first;
            r.acc_out = it->second.second;
            r.has_acceptance = true;
            ++matched_acc;
        } else {
            r.acc_inb = 0.0;
            r.acc_out = 0.0;
            r.has_acceptance = false;
            ++unmatched_acc;
        }

        r.acc_corr_yield = getd(cols, H, "acceptance corrected yield, ep->epg, exp");

        r.Frad = getd(cols, H, "Frad");
        r.Fbin = getd(cols, H, "Fbin");
        r.Vbin = getd(cols, H, "bin_volume");

        out.push_back(r);
        ++kept;

        // Row-by-row print (matches your current style)
        std::printf("[lee] row %d xB[%.3f,%.3f] Q2[%.3f,%.3f] t[%.3f,%.3f] phi[%.0f,%.0f] phiavg=%.3f | raw_inb_sum=%.0f raw_out_sum=%.0f raw_combined=%.0f | contam(inb,out)=(%.6f,%.6f) | signal(inb,out)=(%.6f,%.6f) | acc(inb,out)=(%.6f,%.6f) | acc_corr_yld=%.6f | Frad=%.5f Fbin=%.5f Vbin=%.6g\n",
                    kept,
                    r.xBmin, r.xBmax, r.Q2min, r.Q2max,
                    r.tmin,  r.tmax, r.phimin, r.phimax, r.phiavg,
                    r.raw_inb_sum, r.raw_out_sum, r.raw_combined,
                    r.contam_inb, r.contam_out,
                    r.signal_inb, r.signal_out,
                    r.acc_inb, r.acc_out,
                    r.acc_corr_yield,
                    r.Frad, r.Fbin, r.Vbin);
    }

    lee_info("all_bin_v3 kept rows (valid==1): " + std::to_string(kept) + " ; skipped invalid: " + std::to_string(skipped));
    lee_info("acceptance matched: " + std::to_string(matched_acc) + " ; unmatched: " + std::to_string(unmatched_acc));
    return out;
}