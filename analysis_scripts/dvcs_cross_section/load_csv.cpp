// load_csv.cpp
#include "load_csv.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

// trim helpers
static inline std::string ltrim(std::string s) {
    size_t i=0; while (i<s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    return s.substr(i);
}
static inline std::string rtrim(std::string s) {
    if (s.empty()) return s;
    size_t i=s.size(); while (i>0 && std::isspace(static_cast<unsigned char>(s[i-1]))) --i;
    return s.substr(0,i);
}
static inline std::string trim(std::string s) { return rtrim(ltrim(std::move(s))); }

// very simple CSV splitter (assumes no embedded commas inside quoted fields)
static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (char ch : line) {
        if (ch == '"') { in_quotes = !in_quotes; continue; }
        if (ch == ',' && !in_quotes) {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(ch);
        }
    }
    out.push_back(trim(cur));
    return out;
}

// read whole CSV into rows of strings; return header index map
static bool read_csv(const std::string& path,
                     std::vector<std::vector<std::string>>& rows,
                     std::unordered_map<std::string, size_t>& col,
                     bool verbose_open=true) {
    std::ifstream in(path);
    if (!in.is_open()) {
        std::cerr << "[lee][warn] Failed to open " << path << "\n";
        return false;
    }
    if (verbose_open) std::cout << "[lee] Opened " << path << "\n";
    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "[lee][warn] Empty CSV: " << path << "\n";
        return false;
    }
    auto header = split_csv_line(line);
    for (size_t i=0;i<header.size();++i) col[header[i]] = i;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        rows.push_back(split_csv_line(line));
    }
    return true;
}

// get numeric from row by column name; returns def if missing/empty
static double getd(const std::vector<std::string>& row,
                   const std::unordered_map<std::string,size_t>& col,
                   const std::string& name, double def=0.0) {
    auto it = col.find(name);
    if (it == col.end()) return def;
    size_t idx = it->second;
    if (idx >= row.size()) return def;
    const std::string& s = row[idx];
    if (s.empty()) return def;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return def;
    return v;
}

// get integer-like field (0/1) for valid bin
static int geti(const std::vector<std::string>& row,
                const std::unordered_map<std::string,size_t>& col,
                const std::string& name, int def=0) {
    auto it = col.find(name);
    if (it == col.end()) return def;
    size_t idx = it->second;
    if (idx >= row.size()) return def;
    const std::string& s = row[idx];
    if (s.empty()) return def;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return def;
    return static_cast<int>(v);
}

// key builder with rounding (3 decimals) to be robust to small float diffs
static std::string key_from_bin(double xBmin, double xBmax,
                                double Q2min, double Q2max,
                                double tmin,  double tmax,
                                double phiavg) {
    auto fmt = [](double v)->std::string {
        std::ostringstream os; os.setf(std::ios::fixed); os<<std::setprecision(3)<<v; return os.str();
    };
    std::ostringstream k;
    k << "xb:"  << fmt(xBmin) << "-" << fmt(xBmax)
      << "|Q2:" << fmt(Q2min) << "-" << fmt(Q2max)
      << "|t:"  << fmt(tmin)  << "-" << fmt(tmax)
      << "|phi:"<< fmt(phiavg);
    return k.str();
}

// tolerant equality check for a bin (fallback if hash did not match)
static bool same_bin_tol(double a1,double a2,double eps) { return std::fabs(a1-a2) <= eps; }
static bool bin_matches(const LeeRow& a,
                        double xbmin,double xbmax,
                        double q2min,double q2max,
                        double tmin, double tmax,
                        double phia,
                        double eps_lo=1e-4, double eps_hi=5e-3) {
    // try tight first, then looser
    auto ok = [&](double e)->bool {
        return same_bin_tol(a.xBmin,xbmin,e) && same_bin_tol(a.xBmax,xbmax,e) &&
               same_bin_tol(a.Q2min,q2min,e) && same_bin_tol(a.Q2max,q2max,e) &&
               same_bin_tol(a.tmin,tmin,e)   && same_bin_tol(a.tmax,tmax,e)   &&
               same_bin_tol(a.phiavg,phia,e);
    };
    return ok(eps_lo) || ok(eps_hi);
}

struct AccRow {
    double xBmin=0, xBmax=0;
    double Q2min=0, Q2max=0;
    double tmin=0, tmax=0;
    double phimin=0, phimax=0, phiavg=0;
    double acc_inb=0, acc_out=0;
};

// load all acceptance rows + build hash index
static void load_acceptance(const std::string& acc_csv_path,
                            std::vector<AccRow>& acc_rows,
                            std::unordered_multimap<std::string,size_t>& index) {
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string,size_t> col;
    if (!read_csv(acc_csv_path, rows, col)) return;

    for (size_t i=0;i<rows.size();++i) {
        const auto& r = rows[i];
        AccRow a;
        a.xBmin  = getd(r,col,"xBmin");
        a.xBmax  = getd(r,col,"xBmax");
        a.Q2min  = getd(r,col,"Q2min");
        a.Q2max  = getd(r,col,"Q2max");
        a.tmin   = getd(r,col,"t_abs_min");
        a.tmax   = getd(r,col,"t_abs_max");
        a.phimin = getd(r,col,"phimin");
        a.phimax = getd(r,col,"phimax");
        a.phiavg = getd(r,col,"phiavg");
        a.acc_inb= getd(r,col,"acceptance_inb");
        a.acc_out= getd(r,col,"acceptance_out");

        size_t idx = acc_rows.size();
        acc_rows.push_back(a);
        index.emplace(key_from_bin(a.xBmin,a.xBmax,a.Q2min,a.Q2max,a.tmin,a.tmax,a.phiavg), idx);
    }

    std::cout << "[lee] full_acc rows loaded: " << acc_rows.size() << "\n";
}

} // anon

LeeData load_lee_csvs(const std::string& all_bin_csv,
                      const std::string& full_acc_csv,
                      bool verbose) {
    LeeData out;

    // 1) Read acceptance CSV and build index
    std::vector<AccRow> acc_rows;
    std::unordered_multimap<std::string,size_t> acc_index;
    load_acceptance(full_acc_csv, acc_rows, acc_index);

    // 2) Read all_bin_v3
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string,size_t> col;
    if (!read_csv(all_bin_csv, rows, col)) return out;

    size_t kept = 0, skipped_invalid = 0, matched_acc = 0, unmatched_acc = 0;

    for (size_t i=0;i<rows.size();++i) {
        const auto& r = rows[i];

        // filter on "valid bin"
        int valid = geti(r, col, "valid bin", 0);
        if (valid != 1) { ++skipped_invalid; continue; }

        LeeRow lr;
        lr.valid  = true;
        lr.xBmin  = getd(r,col,"xBmin");
        lr.xBmax  = getd(r,col,"xBmax");
        lr.Q2min  = getd(r,col,"Q2min");
        lr.Q2max  = getd(r,col,"Q2max");
        lr.tmin   = getd(r,col,"t_abs_min");
        lr.tmax   = getd(r,col,"t_abs_max");
        lr.phimin = getd(r,col,"phimin");
        lr.phimax = getd(r,col,"phimax");
        lr.phiavg = getd(r,col,"phiavg");

        // raw yields (sum detector topologies) - inbending
        lr.raw_inb_fd_fd = getd(r,col,"raw yield, ep->epg, (FD, FD), exp, inbending");
        lr.raw_inb_cd_fd = getd(r,col,"raw yield, ep->epg, (CD, FD), exp, inbending");
        lr.raw_inb_cd_ft = getd(r,col,"raw yield, ep->epg, (CD, FT), exp, inbending");
        // raw yields - outbending
        lr.raw_out_fd_fd = getd(r,col,"raw yield, ep->epg, (FD, FD), exp, outbending");
        lr.raw_out_cd_fd = getd(r,col,"raw yield, ep->epg, (CD, FD), exp, outbending");
        lr.raw_out_cd_ft = getd(r,col,"raw yield, ep->epg, (CD, FT), exp, outbending");

        lr.raw_inb_sum = lr.raw_inb_fd_fd + lr.raw_inb_cd_fd + lr.raw_inb_cd_ft;
        lr.raw_out_sum = lr.raw_out_fd_fd + lr.raw_out_cd_fd + lr.raw_out_cd_ft;
        lr.raw_combined = lr.raw_inb_sum + lr.raw_out_sum;

        // contamination and signal yields
        lr.contam_inb = getd(r,col,"contamination ratio, inbending");
        lr.contam_out = getd(r,col,"contamination ratio, outbending");
        lr.signal_inb = getd(r,col,"signal yield, ep->epg, exp, inbending");
        lr.signal_out = getd(r,col,"signal yield, ep->epg, exp, outbending");

        // acceptance corrected yield (combined)
        lr.acc_corr_yld = getd(r,col,"acceptance corrected yield, ep->epg, exp");

        // systematics
        lr.Frad = getd(r,col,"Frad",1.0);
        lr.Fbin = getd(r,col,"Fbin",1.0);
        lr.Vbin = getd(r,col,"bin_volume",1.0);

        // --- match acceptance from full_acc.csv ---
        lr.acc_inb = 0.0;
        lr.acc_out = 0.0;

        const std::string k = key_from_bin(lr.xBmin,lr.xBmax,lr.Q2min,lr.Q2max,lr.tmin,lr.tmax,lr.phiavg);
        auto range = acc_index.equal_range(k);

        bool found = false;
        size_t best_idx = std::numeric_limits<size_t>::max();
        double best_dphi = 1e9;

        // 1) First try hash matches (there can be multiple if duplicates)
        for (auto it = range.first; it != range.second; ++it) {
            const AccRow& a = acc_rows[it->second];
            double dphi = std::fabs(a.phiavg - lr.phiavg);
            if (dphi < best_dphi) { best_dphi = dphi; best_idx = it->second; found = true; }
        }

        // 2) Fallback: tolerant scan
        if (!found) {
            for (size_t j=0;j<acc_rows.size();++j) {
                const AccRow& a = acc_rows[j];
                if (bin_matches(lr, a.xBmin,a.xBmax,a.Q2min,a.Q2max,a.tmin,a.tmax,a.phiavg)) {
                    double dphi = std::fabs(a.phiavg - lr.phiavg);
                    if (dphi < best_dphi) { best_dphi = dphi; best_idx = j; found = true; }
                }
            }
        }

        if (found) {
            lr.acc_inb = acc_rows[best_idx].acc_inb;
            lr.acc_out = acc_rows[best_idx].acc_out;
            ++matched_acc;
        } else {
            ++unmatched_acc;
        }

        out.rows.push_back(lr);
        if (verbose) {
            const LeeRow& t = out.rows.back();
            std::cout << "[lee] row " << kept
                      << " xB[" << t.xBmin << "," << t.xBmax << "]"
                      << " Q2[" << t.Q2min << "," << t.Q2max << "]"
                      << " t["  << t.tmin  << "," << t.tmax  << "]"
                      << " phi["<< t.phimin<< "," << t.phimax<< "]"
                      << " phiavg=" << t.phiavg
                      << " | raw_inb_sum=" << t.raw_inb_sum
                      << " raw_out_sum=" << t.raw_out_sum
                      << " raw_combined=" << t.raw_combined
                      << " | contam(inb,out)=(" << t.contam_inb << "," << t.contam_out << ")"
                      << " | signal(inb,out)=(" << t.signal_inb << "," << t.signal_out << ")"
                      << " | acc(inb,out)=(" << t.acc_inb << "," << t.acc_out << ")"
                      << " | acc_corr_yld=" << t.acc_corr_yld
                      << " | Frad=" << t.Frad << " Fbin=" << t.Fbin << " Vbin=" << t.Vbin
                      << "\n";
        }
        ++kept;
    }

    std::cout << "[lee] all_bin_v3 kept rows (valid==1): " << kept
              << " ; skipped invalid: " << skipped_invalid << "\n";
    std::cout << "[lee] acceptance matched: " << matched_acc
              << " ; unmatched: " << unmatched_acc << "\n";
    return out;
}