// load_csv.cpp
#include "load_csv.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

namespace fs = std::filesystem;

// ---------- small helpers ----------
static inline std::string trim(std::string s) {
    auto issp = [](int ch){ return std::isspace(ch); };
    while (!s.empty() && issp((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && issp((unsigned char)s.back()))  s.pop_back();
    return s;
}

static inline double to_double(const std::string& s, double def = 0.0) {
    try {
        size_t idx = 0;
        double v = std::stod(s, &idx);
        return v;
    } catch (...) { return def; }
}

static std::vector<std::string> parse_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool inq = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') {
            if (inq && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"'); // escaped quote
                ++i;
            } else {
                inq = !inq;
            }
        } else if (c == ',' && !inq) {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(trim(cur));
    return out;
}

// normalize header: lowercase and drop spaces/punct we do not care about
static std::string norm(const std::string& s) {
    std::string r;
    r.reserve(s.size());
    for (char c : s) {
        if (std::isalnum((unsigned char)c)) r.push_back((char)std::tolower((unsigned char)c));
    }
    return r;
}

static int find_col(const std::vector<std::string>& headers_norm, const std::string& query_norm) {
    for (int i = 0; i < (int)headers_norm.size(); ++i) {
        if (headers_norm[i] == query_norm) return i;
    }
    // fallback: substring contains
    for (int i = 0; i < (int)headers_norm.size(); ++i) {
        if (headers_norm[i].find(query_norm) != std::string::npos) return i;
    }
    return -1;
}

static double getd(const std::vector<std::string>& row, int idx, double def = 0.0) {
    if (idx < 0 || idx >= (int)row.size()) return def;
    return to_double(row[idx], def);
}

std::string make_bin_key(double xBmin,double xBmax,
                         double Q2min,double Q2max,
                         double tmin,double tmax,
                         double phimin,double phimax) {
    std::ostringstream os;
    os.setf(std::ios::fixed);
    os << std::setprecision(6)
       << xBmin << "," << xBmax << ","
       << Q2min << "," << Q2max << ","
       << tmin  << "," << tmax  << ","
       << phimin << "," << phimax;
    return os.str();
}

// read a CSV file fully (headers + rows)
struct Table {
    std::vector<std::string> headers;
    std::vector<std::string> headers_norm;
    std::vector<std::vector<std::string>> rows;
};

static Table read_csv(const std::string& path) {
    Table t;
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[lee][warn] Could not open " << path << "\n";
        return t;
    }
    std::string line;
    bool first = true;
    while (std::getline(f, line)) {
        // skip empty lines
        bool allspace = true;
        for (char c : line) if (!std::isspace((unsigned char)c)) { allspace = false; break; }
        if (allspace) continue;

        auto cols = parse_csv_line(line);
        if (first) {
            t.headers = cols;
            t.headers_norm.reserve(cols.size());
            for (const auto& h : cols) t.headers_norm.push_back(norm(h));
            first = false;
        } else {
            t.rows.push_back(std::move(cols));
        }
    }
    return t;
}

// map of normalized names we will query
struct ColMap {
    int xBmin=-1, xBmax=-1, Q2min=-1, Q2max=-1, tmin=-1, tmax=-1, phimin=-1, phimax=-1, phiavg=-1;

    int raw_inb_fd_fd=-1, raw_inb_cd_fd=-1, raw_inb_cd_ft=-1;
    int raw_out_fd_fd=-1, raw_out_cd_fd=-1, raw_out_cd_ft=-1;

    int contam_inb=-1, contam_out=-1;
    int signal_inb=-1, signal_out=-1;

    int acc_corr_yield_combined=-1;

    int Frad=-1, Fbin=-1, bin_volume=-1;
    int xsec=-1;
};

// try to bind columns by likely names; we use contains fallback if exact not found
static ColMap bind_all_columns(const std::vector<std::string>& HN) {
    ColMap m;
    // edges
    m.xBmin  = find_col(HN, "xbmin");
    m.xBmax  = find_col(HN, "xbmax");
    m.Q2min  = find_col(HN, "q2min");
    m.Q2max  = find_col(HN, "q2max");
    m.tmin   = find_col(HN, "tabsmin");
    if (m.tmin < 0) m.tmin = find_col(HN, "tabsmingev2"); // just in case
    m.tmax   = find_col(HN, "tabsmax");
    if (m.tmax < 0) m.tmax = find_col(HN, "tabsmaxgev2");
    m.phimin = find_col(HN, "phimin");
    m.phimax = find_col(HN, "phimax");
    m.phiavg = find_col(HN, "phiavg");

    // raw yields (inb)
    m.raw_inb_fd_fd = find_col(HN, "rawyieldepepgfdfdexpinbending");
    m.raw_inb_cd_fd = find_col(HN, "rawyieldepepgcdfdexpinbending");
    m.raw_inb_cd_ft = find_col(HN, "rawyieldepepgcdftexpinbending");
    // raw yields (out)
    m.raw_out_fd_fd = find_col(HN, "rawyieldepepgfdfdexpoutbending");
    m.raw_out_cd_fd = find_col(HN, "rawyieldepepgcdfdexpoutbending");
    m.raw_out_cd_ft = find_col(HN, "rawyieldepepgcdftexpoutbending");

    // contamination and signal
    m.contam_inb    = find_col(HN, "contaminationratioinbending");
    m.contam_out    = find_col(HN, "contaminationratiooutbending");
    m.signal_inb    = find_col(HN, "signalyieldepepgexpinbending");
    m.signal_out    = find_col(HN, "signalyieldepepgexpoutbending");

    // acceptance corrected combined
    m.acc_corr_yield_combined = find_col(HN, "acceptancecorrectedyieldepepgexp");

    // systematics
    m.Frad       = find_col(HN, "frad");
    m.Fbin       = find_col(HN, "fbin");
    m.bin_volume = find_col(HN, "binvolume");

    // optional cross section
    m.xsec       = find_col(HN, "crosssectionsepepgexp");
    if (m.xsec < 0) m.xsec = find_col(HN, "crosssectionsepepg");
    return m;
}

static void bind_acc_columns(const std::vector<std::string>& HN,
                             int& acc_inb, int& acc_out) {
    acc_inb = find_col(HN, "acceptanceinb");
    acc_out = find_col(HN, "acceptanceoutb");
    if (acc_out < 0) acc_out = find_col(HN, "acceptanceout"); // tolerate slight naming diffs
}

// pretty print one row
static void print_row(const LeeRow& r, size_t idx) {
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(6);
    std::cout
        << "[lee] row " << idx
        << " xB[" << r.xBmin << "," << r.xBmax << "]"
        << " Q2[" << r.Q2min << "," << r.Q2max << "]"
        << " t["  << r.t_abs_min << "," << r.t_abs_max << "]"
        << " phi["<< r.phimin << "," << r.phimax << "]"
        << " phiavg=" << r.phiavg
        << " | raw_inb_sum=" << r.raw_inb_sum()
        << " raw_out_sum=" << r.raw_out_sum()
        << " raw_combined=" << r.raw_combined()
        << " | contam(inb,out)=(" << r.contam_ratio_inb << "," << r.contam_ratio_out << ")"
        << " | signal(inb,out)=(" << r.signal_yield_inb << "," << r.signal_yield_out << ")"
        << " | acc(inb,out)=(" << r.acceptance_inb << "," << r.acceptance_out << ")"
        << " | acc_corr_yld=" << r.acc_corr_yield_combined
        << " | Frad=" << r.Frad << " Fbin=" << r.Fbin << " Vbin=" << r.bin_volume
        << "\n";
}

LeeData load_lee_csvs(const std::string& all_csv,
                      const std::string& full_acc_csv,
                      bool verbose) {
    LeeData out;

    // -------- all_bin_v3.csv (main table) --------
    Table A = read_csv(all_csv);
    if (A.headers.empty()) {
        std::cerr << "[lee][err] No headers in " << all_csv << "\n";
        return out;
    }
    ColMap M = bind_all_columns(A.headers_norm);

    for (size_t i = 0; i < A.rows.size(); ++i) {
        const auto& row = A.rows[i];

        LeeRow r;
        r.xBmin = getd(row, M.xBmin);
        r.xBmax = getd(row, M.xBmax);
        r.Q2min = getd(row, M.Q2min);
        r.Q2max = getd(row, M.Q2max);
        r.t_abs_min = getd(row, M.tmin);
        r.t_abs_max = getd(row, M.tmax);
        r.phimin = getd(row, M.phimin);
        r.phimax = getd(row, M.phimax);
        r.phiavg = getd(row, M.phiavg);

        r.raw_inb_fd_fd = getd(row, M.raw_inb_fd_fd);
        r.raw_inb_cd_fd = getd(row, M.raw_inb_cd_fd);
        r.raw_inb_cd_ft = getd(row, M.raw_inb_cd_ft);
        r.raw_out_fd_fd = getd(row, M.raw_out_fd_fd);
        r.raw_out_cd_fd = getd(row, M.raw_out_cd_fd);
        r.raw_out_cd_ft = getd(row, M.raw_out_cd_ft);

        r.contam_ratio_inb = getd(row, M.contam_inb);
        r.contam_ratio_out = getd(row, M.contam_out);
        r.signal_yield_inb = getd(row, M.signal_inb);
        r.signal_yield_out = getd(row, M.signal_out);

        r.acc_corr_yield_combined = getd(row, M.acc_corr_yield_combined);

        r.Frad = getd(row, M.Frad);
        r.Fbin = getd(row, M.Fbin);
        r.bin_volume = getd(row, M.bin_volume);

        r.xsec = getd(row, M.xsec);

        const std::string key = make_bin_key(r.xBmin, r.xBmax, r.Q2min, r.Q2max,
                                             r.t_abs_min, r.t_abs_max, r.phimin, r.phimax);
        out.index_by_key[key] = out.rows.size();
        out.rows.push_back(r);
    }

    // -------- full_acc.csv (acceptance inb/out only) --------
    Table B = read_csv(full_acc_csv);
    if (!B.headers.empty()) {
        int acc_inb_col = -1, acc_out_col = -1;
        bind_acc_columns(B.headers_norm, acc_inb_col, acc_out_col);

        // bind edges in this file too
        ColMap MB = bind_all_columns(B.headers_norm);

        size_t matched = 0, missed = 0;
        for (size_t i = 0; i < B.rows.size(); ++i) {
            const auto& row = B.rows[i];

            double xBmin = getd(row, MB.xBmin);
            double xBmax = getd(row, MB.xBmax);
            double Q2min = getd(row, MB.Q2min);
            double Q2max = getd(row, MB.Q2max);
            double tmin  = getd(row, MB.tmin);
            double tmax  = getd(row, MB.tmax);
            double pmin  = getd(row, MB.phimin);
            double pmax  = getd(row, MB.phimax);

            const std::string key = make_bin_key(xBmin,xBmax,Q2min,Q2max,tmin,tmax,pmin,pmax);
            auto it = out.index_by_key.find(key);
            if (it == out.index_by_key.end()) {
                ++missed;
                continue;
            }
            LeeRow& r = out.rows[it->second];
            if (acc_inb_col >= 0) r.acceptance_inb = getd(row, acc_inb_col);
            if (acc_out_col >= 0) r.acceptance_out = getd(row, acc_out_col);
            ++matched;
        }
        std::cout << "[lee] acceptance merge: matched=" << matched << " missed=" << missed << "\n";
    } else {
        std::cerr << "[lee][warn] No headers in " << full_acc_csv << " (skipping acceptance merge)\n";
    }

    // -------- verbose dump --------
    if (verbose) {
        for (size_t i = 0; i < out.rows.size(); ++i) {
            print_row(out.rows[i], i);
        }
    }
    return out;
}