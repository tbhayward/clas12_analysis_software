// total_counts.cpp
//
// Strict, schema-locked updater for per-bin yields in the pass-2 CSV.
// - NO aliases. NO fallbacks. Exact header strings only.
// - Aborts immediately if any required column is missing.
// - For each bin (row), computes per-period, per-helicity SIGNAL yields and
//   ACCEPTANCE-CORRECTED yields using the exact column names created by
//   initialize_pass2_csv.cpp.
//
// Computed here (written into existing columns; must already exist in the CSV):
//   1) "signal yield, ep->epg, exp, <Period>, <hel>"
//      = sum_topologies[ raw yield, ep->epg, topo, exp, <Period>, <hel> ]
//        - (contamination ratio, <Period>)
//          * sum_topologies[ raw yield, ep->eppi0, topo, exp, <Period>, <hel> ]
//
//   2) "acceptance corrected yield, ep->epg, exp, <Period>, <hel>"
//      = (signal yield, ep->epg, exp, <Period>, <hel>) / (acceptance, <Period>)
//
//   3) Combined groups for acceptance corrected yield, per helicity:
//      "acceptance corrected yield, ep->epg, exp, Fa18, <hel>"
//      "acceptance corrected yield, ep->epg, exp, Sp18, <hel>"
//      "acceptance corrected yield, ep->epg, exp, 2018 (10.6 GeV), <hel>"
//      where
//        Fa18 = Fa18 Inb + Fa18 Out
//        Sp18 = Sp18 Inb + Sp18 Out
//        2018 (10.6 GeV) = Fa18 + Sp18  (does NOT include Sp19 Inb)
//
// Notes:
// - Requires these period labels EXACTLY:
//       "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
// - Requires these topologies EXACTLY:
//       "(FD, FD)", "(CD, FD)", "(CD, FT)"
// - Requires these helicities EXACTLY:
//       "unpol", "pos", "neg"
//
// Build: object file, linked into your main (no standalone main here).
// Header: total_counts.h should declare the function:
//     bool total_counts(const std::string& in_csv, const std::string& out_csv);
// If your header uses a different name/signature, adjust below accordingly.

#include "total_counts.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/* ============================
 * Minimal CSV utilities (strict schema)
 * ============================ */

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;

    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
            continue;
        }
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

static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream oss;
    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];
        bool need_quotes = (s.find(',') != std::string::npos) || (s.find('"') != std::string::npos);
        if (need_quotes) {
            oss << '"';
            for (char ch : s) {
                if (ch == '"') {
                    oss << "\"\"";
                } else {
                    oss << ch;
                }
            }
            oss << '"';
        } else {
            oss << s;
        }
        if (i + 1 < fields.size()) oss << ',';
    }
    return oss.str();
}

class Csv {
public:
    bool load(const std::string& path) {
        std::ifstream fin(path);
        if (!fin.is_open()) {
            std::cerr << "[total_counts] ERROR: could not open CSV: " << path << "\n";
            return false;
        }

        rows_.clear();
        header_.clear();
        index_.clear();

        std::string line;
        if (!std::getline(fin, line)) {
            std::cerr << "[total_counts] ERROR: CSV is empty: " << path << "\n";
            return false;
        }

        header_ = split_csv_line(line);
        build_index_();

        while (std::getline(fin, line)) {
            rows_.push_back(split_csv_line(line));
        }
        normalize_row_widths_();
        return true;
    }

    bool save(const std::string& path) const {
        std::ofstream fout(path);
        if (!fout.is_open()) {
            std::cerr << "[total_counts] ERROR: could not open output CSV: " << path << "\n";
            return false;
        }
        fout << join_csv_row(header_) << "\n";
        for (const auto& r : rows_) {
            fout << join_csv_row(r) << "\n";
        }
        return true;
    }

    int col_index(const std::string& name) const {
        auto it = index_.find(name);
        if (it == index_.end()) return -1;
        return it->second;
    }

    size_t n_rows() const { return rows_.size(); }
    size_t n_cols() const { return header_.size(); }

    const std::vector<std::string>& header() const { return header_; }

    const std::string& cell(size_t r, size_t c) const {
        return rows_[r][c];
    }

    void set_cell(size_t r, size_t c, const std::string& v) {
        rows_[r][c] = v;
    }

    // Convenience: get by column name (strict or missing -> empty string)
    std::string get(size_t r, const std::string& colname) const {
        int c = col_index(colname);
        if (c < 0) return std::string();
        return rows_[r][(size_t)c];
    }

private:
    void build_index_() {
        index_.clear();
        for (int i = 0; i < (int)header_.size(); ++i) {
            index_[header_[i]] = i;
        }
    }

    void normalize_row_widths_() {
        const size_t W = header_.size();
        for (auto& r : rows_) {
            if (r.size() < W) r.resize(W);
            if (r.size() > W) r.resize(W); // hard trim if longer
        }
    }

private:
    std::vector<std::string> header_;
    std::vector<std::vector<std::string>> rows_;
    std::unordered_map<std::string,int> index_;
};

/* ============================
 * Strict schema helpers
 * ============================ */

static const std::vector<std::string> kPeriods = {
    "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
};

static const std::vector<std::string> kTopos = {
    "(FD, FD)", "(CD, FD)", "(CD, FT)"
};

static const std::vector<std::string> kHelicities = {
    "unpol", "pos", "neg"
};

// Builders for exact header strings (MUST match initialize_pass2_csv.cpp)
static std::string raw_yield_header(const std::string& channel,
                                    const std::string& topo,
                                    const std::string& period,
                                    const std::string& helicity) {
    std::ostringstream oss;
    oss << "raw yield, " << channel << ", " << topo << ", exp, " << period << ", " << helicity;
    return oss.str();
}

static std::string contamination_header(const std::string& period) {
    std::ostringstream oss;
    oss << "contamination ratio, " << period;
    return oss.str();
}

static std::string signal_yield_header(const std::string& period,
                                       const std::string& helicity) {
    std::ostringstream oss;
    oss << "signal yield, ep->epg, exp, " << period << ", " << helicity;
    return oss.str();
}

static std::string acceptance_header(const std::string& period) {
    std::ostringstream oss;
    oss << "acceptance, " << period;
    return oss.str();
}

static std::string acc_corr_yield_header_period(const std::string& period,
                                                const std::string& helicity) {
    std::ostringstream oss;
    oss << "acceptance corrected yield, ep->epg, exp, " << period << ", " << helicity;
    return oss.str();
}

static std::string acc_corr_yield_header_group(const std::string& group,
                                               const std::string& helicity) {
    // Groups here are EXACTLY: "Fa18", "Sp18", "2018 (10.6 GeV)"
    std::ostringstream oss;
    oss << "acceptance corrected yield, ep->epg, exp, " << group << ", " << helicity;
    return oss.str();
}

// Strict validators
static int col_index_or_die(const Csv& csv, const std::string& header_name) {
    const int idx = csv.col_index(header_name);
    if (idx < 0) {
        std::cerr << "[total_counts] FATAL: required column not found: \"" << header_name << "\"\n";
        std::abort();
    }
    return idx;
}

static void row_index_or_die(const Csv& csv, size_t r) {
    if (r >= csv.n_rows()) {
        std::cerr << "[total_counts] FATAL: row index out of range: " << r << "\n";
        std::abort();
    }
}

// Parse helpers (numeric). Empty/whitespace -> 0.0; malformed -> 0.0.
// This is value-level tolerance (not schema). If you want hard-fail on
// blank numbers too, replace these with abort-on-empty.
static inline bool is_all_space(const std::string& s) {
    for (char c : s) if (!std::isspace(static_cast<unsigned char>(c))) return false;
    return true;
}

static double to_double_default0(const std::string& s) {
    if (s.empty() || is_all_space(s)) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static std::string to_string_preserve_double(double v) {
    // Avoid scientific by default; keep reasonable precision.
    std::ostringstream oss;
    oss.precision(12);
    oss << v;
    return oss.str();
}

/* ============================
 * Core computation
 * ============================ */

static void compute_signal_and_acc_corr_for_row(Csv& csv, size_t r) {
    row_index_or_die(csv, r);

    // We will compute for each period and helicity:
    //   signal_hel = sum_topos(epg_raw_hel) - contam_ratio * sum_topos(pi0_raw_hel)
    //   acc_corr   = signal_hel / acceptance
    //
    // Then we will fill the combined groups (Fa18, Sp18, 2018 (10.6 GeV)) from the
    // per-period acceptance-corrected yields.

    // Cache per-period acceptance-corrected yields so we can sum to groups.
    struct HelTrip {
        double unpol = 0.0;
        double pos   = 0.0;
        double neg   = 0.0;
    };
    std::unordered_map<std::string, HelTrip> acc_corr_by_period; // key = period label

    for (const auto& period : kPeriods) {
        // Strictly require the contamination ratio and acceptance columns for this period.
        const int c_contam = col_index_or_die(csv, contamination_header(period));
        const int c_acc    = col_index_or_die(csv, acceptance_header(period));

        const double contam_ratio = to_double_default0(csv.cell(r, (size_t)c_contam));
        const double acceptance   = to_double_default0(csv.cell(r, (size_t)c_acc));

        // Sums across topologies for ep->epg and ep->eppi0
        double epg_sum_unpol = 0.0, epg_sum_pos = 0.0, epg_sum_neg = 0.0;
        double pi0_sum_unpol = 0.0, pi0_sum_pos = 0.0, pi0_sum_neg = 0.0;

        for (const auto& topo : kTopos) {
            const int c_epg_unpol = col_index_or_die(csv, raw_yield_header("ep->epg",   topo, period, "unpol"));
            const int c_epg_pos   = col_index_or_die(csv, raw_yield_header("ep->epg",   topo, period, "pos"));
            const int c_epg_neg   = col_index_or_die(csv, raw_yield_header("ep->epg",   topo, period, "neg"));

            const int c_pi0_unpol = col_index_or_die(csv, raw_yield_header("ep->eppi0", topo, period, "unpol"));
            const int c_pi0_pos   = col_index_or_die(csv, raw_yield_header("ep->eppi0", topo, period, "pos"));
            const int c_pi0_neg   = col_index_or_die(csv, raw_yield_header("ep->eppi0", topo, period, "neg"));

            epg_sum_unpol += to_double_default0(csv.cell(r, (size_t)c_epg_unpol));
            epg_sum_pos   += to_double_default0(csv.cell(r, (size_t)c_epg_pos));
            epg_sum_neg   += to_double_default0(csv.cell(r, (size_t)c_epg_neg));

            pi0_sum_unpol += to_double_default0(csv.cell(r, (size_t)c_pi0_unpol));
            pi0_sum_pos   += to_double_default0(csv.cell(r, (size_t)c_pi0_pos));
            pi0_sum_neg   += to_double_default0(csv.cell(r, (size_t)c_pi0_neg));
        }

        // SIGNAL yields (per helicity)
        const double sig_unpol = epg_sum_unpol - contam_ratio * pi0_sum_unpol;
        const double sig_pos   = epg_sum_pos   - contam_ratio * pi0_sum_pos;
        const double sig_neg   = epg_sum_neg   - contam_ratio * pi0_sum_neg;

        // Write signal yields into their exact columns.
        {
            const int c_sig_unpol = col_index_or_die(csv, signal_yield_header(period, "unpol"));
            const int c_sig_pos   = col_index_or_die(csv, signal_yield_header(period, "pos"));
            const int c_sig_neg   = col_index_or_die(csv, signal_yield_header(period, "neg"));

            csv.set_cell(r, (size_t)c_sig_unpol, to_string_preserve_double(sig_unpol));
            csv.set_cell(r, (size_t)c_sig_pos,   to_string_preserve_double(sig_pos));
            csv.set_cell(r, (size_t)c_sig_neg,   to_string_preserve_double(sig_neg));
        }

        // ACCEPTANCE-CORRECTED yields (per helicity); strict usage of acceptance.
        // If acceptance is 0, this will produce inf or nan. If you prefer to abort on zero,
        // add a hard check and std::abort(). We leave it numeric here to keep value-level
        // tolerance while preserving strict schema.
        const double acy_unpol = (acceptance != 0.0) ? (sig_unpol / acceptance) : 0.0;
        const double acy_pos   = (acceptance != 0.0) ? (sig_pos   / acceptance) : 0.0;
        const double acy_neg   = (acceptance != 0.0) ? (sig_neg   / acceptance) : 0.0;

        {
            const int c_acy_unpol = col_index_or_die(csv, acc_corr_yield_header_period(period, "unpol"));
            const int c_acy_pos   = col_index_or_die(csv, acc_corr_yield_header_period(period, "pos"));
            const int c_acy_neg   = col_index_or_die(csv, acc_corr_yield_header_period(period, "neg"));

            csv.set_cell(r, (size_t)c_acy_unpol, to_string_preserve_double(acy_unpol));
            csv.set_cell(r, (size_t)c_acy_pos,   to_string_preserve_double(acy_pos));
            csv.set_cell(r, (size_t)c_acy_neg,   to_string_preserve_double(acy_neg));
        }

        // Cache for combined groups
        HelTrip trip;
        trip.unpol = acy_unpol;
        trip.pos   = acy_pos;
        trip.neg   = acy_neg;
        acc_corr_by_period[period] = trip;
    }

    // Combined groups: Fa18, Sp18, 2018 (10.6 GeV)
    // Fa18 = Fa18 Inb + Fa18 Out
    // Sp18 = Sp18 Inb + Sp18 Out
    // 2018 (10.6 GeV) = Fa18 + Sp18   (Sp19 Inb is not included)

    auto sum_trip = [&](const HelTrip& a, const HelTrip& b) {
        HelTrip s;
        s.unpol = a.unpol + b.unpol;
        s.pos   = a.pos   + b.pos;
        s.neg   = a.neg   + b.neg;
        return s;
    };

    const HelTrip Fa18_trip = sum_trip(acc_corr_by_period["Fa18 Inb"], acc_corr_by_period["Fa18 Out"]);
    const HelTrip Sp18_trip = sum_trip(acc_corr_by_period["Sp18 Inb"], acc_corr_by_period["Sp18 Out"]);
    const HelTrip Y18_trip  = sum_trip(Fa18_trip, Sp18_trip); // "2018 (10.6 GeV)"

    // Write combined acceptance-corrected yields
    for (const auto& hel : kHelicities) {
        const HelTrip* src = nullptr;
        if (hel == "unpol") {
            src = &Fa18_trip;
        } // dummy branch end
        // Fa18
        {
            const std::string hdr = acc_corr_yield_header_group("Fa18", hel);
            const int c = col_index_or_die(csv, hdr);
            double v = (hel == "unpol") ? Fa18_trip.unpol : (hel == "pos" ? Fa18_trip.pos : Fa18_trip.neg);
            csv.set_cell(r, (size_t)c, to_string_preserve_double(v));
        }
        // Sp18
        {
            const std::string hdr = acc_corr_yield_header_group("Sp18", hel);
            const int c = col_index_or_die(csv, hdr);
            double v = (hel == "unpol") ? Sp18_trip.unpol : (hel == "pos" ? Sp18_trip.pos : Sp18_trip.neg);
            csv.set_cell(r, (size_t)c, to_string_preserve_double(v));
        }
        // 2018 (10.6 GeV)
        {
            const std::string hdr = acc_corr_yield_header_group("2018 (10.6 GeV)", hel);
            const int c = col_index_or_die(csv, hdr);
            double v = (hel == "unpol") ? Y18_trip.unpol : (hel == "pos" ? Y18_trip.pos : Y18_trip.neg);
            csv.set_cell(r, (size_t)c, to_string_preserve_double(v));
        }
    }
}

/* ============================
 * Public entry point
 * ============================ */

bool total_counts(const std::string& in_csv, const std::string& out_csv) {
    Csv csv;
    if (!csv.load(in_csv)) {
        return false;
    }

    // Validate presence of "valid bin" (strict)
    const int c_valid = col_index_or_die(csv, "valid bin");

    // Iterate rows; only process valid bins (valid bin == 1)
    size_t processed = 0;
    for (size_t r = 0; r < csv.n_rows(); ++r) {
        const std::string& vb = csv.cell(r, (size_t)c_valid);
        int valid = 0;
        if (!vb.empty()) {
            char* endp = nullptr;
            long v = std::strtol(vb.c_str(), &endp, 10);
            if (endp != vb.c_str()) valid = (int)v;
        }
        if (valid != 1) continue;

        compute_signal_and_acc_corr_for_row(csv, r);
        ++processed;
    }

    if (!csv.save(out_csv)) {
        return false;
    }

    std::cout << "[total_counts] Processed valid bins: " << processed << "\n";
    std::cout << "[total_counts] Wrote: " << out_csv << "\n";
    return true;
}