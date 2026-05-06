// acceptance.cpp
// -----------------------------------------------------------------------------
// DVCS acceptance from already-counted non-radiative MC yields in the CSV.
//
// New refactored behavior:
//   - No ROOT tree loops.
//   - No event-level cuts here.
//   - Uses generated and reconstructed MC yields already saved by total_counts.cpp.
//   - Computes:
//
//       acceptance, <period> = N_rec_current_corrected / N_gen
//
//     where:
//
//       N_gen = generated yield, ep->epg, mc, <period>
//
//       N_rec_current_corrected = sum over topologies:
//         reconstructed current corrected yield, ep->epg, (FD, FD), mc, <period>
//         reconstructed current corrected yield, ep->epg, (CD, FD), mc, <period>
//         reconstructed current corrected yield, ep->epg, (CD, FT), mc, <period>
//
//   - Writes acceptance as:
//
//       (value,stat,sys)
//
//     with sys = 0 for now.
//
//   - Produces per-period acceptance-vs-phi plots under:
//       <out_root_dir>/acceptance/<PeriodDir>/
//
// Notes:
//   - Global cuts, topology cuts, and 3-sigma cuts are intentionally not applied
//     here because the reconstructed MC yield columns were filled upstream after
//     those cuts.
// -----------------------------------------------------------------------------

#include "acceptance.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TAxis.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

static inline void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static inline bool is_finite_positive(double x) {
    return std::isfinite(x) && x > 0.0;
}

static inline bool is_finite_nonnegative(double x) {
    return std::isfinite(x) && x >= 0.0;
}

static std::string sanitize_dir_token(const std::string& s) {
    std::string out;
    out.reserve(s.size());

    for (char c : s) {
        if (c == ' ') {
            out.push_back('_');
        } else {
            out.push_back(c);
        }
    }

    return out;
}

static std::string canonical_period_dir(const std::string& label) {
    if (label == "Fa18 Inb") return "Fa18_Inb";
    if (label == "Fa18 Out") return "Fa18_Out";
    if (label == "Sp19 Inb") return "Sp19_Inb";
    if (label == "Sp18 Inb") return "Sp18_Inb";
    if (label == "Sp18 Out") return "Sp18_Out";

    std::ostringstream ss;
    ss << "[acceptance] FATAL: unknown period label: " << label;
    fatal(ss.str());

    return "";
}

static void mkdir_p(const std::string& path) {
    if (!path.empty()) {
        gSystem->mkdir(path.c_str(), true);
    }
}

static const std::vector<std::string>& periods() {
    static const std::vector<std::string> v = {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb",
        "Sp18 Out"
    };

    return v;
}

static const std::vector<std::string>& topologies() {
    static const std::vector<std::string> v = {
        "(FD, FD)",
        "(CD, FD)",
        "(CD, FT)"
    };

    return v;
}

// -----------------------------------------------------------------------------
// CSV I/O
// -----------------------------------------------------------------------------

struct CSV {
    std::vector<std::string> header;
    std::unordered_map<std::string, int> index;
    std::vector<std::vector<std::string>> rows;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());

    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
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

        const bool quote =
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

        if (quote) {
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

        if (i + 1 < fields.size()) {
            oss << ',';
        }
    }

    return oss.str();
}

static void load_csv_strict(const std::string& path, CSV& csv) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[acceptance] FATAL: cannot open CSV: " << path;
        fatal(ss.str());
    }

    std::string line;

    if (!std::getline(fin, line)) {
        std::ostringstream ss;
        ss << "[acceptance] FATAL: empty CSV: " << path;
        fatal(ss.str());
    }

    csv.header = split_csv_line(line);
    csv.index.clear();

    for (int i = 0; i < (int)csv.header.size(); ++i) {
        const std::string& name = csv.header[(size_t)i];

        if (csv.index.find(name) != csv.index.end()) {
            std::ostringstream ss;
            ss << "[acceptance] FATAL: duplicate CSV column: " << name;
            fatal(ss.str());
        }

        csv.index[name] = i;
    }

    csv.rows.clear();

    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> row = split_csv_line(line);

        if (row.size() < csv.header.size()) {
            row.resize(csv.header.size(), "");
        }

        if (row.size() != csv.header.size()) {
            std::ostringstream ss;
            ss << "[acceptance] FATAL: CSV row width mismatch in " << path
               << ". Row has " << row.size()
               << " cells, header has " << csv.header.size() << ".";
            fatal(ss.str());
        }

        csv.rows.push_back(std::move(row));
    }
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";

    {
        std::ofstream fout(tmp);

        if (!fout.is_open()) {
            std::ostringstream ss;
            ss << "[acceptance] FATAL: cannot write temporary CSV: " << tmp;
            fatal(ss.str());
        }

        fout << join_csv_row(csv.header) << "\n";

        for (const auto& row : csv.rows) {
            if (row.size() != csv.header.size()) {
                std::ostringstream ss;
                ss << "[acceptance] FATAL: row width mismatch during write.";
                fatal(ss.str());
            }

            fout << join_csv_row(row) << "\n";
        }

        fout.close();

        if (!fout) {
            std::ostringstream ss;
            ss << "[acceptance] FATAL: failed writing temporary CSV: " << tmp;
            fatal(ss.str());
        }
    }

    std::error_code ec;
    std::filesystem::rename(tmp, path, ec);

    if (ec) {
        std::remove(path.c_str());
        std::filesystem::rename(tmp, path, ec);

        if (ec) {
            std::ostringstream ss;
            ss << "[acceptance] FATAL: atomic rename failed from "
               << tmp << " to " << path << ": " << ec.message();
            fatal(ss.str());
        }
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);

    if (it == csv.index.end()) {
        std::ostringstream ss;
        ss << "[acceptance] FATAL: missing required CSV column: " << name;
        fatal(ss.str());
    }

    return it->second;
}

static bool parse_first_number(const std::string& s, double& value) {
    value = std::numeric_limits<double>::quiet_NaN();

    if (s.empty()) {
        return false;
    }

    size_t i = 0;

    while (i < s.size()) {
        const char c = s[i];

        const bool maybe_start =
            (c >= '0' && c <= '9') ||
            c == '-' ||
            c == '+' ||
            c == '.';

        if (maybe_start) {
            break;
        }

        ++i;
    }

    if (i >= s.size()) {
        return false;
    }

    char* endp = nullptr;
    value = std::strtod(s.c_str() + i, &endp);

    if (endp == s.c_str() + i) {
        value = std::numeric_limits<double>::quiet_NaN();
        return false;
    }

    return std::isfinite(value);
}

static double cell_value_or_zero(const CSV& csv, int row, int col) {
    if (row < 0 || row >= (int)csv.rows.size()) {
        return 0.0;
    }

    if (col < 0 || col >= (int)csv.header.size()) {
        return 0.0;
    }

    double v = 0.0;

    if (!parse_first_number(csv.rows[(size_t)row][(size_t)col], v)) {
        return 0.0;
    }

    if (!std::isfinite(v)) {
        return 0.0;
    }

    return v;
}

static double cell_value_strict(const CSV& csv,
                                int row,
                                int col,
                                const std::string& label) {
    if (row < 0 || row >= (int)csv.rows.size()) {
        std::ostringstream ss;
        ss << "[acceptance] FATAL: row index out of range for " << label;
        fatal(ss.str());
    }

    if (col < 0 || col >= (int)csv.header.size()) {
        std::ostringstream ss;
        ss << "[acceptance] FATAL: column index out of range for " << label;
        fatal(ss.str());
    }

    double v = 0.0;

    if (!parse_first_number(csv.rows[(size_t)row][(size_t)col], v)) {
        std::ostringstream ss;
        ss << "[acceptance] FATAL: cannot parse numeric cell for " << label
           << " from value '" << csv.rows[(size_t)row][(size_t)col] << "'";
        fatal(ss.str());
    }

    return v;
}

static std::string triple_string(double value, double stat, double sys) {
    std::ostringstream ss;
    ss << std::setprecision(12)
       << "(" << value << "," << stat << "," << sys << ")";
    return ss.str();
}

// -----------------------------------------------------------------------------
// Column naming
// -----------------------------------------------------------------------------

static std::string col_generated(const std::string& period) {
    std::ostringstream ss;
    ss << "generated yield, ep->epg, mc, " << period;
    return ss.str();
}

static std::string col_reconstructed_current_corrected(const std::string& topo,
                                                       const std::string& period) {
    std::ostringstream ss;
    ss << "reconstructed current corrected yield, ep->epg, "
       << topo << ", mc, " << period;
    return ss.str();
}

static std::string col_acceptance(const std::string& period) {
    std::ostringstream ss;
    ss << "acceptance, " << period;
    return ss.str();
}

static std::string col_phiavg(const std::string& period) {
    std::ostringstream ss;
    ss << "phiavg, " << period;
    return ss.str();
}

static std::string col_xBavg(const std::string& period) {
    std::ostringstream ss;
    ss << "xBavg, " << period;
    return ss.str();
}

// -----------------------------------------------------------------------------
// Row bins and computed acceptance
// -----------------------------------------------------------------------------

struct RowBin {
    double xBmin = std::numeric_limits<double>::quiet_NaN();
    double xBmax = std::numeric_limits<double>::quiet_NaN();

    double Q2min = std::numeric_limits<double>::quiet_NaN();
    double Q2max = std::numeric_limits<double>::quiet_NaN();

    double tmin = std::numeric_limits<double>::quiet_NaN();
    double tmax = std::numeric_limits<double>::quiet_NaN();

    double phimin = std::numeric_limits<double>::quiet_NaN();
    double phimax = std::numeric_limits<double>::quiet_NaN();

    bool valid = false;
};

struct AcceptancePoint {
    double value = std::numeric_limits<double>::quiet_NaN();
    double stat = std::numeric_limits<double>::quiet_NaN();
    double sys = 0.0;

    double n_gen = 0.0;
    double n_rec = 0.0;
};

static bool valid_cell_to_bool(const std::string& s) {
    return s == "1" || s == "1.0" || s == "true" || s == "TRUE";
}

static double to_double_strict_cell(const CSV& csv,
                                    int row,
                                    const std::string& colname) {
    const int c = col_strict(csv, colname);
    return cell_value_strict(csv, row, c, colname);
}

static std::vector<RowBin> load_row_bins(const CSV& csv) {
    const int c_xBmin = col_strict(csv, "xBmin");
    const int c_xBmax = col_strict(csv, "xBmax");
    const int c_Q2min = col_strict(csv, "Q2min");
    const int c_Q2max = col_strict(csv, "Q2max");
    const int c_tmin = col_strict(csv, "t_abs_min");
    const int c_tmax = col_strict(csv, "t_abs_max");
    const int c_phimin = col_strict(csv, "phimin");
    const int c_phimax = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    std::vector<RowBin> out;
    out.reserve(csv.rows.size());

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        RowBin b;

        b.xBmin = cell_value_strict(csv, r, c_xBmin, "xBmin");
        b.xBmax = cell_value_strict(csv, r, c_xBmax, "xBmax");
        b.Q2min = cell_value_strict(csv, r, c_Q2min, "Q2min");
        b.Q2max = cell_value_strict(csv, r, c_Q2max, "Q2max");
        b.tmin = cell_value_strict(csv, r, c_tmin, "t_abs_min");
        b.tmax = cell_value_strict(csv, r, c_tmax, "t_abs_max");
        b.phimin = cell_value_strict(csv, r, c_phimin, "phimin");
        b.phimax = cell_value_strict(csv, r, c_phimax, "phimax");
        b.valid = valid_cell_to_bool(csv.rows[(size_t)r][(size_t)c_valid]);

        out.push_back(b);
    }

    return out;
}

static double bin_center(double a, double b) {
    return 0.5 * (a + b);
}

static double wrap_phi_deg(double phi) {
    double x = std::fmod(phi, 360.0);

    if (x < 0.0) {
        x += 360.0;
    }

    if (x >= 360.0) {
        x = std::nextafter(360.0, 0.0);
    }

    return x;
}

static double binomial_stat(double n_rec, double n_gen) {
    if (!is_finite_positive(n_gen)) {
        return 0.0;
    }

    const double a = n_rec / n_gen;

    if (!std::isfinite(a)) {
        return 0.0;
    }

    // With current-corrected reconstructed yields, n_rec can be slightly above
    // n_gen. Clamp only the variance factor, not the acceptance value itself.
    const double a_for_var = std::max(0.0, std::min(1.0, a));
    const double var = a_for_var * (1.0 - a_for_var) / n_gen;

    if (!(var > 0.0)) {
        return 0.0;
    }

    return std::sqrt(var);
}

static std::vector<AcceptancePoint> compute_acceptance_for_period(const CSV& csv,
                                                                  const std::string& period) {
    const int c_gen = col_strict(csv, col_generated(period));
    const int c_acc = col_strict(csv, col_acceptance(period));
    (void)c_acc;

    std::vector<int> rec_cols;
    rec_cols.reserve(topologies().size());

    for (const auto& topo : topologies()) {
        rec_cols.push_back(col_strict(csv, col_reconstructed_current_corrected(topo, period)));
    }

    std::vector<AcceptancePoint> out;
    out.resize(csv.rows.size());

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        const double n_gen = cell_value_or_zero(csv, r, c_gen);

        double n_rec = 0.0;

        for (int c : rec_cols) {
            n_rec += cell_value_or_zero(csv, r, c);
        }

        AcceptancePoint p;
        p.n_gen = n_gen;
        p.n_rec = n_rec;

        if (is_finite_positive(n_gen) && is_finite_nonnegative(n_rec)) {
            p.value = n_rec / n_gen;
            p.stat = binomial_stat(n_rec, n_gen);
            p.sys = 0.0;
        } else {
            p.value = 0.0;
            p.stat = 0.0;
            p.sys = 0.0;
        }

        out[(size_t)r] = p;
    }

    return out;
}

static void write_acceptance_to_csv(CSV& csv,
                                    const std::string& period,
                                    const std::vector<AcceptancePoint>& acc) {
    const int c_acc = col_strict(csv, col_acceptance(period));

    if (acc.size() != csv.rows.size()) {
        fatal("[acceptance] FATAL: acceptance vector size does not match CSV rows.");
    }

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        const AcceptancePoint& p = acc[(size_t)r];
        csv.rows[(size_t)r][(size_t)c_acc] = triple_string(p.value, p.stat, p.sys);
    }
}

// -----------------------------------------------------------------------------
// Plotting
// -----------------------------------------------------------------------------

struct Edge {
    double a = 0.0;
    double b = 0.0;
};

static bool same_edge(const Edge& x, const Edge& y) {
    return x.a == y.a && x.b == y.b;
}

static void add_unique_edge(std::vector<Edge>& edges, const Edge& e) {
    auto it = std::find_if(edges.begin(), edges.end(), [&](const Edge& z) {
        return same_edge(e, z);
    });

    if (it == edges.end()) {
        edges.push_back(e);
    }
}

static int find_edge_index(const std::vector<Edge>& edges, const Edge& e) {
    for (int i = 0; i < (int)edges.size(); ++i) {
        if (same_edge(edges[(size_t)i], e)) {
            return i;
        }
    }

    return -1;
}

static std::vector<Edge> unique_xB_edges(const std::vector<RowBin>& rows) {
    std::vector<Edge> out;

    for (const RowBin& r : rows) {
        if (!r.valid) {
            continue;
        }

        add_unique_edge(out, Edge{r.xBmin, r.xBmax});
    }

    std::sort(out.begin(), out.end(), [](const Edge& a, const Edge& b) {
        if (a.a != b.a) return a.a < b.a;
        return a.b < b.b;
    });

    return out;
}

static double phi_x_value(const CSV& csv,
                          const RowBin& rowbin,
                          int row_index,
                          int c_phiavg) {
    if (c_phiavg >= 0) {
        double v = 0.0;

        if (parse_first_number(csv.rows[(size_t)row_index][(size_t)c_phiavg], v)) {
            return wrap_phi_deg(v);
        }
    }

    return wrap_phi_deg(bin_center(rowbin.phimin, rowbin.phimax));
}

static std::string xB_label_for_canvas(const CSV& csv,
                                       const std::vector<RowBin>& rows,
                                       const std::vector<int>& row_indices,
                                       const std::string& period,
                                       const Edge& xb) {
    const int c_xBavg = csv.index.count(col_xBavg(period)) ? csv.index.at(col_xBavg(period)) : -1;

    if (c_xBavg >= 0) {
        for (int r : row_indices) {
            double v = 0.0;

            if (parse_first_number(csv.rows[(size_t)r][(size_t)c_xBavg], v)) {
                std::ostringstream ss;
                ss << "x_{B}=" << std::fixed << std::setprecision(3) << v;
                return ss.str();
            }
        }
    }

    (void)rows;

    std::ostringstream ss;
    ss << "x_{B} in [" << std::fixed << std::setprecision(2) << xb.a
       << "," << std::fixed << std::setprecision(2) << xb.b << ")";
    return ss.str();
}

static void plot_acceptance_for_period(const CSV& csv,
                                       const std::vector<RowBin>& rows,
                                       const std::string& period,
                                       const std::vector<AcceptancePoint>& acc,
                                       const std::string& out_root_dir) {
    const std::string outdir =
        out_root_dir + "/acceptance/" + canonical_period_dir(period);

    mkdir_p(outdir);

    const std::vector<Edge> xbs = unique_xB_edges(rows);
    const int c_phiavg = csv.index.count(col_phiavg(period)) ? csv.index.at(col_phiavg(period)) : -1;

    for (const Edge& xb : xbs) {
        std::vector<int> row_indices;
        std::vector<Edge> q_edges;
        std::vector<Edge> t_edges;

        for (int r = 0; r < (int)rows.size(); ++r) {
            const RowBin& rb = rows[(size_t)r];

            if (!rb.valid) {
                continue;
            }

            if (!(rb.xBmin == xb.a && rb.xBmax == xb.b)) {
                continue;
            }

            row_indices.push_back(r);
            add_unique_edge(q_edges, Edge{rb.Q2min, rb.Q2max});
            add_unique_edge(t_edges, Edge{rb.tmin, rb.tmax});
        }

        if (row_indices.empty()) {
            continue;
        }

        std::sort(q_edges.begin(), q_edges.end(), [](const Edge& a, const Edge& b) {
            if (a.a != b.a) return a.a < b.a;
            return a.b < b.b;
        });

        std::sort(t_edges.begin(), t_edges.end(), [](const Edge& a, const Edge& b) {
            if (a.a != b.a) return a.a < b.a;
            return a.b < b.b;
        });

        const int ncols = (int)q_edges.size();
        const int nrows = (int)t_edges.size();

        if (ncols <= 0 || nrows <= 0) {
            continue;
        }

        const int W = 300 * ncols + 160;
        const int H = 260 * nrows + 240;

        TCanvas c("c_acceptance", "", W, H);

        TPad* top = new TPad("top", "top", 0.0, 0.90, 1.0, 1.0);
        TPad* grid = new TPad("grid", "grid", 0.0, 0.0, 1.0, 0.90);

        top->SetFillStyle(0);
        grid->SetFillStyle(0);

        top->Draw();
        grid->Draw();

        top->cd();

        TLatex title;
        title.SetNDC(true);
        title.SetTextFont(42);
        title.SetTextSize(0.45);

        const std::string xb_text =
            xB_label_for_canvas(csv, rows, row_indices, period, xb);

        {
            std::ostringstream ss;
            ss << "Acceptance   " << period << "   " << xb_text;
            title.DrawLatex(0.05, 0.35, ss.str().c_str());
        }

        grid->cd();
        grid->Divide(ncols, nrows, 0.0, 0.0);

        for (int it = 0; it < nrows; ++it) {
            for (int iq = 0; iq < ncols; ++iq) {
                const int pad_index = it * ncols + iq + 1;
                grid->cd(pad_index);

                gPad->SetLeftMargin(0.160);
                gPad->SetRightMargin(0.07);
                gPad->SetTopMargin(0.22);
                gPad->SetBottomMargin(0.18);
                gPad->SetGrid(1, 1);
                gPad->SetTickx(1);
                gPad->SetTicky(1);

                std::vector<std::pair<double, AcceptancePoint>> points;

                for (int r : row_indices) {
                    const RowBin& rb = rows[(size_t)r];

                    const int q_idx = find_edge_index(q_edges, Edge{rb.Q2min, rb.Q2max});
                    const int t_idx = find_edge_index(t_edges, Edge{rb.tmin, rb.tmax});

                    if (q_idx != iq || t_idx != it) {
                        continue;
                    }

                    const double x = phi_x_value(csv, rb, r, c_phiavg);
                    points.push_back(std::make_pair(x, acc[(size_t)r]));
                }

                std::sort(points.begin(), points.end(),
                          [](const std::pair<double, AcceptancePoint>& a,
                             const std::pair<double, AcceptancePoint>& b) {
                              return a.first < b.first;
                          });

                TGraphErrors* g = new TGraphErrors();

                double ymax = 0.0;

                for (int i = 0; i < (int)points.size(); ++i) {
                    const double x = points[(size_t)i].first;
                    const AcceptancePoint& p = points[(size_t)i].second;

                    g->SetPoint(i, x, p.value);
                    g->SetPointError(i, 0.0, p.stat);

                    ymax = std::max(ymax, p.value + p.stat);
                }

                if (!(ymax > 0.0)) {
                    ymax = 1.0;
                }

                TH1F* frame = (TH1F*)gPad->DrawFrame(0.0, 0.0, 360.0, 1.20 * ymax);
                frame->SetTitle("");
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Acceptance");
                frame->GetXaxis()->CenterTitle(true);
                frame->GetYaxis()->CenterTitle(true);
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetLabelSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.060);
                frame->GetXaxis()->SetTitleSize(0.070);
                frame->GetYaxis()->SetTitleSize(0.070);
                frame->GetXaxis()->SetTitleOffset(1.05);
                frame->GetYaxis()->SetTitleOffset(1.10);

                g->SetMarkerStyle(20);
                g->SetMarkerSize(0.8);
                g->SetMarkerColor(kBlack);
                g->SetLineColor(kBlack);
                g->SetLineWidth(1);
                g->Draw("PE1 SAME");

                TLatex txt;
                txt.SetNDC(true);
                txt.SetTextFont(42);
                txt.SetTextSize(0.060);

                {
                    std::ostringstream ss;
                    ss << "Q^{2}=" << std::fixed << std::setprecision(2)
                       << bin_center(q_edges[(size_t)iq].a, q_edges[(size_t)iq].b)
                       << "  |t|=" << std::fixed << std::setprecision(2)
                       << bin_center(t_edges[(size_t)it].a, t_edges[(size_t)it].b);
                    txt.DrawLatex(0.12, 0.83, ss.str().c_str());
                }
            }
        }

        const int idx = (int)std::llround(xb.a * 1000.0);

        std::ostringstream fname;
        fname << outdir << "/acceptance_"
              << canonical_period_dir(period)
              << "_xB_" << idx << ".png";

        c.SaveAs(fname.str().c_str());

        delete top;
        delete grid;
    }
}

} // namespace

bool update_acceptance_csv(const std::string& csv_path,
                           const std::map<std::string, TTree*>& genMcTrees,
                           const std::map<std::string, TTree*>& recMcTrees,
                           const std::string& combined_cuts_json,
                           const std::string& global_cuts_json,
                           const std::string& out_root_dir) {
    (void)genMcTrees;
    (void)recMcTrees;
    (void)combined_cuts_json;
    (void)global_cuts_json;

    try {
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gStyle->SetOptStat(0);

        CSV csv;
        load_csv_strict(csv_path, csv);

        const std::vector<RowBin> rows = load_row_bins(csv);

        std::map<std::string, std::vector<AcceptancePoint>> acc_by_period;

        for (const std::string& period : periods()) {
            std::cout << "[acceptance] Computing CSV-only acceptance for "
                      << period << std::endl;

            std::vector<AcceptancePoint> acc =
                compute_acceptance_for_period(csv, period);

            write_acceptance_to_csv(csv, period, acc);
            acc_by_period[period] = std::move(acc);
        }

        write_csv_atomic(csv_path, csv);

        std::cout << "[acceptance] Updated acceptance columns in: "
                  << csv_path << std::endl;

        for (const std::string& period : periods()) {
            plot_acceptance_for_period(csv,
                                       rows,
                                       period,
                                       acc_by_period.at(period),
                                       out_root_dir);
        }

        std::cout << "[acceptance] Acceptance plots written under: "
                  << out_root_dir << "/acceptance" << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}