#include "cut_variation_systematics.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLine.h>
#include <TStyle.h>
#include <TAxis.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Triple {
    double value = std::numeric_limits<double>::quiet_NaN();
    double stat = std::numeric_limits<double>::quiet_NaN();
    double sys = std::numeric_limits<double>::quiet_NaN();
    bool ok = false;
};

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string, size_t> index;
};

struct DiagnosticRow {
    int bin_index = -1;
    std::string bin_name;
    double xb_min = 0.0;
    double xb_max = 0.0;
    double q2_min = 0.0;
    double q2_max = 0.0;
    double t_min = 0.0;
    double t_max = 0.0;
    double phi_min = 0.0;
    double phi_max = 0.0;
    double phi = 0.0;

    Triple nominal;
    Triple excl_loose;
    Triple excl_tight;
    Triple fid_loose;
    Triple fid_tight;

    double excl_delta_loose_raw = 0.0;
    double excl_delta_tight_raw = 0.0;
    double excl_barlow_loose = 0.0;
    double excl_barlow_tight = 0.0;
    bool excl_keep_loose = false;
    bool excl_keep_tight = false;
    double excl_tight_relative_difference = 0.0;
    bool excl_use_loose_only = false;
    double excl_raw_abs = 0.0;
    double excl_final_abs = 0.0;

    double fid_delta_loose_raw = 0.0;
    double fid_delta_tight_raw = 0.0;
    double fid_barlow_loose = 0.0;
    double fid_barlow_tight = 0.0;
    bool fid_keep_loose = false;
    bool fid_keep_tight = false;
    double fid_tight_relative_difference = 0.0;
    bool fid_use_loose_only = false;
    double fid_raw_abs = 0.0;
    double fid_final_abs = 0.0;
};

std::string trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string normalize_csv_field(std::string field) {
    // Remove complete wrapper-quote layers left by older CSV rewrites while
    // preserving quotation marks that are genuinely internal to a field.
    bool changed = true;
    while (changed && field.size() >= 2) {
        changed = false;
        if (field.front() == '"' && field.back() == '"') {
            field = field.substr(1, field.size() - 2);
            std::string unescaped;
            unescaped.reserve(field.size());
            for (size_t i = 0; i < field.size(); ++i) {
                if (field[i] == '"' && i + 1 < field.size() && field[i + 1] == '"') {
                    unescaped.push_back('"');
                    ++i;
                } else {
                    unescaped.push_back(field[i]);
                }
            }
            field.swap(unescaped);
            changed = true;
        }
    }
    return field;
}

std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            out.push_back(normalize_csv_field(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(normalize_csv_field(cur));
    return out;
}

std::string csv_quote(const std::string& s) {
    if (s.find_first_of(",\"\r\n") == std::string::npos) return s;
    std::string q = "\"";
    for (char c : s) {
        if (c == '"') q += "\"\"";
        else q.push_back(c);
    }
    q += '"';
    return q;
}

CsvTable read_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Could not open CSV: " + path);
    CsvTable t;
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("Empty CSV: " + path);
    t.header = split_csv_line(line);
    for (size_t i = 0; i < t.header.size(); ++i) t.index[trim(t.header[i])] = i;
    while (std::getline(in, line)) {
        auto row = split_csv_line(line);
        row.resize(t.header.size());
        t.rows.push_back(std::move(row));
    }
    return t;
}

void write_csv(const CsvTable& t, const std::string& path) {
    fs::create_directories(fs::path(path).parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Could not write CSV: " + path);
    for (size_t i = 0; i < t.header.size(); ++i) {
        if (i) out << ',';
        out << csv_quote(t.header[i]);
    }
    out << '\n';
    for (const auto& row : t.rows) {
        for (size_t i = 0; i < t.header.size(); ++i) {
            if (i) out << ',';
            out << csv_quote(i < row.size() ? row[i] : "");
        }
        out << '\n';
    }
}

size_t require_col(const CsvTable& t, const std::string& name, const std::string& path) {
    const auto it = t.index.find(name);
    if (it == t.index.end()) throw std::runtime_error("Missing column '" + name + "' in " + path);
    return it->second;
}

size_t ensure_col(CsvTable& t, const std::string& name) {
    const auto it = t.index.find(name);
    if (it != t.index.end()) return it->second;
    const size_t idx = t.header.size();
    t.header.push_back(name);
    t.index[name] = idx;
    for (auto& row : t.rows) row.push_back("");
    return idx;
}

double to_double(const std::string& s, double fallback = 0.0) {
    try {
        size_t pos = 0;
        const double v = std::stod(trim(s), &pos);
        return std::isfinite(v) ? v : fallback;
    } catch (...) {
        return fallback;
    }
}

Triple parse_triple(std::string s) {
    Triple out;
    s = normalize_csv_field(trim(s));
    if (s.empty()) return out;
    while (s.size() >= 2 && ((s.front() == '(' && s.back() == ')') ||
                             (s.front() == '[' && s.back() == ']') ||
                             (s.front() == '{' && s.back() == '}'))) {
        s = trim(s.substr(1, s.size() - 2));
    }
    for (char& c : s) if (c == ';') c = ',';
    const auto fields = split_csv_line(s);
    if (fields.empty()) return out;
    out.value = to_double(fields[0], std::numeric_limits<double>::quiet_NaN());
    out.stat = fields.size() > 1 ? to_double(fields[1], 0.0) : 0.0;
    out.sys = fields.size() > 2 ? to_double(fields[2], 0.0) : 0.0;
    out.ok = std::isfinite(out.value) && std::isfinite(out.stat);
    return out;
}

std::string format_number(double v) {
    if (!std::isfinite(v)) return "";
    std::ostringstream ss;
    ss << std::setprecision(12) << v;
    return ss.str();
}

std::string row_key(const CsvTable& t, const std::vector<std::string>& row) {
    static const std::vector<std::string> candidates = {"bin index", "", "Bin Name"};
    for (const auto& name : candidates) {
        const auto it = t.index.find(name);
        if (it != t.index.end() && it->second < row.size() && !trim(row[it->second]).empty()) {
            return trim(row[it->second]);
        }
    }
    throw std::runtime_error("Could not construct a row key; expected bin index or Bin Name.");
}

std::unordered_map<std::string, size_t> build_row_map(const CsvTable& t) {
    std::unordered_map<std::string, size_t> out;
    for (size_t i = 0; i < t.rows.size(); ++i) out[row_key(t, t.rows[i])] = i;
    return out;
}

struct VariationValidationSummary {
    std::size_t nominal_valid = 0;
    std::size_t variation_valid = 0;
    std::vector<std::string> invalid_keys;
};

VariationValidationSummary validate_variation_cross_sections(
    const CsvTable& nominal,
    std::size_t nominal_col,
    const CsvTable& variation,
    std::size_t variation_col,
    const std::unordered_map<std::string, size_t>& variation_rows) {

    VariationValidationSummary summary;
    for (const auto& row : nominal.rows) {
        const Triple nom = parse_triple(row[nominal_col]);
        if (!nom.ok) continue;
        ++summary.nominal_valid;

        const std::string key = row_key(nominal, row);
        const auto it = variation_rows.find(key);
        if (it == variation_rows.end()) {
            if (summary.invalid_keys.size() < 8) summary.invalid_keys.push_back(key + " (missing row)");
            continue;
        }

        const Triple varied = parse_triple(variation.rows[it->second][variation_col]);
        if (varied.ok) {
            ++summary.variation_valid;
        } else if (summary.invalid_keys.size() < 8) {
            summary.invalid_keys.push_back(key + " (unreadable cross-section tuple)");
        }
    }
    return summary;
}

void require_valid_variation(const std::string& label,
                             const VariationValidationSummary& summary) {
    std::cout << "[cut-systematics] " << std::left << std::setw(26) << label
              << " valid bins: " << summary.variation_valid
              << " / " << summary.nominal_valid << '\n';

    if (summary.nominal_valid == 0) {
        throw std::runtime_error("The nominal cross-section column contains no readable tuples.");
    }
    if (summary.variation_valid != summary.nominal_valid) {
        std::ostringstream message;
        message << label << " contains " << summary.variation_valid << " readable cross sections for "
                << summary.nominal_valid << " nominally populated bins. The nominal CSV was not modified.";
        if (!summary.invalid_keys.empty()) {
            message << " Example failures: ";
            for (std::size_t i = 0; i < summary.invalid_keys.size(); ++i) {
                if (i) message << "; ";
                message << summary.invalid_keys[i];
            }
        }
        throw std::runtime_error(message.str());
    }
}

double barlow_value(const Triple& varied, const Triple& nominal) {
    if (!varied.ok || !nominal.ok) return 0.0;
    const double denom2 = std::fabs(varied.stat * varied.stat - nominal.stat * nominal.stat);
    const double diff = std::fabs(varied.value - nominal.value);
    if (denom2 <= 1e-30) return diff <= 1e-30 ? 0.0 : std::numeric_limits<double>::infinity();
    return diff / std::sqrt(denom2);
}

double relative_difference(double delta, double nominal) {
    const double denom = std::fabs(nominal);
    if (denom <= 1e-30) {
        return std::fabs(delta) <= 1e-30
            ? 0.0
            : std::numeric_limits<double>::infinity();
    }
    return std::fabs(delta) / denom;
}

double pass1_cut_systematic_absolute(double loose_delta,
                                     double tight_delta,
                                     bool use_loose_only) {
    if (use_loose_only) return std::fabs(loose_delta);
    return 0.5 * (std::fabs(loose_delta) + std::fabs(tight_delta));
}

std::string xbin_label(double lo, double hi) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2) << lo << " < x_{B} < " << hi;
    return ss.str();
}


double quantile(std::vector<double> values, double q) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(values.begin(), values.end());
    if (q <= 0.0) return values.front();
    if (q >= 1.0) return values.back();

    const double p = q * static_cast<double>(values.size() - 1);
    const size_t i = static_cast<size_t>(std::floor(p));
    const size_t j = static_cast<size_t>(std::ceil(p));
    const double f = p - static_cast<double>(i);
    return values[i] * (1.0 - f) + values[j] * f;
}

double mean_value(const std::vector<double>& values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    double sum = 0.0;
    for (double v : values) sum += v;
    return sum / static_cast<double>(values.size());
}

void style_note_canvas(TCanvas& c, double top_margin = 0.14) {
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.035);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(top_margin);
    c.SetTicks(1, 1);
}

void style_note_axes(TAxis* x, TAxis* y) {
    if (x) {
        x->SetTitleFont(42);
        x->SetLabelFont(42);
        x->SetTitleSize(0.043);
        x->SetLabelSize(0.037);
        x->SetTitleOffset(1.12);
    }
    if (y) {
        y->SetTitleFont(42);
        y->SetLabelFont(42);
        y->SetTitleSize(0.043);
        y->SetLabelSize(0.037);
        y->SetTitleOffset(1.24);
    }
}

void draw_note_title(const char* title, const char* subtitle = nullptr) {
    TLatex text;
    text.SetNDC();
    text.SetTextFont(42);
    text.SetTextAlign(13);

    text.SetTextSize(0.031);
    text.DrawLatex(0.13, 0.965, title);

    if (subtitle && *subtitle) {
        text.SetTextSize(0.024);
        text.DrawLatex(0.13, 0.925, subtitle);
    }
}

std::string kinematic_subtitle(double xb0, double xb1,
                               double q0, double q1,
                               double t0, double t1) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(3)
       << xb0 << " < x_{B} < " << xb1
       << ",  " << q0 << " < Q^{2} < " << q1 << " GeV^{2}"
       << ",  " << t0 << " < |t| < " << t1 << " GeV^{2}";
    return ss.str();
}

void make_analysis_note_plots(const std::vector<DiagnosticRow>& rows,
                              const std::string& outdir) {
    fs::create_directories(outdir);
    gStyle->SetOptStat(0);

    std::vector<double> raw_relative;
    std::vector<double> final_relative;

    for (const auto& r : rows) {
        if (!r.nominal.ok || std::fabs(r.nominal.value) <= 1e-30) continue;
        const double denom = std::fabs(r.nominal.value);
        raw_relative.push_back(100.0 * r.excl_raw_abs / denom);
        final_relative.push_back(100.0 * r.excl_final_abs / denom);
    }

    if (final_relative.empty()) {
        std::cerr << "[cut-systematics] WARNING: no populated bins for "
                  << "analysis-note exclusivity plots.\n";
        return;
    }

    // ------------------------------------------------------------------
    // Compact numerical summary used later when writing the note.
    // ------------------------------------------------------------------
    {
        const std::string summary_path =
            (fs::path(outdir) / "exclusivity_systematic_summary.csv").string();
        std::ofstream out(summary_path);
        if (!out) {
            throw std::runtime_error(
                "Could not write analysis-note summary: " + summary_path);
        }

        const auto loose_retained = std::count_if(
            rows.begin(), rows.end(),
            [](const DiagnosticRow& r) {
                return r.nominal.ok && r.excl_keep_loose;
            });
        const auto tight_retained = std::count_if(
            rows.begin(), rows.end(),
            [](const DiagnosticRow& r) {
                return r.nominal.ok && r.excl_keep_tight;
            });
        const auto loose_only = std::count_if(
            rows.begin(), rows.end(),
            [](const DiagnosticRow& r) {
                return r.nominal.ok && r.excl_use_loose_only;
            });
        const auto populated = std::count_if(
            rows.begin(), rows.end(),
            [](const DiagnosticRow& r) { return r.nominal.ok; });

        out << "quantity,value\n";
        out << "populated bins," << populated << "\n";
        out << "mean final relative systematic (%),"
            << mean_value(final_relative) << "\n";
        out << "median final relative systematic (%),"
            << quantile(final_relative, 0.50) << "\n";
        out << "16th percentile final relative systematic (%),"
            << quantile(final_relative, 0.16) << "\n";
        out << "84th percentile final relative systematic (%),"
            << quantile(final_relative, 0.84) << "\n";
        out << "95th percentile final relative systematic (%),"
            << quantile(final_relative, 0.95) << "\n";
        out << "mean raw relative systematic (%),"
            << mean_value(raw_relative) << "\n";
        out << "loose variation retained by Barlow," << loose_retained << "\n";
        out << "tight variation retained by Barlow," << tight_retained << "\n";
        out << "tight-instability loose-only bins," << loose_only << "\n";
        out << "Barlow threshold,1.0\n";
        out << "tight-instability threshold,0.50\n";
        out << "nominal containment,0.95\n";
        out << "loose containment,0.98\n";
        out << "tight containment,0.90\n";
    }

    // ------------------------------------------------------------------
    // 1) Distribution of relative uncertainty.
    //    The x range is based on the populated-bin distribution so a few
    //    exceptional tails do not compress the informative region.
    // ------------------------------------------------------------------
    {
        const double x99 = std::max(
            quantile(raw_relative, 0.99),
            quantile(final_relative, 0.99));
        const double xmax = std::max(5.0, std::min(100.0, std::ceil(x99 + 2.0)));

        TH1D h_raw("h_excl_note_raw", "", 40, 0.0, xmax);
        TH1D h_final("h_excl_note_final", "", 40, 0.0, xmax);

        for (double v : raw_relative) {
            if (v >= 0.0 && v < xmax) h_raw.Fill(v);
        }
        for (double v : final_relative) {
            if (v >= 0.0 && v < xmax) h_final.Fill(v);
        }

        if (h_raw.Integral() > 0.0) h_raw.Scale(1.0 / h_raw.Integral());
        if (h_final.Integral() > 0.0) h_final.Scale(1.0 / h_final.Integral());

        TCanvas c("c_excl_note_distribution", "", 1080, 760);
        style_note_canvas(c, 0.13);

        h_raw.SetLineColor(kGray + 2);
        h_raw.SetLineWidth(2);
        h_raw.SetLineStyle(2);
        h_raw.SetFillStyle(0);

        h_final.SetLineColor(kBlue + 1);
        h_final.SetLineWidth(3);
        h_final.SetFillStyle(0);

        h_final.SetMaximum(
            1.25 * std::max(h_raw.GetMaximum(), h_final.GetMaximum()));
        h_final.GetXaxis()->SetTitle("Relative exclusivity-selection systematic (%)");
        h_final.GetYaxis()->SetTitle("Fraction of populated bins");
        style_note_axes(h_final.GetXaxis(), h_final.GetYaxis());

        h_final.Draw("HIST");
        h_raw.Draw("HIST SAME");

        TLegend leg(0.57, 0.72, 0.94, 0.86);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.029);
        leg.AddEntry(&h_raw,
                     "Before statistical-consistency filtering", "l");
        leg.AddEntry(&h_final,
                     "Final assigned uncertainty", "l");
        leg.Draw();

        draw_note_title("Distribution of the exclusivity-selection systematic");

        c.SaveAs(
            (fs::path(outdir) /
             "exclusivity_systematic_distribution.png").string().c_str());
    }

    // ------------------------------------------------------------------
    // 2) xB dependence: median and central 68% spread across Q2, |t|, phi.
    // ------------------------------------------------------------------
    {
        struct Bucket {
            double lo = 0.0;
            double hi = 0.0;
            std::vector<double> values;
        };

        std::map<std::pair<double,double>, Bucket> buckets;

        for (const auto& r : rows) {
            if (!r.nominal.ok || std::fabs(r.nominal.value) <= 1e-30) continue;
            auto& b = buckets[{r.xb_min, r.xb_max}];
            b.lo = r.xb_min;
            b.hi = r.xb_max;
            b.values.push_back(
                100.0 * r.excl_final_abs / std::fabs(r.nominal.value));
        }

        TGraphAsymmErrors g;
        int ip = 0;
        double ymax = 0.0;

        for (const auto& kv : buckets) {
            const Bucket& b = kv.second;
            if (b.values.empty()) continue;

            const double x = 0.5 * (b.lo + b.hi);
            const double ex = 0.5 * (b.hi - b.lo);
            const double med = quantile(b.values, 0.50);
            const double p16 = quantile(b.values, 0.16);
            const double p84 = quantile(b.values, 0.84);

            g.SetPoint(ip, x, med);
            g.SetPointError(ip, ex, ex, med - p16, p84 - med);
            ymax = std::max(ymax, p84);
            ++ip;
        }

        TCanvas c("c_excl_note_xb", "", 1080, 760);
        style_note_canvas(c, 0.17);

        TH1D frame("h_excl_note_xb_frame", "", 100, 0.04, 0.62);
        frame.SetMinimum(0.0);
        frame.SetMaximum(std::max(5.0, 1.20 * ymax));
        frame.GetXaxis()->SetTitle("x_{B}");
        frame.GetYaxis()->SetTitle("Exclusivity-selection systematic (%)");
        style_note_axes(frame.GetXaxis(), frame.GetYaxis());
        frame.Draw();

        g.SetMarkerStyle(20);
        g.SetMarkerSize(1.05);
        g.SetMarkerColor(kBlue + 1);
        g.SetLineColor(kBlue + 1);
        g.SetLineWidth(2);
        g.Draw("PZ SAME");

        draw_note_title(
            "Kinematic dependence of the exclusivity-selection systematic",
            "Points: median over populated (Q^{2}, |t|, #phi) bins; bars: central 68% interval");

        c.SaveAs(
            (fs::path(outdir) /
             "exclusivity_systematic_vs_xB.png").string().c_str());
    }

    // ------------------------------------------------------------------
    // 3) Representative phi dependence of the final cross section under
    //    loose / nominal / tight exclusivity selections.
    //
    //    Prefer the representative cell used elsewhere in the note. If it
    //    is unavailable, choose the populated cell closest to its center.
    // ------------------------------------------------------------------
    {
        struct CellKey {
            double xb0, xb1, q0, q1, t0, t1;
            bool operator<(const CellKey& o) const {
                return std::tie(xb0, xb1, q0, q1, t0, t1) <
                       std::tie(o.xb0, o.xb1, o.q0, o.q1, o.t0, o.t1);
            }
        };

        std::map<CellKey, std::vector<const DiagnosticRow*>> cells;
        for (const auto& r : rows) {
            if (!r.nominal.ok || !r.excl_loose.ok || !r.excl_tight.ok) continue;
            cells[{r.xb_min, r.xb_max, r.q2_min, r.q2_max,
                   r.t_min, r.t_max}].push_back(&r);
        }

        const double target_x = 0.5 * (0.204 + 0.268);
        const double target_q = 0.5 * (2.510 + 3.295);
        const double target_t = 0.5 * (0.250 + 0.400);

        auto best = cells.end();
        double best_score = std::numeric_limits<double>::infinity();

        for (auto it = cells.begin(); it != cells.end(); ++it) {
            if (it->second.size() < 6) continue;
            const CellKey& k = it->first;
            const double xc = 0.5 * (k.xb0 + k.xb1);
            const double qc = 0.5 * (k.q0 + k.q1);
            const double tc = 0.5 * (k.t0 + k.t1);

            const double score =
                std::pow((xc - target_x) / 0.08, 2) +
                std::pow((qc - target_q) / 1.0, 2) +
                std::pow((tc - target_t) / 0.25, 2);

            if (score < best_score) {
                best_score = score;
                best = it;
            }
        }

        if (best != cells.end()) {
            auto selected = best->second;
            std::sort(selected.begin(), selected.end(),
                      [](const DiagnosticRow* a, const DiagnosticRow* b) {
                          return a->phi < b->phi;
                      });

            TGraphErrors g_nom;
            TGraphErrors g_loose;
            TGraphErrors g_tight;

            double ymin = std::numeric_limits<double>::infinity();
            double ymax = -std::numeric_limits<double>::infinity();

            for (int i = 0; i < static_cast<int>(selected.size()); ++i) {
                const auto* r = selected[i];
                g_nom.SetPoint(i, r->phi, r->nominal.value);
                g_nom.SetPointError(i, 0.0, r->nominal.stat);

                g_loose.SetPoint(i, r->phi, r->excl_loose.value);
                g_loose.SetPointError(i, 0.0, r->excl_loose.stat);

                g_tight.SetPoint(i, r->phi, r->excl_tight.value);
                g_tight.SetPointError(i, 0.0, r->excl_tight.stat);

                ymin = std::min(
                    ymin,
                    std::min(r->nominal.value,
                             std::min(r->excl_loose.value, r->excl_tight.value)));
                ymax = std::max(
                    ymax,
                    std::max(r->nominal.value,
                             std::max(r->excl_loose.value, r->excl_tight.value)));
            }

            const double span = std::max(1e-12, ymax - ymin);
            const double ylo = std::max(0.0, ymin - 0.15 * span);
            const double yhi = ymax + 0.28 * span;

            TCanvas c("c_excl_note_phi", "", 1080, 780);
            style_note_canvas(c, 0.18);

            TH1D frame("h_excl_note_phi_frame", "", 100, 0.0, 360.0);
            frame.SetMinimum(ylo);
            frame.SetMaximum(yhi);
            frame.GetXaxis()->SetTitle("#phi (deg)");
            frame.GetYaxis()->SetTitle("Cross section [nb/(GeV^{4} deg)]");
            style_note_axes(frame.GetXaxis(), frame.GetYaxis());
            frame.Draw();

            g_nom.SetMarkerStyle(20);
            g_nom.SetMarkerSize(0.95);
            g_nom.SetMarkerColor(kBlack);
            g_nom.SetLineColor(kBlack);
            g_nom.SetLineWidth(2);

            g_loose.SetMarkerStyle(24);
            g_loose.SetMarkerSize(0.95);
            g_loose.SetMarkerColor(kBlue + 1);
            g_loose.SetLineColor(kBlue + 1);
            g_loose.SetLineWidth(2);

            g_tight.SetMarkerStyle(25);
            g_tight.SetMarkerSize(0.95);
            g_tight.SetMarkerColor(kRed + 1);
            g_tight.SetLineColor(kRed + 1);
            g_tight.SetLineWidth(2);

            g_nom.Draw("PZ SAME");
            g_loose.Draw("PZ SAME");
            g_tight.Draw("PZ SAME");

            TLegend leg(0.43, 0.54, 0.82, 0.73);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(0.0225);
            leg.SetMargin(0.15);
            leg.AddEntry(&g_nom,
                         "Nominal selection (95% containment)", "pe");
            leg.AddEntry(&g_loose,
                         "Loose selection (98% containment)", "pe");
            leg.AddEntry(&g_tight,
                         "Tight selection (90% containment)", "pe");
            leg.Draw();

            const CellKey& k = best->first;
            const std::string subtitle =
                kinematic_subtitle(k.xb0, k.xb1, k.q0, k.q1, k.t0, k.t1);
            draw_note_title(
                "Representative effect of the exclusivity-selection variation",
                subtitle.c_str());

            c.SaveAs(
                (fs::path(outdir) /
                 "exclusivity_variation_phi_example.png").string().c_str());
        }
    }

    // ------------------------------------------------------------------
    // 4) Fraction of loose/tight variations that survive Barlow B >= 1.
    //    This is a diagnostic figure; it is useful for documenting how much
    //    of the observed cut dependence exceeds the statistical expectation.
    // ------------------------------------------------------------------
    {
        struct Counts {
            double lo = 0.0;
            double hi = 0.0;
            int total = 0;
            int loose = 0;
            int tight = 0;
        };

        std::map<std::pair<double,double>, Counts> buckets;

        for (const auto& r : rows) {
            if (!r.nominal.ok) continue;
            auto& b = buckets[{r.xb_min, r.xb_max}];
            b.lo = r.xb_min;
            b.hi = r.xb_max;
            ++b.total;
            if (r.excl_keep_loose) ++b.loose;
            if (r.excl_keep_tight) ++b.tight;
        }

        TGraphErrors g_loose;
        TGraphErrors g_tight;
        int ip = 0;

        for (const auto& kv : buckets) {
            const Counts& b = kv.second;
            if (b.total <= 0) continue;

            const double x = 0.5 * (b.lo + b.hi);
            const double ex = 0.5 * (b.hi - b.lo);
            const double fl = static_cast<double>(b.loose) / b.total;
            const double ft = static_cast<double>(b.tight) / b.total;

            g_loose.SetPoint(ip, x, fl);
            g_loose.SetPointError(ip, ex, 0.0);
            g_tight.SetPoint(ip, x, ft);
            g_tight.SetPointError(ip, ex, 0.0);
            ++ip;
        }

        TCanvas c("c_excl_note_barlow", "", 1080, 760);
        style_note_canvas(c, 0.16);

        TH1D frame("h_excl_note_barlow_frame", "", 100, 0.04, 0.62);
        frame.SetMinimum(0.0);
        frame.SetMaximum(1.08);
        frame.GetXaxis()->SetTitle("x_{B}");
        frame.GetYaxis()->SetTitle("Fraction of variations retained");
        style_note_axes(frame.GetXaxis(), frame.GetYaxis());
        frame.Draw();

        g_loose.SetMarkerStyle(24);
        g_loose.SetMarkerSize(1.0);
        g_loose.SetMarkerColor(kBlue + 1);
        g_loose.SetLineColor(kBlue + 1);
        g_loose.SetLineWidth(2);

        g_tight.SetMarkerStyle(25);
        g_tight.SetMarkerSize(1.0);
        g_tight.SetMarkerColor(kRed + 1);
        g_tight.SetLineColor(kRed + 1);
        g_tight.SetLineWidth(2);

        g_loose.Draw("PZ SAME");
        g_tight.Draw("PZ SAME");

        TLegend leg(0.59, 0.72, 0.94, 0.85);
        leg.SetBorderSize(0);
        leg.SetFillStyle(0);
        leg.SetTextFont(42);
        leg.SetTextSize(0.029);
        leg.AddEntry(&g_loose, "Loose selection (98%)", "pe");
        leg.AddEntry(&g_tight, "Tight selection (90%)", "pe");
        leg.Draw();

        draw_note_title(
            "Statistical-consistency filtering of the cut variations",
            "A variation is retained when B #geq 1");

        c.SaveAs(
            (fs::path(outdir) /
             "exclusivity_barlow_retention_vs_xB.png").string().c_str());
    }

    std::cout << "[cut-systematics] Wrote exclusivity analysis-note outputs to "
              << outdir << '\n';
}


void make_cut_kinematic_summary(const std::vector<DiagnosticRow>& rows,
                                const std::string& outdir,
                                bool fiducial) {
    struct AxisSpec {
        const char* title;
        double DiagnosticRow::*lo;
        double DiagnosticRow::*hi;
    };
    const std::array<AxisSpec,4> specs = {{
        {"x_{B}", &DiagnosticRow::xb_min, &DiagnosticRow::xb_max},
        {"Q^{2} (GeV^{2})", &DiagnosticRow::q2_min, &DiagnosticRow::q2_max},
        {"|t| (GeV^{2})", &DiagnosticRow::t_min, &DiagnosticRow::t_max},
        {"#phi (deg)", &DiagnosticRow::phi_min, &DiagnosticRow::phi_max}
    }};

    const std::string stem = fiducial ? "fiducial" : "exclusivity";
    const std::string what = fiducial ? "fiducial-selection" : "exclusivity-selection";

    TCanvas c(("c_" + stem + "_kinematic_summary").c_str(), "", 1450, 1050);
    c.Divide(2,2,0.002,0.002);

    for (int ia=0; ia<4; ++ia) {
        c.cd(ia+1);
        gPad->SetLeftMargin((ia%2==0) ? 0.14 : 0.12);
        gPad->SetRightMargin(0.035);
        gPad->SetBottomMargin((ia>=2) ? 0.15 : 0.12);
        gPad->SetTopMargin((ia<2) ? 0.16 : 0.10);
        gPad->SetTicks(1,1);

        struct Bucket { double lo=0, hi=0; std::vector<double> values; };
        std::map<std::pair<double,double>,Bucket> buckets;

        for (const auto& r : rows) {
            if (!r.nominal.ok || std::fabs(r.nominal.value)<=1e-30) continue;
            double lo = r.*(specs[ia].lo);
            double hi = r.*(specs[ia].hi);
            if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi>lo)) continue;
            if (ia == 3) {
                double phi = std::isfinite(r.phi) ? r.phi : 0.5*(lo+hi);
                while (phi < 0.0) phi += 360.0;
                while (phi >= 360.0) phi -= 360.0;
                const int ib = std::min(11, std::max(0, int(phi / 30.0)));
                lo = 30.0 * ib;
                hi = lo + 30.0;
            }
            auto& b = buckets[{lo,hi}]; b.lo=lo; b.hi=hi;
            const double u = fiducial ? r.fid_final_abs : r.excl_final_abs;
            b.values.push_back(100.0*u/std::fabs(r.nominal.value));
        }

        TGraphAsymmErrors g;
        int ip=0; double xmin=INFINITY, xmax=-INFINITY, ymax=0;
        for (const auto& kv : buckets) {
            const auto& b=kv.second; if (b.values.empty()) continue;
            const double x=.5*(b.lo+b.hi), ex=.5*(b.hi-b.lo);
            const double med=quantile(b.values,.50), p16=quantile(b.values,.16), p84=quantile(b.values,.84);
            g.SetPoint(ip,x,med); g.SetPointError(ip,ex,ex,med-p16,p84-med); ++ip;
            xmin=std::min(xmin,b.lo); xmax=std::max(xmax,b.hi); ymax=std::max(ymax,p84);
        }
        if (!(xmax>xmin)) { xmin=0; xmax=1; }

        TH1D frame(("h_"+stem+"_kinematic_"+std::to_string(ia)).c_str(),"",100,xmin,xmax);
        frame.SetMinimum(0); frame.SetMaximum(std::max(3.0,1.18*ymax));
        frame.GetXaxis()->SetTitle(specs[ia].title);
        frame.GetYaxis()->SetTitle("Systematic uncertainty (%)");
        frame.GetXaxis()->SetTitleSize(.050); frame.GetYaxis()->SetTitleSize(.048);
        frame.GetXaxis()->SetLabelSize(.043); frame.GetYaxis()->SetLabelSize(.042);
        frame.GetXaxis()->SetTitleOffset(1.08); frame.GetYaxis()->SetTitleOffset(1.16);
        // Draw pad-owned copies.  `frame` and `g` are local to this loop
        // iteration, while the canvas is saved only after all four pads have
        // been filled.  Drawing the original stack objects leaves dangling
        // primitive pointers after the iteration ends and produces blank pads.
        frame.DrawCopy();

        g.SetMarkerStyle(20);
        g.SetMarkerSize(.95);
        g.SetMarkerColor(kBlue+1);
        g.SetLineColor(kBlue+1);
        g.SetLineWidth(2);
        g.DrawClone("PZ SAME");

        // Keep panel labels inside the pad, clear of the canvas-level title.
        TLatex p;
        p.SetNDC();
        p.SetTextFont(42);
        p.SetTextSize(.036);
        const std::string lab=std::string("(")+char('a'+ia)+")";
        p.DrawLatex(.18,.80,lab.c_str());
    }

    c.cd(0);
    TLatex t;
    t.SetNDC();
    t.SetTextFont(42);
    t.SetTextAlign(22);
    t.SetTextSize(.021);
    const std::string title="Kinematic dependence of the "+what+" systematic";
    t.DrawLatex(.50,.975,title.c_str());
    t.SetTextSize(.015);
    t.DrawLatex(.50,.946,
                "Median and central 68% range; #phi projection grouped in 30^{#circ} intervals");
    c.SaveAs((fs::path(outdir)/(stem+"_systematic_kinematic_summary.png")).string().c_str());
}

void make_fiducial_analysis_note_plots(const std::vector<DiagnosticRow>& rows,
                                       const std::string& outdir) {
    fs::create_directories(outdir);
    std::vector<double> raw, fin;
    for (const auto& r:rows) {
        if (!r.nominal.ok || std::fabs(r.nominal.value)<=1e-30) continue;
        raw.push_back(100.0*r.fid_raw_abs/std::fabs(r.nominal.value));
        fin.push_back(100.0*r.fid_final_abs/std::fabs(r.nominal.value));
    }
    if (fin.empty()) return;

    { // summary CSV
        std::ofstream o(fs::path(outdir)/"fiducial_systematic_summary.csv");
        o<<"quantity,value\n";
        o<<"populated bins,"<<fin.size()<<"\n";
        o<<"mean final relative systematic (%),"<<mean_value(fin)<<"\n";
        o<<"median final relative systematic (%),"<<quantile(fin,.50)<<"\n";
        o<<"16th percentile final relative systematic (%),"<<quantile(fin,.16)<<"\n";
        o<<"84th percentile final relative systematic (%),"<<quantile(fin,.84)<<"\n";
        o<<"95th percentile final relative systematic (%),"<<quantile(fin,.95)<<"\n";
    }

    { // distribution
        const double xmax=std::max(5.0,std::min(100.0,std::ceil(std::max(quantile(raw,.99),quantile(fin,.99))+2.0)));
        TH1D hr("h_fid_note_raw","",40,0,xmax), hf("h_fid_note_final","",40,0,xmax);
        for(double v:raw) if(v>=0&&v<xmax) hr.Fill(v);
        for(double v:fin) if(v>=0&&v<xmax) hf.Fill(v);
        if(hr.Integral()>0)hr.Scale(1/hr.Integral()); if(hf.Integral()>0)hf.Scale(1/hf.Integral());
        TCanvas c("c_fid_note_dist","",1080,760); style_note_canvas(c,.13);
        hr.SetLineColor(kGray+2);hr.SetLineWidth(2);hr.SetLineStyle(2);hr.SetFillStyle(0);
        hf.SetLineColor(kBlue+1);hf.SetLineWidth(3);hf.SetFillStyle(0);hf.SetMaximum(1.25*std::max(hr.GetMaximum(),hf.GetMaximum()));
        hf.GetXaxis()->SetTitle("Relative fiducial-selection systematic (%)"); hf.GetYaxis()->SetTitle("Fraction of populated bins");
        style_note_axes(hf.GetXaxis(),hf.GetYaxis()); hf.Draw("HIST"); hr.Draw("HIST SAME");
        TLegend l(.55,.72,.94,.86);l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.028);
        l.AddEntry(&hr,"Before statistical-consistency filtering","l");l.AddEntry(&hf,"Final assigned uncertainty","l");l.Draw();
        draw_note_title("Distribution of the fiducial-selection systematic");
        c.SaveAs((fs::path(outdir)/"fiducial_systematic_distribution.png").string().c_str());
    }

    make_cut_kinematic_summary(rows,outdir,true);

    { // representative cross sections vs phi
        struct Key{double x0,x1,q0,q1,t0,t1;bool operator<(const Key&o)const{return std::tie(x0,x1,q0,q1,t0,t1)<std::tie(o.x0,o.x1,o.q0,o.q1,o.t0,o.t1);}};
        std::map<Key,std::vector<const DiagnosticRow*>> cells;
        for(const auto&r:rows)if(r.nominal.ok&&r.fid_loose.ok&&r.fid_tight.ok)cells[{r.xb_min,r.xb_max,r.q2_min,r.q2_max,r.t_min,r.t_max}].push_back(&r);
        const double tx=.5*(.204+.268),tq=.5*(2.510+3.295),tt=.5*(.250+.400);
        auto best=cells.end();double bs=INFINITY;
        for(auto it=cells.begin();it!=cells.end();++it){if(it->second.size()<6)continue;const auto&k=it->first;double s=pow((.5*(k.x0+k.x1)-tx)/.08,2)+pow((.5*(k.q0+k.q1)-tq)/1.0,2)+pow((.5*(k.t0+k.t1)-tt)/.25,2);if(s<bs){bs=s;best=it;}}
        if(best!=cells.end()){
            auto v=best->second;std::sort(v.begin(),v.end(),[](auto*a,auto*b){return a->phi<b->phi;});
            TGraphErrors gn,gl,gt;double ymin=INFINITY,ymax=-INFINITY;int i=0;
            for(auto*r:v){gn.SetPoint(i,r->phi,r->nominal.value);gn.SetPointError(i,0,r->nominal.stat);gl.SetPoint(i,r->phi,r->fid_loose.value);gl.SetPointError(i,0,r->fid_loose.stat);gt.SetPoint(i,r->phi,r->fid_tight.value);gt.SetPointError(i,0,r->fid_tight.stat);ymin=std::min({ymin,r->nominal.value,r->fid_loose.value,r->fid_tight.value});ymax=std::max({ymax,r->nominal.value,r->fid_loose.value,r->fid_tight.value});++i;}
            const double span=std::max(1e-12,ymax-ymin);TCanvas c("c_fid_note_phi","",1080,800);style_note_canvas(c,.20);
            TH1D fr("h_fid_note_phi_frame","",100,0,360);fr.SetMinimum(std::max(0.0,ymin-.15*span));fr.SetMaximum(ymax+.34*span);fr.GetXaxis()->SetTitle("#phi (deg)");fr.GetYaxis()->SetTitle("Cross section [nb/(GeV^{4} deg)]");style_note_axes(fr.GetXaxis(),fr.GetYaxis());fr.Draw();
            gn.SetMarkerStyle(20);gn.SetMarkerColor(kBlack);gn.SetLineColor(kBlack);gl.SetMarkerStyle(24);gl.SetMarkerColor(kBlue+1);gl.SetLineColor(kBlue+1);gt.SetMarkerStyle(25);gt.SetMarkerColor(kRed+1);gt.SetLineColor(kRed+1);gn.Draw("PZ SAME");gl.Draw("PZ SAME");gt.Draw("PZ SAME");
            TLegend l(.42,.54,.82,.73);
            l.SetBorderSize(0);
            l.SetFillStyle(0);
            l.SetTextFont(42);
            l.SetTextSize(.022);
            l.SetMargin(.15);
            l.AddEntry(&gn,"Nominal fiducial selection","pe");
            l.AddEntry(&gl,"Loose variation (#pm 2^{#circ} outward)","pe");
            l.AddEntry(&gt,"Tight variation (#pm 2^{#circ} inward)","pe");
            l.Draw();
            const auto&k=best->first;const std::string sub=kinematic_subtitle(k.x0,k.x1,k.q0,k.q1,k.t0,k.t1);draw_note_title("Representative effect of the fiducial-selection variation",sub.c_str());
            c.SaveAs((fs::path(outdir)/"fiducial_variation_phi_example.png").string().c_str());
        }
    }
}

void write_diagnostics(const std::vector<DiagnosticRow>& rows, const std::string& path) {
    fs::create_directories(fs::path(path).parent_path());
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Could not write diagnostics: " + path);
    out << "bin index,Bin Name,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,phiavg,"
        << "nominal xs,nominal stat,exclusivity loose xs,exclusivity loose stat,"
        << "exclusivity tight xs,exclusivity tight stat,exclusivity loose raw delta,"
        << "exclusivity tight raw delta,exclusivity loose Barlow,exclusivity tight Barlow,"
        << "exclusivity loose retained,exclusivity tight retained,"
        << "exclusivity tight relative difference,exclusivity loose-only prescription,"
        << "exclusivity raw systematic (nb/GeV^4),exclusivity final systematic (nb/GeV^4),"
        << "fiducial loose xs,fiducial loose stat,fiducial tight xs,fiducial tight stat,"
        << "fiducial loose raw delta,fiducial tight raw delta,fiducial loose Barlow,"
        << "fiducial tight Barlow,fiducial loose retained,fiducial tight retained,"
        << "fiducial tight relative difference,fiducial loose-only prescription,"
        << "fiducial raw systematic (nb/GeV^4),fiducial final systematic (nb/GeV^4)\n";
    for (const auto& r : rows) {
        out << r.bin_index << ',' << csv_quote(r.bin_name) << ','
            << r.xb_min << ',' << r.xb_max << ',' << r.q2_min << ',' << r.q2_max << ','
            << r.t_min << ',' << r.t_max << ',' << r.phi << ','
            << r.nominal.value << ',' << r.nominal.stat << ','
            << r.excl_loose.value << ',' << r.excl_loose.stat << ','
            << r.excl_tight.value << ',' << r.excl_tight.stat << ','
            << r.excl_delta_loose_raw << ',' << r.excl_delta_tight_raw << ','
            << r.excl_barlow_loose << ',' << r.excl_barlow_tight << ','
            << (r.excl_keep_loose ? 1 : 0) << ',' << (r.excl_keep_tight ? 1 : 0) << ','
            << r.excl_tight_relative_difference << ',' << (r.excl_use_loose_only ? 1 : 0) << ','
            << r.excl_raw_abs << ',' << r.excl_final_abs << ','
            << r.fid_loose.value << ',' << r.fid_loose.stat << ','
            << r.fid_tight.value << ',' << r.fid_tight.stat << ','
            << r.fid_delta_loose_raw << ',' << r.fid_delta_tight_raw << ','
            << r.fid_barlow_loose << ',' << r.fid_barlow_tight << ','
            << (r.fid_keep_loose ? 1 : 0) << ',' << (r.fid_keep_tight ? 1 : 0) << ','
            << r.fid_tight_relative_difference << ',' << (r.fid_use_loose_only ? 1 : 0) << ','
            << r.fid_raw_abs << ',' << r.fid_final_abs << '\n';
    }
}

void make_plots(const std::vector<DiagnosticRow>& rows, const std::string& outdir) {
    fs::create_directories(outdir);
    gStyle->SetOptStat(0);

    std::map<std::pair<double,double>, std::vector<const DiagnosticRow*>> groups;
    for (const auto& r : rows) {
        if (!r.nominal.ok) continue;
        groups[{r.xb_min, r.xb_max}].push_back(&r);
    }

    for (auto& kv : groups) {
        auto& v = kv.second;
        std::sort(v.begin(), v.end(), [](const DiagnosticRow* a, const DiagnosticRow* b) {
            if (a->q2_min != b->q2_min) return a->q2_min < b->q2_min;
            if (a->t_min != b->t_min) return a->t_min < b->t_min;
            return a->phi < b->phi;
        });

        const int n = static_cast<int>(v.size());
        if (n == 0) continue;
        std::vector<double> x(n), zero(n,0.0), exRaw(n), exFinal(n), fiRaw(n), fiFinal(n);
        std::vector<double> bExL(n), bExT(n), bFiL(n), bFiT(n);
        std::vector<double> ratioExL(n), ratioExT(n), ratioFiL(n), ratioFiT(n);
        for (int i = 0; i < n; ++i) {
            x[i] = i + 1;
            exRaw[i] = v[i]->excl_raw_abs;
            exFinal[i] = v[i]->excl_final_abs;
            fiRaw[i] = v[i]->fid_raw_abs;
            fiFinal[i] = v[i]->fid_final_abs;
            bExL[i] = v[i]->excl_barlow_loose;
            bExT[i] = v[i]->excl_barlow_tight;
            bFiL[i] = v[i]->fid_barlow_loose;
            bFiT[i] = v[i]->fid_barlow_tight;
            const double nom = v[i]->nominal.value;
            ratioExL[i] = (v[i]->excl_loose.ok && std::fabs(nom)>1e-30) ? v[i]->excl_loose.value/nom : 0.0;
            ratioExT[i] = (v[i]->excl_tight.ok && std::fabs(nom)>1e-30) ? v[i]->excl_tight.value/nom : 0.0;
            ratioFiL[i] = (v[i]->fid_loose.ok && std::fabs(nom)>1e-30) ? v[i]->fid_loose.value/nom : 0.0;
            ratioFiT[i] = (v[i]->fid_tight.ok && std::fabs(nom)>1e-30) ? v[i]->fid_tight.value/nom : 0.0;
        }

        const std::string stem = [&]() {
            std::ostringstream ss;
            ss << "xB_" << std::fixed << std::setprecision(2) << kv.first.first
               << "_" << kv.first.second;
            std::string s = ss.str();
            std::replace(s.begin(), s.end(), '.', 'p');
            return s;
        }();

        {
            TCanvas c("c_sys", "", 1800, 1100);
            c.Divide(1,2);
            c.cd(1);
            TGraph g1(n, x.data(), exRaw.data());
            TGraph g2(n, x.data(), exFinal.data());
            g1.SetMarkerStyle(24); g1.SetLineWidth(2);
            g2.SetMarkerStyle(20); g2.SetLineWidth(2);
            g1.SetTitle(("Exclusivity-cut systematic, " + xbin_label(kv.first.first, kv.first.second)).c_str());
            g1.GetXaxis()->SetTitle("Sequential bin number within x_{B} range");
            g1.GetYaxis()->SetTitle("Point-to-point systematic (nb/GeV^{4})");
            g1.Draw("APL"); g2.Draw("PL SAME");
            TLegend l1(0.72,0.76,0.94,0.90); l1.AddEntry(&g1,"Raw prescription","lp"); l1.AddEntry(&g2,"After Barlow","lp"); l1.Draw();
            c.cd(2);
            TGraph f1(n, x.data(), fiRaw.data());
            TGraph f2(n, x.data(), fiFinal.data());
            f1.SetMarkerStyle(24); f1.SetLineWidth(2);
            f2.SetMarkerStyle(20); f2.SetLineWidth(2);
            f1.SetTitle(("Auxiliary-fiducial-cut systematic, " + xbin_label(kv.first.first, kv.first.second)).c_str());
            f1.GetXaxis()->SetTitle("Sequential bin number within x_{B} range");
            f1.GetYaxis()->SetTitle("Point-to-point systematic (nb/GeV^{4})");
            f1.Draw("APL"); f2.Draw("PL SAME");
            TLegend l2(0.72,0.76,0.94,0.90); l2.AddEntry(&f1,"Raw prescription","lp"); l2.AddEntry(&f2,"After Barlow","lp"); l2.Draw();
            c.SaveAs((fs::path(outdir)/(stem+"_systematics.png")).string().c_str());
        }

        {
            TCanvas c("c_barlow", "", 1800, 1100);
            c.Divide(1,2);
            c.cd(1);
            TGraph a(n,x.data(),bExL.data()), b(n,x.data(),bExT.data());
            a.SetMarkerStyle(24); b.SetMarkerStyle(20); a.SetLineWidth(2); b.SetLineWidth(2);
            a.SetTitle(("Exclusivity Barlow values, " + xbin_label(kv.first.first, kv.first.second)).c_str());
            a.GetXaxis()->SetTitle("Sequential bin number within x_{B} range"); a.GetYaxis()->SetTitle("Barlow B");
            a.Draw("APL"); b.Draw("PL SAME");
            TLine line1(1.0,1.0,n,1.0); line1.SetLineStyle(2); line1.Draw();
            TLegend l1(0.70,0.76,0.94,0.90); l1.AddEntry(&a,"95% variation","lp"); l1.AddEntry(&b,"99.99% variation","lp"); l1.Draw();
            c.cd(2);
            TGraph d(n,x.data(),bFiL.data()), e(n,x.data(),bFiT.data());
            d.SetMarkerStyle(24); e.SetMarkerStyle(20); d.SetLineWidth(2); e.SetLineWidth(2);
            d.SetTitle(("Auxiliary-fiducial Barlow values, " + xbin_label(kv.first.first, kv.first.second)).c_str());
            d.GetXaxis()->SetTitle("Sequential bin number within x_{B} range"); d.GetYaxis()->SetTitle("Barlow B");
            d.Draw("APL"); e.Draw("PL SAME");
            TLine line2(1.0,1.0,n,1.0); line2.SetLineStyle(2); line2.Draw();
            TLegend l2(0.70,0.76,0.94,0.90); l2.AddEntry(&d,"Loose angles","lp"); l2.AddEntry(&e,"Tight angles","lp"); l2.Draw();
            c.SaveAs((fs::path(outdir)/(stem+"_barlow.png")).string().c_str());
        }

        {
            TCanvas c("c_ratios", "", 1800, 1100);
            c.Divide(1,2);
            c.cd(1);
            TGraph a(n,x.data(),ratioExL.data()), b(n,x.data(),ratioExT.data());
            a.SetMarkerStyle(24); b.SetMarkerStyle(20); a.SetLineWidth(2); b.SetLineWidth(2);
            a.SetTitle(("Exclusivity variation / nominal, " + xbin_label(kv.first.first, kv.first.second)).c_str());
            a.GetXaxis()->SetTitle("Sequential bin number within x_{B} range"); a.GetYaxis()->SetTitle("#sigma_{var}/#sigma_{nom}");
            a.Draw("APL"); b.Draw("PL SAME");
            TLine line1(1.0,1.0,n,1.0); line1.SetLineStyle(2); line1.Draw();
            TLegend l1(0.70,0.76,0.94,0.90); l1.AddEntry(&a,"95% variation","lp"); l1.AddEntry(&b,"99.99% variation","lp"); l1.Draw();
            c.cd(2);
            TGraph d(n,x.data(),ratioFiL.data()), e(n,x.data(),ratioFiT.data());
            d.SetMarkerStyle(24); e.SetMarkerStyle(20); d.SetLineWidth(2); e.SetLineWidth(2);
            d.SetTitle(("Auxiliary-fiducial variation / nominal, " + xbin_label(kv.first.first, kv.first.second)).c_str());
            d.GetXaxis()->SetTitle("Sequential bin number within x_{B} range"); d.GetYaxis()->SetTitle("#sigma_{var}/#sigma_{nom}");
            d.Draw("APL"); e.Draw("PL SAME");
            TLine line2(1.0,1.0,n,1.0); line2.SetLineStyle(2); line2.Draw();
            TLegend l2(0.70,0.76,0.94,0.90); l2.AddEntry(&d,"Loose angles","lp"); l2.AddEntry(&e,"Tight angles","lp"); l2.Draw();
            c.SaveAs((fs::path(outdir)/(stem+"_cross_section_ratios.png")).string().c_str());
        }
    }
}

} // namespace

bool update_cut_variation_systematics(const CutVariationSystematicsOptions& options) {
    if (!options.enabled) {
        std::cout << "[cut-systematics] Disabled.\n";
        return true;
    }

    try {
        CsvTable nominal = read_csv(options.nominal_csv);
        const CsvTable ex_loose = read_csv(options.exclusivity_loose_csv);
        const CsvTable ex_tight = read_csv(options.exclusivity_tight_csv);
        const CsvTable fi_loose = read_csv(options.fiducial_loose_csv);
        const CsvTable fi_tight = read_csv(options.fiducial_tight_csv);

        const size_t c_nom = require_col(nominal, options.cross_section_column, options.nominal_csv);
        const size_t c_exl = require_col(ex_loose, options.cross_section_column, options.exclusivity_loose_csv);
        const size_t c_ext = require_col(ex_tight, options.cross_section_column, options.exclusivity_tight_csv);
        const size_t c_fil = require_col(fi_loose, options.cross_section_column, options.fiducial_loose_csv);
        const size_t c_fit = require_col(fi_tight, options.cross_section_column, options.fiducial_tight_csv);

        const auto map_exl = build_row_map(ex_loose);
        const auto map_ext = build_row_map(ex_tight);
        const auto map_fil = build_row_map(fi_loose);
        const auto map_fit = build_row_map(fi_tight);

        // Validate all variation inputs before adding columns or calculating systematics
        // values or writing a backup. Missing/malformed tuples must never be
        // interpreted as zero differences or universal Barlow rejection.
        require_valid_variation(
            "exclusivity loose",
            validate_variation_cross_sections(nominal, c_nom, ex_loose, c_exl, map_exl));
        require_valid_variation(
            "exclusivity tight",
            validate_variation_cross_sections(nominal, c_nom, ex_tight, c_ext, map_ext));
        require_valid_variation(
            "fiducial loose",
            validate_variation_cross_sections(nominal, c_nom, fi_loose, c_fil, map_fil));
        require_valid_variation(
            "fiducial tight",
            validate_variation_cross_sections(nominal, c_nom, fi_tight, c_fit, map_fit));

        const size_t c_ex_raw = ensure_col(nominal, "Syst. err (exclusivity cuts, raw)");
        const size_t c_ex_fin = ensure_col(nominal, "Syst. err (exclusivity cuts)");
        const size_t c_fi_raw = ensure_col(nominal, "Syst. err (fiducial cuts, raw)");
        const size_t c_fi_fin = ensure_col(nominal, "Syst. err (fiducial cuts)");
        const size_t c_pi0 = ensure_col(nominal, "Syst. err (pi0 subtraction)");
        const size_t c_acc = ensure_col(nominal, "Syst. err (Acceptance)");
        const size_t c_frad = ensure_col(nominal, "Syst.err (Frad)");
        const size_t c_fbin = ensure_col(nominal, "Syst.err (Fbin)");
        const size_t c_ptp = ensure_col(nominal, "Syst. err (point-to-point total)");

        const auto col_optional = [&](const std::string& name) -> int {
            const auto it = nominal.index.find(name);
            return it == nominal.index.end() ? -1 : static_cast<int>(it->second);
        };
        const int c_bin = col_optional("bin index");
        const int c_name = col_optional("Bin Name");
        const int c_xlo = col_optional("xBmin");
        const int c_xhi = col_optional("xBmax");
        const int c_qlo = col_optional("Q2min");
        const int c_qhi = col_optional("Q2max");
        const int c_tlo = col_optional("t_abs_min");
        const int c_thi = col_optional("t_abs_max");
        const int c_plo = col_optional("phimin");
        const int c_phi_hi = col_optional("phimax");
        int c_phi = col_optional("phiavg, 10.6 GeV");
        if (c_phi < 0) c_phi = col_optional("phiavg");

        std::vector<DiagnosticRow> diagnostics;
        diagnostics.reserve(nominal.rows.size());

        for (size_t i = 0; i < nominal.rows.size(); ++i) {
            auto& row = nominal.rows[i];
            const std::string key = row_key(nominal, row);
            const auto a = map_exl.find(key), b = map_ext.find(key), c = map_fil.find(key), d = map_fit.find(key);
            if (a == map_exl.end() || b == map_ext.end() || c == map_fil.end() || d == map_fit.end()) {
                throw std::runtime_error("Variation CSV row mismatch for key '" + key + "'.");
            }

            DiagnosticRow r;
            r.bin_index = c_bin >= 0 ? static_cast<int>(to_double(row[c_bin], static_cast<double>(i + 1))) : static_cast<int>(i + 1);
            r.bin_name = c_name >= 0 ? row[c_name] : key;
            r.xb_min = c_xlo >= 0 ? to_double(row[c_xlo]) : 0.0;
            r.xb_max = c_xhi >= 0 ? to_double(row[c_xhi]) : 0.0;
            r.q2_min = c_qlo >= 0 ? to_double(row[c_qlo]) : 0.0;
            r.q2_max = c_qhi >= 0 ? to_double(row[c_qhi]) : 0.0;
            r.t_min = c_tlo >= 0 ? to_double(row[c_tlo]) : 0.0;
            r.t_max = c_thi >= 0 ? to_double(row[c_thi]) : 0.0;
            r.phi_min = c_plo >= 0 ? to_double(row[c_plo]) : 0.0;
            r.phi_max = c_phi_hi >= 0 ? to_double(row[c_phi_hi]) : 0.0;
            r.phi = c_phi >= 0 ? to_double(row[c_phi]) : 0.5*(r.phi_min+r.phi_max);
            if (!std::isfinite(r.phi)) r.phi = 0.5*(r.phi_min+r.phi_max);

            r.nominal = parse_triple(row[c_nom]);
            r.excl_loose = parse_triple(ex_loose.rows[a->second][c_exl]);
            r.excl_tight = parse_triple(ex_tight.rows[b->second][c_ext]);
            r.fid_loose = parse_triple(fi_loose.rows[c->second][c_fil]);
            r.fid_tight = parse_triple(fi_tight.rows[d->second][c_fit]);

            if (r.nominal.ok) {
                // The preflight validation guarantees that all four varied
                // tuples are readable for every nominally populated bin.
                r.excl_delta_loose_raw = r.excl_loose.value - r.nominal.value;
                r.excl_delta_tight_raw = r.excl_tight.value - r.nominal.value;
                r.excl_barlow_loose = barlow_value(r.excl_loose, r.nominal);
                r.excl_barlow_tight = barlow_value(r.excl_tight, r.nominal);
                r.excl_keep_loose = !options.apply_barlow || r.excl_barlow_loose >= options.barlow_threshold;
                r.excl_keep_tight = !options.apply_barlow || r.excl_barlow_tight >= options.barlow_threshold;
                r.excl_tight_relative_difference =
                    relative_difference(r.excl_delta_tight_raw, r.nominal.value);
                r.excl_use_loose_only = options.use_pass1_tight_instability_rule &&
                    r.excl_tight_relative_difference > options.tight_relative_difference_threshold;
                r.excl_raw_abs = pass1_cut_systematic_absolute(
                    r.excl_delta_loose_raw,
                    r.excl_delta_tight_raw,
                    r.excl_use_loose_only);
                r.excl_final_abs = pass1_cut_systematic_absolute(
                    r.excl_keep_loose ? r.excl_delta_loose_raw : 0.0,
                    r.excl_keep_tight ? r.excl_delta_tight_raw : 0.0,
                    r.excl_use_loose_only);

                r.fid_delta_loose_raw = r.fid_loose.value - r.nominal.value;
                r.fid_delta_tight_raw = r.fid_tight.value - r.nominal.value;
                r.fid_barlow_loose = barlow_value(r.fid_loose, r.nominal);
                r.fid_barlow_tight = barlow_value(r.fid_tight, r.nominal);
                r.fid_keep_loose = !options.apply_barlow || r.fid_barlow_loose >= options.barlow_threshold;
                r.fid_keep_tight = !options.apply_barlow || r.fid_barlow_tight >= options.barlow_threshold;
                r.fid_tight_relative_difference =
                    relative_difference(r.fid_delta_tight_raw, r.nominal.value);
                r.fid_use_loose_only = options.use_pass1_tight_instability_rule &&
                    r.fid_tight_relative_difference > options.tight_relative_difference_threshold;
                r.fid_raw_abs = pass1_cut_systematic_absolute(
                    r.fid_delta_loose_raw,
                    r.fid_delta_tight_raw,
                    r.fid_use_loose_only);
                r.fid_final_abs = pass1_cut_systematic_absolute(
                    r.fid_keep_loose ? r.fid_delta_loose_raw : 0.0,
                    r.fid_keep_tight ? r.fid_delta_tight_raw : 0.0,
                    r.fid_use_loose_only);
            }

            row[c_ex_raw] = format_number(r.excl_raw_abs);
            row[c_ex_fin] = format_number(r.excl_final_abs);
            row[c_fi_raw] = format_number(r.fid_raw_abs);
            row[c_fi_fin] = format_number(r.fid_final_abs);

            const double s_pi0 = to_double(row[c_pi0], 0.0);
            const double s_acc = to_double(row[c_acc], 0.0);
            const double s_frad = to_double(row[c_frad], 0.0);
            const double s_fbin = to_double(row[c_fbin], 0.0);
            const double s_ptp = std::sqrt(s_pi0*s_pi0 + s_acc*s_acc +
                                           s_frad*s_frad + s_fbin*s_fbin +
                                           r.excl_final_abs*r.excl_final_abs +
                                           r.fid_final_abs*r.fid_final_abs);
            row[c_ptp] = format_number(s_ptp);
            diagnostics.push_back(r);
        }

        const std::string backup = options.nominal_csv + ".before_cut_systematics";
        fs::copy_file(options.nominal_csv, backup, fs::copy_options::overwrite_existing);
        write_csv(nominal, options.nominal_csv);

        if (options.write_diagnostic_csv) {
            write_diagnostics(diagnostics, (fs::path(options.output_dir)/"cut_variation_diagnostics.csv").string());
        }
        if (options.make_plots) {
            make_plots(diagnostics, (fs::path(options.output_dir)/"plots").string());
        }
        if (options.make_analysis_note_plots) {
            const std::string excl_note_dir =
                (fs::path(options.output_dir) / "analysis_note/exclusivity_selection").string();
            const std::string fid_note_dir =
                (fs::path(options.output_dir) / "analysis_note/fiducial_selection").string();

            make_analysis_note_plots(diagnostics, excl_note_dir);
            make_cut_kinematic_summary(diagnostics, excl_note_dir, false);
            make_fiducial_analysis_note_plots(diagnostics, fid_note_dir);
        }

        const auto excl_loose_only_count = std::count_if(
            diagnostics.begin(), diagnostics.end(),
            [](const DiagnosticRow& r) { return r.excl_use_loose_only; });
        const auto fid_loose_only_count = std::count_if(
            diagnostics.begin(), diagnostics.end(),
            [](const DiagnosticRow& r) { return r.fid_use_loose_only; });
        std::cout << "[cut-systematics] Pass-1 tight-instability rule (threshold "
                  << 100.0 * options.tight_relative_difference_threshold
                  << "%): exclusivity loose-only in " << excl_loose_only_count
                  << " bins; fiducial loose-only in " << fid_loose_only_count << " bins.\n";

        std::cout << "[cut-systematics] Updated " << options.nominal_csv << " with raw and Barlow-filtered "
                  << "exclusivity/fiducial point-to-point systematics for " << diagnostics.size() << " bins.\n";
        return true;
    } catch (const std::exception& e) {
        std::cerr << "[cut-systematics] FATAL: " << e.what() << '\n';
        return false;
    }
}
