// norm_cross_sections.cpp
// -----------------------------------------------------------------------------
// Apply overall-normalization factors to already-computed DVCS cross sections
// and make the same plots as cross_sections.cpp, but using the normalized values.
//
// Convention (from overall_normalization_study.cpp):
//   "norm, <Label>" is defined as (xs / BH)
//
// Therefore BH-normalized cross sections must DIVIDE:
//   sigma_normed = sigma_raw / N
//
// Where:
//   N is read from:   "norm, <Label>"   (scalar or tuple3)
//   sigma_raw from:  "cross sections, ep->epg, exp, <Label>, <hel>"
//   output written:  "normed cross sections, ep->epg, exp, <Label>, <hel>"
//
// Error propagation (independent):
//   y = x / N
//   stat_y = sqrt( (stat_x / N)^2 + (x * stat_N / N^2)^2 )
//   sys_y  = sqrt( (sys_x  / N)^2 + (x * sys_N  / N^2)^2 )
//
// Plotting overlays BH / KM / VGG curves loaded from:
//   theory_json_root/<EnergyDir>/xs_phi_all.json
// keyed by CSV row index (row number in the CSV, header excluded).
//
// Additional ratio plotting:
//   For each xB bin, make separate ratio canvases for:
//     - unpolarized data
//     - positive-helicity data
//     - negative-helicity data
//
//   Each ratio subplot overlays data/model ratios with respect to:
//     - KM15
//     - BH
//     - VGG
//
// Ratio convention:
//   ratio = data / model
//   ratio_stat = stat_data / model
//
// The model uncertainty is not included. The ratio uncertainty is propagated
// only from the statistical uncertainty on the normalized data point.
// -----------------------------------------------------------------------------

#include "norm_cross_sections.h"

#include <algorithm>
#include <cctype>
#include <cmath>
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

#include <nlohmann/json.hpp>

// ROOT
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TString.h>
#include <TH1.h>
#include <TSystem.h>
#include <TLine.h>

namespace fs = std::filesystem;
using json   = nlohmann::json;

using Range = std::pair<double, double>;

namespace {

struct Triple {
    double value;
    double stat;
    double sys;
};

// -----------------------------------------------------------------------------
// Label helpers (match cross_sections.cpp behavior / your conventions)
// -----------------------------------------------------------------------------

static std::string canonical_period_dir(const std::string &label) {
    if (label == "Fa18 Inb")      return "Fa18_Inb";
    if (label == "Fa18 Out")      return "Fa18_Out";
    if (label == "Fa18 Inb Supp") return "Fa18_Inb_Supp";
    if (label == "Sp18 Inb")      return "Sp18_Inb";
    if (label == "Sp18 Out")      return "Sp18_Out";
    if (label == "Sp19 Inb")      return "Sp19_Inb";
    if (label == "Fa18")          return "Fa18";
    if (label == "Sp18")          return "Sp18";
    if (label == "10.6 GeV")      return "10.6_GeV";
    if (label == "10.2 GeV")      return "10.2_GeV";

    std::string out = label;
    std::replace(out.begin(), out.end(), ' ', '_');
    return out;
}

static std::string theory_energy_label_for(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") return "10.2 GeV";
    return "10.6 GeV";
}

// -----------------------------------------------------------------------------
// CSV helpers (robust for quoted fields with escaped quotes)
// -----------------------------------------------------------------------------

static std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;

    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;

    return s.substr(b, e - b);
}

static std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (in_quotes) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                field.push_back(c);
            }
        } else {
            if (c == '"') {
                in_quotes = true;
            } else if (c == ',') {
                out.push_back(field);
                field.clear();
            } else {
                field.push_back(c);
            }
        }
    }

    out.push_back(field);
    return out;
}

static std::string quote_if_needed(const std::string &s) {
    bool need = false;

    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            need = true;
            break;
        }
    }

    if (!need) return s;

    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');

    for (char c : s) {
        if (c == '"') {
            out.push_back('"');
            out.push_back('"');
        } else {
            out.push_back(c);
        }
    }

    out.push_back('"');
    return out;
}

static std::string join_csv_line(const std::vector<std::string> &fields) {
    std::ostringstream oss;

    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) oss << ",";
        oss << quote_if_needed(fields[i]);
    }

    return oss.str();
}

static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(header[i]) == target) {
            return (int)i;
        }
    }

    return -1;
}

static int find_col_required(const std::vector<std::string> &header,
                             const std::string &target) {
    const int idx = find_col_optional(header, target);

    if (idx < 0) {
        throw std::runtime_error("Missing required column: \"" + target + "\"");
    }

    return idx;
}

// -----------------------------------------------------------------------------
// Tuple parsing/formatting
// -----------------------------------------------------------------------------

static Triple parse_tuple3(const std::string &cell_raw) {
    Triple out{0.0, 0.0, 0.0};

    std::string s = trim(cell_raw);
    if (s.empty()) return out;

    if (!s.empty() && s.front() == '(' && s.back() == ')') {
        s = trim(s.substr(1, s.size() - 2));
        if (s.empty()) return out;
    }

    std::vector<std::string> parts;
    std::string token;

    for (size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];

        if (c == ',') {
            parts.push_back(trim(token));
            token.clear();
        } else {
            token.push_back(c);
        }
    }

    if (!token.empty() || s.find(',') != std::string::npos) {
        parts.push_back(trim(token));
    }

    auto to_double_or_zero = [](const std::string &str) -> double {
        if (str.empty()) return 0.0;
        return std::atof(str.c_str());
    };

    if (!parts.empty())    out.value = to_double_or_zero(parts[0]);
    if (parts.size() > 1U) out.stat  = to_double_or_zero(parts[1]);
    if (parts.size() > 2U) out.sys   = to_double_or_zero(parts[2]);

    return out;
}

static bool cell_looks_like_tuple3(const std::string &cell_raw) {
    const std::string s = trim(cell_raw);

    if (s.empty()) return false;
    if (s.front() == '(' && s.back() == ')') return true;
    if (s.find(',') != std::string::npos) return true;

    return false;
}

static bool parse_norm_cell(const std::string &cell_raw, Triple *outN) {
    if (!outN) return false;

    *outN = Triple{0.0, 0.0, 0.0};

    const std::string s = trim(cell_raw);
    if (s.empty()) return false;

    if (cell_looks_like_tuple3(s)) {
        const Triple t = parse_tuple3(s);

        if (!(t.value > 0.0) || !std::isfinite(t.value)) return false;

        *outN = t;
        return true;
    }

    const double N = std::atof(s.c_str());

    if (!(N > 0.0) || !std::isfinite(N)) return false;

    outN->value = N;
    outN->stat  = 0.0;
    outN->sys   = 0.0;

    return true;
}

static std::string tuple3_to_cell(double value, double stat, double sys) {
    std::ostringstream oss;

    oss << "("
        << std::setprecision(10) << value << ", "
        << std::setprecision(10) << stat  << ", "
        << std::setprecision(10) << sys   << ")";

    return oss.str();
}

static void ensure_dir(const fs::path &p) {
    std::error_code ec;

    if (fs::exists(p, ec)) return;

    fs::create_directories(p, ec);
}

// -----------------------------------------------------------------------------
// Theory loading (same JSON structure as cross_sections.cpp generator)
// -----------------------------------------------------------------------------

struct TheoryCurves {
    std::vector<double> phi_deg;

    std::vector<double> bh_unpol;
    std::vector<double> bh_pos;
    std::vector<double> bh_neg;

    std::vector<double> km_unpol;
    std::vector<double> km_pos;
    std::vector<double> km_neg;

    std::vector<double> vgg_unpol;
    std::vector<double> vgg_pos;
    std::vector<double> vgg_neg;
};

static std::map<size_t, TheoryCurves>
load_theory_for_label(const std::string &label,
                      const std::string &theory_root) {
    std::map<size_t, TheoryCurves> out;

    const std::string energy_label = theory_energy_label_for(label);
    const fs::path dir  = fs::path(theory_root) / canonical_period_dir(energy_label);
    const fs::path file = dir / "xs_phi_all.json";

    if (!fs::exists(file)) {
        std::cerr << "[norm_cross_sections] WARNING: no theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    std::ifstream ifs(file);

    if (!ifs) {
        std::cerr << "[norm_cross_sections] WARNING: cannot open theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    json j;

    try {
        ifs >> j;
    } catch (...) {
        std::cerr << "[norm_cross_sections] WARNING: malformed theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    const std::vector<double> phi_deg = j.value("phi_deg", std::vector<double>{});

    if (phi_deg.empty()) {
        std::cerr << "[norm_cross_sections] WARNING: theory JSON has empty phi_deg at "
                  << file.string() << "\n";
        return out;
    }

    if (!j.contains("rows") || !j["rows"].is_object()) {
        std::cerr << "[norm_cross_sections] WARNING: theory JSON has no rows object at "
                  << file.string() << "\n";
        return out;
    }

    for (auto it = j["rows"].begin(); it != j["rows"].end(); ++it) {
        const std::string row_key = it.key();
        const json &cell = it.value();

        size_t row_index = 0;

        try {
            row_index = (size_t)std::stoul(row_key);
        } catch (...) {
            continue;
        }

        TheoryCurves tc;
        tc.phi_deg = phi_deg;

        tc.bh_unpol  = cell["BH"].value("unpol", std::vector<double>{});
        tc.bh_pos    = cell["BH"].value("pos",   std::vector<double>{});
        tc.bh_neg    = cell["BH"].value("neg",   std::vector<double>{});

        tc.km_unpol  = cell["KM"].value("unpol", std::vector<double>{});
        tc.km_pos    = cell["KM"].value("pos",   std::vector<double>{});
        tc.km_neg    = cell["KM"].value("neg",   std::vector<double>{});

        tc.vgg_unpol = cell["VGG"].value("unpol", std::vector<double>{});
        tc.vgg_pos   = cell["VGG"].value("pos",   std::vector<double>{});
        tc.vgg_neg   = cell["VGG"].value("neg",   std::vector<double>{});

        if (!tc.phi_deg.empty()) {
            out[row_index] = std::move(tc);
        }
    }

    std::cout << "[norm_cross_sections] Loaded theory for label \"" << label
              << "\" (energy " << energy_label << ") rows=" << out.size()
              << " from " << file.string() << "\n";

    return out;
}

// -----------------------------------------------------------------------------
// Normalization math
// -----------------------------------------------------------------------------

static void normalize_one_sigma_cell(const Triple &sigma_raw,
                                     const Triple &N,
                                     Triple *sigma_out) {
    if (!sigma_out) return;

    *sigma_out = Triple{0.0, 0.0, 0.0};

    if (!(sigma_raw.value > 0.0) || !(N.value > 0.0)) return;

    const double x  = sigma_raw.value;
    const double Nx = N.value;

    const double y = x / Nx;

    const double stat_x = sigma_raw.stat;
    const double sys_x  = sigma_raw.sys;
    const double stat_N = N.stat;
    const double sys_N  = N.sys;

    const double stat_y = std::sqrt(
        (stat_x / Nx) * (stat_x / Nx) +
        (x * stat_N / (Nx * Nx)) * (x * stat_N / (Nx * Nx))
    );

    const double sys_y = std::sqrt(
        (sys_x / Nx) * (sys_x / Nx) +
        (x * sys_N / (Nx * Nx)) * (x * sys_N / (Nx * Nx))
    );

    sigma_out->value = y;
    sigma_out->stat  = stat_y;
    sigma_out->sys   = sys_y;
}

// -----------------------------------------------------------------------------
// Plotting structures (mirrors cross_sections style)
// -----------------------------------------------------------------------------

struct Point {
    double phi;
    double xs;
    double xs_err; // using stat as error bars (matches cross_sections.cpp style)
};

struct RatioPoint {
    double phi;
    double ratio;
    double ratio_err;
};

struct BinData {
    std::vector<Point> unpol;
    std::vector<Point> pos;
    std::vector<Point> neg;

    size_t theory_row = 0;
    bool   have_theory_row = false;
};

using QTKey = std::pair<Range, Range>; // (Q2, |t|)

struct XSGroupByXB {
    std::map<QTKey, BinData> bins;
    int xb_index = -1;
};

enum class XSecPanelMode {
    All,
    UnpolOnly,
    PosOnly,
    NegOnly
};

enum class RatioModel {
    KM15,
    BH,
    VGG
};

static std::string ratio_model_name(RatioModel model) {
    if (model == RatioModel::KM15) return "KM15";
    if (model == RatioModel::BH)   return "BH";
    return "VGG";
}

static int ratio_model_color(RatioModel model) {
    if (model == RatioModel::KM15) return kMagenta + 1;
    if (model == RatioModel::BH)   return kGreen + 2;
    return kOrange + 7;
}

static std::string ratio_mode_prefix(XSecPanelMode mode) {
    if (mode == XSecPanelMode::UnpolOnly) return "unpol";
    if (mode == XSecPanelMode::PosOnly)   return "pos";
    if (mode == XSecPanelMode::NegOnly)   return "neg";
    return "all";
}

static std::string ratio_mode_label(XSecPanelMode mode) {
    if (mode == XSecPanelMode::UnpolOnly) return "unpolarized";
    if (mode == XSecPanelMode::PosOnly)   return "+ helicity";
    if (mode == XSecPanelMode::NegOnly)   return "- helicity";
    return "all helicities";
}

static const std::vector<Point> *points_for_mode(const BinData &bin,
                                                 XSecPanelMode mode) {
    if (mode == XSecPanelMode::UnpolOnly) return &(bin.unpol);
    if (mode == XSecPanelMode::PosOnly)   return &(bin.pos);
    if (mode == XSecPanelMode::NegOnly)   return &(bin.neg);

    return nullptr;
}

static std::pair<double, double> compute_yrange_for_bin(
    const BinData *bin,
    const std::map<size_t, TheoryCurves> &theory,
    XSecPanelMode mode) {

    double ymin = std::numeric_limits<double>::max();
    double ymax = 0.0;

    auto update_from_points = [&](const std::vector<Point> &v) {
        for (const auto &p : v) {
            if (p.xs > 0.0) {
                const double ylow  = std::max(1e-12, p.xs - p.xs_err);
                const double yhigh = p.xs + p.xs_err;

                if (ylow  > 0.0 && ylow  < ymin) ymin = ylow;
                if (yhigh > ymax) ymax = yhigh;
            }
        }
    };

    auto update_from_curve = [&](const std::vector<double> &ys) {
        for (double y : ys) {
            if (y <= 0.0) continue;
            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }
    };

    if (bin) {
        if (mode == XSecPanelMode::All || mode == XSecPanelMode::UnpolOnly) {
            update_from_points(bin->unpol);
        }

        if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly) {
            update_from_points(bin->pos);
        }

        if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly) {
            update_from_points(bin->neg);
        }

        if (bin->have_theory_row) {
            const auto it_th = theory.find(bin->theory_row);

            if (it_th != theory.end()) {
                const TheoryCurves &tc = it_th->second;

                if (mode == XSecPanelMode::All) {
                    update_from_curve(tc.bh_unpol);

                    update_from_curve(tc.km_unpol);
                    update_from_curve(tc.km_pos);
                    update_from_curve(tc.km_neg);

                    update_from_curve(tc.vgg_unpol);
                    update_from_curve(tc.vgg_pos);
                    update_from_curve(tc.vgg_neg);
                } else if (mode == XSecPanelMode::UnpolOnly) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_unpol);
                    update_from_curve(tc.vgg_unpol);
                } else if (mode == XSecPanelMode::PosOnly) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_pos);
                    update_from_curve(tc.vgg_pos);
                } else if (mode == XSecPanelMode::NegOnly) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_neg);
                    update_from_curve(tc.vgg_neg);
                }
            }
        }
    }

    if (ymax <= 0.0 || !std::isfinite(ymax)) {
        ymin = 1e-4;
        ymax = 1.0;
    } else {
        if (ymin <= 0.0 || !std::isfinite(ymin)) {
            ymin = ymax * 1e-3;
        }

        const double logmin = std::pow(10.0, std::floor(std::log10(ymin)));
        const double logmax = std::pow(10.0, std::ceil (std::log10(ymax)));

        ymin = std::max(1e-4, logmin);
        ymax = logmax;
    }

    return std::make_pair(ymin, 1.2 * ymax);
}

static TGraphErrors *make_xsec_graph(const std::vector<Point> &v,
                                     int mstyle,
                                     int mcolor) {
    if (v.empty()) return nullptr;

    const int N = (int)v.size();

    std::vector<double> x(N);
    std::vector<double> y(N);
    std::vector<double> ex(N);
    std::vector<double> ey(N);

    for (int i = 0; i < N; ++i) {
        x[i]  = v[i].phi;
        y[i]  = v[i].xs;
        ex[i] = 0.0;
        ey[i] = v[i].xs_err;
    }

    TGraphErrors *g = new TGraphErrors(N, x.data(), y.data(), ex.data(), ey.data());

    g->SetMarkerStyle(mstyle);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
    g->SetLineColor(mcolor);
    g->SetMarkerColor(mcolor);

    return g;
}

static TGraphErrors *make_ratio_graph(const std::vector<RatioPoint> &v,
                                      int mstyle,
                                      int mcolor) {
    if (v.empty()) return nullptr;

    const int N = (int)v.size();

    std::vector<double> x(N);
    std::vector<double> y(N);
    std::vector<double> ex(N);
    std::vector<double> ey(N);

    for (int i = 0; i < N; ++i) {
        x[i]  = v[i].phi;
        y[i]  = v[i].ratio;
        ex[i] = 0.0;
        ey[i] = v[i].ratio_err;
    }

    TGraphErrors *g = new TGraphErrors(N, x.data(), y.data(), ex.data(), ey.data());

    g->SetMarkerStyle(mstyle);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
    g->SetLineColor(mcolor);
    g->SetMarkerColor(mcolor);

    return g;
}

static double interp_curve_at_phi(const std::vector<double> &x,
                                  const std::vector<double> &y,
                                  double phi) {
    if (x.empty() || y.empty()) return 0.0;
    if (x.size() != y.size()) return 0.0;

    std::vector<std::pair<double, double> > xy;
    xy.reserve(x.size());

    for (size_t i = 0; i < x.size(); ++i) {
        if (std::isfinite(x[i]) && std::isfinite(y[i]) && y[i] > 0.0) {
            xy.push_back(std::make_pair(x[i], y[i]));
        }
    }

    if (xy.empty()) return 0.0;

    std::sort(xy.begin(), xy.end(),
              [](const std::pair<double, double> &a,
                 const std::pair<double, double> &b) {
                  return a.first < b.first;
              });

    if (xy.size() == 1U) {
        return xy.front().second;
    }

    if (phi <= xy.front().first) {
        return xy.front().second;
    }

    if (phi >= xy.back().first) {
        return xy.back().second;
    }

    for (size_t i = 1; i < xy.size(); ++i) {
        const double x0 = xy[i - 1].first;
        const double x1 = xy[i].first;

        if (phi >= x0 && phi <= x1) {
            const double y0 = xy[i - 1].second;
            const double y1 = xy[i].second;

            if (x1 == x0) return y0;

            const double f = (phi - x0) / (x1 - x0);
            return y0 + f * (y1 - y0);
        }
    }

    return 0.0;
}

static const std::vector<double> *select_model_curve(const TheoryCurves &tc,
                                                     RatioModel model,
                                                     XSecPanelMode hel_mode) {
    if (model == RatioModel::BH) {
        return &(tc.bh_unpol);
    }

    if (model == RatioModel::KM15) {
        if (hel_mode == XSecPanelMode::UnpolOnly) return &(tc.km_unpol);
        if (hel_mode == XSecPanelMode::PosOnly)   return &(tc.km_pos);
        if (hel_mode == XSecPanelMode::NegOnly)   return &(tc.km_neg);

        return nullptr;
    }

    if (model == RatioModel::VGG) {
        if (hel_mode == XSecPanelMode::UnpolOnly) return &(tc.vgg_unpol);
        if (hel_mode == XSecPanelMode::PosOnly)   return &(tc.vgg_pos);
        if (hel_mode == XSecPanelMode::NegOnly)   return &(tc.vgg_neg);

        return nullptr;
    }

    return nullptr;
}

static std::vector<RatioPoint> build_ratios_for_points(
    const std::vector<Point> &points,
    const TheoryCurves &tc,
    RatioModel model,
    XSecPanelMode hel_mode) {

    std::vector<RatioPoint> ratios;

    const std::vector<double> *curve = select_model_curve(tc, model, hel_mode);

    if (!curve) return ratios;
    if (curve->empty()) return ratios;
    if (tc.phi_deg.size() != curve->size()) return ratios;

    for (const auto &p : points) {
        if (!(p.xs > 0.0) || !std::isfinite(p.xs)) continue;
        if (!(p.xs_err >= 0.0) || !std::isfinite(p.xs_err)) continue;

        const double model_value = interp_curve_at_phi(tc.phi_deg, *curve, p.phi);

        if (!(model_value > 0.0) || !std::isfinite(model_value)) continue;

        RatioPoint rp;
        rp.phi       = p.phi;
        rp.ratio     = p.xs / model_value;
        rp.ratio_err = p.xs_err / model_value;

        if (!std::isfinite(rp.ratio)) continue;
        if (!std::isfinite(rp.ratio_err)) continue;

        ratios.push_back(rp);
    }

    std::sort(ratios.begin(), ratios.end(),
              [](const RatioPoint &a, const RatioPoint &b) {
                  return a.phi < b.phi;
              });

    return ratios;
}

static std::pair<double, double> compute_ratio_yrange_for_bin_and_mode(
    const BinData *bin,
    const std::map<size_t, TheoryCurves> &theory,
    XSecPanelMode mode) {

    double ymin = std::numeric_limits<double>::max();
    double ymax = 0.0;

    if (bin && bin->have_theory_row) {
        const auto it_th = theory.find(bin->theory_row);

        if (it_th != theory.end()) {
            const TheoryCurves &tc = it_th->second;
            const std::vector<Point> *points = points_for_mode(*bin, mode);

            if (points) {
                const std::vector<RatioModel> models = {
                    RatioModel::KM15,
                    RatioModel::BH,
                    RatioModel::VGG
                };

                for (RatioModel model : models) {
                    const std::vector<RatioPoint> ratios =
                        build_ratios_for_points(*points, tc, model, mode);

                    for (const auto &p : ratios) {
                        if (!std::isfinite(p.ratio)) continue;
                        if (!std::isfinite(p.ratio_err)) continue;

                        const double lo = std::max(0.0, p.ratio - p.ratio_err);
                        const double hi = p.ratio + p.ratio_err;

                        if (lo < ymin) ymin = lo;
                        if (hi > ymax) ymax = hi;
                    }
                }
            }
        }
    }

    if (!(ymax > 0.0) || !std::isfinite(ymax)) {
        return std::make_pair(0.0, 2.0);
    }

    if (!std::isfinite(ymin) || ymin == std::numeric_limits<double>::max()) {
        ymin = 0.0;
    }

    ymin = std::max(0.0, ymin * 0.80);
    ymax = ymax * 1.20;

    if (ymax < 2.0) {
        ymax = 2.0;
    }

    return std::make_pair(ymin, ymax);
}

// -----------------------------------------------------------------------------
// Canvas builder (one xB bin, one cross-section view mode)
// -----------------------------------------------------------------------------

static void make_normed_xsec_canvas_for_mode(
    const std::string &label,
    const Range &xb_range,
    const XSGroupByXB &group,
    const std::vector<Range> &q2_slice,
    const std::vector<Range> &t_slice,
    const std::map<size_t, TheoryCurves> &theory,
    const fs::path &outdir,
    int xb_idx_for_name,
    XSecPanelMode mode,
    int ncols,
    int nrows,
    int nPads) {

    (void)nPads;

    const auto &bins_for_xB = group.bins;

    if (bins_for_xB.empty()) return;

    const int W = 300 * ncols + 160;
    const int H = 260 * nrows + 240;

    TCanvas *c = new TCanvas("c_normed_xsec", "c_normed_xsec", W, H);

    TPad *pTop = new TPad("pTop", "pTop", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad *pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.78);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.06);

    std::ostringstream tit;
    tit << "Normed cross sections, ep #rightarrow ep#gamma   " << label
        << "   x_{B} in ("
        << std::fixed << std::setprecision(3)
        << xb_range.first << ", " << xb_range.second << ")";

    if      (mode == XSecPanelMode::UnpolOnly) tit << "   (unpolarized only)";
    else if (mode == XSecPanelMode::PosOnly)   tit << "   (+ helicity only)";
    else if (mode == XSecPanelMode::NegOnly)   tit << "   (- helicity only)";

    head.DrawLatex(0.5, 0.86, tit.str().c_str());

    TGraphErrors dummy_unpol;
    TGraphErrors dummy_pos;
    TGraphErrors dummy_neg;

    dummy_unpol.SetMarkerStyle(20);
    dummy_unpol.SetMarkerSize(1.0);
    dummy_unpol.SetLineWidth(2);
    dummy_unpol.SetMarkerColor(kBlack);
    dummy_unpol.SetLineColor(kBlack);

    dummy_pos.SetMarkerStyle(24);
    dummy_pos.SetMarkerSize(1.0);
    dummy_pos.SetLineWidth(2);
    dummy_pos.SetMarkerColor(kRed + 1);
    dummy_pos.SetLineColor(kRed + 1);

    dummy_neg.SetMarkerStyle(25);
    dummy_neg.SetMarkerSize(1.0);
    dummy_neg.SetLineWidth(2);
    dummy_neg.SetMarkerColor(kBlue + 1);
    dummy_neg.SetLineColor(kBlue + 1);

    TGraph dummy_bh;
    TGraph dummy_km_unpol;
    TGraph dummy_km_pos;
    TGraph dummy_km_neg;
    TGraph dummy_vgg_unpol;
    TGraph dummy_vgg_pos;
    TGraph dummy_vgg_neg;

    dummy_bh.SetLineWidth(2);
    dummy_bh.SetLineStyle(2);
    dummy_bh.SetLineColor(kGreen + 2);

    dummy_km_unpol.SetLineWidth(2);
    dummy_km_unpol.SetLineStyle(1);
    dummy_km_unpol.SetLineColor(kMagenta + 1);

    dummy_km_pos.SetLineWidth(2);
    dummy_km_pos.SetLineStyle(2);
    dummy_km_pos.SetLineColor(kMagenta + 1);

    dummy_km_neg.SetLineWidth(2);
    dummy_km_neg.SetLineStyle(3);
    dummy_km_neg.SetLineColor(kMagenta + 1);

    dummy_vgg_unpol.SetLineWidth(2);
    dummy_vgg_unpol.SetLineStyle(1);
    dummy_vgg_unpol.SetLineColor(kOrange + 7);

    dummy_vgg_pos.SetLineWidth(2);
    dummy_vgg_pos.SetLineStyle(2);
    dummy_vgg_pos.SetLineColor(kOrange + 7);

    dummy_vgg_neg.SetLineWidth(2);
    dummy_vgg_neg.SetLineStyle(3);
    dummy_vgg_neg.SetLineColor(kOrange + 7);

    TLegend *legData = new TLegend(0.02, 0.10, 0.32, 0.75);
    legData->SetBorderSize(1);
    legData->SetLineColor(kBlack);
    legData->SetFillColor(kWhite);
    legData->SetFillStyle(1001);
    legData->SetTextFont(42);
    legData->SetTextSize(0.045);

    TLegend *legKM = new TLegend(0.35, 0.10, 0.65, 0.75);
    legKM->SetBorderSize(1);
    legKM->SetLineColor(kBlack);
    legKM->SetFillColor(kWhite);
    legKM->SetFillStyle(1001);
    legKM->SetTextFont(42);
    legKM->SetTextSize(0.045);

    TLegend *legVGG = new TLegend(0.68, 0.10, 0.98, 0.75);
    legVGG->SetBorderSize(1);
    legVGG->SetLineColor(kBlack);
    legVGG->SetFillColor(kWhite);
    legVGG->SetFillStyle(1001);
    legVGG->SetTextFont(42);
    legVGG->SetTextSize(0.045);

    if (mode == XSecPanelMode::All) {
        legData->AddEntry(&dummy_unpol, "data unpolarized (normed)", "lep");
        legData->AddEntry(&dummy_pos,   "data + helicity (normed)", "lep");
        legData->AddEntry(&dummy_neg,   "data - helicity (normed)", "lep");

        legKM->AddEntry(&dummy_bh,       "BH unpolarized", "l");
        legKM->AddEntry(&dummy_km_unpol, "KM unpolarized", "l");
        legKM->AddEntry(&dummy_km_pos,   "KM + helicity",  "l");
        legKM->AddEntry(&dummy_km_neg,   "KM - helicity",  "l");

        legVGG->AddEntry(&dummy_vgg_unpol, "VGG unpolarized", "l");
        legVGG->AddEntry(&dummy_vgg_pos,   "VGG + helicity",  "l");
        legVGG->AddEntry(&dummy_vgg_neg,   "VGG - helicity",  "l");
    } else if (mode == XSecPanelMode::UnpolOnly) {
        legData->AddEntry(&dummy_unpol, "data unpolarized (normed)", "lep");

        legKM->AddEntry(&dummy_bh,       "BH unpolarized", "l");
        legKM->AddEntry(&dummy_km_unpol, "KM unpolarized", "l");

        legVGG->AddEntry(&dummy_vgg_unpol, "VGG unpolarized", "l");
    } else if (mode == XSecPanelMode::PosOnly) {
        legData->AddEntry(&dummy_pos, "data + helicity (normed)", "lep");

        legKM->AddEntry(&dummy_bh,     "BH",            "l");
        legKM->AddEntry(&dummy_km_pos, "KM + helicity", "l");

        legVGG->AddEntry(&dummy_vgg_pos, "VGG + helicity", "l");
    } else if (mode == XSecPanelMode::NegOnly) {
        legData->AddEntry(&dummy_neg, "data - helicity (normed)", "lep");

        legKM->AddEntry(&dummy_bh,     "BH",            "l");
        legKM->AddEntry(&dummy_km_neg, "KM - helicity", "l");

        legVGG->AddEntry(&dummy_vgg_neg, "VGG - helicity", "l");
    }

    legData->Draw();
    legKM->Draw();
    legVGG->Draw();

    for (int r = 0; r < nrows; ++r) {
        const Range &t_range = t_slice[r];

        for (int cc = 0; cc < ncols; ++cc) {
            const Range &q2_range = q2_slice[cc];

            pGrid->cd(r * ncols + cc + 1);

            gPad->SetGrid(1, 1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.07);
            gPad->SetLogy(true);

            QTKey key(q2_range, t_range);

            auto it_bin = bins_for_xB.find(key);

            const BinData *bin_ptr = nullptr;
            if (it_bin != bins_for_xB.end()) {
                bin_ptr = &(it_bin->second);
            }

            double ymin_canvas = 1e-4;
            double ymax_canvas = 1.0;

            {
                auto yr = compute_yrange_for_bin(bin_ptr, theory, mode);
                ymin_canvas = yr.first;
                ymax_canvas = yr.second;
            }

            TH1 *frame = gPad->DrawFrame(0.0, ymin_canvas, 360.0, ymax_canvas);

            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("d^{4}#sigma_{norm} / (dx_{B} dQ^{2} d|t| d#phi)");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.070);
            frame->GetYaxis()->SetTitleSize(0.070);
            frame->GetXaxis()->SetLabelSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.060);
            frame->GetXaxis()->SetTitleOffset(1.00);
            frame->GetYaxis()->SetTitleOffset(1.15);

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.060);
            lab.SetTextAlign(11);
            lab.SetTextFont(42);

            lab.DrawLatex(
                0.12, 0.83,
                Form("Q^{2}=(%.2f, %.2f)  |t|=(%.2f, %.2f)",
                     q2_range.first, q2_range.second,
                     t_range.first,  t_range.second)
            );

            if (!bin_ptr) continue;

            BinData bin = *bin_ptr;

            auto sort_by_phi = [](std::vector<Point> &v) {
                std::sort(v.begin(), v.end(),
                          [](const Point &a, const Point &b) {
                              return a.phi < b.phi;
                          });
            };

            sort_by_phi(bin.unpol);
            sort_by_phi(bin.pos);
            sort_by_phi(bin.neg);

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::UnpolOnly) {
                TGraphErrors *g_unpol = make_xsec_graph(bin.unpol, 20, kBlack);
                if (g_unpol) g_unpol->Draw("P SAME");
            }

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly) {
                TGraphErrors *g_pos = make_xsec_graph(bin.pos, 24, kRed + 1);
                if (g_pos) g_pos->Draw("P SAME");
            }

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly) {
                TGraphErrors *g_neg = make_xsec_graph(bin.neg, 25, kBlue + 1);
                if (g_neg) g_neg->Draw("P SAME");
            }

            if (bin.have_theory_row) {
                auto it_th = theory.find(bin.theory_row);

                if (it_th != theory.end()) {
                    const TheoryCurves &tc = it_th->second;

                    auto draw_curve = [&](const std::vector<double> &ys,
                                          int lstyle,
                                          int lcolor) {
                        if (tc.phi_deg.size() != ys.size() || ys.empty()) return;

                        const int M = (int)ys.size();

                        std::vector<double> xp(M);
                        std::vector<double> yp(M);

                        for (int i = 0; i < M; ++i) {
                            xp[i] = tc.phi_deg[i];
                            yp[i] = ys[i];
                        }

                        TGraph *gth = new TGraph(M, xp.data(), yp.data());

                        gth->SetLineStyle(lstyle);
                        gth->SetLineWidth(2);
                        gth->SetLineColor(lcolor);
                        gth->Draw("L SAME");
                    };

                    if (mode == XSecPanelMode::All) {
                        draw_curve(tc.bh_unpol,  2, kGreen + 2);

                        draw_curve(tc.km_unpol,  1, kMagenta + 1);
                        draw_curve(tc.km_pos,    2, kMagenta + 1);
                        draw_curve(tc.km_neg,    3, kMagenta + 1);

                        draw_curve(tc.vgg_unpol, 1, kOrange + 7);
                        draw_curve(tc.vgg_pos,   2, kOrange + 7);
                        draw_curve(tc.vgg_neg,   3, kOrange + 7);
                    } else if (mode == XSecPanelMode::UnpolOnly) {
                        draw_curve(tc.bh_unpol,  2, kGreen + 2);
                        draw_curve(tc.km_unpol,  1, kMagenta + 1);
                        draw_curve(tc.vgg_unpol, 1, kOrange + 7);
                    } else if (mode == XSecPanelMode::PosOnly) {
                        draw_curve(tc.bh_unpol,  2, kGreen + 2);
                        draw_curve(tc.km_pos,    2, kMagenta + 1);
                        draw_curve(tc.vgg_pos,   2, kOrange + 7);
                    } else if (mode == XSecPanelMode::NegOnly) {
                        draw_curve(tc.bh_unpol,  2, kGreen + 2);
                        draw_curve(tc.km_neg,    3, kMagenta + 1);
                        draw_curve(tc.vgg_neg,   3, kOrange + 7);
                    }
                }
            }
        }
    }

    std::ostringstream fname;

    fname << "normed_cross_sections_";

    if      (mode == XSecPanelMode::UnpolOnly) fname << "unpol_";
    else if (mode == XSecPanelMode::PosOnly)   fname << "pos_";
    else if (mode == XSecPanelMode::NegOnly)   fname << "neg_";

    fname << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    fs::path outpath = outdir / fname.str();

    c->SaveAs(outpath.string().c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Canvas builder (one xB bin, one helicity mode, combined KM15/BH/VGG ratios)
// -----------------------------------------------------------------------------

static void make_ratio_canvas_for_mode(
    const std::string &label,
    const Range &xb_range,
    const XSGroupByXB &group,
    const std::vector<Range> &q2_slice,
    const std::vector<Range> &t_slice,
    const std::map<size_t, TheoryCurves> &theory,
    const fs::path &outdir,
    int xb_idx_for_name,
    XSecPanelMode mode,
    int ncols,
    int nrows,
    int nPads) {

    (void)nPads;

    if (mode == XSecPanelMode::All) return;

    const auto &bins_for_xB = group.bins;

    if (bins_for_xB.empty()) return;

    const int W = 300 * ncols + 160;
    const int H = 260 * nrows + 240;

    TCanvas *c = new TCanvas("c_normed_xsec_ratios", "c_normed_xsec_ratios", W, H);

    TPad *pTop = new TPad("pTopRatio", "pTopRatio", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad *pGrid = new TPad("pGridRatio", "pGridRatio", 0.0, 0.00, 1.0, 0.78);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.06);

    std::ostringstream tit;

    tit << "Normed cross-section data/model ratios, ep #rightarrow ep#gamma   " << label
        << "   " << ratio_mode_label(mode)
        << "   x_{B} in ("
        << std::fixed << std::setprecision(3)
        << xb_range.first << ", " << xb_range.second << ")";

    head.DrawLatex(0.5, 0.86, tit.str().c_str());

    TGraphErrors dummy_km15;
    TGraphErrors dummy_bh;
    TGraphErrors dummy_vgg;

    dummy_km15.SetMarkerStyle(20);
    dummy_km15.SetMarkerSize(1.0);
    dummy_km15.SetMarkerColor(ratio_model_color(RatioModel::KM15));
    dummy_km15.SetLineColor(ratio_model_color(RatioModel::KM15));
    dummy_km15.SetLineWidth(2);

    dummy_bh.SetMarkerStyle(20);
    dummy_bh.SetMarkerSize(1.0);
    dummy_bh.SetMarkerColor(ratio_model_color(RatioModel::BH));
    dummy_bh.SetLineColor(ratio_model_color(RatioModel::BH));
    dummy_bh.SetLineWidth(2);

    dummy_vgg.SetMarkerStyle(20);
    dummy_vgg.SetMarkerSize(1.0);
    dummy_vgg.SetMarkerColor(ratio_model_color(RatioModel::VGG));
    dummy_vgg.SetLineColor(ratio_model_color(RatioModel::VGG));
    dummy_vgg.SetLineWidth(2);

    TLegend *legModel = new TLegend(0.34, 0.14, 0.66, 0.72);
    legModel->SetBorderSize(1);
    legModel->SetLineColor(kBlack);
    legModel->SetFillColor(kWhite);
    legModel->SetFillStyle(1001);
    legModel->SetTextFont(42);
    legModel->SetTextSize(0.045);
    legModel->AddEntry(&dummy_km15, "data / KM15", "lep");
    legModel->AddEntry(&dummy_bh,   "data / BH",   "lep");
    legModel->AddEntry(&dummy_vgg,  "data / VGG",  "lep");
    legModel->Draw();

    const std::vector<RatioModel> models = {
        RatioModel::KM15,
        RatioModel::BH,
        RatioModel::VGG
    };

    for (int r = 0; r < nrows; ++r) {
        const Range &t_range = t_slice[r];

        for (int cc = 0; cc < ncols; ++cc) {
            const Range &q2_range = q2_slice[cc];

            pGrid->cd(r * ncols + cc + 1);

            gPad->SetGrid(1, 1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.07);
            gPad->SetLogy(false);

            QTKey key(q2_range, t_range);

            auto it_bin = bins_for_xB.find(key);

            const BinData *bin_ptr = nullptr;
            if (it_bin != bins_for_xB.end()) {
                bin_ptr = &(it_bin->second);
            }

            auto yr = compute_ratio_yrange_for_bin_and_mode(bin_ptr, theory, mode);
            const double ymin_canvas = yr.first;
            const double ymax_canvas = yr.second;

            TH1 *frame = gPad->DrawFrame(0.0, ymin_canvas, 360.0, ymax_canvas);

            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("data / model");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.070);
            frame->GetYaxis()->SetTitleSize(0.070);
            frame->GetXaxis()->SetLabelSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.060);
            frame->GetXaxis()->SetTitleOffset(1.00);
            frame->GetYaxis()->SetTitleOffset(1.15);

            TLine *unity = new TLine(0.0, 1.0, 360.0, 1.0);
            unity->SetLineColor(kBlack);
            unity->SetLineStyle(1);
            unity->SetLineWidth(2);
            unity->Draw("SAME");

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.060);
            lab.SetTextAlign(11);
            lab.SetTextFont(42);

            lab.DrawLatex(
                0.12, 0.83,
                Form("Q^{2}=(%.2f, %.2f)  |t|=(%.2f, %.2f)",
                     q2_range.first, q2_range.second,
                     t_range.first,  t_range.second)
            );

            if (!bin_ptr) continue;
            if (!bin_ptr->have_theory_row) continue;

            const auto it_th = theory.find(bin_ptr->theory_row);

            if (it_th == theory.end()) continue;

            const TheoryCurves &tc = it_th->second;
            const std::vector<Point> *points = points_for_mode(*bin_ptr, mode);

            if (!points) continue;

            for (RatioModel model : models) {
                const int color = ratio_model_color(model);

                const std::vector<RatioPoint> ratios =
                    build_ratios_for_points(*points, tc, model, mode);

                TGraphErrors *g = make_ratio_graph(ratios, 20, color);

                if (g) g->Draw("P SAME");
            }
        }
    }

    std::ostringstream fname;

    fname << "normed_cross_sections_ratios_"
          << ratio_mode_prefix(mode) << "_"
          << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    fs::path outpath = outdir / fname.str();

    c->SaveAs(outpath.string().c_str());

    delete c;
}

} // end anonymous namespace

// -----------------------------------------------------------------------------
// Public API: wrapper overload (fixes your original main.cpp call)
// -----------------------------------------------------------------------------

bool update_normed_cross_sections_csv(const std::string &csv_main) {
    const std::vector<std::string> labels_to_update = {
        "Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out", "Sp19 Inb",
        "Fa18", "Sp18", "10.6 GeV"
    };

    return update_normed_cross_sections_csv(csv_main, labels_to_update);
}

// -----------------------------------------------------------------------------
// Update normalized cross sections in the CSV
// -----------------------------------------------------------------------------

bool update_normed_cross_sections_csv(const std::string &csv_main,
                                      const std::vector<std::string> &labels_to_update) {
    std::ifstream ifs(csv_main);

    if (!ifs) {
        std::cerr << "[norm_cross_sections] ERROR: cannot open " << csv_main
                  << " for reading.\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }

    ifs.close();

    if (lines.empty()) {
        std::cerr << "[norm_cross_sections] ERROR: CSV " << csv_main
                  << " is empty.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    struct LabelCols {
        int c_norm = -1;

        int c_raw_unpol = -1;
        int c_raw_pos   = -1;
        int c_raw_neg   = -1;

        int c_out_unpol = -1;
        int c_out_pos   = -1;
        int c_out_neg   = -1;
    };

    std::map<std::string, LabelCols> cols;

    try {
        for (const auto &label : labels_to_update) {
            LabelCols lc;

            lc.c_norm = find_col_required(header, "norm, " + label);

            lc.c_raw_unpol = find_col_required(header, "cross sections, ep->epg, exp, " + label + ", unpol");
            lc.c_raw_pos   = find_col_required(header, "cross sections, ep->epg, exp, " + label + ", pos");
            lc.c_raw_neg   = find_col_required(header, "cross sections, ep->epg, exp, " + label + ", neg");

            lc.c_out_unpol = find_col_required(header, "normed cross sections, ep->epg, exp, " + label + ", unpol");
            lc.c_out_pos   = find_col_required(header, "normed cross sections, ep->epg, exp, " + label + ", pos");
            lc.c_out_neg   = find_col_required(header, "normed cross sections, ep->epg, exp, " + label + ", neg");

            cols[label] = lc;
        }
    } catch (const std::exception &e) {
        std::cerr << "[norm_cross_sections] FATAL: " << e.what() << "\n";
        return false;
    }

    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);

    const size_t n_data_rows = (lines.size() > 0 ? lines.size() - 1 : 0);

    std::cout << "[norm_cross_sections] update_normed_cross_sections_csv: data rows = "
              << n_data_rows << "\n";

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            out_lines.push_back(lines[row]);
            continue;
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);

        if (fields.size() != header.size()) {
            std::ostringstream oss;

            oss << "Row " << row << " has " << fields.size()
                << " columns, expected " << header.size();

            std::cerr << "[norm_cross_sections] FATAL: " << oss.str() << "\n";

            return false;
        }

        for (const auto &kv : cols) {
            const LabelCols &lc = kv.second;

            Triple N;

            if (!parse_norm_cell(fields[lc.c_norm], &N)) {
                continue;
            }

            auto do_one = [&](int c_raw, int c_out) {
                const std::string raw_cell = trim(fields[c_raw]);

                if (raw_cell.empty()) return;

                const Triple sigma_raw = parse_tuple3(raw_cell);

                if (!(sigma_raw.value > 0.0) || !std::isfinite(sigma_raw.value)) return;

                Triple sigma_norm;

                normalize_one_sigma_cell(sigma_raw, N, &sigma_norm);

                if (!(sigma_norm.value > 0.0) || !std::isfinite(sigma_norm.value)) return;

                fields[c_out] = tuple3_to_cell(sigma_norm.value, sigma_norm.stat, sigma_norm.sys);
            };

            do_one(lc.c_raw_unpol, lc.c_out_unpol);
            do_one(lc.c_raw_pos,   lc.c_out_pos);
            do_one(lc.c_raw_neg,   lc.c_out_neg);
        }

        out_lines.push_back(join_csv_line(fields));
    }

    fs::path csv_path(csv_main);
    fs::path tmp_path = csv_path;

    tmp_path += ".tmp";

    std::ofstream ofs(tmp_path);

    if (!ofs) {
        std::cerr << "[norm_cross_sections] ERROR: cannot open " << tmp_path
                  << " for writing.\n";
        return false;
    }

    for (const auto &lout : out_lines) {
        ofs << lout << "\n";
    }

    ofs.close();

    std::error_code ec;

    fs::rename(tmp_path, csv_path, ec);

    if (ec) {
        std::cerr << "[norm_cross_sections] ERROR: rename failed: " << ec.message() << "\n";
        return false;
    }

    std::cout << "[norm_cross_sections] Updated CSV with normed cross sections: "
              << csv_main << "\n";

    return true;
}

// -----------------------------------------------------------------------------
// Plotting normalized cross sections (with theory curves)
// -----------------------------------------------------------------------------

bool plot_normed_cross_sections_for_label(const std::string &csv_main,
                                          const std::string &label,
                                          const std::string &theory_json_root,
                                          const std::string &out_root_dir) {
    std::ifstream ifs(csv_main);

    if (!ifs) {
        std::cerr << "[norm_cross_sections] ERROR: cannot open " << csv_main
                  << " for plotting.\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }

    ifs.close();

    if (lines.empty()) {
        std::cerr << "[norm_cross_sections] ERROR: CSV " << csv_main
                  << " is empty in plot_normed_cross_sections_for_label.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    int c_xb_min = -1;
    int c_xb_max = -1;
    int c_q2_min = -1;
    int c_q2_max = -1;
    int c_t_min  = -1;
    int c_t_max  = -1;

    try {
        c_xb_min = find_col_required(header, "xBmin");
        c_xb_max = find_col_required(header, "xBmax");
        c_q2_min = find_col_required(header, "Q2min");
        c_q2_max = find_col_required(header, "Q2max");
        c_t_min  = find_col_required(header, "t_abs_min");
        c_t_max  = find_col_required(header, "t_abs_max");
    } catch (const std::exception &e) {
        std::cerr << "[norm_cross_sections] FATAL: " << e.what()
                  << " in plot_normed_cross_sections_for_label.\n";
        return false;
    }

    const int c_xb_idx = find_col_optional(header, "xB index");

    int c_phiavg = find_col_optional(header, "phiavg, " + label);
    int c_phimin = -1;
    int c_phimax = -1;

    if (c_phiavg < 0) {
        c_phimin = find_col_optional(header, "phimin");
        c_phimax = find_col_optional(header, "phimax");

        if (c_phimin < 0 || c_phimax < 0) {
            std::cerr << "[norm_cross_sections] FATAL: no phiavg or phimin/phimax "
                      << "available for label " << label << ".\n";
            return false;
        }
    }

    const int c_xs_unpol = find_col_optional(
        header, "normed cross sections, ep->epg, exp, " + label + ", unpol");
    const int c_xs_pos   = find_col_optional(
        header, "normed cross sections, ep->epg, exp, " + label + ", pos");
    const int c_xs_neg   = find_col_optional(
        header, "normed cross sections, ep->epg, exp, " + label + ", neg");

    if (c_xs_unpol < 0 && c_xs_pos < 0 && c_xs_neg < 0) {
        std::cerr << "[norm_cross_sections] INFO: no normed cross section columns for label "
                  << label << "; nothing to plot.\n";
        return true;
    }

    if (c_xs_unpol < 0 || c_xs_pos < 0 || c_xs_neg < 0) {
        std::cerr << "[norm_cross_sections] FATAL: incomplete normed cross section columns "
                  << "for label " << label << " (need unpol/pos/neg).\n";
        return false;
    }

    std::map<Range, XSGroupByXB> by_xb;

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) continue;

        std::vector<std::string> fields = split_csv_line(lines[row]);

        if (fields.size() != header.size()) continue;

        const double xbmin = std::atof(trim(fields[c_xb_min]).c_str());
        const double xbmax = std::atof(trim(fields[c_xb_max]).c_str());
        const double q2min = std::atof(trim(fields[c_q2_min]).c_str());
        const double q2max = std::atof(trim(fields[c_q2_max]).c_str());
        const double tmin  = std::atof(trim(fields[c_t_min]).c_str());
        const double tmax  = std::atof(trim(fields[c_t_max]).c_str());

        Range xb_range(xbmin, xbmax);
        Range q2_range(q2min, q2max);
        Range t_range(tmin,  tmax);

        double phi = 0.0;

        if (c_phiavg >= 0) {
            phi = std::atof(trim(fields[c_phiavg]).c_str());
        } else {
            const double pmin = std::atof(trim(fields[c_phimin]).c_str());
            const double pmax = std::atof(trim(fields[c_phimax]).c_str());

            phi = 0.5 * (pmin + pmax);
        }

        const Triple xs_unpol = parse_tuple3(fields[c_xs_unpol]);
        const Triple xs_pos   = parse_tuple3(fields[c_xs_pos]);
        const Triple xs_neg   = parse_tuple3(fields[c_xs_neg]);

        if (xs_unpol.value <= 0.0 && xs_pos.value <= 0.0 && xs_neg.value <= 0.0) {
            continue;
        }

        XSGroupByXB &group = by_xb[xb_range];

        if (group.xb_index < 0 && c_xb_idx >= 0) {
            group.xb_index = std::atoi(trim(fields[c_xb_idx]).c_str());
        }

        QTKey key(q2_range, t_range);
        BinData &bin = group.bins[key];

        if (!bin.have_theory_row) {
            bin.theory_row      = row;
            bin.have_theory_row = true;
        }

        auto add_point = [&](const Triple &xs, std::vector<Point> &vec) {
            if (!(xs.value > 0.0) || !std::isfinite(xs.value)) return;

            Point p;
            p.phi    = phi;
            p.xs     = xs.value;
            p.xs_err = xs.stat;

            vec.push_back(p);
        };

        add_point(xs_unpol, bin.unpol);
        add_point(xs_pos,   bin.pos);
        add_point(xs_neg,   bin.neg);
    }

    if (by_xb.empty()) {
        std::cerr << "[norm_cross_sections] WARNING: no normed xsec data found for label "
                  << label << " to plot.\n";
        return true;
    }

    const std::map<size_t, TheoryCurves> theory =
        load_theory_for_label(label, theory_json_root);

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    fs::path outdir = fs::path(out_root_dir) / canonical_period_dir(label);

    ensure_dir(outdir);

    int xb_canvas_counter = 0;

    for (const auto &kv_xb : by_xb) {
        const Range &xb_range = kv_xb.first;
        const XSGroupByXB &group = kv_xb.second;
        const auto &bins_for_xB = group.bins;

        if (bins_for_xB.empty()) continue;

        std::set<Range> q2_set;
        std::set<Range> t_set;

        for (const auto &kv : bins_for_xB) {
            const QTKey &qt = kv.first;

            q2_set.insert(qt.first);
            t_set.insert(qt.second);
        }

        std::vector<Range> q2_slice(q2_set.begin(), q2_set.end());
        std::vector<Range> t_slice(t_set.begin(), t_set.end());

        if (q2_slice.empty() || t_slice.empty()) continue;

        const int ncols = (int)q2_slice.size();
        const int nrows = (int)t_slice.size();
        const int nPads = ncols * nrows;

        const int xb_idx_for_name =
            (group.xb_index >= 0 ? group.xb_index : xb_canvas_counter);

        make_normed_xsec_canvas_for_mode(label, xb_range, group,
                                         q2_slice, t_slice,
                                         theory, outdir,
                                         xb_idx_for_name,
                                         XSecPanelMode::All,
                                         ncols, nrows, nPads);

        make_normed_xsec_canvas_for_mode(label, xb_range, group,
                                         q2_slice, t_slice,
                                         theory, outdir,
                                         xb_idx_for_name,
                                         XSecPanelMode::UnpolOnly,
                                         ncols, nrows, nPads);

        make_normed_xsec_canvas_for_mode(label, xb_range, group,
                                         q2_slice, t_slice,
                                         theory, outdir,
                                         xb_idx_for_name,
                                         XSecPanelMode::PosOnly,
                                         ncols, nrows, nPads);

        make_normed_xsec_canvas_for_mode(label, xb_range, group,
                                         q2_slice, t_slice,
                                         theory, outdir,
                                         xb_idx_for_name,
                                         XSecPanelMode::NegOnly,
                                         ncols, nrows, nPads);

        make_ratio_canvas_for_mode(label, xb_range, group,
                                   q2_slice, t_slice,
                                   theory, outdir,
                                   xb_idx_for_name,
                                   XSecPanelMode::UnpolOnly,
                                   ncols, nrows, nPads);

        make_ratio_canvas_for_mode(label, xb_range, group,
                                   q2_slice, t_slice,
                                   theory, outdir,
                                   xb_idx_for_name,
                                   XSecPanelMode::PosOnly,
                                   ncols, nrows, nPads);

        make_ratio_canvas_for_mode(label, xb_range, group,
                                   q2_slice, t_slice,
                                   theory, outdir,
                                   xb_idx_for_name,
                                   XSecPanelMode::NegOnly,
                                   ncols, nrows, nPads);

        ++xb_canvas_counter;
    }

    std::cout << "[norm_cross_sections] Plotted normed cross sections for label "
              << label << " into " << outdir.string() << "\n";

    return true;
}