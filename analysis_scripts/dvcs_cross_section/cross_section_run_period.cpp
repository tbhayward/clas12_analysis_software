#include "cross_section_run_period_consistency.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TMath.h>
#include <TROOT.h>
#include <TGraph.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

static inline void info_rpc(const std::string& s) {
    std::cout << "[run-period-consistency] " << s << std::endl;
}

static inline void fatal_rpc(const std::string& s) {
    std::cerr << "[run-period-consistency][FATAL] " << s << std::endl;
    throw std::runtime_error(s);
}

static inline std::string trim_rpc(const std::string& s) {
    size_t i0 = 0;
    while (i0 < s.size() && std::isspace((unsigned char)s[i0])) {
        ++i0;
    }
    size_t i1 = s.size();
    while (i1 > i0 && std::isspace((unsigned char)s[i1 - 1])) {
        --i1;
    }
    return s.substr(i0, i1 - i0);
}

static std::vector<std::string> split_csv_line_rpc(const std::string& line) {
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
    // endfor
    out.push_back(cur);
    return out;
}

static std::unordered_map<std::string, int>
build_header_index_rpc(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> out;
    for (int i = 0; i < (int)header.size(); ++i) {
        out[header[i]] = i;
    }
    // endfor
    return out;
}

static const std::string& get_col_ref_rpc(const std::vector<std::string>& row,
                                          const std::unordered_map<std::string, int>& idx,
                                          const std::string& name) {
    auto it = idx.find(name);
    if (it == idx.end()) {
        static const std::string empty;
        return empty;
    }
    const int j = it->second;
    if (j < 0 || j >= (int)row.size()) {
        static const std::string empty;
        return empty;
    }
    return row[j];
}

static double ToDouble_rpc(const std::string& s) {
    if (s.empty()) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static int ToInt_rpc(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static void require_columns_rpc(const std::unordered_map<std::string, int>& idx,
                                const std::vector<std::string>& cols,
                                const std::string& which_csv) {
    std::vector<std::string> missing;
    for (const auto& c : cols) {
        if (idx.find(c) == idx.end()) {
            missing.push_back(c);
        }
        // endfor
    }
    // endfor

    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "Missing columns in " << which_csv << ": ";
        for (size_t i = 0; i < missing.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << '"' << missing[i] << '"';
        }
        fatal_rpc(oss.str());
    }
}

struct TripletValue_rpc {
    double value    = 0.0;
    double stat_err = 0.0;
    double syst_err = 0.0;
    bool   valid    = false;
};

static TripletValue_rpc parse_triplet_rpc(const std::string& s) {
    TripletValue_rpc tv;
    if (s.empty()) return tv;

    const size_t p0 = s.find('(');
    const size_t p1 = s.find(')');
    if (p0 == std::string::npos || p1 == std::string::npos || p1 <= p0) {
        return tv;
    }
    // endif

    const std::string inner = s.substr(p0 + 1, p1 - p0 - 1);
    std::vector<std::string> parts;
    std::string cur;

    for (char c : inner) {
        if (c == ',') {
            parts.push_back(trim_rpc(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    // endfor
    parts.push_back(trim_rpc(cur));

    if (parts.size() >= 2) {
        tv.value    = ToDouble_rpc(parts[0]);
        tv.stat_err = ToDouble_rpc(parts[1]);
        tv.syst_err = (parts.size() >= 3) ? ToDouble_rpc(parts[2]) : 0.0;
        tv.valid    = true;
    }
    // endif

    return tv;
}

struct PeriodDef_rpc {
    std::string label;
    std::string short_tag;
    int color;
    int marker;
};

static const std::vector<PeriodDef_rpc>& all_period_defs_rpc() {
    static const std::vector<PeriodDef_rpc> defs = {
        {"Fa18 Inb", "fa18_inb", kBlack,      20},
        {"Fa18 Out", "fa18_out", kRed + 1,    24},
        {"Sp19 Inb", "sp19_inb", kBlue + 1,   21},
        {"Sp18 Inb", "sp18_inb", kGreen + 2,  25},
        {"Sp18 Out", "sp18_out", kMagenta+2,  33}
    };
    return defs;
}

static std::vector<PeriodDef_rpc> regular_period_defs_rpc(const std::string& helicity) {
    if (helicity == "unpol") {
        return {
            {"Fa18 Inb", "fa18_inb", kBlack,      20},
            {"Fa18 Out", "fa18_out", kRed + 1,    24},
            {"Sp18 Inb", "sp18_inb", kGreen + 2,  25},
            {"Sp18 Out", "sp18_out", kMagenta+2,  33}
        };
    }
    // endif

    if (helicity == "pos" || helicity == "neg") {
        return {
            {"Fa18 Inb", "fa18_inb", kBlack,   20},
            {"Fa18 Out", "fa18_out", kRed + 1, 24}
        };
    }
    // endif

    fatal_rpc("Unknown helicity in regular_period_defs_rpc: " + helicity);
    return {};
}

static std::vector<PeriodDef_rpc> pairable_period_defs_rpc(const std::string& helicity) {
    if (helicity == "unpol") {
        return all_period_defs_rpc();
    }
    // endif

    if (helicity == "pos" || helicity == "neg") {
        return {
            {"Fa18 Inb", "fa18_inb", kBlack,    20},
            {"Fa18 Out", "fa18_out", kRed + 1,  24},
            {"Sp19 Inb", "sp19_inb", kBlue + 1, 21}
        };
    }
    // endif

    fatal_rpc("Unknown helicity in pairable_period_defs_rpc: " + helicity);
    return {};
}

struct PairMode_rpc {
    std::string dir_name;
    std::string label;
    std::vector<PeriodDef_rpc> periods;
};

static std::vector<PairMode_rpc> build_pair_modes_rpc(const std::string& helicity) {
    const std::vector<PeriodDef_rpc> defs = pairable_period_defs_rpc(helicity);
    std::vector<PairMode_rpc> out;

    for (size_t i = 0; i < defs.size(); ++i) {
        for (size_t j = i + 1; j < defs.size(); ++j) {
            PairMode_rpc pm;
            pm.dir_name = defs[i].short_tag + "_vs_" + defs[j].short_tag;
            pm.label = defs[i].label + " vs " + defs[j].label;
            pm.periods = {defs[i], defs[j]};
            out.push_back(pm);
        }
        // endfor
    }
    // endfor

    return out;
}

struct Row_rpc {
    int    bin_index = 0;
    double xBmin     = 0.0;
    double xBmax     = 0.0;
    double Q2min     = 0.0;
    double Q2max     = 0.0;
    double tmin      = 0.0;
    double tmax      = 0.0;
    double xBavg     = 0.0;
    double Q2avg     = 0.0;
    double tavg      = 0.0;
    double phiavg    = 0.0;
    std::unordered_map<std::string, TripletValue_rpc> xs_by_col;
};

struct AxisSets_rpc {
    std::vector<std::pair<double, double>> xB;
    std::map<int, std::vector<std::pair<double, double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double, double>>> t_by_ix;
};

static AxisSets_rpc build_axes_rpc(const std::vector<Row_rpc>& rows) {
    std::set<std::pair<double, double>> xbset;
    std::map<std::pair<double, double>, std::set<std::pair<double, double>>> q2set_by_xb;
    std::map<std::pair<double, double>, std::set<std::pair<double, double>>> tset_by_xb;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert({r.Q2min, r.Q2max});
        tset_by_xb[xb].insert({r.tmin, r.tmax});
    }
    // endfor

    AxisSets_rpc ax;
    ax.xB.assign(xbset.begin(), xbset.end());
    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& xb = ax.xB[ix];
        ax.Q2_by_ix[ix] = {q2set_by_xb[xb].begin(), q2set_by_xb[xb].end()};
        ax.t_by_ix[ix]  = {tset_by_xb[xb].begin(),  tset_by_xb[xb].end()};
    }
    // endfor
    return ax;
}

static int find_index_rpc(const std::pair<double, double>& r,
                          const std::vector<std::pair<double, double>>& v) {
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] == r) return i;
    }
    // endfor
    return -1;
}

struct PanelPoint_rpc {
    double phi   = 0.0;
    double ratio = 0.0;
    double err   = 0.0;
};

struct PanelSeries_rpc {
    std::vector<double> phi;
    std::vector<double> val;
    std::vector<double> err;
};

struct SpreadBandPoint_rpc {
    double phi               = 0.0;
    double lo                = 0.0;
    double hi                = 0.0;
    double spread            = 0.0;
    double spread_err        = 0.0;
    double weight            = 0.0;
    double weighted_mean_raw = 0.0;
};

struct PanelChi2_rpc {
    double chi2                = 0.0;
    int    ndof                = 0;
    int    phi_rows_used       = 0;
    int    points_used         = 0;
    int    rows_total          = 0;

    double weighted_spread     = 0.0;
    double sum_w_spread        = 0.0;
    double sum_wr_spread       = 0.0;

    double reduced_chi2() const {
        return (ndof > 0) ? (chi2 / (double)ndof) : 0.0;
    }
};

struct PanelBundle_rpc {
    std::map<std::string, PanelSeries_rpc> series_by_period;
    std::map<std::string, std::vector<PanelPoint_rpc>> ratio_series_by_period;
    std::vector<SpreadBandPoint_rpc> spread_band;
    double Q2avg = 0.0;
    double tavg  = 0.0;
    PanelChi2_rpc chi2;
};

struct GlobalSpreadAccumulator_rpc {
    double sum_w_spread  = 0.0;
    double sum_wr_spread = 0.0;
    int    phi_rows_used = 0;
    int    panels_used   = 0;

    double mean() const {
        return (sum_w_spread > 0.0) ? (sum_wr_spread / sum_w_spread) : 0.0;
    }
};

static std::string xs_col_rpc(const std::string& period_label,
                              const std::string& helicity) {
    return "normed cross sections, ep->epg, exp, " + period_label + ", " + helicity;
}

static bool has_helicity_resolved_normed_cross_sections_rpc(const std::string& period_label) {
    if (period_label == "Sp18 Inb" || period_label == "Sp18 Out") {
        return false;
    }
    // endif

    return true;
}

static std::vector<PeriodDef_rpc> period_defs_for_required_columns_rpc(const std::string& helicity) {
    if (helicity == "unpol") {
        return all_period_defs_rpc();
    }
    // endif

    if (helicity == "pos" || helicity == "neg") {
        std::vector<PeriodDef_rpc> out;
        for (const auto& p : all_period_defs_rpc()) {
            if (!has_helicity_resolved_normed_cross_sections_rpc(p.label)) {
                continue;
            }
            // endif
            out.push_back(p);
        }
        // endfor
        return out;
    }
    // endif

    fatal_rpc("Unknown helicity in period_defs_for_required_columns_rpc: " + helicity);
    return {};
}

static std::vector<std::string> all_required_xs_columns_rpc() {
    std::vector<std::string> cols;
    const std::vector<std::string> helicities = {"unpol", "pos", "neg"};

    for (const auto& h : helicities) {
        const std::vector<PeriodDef_rpc> defs = period_defs_for_required_columns_rpc(h);
        for (const auto& p : defs) {
            cols.push_back(xs_col_rpc(p.label, h));
        }
        // endfor
    }
    // endfor

    return cols;
}

static std::vector<Row_rpc> load_rows_rpc(const std::string& hayward_csv_path) {
    std::ifstream fin(hayward_csv_path);
    if (!fin.is_open()) {
        fatal_rpc("Cannot open CSV: " + hayward_csv_path);
    }
    // endif

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal_rpc("CSV appears empty: " + hayward_csv_path);
    }
    // endif

    const std::vector<std::string> header = split_csv_line_rpc(header_line);
    const auto H = build_header_index_rpc(header);

    std::vector<std::string> required = {
        "bin index",
        "valid bin",
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "xBavg, 10.6 GeV",
        "Q2avg, 10.6 GeV",
        "t_abs_avg, 10.6 GeV",
        "phiavg, 10.6 GeV"
    };

    const auto xs_cols = all_required_xs_columns_rpc();
    required.insert(required.end(), xs_cols.begin(), xs_cols.end());
    require_columns_rpc(H, required, "Hayward CSV");

    std::vector<Row_rpc> rows;
    std::string line;
    int input_rows = 0;
    int kept_rows  = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        const std::vector<std::string> cols = split_csv_line_rpc(line);
        if (ToInt_rpc(get_col_ref_rpc(cols, H, "valid bin")) != 1) {
            continue;
        }
        // endif

        Row_rpc r;
        r.bin_index = ToInt_rpc(get_col_ref_rpc(cols, H, "bin index"));
        r.xBmin     = ToDouble_rpc(get_col_ref_rpc(cols, H, "xBmin"));
        r.xBmax     = ToDouble_rpc(get_col_ref_rpc(cols, H, "xBmax"));
        r.Q2min     = ToDouble_rpc(get_col_ref_rpc(cols, H, "Q2min"));
        r.Q2max     = ToDouble_rpc(get_col_ref_rpc(cols, H, "Q2max"));
        r.tmin      = ToDouble_rpc(get_col_ref_rpc(cols, H, "t_abs_min"));
        r.tmax      = ToDouble_rpc(get_col_ref_rpc(cols, H, "t_abs_max"));
        r.xBavg     = ToDouble_rpc(get_col_ref_rpc(cols, H, "xBavg, 10.6 GeV"));
        r.Q2avg     = ToDouble_rpc(get_col_ref_rpc(cols, H, "Q2avg, 10.6 GeV"));
        r.tavg      = ToDouble_rpc(get_col_ref_rpc(cols, H, "t_abs_avg, 10.6 GeV"));
        r.phiavg    = ToDouble_rpc(get_col_ref_rpc(cols, H, "phiavg, 10.6 GeV"));

        bool has_any = false;
        for (const auto& cname : xs_cols) {
            TripletValue_rpc tv = parse_triplet_rpc(get_col_ref_rpc(cols, H, cname));
            if (tv.valid) has_any = true;
            r.xs_by_col[cname] = tv;
        }
        // endfor

        if (!has_any) continue;
        rows.push_back(r);
        ++kept_rows;
    }
    // endwhile

    info_rpc("CSV rows read: " + std::to_string(input_rows));
    info_rpc("Valid rows kept: " + std::to_string(kept_rows));
    return rows;
}

static TGraphErrors* graph_pe1_rpc(const std::vector<double>& X,
                                   const std::vector<double>& Y,
                                   const std::vector<double>& EY,
                                   int markerStyle,
                                   int color) {
    if (X.empty()) return nullptr;

    std::vector<double> ex(X.size(), 0.0);
    auto* g = new TGraphErrors((int)X.size(),
                               const_cast<double*>(X.data()),
                               const_cast<double*>(Y.data()),
                               ex.data(),
                               const_cast<double*>(EY.data()));
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerSize(0.95);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);
    g->Draw("PE1 SAME");
    return g;
}

struct WeightedMeanInfo_rpc {
    double mean    = 0.0;
    double err     = 0.0;
    double sum_w   = 0.0;
    bool   valid   = false;
};

static WeightedMeanInfo_rpc weighted_mean_rpc(const std::vector<double>& values,
                                              const std::vector<double>& errs) {
    WeightedMeanInfo_rpc out;
    if (values.size() != errs.size() || values.empty()) return out;

    double sum_w  = 0.0;
    double sum_wx = 0.0;

    for (size_t i = 0; i < values.size(); ++i) {
        const double sigma = errs[i];
        if (sigma <= 0.0) continue;
        const double w = 1.0 / (sigma * sigma);
        sum_w  += w;
        sum_wx += w * values[i];
    }
    // endfor

    if (sum_w <= 0.0) return out;

    out.mean  = sum_wx / sum_w;
    out.err   = std::sqrt(1.0 / sum_w);
    out.sum_w = sum_w;
    out.valid = true;
    return out;
}

static double ratio_err_to_weighted_mean_rpc(double xi,
                                             double si,
                                             double weighted_mean,
                                             double wi_over_sumw) {
    if (weighted_mean == 0.0) return 0.0;
    if (si <= 0.0) return 0.0;

    const double dr_dxi = (weighted_mean - xi * wi_over_sumw) / (weighted_mean * weighted_mean);
    return std::fabs(dr_dxi) * si;
}

static PanelBundle_rpc make_panel_bundle_rpc(const std::vector<Row_rpc>& rows,
                                             int ix,
                                             int iQ,
                                             int it,
                                             const AxisSets_rpc& ax,
                                             const std::string& helicity,
                                             const std::vector<PeriodDef_rpc>& periods) {
    PanelBundle_rpc out;

    for (const auto& p : periods) {
        out.series_by_period[p.label] = PanelSeries_rpc();
        out.ratio_series_by_period[p.label] = std::vector<PanelPoint_rpc>();
    }
    // endfor

    for (const auto& r : rows) {
        const int ix_row = find_index_rpc({r.xBmin, r.xBmax}, ax.xB);
        if (ix_row != ix) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);
        const int iQ_row = find_index_rpc({r.Q2min, r.Q2max}, Q2s);
        const int it_row = find_index_rpc({r.tmin, r.tmax}, Ts);
        if (iQ_row != iQ || it_row != it) continue;

        ++out.chi2.rows_total;

        struct OneVal_rpc {
            std::string label;
            double value;
            double err;
        };

        std::vector<OneVal_rpc> valid_points;

        for (const auto& p : periods) {
            const std::string cname = xs_col_rpc(p.label, helicity);
            auto itx = r.xs_by_col.find(cname);
            if (itx == r.xs_by_col.end()) continue;

            const TripletValue_rpc& tv = itx->second;
            if (!tv.valid || tv.value <= 0.0 || tv.stat_err <= 0.0) continue;

            PanelSeries_rpc& s = out.series_by_period[p.label];
            s.phi.push_back(r.phiavg);
            s.val.push_back(tv.value);
            s.err.push_back(tv.stat_err);

            valid_points.push_back({p.label, tv.value, tv.stat_err});
        }
        // endfor

        const int n = (int)valid_points.size();
        if (n < 2) continue;

        std::vector<double> values;
        std::vector<double> errs;
        values.reserve(valid_points.size());
        errs.reserve(valid_points.size());

        for (const auto& vp : valid_points) {
            values.push_back(vp.value);
            errs.push_back(vp.err);
        }
        // endfor

        const WeightedMeanInfo_rpc wm = weighted_mean_rpc(values, errs);
        if (!wm.valid || wm.mean <= 0.0 || wm.sum_w <= 0.0) continue;

        double q2avg_sum = out.Q2avg * (double)out.chi2.phi_rows_used;
        double tavg_sum  = out.tavg  * (double)out.chi2.phi_rows_used;
        out.Q2avg = (q2avg_sum + r.Q2avg) / (double)(out.chi2.phi_rows_used + 1);
        out.tavg  = (tavg_sum  + r.tavg ) / (double)(out.chi2.phi_rows_used + 1);

        double min_ratio =  1.0e300;
        double max_ratio = -1.0e300;
        double chi2_row  =  0.0;
        double sum_wr    =  0.0;
        double sum_wrr2  =  0.0;

        for (const auto& vp : valid_points) {
            const double sigma = vp.err;
            if (sigma <= 0.0) continue;

            const double wi = 1.0 / (sigma * sigma);
            const double wi_over_sumw = wi / wm.sum_w;

            const double pull = (vp.value - wm.mean) / sigma;
            chi2_row += pull * pull;

            const double ratio = vp.value / wm.mean;
            const double ratio_err = ratio_err_to_weighted_mean_rpc(vp.value,
                                                                    vp.err,
                                                                    wm.mean,
                                                                    wi_over_sumw);

            PanelPoint_rpc rp;
            rp.phi   = r.phiavg;
            rp.ratio = ratio;
            rp.err   = ratio_err;
            out.ratio_series_by_period[vp.label].push_back(rp);

            sum_wr   += wi * ratio;
            sum_wrr2 += wi * ratio * ratio;

            min_ratio = std::min(min_ratio, ratio);
            max_ratio = std::max(max_ratio, ratio);
        }
        // endfor

        if (min_ratio >= 1.0e299 || max_ratio <= -1.0e299) continue;

        const double ratio_mean = sum_wr / wm.sum_w;
        double variance = sum_wrr2 / wm.sum_w - ratio_mean * ratio_mean;
        if (variance < 0.0) variance = 0.0;
        const double spread = std::sqrt(variance);

        SpreadBandPoint_rpc bp;
        bp.phi               = r.phiavg;
        bp.lo                = min_ratio;
        bp.hi                = max_ratio;
        bp.spread            = spread;
        bp.spread_err        = wm.err / wm.mean;
        bp.weight            = wm.sum_w;
        bp.weighted_mean_raw = wm.mean;
        out.spread_band.push_back(bp);

        out.chi2.sum_w_spread  += bp.weight;
        out.chi2.sum_wr_spread += bp.weight * bp.spread;
        out.chi2.weighted_spread =
            out.chi2.sum_wr_spread / out.chi2.sum_w_spread;

        out.chi2.chi2 += chi2_row;
        out.chi2.ndof += (n - 1);
        ++out.chi2.phi_rows_used;
        out.chi2.points_used += n;
    }
    // endfor

    auto sort_one = [](PanelSeries_rpc& s) {
        if (s.phi.empty()) return;

        std::vector<size_t> idx(s.phi.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
            return s.phi[a] < s.phi[b];
        });

        std::vector<double> phi2, val2, err2;
        phi2.reserve(idx.size());
        val2.reserve(idx.size());
        err2.reserve(idx.size());
        for (size_t i : idx) {
            phi2.push_back(s.phi[i]);
            val2.push_back(s.val[i]);
            err2.push_back(s.err[i]);
        }
        // endfor

        s.phi = std::move(phi2);
        s.val = std::move(val2);
        s.err = std::move(err2);
    };

    auto sort_ratio_one = [](std::vector<PanelPoint_rpc>& pts) {
        std::sort(pts.begin(), pts.end(), [](const PanelPoint_rpc& a, const PanelPoint_rpc& b) {
            return a.phi < b.phi;
        });
    };

    for (auto& kv : out.series_by_period) {
        sort_one(kv.second);
    }
    // endfor

    for (auto& kv : out.ratio_series_by_period) {
        sort_ratio_one(kv.second);
    }
    // endfor

    std::sort(out.spread_band.begin(), out.spread_band.end(), [](const SpreadBandPoint_rpc& a,
                                                                 const SpreadBandPoint_rpc& b) {
        return a.phi < b.phi;
    });

    return out;
}

static bool panel_has_content_rpc(const PanelBundle_rpc& pb) {
    return (pb.chi2.phi_rows_used > 0 && pb.chi2.sum_w_spread > 0.0 && !pb.spread_band.empty());
}

static int count_nonempty_panels_rpc(const std::map<std::tuple<int, int>, PanelBundle_rpc>& panel_map) {
    int count = 0;
    for (const auto& kv : panel_map) {
        if (panel_has_content_rpc(kv.second)) {
            ++count;
        }
    }
    // endfor
    return count;
}

static double weighted_spread_for_xb_rpc(const std::map<std::tuple<int, int>, PanelBundle_rpc>& panel_map) {
    double sum_w  = 0.0;
    double sum_ws = 0.0;

    for (const auto& kv : panel_map) {
        const PanelBundle_rpc& pb = kv.second;
        if (!panel_has_content_rpc(pb)) continue;
        sum_w  += pb.chi2.sum_w_spread;
        sum_ws += pb.chi2.sum_wr_spread;
    }
    // endfor

    return (sum_w > 0.0) ? (sum_ws / sum_w) : 0.0;
}

static void accumulate_global_spread_from_panel_map_rpc(
    const std::map<std::tuple<int, int>, PanelBundle_rpc>& panel_map,
    GlobalSpreadAccumulator_rpc& acc) {

    for (const auto& kv : panel_map) {
        const PanelBundle_rpc& pb = kv.second;
        if (!panel_has_content_rpc(pb)) {
            continue;
        }
        // endif

        acc.sum_w_spread  += pb.chi2.sum_w_spread;
        acc.sum_wr_spread += pb.chi2.sum_wr_spread;
        acc.phi_rows_used += pb.chi2.phi_rows_used;
        acc.panels_used   += 1;
    }
    // endfor
}

static void draw_overlay_canvas_rpc(const std::vector<Row_rpc>& rows,
                                    const AxisSets_rpc& ax,
                                    int ix,
                                    const std::string& helicity,
                                    const std::vector<PeriodDef_rpc>& periods,
                                    const std::string& top_label,
                                    const std::string& out_png,
                                    std::vector<double>& reduced_chi2_values,
                                    std::ofstream& summary_out,
                                    GlobalSpreadAccumulator_rpc& global_acc) {
    const auto& Q2s = ax.Q2_by_ix.at(ix);
    const auto& Ts  = ax.t_by_ix.at(ix);
    const int ncols = (int)Q2s.size();
    const int nrows = (int)Ts.size();
    if (ncols == 0 || nrows == 0) return;

    const int W = 320 * ncols + 240;
    const int H = 260 * nrows + 260;

    std::map<std::tuple<int, int>, PanelBundle_rpc> panel_map;

    for (int it = 0; it < nrows; ++it) {
        for (int iQ = 0; iQ < ncols; ++iQ) {
            PanelBundle_rpc pb = make_panel_bundle_rpc(rows, ix, iQ, it, ax, helicity, periods);
            panel_map[{iQ, it}] = pb;

            const PanelChi2_rpc& chi = pb.chi2;
            summary_out << std::fixed << std::setprecision(6)
                        << "xBmin=" << ax.xB[ix].first << "  "
                        << "xBmax=" << ax.xB[ix].second << "  "
                        << "Q2min=" << Q2s[iQ].first << "  "
                        << "Q2max=" << Q2s[iQ].second << "  "
                        << "tmin=" << Ts[it].first << "  "
                        << "tmax=" << Ts[it].second << "  "
                        << "chi2=" << chi.chi2 << "  "
                        << "ndof=" << chi.ndof << "  "
                        << "chi2_ndf=" << chi.reduced_chi2() << "  "
                        << "weighted_spread=" << chi.weighted_spread << "  "
                        << "phi_rows_used=" << chi.phi_rows_used << "  "
                        << "points_used=" << chi.points_used << "  "
                        << "rows_total=" << chi.rows_total
                        << "\n";

            if (chi.ndof > 0) {
                reduced_chi2_values.push_back(chi.reduced_chi2());
            }
            // endif
        }
        // endfor
    }
    // endfor

    const int nonempty_panels = count_nonempty_panels_rpc(panel_map);
    if (nonempty_panels <= 0) {
        summary_out << "Skipped empty overlay canvas: " << out_png << "\n";
        return;
    }
    // endif

    accumulate_global_spread_from_panel_map_rpc(panel_map, global_acc);

    const double global_weighted_spread = weighted_spread_for_xb_rpc(panel_map);

    const std::string cname = fs::path(out_png).filename().string();
    TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), W, H);
    c->SetFillColor(kWhite);
    c->SetFrameFillColor(kWhite);
    c->cd();

    TPad* pTop = new TPad(Form("pTop_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                          Form("pTop_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                          0.0, 0.90, 1.0, 1.0);
    pTop->SetFillColor(kWhite);
    pTop->SetFrameFillColor(kWhite);
    pTop->SetFillStyle(1001);
    pTop->SetBorderSize(0);
    pTop->Draw();
    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.14);
    head.DrawLatex(0.50, 0.55,
        Form("%s (%s)   x_{B} #in [%.3g, %.3g]   Global weighted mean spread = %.3f",
             top_label.c_str(),
             helicity.c_str(),
             ax.xB[ix].first,
             ax.xB[ix].second,
             global_weighted_spread));

    TLegend* leg = new TLegend(0.03, 0.02, 0.97, 0.48);
    leg->SetNColumns((int)periods.size());
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.16);

    std::vector<TMarker*> legend_markers;
    for (const auto& p : periods) {
        TMarker* m = new TMarker(0.0, 0.0, p.marker);
        m->SetMarkerColor(p.color);
        legend_markers.push_back(m);
        leg->AddEntry(m, p.label.c_str(), "p");
    }
    // endfor
    leg->Draw();

    c->cd();
    TPad* pGrid = new TPad(Form("pGrid_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                           Form("pGrid_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                           0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillColor(kWhite);
    pGrid->SetFrameFillColor(kWhite);
    pGrid->SetFillStyle(1001);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.00, 0.00);

    for (int it = 0; it < nrows; ++it) {
        for (int iQ = 0; iQ < ncols; ++iQ) {
            pGrid->cd(it * ncols + iQ + 1);
            gPad->SetFillColor(kWhite);
            gPad->SetFrameFillColor(kWhite);
            gPad->SetTicks(1, 1);
            gPad->SetTopMargin(0.14);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);
            gPad->SetLogy(0);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 2.0);
            frame->SetLineColor(kBlack);
            frame->SetTitle("");
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Run period / weighted mean");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetXaxis()->SetLabelSize(0.050);
            frame->GetYaxis()->SetLabelSize(0.050);
            frame->GetXaxis()->SetTitleOffset(1.10);
            frame->GetYaxis()->SetTitleOffset(1.55);

            const PanelBundle_rpc& pb = panel_map.at({iQ, it});
            if (!panel_has_content_rpc(pb)) {
                continue;
            }
            // endif

            const double redchi2 = pb.chi2.reduced_chi2();

            frame->SetTitleFont(42, "t");
            frame->SetTitleSize(0.040, "t");
            frame->SetTitle(Form("Q^{2}=%.2f, -t=%.2f, #chi^{2}/ndf=%s, <spread>_{w}=%.2f",
                                 pb.Q2avg,
                                 pb.tavg,
                                 (pb.chi2.ndof > 0 ? Form("%.2f", redchi2) : "n/a"),
                                 pb.chi2.weighted_spread));

            if (!pb.spread_band.empty()) {
                std::vector<double> xb, yb;
                xb.reserve(2 * pb.spread_band.size());
                yb.reserve(2 * pb.spread_band.size());

                for (const auto& bp : pb.spread_band) {
                    xb.push_back(bp.phi);
                    yb.push_back(bp.hi);
                }
                // endfor

                for (int ib = (int)pb.spread_band.size() - 1; ib >= 0; --ib) {
                    xb.push_back(pb.spread_band[ib].phi);
                    yb.push_back(pb.spread_band[ib].lo);
                }
                // endfor

                TGraph* gband = new TGraph((int)xb.size(), xb.data(), yb.data());
                gband->SetFillColorAlpha(kGray + 1, 0.25);
                gband->SetLineColor(kGray + 1);
                gband->SetLineWidth(1);
                gband->Draw("F SAME");
            }
            // endif

            TLine* unity = new TLine(0.0, 1.0, 360.0, 1.0);
            unity->SetLineColor(kGray + 2);
            unity->SetLineStyle(2);
            unity->SetLineWidth(1);
            unity->Draw("SAME");

            for (const auto& p : periods) {
                const auto itp = pb.ratio_series_by_period.find(p.label);
                if (itp == pb.ratio_series_by_period.end()) continue;

                const auto& pts = itp->second;
                if (pts.empty()) continue;

                std::vector<double> x, y, ey;
                x.reserve(pts.size());
                y.reserve(pts.size());
                ey.reserve(pts.size());

                for (const auto& pt : pts) {
                    x.push_back(pt.phi);
                    y.push_back(pt.ratio);
                    ey.push_back(pt.err);
                }
                // endfor

                graph_pe1_rpc(x, y, ey, p.marker, p.color);
            }
            // endfor
        }
        // endfor
    }
    // endfor

    c->SaveAs(out_png.c_str());

    delete leg;
    for (auto* m : legend_markers) delete m;
    delete c;
}

static void draw_half_range_canvas_rpc(const std::vector<Row_rpc>& rows,
                                       const AxisSets_rpc& ax,
                                       int ix,
                                       const std::string& helicity,
                                       const std::vector<PeriodDef_rpc>& periods,
                                       const std::string& top_label,
                                       const std::string& out_png,
                                       std::ofstream& summary_out) {
    const auto& Q2s = ax.Q2_by_ix.at(ix);
    const auto& Ts  = ax.t_by_ix.at(ix);
    const int ncols = (int)Q2s.size();
    const int nrows = (int)Ts.size();
    if (ncols == 0 || nrows == 0) return;

    const int W = 320 * ncols + 240;
    const int H = 260 * nrows + 260;

    std::map<std::tuple<int, int>, PanelBundle_rpc> panel_map;
    double ymax = 0.0;

    for (int it = 0; it < nrows; ++it) {
        for (int iQ = 0; iQ < ncols; ++iQ) {
            PanelBundle_rpc pb = make_panel_bundle_rpc(rows, ix, iQ, it, ax, helicity, periods);
            panel_map[{iQ, it}] = pb;

            if (!panel_has_content_rpc(pb)) continue;

            for (const auto& bp : pb.spread_band) {
                ymax = std::max(ymax, bp.spread + bp.spread_err);
            }
            // endfor
        }
        // endfor
    }
    // endfor

    const int nonempty_panels = count_nonempty_panels_rpc(panel_map);
    if (nonempty_panels <= 0) {
        summary_out << "Skipped empty spread canvas: " << out_png << "\n";
        return;
    }
    // endif

    if (ymax <= 0.0) ymax = 0.25;
    ymax *= 1.20;

    const double global_weighted_spread = weighted_spread_for_xb_rpc(panel_map);

    const std::string cname = fs::path(out_png).filename().string();
    TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), W, H);
    c->SetFillColor(kWhite);
    c->SetFrameFillColor(kWhite);
    c->cd();

    TPad* pTop = new TPad(Form("pTop_hr_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                          Form("pTop_hr_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                          0.0, 0.90, 1.0, 1.0);
    pTop->SetFillColor(kWhite);
    pTop->SetFrameFillColor(kWhite);
    pTop->SetFillStyle(1001);
    pTop->SetBorderSize(0);
    pTop->Draw();
    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.14);
    head.DrawLatex(0.50, 0.55,
        Form("%s (%s)   x_{B} #in [%.3g, %.3g]   Global weighted mean spread = %.3f",
             top_label.c_str(),
             helicity.c_str(),
             ax.xB[ix].first,
             ax.xB[ix].second,
             global_weighted_spread));

    c->cd();
    TPad* pGrid = new TPad(Form("pGrid_hr_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                           Form("pGrid_hr_rpc_%d_%s_%zu", ix, helicity.c_str(), std::hash<std::string>{}(top_label)),
                           0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillColor(kWhite);
    pGrid->SetFrameFillColor(kWhite);
    pGrid->SetFillStyle(1001);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.00, 0.00);

    for (int it = 0; it < nrows; ++it) {
        for (int iQ = 0; iQ < ncols; ++iQ) {
            pGrid->cd(it * ncols + iQ + 1);
            gPad->SetFillColor(kWhite);
            gPad->SetFrameFillColor(kWhite);
            gPad->SetTicks(1, 1);
            gPad->SetTopMargin(0.14);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);
            gPad->SetLogy(0);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, ymax);
            frame->SetLineColor(kBlack);
            frame->SetTitle("");
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Weighted spread");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetXaxis()->SetLabelSize(0.050);
            frame->GetYaxis()->SetLabelSize(0.050);
            frame->GetXaxis()->SetTitleOffset(1.10);
            frame->GetYaxis()->SetTitleOffset(1.55);

            const PanelBundle_rpc& pb = panel_map.at({iQ, it});
            if (!panel_has_content_rpc(pb)) {
                continue;
            }
            // endif

            frame->SetTitleFont(42, "t");
            frame->SetTitleSize(0.040, "t");
            frame->SetTitle(Form("Q^{2}=%.2f, -t=%.2f, <spread>_{w}=%.2f",
                                 pb.Q2avg,
                                 pb.tavg,
                                 pb.chi2.weighted_spread));

            if (!pb.spread_band.empty()) {
                std::vector<double> x, y, ey;
                x.reserve(pb.spread_band.size());
                y.reserve(pb.spread_band.size());
                ey.reserve(pb.spread_band.size());

                for (const auto& bp : pb.spread_band) {
                    x.push_back(bp.phi);
                    y.push_back(bp.spread);
                    ey.push_back(bp.spread_err);
                }
                // endfor

                std::vector<double> ex(x.size(), 0.0);
                TGraphErrors* g = new TGraphErrors((int)x.size(),
                                                   x.data(),
                                                   y.data(),
                                                   ex.data(),
                                                   ey.data());
                g->SetMarkerStyle(20);
                g->SetMarkerSize(0.80);
                g->SetMarkerColor(kBlack);
                g->SetLineColor(kBlack);
                g->SetLineWidth(1);
                g->Draw("PE1 SAME");
            }
            // endif
        }
        // endfor
    }
    // endfor

    c->SaveAs(out_png.c_str());
    summary_out << "Saved spread canvas: " << out_png << "\n";

    delete c;
}

static Double_t reduced_chi2_pdf_rpc(Double_t* x, Double_t* par) {
    const double xx  = x[0];
    const double amp = par[0];
    const double nu  = par[1];

    if (xx <= 0.0 || nu <= 0.0) return 0.0;

    const double half_nu = 0.5 * nu;
    const double norm = std::pow(0.5 * nu, half_nu) / TMath::Gamma(half_nu);
    return amp * norm * std::pow(xx, half_nu - 1.0) * std::exp(-0.5 * nu * xx);
}

static void draw_reduced_chi2_hist_rpc(const std::vector<double>& reduced_chi2_values,
                                       const std::string& helicity,
                                       const std::string& out_png,
                                       std::ofstream& summary_out) {
    if (reduced_chi2_values.empty()) {
        summary_out << "\nNo panels with positive ndf were found for helicity "
                    << helicity << ".\n";
        return;
    }
    // endif

    const double xmax = 50.0;
    const int nbins = 50;
    TH1D* h = new TH1D(Form("h_redchi2_%s_%lld",
                            helicity.c_str(),
                            (long long)std::hash<std::string>{}(out_png)),
                       "",
                       nbins,
                       0.0,
                       xmax);
    for (double v : reduced_chi2_values) {
        h->Fill(v);
    }
    // endfor

    TF1* f = new TF1(Form("f_redchi2_%s_%lld",
                          helicity.c_str(),
                          (long long)std::hash<std::string>{}(out_png)),
                     reduced_chi2_pdf_rpc,
                     0.0,
                     xmax,
                     2);
    f->SetParName(0, "A");
    f->SetParName(1, "#nu");
    f->SetParameters(h->GetMaximum(), 4.0);
    f->SetParLimits(0, 0.0, 1.0e9);
    f->SetParLimits(1, 0.1, 200.0);

    TCanvas* c = new TCanvas(Form("c_redchi2_%s_%lld",
                                  helicity.c_str(),
                                  (long long)std::hash<std::string>{}(out_png)),
                             "",
                             900,
                             700);
    c->SetFillColor(kWhite);
    c->SetFrameFillColor(kWhite);
    c->cd();
    gPad->SetFillColor(kWhite);
    gPad->SetFrameFillColor(kWhite);
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.06);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);

    h->SetLineColor(kBlack);
    h->SetLineWidth(2);
    h->GetXaxis()->SetTitle("Reduced #chi^{2}");
    h->GetYaxis()->SetTitle("Panels");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.050);
    h->GetYaxis()->SetTitleSize(0.050);
    h->GetXaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetLabelSize(0.045);
    h->GetYaxis()->SetTitleOffset(1.25);
    h->Draw("HIST");

    const int fit_status = h->Fit(f, "RQ0");
    f->SetLineColor(kRed + 1);
    f->SetLineWidth(2);
    f->Draw("SAME");

    const double mean = h->GetMean();
    const double rms  = h->GetRMS();
    const double nu   = f->GetParameter(1);
    const double nu_e = f->GetParError(1);

    TLatex lab;
    lab.SetNDC();
    lab.SetTextFont(42);
    lab.SetTextAlign(13);
    lab.SetTextSize(0.040);
    lab.DrawLatex(0.58, 0.90,
        Form("Reduced #chi^{2} distribution (%s)", helicity.c_str()));
    lab.DrawLatex(0.58, 0.84,
        Form("Entries = %d", (int)reduced_chi2_values.size()));
    lab.DrawLatex(0.58, 0.78,
        Form("Mean = %.3f", mean));
    lab.DrawLatex(0.58, 0.72,
        Form("RMS = %.3f", rms));
    lab.DrawLatex(0.58, 0.66,
        Form("Fit status = %d", fit_status));
    lab.DrawLatex(0.58, 0.60,
        Form("Effective #nu = %.3f #pm %.3f", nu, nu_e));

    c->SaveAs(out_png.c_str());

    summary_out << "\nReduced-chi2 histogram summary\n";
    summary_out << "entries=" << reduced_chi2_values.size() << "\n";
    summary_out << std::fixed << std::setprecision(6)
                << "hist_xmax=" << xmax << "\n"
                << "mean=" << mean << "\n"
                << "rms=" << rms << "\n"
                << "fit_status=" << fit_status << "\n"
                << "fit_nu=" << nu << "\n"
                << "fit_nu_err=" << nu_e << "\n"
                << "fit_model=reduced_chi2_pdf_effective_nu\n";

    delete f;
    delete h;
    delete c;
}

static void print_mode_global_summary_rpc(const std::string& mode_label,
                                          const std::string& helicity,
                                          const std::vector<PeriodDef_rpc>& periods,
                                          const GlobalSpreadAccumulator_rpc& acc,
                                          std::ofstream& summary_out) {
    std::ostringstream plist;
    for (size_t i = 0; i < periods.size(); ++i) {
        if (i > 0) plist << ", ";
        plist << periods[i].label;
    }
    // endfor

    std::ostringstream msg;
    msg << std::fixed << std::setprecision(6)
        << "GLOBAL SUMMARY | mode = " << mode_label
        << " | helicity = " << helicity
        << " | periods = [" << plist.str() << "]"
        << " | global_weighted_mean_spread = " << acc.mean()
        << " | total_phi_rows_used = " << acc.phi_rows_used
        << " | panels_used = " << acc.panels_used
        << " | sum_row_weights = " << acc.sum_w_spread;

    info_rpc(msg.str());

    summary_out << "\nGLOBAL SUMMARY\n";
    summary_out << "mode=" << mode_label << "\n";
    summary_out << "helicity=" << helicity << "\n";
    summary_out << "periods=[" << plist.str() << "]\n";
    summary_out << std::fixed << std::setprecision(6)
                << "global_weighted_mean_spread=" << acc.mean() << "\n"
                << "total_phi_rows_used=" << acc.phi_rows_used << "\n"
                << "panels_used=" << acc.panels_used << "\n"
                << "sum_row_weights=" << acc.sum_w_spread << "\n\n";
}

static void run_one_helicity_mode_rpc(const std::vector<Row_rpc>& rows,
                                      const AxisSets_rpc& ax,
                                      const std::string& helicity,
                                      const std::vector<PeriodDef_rpc>& periods,
                                      const std::string& mode_label,
                                      const std::string& output_base_dir) {
    const fs::path base_dir = fs::path(output_base_dir) / helicity;
    const fs::path overlay_dir = base_dir / "overlays";
    fs::create_directories(overlay_dir);

    std::ofstream summary_out((base_dir / "reduced_chi2_summary.txt").string());
    if (!summary_out.is_open()) {
        fatal_rpc("Could not open summary text file for writing in " + base_dir.string());
    }
    // endif

    summary_out << "Run-period consistency summary for mode = " << mode_label
                << " helicity = " << helicity << "\n";

    summary_out << "Included periods: ";
    for (size_t i = 0; i < periods.size(); ++i) {
        if (i > 0) summary_out << ", ";
        summary_out << periods[i].label;
    }
    // endfor
    summary_out << "\n";

    summary_out << "Definition: for each phi row, a weighted mean cross section is computed across the available included run periods using inverse-variance weights w_i = 1/sigma_i^2 from the statistical uncertainties of the corresponding normed cross-section tuples.\n";
    summary_out << "Overlay plots show the ratio of each included run period to that weighted mean. The plotted statistical uncertainty is propagated through r_i = x_i / mu_w, retaining the x_i contribution to mu_w.\n";
    summary_out << "For each phi row, the plotted spread is the inverse-variance-weighted RMS spread in ratio space: spread_phi = sqrt[ sum_i w_i (r_i - <r>_w)^2 / sum_i w_i ].\n";
    summary_out << "The reported panel value <spread>_w is the weighted mean of spread_phi across phi, with row weight W_phi = sum_i w_i.\n";
    summary_out << "Empty canvases are now skipped entirely, and empty subpads are never used in the global summary.\n\n";

    std::vector<double> reduced_chi2_values;
    GlobalSpreadAccumulator_rpc global_acc;

    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const fs::path out_png1 = overlay_dir / Form("period_consistency_xB_%d.png", ix);
        draw_overlay_canvas_rpc(rows,
                                ax,
                                ix,
                                helicity,
                                periods,
                                mode_label,
                                out_png1.string(),
                                reduced_chi2_values,
                                summary_out,
                                global_acc);
        info_rpc("Processed: " + out_png1.string());

        const fs::path out_png2 = overlay_dir / Form("period_consistency_half_range_xB_%d.png", ix);
        draw_half_range_canvas_rpc(rows,
                                   ax,
                                   ix,
                                   helicity,
                                   periods,
                                   mode_label,
                                   out_png2.string(),
                                   summary_out);
        info_rpc("Processed: " + out_png2.string());
    }
    // endfor

    print_mode_global_summary_rpc(mode_label, helicity, periods, global_acc, summary_out);

    const fs::path hist_png = base_dir / "reduced_chi2_distribution.png";
    draw_reduced_chi2_hist_rpc(reduced_chi2_values, helicity, hist_png.string(), summary_out);
    info_rpc("Saved: " + hist_png.string());
}

static void run_regular_mode_rpc(const std::vector<Row_rpc>& rows,
                                 const AxisSets_rpc& ax,
                                 const std::string& output_base_dir) {
    const fs::path regular_dir = fs::path(output_base_dir) / "regular";
    fs::create_directories(regular_dir);

    run_one_helicity_mode_rpc(rows,
                              ax,
                              "unpol",
                              regular_period_defs_rpc("unpol"),
                              "Regular run-period consistency",
                              regular_dir.string());

    run_one_helicity_mode_rpc(rows,
                              ax,
                              "pos",
                              regular_period_defs_rpc("pos"),
                              "Regular run-period consistency",
                              regular_dir.string());

    run_one_helicity_mode_rpc(rows,
                              ax,
                              "neg",
                              regular_period_defs_rpc("neg"),
                              "Regular run-period consistency",
                              regular_dir.string());
}

static void run_all_pair_modes_rpc(const std::vector<Row_rpc>& rows,
                                   const AxisSets_rpc& ax,
                                   const std::string& output_base_dir) {
    const fs::path pairs_dir = fs::path(output_base_dir) / "pairs";
    fs::create_directories(pairs_dir);

    const std::vector<std::string> helicities = {"unpol", "pos", "neg"};

    for (const auto& helicity : helicities) {
        const std::vector<PairMode_rpc> pair_modes = build_pair_modes_rpc(helicity);

        for (const auto& pm : pair_modes) {
            const fs::path mode_dir = pairs_dir / pm.dir_name;
            fs::create_directories(mode_dir);

            run_one_helicity_mode_rpc(rows,
                                      ax,
                                      helicity,
                                      pm.periods,
                                      pm.label,
                                      mode_dir.string());
        }
        // endfor
    }
    // endfor
}

} // anonymous namespace

void plot_cross_section_run_period_consistency(const std::string& hayward_csv_path,
                                               const std::string& output_base_dir) {
    fs::create_directories(output_base_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetCanvasColor(kWhite);
    gStyle->SetPadColor(kWhite);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetStatColor(kWhite);
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "");
    gStyle->SetTitleSize(0.050, "");

    const std::vector<Row_rpc> rows = load_rows_rpc(hayward_csv_path);
    const AxisSets_rpc ax = build_axes_rpc(rows);

    info_rpc("Number of xB bins found: " + std::to_string(ax.xB.size()));

    run_regular_mode_rpc(rows, ax, output_base_dir);
    run_all_pair_modes_rpc(rows, ax, output_base_dir);
}