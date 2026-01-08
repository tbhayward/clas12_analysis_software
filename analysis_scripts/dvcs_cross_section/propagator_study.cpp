// propagator_study.cpp
// -----------------------------------------------------------------------------
// Study of dvcsgen propagator (P1 / ycol-mirror) cut efficiency in DATA.
//
// For each (xB, Q2, |t|, phi) bin defined by dvcs_pass2_analysis.csv, we count:
//
//   N_noP1   : events passing global DVCS cuts with the P1 cut DISABLED
//   N_withP1 : events passing global DVCS cuts with the P1 cut ENABLED
//
// and plot the ratio vs phi:
//
//   eff(phi) = N_withP1 / N_noP1
//
// Notes / constraints:
// - Applies 3-sigma exclusivity cuts for DATA using combined_cuts.json, topology-by-topology.
// - Fails fast on missing CSV columns, missing period keys, missing required branches,
//   or malformed cut JSON.
// - This module only produces plots; it does not write to CSV.
// -----------------------------------------------------------------------------

#include "propagator_study.h"

#include "global_cuts.h"

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
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

// ROOT
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TH1.h>
#include <TLine.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace propagator_study {

using Range = std::pair<double, double>;

static constexpr double kPi = 3.14159265358979323846;

static void fatal(const std::string& msg) {
    throw std::runtime_error(std::string("[propagator_study] FATAL: ") + msg);
}

static std::string canonical_period_dir(const std::string& label) {
    if (label == "Fa18 Inb") return "Fa18_Inb";
    if (label == "Fa18 Out") return "Fa18_Out";
    if (label == "Sp18 Inb") return "Sp18_Inb";
    if (label == "Sp18 Out") return "Sp18_Out";
    if (label == "Sp19 Inb") return "Sp19_Inb";
    if (label == "Fa18")     return "Fa18";
    if (label == "Sp18")     return "Sp18";
    if (label == "10.6 GeV") return "10.6_GeV";

    fatal("Unknown label passed to canonical_period_dir: \"" + label + "\"");
    return "UNREACHABLE";
}

static void ensure_dir(const fs::path& p) {
    if (!fs::exists(p)) {
        fs::create_directories(p);
    }
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            in_quotes = !in_quotes;
            field.push_back(c);
        } else if (c == ',' && !in_quotes) {
            out.push_back(field);
            field.clear();
        } else {
            field.push_back(c);
        }
    }
    out.push_back(field);
    return out;
}

static std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) {
        ++b;
    }
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) {
        --e;
    }
    return s.substr(b, e - b);
}

static std::string unquote(const std::string& s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        const std::string inner = s.substr(1, s.size() - 2);
        std::string out;
        out.reserve(inner.size());
        for (size_t i = 0; i < inner.size(); ++i) {
            if (inner[i] == '"' && i + 1 < inner.size() && inner[i + 1] == '"') {
                out.push_back('"');
                ++i;
            } else {
                out.push_back(inner[i]);
            }
        }
        return out;
    }
    return s;
}

static int find_col(const std::vector<std::string>& header,
                    const std::string& target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        }
    }
    fatal("Missing required CSV column: \"" + target + "\"");
    return -1;
}

static int find_col_optional(const std::vector<std::string>& header,
                             const std::string& target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        }
    }
    return -1;
}

// Wrap angle into [0, 360)
static double wrap_phi_deg(double phi_deg) {
    double x = std::fmod(phi_deg, 360.0);
    if (x < 0.0) {
        x += 360.0;
    }
    return x;
}

static bool phi_in_range_deg(double phi_deg, double phimin, double phimax) {
    // Support wrap-around bins, e.g. [350, 10]
    phi_deg = wrap_phi_deg(phi_deg);
    phimin  = wrap_phi_deg(phimin);
    phimax  = wrap_phi_deg(phimax);

    if (phimin <= phimax) {
        return (phi_deg >= phimin && phi_deg < phimax);
    }
    return (phi_deg >= phimin || phi_deg < phimax);
}

static double phi_center_deg(double phimin, double phimax) {
    const double a = wrap_phi_deg(phimin);
    const double b = wrap_phi_deg(phimax);

    if (a <= b) {
        return 0.5 * (a + b);
    }
    const double bb = b + 360.0;
    const double c = 0.5 * (a + bb);
    return wrap_phi_deg(c);
}

static bool value_in_range(double v, const Range& r) {
    return (v >= r.first && v < r.second);
}

static std::string topology_from_detectors(int detector1, int detector2) {
    // Deterministic mapping (must match the rest of the analysis):
    // 1 = FD, 2 = CD, 3 = FT
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 3) return "CD_FT";
    return "UNKNOWN";
}

// -----------------------------------------------------------------------------
// 3-sigma cut map loader (combined_cuts.json) NEW SCHEMA
//
// Example key:
//   "DVCS_Fa18_Inb_CD_FD": {
//        "data": { "Emiss2": { "mean": ..., "std": ... }, ... },
//        "mc":   { ... }
//   }
//
// This module is DATA-only, so we load only the "data" branch.
// -----------------------------------------------------------------------------

struct Cut1D {
    double mean = 0.0;
    double std = 0.0;
};

using SigmaCut = std::map<std::string, Cut1D>;

// Index by (period_key, topo) -> SigmaCut
struct CutDB {
    std::map<std::string, std::map<std::string, SigmaCut>> cuts; // cuts[period_key][topo]
};

static bool read_mean_std_strict(const json& j, double& mean, double& stdv) {
    if (!j.is_object()) {
        return false;
    }
    if (!j.contains("mean") || !j.contains("std")) {
        return false;
    }
    if (!j["mean"].is_number() || !j["std"].is_number()) {
        return false;
    }
    mean = j["mean"].get<double>();
    stdv = j["std"].get<double>();
    return true;
}

static std::string dvcs_period_prefix_from_period_key(const std::string& period_key) {
    if (period_key == "fa18_inb") return "DVCS_Fa18_Inb";
    if (period_key == "fa18_out") return "DVCS_Fa18_Out";
    if (period_key == "sp18_inb") return "DVCS_Sp18_Inb";
    if (period_key == "sp18_out") return "DVCS_Sp18_Out";
    if (period_key == "sp19_inb") return "DVCS_Sp19_Inb";
    fatal("Unknown period_key in dvcs_period_prefix_from_period_key: \"" + period_key + "\"");
    return "UNREACHABLE";
}

static CutDB load_combined_cuts_json_newschema(const std::string& path,
                                               const std::vector<std::string>& required_period_keys) {
    std::ifstream ifs(path);
    if (!ifs) {
        fatal("Cannot open combined cuts JSON: " + path);
    }

    json j;
    try {
        ifs >> j;
    } catch (...) {
        fatal("Malformed JSON in combined cuts file: " + path);
    }

    if (!j.is_object()) {
        fatal("combined cuts JSON root is not an object: " + path);
    }

    const std::vector<std::string> required_topos = {"FD_FD", "CD_FD", "CD_FT"};

    CutDB db;

    // Strictly require that each (period_key, topo) DVCS key exists, and has a "data" object
    // with at least one variable containing numeric {mean, std}.
    for (const auto& pk : required_period_keys) {
        const std::string prefix = dvcs_period_prefix_from_period_key(pk);

        for (const auto& topo : required_topos) {
            const std::string full_key = prefix + "_" + topo;

            if (!j.contains(full_key)) {
                fatal("combined cuts JSON missing required key: \"" + full_key + "\"");
            }
            const json& entry = j[full_key];
            if (!entry.is_object()) {
                fatal("combined cuts JSON key \"" + full_key + "\" is not an object");
            }
            if (!entry.contains("data")) {
                fatal("combined cuts JSON key \"" + full_key + "\" missing required \"data\" object");
            }
            const json& data_obj = entry["data"];
            if (!data_obj.is_object()) {
                fatal("combined cuts JSON key \"" + full_key + "\" -> \"data\" is not an object");
            }

            SigmaCut sc;
            for (auto it = data_obj.begin(); it != data_obj.end(); ++it) {
                const std::string var = it.key();
                const json& spec = it.value();

                double mean = 0.0;
                double stdv = 0.0;
                if (!read_mean_std_strict(spec, mean, stdv)) {
                    fatal("Cut spec must contain numeric keys {mean, std} for key=" + full_key + " var=" + var);
                }
                if (!(stdv > 0.0) || !std::isfinite(stdv)) {
                    fatal("Non-positive or non-finite std for key=" + full_key + " var=" + var);
                }
                sc[var] = Cut1D{mean, stdv};
            }

            if (sc.empty()) {
                fatal("combined cuts JSON key \"" + full_key + "\" contains empty \"data\" cut map");
            }

            db.cuts[pk][topo] = sc;
        }
    }

    std::cout << "[propagator_study] Loaded combined cuts (new schema) from " << path << "\n";
    for (const auto& pk : required_period_keys) {
        std::cout << "  period_key=" << pk << ":\n";
        for (const auto& topo : required_topos) {
            const auto itp = db.cuts.find(pk);
            if (itp == db.cuts.end()) {
                fatal("Internal error: missing period_key after load: " + pk);
            }
            const auto itt = itp->second.find(topo);
            if (itt == itp->second.end()) {
                fatal("Internal error: missing topo after load: " + pk + " topo=" + topo);
            }
            std::cout << "    " << topo << " vars=" << itt->second.size() << "\n";
        }
    }

    return db;
}

static bool within_3sigma_value(double v, const Cut1D& c) {
    const double d = v - c.mean;
    return (std::fabs(d) <= 3.0 * c.std);
}

// -----------------------------------------------------------------------------
// CSV bin definitions and LUT
// -----------------------------------------------------------------------------

struct RowBin {
    Range xb;
    Range q2;
    Range t_abs;
    double phimin_deg = 0.0;
    double phimax_deg = 0.0;
    int xb_index = -1;
};

struct PhiBinRef {
    double phimin_deg = 0.0;
    double phimax_deg = 0.0;
    size_t row_index = 0;
};

struct LUT {
    std::vector<Range> xb_bins;

    struct XBNode {
        std::vector<Range> q2_bins;
        std::vector<Range> t_bins;
        std::vector<std::vector<std::vector<PhiBinRef>>> phi_bins; // [iq][it][iphi]
        int xb_index_for_name = -1;
    };

    std::vector<XBNode> nodes;
};

static LUT build_lut_from_bins(const std::vector<RowBin>& bins) {
    LUT lut;

    // Unique xB ranges
    std::set<Range> xb_set;
    for (const auto& b : bins) {
        xb_set.insert(b.xb);
    }
    lut.xb_bins.assign(xb_set.begin(), xb_set.end());
    lut.nodes.resize(lut.xb_bins.size());

    for (size_t ixb = 0; ixb < lut.xb_bins.size(); ++ixb) {
        const Range& xbR = lut.xb_bins[ixb];

        std::set<Range> q2_set;
        std::set<Range> t_set;
        int xb_idx_name = -1;

        for (const auto& b : bins) {
            if (b.xb == xbR) {
                q2_set.insert(b.q2);
                t_set.insert(b.t_abs);
                if (xb_idx_name < 0 && b.xb_index >= 0) {
                    xb_idx_name = b.xb_index;
                }
            }
        }

        LUT::XBNode node;
        node.q2_bins.assign(q2_set.begin(), q2_set.end());
        node.t_bins.assign(t_set.begin(), t_set.end());
        node.xb_index_for_name = xb_idx_name;

        const size_t nq = node.q2_bins.size();
        const size_t nt = node.t_bins.size();

        if (nq == 0 || nt == 0) {
            fatal("LUT build produced empty q2 or t bins for an xB slice");
        }

        node.phi_bins.resize(nq);
        for (size_t iq = 0; iq < nq; ++iq) {
            node.phi_bins[iq].resize(nt);
        }

        // The CSV binning can be ragged: not every (Q2,|t|) pair exists for a given xB.
        // Therefore we do NOT fatal on empty phi lists; we leave those cells empty.
        size_t n_empty_cells = 0;
        size_t n_filled_cells = 0;

        for (size_t iq = 0; iq < nq; ++iq) {
            for (size_t it = 0; it < nt; ++it) {
                const Range& q2R = node.q2_bins[iq];
                const Range& tR  = node.t_bins[it];

                std::vector<PhiBinRef> pb;
                for (size_t r = 0; r < bins.size(); ++r) {
                    const auto& b = bins[r];
                    if (b.xb == xbR && b.q2 == q2R && b.t_abs == tR) {
                        PhiBinRef ref;
                        ref.phimin_deg = b.phimin_deg;
                        ref.phimax_deg = b.phimax_deg;
                        ref.row_index  = r;
                        pb.push_back(ref);
                    }
                }

                if (pb.empty()) {
                    ++n_empty_cells;
                    continue;
                }

                std::sort(pb.begin(), pb.end(),
                          [](const PhiBinRef& a, const PhiBinRef& b) {
                              return wrap_phi_deg(a.phimin_deg) < wrap_phi_deg(b.phimin_deg);
                          });

                node.phi_bins[iq][it] = std::move(pb);
                ++n_filled_cells;
            }
        }

        std::cout << "[propagator_study] LUT xB slice ("
                  << std::fixed << std::setprecision(3) << xbR.first << ", " << xbR.second
                  << "): grid=(" << nq << " Q2 bins x " << nt << " |t| bins)"
                  << ", filled_cells=" << n_filled_cells
                  << ", empty_cells=" << n_empty_cells << "\n";

        lut.nodes[ixb] = std::move(node);
    }

    return lut;
}

static bool find_bin_index(const LUT& lut,
                           double xB,
                           double Q2,
                           double t_abs,
                           double phi_deg,
                           size_t& out_row_index) {
    int ixb = -1;
    for (size_t i = 0; i < lut.xb_bins.size(); ++i) {
        if (value_in_range(xB, lut.xb_bins[i])) {
            ixb = (int)i;
            break;
        }
    }
    if (ixb < 0) {
        return false;
    }

    const auto& node = lut.nodes[(size_t)ixb];

    int iq = -1;
    for (size_t i = 0; i < node.q2_bins.size(); ++i) {
        if (value_in_range(Q2, node.q2_bins[i])) {
            iq = (int)i;
            break;
        }
    }
    if (iq < 0) {
        return false;
    }

    int it = -1;
    for (size_t i = 0; i < node.t_bins.size(); ++i) {
        if (value_in_range(t_abs, node.t_bins[i])) {
            it = (int)i;
            break;
        }
    }
    if (it < 0) {
        return false;
    }

    const auto& pb = node.phi_bins[(size_t)iq][(size_t)it];
    if (pb.empty()) {
        return false;
    }

    phi_deg = wrap_phi_deg(phi_deg);

    for (const auto& ref : pb) {
        if (phi_in_range_deg(phi_deg, ref.phimin_deg, ref.phimax_deg)) {
            out_row_index = ref.row_index;
            return true;
        }
    }

    return false;
}

// -----------------------------------------------------------------------------
// Plotting helpers
// -----------------------------------------------------------------------------

struct Point {
    double phi = 0.0;
    double y = 0.0;
    double yerr = 0.0;
};

static TGraphErrors* make_graph(const std::vector<Point>& v) {
    if (v.empty()) {
        return nullptr;
    }

    const int N = (int)v.size();
    std::vector<double> x(N), y(N), ex(N), ey(N);
    for (int i = 0; i < N; ++i) {
        x[i]  = v[i].phi;
        y[i]  = v[i].y;
        ex[i] = 0.0;
        ey[i] = v[i].yerr;
    }

    TGraphErrors* g = new TGraphErrors(N, x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(1);
    g->SetLineColor(kBlack);
    g->SetMarkerColor(kBlack);
    return g;
}

static void set_root_style() {
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
}

static std::string fmt_cell_label(double q2min, double q2max, double tmin, double tmax) {
    std::ostringstream oss;
    oss << "Q^{2} in (" << std::fixed << std::setprecision(2) << q2min << ", " << q2max
        << "), |t| in (" << std::fixed << std::setprecision(2) << tmin << ", " << tmax << ")";
    return oss.str();
}

static void make_canvas_for_xb(const std::string& label,
                               const Range& xb_range,
                               const LUT& lut,
                               const std::vector<RowBin>& bins,
                               const std::vector<double>& N_noP1,
                               const std::vector<double>& N_withP1,
                               const fs::path& outdir,
                               int xb_idx_for_name,
                               double p1_threshold) {
    int ixb = -1;
    for (size_t i = 0; i < lut.xb_bins.size(); ++i) {
        if (lut.xb_bins[i] == xb_range) {
            ixb = (int)i;
            break;
        }
    }
    if (ixb < 0) {
        return;
    }

    const auto& node = lut.nodes[(size_t)ixb];

    const int ncols = (int)node.q2_bins.size();
    const int nrows = (int)node.t_bins.size();
    if (ncols <= 0 || nrows <= 0) {
        return;
    }

    const int nPads = ncols * nrows;

    int W = 300 * ncols + 160;
    int H = 260 * nrows + 240;
    if (W < 1200) {
        W = 1200;
    }
    if (H < 900) {
        H = 900;
    }

    double titleSize = 0.18;
    double legendTextSize = 0.11;
    double cellLabelSize = 0.070;

    if (nPads <= 4) {
        titleSize = 0.14;
        legendTextSize = 0.09;
        cellLabelSize = 0.060;
    }
    if (nPads == 1) {
        titleSize = 0.12;
        legendTextSize = 0.085;
        cellLabelSize = 0.055;
    }

    titleSize *= 0.5;
    legendTextSize *= 0.5;

    std::ostringstream cname;
    cname << "c_prop_" << canonical_period_dir(label) << "_xB" << xb_idx_for_name;

    TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

    TPad* pTop = new TPad("pTop", "pTop", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.78);
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
    head.SetTextSize(titleSize);

    std::ostringstream tit;
    tit << "Propagator cut efficiency (P1 > "
        << std::fixed << std::setprecision(3) << p1_threshold
        << "), ep #rightarrow ep#gamma   " << label
        << "   x_{B} in ("
        << std::fixed << std::setprecision(3)
        << xb_range.first << ", " << xb_range.second << ")";

    head.DrawLatex(0.5, 0.86, tit.str().c_str());

    TGraphErrors dummy;
    dummy.SetMarkerStyle(20);
    dummy.SetMarkerSize(1.0);
    dummy.SetLineWidth(1);
    dummy.SetMarkerColor(kBlack);
    dummy.SetLineColor(kBlack);

    TLegend* leg = new TLegend(0.02, 0.10, 0.40, 0.70);
    leg->SetBorderSize(1);
    leg->SetLineColor(kBlack);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(legendTextSize);
    leg->AddEntry(&dummy, "N(with P1 cut) / N(no P1 cut)", "lep");
    leg->Draw();

    for (int r = 0; r < nrows; ++r) {
        const Range& t_range = node.t_bins[(size_t)r];
        for (int cc = 0; cc < ncols; ++cc) {
            const Range& q2_range = node.q2_bins[(size_t)cc];

            pGrid->cd(r * ncols + cc + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.160);
            gPad->SetRightMargin(0.07);

            // Requested y-axis range:
            TH1* frame = gPad->DrawFrame(0.0, 0.4, 360.0, 1.6);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Fraction");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetXaxis()->SetLabelSize(0.048);
            frame->GetYaxis()->SetLabelSize(0.048);
            frame->GetXaxis()->SetTitleOffset(1.10);
            frame->GetYaxis()->SetTitleOffset(1.35);

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(cellLabelSize);
            lab.SetTextAlign(11);
            lab.SetTextFont(42);

            const std::string cell = fmt_cell_label(q2_range.first, q2_range.second,
                                                    t_range.first, t_range.second);
            lab.DrawLatex(0.12, 0.93, cell.c_str());

            const auto& pb = node.phi_bins[(size_t)cc][(size_t)r];
            if (pb.empty()) {
                continue;
            }

            std::vector<Point> pts;
            pts.reserve(pb.size());

            for (const auto& ref : pb) {
                const size_t ridx = ref.row_index;
                if (ridx >= bins.size()) {
                    fatal("Internal error: phi-bin row_index exceeds bins size");
                }

                const double denom = N_noP1[ridx];
                const double numer = N_withP1[ridx];

                if (denom <= 0.0) {
                    continue;
                }

                const double f = numer / denom;

                // Binomial uncertainty: sqrt(f(1-f)/N) with guard
                double err = 0.0;
                const double arg = f * (1.0 - f) / denom;
                if (arg > 0.0) {
                    err = std::sqrt(arg);
                }

                const double phi_c = phi_center_deg(ref.phimin_deg, ref.phimax_deg);

                // Requested: print out ratios for each plotted point
                std::cout << std::fixed << std::setprecision(6)
                          << "[propagator_study] POINT "
                          << "label=\"" << label << "\" "
                          << "xB=(" << xb_range.first << "," << xb_range.second << ") "
                          << "Q2=(" << q2_range.first << "," << q2_range.second << ") "
                          << "|t|=(" << t_range.first << "," << t_range.second << ") "
                          << "phi=" << phi_c << " "
                          << "numer=" << numer << " denom=" << denom
                          << " eff=" << f << " err=" << err
                          << "\n";

                Point p;
                p.phi  = phi_c;
                p.y    = f;
                p.yerr = err;
                pts.push_back(p);
            }

            std::sort(pts.begin(), pts.end(),
                      [](const Point& a, const Point& b) {
                          return a.phi < b.phi;
                      });

            TGraphErrors* g = make_graph(pts);
            if (g != nullptr) {
                g->Draw("P SAME");
            }

            TLine* l1 = new TLine(0.0, 1.0, 360.0, 1.0);
            l1->SetLineStyle(2);
            l1->SetLineWidth(1);
            l1->Draw("SAME");
        }
    }

    std::ostringstream fname;
    fname << "propagator_study_" << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    const fs::path outpath = outdir / fname.str();
    c->SaveAs(outpath.string().c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Tree reading / counting
// -----------------------------------------------------------------------------

struct BoundVarD {
    std::string name;
    std::unique_ptr<TTreeReaderValue<double>> val;
};

static bool tree_has_branch(TTree* t, const std::string& bname) {
    if (t == nullptr) {
        return false;
    }
    return (t->GetBranch(bname.c_str()) != nullptr);
}

static void require_branch(TTree* t, const std::string& bname, const std::string& context) {
    if (!tree_has_branch(t, bname)) {
        std::ostringstream oss;
        oss << "Missing required branch \"" << bname << "\" in tree "
            << (t ? t->GetName() : "(null)") << " (" << context << ")";
        fatal(oss.str());
    }
}

static void accumulate_counts_for_period(const std::string& period_key,
                                         const TreeVec& trees,
                                         const CutDB& cutdb,
                                         const LUT& lut,
                                         const std::vector<RowBin>& bins,
                                         std::vector<double>& N_noP1,
                                         std::vector<double>& N_withP1,
                                         const GlobalCutConfig& cfg_noP1,
                                         const GlobalCutConfig& cfg_withP1) {
    if (trees.empty()) {
        fatal("No trees provided for period_key=" + period_key);
    }

    for (TTree* tree : trees) {
        if (tree == nullptr) {
            continue;
        }

        // Required branches for binning and global cuts
        require_branch(tree, "x",              period_key);
        require_branch(tree, "Q2",             period_key);
        require_branch(tree, "t1",             period_key);
        require_branch(tree, "phi2",           period_key);
        require_branch(tree, "open_angle_ep2", period_key);
        require_branch(tree, "pTmiss",         period_key);
        require_branch(tree, "detector1",      period_key);
        require_branch(tree, "detector2",      period_key);
        require_branch(tree, "runnum",         period_key);

        // Required for P1 computation inside passes_global_cuts(...)
        require_branch(tree, "e_p",      period_key);
        require_branch(tree, "e_theta",  period_key);
        require_branch(tree, "e_phi",    period_key);
        require_branch(tree, "p2_p",     period_key);
        require_branch(tree, "p2_theta", period_key);
        require_branch(tree, "p2_phi",   period_key);

        // Optional exclusivity variables for 3-sigma cuts
        const std::vector<std::string> opt_vars = {
            "Emiss2",
            "Mx2",
            "Mx2_1",
            "Mx2_2",
            "theta_gamma_gamma",
            "xF",
            "pTmiss"
        };

        TTreeReader reader(tree);

        TTreeReaderValue<double> x(reader, "x");
        TTreeReaderValue<double> Q2(reader, "Q2");
        TTreeReaderValue<double> t1(reader, "t1");
        TTreeReaderValue<double> phi2(reader, "phi2");
        TTreeReaderValue<double> open_angle_ep2(reader, "open_angle_ep2");
        TTreeReaderValue<double> pTmiss(reader, "pTmiss");
        TTreeReaderValue<int> detector1(reader, "detector1");
        TTreeReaderValue<int> detector2(reader, "detector2");
        TTreeReaderValue<int> runnum(reader, "runnum");

        TTreeReaderValue<double> e_p(reader, "e_p");
        TTreeReaderValue<double> e_theta(reader, "e_theta");
        TTreeReaderValue<double> e_phi(reader, "e_phi");
        TTreeReaderValue<double> p2_p(reader, "p2_p");
        TTreeReaderValue<double> p2_theta(reader, "p2_theta");
        TTreeReaderValue<double> p2_phi(reader, "p2_phi");

        std::vector<BoundVarD> opt_bound;
        for (const auto& v : opt_vars) {
            if (tree_has_branch(tree, v)) {
                BoundVarD b;
                b.name = v;
                b.val.reset(new TTreeReaderValue<double>(reader, v.c_str()));
                opt_bound.push_back(std::move(b));
            }
        }

        const size_t n_bins = bins.size();
        if (N_noP1.size() != n_bins || N_withP1.size() != n_bins) {
            fatal("Internal error: count arrays size mismatch vs bins");
        }

        // Diagnostics for UNKNOWN topology: count observed (d1,d2) pairs that map to UNKNOWN
        std::map<std::pair<int,int>, long long> unknown_pairs;
        long long n_unknown_topo = 0;

        while (reader.Next()) {
            const int rnum = *runnum;
            if (is_excluded_run(rnum, cfg_withP1)) {
                continue;
            }

            const double t_abs = -(*t1);
            if (!(t_abs > 0.0) || !std::isfinite(t_abs)) {
                continue;
            }

            double phi_deg = (*phi2) * 180.0 / kPi;
            phi_deg = wrap_phi_deg(phi_deg);

            bool pass_noP1 = false;
            bool pass_withP1 = false;

            try {
                pass_noP1 = passes_global_cuts(*t1,
                                               *open_angle_ep2,
                                               *pTmiss,
                                               period_key,
                                               *e_p, *e_theta, *e_phi,
                                               *p2_p, *p2_theta, *p2_phi,
                                               cfg_noP1);

                if (!pass_noP1) {
                    continue;
                }

                pass_withP1 = passes_global_cuts(*t1,
                                                 *open_angle_ep2,
                                                 *pTmiss,
                                                 period_key,
                                                 *e_p, *e_theta, *e_phi,
                                                 *p2_p, *p2_theta, *p2_phi,
                                                 cfg_withP1);
            } catch (const std::exception& e) {
                std::ostringstream oss;
                oss << "passes_global_cuts threw exception for period_key=" << period_key
                    << " tree=" << tree->GetName()
                    << " (what=" << e.what() << ")";
                fatal(oss.str());
            }

            // Topology and 3-sigma cuts are only relevant after baseline selection passes.
            const int d1 = *detector1;
            const int d2 = *detector2;
            const std::string topo = topology_from_detectors(d1, d2);
            if (topo == "UNKNOWN") {
                ++n_unknown_topo;
                unknown_pairs[std::make_pair(d1, d2)] += 1;
                continue;
            }

            auto itp = cutdb.cuts.find(period_key);
            if (itp == cutdb.cuts.end()) {
                fatal("CutDB missing period_key=\"" + period_key + "\"");
            }
            auto itt = itp->second.find(topo);
            if (itt == itp->second.end()) {
                fatal("CutDB missing topo=\"" + topo + "\" for period_key=\"" + period_key + "\"");
            }
            const SigmaCut& sc = itt->second;

            bool pass_3s = true;
            for (const auto& bv : opt_bound) {
                auto itv = sc.find(bv.name);
                if (itv == sc.end()) {
                    continue;
                }
                const double vv = **(bv.val);
                if (!within_3sigma_value(vv, itv->second)) {
                    pass_3s = false;
                    break;
                }
            }
            if (!pass_3s) {
                continue;
            }

            size_t ridx = 0;
            if (!find_bin_index(lut, *x, *Q2, t_abs, phi_deg, ridx)) {
                continue;
            }

            N_noP1[ridx] += 1.0;
            if (pass_withP1) {
                N_withP1[ridx] += 1.0;
            }
        }

        if (n_unknown_topo > 0) {
            std::cout << "[propagator_study] INFO: period_key=" << period_key
                      << " tree=" << tree->GetName()
                      << " skipped " << n_unknown_topo
                      << " events with UNKNOWN topology (after baseline noP1 cuts)\n";

            // Print the most common unknown (d1,d2) pairs (top 10)
            std::vector<std::pair<std::pair<int,int>, long long>> v;
            v.reserve(unknown_pairs.size());
            for (const auto& kv : unknown_pairs) {
                v.push_back(kv);
            }
            std::sort(v.begin(), v.end(),
                      [](const auto& a, const auto& b) {
                          return a.second > b.second;
                      });

            const size_t nprint = std::min<size_t>(10, v.size());
            std::cout << "[propagator_study] INFO: Most common UNKNOWN (detector1, detector2) pairs:\n";
            for (size_t i = 0; i < nprint; ++i) {
                std::cout << "  (" << v[i].first.first << ", " << v[i].first.second << ") : " << v[i].second << "\n";
            }
        }
    }
}

// -----------------------------------------------------------------------------
// CSV bin loader
// -----------------------------------------------------------------------------

static std::vector<RowBin> load_bins_from_csv(const std::string& csv_main) {
    std::ifstream ifs(csv_main);
    if (!ifs) {
        fatal("Cannot open CSV: " + csv_main);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }
    ifs.close();

    if (lines.empty()) {
        fatal("CSV is empty: " + csv_main);
    }

    const std::vector<std::string> header = split_csv_line(lines[0]);

    const int c_xb_min  = find_col(header, "xBmin");
    const int c_xb_max  = find_col(header, "xBmax");
    const int c_q2_min  = find_col(header, "Q2min");
    const int c_q2_max  = find_col(header, "Q2max");
    const int c_t_min   = find_col(header, "t_abs_min");
    const int c_t_max   = find_col(header, "t_abs_max");
    const int c_phimin  = find_col(header, "phimin");
    const int c_phimax  = find_col(header, "phimax");
    const int c_xb_idx  = find_col_optional(header, "xB index");

    std::vector<RowBin> bins;
    bins.reserve(lines.size() > 1 ? lines.size() - 1 : 0);

    for (size_t r = 1; r < lines.size(); ++r) {
        if (lines[r].empty()) {
            continue;
        }
        const std::vector<std::string> fields = split_csv_line(lines[r]);
        if (fields.size() != header.size()) {
            continue;
        }

        RowBin b;

        b.xb.first  = std::atof(trim(unquote(fields[c_xb_min])).c_str());
        b.xb.second = std::atof(trim(unquote(fields[c_xb_max])).c_str());
        b.q2.first  = std::atof(trim(unquote(fields[c_q2_min])).c_str());
        b.q2.second = std::atof(trim(unquote(fields[c_q2_max])).c_str());
        b.t_abs.first  = std::atof(trim(unquote(fields[c_t_min])).c_str());
        b.t_abs.second = std::atof(trim(unquote(fields[c_t_max])).c_str());

        b.phimin_deg = std::atof(trim(unquote(fields[c_phimin])).c_str());
        b.phimax_deg = std::atof(trim(unquote(fields[c_phimax])).c_str());

        if (c_xb_idx >= 0) {
            b.xb_index = std::atoi(trim(unquote(fields[c_xb_idx])).c_str());
        }

        if (!(b.xb.first < b.xb.second)) {
            continue;
        }
        if (!(b.q2.first < b.q2.second)) {
            continue;
        }
        if (!(b.t_abs.first < b.t_abs.second)) {
            continue;
        }

        bins.push_back(b);
    }

    if (bins.empty()) {
        fatal("No valid bin rows parsed from CSV: " + csv_main);
    }

    std::cout << "[propagator_study] Loaded " << bins.size()
              << " bin-rows from CSV: " << csv_main << "\n";

    return bins;
}

static void require_period_keys(const TreeMap& m,
                                const std::vector<std::string>& keys) {
    for (const auto& k : keys) {
        auto it = m.find(k);
        if (it == m.end()) {
            fatal("data_trees_by_period missing required key: \"" + k + "\"");
        }
        if (it->second.empty()) {
            fatal("data_trees_by_period has empty tree vector for key: \"" + k + "\"");
        }
    }
}

// -----------------------------------------------------------------------------
// Public entry point
// -----------------------------------------------------------------------------

bool run_propagator_study(const std::string& csv_main,
                          const TreeMap& data_trees_by_period,
                          const std::string& combined_cuts_json,
                          const std::string& out_root_dir) {
    set_root_style();

    const std::vector<std::string> base_period_keys = {
        "fa18_inb", "fa18_out", "sp18_inb", "sp18_out", "sp19_inb"
    };
    require_period_keys(data_trees_by_period, base_period_keys);

    const std::vector<RowBin> bins = load_bins_from_csv(csv_main);
    const LUT lut = build_lut_from_bins(bins);

    const CutDB cutdb = load_combined_cuts_json_newschema(combined_cuts_json, base_period_keys);

    GlobalCutConfig cfg_withP1 = default_global_cuts();
    GlobalCutConfig cfg_noP1   = default_global_cuts();

    // Disable P1/ycol-mirror cut in cfg_noP1 only.
    cfg_noP1.enable_dvcsgen_ycol_cut = false;

    const double p1_threshold = cfg_withP1.dvcsgen_ycol_cut;

    struct LabelSpec {
        std::string label;
        std::vector<std::string> members;
    };

    const std::vector<LabelSpec> labels = {
        {"Fa18 Inb", {"fa18_inb"}},
        {"Fa18 Out", {"fa18_out"}},
        {"Sp18 Inb", {"sp18_inb"}},
        {"Sp18 Out", {"sp18_out"}},
        {"Sp19 Inb", {"sp19_inb"}},
        {"Fa18",     {"fa18_inb", "fa18_out"}},
        {"Sp18",     {"sp18_inb", "sp18_out"}},
        {"10.6 GeV", {"fa18_inb", "fa18_out", "sp18_inb", "sp18_out"}}
    };

    for (const auto& L : labels) {
        std::vector<double> N_noP1(bins.size(), 0.0);
        std::vector<double> N_withP1(bins.size(), 0.0);

        std::cout << "[propagator_study] Accumulating counts for label=\""
                  << L.label << "\" using members={";
        for (size_t i = 0; i < L.members.size(); ++i) {
            if (i > 0) {
                std::cout << ", ";
            }
            std::cout << L.members[i];
        }
        std::cout << "}\n";

        for (const auto& period_key : L.members) {
            auto it = data_trees_by_period.find(period_key);
            if (it == data_trees_by_period.end()) {
                fatal("Missing member period_key=\"" + period_key + "\" for label=\"" + L.label + "\"");
            }

            accumulate_counts_for_period(period_key,
                                         it->second,
                                         cutdb,
                                         lut,
                                         bins,
                                         N_noP1,
                                         N_withP1,
                                         cfg_noP1,
                                         cfg_withP1);
        }

        const fs::path outdir = fs::path(out_root_dir) / canonical_period_dir(L.label);
        ensure_dir(outdir);

        int xb_canvas_counter = 0;
        for (size_t ixb = 0; ixb < lut.xb_bins.size(); ++ixb) {
            const Range& xbR = lut.xb_bins[ixb];
            const auto& node = lut.nodes[ixb];

            int xb_idx_for_name = node.xb_index_for_name;
            if (xb_idx_for_name < 0) {
                xb_idx_for_name = xb_canvas_counter;
            }

            make_canvas_for_xb(L.label,
                               xbR,
                               lut,
                               bins,
                               N_noP1,
                               N_withP1,
                               outdir,
                               xb_idx_for_name,
                               p1_threshold);

            ++xb_canvas_counter;
        }

        std::cout << "[propagator_study] Wrote plots for label=\""
                  << L.label << "\" into " << outdir.string() << "\n";
    }

    return true;
}

} // namespace propagator_study