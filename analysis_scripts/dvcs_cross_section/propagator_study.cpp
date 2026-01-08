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
    std::string out = label;
    std::replace(out.begin(), out.end(), ' ', '_');
    return out;
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
        char c = line[i];
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
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
}

static std::string unquote(const std::string& s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        std::string inner = s.substr(1, s.size() - 2);
        std::string out;
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
    if (x < 0.0) x += 360.0;
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
    // wrap-around: accept if phi >= phimin OR phi < phimax
    return (phi_deg >= phimin || phi_deg < phimax);
}

static double phi_center_deg(double phimin, double phimax) {
    double a = wrap_phi_deg(phimin);
    double b = wrap_phi_deg(phimax);

    if (a <= b) {
        return 0.5 * (a + b);
    }
    // wrap-around: treat b as b+360
    double bb = b + 360.0;
    double c = 0.5 * (a + bb);
    return wrap_phi_deg(c);
}

static bool value_in_range(double v, const Range& r) {
    return (v >= r.first && v < r.second);
}

static std::string topology_from_detectors(int detector1, int detector2) {
    // Deterministic mapping:
    // 1 = FD, 2 = CD, 3 = FT (assumed consistent with existing analysis)
    if (detector1 == 1 && detector2 == 1) return "FD_FD";
    if (detector1 == 2 && detector2 == 1) return "CD_FD";
    if (detector1 == 2 && detector2 == 3) return "CD_FT";
    return "UNKNOWN";
}

// -----------------------------------------------------------------------------
// 3-sigma cut map loader (combined_cuts.json)
// -----------------------------------------------------------------------------

struct Cut1D {
    double mu = 0.0;
    double sigma = 0.0;
};

using SigmaCut   = std::map<std::string, Cut1D>;
using TopoCutMap = std::map<std::string, SigmaCut>;

static bool read_mu_sigma(const json& j, double& mu, double& sigma) {
    // Accept a small set of deterministic key spellings.
    // This is not a fallback for analysis logic; it is a strict parser for the cut file schema.
    if (!j.is_object()) return false;

    bool have_mu = false;
    bool have_sigma = false;

    if (j.contains("mu")) {
        mu = j["mu"].get<double>();
        have_mu = true;
    } else if (j.contains("mean")) {
        mu = j["mean"].get<double>();
        have_mu = true;
    }

    if (j.contains("sigma")) {
        sigma = j["sigma"].get<double>();
        have_sigma = true;
    } else if (j.contains("std")) {
        sigma = j["std"].get<double>();
        have_sigma = true;
    }

    return have_mu && have_sigma;
}

static TopoCutMap load_combined_cuts_json(const std::string& path) {
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

    // Expect top-level to contain topology objects (FD_FD, CD_FD, CD_FT).
    TopoCutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string topo = it.key();
        const json& topo_obj = it.value();
        if (!topo_obj.is_object()) continue;

        SigmaCut sc;
        for (auto iv = topo_obj.begin(); iv != topo_obj.end(); ++iv) {
            const std::string var = iv.key();
            const json& spec = iv.value();
            double mu = 0.0, sigma = 0.0;
            if (!read_mu_sigma(spec, mu, sigma)) {
                continue;
            }
            if (!(sigma > 0.0) || !std::isfinite(sigma)) {
                fatal("Non-positive or non-finite sigma for topo=" + topo + " var=" + var);
            }
            sc[var] = Cut1D{mu, sigma};
        }

        if (!sc.empty()) {
            out[topo] = sc;
        }
    }

    if (out.empty()) {
        fatal("combined cuts JSON produced an empty TopoCutMap: " + path);
    }

    // Require the canonical DVCS topologies to exist.
    const std::vector<std::string> required_topos = {"FD_FD", "CD_FD", "CD_FT"};
    for (const auto& t : required_topos) {
        if (out.find(t) == out.end()) {
            fatal("combined cuts JSON missing required topology key: " + t);
        }
    }

    std::cout << "[propagator_study] Loaded combined cuts from " << path
              << " with topologies: ";
    for (auto it = out.begin(); it != out.end(); ++it) {
        if (it != out.begin()) std::cout << ", ";
        std::cout << it->first << "(" << it->second.size() << " vars)";
    }
    std::cout << "\n";

    return out;
}

static bool within_3sigma_value(double v, const Cut1D& c) {
    const double d = v - c.mu;
    return (std::fabs(d) <= 3.0 * c.sigma);
}

// -----------------------------------------------------------------------------
// CSV bin definitions
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
    size_t row_index = 0; // index into bins vector
};

struct LUT {
    std::vector<Range> xb_bins;

    struct XBNode {
        std::vector<Range> q2_bins;
        std::vector<Range> t_bins;
        // phi bins per (q2_idx, t_idx)
        std::vector<std::vector<std::vector<PhiBinRef>>> phi_bins;
        int xb_index_for_name = -1;
    };

    std::vector<XBNode> nodes;
};

static LUT build_lut_from_bins(const std::vector<RowBin>& bins) {
    LUT lut;

    // Unique xB ranges
    std::set<Range> xb_set;
    for (const auto& b : bins) xb_set.insert(b.xb);
    lut.xb_bins.assign(xb_set.begin(), xb_set.end());
    lut.nodes.resize(lut.xb_bins.size());

    for (size_t ixb = 0; ixb < lut.xb_bins.size(); ++ixb) {
        const Range& xbR = lut.xb_bins[ixb];

        std::set<Range> q2_set, t_set;
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
        node.phi_bins.resize(nq);
        for (size_t iq = 0; iq < nq; ++iq) {
            node.phi_bins[iq].resize(nt);
        }

        // Fill phi bins for each (q2,t)
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

                std::sort(pb.begin(), pb.end(),
                          [](const PhiBinRef& a, const PhiBinRef& b) {
                              return wrap_phi_deg(a.phimin_deg) < wrap_phi_deg(b.phimin_deg);
                          });

                node.phi_bins[iq][it] = std::move(pb);
            }
        }

        lut.nodes[ixb] = std::move(node);
    }

    return lut;
}

static bool find_bin_index(const LUT& lut,
                           const std::vector<RowBin>& bins,
                           double xB,
                           double Q2,
                           double t_abs,
                           double phi_deg,
                           size_t& out_row_index) {
    (void)bins;

    // xB
    int ixb = -1;
    for (size_t i = 0; i < lut.xb_bins.size(); ++i) {
        if (value_in_range(xB, lut.xb_bins[i])) {
            ixb = (int)i;
            break;
        }
    }
    if (ixb < 0) return false;

    const auto& node = lut.nodes[(size_t)ixb];

    int iq = -1;
    for (size_t i = 0; i < node.q2_bins.size(); ++i) {
        if (value_in_range(Q2, node.q2_bins[i])) {
            iq = (int)i;
            break;
        }
    }
    if (iq < 0) return false;

    int it = -1;
    for (size_t i = 0; i < node.t_bins.size(); ++i) {
        if (value_in_range(t_abs, node.t_bins[i])) {
            it = (int)i;
            break;
        }
    }
    if (it < 0) return false;

    const auto& pb = node.phi_bins[(size_t)iq][(size_t)it];
    if (pb.empty()) return false;

    for (const auto& ref : pb) {
        if (phi_in_range_deg(phi_deg, ref.phimin_deg, ref.phimax_deg)) {
            out_row_index = ref.row_index;
            return true;
        }
    }

    return false;
}

// -----------------------------------------------------------------------------
// Plotting
// -----------------------------------------------------------------------------

struct Point {
    double phi = 0.0;
    double y = 0.0;
    double yerr = 0.0;
};

static TGraphErrors* make_graph(const std::vector<Point>& v) {
    if (v.empty()) return nullptr;
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

static void make_canvas_for_xb(const std::string& label,
                               const Range& xb_range,
                               const LUT& lut,
                               const std::vector<RowBin>& bins,
                               const std::vector<double>& N_noP1,
                               const std::vector<double>& N_withP1,
                               const fs::path& outdir,
                               int xb_idx_for_name,
                               double p1_threshold) {
    // Identify q2 and t slices for this xB
    int ixb = -1;
    for (size_t i = 0; i < lut.xb_bins.size(); ++i) {
        if (lut.xb_bins[i] == xb_range) {
            ixb = (int)i;
            break;
        }
    }
    if (ixb < 0) return;

    const auto& node = lut.nodes[(size_t)ixb];
    const int ncols = (int)node.q2_bins.size();
    const int nrows = (int)node.t_bins.size();
    if (ncols <= 0 || nrows <= 0) return;

    const int nPads = ncols * nrows;

    int W = 300 * ncols + 160;
    int H = 260 * nrows + 240;
    if (W < 1200) W = 1200;
    if (H < 900)  H = 900;

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

    // Legend in top pad
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

    // Draw pads
    for (int r = 0; r < nrows; ++r) {
        const Range& t_range = node.t_bins[(size_t)r];
        for (int cc = 0; cc < ncols; ++cc) {
            const Range& q2_range = node.q2_bins[(size_t)cc];

            pGrid->cd(r * ncols + cc + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.10);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 1.05);
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
            lab.DrawLatex(
                0.14, 0.93,
                Form("Q^{2} in (%.2f, %.2f), |t| in (%.2f, %.2f)",
                     q2_range.first, q2_range.second,
                     t_range.first,  t_range.second)
            );

            // Collect phi points in this (xB,q2,t)
            const auto& pb = node.phi_bins[(size_t)cc][(size_t)r];
            if (pb.empty()) continue;

            std::vector<Point> pts;
            pts.reserve(pb.size());

            for (const auto& ref : pb) {
                const size_t ridx = ref.row_index;
                if (ridx >= bins.size()) continue;

                const double denom = N_noP1[ridx];
                const double numer = N_withP1[ridx];

                if (denom <= 0.0) continue;

                double f = numer / denom;
                if (f < 0.0) f = 0.0;
                if (f > 1.0) f = 1.0;

                double err = 0.0;
                // Binomial uncertainty: sqrt(f(1-f)/N)
                err = std::sqrt(std::max(0.0, f * (1.0 - f) / denom));

                Point p;
                p.phi  = phi_center_deg(ref.phimin_deg, ref.phimax_deg);
                p.y    = f;
                p.yerr = err;
                pts.push_back(p);
            }

            std::sort(pts.begin(), pts.end(),
                      [](const Point& a, const Point& b) { return a.phi < b.phi; });

            TGraphErrors* g = make_graph(pts);
            if (g) g->Draw("P SAME");

            // Reference line at 1.0
            TLine* l1 = new TLine(0.0, 1.0, 360.0, 1.0);
            l1->SetLineStyle(2);
            l1->SetLineWidth(1);
            l1->Draw("SAME");
        }
    }

    std::ostringstream fname;
    fname << "propagator_study_" << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    fs::path outpath = outdir / fname.str();
    c->SaveAs(outpath.string().c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Counting pass (one period key, one set of trees)
// -----------------------------------------------------------------------------

struct BoundVarD {
    std::string name;
    std::unique_ptr<TTreeReaderValue<double>> val;
};

static bool tree_has_branch(TTree* t, const std::string& bname) {
    if (!t) return false;
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
                                         const TopoCutMap& cutmap,
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
        if (!tree) continue;

        // Required branches for binning and global cuts
        require_branch(tree, "x",               period_key);
        require_branch(tree, "Q2",              period_key);
        require_branch(tree, "t1",              period_key);
        require_branch(tree, "phi2",            period_key);
        require_branch(tree, "open_angle_ep2",  period_key);
        require_branch(tree, "pTmiss",          period_key);
        require_branch(tree, "detector1",       period_key);
        require_branch(tree, "detector2",       period_key);
        require_branch(tree, "runnum",          period_key);

        // Required for P1 computation inside passes_global_cuts(...)
        require_branch(tree, "e_p",      period_key);
        require_branch(tree, "e_theta",  period_key);
        require_branch(tree, "e_phi",    period_key);
        require_branch(tree, "p2_p",     period_key);
        require_branch(tree, "p2_theta", period_key);
        require_branch(tree, "p2_phi",   period_key);

        // Optional exclusivity variables for 3-sigma (apply if both present in tree and in cutmap for the topology)
        const std::vector<std::string> opt_vars = {
            "Emiss2",
            "Mx2",
            "Mx2_1",
            "Mx2_2",
            "theta_gamma_gamma",
            "xF",
            "pTmiss" // already required, but also commonly in cutmap
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

        size_t n_unknown_topo = 0;

        while (reader.Next()) {
            const int rnum = *runnum;
            if (is_excluded_run(rnum, cfg_withP1)) {
                continue;
            }

            const int d1 = *detector1;
            const int d2 = *detector2;
            const std::string topo = topology_from_detectors(d1, d2);
            if (topo == "UNKNOWN") {
                ++n_unknown_topo;
                continue;
            }

            auto it_topo = cutmap.find(topo);
            if (it_topo == cutmap.end()) {
                fatal("Cut map missing topology \"" + topo + "\"");
            }
            const SigmaCut& sc = it_topo->second;

            // 3-sigma cuts: apply for variables that are both bound and present in sc
            bool pass_3s = true;
            for (const auto& bv : opt_bound) {
                auto itv = sc.find(bv.name);
                if (itv == sc.end()) continue;
                const double val = **(bv.val);
                if (!within_3sigma_value(val, itv->second)) {
                    pass_3s = false;
                    break;
                }
            }
            if (!pass_3s) continue;

            const double t_abs = -(*t1);
            if (!(t_abs > 0.0) || !std::isfinite(t_abs)) continue;

            double phi_deg = (*phi2) * 180.0 / M_PI;
            phi_deg = wrap_phi_deg(phi_deg);

            // Global cuts without P1
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

                if (!pass_noP1) continue;

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

            // Bin lookup
            size_t ridx = 0;
            if (!find_bin_index(lut, bins, *x, *Q2, t_abs, phi_deg, ridx)) {
                continue;
            }

            // Denominator always increments here
            N_noP1[ridx] += 1.0;
            if (pass_withP1) {
                N_withP1[ridx] += 1.0;
            }
        }

        if (n_unknown_topo > 0) {
            std::cout << "[propagator_study] INFO: period_key=" << period_key
                      << " tree=" << tree->GetName()
                      << " skipped " << n_unknown_topo
                      << " events with UNKNOWN topology\n";
        }
    } // end tree loop
}

// -----------------------------------------------------------------------------
// Main entry point
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

    std::vector<std::string> header = split_csv_line(lines[0]);

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
        if (lines[r].empty()) continue;
        std::vector<std::string> fields = split_csv_line(lines[r]);
        if (fields.size() != header.size()) continue;

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

        if (!(b.xb.first < b.xb.second)) continue;
        if (!(b.q2.first < b.q2.second)) continue;
        if (!(b.t_abs.first < b.t_abs.second)) continue;

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

bool run_propagator_study(const std::string& csv_main,
                          const TreeMap& data_trees_by_period,
                          const std::string& combined_cuts_json,
                          const std::string& out_root_dir) {
    set_root_style();

    // Require the five base periods as inputs
    const std::vector<std::string> base_period_keys = {
        "fa18_inb", "fa18_out", "sp18_inb", "sp18_out", "sp19_inb"
    };
    require_period_keys(data_trees_by_period, base_period_keys);

    // Load bins and LUT
    const std::vector<RowBin> bins = load_bins_from_csv(csv_main);
    const LUT lut = build_lut_from_bins(bins);

    // Load 3-sigma cuts
    const TopoCutMap cutmap = load_combined_cuts_json(combined_cuts_json);

    // Prepare global cut configs:
    //   cfg_noP1: same as default but P1 cut disabled
    //   cfg_withP1: default (P1 cut enabled)
    GlobalCutConfig cfg_withP1 = default_global_cuts();
    GlobalCutConfig cfg_noP1   = default_global_cuts();
    cfg_noP1.enable_dvcsgen_ycol_cut = false;

    const double p1_threshold = cfg_withP1.dvcsgen_ycol_cut;

    // Labels to produce (display labels)
    struct LabelSpec {
        std::string label;
        std::vector<std::string> members; // period keys
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

    // For each label, accumulate counts and then plot
    for (const auto& L : labels) {
        std::vector<double> N_noP1(bins.size(), 0.0);
        std::vector<double> N_withP1(bins.size(), 0.0);

        std::cout << "[propagator_study] Accumulating counts for label=\""
                  << L.label << "\" using members={";
        for (size_t i = 0; i < L.members.size(); ++i) {
            if (i > 0) std::cout << ", ";
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
                                         cutmap,
                                         lut,
                                         bins,
                                         N_noP1,
                                         N_withP1,
                                         cfg_noP1,
                                         cfg_withP1);
        }

        // Output directory for this label
        fs::path outdir = fs::path(out_root_dir) / canonical_period_dir(L.label);
        ensure_dir(outdir);

        // Iterate xB canvases and make plots
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