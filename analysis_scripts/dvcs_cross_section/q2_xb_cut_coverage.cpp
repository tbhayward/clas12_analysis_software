#include "q2_xb_cut_coverage.h"

#include "exclusivity_cuts.h"
#include "global_cuts.h"

#include <nlohmann/json.hpp>

#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TDataType.h>
#include <TVirtualPad.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

static constexpr double kXBMin = 0.0;
static constexpr double kXBMax = 0.7;
static constexpr double kQ2Min = 0.0;
static constexpr double kQ2Max = 7.0;

static constexpr int kInbendingColor = kRed + 1;
static constexpr int kOutbendingColor = kBlue + 1;
static constexpr int kBinLineColor = kSpring - 5;

static std::string topo_to_key(Topology topo) {
    switch (topo) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }

    throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Unknown topology enum value.");
}

static bool event_topology_from_detectors(int detector1, int detector2, Topology& topo) {
    if (detector1 == 1 && detector2 == 1) {
        topo = Topology::FD_FD;
        return true;
    }
    // endif

    if (detector1 == 2 && detector2 == 1) {
        topo = Topology::CD_FD;
        return true;
    }
    // endif

    if (detector1 == 2 && detector2 == 0) {
        topo = Topology::CD_FT;
        return true;
    }
    // endif

    return false;
}

static std::string period_code_dvcs_from_label(const std::string& run_tag) {
    std::string nice = run_tag;

    if (!nice.empty()) {
        nice[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(nice[0])));
    }
    // endif

    for (size_t i = 0; i + 1 < nice.size(); ++i) {
        if (nice[i] == '_') {
            nice[i + 1] = static_cast<char>(std::toupper(static_cast<unsigned char>(nice[i + 1])));
        }
        // endif
    }
    // endfor

    return std::string("DVCS_") + nice;
}

static bool within_3sigma(double val, const Stats& stats) {
    return (val >= stats.mean - 3.0 * stats.std) && (val <= stats.mean + 3.0 * stats.std);
}

static bool passes_3sigma_cuts(
    const std::map<std::string, Stats>& cuts,
    const std::map<std::string, double>& values
) {
    for (const auto& kv : cuts) {
        const std::string& var = kv.first;
        const Stats& stats = kv.second;

        const auto it = values.find(var);
        if (it == values.end()) {
            continue;
        }
        // endif

        if (!within_3sigma(it->second, stats)) {
            return false;
        }
        // endif
    }
    // endfor

    return true;
}

static std::vector<std::string> parse_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
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
            // endif
        } else if (c == ',' && !in_quotes) {
            fields.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
        // endif
    }
    // endfor

    fields.push_back(cur);
    return fields;
}

static double parse_required_double(const std::string& text, const std::string& context) {
    if (text.empty()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Empty numeric field while parsing " + context + ".");
    }
    // endif

    char* endptr = nullptr;
    const double val = std::strtod(text.c_str(), &endptr);

    if (endptr == text.c_str()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Could not parse numeric field '" + text + "' while parsing " + context + ".");
    }
    // endif

    return val;
}

static int find_required_column(const std::vector<std::string>& header, const std::string& name) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) {
            return static_cast<int>(i);
        }
        // endif
    }
    // endfor

    throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Required CSV column is missing: " + name);
}

struct BinGrid {
    std::vector<double> xb_edges;
    std::vector<double> q2_edges;
};

static void add_unique_edge(std::vector<double>& edges, double val) {
    const double rounded = std::round(val * 1000000.0) / 1000000.0;

    for (double existing : edges) {
        if (std::abs(existing - rounded) < 1.0e-9) {
            return;
        }
        // endif
    }
    // endfor

    edges.push_back(rounded);
}

static BinGrid load_bin_grid_from_csv(const std::string& csv_path) {
    std::ifstream fin(csv_path);
    if (!fin) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Could not open CSV: " + csv_path);
    }
    // endif

    std::string line;
    if (!std::getline(fin, line)) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: CSV is empty: " + csv_path);
    }
    // endif

    const std::vector<std::string> header = parse_csv_line(line);

    const int ixmin = find_required_column(header, "xBmin");
    const int ixmax = find_required_column(header, "xBmax");
    const int iqmin = find_required_column(header, "Q2min");
    const int iqmax = find_required_column(header, "Q2max");

    BinGrid grid;

    Long64_t row_index = 1;
    while (std::getline(fin, line)) {
        ++row_index;

        if (line.empty()) {
            continue;
        }
        // endif

        const std::vector<std::string> row = parse_csv_line(line);

        const int needed = std::max(std::max(ixmin, ixmax), std::max(iqmin, iqmax));
        if (static_cast<int>(row.size()) <= needed) {
            std::ostringstream ss;
            ss << "[q2_xb_cut_coverage] FATAL: CSV row " << row_index
               << " has too few columns.";
            throw std::runtime_error(ss.str());
        }
        // endif

        const double xbmin = parse_required_double(row[ixmin], "xBmin");
        const double xbmax = parse_required_double(row[ixmax], "xBmax");
        const double q2min = parse_required_double(row[iqmin], "Q2min");
        const double q2max = parse_required_double(row[iqmax], "Q2max");

        add_unique_edge(grid.xb_edges, xbmin);
        add_unique_edge(grid.xb_edges, xbmax);
        add_unique_edge(grid.q2_edges, q2min);
        add_unique_edge(grid.q2_edges, q2max);
    }

    std::sort(grid.xb_edges.begin(), grid.xb_edges.end());
    std::sort(grid.q2_edges.begin(), grid.q2_edges.end());

    if (grid.xb_edges.empty()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: No xB bin edges found in CSV.");
    }
    // endif

    if (grid.q2_edges.empty()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: No Q2 bin edges found in CSV.");
    }
    // endif

    return grid;
}

struct BranchBinder {
    int runnum = 0;
    int detector1 = 0;
    int detector2 = 0;

    double x = 0.0;
    double Q2 = 0.0;

    double t1 = 0.0;
    double open_angle_ep2 = 0.0;
    double pTmiss = 0.0;

    double Emiss2 = 0.0;
    double Mx2 = 0.0;
    double Mx2_1 = 0.0;
    double Mx2_2 = 0.0;
    double xF = 0.0;
    double Delta_phi = 0.0;
    double theta_gamma_gamma = 0.0;

    double e_p = 0.0;
    double e_theta = 0.0;
    double e_phi = 0.0;

    double p2_p = 0.0;
    double p2_theta = 0.0;
    double p2_phi = 0.0;

    bool has_runnum = false;
    bool has_detector1 = false;
    bool has_detector2 = false;

    bool has_x = false;
    bool has_Q2 = false;

    bool has_t1 = false;
    bool has_open_angle_ep2 = false;
    bool has_pTmiss = false;

    bool has_Emiss2 = false;
    bool has_Mx2 = false;
    bool has_Mx2_1 = false;
    bool has_Mx2_2 = false;
    bool has_xF = false;
    bool has_Delta_phi = false;
    bool has_theta_gamma_gamma = false;

    bool has_e_p = false;
    bool has_e_theta = false;
    bool has_e_phi = false;

    bool has_p2_p = false;
    bool has_p2_theta = false;
    bool has_p2_phi = false;

    void bind(TTree* tree, const std::string& tree_key) {
        if (!tree) {
            throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Null TTree for key: " + tree_key);
        }
        // endif

        tree->SetBranchStatus("*", 0);

        auto enable = [&](const char* name) {
            if (tree->GetBranch(name)) {
                tree->SetBranchStatus(name, 1);
            }
            // endif
        };

        enable("runnum");
        enable("detector1");
        enable("detector2");

        enable("x");
        enable("Q2");

        enable("t1");
        enable("open_angle_ep2");
        enable("pTmiss");

        enable("Emiss2");
        enable("Mx2");
        enable("Mx2_1");
        enable("Mx2_2");
        enable("xF");
        enable("Delta_phi");
        enable("theta_gamma_gamma");

        enable("e_p");
        enable("e_theta");
        enable("e_phi");

        enable("p2_p");
        enable("p2_theta");
        enable("p2_phi");

        auto bind_int = [&](const char* name, int* addr, bool& has_branch) {
            if (tree->GetBranch(name)) {
                tree->SetBranchAddress(name, addr, nullptr, nullptr, kInt_t, false);
                has_branch = true;
            }
            // endif
        };

        auto bind_double = [&](const char* name, double* addr, bool& has_branch) {
            if (tree->GetBranch(name)) {
                tree->SetBranchAddress(name, addr, nullptr, nullptr, kDouble_t, false);
                has_branch = true;
            }
            // endif
        };

        bind_int("runnum", &runnum, has_runnum);
        bind_int("detector1", &detector1, has_detector1);
        bind_int("detector2", &detector2, has_detector2);

        bind_double("x", &x, has_x);
        bind_double("Q2", &Q2, has_Q2);

        bind_double("t1", &t1, has_t1);
        bind_double("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bind_double("pTmiss", &pTmiss, has_pTmiss);

        bind_double("Emiss2", &Emiss2, has_Emiss2);
        bind_double("Mx2", &Mx2, has_Mx2);
        bind_double("Mx2_1", &Mx2_1, has_Mx2_1);
        bind_double("Mx2_2", &Mx2_2, has_Mx2_2);
        bind_double("xF", &xF, has_xF);
        bind_double("Delta_phi", &Delta_phi, has_Delta_phi);
        bind_double("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);

        bind_double("e_p", &e_p, has_e_p);
        bind_double("e_theta", &e_theta, has_e_theta);
        bind_double("e_phi", &e_phi, has_e_phi);

        bind_double("p2_p", &p2_p, has_p2_p);
        bind_double("p2_theta", &p2_theta, has_p2_theta);
        bind_double("p2_phi", &p2_phi, has_p2_phi);

        std::vector<std::string> missing;

        if (!has_runnum) missing.push_back("runnum");
        if (!has_detector1) missing.push_back("detector1");
        if (!has_detector2) missing.push_back("detector2");

        if (!has_x) missing.push_back("x");
        if (!has_Q2) missing.push_back("Q2");

        if (!has_t1) missing.push_back("t1");
        if (!has_open_angle_ep2) missing.push_back("open_angle_ep2");
        if (!has_pTmiss) missing.push_back("pTmiss");

        if (!has_Emiss2) missing.push_back("Emiss2");
        if (!has_Mx2) missing.push_back("Mx2");
        if (!has_Mx2_1) missing.push_back("Mx2_1");
        if (!has_Mx2_2) missing.push_back("Mx2_2");
        if (!has_xF) missing.push_back("xF");
        if (!has_Delta_phi) missing.push_back("Delta_phi");
        if (!has_theta_gamma_gamma) missing.push_back("theta_gamma_gamma");

        if (!missing.empty()) {
            std::ostringstream ss;
            ss << "[q2_xb_cut_coverage] FATAL: Tree '" << tree_key
               << "' is missing required branch(es): ";

            for (size_t i = 0; i < missing.size(); ++i) {
                if (i > 0) {
                    ss << ", ";
                }
                // endif

                ss << missing[i];
            }
            // endfor

            throw std::runtime_error(ss.str());
        }
        // endif
    }

    std::map<std::string, double> values_map() const {
        std::map<std::string, double> values;

        values["Delta_phi"] = Delta_phi;
        values["theta_gamma_gamma"] = theta_gamma_gamma;
        values["pTmiss"] = pTmiss;
        values["xF"] = xF;
        values["Emiss2"] = Emiss2;
        values["Mx2"] = Mx2;
        values["Mx2_1"] = Mx2_1;
        values["Mx2_2"] = Mx2_2;

        return values;
    }
};

static bool passes_global_for_topology(
    const BranchBinder& b,
    Topology topo,
    const std::string& period_label
) {
    GlobalCutConfig gcfg = default_global_cuts();
    gcfg.enable_topology_filter = true;

    if (topo == Topology::FD_FD) {
        gcfg.required_detector1 = 1;
        gcfg.required_detector2 = 1;
    } else if (topo == Topology::CD_FD) {
        gcfg.required_detector1 = 2;
        gcfg.required_detector2 = 1;
    } else if (topo == Topology::CD_FT) {
        gcfg.required_detector1 = 2;
        gcfg.required_detector2 = 0;
    } else {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Unknown topology in passes_global_for_topology.");
    }
    // endif

    if (gcfg.enable_dvcsgen_ycol_cut) {
        std::vector<std::string> missing;

        if (!b.has_e_p) missing.push_back("e_p");
        if (!b.has_e_theta) missing.push_back("e_theta");
        if (!b.has_e_phi) missing.push_back("e_phi");
        if (!b.has_p2_p) missing.push_back("p2_p");
        if (!b.has_p2_theta) missing.push_back("p2_theta");
        if (!b.has_p2_phi) missing.push_back("p2_phi");

        if (!missing.empty()) {
            std::ostringstream ss;
            ss << "[q2_xb_cut_coverage] FATAL: dvcsgen ycol cut is enabled, but tree for period_label='"
               << period_label << "' is missing required branch(es): ";

            for (size_t i = 0; i < missing.size(); ++i) {
                if (i > 0) {
                    ss << ", ";
                }
                // endif

                ss << missing[i];
            }
            // endfor

            throw std::runtime_error(ss.str());
        }
        // endif

        return passes_global_cuts(
            b.t1,
            b.open_angle_ep2,
            b.pTmiss,
            b.detector1,
            b.detector2,
            period_label,
            b.e_p,
            b.e_theta,
            b.e_phi,
            b.p2_p,
            b.p2_theta,
            b.p2_phi,
            gcfg
        );
    }
    // endif

    return passes_global_cuts(
        b.t1,
        b.open_angle_ep2,
        b.pTmiss,
        b.detector1,
        b.detector2,
        gcfg
    );
}

struct PeriodInput {
    std::string tree_key;
    std::string period_label;
    std::string display_label;
    bool is_inbending = true;
};

struct PointSet {
    std::vector<double> x;
    std::vector<double> q2;
    bool is_inbending = true;
    std::string display_label;
};

struct DrawnObjects {
    std::vector<TH2D*> frames;
    std::vector<TGraph*> graphs;
    std::vector<TLine*> lines;
    std::vector<TLatex*> latex;
};

using CutLookup = std::map<std::string, std::map<Topology, CutDict>>;

static Stats read_stats_from_json(const nlohmann::json& jstats, const std::string& key) {
    if (!jstats.contains("mean") || !jstats.contains("std")) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Cut stats for key '" + key + "' must contain mean and std.");
    }
    // endif

    Stats stats;
    stats.mean = jstats.at("mean").get<double>();
    stats.std = jstats.at("std").get<double>();

    return stats;
}

static CutLookup load_data_cuts(const std::string& json_path) {
    std::ifstream fin(json_path);
    if (!fin) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Could not open combined cuts JSON: " + json_path);
    }
    // endif

    nlohmann::json j;
    fin >> j;

    CutLookup cuts;

    const std::vector<std::string> period_labels = {
        "fa18_inb",
        "fa18_out",
        "sp18_inb",
        "sp18_out",
        "sp19_inb"
    };

    const std::vector<Topology> topologies = {
        Topology::FD_FD,
        Topology::CD_FD,
        Topology::CD_FT
    };

    for (const std::string& period_label : period_labels) {
        const std::string period_code = period_code_dvcs_from_label(period_label);

        for (Topology topo : topologies) {
            const std::string topo_key = topo_to_key(topo);
            const std::string full_key = period_code + "_" + topo_key;

            if (!j.contains(full_key)) {
                throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Missing cuts key in combined cuts JSON: " + full_key);
            }
            // endif

            const nlohmann::json& jobj = j.at(full_key);

            if (!jobj.contains("data")) {
                throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Missing data cuts object for key: " + full_key);
            }
            // endif

            const nlohmann::json& jdata = jobj.at("data");

            CutDict cd;

            for (auto it = jdata.begin(); it != jdata.end(); ++it) {
                const std::string var = it.key();
                cd.data[var] = read_stats_from_json(it.value(), full_key + ".data." + var);
            }
            // endfor

            cuts[period_label][topo] = cd;
        }
        // endfor
    }
    // endfor

    return cuts;
}

static TTree* get_required_tree(
    const std::map<std::string, TTree*>& trees,
    const std::string& key
) {
    const auto it = trees.find(key);

    if (it == trees.end() || !(it->second)) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Required DVCS data tree is missing: " + key);
    }
    // endif

    return it->second;
}

static PointSet collect_points_for_period(
    TTree* tree,
    const PeriodInput& period,
    const CutLookup& cuts
) {
    const auto period_cut_it = cuts.find(period.period_label);
    if (period_cut_it == cuts.end()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: No cut set found for period label: " + period.period_label);
    }
    // endif

    PointSet points;
    points.is_inbending = period.is_inbending;
    points.display_label = period.display_label;

    BranchBinder b;
    b.bind(tree, period.tree_key);

    const Long64_t nentries = tree->GetEntries();

    Long64_t accepted = 0;
    Long64_t in_range = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (is_excluded_run(b.runnum)) {
            continue;
        }
        // endif

        Topology topo = Topology::FD_FD;
        if (!event_topology_from_detectors(b.detector1, b.detector2, topo)) {
            continue;
        }
        // endif

        const auto topo_cut_it = period_cut_it->second.find(topo);
        if (topo_cut_it == period_cut_it->second.end()) {
            throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Missing topology cuts for period label: " + period.period_label);
        }
        // endif

        if (!passes_global_for_topology(b, topo, period.period_label)) {
            continue;
        }
        // endif

        const std::map<std::string, double> values = b.values_map();

        if (!passes_3sigma_cuts(topo_cut_it->second.data, values)) {
            continue;
        }
        // endif

        ++accepted;

        if (b.x < kXBMin || b.x > kXBMax || b.Q2 < kQ2Min || b.Q2 > kQ2Max) {
            continue;
        }
        // endif

        ++in_range;

        points.x.push_back(b.x);
        points.q2.push_back(b.Q2);
    }
    // endfor

    std::cout << "[q2_xb_cut_coverage] " << period.tree_key
              << " entries=" << nentries
              << " accepted_after_global_and_3sigma=" << accepted
              << " plotted_in_range=" << in_range
              << "\n";

    return points;
}

static TGraph* make_graph(
    const std::vector<double>& xs,
    const std::vector<double>& ys,
    int color
) {
    auto* graph = new TGraph(static_cast<int>(xs.size()));

    for (int i = 0; i < static_cast<int>(xs.size()); ++i) {
        graph->SetPoint(i, xs[static_cast<size_t>(i)], ys[static_cast<size_t>(i)]);
    }
    // endfor

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.10);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);

    return graph;
}

static std::string format_axis_value(double x) {
    std::ostringstream ss;

    if (std::abs(x - std::round(x)) < 1.0e-9) {
        ss << static_cast<int>(std::round(x));
        return ss.str();
    }
    // endif

    ss << std::fixed << std::setprecision(2) << x;

    std::string out = ss.str();

    while (!out.empty() && out.back() == '0') {
        out.pop_back();
    }

    if (!out.empty() && out.back() == '.') {
        out.pop_back();
    }
    // endif

    return out;
}

static TH2D* make_frame(const std::string& name) {
    auto* frame = new TH2D(
        name.c_str(),
        "",
        100,
        kXBMin,
        kXBMax,
        100,
        kQ2Min,
        kQ2Max
    );

    frame->SetDirectory(nullptr);

    frame->GetXaxis()->SetTitle("");
    frame->GetYaxis()->SetTitle("");

    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetYaxis()->SetLabelSize(0.0);

    frame->GetXaxis()->SetNdivisions(0);
    frame->GetYaxis()->SetNdivisions(0);

    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetYaxis()->SetTickLength(0.0);

    return frame;
}

static void draw_bin_boundary_lines_on_top(const BinGrid& grid, DrawnObjects& drawn) {
    for (double xb : grid.xb_edges) {
        if (xb < kXBMin || xb > kXBMax) {
            continue;
        }
        // endif

        auto* line = new TLine(xb, kQ2Min, xb, kQ2Max);
        line->SetLineColor(kBinLineColor);
        line->SetLineStyle(1);
        line->SetLineWidth(1);
        line->Draw("same");

        drawn.lines.push_back(line);
    }
    // endfor

    for (double q2 : grid.q2_edges) {
        if (q2 < kQ2Min || q2 > kQ2Max) {
            continue;
        }
        // endif

        auto* line = new TLine(kXBMin, q2, kXBMax, q2);
        line->SetLineColor(kBinLineColor);
        line->SetLineStyle(1);
        line->SetLineWidth(1);
        line->Draw("same");

        drawn.lines.push_back(line);
    }
    // endfor
}

static void draw_custom_bin_edge_labels(const BinGrid& grid, DrawnObjects& drawn) {
    const double left = gPad->GetLeftMargin();
    const double right = gPad->GetRightMargin();
    const double bottom = gPad->GetBottomMargin();
    const double top = gPad->GetTopMargin();

    auto* xlatex = new TLatex();
    xlatex->SetNDC(true);
    xlatex->SetTextAlign(23);
    xlatex->SetTextSize(0.030);

    for (double xb : grid.xb_edges) {
        if (xb < kXBMin || xb > kXBMax) {
            continue;
        }
        // endif

        const double x_ndc = left + (1.0 - left - right) * (xb - kXBMin) / (kXBMax - kXBMin);
        const double y_ndc = bottom * 0.58;

        xlatex->DrawLatex(x_ndc, y_ndc, format_axis_value(xb).c_str());
    }
    // endfor

    drawn.latex.push_back(xlatex);

    auto* ylatex = new TLatex();
    ylatex->SetNDC(true);
    ylatex->SetTextAlign(32);
    ylatex->SetTextSize(0.030);

    for (double q2 : grid.q2_edges) {
        if (q2 < kQ2Min || q2 > kQ2Max) {
            continue;
        }
        // endif

        const double x_ndc = left * 0.74;
        const double y_ndc = bottom + (1.0 - bottom - top) * (q2 - kQ2Min) / (kQ2Max - kQ2Min);

        ylatex->DrawLatex(x_ndc, y_ndc, format_axis_value(q2).c_str());
    }
    // endfor

    drawn.latex.push_back(ylatex);
}

static void draw_axis_titles(DrawnObjects& drawn) {
    const double left = gPad->GetLeftMargin();
    const double bottom = gPad->GetBottomMargin();

    auto* x_title = new TLatex();
    x_title->SetNDC(true);
    x_title->SetTextAlign(22);
    x_title->SetTextSize(0.046);
    x_title->DrawLatex(0.5, bottom * 0.20, "x_{B}");
    drawn.latex.push_back(x_title);

    auto* y_title = new TLatex();
    y_title->SetNDC(true);
    y_title->SetTextAlign(22);
    y_title->SetTextSize(0.046);
    y_title->SetTextAngle(90.0);
    y_title->DrawLatex(left * 0.24, 0.5, "Q^{2} (GeV^{2})");
    drawn.latex.push_back(y_title);
}

static void draw_period_title(const PointSet& points, DrawnObjects& drawn) {
    auto* title = new TLatex();
    title->SetNDC(true);
    title->SetTextAlign(22);
    title->SetTextSize(0.062);
    title->DrawLatex(0.50, 0.965, points.display_label.c_str());
    drawn.latex.push_back(title);

    auto* subtitle = new TLatex();
    subtitle->SetNDC(true);
    subtitle->SetTextAlign(22);
    subtitle->SetTextSize(0.040);

    const char* bending_label = points.is_inbending ? "Inbending DATA" : "Outbending DATA";
    subtitle->DrawLatex(0.50, 0.920, bending_label);
    drawn.latex.push_back(subtitle);
}

static void draw_panel(
    const std::string& frame_name,
    const PointSet& points,
    const BinGrid& grid,
    DrawnObjects& drawn
) {
    TH2D* frame = make_frame(frame_name);
    frame->Draw("axis");
    drawn.frames.push_back(frame);

    const int color = points.is_inbending ? kInbendingColor : kOutbendingColor;
    TGraph* graph = make_graph(points.x, points.q2, color);

    if (graph->GetN() > 0) {
        graph->Draw("P same");
    }
    // endif

    drawn.graphs.push_back(graph);

    draw_bin_boundary_lines_on_top(grid, drawn);

    frame->Draw("axis same");

    draw_custom_bin_edge_labels(grid, drawn);
    draw_axis_titles(drawn);
    draw_period_title(points, drawn);
}

static void draw_empty_panel(DrawnObjects& drawn) {
    gPad->SetFillColor(kWhite);
    gPad->SetFrameFillColor(kWhite);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.16);
    gPad->SetBottomMargin(0.22);

    auto* latex = new TLatex();
    latex->SetNDC(true);
    latex->SetTextAlign(22);
    latex->SetTextSize(0.055);
    latex->DrawLatex(0.50, 0.50, "No outbending Sp19 period");
    drawn.latex.push_back(latex);
}

static void setup_pad_margins() {
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.16);
    gPad->SetBottomMargin(0.22);
    gPad->SetGrid(0, 0);
    gPad->SetTicks(1, 1);
}

static bool file_exists_nonempty(const std::string& path) {
    std::ifstream fin(path, std::ios::binary | std::ios::ate);
    if (!fin) {
        return false;
    }
    // endif

    return fin.tellg() > 0;
}

static std::string shell_quote(const std::string& s) {
    std::string out = "'";

    for (char c : s) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out += c;
        }
        // endif
    }
    // endfor

    out += "'";
    return out;
}

static void write_png_from_pdf_or_canvas(
    TCanvas* canvas,
    const std::string& out_pdf,
    const std::string& out_png
) {
    canvas->cd();
    canvas->Modified();
    canvas->Update();
    gSystem->ProcessEvents();

    const std::string gs_cmd =
        "gs -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r200 "
        "-dTextAlphaBits=4 -dGraphicsAlphaBits=4 "
        "-sOutputFile=" + shell_quote(out_png) + " " + shell_quote(out_pdf);

    const int gs_status = std::system(gs_cmd.c_str());

    if (gs_status == 0 && file_exists_nonempty(out_png)) {
        std::cout << "[q2_xb_cut_coverage] PNG rasterized from PDF using Ghostscript.\n";
        return;
    }
    // endif

    std::cerr << "[q2_xb_cut_coverage] WARNING: Ghostscript PNG rasterization failed. "
              << "Falling back to ROOT canvas PNG output.\n";

    canvas->Print(out_png.c_str());
}

} // namespace

bool plot_q2_xb_cut_coverage(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& pass2_csv_path,
    const std::string& combined_cuts_json_path,
    const std::string& output_dir
) {
    try {
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);
        gStyle->SetTitleBorderSize(0);
        gStyle->SetTitleFillColor(0);
        gStyle->SetTitleFontSize(0.0);

        if (!output_dir.empty()) {
            gSystem->mkdir(output_dir.c_str(), true);
        }
        // endif

        const BinGrid grid = load_bin_grid_from_csv(pass2_csv_path);
        const CutLookup cuts = load_data_cuts(combined_cuts_json_path);

        const std::vector<PeriodInput> periods = {
            {"DVCS_Fa18_inb", "fa18_inb", "Fa18 Inb", true},
            {"DVCS_Fa18_out", "fa18_out", "Fa18 Out", false},
            {"DVCS_Sp18_inb", "sp18_inb", "Sp18 Inb", true},
            {"DVCS_Sp18_out", "sp18_out", "Sp18 Out", false},
            {"DVCS_Sp19_inb", "sp19_inb", "Sp19 Inb", true}
        };

        std::vector<PointSet> all_points;
        all_points.reserve(periods.size());

        for (const PeriodInput& period : periods) {
            TTree* tree = get_required_tree(dvcsDataTrees, period.tree_key);
            all_points.push_back(collect_points_for_period(tree, period, cuts));
        }
        // endfor

        DrawnObjects drawn;

        auto* canvas = new TCanvas(
            "c_q2_xb_cut_coverage",
            "Q2 vs xB coverage after global and 3sigma cuts",
            2400,
            1700
        );

        canvas->Divide(2, 3, 0.025, 0.035);

        for (int i = 0; i < static_cast<int>(all_points.size()); ++i) {
            canvas->cd(i + 1);
            setup_pad_margins();

            std::ostringstream frame_name;
            frame_name << "h_q2_xb_frame_" << i;

            draw_panel(
                frame_name.str(),
                all_points[static_cast<size_t>(i)],
                grid,
                drawn
            );

            gPad->Modified();
            gPad->Update();
        }
        // endfor

        canvas->cd(6);
        draw_empty_panel(drawn);
        gPad->Modified();
        gPad->Update();

        canvas->cd();
        canvas->Modified();
        canvas->Update();
        gSystem->ProcessEvents();

        const std::string out_pdf = output_dir + "/q2_xb_cut_coverage_by_period.pdf";
        const std::string out_png = output_dir + "/q2_xb_cut_coverage_by_period.png";

        canvas->SaveAs(out_pdf.c_str());
        write_png_from_pdf_or_canvas(canvas, out_pdf, out_png);

        std::cout << "[q2_xb_cut_coverage] Wrote " << out_pdf << "\n";
        std::cout << "[q2_xb_cut_coverage] Wrote " << out_png << "\n";

        delete canvas;

        for (TH2D* obj : drawn.frames) {
            delete obj;
        }
        // endfor

        for (TGraph* obj : drawn.graphs) {
            delete obj;
        }
        // endfor

        for (TLine* obj : drawn.lines) {
            delete obj;
        }
        // endfor

        for (TLatex* obj : drawn.latex) {
            delete obj;
        }
        // endfor

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return false;
    }
    // endtry
}