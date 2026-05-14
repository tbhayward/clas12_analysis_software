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
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TDataType.h>

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

    if (detector1 == 2 && detector2 == 1) {
        topo = Topology::CD_FD;
        return true;
    }

    if (detector1 == 2 && detector2 == 0) {
        topo = Topology::CD_FT;
        return true;
    }

    return false;
}

static std::string period_code_dvcs_from_label(const std::string& run_tag) {
    std::string nice = run_tag;

    if (!nice.empty()) {
        nice[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(nice[0])));
    }

    for (size_t i = 0; i + 1 < nice.size(); ++i) {
        if (nice[i] == '_') {
            nice[i + 1] = static_cast<char>(std::toupper(static_cast<unsigned char>(nice[i + 1])));
        }
    }

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

        if (!within_3sigma(it->second, stats)) {
            return false;
        }
    }

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
        } else if (c == ',' && !in_quotes) {
            fields.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    fields.push_back(cur);
    return fields;
}

static double parse_required_double(const std::string& text, const std::string& context) {
    if (text.empty()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Empty numeric field while parsing " + context + ".");
    }

    char* endptr = nullptr;
    const double val = std::strtod(text.c_str(), &endptr);

    if (endptr == text.c_str()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Could not parse numeric field '" + text + "' while parsing " + context + ".");
    }

    return val;
}

static int find_required_column(const std::vector<std::string>& header, const std::string& name) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) {
            return static_cast<int>(i);
        }
    }

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
    }

    edges.push_back(rounded);
}

static BinGrid load_bin_grid_from_csv(const std::string& csv_path) {
    std::ifstream fin(csv_path);
    if (!fin) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Could not open CSV: " + csv_path);
    }

    std::string line;
    if (!std::getline(fin, line)) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: CSV is empty: " + csv_path);
    }

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

        const std::vector<std::string> row = parse_csv_line(line);

        const int needed = std::max(std::max(ixmin, ixmax), std::max(iqmin, iqmax));
        if (static_cast<int>(row.size()) <= needed) {
            std::ostringstream ss;
            ss << "[q2_xb_cut_coverage] FATAL: CSV row " << row_index
               << " has too few columns.";
            throw std::runtime_error(ss.str());
        }

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

    if (grid.q2_edges.empty()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: No Q2 bin edges found in CSV.");
    }

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

        tree->SetBranchStatus("*", 0);

        auto enable = [&](const char* name) {
            if (tree->GetBranch(name)) {
                tree->SetBranchStatus(name, 1);
            }
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
        };

        auto bind_double = [&](const char* name, double* addr, bool& has_branch) {
            if (tree->GetBranch(name)) {
                tree->SetBranchAddress(name, addr, nullptr, nullptr, kDouble_t, false);
                has_branch = true;
            }
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

                ss << missing[i];
            }

            throw std::runtime_error(ss.str());
        }
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

                ss << missing[i];
            }

            throw std::runtime_error(ss.str());
        }

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
    bool is_inbending = true;
};

struct PointSet {
    std::vector<double> inb_x;
    std::vector<double> inb_q2;

    std::vector<double> out_x;
    std::vector<double> out_q2;
};

using CutLookup = std::map<std::string, std::map<Topology, CutDict>>;

static Stats read_stats_from_json(const nlohmann::json& jstats, const std::string& key) {
    if (!jstats.contains("mean") || !jstats.contains("std")) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Cut stats for key '" + key + "' must contain mean and std.");
    }

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

            const nlohmann::json& jobj = j.at(full_key);

            if (!jobj.contains("data")) {
                throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Missing data cuts object for key: " + full_key);
            }

            const nlohmann::json& jdata = jobj.at("data");

            CutDict cd;

            for (auto it = jdata.begin(); it != jdata.end(); ++it) {
                const std::string var = it.key();
                cd.data[var] = read_stats_from_json(it.value(), full_key + ".data." + var);
            }

            cuts[period_label][topo] = cd;
        }
    }

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

    return it->second;
}

static void collect_points_for_period(
    TTree* tree,
    const std::string& tree_key,
    const std::string& period_label,
    bool is_inbending,
    const CutLookup& cuts,
    PointSet& points
) {
    const auto period_cut_it = cuts.find(period_label);
    if (period_cut_it == cuts.end()) {
        throw std::runtime_error("[q2_xb_cut_coverage] FATAL: No cut set found for period label: " + period_label);
    }

    BranchBinder b;
    b.bind(tree, tree_key);

    const Long64_t nentries = tree->GetEntries();

    Long64_t accepted = 0;
    Long64_t in_range = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (is_excluded_run(b.runnum)) {
            continue;
        }

        Topology topo = Topology::FD_FD;
        if (!event_topology_from_detectors(b.detector1, b.detector2, topo)) {
            continue;
        }

        const auto topo_cut_it = period_cut_it->second.find(topo);
        if (topo_cut_it == period_cut_it->second.end()) {
            throw std::runtime_error("[q2_xb_cut_coverage] FATAL: Missing topology cuts for period label: " + period_label);
        }

        if (!passes_global_for_topology(b, topo, period_label)) {
            continue;
        }

        const std::map<std::string, double> values = b.values_map();

        if (!passes_3sigma_cuts(topo_cut_it->second.data, values)) {
            continue;
        }

        ++accepted;

        if (b.x < kXBMin || b.x > kXBMax || b.Q2 < kQ2Min || b.Q2 > kQ2Max) {
            continue;
        }

        ++in_range;

        if (is_inbending) {
            points.inb_x.push_back(b.x);
            points.inb_q2.push_back(b.Q2);
        } else {
            points.out_x.push_back(b.x);
            points.out_q2.push_back(b.Q2);
        }
    }

    std::cout << "[q2_xb_cut_coverage] " << tree_key
              << " entries=" << nentries
              << " accepted_after_global_and_3sigma=" << accepted
              << " plotted_in_range=" << in_range
              << "\n";
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

    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.18);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);

    return graph;
}

static void draw_grid_lines(const BinGrid& grid) {
    for (double xb : grid.xb_edges) {
        if (xb < kXBMin || xb > kXBMax) {
            continue;
        }

        auto* line = new TLine(xb, kQ2Min, xb, kQ2Max);
        line->SetLineColor(kGray + 1);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("same");
    }

    for (double q2 : grid.q2_edges) {
        if (q2 < kQ2Min || q2 > kQ2Max) {
            continue;
        }

        auto* line = new TLine(kXBMin, q2, kXBMax, q2);
        line->SetLineColor(kGray + 1);
        line->SetLineStyle(2);
        line->SetLineWidth(1);
        line->Draw("same");
    }
}

static TH2D* make_frame(const std::string& name, const std::string& title) {
    auto* frame = new TH2D(
        name.c_str(),
        title.c_str(),
        100,
        kXBMin,
        kXBMax,
        100,
        kQ2Min,
        kQ2Max
    );

    frame->SetDirectory(nullptr);
    frame->GetXaxis()->SetTitle("x_{B}");
    frame->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleOffset(0.95);
    frame->GetYaxis()->SetTitleOffset(1.20);

    return frame;
}

static void draw_panel(
    const std::string& frame_name,
    const std::string& title,
    const PointSet& points,
    const BinGrid& grid,
    bool draw_outbending
) {
    TH2D* frame = make_frame(frame_name, title);
    frame->Draw("axis");

    draw_grid_lines(grid);

    TGraph* g_inb = make_graph(points.inb_x, points.inb_q2, kRed + 1);
    TGraph* g_out = make_graph(points.out_x, points.out_q2, kBlue + 1);

    if (g_inb->GetN() > 0) {
        g_inb->Draw("P same");
    }

    if (draw_outbending && g_out->GetN() > 0) {
        g_out->Draw("P same");
    }

    frame->Draw("axis same");

    auto* legend = new TLegend(0.16, 0.78, 0.52, 0.91);
    legend->SetFillColor(kWhite);
    legend->SetFillStyle(1001);
    legend->SetBorderSize(1);
    legend->SetTextSize(0.040);

    legend->AddEntry(g_inb, "Inbending DATA", "p");

    if (draw_outbending) {
        legend->AddEntry(g_out, "Outbending DATA", "p");
    }

    legend->Draw("same");
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

        if (!output_dir.empty()) {
            gSystem->mkdir(output_dir.c_str(), true);
        }

        const BinGrid grid = load_bin_grid_from_csv(pass2_csv_path);
        const CutLookup cuts = load_data_cuts(combined_cuts_json_path);

        PointSet fa_sp_18_points;
        PointSet sp19_points;

        const std::vector<PeriodInput> left_periods = {
            {"DVCS_Fa18_inb",  "fa18_inb",  true},
            {"DVCS_Sp18_inb",  "sp18_inb",  true},
            {"DVCS_Fa18_out",  "fa18_out",  false},
            {"DVCS_Sp18_out",  "sp18_out",  false}
        };

        const std::vector<PeriodInput> right_periods = {
            {"DVCS_Sp19_inb",  "sp19_inb",  true}
        };

        for (const PeriodInput& period : left_periods) {
            TTree* tree = get_required_tree(dvcsDataTrees, period.tree_key);

            collect_points_for_period(
                tree,
                period.tree_key,
                period.period_label,
                period.is_inbending,
                cuts,
                fa_sp_18_points
            );
        }

        for (const PeriodInput& period : right_periods) {
            TTree* tree = get_required_tree(dvcsDataTrees, period.tree_key);

            collect_points_for_period(
                tree,
                period.tree_key,
                period.period_label,
                period.is_inbending,
                cuts,
                sp19_points
            );
        }

        auto* canvas = new TCanvas(
            "c_q2_xb_cut_coverage",
            "Q2 vs xB coverage after global and 3sigma cuts",
            1600,
            750
        );

        canvas->Divide(2, 1);

        canvas->cd(1);
        gPad->SetLeftMargin(0.13);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.13);
        gPad->SetGrid(0, 0);

        draw_panel(
            "h_q2_xb_frame_fa18_sp18",
            "Fa18 + Sp18 DATA after global and 3#sigma DVCS cuts",
            fa_sp_18_points,
            grid,
            true
        );

        canvas->cd(2);
        gPad->SetLeftMargin(0.13);
        gPad->SetRightMargin(0.04);
        gPad->SetTopMargin(0.10);
        gPad->SetBottomMargin(0.13);
        gPad->SetGrid(0, 0);

        draw_panel(
            "h_q2_xb_frame_sp19",
            "Sp19 Inb DATA after global and 3#sigma DVCS cuts",
            sp19_points,
            grid,
            false
        );

        canvas->cd();

        const std::string out_png = output_dir + "/q2_xb_cut_coverage.png";
        const std::string out_pdf = output_dir + "/q2_xb_cut_coverage.pdf";

        canvas->SaveAs(out_png.c_str());
        canvas->SaveAs(out_pdf.c_str());

        std::cout << "[q2_xb_cut_coverage] Wrote " << out_png << "\n";
        std::cout << "[q2_xb_cut_coverage] Wrote " << out_pdf << "\n";

        delete canvas;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return false;
    }
}