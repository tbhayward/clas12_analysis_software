// total_counts.cpp
// -----------------------------------------------------------------------------
// Fill raw data counts and MC counts into the pass-2 CSV.
//
// DATA columns filled:
//   raw yield, ep->epg,   <topo label>, exp, <period>, unpol
//   raw yield, ep->eppi0, <topo label>, exp, <period>, unpol
//   raw yield, ep->epg,   <topo label>, exp, <period>, pos/neg for helicity-qualified periods only
//   raw yield, ep->eppi0, <topo label>, exp, <period>, pos/neg for helicity-qualified periods only
//
// MC columns filled:
//   generated yield, ep->epg, mc, <period>
//   reconstructed yield, ep->epg, mc, <period>
//   reconstructed yield, ep->epg, <topo label>, mc, <period>
//
//   generated yield, ep->eppi0, mc, <period>
//   reconstructed yield, ep->eppi0, mc, <period>
//   reconstructed yield, ep->eppi0, <topo label>, mc, <period>
//
//   reconstructed yield, ep->eppi0->epg, mc, <period>
//   reconstructed yield, ep->eppi0->epg, <topo label>, mc, <period>
//
// In addition to the untouched unit-weight columns above, the nominal workflow
// also fills:
//   normalized raw yield, ...                         (DATA)
//   reconstructed current corrected yield, ...       (reconstructed MC)
//
// Those are accumulated event-by-event with the regional current-response model
// produced upstream by current_dependence.cpp. Raw totals remain unit-weighted.
//
// Speed/stability notes:
//   - Parallelized over independent ROOT trees/work items.
//   - Hard cap of seven workers.
//   - ROOT branch binding is mutex-protected.
//   - Each worker writes to local storage only; results are merged after each tree.
//   - Fast row lookup avoids scanning every CSV row for every event.
// -----------------------------------------------------------------------------

#include "total_counts.h"

#include "periods.h"
#include "global_cuts.h"

// ROOT
#include <TTree.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TLine.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH1F.h>
#include <TAxis.h>

// JSON
#include <nlohmann/json.hpp>

// C++ stdlib
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

static constexpr double PI      = 3.14159265358979323846;
static constexpr double RAD2DEG = 180.0 / PI;

enum class Channel {
    DVCS,
    EPPI0,
    EPPI0_BKG_AS_DVCS
};

enum class SampleKind {
    DATA,
    MC_GEN,
    MC_REC
};

struct ChannelConfig {
    Channel channel = Channel::DVCS;
    std::string csv_channel;
    std::string cut_prefix;
    std::string plot_subdir;
    std::string plot_file_token;
    std::string title_label;
    bool uses_dvcs_topology_cuts = true;
};

static ChannelConfig dvcs_config() {
    ChannelConfig cfg;
    cfg.channel = Channel::DVCS;
    cfg.csv_channel = "ep->epg";
    cfg.cut_prefix = "DVCS";
    cfg.plot_subdir = "DVCS";
    cfg.plot_file_token = "";
    cfg.title_label = "ep #rightarrow ep#gamma";
    cfg.uses_dvcs_topology_cuts = true;
    return cfg;
}

static ChannelConfig eppi0_config() {
    ChannelConfig cfg;
    cfg.channel = Channel::EPPI0;
    cfg.csv_channel = "ep->eppi0";
    cfg.cut_prefix = "eppi0";
    cfg.plot_subdir = "eppi0";
    cfg.plot_file_token = "eppi0";
    cfg.title_label = "ep #rightarrow ep#pi_{0}";
    cfg.uses_dvcs_topology_cuts = false;
    return cfg;
}

static ChannelConfig eppi0_bkg_as_dvcs_config() {
    ChannelConfig cfg;
    cfg.channel = Channel::EPPI0_BKG_AS_DVCS;
    cfg.csv_channel = "ep->eppi0->epg";
    cfg.cut_prefix = "DVCS";
    cfg.plot_subdir = "eppi0_bkg_as_dvcs";
    cfg.plot_file_token = "eppi0_bkg_as_dvcs";
    cfg.title_label = "ep#pi_{0} #rightarrow ep#gamma selection";
    cfg.uses_dvcs_topology_cuts = true;
    return cfg;
}

struct WorkConfig {
    ChannelConfig channel_cfg;
    SampleKind sample_kind = SampleKind::DATA;
    bool write_data_raw_columns = false;
    bool write_mc_generated_columns = false;
    bool write_mc_reconstructed_columns = false;
    bool make_plots = false;
};


static inline void fatal(const std::string& msg);

static constexpr int kCurrentRegionCount = 7;

static const std::array<std::string, kCurrentRegionCount>& current_region_names() {
    static const std::array<std::string, kCurrentRegionCount> names = {
        "FT", "S1", "S2", "S3", "S4", "S5", "S6"
    };
    return names;
}

static int current_region_index(int detector2, double p2_phi_rad, bool has_p2_phi) {
    if (detector2 == 0) return 0;
    if (detector2 != 1 || !has_p2_phi || !std::isfinite(p2_phi_rad)) return -1;
    double phi = std::fmod(p2_phi_rad * RAD2DEG, 360.0);
    if (phi < 0.0) phi += 360.0;
    int sector = 0;
    if (phi >= 330.0 || phi < 30.0) sector = 1;
    else if (phi < 90.0) sector = 2;
    else if (phi < 150.0) sector = 3;
    else if (phi < 210.0) sector = 4;
    else if (phi < 270.0) sector = 5;
    else if (phi < 330.0) sector = 6;
    return (sector >= 1 && sector <= 6) ? sector : -1;
}

struct CurrentResponseEntry {
    bool valid = false;
    double parameter = std::numeric_limits<double>::quiet_NaN();
    double parameter_stat = 0.0;

    // Optional centered polar-angle term for DATA:
    //   s(theta) = parameter + angular_gradient*(theta-angular_center).
    // The fitted centering makes the regional baseline and common gradient
    // effectively orthogonal, so their calibration uncertainties can be
    // propagated as separate correlated nuisances.
    bool has_angular_model = false;
    std::string angular_variable;
    double angular_center = 0.0;
    double angular_gradient = 0.0;
    double angular_gradient_stat = 0.0;
};

struct CurrentResponseModel {
    bool valid = false;
    std::map<std::string, std::map<std::string, std::array<CurrentResponseEntry, kCurrentRegionCount>>> data;
    std::map<std::string, std::map<std::string, std::array<CurrentResponseEntry, kCurrentRegionCount>>> mc;
    std::map<std::string, std::unordered_map<int, int>> run_current_nA;
    std::map<std::string, std::unordered_set<int>> excluded_data_runs;
};

static CurrentResponseModel load_current_response_model(const std::string& path) {
    CurrentResponseModel model;
    std::ifstream fin(path);
    if (!fin.is_open()) {
        fatal("[total_counts] FATAL: cannot open current-response model JSON: " + path);
    }
    nlohmann::json j;
    fin >> j;

    auto load_block = [&](const char* sample,
                          std::map<std::string, std::map<std::string, std::array<CurrentResponseEntry, kCurrentRegionCount>>>& dst,
                          const char* value_key,
                          const char* stat_key) {
        if (!j.contains(sample) || !j[sample].is_object()) return;
        for (auto ch = j[sample].begin(); ch != j[sample].end(); ++ch) {
            for (auto per = ch.value().begin(); per != ch.value().end(); ++per) {
                std::array<CurrentResponseEntry, kCurrentRegionCount> arr{};
                for (int ir = 0; ir < kCurrentRegionCount; ++ir) {
                    const std::string& rn = current_region_names()[ir];
                    if (!per.value().contains(rn)) continue;
                    const auto& e = per.value()[rn];
                    if (!e.contains(value_key)) continue;
                    const double v = e[value_key].get<double>();
                    const double se = e.contains(stat_key) ? e[stat_key].get<double>() : 0.0;
                    if (std::isfinite(v) && std::isfinite(se) && se >= 0.0) {
                        arr[ir].valid = true;
                        arr[ir].parameter = v;
                        arr[ir].parameter_stat = se;

                        if (std::string(sample) == "data" &&
                            e.contains("angular_model") &&
                            e["angular_model"].is_object()) {
                            const auto& a = e["angular_model"];
                            const std::string variable =
                                a.value("variable", std::string());
                            const double center =
                                a.value("center", std::numeric_limits<double>::quiet_NaN());
                            const double gradient =
                                a.value("gradient_per_nA_per_deg",
                                        std::numeric_limits<double>::quiet_NaN());
                            const double gradient_stat =
                                a.value("gradient_stat",
                                        std::numeric_limits<double>::quiet_NaN());
                            if (variable == "e_theta" &&
                                std::isfinite(center) &&
                                std::isfinite(gradient) &&
                                std::isfinite(gradient_stat) &&
                                gradient_stat >= 0.0) {
                                arr[ir].has_angular_model = true;
                                arr[ir].angular_variable = variable;
                                arr[ir].angular_center = center;
                                arr[ir].angular_gradient = gradient;
                                arr[ir].angular_gradient_stat = gradient_stat;
                            } else {
                                fatal("[total_counts] FATAL: malformed DATA angular current-response model in " + path);
                            }
                        }
                    }
                }
                dst[ch.key()][per.key()] = arr;
            }
        }
    };

    load_block("data", model.data, "relative_slope_per_nA", "relative_slope_stat");
    load_block("mc", model.mc, "reference_factor", "reference_factor_stat");

    if (j.contains("run_current_nA") && j["run_current_nA"].is_object()) {
        for (auto per = j["run_current_nA"].begin(); per != j["run_current_nA"].end(); ++per) {
            for (auto rr = per.value().begin(); rr != per.value().end(); ++rr) {
                model.run_current_nA[per.key()][std::stoi(rr.key())] = rr.value().get<int>();
            }
        }
    }

    if (j.contains("excluded_data_runs") && j["excluded_data_runs"].is_object()) {
        for (auto per = j["excluded_data_runs"].begin(); per != j["excluded_data_runs"].end(); ++per) {
            if (!per.value().is_object()) continue;
            for (auto rr = per.value().begin(); rr != per.value().end(); ++rr) {
                model.excluded_data_runs[per.key()].insert(std::stoi(rr.key()));
            }
        }
    }

    model.valid = !model.data.empty() && !model.mc.empty() && !model.run_current_nA.empty();
    if (!model.valid) {
        fatal("[total_counts] FATAL: current-response model is incomplete: " + path);
    }
    return model;
}

static void validate_final_current_response_prescription(
    const CurrentResponseModel& model,
    bool require_sp18_out_epg_e_theta_current_model) {

    if (!require_sp18_out_epg_e_theta_current_model) return;

    const auto ic = model.data.find("ep->epg");
    if (ic == model.data.end()) {
        fatal("[total_counts] FATAL: final current-response prescription requires ep->epg DATA calibration.");
    }
    const auto ip = ic->second.find("Sp18 Out");
    if (ip == ic->second.end()) {
        fatal("[total_counts] FATAL: final current-response prescription requires Sp18 Out ep->epg DATA calibration.");
    }

    for (int ir = 0; ir < kCurrentRegionCount; ++ir) {
        const CurrentResponseEntry& e = ip->second[ir];
        if (!e.valid || !e.has_angular_model ||
            e.angular_variable != "e_theta" ||
            !std::isfinite(e.angular_center) ||
            !std::isfinite(e.angular_gradient) ||
            !std::isfinite(e.angular_gradient_stat)) {
            fatal("[total_counts] FATAL: final current-response prescription requires a valid centered linear e_theta model for Sp18 Out ep->epg region '" +
                  current_region_names()[ir] + "'. Regenerate current_response_model.json with use_sp18_out_e_theta_response_model=true.");
        }
    }

    std::cout << "[total_counts] Final current-response prescription validated: "
              << "Sp18 Out ep->epg uses regional + centered linear e_theta DATA response; "
              << "all other angular extensions remain disabled." << std::endl;
}

struct EventCurrentWeight {
    double weight = 1.0;
    double derivative = 0.0;  // derivative wrt the regional baseline parameter
    double parameter_stat = 0.0;
    double angular_derivative = 0.0; // derivative wrt common angular gradient
    double angular_parameter_stat = 0.0;
    int region = -1;
    bool skip_corrected = false;
};

static const CurrentResponseEntry* find_current_response_entry(
    const std::map<std::string, std::map<std::string, std::array<CurrentResponseEntry, kCurrentRegionCount>>>& table,
    const std::string& channel,
    const std::string& period,
    int region) {
    auto ic = table.find(channel);
    if (ic == table.end()) return nullptr;
    auto ip = ic->second.find(period);
    if (ip == ic->second.end() || region < 0 || region >= kCurrentRegionCount) return nullptr;
    return &ip->second[region];
}

static inline bool env_flag(const char* name) {
    return (std::getenv(name) != nullptr);
}

static inline void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static inline std::string to_lower_ascii(std::string s) {
    for (char& c : s) {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

static inline double wrap_phi_deg(double phi_deg) {
    double p = std::fmod(phi_deg, 360.0);

    if (p < 0.0) {
        p += 360.0;
    }

    if (p >= 360.0) {
        p = std::nextafter(360.0, 0.0);
    }

    return p;
}


static inline double delta_phi_rad_from_two_phi(double phi_a, double phi_b) {
    double d = std::fmod(phi_a - phi_b, 2.0 * PI);

    if (d <= -PI) {
        d += 2.0 * PI;
    }

    if (d > PI) {
        d -= 2.0 * PI;
    }

    return std::fabs(d);
}

static inline bool in_range(double v, double a, double b) {
    return (v >= a) && (v < b);
}

static inline bool row_accepts_phi(double phi_deg, double pmin_deg, double pmax_deg) {
    if (pmax_deg > pmin_deg) {
        return in_range(phi_deg, pmin_deg, pmax_deg);
    }

    return (phi_deg >= pmin_deg) || (phi_deg < pmax_deg);
}

enum class TopologyIndex : int {
    FD_FD = 0,
    CD_FD = 1,
    CD_FT = 2,
    INVALID = -1
};

static inline TopologyIndex topology_index(int det1, int det2) {
    if (det1 == 1 && det2 == 1) {
        return TopologyIndex::FD_FD;
    }

    if (det1 == 2 && det2 == 1) {
        return TopologyIndex::CD_FD;
    }

    if (det1 == 2 && det2 == 0) {
        return TopologyIndex::CD_FT;
    }

    return TopologyIndex::INVALID;
}

static inline const std::string& topology_name(TopologyIndex topo) {
    static const std::string kFdFd = "FD_FD";
    static const std::string kCdFd = "CD_FD";
    static const std::string kCdFt = "CD_FT";
    static const std::string kInvalid;

    switch (topo) {
        case TopologyIndex::FD_FD: return kFdFd;
        case TopologyIndex::CD_FD: return kCdFd;
        case TopologyIndex::CD_FT: return kCdFt;
        default:                   return kInvalid;
    }
}

static inline std::string topo_label_for_csv(const std::string& topoDir) {
    if (topoDir == "FD_FD") {
        return "(FD, FD)";
    }

    if (topoDir == "CD_FD") {
        return "(CD, FD)";
    }

    if (topoDir == "CD_FT") {
        return "(CD, FT)";
    }

    return "";
}

static inline std::string canonical_period_dir(const std::string& label) {
    if (label == "Fa18 Inb") {
        return "Fa18_Inb";
    }

    if (label == "Fa18 Out") {
        return "Fa18_Out";
    }

    if (label == "Fa18 Inb Supp") {
        return "Fa18_Inb_Supp";
    }

    if (label == "Sp18 Inb") {
        return "Sp18_Inb";
    }

    if (label == "Sp18 Out") {
        return "Sp18_Out";
    }

    if (label == "Sp19 Inb") {
        return "Sp19_Inb";
    }

    std::ostringstream ss;
    ss << "[total_counts] FATAL: unknown period label for canonical_period_dir: '"
       << label << "'";
    fatal(ss.str());

    return "";
}

static inline std::string canonical_group_dir(const std::string& label) {
    if (label == "Fa18") {
        return "Fa18";
    }

    if (label == "Sp18") {
        return "Sp18";
    }

    if (label == "10.6 GeV") {
        return "10.6 GeV";
    }

    std::ostringstream ss;
    ss << "[total_counts] FATAL: unknown group label for canonical_group_dir: '"
       << label << "'";
    fatal(ss.str());

    return "";
}

static inline bool is_combined_group_label(const std::string& label) {
    return (label == "Fa18" || label == "Sp18" || label == "10.6 GeV");
}

static inline bool is_supplemental_label(const std::string& label) {
    return (label == "Fa18 Inb Supp");
}

static inline bool should_skip_csv_for_label(const std::string& label) {
    if (is_supplemental_label(label)) {
        return true;
    }

    if (is_combined_group_label(label)) {
        return true;
    }

    return false;
}

static inline bool has_helicity_resolved_data_columns(const std::string& label) {
    if (label == "Sp18 Inb") {
        return false;
    }

    if (label == "Sp18 Out") {
        return false;
    }

    return true;
}

static inline bool ends_with_path_component(const std::string& path,
                                                const std::string& component) {
    if (path == component) {
        return true;
    }

    if (path.size() <= component.size()) {
        return false;
    }

    const size_t start = path.size() - component.size();
    return path.compare(start, component.size(), component) == 0 &&
           (start == 0 || path[start - 1] == '/');
}

static inline std::string normalize_total_counts_root(const std::string& out_root_dir) {
    if (out_root_dir.empty()) {
        return "output/total_counts_plots";
    }

    if (ends_with_path_component(out_root_dir, "total_counts_plots")) {
        return out_root_dir;
    }

    return out_root_dir + "/total_counts_plots";
}

static inline std::string out_root_for_label(const ChannelConfig& channel_cfg,
                                             const std::string& label,
                                             const std::string& out_root_dir) {
    std::string base = normalize_total_counts_root(out_root_dir);

    if (channel_cfg.plot_subdir.empty()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: empty plot_subdir for channel "
           << channel_cfg.csv_channel;
        fatal(ss.str());
    }

    base += "/" + channel_cfg.plot_subdir;

    if (is_combined_group_label(label)) {
        return base + "/" + canonical_group_dir(label);
    }

    return base + "/" + canonical_period_dir(label);
}

static inline void mkdir_p(const std::string& path) {
    if (path.empty()) {
        return;
    }

    gSystem->mkdir(path.c_str(), true);
}

// -----------------------------------------------------------------------------
// Period key parsing
// -----------------------------------------------------------------------------

struct PeriodTags {
    std::string tree_key;
    std::string period_display;
    std::string period_label;
    std::string period_code;
};

static inline PeriodTags parse_period_tags_from_tree_key(const std::string& tree_key) {
    const std::string s = to_lower_ascii(tree_key);

    PeriodTags t;
    t.tree_key = tree_key;

    auto has = [&](const char* sub) {
        return (s.find(sub) != std::string::npos);
    };

    if (has("fa18") && has("inb") && has("supp")) {
        t.period_display = "Fa18 Inb Supp";
        t.period_label   = "fa18_inb";
        t.period_code    = "Fa18_Inb_Supp";
        return t;
    }

    if (has("fa18") && has("inb")) {
        t.period_display = "Fa18 Inb";
        t.period_label   = "fa18_inb";
        t.period_code    = "Fa18_Inb";
        return t;
    }

    if (has("fa18") && has("out")) {
        t.period_display = "Fa18 Out";
        t.period_label   = "fa18_out";
        t.period_code    = "Fa18_Out";
        return t;
    }

    if (has("sp18") && has("inb")) {
        t.period_display = "Sp18 Inb";
        t.period_label   = "sp18_inb";
        t.period_code    = "Sp18_Inb";
        return t;
    }

    if (has("sp18") && has("out")) {
        t.period_display = "Sp18 Out";
        t.period_label   = "sp18_out";
        t.period_code    = "Sp18_Out";
        return t;
    }

    if (has("sp19") && has("inb")) {
        t.period_display = "Sp19 Inb";
        t.period_label   = "sp19_inb";
        t.period_code    = "Sp19_Inb";
        return t;
    }

    std::ostringstream ss;
    ss << "[total_counts] FATAL: cannot map tree_key '" << tree_key
       << "' to PeriodTags.";
    fatal(ss.str());

    return t;
}

static inline std::string combined_cuts_key(const ChannelConfig& channel_cfg,
                                            const PeriodTags& tags,
                                            const std::string& topoDir) {
    return channel_cfg.cut_prefix + "_" + tags.period_code + "_" + topoDir;
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

    bool inq = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            if (inq && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                inq = !inq;
            }
        } else if (c == ',' && !inq) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    out.push_back(cur);
    return out;
}

static bool load_csv(const std::string& path, CSV& csv) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cannot open CSV: " << path;
        fatal(ss.str());
    }

    std::string line;

    if (!std::getline(fin, line)) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: empty CSV: " << path;
        fatal(ss.str());
    }

    csv.header = split_csv_line(line);
    csv.index.clear();

    for (int i = 0; i < (int)csv.header.size(); ++i) {
        if (csv.index.find(csv.header[i]) != csv.index.end()) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: duplicate CSV column: '" << csv.header[i] << "'";
            fatal(ss.str());
        }

        csv.index[csv.header[i]] = i;
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
            ss << "[total_counts] FATAL: CSV row width mismatch while loading '"
               << path << "'. Row has " << row.size()
               << " cells, header has " << csv.header.size() << ".";
            fatal(ss.str());
        }

        csv.rows.push_back(std::move(row));
    }

    return true;
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";

    std::ofstream fout(tmp);

    if (!fout.is_open()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cannot write temp CSV: " << tmp;
        fatal(ss.str());
    }

    auto write_cell = [&](const std::string& s) {
        const bool needq =
            (s.find(',') != std::string::npos) ||
            (s.find('"') != std::string::npos) ||
            (s.find('\n') != std::string::npos) ||
            (s.find('\r') != std::string::npos);

        if (!needq) {
            fout << s;
            return;
        }

        fout << '"';

        for (char ch : s) {
            if (ch == '"') {
                fout << "\"\"";
            } else {
                fout << ch;
            }
        }

        fout << '"';
    };

    for (size_t i = 0; i < csv.header.size(); ++i) {
        write_cell(csv.header[i]);

        if (i + 1 < csv.header.size()) {
            fout << ',';
        }
    }

    fout << "\n";

    for (const auto& row : csv.rows) {
        if (row.size() != csv.header.size()) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: CSV row width mismatch during write. Row has "
               << row.size() << " cells, header has " << csv.header.size() << ".";
            fatal(ss.str());
        }

        for (size_t i = 0; i < row.size(); ++i) {
            write_cell(row[i]);

            if (i + 1 < row.size()) {
                fout << ',';
            }
        }

        fout << "\n";
    }

    fout.close();

    if (!fout) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: write failed for temp CSV: " << tmp;
        fatal(ss.str());
    }

    (void)std::remove(path.c_str());

    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: rename failed from '" << tmp
           << "' to '" << path << "'";
        fatal(ss.str());
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);

    if (it == csv.index.end()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: missing required CSV column: '" << name << "'";
        fatal(ss.str());
    }

    return it->second;
}

static inline std::string col_data_raw_counts(const ChannelConfig& channel_cfg,
                                              const std::string& topo_label,
                                              const std::string& period_display,
                                              const std::string& helicity) {
    return std::string("raw yield, ") + channel_cfg.csv_channel + ", " +
           topo_label + ", exp, " + period_display + ", " + helicity;
}

static inline std::string col_mc_generated(const ChannelConfig& channel_cfg,
                                           const std::string& period_display) {
    return std::string("generated yield, ") + channel_cfg.csv_channel +
           ", mc, " + period_display;
}

static inline std::string col_mc_reconstructed_total(const ChannelConfig& channel_cfg,
                                                     const std::string& period_display) {
    return std::string("reconstructed yield, ") + channel_cfg.csv_channel +
           ", mc, " + period_display;
}

static inline std::string col_mc_reconstructed_topo(const ChannelConfig& channel_cfg,
                                                    const std::string& topo_label,
                                                    const std::string& period_display) {
    return std::string("reconstructed yield, ") + channel_cfg.csv_channel + ", " +
           topo_label + ", mc, " + period_display;
}

static inline std::string col_data_normalized_counts(const ChannelConfig& channel_cfg,
                                                     const std::string& topo_label,
                                                     const std::string& period_display,
                                                     const std::string& helicity) {
    return std::string("normalized raw yield, ") + channel_cfg.csv_channel + ", " +
           topo_label + ", exp, " + period_display + ", " + helicity;
}

static inline std::string col_mc_current_corrected_total(const ChannelConfig& channel_cfg,
                                                         const std::string& period_display) {
    return std::string("reconstructed current corrected yield, ") + channel_cfg.csv_channel +
           ", mc, " + period_display;
}

static inline std::string col_mc_current_corrected_topo(const ChannelConfig& channel_cfg,
                                                        const std::string& topo_label,
                                                        const std::string& period_display) {
    return std::string("reconstructed current corrected yield, ") + channel_cfg.csv_channel + ", " +
           topo_label + ", mc, " + period_display;
}

// -----------------------------------------------------------------------------
// Row bins and fast lookup
// -----------------------------------------------------------------------------

struct RowBin {
    double xBmin = std::numeric_limits<double>::quiet_NaN();
    double xBmax = std::numeric_limits<double>::quiet_NaN();
    double Q2min = std::numeric_limits<double>::quiet_NaN();
    double Q2max = std::numeric_limits<double>::quiet_NaN();
    double tmin  = std::numeric_limits<double>::quiet_NaN();
    double tmax  = std::numeric_limits<double>::quiet_NaN();
    double pmin  = std::numeric_limits<double>::quiet_NaN();
    double pmax  = std::numeric_limits<double>::quiet_NaN();
    bool valid = false;
};

struct AxisBin {
    double min = 0.0;
    double max = 0.0;
};

struct FastBinning {
    std::vector<AxisBin> xbins;
    std::vector<AxisBin> qbins;
    std::vector<AxisBin> tbins;

    std::vector<std::vector<std::vector<std::vector<int>>>> rows_by_xqt;
};

static inline double to_double_strict(const std::string& s, const std::string& what) {
    if (s.empty()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: empty numeric cell for '" << what << "'";
        fatal(ss.str());
    }

    char* e = nullptr;
    const double v = std::strtod(s.c_str(), &e);

    if (e == s.c_str()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: parse failure for '" << what
           << "' value '" << s << "'";
        fatal(ss.str());
    }

    return v;
}

static inline bool to_bool_valid(const std::string& s) {
    return (s == "1" || s == "1.0" || s == "true" || s == "TRUE");
}

static std::vector<RowBin> load_row_bins_from_csv(const CSV& csv) {
    const int c_xBmin = col_strict(csv, "xBmin");
    const int c_xBmax = col_strict(csv, "xBmax");
    const int c_Q2min = col_strict(csv, "Q2min");
    const int c_Q2max = col_strict(csv, "Q2max");
    const int c_tmin  = col_strict(csv, "t_abs_min");
    const int c_tmax  = col_strict(csv, "t_abs_max");
    const int c_pmin  = col_strict(csv, "phimin");
    const int c_pmax  = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    std::vector<RowBin> rows;
    rows.reserve(csv.rows.size());

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        const auto& row = csv.rows[r];

        RowBin b;
        b.xBmin = to_double_strict(row[c_xBmin], "xBmin");
        b.xBmax = to_double_strict(row[c_xBmax], "xBmax");
        b.Q2min = to_double_strict(row[c_Q2min], "Q2min");
        b.Q2max = to_double_strict(row[c_Q2max], "Q2max");
        b.tmin  = to_double_strict(row[c_tmin],  "t_abs_min");
        b.tmax  = to_double_strict(row[c_tmax],  "t_abs_max");
        b.pmin  = to_double_strict(row[c_pmin],  "phimin");
        b.pmax  = to_double_strict(row[c_pmax],  "phimax");
        b.valid = to_bool_valid(row[c_valid]);

        rows.push_back(b);
    }

    return rows;
}

static inline bool axis_bin_equal(const AxisBin& a, const AxisBin& b) {
    return (a.min == b.min && a.max == b.max);
}

static void add_unique_axis_bin(std::vector<AxisBin>& bins, double minv, double maxv) {
    AxisBin b;
    b.min = minv;
    b.max = maxv;

    auto it = std::find_if(bins.begin(), bins.end(), [&](const AxisBin& x) {
        return axis_bin_equal(x, b);
    });

    if (it == bins.end()) {
        bins.push_back(b);
    }
}

static void sort_axis_bins(std::vector<AxisBin>& bins) {
    std::sort(bins.begin(), bins.end(), [](const AxisBin& a, const AxisBin& b) {
        if (a.min != b.min) {
            return a.min < b.min;
        }

        return a.max < b.max;
    });
}

static int find_axis_bin_index(const std::vector<AxisBin>& bins, double value) {
    for (int i = 0; i < (int)bins.size(); ++i) {
        if (value >= bins[i].min && value < bins[i].max) {
            return i;
        }
    }

    return -1;
}

static int find_axis_bin_exact(const std::vector<AxisBin>& bins, double minv, double maxv) {
    for (int i = 0; i < (int)bins.size(); ++i) {
        if (bins[i].min == minv && bins[i].max == maxv) {
            return i;
        }
    }

    return -1;
}

static FastBinning build_fast_binning(const std::vector<RowBin>& rows) {
    FastBinning fb;

    for (const auto& r : rows) {
        if (!r.valid) {
            continue;
        }

        add_unique_axis_bin(fb.xbins, r.xBmin, r.xBmax);
        add_unique_axis_bin(fb.qbins, r.Q2min, r.Q2max);
        add_unique_axis_bin(fb.tbins, r.tmin,  r.tmax);
    }

    sort_axis_bins(fb.xbins);
    sort_axis_bins(fb.qbins);
    sort_axis_bins(fb.tbins);

    fb.rows_by_xqt.resize(fb.xbins.size());

    for (size_t ix = 0; ix < fb.xbins.size(); ++ix) {
        fb.rows_by_xqt[ix].resize(fb.qbins.size());

        for (size_t iq = 0; iq < fb.qbins.size(); ++iq) {
            fb.rows_by_xqt[ix][iq].resize(fb.tbins.size());
        }
    }

    for (int r = 0; r < (int)rows.size(); ++r) {
        const RowBin& row = rows[r];

        if (!row.valid) {
            continue;
        }

        const int ix = find_axis_bin_exact(fb.xbins, row.xBmin, row.xBmax);
        const int iq = find_axis_bin_exact(fb.qbins, row.Q2min, row.Q2max);
        const int it = find_axis_bin_exact(fb.tbins, row.tmin,  row.tmax);

        if (ix < 0 || iq < 0 || it < 0) {
            fatal("[total_counts] FATAL: failed to build fast row bin lookup.");
        }

        fb.rows_by_xqt[ix][iq][it].push_back(r);
    }

    std::cout << "[total_counts] Fast bin lookup built with "
              << fb.xbins.size() << " xB bins, "
              << fb.qbins.size() << " Q2 bins, "
              << fb.tbins.size() << " |t| bins."
              << std::endl;

    return fb;
}

// -----------------------------------------------------------------------------
// Combined 3-exclusivity cuts loader
// -----------------------------------------------------------------------------

struct SigmaStats {
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std  = std::numeric_limits<double>::quiet_NaN();
    double cut_low = std::numeric_limits<double>::quiet_NaN();
    double cut_high = std::numeric_limits<double>::quiet_NaN();
    double quantile = 0.0;
    std::string mode = "symmetric_3sigma";
};

using CutVarMap = std::unordered_map<std::string, SigmaStats>;
using TopoCutMap = std::unordered_map<std::string, CutVarMap>;

static inline bool within_cut_window(double v, const SigmaStats& s) {
    if (!std::isfinite(v)) {
        return false;
    }

    if (s.mode == "upper_quantile") {
        if (!std::isfinite(s.cut_high)) {
            return true;
        }
        return v <= s.cut_high;
    }

    double lo = s.cut_low;
    double hi = s.cut_high;

    if (!(std::isfinite(lo) && std::isfinite(hi)) || hi <= lo) {
        if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) {
            return true;
        }
        lo = s.mean - 3.0 * s.std;
        hi = s.mean + 3.0 * s.std;
    }

    return (v >= lo && v <= hi);
}

static TopoCutMap load_combined_cuts(const std::string& combined_cuts_json,
                                     const std::string& sample_key) {
    std::ifstream fin(combined_cuts_json);

    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cannot open combined cuts JSON: "
           << combined_cuts_json;
        fatal(ss.str());
    }

    nlohmann::json j;

    try {
        fin >> j;
    } catch (const std::exception& e) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: JSON parse failed for "
           << combined_cuts_json << " : " << e.what();
        fatal(ss.str());
    }

    if (!j.is_object()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: combined cuts JSON is not an object: "
           << combined_cuts_json;
        fatal(ss.str());
    }

    TopoCutMap out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const auto& block = it.value();

        if (!block.is_object()) {
            continue;
        }

        if (!block.contains(sample_key)) {
            continue;
        }

        const auto& data = block[sample_key];

        if (!data.is_object()) {
            continue;
        }

        CutVarMap vm;

        for (auto vit = data.begin(); vit != data.end(); ++vit) {
            const std::string var = vit.key();
            const auto& stats = vit.value();

            if (!stats.is_object()) {
                continue;
            }

            if (!stats.contains("mean") || !stats.contains("std")) {
                continue;
            }

            SigmaStats s;

            try {
                s.mean = stats["mean"].get<double>();
                s.std  = stats["std"].get<double>();

                if (stats.contains("cut_low")) {
                    s.cut_low = stats["cut_low"].get<double>();
                }
                if (stats.contains("cut_high")) {
                    s.cut_high = stats["cut_high"].get<double>();
                }
                if (stats.contains("quantile")) {
                    s.quantile = stats["quantile"].get<double>();
                }
                if (stats.contains("mode")) {
                    s.mode = stats["mode"].get<std::string>();
                }
            } catch (...) {
                continue;
            }

            if (!std::isfinite(s.cut_low) || !std::isfinite(s.cut_high) || s.cut_high <= s.cut_low) {
                if (std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0) {
                    s.cut_low = s.mean - 3.0 * s.std;
                    s.cut_high = s.mean + 3.0 * s.std;
                }
            }

            if (std::isfinite(s.cut_high)) {
                vm.emplace(var, s);
            }
        }

        if (!vm.empty()) {
            out.emplace(key, std::move(vm));
        }
    }

    std::cout << "[total_counts] Loaded " << sample_key << " exclusivity cuts for "
              << out.size() << " topology keys from " << combined_cuts_json
              << std::endl;

    return out;
}

// -----------------------------------------------------------------------------
// Branch binder
// -----------------------------------------------------------------------------

static std::mutex g_root_bind_mutex;

struct BranchBinder {
    int runnum = 0; bool has_runnum = false;

    int detector1 = 0; bool has_detector1 = false;
    int detector2 = 0; bool has_detector2 = false;

    int helicity = 0; bool has_helicity = false;

    double x = 0.0;     bool has_x = false;
    double Q2 = 0.0;    bool has_Q2 = false;
    double t1 = 0.0;    bool has_t1 = false;
    double phi2 = 0.0;  bool has_phi2 = false;
    double open_angle_ep2 = 0.0; bool has_open_angle = false;

    double Delta_phi = 0.0;         bool has_Delta_phi = false;
    double theta = 0.0;             bool has_theta = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0;     bool has_theta_pi0_pi0 = false;
    double pTmiss = 0.0;            bool has_pTmiss = false;
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;

    double e_p = 0.0;       bool has_e_p = false;
    double e_theta = 0.0;   bool has_e_theta = false;
    double e_phi = 0.0;     bool has_e_phi = false;

    double p1_theta = 0.0;  bool has_p1_theta = false;
    double p1_phi = 0.0;    bool has_p1_phi = false;

    double p2_p = 0.0;      bool has_p2_p = false;
    double p2_theta = 0.0;  bool has_p2_theta = false;
    double p2_phi = 0.0;    bool has_p2_phi = false;

    void bind(TTree* t, const WorkConfig& work_cfg) {
        if (!t) {
            return;
        }

        std::lock_guard<std::mutex> lock(g_root_bind_mutex);

        t->SetBranchStatus("*", 0);

        // Use a modest per-tree read cache. Only explicitly enabled branches are
        // added below, so generated MC does not spend I/O on reconstructed-only
        // quantities. With seven workers this remains bounded at about 112 MiB.
        static constexpr Long64_t kTreeCacheBytes = 16LL * 1024LL * 1024LL;
        t->SetCacheSize(kTreeCacheBytes);

        auto ena = [&](const char* n) {
            if (t->GetBranch(n)) {
                t->SetBranchStatus(n, 1);
                t->AddBranchToCache(n, true);
            }
        };

        const bool is_gen = (work_cfg.sample_kind == SampleKind::MC_GEN);
        const bool is_data = (work_cfg.sample_kind == SampleKind::DATA);

        // Every sample needs only these four branches for CSV-bin matching.
        ena("x");
        ena("Q2");
        ena("t1");
        ena("phi2");

        if (!is_gen) {
            const GlobalCutConfig& cfg = default_global_cuts();

            ena("runnum");
            ena("detector1");
            ena("detector2");

            // Needed by the optional Sp18-Out centered theta_e current
            // response, independent of whether other global cuts require it.
            ena("e_theta");

            if (is_data) {
                ena("helicity");
            }

            // Global cuts always require open_angle_ep2. pTmiss is required by
            // its global cut and is also one of the exclusivity variables.
            ena("open_angle_ep2");
            ena("pTmiss");

            // Production exclusivity variables emitted by the Python
            // optimization. The direct-pi0 JSON retains the logical
            // theta_gamma_gamma key, but its physical branch is
            // theta_pi0_pi0.
            ena("Delta_phi");
            ena("theta");
            ena("pTmiss");
            ena("Emiss2");
            ena("Mx2");
            ena("Mx2_2");

            if (work_cfg.channel_cfg.channel == Channel::EPPI0) {
                ena("theta_pi0_pi0");
            } else {
                ena("theta_gamma_gamma");
            }

            // Delta_phi falls back to p1_phi/p2_phi if the direct branch is
            // absent. These phi branches are also required by sector cuts.
            ena("p1_phi");
            ena("p2_phi");

            if (global_cuts_require_sector_phi(cfg) ||
                global_cuts_require_auxiliary_kinematics(cfg) ||
                cfg.enable_dvcsgen_ycol_cut) {
                ena("e_phi");
            }

            if (global_cuts_require_auxiliary_kinematics(cfg) ||
                cfg.enable_dvcsgen_ycol_cut) {
                ena("e_p");
                ena("e_theta");
                ena("p2_p");
                ena("p2_theta");
            }

            if (global_cuts_require_auxiliary_kinematics(cfg)) {
                ena("p1_theta");
            }
        }

        t->StopCacheLearningPhase();

        auto bI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n) && t->GetBranchStatus(n)) {
                t->SetBranchAddress(n, a);
                f = true;
            }
        };

        auto bD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n) && t->GetBranchStatus(n)) {
                t->SetBranchAddress(n, a);
                f = true;
            }
        };

        bI("runnum", &runnum, has_runnum);

        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);
        bI("helicity",  &helicity,  has_helicity);

        bD("x",    &x,    has_x);
        bD("Q2",   &Q2,   has_Q2);
        bD("t1",   &t1,   has_t1);
        bD("phi2", &phi2, has_phi2);
        bD("Delta_phi", &Delta_phi, has_Delta_phi);

        bD("open_angle_ep2", &open_angle_ep2, has_open_angle);

        bD("theta",             &theta,             has_theta);
        bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bD("theta_pi0_pi0",     &theta_pi0_pi0,     has_theta_pi0_pi0);
        bD("pTmiss",            &pTmiss,            has_pTmiss);
        bD("Emiss2",            &Emiss2,            has_Emiss2);
        bD("Mx2",               &Mx2,               has_Mx2);
        bD("Mx2_2",             &Mx2_2,             has_Mx2_2);

        bD("e_p",     &e_p,     has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi",   &e_phi,   has_e_phi);
        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi",  &p1_phi,  has_p1_phi);

        bD("p2_p",     &p2_p,     has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi",   &p2_phi,   has_p2_phi);
    }

    bool ready_for_generated_matching() const {
        return has_x && has_Q2 && has_t1 && has_phi2;
    }

    bool ready_for_reconstructed_matching() const {
        return has_detector1 && has_detector2 && has_x && has_Q2 && has_t1 && has_phi2;
    }

    double phi_deg() const {
        return wrap_phi_deg(phi2 * RAD2DEG);
    }

    double t_abs() const {
        return std::fabs(t1);
    }

    double delta_phi_value(bool& has_val) const {
        if (has_Delta_phi) {
            has_val = true;
            return Delta_phi;
        }

        if (has_p1_phi && has_p2_phi) {
            has_val = true;
            return delta_phi_rad_from_two_phi(p1_phi, p2_phi);
        }

        has_val = false;
        return 0.0;
    }
};

// -----------------------------------------------------------------------------
// Counts accumulation
// -----------------------------------------------------------------------------

struct HelCounts {
    double unpol = 0.0;
    double pos   = 0.0;
    double neg   = 0.0;
};

using RowCounts = std::unordered_map<int, HelCounts>;

struct WeightedHelCounts {
    HelCounts sumw;
    HelCounts sumw2;
    std::array<HelCounts, kCurrentRegionCount> derivative_by_region{};
    std::array<double, kCurrentRegionCount> parameter_stat{};
    HelCounts angular_derivative;
    double angular_parameter_stat = 0.0;
};

using WeightedRowCounts = std::unordered_map<int, WeightedHelCounts>;

struct CutFlowSummary {
    long long entries = 0;
    long long valid_topology = 0;
    long long global_pass = 0;
    long long sigma_pass = 0;
    long long matched = 0;

    std::unordered_map<std::string, long long> topology_entries;
    std::unordered_map<std::string, long long> topology_global_pass;
    std::unordered_map<std::string, long long> topology_sigma_pass;
    std::unordered_map<std::string, long long> topology_matched;

    // Diagnostic only. These are evaluated after topology + global cuts.
    // sigma_single_pass[var] = number of events passing that one exclusivity cut alone.
    // sigma_cumulative_pass[var] = number of events surviving the ordered cumulative
    // sequence up to and including that variable.
    std::unordered_map<std::string, long long> sigma_single_pass;
    std::unordered_map<std::string, long long> sigma_cumulative_pass;
    std::unordered_map<std::string, std::unordered_map<std::string, long long>> topology_sigma_single_pass;
    std::unordered_map<std::string, std::unordered_map<std::string, long long>> topology_sigma_cumulative_pass;
};

struct WorkCounts {
    RowCounts total_counts;
    std::unordered_map<std::string, RowCounts> topo_counts;
    WeightedRowCounts corrected_total_counts;
    std::unordered_map<std::string, WeightedRowCounts> corrected_topo_counts;
    CutFlowSummary flow;
};

static constexpr std::size_t kProductionVariableCount = 7;

static const std::array<std::string, kProductionVariableCount>&
production_variable_order() {
    static const std::array<std::string, kProductionVariableCount> vars = {
        "Delta_phi",
        "theta",
        "theta_gamma_gamma",
        "pTmiss",
        "Emiss2",
        "Mx2",
        "Mx2_2"
    };
    return vars;
}

static inline double branch_value_for_production_variable(
    const BranchBinder& b,
    const std::string& variable,
    bool use_pi0_angle,
    bool& has_val) {

    has_val = true;

    if (variable == "Delta_phi") {
        return b.delta_phi_value(has_val);
    }
    if (variable == "theta") {
        has_val = b.has_theta;
        return b.theta;
    }
    if (variable == "theta_gamma_gamma") {
        if (use_pi0_angle) {
            has_val = b.has_theta_pi0_pi0;
            return b.theta_pi0_pi0;
        }
        has_val = b.has_theta_gamma_gamma;
        return b.theta_gamma_gamma;
    }
    if (variable == "pTmiss") {
        has_val = b.has_pTmiss;
        return b.pTmiss;
    }
    if (variable == "Emiss2") {
        has_val = b.has_Emiss2;
        return b.Emiss2;
    }
    if (variable == "Mx2") {
        has_val = b.has_Mx2;
        return b.Mx2;
    }
    if (variable == "Mx2_2") {
        has_val = b.has_Mx2_2;
        return b.Mx2_2;
    }

    has_val = false;
    return 0.0;
}

struct CompiledSigmaPlan {
    std::array<const SigmaStats*, kProductionVariableCount> stats{};
    std::array<std::string, kProductionVariableCount> names{};
    std::string key;
    bool use_pi0_angle = false;
};

static CompiledSigmaPlan compile_sigma_plan(const ChannelConfig& channel_cfg,
                                            const TopoCutMap& cuts,
                                            const std::string& key) {
    CompiledSigmaPlan plan;
    plan.key = key;
    plan.use_pi0_angle = (channel_cfg.channel == Channel::EPPI0);
    plan.names = production_variable_order();

    const auto it = cuts.find(key);
    if (it == cuts.end()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: missing exclusivity-cut key in "
           << "combined_cuts.json: '" << key << "'";
        fatal(ss.str());
    }

    const CutVarMap& vm = it->second;
    for (const auto& entry : vm) {
        const auto allowed = std::find(
            plan.names.begin(), plan.names.end(), entry.first);
        if (allowed == plan.names.end()) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: cut key '" << key
               << "' contains unsupported production variable '"
               << entry.first << "'. Expected only Delta_phi, theta, "
               << "theta_gamma_gamma, pTmiss, Emiss2, Mx2 and Mx2_2.";
            fatal(ss.str());
        }
    }

    for (std::size_t i = 0; i < kProductionVariableCount; ++i) {
        const auto iv = vm.find(plan.names[i]);
        if (iv != vm.end()) {
            plan.stats[i] = &iv->second;
        }
    }

    // Mx2 is the forced first production cut in the Python optimizer and must
    // therefore always be present. The remaining variables are optional
    // because the optimizer stops when an additional cut no longer improves
    // S/sqrt(S+B) with positive pi0 discrimination.
    const auto mx2 = vm.find("Mx2");
    if (mx2 == vm.end()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: cut key '" << key
           << "' is missing the mandatory first production cut 'Mx2'.";
        fatal(ss.str());
    }

    return plan;
}

struct DenseSigmaDiagnostics {
    std::array<long long, kProductionVariableCount> single{};
    std::array<long long, kProductionVariableCount> cumulative{};
    std::array<std::array<long long, kProductionVariableCount>, 3> topo_single{};
    std::array<std::array<long long, kProductionVariableCount>, 3> topo_cumulative{};
};

static inline bool fill_sigma_cut_diagnostics_compiled(
    const CompiledSigmaPlan& plan,
    const BranchBinder& b,
    int topo_index_value,
    DenseSigmaDiagnostics& diag) {

    bool cumulative_ok = true;
    bool all_ok = true;

    for (std::size_t i = 0; i < kProductionVariableCount; ++i) {
        const SigmaStats* cut = plan.stats[i];

        // A null entry means the Python production optimizer did not accept
        // this variable. It is not a cut and must not alter the cumulative
        // diagnostics.
        if (cut == nullptr) {
            continue;
        }

        bool has_val = false;
        const double val = branch_value_for_production_variable(
            b, plan.names[i], plan.use_pi0_angle, has_val);

        if (!has_val) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: cut key '" << plan.key
               << "' requires variable '" << plan.names[i]
               << "', but the corresponding branch is missing in this tree.";
            if (plan.use_pi0_angle &&
                plan.names[i] == "theta_gamma_gamma") {
                ss << " Direct eppi0 samples require branch "
                   << "'theta_pi0_pi0' for this logical cut.";
            }
            fatal(ss.str());
        }

        const bool pass_this = within_cut_window(val, *cut);
        if (pass_this) {
            ++diag.single[i];
            ++diag.topo_single[topo_index_value][i];
        } else {
            all_ok = false;
        }

        if (cumulative_ok && pass_this) {
            ++diag.cumulative[i];
            ++diag.topo_cumulative[topo_index_value][i];
        } else {
            cumulative_ok = false;
        }
    }

    return all_ok;
}

static void materialize_sigma_diagnostics(
    const ChannelConfig&,
    const DenseSigmaDiagnostics& diag,
    CutFlowSummary& flow) {

    static const std::array<TopologyIndex, 3> kTopologies = {
        TopologyIndex::FD_FD,
        TopologyIndex::CD_FD,
        TopologyIndex::CD_FT
    };

    const auto& vars = production_variable_order();

    for (std::size_t i = 0; i < kProductionVariableCount; ++i) {
        const std::string& var = vars[i];

        if (diag.single[i] != 0) {
            flow.sigma_single_pass[var] = diag.single[i];
        }
        if (diag.cumulative[i] != 0) {
            flow.sigma_cumulative_pass[var] = diag.cumulative[i];
        }

        for (int ti = 0; ti < 3; ++ti) {
            const std::string& topo = topology_name(kTopologies[ti]);
            if (diag.topo_single[ti][i] != 0) {
                flow.topology_sigma_single_pass[topo][var] =
                    diag.topo_single[ti][i];
            }
            if (diag.topo_cumulative[ti][i] != 0) {
                flow.topology_sigma_cumulative_pass[topo][var] =
                    diag.topo_cumulative[ti][i];
            }
        }
    }
}

static inline bool passes_global_cuts_dispatch(const BranchBinder& b,
                                               const std::string& period_label) {
    const GlobalCutConfig& cfg = default_global_cuts();

    if (!(b.has_t1 && b.has_open_angle)) return false;
    if (cfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (b.has_runnum && is_excluded_run(b.runnum)) return false;

    if (cfg.enable_topology_filter || global_cuts_require_sector_phi(cfg) || cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_detector1 && b.has_detector2)) {
            fatal("[total_counts] FATAL: topology/sector/global-ycol selection requires detector1/detector2 branches.");
        }
    }

    if (global_cuts_require_sector_phi(cfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            fatal("[total_counts] FATAL: sector selection requires e_phi, p1_phi, and p2_phi branches.");
        }
    }

    if (global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[total_counts] FATAL: auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta, p2_phi branches.");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (cfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[total_counts] FATAL: dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi branches.");
        }

        if (global_cuts_require_sector_phi(cfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      period_label,
                                      b.e_p, b.e_theta, b.e_phi, b.p1_theta, b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      cfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
    }

    if (global_cuts_require_sector_phi(cfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  cfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              cfg);
}

static inline void add_count(HelCounts& h, bool split_helicity, int helicity) {
    if (!split_helicity) {
        h.unpol += 1.0;
        return;
    }

    if (helicity > 0) {
        h.pos += 1.0;
    } else if (helicity < 0) {
        h.neg += 1.0;
    } else {
        h.unpol += 1.0;
    }
}


static inline void add_weighted_count(WeightedHelCounts& h,
                                      bool split_helicity,
                                      int helicity,
                                      const EventCurrentWeight& cw) {
    if (cw.region < 0 || cw.region >= kCurrentRegionCount) return;
    auto add_one = [&](double& sw, double& sw2, double& deriv, double& aderiv) {
        sw += cw.weight;
        sw2 += cw.weight * cw.weight;
        deriv += cw.derivative;
        aderiv += cw.angular_derivative;
    };
    h.parameter_stat[cw.region] = cw.parameter_stat;
    h.angular_parameter_stat =
        std::max(h.angular_parameter_stat, cw.angular_parameter_stat);
    if (!split_helicity) {
        add_one(h.sumw.unpol, h.sumw2.unpol,
                h.derivative_by_region[cw.region].unpol,
                h.angular_derivative.unpol);
    } else if (helicity > 0) {
        add_one(h.sumw.pos, h.sumw2.pos,
                h.derivative_by_region[cw.region].pos,
                h.angular_derivative.pos);
    } else if (helicity < 0) {
        add_one(h.sumw.neg, h.sumw2.neg,
                h.derivative_by_region[cw.region].neg,
                h.angular_derivative.neg);
    } else {
        add_one(h.sumw.unpol, h.sumw2.unpol,
                h.derivative_by_region[cw.region].unpol,
                h.angular_derivative.unpol);
    }
}

static EventCurrentWeight event_current_weight(const CurrentResponseModel& model,
                                               const WorkConfig& work_cfg,
                                               const PeriodTags& tags,
                                               const BranchBinder& b,
                                               bool use_epg_mc_factor_for_eppi0_bkg) {
    EventCurrentWeight out;
    const int region = current_region_index(b.detector2, b.p2_phi, b.has_p2_phi);
    if (region < 0) {
        fatal("[total_counts] FATAL: cannot determine FT/FD sector for current correction in tree '" + tags.tree_key + "'.");
    }
    out.region = region;

    if (work_cfg.sample_kind == SampleKind::DATA) {
        const CurrentResponseEntry* e = find_current_response_entry(
            model.data, work_cfg.channel_cfg.csv_channel, tags.period_display, region);
        if (!e || !e->valid) {
            fatal("[total_counts] FATAL: missing DATA regional current response for channel='" +
                  work_cfg.channel_cfg.csv_channel + "' period='" + tags.period_display +
                  "' region='" + current_region_names()[region] + "'.");
        }
        auto ip = model.run_current_nA.find(tags.period_display);
        if (ip == model.run_current_nA.end()) {
            fatal("[total_counts] FATAL: no run-current map for period '" + tags.period_display + "'.");
        }
        auto ir = ip->second.find(b.runnum);
        if (ir == ip->second.end()) {
            const auto iex = model.excluded_data_runs.find(tags.period_display);
            const bool explicitly_excluded =
                (iex != model.excluded_data_runs.end() && iex->second.count(b.runnum) != 0);
            if (explicitly_excluded) {
                out.skip_corrected = true;
                out.weight = 0.0;
                out.derivative = 0.0;
                out.parameter_stat = 0.0;
                return out;
            }
            fatal("[total_counts] FATAL: run " + std::to_string(b.runnum) +
                  " has no current assignment in response model for period '" + tags.period_display +
                  "' and is not explicitly marked as an unusable/zero-charge run.");
        }
        const double I = (double)ir->second;
        double effective_slope = e->parameter;
        double angular_offset = 0.0;
        if (e->has_angular_model) {
            if (e->angular_variable != "e_theta" || !b.has_e_theta ||
                !std::isfinite(b.e_theta)) {
                fatal("[total_counts] FATAL: DATA angular current response requires a finite e_theta branch.");
            }
            const double theta_deg = b.e_theta * RAD2DEG;
            angular_offset = theta_deg - e->angular_center;
            effective_slope += e->angular_gradient * angular_offset;
        }

        const double R = 1.0 + effective_slope * I;
        if (!(std::isfinite(R) && R > 0.0)) {
            fatal("[total_counts] FATAL: non-positive DATA current response encountered.");
        }
        out.weight = 1.0 / R;
        out.derivative = -I / (R * R);
        out.parameter_stat = e->parameter_stat;
        if (e->has_angular_model) {
            out.angular_derivative = -I * angular_offset / (R * R);
            out.angular_parameter_stat = e->angular_gradient_stat;
        }
        return out;
    }

    if (work_cfg.sample_kind == SampleKind::MC_REC) {
        std::string response_channel = work_cfg.channel_cfg.csv_channel;
        if (work_cfg.channel_cfg.channel == Channel::EPPI0_BKG_AS_DVCS) {
            response_channel = use_epg_mc_factor_for_eppi0_bkg ? "ep->epg" : "ep->eppi0";
        }
        const CurrentResponseEntry* e = find_current_response_entry(
            model.mc, response_channel, tags.period_display, region);
        if (!e || !e->valid || !(e->parameter > 0.0)) {
            fatal("[total_counts] FATAL: missing MC regional current response for channel='" +
                  response_channel + "' period='" + tags.period_display +
                  "' region='" + current_region_names()[region] + "'.");
        }
        const double f = e->parameter;
        out.weight = 1.0 / f;
        out.derivative = -1.0 / (f * f);
        out.parameter_stat = e->parameter_stat;
    }
    return out;
}

static WorkCounts accumulate_counts_for_tree(const WorkConfig& work_cfg,
                                             const PeriodTags& tags,
                                             TTree* tree,
                                             const std::vector<RowBin>& rows,
                                             const FastBinning& fast_bins,
                                             const TopoCutMap& sigma_cuts,
                                             bool trace_matches,
                                             const CurrentResponseModel* current_model,
                                             bool use_epg_mc_current_factor_for_eppi0_bkg) {
    WorkCounts out;

    if (!tree) {
        return out;
    }

    BranchBinder b;
    b.bind(tree, work_cfg);

    const bool is_gen = (work_cfg.sample_kind == SampleKind::MC_GEN);
    const bool is_data = (work_cfg.sample_kind == SampleKind::DATA);
    const bool split_helicity = is_data;

    if (is_gen) {
        if (!b.ready_for_generated_matching()) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: missing required generated-MC matching branches in tree for '"
               << tags.tree_key << "'. Required: x, Q2, t1, phi2.";
            fatal(ss.str());
        }
    } else {
        if (!b.ready_for_reconstructed_matching()) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: missing required reconstructed matching branches in tree for '"
               << tags.tree_key << "'. Required: detector1, detector2, x, Q2, t1, phi2.";
            fatal(ss.str());
        }

        if (is_data && !b.has_helicity) {
            std::ostringstream ss;
            ss << "[total_counts] FATAL: missing required branch 'helicity' in data tree for '"
               << tags.tree_key << "'.";
            fatal(ss.str());
        }
    }

    // Dense row-indexed storage avoids two unordered-map lookups for every
    // matched event. The public/downstream representation is rebuilt once at
    // the end of the tree, so all existing interfaces and CSV behavior remain
    // unchanged.
    std::vector<HelCounts> total_dense(rows.size());
    std::array<std::vector<HelCounts>, 3> topo_dense;
    std::vector<WeightedHelCounts> corrected_total_dense;
    std::array<std::vector<WeightedHelCounts>, 3> corrected_topo_dense;
    const bool apply_current_weights = (!is_gen && current_model != nullptr && tags.period_display != "Fa18 Inb Supp");
    if (apply_current_weights) corrected_total_dense.resize(rows.size());
    if (!is_gen) {
        for (auto& v : topo_dense) v.resize(rows.size());
        if (apply_current_weights) {
            for (auto& v : corrected_topo_dense) v.resize(rows.size());
        }
    }

    // The sigma-cut key depends only on channel, period, and topology. Build
    // the three strings once per tree instead of once per accepted event.
    std::array<std::string, 3> sigma_keys;
    std::array<CompiledSigmaPlan, 3> sigma_plans;
    std::array<bool, 3> sigma_plan_ready{{false, false, false}};
    DenseSigmaDiagnostics sigma_diag;
    if (!is_gen) {
        sigma_keys[0] = combined_cuts_key(work_cfg.channel_cfg, tags, topology_name(TopologyIndex::FD_FD));
        sigma_keys[1] = combined_cuts_key(work_cfg.channel_cfg, tags, topology_name(TopologyIndex::CD_FD));
        sigma_keys[2] = combined_cuts_key(work_cfg.channel_cfg, tags, topology_name(TopologyIndex::CD_FT));
    }

    const Long64_t N = tree->GetEntries();
    const bool dbg = env_flag("TOTAL_COUNTS_DEBUG");

    long long n_global_pass = 0;
    long long n_sigma_pass  = 0;
    long long n_used        = 0;
    long long n_current_excluded = 0;
    std::unordered_set<int> current_excluded_runs_seen;

    out.flow.entries = (long long)N;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        TopologyIndex topo = TopologyIndex::INVALID;
        const std::string* topo_dir_ptr = nullptr;
        int topo_idx = -1;

        if (!is_gen) {
            topo = topology_index(b.detector1, b.detector2);

            if (topo == TopologyIndex::INVALID) {
                continue;
            }

            topo_idx = static_cast<int>(topo);
            const std::string& topoDir = topology_name(topo);
            topo_dir_ptr = &topoDir;

            ++out.flow.valid_topology;
            ++out.flow.topology_entries[topoDir];

            if (!passes_global_cuts_dispatch(b, tags.period_label)) {
                continue;
            }

            ++n_global_pass;
            ++out.flow.global_pass;
            ++out.flow.topology_global_pass[topoDir];

            if (!sigma_plan_ready[topo_idx]) {
                sigma_plans[topo_idx] = compile_sigma_plan(work_cfg.channel_cfg,
                                                           sigma_cuts,
                                                           sigma_keys[topo_idx]);
                sigma_plan_ready[topo_idx] = true;
            }

            if (!fill_sigma_cut_diagnostics_compiled(sigma_plans[topo_idx],
                                                     b,
                                                     topo_idx,
                                                     sigma_diag)) {
                continue;
            }

            ++n_sigma_pass;
            ++out.flow.sigma_pass;
            ++out.flow.topology_sigma_pass[topoDir];
        }

        EventCurrentWeight current_weight;
        if (apply_current_weights) {
            current_weight = event_current_weight(*current_model, work_cfg, tags, b,
                                                  use_epg_mc_current_factor_for_eppi0_bkg);
            if (current_weight.skip_corrected) {
                ++n_current_excluded;
                current_excluded_runs_seen.insert(b.runnum);
            }
        }

        const double phi_deg = b.phi_deg();
        const double tabs = b.t_abs();

        const int ix = find_axis_bin_index(fast_bins.xbins, b.x);
        if (ix < 0) {
            continue;
        }

        const int iq = find_axis_bin_index(fast_bins.qbins, b.Q2);
        if (iq < 0) {
            continue;
        }

        const int it = find_axis_bin_index(fast_bins.tbins, tabs);
        if (it < 0) {
            continue;
        }

        const std::vector<int>& candidate_rows = fast_bins.rows_by_xqt[ix][iq][it];

        bool matched_any = false;

        for (int r : candidate_rows) {
            const RowBin& w = rows[r];

            if (!row_accepts_phi(phi_deg, w.pmin, w.pmax)) {
                continue;
            }

            add_count(total_dense[r], split_helicity, b.helicity);

            if (!is_gen) {
                add_count(topo_dense[topo_idx][r], split_helicity, b.helicity);
            }
            if (apply_current_weights && !current_weight.skip_corrected) {
                add_weighted_count(corrected_total_dense[r], split_helicity, b.helicity, current_weight);
                add_weighted_count(corrected_topo_dense[topo_idx][r], split_helicity, b.helicity, current_weight);
            }

            matched_any = true;

            if (trace_matches) {
                std::cout << "[total_counts][TRACE] channel=" << work_cfg.channel_cfg.csv_channel
                          << " sample=" << (is_gen ? "gen" : (is_data ? "data" : "rec"))
                          << " tree=" << tags.tree_key
                          << " topo=" << (is_gen ? std::string("GEN") : *topo_dir_ptr)
                          << " row=" << r
                          << " x=" << b.x
                          << " Q2=" << b.Q2
                          << " |t|=" << tabs
                          << " phi(deg)=" << phi_deg
                          << " hel=" << b.helicity
                          << std::endl;
            }
        }

        if (matched_any) {
            ++n_used;
            ++out.flow.matched;

            if (!is_gen) {
                ++out.flow.topology_matched[*topo_dir_ptr];
            }
        }

        if (dbg && i < 3) {
            std::cout << "[total_counts][DEBUG] channel=" << work_cfg.channel_cfg.csv_channel
                      << " tree=" << tags.tree_key
                      << " i=" << (long long)i
                      << " sample=" << (is_gen ? "gen" : (is_data ? "data" : "rec"))
                      << " topo=" << (is_gen ? std::string("GEN") : *topo_dir_ptr)
                      << " hel=" << b.helicity
                      << " x=" << b.x
                      << " Q2=" << b.Q2
                      << " t1=" << b.t1
                      << " phi2(rad)=" << b.phi2
                      << " phi(deg)=" << phi_deg
                      << std::endl;
        }
    }

    if (!is_gen) {
        materialize_sigma_diagnostics(work_cfg.channel_cfg, sigma_diag, out.flow);
    }

    auto nonzero = [](const HelCounts& h) {
        return h.unpol != 0.0 || h.pos != 0.0 || h.neg != 0.0;
    };

    for (int r = 0; r < (int)rows.size(); ++r) {
        if (nonzero(total_dense[r])) {
            out.total_counts.emplace(r, total_dense[r]);
        }
    }

    if (apply_current_weights) {
        for (int r = 0; r < (int)rows.size(); ++r) {
            const WeightedHelCounts& h = corrected_total_dense[r];
            if (h.sumw.unpol != 0.0 || h.sumw.pos != 0.0 || h.sumw.neg != 0.0) {
                out.corrected_total_counts.emplace(r, h);
            }
        }
    }

    if (!is_gen) {
        static const std::array<TopologyIndex, 3> kTopologies = {
            TopologyIndex::FD_FD,
            TopologyIndex::CD_FD,
            TopologyIndex::CD_FT
        };

        for (int ti = 0; ti < 3; ++ti) {
            RowCounts dst;
            WeightedRowCounts corrected_dst;
            for (int r = 0; r < (int)rows.size(); ++r) {
                if (nonzero(topo_dense[ti][r])) {
                    dst.emplace(r, topo_dense[ti][r]);
                }
                if (apply_current_weights) {
                    const WeightedHelCounts& wh = corrected_topo_dense[ti][r];
                    if (wh.sumw.unpol != 0.0 || wh.sumw.pos != 0.0 || wh.sumw.neg != 0.0) {
                        corrected_dst.emplace(r, wh);
                    }
                }
            }

            // Preserve the original behavior: a topology key exists only if at
            // least one event matched a CSV row for that topology.
            if (!dst.empty()) {
                out.topo_counts.emplace(topology_name(kTopologies[ti]), std::move(dst));
            }
            if (!corrected_dst.empty()) {
                out.corrected_topo_counts.emplace(topology_name(kTopologies[ti]), std::move(corrected_dst));
            }
        }
    }

    std::cout << "[total_counts] channel=" << work_cfg.channel_cfg.csv_channel
              << " sample=" << (is_gen ? "gen" : (is_data ? "data" : "rec"))
              << " tree=" << tags.tree_key
              << " entries=" << (long long)N
              << " global_pass=" << n_global_pass
              << " sig_pass=" << n_sigma_pass
              << " matched=" << n_used
              << " current_excluded=" << n_current_excluded
              << std::endl;

    if (!current_excluded_runs_seen.empty()) {
        std::ostringstream ss;
        bool first = true;
        for (int run : current_excluded_runs_seen) {
            if (!first) ss << ",";
            ss << run;
            first = false;
        }
        std::cout << "[total_counts] Skipped from current-corrected accumulator only: "
                  << n_current_excluded << " accepted event(s) from explicitly unusable/zero-charge run(s): "
                  << ss.str() << ". Raw yield-total bookkeeping is unchanged." << std::endl;
    }

    return out;
}

// -----------------------------------------------------------------------------
// Formatting and writing counts
// -----------------------------------------------------------------------------

static inline std::string fmt0(double v) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(0) << v;
    return oss.str();
}

static inline std::string fmt_count_triple(double n) {
    if (!(std::isfinite(n) && n >= 0.0)) {
        return "";
    }

    const double stat = (n > 0.0) ? std::sqrt(n) : 0.0;

    std::ostringstream oss;
    oss << std::setprecision(12)
        << "(" << n << "," << stat << ",0)";
    return oss.str();
}

static RowCounts sum_row_counts(const RowCounts& a, const RowCounts& b) {
    RowCounts out = a;

    for (const auto& kv : b) {
        const int row = kv.first;
        const HelCounts& h = kv.second;

        HelCounts& o = out[row];

        o.unpol += h.unpol;
        o.pos += h.pos;
        o.neg += h.neg;
    }

    return out;
}


static WeightedRowCounts sum_weighted_row_counts(const WeightedRowCounts& a,
                                                 const WeightedRowCounts& b) {
    WeightedRowCounts out = a;
    auto add_hel = [](HelCounts& dst, const HelCounts& src) {
        dst.unpol += src.unpol;
        dst.pos += src.pos;
        dst.neg += src.neg;
    };
    for (const auto& kv : b) {
        WeightedHelCounts& o = out[kv.first];
        const WeightedHelCounts& h = kv.second;
        add_hel(o.sumw, h.sumw);
        add_hel(o.sumw2, h.sumw2);
        for (int ir = 0; ir < kCurrentRegionCount; ++ir) {
            add_hel(o.derivative_by_region[ir], h.derivative_by_region[ir]);
            o.parameter_stat[ir] = std::max(o.parameter_stat[ir], h.parameter_stat[ir]);
        }
        add_hel(o.angular_derivative, h.angular_derivative);
        o.angular_parameter_stat =
            std::max(o.angular_parameter_stat, h.angular_parameter_stat);
    }
    return out;
}

static double weighted_stat_variance(const WeightedHelCounts& h,
                                     char component) {
    auto comp = [&](const HelCounts& x) -> double {
        if (component == 'p') return x.pos;
        if (component == 'n') return x.neg;
        return x.unpol;
    };

    // Event-counting statistics only.  The fitted current-response parameters
    // are calibration quantities shared by many events and many analysis bins;
    // their uncertainties are therefore correlated systematics, not additional
    // Poisson/statistical variance.
    return std::max(0.0, comp(h.sumw2));
}

static std::string fmt_weighted_triple(double value, double variance) {
    if (!(std::isfinite(value) && value >= 0.0 && std::isfinite(variance) && variance >= 0.0)) return "";
    std::ostringstream oss;
    oss << std::setprecision(12) << "(" << value << "," << std::sqrt(variance) << ",0)";
    return oss.str();
}


static double weighted_total_value(const WeightedHelCounts& h) {
    return h.sumw.unpol + h.sumw.pos + h.sumw.neg;
}

static double weighted_total_variance(const WeightedHelCounts& h) {
    // Event-counting statistics only.  Current-response fit uncertainties are
    // exported separately by write_current_nuisance_response_csv().
    return std::max(
        0.0,
        h.sumw2.unpol + h.sumw2.pos + h.sumw2.neg
    );
}

static void add_count_map(std::unordered_map<std::string, long long>& dst,
                          const std::unordered_map<std::string, long long>& src) {
    for (const auto& kv : src) {
        dst[kv.first] += kv.second;
    }
}

static void add_nested_count_map(
    std::unordered_map<std::string, std::unordered_map<std::string, long long>>& dst,
    const std::unordered_map<std::string, std::unordered_map<std::string, long long>>& src) {
    for (const auto& outer : src) {
        for (const auto& inner : outer.second) {
            dst[outer.first][inner.first] += inner.second;
        }
    }
}

static CutFlowSummary sum_cut_flow(const CutFlowSummary& a, const CutFlowSummary& b) {
    CutFlowSummary out = a;

    out.entries += b.entries;
    out.valid_topology += b.valid_topology;
    out.global_pass += b.global_pass;
    out.sigma_pass += b.sigma_pass;
    out.matched += b.matched;

    add_count_map(out.topology_entries, b.topology_entries);
    add_count_map(out.topology_global_pass, b.topology_global_pass);
    add_count_map(out.topology_sigma_pass, b.topology_sigma_pass);
    add_count_map(out.topology_matched, b.topology_matched);

    add_count_map(out.sigma_single_pass, b.sigma_single_pass);
    add_count_map(out.sigma_cumulative_pass, b.sigma_cumulative_pass);
    add_nested_count_map(out.topology_sigma_single_pass, b.topology_sigma_single_pass);
    add_nested_count_map(out.topology_sigma_cumulative_pass, b.topology_sigma_cumulative_pass);

    return out;
}

static double row_counts_total(const RowCounts& rc) {
    double total = 0.0;

    for (const auto& kv : rc) {
        const HelCounts& h = kv.second;
        total += h.unpol + h.pos + h.neg;
    }

    return total;
}

struct CountCollection {
    WorkConfig work_cfg;
    std::unordered_map<std::string, RowCounts> total_by_period;
    std::unordered_map<std::string, std::unordered_map<std::string, RowCounts>> topo_by_period;
    std::unordered_map<std::string, WeightedRowCounts> corrected_total_by_period;
    std::unordered_map<std::string, std::unordered_map<std::string, WeightedRowCounts>> corrected_topo_by_period;
    std::unordered_map<std::string, CutFlowSummary> flow_by_period;
};

static std::string collection_key(const WorkConfig& cfg) {
    std::ostringstream ss;
    ss << cfg.channel_cfg.csv_channel << "::";

    if (cfg.sample_kind == SampleKind::DATA) {
        ss << "data";
    } else if (cfg.sample_kind == SampleKind::MC_GEN) {
        ss << "gen";
    } else {
        ss << "rec";
    }

    return ss.str();
}


static const std::vector<std::string>& diagnostic_period_order() {
    static const std::vector<std::string> periods = {
        "Fa18 Inb",
        "Sp18 Inb",
        "Fa18 Out",
        "Sp18 Out",
        "Sp19 Inb"
    };

    return periods;
}

static const std::vector<std::string>& diagnostic_topology_order() {
    static const std::vector<std::string> topologies = {
        "FD_FD",
        "CD_FD",
        "CD_FT"
    };

    return topologies;
}

static double safe_ratio(double numerator, double denominator) {
    if (!(std::isfinite(numerator) && std::isfinite(denominator)) || denominator == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return numerator / denominator;
}

static std::string fmt_diag(double v, int precision = 6) {
    if (!std::isfinite(v)) {
        return "nan";
    }

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << v;
    return ss.str();
}

static long long flow_stage_value(const CutFlowSummary& f, const std::string& stage) {
    if (stage == "entries") return f.entries;
    if (stage == "topology") return f.valid_topology;
    if (stage == "global") return f.global_pass;
    if (stage == "sigma") return f.sigma_pass;
    if (stage == "matched") return f.matched;
    return 0;
}

static long long flow_topology_stage_value(const CutFlowSummary& f,
                                           const std::string& topo,
                                           const std::string& stage) {
    const std::unordered_map<std::string, long long>* m = nullptr;

    if (stage == "topology") {
        m = &f.topology_entries;
    } else if (stage == "global") {
        m = &f.topology_global_pass;
    } else if (stage == "sigma") {
        m = &f.topology_sigma_pass;
    } else if (stage == "matched") {
        m = &f.topology_matched;
    } else {
        return 0;
    }

    auto it = m->find(topo);
    return (it == m->end()) ? 0 : it->second;
}

static const CountCollection* find_collection_by_channel_and_kind(
    const std::map<std::string, CountCollection>& collections,
    const std::string& csv_channel,
    SampleKind kind) {

    for (const auto& kv : collections) {
        const CountCollection& C = kv.second;

        if (C.work_cfg.channel_cfg.csv_channel == csv_channel &&
            C.work_cfg.sample_kind == kind) {
            return &C;
        }
    }

    return nullptr;
}

static double period_row_count_total(const CountCollection* C,
                                     const std::string& period) {
    if (!C) return std::numeric_limits<double>::quiet_NaN();

    auto it = C->total_by_period.find(period);
    if (it == C->total_by_period.end()) return 0.0;

    return row_counts_total(it->second);
}

static double period_topology_row_count_total(const CountCollection* C,
                                              const std::string& period,
                                              const std::string& topo) {
    if (!C) return std::numeric_limits<double>::quiet_NaN();

    auto ip = C->topo_by_period.find(period);
    if (ip == C->topo_by_period.end()) return 0.0;

    auto it = ip->second.find(topo);
    if (it == ip->second.end()) return 0.0;

    return row_counts_total(it->second);
}

static const CutFlowSummary* period_flow(const CountCollection& C,
                                         const std::string& period) {
    auto it = C.flow_by_period.find(period);
    if (it == C.flow_by_period.end()) return nullptr;
    return &(it->second);
}

static void print_period_ratio_line(const CountCollection& recC,
                                    const CountCollection* genC,
                                    const std::string& numerator,
                                    const std::string& denominator) {
    const CutFlowSummary* fn = period_flow(recC, numerator);
    const CutFlowSummary* fd = period_flow(recC, denominator);

    if (!fn || !fd) {
        std::cout << "[total_counts][REC-MC-RATIO] channel="
                  << recC.work_cfg.channel_cfg.csv_channel
                  << " ratio=" << numerator << "/" << denominator
                  << " missing_period_flow" << std::endl;
        return;
    }

    const double gen_num = period_row_count_total(genC, numerator);
    const double gen_den = period_row_count_total(genC, denominator);
    const double rec_num = period_row_count_total(&recC, numerator);
    const double rec_den = period_row_count_total(&recC, denominator);
    const double acc_num = safe_ratio(rec_num, gen_num);
    const double acc_den = safe_ratio(rec_den, gen_den);

    std::cout << "[total_counts][REC-MC-RATIO] channel="
              << recC.work_cfg.channel_cfg.csv_channel
              << " ratio=" << numerator << "/" << denominator
              << " entries=" << fmt_diag(safe_ratio((double)fn->entries, (double)fd->entries))
              << " topology=" << fmt_diag(safe_ratio((double)fn->valid_topology, (double)fd->valid_topology))
              << " global=" << fmt_diag(safe_ratio((double)fn->global_pass, (double)fd->global_pass))
              << " sigma=" << fmt_diag(safe_ratio((double)fn->sigma_pass, (double)fd->sigma_pass))
              << " matched=" << fmt_diag(safe_ratio((double)fn->matched, (double)fd->matched))
              << " gen_matched=" << fmt_diag(safe_ratio(gen_num, gen_den))
              << " csv_acceptance_like_double_input=" << fmt_diag(safe_ratio(acc_num, acc_den))
              << " acc_num=" << fmt_diag(acc_num)
              << " acc_den=" << fmt_diag(acc_den)
              << std::endl;
}

static void print_topology_ratio_line(const CountCollection& recC,
                                      const CountCollection* genC,
                                      const std::string& numerator,
                                      const std::string& denominator,
                                      const std::string& topo) {
    const CutFlowSummary* fn = period_flow(recC, numerator);
    const CutFlowSummary* fd = period_flow(recC, denominator);

    if (!fn || !fd) return;

    const double gen_num = period_row_count_total(genC, numerator);
    const double gen_den = period_row_count_total(genC, denominator);
    const double rec_num = period_topology_row_count_total(&recC, numerator, topo);
    const double rec_den = period_topology_row_count_total(&recC, denominator, topo);
    const double acc_num = safe_ratio(rec_num, gen_num);
    const double acc_den = safe_ratio(rec_den, gen_den);

    std::cout << "[total_counts][REC-MC-TOPO-RATIO] channel="
              << recC.work_cfg.channel_cfg.csv_channel
              << " topo=" << topo
              << " ratio=" << numerator << "/" << denominator
              << " topology_entries=" << fmt_diag(safe_ratio(
                     (double)flow_topology_stage_value(*fn, topo, "topology"),
                     (double)flow_topology_stage_value(*fd, topo, "topology")))
              << " global=" << fmt_diag(safe_ratio(
                     (double)flow_topology_stage_value(*fn, topo, "global"),
                     (double)flow_topology_stage_value(*fd, topo, "global")))
              << " sigma=" << fmt_diag(safe_ratio(
                     (double)flow_topology_stage_value(*fn, topo, "sigma"),
                     (double)flow_topology_stage_value(*fd, topo, "sigma")))
              << " matched=" << fmt_diag(safe_ratio(
                     (double)flow_topology_stage_value(*fn, topo, "matched"),
                     (double)flow_topology_stage_value(*fd, topo, "matched")))
              << " rec_topo_over_gen_ratio=" << fmt_diag(safe_ratio(acc_num, acc_den))
              << " acc_num=" << fmt_diag(acc_num)
              << " acc_den=" << fmt_diag(acc_den)
              << std::endl;
}

static long long flow_sigma_var_value(const CutFlowSummary& f,
                                      const std::string& mode,
                                      const std::string& var) {
    const std::unordered_map<std::string, long long>* m = nullptr;

    if (mode == "single") {
        m = &f.sigma_single_pass;
    } else if (mode == "cumulative") {
        m = &f.sigma_cumulative_pass;
    } else {
        return 0;
    }

    auto it = m->find(var);
    return (it == m->end()) ? 0 : it->second;
}

static long long flow_topology_sigma_var_value(const CutFlowSummary& f,
                                               const std::string& topo,
                                               const std::string& mode,
                                               const std::string& var) {
    const std::unordered_map<std::string, std::unordered_map<std::string, long long>>* outer = nullptr;

    if (mode == "single") {
        outer = &f.topology_sigma_single_pass;
    } else if (mode == "cumulative") {
        outer = &f.topology_sigma_cumulative_pass;
    } else {
        return 0;
    }

    auto ito = outer->find(topo);

    if (ito == outer->end()) {
        return 0;
    }

    auto iti = ito->second.find(var);
    return (iti == ito->second.end()) ? 0 : iti->second;
}

static void print_sigma_variable_period_lines(const CountCollection& recC,
                                              const CutFlowSummary& f,
                                              const std::string& period,
                                              const std::string& topo) {
    const auto& vars = production_variable_order();

    const double denominator = (topo == "ALL")
        ? (double)f.global_pass
        : (double)flow_topology_stage_value(f, topo, "global");

    for (const std::string& var : vars) {
        const long long single_count = (topo == "ALL")
            ? flow_sigma_var_value(f, "single", var)
            : flow_topology_sigma_var_value(f, topo, "single", var);

        const long long cumulative_count = (topo == "ALL")
            ? flow_sigma_var_value(f, "cumulative", var)
            : flow_topology_sigma_var_value(f, topo, "cumulative", var);

        std::cout << "[total_counts][REC-MC-SIGMA-CUT] channel="
                  << recC.work_cfg.channel_cfg.csv_channel
                  << " period=" << period
                  << " topo=" << topo
                  << " var=" << var
                  << " single=" << single_count
                  << " cumulative=" << cumulative_count
                  << " single/global=" << fmt_diag(safe_ratio((double)single_count, denominator))
                  << " cumulative/global=" << fmt_diag(safe_ratio((double)cumulative_count, denominator))
                  << std::endl;
    }
}

static void print_sigma_variable_ratio_lines(const CountCollection& recC,
                                             const std::string& topo) {
    const CutFlowSummary* f_si = period_flow(recC, "Sp18 Inb");
    const CutFlowSummary* f_fi = period_flow(recC, "Fa18 Inb");
    const CutFlowSummary* f_so = period_flow(recC, "Sp18 Out");
    const CutFlowSummary* f_fo = period_flow(recC, "Fa18 Out");

    if (!(f_si && f_fi && f_so && f_fo)) {
        return;
    }

    const auto& vars = production_variable_order();

    auto denom = [&](const CutFlowSummary& f)->double {
        return (topo == "ALL")
            ? (double)f.global_pass
            : (double)flow_topology_stage_value(f, topo, "global");
    };

    auto val = [&](const CutFlowSummary& f, const std::string& mode, const std::string& var)->double {
        const long long n = (topo == "ALL")
            ? flow_sigma_var_value(f, mode, var)
            : flow_topology_sigma_var_value(f, topo, mode, var);

        return safe_ratio((double)n, denom(f));
    };

    for (const std::string& var : vars) {
        for (const std::string& mode : {std::string("single"), std::string("cumulative")}) {
            const double rinb = safe_ratio(val(*f_si, mode, var), val(*f_fi, mode, var));
            const double rout = safe_ratio(val(*f_so, mode, var), val(*f_fo, mode, var));

            std::cout << "[total_counts][REC-MC-SIGMA-RATIO] channel="
                      << recC.work_cfg.channel_cfg.csv_channel
                      << " topo=" << topo
                      << " var=" << var
                      << " mode=" << mode
                      << " Sp18Inb_over_Fa18Inb=" << fmt_diag(rinb)
                      << " Sp18Out_over_Fa18Out=" << fmt_diag(rout)
                      << " double_ratio=" << fmt_diag(safe_ratio(rinb, rout))
                      << std::endl;
        }
    }
}

static void write_reconstructed_mc_survival_csv(const std::map<std::string, CountCollection>& collections,
                                                const std::string& out_root_dir) {
    const std::string root = normalize_total_counts_root(out_root_dir);
    mkdir_p(root);

    const std::string path = root + "/reconstructed_mc_cutflow_diagnostics.csv";
    std::ofstream out(path);

    if (!out) {
        std::cerr << "[total_counts] WARNING: could not write diagnostic CSV: " << path << std::endl;
        return;
    }

    out << "channel,period,topology,stage,variable,mode,count,denominator,fraction,"
        << "Sp18Inb_over_Fa18Inb,Sp18Out_over_Fa18Out,double_ratio\n";

    auto write_row = [&](const std::string& channel,
                         const std::string& period,
                         const std::string& topo,
                         const std::string& stage,
                         const std::string& var,
                         const std::string& mode,
                         double count,
                         double denominator,
                         double rinb,
                         double rout) {
        out << channel << "," << period << "," << topo << ","
            << stage << "," << var << "," << mode << ","
            << fmt_diag(count) << "," << fmt_diag(denominator) << ","
            << fmt_diag(safe_ratio(count, denominator)) << ","
            << fmt_diag(rinb) << "," << fmt_diag(rout) << ","
            << fmt_diag(safe_ratio(rinb, rout)) << "\n";
    };

    for (const auto& kv : collections) {
        const CountCollection& recC = kv.second;

        if (recC.work_cfg.sample_kind != SampleKind::MC_REC) {
            continue;
        }

        for (const std::string& period : diagnostic_period_order()) {
            const CutFlowSummary* f = period_flow(recC, period);

            if (!f) {
                continue;
            }

            for (const std::string& topo : std::vector<std::string>{"ALL", "FD_FD", "CD_FD", "CD_FT"}) {
                const double den = (topo == "ALL")
                    ? (double)f->global_pass
                    : (double)flow_topology_stage_value(*f, topo, "global");

                const long long sigma_count = (topo == "ALL")
                    ? f->sigma_pass
                    : flow_topology_stage_value(*f, topo, "sigma");

                write_row(recC.work_cfg.channel_cfg.csv_channel, period, topo,
                          "sigma_all", "all", "all", sigma_count, den,
                          std::numeric_limits<double>::quiet_NaN(),
                          std::numeric_limits<double>::quiet_NaN());

                for (const std::string& var : production_variable_order()) {
                    const long long single_count = (topo == "ALL")
                        ? flow_sigma_var_value(*f, "single", var)
                        : flow_topology_sigma_var_value(*f, topo, "single", var);
                    const long long cumulative_count = (topo == "ALL")
                        ? flow_sigma_var_value(*f, "cumulative", var)
                        : flow_topology_sigma_var_value(*f, topo, "cumulative", var);

                    write_row(recC.work_cfg.channel_cfg.csv_channel, period, topo,
                              "sigma_var", var, "single", single_count, den,
                              std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());

                    write_row(recC.work_cfg.channel_cfg.csv_channel, period, topo,
                              "sigma_var", var, "cumulative", cumulative_count, den,
                              std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN());
                }
            }
        }

        for (const std::string& topo : std::vector<std::string>{"ALL", "FD_FD", "CD_FD", "CD_FT"}) {
            const CutFlowSummary* f_si = period_flow(recC, "Sp18 Inb");
            const CutFlowSummary* f_fi = period_flow(recC, "Fa18 Inb");
            const CutFlowSummary* f_so = period_flow(recC, "Sp18 Out");
            const CutFlowSummary* f_fo = period_flow(recC, "Fa18 Out");

            if (!(f_si && f_fi && f_so && f_fo)) {
                continue;
            }

            auto den = [&](const CutFlowSummary& f)->double {
                return (topo == "ALL")
                    ? (double)f.global_pass
                    : (double)flow_topology_stage_value(f, topo, "global");
            };

            auto frac = [&](const CutFlowSummary& f, const std::string& mode, const std::string& var)->double {
                const long long n = (mode == "all")
                    ? ((topo == "ALL") ? f.sigma_pass : flow_topology_stage_value(f, topo, "sigma"))
                    : ((topo == "ALL")
                        ? flow_sigma_var_value(f, mode, var)
                        : flow_topology_sigma_var_value(f, topo, mode, var));

                return safe_ratio((double)n, den(f));
            };

            for (const std::string& var : production_variable_order()) {
                for (const std::string& mode : {std::string("single"), std::string("cumulative")}) {
                    const double rinb = safe_ratio(frac(*f_si, mode, var), frac(*f_fi, mode, var));
                    const double rout = safe_ratio(frac(*f_so, mode, var), frac(*f_fo, mode, var));

                    write_row(recC.work_cfg.channel_cfg.csv_channel, "RATIO", topo,
                              "sigma_var", var, mode,
                              std::numeric_limits<double>::quiet_NaN(),
                              std::numeric_limits<double>::quiet_NaN(),
                              rinb, rout);
                }
            }

            const double rinb_all = safe_ratio(frac(*f_si, "all", "all"), frac(*f_fi, "all", "all"));
            const double rout_all = safe_ratio(frac(*f_so, "all", "all"), frac(*f_fo, "all", "all"));

            write_row(recC.work_cfg.channel_cfg.csv_channel, "RATIO", topo,
                      "sigma_all", "all", "all",
                      std::numeric_limits<double>::quiet_NaN(),
                      std::numeric_limits<double>::quiet_NaN(),
                      rinb_all, rout_all);
        }
    }

    std::cout << "[total_counts][REC-MC-SURVIVAL] Wrote diagnostic CSV: " << path << std::endl;
}

static void print_reconstructed_mc_survival_summary(
    const std::map<std::string, CountCollection>& collections) {

    std::cout << "\n[total_counts][REC-MC-SURVIVAL] =====================================================" << std::endl;
    std::cout << "[total_counts][REC-MC-SURVIVAL] Reconstructed-MC cut-flow diagnostic." << std::endl;
    std::cout << "[total_counts][REC-MC-SURVIVAL] Key stages: entries -> valid topology -> global cuts -> 3exclusivity cuts -> matched CSV bin." << std::endl;
    std::cout << "[total_counts][REC-MC-SURVIVAL] The acceptance-like number here is final reconstructed matched counts divided by generated matched counts." << std::endl;

    for (const auto& kv : collections) {
        const CountCollection& recC = kv.second;

        if (recC.work_cfg.sample_kind != SampleKind::MC_REC) {
            continue;
        }

        const CountCollection* genC = find_collection_by_channel_and_kind(
            collections,
            recC.work_cfg.channel_cfg.csv_channel,
            SampleKind::MC_GEN);

        std::cout << "\n[total_counts][REC-MC-SURVIVAL] channel="
                  << recC.work_cfg.channel_cfg.csv_channel << std::endl;

        for (const std::string& period : diagnostic_period_order()) {
            const CutFlowSummary* f = period_flow(recC, period);

            if (!f) {
                continue;
            }

            const double gen_matched = period_row_count_total(genC, period);
            const double rec_matched = period_row_count_total(&recC, period);
            const double acc_like = safe_ratio(rec_matched, gen_matched);

            std::cout << "[total_counts][REC-MC-SURVIVAL] period=" << std::setw(9) << period
                      << " entries=" << f->entries
                      << " topology=" << f->valid_topology
                      << " global=" << f->global_pass
                      << " sigma=" << f->sigma_pass
                      << " matched=" << f->matched
                      << " f_topology=" << fmt_diag(safe_ratio((double)f->valid_topology, (double)f->entries))
                      << " f_global/topology=" << fmt_diag(safe_ratio((double)f->global_pass, (double)f->valid_topology))
                      << " f_sigma/global=" << fmt_diag(safe_ratio((double)f->sigma_pass, (double)f->global_pass))
                      << " f_matched/sigma=" << fmt_diag(safe_ratio((double)f->matched, (double)f->sigma_pass))
                      << " gen_matched=" << fmt_diag(gen_matched, 0)
                      << " rec_matched=" << fmt_diag(rec_matched, 0)
                      << " rec/gen=" << fmt_diag(acc_like)
                      << std::endl;

            for (const std::string& topo : diagnostic_topology_order()) {
                const long long te = flow_topology_stage_value(*f, topo, "topology");
                const long long tg = flow_topology_stage_value(*f, topo, "global");
                const long long ts = flow_topology_stage_value(*f, topo, "sigma");
                const long long tm = flow_topology_stage_value(*f, topo, "matched");
                const double topo_rec_matched = period_topology_row_count_total(&recC, period, topo);
                const double topo_acc_like = safe_ratio(topo_rec_matched, gen_matched);

                std::cout << "[total_counts][REC-MC-SURVIVAL-TOPO] period=" << std::setw(9) << period
                          << " topo=" << topo
                          << " topology=" << te
                          << " global=" << tg
                          << " sigma=" << ts
                          << " matched=" << tm
                          << " f_global/topology=" << fmt_diag(safe_ratio((double)tg, (double)te))
                          << " f_sigma/global=" << fmt_diag(safe_ratio((double)ts, (double)tg))
                          << " f_matched/sigma=" << fmt_diag(safe_ratio((double)tm, (double)ts))
                          << " rec_topo_matched=" << fmt_diag(topo_rec_matched, 0)
                          << " rec_topo/gen=" << fmt_diag(topo_acc_like)
                          << std::endl;
            }
        }

        for (const std::string& period : diagnostic_period_order()) {
            const CutFlowSummary* f = period_flow(recC, period);

            if (!f) {
                continue;
            }

            print_sigma_variable_period_lines(recC, *f, period, "ALL");

            for (const std::string& topo : diagnostic_topology_order()) {
                print_sigma_variable_period_lines(recC, *f, period, topo);
            }
        }

        for (const std::string& topo : std::vector<std::string>{"ALL", "FD_FD", "CD_FD", "CD_FT"}) {
            print_sigma_variable_ratio_lines(recC, topo);
        }

        print_period_ratio_line(recC, genC, "Sp18 Inb", "Fa18 Inb");
        print_period_ratio_line(recC, genC, "Sp18 Out", "Fa18 Out");

        const CutFlowSummary* f_si = period_flow(recC, "Sp18 Inb");
        const CutFlowSummary* f_fi = period_flow(recC, "Fa18 Inb");
        const CutFlowSummary* f_so = period_flow(recC, "Sp18 Out");
        const CutFlowSummary* f_fo = period_flow(recC, "Fa18 Out");

        if (f_si && f_fi && f_so && f_fo) {
            const std::vector<std::string> stages = {"entries", "topology", "global", "sigma", "matched"};

            std::cout << "[total_counts][REC-MC-DOUBLE-RATIO] channel="
                      << recC.work_cfg.channel_cfg.csv_channel;

            for (const std::string& stage : stages) {
                const double rinb = safe_ratio((double)flow_stage_value(*f_si, stage),
                                               (double)flow_stage_value(*f_fi, stage));
                const double rout = safe_ratio((double)flow_stage_value(*f_so, stage),
                                               (double)flow_stage_value(*f_fo, stage));

                std::cout << " " << stage << "=" << fmt_diag(safe_ratio(rinb, rout));
            }

            const double acc_si = safe_ratio(period_row_count_total(&recC, "Sp18 Inb"), period_row_count_total(genC, "Sp18 Inb"));
            const double acc_fi = safe_ratio(period_row_count_total(&recC, "Fa18 Inb"), period_row_count_total(genC, "Fa18 Inb"));
            const double acc_so = safe_ratio(period_row_count_total(&recC, "Sp18 Out"), period_row_count_total(genC, "Sp18 Out"));
            const double acc_fo = safe_ratio(period_row_count_total(&recC, "Fa18 Out"), period_row_count_total(genC, "Fa18 Out"));

            std::cout << " acceptance_like=" << fmt_diag(safe_ratio(safe_ratio(acc_si, acc_fi), safe_ratio(acc_so, acc_fo)))
                      << std::endl;
        }

        for (const std::string& topo : diagnostic_topology_order()) {
            print_topology_ratio_line(recC, genC, "Sp18 Inb", "Fa18 Inb", topo);
            print_topology_ratio_line(recC, genC, "Sp18 Out", "Fa18 Out", topo);
        }
    }

    std::cout << "[total_counts][REC-MC-SURVIVAL] =====================================================\n" << std::endl;
}



static std::string csv_escape(const std::string& value) {
    if (value.find_first_of(",\"\n\\r") == std::string::npos) {
        return value;
    }
    std::string out = "\"";
    for (char c : value) {
        if (c == '"') out += "\"\"";
        else out.push_back(c);
    }
    out += "\"";
    return out;
}

static std::string current_nuisance_source_period(
    const CountCollection& C,
    const std::string& period) {

    // Sp19 DATA uses the Fa18 Inb response calibration.  Give the transferred
    // parameter the same nuisance key so its correlation with the Fa18
    // calibration is retained explicitly.
    if (C.work_cfg.sample_kind == SampleKind::DATA &&
        period == "Sp19 Inb") {
        return "Fa18 Inb";
    }

    return period;
}

static std::string current_nuisance_response_channel(
    const CountCollection& C,
    bool use_epg_mc_current_factor_for_eppi0_bkg) {

    if (C.work_cfg.sample_kind == SampleKind::MC_REC &&
        C.work_cfg.channel_cfg.channel == Channel::EPPI0_BKG_AS_DVCS &&
        use_epg_mc_current_factor_for_eppi0_bkg) {
        return "ep->epg";
    }

    return C.work_cfg.channel_cfg.csv_channel;
}

static std::string current_sample_token(const CountCollection& C) {
    const std::string channel = C.work_cfg.channel_cfg.csv_channel;

    if (C.work_cfg.sample_kind == SampleKind::DATA) {
        if (channel == "ep->epg") return "epg_data";
        if (channel == "ep->eppi0") return "eppi0_data";
        return channel + "_data";
    }

    if (C.work_cfg.sample_kind == SampleKind::MC_GEN) {
        if (channel == "ep->epg") return "epg_gen";
        if (channel == "ep->eppi0") return "eppi0_gen";
        return channel + "_gen";
    }

    if (channel == "ep->epg") return "epg_rec";
    if (channel == "ep->eppi0") return "eppi0_rec";
    if (channel == "ep->eppi0->epg") return "misid_rec";
    return channel + "_rec";
}

static void write_current_nuisance_response_csv(
    const std::map<std::string, CountCollection>& collections,
    const std::string& path,
    bool use_epg_mc_current_factor_for_eppi0_bkg) {

    const std::filesystem::path out_path(path);
    if (out_path.has_parent_path()) {
        std::filesystem::create_directories(out_path.parent_path());
    }

    std::ofstream out(path);
    if (!out.is_open()) {
        fatal("[total_counts] FATAL: cannot write current nuisance response CSV: " + path);
    }

    out << "row,period,sample,nuisance,nominal_yield,delta_yield_1sigma\n";
    out << std::setprecision(14);

    for (const auto& collection_kv : collections) {
        const CountCollection& C = collection_kv.second;
        const std::string sample = current_sample_token(C);

        // Generated MC has no current-response nuisance, but its nominal yield
        // is needed downstream to reconstruct S/A for the systematic study.
        if (C.work_cfg.sample_kind == SampleKind::MC_GEN) {
            for (const auto& per_kv : C.total_by_period) {
                const std::string& period = per_kv.first;
                if (should_skip_csv_for_label(period)) continue;

                for (const auto& row_kv : per_kv.second) {
                    const HelCounts& h = row_kv.second;
                    const double nominal = h.unpol + h.pos + h.neg;
                    out << row_kv.first << ','
                        << csv_escape(period) << ','
                        << csv_escape(sample) << ','
                        << csv_escape("NOMINAL") << ','
                        << nominal << ",0\n";
                }
            }
            continue;
        }

        for (const auto& per_kv : C.corrected_total_by_period) {
            const std::string& period = per_kv.first;
            if (should_skip_csv_for_label(period)) continue;

            const std::string source_period =
                current_nuisance_source_period(C, period);

            const std::string response_channel =
                current_nuisance_response_channel(
                    C,
                    use_epg_mc_current_factor_for_eppi0_bkg);

            const std::string sample_kind =
                C.work_cfg.sample_kind == SampleKind::DATA ? "data" : "mc";

            for (const auto& row_kv : per_kv.second) {
                const int row = row_kv.first;
                const WeightedHelCounts& h = row_kv.second;
                const double nominal = weighted_total_value(h);

                // Always write one nominal record.
                out << row << ','
                    << csv_escape(period) << ','
                    << csv_escape(sample) << ','
                    << csv_escape("NOMINAL") << ','
                    << nominal << ",0\n";

                for (int ir = 0; ir < kCurrentRegionCount; ++ir) {
                    const HelCounts& d = h.derivative_by_region[ir];
                    const double deriv = d.unpol + d.pos + d.neg;
                    const double sigma_parameter = h.parameter_stat[ir];

                    if (!(std::isfinite(deriv) &&
                          std::isfinite(sigma_parameter) &&
                          sigma_parameter > 0.0)) {
                        continue;
                    }

                    const double delta = deriv * sigma_parameter;
                    if (delta == 0.0) continue;

                    std::ostringstream nuisance;
                    nuisance << sample_kind << ':'
                             << response_channel << ':'
                             << source_period << ':'
                             << current_region_names()[ir];

                    out << row << ','
                        << csv_escape(period) << ','
                        << csv_escape(sample) << ','
                        << csv_escape(nuisance.str()) << ','
                        << nominal << ','
                        << delta << '\n';
                }

                const double angular_derivative =
                    h.angular_derivative.unpol +
                    h.angular_derivative.pos +
                    h.angular_derivative.neg;

                if (C.work_cfg.sample_kind == SampleKind::DATA &&
                    std::isfinite(angular_derivative) &&
                    std::isfinite(h.angular_parameter_stat) &&
                    h.angular_parameter_stat > 0.0 &&
                    angular_derivative != 0.0) {

                    const double delta =
                        angular_derivative * h.angular_parameter_stat;

                    std::ostringstream nuisance;
                    nuisance << "data:" << response_channel << ':'
                             << source_period << ":theta_e_gradient";

                    out << row << ','
                        << csv_escape(period) << ','
                        << csv_escape(sample) << ','
                        << csv_escape(nuisance.str()) << ','
                        << nominal << ','
                        << delta << '\n';
                }
            }
        }
    }

    std::cout << "[total_counts] Wrote correlated current-response nuisance "
              << "yield responses to " << path << std::endl;
}

static void write_collection_to_csv(CSV& csv,
                                    const CountCollection& C) {
    const WorkConfig& cfg = C.work_cfg;

    for (const auto& kvp : C.total_by_period) {
        const std::string& period_display = kvp.first;
        const RowCounts& total_counts = kvp.second;

        if (should_skip_csv_for_label(period_display)) {
            continue;
        }

        if (cfg.write_mc_generated_columns) {
            const int c = col_strict(csv, col_mc_generated(cfg.channel_cfg, period_display));

            for (const auto& row_kv : total_counts) {
                const int r = row_kv.first;
                const HelCounts& h = row_kv.second;

                if (r < 0 || r >= (int)csv.rows.size()) {
                    fatal("[total_counts] FATAL: row index out of range while writing generated MC.");
                }

                csv.rows[r][c] = fmt_count_triple(h.unpol + h.pos + h.neg);
            }
        }

        if (cfg.write_mc_reconstructed_columns) {
            const int c = col_strict(csv, col_mc_reconstructed_total(cfg.channel_cfg, period_display));

            for (const auto& row_kv : total_counts) {
                const int r = row_kv.first;
                const HelCounts& h = row_kv.second;

                if (r < 0 || r >= (int)csv.rows.size()) {
                    fatal("[total_counts] FATAL: row index out of range while writing reconstructed MC total.");
                }

                csv.rows[r][c] = fmt_count_triple(h.unpol + h.pos + h.neg);
            }
        }
    }

    for (const auto& kvp : C.topo_by_period) {
        const std::string& period_display = kvp.first;
        const auto& topo_map = kvp.second;

        if (should_skip_csv_for_label(period_display)) {
            continue;
        }

        for (const auto& kt : topo_map) {
            const std::string& topoDir = kt.first;
            const RowCounts& rc = kt.second;

            const std::string topoLabel = topo_label_for_csv(topoDir);

            if (topoLabel.empty()) {
                std::ostringstream ss;
                ss << "[total_counts] FATAL: cannot map topoDir '" << topoDir
                   << "' to CSV topology label.";
                fatal(ss.str());
            }

            if (cfg.write_data_raw_columns) {
                const int c_unpol = col_strict(csv, col_data_raw_counts(cfg.channel_cfg, topoLabel, period_display, "unpol"));
                const bool write_helicity_resolved = has_helicity_resolved_data_columns(period_display);
                const int c_pos = write_helicity_resolved
                                  ? col_strict(csv, col_data_raw_counts(cfg.channel_cfg, topoLabel, period_display, "pos"))
                                  : -1;
                const int c_neg = write_helicity_resolved
                                  ? col_strict(csv, col_data_raw_counts(cfg.channel_cfg, topoLabel, period_display, "neg"))
                                  : -1;

                for (const auto& row_kv : rc) {
                    const int r = row_kv.first;
                    const HelCounts& h = row_kv.second;

                    if (r < 0 || r >= (int)csv.rows.size()) {
                        fatal("[total_counts] FATAL: row index out of range while writing data raw counts.");
                    }

                    const double unpol = h.pos + h.neg + h.unpol;

                    csv.rows[r][c_unpol] = fmt_count_triple(unpol);

                    if (write_helicity_resolved) {
                        csv.rows[r][c_pos] = fmt_count_triple(h.pos);
                        csv.rows[r][c_neg] = fmt_count_triple(h.neg);
                    }
                }
            }

            if (cfg.write_mc_reconstructed_columns) {
                const int c_topo = col_strict(csv, col_mc_reconstructed_topo(cfg.channel_cfg, topoLabel, period_display));

                for (const auto& row_kv : rc) {
                    const int r = row_kv.first;
                    const HelCounts& h = row_kv.second;

                    if (r < 0 || r >= (int)csv.rows.size()) {
                        fatal("[total_counts] FATAL: row index out of range while writing reconstructed MC topology counts.");
                    }

                    csv.rows[r][c_topo] = fmt_count_triple(h.unpol + h.pos + h.neg);
                }
            }
        }
    }


    // Event-level current-corrected outputs. Raw/unit-weight columns above are
    // intentionally left unchanged so existing yield-total note material is
    // exactly reproducible.
    if (cfg.write_mc_reconstructed_columns) {
        for (const auto& kvp : C.corrected_total_by_period) {
            const std::string& period_display = kvp.first;
            if (should_skip_csv_for_label(period_display)) continue;
            const int c = col_strict(csv, col_mc_current_corrected_total(cfg.channel_cfg, period_display));
            for (const auto& row_kv : kvp.second) {
                const int r = row_kv.first;
                const WeightedHelCounts& h = row_kv.second;
                if (r < 0 || r >= (int)csv.rows.size()) fatal("[total_counts] FATAL: corrected MC row index out of range.");
                csv.rows[r][c] = fmt_weighted_triple(weighted_total_value(h), weighted_total_variance(h));
            }
        }
    }

    for (const auto& kvp : C.corrected_topo_by_period) {
        const std::string& period_display = kvp.first;
        if (should_skip_csv_for_label(period_display)) continue;
        const bool write_helicity_resolved = has_helicity_resolved_data_columns(period_display);

        for (const auto& kt : kvp.second) {
            const std::string topoLabel = topo_label_for_csv(kt.first);
            if (topoLabel.empty()) fatal("[total_counts] FATAL: cannot map corrected topology to CSV label.");

            if (cfg.write_data_raw_columns) {
                const int c_unpol = col_strict(csv, col_data_normalized_counts(cfg.channel_cfg, topoLabel, period_display, "unpol"));
                const int c_pos = write_helicity_resolved ? col_strict(csv, col_data_normalized_counts(cfg.channel_cfg, topoLabel, period_display, "pos")) : -1;
                const int c_neg = write_helicity_resolved ? col_strict(csv, col_data_normalized_counts(cfg.channel_cfg, topoLabel, period_display, "neg")) : -1;

                for (const auto& row_kv : kt.second) {
                    const int r = row_kv.first;
                    const WeightedHelCounts& h = row_kv.second;
                    if (r < 0 || r >= (int)csv.rows.size()) fatal("[total_counts] FATAL: corrected DATA row index out of range.");
                    csv.rows[r][c_unpol] = fmt_weighted_triple(weighted_total_value(h), weighted_total_variance(h));
                    if (write_helicity_resolved) {
                        csv.rows[r][c_pos] = fmt_weighted_triple(h.sumw.pos, weighted_stat_variance(h, 'p'));
                        csv.rows[r][c_neg] = fmt_weighted_triple(h.sumw.neg, weighted_stat_variance(h, 'n'));
                    }
                }
            }

            if (cfg.write_mc_reconstructed_columns) {
                const int c = col_strict(csv, col_mc_current_corrected_topo(cfg.channel_cfg, topoLabel, period_display));
                for (const auto& row_kv : kt.second) {
                    const int r = row_kv.first;
                    const WeightedHelCounts& h = row_kv.second;
                    if (r < 0 || r >= (int)csv.rows.size()) fatal("[total_counts] FATAL: corrected MC topology row index out of range.");
                    csv.rows[r][c] = fmt_weighted_triple(weighted_total_value(h), weighted_total_variance(h));
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
// Plotting data collections
// -----------------------------------------------------------------------------

struct RowPoint {
    double phi_x = 0.0;
    double pos = 0.0;
    double neg = 0.0;
    double pos_err = 0.0;
    double neg_err = 0.0;
    bool valid = false;
};

static inline double bin_center(double a, double b) {
    return 0.5 * (a + b);
}

static inline double poisson_err(double n) {
    return (n > 0.0) ? std::sqrt(n) : 0.0;
}

static inline bool cell_is_number(const std::string& s) {
    if (s.empty()) {
        return false;
    }

    char* e = nullptr;
    (void)std::strtod(s.c_str(), &e);

    return (e != s.c_str());
}

static inline double cell_to_double_or_nan(const std::string& s) {
    if (!cell_is_number(s)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::strtod(s.c_str(), nullptr);
}

static std::string col_phiavg_for_period(const std::string& period_display) {
    return std::string("phiavg, ") + period_display;
}

static std::string col_xBavg_for_period(const std::string& period_display) {
    return std::string("xBavg, ") + period_display;
}

static void plot_one_xB_canvas(const ChannelConfig& channel_cfg,
                               const std::string& outdir,
                               const std::string& label_for_title,
                               const std::string& topoDir,
                               const std::string& topoLabel,
                               const std::string& xB_text,
                               double xBmin_for_file,
                               const std::vector<RowBin>& rows,
                               const std::vector<int>& row_indices,
                               const std::vector<RowPoint>& points) {
    if (row_indices.empty()) {
        return;
    }

    struct Edge {
        double a;
        double b;
    };

    std::vector<Edge> Q2bins;
    std::vector<Edge> tbins;

    auto edge_eq = [](const Edge& u, const Edge& v) {
        return (u.a == v.a && u.b == v.b);
    };

    for (int ridx : row_indices) {
        const RowBin& r = rows[ridx];

        Edge q{r.Q2min, r.Q2max};
        Edge t{r.tmin,  r.tmax};

        if (std::find_if(Q2bins.begin(), Q2bins.end(),
                         [&](const Edge& e){ return edge_eq(e, q); }) == Q2bins.end()) {
            Q2bins.push_back(q);
        }

        if (std::find_if(tbins.begin(), tbins.end(),
                         [&](const Edge& e){ return edge_eq(e, t); }) == tbins.end()) {
            tbins.push_back(t);
        }
    }

    auto sort_edges = [](std::vector<Edge>& v) {
        std::sort(v.begin(), v.end(), [](const Edge& p, const Edge& q) {
            if (p.a != q.a) {
                return p.a < q.a;
            }

            return p.b < q.b;
        });
    };

    sort_edges(Q2bins);
    sort_edges(tbins);

    const int ncols = (int)Q2bins.size();
    const int nrows = (int)tbins.size();

    if (ncols <= 0 || nrows <= 0) {
        return;
    }

    const int W = 300 * ncols + 160;
    const int H = 260 * nrows + 240;

    TCanvas c("c_total_counts", "", W, H);

    TPad* top = new TPad("top", "top", 0.0, 0.90, 1.0, 1.0);
    TPad* grid = new TPad("grid", "grid", 0.0, 0.0, 1.0, 0.90);

    top->SetFillStyle(0);
    grid->SetFillStyle(0);

    top->Draw();
    grid->Draw();

    top->cd();

    TLatex t;
    t.SetNDC(true);
    t.SetTextFont(42);
    t.SetTextSize(0.45);

    {
        std::ostringstream ss;
        ss << channel_cfg.title_label << "   "
           << label_for_title << "   "
           << xB_text << "   "
           << topoLabel;
        t.DrawLatex(0.05, 0.35, ss.str().c_str());
    }

    grid->cd();
    grid->Divide(ncols, nrows, 0.0, 0.0);

    auto find_q = [&](double a, double b) {
        for (int i = 0; i < (int)Q2bins.size(); ++i) {
            if (Q2bins[i].a == a && Q2bins[i].b == b) {
                return i;
            }
        }

        return -1;
    };

    auto find_t = [&](double a, double b) {
        for (int i = 0; i < (int)tbins.size(); ++i) {
            if (tbins[i].a == a && tbins[i].b == b) {
                return i;
            }
        }

        return -1;
    };

    for (int it = 0; it < nrows; ++it) {
        for (int iq = 0; iq < ncols; ++iq) {
            const int pad_idx = it * ncols + iq + 1;
            grid->cd(pad_idx);

            gPad->SetLeftMargin(0.160);
            gPad->SetRightMargin(0.07);
            gPad->SetTopMargin(0.22);
            gPad->SetBottomMargin(0.18);
            gPad->SetGrid(1, 1);
            gPad->SetTickx(1);
            gPad->SetTicky(1);

            std::vector<RowPoint> cell;
            cell.reserve(row_indices.size());

            for (int ridx : row_indices) {
                const RowBin& r = rows[ridx];

                const int q = find_q(r.Q2min, r.Q2max);
                const int tt = find_t(r.tmin, r.tmax);

                if (q != iq || tt != it) {
                    continue;
                }

                const RowPoint& p = points[ridx];

                if (!p.valid) {
                    continue;
                }

                cell.push_back(p);
            }

            std::sort(cell.begin(), cell.end(), [](const RowPoint& a, const RowPoint& b) {
                return a.phi_x < b.phi_x;
            });

            const int N = (int)cell.size();

            TGraphErrors* gpos = new TGraphErrors();
            TGraphErrors* gneg = new TGraphErrors();

            gpos->SetMarkerStyle(24);
            gpos->SetMarkerColor(kRed);
            gpos->SetLineColor(kRed);
            gpos->SetLineWidth(1);

            gneg->SetMarkerStyle(20);
            gneg->SetMarkerColor(kBlue);
            gneg->SetLineColor(kBlue);
            gneg->SetLineWidth(1);

            double ymax = 0.0;

            for (int i = 0; i < N; ++i) {
                const double x = cell[i].phi_x;

                gpos->SetPoint(i, x, cell[i].pos);
                gpos->SetPointError(i, 0.0, cell[i].pos_err);

                gneg->SetPoint(i, x, cell[i].neg);
                gneg->SetPointError(i, 0.0, cell[i].neg_err);

                ymax = std::max(ymax,
                                std::max(cell[i].pos + cell[i].pos_err,
                                         cell[i].neg + cell[i].neg_err));
            }

            TH1F* frame = (TH1F*)gPad->DrawFrame(0.0, 0.0, 360.0, std::max(1.0, 1.15 * ymax));

            frame->SetTitle("");
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Counts");

            frame->GetXaxis()->CenterTitle(true);
            frame->GetYaxis()->CenterTitle(true);

            frame->GetXaxis()->SetNdivisions(505);

            frame->GetXaxis()->SetLabelSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.060);
            frame->GetXaxis()->SetTitleSize(0.070);
            frame->GetYaxis()->SetTitleSize(0.070);

            frame->GetXaxis()->SetTitleOffset(1.05);
            frame->GetYaxis()->SetTitleOffset(1.10);

            gpos->Draw("PE1 SAME");
            gneg->Draw("PE1 SAME");

            TLegend* leg = new TLegend(0.60, 0.73, 0.93, 0.92);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetBorderSize(1);
            leg->SetTextFont(42);
            leg->SetTextSize(0.055);
            leg->AddEntry(gpos, "+ helicity", "p");
            leg->AddEntry(gneg, "- helicity", "p");
            leg->Draw();

            TLatex txt;
            txt.SetNDC(true);
            txt.SetTextFont(42);
            txt.SetTextSize(0.060);

            const double Q2c = bin_center(Q2bins[iq].a, Q2bins[iq].b);
            const double tc  = bin_center(tbins[it].a,  tbins[it].b);

            {
                std::ostringstream ss;
                ss << "Q^{2}=" << std::fixed << std::setprecision(2) << Q2c
                   << "  |t|=" << std::fixed << std::setprecision(2) << tc;
                txt.DrawLatex(0.12, 0.83, ss.str().c_str());
            }
        }
    }

    mkdir_p(outdir);

    const int idx = (int)std::llround(xBmin_for_file * 1000.0);

    std::ostringstream fname;
    fname << outdir << "/plot_total_counts_";

    if (!channel_cfg.plot_file_token.empty()) {
        fname << channel_cfg.plot_file_token << "_";
    }

    fname << (is_combined_group_label(label_for_title)
              ? canonical_group_dir(label_for_title)
              : canonical_period_dir(label_for_title))
          << "_" << topoDir << "_xB_" << idx << ".png";

    c.SaveAs(fname.str().c_str());

    delete top;
    delete grid;
}

static void make_plots_for_label_and_topo(const ChannelConfig& channel_cfg,
                                          const std::string& label,
                                          const std::string& topoDir,
                                          const std::string& out_root_dir,
                                          const CSV& csv,
                                          const std::vector<RowBin>& rows,
                                          const RowCounts& row_counts_for_topo) {
    struct XBEdge {
        double a;
        double b;
    };

    std::vector<XBEdge> xbs;

    for (int r = 0; r < (int)rows.size(); ++r) {
        const RowBin& w = rows[r];

        if (!w.valid) {
            continue;
        }

        XBEdge e{w.xBmin, w.xBmax};

        auto it = std::find_if(xbs.begin(), xbs.end(), [&](const XBEdge& z) {
            return (z.a == e.a && z.b == e.b);
        });

        if (it == xbs.end()) {
            xbs.push_back(e);
        }
    }

    std::sort(xbs.begin(), xbs.end(), [](const XBEdge& p, const XBEdge& q) {
        if (p.a != q.a) {
            return p.a < q.a;
        }

        return p.b < q.b;
    });

    const bool use_phiavg = (!is_combined_group_label(label));

    int c_phiavg = -1;

    if (use_phiavg) {
        const std::string name = col_phiavg_for_period(label);
        auto it = csv.index.find(name);

        if (it != csv.index.end()) {
            c_phiavg = it->second;
        }
    }

    int c_xBavg = -1;

    if (!is_combined_group_label(label)) {
        const std::string name = col_xBavg_for_period(label);
        auto it = csv.index.find(name);

        if (it != csv.index.end()) {
            c_xBavg = it->second;
        }
    }

    const std::string topoLabel = topo_label_for_csv(topoDir);

    if (topoLabel.empty()) {
        std::ostringstream ss;
        ss << "[total_counts] FATAL: unknown topoDir '" << topoDir << "'";
        fatal(ss.str());
    }

    const std::string outdir = out_root_for_label(channel_cfg, label, out_root_dir) + "/" + topoDir;
    mkdir_p(outdir);

    for (const auto& xb : xbs) {
        std::vector<int> row_indices;
        row_indices.reserve(256);

        std::vector<RowPoint> points;
        points.resize(rows.size());

        for (int r = 0; r < (int)rows.size(); ++r) {
            const RowBin& w = rows[r];

            if (!w.valid) {
                continue;
            }

            if (!(w.xBmin == xb.a && w.xBmax == xb.b)) {
                continue;
            }

            row_indices.push_back(r);

            RowPoint p;

            if (c_phiavg >= 0) {
                const std::string& cell = csv.rows[r][c_phiavg];
                const double v = cell_to_double_or_nan(cell);

                if (std::isfinite(v)) {
                    p.phi_x = wrap_phi_deg(v);
                } else {
                    p.phi_x = wrap_phi_deg(bin_center(w.pmin, w.pmax));
                }
            } else {
                p.phi_x = wrap_phi_deg(bin_center(w.pmin, w.pmax));
            }

            auto itc = row_counts_for_topo.find(r);

            double pos = 0.0;
            double neg = 0.0;

            if (itc != row_counts_for_topo.end()) {
                pos = itc->second.pos;
                neg = itc->second.neg;
            }

            p.pos = pos;
            p.neg = neg;
            p.pos_err = poisson_err(pos);
            p.neg_err = poisson_err(neg);
            p.valid = true;

            points[r] = p;
        }

        std::ostringstream xbtxt;

        if (c_xBavg >= 0) {
            double xBavg = std::numeric_limits<double>::quiet_NaN();

            for (int r : row_indices) {
                const double v = cell_to_double_or_nan(csv.rows[r][c_xBavg]);

                if (std::isfinite(v)) {
                    xBavg = v;
                    break;
                }
            }

            if (std::isfinite(xBavg)) {
                xbtxt << "x_{B}=" << std::fixed << std::setprecision(3) << xBavg;
            } else {
                xbtxt << "x_{B} in [" << std::fixed << std::setprecision(2) << xb.a
                      << "," << std::fixed << std::setprecision(2) << xb.b << ")";
            }
        } else {
            xbtxt << "x_{B} in [" << std::fixed << std::setprecision(2) << xb.a
                  << "," << std::fixed << std::setprecision(2) << xb.b << ")";
        }

        plot_one_xB_canvas(channel_cfg,
                           outdir,
                           label,
                           topoDir,
                           topoLabel,
                           xbtxt.str(),
                           xb.a,
                           rows,
                           row_indices,
                           points);
    }
}

// -----------------------------------------------------------------------------
// Work item construction
// -----------------------------------------------------------------------------

struct WorkItem {
    WorkConfig work_cfg;
    PeriodTags tags;
    TTree* tree = nullptr;
};

static std::vector<WorkItem> build_work_items_for_map(const WorkConfig& work_cfg,
                                                      const std::map<std::string, TTree*>& trees) {
    std::vector<WorkItem> out;

    for (const auto& kv : trees) {
        if (!kv.second) {
            continue;
        }

        WorkItem w;
        w.work_cfg = work_cfg;
        w.tags = parse_period_tags_from_tree_key(kv.first);
        w.tree = kv.second;

        out.push_back(w);
    }

    return out;
}

static void append_items(std::vector<WorkItem>& out,
                         const std::vector<WorkItem>& add) {
    out.insert(out.end(), add.begin(), add.end());
}

static bool looks_like_nobkg_dvcs_mc_key(const std::string& key) {
    const std::string k = to_lower_ascii(key);

    // Preferred explicit tag.
    if (k.find("nobkg") != std::string::npos) {
        return true;
    }

    // Backward-compatible with the current load_trees.cpp convention: the
    // no-background dvcsgen files are loaded into the current-study maps and
    // tagged as 0nA rather than nobkg.
    if (k.find("0na") != std::string::npos) {
        return true;
    }

    return false;
}

static std::map<std::string, TTree*> filter_nobkg_tree_map(
    const std::map<std::string, TTree*>& in) {

    std::map<std::string, TTree*> out;

    for (const auto& kv : in) {
        if (!kv.second) {
            continue;
        }

        if (looks_like_nobkg_dvcs_mc_key(kv.first)) {
            out[kv.first] = kv.second;
        }
    }

    return out;
}


// -----------------------------------------------------------------------------
// Compact analysis-note raw-yield summaries
// -----------------------------------------------------------------------------

static double sum_hel_counts(const HelCounts& h) {
    return h.unpol + h.pos + h.neg;
}

static double sum_row_counts_all(const RowCounts& rc) {
    double s = 0.0;
    for (const auto& kv : rc) s += sum_hel_counts(kv.second);
    return s;
}

static int populated_row_count(const RowCounts& rc) {
    int n = 0;
    for (const auto& kv : rc) if (sum_hel_counts(kv.second) > 0.0) ++n;
    return n;
}

static const std::vector<std::string>& note_period_order() {
    static const std::vector<std::string> p = {
        "Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"
    };
    return p;
}

static const std::array<int,3>& note_topology_colors() {
    static const std::array<int,3> c = {{kBlue+1, kRed+1, kGreen+2}};
    return c;
}

static const std::array<int,5>& note_period_colors() {
    static const std::array<int,5> c = {{kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+7}};
    return c;
}

static std::string note_channel_token(const ChannelConfig& cfg) {
    if (cfg.channel == Channel::DVCS) return "dvcs";
    if (cfg.channel == Channel::EPPI0) return "eppi0";
    return "other";
}

static std::string note_channel_axis_label(const ChannelConfig& cfg) {
    if (cfg.channel == Channel::DVCS) return "e'p'#gamma candidate events";
    if (cfg.channel == Channel::EPPI0) return "e'p'#pi^{0} events";
    return "events";
}

static void write_note_yield_table(const CountCollection& C,
                                   const std::string& outdir) {
    gSystem->mkdir(outdir.c_str(), true);
    const std::string token = note_channel_token(C.work_cfg.channel_cfg);
    const std::string path = outdir + "/raw_yield_summary_" + token + ".csv";
    std::ofstream out(path);
    if (!out.is_open()) {
        std::cerr << "[total_counts] WARNING: cannot write analysis-note summary: " << path << std::endl;
        return;
    }
    out << "period,topology,raw_count,period_fraction_percent,populated_4d_bins\n";
    for (const std::string& period : note_period_order()) {
        auto ip = C.topo_by_period.find(period);
        if (ip == C.topo_by_period.end()) continue;
        double total = 0.0;
        for (const std::string& topo : diagnostic_topology_order()) {
            auto it = ip->second.find(topo);
            if (it != ip->second.end()) total += sum_row_counts_all(it->second);
        }
        for (const std::string& topo : diagnostic_topology_order()) {
            auto it = ip->second.find(topo);
            const double n = (it == ip->second.end()) ? 0.0 : sum_row_counts_all(it->second);
            const int populated = (it == ip->second.end()) ? 0 : populated_row_count(it->second);
            out << period << ',' << topo << ',' << std::fixed << std::setprecision(0) << n << ','
                << std::setprecision(3) << (total > 0.0 ? 100.0*n/total : 0.0) << ',' << populated << '\n';
        }
        out << period << ",ALL," << std::fixed << std::setprecision(0) << total << ",100.000,";
        auto ittot=C.total_by_period.find(period);
        out << (ittot==C.total_by_period.end()?0:populated_row_count(ittot->second)) << '\n';
    }
    std::cout << "[total_counts] Analysis-note raw-yield table written to: " << path << std::endl;
}

static void make_note_yield_summary_plot(const CountCollection& C,
                                         const std::vector<RowBin>& rows,
                                         const std::string& outdir) {
    const ChannelConfig& cfg = C.work_cfg.channel_cfg;
    gSystem->mkdir(outdir.c_str(), true);
    gStyle->SetOptStat(0);

    // Unique xB bins in production order.
    std::vector<AxisBin> xbins;
    for (const RowBin& r : rows) {
        if (!r.valid) continue;
        AxisBin b{r.xBmin,r.xBmax};
        auto it=std::find_if(xbins.begin(),xbins.end(),[&](const AxisBin& z){return z.min==b.min&&z.max==b.max;});
        if(it==xbins.end()) xbins.push_back(b);
    }
    std::sort(xbins.begin(),xbins.end(),[](const AxisBin&a,const AxisBin&b){return a.min<b.min;});

    TCanvas c(("c_note_yields_"+note_channel_token(cfg)).c_str(),"",1600,720);
    c.Divide(2,1,0.015,0.015);

    // Left: grouped topology totals for each period.
    c.cd(1);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.18); gPad->SetTopMargin(0.08);
    const int np = static_cast<int>(note_period_order().size());
    double ymax=0.0;
    std::array<std::vector<double>,3> topo_totals;
    for (int it=0;it<3;++it) topo_totals[it].assign(np,0.0);
    for(int iper=0;iper<np;++iper){
        auto ip=C.topo_by_period.find(note_period_order()[iper]);
        if(ip==C.topo_by_period.end()) continue;
        for(int it=0;it<3;++it){
            auto ir=ip->second.find(diagnostic_topology_order()[it]);
            if(ir!=ip->second.end()) topo_totals[it][iper]=sum_row_counts_all(ir->second);
            ymax=std::max(ymax,topo_totals[it][iper]);
        }
    }
    TH1F frame("frame_note_period_totals","",np,0.5,np+0.5);
    frame.SetMinimum(0.0); frame.SetMaximum(ymax>0?1.28*ymax:1.0);
    frame.GetYaxis()->SetTitle(("raw " + note_channel_axis_label(cfg)).c_str());
    frame.GetXaxis()->SetTitle("run period");
    for(int i=0;i<np;++i) frame.GetXaxis()->SetBinLabel(i+1,note_period_order()[i].c_str());
    frame.GetXaxis()->SetLabelSize(0.045); frame.GetYaxis()->SetLabelSize(0.042);
    frame.GetXaxis()->SetTitleSize(0.050); frame.GetYaxis()->SetTitleSize(0.048);
    frame.GetYaxis()->SetTitleOffset(1.45); frame.Draw("AXIS");
    TLegend leg1(0.54,0.73,0.93,0.91); leg1.SetBorderSize(0); leg1.SetFillStyle(0); leg1.SetTextSize(0.040);
    std::vector<TH1F*> bars;
    for(int it=0;it<3;++it){
        TH1F* h=new TH1F(("h_note_topo_"+std::to_string(it)).c_str(),"",np,0.5,np+0.5);
        h->SetDirectory(nullptr); h->SetFillColor(note_topology_colors()[it]); h->SetLineColor(note_topology_colors()[it]);
        h->SetBarWidth(0.24); h->SetBarOffset(0.12+0.27*it);
        for(int iper=0;iper<np;++iper) h->SetBinContent(iper+1,topo_totals[it][iper]);
        h->Draw("BAR SAME");
        std::string lbl=topo_label_for_csv(diagnostic_topology_order()[it]);
        leg1.AddEntry(h,lbl.c_str(),"f"); bars.push_back(h);
    }
    leg1.Draw();

    // Right: raw yield projected onto xB, summed over Q2, -t, phi, and topology.
    c.cd(2);
    gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.04); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.08);
    double yprojmax=0.0;
    std::array<std::vector<double>,5> xs,ys;
    for(int iper=0;iper<np;++iper){
        xs[iper].resize(xbins.size()); ys[iper].assign(xbins.size(),0.0);
        auto irc=C.total_by_period.find(note_period_order()[iper]);
        for(size_t ix=0;ix<xbins.size();++ix) xs[iper][ix]=0.5*(xbins[ix].min+xbins[ix].max);
        if(irc==C.total_by_period.end()) continue;
        for(const auto& kv:irc->second){
            const int r=kv.first; if(r<0||r>=static_cast<int>(rows.size())||!rows[r].valid) continue;
            for(size_t ix=0;ix<xbins.size();++ix){
                if(rows[r].xBmin==xbins[ix].min && rows[r].xBmax==xbins[ix].max){ ys[iper][ix]+=sum_hel_counts(kv.second); break; }
            }
        }
        for(double v:ys[iper]) yprojmax=std::max(yprojmax,v);
    }
    const double xmin=xbins.empty()?0.0:xbins.front().min;
    const double xmax=xbins.empty()?1.0:xbins.back().max;
    TH1F frame2("frame_note_xb_projection","",100,xmin,xmax);
    frame2.SetMinimum(0.0); frame2.SetMaximum(yprojmax>0?1.25*yprojmax:1.0);
    frame2.GetXaxis()->SetTitle("x_{B}"); frame2.GetYaxis()->SetTitle(("raw " + note_channel_axis_label(cfg)).c_str());
    frame2.GetXaxis()->SetTitleSize(0.050); frame2.GetYaxis()->SetTitleSize(0.048); frame2.GetYaxis()->SetTitleOffset(1.45);
    frame2.GetXaxis()->SetLabelSize(0.042); frame2.GetYaxis()->SetLabelSize(0.042); frame2.Draw("AXIS");
    TLegend leg2(0.58,0.67,0.93,0.91); leg2.SetBorderSize(0); leg2.SetFillStyle(0); leg2.SetTextSize(0.038);
    std::vector<TGraphErrors*> graphs;
    for(int iper=0;iper<np;++iper){
        if(xs[iper].empty()) continue;
        TGraphErrors* g=new TGraphErrors(static_cast<int>(xs[iper].size()));
        for(int ix=0;ix<static_cast<int>(xs[iper].size());++ix){ g->SetPoint(ix,xs[iper][ix],ys[iper][ix]); g->SetPointError(ix,0.0,std::sqrt(std::max(0.0,ys[iper][ix]))); }
        g->SetLineColor(note_period_colors()[iper]); g->SetMarkerColor(note_period_colors()[iper]); g->SetMarkerStyle(20+iper); g->SetMarkerSize(0.9); g->SetLineWidth(2);
        g->Draw("LP SAME"); leg2.AddEntry(g,note_period_order()[iper].c_str(),"lp"); graphs.push_back(g);
    }
    leg2.Draw();

    const std::string path=outdir+"/raw_yield_summary_"+note_channel_token(cfg)+".png";
    c.SaveAs(path.c_str());
    std::cout << "[total_counts] Analysis-note raw-yield plot written to: " << path << std::endl;
    for(TH1F* h:bars) delete h;
    for(TGraphErrors* g:graphs) delete g;
}

static void write_analysis_note_raw_yield_outputs(const std::map<std::string, CountCollection>& collections,
                                                  const std::vector<RowBin>& rows,
                                                  const std::string& out_root_dir) {
    const std::string outdir=normalize_total_counts_root(out_root_dir)+"/analysis_note";
    for(const auto& kv:collections){
        const CountCollection& C=kv.second;
        if(C.work_cfg.sample_kind!=SampleKind::DATA) continue;
        if(C.work_cfg.channel_cfg.channel!=Channel::DVCS && C.work_cfg.channel_cfg.channel!=Channel::EPPI0) continue;
        write_note_yield_table(C,outdir);
        make_note_yield_summary_plot(C,rows,outdir);
    }
}

// -----------------------------------------------------------------------------
// Public entry
// -----------------------------------------------------------------------------

} // namespace

bool update_total_counts_csv(const std::string& csv_path,
                             const std::map<std::string, TTree*>& dvcsDataTrees,
                             const std::map<std::string, TTree*>& eppi0DataTrees,
                             const std::map<std::string, TTree*>& dvcsGenMcTrees,
                             const std::map<std::string, TTree*>& dvcsRecMcTrees,
                             const std::map<std::string, TTree*>& eppi0GenMcTrees,
                             const std::map<std::string, TTree*>& eppi0RecMcTrees,
                             const std::map<std::string, TTree*>& eppi0BkgTrees,
                             const std::string& combined_cuts_json,
                             const std::string& out_root_dir,
                             int max_workers,
                             const TotalCountsOptions& options,
                             const std::map<std::string, TTree*>& dvcsNoBkgGenMcTrees,
                             const std::map<std::string, TTree*>& dvcsNoBkgRecMcTrees) {
    try {
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gStyle->SetOptStat(0);

        const bool trace_matches = env_flag("TOTAL_COUNTS_TRACE_MATCHES");

        CSV csv;
        load_csv(csv_path, csv);

        const std::vector<RowBin> rows = load_row_bins_from_csv(csv);
        const FastBinning fast_bins = build_fast_binning(rows);

        const TopoCutMap sigma_cuts_data = load_combined_cuts(combined_cuts_json, "data");
        const TopoCutMap sigma_cuts_mc   = load_combined_cuts(combined_cuts_json, "mc");

        CurrentResponseModel current_model;
        const CurrentResponseModel* current_model_ptr = nullptr;
        if (options.apply_event_level_current_correction) {
            current_model = load_current_response_model(options.current_response_model_json);
            validate_final_current_response_prescription(
                current_model,
                options.require_sp18_out_epg_e_theta_current_model);
            current_model_ptr = &current_model;
            std::cout << "[total_counts] Event-level regional current correction enabled from "
                      << options.current_response_model_json << std::endl;
        }

        WorkConfig dvcs_data;
        dvcs_data.channel_cfg = dvcs_config();
        dvcs_data.sample_kind = SampleKind::DATA;
        dvcs_data.write_data_raw_columns = true;
        dvcs_data.make_plots = options.make_plots;

        WorkConfig eppi0_data;
        eppi0_data.channel_cfg = eppi0_config();
        eppi0_data.sample_kind = SampleKind::DATA;
        eppi0_data.write_data_raw_columns = true;
        eppi0_data.make_plots = options.make_plots;

        WorkConfig dvcs_gen;
        dvcs_gen.channel_cfg = dvcs_config();
        dvcs_gen.sample_kind = SampleKind::MC_GEN;
        dvcs_gen.write_mc_generated_columns = true;

        WorkConfig dvcs_rec;
        dvcs_rec.channel_cfg = dvcs_config();
        dvcs_rec.sample_kind = SampleKind::MC_REC;
        dvcs_rec.write_mc_reconstructed_columns = true;

        WorkConfig eppi0_gen;
        eppi0_gen.channel_cfg = eppi0_config();
        eppi0_gen.sample_kind = SampleKind::MC_GEN;
        eppi0_gen.write_mc_generated_columns = true;

        WorkConfig eppi0_rec;
        eppi0_rec.channel_cfg = eppi0_config();
        eppi0_rec.sample_kind = SampleKind::MC_REC;
        eppi0_rec.write_mc_reconstructed_columns = true;

        WorkConfig eppi0_bkg_rec;
        eppi0_bkg_rec.channel_cfg = eppi0_bkg_as_dvcs_config();
        eppi0_bkg_rec.sample_kind = SampleKind::MC_REC;
        eppi0_bkg_rec.write_mc_reconstructed_columns = true;

        std::map<std::string, TTree*> dvcs_gen_for_counts = dvcsGenMcTrees;
        std::map<std::string, TTree*> dvcs_rec_for_counts = dvcsRecMcTrees;

        if (options.use_nobkg_dvcs_mc_counts) {
            const std::map<std::string, TTree*>& gen_source =
                dvcsNoBkgGenMcTrees.empty() ? dvcsGenMcTrees : dvcsNoBkgGenMcTrees;
            const std::map<std::string, TTree*>& rec_source =
                dvcsNoBkgRecMcTrees.empty() ? dvcsRecMcTrees : dvcsNoBkgRecMcTrees;

            dvcs_gen_for_counts = filter_nobkg_tree_map(gen_source);
            dvcs_rec_for_counts = filter_nobkg_tree_map(rec_source);

            if (dvcs_gen_for_counts.empty() || dvcs_rec_for_counts.empty()) {
                std::cerr << "[total_counts] FATAL: no-background DVCS MC override was enabled, "
                          << "but no no-background generated/reconstructed DVCS trees were found. "
                          << "Accepted key tags are explicit 'nobkg' or the existing current-study "
                          << "'0nA' convention." << std::endl;
                std::cerr << "[total_counts] Available generated DVCS override/source keys:" << std::endl;
                for (const auto& kv : gen_source) {
                    std::cerr << "  gen key: " << kv.first << std::endl;
                }
                std::cerr << "[total_counts] Available reconstructed DVCS override/source keys:" << std::endl;
                for (const auto& kv : rec_source) {
                    std::cerr << "  rec key: " << kv.first << std::endl;
                }
                fatal("[total_counts] FATAL: no-background DVCS MC override has no usable input trees.");
            }

            std::cout << "[total_counts] No-background DVCS MC override enabled: "
                      << "using " << dvcs_gen_for_counts.size()
                      << " generated and " << dvcs_rec_for_counts.size()
                      << " reconstructed no-background ep->epg MC tree(s) for DVCS MC counts "
                      << "(accepted tags: 'nobkg' or '0nA'). "
                      << "Data, ep->eppi0 MC, and ep->eppi0->epg background MC are unchanged."
                      << std::endl;
        }

        std::vector<WorkItem> work_items;

        append_items(work_items, build_work_items_for_map(dvcs_data, dvcsDataTrees));
        append_items(work_items, build_work_items_for_map(eppi0_data, eppi0DataTrees));
        append_items(work_items, build_work_items_for_map(dvcs_gen, dvcs_gen_for_counts));
        append_items(work_items, build_work_items_for_map(dvcs_rec, dvcs_rec_for_counts));
        append_items(work_items, build_work_items_for_map(eppi0_gen, eppi0GenMcTrees));
        append_items(work_items, build_work_items_for_map(eppi0_rec, eppi0RecMcTrees));
        append_items(work_items, build_work_items_for_map(eppi0_bkg_rec, eppi0BkgTrees));

        if (work_items.empty()) {
            fatal("[total_counts] FATAL: no trees available for total_counts.");
        }

        std::cout << "[total_counts] Will process " << work_items.size()
                  << " tree work item(s)." << std::endl;

        std::map<std::string, CountCollection> collections;

        auto ensure_collection = [&](const WorkConfig& cfg) {
            const std::string key = collection_key(cfg);

            if (collections.find(key) == collections.end()) {
                CountCollection C;
                C.work_cfg = cfg;
                collections[key] = C;
            }
        };

        for (const auto& item : work_items) {
            ensure_collection(item.work_cfg);
        }

        std::mutex merge_mutex;

        int nth = std::max(1, std::min(7, max_workers));
        nth = std::min(nth, (int)work_items.size());

        std::cout << "[total_counts] Using " << nth
                  << " worker thread(s), capped at 7." << std::endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
        for (int i = 0; i < (int)work_items.size(); ++i) {
            const WorkItem& w = work_items[i];

            const bool use_data_cuts = (w.work_cfg.sample_kind == SampleKind::DATA);
            const TopoCutMap& cuts = use_data_cuts ? sigma_cuts_data : sigma_cuts_mc;

            const WorkCounts counts =
                accumulate_counts_for_tree(w.work_cfg,
                                           w.tags,
                                           w.tree,
                                           rows,
                                           fast_bins,
                                           cuts,
                                           trace_matches,
                                           current_model_ptr,
                                           options.use_epg_mc_current_factor_for_eppi0_bkg);

            std::lock_guard<std::mutex> lock(merge_mutex);

            CountCollection& C = collections[collection_key(w.work_cfg)];

            C.total_by_period[w.tags.period_display] =
                sum_row_counts(C.total_by_period[w.tags.period_display],
                               counts.total_counts);

            if (!counts.corrected_total_counts.empty()) {
                C.corrected_total_by_period[w.tags.period_display] =
                    sum_weighted_row_counts(C.corrected_total_by_period[w.tags.period_display],
                                            counts.corrected_total_counts);
            }

            C.flow_by_period[w.tags.period_display] =
                sum_cut_flow(C.flow_by_period[w.tags.period_display],
                             counts.flow);

            for (const auto& kv : counts.topo_counts) {
                const std::string& topoDir = kv.first;
                const RowCounts& rc = kv.second;

                C.topo_by_period[w.tags.period_display][topoDir] =
                    sum_row_counts(C.topo_by_period[w.tags.period_display][topoDir],
                                   rc);
            }

            for (const auto& kv : counts.corrected_topo_counts) {
                const std::string& topoDir = kv.first;
                const WeightedRowCounts& rc = kv.second;
                C.corrected_topo_by_period[w.tags.period_display][topoDir] =
                    sum_weighted_row_counts(C.corrected_topo_by_period[w.tags.period_display][topoDir], rc);
            }
        }

        print_reconstructed_mc_survival_summary(collections);
        write_reconstructed_mc_survival_csv(collections, out_root_dir);

        if (options.write_current_nuisance_responses) {
            write_current_nuisance_response_csv(
                collections,
                options.current_nuisance_response_csv,
                options.use_epg_mc_current_factor_for_eppi0_bkg);
        }

        for (const auto& kv : collections) {
            write_collection_to_csv(csv, kv.second);
        }

        write_csv_atomic(csv_path, csv);

        std::cout << "[total_counts] Updated data and MC count columns in: "
                  << csv_path << std::endl;

        if (options.make_note_outputs) {
            write_analysis_note_raw_yield_outputs(collections, rows, out_root_dir);
        }

        for (const auto& kv : collections) {
            const CountCollection& C = kv.second;

            if (!C.work_cfg.make_plots) {
                continue;
            }

            for (const auto& kvp : C.topo_by_period) {
                const std::string& period_display = kvp.first;
                const auto& topoMap = kvp.second;

                for (const auto& kt : topoMap) {
                    make_plots_for_label_and_topo(C.work_cfg.channel_cfg,
                                                  period_display,
                                                  kt.first,
                                                  out_root_dir,
                                                  csv,
                                                  rows,
                                                  kt.second);
                }
            }
        }

        if (options.make_plots) {
            std::cout << "[total_counts] Plots written under: "
                      << normalize_total_counts_root(out_root_dir) << std::endl;
        } else {
            std::cout << "[total_counts] Per-bin total-count plots disabled for this run."
                      << std::endl;
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}

// =============================================================================
// Pass-2 acceptance model-dependence study by iterative reweighting
// =============================================================================

namespace {

struct ARWEvent {
    double x = 0.0;
    double q2 = 0.0;
    double tabs = 0.0;
    double phi = 0.0;
    double base_weight = 1.0;
    int row = -1;
    int subcell = -1;
};

struct ARWFineAxis {
    std::vector<double> edges;
    bool periodic = false;

    int bin(double value) const {
        if (edges.size() < 2 || !std::isfinite(value)) return -1;
        double v = value;
        if (periodic) {
            v = std::fmod(v, 360.0);
            if (v < 0.0) v += 360.0;
        }
        if (v < edges.front() || v > edges.back()) return -1;
        if (v == edges.back()) return (int)edges.size()-2;
        auto it = std::upper_bound(edges.begin(),edges.end(),v);
        const int ib = (int)std::distance(edges.begin(),it)-1;
        return (ib>=0 && ib<(int)edges.size()-1) ? ib : -1;
    }

    double center(int ib) const {
        return 0.5*(edges[(size_t)ib]+edges[(size_t)ib+1]);
    }

    int nbins() const { return std::max(0,(int)edges.size()-1); }
};

struct ARWWeightModel {
    std::array<ARWFineAxis,4> axes;
    std::array<std::vector<double>,4> factors;

    double value(const ARWEvent& e) const {
        const std::array<double,4> v={{e.x,e.q2,e.tabs,e.phi}};
        double w=1.0;
        for(int iv=0;iv<4;++iv){
            const int ib=axes[(size_t)iv].bin(v[(size_t)iv]);
            if(ib<0 || ib>=(int)factors[(size_t)iv].size()) continue;
            w*=factors[(size_t)iv][(size_t)ib];
        }
        return w;
    }
};

static double arw_tuple_first(const std::string& raw) {
    std::string s=raw;
    s.erase(std::remove_if(s.begin(),s.end(),
                           [](unsigned char c){return std::isspace(c);}),
            s.end());
    if(s.empty()) return std::numeric_limits<double>::quiet_NaN();
    if(s.front()=='('){
        const size_t comma=s.find(',');
        const std::string token=s.substr(1,comma==std::string::npos
                                            ? std::string::npos
                                            : comma-1);
        char* end=nullptr;
        const double v=std::strtod(token.c_str(),&end);
        return end==token.c_str()
            ? std::numeric_limits<double>::quiet_NaN() : v;
    }
    char* end=nullptr;
    const double v=std::strtod(s.c_str(),&end);
    return end==s.c_str()
        ? std::numeric_limits<double>::quiet_NaN() : v;
}

static int arw_find_row(const BranchBinder& b,
                        const std::vector<RowBin>& rows,
                        const FastBinning& fast_bins) {
    const int ix=find_axis_bin_index(fast_bins.xbins,b.x);
    const int iq=find_axis_bin_index(fast_bins.qbins,b.Q2);
    const int it=find_axis_bin_index(fast_bins.tbins,b.t_abs());
    if(ix<0||iq<0||it<0) return -1;
    const double phi=b.phi_deg();
    for(int r:fast_bins.rows_by_xqt[ix][iq][it]){
        if(!rows[(size_t)r].valid) continue;
        if(row_accepts_phi(phi,rows[(size_t)r].pmin,rows[(size_t)r].pmax))
            return r;
    }
    return -1;
}

static int arw_subcell(const ARWEvent& e,const RowBin& r) {
    const double xm=0.5*(r.xBmin+r.xBmax);
    const double qm=0.5*(r.Q2min+r.Q2max);
    const double tm=0.5*(r.tmin+r.tmax);
    const double pm=0.5*(r.pmin+r.pmax);
    const int ix=e.x>=xm;
    const int iq=e.q2>=qm;
    const int it=e.tabs>=tm;
    const int ip=e.phi>=pm;
    return (((ix*2)+iq)*2+it)*2+ip;
}

static std::array<ARWFineAxis,4> arw_build_axes(
    const std::vector<RowBin>& rows) {

    std::array<std::vector<double>,4> raw;
    for(const auto& r:rows){
        if(!r.valid) continue;
        const std::array<std::pair<double,double>,4> ranges={{
            {r.xBmin,r.xBmax},{r.Q2min,r.Q2max},
            {r.tmin,r.tmax},{r.pmin,r.pmax}
        }};
        for(int iv=0;iv<4;++iv){
            const double a=ranges[(size_t)iv].first;
            const double b=ranges[(size_t)iv].second;
            raw[(size_t)iv].push_back(a);
            raw[(size_t)iv].push_back(0.5*(a+b));
            raw[(size_t)iv].push_back(b);
        }
    }

    std::array<ARWFineAxis,4> out;
    for(int iv=0;iv<4;++iv){
        auto& v=raw[(size_t)iv];
        std::sort(v.begin(),v.end());
        std::vector<double> u;
        for(double x:v){
            if(u.empty() || std::fabs(x-u.back())>1e-10) u.push_back(x);
        }
        out[(size_t)iv].edges=u;
        out[(size_t)iv].periodic=(iv==3);
    }
    return out;
}

static TTree* arw_tree_for_period(
    const std::map<std::string,TTree*>& trees,
    const std::string& period) {
    for(const auto& kv:trees){
        if(!kv.second) continue;
        try{
            if(parse_period_tags_from_tree_key(kv.first).period_display==period)
                return kv.second;
        }catch(...){}
    }
    return nullptr;
}

static std::vector<ARWEvent> arw_collect(
    TTree* tree,
    const WorkConfig& work_cfg,
    const PeriodTags& tags,
    const std::vector<RowBin>& rows,
    const FastBinning& fast_bins,
    const TopoCutMap& sigma_cuts,
    const CurrentResponseModel* current_model,
    const CSV& csv,
    const AcceptanceReweightingOptions& options,
    bool is_generated) {

    std::vector<ARWEvent> out;
    if(!tree) return out;

    BranchBinder b;
    b.bind(tree,work_cfg);

    std::array<std::string,3> sigma_keys;
    std::array<CompiledSigmaPlan,3> sigma_plans;
    std::array<bool,3> sigma_ready{{false,false,false}};
    if(!is_generated){
        sigma_keys[0]=combined_cuts_key(
            work_cfg.channel_cfg,tags,topology_name(TopologyIndex::FD_FD));
        sigma_keys[1]=combined_cuts_key(
            work_cfg.channel_cfg,tags,topology_name(TopologyIndex::CD_FD));
        sigma_keys[2]=combined_cuts_key(
            work_cfg.channel_cfg,tags,topology_name(TopologyIndex::CD_FT));
    }

    DenseSigmaDiagnostics dummy_diag;

    int c_cont=-1;
    if(work_cfg.sample_kind==SampleKind::DATA &&
       options.apply_pi0_signal_fraction_to_data){
        const std::string col="contamination ratio, "+tags.period_display;
        auto it=csv.index.find(col);
        if(it!=csv.index.end()) c_cont=it->second;
    }

    const Long64_t N=tree->GetEntries();
    out.reserve((size_t)std::min<Long64_t>(N,3000000));

    for(Long64_t i=0;i<N;++i){
        tree->GetEntry(i);

        if(!is_generated){
            const TopologyIndex topo=topology_index(b.detector1,b.detector2);
            if(topo==TopologyIndex::INVALID) continue;
            const int ti=(int)topo;

            if(!passes_global_cuts_dispatch(b,tags.period_label)) continue;

            if(!sigma_ready[(size_t)ti]){
                sigma_plans[(size_t)ti]=compile_sigma_plan(
                    work_cfg.channel_cfg,sigma_cuts,sigma_keys[(size_t)ti]);
                sigma_ready[(size_t)ti]=true;
            }
            if(!fill_sigma_cut_diagnostics_compiled(
                    sigma_plans[(size_t)ti],b,ti,dummy_diag)) continue;
        }

        const int row=arw_find_row(b,rows,fast_bins);
        if(row<0) continue;

        ARWEvent e;
        e.x=b.x;
        e.q2=b.Q2;
        e.tabs=b.t_abs();
        e.phi=b.phi_deg();
        e.row=row;
        e.base_weight=1.0;

        if(!is_generated && options.apply_current_correction &&
           current_model!=nullptr){
            const EventCurrentWeight cw=event_current_weight(
                *current_model,work_cfg,tags,b,true);
            if(cw.skip_corrected) continue;
            e.base_weight*=cw.weight;
        }

        if(work_cfg.sample_kind==SampleKind::DATA && c_cont>=0){
            const double c=arw_tuple_first(
                csv.rows[(size_t)row][(size_t)c_cont]);
            if(std::isfinite(c))
                e.base_weight*=std::max(0.0,std::min(1.0,1.0-c));
        }

        e.subcell=arw_subcell(e,rows[(size_t)row]);
        out.push_back(e);
    }

    return out;
}

static std::vector<double> arw_hist(
    const std::vector<ARWEvent>& ev,
    const ARWWeightModel& model,
    int iv,
    bool apply_model) {

    const int n=model.axes[(size_t)iv].nbins();
    std::vector<double> h((size_t)n,0.0);
    for(const auto& e:ev){
        const std::array<double,4> v={{e.x,e.q2,e.tabs,e.phi}};
        const int ib=model.axes[(size_t)iv].bin(v[(size_t)iv]);
        if(ib<0) continue;
        double w=e.base_weight;
        if(apply_model) w*=model.value(e);
        h[(size_t)ib]+=w;
    }
    double sum=0.0;
    for(double x:h) sum+=x;
    if(sum>0.0) for(double& x:h) x/=sum;
    return h;
}

static double arw_shape_distance(const std::vector<double>& a,
                                 const std::vector<double>& b) {
    if(a.size()!=b.size()) return 1.0;
    double d=0.0;
    for(size_t i=0;i<a.size();++i) d+=std::fabs(a[i]-b[i]);
    return 0.5*d;
}

static std::vector<double> arw_smoothed_ratio(
    const std::vector<double>& data,
    const std::vector<double>& mc,
    bool periodic,
    int min_entries_proxy,
    const std::vector<double>& data_raw,
    const std::vector<double>& mc_raw,
    const AcceptanceReweightingOptions& options) {

    const int n=(int)data.size();
    std::vector<double> r((size_t)n,1.0);

    for(int i=0;i<n;++i){
        if(i>=(int)data_raw.size()||i>=(int)mc_raw.size()) continue;
        if(data_raw[(size_t)i]<min_entries_proxy ||
           mc_raw[(size_t)i]<min_entries_proxy) continue;
        if(mc[(size_t)i]<=0.0 || data[(size_t)i]<=0.0) continue;
        r[(size_t)i]=data[(size_t)i]/mc[(size_t)i];
    }

    std::vector<double> s=r;
    for(int i=0;i<n;++i){
        int im=i-1,ip=i+1;
        if(periodic){
            if(im<0) im=n-1;
            if(ip>=n) ip=0;
        }
        double num=2.0*r[(size_t)i],den=2.0;
        if(im>=0){num+=r[(size_t)im];den+=1.0;}
        if(ip<n){num+=r[(size_t)ip];den+=1.0;}
        s[(size_t)i]=num/den;
    }

    for(double& x:s){
        x=std::pow(std::max(1e-6,x),options.damping_power);
        x=std::max(options.per_iteration_weight_min,
                   std::min(options.per_iteration_weight_max,x));
    }
    return s;
}

static std::vector<double> arw_raw_counts(
    const std::vector<ARWEvent>& ev,
    const ARWWeightModel& model,
    int iv,
    bool apply_model) {
    const int n=model.axes[(size_t)iv].nbins();
    std::vector<double> h((size_t)n,0.0);
    for(const auto& e:ev){
        const std::array<double,4> v={{e.x,e.q2,e.tabs,e.phi}};
        const int ib=model.axes[(size_t)iv].bin(v[(size_t)iv]);
        if(ib<0) continue;
        double w=e.base_weight;
        if(apply_model) w*=model.value(e);
        h[(size_t)ib]+=w;
    }
    return h;
}

struct ARWFitResult {
    ARWWeightModel model;
    std::array<double,4> before_distance{{0,0,0,0}};
    std::array<double,4> after_distance{{0,0,0,0}};
    int iterations=0;
};

static ARWFitResult arw_fit(
    const std::vector<ARWEvent>& data,
    const std::vector<ARWEvent>& rec,
    const std::array<ARWFineAxis,4>& axes,
    const AcceptanceReweightingOptions& options) {

    ARWFitResult result;
    result.model.axes=axes;
    for(int iv=0;iv<4;++iv)
        result.model.factors[(size_t)iv].assign(
            (size_t)axes[(size_t)iv].nbins(),1.0);

    for(int iv=0;iv<4;++iv){
        result.before_distance[(size_t)iv]=arw_shape_distance(
            arw_hist(data,result.model,iv,false),
            arw_hist(rec,result.model,iv,false));
    }

    double previous=std::numeric_limits<double>::infinity();

    for(int iter=0;iter<options.max_iterations;++iter){
        for(int iv=0;iv<4;++iv){
            const auto dh=arw_hist(data,result.model,iv,false);
            const auto mh=arw_hist(rec,result.model,iv,true);
            const auto dr=arw_raw_counts(data,result.model,iv,false);
            const auto mr=arw_raw_counts(rec,result.model,iv,true);
            const auto update=arw_smoothed_ratio(
                dh,mh,axes[(size_t)iv].periodic,
                options.minimum_entries_per_fine_bin,dr,mr,options);

            auto& f=result.model.factors[(size_t)iv];
            for(size_t ib=0;ib<f.size();++ib){
                f[ib]*=update[ib];
                f[ib]=std::max(options.cumulative_weight_min,
                               std::min(options.cumulative_weight_max,f[ib]));
            }
        }

        double maxd=0.0;
        for(int iv=0;iv<4;++iv){
            const double d=arw_shape_distance(
                arw_hist(data,result.model,iv,false),
                arw_hist(rec,result.model,iv,true));
            maxd=std::max(maxd,d);
        }

        result.iterations=iter+1;
        if(maxd<options.convergence_shape_distance ||
           previous-maxd<options.convergence_improvement)
            break;
        previous=maxd;
    }

    for(int iv=0;iv<4;++iv){
        result.after_distance[(size_t)iv]=arw_shape_distance(
            arw_hist(data,result.model,iv,false),
            arw_hist(rec,result.model,iv,true));
    }
    return result;
}

static std::vector<double> arw_acceptance(
    const std::vector<ARWEvent>& gen,
    const std::vector<ARWEvent>& rec,
    const ARWWeightModel* model,
    size_t nrows) {

    std::vector<double> ng(nrows,0.0),nr(nrows,0.0);
    for(const auto& e:gen){
        double w=e.base_weight;
        if(model) w*=model->value(e);
        if(e.row>=0 && (size_t)e.row<nrows) ng[(size_t)e.row]+=w;
    }
    for(const auto& e:rec){
        double w=e.base_weight;
        if(model) w*=model->value(e);
        if(e.row>=0 && (size_t)e.row<nrows) nr[(size_t)e.row]+=w;
    }

    std::vector<double> a(nrows,std::numeric_limits<double>::quiet_NaN());
    for(size_t r=0;r<nrows;++r)
        if(ng[r]>0.0) a[r]=nr[r]/ng[r];
    return a;
}

struct ARWBHCell {
    bool valid=false;
    double target_mass=0.0;
};

static std::map<std::pair<int,int>,ARWBHCell> arw_load_bh_grid(
    const std::string& path,const std::string& energy_tag) {

    std::map<std::pair<int,int>,ARWBHCell> out;
    std::ifstream in(path);
    if(!in.is_open()) return out;

    std::string line;
    if(!std::getline(in,line)) return out;
    const auto head=split_csv_line(line);
    std::unordered_map<std::string,int> idx;
    for(int i=0;i<(int)head.size();++i) idx[head[(size_t)i]]=i;

    auto col=[&](const char* n)->int{
        auto it=idx.find(n);
        return it==idx.end()?-1:it->second;
    };
    const int ce=col("energy_tag"),cr=col("row"),cs=col("subcell");
    const int cb=col("bh_xs"),cv=col("volume");
    if(ce<0||cr<0||cs<0||cb<0||cv<0) return out;

    while(std::getline(in,line)){
        if(line.empty()) continue;
        const auto v=split_csv_line(line);
        if((int)v.size()<=std::max({ce,cr,cs,cb,cv})) continue;
        if(v[(size_t)ce]!=energy_tag) continue;
        const int r=std::atoi(v[(size_t)cr].c_str());
        const int s=std::atoi(v[(size_t)cs].c_str());
        const double bh=std::atof(v[(size_t)cb].c_str());
        const double vol=std::atof(v[(size_t)cv].c_str());
        if(!(std::isfinite(bh)&&bh>0.0&&std::isfinite(vol)&&vol>0.0))
            continue;
        out[{r,s}]={true,bh*vol};
    }
    return out;
}

static std::vector<double> arw_bh_event_weights(
    const std::vector<ARWEvent>& gen,
    const std::map<std::pair<int,int>,ARWBHCell>& grid) {

    std::map<std::pair<int,int>,double> counts;
    for(const auto& e:gen) counts[{e.row,e.subcell}]+=e.base_weight;

    std::vector<double> weights(gen.size(),1.0);
    std::map<std::pair<int,int>,double> raw;
    for(const auto& kv:grid){
        const double n=counts[kv.first];
        if(n>0.0 && kv.second.valid)
            raw[kv.first]=kv.second.target_mass/n;
    }

    // Normalize the average generated-event weight to one.
    double sw=0.0,sn=0.0;
    for(const auto& e:gen){
        auto it=raw.find({e.row,e.subcell});
        if(it==raw.end()) continue;
        sw+=e.base_weight*it->second;
        sn+=e.base_weight;
    }
    const double norm=(sw>0.0&&sn>0.0)?sn/sw:1.0;
    for(auto& kv:raw) kv.second*=norm;

    for(size_t i=0;i<gen.size();++i){
        auto it=raw.find({gen[i].row,gen[i].subcell});
        if(it!=raw.end()) weights[i]=it->second;
    }
    return weights;
}

static std::vector<double> arw_acceptance_bh(
    const std::vector<ARWEvent>& gen,
    const std::vector<ARWEvent>& rec,
    const std::map<std::pair<int,int>,ARWBHCell>& grid,
    size_t nrows) {

    std::map<std::pair<int,int>,double> gen_count;
    for(const auto& e:gen)
        gen_count[{e.row,e.subcell}]+=e.base_weight;

    std::map<std::pair<int,int>,double> cell_weight;
    for(const auto& kv:grid){
        const double n=gen_count[kv.first];
        if(n>0.0 && kv.second.valid)
            cell_weight[kv.first]=kv.second.target_mass/n;
    }

    double sw=0.0,sn=0.0;
    for(const auto& e:gen){
        auto it=cell_weight.find({e.row,e.subcell});
        if(it==cell_weight.end()) continue;
        sw+=e.base_weight*it->second;
        sn+=e.base_weight;
    }
    const double norm=(sw>0.0&&sn>0.0)?sn/sw:1.0;
    for(auto& kv:cell_weight) kv.second*=norm;

    std::vector<double> ng(nrows,0.0),nr(nrows,0.0);
    for(const auto& e:gen){
        double w=e.base_weight;
        auto it=cell_weight.find({e.row,e.subcell});
        if(it!=cell_weight.end()) w*=it->second;
        if(e.row>=0&&(size_t)e.row<nrows) ng[(size_t)e.row]+=w;
    }
    for(const auto& e:rec){
        double w=e.base_weight;
        auto it=cell_weight.find({e.row,e.subcell});
        if(it!=cell_weight.end()) w*=it->second;
        if(e.row>=0&&(size_t)e.row<nrows) nr[(size_t)e.row]+=w;
    }

    std::vector<double> a(nrows,std::numeric_limits<double>::quiet_NaN());
    for(size_t r=0;r<nrows;++r)
        if(ng[r]>0.0) a[r]=nr[r]/ng[r];
    return a;
}

static double arw_population_stddev(const std::vector<double>& values) {
    std::vector<double> v;
    for(double x:values) if(std::isfinite(x)) v.push_back(x);
    if(v.size()<2) return std::numeric_limits<double>::quiet_NaN();
    double m=0.0;
    for(double x:v) m+=x;
    m/=v.size();
    double s=0.0;
    for(double x:v) s+=(x-m)*(x-m);
    return std::sqrt(s/v.size());
}

static int arw_ensure_column(CSV& csv,const std::string& name) {
    auto it=csv.index.find(name);
    if(it!=csv.index.end()) return it->second;
    const int idx=(int)csv.header.size();
    csv.header.push_back(name);
    csv.index[name]=idx;
    for(auto& row:csv.rows) row.push_back("");
    return idx;
}

static void arw_write_closure_canvas(
    const std::string& path,
    const std::string& period,
    const std::vector<ARWEvent>& data,
    const std::vector<ARWEvent>& rec,
    const ARWFitResult& fit) {

    TCanvas c(("c_arw_"+period).c_str(),"",1200,900);
    c.Divide(2,2,0.002,0.002);
    const std::array<std::string,4> xt={{
        "x_{B}","Q^{2} (GeV^{2})","|t| (GeV^{2})","#phi (deg)"
    }};

    for(int iv=0;iv<4;++iv){
        c.cd(iv+1);
        gPad->SetLeftMargin(.14);
        gPad->SetRightMargin(.04);
        gPad->SetBottomMargin(.14);
        gPad->SetTopMargin(.16);
        gPad->SetTicks(1,1);

        const auto dh=arw_hist(data,fit.model,iv,false);
        const ARWWeightModel unit=[&](){
            ARWWeightModel u=fit.model;
            for(auto& f:u.factors) std::fill(f.begin(),f.end(),1.0);
            return u;
        }();
        const auto mb=arw_hist(rec,unit,iv,false);
        const auto ma=arw_hist(rec,fit.model,iv,true);

        const auto& ax=fit.model.axes[(size_t)iv];
        TH1D frame(("h_arw_frame_"+std::to_string(iv)).c_str(),"",
                   std::max(1,ax.nbins()),
                   ax.edges.front(),ax.edges.back());
        frame.SetStats(0);
        double ymax=.0;
        for(double y:dh) ymax=std::max(ymax,y);
        for(double y:mb) ymax=std::max(ymax,y);
        for(double y:ma) ymax=std::max(ymax,y);
        frame.SetMinimum(0.0);
        frame.SetMaximum(std::max(0.02,1.25*ymax));
        frame.GetXaxis()->SetTitle(xt[(size_t)iv].c_str());
        frame.GetYaxis()->SetTitle("Normalized event fraction");
        frame.GetXaxis()->SetTitleSize(.046);
        frame.GetYaxis()->SetTitleSize(.043);
        frame.GetXaxis()->SetLabelSize(.038);
        frame.GetYaxis()->SetLabelSize(.036);
        frame.GetYaxis()->SetTitleOffset(1.35);
        frame.DrawCopy();

        TGraph gd,gb,ga;
        for(int ib=0;ib<ax.nbins();++ib){
            const double x=ax.center(ib);
            gd.SetPoint(ib,x,dh[(size_t)ib]);
            gb.SetPoint(ib,x,mb[(size_t)ib]);
            ga.SetPoint(ib,x,ma[(size_t)ib]);
        }
        gd.SetMarkerStyle(20); gd.SetMarkerColor(kBlack);
        gd.SetLineColor(kBlack); gd.SetLineWidth(2);
        gb.SetMarkerStyle(24); gb.SetMarkerColor(kRed+1);
        gb.SetLineColor(kRed+1); gb.SetLineWidth(2);
        ga.SetMarkerStyle(25); ga.SetMarkerColor(kMagenta+2);
        ga.SetLineColor(kMagenta+2); ga.SetLineWidth(3);
        gd.DrawClone("LP SAME");
        gb.DrawClone("LP SAME");
        ga.DrawClone("LP SAME");

        TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.030);
        std::ostringstream ss;
        ss<<"D: "<<std::fixed<<std::setprecision(3)
          <<fit.before_distance[(size_t)iv]<<" #rightarrow "
          <<fit.after_distance[(size_t)iv];
        t.DrawLatex(.18,.78,ss.str().c_str());

        if(iv==0){
            TLegend l(.43,.62,.93,.80);
            l.SetBorderSize(0);l.SetFillStyle(0);
            l.SetTextFont(42);l.SetTextSize(.026);
            l.AddEntry(&gd,"DATA signal estimate","lp");
            l.AddEntry(&gb,"nominal reconstructed MC","lp");
            l.AddEntry(&ga,"iteratively reweighted MC","lp");
            l.DrawClone();
        }
    }

    c.cd(0);
    TLatex title;title.SetNDC();title.SetTextAlign(22);
    title.SetTextFont(42);title.SetTextSize(.022);
    title.DrawLatex(.50,.978,
        ("Iterative acceptance reweighting closure: "+period).c_str());
    title.SetTextSize(.0155);
    std::ostringstream ss;
    ss<<"Normalized shapes after nominal selection; "
      <<fit.iterations<<" iteration(s)";
    title.DrawLatex(.50,.952,ss.str().c_str());
    c.SaveAs(path.c_str());
}

// -----------------------------------------------------------------------------
// Synthetic-MC closure validation for reconstructed -> generated weight transfer
// -----------------------------------------------------------------------------

static double arw_test_weight(const ARWEvent& e,int test_id) {
    const double zx=std::max(-1.5,std::min(1.5,(e.x-0.30)/0.22));
    const double zq=std::max(-1.5,std::min(1.5,(e.q2-3.0)/2.0));
    const double zt=std::max(-1.5,std::min(1.5,(e.tabs-0.45)/0.38));
    const double ph=e.phi*3.14159265358979323846/180.0;

    double logw=0.0;
    switch(test_id){
        case 0: logw=+0.45*zx; break;
        case 1: logw=-0.40*zq; break;
        case 2: logw=+0.45*zt; break;
        case 3: logw=+0.35*std::cos(ph)+0.12*std::cos(2.0*ph); break;
        default:
            logw=+0.28*zx-0.22*zq+0.26*zt
                +0.24*std::cos(ph)+0.08*std::cos(2.0*ph);
            break;
    }
    return std::max(0.45,std::min(2.20,std::exp(logw)));
}

static const std::array<std::string,5>& arw_test_names() {
    static const std::array<std::string,5> names={{
        "xB tilt","Q2 tilt","|t| tilt","phi modulation","combined"
    }};
    return names;
}

static std::vector<ARWEvent> arw_make_pseudodata(
    const std::vector<ARWEvent>& rec,int test_id) {
    std::vector<ARWEvent> out=rec;
    for(auto& e:out) e.base_weight*=arw_test_weight(e,test_id);
    return out;
}

static std::vector<double> arw_acceptance_known_weight(
    const std::vector<ARWEvent>& gen,
    const std::vector<ARWEvent>& rec,
    int test_id,
    size_t nrows) {

    std::vector<double> ng(nrows,0.0),nr(nrows,0.0);
    for(const auto& e:gen){
        if(e.row<0||(size_t)e.row>=nrows) continue;
        ng[(size_t)e.row]+=e.base_weight*arw_test_weight(e,test_id);
    }
    for(const auto& e:rec){
        if(e.row<0||(size_t)e.row>=nrows) continue;
        nr[(size_t)e.row]+=e.base_weight*arw_test_weight(e,test_id);
    }
    std::vector<double> a(nrows,std::numeric_limits<double>::quiet_NaN());
    for(size_t r=0;r<nrows;++r)
        if(ng[r]>0.0) a[r]=nr[r]/ng[r];
    return a;
}

struct ARWSyntheticPoint {
    std::string period;
    std::string test;
    int row=-1;
    double nominal=std::numeric_limits<double>::quiet_NaN();
    double exact=std::numeric_limits<double>::quiet_NaN();
    double recovered=std::numeric_limits<double>::quiet_NaN();
    double exact_shift=std::numeric_limits<double>::quiet_NaN();
    double recovered_shift=std::numeric_limits<double>::quiet_NaN();
    double closure_bias=std::numeric_limits<double>::quiet_NaN();
};

static std::vector<ARWSyntheticPoint> arw_run_synthetic_closure(
    const std::string& period,
    const std::vector<ARWEvent>& gen,
    const std::vector<ARWEvent>& rec,
    const std::array<ARWFineAxis,4>& axes,
    const AcceptanceReweightingOptions& options,
    size_t nrows) {

    std::vector<ARWSyntheticPoint> out;
    const auto nominal=arw_acceptance(gen,rec,nullptr,nrows);

    for(int itest=0;itest<5;++itest){
        const auto pseudo=arw_make_pseudodata(rec,itest);
        const ARWFitResult fit=arw_fit(pseudo,rec,axes,options);
        const auto exact=arw_acceptance_known_weight(gen,rec,itest,nrows);
        const auto recovered=arw_acceptance(gen,rec,&fit.model,nrows);

        for(size_t r=0;r<nrows;++r){
            const double a0=nominal[r],ae=exact[r],ar=recovered[r];
            if(!(std::isfinite(a0)&&a0>0.0&&std::isfinite(ae)&&ae>0.0&&
                 std::isfinite(ar)&&ar>0.0)) continue;
            ARWSyntheticPoint p;
            p.period=period;
            p.test=arw_test_names()[(size_t)itest];
            p.row=(int)r;
            p.nominal=a0;
            p.exact=ae;
            p.recovered=ar;
            p.exact_shift=std::fabs(ae-a0)/a0;
            p.recovered_shift=std::fabs(ar-a0)/a0;
            p.closure_bias=std::fabs(ar-ae)/a0;
            out.push_back(std::move(p));
        }
    }
    return out;
}

static void arw_write_synthetic_csv(
    const std::string& path,
    const std::vector<ARWSyntheticPoint>& points) {
    std::ofstream out(path);
    out<<"period,test,row,acceptance_nominal,acceptance_exact,"
          "acceptance_recovered,exact_shift_frac,recovered_shift_frac,"
          "closure_bias_frac\\n";
    for(const auto& p:points)
        out<<p.period<<','<<p.test<<','<<p.row<<','
           <<p.nominal<<','<<p.exact<<','<<p.recovered<<','
           <<p.exact_shift<<','<<p.recovered_shift<<','
           <<p.closure_bias<<'\\n';
}

static double arw_quantile(std::vector<double> v,double q) {
    v.erase(std::remove_if(v.begin(),v.end(),
                           [](double x){return !std::isfinite(x);}),v.end());
    if(v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(),v.end());
    const double x=q*(v.size()-1);
    const size_t i=(size_t)std::floor(x),j=std::min(v.size()-1,i+1);
    return v[i]+(x-i)*(v[j]-v[i]);
}

static void arw_write_synthetic_canvas(
    const std::string& path,
    const std::vector<ARWSyntheticPoint>& points) {

    if(points.empty()) return;
    TCanvas c("c_arw_synthetic","",1250,900);
    c.Divide(2,2,.002,.002);

    std::vector<double> exact,recovered,bias,ratio;
    for(const auto& p:points){
        exact.push_back(100.0*p.exact_shift);
        recovered.push_back(100.0*p.recovered_shift);
        bias.push_back(100.0*p.closure_bias);
        if(p.exact_shift>0.002)
            ratio.push_back(p.closure_bias/p.exact_shift);
    }

    // (a) recovered vs exact acceptance shift.
    c.cd(1);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);
    gPad->SetTopMargin(.14);gPad->SetRightMargin(.04);gPad->SetTicks(1,1);
    double mx=0.0;for(double x:exact)mx=std::max(mx,x);
    for(double y:recovered)mx=std::max(mx,y);mx=std::max(1.0,1.10*mx);
    TH1D f1("h_arw_syn_frame1","",100,0,mx);f1.SetStats(0);
    f1.SetMinimum(0);f1.SetMaximum(mx);
    f1.GetXaxis()->SetTitle("Exact acceptance shift (%)");
    f1.GetYaxis()->SetTitle("Recovered acceptance shift (%)");
    f1.GetYaxis()->SetTitleOffset(1.35);f1.DrawCopy();
    TGraph g;for(int i=0;i<(int)exact.size();++i)g.SetPoint(i,exact[i],recovered[i]);
    g.SetMarkerStyle(20);g.SetMarkerSize(.45);g.SetMarkerColor(kBlue+1);g.DrawClone("P SAME");
    TLine one(0,0,mx,mx);one.SetLineStyle(2);one.SetLineWidth(2);one.DrawClone();
    TLatex lab;lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.032);lab.DrawLatex(.18,.84,"(a)");

    // (b) closure bias distribution.
    c.cd(2);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);
    gPad->SetTopMargin(.14);gPad->SetRightMargin(.04);gPad->SetTicks(1,1);
    const double b95=arw_quantile(bias,.95);const double bmax=std::max(0.5,1.35*b95);
    TH1D hb("h_arw_syn_bias","",60,0,bmax);hb.SetStats(0);
    for(double x:bias) if(x<=bmax) hb.Fill(x);
    if(hb.Integral()>0)hb.Scale(1.0/hb.Integral());
    hb.GetXaxis()->SetTitle("|A_{recovered}-A_{exact}|/A_{0} (%)");
    hb.GetYaxis()->SetTitle("Fraction of closure bins");hb.GetYaxis()->SetTitleOffset(1.35);
    hb.SetLineWidth(3);hb.SetLineColor(kMagenta+2);hb.DrawCopy("HIST");
    lab.DrawLatex(.18,.84,"(b)");

    // (c) bias relative to imposed exact shift.
    c.cd(3);gPad->SetLeftMargin(.14);gPad->SetBottomMargin(.14);
    gPad->SetTopMargin(.14);gPad->SetRightMargin(.04);gPad->SetTicks(1,1);
    const double r95=arw_quantile(ratio,.95);const double rmax=std::max(.25,1.30*r95);
    TH1D hr("h_arw_syn_ratio","",60,0,rmax);hr.SetStats(0);
    for(double x:ratio) if(x<=rmax)hr.Fill(x);
    if(hr.Integral()>0)hr.Scale(1.0/hr.Integral());
    hr.GetXaxis()->SetTitle("Closure bias / imposed acceptance shift");
    hr.GetYaxis()->SetTitle("Fraction of closure bins");hr.GetYaxis()->SetTitleOffset(1.35);
    hr.SetLineWidth(3);hr.SetLineColor(kGreen+2);hr.DrawCopy("HIST");
    lab.DrawLatex(.18,.84,"(c)");

    // (d) median closure bias by test shape.
    c.cd(4);gPad->SetLeftMargin(.17);gPad->SetBottomMargin(.20);
    gPad->SetTopMargin(.14);gPad->SetRightMargin(.04);gPad->SetTicks(1,1);
    TH1D ht("h_arw_syn_tests","",5,0.5,5.5);ht.SetStats(0);
    for(int it=0;it<5;++it){
        std::vector<double> vb;
        for(const auto& p:points)if(p.test==arw_test_names()[(size_t)it])vb.push_back(100.0*p.closure_bias);
        ht.SetBinContent(it+1,arw_quantile(vb,.50));
        ht.GetXaxis()->SetBinLabel(it+1,arw_test_names()[(size_t)it].c_str());
    }
    ht.GetXaxis()->SetTitle("Synthetic distortion");
    ht.GetYaxis()->SetTitle("Median closure bias (%)");ht.GetYaxis()->SetTitleOffset(1.50);
    ht.SetMinimum(0.0);ht.SetMaximum(std::max(0.5,1.35*ht.GetMaximum()));
    ht.SetMarkerStyle(20);ht.SetMarkerSize(1.2);ht.SetLineWidth(2);ht.DrawCopy("P");
    ht.GetXaxis()->LabelsOption("v");lab.DrawLatex(.20,.84,"(d)");

    c.cd(0);TLatex title;title.SetNDC();title.SetTextAlign(22);title.SetTextFont(42);
    title.SetTextSize(.022);title.DrawLatex(.50,.978,
        "Synthetic closure of reconstructed-to-generated acceptance reweighting");
    title.SetTextSize(.0155);title.DrawLatex(.50,.952,
        "Known smooth MC distortions; analytic-weight acceptance compared with iterative recovery");
    c.SaveAs(path.c_str());
}

static void arw_write_candidate_distribution(
    const std::string& path,
    const std::vector<double>& f106,
    const std::vector<double>& f102) {
    TCanvas c("c_arw_candidate_dist","",1050,720);
    c.SetLeftMargin(.13);c.SetRightMargin(.04);c.SetBottomMargin(.13);c.SetTopMargin(.18);c.SetTicks(1,1);
    std::vector<double> all=f106;all.insert(all.end(),f102.begin(),f102.end());
    const double p95=arw_quantile(all,.95);const double xmax=std::max(2.0,1.35*100.0*p95);
    TH1D h1("h_arw_cand106","",60,0,xmax),h2("h_arw_cand102","",60,0,xmax);
    h1.SetStats(0);h2.SetStats(0);
    for(double x:f106)if(std::isfinite(x)&&100*x<=xmax)h1.Fill(100*x);
    for(double x:f102)if(std::isfinite(x)&&100*x<=xmax)h2.Fill(100*x);
    if(h1.Integral()>0)h1.Scale(1.0/h1.Integral());if(h2.Integral()>0)h2.Scale(1.0/h2.Integral());
    double ymax=1.2*std::max(h1.GetMaximum(),h2.GetMaximum());
    TH1D fr("h_arw_cand_frame","",60,0,xmax);fr.SetStats(0);fr.SetMinimum(0);fr.SetMaximum(ymax);
    fr.GetXaxis()->SetTitle("DATA-reweighting acceptance systematic (%)");
    fr.GetYaxis()->SetTitle("Fraction of populated bins");fr.GetYaxis()->SetTitleOffset(1.25);fr.DrawCopy();
    h1.SetLineColor(kBlue+1);h1.SetLineWidth(3);h2.SetLineColor(kRed+1);h2.SetLineWidth(3);h2.SetLineStyle(2);
    h1.DrawClone("HIST SAME");h2.DrawClone("HIST SAME");
    TLegend l(.58,.66,.92,.80);l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.029);
    l.AddEntry(&h1,"10.6 GeV combined","l");l.AddEntry(&h2,"10.2 GeV Sp19 Inb","l");l.DrawClone();
    TLatex t;t.SetNDC();t.SetTextFont(42);t.SetTextSize(.030);
    std::ostringstream ss;ss<<"median: "<<std::fixed<<std::setprecision(2)<<100*arw_quantile(f106,.50)
                            <<"% (10.6), "<<100*arw_quantile(f102,.50)<<"% (10.2)";
    t.DrawLatex(.16,.86,ss.str().c_str());
    t.SetTextAlign(22);t.SetTextSize(.022);t.DrawLatex(.50,.965,"Pass-2 DATA-driven acceptance model uncertainty");
    c.SaveAs(path.c_str());
}

static void arw_write_candidate_kinematic_summary(
    const std::string& path,
    const std::vector<RowBin>& rows,
    const CSV& csv,
    int c106,
    int c102) {

    struct VSpec { const char* title; int which; };
    const std::array<VSpec,4> vs={{
        {"x_{B}",0},{"Q^{2} (GeV^{2})",1},
        {"|t| (GeV^{2})",2},{"#phi (deg)",3}
    }};

    TCanvas c("c_arw_candidate_kin","",1250,900);
    c.Divide(2,2,.002,.002);

    for(int iv=0;iv<4;++iv){
        std::map<std::pair<double,double>,std::vector<double>> a,b;
        for(size_t r=0;r<rows.size()&&r<csv.rows.size();++r){
            if(!rows[r].valid) continue;
            std::pair<double,double> key;
            if(iv==0) key={rows[r].xBmin,rows[r].xBmax};
            else if(iv==1) key={rows[r].Q2min,rows[r].Q2max};
            else if(iv==2) key={rows[r].tmin,rows[r].tmax};
            else key={rows[r].pmin,rows[r].pmax};
            const double f1=arw_tuple_first(csv.rows[r][(size_t)c106]);
            const double f2=arw_tuple_first(csv.rows[r][(size_t)c102]);
            if(std::isfinite(f1)) a[key].push_back(100.0*f1);
            if(std::isfinite(f2)) b[key].push_back(100.0*f2);
        }

        std::vector<double> xa,ya,xb,yb;
        double ymax=0.0;
        for(const auto& kv:a){
            const double y=arw_quantile(kv.second,.50);
            xa.push_back(.5*(kv.first.first+kv.first.second));ya.push_back(y);ymax=std::max(ymax,y);
        }
        for(const auto& kv:b){
            const double y=arw_quantile(kv.second,.50);
            xb.push_back(.5*(kv.first.first+kv.first.second));yb.push_back(y);ymax=std::max(ymax,y);
        }
        if(xa.empty()&&xb.empty()) continue;
        double xmin=1e9,xmax=-1e9;
        for(double x:xa){xmin=std::min(xmin,x);xmax=std::max(xmax,x);}for(double x:xb){xmin=std::min(xmin,x);xmax=std::max(xmax,x);}
        const double dx=(xmax>xmin)?.04*(xmax-xmin):.5;xmin-=dx;xmax+=dx;

        c.cd(iv+1);gPad->SetLeftMargin(.14);gPad->SetRightMargin(.04);
        gPad->SetBottomMargin(.14);gPad->SetTopMargin(.16);gPad->SetTicks(1,1);
        TH1D fr(("h_arw_cand_kin_"+std::to_string(iv)).c_str(),"",100,xmin,xmax);
        fr.SetStats(0);fr.SetMinimum(0);fr.SetMaximum(std::max(1.0,1.25*ymax));
        fr.GetXaxis()->SetTitle(vs[(size_t)iv].title);
        fr.GetYaxis()->SetTitle("Median acceptance systematic (%)");fr.GetYaxis()->SetTitleOffset(1.35);fr.DrawCopy();
        TGraph g1((int)xa.size(),xa.data(),ya.data()),g2((int)xb.size(),xb.data(),yb.data());
        g1.SetMarkerStyle(20);g1.SetMarkerColor(kBlue+1);g1.SetLineColor(kBlue+1);g1.SetLineWidth(3);
        g2.SetMarkerStyle(24);g2.SetMarkerColor(kRed+1);g2.SetLineColor(kRed+1);g2.SetLineWidth(3);
        if(!xa.empty())g1.DrawClone("LP SAME");if(!xb.empty())g2.DrawClone("LP SAME");
        TLatex lab;lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.032);std::string pn="("+std::string(1,char('a'+iv))+")";lab.DrawLatex(.18,.84,pn.c_str());
        if(iv==0){TLegend l(.53,.65,.92,.80);l.SetBorderSize(0);l.SetFillStyle(0);l.SetTextFont(42);l.SetTextSize(.027);l.AddEntry(&g1,"10.6 GeV combined","lp");l.AddEntry(&g2,"10.2 GeV Sp19 Inb","lp");l.DrawClone();}
    }
    c.cd(0);TLatex t;t.SetNDC();t.SetTextAlign(22);t.SetTextFont(42);t.SetTextSize(.022);
    t.DrawLatex(.50,.978,"Kinematic dependence of DATA-reweighted acceptance systematic");
    t.SetTextSize(.0155);t.DrawLatex(.50,.952,"Full nominal-to-reweighted acceptance excursion; medians in each kinematic interval");
    c.SaveAs(path.c_str());
}


} // namespace

bool run_acceptance_reweighting_study(
    const std::string& csv_path,
    const std::map<std::string,TTree*>& dvcsDataTrees,
    const std::map<std::string,TTree*>& dvcsGenMcTrees,
    const std::map<std::string,TTree*>& dvcsRecMcTrees,
    const AcceptanceReweightingOptions& options) {

    try{
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gStyle->SetOptStat(0);
        std::filesystem::create_directories(options.output_dir);

        CSV csv;
        load_csv(csv_path,csv);
        const auto rows=load_row_bins_from_csv(csv);
        const auto fast_bins=build_fast_binning(rows);
        const auto axes=arw_build_axes(rows);
        const TopoCutMap cuts_data=
            load_combined_cuts(options.combined_cuts_json,"data");
        const TopoCutMap cuts_mc=
            load_combined_cuts(options.combined_cuts_json,"mc");

        CurrentResponseModel current_model;
        const CurrentResponseModel* current_ptr=nullptr;
        if(options.apply_current_correction){
            current_model=load_current_response_model(
                options.current_response_model_json);
            current_ptr=&current_model;
        }

        if(options.enable_bh_reweighting &&
           options.build_bh_grid_if_missing &&
           !std::filesystem::exists(options.bh_grid_csv)){
            std::ostringstream cmd;
            cmd<<"python3 "<<options.bh_grid_script
               <<" --analysis-csv "<<csv_path
               <<" --output "<<options.bh_grid_csv
               <<" --workers "<<options.bh_grid_workers
               <<" --subdivisions 2";
            std::cout<<"[acceptance-reweighting] Building pure-BH subcell grid.\n";
            const int rc=std::system(cmd.str().c_str());
            if(rc!=0)
                std::cerr<<"[acceptance-reweighting] WARNING: BH-grid helper "
                         <<"returned "<<rc<<"; continuing data-driven only.\n";
        }

        WorkConfig data_cfg; data_cfg.channel_cfg=dvcs_config();
        data_cfg.sample_kind=SampleKind::DATA;
        WorkConfig gen_cfg; gen_cfg.channel_cfg=dvcs_config();
        gen_cfg.sample_kind=SampleKind::MC_GEN;
        WorkConfig rec_cfg; rec_cfg.channel_cfg=dvcs_config();
        rec_cfg.sample_kind=SampleKind::MC_REC;

        const std::array<std::string,5> periods={{
            "Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out","Sp19 Inb"
        }};

        const int c_candidate_106=arw_ensure_column(
            csv,"acceptance reweighting candidate sys frac, 10.6 GeV");
        const int c_candidate_102=arw_ensure_column(
            csv,"acceptance reweighting candidate sys frac, Sp19 Inb");
        const int c_data_106=arw_ensure_column(
            csv,"acceptance reweighting data-driven sys frac, 10.6 GeV");
        const int c_data_102=arw_ensure_column(
            csv,"acceptance reweighting data-driven sys frac, Sp19 Inb");

        struct PeriodResult {
            std::vector<double> nominal,data_rw,candidate;
        };
        std::map<std::string,PeriodResult> results;

        std::ofstream summary(
            std::filesystem::path(options.output_dir)/
            "acceptance_reweighting_per_period_summary.csv");
        summary<<"period,row,acceptance_nominal,acceptance_data_reweighted,"
               <<"relative_acceptance_shift,candidate_cross_section_fraction\n";

        std::vector<ARWSyntheticPoint> synthetic_points;

        for(const auto& period:periods){
            TTree* dt=arw_tree_for_period(dvcsDataTrees,period);
            TTree* gt=arw_tree_for_period(dvcsGenMcTrees,period);
            TTree* rt=arw_tree_for_period(dvcsRecMcTrees,period);
            if(!dt||!gt||!rt){
                std::cerr<<"[acceptance-reweighting] WARNING: missing tree(s) "
                         <<"for "<<period<<".\n";
                continue;
            }

            const PeriodTags dtag=parse_period_tags_from_tree_key(
                [&](){
                    for(const auto& kv:dvcsDataTrees)
                        if(kv.second==dt) return kv.first;
                    return std::string();
                }());
            const PeriodTags gtag=parse_period_tags_from_tree_key(
                [&](){
                    for(const auto& kv:dvcsGenMcTrees)
                        if(kv.second==gt) return kv.first;
                    return std::string();
                }());
            const PeriodTags rtag=parse_period_tags_from_tree_key(
                [&](){
                    for(const auto& kv:dvcsRecMcTrees)
                        if(kv.second==rt) return kv.first;
                    return std::string();
                }());

            std::cout<<"[acceptance-reweighting] Collecting nominal selected "
                     <<"events for "<<period<<"...\n";
            const auto data=arw_collect(
                dt,data_cfg,dtag,rows,fast_bins,cuts_data,current_ptr,
                csv,options,false);
            const auto rec=arw_collect(
                rt,rec_cfg,rtag,rows,fast_bins,cuts_mc,current_ptr,
                csv,options,false);
            const auto gen=arw_collect(
                gt,gen_cfg,gtag,rows,fast_bins,cuts_mc,nullptr,
                csv,options,true);

            std::cout<<"[acceptance-reweighting] "<<period
                     <<" selected: data="<<data.size()
                     <<" recMC="<<rec.size()
                     <<" genMC="<<gen.size()<<"\n";

            const ARWFitResult fit=arw_fit(data,rec,axes,options);
            arw_write_closure_canvas(
                (std::filesystem::path(options.output_dir)/
                 ("shape_closure_"+dtag.period_code+".png")).string(),
                period,data,rec,fit);

            PeriodResult pr;
            pr.nominal=arw_acceptance(
                gen,rec,nullptr,csv.rows.size());
            pr.data_rw=arw_acceptance(
                gen,rec,&fit.model,csv.rows.size());

            pr.candidate.assign(csv.rows.size(),
                std::numeric_limits<double>::quiet_NaN());

            for(size_t r=0;r<csv.rows.size();++r){
                const double a0=pr.nominal[r];
                const double ad=pr.data_rw[r];
                if(!(std::isfinite(a0)&&a0>0.0&&std::isfinite(ad)&&ad>0.0))
                    continue;

                // Final candidate definition under study: full DATA-reweighting
                // excursion relative to the nominal acceptance.  Pure BH is not
                // included because it is a stress-test shape, not our estimate
                // of the physical DVCS+BH event distribution.
                const double rd_acceptance=std::fabs(ad-a0)/a0;
                const double rd_cross_section=std::fabs(a0/ad-1.0);
                pr.candidate[r]=rd_cross_section;

                summary<<period<<','<<r<<','<<a0<<','<<ad<<','
                       <<rd_acceptance<<','<<rd_cross_section<<'\n';
            }

            std::cout<<"[acceptance-reweighting] Running synthetic transfer "
                     <<"closure for "<<period<<"...\n";
            const auto syn=arw_run_synthetic_closure(
                period,gen,rec,axes,options,csv.rows.size());
            synthetic_points.insert(synthetic_points.end(),syn.begin(),syn.end());

            results[period]=std::move(pr);
        }

        // Effective 10.6-GeV candidate: weight period-level acceptance shifts by
        // the nominal acceptance-corrected unpolarized yield contribution.
        for(size_t r=0;r<csv.rows.size();++r){
            double ytot=0.0;
            double delta_data=0.0;
            double nominal_combined=0.0;
            bool have_any=false;

            for(const std::string period:
                {"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"}){
                auto ir=results.find(period);
                if(ir==results.end()) continue;

                const std::string ycol=
                    "acceptance corrected yield, ep->epg, exp, "+
                    period+", unpol";
                auto iy=csv.index.find(ycol);
                if(iy==csv.index.end()) continue;
                const double y=arw_tuple_first(
                    csv.rows[r][(size_t)iy->second]);
                if(!(std::isfinite(y)&&y>=0.0)) continue;

                const double a0=ir->second.nominal[r];
                const double ad=ir->second.data_rw[r];
                if(!(std::isfinite(a0)&&a0>0.0&&std::isfinite(ad)&&ad>0.0))
                    continue;

                ytot+=y;
                nominal_combined+=y;
                delta_data+=y*(a0/ad);
                have_any=true;
            }

            if(have_any && ytot>0.0){
                const double rd=std::fabs(delta_data-nominal_combined)/ytot;
                csv.rows[r][(size_t)c_data_106]=std::to_string(rd);
                csv.rows[r][(size_t)c_candidate_106]=std::to_string(rd);
            }

            auto sp=results.find("Sp19 Inb");
            if(sp!=results.end()){
                const double a0=sp->second.nominal[r];
                const double ad=sp->second.data_rw[r];
                if(std::isfinite(a0)&&a0>0.0&&std::isfinite(ad)&&ad>0.0){
                    const double rd=std::fabs(a0/ad-1.0);
                    csv.rows[r][(size_t)c_data_102]=std::to_string(rd);
                    csv.rows[r][(size_t)c_candidate_102]=std::to_string(rd);
                }
            }
        }

        // Synthetic closure and candidate-distribution diagnostics.
        arw_write_synthetic_csv(
            (std::filesystem::path(options.output_dir)/
             "synthetic_reweighting_closure.csv").string(),
            synthetic_points);
        arw_write_synthetic_canvas(
            (std::filesystem::path(options.output_dir)/
             "synthetic_reweighting_closure_summary.png").string(),
            synthetic_points);

        std::vector<double> candidate106,candidate102;
        for(const auto& row:csv.rows){
            const double a=arw_tuple_first(row[(size_t)c_candidate_106]);
            const double b=arw_tuple_first(row[(size_t)c_candidate_102]);
            if(std::isfinite(a)) candidate106.push_back(a);
            if(std::isfinite(b)) candidate102.push_back(b);
        }
        arw_write_candidate_distribution(
            (std::filesystem::path(options.output_dir)/
             "acceptance_data_reweighting_systematic_distribution.png").string(),
            candidate106,candidate102);
        arw_write_candidate_kinematic_summary(
            (std::filesystem::path(options.output_dir)/
             "acceptance_data_reweighting_kinematic_summary.png").string(),
            rows,csv,c_candidate_106,c_candidate_102);

        if(!synthetic_points.empty()){
            std::vector<double> vbias,vratio;
            for(const auto& p:synthetic_points){
                vbias.push_back(p.closure_bias);
                if(p.exact_shift>0.002) vratio.push_back(p.closure_bias/p.exact_shift);
            }
            std::cout<<"[acceptance-reweighting] Synthetic closure: median bias="
                     <<100.0*arw_quantile(vbias,.50)<<"%, 95th percentile="
                     <<100.0*arw_quantile(vbias,.95)<<"%, median bias/imposed-shift="
                     <<arw_quantile(vratio,.50)<<".\n";
        }

        if(options.install_candidate_as_production_systematic){
            auto ia=csv.index.find("Syst. err (Acceptance)");
            auto ix=csv.index.find(
                "normed cross sections, ep->epg, exp, 10.6 GeV, unpol");
            if(ia==csv.index.end()||ix==csv.index.end())
                throw std::runtime_error(
                    "Cannot install acceptance candidate: required columns missing.");
            for(size_t r=0;r<csv.rows.size();++r){
                const double frac=arw_tuple_first(
                    csv.rows[r][(size_t)c_candidate_106]);
                const double xs=arw_tuple_first(
                    csv.rows[r][(size_t)ix->second]);
                if(std::isfinite(frac)&&std::isfinite(xs))
                    csv.rows[r][(size_t)ia->second]=std::to_string(
                        std::fabs(xs)*frac);
            }
        }

        write_csv_atomic(csv_path,csv);

        std::cout<<"[acceptance-reweighting] Study complete. Candidate columns "
                 <<"were written to "<<csv_path<<".\n";
        std::cout<<"[acceptance-reweighting] Production acceptance systematic "
                 <<(options.install_candidate_as_production_systematic
                    ?"WAS":"was NOT")
                 <<" replaced.\n";
        return true;
    }catch(const std::exception& e){
        std::cerr<<"[acceptance-reweighting] ERROR: "<<e.what()<<"\n";
        return false;
    }
}
