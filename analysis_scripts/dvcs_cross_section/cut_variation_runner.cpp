#include "cut_variation_runner.h"

#include "acceptance.h"
#include "cross_sections.h"
#include "cut_variation_systematics.h"
#include "exclusivity_cuts.h"
#include "global_cuts.h"
#include "norm_cross_sections.h"
#include "pi0_contamination.h"
#include "pi0_corrected_counts.h"
#include "total_counts.h"
#include "unfolding.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string, std::size_t> index;
};

struct TripleCell {
    bool ok = false;
    double value = 0.0;
    double stat = 0.0;
    double sys = 0.0;
};

std::vector<std::string> split_csv(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                quoted = !quoted;
            }
        } else if (c == ',' && !quoted) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

std::string csv_escape(const std::string& s) {
    if (s.find_first_of(",\"\n\r") == std::string::npos) return s;
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out.push_back(c);
    }
    out += "\"";
    return out;
}

CsvTable read_csv(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open CSV: " + path);
    CsvTable t;
    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("empty CSV: " + path);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    t.header = split_csv(line);
    for (std::size_t i = 0; i < t.header.size(); ++i) t.index[t.header[i]] = i;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        auto row = split_csv(line);
        row.resize(t.header.size());
        t.rows.push_back(std::move(row));
    }
    return t;
}

void write_csv(const std::string& path, const CsvTable& t) {
    const std::string tmp = path + ".tmp";
    std::ofstream out(tmp);
    if (!out) throw std::runtime_error("cannot write CSV: " + tmp);
    for (std::size_t i = 0; i < t.header.size(); ++i) {
        if (i) out << ',';
        out << csv_escape(t.header[i]);
    }
    out << '\n';
    for (const auto& row : t.rows) {
        for (std::size_t i = 0; i < t.header.size(); ++i) {
            if (i) out << ',';
            out << csv_escape(i < row.size() ? row[i] : std::string());
        }
        out << '\n';
    }
    out.close();
    fs::rename(tmp, path);
}

TripleCell parse_triple(const std::string& cell) {
    TripleCell v;
    std::string s = cell;
    s.erase(std::remove_if(s.begin(), s.end(), [](unsigned char c){ return std::isspace(c); }), s.end());
    if (s.empty()) return v;
    if (s.front() == '(' && s.back() == ')') s = s.substr(1, s.size() - 2);
    std::replace(s.begin(), s.end(), ';', ',');
    std::stringstream ss(s);
    std::string token;
    std::vector<double> vals;
    while (std::getline(ss, token, ',')) {
        try { vals.push_back(std::stod(token)); } catch (...) { return v; }
    }
    if (vals.empty()) return v;
    v.ok = std::isfinite(vals[0]);
    v.value = vals[0];
    v.stat = vals.size() > 1 ? vals[1] : 0.0;
    v.sys = vals.size() > 2 ? vals[2] : 0.0;
    return v;
}

std::string format_triple(double value, double stat, double sys) {
    std::ostringstream os;
    os << '(' << std::setprecision(12) << value << ',' << stat << ',' << sys << ')';
    return os.str();
}

bool starts_with(const std::string& s, const std::string& prefix) {
    return s.rfind(prefix, 0) == 0;
}

// The nominal CSV contains the exact effective correction applied in every bin.
// Preserve that correction under a cut variation by inferring target/source from
// the nominal CSV, then applying the same multiplicative factor to the varied source.
void reapply_nominal_fixed_corrections(const std::string& nominal_path,
                                       const std::string& varied_path) {
    const CsvTable nominal = read_csv(nominal_path);
    CsvTable varied = read_csv(varied_path);
    if (nominal.header != varied.header || nominal.rows.size() != varied.rows.size()) {
        throw std::runtime_error("nominal/variation CSV layouts differ");
    }

    for (std::size_t target_col = 0; target_col < varied.header.size(); ++target_col) {
        const std::string& target = varied.header[target_col];
        std::string source;
        if (starts_with(target, "normalized raw yield, ")) {
            source = "raw yield, " + target.substr(std::string("normalized raw yield, ").size());
        } else if (starts_with(target, "reconstructed current corrected yield, ")) {
            source = "reconstructed yield, " + target.substr(std::string("reconstructed current corrected yield, ").size());
        } else {
            continue;
        }

        auto it = varied.index.find(source);
        if (it == varied.index.end()) continue;
        const std::size_t source_col = it->second;

        for (std::size_t r = 0; r < varied.rows.size(); ++r) {
            const TripleCell nsrc = parse_triple(nominal.rows[r][source_col]);
            const TripleCell ntgt = parse_triple(nominal.rows[r][target_col]);
            const TripleCell vsrc = parse_triple(varied.rows[r][source_col]);
            if (!nsrc.ok || !ntgt.ok || !vsrc.ok || std::abs(nsrc.value) < 1e-30) {
                varied.rows[r][target_col].clear();
                continue;
            }
            const double scale = ntgt.value / nsrc.value;
            varied.rows[r][target_col] = format_triple(
                vsrc.value * scale,
                std::abs(scale) * vsrc.stat,
                std::abs(scale) * vsrc.sys);
        }
    }
    write_csv(varied_path, varied);
}

struct VariationSpec {
    std::string name;
    std::string csv_name;
    double sigma = 3.0;
    double quantile = 0.99;
    int fiducial_direction = 0; // -1 loose, 0 nominal, +1 tight
};

GlobalCutConfig varied_global_config(const GlobalCutConfig& nominal, int direction) {
    GlobalCutConfig cfg = nominal;
    if (direction == 0) return cfg;
    const double d = 2.0 * static_cast<double>(direction);
    cfg.auxiliary_e_theta_min_deg += d;
    cfg.auxiliary_e_theta_max_deg -= d;
    cfg.auxiliary_fd_proton_theta_min_deg += d;
    cfg.auxiliary_fd_proton_theta_max_deg -= d;
    cfg.auxiliary_fd_photon_theta_min_deg += d;
    cfg.auxiliary_fd_photon_theta_max_deg -= d;
    cfg.auxiliary_cd_proton_theta_min_deg += d;
    cfg.auxiliary_cd_proton_theta_max_deg -= d;
    return cfg;
}

bool produce_variation(
    const VariationSpec& spec,
    const AutomaticCutVariationOptions& options,
    const GlobalCutConfig& nominal_global_cuts,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::map<std::string, TTree*>& currentStudyGenMcTrees,
    const std::map<std::string, TTree*>& currentStudyRecMcTrees,
    bool use_nobkg_dvcs_mc_for_acceptance,
    bool use_epg_mc_current_factor_for_eppi0_bkg) {

    const fs::path base = fs::path(options.output_dir) / "variations" / spec.name;
    const fs::path json_dir = base / "jsons";
    const fs::path plots_dir = base / "exclusivity_plots";
    const fs::path stage_out = base / "analysis_output";
    const fs::path csv_dir = fs::path(options.output_dir) / "csv";
    const fs::path csv_path = csv_dir / spec.csv_name;
    fs::create_directories(json_dir);
    fs::create_directories(plots_dir);
    fs::create_directories(stage_out);
    fs::create_directories(csv_dir);
    fs::copy_file(options.nominal_csv, csv_path, fs::copy_options::overwrite_existing);

    const GlobalCutConfig cfg = varied_global_config(nominal_global_cuts, spec.fiducial_direction);
    set_default_global_cuts(cfg);
    write_global_cuts_config_json(json_dir.string());

    ExclusivityDiagnosticConfig excl_cfg;
    excl_cfg.symmetric_sigma_multiplier = spec.sigma;
    excl_cfg.upper_tail_quantile = spec.quantile;
    excl_cfg.enable = false;
    runAllExclusivityCuts(dataTrees, recMcTrees, eppi0DataTrees, eppi0RecMcTrees,
                          json_dir.string(),
                          plots_dir.string(),
                          options.max_workers, excl_cfg);

    const std::string cuts_json = (json_dir / "combined_cuts.json").string();
    TotalCountsOptions count_opts;
    count_opts.use_nobkg_dvcs_mc_counts = use_nobkg_dvcs_mc_for_acceptance;
    if (!update_total_counts_csv(csv_path.string(), dataTrees, eppi0DataTrees,
                                 genMcTrees, recMcTrees, eppi0GenMcTrees,
                                 eppi0RecMcTrees, eppi0BkgTrees, cuts_json,
                                 stage_out.string(), options.max_workers, count_opts,
                                 currentStudyGenMcTrees, currentStudyRecMcTrees)) return false;

    reapply_nominal_fixed_corrections(options.nominal_csv, csv_path.string());

    Pi0ContaminationOptions pi0_opts;
    pi0_opts.use_epg_mc_current_factor_for_eppi0_bkg = use_epg_mc_current_factor_for_eppi0_bkg;
    if (!compute_pi0_contamination_overall(dataTrees, eppi0DataTrees, eppi0RecMcTrees,
                                           eppi0BkgTrees, cuts_json, csv_path.string(),
                                           stage_out.string(), 1, pi0_opts)) return false;
    if (!update_pi0_corrected_counts_csv(csv_path.string(), stage_out.string())) return false;
    if (!update_acceptance_csv(csv_path.string(), genMcTrees, recMcTrees, cuts_json,
                               (json_dir / "global_cuts_config.json").string(),
                               stage_out.string())) return false;
    if (!update_unfolded_yields_csv(csv_path.string(), stage_out.string())) return false;

    LumiBuildOptions lumi_opts;
    lumi_opts.use_second_column_charge_for_all_unpolarized = true;
    lumi_opts.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized = false;
    lumi_opts.columns_3_to_5_charge_sum_scale = 1.025;
    const LumiMap lumi = build_lumi_map(lumi_opts);
    if (!compute_cross_sections(csv_path.string(), lumi)) return false;
    if (!update_normed_cross_sections_csv(csv_path.string())) return false;

    std::cout << "[cut-variation-runner] Completed " << spec.name
              << " -> " << csv_path.string() << std::endl;
    return true;
}

} // namespace

bool run_automatic_cut_variation_systematics(
    const AutomaticCutVariationOptions& options,
    const GlobalCutConfig& nominal_global_cuts,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0GenMcTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::map<std::string, TTree*>& currentStudyGenMcTrees,
    const std::map<std::string, TTree*>& currentStudyRecMcTrees,
    bool use_nobkg_dvcs_mc_for_acceptance,
    bool use_epg_mc_current_factor_for_eppi0_bkg) {

    if (!options.enabled) return true;
    try {
        fs::create_directories(options.output_dir);
        fs::copy_file(options.nominal_csv,
                      fs::path(options.output_dir) / "nominal_csv_before_cut_systematics.csv",
                      fs::copy_options::overwrite_existing);

        const std::vector<VariationSpec> variations = {
            {"exclusivity_95", "exclusivity_95.csv", 2.0, 0.95, 0},
            {"exclusivity_99p99", "exclusivity_99p99.csv", 4.0, 0.9999, 0},
            {"fiducial_loose", "fiducial_loose.csv", 3.0, 0.99, -1},
            {"fiducial_tight", "fiducial_tight.csv", 3.0, 0.99, +1}
        };

        for (const auto& v : variations) {
            if (!produce_variation(v, options, nominal_global_cuts,
                                   dataTrees, genMcTrees, recMcTrees,
                                   eppi0DataTrees, eppi0GenMcTrees,
                                   eppi0RecMcTrees, eppi0BkgTrees,
                                   currentStudyGenMcTrees, currentStudyRecMcTrees,
                                   use_nobkg_dvcs_mc_for_acceptance,
                                   use_epg_mc_current_factor_for_eppi0_bkg)) {
                set_default_global_cuts(nominal_global_cuts);
                return false;
            }
        }
        set_default_global_cuts(nominal_global_cuts);

        CutVariationSystematicsOptions syst;
        syst.enabled = true;
        syst.apply_barlow = true;
        syst.barlow_threshold = 1.0;
        syst.make_plots = options.make_final_diagnostic_plots;
        syst.write_diagnostic_csv = true;
        syst.nominal_csv = options.nominal_csv;
        syst.exclusivity_loose_csv = (fs::path(options.output_dir) / "csv/exclusivity_95.csv").string();
        syst.exclusivity_tight_csv = (fs::path(options.output_dir) / "csv/exclusivity_99p99.csv").string();
        syst.fiducial_loose_csv = (fs::path(options.output_dir) / "csv/fiducial_loose.csv").string();
        syst.fiducial_tight_csv = (fs::path(options.output_dir) / "csv/fiducial_tight.csv").string();
        syst.output_dir = options.output_dir;
        return update_cut_variation_systematics(syst);
    } catch (const std::exception& e) {
        set_default_global_cuts(nominal_global_cuts);
        std::cerr << "[cut-variation-runner] ERROR: " << e.what() << std::endl;
        return false;
    }
}
