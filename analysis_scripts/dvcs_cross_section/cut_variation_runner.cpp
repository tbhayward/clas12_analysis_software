#include "cut_variation_runner.h"

#include "acceptance.h"
#include "bsa.h"
#include "cross_sections.h"
#include "cut_variation_systematics.h"
#include "python_exclusivity_runner.h"
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
#include <chrono>
#include <iomanip>
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

std::string normalize_csv_field(std::string field) {
    // A correctly parsed CSV field must not retain the quote characters used
    // solely to protect commas. Older variation files could contain one or
    // more extra wrapper-quote layers, e.g. """(1,2,0)""". Remove only
    // complete matching wrapper layers and collapse doubled quotes.
    bool changed = true;
    while (changed && field.size() >= 2) {
        changed = false;
        if (field.front() == '"' && field.back() == '"') {
            field = field.substr(1, field.size() - 2);
            std::string unescaped;
            unescaped.reserve(field.size());
            for (std::size_t i = 0; i < field.size(); ++i) {
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
            out.push_back(normalize_csv_field(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(normalize_csv_field(cur));
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
    s = normalize_csv_field(s);
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

std::size_t ensure_column(CsvTable& t, const std::string& name) {
    auto it = t.index.find(name);
    if (it != t.index.end()) return it->second;
    const std::size_t idx = t.header.size();
    t.header.push_back(name);
    t.index[name] = idx;
    for (auto& row : t.rows) row.push_back("");
    return idx;
}

void update_bsa_cut_systematics(const AutomaticCutVariationOptions& options) {
    CsvTable nominal = read_csv(options.nominal_csv);
    const CsvTable excl_loose = read_csv((fs::path(options.output_dir) / "csv/exclusivity_loose_98.csv").string());
    const CsvTable excl_tight = read_csv((fs::path(options.output_dir) / "csv/exclusivity_tight_90.csv").string());
    const CsvTable fid_loose = read_csv((fs::path(options.output_dir) / "csv/fiducial_loose.csv").string());
    const CsvTable fid_tight = read_csv((fs::path(options.output_dir) / "csv/fiducial_tight.csv").string());
    if (nominal.rows.size() != excl_loose.rows.size() || nominal.rows.size() != excl_tight.rows.size() ||
        nominal.rows.size() != fid_loose.rows.size() || nominal.rows.size() != fid_tight.rows.size()) {
        throw std::runtime_error("BSA cut-systematic CSV row counts differ");
    } //endif

    fs::create_directories(fs::path(options.output_dir) / "bsa");
    std::ofstream diag(fs::path(options.output_dir) / "bsa/bsa_cut_variation_diagnostics.csv");
    diag << "group,row,A_nominal,stat_nominal,A_excl_loose,stat_excl_loose,A_excl_tight,stat_excl_tight,"
            "A_fid_loose,stat_fid_loose,A_fid_tight,stat_fid_tight,"
            "delta_excl_loose,delta_excl_tight,delta_fid_loose,delta_fid_tight,"
            "z_indep_excl_loose,z_indep_excl_tight,z_indep_fid_loose,z_indep_fid_tight,"
            "exclusivity_sys,fiducial_sys,total_cut_sys\n";

    const std::vector<std::string> groups = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out", "Fa18", "Sp18", "10.6 GeV"
    };
    for (const std::string& group : groups) {
        const std::string base = "BSA, counts, " + group;
        auto get_idx = [&](const CsvTable& t) -> std::size_t {
            auto it = t.index.find(base);
            if (it == t.index.end()) throw std::runtime_error("missing BSA column in variation CSV: " + base);
            return it->second;
        };
        const std::size_t in = get_idx(nominal);
        const std::size_t iel = get_idx(excl_loose);
        const std::size_t iet = get_idx(excl_tight);
        const std::size_t ifl = get_idx(fid_loose);
        const std::size_t ift = get_idx(fid_tight);
        const std::size_t cex = ensure_column(nominal, base + ", exclusivity sys");
        const std::size_t cfi = ensure_column(nominal, base + ", fiducial sys");
        const std::size_t ctot = ensure_column(nominal, base + ", cut sys");

        for (std::size_t r = 0; r < nominal.rows.size(); ++r) {
            const TripleCell n = parse_triple(nominal.rows[r][in]);
            const TripleCell el = parse_triple(excl_loose.rows[r][iel]);
            const TripleCell et = parse_triple(excl_tight.rows[r][iet]);
            const TripleCell fl = parse_triple(fid_loose.rows[r][ifl]);
            const TripleCell ft = parse_triple(fid_tight.rows[r][ift]);
            if (!n.ok || !el.ok || !et.ok || !fl.ok || !ft.ok) {
                nominal.rows[r][cex].clear(); nominal.rows[r][cfi].clear(); nominal.rows[r][ctot].clear();
                continue;
            } //endif
            // BSA crosses zero, so use absolute A_LU changes; no relative-difference
            // instability criterion is applied. The symmetric loose/tight average
            // mirrors the pass-1 cut prescription without dividing by A_LU.
            const double sex = 0.5 * (std::abs(el.value - n.value) + std::abs(et.value - n.value));
            const double sfi = 0.5 * (std::abs(fl.value - n.value) + std::abs(ft.value - n.value));
            const double stot = std::hypot(sex, sfi);
            nominal.rows[r][cex] = std::to_string(sex);
            nominal.rows[r][cfi] = std::to_string(sfi);
            nominal.rows[r][ctot] = std::to_string(stot);
            const double del = el.value - n.value;
            const double det = et.value - n.value;
            const double dfl = fl.value - n.value;
            const double dft = ft.value - n.value;
            auto z_independent = [](double delta, double s0, double s1) {
                // This deliberately treats the samples as independent, so it is
                // only a reference scale.  Nominal/varied samples overlap and
                // therefore have positive statistical covariance; the production
                // systematic remains the raw pass-1-style BSA displacement above.
                const double den = std::hypot(s0, s1);
                return den > 0.0 ? delta / den : 0.0;
            };
            diag << group << ',' << r << ',' << n.value << ',' << n.stat << ','
                 << el.value << ',' << el.stat << ',' << et.value << ',' << et.stat << ','
                 << fl.value << ',' << fl.stat << ',' << ft.value << ',' << ft.stat << ','
                 << del << ',' << det << ',' << dfl << ',' << dft << ','
                 << z_independent(del,n.stat,el.stat) << ',' << z_independent(det,n.stat,et.stat) << ','
                 << z_independent(dfl,n.stat,fl.stat) << ',' << z_independent(dft,n.stat,ft.stat) << ','
                 << sex << ',' << sfi << ',' << stot << '\n';
        } //endfor
    } //endfor
    write_csv(options.nominal_csv, nominal);
    std::cout << "[cut-variation-runner] Added absolute BSA exclusivity/fiducial cut systematics to "
              << options.nominal_csv << "\n";
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

enum class ExclusivityCutMode {
    ProductionTight90,
    ProductionLoose98,
    RefitNominal95
};

struct VariationSpec {
    std::string name;
    std::string csv_name;
    int fiducial_direction = 0; // -1 loose, 0 nominal, +1 tight
    ExclusivityCutMode exclusivity_mode =
        ExclusivityCutMode::RefitNominal95;
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


bool copy_file_checked(const fs::path& source,
                       const fs::path& destination,
                       const std::string& label) {
    std::error_code error;
    if (!fs::is_regular_file(source, error) || error) {
        throw std::runtime_error(
            "missing " + label + ": " + source.string());
    }

    fs::create_directories(destination.parent_path(), error);
    if (error) {
        throw std::runtime_error(
            "cannot create directory " + destination.parent_path().string() +
            ": " + error.message());
    }

    fs::copy_file(source, destination,
                  fs::copy_options::overwrite_existing, error);
    if (error) {
        throw std::runtime_error(
            "cannot copy " + label + " from " + source.string() +
            " to " + destination.string() + ": " + error.message());
    }
    return true;
}

std::string prepare_exclusivity_json(
    const VariationSpec& spec,
    const AutomaticCutVariationOptions& options,
    const GlobalCutConfig& cfg,
    const fs::path& base,
    const fs::path& json_dir) {

    set_default_global_cuts(cfg);
    write_global_cuts_config_json(json_dir.string(), cfg);

    if (spec.exclusivity_mode ==
        ExclusivityCutMode::ProductionTight90) {
        const fs::path source =
            fs::path(options.production_cuts_dir) /
            "combined_cuts_90.json";
        const fs::path destination =
            json_dir / "combined_cuts.json";
        copy_file_checked(source, destination,
                          "production 90% exclusivity JSON");
        return destination.string();
    }

    if (spec.exclusivity_mode ==
        ExclusivityCutMode::ProductionLoose98) {
        const fs::path source =
            fs::path(options.production_cuts_dir) /
            "combined_cuts_98.json";
        const fs::path destination =
            json_dir / "combined_cuts.json";
        copy_file_checked(source, destination,
                          "production 98% exclusivity JSON");
        return destination.string();
    }

    PythonExclusivityOptions python_options;
    python_options.enabled = true;
    python_options.force_rerun = true;
    python_options.python_executable = options.python_executable;
    python_options.script_path = options.python_script_path;
    python_options.global_cuts_json =
        (json_dir / "global_cuts_config.json").string();
    python_options.output_directory =
        (base / "python_exclusivity_fit").string();
    python_options.install_directory = json_dir.string();
    python_options.workers =
        std::max(1, std::min(options.max_workers, 7));
    python_options.tight_containment =
        options.tight_containment;
    python_options.nominal_containment =
        options.nominal_containment;
    python_options.loose_containment =
        options.loose_containment;

    if (!run_python_exclusivity_analysis(python_options)) {
        throw std::runtime_error(
            "Python exclusivity refit failed for variation " +
            spec.name);
    }

    const fs::path nominal_json =
        json_dir / "combined_cuts.json";
    if (!fs::is_regular_file(nominal_json)) {
        throw std::runtime_error(
            "Python exclusivity refit did not install " +
            nominal_json.string());
    }
    return nominal_json.string();
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

    using Clock = std::chrono::steady_clock;
    const auto variation_t0 = Clock::now();
    auto elapsed = [](const Clock::time_point& t0) {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };
    auto stage_start = [&](const std::string& name) {
        std::cout << "[cut-variation-runner]   [" << spec.name << "] "
                  << name << " started." << std::endl;
        return Clock::now();
    };
    auto stage_done = [&](const std::string& name, const Clock::time_point& t0) {
        std::cout << "[cut-variation-runner]   [" << spec.name << "] "
                  << name << " completed in " << std::fixed << std::setprecision(1)
                  << elapsed(t0) << " s." << std::defaultfloat << std::setprecision(6)
                  << std::endl;
    };

    const fs::path base = fs::path(options.output_dir) / "variations" / spec.name;
    const fs::path json_dir = base / "jsons";
    const fs::path stage_out = base / "analysis_output";
    const fs::path csv_dir = fs::path(options.output_dir) / "csv";
    const fs::path csv_path = csv_dir / spec.csv_name;
    fs::create_directories(json_dir);
    fs::create_directories(stage_out);
    fs::create_directories(csv_dir);
    fs::copy_file(options.nominal_csv, csv_path, fs::copy_options::overwrite_existing);

    const GlobalCutConfig cfg =
        varied_global_config(nominal_global_cuts,
                             spec.fiducial_direction);

    const auto cuts_t0 = stage_start("exclusivity-cut preparation/refit");
    const std::string cuts_json =
        prepare_exclusivity_json(spec, options, cfg, base, json_dir);
    stage_done("exclusivity-cut preparation/refit", cuts_t0);
    TotalCountsOptions count_opts;
    count_opts.use_nobkg_dvcs_mc_counts = use_nobkg_dvcs_mc_for_acceptance;
    count_opts.make_plots = false;
    count_opts.make_note_outputs = false;
    const auto counts_t0 = stage_start("total_counts");
    if (!update_total_counts_csv(csv_path.string(), dataTrees, eppi0DataTrees,
                                 genMcTrees, recMcTrees, eppi0GenMcTrees,
                                 eppi0RecMcTrees, eppi0BkgTrees, cuts_json,
                                 stage_out.string(), options.max_workers, count_opts,
                                 currentStudyGenMcTrees, currentStudyRecMcTrees)) return false;
    stage_done("total_counts", counts_t0);

    // Current-corrected DATA and reconstructed-MC columns are now filled
    // event-by-event by total_counts.cpp for the varied selection itself.
    // Do not rescale them from the nominal CSV; that would erase the regional
    // composition change induced by the cut variation.

    Pi0ContaminationOptions pi0_opts;
    pi0_opts.use_epg_mc_current_factor_for_eppi0_bkg = use_epg_mc_current_factor_for_eppi0_bkg;
    const auto pi0_t0 = stage_start("pi0 contamination");
    if (!compute_pi0_contamination_overall(dataTrees, eppi0DataTrees, eppi0RecMcTrees,
                                           eppi0BkgTrees, cuts_json, csv_path.string(),
                                           stage_out.string(), 1, pi0_opts)) return false;
    if (!update_pi0_corrected_counts_csv(csv_path.string(), stage_out.string())) return false;
    stage_done("pi0 contamination", pi0_t0);

    // BSA must be recomputed from the varied event selection itself.  This is
    // intentionally done before acceptance/radiative/bin-centering stages,
    // because those helicity-independent corrections are not part of A_LU.
    BSAOptions bsa_opts;
    bsa_opts.csv_path = csv_path.string();
    bsa_opts.combined_cuts_json = cuts_json;
    bsa_opts.output_root = (stage_out / "bsa").string();
    bsa_opts.enable_pi0_subtraction = true;
    bsa_opts.pi0_leakage_relative_uncertainty = 0.10;
    bsa_opts.make_plots = false;
    bsa_opts.make_photon_topology_study = false;
    bsa_opts.make_helicity_scrambling_study = false;
    bsa_opts.max_workers = options.max_workers;
    const auto bsa_t0 = stage_start("BSA");
    if (!update_bsa_counts_csv(dataTrees, eppi0DataTrees, bsa_opts)) return false;
    stage_done("BSA", bsa_t0);

    const auto acc_t0 = stage_start("acceptance");
    if (!update_acceptance_csv(csv_path.string(), genMcTrees, recMcTrees, cuts_json,
                               (json_dir / "global_cuts_config.json").string(),
                               stage_out.string())) return false;
    stage_done("acceptance", acc_t0);

    const auto unfold_t0 = stage_start("unfolding");
    if (!update_unfolded_yields_csv(csv_path.string(), stage_out.string())) return false;
    stage_done("unfolding", unfold_t0);

    LumiBuildOptions lumi_opts;
    lumi_opts.charge_csv_path =
        "imports/integrated_luminosity/global.csv";
    const LumiMap lumi = build_lumi_map(lumi_opts);
    const auto xsec_t0 = stage_start("cross sections");
    if (!compute_cross_sections(csv_path.string(), lumi)) return false;
    if (!update_normed_cross_sections_csv(csv_path.string())) return false;
    stage_done("cross sections", xsec_t0);

    std::cout << "[cut-variation-runner] Completed " << spec.name
              << " in " << std::fixed << std::setprecision(1)
              << elapsed(variation_t0) << " s -> " << csv_path.string()
              << std::defaultfloat << std::setprecision(6) << std::endl;
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
        if (!(options.tight_containment > 0.0 &&
              options.tight_containment <
                  options.nominal_containment &&
              options.nominal_containment <
                  options.loose_containment &&
              options.loose_containment < 1.0)) {
            throw std::runtime_error(
                "containments must satisfy 0 < tight < nominal < loose < 1");
        }
        if (options.max_workers < 1) {
            throw std::runtime_error(
                "max_workers must be at least 1");
        }

        fs::create_directories(options.output_dir);
        fs::copy_file(options.nominal_csv,
                      fs::path(options.output_dir) / "nominal_csv_before_cut_systematics.csv",
                      fs::copy_options::overwrite_existing);

        const std::vector<VariationSpec> variations = {
            {"exclusivity_loose_98",
             "exclusivity_loose_98.csv",
             0,
             ExclusivityCutMode::ProductionLoose98},
            {"exclusivity_tight_90",
             "exclusivity_tight_90.csv",
             0,
             ExclusivityCutMode::ProductionTight90},
            {"fiducial_loose",
             "fiducial_loose.csv",
             -1,
             ExclusivityCutMode::RefitNominal95},
            {"fiducial_tight",
             "fiducial_tight.csv",
             +1,
             ExclusivityCutMode::RefitNominal95}
        };

        for (size_t iv = 0; iv < variations.size(); ++iv) {
            const auto& v = variations[iv];
            std::cout << "[cut-variation-runner] Variation " << (iv + 1)
                      << "/" << variations.size() << ": " << v.name
                      << std::endl;
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

        // BSA-specific cut systematics use absolute A_LU differences because
        // the asymmetry crosses zero. They are intentionally kept separate from
        // the cross-section Barlow/relative-systematic machinery below.
        update_bsa_cut_systematics(options);

        CutVariationSystematicsOptions syst;
        syst.enabled = true;
        syst.apply_barlow = true;
        syst.barlow_threshold = 1.0;
        syst.use_pass1_tight_instability_rule =
            options.use_pass1_tight_instability_rule;
        syst.tight_relative_difference_threshold =
            options.tight_relative_difference_threshold;
        syst.make_plots = options.make_final_diagnostic_plots;
        syst.write_diagnostic_csv = true;
        syst.nominal_csv = options.nominal_csv;
        syst.exclusivity_loose_csv =
            (fs::path(options.output_dir) /
             "csv/exclusivity_loose_98.csv").string();
        syst.exclusivity_tight_csv =
            (fs::path(options.output_dir) /
             "csv/exclusivity_tight_90.csv").string();
        syst.fiducial_loose_csv = (fs::path(options.output_dir) / "csv/fiducial_loose.csv").string();
        syst.fiducial_tight_csv = (fs::path(options.output_dir) / "csv/fiducial_tight.csv").string();
        syst.output_dir = options.output_dir;
        std::cout << "[cut-variation-runner] All " << variations.size()
                  << " variation CSVs complete; building Barlow/systematic summary."
                  << std::endl;
        const auto syst_t0 = std::chrono::steady_clock::now();
        const bool ok = update_cut_variation_systematics(syst);
        const double syst_sec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - syst_t0).count();
        std::cout << "[cut-variation-runner] Barlow/systematic summary "
                  << (ok ? "completed" : "failed")
                  << " in " << std::fixed << std::setprecision(1)
                  << syst_sec << " s." << std::defaultfloat << std::setprecision(6)
                  << std::endl;
        return ok;
    } catch (const std::exception& e) {
        set_default_global_cuts(nominal_global_cuts);
        std::cerr << "[cut-variation-runner] ERROR: " << e.what() << std::endl;
        return false;
    }
}
