// main_systematics.cpp
// -----------------------------------------------------------------------------
// Systematic-uncertainty driver for the DVCS pass-2 CSV.
//
// This executable is intentionally CSV-only.  It does not load ROOT trees.  The
// individual systematic studies added here should read the quantities already
// written to output/csvs/dvcs_pass2_analysis.csv, update the appropriate CSV
// columns, and fail fast if required input or output columns are missing.
//
// Current shell behavior:
//   - create output/systematics;
//   - load the pass-2 CSV header and rows;
//   - verify that the normed combined cross-section columns exist;
//   - verify that the newly added combination-systematic columns exist;
//   - make a backup copy of the CSV before future update stages are added.
//
// Future stages should be added as explicit function calls in main(), following
// the same pattern used by main.cpp and cross_check_lee_main.cpp.
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "combination_systematics.h"
#include "correlated_scale_systematics.h"
#include "pass1_systematics_import.h"
#include "pi0_systematics.h"
#include "combination_point_to_point_systematics.h"
#include "sp19_inb_energy_scaling_systematics.h"
#include "run_period_consistency_systematics.h"
#include "systematic_projection_plots.h"
#include "radiative_systematic_plots.h"
#include "bin_centering_systematic_plots.h"

namespace fs = std::filesystem;

namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rows;
    std::unordered_map<std::string, int> index;
};

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());

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
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    out.push_back(cur);
    return out;
}

static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream oss;

    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];

        const bool need_quotes =
            (s.find(',') != std::string::npos) ||
            (s.find('"') != std::string::npos) ||
            (s.find('\n') != std::string::npos) ||
            (s.find('\r') != std::string::npos);

        if (need_quotes) {
            oss << '"';

            for (const char ch : s) {
                if (ch == '"') {
                    oss << "\"\"";
                } else {
                    oss << ch;
                }
            }

            oss << '"';
        } else {
            oss << s;
        }

        if (i + 1 < fields.size()) {
            oss << ',';
        }
    }

    return oss.str();
}

static std::unordered_map<std::string, int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> idx;

    for (int i = 0; i < (int)header.size(); ++i) {
        idx[header[(size_t)i]] = i;
    }

    return idx;
}

static void assert_no_duplicate_columns(const std::vector<std::string>& header) {
    std::set<std::string> seen;
    std::set<std::string> duplicates;

    for (const auto& name : header) {
        if (!seen.insert(name).second) {
            duplicates.insert(name);
        }
    }

    if (!duplicates.empty()) {
        std::ostringstream msg;
        msg << "CSV header contains duplicate columns:";

        for (const auto& name : duplicates) {
            msg << "\n  - " << name;
        }

        throw std::runtime_error(msg.str());
    }
}

static CsvTable read_csv_or_throw(const std::string& path) {
    std::ifstream fin(path);

    if (!fin.is_open()) {
        throw std::runtime_error("Could not open CSV: " + path);
    }

    CsvTable table;

    std::string line;
    if (!std::getline(fin, line)) {
        throw std::runtime_error("CSV is empty: " + path);
    }

    table.header = split_csv_line(line);
    assert_no_duplicate_columns(table.header);
    table.index = build_header_index(table.header);

    while (std::getline(fin, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> row = split_csv_line(line);

        if (row.size() < table.header.size()) {
            row.resize(table.header.size());
        }

        if (row.size() != table.header.size()) {
            std::ostringstream msg;
            msg << "CSV row has " << row.size()
                << " columns, but header has " << table.header.size()
                << " columns. File: " << path;
            throw std::runtime_error(msg.str());
        }

        table.rows.push_back(std::move(row));
    }

    return table;
}

static void write_csv_or_throw(const std::string& path, const CsvTable& table) {
    const std::string tmp_path = path + ".tmp";
    std::ofstream fout(tmp_path);

    if (!fout.is_open()) {
        throw std::runtime_error("Could not open temporary CSV for writing: " + tmp_path);
    }

    fout << join_csv_row(table.header) << "\n";

    for (const auto& row : table.rows) {
        if (row.size() != table.header.size()) {
            throw std::runtime_error("Internal CSV row/header size mismatch before writing.");
        }

        fout << join_csv_row(row) << "\n";
    }

    fout.close();

    if (!fout) {
        throw std::runtime_error("Failed while writing temporary CSV: " + tmp_path);
    }

    fs::rename(tmp_path, path);
}

static void require_columns(const CsvTable& table,
                            const std::vector<std::string>& required,
                            const std::string& context) {
    std::vector<std::string> missing;

    for (const auto& name : required) {
        if (table.index.find(name) == table.index.end()) {
            missing.push_back(name);
        }
    }

    if (!missing.empty()) {
        std::ostringstream msg;
        msg << "Missing required columns for " << context << ":";

        for (const auto& name : missing) {
            msg << "\n  - " << name;
        }

        throw std::runtime_error(msg.str());
    }
}

static std::vector<std::string> combined_normed_cross_section_columns() {
    return {
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol",
        "normed cross sections, ep->epg, exp, Fa18, unpol",
        "normed cross sections, ep->epg, exp, Fa18, pos",
        "normed cross sections, ep->epg, exp, Fa18, neg",
        "normed cross sections, ep->epg, exp, Sp18, unpol",
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol"
    };
}

static std::vector<std::string> combination_systematic_columns() {
    std::vector<std::string> out;
    const std::vector<std::pair<std::string, std::string> > xs_targets = {
        {"10.6 GeV", "unpol"}, {"Fa18", "unpol"}, {"Fa18", "pos"},
        {"Fa18", "neg"}, {"Sp18", "unpol"}, {"Sp19 Inb", "unpol"}
    };
    for (const auto& target : xs_targets) {
        const std::string base = "normed cross sections, ep->epg, exp, " + target.first + ", " + target.second;
        out.push_back(base + ", combination sys");
        out.push_back(base + ", target thickness and charge sys");
        out.push_back(base + ", total scale sys");
    }
    for (const auto& label : std::vector<std::string>{"10.6 GeV", "Fa18", "Sp18", "Sp19 Inb"}) {
        const std::string base = "BSA, counts, " + label;
        out.push_back(base + ", combination sys");
        out.push_back(base + ", beam polarization sys");
        out.push_back(base + ", total scale sys");
    }
    return out;
}

static std::vector<std::string> pass1_fixed_systematic_columns() {
    return {
        "Syst. err (pi0 subtraction)",
        "Syst. err (Acceptance)",
        "Syst.err (Frad)",
        "Syst.err (Fbin)",
        "Syst. err (point-to-point total)"
    };
}

static void validate_combination_systematics_schema(const CsvTable& table) {
    require_columns(table,
                    combined_normed_cross_section_columns(),
                    "combined normed cross-section systematics inputs");

    require_columns(table,
                    combination_systematic_columns(),
                    "combination systematic outputs");

    require_columns(table,
                    pass1_fixed_systematic_columns(),
                    "pass-1 fixed systematic outputs");
}

static void ensure_column(CsvTable& table,
                          const std::string& column,
                          const std::string& initial_value = "") {
    if (table.index.find(column) != table.index.end()) {
        return;
    }

    table.index[column] = (int)table.header.size();
    table.header.push_back(column);

    for (auto& row : table.rows) {
        row.push_back(initial_value);
    }
}

static bool ensure_systematics_output_columns(const std::string& csv_path) {
    CsvTable table = read_csv_or_throw(csv_path);

    size_t n_added = 0;
    auto add_if_missing = [&](const std::string& col) {
        if (table.index.find(col) == table.index.end()) {
            ensure_column(table, col);
            ++n_added;
        }
    };

    for (const auto& col : pass1_fixed_systematic_columns()) {
        add_if_missing(col);
    }
    for (const auto& col : combination_systematic_columns()) {
        add_if_missing(col);
    }

    if (n_added > 0) {
        write_csv_or_throw(csv_path, table);
        std::cout << "[systematics] Added " << n_added
                  << " missing systematics output columns to " << csv_path << "\n";
    }

    return n_added > 0;
}

static std::string trim_copy(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;
    return s.substr(b, e - b);
}

static double scalar_value(const std::string& raw) {
    const std::string s = trim_copy(raw);
    if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* end = nullptr;
    const double v = std::strtod(s.c_str(), &end);
    return (end == s.c_str()) ? std::numeric_limits<double>::quiet_NaN() : v;
}

static double tuple_first_value(const std::string& raw) {
    std::string s = trim_copy(raw);
    while (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        s = trim_copy(s.substr(1, s.size() - 2));
    if (s.size() < 3 || s.front() != '(' || s.back() != ')')
        return std::numeric_limits<double>::quiet_NaN();
    s = s.substr(1, s.size() - 2);
    const size_t comma = s.find(',');
    return scalar_value(comma == std::string::npos ? s : s.substr(0, comma));
}

static std::string format_scalar(double v) {
    if (!std::isfinite(v)) return std::string();
    std::ostringstream ss;
    ss << std::setprecision(12) << v;
    return ss.str();
}


static double quantile_copy(std::vector<double> v, double q) {
    v.erase(std::remove_if(v.begin(), v.end(),
                           [](double x){ return !std::isfinite(x); }), v.end());
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    if (q <= 0.0) return v.front();
    if (q >= 1.0) return v.back();
    const double pos = q * (double)(v.size() - 1);
    const size_t lo = (size_t)std::floor(pos);
    const size_t hi = (size_t)std::ceil(pos);
    const double f = pos - (double)lo;
    return (1.0 - f) * v[lo] + f * v[hi];
}

// Recover the combined 10.6-GeV acceptance candidate from the already-written
// acceptance-reweighting diagnostics.  This is intentionally CSV-only: if an
// expensive acceptance run finished but its final 10.6-GeV columns were left
// blank (for example because an older CSV contained an extra wrapper-quote
// layer around tuple-valued yield cells), we can reconstruct exactly the same
// combination without rerunning any ROOT event loops.
static bool recover_acceptance_reweighting_from_diagnostics(CsvTable& t) {
    const std::string summary_path =
        "output/systematics/acceptance_reweighting/acceptance_reweighting_per_period_summary.csv";
    const std::string closure_path =
        "output/systematics/acceptance_reweighting/synthetic_reweighting_closure.csv";

    if (!fs::exists(summary_path) || !fs::exists(closure_path)) return false;

    const CsvTable ps = read_csv_or_throw(summary_path);
    const CsvTable sc = read_csv_or_throw(closure_path);

    require_columns(ps, {"period","row","acceptance_nominal","acceptance_data_reweighted"},
                    "acceptance-reweighting period summary recovery");
    require_columns(sc, {"period","row","closure_bias_frac"},
                    "acceptance-reweighting closure recovery");

    const std::array<std::string,4> p10 = {{
        "Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"
    }};

    struct Pair { double a0=std::numeric_limits<double>::quiet_NaN();
                  double ad=std::numeric_limits<double>::quiet_NaN(); };
    std::map<std::string,std::unordered_map<size_t,Pair>> acc;
    std::map<std::string,std::unordered_map<size_t,std::vector<double>>> biases;

    const int ips_per=ps.index.at("period");
    const int ips_row=ps.index.at("row");
    const int ips_a0=ps.index.at("acceptance_nominal");
    const int ips_ad=ps.index.at("acceptance_data_reweighted");
    for (const auto& r: ps.rows) {
        const std::string per=trim_copy(r[(size_t)ips_per]);
        const double rr=scalar_value(r[(size_t)ips_row]);
        if (!std::isfinite(rr) || rr < 0.0) continue;
        Pair x;
        x.a0=scalar_value(r[(size_t)ips_a0]);
        x.ad=scalar_value(r[(size_t)ips_ad]);
        if (std::isfinite(x.a0) && x.a0>0.0 && std::isfinite(x.ad) && x.ad>0.0)
            acc[per][(size_t)std::llround(rr)]=x;
    }

    const int isc_per=sc.index.at("period");
    const int isc_row=sc.index.at("row");
    const int isc_b=sc.index.at("closure_bias_frac");
    for (const auto& r: sc.rows) {
        const std::string per=trim_copy(r[(size_t)isc_per]);
        const double rr=scalar_value(r[(size_t)isc_row]);
        const double b=scalar_value(r[(size_t)isc_b]);
        if (!std::isfinite(rr) || rr<0.0 || !std::isfinite(b)) continue;
        biases[per][(size_t)std::llround(rr)].push_back(b);
    }

    const int c_data10 = t.index.at("acceptance reweighting data-driven sys frac, 10.6 GeV");
    const int c_transfer10 = t.index.at("acceptance reweighting transfer closure sys frac, 10.6 GeV");
    const int c_cons10 = t.index.at("acceptance reweighting conservative sys frac, 10.6 GeV");
    const int c_cand10 = t.index.at("acceptance reweighting candidate sys frac, 10.6 GeV");

    size_t n=0;
    for (size_t i=0;i<t.rows.size();++i) {
        // Preserve a valid value from the tree-based stage.
        if (std::isfinite(scalar_value(t.rows[i][(size_t)c_cons10]))) continue;

        double ytot=0.0, weighted_ratio=0.0, trmax=0.0;
        bool have_transfer=false;
        for (const auto& per:p10) {
            auto ia=acc[per].find(i);
            if (ia==acc[per].end()) continue;
            const std::string yc="acceptance corrected yield, ep->epg, exp, "+per+", unpol";
            auto iy=t.index.find(yc);
            if (iy==t.index.end()) continue;
            const double y=tuple_first_value(t.rows[i][(size_t)iy->second]);
            if (!std::isfinite(y) || y<0.0) continue;
            ytot += y;
            weighted_ratio += y*(ia->second.a0/ia->second.ad);
            auto ibp=biases.find(per);
            if (ibp!=biases.end()) {
                auto ib=ibp->second.find(i);
                if (ib!=ibp->second.end() && !ib->second.empty()) {
                    const double tr=quantile_copy(ib->second,0.95);
                    if (std::isfinite(tr)) { trmax=std::max(trmax,tr); have_transfer=true; }
                }
            }
        }
        if (ytot<=0.0) continue;
        const double ddata=std::fabs(weighted_ratio/ytot - 1.0);
        const double tr=have_transfer ? trmax : 0.0;
        const double cons=std::hypot(ddata,tr);
        t.rows[i][(size_t)c_data10]=format_scalar(ddata);
        t.rows[i][(size_t)c_transfer10]=format_scalar(tr);
        t.rows[i][(size_t)c_cons10]=format_scalar(cons);
        t.rows[i][(size_t)c_cand10]=format_scalar(cons);
        ++n;
    }
    if (n>0)
        std::cout << "[acceptance-systematics] Recovered " << n
                  << " combined 10.6-GeV acceptance candidates from the completed "
                  << "acceptance diagnostic files; no event-loop rerun was required.\n";
    return n>0;
}

// Promote the reviewed acceptance-reweighting candidate to the production
// point-to-point systematic.  The tree-based study writes fractional
// uncertainties for 10.6 GeV and Sp19 separately.  Here we convert those
// fractions to absolute cross-section uncertainties and rebuild the total
// point-to-point uncertainty with the final pass-2 components.
static bool install_acceptance_reweighting_systematic(const std::string& csv_path) {
    CsvTable t = read_csv_or_throw(csv_path);

    recover_acceptance_reweighting_from_diagnostics(t);

    const std::string c_frac10 =
        "acceptance reweighting conservative sys frac, 10.6 GeV";
    const std::string c_fracsp =
        "acceptance reweighting conservative sys frac, Sp19 Inb";
    const std::string c_xs10 =
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol";
    const std::string c_xssp =
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol";

    require_columns(t, {c_frac10, c_fracsp, c_xs10, c_xssp,
                       "Syst. err (pi0 subtraction)",
                       "Syst. err (Acceptance)",
                       "Syst.err (Frad)",
                       "Syst.err (Fbin)",
                       "Syst. err (exclusivity cuts)",
                       "Syst. err (fiducial cuts)",
                       "Syst. err (point-to-point total)",
                       "Syst. err (pi0 subtraction), Sp19 Inb (10.2 GeV)",
                       "Syst.err (Frad), Sp19 Inb (10.2 GeV)"},
                    "acceptance-reweighting production assignment");

    ensure_column(t, "Syst. err (Acceptance), Sp19 Inb (10.2 GeV)");
    ensure_column(t, "Syst. err (point-to-point total), Sp19 Inb (10.2 GeV)");

    const int i_frac10 = t.index.at(c_frac10);
    const int i_fracsp = t.index.at(c_fracsp);
    const int i_xs10 = t.index.at(c_xs10);
    const int i_xssp = t.index.at(c_xssp);
    const int i_acc10 = t.index.at("Syst. err (Acceptance)");
    const int i_accsp = t.index.at("Syst. err (Acceptance), Sp19 Inb (10.2 GeV)");
    const int i_ptp10 = t.index.at("Syst. err (point-to-point total)");
    const int i_ptpsp = t.index.at("Syst. err (point-to-point total), Sp19 Inb (10.2 GeV)");

    const std::vector<std::string> components10 = {
        "Syst. err (pi0 subtraction)",
        "Syst. err (Acceptance)",
        "Syst.err (Frad)",
        "Syst.err (Fbin)",
        "Syst. err (exclusivity cuts)",
        "Syst. err (fiducial cuts)"
    };

    size_t n10 = 0, nsp = 0;
    std::vector<double> frac10_values, fracsp_values;

    for (auto& row : t.rows) {
        const double xs10 = tuple_first_value(row[(size_t)i_xs10]);
        const double xssp = tuple_first_value(row[(size_t)i_xssp]);
        const double f10 = scalar_value(row[(size_t)i_frac10]);
        const double fsp = scalar_value(row[(size_t)i_fracsp]);

        if (std::isfinite(xs10) && std::isfinite(f10) && f10 >= 0.0) {
            row[(size_t)i_acc10] = format_scalar(std::fabs(xs10) * f10);
            frac10_values.push_back(f10);
            ++n10;
        } else {
            row[(size_t)i_acc10].clear();
        }

        if (std::isfinite(xssp) && std::isfinite(fsp) && fsp >= 0.0) {
            row[(size_t)i_accsp] = format_scalar(std::fabs(xssp) * fsp);
            fracsp_values.push_back(fsp);
            ++nsp;
        } else {
            row[(size_t)i_accsp].clear();
        }

        // Final 10.6-GeV point-to-point total from the six production terms.
        double sum10 = 0.0;
        bool ok10 = true;
        for (const auto& col : components10) {
            const double e = scalar_value(row[(size_t)t.index.at(col)]);
            if (!std::isfinite(e) || e < 0.0) {
                ok10 = false;
                break;
            }
            sum10 += e * e;
        }
        row[(size_t)i_ptp10] = ok10 ? format_scalar(std::sqrt(sum10)) : std::string();

        // Sp19 has dedicated pi0, acceptance and radiative terms.  Fbin,
        // exclusivity and fiducial are transferred using their 10.6-GeV
        // bin-wise fractions, matching the existing pass-2 prescription.
        bool oksp = std::isfinite(xssp) && std::isfinite(xs10) && std::fabs(xs10) > 0.0;
        double sumsp = 0.0;
        if (oksp) {
            for (const auto& col : std::vector<std::string>{
                    "Syst. err (pi0 subtraction), Sp19 Inb (10.2 GeV)",
                    "Syst. err (Acceptance), Sp19 Inb (10.2 GeV)",
                    "Syst.err (Frad), Sp19 Inb (10.2 GeV)"}) {
                const double e = scalar_value(row[(size_t)t.index.at(col)]);
                if (!std::isfinite(e) || e < 0.0) { oksp = false; break; }
                sumsp += e * e;
            }
        }
        if (oksp) {
            for (const auto& col : std::vector<std::string>{
                    "Syst.err (Fbin)",
                    "Syst. err (exclusivity cuts)",
                    "Syst. err (fiducial cuts)"}) {
                const double e10 = scalar_value(row[(size_t)t.index.at(col)]);
                if (!std::isfinite(e10) || e10 < 0.0) { oksp = false; break; }
                const double frac = std::fabs(e10 / xs10);
                const double esp = std::fabs(xssp) * frac;
                sumsp += esp * esp;
            }
        }
        row[(size_t)i_ptpsp] = oksp ? format_scalar(std::sqrt(sumsp)) : std::string();
    }

    write_csv_or_throw(csv_path, t);

    auto median_fraction = [](std::vector<double> v) {
        v.erase(std::remove_if(v.begin(), v.end(),
                               [](double x){ return !std::isfinite(x); }), v.end());
        if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
        std::sort(v.begin(), v.end());
        const size_t n = v.size();
        return n % 2 ? v[n/2] : 0.5*(v[n/2 - 1] + v[n/2]);
    };

    std::cout << "[acceptance-systematics] Installed conservative pass-2 acceptance "
              << "uncertainty in " << n10 << " 10.6-GeV bins and " << nsp
              << " Sp19 bins. Median fractions: "
              << 100.0*median_fraction(frac10_values) << "% (10.6), "
              << 100.0*median_fraction(fracsp_values) << "% (Sp19).\n";
    std::cout << "[acceptance-systematics] Recomputed dedicated 10.6-GeV and Sp19 "
              << "point-to-point totals.\n";
    return true;
}

static void make_output_dirs() {
    fs::create_directories("output");
    fs::create_directories("output/csvs");
    fs::create_directories("output/systematics");
}

static void backup_csv_or_throw(const std::string& csv_path,
                                const std::string& backup_path) {
    fs::copy_file(csv_path,
                  backup_path,
                  fs::copy_options::overwrite_existing);
}

static bool validate_systematics_csv_schema(const std::string& csv_path) {
    CsvTable table = read_csv_or_throw(csv_path);

    validate_combination_systematics_schema(table);

    std::cout << "[systematics] CSV rows loaded: " << table.rows.size() << "\n";
    std::cout << "[systematics] CSV columns loaded: " << table.header.size() << "\n";
    std::cout << "[systematics] Schema validation passed.\n";

    return true;
}

} // namespace

int main(int argc, char* argv[]) {
    std::cout << "[systematics] main_systematics starting...\n";

    const std::string csv_main =
        (argc >= 2) ? std::string(argv[1]) : std::string("output/csvs/dvcs_pass2_analysis.csv");

    try {
        make_output_dirs();

        const std::string backup_path =
            "output/csvs/dvcs_pass2_analysis_backup_systematics.csv";

        backup_csv_or_throw(csv_main, backup_path);
        std::cout << "[systematics] Backed up CSV to " << backup_path << "\n";

        ensure_systematics_output_columns(csv_main);

        if (!validate_systematics_csv_schema(csv_main)) {
            std::cerr << "[systematics] FATAL: validate_systematics_csv_schema failed.\n";
            return 1;
        }

        const std::string pass1_systematics_path =
            (argc >= 3) ? std::string(argv[2]) : std::string("imports/pass1_systematic_summary.csv");

        if (!import_pass1_systematics(csv_main, pass1_systematics_path)) {
            std::cerr << "[systematics] FATAL: import_pass1_systematics failed.\n";
            return 1;
        }

        // Recalculate the pi0-subtraction systematic from the pass-2
        // contamination itself.  The 7.2% relative uncertainty on the
        // acceptance-corrected background is inherited from the dedicated
        // pass-1 variation study, while the actual background fraction is
        // determined from the present pass-2 yields.
        if (!pi0_systematics(
                csv_main,
                pass1_systematics_path,
                "imports/all_bin_v3.csv",
                "output/systematics/pi0_systematics")) {
            std::cerr << "[systematics] FATAL: pi0_systematics failed.\n";
            return 1;
        }

        if (!install_acceptance_reweighting_systematic(csv_main)) {
            std::cerr << "[systematics] FATAL: install_acceptance_reweighting_systematic failed.\n";
            return 1;
        }

        if (!combination_systematics(csv_main, "output/systematics")) {
            std::cerr << "[systematics] FATAL: combination_systematics failed.\n";
            return 1;
        }

        // Build the authoritative final correlated-scale category only after
        // the current-systematic columns and the legacy run-period diagnostics
        // are available.  The charge/target-thickness normalization is
        // kept separate and is not folded into this nuisance.
        CorrelatedScaleSystematicsOptions correlated_opts;
        correlated_opts.theta_bin_width_deg = 4.0;
        correlated_opts.min_ratio_points_per_period = 25;
        correlated_opts.uncorrelated_normalization_fraction = 0.021633307652784; // sqrt(0.012^2 + 0.018^2) = 2.16%
        correlated_opts.output_dir =
            "output/systematics/correlated_scale";

        if (!correlated_scale_systematics(csv_main, correlated_opts)) {
            std::cerr
                << "[systematics] FATAL: correlated_scale_systematics failed.\n";
            return 1;
        }

        if (!make_radiative_systematic_analysis_note_plots(
                csv_main,
                "output/systematics/analysis_note/radiative_corrections")) {
            std::cerr << "[systematics] FATAL: "
                      << "make_radiative_systematic_analysis_note_plots failed.\n";
            return 1;
        }

        // Duplicate the radiative-systematic figures into the familiar
        // radiative-corrections analysis-note directory as well as the
        // canonical output/systematics tree.
        make_radiative_systematic_analysis_note_plots(
            csv_main,
            "output/radiative_corrections/analysis_note");

        if (!make_bin_centering_systematic_analysis_note_plots(
                pass1_systematics_path,
                "imports/all_bin_v3.csv",
                "output/systematics/analysis_note/bin_centering")) {
            std::cerr << "[systematics] FATAL: "
                      << "make_bin_centering_systematic_analysis_note_plots failed.\n";
            return 1;
        }

        // Likewise keep the bin-centering systematic figures beside the
        // bin-centering correction figures for analysis-note packaging.
        make_bin_centering_systematic_analysis_note_plots(
            pass1_systematics_path,
            "imports/all_bin_v3.csv",
            "output/bin_centering_plots/analysis_note");

        if (!sp19_inb_energy_scaling_systematics(csv_main, "output/systematics")) {
            std::cerr << "[systematics] FATAL: sp19_inb_energy_scaling_systematics failed.\n";
            return 1;
        }

        if (!combination_point_to_point_systematics(csv_main, "output/systematics")) {
            std::cerr << "[systematics] FATAL: combination_point_to_point_systematics failed.\n";
            return 1;
        }

        if (!make_systematic_projection_plots(
                csv_main,
                "output/systematics/point_to_point_projections")) {
            std::cerr << "[systematics] FATAL: make_systematic_projection_plots failed.\n";
            return 1;
        }

        // if (!run_period_consistency_systematics(csv_main, "output/systematics", false)) {
        //     std::cerr << "[systematics] FATAL: run_period_consistency_systematics failed.\n";
        //     return 1;
        // }
        
    } catch (const std::exception& e) {
        std::cerr << "[systematics] FATAL: " << e.what() << "\n";
        return 1;
    }

    std::cout << "[systematics] main_systematics complete.\n";
    return 0;
}
