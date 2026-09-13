// photon_efficiency_preliminary_correction.cpp
//
// Extremely preliminary, stand-alone post-processing test for the DVCS pass-2
// CSV.  This is intentionally NOT the production photon-efficiency correction.
//
// Purpose:
//   Apply the current integrated photon data/MC efficiency ratios as a simple
//   topology- and photon-energy-weighted multiplicative correction to an
//   already-produced pass-2 CSV, so the corrected cross sections can be fed
//   immediately to compare_world_dvcs_cross_sections.py.
//
// Current convention:
//     R = epsilon_data / epsilon_MC
//     sigma_corrected = sigma_current * (epsilon_MC / epsilon_data) = sigma/R
//
// Photon detector by topology:
//     (FD, FD) -> photon FD
//     (CD, FD) -> photon FD
//     (CD, FT) -> photon FT
//
// Photon energy:
//   The tag-and-probe correction is parameterized versus the PROBE photon
//   energy.  For this exploratory application to a measured DVCS photon, use
//   the exclusive-photon lab energy reconstructed from the mean bin kinematics:
//
//       E_gamma = nu - |t|/(2 Mp)
//               = Q2/(2 Mp xB) - |t|/(2 Mp)
//
//   This avoids requiring an event-level photon-energy column in the final CSV.
//
// IMPORTANT APPROXIMATION:
//   This program acts POST HOC on the final cross section.  It estimates an
//   effective correction from the topology composition of the current-corrected
//   selected yields.  The production implementation should eventually apply
//   the photon correction upstream at event/yield level before pi0 subtraction
//   and acceptance unfolding.
//
// Build:
//   g++ -O2 -std=c++17 photon_efficiency_preliminary_correction.cpp \
//       -o photon_efficiency_preliminary_correction
//
// Typical exploratory use:
//   ./photon_efficiency_preliminary_correction
//
// Then compare the corrected combined 10.6-GeV CSV:
//   python external_scripts/compare_world_dvcs_cross_sections.py \
//       --pass2-file output/csvs/dvcs_pass2_analysis_photon_eff_prelim.csv \
//       --outdir ../output/world_dvcs_cross_section_comparison_photon_eff_prelim
//
// The constants below are deliberately centralized so later detector-, period-,
// sector-, theta-, and energy-dependent maps can replace them cleanly.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Configuration: CURRENT EXTREMELY PRELIMINARY integrated ratios.
// Update only this block for quick future tests.
// -----------------------------------------------------------------------------

constexpr double MP_GEV = 0.9382720813;
constexpr double ENERGY_SPLIT_GEV = 2.0;

// R = epsilon_data / epsilon_MC.
//
// FD values: current standard / validated core method.
// FT values: current robust sparse-bin Poisson/core estimate.
//
// These should be treated as exploratory only.
struct EfficiencyRatios {
    double fd_low  = 1.0304;
    double fd_high = 0.7024;
    double ft_low  = 0.5927;
    double ft_high = 0.3506;
};

const EfficiencyRatios EFF{};

// Default paths.
const std::string DEFAULT_INPUT =
    "output/csvs/dvcs_pass2_analysis.csv";
const std::string DEFAULT_OUTPUT =
    "output/csvs/dvcs_pass2_analysis_photon_eff_prelim.csv";

// Correct both the requested Fa18 Inb measurement and the combined 10.6-GeV
// quantity used by compare_world_dvcs_cross_sections.py.
const std::vector<std::string> DEFAULT_LABELS = {
    "Fa18 Inb",
    "10.6 GeV"
};

// Period membership for combined labels.  The same preliminary efficiency map
// is deliberately transferred to all periods in this quick sensitivity test.
const std::map<std::string, std::vector<std::string>> LABEL_MEMBERS = {
    {"Fa18 Inb", {"Fa18 Inb"}},
    {"Fa18", {"Fa18 Inb", "Fa18 Out"}},
    {"10.6 GeV", {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"}}
};

// -----------------------------------------------------------------------------
// Minimal robust CSV utilities.
// -----------------------------------------------------------------------------

std::vector<std::string> parse_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool quoted = false;

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (quoted) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    cur.push_back('"');
                    ++i;
                } else {
                    quoted = false;
                }
            } else {
                cur.push_back(c);
            }
        } else {
            if (c == '"') {
                quoted = true;
            } else if (c == ',') {
                out.push_back(cur);
                cur.clear();
            } else {
                cur.push_back(c);
            }
        }
    }

    out.push_back(cur);
    return out;
}

std::string csv_escape(const std::string& s) {
    const bool need =
        s.find(',') != std::string::npos ||
        s.find('"') != std::string::npos ||
        s.find('\n') != std::string::npos ||
        s.find('\r') != std::string::npos;

    if (!need) return s;

    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out.push_back(c);
    }
    out += '"';
    return out;
}

void write_csv_row(std::ostream& os, const std::vector<std::string>& row) {
    for (std::size_t i = 0; i < row.size(); ++i) {
        if (i) os << ',';
        os << csv_escape(row[i]);
    }
    os << '\n';
}

std::map<std::string, std::size_t>
make_header_map(const std::vector<std::string>& header) {
    std::map<std::string, std::size_t> m;
    for (std::size_t i = 0; i < header.size(); ++i)
        m[header[i]] = i;
    return m;
}

double parse_double(const std::string& s, double fallback = 0.0) {
    if (s.empty()) return fallback;
    try {
        std::size_t n = 0;
        const double x = std::stod(s, &n);
        if (n == 0 || !std::isfinite(x)) return fallback;
        return x;
    } catch (...) {
        return fallback;
    }
}

std::string fmt(double x) {
    std::ostringstream ss;
    ss << std::setprecision(12) << x;
    return ss.str();
}

// The cross-section CSV stores entries such as "(12.34, 0.56)".
// Scale both central value and statistical uncertainty.
bool scale_tuple_cell(std::string& cell, double factor) {
    std::string s = cell;
    s.erase(std::remove_if(s.begin(), s.end(),
                           [](unsigned char c){ return std::isspace(c); }),
            s.end());

    if (s.size() < 5 || s.front() != '(' || s.back() != ')')
        return false;

    s = s.substr(1, s.size() - 2);
    const auto comma = s.find(',');
    if (comma == std::string::npos) return false;

    const double v = parse_double(s.substr(0, comma), NAN);
    const double e = parse_double(s.substr(comma + 1), NAN);
    if (!std::isfinite(v) || !std::isfinite(e)) return false;

    std::ostringstream out;
    out << "(" << std::setprecision(12) << factor * v
        << ", " << std::setprecision(12) << std::fabs(factor) * e << ")";
    cell = out.str();
    return true;
}

bool has_col(const std::map<std::string, std::size_t>& H,
             const std::string& name) {
    return H.find(name) != H.end();
}

double getd(const std::vector<std::string>& row,
            const std::map<std::string, std::size_t>& H,
            const std::string& name,
            double fallback = 0.0) {
    auto it = H.find(name);
    if (it == H.end() || it->second >= row.size()) return fallback;
    return parse_double(row[it->second], fallback);
}

// -----------------------------------------------------------------------------
// Physics helpers.
// -----------------------------------------------------------------------------

double photon_energy_from_bin(double xB, double Q2, double t_abs) {
    if (!(xB > 0.0) || !(Q2 > 0.0) || !(t_abs >= 0.0))
        return NAN;

    const double nu = Q2 / (2.0 * MP_GEV * xB);
    const double Eg = nu - t_abs / (2.0 * MP_GEV);
    return (std::isfinite(Eg) && Eg > 0.0) ? Eg : NAN;
}

double ratio_data_over_mc(bool photon_is_ft, double Egamma) {
    const bool high = Egamma >= ENERGY_SPLIT_GEV;

    if (photon_is_ft)
        return high ? EFF.ft_high : EFF.ft_low;

    return high ? EFF.fd_high : EFF.fd_low;
}

double correction_mc_over_data(bool photon_is_ft, double Egamma) {
    const double r = ratio_data_over_mc(photon_is_ft, Egamma);
    return (r > 0.0) ? 1.0 / r : 1.0;
}

struct EffectiveCorrection {
    bool valid = false;
    double correction = 1.0;

    double fd_weight = 0.0;
    double ft_weight = 0.0;
    double total_weight = 0.0;

    double fd_corrected_weight = 0.0;
    double ft_corrected_weight = 0.0;

    double mean_Egamma = NAN;
};

// Estimate the effective correction using the topology-resolved,
// current-corrected selected epgamma yields already stored in the CSV.
EffectiveCorrection effective_correction(
    const std::vector<std::string>& row,
    const std::map<std::string, std::size_t>& H,
    const std::string& label)
{
    EffectiveCorrection out;

    auto mit = LABEL_MEMBERS.find(label);
    if (mit == LABEL_MEMBERS.end()) return out;

    double sum_Ew = 0.0;
    double sum_w  = 0.0;
    double corrected = 0.0;

    for (const std::string& period : mit->second) {
        // Prefer period-specific mean kinematics.
        double xB = getd(row, H, "xBavg, " + period, NAN);
        double Q2 = getd(row, H, "Q2avg, " + period, NAN);
        double tt = getd(row, H, "t_abs_avg, " + period, NAN);

        // Fall back to target-label means if a period did not contribute enough
        // events to have its own finite mean.
        if (!std::isfinite(xB)) xB = getd(row, H, "xBavg, " + label, NAN);
        if (!std::isfinite(Q2)) Q2 = getd(row, H, "Q2avg, " + label, NAN);
        if (!std::isfinite(tt)) tt = getd(row, H, "t_abs_avg, " + label, NAN);

        const double Eg = photon_energy_from_bin(xB, Q2, tt);
        if (!std::isfinite(Eg)) continue;

        const std::string prefix =
            "normalized raw yield, ep->epg, ";

        const double y_fd_fd = std::max(
            0.0,
            getd(row, H,
                 prefix + "(FD, FD), exp, " + period + ", unpol",
                 0.0));

        const double y_cd_fd = std::max(
            0.0,
            getd(row, H,
                 prefix + "(CD, FD), exp, " + period + ", unpol",
                 0.0));

        const double y_cd_ft = std::max(
            0.0,
            getd(row, H,
                 prefix + "(CD, FT), exp, " + period + ", unpol",
                 0.0));

        const double y_fd = y_fd_fd + y_cd_fd;
        const double y_ft = y_cd_ft;

        const double c_fd = correction_mc_over_data(false, Eg);
        const double c_ft = correction_mc_over_data(true,  Eg);

        out.fd_weight += y_fd;
        out.ft_weight += y_ft;
        out.fd_corrected_weight += y_fd * c_fd;
        out.ft_corrected_weight += y_ft * c_ft;

        corrected += y_fd * c_fd + y_ft * c_ft;

        const double yw = y_fd + y_ft;
        sum_Ew += yw * Eg;
        sum_w  += yw;
    }

    out.total_weight = out.fd_weight + out.ft_weight;

    if (out.total_weight > 0.0) {
        out.valid = true;
        out.correction = corrected / out.total_weight;
        out.mean_Egamma = (sum_w > 0.0) ? sum_Ew / sum_w : NAN;
    }

    return out;
}

// Add a new output column if absent and return its index.
std::size_t ensure_column(
    std::vector<std::string>& header,
    std::map<std::string, std::size_t>& H,
    const std::string& name)
{
    auto it = H.find(name);
    if (it != H.end()) return it->second;

    const std::size_t idx = header.size();
    header.push_back(name);
    H[name] = idx;
    return idx;
}

void resize_row(std::vector<std::string>& row, std::size_t n) {
    if (row.size() < n) row.resize(n);
}

// Scale all absolute uncertainty columns that the world-data comparison
// interprets in the same units as the cross section. Fractional systematics are
// deliberately left unchanged.
void scale_absolute_uncertainty_if_present(
    std::vector<std::string>& row,
    const std::map<std::string, std::size_t>& H,
    const std::string& col,
    double factor)
{
    auto it = H.find(col);
    if (it == H.end() || it->second >= row.size()) return;

    const double x = parse_double(row[it->second], NAN);
    if (std::isfinite(x))
        row[it->second] = fmt(std::fabs(factor) * x);
}

} // namespace

int main(int argc, char** argv) {
    std::string input  = DEFAULT_INPUT;
    std::string output = DEFAULT_OUTPUT;

    if (argc >= 2) input = argv[1];
    if (argc >= 3) output = argv[2];

    std::ifstream in(input);
    if (!in) {
        std::cerr << "[photon-eff-prelim] ERROR: cannot open input CSV: "
                  << input << "\n";
        return 1;
    }

    std::string line;
    if (!std::getline(in, line)) {
        std::cerr << "[photon-eff-prelim] ERROR: empty input CSV\n";
        return 1;
    }

    std::vector<std::string> header = parse_csv_line(line);
    std::map<std::string, std::size_t> H = make_header_map(header);

    // Validate topology-yield inputs before touching anything.
    for (const std::string& period :
         std::vector<std::string>{"Fa18 Inb","Fa18 Out","Sp18 Inb","Sp18 Out"})
    {
        for (const std::string& topo :
             std::vector<std::string>{"(FD, FD)","(CD, FD)","(CD, FT)"})
        {
            const std::string c =
                "normalized raw yield, ep->epg, " + topo +
                ", exp, " + period + ", unpol";

            if (!has_col(H, c)) {
                std::cerr
                    << "[photon-eff-prelim] ERROR: required column missing:\n  "
                    << c << "\n";
                return 2;
            }
        }
    }

    // Audit columns.
    struct AuditCols {
        std::size_t correction;
        std::size_t data_mc_ratio_effective;
        std::size_t fd_fraction;
        std::size_t ft_fraction;
        std::size_t mean_Egamma;
    };

    std::map<std::string, AuditCols> audit;

    for (const std::string& label : DEFAULT_LABELS) {
        AuditCols a;
        a.correction = ensure_column(
            header, H,
            "photon eff prelim correction mc/data, " + label);

        a.data_mc_ratio_effective = ensure_column(
            header, H,
            "photon eff prelim effective data/mc, " + label);

        a.fd_fraction = ensure_column(
            header, H,
            "photon eff prelim FD photon yield fraction, " + label);

        a.ft_fraction = ensure_column(
            header, H,
            "photon eff prelim FT photon yield fraction, " + label);

        a.mean_Egamma = ensure_column(
            header, H,
            "photon eff prelim mean Egamma, " + label);

        audit[label] = a;

        // Preserve uncorrected XS cells for one-to-one auditing.
        for (const std::string& prefix :
             std::vector<std::string>{"cross sections, ", "normed cross sections, "})
        {
            const std::string original =
                prefix + "ep->epg, exp, " + label + ", unpol";

            if (has_col(H, original)) {
                ensure_column(
                    header, H,
                    "before photon eff prelim, " + original);
            }
        }
    }

    std::ofstream out(output);
    if (!out) {
        std::cerr << "[photon-eff-prelim] ERROR: cannot open output CSV: "
                  << output << "\n";
        return 3;
    }

    write_csv_row(out, header);

    std::size_t nrow = 0;
    std::size_t ncorr = 0;
    double min_corr = 1e99;
    double max_corr = -1e99;
    double sum_corr = 0.0;

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        std::vector<std::string> row = parse_csv_line(line);
        resize_row(row, header.size());
        ++nrow;

        for (const std::string& label : DEFAULT_LABELS) {
            const EffectiveCorrection ec =
                effective_correction(row, H, label);

            const AuditCols& a = audit.at(label);

            if (!ec.valid) {
                row[a.correction] = "1";
                row[a.data_mc_ratio_effective] = "1";
                row[a.fd_fraction] = "0";
                row[a.ft_fraction] = "0";
                row[a.mean_Egamma] = "";
                continue;
            }

            row[a.correction] = fmt(ec.correction);
            row[a.data_mc_ratio_effective] =
                fmt(ec.correction > 0.0 ? 1.0 / ec.correction : 1.0);

            row[a.fd_fraction] =
                fmt(ec.fd_weight / ec.total_weight);

            row[a.ft_fraction] =
                fmt(ec.ft_weight / ec.total_weight);

            row[a.mean_Egamma] =
                std::isfinite(ec.mean_Egamma) ? fmt(ec.mean_Egamma) : "";

            // Scale raw + normed final cross-section cells for this target label.
            for (const std::string& prefix :
                 std::vector<std::string>{"cross sections, ", "normed cross sections, "})
            {
                const std::string xscol =
                    prefix + "ep->epg, exp, " + label + ", unpol";

                auto xit = H.find(xscol);
                if (xit == H.end()) continue;

                const std::string backup =
                    "before photon eff prelim, " + xscol;

                auto bit = H.find(backup);
                if (bit != H.end())
                    row[bit->second] = row[xit->second];

                scale_tuple_cell(row[xit->second], ec.correction);
            }

            // For the finalized 10.6-GeV world comparison, scale absolute
            // point-to-point errors along with the cross section. Fractional
            // correlated/overall normalization systematics remain fractions.
            if (label == "10.6 GeV") {
                scale_absolute_uncertainty_if_present(
                    row, H, "Syst. err (point-to-point total)", ec.correction);
            }

            ++ncorr;
            min_corr = std::min(min_corr, ec.correction);
            max_corr = std::max(max_corr, ec.correction);
            sum_corr += ec.correction;
        }

        write_csv_row(out, row);
    }

    std::cout << "\n[photon-eff-prelim] wrote: " << output << "\n";
    std::cout << "[photon-eff-prelim] input rows: " << nrow << "\n";

    if (ncorr > 0) {
        std::cout
            << "[photon-eff-prelim] effective correction range = "
            << min_corr << " .. " << max_corr
            << ", mean = " << sum_corr / double(ncorr) << "\n";
    }

    std::cout << "\nCurrent exploratory epsilon_data/epsilon_MC constants:\n";
    std::cout << "  FD E<2  : " << EFF.fd_low
              << "  -> sigma multiplier " << 1.0/EFF.fd_low << "\n";
    std::cout << "  FD E>=2 : " << EFF.fd_high
              << "  -> sigma multiplier " << 1.0/EFF.fd_high << "\n";
    std::cout << "  FT E<2  : " << EFF.ft_low
              << "  -> sigma multiplier " << 1.0/EFF.ft_low << "\n";
    std::cout << "  FT E>=2 : " << EFF.ft_high
              << "  -> sigma multiplier " << 1.0/EFF.ft_high << "\n";

    std::cout
        << "\nWARNING: this is a POST-HOC sensitivity test only. "
        << "Do not use this CSV as the production correction.\n";

    return 0;
}
