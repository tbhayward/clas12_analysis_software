#include "branch_data_mc_comparison.h"
#include "branch_plot_config.h"
#include "exclusivity_cuts.h"
#include "global_cuts.h"

#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TROOT.h>
#include <TH1.h>
#include <TString.h>
#include <TDataType.h>
#include <TObjArray.h>
#include <TMath.h>

#include <nlohmann/json.hpp>

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
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

static constexpr bool kSkipGlobalCuts = false;
static constexpr bool kSkipExclusivityCuts = false;
static constexpr int kNClasdSectors = 6;
static constexpr double kPi = 3.14159265358979323846;
static constexpr double kRad2Deg = 180.0 / kPi;

struct PeriodDef {
    std::string pretty;
    std::string period_label;
    std::string data_key;
    std::string mc_key;
};

struct RangeInfo {
    double xmin = 0.0;
    double xmax = 1.0;
    int nbins = 100;
    bool is_integer_like = false;
    bool is_bool_like = false;
    bool from_config = false;
};

struct RowBin {
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;
    double phimin = 0.0;
    double phimax = 0.0;
    bool valid = false;
};

struct PeriodCorrections {
    double current_exp = 1.0;
    double poly[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
};

struct CsvInfo {
    std::vector<std::string> header;
    std::unordered_map<std::string, int> col;
    std::vector<std::vector<std::string> > rows;
    std::vector<RowBin> bins;
    std::map<std::string, PeriodCorrections> dvcs_corr;
    std::map<std::string, PeriodCorrections> eppi0_corr;
};

struct ChannelHistStore {
    std::vector<std::string> branch_names;
    std::map<std::string, RangeInfo> range_by_branch;
    std::map<std::string, std::vector<TH1D*> > data_hists_by_branch;
    std::map<std::string, std::vector<TH1D*> > mc_hists_by_branch;
    std::map<std::string, std::vector<std::vector<TH1D*> > > data_sector_hists_by_branch;
    std::map<std::string, std::vector<std::vector<TH1D*> > > mc_sector_hists_by_branch;
};

static std::mutex g_root_bind_mutex;

// -------------------- string / CSV helpers --------------------

static std::string sanitizeName(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_') out.push_back(c);
        else out.push_back('_');
    }
    return out;
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
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
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    out.push_back(cur);
    return out;
}

static double parse_double_or_nan(const std::string& s) {
    if (s.empty()) return std::numeric_limits<double>::quiet_NaN();
    char* endp = nullptr;
    const double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return std::numeric_limits<double>::quiet_NaN();
    return v;
}

static bool parse_tuple_numbers(const std::string& cell, std::vector<double>& out) {
    out.clear();
    std::string s = cell;
    for (char& c : s) {
        if (c == '(' || c == ')' || c == '[' || c == ']') c = ' ';
    }
    std::vector<std::string> parts = split_csv_line(s);
    for (const auto& p : parts) {
        const double v = parse_double_or_nan(p);
        if (std::isfinite(v)) out.push_back(v);
    }
    return !out.empty();
}

static int col_strict(const CsvInfo& csv, const std::string& name) {
    auto it = csv.col.find(name);
    if (it == csv.col.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing required CSV column: " << name;
        throw std::runtime_error(ss.str());
    }
    return it->second;
}

static std::string current_eff_col(const std::string& channel,
                                   const std::string& period) {
    return "current efficiency factor, " + channel + ", exp, " + period;
}

static std::string poly_col(const std::string& period) {
    return "eppi0 cross-section normalization polynomial, ep->eppi0, data_over_mc, " + period;
}

static std::string contam_col(const std::string& period) {
    return "contamination ratio, " + period;
}

static bool is_true_valid(const std::string& s) {
    return (s == "1" || s == "1.0" || s == "true" || s == "TRUE");
}

static double first_tuple_value_or_default(const std::string& cell, double def) {
    std::vector<double> vals;
    if (!parse_tuple_numbers(cell, vals)) return def;
    if (vals.empty() || !std::isfinite(vals[0])) return def;
    return vals[0];
}

static PeriodCorrections read_period_corrections(const CsvInfo& csv,
                                                 const std::string& channel,
                                                 const std::string& period) {
    PeriodCorrections c;

    const int c_eff = col_strict(csv, current_eff_col(channel, period));
    const int c_poly = col_strict(csv, poly_col(period));

    bool have_eff = false;
    bool have_poly = false;

    for (const auto& row : csv.rows) {
        if (!have_eff) {
            const double v = first_tuple_value_or_default(row[c_eff], std::numeric_limits<double>::quiet_NaN());
            if (std::isfinite(v) && v > 0.0) {
                c.current_exp = v;
                have_eff = true;
            }
        }

        if (!have_poly) {
            std::vector<double> coeffs;
            if (parse_tuple_numbers(row[c_poly], coeffs) && coeffs.size() >= 5) {
                bool ok = true;
                for (int i = 0; i < 5; ++i) ok = ok && std::isfinite(coeffs[(size_t)i]);
                if (ok) {
                    for (int i = 0; i < 5; ++i) c.poly[i] = coeffs[(size_t)i];
                    have_poly = true;
                }
            }
        }

        if (have_eff && have_poly) break;
    }

    if (!have_eff) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: no valid current-efficiency factor found for "
           << channel << " period " << period;
        throw std::runtime_error(ss.str());
    }

    if (!have_poly) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: no valid eppi0 normalization polynomial found for period "
           << period;
        throw std::runtime_error(ss.str());
    }

    return c;
}

static CsvInfo load_csv_info(const std::string& csv_path) {
    CsvInfo csv;

    std::ifstream fin(csv_path);
    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: cannot open CSV: " << csv_path;
        throw std::runtime_error(ss.str());
    }

    std::string line;
    if (!std::getline(fin, line)) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: empty CSV: " << csv_path;
        throw std::runtime_error(ss.str());
    }

    csv.header = split_csv_line(line);
    for (int i = 0; i < (int)csv.header.size(); ++i) {
        if (csv.col.find(csv.header[(size_t)i]) != csv.col.end()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: duplicate CSV column: " << csv.header[(size_t)i];
            throw std::runtime_error(ss.str());
        }
        csv.col[csv.header[(size_t)i]] = i;
    }

    while (std::getline(fin, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < csv.header.size()) row.resize(csv.header.size(), "");
        if (row.size() != csv.header.size()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: CSV row width mismatch in " << csv_path;
            throw std::runtime_error(ss.str());
        }
        csv.rows.push_back(std::move(row));
    }

    const int c_xmin = col_strict(csv, "xBmin");
    const int c_xmax = col_strict(csv, "xBmax");
    const int c_Qmin = col_strict(csv, "Q2min");
    const int c_Qmax = col_strict(csv, "Q2max");
    const int c_tmin = col_strict(csv, "t_abs_min");
    const int c_tmax = col_strict(csv, "t_abs_max");
    const int c_pmin = col_strict(csv, "phimin");
    const int c_pmax = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    csv.bins.reserve(csv.rows.size());
    for (const auto& row : csv.rows) {
        RowBin b;
        b.xBmin = parse_double_or_nan(row[c_xmin]);
        b.xBmax = parse_double_or_nan(row[c_xmax]);
        b.Q2min = parse_double_or_nan(row[c_Qmin]);
        b.Q2max = parse_double_or_nan(row[c_Qmax]);
        b.tmin = parse_double_or_nan(row[c_tmin]);
        b.tmax = parse_double_or_nan(row[c_tmax]);
        b.phimin = parse_double_or_nan(row[c_pmin]);
        b.phimax = parse_double_or_nan(row[c_pmax]);
        b.valid = is_true_valid(row[c_valid]);
        csv.bins.push_back(b);
    }

    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };

    for (const auto& p : periods) {
        csv.dvcs_corr[p] = read_period_corrections(csv, "ep->epg", p);
        csv.eppi0_corr[p] = read_period_corrections(csv, "ep->eppi0", p);
        (void)col_strict(csv, contam_col(p));
    }

    std::cout << "[branch_data_mc_comparison] Loaded CSV corrections and "
              << csv.bins.size() << " valid-bin rows from " << csv_path << std::endl;

    return csv;
}

// -------------------- labels and mapping --------------------

static std::string channelToStr(Channel ch) {
    return (ch == Channel::DVCS) ? "dvcs" : "eppi0";
}

static std::string channelCsvName(Channel ch) {
    return (ch == Channel::DVCS) ? "ep->epg" : "ep->eppi0";
}

static std::string channelPretty(Channel ch) {
    return (ch == Channel::DVCS) ? "DVCS" : "eppi0";
}

static std::string topoToKey(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

static std::string periodCode(Channel ch, const std::string& runTag) {
    const std::string prefix = (ch == Channel::DVCS) ? "DVCS_" : "eppi0_";

    if (runTag == "fa18_inb") return prefix + "Fa18_Inb";
    if (runTag == "fa18_out") return prefix + "Fa18_Out";
    if (runTag == "sp18_inb") return prefix + "Sp18_Inb";
    if (runTag == "sp18_out") return prefix + "Sp18_Out";
    if (runTag == "sp19_inb") return prefix + "Sp19_Inb";

    std::ostringstream ss;
    ss << "[branch_data_mc_comparison] FATAL: unknown runTag \"" << runTag
       << "\" in periodCode().";
    throw std::runtime_error(ss.str());
}

static bool topologyFromDetectors(int detector1, int detector2, Topology& topo_out) {
    if (detector1 == 1 && detector2 == 1) {
        topo_out = Topology::FD_FD;
        return true;
    }
    if (detector1 == 2 && detector2 == 1) {
        topo_out = Topology::CD_FD;
        return true;
    }
    if (detector1 == 2 && detector2 == 0) {
        topo_out = Topology::CD_FT;
        return true;
    }
    return false;
}

static bool within3Sigma(double val, const Stats& s) {
    if (!(std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0)) return true;
    return (val >= s.mean - 3.0 * s.std) && (val <= s.mean + 3.0 * s.std);
}

static bool passes3SigmaCuts(const std::map<std::string, Stats>& cuts,
                             const std::map<std::string, double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: missing value for 3-sigma cut variable "
               << kv.first;
            throw std::runtime_error(ss.str());
        }
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

static std::vector<PeriodDef> getDvcsPeriods() {
    return {
        {"Sp18 Inb", "sp18_inb", "DVCS_Sp18_inb", "DVCS_Sp18_inb_rec"},
        {"Sp18 Out", "sp18_out", "DVCS_Sp18_out", "DVCS_Sp18_out_rec"},
        {"Fa18 Inb", "fa18_inb", "DVCS_Fa18_inb", "DVCS_Fa18_inb_rec"},
        {"Fa18 Out", "fa18_out", "DVCS_Fa18_out", "DVCS_Fa18_out_rec"},
        {"Sp19 Inb", "sp19_inb", "DVCS_Sp19_inb", "DVCS_Sp19_inb_rec"}
    };
}

static std::vector<PeriodDef> getEppi0Periods() {
    return {
        {"Sp18 Inb", "sp18_inb", "DVCS_Sp18_inb_eppi0", "DVCS_Sp18_inb_rec_mc"},
        {"Sp18 Out", "sp18_out", "DVCS_Sp18_out_eppi0", "DVCS_Sp18_out_rec_mc"},
        {"Fa18 Inb", "fa18_inb", "DVCS_Fa18_inb_eppi0", "DVCS_Fa18_inb_rec_mc"},
        {"Fa18 Out", "fa18_out", "DVCS_Fa18_out_eppi0", "DVCS_Fa18_out_rec_mc"},
        {"Sp19 Inb", "sp19_inb", "DVCS_Sp19_inb_eppi0", "DVCS_Sp19_inb_rec_mc"}
    };
}

static bool shouldConvertRadiansToDegrees(const std::string& branch_name) {
    return (
        branch_name == "e_theta" ||
        branch_name == "p1_theta" ||
        branch_name == "p2_theta"
    );
}

static double convertValueForPlot(const std::string& branch_name, double value) {
    if (shouldConvertRadiansToDegrees(branch_name)) {
        return value * kRad2Deg;
    }
    return value;
}

static std::string axisTitleForBranch(const std::string& branch_name) {
    if (shouldConvertRadiansToDegrees(branch_name)) {
        return branch_name + " (deg)";
    }
    return branch_name;
}

static double wrapPhiDeg(double phi_deg) {
    double out = std::fmod(phi_deg, 360.0);
    if (out < 0.0) out += 360.0;
    if (out >= 360.0) out = std::nextafter(360.0, 0.0);
    return out;
}

static bool in_range(double v, double a, double b) {
    return (v >= a && v < b);
}

static bool row_accepts_phi(double phi, double a, double b) {
    if (b > a) return in_range(phi, a, b);
    return (phi >= a || phi < b);
}

static int find_csv_row_for_event(const CsvInfo& csv,
                                  double x,
                                  double Q2,
                                  double t_abs,
                                  double phi_deg) {
    for (int i = 0; i < (int)csv.bins.size(); ++i) {
        const RowBin& r = csv.bins[(size_t)i];
        if (!r.valid) continue;
        if (!in_range(x, r.xBmin, r.xBmax)) continue;
        if (!in_range(Q2, r.Q2min, r.Q2max)) continue;
        if (!in_range(t_abs, r.tmin, r.tmax)) continue;
        if (!row_accepts_phi(phi_deg, r.phimin, r.phimax)) continue;
        return i;
    }
    return -1;
}

static int getClas12SectorFromPhiDeg(double phi_deg) {
    const double phi = wrapPhiDeg(phi_deg);

    if ((phi >= 330.0 && phi < 360.0) || (phi >= 0.0 && phi < 30.0)) return 1;
    if (phi >= 30.0  && phi < 90.0)  return 2;
    if (phi >= 90.0  && phi < 150.0) return 3;
    if (phi >= 150.0 && phi < 210.0) return 4;
    if (phi >= 210.0 && phi < 270.0) return 5;
    if (phi >= 270.0 && phi < 330.0) return 6;

    return 0;
}

static std::string sectorSourcePhiBranch(const std::string& branch_name) {
    if (branch_name.rfind("p1_", 0) == 0) return "p1_phi";
    if (branch_name.rfind("p2_", 0) == 0) return "p2_phi";
    return "e_phi";
}

static std::string cutModeText() {
    if (kSkipGlobalCuts && kSkipExclusivityCuts) return "No global cuts, no exclusivity cuts";
    if (kSkipGlobalCuts && !kSkipExclusivityCuts) return "No global cuts, topology-matched 3#sigma exclusivity cuts";
    if (!kSkipGlobalCuts && kSkipExclusivityCuts) return "Global cuts only";
    return "Global cuts + topology-matched 3#sigma exclusivity cuts";
}

static double eval_poly4(const double coeffs[5], double x) {
    return coeffs[0] + x * (coeffs[1] + x * (coeffs[2] + x * (coeffs[3] + x * coeffs[4])));
}

static double event_weight(Channel ch,
                           const PeriodDef& period,
                           const CsvInfo& csv,
                           int csv_row,
                           double p1_theta_rad) {
    const std::map<std::string, PeriodCorrections>& corr_map =
        (ch == Channel::DVCS) ? csv.dvcs_corr : csv.eppi0_corr;

    auto it_corr = corr_map.find(period.pretty);
    if (it_corr == corr_map.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing corrections for period " << period.pretty;
        throw std::runtime_error(ss.str());
    }

    const PeriodCorrections& corr = it_corr->second;
    const double theta_deg = p1_theta_rad * kRad2Deg;
    double R = eval_poly4(corr.poly, theta_deg);

    if (!(std::isfinite(R) && R > 0.0)) {
        return 0.0;
    }

    double w = 1.0 / (corr.current_exp * R);

    if (ch == Channel::DVCS) {
        if (csv_row < 0 || csv_row >= (int)csv.rows.size()) return 0.0;
        const int c_contam = col_strict(csv, contam_col(period.pretty));
        const double contam = first_tuple_value_or_default(csv.rows[(size_t)csv_row][c_contam], 0.0);
        if (std::isfinite(contam)) {
            w *= (1.0 - contam);
        }
    }

    if (!std::isfinite(w) || w < 0.0) return 0.0;
    return w;
}

// -------------------- JSON cuts --------------------

static std::map<std::string, CutDict> loadCombinedCutsJson(const std::string& path) {
    std::cout << "[branch_data_mc_comparison] combined_cuts_json path = " << path << std::endl;

    std::ifstream fin(path);
    if (!fin.is_open()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: cannot open combined cuts JSON: " << path;
        throw std::runtime_error(ss.str());
    }

    nlohmann::json j;
    try {
        fin >> j;
    } catch (const std::exception& e) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: failed parsing " << path << ": " << e.what();
        throw std::runtime_error(ss.str());
    }

    std::map<std::string, CutDict> out;

    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string key = it.key();
        const auto& block = it.value();
        if (!block.is_object()) continue;

        CutDict cd;

        auto fill = [&](const char* sample, std::map<std::string, Stats>& target) {
            if (!block.contains(sample) || !block[sample].is_object()) return;
            for (auto vit = block[sample].begin(); vit != block[sample].end(); ++vit) {
                const auto& obj = vit.value();
                if (!obj.is_object()) continue;
                if (!obj.contains("mean") || !obj.contains("std")) continue;
                Stats s;
                s.mean = obj["mean"].get<double>();
                s.std = obj["std"].get<double>();
                target[vit.key()] = s;
            }
        };

        fill("data", cd.data);
        fill("mc", cd.mc);

        out[key] = cd;
    }

    return out;
}

// -------------------- branch inspection --------------------

static bool isSupportedScalarNumericLeaf(const TLeaf* leaf) {
    if (!leaf) return false;
    if (leaf->GetLenStatic() != 1) return false;
    if (leaf->GetLeafCount() != nullptr) return false;

    const std::string t = leaf->GetTypeName();

    return (
        t == "Bool_t"   ||
        t == "Char_t"   ||
        t == "UChar_t"  ||
        t == "Short_t"  ||
        t == "UShort_t" ||
        t == "Int_t"    ||
        t == "UInt_t"   ||
        t == "Long_t"   ||
        t == "ULong_t"  ||
        t == "Long64_t" ||
        t == "ULong64_t"||
        t == "Float_t"  ||
        t == "Double_t"
    );
}

static bool isIntegerLikeType(const TLeaf* leaf) {
    if (!leaf) return false;
    const std::string t = leaf->GetTypeName();
    return (
        t == "Bool_t"   ||
        t == "Char_t"   ||
        t == "UChar_t"  ||
        t == "Short_t"  ||
        t == "UShort_t" ||
        t == "Int_t"    ||
        t == "UInt_t"   ||
        t == "Long_t"   ||
        t == "ULong_t"  ||
        t == "Long64_t" ||
        t == "ULong64_t"
    );
}

static bool isBoolType(const TLeaf* leaf) {
    if (!leaf) return false;
    return std::string(leaf->GetTypeName()) == "Bool_t";
}

static TLeaf* getScalarNumericLeaf(TTree* tree, const std::string& branch_name) {
    if (!tree) return nullptr;

    TBranch* br = tree->GetBranch(branch_name.c_str());
    if (!br) return nullptr;

    TLeaf* leaf = br->GetLeaf(branch_name.c_str());
    if (!leaf) leaf = tree->GetLeaf(branch_name.c_str());
    if (!leaf) return nullptr;
    if (!isSupportedScalarNumericLeaf(leaf)) return nullptr;
    return leaf;
}

static std::set<std::string> getSupportedScalarNumericBranchNames(TTree* tree) {
    std::set<std::string> out;
    if (!tree) return out;

    TObjArray* branches = tree->GetListOfBranches();
    if (!branches) return out;

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch* br = dynamic_cast<TBranch*>(branches->At(i));
        if (!br) continue;

        const std::string name = br->GetName();
        TLeaf* leaf = br->GetLeaf(name.c_str());
        if (!leaf) leaf = tree->GetLeaf(name.c_str());
        if (!leaf) continue;

        if (isSupportedScalarNumericLeaf(leaf)) {
            out.insert(name);
        }
    }

    return out;
}

// -------------------- event binder for cuts and row matching --------------------

struct BranchBinder {
    int runnum = 0;       bool has_runnum = false;
    int detector1 = 0;    bool has_detector1 = false;
    int detector2 = 0;    bool has_detector2 = false;

    double x = 0.0;       bool has_x = false;
    double Q2 = 0.0;      bool has_Q2 = false;
    double phi2 = 0.0;    bool has_phi2 = false;

    double t1 = 0.0;                bool has_t1 = false;
    double open_angle_ep2 = 0.0;    bool has_open_angle_ep2 = false;
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_1 = 0.0;             bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;
    double pTmiss = 0.0;            bool has_pTmiss = false;
    double xF = 0.0;                bool has_xF = false;
    double Delta_phi = 0.0;         bool has_Delta_phi = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0;     bool has_theta_pi0_pi0 = false;

    double e_p = 0.0;       bool has_e_p = false;
    double e_theta = 0.0;   bool has_e_theta = false;
    double e_phi = 0.0;     bool has_e_phi = false;

    double p1_theta = 0.0;  bool has_p1_theta = false;
    double p1_phi = 0.0;    bool has_p1_phi = false;

    double p2_p = 0.0;      bool has_p2_p = false;
    double p2_theta = 0.0;  bool has_p2_theta = false;
    double p2_phi = 0.0;    bool has_p2_phi = false;

    void bindCutBranches(TTree* t, Channel ch) {
        if (!t) return;

        auto ena = [&](const char* n) {
            if (t->GetBranch(n)) t->SetBranchStatus(n, 1);
        };

        ena("runnum");
        ena("detector1");
        ena("detector2");
        ena("x");
        ena("Q2");
        ena("phi2");
        ena("t1");
        ena("open_angle_ep2");
        ena("Emiss2");
        ena("Mx2");
        ena("Mx2_1");
        ena("Mx2_2");
        ena("pTmiss");
        ena("xF");
        ena("Delta_phi");
        if (ch == Channel::DVCS) ena("theta_gamma_gamma");
        else                     ena("theta_pi0_pi0");

        ena("e_p");
        ena("e_theta");
        ena("e_phi");
        ena("p1_theta");
        ena("p1_phi");
        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        auto bI = [&](const char* n, int* a, bool& f) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, a, nullptr, nullptr, kInt_t, false);
                f = true;
            }
        };

        auto bD = [&](const char* n, double* a, bool& f) {
            if (t->GetBranch(n)) {
                t->SetBranchAddress(n, a, nullptr, nullptr, kDouble_t, false);
                f = true;
            }
        };

        bI("runnum",    &runnum,    has_runnum);
        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);

        bD("x", &x, has_x);
        bD("Q2", &Q2, has_Q2);
        bD("phi2", &phi2, has_phi2);
        bD("t1", &t1, has_t1);
        bD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bD("Emiss2", &Emiss2, has_Emiss2);
        bD("Mx2", &Mx2, has_Mx2);
        bD("Mx2_1", &Mx2_1, has_Mx2_1);
        bD("Mx2_2", &Mx2_2, has_Mx2_2);
        bD("pTmiss", &pTmiss, has_pTmiss);
        bD("xF", &xF, has_xF);
        bD("Delta_phi", &Delta_phi, has_Delta_phi);

        if (ch == Channel::DVCS) bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        else                     bD("theta_pi0_pi0",     &theta_pi0_pi0,     has_theta_pi0_pi0);

        bD("e_p", &e_p, has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi", &e_phi, has_e_phi);
        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi", &p1_phi, has_p1_phi);
        bD("p2_p", &p2_p, has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi", &p2_phi, has_p2_phi);
    }

    std::map<std::string, double> valuesMap(Channel ch) const {
        std::map<std::string, double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
        if (ch == Channel::DVCS) {
            if (has_theta_gamma_gamma) m["theta_gamma_gamma"] = theta_gamma_gamma;
        } else {
            if (has_theta_pi0_pi0) m["theta_pi0_pi0"] = theta_pi0_pi0;
        }
        if (has_pTmiss) m["pTmiss"] = pTmiss;
        if (has_xF) m["xF"] = xF;
        if (has_Emiss2) m["Emiss2"] = Emiss2;
        if (has_Mx2) m["Mx2"] = Mx2;
        if (has_Mx2_1) m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2) m["Mx2_2"] = Mx2_2;
        return m;
    }
};

// -------------------- tree and branch validation --------------------

static void requireAllPeriodTreesPresent(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees)
{
    std::vector<std::string> missing;

    for (const auto& p : periods) {
        auto itd = dataTrees.find(p.data_key);
        if (itd == dataTrees.end() || itd->second == nullptr) {
            missing.push_back("data:" + p.data_key);
        }

        auto itm = recMcTrees.find(p.mc_key);
        if (itm == recMcTrees.end() || itm->second == nullptr) {
            missing.push_back("rec_mc:" + p.mc_key);
        }
    }

    if (!missing.empty()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] Missing required trees for channel "
           << channel_name << ": ";
        for (size_t i = 0; i < missing.size(); ++i) {
            ss << missing[i];
            if (i + 1 < missing.size()) ss << ", ";
        }
        throw std::runtime_error(ss.str());
    }
}

static std::vector<std::string> getCommonBranchesAcrossPeriods(
    const std::string& channel_name,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees)
{
    requireAllPeriodTreesPresent(channel_name, periods, dataTrees, recMcTrees);

    bool first = true;
    std::set<std::string> common;

    for (const auto& p : periods) {
        TTree* dt = dataTrees.at(p.data_key);
        TTree* mt = recMcTrees.at(p.mc_key);

        std::set<std::string> sdata = getSupportedScalarNumericBranchNames(dt);
        std::set<std::string> smc   = getSupportedScalarNumericBranchNames(mt);

        std::set<std::string> both;
        std::set_intersection(
            sdata.begin(), sdata.end(),
            smc.begin(), smc.end(),
            std::inserter(both, both.begin())
        );

        if (first) {
            common = both;
            first = false;
        } else {
            std::set<std::string> next_common;
            std::set_intersection(
                common.begin(), common.end(),
                both.begin(), both.end(),
                std::inserter(next_common, next_common.begin())
            );
            common.swap(next_common);
        }
    }

    std::vector<std::string> out;
    for (const auto& b : common) {
        const BranchHistConfig* cfg = findBranchPlotConfig(b);
        if (cfg && !cfg->enabled) continue;
        out.push_back(b);
    }

    std::sort(out.begin(), out.end());
    return out;
}

// -------------------- histogram range logic --------------------

static RangeInfo determineAutoRangeForBranch(
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& branch_name)
{
    RangeInfo info;

    bool initialized = false;
    double global_min = 0.0;
    double global_max = 0.0;
    bool integer_like = true;
    bool bool_like = true;

    for (const auto& p : periods) {
        TTree* dt = dataTrees.at(p.data_key);
        TTree* mt = recMcTrees.at(p.mc_key);

        TLeaf* dleaf = getScalarNumericLeaf(dt, branch_name);
        TLeaf* mleaf = getScalarNumericLeaf(mt, branch_name);

        if (!dleaf || !mleaf) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] Branch " << branch_name
               << " is not available as scalar numeric in all required trees.";
            throw std::runtime_error(ss.str());
        }

        integer_like = integer_like && isIntegerLikeType(dleaf) && isIntegerLikeType(mleaf);
        bool_like = bool_like && isBoolType(dleaf) && isBoolType(mleaf);

        const double dmin_raw = dt->GetMinimum(branch_name.c_str());
        const double dmax_raw = dt->GetMaximum(branch_name.c_str());
        const double mmin_raw = mt->GetMinimum(branch_name.c_str());
        const double mmax_raw = mt->GetMaximum(branch_name.c_str());

        const double dmin = convertValueForPlot(branch_name, dmin_raw);
        const double dmax = convertValueForPlot(branch_name, dmax_raw);
        const double mmin = convertValueForPlot(branch_name, mmin_raw);
        const double mmax = convertValueForPlot(branch_name, mmax_raw);

        if (!initialized) {
            global_min = std::min(dmin, mmin);
            global_max = std::max(dmax, mmax);
            initialized = true;
        } else {
            global_min = std::min(global_min, std::min(dmin, mmin));
            global_max = std::max(global_max, std::max(dmax, mmax));
        }
    }

    info.is_integer_like = integer_like;
    info.is_bool_like = bool_like;
    info.from_config = false;

    if (!initialized) {
        info.xmin = 0.0;
        info.xmax = 1.0;
        info.nbins = 100;
        return info;
    }

    if (bool_like) {
        info.xmin = -0.5;
        info.xmax = 1.5;
        info.nbins = 2;
        return info;
    }

    if (integer_like) {
        const long long imin = static_cast<long long>(std::floor(global_min));
        const long long imax = static_cast<long long>(std::ceil(global_max));
        const long long span = imax - imin + 1LL;

        if (span >= 1LL && span <= 200LL) {
            info.xmin = static_cast<double>(imin) - 0.5;
            info.xmax = static_cast<double>(imax) + 0.5;
            info.nbins = static_cast<int>(span);
            return info;
        }
    }

    if (global_min == global_max) {
        double delta = 1.0;
        if (std::abs(global_min) > 0.0) delta = 0.05 * std::abs(global_min);
        info.xmin = global_min - delta;
        info.xmax = global_max + delta;
        info.nbins = 100;
        return info;
    }

    const double width = global_max - global_min;
    const double pad = 0.05 * width;

    info.xmin = global_min - pad;
    info.xmax = global_max + pad;
    info.nbins = 100;

    return info;
}

static RangeInfo determineRangeForBranch(
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::string& branch_name)
{
    const BranchHistConfig* cfg = findBranchPlotConfig(branch_name);
    if (cfg) {
        if (!cfg->enabled) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] Branch " << branch_name
               << " is explicitly disabled in branch_plot_config.";
            throw std::runtime_error(ss.str());
        }
        if (!(cfg->xhigh > cfg->xlow) || cfg->nbins <= 0) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] Invalid configured histogram range for branch "
               << branch_name << ": nbins=" << cfg->nbins
               << " xlow=" << cfg->xlow
               << " xhigh=" << cfg->xhigh;
            throw std::runtime_error(ss.str());
        }

        RangeInfo info;
        info.xmin = cfg->xlow;
        info.xmax = cfg->xhigh;
        info.nbins = cfg->nbins;
        info.from_config = true;
        return info;
    }

    return determineAutoRangeForBranch(periods, dataTrees, recMcTrees, branch_name);
}

// -------------------- one-pass filling --------------------

static bool passesGlobalCutsForBinder(const BranchBinder& b,
                                      const std::string& period_label)
{
    if (kSkipGlobalCuts) return true;

    if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) return false;
    if (!(b.has_detector1 && b.has_detector2)) return false;

    GlobalCutConfig gcfg = default_global_cuts();
    gcfg.enable_topology_filter = false;

    if (global_cuts_require_sector_phi(gcfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            throw std::runtime_error("[branch_data_mc_comparison] FATAL: sector selection requires e_phi, p1_phi, p2_phi.");
        }
    }

    if (gcfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            throw std::runtime_error("[branch_data_mc_comparison] FATAL: dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.");
        }

        if (global_cuts_require_sector_phi(gcfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      period_label,
                                      b.e_p, b.e_theta, b.e_phi, b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      gcfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  gcfg);
    }

    if (global_cuts_require_sector_phi(gcfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  gcfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              gcfg);
}

static void enablePlotBranches(TTree* tree, const std::vector<std::string>& branch_names) {
    for (const auto& name : branch_names) {
        if (tree->GetBranch(name.c_str())) {
            tree->SetBranchStatus(name.c_str(), 1);
        }
    }
}

static std::map<std::string, TLeaf*> buildPlotLeafMap(
    TTree* tree,
    const std::vector<std::string>& branch_names)
{
    std::map<std::string, TLeaf*> leaf_map;
    for (const auto& name : branch_names) {
        TLeaf* leaf = getScalarNumericLeaf(tree, name);
        if (!leaf) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: expected scalar numeric leaf for branch "
               << name << " but did not find one while building plot leaf map.";
            throw std::runtime_error(ss.str());
        }
        leaf_map[name] = leaf;
    }
    return leaf_map;
}

static int determineSectorForBranch(const std::string& branch_name, const BranchBinder& b) {
    const std::string sector_phi_branch = sectorSourcePhiBranch(branch_name);

    double phi_rad = 0.0;
    bool have_phi = false;

    if (sector_phi_branch == "e_phi") {
        phi_rad = b.e_phi;
        have_phi = b.has_e_phi;
    } else if (sector_phi_branch == "p1_phi") {
        phi_rad = b.p1_phi;
        have_phi = b.has_p1_phi;
    } else if (sector_phi_branch == "p2_phi") {
        phi_rad = b.p2_phi;
        have_phi = b.has_p2_phi;
    } else {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: unrecognized sector phi source branch "
           << sector_phi_branch << " for plot branch " << branch_name;
        throw std::runtime_error(ss.str());
    }

    if (!have_phi) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing required phi branch "
           << sector_phi_branch << " needed to sectorize plot branch " << branch_name;
        throw std::runtime_error(ss.str());
    }

    const double phi_deg = phi_rad * kRad2Deg;
    const int sector = getClas12SectorFromPhiDeg(phi_deg);

    if (sector < 1 || sector > kNClasdSectors) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: computed invalid sector " << sector
           << " from " << sector_phi_branch << " for plot branch " << branch_name;
        throw std::runtime_error(ss.str());
    }

    return sector;
}

static void fillHistogramsForTreeSinglePass(
    TTree* tree,
    Channel ch,
    const PeriodDef& period,
    const std::vector<std::string>& branch_names,
    const std::map<std::string, TLeaf*>& leaf_map,
    const std::map<std::string, CutDict>& combinedCuts,
    const CsvInfo& csv,
    bool use_data_cuts,
    bool is_data,
    std::map<std::string, TH1D*>& hist_by_branch,
    std::map<std::string, std::vector<TH1D*> >& sector_hists_by_branch)
{
    if (!tree) {
        throw std::runtime_error("[branch_data_mc_comparison] Null TTree passed to fillHistogramsForTreeSinglePass.");
    }

    BranchBinder b;
    b.bindCutBranches(tree, ch);

    GlobalCutConfig gcfg = default_global_cuts();

    const Long64_t nentries = tree->GetEntries();
    long long n_pass = 0;
    long long n_row = 0;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        if (!kSkipGlobalCuts) {
            if (b.has_runnum && is_excluded_run(b.runnum, gcfg)) continue;
        }

        if (!(b.has_detector1 && b.has_detector2)) continue;

        Topology topo;
        if (!topologyFromDetectors(b.detector1, b.detector2, topo)) continue;

        if (!passesGlobalCutsForBinder(b, period.period_label)) continue;

        if (!kSkipExclusivityCuts) {
            const std::string cut_key = periodCode(ch, period.period_label) + "_" + topoToKey(topo);
            auto itCuts = combinedCuts.find(cut_key);
            if (itCuts == combinedCuts.end()) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: missing combined cuts entry for key "
                   << cut_key;
                throw std::runtime_error(ss.str());
            }

            const std::map<std::string, double> vals = b.valuesMap(ch);
            const std::map<std::string, Stats>& cut_map = use_data_cuts ? itCuts->second.data : itCuts->second.mc;

            if (!passes3SigmaCuts(cut_map, vals)) continue;
        }

        if (!(b.has_x && b.has_Q2 && b.has_t1 && b.has_phi2)) {
            throw std::runtime_error("[branch_data_mc_comparison] FATAL: missing x/Q2/t1/phi2 needed for CSV-bin correction lookup.");
        }

        const int csv_row = find_csv_row_for_event(csv,
                                                   b.x,
                                                   b.Q2,
                                                   std::fabs(b.t1),
                                                   wrapPhiDeg(b.phi2 * kRad2Deg));
        if (csv_row < 0) continue;
        ++n_row;

        double weight = 1.0;
        if (is_data) {
            if (!b.has_p1_theta) {
                throw std::runtime_error("[branch_data_mc_comparison] FATAL: missing p1_theta needed for eppi0 normalization polynomial weight.");
            }
            weight = event_weight(ch, period, csv, csv_row, b.p1_theta);
            if (!(weight > 0.0)) continue;
        }

        ++n_pass;

        for (const auto& branch_name : branch_names) {
            auto itLeaf = leaf_map.find(branch_name);
            if (itLeaf == leaf_map.end() || itLeaf->second == nullptr) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: missing plot leaf during fill for branch "
                   << branch_name;
                throw std::runtime_error(ss.str());
            }

            auto itHist = hist_by_branch.find(branch_name);
            if (itHist == hist_by_branch.end() || itHist->second == nullptr) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: missing histogram during fill for branch "
                   << branch_name;
                throw std::runtime_error(ss.str());
            }

            auto itSectorVec = sector_hists_by_branch.find(branch_name);
            if (itSectorVec == sector_hists_by_branch.end() || static_cast<int>(itSectorVec->second.size()) != kNClasdSectors) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: missing sector histogram vector during fill for branch "
                   << branch_name;
                throw std::runtime_error(ss.str());
            }

            const double raw_value = itLeaf->second->GetValue(0);
            const double plot_value = convertValueForPlot(branch_name, raw_value);
            itHist->second->Fill(plot_value, weight);

            const int sector = determineSectorForBranch(branch_name, b);
            TH1D* h_sector = itSectorVec->second[sector - 1];
            if (!h_sector) {
                std::ostringstream ss;
                ss << "[branch_data_mc_comparison] FATAL: null sector histogram for branch "
                   << branch_name << " sector " << sector;
                throw std::runtime_error(ss.str());
            }
            h_sector->Fill(plot_value, weight);
        }
    }

    std::cout << "[branch_data_mc_comparison] Filled "
              << (is_data ? "data" : "rec MC") << " channel=" << channelToStr(ch)
              << " period=" << period.pretty
              << " entries=" << (long long)nentries
              << " csv_bin_matched=" << n_row
              << " filled=" << n_pass
              << std::endl;
}

static void normalizeHist(TH1D* h) {
    if (!h) return;
    const double integral = h->Integral();
    if (integral > 0.0) h->Scale(1.0 / integral);
}

static void styleHistogram(TH1D* h, int color, int marker_style) {
    if (!h) return;
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(marker_style);
    h->SetMarkerSize(0.75);
    h->SetLineWidth(2);
    h->SetStats(0);
}

static std::string makeCanvasTitle(const std::string& channel_name, const std::string& branch_name) {
    if (channel_name == "eppi0") return "e p #rightarrow e p #pi_{0} : " + branch_name;
    return "DVCS : " + branch_name;
}

// -------------------- histogram creation --------------------

static ChannelHistStore createHistStore(
    Channel ch,
    const std::vector<PeriodDef>& periods,
    const std::vector<std::string>& branch_names,
    const std::map<std::string, RangeInfo>& range_by_branch)
{
    ChannelHistStore store;
    store.branch_names = branch_names;
    store.range_by_branch = range_by_branch;

    for (const auto& branch_name : branch_names) {
        const auto itRange = range_by_branch.find(branch_name);
        if (itRange == range_by_branch.end()) {
            std::ostringstream ss;
            ss << "[branch_data_mc_comparison] FATAL: missing range info for branch "
               << branch_name;
            throw std::runtime_error(ss.str());
        }

        const RangeInfo& rinfo = itRange->second;

        std::vector<TH1D*> data_vec;
        std::vector<TH1D*> mc_vec;
        std::vector<std::vector<TH1D*> > data_sector_vec;
        std::vector<std::vector<TH1D*> > mc_sector_vec;

        data_vec.reserve(periods.size());
        mc_vec.reserve(periods.size());
        data_sector_vec.reserve(periods.size());
        mc_sector_vec.reserve(periods.size());

        for (size_t ip = 0; ip < periods.size(); ++ip) {
            const PeriodDef& p = periods[ip];

            TH1D* hd = new TH1D(
                ("h_data_" + sanitizeName(channelToStr(ch)) + "_" +
                 sanitizeName(p.period_label) + "_" + sanitizeName(branch_name)).c_str(),
                "",
                rinfo.nbins,
                rinfo.xmin,
                rinfo.xmax
            );
            hd->SetDirectory(nullptr);
            hd->Sumw2();

            TH1D* hm = new TH1D(
                ("h_mc_" + sanitizeName(channelToStr(ch)) + "_" +
                 sanitizeName(p.period_label) + "_" + sanitizeName(branch_name)).c_str(),
                "",
                rinfo.nbins,
                rinfo.xmin,
                rinfo.xmax
            );
            hm->SetDirectory(nullptr);
            hm->Sumw2();

            std::vector<TH1D*> data_sector_period_vec;
            std::vector<TH1D*> mc_sector_period_vec;
            data_sector_period_vec.reserve(kNClasdSectors);
            mc_sector_period_vec.reserve(kNClasdSectors);

            for (int isec = 1; isec <= kNClasdSectors; ++isec) {
                TH1D* hd_sec = new TH1D(
                    ("h_data_" + sanitizeName(channelToStr(ch)) + "_" +
                     sanitizeName(p.period_label) + "_" + sanitizeName(branch_name) +
                     "_sector_" + std::to_string(isec)).c_str(),
                    "",
                    rinfo.nbins,
                    rinfo.xmin,
                    rinfo.xmax
                );
                hd_sec->SetDirectory(nullptr);
                hd_sec->Sumw2();

                TH1D* hm_sec = new TH1D(
                    ("h_mc_" + sanitizeName(channelToStr(ch)) + "_" +
                     sanitizeName(p.period_label) + "_" + sanitizeName(branch_name) +
                     "_sector_" + std::to_string(isec)).c_str(),
                    "",
                    rinfo.nbins,
                    rinfo.xmin,
                    rinfo.xmax
                );
                hm_sec->SetDirectory(nullptr);
                hm_sec->Sumw2();

                data_sector_period_vec.push_back(hd_sec);
                mc_sector_period_vec.push_back(hm_sec);
            }

            data_vec.push_back(hd);
            mc_vec.push_back(hm);
            data_sector_vec.push_back(data_sector_period_vec);
            mc_sector_vec.push_back(mc_sector_period_vec);
        }

        store.data_hists_by_branch[branch_name] = data_vec;
        store.mc_hists_by_branch[branch_name] = mc_vec;
        store.data_sector_hists_by_branch[branch_name] = data_sector_vec;
        store.mc_sector_hists_by_branch[branch_name] = mc_sector_vec;
    }

    return store;
}

static void deleteHistStore(ChannelHistStore& store) {
    for (auto& kv : store.data_hists_by_branch) {
        for (TH1D* h : kv.second) delete h;
    }
    for (auto& kv : store.mc_hists_by_branch) {
        for (TH1D* h : kv.second) delete h;
    }
    for (auto& kv : store.data_sector_hists_by_branch) {
        for (auto& per_period : kv.second) {
            for (TH1D* h : per_period) delete h;
        }
    }
    for (auto& kv : store.mc_sector_hists_by_branch) {
        for (auto& per_period : kv.second) {
            for (TH1D* h : per_period) delete h;
        }
    }

    store.data_hists_by_branch.clear();
    store.mc_hists_by_branch.clear();
    store.data_sector_hists_by_branch.clear();
    store.mc_sector_hists_by_branch.clear();
}

// -------------------- plotting --------------------

static void saveOneBranchCanvas(
    Channel ch,
    const std::vector<PeriodDef>& periods,
    const ChannelHistStore& store,
    const std::string& branch_name,
    const std::string& channel_out_dir)
{
    auto itRange = store.range_by_branch.find(branch_name);
    if (itRange == store.range_by_branch.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing range for branch "
           << branch_name << " during saveOneBranchCanvas";
        throw std::runtime_error(ss.str());
    }
    const RangeInfo& rinfo = itRange->second;

    auto itData = store.data_hists_by_branch.find(branch_name);
    auto itMc   = store.mc_hists_by_branch.find(branch_name);
    if (itData == store.data_hists_by_branch.end() || itMc == store.mc_hists_by_branch.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing hist vectors for branch "
           << branch_name;
        throw std::runtime_error(ss.str());
    }

    std::vector<TH1D*> data_hists = itData->second;
    std::vector<TH1D*> mc_hists   = itMc->second;

    double global_ymax = 0.0;

    for (size_t i = 0; i < periods.size(); ++i) {
        normalizeHist(data_hists[i]);
        normalizeHist(mc_hists[i]);

        styleHistogram(data_hists[i], kBlue, 20);
        styleHistogram(mc_hists[i],   kRed,  24);

        global_ymax = std::max(global_ymax, data_hists[i]->GetMaximum());
        global_ymax = std::max(global_ymax, mc_hists[i]->GetMaximum());
    }

    if (global_ymax <= 0.0) global_ymax = 1.0;

    TCanvas c(
        ("c_" + sanitizeName(channelToStr(ch)) + "_" + sanitizeName(branch_name)).c_str(),
        "",
        2100,
        1200
    );
    c.Divide(3, 2, 0.002, 0.002);

    for (size_t i = 0; i < periods.size(); ++i) {
        c.cd(static_cast<int>(i) + 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.12);
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        TH1D* hd = data_hists[i];
        TH1D* hm = mc_hists[i];

        hd->SetTitle("");
        hd->GetXaxis()->SetTitle(axisTitleForBranch(branch_name).c_str());
        hd->GetYaxis()->SetTitle("Normalized counts");
        hd->GetXaxis()->SetTitleOffset(1.10);
        hd->GetYaxis()->SetTitleOffset(1.70);
        hd->SetMaximum(1.25 * global_ymax);

        hd->Draw("HIST");
        hm->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.54, 0.66, 0.92, 0.84);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(hd, "Data, corrected", "l");
        leg->AddEntry(hm, "Reconstructed MC", "l");
        leg->Draw();

        TLatex period_latex;
        period_latex.SetNDC();
        period_latex.SetTextFont(42);
        period_latex.SetTextSize(0.055);
        period_latex.SetTextAlign(13);
        period_latex.DrawLatex(0.18, 0.92, periods[i].pretty.c_str());
    }

    c.cd(6);
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.05);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAlign(22);
    latex.SetTextFont(42);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.50, 0.75, makeCanvasTitle(channelToStr(ch), branch_name).c_str());
    latex.SetTextSize(0.038);
    latex.DrawLatex(0.50, 0.60, cutModeText().c_str());
    latex.DrawLatex(0.50, 0.51, "Data weighted by current-efficiency and e#pi_{0} normalization corrections");
    if (ch == Channel::DVCS) latex.DrawLatex(0.50, 0.43, "DVCS data also weighted by row-level (1 - #pi_{0} contamination)");
    else                     latex.DrawLatex(0.50, 0.43, "e#pi_{0} data corrected with the e#pi_{0} normalization polynomial");
    latex.DrawLatex(0.50, 0.35, "Each histogram normalized to its own integral");

    std::ostringstream ss;
    ss << "Range: [" << rinfo.xmin << ", " << rinfo.xmax << "], bins = " << rinfo.nbins;
    latex.DrawLatex(0.50, 0.27, ss.str().c_str());

    if (rinfo.from_config) latex.DrawLatex(0.50, 0.19, "Using explicit branch_plot_config range");
    else                   latex.DrawLatex(0.50, 0.19, "Using automatic branch range");

    const std::string out_name =
        channel_out_dir + "/" + sanitizeName(branch_name) + "_data_vs_rec_mc.png";

    c.SaveAs(out_name.c_str());
}

static void saveOneBranchSectorCanvas(
    Channel ch,
    const PeriodDef& period,
    size_t period_index,
    const ChannelHistStore& store,
    const std::string& branch_name,
    const std::string& sector_out_dir)
{
    auto itRange = store.range_by_branch.find(branch_name);
    if (itRange == store.range_by_branch.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing range for branch "
           << branch_name << " during saveOneBranchSectorCanvas";
        throw std::runtime_error(ss.str());
    }
    const RangeInfo& rinfo = itRange->second;

    auto itData = store.data_sector_hists_by_branch.find(branch_name);
    auto itMc   = store.mc_sector_hists_by_branch.find(branch_name);
    if (itData == store.data_sector_hists_by_branch.end() || itMc == store.mc_sector_hists_by_branch.end()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: missing sector hist vectors for branch "
           << branch_name;
        throw std::runtime_error(ss.str());
    }

    if (period_index >= itData->second.size() || period_index >= itMc->second.size()) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: invalid period_index " << period_index
           << " for branch " << branch_name << " during saveOneBranchSectorCanvas";
        throw std::runtime_error(ss.str());
    }

    std::vector<TH1D*> data_hists = itData->second[period_index];
    std::vector<TH1D*> mc_hists   = itMc->second[period_index];

    if (static_cast<int>(data_hists.size()) != kNClasdSectors ||
        static_cast<int>(mc_hists.size())   != kNClasdSectors) {
        std::ostringstream ss;
        ss << "[branch_data_mc_comparison] FATAL: sector histogram size mismatch for branch "
           << branch_name << " period " << period.period_label;
        throw std::runtime_error(ss.str());
    }

    double global_ymax = 0.0;

    for (int isec = 0; isec < kNClasdSectors; ++isec) {
        normalizeHist(data_hists[isec]);
        normalizeHist(mc_hists[isec]);

        styleHistogram(data_hists[isec], kBlue, 20);
        styleHistogram(mc_hists[isec],   kRed,  24);

        global_ymax = std::max(global_ymax, data_hists[isec]->GetMaximum());
        global_ymax = std::max(global_ymax, mc_hists[isec]->GetMaximum());
    }

    if (global_ymax <= 0.0) global_ymax = 1.0;

    TCanvas c(
        ("c_sector_" + sanitizeName(channelToStr(ch)) + "_" +
         sanitizeName(period.period_label) + "_" + sanitizeName(branch_name)).c_str(),
        "",
        2100,
        1200
    );
    c.Divide(3, 2, 0.002, 0.002);

    for (int isec = 0; isec < kNClasdSectors; ++isec) {
        c.cd(isec + 1);
        gPad->SetLeftMargin(0.16);
        gPad->SetRightMargin(0.06);
        gPad->SetBottomMargin(0.13);
        gPad->SetTopMargin(0.12);
        gPad->SetTickx(1);
        gPad->SetTicky(1);

        TH1D* hd = data_hists[isec];
        TH1D* hm = mc_hists[isec];

        hd->SetTitle("");
        hd->GetXaxis()->SetTitle(axisTitleForBranch(branch_name).c_str());
        hd->GetYaxis()->SetTitle("Normalized counts");
        hd->GetXaxis()->SetTitleOffset(1.10);
        hd->GetYaxis()->SetTitleOffset(1.70);
        hd->SetMaximum(1.25 * global_ymax);

        hd->Draw("HIST");
        hm->Draw("HIST SAME");

        TLegend* leg = new TLegend(0.54, 0.66, 0.92, 0.84);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetTextFont(42);
        leg->SetTextSize(0.035);
        leg->AddEntry(hd, "Data, corrected", "l");
        leg->AddEntry(hm, "Reconstructed MC", "l");
        leg->Draw();

        TLatex pad_latex;
        pad_latex.SetNDC();
        pad_latex.SetTextFont(42);
        pad_latex.SetTextSize(0.050);
        pad_latex.SetTextAlign(13);
        pad_latex.DrawLatex(0.18, 0.92, ("Sector " + std::to_string(isec + 1)).c_str());
    }

    c.cd(6);
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.05);
    gPad->SetTopMargin(0.05);

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAlign(22);
    latex.SetTextFont(42);
    latex.SetTextSize(0.050);
    latex.DrawLatex(0.50, 0.74, (makeCanvasTitle(channelToStr(ch), branch_name) + "   " + period.pretty).c_str());
    latex.SetTextSize(0.038);
    latex.DrawLatex(0.50, 0.60, "Sector-by-sector");
    latex.DrawLatex(0.50, 0.50, cutModeText().c_str());
    latex.DrawLatex(0.50, 0.42, "Data are corrected with CSV current-efficiency and e#pi_{0} normalization weights");

    std::ostringstream ss;
    ss << "Range: [" << rinfo.xmin << ", " << rinfo.xmax << "], bins = " << rinfo.nbins;
    latex.DrawLatex(0.50, 0.32, ss.str().c_str());

    const std::string out_name =
        sector_out_dir + "/" + sanitizeName(branch_name) + "_" +
        sanitizeName(period.period_label) + "_sector_by_sector_data_vs_rec_mc.png";

    c.SaveAs(out_name.c_str());
}

// -------------------- channel driver --------------------

static void fill_one_period_pair(
    Channel ch,
    const PeriodDef& p,
    size_t ip,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::vector<std::string>& branch_names,
    const std::map<std::string, CutDict>& combinedCuts,
    const CsvInfo& csv,
    ChannelHistStore& store)
{
    {
        TTree* tree = dataTrees.at(p.data_key);
        {
            std::lock_guard<std::mutex> lock(g_root_bind_mutex);
            tree->SetBranchStatus("*", 0);
            enablePlotBranches(tree, branch_names);
        }

        std::map<std::string, TLeaf*> leaf_map = buildPlotLeafMap(tree, branch_names);

        std::map<std::string, TH1D*> hist_by_branch;
        std::map<std::string, std::vector<TH1D*> > sector_hist_by_branch;

        for (const auto& branch_name : branch_names) {
            hist_by_branch[branch_name] = store.data_hists_by_branch.at(branch_name)[ip];
            sector_hist_by_branch[branch_name] = store.data_sector_hists_by_branch.at(branch_name)[ip];
        }

        fillHistogramsForTreeSinglePass(
            tree,
            ch,
            p,
            branch_names,
            leaf_map,
            combinedCuts,
            csv,
            true,
            true,
            hist_by_branch,
            sector_hist_by_branch
        );

        std::lock_guard<std::mutex> lock(g_root_bind_mutex);
        tree->SetBranchStatus("*", 1);
    }

    {
        TTree* tree = recMcTrees.at(p.mc_key);
        {
            std::lock_guard<std::mutex> lock(g_root_bind_mutex);
            tree->SetBranchStatus("*", 0);
            enablePlotBranches(tree, branch_names);
        }

        std::map<std::string, TLeaf*> leaf_map = buildPlotLeafMap(tree, branch_names);

        std::map<std::string, TH1D*> hist_by_branch;
        std::map<std::string, std::vector<TH1D*> > sector_hist_by_branch;

        for (const auto& branch_name : branch_names) {
            hist_by_branch[branch_name] = store.mc_hists_by_branch.at(branch_name)[ip];
            sector_hist_by_branch[branch_name] = store.mc_sector_hists_by_branch.at(branch_name)[ip];
        }

        fillHistogramsForTreeSinglePass(
            tree,
            ch,
            p,
            branch_names,
            leaf_map,
            combinedCuts,
            csv,
            false,
            false,
            hist_by_branch,
            sector_hist_by_branch
        );

        std::lock_guard<std::mutex> lock(g_root_bind_mutex);
        tree->SetBranchStatus("*", 1);
    }
}

static void runChannelComparisons(
    Channel ch,
    const std::vector<PeriodDef>& periods,
    const std::map<std::string, TTree*>& dataTrees,
    const std::map<std::string, TTree*>& recMcTrees,
    const std::map<std::string, CutDict>& combinedCuts,
    const CsvInfo& csv,
    const std::string& outPlotDir,
    int max_workers)
{
    const std::string channel_name = channelToStr(ch);

    requireAllPeriodTreesPresent(channel_name, periods, dataTrees, recMcTrees);

    const std::vector<std::string> branch_names =
        getCommonBranchesAcrossPeriods(channel_name, periods, dataTrees, recMcTrees);

    if (branch_names.empty()) {
        std::cout << "[branch_data_mc_comparison] No common scalar numeric branches found for "
                  << channel_name << ". Skipping." << std::endl;
        return;
    }

    std::map<std::string, RangeInfo> range_by_branch;
    for (const auto& branch_name : branch_names) {
        range_by_branch[branch_name] = determineRangeForBranch(periods, dataTrees, recMcTrees, branch_name);
    }

    ChannelHistStore store = createHistStore(ch, periods, branch_names, range_by_branch);

    const std::string channel_out_dir = outPlotDir + "/" + channel_name;
    const std::string sector_out_dir = channel_out_dir + "/sector_by_sector";
    std::filesystem::create_directories(channel_out_dir);
    std::filesystem::create_directories(sector_out_dir);

    std::cout << "[branch_data_mc_comparison] " << channel_name
              << ": found " << branch_names.size()
              << " common scalar numeric branches." << std::endl;

    int nth = std::max(1, std::min(5, max_workers));
    nth = std::min(nth, std::max(1, (int)periods.size()));

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int ip = 0; ip < (int)periods.size(); ++ip) {
        fill_one_period_pair(ch,
                             periods[(size_t)ip],
                             (size_t)ip,
                             dataTrees,
                             recMcTrees,
                             branch_names,
                             combinedCuts,
                             csv,
                             store);
    }

    for (size_t i = 0; i < branch_names.size(); ++i) {
        const std::string& branch_name = branch_names[i];
        std::cout << "[branch_data_mc_comparison] " << channel_name
                  << " : (" << (i + 1) << "/" << branch_names.size() << ") "
                  << branch_name << std::endl;

        try {
            saveOneBranchCanvas(
                ch,
                periods,
                store,
                branch_name,
                channel_out_dir
            );
        } catch (const std::exception& e) {
            std::cerr << "[branch_data_mc_comparison] WARNING: skipping overall branch plot "
                      << branch_name << " because: " << e.what() << std::endl;
        }

        for (size_t ip = 0; ip < periods.size(); ++ip) {
            try {
                saveOneBranchSectorCanvas(
                    ch,
                    periods[ip],
                    ip,
                    store,
                    branch_name,
                    sector_out_dir
                );
            } catch (const std::exception& e) {
                std::cerr << "[branch_data_mc_comparison] WARNING: skipping sector plot for branch "
                          << branch_name << ", period " << periods[ip].period_label
                          << " because: " << e.what() << std::endl;
            }
        }
    }

    deleteHistStore(store);
}

} // namespace

bool runAllBranchDataMcComparisons(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& outPlotDir,
    int max_workers)
{
    try {
        ROOT::EnableThreadSafety();
        TH1::AddDirectory(kFALSE);
        gROOT->SetBatch(kTRUE);
        gStyle->SetOptStat(0);

        const CsvInfo csv = load_csv_info(csv_path);

        std::map<std::string, CutDict> combinedCuts;
        if (!kSkipExclusivityCuts) {
            combinedCuts = loadCombinedCutsJson(combined_cuts_json);
        }

        const std::string base_out_dir = outPlotDir + "/branch_data_mc_comparisons";
        std::filesystem::create_directories(base_out_dir);

        try {
            runChannelComparisons(
                Channel::DVCS,
                getDvcsPeriods(),
                dvcsDataTrees,
                dvcsRecMcTrees,
                combinedCuts,
                csv,
                base_out_dir,
                max_workers
            );
        } catch (const std::exception& e) {
            std::cerr << "[branch_data_mc_comparison] DVCS skipped: " << e.what() << std::endl;
        }

        try {
            runChannelComparisons(
                Channel::EPPI0,
                getEppi0Periods(),
                eppi0DataTrees,
                eppi0RecMcTrees,
                combinedCuts,
                csv,
                base_out_dir,
                max_workers
            );
        } catch (const std::exception& e) {
            std::cerr << "[branch_data_mc_comparison] eppi0 skipped: " << e.what() << std::endl;
        }

        std::cout << "[branch_data_mc_comparison] Done." << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}
