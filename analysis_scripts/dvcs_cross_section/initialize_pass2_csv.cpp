// initialize_pass2_csv.cpp

#include "initialize_pass2_csv.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/**
 * OVERVIEW
 * --------
 * 1) Read header row from Lee's CSV and build a "name -> index" lookup.
 * 2) Build the new pass-2 analysis header.
 * 3) Validate that the new header contains no duplicate column names.
 * 4) Iterate Lee's rows, keeping only rows where "valid bin" == 1.
 * 5) For each kept row, populate a new row vector sized to the new header.
 *    - Copy requested columns from Lee.
 *    - Leave all analysis-produced columns blank.
 * 6) Write the new CSV file.
 *
 * IMPORTANT SCHEMA NOTES
 * ----------------------
 * - The current-efficiency correction is represented upstream of the final
 *   cross-section calculation. cross_sections.cpp should not read or apply
 *   imports/efficiency.json after this refactor.
 *
 * - Current-efficiency factors are stored as one data factor and one MC factor
 *   per channel and period. Each cell is intended to hold a two-tuple:
 *
 *       (value,stat)
 *
 * - The eppi0 cross-section normalization factor is represented separately
 *   from the current-efficiency correction. The fitted eppi0 ratio is:
 *
 *       R_pi0(theta_p) = N_data(theta_p) / N_MC(theta_p)
 *
 *   The DVCS and eppi0 data yield corrections later divide event weights by
 *   this factor.
 *
 * - Generated and reconstructed MC yields are written to the CSV so that
 *   acceptance and normalization steps can be computed algebraically from the
 *   CSV rather than by re-looping over the ROOT trees.
 *
 * - Reconstructed-current-corrected MC columns are created here but should be
 *   filled later, after current-efficiency factors have been determined.
 */

/* ============================
 * Minimal CSV utilities
 * ============================ */

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
            for (char ch : s) {
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

/* ============================
 * Header lookup helpers
 * ============================ */

static std::unordered_map<std::string, int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> m;

    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }

    return m;
}

static std::string get_col(const std::vector<std::string>& row,
                           const std::unordered_map<std::string, int>& idx,
                           const std::string& name) {
    const auto it = idx.find(name);

    if (it == idx.end()) {
        return std::string();
    }

    const int j = it->second;

    if (j < 0 || j >= (int)row.size()) {
        return std::string();
    }

    return row[j];
}

static int ToInteger(const std::string& s) {
    if (s.empty()) {
        return 0;
    }

    char* endp = nullptr;
    const long v = std::strtol(s.c_str(), &endp, 10);

    if (endp == s.c_str()) {
        return 0;
    }

    return (int)v;
}

static bool header_has_duplicates(const std::vector<std::string>& header,
                                  std::vector<std::string>& duplicates) {
    std::set<std::string> seen;
    std::set<std::string> dup_set;

    for (const auto& name : header) {
        if (!seen.insert(name).second) {
            dup_set.insert(name);
        }
    }

    duplicates.assign(dup_set.begin(), dup_set.end());
    return !duplicates.empty();
}

/* ============================
 * Shared schema lists
 * ============================ */

static const std::vector<std::string>& base_periods() {
    static const std::vector<std::string> v = {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb",
        "Sp18 Out"
    };
    return v;
}

static const std::vector<std::string>& avg_groups() {
    static const std::vector<std::string> v = {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb",
        "Sp18 Out",
        "Fa18",
        "Sp18",
        "10.6 GeV"
    };
    return v;
}

static const std::vector<std::string>& cross_section_groups() {
    static const std::vector<std::string> v = {
        "10.6 GeV",
        "Fa18",
        "Sp18"
    };
    return v;
}

static const std::vector<std::string>& acceptance_corrected_groups() {
    static const std::vector<std::string> v = {
        "Fa18",
        "Sp18",
        "2018 (10.6 GeV)"
    };
    return v;
}

static const std::vector<std::string>& helicities() {
    static const std::vector<std::string> v = {
        "unpol",
        "pos",
        "neg"
    };
    return v;
}

static const std::vector<std::string>& topologies() {
    static const std::vector<std::string> v = {
        "(FD, FD)",
        "(CD, FD)",
        "(CD, FT)"
    };
    return v;
}

static const std::vector<std::string>& physics_channels() {
    static const std::vector<std::string> v = {
        "ep->epg",
        "ep->eppi0"
    };
    return v;
}

/* ============================
 * Build the NEW header schema
 * ============================ */

static void add_grouped_avg_columns(std::vector<std::string>& H,
                                    const std::string& base) {
    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << base << ", " << g;
        H.push_back(name.str());
    }
}

static void add_theta_group_columns_after_phiavg_10p6(std::vector<std::string>& H) {
    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << "e_theta, " << g;
        H.push_back(name.str());
    }

    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << "p_theta, " << g;
        H.push_back(name.str());
    }

    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << "g_theta, " << g;
        H.push_back(name.str());
    }
}

static void add_bin_definition_columns(std::vector<std::string>& H) {
    H.push_back("bin index");
    H.push_back("Bin Name");

    H.push_back("xBmin");
    H.push_back("xBmax");
    add_grouped_avg_columns(H, "xBavg");

    H.push_back("Q2min");
    H.push_back("Q2max");
    add_grouped_avg_columns(H, "Q2avg");

    H.push_back("t_abs_min");
    H.push_back("t_abs_max");
    add_grouped_avg_columns(H, "t_abs_avg");

    H.push_back("phimin");
    H.push_back("phimax");
    add_grouped_avg_columns(H, "phiavg");

    add_theta_group_columns_after_phiavg_10p6(H);

    H.push_back("tmin");
    H.push_back("tcol");
    H.push_back("P1");
    H.push_back("P2");
    H.push_back("maxP1P2");
    H.push_back("minP1P2");
    H.push_back("intP1P2");
}

static void add_raw_yield_columns_for_channel(std::vector<std::string>& H,
                                              const std::string& channel) {
    for (const auto& per : base_periods()) {
        for (const auto& topo : topologies()) {
            for (const auto& hel : helicities()) {
                std::ostringstream name;
                name << "raw yield, " << channel << ", " << topo
                     << ", exp, " << per << ", " << hel;
                H.push_back(name.str());
            }
        }
    }
}

static void add_raw_yield_columns(std::vector<std::string>& H) {
    add_raw_yield_columns_for_channel(H, "ep->epg");
    add_raw_yield_columns_for_channel(H, "ep->eppi0");
}

static void add_current_efficiency_columns(std::vector<std::string>& H) {
    for (const auto& channel : physics_channels()) {
        for (const auto& per : base_periods()) {
            std::ostringstream name;
            name << "current efficiency factor, " << channel << ", exp, " << per;
            H.push_back(name.str());
        }

        for (const auto& per : base_periods()) {
            std::ostringstream name;
            name << "current efficiency factor, " << channel << ", mc, " << per;
            H.push_back(name.str());
        }
    }
}

static std::vector<std::string> eppi0_normalization_regions() {
    return {"Sector 1", "Sector 2", "Sector 3", "Sector 4", "Sector 5", "Sector 6", "CD"};
}

static void add_eppi0_normalization_columns(std::vector<std::string>& H) {
    for (const auto& per : base_periods()) {
        for (const auto& region : eppi0_normalization_regions()) {
            std::ostringstream name;
            name << "eppi0 cross-section normalization factor, ep->eppi0, data_over_mc, "
                 << region << ", " << per;
            H.push_back(name.str());
        }
    }

    for (const auto& per : base_periods()) {
        for (const auto& region : eppi0_normalization_regions()) {
            std::ostringstream name;
            name << "eppi0 cross-section normalization cubic, ep->eppi0, data_over_mc, "
                 << region << ", " << per;
            H.push_back(name.str());
        }
    }
}

static void add_normalized_raw_yield_columns_for_channel(std::vector<std::string>& H,
                                                         const std::string& channel) {
    for (const auto& per : base_periods()) {
        for (const auto& topo : topologies()) {
            for (const auto& hel : helicities()) {
                std::ostringstream name;
                name << "normalized raw yield, " << channel << ", " << topo
                     << ", exp, " << per << ", " << hel;
                H.push_back(name.str());
            }
        }
    }
}

static void add_normalized_raw_yield_columns(std::vector<std::string>& H) {
    add_normalized_raw_yield_columns_for_channel(H, "ep->epg");
    add_normalized_raw_yield_columns_for_channel(H, "ep->eppi0");
}

static void add_mc_yield_columns_for_channel(std::vector<std::string>& H,
                                             const std::string& channel) {
    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "generated yield, " << channel << ", mc, " << per;
        H.push_back(name.str());
    }

    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "reconstructed yield, " << channel << ", mc, " << per;
        H.push_back(name.str());
    }

    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "reconstructed current corrected yield, " << channel << ", mc, " << per;
        H.push_back(name.str());
    }

    for (const auto& per : base_periods()) {
        for (const auto& topo : topologies()) {
            std::ostringstream name;
            name << "reconstructed yield, " << channel << ", " << topo
                 << ", mc, " << per;
            H.push_back(name.str());
        }
    }

    for (const auto& per : base_periods()) {
        for (const auto& topo : topologies()) {
            std::ostringstream name;
            name << "reconstructed current corrected yield, " << channel << ", "
                 << topo << ", mc, " << per;
            H.push_back(name.str());
        }
    }
}

static void add_eppi0_background_mc_yield_columns(std::vector<std::string>& H) {
    const std::string channel = "ep->eppi0->epg";

    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "reconstructed yield, " << channel << ", mc, " << per;
        H.push_back(name.str());
    }

    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "reconstructed current corrected yield, " << channel << ", mc, " << per;
        H.push_back(name.str());
    }

    for (const auto& per : base_periods()) {
        for (const auto& topo : topologies()) {
            std::ostringstream name;
            name << "reconstructed yield, " << channel << ", " << topo
                 << ", mc, " << per;
            H.push_back(name.str());
        }
    }

    for (const auto& per : base_periods()) {
        for (const auto& topo : topologies()) {
            std::ostringstream name;
            name << "reconstructed current corrected yield, " << channel << ", "
                 << topo << ", mc, " << per;
            H.push_back(name.str());
        }
    }
}

static void add_mc_yield_columns(std::vector<std::string>& H) {
    add_mc_yield_columns_for_channel(H, "ep->epg");
    add_mc_yield_columns_for_channel(H, "ep->eppi0");
    add_eppi0_background_mc_yield_columns(H);
}

static void add_contamination_columns(std::vector<std::string>& H) {
    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "contamination ratio, " << per;
        H.push_back(name.str());
    }
}

static void add_signal_yield_columns(std::vector<std::string>& H) {
    for (const auto& per : base_periods()) {
        for (const auto& hel : helicities()) {
            std::ostringstream name;
            name << "signal yield, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }
}

static void add_acceptance_columns(std::vector<std::string>& H) {
    for (const auto& per : base_periods()) {
        std::ostringstream name;
        name << "acceptance, " << per;
        H.push_back(name.str());
    }
}

static void add_acc_corrected_yield_columns(std::vector<std::string>& H) {
    for (const auto& per : base_periods()) {
        for (const auto& hel : helicities()) {
            std::ostringstream name;
            name << "acceptance corrected yield, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }

    for (const auto& grp : acceptance_corrected_groups()) {
        for (const auto& hel : helicities()) {
            std::ostringstream name;
            name << "acceptance corrected yield, ep->epg, exp, " << grp << ", " << hel;
            H.push_back(name.str());
        }
    }
}

static void add_frad_fbin_binvol_columns(std::vector<std::string>& H) {
    H.push_back("Frad, 10.6 GeV");
    H.push_back("Frad, 10.2 GeV");
    H.push_back("Fbin, 10.6 GeV");
    H.push_back("Fbin, 10.2 GeV");
    H.push_back("bin_volume, 10.6 GeV");
    H.push_back("bin_volume, 10.2 GeV");
}

static void add_luminosity_columns(std::vector<std::string>& H) {
    H.push_back("integrated luminosity, Fa18 Inb (nC)");
    H.push_back("integrated luminosity, Fa18 Out (nC)");
    H.push_back("integrated luminosity, Sp19 Inb (nC)");
    H.push_back("integrated luminosity, Sp18 Inb (nC)");
    H.push_back("integrated luminosity, Sp18 Out (nC)");
    H.push_back("integrated luminosity, Fa18 (nC)");
    H.push_back("integrated luminosity, Sp18 (nC)");
    H.push_back("integrated luminosity, 10.6 GeV (nC)");
}

static void add_norm_columns(std::vector<std::string>& H) {
    H.push_back("norm, Fa18 Inb");
    H.push_back("norm, Fa18 Out");
    H.push_back("norm, Sp19 Inb");
    H.push_back("norm, Sp18 Inb");
    H.push_back("norm, Sp18 Out");
    H.push_back("norm, Fa18");
    H.push_back("norm, Sp18");
    H.push_back("norm, 10.6 GeV");
}

static void add_cross_section_columns_without_normed(std::vector<std::string>& H,
                                                     const std::string& prefix) {
    for (const auto& per : base_periods()) {
        for (const auto& hel : helicities()) {
            std::ostringstream name;
            name << prefix << "cross sections, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }

    for (const auto& grp : cross_section_groups()) {
        for (const auto& hel : helicities()) {
            std::ostringstream name;
            name << prefix << "cross sections, ep->epg, exp, " << grp << ", " << hel;
            H.push_back(name.str());
        }
    }
}

static void add_cross_section_columns(std::vector<std::string>& H) {
    add_cross_section_columns_without_normed(H, "");
    add_norm_columns(H);
    add_cross_section_columns_without_normed(H, "normed ");
}

static void add_bsa_columns(std::vector<std::string>& H) {
    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << "BSA, counts, " << g;
        H.push_back(name.str());
    }

    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << "BSA, sigma, " << g;
        H.push_back(name.str());
    }
}

static void add_valid_and_prefactor_columns(std::vector<std::string>& H) {
    H.push_back("valid bin");

    H.push_back("cross section prefactor, 10.6 GeV");
    H.push_back("cross section prefactor, 10.2 GeV");

    H.push_back("BH2 prefactor, 10.6 GeV");
    H.push_back("BH2 prefactor, 10.2 GeV");

    H.push_back("BH2 c0, 10.6 GeV");
    H.push_back("BH2 c0, 10.2 GeV");

    H.push_back("BH2 c1, 10.6 GeV");
    H.push_back("BH2 c1, 10.2 GeV");

    H.push_back("BH2 c2, 10.6 GeV");
    H.push_back("BH2 c2, 10.2 GeV");

    H.push_back("Int prefactor, 10.6 GeV");
    H.push_back("Int prefactor, 10.2 GeV");

    H.push_back("Int c0 (KM15), 10.6 GeV");
    H.push_back("Int c0 (KM15), 10.2 GeV");

    H.push_back("Int c1 (KM15), 10.6 GeV");
    H.push_back("Int c1 (KM15), 10.2 GeV");

    H.push_back("Int c2 (KM15), 10.6 GeV");
    H.push_back("Int c2 (KM15), 10.2 GeV");

    H.push_back("DVCS2 prefactor, 10.6 GeV");
    H.push_back("DVCS2 prefactor, 10.2 GeV");

    H.push_back("DVCS2 c0 (KM15), 10.6 GeV");
    H.push_back("DVCS2 c0 (KM15), 10.2 GeV");

    H.push_back("DVCS2 c1 (KM15), 10.6 GeV");
    H.push_back("DVCS2 c1 (KM15), 10.2 GeV");
}

static std::vector<std::string> build_new_header() {
    std::vector<std::string> H;
    H.reserve(1050);

    add_bin_definition_columns(H);

    add_raw_yield_columns(H);

    add_current_efficiency_columns(H);
    add_eppi0_normalization_columns(H);
    add_normalized_raw_yield_columns(H);

    add_mc_yield_columns(H);

    add_contamination_columns(H);
    add_signal_yield_columns(H);
    add_acceptance_columns(H);
    add_acc_corrected_yield_columns(H);

    add_frad_fbin_binvol_columns(H);
    add_luminosity_columns(H);
    add_cross_section_columns(H);
    add_bsa_columns(H);
    add_valid_and_prefactor_columns(H);

    return H;
}

/* ============================
 * Copy plan: which Lee columns to copy
 * ============================ */

struct CopyMap {
    std::vector<std::pair<std::string, std::string>> pairs;
};

static CopyMap build_copy_map() {
    CopyMap M;

    M.pairs.push_back({"bin index", "bin index"});
    M.pairs.push_back({"Bin Name", "Bin Name"});

    M.pairs.push_back({"xBmin", "xBmin"});
    M.pairs.push_back({"xBmax", "xBmax"});

    M.pairs.push_back({"Q2min", "Q2min"});
    M.pairs.push_back({"Q2max", "Q2max"});

    M.pairs.push_back({"t_abs_min", "t_abs_min"});
    M.pairs.push_back({"t_abs_max", "t_abs_max"});

    M.pairs.push_back({"phimin", "phimin"});
    M.pairs.push_back({"phimax", "phimax"});

    M.pairs.push_back({"tmin", "tmin"});
    M.pairs.push_back({"tcol", "tcol"});
    M.pairs.push_back({"P1", "P1"});
    M.pairs.push_back({"P2", "P2"});
    M.pairs.push_back({"maxP1P2", "maxP1P2"});
    M.pairs.push_back({"minP1P2", "minP1P2"});
    M.pairs.push_back({"intP1P2", "intP1P2"});

    M.pairs.push_back({"valid bin", "valid bin"});

    M.pairs.push_back({"cross section prefactor, 10.6 GeV", "cross section prefactor"});

    M.pairs.push_back({"BH2 prefactor, 10.6 GeV", "BH2 prefactor"});
    M.pairs.push_back({"BH2 c0, 10.6 GeV", "BH2 c0"});
    M.pairs.push_back({"BH2 c1, 10.6 GeV", "BH2 c1"});
    M.pairs.push_back({"BH2 c2, 10.6 GeV", "BH2 c2"});

    M.pairs.push_back({"Int prefactor, 10.6 GeV", "Int prefactor"});
    M.pairs.push_back({"Int c0 (KM15), 10.6 GeV", "Int c0 (KM15)"});
    M.pairs.push_back({"Int c1 (KM15), 10.6 GeV", "Int c1 (KM15)"});
    M.pairs.push_back({"Int c2 (KM15), 10.6 GeV", "Int c2 (KM15)"});

    M.pairs.push_back({"DVCS2 prefactor, 10.6 GeV", "DVCS2 prefactor"});
    M.pairs.push_back({"DVCS2 c0 (KM15), 10.6 GeV", "DVCS2 c0 (KM15)"});
    M.pairs.push_back({"DVCS2 c1 (KM15), 10.6 GeV", "DVCS2 c1 (KM15)"});

    return M;
}

/* ============================
 * Main initializer
 * ============================ */

bool initialize_pass2_csv(const std::string& lee_csv_path,
                          const std::string& out_csv_path) {
    std::ifstream fin(lee_csv_path);

    if (!fin.is_open()) {
        std::cerr << "[init] ERROR: Could not open Lee CSV: "
                  << lee_csv_path << std::endl;
        return false;
    }

    std::string lee_header_line;

    if (!std::getline(fin, lee_header_line)) {
        std::cerr << "[init] ERROR: Lee CSV appears to be empty: "
                  << lee_csv_path << std::endl;
        return false;
    }

    std::vector<std::string> lee_header = split_csv_line(lee_header_line);

    if (!lee_header.empty()) {
        std::string first = lee_header[0];

        bool unnamed = first.empty();

        if (!unnamed) {
            std::string low;
            low.reserve(first.size());

            for (char c : first) {
                low.push_back((char)std::tolower((unsigned char)c));
            }

            if (low.find("unnamed") != std::string::npos) {
                unnamed = true;
            }
        }

        if (unnamed) {
            lee_header[0] = "bin index";
        }
    }

    const std::unordered_map<std::string, int> lee_index =
        build_header_index(lee_header);

    std::vector<std::string> new_header = build_new_header();

    std::vector<std::string> duplicates;

    if (header_has_duplicates(new_header, duplicates)) {
        std::cerr << "[init] FATAL: New CSV header contains duplicate column names:\n";

        for (const auto& d : duplicates) {
            std::cerr << "  - " << d << "\n";
        }

        return false;
    }

    const std::unordered_map<std::string, int> new_index =
        build_header_index(new_header);

    std::ofstream fout(out_csv_path);

    if (!fout.is_open()) {
        std::cerr << "[init] ERROR: Could not open output CSV for writing: "
                  << out_csv_path << std::endl;
        return false;
    }

    fout << join_csv_row(new_header) << "\n";

    const CopyMap cmap = build_copy_map();

    int input_rows = 0;
    int kept_rows = 0;

    std::string line;

    while (std::getline(fin, line)) {
        ++input_rows;

        if (line.empty()) {
            continue;
        }

        std::vector<std::string> cols = split_csv_line(line);

        const std::string valid_s = get_col(cols, lee_index, "valid bin");
        const int valid = ToInteger(valid_s);

        if (valid != 1) {
            continue;
        }

        std::vector<std::string> out_row(new_header.size(), std::string());

        for (const auto& kv : cmap.pairs) {
            const std::string& new_name = kv.first;
            const std::string& old_name = kv.second;

            const auto it_new = new_index.find(new_name);

            if (it_new == new_index.end()) {
                std::cerr << "[init] FATAL: Copy-map target column is missing from new header: "
                          << new_name << std::endl;
                return false;
            }

            const int new_col = it_new->second;

            if (new_col < 0 || new_col >= (int)out_row.size()) {
                std::cerr << "[init] FATAL: Internal new-column index out of range for column: "
                          << new_name << std::endl;
                return false;
            }

            out_row[(size_t)new_col] = get_col(cols, lee_index, old_name);
        }

        fout << join_csv_row(out_row) << "\n";
        ++kept_rows;
    }

    std::cout << "[init] Lee rows read: " << input_rows << "\n";
    std::cout << "[init] Valid rows kept (valid bin == 1): " << kept_rows << "\n";
    std::cout << "[init] Output columns: " << new_header.size() << "\n";
    std::cout << "[init] Wrote output: " << out_csv_path << "\n";

    return true;
}