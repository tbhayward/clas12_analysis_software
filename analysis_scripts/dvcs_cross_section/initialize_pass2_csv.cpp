#include "initialize_pass2_csv.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <iterator>

/**
 * OVERVIEW
 * --------
 * 1) Read header row from Lee's CSV and build a "name -> index" lookup.
 * 2) Build the new header (the extended schema you specified).
 * 3) Iterate Lee's rows, keep only those with "valid bin" == 1.
 * 4) For each kept row, populate a new row vector sized to the new header.
 *    - Copy requested columns from Lee.
 *    - Leave all other analysis-produced columns blank for now.
 * 5) Write the new CSV file.
 *
 * GROUPED AVERAGES
 * ----------------
 * We only keep period/combined-group averages (no singletons). For each variable
 * (xB, Q2, t_abs, phi) the eight grouped averages appear right after the matching
 * "...max" column.
 *
 * ADDITIONS (THIS UPDATE)
 * -----------------------
 * 1) After "phiavg, 10.6 GeV", add grouped theta columns:
 *      - e_theta, <group> for specific groups (as requested, note Sp18 Inb omitted)
 *      - p_theta, <group> for groups (including Sp18 Inb)
 *      - g_theta, <group> for groups (including Sp18 Inb)
 *
 * 2) After "cross sections, ep->epg, exp, Sp18, neg", add:
 *      - norm, <group> for {Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out, Fa18, Sp18, 10.6 GeV}
 *      - then add "normed " + <each existing cross section column> (same list/order as cross sections)
 */

/* ============================
 * Minimal CSV utilities
 * ============================ */

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;

    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
            continue;
        }
        if (c == ',' && !in_quotes) {
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
        bool need_quotes = (s.find(',') != std::string::npos) || (s.find('"') != std::string::npos);
        if (need_quotes) {
            oss << '"';
            for (char ch : s) {
                if (ch == '"') oss << "\"\""; else oss << ch;
            }
            oss << '"';
        } else {
            oss << s;
        }
        if (i + 1 < fields.size()) oss << ',';
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
                           const std::unordered_map<std::string,int>& idx,
                           const std::string& name) {
    auto it = idx.find(name);
    if (it == idx.end()) return std::string();
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return std::string();
    return row[j];
}

static int ToInteger(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

/* ============================
 * Triplet formatter (unused here)
 * ============================ */

static std::string make_triplet(double v, double stat, double syst) {
    std::ostringstream oss;
    oss << v << "," << stat << "," << syst;
    return oss.str();
}

/* ============================
 * Build the NEW header (schema)
 * ============================ */

static const std::vector<std::string>& avg_groups() {
    static const std::vector<std::string> g = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out",
        "Fa18", "Sp18", "10.6 GeV"
    };
    return g;
}

static void add_grouped_avg_columns(std::vector<std::string>& H, const std::string& base) {
    for (const auto& g : avg_groups()) {
        std::ostringstream name;
        name << base << ", " << g;
        H.push_back(name.str());
    }
}

static void add_theta_group_columns_after_phiavg_10p6(std::vector<std::string>& H) {
    // Inserted immediately after "phiavg, 10.6 GeV" (which is last in avg_groups()).

    // e_theta sequence EXACTLY as requested (note: Sp18 Inb intentionally omitted)
    H.push_back("e_theta, Fa18 Inb");
    H.push_back("e_theta, Fa18 Out");
    H.push_back("e_theta, Sp19 Inb");
    H.push_back("e_theta, Sp18 Out");
    H.push_back("e_theta, Fa18");
    H.push_back("e_theta, Sp18");
    H.push_back("e_theta, 10.6 GeV");

    // p_theta sequence EXACTLY as requested
    H.push_back("p_theta, Fa18 Inb");
    H.push_back("p_theta, Fa18 Out");
    H.push_back("p_theta, Sp19 Inb");
    H.push_back("p_theta, Sp18 Inb");
    H.push_back("p_theta, Sp18 Out");
    H.push_back("p_theta, Fa18");
    H.push_back("p_theta, Sp18");
    H.push_back("p_theta, 10.6 GeV");

    // g_theta sequence EXACTLY as requested
    H.push_back("g_theta, Fa18 Inb");
    H.push_back("g_theta, Fa18 Out");
    H.push_back("g_theta, Sp19 Inb");
    H.push_back("g_theta, Sp18 Inb");
    H.push_back("g_theta, Sp18 Out");
    H.push_back("g_theta, Fa18");
    H.push_back("g_theta, Sp18");
    H.push_back("g_theta, 10.6 GeV");
}

static void add_bin_definition_columns(std::vector<std::string>& H) {
    // Bin-definition columns (copied where specified)
    H.push_back("bin index");
    H.push_back("Bin Name");

    // xB
    H.push_back("xBmin");
    H.push_back("xBmax");
    add_grouped_avg_columns(H, "xBavg"); // immediately after xBmax

    // Q2
    H.push_back("Q2min");
    H.push_back("Q2max");
    add_grouped_avg_columns(H, "Q2avg"); // immediately after Q2max

    // t_abs
    H.push_back("t_abs_min");
    H.push_back("t_abs_max");
    add_grouped_avg_columns(H, "t_abs_avg"); // immediately after t_abs_max

    // phi
    H.push_back("phimin");
    H.push_back("phimax");
    add_grouped_avg_columns(H, "phiavg"); // immediately after phimax

    // NEW: insert theta columns immediately after "phiavg, 10.6 GeV"
    add_theta_group_columns_after_phiavg_10p6(H);

    // remaining bin-definition helpers
    H.push_back("tmin");
    H.push_back("tcol");
    H.push_back("P1");
    H.push_back("P2");
    H.push_back("maxP1P2");
    H.push_back("minP1P2");
    H.push_back("intP1P2");
}

static void add_raw_yield_columns(std::vector<std::string>& H) {
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    const std::vector<std::string> topos = { "(FD, FD)", "(CD, FD)", "(CD, FT)" };

    // ep->epg
    for (const auto& per : periods) {
        for (const auto& topo : topos) {
            for (const auto& hel : helicities) {
                std::ostringstream name;
                name << "raw yield, ep->epg, " << topo << ", exp, " << per << ", " << hel;
                H.push_back(name.str());
            }
        }
    }
    // ep->eppi0
    for (const auto& per : periods) {
        for (const auto& topo : topos) {
            for (const auto& hel : helicities) {
                std::ostringstream name;
                name << "raw yield, ep->eppi0, " << topo << ", exp, " << per << ", " << hel;
                H.push_back(name.str());
            }
        }
    }
}

static void add_contamination_columns(std::vector<std::string>& H) {
    const std::vector<std::string> per = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    for (const auto& p : per) {
        std::ostringstream name;
        name << "contamination ratio, " << p;
        H.push_back(name.str());
    }
}

static void add_signal_yield_columns(std::vector<std::string>& H) {
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    for (const auto& per : periods) {
        for (const auto& hel : helicities) {
            std::ostringstream name;
            name << "signal yield, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }
}

static void add_acceptance_columns(std::vector<std::string>& H) {
    const std::vector<std::string> per = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    for (const auto& p : per) {
        std::ostringstream name;
        name << "acceptance, " << p;
        H.push_back(name.str());
    }
}

static void add_acc_corrected_yield_columns(std::vector<std::string>& H) {
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    const std::vector<std::string> combined = {
        "Fa18", "Sp18", "2018 (10.6 GeV)"
    };

    // Per period
    for (const auto& per : periods) {
        for (const auto& hel : helicities) {
            std::ostringstream name;
            name << "acceptance corrected yield, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }
    // Combined
    for (const auto& grp : combined) {
        for (const auto& hel : helicities) {
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
    // EXACT insertion sequence requested
    H.push_back("norm, Fa18 Inb");
    H.push_back("norm, Fa18 Out");
    H.push_back("norm, Sp19 Inb");
    H.push_back("norm, Sp18 Inb");
    H.push_back("norm, Sp18 Out");
    H.push_back("norm, Fa18");
    H.push_back("norm, Sp18");
    H.push_back("norm, 10.6 GeV");
}

static void add_normed_cross_section_columns(std::vector<std::string>& H) {
    // Repeat the SAME cross section columns and SAME ordering, but with "normed " prefixed.
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    const std::vector<std::string> combined = {
        "10.6 GeV", "Fa18", "Sp18"
    };

    for (const auto& per : periods) {
        for (const auto& hel : helicities) {
            std::ostringstream name;
            name << "normed cross sections, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }
    for (const auto& grp : combined) {
        for (const auto& hel : helicities) {
            std::ostringstream name;
            name << "normed cross sections, ep->epg, exp, " << grp << ", " << hel;
            H.push_back(name.str());
        }
    }
}

static void add_cross_section_columns(std::vector<std::string>& H) {
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    const std::vector<std::string> combined = {
        "10.6 GeV", "Fa18", "Sp18"
    };

    // Per period
    for (const auto& per : periods) {
        for (const auto& hel : helicities) {
            std::ostringstream name;
            name << "cross sections, ep->epg, exp, " << per << ", " << hel;
            H.push_back(name.str());
        }
    }
    // Combined
    for (const auto& grp : combined) {
        for (const auto& hel : helicities) {
            std::ostringstream name;
            name << "cross sections, ep->epg, exp, " << grp << ", " << hel;
            H.push_back(name.str());
        }
    }

    // NEW: insert norm columns immediately after "cross sections, ep->epg, exp, Sp18, neg"
    // Since the above loops end on (grp="Sp18", hel="neg"), appending here is exactly "right after".
    add_norm_columns(H);

    // NEW: immediately after norm columns, append "normed" copies of all cross section columns.
    add_normed_cross_section_columns(H);
}

static void add_bsa_columns(std::vector<std::string>& H) {
    const std::vector<std::string> counts_groups = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out",
        "Fa18", "Sp18", "10.6 GeV"
    };
    for (const auto& g : counts_groups) {
        std::ostringstream name;
        name << "BSA, counts, " << g;
        H.push_back(name.str());
    }
    for (const auto& g : counts_groups) {
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
    H.reserve(650);

    add_bin_definition_columns(H);  // includes grouped avg columns + NEW theta group columns after phiavg, 10.6 GeV
    add_raw_yield_columns(H);
    add_contamination_columns(H);
    add_signal_yield_columns(H);
    add_acceptance_columns(H);
    add_acc_corrected_yield_columns(H);
    add_frad_fbin_binvol_columns(H);
    add_luminosity_columns(H);
    add_cross_section_columns(H);   // includes NEW norm columns + NEW normed cross section columns
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

    // Bin-definition copies (no singleton avgs in new schema)
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

    // valid bin
    M.pairs.push_back({"valid bin", "valid bin"});

    // Prefactors and coefficients (copy into 10.6 GeV versions)
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
    // 1) Open Lee's CSV
    std::ifstream fin(lee_csv_path);
    if (!fin.is_open()) {
        std::cerr << "[init] ERROR: Could not open Lee CSV: " << lee_csv_path << std::endl;
        return false;
    }

    // 2) Read header row from Lee
    std::string lee_header_line;
    if (!std::getline(fin, lee_header_line)) {
        std::cerr << "[init] ERROR: Lee CSV appears to be empty: " << lee_csv_path << std::endl;
        return false;
    }
    std::vector<std::string> lee_header = split_csv_line(lee_header_line);

    // Normalize first column name if blank/unnamed: call it "bin index"
    if (!lee_header.empty()) {
        std::string first = lee_header[0];
        bool unnamed = first.empty();
        if (!unnamed) {
            std::string low;
            low.reserve(first.size());
            for (char c : first) low.push_back((char)std::tolower((unsigned char)c));
            if (low.find("unnamed") != std::string::npos) unnamed = true;
        }
        if (unnamed) lee_header[0] = "bin index";
    }

    // 3) Build Lee's header index
    std::unordered_map<std::string,int> L = build_header_index(lee_header);

    // 4) Build the new header schema
    std::vector<std::string> new_header = build_new_header();

    // 5) Open output CSV and write new header
    std::ofstream fout(out_csv_path);
    if (!fout.is_open()) {
        std::cerr << "[init] ERROR: Could not open output CSV for writing: " << out_csv_path << std::endl;
        return false;
    }
    fout << join_csv_row(new_header) << "\n";

    // 6) Prepare copy map
    CopyMap cmap = build_copy_map();

    // 7) Read Lee rows, filter on "valid bin" == 1, and write new rows
    int input_rows = 0;
    int kept_rows  = 0;

    std::string line;
    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line(line);

        // Fetch "valid bin" from Lee
        std::string valid_s = get_col(cols, L, "valid bin");
        int valid = ToInteger(valid_s);

        // Only keep valid bins
        if (valid != 1) continue;

        // Create a blank output row with the right number of columns
        std::vector<std::string> out_row(new_header.size(), std::string());

        // Copy fields defined in the copy map
        for (const auto& kv : cmap.pairs) {
            const std::string& new_name = kv.first;
            const std::string& old_name = kv.second;

            auto it_new = std::find(new_header.begin(), new_header.end(), new_name);
            if (it_new == new_header.end()) {
                continue;
            }
            size_t new_idx = (size_t)std::distance(new_header.begin(), it_new);

            std::string v = get_col(cols, L, old_name);
            out_row[new_idx] = v;
        }

        // All analysis-produced fields (including grouped averages, theta columns, norms, and normed cross sections) remain blank at init
        fout << join_csv_row(out_row) << "\n";
        ++kept_rows;
    }

    std::cout << "[init] Lee rows read: " << input_rows << "\n";
    std::cout << "[init] Valid rows kept (valid bin == 1): " << kept_rows << "\n";
    std::cout << "[init] Wrote output: " << out_csv_path << "\n";

    return true;
}