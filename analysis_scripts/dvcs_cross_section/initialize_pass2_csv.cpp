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

/**
 *
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
 * CSV FORMAT
 * ----------
 * We only need a simple CSV parser that can handle quoted fields and commas.
 * We'll implement a minimal split function that respects double quotes.
 *
 * TRIPLE FIELDS
 * -------------
 * Certain fields will later store triplets "value,stat,syst" inside a single CSV cell.
 * In this initializer, we leave those cells blank (empty string). Downstream code will
 * fill them with quoted triples (e.g. "12.3,0.5,0.7"). We provide a simple formatter below.
 */

/* ============================
 * Minimal CSV utilities
 * ============================ */

/**
 * split_csv_line
 *
 * Splits a CSV line into fields with very basic quote handling.
 * - If we are inside double quotes, commas are treated as literal commas.
 * - Quotes themselves are not included in the returned text.
 *
 */
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

/**
 * join_csv_row
 *
 * Writes a vector of fields to a CSV line.
 * - If a field contains a comma or a quote, we wrap it in quotes and
 *   escape any embedded quotes by doubling them.
 */
static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream oss;
    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];
        bool need_quotes = (s.find(',') != std::string::npos) || (s.find('"') != std::string::npos);
        if (need_quotes) {
            oss << '"';
            for (char c : s) {
                if (c == '"') oss << "\"\""; else oss << c;
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

/**
 * Build a case-sensitive name->index map for Lee's header.
 * We use exact names because Lee's CSV has stable column titles.
 */
static std::unordered_map<std::string, int>
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> m;
    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }
    return m;
}

/**
 * get_col
 *
 * Safe getter that returns empty string if a named column is missing
 * or out-of-range in the given row. We do not throw; we leave cells blank.
 */
static std::string get_col(const std::vector<std::string>& row,
                           const std::unordered_map<std::string,int>& idx,
                           const std::string& name) {
    auto it = idx.find(name);
    if (it == idx.end()) return std::string();
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) return std::string();
    return row[j];
}

/**
 * ToInteger
 *
 * Converts a string to integer. Returns 0 if conversion fails.
 * Used only for "valid bin" filtering.
 */
static int ToInteger(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

/* ============================
 * Triplet formatter
 * ============================
 *
 * We will leave all triplet fields blank at initialization.
 * This helper shows how we will format triplets later:
 *
 *   make_triplet(12.3, 0.5, 0.7)  ->  "12.3,0.5,0.7"
 *
 * If you want to populate a single-number-only field as a triplet,
 * you can do: make_triplet(value, NAN, NAN) or just leave it blank here.
 */
static std::string make_triplet(double v, double stat, double syst) {
    std::ostringstream oss;
    oss << v << "," << stat << "," << syst;
    return oss.str();
}

/* ============================
 * Build the NEW header (schema)
 * ============================
 *
 * We construct the header in a very explicit, readable way.
 * Repeated families of columns (periods, helicities, topologies, energies)
 * are created by simple loops for clarity and to avoid typos.
 */

static void add_bin_definition_columns(std::vector<std::string>& H) {
    // Bin-definition columns (copied where specified)
    H.push_back("bin index");
    H.push_back("Bin Name");
    H.push_back("xBmin");
    H.push_back("xBmax");
    H.push_back("xBavg");      // blank at init
    H.push_back("Q2min");
    H.push_back("Q2max");
    H.push_back("Q2avg");      // blank at init
    H.push_back("t_abs_min");
    H.push_back("t_abs_max");
    H.push_back("t_abs_avg");  // blank at init
    H.push_back("phimin");
    H.push_back("phimax");
    H.push_back("phiavg");     // blank at init
    H.push_back("tmin");
    H.push_back("tcol");
    H.push_back("P1");
    H.push_back("P2");
    H.push_back("maxP1P2");
    H.push_back("minP1P2");
    H.push_back("intP1P2");
}

static void add_raw_yield_columns(std::vector<std::string>& H) {
    // Periods we will support
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    // Helicity states for split
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    // Topologies
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
    // Triplets; shared for helicities
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
    // Helicity-split signal yield for ep->epg
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
    // Triplets; shared for helicities
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
    // Acceptance-corrected yield, ep->epg, helicity split, per period and combined groups
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
    // Scalars; by beam energy
    H.push_back("Frad, 10.6 GeV");
    H.push_back("Frad, 10.2 GeV");
    H.push_back("Fbin, 10.6 GeV");
    H.push_back("Fbin, 10.2 GeV");
    H.push_back("bin_volume, 10.6 GeV");
    H.push_back("bin_volume, 10.2 GeV");
}

static void add_luminosity_columns(std::vector<std::string>& H) {
    // Scalars; per period and combined
    H.push_back("integrated luminosity, Fa18 Inb (nC)");
    H.push_back("integrated luminosity, Fa18 Out (nC)");
    H.push_back("integrated luminosity, Sp19 Inb (nC)");
    H.push_back("integrated luminosity, Sp18 Inb (nC)");
    H.push_back("integrated luminosity, Sp18 Out (nC)");
    H.push_back("integrated luminosity, Fa18 (nC)");
    H.push_back("integrated luminosity, Sp18 (nC)");
    H.push_back("integrated luminosity, 10.6 GeV (nC)");
}

static void add_cross_section_columns(std::vector<std::string>& H) {
    // Cross sections, ep->epg, helicity split, per period and combined groups
    const std::vector<std::string> periods = {
        "Fa18 Inb", "Fa18 Out", "Sp19 Inb", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> helicities = { "unpol", "pos", "neg" };
    const std::vector<std::string> combined = {
        "Fa18", "Sp18", "10.6 GeV"
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
}

static void add_bsa_columns(std::vector<std::string>& H) {
    // Triplets; counts and sigma flavors, by group
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
    // valid bin
    H.push_back("valid bin");

    // Kinematic prefactors and Fourier coefficients
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
    H.reserve(300); // generous

    add_bin_definition_columns(H);
    add_raw_yield_columns(H);
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
 * ============================
 *
 * We list the exact column names expected in Lee's CSV for the fields we COPY.
 * Everything else in the new schema will be left blank.
 */

struct CopyMap {
    // Map new-name -> old-name
    std::vector<std::pair<std::string, std::string>> pairs;
};

static CopyMap build_copy_map() {
    CopyMap M;

    // Bin-definition copies
    M.pairs.push_back({"bin index", "bin index"});
    M.pairs.push_back({"Bin Name", "Bin Name"});
    M.pairs.push_back({"xBmin", "xBmin"});
    M.pairs.push_back({"xBmax", "xBmax"});
    // xBavg blank by design
    M.pairs.push_back({"Q2min", "Q2min"});
    M.pairs.push_back({"Q2max", "Q2max"});
    // Q2avg blank by design
    M.pairs.push_back({"t_abs_min", "t_abs_min"});
    M.pairs.push_back({"t_abs_max", "t_abs_max"});
    // t_abs_avg blank by design
    M.pairs.push_back({"phimin", "phimin"});
    M.pairs.push_back({"phimax", "phimax"});
    // phiavg blank by design
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
            // Some CSV tools export "Unnamed: 0" for the first column; treat that as "bin index".
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

            // Find index of new_name in new_header
            auto it_new = std::find(new_header.begin(), new_header.end(), new_name);
            if (it_new == new_header.end()) {
                // Should not happen; skip if we made a typo
                continue;
            }
            size_t new_idx = (size_t)std::distance(new_header.begin(), it_new);

            // Extract from Lee row
            std::string v = get_col(cols, L, old_name);
            out_row[new_idx] = v; // scalar copy as-is
        }

        // All other fields remain blank at initialization stage.
        // Triplet-type cells will be filled later as "value,stat,syst" by downstream code.

        fout << join_csv_row(out_row) << "\n";
        ++kept_rows;
    }

    std::cout << "[init] Lee rows read: " << input_rows << "\n";
    std::cout << "[init] Valid rows kept (valid bin == 1): " << kept_rows << "\n";
    std::cout << "[init] Wrote output: " << out_csv_path << "\n";

    return true;
}