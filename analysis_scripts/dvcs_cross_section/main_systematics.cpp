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
#include <cctype>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include "run_period_consistency_systematics.h"

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
    return {
        "normed cross sections, ep->epg, exp, 10.6 GeV, unpol, combination sys",
        "normed cross sections, ep->epg, exp, Fa18, unpol, combination sys",
        "normed cross sections, ep->epg, exp, Fa18, pos, combination sys",
        "normed cross sections, ep->epg, exp, Fa18, neg, combination sys",
        "normed cross sections, ep->epg, exp, Sp18, unpol, combination sys",
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol, combination sys"
    };
}

static void validate_combination_systematics_schema(const CsvTable& table) {
    require_columns(table,
                    combined_normed_cross_section_columns(),
                    "combined normed cross-section systematics inputs");

    require_columns(table,
                    combination_systematic_columns(),
                    "combined normed cross-section systematic outputs");
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

static bool update_combination_systematics_csv_shell(const std::string& csv_path) {
    CsvTable table = read_csv_or_throw(csv_path);

    validate_combination_systematics_schema(table);

    std::cout << "[systematics] CSV rows loaded: " << table.rows.size() << "\n";
    std::cout << "[systematics] CSV columns loaded: " << table.header.size() << "\n";
    std::cout << "[systematics] Combination-systematics schema validation passed.\n";
    std::cout << "[systematics] No combination-systematic values are written yet.\n";

    // Keep the writer compiled and ready for future update functions.  This is
    // intentionally disabled until the first real systematic calculation is
    // added, so running the shell does not rewrite the CSV.
    // write_csv_or_throw(csv_path, table);

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

        if (!update_combination_systematics_csv_shell(csv_main)) {
            std::cerr << "[systematics] FATAL: update_combination_systematics_csv_shell failed.\n";
            return 1;
        }

        if (!run_period_consistency_systematics(csv_main, "output/systematics")) {
            std::cerr << "[systematics] FATAL: run_period_consistency_systematics failed.\n";
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "[systematics] FATAL: " << e.what() << "\n";
        return 1;
    }

    std::cout << "[systematics] main_systematics complete.\n";
    return 0;
}