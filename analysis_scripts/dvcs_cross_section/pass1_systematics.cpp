#include "pass1_systematics_import.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct CsvTable {
    std::vector<std::string> header;
    std::vector<std::vector<std::string> > rows;
    std::unordered_map<std::string, int> index;
};

static std::string trim(const std::string& s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) {
        ++b;
    }
    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) {
        --e;
    }
    return s.substr(b, e - b);
}

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
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

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

static void assert_no_duplicate_columns(const std::vector<std::string>& header,
                                        const std::string& path) {
    std::set<std::string> seen;
    std::set<std::string> duplicates;
    for (const auto& name : header) {
        if (!seen.insert(name).second) {
            duplicates.insert(name);
        }
    }
    if (!duplicates.empty()) {
        std::ostringstream msg;
        msg << "CSV header contains duplicate columns in " << path << ":";
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
    if (!table.header.empty()) {
        std::string first = table.header[0];
        std::string low;
        low.reserve(first.size());
        for (const char c : first) {
            low.push_back((char)std::tolower((unsigned char)c));
        }
        if (first.empty() || low.find("unnamed") != std::string::npos) {
            table.header[0] = "bin index";
        }
    }

    assert_no_duplicate_columns(table.header, path);
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
            msg << "CSV row has " << row.size() << " columns but header has "
                << table.header.size() << " columns in " << path;
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
            throw std::runtime_error("Internal CSV row/header size mismatch while writing " + path);
        }
        fout << join_csv_row(row) << "\n";
    }

    fout.close();
    if (!fout) {
        throw std::runtime_error("Failed while writing temporary CSV: " + tmp_path);
    }

    std::remove(path.c_str());
    if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
        throw std::runtime_error("Could not replace CSV with temporary file: " + path);
    }
}

static bool has_column(const CsvTable& table, const std::string& name) {
    return table.index.find(name) != table.index.end();
}

static int column_or_throw(const CsvTable& table,
                           const std::string& name,
                           const std::string& context) {
    const auto it = table.index.find(name);
    if (it == table.index.end()) {
        throw std::runtime_error("Missing required column for " + context + ": " + name);
    }
    return it->second;
}

static std::string get_cell(const std::vector<std::string>& row,
                            const CsvTable& table,
                            const std::string& name) {
    const auto it = table.index.find(name);
    if (it == table.index.end()) {
        return std::string();
    }
    const int col = it->second;
    if (col < 0 || col >= (int)row.size()) {
        return std::string();
    }
    return row[(size_t)col];
}

static double to_double(const std::string& s) {
    const std::string t = trim(s);
    if (t.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    char* endp = nullptr;
    const double v = std::strtod(t.c_str(), &endp);
    if (endp == t.c_str()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return v;
}

static int to_int(const std::string& s) {
    const std::string t = trim(s);
    if (t.empty()) {
        return 0;
    }
    char* endp = nullptr;
    const long v = std::strtol(t.c_str(), &endp, 10);
    if (endp == t.c_str()) {
        return 0;
    }
    return (int)v;
}

static std::string format_double(double value) {
    if (!std::isfinite(value)) {
        return std::string();
    }
    std::ostringstream oss;
    oss << std::setprecision(12) << value;
    return oss.str();
}

static const std::vector<std::string>& pass1_systematic_component_columns() {
    static const std::vector<std::string> cols = {
        "Syst. err (pi0 subtraction)",
        "Syst. err (Acceptance)",
        "Syst.err (Frad)",
        "Syst.err (Fbin)"
    };
    return cols;
}

static const std::string& pass1_systematic_total_column() {
    static const std::string col = "Syst. err (point-to-point total)";
    return col;
}

static const std::vector<std::string>& pass1_systematic_destination_columns() {
    static const std::vector<std::string> cols = {
        "Syst. err (pi0 subtraction)",
        "Syst. err (Acceptance)",
        "Syst.err (Frad)",
        "Syst.err (Fbin)",
        "Syst. err (point-to-point total)"
    };
    return cols;
}

static double quadrature_total_from_components(const std::map<std::string, std::string>& values) {
    double sum2 = 0.0;
    bool any = false;
    for (const auto& col : pass1_systematic_component_columns()) {
        const auto it = values.find(col);
        if (it == values.end()) {
            continue;
        }
        const double v = to_double(it->second);
        if (!std::isfinite(v)) {
            continue;
        }
        sum2 += v * v;
        any = true;
    }
    return any ? std::sqrt(sum2) : std::numeric_limits<double>::quiet_NaN();
}

struct Pass1SystValues {
    bool used = false;
    std::map<std::string, std::string> values;
};

static std::string make_boundary_key(double xb_min,
                                     double xb_max,
                                     double q2_min,
                                     double q2_max,
                                     double t_min,
                                     double t_max,
                                     double phi_min,
                                     double phi_max) {
    std::ostringstream key;
    key << std::fixed << std::setprecision(6)
        << xb_min << '|'
        << xb_max << '|'
        << q2_min << '|'
        << q2_max << '|'
        << t_min << '|'
        << t_max << '|'
        << phi_min << '|'
        << phi_max;
    return key.str();
}

static std::string make_boundary_key_from_pass1_row(const CsvTable& pass1,
                                                    const std::vector<std::string>& row) {
    return make_boundary_key(to_double(get_cell(row, pass1, "xBmin")),
                             to_double(get_cell(row, pass1, "xBmax")),
                             to_double(get_cell(row, pass1, "Q2min")),
                             to_double(get_cell(row, pass1, "Q2max")),
                             to_double(get_cell(row, pass1, "tmin")),
                             to_double(get_cell(row, pass1, "tmax")),
                             to_double(get_cell(row, pass1, "phimin")),
                             to_double(get_cell(row, pass1, "phimax")));
}

static std::string make_boundary_key_from_pass2_row(const CsvTable& pass2,
                                                    const std::vector<std::string>& row) {
    return make_boundary_key(to_double(get_cell(row, pass2, "xBmin")),
                             to_double(get_cell(row, pass2, "xBmax")),
                             to_double(get_cell(row, pass2, "Q2min")),
                             to_double(get_cell(row, pass2, "Q2max")),
                             to_double(get_cell(row, pass2, "t_abs_min")),
                             to_double(get_cell(row, pass2, "t_abs_max")),
                             to_double(get_cell(row, pass2, "phimin")),
                             to_double(get_cell(row, pass2, "phimax")));
}

static void validate_schema(const CsvTable& pass2,
                            const CsvTable& pass1) {
    column_or_throw(pass2, "bin index", "pass-2 systematic import");
    column_or_throw(pass2, "xBmin", "pass-2 systematic import");
    column_or_throw(pass2, "xBmax", "pass-2 systematic import");
    column_or_throw(pass2, "Q2min", "pass-2 systematic import");
    column_or_throw(pass2, "Q2max", "pass-2 systematic import");
    column_or_throw(pass2, "t_abs_min", "pass-2 systematic import");
    column_or_throw(pass2, "t_abs_max", "pass-2 systematic import");
    column_or_throw(pass2, "phimin", "pass-2 systematic import");
    column_or_throw(pass2, "phimax", "pass-2 systematic import");

    column_or_throw(pass1, "bin index", "pass-1 systematic import");
    column_or_throw(pass1, "xBmin", "pass-1 systematic import");
    column_or_throw(pass1, "xBmax", "pass-1 systematic import");
    column_or_throw(pass1, "Q2min", "pass-1 systematic import");
    column_or_throw(pass1, "Q2max", "pass-1 systematic import");
    column_or_throw(pass1, "tmin", "pass-1 systematic import");
    column_or_throw(pass1, "tmax", "pass-1 systematic import");
    column_or_throw(pass1, "phimin", "pass-1 systematic import");
    column_or_throw(pass1, "phimax", "pass-1 systematic import");

    for (const auto& col : pass1_systematic_component_columns()) {
        column_or_throw(pass1, col, "pass-1 systematic source");
    }
    for (const auto& col : pass1_systematic_destination_columns()) {
        column_or_throw(pass2, col, "pass-2 systematic destination");
    }
}

} // namespace

bool import_pass1_systematics(const std::string& csv_path,
                              const std::string& pass1_summary_path) {
    try {
        CsvTable pass2 = read_csv_or_throw(csv_path);
        const CsvTable pass1 = read_csv_or_throw(pass1_summary_path);

        validate_schema(pass2, pass1);

        std::map<int, Pass1SystValues> by_bin_index;
        std::map<std::string, Pass1SystValues> by_boundary_key;

        for (const auto& row : pass1.rows) {
            Pass1SystValues values;
            for (const auto& col : pass1_systematic_component_columns()) {
                values.values[col] = get_cell(row, pass1, col);
            }
            values.values[pass1_systematic_total_column()] =
                format_double(quadrature_total_from_components(values.values));

            const int bin_index = to_int(get_cell(row, pass1, "bin index"));
            if (bin_index > 0) {
                by_bin_index[bin_index] = values;
            }

            const std::string key = make_boundary_key_from_pass1_row(pass1, row);
            if (!key.empty()) {
                by_boundary_key[key] = values;
            }
        }

        int matched_by_index = 0;
        int matched_by_boundaries = 0;
        int unmatched = 0;

        for (auto& row : pass2.rows) {
            const int bin_index = to_int(get_cell(row, pass2, "bin index"));

            const Pass1SystValues* values = nullptr;
            auto it_index = by_bin_index.find(bin_index);
            if (it_index != by_bin_index.end()) {
                values = &it_index->second;
                ++matched_by_index;
            } else {
                const std::string key = make_boundary_key_from_pass2_row(pass2, row);
                auto it_key = by_boundary_key.find(key);
                if (it_key != by_boundary_key.end()) {
                    values = &it_key->second;
                    ++matched_by_boundaries;
                }
            }

            if (!values) {
                ++unmatched;
                continue;
            }

            for (const auto& col : pass1_systematic_destination_columns()) {
                const int j = column_or_throw(pass2, col, "pass-2 systematic destination");
                const auto it_value = values->values.find(col);
                row[(size_t)j] = (it_value == values->values.end()) ? std::string() : it_value->second;
            }
        }

        write_csv_or_throw(csv_path, pass2);

        std::cout << "[pass1-systematics] Imported fixed pass-1 systematic components from "
                  << pass1_summary_path << "\n";
        std::cout << "[pass1-systematics] Rows matched by bin index: " << matched_by_index
                  << ", by boundaries: " << matched_by_boundaries
                  << ", unmatched: " << unmatched << "\n";
        std::cout << "[pass1-systematics] Filled columns:";
        for (const auto& col : pass1_systematic_destination_columns()) {
            std::cout << "\n  - " << col;
        }
        std::cout << "\n";

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[pass1-systematics] ERROR: " << e.what() << "\n";
        return false;
    }
}
